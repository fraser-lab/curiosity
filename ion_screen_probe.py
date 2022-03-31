from probe_core import Probe
from iotbx.pdb import common_residue_names_water as WATER_RES_NAMES
from libtbx import easy_pickle
from scitbx.array_family import flex
from scitbx import matrix
import numpy as np
import time
import os

class AmIAnIonML(Probe):
  """Probe derived from Am I an Ion class for purposes of validation, but functioning
  completely differently. We'll locate nearest two neighbors to define a coordinate system
  against which to normalize. (Reject if there aren't at least two atoms within 5 Ang.)
  Then sample experimental density on a grid in that coordinate system. Print these
  values to a log file during training, or analyze immediately during use of a
  trained classifier."""
  def __init__(self, expedition):
    self.expedition = expedition
    self.experiment_type = None
    self.applicable_structure_types = None # allow all
    metals = ["Mg", "MG", "Fe", "FE", "Fe2", "FE2", "Fe3", "FE3", "Zn", "ZN", "Zn2", "ZN2"]
    metal_elements = ["MG", "FE", "ZN"]
    self.applicable_residue_list = WATER_RES_NAMES + metals
    self.applicable_coordination_partners = ["O", "N", "S"] + metal_elements
      # also check ions to see if they are ions, so we won't have to rename them as
      # waters in order to assess precision and recall
    self.water_and_ion_names = set([item.lower() for item in self.applicable_residue_list])

  def validate_expedition(self):
    for etype in self.expedition.maps:
      if 'expt' in self.expedition.maps[etype]:
        # right now we just grab whatever is the first map to meet all reqs
        self.experiment_type = etype
        self.setup_manager()
        return True
    return False

  def setup_manager(self):
    self.map_manager = self.expedition.maps[self.experiment_type]['expt']
    here = os.path.abspath(os.path.dirname(__file__))
    self.classifier = easy_pickle.load(os.path.join(here, "ml", "ion_detector_expt_ours_norm_tuned.pkl"))
    tmp = easy_pickle.load(os.path.join(here, "ml", "ion_lookup_expt_ours_norm.pkl"))
    self.lookup = {tmp[key]:key for key in tmp}
    # TODO for classifier:
    # - full path to pkl
    # - bundle release 0.1 and give it a zenodo doi
    # - update README to make clear this needs to run with py3
    print ("... completed probe setup at {timestr}".format(timestr=time.asctime()))

  def get_basis_set(self, position, radius=3):
    neighboring_atoms = self.expedition.locate_neighbors(position)
    coordinating_positions = [matrix.col(n.xyz) for n in neighboring_atoms \
      if n.element.upper() in self.applicable_coordination_partners]
    here = matrix.col(position)
    dists = [(here - other).length() for other in coordinating_positions]
    order = [dists.index(d) for d in sorted(dists)]
    sorted_positions = [coordinating_positions[i] for i in order if dists[i] <= radius]
    if len(sorted_positions) < 2:
      return # not enough information
    unit_vec_1 = matrix.col(sorted_positions[0]).normalize()
    unit_vec_3 = unit_vec_1.cross(matrix.col(sorted_positions[1])).normalize()
    unit_vec_2 = unit_vec_3.cross(unit_vec_1).normalize()
    return [unit_vec_1, unit_vec_2, unit_vec_3]

  def get_map_density_grid(self, origin, basis_set, grid_spacing, sampling_radius, map_manager):
    sampling_points = int(sampling_radius / grid_spacing)
    n_points_on_edge = sampling_points * 2 + 1
    grid_coords = flex.vec3_double(n_points_on_edge**3)
    sampling = range(-sampling_points, sampling_points+1)
    for i in sampling:
      for j in sampling:
        for k in sampling:
          position = matrix.col(origin) \
                      + grid_spacing * i * basis_set[0] \
                      + grid_spacing * j * basis_set[1] \
                      + grid_spacing * k * basis_set[2]
          index = i*n_points_on_edge**2 + j*n_points_on_edge + k
          grid_coords[index] = tuple(position)
    return map_manager.density_at_sites_cart(grid_coords)

  def probe_at_residue(self, residue, chain_id, resid, resname, struct_type, training=False):
    atoms = [atom for atom in residue.atoms() if atom.element.strip() != "H"]
    # we don't match to "O" just in case we are allowing ions to be re-identified here
    if not len(atoms) == 1:
      print ("too many atoms in candidate water")
      return
    if not resname.lower() in self.water_and_ion_names:
      print("{resname} does not appear to be a water".format(
        resname=resname))
      return
    position = atoms[0].xyz
    basis_set = self.get_basis_set(position)
    print("... determined orientation-normalized basis set at {timestr}".format(timestr=time.asctime()))
    if basis_set is None:
      return # not enough information
    density_grid = list(self.get_map_density_grid(position,
                                                  basis_set,
                                                  grid_spacing=0.5,
                                                  sampling_radius=3,
                                                  map_manager=self.map_manager))
    density_scale_factor = 1/(sum(density_grid)/len(density_grid))
    density_grid_normalized = [d*density_scale_factor for d in density_grid]
    print("... fetched grid of map densities at {timestr}".format(timestr=time.asctime()))
    if training:
      printable = map(str, density_grid_normalized)
      print ("... probed one water at {timestr}".format(timestr=time.asctime()))
      return ("density grid at water: " + " ".join(printable) + "\n")
    else:
      choice = self.lookup[self.classifier.predict(np.asarray(density_grid_normalized).reshape(1,-1))[0]]
      print ("... probed one water at {timestr}".format(timestr=time.asctime()))
      if not "ion" in self.expedition.n_true_positives.keys():
        self.expedition.n_true_positives["ion"] = 0
        self.expedition.n_true_negatives["ion"] = 0
        self.expedition.n_false_positives["ion"] = 0
        self.expedition.n_false_negatives["ion"] = 0
      if choice.upper() != resname.upper():
        if resname.upper() in WATER_RES_NAMES:
          self.expedition.n_false_positives["ion"] += 1
        else:
          self.expedition.n_false_negatives["ion"] += 1
        return (choice)
      else:
        if resname.upper() in WATER_RES_NAMES:
          self.expedition.n_true_negatives["ion"] += 1
        else:
          self.expedition.n_true_positives["ion"] += 1
        return