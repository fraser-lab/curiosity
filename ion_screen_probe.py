from probe_core import Probe
from iotbx.pdb import common_residue_names_water as WATER_RES_NAMES
from libtbx import easy_pickle
from scitbx.array_family import flex
from scitbx import matrix
import numpy as np
import time

class AmIAnIon(Probe):
  """Determine if a water might be better modeled as an ion. Initialize once per
  model only, so we can do some preliminary work outside the probe method."""
  def __init__(self, expedition):
    self.expedition = expedition
    self.experiment_type = None
    self.applicable_structure_types = None # allow all
    self.applicable_residue_list = ["HOH","WAT","OH","DOD","D3O"]

  def validate_expedition(self):
    for etype in self.expedition.managers:
      if 'fmodel' in self.expedition.managers[etype] and \
        'expt' in self.expedition.maps[etype]:
        # right now we just grab whatever is the first map to meet all reqs
        self.experiment_type = etype
        self.setup_manager()
        return True
    return False

  def setup_manager(self):
    from mmtbx.ions.identify import create_manager, ion_identification_phil_str
    geometry_manager = self.expedition.model.get_restraints_manager().geometry
    params = parse(ion_identification_phil_str, process_includes=True).extract()
    self.manager = create_manager(self.expedition.hier,
                                  geometry_manager,
                                  self.expedition.managers[self.experiment_type]['fmodel'],
                                  wavelength=0,
                                  params=params)
    print ("... completed probe setup at {timestr}".format(timestr=time.asctime()))

  def probe_at_residue(self, residue, chain_id, resid, resname, struct_type):
    from mmtbx.ions.identify import AtomProperties
    atoms = [atom for atom in residue.atoms() if atom.name.strip() == "O"]
    if not len(atoms) == 1:
      print ("too many atoms in candidate water")
      return
    if not resname in WATER_RES_NAMES:
      print("{resname} does not appear to be a water".format(
        resname=resname))
      return
    waters = self.manager.water_selection()
    self.manager.atoms_to_props = dict((i_seq, AtomProperties(i_seq, self.manager)) for i_seq in waters)
    water_result = self.manager.analyze_water(atoms[0].i_seq)
    print ("... completed analysis of one water molecule at {timestr}".format(timestr=time.asctime()))
    if water_result:
      accepted = water_result.matching_candidates
      try:
        accepted = [a[0].element for a in accepted]
        print("accepted ions are logged as tuples") # FIXME
      except TypeError:
        accepted = [a.element for a in accepted]
        print("accepted ions are logged as objects") # FIXME
      rejected = [r[0].element for r in water_result.rejected_candidates]
      return ("ions accepted", accepted, "ions rejected", rejected)
    return

class AmIAnIonML(AmIAnIon):
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

  def setup_manager(self):
    # no manager at this time -- use gridded atoms to locate neighbors
    # # prepare an array of all water positions in the entire model
    # self.waters_x = flex.double()
    # self.waters_y = flex.double()
    # self.waters_z = flex.double()
    # for chain in self.expedition.hier.chains():
    #   for residue in chain.residue_groups():
    #     resname = residue.unique_resnames()[0].strip()
    #     if not resname.lower() in self.water_and_ion_names: continue
    #     oxygens = [atom for atom in residue.atoms() if atom.name.strip() == "O"]
    #     if not len(oxygens) == 1: continue
    #     position = oxygens[0].xyz
    #     self.waters_x.append(position[0])
    #     self.waters_y.append(position[1])
    #     self.waters_z.append(position[2])
    self.map_manager = self.expedition.maps[self.experiment_type]['expt']
    self.classifier = easy_pickle.load("mlp_classifier.pkl")
    # TODO for classifier:
    # - full path to pkl
    # - label file appropriately (e.g. classifier_waters_and_ions)
    # - bundle release 0.1 and give it a zenodo doi
    # - update README to make clear this needs to run with py3
    print ("... completed probe setup at {timestr}".format(timestr=time.asctime()))

  def get_basis_set(self, position, radius=3):
    # x,y,z = position
    # diffx = self.waters_x - x
    # diffy = self.waters_y - y
    # diffz = self.waters_z - z
    # nearby = diffx < radius and diffx > -1*radius and \
    #          diffy < radius and diffy > -1*radius and \
    #          diffz < radius and diffz > -1*radius
    #          # selects in a cube, not a sphere, for computational efficiency
    #          # also, using flex doubles instead of vec3_doubles at this stage
    #          # so that we can use overloaded vectorized operators on them
    # nearx = self.waters_x.select(nearby)
    # neary = self.waters_y.select(nearby)
    # nearz = self.waters_z.select(nearby)
    # diffx = diffx.select(nearby)
    # diffy = diffy.select(nearby)
    # diffz = diffz.select(nearby)
    # sq_dists = (diffx*diffx + diffy*diffy + diffz*diffz)
    # order = [list(sq_dists).index(d) for d in sorted(sq_dists)]
    # sorted_positions = flex.vec3_double([(nearx[i],neary[i],nearz[i]) for i in order if sq_dists[i] <= radius**2])
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
    density_grid = self.get_map_density_grid(position,
                                             basis_set,
                                             grid_spacing=0.5,
                                             sampling_radius=3,
                                             map_manager=self.map_manager)
    print("... fetched grid of map densities at {timestr}".format(timestr=time.asctime()))
    if training:
      printable = map(str, list(density_grid))
      print ("... probed one water at {timestr}".format(timestr=time.asctime()))
      return ("density grid at water: " + " ".join(printable) + "\n")
    else:
      choice = self.lookup[self.classifier.predict(np.asarray(density_grid).reshape(1,-1))[0]]
      print ("... probed one water at {timestr}".format(timestr=time.asctime()))
      if choice.upper() != resname.upper():
        return (choice)