from probe_core import Probe
from libtbx import easy_pickle
from scitbx.array_family import flex
from scitbx import matrix
from libtbx.utils import Sorry
import numpy as np
import time
import os

modifs_on_A = ["A", "MA6", "6MZ", "2MA", "A2M", "2M8"]
modifs_on_DA = ["DA"]
modifs_on_G = ["G", "G7M", "2MG", "M2G", "1MG", "OMG"]
modifs_on_DG = ["DG"]
purines = modifs_on_A + modifs_on_DA + modifs_on_G + modifs_on_DG
modifs_on_U = ["U", "UR3", "5MU", "OMU", "PSU"]
modifs_on_DT = ["DT"]
modifs_on_C = ["C", "5MC", "OMC", "40C"]
modifs_on_DC = ["DC"]
pyrimidines = modifs_on_U + modifs_on_DT + modifs_on_C + modifs_on_DC
nucleotides = purines + pyrimidines

def locate_atom_by_name(residue, name):
  """find the first copy of an atom by name on a hierarchical
  residue object"""
  atoms = residue.atom_groups()[0].atoms()
  for i in range(len(atoms)):
    if atoms[i].name.strip() == name.strip():
      return atoms[i]
  raise Sorry("couldn't locate requested atom")

def get_atom_positions(reference_atoms_tup, residue):
  reference_atoms = [
    locate_atom_by_name(residue, atom) for atom in reference_atoms_tup]
  reference_positions = [matrix.col(atom.xyz) for atom in reference_atoms]
  return reference_positions

class IsModifiedNucleotide(Probe):
  """Templated of AmIAnIonML, but applied to searching for posttranscriptional
  modifications (updated from qPTxM to use neural net in place of random
  forest classifier)."""
  def __init__(self, expedition):
    self.expedition = expedition
    self.experiment_type = None
    # self.applicable_structure_types = "nucleic acid"
    self.applicable_structure_types = None
    self.applicable_residue_list = nucleotides #["A","U","C","G","DA","DT","DC","DG"]

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
    self.classifier_A = easy_pickle.load(os.path.join(here, "ml", "modif_detector_A_tuned.pkl"))
    self.classifier_U = easy_pickle.load(os.path.join(here, "ml", "modif_detector_U_tuned.pkl"))
    self.classifier_C = easy_pickle.load(os.path.join(here, "ml", "modif_detector_C_tuned.pkl"))
    self.classifier_G = easy_pickle.load(os.path.join(here, "ml", "modif_detector_G_tuned.pkl"))
    # lookups are resname:int so we will need to reverse them each
    tmp = easy_pickle.load(os.path.join(here, "ml", "A_lookup.pkl"))
    self.lookup_A = {tmp[key]:key for key in tmp}
    tmp = easy_pickle.load(os.path.join(here, "ml", "U_lookup.pkl"))
    self.lookup_U = {tmp[key]:key for key in tmp}
    tmp = easy_pickle.load(os.path.join(here, "ml", "C_lookup.pkl"))
    self.lookup_C = {tmp[key]:key for key in tmp}
    tmp = easy_pickle.load(os.path.join(here, "ml", "G_lookup.pkl"))
    self.lookup_G = {tmp[key]:key for key in tmp}
    print ("... completed probe setup at {timestr}".format(timestr=time.asctime()))

  def get_type_and_origin_and_basis_set_nucleobase(self, residue, resname):
    """Define frame of reference for sampling map around this residue."""
    if resname in purines:#("A", "DA", "G", "DG"):
      reference_atoms_tuple = ("C4", "C5", "N3")
      base_type = "purine"
    elif resname in pyrimidines:#("U", "DT", "C", "DC"):
      reference_atoms_tuple = ("C2", "N3", "N1")
      base_type = "pyrimidine"
    else:
      raise Sorry("Nucleotide not recognized: {n}.".format(n=resname))
    ref1, ref2, ref3 = get_atom_positions(reference_atoms_tuple, residue)
    origin = (ref1 + ref2)/2.
    dir1 = ref2 - ref1
    unit_vec_1 = dir1.normalize()
    dir2 = ref3 - ref1
    unit_vec_3 = unit_vec_1.cross(matrix.col(dir2)).normalize()
    unit_vec_2 = unit_vec_3.cross(unit_vec_1).normalize()
    basis_set = [unit_vec_1, unit_vec_2, unit_vec_3]
    return (base_type, origin, basis_set)

  def get_origin_and_basis_set_sugar(self, residue):
    """Define frame of reference for sampling map around this residue."""
    reference_atoms_tuple = ("C3'", "C2'", "O3'")
    ref1, ref2, ref3 = get_atom_positions(reference_atoms_tuple, residue)
    dir1 = ref2 - ref1
    unit_vec_1 = dir1.normalize()
    dir2 = ref3 - ref1
    unit_vec_3 = unit_vec_1.cross(matrix.col(dir2)).normalize()
    unit_vec_2 = unit_vec_3.cross(unit_vec_1).normalize()
    basis_set = [unit_vec_1, unit_vec_2, unit_vec_3]
    return (ref1, basis_set)

  def get_map_density_grid(self, origin, basis_set, part):
    unit1, unit2, unit3 = basis_set
    # if part == "purine":
    #   dir1_min = -4*unit1
    #   dir1_max = 4.5*unit1
    #   dir2_min = -4.5*unit2
    #   dir2_max = 6*unit2
    #   dir3_min = -1*unit3
    #   dir3_max = unit3
    # elif part == "pyrimidine":
    #   dir1_min = -3.5*unit1
    #   dir1_max = 4.5*unit1
    #   dir2_min = -4.5*unit2
    #   dir2_max = 4.5*unit2
    #   dir3_min = -1*unit3
    #   dir3_max = unit3
    # elif part == "sugar":
    #   dir1_min = -2.5*unit1
    #   dir1_max = 4.5*unit1
    #   dir2_min = 0
    #   dir2_max = 3.5*unit2
    #   dir3_min = -2*unit3
    #   dir3_max = 3*unit3
    # else:
    #   raise NotImplementedError("Don't recognize this part")
    # grid_coords = flex.vec3_double()
    # for shift1 in range(dir1_min, dir1_max, 0.5):
    #   for shift2 in range(dir2_min, dir2_max, 0.5):
    #     for shift3 in range(dir3_min, dir3_max, 0.5):
    #       grid_coords.append(origin + shift1 + shift2 + shift3)
    # if part == "purine":
    #   dir1_min = -4
    #   dir1_max = 4.5
    #   dir2_min = -4.5
    #   dir2_max = 6
    #   dir3_min = -1
    #   dir3_max = 1
    # elif part == "pyrimidine":
    #   dir1_min = -3.5
    #   dir1_max = 4.5
    #   dir2_min = -4.5
    #   dir2_max = 4.5
    #   dir3_min = -1
    #   dir3_max = 1
    # elif part == "sugar":
    #   dir1_min = -2.5
    #   dir1_max = 4.5
    #   dir2_min = 0
    #   dir2_max = 3.5
    #   dir3_min = -2
    #   dir3_max = 3
    # else:
    #   raise NotImplementedError("Don't recognize this part")
    # grid_coords = flex.vec3_double()
    # for shift1 in range(dir1_min, dir1_max, 0.5):
    #   for shift2 in range(dir2_min, dir2_max, 0.5):
    #     for shift3 in range(dir3_min, dir3_max, 0.5):
    #       grid_coords.append(origin + shift1*unit1 + shift2*unit2 + shift3*unit3)
    if part == "purine":
      dir1_min = -8
      dir1_max = 9
      dir2_min = -9
      dir2_max = 12
      dir3_min = -2
      dir3_max = 2
    elif part == "pyrimidine":
      dir1_min = -7
      dir1_max = 9
      dir2_min = -9
      dir2_max = 9
      dir3_min = -2
      dir3_max = 2
    elif part == "sugar":
      dir1_min = -5
      dir1_max = 9
      dir2_min = 0
      dir2_max = 7
      dir3_min = -4
      dir3_max = 6
    else:
      raise NotImplementedError("Don't recognize this part")
    grid_coords = flex.vec3_double()
    for shift1 in range(dir1_min, dir1_max):
      for shift2 in range(dir2_min, dir2_max):
        for shift3 in range(dir3_min, dir3_max):
          grid_coords.append(origin + .5*shift1*unit1 + .5*shift2*unit2 + .5*shift3*unit3)
    return self.map_manager.density_at_sites_cart(grid_coords)

  def probe_at_residue(self, residue, chain_id, resid, resname, struct_type, training=False):
    base_type, base_origin, base_basis_set = \
      self.get_type_and_origin_and_basis_set_nucleobase(residue, resname)
    sugar_origin, sugar_basis_set = \
      self.get_origin_and_basis_set_sugar(residue)
    print("... determined orientation-normalized basis sets at {timestr}".format(timestr=time.asctime()))
    map_values_base = self.get_map_density_grid(base_origin, base_basis_set, base_type)
    map_values_sugar = self.get_map_density_grid(sugar_origin, sugar_basis_set, "sugar")
    print("... fetched grid of map densities at {timestr}".format(timestr=time.asctime()))
    density_grid = list(map_values_base) + list(map_values_sugar)
    density_scale_factor = 1/(sum(density_grid)/len(density_grid))
    density_grid_normalized = [d*density_scale_factor for d in density_grid]
    if training:
      printable = map(str, density_grid_normalized)
      print ("... probed one nucleotide at {timestr}".format(timestr=time.asctime()))
      return ("density grid at nucleotide: " + " ".join(printable) + "\n")
    else:
      if resname in modifs_on_A:
        classif = self.classifier_A
        lookup = self.lookup_A
      elif resname in modifs_on_U:
        classif = self.classifier_U
        lookup = self.lookup_U
      elif resname in modifs_on_C:
        classif = self.classifier_C
        lookup = self.lookup_C
      elif resname in modifs_on_G:
        classif = self.classifier_G
        lookup = self.lookup_G
      else:
        print ("we don't have a classifier for this nucleotide yet")
        return
      choice = lookup[classif.predict(np.asarray(density_grid_normalized).reshape(1,-1))[0]]
      print ("... probed one water at {timestr}".format(timestr=time.asctime()))
      if not "ptxm" in self.expedition.n_true_positives.keys():
        self.expedition.n_true_positives["ptxm"] = 0
        self.expedition.n_false_positives["ptxm"] = 0
        self.expedition.n_true_negatives["ptxm"] = 0
        self.expedition.n_false_negatives["ptxm"] = 0
      if choice.upper() != resname.upper():
        if resname.upper() in ["A", "U", "C", "G"]:
          self.expedition.n_false_positives["ptxm"] += 1
          # if the model is already correct, we've just discovered modifications that shouldn't actually be there
        else:
          self.expedition.n_false_negatives["ptxm"] += 1
          # if the model is already correct, we've just failed to discover a modification that should be modeled
        return (choice)
      else:
        if resname.upper() in ["A", "U", "C", "G"]:
          self.expedition.n_true_negatives["ptxm"] += 1
          # model and probe agree there should be no modification here
        else:
          self.expedition.n_true_positives["ptxm"] += 1
          # model and probe agree there should be a modification here
        return