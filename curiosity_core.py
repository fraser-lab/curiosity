from __future__ import division
from scitbx.array_family import flex
# from iotbx.data_manager import DataManager
from mmtbx.masks import manager as masks_manager
from mmtbx.f_model import manager as fmodel_manager
from load_probes import probe_collection
import time

# Curiosity: a tool for exploring electron density or eletrostatic potential maps.

class Expedition(object):
  """Step through every modeled residue and search for recognized features in the
  electron density or electrostatic potential map. The search should be abstracted
  away such that instances of this class know how to begin a search but do not
  dictate what that entails. This class describes how to traverse the model and
  how to keep track of any discoveries."""
  def __init__(self, map_model_manager_in, params, fmodel=None, fobs=None):
    self.mmm = map_model_manager_in
    # self.fobs = fobs
    # self.fmodel = fmodel
    self.params = params
    self.discoveries = {}
    self.n_tested = 0
    self.n_true_positives = {}
    self.n_true_negatives = {}
    self.n_false_positives = {}
    self.n_false_negatives = {}
    self.model = self.mmm.model()
    self.hier = self.model.get_hierarchy()
    print ("... constructed model hierarchy at {timestr}".format(timestr=time.asctime()))
    # only needed when using the non-ML, Phenix code for identifying ions
    # self.model.process(make_restraints=True)
    # print ("... processed input model with restraints at {timestr}".format(timestr=time.asctime()))
    self.grid_atoms()
    print ("... divided model atoms by grid position at {timestr}".format(timestr=time.asctime()))
    self.inventory_maps()
    print ("... inventoried maps at {timestr}".format(timestr=time.asctime()))
    self.inventory_probes()
    print ("... inventoried probes at {timestr}".format(timestr=time.asctime()))
    # self.setup_fmodel()

  # def setup_fmodel(self):


    # if not self.params.calculated_map_file:
      # compute fmodel
      # from mmtbx.utils import process_command_line_args
      # process_command_line_args(
      #   args=args, log=log, master_params=fmodel_from_xray_structure_master_params)
      # fmodel_miller_arrays = processed_args.reflection_files[0].as_miller_arrays()

  def grid_atoms(self):
    """Prepare an accessor for model atoms by position -- this will allow us to avoid
    screening the entire model every time we need to know what's nearby. Atoms are
    organized by blocks of 10x10x10 Angstroms in nested dictionaries."""
    grid = {}
    for chain in self.hier.chains():
      for residue in chain.residue_groups():
        for atom in residue.atoms():
          x,y,z = atom.xyz
          xx,yy,zz = map(int, (x//10,y//10,z//10))
          if not xx in grid.keys():
            grid[xx] = {}
          if not yy in grid[xx].keys():
            grid[xx][yy] = {}
          if not zz in grid[xx][yy].keys():
            grid[xx][yy][zz] = []
          grid[xx][yy][zz].append(atom)
    self.gridded_atoms = grid

  def locate_neighbors(self, position):
    x,y,z = position
    xx,yy,zz = map(int, (x//10,y//10,z//10))
    neighbors = self.gridded_atoms[xx][yy][zz][::] # DO NOT modify in place!
    # depending on which edges we're nearer, add the neighbors on the adjoining grid blocks as well
    def which_neighbors(dimension=1):
      if dimension == 1:
        direction = 1 if (x - xx >= 5) else -1
        return self.gridded_atoms[xx + direction][yy][zz]
      elif dimension == 2:
        direction = 1 if (y - yy >= 5) else -1
        return self.gridded_atoms[xx][yy + direction][zz]
      else:
        direction = 1 if (z - zz >= 5) else -1
        return self.gridded_atoms[xx][yy][zz + direction]
    for dimension in range(3):
      try:
        neighbors.extend(which_neighbors(dimension))
      except KeyError:
        continue
    return neighbors

  def inventory_maps(self):
    """Take an inventory of the types of maps available in the map_model_manager."""
    self.maps = {}
    self.managers = {}
    map_ids = self.mmm.map_id_list()
    for label in map_ids:
      if label == 'map_manager': continue
      try:
        etype, map_type = label.split("_")
      except ValueError:
        continue # could be 'mask'
      if etype not in self.maps.keys():
        self.maps[etype] = {}
      self.maps[etype][map_type] = self.mmm.get_map_manager_by_id(label)
    self.xray_structure = self.mmm.xray_structure()
    self.mmm.create_mask_around_atoms()
    basic_mask_manager = self.mmm._map_dict['mask']
    mask_coeffs = basic_mask_manager.map_as_fourier_coefficients()
    self.mask_manager = masks_manager(mask_coeffs,
                                      xray_structure=self.xray_structure)
    for etype in self.maps:
      if 'expt' in self.maps[etype].keys() and \
        'fmodel' in self.maps[etype].keys():
        fobs_coeffs = self.maps[etype]['expt'].map_as_fourier_coefficients().intensities()
        fmodel_coeffs = self.maps[etype]['fmodel'].map_as_fourier_coefficients()
        fmodel_mgr = fmodel_manager(f_obs          = fobs_coeffs,
                                    f_calc         = fmodel_coeffs,
                                    f_mask         = mask_coeffs,
                                    mask_manager   = self.mask_manager,
                                    xray_structure = self.xray_structure)
        self.managers[etype] = {}
        self.managers[etype]['fmodel'] = fmodel_mgr

  def experiments(self):
    """Return all available experiment types."""
    return self.maps.keys()

  def inventory_probes(self):
    """Check what tests we can run on this expedition."""
    self.probes = []
    for probe_class in probe_collection:
      probe = probe_class(self)
      if probe.validate_expedition():
        self.probes.append(probe)

  def add_map(self, mm, expt_type, map_type):
    """Update the map_model_manager with an additional map_manager."""
    mm.labels = [map_type]
    self.mmm.add_map_manager_by_id(mm, "_".join([expt_type, map_type]))
    # if we want to be able to add a map from file, we will need to
    # construct a map manager for it before passing it in
    self.inventory_maps()
    self.inventory_probes()

  def walk(self):
    """Traverse the model to access each residue. Heavily borrowing from qPTxM."""
    print ("Began main analysis (walk function) at {timestr}".format(timestr=time.asctime()))
    for chain in self.hier.chains():
      chain_id = chain.id.strip()
      if chain_id not in self.discoveries.keys():
        self.discoveries[chain_id] = []
      struct_type = "protein" if chain.is_protein() else "na" # nucleic acid
      # TODO: get rid of struct_type once we have a mechanism for determining
      # it from the resiud object itself (lookup tables)
      # below, we use the residue_groups() instead of the residues() accessor so
      # we don't lose the connection to the parent model. It's a little less
      # convenient to use but being able to copy residues is worth it. Also note:
      # same accessor for protein and nucleic acids.
      for residue in chain.residue_groups():
        resid = residue.resid().strip()
        resname = residue.unique_resnames()[0].strip()
        print ([a.i_seq for a in residue.atoms()])
        result_list = self.step(residue, chain_id, resid, resname, struct_type)
        self.n_tested += 1
        print ("... completed one step at {timestr}".format(timestr=time.asctime()))
        if result_list:
          self.discoveries[chain_id].append((residue, result_list))

  def step(self, residue, chain_id, resid, resname, struct_type):
    """Explore the density at a given residue."""
    result_list = []
    for probe in self.probes:
      if probe.validate_structure_type(struct_type) and probe.validate_residue_type(residue):
        result = probe.probe_at_residue(residue, chain_id, resid, resname, struct_type)
        if result is not None:
          result_list.append(result)
    return result_list

  def log_summary(self, outfile):
    """Log accuracy metrics if we presume the model to be ground truth."""
    with open(outfile, "w") as out:
      for probe_type in self.n_true_positives.keys():
        out.write("Results for probing {p}:\n".format(p=probe_type))
        out.write("Total number tested: {n}\n".format(n=self.n_tested))
        tp = self.n_true_positives[probe_type]
        tn = self.n_true_negatives[probe_type]
        fp = self.n_false_positives[probe_type]
        fn = self.n_false_negatives[probe_type]
        out.write("True positives: {n}\n".format(n=tp))
        out.write("True negatives: {n}\n".format(n=tn))
        out.write("False positives: {n}\n".format(n=fp))
        out.write("False negatives: {n}\n".format(n=fn))
        try:
          precision = tp/(tp+fp)
        except ZeroDivisionError:
          precision = "[undefined - no true positives or false positives identified]"
        try:
          recall = tp/(tp+fn)
        except ZeroDivisionError:
          recall = "[undefined = no true positives or false negatives identified]"
        # just in case we've screwed up record keeping somewhere,
        # we won't assume n_tested is correct for these calculations
        try:
          accuracy = (tp+tn)/(tp+tn+fp+fn)
        except ZeroDivisionError:
          accuracy = "[undefined]"
        try:
          f1 = 2*precision*recall/(precision + recall)
        except TypeError:
          f1 = "[undefined]"
        out.write("Precision: {p}\n".format(p=precision))
        out.write("Recall: {r}\n".format(r=recall))
        out.write("F1-score: {f1}\n".format(f1=f1))
        out.write("Overall accuracy: {a}\n\n".format(a=accuracy))
    return outfile

  def log_results(self, outfile):
    """Log discoveries to file. For now just simple strings."""
    with open(outfile, "w") as out:
      for chain in self.discoveries:
        for (residue, disco_list) in self.discoveries[chain]:
          # TODO: use functionality from pdb parsing tools instead
          resid = residue.resid().strip()
          resname = residue.unique_resnames()[0].strip()
          out.write(" ".join([chain, resid, resname] + disco_list) + "\n")
    return outfile

  # def apply_selected_discoveries(self, discoveries_dict):
    # """Given a dictionary of discoveries (proposed changes to the
    # model), make these changes and return the modified model."""
    # modified_model = self.hier.deepcopy()
    # Will need to be able to identify a modification function for each probe
    # and be able to tell which probe to use for each recorded discovery.
