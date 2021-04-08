from __future__ import division
from iotbx.data_manager import DataManager
from probes import probe_collection

# Curiosity: a tool for exploring electron density or eletrostatic potential maps.

class Expedition(object):
  """Step through every modeled residue and search for recognized features in the
  electron density or electrostatic potential map. The search should be abstracted
  away such that instances of this class know how to begin a search but do not
  dictate what that entails. This class describes how to traverse the model and
  how to keep track of any discoveries."""
  def __init__(self, map_model_manager_in, params):
    self.mmm = map_model_manager_in
    self.params = params
    self.discoveries = {}
    self.hier = self.mmm.model().get_hierarchy()
    self.inventory_maps()
    self.inventory_probes()

  def inventory_maps(self):
    """Take an inventory of the types of maps available in the map_model_manager."""
    self.maps = {'expt':{}, 'diff':{}, 'fmodel':{}}
    self.experiments = []
    map_managers = self.mmm.map_managers()
    for mm in map_managers:
      etype = mm.experiment_type()
      if etype not in experiments:
        experiments.append(etype)
      if 'expt' in mm.labels:
        self.maps['expt'][etype] = mm
        continue
      if 'diff' in mm.labels:
        self.maps['diff'][etype] = mm
        continue
      if 'fmodel' in mm.labels:
        mm.set_experiment_type('fmodel')
        self.maps['fmodel']['fmodel'] = mm

  def inventory_probes(self):
    """Check what tests we can run on this expedition."""
    self.probes = []
    for probe_class in probe_collection:
      probe = probe_class(self)
      if probe.validate_expedition():
        self.probes.append(probe)

  def add_map(self, newmap, tag):
    """Update the map_model_manager with an additional map_manager."""
    self.mmm.add_map_manager_by_id(newmap, tag)
    # if we want to be able to add a map from file, we will need to
    # construct a map manager for it before passing it in
    self.inventory_maps()

  def walk(self):
    """Traverse the model to access each residue. Heavily borrowing from qPTxM."""
    for chain in self.hier.chains():
      chain_id = chain.id.strip()
      discoveries[chain_id] = []
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
        result_list = self.step(residue, struct_type=struct_type)
        if result_list:
          discoveries.append((residue, result_list))

  def step(self, residue, struct_type="protein"):
    """Explore the density at a given residue."""
    result_list = []
    for probe in self.probes:
      result = probe.probe_at_residue(residue, struct_type)
      if result is not None:
        result_list.append(result)
    return result_list

