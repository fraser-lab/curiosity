from __future__ import division

# Curiosity: a tool for exploring electron density or eletrostatic potential maps.

class ExploreModel(object):
  """Step through every modeled residue and search for recognized features in the
  electron density or electrostatic potential map. The search should be abstracted
  away such that instances of this class know how to begin a search but do not
  dictate what that entails. This class describes how to traverse the model and
  how to keep track of any discoveries."""
  def __init__(self, hierarchical_model_in, map_collection_in, params):
    self.hier = hierarchical_model_in # cctbx hierarchical model object
    self.maps = map_collection_in # should be passing this en masse to search step
    self.params = params
    self.discoveries = {}

  def walk(self):
    """Traverse the model to access each residue. Heavily borrowing from qPTxM."""
    for chain in self.hier.chains():
      chain_id = chain.id.strip()
      struct_type = "protein" if chain.is_protein() else "na" # nucleic acid
      i = 0 # index of the residue in the chain object
      # below, we use the residue_groups() instead of the residues() accessor so
      # we don't lose the connection to the parent model. It's a little less
      # convenient to use but being able to copy residues is worth it. Also note:
      # same accessor for protein and nucleic acids.
      for residue in chain.residue_groups():
        resid = residue.resid().strip()
        resname = residue.unique_resnames()[0].strip()
        self.step(residue, i, type=struct_type)
        i += 1

  def step(self, residue, i, type="protein"):
    """Explore the density at a given residue."""
    # TODO: probably should start with a CC of model to map
    # TODO: note any strong difference density
    # TODO: if the residue type is recognized and we have any features we know
    # to look for at that type of residue, look for them
    # TODO: keep track of results (use self.discoveries, include position in the
    # model and a label for what we're looking at, to use with Coot plugin)
    pass

