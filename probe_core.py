# from cctbx import geometry_restraints
from libtbx.phil import parse
from iotbx.pdb import common_residue_names_water as WATER_RES_NAMES

class Probe(object):
  """Scaffold what it means to probe a map at a residue depending on what we're
  looking for and what kind of map(s) we have available. Overwrite if we are
  placing any conditions on structure type (protein, nucleic acid) or residue
  type (e.g. only applicable to histidines)."""
  def __init__(self, expedition):
    self.expedition = expedition
    self.applicable_structure_types = None # allow all
    self.applicable_residue_list = None # allow all

  def validate_expedition(self):
    """Overwrite this in subclasses where there are conditions on which types
    of maps will be needed for this test."""
    return True

  def validate_structure_type(self, structure_type):
    """If specific structure types have been specified, ensure the residue
    we are examining is of one of those types."""
    if self.applicable_structure_types is None:
      return True
    return True if structure_type in self.applicable_structure_types else False

  def validate_residue_type(self, residue):
    """If specific residue types have been specified, ensure the residue
    we are examining is of one of those types."""
    if self.applicable_residue_list is None:
      return True
    return True if residue.unique_resnames()[0].strip() in self.applicable_residue_list else False

  def probe_at_residue(self, residue, chain_id, resid, resname, struct_type):
    """Overwrite this with how the probe should explore the map(s) at the residue
    and return a results object that can be meaningfully interpreted later. This
    method should run validate_structure_type and/or validate_residue_type if
    applicable."""
    pass

  def interpret_results(self, result):
    """Overwrite this with how the result object produced by probe_at_residue
    should be interpreted when the expedition is finished and reporting back."""
    pass
