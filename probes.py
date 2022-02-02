# from cctbx import geometry_restraints
from libtbx.phil import parse
from libtbx import easy_pickle
from scitbx.array_family import flex
from scitbx import matrix
from iotbx.pdb import common_residue_names_water as WATER_RES_NAMES
import time
import numpy as np

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

class MapModelCC(Probe):
  """Given any type of map, get the CC between the residue and map. This class
  should not be used directly -- use a subclass where validate_expedition and
  probe_at_residue are defined."""
  def probe_at_residue_by_expt(self, residue, chain_id, resid, resname, structure_type, expt_type, map_type):
    # get mask around residue
    box_mmm = self.expedition.mmm.extract_all_maps_around_model(
      selection_string="chain {chain} and resname {resname} and resseq {resid}".format(chain=chain_id, resname=resname, resid=resid)
    )
    box_mmm.create_mask_around_atoms()
    box_mmm.apply_mask_to_maps()
    cc = box_mmm.map_map_cc('{type}_expt','{type}_fmodel'.format(type=expt_type))
    return cc
  def validate_expedition(self, expt_type):
    if not expt_type in self.expedition.maps.keys() or \
      not 'expt' in self.expedition.maps[expt_type] or \
      not 'fmodel' in self.expedition.maps[expt_type]:
      return False
    return True

class CryoMapModelCC(MapModelCC):
  """MapModelCC for cryoEM maps"""
  def probe_at_residue(self, residue, chain_id, resid, resname, structure_type):
    MapModelCC.probe_at_residue_by_expt(self, residue, chain_id, resid, resname, structure_type, 'electron', 'expt')
  def validate_expedition(self):
    return MapModelCC.validate_expedition(self, 'electron')

class XrayMapModelCC(MapModelCC):
  """MapModelCC for X-ray maps"""
  def probe_at_residue(self, residue, chain_id, resid, resname, structure_type):
    MapModelCC.probe_at_residue_by_expt(self, residue, chain_id, resid, resname, structure_type, 'xray', 'expt')
  def validate_expedition(self):
    return MapModelCC.validate_expedition(self, 'xray')

class NeutronMapModelCC(MapModelCC):
  """MapModelCC for neutron maps"""
  def probe_at_residue(self, residue, chain_id, resid, resname, structure_type):
    MapModelCC.probe_at_residue_by_expt(self, residue, chain_id, resid, resname, structure_type, 'neutron', 'expt')
  def validate_expedition(self):
    return MapModelCC.validate_expedition(self, 'xray')

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




