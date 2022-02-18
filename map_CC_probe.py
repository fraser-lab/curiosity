from probe_core import Probe

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

