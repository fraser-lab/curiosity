from probe_core import Probe
import time
from scitbx.array_family import flex

class CompareDENQDensities(Probe):
  """For analysis of a ribosome where some asparagine and glutamine
  residues are expected to be aspartate and glutamate instead,
  analyze density at these four types of residues and store to a
  data object to compare between two structures."""
  def __init__(self, expedition):
      self.expedition = expedition
      self.applicable_structure_types = 'protein'
      self.applicable_residue_list = ['GLU', 'ASP', 'GLN', 'ASN', 'LEU']

  def validate_expedition(self):
    for etype in self.expedition.managers:
      if 'fmodel' in self.expedition.managers[etype] and \
        'expt' in self.expedition.maps[etype]:
        # we will use any experimental map we can get a hold of
        self.map_manager = self.expedition.maps[etype]['expt']
        return True
    return False

  def probe_at_residue(self, residue, chain_id, resid, resname, struct_type):
    carboxylate_amide_atoms = ["OD1", "OD2", "ND1", "ND2", "OE1", "OE2", "NE1", "NE2", "CD1", "CD2"]
    coordinates = flex.vec3_double([atom.xyz for atom in residue.atoms() \
      if atom.name.strip() in carboxylate_amide_atoms])
    densities = self.map_manager.density_at_sites_cart(coordinates)
    print ("... probed one {res} residue at {timestr}".format(res=resname, timestr=time.asctime()))
    # return ("densities at {}: ".format(res=resname) + " ".join(map(str, list(densities))))
    return ("total density at {res}: ".format(res=resname) + str(sum(densities)))