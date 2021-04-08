from __future__ import division

# Curiosity: a tool for exploring electron density or eletrostatic potential maps.
#
# This program accepts as input a model (PDB or mMCIF) and a map (MRC, CCP4 or MTZ).
# If supplying structure factors, curiosity will compute a calculated map and a
# difference map; otherwise these should be supplied as well. The program will look
# for unmodeled features matching particular criteria and produce a script for
# jumping directly to each of these features in Coot. This program expands on its
# predecessor qPTxM, for quantification of post-transcriptoinal modifications, and
# reuses much of its infrastructure. Also like qPTxM, it uses the computational
# crystallographic toolbox (cctbx) for handling maps and models.

import libtbx.phil
from libtbx.utils import Sorry, Usage

master_phil_str = """
model_file = None
  .type = path
  .help = "Path to the model to be analyzed, in pdb or mmcif format"
map_file = None
  .type = path
  .help = "Path to the map into which the model was built, in mrc or ccp4 format"
mtz_file = None
  .type = path
  .help = "Path to an mtz containing structure factors"
difference_map_file = None
  .type = path
  .help = "Path to a difference map to use directly (if not calculating one)."
calculated_map_file = None
  .type = path
  .help = "Path to a calculated map to use directly (if not calculating one)."
experiment = *electron xray neutron
  .type = choice
  .help = "Type of experiment (dictating the scattering table to use for the"
  .help = "calculated map)"
d_min = 3
  .type = float
  .help = "Estimated global resolution in Angstroms"
"""

master_phil = libtbx.phil.parse(master_phil_str)

helpstr = """
Available parameters:
""" + master_phil_str
# TODO: useful helpstring

def validate_params(params):
  import os
  if not params.model_file:
    raise Sorry("A molecular model (.pdb or .cif) is required.")
  if (not params.map_file) and (not params.mtz_file):
    raise Sorry("Please supply a map or structure factors.")
  if (not params.mtz_file and not (params.difference_map_file and params.calculated_map_file)):
    raise Sorry("Please supply either structure factors (.mtz) or calculated and difference map files (.map, .mrc or .ccp4).")
  if params.d_min is None:
    raise Sorry("Please supply an estimated global resolution d_min in Angstroms.")

def run(args):
  if not args or "-h" in args or "--help" in args:
    raise Usage(helpstr)
  from iotbx import file_reader
  import iotbx.phil
  cmdline = iotbx.phil.process_command_line_with_files(
    args=args,
    master_phil=master_phil,
    )
  if cmdline.unused_args:
    print "\nEncountered unrecognized parameters:", str(cmdline.unused_args), "\n"
    return
  params = cmdline.work.extract()
  validate_params(params)
  with open("params.out", "wb") as outf:
    outf.write(cmdline.work.as_str())
  # process the model
  model_in = file_reader.any_file(params.model_file, force_type="pdb")
  model_in.check_file_type("pdb")
  hier_model = model_in.file_object.construct_hierarchy()
  # TODO: use map_model_manager instead
  # TODO: probably want to add hydrogens if not present
  # process the map(s)
  def get_mtz_file_object(mtz_path):
    if not mtz_path: return
    pass #TODO load mtz and produce calc and diff maps
  def get_map_file_object(map_path):
    if not map_path: return
    map_in = file_reader.any_file(map_path, force_type="ccp4_map")
    map_in.check_file_type("ccp4_map")
    return map_in.file_object
  # TODO: use map_model_manager for the maps as well
  # TODO: get a matching fmodel for each experiment type
  # WISHLIST: fmodel by local resolution in the matching map
  # TODO: calculate a difference map for each experimental and fmodel pair
  # TODO: make sure everything makes it into the map_model_manager with accurate
  # and informative labels
  # TODO: run it
  # TODO: write discoveries to file
  # TODO: write Coot script to look through discoveries, labeled appropriately
  print "Analysis complete!"

if __name__=="__main__":
  import sys
  run(sys.argv[1:])
