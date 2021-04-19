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
  if params.mtz_file:
    raise NotImplementedError("Processing of mtz files not yet implemented.")
  if (not params.difference_map_file or not params.calculated_map_file):
    raise NotImplementedError("Generation of fmodel and difference maps not yet implemented. Please supply both maps.")
  if params.d_min is None:
    raise Sorry("Please supply an estimated global resolution d_min in Angstroms.")

def run(args):
  if not args or "-h" in args or "--help" in args:
    raise Usage(helpstr)
  from iotbx import file_reader
  import iotbx.phil
  import phenix_util
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
  # process the inputs and construct a map_model_manager
  from iotbx.data_manager import DataManager
  dm = DataManager()
  mmm = dm.get_map_model_manager(model_file=params.model_file, map_files=params.map_file)
  # TODO: probably want to add hydrogens to the model if not present
  expt_mm = mmm.map_managers()[0]
  expt_mm.labels = ['expt']
  mmm.add_map_manager_by_id(expt_mm, '%s_expt' % params.experiment)
  # TODO: calculate fmodel map if not supplied
  fmodel_mm = file_reader.any_file(params.calculated_map_file).file_object
  fmodel_mm.labels = ['fmodel']
  mmm.add_map_manager_by_id(fmodel_mm, '%s_fmodel' % params.experiment)
  # TODO: calculate difference map if not supplied
  diff_mm = file_reader.any_file(params.difference_map_file).file_object
  diff_mm.labels = ['diff']
  mmm.add_map_manager_by_id(diff_mm, '%s_diff' % params.experiment)
  # TODO: repeat for additional supplied maps from other types of experiments (refactor phil too)
  # WISHLIST: fmodel by local resolution in the matching map
  from curiosity_core import Expedition
  curiosity_expedition = Expedition(mmm, params)
  curiosity_expedition.walk()
  # TODO: write discoveries to file
  # TODO: write Coot script to look through discoveries, labeled appropriately
  print "Analysis complete!"

if __name__=="__main__":
  import sys
  run(sys.argv[1:])
