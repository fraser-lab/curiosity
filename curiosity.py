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
# from scitbx.array_family import flex
from iotbx.reflection_file_reader import any_reflection_file
import time

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
fmodel_file = None
  .type = path
  .help = "Path to a reflection file of calculated structure factors from the model."
experiment = *electron xray neutron
  .type = choice
  .help = "Type of experiment (dictating the scattering table to use for the"
  .help = "calculated map)"
d_min = 3
  .type = float
  .help = "Estimated global resolution in Angstroms"
"""

master_phil = libtbx.phil.parse(master_phil_str)

helpstr = """Curiosity: a tool for exploring crystallographic and cryoEM maps.

My job is to look over your model and find things nearby that aren't modeled
yet, or maybe are modeled wrong. I use the model to look at what's already
there but mainly my analysis is of the map itself. I hope you find me useful!

Some suggestions for best results:
- Big maps == big memory requirements. Especially for EM maps, please run
  phenix.map_box first. Some things won't run at all without this step.
- For finding ions: run phenix.douse before curiosity to make sure you have
  a full set of water candidates to start from. You can choose to set
  keep_input_water=True if you want to preserve positions you've already
  modeled or refined.

Available parameters:
""" + master_phil_str
# TODO: useful helpstring

def validate_params(params):
  import os
  if not params.model_file:
    raise Sorry("A molecular model (.pdb or .cif) is required.")
  if (not params.map_file) and (not params.mtz_file):
    raise Sorry("Please supply a map or structure factors.")
  # if (not params.mtz_file and not (params.difference_map_file and params.calculated_map_file)):
  #   raise Sorry("Please supply either structure factors (.mtz) or calculated and difference map files (.map, .mrc or .ccp4).")
  # if not params.mtz_file:
  #   raise Sorry("Please supply an mtz with f_obs.")
  if (not params.difference_map_file or not params.calculated_map_file):
    raise NotImplementedError("Generation of fmodel and difference maps not yet implemented. Please supply both maps.")
  if params.d_min is None:
    raise Sorry("Please supply an estimated global resolution d_min in Angstroms.")
  print ("Finished validating parameters at {timestr}".format(timestr=time.asctime()))

def run(args):
  if not args or "-h" in args or "--help" in args:
    raise Usage(helpstr)
  # from iotbx import file_reader
  import iotbx.phil
  # import phenix_util
  cmdline = iotbx.phil.process_command_line_with_files(
    args=args,
    master_phil=master_phil,
    )
  if cmdline.unused_args:
    print ("\nEncountered unrecognized parameters:", str(cmdline.unused_args), "\n")
    return
  params = cmdline.work.extract()
  validate_params(params)
  with open("params.out", "w") as outf:
    outf.write(cmdline.work.as_str())
  print ("Finished parsing parameters at {timestr}".format(timestr=time.asctime()))
  # process the inputs and construct a map_model_manager
  from iotbx.data_manager import DataManager
  dm = DataManager()
  mmm = dm.get_map_model_manager(
    model_file=params.model_file, map_files=params.map_file)
  # TODO: probably want to add hydrogens to the model if not present
  expt_mm = mmm.map_managers()[0]
  expt_mm.labels = ['expt']
  mmm.add_map_manager_by_id(expt_mm, '%s_expt' % params.experiment)
  print ("... added experimental map to map_model_manager at {timestr}".format(timestr=time.asctime()))
  # TODO: generate fobs from map data if mtz not provided
  # expt_millers = expt_mm.map_as_fourier_coefficients()
  if params.mtz_file:
    fobs_arrays = any_reflection_file(params.mtz_file).as_miller_arrays()
    for arr in fobs_arrays:
      labels = arr.info().labels
      if 'FC' in labels: break # FIXME this is terrible
    fobs_array = arr.as_amplitude_array()
    print ("... processed mtz file at {timestr}".format(timestr=time.asctime()))
  else:
    fobs_array = None
  # import pdb; pdb.set_trace()
  # TODO: handle properly and select correct arrays
  # TODO: calculate fmodel map if not supplied
  # TODO: calculate fmodel miller array from map if map but not mtz supplied
  # fmodel_obj = any_reflection_file(params.fmodel_file)
  # # import pdb; pdb.set_trace()
  if params.fmodel_file:
    fmodel_miller_array = any_reflection_file(params.fmodel_file).as_miller_arrays()[0]
    print ("... processed fmodel file at {timestr}".format(timestr=time.asctime()))
  else:
    fmodel_miller_array = None
    # # TODO: check for multiple, select the right one
    # fmodel_miller_array = fmodel_miller_array.map_to_asu().customized_copy(
    #     data = flex.double(fmodel_miller_array.data().size(), 1))
    # fmodel_mm = file_reader.any_file(params.calculated_map_file).file_object
  if params.calculated_map_file:
    fmodel_mm = dm.get_map_model_manager(model_file=params.model_file, map_files=params.calculated_map_file).map_manager()
    fmodel_mm.labels = ['fmodel']
    mmm.add_map_manager_by_id(fmodel_mm, '%s_fmodel' % params.experiment)
    print ("... added fmodel map to map_model_manager at {timestr}".format(timestr=time.asctime()))
  else:
    # TODO: calculate difference map if not supplied
    # diff_mm = file_reader.any_file(params.difference_map_file).file_object
    fmodel_mm = None
  if params.difference_map_file:
    diff_mm = dm.get_map_model_manager(model_file=params.model_file, map_files=params.difference_map_file).map_manager()
    diff_mm.labels = ['diff']
    mmm.add_map_manager_by_id(diff_mm, '%s_diff' % params.experiment)
    print ("... added difference map to map_model_manager at {timestr}".format(timestr=time.asctime()))
  else:
    diff_mm = None
    # TODO: repeat for additional supplied maps from other types of experiments (refactor phil too)
    # WISHLIST: fmodel by local resolution in the matching map
  from curiosity_core import Expedition
  curiosity_expedition = Expedition(mmm,
                                    params,
                                    fobs=fobs_array,
                                    fmodel=fmodel_miller_array)
  print ("Finished initializing curiosity expedition at {timestr}".format(timestr=time.asctime()))
  curiosity_expedition.walk()
  print ("Finished expedition analysis at {timestr}".format(timestr=time.asctime()))
  # print ("Discoveries:")
  # print(curiosity_expedition.discoveries)
  log = curiosity_expedition.log_results("{model}_curiosity_results.out".format(model=params.model_file))
  print("Wrote discoveries to {log}".format(log=log))
  summary = curiosity_expedition.log_summary("{model}_curiosity_summary.out".format(model=params.model_file))
  print("Wrote summary to {summary}".format(summary=summary))
  # TODO: write Coot script to look through discoveries, labeled appropriately
  # curiosity_expedition.apply_selected_discoveries(curiosity_expedition.discoveries)
  # print("Applied all discovered modifications to model.")
  # import pdb; pdb.set_trace()
  print ("Analysis complete!")

if __name__=="__main__":
  import sys
  print ("Beginning execution of curiosity.py at {timestr}".format(timestr=time.asctime()))
  run(sys.argv[1:])
