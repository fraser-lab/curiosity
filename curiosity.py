from codecs import ignore_errors

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
from libtbx import easy_run
import os, shutil
import time

master_phil_str = """
model_file = None
  .type = path
  .help = "Path to the model to be analyzed, in pdb or mmcif format"
map_file = None
  .type = path
  .help = "Path to the map into which the model was built, in mrc or ccp4 format"
refls_file = None
  .type = path
  .help = "Path to an mtz or sf-cif containing structure factors"
refls_labels = None
  .type = str
  .help = "Column labels for experimental measurements in sf-cif file"
rfree_labels = None
  .type = str
  .help = "Column labels for R-free labels in sf-cif file"
expt_mtz_labels = None
  .type = str
  .help = "Column labels for 2mFo-DFc map in mtz file"
diff_mtz_labels = None
  .type = str
  .help = "Column labels for mFo-DFc map in mtz file"
difference_map_file = None
  .type = path
  .help = "Path to a difference map to use directly (if not calculating one)."
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
  if (not params.map_file) and (not params.refls_file):
    raise Sorry("Please supply a map or structure factors.")
  if (not params.refls_file and not (params.difference_map_file) and not params.experiment == "electron"):
    raise Sorry("Please supply either structure factors (.mtz or .cif) or experimental and difference map files (.map, .mrc or .ccp4).")
  if params.d_min is None:
    raise Sorry("Please supply an estimated global resolution d_min in Angstroms.")
  print ("Finished validating parameters at {timestr}".format(timestr=time.asctime()))

def run(args):
  if not args or "-h" in args or "--help" in args:
    raise Usage(helpstr)
  import iotbx.phil
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
  # generate any missing maps from reflection tables if supplied
  if not params.map_file or not params.difference_map_file:
    if params.experiment in ["xray", "neutron"]:
      # TODO: someday we'll have to account for electron diffraction in addition to electron microscopy
      assert params.refls_file, "Sorry: cannot proceed without either all required maps or structure factors"
      if os.path.splitext(params.refls_file)[-1] == ".mtz" and params.expt_mtz_labels and params.diff_mtz_labels:
        mtz_path = params.refls_file
        command = "phenix.mtz2map {model} {mtz} labels={labels}".format(
          model=params.model_file, mtz=mtz_path, labels=params.expt_mtz_labels)
        easy_run.fully_buffered(command=command).raise_if_errors()
        tmp = "{model}_1.ccp4".format(model=os.path.splitext(params.model_file)[0])
        expt_map_path = "{model}_expt.ccp4".format(model=os.path.splitext(params.model_file)[0])
        os.rename(tmp, expt_map_path)
        command = "phenix.mtz2map {model} {mtz} labels={labels}".format(
          model=params.model_file, mtz=mtz_path, labels=params.diff_mtz_labels)
        easy_run.fully_buffered(command=command).raise_if_errors()
        tmp = "{model}_1.ccp4".format(model=os.path.splitext(params.model_file)[0])
        diff_map_path = "{model}_diff.ccp4".format(model=os.path.splitext(params.model_file)[0])
        os.rename(tmp, diff_map_path)
      else:
        # generate maps from structure factors
        scattering_table_lookup = {"xray":"n_gaussian","electron":"electron","neutron":"neutron"}
        command = "phenix.maps {model} {refls} maps.input.reflection_data.labels={refls_labels} maps.input.reflection_data.r_free_flags.label={rfree_labels} scattering_table={stable}".format(model=params.model_file, refls=params.refls_file, refls_labels=params.refls_labels, rfree_labels=params.rfree_labels, stable=scattering_table_lookup[params.experiment])
        easy_run.fully_buffered(command=command).raise_if_errors()
        mtz_path = "{model}_map_coeffs.mtz".format(model=os.path.splitext(params.model_file)[0])
        # now write out the maps
        command = "phenix.mtz2map {model} {mtz}".format(model=params.model_file, mtz=mtz_path)
        easy_run.fully_buffered(command=command).raise_if_errors()
        expt_map_path = "{model}_map_coeffs_2mFo-DFc.ccp4".format(model=os.path.splitext(params.model_file)[0])
        diff_map_path = "{model}_map_coeffs_mFo-DFc.ccp4".format(model=os.path.splitext(params.model_file)[0])
      assert(os.path.exists(expt_map_path)), "Sorry: could not generate 2mFo-DFc map from reflections file"
      assert(os.path.exists(diff_map_path)), "Sorry: could not generate mFo-DFc map from reflections file"
    else: # electron map -- we can generate the difference map in real space
      assert params.map_file, "Sorry: cannot proceed without experimental map"
      # generate boxed maps for both original and difference maps -- must match
      ext = os.path.splitext(params.model_file)[-1]
      shutil.copyfile(params.model_file, "expt{ext}".format(ext=ext))
      command = "phenix.map_box expt{ext} {map}".format(map=params.map_file, ext=ext)
      easy_run.fully_buffered(command=command).raise_if_errors()
      expt_map_path = "expt_box.ccp4"
      if params.difference_map_file:
        unboxed_diff_map_path = params.difference_map_file
      else:
        command = "phenix.real_space_diff_map {model} {map} resolution={res}".format(
          model=params.model_file, map=expt_map_path, res=params.d_min, ext=ext)
        easy_run.fully_buffered(command=command).raise_if_errors()
        unboxed_diff_map_path = "map_model_difference_1.ccp4"
      shutil.copyfile(params.model_file, "diff{ext}".format(ext=ext))
      command = "phenix.map_box diff{ext} {map}".format(map=unboxed_diff_map_path, ext=ext)
      easy_run.fully_buffered(command=command).raise_if_errors()
      diff_map_path = "diff_box.ccp4"
      os.remove("expt{ext}".format(ext=ext)) # model copies only used to force generation of boxed maps with associated names
      os.remove("diff{ext}".format(ext=ext))
  else:
    expt_map_path = params.map_file
    diff_map_path = params.difference_map_file
  # load experimental map
  mmm = dm.get_map_model_manager(
    model_file=params.model_file, map_files=expt_map_path, ignore_symmetry_conflicts=True)
  # TODO: probably want to add hydrogens to the model if not present
  expt_mm = mmm.map_managers()[0]
  expt_mm.labels = ['expt']
  mmm.add_map_manager_by_id(expt_mm, '%s_expt' % params.experiment)
  print ("... added experimental map to map_model_manager at {timestr}".format(timestr=time.asctime()))
  # load difference map
  diff_mm = dm.get_map_model_manager(model_file=params.model_file, map_files=diff_map_path, ignore_symmetry_conflicts=True).map_manager()
  diff_mm.labels = ['diff']
  mmm.add_map_manager_by_id(diff_mm, '%s_diff' % params.experiment)
  print ("... added difference map to map_model_manager at {timestr}".format(timestr=time.asctime()))
  # TODO: repeat for additional supplied maps from other types of experiments (refactor phil too)
  # WISHLIST: fmodel by local resolution in the matching map
  from curiosity_core import Expedition
  curiosity_expedition = Expedition(mmm, params)
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
  start = time.time()
  print ("Beginning execution of curiosity.py at {timestr}".format(timestr=time.asctime()))
  run(sys.argv[1:])
  end = time.time()
  duration = end - start
  seconds = duration % 60
  minutes = (duration // 60) % 3600
  hours = (duration // 3600) % (24*3600)
  days = duration // (24*3600)
  timestr = "{:.2f} seconds".format(seconds)
  if minutes: timestr = "{m} minutes, ".format(m=minutes) + timestr
  if hours: timestr = "{h} hours, ".format(h=hours) + timestr
  if days: timestr = "{d} days, ".format(d=days) + timestr
  print("Program execution took {timestr}.".format(timestr=timestr))
