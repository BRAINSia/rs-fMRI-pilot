#!/usr/bin/env python
"""
Usage:
  restingState.py [-h | --help]
  restingState.py [options] [-g | -b] (-n NAME | --name NAME) SESSION...

Arguments:
  -n NAME, --name NAME        experiment name, format: 'YYYYMMDD_<experiment>'
  SESSION                     one or more session IDs

Options:
  -h, --help                  show this help message and exit
  -d, --debug                 run in debug mode
  -f, --force                 force DataSink rewriting
  -p, --preprocess            run preprocessing pipeline (not Cleveland)
  -g, --maskGM                global signal regression by masking gray matter (cannot be combined with -b)
  -b, --maskWB                whole brain (cannot be combined with -g)
  -w, --maskSeeds             mask white matter from seeds
  -F FORMAT, --format FORMAT  output format, values: afni, nifti, nifti_gz [default: nifti]

Note: Logs are written to the $PWD/$USER/logs directory

Example:
  restingState.py -n my_new_experiment 0001 0002 0003
  restingState.py -fgF nifti_gz --name new_experiment_w_gray_mask 00001 00002 00003
  restingState.py -bn my_new_experiment_brain_mask --format afni 0001 0002 0003

"""
from nipype import config, logging

# TODO: Modify virtualenv to include formatFMRI.sh
import os
import re
import sys

from nipype.interfaces.freesurfer.preprocess import *
from nipype.interfaces.utility import Function, IdentityInterface, Rename, Select
from nipype.interfaces.afni.preprocess import Copy
import nipype.pipeline.engine as pipe
import numpy

import afninodes
import dataio
from utilities import *
import registrationWorkflow
import preprocessWorkflow
import nuisanceWorkflow
import seedWorkflow


def pipeline(args):
    if args['debug']:
        config.enable_debug_mode()
        logging.update_logging(config)

    # CONSTANTS
    sessionID = args['session']
    outputType = args['format'].upper()
    fOutputType = args['freesurfer']
    preprocessOn = args['preprocess']
    maskGM = args['maskgm']
    maskWholeBrain = args['maskwb']
    maskWhiteMatterFromSeeds = args['maskseeds']
    REWRITE_DATASINKS = args['force']  # Toggle REWRITE_DATASINKS per command line flag
    # print args['name']
    CACHE_DIR = args['name'] + "_Cache"  # Cache directory
    RESULTS_DIR = args['name'] + "_Results"  # Universal datasink directory
    t1_experiment = "20141001_PREDICTHD_long_Results"  #"20130729_PREDICT_Results"
    nacAtlasFile = "/Shared/paulsen/Experiments/rsFMRI/rs-fMRI-pilot/ReferenceAtlas/template_t1.nii.gz"
    nacWholeBrainFile = "/Shared/paulsen/Experiments/rsFMRI/rs-fMRI-pilot/ReferenceAtlas/template_brain.nii.gz"
    nacAtlasLabel = "/Shared/paulsen/Experiments/rsFMRI/rs-fMRI-pilot/ReferenceAtlas/template_nac_labels.nii.gz"
    nacResampleResolution = (2.0, 2.0, 2.0)
    downsampledNACfilename = 'downsampledNACatlas.nii.gz'
    # if preprocessing and nuisance is needed
    args['csf_input'] = nacAtlasFile
    args['wb_input'] = nacWholeBrainFile

    preproc = pipe.Workflow(updatehash=True, name=args['name'])
    userLogDir = os.path.join(os.getcwd(), os.environ['USER'], 'logs')
    if not os.path.isdir(userLogDir):
        os.makedirs(userLogDir, 0777)
    config.update_config({'logging': {'log_directory':userLogDir}})
    logging.update_logging(config)
    preproc.base_dir = os.path.abspath("/Shared/sinapse/CACHE")

    sessions = pipe.Node(interface=IdentityInterface(fields=['session_id']), name='sessionIDs')
    sessions.iterables = ('session_id', sessionID)
    downsampleAtlas = pipe.Node(interface=Function(function=resampleImage,
                                                   input_names=['inputVolume', 'outputVolume', 'resolution'],
                                                   output_names=['outputVolume']),
                                name="downsampleAtlas")
    downsampleAtlas.inputs.inputVolume = nacAtlasFile
    downsampleAtlas.inputs.outputVolume = downsampledNACfilename
    downsampleAtlas.inputs.resolution = [int(x) for x in nacResampleResolution]

    # HACK: Remove node from pipeline until Nipype/AFNI file copy issue is resolved
    # fmri_DataSink = pipe.Node(interface=DataSink(), name="fmri_DataSink")
    # fmri_DataSink.overwrite = REWRITE_DATASINKS
    # Output to: /Shared/paulsen/Experiments/YYYYMMDD_<experiment>_Results/fmri
    # fmri_DataSink.inputs.base_directory = os.path.join(preproc.base_dir, RESULTS_DIR, 'fmri')
    # fmri_DataSink.inputs.substitutions = [('to_3D_out+orig', 'to3D')]
    # fmri_DataSink.inputs.parameterization = False
    #
    # preproc.connect([(sessions, fmri_DataSink, [('session_id', 'container')])])
    # END HACK

    rgwf = registrationWorkflow.workflow(t1_experiment, outputType, name="registration_wkfl")
    preproc.connect([(sessions, rgwf, [('session_id', "inputs.session_id")])])

    detrend = afninodes.detrendnode(outputType, 'afni3Ddetrend')
    # define grabber
    site = "*"
    subject = "*"
    if preprocessOn:
        grabber = dataio.iowaGrabber(t1_experiment, site, subject, maskGM, maskWholeBrain)
        preproc.connect([(sessions, grabber, [('session_id', 'session_id')]),
                         (grabber, rgwf,     [('t1_File', 'inputs.t1')])])

        preprocessing = preprocessWorkflow.prepWorkflow(skipCount=6, outputType=outputType)
        args.pop('name')  # HACK: prevent name conflict with nuisance workflow
        nuisance = nuisanceWorkflow.workflow(outputType=outputType, **args)
        preproc.connect([(grabber, preprocessing,  [('fmri_dicom_dir', 'to_3D.infolder'),
                                                    ('fmri_dicom_dir', 'formatFMRINode.dicomDirectory')]),
                         (grabber, nuisance,       [('whmFile', 'wm.warpWMtoFMRI.input_image')]),
                         (preprocessing, rgwf,     [('merge.out_file', 'inputs.fmri'),  # 7
                                                    ('automask.out_file', 'tstat.mask_file')]),  # *optional*
                         (rgwf, nuisance,          [('outputs.fmri_reference', 'csf.warpCSFtoFMRI.reference_image'),  # CSF
                                                    ('outputs.nac2fmri_list', 'csf.warpCSFtoFMRI.transforms'),
                                                    ('outputs.fmri_reference', 'wm.warpWMtoFMRI.reference_image'),    # WM
                                                    ('outputs.t12fmri_list', 'wm.warpWMtoFMRI.transforms')]),
                        ])
        if maskGM:
            preproc.connect([(grabber, nuisance,       [('gryFile', 'gm.warpGMtoFMRI.input_image')]),
                             (rgwf, nuisance,          [('outputs.fmri_reference', 'gm.warpGMtoFMRI.reference_image'),
                                                        ('outputs.t12fmri_list', 'gm.warpGMtoFMRI.transforms')]),
                             (preprocessing, nuisance, [('calc.out_file', 'gm.afni3DmaskAve_grm.in_file'),
                                                        ('volreg.oned_file', 'afni3Ddeconvolve.stim_file_4')])])
        elif maskWholeBrain:
            preproc.connect([(rgwf, nuisance,          [('outputs.fmri_reference', 'wb.warpCSFtoFMRI.reference_image'),
                                                        ('outputs.nac2fmri_list', 'wb.warpCSFtoFMRI.transforms')]),
                             (preprocessing, nuisance, [('calc.out_file', 'wb.afni3DmaskAve_whole.in_file'),
                                                        ('volreg.oned_file', 'afni3Ddeconvolve.stim_file_4')])])
        else:
            preproc.connect([(preprocessing, nuisance, [('volreg.oned_file', 'afni3Ddeconvolve.stim_file_3')])])

        preproc.connect([(preprocessing, nuisance, [('calc.out_file', 'wm.afni3DmaskAve_wm.in_file'),
                                                    ('calc.out_file', 'csf.afni3DmaskAve_csf.in_file'),
                                                    ('calc.out_file', 'afni3Ddeconvolve.in_file')]),
                         (nuisance, detrend,       [('afni3Ddeconvolve.out_errts', 'in_file')])])  # 13
    else:
        cleveland_grabber = dataio.clevelandGrabber()
        grabber = dataio.autoworkupGrabber(t1_experiment, site, subject)
        converter = pipe.Node(interface=Copy(), name='converter')  # Convert ANALYZE to AFNI

        preproc.connect([(sessions, grabber,            [('session_id', 'session_id')]),
                         (grabber, rgwf,                [('t1_File', 'inputs.t1')]),
                         (sessions, cleveland_grabber,  [('session_id', 'session_id')]),
                         (cleveland_grabber, converter, [('fmriHdr', 'in_file')]),
                         (converter, rgwf,              [('out_file', 'inputs.fmri')]),
                         (converter, detrend,           [('out_file', 'in_file')]),  # in fMRI_space
                        ])

    t1_wf = registrationWorkflow.t1Workflow()
    babc_wf = registrationWorkflow.babcWorkflow()
    # HACK: No EPI
    # epi_wf = registrationWorkflow.epiWorkflow()
    lb_wf = registrationWorkflow.labelWorkflow()
    seed_wf = registrationWorkflow.seedWorkflow()
    fourier = afninodes.fouriernode(outputType, 'fourier') # Fourier is the last NIFTI file format in the AFNI pipeline

    preproc.connect([(detrend, fourier,       [('out_file', 'in_file')]), # Per Dawei, bandpass after running 3dDetrend
                     (grabber, t1_wf,         [('t1_File', 'warpT1toFMRI.input_image')]),
                     (rgwf, t1_wf,            [('outputs.fmri_reference', 'warpT1toFMRI.reference_image'),  # T1
                                               ('outputs.t12fmri_list', 'warpT1toFMRI.transforms')]),
                     (grabber, babc_wf,       [('csfFile', 'warpBABCtoFMRI.input_image')]),
                     (rgwf, babc_wf,          [('outputs.fmri_reference', 'warpBABCtoFMRI.reference_image'),  # Labels
                                               ('outputs.t12fmri_list', 'warpBABCtoFMRI.transforms')]),
                     # HACK: No EPI
                     # (downsampleAtlas, epi_wf, [('outputVolume', 'warpEPItoNAC.reference_image')]),
                     # (rgwf, epi_wf,         [('outputs.fmri2nac_list', 'warpEPItoNAC.transforms')]),
                     # (fourier, epi_wf,      [('out_file', 'warpEPItoNAC.input_image')]),
                     # END HACK
                     (downsampleAtlas, lb_wf, [('outputVolume', 'warpLabeltoNAC.reference_image')]),
                     (rgwf, lb_wf,            [('outputs.fmri2nac_list', 'warpLabeltoNAC.transforms')]),
                     (t1_wf, seed_wf,         [('warpT1toFMRI.output_image', 'warpSeedtoFMRI.reference_image')]),
                     (rgwf, seed_wf,          [('outputs.nac2fmri_list', 'warpSeedtoFMRI.transforms')]),
                     ])

    renameMasks = pipe.Node(interface=Rename(format_string='%(label)s_mask'), name='renameMasksAtlas')
    renameMasks.inputs.keep_ext = True

    atlas_DataSink = dataio.datasink(base_directory=preproc.base_dir, container=RESULTS_DIR,
                                     overwrite=REWRITE_DATASINKS, name="atlas_DataSink")

    preproc.connect([(renameMasks, atlas_DataSink,     [('out_file', 'Atlas')]),
                     (downsampleAtlas, atlas_DataSink, [('outputVolume', 'Atlas.@resampled')]),
                    ])

    renameMasks2 = pipe.Node(interface=Rename(format_string='%(session)s_%(label)s_mask'), name='renameMasksFMRI')
    renameMasks2.inputs.keep_ext = True
    preproc.connect(sessions, 'session_id', renameMasks2, 'session')

    clipSeedWithVentriclesNode = pipe.Node(interface=Function(function=clipSeedWithVentricles,
                                           input_names=['unclipped_seed_fn', 'fmriBABCSeg_fn', 'desired_out_seed_fn'],
                                           output_names=['clipped_seed_fn']),
                                           name='clipSeedWithVentriclesNode')
    clipSeedWithVentriclesNode.inputs.desired_out_seed_fn = "clipped_seed.nii.gz"

    preproc.connect(seed_wf, 'warpSeedtoFMRI.output_image', clipSeedWithVentriclesNode, 'unclipped_seed_fn')
    preproc.connect(babc_wf, 'warpBABCtoFMRI.output_image', clipSeedWithVentriclesNode, 'fmriBABCSeg_fn')
    if not maskWhiteMatterFromSeeds:
        preproc.connect(clipSeedWithVentriclesNode, 'clipped_seed_fn', renameMasks2, 'in_file')
    else:
        clipSeedWithWhiteMatterNode = pipe.Node(interface=Function(function=clipSeedWithWhiteMatter,
                                                                   input_names=['seed', 'mask', 'outfile'],
                                                                   output_names=['outfile']),
                                                name='clipSeedWithWhiteMatterNode')
        clipSeedWithWhiteMatterNode.inputs.outfile = 'clipped_wm_seed.nii.gz'
        preproc.connect(babc_wf, 'warpBABCtoFMRI.output_image', clipSeedWithWhiteMatterNode, 'mask')
        preproc.connect(clipSeedWithVentriclesNode, 'clipped_seed_fn', clipSeedWithWhiteMatterNode, 'seed')
        preproc.connect(clipSeedWithWhiteMatterNode, 'outfile', renameMasks2, 'in_file')
    # Labels are iterated over, so we need a seperate datasink to avoid overwriting any preprocessing
    # results when the labels are iterated (e.g. To3d output)
    # Write out to: /Shared/sinapse/CACHE/YYYYMMDD_<experiment>_Results/EPI
    fmri_label_DataSink = dataio.datasink(os.path.join(preproc.base_dir, RESULTS_DIR), container='EPI',
                                          name='fmri_label_DataSink', overwrite=REWRITE_DATASINKS)
    preproc.connect(sessions, 'session_id', fmri_label_DataSink, 'container')
    preproc.connect(renameMasks2, 'out_file', fmri_label_DataSink, 'masks')
    preproc.connect(fourier,'out_file', fmri_label_DataSink, 'masks.@bandpass')

    roiMedian = afninodes.maskavenode('AFNI_1D', 'afni_roiMedian', '-mrange 1 1')
    preproc.connect(renameMasks2, 'out_file', roiMedian, 'mask')
    preproc.connect(fourier, 'out_file', roiMedian, 'in_file')

    correlate = afninodes.fimnode('Correlation', 'afni_correlate')
    preproc.connect(roiMedian, 'out_file', correlate, 'ideal_file')
    preproc.connect(fourier, 'out_file', correlate, 'in_file')

    regionLogCalc = afninodes.logcalcnode(outputType, 'afni_regionLogCalc')
    preproc.connect(correlate, 'out_file', regionLogCalc, 'in_file_a')

    renameZscore = pipe.Node(interface=Rename(format_string="%(session)s_%(label)s_zscore"), name='renameZscore')
    renameZscore.inputs.keep_ext = True
    preproc.connect(sessions, 'session_id', renameZscore, 'session')
    preproc.connect(regionLogCalc, 'out_file', renameZscore, 'in_file')
    preproc.connect(renameZscore, 'out_file', fmri_label_DataSink, 'zscores')

    # Move z values back into NAC atlas space
    # preproc.connect(downsampleAtlas, 'outputVolume', lb_wf, 'warpLabeltoNAC.reference_image')
    preproc.connect(regionLogCalc, 'out_file', lb_wf, 'warpLabeltoNAC.input_image')

    renameZscore2 = pipe.Node(interface=Rename(format_string="%(session)s_%(label)s_result"), name='renameZscore2')
    renameZscore2.inputs.keep_ext = True
    preproc.connect(sessions, 'session_id', renameZscore2, 'session')
    preproc.connect(lb_wf, 'warpLabeltoNAC.output_image', renameZscore2, 'in_file')
    preproc.connect(renameZscore2, 'out_file', atlas_DataSink, 'Atlas.@zscore')
    # Connect seed subworkflow
    seedSubflow = seedWorkflow.workflow(outputType='NIFTI_GZ', name='seed_wkfl')
    preproc.connect([(downsampleAtlas, seedSubflow,    [('outputVolume', 'afni3Dcalc_seeds.in_file_a')]),
                     (seedSubflow, renameMasks,        [('afni3Dcalc_seeds.out_file', 'in_file'),
                                                        ('selectLabel.out', 'label')]),
                     (seedSubflow, renameMasks2,       [('selectLabel.out', 'label')]),
                     (seedSubflow, renameZscore,       [('selectLabel.out', 'label')]),
                     (seedSubflow, renameZscore2,      [('selectLabel.out', 'label')]),
                     (seedSubflow, seed_wf,               [('afni3Dcalc_seeds.out_file', 'warpSeedtoFMRI.input_image')])
                    ])
    if args['debug']:
        preproc.write_hierarchical_dotfile(dotfilename='debug_hier.dot')
        preproc.write_graph(dotfilename='debug_orig.dot', graph2use="orig",  #='hierarchical',
                            format='png', simple_form=False)
        # preproc.run()
        # Run restingState on the all threads
        # Setup environment for CPU load balancing of ITK based programs.
        # --------
        # import multiprocessing
        # total_CPUS = 11  # multiprocessing.cpu_count()
        # preproc.run(plugin='MultiProc', plugin_args={'n_proc': total_CPUS})
        # --------
        # Run restingState on the local cluster
        # preproc.run(plugin='SGE', plugin_args={'template': os.path.join(os.getcwd(), 'ENV/bin/activate'),
        #                                        'qsub_args': '-S /bin/bash -cwd'})
    else:
        preproc.write_graph(dotfilename=args['name'] + '.dot')
        import multiprocessing
        # Setup environment for CPU load balancing of ITK based programs.
        total_CPUS = multiprocessing.cpu_count()
        preproc.run(plugin='MultiProc', plugin_args={'n_proc': total_CPUS})
    return 0

if __name__ == '__main__':
    from docopt import docopt
    args = docopt(__doc__, version='1.0')
    keys = args.keys()
    for key in keys:
        # Return dictionary with lowercase keys, without leading "-"'s
        value = args.pop(key)
        key = key.lstrip('-')
        args[key.lower()] = value
        if key == 'name':  # check that experiment naming convention is kept
            if re.match(r"20[0-9]{6}_", value) is None:
                raise RuntimeError("Experiment name must begin with 'YYYYMMDD_'")
    freesurferOutputTypes = {"nifti_gz": "niigz",
                             "afni": "afni",
                             "nifti": "nii"}
    args['freesurfer'] = freesurferOutputTypes[args['format']]
    outvalue = pipeline(args)
    sys.exit(outvalue)
