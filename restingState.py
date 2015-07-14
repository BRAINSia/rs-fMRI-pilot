#!/usr/bin/env python
"""
Usage:
  restingState.py [-h | --help]
  restingState.py [options] [-g | -b] (-n NAME | --name NAME) (-s SEEDS | --seeds SEEDS) SESSION...

Arguments:
  -n NAME, --name NAME        experiment name, format: 'YYYYMMDD_<experiment>'
  -s SEEDS, --seeds SEEDS     seed fiducial file in .fcsv format
  SESSION                     one or more session IDs

Options:
  -h, --help                  show this help message and exit
  -b, --maskWB                whole brain (requires -p)
  -d, --debug                 run in debug mode
  -f, --force                 force DataSink rewriting
  -F FORMAT, --format FORMAT  output format, values: afni, nifti, nifti_gz [default: nifti]
  -g, --maskGM                global signal regression by masking gray matter (requires -p)
  -p, --preprocess            run preprocessing pipeline (not Cleveland)
  -P, --plot                  write out the workflow graphs to /tmp and exit
  -w, --maskSeeds             mask white matter from seeds

Note: Logs are written to the /tmp directory

Example:
  restingState.py -n my_new_experiment 0001 0002 0003
  restingState.py -fgF nifti_gz --name new_experiment_w_gray_mask 00001 00002 00003
  restingState.py -bn my_new_experiment_brain_mask --format afni 0001 0002 0003

"""
from nipype import config, logging

import os
import re
import sys

# from nipype.interfaces.freesurfer.preprocess import *
from nipype.interfaces.utility import Function, IdentityInterface, Rename
from afni.preprocess import Copy
import nipype.pipeline.engine as pipe

import afninodes
import dataio
from utilities import resampleImage, clipSeedWithVentricles, clipSeedWithWhiteMatter
import registrationWorkflow
import preprocessWorkflow
import nuisanceWorkflow
import seedWorkflow


def makeSupportDir(name, suffix):
    """
    Create /tmp/name.suffix, if it d.n.e. This prevents permission collisions with different users
    running the same pipeline
    """
    retval = os.path.abspath(os.path.join(os.path.sep, 'tmp', name + '.' + suffix))
    if not os.path.isdir(retval):
        os.makedirs(retval, 0777)
    return retval


def pipeline(args):
    if args['debug']:
        config.enable_debug_mode()
    config.update_config({'logging': {'log_directory':makeSupportDir(args['name'], "logs")}})
    logging.update_logging(config)

    # CONSTANTS
    sessionID = args['session']
    outputType = args['format'].upper()
    fOutputType = args['freesurfer']
    preprocessOn = args['preprocess']
    maskGM = args['maskgm']
    maskWholeBrain = args['maskwb']
    maskWhiteMatterFromSeeds = args['maskseeds']
    # print args['name']
    t1_experiment = "20141001_PREDICTHD_long_Results"  #"20130729_PREDICT_Results"
    atlasFile = "/Shared/paulsen/Experiments/rsFMRI/rs-fMRI-pilot/ReferenceAtlas/template_t1.nii.gz"
    wholeBrainFile = "/Shared/paulsen/Experiments/rsFMRI/rs-fMRI-pilot/ReferenceAtlas/template_brain.nii.gz"
    atlasLabel = "/Shared/paulsen/Experiments/rsFMRI/rs-fMRI-pilot/ReferenceAtlas/template_nac_labels.nii.gz"
    resampleResolution = (2.0, 2.0, 2.0)
    downsampledfilename = 'downsampled_atlas.nii.gz'
    # if preprocessing and nuisance is needed
    # args['csf_input'] = atlasFile
    # args['wb_input'] = wholeBrainFile

    master = pipe.Workflow(name=args['name'] + "_CACHE")
    master.base_dir = os.path.abspath("/Shared/sinapse/CACHE")

    sessions = pipe.Node(interface=IdentityInterface(fields=['session_id']), name='sessionIDs')
    sessions.iterables = ('session_id', sessionID)
    downsampleAtlas = pipe.Node(interface=Function(function=resampleImage,
                                                   input_names=['inputVolume', 'outputVolume', 'resolution'],
                                                   output_names=['outputVolume']),
                                name="downsampleAtlas")
    downsampleAtlas.inputs.inputVolume = atlasFile
    downsampleAtlas.inputs.outputVolume = downsampledfilename
    downsampleAtlas.inputs.resolution = [int(x) for x in resampleResolution]

    # HACK: Remove node from pipeline until Nipype/AFNI file copy issue is resolved
    # fmri_DataSink = pipe.Node(interface=DataSink(), name="fmri_DataSink")
    # fmri_DataSink.overwrite = REWRITE_DATASINKS
    # Output to: /Shared/paulsen/Experiments/YYYYMMDD_<experiment>_Results/fmri
    # fmri_DataSink.inputs.base_directory = os.path.join(master.base_dir, RESULTS_DIR, 'fmri')
    # fmri_DataSink.inputs.substitutions = [('to_3D_out+orig', 'to3D')]
    # fmri_DataSink.inputs.parameterization = False
    #
    # master.connect([(sessions, fmri_DataSink, [('session_id', 'container')])])
    # END HACK

    registration = registrationWorkflow.workflow(t1_experiment, outputType, name="registration_wkfl")
    master.connect([(sessions, registration, [('session_id', "inputs.session_id")])])

    detrend = afninodes.detrendnode(outputType, 'afni3Ddetrend')
    # define grabber
    site = "*"
    subject = "*"
    if preprocessOn:
        grabber = dataio.iowaGrabber(t1_experiment, site, subject, maskGM, maskWholeBrain)
        master.connect([(sessions, grabber, [('session_id', 'session_id')]),
                         (grabber, registration,     [('t1_File', 'inputs.t1')])])
        # Why isn't preprocessWorkflow.workflow() used instead? It would avoid most of the nuisance connections here...
        preprocessing = preprocessWorkflow.prepWorkflow(skipCount=6, outputType=outputType)
        name = args.pop('name')  # HACK: prevent name conflict with nuisance workflow
        nuisance = nuisanceWorkflow.workflow(outputType=outputType, **args)
        args['name'] = name  # END HACK
        master.connect([(grabber, preprocessing,      [('fmri_dicom_dir', 'to_3D.infolder'),
                                                       ('fmri_dicom_dir', 'formatFMRINode.dicomDirectory')]),
                        (grabber, nuisance,           [('whmFile', 'wm.warpWMtoFMRI.input_image')]),
                        (preprocessing, registration, [('merge.out_file', 'inputs.fmri'),  # 7
                                                       ('automask.out_file', 'tstat.mask_file')]),  # *optional*
                        (registration, nuisance,      [('outputs.fmri_reference', 'csf.warpCSFtoFMRI.reference_image'),  # CSF
                                                       ('outputs.nac2fmri_list', 'csf.warpCSFtoFMRI.transforms'),
                                                       ('outputs.fmri_reference', 'wm.warpWMtoFMRI.reference_image'),    # WM
                                                       ('outputs.t12fmri_list', 'wm.warpWMtoFMRI.transforms')]),
                        ])
        warpCSFtoFMRInode = nuisance.get_node('csf').get_node('warpCSFtoFMRI')
        warpCSFtoFMRInode.input_file = atlasFile
        if maskGM:
            master.connect([(grabber, nuisance,       [('gryFile', 'gm.warpGMtoFMRI.input_image')]),
                            (registration, nuisance,  [('outputs.fmri_reference', 'gm.warpGMtoFMRI.reference_image'),
                                                       ('outputs.t12fmri_list', 'gm.warpGMtoFMRI.transforms')]),
                            (preprocessing, nuisance, [('calc.out_file', 'gm.afni3DmaskAve_grm.in_file'),
                                                       ('volreg.oned_file', 'afni3Ddeconvolve.stim_file_4')])])
        elif maskWholeBrain:
            master.connect([(registration, nuisance,  [('outputs.fmri_reference', 'wb.warpBraintoFMRI.reference_image'),
                                                       ('outputs.nac2fmri_list', 'wb.warpBraintoFMRI.transforms')]),
                            (preprocessing, nuisance, [('calc.out_file', 'wb.afni3DmaskAve_whole.in_file'),
                                                       ('volreg.oned_file', 'afni3Ddeconvolve.stim_file_4')])])
            warpBraintoFMRInode = nuisance.get_node('wb').get_node('warpBraintoFMRI')
            warpBraintoFMRInode.input_mask = wholeBrainFile
        else:
            master.connect([(preprocessing, nuisance, [('volreg.oned_file', 'afni3Ddeconvolve.stim_file_3')])])

        master.connect([(preprocessing, nuisance, [('calc.out_file', 'wm.afni3DmaskAve_wm.in_file'),
                                                   ('calc.out_file', 'csf.afni3DmaskAve_csf.in_file'),
                                                   ('calc.out_file', 'afni3Ddeconvolve.in_file')]),
                        (nuisance, detrend,       [('afni3Ddeconvolve.out_errts', 'in_file')])])  # 13
    else:
        cleveland_grabber = dataio.clevelandGrabber()
        grabber = dataio.autoworkupGrabber(t1_experiment, site, subject)
        converter = pipe.Node(interface=Copy(), name='converter')  # Convert ANALYZE to AFNI

        master.connect([(sessions, grabber,            [('session_id', 'session_id')]),
                         (grabber, registration,        [('t1_File', 'inputs.t1')]),
                         (sessions, cleveland_grabber,  [('session_id', 'session_id')]),
                         (cleveland_grabber, converter, [('fmriHdr', 'in_file')]),
                         (converter, registration,      [('out_file', 'inputs.fmri')]),
                         (converter, detrend,           [('out_file', 'in_file')]),  # in fMRI_space
                        ])

    t1_wf = registrationWorkflow.t1Workflow()
    babc_wf = registrationWorkflow.babcWorkflow()
    # HACK: No EPI
    # epi_wf = registrationWorkflow.epiWorkflow()
    lb_wf = registrationWorkflow.labelWorkflow()
    seed_wf = registrationWorkflow.seedWorkflow()
    bandpass = afninodes.fouriernode(outputType, 'fourier') # Fourier is the last NIFTI file format in the AFNI pipeline

    master.connect([(detrend, bandpass,       [('out_file', 'in_file')]), # Per Dawei, bandpass after running 3dDetrend
                     (grabber, t1_wf,         [('t1_File', 'warpT1toFMRI.input_image')]),
                     (registration, t1_wf,    [('outputs.fmri_reference', 'warpT1toFMRI.reference_image'),  # T1
                                               ('outputs.t12fmri_list', 'warpT1toFMRI.transforms')]),
                     (grabber, babc_wf,       [('csfFile', 'warpBABCtoFMRI.input_image')]),
                     (registration, babc_wf,  [('outputs.fmri_reference', 'warpBABCtoFMRI.reference_image'),  # Labels
                                               ('outputs.t12fmri_list', 'warpBABCtoFMRI.transforms')]),
                     # HACK: No EPI
                     # (downsampleAtlas, epi_wf, [('outputVolume', 'warpEPItoNAC.reference_image')]),
                     # (registration, epi_wf,    [('outputs.fmri2nac_list', 'warpEPItoNAC.transforms')]),
                     # (bandpass, epi_wf,         [('out_file', 'warpEPItoNAC.input_image')]),
                     # END HACK
                     (downsampleAtlas, lb_wf, [('outputVolume', 'warpLabeltoNAC.reference_image')]),
                     (registration, lb_wf,    [('outputs.fmri2nac_list', 'warpLabeltoNAC.transforms')]),
                     (t1_wf, seed_wf,         [('warpT1toFMRI.output_image', 'warpSeedtoFMRI.reference_image')]),
                     (registration, seed_wf,  [('outputs.nac2fmri_list', 'warpSeedtoFMRI.transforms')]),
                     ])

    renameMasks = pipe.Node(interface=Rename(format_string='%(label)s_mask'), name='renameMasksAtlas')
    renameMasks.inputs.keep_ext = True
    atlas_DataSink = dataio.atlasSink(base_directory=master.base_dir, **args)
    master.connect([(renameMasks, atlas_DataSink,     [('out_file', 'Atlas')]),
                    (downsampleAtlas, atlas_DataSink, [('outputVolume', 'Atlas.@resampled')]),
                    ])

    renameMasks2 = pipe.Node(interface=Rename(format_string='%(session)s_%(label)s_mask'), name='renameMasksFMRI')
    renameMasks2.inputs.keep_ext = True
    master.connect(sessions, 'session_id', renameMasks2, 'session')

    clipSeedWithVentriclesNode = pipe.Node(interface=Function(function=clipSeedWithVentricles,
                                           input_names=['seed', 'label', 'outfile'],
                                           output_names=['clipped_seed_fn']),
                                           name='clipSeedWithVentriclesNode')
    clipSeedWithVentriclesNode.inputs.outfile = "clipped_seed.nii.gz"

    master.connect(seed_wf, 'warpSeedtoFMRI.output_image', clipSeedWithVentriclesNode, 'seed')
    master.connect(babc_wf, 'warpBABCtoFMRI.output_image', clipSeedWithVentriclesNode, 'label')
    if not maskWhiteMatterFromSeeds:
        master.connect(clipSeedWithVentriclesNode, 'clipped_seed_fn', renameMasks2, 'in_file')
    else:
        clipSeedWithWhiteMatterNode = pipe.Node(interface=Function(function=clipSeedWithWhiteMatter,
                                                                   input_names=['seed', 'mask', 'outfile'],
                                                                   output_names=['outfile']),
                                                name='clipSeedWithWhiteMatterNode')
        clipSeedWithWhiteMatterNode.inputs.outfile = 'clipped_wm_seed.nii.gz'
        master.connect(babc_wf, 'warpBABCtoFMRI.output_image', clipSeedWithWhiteMatterNode, 'mask')
        master.connect(clipSeedWithVentriclesNode, 'clipped_seed_fn', clipSeedWithWhiteMatterNode, 'seed')
        master.connect(clipSeedWithWhiteMatterNode, 'outfile', renameMasks2, 'in_file')
    # Labels are iterated over, so we need a seperate datasink to avoid overwriting any preprocessing
    # results when the labels are iterated (e.g. To3d output)
    # Write out to: /Shared/sinapse/CACHE/YYYYMMDD_<experiment>_Results/<SESSION>
    fmri_label_DataSink = dataio.fmriSink(master.base_dir, **args)
    master.connect(sessions, 'session_id', fmri_label_DataSink, 'container')
    master.connect(renameMasks2, 'out_file', fmri_label_DataSink, 'masks')
    master.connect(bandpass,'out_file', fmri_label_DataSink, 'masks.@bandpass')

    roiMedian = afninodes.maskavenode('AFNI_1D', 'afni_roiMedian', '-mrange 1 1')
    master.connect(renameMasks2, 'out_file', roiMedian, 'mask')
    master.connect(bandpass, 'out_file', roiMedian, 'in_file')

    correlate = afninodes.fimnode('Correlation', 'afni_correlate')
    master.connect(roiMedian, 'out_file', correlate, 'ideal_file')
    master.connect(bandpass, 'out_file', correlate, 'in_file')

    regionLogCalc = afninodes.logcalcnode(outputType, 'afni_regionLogCalc')
    master.connect(correlate, 'out_file', regionLogCalc, 'in_file_a')

    renameZscore = pipe.Node(interface=Rename(format_string="%(session)s_%(label)s_zscore"), name='renameZscore')
    renameZscore.inputs.keep_ext = True
    master.connect(sessions, 'session_id', renameZscore, 'session')
    master.connect(regionLogCalc, 'out_file', renameZscore, 'in_file')
    master.connect(renameZscore, 'out_file', fmri_label_DataSink, 'zscores')

    # Move z values back into NAC atlas space
    # master.connect(downsampleAtlas, 'outputVolume', lb_wf, 'warpLabeltoNAC.reference_image')
    master.connect(regionLogCalc, 'out_file', lb_wf, 'warpLabeltoNAC.input_image')

    renameZscore2 = pipe.Node(interface=Rename(format_string="%(session)s_%(label)s_result"), name='renameZscore2')
    renameZscore2.inputs.keep_ext = True
    master.connect(sessions, 'session_id', renameZscore2, 'session')
    master.connect(lb_wf, 'warpLabeltoNAC.output_image', renameZscore2, 'in_file')
    master.connect(renameZscore2, 'out_file', atlas_DataSink, 'Atlas.@zscore')
    # Connect seed subworkflow
    seedSubflow = seedWorkflow.workflow(args['seeds'], outputType='NIFTI_GZ', name='seed_wkfl')
    master.connect([(downsampleAtlas, seedSubflow,    [('outputVolume', 'afni3Dcalc_seeds.in_file_a')]),
                     (seedSubflow, renameMasks,        [('afni3Dcalc_seeds.out_file', 'in_file'),
                                                        ('selectLabel.out', 'label')]),
                     (seedSubflow, renameMasks2,       [('selectLabel.out', 'label')]),
                     (seedSubflow, renameZscore,       [('selectLabel.out', 'label')]),
                     (seedSubflow, renameZscore2,      [('selectLabel.out', 'label')]),
                     (seedSubflow, seed_wf,               [('afni3Dcalc_seeds.out_file', 'warpSeedtoFMRI.input_image')])
                    ])
    imageDir = makeSupportDir(args['name'], "images")
    if args['plot']:
        registration.write_graph(dotfilename=os.path.join(imageDir, 'register.dot'), graph2use='orig', format='png',
                                 simple_form=False)
        if preprocessOn:
            preprocessing.write_graph(dotfilename=os.path.join(imageDir, 'preprocess.dot'), graph2use='orig', format='png',
                                      simple_form=False)
            nuisance.write_graph(dotfilename=os.path.join(imageDir, 'nuisance.dot'), graph2use='orig', format='png',
                                 simple_form=False)
        seedSubflow.write_graph(dotfilename=os.path.join(imageDir, 'seed.dot'), graph2use='orig', format='png',
                                 simple_form=False)
        master.write_graph(dotfilename=os.path.join(imageDir, 'master.dot'), graph2use="orig", format='png', simple_form=False)
    elif args['debug']:
        try:
            master.run()  #updatehash=True)
            # Run restingState on the all threads
            # Setup environment for CPU load balancing of ITK based programs.
            # --------
            # import multiprocessing
            # total_CPUS = 10  # multiprocessing.cpu_count()
            # master.run(plugin='MultiProc', plugin_args={'n_proc': total_CPUS})  #, updatehash=True)
            # --------
            # Run restingState on the local cluster
            # master.run(plugin='SGE', plugin_args={'template': os.path.join(os.getcwd(), 'ENV/bin/activate'),
            #                                        'qsub_args': '-S /bin/bash -cwd'})  #, updatehash=True)
        except:
            pass
        master.name = "master"  # HACK: Bug in Graphviz for nodes beginning with numbers
        master.write_graph(dotfilename=os.path.join(imageDir, 'debug_hier.dot'), graph2use="colored", format='png')
        master.write_graph(dotfilename=os.path.join(imageDir, 'debug_orig.dot'), graph2use="flat", format='png')
    else:
        #HACK commented out for now
        #master.write_graph(dotfilename='images/' + args['name'] + '.dot')
        import multiprocessing
        # Setup environment for CPU load balancing of ITK based programs.
        total_CPUS = multiprocessing.cpu_count()
        master.run(plugin='MultiProc', plugin_args={'n_proc': total_CPUS})  #, updatehash=True)
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
            if value is None:  # no flags were given
                raise RuntimeWarning("Must specify a value for --name or provide the --help flag")
                sys.exit(1)
            if re.match(r"20[0-9]{6}_", value) is None:
                raise RuntimeError("Experiment name must begin with 'YYYYMMDD_'")
    freesurferOutputTypes = {"nifti_gz": "niigz",
                             "afni": "afni",
                             "nifti": "nii"}
    args['freesurfer'] = freesurferOutputTypes[args['format']]
    # Docopt doesn't recognize 'optional or', e.g. "[-p [-g | -b]]"
    if args['maskgm'] or args['maskwb']:
        assert args['preprocess'], "-g and -b flags must accompany -p"
    outvalue = pipeline(args)
    sys.exit(outvalue)
