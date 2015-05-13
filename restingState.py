#!/usr/bin/env python
"""
Usage:
  restingState.py [-h | --help]
  restingState.py [-f | --force] [-p] [-g | -b] [-w] [-F FORMAT | --format=FORMAT] (-n NAME| --name=NAME) SESSION...

Arguments:
  -n NAME, --name=NAME        experiment name
  SESSION                     one or more session IDs

Optional:
  -h, --help                  show this help message and exit
  -f, --force                 force DataSink rewriting
  -p                          run preprocessing pipeline (not Cleveland)
  -g                          global signal regression by masking gray matter
  -b                          whole brain
  -w                          mask white matter from seeds
  -F FORMAT, --format=FORMAT  output format, values: afni, nifti, nifti_gz [default: nifti]

Example:
  restingState.py -n my_new_experiment 0001 0002 0003
  restingState.py -fgF nifti_gz --name new_experiment_w_gray_mask 00001 00002 00003
  restingState.py -bn my_new_experiment_brain_mask --format afni 0001 0002 0003

"""
# TODO: Modify virtualenv to include formatFMRI.sh
import os
import sys

from nipype.interfaces.freesurfer.preprocess import *
from nipype.interfaces.utility import Function, IdentityInterface, Rename, Select
import nipype.pipeline.engine as pipe
import numpy

import afninodes
import dataio
import registrationWorkflow
import preprocessWorkflow
import seedWorkflow

from utilities import *


def pipeline(args):
    # CONSTANTS
    sessionID = args['session']
    outputType = args['format'].upper()
    fOutputType = args['freesurfer']
    preprocessOn = args['p']
    maskGM = args['g']
    maskWholeBrain = args['b']
    maskWhiteMatterFromSeeds = args['w']
    REWRITE_DATASINKS = args['force']  # Toggle REWRITE_DATASINKS per command line flag
    # print args['name']
    CACHE_DIR = "workflow_" + args['name']  # Cache directory
    RESULTS_DIR = args['name'] + "_Results"  # Universal datasink directory
    t1_experiment = "20130729_PREDICT_Results"
    nacAtlasFile = "/Shared/paulsen/Experiments/rsFMRI/rs-fMRI-pilot/ReferenceAtlas/template_t1.nii.gz"
    nacWholeBrainFile = "/Shared/paulsen/Experiments/rsFMRI/rs-fMRI-pilot/ReferenceAtlas/template_brain.nii.gz"
    nacAtlasLabel = "/Shared/paulsen/Experiments/rsFMRI/rs-fMRI-pilot/ReferenceAtlas/template_nac_labels.nii.gz"
    nacResampleResolution = (2.0, 2.0, 2.0)
    downsampledNACfilename = 'downsampledNACatlas.nii.gz'
    # if preprocessing and nuisance is needed
    args['csf_input'] = nacAtlasFile
    args['wb_input'] = nacWholeBrainFile

    preproc = pipe.Workflow(updatehash=True, name=CACHE_DIR)
    preproc.base_dir = os.getcwd()

    # HACK: Remove node from pipeline until Nipype/AFNI file copy issue is resolved
    # fmri_DataSink = pipe.Node(interface=DataSink(), name="fmri_DataSink")
    # fmri_DataSink.overwrite = REWRITE_DATASINKS
    # fmri_DataSink.inputs.base_directory = os.path.join(preproc.base_dir, RESULTS_DIR, 'fmri')
    # '/Shared/paulsen/Experiments/20130417_rsfMRI_Results'
    # fmri_DataSink.inputs.substitutions = [('to_3D_out+orig', 'to3D')]
    # fmri_DataSink.inputs.parameterization = False
    # END HACK

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
    # preproc.connect(sessions, 'session_id', fmri_DataSink, 'container')
    # END HACK
    rgwf = registrationWorkflow.workflow(t1_experiment, outputType, name="registration_wkfl")
    preproc.connect([(sessions,        rgwf, [('session_id', "transformGrabber.session_id")]),
                     (downsampleAtlas, rgwf, [('outputVolume', 'warpEPItoNAC.reference_image')]),
                     ])
    # define grabber
    if preprocessOn:
        site = "*"
        infields = ['session_id']
        tissuecls = 'Experiments/{experiment}/{site}/*/%s/TissueClassify/%s.nii.gz'.format(site=site, experiment=t1_experiment)
        posterior = 'Experiments/{experiment}/{site}/*/%s/ACCUMULATED_POSTERIORS/POSTERIOR_%s_TOTAL.nii.gz'.format(
            site=site, experiment=t1_experiment)
        field_template = dict(fmri_dicom_dir='MRx/{site}/*/%s/%s/%s/*'.format(site=site),
                              csfFile=tissuecls,
                              whmFile=posterior,
                              t1_File=tissuecls)
        template_args = dict(fmri_dicom_dir=[['session_id', 'ANONRAW', 'FMRI_RestingStateConnectivity']],
                             csfFile=[['session_id', 'fixed_brainlabels_seg']],
                             whmFile=[['session_id', 'WM']],
                             t1_File=[['session_id', 't1_average_BRAINSABC']])
        if maskGM:
            field_template['gryFile'] = posterior
            template_args['gryFile'] = [['session_id', 'GM']]
            # For cerebrum gray matter ONLY:
            # field_template['gryFile'] = tissuecls
            # template_args['gryFile'] = [['session_id', 'POSTERIOR_SURFGM']]
        elif maskWholeBrain:
            pass  # No need for grabber, we're using NAC-atlas file
        grabber = dataio.iowaGrabber(infields, field_template, template_args, base_directory="/Shared/paulsen")
    else:
        infields = ['subject_id', 'session_id', 'year', 'day']
        outfields = ['fmriHdr', 't1_File']
        grabber = dataio.clevelandGrabber(infields, outfields,
                                          base_directory='/Shared/paulsen/Experiments',
                                          experiment=t1_experiment)
    preproc.connect(sessions, 'session_id', grabber, 'session_id')

    ##############
    # FMRI space #
    ##############
    # HACK: Remove node from pipeline until Nipype/AFNI file copy issue is resolved
    # renameTo3D = pipe.Node(Rename(format_string='%(session)s_to3D'), name='renameTo3D')
    # renameTo3D.inputs.keep_ext = True
    # preproc.connect(to_3D, 'out_file', renameTo3D, 'in_file')
    # preproc.connect(sessions, 'session_id', renameTo3D, 'session')
    # preproc.connect(renameTo3D, 'out_file', fmri_DataSink, '@To3D')
    # END HACK
    detrend = afninodes.detrendnode(outputType, 'afni3Ddetrend')
    fourier = afninodes.fouriernode(outputType, 'fourier') # Fourier is the last NIFTI file format in the AFNI pipeline
    preproc.add_nodes([detrend, fourier])  #DEBUG
    preproc.write_graph(dotfilename='stage1_h.dot', graph2use='hierarchical', format='png', simple_form=False)  #DEBUG
    if preprocessOn:
        args.pop('name')
        preprocessing = preprocessWorkflow.workflow(skipCount=6, outputType=outputType, name="preprocessing_wkfl", **args)
        preproc.connect([(grabber, preprocessing,  [('fmri_dicom_dir', 'prep.to_3D.infolder'),
                                                    ('fmri_dicom_dir', 'prep.formatFMRINode.dicomDirectory')]),
                         (preprocessing, rgwf,     [('prep.merge.out_file', 'tstat.in_file'),  # 7
                                                    ('prep.automask.out_file', 'tstat.mask_file')]),
                                                    ])
        preproc.write_graph(dotfilename='stage2a_h.dot', graph2use='hierarchical', format='png', simple_form=False)  #DEBUG
        preproc.connect([(rgwf, preprocessing,     [('tstat.out_file', 'nuisance.csf.warpCSFtoFMRI.reference_image'),  # CSF
                                                    ('nactoFMRI_list.out', 'nuisance.csf.warpCSFtoFMRI.transforms'),
                                                    ('tstat.out_file', 'nuisance.wm.warpWMtoFMRI.reference_image'),    # WM
                                                    ('t1toFMRI_list.out', 'nuisance.wm.warpWMtoFMRI.transforms'),
                                                    # ('tstat.out_file', 'nuisance.t1.warpCSFtoFMRI.reference_image'),  # T1
                                                    # ('t1toFMRI_list.out', 'nuisance.t1.warpCSFtoFMRI.transforms'),
                                                    # ('tstat.out_file', 'nuisance.babc.warpCSFtoFMRI.reference_image'),  # Labels
                                                    # ('t1toFMRI_list.out', 'nuisance.babc.warpCSFtoFMRI.transforms')
                                                    ]),

                        ])
        preproc.write_graph(dotfilename='stage2_h.dot', graph2use='hierarchical', format='png', simple_form=False)  #DEBUG
        if maskGM:
            preproc.connect([(grabber, preprocessing,  [('gryFile', 'nuisance.gm.warpGMtoFMRI.input_image')]), ])
            preproc.write_graph(dotfilename='stage3a_h.dot', graph2use='hierarchical', format='png', simple_form=False)  #DEBUG
            preproc.connect([(rgwf,    preprocessing,  [('tstat.out_file', 'nuisance.gm.warpGMtoFMRI.reference_image'),
                                                   ('t1toFMRI_list.out', 'nuisance.gm.warpGMtoFMRI.transforms')]),
                            ])
            preproc.write_graph(dotfilename='stage3b_h.dot', graph2use='hierarchical', format='png', simple_form=False)  #DEBUG
        elif maskWholeBrain:
            preproc.connect([(rgwf, nuisance, [('tstat.out_file', 'wb.warpCSFtoFMRI.reference_image'),
                                               ('nactoFMRI_list.out', 'wb.warpCSFtoFMRI.transforms')]),
                            ])
            preproc.write_graph(dotfilename='stage3b_h.dot', graph2use='hierarchical', format='png', simple_form=False)  #DEBUG
        preproc.connect(preprocessing, 'nuisance.afni3Ddeconvolve.out_errts', detrend, 'in_file')  # 13

    preproc.connect([(grabber, rgwf,           [('t1_File', 'brainsFit.fixedVolume'),
                                                ('t1_File', 'warpT1toFMRI.input_image'),
                                                ('csfFile', 'warpBABCtoFMRI.input_image')]),
                     (grabber, nuisance,       [('whmFile', 'wm.warpWMtoFMRI.input_image')]),
                     (detrend, fourier,        [('out_file', 'in_file')]), # Per Dawei, bandpass after running 3dDetrend
                     (fourier, rgwf,           [('out_file', 'warpEPItoNAC.input_image')]),
                     ])

    preproc.write_graph(dotfilename='stage3_h.dot', graph2use='hierarchical', format='png', simple_form=False)  #DEBUG

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

    # TODO: seedwkfl = seedWorkflow.workflow()

    clipSeedWithVentriclesNode = pipe.Node(interface=Function(function=clipSeedWithVentricles,
                                           input_names=['unclipped_seed_fn', 'fmriBABCSeg_fn', 'desired_out_seed_fn'],
                                           output_names=['clipped_seed_fn']),
                                           name='clipSeedWithVentriclesNode')
    clipSeedWithVentriclesNode.inputs.desired_out_seed_fn = "clipped_seed.nii.gz"

    preproc.connect(rgwf, 'warpSeedtoFMRI.output_image', clipSeedWithVentriclesNode, 'unclipped_seed_fn')
    preproc.connect(rgwf, 'warpBABCtoFMRI.output_image', clipSeedWithVentriclesNode, 'fmriBABCSeg_fn')
    if not maskWhiteMatterFromSeeds:
        preproc.connect(clipSeedWithVentriclesNode, 'clipped_seed_fn', renameMasks2, 'in_file')
    else:
        clipSeedWithWhiteMatterNode = pipe.Node(interface=Function(function=clipSeedWithWhiteMatter,
                                                                   input_names=['seed', 'mask', 'outfile'],
                                                                   output_names=['outfile']),
                                                name='clipSeedWithWhiteMatterNode')
        clipSeedWithWhiteMatterNode.inputs.outfile = 'clipped_wm_seed.nii.gz'
        preproc.connect(rgwf, 'warpBABCtoFMRI.output_image', clipSeedWithWhiteMatterNode, 'mask')
        preproc.connect(clipSeedWithVentriclesNode, 'clipped_seed_fn', clipSeedWithWhiteMatterNode, 'seed')
        preproc.connect(clipSeedWithWhiteMatterNode, 'outfile', renameMasks2, 'in_file')
    preproc.write_graph(dotfilename='stage4_h.dot', graph2use='hierarchical', format='png', simple_form=False)  #DEBUG
    # Labels are iterated over, so we need a seperate datasink to avoid overwriting any preprocessing
    # results when the labels are iterated (e.g. To3d output)
    fmri_label_DataSink = dataio.datasink(os.path.join(preproc.base_dir, RESULTS_DIR), container='EPI',
                                          name='fmri_label_DataSink', overwrite=REWRITE_DATASINKS)
    # '/Shared/paulsen/Experiments/20130417_rsfMRI_Results/EPI'
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
    preproc.connect(downsampleAtlas, 'outputVolume', rgwf, 'warpLabeltoNAC.reference_image')
    preproc.connect(regionLogCalc, 'out_file', rgwf, 'warpLabeltoNAC.input_image')

    renameZscore2 = pipe.Node(interface=Rename(format_string="%(session)s_%(label)s_result"), name='renameZscore2')
    renameZscore2.inputs.keep_ext = True
    preproc.connect(sessions, 'session_id', renameZscore2, 'session')
    preproc.connect(rgwf, 'warpLabeltoNAC.output_image', renameZscore2, 'in_file')
    preproc.connect(renameZscore2, 'out_file', atlas_DataSink, 'Atlas.@zscore')
    # Connect seed subworkflow
    seedSubflow = seedWorkflow.workflow(outputType='NIFTI_GZ', name='seed_wkfl')
    preproc.connect([(downsampleAtlas, seedSubflow,    [('outputVolume', 'afni3Dcalc_seeds.in_file_a')]),
                     (seedSubflow, renameMasks,        [('afni3Dcalc_seeds.out_file', 'in_file'),
                                                        ('selectLabel.out', 'label')]),
                     (seedSubflow, renameMasks2,       [('selectLabel.out', 'label')]),
                     (seedSubflow, renameZscore,       [('selectLabel.out', 'label')]),
                     (seedSubflow, renameZscore2,      [('selectLabel.out', 'label')]),
                     (seedSubflow, rgwf,               [('afni3Dcalc_seeds.out_file', 'warpSeedtoFMRI.input_image')])
                    ])

    # preproc.write_graph()
    # preproc.write_hierarchical_dotfile(dotfilename='dave.dot')
    if os.environ['USER'] == 'dmwelch' and True:
        preproc.write_graph(dotfilename='preprocessing.dot', graph2use='hierarchical', format='png', simple_form=False)
        return 0
        # Run restingState on the local cluster
        preproc.run(plugin='SGE', plugin_args={'template': os.path.join(os.getcwd(), 'ENV/bin/activate'),
                                               'qsub_args': '-S /bin/bash -cwd'})
    else:
        import multiprocessing
        # Setup environment for CPU load balancing of ITK based programs.
        total_CPUS = multiprocessing.cpu_count()
        preproc.run(plugin='MultiProc', plugin_args={'n_proc': total_CPUS})

if __name__ == '__main__':
    from docopt import docopt
    args = docopt(__doc__, version='1.0')
    keys = args.keys()
    for key in keys:
        # Return dictionary with lowercase keys, without leading "-"'s
        value = args.pop(key)
        key = key.lstrip('-')
        args[key.lower()] = value
    freesurferOutputTypes = {"nifti_gz": "niigz",
                             "afni": "afni",
                             "nifti": "nii"}
    args['freesurfer'] = freesurferOutputTypes[args['format']]
    outvalue = pipeline(args)
    sys.exit(outvalue)
