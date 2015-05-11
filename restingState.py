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

import SEMTools as sem
from nipype.interfaces.ants.registration import Registration
from nipype.interfaces.ants import ApplyTransforms
from nipype.interfaces.freesurfer.preprocess import *
from nipype.interfaces.utility import Function, IdentityInterface, Rename, Select
from nipype.interfaces.utility import Merge as Merger
0import nipype.pipeline.engine as pipe
import numpy

import dataio
import registrationWorkflow
import preprocessWorkflow
import nuisanceWorkflow
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
    nacAtlasLabel = "/Shared/paulsen/Experiments/rsFMRI/rs-fMRI-pilot/ReferenceAtlas/template_nac_labels.nii.gz"
    nacResampleResolution = (2.0, 2.0, 2.0)
    downsampledNACfilename = 'downsampledNACatlas.nii.gz'

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

    # HACK: Remove node from pipeline until Nipype/AFNI file copy issue is resolved
    # preproc.connect(sessions, 'session_id', fmri_DataSink, 'container')
    # END HACK
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
        infields = ['subject_id', 'session_id', 'year', 'day'],
        outfields = ['fmriHdr', 't1_File']
        grabber = dataio.clevelandGrabber(infields, outfields,
                                          base_directory='/Shared/paulsen/Experiments',
                                          experiment=t1_experiment)

    preproc.connect(sessions, 'session_id', grabber, 'session_id')

    formatFMRINode = pipe.Node(interface=Function(function=formatFMRI,
                                                  input_names=['dicomDirectory'],
                                                  output_names=['modality', 'numberOfSlices',
                                                                'numberOfFiles',
                                                                'repetitionTime', 'sliceOrder']),
                               name='formatFMRINode')
    preproc.connect(grabber, 'fmri_dicom_dir', formatFMRINode, 'dicomDirectory')

    to_3D_str = afninodes.to3dstrnode('strCreate')
    to_3D = afninodes.to3dnode('afniTo3D')
    preproc.connect([(formatFMRINode, to_3D_str, [('numberOfSlices', 'slices'),
                                                  ('numberOfFiles', 'volumes'),
                                                  ('repetitionTime', to_3D_str, 'repTime'),
                                                  ('sliceOrder', to_3D_str, 'order')
                                                  ]),
                     (to_3D_str, to_3D, [('funcparams', 'funcparams')])
                    ])
    ##############
    # FMRI space #
    ##############
    preproc.connect(grabber, 'fmri_dicom_dir', to_3D, 'infolder')

    # HACK: Remove node from pipeline until Nipype/AFNI file copy issue is resolved
    # renameTo3D = pipe.Node(Rename(format_string='%(session)s_to3D'), name='renameTo3D')
    # renameTo3D.inputs.keep_ext = True
    # preproc.connect(to_3D, 'out_file', renameTo3D, 'in_file')
    # preproc.connect(sessions, 'session_id', renameTo3D, 'session')
    # preproc.connect(renameTo3D, 'out_file', fmri_DataSink, '@To3D')
    # END HACK

    rgwf = registrationWorkflow.workflow(t1_experiment, outputType, name="registration_wkfl")
    downsampleAtlas = pipe.Node(interface=Function(function=resampleImage,
                                                   input_names=['inputVolume', 'outputVolume', 'resolution'],
                                                   output_names=['outputVolume']),
                                name="downsampleAtlas")
    downsampleAtlas.inputs.inputVolume = nacAtlasFile
    downsampleAtlas.inputs.outputVolume = downsampledNACfilename
    downsampleAtlas.inputs.resolution = [int(x) for x in nacResampleResolution]

    preproc.connect(sessions, 'session_id', rgwf, "transformGrabber.session_id")
    preproc.connect(grabber, 't1_File', rgwf.bFit, 'fixedVolume')
    preproc.connect(grabber, 't1_File', rgwf.warpT1ToFMRI, 'input_image')
    preproc.connect(grabber, 'csfFile', rgwf.warpBABCSegToFMRI, 'input_image')
    preproc.connect(downsampleAtlas, 'outputVolume', rgwf.fmriToNAC_epi, 'reference_image')

    if preprocessOn:
        preprocessing = preprocessWorkflow.workflow(skipCount=6, name="preprocessing_wkfl", outputType=outputType)
        nuisance = nuisanceWorkflow.workflow(nacAtlasLabel, outputType, name="nuisance_wkfl")

        preproc.connect([(to_3D, preprocessing, [('out_file', 'refit.in_file')]),     # 1a
                         (preprocessing, rgwf,  [('merge.out_file', 'tstat.in_file'),  # 7 ### 'mean_file' -> 'in_file'])
                                                 ('automask.out_file', 'tstat.mask_file')])
                         (rgwf, rgwf, [('0_List_forwardTransformNACToFMRI.out', 'antsApplyTransformsCSF.transforms'),
                                       ('afni3DtStat.out_file', 'antsApplyTransformsCSF.reference_image'),
                                       ('0_List_forwardTransformT1ToFMRI.out', 'antsApplyTransformsWHM.transforms'),
                                       ('afni3DtStat.out_file', 'antsApplyTransformsWHM.reference_image')]),
                         (rgwf, nuisance, [('antsApplyTransformWHM.output_image', 'wmMask.input_file'),
                                           ('antsApplyTransformCSF.output_image', 'afni3DmaskAve_csf.mask')]),
                         (nuisance, rgwf, [('csfMask.output_file', 'antsApplyTransformsCSF.input_image')]),
                         (preprocessing, nuisance, [('afni3Dcalc.out_file', 'afni3DmaskAve_wm.in_file'),
                                                    ('afni3Dcalc.out_file', 'afni3Ddeconvolve.in_file')]),
                         (grabber, rgwf,  [('whmFile', 'antsApplyTransformWHM.input_image')]),
                        ])

        if maskGM:
            preproc.connect([(rgwf, rgwf, [('0_List_forwardTransformT1ToFMRI.out', 'antsApplyTransformsGRM.transforms'),
                                           ('afni3DtStat.out_file', 'antsApplyTransformsGRM.reference_image')]),
                             (rgwf, nuisance, [('antsApplyTransformGRM.output_image', 'grmMask.input_file')]),
                             (preprocessing, nuisance, [('afni3Dcalc.out_file', 'afni3DmaskAve_grm.in_file')])
                             (grabber, rgwf,  [('gryFile', 'antsApplyTransformGRM.input_image')])
                            ])
        elif maskWholeBrain:
            nacWholeBrainFile = "/Shared/paulsen/Experiments/rsFMRI/rs-fMRI-pilot/ReferenceAtlas/template_brain.nii.gz"
            node = rgwf.get_node('antsApplyTransformsWholeBrain')
            node.inputs.input_image = nacWholeBrainFile
            preproc.connect([(rgwf, rgwf, [('0_List_forwardTransformT1ToFMRI.out', 'antsApplyTransformsWholeBrain.transforms'),
                                           ('afni3DtStat.out_file', 'antsApplyTransformsWholeBrain.reference_image')]),
                             (rgwf, nuisance, [('antsApplyTransformGRM.output_image', 'wholeBrainMask.input_file')]),
                             (preprocessing, nuisance, [('afni3Dcalc.out_file', 'afni3DmaskAve_whole.in_file')])
                            ])
        preproc.connect(deconvolve, 'out_errts', detrend, 'in_file')  # 13
    # Per Dawei, bandpass after running 3dDetrend
    detrend = afninodes.detrendnode(outputType, 'afni3Ddetrend')
    preproc.connect(detrend, 'out_file', fourier, 'in_file')

    renameMasks = pipe.Node(interface=Rename(format_string='%(label)s_mask'), name='renameMasksAtlas')
    renameMasks.inputs.keep_ext = True

    atlas_DataSink = dataio.datasink(base_directory=preproc.base_dir, container=RESULTS_DIR,
                                     overwrite=REWRITE_DATASINKS, name="atlas_DataSink")


    preproc.connect([(renameMasks, atlas_DataSink,     [('out_file', 'Atlas')]),
                     (downsampleAtlas, atlas_DataSink, [('outputVolume', 'Atlas.@resampled')]),
                     (warpT1ToFMRI, nacToFMRI,         [('output_image', 'reference_image')]),
                    ])

    renameMasks2 = pipe.Node(interface=Rename(format_string='%(session)s_%(label)s_mask'), name='renameMasksFMRI')
    renameMasks2.inputs.keep_ext = True
    preproc.connect(sessions, 'session_id', renameMasks2, 'session')

    seedwkfl = seedWorkflow.workflow()  # TODO

    clipSeedWithVentriclesNode = pipe.Node(interface=Function(function=clipSeedWithVentricles,
                                           input_names=['unclipped_seed_fn', 'fmriBABCSeg_fn', 'desired_out_seed_fn'],
                                           output_names=['clipped_seed_fn']),
                                           name='clipSeedWithVentriclesNode')
    clipSeedWithVentriclesNode.inputs.desired_out_seed_fn = "clipped_seed.nii.gz"

    preproc.connect(rgwf.nacToFMRI, 'output_image', clipSeedWithVentriclesNode, 'unclipped_seed_fn')
    preproc.connect(rgwf.warpBABCSegToFMRI, 'output_image', clipSeedWithVentriclesNode, 'fmriBABCSeg_fn')
    if not maskWhiteMatterFromSeeds:
        preproc.connect(clipSeedWithVentriclesNode, 'clipped_seed_fn', renameMasks2, 'in_file')
    else:
        clipSeedWithWhiteMatterNode = pipe.Node(interface=Function(function=clipSeedWithWhiteMatter,
                                                                   input_names=['seed', 'mask', 'outfile'],
                                                                   output_names=['outfile']),
                                                name='clipSeedWithWhiteMatterNode')
        clipSeedWithWhiteMatterNode.inputs.outfile = 'clipped_wm_seed.nii.gz'
        preproc.connect(rgwf.warpBABCSegToFMRI, 'output_image', clipSeedWithWhiteMatterNode, 'mask')
        preproc.connect(clipSeedWithVentriclesNode, 'clipped_seed_fn', clipSeedWithWhiteMatterNode, 'seed')
        preproc.connect(clipSeedWithWhiteMatterNode, 'outfile', renameMasks2, 'in_file')

    # Labels are iterated over, so we need a seperate datasink to avoid overwriting any preprocessing
    # results when the labels are iterated (e.g. To3d output)
    fmri_label_DataSink = dataio.datasink(os.path.join(preproc.base_dir, RESULTS_DIR), 'EPI'
                                          'fmri_label_DataSink', REWRITE_DATASINKS)
    # '/Shared/paulsen/Experiments/20130417_rsfMRI_Results/EPI'
    preproc.connect(sessions, 'session_id', fmri_label_DataSink, 'container')
    preproc.connect(renameMasks2, 'out_file', fmri_label_DataSink, 'masks')
    preproc.connect(rgwf.fourier, 'out_file', fmri_label_DataSink, 'masks.@bandpass')

    roiMedian = maskavenode('AFNI_1D', 'afni_roiMedian', '-mrange 1 1')
    preproc.connect(renameMasks2, 'out_file', roiMedian, 'mask')
    preproc.connect(rgwf.fourier, 'out_file', roiMedian, 'in_file')

    correlate = afninodes.fimnode('Correlation', 'afni_correlate')
    preproc.connect(roiMedian, 'out_file', correlate, 'ideal_file')
    preproc.connect(rgwf.fourier, 'out_file', correlate, 'in_file')

    regionLogCalc = afninodes.logcalcnode(outputType, 'afni_regionLogCalc')
    preproc.connect(correlate, 'out_file', regionLogCalc, 'in_file_a')

    renameZscore = pipe.Node(interface=Rename(format_string="%(session)s_%(label)s_zscore"), name='renameZscore')
    renameZscore.inputs.keep_ext = True
    preproc.connect(sessions, 'session_id', renameZscore, 'session')
    preproc.connect(regionLogCalc, 'out_file', renameZscore, 'in_file')
    preproc.connect(renameZscore, 'out_file', fmri_label_DataSink, 'zscores')

    # Move z values back into NAC atlas space
    fmriToNAC_label = fmriToNAC_epi.clone(name='fmriToNac_label')
    fmriToNAC_label.inputs.interpolation = 'Linear'
    fmriToNAC_label.inputs.invert_transform_flags = [False, False]
    preproc.connect(downsampleAtlas, 'outputVolume', fmriToNAC_label, 'reference_image')
    preproc.connect(regionLogCalc, 'out_file', fmriToNAC_label, 'input_image')
    preproc.connect(rgwf.reverseTransform, 'out', fmriToNAC_label, 'transforms')

    renameZscore2 = pipe.Node(interface=Rename(format_string="%(session)s_%(label)s_result"), name='renameZscore2')
    renameZscore2.inputs.keep_ext = True
    preproc.connect(sessions, 'session_id', renameZscore2, 'session')
    preproc.connect(fmriToNAC_label, 'output_image', renameZscore2, 'in_file')
    preproc.connect(renameZscore2, 'out_file', atlas_DataSink, 'Atlas.@zscore')

    # Connect seed subworkflow
    seedSubflow = seedWorkflow.workflow(name='seed_wkfl')
    preproc.connect([(downsampleAtlas, seedSubflow,    [('outputVolume', 'afni3Dcalc_seeds.in_file_a')]),
                     (seedSubflow, renameMasks,        [('afni3Dcalc_seeds.out_file', 'in_file'),
                                                        ('selectLabel.out', 'label')]),
                     (seedSubflow, renameMasks2,       [('afni3Dcalc_seeds.out_file', 'in_file'),
                                                        ('selectLabel.out', 'label')]),
                     (seedSubflow, renameZscore,       [('selectLabel.out', 'label')]),
                     (seedSubflow, renameZscore2,      [('selectLabel.out', 'label')]),
                     (seedSubflow, rgwf.nacToFMRI,     [('afni3Dcalc_seeds.out_file', 'input_image')])
                    ])

    preproc.write_graph()
    # preproc.write_hierarchical_dotfile(dotfilename='dave.dot')
    if os.environ['USER'] == 'dmwelch' and False:
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
