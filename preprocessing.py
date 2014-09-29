#!/usr/bin/env python
"""
Usage:
  preprocessing.py [-h | --help]
  preprocessing.py [-f | --force] [-g | -b] [-w] [-F FORMAT | --format=FORMAT] (-n NAME| --name=NAME) SESSION...

Arguments:
  -n NAME, --name=NAME        experiment name
  SESSION                     one or more session IDs

Optional:
  -h, --help                  show this help message and exit
  -f, --force                 force DataSink rewriting
  -g                          global signal regression by masking gray matter
  -b                          whole brain
  -w                          mask white matter from seeds
  -F FORMAT, --format=FORMAT  output format, values: afni, nifti, nifti_gz [default: nifti]

Example:
  preprocessing.py -n my_new_experiment 0001 0002 0003
  preprocessing.py -fgF nifti_gz --name new_experiment_w_gray_mask 00001 00002 00003
  preprocessing.py -bn my_new_experiment_brain_mask --format afni 0001 0002 0003

"""
# TODO: Modify virtualenv to include formatFMRI.sh
import os
import sys

import SEMTools as sem
from nipype.interfaces.ants.registration import Registration
from nipype.interfaces.ants import ApplyTransforms
from nipype.interfaces.afni.preprocess import *
from nipype.interfaces.freesurfer.preprocess import *
from nipype.interfaces.io import DataSink, DataGrabber
from nipype.interfaces.utility import Function, IdentityInterface, Rename, Select
from nipype.interfaces.utility import Merge as Merger
import nipype.pipeline.engine as pipe
import numpy

from utilities import *


def pipeline(args):
    sessionID = args['session']
    outputType = args['format'].upper()
    fOutputType = args['freesurfer']
    maskGM = args['g']
    maskWholeBrain = args['b']
    maskWhiteMatterFromSeeds = args['w']
    print args['name']
    CACHE_DIR = "workflow_" + args['name']  # Cache directory
    RESULTS_DIR = args['name'] + "_Results"  # Universal datasink directory
    preproc = pipe.Workflow(updatehash=True, name=CACHE_DIR)
    preproc.base_dir = os.getcwd()
    REWRITE_DATASINKS = args['force']  # Toggle REWRITE_DATASINKS per command line flag

    fmri_DataSink = pipe.Node(interface=DataSink(), name="fmri_DataSink")
    fmri_DataSink.overwrite = REWRITE_DATASINKS
    # HACK: Remove node from pipeline until Nipype/AFNI file copy issue is resolved
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
    outfields = ['fmri_dicom_dir', 't1_File', 'csfFile', 'whmFile']
    if maskGM:
        outfields.append('gryFile')
    elif maskWholeBrain:
        pass  # No need for grabber, we're using NAC-atlas file
    grabber = pipe.Node(interface=DataGrabber(infields=['session_id'],
                                              outfields=outfields), name='dataGrabber')
    grabber.inputs.base_directory = '/Shared/paulsen'  # TODO:  Replace with /Shared/paulsen
    grabber.inputs.template = '*'
    site = "*"
    experiment = '20130729_PREDICT_Results'
    fmriRegex = 'MRx/{site}/*/%s/%s/%s/*'.format(site=site)
    probRegex = 'Experiments/{experiment}/{site}/*/%s/%s/%s'.format(site=site, experiment=experiment)
    field_template = dict(fmri_dicom_dir=fmriRegex,
                          csfFile=probRegex,
                          whmFile=probRegex)
    template_args = dict(fmri_dicom_dir=[['session_id', 'ANONRAW', 'FMRI_RestingStateConnectivity']],
                         csfFile=[['session_id', 'TissueClassify', 'fixed_brainlabels_seg.nii.gz']],
                         whmFile=[['session_id', 'ACCUMULATED_POSTERIORS', 'POSTERIOR_WM_TOTAL.nii.gz']])
    field_template['t1_File'] = probRegex
    template_args['t1_File'] = [['session_id', 'TissueClassify', 't1_average_BRAINSABC.nii.gz']]
    grabber.inputs.field_template = field_template
    grabber.inputs.template_args = template_args
    preproc.connect(sessions, 'session_id', grabber, 'session_id')
    # CONSTANTS
    nacAtlasFile = "/Shared/paulsen/Experiments/rsFMRI/rs-fMRI-pilot/ReferenceAtlas/template_t1.nii.gz"
    nacAtlasLabel = "/Shared/paulsen/Experiments/rsFMRI/rs-fMRI-pilot/ReferenceAtlas/template_nac_labels.nii.gz"
    nacResampleResolution = (2.0, 2.0, 2.0)
    downsampledNACfilename = 'downsampledNACatlas.nii.gz'
    # WORKFLOW
    transformGrabber = pipe.Node(interface=DataGrabber(infields=['session_id'],
                                                       outfields=['atlasToSessionTransform',
                                                                  'sessionToAtlasTransform']),
                                 name='transformGrabber')
    transformGrabber.inputs.base_directory = '/Shared/paulsen/Experiments/{experiment}/SubjectToAtlasWarped'.format(
        experiment=experiment)
    transformGrabber.inputs.template = '*'
    transformRegex = '%s/AtlasToSubject_%sComposite.h5'
    transformGrabber.inputs.field_template = dict(atlasToSessionTransform=transformRegex,
                                                  sessionToAtlasTransform=transformRegex)
    transformGrabber.inputs.template_args = dict(atlasToSessionTransform=[['session_id', '']],
                                                 sessionToAtlasTransform=[['session_id', 'Inverse']])
    preproc.connect(sessions, 'session_id', transformGrabber, 'session_id')

    # NOTE: antsApplyTransforms takes transforms in REVERSE order!!!
    # Concatenate transforms: NAC -> T1 -> fMRI
    forwardTransform = pipe.Node(Merger(2), name='0_List_forwardTransformNACToFMRI')
    preproc.connect(transformGrabber, 'atlasToSessionTransform', forwardTransform, 'in2')

    # Concatenate transforms: NAC <- T1 <- fMRI
    reverseTransform = pipe.Node(Merger(2), name='0_List_reverseTransformFMRIToNAC')
    preproc.connect(transformGrabber, 'sessionToAtlasTransform', reverseTransform, 'in1')

    formatFMRINode = pipe.Node(interface=Function(function=formatFMRI,
                                                  input_names=['dicomDirectory'],
                                                  output_names=['modality', 'numberOfSlices',
                                                                'numberOfFiles',
                                                                'repetitionTime', 'sliceOrder']),
                               name='formatFMRINode')
    preproc.connect(grabber, 'fmri_dicom_dir', formatFMRINode, 'dicomDirectory')

    to_3D_str = pipe.Node(interface=Function(function=strCreate,
                                             input_names=['time', 'slices', 'volumes', 'repTime', 'order'],
                                             output_names=['out_string']),
                          name='strCreate')
    to_3D_str.inputs.time = 'zt'
    preproc.connect(formatFMRINode, 'numberOfSlices', to_3D_str, 'slices')  # 1
    preproc.connect(formatFMRINode, 'numberOfFiles', to_3D_str, 'volumes')
    preproc.connect(formatFMRINode, 'repetitionTime', to_3D_str, 'repTime')
    preproc.connect(formatFMRINode, 'sliceOrder', to_3D_str, 'order')

    ##############
    # FMRI space #
    ##############

    to_3D = pipe.Node(interface=To3D(), name='afniTo3D')
    to_3D.inputs.outputtype = 'AFNI'
    to_3D.inputs.datum = 'short'
    to_3D.inputs.filetype = 'epan'
    to_3D.inputs.prefix = 'to_3D_out'
    preproc.connect(grabber, 'fmri_dicom_dir', to_3D, 'infolder')
    preproc.connect(to_3D_str, 'out_string', to_3D, 'funcparams')

    # HACK: Remove node from pipeline until Nipype/AFNI file copy issue is resolved
    # renameTo3D = pipe.Node(Rename(format_string='%(session)s_to3D'), name='renameTo3D')
    # renameTo3D.inputs.keep_ext = True
    # preproc.connect(to_3D, 'out_file', renameTo3D, 'in_file')
    # preproc.connect(sessions, 'session_id', renameTo3D, 'session')
    # preproc.connect(renameTo3D, 'out_file', fmri_DataSink, '@To3D')
    # END HACK

    # 1a
    refit = pipe.Node(interface=Refit(), name='afni3Drefit')
    refit.inputs.outputtype = 'AFNI'
    refit.inputs.deoblique = True
    preproc.connect(to_3D, 'out_file', refit, 'in_file')

    # 2
    skipCount = 6

    def strToIntMinusOne(string):
        return int(string) - 1

    despike = pipe.Node(interface=Despike(), name='afni3Ddespike')
    despike.inputs.outputtype = outputType
    # Since we want to remove the first four from the output,
    # we should ignore them for the despike analysis, right?
    despike.inputs.start = skipCount
    despike.inputs.ignore = skipCount
    preproc.connect(refit, 'out_file', despike, 'in_file')
    preproc.connect(formatFMRINode, ('numberOfFiles', strToIntMinusOne), despike, 'end')

    # 3
    volreg = pipe.Node(interface=Volreg(), name='afni3DvolReg')
    volreg.inputs.outputtype = outputType
    volreg.inputs.timeshift = False  # 0
    volreg.inputs.zpad = 3
    volreg.inputs.interp = 'cubic'
    volreg.inputs.maxite = 50
    volreg.inputs.thresh = 0.001
    volreg.inputs.rot_thresh = 0.001
    volreg.inputs.delta = 0.1
    volreg.inputs.final = 'Fourier'
    volreg.inputs.twopass = True
    volreg.inputs.twodup = True
    volreg.inputs.coarse = [2, 2]
    volreg.inputs.coarserot = True
    volreg.inputs.base = 9
    volreg.inputs.oned_file = 'volReg.1D'
    preproc.connect(despike, 'out_file', volreg, 'in_file')  # 3

    # 4
    zpad = pipe.Node(interface=Zeropad(), name='afniZeropad')
    zpad.inputs.plane = 'IS'
    zpad.inputs.numberOfPlanes = 44
    zpad.inputs.is_mm = False
    preproc.connect(volreg, 'out_file', zpad, 'in_file')  # 4

    # 5
    merge = pipe.Node(interface=Merge(), name='afni3Dmerge')
    merge.inputs.outputtype = outputType
    merge.inputs.blurfwhm = 6
    merge.inputs.doall = True
    merge.inputs.args = '-1noneg -1clip 100'  # TODO
    # TODO: implement
    # merge.inputs.onenoneg = True
    # merge.inputs.oneclip = 100
    # 5 ### NOTE: ONLY ONE FILE TO in_files. JV changed blurr back to 6 on 5/14/2013.
    preproc.connect(zpad, 'out_file', merge, 'in_files')

    # 6
    automask = pipe.Node(interface=Automask(), name='afni3Dautomask')
    automask.inputs.outputtype = outputType
    automask.inputs.dilate = 1
    preproc.connect(merge, 'out_file', automask, 'in_file')  # 6

    # 7
    tstat = pipe.Node(interface=TStat(), name='afni3DtStat')
    tstat.inputs.outputtype = outputType
    tstat.inputs.args = '-mean'  # TODO
    preproc.connect(merge, 'out_file', tstat, 'in_file')  # 7 ### 'mean_file' -> 'in_file'
    preproc.connect(automask, 'out_file', tstat, 'mask_file')

    calc = pipe.Node(interface=Calc(letters=['a', 'b']), name='afni3Dcalc')
    calc.inputs.outputtype = outputType
    calc.inputs.expr = "a * b"
    preproc.connect(merge, 'out_file', calc, 'in_file_a')
    preproc.connect(automask, 'out_file', calc, 'in_file_b')

    # 8
    fourier = pipe.Node(interface=Fourier(), name='afni3Dfourier')
    fourier.inputs.outputtype = outputType
    fourier.inputs.highpass = 0.011
    fourier.inputs.lowpass = 0.1
    # fourier.inputs.args = '-retrend' # removed on 9/20/13
    fourier.inputs.outputtype = outputType

    # 9
    bFit = pipe.Node(interface=sem.BRAINSFit(), name='brainsFit')
    bFit.inputs.initializeTransformMode = 'useCenterOfHeadAlign'
    bFit.inputs.maskProcessingMode = 'ROIAUTO'
    bFit.inputs.ROIAutoDilateSize = 10.0
    bFit.inputs.useRigid = True
    bFit.inputs.costMetric = 'MMI'  # (default)
    bFit.inputs.numberOfSamples = 100000  # (default)
    bFit.inputs.outputTransform = True
    ##########################################################################
    # The registration of the T1 -> fMRI is good _most_ of the time, but sometimes it hits a local minimum during optimization. #
    # To prevent this, we set the fMRI image as the moving image and the T1 as the fixed.  When we run ApplyTransforms, we pass #
    # it the inverse rigid transform file (or give it the 'inverse' flag)                                                       #
    ##########################################################################

    preproc.connect(tstat, 'out_file', bFit, 'movingVolume')
    preproc.connect(grabber, 't1_File', bFit, 'fixedVolume')

    forwardTransformT1ToFMRI = pipe.Node(Merger(1), name='0_List_forwardTransformT1ToFMRI')
    preproc.connect(bFit, 'outputTransform', forwardTransformT1ToFMRI, 'in1')

    # preprocessing_part2.sh
    warpT1ToFMRI = pipe.Node(interface=ApplyTransforms(), iterfield=['input_image'],
                             name='antsApplyTransformsT1')
    warpT1ToFMRI.inputs.interpolation = 'Linear'  # Validation test value: NearestNeighbor
    warpT1ToFMRI.inputs.invert_transform_flags = [True]
    preproc.connect(tstat, 'out_file', warpT1ToFMRI, 'reference_image')
    preproc.connect(grabber, 't1_File', warpT1ToFMRI, 'input_image')
    preproc.connect(forwardTransformT1ToFMRI, 'out', warpT1ToFMRI, 'transforms')

    preproc.connect(bFit, 'outputTransform', forwardTransform, 'in1')
    preproc.connect(bFit, 'outputTransform', reverseTransform, 'in2')

    warpBABCSegToFMRI = pipe.Node(interface=ApplyTransforms(), iterfield=['input_image'],
                                  name='warpBABCSegToFMRI')
    warpBABCSegToFMRI.inputs.interpolation = 'NearestNeighbor'
    warpBABCSegToFMRI.inputs.invert_transform_flags = [True]
    preproc.connect(tstat, 'out_file', warpBABCSegToFMRI, 'reference_image')
    preproc.connect(grabber, 'csfFile', warpBABCSegToFMRI, 'input_image')
    preproc.connect(forwardTransformT1ToFMRI, 'out', warpBABCSegToFMRI, 'transforms')

    # 10
    csfmask = pipe.Node(interface=Function(function=generateTissueMask,
                                           input_names=['input_file', 'fileName', 'low',
                                                        'high', 'erodeFlag', 'binary', 'largest'],
                                           output_names=['output_file']),
                        name='csfMask')
    csfmask.inputs.fileName = 'csfMask.nii'
    csfmask.inputs.low = 3
    csfmask.inputs.high = 42
    csfmask.inputs.erodeFlag = False
    csfmask.inputs.binary = True
    csfmask.inputs.largest = False
    csfmask.inputs.input_file = nacAtlasLabel

    warpAtlasCSF = warpT1ToFMRI.clone('antsApplyTransformCSF')
    warpAtlasCSF.inputs.invert_transform_flags = [True, False]
    preproc.connect(csfmask, 'output_file', warpAtlasCSF, 'input_image')  # 10
    preproc.connect(forwardTransform, 'out', warpAtlasCSF, 'transforms')
    preproc.connect(tstat, 'out_file', warpAtlasCSF, 'reference_image')

    warpWHM = warpT1ToFMRI.clone('antsApplyTransformWHM')
    preproc.connect(grabber, 'whmFile', warpWHM, 'input_image')
    preproc.connect(tstat, 'out_file', warpWHM, 'reference_image')
    preproc.connect(forwardTransformT1ToFMRI, 'out', warpWHM, 'transforms')

    # 11
    wmmask = pipe.Node(interface=Function(function=generateTissueMask,
                                          input_names=['input_file', 'fileName', 'low',
                                                       'high', 'erodeFlag', 'binary', 'largest'],
                                          output_names=['output_file']),
                       name='wmMask')
    wmmask.inputs.fileName = 'whiteMatterMask.nii'
    wmmask.inputs.low = 0.99
    wmmask.inputs.high = 1.0
    wmmask.inputs.erodeFlag = True
    wmmask.inputs.binary = False
    wmmask.inputs.largest = False
    preproc.connect(warpWHM, 'output_image', wmmask, 'input_file')  # 11

    csfAvg = pipe.Node(interface=Maskave(), name='afni3DmaskAve_csf')
    csfAvg.inputs.outputtype = 'AFNI_1D'  # outputType
    csfAvg.inputs.args = '-median'  # TODO
    csfAvg.inputs.quiet = True
    preproc.connect(warpAtlasCSF, 'output_image', csfAvg, 'mask')

    wmAvg = pipe.Node(interface=Maskave(), name='afni3DmaskAve_wm')
    wmAvg.inputs.outputtype = 'AFNI_1D'  # outputType
    wmAvg.inputs.args = '-median'  # TODO
    wmAvg.inputs.quiet = True
    preproc.connect(wmmask, 'output_file', wmAvg, 'mask')

    preproc.connect(calc, 'out_file', csfAvg, 'in_file')
    preproc.connect(calc, 'out_file', wmAvg, 'in_file')

    if maskGM:
        #------------------------------ GRAY MATTER MASK ------------------------------
        grabber.inputs.field_template['gryFile'] = probRegex
        grabber.inputs.template_args['gryFile'] = [
            ['session_id', 'ACCUMULATED_POSTERIORS', 'POSTERIOR_GM_TOTAL.nii.gz']]
        # For cerebrum gray matter ONLY:
        # grabber.inputs.template_args['gryFile'] = [['session_id', 'TissueClassify', 'POSTERIOR_SURFGM.nii.gz']]
        warpGRM = warpT1ToFMRI.clone('antsApplyTransformGRM')
        warpGRM.inputs.invert_transform_flags = [True]
        preproc.connect(grabber, 'gryFile', warpGRM, 'input_image')
        preproc.connect(tstat, 'out_file', warpGRM, 'reference_image')
        preproc.connect(forwardTransformT1ToFMRI, 'out', warpGRM, 'transforms')

        grmmask = pipe.Node(interface=Function(function=generateTissueMask,
                                               input_names=['input_file', 'fileName', 'low',
                                                            'high', 'erodeFlag', 'binary', 'largest'],
                                               output_names=['output_file']),
                            name='grmMask')
        grmmask.inputs.fileName = 'grmMask.nii'
        grmmask.inputs.low = 0.99
        grmmask.inputs.high = 1.0
        grmmask.inputs.erodeFlag = False
        grmmask.inputs.binary = False
        grmmask.inputs.largest = False
        preproc.connect(warpGRM, 'output_image', grmmask, 'input_file')

        grmAvg = pipe.Node(interface=Maskave(), name='afni3DmaskAve_grm')
        grmAvg.inputs.outputtype = 'AFNI_1D'  # outputType
        grmAvg.inputs.args = '-median'  # TODO
        grmAvg.inputs.quiet = True
        preproc.connect(grmmask, 'output_file', grmAvg, 'mask')
        preproc.connect(calc, 'out_file', grmAvg, 'in_file')
        # 12
        deconvolve = pipe.Node(interface=Deconvolve(fileCount=4, seriesCount=9), name='afni3Ddeconvolve')
        # deconvolve.inputs.outputtype = outputType # 'AFNI'
        deconvolve.inputs.ignoreWarnings = 10
        deconvolve.inputs.nullHypothesisPolynomialDegree = 1
        deconvolve.inputs.full_first = True
        deconvolve.inputs.is_float = True
        deconvolve.inputs.tout = True
        deconvolve.inputs.rout = True
        deconvolve.inputs.fout = True
        deconvolve.inputs.bucket = 'Rest_bp_Decon'
        deconvolve.inputs.fitts = 'full_fitts_Decon'
        deconvolve.inputs.errts = 'errts_Decon'
        deconvolve.inputs.stim_label_1 = "Median_CSF"
        deconvolve.inputs.is_stim_base_1 = False
        deconvolve.inputs.stim_label_2 = "Median_WM"
        deconvolve.inputs.is_stim_base_2 = False
        deconvolve.inputs.stim_label_3 = "Median_GM"
        deconvolve.inputs.is_stim_base_3 = False
        deconvolve.inputs.stim_label_4 = "roll"
        deconvolve.inputs.is_stim_base_4 = True
        deconvolve.inputs.stim_label_5 = 'pitch'
        deconvolve.inputs.is_stim_base_5 = True
        deconvolve.inputs.stim_label_6 = 'yaw'
        deconvolve.inputs.is_stim_base_6 = True
        deconvolve.inputs.stim_label_7 = 'dS'
        deconvolve.inputs.is_stim_base_7 = True
        deconvolve.inputs.stim_label_8 = 'dL'
        deconvolve.inputs.is_stim_base_8 = True
        deconvolve.inputs.stim_label_9 = 'dP'
        deconvolve.inputs.is_stim_base_9 = True
        preproc.connect(csfAvg, 'out_file', deconvolve, 'stim_file_1')
        preproc.connect(wmAvg, 'out_file', deconvolve, 'stim_file_2')
        preproc.connect(grmAvg, 'out_file', deconvolve, 'stim_file_3')
        preproc.connect(volreg, 'oned_file', deconvolve, 'stim_file_4')

    elif maskWholeBrain:
        # Mask the whole brain
        nacWholeBrainFile = "/Shared/paulsen/Experiments/rsFMRI/rs-fMRI-pilot/ReferenceAtlas/template_brain.nii.gz"

        warpWholeMask = warpT1ToFMRI.clone('antsApplyTransformWholeBrain')
        warpWholeMask.inputs.invert_transform_flags = [True, False]
        warpWholeMask.inputs.input_image = nacWholeBrainFile
        preproc.connect(forwardTransform, 'out', warpWholeMask, 'transforms')
        preproc.connect(tstat, 'out_file', warpWholeMask, 'reference_image')

        wholeBrainMask = wmmask.clone('wholeBrainMask')
        wholeBrainMask.inputs.fileName = 'wholeBrainMask.nii'
        # Invert the mask?
        wholeBrainMask.inputs.low = 0.5   # 0.0
        wholeBrainMask.inputs.high = 1.0  # 0.5
        wholeBrainMask.inputs.erodeFlag = False
        wholeBrainMask.inputs.binary = False
        wholeBrainMask.inputs.largest = True
        preproc.connect(warpWholeMask, 'output_image', wholeBrainMask, 'input_file')

        wholeMaskAvg = pipe.Node(interface=Maskave(), name='afni3DmaskAve_whole')
        wholeMaskAvg.inputs.outputtype = 'AFNI_1D'  # outputType
        wholeMaskAvg.inputs.args = '-median'  # TODO
        wholeMaskAvg.inputs.quiet = True
        preproc.connect(wholeBrainMask, 'output_file', wholeMaskAvg, 'mask')
        preproc.connect(calc, 'out_file', wholeMaskAvg, 'in_file')
        # 12
        deconvolve = pipe.Node(interface=Deconvolve(fileCount=4, seriesCount=9), name='afni3Ddeconvolve')
        # deconvolve.inputs.outputtype = outputType # 'AFNI'
        deconvolve.inputs.ignoreWarnings = 10
        deconvolve.inputs.nullHypothesisPolynomialDegree = 1
        deconvolve.inputs.full_first = True
        deconvolve.inputs.is_float = True
        deconvolve.inputs.tout = True
        deconvolve.inputs.rout = True
        deconvolve.inputs.fout = True
        deconvolve.inputs.bucket = 'Rest_bp_Decon'
        deconvolve.inputs.fitts = 'full_fitts_Decon'
        deconvolve.inputs.errts = 'errts_Decon'
        deconvolve.inputs.stim_label_1 = "Median_CSF"
        deconvolve.inputs.is_stim_base_1 = False
        deconvolve.inputs.stim_label_2 = "Median_WM"
        deconvolve.inputs.is_stim_base_2 = False
        deconvolve.inputs.stim_label_3 = "Median_WholeBrain"
        deconvolve.inputs.is_stim_base_3 = False
        deconvolve.inputs.stim_label_4 = "roll"
        deconvolve.inputs.is_stim_base_4 = True
        deconvolve.inputs.stim_label_5 = 'pitch'
        deconvolve.inputs.is_stim_base_5 = True
        deconvolve.inputs.stim_label_6 = 'yaw'
        deconvolve.inputs.is_stim_base_6 = True
        deconvolve.inputs.stim_label_7 = 'dS'
        deconvolve.inputs.is_stim_base_7 = True
        deconvolve.inputs.stim_label_8 = 'dL'
        deconvolve.inputs.is_stim_base_8 = True
        deconvolve.inputs.stim_label_9 = 'dP'
        deconvolve.inputs.is_stim_base_9 = True
        preproc.connect(csfAvg, 'out_file', deconvolve, 'stim_file_1')
        preproc.connect(wmAvg, 'out_file', deconvolve, 'stim_file_2')
        preproc.connect(wholeMaskAvg, 'out_file', deconvolve, 'stim_file_3')
        # preproc.connect(wholeBrainMask, 'out_file', deconvolve, 'mask')
        preproc.connect(volreg, 'oned_file', deconvolve, 'stim_file_4')

    else:
        deconvolve = pipe.Node(interface=Deconvolve(fileCount=3, seriesCount=8), name='afni3Ddeconvolve')
        # deconvolve.inputs.outputtype = outputType # 'AFNI'
        deconvolve.inputs.ignoreWarnings = 10
        deconvolve.inputs.nullHypothesisPolynomialDegree = 1
        deconvolve.inputs.full_first = True
        deconvolve.inputs.is_float = True
        deconvolve.inputs.tout = True
        deconvolve.inputs.rout = True
        deconvolve.inputs.fout = True
        deconvolve.inputs.bucket = 'Rest_bp_Decon'
        deconvolve.inputs.fitts = 'full_fitts_Decon'
        deconvolve.inputs.errts = 'errts_Decon'
        deconvolve.inputs.stim_label_1 = "Median_CSF"
        deconvolve.inputs.is_stim_base_1 = False
        deconvolve.inputs.stim_label_2 = "Median_WM"
        deconvolve.inputs.is_stim_base_2 = False
        deconvolve.inputs.stim_label_3 = "roll"
        deconvolve.inputs.is_stim_base_3 = True
        deconvolve.inputs.stim_label_4 = 'pitch'
        deconvolve.inputs.is_stim_base_4 = True
        deconvolve.inputs.stim_label_5 = 'yaw'
        deconvolve.inputs.is_stim_base_5 = True
        deconvolve.inputs.stim_label_6 = 'dS'
        deconvolve.inputs.is_stim_base_6 = True
        deconvolve.inputs.stim_label_7 = 'dL'
        deconvolve.inputs.is_stim_base_7 = True
        deconvolve.inputs.stim_label_8 = 'dP'
        deconvolve.inputs.is_stim_base_8 = True
        preproc.connect(csfAvg, 'out_file', deconvolve, 'stim_file_1')
        preproc.connect(wmAvg, 'out_file', deconvolve, 'stim_file_2')
        preproc.connect(volreg, 'oned_file', deconvolve, 'stim_file_3')
    # 13
    detrend = pipe.Node(interface=Detrend(), name='afni3Ddetrend')
    detrend.inputs.outputtype = outputType  # 'AFNI'
    detrend.inputs.suffix = '_dt'
    detrend.inputs.args = '-polort 3'  # TODO
    preproc.connect(deconvolve, 'out_errts', detrend, 'in_file')  # 13

    # Per Dawei, bandpass after running 3dDetrend
    preproc.connect(calc, 'out_file', deconvolve, 'in_file')
    preproc.connect(detrend, 'out_file', fourier, 'in_file')

    downsampleAtlas = pipe.Node(interface=Function(function=resampleImage,
                                                   input_names=['inputVolume', 'outputVolume', 'resolution'],
                                                   output_names=['outputVolume']),
                                name="downsampleAtlas")
    downsampleAtlas.inputs.inputVolume = nacAtlasFile
    downsampleAtlas.inputs.outputVolume = downsampledNACfilename
    downsampleAtlas.inputs.resolution = [int(x) for x in nacResampleResolution]

    fmriToNAC_epi = pipe.Node(interface=ApplyTransforms(), name='fmriToNac_epi')
    fmriToNAC_epi.inputs.interpolation = 'Linear'
    fmriToNAC_epi.inputs.invert_transform_flags = [False, False]
    # Fourier is the last NIFTI file format in the AFNI pipeline
    preproc.connect(fourier, 'out_file', fmriToNAC_epi, 'input_image')
    preproc.connect(downsampleAtlas, 'outputVolume', fmriToNAC_epi, 'reference_image')
    preproc.connect(reverseTransform, 'out', fmriToNAC_epi, 'transforms')

    # Create seed points
    labels, seeds = getAtlasPoints('seeds.fcsv')

    seedsIdentity = pipe.Node(interface=IdentityInterface(fields=['index']),
                              name='seedsIdentity')
    seedsIdentity.iterables = ('index', range(len(labels)))

    selectSeed = pipe.Node(interface=Select(), name='selectSeed')
    selectSeed.inputs.inlist = seeds
    preproc.connect(seedsIdentity, 'index', selectSeed, 'index')

    selectLabel = pipe.Node(interface=Select(), name='selectLabel')
    selectLabel.inputs.inlist = labels
    preproc.connect(seedsIdentity, 'index', selectLabel, 'index')

    points = pipe.Node(interface=Function(function=createSphereExpression,
                                          input_names=['coordinates', 'radius'],
                                          output_names=['expression']),
                       name='createSphereExpression')
    preproc.connect(selectSeed, 'out', points, 'coordinates')

    spheres = pipe.Node(interface=Calc(letters=['a']),
                        name='afni3Dcalc_seeds')
    spheres.inputs.outputtype = outputType
    preproc.connect(downsampleAtlas, 'outputVolume', spheres, 'in_file_a')
    spheres.inputs.args = '-nscale'

    preproc.connect(points, 'expression', spheres, 'expr')

    renameMasks = pipe.Node(interface=Rename(format_string='%(label)s_mask'), name='renameMasksAtlas')
    renameMasks.inputs.keep_ext = True
    preproc.connect(selectLabel, 'out', renameMasks, 'label')
    preproc.connect(spheres, 'out_file', renameMasks, 'in_file')

    atlas_DataSink = pipe.Node(interface=DataSink(), name="atlas_DataSink")
    atlas_DataSink.inputs.base_directory = preproc.base_dir  # '/Shared/paulsen/Experiments/20130417_rsfMRI_Results'
    atlas_DataSink.inputs.container = RESULTS_DIR
    atlas_DataSink.inputs.parameterization = False
    atlas_DataSink.overwrite = REWRITE_DATASINKS
    preproc.connect(renameMasks, 'out_file', atlas_DataSink, 'Atlas')
    preproc.connect(downsampleAtlas, 'outputVolume', atlas_DataSink, 'Atlas.@resampled')

    # Warp seed output to FMRI--this will be the input into the clipping code. JV -- Dec. 17 2013
    nacToFMRI = pipe.Node(interface=ApplyTransforms(), name="warpNACseedsToFMRI")
    nacToFMRI.inputs.interpolation = 'NearestNeighbor'
    nacToFMRI.inputs.invert_transform_flags = [True, False]
    preproc.connect(spheres, 'out_file', nacToFMRI, 'input_image')
    preproc.connect(warpT1ToFMRI, 'output_image', nacToFMRI, 'reference_image')
    preproc.connect(forwardTransform, 'out', nacToFMRI, 'transforms')

    renameMasks2 = pipe.Node(interface=Rename(format_string='%(session)s_%(label)s_mask'), name='renameMasksFMRI')
    renameMasks2.inputs.keep_ext = True
    preproc.connect(selectLabel, 'out', renameMasks2, 'label')
    preproc.connect(sessions, 'session_id', renameMasks2, 'session')

    clipSeedWithVentriclesNode = pipe.Node(interface=Function(function=clipSeedWithVentricles,
                                           input_names=['unclipped_seed_fn', 'fmriBABCSeg_fn', 'desired_out_seed_fn'],
                                           output_names=['clipped_seed_fn']),
                                           name='clipSeedWithVentriclesNode')

    clipSeedWithVentriclesNode.inputs.desired_out_seed_fn = "clipped_seed.nii.gz"

    preproc.connect(nacToFMRI, 'output_image', clipSeedWithVentriclesNode, 'unclipped_seed_fn')
    preproc.connect(warpBABCSegToFMRI, 'output_image', clipSeedWithVentriclesNode, 'fmriBABCSeg_fn')
    if not maskWhiteMatterFromSeeds:
        preproc.connect(clipSeedWithVentriclesNode, 'clipped_seed_fn', renameMasks2, 'in_file')
    else:
        clipSeedWithWhiteMatterNode = pipe.Node(interface=Function(function=clipSeedWithWhiteMatter,
                                                                   input_names=['seed', 'mask', 'outfile'],
                                                                   output_names=['outfile']),
                                                name='clipSeedWithWhiteMatterNode')
        clipSeedWithWhiteMatterNode.inputs.outfile = 'clipped_wm_seed.nii.gz'
        preproc.connect(warpBABCSegToFMRI, 'output_image', clipSeedWithWhiteMatterNode, 'mask')
        preproc.connect(clipSeedWithVentriclesNode, 'clipped_seed_fn', clipSeedWithWhiteMatterNode, 'seed')
        preproc.connect(clipSeedWithWhiteMatterNode, 'outfile', renameMasks2, 'in_file')

    # Labels are iterated over, so we need a seperate datasink to avoid overwriting any preprocessing
    # results when the labels are iterated (e.g. To3d output)
    fmri_label_DataSink = pipe.Node(interface=DataSink(), name='fmri_label_DataSink')
    # '/Shared/paulsen/Experiments/20130417_rsfMRI_Results/EPI'
    fmri_label_DataSink.inputs.base_directory = os.path.join(preproc.base_dir, RESULTS_DIR, 'EPI')
    fmri_label_DataSink.inputs.parameterization = False
    fmri_label_DataSink.overwrite = REWRITE_DATASINKS
    preproc.connect(sessions, 'session_id', fmri_label_DataSink, 'container')
    preproc.connect(renameMasks2, 'out_file', fmri_label_DataSink, 'masks')
    preproc.connect(fourier, 'out_file', fmri_label_DataSink, 'masks.@bandpass')

    roiMedian = pipe.Node(interface=Maskave(), name='afni_roiMedian')
    roiMedian.inputs.outputtype = 'AFNI_1D'
    roiMedian.inputs.args = '-median -mrange 1 1'  # TODO
    roiMedian.inputs.quiet = True
    preproc.connect(renameMasks2, 'out_file', roiMedian, 'mask')
    preproc.connect(fourier, 'out_file', roiMedian, 'in_file')

    correlate = pipe.Node(interface=Fim(), name='afni_correlate')
    correlate.inputs.out = 'Correlation'
    preproc.connect(roiMedian, 'out_file', correlate, 'ideal_file')
    preproc.connect(fourier, 'out_file', correlate, 'in_file')

    regionLogCalc = pipe.Node(interface=Calc(letters=['a']), name='afni_regionLogCalc')
    regionLogCalc.inputs.outputtype = outputType
    regionLogCalc.inputs.expr = 'log((1+a)/(1-a))/2'
    preproc.connect(correlate, 'out_file', regionLogCalc, 'in_file_a')

    renameZscore = pipe.Node(interface=Rename(format_string="%(session)s_%(label)s_zscore"), name='renameZscore')
    renameZscore.inputs.keep_ext = True
    preproc.connect(sessions, 'session_id', renameZscore, 'session')
    preproc.connect(selectLabel, 'out', renameZscore, 'label')
    preproc.connect(regionLogCalc, 'out_file', renameZscore, 'in_file')
    preproc.connect(renameZscore, 'out_file', fmri_label_DataSink, 'zscores')

    # Move z values back into NAC atlas space
    fmriToNAC_label = fmriToNAC_epi.clone(name='fmriToNac_label')
    fmriToNAC_label.inputs.interpolation = 'Linear'
    fmriToNAC_label.inputs.invert_transform_flags = [False, False]
    preproc.connect(downsampleAtlas, 'outputVolume', fmriToNAC_label, 'reference_image')
    preproc.connect(regionLogCalc, 'out_file', fmriToNAC_label, 'input_image')
    preproc.connect(reverseTransform, 'out', fmriToNAC_label, 'transforms')

    renameZscore2 = pipe.Node(interface=Rename(format_string="%(session)s_%(label)s_result"), name='renameZscore2')
    renameZscore2.inputs.keep_ext = True
    preproc.connect(sessions, 'session_id', renameZscore2, 'session')
    preproc.connect(selectLabel, 'out', renameZscore2, 'label')
    preproc.connect(fmriToNAC_label, 'output_image', renameZscore2, 'in_file')
    preproc.connect(renameZscore2, 'out_file', atlas_DataSink, 'Atlas.@zscore')

    preproc.write_graph()
    # preproc.write_hierarchical_dotfile(dotfilename='dave.dot')
    if os.environ['USER'] == 'dmwelch' and False:
        # Run preprocessing on the local cluster
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
