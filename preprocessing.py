#!/usr/bin/env python
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

### TODO: Modify virtualenv to include formatFMRI.sh

import argparse
import os
import sys

import SEMTools as sem
from nipype.interfaces.ants.registration import Registration
from nipype.interfaces.ants.resampling import ApplyTransforms
from nipype.interfaces.afni.preprocess import *
from nipype.interfaces.freesurfer.preprocess import *
from nipype.interfaces.io import DataSink, DataGrabber
from nipype.interfaces.utility import Function, IdentityInterface, Rename, Select
from nipype.interfaces.utility import Merge as Merger
import nipype.pipeline.engine as pipe
import numpy

from utilities import *

def pipeline(args):
    sessionID = args.sessionID
    outputType = args.outputType
    fOutputType = args.fOutputType
    print args.name
    preproc = pipe.Workflow(updatehash=True, name="workflow_" + args.name)
    if  args.pipelineOption == 'iowa':
        preproc.name += "_iowa"
    else:
        preproc.name += "_csail"
    preproc.base_dir = os.getcwd()

    fmri_DataSink = pipe.Node(interface=DataSink(), name="fmri_DataSink")
    ### HACK: Remove node from pipeline until Nipype/AFNI file copy issue is resolved
    # fmri_DataSink.inputs.base_directory = os.path.join(preproc.base_dir, 'results', 'fmri') # '/paulsen/Experiments/20130417_rsfMRI_Results'
    # fmri_DataSink.inputs.substitutions = [('to_3D_out+orig', 'to3D')]
    # fmri_DataSink.inputs.parameterization = False
    ### END HACK

    sessions = pipe.Node(interface=IdentityInterface(fields=['session_id']), name='sessionIDs')
    sessions.iterables = ('session_id', sessionID)

    ### HACK: Remove node from pipeline until Nipype/AFNI file copy issue is resolved
    # preproc.connect(sessions, 'session_id', fmri_DataSink, 'container')
    ### END HACK

    grabber = pipe.Node(interface=DataGrabber(infields=['session_id'],
                                              outfields=['fmri_dicom_dir', 't1_File',
                                                         'faparc_File', 'faparc2009_File', 'csfFile',
                                                         'whmFile']), name='dataGrabber')
    grabber.inputs.base_directory = '/paulsen'
    grabber.inputs.template = '*'
    site = 'FMRI_HD_024' ### HACK: site hardcoded
    fmriRegex = 'MRx/{site}/*/%s/%s/%s/*'.format(site=site)
    fS_Regex = 'Experiments/20120722_JOY_DWI/{site}/*/%s/%s/%s_*_%s_FS/%s/%s'.format(site=site)
    probRegex = 'Experiments/20120801.SubjectOrganized_Results/{site}/*/%s/%s/%s'.format(site=site)
    grabber.inputs.field_template = dict(fmri_dicom_dir=fmriRegex,
                                         t1_File=fS_Regex,
                                         faparc_File=fS_Regex,
                                         faparc2009_File=fS_Regex,
                                         csfFile=probRegex,
                                         whmFile=probRegex)
    grabber.inputs.template_args = dict(fmri_dicom_dir=[['session_id', 'ANONRAW',
                                                         'FMRI_RestingStateConnectivity']],
                                        t1_File=[['session_id', '10_AUTO.NN3Tv20110419',
                                                  'JOY_v51_2011', 'session_id', 'mri', 'brain.mgz']],
                                        faparc_File=[['session_id', '10_AUTO.NN3Tv20110419',
                                                  'JOY_v51_2011', 'session_id', 'mri_nifti',
                                                  'aparc+aseg.nii.gz']],
                                        faparc2009_File=[['session_id', '10_AUTO.NN3Tv20110419',
                                                  'JOY_v51_2011', 'session_id', 'mri_nifti',
                                                  'aparc.a2009s+aseg.nii.gz']],
                                        csfFile=[['session_id', 'ACCUMULATED_POSTERIORS',
                                                   'POSTERIOR_CSF_TOTAL.nii.gz']],
                                        whmFile=[['session_id', 'ACCUMULATED_POSTERIORS',
                                                  'POSTERIOR_WM_TOTAL.nii.gz']])
    preproc.connect(sessions, 'session_id', grabber, 'session_id')

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
    preproc.connect(formatFMRINode, 'numberOfSlices', to_3D_str, 'slices')                 #1
    preproc.connect(formatFMRINode, 'numberOfFiles', to_3D_str, 'volumes')
    preproc.connect(formatFMRINode, 'repetitionTime', to_3D_str, 'repTime')
    preproc.connect(formatFMRINode, 'sliceOrder', to_3D_str, 'order')

    to_3D = pipe.Node(interface=To3D(), name='afniTo3D')
    to_3D.inputs.outputtype = 'AFNI'
    to_3D.inputs.datum = 'short'
    to_3D.inputs.filetype = 'epan'
    to_3D.inputs.prefix = 'to_3D_out'
    preproc.connect(grabber, 'fmri_dicom_dir', to_3D, 'infolder')
    preproc.connect(to_3D_str, 'out_string', to_3D, 'funcparams')

    ### HACK: Remove node from pipeline until Nipype/AFNI file copy issue is resolved
    # renameTo3D = pipe.Node(Rename(format_string='%(session)s_to3D'), name='renameTo3D')
    # renameTo3D.inputs.keep_ext = True
    # preproc.connect(to_3D, 'out_file', renameTo3D, 'in_file')
    # preproc.connect(sessions, 'session_id', renameTo3D, 'session')
    # preproc.connect(renameTo3D, 'out_file', fmri_DataSink, '@To3D')
    ### END HACK

    #1a
    refit = pipe.Node(interface=Refit(), name='afni3Drefit')
    refit.inputs.outputtype = 'AFNI'
    refit.inputs.deoblique = True
    preproc.connect(to_3D, 'out_file', refit, 'in_file')

    #2
    skipCount = 4
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

    #3
    volreg = pipe.Node(interface=Volreg(), name='afni3DvolReg')
    volreg.inputs.outputtype = outputType
    volreg.inputs.timeshift = False # 0
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
    preproc.connect(despike, 'out_file', volreg, 'in_file')               #3

    #4
    zpad = pipe.Node(interface=Zeropad(), name='afniZeropad')
    zpad.inputs.plane = 'IS'
    zpad.inputs.numberOfPlanes = 44
    zpad.inputs.is_mm = False
    preproc.connect(volreg, 'out_file', zpad, 'in_file')                  #4

    #5
    merge = pipe.Node(interface=Merge(), name='afni3Dmerge')
    merge.inputs.outputtype = outputType
    merge.inputs.blurfwhm = 6
    merge.inputs.doall = True
    merge.inputs.args = '-1noneg -1clip 100' # TODO
    # TODO: implement
    # merge.inputs.onenoneg = True
    # merge.inputs.oneclip = 100
    preproc.connect(zpad, 'out_file', merge, 'in_files')                  #5 ### NOTE: ONLY ONE FILE TO in_files

    #6
    automask = pipe.Node(interface=Automask(), name='afni3Dautomask')
    automask.inputs.outputtype = outputType
    automask.inputs.dilate = 1
    preproc.connect(merge, 'out_file', automask, 'in_file')               #6

    #7
    tstat = pipe.Node(interface=TStat(), name='afni3DtStat')
    tstat.inputs.outputtype = outputType
    tstat.inputs.args = '-mean' # TODO
    preproc.connect(merge, 'out_file', tstat, 'in_file')                  #7 ### 'mean_file' -> 'in_file'
    preproc.connect(automask, 'out_file', tstat, 'mask_file')

    calc = pipe.Node(interface=Calc(letters=['a', 'b']), name='afni3Dcalc')
    calc.inputs.outputtype = outputType
    calc.inputs.expr = "a * b"
    preproc.connect(merge, 'out_file', calc, 'in_file_a')
    preproc.connect(automask, 'out_file', calc, 'in_file_b')

    #8
    fourier = pipe.Node(interface=Fourier(), name='afni3Dfourier')
    fourier.inputs.outputtype = outputType
    fourier.inputs.highpass = 0.011
    fourier.inputs.lowpass = 0.1
    fourier.inputs.args = '-retrend' # TODO
    fourier.inputs.outputtype = outputType
    preproc.connect(calc, 'out_file', fourier, 'in_file')

    #9
    converter = pipe.Node(interface=MRIConvert(), name='freesurferMRIconvert')
    converter.inputs.out_type = fOutputType
    converter.inputs.in_type = 'mgz'
    converter.inputs.force_ras = True
    converter.inputs.out_datatype = 'short'

    convertT1 = converter.clone('T1Converter')
    preproc.connect(grabber, 't1_File', convertT1, 'in_file')

    convertF1 = converter.clone('F1Converter')
    convertF2 = converter.clone('F2Converter')

    bFit = pipe.Node(interface=sem.BRAINSFit(), name='brainsFit')
    bFit.inputs.initializeTransformMode = 'useCenterOfHeadAlign'
    bFit.inputs.maskProcessingMode = 'ROIAUTO'
    bFit.inputs.ROIAutoDilateSize = 10.0
    bFit.inputs.useRigid = True
    bFit.inputs.costMetric = 'MMI' # (default)
    bFit.inputs.numberOfSamples = 100000 # (default)
    bFit.inputs.outputTransform = True
    preproc.connect(convertT1, 'out_file', bFit, 'movingVolume')
    preproc.connect(tstat, 'out_file', bFit, 'fixedVolume')

    ### preprocessing_part2.sh
    warpT1ToFMRI = pipe.Node(interface=ApplyTransforms(), iterfield=['input_image'],
                             name='antsApplyTransformsT1')
    warpT1ToFMRI.inputs.interpolation='Linear' # Validation test value: NearestNeighbor
    preproc.connect(convertT1, 'out_file', warpT1ToFMRI, 'input_image') # connected to brain.nii NOT brain.mgz
    preproc.connect(tstat, 'out_file', warpT1ToFMRI, 'reference_image')
    preproc.connect([(bFit, warpT1ToFMRI, [(('outputTransform', makeList), 'transforms')])])

    # warpT1ToFMRI = warpNN.clone('antsApplyTransformT1')
    warpCSF = warpT1ToFMRI.clone('antsApplyTransformCSF')
    preproc.connect(grabber, 'csfFile', warpCSF, 'input_image')
    preproc.connect(tstat, 'out_file', warpCSF, 'reference_image')
    preproc.connect([(bFit, warpCSF, [(('outputTransform', makeList), 'transforms')])])

    warpWHM = warpT1ToFMRI.clone('antsApplyTransformWHM')
    preproc.connect(grabber, 'whmFile', warpWHM, 'input_image')
    preproc.connect(tstat, 'out_file', warpWHM, 'reference_image')
    preproc.connect([(bFit, warpWHM, [(('outputTransform', makeList), 'transforms')])])

    #10
    csfmask = pipe.Node(interface=Function(function=generateTissueMask,
                                           input_names=['input_file','low', 'high', 'erodeFlag'],
                                           output_names=['output_file']),
                        name='csfMask')
    csfmask.inputs.low = 0.99
    csfmask.inputs.high = 1.0
    csfmask.inputs.erodeFlag = False
    preproc.connect(warpCSF, 'output_image', csfmask, 'input_file')       #10

    #11
    wmmask = pipe.Node(interface=Function(function=generateTissueMask,
                                           input_names=['input_file','low', 'high', 'erodeFlag'],
                                           output_names=['output_file']),
                        name='wmMask')
    wmmask.inputs.low = 0.99
    wmmask.inputs.high = 1.0
    wmmask.inputs.erodeFlag = True
    preproc.connect(warpWHM, 'output_image', wmmask, 'input_file')        #11

    csfAvg = pipe.Node(interface=Maskave(), name='afni3DmaskAve_csf')
    csfAvg.inputs.outputtype = 'AFNI_1D' #outputType
    csfAvg.inputs.args = '-median' # TODO
    csfAvg.inputs.quiet = True
    preproc.connect(csfmask, 'output_file', csfAvg, 'mask')
    preproc.connect(fourier, 'out_file', csfAvg, 'in_file')

    wmAvg = pipe.Node(interface=Maskave(), name='afni3DmaskAve_wm')
    wmAvg.inputs.outputtype = 'AFNI_1D' #outputType
    wmAvg.inputs.args = '-median' # TODO
    wmAvg.inputs.quiet = True
    preproc.connect(wmmask, 'output_file', wmAvg, 'mask')
    preproc.connect(fourier, 'out_file', wmAvg, 'in_file')

    #12
    deconvolve = pipe.Node(interface=Deconvolve(fileCount=3, seriesCount=8), name='afni3Ddeconvolve')
    # deconvolve.inputs.outputtype = outputType # 'AFNI'
    deconvolve.inputs.ignoreWarnings = 4
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
    preproc.connect(fourier, 'out_file', deconvolve, 'in_file')
    preproc.connect(csfAvg, 'out_file', deconvolve, 'stim_file_1')
    preproc.connect(wmAvg, 'out_file', deconvolve, 'stim_file_2')
    preproc.connect(volreg, 'oned_file', deconvolve, 'stim_file_3')

    #13
    detrend = pipe.Node(interface=Detrend(), name='afni3Ddetrend')
    detrend.inputs.outputtype = outputType # 'AFNI'
    detrend.inputs.suffix = '_dt'
    detrend.inputs.args = '-polort 3' # TODO
    preproc.connect(deconvolve, 'out_errts', detrend, 'in_file')          #13

    if args.pipelineOption == 'iowa':
        nacAtlasFile = "/IPLlinux/ipldev/scratch/welchdm/bld/rs-fMRI/BSA/ReferenceAtlas-build/Atlas/Atlas_20130106/template_t1.nii.gz"

        transformGrabber = pipe.Node(interface=DataGrabber(infields=['session_id'],
                                                           outfields=['atlasToSessionTransform',
                                                                      'sessionToAtlasTransform']),
                                     name='transformGrabber')
        transformGrabber.inputs.base_directory = '/paulsen/Experiments/20130202_PREDICTHD_Results/SubjectToAtlasWarped'
        transformGrabber.inputs.template = '*'
        transformRegex = '%s/AtlasToSubject_%sComposite.h5'
        transformGrabber.inputs.field_template = dict(atlasToSessionTransform=transformRegex,
                                                      sessionToAtlasTransform=transformRegex)
        transformGrabber.inputs.template_args = dict(atlasToSessionTransform=[['session_id', '']],
                                                     sessionToAtlasTransform=[['session_id', 'Inverse']])

        preproc.connect(sessions, 'session_id', transformGrabber, 'session_id')

        # Concatenate transforms: NAC -> T1 -> fMRI
        concatTransformsNode = pipe.Node(interface=Function(function=concatTransforms,
                                                            input_names=['transform1', 'transform2'],
                                                            output_names=['transformList']),
                                         name='concatTransforms')
        preproc.connect(bFit, 'outputTransform', concatTransformsNode, 'transform1')
        preproc.connect(transformGrabber, 'atlasToSessionTransform', concatTransformsNode, 'transform2')

        # Concatenate transforms: NAC <- T1 <- fMRI
        invertConcatTransformsNode = concatTransformsNode.clone(name='invertConcatTransforms')
        preproc.connect(bFit, 'outputTransform', invertConcatTransformsNode, 'transform1')
        preproc.connect(transformGrabber, 'sessionToAtlasTransform', invertConcatTransformsNode, 'transform2')

        fmriToNAC_epi = pipe.Node(interface=ApplyTransforms(), name='fmriToNac_epi')
        fmriToNAC_epi.inputs.interpolation = 'Linear'
        fmriToNAC_epi.inputs.invert_transform_flags = [True, False]
        fmriToNAC_epi.inputs.reference_image = nacAtlasFile
        preproc.connect(detrend, 'out_file', fmriToNAC_epi, 'input_image') # Detrend is the last NIFTI file format in the AFNI pipeline
        preproc.connect(invertConcatTransformsNode, 'transformList', fmriToNAC_epi, 'transforms')

        ### Create seed points
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
        spheres.inputs.in_file_a = nacAtlasFile
        spheres.inputs.args = '-nscale'

        preproc.connect(points, 'expression', spheres, 'expr')

        renameMasks = pipe.Node(interface=Rename(format_string='%(label)s_mask'), name='renameMasksAtlas')
        renameMasks.inputs.keep_ext = True
        preproc.connect(selectLabel, 'out', renameMasks, 'label')
        preproc.connect(spheres, 'out_file', renameMasks, 'in_file')

        atlas_DataSink = pipe.Node(interface=DataSink(), name="atlas_DataSink")
        atlas_DataSink.inputs.base_directory = preproc.base_dir # '/paulsen/Experiments/20130417_rsfMRI_Results'
        atlas_DataSink.inputs.container = 'results'
        atlas_DataSink.inputs.parameterization = False
        preproc.connect(renameMasks, 'out_file', atlas_DataSink, 'Atlas')

        # Warp seed output to FMRI
        nacToFMRI = pipe.Node(interface=ApplyTransforms(), name="nacToFMRI")
        nacToFMRI.inputs.interpolation = 'NearestNeighbor'
        preproc.connect(spheres, 'out_file', nacToFMRI, 'input_image')
        preproc.connect(warpT1ToFMRI, 'output_image', nacToFMRI, 'reference_image')
        preproc.connect(concatTransformsNode, 'transformList', nacToFMRI, 'transforms')

        ### TEST
        fmriToNAC_test = fmriToNAC_epi.clone(name='fmriToNac_test')
        fmriToNAC_test.inputs.interpolation = 'NearestNeighbor'
        fmriToNAC_test.inputs.invert_transform_flags = [True, False]
        fmriToNAC_test.inputs.reference_image = nacAtlasFile
        preproc.connect(nacToFMRI, 'output_image', fmriToNAC_test, 'input_image')
        preproc.connect(invertConcatTransformsNode, 'transformList', fmriToNAC_test, 'transforms')
        ### END TEST

        renameMasks2 = pipe.Node(interface=Rename(format_string='%(session)s_%(label)s_mask'), name='renameMasksFMRI')
        renameMasks2.inputs.keep_ext = True
        preproc.connect(selectLabel, 'out', renameMasks2, 'label')
        preproc.connect(sessions, 'session_id', renameMasks2, 'session')
        preproc.connect(nacToFMRI, 'output_image', renameMasks2, 'in_file')

        # Labels are iterated over, so we need a seperate datasink to avoid overwriting any preprocessing
        # results when the labels are iterated (e.g. To3d output)
        fmri_label_DataSink = fmri_DataSink.clone(name='fmri_label_DataSink')
        fmri_label_DataSink.inputs.base_directory = os.path.join(preproc.base_dir, 'Results', 'EPI') # '/paulsen/Experiments/20130417_rsfMRI_Results/EPI'
        fmri_label_DataSink.inputs.parameterization = False
        preproc.connect(sessions, 'session_id', fmri_label_DataSink, 'container')
        preproc.connect(renameMasks2, 'out_file', fmri_label_DataSink, 'masks')
        preproc.connect(detrend, 'out_file', fmri_label_DataSink, 'masks.@detrend')

        roiMedian = pipe.Node(interface=Maskave(), name='afni_roiMedian')
        roiMedian.inputs.outputtype = 'AFNI_1D'
        roiMedian.inputs.args = '-median -mrange 1 1'  ### TODO
        roiMedian.inputs.quiet = True
        preproc.connect(nacToFMRI, 'output_image', roiMedian, 'mask')
        preproc.connect(detrend, 'out_file', roiMedian, 'in_file')

        correlate = pipe.Node(interface=Fim(), name='afni_correlate')
        correlate.inputs.out = 'Correlation'
        preproc.connect(roiMedian, 'out_file', correlate, 'ideal_file')
        preproc.connect(deconvolve, 'out_errts', correlate, 'in_file')

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

        ### Move z values back into NAC atlas space
        fmriToNAC_label = fmriToNAC_epi.clone(name='fmriToNac_label')
        fmriToNAC_label.inputs.interpolation = 'Linear'
        fmriToNAC_label.inputs.invert_transform_flags = [True, False]
        fmriToNAC_label.inputs.reference_image = nacAtlasFile
        preproc.connect(regionLogCalc, 'out_file', fmriToNAC_label, 'input_image')
        preproc.connect(invertConcatTransformsNode, 'transformList', fmriToNAC_label, 'transforms')

        renameZscore2 = pipe.Node(interface=Rename(format_string="%(session)s_%(label)s_result"), name='renameZscore2')
        renameZscore2.inputs.keep_ext = True
        preproc.connect(sessions, 'session_id', renameZscore2, 'session')
        preproc.connect(selectLabel, 'out', renameZscore2, 'label')
        preproc.connect(fmriToNAC_label, 'output_image', renameZscore2, 'in_file')
        preproc.connect(renameZscore2, 'out_file', atlas_DataSink, 'atlas.@zscore')

    elif args.pipelineOption == 'csail':

        warpFS1 = pipe.Node(interface=ApplyTransforms(), name='antsApplyTransformsFS1')
        warpFS1.inputs.interpolation='MultiLabel'
        preproc.connect([(bFit, warpFS1, [(('outputTransform', makeList), 'transforms')])])
        preproc.connect(grabber, 'faparc_File', warpFS1, 'input_image')
        preproc.connect(tstat, 'out_file', warpFS1, 'reference_image')

        warpFS2 = pipe.Node(interface=ApplyTransforms(), name='antsApplyTransformsFS2')
        warpFS2.inputs.interpolation='MultiLabel'
        preproc.connect([(bFit, warpFS2, [(('outputTransform', makeList), 'transforms')])])
        preproc.connect(grabber, 'faparc2009_File', warpFS2, 'input_image')
        preproc.connect(tstat, 'out_file', warpFS2, 'reference_image')

        ###        2) Binary threshold label image for each label
        ###        3) Calculate the average fMRI value per label
        roiStats1 = pipe.Node(interface=ROIStats(), name='afni3DroiStats1')
        roiStats1.inputs.outputtype = 'AFNI_1D'
        roiStats1.inputs.args = '-nzmedian -nomeanout'
        preproc.connect(warpFS1, 'output_image', roiStats1, 'mask')
        preproc.connect(detrend, 'out_file', roiStats1, 'in_file')

        roiStats2 = roiStats1.clone(name='afni3DroiStats2')
        preproc.connect(warpFS2, 'output_image', roiStats2, 'mask')
        preproc.connect(detrend, 'out_file', roiStats2, 'in_file')

        ###         4) Output label covariance matrix to a file
        writeFile1 = pipe.Node(interface=Function(function=writeMat,
                                                  input_names=['in_file'],
                                                  output_names=['corr_file', 'label_file', 'raw_file']),
                               name='writeMatFile1')
        preproc.connect(roiStats1, 'stats', writeFile1, 'in_file')

        writeFile2 = writeFile1.clone(name='writeMatFile2')
        preproc.connect(roiStats2, 'stats', writeFile2, 'in_file')

        sinker = pipe.Node(interface=DataSink(), name="rsDataSink")
        sinker.inputs.base_directory = '/paulsen/Experiments/rsFMRI-test/csail_Results'
        find_matlab = 'errts_Decon+orig_dtroiStat'
        repl_matlab = 'freesurfer'
        sinker.inputs.substitutions = [(find_matlab, repl_matlab)]
        sinker.inputs.regexp_substitutions = [('([a-zA-Z0-9_]*)/freesurfer', 'freesurfer')]
        def createStr(value):
            return str(value)

        preproc.connect(sessions, ('session_id', createStr), sinker, 'container')
        preproc.connect(writeFile1, 'corr_file', sinker, 'aparc.@corr')
        preproc.connect(writeFile1, 'label_file', sinker, 'aparc.@label')
        preproc.connect(writeFile1, 'raw_file', sinker, 'aparc.@raw')

        preproc.connect(writeFile2, 'corr_file', sinker, 'aparca2009s.@corr')
        preproc.connect(writeFile2, 'label_file', sinker, 'aparca2009s.@label')
        preproc.connect(writeFile2, 'raw_file', sinker, 'aparca2009s.@raw')

    preproc.write_graph()
    # preproc.write_hierarchical_dotfile(dotfilename='dave.dot')
    preproc.run(plugin='MultiProc', plugin_args={'n_proc':12})

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Preprocessing script for resting state fMRI')
    parser.add_argument('-t', '--filetype', action='store', dest='outputType', type=str, required=False,
                        default='NIFTI', help='File type for outputs.  Values: "NIFTI_GZ", "AFNI", "NIFTI" (default)')
    parser.add_argument('-n','--name', action='store', dest='name', type=str, required=True,
                        help='Name (required)')
    parser.add_argument('-s', '--session', action='store', dest='sessionID', type=str, required=True,
                        nargs='*', help='list of session IDs (required)')
    parser.add_argument('-p', '--pipelineOption', action='store', dest='pOption', type=str, required=False, default='csail',
                        help='Specify the pipeline to execute.  Values: "IOWA", "CSAIL"(default)')
    args = parser.parse_args()
    freesurferOutputTypes = {"NIFTI_GZ" : "niigz",
                             "AFNI" : "afni",
                             "NIFTI" : "nii"}
    pipelineOptions = {"CSAIL" : "csail",
                       "IOWA" : "iowa"}
    args.fOutputType = freesurferOutputTypes[args.outputType]
    args.pipelineOption = pipelineOptions[args.pOption]
    outvalue = pipeline(args)
    sys.exit(outvalue)
