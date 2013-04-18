#!/usr/bin/env python
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

import argparse
import os
import sys

import SEMTools as sem
from nipype.interfaces.ants.registration import Registration
from nipype.interfaces.ants.resampling import ApplyTransforms
from nipype.interfaces.afni.preprocess import *
from nipype.interfaces.freesurfer.preprocess import *
from nipype.interfaces.io import DataSink, DataGrabber
from nipype.interfaces.utility import Function, IdentityInterface
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

    ## preprocessing_Nov1_update.sh
    sessions = pipe.Node(interface=IdentityInterface(fields=['session_id']), name='sessionIDs')
    sessions.iterables = ('session_id', sessionID)

    grabber = pipe.Node(interface=DataGrabber(infields=['session_id'],
                                              outfields=['fmri_dicom_dir', 't1_File',
                                                         'faparc_File', 'faparc2009_File', 'csfFile',
                                                         'whmFile']), name='dataGrabber')
    grabber.inputs.base_directory = '/paulsen'
    grabber.inputs.template = '*'
    fmriRegex = 'MRx/*/*/%s/%s/%s/*'
    fS_Regex = 'Experiments/20120722_JOY_DWI/*/*/%s/%s/%s_*_%s_FS/%s/%s'
    probRegex = 'Experiments/20130729_PREDICT_Results/*/*/%s/%s/%s'
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
    # if len(sessionID) > 1:
    preproc.connect(sessions, 'session_id', grabber, 'session_id')
    # else:
    #     grabber.inputs.session_id = sessionID[0]

    dwiConverter = pipe.Node(interface=sem.DWIConvert(), name='dwiConverter')
    dwiConverter.inputs.conversionMode = 'DicomToNrrd'
    dwiConverter.inputs.outputDirectory = '.'  # Cache directory
    dwiConverter.inputs.outputVolume = 'converted.nrrd'
    preproc.connect(grabber, 'fmri_dicom_dir', dwiConverter, 'inputDicomDirectory')

    grep = pipe.Node(interface=Function(function=readNrrdHeader,
                                        input_names=['fileName'],
                                        output_names=['slices', 'volumes']),
                                        name='nrrdGrep')
    preproc.connect(dwiConverter, 'outputVolume', grep, 'fileName')

    dicom = pipe.Node(interface=Function(function=dicomRead, input_names=['infolder'],
                                         output_names=['repTime']), name='dicomRead')
    preproc.connect(grabber, 'fmri_dicom_dir', dicom, 'infolder')

    #1
    to_3D_str = pipe.Node(interface=Function(function=strCreate,
                                             input_names=['slices', 'volumes', 'repTime'],
                                             output_names=['out_string']),
                          name='strCreate')
    preproc.connect(grep, 'slices', to_3D_str, 'slices')                 #1
    preproc.connect(grep, 'volumes', to_3D_str, 'volumes')
    preproc.connect(dicom, 'repTime', to_3D_str, 'repTime')

    to_3D = pipe.Node(interface=To3D(), name='afniTo3D')
    to_3D.inputs.outputtype = 'AFNI'
    to_3D.inputs.datum = 'short'
    to_3D.inputs.filetype = 'epan'
    to_3D.inputs.prefix = 'to_3D_out'
    preproc.connect(grabber, 'fmri_dicom_dir', to_3D, 'infolder')
    preproc.connect(to_3D_str, 'out_string', to_3D, 'funcparams')

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
    preproc.connect(grep, ('volumes', strToIntMinusOne), despike, 'end')

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
    warpT1 = pipe.Node(interface=ApplyTransforms(), iterfield=['input_image'],
                          name='antsApplyTransformsT1')
    warpT1.inputs.interpolation='Linear' # Validation test value: NearestNeighbor
    preproc.connect(convertT1, 'out_file', warpT1, 'input_image') # connected to brain.nii NOT brain.mgz
    preproc.connect(tstat, 'out_file', warpT1, 'reference_image')
    preproc.connect([(bFit, warpT1, [(('outputTransform', makeList), 'transforms')])])

    # warpT1 = warpNN.clone('antsApplyTransformT1')
    warpCSF = warpT1.clone('antsApplyTransformCSF')
    preproc.connect(grabber, 'csfFile', warpCSF, 'input_image')
    preproc.connect(tstat, 'out_file', warpCSF, 'reference_image')
    preproc.connect([(bFit, warpCSF, [(('outputTransform', makeList), 'transforms')])])

    warpWHM = warpT1.clone('antsApplyTransformWHM')
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
        ### Create seed points
        talairach, nac = getAtlasPoints('seeds.txt')
        nacAtlasFile = "/IPLlinux/ipldev/scratch/welchdm/bld/rs-fMRI/BSA/ReferenceAtlas-build/Atlas/Atlas_20130106/template_t1.nii.gz"
        points = pipe.Node(interface=Function(function=createSphereExpression, input_names=['coordinates', 'radius'],
                                              output_names=['expression']), name='createSphereExpression')
        points.iterables = ('coordinates', nac)

        spheres = pipe.Node(interface=Calc(letters=['a']), name='afni3Dcalc_seeds')
        spheres.inputs.outputtype = outputType
        spheres.inputs.in_file_a = nacAtlasFile
        spheres.inputs.args = '-nscale'

        preproc.connect(points, 'expression', spheres, 'expr')

        seed_DataSink = pipe.Node(interface=DataSink(), name="seed_DataSink")
        seed_DataSink.inputs.base_directory = '/paulsen/Experiments/20130417_rsfMRI_Results'

        # Warp seed output to detrend
        transformGrabber = pipe.Node(interface=DataGrabber(infields=['session_id'],
                                                           outfields=['atlasToSessionTransform',
                                                                      'sessionToAtlasTransform']),
                                     name='transformGrabber')
        transformGrabber.inputs.base_directory = '/paulsen/Experiments/20130202_PREDICTHD_Results/SubjectToAtlasWarped'
        Transformgrabber.inputs.template = '*'
        transformRegex = '%s/AtlasToSubject_%sComposite.h5'
        transformGrabber.inputs.field_template = dict(atlasToSessionTransform=transformRegex,
                                                      sessionToAtlasTransform=transformRegex)
        transformGrabber.inputs.template_args = dict(atlasToSessionTransform=[['session_id', '']],
                                                     sessionToAtlasTransform=[['session_id', 'Inverse']])

        preproc.connect(sessions, 'session_id', transformGrabber, 'session_id')

        nacToFMRI = pipe.Node(interface=ApplyTransforms(), name="antsApplyTransformsNAC")
        nacToFMRI.inputs.interpolation = 'Linear'

        preproc.connect(spheres, 'out_file', nacToFMRI, 'input_image')
        preproc.connect(warpT1, 'output_image', nacToFMRI, 'reference_image')

        # Apply two transforms: NAC -> T1 -> fMRI
        concatTransformsNode = pipe.Node(interface=Function(function=concatTransforms,
                                                            input_names=['transform1', 'transform2'],
                                                            output_names=['transformList']),
                                         name='concatTransforms')
        preproc.connect(bFit, 'outputTransform', concatTransformsNode, 'transform1')
        preproc.connect(transformGrabber, 'atlasToSessionTransform', concatTransformsNode, 'transform2')
        preproc.connect(concatTransformsNode, 'transformList', nacToFMRI, 'transforms')

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

        ### Move z values back into NAC atlas space

        brainsFitFMRItoT1 = bFit.clone(name='brainsFitFMRItoT1')
        warpT1toAtlas = warpT1.clone(name='warpT1toAtlas')

        preproc.connect(convertT1, 'out_file', brainsFitFMRItoT1, 'fixedVolume')
        preproc.connect(regionLogCalc, 'out_file', brainsFitFMRItoT1, 'movingVolume')
        preproc.connect(brainsFitFMRItoT1, 'outputVolume', warpT1toAtlas, 'input_image')
        warpT1toAtlas.inputs.reference_image = nacAtlasFile
        preproc.connect([transformGrabber, warpT1toAtlas],
                        [(('sessionToAtlasTransform', makeList), 'transforms')])


        ## def createStr(value):
        ##     return str(value)

        preproc.connect(sessions, ('session_id', createStr), sinker, 'container')
        preproc.connect(regionLogCalc, 'out_file', sinker, 'aparc.@logCalc')

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

        ###TODO: DataSink
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
    preproc.write_hierarchical_dotfile(dotfilename='dave.dot')
    #preproc.run()

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
