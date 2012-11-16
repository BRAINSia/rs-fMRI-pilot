#!/usr/bin/env python
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

import argparse
import os
import sys

### HACK: Set PYTHONPATH, DYLD_LIBRARY_PATH, FREESURFER_HOME, and SUBJECTS_DIR for Athena ###
### export PATH=$PATH:/opt/afni-OSX_10.7:/opt/freesurfer_v4.5.0-full/bin ###
os.environ['FREESURFER_HOME'] = '/opt/freesurfer_v4.5.0-full'
os.environ['SUBJECTS_DIR'] = '/paulsen/MRx'
try:
    old_dyld = os.environ['DYLD_LIBRARY_PATH']
except KeyError:
    old_dyld = ''
old_path = os.environ['PATH']
try:
    old_fallback = os.environ['DYLD_FALLBACK_LIBRARY_PATH']
except KeyError:
    old_fallback = ''

os.environ['PATH'] = '/Volumes/scratch/welchdm/bld/BSA-20121108/bin:' + \
    '/opt/afni-OSX_10.7:' + '/opt/freesurfer_v4.5.0-full/bin:' + old_path

os.environ['DYLD_LIBRARY_PATH'] = '/Volumes/scratch/welchdm/bld/BSA-20121108/lib:' + \
    '/Volumes/scratch/welchdm/bld/BSA-20121108/bin:' + \
    '/opt/freesurfer_v4.5.0-full/lib:' + old_dyld

os.environ['DYLD_FALLBACK_LIBRARY_PATH'] = '/opt/afni-OSX_10.7:' + old_fallback

sys.path.insert(1, '/Volumes/scratch/welchdm/src/nipype/nipype')
sys.path.insert(2, '/Volumes/scratch/welchdm/bld/BSA-20121108/SimpleITK-build/lib')
sys.path.insert(3, '/Volumes/scratch/welchdm/src/BRAINSStandAlone/AutoWorkup')
### END HACK ###

from SEMTools import BRAINSFit
from nipype.interfaces.ants.registration import Registration
from nipype.interfaces.ants.resampling import ApplyTransforms
from nipype.interfaces.afni.preprocess import *
from nipype.interfaces.freesurfer.preprocess import *
from nipype.interfaces.io import DataSink, DataGrabber
from nipype.interfaces.utility import Function, IdentityInterface
import nipype.pipeline.engine as pipe
import numpy

from convert import *
from utilities import *
from writeMat import *
### from nodes import *  ### TODO!!!


def pipeline(args):
    sessionID = args.sessionID
    outputType = args.outputType
    fOutputType = args.fOutputType
    preproc = pipe.Workflow(name='rs_fmri_workflow')
    preproc.base_dir = os.getcwd()
    #--------------------------------------------------------------------------------------
    sessions = pipe.Node(interface=IdentityInterface(fields=['session_id']), name='sessionIDs')
    sessions.iterables = ('session_id', sessionID)
    #--------------------------------------------------------------------------------------
    grabber = pipe.Node(interface=DataGrabber(infields=['session_id'],
                                              outfields=['fmri_dicom_dir', 't1_File', 'f1_File',
                                                         'f2_File', 'csfFile', 'whmFile']),
                        name='dataGrabber')
    grabber.inputs.base_directory = '/paulsen'
    grabber.inputs.template = '*'
    fmriRegex = 'MRx/FMRI_HD_120/*/%s/%s/%s/*'
    fS_Regex = 'Experiments/20120722_JOY_DWI/FMRI_HD_120/*/%s/%s/%s_*_%s_FS/%s/%s'
    probRegex = 'Experiments/20120801.SubjectOrganized_Results/FMRI_HD_120/*/%s/%s/%s.nii.gz'
    grabber.inputs.field_template = dict(fmri_dicom_dir=fmriRegex,
                                         t1_File=fS_Regex,
                                         f1_File=fS_Regex,
                                         f2_File=fS_Regex,
                                         csfFile=probRegex,
                                         whmFile=probRegex)
    grabber.inputs.template_args = dict(fmri_dicom_dir=[['session_id', 'ANONRAW',
                                                         'FMRI_RestingStateConnectivity']],
                                        t1_File=[['session_id', '10_AUTO.NN3Tv20110419',
                                                  'JOY_v51_2011', 'session_id', 'mri', 'brain.mgz']],
                                        f1_File=[['session_id', '10_AUTO.NN3Tv20110419',
                                                  'JOY_v51_2011', 'session_id', 'mri_nifti',
                                                  'aparc+aseg.nii.gz']],
                                        f2_File=[['session_id', '10_AUTO.NN3Tv20110419',
                                                  'JOY_v51_2011', 'session_id', 'mri_nifti',
                                                  'aparc.a2009s+aseg.nii.gz']],
                                        csfFile=[['session_id', 'ACCUMULATED_POSTERIORS',
                                                   'POSTERIOR_CSF_TOTAL']],
                                        whmFile=[['session_id', 'ACCUMULATED_POSTERIORS',
                                                  'POSTERIOR_WM_TOTAL']])
    #--------------------------------------------------------------------------------------
    dicomNrrd = pipe.Node(interface=DWIconvert(), name='dicomToNrrd')
    dicomNrrd.inputs.conversionMode = 'DicomToNrrd'
    dicomNrrd.inputs.outputVolume = 'converted.nrrd'
    #--------------------------------------------------------------------------------------
    grep = pipe.Node(interface=Function(function=readNrrdHeader,
                                        input_names=['fileName'],
                                        output_names=['slices', 'volumes']),
                                        name='nrrdGrep')
    #--------------------------------------------------------------------------------------
    to_3D_str = pipe.Node(interface=Function(function=strCreate,
                                             input_names=['slices', 'volumes', 'repTime'],
                                             output_names=['out_string']),
                                             name='strCreate')
    #--------------------------------------------------------------------------------------
    dicom = pipe.Node(interface=Function(function=dicomRead, input_names=['infolder'],
                                         output_names=['repTime']), name='dicomRead')
    #--------------------------------------------------------------------------------------
    to_3D = pipe.Node(interface=To3D(), name='afniTo3D')
    to_3D.inputs.outputtype = 'AFNI'
    to_3D.inputs.datum = 'short'
    to_3D.inputs.filetype = 'epan'
    to_3D.inputs.prefix = 'to_3D_out'
    #--------------------------------------------------------------------------------------
    refit = pipe.Node(interface=Refit(), name='afni3Drefit')
    refit.inputs.outputtype = 'AFNI'
    refit.inputs.deoblique = True
    #--------------------------------------------------------------------------------------
    despike = pipe.Node(interface=Despike(), name='afni3Ddespike')
    despike.inputs.outputtype = outputType
    despike.inputs.ignore = 4
    #--------------------------------------------------------------------------------------
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
    #--------------------------------------------------------------------------------------
    zpad = pipe.Node(interface=Zeropad(), name='afniZeropad')
    zpad.inputs.plane = 'IS'
    zpad.inputs.numberOfPlanes = 44
    zpad.inputs.is_mm = False
    #--------------------------------------------------------------------------------------
    tstat = pipe.Node(interface=TStat(), name='afni3DtStat')
    tstat.inputs.outputtype = outputType
    tstat.inputs.args = '-mean' # TODO
    #--------------------------------------------------------------------------------------
    calc = pipe.Node(interface=Calc(letters=['a', 'b','c']), name='afni3Dcalc')
    calc.inputs.outputtype = outputType
    calc.inputs.expr = "((a - b) + 1000) * c"
    #--------------------------------------------------------------------------------------
    fourier = pipe.Node(interface=Fourier(), name='afni3Dfourier')
    fourier.inputs.outputtype = outputType
    fourier.inputs.highpass = 0.011
    fourier.inputs.lowpass = 0.1
    fourier.inputs.args = '-retrend' # TODO
    fourier.inputs.outputtype = outputType
    #--------------------------------------------------------------------------------------
    merge = pipe.Node(interface=Merge(), name='afni3Dmerge')
    merge.inputs.outputtype = outputType
    merge.inputs.blurfwhm = 6
    merge.inputs.doall = True
    merge.inputs.args = '-1noneg -1clip 100' # TODO
    # TODO: implement
    # merge.inputs.onenoneg = True
    # merge.inputs.oneclip = 100
    #--------------------------------------------------------------------------------------
    automask = pipe.Node(interface=Automask(), name='afni3Dautomask')
    automask.inputs.outputtype = outputType
    automask.inputs.dilate = 1
    #--------------------------------------------------------------------------------------
    converter = pipe.Node(interface=MRIConvert(), name='freesurferMRIconvert')
    converter.inputs.out_type = fOutputType
    converter.inputs.in_type = 'mgz'
    converter.inputs.force_ras = True

    convertT1 = converter.clone('T1Converter')
    convertF1 = converter.clone('F1Converter')
    convertF2 = converter.clone('F2Converter')
    #--------------------------------------------------------------------------------------
    bFit = pipe.Node(interface=BRAINSFit(), name='brainsFit')
    bFit.inputs.initializeTransformMode = 'useCenterOfHeadAlign'
    bFit.inputs.maskProcessingMode = 'ROIAUTO'
    bFit.inputs.ROIAutoDilateSize = 10.0
    bFit.inputs.useRigid = True
    bFit.inputs.costMetric = 'MMI' # (default)
    bFit.inputs.numberOfSamples = 100000 # (default)
    bFit.inputs.outputTransform = True
    #--------------------------------------------------------------------------------------
    warpT1 = pipe.MapNode(interface=ApplyTransforms(), iterfield=['input_image'],
                          name='antsApplyTransformsT1')
    warpT1.inputs.interpolation='NearestNeighbor'

    # warpT1 = warpNN.clone('antsApplyTransformT1')
    warpCSF = warpT1.clone('antsApplyTransformCSF')
    warpWHM = warpT1.clone('antsApplyTransformWHM')
    #--------------------------------------------------------------------------------------
    warpFS1 = pipe.Node(interface=ApplyTransforms(), name='antsApplyTransformsFS1')
    warpFS1.inputs.interpolation='MultiLabel'

    warpFS2 = pipe.Node(interface=ApplyTransforms(), name='antsApplyTransformsFS2')
    warpFS2.inputs.interpolation='MultiLabel'
    #--------------------------------------------------------------------------------------
    csfmask = pipe.Node(interface=Function(function=generateTissueMask,
                                           input_names=['input_file','low', 'high', 'erodeFlag'],
                                           output_names=['output_file']),
                        name='csfMask')
    csfmask.inputs.low = 0.99
    csfmask.inputs.high = 1.0
    csfmask.inputs.erodeFlag = False
    #--------------------------------------------------------------------------------------
    wmmask = pipe.Node(interface=Function(function=generateTissueMask,
                                           input_names=['input_file','low', 'high', 'erodeFlag'],
                                           output_names=['output_file']),
                        name='wmMask')
    wmmask.inputs.low = 0.99
    wmmask.inputs.high = 1.0
    wmmask.inputs.erodeFlag = True
    #--------------------------------------------------------------------------------------
    csfAvg = pipe.Node(interface=Maskave(), name='afni3DmaskAve_csf')
    csfAvg.inputs.outputtype = outputType
    csfAvg.inputs.args = '-median' # TODO
    csfAvg.inputs.quiet = True
    #--------------------------------------------------------------------------------------
    wmAvg = pipe.Node(interface=Maskave(), name='afni3DmaskAve_wm')
    wmAvg.inputs.outputtype = outputType
    wmAvg.inputs.args = '-median' # TODO
    wmAvg.inputs.quiet = True
    #--------------------------------------------------------------------------------------
    deconvolve = pipe.Node(interface=Deconvolve(fileCount=3, seriesCount=8), name='afni3Ddeconvolve')
    deconvolve.inputs.outputtype = outputType
    deconvolve.inputs.ignoreWarnings = 4
    deconvolve.inputs.nullHypothesisPolynomialDegree = 1
    deconvolve.inputs.full_first = True
    deconvolve.inputs.is_float = True
    deconvolve.inputs.tout = True
    deconvolve.inputs.rout = True
    deconvolve.inputs.fout = True
    deconvolve.inputs.bucket = 'Rest_bp_Decon+tlr'
    deconvolve.inputs.fitts = 'full_fitts_Decon+tlrc'
    deconvolve.inputs.errts = 'errts_Decon+tlrc'
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
    #--------------------------------------------------------------------------------------
    detrend = pipe.Node(interface=Detrend(), name='afni3Ddetrend')
    detrend.inputs.outputtype = outputType
    detrend.inputs.out_file = 'errts_Decon_dt+tlrc'
    detrend.inputs.args = '-polort 3' # TODO
    #--------------------------------------------------------------------------------------
    # TODO: Iterate over this node using the args = '-median -mrange %d %d' function...
    regionAvg = pipe.Node(interface=Maskave(), name='afni3DmaskAve_region')
    regionAvg.inputs.outputtype = outputType
    regionAvg.inputs.args = '-median' # TODO
    regionAvg.inputs.quiet = True
    #--------------------------------------------------------------------------------------
    # writeFile = pipe.Node(interface=Function(function=writeMat,
    #     input_names=['data', 'fileName', 'key']),
    #                      name='writeMatFile')
    #--------------------------------------------------------------------------------------
    # Connect pipeline
    if len(sessionID) > 1:
        preproc.connect(sessions, 'session_id', grabber, 'session_id')
    else:
        grabber.inputs.session_id = sessionID[0]
    preproc.connect(grabber, 'fmri_dicom_dir', to_3D, 'infolder')
    preproc.connect(grabber, 'fmri_dicom_dir', dicom, 'infolder')
    preproc.connect(grabber, 'fmri_dicom_dir', dicomNrrd, 'inputDicomDirectory')
    preproc.connect(dicomNrrd, 'outputVolume', grep, 'fileName')
    ### preprocessing_Nov1_update.sh
    preproc.connect(grep, 'slices', to_3D_str, 'slices')                 #1
    preproc.connect(grep, 'volumes', to_3D_str, 'volumes')
    preproc.connect(dicom, 'repTime', to_3D_str, 'repTime')
    preproc.connect(to_3D_str, 'out_string', to_3D, 'funcparams')
    preproc.connect(to_3D, 'out_file', refit, 'in_file')                  #1a
    preproc.connect(refit, 'out_file', despike, 'in_file')                #2
    preproc.connect(despike, 'out_file', volreg, 'in_file')               #3
    preproc.connect(volreg, 'out_file', zpad, 'in_file')                  #4
    preproc.connect(zpad, 'out_file', merge, 'in_files')                  #5 ### NOTE: ONLY ONE FILE TO in_files
    preproc.connect(merge, 'out_file', automask, 'in_file')               #6
    preproc.connect(merge, 'out_file', tstat, 'in_file')                  #7 ### 'mean_file' -> 'in_file'
    preproc.connect(merge, 'out_file', calc, 'in_file_a')
    preproc.connect(automask, 'out_file', tstat, 'mask_file')
    preproc.connect(tstat, 'out_file', calc, 'in_file_b')
    preproc.connect(automask, 'out_file', calc, 'in_file_c')
    preproc.connect(calc, 'out_file', fourier, 'in_file')                 #8
    #### Freesurfer section
    preproc.connect(grabber, 't1_File', convertT1, 'in_file')             #9
    preproc.connect(convertT1, 'out_file', bFit, 'movingVolume')
    preproc.connect(tstat, 'out_file', bFit, 'fixedVolume')
    ### preprocessing_part2.sh
    ### These will change to Linear later...
    preproc.connect([(bFit, warpT1, [(('outputTransform', makeList), 'transforms')])])
    preproc.connect(convertT1, 'out_file', warpT1, 'input_image') # connected to brain.nii NOT brain.mgz
    preproc.connect(tstat, 'out_file', warpT1, 'reference_image')

    preproc.connect([(bFit, warpCSF, [(('outputTransform', makeList), 'transforms')])])
    preproc.connect(grabber, 'csfFile', warpCSF, 'input_image')
    preproc.connect(tstat, 'out_file', warpCSF, 'reference_image')
    preproc.connect([(bFit, warpWHM, [(('outputTransform', makeList), 'transforms')])])
    preproc.connect(grabber, 'whmFile', warpWHM, 'input_image')
    preproc.connect(tstat, 'out_file', warpWHM, 'reference_image')

    preproc.connect(warpCSF, 'output_image', csfmask, 'input_file')       #10
    preproc.connect(warpWHM, 'output_image', wmmask, 'input_file')        #11
    ### TODO: Extract ROIs
    ###         1) Register labels with fMRI
    ## Register multilabel files
    preproc.connect([(bFit, warpFS1, [(('outputTransform', makeList), 'transforms')])])
    preproc.connect(grabber, 'f1_File', warpFS1, 'input_image')
    preproc.connect(tstat, 'out_file', warpFS1, 'reference_image')
    preproc.connect([(bFit, warpFS2, [(('outputTransform', makeList), 'transforms')])])
    preproc.connect(grabber, 'f2_File', warpFS2, 'input_image')
    preproc.connect(tstat, 'out_file', warpFS2, 'reference_image')

    preproc.connect(csfmask, 'output_file', csfAvg, 'mask')
    preproc.connect(tstat, 'out_file', csfAvg, 'in_file')
    preproc.connect(wmmask, 'output_file', wmAvg, 'mask')
    preproc.connect(tstat, 'out_file', wmAvg, 'in_file')
    preproc.connect(fourier, 'out_file', deconvolve, 'in_file')           #12
    preproc.connect(csfAvg, 'out_file', deconvolve, 'stim_file_1')
    preproc.connect(wmAvg, 'out_file', deconvolve, 'stim_file_2')
    preproc.connect(volreg, 'oned_file', deconvolve, 'stim_file_3')
    preproc.connect(deconvolve, 'out_errts', detrend, 'in_file')          #13
    ###         2) Binary threshold label image for each label
    ###         3) Calculate the average fMRI value per label
    ###         4) Output label and average to a dictionary?
    ### TODO: create NxN matrix files
    ### TODO: DataSink
    preproc.write_graph()
    preproc.write_hierarchical_dotfile(dotfilename='dave.dot')
    preproc.run()

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Preprocessing script for resting state fMRI')
    parser.add_argument('-t', '--filetype', action='store', dest='outputType', type=str, required=False,
                        default='NIFTI', help='File type for outputs.  Values: "NIFTI_GZ", "AFNI", "NIFTI" (default)')
    parser.add_argument('-n','--name', action='store', dest='name', type=str, required=True,
                        help='Name (required)')
    parser.add_argument('-s', '--session', action='store', dest='sessionID', type=int, required=True,
                        nargs='*', help='list of session IDs (required)')
    args = parser.parse_args()
    freesurferOutputTypes = {"NIFTI_GZ" : "niigz",
                             "AFNI" : "afni",
                             "NIFTI" : "nii"}
    args.fOutputType = freesurferOutputTypes[args.outputType]
    outvalue = pipeline(args)
    ### HACK ###
    os.environ['FREESURFER_HOME'] = ''
    os.environ['SUBJECTS_DIR'] = ''
    os.environ['PATH'] = old_path
    os.environ['DYLD_LIBRARY_PATH'] = old_dyld
    os.environ['DYLD_FALLBACK_LIBRARY_PATH'] = old_fallback
    sys.path.pop(3)
    sys.path.pop(2)
    sys.path.pop(1)
    ### END HACK ###
    sys.exit(outvalue)
