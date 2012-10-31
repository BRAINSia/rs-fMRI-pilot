#!/usr/bin/env python
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

import argparse
import os
import sys

### HACK: Set PYTHONPATH, DYLD_LIBRARY_PATH, FREESURFER_HOME, and SUBJECTS_DIR for Athena ###
os.environ['FREESURFER_HOME'] = '/opt/freesurfer'
os.environ['SUBJECTS_DIR'] = '/paulsen/MRx'
try:
    old_path = os.environ['DYLD_LIBRARY_PATH']
except KeyError:
    old_path = ''
try:
    old_fallback = os.environ['DYLD_FALLBACK_LIBRARY_PATH']
except KeyError:
    old_fallback = ''
os.environ['DYLD_LIBRARY_PATH'] = '/Volumes/scratch/welchdm/local/lib:' + old_path
os.environ['DYLD_FALLBACK_LIBRARY_PATH'] = '/opt/afni-OSX_10.7:' + old_fallback
sys.path.insert(1, '/Volumes/scratch/welchdm/src/nipype')
### END HACK ###

import nipype.pipeline.engine as pipe
from nipype.interfaces.io import DataSink, DataGrabber
from nipype.interfaces.utility import Function, IdentityInterface
from nipype.interfaces.afni.preprocess import *
from nipype.interfaces.freesurfer.preprocess import *

from convert import *
from nodes import *  ### TODO!!!

def readNrrdHeader(fileName):
    """
    Grep out the number of slices and volumes from the Nrrd file created by DWIConvert
    """
    import re
    fID = open(fileName, 'r')
    try:
        for line in fID.readlines():
            search = re.search('(?<=sizes:\s)(?:[0-9]*)\s(?:[0-9])*\s(?P<slices>[0-9]*)\s(?P<volumes>[0-9]*)', line)
            if search is not None:
                slices = search.groupdict()['slices']
                volumes = search.groupdict()['volumes']
                return slices, volumes
    finally:
        fID.close()
    raise IOError('Nrrd file %s has regex match for slice and volume numbers!' % fileName)

def strCreate(slices, volumes, repTime):
    # return "-time:zt %s %s %s alt-z" % (slices, volumes, repTime)
    return "-time:zt %s %s %s FROM_IMAGE" % (slices, volumes, repTime)

def dicomRead(infolder):
    """
    Use PyDicom to get the DICOM repetition time
    """
    import os
    import dicom
    infolder = os.path.abspath(infolder)
    for ii in os.listdir(infolder):
        meta = dicom.read_file(os.path.join(infolder, ii))
        if "RepetitionTime" in meta:
            return meta.data_element('RepetitionTime').value
    raise Exception('None of the dicom files in %s contain \
                     "Repetition Time" where PyDicom can find it!' % infolder)

def generateTissueMask(input_file, low=0.0, high=1.0, erodeFlag=False):
    """
    Using posterior tissue probability file, threshold by
    a low and/or high value, erode by 2mm, and take the ceil()
    Return the result as a binary mask file
    """
    import os
    import SimpleITK as sitk
    reader = sitk.ImageFileReader()
    reader.SetFileName(input_file)
    image = reader.Execute()
    ## Compute brain mask
    binaryMask = sitk.BinaryThreshold(image, low, high, 1, 0)
    if erodeFlag:
        radiusMM = 2
        kernel = sitk.BinaryErodeImageFilter() # kernel = sitk.ErodeObjectMorphologyImageFilter()
        kernel.SetKernelType(kernel.Ball)
        kernel.SetKernelRadius(radiusMM)
        erodedMask = sitk.ErodeObjectMorphology(binaryMask, kernel.GetKernelType())
        connected = sitk.ConnectedComponent(erodedMask)
    else:
        connected = sitk.ConnectedComponent(binaryMask)
    relabeled = sitk.RelabelComponent(connected)
    maskOnly = sitk.BinaryThreshold(relabeled, 1.0, 1.0, 1, 0)
    # maskOnly = getLargestConnectedRegion(erodedMask)
    writer = sitk.ImageFileWriter()
    writer.SetFileName(os.path.join(os.getcwd(), 'maskOnly.nii.gz'))
    writer.Execute(maskOnly)
    return writer.GetFileName()

def pipeline(args):
    sessionID = args.sessionID
    outputType = args.outputType
    fOutputType = args.fOutputType
    preproc = pipe.Workflow(name='rs_fmri_preprocessing')
    preproc.base_dir = os.getcwd()
    #--------------------------------------------------------------------------------------
    sessions = pipe.Node(interface=IdentityInterface(fields=['session_id']), name='sessionIDs')
    sessions.iterables = ('session_id', sessionID)
    #--------------------------------------------------------------------------------------
    grabber = pipe.Node(interface=DataGrabber(infields=['session_id'],
                                              outfields=['fmri_dicom_dir', 't1_file',
                                                         'csf_file', 'wm_file']),
                        name='dataGrabber')
    grabber.inputs.base_directory = '/paulsen'
    grabber.inputs.template = '*'
    grabber.inputs.field_template = dict(fmri_dicom_dir='MRx/FMRI_HD_120/*/%s/%s/%s/*',
                                         t1_file='Experiments/20120722_JOY_DWI/FMRI_HD_120/*/%s/%s/%s_*_%s_FS/mri/T1.mgz',
                                         csf_file='Experiments/20120801.SubjectOrganized_Results/FMRI_HD_120/*/%s/%s/%s.nii.gz',
                                         wm_file='Experiments/20120801.SubjectOrganized_Results/FMRI_HD_120/*/%s/%s/%s.nii.gz')
    grabber.inputs.template_args = dict(fmri_dicom_dir=[['session_id', 'ANONRAW', 'FMRI_RestingStateConnectivity']],
                                        t1_file=[['session_id', '10_AUTO.NN3Tv20110419','JOY_v51_2011', 'session_id']],
                                        csf_file=[['session_id','ACCUMULATED_POSTERIORS', 'POSTERIOR_CSF_TOTAL']],
                                        wm_file=[['session_id','ACCUMULATED_POSTERIORS', 'POSTERIOR_WM_TOTAL']])
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
    to_3D.inputs.datum = 'float'
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
    despike.inputs.in_file = os.path.join(os.getcwd(), 'rs_fmri_preprocessing',
                                          'afni3Drefit', 'to_3D_out+orig.BRIK') # TODO
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
    tstat = pipe.Node(interface=TStat(), name='afni3DtStat')
    tstat.inputs.outputtype = outputType
    tstat.inputs.args = '-mean' # TODO
    #--------------------------------------------------------------------------------------
    calc = pipe.Node(interface=Calc(), name='afni3Dcalc')
    calc.inputs.outputtype = outputType
    calc.inputs.expr = "'(a - b) + 1000'"
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
    #--------------------------------------------------------------------------------------
    skullstrip = pipe.Node(interface=SkullStrip(), name='afni3DskullStrip')
    skullstrip.inputs.outputtype = outputType
    #--------------------------------------------------------------------------------------
    register1 = pipe.Node(interface=Allineate(), name='afni3Dallineate1')
    register1.inputs.outputtype = outputType
    register1.inputs.warp = 'shift_rotate'
    register1.inputs.cost = 'mutualinfo'
    register1.inputs.cmass = True
    register1.inputs.interp = 'triquintic'
    register1.inputs.final = 'triquintic'
    register1.inputs.onedmatrix_save = True
    ### HACK
    register1.inputs.out_file = 'register1.nii'
    ### END HACK
    #--------------------------------------------------------------------------------------
    register2 = pipe.Node(interface=Allineate(), name='afni3Dallineate2')
    register2.inputs.outputtype = outputType
    register2.inputs.final = 'triquintic'
    ### HACK
    register2.inputs.out_file = 'register2.nii'
    ### END HACK
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
    wmmask.inputs.low = 1.0
    wmmask.inputs.erodeFlag = True
    #--------------------------------------------------------------------------------------
    csfAvg = pipe.Node(interface=Maskave(), name='afni3DmaskAve_csf')
    csfAvg.inputs.outputtype = outputType
    csfAvg.inputs.args = '-median -mrange 1 1' # TODO
    csfAvg.inputs.quiet = True
    #--------------------------------------------------------------------------------------
    wmAvg = pipe.Node(interface=Maskave(), name='afni3DmaskAve_wm')
    wmAvg.inputs.outputtype = outputType
    wmAvg.inputs.args = '-median -mrange 1 2' # TODO
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
    #--------------------------------------------------------------------------------------
    detrend = pipe.Node(interface=Detrend(), name='afni3Ddetrend')
    detrend.inputs.outputtype = outputType
    detrend.inputs.out_file = 'errts_Decon_dt+tlrc'
    detrend.inputs.args = '-polort 3' # TODO
    #--------------------------------------------------------------------------------------
    regionAvg = pipe.Node(interface=Maskave(), name='afni3DmaskAve_region')
    regionAvg.inputs.outputtype = outputType
    regionAvg.inputs.args = '-median -mrange 1 1' # TODO
    regionAvg.inputs.quiet = True
    #--------------------------------------------------------------------------------------
    # Connect pipeline
    if len(sessionID) > 1:
        preproc.connect(sessions, 'session_id', grabber, 'session_id')
    else:
        grabber.inputs.session_id = sessionID[0]
    #1.
    preproc.connect(grabber, 'fmri_dicom_dir', to_3D, 'infolder')
    preproc.connect(grabber, 'fmri_dicom_dir', dicom, 'infolder')
    preproc.connect(grabber, 'fmri_dicom_dir', dicomNrrd, 'inputDicomDirectory')
    preproc.connect(dicomNrrd, 'outputVolume', grep, 'fileName')
    preproc.connect(grep, 'slices', to_3D_str, 'slices')
    preproc.connect(grep, 'volumes', to_3D_str, 'volumes')
    preproc.connect(dicom, 'repTime', to_3D_str, 'repTime')
    preproc.connect(to_3D_str, 'out_string', to_3D, 'funcparams')
    preproc.connect(to_3D, 'out_file', refit, 'in_file')
    #2.
    preproc.connect(refit, 'out_file', despike, 'in_file')
    #3.
    preproc.connect(despike, 'out_file', volreg, 'in_file')
    #4.
    preproc.connect(volreg, 'out_file', tstat, 'in_file')
    preproc.connect(volreg, 'out_file', calc, 'in_file_a')
    preproc.connect(tstat, 'out_file', calc, 'in_file_b')
    #5.
    preproc.connect(calc, 'out_file', fourier, 'in_file')
    #6.
    preproc.connect(fourier, 'out_file', merge, 'in_files')
    #7.
    preproc.connect(grabber, 't1_file', converter, 'in_file')
    preproc.connect(converter, 'out_file', skullstrip, 'in_file')
    preproc.connect(skullstrip, 'out_file', register1, 'base_file')
    # preproc.connect(converter, 'out_file', register1, 'master_file')
    preproc.connect(calc, 'out_file', register1, 'in_file')
    preproc.connect(register1, 'onedmatrix_out', register2, 'onedmatrix')
    preproc.connect(fourier, 'out_file', automask, 'in_file')
    preproc.connect(automask, 'out_file', register2, 'in_file')
    #8.
    preproc.connect(grabber, 'csf_file', csfmask, 'input_file')
    preproc.connect(grabber, 'wm_file', wmmask, 'input_file')
    #9.
    preproc.connect(csfmask, 'output_file', csfAvg, 'mask')
    preproc.connect(fourier, 'out_file', csfAvg, 'in_file')
    preproc.connect(wmmask, 'output_file', wmAvg, 'mask')
    preproc.connect(fourier, 'out_file', wmAvg, 'in_file')
    preproc.connect(fourier, 'out_file', deconvolve, 'in_file')
    preproc.connect(automask, 'out_file', deconvolve, 'mask') ### VERIFY WITH JATIN!!!
    preproc.connect(csfAvg, 'out_file', deconvolve, 'stim_file_1')
    deconvolve.inputs.stim_label_1 = "Median_CSF"
    preproc.connect(wmAvg, 'out_file', deconvolve, 'stim_file_2')
    deconvolve.inputs.stim_label_2 = "Median_WM"
    preproc.connect(volreg, 'oned_file', deconvolve, 'stim_file_3')
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
    preproc.connect(deconvolve, 'out_errts', detrend, 'in_file')
    ### TODO: DataSink
    preproc.write_graph()
    preproc.write_hierarchical_dotfile(dotfilename='dave.dot')
    # preproc.run()

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
    os.environ['DYLD_LIBRARY_PATH'] = old_path
    os.environ['DYLD_FALLBACK_LIBRARY_PATH'] = old_fallback
    sys.path.pop(1)
    ### END HACK ###
    sys.exit(outvalue)
