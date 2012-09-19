#!/usr/bin/env python
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

import argparse
import os
import sys

import nipype.pipeline.engine as pipe
from nipype.interfaces.io import DataSink, DataGrabber
from nipype.interfaces.utility import Function, IdentityInterface
from nipype.interfaces.afni.preprocess import *
from nipype.interfaces.freesurfer.preprocess import *

from .convert import *

def readNrrdHeader(fileName):
    """
    Grep out the number of slices and volumes from the Nrrd file created by DWIConvert
    """
    import re
    fID = open(fileName, 'r')
    try:
        for line in fID.readlines():
            search = re.search('(?<=sizes:\s)(?:[0-9]*)\s(?:[0-9])*\s(?P<slices>[0-9]*)\s(?P<volumes>[0-9]*)')
            if search is not None:
                slices = search.groupdict()['slices']
                volumes = search.groupdict()['volumes']
                return slices, volumes
    finally:
        fID.close()
    raise IOError('Nrrd file %s has regex match for slice and volume numbers!' % fileName)

def strCreate(slices, volumes, repTime):
    return "-time:zt %s %s %d alt-z" % (slices, volumes, repTime)

def dicomRead(infolder):
    """
    Use PyDicom to get the DICOM repetition time
    """
    import os
    import dicom
    for ii in os.listdir(os.path.abspath(infolder)):
        meta = dicom.read_file(ii)
        if "RepetitionTime" in meta:
            return meta.data_element('RepetitionTime')
    raise Exception('None of the dicom files in %s contain \
                     "Repetition Time" where PyDicom can find it!' % infolder)

def generateTissueMask(inputLabelMap, low=0.0, high=1.0):
    """
    Using posterior tissue probability file, threshold by
    a low and/or high value, erode by 2mm, and take the ceil()
    Return the result as a binary mask file
    """
    import SimpleITK as sitk

    def getLargestConnectedRegion(mask):
        """
        Labels connected regions in order from largest (#1) to smallest,
        then thresholds to only return the largest region
        """
        connected = sitk.ConnectedComponent(mask)
        relabeled = sitk.RelabelComponent(connected)
        largestLabel = sitk.BinaryThreshold(relabeled, 1, 1, 1, 0)
        return largestLabel

    ## Compute brain mask
    binaryMask = sitk.BinaryThreshold(input_file, low, high)
    kernel = sitk.ErodeObjectMorphology.KernelType('Ball')
    radiusMM = 2
    # Does this take the ceil() of the erosion?
    erodedMask = sitk.ErodeObjectMorphology(binaryMask, radiusMM, kernel, 1, 0)
    maskOnly = getLargestConnectedRegion(erodedMask)
    return maskOnly

def pipeline(**kwds):
    preproc = pipe.Workflow()

    # filePaths = glob.glob('/paulsen/MRx/FMRI_HD_*/*/%s/ANONRAW/FMRI_RestingStateConnectivity/*/*.IMA' % sessionID)
    # if len(filePath) == 0:
    #     filePaths = glob.glob('/paulsen/MRx/FMRI_HD_*/*/%s/ANONRAW/FMRI_RestingStateConnectivity/*/*.dcm' % sessionID)
    # header = dicom.read_file(filePaths[0])
    # repetitionTime = header.RepetitionTime

    sessions = pipe.Node(interface=IdentityInterface(fields=['session_id']), name='sessionIDs')
    sessions.iterables = ('session_id', sessionID)

    grabber = pipe.Node(interface=DataGrabber(infields=['session_id'],
                                              outfields=['fmri_dicom_dir', 't1_file',
                                                         'csf_file', 'wm_file']),
        name='dataGrabber')
    grabber.inputs.base_directory = '/paulsen'
    grabber.inputs.template = '*'
    grabber.inputs.field_template = dict(fmri_dicom_dir='MRx/FMRI_HD_120/*/%s/%s/%s/*'
                                         t1_file='Experiments/20120722_JOY_DWI/FMRI_HD_120/*/%s/%s/%s_*_%s/T1.mgz'
                                         csf_file='Experiments/20120801.SubjectOrganized_Results/PHD_120/*/%s/%s/%s.nii.gz',
                                         wm_file= 'Experiments/20120801.SubjectOrganized_Results/PHD_120/*/%s/%s/%s.nii.gz')
    grabber.inputs.template_args = dict(fmri_dicom_dir=[['session_id', 'ANONRAW', 'FMRI_RestingStateConnectivity']],
                                         t1_file=[['session_id', '10_AUTO.NN3Tv20110419','JOY_v51_2011', 'session_id']],
                                         csf_file=[['session_id','ACCUMULATED_POSTERIORS', 'POSTERIOR_CSF_TOTAL']],
                                         wm_file=[['session_id','ACCUMULATED_POSTERIORS', 'POSTERIOR_WM_TOTAL']])

    # fMRIgrabber = pipe.Node(interface=DataGrabber(infields=['session_id'],
    #                                               outfields=['outfolder']),
    #                         name='fMRIdata')
    # fMRIgrabber.inputs.base_directory = '/paulsen/MRx/FMRI_HD_120'
    # fMRIgrabber.inputs.template = '*/%s/ANONRAW/FMRI_RestingStateConnectivity/*'
    # fMRIgrabber.inputs.sorted = True
    # fMRIgrabber.inputs.template_args['outfolder'] = [['session_id', 'session_id']]

    # t1grabber = pipe.Node(interface=DataGrabber(infields=['session_id'],
    #                                             outfields=['out_file']),
    #                       name='T1data')
    # t1grabber.inputs.base_directory = '/paulsen/Experiments/20120722_JOY_DWI/FMRI_HD_120'
    # t1grabber.inputs.template = '*/%s/10_AUTO.NN3Tv20110419/JOY_v51_2011_*_%s/T1.mgz'
    # t1grabber.inputs.template_args['out_file'] = [['session_id', 'session_id']]

    # posteriorGrabber = pipe.Node(interface=DataGrabber(infields=['session_id', 'posteriors'],
    #                                                    outfields=['csf_file', 'wm_file']),
    #                              name='posteriorData')
    # t1grabber.inputs.base_directory = '/paulsen/Experiments/20120801.SubjectOrganized_Results/PHD_120'
    # t1grabber.inputs.template = '*/%s/ACCUMULATED_POSTERIORS/%s'
    # t1grabber.inputs.template_args['csf_file'] = [['session_id', 'POSTERIOR_CSF_TOTAL.nii.gz']]
    # t1grabber.inputs.template_args['wm_file'] = [['session_id', 'POSTERIOR_WM_TOTAL.nii.gz']]

    if len(sessionID) > 1:
        preprocess.connect(sessions, 'session_id', grabber, 'session_id')
        preprocess.connect(sessions, 'session_id', grabber, 'session_id')
    else:
        grabber.inputs.session_id = sessionID[0]
        grabber.inputs.session_id = sessionID[0]

    dicomNrrd = pipe.Node(interface=DWIconvert(), name='dicomToNrrd')
    dicomNrrd.inputs.conversionMode = 'DicomToNrrd'

    grep = pipe.Node(interface=Function(function=readNrrdHeader(),
                                        input_names=['fileName'],
                                        output_names=['slices, volumes']),
                                        name='nrrdGrep')

    strCreate = pipe.Node(interface=Function(function=strCreate(),
                                             input_names=['slices', 'volumes', 'repTime'],
                                             output_names=['out_string']),
                                             name='strCreate')

    dicom = pipe.Node(interface=Function(function=dicomRead(),
                                         input_names=['infolder'],
                                         output_names=['repTime']),
                                         name='dicomRead')

    to_3D = pipe.Node(interface=To3D(), name='afniTo3D')
    to_3D.inputs.outputtype = outputType
    to_3D.inputs.args = "-orient RAS" # TODO
    to_3D.inputs.filetype = 'epan'

    refit = pipe.Node(interface=Refit(), name='afni3Drefit')
    refit.inputs.outputtype = outputType
    refit.inputs.deoblique = True

    despike = pipe.Node(interface=Despike(), name='afni3Ddespike')
    despike.inputs.outputtype = outputType
    despike.inputs.args = '-ignore 4'

    volreg = pipe.Node(interface=Volreg(), name='afni3DvolReg')
    volreg.inputs.outputtype = outputType
    volreg.inputs.timeshift = False # 0
    volreg.inputs.zpad = 3
    volreg.inputs.args = '-cubic -maxite 50 -x_thresh 0.001 -rot_thresh 0.001 -delta 0.1 -final Fourier \
                          -twopass -twodup -coarse 2 2 -coarserot -base 9 '
    ### TODO: implement in Nipype
    # volreg.inputs.cubic = True
    # volreg.inpus.maxite = 50
    # volreg.inpus.x_thresh = 0.001
    # volreg.inpus.rot_thresh = 0.001
    # volreg.inpus.delta = 0.1
    # volreg.inpus.final = 'Fourier'
    # volreg.inpus.twopass = True
    # volreg.inpus.twodup = True
    # volreg.inpus.coarse = [2, 2]
    # volreg.inpus.coarserot = True
    # volreg.inpus.base = 9

    tstat = pipe.Node(interface=TStat(), name='afni3DtStat')
    tstat.inputs.outputtype = outputType
    tstat.inputs.args = '-mean'

    calc = pipe.Node(interface=Calc(), name='afni3Dcalc')
    calc.inputs.outputtype = outputType
    calc.inputs.expr = '(a - b) + 1000'

    fourier = pipe.Node(interface=Fourier(), name='afni3Dfourier')
    fourier.inputs.outputtype = outputType
    fourier.inputs.highpass = 0.011
    fourier.inputs.lowpass = 0.1
    fourier.inputs.args = '-retrend'
    fourier.inputs.outputtype = outputType

    merge = pipe.Node(interface=Merge(), name='afni3Dmerge')
    merge.inputs.outputtype = outputType
    merge.inputs.blurfwhm = 6
    merge.inputs.doall = True
    merge.inputs.args = '-1noneg -1clip 100'

    converter = pipe.Node(interface=MRIConvert(), name='freesurferMRIconvert')
    converter.inputs.out_type = fOutputType
    converter.inputs.in_orientation = 'RAS'
    converter.inputs.in_type = 'mgz'
    converter.inputs.out_orientation = 'RAS'

    register1 = pipe.Node(interface=Allineate(), name='afni3Dallineate1')
    register1.inputs.outputtype = outputType
    register1.inputs.warp = 'shift_rotate'
    register1.inputs.cost = 'mutualinfo'
    register1.inputs.cmass = True
    register1.inputs.interp = 'triquintic'
    register1.inputs.final = 'triquintic'

    register2 = pipe.Node(interface=Allineate(), name='afni3Dallineate2')
    register2.inputs.outputtype = outputType
    register2.inputs.final = 'triquintic'

    csfmask = pipe.Node(interface=Function(function=generateTissueMask,
                                           input_names=['input_file','low', 'high'],
                                           output_names=['output_file']),
                        name='csfMask')
    csfmask.inputs.low = 0.99

    wmmask = pipe.Node(interface=Function(function=generateTissueMask,
                                           input_names=['input_file','low', 'high'],
                                           output_names=['output_file']),
                        name='wmMask')
    wmmask.inputs.low = 0.99

    # Connect pipeline
    #1.
    preproc.connect(grabber, 'fmri_dicom_dir', to_3D, 'infolder')
    preproc.connect(grabber, 'fmri_dicom_dir', dicom, 'infolder')
    preproc.connect(grabber, 'fmri_dicom_dir', dicomNrrd, 'inputDicomDirectory')
    preproc.connect(dicomNrrd, 'outputVolume', grep, 'fileName')
    preproc.connect(grep, 'slices', strCreate, 'slices')
    preproc.connect(grep, 'volumes', strCreate, 'volumes')
    preproc.connect(dicom, 'repTime', strCreate, 'repTime')
    preproc.connect(strCreate, 'out_string', to_3D, 'funcparams')
    preproc.connect(to_3D, 'out_file', refit, 'in_file')
    #2.
    preproc.connect(to_3D, 'out_file', despike, 'in_file')
    #3.
    preproc.connect(despike, 'out_file', volreg, 'in_file')
    #4.
    preproc.connect(volreg, 'out_file', tstat, 'in_file')
    preproc.connect(volreg, 'out_file', calc, 'in_file_a')
    preproc.connect(tstat, 'out_file', calc, 'in_file_b')
    #5.
    preproc.connect(calc, 'out_file', fourier, 'in_file')
    #6.
    preproc.connect(fourier, 'out_file', merge, 'in_file')
    #7.
    preproc.connect(grabber, 't1_file', converter, 'in_file')
    #8.
    preproc.connect(converter, 'out_file', register1, 'base_file')
    preproc.connect(calc, 'out_file', register1, 'in_file')
    preproc.connect(converter, 'out_file', register1, 'master_file')
    preproc.connect(register1, 'onedmatrix', register2, 'onedmatrix')
    preproc.connect(fourier, 'out_file', register2, 'in_file')
    #9.
    preproc.connect(grabber, 'csf_file', csfmask, 'input_file')
    preproc.connect(grabber, 'wm_file', wmmask, 'input_file')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Preprocessing script for resting state fMRI')
    parser.add_argument('-t', '--filetype', action='store', dest='outputType', type=str, required=False,
                        default='NIFTI', help='File type for outputs.  Values: "NIFTI_GZ", "AFNI", "NIFTI" (default)'
    parser.add_argument('-n','--name', action='store', dest='name', type=str, required=True,
                        help='Name (required)')
    parser.add_argument('-s', '--session', action='store', dest='sessionID', type=list, required=True,
                        help='list of session IDs (required)')
    args = parser.parse_args()
    freesurferOutputTypes = {"NIFTI_GZ" : "niigz",
                             "AFNI" : "afni",
                             "NIFTI" : "nii"}
    args[fOutputType] = freesurferOutputTypes[args['outputType']]

    sys.exit(pipeline(args))
