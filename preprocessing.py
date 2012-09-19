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
import dicom

def dicomToNrrd(fileName):
    import shlex
    import subprocess
    _cmd = '/ipldev/sharedopt/20120722/Darwin_i386/DicomToNrrdConverter/DWIConvert-build/DWIConvert'
    args = '--conversionMode DicomToNrrd --outputDirectory %s --inputDicomDirectory %s' % () # TODO: pass relevant args
    command = shlex.split(" ".join(_cmd, args))
    subprocess.call(command, shell=True)
    return outputFiles # TODO: get the path to the outfile

def readNrrdHeader(fileName):
    headerName = dicomToNrrd(fileName)
    # TODO: grep header for numberOfSlices and numberOfVolumes
    return numberOfSlices, numberOfVolumes

def generateTissueMask(inputLabelMap, low=0.0, high=1.0):
    import SimpleITK as sitk

    #GeodesicActiveContourLevelSet
    def getLargestConnectedRegion(input_file, low=0, high=100):
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
    csfMask = sitk.ErodeObjectMorphology(binaryMask, radiusMM, kernel, 1, 0)
    ## erodedMask = sitk.ErodeObjectMorphology(binaryMask, 2)
    ## csfMask = sitk.BinaryThreshold(input_file * erodedMask, ### low, ### high, 1, 0)
    csfMaskOnly = getLargestConnectedRegion(csfMask)
    csfEroded = sitk.ErodeObjectMorphology(csfMaskOnly, 1)
    # Binary dilation
    csfDilated = sitk.DilateObjectMorphology(csfEroded, 5)
    csfMaskFinal = getLargestConnectedRegion(csfDilated)
    return csfMaskFinal * csfMask


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
                                              outfields=['fmri_folder', 't1_file',
                                                         'csf_file', 'wm_file']),
        name='dataGrabber')
    grabber.inputs.base_directory = '/paulsen'
    grabber.inputs.template = '*'
    ### TODO: rename 'outfolder' to 'fmri_folder'
    ### TODO: rename 'out_file' from t1grabber to 't1_file'
    grabber.inputs.field_template = dict(fmri_folder='MRx/FMRI_HD_120/*/%s/%s/%s/*'
                                         t1_file='Experiments/20120722_JOY_DWI/FMRI_HD_120/*/%s/%s/%s_*_%s/T1.mgz'
                                         csf_file='Experiments/20120801.SubjectOrganized_Results/PHD_120/*/%s/%s/POSTERIOR_CSF_TOTAL.nii.gz',
                                         wm_file= 'Experiments/20120801.SubjectOrganized_Results/PHD_120/*/%s/%s/POSTERIOR_WM_TOTAL.nii.gz')
    grabber.inputs.template_args = dict(fmri_folder=[['session_id', 'ANONRAW', 'FMRI_RestingStateConnectivity']],
                                         t1_file=[['session_id', '10_AUTO.NN3Tv20110419','JOY_v51_2011', 'session_id']],
                                         csf_file=[['session_id','ACCUMULATED_POSTERIORS']],
                                         wm_file=[['session_id','ACCUMULATED_POSTERIORS']])

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

    dicomNrrd = pipe.Node(interface=Function(), name='dicomToNrrd') # TODO: Flush this out

    to_3D = pipe.Node(interface=To3D(), name='afniTo3D')
    to_3D.inputs.outputtype = outputType
    to_3D.inputs.infolder = # TODO: dataGrabber node
    to_3D.inputs.funcparams = "-time:zt ${numberOfSlices} ${numberOfVolumes} ${repetitionTime} alt-z" # TODO: connect to a string creation node
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
    # volreg.outputs.onedfile = 'Rest_mt.1D' # TODO: Verify with Jatin if should connect to a dataGrabber?
    volreg.inputs.args = '-cubic -maxite 50 -x_thresh 0.001 -rot_thresh 0.001 -delta 0.1 -final Fourier \
                          -twopass -twodup -coarse 2 2 -coarserot -base 9 '

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

    convert = pipe.Node(interface=MRIConvert(), name='freesurferMRIconvert')
    convert.inputs.in_orientation = 'RAS'
    convert.inputs.in_type = 'mgz'
    convert.inputs.out_orientation = 'RAS'
    convert.inputs.out_type = fOutputType

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
    preproc.connect(grabber, 'fmir_folder', to_3D, 'infolder')
    preproc.connect(to_3D, 'out_file', refit, 'in_file')
    #2.
    preproc.connect(to_3D, 'out_file', despike, 'in_file') # TODO: connect refit to despike???
    #3.
    preproc.connect(despike, 'out_file', volreg, 'in_file')
    #4.
    preproc.connect(volreg, 'out_file', tstat, 'in_file')
    preproc.connect(volreg, 'out_file', calc, 'in_file_a')
    preproc.connect(tstat, 'out_file', calc, 'in_file_b')
    #5.
    preproc.connect(calc, 'out_file', fourier, 'in_file')
    #6.
    preproc.connect(fourier, '_file', merge, 'in_file')
    #7.
    preproc.connect(grabber, 't1_file', convert, 'in_file')
    #8.
    preproc.connect(convert, 'out_file', register1, 'base_file')
    preproc.connect(calc, 'out_file', register1, 'in_file')
    preproc.connect(convert, 'out_file', register1, 'master_file')
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
