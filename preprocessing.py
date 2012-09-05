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

def pipeline(**kwds):
    preproc = pipe.Workflow()

    # filePaths = glob.glob('/paulsen/MRx/FMRI_HD_*/*/%s/ANONRAW/FMRI_RestingStateConnectivity/*/*.IMA' % sessionID)
    # if len(filePath) == 0:
    #     filePaths = glob.glob('/paulsen/MRx/FMRI_HD_*/*/%s/ANONRAW/FMRI_RestingStateConnectivity/*/*.dcm' % sessionID)
    # header = dicom.read_file(filePaths[0])
    # repetitionTime = header.RepetitionTime

    sessions = pipe.Node(interface=IdentityInterface(fields=['session_id']), name='sessionIDs')
    sessions.iterables = ('session_id', sessionID)

    fMRIgrabber = pipe.Node(interface=DataGrabber(infields=['session_id'], outfields=['outfolder']),
                            name='fMRIdata')
    fMRIgrabber.inputs.base_directory = '/paulsen/MRx/FMRI_HD_120'
    fMRIgrabber.inputs.template = '*/%s/ANONRAW/FMRI_RestingStateConnectivity/*'
    fMRIgrabber.inputs.template_args['outfolder'] = [['session_id', 'session_id']]

    t1grabber = pipe.Node(interface=DataGrabber(infields=['session_id'], outfields=['out_file']),
                          name='T1data')
    t1grabber.inputs.base_directory = '/paulsen/Experiments/20120722_JOY_DWI/FMRI_HD_120'
    t1grabber.inputs.template = '*/%s/10_AUTO.NN3Tv20110419/JOY_v51_2011_*_%s/T1.mgz'
    t1grabber.inputs.template_args['out_file'] = [['session_id', 'session_id']]

    if len(sessionID) > 1:
        preprocess.connect(sessions, 'session_id', fMRIgrabber, 'session_id')
        preprocess.connect(sessions, 'session_id', t1grabber, 'session_id')
    else:
        fMRIgrabber.inputs.session_id = sessionID[0]
        t1grabber.inputs.session_id = sessionID[0]

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
    despike.inputs.args = '-prefix Rest_ds -ignore 4' # -ignore 4

    volreg = pipe.Node(interface=Volreg(), name='afni3DvolReg')
    volreg.inputs.outputtype = outputType
    volreg.inputs.timeshift = False # 0
    volreg.inputs.zpad = 3
    # volreg.inputs.onedfile = 'Rest_mt.1D' # TODO: Verify with Jatin if should connect to a dataGrabber?
    volreg.inputs.args = '-cubic -maxite 50 -x_thresh 0.001 -rot_thresh 0.001 -delta 0.1 -final Fourier \
                          -twopass -twodup -coarse 2 2 -coarserot -base 9 '

    tstat = pipe.Node(interface=TStat(), name='afni3DtStat')
    tstat.inputs.outputtype = outputType
    tstat.inputs.args = '-mean'

    calc = pipe.Node(interface=Calc(), name='afni3Dcalc')
    calc.inputs.outputtype = outputType
    calc.inputs.expr = 'a-b'

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

    ### TODO: check/implement this node ###
    zeroPad = pipe.Node(interface=Zeropad(), name='afni3DzeroPad')
    zeroPad.inputs.direction_list = ['IS']
    zeroPad.inputs.pad_count_list = [44]
    zeroPad.inputs.mm = False

    convert = pipe.Node(interface=MRIConvert(), name='freesurferMRIconvert')
    convert.inputs.in_orientation = 'RAS'
    convert.inputs.in_type = 'mgz'
    convert.inputs.out_orientation = 'RAS'
    convert.inputs.out_type = fOutputType

    ### TODO ###
    register = pipe.Node(interface=????(), name='??????')

    # Connect pipeline
    #1.
    preproc.connect(fMRIgrabber, 'outfolder', to_3D, 'infolder')
    preproc.connect(to_3D, 'out_file', refit, 'in_file')
    #2.
    preproc.connect(to_3D, 'out_file', despike, 'in_file') # TODO: connect refit to despike???
    #3.
    preproc.connect(despike, 'out_file', volreg, 'in_file')
    preproc.connect('?????', 'out_file', volreg, 'onedfile') # TODO: Verify with Jatin if should connect to a dataGrabber?
    #4.
    preproc.connect(volreg, 'out_file', tstat, 'in_file')
    preproc.connect(volreg, 'out_file', calc, 'in_file_a')
    preproc.connect(tstat, 'out_file', calc, 'in_file_b')
    #5.
    preproc.connect(volreg, 'out_file', fourier, 'in_file') # TODO: Verify with Jatin - where does 'Rest_zpad' come from?
    #6.
    preproc.connect(fourier, 'out_file', merge, 'in_file')
    #7.
    preproc.connect(calc, 'out_file', zeroPad, 'in_file')
    preproc.connect(zeroPad, 'out_file', ???, '???') # TODO: Verify with Jatin - 'Rest_zpad' created here?
    #8.
    preproc.connect(t1grabber, 'out_file', convert, 'in_file')
    ### TODO: Decide on registration method and direction ###
    preproc.connect(convert, 'out_file', register, 'fixed') # Do we want to use antsRegistration or BRAINSFit?
    preproc.connect(zeroPad, 'out_file', register, 'moving') # Should fixed be the T1 or the fMRI data?  Does it matter?

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
