# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
import nipype.pipeline.engine as pipe
from nipype.interfaces.io import DataSink, DataGrabber
from nipype.interfaces.utility import Function, IdentityInterface
from nipype.interfaces.afni.preprocess import *
from nipype.interfaces.freesurfer.preprocess import *
from nipype.interfaces.ants import ApplyTransforms


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


def applyTransformNode(name, transform, **kwargs):
    """ input fields are kwargs e.g. 'interpolation', 'invert_transform_flags', etc. """
    kwargs.setdefault('interpolation', 'Linear')
    if transform == 'fmri2nac':
        kwargs['invert_transform_flags'] = [False, False]
    elif transform == 't12fmri':
        kwargs['invert_transform_flags'] = [True]
    elif transform == 'nac2fmri':
        kwargs['invert_transform_flags'] = [True, False]
    else:
        pass
    # NOTE: antsApplyTransforms takes transforms in REVERSE order!!!
    node = pipe.Node(interface=ApplyTransforms(), iterfield=['input_image'], name=name)
    for k, v in kwargs.items():
        setattr(node.inputs, k, v)
    return node


def sessionNode(outputType):
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
    grabber.inputs.base_directory = '/Shared/paulsen'
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
