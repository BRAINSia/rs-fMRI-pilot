import nipype.pipeline.engine as pipe
from afni.preprocess import *
from nipype.interfaces.utility import Function


def strCreate(time, slices, volumes, repTime, order):
    """ Create a funcparams string for To3D() """
    return "%s %s %s %s %s" % (time, slices, volumes, repTime, order.strip())


def to3dstrnode(name):
    node = pipe.Node(interface=Function(function=strCreate,
                                        input_names=['time', 'slices', 'volumes', 'repTime', 'order'],
                                        output_names=['funcparams']),
                          name=name)
    node.inputs.time = 'zt'
    return node


def to3dnode(name):
    """ Use 3dAFNIto3D() to convert AFNI data to 3D file """
    node = pipe.Node(interface=To3D(), name=name)
    node.inputs.outputtype = 'AFNI'
    node.inputs.datum = 'short'
    node.inputs.filetype = 'epan'
    node.inputs.prefix = 'to_3D_out'
    # node.inputs.funcparams = strCreate()
    return node


def refitnode(name):
    """ Use 3dRefit() to replace transformation matrix in header with cardinal matrix """
    node = pipe.Node(interface=Refit(), name=name)
    node.inputs.outputtype = 'AFNI'
    node.inputs.deoblique = True
    return node


def despikenode(outputType, skipCount, name):
    """ Use 3dDespike() to remove data 'spikes' """
    node = pipe.Node(interface=Despike(), name=name)
    node.inputs.outputtype = outputType
    node.inputs.start = skipCount
    node.inputs.ignore = skipCount
    return node


def volregnode(outputType, name):
    """ Use 3dvolreg() to register sub-bricks to base brick """
    node = pipe.Node(interface=Volreg(), name=name)
    node.inputs.outputtype = outputType
    node.inputs.timeshift = False  # 0
    node.inputs.zpad = 3
    node.inputs.interp = 'cubic'
    node.inputs.maxite = 50
    node.inputs.thresh = 0.001
    node.inputs.rot_thresh = 0.001
    node.inputs.delta = 0.1
    node.inputs.final = 'Fourier'
    node.inputs.twopass = True
    node.inputs.twodup = True
    node.inputs.coarse = [2, 2]
    node.inputs.coarserot = True
    node.inputs.base = 9
    node.inputs.oned_file = 'volReg.1D'
    return node


def zeropadnode(name):
    """ Use 3dZeropad() to pad the image along the inferior-superior axis s.t. the final dataset has 44 planes """
    node = pipe.Node(interface=Zeropad(), name=name)
    node.inputs.plane = 'IS'
    node.inputs.numberOfPlanes = 44
    node.inputs.is_mm = False
    return node


def mergenode(outputType, name):
    """
    Use 3dmerge() to modify all sub-bricks in the dataset:
      * set all negative values -> 0,
      * set all values over 100 -> 0, and
      * apply a gaussian blur of 6mm on the data
    """
    node = pipe.Node(interface=Merge(), name=name)
    node.inputs.outputtype = outputType
    node.inputs.blurfwhm = 6
    node.inputs.doall = True
    # TODO: implement
    # node.inputs.onenoneg = True
    # node.inputs.oneclip = 100
    node.inputs.args = '-1noneg -1clip 100'
    # END TODO
    return node


def automasknode(outputType, name):
    """ Use 3dAutomask() to create a brain-only mask with 1mm dilation """
    node = pipe.Node(interface=Automask(), name=name)
    node.inputs.outputtype = outputType
    node.inputs.dilate = 1
    return node


def tstatnode(outputType, name):
    """ Use 3dTstat() to compute voxel-wise mean """
    node = pipe.Node(interface=TStat(), name=name)
    node.inputs.outputtype = outputType
    node.inputs.args = '-mean'  # TODO
    # node.inputs.suffix = '_tstat'
    return node


def multiplynode(outputType, name):
    """ Use 3dCalc() to multiply two images together """
    node = pipe.Node(interface=Calc(letters=['a', 'b']), name=name)
    node.inputs.outputtype = outputType
    node.inputs.expr = "a * b"
    return node


def fouriernode(outputType, name):
    """
    Use 3dFourier() to apply a bandpass filter on the data (high=0.011, low=0.1).
    This should be done ***after*** detrending!
    """
    node = pipe.Node(interface=Fourier(), name=name)
    node.inputs.outputtype = outputType
    node.inputs.highpass = 0.011
    node.inputs.lowpass = 0.1
    # node.inputs.args = '-retrend' # removed on 9/20/13
    return node


def maskavenode(outputType, name, args=''):
    """ Use 3dmaskave() to compute the median voxel value for the 3D data set """
    node = pipe.Node(interface=Maskave(), name=name)
    node.inputs.outputtype = outputType
    node.inputs.args = '-median ' + args  # TODO
    node.inputs.quiet = True
    return node


def deconvolvenode(labels, name):
    """
    Use 3dDeconvolve() to compute the linear regression of the data with varying stimulus inputs determined by the
    labels provided at runtime.  Output needed for 3dDetrend()
     """
    defaults = ("roll", "pitch", "yaw", "dS", "dL", "dP")
    labels = tuple(labels)
    all_labels = labels + defaults
    node = pipe.Node(interface=Deconvolve(fileCount=(len(labels) + 1), seriesCount=len(all_labels)), name=name)
    # node.inputs.outputtype = outputType # 'AFNI'
    node.inputs.ignoreWarnings = 10
    node.inputs.nullHypothesisPolynomialDegree = 1
    node.inputs.full_first = True
    node.inputs.is_float = True
    node.inputs.tout = True
    node.inputs.rout = True
    node.inputs.fout = True
    node.inputs.bucket = 'Rest_bp_Decon'
    node.inputs.fitts = 'full_fitts_Decon'
    node.inputs.errts = 'errts_Decon'
    # tuple is zero-based indexed, while Deconvolve uses one-based indexing, so we use enumerate with start=1
    for index, value in enumerate(labels, start=1):
        stim_key = "stim_label_{0}".format(index)
        setattr(node.inputs, stim_key, value)

        base_key = "is_stim_base_{0}".format(index)
        if index <= (len(labels) + 1):
            setattr(node.inputs, base_key, False)
        else:
            setattr(node.inputs, base_key, True)
    return node


def detrendnode(outputType, name):
    """
    Use 3dDetrend() to detrend the data.
    This should always be done ***before*** applying a bandpass filter
    """
    node = pipe.Node(interface=Detrend(), name=name)
    node.inputs.outputtype = outputType
    node.inputs.suffix = '_dt'
    node.inputs.args = '-polort 3'  # TODO
    return node


def fimnode(out, name):
    """
    Use 3dfim+() to calculate the cross-correlation of an ideal reference waveform with the measured FMRI
    time series for each voxel
    """
    node = pipe.Node(interface=Fim(), name=name)
    node.inputs.out = out
    return node


def logcalcnode(outputType, name):
    """
    Use 3dCalc() to calculate math::

        \frac{log(\frac{1 + a}{1 - a})}{2}

    per voxel.
    """
    node = pipe.Node(interface=Calc(letters=['a']), name=name)
    node.inputs.outputtype = outputType
    node.inputs.expr = 'log((1+a)/(1-a))/2'
    return node
