import nipype.pipeline.engine as pipe
from nipype.interfaces.utility import Merge as Merger

import dataio
import afninodes


def transformNode(name):
    """ BRAINSFit node """
    bFit = pipe.Node(interface=sem.BRAINSFit(), name='brainsFit')
    bFit.inputs.initializeTransformMode = 'useCenterOfHeadAlign'
    bFit.inputs.maskProcessingMode = 'ROIAUTO'
    bFit.inputs.ROIAutoDilateSize = 10.0
    bFit.inputs.useRigid = True
    bFit.inputs.costMetric = 'MMI'  # (default)
    bFit.inputs.numberOfSamples = 100000  # (default)
    bFit.inputs.outputTransform = True
    return bFit


def mergeNode(count, name):
    pass


def applyTransformsNode(name, **kwargs):
    """ input fields are kwargs e.g. 'interpolation', invert_transform_flags', etc. """
    node = pipe.Node(interface=ApplyTransforms(), iterfield=['input_image'], name=name)
    for k, v in kwargs:
        setattr(node.inputs, k, v)
    return node


def workflow(t1_experiment, outputType, name):
    """ registration workflow

    Connections needed:
      (in)
      tstat.in_file,
      tstat.mask_file,
      transformGrabber.session_id,
      bFit.fixedVolume

      (out)
      tstat.out_file
      forwardTransformT1ToFMRI.out

    ##########################################################################
    # The registration of the T1 -> fMRI is good _most_ of the time, but     #
    # sometimes it hits a local minimum during optimization. To prevent this,#
    # we set the fMRI image as the moving image and the T1 as the fixed.     #
    # When we run ApplyTransforms, we pass it the inverse rigid transform    #
    # file (or give it the 'inverse' flag)                                   #
    ##########################################################################
    """

    register = pipe.Workflow(updatehash=True, name=name)

    transformGrabber = dataio.transformGrabber(experiment=t1_experiment)
    # Concatenate transforms
    # NOTE: antsApplyTransforms takes transforms in REVERSE order!!!
    forwardTransform = pipe.Node(Merger(2), name='0_List_forwardTransformNACToFMRI')  # NAC -> T1 -> fMRI
    forwardTransformT1ToFMRI = pipe.Node(Merger(1), name='0_List_forwardTransformT1ToFMRI')  # T1 -> fMRI
    reverseTransform = pipe.Node(Merger(2), name='0_List_reverseTransformFMRIToNAC')  # NAC <- T1 <- fMRI

    tstat = afninodes.tstatnode(outputType, 'afni3DtStat')
    fourier = afninodes.fouriernode(outputType, 'afni3Dfourier') # Fourier is the last NIFTI file format in the AFNI pipeline

    bFit = transformNode(name='brainsFit')

    warpT1ToFMRI = applyTransformNode(name='antsApplyTransformsT1',
                                      interpolation='Linear',   # Validation test value: NearestNeighbor
                                      invert_transform_flags=[True])

    warpWHM = warpT1ToFMRI.clone('antsApplyTransformWHM')
    warpGRM = warpT1ToFMRI.clone('antsApplyTransformGRM')  # maskGM

    warpBABCSegToFMRI = applyTransformNode(name='warpBABCSegToFMRI',
                                           interpolation='NearestNeighbor',
                                           invert_transform_flags=[True])

    warpAtlasCSF = applyTransformNode(name='antsApplyTransformsCSF',
                                      interpolation='Linear',
                                      invert_transform_flags = [True, False])

    warpWholeMask = applyTransformNode(name='antsApplyTransformsWholeBrain',
                                       interpolation='Linear',
                                       invert_transform_flags = [True, False])  # maskWholeBrain

    nacToFMRI = applyTransformNode(name="warpNACseedsToFMRI",   # Warp seed output to FMRI
                                   interpolation='NearestNeighbor',
                                   invert_transform_flags=[True, False])

    fmriToNAC_epi = applyTransformNode(name='fmriToNac_epi',
                                       interpolation = 'Linear',
                                       invert_transform_flags = [False, False])

    register.connect([(tstat, bFit,                                 [('out_file', 'movingVolume')]),
                      # pass images to the warpers
                      (tstat, warpT1ToFMRI,                         [('out_file', 'reference_image')]),
                      (tstat, warpBABCSegToFMRI,                    [('out_file', 'reference_image')]),
                      (fourier, fmriToNAC_epi,                      [('out_file', 'input_image')]),
                      # merge the forward transform
                      (bFit, forwardTransform,                      [('outputTransform', 'in1')]),
                      (transformGrabber, forwardTransform,          [('atlasToSessionTransform', 'in2')]),
                      # merge the partial forward transform
                      (bFit, forwardTransformT1ToFMRI,              [('outputTransform', 'in1')]),
                      # merge the reverse transform
                      (transformGrabber, reverseTransform,          [('sessionToAtlasTransform', 'in1')]),
                      (bFit, reverseTransform,                      [('outputTransform', 'in2')]),
                      # pass the transforms to the warpers
                      (forwardTransformT1ToFMRI, warpT1ToFMRI,      [('out', 'transforms')]),
                      (forwardTransformT1ToFMRI, warpBABCSegToFMRI, [('out', 'transforms')]),
                      (forwardTransform, nacToFMRI,                 [('out', 'transforms')]),
                      (reverseTransform, fmriToNAC_epi,             [('out', 'transforms')]),
                      ])
    return register
