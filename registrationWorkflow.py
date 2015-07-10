import nipype.pipeline.engine as pipe
from nipype.interfaces.utility import Merge, IdentityInterface
from nipype.interfaces.ants import ApplyTransforms
# from nipype.interfaces.ants.registration import Registration

import SEMTools as sem

import dataio
import afninodes


def transformNode(name):
    """ BRAINSFit node """
    bFit = pipe.Node(interface=sem.BRAINSFit(), name=name)
    bFit.inputs.initializeTransformMode = 'useCenterOfHeadAlign'
    bFit.inputs.maskProcessingMode = 'ROIAUTO'
    bFit.inputs.ROIAutoDilateSize = 10.0
    bFit.inputs.useRigid = True
    bFit.inputs.costMetric = 'MMI'  # (default)
    bFit.inputs.numberOfSamples = 100000  # (default)
    bFit.inputs.outputTransform = True
    return bFit


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


def t1Workflow():
    wkfl = pipe.Workflow(name='t1_wf')
    warp = applyTransformNode(name='warpT1toFMRI', transform='t12fmri')  # Validate using 'NearestNeighbor'
    wkfl.add_nodes([warp])
    return wkfl


def babcWorkflow():
    wkfl = pipe.Workflow(name='babc_wf')
    warp = applyTransformNode(name='warpBABCtoFMRI', transform='t12fmri', interpolation='NearestNeighbor')
    wkfl.add_nodes([warp])
    return wkfl


def epiWorkflow():
    wkfl = pipe.Workflow(name='epi_wf')
    warp = applyTransformNode(name='warpEPItoNAC', transform='fmri2nac')
    wkfl.add_nodes([warp])
    return wkfl


def labelWorkflow():
    wkfl = pipe.Workflow(name='lb_wf')
    warp = applyTransformNode(name='warpLabeltoNAC', transform='fmri2nac')
    wkfl.add_nodes([warp])
    return wkfl


def seedWorkflow():
    wkfl = pipe.Workflow(name='seed_wf')
    warp = applyTransformNode(name="warpSeedtoFMRI", transform='nac2fmri', interpolation='NearestNeighbor')
    wkfl.add_nodes([warp])
    return wkfl


def workflow(t1_experiment, outputType, name):
    """ registration workflow

    Connections needed:
      (in)
      tstat.in_file,
      tstat.mask_file, **optional**
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
    register = pipe.Workflow(name=name)

    # Nodes
    in_fields = ['fmri', 't1', 'session_id']
    inputnode = pipe.Node(interface=IdentityInterface(fields=in_fields), name='inputs')
    bFit = transformNode(name='brainsFit')
    tstat = afninodes.tstatnode(outputType, 'tstat')  # fMRI reference image
    txGrabber = dataio.transformGrabber(experiment=t1_experiment)
    # NAC -> T1 -> fMRI
    nactoFMRI_list = pipe.Node(Merge(2), name='nactoFMRI_list')
    register.connect([(bFit,      nactoFMRI_list, [('outputTransform', 'in1')]),
                      (txGrabber, nactoFMRI_list, [('atlasToSessionTransform', 'in2')])])
    # T1 -> fMRI
    t1toFMRI_list = pipe.Node(Merge(1), name='t1toFMRI_list')
    register.connect([(bFit,      t1toFMRI_list,  [('outputTransform', 'in1')])])
    # fMRI -> T1 -> NAC
    fmritoNAC_list = pipe.Node(Merge(2), name='fmritoNAC_list')
    register.connect([(txGrabber, fmritoNAC_list, [('sessionToAtlasTransform', 'in1')]),
                      (bFit,      fmritoNAC_list, [('outputTransform', 'in2')])])

    out_fields = ['fmri_reference', 't1_reference', 'nac2fmri_list', 't12fmri_list', 'fmri2nac_list']
    outputnode = pipe.Node(interface=IdentityInterface(fields=out_fields), name='outputs')

    register.connect([(inputnode, bFit,           [('t1', 'fixedVolume')]),
                      (inputnode, tstat,          [('fmri', 'in_file')]),
                      (inputnode, txGrabber,      [('session_id', 'session_id')]),
                      (tstat,     bFit,           [('out_file', 'movingVolume')]),
                      (inputnode, outputnode,     [('t1', 't1_reference')]),
                      (tstat, outputnode,         [('out_file', 'fmri_reference')]),
                      (nactoFMRI_list, outputnode,[('out', 'nac2fmri_list')]),
                      (t1toFMRI_list, outputnode, [('out', 't12fmri_list')]),
                      (fmritoNAC_list, outputnode,[('out', 'fmri2nac_list')]),
                      ])
    return register
