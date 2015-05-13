import nipype.pipeline.engine as pipe
from nipype.interfaces.utility import Merge
from nipype.interfaces.ants.registration import Registration
from nipype.interfaces.ants import ApplyTransforms

import SEMTools as sem

import dataio
import afninodes
from nodes import applyTransformNode


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


def mergeNode(count, name):
    pass


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

    # Nodes
    bFit = transformNode(name='brainsFit')
    tstat = afninodes.tstatnode(outputType, 'tstat')  # fMRI reference image
    transformGrabber = dataio.transformGrabber(experiment=t1_experiment)

    # NAC -> T1 -> fMRI
    nactoFMRI_list = pipe.Node(Merge(2), name='nactoFMRI_list')
    warpSeedtoFMRI  = applyTransformNode(name="warpSeedtoFMRI", transform='nac2fmri', interpolation='NearestNeighbor')
    # T1 -> fMRI
    t1toFMRI_list = pipe.Node(Merge(1), name='t1toFMRI_list')
    warpT1toFMRI   = applyTransformNode(name='warpT1toFMRI', transform='t12fmri')  # Validate using 'NearestNeighbor'
    warpBABCtoFMRI = applyTransformNode(name='warpBABCtoFMRI', transform='t12fmri', interpolation='NearestNeighbor')
    # fMRI -> T1 -> NAC
    fmritoNAC_list = pipe.Node(Merge(2), name='fmritoNAC_list')
    warpEPItoNAC = applyTransformNode(name='warpEPItoNAC', transform='fmri2nac')
    warpLabeltoNAC  = applyTransformNode(name='warpLabeltoNAC', transform='fmri2nac')

    register.connect([(tstat, bFit,                     [('out_file', 'movingVolume')]),
                      # (dataGrabber, bFit              [('t1_file', 'fixedVolume')]),
                      # NAC -> T1 -> fMRI
                      ###################
                      (bFit,             nactoFMRI_list, [('outputTransform', 'in1')]),
                      (transformGrabber, nactoFMRI_list, [('atlasToSessionTransform', 'in2')]),
                      # seeds
                      (warpT1toFMRI,   warpSeedtoFMRI,   [('output_image', 'reference_image')]),
                      (nactoFMRI_list, warpSeedtoFMRI,   [('out', 'transforms')]),
                      # T1 => fMRI
                      ############
                      (bFit, t1toFMRI_list,              [('outputTransform', 'in1')]),
                      # T1
                      (tstat,         warpT1toFMRI,      [('out_file', 'reference_image')]),
                      (t1toFMRI_list, warpT1toFMRI,      [('out', 'transforms')]),
                      # BABC labels
                      (tstat,         warpBABCtoFMRI,    [('out_file', 'reference_image')]),
                      (t1toFMRI_list, warpBABCtoFMRI,    [('out', 'transforms')]),
                      # fMRI => NAC
                      #############
                      (transformGrabber, fmritoNAC_list, [('sessionToAtlasTransform', 'in1')]),
                      (bFit,             fmritoNAC_list, [('outputTransform', 'in2')]),
                      # session-specific seeds
                      (fmritoNAC_list, warpEPItoNAC,     [('out', 'transforms')]),
                      (fmritoNAC_list, warpLabeltoNAC,   [('out', 'transforms')]),
                      ])

    register.write_graph(dotfilename='register.dot', graph2use='orig', format='png', simple_form=False)
    return register
