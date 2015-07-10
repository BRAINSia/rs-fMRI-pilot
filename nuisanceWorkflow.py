import nipype.pipeline.engine as pipe
from nipype.interfaces.utility import Function  #, IdentityInterface

import afninodes
from utilities import generateTissueMask
from registrationWorkflow import applyTransformNode


def maskNode(name, **inputkwargs):
    binaries = ['erode', 'binary', 'largest']
    inputs = ['input_file', 'fileName', 'low', 'high'] + binaries
    outputs = ['output_file']
    node = pipe.Node(interface=Function(function=generateTissueMask, input_names=inputs, output_names=outputs),
                     name=name)
    flags = inputkwargs.pop('flags', [])
    for k in binaries:
        if k in flags:
            setattr(node.inputs, k, True)
        else:
            setattr(node.inputs, k, False)

    for k, v in inputkwargs.items():
        setattr(node.inputs, k, v)
    return node


def csfWorkflow(input_file):
    # Nodes
    warp = applyTransformNode(name='warpCSFtoFMRI', transform='nac2fmri')
    mask = maskNode(name='csfMask', fileName='csfMask.nii', low=3, high=42, flags=['binary'], input_file=input_file)
    avg = afninodes.maskavenode('AFNI_1D', 'afni3DmaskAve_csf')
    # Pipeline
    csf = pipe.Workflow(name='csf')
    csf.connect([(warp, mask, [('output_image', 'input_file')]),
                 (mask, avg,  [('output_file', 'mask')]),
                ])
    return csf


def wmWorkflow():
    # Nodes
    warp = applyTransformNode(name='warpWMtoFMRI', transform='t12fmri')
    mask = maskNode(name='wmMask', fileName='whiteMatterMask.nii', low=0.99, high=1.0, flags=['erode'])
    avg = afninodes.maskavenode('AFNI_1D', 'afni3DmaskAve_wm')
    # Pipeline
    wm = pipe.Workflow(name='wm')
    wm.connect([(warp, mask, [('output_image', 'input_file')]),
                (mask, avg,  [('output_file', 'mask')]),
               ])
    return wm


def gmWorkflow():
    # Nodes
    warp = applyTransformNode(name='warpGMtoFMRI', transform='t12fmri')
    mask = maskNode(name='grmMask', fileName='grmMask.nii', low=0.99, high=1.0)
    avg = afninodes.maskavenode('AFNI_1D', 'afni3DmaskAve_grm')
    # Pipeline
    gm = pipe.Workflow(name='gm')
    gm.connect([(warp, mask, [('output_image', 'input_file')]),
                (mask, avg,  [('output_file', 'mask')]),
               ])
    return gm


def wbWorkflow(input_mask):
    # Nodes
    warp = applyTransformNode(name='warpBraintoFMRI', transform='nac2fmri')
    warp.inputs.input_image = input_mask
    mask = maskNode(name='wholeBrainMask', fileName='wholeBrainMask.nii', low=0.5, high=1.0, flags=['largest'])
    avg = afninodes.maskavenode('AFNI_1D', name='afni3DmaskAve_whole')
    # Pipeline
    wb = pipe.Workflow(name='wb')
    wb.connect([(warp, mask, [('output_image', 'input_file')]),
                (mask, avg,  [('output_file', 'mask')]),
               ])
    return wb


def workflow(outputType, name="nuisance", **kwargs):
    # Nodes
    # in_fields = ["oned_file", "nactoFMRI", "t1toFMRI", "avg_file"]
    # inputnode = pipe.Node(interface=IdentityInterface(fields=in_fields), name="inputs")
    csf_sub = csfWorkflow(kwargs['csf_input'])
    wm_sub = wmWorkflow()

    nuisance = pipe.Workflow(name=name)
    # nuisance.connect([(inputnode, csf_sub, [('avg_file', 'afni3DmaskAve_csf.in_file')

    if kwargs['maskgm']:
        gm_sub = gmWorkflow()  # Mask gray matter
        deconvolve = afninodes.deconvolvenode(("Median_CSF", "Median_WM", "Median_GM"), "afni3Ddeconvolve")
        nuisance.connect(gm_sub, 'afni3DmaskAve_grm.out_file', deconvolve, 'stim_file_3')
    elif kwargs['maskwb']:
        wb_sub = wbWorkflow(kwargs['wb_input'])  # Mask the whole brain
        deconvolve = afninodes.deconvolvenode(("Median_CSF", "Median_WM", "Median_WholeBrain"), "afni3Ddeconvolve")
        nuisance.connect(wb_sub, 'afni3DmaskAve_whole.out_file', deconvolve, 'stim_file_3')
    else:
        deconvolve = afninodes.deconvolvenode(("Median_CSF", "Median_WM"), "afni3Ddeconvolve")
    nuisance.connect(csf_sub, 'afni3DmaskAve_csf.out_file', deconvolve, 'stim_file_1')
    nuisance.connect(wm_sub, 'afni3DmaskAve_wm.out_file', deconvolve, 'stim_file_2')
    return nuisance
