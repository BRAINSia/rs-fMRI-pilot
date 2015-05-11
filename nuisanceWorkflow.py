import nipype.pipeline.engine as pipe
from nipype.interfaces.utility import Merge as Merger

import dataio
import afninodes
from utilities import generateTissueMask


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


def workflow(nacAtlasLabel, outputType, name):

    preproc = pipe.Workflow(updatehash=True, name=name)

    csfmask = maskNode(name='csfMask', fileName='csfMask.nii', low=3, high=42, flags=['binary'], input_file=nacAtlasLabel)
    csfAvg = afninodes.maskavenode('AFNI_1D', 'afni3DmaskAve_csf')

    wmmask = maskNode(name='wmMask', fileName='whiteMatterMask.nii', low=0.99, high=1.0, flags=['erode'])
    wmAvg = afninodes.maskavenode('AFNI_1D', 'afni3DmaskAve_wm')
    preproc.connect(wmmask, 'out_file', wmAvg, 'mask')

    if maskGM:
        #------------------------------ GRAY MATTER MASK ------------------------------
        grmmask = maskNode(name='grmMask', fileName='grmMask.nii', low=0.99, high=1.0)
                           largest=False)

        grmAvg = maskavenode('AFNI_1D', 'afni3DmaskAve_grm')
        preproc.connect(grmmask, 'output_file', grmAvg, 'mask')
        # 12
        deconvolve = afninodes.deconvolvenode(("Median_CSF", "Median_WM", "Median_GM"), "afni3Ddeconvolve")
        preproc.connect(csfAvg, 'out_file', deconvolve, 'stim_file_1')
        preproc.connect(wmAvg, 'out_file', deconvolve, 'stim_file_2')
        preproc.connect(grmAvg, 'out_file', deconvolve, 'stim_file_3')
        preproc.connect(volreg, 'oned_file', deconvolve, 'stim_file_4')

    elif maskWholeBrain:
        # Mask the whole brain
        nacWholeBrainFile = "/Shared/paulsen/Experiments/rsFMRI/rs-fMRI-pilot/ReferenceAtlas/template_brain.nii.gz"

        warpWholeMask = rgwf.warpT1ToFMRI.clone('antsApplyTransformWholeBrain')
        warpWholeMask.inputs.invert_transform_flags = [True, False]
        warpWholeMask.inputs.input_image = nacWholeBrainFile
        preproc.connect(rgwf.forwardTransform, 'out', warpWholeMask, 'transforms')
        preproc.connect(rgwf.tstat, 'out_file', warpWholeMask, 'reference_image')

        wholeBrainMask = maskNode(name='wholeBrainMask', fileName='wholeBrainMask.nii', low=0.5, high=1.0, flags=['largest'])
        preproc.connect(warpWholeMask, 'output_image', wholeBrainMask, 'input_file')

        wholeMaskAvg = pipe.Node(interface=Maskave(), name='afni3DmaskAve_whole')
        wholeMaskAvg.inputs.outputtype = 'AFNI_1D'  # outputType
        wholeMaskAvg.inputs.args = '-median'  # TODO
        wholeMaskAvg.inputs.quiet = True
        preproc.connect(wholeBrainMask, 'output_file', wholeMaskAvg, 'mask')
        preproc.connect(calc, 'out_file', wholeMaskAvg, 'in_file')
        # 12
        deconvolve = afninodes.deconvolvenode(("Median_CSF", "Median_WM", "Median_WholeBrain"), "afni3Ddeconvolve")
        preproc.connect(csfAvg, 'out_file', deconvolve, 'stim_file_1')
        preproc.connect(wmAvg, 'out_file', deconvolve, 'stim_file_2')
        preproc.connect(wholeMaskAvg, 'out_file', deconvolve, 'stim_file_3')
        # preproc.connect(wholeBrainMask, 'out_file', deconvolve, 'mask')
        preproc.connect(volreg, 'oned_file', deconvolve, 'stim_file_4')

    else:
        deconvolve = afninodes.deconvolvenode(("Median_CSF", "Median_WM"), "afni3Ddeconvolve")
        preproc.connect(csfAvg, 'out_file', deconvolve, 'stim_file_1')
        preproc.connect(wmAvg, 'out_file', deconvolve, 'stim_file_2')
        preproc.connect(volreg, 'oned_file', deconvolve, 'stim_file_3')

    return preproc
