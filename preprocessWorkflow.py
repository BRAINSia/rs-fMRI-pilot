import nipype.pipeline.engine as pipe
from nipype.interfaces.utility import Function, IdentityInterface

import afninodes
import utilities
import nuisanceWorkflow


def prepWorkflow(skipCount, outputType, name="prep"):
    preprocessing = pipe.Workflow(name=name)

    # Nodes
    # inputnode = pipe.Node(interface=IdentityInterface(fields=in_fields), name='inputs')
    format_out = ['modality', 'numberOfSlices', 'numberOfFiles', 'repetitionTime', 'sliceOrder']
    formatFMRINode = pipe.Node(interface=Function(function=utilities.formatFMRI,
                                                  input_names=['dicomDirectory'],
                                                  output_names=format_out),
                               name='formatFMRINode')

    to_3D_str = afninodes.to3dstrnode('strCreate')
    to_3D = afninodes.to3dnode('to_3D')
    refit = afninodes.refitnode('refit')
    despike = afninodes.despikenode(outputType, skipCount, 'despike')
    volreg = afninodes.volregnode(outputType, 'volreg')
    zeropad = afninodes.zeropadnode('zeropad')
    merge = afninodes.mergenode(outputType, 'merge')
    automask = afninodes.automasknode(outputType, 'automask')
    calc = afninodes.multiplynode(outputType, 'calc')

    def strToIntMinusOne(string):
        return int(string) - 1

    preprocessing.connect([(formatFMRINode, to_3D_str, [('numberOfSlices', 'slices'),
                                                        ('numberOfFiles', 'volumes'),
                                                        ('repetitionTime', 'repTime'),
                                                        ('sliceOrder', 'order')]),
                           (to_3D_str, to_3D,          [('funcparams', 'funcparams')]),
                           (formatFMRINode, despike,   [(('numberOfFiles', strToIntMinusOne), 'end')]),
                           (to_3D, refit,              [('out_file', 'in_file')]),    # 1a
                           (refit, despike,            [('out_file', 'in_file')]),    # 2
                           (despike, volreg,           [('out_file', 'in_file')]),    # 3
                           (volreg, zeropad,           [('out_file', 'in_file')]),    # 4
                           (zeropad, merge,            [('out_file', 'in_files')]),   # 5
                           (merge, automask,           [('out_file', 'in_file')]),    # 6
                           (merge, calc,               [('out_file', 'in_file_a')]),  # 7
                           (automask, calc,            [('out_file', 'in_file_b')]),  # 8
                           ])
    return preprocessing


def workflow(skipCount, outputType, name, **kwargs):
    """ preprocessing workflow

    Connections needed:
      (in)
      formatFMRI.dicomDirectory
      to_3D.infolder

      (out)
      merge.out_file
      automask.out_file
      calc.out_file

    """
    # Nodes
    prep = prepWorkflow(skipCount, outputType, name="prep")
    nuisance = nuisanceWorkflow.workflow(outputType=outputType, name="nuisance", **kwargs)

    master = pipe.Workflow(name=name)

    master.connect([(prep, nuisance, [('calc.out_file', 'wm.afni3DmaskAve_wm.in_file'),
                                      ('calc.out_file', 'csf.afni3DmaskAve_csf.in_file'),
                                      ('calc.out_file', 'afni3Ddeconvolve.in_file')])])
    if kwargs['maskgm']:
        master.connect([(prep, nuisance, [('calc.out_file', 'gm.afni3DmaskAve_grm.in_file'),
                                          ('volreg.oned_file', 'afni3Ddeconvolve.stim_file_4')])])
    elif kwargs['maskwb']:
        master.connect([(prep, nuisance, [('calc.out_file', 'wb.afni3DmaskAve_whole.in_file'),
                                          ('volreg.oned_file', 'afni3Ddeconvolve.stim_file_4')])])
    else:
        master.connect([(prep, nuisance, [('volreg.oned_file', 'afni3Ddeconvolve.stim_file_3')])])
    return master
