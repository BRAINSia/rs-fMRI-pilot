#!/usr/bin/env python
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

import os

from nipype.interfaces.base import (CommandLine, CommandLineInputSpec,
                                    InputMultiPath, traits, TraitedSpec,
                                    OutputMultiPath, isdefined,
                                    File, Directory)
from nipype.utils.filemanip import split_filename

class DWIconvertInputSpec(CommandLineInputSpec):
    conversionMode = traits.Enum('DicomToNrrd', 'DicomToFSL', 'NrrdToFSL', 'FSLToNrrd',
                                 argstr='--conversionMode %s', usedefault=True,
                                 desc='Determine which conversion to perform')
    outputVolume = traits.Str(argstr='--outputVolume %s')
    outputDirectory = Directory(os.path.getcwd(), exists=True, argstr='--outputDirectory %s',
                                usedefault=True, desc='Directory holding the output NRRD format')
    inputDicomDirectory = Directory(exists=True, argstr='--inputDicomDirectory %s',
                                    usedefault=False, desc='Directory holding Dicom series')


class DWIconvertOutputSpec(TraitedSpec):
    outputVolume = File(exists=True)


class DWIconvert(CommandLine):
    _cmd = 'ipldev/sharedopt/20120722/Darwin_i386/DicomToNrrdConverter/DWIConvert-build/DWIConvert'
    input_spec = DWIconvertInputSpec
    output_spec = DWIconvertOutputSpec

    def _format_arg(self, opt, spec, val):
        return super(convert, self)._format_arg(self, opt, spec, val)

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['outputVolume'] = os.path.join(self.inputs.outputDirectory, self.inputs.outputVolume)
        return outputs
