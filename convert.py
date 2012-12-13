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
    outputVolume = traits.Str('converted_dwi', argstr='%s', usedefault=True)
    outputDirectory = Directory(argstr='%s', desc='Directory holding the output file')
    inputDicomDirectory = Directory(exists=True, argstr='--inputDicomDirectory %s',
                                    usedefault=False, desc='Directory holding Dicom series')


class DWIconvertOutputSpec(TraitedSpec):
    outputVolume = File(exists=True)


class DWIconvert(CommandLine):
    _cmd = '/ipldev/sharedopt/20120722/Darwin_i386/DicomToNrrdConverter/DWIConvert-build/DWIConvert'
    input_spec = DWIconvertInputSpec
    output_spec = DWIconvertOutputSpec

    def _format_arg(self, opt, spec, val):
        if opt == 'outputVolume':
            return '--outputVolume %s' % self._gen_filename(val)
        elif opt == 'outputDirectory':
            return '--outputDirectory %s' % os.path.abspath(val)
        return super(DWIconvert, self)._format_arg(opt, spec, val)

    def _gen_filename(self, name):
        if len(name.split('.')) <= 1:
            if self.inputs.conversionMode == 'DicomToNrrd':
                name = name + '.nrrd'
            else:
                return super(DWIConvert, self)._gen_filename(self, name)
        return name

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['outputVolume'] = os.path.abspath(os.path.join(self.inputs.outputDirectory,
                                               self._gen_filename(self.inputs.outputVolume)))
        return outputs
