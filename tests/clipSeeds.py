#!/usr/bin/env python

import os
import sys

import SEMTools as sem
from nipype.interfaces.ants.registration import Registration
from nipype.interfaces.ants import ApplyTransforms
from nipype.interfaces.afni.preprocess import *
from nipype.interfaces.freesurfer.preprocess import *
from nipype.interfaces.io import DataSink, DataGrabber
from nipype.interfaces.utility import Function, IdentityInterface, Rename, Select
from nipype.interfaces.utility import Merge as Merger
import nipype.pipeline.engine as pipe
import numpy

from utilities import *

def writeSeedFiles():
    CACHE_DIR = 'seeds_CACHE'
    RESULTS_DIR = 'seeds'
    REWRITE_DATASINKS = True
    nacAtlasFile = "/Shared/paulsen/Experiments/rsFMRI/rs-fMRI-pilot/ReferenceAtlas/template_t1.nii.gz"
    nacAtlasLabel = "/Shared/paulsen/Experiments/rsFMRI/rs-fMRI-pilot/ReferenceAtlas/template_nac_labels.nii.gz"
    nacResampleResolution = (2.0, 2.0, 2.0)
    downsampledNACfilename = 'downsampledNACatlas.nii.gz'

    preproc = pipe.Workflow(updatehash=True, name=CACHE_DIR)
    preproc.base_dir = os.getcwd()

    labels, seeds = getAtlasPoints('/Shared/paulsen/Experiments/rsFMRI/rs-fMRI-pilot/seeds.fcsv')

    seedsIdentity = pipe.Node(interface=IdentityInterface(fields=['index']),
                              name='seedsIdentity')
    seedsIdentity.iterables = ('index', range(len(labels)))

    selectSeed = pipe.Node(interface=Select(), name='selectSeed')
    selectSeed.inputs.inlist = seeds
    preproc.connect(seedsIdentity, 'index', selectSeed, 'index')

    selectLabel = pipe.Node(interface=Select(), name='selectLabel')
    selectLabel.inputs.inlist = labels
    preproc.connect(seedsIdentity, 'index', selectLabel, 'index')

    points = pipe.Node(interface=Function(function=createSphereExpression,
                                          input_names=['coordinates', 'radius'],
                                          output_names=['expression']),
                       name='createSphereExpression')
    preproc.connect(selectSeed, 'out', points, 'coordinates')

    downsampleAtlas = pipe.Node(interface=Function(function=resampleImage,
                                                   input_names=['inputVolume', 'outputVolume', 'resolution'],
                                                   output_names=['outputVolume']),
                                name="downsampleAtlas")
    downsampleAtlas.inputs.inputVolume = nacAtlasFile
    downsampleAtlas.inputs.outputVolume = downsampledNACfilename
    downsampleAtlas.inputs.resolution = [int(x) for x in nacResampleResolution]

    spheres = pipe.Node(interface=Calc(letters=['a']),
                        name='afni3Dcalc_seeds')
    spheres.inputs.outputtype = 'NIFTI_GZ'
    preproc.connect(downsampleAtlas, 'outputVolume', spheres, 'in_file_a')
    spheres.inputs.args = '-nscale'

    preproc.connect(points, 'expression', spheres, 'expr')

    renameMasks = pipe.Node(interface=Rename(format_string='%(label)s_mask'), name='renameMasksAtlas')
    renameMasks.inputs.keep_ext = True
    preproc.connect(selectLabel, 'out', renameMasks, 'label')
    preproc.connect(spheres, 'out_file', renameMasks, 'in_file')

    atlas_DataSink = pipe.Node(interface=DataSink(), name="atlas_DataSink")
    atlas_DataSink.inputs.base_directory = preproc.base_dir  # '/Shared/paulsen/Experiments/20130417_rsfMRI_Results'
    atlas_DataSink.inputs.container = RESULTS_DIR
    atlas_DataSink.inputs.parameterization = False
    atlas_DataSink.overwrite = REWRITE_DATASINKS
    preproc.connect(renameMasks, 'out_file', atlas_DataSink, 'Atlas')
    preproc.connect(downsampleAtlas, 'outputVolume', atlas_DataSink, 'Atlas.@resampled')
    preproc.run()


def test_pipeline():
    maskWhiteMatterFromSeeds = True
    CACHE_DIR = 'seed_Test_cache'
    preproc = pipe.Workflow(updatehash=True, name=CACHE_DIR)
    preproc.base_dir = os.getcwd()

    labelFile = "/Shared/paulsen/Experiments/rsFMRI/rs-fMRI-pilot/ReferenceAtlas/template_ABC_labels.nii.gz"
    downsampledLabel = 'template_downsampled.nii.gz'
    nacResampleResolution = (2.0, 2.0, 2.0)
    downsampleAtlas = pipe.Node(interface=Function(function=resampleImage,
                                                   input_names=['inputVolume', 'outputVolume', 'resolution'],
                                                   output_names=['outputVolume']),
                                name="downsampleAtlas")
    downsampleAtlas.inputs.inputVolume = labelFile
    downsampleAtlas.inputs.outputVolume = downsampledLabel
    downsampleAtlas.inputs.resolution = [int(x) for x in nacResampleResolution]

    grabber = pipe.Node(interface=DataGrabber(infields=['label'], outfields=['seedfile']), name='dataGrabber')
    grabber.inputs.base_directory = os.path.abspath('seeds/Atlas')
    grabber.inputs.template = '%s_mask.nii.gz'

    labels, seeds = getAtlasPoints('/Shared/paulsen/Experiments/rsFMRI/rs-fMRI-pilot/seeds.fcsv')
    print labels
    seedsIdentity = pipe.Node(interface=IdentityInterface(fields=['index']),
                              name='seedsIdentity')
    seedsIdentity.iterables = ('index', range(len(labels)))

    labelNode = pipe.Node(interface=IdentityInterface(fields=['label']), name='labelNode')
    labelNode.iterables = ('label', labels)
    preproc.connect(labelNode, 'label', grabber, 'label')

    # import glob
    # seedfiles = [os.path.abspath(x) for x in glob.glob('seeds/Atlas/*mask.nii.gz')]
    # infiles.iterables = ('infile', seedfiles)

    atlas_DataSink = pipe.Node(interface=DataSink(), name="atlas_DataSink")
    atlas_DataSink.inputs.base_directory = preproc.base_dir  # '/Shared/paulsen/Experiments/20130417_rsfMRI_Results'
    atlas_DataSink.inputs.container = 'seed_Test'
    atlas_DataSink.inputs.parameterization = False
    atlas_DataSink.overwrite = True

    clipSeedWithVentriclesNode = pipe.Node(interface=Function(function=clipSeedWithVentricles,
                                           input_names=['unclipped_seed_fn', 'fmriBABCSeg_fn', 'desired_out_seed_fn'],
                                           output_names=['clipped_seed_fn']),
                                           name='clipSeedWithVentriclesNode')
    def create_label(label, _type):
        return '{0}_{1}_clipped.nii.gz'.format(label, _type)
    preproc.connect([(labelNode, clipSeedWithVentriclesNode, [(('label', create_label, 'csf'), 'desired_out_seed_fn')])])
    preproc.connect(downsampleAtlas, 'outputVolume', clipSeedWithVentriclesNode, 'fmriBABCSeg_fn')

    preproc.connect(grabber, 'seedfile', clipSeedWithVentriclesNode, 'unclipped_seed_fn')
    # preproc.connect(infiles, 'infile', clipSeedWithVentriclesNode, 'desired_out_seed_fn')
    preproc.connect(clipSeedWithVentriclesNode, 'clipped_seed_fn', atlas_DataSink, 'CSF')

    if maskWhiteMatterFromSeeds:
        clipSeedWithWhiteMatterNode = pipe.Node(interface=Function(function=clipSeedWithWhiteMatter,
                                                                   input_names=['seed', 'mask', 'outfile'],
                                                                   output_names=['outfile']),
                                                name='clipSeedWithWhiteMatterNode')
        preproc.connect([(labelNode, clipSeedWithWhiteMatterNode, [(('label', create_label, 'wm'), 'outfile')])])
        preproc.connect(downsampleAtlas, 'outputVolume', clipSeedWithWhiteMatterNode, 'mask')

        preproc.connect(clipSeedWithVentriclesNode, 'clipped_seed_fn', clipSeedWithWhiteMatterNode, 'seed')
        preproc.connect(clipSeedWithWhiteMatterNode, 'outfile', atlas_DataSink, 'WM')

    preproc.run()


if __name__ == '__main__':
    if not os.path.exists('seeds/Atlas'):
        writeSeedFiles()
    test_pipeline()
