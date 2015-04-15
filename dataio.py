import nipype.pipeline.engine as pipe
from nipype.interfaces.io import DataSink, DataGrabber

def iowaGrabber(infields, template, template_args, base_directory):
    assert template.keys() == template_args.keys(), "Template keys do not match!"
    grabber = pipe.Node(interface=DataGrabber(infields=infields, outfields=template.keys()), name='dataGrabber')
    grabber.inputs.base_directory = base_directory
    grabber.inputs.template = '*'
    grabber.inputs.field_template = template
    grabber.inputs.template_args = template_args
    return grabber


def clevelandGrabber(infields, outfields, base_directory, experiment):
    grabber = pipe.Node(interface=DataGrabber(infields=infields,
                                              outfields=outfields), name='dataGrabber')
    grabber.inputs.base_directory = base_directory
    grabber.inputs.template = '*'

    fmriRegex = '{experiment}/%s/Y%s/day%s/*.hdr'.format(experiment="20150409_ClevelandPreProcessedfMRI")
    t1Regex = '{experiment}/FMRI_HD_120/%s/%s/TissueClassify/ti_average_BRAINSABC.nii.gz'.format(experiment=experiment)
    field_template = dict(fmriHdr=fmriRegex,
                          t1_File=t1Regex)
    template_args = dict(fmriHdr=[['subject_id', 'year', 'day']],
                         t1_File=[['subject_id', 'session_id']])
    grabber.inputs.field_template = field_template
    grabber.inputs.template_args = template_args
    return grabber


def transformGrabber(experiment):
    grabber = pipe.Node(interface=DataGrabber(infields=['session_id'],
                                              outfields=['atlasToSessionTransform', 'sessionToAtlasTransform']),
                        name='transformGrabber')
    grabber.inputs.base_directory = '/Shared/paulsen/Experiments/{experiment}/SubjectToAtlasWarped'.format(experiment=experiment)
    grabber.inputs.template = '*'
    transformRegex = '%s/AtlasToSubject_%sComposite.h5'
    grabber.inputs.field_template = dict(atlasToSessionTransform=transformRegex,
                                         sessionToAtlasTransform=transformRegex)
    grabber.inputs.template_args = dict(atlasToSessionTransform=[['session_id', '']],
                                        sessionToAtlasTransform=[['session_id', 'Inverse']])
    return grabber


def datasink(base_directory, container, name, overwrite=False):
    output = pipe.Node(interface=DataSink(), name=name)
    output.inputs.base_directory = base_directory
    output.inputs.container = container
    output.inputs.parameterization = False
    output.overwrite = overwrite
    return output
