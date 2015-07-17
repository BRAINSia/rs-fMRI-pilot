import os.path

import nipype.pipeline.engine as pipe
from nipype.interfaces.io import DataSink, DataGrabber

# Templates
autoworkup = 'Experiments/{experiment}/{site}/{subject}/%s/{result}/{fname}.nii.gz'
raw_dicom = 'MRx/{site}/{subject}/%s/%s/%s/*'
cleveland = '20150409_ClevelandPreProcessedfMRI/%s/*.hdr'
transform = '{experiment}/SubjectToAtlasWarped/%s/%s.h5'


def datagrabber(infields, template, field_template, template_args, base_directory, name):
    assert field_template.keys() == template_args.keys(), "Template keys do not match!"
    grabber = pipe.Node(interface=DataGrabber(infields=infields, outfields=field_template.keys()), name=name)
    grabber.inputs.sort_filelist = False
    grabber.inputs.base_directory = base_directory
    grabber.inputs.template = template
    grabber.inputs.field_template = field_template
    grabber.inputs.template_args = template_args
    return grabber


def datasink(base_directory, name, container=None, overwrite=False):
    output = pipe.Node(interface=DataSink(), name=name)
    output.inputs.base_directory = base_directory
    if container is not None:
        output.inputs.container = container
    output.inputs.parameterization = False
    output.overwrite = overwrite
    return output


def iowaGrabber(experiment, site, subject, maskGM=False, maskWholeBrain=False, **kwargs):
    kwargs["name"] = "iowaGrabber"
    kwargs["field_template"] = dict()
    kwargs["template_args"] = dict()
    # get white matter
    posterior = autoworkup.format(experiment=experiment, site=site, subject=subject, result='ACCUMULATED_POSTERIORS',
                                  fname='POSTERIOR_%s_TOTAL')
    kwargs["field_template"]['whmFile'] = posterior
    kwargs["template_args"]['whmFile'] = [['session_id', 'WM']]
    # get fMRI DICOMs
    kwargs["field_template"]['fmri_dicom_dir'] = raw_dicom.format(site=site, subject=subject)
    kwargs["template_args"]['fmri_dicom_dir'] = [['session_id', 'ANONRAW', 'FMRI_RestingStateConnectivity']]
    if maskGM:
        kwargs["field_template"]['gryFile'] = posterior
        kwargs["template_args"]['gryFile'] = [['session_id', 'GM']]
        # For cerebrum gray matter ONLY:
        # kwargs["field_template"]['gryFile'] = tissuecls
        # kwargs["template_args"]['gryFile'] = [['session_id', 'POSTERIOR_SURFGM']]
    elif maskWholeBrain:
        pass  # No need for grabber, we're using NAC-atlas file
    return autoworkupGrabber(experiment, site, subject, **kwargs)


def autoworkupGrabber(experiment, site, subject, **kwargs):
    kwargs.setdefault("name", "t1Grabber")
    kwargs.setdefault("base_directory", "/Shared/paulsen")  #Largest common path for both raw_dicom and autoworkup templates
    kwargs.setdefault("template", "*")
    kwargs.setdefault("infields", ['session_id'])
    ftemplate = autoworkup.format(experiment=experiment, site=site, subject=subject, result="TissueClassify", fname="%s")
    if kwargs.has_key("field_template"):  # called from iowaGrabber
        kwargs["field_template"]["t1_File"] = ftemplate
        kwargs["field_template"]["csfFile"] = ftemplate
        kwargs["template_args"]["t1_File"] = [['session_id', 't1_average_BRAINSABC']]
        kwargs["template_args"]["csfFile"] = [['session_id', 'complete_brainlabels_seg']]  #'fixed_brainlabels_seg']])
    else:
        kwargs["field_template"] = dict(t1_File=ftemplate, csfFile=ftemplate)
        kwargs["template_args"] = dict(t1_File=[['session_id', 't1_average_BRAINSABC']],
                                       csfFile=[['session_id', 'complete_brainlabels_seg']])  #'fixed_brainlabels_seg']])
    return datagrabber(**kwargs)


def clevelandGrabber(**kwargs):
    kwargs.setdefault("name", "clevelandGrabber")
    kwargs.setdefault("base_directory", "/Shared/paulsen/Experiments")
    kwargs.setdefault("infields", ['session_id'])
    kwargs.setdefault("template", "*")
    kwargs.setdefault("field_template", dict(fmriHdr=cleveland))
    kwargs.setdefault("template_args", dict(fmriHdr=[['session_id']]) )
    return datagrabber(**kwargs)


def transformGrabber(experiment):
    grabber = pipe.Node(interface=DataGrabber(infields=["session_id"],
                                              outfields=["atlasToSessionTransform", "sessionToAtlasTransform"]),
                        name="transformGrabber")
    grabber.inputs.sort_filelist = False
    grabber.inputs.base_directory = "/Shared/paulsen/Experiments"
    grabber.inputs.template = "*"
    transformRegex = transform.format(experiment=experiment)
    grabber.inputs.field_template = dict(atlasToSessionTransform=transformRegex,
                                         sessionToAtlasTransform=transformRegex)
    grabber.inputs.template_args = dict(atlasToSessionTransform=[["session_id", "AtlasToSubjectPostBABC_SyNComposite"]],
                                        sessionToAtlasTransform=[["session_id", "AtlasToSubjectPostBABC_SyNInverseComposite"]])
    return grabber


def atlasSink(**kwargs):
    kwargs['container'] = kwargs.pop('name') + "_Results"  # The 'name' from the command line
    kwargs['name'] = "atlas_DataSink"
    kwargs['overwrite'] = kwargs.pop('force')
    keys = kwargs.keys()
    for key in keys:
        if key not in ['base_directory', 'container', 'name', 'overwrite']:
            kwargs.pop(key)
    return datasink(**kwargs)


def fmriSink(result_directory, **kwargs):
    # 'container' must be reserved for <SESSION>, so we set the base_directory to the
    kwargs['base_directory'] = os.path.join(result_directory, kwargs['name'] + "_Results")
    kwargs['name'] = "fmri_label_DataSink"
    kwargs['overwrite'] = kwargs.pop('force')
    keys = kwargs.keys()
    for key in keys:
        if key not in ['base_directory', 'name', 'overwrite']:
            kwargs.pop(key)

    return datasink(**kwargs)
