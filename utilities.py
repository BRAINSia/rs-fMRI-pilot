#!/usr/bin/env python
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

def makeList(input):
    """
    Hack to avoid making a MapNode for a list of length == 1
    """
    return [input]

def splitList(in_list):
    t1_out = in_list[0]
    label1_out = in_lits[1]
    label2_out = in_lits[2]
    return t1_out, label1_out, label2_out

def readNrrdHeader(fileName):
    """
    Grep out the number of slices and volumes from the Nrrd file created by DWIConvert
    """
    import re
    fID = open(fileName, 'r')
    try:
        for line in fID.readlines():
            search = re.search('(?<=sizes:\s)(?:[0-9]*)\s(?:[0-9])*\s(?P<slices>[0-9]*)\s(?P<volumes>[0-9]*)', line)
            if search is not None:
                slices = search.groupdict()['slices']
                volumes = search.groupdict()['volumes']
                return slices, volumes
    finally:
        fID.close()
    raise IOError('Nrrd file %s has regex match for slice and volume numbers!' % fileName)

def strCreate(slices, volumes, repTime):
    """ Removed FROM_IMAGE string in favor of hard-coded alt+z2 in afni.preprocess.py """
    return "%s %s %s" % (slices, volumes, repTime)

def dicomRead(infolder):
    """
    Use PyDicom to get the DICOM repetition time
    """
    import os
    import dicom
    infolder = os.path.abspath(infolder)
    for ii in os.listdir(infolder):
        meta = dicom.read_file(os.path.join(infolder, ii))
        if "RepetitionTime" in meta:
            return meta.data_element('RepetitionTime').value
    raise Exception('None of the dicom files in %s contain \
                     "Repetition Time" where PyDicom can find it!' % infolder)

def generateTissueMask(input_file, low=0.0, high=1.0, erodeFlag=True):
    """
    Using posterior tissue probability file, threshold by
    a low and/or high value, erode by 2mm, and take the ceil()
    Return the result as a binary mask file
    """
    import os
    import SimpleITK as sitk
    assert (isinstance(input_file[0], str)), input_file
    image = sitk.ReadImage(input_file[0])
    ## Compute brain mask
    inValue = 1
    outValue = 0
    binaryMask = sitk.BinaryThreshold(image, low, high, inValue, outValue)
    if erodeFlag:
        fileName = 'whiteMatterMask.nii'
        radiusMM = 1
        erodedMask = sitk.BinaryErode(binaryMask, radiusMM)
        sitk.WriteImage(erodedMask, os.path.abspath('eroded_' + fileName))
        connected = sitk.ConnectedComponent(erodedMask)
        sortedComp = sitk.RelabelComponent(connected)
        maskOnly = sitk.BinaryThreshold(sortedComp, 0.99, 1.0, inValue, outValue)
    else:
        fileName = 'csfMask.nii'
        connected = sitk.ConnectedComponent(binaryMask)
        sortedComp = sitk.RelabelComponent(connected)
        #        maskOnly = sitk.BinaryThreshold(sortedComp, 0.0, 1.0, inValue, outValue)
        maskOnly = sitk.BinaryThreshold(sortedComp, 0.99, 1.0, inValue, outValue)
    sitk.WriteImage(binaryMask, os.path.abspath('binary_' + fileName))
    sitk.WriteImage(connected, os.path.abspath('connected_' + fileName))
    sitk.WriteImage(sortedComp, os.path.abspath('sorted_' + fileName))
    sitk.WriteImage(maskOnly, os.path.abspath(fileName))
    return os.path.abspath(fileName)
