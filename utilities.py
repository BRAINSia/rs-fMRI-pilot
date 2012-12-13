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
    image = sitk.ReadImage(input_file)
    ## Compute brain mask
    inValue = 1
    outValue = 0
    binaryMask = sitk.BinaryThreshold(image, low, high)
    if erodeFlag:
        fileName = 'whiteMatterMask.nii'
        radiusMM = 1
        erodedMask = sitk.BinaryErode(binaryMask, radiusMM)
        sitk.WriteImage(erodedMask, os.path.abspath('eroded_' + fileName))
        connected = sitk.ConnectedComponent(erodedMask)
        sortedComp = sitk.RelabelComponent(connected, 10) # HACK
        maskOnly = sitk.BinaryThreshold(sortedComp, 1, 1)
    else:
        fileName = 'csfMask.nii'
        connected = sitk.ConnectedComponent(binaryMask)
        sortedComp = sitk.RelabelComponent(connected)
        #        maskOnly = sitk.BinaryThreshold(sortedComp, 0.0, 1.0, inValue, outValue)
        maskOnly = sitk.BinaryThreshold(sortedComp, 2.0, 2.0)
    sitk.WriteImage(binaryMask, os.path.abspath('binary_' + fileName))
    sitk.WriteImage(connected, os.path.abspath('connected_' + fileName))
    sitk.WriteImage(sortedComp, os.path.abspath('sorted_' + fileName))
    sitk.WriteImage(maskOnly, os.path.abspath('final_' + fileName))
    return os.path.abspath('final_' + fileName)

def getLabelList(label_file, arg_template):
    import os
    import SimpleITK as sitk
    assert os.path.exists(label_file)
    label = sitk.ReadImage(label_file)
    labelStat = sitk.LabelStatisticsImageFilter()
    casted = sitk.Cast(label, sitk.sitkUInt16)
    labelStat.Execute(casted, casted)
    labelList = labelStat.GetValidLabels()
    arg_str = []
    for labelCode in labelList:
        if labelCode > 0:
            arg_str.append(arg_template.format(labelCode))
    return arg_str, labelList

def writeMat(in_file):
    """
    Pass output of afni.preprocess.ROIStat() to be written as .mat file

    Nota bene: ROIStats() should only output ONE AND ONLY ONE
    Based on code on Stackflow:
           (http://stackoverflow.com/questions/1095265/matrix-from-python-to-matlab)
    """
    import os
    import numpy as np
    from scipy.io import savemat
    with open(in_file, 'r') as fID:
        count = 0
        timepoint = []
        row_data = []
        for valueStr in fID.readlines():
            value_list = valueStr.split()
            if count == 0:
                column_header = value_list[2:]
                valid_labels = [item.strip('NZMed_') for item in column_header]
                columns = len(column_header)
            else:
                timepoint.append(value_list[1].split('[')[0])
                row_data.append(value_list[2:])
            count += 1
        rows = len(timepoint)
        data = np.array(row_data, dtype='float', ndmin=2)
        labels = np.array(valid_labels, dtype='int', ndmin=1)
    corr = np.corrcoef(data)
    temp, _ = os.path.basename(in_file).split('.')
    savemat(file_name=temp + '_corr', mdict={'data':corr}, appendmat=True, oned_as='row')
    savemat(file_name=temp + '_labels', mdict={'data':labels}, appendmat=True, oned_as='row')
    savemat(file_name=temp + '_raw', mdict={'data':data}, appendmat=True, oned_as='row')
    return temp + '_corr.mat'

def generateMatSubstitution(in_file, session):
    import os.path
    oldString = os.path.basename(in_file)
    replaceString = "{session}_{suffix}".format(session=session,
                                                suffix=oldString.split('roiStat_')[1])
    return [(oldString, replaceString)]
