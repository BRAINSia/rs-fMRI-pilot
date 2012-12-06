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

# def getCovarianceMatrix(fmri_file, label_file):
#     """
#     Read in a 4D nifti file and a Freesurfer label file and compute the median value
#     for all regions in all volumes and compute the covariance matrix for those values
#     """
#     import SimpleITK as sitk
#     def seperate_volumes(filename):
#         """
#         Write the seperate volumes of a 4D (3D + time) file and return a list of the filenames
#         """
#         from nibabel import nifti1 as nifti
#         import numpy as np
#         import os.path
#         assert isinstance(index, int)
#         image = nifti.load(filename)
#         fourD = image.get_data()
#         fileList = []
#         for index in range(len(fourD[...,:])):
#             header = image.get_header()
#             affine = header.get_best_affine()
#             threeD = fourD[...,index]
#             newVolume = nifti.Nifti1Image(data=threeD, affine=affine, header=header)
#             newVolume.update_header()
#             newFile = 'threeD_%d.nii' % index
#             newVolume.set_filename(newFile)
#             newVolume.to_files()
#             assert os.path.exists(os.path.abspath(newFile))
#             fileList.append(os.path.abspath(newFile))
#         return fileList
#     def writeMat(out_file, labelList, fileNames):
#         """
#         Pass an array of data to be written to a .mat file
#         Based on code on Stackflow:
#                (http://stackoverflow.com/questions/1095265/matrix-from-python-to-matlab)
#         """
#         import os
#         import numpy as np
#         from scipy.io import savemat
#         columns = len(labelList)
#         rows = len(fileNames)
#         data = np.array(np.empty((rows, columns)), ndmin=2)
#         index = 0
#         for fileName in fileNames:
#             with open(fileName, 'r') as fID:
#                 valueStr = fID.readlines()
#                 temp = [float(item.rstrip('\n')) for item in valueStr]
#                 np.put(data, index, np.array(temp, ndmin=2))
#                 index += columns
#         covar = np.cov(data)
#         if out_file[-4:] == '.mat':
#             out_file = out_file[:-4]
#         fileName = os.path.abspath(out_file)
#         for fname in [fileName, fileName + 'labels']:
#             if isinstance(covar, dict):
#                 savemat(file_name=fname, mdict=covar, appendmat=True, oned_as='row')
#             savemat(file_name=fname, mdict={'data':covar}, appendmat=True, oned_as='row')
#         fileName = fileName + '.mat'
#         return fileName
#     label = sitk.Cast(sitk.ReadImage(label_file), sitk.sitkUInt16)
#     labelStat = sitk.LabelStatisticsImageFilter()
#     volArray = [[]]
#     for volFile in seperate_volumes(fmri_file):
#         newVol = sitk.ReadImage(volFile)
#         newVol.GetPixelIDTypeAsString()
#         labelStat.Execute(newVol, label)
#         labelList = []
#         for lb in labelStat.GetValidLabels():
#             labelList.append(labelStat.GetMedian(lb))
#         volArray.append(labelList)
#     import numpy as np
#     rows =
#     columns = len(labelList)
#     data = np.array(np.empty((rows, columns)), ndmin=2)
#     index = 0
#     for row in range(rows):
#         numpy.put( data, index, np.array(labelList, ndmin=2))
#         index += columns
#     frow, fcolumn =  data.shape
#     covar = np.cov(data)

def writeMat(out_file, fileNames, volumeCount, skippedVolCount):
    """
    Pass an array of data to be written to a .mat file
    Based on code on Stackflow:
           (http://stackoverflow.com/questions/1095265/matrix-from-python-to-matlab)
    """
    import os
    import numpy as np
    from scipy.io import savemat
    rows = len(fileNames)
    columns = int(volumeCount) - skippedVolCount # - 4 ## HACK: "- 4""
    data = np.array(np.empty((rows, columns)), ndmin=2)
    index = 0
    for fileName in fileNames:
        with open(fileName, 'r') as fID:
            valueStr = fID.readlines()
            temp = [item.rstrip('\n') for item in valueStr]
            np.put(data, range(index, index + columns), temp)
            for ii in range(len(temp)):
                assert float(temp[ii]) == data.flat[index+ii], \
                  "Not equal!, Index %d, temp %f, data %f" % (ii, float(temp[ii]), data.flat[index+ii])
            index += columns
    # covar = np.cov(data.T)
    corr = np.corrcoef(data.T)
    if out_file[-4:] == '.mat':
        out_file = out_file[:-4]
    fileName = os.path.abspath(out_file)
    valid_labels = [2, 4, 5, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 18, 24, 26, 28, 30, 31, 41, 43, 44, 46, 47, 49, 50, 51, 52,
                    53, 54, 58, 60, 62, 63, 77, 85, 251, 252, 253, 254, 255, 1000, 1001, 1002, 1003, 1005, 1006, 1007, 1008,
                    1009, 1010, 1011, 1012, 1013, 1014, 1015, 1016, 1017, 1018, 1019, 1020, 1021, 1022, 1023, 1024, 1025, 1026,
                    1027, 1028, 1029, 1030, 1031, 1032, 1033, 1034, 1035, 2000, 2001, 2002, 2003, 2005, 2006, 2007, 2008, 2009,
                    2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026, 2027,
                    2028, 2029, 2030, 2031, 2032, 2033, 2034, 2035]
    # for fname in [fileName]:
        # if isinstance(corr, dict):
        #     savemat(file_name=fname, mdict=corr, appendmat=True, oned_as='row')
    savemat(file_name=fileName, mdict={'data':corr}, appendmat=True, oned_as='row')
    savemat(file_name='raw_values', mdict={'data':data}, appendmat=True, oned_as='row')
    savemat(file_name='validLabels', mdict={'data':valid_labels}, appendmat=True, oned_as='row')

    fileName = fileName + '.mat'
    return fileName
