#!/usr/bin/env python
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

def writeMat(data, fileName, key='data'):
    """
    Pass an array of data to be written to a .mat file
    Based on code on Stackflow:
           (http://stackoverflow.com/questions/1095265/matrix-from-python-to-matlab)
    """
    import os
    from scipy.io import savemat

    appendmat = True
    if fileName[-4:] == '.mat':
        appendmat = False
    fileName = os.path.abspath(fileName)
    if isinstance(data, dict):
        savemat(file_name=fileName, mdict=data, appendmat=appendmat, oned_as='row')
    savemat(file_name=fileName, mdict={key:data}, appendmat=appendmat, oned_as='row')

def registerLabelsTofMRI(label_file, fMRI_file):
    import os
    import SimpleITK as sitk

    labelImage = sitk.ReadImage(label_file)
    labelStat = sitk.LabelStatistics(labelImage)
    labelMin = labelStat.GetMinimum()
    labelMax = labelStat.GetMaximum()
    fmriImage = sitk.ReadImage(fMRI_file)
    size = new_size = list(fmriImage.GetSize())
    new_size[3] = 0
    for vIndex in range(size[3]):
        index = [0, 0, 0, vIndex]
        extractor = sitk.ExtractImageFilter()
        extractor.SetSize(new_size)
        extractor.SetIndex(index)




def createDataDictionary():
    pass
