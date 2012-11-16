#!/usr/bin/env python
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

def writeMat(out_file, labelList, fileNames):
    """
    Pass an array of data to be written to a .mat file
    Based on code on Stackflow:
           (http://stackoverflow.com/questions/1095265/matrix-from-python-to-matlab)
    """
    import os
    import numpy as np
    from scipy.io import savemat
    columns = len(labelList)
    rows = len(fileNames)
    data = np.array(np.empty((rows, columns)), ndmin=2)
    index = 0
    for fileName in fileNames:
        with open(fileName, 'r') as fID:
            valueStr = fID.readlines()
            temp = [float(item.rstrip('\n')) for item in valueStr]
            np.put(data, index, np.array(temp, ndmin=2))
            index += columns
    covar = np.cov(data)
    if out_file[-4:] == '.mat':
        out_file = out_file[:-4]
    fileName = os.path.abspath(out_file)
    for fname in [fileName, fileName + 'labels']:
        if isinstance(covar, dict):
            savemat(file_name=fname, mdict=covar, appendmat=True, oned_as='row')
        savemat(file_name=fname, mdict={'data':covar}, appendmat=True, oned_as='row')
    fileName = fileName + '.mat'
    return fileName
