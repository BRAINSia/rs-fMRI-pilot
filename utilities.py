#!/usr/bin/env python
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

def makeList(input):
    """
    Hack to avoid making a MapNode for a list of length == 1
    """
    return [input]

def concatTransforms(transform1, transform2):
    """
    Hack to put two transforms into a list
    """
    return [transform1, transform2]

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

def strCreate(time, slices, volumes, repTime, order):
    """
    Create a funcparams string for To3D()
    """
    return "%s %s %s %s %s" % (time, slices, volumes, repTime, order.strip())

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
        sitk.WriteImage(binaryMask, os.path.abspath('binary_' + fileName))
        sitk.WriteImage(connected, os.path.abspath('connected_' + fileName))
        sitk.WriteImage(sortedComp, os.path.abspath('sorted_' + fileName))
        sitk.WriteImage(maskOnly, os.path.abspath('final_' + fileName))
    else:
        fileName = 'csfMask.nii'
        # csfLabels = [4,23]
        image1 = sitk.BinaryThreshold(image, 3, 4)
        image2 = sitk.BinaryThreshold(image, 42, 43)
        final_image = image1 + image2
        sitk.WriteImage(final_image, os.path.abspath('final_' + fileName))
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

    Nota bene: ROIStats() should only output ONE AND ONLY ONE file
        The format of the input file is:
        --------------------------
        File  Sub-brick  NZMed_2  NZMed_4 ...
        <filename>  0[#4] -0.337481	 0.358179  ...
        <filename>  0[#4] 1.524675	 0.230819  ...
        ...
        --------------------------

    Based on code from Stackflow:
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
                # Header line
                column_header = value_list[2:]
                regionsLabel = [item.strip('NZMed_') for item in column_header]
                # columns = len(column_header)
                count += 1
            else:
                # timepoint.append(value_list[1].split('[')[0]) # Sub-brick
                row_data.append(value_list[2:])
        #rows = len(timepoint)
        data = np.array(row_data, dtype='float', ndmin=2)
        labels = np.array(regionsLabel, dtype='int', ndmin=1)
    corr = np.corrcoef(data, rowvar=0) # Set rowvar > 0 to correlate timepoints, rowvar = 0 for regions
    temp, _ = os.path.basename(in_file).split('.')
    corrFile = os.path.abspath(temp + '_corr.mat')
    labelFile = os.path.abspath(temp + '_labels.mat')
    rawFile = os.path.abspath(temp + '_raw.mat')
    savemat(file_name=corrFile, mdict={'data':corr}, appendmat=True, oned_as='row')
    savemat(file_name=labelFile, mdict={'data':labels}, appendmat=True, oned_as='row')
    savemat(file_name=rawFile, mdict={'data':data}, appendmat=True, oned_as='row')
    return corrFile, labelFile, rawFile

def generateMatSubstitution(in_file, session):
    import os.path
    oldString = os.path.basename(in_file)
    replaceString = "{session}_{suffix}".format(session=session,
                                                suffix=oldString.split('roiStat_')[1])
    return [(oldString, replaceString)]

def getAtlasPoints(filename):
    """ Assumes file with text header:

    """
    from csv import reader, DictReader

    class FiducialReader(DictReader):
        def __init__(self, fid, commentchar='#', *args, **kwds):
            if issubclass(DictReader, object):
                super(DictReader, self).__init__(fid, *args, **kwds)
            else:
                DictReader.__init__(self, fid, *args, **kwds)
            self.commentchar = commentchar
            self.leadingfield = self.commentchar + 'label'

        def __iter__(self):
            return self

        @property
        def fieldnames(self):
            while self._fieldnames is None or self._fieldnames[0] != self.leadingfield:
                try:
                    self._fieldnames = self.reader.next()
                except StopIteration:
                    pass
            self.line_num = self.reader.line_num
            return self._fieldnames

        @fieldnames.setter
        def fieldnames(self, value):
            self._fieldnames = values

        def next(self):
            if self.line_num == 0:
               # Used only for its side effect.
               self.fieldnames
            row = self.reader.next()
            self.line_num = self.reader.line_num
            # unlike the basic reader, we prefer not to return blanks,
            # because we will typically wind up with a dict full of None
            # values
            # also, if the line begins with a comment character we
            # shouldn't return it either
            while row == [] or row[0][0] == self.commentchar:
                row = self.reader.next()
            d = dict(zip(self.fieldnames, row))
            lf = len(self.fieldnames)
            lr = len(row)
            if lf < lr:
                d[self.restkey] = row[lf:]
            elif lf > lr:
                for key in self.fieldnames[lr:]:
                    d[key] = self.restval
            return d


    with open(filename, 'r') as fid:
        fcsvreader = FiducialReader(fid)
        labels = []
        nac = []
        for line in fcsvreader:
            labels.append(line['#label'])
            nac.append((float(line['x']), float(line['y']), float(line['z'])))
    return labels, nac


def createSphereExpression(coordinates, radius=5):
    """
    For left DLPFC with Tailarach coordinates (-42, 30, 24) and 5mm radius sphere:

    3dcalc -a TT_icbm452+tlrc.nii \
           -expr 'step(25-(x-42)*(x-42)-(y-(-30))*(y-(-30))-(z-24)*(z-24))' \
           -prefix left_dlPFC3+tlrc.nii

    ### Nota Bene: x- and y- coordinates change sign in expression

    Return: 'step(25-(x-42)*(x-42)-(y+30)*(y+30)-(z-24)*(z-24))'

    Test: (-)
    """
    coordinates = tuple(coordinates)
    expression = 'step(%d-' % radius**2
    axes = ('x', 'y', 'z')
    for index in range(len(axes)):
        if axes[index] == 'z':
            value = int(coordinates[index]) * (-1)
            nextChar = ')'
        else:
            # We are reversing the y-coordinate b/c AFNI's RAS space is NOT DICOM-compliant.
            # It is actually "small-centric", not "large-centric"
            # if axes[index] == 'y':
            #     value = int(coordinates[index]) * (-1)
            # else:
            #     value = int(coordinates[index])
            value = int(coordinates[index]) * (-1)
            nextChar = '-'
        if value >= 0:
            expression += '({0}+{1})*({0}+{1})'.format(axes[index], value)
        else:
            # abs() to avoid printing "--"
            expression += '({0}-{1})*({0}-{1})'.format(axes[index], abs(value))
        expression += nextChar
    return expression


def formatFMRI(dicomDirectory):
    """

    :param dicomDirectory:
    """
    import subprocess

    outputs = subprocess.check_output(['formatFMRI.sh', dicomDirectory], stderr=subprocess.STDOUT).split(" ")
    # outputs = cmd.stdout.read().split(" ")
    # errors = cmd.stderr.read().split(" ")
    sliceOrder = outputs.pop()
    repetitionTime = outputs.pop()
    numberOfFiles = outputs.pop()
    numberOfSlices = outputs.pop()
    modality = outputs.pop()

    return modality, numberOfSlices, numberOfFiles, repetitionTime, sliceOrder

def resampleImage(inputVolume, outputVolume, resolution=(1.0, 1.0, 1.0)):
    """
    Resample the image using Linear interpolation (and identity transform) to the given resolution and write out the file

    """
    import os

    import SimpleITK as sitk

    image = sitk.ReadImage(inputVolume)
    old_size = list(image.GetSize())
    old_res = list(image.GetSpacing())
    new_res = list(resolution)
    new_size = [ int((old_size[i] * old_res[i]) / new_res[i]) for i in range(len(new_res))]

    resampler = sitk.ResampleImageFilter()
    resampler.SetSize(tuple(new_size))
    resampler.SetTransform(sitk.Transform(3, sitk.sitkIdentity))
    resampler.SetInterpolator(sitk.sitkLinear)
    resampler.SetOutputOrigin(image.GetOrigin())
    resampler.SetOutputSpacing(tuple(new_res))
    resampler.SetOutputDirection(image.GetDirection())
    outImage = resampler.Execute(image)

    felements = outputVolume.split('.')
    if len(felements) == 1:
        outputVolume += '.nii.gz'
    outputVolume = os.path.join(os.getcwd(), outputVolume)
    sitk.WriteImage(outImage, outputVolume)
    return outputVolume
