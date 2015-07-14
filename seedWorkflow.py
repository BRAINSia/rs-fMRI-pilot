from csv import DictReader

from nipype.interfaces.utility import Function, IdentityInterface, Select
import nipype.pipeline.engine as pipe

from afni.preprocess import Calc


class FiducialReader(DictReader):
    """ Assumes file with text header """
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
    def fieldnames(self, values):
        self._fieldnames = values

    def next(self):
        if self.line_num == 0:
            # Used only for its side effect.
            self.fieldnames
        row = self.reader.next()
        self.line_num = self.reader.line_num
        # Unlike the basic reader, we prefer not to return blanks because we will
        # typically wind up with a dict full of None values.
        # If the line begins with a comment character we shouldn't return it either.
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


def getAtlasPoints(filename):
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

    >>> createSphereExpression((-42, 30, 24), radius=5)
    step(25-(x-42)*(x-42)-(y+30)*(y+30)-(z-24)*(z-24))

    """
    coordinates = tuple(coordinates)
    expression = 'step(%d-' % radius ** 2
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
            expression += '({0}-{1})*({0}-{1})'.format(axes[index], abs(value))  # abs() avoids printing "--"
        expression += nextChar
    return expression


def workflow(filename, outputType, name):
    seedworkflow = pipe.Workflow(name=name)
    labels, seeds = getAtlasPoints(filename)  # Create seed points

    seedsIdentity = pipe.Node(interface=IdentityInterface(fields=['index']), name='seedsIdentity')
    seedsIdentity.iterables = ('index', range(len(labels)))

    selectSeed = pipe.Node(interface=Select(), name='selectSeed')
    selectSeed.inputs.inlist = seeds

    selectLabel = pipe.Node(interface=Select(), name='selectLabel')
    selectLabel.inputs.inlist = labels

    points = pipe.Node(interface=Function(function=createSphereExpression,
                                          input_names=['coordinates', 'radius'],
                                          output_names=['expression']),
                       name='createSphereExpression')

    spheres = pipe.Node(interface=Calc(letters=['a']), name='afni3Dcalc_seeds')
    spheres.inputs.outputtype = outputType
    spheres.inputs.args = '-nscale'

    seedworkflow.connect([(seedsIdentity, selectSeed, [('index', 'index')]),
                      (seedsIdentity, selectLabel, [('index', 'index')]),
                      (selectSeed, points, [('out', 'coordinates')]),
                      (points, spheres, [('expression', 'expr')])])
    return seedworkflow
