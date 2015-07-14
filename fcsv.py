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
