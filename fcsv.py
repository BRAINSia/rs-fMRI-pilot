from csv import reader, DictReader

class FiducialReader(DictReader):
    def __init__(self, filename, fieldnames=None, restkey=None, restval=None,
                 dialect='excel', commentchar='#', *args, **kwds):
        self._fieldnames = fieldnames   # list of keys for the dict
        self.restkey = restkey          # key to catch long rows
        self.restval = restval          # default value for short rows
        self.reader = reader(filename, dialect, *args, **kwds)
        self.dialect = dialect
        self.line_num = 0
        self.comment_line_num = -1
        super(DictReader, self).__init__(*args, **kwds)

    def __iter__(self):
        return self

    @property
    def fieldnames(self):
        if self._fieldnames is None:
            try:
                self._fieldnames = self.reader.next()
            except StopIteration:
                pass
        self.line_num = self.reader.line_num
        return self._fieldnames

    @fieldnames.setter
    def fieldnames(self, value):
        self._fieldnames = value

    def next(self):
        if self.line_num == 0:
           # Used only for its side effect.
           self.fieldnames
        row = self.reader.next()
        assert self.comment_line_num < self.line_num, 'First line MUST be a commented line!'
        temp = []
        if row[0][0] = self.commentchar:
            temp = row
            self.comment_line_num += 1

        else:


        self.line_num = self.reader.line_num

        # unlike the basic reader, we prefer not to return blanks,
        # because we will typically wind up with a dict full of None
        # values
        while row == []:
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
