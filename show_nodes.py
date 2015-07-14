#!/bin/env python

import cPickle
import gzip
from pprint import pprint
import sys

from nipype.interfaces.base import InterfaceResult, Interface

if __name__ == '__main__':
    infile = sys.argv[1]
    with gzip.open(infile) as fid:
        node = cPickle.load(fid)
    if isinstance(node, InterfaceResult):
        for key in ['inputs', 'interface', 'outputs', 'provenance', 'runtime', 'version']:
            print "---- %s ----" % key
            pprint(eval("node.%s" % key))
    elif isinstance(node, Interface):
        print "This is an interactive node - you should load this in iPython to explore"
    else:
        pprint(node)
    sys.exit(0)
