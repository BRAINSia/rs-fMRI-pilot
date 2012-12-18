#!/usr/bin/env python
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

import fnmatch
import glob
import os


def findMatches(filename, pattern = '*', index = 5):
    with open(filename, 'w+') as fid:
        for path in glob.iglob(pattern):
            parts = path.split(os.sep)
            fid.write('{0}\n'.format(parts[index]))
    with open(filename, 'r') as fid2:
        raw_sessions = fid2.readlines()
    sessions = []
    for item in raw_sessions:
        sessions.append(str(item.rstrip('\n')))
    return set(sessions)


if __name__ == '__main__':
    mrx = findMatches('dicom.txt', '/paulsen/MRx/FMRI_HD_120/*/*/ANONRAW/FMRI_RestingStateConnectivity', 5)
    print "Found all DICOMs"
    print mrx
    fs = findMatches('freesurf.txt', '/paulsen/Experiments/20120722_JOY_DWI/FMRI_HD_120/*/*/10_AUTO.NN3Tv20110419/*_FS/mri', 6)
    print "Found all brain.mgz"
    print fs
    canidates = mrx.intersection(fs)
    post = findMatches('posterior.txt',
                       '/paulsen/Experiments/20120801.SubjectOrganized_Results/FMRI_HD_120/*/*/ACCUMULATED_POSTERIORS', 6)
    print "Found all posteriors"
    print post
    usable = canidates.intersection(post)
    print "\n\n\nThis is the usable list: \n\n\n"
    print usable
