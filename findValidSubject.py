#!/usr/bin/env python
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

import fnmatch
import os

def findMatches(baseDir = '/', pattern = '*', index = 5):
    sessions = set()
    for root, dirs, files in os.walk(baseDir):
        for directory in dirs:
            full = os.path.join(root, directory)
            if fnmatch.fnmatch(full, pattern):
                # print "Match: \t%s \n\t %s\n\n" %(full, pattern)
                parts = full.split(os.sep)
                # print parts[index]
                sessions.add(parts[index])
    return sessions

if __name__ == '__main__':
    mrx = findMatches('/paulsen/MRx/FMRI_HD_120', '*/FMRI_HD_120/*/*/ANONRAW/FMRI_RestingStateConnectivity', 5)
    fs = findMatches('/paulsen/Experiments/20120722_JOY_DWI',
                     '*/FMRI_HD_120/*/*/10_AUTO.NN3Tv20110419/*_FS/mri', 6)
    canidates = mrx.intersection(fs)
    post = findMatches('/paulsen/Experiments/20120801.SubjectOrganized_Results',
                       '*/FMRI_HD_120/*/*/ACCUMULATED_POSTERIORS', 6)
    usable = canidates.intersection(post)
    print "\n\n\nThis is the usable list: \n\n\n"
    print usable
