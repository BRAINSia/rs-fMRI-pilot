#!/usr/bin/env python
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

import fnmatch
import os

def findMatches(baseDir = '/paulsen/MRx', pattern = 'FMRI_HD_120/*/*'):
    sessions = set()
    for root, dirs, files in os.walk(baseDir):
        for directory in dirs:
            if fnmatch.fnmatch(directory, pattern):
                parts = directory.split(os.sep)
                sessions.add(parts[4])
    return sessions

if __name__ == '__main__':
    mrx = findMatches()
    fs = findMatches('/paulsen/Experiments/20120722_JOY_DWI',
                     'FMRI_HD_120/*/*/10_AUTO.NN3Tv20110419/*_FS/mri')
    canidates = mrx.intersection(fs)
    post = findMatches('/paulsen/Experiments/20120801.SubjectOrganized_Results',
                       'FMRI_HD_120/*/*/ACCUMULATED_POSTERIORS')
    usable = canidates.intersection(post)
    print usable
