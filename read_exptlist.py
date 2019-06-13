# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.4'
#       jupytext_version: 1.1.6
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

import numpy as np


def read_exptlist(filename,lines_per_expt=3):
    with open('exptlist.txt',mode='r') as f:
        content = f.read().splitlines()
    foldname = []
    filenames = []
    for iline,line in enumerate(content):
        if np.mod(iline,lines_per_expt)==0:
            foldname.append(line)
        elif np.mod(iline,lines_per_expt)==1:
            filenames.append([int(x) for x in line.split(',')])
    return foldname,filenames


