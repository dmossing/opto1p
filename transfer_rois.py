#!/usr/bin/env python

import sys
import os
import numpy as np
#sys.path.insert(0, '/home/mossing/code/adesnal')
import run_pipeline_tiffs as rpt
import read_exptlist as re

matfile_fold = '/home/mossing/modulation/matfiles/'
suite2p_fold = '/home/mossing/data1/suite2P/results/'

def save_meanImg(datafold):
    vars_of_interest = ['meanImg','meanImg_chan2','meanImg_chan2_corrected','meanImgE']
    nplanes = 4
    planefolds = [datafold + '/suite2p/plane%d/' % iplane for iplane in range(nplanes)] 
    for fold in planefolds:
        ops = np.load(fold+'ops.npy',allow_pickle=True)[()]
        for var in vars_of_interest:
            if var in ops:
                np.save(fold+var+'.npy',ops[var])

def run(exptfilename,fileline=(1,2)):
    
    foldname = []
    filenames = []
    foldname,filenames = re.read_exptlist(exptfilename,lines_per_expt=4,fileline=fileline)
    print(fileline)
    print(filenames)
    print(foldname)
    
#    for i in range(len(foldname)):
#        # do matlab stuff to save cropping rectangles
#        print('now saving a cropping rectangle for ' + foldname[i])
    
    for i in range(len(foldname)):

        fileparts = foldname[i].split('/') 
        date = fileparts[0]
        animalid = fileparts[1]
        expt_ids = [str(x) for x in filenames[i]]
        subfold = '_'.join(expt_ids)

        thisfold = suite2p_fold + animalid + '/' + date + '/' + subfold + '/'

        save_meanImg(thisfold)

        matlab_cmd = '"' + "s2p_output_to_opto_corrected_rois('" + thisfold + "','datafold','" + matfile_fold + "'); exit;" + '"'

        print(matlab_cmd)
        os.system('matlab -r ' + matlab_cmd)


if __name__ == "__main__":
    if len(sys.argv)>2:
        fileline = (int(sys.argv[2]),)
    else:
        fileline = (1,2)
    run(sys.argv[1],fileline=fileline)
