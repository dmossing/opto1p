#!/usr/bin/env python

import sys
import os

with open('path_locations.txt','r') as f:
    path_locs = f.read().splitlines()
for loc in path_locs:    
    print(loc)
    sys.path.insert(0,loc)
import run_pipeline_tiffs as rpt
import read_exptlist as re

delete_raw = True

def run(exptfilename,diameter=15,sbx_fold='/home/mossing/modulation/2P/',suite2p_fold='/home/mossing/data1/suite2P/',fast_disk='/home/mossing/data_ssd/suite2P/bin',matfile_fold='/home/mossing/modulation/matfiles/'):

    raw_fold = suite2p_fold + 'raw/'
    result_fold = suite2p_fold + 'results/'

    foldname = []
    filenames = []
    foldname,filenames_1ch = re.read_exptlist(exptfilename,lines_per_expt=4,fileline=(1,2))  # should this be (1,2) ? 
    _,filenames_2ch = re.read_exptlist(exptfilename,lines_per_expt=4,fileline=2)
    
    for i in range(len(foldname)):

        fileparts = foldname[i].split('/') 
        date = fileparts[0]
        animalid = fileparts[1]
        expt_ids_1ch = [str(x) for x in filenames_1ch[i]]
        expt_ids_2ch = [str(x) for x in filenames_2ch[i]]

        try:
            shutil.rmtree(fast_disk+'/suite2p')
            os.system('rm -rf '+fast_disk+'/suite2p')
            print('fast disk contents deleted')
        except:
            print('fast disk location empty')
        
        matlab_options = "options.green_only = 1; options.targetfold = '" + raw_fold + "'; options.data_foldbase = '" + sbx_fold + "'; options.matfile_foldbase = '" + matfile_fold + "';"
        matlab_cmd = '"' + matlab_options + "gen_2channel_tiffs('" + foldname[i] + "'," + str(filenames_1ch[i]) + ",options); exit" + '"'
        print(matlab_cmd)
        os.system('matlab -r ' + matlab_cmd)

        try:
            shutil.rmtree(fast_disk+'/suite2p')
            os.system('rm -rf '+fast_disk+'/suite2p')
            print('fast disk contents deleted')
        except:
            print('fast disk location empty')

        rpt.process_data(animalid,date,expt_ids_1ch,delete_raw=delete_raw,raw_base=raw_fold,result_base=result_fold,diameter=diameter,nchannels=1,fast_disk=fast_disk)

        matlab_options = "addpath(genpath('~/adesnal')); options.green_only = 0; options.targetfold = '" + raw_fold + "'; options.data_foldbase = '" + sbx_fold + "'; options.matfile_foldbase = '" + matfile_fold + "';"
        matlab_cmd = '"' + matlab_options + "gen_2channel_tiffs('" + foldname[i] + "'," + str(filenames_2ch[i]) + ",options); exit" + '"'
        print(matlab_cmd)
        os.system('matlab -r ' + matlab_cmd)
    
        rpt.add_2ch_data(animalid,date,expt_ids_1ch,expt_ids_2ch,raw_base=raw_fold,result_base=result_fold,delete_raw=delete_raw,diameter=diameter,fast_disk=fast_disk)

        try:
            shutil.rmtree(fast_disk+'/suite2p')
            os.system('rm -rf '+fast_disk+'/suite2p')
            print('fast disk contents deleted')
        except:
            print('fast disk location empty')

if __name__ == "__main__":
    if len(sys.argv)==5:
        location = sys.argv[4]
        if location == 'big-boi':
            suite2p_fold  = '/home/mossing/data1/suite2P/'
            fast_disk = '/home/mossing/data_ssd/suite2P/bin/'
            matfile_fold = '/home/mossing/modulation/matfiles/'
            sbx_fold = sys.argv[3]
        elif (location == 'cluster') or (location == 'savio'):
            suite2p_fold  = '/global/scratch/mossing/2Pdata/suite2P/'
            fast_disk = '/global/scratch/mossing/2Pdata/suite2P/bin/'
            sbx_fold = '/global/scratch/mossing/2Pdata/'
            matfile_fold = '/global/scratch/mossing/matfiles/'
        run(sys.argv[1],diameter=sys.argv[2],sbx_fold=sbx_fold,suite2p_fold=suite2p_fold,fast_disk=fast_disk,matfile_fold=matfile_fold)
    elif len(sys.argv)==4:
        run(sys.argv[1],diameter=sys.argv[2],sbx_fold=sys.argv[3])
    elif len(sys.argv)==3:
        run(sys.argv[1],diameter=sys.argv[2])
    else:
        run(sys.argv[1])
