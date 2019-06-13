#!/usr/bin/env python

import sys
import os
sys.path.insert(0, '/home/mossing/code/adesnal')
import run_pipeline_tiffs as rpt

suite2p_base = '/home/mossing/data_ssd/suite2P/'
result_base = suite2p_base + 'results/'
raw_base = suite2p_base + 'raw/'

foldname = []
filenames = []

foldname.append('190521/M0092/')
filenames.append([1,2,3,4])

foldname.append('190525/M0092/')
filenames.append([1,2,3,4])

foldname.append('190530/M0092/')
filenames.append([1,2,3])

for i in range(len(foldname)):
    
    matlab_cmd = '"' + "gen_2channel_tiffs('" + foldname[i] + "'," + str(filenames[i]) + ",1); exit" + '"'
    print(matlab_cmd)
    os.system('matlab -r ' + matlab_cmd)

    fileparts = foldname[i].split('/') 
    date = fileparts[0]
    animalid = fileparts[1]
    expt_ids = [str(x) for x in filenames[i]]
    rpt.process_data(animalid,date,expt_ids,nchannels=1,delete_raw=True,result_base=result_base,raw_base=raw_base)
