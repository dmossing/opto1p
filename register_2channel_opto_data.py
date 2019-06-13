#!/usr/bin/env python

import sys
import os
sys.path.insert(0, '/home/mossing/code/adesnal')
import run_pipeline_tiffs as rpt

foldname = []
filenames = []

#foldname.append('190301/M9835/')
#filenames.append([1,888])
#
#foldname.append('190307/M9835/')
#filenames.append([1,999])
#
foldname.append('190318/M10338/')
filenames.append([2,999])

foldname.append('190320/M10365/')
filenames.append([1,999])

foldname.append('190410/M10368/')
filenames.append([2])

foldname.append('190411/M0002/')
filenames.append([4])

foldname.append('190501/M0094/')
filenames.append([4])

for i in range(len(foldname)):
    
    matlab_cmd = '"' + "gen_2channel_tiffs('" + foldname[i] + "'," + str(filenames[i]) + "); exit" + '"'
    print(matlab_cmd)
    os.system('matlab -r ' + matlab_cmd)

    fileparts = foldname[i].split('/') 
    date = fileparts[0]
    animalid = fileparts[1]
    expt_ids = [str(x) for x in filenames[i]]
    rpt.process_data(animalid,date,expt_ids,nchannels=2,delete_raw=True)
