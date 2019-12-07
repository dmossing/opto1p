#!/usr/bin/env python

import sys
import os
sys.path.insert(0, '/home/mossing/code/adesnal')
import run_pipeline_tiffs as rpt
import read_exptlist as re

suite2p_fold = '/home/mossing/data1/suite2P/'
sbx_fold = '/home/mossing/modulation/2P/'
fast_disk = '/home/mossing/data_ssd/suite2P/bin'
delete_raw = True

raw_fold = suite2p_fold + 'raw/'
result_fold = suite2p_fold + 'results/'

def run(exptfilename,diameter=15):
    
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

        matlab_options = "options.green_only = 0; options.targetfold = '" + raw_fold + "'; options.data_foldbase = '" + sbx_fold + "'; "
        matlab_cmd = '"' + matlab_options + "gen_2channel_tiffs('" + foldname[i] + "'," + str(filenames_2ch[i]) + ",options); exit" + '"'
        print(matlab_cmd)
        os.system('matlab -r ' + matlab_cmd)
    
        rpt.add_2ch_data(animalid,date,expt_ids_1ch,expt_ids_2ch,raw_base=raw_fold,result_base=result_fold,delete_raw=delete_raw,diameter=diameter,fast_disk=fast_disk)

        try:
            shutil.rmtree(fast_disk+'/suite2p')
            print('fast disk contents deleted')
        except:
            print('fast disk location empty')

if __name__ == "__main__":
    if len(sys.argv)>2:
        run(sys.argv[1],diameter=sys.argv[2])
    else:
        run(sys.argv[1])
