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

# +
import os
import scipy.io as sio
import matplotlib.pyplot as plt
# %matplotlib notebook
import numpy as np
import h5py
from oasis.functions import deconvolve
from oasis import oasisAR1, oasisAR2
import pyute as ut

from importlib import reload
reload(ut)
import scipy.ndimage.filters as sfi
import scipy.stats as sst
import scipy.ndimage.measurements as snm
from mpl_toolkits.mplot3d import Axes3D
import size_contrast_opto_analysis as scoa
reload(scoa)
import retinotopy_analysis as rt
reload(rt)
import naka_rushton_analysis as nra
import pdb
import size_contrast_analysis as sca
import skimage.segmentation as sks
import skimage.morphology as skm


# +
def tack_on(thisfold,thisfile,rg=(2,-10),criterion=lambda x:np.abs(x)>100,datafoldbase=None,stimfoldbase=None):
    folds.append(thisfold)
    files.append(thisfile)
    rgs.append(rg)
    criteria.append(criterion)
    datafoldbases.append(datafoldbase)
    stimfoldbases.append(stimfoldbase)
    
folds = []
files = []
rgs = []
criteria = []
datafoldbases = []
stimfoldbases = []

# thisfold = '181030/M9826/'
# thisfile = 'M9826_050_'
# retnumber = '001'
# thisrg = (1,-10)
# tack_on(thisfold,thisfile+retnumber,thisrg,criterion=lambda x:np.abs(x)<100)

# thisfold = '181109/M9826/'
# thisfile = 'M9826_065_'
# retnumber = '004'
# thisrg = (1,-10)
# tack_on(thisfold,thisfile+retnumber,criterion=lambda x:np.abs(x)>100,rg=thisrg)

# thisfold = '181115/M9826/'
# thisfile = 'M9826_100_'
# retnumber = '003'
# thisrg = (1,-10)
# tack_on(thisfold,thisfile+retnumber,criterion=lambda x:np.abs(x)>100,rg=thisrg)

# thisfold = '190221/M9835/'
# thisfile = 'M9835_115_'
# retnumber = '002'
# thisrg = (1,-10)
# datafoldbase='/media/mossing/data_ssd/data/2P/green_only/'
# stimfoldbase='/home/mossing/modulation/visual_stim/'
# tack_on(thisfold,thisfile+retnumber,criterion=lambda x:np.abs(x)>100,rg=thisrg,datafoldbase=datafoldbase,stimfoldbase=stimfoldbase)

# thisfold = '190223/M9835/'
# thisfile = 'M9835_110_'
# retnumber = '002'
# thisrg = (1,-10)
# datafoldbase='/media/mossing/data_ssd/data/2P/'
# stimfoldbase='/home/mossing/modulation/visual_stim/'
# tack_on(thisfold,thisfile+retnumber,criterion=lambda x:np.abs(x)>0,rg=thisrg,datafoldbase=datafoldbase,stimfoldbase=stimfoldbase)

#thisfold = '190227/M9835/' # this one shows confirmation of VIP suppression in a few cells (~ all 3 that were clearly visible)
#thisfile = 'M9835_250_'
#retnumber = '003'
#thisrg = (1,-10)
#datafoldbase='/media/mossing/data_ssd/data/2P/'
#stimfoldbase='/home/mossing/modulation/visual_stim/'
#tack_on(thisfold,thisfile+retnumber,criterion=lambda x:np.abs(x)>0,rg=thisrg,datafoldbase=datafoldbase,stimfoldbase=stimfoldbase)

thisfold = '190301/M9835/'
thisfile = 'M9835_175_'
retnumber = '001'
thisrg = (0,-10)
datafoldbase='/media/mossing/backup_1/data/2P/'
stimfoldbase='/home/mossing/modulation/visual_stim/'
tack_on(thisfold,thisfile+retnumber,criterion=lambda x:np.abs(x)>0,rg=thisrg,datafoldbase=datafoldbase,stimfoldbase=stimfoldbase)

thisfold = '190307/M9835/'
thisfile = 'M9835_095_'
retnumber = '001'
thisrg = (1,-11)
datafoldbase='/media/mossing/backup_1/data/2P/'
stimfoldbase='/home/mossing/modulation/visual_stim/'
tack_on(thisfold,thisfile+retnumber,criterion=lambda x:np.abs(x)>0,rg=thisrg,datafoldbase=datafoldbase,stimfoldbase=stimfoldbase)

thisfold = '190314/M10388/'
thisfile = 'M10388_140_'
retnumber = '001'
thisrg = (1,-11)
datafoldbase='/media/mossing/backup_1/data/2P/'
stimfoldbase='/home/mossing/modulation/visual_stim/'
tack_on(thisfold,thisfile+retnumber,criterion=lambda x:np.abs(x)>0,rg=thisrg,datafoldbase=datafoldbase,stimfoldbase=stimfoldbase)

thisfold = '190318/M10338/'
thisfile = 'M10338_180_'
retnumber = '002'
thisrg = (1,-10)
datafoldbase='/media/mossing/backup_0/data/2P/'
stimfoldbase='/home/mossing/modulation/visual_stim/'
tack_on(thisfold,thisfile+retnumber,criterion=lambda x:np.abs(x)>0,rg=thisrg,datafoldbase=datafoldbase,stimfoldbase=stimfoldbase)

thisfold = '190320/M10365/'
thisfile = 'M10365_125_'
retnumber = '001'
thisrg = (0,-11)
datafoldbase='/media/mossing/backup_0/data/2P/'
stimfoldbase='/home/mossing/modulation/visual_stim/'
tack_on(thisfold,thisfile+retnumber,criterion=lambda x:np.abs(x)>0,rg=thisrg,datafoldbase=datafoldbase,stimfoldbase=stimfoldbase)

thisfold = '190410/M10368/'
thisfile = 'M10368_080_'
retnumber = '002'
thisrg = (1,-10)
datafoldbase='/media/mossing/backup_1/data/2P/'
stimfoldbase='/home/mossing/modulation/visual_stim/'
tack_on(thisfold,thisfile+retnumber,criterion=lambda x:np.abs(x)>0,rg=thisrg,datafoldbase=datafoldbase,stimfoldbase=stimfoldbase)

thisfold = '190411/M0002/'
thisfile = 'M0002_100_'
retnumber = '004'
thisrg = (1,-10)
datafoldbase='/media/mossing/backup_0/data/2P/'
stimfoldbase='/home/mossing/modulation/visual_stim/'
tack_on(thisfold,thisfile+retnumber,criterion=lambda x:np.abs(x)>0,rg=thisrg,datafoldbase=datafoldbase,stimfoldbase=stimfoldbase)

thisfold = '190501/M0094/'
thisfile = 'M0094_130_'
retnumber = '002'
thisrg = (0,-11)
datafoldbase='/media/mossing/backup_0/data/2P/'
stimfoldbase='/home/mossing/modulation/visual_stim/'
tack_on(thisfold,thisfile+retnumber,criterion=lambda x:np.abs(x)>0,rg=thisrg,datafoldbase=datafoldbase,stimfoldbase=stimfoldbase)

reload(rt)
reload(ut)
ret,paramdict,pval,trialrun,has_inverse,nbydepth,proc = rt.analyze_everything(folds,files,rgs,criteria,datafoldbase=datafoldbases,stimfoldbase=stimfoldbases)
# -

keylist = list(ret.keys())

scoa.add_data_struct_h5('vip_halo_data_struct.hdf5',cell_type='synapsin',keylist=keylist,proc=proc)

k = 0
# sio.savemat('vip_roi_ctrs.mat',{'ctrx':paramdict[keylist[k]]['xo'][()],'ctry':paramdict[keylist[k]]['yo'][()]})

for k in range(len(keylist)):
    proc[keylist[k]]['paramdict'] = paramdict[keylist[k]]
    proc[keylist[k]]['trialrun'] = trialrun[keylist[k]]
    proc[keylist[k]]['pval'] = pval[keylist[k]]
    proc[keylist[k]]['ret'] = ret[keylist[k]]
    proc[keylist[k]]['nbydepth'] = nbydepth[keylist[k]]

# +
# data_struct = rt.gen_full_data_struct(keylist=keylist,proc=proc,nbefore=8,nafter=8)

# +
# sio.savemat('vip_retinotopy_data_struct.mat',data_struct)
# -

k = 0
plt.figure()
plt.plot(np.nanmean(np.nanmean(proc[keylist[k]]['strialwise'],0),0))
for i in range(int(np.floor(proc[keylist[k]]['strialwise'].shape[0]/100))):
    plt.plot(np.nanmean(proc[keylist[k]]['strialwise'][i*100:(i+1)*100],1).T,alpha=1e-2,c='b')
plt.ylim((0,0.1))

plt.figure()
k = 0
ut.imshow_in_rows(ret[keylist[k]][:100])

# +
folds = []
files = []
rets = []
adjust_fns = []
rgs = []
criteria = []
datafoldbases = []
stimfoldbases = []

def tack_on(thisfold,thisfile,retnumber,frame_adjust=None,rg=(1,0),criterion=lambda x: np.abs(x)>100,datafoldbase=None,stimfoldbase=None):
    folds.append(thisfold)
    files.append(thisfile)
    rets.append(retnumber)
    adjust_fns.append(frame_adjust)
    rgs.append(rg)
    criteria.append(criterion)
    datafoldbases.append(datafoldbase)
    stimfoldbases.append(stimfoldbase)

# thisfold = '181030/M9826/'
# thisfile = 'M9826_050_002'
# retnumber = '001'
# criterion = lambda x: np.abs(x)<100
# frame_adjust = lambda x: np.hstack((x[0::4][:,np.newaxis],x[3::4][:,np.newaxis])).flatten()
# rg = (1,0)
# tack_on(thisfold,thisfile,retnumber,rg=rg,criterion=criterion,frame_adjust=frame_adjust)

# thisfold = '181109/M9826/'
# thisfile = 'M9826_065_005'
# retnumber = '004'
# criterion = lambda x: np.abs(x)>100
# frame_adjust = lambda x: np.hstack((x[0::4][:,np.newaxis],x[3::4][:,np.newaxis])).flatten()
# rg = (1,0)
# tack_on(thisfold,thisfile,retnumber,rg=rg,criterion=criterion,frame_adjust=frame_adjust)

# thisfold = '181115/M9826/'
# thisfile = 'M9826_100_004'
# retnumber = '003'
# criterion = lambda x: np.abs(x)>100
# frame_adjust = lambda x: np.hstack((x[0::4][:,np.newaxis],x[3::4][:,np.newaxis])).flatten()
# rg = (1,0)
# tack_on(thisfold,thisfile,retnumber,rg=rg,criterion=criterion,frame_adjust=frame_adjust)

# thisfold = '190223/M9835/'
# thisfile = 'M9835_110_003'
# thisrg = (1,-1)
# retnumber = '002'
# datafoldbase='/media/mossing/data_ssd/data/2P/'
# stimfoldbase='/home/mossing/modulation/visual_stim/'
# frame_adjust = lambda x: np.hstack((x[0::4][:,np.newaxis],x[3::4][:,np.newaxis])).flatten()
# tack_on(thisfold,thisfile,retnumber,frame_adjust=frame_adjust,criterion=lambda x:np.abs(x)>100,rg=thisrg,datafoldbase=datafoldbase,stimfoldbase=stimfoldbase)

#thisfold = '190227/M9835/'
#thisfile = 'M9835_250_004'
#thisrg = (1,0)
#retnumber = '003'
#datafoldbase='/media/mossing/data_ssd/data/2P/'
#stimfoldbase='/home/mossing/modulation/visual_stim/'
#frame_adjust = lambda x: np.hstack((x[0::4][:,np.newaxis],x[3::4][:,np.newaxis])).flatten()
#tack_on(thisfold,thisfile,retnumber,frame_adjust=frame_adjust,criterion=lambda x:np.abs(x)>100,rg=thisrg,datafoldbase=datafoldbase,stimfoldbase=stimfoldbase)

thisfold = '190301/M9835/'
thisfile = 'M9835_175_002'
thisrg = (1,-1)
retnumber = '001'
datafoldbase='/media/mossing/backup_1/data/2P/'
stimfoldbase='/home/mossing/modulation/visual_stim/'
frame_adjust = lambda x: np.hstack((x[0::4][:,np.newaxis],x[3::4][:,np.newaxis])).flatten()
tack_on(thisfold,thisfile,retnumber,frame_adjust=frame_adjust,criterion=lambda x:np.abs(x)>100,rg=thisrg,datafoldbase=datafoldbase,stimfoldbase=stimfoldbase)

thisfold = '190318/M10338/'
thisfile = 'M10338_180_003'
thisrg = (1,0)
retnumber = '002'
datafoldbase='/media/mossing/backup_0/data/2P/'
stimfoldbase='/home/mossing/modulation/visual_stim/'
frame_adjust = lambda x: np.hstack((x[0::4][:,np.newaxis],x[3::4][:,np.newaxis])).flatten()
tack_on(thisfold,thisfile,retnumber,frame_adjust=frame_adjust,criterion=lambda x:np.abs(x)>100,rg=thisrg,datafoldbase=datafoldbase,stimfoldbase=stimfoldbase)

thisfold = '190320/M10365/'
thisfile = 'M10365_125_002'
thisrg = (1,0)
retnumber = '001'
datafoldbase='/media/mossing/backup_0/data/2P/'
stimfoldbase='/home/mossing/modulation/visual_stim/'
frame_adjust = lambda x: np.hstack((x[0::4][:,np.newaxis],x[3::4][:,np.newaxis])).flatten()
tack_on(thisfold,thisfile,retnumber,frame_adjust=frame_adjust,criterion=lambda x:np.abs(x)>100,rg=thisrg,datafoldbase=datafoldbase,stimfoldbase=stimfoldbase)

thisfold = '190410/M10368/'
thisfile = 'M10368_080_003'
thisrg = (1,-1)
retnumber = '002'
datafoldbase='/media/mossing/backup_1/data/2P/'
stimfoldbase='/home/mossing/modulation/visual_stim/'
frame_adjust = lambda x: np.hstack((x[0::4][:,np.newaxis],x[3::4][:,np.newaxis])).flatten()
tack_on(thisfold,thisfile,retnumber,frame_adjust=frame_adjust,criterion=lambda x:np.abs(x)>100,rg=thisrg,datafoldbase=datafoldbase,stimfoldbase=stimfoldbase)

thisfold = '190411/M0002/'
thisfile = 'M0002_100_005'
thisrg = (1,0)
retnumber = '004'
datafoldbase='/media/mossing/backup_0/data/2P/'
stimfoldbase='/home/mossing/modulation/visual_stim/'
frame_adjust = lambda x: np.hstack((x[0::4][:,np.newaxis],x[3::4][:,np.newaxis])).flatten()
tack_on(thisfold,thisfile,retnumber,frame_adjust=frame_adjust,criterion=lambda x:np.abs(x)>100,rg=thisrg,datafoldbase=datafoldbase,stimfoldbase=stimfoldbase)

thisfold = '190501/M0094/'
thisfile = 'M0094_130_003'
thisrg = (1,-1)
retnumber = '002'
datafoldbase='/media/mossing/backup_0/data/2P/'
stimfoldbase='/home/mossing/modulation/visual_stim/'
frame_adjust = lambda x: np.hstack((x[0::4][:,np.newaxis],x[3::4][:,np.newaxis])).flatten()
tack_on(thisfold,thisfile,retnumber,frame_adjust=frame_adjust,criterion=lambda x:np.abs(x)>100,rg=thisrg,datafoldbase=datafoldbase,stimfoldbase=stimfoldbase)


reload(ut)
reload(scoa)
soriavg,strialavg,lb,ub,pval,nbydepth,spont,ret_vars,Smean_stat = scoa.analyze_everything_by_criterion(folds,files,rets,adjust_fns,rgs,criteria=criteria,datafoldbase=datafoldbases,stimfoldbase=stimfoldbases,procname='size_contrast_opto_proc.hdf5')
# -

reload(scoa)
with h5py.File('size_contrast_opto_proc.hdf5',mode='r') as proc:
    keylist = list([x for x in proc.keys() if not x[0]=='_'])
    scoa.add_data_struct_h5('vip_halo_data_struct.hdf5',cell_type='synapsin', keylist=keylist, frame_rate_dict=None, proc=proc, nbefore=8, nafter=8)

pdb.pm()

keylist = list(proc.keys())

arglist = ['runtrial','trialrun','angle','size','contrast','light','trialwise','strialwise','dtrialwise','trialwise_t_offset']
axlist = [0,0,0,0,0,0,1,1,1,1]

proc[keylist[1]][1].keys()

for key in keylist:
    print(ret_vars[key].keys())

# +
ontarget_ret_lax = {}
keylist = list(proc.keys())
for key in keylist:
#     ontarget_ret_lax[key] = np.nanmin(pval[key],1)<0.05/np.logical_not(np.isnan(pval[key])).sum(1) #np.ones((soriavg[key].shape[0],),dtype='bool')
    ontarget_ret_lax[key],_ = rt.ontarget_by_retinotopy(ret_vars[key],ctr=ret_vars[key]['position'],rg=10)

offtarget_ret_lax = {}
keylist = list(proc.keys())
for key in keylist:
#     ontarget_ret_lax[key] = np.nanmin(pval[key],1)<0.05/np.logical_not(np.isnan(pval[key])).sum(1) #np.ones((soriavg[key].shape[0],),dtype='bool')
    offtarget_ret_lax[key],_ = rt.ontarget_by_retinotopy(ret_vars[key],ctr=ret_vars[key]['position'],rg=20)
    offtarget_ret_lax[key][ontarget_ret_lax[key]] = False
# -

for k in range(len(keylist)):
    plt.figure()
    plt.plot(proc[keylist[k]][1]['strialwise'].mean(0).mean(0))

reload(scoa)
data_struct = scoa.gen_full_data_struct(cell_type='PyrL23', keylist=keylist, frame_rate_dict=None, proc=proc, ret_vars=ret_vars, nbefore=4, nafter=4)

sio.savemat('vip_opto_data_struct.mat',data_struct)

data_struct.keys()

pval_size_contrast = {}
for k,key in enumerate(data_struct.keys()):
    pval_size_contrast[key] = pval[keylist[k]]

sio.savemat('pval_size_contrast.mat',pval_size_contrast)


