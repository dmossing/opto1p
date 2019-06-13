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
# -

with h5py.File('/media/mossing/backup_1/data/2P/190227/M9835/ot/M9835_250_004_ot_000.rois',mode='r') as matfile:
    redratio = matfile['redratio'][:][:,0]
    msk = matfile['msk'][:].transpose((0,2,1))
    ctr = matfile['ctr'][:]

arglist = ['runtrial','trialrun','angle','size','contrast','light','trialwise','strialwise','dtrialwise','trialwise_t_offset']
axlist = [0,0,0,0,0,0,1,1,1,1]

data_struct = sio.loadmat('vip_opto_data_struct.mat',squeeze_me=True)

data_struct.keys()
keylist = [key for key in data_struct.keys() if not key[0]=='_']

data_struct[keylist[0]][()].dtype

iexpt = 0
key = keylist[iexpt]
data_struct[key][()]['stimulus_id']

pval_size_contrast = sio.loadmat('pval_size_contrast.mat',squeeze_me=True)

tuning = [None]*len(keylist)
strialavg = [None]*len(keylist)
ontarget_ret_lax = [None]*len(keylist)
bottom_3 = [None]*len(keylist)
run_cutoff = 7
for iexpt in range(len(keylist)):
    key = keylist[iexpt]
    data = data_struct[key][()]['decon']
    nbefore = data_struct[key][()]['nbefore']
    nafter = data_struct[key][()]['nafter']
    usize = data_struct[key][()]['stimulus_size_deg']
    ucontrast = data_struct[key][()]['stimulus_contrast']
    uangle = data_struct[key][()]['stimulus_direction']
    ulight = data_struct[key][()]['stimulus_light']
    size = usize[data_struct[key][()]['stimulus_id'][0]]
    contrast = ucontrast[data_struct[key][()]['stimulus_id'][1]]
    angle = uangle[data_struct[key][()]['stimulus_id'][2]]
    light = ulight[data_struct[key][()]['stimulus_id'][3]]
    trialrun = data_struct[key][()]['running_speed_cm_s']<run_cutoff
    ontarget_ret_lax[iexpt] = np.logical_and(data_struct[key][()]['rf_distance_deg']<5,data_struct[key][()]['rf_mapping_pval']<0.05)
    tuning[iexpt] = np.zeros((data.shape[0],len(usize),len(ucontrast),2)) #,data.shape[-1]))
    strialavg[iexpt] = np.zeros((data.shape[0],len(usize),len(ucontrast),len(uangle),2))#,data.shape[-1]))
    for k in range(2):
        for i in range(len(usize)):
            for j in range(len(ucontrast)):
                lkat = np.logical_and(np.logical_and(np.logical_and(size==usize[i],light==k),contrast==ucontrast[j]),np.logical_not(trialrun))
                tuning[iexpt][:,i,j,k] = data[:,lkat].mean(1) # for expts before and on 2/28, sign is reversed
                for w in range(len(uangle)):
                    strialavg[iexpt][:,i,j,w,k] = data[:,np.logical_and(lkat,angle==uangle[w])].mean(1) # for expts before and on 2/28, sign is reversed
            # convention here will be that 0 index is light off, 1 index is light on
    tuning[iexpt] = tuning[iexpt]/tuning[iexpt][:,:,:,:].max(1).max(1).max(1)[:,np.newaxis,np.newaxis,np.newaxis] #,np.newaxis] # ,nbefore:-nafter]
#     bottom_3[iexpt] = np.zeros((tuning[iexpt].shape[0],),dtype='bool')
#     bottom_3[iexpt][:int(np.cumsum(nbydepth[key])[-2])] = True

ops = np.load('/media/mossing/backup_1/data/suite2P/results/M9835/190227/3_4/suite2p/ops1.npy')

mred = ops[0]['meanImg_chan2_corrected']
mgreen = ops[0]['meanImg_chan2']

plt.figure()
plt.imshow(ops[0]['sdmov'])

plt.figure()
plt.imshow(msk.sum(0)*mred)
plt.figure()
plt.imshow(mred)
plt.figure()
plt.imshow(mgreen)#*(msk.sum(0)>0))


def combineRG(r,g,sat_by=(2,1)):
    R = sat_by[0]*r/r.max()
    G = sat_by[1]*g/g.max()
    B = np.zeros_like(r)
    return np.concatenate((R[:,:,np.newaxis],G[:,:,np.newaxis],B[:,:,np.newaxis]),axis=2)


red_green = combineRG(mred,mgreen,sat_by=(2,0.75))

np.where(redratio>0.65)[0]

plt.figure()
# for i in range(ctr[redratio>0.65].shape[0]):
#     plt.text(ctr[redratio>0.65][i,1],ctr[redratio>0.65][i,0],str(i))
ind = 56
img = sks.mark_boundaries(red_green,msk[ind].astype('int'),color=(1,1,1))
plt.imshow(img)
# plt.text(ctr[ind,1],ctr[ind,0],str(ind))
plt.axis('off')

np.where(redratio>0.65)[0][9]

tuning[0].shape

plt.figure(figsize=(5,15))
red_candidates = red_candidates-1 # convert from matlab indexing
for i in range(tuning.shape[0]): #(red_candidates.size):
    plt.subplot(38,10,i+1)
    ind = i# red_candidates[i]
    plt.plot(tuning[ind,:,0].T,c='k',alpha=0.2)
    plt.plot(tuning[ind,:,1].T,c='r',alpha=0.2)
    plt.axis('off')

np.where(redratio>0.65)[0]

plt.figure(figsize=(8,2))
ind = 10
for i in range(8):
    plt.subplot(1,8,i+1)
    plt.plot(tuning[ind,i,0],c='k',alpha=0.5)
    plt.plot(tuning[ind,i,1],c='r',alpha=1)
    plt.ylim((-1,5))
    plt.axis('off')
plt.figure()
img = sks.mark_boundaries(red_green,skm.binary_dilation(msk[ind]).astype('int'),color=(1,1,1))
plt.imshow(img)

redratio[69]

plt.figure(figsize=(8,2))
# actually vip negative
ind = 56
for i in range(8):
    plt.subplot(1,8,i+1)
    plt.plot(tuning[ind,i,0],c='k',alpha=0.5)
    plt.plot(tuning[ind,i,1],c='r',alpha=1)
    plt.ylim((-1,5))
    plt.axis('off')
# plt.savefig('representative_putative_vip_cell_56.png')
# plt.savefig('representative_putative_vip_cell_56.pdf')
plt.figure()
img = sks.mark_boundaries(red_green,skm.binary_dilation(msk[ind]).astype('int'),color=(1,1,1))
plt.imshow(img)
plt.axis('off')

plt.figure(figsize=(8,2))
ind = 51
for i in range(8):
    plt.subplot(1,8,i+1)
    plt.plot(tuning[ind,i,0],c='k',alpha=0.5)
    plt.plot(tuning[ind,i,1],c='r',alpha=1)
    plt.ylim((-1,5))
    plt.axis('off')
plt.figure()
img = sks.mark_boundaries(red_green,skm.binary_dilation(msk[ind]).astype('int'),color=(1,1,1))
plt.imshow(img)
plt.axis('off')

plt.figure(figsize=(8,2))
ind = 69
for i in range(8):
    plt.subplot(1,8,i+1)
    plt.plot(tuning[ind,i,0],c='k',alpha=0.5)
    plt.plot(tuning[ind,i,1],c='r',alpha=1)
    plt.ylim((-1,3))
    plt.axis('off')
plt.figure()
img = sks.mark_boundaries(red_green,skm.binary_dilation(msk[ind]).astype('int'),color=(1,1,1))
plt.imshow(img)
plt.axis('off')

plt.figure(figsize=(8,2))
ind = 133
for i in range(8):
    plt.subplot(1,8,i+1)
    plt.plot(tuning[ind,i,0],c='k',alpha=0.5)
    plt.plot(tuning[ind,i,1],c='r',alpha=1)
    plt.ylim((-1,3))
    plt.axis('off')
plt.figure()
img = sks.mark_boundaries(red_green,skm.binary_dilation(msk[ind]).astype('int'),color=(1,1,1))
plt.imshow(img)
plt.axis('off')

plt.figure(figsize=(8,2))
ind = 79
for i in range(8):
    plt.subplot(1,8,i+1)
    plt.plot(tuning[ind,i,0],c='k',alpha=0.5)
    plt.plot(tuning[ind,i,1],c='r',alpha=1)
    plt.ylim((-1,3))
    plt.axis('off')
plt.figure()
img = sks.mark_boundaries(red_green,skm.binary_dilation(msk[ind]).astype('int'),color=(1,1,1))
plt.imshow(img)
plt.axis('off')

plt.figure(figsize=(8,2))
ind = 91
for i in range(8):
    plt.subplot(1,8,i+1)
    plt.plot(tuning[ind,i,0],c='k',alpha=0.5)
    plt.plot(tuning[ind,i,1],c='r',alpha=1)
    plt.ylim((-1,3))
    plt.axis('off')
plt.figure()
img = sks.mark_boundaries(red_green,skm.binary_dilation(msk[ind]).astype('int'),color=(1,1,1))
plt.imshow(img)
plt.axis('off')

plt.figure(figsize=(8,2))
ind = 347
for i in range(8):
    plt.subplot(1,8,i+1)
    plt.plot(tuning[ind,i,0],c='k',alpha=0.5)
    plt.plot(tuning[ind,i,1],c='r',alpha=1)
    plt.ylim((-1,3))
    plt.axis('off')
plt.figure()
img = sks.mark_boundaries(red_green,skm.binary_dilation(msk[ind]).astype('int'),color=(1,1,1))
plt.imshow(img)
plt.axis('off')

inds = (10,51,56,69,79,91,133,347,370)
for ind in inds:
    plt.figure(figsize=(8,2))
    for i in range(8):
        plt.subplot(1,8,i+1)
        plt.plot(tuning[ind,i,0],c='k',alpha=0.5)
        plt.plot(tuning[ind,i,1],c='r',alpha=1)
        plt.ylim((-1,3))
        plt.axis('off')
        
    plt.figure()
    img = sks.mark_boundaries(red_green,skm.binary_dilation(msk[ind]).astype('int'),color=(1,1,1))
    c = ctr[ind].astype('int')
    plt.imshow(img[c[0]-25:c[0]+25,c[1]-25:c[1]+25])
    plt.axis('off')

np.where(redratio>0.65)[0]

sdmov = ops[0]['sdmov']

plt.figure()
plt.imshow(sdmov)

ops[0].keys()

ops[0]['maxregshift']

offsetx = 12
offsety = 15
slicey = slice(offsety,offsety+sdmov.shape[0])
slicex = slice(offsetx,offsetx+sdmov.shape[1])
mean_sd = combineRG(mgreen[slicey,slicex],sdmov,sat_by=(1,1))
red_green_sd = combineRG(mred[slicey,slicex],sdmov-sdmov.min(),sat_by=(2,3))

plt.figure()
plt.imshow(red_green_sd)

(redratio>0.5).sum()

plt.figure(figsize=(15,7.5))
for i,ind in enumerate(np.where(redratio>0.5)[0]):
    plt.subplot(5,10,i+1)
#     img = sks.mark_boundaries(red_green,skm.binary_dilation(msk[ind]).astype('int'),color=(1,1,1))
    img = sks.mark_boundaries(red_green_sd,skm.binary_dilation(msk[ind,slicey,slicex]).astype('int'),color=(1,1,1))
#     c = ctr[ind].astype('int')
    c = (ctr[ind]-np.array((offsety,offsetx))).astype('int')
    plt.imshow(img[c[0]-25:c[0]+25,c[1]-25:c[1]+25])
    plt.axis('off')
