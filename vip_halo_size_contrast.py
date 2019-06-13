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
# with h5py.File('/media/mossing/data_ssd/data/2P/190227/M9835/ot/M9835_250_004_ot_000.rois',mode='r') as matfile:
#     redratio = matfile['redratio'][:][:,0]
#     msk = matfile['msk'][:].transpose((0,2,1))
#     ctr = matfile['ctr'][:]
# -

arglist = ['runtrial','trialrun','angle','size','contrast','light','trialwise','strialwise','dtrialwise','trialwise_t_offset']
axlist = [0,0,0,0,0,0,1,1,1,1]

proc[keylist[1]][0].keys()

# +
# plt.figure()
# plt.hist(ret_vars[keylist[0]]['paramdict_normal']['amplitude'][()][ret_vars[keylist[0]]['pval_ret']<0.05],bins=np.linspace(-1,1,100))

# +
#ontarget_ret_lax = {}
#keylist = list(proc.keys())
#for key in keylist:
##     ontarget_ret_lax[key] = np.nanmin(pval[key],1)<0.05/np.logical_not(np.isnan(pval[key])).sum(1) #np.ones((soriavg[key].shape[0],),dtype='bool')
#    ontarget_ret_lax[key],_ = rt.ontarget_by_retinotopy(ret_vars[key],ctr=ret_vars[key]['position'],rg=10)

#offtarget_ret_lax = {}
#keylist = list(proc.keys())
#for key in keylist:
##     ontarget_ret_lax[key] = np.nanmin(pval[key],1)<0.05/np.logical_not(np.isnan(pval[key])).sum(1) #np.ones((soriavg[key].shape[0],),dtype='bool')
#    offtarget_ret_lax[key],_ = rt.ontarget_by_retinotopy(ret_vars[key],ctr=ret_vars[key]['position'],rg=20)
#    offtarget_ret_lax[key][ontarget_ret_lax[key]] = False
# -

data_struct = sio.loadmat('vip_opto_data_struct.mat',squeeze_me=True)

# +
ontarget_ret_lax = {}
keylist = [x for x in list(data_struct.keys()) if not x[0]=='_']
for key in keylist:
#     ontarget_ret_lax[key] = np.nanmin(pval[key],1)<0.05/np.logical_not(np.isnan(pval[key])).sum(1) #np.ones((soriavg[key].shape[0],),dtype='bool')
    ontarget_ret_lax[key] = np.logical_and(data_struct[key][()]['rf_mapping_pval']<0.05,data_struct[key][()]['rf_distance_deg']<5)
    ontarget_ret_lax[key] = np.ones(ontarget_ret_lax[key].shape,dtype='bool')
    
offtarget_ret_lax = {}
for key in keylist:
#     ontarget_ret_lax[key] = np.nanmin(pval[key],1)<0.05/np.logical_not(np.isnan(pval[key])).sum(1) #np.ones((soriavg[key].shape[0],),dtype='bool')
    offtarget_ret_lax[key] = np.logical_and(data_struct[key][()]['rf_mapping_pval']<0.05,data_struct[key][()]['rf_distance_deg']<20)
    offtarget_ret_lax[key][ontarget_ret_lax[key]] = False
# -



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



keylist

for k in range(len(tuning)):
#     data = tuning[k][np.logical_and(bottom_3[k],ontarget_ret_lax[keylist[k]])][:,:,:,:,nbefore:-nafter].mean(0).mean(-1)
    data = tuning[k][ontarget_ret_lax[k]].mean(0) #[:,:,:,nbefore:-nafter].mean(-1)
    plt.figure()
    plt.subplot(1,2,1)
    plt.imshow(data[:,:,0],vmax=data.max())
    plt.subplot(1,2,2)
    plt.imshow(data[:,:,1],vmax=data.max())

plt.figure()
data = np.concatenate([tuning[k][ontarget_ret_lax[k]] for k in (0,2,3,4,5)],axis=0)
for s in range(4):
    this_data = np.nanmean(data[np.argmax(data[:,:,-1,0],axis=1)==s],axis=0)
    for ll in range(2):
        plt.subplot(2,4,ll*4+s+1)
        plt.imshow(this_data[:,:,ll],vmax=this_data.max(),interpolation='bicubic')
        plt.axis('off')

si = 1 - data[:,-1,:,:]/np.max(data,axis=1)

plt.figure()
plt.scatter(si[:,1,0],si[:,-1,0])
plt.plot((0,1),(0,1),c='r')
plt.xlabel('low contrast SI')
plt.ylabel('high contrast SI')

plt.figure()
bins = np.linspace(0,1,100)
h0hi = plt.hist(si[:,-1,0],bins=bins)
h1hi = plt.hist(si[:,-1,1],bins=bins)

plt.figure()
bins = np.linspace(0,1,100)
h0lo = plt.hist(si[:,1,0],bins=bins)
h1lo = plt.hist(si[:,1,1],bins=bins)

for y in (0,1,2):
    data = np.concatenate([tuning[k][ontarget_ret_lax[k]] for k in (y,)],axis=0)
    si = 1 - data[:,-1,:,:]/np.max(data,axis=1)
    plt.figure(figsize=(10,2))
    for i in range(6):
        plt.subplot(1,6,i+1)
        h0lo = plt.hist(si[:,i,0],bins=bins)
        h1lo = plt.hist(si[:,i,1],bins=bins)
        plt.cla()
        plt.plot(bins[1:],np.cumsum(h0lo[0]))
        plt.plot(bins[1:],np.cumsum(h1lo[0]))

# tuning[k][np.logical_and(bottom_3[k],ontarget_ret_lax[keylist[k]])].shape

k = 0
ind = 12
data = tuning[k][ontarget_ret_lax[k]][ind]
plt.figure()
plt.subplot(1,2,1)
plt.imshow(data[:,:,0],vmax=data.max())
plt.subplot(1,2,2)
plt.imshow(data[:,:,1],vmax=data.max())

tuning[k].shape

k = 2
ut.imshow_in_pairs(tuning[k][ontarget_ret_lax[k]][:,:,:,0],tuning[k][ontarget_ret_lax[k]][:,:,:,1])

k = 0
ut.imshow_in_pairs(tuning[k][ontarget_ret_lax[keylist[k]]][:,:,:,0,nbefore:-nafter].mean(-1),tuning[k][ontarget_ret_lax[keylist[k]]][:,:,:,1,nbefore:-nafter].mean(-1))

# os.mkdir('single_neuron_vip_halo_running')
k = 0
for ind in range(ontarget_ret_lax[k].sum()):
    data = tuning[k][ontarget_ret_lax[k]][ind]
    fig = plt.figure(figsize=(10,3))
    for s in range(4):
        plt.subplot(1,4,s+1)
        plt.title(str(int(np.round(usize[s])))+'$^o$')
        plt.plot(100*ucontrast,data[s,:,0],c='k')#,pct=(2.5,97.5))
        plt.plot(100*ucontrast,data[s,:,1],c='r')
        plt.xlabel('contrast (%)')
        plt.ylim((0,1.1*data.max()))
        if s==0:
            plt.ylabel('event rate / max event rate')
            plt.legend(['lights off','lights on'])
        else:
            plt.gca().get_yaxis().set_ticklabels([])
    fig.tight_layout()
#     plt.savefig('single_neuron_vip_halo_running/%03d.png' % ind)

plt.figure()
plt.plot(np.nanmean(np.nanmean(np.nanmean(np.nanmean(tuning[0],0),0),0),0))

k = 0
tuning[k][np.logical_and(bottom_3[k],ontarget_ret_lax[keylist[k]])][:,:,:,:,nbefore:-nafter].mean(-1).shape
# ontarget_ret_laxs = ontarget_ret_lax.copy()
# ontarget_ret_laxs[keylist[1]] = np.ones(ontarget_ret_lax[keylist[1]].shape,dtype='bool')

tuning[k][ontarget_ret_lax[k]].shape

usize0 = np.concatenate(((0,),usize))
for k in (3,): #range(len(tuning)):
    #data = np.concatenate([np.nanmean(tuning[k][ontarget_ret_lax[k]][:,:,:,:,nbefore:-nafter],-1) for k in (0,1,2)],axis=0).mean(0)
    data = np.concatenate([tuning[k][ontarget_ret_lax[k]] for k in (k,)],axis=0).mean(0)
    plt.figure(figsize=(10,3))
    for c in range(6):
        plt.subplot(1,6,c+1)
        plt.plot(usize0,np.concatenate(((np.nanmean(data[:,0,0]),),data[:,c,0])))
        plt.plot(usize0,np.concatenate(((np.nanmean(data[:,0,1]),),data[:,c,1])))
        plt.ylim((np.nanmin(data),1.1*np.nanmax(data)))
        plt.axis('off')

    plt.figure(figsize=(10,3))
    for s in range(4):
        plt.subplot(1,4,s+1)
        plt.plot(ucontrast,data[s,:,0])
        plt.plot(ucontrast,data[s,:,1])
        plt.ylim((0,np.nanmax(data)))
        plt.axis('off')

# +
# ontarget_temp = [None]*3
# for k in range(3):
# #     ontarget_temp[k] = tuning[k][:,0,-1,0]>tuning[k][:,:,0,0].mean(1)
#     ontarget_temp[k] = pval_size_contrast[keylist[k]][1].min(1)<0.05

# +
# for k in range(len(keylist)):
#     print((pval_size_contrast[keylist[k]][1].min(1)<0.05/4).sum())

# +
# [(pval_size_contrast[keylist[k]][0].min(1)<0.05).sum() for k in range(3)]
# -

usize0 = np.concatenate(((0,),usize))
for kk in range(1):
    #data = np.concatenate([np.nanmean(tuning[k][ontarget_ret_lax[k]][:,:,:,:,nbefore:-nafter],-1) for k in (0,1,2)],axis=0).mean(0)
    data = np.concatenate([tuning[k][ontarget_ret_lax[k]] for k in (0,2,3,5)],axis=0).mean(0)
    plt.figure(figsize=(10,3))
    for c in range(6):
        plt.subplot(1,6,c+1)
        plt.plot(usize0,np.concatenate(((np.nanmean(data[:,0,0]),),data[:,c,0])))
        plt.plot(usize0,np.concatenate(((np.nanmean(data[:,0,1]),),data[:,c,1])))
        plt.ylim((np.nanmin(data),1.1*np.nanmax(data)))
        plt.axis('off')

    plt.figure(figsize=(10,3))
    for s in range(4):
        plt.subplot(1,4,s+1)
        plt.plot(ucontrast,data[s,:,0])
        plt.plot(ucontrast,data[s,:,1])
        plt.ylim((0,np.nanmax(data)))
        plt.axis('off')

for klist in ((0,),(1,),(2,)):
    #data = np.concatenate([np.nanmean(tuning[k][ontarget_ret_lax[keylist[k]]][:,:,:,:,nbefore:-nafter],-1) for k in klist],axis=0).mean(0)
    data = np.concatenate([tuning[k][ontarget_temp[k]] for k in klist],axis=0).mean(0)
    plt.figure(figsize=(10,3))
    for c in range(6):
        plt.subplot(1,6,c+1)
        plt.plot(usize0,np.concatenate(((np.nanmean(data[:,0,0]),),data[:,c,0])))
        plt.plot(usize0,np.concatenate(((np.nanmean(data[:,0,1]),),data[:,c,1])))
        plt.ylim((np.nanmin(data),1.1*np.nanmax(data)))
        plt.axis('off')

    plt.figure(figsize=(10,3))
    for s in range(4):
        plt.subplot(1,4,s+1)
        plt.plot(ucontrast,data[s,:,0])
        plt.plot(ucontrast,data[s,:,1])
        plt.ylim((0,np.nanmax(data)))
        plt.axis('off')

plt.figure(figsize=(10,3))
for s in range(2):
    plt.subplot(1,2,s+1)
    plt.plot(ucontrast,data[2*s,:,0])
    plt.plot(ucontrast,data[2*s,:,1])
    plt.ylim((0,data.max()))
    plt.axis('off')

# +
data = tuning[k][ontarget_ret_lax[keylist[k]]][:,:,:,:,nbefore:-nafter].mean(0).mean(-1)
plt.figure(figsize=(10,3))
for c in range(6):
    plt.subplot(1,6,c+1)
    plt.plot(usize0,np.concatenate(((data[:,0,0].mean(),),data[:,c,0])))
    plt.plot(usize0,np.concatenate(((data[:,0,1].mean(),),data[:,c,1])))
    plt.ylim((data.min(),1.1*data.max()))
    plt.axis('off')
    
plt.figure(figsize=(10,3))
for s in range(4):
    plt.subplot(1,4,s+1)
    plt.plot(ucontrast,data[s,:,0])
    plt.plot(ucontrast,data[s,:,1])
    plt.ylim((0,data.max()))
    plt.axis('off')

# +
data = tuning[k][ontarget_ret_lax[keylist[k]]][:,:,:,:,nbefore:-nafter].mean(0).mean(-1)
# data = np.concatenate([np.nanmean(tuning[k][np.logical_and(bottom_3[k],ontarget_ret_lax[keylist[k]])][:,:,:,:,nbefore:-nafter],-1) for k in range(1,len(keylist))],axis=0).mean(0)

plt.figure(figsize=(6,3))
for s in (0,2):
    plt.subplot(1,2,1)
    plt.title('lights off, on target')
    plt.plot(ucontrast,data[s,:,0]) # ucontrast,
    plt.ylim((0,data.max()))
    plt.axis('off')
    plt.subplot(1,2,2)
    plt.title('lights on')
    plt.plot(ucontrast,data[s,:,1]) # ucontrast,
    plt.ylim((0,data.max()))
    plt.axis('off')
plt.subplot(1,2,2)
plt.legend(['8 deg','22 deg'])
# plt.savefig('vip_halo_ontarget_contrast.png')

# data = tuning[offtarget_ret_lax[keylist[k]]][:,:,:,:,nbefore:-nafter].mean(0).mean(-1)

# plt.figure(figsize=(6,3))
# for s in (0,2):
#     plt.subplot(1,2,1)
#     plt.title('lights off, off target')
#     plt.plot(data[s,:,0]) # ucontrast,
#     plt.ylim((0,data.max()))
#     plt.axis('off')
#     plt.subplot(1,2,2)
#     plt.title('lights on')
#     plt.plot(data[s,:,1]) # ucontrast,
#     plt.ylim((0,data.max()))
#     plt.axis('off')
# plt.subplot(1,2,2)
# plt.legend(['8 deg','22 deg'])
# plt.savefig('vip_halo_offtarget_contrast.png')

# +
data = tuning[k][ontarget_ret_lax[keylist[k]]][:,:,:,:,nbefore:-nafter].mean(0).mean(-1)
# data = np.concatenate([np.nanmean(tuning[k][np.logical_and(bottom_3[k],ontarget_ret_lax[keylist[k]])][:,:,:,:,nbefore:-nafter],-1) for k in range(1,len(keylist))],axis=0).mean(0)

plt.figure(figsize=(6,3))
for s in (1,3,5):
    plt.subplot(1,2,1)
    plt.title('lights off, on target')
    plt.plot(usize,data[:,s,0].T) # ucontrast,
    plt.ylim((0,data.max()))
    plt.axis('off')
    plt.subplot(1,2,2)
    plt.title('lights on')
    plt.plot(usize,data[:,s,1].T) # ucontrast,
    plt.ylim((0,data.max()))
    plt.axis('off')
#plt.subplot(1,2,2)
#plt.legend(['8 deg','22 deg'])
# plt.savefig('vip_halo_ontarget_contrast.png')

# data = tuning[offtarget_ret_lax[keylist[k]]][:,:,:,:,nbefore:-nafter].mean(0).mean(-1)

# plt.figure(figsize=(6,3))
# for s in (0,2):
#     plt.subplot(1,2,1)
#     plt.title('lights off, off target')
#     plt.plot(data[s,:,0]) # ucontrast,
#     plt.ylim((0,data.max()))
#     plt.axis('off')
#     plt.subplot(1,2,2)
#     plt.title('lights on')
#     plt.plot(data[s,:,1]) # ucontrast,
#     plt.ylim((0,data.max()))
#     plt.axis('off')
# plt.subplot(1,2,2)
# plt.legend(['8 deg','22 deg'])
# plt.savefig('vip_halo_offtarget_contrast.png')
# -

data.shape

k = 1
data = np.concatenate([np.nanmean(tuning[k][:,:,:,:,nbefore:-nafter],-1) for k in range(1,len(keylist))],axis=0)
data = data/data.max(-1).max(-1).max(-1)[:,np.newaxis,np.newaxis,np.newaxis]
ut.plot_in_rows(data[:,::2,:,0].transpose((0,2,1)))

ontarget_ret_lax[keylist[1]] = ontarget

import naka_rushton_analysis as nra

reload(ut)
reload(nra)

kk = (0,)
data = np.concatenate([tuning[k][ontarget_temp[k]] for k in kk],axis=0)
c = np.array((0,6,12,25,50,100))
for ind in (17,):
    plt.figure(figsize=(6,3))
    for ll in range(2):
        plt.subplot(1,2,ll+1)
        plt.plot(c,data[ind,:,:,ll].T)
        plt.ylim((0,1))

kk = (2,)
data = np.concatenate([tuning[k][ontarget_temp[k]] for k in kk],axis=0)
c = np.array((0,6,12,25,50,100))
for ind in (15,34,49):
    plt.figure(figsize=(6,3))
    for ll in range(2):
        plt.subplot(1,2,ll+1)
        plt.plot(c,data[ind,:,:,ll].T)
        plt.ylim((0,1))

foldname = 'all_size_contrast_lights_off_animals'
if not os.path.exists(foldname):
    os.mkdir(foldname)
for ko in range(3):
    kk = (ko,)
    data = np.concatenate([tuning[k][ontarget_temp[k]] for k in kk],axis=0)
    ut.imshow_in_rows(data[:,:,:,0])
    plt.savefig(foldname+'/all_size_contrast_lights_off'+str(kk[0])+'.pdf')

kk = (0,)
foldname = 'size_contrast_example_neurons/'
if not os.path.exists(foldname):
    os.mkdir(foldname)
data = np.concatenate([tuning[k][ontarget_temp[k]] for k in kk],axis=0)
for ind in (3,60):
    plt.figure()
    plt.imshow(data[ind,:,:,0])
    plt.yticks(np.arange(4),[8,13,22,36])
    plt.xticks(np.arange(6),[0,6,12,25,50,100])
    plt.ylabel('size ($^o$)')
    plt.xlabel('contrast (%)')
    plt.savefig(foldname+str(kk[0])+'_'+str(ind)+'.pdf')

kk = (2,)
data = np.concatenate([tuning[k][ontarget_temp[k]] for k in kk],axis=0)
ut.imshow_in_rows(data[:,:,:,0])

ut.imshow_in_pairs(data[:,:,:,0],data[:,:,:,1])

# +
    # data = tuning[ontarget_ret_lax[keylist[k]]][:,:,:,:,nbefore:-nafter].mean(-1)
#     data = np.concatenate([np.nanmean(tuning[k][ontarget_temp[k]][:,:,:,:,nbefore:-nafter],-1) for k in range(len(keylist))],axis=0)
#     data = np.concatenate([np.nanmean(tuning[k][:,:,:,:,nbefore:-nafter],-1) for k in range(len(keylist))],axis=0)
    if not os.path.exists('vip_halo_size_tuning_animals'):
        os.mkdir('vip_halo_size_tuning_animals')
    for kk in [(k,) for k in range(len(keylist))]:
        data = np.concatenate([tuning[k][ontarget_temp[k]] for k in kk],axis=0)
        data = data/data.max(-1).max(-1).max(-1)[:,np.newaxis,np.newaxis,np.newaxis]
        mn = data.mean(0)
        lb,ub = ut.bootstrap(data,fn=np.mean,axis=0,pct=(2.5,97.5))
        gray = np.tile(data[:,:,0].mean(1)[:,np.newaxis,np.newaxis],(1,1,data.shape[2],1))
        # data = data.mean(0)
        fig = plt.figure(figsize=(12,3))
        for c in range(6):
            plt.subplot(1,6,c+1)
            plt.title(str(int(np.round(100*ucontrast[c])))+'%')
            ut.plot_bootstrapped_errorbars_hillel(usize0,np.concatenate((gray,data),axis=1)[:,:,c].transpose((0,2,1)),colors=['k','r'],pct=(2.5,97.5),markersize=10)
    #         plt.ylim((0,1))
            if c==0:
                plt.ylabel('event rate / max event rate')
                plt.legend(['lights off','lights on'])
            else:
                a = 0
#                 plt.gca().get_yaxis().set_ticklabels([])
            plt.xlabel('size (deg)')
        fig.tight_layout()
        plt.savefig('vip_halo_size_tuning_animals/size_tuning_by_contrast_vip_halo_running_' + str(kk[0]) + '.pdf')

# #     fig = plt.figure(figsize=(10,3))
#     for s in range(4):
# #         plt.subplot(1,4,s+1)
#         plt.figure()
#         plt.title(str(int(np.round(usize[s])))+'$^o$')
#         plt.plot((0,0),(0.5,0.5),c='k')
#         plt.plot((0,0),(0.5,0.5),c='r')
#         ut.plot_bootstrapped_errorbars_hillel(100*ucontrast,data[:,s,:].transpose((0,2,1)),colors=['k','r'],pct=(2.5,97.5),markersize=10,linewidth=0)
#         c = 100*np.logspace(-20,0,500)
#         params = nra.fit_opt_params(100*ucontrast,np.nanmean(data[:,s,:].transpose((0,2,1)),0)[0])
#         plt.plot(c,nra.naka_rushton(c,params),'k')
#         params = nra.fit_opt_params(100*ucontrast,np.nanmean(data[:,s,:].transpose((0,2,1)),0)[1])
#         plt.plot(c,nra.naka_rushton(c,params),'r')
#         plt.xlabel('contrast (%)')
# #         plt.ylim((0,1))
# #         if s==0:
#         plt.ylabel('event rate / max event rate')
#         plt.legend(['lights off','lights on'])
#         plt.ylim((0.15,0.75))
# #         plt.savefig('opto_comparison_'+str(s)+'_contrast_naka_rushton.pdf')
# #         else:
# #             plt.gca().get_yaxis().set_ticklabels([])
#     fig.tight_layout()
#     plt.savefig('contrast_tuning_by_size_vip_halo_running_1.png')

# -

pdb.pm()

data.shape

[np.nanmean(tuning[k][ontarget_ret_lax[keylist[k]]][:,:,:,:,nbefore:-nafter],-1).shape[0] for k in range(len(keylist))]

    reload(ut)
    if not os.path.exists('vip_halo_running'):
        os.mkdir('vip_halo_running')
    data = np.concatenate([tuning[k][ontarget_ret_lax[k]] for k in (3,)],axis=0) #(0,2,3,5)],axis=0)
#     data = np.concatenate([np.nanmean(tuning[k][:,:,:,:,nbefore:-nafter],-1) for k in range(len(keylist))],axis=0)
    data = data/data.max(-1).max(-1).max(-1)[:,np.newaxis,np.newaxis,np.newaxis]
    mn = data.mean(0)
    lb,ub = ut.bootstrap(data,fn=np.mean,axis=0,pct=(2.5,97.5))
    gray = np.tile(data[:,:,0].mean(1)[:,np.newaxis,np.newaxis],(1,1,data.shape[2],1))
    # data = data.mean(0)
    for c in range(6):
        fig = plt.figure()
        plt.title(str(int(np.round(100*ucontrast[c])))+'%')
        ut.plot_bootstrapped_errorbars_hillel(usize0,np.concatenate((gray,data),axis=1)[:,:,c].transpose((0,2,1)),colors=[np.array((0,0,0,(5+1)/6)),np.array((1,0,0,(5+1)/6))],pct=(2.5,97.5),norm_to_max=False)
#         plt.ylim((0,1))
        plt.ylabel('event rate / max event rate')
        plt.legend(['lights off','lights on'])
        plt.xlabel('size (deg)')
        fig.tight_layout()
        plt.savefig('vip_halo_running/size_tuning_by_contrast_'+str(c)+'.pdf')


    reload(ut)
    if not os.path.exists('vip_halo_running'):
        os.mkdir('vip_halo_running')
    data = np.concatenate([tuning[k][ontarget_ret_lax[k]] for k in (3,)],axis=0) #(0,2,3,5)],axis=0)
#     data = np.concatenate([np.nanmean(tuning[k][:,:,:,:,nbefore:-nafter],-1) for k in range(len(keylist))],axis=0)
    data = data/data.max(-1).max(-1).max(-1)[:,np.newaxis,np.newaxis,np.newaxis]
    mn = data.mean(0)
    lb,ub = ut.bootstrap(data,fn=np.mean,axis=0,pct=(2.5,97.5))
    gray = np.tile(data[:,:,0].mean(1)[:,np.newaxis,np.newaxis],(1,1,data.shape[2],1))
    # data = data.mean(0)
    for c in range(6):
        fig = plt.figure()
        plt.title(str(int(np.round(100*ucontrast[c])))+'%')
        ut.plot_bootstrapped_errorbars_hillel(usize0,np.concatenate((gray,data),axis=1)[:,:,c].transpose((0,2,1)),colors=[np.array((0,0,0,(5+1)/6)),np.array((1,0,0,(5+1)/6))],pct=(2.5,97.5),norm_to_max=True)
#         plt.ylim((0,1))
        plt.ylabel('evoked event rate / max event rate')
        plt.legend(['lights off','lights on'])
        plt.xlabel('size (deg)')
        fig.tight_layout()
        plt.savefig('vip_halo_running/size_tuning_by_contrast_'+str(c)+'_peak_norm.pdf')


np.concatenate((gray,data),axis=1)[:,:,c].transpose((0,2,1)).shape

    # data = tuning[ontarget_ret_lax[keylist[k]]][:,:,:,:,nbefore:-nafter].mean(-1)
    data = np.concatenate([np.nanmean(tuning[k][ontarget[keylist[k]]][:,:,:,:,nbefore:-nafter],-1) for k in range(len(keylist))],axis=0)
#     data = np.concatenate([np.nanmean(tuning[k][:,:,:,:,nbefore:-nafter],-1) for k in range(len(keylist))],axis=0)
    data = data/data.max(-1).max(-1).max(-1)[:,np.newaxis,np.newaxis,np.newaxis]
    mn = data.mean(0)
    lb,ub = ut.bootstrap(data,fn=np.mean,axis=0,pct=(2.5,97.5))
    gray = np.tile(data[:,:,0].mean(1)[:,np.newaxis,np.newaxis],(1,1,data.shape[2],1))
    pickout = np.array((1,3,5))
    # data = data.mean(0)
    fig = plt.figure()#figsize=(12,3))
    for ll in (0,1):
        plt.subplot(1,2,ll+1)
        for c in pickout:
            #plt.subplot(1,6,c+1)
            plt.title(['lights off','lights on'][ll])
            ut.plot_bootstrapped_errorbars_hillel(usize0,np.concatenate((gray,data),axis=1)[:,:,c].transpose((0,2,1))[:,ll:ll+1],colors=[(ll,0,0,(c+1)/6)],pct=(2.5,97.5))
            plt.ylim((0,1))
        
            if c==0:
                plt.ylabel('event rate / max event rate')
            else:
                plt.gca().get_yaxis().set_ticklabels([])
        plt.xlabel('size ($^o$)')
        plt.legend([str(int(x))+'%' for x in 100*ucontrast[pickout]])
    fig.tight_layout()
    plt.savefig('size_tuning_by_contrast_vip_halo_running.pdf')


fig = plt.figure(figsize=(6,3))
for s in range(2):
    plt.subplot(1,2,s+1)
    plt.title(str(int(np.round(usize[2*s])))+'$^o$')
    ut.plot_bootstrapped_errorbars_hillel(np.arange(6),data[:,2*s,:].transpose((0,2,1)),colors=['k','r'],pct=(2.5,97.5))
    plt.xticks(np.arange(6),[str(x) for x in (0,6,12,25,50,100)])
    plt.xlabel('contrast (%)')
    plt.ylim((0,1))
    if s==0:
        plt.ylabel('event rate / max event rate')
        plt.legend(['lights off','lights on'])
    else:
        plt.gca().get_yaxis().set_ticklabels([])
fig.tight_layout()
# plt.savefig('vip_halo_8_22_by_size.pdf')
# plt.savefig('vip_halo_8_22_by_size.png')

pcutoff = 0.05
ontarget = ontarget_ret_lax.copy()
ontarget[keylist[1]] = pval[keylist[1]][1].min(1) < pcutoff/4

fig = plt.figure(figsize=(6,3))
for ll in range(2):
    plt.subplot(1,2,ll+1)
    if ll==0:
        plt.title('lights off')
    else:
        plt.title('lights on')
    ut.plot_bootstrapped_errorbars_hillel(100*ucontrast,data[:,::3,:,ll],pct=(2.5,97.5))
#     plt.xticks(ucontrast,[str(x) for x in (0,6,12,25,50,100)])
    plt.xlabel('contrast (%)')
    plt.ylim((0,1))
    if ll==0:
        plt.ylabel('event rate / max event rate')
        plt.legend(['8$^o$','36$^o$'])
    else:
        plt.gca().get_yaxis().set_ticklabels([])
fig.tight_layout()
plt.savefig('vip_halo_8_36_by_opto.pdf')
# plt.savefig('vip_halo_8_22_by_opto.png')

# +
data = tuning[ontarget_ret_lax[keylist[0]]][:,:,:,:,nbefore:-nafter].mean(-1)
data = data/data.max(-1).max(-1).max(-1)[:,np.newaxis,np.newaxis,np.newaxis]
# mn = data.mean(0)
# lb,ub = ut.bootstrap(data,fn=np.mean,axis=0,pct=(2.5,97.5))
gray = np.tile(data[:,:,0].mean(1)[:,np.newaxis,np.newaxis],(1,1,data.shape[2],1))
# # data = data.mean(0)
fig = plt.figure(figsize=(12,3))
for c in range(6):
    plt.subplot(1,6,c+1)
    plt.title(str(int(np.round(100*ucontrast[c])))+'%')
    ut.plot_bootstrapped_errorbars_hillel(usize0,np.concatenate((gray,data),axis=1)[:,:,c:c+1,0].transpose((0,2,1))-np.concatenate((gray,data),axis=1)[:,:,c:c+1,1].transpose((0,2,1)),pct=(2.5,97.5))
    plt.ylim((0,1))
    if c==0:
        plt.ylabel('event rate / max event rate')
        plt.legend(['lights off','lights on'])
    else:
        plt.gca().get_yaxis().set_ticklabels([])
    plt.xlabel('size (deg)')
fig.tight_layout()
# # plt.savefig('size_tuning_by_contrast_vip_halo.png')

    
fig = plt.figure(figsize=(10,3))
for s in range(4):
    plt.subplot(1,4,s+1)
    plt.title(str(int(np.round(usize[s])))+'$^o$')
    ut.plot_bootstrapped_errorbars_hillel(100*ucontrast,data[:,s:s+1,:,0]-data[:,s:s+1,:,1],pct=(2.5,97.5))
    plt.xlabel('contrast (%)')
    plt.ylim((0,1))
    if s==0:
        plt.ylabel('event rate / max event rate')
#         plt.legend(['lights off','lights on'])
    else:
        plt.gca().get_yaxis().set_ticklabels([])
fig.tight_layout()
# -

data[:,s,:,0].shape

np.concatenate((gray,data),axis=1).transpose((0,3,1,2)).shape

# +
data = tuning[ontarget_ret_lax[keylist[k]]][:,:,:,:,nbefore:-nafter].mean(-1)
data = data/data.max(-1).max(-1).max(-1)[:,np.newaxis,np.newaxis,np.newaxis]
mn = data.mean(0)
lb,ub = ut.bootstrap(data,fn=np.mean,axis=0,pct=(2.5,97.5))
# data = data.mean(0)
plt.figure(figsize=(10,3))
for c in range(6):
    plt.subplot(1,6,c+1)
    plt.title(str(int(np.round(100*ucontrast[c])))+'%')
#     plt.plot(usize0,np.concatenate(((mn[:,0,0].mean(),),mn[:,c,0])))
    plt.fill_between(usize0,np.concatenate(((lb[:,0,0].mean(),),lb[:,c,0])),np.concatenate(((ub[:,0,0].mean(),),ub[:,c,0])),color='k')
#     plt.plot(usize0,np.concatenate(((mn[:,0,1].mean(),),mn[:,c,1])))
    plt.fill_between(usize0,np.concatenate(((lb[:,0,1].mean(),),lb[:,c,1])),np.concatenate(((ub[:,0,1].mean(),),ub[:,c,1])),color='r')
    plt.ylim((lb.min()-0.05,ub.max()+0.05))
    plt.axis('off')
plt.subplot(1,6,1)
plt.axis('on')
plt.xlabel('size (deg)')
plt.ylabel('event rate / max event rate')
    
plt.figure(figsize=(10,3))
for s in range(4):
    plt.subplot(1,4,s+1)
    plt.title(str(int(np.round(usize[s])))+'$^o$')
#     plt.plot(ucontrast,mn[s,:,0])
    plt.fill_between(100*ucontrast,lb[s,:,0],ub[s,:,0],color='k')
    plt.fill_between(100*ucontrast,lb[s,:,1],ub[s,:,1],color='r')
    plt.ylim((lb.min()-0.05,ub.max()+0.05))
#     plt.axis('off')
# plt.subplot(1,4,1)
# plt.axis('on')
    plt.xlabel('contrast (%)')
    if s==0:
        plt.ylabel('event rate / max event rate')
    else:
        plt.yticks
# -

plt.figure()
colors = plt.cm.viridis(np.linspace(0,1,4))
for i in range(4):
    plt.scatter(mn[i,:,0].flatten(),mn[i,:,1].flatten(),c=colors[i])
    plt.plot(mn[i,:,0].flatten(),mn[i,:,1].flatten(),c=colors[i])
plt.xlabel('event rate/max event rate, lights off')
plt.ylabel('event rate/max event rate, lights on')
# plt.axis('equal')
plt.legend([str(int(np.round(s)))+'$^o$' for s in usize])
plt.plot(mn[:,:,0].flatten(),mn[:,:,0].flatten(),c='k')
plt.title('Contrast response functions with VIP silencing')
# plt.savefig('crf_vip_halo_lights_on_off.png')

plt.figure()
# plt.plot(mn[:,:,0].flatten(),mn[:,:,0].flatten(),c='k')
colors = plt.cm.viridis(np.linspace(0,1,6))
for i in range(6):
    plt.scatter(mn[:,i,0].flatten(),mn[:,i,1].flatten(),c=colors[i][np.newaxis])
    plt.plot(mn[:,i,0].flatten(),mn[:,i,1].flatten(),c=colors[i])
plt.plot((0,1),(0,1))
plt.xlabel('PC event rate, lights off')
plt.ylabel('PC event rate, lights on')
# plt.axis('equal')
plt.legend([str(int(np.round(100*s)))+'%' for s in ucontrast])
# plt.title('Size tuning functions with VIP silencing')

colors[i].shape

mn.shape

# +
data = tuning[ontarget_ret_lax[keylist[0]]][:,:,:,:,nbefore:-nafter].mean(-1)
# data = data/data.max(1).max(1)[:,np.newaxis,np.newaxis]
data = (data-data.min(1).min(1)[:,np.newaxis,np.newaxis])/(data.max(1).max(1)-data.min(1).min(1))[:,np.newaxis,np.newaxis]
mn = data.mean(0)
lb,ub = ut.bootstrap(data,fn=np.mean,axis=0,pct=(2.5,97.5))
# data = data.mean(0)
plt.figure(figsize=(10,3))
for c in range(6):
    plt.subplot(1,6,c+1)
#     plt.plot(usize0,np.concatenate(((mn[:,0,0].mean(),),mn[:,c,0])))
    plt.fill_between(usize0,-lb[:,0,0].mean()+np.concatenate(((lb[:,0,0].mean(),),lb[:,c,0])),-lb[:,0,0].mean()+np.concatenate(((ub[:,0,0].mean(),),ub[:,c,0])))
#     plt.plot(usize0,np.concatenate(((mn[:,0,1].mean(),),mn[:,c,1])))
    plt.fill_between(usize0,-lb[:,0,1].mean()+np.concatenate(((lb[:,0,1].mean(),),lb[:,c,1])),-lb[:,0,1].mean()+np.concatenate(((ub[:,0,1].mean(),),ub[:,c,1])))
    plt.ylim((0,1))
    plt.axis('off')
    
plt.figure(figsize=(10,3))
for s in range(4):
    plt.subplot(1,4,s+1)
#     plt.plot(ucontrast,mn[s,:,0])
    plt.fill_between(ucontrast,lb[s,:,0]-lb[s,0,0],ub[s,:,0]-lb[s,0,0])
    plt.fill_between(ucontrast,lb[s,:,1]-lb[s,0,1],ub[s,:,1]-lb[s,0,1])
    plt.ylim((0,1))
    plt.axis('off')
# -

(lb[s,:,1]-lb[s,0,1]).shape

plt.figure(figsize=(8,4))
plt.subplot(1,2,1)
for s in range(4):
    plt.plot(ucontrast,data[s,:,0])
    plt.ylim((0,data.max()))
    plt.axis('off')
plt.subplot(1,2,2)
for s in range(4):
    plt.plot(ucontrast,data[s,:,1])
    plt.ylim((0,data.max()))
    plt.axis('off')
#     plt.plot(data[s,:,1])

plt.figure(figsize=(6,3))
for c in range(6):
    plt.subplot(1,2,1)
    plt.plot(usize0,np.concatenate(((data[:,0,0].mean(),),data[:,c,0])))
    plt.ylim((data.min(),1.1*data.max()))
    plt.axis('off')
    plt.subplot(1,2,2)
    plt.plot(usize0,np.concatenate(((data[:,0,1].mean(),),data[:,c,1])))
    plt.ylim((data.min(),1.1*data.max()))
    plt.axis('off')

data[s,:,0].max()

plt.figure(figsize=(8,4))
plt.subplot(1,2,1)
for s in range(4):
    plt.plot(ucontrast,ut.norm01(data[s,:,0]))
    plt.ylim((0,1))
    plt.axis('off')
plt.subplot(1,2,2)
for s in range(4):
    plt.plot(ucontrast,ut.norm01(data[s,:,1]))
    plt.ylim((0,1))
    plt.axis('off')

(ret_vars[keylist[0]]['pval_ret']<0.01).sum()


# +
def symmetrize(arr):
    return np.concatenate((arr,arr[0:1]))

os.mkdir('size_contrast_vip_halo')

for ind in range(strialavg.shape[0]):
    if ret_vars[keylist[0]]['pval_ret'][ind]<0.01:
        fig = plt.figure()
        nsize = strialavg.shape[1]
        ncontrast = strialavg.shape[2]
        nangle = strialavg.shape[3]
        data = strialavg[ind,:,:,:,:,nbefore:-nafter].mean(-1)
        for i in range(nsize):
            for j in range(ncontrast):
                plt.subplot(nsize,ncontrast,i*ncontrast+j+1)
                plt.fill_between(np.arange(nangle+1),symmetrize(data[i,j,:,0]),alpha=0.5)
                plt.plot(np.arange(nangle+1),symmetrize(data[i,j,:,0]))
                plt.ylim(data.min(),data.max())
                plt.fill_between(np.arange(nangle+1),symmetrize(data[i,j,:,1]),alpha=0.5)
                plt.plot(np.arange(nangle+1),symmetrize(data[i,j,:,1]))
                plt.axis('off')
    #     fig.patch.set_facecolor('white')
        plt.savefig('size_contrast_vip_halo/%04d.tif' % ind, transparent=False)
# -

plt.figure(figsize=(4,1))
for s in range(4):
    plt.subplot(1,4,s+1)
    plt.plot(data[s,:,0])
    plt.plot(data[s,:,1])
    plt.ylim((0,data.max()))
    plt.axis('off')

c = 2
plt.figure()
plt.plot(usize,data[:,c,0])
plt.plot(usize,data[:,c,1])
plt.ylim((data.min(),0.8))

plt.figure(figsize=(5,15))
for ind in range(tuning.shape[0]):
    plt.subplot(38,10,ind+1)
    plt.imshow(tuning[ind,:,:,0,nbefore:-nafter].mean(-1))
    #plt.plot(tuning[ind,:,1,nbefore:-nafter].mean(-1))
    plt.axis('off')

plt.figure(figsize=(5,15))
for ind in range(tuning.shape[0]):
    plt.subplot(38,10,ind+1)
    plt.plot(tuning[ind,:,0,nbefore:-nafter].mean(-1))
    plt.plot(tuning[ind,:,1,nbefore:-nafter].mean(-1))
    plt.axis('off')

plt.figure()#figsize=(5,15))
for ind in range(tuning[redratio>0.65].shape[0]):
    plt.subplot(3,10,ind+1)
    plt.plot(tuning[np.where(redratio>0.65)[0]][ind,:,0,nbefore:-nafter].mean(-1))
    plt.plot(tuning[np.where(redratio>0.65)[0]][ind,:,1,nbefore:-nafter].mean(-1))
    plt.axis('off')

ops = np.load('/media/mossing/backup_1/data/suite2P/results/M9835/190227/3_4/suite2p/ops1.npy')

mred = ops[0]['meanImg_chan2_corrected']
mgreen = ops[0]['meanImg_chan2']

ops[0].keys()

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

tuning.shape

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

for ind in range(10):
    plt.figure()
    for i,rs in enumerate((0,1)):
        data = soriavg[key][rs][ontarget_ret_lax[key]][ind]
        plt.subplot(2,2,2*i+1)
        plt.imshow(data[:,:,0],vmax=data.max())
        plt.subplot(2,2,2*i+2)
        plt.imshow(data[:,:,1],vmax=data.max())



# +
rs = 0
data = soriavg[key][rs][ontarget_ret_lax[key]].mean(0)
plt.figure()
plt.subplot(1,2,1)
plt.imshow(data[:,:,0],vmax=data.max(),vmin=data.min())
plt.subplot(1,2,2)
plt.imshow(data[:,:,1],vmax=data.max(),vmin=data.min())

plt.figure()
plt.plot(data[0,:,0])
plt.plot(data[0,:,1])

plt.figure()
plt.plot(data[1,:,0])
plt.plot(data[1,:,1])

plt.figure()
plt.plot(data[2,:,0])
plt.plot(data[2,:,1])
# -



# +
plt.figure()
ind = 0
plt.subplot(1,2,1)
plt.imshow(soriavg[key][runstatus][ind][:,:,0],vmax=soriavg[key][1][ind].max())
plt.subplot(1,2,2)
plt.imshow(soriavg[key][runstatus][ind][:,:,1],vmax=soriavg[key][1][ind].max())

plt.figure()
ind = 6
plt.subplot(1,2,1)
plt.imshow(soriavg[key][runstatus][ind][:,:,0],vmax=soriavg[key][1][ind].max())
plt.subplot(1,2,2)
plt.imshow(soriavg[key][runstatus][ind][:,:,1],vmax=soriavg[key][1][ind].max())

plt.figure()
ind = 11
plt.subplot(1,2,1)
plt.imshow(soriavg[key][runstatus][ind][:,:,0],vmax=soriavg[key][1][ind].max())
plt.subplot(1,2,2)
plt.imshow(soriavg[key][runstatus][ind][:,:,1],vmax=soriavg[key][1][ind].max())

plt.figure()
ind = 15
plt.subplot(1,2,1)
plt.imshow(soriavg[key][runstatus][ind][:,:,0],vmax=soriavg[key][1][ind].max())
plt.subplot(1,2,2)
plt.imshow(soriavg[key][runstatus][ind][:,:,1],vmax=soriavg[key][1][ind].max())

plt.figure()
ind = 16
plt.subplot(1,2,1)
plt.imshow(soriavg[key][runstatus][ind][:,:,0],vmax=soriavg[key][1][ind].max())
plt.subplot(1,2,2)
plt.imshow(soriavg[key][runstatus][ind][:,:,1],vmax=soriavg[key][1][ind].max())

plt.figure()
ind = 70
plt.subplot(1,2,1)
plt.imshow(soriavg[key][runstatus][ind][:,:,0],vmax=soriavg[key][1][ind].max())
plt.subplot(1,2,2)
plt.imshow(soriavg[key][runstatus][ind][:,:,1],vmax=soriavg[key][1][ind].max())
# -

thisfold = keylist[1]
strialwise = proc[thisfold]['strialwise']
trialwise = proc[thisfold]['trialwise']
contrast = proc[thisfold]['contrast']
size = proc[thisfold]['size']
trialrun = proc[thisfold]['trialrun']
dfof = proc[thisfold]['dfof']
frame = proc[thisfold]['frame']
angle = proc[thisfold]['angle']
light = proc[thisfold]['light']

trialwise.shape


def compute_oriavg(proc):
#     keylist = list(proc.keys())
    oriavg = {}
    for thisfold in proc.keys():
        
        trialwise = proc[thisfold]['trialwise']
        size = proc[thisfold]['size']
        contrast = proc[thisfold]['contrast']
        light = proc[thisfold]['light']
        trialrun = proc[thisfold]['trialrun']
        
        usize = np.unique(size)
        ucontrast = np.unique(contrast)
        ulight = np.unique(light)
        
        oriavg[thisfold] = np.zeros((trialwise.shape[0],len(usize),len(ucontrast),len(ulight)))
        for i,s in enumerate(usize):
            for j,c in enumerate(ucontrast):
                for k,ell in enumerate(ulight):
                    oriavg[thisfold][:,i,j,k] = np.nanmean(np.nanmean(trialwise[:,np.logical_and(np.logical_and(np.logical_and(size==s,contrast==c),trialrun),light==ell),5:],axis=1),axis=-1)
        
    return oriavg


a = np.array((np.arange(7),np.arange(7,14)))
a.flatten()

lkat = [None]*2
for k in range(2):
    smallest_peak = np.argmax(soriavg[thisfold][1-k][:,:,:,0].reshape((soriavg[thisfold][1].shape[0],-1)),1)==6
    second_smallest_peak = np.argmax(soriavg[thisfold][1-k][:,:,:,0].reshape((soriavg[thisfold][1].shape[0],-1)),1)==13
    third_smallest_peak = np.argmax(soriavg[thisfold][1-k][:,:,:,0].reshape((soriavg[thisfold][1].shape[0],-1)),1)==13
    lkat[k] = np.logical_or(np.logical_or(smallest_peak,second_smallest_peak),third_smallest_peak)
ut.imshow_in_rows(soriavg[thisfold][0][:,:,:,0][lkat[0]])
ut.imshow_in_rows(soriavg[thisfold][0][:,:,:,0][lkat[1]])

reload(sca)

# +
# on target based on size/contrast responses
ucontrast = {}
to_take_c = {}
intersect = np.array(())
for key in keylist:
    gdind = sca.find_gdind(soriavg[key])
    ucontrast[key] = np.unique(proc[key][gdind]['contrast'])
for key in keylist:
    if not intersect.size:
        intersect = ucontrast[key]
    else:
        intersect = np.intersect1d(intersect,ucontrast[key])
for key in keylist:
    to_take_c[key] = np.where(np.in1d(ucontrast[key],intersect))[0]
    
c = 100*intersect
    
usize = {}
to_take_s = {}
intersect = np.array(())
for key in keylist:
    gdind = sca.find_gdind(soriavg[key])
    usize[key] = np.unique(proc[key][gdind]['size'])
for key in keylist:
    if not intersect.size:
        intersect = usize[key]
    else:
        intersect = np.intersect1d(intersect,usize[key])
for key in keylist:
    to_take_s[key] = np.where(np.in1d(usize[key],intersect))[0]

s = intersect
s0 = np.concatenate(((0,),s))

snorm = {}
snorm[thisfold] = [None]*2
for k in range(2):
    lkatdict = {}
    lkatdict[thisfold] = lkat[k]
    snorm[thisfold][k] = sca.get_norm_curves_row_col(soriavg,lkatdict,sizes=to_take_s,contrasts=to_take_c,append_gray=True)[k]
# -

soriavg[thisfold][k][lkat[k]].shape

k = 0
ut.imshow_in_rows(soriavg[thisfold][k][lkat[k]][:,:,:,0])

k = 0
plt.figure()
cap = -2
for i,ind in enumerate(range(lkat[k].sum())):
    plt.subplot(7,10,i+1)
    for lind in (0,1):
        mn_tgt = sca.append_gray(soriavg[keylist[0]][0][lkat[k]][ind][:,:,lind])[:,cap]
        lb_tgt = sca.append_gray(lb[keylist[0]][0][lkat[k]][ind].mean(0)[:,:,lind])[:,cap]
        ub_tgt = sca.append_gray(ub[keylist[0]][0][lkat[k]][ind].mean(0)[:,:,lind])[:,cap]
        ut.plot_errorbar_hillel(s0,mn_tgt.T,lb_tgt.T,ub_tgt.T,markersize=20)#,colors=colors[lind])
#         plt.ylabel('event rate')
#         plt.xlabel('size ($^o$)')
#         plt.legend([str(int(np.round(x)))+'%' for x in c[2:cap+1]])
#         plt.ylim((-0.05,2.4))
    plt.axis('off')

np.where(lkat[k])[0][ind]

135+185

nbydepth[thisfold]

lb[keylist[0]][0][lkat[k]][ind].mean(0)[:,:,lind].shape

# +
# ind = 10
# plt.figure()
# plt.plot(s,soriavg[thisfold][k][lkat[k]][ind,:,:,0])
reload(sca)
# ind = 6
ind = 40
# ind = 63
# ind = 24
k = 0
cs = ['k','r']

colors = plt.cm.viridis(np.linspace(0,1,2))
plt.figure(figsize=(10,2))
for cap in range(2,7):
    plt.subplot(1,5,cap-1)
    for lind in (0,1):
        mn_tgt = sca.append_gray(soriavg[keylist[0]][0][lkat[k]][ind][:,:,lind])[:,cap]
        lb_tgt = sca.append_gray(lb[keylist[0]][0][lkat[k]][ind].mean(0)[:,:,lind])[:,cap]
        ub_tgt = sca.append_gray(ub[keylist[0]][0][lkat[k]][ind].mean(0)[:,:,lind])[:,cap]
        ut.plot_errorbar_hillel(s0,mn_tgt.T,lb_tgt.T,ub_tgt.T,c=cs[lind])#,colors=colors[lind])
        plt.ylim((-0.05,0.8))
    plt.title(str(int(np.round(c[cap])))+'% contrast')

plt.subplot(1,5,1)
plt.ylabel('event rate')
plt.xlabel('size ($^o$)')
# plt.legend([str(int(np.round(x)))+'%' for x in c[2:cap+1]])
plt.savefig('vip_halo_results/example_size_by_contrast/ex1_size_tuning'+str(cap-2)+'.pdf')

# +
c = np.array((0,3,6,12,25,50,100))
s = np.logspace(np.log10(5),np.log10(36.5),5)
thisfold = keylist[0]
labels = ['running','non_running']

for k in range(2):
    plt.figure(figsize=(10,2))
    for i in (0,):
        plt.subplot(1,5,i+1)
        ut.plot_bootstrapped_errorbars_hillel(c,soriavg[thisfold][k][lkat[k]][:,i,:,:].transpose((0,2,1)),pct=(16,84),colors=cs)
        plt.ylim((0,1))
        plt.xlabel('contrast (%)')
        plt.ylabel('event rate')
        plt.title(str(int(np.round(s[i])))+'$^o$ size')
    for i in range(1,soriavg[thisfold][k].shape[1]):
        plt.subplot(1,5,i+1)
        ut.plot_bootstrapped_errorbars_hillel(c,soriavg[thisfold][k][lkat[k]][:,i,:,:].transpose((0,2,1)),pct=(16,84),colors=cs)
        plt.ylim((0,1))
        plt.title(str(int(np.round(s[i])))+'$^o$ size')
        plt.gca().get_yaxis().set_visible(False)
    plt.savefig('vip_halo_results/contrast_by_size_non_normalized_'+labels[k]+'.pdf')

    plt.figure(figsize=(10,2))
    for i in (0,):
        plt.subplot(1,5,i+1)
        ut.plot_bootstrapped_errorbars_hillel(s,soriavg[thisfold][k][lkat[k]][:,:,2+i,:].transpose((0,2,1)),pct=(16,84),colors=cs)
        plt.ylim((0,1))
        plt.xlabel('size ($^o$)')
        plt.ylabel('event rate')
        plt.title(str(int(np.round(c[i+2])))+'% contrast')
    for i in range(1,5):
        plt.subplot(1,5,i+1)
        ut.plot_bootstrapped_errorbars_hillel(s,soriavg[thisfold][k][lkat[k]][:,:,2+i,:].transpose((0,2,1)),pct=(16,84),colors=cs)
        plt.ylim((0,1))
        plt.gca().get_yaxis().set_visible(False)
        plt.title(str(int(np.round(c[i+2])))+'% contrast')
    plt.savefig('vip_halo_results/size_by_contrast_non_normalized_'+labels[k]+'.pdf')

# +
c = np.array((0,3,6,12,25,50,100))
s = np.logspace(np.log10(5),np.log10(36.5),5)
thisfold = keylist[0]

labels = ['running','non_running']

for k in range(2):
    plt.figure(figsize=(10,2))
    for i in (0,):
        plt.subplot(1,5,i+1)
        ut.plot_bootstrapped_errorbars_hillel(c,snorm[thisfold][k][:,1+i,:,:].transpose((0,2,1)),pct=(16,84))
        plt.ylim((0,1))
        plt.xlabel('contrast (%)')
        plt.ylabel('event rate/max event rate')
        plt.title(str(int(np.round(s[i])))+'$^o$ size')
    for i in range(1,soriavg[thisfold][k].shape[1]):
        plt.subplot(1,5,i+1)
        ut.plot_bootstrapped_errorbars_hillel(c,snorm[thisfold][k][:,1+i,:,:].transpose((0,2,1)),pct=(16,84))
        plt.ylim((0,1))
        plt.title(str(int(np.round(s[i])))+'$^o$ size')
        plt.gca().get_yaxis().set_visible(False)
    plt.savefig('vip_halo_results/contrast_by_size_normalized_'+labels[k]+'.pdf')

    plt.figure(figsize=(10,2))
    for i in (0,):
        plt.subplot(1,5,i+1)
        ut.plot_bootstrapped_errorbars_hillel(s0,snorm[thisfold][k][:,:,2+i,:].transpose((0,2,1)),pct=(16,84))
        plt.ylim((0,1))
        plt.xlabel('size ($^o$)')
        plt.ylabel('event rate/max event rate')
        plt.title(str(int(np.round(c[i+2])))+'% contrast')
    for i in range(1,5):
        plt.subplot(1,5,i+1)
        ut.plot_bootstrapped_errorbars_hillel(s0,snorm[thisfold][k][:,:,2+i,:].transpose((0,2,1)),pct=(16,84))
        plt.ylim((0,1))
        plt.gca().get_yaxis().set_visible(False)
        plt.title(str(int(np.round(c[i+2])))+'% contrast')
    plt.savefig('vip_halo_results/size_by_contrast_normalized_'+labels[k]+'.pdf')


# c = np.array((0,3,6,12,25,50,100))
# s = np.logspace(np.log10(5),np.log10(36.5),5)
# thisfold = keylist[0]

# for k in range(2):
#     plt.figure(figsize=(10,2))
#     for i in range(snorm[thisfold][k].shape[1]-1):
#         plt.subplot(1,5,i+1)
#         ut.plot_bootstrapped_errorbars_hillel(c,snorm[thisfold][k][:,1+i,:,:].transpose((0,2,1)),pct=(2.5,97.5))
#         plt.ylim((0,1))

#     plt.figure(figsize=(10,2))
#     for i in range(5):
#         plt.subplot(1,5,i+1)
#         ut.plot_bootstrapped_errorbars_hillel(s0,snorm[thisfold][k][:,:,2+i,:].transpose((0,2,1)),pct=(2.5,97.5))
#         plt.ylim((0,1))
# -

os.mkdir('vip_halo_results')

# +
c = np.array((0,3,6,12,25,50,100))
s = np.logspace(np.log10(5),np.log10(36.5),5)
thisfold = keylist[0]

plt.figure(figsize=(10,2))
for i in range(soriavg[thisfold][0].shape[1]):
    plt.subplot(1,5,i+1)
    ut.plot_bootstrapped_errorbars_hillel(c,soriavg[thisfold][1][lkat][:,i,:,:].transpose((0,2,1)),pct=(2.5,97.5))
    plt.ylim((0,1))
    
plt.figure(figsize=(10,2))
for i in range(5):
    plt.subplot(1,5,i+1)
    ut.plot_bootstrapped_errorbars_hillel(s,soriavg[thisfold][1][lkat][:,:,2+i,:].transpose((0,2,1)),pct=(2.5,97.5))
    plt.ylim((0,1))

# +
plt.figure()
for i,thisfold in enumerate(keylist):
    plt.subplot(1,3,i+1)
    for k in range(soriavg[thisfold].shape[1]):
        plt.plot(soriavg[thisfold][ontarget_ret_lax[thisfold]][:,k,:,0].mean(0))
        
plt.figure()
for i,thisfold in enumerate(keylist):
    plt.subplot(1,3,i+1)
    for k in range(soriavg[thisfold].shape[1]):
        plt.plot(soriavg[thisfold][ontarget_ret_lax[thisfold]][:,k,:,1].mean(0))
# -

thisfold = keylist[0]
plt.figure()
plt.plot(soriavg[thisfold][ontarget_ret_lax[thisfold]][:,3,:,0].mean(0))
plt.plot(soriavg[thisfold][ontarget_ret_lax[thisfold]][:,3,:,1].mean(0))

ontarget_ret_lax[thisfold].sum()

strialavg[thisfold].shape

strialwise.mean(0).mean(0)

plt.figure()
plt.plot(np.nanmean(trialwise[:,light==0],0).mean(0))
plt.plot(np.nanmean(trialwise[:,light>0],0).mean(0))

plt.figure(figsize=(12,6))
for j in range(5):
    plt.subplot(2,5,j+1)
    for i in range(7):
        plt.plot(np.nanmean(np.nanmean(strialavg[thisfold][ontarget_ret_lax[thisfold]][:,:,j,i,0],0),0))
    plt.ylim((0,10))
for j in range(5):
    plt.subplot(2,5,j+6)
    for i in range(7):
        plt.plot(np.nanmean(np.nanmean(strialavg[thisfold][ontarget_ret_lax[thisfold]][:,:,j,i,1],0),0))
    plt.ylim((0,10))

snorm = scoa.get_norm_curves(oriavg,ontarget_ret_lax)

snorm.shape

plt.figure()
plt.plot(snorm.mean(0)[:,:,1].T)

plt.figure()
thisfold = keylist[1]
plt.plot(np.nanmean(oriavg[thisfold][ontarget_ret_lax[thisfold]],0)[:,:,1].T)

c = np.array((0,0.03,0.06,0.12,0.25,0.5,1))*100

oria

# +
plt.figure(figsize=(5,5))
thisfold = keylist[1]
for i in range(25,35):
    plt.subplot(5,5,i-25+1)
    plt.plot(s,oriavg[thisfold][ontarget_ret_lax[thisfold]][i,:,:,0])
    plt.axis('off')
    
plt.figure(figsize=(5,5))
thisfold = keylist[1]
for i in range(25,35):
    plt.subplot(5,5,i-25+1)
    plt.plot(s,oriavg[thisfold][ontarget_ret_lax[thisfold]][i,:,:,1])
    plt.axis('off')
# -

colors = plt.cm.viridis(np.linspace(0,1,7))

colors

gd = (5,11,12,13,24,28)
plt.figure(figsize=(10,4))
for i,g in enumerate(gd):
    plt.subplot(2,6,i+1)
    for j in range(7):
        plt.plot(s,oriavg[thisfold][ontarget_ret_lax[thisfold]][g,:,j,0],c=colors[j])
    plt.axis('off')
    plt.subplot(2,6,6+i+1)
    for j in range(7):
        plt.plot(s,oriavg[thisfold][ontarget_ret_lax[thisfold]][g,:,j,1],c=colors[j])
    plt.axis('off')
plt.subplot(2,6,7)
plt.axis('on')
plt.xlabel('size (deg.)')
plt.yticks([])
plt.savefig('vip_halo_size_contrast_pilot.pdf')

gd = (5,11,12,13,24,28)
plt.figure(figsize=(10,4))
for i,g in enumerate(gd):
    plt.subplot(2,6,i+1)
    for j in range(5):
        plt.plot(c,oriavg[thisfold][ontarget_ret_lax[thisfold]][g,j,:,0],c=colors[j])
    plt.axis('off')
    plt.subplot(2,6,6+i+1)
    for j in range(5):
        plt.plot(c,oriavg[thisfold][ontarget_ret_lax[thisfold]][g,j,:,1],c=colors[j])
    plt.axis('off')
plt.subplot(2,6,7)
plt.axis('on')
plt.xlabel('contrast (%)')
plt.yticks([])
plt.savefig('vip_halo_contrast_size_pilot.pdf')

plt.figure(figsize=(9,6))
thisfold = keylist[1]
X = np.nanmean(np.nanmean(np.nanmean(strialavg[thisfold][ontarget_ret_lax[thisfold]][:,:,:,:,:,5:10],0),0),-1)
for j in range(5):
    plt.subplot(1,2,1)
    plt.plot(c,X[j,:,0])
#     plt.ylim((0,10))
for j in range(5):
    plt.subplot(1,2,2)
    plt.plot(c,X[j,:,1])
#     plt.ylim((0,10))

np.logspace(np.log10(5),np.log10(36.5),5)


def norm_each(X):
    return X/np.nanmax(X.reshape(X.shape[0],-1),axis=1)[:,np.newaxis,np.newaxis,np.newaxis]


strialavg[thisfold][ontarget_ret_lax[thisfold]][:,:,:,:,:,5:10].shape

plt.figure(figsize=(9,6))
thisfold = keylist[1]
Y = np.nanmean(np.nanmean(strialavg[thisfold][ontarget_ret_lax[thisfold]][:,:,:,:,:,5:10],1),-1)
snorm = norm_each(Y)
X = snorm
x = np.nanmean(X,0)
s = np.logspace(np.log10(5),np.log10(36.5),5)
for i in range(7):
    plt.subplot(1,2,1)
    plt.plot(s,x[:,i,0])
    plt.ylim((0,1))
for i in range(7):
    plt.subplot(1,2,2)
    plt.plot(s,x[:,i,1])
    plt.ylim((0,1))

strialavg[thisfold][ontarget_ret_lax[thisfold]].shape

plt.figure(figsize=(9,6))
for j in range(5):
    plt.subplot(1,2,1)
    plt.plot(c,strialavg[thisfold][ontarget_ret_lax[thisfold]][:,:,j,:,0,5:10].mean(0).mean(0).mean(-1))
    plt.ylim((0,2))
for j in range(5):
    plt.subplot(1,2,2)
    plt.plot(c,strialavg[thisfold][ontarget_ret_lax[thisfold]][:,:,j,:,1,5:10].mean(0).mean(0).mean(-1))
    plt.ylim((0,2))

plt.figure(figsize=(9,6))
for j in range(5):
    plt.subplot(1,2,1)
    plt.plot(c,strialavg[thisfold][ontarget_ret_lax[thisfold]][:,:,j,:,0,5:10].mean(0).mean(0).mean(-1)-strialavg[thisfold][ontarget_ret_lax[thisfold]][:,:,j,:,1,5:10].mean(0).mean(0).mean(-1))
    plt.ylim((0,2))

x = strialavg[thisfold][ontarget_ret_lax[thisfold]][:,:,:,:,:,5:10].mean(0).mean(0).mean(-1)

plt.figure(figsize=(10,5))
plt.subplot(1,2,1)
for j in range(7):
    plt.scatter(x[:,j,0].flatten(),x[:,j,1].flatten())
plt.plot(x[:,:,0].flatten(),x[:,:,0].flatten(),'r')
plt.subplot(1,2,2)
for i in range(5):
    plt.scatter(x[i,:,0].flatten(),x[i,:,1].flatten())
plt.plot(x[:,:,0].flatten(),x[:,:,0].flatten(),'r')

plt.figure(figsize=(9,6))
for j in range(3):
    plt.subplot(1,2,1)
    plt.plot(np.nanmean(np.nanmean(np.nanmean(snorm[:,:,j,:,0,5:10],-1),0),0))
    plt.ylim((0,5))
for j in range(3):
    plt.subplot(1,2,2)
    plt.plot(snorm[:,:,j,:,1,5:10].mean(0).mean(0).mean(1))
    plt.ylim((0,5))

plt.figure()
ylim = (0,1.5)
plt.subplot(1,4,1)
plt.plot(snorm[0,:,:,0])
plt.ylim(ylim)
plt.subplot(1,4,2)
plt.plot(snorm[1,:,:,0])
plt.ylim(ylim)
plt.subplot(1,4,3)
plt.plot(snorm[0,:,:,1])
plt.ylim(ylim)
plt.subplot(1,4,4)
plt.plot(snorm[1,:,:,1])
plt.ylim(ylim)

plt.figure()
plt.plot(lb[thisfold][:,:,:,:,0].mean(0).mean(0).T)

ut.imshow_in_rows(soriavg[thisfold][ontarget_ret_lax[thisfold]][:,:,:,0])

frame_adjust = lambda x: np.hstack((x[0::4][:,np.newaxis],x[3::4][:,np.newaxis])).flatten()


frame_adjust(np.arange(100))


