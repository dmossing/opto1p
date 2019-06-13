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

import os
import scipy.io as sio
import matplotlib.pyplot as plt
# %matplotlib notebook
import autograd.numpy as np
from autograd import jacobian
import h5py
from oasis.functions import deconvolve
from oasis import oasisAR1, oasisAR2
import pyute as ut
from importlib import reload
reload(ut)
import matplotlib
matplotlib.rcParams['figure.facecolor'] = '1'
import scipy.stats as sst
from mpl_toolkits.mplot3d import Axes3D
import retinotopy_analysis as rt
reload(rt)
import scipy.optimize as sop
import pdb
import direct_fr_estimation_nd as dfre_nd
from direct_fr_estimation_nd import models_nd
reload(dfre_nd)
import pickle as pkl
from oasis.functions import deconvolve, estimate_parameters
import scipy.signal as ssi
import warnings

data_struct = sio.loadmat('vip_opto_data_struct.mat',squeeze_me=True)


def von_mises_function(x,kappa,mu,ht):
    return (np.exp(kappa*np.cos(x-mu)) + ht*np.exp(kappa*np.cos(x-mu-np.pi)))/(np.exp(kappa) + ht*np.exp(-kappa))


def naka_rushton_function(c,A,n,c50,top,bottom):
    return A*(c**n + top)/(c50**n + c**n + bottom)


def contrast_parametric_ori(ctuning,kappa,mu,ht,b): # cr the coupling coefficient to run speed
    x = data_obj.stim_id[2]*np.pi/4
    c = data_obj.stim_id[1]
    ori_multiplier = von_mises_function(x,kappa,mu,ht)
    contrast_multiplier = ctuning[c] - b # naka_rushton_function(c,A,n,c50,top,bottom)
    step_pulse = helper_vars['step_pulse'] # (T,), where T is the number of time points
    trialpart = contrast_multiplier*ori_multiplier
    return step_pulse[np.newaxis,:,np.newaxis]*trialpart[np.newaxis,np.newaxis,:] + b # (1,T,N)


def size_contrast_parametric_ori(sctuning,kappa,mu,ht,b): # cr the coupling coefficient to run speed
    x = data_obj.stim_id[2]*np.pi/4
    c = data_obj.stim_id[1]
    s = data_obj.stim_id[0]
    sctuning = sctuning.reshape((4,6))
    ori_multiplier = von_mises_function(x,kappa,mu,ht)
    size_contrast_multiplier = sctuning[s,c] - b # naka_rushton_function(c,A,n,c50,top,bottom)
    step_pulse = helper_vars['step_pulse'] # (T,), where T is the number of time points
    trialpart = size_contrast_multiplier*ori_multiplier
    return step_pulse[np.newaxis,:,np.newaxis]*trialpart[np.newaxis,np.newaxis,:] + b # (1,T,N)


def size_contrast_parametric_ori_varying_width(sctuning,kappa,mu,ht,b): # cr the coupling coefficient to run speed
    x = data_obj.stim_id[2]*np.pi/4
    c = data_obj.stim_id[1]
    s = data_obj.stim_id[0]
    sctuning = sctuning.reshape((4,6))
    ori_multiplier = von_mises_function(x,kappa[s],mu,ht)
    size_contrast_multiplier = sctuning[s,c] - b # naka_rushton_function(c,A,n,c50,top,bottom)
    step_pulse = helper_vars['step_pulse'] # (T,), where T is the number of time points
    trialpart = size_contrast_multiplier*ori_multiplier
    return step_pulse[np.newaxis,:,np.newaxis]*trialpart[np.newaxis,np.newaxis,:] + b # (1,T,N)


def size_contrast_parametric_ori_varying_width_lights(sctuning,kappa,mu,ht,b): # cr the coupling coefficient to run speed
    x = data_obj.stim_id[2]*np.pi/4
    c = data_obj.stim_id[1]
    s = data_obj.stim_id[0]
    lt = data_obj.stim_id[3]
    sctuning = sctuning.reshape((4,6,2))
    ori_multiplier = von_mises_function(x,kappa[s],mu,ht)
    size_contrast_multiplier = sctuning[s,c,lt] - b # naka_rushton_function(c,A,n,c50,top,bottom)
    step_pulse = helper_vars['step_pulse'] # (T,), where T is the number of time points
    trialpart = size_contrast_multiplier*ori_multiplier
    return step_pulse[np.newaxis,:,np.newaxis]*trialpart[np.newaxis,np.newaxis,:] + b # (1,T,N)


keylist = ['session_190301_M9835']
stimulus_id = data_struct[keylist[0]]['stimulus_id'][()]
nbefore = 4 #data_struct[keylist[0]]['nbefore'][()]
nafter = 4 # data_struct[keylist[0]]['nafter'][()]
running_speed_cm_s = data_struct[keylist[0]]['running_speed_cm_s'][()]
F = data_struct[keylist[0]]['F'][()]
stimlen = F.shape[2]-nbefore-nafter

# +
(ns,nc,nd,nl) = stimulus_id.max(1)+1
sizedict = {'ctuning':nc,'kappa':1,'mu':1,'ht':1,'b':1} #,'mu_nr':2,'sa_nr':1,'phi_nr':1,'sc_nr':1,'A_nr':1,'b_nr':1}
helper_vars = {}
helper_vars['step_pulse'] = np.concatenate((np.zeros((nbefore-1,)),np.ones((stimlen,)),np.zeros((nafter,))))
running_speed_cm_s[np.isnan(running_speed_cm_s)] = 0
helper_vars['runtrial'] = running_speed_cm_s < 7

eps = 5e-3
    
bdc = (eps,np.inf)
bdk = (0,4)
bdm = (-np.inf,np.inf)
bdht = (0,1)
bdb = (eps,np.inf)

bdd = {}

bdd['ctuning'] = [bdc]*nc
#bdd['sctuning'] = [bdc]*nc*ns
bdd['kappa'] = bdk
bdd['mu'] = bdm
bdd['ht'] = bdht
bdd['b'] = bdb
helper_vars['bounds_dict'] = bdd

fn_obj = dfre_nd.models_nd.model_from_rate_function(contrast_parametric_ori,sizedict,helper_vars)
# -

#tr = np.zeros((dtrialwise.shape[0],8))
imax = 50
tr = np.zeros((imax,ns,nc+4))
for isize in range(ns):
    runtrial = np.logical_and(np.logical_and(helper_vars['runtrial'],stimulus_id[0]==isize),stimulus_id[3]==0)
    for ind in range(imax):
        print(ind)
        np.seterr(all='print')
        data_obj = dfre_nd.data_obj(stimulus_id[:,runtrial],F[ind].T[:,runtrial],F[ind].T[:,runtrial],F[ind].flatten(),nbefore,nafter)
        fit_obj = dfre_nd.fit_obj(data_obj,pct_spike=97.5)
        np.seterr(all='warn')
        prestim = data_obj.F[:,:nbefore].mean()
        during = data_obj.F[:,nbefore:-nafter].mean()
        if prestim < during:
            guessA = 1
            guessb = eps
        else:
            guessA = eps
            guessb = 1
        tg = {}
        tg['ctuning'] = np.ones((nc,))*guessA
        tg['kappa'] = 1
        tg['mu'] = 0
        tg['ht'] = 1
        tg['b'] = guessb
        helper_vars['theta_guess_dict'] = tg
        fn_obj = dfre_nd.models_nd.model_from_rate_function(size_contrast_parametric_ori,sizedict,helper_vars) #
        warnings.filterwarnings('error')
        warnings.filterwarnings('ignore',category=DeprecationWarning)
        tr[ind,isize] = fit_obj.fit(fn_obj,factr=1e5,epsilon=1e-8)[0]

ut.plot_in_rows(tr[:50,2,:6])

ut.imshow_in_rows(tr[:50,:,:6])

plt.figure()
for ind in range(40,50):
    kappa = tr[ind,:,6]
    mu = tr[ind,:,7]
    ht = tr[ind,:,8]
    x = np.arange(32)*np.pi/16
    plt.subplot(10,1,ind-40+1)
    plt.imshow(von_mises_function(x[np.newaxis],kappa[:,np.newaxis],mu[:,np.newaxis],ht[:,np.newaxis]))

# +
(ns,nc,nd,nl) = stimulus_id.max(1)+1
#sizedict = {'sctuning':nc*ns,'kappa':1,'mu':1,'ht':1,'b':1} #,'mu_nr':2,'sa_nr':1,'phi_nr':1,'sc_nr':1,'A_nr':1,'b_nr':1}
sizedict = {'sctuning':nc*ns,'kappa':ns,'mu':1,'ht':1,'b':1} #,'mu_nr':2,'sa_nr':1,'phi_nr':1,'sc_nr':1,'A_nr':1,'b_nr':1}
helper_vars = {}
helper_vars['step_pulse'] = np.concatenate((np.zeros((nbefore-1,)),np.ones((stimlen,)),np.zeros((nafter,))))
running_speed_cm_s[np.isnan(running_speed_cm_s)] = 0
helper_vars['runtrial'] = running_speed_cm_s < 7

eps = 5e-3
    
bdc = (eps,np.inf)
bdk = (0,4)
bdm = (-np.inf,np.inf)
bdht = (0,1)
bdb = (eps,np.inf)

bdd = {}

bdd['sctuning'] = [bdc]*nc*ns
bdd['kappa'] = bdk
bdd['mu'] = bdm
bdd['ht'] = bdht
bdd['b'] = bdb
helper_vars['bounds_dict'] = bdd

fn_obj = dfre_nd.models_nd.model_from_rate_function(size_contrast_parametric_ori,sizedict,helper_vars)
# -

imax = 500 #F.shape[0]
tr = np.zeros((imax,2,ns*nc+4))
for ilight in range(2):
    runtrial = np.logical_and(helper_vars['runtrial'],stimulus_id[3]==ilight)
    for ind in range(imax):
        if np.nanmax(F[ind])>0:
            print(ind)
            np.seterr(all='print')
            data_obj = dfre_nd.data_obj(stimulus_id[:,runtrial],F[ind].T[:,runtrial],F[ind].T[:,runtrial],F[ind].flatten(),nbefore,nafter)
            fit_obj = dfre_nd.fit_obj(data_obj,pct_spike=97.5)
            np.seterr(all='warn')
            prestim = data_obj.F[:,:nbefore].mean()
            during = data_obj.F[:,nbefore:-nafter].mean()
            if prestim < during:
                guessA = 1
                guessb = eps
            else:
                guessA = eps
                guessb = 1
            tg = {}
            tg['sctuning'] = np.ones((ns*nc,))*guessA
            tg['kappa'] = 1
            tg['mu'] = 0
            tg['ht'] = 1
            tg['b'] = guessb
            helper_vars['theta_guess_dict'] = tg
            fn_obj = dfre_nd.models_nd.model_from_rate_function(size_contrast_parametric_ori,sizedict,helper_vars) #
            warnings.filterwarnings('error')
            warnings.filterwarnings('ignore',category=DeprecationWarning)
            tr[ind,ilight] = fit_obj.fit(fn_obj,factr=1e5,epsilon=1e-8)[0]

pdb.pm()

ut.imshow_in_rows(tr[50:100,0,:-4].reshape((-1,4,6)))
ut.imshow_in_rows(tr[50:100,1,:-4].reshape((-1,4,6)))

plt.figure()
plt.plot((TR[:imax,:,-1,-1]/TR[:imax,:,-1,-1].max(1)[:,np.newaxis]).T,c='b',alpha=0.02)

TR = tr[:,:,:-4].reshape((-1,2,4,6))
#ind = 1 #putative SST
ind = 31
plt.figure()
plt.subplot(1,2,1)
plt.imshow(TR[ind,0],vmax=TR[ind].max())
plt.subplot(1,2,2)
plt.imshow(TR[ind,1],vmax=TR[ind].max())

# +
(ns,nc,nd,nl) = stimulus_id.max(1)+1
#sizedict = {'sctuning':nc*ns,'kappa':1,'mu':1,'ht':1,'b':1} #,'mu_nr':2,'sa_nr':1,'phi_nr':1,'sc_nr':1,'A_nr':1,'b_nr':1}
sizedict = {'sctuning':nc*ns,'kappa':ns,'mu':1,'ht':1,'b':1} #,'mu_nr':2,'sa_nr':1,'phi_nr':1,'sc_nr':1,'A_nr':1,'b_nr':1}
helper_vars = {}
helper_vars['step_pulse'] = np.concatenate((np.zeros((nbefore-1,)),np.ones((stimlen,)),np.zeros((nafter,))))
running_speed_cm_s[np.isnan(running_speed_cm_s)] = 0
helper_vars['runtrial'] = running_speed_cm_s < 7

eps = 5e-3
    
bdc = (eps,np.inf)
bdk = (0,4)
bdm = (-np.inf,np.inf)
bdht = (0,1)
bdb = (eps,np.inf)

bdd = {}

bdd['sctuning'] = [bdc]*nc*ns
bdd['kappa'] = [bdk]*ns # bdk
bdd['mu'] = bdm
bdd['ht'] = bdht
bdd['b'] = bdb
helper_vars['bounds_dict'] = bdd

fn_obj = dfre_nd.models_nd.model_from_rate_function(size_contrast_parametric_ori_varying_width,sizedict,helper_vars)
# -

imax = F.shape[0]
tr_vw = np.zeros((imax,2,ns*nc+ns+3))
for ilight in range(2):
    runtrial = np.logical_and(helper_vars['runtrial'],stimulus_id[3]==ilight)
    for ind in range(imax):
        if np.nanmax(F[ind])>0:
            print(ind)
            np.seterr(all='print')
            data_obj = dfre_nd.data_obj(stimulus_id[:,runtrial],F[ind].T[:,runtrial],F[ind].T[:,runtrial],F[ind].flatten(),nbefore,nafter)
            fit_obj = dfre_nd.fit_obj(data_obj,pct_spike=97.5)
            np.seterr(all='warn')
            prestim = data_obj.F[:,:nbefore].mean()
            during = data_obj.F[:,nbefore:-nafter].mean()
            if prestim < during:
                guessA = 1
                guessb = eps
            else:
                guessA = eps
                guessb = 1
            tg = {}
            tg['sctuning'] = np.ones((ns*nc,))*guessA
            tg['kappa'] = np.ones((ns,)) # 1
            tg['mu'] = 0
            tg['ht'] = 1
            tg['b'] = guessb
            helper_vars['theta_guess_dict'] = tg
            fn_obj = dfre_nd.models_nd.model_from_rate_function(size_contrast_parametric_ori_varying_width,sizedict,helper_vars) #
            warnings.filterwarnings('error')
            warnings.filterwarnings('ignore',category=DeprecationWarning)
            tr_vw[ind,ilight] = fit_obj.fit(fn_obj,factr=1e5,epsilon=1e-8)[0]

# +
(ns,nc,nd,nl) = stimulus_id.max(1)+1
#sizedict = {'sctuning':nc*ns,'kappa':1,'mu':1,'ht':1,'b':1} #,'mu_nr':2,'sa_nr':1,'phi_nr':1,'sc_nr':1,'A_nr':1,'b_nr':1}
sizedict = {'sctuning':nc*ns*2,'kappa':ns,'mu':1,'ht':1,'b':1} #,'mu_nr':2,'sa_nr':1,'phi_nr':1,'sc_nr':1,'A_nr':1,'b_nr':1}
helper_vars = {}
helper_vars['step_pulse'] = np.concatenate((np.zeros((nbefore-1,)),np.ones((stimlen,)),np.zeros((nafter,))))
running_speed_cm_s[np.isnan(running_speed_cm_s)] = 0
helper_vars['runtrial'] = running_speed_cm_s < 7

eps = 5e-3
    
bdc = (eps,np.inf)
bdk = (0,4)
bdm = (-np.inf,np.inf)
bdht = (0,1)
bdb = (eps,np.inf)

bdd = {}

bdd['sctuning'] = [bdc]*sizedict['sctuning']
bdd['kappa'] = [bdk]*sizedict['kappa'] # bdk
bdd['mu'] = bdm
bdd['ht'] = bdht
bdd['b'] = bdb
helper_vars['bounds_dict'] = bdd

fn_obj = dfre_nd.models_nd.model_from_rate_function(size_contrast_parametric_ori_varying_width,sizedict,helper_vars)
# -

imax = F.shape[0]
tr_l = np.zeros((imax,2*ns*nc+ns+3))
for ilight in range(1):
    runtrial = helper_vars['runtrial']
    for ind in range(imax):
        if np.nanmax(F[ind])>0:
            print(ind)
            np.seterr(all='print')
            data_obj = dfre_nd.data_obj(stimulus_id[:,runtrial],F[ind].T[:,runtrial],F[ind].T[:,runtrial],F[ind].flatten(),nbefore,nafter)
            fit_obj = dfre_nd.fit_obj(data_obj,pct_spike=97.5)
            np.seterr(all='warn')
            prestim = data_obj.F[:,:nbefore].mean()
            during = data_obj.F[:,nbefore:-nafter].mean()
            if prestim < during:
                guessA = 1
                guessb = eps
            else:
                guessA = eps
                guessb = 1
            tg = {}
            tg['sctuning'] = np.ones((ns*nc*2,))*guessA
            tg['kappa'] = np.ones((ns,)) # 1
            tg['mu'] = 0
            tg['ht'] = 1
            tg['b'] = guessb
            helper_vars['theta_guess_dict'] = tg
            fn_obj = dfre_nd.models_nd.model_from_rate_function(size_contrast_parametric_ori_varying_width_lights,sizedict,helper_vars) #
            warnings.filterwarnings('error')
            warnings.filterwarnings('ignore',category=DeprecationWarning)
            tr_l[ind] = fit_obj.fit(fn_obj,factr=1e5,epsilon=1e-8)[0]

plt.figure()
for ind in range(40,50):
    kappa = tr_vw[ind,0,6:10]
    mu = tr[ind,0,10]
    ht = tr[ind,0,11]
    x = np.arange(32)*np.pi/16
    plt.subplot(10,1,ind-40+1)
    plt.imshow(von_mises_function(x[np.newaxis],kappa[:,np.newaxis],mu,ht))

kappa = tr_vw[:100,0,6:10]
fwhm = np.zeros_like(kappa)
fwhm[kappa>0.1] = np.arccos(1-np.log(2)/kappa[kappa>0.1])
plt.figure()
plt.plot(fwhm.T,c='b',alpha=0.1)

(1-np.log(2)/kappa[kappa>0.1]).min()

(1+np.log(2)/kappa[kappa>0.1])

dsr = sio.loadmat('vip_retinotopy_data_struct.mat',squeeze_me=True,struct_as_record=True)

keylist = [x for x in list(dsr.keys()) if not x[0]=='_']
k = 0

stimlen = dsr[keylist[k]]['F'][()].shape[-1]-nbefore-nafter

keylist


# +
def gaussian_rate_function(mu,sa,phi,sc,A,b): # cr the coupling coefficient to run speed
    xx = data_obj.stim_id #helper_vars['xx'] # (2,N), where N is the # of trials. Each is one of Ng grid locations
    step_pulse = helper_vars['step_pulse'] # (T,), where T is the number of time points
    trialpart = gaussian_trialpart_function(xx,mu,sa,phi,sc,A,b)
    #return (1+cr*sdxdt)*(step_pulse[np.newaxis,:,np.newaxis]*trialpart[np.newaxis,np.newaxis,:] + b) # (1,T,N)
    return step_pulse[np.newaxis,:,np.newaxis]*trialpart[np.newaxis,np.newaxis,:] + b # (1,T,N)

def gaussian_trialpart_function(xx,mu,sa,phi,sc,A,b):
    #siginv = np.array(((sa,sb),(sb,sc)))[:,:,0] # (2,2)
    rotmat = np.array(((np.cos(phi[0]),-np.sin(phi[0])),(np.sin(phi[0]),np.cos(phi[0])))) #[:,:,0]
    sigdiag = np.array(((sa[0],0),(0,sc[0]))) #[:,:,0] # (2,2)
    siginv = np.dot(np.dot(rotmat.T,sigdiag),rotmat)
    trialpart = (A-b)*np.exp(np.sum(-0.5*(xx-mu[:,np.newaxis])*np.dot(siginv,xx-mu[:,np.newaxis]),axis=0)) # (2,2) dot (2,N) -> (2,N). result is (N,)
    return trialpart


# +
sizedict = {'mu_r':2,'sa_r':1,'phi_r':1,'sc_r':1,'A_r':1,'b_r':1} #,'mu_nr':2,'sa_nr':1,'phi_nr':1,'sc_nr':1,'A_nr':1,'b_nr':1}
helper_vars = {}

locinds = dsr[keylist[k]]['stimulus_id'][()].T
strialwise = dsr[keylist[k]]['strialwise'][()]
nbefore = dsr[keylist[k]]['nbefore'][()]
nafter = dsr[keylist[k]]['nafter'][()]
F = dsr[keylist[k]]['F'][()]
stimlen = F.shape[-1]-nbefore-nafter
#strialwise = dsr[keylist[k]]['strialwise'][()]


(ny,nx) = locinds.max(1)
helper_vars['step_pulse'] = np.concatenate((np.zeros((nbefore-1,)),np.ones((stimlen,)),np.zeros((nbefore,))))
helper_vars['runtrial'] = dsr[keylist[k]]['running_speed_cm_s'][()] < 7

eps = 5e-3

guessy = (ny-1)/2
guessx = (nx-1)/2
guess_sa = 1/((ny-1)/2)**2
guess_phi = 0
guess_sc = 1/((nx-1)/2)**2
    
bdmux = (-0.5,float(nx)+0.5)
bdmuy = (-0.5,float(ny)+0.5)
bdsiginvii = (np.minimum(guess_sa,guess_sc),4) #(eps,np.inf)
bdphi = (-np.inf,np.inf)
bdA = (eps,np.inf)
bdb = (eps,np.inf)



bdd = {}
bdd['mu_r'] = [bdmuy,bdmux]
bdd['sa_r'] = bdsiginvii
bdd['phi_r'] = bdphi
bdd['sc_r'] = bdsiginvii
bdd['A_r'] = bdA
bdd['b_r'] = bdb
helper_vars['bounds_dict'] = bdd

fn_obj = dfre_nd.models_nd.model_from_rate_function(gaussian_rate_function,sizedict,helper_vars)
# -

reload(dfre_nd)

# imax = 100
tr_ret = np.zeros((F.shape[0],7))
runtrial = helper_vars['runtrial']
for ind in range(F.shape[0]):
    print(ind)
    if not np.nanmax(F[ind])==0:
        np.seterr(all='print')
        data_obj = dfre_nd.data_obj(locinds[:,runtrial],F[ind].T[:,runtrial],strialwise[ind].T[:,runtrial],F[ind].flatten(),nbefore,nafter)
        fit_obj = dfre_nd.fit_obj(data_obj,pct_spike=97.5)
        np.seterr(all='warn')
        prestim = data_obj.F[:,:nbefore].mean()
        during = data_obj.F[:,nbefore:-nafter].mean()
        if prestim < during:
            guessA = 1
            guessb = eps
        else:
            guessA = eps
            guessb = 1
        tg = {}
        tg['mu_r'] = np.array((guessy,guessx))
        tg['sa_r'] = guess_sa
        tg['phi_r'] = guess_phi
        tg['sc_r'] = guess_sc
        tg['A_r'] = guessA
        tg['b_r'] = guessb
        helper_vars['theta_guess_dict'] = tg
        fn_obj = dfre_nd.models_nd.model_from_rate_function(gaussian_rate_function,sizedict,helper_vars) #
        warnings.filterwarnings('error')
        warnings.filterwarnings('ignore',category=DeprecationWarning)
        try:
            tr_ret[ind] = fit_obj.fit(fn_obj,factr=1e5,epsilon=1e-8)[0]
        except:
            print(str(ind) + ' failed')
    else:
        print(str(ind) + ' is all zeros')


ret.shape

plt.figure()
plt.plot(tr_ret[:,0],alpha=0.2)
plt.plot(tr_ret[:,1],alpha=0.2)

plt.figure()
plt.hist(tr_ret[:,0],bins=100)

ontarget = np.logical_and(np.abs(tr_ret[:,0]-5)<2,np.abs(tr_ret[:,1]-5)<2)

ut.imshow_in_rows(tr_l[ontarget,:48:2].reshape((-1,4,6)))

plt.figure()
resp = tr_l[ontarget,:48].reshape((-1,4,6,2))
resp = resp/resp.max(-1).max(-1).max(-1)[:,np.newaxis,np.newaxis,np.newaxis]
mn = resp.mean(0)
plt.subplot(1,2,1)
plt.imshow(mn[:,:,0],vmin=mn.min(),vmax=mn.max())
plt.subplot(1,2,2)
plt.imshow(mn[:,:,1],vmin=mn.min(),vmax=mn.max())

sio.savemat('dfre_params')

ut.imshow_in_pairs(resp[:,:,:,0],resp[:,:,:,1])

pdb.pm()

plt.figure()
plt.plot(F[23].T,c='b',alpha=0.1)

# +
Ng = 9
imax = F.shape[0]
xx,yy = np.meshgrid(np.arange(Ng),np.arange(Ng))
xx = 1+np.concatenate((yy[np.newaxis],xx[np.newaxis]),axis=0).reshape((2,-1))
g = np.zeros((imax,Ng,Ng))

slicedict = {key:value for key,value in zip(fn_obj.paramlist,fn_obj.argslice)}

for i in range(imax):
    A = tr_ret[i][slicedict['A_r']]
    b = tr_ret[i][slicedict['b_r']]
    mu = tr_ret[i][slicedict['mu_r']]
    sa = tr_ret[i][slicedict['sa_r']]
    phi = tr_ret[i][slicedict['phi_r']]
    sc = tr_ret[i][slicedict['sc_r']]
    g[i] = gaussian_trialpart_function(xx,mu,sa,phi,sc,A,b).reshape((Ng,Ng))
# -

ontarget = np.logical_and(np.abs(tr_ret[:,0]-5)<0.5,np.abs(tr_ret[:,1]-5)<0.5)
ut.imshow_in_rows(g[ontarget])


