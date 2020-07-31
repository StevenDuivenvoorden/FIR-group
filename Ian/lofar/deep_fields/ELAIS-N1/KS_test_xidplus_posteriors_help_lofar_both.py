from astropy.io import ascii, fits
import astropy
import pylab as plt

from astropy import wcs
from astropy.table import Table,Column,join,hstack,vstack
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.cosmology import WMAP9 as cosmo
import pymoc
import glob
from time import sleep
import os
from astropy.nddata.utils import Cutout2D
import herschelhelp.cutouts_server as ctts
import scipy

import numpy as np
import xidplus
from xidplus import moc_routines
import pickle
import xidplus.catalogue as cat

import sys
from herschelhelp_internal.utils import inMoc,flux_to_mag
from xidplus.stan_fit import SPIRE

import aplpy
import seaborn as sns
#sns.set(color_codes=True)
import pandas as pd
#sns.set_style("white")
import xidplus.posterior_maps as postmaps
import xidplus.catalogue as cat
from herschelhelp_internal.masterlist import merge_catalogues, nb_merge_dist_plot, specz_merge
from herschelhelp import image_plotting,utils
import pyvo as vo
import glob

import copy
from scipy.stats import ks_2samp



def find_pixels(ra_target,dec_target,radius,priors,posterior):
    
    rep_map_new = xidplus.posterior_maps.replicated_maps(priors,posterior)
    
    bval_map_new = xidplus.posterior_maps.make_Bayesian_pval_maps(priors[0],rep_map_new[0])
    pval_map_new = xidplus.posterior_maps.make_fits_image(priors[0],bval_map_new)
    
    wcs_new = wcs.WCS(pval_map_new[1].header)
    x = np.arange(0,pval_map_new[1].header['NAXIS1'],1)
    y = np.arange(0,pval_map_new[1].header['NAXIS2'],1)
    
    x,y = np.meshgrid(x,y)
    ras,decs = wcs_new.wcs_pix2world(x,y,0,ra_dec_order=True)
    
    pixels = []
    for n in range(len(ra_target)):
        dist = (ras-ra_target[n])**2 + (decs-dec_target[n])**2
        pixel_mask = dist<radius**2
        mask = ~np.isnan(pval_map_new[1].data)
        pixels_mask_flat = pixel_mask[mask].flatten()
        pixels.append(np.arange(0,np.sum(mask),1)[pixels_mask_flat])
    
    return(pixels)

def pval_map_sum(ra_target,dec_target,prior_new,prior_old,posterior_new,posterior_old,radius):
    
    #pixels_newx = prior_new.sx_pix
    #pixels_newy = prior_new.sy_pix
    
    #pixels_oldx = prior_old.sx_pix
    #pixels_oldy = prior_old.sy_pix
    
    #ids = priors_new.ID
    #mask = ids==source_id
    
    rep_map_new = xidplus.posterior_maps.replicated_maps(prior_new,posterior_new)
    bval_map_new = xidplus.posterior_maps.make_Bayesian_pval_maps(prior_new[0],rep_map_new[0])
    pval_map_new = xidplus.posterior_maps.make_fits_image(prior_new[0],bval_map_new)
    
    wcs_new = wcs.WCS(pval_map_new[1].header)
    x = np.arange(0,pval_map_new[1].header['NAXIS1'],1)
    y = np.arange(0,pval_map_new[1].header['NAXIS2'],1)

    x,y = np.meshgrid(x,y)
    ras,decs = wcs_new.wcs_pix2world(x,y,0,ra_dec_order=True)
    tot_pval_new = []
    for n in range(len(ra_target)):
        dist = (ras-ra_target[n])**2 + (decs-dec_target[n])**2
        mask = dist<radius**2
        print('number of pixel within radius is: {}'.format(np.sum(mask)))
        x,y = wcs_new.wcs_world2pix(ras[mask],decs[mask],0,ra_dec_order=True)
        x = x.astype(int)
        y = y.astype(int)
        tot_pval_new.append(np.sum(pval_map_new[1].data[y,x]))
    
    
    rep_map_old = xidplus.posterior_maps.replicated_maps(prior_old,posterior_old)
    bval_map_old = xidplus.posterior_maps.make_Bayesian_pval_maps(prior_old[0],rep_map_old[0])
    pval_map_old = xidplus.posterior_maps.make_fits_image(prior_old[0],bval_map_old)
    
    wcs_old = wcs.WCS(pval_map_old[1].header)
    x = np.arange(0,prior_old[0].imhdu['NAXIS1'],1)
    y = np.arange(0,prior_old[0].imhdu['NAXIS2'],1)
    x,y = np.meshgrid(x,y)
    ras,decs = wcs_old.wcs_pix2world(x,y,0,ra_dec_order=True)
    
    tot_pval_old = []
    for n in range(len(ra_target)):
        dist = (ras-ra_target[n])**2 + (decs-dec_target[n])**2
        mask = dist<radius**2
        print('number of pixel within radius is: {}'.format(np.sum(mask)))
        x,y = wcs_old.wcs_world2pix(ras[mask],decs[mask],0,ra_dec_order=True)
        x = x.astype(int)
        y = y.astype(int)
        tot_pval_old.append(np.sum(pval_map_old[1].data[y,x]))
    
    
    return(tot_pval_new,tot_pval_old)

def pval_summary(pixels,sigma_level,rep_map,prior):
    
    samples = rep_map[pixels,:]
    
    t = np.sum(((samples - prior.sim[pixels, None]) / (np.sqrt(2) * prior.snim[pixels, None])) ** 2.0, axis=0)
    ind_T = t / len(pixels) > sigma_level
    Bayes_pval_res_vals = ind_T.sum()/np.float(rep_map.shape[1])
    return(Bayes_pval_res_vals)

def get_xidplus_posterior(helpids,order,tiles,filepath=None,ras=[],decs=[]):
    # function to find the healpix for the given sources of the given order. You can provid this function a
    # list of help ids and it will return the healpix tile they are in. Alternatively you can provide a list
    # of ras and decs and it will do the same only faster. Finally if you give it the directory that the xidplus
    # posteriors are stored in then it will return you the path the posterior so it can be loaded in.
    
    # you need to provide this function a list of help ids, the order of the tiles desired and a list of tiles 
    # that you have posteriors for
    
    # if you haven't provided ras or decs it converts the help ids into ras and decs
    if (len(ras)==0) & (len(decs)==0):
        ras = []
        decs = []
        
        for helpid in helpids:
            ra,dec = utils.help_id_to_ra_dec(helpid)
            ras.append(ra)
            decs.append(dec)
        ras = np.array(ras)
        decs = np.array(decs)
        
    # finds the healpix tile of the given order that the coordinate is in
    hpxids = xidplus.moc_routines.get_HEALPix_pixels(order,ras,decs,unique=False)
    
    # tells you how many of the sources you gave it have matching XID+ posteriors from
    # the tile lits you gave
    in_tiles = np.array([healid in tiles for healid in hpxids])
    print('of the {} sources given, {} are within the tile list you have'.format(len(ras),np.sum(in_tiles)))
    
    # if a filepath was provided then here it gives you the filepath to each posterior and returns the 
    #hpxids and the filenames
    if filepath!=None:
        files = []
        for hpxid in hpxids:
            file = filepath + 'Tile_{}_{}.pkl'.format(hpxid,order)
            files.append(file)
        return(hpxids,files)
    
    # returns just the healpix ids if no filepath was given
    return(hpxids)

lofar_runs = glob.glob('data/fir/SPIRE/*/*.pkl')

lofar = Table.read('data/data_release/final_cross_match_catalogue-v0.5.fits')
mask = (~np.isnan(lofar['F_SPIRE_250'])) | (~np.isnan(lofar['F_SPIRE_350'])) | (~np.isnan(lofar['F_SPIRE_500']))
lofar = lofar[~mask]

taskid = np.int(os.environ['SGE_TASK_ID'])-1
batch_size = 200
if taskid*batch_size>len(lofar):
    print('Task id is too high. Trying to run code on more sources than exist')
    sys.exit()
ind_low = taskid*batch_size
if taskid*batch_size+batch_size>len(lofar):
    ind_up = len(lofar)
else:
    ind_up = taskid*batch_size+batch_size

ras = lofar['optRA'][ind_low:ind_up]
mask = np.isnan(ras)
ras[mask] = lofar['RA'][ind_low:ind_up][mask]

decs = lofar['optDec'][ind_low:ind_up]
mask = np.isnan(decs)
decs[mask] = lofar['DEC'][ind_low:ind_up][mask]

ids_centre = lofar['Source_Name'][ind_low:ind_up]

dir_list = glob.glob('data/fir/SPIRE/KS_results/KS_lofar_rerun_*.pkl')
num_done = []
for folder in dir_list:
    num_done.append(int(folder.split('_')[-1].replace('.pkl','')))
if taskid in num_done:
    sys.exit()

sources_done = []
both_result = []
lofar_result = []
help_result = []

for n,name in enumerate(ids_centre):
    
    if n%10==0:
        print(n)
        #print('loading new lofar and HELP posterior')
        lofar_file = 'data/fir/SPIRE_no_help/xidplus_run_{}/lofar_xidplus_fir_{}_rerun.pkl'.format(int(n/10),int(n/10))
        priors_lofar,posterior_lofar = xidplus.load(lofar_file)
        rep_map_lofar = postmaps.replicated_maps(priors_lofar,posterior_lofar)
        
        help_file = 'data/fir/SPIRE_no_lofar/xidplus_run_{}/lofar_xidplus_fir_{}_rerun.pkl'.format(int(n/10),int(n/10))
        priors_help,posterior_help = xidplus.load(help_file)
        rep_map_help = postmaps.replicated_maps(priors_help,posterior_help)
    if n%20==0:
        #print('loading new both posterior')
        both_file = 'data/fir/SPIRE/xidplus_run_{}/lofar_xidplus_fir_{}_rerun.pkl'.format(int(n/20),int(n/20))
        priors_both,posterior_both = xidplus.load(both_file)
        rep_map_both = postmaps.replicated_maps(priors_both,posterior_both)
        
    ra_target = ras[n]
    dec_target = decs[n]
    
    pixels_both = find_pixels([ra_target],[dec_target],18/3600,priors_both,posterior_both)
    pixels_lofar = find_pixels([ra_target],[dec_target],18/3600,priors_lofar,posterior_lofar)
    pixels_help = find_pixels([ra_target],[dec_target],18/3600,priors_help,posterior_help)
    
    for n,source in enumerate(pixels_lofar):
        sigmas = np.arange(0,3,0.1)
        y_both = []
        y_lofar = []
        y_help = []
        for sigma in sigmas:
            y_both.append(pval_summary(pixels_both,sigma,rep_map_both[0],priors_both[0]))
            y_lofar.append(pval_summary(pixels_lofar,sigma,rep_map_lofar[0],priors_lofar[0]))
            y_help.append(pval_summary(pixels_help,sigma,rep_map_help[0],priors_help[0]))
    
        both_result.append(y_both)
        lofar_result.append(y_lofar)
        help_result.append(y_help)
    
        plt.plot(sigmas,y_both,label='both')
        plt.plot(sigmas,y_lofar,label='lofar')
        plt.plot(sigmas,y_help,label='help')
        plt.legend()
        plt.show()

        sources_done.append([name,both_file,lofar_file,help_file])

ks_test = []
sign = []

for n in range(len(help_result)):

    ks_test.append(ks_2samp(np.array(both_result[n]),np.array(help_result[n])))
            
    diff = np.array(both_result[n]) - np.array(help_result[n])
    abs_diff = abs(diff)
    sign.append(diff[np.argmax(abs_diff)])
sign = np.array(sign)
mask = sign<0
sign[mask] = -1

mask = sign>0
sign[mask] = 1
            
f = open('data/fir/SPIRE/KS_results/KS_lofar_rerun_{}.pkl'.format(taskid),'wb')
pickle.dump([ks_test,sign,sources_done,both_result,lofar_result,help_result],f)