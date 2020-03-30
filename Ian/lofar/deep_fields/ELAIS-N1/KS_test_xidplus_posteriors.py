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

batch_size = 20


sources_done = []
ks_test = []
sign = []
for file in lofar_runs: 
    print(file)
    priors,posterior = xidplus.load(file)
    rep_map_lofar = postmaps.replicated_maps(priors,posterior)
    
    filename = file

    taskid = int(filename.split('_')[-2])
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

    #SPIRE_cat = cat.create_SPIRE_cat(posterior,priors[0],priors[1],priors[2])
    #SPIRE_cat = Table.read(SPIRE_cat)
    #mask = ['ILTJ' in SPIRE_cat['HELP_ID'][i] for i in range(len(SPIRE_cat))]

    #ids_lofar = SPIRE_cat[mask]['HELP_ID']
    #ras_lofar = SPIRE_cat[mask]['RA']
    #decs_lofar = SPIRE_cat[mask]['Dec']

    filenames = glob.glob('/lustre/scratch/astro/pdh21/ELAIS_N1/SPIRE/output/Tile_*_9.pkl')
    tiles = [int(name.split('/')[-1].split('_')[1]) for name in filenames]
    healpixids,filepaths = get_xidplus_posterior(ids_centre,9,tiles,ras=ras,decs=decs,filepath='/lustre/scratch/astro/pdh21/ELAIS_N1/SPIRE/output/')

    for n,help_run in enumerate(np.unique(filepaths)):
        #print(help_run)
        ra_target = ras[np.where(np.array(filepaths)==help_run)]
        dec_target = decs[np.where(np.array(filepaths)==help_run)]
        ids_target = ids_centre[np.where(np.array(filepaths)==help_run)]
        
        priors_help,posterior_help = xidplus.load(help_run)
        rep_map_help = postmaps.replicated_maps(priors_help,posterior_help)
    
        pixels_lofar = find_pixels(ra_target,dec_target,18/3600,priors,posterior)
        pixels_help = find_pixels(ra_target,dec_target,18/3600,priors_help,posterior_help)
        
        #print(pixels_lofar)

        for n,source in enumerate(pixels_lofar):
            sigmas = np.arange(0,3,0.1)
            y_lofar = []
            y_help = []
            for sigma in sigmas:
                y_lofar.append(pval_summary(pixels_lofar[n],sigma,rep_map_lofar[0],priors[0]))
                y_help.append(pval_summary(pixels_help[n],sigma,rep_map_help[0],priors_help[0]))
    
            #print(ks_2samp(np.array(y_lofar),np.array(y_help)))
            ks_test.append(ks_2samp(np.array(y_lofar),np.array(y_help)))
            
            diff = np.array(y_lofar) - np.array(y_help)
            abs_diff = abs(diff)
            sign.append(diff[np.argmax(abs_diff)])

            sources_done.append([ids_target[n],file,help_run])
            
f = open('data/fir/SPIRE/KS_lofar_rerun.pkl','wb')
pickle.dump([ks_test,sign,sources_done],f)