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
from astropy.convolution import Gaussian2DKernel
#import herschelhelp.cutouts_server as ctts
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
#from herschelhelp_internal.masterlist import merge_catalogues, nb_merge_dist_plot, specz_merge
from herschelhelp import image_plotting,utils
import pyvo as vo
import glob

import copy
#imported to be able to write tables while astropy is broken
from astropy.io import registry
from astropy.table.info import serialize_method_as

lofar_orig = Table.read('../data/data_release/final_component_catalogue-v1.0.fits')

fname = 'data/injected_source_radio_map.fits'
hdulist = fits.open(fname)
radim_header = hdulist[0].header
radim_wcs = wcs.WCS(radim_header).celestial
radim_data = hdulist[0].data #convert to mJy
hdulist.close()
radim_header['NAXIS']=2
radim_header['WCSAXES']=2

fname = '../data/data_release/radio_rms_image.fits'
hdulist = fits.open(fname)
radim_header_err = hdulist[0].header
radim_wcs_err = wcs.WCS(radim_header_err).celestial
radim_err = hdulist[0].data[0][0].astype('f4')*1000 #convert to mJy
hdulist.close()


'''for n in ['3','4']:
    for val in ['RPIX','DELT','UNIT','TYPE','RVAL']:
        try:
            radim_header.remove('C{}{}'.format(val,n))
        except:()
radim_header.remove('NAXIS3')
radim_header.remove('NAXIS4')
radim_header['NAXIS']=2
radim_header['WCSAXES']=2'''

#generate a psf for the first object in the catalogue 

n_size = 101
FWHM_factor = 2*np.sqrt(2*np.log(2))
sig_maj = 4/FWHM_factor#/2.355


prf = Gaussian2DKernel(sig_maj,x_size=n_size,y_size=n_size)
prf.normalize(mode='peak')

ras_old = lofar_orig['RA']
decs_old = lofar_orig['DEC']
ids_old = lofar_orig['Component_Name']

injected_cat = Table.read('data/injected_source_cat.fits')

ras_new = injected_cat['RA_injected']
decs_new = injected_cat['DEC_injected']
ids_new = np.arange(0,len(ras_new),1)

ras_all = np.append(ras_new,ras_old.data.data)
decs_all = np.append(decs_new,decs_old.data.data)
ids_all = np.append(ids_new,ids_old.data.data).astype(str)

n = (np.int(os.environ['SGE_TASK_ID'])-1)

ras = ras_new[n:n+1]
decs = decs_new[n:n+1]
ID = ids_new[n:n+1]

from astropy.coordinates import SkyCoord
from astropy import units as u
c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree)  
import pymoc
moc = pymoc.util.catalog.catalog_to_moc(c,60,15)


#create a cutout of the map around the centre cooridnate to speed up creation of the prior class
size = 120
cutout = Cutout2D(radim_data,c,size,wcs=radim_wcs,copy=True) 
im_data = cutout.data
im_wcs = cutout.wcs
im_header = im_wcs.to_header()
im_header['NAXIS'] = 2
im_header['NAXIS1'] = im_wcs.pixel_shape[0]
im_header['NAXIS2'] = im_wcs.pixel_shape[1]

cutout = Cutout2D(radim_err,c,size,wcs=radim_wcs,copy=True) 
im_err = cutout.data


prior=xidplus.prior(im_data,im_err,im_header,im_header, moc=moc)#Initialise with map, uncertianty map, wcs info and primary header

prior.prior_cat(ras_all,decs_all,'prior_cat',ID=ids_all,moc=moc)#Set input catalogue
prior.prior_bkg(-5.0,5)#Set prior on background (assumes Gaussian pdf with mu and sigma)

psf_dirty,coords,_ = pickle.load(open('data/PSF_coords.pkl','rb'))

prfs = []
pinds = []

for n in range(len(prior.sra)):
    if 'ILTJ' in prior.ID[n]:
        prfs.append(prf_true.array)
        
    else:
        prfs.append(psf_dirty[0])

    pinds.append(np.arange(0,prfs[n].shape[0],1))


prior.set_prfs(prfs,pinds,pinds)

prior.get_pointing_matrix_multiple_prf()
prior.upper_lim_map()

print('fitting '+ str(prior.nsrc)+' sources \n')
print('using ' +  str(prior.snpix) + ' pixels')

from xidplus.numpyro_fit import LOFAR150
fit=LOFAR150.LOFAR_150(prior,num_samples=1000,num_warmup=1000)
samples = fit.get_samples()['src_f']
posterior = xidplus.posterior_numpyro(fit,[prior])

if os.path.exists('data/xidplus_results/injected_sources/xidplus_run_{}'.format(n))==True:()
else:
    os.mkdir('data/xidplus_results/injected_sources/xidplus_run_{}'.format(n))
     
xidplus.save([prior],posterior,'data/xidplus_results/injected_sources/xidplus_run_{}/injected_sources_{}.pkl'.format(n,n))

LOFAR_cat=cat.create_LOFAR_cat(posterior,prior)
LOFAR_cat = Table.read(LOFAR_cat)
     
with serialize_method_as(LOFAR_cat, None):
            registry.write(LOFAR_cat,'data/xidplus_results/injected_sources/xidplus_run_{}/injected_sources_cat_{}.fits'.format(n,n),format='fits',overwrite=True)