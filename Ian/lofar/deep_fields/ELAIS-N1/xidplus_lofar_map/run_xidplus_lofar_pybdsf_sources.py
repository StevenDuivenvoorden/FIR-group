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

lofar_orig = Table.read('../data/data_release/final_cross_match_catalogue-v1.0.fits')

fname = '../data/image_full_ampphase_di_m.NS.psf.fits'
hdulist = fits.open(fname)
radim_header_psf = hdulist[0].header
radim_wcs_psf = wcs.WCS(radim_header_psf).celestial
radim_data_psf = hdulist[0].data[0][0]
radim_pixsize = abs(radim_wcs_psf.wcs.cdelt[0])*3600
hdulist.close()

fname = '../data/data_release/radio_image.fits'
hdulist = fits.open(fname)
radim_header = hdulist[0].header
radim_wcs = wcs.WCS(radim_header).celestial
radim_data = hdulist[0].data[0][0].astype('f4')*1000 #convert to mJy
hdulist.close()
radim_header['NAXIS']=2
radim_header['WCSAXES']=2

fname = '../data/data_release/radio_rms_image.fits'
hdulist = fits.open(fname)
radim_header_err = hdulist[0].header
radim_wcs_err = wcs.WCS(radim_header_err).celestial
radim_err = hdulist[0].data[0][0].astype('f4')*1000 #convert to mJy
hdulist.close()

fname = '../data/data_release/image_full_ampphase_di_m.NS.mask01-2.fits'
hdulist = fits.open(fname)
radim_header_mask = hdulist[0].header
mask_deconv_wcs = wcs.WCS(radim_header_mask).celestial
mask_deconv = hdulist[0].data[0][0]
hdulist.close()

for n in ['3','4']:
    for val in ['RPIX','DELT','UNIT','TYPE','RVAL']:
        try:
            radim_header.remove('C{}{}'.format(val,n))
        except:()
radim_header.remove('NAXIS3')
radim_header.remove('NAXIS4')
radim_header['NAXIS']=2
radim_header['WCSAXES']=2

#generate a psf for the first object in the catalogue 

n_size = 101
FWHM_factor = 2*np.sqrt(2*np.log(2))
sig_maj = 4/FWHM_factor#/2.355


prf = Gaussian2DKernel(sig_maj,x_size=n_size,y_size=n_size)
prf.normalize(mode='peak')



n = (np.int(os.environ['SGE_TASK_ID'])-1)*10

ras = lofar_orig['RA'][n:n+1]
decs = lofar_orig['DEC'][n:n+1]
ID = lofar_orig['Source_Name'][n:n+1]

from astropy.coordinates import SkyCoord
from astropy import units as u
c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree)  
import pymoc
moc_near=pymoc.util.catalog.catalog_to_moc(c,60,15)

moc_mask = inMoc(lofar_orig['RA'],lofar_orig['DEC'],moc_near)
moc = pymoc.util.catalog.catalog_to_moc(c,120,15)


prior=xidplus.prior(radim_data,radim_err,radim_header,radim_header, moc=moc)#Initialise with map, uncertianty map, wcs info and primary header
prior.prior_cat(lofar_orig['RA'][moc_mask],lofar_orig['DEC'][moc_mask],'prior_cat',ID=lofar_orig['Source_Name'][moc_mask],moc=moc)#Set input catalogue
prior.prior_bkg(-5.0,5)#Set prior on background (assumes Gaussian pdf with mu and sigma)

pind = np.arange(0,n_size)#*(1.0/radim_pixsize)
prior.set_prf(prf.array,pind,pind)

prior.get_pointing_matrix()
prior.upper_lim_map()

print('fitting '+ str(prior.nsrc)+' sources \n')
print('using ' +  str(prior.snpix) + ' pixels')

from xidplus.numpyro_fit import LOFAR150
fit=LOFAR150.LOFAR_150(prior,num_samples=1000,num_warmup=1000)
samples = fit.get_samples()['src_f']
posterior = xidplus.posterior_numpyro(fit,[prior])

if os.path.exists('data/xidplus_results/pybdsf/xidplus_run_{}'.format(n))==True:()
else:
    os.mkdir('data/xidplus_results/pybdsf/xidplus_run_{}'.format(n))
     
xidplus.save([prior],posterior,'data/xidplus_results/pybdsf/xidplus_run_{}/lofar_xidplus_fir_{}_rerun.pkl'.format(n,n))

LOFAR_cat=cat.create_LOFAR_cat(posterior,prior,0)
LOFAR_cat = Table.read(LOFAR_cat)
     
with serialize_method_as(LOFAR_cat, None):
            registry.write(LOFAR_cat,'data/xidplus_results/pybdsf/xidplus_run_{}/lofar_xidplus_fir_{}_rerun.fits'.format(n,n),format='fits',overwrite=True)