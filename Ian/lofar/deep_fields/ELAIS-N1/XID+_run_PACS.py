from astropy.io import ascii, fits
import astropy
import pylab as plt
from astropy import wcs
from astropy.table import Table,Column,join,hstack
from astropy.coordinates import SkyCoord
from astropy import units as u
import pymoc
import glob
from time import sleep
import os


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
import pyvo as vo

#imported to be able to write tables while astropy is broken
from astropy.io import registry
from astropy.table.info import serialize_method_as



lofar = Table.read('data/data_release/final_cross_match_catalogue-v0.5.fits')
mask = (~np.isnan(lofar['F_SPIRE_250'])) | (~np.isnan(lofar['F_SPIRE_350'])) | (~np.isnan(lofar['F_SPIRE_500']))
lofar = lofar[~mask]

#remove if not running just the changed sources between v0.6 and v0.7
#lofar = Table.read('data/data_release/EN1_v0.6_v0.7_changedIDs.fits')

print(len(lofar))   


taskid = np.int(os.environ['SGE_TASK_ID'])-1
#taskid=3
print('taskid is: {}'.format(taskid))

batch_size = 20

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

ids = lofar['Source_Name'][ind_low:ind_up]

print('number of sources is: {}'.format(len(ras)))

imfolder='../../../../../HELP/dmu_products/dmu18/dmu18_ELAIS-N1/data/'

im100fits=imfolder + 'ELAIS-N1-100um-img_wgls.fits'#PACS 100 map
nim100fits=imfolder + 'ELAIS-N1-100um-img_noise.fits'#PACS 100 noise map
im160fits=imfolder + 'ELAIS-N1-160um-img_wgls.fits'#PACS 160 map
nim160fits=imfolder + 'ELAIS-N1-160um-img_noise.fits'#PACS 100 noise map

#-----100-------------
hdulist = fits.open(im100fits)
im100phdu=hdulist[0].header
im100hdu=hdulist[0].header
im100=hdulist[0].data
w_100 = wcs.WCS(hdulist[0].header)
pixsize100=3600.0*np.abs(hdulist[0].header['CDELT1']) #pixel size (in arcseconds)
hdulist.close()

hdulist = fits.open(nim100fits)
nim100=hdulist[0].data
hdulist.close()

#-----160-------------
hdulist = fits.open(im160fits)
im160phdu=hdulist[0].header
im160hdu=hdulist[0].header

im160=hdulist[0].data #convert to mJy
w_160 = wcs.WCS(hdulist[0].header)
pixsize160=3600.0*np.abs(hdulist[0].header['CDELT1']) #pixel size (in arcseconds)
hdulist.close()

hdulist = fits.open(nim160fits)
nim160=hdulist[0].data
hdulist.close()

prior_cat = Table.read('data/data_release/xidplus_prior_cat.fits')
#prior_cat = Table.read('data/data_release/xidplus_prior_cat_v0_7.fits')

from astropy.coordinates import SkyCoord
from astropy import units as u
c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree)  
import pymoc
moc=pymoc.util.catalog.catalog_to_moc(c,20,15)

#---prior100--------
prior100=xidplus.prior(im100,nim100,im100phdu,im100hdu, moc=moc)#Initialise with map, uncertianty map, wcs info and primary header
prior100.prior_cat(prior_cat['ra'],prior_cat['dec'],'prior_cat',ID=prior_cat['help_id'])#Set input catalogue
prior100.prior_bkg(0.0,5)#Set prior on background (assumes Gaussian pdf with mu and sigma)

#---prior160--------
prior160=xidplus.prior(im160,nim160,im160phdu,im160hdu, moc=moc)
prior160.prior_cat(prior_cat['ra'],prior_cat['dec'],'prior_cat',ID=prior_cat['help_id'])
prior160.prior_bkg(0.0,5)

pacs100_psf=fits.open(imfolder+'dmu18_PACS_100_PSF_ELAIS-N1_20170720.fits')
pacs160_psf=fits.open(imfolder+'dmu18_PACS_160_PSF_ELAIS-N1_20170720.fits')

centre100=np.long((pacs100_psf[1].header['NAXIS1']-1)/2)
radius100=15
centre160=np.long((pacs160_psf[1].header['NAXIS1']-1)/2)
radius160=25

pind100=np.arange(0,radius100+1+radius100,1)*3600*np.abs(pacs100_psf[1].header['CDELT1'])/pixsize100 #get 100 scale in terms of pixel scale of map
pind160=np.arange(0,radius160+1+radius160,1)*3600*np.abs(pacs160_psf[1].header['CDELT1'])/pixsize160 #get 160 scale in terms of pixel scale of map

prior100.set_prf(pacs100_psf[1].data[centre100-radius100:centre100+radius100+1,centre100-radius100:centre100+radius100+1]/1000.0,
                pind100,pind100)
prior160.set_prf(pacs160_psf[1].data[centre160-radius160:centre160+radius160+1,centre160-radius160:centre160+radius160+1]/1000.0,
                pind160,pind160)

prior100.get_pointing_matrix()
prior160.get_pointing_matrix()

from xidplus.stan_fit import PACS
fit=PACS.all_bands(prior100,prior160,iter=1000)

posterior=xidplus.posterior_stan(fit,[prior100,prior160])

priors = [prior100,prior160]

import xidplus.catalogue as cat
PACS_cat=cat.create_PACS_cat(posterior,priors[0],priors[1])
PACS_cat = Table.read(PACS_cat)

mask = [PACS_cat['help_id'][i] in ids for i in range(len(PACS_cat))]
PACS_cat = PACS_cat[mask]

if os.path.exists('data/fir/PACS/v2/xidplus_run_{}'.format(taskid))==True:()
else:
    os.mkdir('data/fir/PACS/v2/xidplus_run_{}'.format(taskid))

#the next couple of lines are an alternative way to save astropy table since the Table.write method is currently broken
xidplus.save([prior100,prior160],posterior,'data/fir/PACS/v2/xidplus_run_{}/lofar_xidplus_fir_{}_rerun.pkl'.format(taskid,taskid))
with serialize_method_as(PACS_cat, None):
            registry.write(PACS_cat,'data/fir/PACS/v2/xidplus_run_{}/lofar_xidplus_fir_{}_rerun.fits'.format(taskid,taskid),format='fits',overwrite=True)
#Table.write(PACS_cat,'data/fir/PACS/xidplus_run_{}/lofar_xidplus_fir_{}_rerun.fits'.format(taskid,taskid),format='fits',overwrite=True)

