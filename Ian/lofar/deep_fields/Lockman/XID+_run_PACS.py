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
from herschelhelp_internal.masterlist import merge_catalogues, nb_merge_dist_plot, specz_merge
import pyvo as vo


lofar = Table.read('data/data_release/final_cross_match_catalogue-v0.5.fits')
mask = (~np.isnan(lofar['F_SPIRE_250'])) | (~np.isnan(lofar['F_SPIRE_350'])) | (~np.isnan(lofar['F_SPIRE_500']))
lofar = lofar[~mask]

print(len(lofar))   


taskid = np.int(os.environ['SGE_TASK_ID'])-1
#taskid=3
print('taskid is: {}'.format(taskid))

'''dir_list = glob.glob('data/fir/*')
num_done = []
for folder in dir_list:
    num_done.append(int(folder.split('_')[-1]))
if taskid in num_done:
    sys.exit()'''

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



imfolder='../../../../../HELP/dmu_products/dmu18/dmu18_HELP-PACS-maps/data/'

im100fits=imfolder + 'Lockman-SWIRE_PACS100_v0.9.fits'#PACS 100 map
im160fits=imfolder + 'Lockman-SWIRE_PACS160_v0.9.fits'#PACS 160 map


#-----100-------------
hdulist = fits.open(im100fits)
im100phdu=hdulist[0].header
im100hdu=hdulist['image'].header
im100=hdulist['image'].data
w_100 = wcs.WCS(hdulist['image'].header)
pixsize100=3600.0*np.abs(hdulist['image'].header['CDELT1']) #pixel size (in arcseconds)

nim100=hdulist['error'].data
hdulist.close()

#-----160-------------
hdulist = fits.open(im160fits)
im160phdu=hdulist[0].header
im160hdu=hdulist['image'].header

im160=hdulist['image'].data #convert to mJy
w_160 = wcs.WCS(hdulist['image'].header)
pixsize160=3600.0*np.abs(hdulist['image'].header['CDELT1']) #pixel size (in arcseconds)

nim160=hdulist['error'].data
hdulist.close()

prior_cat = Table.read('data/data_release/xidplus_prior_cat_rerun_mips.fits')

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

pacs100_psf=fits.open('../../../../../HELP/dmu_products/dmu18/dmu18_Lockman-SWIRE/data/dmu18_PACS_100_PSF_Lockman-SWIRE_20171214.fits')
pacs160_psf=fits.open('../../../../../HELP/dmu_products/dmu18/dmu18_Lockman-SWIRE/data/dmu18_PACS_160_PSF_Lockman-SWIRE_20171214.fits')

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

prior100.upper_lim_map()
prior160.upper_lim_map()

from xidplus.stan_fit import PACS
fit=PACS.all_bands(prior100,prior160,iter=1000)

posterior=xidplus.posterior_stan(fit,[prior100,prior160])

priors = [prior100,prior160]

import xidplus.catalogue as cat
PACS_cat=cat.create_PACS_cat(posterior,priors[0],priors[1])
PACS_cat = Table.read(PACS_cat)

mask = [PACS_cat['help_id'][i] in ids for i in range(len(PACS_cat))]
PACS_cat = PACS_cat[mask]

if os.path.exists('data/fir/PACS/xidplus_run_{}'.format(taskid))==True:()
else:
    os.mkdir('data/fir/PACS/xidplus_run_{}'.format(taskid))
Table.write(PACS_cat,'data/fir/PACS/xidplus_run_{}/lofar_xidplus_fir_{}_rerun.fits'.format(taskid,taskid),overwrite=True)

xidplus.save([prior100,prior160],posterior,'data/fir/PACS/xidplus_run_{}/lofar_xidplus_fir_{}_rerun.pkl'.format(taskid,taskid))