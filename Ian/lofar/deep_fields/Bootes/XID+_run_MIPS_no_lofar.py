from astropy.io import ascii, fits
import astropy
import pylab as plt
from astropy import wcs
from astropy.table import Table,Column,join,hstack
from astropy.coordinates import SkyCoord
from astropy import units as u
import pymoc
import glob
import time
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

from astropy.io import registry
from astropy.table.info import serialize_method_as

lofar = Table.read('data/data_release/final_cross_match_catalogue-v0.5.fits')
mask = (~np.isnan(lofar['F_MIPS_24'])) 
lofar = lofar[~mask]

#remove if not running just the changed sources between v0.6 and v0.7
#lofar = Table.read('data/data_release/Bootes_v0.6_v0.7_changedIDs.fits')
                        
taskid = np.int(os.environ['SGE_TASK_ID'])-1
print('taskid is: {}'.format(taskid))

'''dir_list = glob.glob('data/fir/MIPS/*')
num_done = []
for folder in dir_list:
    num_done.append(int(folder.split('_')[-1]))
if taskid in num_done:
    sys.exit()'''

batch_size = 50
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

print('reading in prior cat')
prior_cat = Table.read('data/data_release/xidplus_prior_cat_MIPS_rerun.fits')
mask = np.array(['ILTJ' in name for name in prior_cat['help_id']])
prior_cat = prior_cat[~mask]
#prior_cat = Table.read('data/data_release/xidplus_prior_cat_MIPS_v0_7.fits')
MIPS_lower = prior_cat['MIPS_lower']
MIPS_upper = prior_cat['MIPS_upper']

from astropy.coordinates import SkyCoord
from astropy import units as u
c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree)  
import pymoc
moc=pymoc.util.catalog.catalog_to_moc(c,30,15)


#Read in the herschel images
imfolder='../../../../../HELP/dmu_products/dmu17/dmu17_HELP_Legacy_maps/Bootes/data/'

pswfits=imfolder+'wp4_bootes_mips24_map_v1.0.fits.gz'#SPIRE 250 map

print('reading in MIPS map')
MIPS_Map = fits.open(pswfits)

#-----250-------------
hdulist = fits.open(pswfits)
im250phdu=hdulist[0].header
im250hdu=hdulist[1].header

im250=hdulist[1].data #convert to mJy
nim250=hdulist[2].data #convert to mJy
w_250 = wcs.WCS(hdulist[1].header)
hdulist.close()

MIPS_psf=fits.open(imfolder+'dmu17_MIPS_PSF_Bootes_20180129.fits')


#---prior250--------
print('creating prior')
prior250=xidplus.prior(im250,nim250,im250phdu,im250hdu, moc=moc)#Initialise with map, uncertianty map, wcs info and primary header
prior250.prior_cat(prior_cat['ra'],prior_cat['dec'],'prior_cat',flux_lower=MIPS_lower,flux_upper=MIPS_upper,ID=prior_cat['help_id'])#Set input catalogue
prior250.prior_bkg(-5.0,5)#Set prior on background (assumes Gaussian pdf with mu and sigma)

centre=np.long((MIPS_psf[1].header['NAXIS1']-1)/2)
radius=20
prior250.set_prf(MIPS_psf[1].data[centre-radius:centre+radius+1,centre-radius:centre+radius+1]/1.0E6,np.arange(0,41/2.0,0.5),np.arange(0,41/2.0,0.5))#requires PRF as 2d grid, and x and y bins for grid (in pixel scale)

prior250.get_pointing_matrix()


from xidplus.stan_fit import MIPS,SPIRE
fit=MIPS.MIPS_24(prior250,iter=1000)

posterior=xidplus.posterior_stan(fit,[prior250])
priors = [prior250]

import xidplus.catalogue as cat
MIPS_cat=cat.create_MIPS_cat(posterior,priors[0],0)
MIPS_cat = Table.read(MIPS_cat)

mask = [MIPS_cat['help_id'][i] in ids for i in range(len(MIPS_cat))]
MIPS_cat = MIPS_cat[mask]

if os.path.exists('data/fir/MIPS_no_lofar/xidplus_run_{}'.format(taskid))==True:()
else:
    os.mkdir('data/fir/MIPS_no_lofar/xidplus_run_{}'.format(taskid))
    
xidplus.save([prior250],posterior,'data/fir/MIPS_no_lofar/xidplus_run_{}/lofar_xidplus_fir_{}'.format(taskid,taskid))

#the next couple of lines are an alternative way to save astropy table since the Table.write method is currently broken
with serialize_method_as(MIPS_cat, None):
            registry.write(MIPS_cat,'data/fir/MIPS_no_lofar/xidplus_run_{}/lofar_xidplus_fir_{}.fits'.format(taskid,taskid),format='fits',overwrite=True)    
#Table.write(MIPS_cat,'data/fir/MIPS/xidplus_run_{}/lofar_xidplus_fir_{}_rerun.fits'.format(taskid,taskid),overwrite=True)


