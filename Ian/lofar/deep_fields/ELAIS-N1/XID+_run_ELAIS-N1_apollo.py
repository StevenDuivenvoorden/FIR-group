print('started running python script')
import astropy
from astropy.io import ascii, fits
import pylab as plt
from astropy import wcs
from astropy.table import Table,Column,join,hstack
from astropy.coordinates import SkyCoord
from astropy import units as u
import pymoc
import glob
from time import sleep
import os
import sys

print('importing modules.....half way')

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

print('finished importing modules')


taskid = np.int(os.environ['SGE_TASK_ID'])-1
#taskid=3
print('taskid is: {}'.format(taskid))

#only here for rerunning EN1 with the updated catalogue
lofar01 = Table.read('data/data_release/final_cross_match_catalogue-v0.1.fits')
lofar05 = Table.read('data/data_release/final_cross_match_catalogue-v0.5.fits')
new_x_old = join(lofar05,lofar01,join_type='left',keys='Source_Name')

mask = new_x_old['optRA_1']==new_x_old['optRA_2']
lofar = lofar05[~mask]

#Read in the LOFAR data both optical and radio
#lofar = Table.read('data/data_release/final_cross_match_catalogue-v0.1.fits')
mask = (~np.isnan(lofar['F_SPIRE_250'])) | (~np.isnan(lofar['F_SPIRE_350'])) | (~np.isnan(lofar['F_SPIRE_500']))
lofar = lofar[~mask]

print(len(lofar))

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

#print(len(ras))

'''columns = 'ra','dec','help_id','flag_optnir_det','f_mips_24'
masterlist = Table.read('../../../../../HELP/dmu_products/dmu32/dmu32_ELAIS-N1/data/ELAIS-N1_20171020.fits')
help_masterlist = masterlist[columns]
masterlist = 0

lofar_coords = SkyCoord(ra,dec,unit='deg')
help_coords = SkyCoord(help_masterlist['ra'],help_masterlist['dec'],unit='deg')
radius = 2
idx_help, idx_lofar, d2d, d3d = lofar_coords.search_around_sky(
    help_coords, radius*u.arcsec)

prior_mask = (help_masterlist['flag_optnir_det']>=5) & (help_masterlist['f_mips_24']>20)
prior_cat = help_masterlist[prior_mask]
for n,pos in enumerate(ras):
    prior_cat.add_row([ras[n],decs[n],ids[n],-99,np.nan])'''

prior_cat = Table.read('data/data_release/xidplus_prior_cat_rerun.fits')
#xid_rerun = Column(name='XID_rerun',data=np.zeros(len(prior_cat))-99)
#prior_cat.add_column(xid_rerun)

'''
#if there is only one crossmatch within the search radius then match them if the source is in the prior list
#otherwise add the ra and dec of the lofar optical counterpart to the prior list
XID_rerun = []
source_type = []
mask = (help_masterlist[idx_help]['flag_optnir_det']>=5) & (help_masterlist[idx_help]['f_mips_24']>20)
idx_true = idx_help[mask]

if len(d2d) == 0:
    min_d2d = np.nan
else:
    min_d2d = np.min(d2d)

if d2d[mask] != min_d2d:
    source_type.append('nearer_non_prior')
    XID_rerun.append(False)
if np.sum(mask)==0:
    prior_cat.add_row([ra,dec,'lofar',np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,-99,np.nan,False])
    source_type.append('radio_position')
    XID_rerun.append(True)
if np.sum(mask)>1:
    nearer = ''
    if np.min(d2d[mask]) != min_d2d:
        nearer = 'nearer_non_prior'
    for n in range(np.sum(mask)):
        source_type.append('multiple_prior'+nearer)
        XID_rerun.append(False)
if np.sum(mask)==1:
    source_type.append('nearest_prior')
    XID_rerun.append(False)
    if d2d[mask] != np.min(d2d):
        source_type.append('nearer_non_prior')
        XID_rerun.append(False)

print(XID_rerun)

if False in XID_rerun:
    mask = lofar['ra']==ra
    lofar_fir = hstack([lofar[mask],help_masterlist[idx_true]])
    XID_rerun_col = Column(data=XID_rerun,name='XID_rerun',dtype=bool)
    source_type_col = Column(data=source_type,name='source_type',dtype=str)
    lofar_fir.add_columns([XID_rerun_col,source_type_col])
    if os.path.exists('data/fir_v2/xidplus_run_{}'.format(taskid)):()
    else:
        os.mkdir('data/fir_v2/xidplus_run_{}'.format(taskid))
    Table.write(lofar_fir,'data/fir_v2/xidplus_run_{}/lofar_xidplus_fir_{}.fits'.format(taskid,taskid),overwrite=True)
'''


#Read in the herschel images
imfolder='../../../../../HELP/dmu_products/dmu19/dmu19_HELP-SPIRE-maps/data/'

pswfits=imfolder+'ELAIS-N1_SPIRE250_v1.0.fits'#SPIRE 250 map
pmwfits=imfolder+'ELAIS-N1_SPIRE350_v1.0.fits'#SPIRE 350 map
plwfits=imfolder+'ELAIS-N1_SPIRE500_v1.0.fits'#SPIRE 500 map

#-----250-------------
hdulist = fits.open(pswfits)
im250phdu=hdulist[0].header
im250hdu=hdulist['image'].header

im250=hdulist['image'].data*1.0E3 #convert to mJy
nim250=hdulist['error'].data*1.0E3 #convert to mJy
w_250 = wcs.WCS(hdulist['image'].header)
pixsize250=3600.0*w_250.wcs.cd[1,1] #pixel size (in arcseconds)
hdulist.close()
#-----350-------------
hdulist = fits.open(pmwfits)
im350phdu=hdulist[0].header
im350hdu=hdulist['image'].header

im350=hdulist['image'].data*1.0E3 #convert to mJy
nim350=hdulist['error'].data*1.0E3 #convert to mJy
w_350 = wcs.WCS(hdulist['image'].header)
pixsize350=3600.0*w_350.wcs.cd[1,1] #pixel size (in arcseconds)
hdulist.close()
#-----500-------------
hdulist = fits.open(plwfits)
im500phdu=hdulist[0].header
im500hdu=hdulist['image'].header 
im500=hdulist['image'].data*1.0E3 #convert to mJy
nim500=hdulist['error'].data*1.0E3 #convert to mJy
w_500 = wcs.WCS(hdulist['image'].header)
pixsize500=3600.0*w_500.wcs.cd[1,1] #pixel size (in arcseconds)
hdulist.close()

from astropy.coordinates import SkyCoord
from astropy import units as u
c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree)  
import pymoc
moc=pymoc.util.catalog.catalog_to_moc(c,60,15)

#---prior250--------
prior250=xidplus.prior(im250,nim250,im250phdu,im250hdu, moc=moc)#Initialise with map, uncertianty map, wcs info and primary header
prior250.prior_cat(prior_cat['ra'],prior_cat['dec'],'prior_cat',ID=prior_cat['help_id'])#Set input catalogue
prior250.prior_bkg(-5.0,5)#Set prior on background (assumes Gaussian pdf with mu and sigma)
#---prior350--------
prior350=xidplus.prior(im350,nim350,im350phdu,im350hdu, moc=moc)
prior350.prior_cat(prior_cat['ra'],prior_cat['dec'],'prior_cat',ID=prior_cat['help_id'])
prior350.prior_bkg(-5.0,5)

#---prior500--------
prior500=xidplus.prior(im500,nim500,im500phdu,im500hdu, moc=moc)
prior500.prior_cat(prior_cat['ra'],prior_cat['dec'],'prior_cat',ID=prior_cat['help_id'])
prior500.prior_bkg(-5.0,5)

#pixsize array (size of pixels in arcseconds)
pixsize=np.array([pixsize250,pixsize350,pixsize500])
#point response function for the three bands
prfsize=np.array([18.15,25.15,36.3])
#use Gaussian2DKernel to create prf (requires stddev rather than fwhm hence pfwhm/2.355)
from astropy.convolution import Gaussian2DKernel

##---------fit using Gaussian beam-----------------------
prf250=Gaussian2DKernel(prfsize[0]/2.355,x_size=101,y_size=101)
prf250.normalize(mode='peak')
prf350=Gaussian2DKernel(prfsize[1]/2.355,x_size=101,y_size=101)
prf350.normalize(mode='peak')
prf500=Gaussian2DKernel(prfsize[2]/2.355,x_size=101,y_size=101)
prf500.normalize(mode='peak')

pind250=np.arange(0,101,1)*1.0/pixsize[0] #get 250 scale in terms of pixel scale of map
pind350=np.arange(0,101,1)*1.0/pixsize[1] #get 350 scale in terms of pixel scale of map
pind500=np.arange(0,101,1)*1.0/pixsize[2] #get 500 scale in terms of pixel scale of map

prior250.set_prf(prf250.array,pind250,pind250)#requires PRF as 2d grid, and x and y bins for grid (in pixel scale)
prior350.set_prf(prf350.array,pind350,pind350)
prior500.set_prf(prf500.array,pind500,pind500)

prior250.get_pointing_matrix()
prior350.get_pointing_matrix()
prior500.get_pointing_matrix()

prior250.upper_lim_map()
prior350.upper_lim_map()
prior500.upper_lim_map()

print('fitting '+ str(prior250.nsrc)+' sources \n')
print('using ' +  str(prior250.snpix)+', '+ str(prior350.snpix)+' and '+ str(prior500.snpix)+' pixels')

from xidplus.stan_fit import SPIRE
fit=SPIRE.all_bands(prior250,prior350,prior500,iter=1000)

posterior=xidplus.posterior_stan(fit,[prior250,prior350,prior500])
priors = [prior250,prior350,prior500]
import xidplus.catalogue as cat
SPIRE_cat=cat.create_SPIRE_cat(posterior,priors[0],priors[1],priors[2])
SPIRE_cat = Table.read(SPIRE_cat)

#mask = ['ILTJ' in SPIRE_cat['HELP_ID'][i] for i in range(len(SPIRE_cat)) if SPIRE_cat['RA'][i] in ras]
mask = [SPIRE_cat['HELP_ID'][i] in ids for i in range(len(SPIRE_cat))]
SPIRE_cat = SPIRE_cat[mask]

'''mask = SPIRE_cat['HELP_ID']=='lofar'
mask_lofar = lofar['ra']==ra
mask_pcat = prior_cat['help_id']=='lofar'
SPIRE_cat = SPIRE_cat[mask]
SPIRE_cat.add_columns([prior_cat['f_mips_24'][mask_pcat],prior_cat['flag_optnir_det'][mask_pcat]])
lofar_fir = hstack([lofar[mask_lofar],SPIRE_cat])
XID_rerun_col = Column(data=XID_rerun,name='XID_rerun',dtype=bool)
source_type_col = Column(data=source_type,name='source_type',dtype=str)
lofar_fir.add_columns([XID_rerun_col,source_type_col])'''
   
if os.path.exists('data/fir_v10/xidplus_run_{}'.format(taskid))==True:()
else:
    os.mkdir('data/fir_v10/xidplus_run_{}'.format(taskid))
Table.write(SPIRE_cat,'data/fir_v10/xidplus_run_{}/lofar_xidplus_fir_{}_rerun.fits'.format(taskid,taskid),overwrite=True)

xidplus.save([prior250,prior350,prior500],posterior,'data/fir_v10/xidplus_run_{}/lofar_xidplus_fir_{}_rerun'.format(taskid,taskid))
