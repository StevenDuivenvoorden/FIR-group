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

from astropy.io import registry
from astropy.table.info import serialize_method_as


priors,posterior = xidplus.load('data/fir/PACS/xidplus_run_5/lofar_xidplus_fir_5.pkl')
prior100,prior160 = priors


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
xidplus.save([prior100,prior160],posterior,'data/fir/PACS/xidplus_run_{}/lofar_xidplus_fir_{}'.format(taskid,taskid))
#the next couple of lines are an alternative way to save astropy table since the Table.write method is currently broken
with serialize_method_as(PACS_cat, None):
            registry.write(PACS_cat,'data/fir/PACS/xidplus_run_{}/lofar_xidplus_fir_{}.fits'.format(taskid,taskid),format='fits')    
#Table.write(PACS_cat,'data/fir/PACS/xidplus_run_{}/lofar_xidplus_fir_{}_rerun.fits'.format(taskid,taskid),overwrite=True)


