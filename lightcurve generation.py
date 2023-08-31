#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 16:28:38 2023

@author: jameshenderson
"""

import sncosmo
import numpy as np
import matplotlib.pyplot as plt
import warnings
import corner
from astropy.io import ascii
from astropy.table import Table
import astropy
import sntd
import pickle
from test_mockLC_lensed_visual import GLSNe
from test_mockLC_lensed_visual import unresolved_and_individual_lcs
import time

# LIGHTCURVE GENERATION FILE

# Params -- these are not being varied
z_val = 0.3544   # source supernova redshift
t0_val = 0.
x1_val = 0.
c_val = 0.
x0_val = 6e-05
cad = 5. #cadence
# Lensing
nimages = 2  # Multiplicity
mu_1 = 2
mu_2 = 7
mag_list = [mu_1, mu_2]

with open('params1.txt', 'w') as file:
    # Write data to the file
    file.write(f'z: {z_val}, t0: {t0_val}, x1: {x1_val}, c: {c_val}, x0: {x0_val}, cad: {cad}, nimages: {nimages}, mu1: {mu_1}, mu2: {mu_2}')



time_delay_list = [5,10,15,20,25,30]

band_list_ztf = ['ztfg','ztfr', 'ztfi']

# Creating functions for sn_parameters
def x1_func():
    return x1_val

def c_func():
    return c_val

sn_parameters = {'x1':x1_func,'c':c_func}   # i assume you don't need to specify z because specified in other argument?    
    
# Creating zp list in format sntd.createMISN wants
zp_list = [25.] * len(band_list_ztf)

for dt in time_delay_list:
    
    dt_1 = -dt/2 + 20
    dt_2 = dt/2 + 20
    
    myMISN = sntd.createMultiplyImagedSN(sourcename='salt2-extended', snType='Ia', redshift=z_val, bands=band_list_ztf,
                  zp=zp_list, cadence=cad, epochs=35.,time_delays=[dt_1, dt_2], magnifications=mag_list,
      objectName='My Type Ia SN',telescopename='ztf',av_host=False, numImages = nimages, sn_params = sn_parameters) 
    myMISN.plot_object(bands='ztfg')

    # Saving data to astropy tables - this is our created data
    image_1_dat = myMISN.images['image_1']['table']
    image_2_dat = myMISN.images['image_2']['table']
    image_dat = astropy.table.vstack([image_1_dat, image_2_dat])
    print(image_dat)
    
    # Creating MISN instance using data
    # resolved_MISN=sntd.table_factory([image_1_dat,image_2_dat],telescopename='telescope',object_name='example_SN')
    #print(resolved_MISN)
    # print(resolved_MISN.images['image_1']['__dict__']['table'])
    # print(resolved_MISN.images['image_2']['__dict__']['table'])
    
    # Saving image 1
    data_file_name1 = f'td{dt}days-params1-img1'
    # ascii.write(image_1_dat, data_file_name1, overwrite=True)    #writing this lensed lc info to file to be read and fitted to by model 2
    # data = sncosmo.read_lc(data_file_name)   # producing lensed lc data
    
    fit_output_file_path = f'Lightcurve data/{data_file_name1}'
    
    with open(fit_output_file_path, 'wb') as file:
        pickle.dump(image_1_dat, file)
    
    # Saving image 2
    data_file_name2 = f'td{dt}days-params1-img2'
    # ascii.write(image_1_dat, data_file_name1, overwrite=True)    #writing this lensed lc info to file to be read and fitted to by model 2
    # data = sncosmo.read_lc(data_file_name)   # producing lensed lc data
    
    fit_output_file_path = f'Lightcurve data/{data_file_name2}'
    
    with open(fit_output_file_path, 'wb') as file:
        pickle.dump(image_2_dat, file)












