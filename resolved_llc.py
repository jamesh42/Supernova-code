#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 10:26:48 2023

@author: jameshenderson
"""

import sncosmo
import numpy as np
import matplotlib.pyplot as plt
import warnings
import corner
warnings.filterwarnings('ignore') # to make output easier to read on jupyter nb's
from astropy.io import ascii
from astropy.table import Table
#from SD_resolved import GLSNe
#from SD_resolved import unresolved_and_individual_lcs
import sntd

# SNTD method

z_val = 0.3544   # source supernova redshift
t0_val = 0.
x1_val = 0.
c_val = 0.
x0_val = 6e-05
# Lensing
nimages = 2  # Multiplicity
# delta_t = 20
# min_delta_t = 10
# max_delta_t = 30
mu_1 = 2
mu_2 = 7

# dt_val =10

# dt_1 = -dt_val/2
# dt_2 = dt_val/2

band_list = ['ztfg','ztfr', 'ztfi']

zp_list = [25.] * len(band_list)

# dt_1 = 20
# dt_2 = 70

# delta_t_list = [50]

delta_t_list = np.linspace(0,40,5)

delta_t_list = [10]
#%%

# CREATING DATA

def x1_func():
    # presumably this seemingly ridiculous functionality is because in reality you'd only know x1 within some distribution, which you'd describe with a function
    return x1_val

def c_func():
    return c_val

sn_parameters = {'x1':x1_func,'c':c_func}   # i assume you don't need to specify z because specified in other argument?

residual_list_list = []   # list for storing lists of residuals
delay_error_list_list = []


source_model_list = ['salt2-extended', 'hsiao']
#source_model_list = ['salt2-extended']


#%%
for dt_val in delta_t_list:
    # Creating data  
    
    dt_1 = -dt_val/2
    dt_2 = dt_val/2
    
    myMISN = sntd.createMultiplyImagedSN(sourcename='salt2-extended', snType='Ia', redshift=z_val, bands=band_list,
                 zp=[25., 25., 25.], cadence=5., epochs=35.,time_delays=[dt_1, dt_2], magnifications=[mu_1,mu_2],
     objectName='My Type Ia SN',telescopename='ztf',av_host=False, numImages = nimages, sn_params = sn_parameters)
    
    myMISN.plot_object()
    
    # Saving data to astropy tables - this is our created data
    image_1_dat = myMISN.images['image_1']['table']
    image_2_dat = myMISN.images['image_2']['table']
    
   
    # FITTING DATA

    
    # Creating MISN instance using data
    new_MISN=sntd.table_factory([image_1_dat,image_2_dat],telescopename='telescope',object_name='example_SN')
    print(new_MISN)
    
    residual_list = []   # stores residuals for all sources but at same time delay
    delay_error_list = []
 
    for source_model in source_model_list:
        # I am pretty sure here I am just creating the model + data, extracting the data and then recreating the model again, but it's better to be safe and unpack then repack than risk contaminating independent test
        print(source_model)   # IF I REMOVE THIS, THE CODE BREAKS -- UTTERLY BIZARRE
        # Fitting again with SNTD
        
    max_retries = 10
    retry_count = 0
    
    while retry_count < max_retries:
        
        try:
            
            # FOR ANY MODEL EXCEPT SALT2, AMPLITUDE MUST BE USED IN PLACE OF X0
            if source_model == 'salt2-extended':
            
                fitCurves=sntd.fit_data(new_MISN,snType='Ia', models=source_model,bands=band_list,
                                params=['x0','x1','t0','c'],constants={'z':z_val},
                                bounds={'t0':(-25,25),'x1':(-2,2),'c':(-1,1)}, method='parallel',npoints=100)
                fitCurves.plot_object(showFit=True)
                plt.show()
                
                time_delay_SNTD = fitCurves.parallel.time_delays['image_2'] - fitCurves.parallel.time_delays['image_1']
                print(time_delay_SNTD)
                residual = time_delay_SNTD - dt_val   #output time delay minus input time delay
                print('residual: ' + str(residual) + ' days')
                residual_list.append(residual)
                t1_error = fitCurves.images['image_1'].fits.res.errors['t0']
                t2_error = fitCurves.images['image_2'].fits.res.errors['t0']
                # So, assuming that t0 in SNTD is really the time delay
                total_t_error = np.sqrt(t1_error**2 + t2_error**2)
                
                delay_error_list.append(total_t_error)
                
                break
    
                
            else:
            
                fitCurves=sntd.fit_data(new_MISN,snType='Ia', models=source_model,bands=band_list,
                                params=['amplitude','x1','t0','c'],constants={'z':z_val},
                                bounds={'t0':(-25,25),'x1':(-2,2),'c':(-1,1)}, method='parallel',npoints=100)
                fitCurves.plot_object(showFit=True)
                plt.show()
                
                time_delay_SNTD = fitCurves.parallel.time_delays['image_2'] - fitCurves.parallel.time_delays['image_1']
                print(time_delay_SNTD)
                residual = time_delay_SNTD - dt_val
                print('residual: ' + str(residual) + ' days')
                residual_list.append(residual)
                #print('errors!')
                #print(type(fitCurves.images))
                print(fitCurves.images['image_1'].fits.res.errors['t0'])
                print(fitCurves.images['image_2'].fits.res.errors['t0'])
                
                
                print(fitCurves.parallel.time_delays)
                print(fitCurves.parallel.time_delay_errors)   #why are these arrays???
                print('image 1')
                print(fitCurves.images['image_1'].fits.res)
                print('image 2')
                print(fitCurves.images['image_2'].fits.res)
    
                #t1 = fitCurves.images['image_1'].fits.res.param_dict['t0']   # why  on earth would image_1 not work - it turns out not to exist! Maybe by definition set to 0?
                t2 = fitCurves.images['image_2'].fits.res.param_dict['t0']
                #print('xxx')
                #print(time_delay_SNTD - (t2-t1))
                print(t2)
                
                t1_error = fitCurves.images['image_1'].fits.res.errors['t0']    #what is this quantity? turns out it's not what I thought.....
                t2_error = fitCurves.images['image_2'].fits.res.errors['t0']
                print('errors:')
                print(fitCurves.parallel.time_delay_errors)
                # So, assuming that t0 in SNTD is really the time delay
                total_t_error = np.sqrt(t1_error**2 + t2_error**2)
                print(total_t_error)
                
                delay_error_list.append(total_t_error)
                # NOW YOU JUST NEED TO PLOT IT! (ALSO CHECK WHETHER t0 IS THE TIME DELAY!)
                break
                
        except ZeroDivisionError:
            retry_count += 1
            continue
        
        except Exception:
            print(Exception)
            break
            
    if retry_count == max_retries:
        print("Failed after maximum retries.")
    else:
        print("Code executed successfully. Took " + str(retry_count+1) + " tries." )
        
            
            
            
            #27/7/2023
            
            
            
            
    
    residual_list_list.append(residual_list)
    delay_error_list_list.append(delay_error_list)
        
#%%

print(residual_list)
print(residual_list_list)
print(delay_error_list_list)

colour_list = ['blue', 'orange']

plt.title(r'Time delay residuals for various fitting models')
plt.xlabel(r'Original time delay (days)')
plt.ylabel(r'Output time delay residuals (days)')
for i in range(len(delta_t_list)):
    for j in range(len(source_model_list)):
        plt.plot(delta_t_list[i], residual_list_list[i][j], label = source_model_list[j], marker='o', c = colour_list[j])
        plt.errorbar(delta_t_list[i], residual_list_list[i][j], yerr=delay_error_list_list[i][j], label = source_model_list[j], marker='o', c = colour_list[j])
handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys())   #just copied off stack exchange to stop legend repeating itself
#plt.legend()
plt.show()



"""
Questions: 
    
    -what is t1_error = fitCurves.images['image_1'].fits.res.errors['t0']? I have no idea what this, and the corresponding image_2 quantity, is

"""


