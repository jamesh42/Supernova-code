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
from astropy.io import ascii
from astropy.table import Table
import sntd
import pickle
import time

# This is the main code I was using. It saves simulations (the original input parameters, the fitted parameters, and 
# information about the simulation (i.e. which parameters are fitted)) as pickle files, for other scripts to unpickle
# and use without worrying about result reproducability or having to spend lots of time rerunning the simulation.

# For timing simulation
start_time = time.time()

# Parameters
z_val = 0.3544   # source supernova redshift
t0_val = 0.
x1_val = 0.
c_val = 0.
x0_val = 6e-05
# Lensing
nimages = 2  # Multiplicity
mu_1 = 2
mu_2 = 7

# Band lists - for choosing which filters to simulate
band_list_ztf = ['ztfg','ztfr', 'ztfi']
band_list_csp = ['cspjd', 'cspyd', 'csphd']
band_list_full = ['ztfg', 'ztfr', 'ztfi', 'cspjd', 'cspyd', 'csphd']
# For choosing which magnifications to use
magnifications_list = [[2,7]]

# Creating zp list in format sntd.createMISN wants
zp_list = [25.] * len(band_list_ztf)

# Time delay list
delta_t_value_list = [10] # time delays to simulates

#%%

# Creating functions for sn_parameters
def x1_func():
    # presumably this seemingly ridiculous functionality is because in reality you'd only know x1 within some distribution, which you'd describe with a function
    return x1_val

def c_func():
    return c_val

sn_parameters = {'x1':x1_func,'c':c_func}   # i assume you don't need to specify z because specified in other argument?

# List of models to use and band lists, etc.
source_model_list = ['salt2-extended', 'hsiao']
source_model_list = ['salt2-extended'] # can simulate salt2 or salt2 and hsiao
band_list_list = [band_list_ztf, band_list_full] #
band_list_labels = ['ztf', 'both']

# Extinction lists
a_v_val_list = [0.5]   # for simulating multiple extinctions

# Lists for appending information to during simulation
a_v_list = []
delay_residual_list = []
delta_t_list = []
source_model_iterated_list = []   # just list to remember each data point's model used
delay_error_list = []
mu_1_list = []
mu_2_list = []
mag_iteration_list = []
fitted_mag_list = []
mag1_err_list = []
mag2_err_list = []

# Starting simulation, which loops over all the values to fit and use in simulation
for a_v_val in a_v_val_list:
    for mag_list in magnifications_list:
        
        residual_list_list = []   # list for storing lists of residuals
        delay_error_list_list = []
        
        for dt_val in delta_t_value_list:
            
            # Creating data  
            dt_1 = -dt_val/2
            dt_2 = dt_val/2
            
            myMISN = sntd.createMultiplyImagedSN(sourcename='salt2-extended', snType='Ia', redshift=z_val, bands=band_list_ztf,
                         zp=zp_list, cadence=5., epochs=35.,time_delays=[dt_1, dt_2], magnifications=mag_list,
             objectName='My Type Ia SN',telescopename='ztf',av_host=a_v_val, numImages = nimages, sn_params = sn_parameters) 
            myMISN.plot_object()
            
            # Saving data to astropy tables - this is our created data
            image_1_dat = myMISN.images['image_1']['table']
            image_2_dat = myMISN.images['image_2']['table']
           
            # Creating MISN instance using data
            new_MISN=sntd.table_factory([image_1_dat,image_2_dat],telescopename='telescope',object_name='example_SN')
            
            dust = sncosmo.CCM89Dust()
         
            # FITTING DATA
            for source_model in source_model_list:
                # I realise here I am just creating the model + data, extracting the data and then recreating the model again, but it's better to be safe and unpack then repack than risk contaminating independent test
                print(source_model)   # Inexplicably if I remove this, the code can break. So I just left it in here
                
                max_retries = 15
                retry_count = 0
                
                while retry_count < max_retries:   # since some runs will not work, but repeating the run can, I try a maximum of 15 times before admitting defeat
                    
                    try:
                        
                        # FOR ANY MODEL EXCEPT SALT2, AMPLITUDE MUST BE USED IN PLACE OF X0
                        if source_model == 'salt2-extended':
                        
                            fitting_model=sncosmo.Model(source_model,effects=[dust,dust],effect_names=['lens','host'],effect_frames=['free','rest'])    
                        
                            fitCurves=sntd.fit_data(new_MISN,snType='Ia', models=fitting_model,bands=band_list_ztf,   # here is where the fitting goes on, choose band list etc. accordingly
                                            params=['x0','x1','t0','c', 'hostebv'],constants={'z':z_val, 'lensebv': 0., 'hostr_v':2.0},
                                            bounds={'t0':(-25,25),'x1':(-2,2),'c':(-1,1), 'hostebv':(0.02,0.333)}, method='parallel',npoints=100)
                            fitCurves.plot_object(bands='ztfg', showFit=True)
                            plt.show()
                            
                            # time delay info
                            time_delay_SNTD = fitCurves.parallel.time_delays['image_2'] - fitCurves.parallel.time_delays['image_1']   # image 1 should always be 0!
                            delay_residual = time_delay_SNTD - dt_val
                            
                            # fitting magnifications info
                            fit_mag_1 = fitCurves.parallel.magnifications['image_1']
                            fit_mag_2 = fitCurves.parallel.magnifications['image_2']
                            
                            fit_mags = [fit_mag_1, fit_mag_2]   # len 2 array with fitted magnifications for image 1 and 2 respectively
                            fitted_mag_list.append(fit_mags)
                            
                            # storing errors
                            mag_1_errs = fitCurves.parallel.magnification_errors['image_1']   # think this is always 0?
                            mag_2_errs = fitCurves.parallel.magnification_errors['image_2']    
                            mag_1_err = [abs(ele) for ele in mag_1_errs]
                            mag_2_err = [abs(ele) for ele in mag_2_errs]
            
                            del_errors = fitCurves.parallel.time_delay_errors
                            down_err = abs(del_errors['image_2'][0])
                            up_err = del_errors['image_2'][1]
            
                            delay_error = [down_err, up_err]
            
                            delay_error_list.append(delay_error)
                            mag1_err_list.append(mag_1_err)
                            mag2_err_list.append(mag_2_err)
                            
                            # appending info to lists
                            a_v_list.append(a_v_val)
                            delay_residual_list.append(delay_residual)
                            delta_t_list.append(dt_val)
                            source_model_iterated_list.append(source_model)
                            mag_iteration_list.append(mag_list)   #for input magnifications
                            
                            break
                
                            
                        else:  # for other source model, requires different specification since salt2 uses x0 and others use 'amplitude'. Besides that this else: section is identical.
                        
                        
                            fitting_model=sncosmo.Model(source_model,effects=[dust,dust],effect_names=['lens','host'],effect_frames=['free','rest'])    
                        
                            fitCurves=sntd.fit_data(new_MISN,snType='Ia', models=fitting_model,bands=band_list_ztf,
                                            params=['amplitude','x1','t0','c', 'hostebv'],constants={'z':z_val, 'lensebv': 0., 'hostr_v':2.0},
                                            bounds={'t0':(-25,25),'x1':(-2,2),'c':(-1,1), 'hostebv':(0.02,0.333)}, method='parallel',npoints=100)
                            fitCurves.plot_object(showFit=True)
                            plt.show()
                            
                            time_delay_SNTD = fitCurves.parallel.time_delays['image_2'] - fitCurves.parallel.time_delays['image_1']
                            delay_residual = time_delay_SNTD - dt_val
                            
                            fit_mag_1 = fitCurves.parallel.magnifications['image_1']
                            fit_mag_2 = fitCurves.parallel.magnifications['image_2']
                            
                            fit_mags = [fit_mag_1, fit_mag_2]
                            fitted_mag_list.append(fit_mags)
                            
                            mag_1_errs = fitCurves.parallel.magnification_errors['image_1']
                            mag_2_errs = fitCurves.parallel.magnification_errors['image_2']    
                            mag_1_err = [abs(ele) for ele in mag_1_errs]
                            mag_2_err = [abs(ele) for ele in mag_2_errs]
            
                            del_errors = fitCurves.parallel.time_delay_errors
                            down_err = abs(del_errors['image_2'][0])
                            up_err = del_errors['image_2'][1]
            
                            delay_error = [down_err, up_err]
            
                            delay_error_list.append(delay_error)
                            mag1_err_list.append(mag_1_err)
                            mag2_err_list.append(mag_2_err)
                            
                            a_v_list.append(a_v_val)
                            delay_residual_list.append(delay_residual)
                            delta_t_list.append(dt_val)
                            source_model_iterated_list.append(source_model)
                            mag_iteration_list.append(mag_list)   #for input magnifications
                            
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
                    print(f"Code executed successfully. Took {retry_count + 1} tries.")
                        
        
end_time = time.time()
time_taken = end_time - start_time

print(f'Time taken for simulation: {time_taken}')        

#%%

# Storing information for pickling
data = Table()
data['A_v'] = a_v_list
data['Time delay residual'] = delay_residual_list
data['Time delay'] = delta_t_list
data['Source model'] = source_model_iterated_list
data['Time delay errors'] = delay_error_list
data['Input magnifications'] = mag_iteration_list
data['Fitted magnifications'] = fitted_mag_list
data['Image 1 fitted mag errors'] = mag1_err_list   # this should always just be zero
data['Image 2 fitted mag errors'] = mag2_err_list

# File path to save the pickle file
fit_output_file_path = 'Pickle/Full sim tables/misc..pkl' #misc so that data is not overwritten accidentally

# Open the file in binary mode and save the data using pickle
with open(fit_output_file_path, 'wb') as file:
    pickle.dump(data, file)
    
#%%

# The more functional part of the code is above. This last section was for an earlier plot I was making of the magnification residuals versus extinction.
# It is probably worth deleting, but I left it in case this test is to be run again. But as mentioned it is more likely that this code is only of further use
# insofar as it can be used to run simulations and store information to pickle files.

mag_1_inputs = []
mag_1_fits = []
mag_2_inputs = []
mag_2_fits = []

for i, mu_input_list in enumerate(mag_iteration_list):
    mag_1_inputs.append(mu_input_list[0])
    mag_2_inputs.append(mu_input_list[1])
    mag_1_fits.append(fitted_mag_list[i][0])
    mag_2_fits.append(fitted_mag_list[i][1])

mag1_residuals = np.array(mag_1_fits) - np.array(mag_1_inputs)
mag2_residuals = np.array(mag_2_fits) - np.array(mag_2_inputs)

delay_up_errors = []
delay_down_errors = []

print(data)
print(data['Image 1 fitted mag errors'])
print(data['Image 2 fitted mag errors'])

for i, err in enumerate(delay_error_list): #for each data point
    delay_down_errors.append(err[0])
    delay_up_errors.append(err[1])
    
mu1_up_errors = []
mu1_down_errors = []
    
for i, err in enumerate(mag1_err_list): #for each data point
    mu1_down_errors.append(err[0])
    mu1_up_errors.append(err[1])
    
mu2_up_errors = []
mu2_down_errors = []
    
for i, err in enumerate(mag2_err_list): #for each data point
    mu2_down_errors.append(err[0])
    mu2_up_errors.append(err[1])

plt.title(r'Time delay residuals for different magnifications ($\Delta t = 10$)')   #fitted only with salt2 (hsiao data not included in this)
plt.xlabel(r'$A_V$')
plt.ylabel(r'Output $\mu$ residuals')

for i, A_v in enumerate(np.array(data['A_v'])):   # i want a plot 
    
    y_err_delays = [[delay_down_errors[i]], [delay_up_errors[i]]] # this needs to be in this format, ie errorbar only accepts form [[down_error], [up_error]] for a given point (or of course [[down_errorS], [up_errorS]])
    y_err_mu1 = [[mu1_down_errors[i]], [mu1_up_errors[i]]]
    y_err_mu2 = [[mu2_down_errors[i]], [mu2_up_errors[i]]]
    
    if np.array(data['Source model'][i]) == 'salt2-extended': 
        if np.array(data["Input magnifications"][i])[0] == 2 and np.array(data["Input magnifications"][i])[1] == 7:
            plt.scatter(A_v, mag1_residuals[i], label=f'$\mu_1$, input mags: {np.array(data["Input magnifications"][i])}', c = 'red', marker = '*')
            plt.scatter(A_v, mag2_residuals[i], label=f'$\mu_2$; input mags: {np.array(data["Input magnifications"][i])}', c='green', marker = '*')
            plt.errorbar(A_v, mag1_residuals[i], yerr = y_err_mu1, c = 'red')
            plt.errorbar(A_v, mag2_residuals[i], yerr = y_err_mu2, c='green')
        else: 
            plt.scatter(A_v, mag1_residuals[i], label=f'$\mu_1$; input mags: {np.array(data["Input magnifications"][i])}', c = 'red', marker = 's')
            plt.scatter(A_v, mag2_residuals[i], label=f'$\mu_2$; input mags: {np.array(data["Input magnifications"][i])}', c='green', marker = 's')
            plt.errorbar(A_v, mag1_residuals[i], yerr = y_err_mu1, c = 'red')
            plt.errorbar(A_v, mag2_residuals[i], yerr = y_err_mu2, c = 'green')
    else:
        pass
    
handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
legend_offset = (0.5, 0.9)  # (x, y) coordinates relative to the axes
plt.legend(by_label.values(), by_label.keys(), loc='upper center', bbox_to_anchor=legend_offset)   #just copied off stack exchange to stop legend repeating itself
plt.show()



