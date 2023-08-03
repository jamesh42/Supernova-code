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

start_time = time.time()

z_val = 0.3544   # source supernova redshift
t0_val = 0.
x1_val = 0.
c_val = 0.
x0_val = 6e-05
# Lensing
nimages = 2  # Multiplicity
mu_1 = 2
mu_2 = 7

# Band lists
band_list_ztf = ['ztfg','ztfr', 'ztfi']
band_list_csp = ['cspjd', 'cspyd', 'csphd']
band_list_full = ['ztfg', 'ztfr', 'ztfi', 'cspjd', 'cspyd', 'csphd']

# Creating zp list in format sntd.createMISN wants
zp_list = [25.] * len(band_list_full)

# Time delay list
delta_t_list = np.linspace(0,40,5)
delta_t_list = [10,20]
delta_t_list= [10, 15, 20]
delta_t_list = [10,20]
#delta_t_list = [10]
delta_t_value_list = [10]

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
source_model_list = ['salt2-extended']
#source_model_list = ['salt2-extended']
band_list_list = [band_list_ztf, band_list_full]
band_list_labels = ['ztf', 'both']

# Extinction lists
a_v_val_list = np.linspace(0.1,0.5, 2)
a_v_val_list = np.linspace(0.1,0.5,5)

#a_v_val_list = [0.25]
#%

a_v_list = []
delay_residual_list = []
delta_t_list = []
source_model_iterated_list = []   # just list to remember each data point's model used
delay_error_list = []
mu_1_list = []
mu_2_list = []
magnifications_list = [[2,7], [5,5]]
mag_iteration_list = []
fitted_mag_list = []
mag1_err_list = []
mag2_err_list = []

print('xxx')
#a_v_val = 0.25

for a_v_val in a_v_val_list:
    for mag_list in magnifications_list:
        
        residual_list_list = []   # list for storing lists of residuals
        delay_error_list_list = []
        
        for dt_val in delta_t_value_list:
            
            # Creating data  
            
            dt_1 = -dt_val/2
            dt_2 = dt_val/2
            
            myMISN = sntd.createMultiplyImagedSN(sourcename='salt2-extended', snType='Ia', redshift=z_val, bands=band_list_full,
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
                # I am pretty sure here I am just creating the model + data, extracting the data and then recreating the model again, but it's better to be safe and unpack then repack than risk contaminating independent test
                print(source_model)   # IF I REMOVE THIS, THE CODE BREAKS -- UTTERLY BIZARRE
                # Fitting again with SNTD
                
                max_retries = 15
                retry_count = 0
                
                while retry_count < max_retries:
                    
                    try:
                        
                        # FOR ANY MODEL EXCEPT SALT2, AMPLITUDE MUST BE USED IN PLACE OF X0
                        if source_model == 'salt2-extended':
                        
                            fitting_model=sncosmo.Model(source_model,effects=[dust,dust],effect_names=['lens','host'],effect_frames=['free','rest'])    
                        
                            fitCurves=sntd.fit_data(new_MISN,snType='Ia', models=fitting_model,bands=band_list_full,
                                            params=['x0','x1','t0','c', 'hostebv'],constants={'z':z_val, 'lensebv': 0., 'hostr_v':2.0},
                                            bounds={'t0':(-25,25),'x1':(-2,2),'c':(-1,1), 'hostebv':(0.02,0.333)}, method='parallel',npoints=100)
                            fitCurves.plot_object(showFit=True)
                            plt.show()
                            
                            time_delay_SNTD = fitCurves.parallel.time_delays['image_2'] - fitCurves.parallel.time_delays['image_1']   # image 1 should always be 0!
                            delay_residual = time_delay_SNTD - dt_val
                            
                            fit_mag_1 = fitCurves.parallel.magnifications['image_1']
                            fit_mag_2 = fitCurves.parallel.magnifications['image_2']
                            
                            fit_mags = [fit_mag_1, fit_mag_2]   # len 2 array with fitted magnifications for image 1 and 2 respectively
                            fitted_mag_list.append(fit_mags)
                            
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
                            
                            a_v_list.append(a_v_val)
                            delay_residual_list.append(delay_residual)   #for trying new method, 2 denotes new method
                            delta_t_list.append(dt_val)
                            source_model_iterated_list.append(source_model)
                            mag_iteration_list.append(mag_list)   #for input magnifications
                            
                            break
                
                            
                        else:
                        
                        
                            fitting_model=sncosmo.Model(source_model,effects=[dust,dust],effect_names=['lens','host'],effect_frames=['free','rest'])    
                        
                            fitCurves=sntd.fit_data(new_MISN,snType='Ia', models=fitting_model,bands=band_list_full,
                                            params=['amplitude','x1','t0','c', 'hostebv'],constants={'z':z_val, 'lensebv': 0., 'hostr_v':2.0},
                                            bounds={'t0':(-25,25),'x1':(-2,2),'c':(-1,1), 'hostebv':(0.02,0.333)}, method='parallel',npoints=100)
                            fitCurves.plot_object(showFit=True)
                            plt.show()
                            
                            time_delay_SNTD = fitCurves.parallel.time_delays['image_2'] - fitCurves.parallel.time_delays['image_1']
                            delay_residual = time_delay_SNTD - dt_val
                            
                            fit_mag_1 = fitCurves.parallel.magnifications['image_1']
                            fit_mag_2 = fitCurves.parallel.magnifications['image_2']
                            
                            fit_mags = [fit_mag_1, fit_mag_2]   # len 2 array with fitted magnifications for image 1 and 2 respectively
                            fitted_mag_list.append(fit_mags)
                            
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
                            
                            a_v_list.append(a_v_val)
                            delay_residual_list.append(delay_residual)   #for trying new method, 2 denotes new method
                            delta_t_list.append(dt_val)
                            source_model_iterated_list.append(source_model)
                            mag_iteration_list.append(mag_list)   #for input magnifications
                            
                            break
                            
                    except ZeroDivisionError:
                        retry_count += 1   # USE ENUMERATE IN THIS 5:21 of python vid
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

print(f'Time taken: {time_taken}')        

#%%

data = Table()
data['A_v'] = a_v_list # list of times for which we have data
data['Time delay residual'] = delay_residual_list   # corresponding list of simulated fluxes for each time
data['Time delay'] = delta_t_list   # applying systematic error
data['Source model'] = source_model_iterated_list
data['Time delay errors'] = delay_error_list
data['Input magnifications'] = mag_iteration_list
data['Fitted magnifications'] = fitted_mag_list
data['Image 1 fitted mag errors'] = mag1_err_list   # this should always just be zero
data['Image 2 fitted mag errors'] = mag2_err_list

mag_1_inputs = []
mag_1_fits = []
mag_2_inputs = []
mag_2_fits = []

for i, mu_input_list in enumerate(mag_iteration_list):
    mag_1_inputs.append(mu_input_list[0])
    mag_2_inputs.append(mu_input_list[1])
    mag_1_fits.append(fitted_mag_list[i][0])
    mag_2_fits.append(fitted_mag_list[i][1])


print(mag_1_fits)
print(mag_2_fits)
mag1_residuals = np.array(mag_1_fits) - np.array(mag_1_inputs)
mag2_residuals = np.array(mag_2_fits) - np.array(mag_2_inputs)

# I believe this ^ is all the key information from this run, besides input parameters and parameters varied, 
# but I will specify those in the overleaf. This is all the information needed to plot relevant graphs, 
# so this is what will be pickled and then can be re-plotted in any way desired in the future.

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
    
# File path to save the pickle file
fit_output_file_path = 'Pickle/Full sim tables/misc.pkl'  #specify which one!

# Open the file in binary mode and save the data using pickle
with open(fit_output_file_path, 'wb') as file:
    pickle.dump(data, file)

# av_min = a_v_list[0]
# av_max = a_v_list[-1]

plt.title(r'Time delay residuals for different magnifications ($\Delta t = 10$)')   #fitted only with salt2 (hsiao data not included in this)
plt.xlabel(r'$A_V$')
plt.ylabel(r'Output $\mu$ residuals')
# plt.ylim(-1,1)
#plt.ylim(-5,0)
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
# Place the legend slightly above the center
legend_offset = (0.5, 0.9)  # (x, y) coordinates relative to the axes
plt.legend(by_label.values(), by_label.keys(), loc='upper center', bbox_to_anchor=legend_offset)   #just copied off stack exchange to stop legend repeating itself
plt.show()

# So red denotes mu1


# plt.title(r'Time delay residuals for various fitting models')
# plt.xlabel(r'Original time delay (days)')
# plt.ylabel(r'Output time delay residuals (days)')
# plt.ylim(-2,2)
# for i, A_v in enumerate(np.array(data['A_v'])):
    
#     y_err = [[down_errors[i]], [up_errors[i]]] # this needs to be in this format, ie errorbar only accepts form [[down_error], [up_error]] for a given point (or of course [[down_errorS], [up_errorS]])
    
#     if np.array(data['Source model'][i]) == 'salt2-extended': 
#         plt.scatter(np.array(data['Time delay'][i]), np.array(data['Time delay residual'][i]), marker = 'x', c = A_v, vmin = av_min, vmax = av_max, label = 'salt2-extended')
#         plt.errorbar(np.array(data['Time delay'][i]), np.array(data['Time delay residual'][i]), yerr = y_err, color=plt.cm.viridis((A_v - av_min) / (av_max - av_min)))
#     elif np.array(data['Source model'][i]) == 'hsiao': 
#         plt.scatter(np.array(data['Time delay'][i]), np.array(data['Time delay residual'][i]), marker = 'o', c = A_v, vmin = av_min, vmax = av_max, label = 'hsiao')
#         plt.errorbar(np.array(data['Time delay'][i]), np.array(data['Time delay residual'][i]), yerr = y_err, color=plt.cm.viridis((A_v - av_min) / (av_max - av_min)))

# plt.colorbar(label=r'$A_V$')
# handles, labels = plt.gca().get_legend_handles_labels()
# by_label = dict(zip(labels, handles))
# plt.legend(by_label.values(), by_label.keys(), loc='upper right')   #just copied off stack exchange to stop legend repeating itself
# plt.show()

# USE ENUMERATE FOR ALL FUTURE LOOPS
# ALSO USE f STRINGS FOR ALL LABELS ETC.

# plot lightcurves IN OVERLEAF AS KAISEY SAYS, USEFUL TO SEE, COULD IR BE MISSING THE SECOND PEAK AT LARGE TIME DELAYS SINCE 
# FINITE TIME SERIES SO DATA IS TRUNCATED CAUSING LOSS OF INFORMATION AND THEREFORE WORSE FIT




#%%

"""
Questions: 
    
    -what is t1_error = fitCurves.images['image_1'].fits.res.errors['t0']? I have no idea what this, and the corresponding image_2 quantity, is
    - what is relationship between t1_error above and the parallel.time_delay_errors ?
    
    -  #t1 = fitCurves.images['image_1'].fits.res.param_dict['t0']   # why  on earth would image_1 not work - it turns out not to exist! Maybe by definition set to 0?
    - here images.['image_1']..... param dict doesnt exist but it does for image 2, why?
    
""" 


