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
#delta_t_list = [10]

#%%ß

# Creating functions for sn_parameters
def x1_func():
    # presumably this seemingly ridiculous functionality is because in reality you'd only know x1 within some distribution, which you'd describe with a function
    return x1_val

def c_func():
    return c_val

sn_parameters = {'x1':x1_func,'c':c_func}   # i assume you don't need to specify z because specified in other argument?

# List of models to use and band lists, etc.
source_model_list = ['salt2-extended', 'hsiao']
#source_model_list = ['salt2-extended']
band_list_list = [band_list_ztf, band_list_full]
band_list_labels = ['ztf', 'both']

# Extinction lists
a_v_val_list = np.linspace(0.1,0.5, 2)
a_v_val_list = np.linspace(0.1,0.5,5)
#a_v_val_list = [0.1]

# residual_list_list_list = []
# delay_error_list_list_list = []
#%%

a_v_list2 = []
residual_list2 = []
delta_t_list2 = []

N = len(a_v_val_list) * len(source_model_list) * len(delta_t_list)
info_array = np.zeros((N,0))   # but this will surely break if any of the loops fail

source_model_iterated_list = []   # just list to remember each data point's model used
delay_error_list = []

for a_v_val in a_v_val_list:
    
    residual_list_list = []   # list for storing lists of residuals
    delay_error_list_list = []
    
    for dt_val in delta_t_list:
        
        # Creating data  
        
        dt_1 = -dt_val/2
        dt_2 = dt_val/2
        
        myMISN = sntd.createMultiplyImagedSN(sourcename='salt2-extended', snType='Ia', redshift=z_val, bands=band_list_full,
                     zp=zp_list, cadence=5., epochs=35.,time_delays=[dt_1, dt_2], magnifications=[mu_1,mu_2],
         objectName='My Type Ia SN',telescopename='ztf',av_host=a_v_val, numImages = nimages, sn_params = sn_parameters) 
        myMISN.plot_object()
        
        
        # Saving data to astropy tables - this is our created data
        image_1_dat = myMISN.images['image_1']['table']
        image_2_dat = myMISN.images['image_2']['table']
       
        # Creating MISN instance using data
        new_MISN=sntd.table_factory([image_1_dat,image_2_dat],telescopename='telescope',object_name='example_SN')
        
        # residual_list = []   # stores residuals for all sources but at same time delay

     
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
                    
                        fitCurves=sntd.fit_data(new_MISN,snType='Ia', models=source_model,bands=band_list_full,
                                        params=['x0','x1','t0','c'],constants={'z':z_val},
                                        bounds={'t0':(-25,25),'x1':(-2,2),'c':(-1,1)}, method='parallel',npoints=100)
                        fitCurves.plot_object(showFit=True)
                        plt.show()
                        
                        time_delay_SNTD = fitCurves.parallel.time_delays['image_2'] - fitCurves.parallel.time_delays['image_1']
                        residual = time_delay_SNTD - dt_val
                        #residual_list.append(residual)
        
                        del_errors = fitCurves.parallel.time_delay_errors
                        print(del_errors)
                        down_err = abs(del_errors['image_2'][0])
                        up_err = del_errors['image_2'][1]
        
                        delay_error = [down_err, up_err]
        
                        delay_error_list.append(delay_error)
                        
                        a_v_list2.append(a_v_val)
                        residual_list2.append(residual)   #for trying new method, 2 denotes new method
                        delta_t_list2.append(dt_val)
                        source_model_iterated_list.append(source_model)
                        
                        # # File path to save the pickle file
                        # fit_output_file_path = 'Pickle/Specific fits/lc_and_fit.pkl'
                        
                        # fit_data_to_save = {'og_MISN': myMISN, 'specific_fit': fitCurves}
                        
                        # # Open the file in binary mode and save the data using pickle
                        # with open(fit_output_file_path, 'wb') as file:
                        #     pickle.dump(fit_data_to_save, file)
                        
                        # ^ IS THIS REALLY NECESSARY? SURELY IT'S MUCH MORE EFFICIENT JUST TO STORE ALL RELEVANT VARIABLES IN 'data' TABLE AT END, rather than keeping
                        # individual fits? Seems like it is needlessly complicated, it's probably always quicker to just run the sim again than rummage through
                        #pickle and find the exact lc fit you want?
                        
                        break
            
                        
                    else:
                    
                        fitCurves=sntd.fit_data(new_MISN,snType='Ia', models=source_model,bands=band_list_full,
                                        params=['amplitude','x1','t0','c'],constants={'z':z_val},
                                        bounds={'t0':(-25,25),'x1':(-2,2),'c':(-1,1)}, method='parallel',npoints=100)
                        fitCurves.plot_object(showFit=True)
                        plt.show()
                        
                        time_delay_SNTD = fitCurves.parallel.time_delays['image_2'] - fitCurves.parallel.time_delays['image_1']
                        residual = time_delay_SNTD - dt_val
                        #residual_list.append(residual)
        
                        del_errors = fitCurves.parallel.time_delay_errors
                        print(del_errors)
                        down_err = abs(del_errors['image_2'][0])
                        up_err = del_errors['image_2'][1]
                        print(down_err, up_err)
                        
                        delay_error = [down_err, up_err]
        
                        delay_error_list.append(delay_error)
                        
                        a_v_list2.append(a_v_val)
                        residual_list2.append(residual)
                        delta_t_list2.append(dt_val)
                        source_model_iterated_list.append(source_model)
                        
                        # # File path to save the pickle file
                        # fit_output_file_path = 'Pickle/Specific fits/lc_and_fit.pkl'
                        
                        # fit_data_to_save = {'og_MISN': myMISN, 'specific_fit': fitCurves}
                        
                        # # Open the file in binary mode and save the data using pickle
                        # with open(fit_output_file_path, 'wb') as file:
                        #     pickle.dump(fit_data_to_save, file)
                        
                        
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
                print("Code executed successfully. Took " + str(retry_count+1) + " tries." )
                    
        
    #     residual_list_list.append(residual_list)
    #     delay_error_list_list.append(delay_error_list)
    # delay_error_list_list_list.append(delay_error_list_list)   #silly
    # residual_list_list_list.append(residual_list_list)
            
#%%

data = Table()
data['A_v'] = a_v_list2 # list of times for which we have data
data['Time delay residual'] = residual_list2   # corresponding list of simulated fluxes for each time
data['Time delay'] = delta_t_list2   # applying systematic error
data['Source model'] = source_model_iterated_list
data['Time delay errors'] = delay_error_list

# I believe this ^ is all the key information from this run, besides input parameters and parameters varied, 
# but I will specify those in the overleaf. This is all the information needed to plot relevant graphs, 
# so this is what will be pickled and then can be re-plotted in any way desired in the future.

up_errors = []
down_errors = []

for i, err in enumerate(delay_error_list): #for each data point
    down_errors.append(err[0])
    up_errors.append(err[1])
    
print(up_errors)


# File path to save the pickle file
fit_output_file_path = 'Pickle/Full sim tables/misc.pkl'  #specify which one!

# Open the file in binary mode and save the data using pickle
with open(fit_output_file_path, 'wb') as file:
    pickle.dump(data, file)

print(delay_error_list)
print(data)

#print(delay_error_list[0])

av_min = a_v_list2[0]
av_max = a_v_list2[-1]

plt.title(r'Time delay residuals for various fitting models')
plt.xlabel(r'Original time delay (days)')
plt.ylabel(r'Output time delay residuals (days)')
#plt.ylim(-0,2)

for i, A_v in enumerate(np.array(data['A_v'])):
    
    y_err = [[down_errors[i]], [up_errors[i]]] # this needs to be in this format, ie errorbar only accepts form [[down_error], [up_error]] for a given point (or of course [[down_errorS], [up_errorS]])
    
    if np.array(data['Source model'][i]) == 'salt2-extended': 
        plt.scatter(np.array(data['Time delay'][i]), np.array(data['Time delay residual'][i]), marker = 'x', c = A_v, vmin = av_min, vmax = av_max, label = 'salt2-extended')
        plt.errorbar(np.array(data['Time delay'][i]), np.array(data['Time delay residual'][i]), yerr = y_err, color=plt.cm.viridis((A_v - av_min) / (av_max - av_min)))
    elif np.array(data['Source model'][i]) == 'hsiao': 
        plt.scatter(np.array(data['Time delay'][i]), np.array(data['Time delay residual'][i]), marker = 'o', c = A_v, vmin = av_min, vmax = av_max, label = 'hsiao')
        plt.errorbar(np.array(data['Time delay'][i]), np.array(data['Time delay residual'][i]), yerr = y_err, color=plt.cm.viridis((A_v - av_min) / (av_max - av_min)))

plt.colorbar(label=r'$A_V$')
handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys(), loc='upper right')   #just copied off stack exchange to stop legend repeating itself
plt.show()

# USE ENUMERATE FOR ALL FUTURE LOOPS
# ALSO USE f STRINGS FOR ALL LABELS ETC.

# plot lightcurves IN OVERLEAF AS KAISEY SAYS, USEFUL TO SEE, COULD IR BE MISSING THE SECOND PEAK AT LARGE TIME DELAYS SINCE 
# FINITE TIME SERIES SO DATA IS TRUNCATED CAUSING LOSS OF INFORMATION AND THEREFORE WORSE FIT


#%%

# colour_list = ['blue', 'orange', 'red', 'green']

# # Define the color list based on the combinations of k and j
# color_list = {
#     (0, 0): 'blue',
#     (0, 1): 'orange',
#     (1, 0): 'red',
#     (1, 1): 'green',
#     (2,0): 'purple',
#     (2,1): 'cyan',
    
#     # Add more combinations as needed
# }


# plt.title(r'Time delay residuals for various fitting models')
# plt.xlabel(r'Original time delay (days)')
# plt.ylabel(r'Output time delay residuals (days)')
# #.ylim(3,-3)
# color_count = 0
# for k in range(len(a_v_val_list)):
#     for i in range(len(delta_t_list)):
#         for j in range(len(source_model_list)):
#             # Use the color mapping based on k and j
#             color = color_list[(k, j)]
#             #it_label = source_model_list[j] + ' Av=' + str(a_v_val_list[k])
#             it_label = f'{source_model_list[j]}, Av = {a_v_val_list[k]}'
#             plt.plot(delta_t_list[i], residual_list_list_list[k][i][j], label = it_label, marker='o', color = color)
#             y_err = [[delay_error_list_list_list[k][i][j][0]], [delay_error_list_list_list[k][i][j][1]]]   #getting that data points error bar values
#             plt.errorbar(delta_t_list[i], residual_list_list_list[k][i][j], yerr=y_err, label = it_label, marker='o', color = color)
# handles, labels = plt.gca().get_legend_handles_labels()
# by_label = dict(zip(labels, handles))
# plt.legend(by_label.values(), by_label.keys(), loc='upper left')   #just copied off stack exchange to stop legend repeating itself
# #plt.legend()
# plt.show()

# USE ENUMERATE FOR ALL FUTURE LOOPS
# ALSO USE f STRINGS FOR ALL LABELS ETC.





#%%

"""
Questions: 
    
    -what is t1_error = fitCurves.images['image_1'].fits.res.errors['t0']? I have no idea what this, and the corresponding image_2 quantity, is
    - what is relationship between t1_error above and the parallel.time_delay_errors ?
    
    -  #t1 = fitCurves.images['image_1'].fits.res.param_dict['t0']   # why  on earth would image_1 not work - it turns out not to exist! Maybe by definition set to 0?
    - here images.['image_1']..... param dict doesnt exist but it does for image 2, why?
    
""" 


