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
from test_mockLC_lensed_visual import GLSNe
from test_mockLC_lensed_visual import unresolved_and_individual_lcs
import time

start_time = time.time()

# This code was that used to investigate the effect of using priors (bounds) on different fitted parameters.
# In its current state, it fits with x1 and c bounds of [0.2,0.1,0.05], and with no bounds at all, and at the end
# plots all the results to compare bounds with no bounds, and tighter bounds with looser bounds. It can of course
# be simply extended to have bounds on other parameters, and tighter/looser bounds etc..

np.random.seed(2000)

# Params
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

# Band lists
band_list_ztf = ['ztfg','ztfr', 'ztfi']
band_list_full = ['ztfg', 'ztfr', 'ztfi', 'cspjd', 'cspyd', 'csphd']

# Time delay list
mag_list = [mu_1,mu_2]
delta_t_list =[10, 15, 20]

# Creating functions for sn_parameters
def x1_func():
    return x1_val

def c_func():
    return c_val

sn_parameters = {'x1':x1_func,'c':c_func}   # i assume you don't need to specify z because specified in other argument?    
    
# Misc. lists, for storing information as fitting loop iterates
end_delay_residual_list = []   # for storing time delay residuals
bounds_list = [True, False]   #whether or not to use priors/bounds
prior_bounds_used = [] # for keeping track of whether prior bounds were used in a given loop, 1 if yes, 0 if no
dt_list = []
prior_bounds_list = [] # for saving prior bounds themselves to make sure they're sensible
td_error_list = [] # for time delay errors
time_info_list = []
fitted_param_list = []

#%%

# Creating functions for sn_parameters
def x1_func():
    return x1_val

def c_func():
    return c_val

sn_parameters = {'x1':x1_func,'c':c_func}   # i assume you don't need to specify z because specified in other argument?    

x1_bound_list = [0.2,0.1,0.05]
x1_bound_append_list = []

# Begin fitting lightcurves with and without priors
for bounds in bounds_list:
    for i, dt in enumerate(delta_t_list):
        for j, x1_bound in enumerate(x1_bound_list):
            
            # To avoid exceptions
            max_retries = 15
            retry_count = 0
                   
            while retry_count < max_retries:
                try:              
        
                    dt_1 = -dt/2 + 20
                    dt_2 = dt/2 + 20
                    
                    # |||| Importing (resolved) image data ||||
                
                    # Importing data from 'lightcurve generation.py' to ensure all lightcurves used in analysis have same normalisation
                    data_file_name1 = f'td{dt}days-params1-img1-full'
                    output_file_path1 = f'Lightcurve data/{data_file_name1}'
                    
                    data_file_name2 = f'td{dt}days-params1-img2-full'
                    output_file_path2 = f'Lightcurve data/{data_file_name2}'

                    # Open the file in binary mode and load the data using pickle
                    with open(output_file_path1, 'rb') as file:
                        image_1_data = pickle.load(file)
                        
                    with open(output_file_path2, 'rb') as file:
                        image_2_data = pickle.load(file)
                  
                    # Creating MISN instance using data
                    resolved_MISN=sntd.table_factory([image_1_data,image_2_data],telescopename='telescope',object_name='example_SN')

                    c_bound = x1_bound   # setting relation between c and x1 bounds, in this case I made them the same
                    
                    if bounds == True: # fitting with bounds on x1 and c!
                        # Note: I tried to softcode the parameters to fit and constants, but after a long time trying I could not get SNTD to work with it and I could not figure out why
                        # so instead, the params are unfortunately hard-coded in, but it is not such a cumbersome task to change them, just change them manually in this if-else statement.
                        resolved_fitCurves=sntd.fit_data(resolved_MISN,snType='Ia', models='salt2-extended',bands=band_list_ztf,
                                        params=['x0', 't0', 'dt_1', 'dt_2', 'x1', 'c'],constants={'z':z_val, 'mu_1': mu_1, 'mu_2': mu_2},
                                        bounds={'t0':(-25,25), 'dt_1':(-30,10), 'dt_2':(-10,30), 'x1': (-x1_bound,x1_bound), 'c': (-c_bound,c_bound)}, method='parallel',npoints=100)
                        resolved_fitCurves.plot_object(bands='ztfg', showFit=True)
                        plt.show()
                        prior_bounds_used.append(1) # to keep track of when we used bounds and when we didn't
                        
                    else: # fitting without bounds on x1 and c!
                        resolved_fitCurves=sntd.fit_data(resolved_MISN,snType='Ia', models='salt2-extended',bands=band_list_ztf,
                                               params=['x0', 't0', 'dt_1', 'dt_2', 'x1', 'c'],constants={'z':z_val, 'mu_1': mu_1, 'mu_2': mu_2},   # changing order of params entered here doesnt affect fitcurves.model.parameters order
                                               bounds={'t0':(-25,25), 'dt_1':(-30,10), 'dt_2':(-10,30), 'x1': (-3,3), 'c': (-1,1)}, method='parallel',npoints=100)
                        resolved_fitCurves.plot_object(bands='ztfg', showFit=True)
                        plt.show()
                        prior_bounds_used.append(0)
                        
                    # Computing and storing in lists some interesting quantities
                    res_time_delay_SNTD = resolved_fitCurves.parallel.time_delays['image_2'] - resolved_fitCurves.parallel.time_delays['image_1']
                    delay_residual = res_time_delay_SNTD - dt
    
                    end_delay_residual_list.append(delay_residual)
                    dt_list.append(dt)
                    x1_bound_append_list.append(x1_bound)
                    
                    # Finding and storing errors
                    del_errors = resolved_fitCurves.parallel.time_delay_errors
                    down_err = abs(del_errors['image_2'][0])
                    up_err = del_errors['image_2'][1]
    
                    delay_error = [down_err, up_err]
                    td_error_list.append(delay_error)
                    
                    break
        
                except ZeroDivisionError:
                    retry_count += 1
                    continue
                
                except Exception as e:
                    print(f"An error occurred: {e}")
                    break
                    
            if retry_count == max_retries:
                print("Failed after maximum retries.")
            else:
                print(f"Code executed successfully. Took {retry_count + 1} tries.")

        
#%%

# Plotting and analysis

delay_up_errors = []
delay_down_errors = []

for i, err in enumerate(td_error_list): #for each data point
    delay_down_errors.append(err[0])
    delay_up_errors.append(err[1])
    
# Plotting the error on the time delay versus x1, c bounds, with IR data included
plt.title('Time delay error versus x1, c bound (with csp)')
plt.xlabel('x1 bound')
plt.ylabel('Lower + Upper time delay error bar')
for i, dt in enumerate(dt_list):
    if prior_bounds_used[i] == 1:
        if dt == 10:
            plt.scatter(x1_bound_append_list[i], np.array(delay_down_errors[i]) + np.array(delay_up_errors[i]), color = 'blue', label = '10 days')
        if dt == 15:
            plt.scatter(x1_bound_append_list[i], np.array(delay_down_errors[i]) + np.array(delay_up_errors[i]), color = 'red', label = '15 days')
        if dt == 20:
            plt.scatter(x1_bound_append_list[i], np.array(delay_down_errors[i]) + np.array(delay_up_errors[i]), color = 'green', label = '20 days')
    elif prior_bounds_used[i] == 0 and x1_bound_append_list[i] == 0.05:
        if dt == 10:
            plt.scatter(3, np.array(delay_down_errors[i]) + np.array(delay_up_errors[i]), color = 'cyan', label = 'No prior, 10 days')
        if dt == 15:
            plt.scatter(3, np.array(delay_down_errors[i]) + np.array(delay_up_errors[i]), color = 'orange', label = 'No prior, 15 days')
        if dt == 20:
            plt.scatter(3, np.array(delay_down_errors[i]) + np.array(delay_up_errors[i]), color = 'brown', label = 'No prior, 20 days')

handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles)) 
plt.legend(by_label.values(), by_label.keys(), loc='center')
plt.show()
    
# Plotting the actual time delay residual for cases with and without bounds
plt.title('Time delay residuals using x1, c priors (with csp)')
plt.xlabel('Time delay (days)')
plt.ylabel('Time delay residuals (days)')
for i, dt in enumerate(dt_list):
    if prior_bounds_used[i] == 1:
        if x1_bound_append_list[i] == 0.05:
            plt.scatter(dt, end_delay_residual_list[i], label = r'Priors, x1, c [$\pm$ 0.05]', color = 'blue')
            plt.errorbar(dt, end_delay_residual_list[i], yerr = [[delay_down_errors[i]], [delay_up_errors[i]]], linestyle='', color = 'blue')
        elif x1_bound_append_list[i] == 0.1:
            plt.scatter(dt, end_delay_residual_list[i], label = 'Priors, x1, c [$\pm$ 0.1]', color = 'red')
            plt.errorbar(dt, end_delay_residual_list[i], yerr = [[delay_down_errors[i]], [delay_up_errors[i]]], linestyle='', color = 'red')
        else:
            plt.scatter(dt, end_delay_residual_list[i], label = 'Priors, x1, c [$\pm$ 0.2]', color = 'green')
            plt.errorbar(dt, end_delay_residual_list[i], yerr = [[delay_down_errors[i]], [delay_up_errors[i]]], linestyle='', color = 'green')
    if prior_bounds_used[i] == 0:
        if x1_bound_append_list[i] == 0.1:
            plt.scatter(dt, end_delay_residual_list[i], label = r'No priors', color = 'brown')
            plt.errorbar(dt, end_delay_residual_list[i], yerr = [[delay_down_errors[i]], [delay_up_errors[i]]], linestyle='', color = 'brown')

handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles)) 
plt.legend(by_label.values(), by_label.keys(), loc='lower left')
plt.show()




