#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  4 12:38:11 2023

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

# This code is designed to work with saved lightcurves. Lightcurve generation.py saves (pickles) generated supernovae lightcurves
# and this code opens those lightcurves and manipulates them in whichever way desired. The reason I stopped generating
# and fitting lightcurves in the same file was because my lightcurves would change normalisation seemingly randomly and 
# I wanted to make sure I was using the same lightcurve each time for consistency.


#%%
#  Varying model and Av whilst fixing Av=0 in fitting and only simulating with SALT2.

# File path to the pickle file
output_file_path = 'Pickle/Full sim tables/varying_model_and_host_A_v.pkl'

# Open the file in binary mode and load the data using pickle
with open(output_file_path, 'rb') as file:
    loaded_data = pickle.load(file)

print(loaded_data)
# Data can now be used as required.

#%%

#  Varying model and Av but also fitting with Av

# File path to the pickle file
output_file_path = 'Pickle/Full sim tables/varying_model_and_host_a_v_fit_a_v_too.pkl'

# Open the file in binary mode and load the data using pickle
with open(output_file_path, 'rb') as file:
    loaded_data = pickle.load(file)

print(loaded_data)
# Data can now be used as required.


#%%

# Seeing effect of different magnifications on the time delays

# File path to the pickle file
output_file_path = 'Pickle/Full sim tables/mags_vs_time_delays.pkl' 

# Open the file in binary mode and load the data using pickle
with open(output_file_path, 'rb') as file:
    loaded_data = pickle.load(file)

print(loaded_data)
# Data can now be used as required.


#%%

# Residual fit magnifications vs Av (Av still held as a free parameter)

# File path to the pickle file
output_file_path = 'Pickle/Full sim tables/res_mags_vs_av.pkl'

# Open the file in binary mode and load the data using pickle
with open(output_file_path, 'rb') as file:
    data = pickle.load(file)

print(data)
 
# Creating lists to store original and fitted magnifications with
mu1_inputs = [] # mu1 original value
mu1_fits = [] # mu1 fitted value
mu2_inputs = []   # etc.
mu2_fits = []

# Adding relevant data to those lists from pickled file
for i, mu_input_list in enumerate(np.array(data['Input magnifications'])):
    mu1_inputs.append(mu_input_list[0])
    mu2_inputs.append(mu_input_list[1])
    mu1_fits.append(data['Fitted magnifications'][i][0])
    mu2_fits.append(data['Fitted magnifications'][i][1])

# Creating and filling lists of the relevant errors from the data
delay_up_errors = []
delay_down_errors = []
for i, err in enumerate(np.array(data['Time delay errors'])): #for each data point
    delay_down_errors.append(err[0])
    delay_up_errors.append(err[1])
    
mu1_up_errors = []
mu1_down_errors = []
for i, err in enumerate(np.array(data['Image 1 fitted mag errors'])): #for each data point
    mu1_down_errors.append(err[0])
    mu1_up_errors.append(err[1])
    
mu2_up_errors = []
mu2_down_errors = []
for i, err in enumerate(np.array(data['Image 2 fitted mag errors'])): #for each data point
    mu2_down_errors.append(err[0])
    mu2_up_errors.append(err[1])
    
mu1_up_errors = np.array(mu1_up_errors)   
mu1_down_errors = np.array(mu1_down_errors)   
mu2_up_errors = np.array(mu2_up_errors)  
mu2_down_errors = np.array(mu2_down_errors)   

# Taking average errors of up and down error bar, because propagating asymmetric error bars is non-trivial
average_error_mu1 = mu1_up_errors + mu1_down_errors / 2
print(average_error_mu1)
average_error_mu2 = mu2_up_errors + mu2_down_errors / 2
print(average_error_mu2)

# Computing the actual ratios of the fitted and original relative magnifications
fit_ratios = np.array(mu2_fits) / np.array(mu1_fits)
input_ratios = np.array(mu2_inputs) / np.array(mu1_inputs)

# Computing relevant error
fit_error_ratios = fit_ratios * ((average_error_mu1 / mu1_fits)**2 + (average_error_mu2 / mu2_fits)**2)**0.5 # error is plus or minus this value

# Making plot of fitted and original relative magnification ratios versus A_v.
plt.title(r'Relative ratios of magnifications for different $A_V$ ($\Delta t= 10$)')   #fitted only with salt2 (hsiao data not included in this)
plt.xlabel(r'$A_V$')
plt.ylabel(r'$\mu_2/ \mu_1$ ratios')
for i, A_v in enumerate(np.array(data['A_v'])):
    
    y_err_delays = [[delay_down_errors[i]], [delay_up_errors[i]]] # this needs to be in this format, ie errorbar only accepts form [[down_error], [up_error]] for a given point (or of course [[down_errorS], [up_errorS]])
    y_err_mu1 = [[mu1_down_errors[i]], [mu1_up_errors[i]]]
    y_err_mu2 = [[mu2_down_errors[i]], [mu2_up_errors[i]]]
    
    if np.array(data['Source model'][i]) == 'salt2-extended': 
        if np.array(data["Input magnifications"][i])[0] == 2 and np.array(data["Input magnifications"][i])[1] == 7:
            plt.scatter(A_v, fit_ratios[i], label=f'Fit ratios, [2,7]', c = 'red')
            plt.axhline(y=input_ratios[i], label = f'Input ratio, {np.array(data["Input magnifications"][i])}', c='red', linestyle = '--')
            plt.errorbar(A_v, fit_ratios[i], yerr = fit_error_ratios[i], c = 'red')
        else: 
            plt.scatter(A_v, fit_ratios[i], label=f'Fit ratios, {np.array(data["Input magnifications"][i])}', c = 'green')
            plt.axhline(y=input_ratios[i], label = f'Input ratio, {np.array(data["Input magnifications"][i])}', c = 'green', linestyle = '--')
            plt.errorbar(A_v, fit_ratios[i], yerr = fit_error_ratios[i], c = 'green')
    else:
        pass
    
handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys(), loc='right')   #just copied off stack exchange to stop legend repeating itself
plt.show()




















