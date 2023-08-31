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

# Time delay list
delta_t_value = 10
mag_list = [mu_1,mu_2]
delta_t_list =[10, 15, 20]
# delta_t_list = [15]
# delta_t_list = [10]

# Creating functions for sn_parameters
def x1_func():
    return x1_val

def c_func():
    return c_val

sn_parameters = {'x1':x1_func,'c':c_func}   # i assume you don't need to specify z because specified in other argument?    
    
# Generating lensed lightcurve model using Suhail's code
source = GLSNe('salt2-extended', nimages) # choose sn source
model = sncosmo.Model(source) # sncosmo model
model.set(**{'z':z_val})   # zsource?
# Setting model 1 parameters NOT varied over
model.set(**{'x1': x1_val})
model.set(**{'c': c_val})
model.set(**{'x0': x0_val})
model.set(**{'mu_2':mu_2})
model.set(**{'t0':t0_val})
model.set(**{'mu_1':mu_1})   #magnification of the first image

# Misc. lists and stuff
end_delay_residual_list = []   # for storing time delay residuals
bounds_list = [True, False]   #whether or not to use unresolved fit as priors
prior_bounds_used = [] # for keeping track of whether prior bounds were used in a given loop,1 if yes, 0 if no
dt_list = []
prior_bounds_list = [] # for saving prior bounds themselves to make sure they're sensible
td_error_list = []
time_info_list = []
fitted_param_list = []

# Fitting parameter parameters, commented out because code breaks when I soft code these into the fitting functions
# params_to_fit = ['x0', 't0', 'dt_1', 'dt_2', 'x1', 'c', 'mu_1', 'mu_2'] # removed mu's since need to understand lens potential for time delays to be useful anyway
# default_bounds = {'t0':(-25,25), 'dt_1':(-30,0), 'dt_2':(0,30), 'x1': (-3,3), 'c': (-1,1), 'mu_1':(1,4), 'mu_2':(5,8)}
# params_to_fit = ['x0', 't0', 'dt_1', 'dt_2', 'x1', 'c']
# #default_bounds = {'t0':(-25,25), 'dt_1':(-30,0), 'dt_2':(0,30), 'x1': (-3,3), 'c': (-1,1)}
# constants_in_fits = {'z':z_val, 'mu_1': mu_1, 'mu_2': mu_2}

#%%

# Creating functions for sn_parameters
def x1_func():
    return x1_val

def c_func():
    return c_val

sn_parameters = {'x1':x1_func,'c':c_func}   # i assume you don't need to specify z because specified in other argument?    
    

x1_bound_list = [0.2,0.1,0.05]
x1_bound_append_list = []
#x1_bound_list = [0.2]

time_lc1_list=[]
time_lc2_list=[]
flux1_list=[]
flux2_list=[]
flux1_err_list=[]
flux2_err_list=[]

for bounds in bounds_list:
    for i, dt in enumerate(delta_t_list):
        for j, x1_bound in enumerate(x1_bound_list):
            
            # To avoid exceptions
            max_retries = 15
            retry_count = 0
                   
            while retry_count < max_retries:
                try:    
        
                    start = time.time()            
        
                    dt_1 = -dt/2  # Adding 15 to each of these did not solve the problem at all
                    dt_2 = dt/2
                    
                    # |||| Creating resolved image data ||||
                
                    # Creating zp list in format sntd.createMISN wants
                    zp_list = [25.] * len(band_list_ztf)
                
                    data_file_name1 = f'td{dt}days-params1-img1'
                    output_file_path1 = f'Lightcurve data/{data_file_name1}'
                    
                    data_file_name2 = f'td{dt}days-params1-img2'
                    output_file_path2 = f'Lightcurve data/{data_file_name2}'

                    # Open the file in binary mode and load the data using pickle
                    with open(output_file_path1, 'rb') as file:
                        image_1_data = pickle.load(file)
                        
                    with open(output_file_path2, 'rb') as file:
                        image_2_data = pickle.load(file)
                
                    # Saving data to astropy tables - this is our created data
                    # image_1_dat = myMISN.images['image_1']['table']
                    # image_2_dat = myMISN.images['image_2']['table']
                    
                    # Creating MISN instance using data
                    resolved_MISN=sntd.table_factory([image_1_data,image_2_data],telescopename='telescope',object_name='example_SN')
                    ztfg_indices_image_1 = np.where(np.array(resolved_MISN.images['image_1']['__dict__']['table']['band']) == 'ztfg')
                    ztfg_indices_image_2 = np.where(np.array(resolved_MISN.images['image_2']['__dict__']['table']['band']) == 'ztfg')
                    print(ztfg_indices_image_1)
                    print(resolved_MISN.images['image_1']['__dict__']['table'])
                    fluxes1 = np.array(resolved_MISN.images['image_1']['__dict__']['table']['flux'][ztfg_indices_image_1])
                    flux1_err = np.array(resolved_MISN.images['image_1']['__dict__']['table']['fluxerr'][ztfg_indices_image_1])
                    flux2_err = np.array(resolved_MISN.images['image_2']['__dict__']['table']['fluxerr'][ztfg_indices_image_2])
                    fluxes2 = np.array(resolved_MISN.images['image_2']['__dict__']['table']['flux'][ztfg_indices_image_2])
                    time_lc1 = np.array(resolved_MISN.images['image_1']['__dict__']['table']['time'][ztfg_indices_image_1])
                    time_lc2 = np.array(resolved_MISN.images['image_2']['__dict__']['table']['time'][ztfg_indices_image_2])
                    
                                
                    time_lc1_list.append(time_lc1)
                    time_lc2_list.append(time_lc2)
                    flux1_list.append(fluxes1)
                    flux2_list.append(fluxes2)
                    flux1_err_list.append(flux1_err)
                    flux2_err_list.append(flux2_err)
                    
                    # |||| Creating unresolved image data ||||
            
                    model.set(**{'dt_1':dt_1})
                    model.set(**{'dt_2':dt_2})
                    
                    tspace = 'discrete'
                    lc = unresolved_and_individual_lcs(model, savelc=True, nimages = nimages, tgrid=tspace, bands = band_list_ztf, sample_freq = cad)   #keep this on full, ie simulate all and then only use some
                    data_file_name = 'data.txt'
                    ascii.write(lc, data_file_name, overwrite=True)    #writing this lensed lc info to file to be read and fitted to by model 2
                    data = sncosmo.read_lc(data_file_name)   # producing lensed lc data
                    
                    # Setting up model2 (fitting model) using SNTD
                    image_a = sncosmo.Model('salt2-extended')
                    image_b = sncosmo.Model('salt2-extended')
                    unresolved = sntd.unresolvedMISN([image_a,image_b])
                    
                    # Initialising MISN object using unresolved data
                    unresolved_MISN=sntd.table_factory(data,telescopename='Unknown',object_name='Unresolved')   # basically just creates MISN instance, but in a clever way to use your given data
                
                    data_creation_time = time.time()
                
                    # |||| Fitting unresolved data, and saving fitted params as priors (bounds) ||||
                    
                    # Note how we're setting unresolved priors to actual value, this is kind of silly perhaps and defeats the point?
                    unresolved.set_delays([dt_1,dt_2])   # So does that (comment below) mean that these are essentially just priors?
                    unresolved.set_magnifications([mu_1,mu_2])   # I've noticed if you comment out either of these lines, it breaks. So I think you need to provide this info for it to work.
                    unresolved.set(z=z_val)
                    unresolved.set(t0=t0_val)
                    unresolved.set(x0= x0_val)
    
                    # Fitting unresolved lightcurve
                    unresolved_fitCurves=sntd.fit_data(unresolved_MISN,snType='Ia', models=unresolved,bands=band_list_ztf,   # there is another filter ztfg but not shown here because it breaks the code
                                    params=['x0', 't0', 'dt_1', 'dt_2', 'x1', 'c'],constants={'z':z_val, 'mu_1': mu_1, 'mu_2': mu_2},   # changing order of params entered here doesnt affect fitcurves.model.parameters order
                                    bounds={'t0':(-25,25), 'dt_1':(-30,10), 'dt_2':(-10,30), 'x1': (-3,3), 'c': (-1,1)}, method='parallel',npoints=100)
                                    # ||| THESE MUST BE THE SAME AS FOR NON-PRIOR FITTING BOUNDS, CANNOT PUT IN default_bounds bc this breaks it for some reason
                    # Saving results as priors
                   #print(list(zip(unresolved_fitCurves.images['image_1'].fits.model.param_names, unresolved_fitCurves.images['image_1'].fits.model.parameters)))
                    fitted_params = unresolved_fitCurves.images['image_1'].fits.model.parameters   # image 1?? I think image 1 is correct, since there is only one image because it is resolved!?
                    param_names = unresolved_fitCurves.images['image_1'].fits.model.param_names
                    #param_names = unresolved_fitCurves.images['image_2']  think this breaks it because image2 doesnt exist
                    # print('FITTED PARAMS:')
                    # print(fitted_params)
                    # print('Param names')
                    # print(param_names)
                    #hard code it then, will not work in general very well
                    x0 = fitted_params[2]
                    x1 = fitted_params[3]
                    c = fitted_params[4]
                    mu_1 = fitted_params[6]
                    dt_1_unres = fitted_params[5]
                    dt_2_unres = fitted_params[7]
                    mu_2 = fitted_params[8]
                    
                    fitted_param_list.append(fitted_params)
                    
                    bt = 0.01 #'bounds tightness'  # smaller is tighter
                    
                
                    # resolved_prior_bounds = {'t0': (-25,25), 'x1': (x1-bt*x1 - 0.01, x1+bt*x1 + 0.01), \
                    #                           'c': (c-bt*c - 0.01, c+bt*c + 0.01),\
                    #                               'dt_1': (dt_1_unres-bt*dt_1_unres, dt_1_unres+bt*dt_1_unres), 'dt_2': (dt_2_unres-bt*dt_2_unres, dt_2_unres+bt*dt_2_unres)}    
                        
                    #bounds2 = {'t0':(-25,25), 'dt_1':(-30,10), 'dt_2':(-10,30), 'x1': (-x1_bound,x1_bound), 'c': (-1,1)}
                        
                    # print(resolved_prior_bounds)
                    
                    # prior_bounds_list.append(resolved_prior_bounds)
                    
                    fit_unresolved_time = time.time()
                
                    # Fitting RESOLVED data using unresolved results as PRIOR BOUNDS!
                    
                    #default_bounds = {'t0':(-25,25), 'dt_1':(-30,10), 'dt_2':(-10,30), 'x1': (-3,3), 'c': (-1,1)}
                      
                    #resolved_prior_bounds  = {'t0': (-25,25), 'x1': (x1-bt*x1 - 0.01, x1+bt*x1 + 0.01), \
                                              # 'c': (c-bt*c - 0.01, c+bt*c + 0.01)}
                    c_bound = x1_bound
                    
                    if bounds == True:
                        resolved_fitCurves=sntd.fit_data(resolved_MISN,snType='Ia', models='salt2-extended',bands=band_list_ztf,
                                        #params=params_to_fit,constants=constants_in_fits,
                                        params=['x0', 't0', 'dt_1', 'dt_2', 'x1', 'c'],constants={'z':z_val, 'mu_1': mu_1, 'mu_2': mu_2},   # changing order of params entered here doesnt affect fitcurves.model.parameters order
                                        bounds={'t0':(-25,25), 'dt_1':(-30,10), 'dt_2':(-10,30), 'x1': (-x1_bound,x1_bound), 'c': (-1,1)}, method='parallel',npoints=100)
                                        #bounds={'t0':(-25,25), 'dt_1':(-10,0), 'dt_2':(0,10), 'x1': (-3,3), 'c': (-1,1)}, method='parallel',npoints=100)
                        print(bounds)
                        resolved_fitCurves.plot_object(bands='ztfg', showFit=True)
                        plt.show()
                        prior_bounds_used.append(1)
                        
                    else: 
                        resolved_fitCurves=sntd.fit_data(resolved_MISN,snType='Ia', models='salt2-extended',bands=band_list_ztf,
                                               params=['x0', 't0', 'dt_1', 'dt_2', 'x1', 'c'],constants={'z':z_val, 'mu_1': mu_1, 'mu_2': mu_2},   # changing order of params entered here doesnt affect fitcurves.model.parameters order
                                               bounds={'t0':(-25,25), 'dt_1':(-30,10), 'dt_2':(-10,30), 'x1': (-3,3), 'c': (-1,1)}, method='parallel',npoints=100)
                                               # ||| THESE MUST BE THE SAME AS FOR NON-PRIOR FITTING BOUNDS, CANNOT PUT IN default_bounds bc this breaks it for some reason
                        resolved_fitCurves.plot_object(bands='ztfg', showFit=True)
                        plt.show()
                        prior_bounds_used.append(0)
                        
                    res_time_delay_SNTD = resolved_fitCurves.parallel.time_delays['image_2'] - resolved_fitCurves.parallel.time_delays['image_1']
                    delay_residual = res_time_delay_SNTD - dt
                    
                    end_delay_residual_list.append(delay_residual)
                    print(delay_residual)
                    dt_list.append(dt)
                    x1_bound_append_list.append(x1_bound)
                    
                    del_errors = resolved_fitCurves.parallel.time_delay_errors
                    down_err = abs(del_errors['image_2'][0])
                    up_err = del_errors['image_2'][1]
    
                    delay_error = [down_err, up_err]
                    td_error_list.append(delay_error)
                    
                    end_time = time.time()
                    
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
                
            time_create =  data_creation_time - start
            time_fit_unresolved = fit_unresolved_time - data_creation_time
            time_fit_resolved = end_time - fit_unresolved_time
            time_info = f'Data creation took {time_create}s, fitting unresolved took {time_fit_unresolved}, fitting resolved took {time_fit_resolved}'
            print(time_info)
            time_info_list.append(time_info)

            
#%%
            
# For plotting lightcurves to see if they line up

# colour_list = ['blue', 'red', 'green', 'pink']

# print(len(flux1_list))
# print(len(flux2_list))

# plt.title('Lightcurves for two consecutive runs (ZTFG)')
# plt.ylabel('Flux')
# plt.xlabel('Observer frame time (days')
# for i, fluxes1 in enumerate(flux1_list):
#     plt.scatter(time_lc1_list[i], fluxes1, color = colour_list[i], label = 'Image 1')
#     plt.errorbar(time_lc1_list[i], fluxes1, yerr = flux1_err_list[i], color = colour_list[i], linestyle='')
#     plt.scatter(time_lc2_list[i], flux2_list[i], color = colour_list[i], label = 'Image 2', marker = '*')
#     plt.errorbar(time_lc2_list[i], flux2_list[i], yerr = flux2_err_list[i], linestyle='')    
#     print('bbb')

# plt.scatter(time_lc1_list[0], flux1_list[0], color = colour_list[0], label = 'Image 1', marker = 's')
# plt.errorbar(time_lc1_list[0], flux1_list[0], yerr = flux1_err_list[0], color = colour_list[0], linestyle='')
# plt.scatter(time_lc2_list[0], flux2_list[0], color = colour_list[0], label = 'Image 2', marker = 's')
# plt.errorbar(time_lc2_list[0], flux2_list[0], yerr = flux2_err_list[0], linestyle='')    

# handles, labels = plt.gca().get_legend_handles_labels()
# by_label = dict(zip(labels, handles)) 
# plt.legend(by_label.values(), by_label.keys(), loc='upper right')
# # plt.legend()
# plt.show()

#%%

delay_up_errors = []
delay_down_errors = []

# print(td_error_list)

for i, err in enumerate(td_error_list): #for each data point
    delay_down_errors.append(err[0])
    delay_up_errors.append(err[1])
    
print(prior_bounds_used)
    
plt.title('Time delay error versus x1 bound')
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
    elif prior_bounds_used[i] == 0:
    # elif prior_bounds_used[i] == 0 and x1_bound_append_list[i] == 0.2:
        if dt == 10:
            plt.scatter(3, np.array(delay_down_errors[i]) + np.array(delay_up_errors[i]), color = 'cyan', label = 'No prior, 10 days')
        if dt == 15:
            plt.scatter(3, np.array(delay_down_errors[i]) + np.array(delay_up_errors[i]), color = 'orange', label = 'No prior, 15 days')
        if dt == 20:
            plt.scatter(3, np.array(delay_down_errors[i]) + np.array(delay_up_errors[i]), color = 'brown', label = 'No prior, 20 days')

#plt.scatter(x1_bound_append_list, np.array(delay_down_errors) + np.array(delay_up_errors))
handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles)) 
plt.legend(by_label.values(), by_label.keys(), loc='center')
plt.show()
    
print(delay_down_errors)
color_list1 = ['blue', 'green', 'pink']
color_list2 = ['red', 'brown', 'orange']

label_priors = ['No priors', 'Priors']

L = int(0.5*len(dt_list))
plt.title('Time delay residuals using x1 priors')
plt.xlabel('Time delay (days)')
plt.ylabel('Time delay residuals (days)')
for i, dt in enumerate(dt_list):
    if prior_bounds_used[i] == 1:
        if x1_bound_append_list[i] == 0.05:
            #plt.scatter(dt, end_delay_residual_list[i], label = f'{label_priors[prior_bounds_used[i]]}, x1 bound = {x1_bound_append_list[i]}', color = 'blue')
            plt.scatter(dt, end_delay_residual_list[i], label = r'Priors, x1 [$\pm$ 0.05]', color = 'blue')
            plt.errorbar(dt, end_delay_residual_list[i], yerr = [[delay_down_errors[i]], [delay_up_errors[i]]], linestyle='', color = 'blue')
        elif x1_bound_append_list[i] == 0.1:
            plt.scatter(dt, end_delay_residual_list[i], label = 'Priors, x1 [$\pm$ 0.1]', color = 'red')
            plt.errorbar(dt, end_delay_residual_list[i], yerr = [[delay_down_errors[i]], [delay_up_errors[i]]], linestyle='', color = 'red')
        else:
            plt.scatter(dt, end_delay_residual_list[i], label = 'Priors, x1 [$\pm$ 0.2]', color = 'green')
            plt.errorbar(dt, end_delay_residual_list[i], yerr = [[delay_down_errors[i]], [delay_up_errors[i]]], linestyle='', color = 'green')
    if prior_bounds_used[i] == 0:
        if x1_bound_append_list[i] == 0.05:   # but it shouldnt matter because priors are not being used in this (==0), although as it turns out it does matter because of, i guess, intrinsic variation in sim
            plt.scatter(dt, end_delay_residual_list[i], label = r'No priors', color = 'brown')
            plt.errorbar(dt, end_delay_residual_list[i], yerr = [[delay_down_errors[i]], [delay_up_errors[i]]], linestyle='', color = 'brown')
        # elif x1_bound_append_list[i] == 0.1:
        #     plt.scatter(dt, end_delay_residual_list[i], label = 'No priors, x1 [$\pm$ 0.1]', color = 'red', marker = '*')
        #     plt.errorbar(dt, end_delay_residual_list[i], yerr = [[delay_down_errors[i]], [delay_up_errors[i]]], linestyle='', color = 'red')
        # else:
        #     plt.scatter(dt, end_delay_residual_list[i], label = 'No priors, x1 [$\pm$ 0.2]', color = 'green', marker = '*')
        #     plt.errorbar(dt, end_delay_residual_list[i], yerr = [[delay_down_errors[i]], [delay_up_errors[i]]], linestyle='', color = 'green')
            
handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles)) 
plt.legend(by_label.values(), by_label.keys(), loc='upper left')
# plt.legend()
plt.show()




