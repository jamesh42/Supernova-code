#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 16:52:32 2023

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
from test_mockLC_lensed_visual import GLSNe
from test_mockLC_lensed_visual import unresolved_and_individual_lcs
import sntd


#%%

# Generating unresolved data using SD's code

model1_type = 'salt2-extended'

# Input and lensing parameters
# Input
z_val = 0.3544   # source supernova redshift
t0_val = 0
x1_val = 0.
c_val = 0.
x0_val = 6e-05
# Lensing
nimages = 2  # Multiplicity
source = GLSNe(model1_type, nimages) # choose sn source
#dust = sncosmo.CCM89Dust()
delta_t = 20
min_delta_t = 10
max_delta_t = 30
mu_1 = 2
mu_2 = 7

delta_t_list = np.linspace(min_delta_t, max_delta_t, 4)


# param_out_array = np.zeros((len(out_param_vals), 0)) 
# param_in_array = np.zeros((len(out_param_vals), 0)) 

param_out_array = np.zeros((9, 0)) 
param_in_array = np.zeros((9, 0))

#%%

# Fitting and plotting data using SNTD
for dt_val in delta_t_list:
    # Creating lensing
    # Model 1 (input model)
    #model1 = sncosmo.Model(source, effects=[dust], effect_names=['MW'], effect_frames=['obs']) # sncosmo model
    model = sncosmo.Model(source) # sncosmo model
    MWebv = 0.03   # what is this? Milky way ebv?
    #model1.set(**{'MWebv':MWebv})   # what are the **'s?
    model.set(**{'z':z_val})   # zsource?
    # Setting model 1 parameters NOT varied over
    model.set(**{'x1': x1_val})
    model.set(**{'c': c_val})
    model.set(**{'x0': x0_val})
    model.set(**{'mu_2':mu_2})
    model.set(**{'t0':t0_val})
    model.set(**{'mu_1':mu_1})   #magnification of the first image

    data_file_name = 'data.txt'

    dt_1 = -dt_val/2
    dt_2 = dt_val/2

    model.set(**{'dt_1':dt_1})
    model.set(**{'dt_2':dt_2})

    tspace = 'discrete'   # I CHANGED SD'S CODE SO THAT BY DEFAULT THIS WAS 5 DAYS

    # Generating lensed lightcurve using Suhail's code
    lc = unresolved_and_individual_lcs(model, savelc=True, nimages = nimages, tgrid=tspace, sample_freq = 5)   # I modified code to build in sample freq = 5 days and flux error = 0.1 (10%)
    ascii.write(lc, data_file_name, overwrite=True)    #writing this lensed lc info to file to be read and fitted to by model 2
    data = sncosmo.read_lc(data_file_name)   # producing lensed lc data

    # Setting up model2 (fitting model) using SNTD
    image_a = sncosmo.Model('salt2-extended')
    image_b = sncosmo.Model('salt2-extended')
    unresolved=sntd.unresolvedMISN([image_a,image_b])
    
    # Setting parameters
    unresolved.set_delays([-5,20])   # So does that (comment below) mean that these are essentially just priors?
    unresolved.set_magnifications([10,1])   # I've noticed if you comment out either of these lines, it breaks. So I think you need to provide this info for it to work.
    unresolved.set(z=z_val)
    unresolved.set(t0=t0_val)
    
    # Initialising MISN object using unresolved data
    new_MISN=sntd.table_factory(data,telescopename='Unknown',object_name='Unresolved')   # basically just creates MISN instance, but in a clever way to use your given data
    
    # Fitting
    fitCurves=sntd.fit_data(new_MISN,snType='Ia', models=unresolved,bands=['ztfr', 'ztfi'],   # there is another filter ztfg but not shown here because it breaks the code
                    params=['x0','x1','t0','c','mu_1', 'dt_2','mu_2', 'dt_1'],constants={'z':z_val},   # changing order of params entered here doesnt affect fitcurves.model.parameters order
                    bounds={'t0':(-25,25),'x1':(-2,2),'c':(-1,1), 'mu_1': (1,3), 'dt_2': (10,30), 'dt_1': (-10,-30),'mu_2': (6,8)})
    
    print(list(zip(fitCurves.images['image_1'].fits.model.param_names, fitCurves.images['image_1'].fits.model.parameters)))
    
    # Collecting output parameters
    param_out_dict = list(zip(fitCurves.images['image_1'].fits.model.param_names, fitCurves.images['image_1'].fits.model.parameters))   # I don't understand this code so am just copying it from readthedocs

    out_param_vals = []
    for i in range(len(param_out_dict)):   #collecting output values
        out_param_vals.append(param_out_dict[i][1])
        
    out_param_names = []
    for i in range(len(param_out_dict)):    # collecting output names
        out_param_names.append(param_out_dict[i][0])
        
    # Collecting input parameters
    in_param_vals = model.parameters
    in_param_names = model.param_names
    
    # Code to ensure that SD input parameters are in same order as SNTD's fitting parameters, so that they can be compared
    for i in range(len(out_param_names)):   # setting SD's list to be like SNTD's
        if in_param_names[i] == out_param_names[i]:
            continue
        elif in_param_names[i] != out_param_names[i]:
            to_swap_index = in_param_names.index(str(out_param_names[i]))    #index in param_in_param_names of corresponding ith component of out_param_names
            in_param_names[i], in_param_names[to_swap_index] = in_param_names[to_swap_index], in_param_names[i]   # swapping param names
            in_param_vals[i], in_param_vals[to_swap_index] = in_param_vals[to_swap_index], in_param_vals[i]   # swapping the actual values! (which was the point of this exercise)
        
    # Or it was never broken in the first place? So the problem for some reason arises only when
    # you run the whole code. If you run just this cell having already run the whole code, the issue
    # somehow does not exist. I don't see why, but I don't think it's too important. The point is, either way, my swap code fixes it
        

    # Making new in and out vals numpy arrays
    out_param_vals = np.array(out_param_vals)
    in_param_vals = np.array(in_param_vals)
    

    param_in_array = np.concatenate((param_in_array, in_param_vals.reshape(-1,1)), axis=1) 
    
    param_out_array = np.concatenate((param_out_array, out_param_vals.reshape(-1,1)), axis=1) 
    
    input_overlay = in_param_vals[1:]
    print(input_overlay)
    
    # Plotting
    fig1 = fitCurves.plot_object(showFit=True,plot_unresolved=True)
    fig2 = fitCurves.plot_fit()    # gives corner plot of probabilities

    # Extract the axes
    ndim  = 8 # must be equal to number of graphs    
    axes = np.array(fig2.axes).reshape((ndim, ndim))
    
    # Loop over the diagonal
    for i in range(ndim):
        ax = axes[i, i]
        ax.axvline(input_overlay[i], color="g")
    
    # Loop over the histograms
    for yi in range(ndim):
        for xi in range(yi):
            ax = axes[yi, xi]
            ax.axvline(input_overlay[xi], color="g")
            ax.axhline(input_overlay[yi], color="g")
            ax.plot(input_overlay[xi], input_overlay[yi], "sg")
    
    #plt.show()
    fig1.savefig('double_lc.png')
    fig2.savefig('cornerdouble.png')



#%% 

# Residuals
# Computing residuals - this could be made shorter I'm sure, but would involve some tricky formatting loops
z_res = abs(param_out_array[0] - param_in_array[0]) #z residuals
t0_res = abs(param_out_array[1] - param_in_array[1])
x0_res = abs(param_out_array[2] - param_in_array[2]) #havent tested all of these indices btw
x1_res = abs(param_out_array[3] - param_in_array[3])
dt1_res = abs(param_out_array[4] - param_in_array[4])
mu1_res = abs(param_out_array[5] - param_in_array[5])
dt2_res = abs(param_out_array[6] - param_in_array[6])
mu2_res = abs(param_out_array[7] - param_in_array[7])

print(param_in_array)
print(param_out_array)

for i in range(len(out_param_vals)):
    plt.title(str(in_param_names[i] + ' residuals'))
    plt.xlabel(r'$\Delta t$')
    plt.ylabel(str(in_param_names[i] + ' residuals'))
    plt.scatter(delta_t_list, abs(param_out_array[i] - param_in_array[i]))
    #plt.colorbar(label=r'$\mu_1$')
    plt.plot()
    plt.show()
    #plt.savefig('Plot folder/x_1vsdt.png', dpi=2000)



#chi = fitCurves.series_chisq

# Is this only image 1? I suppose the plot says only image 1 is fitted, and so I'd imagine that this is really for the whole fit
chi = fitCurves.images['image_1'].chisq
reduced_chi = fitCurves.images['image_1'].reduced_chisq    # THIS IS HOW YOU OBTAIN IT. 
#FOR SOME REASON, fitCurves.series_chisq doesn't work, even though series_chisq IS a property of the MISN class too
print(chi)
print(reduced_chi)
# print('x')



"""
Questions:
    
- do the 'priors' matter? Comparison video on phone from set_delays (-5,5) to (-5,20) and set_magnifications from i think (2,5) to (10,1)
 - CLEARLY THEY ARE SOMEWHAT SENSITIVE TO IT

"""




