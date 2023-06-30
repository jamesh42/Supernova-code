#!/usr/bin/env python
# coding: utf-8

#!/usr/bin/env python
# coding: utf-8

# In[44]:

import sncosmo
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore') # to make output easier to read on jupyter nb's
from astropy.io import ascii
from astropy.table import Table
from test_mockLC_lensed_visual import GLSNe
from test_mockLC_lensed_visual import unresolved_and_individual_lcs

#NOTE: ONLY THING I CHANGED IS I COMMENTED OUT FIRST TWO IMPORTS BECAUSE I COULDN'T GET THESE TO WORK

# By convention, in this simulation, image 2 is held fixed and image 1 is varied over





# Models to use
model1_type = 'salt2-extended'
model2_type = 'salt2-extended'  #for lightcurves.py (not lightcurves_multmodel), these must both be salt2 to work.

# Input and lensing parameters
# Input
z_val = 0.03   # source supernova redshift
t0_val = 0
x1_val = 0.
c_val = 0.
x0_val = 6e-05
# Lensing
nimages = 2  # Multiplicity
source = GLSNe(model1_type, nimages) # choose sn source
dust = sncosmo.CCM89Dust()
delta_t_min = 0
delta_t_max = 15 # preferably difference is even number?
mu_1_min = 1
mu_1_max = 10
mu_2 = 10

# Initialising parameter space to vary over
param_space_dim = 2

delta_t_list = np.linspace(delta_t_min, delta_t_max, param_space_dim)
mu_1_list = np.linspace(mu_1_min, mu_1_max, param_space_dim)

# Initialising models
# Model 1 (input model)
model1 = sncosmo.Model(source, effects=[dust], effect_names=['MW'], effect_frames=['obs']) # sncosmo model
MWebv = 0.03   # what is this? Milky way ebv?
model1.set(**{'MWebv':MWebv})   # what are the **'s?
model1.set(**{'z':z_val})   # zsource?
# Setting model 1 parameters NOT varied over
model1.set(**{'x1': x1_val})
model1.set(**{'c': c_val})
model1.set(**{'x0': x0_val})
model1.set(**{'mu_2':mu_2})
model1.set(**{'t0':t0_val})

model1_num_params = len(model1.param_names)

# Model 2 (output model)
model2 = sncosmo.Model(source=model2_type)
model2.set(z=z_val)  #we know this quantity
model2.set(t0=0)   # this fixes nans, but i dont know why

model2_num_params = len(model2.param_names)

params_to_fit = ['t0', 'x0', 'x1', 'c'] # specifies the parameters to vary in the fitting
# SHOULD I REALLY BE FIXING T0 HERE ETC.? PLUS DUST!

# Misc.

# Creating time list (USELESS RIGHT NOW)
start_time = -15  #days, what should this be?
end_time = 40
sample_freq = 5 #days, to ensure data is realistic and observation taken every ~5 days
# HOW DO I INCORPORATE THIS INTO LENSED SNE CODE?

no_time_points = (end_time - start_time) / sample_freq
no_time_points = int(no_time_points)
time_list = np.linspace(start_time, end_time, no_time_points)

# Some observational params
frac_err = 0.1 #10%?
zp_val = 25. 
data_file_name = 'data.txt'

# Initialising arrays to store input and output parameters

param_in_array = np.zeros((model1_num_params, 0))    #array with input parameters for each loop/iteration in each column

param_out_array = np.zeros((model2_num_params, 0)) 
param_err_array = np.zeros((len(params_to_fit), 0))  #since only the parameters being fitted have associated uncertainties


colour_list = ['red', 'green', 'blue', 'orange', 'purple']  #add more and then trim

#%%

# Iterating over parameter space
for val_delta_t in delta_t_list:
    for mu_1_val in mu_1_list:
        # Step 1: Setting up lensed SN model
        
        # Setting model 1 parameters that ARE varied over

        dt_1 = -val_delta_t/2
        dt_2 = val_delta_t/2
        
        model1.set(**{'mu_1':mu_1_val})   #magnification of the first image
        model1.set(**{'dt_1':dt_1})
        model1.set(**{'dt_2':dt_2})

        # Saving model1 params to array
        param_in_array = np.concatenate((param_in_array, model1.parameters.reshape(-1,1)), axis=1)   #reshape makes it a column instead of a row, which is necessary to concatenate it

        #tspace = 'discrete'
        tspace = 'continuous'   # I don't know what this means
        
        # Generating lensed lightcurve using Suhail's code
        lc = unresolved_and_individual_lcs(model1, savelc=True, nimages = nimages, tgrid=tspace)
        ascii.write(lc, data_file_name, overwrite=True)    #writing this lensed lc info to file to be read and fitted to by model 2
        data = sncosmo.read_lc(data_file_name)   # producing lensed lc data
        
        # Step 2: Taking produced data and using SALT2 model to produce line of best fit
        
        # Run the fit, using pre-defined model2
        try:
            result, fitted_model = sncosmo.fit_lc(
                data, model2,
                params_to_fit, modelcov = True)  ## parameters of model to vary + bounds on parameters (if any)
            
            # Collecting output values for later residuals
            param_out_array = np.concatenate((param_out_array, fitted_model.parameters.reshape(-1,1)), axis=1)

            # Errors
            errors = result.errors.values()   # extracting error values from error dictionary
            errors = np.fromiter(errors, dtype=float)    # making into numpy array
            param_err_array = np.concatenate((param_err_array, errors.reshape(-1,1)), axis=1) # storing
            
            # Printing fit data (more in sncosmo fitting on website)
            #print("chi^2 value at minimum:", result.chisq)

            # Plotting fitted lightcurve
            #fig2_name = 'output_fig_x1' + str(round(val_x1,2)) + '_c' + str(round(val_c,2)) + '.png'   #params produced by data
            #fig2 = sncosmo.plot_lc(data, model=fitted_model, errors=result.errors)   #wrapper for matplotlib.pyplot as plt
            #fig2.savefig('Lensed Lightcurve folder/' + fig2_name)
        
        except Exception as e:
            print("Error for delta_t = " + str(val_delta_t) + ' and mu_1 = ' + str(mu_1_val))
            print(e)
            
        finally:
            continue
        
#%%

# Code for comparing input to output parameters
delta_t_list = param_in_array[4] - param_in_array[2]    #dt2 - dt1
mu_1_list = param_in_array[3]

# Computing residuals - this could be made shorter I'm sure, but would involve some tricky formatting loops
z_res = param_out_array[0] - param_in_array[0] #z residuals
t0_res = param_out_array[1] - param_in_array[1]
x0_res = param_out_array[2] - param_in_array[6] #havent tested all of these indices btw
x1_res = param_out_array[3] - param_in_array[7]
c_res = param_out_array[4] - param_in_array[8]

# Plotting residuals
fig, (ax1, ax2) = plt.subplots(1,2, figsize = (10,3.5))

# Plotting x1 residuals
ax1.set_title(r'$x_1$ input vs output residuals with varying time delay $\Delta t$')
ax1.set_xlabel(r'$\Delta t$')
ax1.set_ylabel(r'$x_1$ residual ($x_1$(output) - $x_1$(input))')
ax1.errorbar(delta_t_list, x1_res, yerr=param_err_array[2], fmt='o')

# Plotting c residuals
ax2.set_title(r'$x_1$ input vs output residuals with varying magnification of image 1')
ax2.set_xlabel(r'$\mu_1$')
ax2.set_ylabel(r'$x_1$ residual ($x_1$(output) - $c$(input))')
ax2.errorbar(mu_1_list, x1_res, yerr=param_err_array[2], fmt='o')

plt.tight_layout(pad=4)   # stops plots from overlapping somehow

fig.savefig('Plot folder/lensingresiduals.png', dpi = 2000)   # change this



















# In[51]:

# New questions:
"""
- seriously get to the bottom of why the residuals would ever not be zero (besides imperfect sampling)

- why does fixing t0 fix everything? and is that even fair, when the point of the project is to calculate
the time delay (i.e t0 for the different lightcurves)?

"""


# In[52]:


# misc comments
"""
- why doesnt reducing sample from 5 days to 1 day fix divergence in input output parameters? or, more generally,
why do plots disagree at all?


-should i have a silly amount of plots like I do? Do I really need to record all of them as I have done?

# just worry about salt2 now, come back to other models after lensing done

"""

#%%
"""
Answered questions:
    - what magnitude system should we use (ab, vega, etc.) and why? just use ab cos convention
    
    - what is zp? why would you use it? just so you can use your own arbitary non-built in unit system? zero-point, useful for comparing two different telescopes with different y-axis constants, e.g. ground based to space nased, you need to accoumt for 
    different instruments/telescopes

    - does fitting fit using every filter or just one? SED has all wavelengths so fits all at same time

"""
