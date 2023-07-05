#!/usr/bin/env python
# coding: utf-8

#!/usr/bin/env python
# coding: utf-8

# In[44]:

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

#NOTE: ONLY THING I CHANGED IS I COMMENTED OUT FIRST TWO IMPORTS BECAUSE I COULDN'T GET THESE TO WORK

# By convention, in this simulation, image 2 is held fixed and image 1 is varied over





# Models to use
model1_type = 'salt2-extended'
model2_type = 'salt2-extended'  #for lightcurves.py (not lightcurves_multmodel), these must both be salt2 to work.

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
delta_t_min = 0
delta_t_max = 15 # preferably difference is even number?
mu_1_min = 1
mu_1_min = 0
mu_1_max = 10
mu_2 = 5

# Initialising parameter space to vary over
param_space_dim = 40

delta_t_list = np.linspace(delta_t_min, delta_t_max, param_space_dim)
mu_1_list = np.linspace(mu_1_min, mu_1_max, param_space_dim)

# Initialising models
# Model 1 (input model)
#model1 = sncosmo.Model(source, effects=[dust], effect_names=['MW'], effect_frames=['obs']) # sncosmo model
model1 = sncosmo.Model(source) # sncosmo model
MWebv = 0.03   # what is this? Milky way ebv?
#model1.set(**{'MWebv':MWebv})   # what are the **'s?
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
#DUST?

# Misc.

# Some observational params
frac_err = 0.1 #10%?
zp_val = 25. 
data_file_name = 'data.txt'

# Initialising arrays to store input and output parameters

param_in_array = np.zeros((model1_num_params, 0))    #array with input parameters for each loop/iteration in each column

param_out_array = np.zeros((model2_num_params, 0)) 
param_err_array = np.zeros((len(params_to_fit), 0))  #since only the parameters being fitted have associated uncertainties

# Storing chisq etc.
extra_info_to_store = 2
info_array = np.zeros((extra_info_to_store, 0))

loop_no = 0

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

        tspace = 'discrete'   # I CHANGED SD'S CODE SO THAT BY DEFAULT THIS WAS 5 DAYS

        # Generating lensed lightcurve using Suhail's code
        lc = unresolved_and_individual_lcs(model1, savelc=True, nimages = nimages, tgrid=tspace)   # I modified code to build in sample freq = 5 days and flux error = 0.1 (10%)
        ascii.write(lc, data_file_name, overwrite=True)    #writing this lensed lc info to file to be read and fitted to by model 2
        data = sncosmo.read_lc(data_file_name)   # producing lensed lc data
        
        # Step 2: Taking produced data and using SALT2 model to produce line of best fit
        
        loop_no += 1
        print(loop_no)
        
        # Run the fit, using pre-defined model2
        try:
            result, fitted_model = sncosmo.fit_lc(
                data, model2,
                params_to_fit, modelcov = True)  ## parameters of model to vary + bounds on parameters (if any)
            
            # Saving model1 params to array
            param_in_array = np.concatenate((param_in_array, model1.parameters.reshape(-1,1)), axis=1)   #reshape makes it a column instead of a row, which is necessary to concatenate it
            
            # Collecting model2 output values for later residuals
            param_out_array = np.concatenate((param_out_array, fitted_model.parameters.reshape(-1,1)), axis=1)

            # Storing dof and chisq
            info_to_store = [result.chisq, result.ndof]
            info_to_store = np.array(info_to_store)
            info_array = np.concatenate((info_array, info_to_store.reshape(-1,1)), axis=1)   # so info array has result.chisq in 1st row and result.dnof in 2nd

            # Errors
            errors = result.errors.values()   # extracting error values from error dictionary
            errors = np.fromiter(errors, dtype=float)    # making into numpy array
            param_err_array = np.concatenate((param_err_array, errors.reshape(-1,1)), axis=1) # storing

            # Plotting fitted lightcurve
            #fig2_name = 'output_fig_x1' + str(round(val_x1,2)) + '_c' + str(round(val_c,2)) + '.png'   #params produced by data
            #fig2_name = 'test'
            # fig2 = sncosmo.plot_lc(data, model=fitted_model, errors=result.errors)   #wrapper for matplotlib.pyplot as plt
            # fig2.savefig('Lensed Lightcurve folder/' + fig2_name)
            # plt.show()]
            
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
z_res = abs(param_out_array[0] - param_in_array[0]) #z residuals
t0_res = abs(param_out_array[1] - param_in_array[1])
x0_res = abs(param_out_array[2] - param_in_array[6]) #havent tested all of these indices btw
x1_res = abs(param_out_array[3] - param_in_array[7])
c_res = abs(param_out_array[4] - param_in_array[8])

c_out = param_out_array[4]

plt.show()   # getting rid of any graphs from above

# Plotting residuals
#fig, (ax1, ax2) = plt.subplots(1,2, figsize = (10,3.5))

# Plotting x1 residuals
plt.title(r'$x_1$ input vs output residuals with time delay $\Delta t$')
plt.xlabel(r'$\Delta t$')
plt.ylabel(r'$x_1$ residual')
#plt.errorbar(delta_t_list, x1_res, yerr=param_err_array[2], fmt='.')
plt.scatter(delta_t_list, x1_res, c = mu_1_list, cmap = 'Reds')
plt.colorbar(label=r'$\mu_1$')
plt.savefig('Plot folder/x_1vsdt.png', dpi=2000)

plt.show()

plt.title(r'$x_1$ input vs output residuals with magnification of image 1')
plt.xlabel(r'$\mu_1$')
plt.ylabel(r'$x_1$ residual')
#plt.errorbar(mu_1_list, x1_res, yerr=param_err_array[2], fmt='.')
plt.scatter(mu_1_list, x1_res, c = delta_t_list, cmap = 'Greens') #which works when we have defined dtlist st it is in order from smallest to largest value - actually works by default
plt.colorbar(label=r'$\Delta t$')
plt.savefig('Plot folder/x_1vsmu1.png', dpi=2000)

plt.show()

# Plotting c residuals
plt.title(r'$c$ input vs output residuals with time delay $\Delta t$')
plt.xlabel(r'$\Delta t$')
plt.ylabel(r'$c$ residual')
#plt.errorbar(delta_t_list, x1_res, yerr=param_err_array[2], fmt='.')
plt.scatter(delta_t_list, c_res, c = mu_1_list, cmap = 'Reds')
plt.colorbar(label=r'$\mu_1$')
plt.savefig('Plot folder/cvsdt.png', dpi=2000)

plt.show()

plt.title(r'c$ input vs output residuals with magnification of image 1')
plt.xlabel(r'$\mu_1$')
plt.ylabel(r'$c$ residual')
#plt.errorbar(mu_1_list, x1_res, yerr=param_err_array[2], fmt='.')
plt.scatter(mu_1_list, c_res, c = delta_t_list, cmap = 'Greens') #which works when we have defined dtlist st it is in order from smallest to largest value - actually works by default
plt.colorbar(label=r'$\Delta t$')
plt.savefig('Plot folder/cvsmu1.png', dpi=2000)

plt.show()


#fig.savefig('Plot folder/lensingresiduals.png', dpi = 2000)   # change this

#%%

# Initialising various arrays to make into corner plots
corner_array1 = np.zeros([5, len(x1_res)])
corner_array2 = np.zeros([3, len(x1_res)])

reduced_chi2 = info_array[0] / info_array[1]

# In its untransposed form, corner array has a row for each variable, and a column for each iteration

corner_array1[0] = param_out_array[3]   #x1 output
corner_array1[1] = param_out_array[4]   #c output
corner_array1[2] = reduced_chi2
corner_array1[3] = delta_t_list
corner_array1[4] = mu_1_list

corner_array2[0] = reduced_chi2
corner_array2[1] = mu_1_list
corner_array2[2] = delta_t_list

# Pearson stat
corr_coef1 = np.corrcoef(corner_array1)

print('Correlation coefficients: ' + str(corr_coef1))

figure = corner.corner(corner_array1.T, labels=[r'$x_1$ output', 'c', r'Reduced $\chi^2$', r'$\Delta t$', r'$\mu_1$'])
figure.savefig('corner.png')

# Note: "I don't believe that you want more dimensions than samples!" error appears when you have more columns (dimensions) than rows (samples), that is, you need to loop over many iterations otherwise corner will just refuse to work...

figure = corner.corner(corner_array2.T, labels=[r'Reduced $\chi^2$', r'$\mu_1$', r'$\Delta t$'])
figure.savefig('corner.png')



# In[51]:

# New questions:
"""
- Should I be making residuals absolute magnitudes or not?
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

- - why does fixing t0 fix everything? and is that even fair, when the point of the project is to calculate
the time delay (i.e t0 for the different lightcurves)?
"""
