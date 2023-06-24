#!/usr/bin/env python
# coding: utf-8

# In[38]:


import sncosmo
import numpy as np
import matplotlib.pyplot as plt


# In[39]:


# Inputting intrinsic parameters

z_val = 0.03
t0_val = 0
x1_list = np.linspace(-3, 3, 5)
c_list = np.linspace(-0.3, 0.3, 5)   #param space

# Peak_abs_mag calibration params
mag_filter = 'bessellb'
abs_mag = -19.1
mag_sys = 'ab'


# In[40]:


# Creating time list
start_time = -40   #days
end_time = 40
sample_freq = 5 #days, to ensure data is realistic and observation taken every ~5 days

no_time_points = (end_time - start_time) / sample_freq
no_time_points = int(no_time_points)
time_list = np.linspace(-40, 40, no_time_points)


# In[41]:


# Producing and organising data for plotting

import warnings
warnings.filterwarnings('ignore') # to make output easier to read

from astropy.io import ascii
from astropy.table import Table

# Some params
frac_err = 0.1 #10%?
zp_val = 25. #why?
filter_list = ['sdssg', 'sdssr', 'sdssi', 'sdssz'] # code breaks for sdssz
filter_list = ['sdssg', 'sdssr', 'sdssi']
data_file_name = 'data.txt'

print("Other input parameters are z = " + str(z_val) + " and t0 = " + str(t0_val) + " with absolute magnitude in " \
      + mag_filter + " of " + str(abs_mag))

success_count = 0 

# Iterating over parameter space
for val_x1 in x1_list:
    for val_c in c_list:
        try:
            # Step 1: Putting in input parameters and producing data acc. to SALT2 model
            model1 = sncosmo.Model(source='salt2')   #which is now an instance of the class Model
            model1.set(z=z_val, t0=t0_val, x1 = val_x1, c = val_c) 
            model1.set_source_peakabsmag(abs_mag, mag_filter, mag_sys) # readthedocs says do this after setting other parameters

#             print("Iteration input parameters")
#             print(model1.parameters)

            # Initialising lists that will be used to fill table
            flux_list = []
            filter_list_list = []
            time_list_list = list(time_list) * len(filter_list)   #must make list again (not np array) so that * len(filter_list) works

            # Initialising data table
            data = Table()

            # Filling lists
            for filter in filter_list: #computing flux in each filter, then appending to total list (since plot_lc can read all data at once)
                try:
                    # Creating synthetic photometry
                    mod1_flux = model1.bandflux(filter, time_list, zp=25., zpsys='ab')
                    for flux in mod1_flux:   #for flux (at a given time)
                        flux_list.append(flux)
                        filter_list_list.append(str(filter))   #noting corresponding filter
                        
                except Exception as e:
                    #print(e)
                    pass
                
                else:
                    #print("Filter OK")
                    pass
                
                finally:
                    continue

            flux_err_list = np.array(flux_list) * frac_err   # applying systematic error
            flux_err_list = np.abs(flux_err_list)    #fixes flux error issue

            # Using ascii in Astropy to create data set
            data['time'] = time_list_list #dtype?  # list of times for which we have data
            data['flux'] = flux_list   # corresponding list of simulated fluxes for each time
            data['fluxerr'] = flux_err_list   # applying systematic error
            data['band'] = filter_list_list   # list of filter name for each data point (which in this case is the same, always)
            data['zp'] = [zp_val] * no_time_points * len(filter_list)   # code requires this so whatever
            data['zpsys'] = ['ab'] * no_time_points * len(filter_list)  # ^
            ascii.write(data, data_file_name, overwrite=True)   # overwriting is okay bc could have used dataframe anyway

            # Finally compiling data file for this iteration
            data = sncosmo.read_lc(data_file_name)   #read_lc is sncosmo's method for reading lightcurve data 

            # Plotting!
            print("---")
            print(model1.parameters)
            
            #sncosmo.plot_lc(data, model=model1)

            # Step 2: Taking produced data and using SALT2 model to produce line of best fit, to ensure that line of best fit 
            # corresponds to original input parameters (we "recover" the input parameters)

            # Initialising new model (to ensure independent test)
            model2 = sncosmo.Model(source='salt2')

            # Run the fit (does this fit over every filter or just one? or does it not matter?)
            result, fitted_model = sncosmo.fit_lc(
                data, model2,
                ['z', 't0', 'x0','x1', 'c'],  # parameters of model to vary
                bounds={'z':(0.1, 0.7)}, modelcov = True)  # bounds on parameters (if any) - how do i not do any

            #print("Iteration output parameters")
            #print(fitted_model.parameters)


            # Printing fit data

    #         print("Number of chi^2 function calls:", result.ncall)
    #         print("Number of degrees of freedom in fit:", result.ndof)
    #         print("chi^2 value at minimum:", result.chisq)
    #         print("model parameters:", result.param_names)
    #         print("best-fit values:", result.parameters)
    #         print("The result contains the following attributes:\n", result.keys())

            #print(data)

            # Plotting lightcurve!

            #sncosmo.plot_lc(data, model=fitted_model, errors=result.errors)   #wrapper for matplotlib.pyplot as plt

        except Exception as e:
            print("Error for x1 = " + str(val_x1) + ' and c = ' + str(val_c))
            print(e)
        
        else:
            print("Success for x1 = " + str(val_x1) + ' and c = ' + str(val_c))
            success_count +=1

        finally:
            continue

print(success_count)


# In[42]:


# Code for comparing input to output parameters (it does not work yet, I'll come back to it)

old_x0 = model1

print(model1.param_names)
print(model1.parameters)

# input parameters used to produce data
z_dat = model1.parameters[0]
t0_dat = model1.parameters[0]
x0_dat = model1.parameters[2]
x1_dat = model1.parameters[3]
c_dat = model1.parameters[4]

# parameters produced by using SALT2 to be fit data produced originally
z_fit = model2.parameters[0]
t0_fit = model2.parameters[0]
x0_fit = model2.parameters[2]
x1_fit = model2.parameters[3]
c_fit = model2.parameters[4]

plt.plot(z_dat, z_fit)


# In[43]:


# New questions:
"""
- what magnitude system should we use (ab, vega, etc.) and why?

- what is zp? why would you use it? just so you can use your own arbitary non-built in unit system?

- is there any reason why you would ever need both bandmag AND bandflux? Or just one or the other

- why does it plot everything twice in some cases?

- bandpass 'sdssg' [3620, .., 5620] outside spectral range [3700, .., 17020] - figure out what this error code means

- does fitting fit using every filter or just one?
"""

