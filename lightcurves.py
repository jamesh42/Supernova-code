#!/usr/bin/env python
# coding: utf-8

# In[44]:

import sncosmo
import numpy as np
import matplotlib.pyplot as plt


# Inputting intrinsic parameters

z_val = 0.03
t0_val = 0

# Varying parameters
x1_low = -3
x1_high = 3
c_low = -0.3
c_high = 0.3

x1_list = np.linspace(x1_low, x1_high, 10)
c_list = np.linspace(c_low, c_high, 10)   #param space

# Peak_abs_mag calibration params
mag_filter = 'bessellb'
abs_mag = -19.1
mag_sys = 'ab'


# Creating time list
start_time = -15  #days, what should this be?
end_time = 40
sample_freq = 5 #days, to ensure data is realistic and observation taken every ~5 days

no_time_points = (end_time - start_time) / sample_freq
no_time_points = int(no_time_points)
time_list = np.linspace(start_time, end_time, no_time_points)


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

# Initialising arrays to store input and output parameters for plotting residuals later
z_dat_list = []
t0_dat_list = []
x0_dat_list = []
x1_dat_list = []
c_dat_list = []

z_fit_list = []
t0_fit_list = []
x0_fit_list = []
x1_fit_list = []
c_fit_list = []

x1_err_list = []
c_err_list = []

# In[48]:


# Iterating over parameter space
for val_x1 in x1_list:
    for val_c in c_list:
        try:
            # Step 1: Putting in input parameters and producing data acc. to SALT2 model
            model1 = sncosmo.Model(source='salt2')   #which is now an instance of the class Model
            model1.set(z=z_val, t0=t0_val, x1 = val_x1, c = val_c) 
            model1.set_source_peakabsmag(abs_mag, mag_filter, mag_sys) # readthedocs says do this after setting other parameters

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
            
            # input parameters used to produce data
            z_dat = model1.parameters[0]
            t0_dat = model1.parameters[1]
            x0_dat = model1.parameters[2]
            x1_dat = model1.parameters[3]
            c_dat = model1.parameters[4]
            
            z_dat_list.append(z_dat)
            t0_dat_list.append(t0_dat)
            x0_dat_list.append(x0_dat)
            x1_dat_list.append(x1_dat)
            c_dat_list.append(c_dat)
            
            # Plotting!
            
#             fig1 = sncosmo.plot_lc(data, model=model1)
#             fig1_name = 'input_fig_x1' + str(round(val_x1,2)) + '_c' + str(round(val_c,2)) + '.png'  #data produced by input params
#             fig1.savefig('Lightcurve folder/' + fig1_name)

            # Step 2: Taking produced data and using SALT2 model to produce line of best fit, to ensure that line of best fit 
            # corresponds to original input parameters (we "recover" the input parameters)

            # Initialising new model (to ensure independent test)
            model2 = sncosmo.Model(source='salt2')

            # Run the fit (does this fit over every filter or just one? or does it not matter?)
            result, fitted_model = sncosmo.fit_lc(
                data, model2,
                ['t0', 'x0', 'x1', 'c'], modelcov = True)  ## parameters of model to vary + bounds on parameters (if any) - how do i not do any
            
            # Errors
            x1_err = result.errors['x1']
            c_err = result.errors['c']
            x1_err_list.append(x1_err)
            c_err_list.append(c_err)
            
            # Collecting input values for later residuals
            z_fit = fitted_model.parameters[0]
            t0_fit = fitted_model.parameters[1]
            x0_fit = fitted_model.parameters[2]
            x1_fit = fitted_model.parameters[3]
            c_fit = fitted_model.parameters[4]

            z_fit_list.append(z_fit)
            t0_fit_list.append(t0_fit)
            x0_fit_list.append(x0_fit)
            x1_fit_list.append(x1_fit)
            c_fit_list.append(c_fit)
            
            # Printing fit data

            # print("Number of chi^2 function calls:", result.ncall)
            # print("Number of degrees of freedom in fit:", result.ndof)
            # print("chi^2 value at minimum:", result.chisq)
            # print("model parameters:", result.param_names)
            # print("best-fit values:", result.parameters)
            # print("The result contains the following attributes:\n", result.keys())

            # Plotting lightcurve!
#             fig2_name = 'output_fig_x1' + str(round(val_x1,2)) + '_c' + str(round(val_c,2)) + '.png'   #params produced by data
            #fig2 = sncosmo.plot_lc(data, model=fitted_model, errors=result.errors)   #wrapper for matplotlib.pyplot as plt
            #fig2.savefig('Lightcurve folder/' + fig2_name)
            
        except Exception as e:
            print("Error for x1 = " + str(val_x1) + ' and c = ' + str(val_c))
            print(e)
        
        else:
            print("Success for x1 = " + str(val_x1) + ' and c = ' + str(val_c))
            success_count +=1

        finally:
            continue

print(success_count)

# Code for comparing input to output parameters

# Converting lists to arrays
z_dat_list = np.array(z_dat_list)
t0_dat_list = np.array(t0_dat_list)
x0_dat_list = np.array(x0_dat_list)
x1_dat_list = np.array(x1_dat_list)
c_dat_list = np.array(c_dat_list)

z_fit_list = np.array(z_fit_list)
t0_fit_list = np.array(t0_fit_list)
x0_fit_list = np.array(x0_fit_list)
x1_fit_list = np.array(x1_fit_list)
c_fit_list = np.array(c_fit_list)

# Residuals!
z_res = z_fit_list - z_dat_list #z-residuals
t0_res = t0_fit_list - t0_dat_list 
x0_res = x0_fit_list - x0_dat_list 
x1_res = x1_fit_list - x1_dat_list 
c_res = c_fit_list - c_dat_list

print(x1_err_list)
print(c_err_list)


# In[50]:

# Plotting residuals
fig, (ax1, ax2) = plt.subplots(1,2, figsize = (10,3.5))

# Plotting x1 residuals
ax1.set_title(r'$x_1$ Input vs Output residuals')
ax1.set_xlabel(r'$x_1$')
ax1.set_ylabel(r'$x_1$ residual ($x_1$(output) - $x_1$(input))')
ax1.errorbar(x1_dat_list, x1_res, yerr=x1_err_list, fmt='o')

# Plotting c residuals
ax2.set_title(r'$c$ Input vs Output residuals')
ax2.set_xlabel(r'$c$')
ax2.set_ylabel(r'$c$ residual ($c$(output) - $c$(input))')
#ax2.scatter(c_dat_list, c_res)
ax2.errorbar(c_dat_list, c_res, yerr=c_err_list, fmt='o')

plt.tight_layout(pad=3)   # stops plots from overlapping somehow

fig.savefig('Plot folder/cx1residuals.png', dpi = 2000)

# In[51]:


# New questions:
"""
- what magnitude system should we use (ab, vega, etc.) and why?

- what is zp? why would you use it? just so you can use your own arbitary non-built in unit system?

- is there any reason why you would ever need both bandmag AND bandflux? Or just one or the other

- why does it plot everything twice in some cases?

- bandpass 'sdssg' [3620, .., 5620] outside spectral range [3700, .., 17020] - figure out what this error code means

- does fitting fit using every filter or just one?

- what to do about sdssz breaking code....?
"""


# In[52]:


# misc comments
"""
- why doesnt reducing sample from 5 days to 1 day fix divergence in input output parameters? or, more generally,
why do plots disagree at all?


-should i have a silly amount of plots like I do? Do I really need to record all of them as I have done?

"""

