#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 10:26:48 2023

@author: jameshenderson
"""

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
from test_mockLC_lensed_visual_resolved import GLSNe
from test_mockLC_lensed_visual_resolved import unresolved_and_individual_lcs
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

#delta_t_list = np.linspace(min_delta_t, max_delta_t, 4)
delta_t_list = [10,10,10,10,10,10,10,10,10]
#delta_t_list = delta_t_list = np.linspace(10,30, 15 )

delta_t_list = [10, 20, 30]

param_out_array = np.zeros((9, 0)) 
param_in_array = np.zeros((9, 0))

params_to_fit = ['x0', 't0', 'dt_1', 'dt_2']
# For reference, param_out_names order is: z,t0,x0,x1,c,dt1,mu1,dt2,mu2
params_to_hold_constant = {'z':z_val} 
params_to_hold_constant = {'z':z_val, 'mu_1': mu_1, 'mu_2': mu_2, 'x0':x0_val, 'x1':x1_val, 'c':c_val} 

# td_list_list = [[-5,5], [-20,0],[-15,15]] # time delay list
# mu_list_list = [[2,7], [7,2], [8,8]]

band_list_ztf =['ztfg', 'ztfr', 'ztfi']
band_list_csp = ['cspjd', 'cspyd', 'csphd']
band_list_full = ['ztfg', 'ztfr', 'ztfi', 'cspjd', 'cspyd', 'csphd']

band_list_list = [band_list_ztf, band_list_csp, band_list_full]

cadence_list = [0.1,1,5,10]

#%%



# Fitting and plotting data using SNTD
for cadence in cadence_list:
#for band_list in band_list_list:
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
        lc = unresolved_and_individual_lcs(model, savelc=True, nimages = nimages, tgrid=tspace, bands = band_list_full, sample_freq = cadence)   #keep this on full, ie simulate all and then only use some
        ascii.write(lc, data_file_name, overwrite=True)    #writing this lensed lc info to file to be read and fitted to by model 2
        data = sncosmo.read_lc(data_file_name)   # producing lensed lc data
    
        # Setting up model2 (fitting model) using SNTD
        image_a = sncosmo.Model('salt2-extended')
        image_b = sncosmo.Model('salt2-extended')
        unresolved=sntd.unresolvedMISN([image_a,image_b])
        
        # Setting parameters
        
        # Setting random delays to see how much it affects final results
        
        #unresolved.set_delays([np.random.random_integers(-20,0), np.random.random_integers(0,20)])
        #unresolved.set_magnifications([np.random.random_integers(0,15), np.random.random_integers(0,15)])
        
        
        unresolved.set_delays([-5,5])   # So does that (comment below) mean that these are essentially just priors?
        #unresolved.set_delays(td_val_list)
        unresolved.set_magnifications([2,7])   # I've noticed if you comment out either of these lines, it breaks. So I think you need to provide this info for it to work.
        #unresolved.set_magnifications(mu_val_list) 
        unresolved.set(z=z_val)
        unresolved.set(t0=t0_val)
        unresolved.set(x0= x0_val)
        
        # Initialising MISN object using unresolved data
        new_MISN=sntd.table_factory(data,telescopename='Unknown',object_name='Unresolved')   # basically just creates MISN instance, but in a clever way to use your given data
        
        fitCurves=sntd.fit_data(new_MISN,snType='Ia', models=unresolved,bands=band_list_full,   # there is another filter ztfg but not shown here because it breaks the code
                        #params=params_to_fit,constants=params_to_hold_constant,   # changing order of params entered here doesnt affect fitcurves.model.parameters order
                        params=params_to_fit,constants={'z':z_val, 'x1':x1_val, 'c':c_val, 'mu_1': mu_1, 'mu_2': mu_2} ,   # changing order of params entered here doesnt affect fitcurves.model.parameters order
                        #bounds={'t0':(-25,25),'x1':(-2,2),'c':(-1,1), 'mu_1': (1,3), 'dt_2': (10,30), 'dt_1': (-10,-30),'mu_2': (6,8)})
                        #bounds=fit_bounds)   #works for the first loop, breaks thereafter? Completely bizarre error
                        bounds={'t0':(-25,25), 'dt_1':(-30,0), 'dt_2':(0,30)}, method='parallel',npoints=100)
        
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
    
        # Making new in and out vals numpy arrays
        out_param_vals = np.array(out_param_vals)
        in_param_vals = np.array(in_param_vals)
        
    
        param_in_array = np.concatenate((param_in_array, in_param_vals.reshape(-1,1)), axis=1) 
        
        param_out_array = np.concatenate((param_out_array, out_param_vals.reshape(-1,1)), axis=1) 
        
        # Computing input overlay
        input_overlay_names = []  # ie names of input parameters being fit
        input_overlay_vals = [] # ie input values of parameters being fit
        
        input_overlay = []   # list of corresponding input parameters for fitted parameters
        # Loop for finding parameters to overlay
        for i in range(len(out_param_names)):
            for j in range(len(params_to_fit)):
                if out_param_names[i] == params_to_fit[j]:   #ie when you've found a match between location of input/output (theyre the same order) parameter and parameter to fit
                    input_overlay.append(in_param_vals[i])   # keep that parameter to fit's input parameter
            else:
                continue
        
        # Plotting
        fig1 = fitCurves.plot_object(showFit=True,plot_unresolved=True)
        fig2 = fitCurves.plot_fit()    # gives corner plot of probabilities
    
        # Extract the axes
        ndim  = 4  # must be equal to number of graphs, couldnt soft code because of weird error with defining bounds etc. outside of loop    
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
        
        plt.show()
        fig1.savefig('double_lc.png')
        fig2.savefig('cornerdouble.png')
    
        res = fitCurves.images['image_1'].fits.res
            
        print(res)
        samples = res.samples
        samples = res.samples.T #making it so that columns become rows, ie rows
        print('samples:')
        print(samples)
        fig3 = corner.hist2d(
                samples[2], samples[3])
        plt.xlabel('dt_1')
        plt.ylabel('dt_2')
        
        print(input_overlay[2:])
        
        plt.plot(input_overlay[2], input_overlay[3], color='green')
        plt.show()
        #fig3.savefig('cornerhist.png')
        
        chi = fitCurves.images['image_1'].chisq
        reduced_chi = fitCurves.images['image_1'].reduced_chisq    # THIS IS HOW YOU OBTAIN IT. 
        #FOR SOME REASON, fitCurves.series_chisq doesn't work, even though series_chisq IS a property of the MISN class too
        print(chi)
        print(reduced_chi)
        # print('x')
    
    

#%% 

# Residuals

label_list = ['ztf', 'csp', 'both']
label_list = ['0.1 day cadence', '1 day cadence', '5 day cadence', '10 day cadence']
colour_list = ['red','blue','green', 'purple']

time_len = len(delta_t_list)

for i in range(len(out_param_vals)):
    plt.title(str(in_param_names[i] + ' residuals'))
    plt.xlabel(r'$\Delta t$')
    plt.ylabel(str(in_param_names[i] + ' residuals'))
    # for j in range(len(param_in_array[0])):
        # plt.scatter(delta_t_list_2[0], (param_out_array[i][j] - param_in_array[i][j]), c=colour_list[j], label=label_list[j])
    for j in range(len(label_list)):
        plt.scatter(delta_t_list, (param_out_array[i][time_len*j:time_len*(j+1)] - param_in_array[i][time_len*j:time_len*(j+1)]), c=colour_list[j], label=label_list[j])   #only works with sim loops as they are
    #plt.colorbar(label=r'$\mu_1$')
    plt.legend()
    plt.plot()
    #plt.legend()
    plt.show()
    #plt.savefig('Plot folder/x_1vsdt.png', dpi=2000)


# Is this only image 1? I suppose the plot says only image 1 is fitted, and so I'd imagine that this is really for the whole fit
chi = fitCurves.images['image_1'].chisq
reduced_chi = fitCurves.images['image_1'].reduced_chisq    # THIS IS HOW YOU OBTAIN IT. 
#FOR SOME REASON, fitCurves.series_chisq doesn't work, even though series_chisq IS a property of the MISN class too
print(chi)
print(reduced_chi)
# print('x')



"""
Questions:
    
- why does SNTD need to exist for RESOLVED images? Would the process not be really easy? Just subtract the time for one peak from the other?
    
- do the 'priors' matter? Comparison video on phone from set_delays (-5,5) to (-5,20) and set_magnifications from i think (2,5) to (10,1)
 - CLEARLY THEY ARE SOMEWHAT SENSITIVE TO IT, initial guesses
 
 - why is t0 real range not shown

- plot together
fix everything except time delays

r, i, g,z j nad h, y g and h

Essentially just do with time delay parameters, just do with parameters he said originally, see results,
try to figure out how important initial guesses are, put corner plots and lightcurves on the same figure in overleaf to avoid confusion,
then add the other filters to see how effective these are at reducing errors on time delays etc.

Apparently maybe nobody has tested the filters' effects on this analysis before, so that would be very interesting'

"""
"""
Explain microlensing maps

SNTD microlensing treatment -- Stephen thought this was strange

Fitting a gaussian process- what is this?

Read SNTD paper - how does SNTD actually work?

Presumably also read SNCosmo paper too

Less important:
    - what is a truncated gaussian compared to what is gaussian? What is gaussian fitting in general?

    - remind self of what the SNe parameters are, Rv, Av, mu etc.
"""





