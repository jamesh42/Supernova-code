#!/usr/bin/env python
# coding: utf-8

# ## Notebook to use an SN model to get the lightcurves for a lensed SN
# 
# ## Add the lensing parameters 
# 
# ## Fitting LC's to get SALT2 parameters

# In[37]:


#from ztfquery import fritz
#from ztftarget import target
import sncosmo
import numpy as np
import matplotlib.pyplot as plt 

from astropy.table import Table


# In[39]:


class GLSNe(sncosmo.Source):
    '''
    sntype = "salt2"
    nimages = 2
    '''

    def __init__(self, sntype="salt2", nimages=2, name=None, version=None):
        self.name = name
        self.version = version
        self._sntype = sntype
        self._source = sncosmo.get_source(sntype)
        self._nimages = nimages

        # lensed parameters
        self._parameters = list(np.concatenate([[0, 1] for k in range(1, nimages + 1)]))
        self._param_names = list(np.concatenate([['dt_%i' % k, 'mu_%i' % k] for k in range(1, nimages + 1)]))
        self.param_names_latex = list(np.concatenate([['dt_%i' % k, 'mu_%i' % k] for k in range(1, nimages + 1)]))

        # Add SN parameters
        self._parameters.extend(self._source._parameters)
        self._param_names.extend(self._source._param_names)
        self.param_names_latex.extend(self._source.param_names_latex)

        self._current_parameters = self._parameters.copy()

    def minwave(self):
        return self._source.minwave()

    def maxwave(self):
        return self._source.maxwave()

    def minphase(self):
        return self._source.minphase()

    def maxphase(self):
        return self._source.maxphase()

    def update_param(self):
        param_tmp = list(self._parameters[self._nimages * 2:])
        for n_ in self._source._param_names:
            self._source.set(**{n_: param_tmp.pop(0)})

        self._current_parameters = self._parameters.copy()

    def _flux(self, phase, wave):
        if np.any(self._current_parameters != self._parameters):
            self.update_param()

        out = np.zeros((len(phase), len(wave)))

        # mask = ((phase >= self._sources[0].minphase()) &
        #        (phase <= self._sources[0].maxphase()))
        # out[mask, :] = self._source._flux(phase[mask], wave)

        for k in range(0, self._nimages * 2, 2):
            if k==0:
                dt = self._parameters[k]
                mu = self._parameters[k + 1]
                out[:] = self._source._flux(phase - dt, wave) * mu
            else:
                dt = self._parameters[k]
                mu = self._parameters[k + 1]
                out[:] += self._source._flux(phase - dt, wave) * mu

        return out


# In[45]:


def unresolved_and_individual_lcs(unresolvedmodel, savelc = False, tgrid='continuous', snmodel='salt2-extended', nimages=2, bands=['ztfg', 'ztfr', 'ztfi'],
                                  plotunresolved=True, plotimages=True):
    '''Input model and output the plot with the unresolve lc in solid and individual in dashed lines.
    All bands in the same panel'''

    import matplotlib.pyplot as plt
    import numpy as np
    import sncosmo

    fig, ax = plt.subplots()
    
    if tgrid == 'discrete':
        time = np.linspace(unresolvedmodel.mintime(), unresolvedmodel.maxtime(), 10)
    else:
        time = np.linspace(unresolvedmodel.mintime(), unresolvedmodel.maxtime(), 100)
    # unresolved
    b_arr = np.concatenate([['ztfg' for ii in range(len(time))], ['ztfr' for ii in range(len(time))], 
                   ['ztfi' for ii in range(len(time))]])
    time_arr = np.concatenate([time, time, time])
    if plotunresolved:
        for b in bands:
            if b=='ztfg': color = 'green'
            elif b=='ztfr': color='red'
            elif b=='ztfi': color='orange'
            else: color='grey'
            y = unresolvedmodel.bandmag(b, 'ab', time)
            mask = (~np.isnan(y)) & (~np.isinf(y))
            ax.plot(time[mask], y[mask], '-', color=color, alpha=0.5, label=b)
        
        if savelc: 
            lc_tab = Table()
            y = unresolvedmodel.bandmag(b_arr, 'ab', time_arr)
            mask = (~np.isnan(y)) & (~np.isinf(y))
            lc_tab['time'] = time_arr[mask]
            lc_tab['band'] = b_arr[mask]
            fl =  pow(10, -0.4 * (y[mask] - 25))
            lc_tab['flux'] = fl
            lc_tab['flux_err'] = 0.1 * fl
            lc_tab['zp']  = [25 for ii in range(len(time_arr[mask]))]
            lc_tab['zpsys'] = ['ab' for ii in range(len(time_arr[mask]))]
            
    # individual images lcs
    if plotimages:
        if "salt2" in snmodel:
            dict_param = {'z': unresolvedmodel.parameters[np.where(np.array(unresolvedmodel.param_names) == 'z')],
                          't0': unresolvedmodel.parameters[np.where(np.array(unresolvedmodel.param_names) == 't0')],
                          'x0': unresolvedmodel.parameters[np.where(np.array(unresolvedmodel.param_names) == 'x0')],
                          'x1': unresolvedmodel.parameters[np.where(np.array(unresolvedmodel.param_names) == 'x1')],
                          'c': unresolvedmodel.parameters[np.where(np.array(unresolvedmodel.param_names) == 'c')]
                          }
        else:
            print('Error no dict_param')

        if 'MWebv' in unresolvedmodel.param_names:
            # Add MW extinction
            print('Add MW extinction MWebv=', float(unresolvedmodel.parameters[np.where(np.array(unresolvedmodel.param_names) == 'MWebv')]))
            model = sncosmo.Model(snmodel)
            dust = sncosmo.CCM89Dust()
            model = sncosmo.Model(snmodel, effects=[dust], effect_names=['MW'], effect_frames=['obs'])
            model.set(**{'MWebv':float(unresolvedmodel.parameters[np.where(np.array(unresolvedmodel.param_names) == 'MWebv')])})
        else:
            model = sncosmo.Model(snmodel)

        model.set(**dict_param)

        for n in range(nimages):
            for b in bands:
                if b == 'ztfg':
                    color = 'green'
                elif b == 'ztfr':
                    color = 'red'
                elif b == 'ztfi':
                    color = 'orange'
                else:
                    color = 'grey'
                dt = unresolvedmodel.parameters[np.where(np.array(unresolvedmodel.param_names) == 'dt_'+str(n+1))]
                mu = unresolvedmodel.parameters[np.where(np.array(unresolvedmodel.param_names) == 'mu_'+str(n+1))]
                model.set(**{'x0': dict_param['x0']*mu})
                y = model.bandmag(b, 'ab', time)
                mask = (~np.isnan(y)) & (~np.isinf(y))
                if tgrid == 'discrete':
                    ax.plot(time[mask]+dt, y[mask], color=color, marker='s', alpha=0.25)    
                else:
                    ax.plot(time[mask]+dt, y[mask], '--', color=color, alpha=0.25)


            
    ax.legend(loc=0)
    ax.invert_yaxis()
    ax.set_xlabel('time')
    ax.set_ylabel('mag')
    if savelc:
        return lc_tab
    else:
        return fig, ax



# In[49]:


nimages = 2  # Multiplicity
source = GLSNe("salt2-extended", nimages) # choose sn source
dust = sncosmo.CCM89Dust()
model = sncosmo.Model(source, effects=[dust], effect_names=['MW'] , effect_frames=['obs']) # sncosmo model
MWebv = 0.03
model.set(**{'MWebv':MWebv})
zsource = 0.3544
model.set(**{'z':zsource})
#model.set_source_peakabsmag(-19.1, "bessellb", 'ab')

model.set(**{'x1': 0.})
model.set(**{'c': 0.})
model.set(**{'x0': 6e-05})
model.set(**{'mu_1':9})
model.set(**{'mu_2':5})
model.set(**{'dt_1':-10})
model.set(**{'dt_2':10})
model.set(**{'t0':5})
#tspace = 'discrete'
tspace = 'continuous'

lc = unresolved_and_individual_lcs(model, savelc=True, nimages = nimages, tgrid=tspace)
print(lc)
sncosmo.plot_lc(lc)


# In[ ]:




