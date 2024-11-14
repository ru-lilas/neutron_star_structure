import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from astropy import units as u
from astropy import constants as const
import numpy as np

G, M, R = const.G.cgs, const.M_sun.cgs, const.R_sun.cgs

fontsize_default = 14
fontsize_tick = 14

filename = './runge_kutta_4_output.csv'
eos_name = 'APR4'
data = pd.read_csv(filename, skip_blank_lines=False)


# data convertion (You can rewrite as needed)
data['m']   = data['y1'] / M
data['rho'] = data['y0']
data['r']   = data['t']

# data processing (DON'T REWRITE VARIABLE NAME)
value   = data[['r', 'rho', 'm', 'P', 'H']]
k_0     = data[['t','k1[0]', 'k2[0]', 'k3[0]', 'k4[0]']]
k_1     = data[['t','k1[1]', 'k2[1]', 'k3[1]', 'k4[1]']]

print(value.head())
filepath = './value_NS_' + eos_name + '.csv'
value.to_csv(filepath, float_format='%.10e', index=False)

fig, ax = plt.subplots(3,2, figsize=(16,24))
ax[0,0].set_xscale('log')
ax[0,0].set_yscale('log')
ax[0,0].minorticks_on()
ax[0,0].plot(value['r'], value['rho'], color='#FF4B00', label='APR4', marker='o', linestyle='')
ax[0,0].tick_params(top=True, right=True, direction='in', which='major')
ax[0,0].tick_params(top=True, right=True, direction='in', which='minor')
ax[0,0].set_xlabel(r'$r$ (cm)', fontsize= fontsize_default)
ax[0,0].set_ylabel(r'$\rho$ (cc)', fontsize= fontsize_default)

ax[0,1].set_xscale('linear')
ax[0,1].set_yscale('linear')
ax[0,1].minorticks_on()
ax[0,1].plot(value['r'], value['m'], color='#FF4B00', label='APR4', marker='o', linestyle='')
ax[0,1].tick_params(top=True, right=True, direction='in', which='major')
ax[0,1].tick_params(top=True, right=True, direction='in', which='minor')
ax[0,1].set_xlabel(r'$r$ (cm)', fontsize= fontsize_default)
ax[0,1].set_ylabel(r'$m$ / $M_\odot$', fontsize= fontsize_default)

ax[1,0].set_xscale('log')
ax[1,0].set_yscale('log')
ax[1,0].minorticks_on()
ax[1,0].plot(value['rho'], value['P'], color='#FF4B00', label='APR4', marker='o', linestyle='')
ax[1,0].tick_params(top=True, right=True, direction='in', which='major')
ax[1,0].tick_params(top=True, right=True, direction='in', which='minor')
ax[1,0].set_xlabel(r'$\rho$ (cc)', fontsize= fontsize_default)
ax[1,0].set_ylabel(r'$P$ (dyn cm$^{-2}$)', fontsize= fontsize_default)

ax[1,1].set_xscale('log')
ax[1,1].set_yscale('log')
ax[1,1].minorticks_on()
ax[1,1].plot(value['r'], value['H'], color='#FF4B00', label='APR4', marker='o', linestyle='')
ax[1,1].tick_params(top=True, right=True, direction='in', which='major')
ax[1,1].tick_params(top=True, right=True, direction='in', which='minor')
ax[1,1].set_xlabel(r'$r$ (cm)', fontsize= fontsize_default)
ax[1,1].set_ylabel(r'Scale height $H$ (cm)', fontsize= fontsize_default)

ax[2,0].set_xscale('log')
ax[2,0].set_yscale('log')
ax[2,0].minorticks_on()
ax[2,0].plot(value['r'], value['P'], color='#FF4B00', label='APR4', marker='o', linestyle='')
ax[2,0].tick_params(top=True, right=True, direction='in', which='major')
ax[2,0].tick_params(top=True, right=True, direction='in', which='minor')
ax[2,0].set_xlabel(r'$r$ (cm)', fontsize= fontsize_default)
ax[2,0].set_ylabel(r'$P$ (dyn cm$^{-2}$)', fontsize= fontsize_default)
plt.savefig("process.png")