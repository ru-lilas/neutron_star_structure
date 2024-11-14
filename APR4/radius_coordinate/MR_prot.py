import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from astropy import units as u
from astropy import constants as const
import numpy as np

G, M, R = const.G.cgs, const.M_sun.cgs, const.R_sun.cgs

fontsize_default = 16
fontsize_tick = 14

eos_name = 'APR4'
filename = './NS_MR_relation_'+ eos_name +'.csv'
data = pd.read_csv(filename, skip_blank_lines=False)


# data convertion (You can rewrite as needed)
data['m']   = data['y[1]']/M
data['rho'] = data['y[0]']
data['r']   = data['t']/10**5


fig, ax = plt.subplots(1,2,figsize=(12,6))
ax[0].set_xscale('linear')
ax[0].set_yscale('linear')
ax[0].minorticks_on()
ax[0].plot(data['rho'], data['m'], color='#FF4B00', label=eos_name, marker='o', linestyle='')
ax[0].tick_params(top=True, right=True, direction='in', which='major')
ax[0].tick_params(top=True, right=True, direction='in', which='minor')
ax[0].set_xlim(0,2.5e15)
ax[0].set_xlabel(r'$\rho$ (cc)', fontsize= fontsize_default)
ax[0].set_ylabel(r'$M$ ($M_\odot$)', fontsize= fontsize_default)

ax[1].set_xscale('linear')
ax[1].set_yscale('linear')
ax[1].minorticks_on()
ax[1].plot(data['r'], data['m'], color='#FF4B00', label=eos_name, marker='o', linestyle='')
ax[1].tick_params(top=True, right=True, direction='in', which='major')
ax[1].tick_params(top=True, right=True, direction='in', which='minor')
ax[1].set_xlim(7,20)
ax[1].set_xlabel(r'$R$ (km)', fontsize= fontsize_default)
ax[1].set_ylabel(r'$M$ ($M_\odot$)', fontsize= fontsize_default)
# plt.gca().xaxis.get_major_formatter().set_useOffset(False)
plt.savefig('MR_relation_'+eos_name+'.png')