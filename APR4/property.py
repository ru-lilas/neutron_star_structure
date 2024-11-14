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
# data = pd.read_csv(filename)
data = pd.read_csv(filename, skip_blank_lines=False)

#print(data.applymap(lambda x: isinstance(x, str)))
blanck_indice = data[data.isnull().all(axis=1)].index
# print(blanck_indice)
#data = data[blanck_indice[100]+2:blanck_indice[101]]
#data = data.apply(pd.to_numeric, errors='coerce')

#print(data)


data['mass'] = data['m'] / M.value
data['P'] = 10**data['logP']
# data['rho'] = data['y0']
# data['r'] = data['y1']
label_R = r"$R$ (cm)"

print(data.head())

fig, ax = plt.subplots(2,2, figsize=(16,16))

ax[0,0].set_xscale('log')
ax[0,0].set_yscale('log')
ax[0,0].minorticks_on()
ax[0,0].plot(data['rho'], data['P'], color='#FF4B00', label='APR4', marker='o', linestyle='')
ax[0,0].tick_params(top=True, right=True, direction='in', which='major')
ax[0,0].tick_params(top=True, right=True, direction='in', which='minor')
ax[0,0].set_xlabel(r'$\rho$(g cm-2)', fontsize= fontsize_default)
ax[0,0].set_ylabel(r'$P$ (dyn cm-2)', fontsize= fontsize_default)
# ax[2,1].legend(loc='upper right')

ax[0,1].set_xscale('log')
ax[0,1].set_yscale('log')
ax[0,1].minorticks_on()
ax[0,1].plot(data['r'], data['rho'], color='#FF4B00', label='APR4', marker='o', linestyle='')
ax[0,1].tick_params(top=True, right=True, direction='in', which='major')
ax[0,1].tick_params(top=True, right=True, direction='in', which='minor')
ax[0,1].set_xlabel(r'$R$(cm)', fontsize= fontsize_default)
ax[0,1].set_ylabel(r'$\rho$', fontsize= fontsize_default)
# ax[0,1].legend(loc='upper right')

ax[1,0].set_xscale('log')
ax[1,0].set_yscale('log')
ax[1,0].minorticks_on()
ax[1,0].plot(data['r'], data['P'], color='#FF4B00', label='APR4', marker='o', linestyle='')
ax[1,0].tick_params(top=True, right=True, direction='in', which='major')
ax[1,0].tick_params(top=True, right=True, direction='in', which='minor')
ax[1,0].set_xlabel(label_R, fontsize= fontsize_default)
ax[1,0].set_ylabel(r'$P$ (dyn cm-2)', fontsize= fontsize_default)
# ax[2,0].legend(loc='upper right')

ax[1,1].set_xscale('log')
ax[1,1].set_yscale('linear')
ax[1,1].minorticks_on()
ax[1,1].plot(data['r'], data['mass'], color='#FF4B00', label='APR4', marker='o', linestyle='')
ax[1,1].tick_params(top=True, right=True, direction='in', which='major')
ax[1,1].tick_params(top=True, right=True, direction='in', which='minor')
ax[1,1].set_xlabel(label_R, fontsize= fontsize_default)
ax[1,1].set_ylabel(r'Enclosed Mass /$\,M_{\odot}$', fontsize= fontsize_default)
# ax[2,1].legend(loc='upper right')
plt.savefig("enclosed_mass.png")