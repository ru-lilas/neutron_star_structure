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
data = pd.read_csv(filename)


# data convertion (You can rewrite as needed)
data['m']   = data['t']
data['rho'] = data['y0']
data['r']   = data['y1']

# data processing (DON'T REWRITE VARIABLE NAME)
value   = data[['r', 'rho', 'm']]
k_0     = data[['t','k1[0]', 'k2[0]', 'k3[0]', 'k4[0]']]
k_1     = data[['t','k1[1]', 'k2[1]', 'k3[1]', 'k4[1]']]

# colour code
RED = '#FF4B00'
GRN = '#03AF7A'
BLE = '#005AFF'
CYN = '#4DC4FF'

print(value.head())
filepath = './value_NS_' + eos_name + '.csv'
value.to_csv(filepath, float_format='%.10e', index=False)

fig, ax = plt.subplots(2,2, figsize=(16,16))
ax[0,0].set_xscale('log')
ax[0,0].set_yscale('log')
ax[0,0].minorticks_on()
ax[0,0].plot(value['m'], value['rho'], color = RED, label='APR4', marker='o', linestyle='')
ax[0,0].tick_params(top=True, right=True, direction='in', which='major')
ax[0,0].tick_params(top=True, right=True, direction='in', which='minor')
ax[0,0].set_xlabel(r'$m$ (g)', fontsize= fontsize_default)
ax[0,0].set_ylabel(r'$\rho$ (g/cm$^{3}$)', fontsize= fontsize_default)

ax[0,1].set_xscale('log')
ax[0,1].set_yscale('log')
ax[0,1].minorticks_on()
ax[0,1].plot(value['m'], value['r'], color = GRN, label='APR4', marker='o', linestyle='')
ax[0,1].tick_params(top=True, right=True, direction='in', which='major')
ax[0,1].tick_params(top=True, right=True, direction='in', which='minor')
ax[0,1].set_xlabel(r'$m$ (g)', fontsize= fontsize_default)
ax[0,1].set_ylabel(r'$r$ (cm)', fontsize= fontsize_default)
plt.savefig("physical_value.png")

fig, ax = plt.subplots(1,2,figsize=(16,8))
ax[0].set_xscale('log')
ax[1].set_xscale('log')
ax[0].set_yscale('log')
ax[1].set_yscale('log')
ax[0].minorticks_on()
ax[1].minorticks_on()
ax[0].plot(data['t'], abs(data['k1[0]']), color = RED, marker='', linestyle='-', label = r'$k_1$')
ax[0].plot(data['t'], abs(data['k2[0]']), color = GRN, marker='', linestyle='-', label = r'$k_2$')
ax[0].plot(data['t'], abs(data['k3[0]']), color = BLE, marker='', linestyle='-', label = r'$k_3$')
ax[0].plot(data['t'], abs(data['k4[0]']), color = CYN, marker='', linestyle='-', label = r'$k_4$')
ax[1].plot(data['t'], data['k1[1]'], color = RED, marker='', linestyle='-', label = r'$k_1$')
ax[1].plot(data['t'], data['k2[1]'], color = GRN, marker='', linestyle='-', label = r'$k_2$')
ax[1].plot(data['t'], data['k3[1]'], color = BLE, marker='', linestyle='-', label = r'$k_3$')
ax[1].plot(data['t'], data['k4[1]'], color = CYN, marker='', linestyle='-', label = r'$k_4$')
ax[0].tick_params(top=True, right=True, direction='in', which='major')
ax[0].tick_params(top=True, right=True, direction='in', which='minor')
ax[0].set_xlabel(r'$t$', fontsize= fontsize_default)
ax[1].tick_params(top=True, right=True, direction='in', which='major')
ax[1].tick_params(top=True, right=True, direction='in', which='minor')
ax[1].set_xlabel(r'$t$', fontsize= fontsize_default)
ax[0].set_ylabel(r'$k$[0]', fontsize= fontsize_default)
ax[1].set_ylabel(r'$k$[1]', fontsize= fontsize_default)
ax[0].legend()
ax[1].legend()
plt.savefig('t_k[1].png')