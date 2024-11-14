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

filename = './APR4/MRrelation_APR4.csv'
data_APR4 = pd.read_csv(filename)
data_APR4['mass'] =data_APR4['mass'] / M.value
data_APR4['r'] = data_APR4['r'] / 10**5

filename = './ALF2/MRrelation_ALF2.csv'
data_ALF2 = pd.read_csv(filename)
data_ALF2['mass'] =data_ALF2['mass'] / M.value
data_ALF2['r'] = data_ALF2['r'] / 10**5

filename = './H4/MRrelation_H4.csv'
data_H4 = pd.read_csv(filename)
data_H4['mass'] =data_H4['mass'] / M.value
data_H4['r'] = data_H4['r'] / 10**5

filename = './MS1/MRrelation_MS1.csv'
data_MS1 = pd.read_csv(filename)
data_MS1['mass'] =data_MS1['mass'] / M.value
data_MS1['r'] = data_MS1['r'] / 10**5

fig, (ax1,ax2) = plt.subplots(1,2, figsize=(12,6))

# ax1.grid(True)
ax1.set_xscale('linear')
ax1.set_yscale('linear')
ax1.minorticks_on()
ax1.xaxis.set_major_locator(ticker.MultipleLocator(1.0e15))
ax1.xaxis.set_minor_locator(ticker.MultipleLocator(1.0e14))
ax1.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
ax1.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
ax1.plot(data_APR4['rho'], data_APR4['mass'], color='#FF4B00', label='APR4', marker='o')
# ax1.plot(data_ALF2['rho'], data_ALF2['mass'], color='#03AF7A', linestyle='--', label='ALF2')
# ax1.plot(data_H4['rho'], data_H4['mass'], color='#005AFF', linestyle=':', label='H4')
# ax1.plot(data_MS1['rho'], data_MS1['mass'], color='#4DC4FF', linestyle='-.', label='MS1')
ax1.tick_params(top=True, right=True, direction='in', which='major')
ax1.tick_params(top=True, right=True, direction='in', which='minor')
ax1.set_xlim(0,2.5e15)
ax1.set_ylim(0,3.0)
ax1.set_xticks([0,1.0e15,2.0e15])
ax1.set_xlabel(r'$\rho$ (g cm$^{-3}$)', fontsize= fontsize_default)
ax1.set_ylabel(r'$M\,/\,M_{\odot}$', fontsize= fontsize_default)
ax1.legend(loc='lower right')

ax2.set_xscale('linear')
ax2.set_yscale('linear')
ax1.minorticks_on()
ax2.xaxis.set_major_locator(ticker.MultipleLocator(5))
ax2.xaxis.set_minor_locator(ticker.MultipleLocator(1))
ax2.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
ax2.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
ax2.plot(data_APR4['r'], data_APR4['mass'], color='#FF4B00', label='APR4')
# ax2.plot(data_ALF2['r'], data_ALF2['mass'], color='#03AF7A', linestyle='--', label='ALF2')
# ax2.plot(data_H4['r'], data_H4['mass'], color='#005AFF', linestyle=':', label='H4')
# ax2.plot(data_MS1['r'], data_MS1['mass'], color='#4DC4FF', linestyle='-.', label='MS1')
ax2.tick_params(top=True, right=True, direction='in', which='major')
ax2.tick_params(top=True, right=True, direction='in', which='minor')
ax2.set_xlim(7,20)
ax2.set_ylim(0,3.0)
ax2.set_xlabel(r'$R$ (km)', fontsize= fontsize_default)
ax2.set_ylabel(r'$M\,/\,M_{\odot}$', fontsize= fontsize_default)
ax2.legend(loc='upper right')

plt.savefig('MRrelation.png')