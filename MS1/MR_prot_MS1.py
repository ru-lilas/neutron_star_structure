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

RED = '#FF4B00'
BLUE = '#005AFF'
GREEN = '#03AF7A'
CYAN = '#4DC4FF'

eos_name = ['APR4', 'ALF2', 'H4', 'MS1']
# filename = ['./'+eos_name[i]+'/NS_MR_relation_'+ eos_name[i] +'.csv' for i in range(4)]
filename = './NS_profile_'+ eos_name[3] +'.csv'
# data_APR4 = pd.read_csv(filename, skip_blank_lines=False)
# data_ALF2 = pd.read_csv(filename, skip_blank_lines=False)
# data_APR4 = pd.read_csv(filename[0], skip_blank_lines=False)
# data_ALF2 = pd.read_csv(filename[1], skip_blank_lines=False)
# data_H4 = pd.read_csv(filename[2], skip_blank_lines=False)
# data_H4 = pd.read_csv(filename, skip_blank_lines=False)
data_MS1 = pd.read_csv(filename, skip_blank_lines=False)
# data_MS1 = pd.read_csv(filename[3], skip_blank_lines=False)

# data convertion (You can rewrite as needed)
# data_APR4['m']   = data_APR4['mass']/M
# data_APR4['rho'] = data_APR4['central_density']
# data_APR4['r']   = data_APR4['radius']/10**5

# data_ALF2['m']   = data_ALF2['mass']/M
# data_ALF2['rho'] = data_ALF2['central_density']
# data_ALF2['r']   = data_ALF2['radius']/10**5

# data_H4['m']   = data_H4['mass']/M
# data_H4['rho'] = data_H4['central_density']
# data_H4['r']   = data_H4['radius']/10**5

data_MS1['m']   = data_MS1['mass']/M
data_MS1['rho'] = data_MS1['central_density']
data_MS1['r']   = data_MS1['radius']/10**5

fig, ax = plt.subplots(1,2,figsize=(24,8))
ax[0].set_xscale('linear')
ax[0].set_yscale('linear')
ax[0].minorticks_on()

# ax[0].plot(data_APR4['rho'] , data_APR4['m'] , color=RED, label=eos_name[0], marker='o', linestyle='')
# ax[0].plot(data_ALF2['rho'] , data_ALF2['m'] , color=GREEN, label=eos_name[1], marker='o', linestyle='')
# ax[0].plot(data_H4['rho']   , data_H4['m']   , color=BLUE, label=eos_name[2], marker='o', linestyle='')
ax[0].plot(data_MS1['rho']  , data_MS1['m']  , color=CYAN, label=eos_name[3], marker='o', linestyle='')

ax[0].tick_params(top=True, right=True, direction='in', which='major')
ax[0].tick_params(top=True, right=True, direction='in', which='minor')
ax[0].set_xlim(0,2.5e15)
ax[0].set_ylim(0,3)
ax[0].set_xlabel(r'$\rho_c$ (g cm$^{-3}$)', fontsize= fontsize_default)
ax[0].set_ylabel(r'$M$ ($M_\odot$)', fontsize= fontsize_default)
ax[0].legend()

ax[1].set_xscale('linear')
ax[1].set_yscale('linear')
ax[1].minorticks_on()

# ax[1].plot(data_APR4['r']   , data_APR4['m'], color=RED, label=eos_name[0], marker='o', linestyle='')
# ax[1].plot(data_ALF2['r']   , data_ALF2['m'], color=GREEN, label=eos_name[1], marker='o', linestyle='')
# ax[1].plot(data_H4['r']     , data_H4['m']  , color=BLUE, label=eos_name[2], marker='o', linestyle='')
ax[1].plot(data_MS1['r']    , data_MS1['m'] , color=CYAN, label=eos_name[3], marker='o', linestyle='')

ax[1].tick_params(top=True, right=True, direction='in', which='major')
ax[1].tick_params(top=True, right=True, direction='in', which='minor')
ax[1].set_xlim(7,20)
ax[1].set_ylim(0,3)
ax[1].set_xlabel(r'$R$ (km)', fontsize= fontsize_default)
ax[1].set_ylabel(r'$M$ ($M_\odot$)', fontsize= fontsize_default)
ax[1].legend()
# plt.gca().xaxis.get_major_formatter().set_useOffset(False)
plt.savefig('./MR_relation'+eos_name[3]+'.png')