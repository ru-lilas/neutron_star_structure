#===Equation of State MS1=========================================================================#
import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
import sys
import pandas as pd
import matplotlib.pyplot as plt
from astropy import units as u
from astropy import constants as const
import numpy as np
EOS_NAME = 'MS1'
EOS_TYPE=3 # APR4=0, ALF2=1, H4=2, MS1=3

#===Physical Constants============================================================================#
G, M, R, c = const.G.cgs, const.M_sun.cgs, const.R_sun.cgs, const.c.cgs
cc = u.g/(u.cm)**3
pressure = u.erg/(u.cm)**3

#===matplotlib settings===========================================================================#
fontsize_default = 14
fontsize_tick = 14

#===Reading EOS profile===========================================================================#
filename = "../EOS_profile.csv"
data = pd.read_csv(filename)
print(data.head())

#===Initialize array==============================================================================#
index = data['EOS'].tolist()
gamma = [1.35692395,0,0,0]
kappa = [3.99873692e-8*c.value**2,0,0,0]
# kappa = [3.99873692e-8,0,0,0]
boundary_rho = [5.0e13, 0, 10**(14.7), 1.0e15, 3.0e15]

#===Define calculation of kappa and rho=========================================================#
def kappa_value(piecewise_i):
    if type(piecewise_i) != int:
        sys.exit('Error: variable "piecewise" is not int type.')
    elif piecewise_i == 1:
        value = 10**(data.iloc[EOS_TYPE]['logP2']) * boundary_rho[2]**(-gamma[1])
        return value
    elif piecewise_i in [2,3]:
        value = kappa[piecewise_i - 1]*boundary_rho[piecewise_i]**(gamma[piecewise_i-1] - gamma[piecewise_i])
        return value
    else:
        sys.exit(f'Error: Invalid piecewise_i = {piecewise_i}')

def a_constant(piecewise_i):
    if type(piecewise_i) != int:
        sys.exit('Error: "piecewise_i" is not int type.')
    elif piecewise_i == 0:
        value = 0
        return value
    elif piecewise_i in [1,2,3]:
        value = \
            a_constant(piecewise_i-1) \
            + kappa[piecewise_i-1]*boundary_rho[piecewise_i-1]**(gamma[piecewise_i-1]-1) / ((gamma[piecewise_i-1] - 1)*(c.value)**2)\
            - kappa[piecewise_i]*boundary_rho[piecewise_i-1]**(gamma[piecewise_i]-1)  / ((gamma[piecewise_i] - 1)*(c.value)**2)
        return value
    else:
        sys.exit(f'Error: Invalid piecewise_i = {piecewise_i}')

def rho1(kappa1, gamma1):
    kappa0 = kappa[0]
    gamma0 = gamma[0]
    return (kappa0/kappa1)**(1/(gamma1-gamma0))

def EOS_pressure(kappa, rho, gamma):
    return kappa*rho**gamma
#===APR4 calculation=============================================================================#
gamma[1:] = [data.iloc[EOS_TYPE][f'gamma{i}'] for i in range(1,4)]
# kappa[1:] = [kappa_value(i) for i in range(1,4)]
kappa[1] = kappa_value(1)
kappa[2] = kappa_value(2)
kappa[3] = kappa_value(3)
boundary_rho[1] = rho1(kappa[1], gamma[1])
a_parameter = [a_constant(i) for i in range(4)]
#===dataframe=====================================================================================#
EOSdata = [kappa, gamma, boundary_rho[:3], a_parameter]
data = pd.DataFrame(data=EOSdata, index=['kappa', 'gamma', 'density_at_boundary', 'a_parameter'], columns=[0,1,2,3])
data.to_csv(EOS_NAME+'.csv', float_format='%.10e', index='False')
#===plotting======================================================================================#
fig, ax = plt.subplots()
#---plot settings---------------------------------------------------------------------------------#
ax.tick_params(top=True, right=True, direction='in', which='minor')
ax.tick_params(top=True, right=True, direction='in', which='major')
ax.set_xlabel(r'$\rho$ (g/cm3)', fontsize= fontsize_default)
ax.set_ylabel(r'$P$ (dyn/cm2)', fontsize= fontsize_default)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(5.0e13, 3.0e15)
ax.set_ylim(1.0e32, 1.0e37)
#-------------------------------------------------------------------------------------------------#
x = [np.geomspace(boundary_rho[i], boundary_rho[i+1], 100) for i in range(4)] # rho range
y = [EOS_pressure(kappa[i], x[i], gamma[i]) for i in range(4)] # pressure range

# [ax.plot(x[i], y[i], color='#FF4B00') for i in range(4)]
ax.plot(x[0], y[0], color = '#FF4B00')
ax.plot(x[1], y[1], color = '#FF4B00', linestyle='--')
ax.plot(x[2], y[2], color = '#FF4B00', linestyle=':')
ax.plot(x[3], y[3], color = '#FF4B00', linestyle='-.')

plt.savefig('./eos_'+EOS_NAME+'.png')