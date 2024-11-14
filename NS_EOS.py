#===Equation of State=============================================================================#
import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
import pandas as pd
import matplotlib.pyplot as plt
from astropy import units as u
from astropy import constants as const
import numpy as np

#===Physical Constants============================================================================#
G, M, R, c = const.G.cgs, const.M_sun.cgs, const.R_sun.cgs, const.c.cgs
cc = u.g/(u.cm)**3
pressure = u.erg/(u.cm)**3
#===matplotlib settings===========================================================================#
fontsize_default = 14
fontsize_tick = 14

filename = "./EOS_profile.csv"
data = pd.read_csv(filename)
print(data.head())
index = data['EOS'].tolist()
gamma = [1.35692395,0,0,0]
kappa = [3.99873692e-8*c.value**2,0,0,0]
boundary_rho = [5.0e13, 0, 10**(14.7), 1.0e15, 3.0e15]
boundary_pressure = [0,0,0,0]
#===mk free parameter list========================================================================#
EOS_APR4, EOS_ALF2, EOS_H4, EOS_MS1 = 0,1,2,3

#print(APR4_list)

def P2(logP2):
    return (10**logP2)

def kappa1(P2, gamma1):
    return P2 * boundary_rho[2]**(-gamma1)

def kappa2(kappa1, gamma1, gamma2):
    rho2 = boundary_rho[2]
    return kappa1*rho2**(gamma1-gamma2)

def kappa3(kappa2, gamma2, gamma3):
    rho3 = boundary_rho[3]
    return kappa2*rho3**(gamma2-gamma3)

def rho1(kappa1, gamma1):
    kappa0 = kappa[0]
    gamma0 = gamma[0]
    return (kappa0/kappa1)**(1/(gamma1-gamma0))
#===APR4==========================================================================================#
gamma[1:] = [data.iloc[EOS_APR4][f'gamma{i}'] for i in range(1,4)]
boundary_pressure[2] = P2(data.iloc[EOS_APR4]['logP2'])
kappa[1] = kappa1(boundary_pressure[2], gamma[1])
kappa[2] = kappa2(kappa[1], gamma[1], gamma[2])
kappa[3] = kappa3(kappa[2], gamma[2], gamma[3])
boundary_rho[1] = rho1(kappa[1], gamma[1])

#---add parameter into dataframe------------------------------------#
# EOSdata = [kappa, gamma, boundary_rho[1:]]
# data = pd.DataFrame(data=EOSdata, index=['kappa', 'gamma', 'density_at_boundary'], columns=[0,1,2,3])
# data.to_csv('APR4.csv')
def EOS_pressure(kappa, rho, gamma):
    return kappa*rho**gamma

#===plotting====================================================#
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

rhorange = [np.geomspace(boundary_rho[i], boundary_rho[i+1], 100) for i in range(4)]
P = [EOS_pressure(kappa[i], rhorange[i], gamma[i]) for i in range(4)]

#===ALF2==========================================================================================#
gamma[1:] = [data.iloc[EOS_ALF2][f'gamma{i}'] for i in range(1,4)]
boundary_pressure[2] = P2(data.iloc[EOS_ALF2]['logP2'])
kappa[1] = kappa1(boundary_pressure[2], gamma[1])
kappa[2] = kappa2(kappa[1], gamma[1], gamma[2])
kappa[3] = kappa3(kappa[2], gamma[2], gamma[3])
boundary_rho[1] = rho1(kappa[1], gamma[1])
rhorange_ALF2 = [np.geomspace(boundary_rho[i], boundary_rho[i+1], 100) for i in range(4)]
P_ALF2 = [EOS_pressure(kappa[i], rhorange_ALF2[i], gamma[i]) for i in range(4)]
#=================================================================================================#

#===H4============================================================================================#
gamma[1:] = [data.iloc[EOS_H4][f'gamma{i}'] for i in range(1,4)]
boundary_pressure[2] = P2(data.iloc[EOS_H4]['logP2'])
kappa[1] = kappa1(boundary_pressure[2], gamma[1])
kappa[2] = kappa2(kappa[1], gamma[1], gamma[2])
kappa[3] = kappa3(kappa[2], gamma[2], gamma[3])
boundary_rho[1] = rho1(kappa[1], gamma[1])
rhorange_H4 = [np.geomspace(boundary_rho[i], boundary_rho[i+1], 100) for i in range(4)]
P_H4 = [EOS_pressure(kappa[i], rhorange_H4[i], gamma[i]) for i in range(4)]
#=================================================================================================#

#===MS1===========================================================================================#
gamma[1:] = [data.iloc[EOS_MS1][f'gamma{i}'] for i in range(1,4)]
boundary_pressure[2] = P2(data.iloc[EOS_MS1]['logP2'])
kappa[1] = kappa1(boundary_pressure[2], gamma[1])
kappa[2] = kappa2(kappa[1], gamma[1], gamma[2])
kappa[3] = kappa3(kappa[2], gamma[2], gamma[3])
boundary_rho[1] = rho1(kappa[1], gamma[1])
rhorange_MS1 = [np.geomspace(boundary_rho[i], boundary_rho[i+1], 100) for i in range(4)]
P_MS1 = [EOS_pressure(kappa[i], rhorange_MS1[i], gamma[i]) for i in range(4)]
#=================================================================================================#
#print(rhorange)
#print(P)

[ax.plot(rhorange[i], P[i], color='#FF4B00') for i in range(4)]
[ax.plot(rhorange_ALF2[i], P_ALF2[i], color='#03AF7A', linestyle='--') for i in range(4)]
[ax.plot(rhorange_H4[i], P_H4[i], color='#005AFF', linestyle=':') for i in range(4)]
[ax.plot(rhorange_MS1[i], P_MS1[i], color='#4DC4FF', linestyle='-.') for i in range(4)]
#ax.legend()
plt.savefig('eos.png')