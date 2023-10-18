# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 15:15:04 2020

@author: I Love Python
"""


# import
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import os
# Prepare
# Path
main_path = r'C:\Users\n8748\Host\Class Report & Work\108-2_AtmSci2002_大氣熱力學\CA\CA8.5'
ERA5_path = r'C:\Users\n8748\Host\Class Report & Work\108-2_AtmSci2002_大氣熱力學\CA\CA5\ERA5'
# Read .nc file
os.chdir(ERA5_path)
ERA5_Bengkulu = nc.Dataset('Bengkulu_2019_pl_tq.nc')

# Constants
    # Physical
Rd = 287.05 # Air constant [J/kg*K]
g = 9.81 # Gravity accelation
Lv = 2.5E6 # J/kg
Rw = 461.5 # J/(kg·K)
Cp = 1004
    # Positions(index)
lat = 2
lon = 2

# Data Handling

    # Def fucntion
    
def Tv(Temp,MR):
    res = Temp*(1+0.608*MR)
    res = np.flip(res)
    return res
def Alt(Temp_v,Pressure):
    # Pre-
    Pres = np.flip(Pressure)
    res = np.array([])
    # I.C.
    res = np.append(res,0)
    # Append
    for i in range(len(Pres)-1):
        temp = res[-1] + Rd*Temp_v[i]/g*np.log(Pres[i]/Pres[i+1])
        res = np.append(res,temp)
    return res
def Qvs(Pres,Temp):
    e_s0 = 6.11 # hPa
    T_0 = 273 # K
    epsilon = Rd/Rw
    e_s = e_s0 * np.exp((Lv/Rw)*((1/T_0)-(1/Temp)))
    qvs = (epsilon*e_s)/(Pres - (1-epsilon)*e_s)
    return qvs
def T_vp(IC_Temp,IC_MR,Pressure,Alt):
    # Pre-
    res = np.array([])
    Pres = np.flip(Pressure)
    q = np.array([])
    # I.C.
    res = np.append(res,IC_Temp)
    q = np.append(q,IC_MR)
    # Append
    for i in range(len(Pres)-1):
        if IC_MR < Qvs(Pres[i],res[i]):
            temp = res[i] + (-g/Cp)*(Alt[i+1]-Alt[i])
            res = np.append(res,temp)
            q = np.append(q,IC_MR)
        else:
            Cp_star = Cp*(1 + (Lv**2*Qvs(Pres[i],res[i]))/(Cp*Rw*res[-1]**2))
            temp = res[i] + (-g/Cp_star)*(Alt[i+1]-Alt[i])
            q = np.append(q,Qvs(Pres[i],res[i]))
            res = np.append(res,temp)
    res = np.flip(Tv(res,q))
    return res
def CAPE(Tvp,Tve,Alt):
    CAPE = 0
    for i in np.where(Tvp >= Tve)[0][:-1]:
        dz = Alt[i+1] - Alt[i]
        CAPE += ((Tvp[i] - Tve[i])/Tve[i])*g*dz
    return CAPE
def RH_surface(IC_MR,IC_Temp):
    qv = IC_MR
    qvs = Qvs(1000, IC_Temp)
    RH = qv/qvs
    return RH
    # Calculation
    
    # Q2
Time = np.linspace(1,12,1460)
CAPE_yr = np.array([])

for i in range(1460):
    Pres = ERA5_Bengkulu.variables['level'][10:]
    Temp = ERA5_Bengkulu.variables['t'][i,:,2,2][10:]
    MR = ERA5_Bengkulu.variables['q'][i,:,2,2][10:]
    Tve = Tv(Temp,MR)
    Height = Alt(Tve,Pres)
    Tvp = T_vp(Temp[-1],MR[-1],Pres,Height)
    temp = CAPE(Tvp,Tve,Height)
    CAPE_yr = np.append(CAPE_yr,temp)
ave_Q2 = [np.mean(CAPE_yr)]*1460
    
    # Q3
Time = np.linspace(1,12,1460)
CAPE_yr_3K = np.array([])
for i in range(1460):
    Pres = ERA5_Bengkulu.variables['level'][10:]
    Temp = ERA5_Bengkulu.variables['t'][i,:,2,2][10:]
    MR = ERA5_Bengkulu.variables['q'][i,:,2,2][10:]
    Tve = Tv(Temp,MR)
    Height = Alt(Tve,Pres)
    Tvp = T_vp(Temp[-1]+3,MR[-1],Pres,Height)
    temp = CAPE(Tvp,Tve,Height)
    CAPE_yr_3K = np.append(CAPE_yr_3K,temp)
ave_Q3 = [np.mean(CAPE_yr_3K)]*1460
    
    # Q4
Time = np.linspace(1,12,1460)
CAPE_yr_3K_RH = np.array([])
for i in range(1460):
    Pres = ERA5_Bengkulu.variables['level'][10:]
    Temp = ERA5_Bengkulu.variables['t'][i,:,2,2][10:]
    MR = ERA5_Bengkulu.variables['q'][i,:,2,2][10:]
    Tve = Tv(Temp,MR)
    Height = Alt(Tve,Pres)
    RH = RH_surface(MR[-1], Temp[-1])
    Fix_MR = RH * Qvs(Pres[-1], Temp[-1]+3)
    Tvp = T_vp(Temp[-1]+3,Fix_MR,Pres,Height)
    temp = CAPE(Tvp,Tve,Height)
    CAPE_yr_3K_RH = np.append(CAPE_yr_3K_RH,temp)
ave_Q4 = [np.mean(CAPE_yr_3K_RH)]*1460
    
# Draw
os.chdir(main_path)

# Fig.1
plt.figure(num = 1, figsize = (24,15))
    # Plot
plt.plot(Time,CAPE_yr,label = 'CAPE')
plt.plot(Time,ave_Q2,linewidth=4,linestyle = '--')
    # X-axis
plt.xlabel('Months',fontsize = 24)
plt.xlim(1,12)
plt.xticks([1,2,3,4,5,6,7,8,9,10,11,12],['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'],fontsize = 14)
    # Y-axis
plt.ylabel('Energy [J/Kg]',fontsize = 24)
plt.yticks(fontsize=14)
    # Legend & Title
plt.legend(prop={'size': 24})
plt.title('Series of the 6-hourly CAPE',fontsize = 34)
    # Save & Show
plt.savefig('CA8.5-Q2.png',dpi = 300)
plt.show()

# Fig.2
plt.figure(num = 2, figsize = (24,15))
    # Plot
plt.plot(Time,CAPE_yr_3K,label = 'CAPE , +3K')
plt.plot(Time,ave_Q3,linewidth=4,linestyle = '--')
    # X-axis
plt.xlabel('Months',fontsize = 24)
plt.xlim(1,12)
plt.xticks([1,2,3,4,5,6,7,8,9,10,11,12],['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'],fontsize = 14)
    # Y-axis
plt.ylabel('Energy [J/Kg]',fontsize = 24)
plt.ylim(10000,24000)
plt.yticks(fontsize=14)
    # Legend & Title
plt.legend(prop={'size': 24})
plt.title('Series of the 6-hourly CAPE, surface warmer',fontsize = 34)
    # Save & Show
plt.savefig('CA8.5-Q3.png',dpi = 300)
plt.show()

# Fig.3
plt.figure(num = 3, figsize = (24,15))
    # Plot
plt.plot(Time,CAPE_yr_3K_RH,label = 'CAPE, +3K, RH')
plt.plot(Time,ave_Q4,linewidth=4,linestyle = '--')
    # X-axis
plt.xlabel('Months',fontsize = 24)
plt.xlim(1,12)
plt.xticks([1,2,3,4,5,6,7,8,9,10,11,12],['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'],fontsize = 14)
    # Y-axis
plt.ylabel('Energy [J/Kg]',fontsize = 24)
plt.ylim(17000,27000)
plt.yticks(fontsize=14)
    # Legend & Title
plt.legend(prop={'size': 24})
plt.title('Series of the 6-hourly CAPE, surface warmer , same RH',fontsize = 34)
    # Save & Show
plt.savefig('CA8.5-Q4.png',dpi = 300)
plt.show()
    
    
    


        