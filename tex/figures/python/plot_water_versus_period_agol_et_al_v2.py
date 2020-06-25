#!/usr/bin/env python
#-*- coding:Utf-8 -*-

from	scipy.ndimage.filters                import	gaussian_filter1d	
from	scipy.ndimage.filters                import	convolve1d
from	scipy.ndimage.filters                import	uniform_filter1d
from numpy import *		
import  numpy                 as        np
import  matplotlib.pyplot     as        mpl
import  math
from math import log
from matplotlib.colors import LogNorm
import pylab
from matplotlib.patches import Ellipse


mpl.figure(1,figsize=(10,6))

ax1 = mpl.gca() # period 'axis'
ax2 = ax1.twiny() # semi-major 'axis'

# Period of planets (in days)
P_planets=[1.510826,2.421938,4.049218,6.101013,9.207541,12.352445,18.772863]

#CMF50

CMF50_T1b=[3e-5,3e-4,1e-3]
CMF50_T1c=[3.e-5,2.e-4,6.e-4]
CMF50_T1d=[5.e-5,3e-5]
CMF50_T1e=[6.3e-2,8.01e-2,10.29e-2]
CMF50_T1f=[8.71e-2,10.68e-2,12.64e-2]
CMF50_T1g=[11.04e-2,12.85e-2,15.10e-2]
CMF50_T1h=[7.17e-2,11.94e-2,16.71e-2]

ax1.errorbar(P_planets[0]-0.15,CMF50_T1b[1],yerr=[[CMF50_T1b[1]-CMF50_T1b[0]],[CMF50_T1b[2]-CMF50_T1b[1]]],color='silver',alpha=1.0,lw=2)
ax1.errorbar([-1],[-1],[1], color='silver',alpha=1.0,label=r'50$\%$',lw=2,marker='o',linestyle="none")
ax1.plot(P_planets[0]-0.15,CMF50_T1b[1],'o',color='silver')

ax1.errorbar(P_planets[1]-0.15,CMF50_T1c[1],yerr=[[CMF50_T1c[1]-CMF50_T1c[0]],[CMF50_T1c[2]-CMF50_T1c[1]]],color='silver',alpha=1.0,lw=2)
ax1.plot(P_planets[1]-0.15,CMF50_T1c[1],'o',color='silver')

ax1.errorbar(P_planets[2]-0.15,CMF50_T1d[0],yerr=[[CMF50_T1d[0]-5.e-6],[0]],uplims=False,lolims=True,color='silver',alpha=1.0,lw=2)

ax1.errorbar(P_planets[3]-0.15,CMF50_T1e[1],yerr=[[CMF50_T1e[1]-CMF50_T1e[0]],[CMF50_T1e[2]-CMF50_T1e[1]]],color='silver',alpha=1.0,lw=2)
ax1.plot(P_planets[3]-0.15,CMF50_T1e[1],'o',color='silver')

ax1.errorbar(P_planets[4]-0.15,CMF50_T1f[1],yerr=[[CMF50_T1f[1]-CMF50_T1f[0]],[CMF50_T1f[2]-CMF50_T1f[1]]],color='silver',alpha=1.0,lw=2)
ax1.plot(P_planets[4]-0.15,CMF50_T1f[1],'o',color='silver')

ax1.errorbar(P_planets[5]-0.15,CMF50_T1g[1],yerr=[[CMF50_T1g[1]-CMF50_T1g[0]],[CMF50_T1g[2]-CMF50_T1g[1]]],color='silver',alpha=1.0,lw=2)
ax1.plot(P_planets[5]-0.15,CMF50_T1g[1],'o',color='silver')

ax1.errorbar(P_planets[6]-0.15,CMF50_T1h[1],yerr=[[CMF50_T1h[1]-CMF50_T1h[0]],[CMF50_T1h[2]-CMF50_T1h[1]]],color='silver',alpha=1.0,lw=2)
ax1.plot(P_planets[6]-0.15,CMF50_T1h[1],'o',color='silver')


#CMF32

CMF32_T1b=[5.e-5,3e-5]
CMF32_T1c=[2.e-5,1e-5]
CMF32_T1d=[1.e-5,0.5e-5]
CMF32_T1e=[1.43/100.,2.88/100.,4.31/100.]
CMF32_T1f=[3.05/100.,4.39/100.,6.18/100.]
CMF32_T1g=[4.59/100.,6.23/100.,8.29/100.]
CMF32_T1h=[2.65/100.,5.97/100.,10.60/100.]

ax1.errorbar(P_planets[0]-0.05,CMF32_T1b[0],yerr=[[CMF32_T1b[0]-5.e-6],[0]],uplims=False,lolims=True,color='steelblue',alpha=1.0,lw=2)
ax1.errorbar([-1],[-1],[1], color='steelblue',alpha=1.0,label=r'32.5$\%$',lw=2,marker='o',linestyle="none")

ax1.errorbar(P_planets[1]-0.05,CMF32_T1c[0],yerr=[[CMF32_T1c[0]-5.e-6],[0]],uplims=False,lolims=True,color='steelblue',alpha=1.0,lw=2)

ax1.errorbar(P_planets[2]-0.05,CMF32_T1d[0],yerr=[[CMF32_T1d[0]-5.e-6],[0]],uplims=False,lolims=True,color='steelblue',alpha=1.0,lw=2)

ax1.errorbar(P_planets[3]-0.05,CMF32_T1e[1],yerr=[[CMF32_T1e[1]-CMF32_T1e[0]],[CMF32_T1e[2]-CMF32_T1e[1]]],color='steelblue',alpha=1.0,lw=2)
ax1.plot(P_planets[3]-0.05,CMF32_T1e[1],'o',color='steelblue')

ax1.errorbar(P_planets[4]-0.05,CMF32_T1f[1],yerr=[[CMF32_T1f[1]-CMF32_T1f[0]],[CMF32_T1f[2]-CMF32_T1f[1]]],color='steelblue',alpha=1.0,lw=2)
ax1.plot(P_planets[4]-0.05,CMF32_T1f[1],'o',color='steelblue')

ax1.errorbar(P_planets[5]-0.05,CMF32_T1g[1],yerr=[[CMF32_T1g[1]-CMF32_T1g[0]],[CMF32_T1g[2]-CMF32_T1g[1]]],color='steelblue',alpha=1.0,lw=2)
ax1.plot(P_planets[5]-0.05,CMF32_T1g[1],'o',color='steelblue')

ax1.errorbar(P_planets[6]-0.05,CMF32_T1h[1],yerr=[[CMF32_T1h[1]-CMF32_T1h[0]],[CMF32_T1h[2]-CMF32_T1h[1]]],color='steelblue',alpha=1.0,lw=2)
ax1.plot(P_planets[6]-0.05,CMF32_T1h[1],'o',color='steelblue')

#CMF25

CMF25_T1b=[1.e-5,0.5e-5]
CMF25_T1c=[1.e-5,0.5e-5]
CMF25_T1d=[1.e-5,0.5e-5]
CMF25_T1e=[2.51e-2]
CMF25_T1f=[1.17e-2,2.35e-2,3.91e-2]
CMF25_T1g=[2.56e-2,3.65e-2,5.47e-2]
CMF25_T1h=[0.61e-2,3.63e-2,7.87e-2]

ax1.errorbar(P_planets[0]+0.05,CMF25_T1b[0],yerr=[[CMF25_T1b[0]-5.e-6],[0]],uplims=False,lolims=True,color='orange',alpha=1.0,lw=2)
ax1.errorbar([-1],[-1],[1], color='orange',alpha=1.0,label=r'25$\%$',lw=2,marker='o',linestyle="none")

ax1.errorbar(P_planets[1]+0.05,CMF25_T1c[0],yerr=[[CMF25_T1c[0]-5.e-6],[0]],uplims=False,lolims=True,color='orange',alpha=1.0,lw=2)

ax1.errorbar(P_planets[2]+0.05,CMF25_T1d[0],yerr=[[CMF25_T1d[0]-5.e-6],[0]],uplims=False,lolims=True,color='orange',alpha=1.0,lw=2)

ax1.errorbar(P_planets[3]+0.05,CMF25_T1e[0],yerr=[[CMF25_T1e[0]-5.e-6],[0]],uplims=False,lolims=True,color='orange',alpha=1.0,lw=2)

ax1.errorbar(P_planets[4]+0.05,CMF25_T1f[1],yerr=[[CMF25_T1f[1]-CMF25_T1f[0]],[CMF25_T1f[2]-CMF25_T1f[1]]],color='orange',alpha=1.0,lw=2)
ax1.plot(P_planets[4]+0.05,CMF25_T1f[1],'o',color='orange')

ax1.errorbar(P_planets[5]+0.05,CMF25_T1g[1],yerr=[[CMF25_T1g[1]-CMF25_T1g[0]],[CMF25_T1g[2]-CMF25_T1g[1]]],color='orange',alpha=1.0,lw=2)
ax1.plot(P_planets[5]+0.05,CMF25_T1g[1],'o',color='orange')

ax1.errorbar(P_planets[6]+0.05,CMF25_T1h[1],yerr=[[CMF25_T1h[1]-CMF25_T1h[0]],[CMF25_T1h[2]-CMF25_T1h[1]]],color='orange',alpha=1.0,lw=2)
ax1.plot(P_planets[6]+0.05,CMF25_T1h[1],'o',color='orange')

#CMF18

CMF18_T1b=[1.e-5]
CMF18_T1c=[1.e-5]
CMF18_T1d=[1.e-5]
CMF18_T1e=[0.36e-2]
CMF18_T1f=[1.65e-2]
CMF18_T1g=[0.30e-2,1.51e-2,3.15e-2]
CMF18_T1h=[5.47e-2]

ax1.errorbar(P_planets[0]+0.15,CMF18_T1b[0],yerr=[[CMF18_T1b[0]-5.e-6],[0]],uplims=False,lolims=True,color='red',alpha=1.0,lw=2)
ax1.errorbar([-1],[-1],[1], color='red',alpha=1.0,label=r'18$\%$',lw=2,marker='o',linestyle="none")

ax1.errorbar(P_planets[1]+0.15,CMF18_T1c[0],yerr=[[CMF18_T1c[0]-5.e-6],[0]],uplims=False,lolims=True,color='red',alpha=1.0,lw=2)

ax1.errorbar(P_planets[2]+0.15,CMF18_T1d[0],yerr=[[CMF18_T1d[0]-5.e-6],[0]],uplims=False,lolims=True,color='red',alpha=1.0,lw=2)

ax1.errorbar(P_planets[3]+0.15,CMF18_T1e[0],yerr=[[CMF18_T1e[0]-5.e-6],[0]],uplims=False,lolims=True,color='red',alpha=1.0,lw=2)

ax1.errorbar(P_planets[4]+0.15,CMF18_T1f[0],yerr=[[CMF18_T1f[0]-5.e-6],[0]],uplims=False,lolims=True,color='red',alpha=1.0,lw=2)

ax1.errorbar(P_planets[5]+0.15,CMF18_T1g[1],yerr=[[CMF18_T1g[1]-CMF18_T1g[0]],[CMF18_T1g[2]-CMF18_T1g[1]]],color='red',alpha=1.0,lw=2)
ax1.plot(P_planets[5]+0.15,CMF18_T1g[1],'o',color='red')

ax1.errorbar(P_planets[6]+0.15,CMF18_T1h[0],yerr=[[CMF18_T1h[0]-5.e-6],[0]],uplims=False,lolims=True,color='red',alpha=1.0,lw=2)

ax1.plot([-0.5,30.],[0.023e-2,0.023e-2],'--',color='black') # earth
ax1.text(10.5, 1e-4, 'Earth water ocean content', fontsize=15)

#ax1.plot([-0.5,30.],[0.023e-2*0.00001,0.023e-2*0.00001],'--',color='black') # venus *0.00001
#ax1.plot([-0.5,30.],[0.023e-2*0.01,0.023e-2*0.01],'--',color='black') # mars *0.01

# PERIOD AXIS
ax1.set_xlim(0.,20.)
ax1.set_ylim(1e-6,0.4)
ax1.set_ylabel(r'Water mass fraction')
ax1.set_xlabel(r'Period (days)')

# SEMI-MAJOR AXIS
T2a3_cste=(1.510826)**2./(0.01154)**3.
ax2.set_xlabel("Semi-major axis (AU)")
ax2.set_xlim(0.,20.)
ax2.set_xticks([(T2a3_cste*0.01**3.)**(1./2.),(T2a3_cste*0.02**3.)**(1./2.),(T2a3_cste*0.03**3.)**(1./2.),(T2a3_cste*0.04**3.)**(1./2.),(T2a3_cste*0.05**3.)**(1./2.),(T2a3_cste*0.06**3.)**(1./2.)])
ax2.set_xticklabels(['0.01','0.02','0.03','0.04','0.05','0.06'])

ax1.legend(title="Iron Mass Fraction",loc=2,numpoints=1)
ax1.semilogy()

mpl.savefig('figure_water_content_versus_period_agol_et_al.png',dpi=500)
mpl.show()
