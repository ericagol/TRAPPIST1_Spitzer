#!/usr/bin/env python
#-*- coding:Utf-8 -*-

from	scipy.ndimage.filters                import	gaussian_filter1d	
from	scipy.ndimage.filters                import	convolve1d
from	scipy.ndimage.filters                import	uniform_filter1d
from    numpy import *		
import  numpy                 as        np
import  matplotlib.pyplot     as        mpl
import  math
from    math import log
from    matplotlib.colors import LogNorm
import  pylab
from    matplotlib.patches import Ellipse
from    scipy.interpolate import interp1d
from    scipy import stats  

T1b = loadtxt('../../../data/POSTERIOR_NORM_DENSITY/out_T1b_norm')
T1c = loadtxt('../../../data/POSTERIOR_NORM_DENSITY/out_T1c_norm')
T1d = loadtxt('../../../data/POSTERIOR_NORM_DENSITY/out_T1d_norm')
T1e = loadtxt('../../../data/POSTERIOR_NORM_DENSITY/out_T1e_norm')
T1f = loadtxt('../../../data/POSTERIOR_NORM_DENSITY/out_T1f_norm')
T1g = loadtxt('../../../data/POSTERIOR_NORM_DENSITY/out_T1g_norm')
T1h = loadtxt('../../../data/POSTERIOR_NORM_DENSITY/out_T1h_norm')


mpl.figure(1)
resol=50

mpl.hist(T1b[:,1],resol,alpha=1.0,histtype=u'step',density=True,color='blue',lw=1.3)
mpl.hist(T1c[:,1],resol,alpha=1.0,histtype=u'step',density=True,color='darkorange',lw=1.3)
mpl.hist(T1d[:,1],resol,alpha=1.0,histtype=u'step',density=True,color='limegreen',lw=1.3)
mpl.hist(T1e[:,1],resol,alpha=1.0,histtype=u'step',density=True,color='red',lw=1.3)
mpl.hist(T1f[:,1],resol,alpha=1.0,histtype=u'step',density=True,color='purple',lw=1.3)
mpl.hist(T1g[:,1],resol,alpha=1.0,histtype=u'step',density=True,color='saddlebrown',lw=1.3)
mpl.hist(T1h[:,1],resol,alpha=1.0,histtype=u'step',density=True,color='fuchsia',lw=1.3)

# this is 'dummy' lines to have the label rights
mpl.plot([-2,-1],[-2,-1],alpha=1.0,color='blue',lw=1.3,label='b')
mpl.plot([-2,-1],[-2,-1],alpha=1.0,color='darkorange',lw=1.3,label='c')
mpl.plot([-2,-1],[-2,-1],alpha=1.0,color='limegreen',lw=1.3,label='d')
mpl.plot([-2,-1],[-2,-1],alpha=1.0,color='red',lw=1.3,label='e')
mpl.plot([-2,-1],[-2,-1],alpha=1.0,color='purple',lw=1.3,label='f')
mpl.plot([-2,-1],[-2,-1],alpha=1.0,color='saddlebrown',lw=1.3,label='g')
mpl.plot([-2,-1],[-2,-1],alpha=1.0,color='fuchsia',lw=1.3,label='h')


#mpl.xlim(0.75,1.15)
mpl.xlim(0.75,1.25)
mpl.ylim(0.,16.)
mpl.ylabel(r'Probability Density Function')
mpl.xlabel(r'Normalized density (by a 20$\%$ Fe, 80$\%$ MgSiO$_3$ interior)')
mpl.legend()
mpl.savefig('../plot_normalized_density_agol_et_al.pdf',dpi=300)
