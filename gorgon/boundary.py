# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 16:52:03 2016

@author: lm2410
"""
import scipy.constants as const
import scipy.signal as sig
import numpy as np

def find_BS(x,rho,drho_min=0.05,debug=False):
    x_BS,i_BS = find_BS_cell(x,rho,drho_min,debug=debug)
    if(i_BS==-99):
        return x[0]
    else:
        return BS_sg_parabolic(x,rho,i_BS,debug=debug)

def gradient(y,x):
#     dy = np.diff(y)
#     xc = 0.5*(x[1:]+x[:-1])
    
    dy = y[2:]-y[1:-1]
    dy2 = y[2:]-y[:-2]
    dy = (27*dy - dy2)/24
    xc = 0.5*(x[2:]+x[1:-1])
    
    #drho = np.gradient(rho)
    #xc = x
    return dy,xc

def find_BS_cell(x,rho,drho_min,debug=False):
    from scipy.signal import argrelmax
    if(debug):
        import matplotlib.pyplot as plt
        
    ''' Takes values along the x line '''
    # Calculate mass density gradient
    drho,xc = gradient(rho,x)
    
    ''' Normalise drho'''
    C = drho.max()
    drho = drho/C
    
    ''' Plot '''
    if(debug):
        fig, ax = plt.subplots(2, 1, sharex=True)
        ax[0].plot(x, rho, '.-')
        ax[1].plot(xc, drho, '.-')
    
    ''' Find negative peaks '''
    pks = argrelmax(drho)[0]
    if(pks.shape[0] == 0):
        print('BS: No peaks found')
        return np.nan,-99
    
    if(debug):
        for i in pks:
            for axi in ax:
                axi.axvline(x=xc[i], linestyle='dashed', color='k')
            
    ''' Sort by descending distance '''
    s = np.argsort(xc[pks])
    pks = pks[s]
    
    ''' Remove insignificant peaks'''
    el = drho[pks]>drho_min
    pks = pks[el]
    
    if(pks.shape[0] == 0):
        print('BS: No significant peaks found')
        return np.nan,-99
    
    if(debug):
        for i in pks:
            for axi in ax:
                axi.axvline(x=xc[i], linestyle='dashed', color='b')
        ax[1].axhline(y=drho_min, linestyle='dashed', color='k')
    
    ''' Pick furthest shock '''
    i_BS = pks[0]
    
    if(debug):
        for axi in ax:
            axi.axvline(x=xc[i_BS], color='r')
    
        print(pks[0],xc[pks[0]])    
    return xc[i_BS],i_BS

def BS_sg_parabolic(x,rho,i_BS,debug=False):
    if(debug):
        import matplotlib.pyplot as plt
    ''' Calculate gradient again'''
    dy,xc = gradient(rho,x)
    
    if(debug):
        ''' Plot '''
        fig, ax = plt.subplots(2, 1, sharex=True)
        ax[0].plot(x,rho,'.-')
        ax[0].set_ylabel('f')
        ax[1].plot(xc,dy,'.-')
        ax[1].set_ylabel('$\partial_x$ f')
    
    ''' Select points around x_BS_grid'''
    di = 1
    
    dyi = dy[i_BS-di:i_BS+di+1]
    xci = xc[i_BS-di:i_BS+di+1]
    
    if(debug):
        ax[1].plot(xci,dyi,'r.')
    
    ''' Parabolic estimate '''
    x2 = xci*xci
    
    a = dyi[1]*(xci[2]-xci[0])-dyi[0]*(xci[2]-xci[1])-dyi[2]*(xci[1]-xci[0])
    a /= ( x2[0]*(xci[1]-xci[2]) - x2[2]*(xci[1]-xci[0]) - x2[1]*(xci[0]-xci[2]) )

    b = dyi[1]-dyi[0]+a*(x2[0]-x2[1])
    b /= (xci[1]-xci[0])
    
    x_BS = -0.5*b/a
    
    if(debug):
        ''' Plot '''
        for axi in ax:
            axi.axvline(x=x_BS,linestyle='dashed',color='r')
            axi.axvline(x=xci[1],linestyle='dashed',color='g')
        xi = np.linspace(xci[0],xci[-1])
        c = dyi[0] - a*xci[0]**2-b*xci[0]
        dyi = a*xi**2+b*xi+c
        ax[1].plot(xi,dyi,'r')
        ax[0].set_xlim(xc[i_BS-5],xc[i_BS+5])
        plt.show()
        
    return x_BS

def find_MP(x, j, rho, debug=False):
    
    dx = x[1]-x[0]
    if debug:
        print('Input array shapes: x {}, j {}, rho {}'.format(x.shape, j.shape, rho.shape))
    drho = np.gradient(rho, dx)
    
    # Find Peaks
    pks = sig.argrelmax(j)[0]
    if debug:
        print('No. of peaks in j: {}'.format(pks.shape))
    
    # Ensure mass decreases over MP
    el_rho = drho[pks]>=0
    pks = pks[el_rho]
    if debug:
        print('No. of peaks after rho test {}'.format(pks.shape))
    
    
    # Choose largest current density
    #s = np.argsort(j[pks])
    s = np.argsort(drho[pks])
    pks = pks[s]
    
    i_MP = pks[-1]
        
    # Estimate sub-grid position using parabola
    x_MP = parabolic_estimate(x[i_MP-1:i_MP+2], j[i_MP-1:i_MP+2])
    j_MP = j[i_MP]
    
    return x_MP, j_MP
	
def find_MP_old(x, j, rho, B, P_T, x_BS = None, show_errors=True):    
    #j = line['j']
    #rho = line['rho']
    #B = line['B']
    #P_T = line['P_T']
        
    # Calculate 
    i = sig.argrelmax(j)[0]

    if len(i) == 0:
        return (np.nan,np.nan)
            
    ''' Sort by descending distance '''
    s = np.argsort(x[i])
    i = i[s]
        
    #''' Sort by largest j peak '''
    #s = np.argsort(j[i])
    #s = s[::-1]
    #i = i[s]
    
    ''' Remove values near edges '''
    el = np.logical_and(2<i,i<len(j)-2)
    i = i[el]
    
    ''' Test other properties:
         - B increases across the magnetopause
         - rho, mom decreases across the magnetopause '''
         
    i0 = i-2
    i1 = i+2
         
    el_B = B[i0]<B[i1]
    el_rho = rho[i0]>rho[i1]
    el_Pbal = P_T[i0]>P_T[i1]
    el_rho2 = rho[i]>1e-23
    
    if( x_BS is None):  # Remove peaks outside bow shock
        el_BS = True
    else:
        el_BS = np.abs(x[i])<np.abs(x_BS-1.) 
    
    el = np.logical_and(el_B,el_rho)
    el = np.logical_and(el,el_rho2)
    el = np.logical_and(el,el_Pbal)
    el = np.logical_and(el,el_BS)
    
    ''' Store value '''
    i = i[el]
    
    if len(i) == 0:
        if(show_errors): print('Error: No MP found')
        return np.nan #(np.nan,np.nan)
        
    i = i[0]    #If multiple values, take the largest gradient element: i[0] (since its sorted)
    
    ''' Sub-grid BS position using parabola '''
    x_MP_SG = parabolic_estimate(x[i-1:i+2],j[i-1:i+2])
    x_MP = x[i]
        
    return x_MP_SG
    
def parabolic_estimate(x,y):
    
    x2 = x*x
    
    a = y[1]*(x[2]-x[0])-y[0]*(x[2]-x[1])-y[2]*(x[1]-x[0])
    a /= ( x2[0]*(x[1]-x[2]) - x2[2]*(x[1]-x[0]) - x2[1]*(x[0]-x[2]) )

    b = y[1]-y[0]+a*(x2[0]-x2[1])
    b /= (x[1]-x[0])
    
    return -b/(2*a)
    
    