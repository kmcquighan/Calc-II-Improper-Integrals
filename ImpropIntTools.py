# -*- coding: utf-8 -*-
"""
Copyright Kelly McQuighan 2016

These tools can be used to visualize different improper integrals.
"""

from matplotlib import pyplot as plt
import numpy as np
from numpy import *
import scipy.integrate as inte
import matplotlib as mpl
mpl.rcParams['font.size'] = 17
import scipy.interpolate as interp

def plots(func, a,b,n,method,ax):
    
    ax.axvline(0.,color='#666666',linewidth=1)
    ax.axhline(0.,color='#666666',linewidth=1)
    if (a>0):
        xlarge = np.linspace(0.,1.03*b,1000)
    elif (b<0):
        xlarge = np.linspace(1.03*a,0.,1000)
    else:
        xlarge = np.linspace(1.03*a,1.03*b,1000)
    flarge = func(xlarge)
    ax.plot(xlarge,flarge,'b', linewidth=5)
    
    ax.set_xlim([xlarge[0], xlarge[999]])
    smallest = np.min(flarge)
    largest = np.max(flarge)
    
    dx = (b-a)/n
    
    xs = np.linspace(a,b,n+1)
    fxs = func(xs)
    
    if method.lower()=='left':
        for i in range(n):
            points = [[xs[i], 0], [xs[i], fxs[i]], [xs[i+1], fxs[i]], [xs[i+1],0]]
            poly = plt.Polygon(points, fc='g', edgecolor='g', alpha=0.3, linewidth=3)
            ax.add_patch(poly)
    elif method.lower()=='right':
        for i in range(n):
            points = [[xs[i], 0], [xs[i], fxs[i+1]], [xs[i+1], fxs[i+1]], [xs[i+1],0]]
            poly = plt.Polygon(points, fc='g', edgecolor='g', alpha=0.3, linewidth=3)
            ax.add_patch(poly)
    elif method.lower()=='midpoint':
        x = np.linspace(a+dx/2.,b-dx/2.,n)
        fx = func(x)
        for i in range(n):
            points = [[xs[i], 0], [xs[i], fx[i]], [xs[i+1], fx[i]], [xs[i+1],0]]
            poly = plt.Polygon(points, fc='g', edgecolor='g', alpha=0.3, linewidth=3)
            ax.add_patch(poly)
    elif method.lower()=='trapezoid':
        for i in range(n):
            points = [[xs[i], 0], [xs[i], fxs[i]], [xs[i+1], fxs[i+1]], [xs[i+1],0]]
            poly = plt.Polygon(points, fc='g', edgecolor='g', alpha=0.3, linewidth=3)
            ax.add_patch(poly)  
    elif method.lower()=='simpson':
        # note: this implementation keeps the number of grid points the same         
        for i in range(0,n,2):
            lag = interp.lagrange([xs[i], xs[i+1], xs[i+2]], [fxs[i], fxs[i+1], fxs[i+2]])
            section = np.linspace(xs[i], xs[i+2], 100)
            fsec = lag(section)
            ax.fill_between(section,fsec, facecolor='g', edgecolor='g', alpha=0.3, linewidth=3)
        
        x_mid = np.linspace(a+dx,b-dx,n/2)
        vert = np.ones(100)   
        for the_x in x_mid:
            ax.plot(the_x*vert,np.linspace(0,func(the_x),100),'g--', linewidth=3, alpha=0.5)
        
    else:
        print ('ERROR: You have not specified a valid method. Please check for typos.')
    
    if smallest>0:
        ax.set_ylim([0,1.1*largest])
    elif largest<0:
        ax.set_ylim([1.1*smallest,0])
    else:
        ax.set_ylim([1.1*smallest, 1.1*largest])
    ax.set_xlabel('x')
    ax.set_ylabel('f(x)')
    
def evalArea(func,a,b,n, method):
    
    dx = (b-a)/n
    
    if method.lower()=='left':
        x = np.linspace(a,b-dx,n)
        fx = func(x)
        area = np.sum(fx)*dx
    elif method.lower()=='right':
        x = np.linspace(a+dx,b,n)
        fx = func(x)
        area = np.sum(fx)*dx
    elif method.lower()=='midpoint':
        x = np.linspace(a+dx/2.,b-dx/2.,n)
        fx = func(x)
        area = np.sum(fx)*dx
    elif method.lower()=='trapezoid':
        x = np.linspace(a+dx,b-dx,n-1)
        fx = func(x)
        area = dx*(0.5*func(a)+0.5*func(b)+np.sum(fx))
    elif method.lower()=='simpson':
        x_mid = np.linspace(a+dx,b-dx,n/2)
        x_trap = np.linspace(a+2*dx,b-2*dx,n/2-1)
        fx_mid = func(x_mid)
        fx_trap = func(x_trap)
        area = dx/3.*(4*np.sum(fx_mid)+func(a)+func(b)+2*np.sum(fx_trap))
    else:
        print ('ERROR: You have not specified a valid method. Please check for typos.')
    
    return area

def plotIntegralApprox(f,a,b,n):
    
    a = float(a)
    b = float(b)
    assert (a<=b) # "ERROR! a must be less than b!"
    nmax = 100

    func = eval("lambda x: " + f)
    #I = inte.quad(func, a, b)[0]
        
    fig = plt.figure(figsize=(20, 6))
       
    ax1 = fig.add_subplot(1,2,1)
    ax2 = fig.add_subplot(1,2,2)
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    
    plots(func,a,b,n,"left",ax1)
    
    n_all = range(nmax+1)
    ax2.set_xlim([0,nmax])
    
    idx = n
    
    I = np.zeros(nmax+1)
    for i in range(1,nmax+1):
        I[i] = evalArea(func,a,b,n_all[i],"left")
    ax2.plot(n_all[1:idx+1],I[1:idx+1],'ro',markersize=10)
    ax2.plot(n_all[idx], I[idx],'go', markersize=13)
    
    if (min(I)==0.):
        ax2.set_ylim([0,1.05*max(I)])
    elif (max(I)==0.):
        ax2.set_ylim([1.05*min(I),0])
    else:
        ax2.set_ylim([1.05*min(I),1.05*max(I)])
   
    #ax1.axhline(0.,linewidth=1,color='k')
    ax2.axhline(0.,linewidth=1,color='k')
    
    ax1.set_xlabel('x', fontsize=36)
    ax1.set_title('f(x)', fontsize=36)
    ax2.set_xlabel('n', fontsize=36)
    ax2.set_title(r'$\sum_{i=1}^n f(x_i)\Delta x$',fontsize=36, y=1.1)
    
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    plt.suptitle('f(x) = '+f+ ' for $x\in$ [%.0f, %.0f]' %(a,b), fontsize=36, y=1.1)
    mpl.rc('xtick', labelsize=20) 
    mpl.rc('ytick', labelsize=20)
    
def plotINFInterval(f,a,b):
    
    a = float(a)
    assert (a<=b) # "ERROR! a must be less than b!"
    amin = 0.2
    bmax = 100

    func = eval("lambda x: " + f)
        
    fig = plt.figure(figsize=(20, 6))
       
    ax1 = fig.add_subplot(1,2,1)
    ax2 = fig.add_subplot(1,2,2)
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    
    x = np.linspace(a,bmax,1000)
    x_all = np.linspace(amin,1.03*bmax,1000)
    y = func(x)
    y_all = func(x_all)
    ax1.set_xlim([0,1.03*bmax])
    ax2.set_xlim([0,1.03*bmax])
    
    #dx = (bmax-a)/1000.
    idx = np.where(x>=b)[0][0]
    
    ax1.plot(x_all,y_all,'b',linewidth=5)
    ax1.fill_between(x[:idx+1],y[:idx+1], facecolor='g', edgecolor='g', alpha=0.3, linewidth=3)
    
    if (min(y_all)>0.):
        ax1.set_ylim([0,1.03*max(y_all)])
    elif (max(y_all)<0.):
        ax1.set_ylim([1.03*min(y_all),0.])
    else:
        ax1.set_ylim([1.03*min(y_all),1.03*max(y_all)])
    
    I = np.zeros(1000)
    for i in range(1000):
        I[i] = inte.quad(func,a,x[i])[0]
    ax2.plot(x[:idx],I[:idx],'r',linewidth=5)
    ax2.plot(x[idx], I[idx],'go', markersize=13)
    
    if (min(I)==0.):
        ax2.set_ylim([0.,1.03*max(I)])
    elif (max(I)==0.):
        ax2.set_ylim([1.03*min(I),0.])
    else:
        ax2.set_ylim([1.03*min(I),1.03*max(I)])
   
    ax1.axhline(0.,linewidth=1,color='k')
    ax2.axhline(0.,linewidth=1,color='k')
    
    ax1.set_xlabel('x', fontsize=36)
    ax1.set_title('f(x)', fontsize=36)
    ax2.set_xlabel('b', fontsize=36)
    ax2.set_title(r'$\int_a^b f(x)dx$',fontsize=36, y=1.1)
    
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    plt.suptitle('f(x) = '+f+ ' for $x\in$ [%.0f, %.0f]' %(a,b), fontsize=36, y=1.1)
    mpl.rc('xtick', labelsize=20) 
    mpl.rc('ytick', labelsize=20)

def plotINFIntegrand(g,m,b):
    
    b = float(b)
    a = 10**m
    mmin = -6
    
    func = eval("lambda x: " + g)
        
    fig = plt.figure(figsize=(20, 12))
       
    ax1 = fig.add_subplot(1,2,1)
    ax2 = fig.add_subplot(2,2,2)
    ax3 = fig.add_subplot(2,2,4)
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    
    #x = np.linspace(amin,b,1000)
    x = np.logspace(mmin,np.log10(b),1000)
    x_all = np.linspace(0.05,1.03*b,1000)
    y = func(x)
    y_all = func(x_all)
    ax1.set_xlim([0,1.03*b])
    ax2.set_xlim([-0.01,1.03*b])
    ax3.set_xlim([1.01*mmin,np.log10(1.2*b)])
    
    idx = np.where(x>=a)[0][0]
    
    ax1.plot(x_all,y_all,'b',linewidth=5)
    ax1.fill_between(x[idx:],y[idx:], facecolor='g', edgecolor='g', alpha=0.3, linewidth=3)
    if (min(y_all)>0.):
        ax1.set_ylim([0,max(y_all)])
    elif (max(y_all)<0.):
        ax1.set_ylim([min(y_all),0.])
    else:
        ax1.set_ylim([min(y_all),max(y_all)])
    
    I = np.zeros(1000)
    for i in range(1000):
        I[i] = inte.quad(func,x[i],b)[0]
    ax2.plot(x[idx:],I[idx:],'r',linewidth=5)
    ax2.plot(x[idx], I[idx],'go', markersize=13)
    
    ax3.plot(np.log10(x[idx:]),I[idx:],'r',linewidth=5)
    ax3.plot(np.log10(x[idx]), I[idx],'go', markersize=13)

    if (min(I)==0.):
        ax2.set_ylim([0.,1.03*max(I)])
        ax3.set_ylim([0.,1.03*max(I)])
    elif (max(I)==0.):
        ax2.set_ylim([1.03*min(I),0.])
        ax3.set_ylim([1.03*min(I),-0.])
    else:
        ax2.set_ylim([1.03*min(I),1.03*max(I)])
        ax3.set_ylim([1.03*min(I),1.03*max(I)])

    ax1.axhline(0.,linewidth=1,color='k')
    ax2.axhline(0.,linewidth=1,color='k')
    ax3.axhline(0.,linewidth=1,color='k')

    ax1.set_xlabel('x', fontsize=36)
    ax1.set_title('g(x)', fontsize=36)
    ax2.set_xlabel(r'$a$', fontsize=36)
    ax2.set_title(r'$\int_a^b g(x)dx$',fontsize=36, y=1.1)
    ax3.set_xlabel(r'$\log_{10}(a)$', fontsize=36)
    
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    plt.suptitle('g(x) = '+g+ ' for $x\in$ [%.2e, %.2f]' %(a,b), fontsize=36, y=1.1)
    mpl.rc('xtick', labelsize=20) 
    mpl.rc('ytick', labelsize=20)