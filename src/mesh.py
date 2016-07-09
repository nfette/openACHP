# -*- coding: utf-8 -*-
"""
Created on Fri Jul 08 21:01:32 2016

@author: nfette
"""

import numpy as np
import scipy.interpolate

def mesh1D(fun,x0,x1,interp):
    x = np.linspace(x0,x1)
    y = fun(x)
    i = interp(x,y)
    return x,y,i
    
def smooth(fun,x,y,i,rhomin,n):
    # The y-error on an interval Dx is c^2 = (1/2 * Dx)^2 * |d2y/dx2|.
    # To keep constant error over the curve, this should be constant.
    # Define density rho = 1 / Dx. Then
    # rho = |d2y/dx2| / 2c
    dd = i.derivative(2)
    rho = np.vectorize(lambda(x):max(np.abs(dd(x)),rhomin))(x)
    density = scipy.interpolate.Akima1DInterpolator(x,rho)
    intd = density.antiderivative()
    yy = intd(x)
    yy *= (n - 1) / (yy.max() - yy[0])
    newmap = scipy.interpolate.Akima1DInterpolator(yy,x)
    xnew = newmap(np.arange(n))
    ynew = fun(xnew)
    return xnew, ynew
    
if __name__ == "__main__":
    f= np.exp
    x,y,i=mesh1D(f,0,2,scipy.interpolate.Akima1DInterpolator)
    xnew,ynew = smooth(f,x,y,i,0.01,50)
    
    
    import matplotlib.pyplot as plt
    plt.plot(x,y)
    #plt.plot(x,yy)
    plt.plot(xnew,ynew,'.')

    
    