# -*- coding: utf-8 -*-
"""
Created on Tue Feb 17 16:17:43 2015

@author: nfette
"""

import ammonia_props
from hw2_1 import CelsiusToKelvin as C2K

if __name__ == "__main__":
    # Example 3.1
    myprops = ammonia_props.AmmoniaProps()
    f1 = myprops.props('TPx')
    T, P, x = C2K(20.), 10., 0.5 # K, bar, dim
    
    state = f1(T, P, x)
    print state
    print "g = {}".format(state.gibbs())
    print f1.dgdxetc(T=T,P=P,x=x)
    
    print f1.getOutputUnits()
    print
    
    # Exercise 3.2
    f2 = myprops.props('TPQu')
    T, P = C2K(100.), 10. # K, bar
    
    for (Qu,label) in [(0,'sat liquid'),(1,'sat vapor'),(0.999,'almost sat vapor')]:
        print label
        state = f2(T, P, Qu)
        print state
        print "g = {}".format(state.gibbs())
        x = state.x
        print f1.dgdxetc(T=T,P=P,x=x)
        print
        