# -*- coding: utf-8 -*-
"""
Created on Fri Mar 04 09:54:22 2016

@author: nfette
"""
import CoolProp.CoolProp as CP
CP.PropsSI('P','T',373,'Q',0,"INCOMP::LiBr[0.6]")
amm = lambda(x):'REFPROP::water[{}]&ammonia[{}]'.format(1-x,x)
f1 = amm(0.302385)
print f1
print CP.PropsSI('D','T',300,'P',101325,f1)

print CP.PropsSI('C','T',300,'P',101325,'INCOMP::LiBr[0.23]')
Psat = CP.PropsSI('P','T',300,'Q',0,'INCOMP::LiBr[0.23]')
print Psat
