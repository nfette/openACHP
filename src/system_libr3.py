# -*- coding: utf-8 -*-
"""
Created on Sat Apr 30 03:26:36 2016

This file solves a thermal system with a chiller.

@author: nfette
"""

import libr3

class system1:
    def __init__(self,UA):
        self.UA = UA
        self.chiller = libr3.ChillerLiBr1()
        self.chiller.iterate1()
    def calcUA(self):
        pass
    def optimize(self):
        pass
    
if __name__ == "__main__":
    sys = system1(100)
    sys.optimize()
    print sys.chiller