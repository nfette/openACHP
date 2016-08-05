# -*- coding: utf-8 -*-
"""
Created on Sun Dec 14 19:10:32 2014

@author: nfette

External procedures are documented in EES help, or see external_procedures.htm.
You can import a C-type DLL using Python ctypes module. The EES interface
uses one function with several modes to return info about the function itself,
and then the values you are evaluating. This module creates an interface
to set up the function calls in more Pythonic syntax.

I'm a little confused. enum works in Python(x,y) for Windows, which I thought
was strictly python 2.7. However, on command line linux, enum works only
with python3 (supported by enum documentation).
"""

from __future__ import print_function
import ctypes
import enum
import os.path

@enum.unique
class mode(enum.IntEnum):
    getCallFormat = -1
    getInputUnits = -2
    getOutputUnits = -3
    call = 0

#class EesExternalProcedure (ctypes.Structure):
class EesParamRec (ctypes.Structure):
    pass

EesParamRec._fields_ = [("value",ctypes.c_double),
                        ("next",ctypes.POINTER(EesParamRec))]

def EesParamRec2List(rec):
    res = []
    res.append(rec.value)
    while (rec.next):
        rec = rec.next.contents
        res.append(rec.value)
    return res

def List2EesParamRec(listin):
    # Create all the nodes, defaulting to null next pointer
    recs = map(EesParamRec, listin)
    # Now set the next pointers
    for i in range(len(recs)-1):
        recs[i].next = ctypes.pointer(recs[i+1])
    if not recs:
        # "The procedure will have one or more inputs." You messed up.
        pass
    return recs[0]

EesStringData = ctypes.c_char * 256
#EesStringData = ctypes.c_char_p # does not work

class EES_DLP:
    def __init__(self, path):
        self.path = path
        self.mydll = ctypes.WinDLL(self.path)
        #myDLL=ctypes.cdll.LoadLibrary(self.path)
        # The function has the same name as the file, per EES documentation
        self.name = os.path.splitext(os.path.basename(path))[0].upper()
        #print(self.name)
        self.func = self.mydll[self.name]
        self.func.argtypes=[EesStringData, ctypes.POINTER(ctypes.c_int),
                       ctypes.POINTER(EesParamRec), ctypes.POINTER(EesParamRec)]

    def wrapper(self, s, intmode, inarglist):
        strdata = EesStringData(" ")
        strdata.raw = "{:256}".format(s)
        intmode = ctypes.c_int(intmode)
        inargs = List2EesParamRec(inarglist)
        outargs = List2EesParamRec(range(8))
        self.func(strdata, ctypes.byref(intmode),
            ctypes.byref(inargs), ctypes.byref(outargs))
        outarglist = EesParamRec2List(outargs)
        return strdata.value, outarglist
        
    def getCallFormat(self,S="",inarglist=[0]):
        callFormat,_ = self.wrapper(S,mode.getCallFormat,inarglist)
        invars,outvars = callFormat.split('(')[1].split(')')[0].split(':')
        invars = map(str.strip,invars.split(','))
        outvars = map(str.strip,outvars.split(','))
        return callFormat, invars, outvars
    def getInputUnits(self,S="",inarglist=[0]):
        return self.wrapper(S,mode.getInputUnits,inarglist)[0].split(',')
    def getOutputUnits(self,S="",inarglist=[0]):
        return self.wrapper(S,mode.getOutputUnits,inarglist)[0].split(',')
    def call(self,S,inarglist):
        return self.wrapper(S,mode.call,inarglist)
        
if __name__ == "__main__":
    myDLL = EES_DLP(r'C:\EES32\Userlib\EES_System\nh3h2o.dlp')
    
    callFormat, invars, outvars = myDLL.getCallFormat()
    print(callFormat)
    print("Invars, outvars:{},{}".format(invars,outvars))
    
    print("Inputs:")
    inarglist = [123, 450, 10, 0.5]
    inunits = myDLL.getInputUnits(inarglist=inarglist)
    for p in zip(invars, inarglist, inunits): print(p)
    
    print("Outputs:")
    outunits = myDLL.getOutputUnits()
    
    inarglist = [123, 450, 10, 0.5]
    s0,outarglist = myDLL.call("", inarglist)
    #print outarglist
    for p in zip(outvars, outarglist, outunits): print(p)
    
    callFormat, invars, outvars = myDLL.getCallFormat()
    print(callFormat)
    print("Invars, outvars:{},{}".format(invars,outvars))
    
    print("Inputs:")
    inarglist = [123, 450, 10, 0.5]
    inunits = myDLL.getInputUnits(inarglist=inarglist)
    for p in zip(invars, inarglist, inunits): print(p)
    
    print("Outputs:")
    outunits = myDLL.getOutputUnits()
    
    inarglist = [123, 450, 10, 0.5]
    s0,outarglist = myDLL.call("", inarglist)
    #print outarglist
    for p in zip(outvars, outarglist, outunits): print(p)
