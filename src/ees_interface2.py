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

class wrappedProcedure:
    def __init__(self,func):
        self.func = func
        self.func.argtypes=[EesStringData, ctypes.POINTER(ctypes.c_int),
                       ctypes.POINTER(EesParamRec), ctypes.POINTER(EesParamRec)]

    def wrapper(self, s, intmode, inarglist):
        strdata = EesStringData(" ")
        strdata.raw = "{:256}".format(s)
        intmode = ctypes.c_int(intmode)
        #print(intmode)
        inargs = List2EesParamRec(inarglist)
        outargs = List2EesParamRec(range(5))
        self.func(strdata, ctypes.byref(intmode),
            ctypes.byref(inargs), ctypes.byref(outargs))
        #print(intmode)
        #print(outargs)
        outarglist = EesParamRec2List(outargs)
        #print(strdata.value)
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
    
class wrappedFunction:
    def __init__(self,func):
        self.func = func
        self.func.argtypes=[EesStringData, ctypes.POINTER(ctypes.c_int),
                       ctypes.POINTER(EesParamRec)]
        self.func.restype=ctypes.c_double

    def wrapper(self, s, intmode, inarglist):
        strdata = EesStringData(" ")
        strdata.raw = "{:256}".format(s)
        intmode = ctypes.c_int(intmode)
        #print(intmode)
        inargs = List2EesParamRec(inarglist)
        outargs = self.func(strdata, ctypes.byref(intmode),
            ctypes.byref(inargs))
        #print(intmode)
        #print(type(outargs))
        #print(strdata.value)
        return strdata.value, outargs
        
    def getCallFormat(self,S="",inarglist=[0]):
        callFormat,_ = self.wrapper(S,mode.getCallFormat,inarglist)
        invars = callFormat.split('(')[1].split(')')[0]
        invars = map(str.strip,invars.split(','))
        return callFormat, invars
    def getInputUnits(self,S="",inarglist=[0]):
        return self.wrapper(S,mode.getInputUnits,inarglist)[0].split(',')
    def getOutputUnits(self,S="",inarglist=[0]):
        return self.wrapper(S,mode.getOutputUnits,inarglist)[0].split(',')
    def call(self,S,inarglist):
        return self.wrapper(S,mode.call,inarglist)

class EES_DLL:
    def __init__(self, path):
        self.path = path
        self.name = os.path.splitext(os.path.basename(path))[0].upper()
        #print(self.name)
        
        try:
            # This works for older libraries like LiBr.DLL        
            self.mydll = ctypes.WinDLL(self.path)
            self.setupNames()
        except:
            # This works for newer libraries like SSCLiBr.DLL
            self.mydll = ctypes.cdll.LoadLibrary(self.path)
            #self.mydll = ctypes.CDLL(self.path)
            self.setupNames()
        
        # Wrap the functions.
        self.func = {name:wrappedFunction(self.mydll[name]) for name in self.getDLFnames()}
        self.proc = {name:wrappedProcedure(self.mydll[name]) for name in self.getDLPnames()}
        
    def setupNames(self):
        # The function names are returned by some interface functions.
        for i in ['DLFNames','DLPNames','FDLNames']:
            self.mydll['DLFNames'].argtypes = [ctypes.c_char_p]
            self.mydll['DLFNames'].restype = None
        funcnames = self.getDLPnames()+self.getDLFnames()+self.getFDLnames()
        #print(funcnames)
        
    def getDLFnames(self):
        
        strdata = EesStringData(" ")
        strdata.raw = "{:256}".format("")
        self.mydll['DLFNames'](strdata)
        newnames = strdata.value.strip().split(',')
        return [n for n in newnames if len(n) > 0]
    def getDLPnames(self):
        strdata = EesStringData(" ")
        strdata.raw = "{:256}".format("")
        self.mydll['DLPNames'](strdata)
        newnames = strdata.value.strip().split(',')
        return [n for n in newnames if len(n) > 0]
    def getFDLnames(self):
        strdata = EesStringData(" ")
        strdata.raw = "{:256}".format("")
        self.mydll['FDLNames'](strdata)
        newnames = strdata.value.strip().split(',')
        return [n for n in newnames if len(n) > 0]
        
        
if __name__ == "__main__":
    # LiBr?
    LiBr_path = r'C:\EES32\Userlib\Libr\LIBR.dll'
    myDLL = EES_DLL(LiBr_path)
    qlibr = myDLL.proc['Q_LIBR']
    
    callFormat, invars, outvars = qlibr.getCallFormat()
    print(callFormat)
    print("Invars, outvars:{},{}".format(invars,outvars))
    
    inarglist = [154.5, 0.673, 62.5, 2]
    print("In values: {}".format(inarglist))
    
    inunits = qlibr.getInputUnits(inarglist=inarglist)
    #for p in zip(invars, inarglist, inunits): print(p)
    print("Inputs units: {}".format(inunits))
    
    outunits = qlibr.getOutputUnits()
    print("Output units: {}".format(outunits))

    inarglist = [154.5, 0.673, 62.5, 2]
    s0,outarglist = qlibr.call("", inarglist)
    print("Output values: {}".format(outarglist))
    #for p in zip(outvars, outarglist, outunits): print(p)

    tlibr = myDLL.func['T_LIBR']
    
    callFormat, invars = tlibr.getCallFormat()
    print(callFormat)
    print("Invars: {}".format(invars))
    
    inarglist = [100, 50, 2]
    print("In values: {}".format(inarglist))
    
    inunits = tlibr.getInputUnits(inarglist=inarglist)
    #for p in zip(invars, inarglist, inunits): print(p)
    print("Inputs units: {}".format(inunits))
    
    outunits = tlibr.getOutputUnits()
    print("Output units: {}".format(outunits))

    inarglist = [110, 50, 2]
    s0,outarglist = tlibr.call("", inarglist)
    print("Output values: {}".format(outarglist))
    
    LiBrSSC_path = r'C:\EES32\Userlib\Libr\SSCLiBr.dll'
    SSC_DLL = EES_DLL(LiBrSSC_path)
    ssclibrh = SSC_DLL.func['LiBrSSCh']
    inarglist=[20,0.5]
    inunits = ssclibrh.getInputUnits(inarglist=inarglist)
    print("Inputs units: {}".format(inunits))
    outunits = ssclibrh.getOutputUnits()
    print("Output units: {}".format(outunits))
    s0,outarglist = ssclibrh.call("",inarglist)
    print(s0)
    print(outarglist)
    