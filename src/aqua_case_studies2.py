# -*- encoding: utf-8
import numpy as np
import system_aqua1
from hw2_1 import CelsiusToKelvin as C2K
import timeit
import matplotlib.pyplot as plt

def playback(method,folder,U=100,rejectTT=np.arange(20,61,5)):
    hashes = np.zeros_like(rejectTT, dtype='<U32')
    data = {}
    QQmax = []
    fig1 = plt.figure()
    fig2 = plt.figure()

    for i,rejectT in enumerate(rejectTT):
        rT = C2K(rejectT)
        xB = np.array([400,1, rT,3, rT,5, 285,4, rT,0.15], dtype=np.double)
        # m_rich, T_evap, T_cond, T_rect, T_abs_outlet, T_gen_outlet = x        
        xC0 = np.array([0.05, 278.45, rT+7, rT+8, rT+5, 395.15])
        h = system_aqua1.hasher(xB)
        hashes[i] = h
        def func():
            
            #data[h]=system_aqua1.makeOrGetProblemForBoundary(
            #    xB,U,xC0,method="COBYLA",options=dict(rhobeg=0.01),
            #    folder='data')
            
            #data[h]=system_aqua1.makeOrGetProblemForBoundary(
            #    xB,U,xC0, folder='data2')
            
            #data[h]=system_aqua1.makeOrGetProblemForBoundary(
            #    xB,U,xC0,method="Nelder-Mead", folder='data3')
    
    #        data[h]=system_aqua1.makeOrGetProblemForBoundary(
    #            xB,U,xC0,method="Powell", folder='data4')
    
            #data[h]=system_aqua1.makeOrGetProblemForBoundary(
            #    xB,U,xC0,method="CG", folder='data5')
    
            # This looks promising
            #data[h]=system_aqua1.makeOrGetProblemForBoundary(
            #    xB,U,xC0,method="BFGS",folder='data6')
            
            # Does not work at all, because Jacobian is required
            #data[h]=system_aqua1.makeOrGetProblemForBoundary(
            #    xB,U,xC0,method="Newton-CG",folder='data7')
            
            #data[h]=system_aqua1.makeOrGetProblemForBoundary(
            #    xB,U,xC0,method="L-BFGS-B", folder='data8')
            
            #data[h]=system_aqua1.makeOrGetProblemForBoundary(
            #    xB,U,xC0,method="TNC", folder='data9')
    
            # We should try this with and without constraints
            #data[h]=system_aqua1.makeOrGetProblemForBoundary(
            #    xB,U,xC0,method="SLSQP", folder='data10')
    
            # This requires the Jacobian, so does not work for me
            #data[h]=system_aqua1.makeOrGetProblemForBoundary(
            #    xB,U,xC0,method="dogleg", folder='data11')
            
            # This requires the Jacobian, so does not work for me
            #data[h]=system_aqua1.makeOrGetProblemForBoundary(
            #    xB,U,xC0,method="trust-ncg", folder='data12')
            
            data[h] = system_aqua1.makeOrGetProblemForBoundary(
                xB, U, xC0, method=method, folder=folder)
            
    
        if False:
            t = timeit.timeit('func()', number=1, globals=globals())
            xB,bdry,p,opt,err = data[h]
            print('Done in ', t, 's')
        else:
            func()
            xB,bdry,p,opt,err = data[h]
        
        # Construct structured arrays from the data
        pi = np.array(p.input)
        d=np.dtype([('Q','f'),('cons','(11,)f')])
        po=np.empty(len(p.input),dtype=d)
        po.fill(np.nan)
        for j,ix in enumerate(p.input):
            try:
                po[j] = p.lookup(ix)
            except Exception as e:
                break
        
            
        # Mask the output by validity and find the maximum
        mask = (po['cons'] >= 0).all(axis=1)
        masknan = np.zeros_like(po['Q'])
        masknan.fill(np.nan)
        masknan[mask] = 1
        try:
            Qmax = po['Q'][mask].max()
        except:
            Qmax = 0
        QQmax.append(Qmax)
    
        
        if True:
        
            # Plot masked objective history
            plt.figure()
            plt.plot(po['Q'] * masknan)
            plt.xlabel('iteration')
            plt.ylabel('Cooling capacity')
            plt.title('Objective max = {}'.format(Qmax))
        
        
        if False:
            # Plot constraints
            fig,axs = plt.subplots(11,1,figsize=(6,20))
            axs[0].set_title('Constraints')
            for j,ax in enumerate(axs):
                ax.plot(po['cons'][:,j],label=str(j))
                ax.legend()
            ax.set_xlabel('iteration')
        
        #xopt = array([   0.50686621,  277.35844879,  313.01669312,  313.78839   ,
        #    309.54452247,  376.29152241])
        
        if False:
            # Flat plot of input variables
            plt.figure()
            for j in range(6):
                plt.plot(pi[:,j]/pi[0,j],label=str(j))
                
            plt.legend(loc='best')
            plt.grid()
            plt.xlabel('iteration')
            plt.ylabel('Normalized input variable')
        
    plt.figure()
    plt.plot(rejectTT,QQmax)
    plt.xlabel('Rejectiont temperature')
    plt.ylabel('Optimal cooling capacity')
    plt.show()

for method,folder in [('COBYLA','data'),
                      (None,'data2'),
                      ('Nelder-Mead','data3'),
                      ('Powell','data4'),
                      ('CG','data5'),
                      ('BFGS','data6'),
                      ('Newton-CG','data7'),
                      ('L-BFGS-B','data8'),
                      ('TNC','data9'),
                      ('SLSQP','data10'),
                      ('dogleg','data11'),
                      ('trust-ncg','data12')]:
        pass
for method,folder in [('COBYLA','data'),]:
    playback(method,folder)
