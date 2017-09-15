# -*- coding: utf-8 -*-
"""
Created on Sun Jul 10 23:39:27 2016

@author: nfette
"""

import numpy as np
import tabulate
from HRHX_integral_model import plotFlow as plotFlow, \
    counterflow_integrator as hxci, \
    streamExample1 as se1
import ammonia1
from scipy.optimize import minimize, basinhopping, brute
from scipy.special import expit
import hashlib
import pickle


def makeBoundary(x):
    """Given a iterable x that represents 5 stream inlets as T,m pairs,
    returns the Boundary object representing that.
    
    t0,m0,t1,m1,t2,m2,t3,m3,t4,m4 = x
    
    Note that boundary streams are
    
    0
        heat
    1
        absorberReject
    2
        condenserReject
    3
        cold
    4
        rectifierReject
    """
    t0, m0, t1, m1, t2, m2, t3, m3, t4, m4 = x
    # Units: K, kg/s, kW/kg-K
    cp = 4.179
    return Boundary(
        heat=se1(t0, m0, cp),
        absorberReject=se1(t1, m1, cp),
        condReject=se1(t2, m2, cp),
        cold=se1(t3, m3, cp),
        rectifierReject=se1(t4, m4, cp))


def makeChiller(x):
    """Given a iterable x that represents a chiller state,
    return a chiller object.
    
    Args
    ====
    
    x[0]
        m_rich (kg/s)
    x[1]
        T_evap (K)
    x[2]
        T_cond (K)
    x[3]
        T_rect (K)
    x[4]
        T_abs_outlet (K)
    x[5]    
        T_gen_outlet (K)
    """
    chiller = ammonia1.AmmoniaChiller()
    m_rich, T_evap, T_cond, T_rect, T_abs_outlet, T_gen_outlet = x
    chiller.update(m_rich=m_rich,
                   T_evap=T_evap,
                   T_cond=T_cond,
                   T_rect=T_rect,
                   T_abs_outlet=T_abs_outlet,
                   T_gen_outlet=T_gen_outlet)
    return chiller


class Boundary(object):
    def __init__(self, heat, absorberReject, condReject,
                 cold, rectifierReject):
        self.heat = heat
        self.absorberReject = absorberReject
        self.condReject = condReject
        self.cold = cold
        self.rectifierReject = rectifierReject

    def __repr__(self):
        import tabulate
        result = []
        for name in "heat absorberReject condReject cold rectifierReject".split():
            stream = self.__getattribute__(name)
            result.append((name, stream.T_inlet, stream.mdot, stream.cp))
        return tabulate.tabulate(result, "stream T_inlet mdot cp".split())

    def _repr_html_(self):
        import tabulate
        result = []
        for name in "heat absorberReject condReject cold rectifierReject".split():
            stream = self.__getattribute__(name)
            result.append((name, stream.T_inlet, stream.mdot, stream.cp))
        return tabulate.tabulate(result,
                                 "stream T_inlet mdot cp".split(),
                                 tablefmt='html')


class System(object):
    def __init__(self, boundary, chiller):
        self.boundary = boundary
        self.chiller = chiller

        self.hxs = {}
        self.hxs['gen'] = hxci(chiller.getGeneratorStream(), boundary.heat)
        self.hxs['rect'] = hxci(boundary.rectifierReject, chiller.getRectifierStream())
        self.hxs['abs'] = hxci(boundary.absorberReject, chiller.getAbsorberStream())
        self.hxs['cond'] = hxci(boundary.condReject, chiller.getCondenserStream())
        self.hxs['evap'] = hxci(chiller.getEvaporatorStream(), boundary.cold)
        self.Q = {'gen': chiller.Q_gen,
                  'rect': -chiller.Q_reflux,
                  'abs': -chiller.Q_abs,
                  'cond': -chiller.Q_cond,
                  'evap': chiller.Q_evap}

        self.totalUA = 0
        self.data = []
        for name in self.hxs:
            self.hxs[name].calcQmax()
            deltaT, epsilon, UA = np.inf, np.inf, np.inf
            try:
                deltaT, epsilon, UA = self.hxs[name].calcUA2(self.Q[name])
            except ValueError as e:
                pass
            self.data.append((name, deltaT, epsilon, UA, self.Q[name]))
            self.totalUA += UA

    def __repr__(self):
        result = "{}\ntotalUA = {}".format(
            tabulate.tabulate(self.data, "name deltaT epsilon UA Q".split()),
            self.totalUA)
        return result

    def _repr_html_(self):
        total_line = [["total", 0, 0, self.totalUA, 0]]
        table = tabulate.tabulate(self.data + total_line,
                                  "name deltaT epsilon UA Q".split(),
                                  tablefmt='html')
        return """{}""".format(table)

    def display(self):
        import matplotlib.pyplot as plt
        for name in self.hxs:
            plotFlow(self.hxs[name], None, self.Q[name])
            plt.title(name)
        """
        Q = [sys.chiller.Q_gen,sys.chiller.Q_abs,sys.chiller.Q_cond,
             sys.chiller.Q_evap,sys.chiller.Q_reflux]
        UA = [sys.genUA,sys.absUA,sys.condUA,sys.evapUA,sys.rectUA]
        component = "Generator Absorber Condenser Evaporator Rectifier".split()
        width = 1

        plt.figure()
        bar1=plt.bar(range(5),Q,width,tick_label=component)
        plt.figure()
        bar2=plt.bar(range(5),UA,width,tick_label=component)
        """

class Problem(object):
    r"""
    A wrapper for aqua-water absorption system plus other features.
    
    * Constraint evaluation
    * Storage of previously evaluated systems for faster lookup of repeats
    
    If ``hardConstraints=True``, then constraints will be applied to modify the
    value returned by ``objective()`` using this formula (no reference, just
    something I invented)::
        
        for each constraint:
            objective *= expit(constraint/mu)
    
    where each :math:`constraint \ge 0` in the feasible region following the
    convention of ``scipy.optimize.minimize`` constraints. This mechanism is
    not great because on the boundary, the modifier is 0.5, so the constraint
    can be violated.
    
    Smarter constraint mechanisms are
    suggested by Michael Wetter in the GenOpt manual, called boundary and
    penalty functions. In this case each :math:`g(x) \le 0`, and 
    
    .. math::
        
        \tilde{f}(x,\mu) &=& f(x) + \mu\frac{1}{\sum g^i(x)}
        
        \tilde{f}(x,\mu) &=& f(x) + \mu\sum{\max(0,g^i(x))^2}
    
    The boundary functions explode at the boundary to encourage staying within
    the interior of the feasible region. The penalty functions allow to escape
    the feasible region, which accomodates equality constraints.
        
    Inputs
    ======
    bdry : Boundary
        The object listing external streams that are the boundary conditions
    UAgoal : float
        The UAgoal
    constraintMode : string, optional
        How to handle constraints.
        
        None (default):
            Defaults = None (then constraints are available to scipy.optimize as a
            list).
        'expit':
            multiply the objective by a scalar that tapers from 1
            when constraints are satisfied to 0 when they are violated.
        'genopt1':
            Use boundary functions for hard constraints (valid cycle and each
            UA > 0), and penalty functions for soft constraints (total UA goal).
            No scaling is applied.
        
    mu : float, optional
        Scale over which constraints send the objective to zero, eg.
        ``objective *= expit(constraint/mu)``. Default = 0.1.
        
    Attributes
    ==========
    Ncons : int
        The number of constraints
    constraints : list
        Each item is a dictionary having 'type','fun','args' keys.
        The list is empty if Problem was initialized with hardConstraints=True.
    input : list
        Each item is the vector received from minimize, in order.
    output : dict
        Each keys is the hash of an input vector; values are (Q,cons) where
        Q is the cooling capacity and cons are the raw constraint values.
    """

    def __init__(self, bdry, UAgoal, constraintMode=None, mu=0.1):
        self.bdry = bdry
        self.UAgoal = UAgoal
        self.constraintMode = constraintMode
        self.mu = 0.1
        self.input = []
        self.output = dict()
        self.Ncons = 11
        if constraintMode == None:
            # Soft constraints mode: this is sent to minimizer
            self.constraints = [{'type': 'ineq',
                                 'fun': self.constraint,
                                 'args': (i,)} for i in range(self.Ncons)]
        else:
            # Hard constraint mode: objective function modified
            self.constraints = []



    def objective(self, x, stepNumber=1):
        """Return the value of the objective function to minimize. This is
        :math:`-Q_{evap}` modified depending on constraintMode.
        """
        Q, cons = self.lookup(x, printing=True)
        if self.constraintMode == 'expit':
            mu = self.mu
            for c in cons:
                Q *= expit(c / mu)
            return -Q
        elif self.constraintMode == 'genopt1':
            # Each cons[i] > 0 by scipy convention,
            # So need to change the sign of the second term.

            # GenOpt eqn. 8.6
            barrier = 0
            for c in cons[:10]:
                g = -c
                barrier += g
            # Eqn. 8.7
            mu1 = stepNumber ** (-2)

            # GenOpt eqn. 8.8
            penalty = 0
            for c in cons[10:]:
                g = -c
                penalty += max(0, g) ** 2
            # Eqn 8.9
            mu2 = stepNumber ** (2)

            f_tilde = -Q + mu1 * barrier + mu2 * penalty
            return f_tilde
        else:
            return -Q

    def constraint(self, x, *args):
        i, = args
        Q, cons = self.lookup(x)
        return cons[i]

    def lookup(self, x, printing=False):
        # print("Looking up {}".format(x))
        h = hasher(x)
        if h in self.output:
            return self.output[h]
        else:
            self.input.append(x.copy())
            cons = [x[0],
                    x[2] - x[1],
                    x[3] - x[2],
                    x[5] - x[3],
                    x[5] - x[4]]
            try:
                sys = System(self.bdry, makeChiller(x))
                Q = sys.chiller.Q_evap
                for name, deltaT, epsilon, UA, Qhx in sys.data:
                    cons.append(deltaT)
                cons.append(self.UAgoal - sys.totalUA)
            except Exception as e:
                print("At {} caught {}".format(h, e))
                Q = 0
                while len(cons) < self.Ncons:
                    cons.append(-1)
            self.output[h] = (Q, cons)
            if printing:
                print("@ hash {:15}, maxcv = {:12.5f}, Q = {:12.5f}".format(h, -min(cons), Q))

            return Q, cons


def systemConstraints(sys, UAgoal):
    T_evap = sys.chiller.refrig_evap_outlet.T
    T_cond = sys.chiller.refrig_cond_outlet.T
    T_rect = sys.chiller.refrig_rect_outlet.T
    T_gen_outlet = sys.chiller.weak_gen_outlet.T
    T_abs_outlet = sys.chiller.rich_abs_outlet.T
    hardConstraints = [T_cond - T_evap,
                       T_rect - T_cond,
                       T_gen_outlet - T_rect,
                       T_gen_outlet - T_abs_outlet]
    # name, deltaT, epsilon, UA, Q
    deltaTT = [hx[1] for hx in sys.data]
    for deltaT in deltaTT:
        hardConstraints.append(deltaT)

    softConstraints = [UAgoal - sys.totalUA]

    # GenOpt eqn. 8.6
    barrier = 0
    for c in hardConstraints:
        g = -c
        barrier += g
    # Eqn. 8.7
    # mu1 = stepNumber ** (-2)

    # GenOpt eqn. 8.8
    penalty = 0
    for c in softConstraints:
        g = -c
        penalty += max(0, g) ** 2
    # Eqn 8.9
    # mu2 = stepNumber ** (2)

    return barrier, penalty


def main():
    # Boundary
    xB0 = [400, 1,
           305, 3,
           305, 5,
           285, 4,
           305, 0.15]
    bdry = makeBoundary(xB0)
    print(bdry)
    # Chiller
    xC0 = np.array((0.40, 278.45, 311.85, 313.65, 310.15, 374.15))
    xC0 = np.array([0.51284472, 277.97717012, 312.16427764, 313.6952877,
                    310.24856734, 374.14020482])

    ch = makeChiller(xC0)
    print(ch)
    # System
    sys = System(bdry, ch)
    print(sys)
    # sys.display()
    return bdry, xC0, ch, sys


def hasher(x):
    if not isinstance(x, np.ndarray):
        x = np.array(x, dtype=np.double)
    # The built-in hash is randomly seeded = bad!
    # It returns an int.
    # x.flags.writeable = False
    # h = hash(y.data.tobytes())
    # x.flags.writeable = True
    # The shortest digest in hashlib is md5 at 16 bytes. Good enough.
    # Can output as hex string (32 chars).
    h = hashlib.md5(x.astype(np.double)).hexdigest()
    return h


def makeOrGetProblemForBoundary(xB, U, xC, method=None, options=None, folder='data', create=False):
    """Wrap the minimizer and problem with disk storage.
    If the given boundary constraint has been tried before, the data file
    in ``folder`` will be loaded to save time.
    
    Inputs
    ======
    xB : array
        Input to makeBoundary
    UA : float
        Goal for the UA value
    xC : array
        Initial input to makeChiller (ignored when loading from file)
    method : string, optional
        Passed directly to minimize
    options : dict, optional
        Passed directly to minimize
    folder : string, optional
        Path to the data folder relative to ../ (defaults to 'data')
    create : bool, optional
        Whether optimization should proceed if the problem is not found on disk
        (defaults to False).
    
    Returns
    =======
    xB : array
        Echo the boundary input back at you!
    bdry : Boundary
        The actual boundary object (so many streams)
    p : Problem
        Object with reference to boundary, list of input points,
        and dictionary of outputs
    opt : OptimizeResult or None
        The result if one was reached, or None if an error occured
    err : Exception
        The exception, if one was caught during optimization. Errors are caught
        so that data may be saved.
    """
    h = hasher(xB)
    fname = '../{}/system_aqua1_{}_{}.pkl'.format(folder, U, h)
    try:
        with open(fname, 'rb') as f:
            data = pickle.load(f)
        print("Loaded the case from ", fname)
        xB, bdry, p, opt, err = data
    except FileNotFoundError:
        print("Running the case for ", fname)
        # Create the data
        opt, err = None, None
        bdry = makeBoundary(xB)
        p = Problem(bdry, U, constraintMode='expit')

        if create:
            try:
                opt = minimize(p.objective,
                               xC,
                               constraints=p.constraints,
                               method=method,
                               options=options)
                print(opt)
            except:
                import sys
                c, err, tb = sys.exc_info()
                print(err)

            print("Saving to ", fname)
            # Now store the data
            data = xB, bdry, p, opt, err
            with open(fname, 'wb') as f:
                pickle.dump(data, f)

    return xB, bdry, p, opt, err


def main1():
    """Script to try out your options for optimization"""

    # bdry, xC0, ch, sys = main()
    U = 100
    xB = [400, 1, 305, 3, 305, 5, 285, 4, 305, 0.15]
    xC = np.array([0.51284472, 277.97717012, 312.16427764, 313.6952877,
                   310.24856734, 374.14020482])
    # Moved to aqua_case_studies.py
    # xB,bdry,p,opt,err = makeOrGetProblemForBoundary(xB,U,xC)


    # opt = scipy.optimize.minimize(p.objective, xC0, constraints=p.constraints,
    #    method="COBYLA",
    #    options={"rhobeg":0.05})
    # "maxiter":100
    # "catol":0.1

    bdry = makeBoundary(xB)
    # p = Problem(bdry, U, hardConstraints=True)
    p = Problem(bdry, U, constraintMode='expit')
    opt = basinhopping(p.objective, xC, niter=3,
                       minimizer_kwargs=dict(options=dict(maxiter=2)))
    # tol = 1
    # for i in range(5):
    #    tol *= 0.25
    #    hi,low = 1 - tol, 1 + tol
    #    ranges = [(xi * low, xi * hi) for xi in xC]
    #    xC = brute(p.objective,ranges,Ns=3)
    #    print(xC)

    print(opt)
    xC1 = opt.x
    ch1 = makeChiller(xC1)
    print(ch1)
    sys1 = System(bdry, ch1)
    print(sys1)
    sys1.display()

    d = np.dtype([('Q', 'f'), ('cons', '(11,)f')])
    po = np.empty(len(p.input), dtype=d)
    po.fill(np.nan)
    for i, ix in enumerate(p.input):
        po[i] = p.lookup(ix)
    import matplotlib.pyplot as plt
    plt.figure()
    plt.plot(po['Q'])
    plt.xlabel('iteration')
    plt.ylabel('Cooling capacity')
    plt.title('Objective')

    fig, axs = plt.subplots(11, 1, figsize=(6, 20))
    axs[0].set_title('Constraints')
    for i, ax in enumerate(axs):
        ax.plot(po['cons'][:, i], label=str(i))
        ax.legend()
    ax.set_xlabel('iteration')

    # xopt = array([   0.50686621,  277.35844879,  313.01669312,  313.78839   ,
    #    309.54452247,  376.29152241])

    pi = np.array(p.input)
    plt.figure()
    for i in range(6):
        plt.plot(pi[:, i] / pi[0, i], label=str(i))

    plt.legend(loc='best')
    plt.grid()
    plt.xlabel('iteration')
    plt.ylabel('Normalized input variable')


def main2():
    """Script to run with GenOpt, using a stupid barrier function"""

    import sys  # command args
    import traceback  # display errors into log
    import json  # read input file with standard, structured format
    from decimal import Decimal  # convert 'inf' to 'Infinity' etc

    infile, outfile, logfile = sys.argv[1:4]

    with open(logfile, 'w') as log:
        print('infile="{}"'.format(infile), file=log)
        print('outfile="{}"'.format(outfile), file=log)
        print('logfile="{}"'.format(logfile), file=log)

        try:
            with open(infile) as fin:
                data = json.load(fin)
            stepNumber = data['stepNumber']
            UAgoal = data['UAgoal']
            xB = data['xB']
            xC = data['xC']
            print('UAgoal={}'.format(UAgoal), file=log)
            print('xB={}'.format(xB), file=log)
            print('xC={}'.format(xC), file=log)
            bdry = makeBoundary(xB)
            chiller = makeChiller(xC)
            sys = System(bdry, chiller)
            print(sys, file=log)

            with open(outfile, 'w') as fout:
                print("Q_evap=", chiller.Q_evap, file=fout)
                sum_g_barrier = 0
                for row in sys.data:
                    name, deltaT, epsilon, UA, Qhx = row
                    print("UA_{}={}".format(name, Decimal(UA)), file=fout)
                    print("deltaT_{}={}".format(name, Decimal(deltaT)), file=fout)
                    sum_g_barrier += deltaT
                print("sum_g_barrier={}".format(Decimal(sum_g_barrier)), file=fout)
                mu = 3 / stepNumber
                print("mu={}".format(mu), file=fout)

        except:
            traceback.print_exc(file=log)


def main3(infile, outfile, logfile):
    """Script to run with GenOpt, using more intelligent constraints."""

    import traceback  # display errors into log
    import json  # read input file with standard, structured format
    from decimal import Decimal  # convert 'inf' to 'Infinity' etc

    with open(logfile, 'w') as log:
        print('infile="{}"'.format(infile), file=log)
        print('outfile="{}"'.format(outfile), file=log)
        print('logfile="{}"'.format(logfile), file=log)

        try:
            with open(infile) as fin:
                data = json.load(fin)
            stepNumber = data['stepNumber']
            UAgoal = data['UAgoal']
            xB = data['xB']
            xC = data['xC']
            print('UAgoal={}'.format(UAgoal), file=log)
            print('xB={}'.format(xB), file=log)
            print('xC={}'.format(xC), file=log)
            bdry = makeBoundary(xB)
            chiller = makeChiller(xC)
            sys = System(bdry, chiller)
            print(sys, file=log)
            barrier, penalty = systemConstraints(sys, UAgoal)

            with open(outfile, 'w') as fout:
                if np.isinf(barrier) or np.isinf(penalty):
                    print("Q_evap=", -1, file=fout)
                    print("barrier=", 0, file=fout)
                    print("penalty=", 0, file=fout)
                else:
                    print("Q_evap=", chiller.Q_evap, file=fout)
                    print("barrier=", barrier, file=fout)
                    print("penalty=", penalty, file=fout)

        except:
            traceback.print_exc(file=log)


if __name__ == "__main__":
    import sys  # command args
    main3(*sys.argv)
