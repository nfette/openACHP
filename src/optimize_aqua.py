import numpy
import ammonia1
import system_aqua1
import scipy.optimize
import aqua_chiller_spec1
import pandas


def saturate(x, bottom=-numpy.inf, top=0):
    a_bottom = numpy.empty_like(x)
    a_top = numpy.empty_like(x)
    a_bottom.fill(bottom)
    a_top.fill(top)
    return numpy.minimum(a_top,
                         numpy.maximum(a_bottom,
                                       x))


def barrier1(c, length_scale):
    """The domain for B is the feasible set only.
    We want B large near boundary, but small elsewhere.
    Feasible means c > 0 and g < 0. Please input c.
    """
    g = numpy.negative(c)
    b = numpy.exp(g / length_scale)
    return numpy.sum(b)


def decay1(step_number, initial_value = 1., rate = 1.):
    """A decaying function to scale barrier functions.
    step_number: as this increases, I decay.
    initial_value: value when step_number = 0.
    decay_rate: how fast to decay.
    """
    # Magnitude tends to zero, slowly
    mu_B = initial_value * numpy.exp(-rate * step_number)
    return mu_B


def penalty1(c, length_scale):
    """We want P = 0 for feasible, P > 0 elsewhere
    Feasible means c > 0 and g < 0.
    """
    g = numpy.negative(c)
    g = saturate(g,bottom=0,top=numpy.inf)
    p = (g / length_scale)**2
    return numpy.sum(p)


def grow1(step_number, initial_value = 1., rate = 1.):
    """A growing function to scale penalty functions."""
    # Magnitude tends to infinite, bit more quickly
    mu_P = initial_value * numpy.exp(rate * step_number)
    return mu_P


class Problem_A:
    def __init__(self, bdry, hx_delta_t_required=1.0, mu=0.1):
        self.bdry = bdry
        self.hx_delta_t_required = hx_delta_t_required
        self.mu = mu
        self.Ncons = 14
        self.n_calls = 0
        # Soft constraints mode: this is sent to minimizer
        self.constraints_each = [{'type': 'ineq',
                             'fun': self.constraint,
                             'args': (i,)
                            } for i in range(self.Ncons)]
        self.constraints_all = {'type': 'ineq',
                                'fun': self.constraint,
                               }

    def objective(self, x):
        step_number = numpy.floor(self.n_calls / 7)
        self.n_calls += 1
        print("[Problem_A.objective] Input vector:")
        print(x)

        Q,B,P = 0.,0.,0.
        try:
            spec = aqua_chiller_spec1.makeSpec(*x)
            a = aqua_chiller_spec1.AquaChillerSpec1(spec,False)
            ch = aqua_chiller_spec1.makeChiller(a.mapped)
            sys = system_aqua1.System(self.bdry, ch)

            # Barriers
            # Magnitude tends to zero, slowly
            mu_B = 1000. * numpy.exp(-0.1 * step_number)
            length_scale_b = 1
            # Or ... magnitude fixed, but shape changes
            mu_B = 1000
            length_scale_b = 1 * numpy.exp(-0.1 * step_number)

            # These are zero at the boundary ...
            barriers = [ch.check_rectifier_delta_T]
            B = mu_B * barrier1(barriers,length_scale_b)

            # Penalties
            # Magnitude tends to infinite
            mu_P = 1 * numpy.exp(0.3 * step_number)
            penalties = (spec - a.mapped) ** 2
            P = mu_P * penalty1(penalties,1)

            Q = sys.chiller.Q_evap
            
            print("Chiller Q = ",Q)
                  
        except KeyboardInterrupt as e:
            raise e
        except:
            Q = numpy.inf

        print(self.n_calls, step_number, Q, B, P, "\n", flush=True)
        return -Q + B + P

    def constraint(self, x, *args):
        print("[Problem_A.constraint] Input vector:")
        print(x)
        spec = aqua_chiller_spec1.makeSpec(*x)
        a = aqua_chiller_spec1.AquaChillerSpec1(spec,False)
        print("Mapped chiller spec:")
        print(a.mapped)
        ch = aqua_chiller_spec1.makeChiller(a.mapped)
        sys = system_aqua1.System(self.bdry, ch)
        cons = pandas.concat([a.C, sys.df.deltaT - self.hx_delta_t_required])
        if len(args) > 0:
            i, = args
            return cons[i]
        else:
            return cons

    def callback(self, x):
        print("Did an iteration at ", x)


class Problem_B:
    def __init__(self, bdry, Q_goal, mu=0.1):
        self.bdry = bdry
        self.Q_goal = Q_goal
        self.mu = mu
        self.Ncons = 7
        self.n_calls = 0
        # Soft constraints mode: this is sent to minimizer
        self.constraints = [{'type': 'ineq',
                             'fun': self.constraint,
                             'args': (i,)
                            } for i in range(self.Ncons)]
            
    def objective_raw(self, xC):
        try:
            ch = system_aqua1.makeChiller(xC)
            sys = system_aqua1.System(self.bdry, ch)
            UA = sys.totalUA
        except:
            UA = numpy.nan
        return UA

        
    def objective(self, xC):
        step_number = numpy.floor(self.n_calls / 7)
        self.n_calls += 1
        #print(xC,flush=True)
        UA,B,P = 0.,0.,0.
        try:
            ch = system_aqua1.makeChiller(xC)
            sys = system_aqua1.System(self.bdry, ch)
            
            # Barriers
            # Magnitude tends to zero, slowly
            mu_B = 1000. * numpy.exp(-0.1 * step_number)
            length_scale_b = 1
            # Or ... magnitude fixed, but shape changes
            mu_B = 1000
            length_scale_b = 1 * numpy.exp(-0.1 * step_number)
            
            # These are zero at the boundary ...
            barriers = [ch.check_rectifier_delta_T] \
                       + [deltaT
                          for name, deltaT, epsilon, UA, Qhx in sys.data]
            B = mu_B * barrier1(barriers,length_scale_b)
            
            # Penalties
            # Magnitude tends to infinite
            mu_P = 1 * numpy.exp(0.3 * step_number)
            penalties = [ch.Q_evap - self.Q_goal]
            P = mu_P * penalty1(penalties,1)
            
            UA = sys.totalUA
        except KeyboardInterrupt as e:
            raise e
        except:
            UA = numpy.inf
        
        # print(self.n_calls, step_number, UA, B, P, "\n", flush=True)
        return UA + B + P
    
    def constraint(self, x, *args):
        cons = [x[0] - 0.1,
                1. - x[0],
                x[2] - x[1] - 1.0,
                x[3] - x[2] - 0.1,
                x[4] - x[1] - 10.0,
                x[5] - x[3] - 1.0,
                x[5] - x[4] - 1.0]
        if len(args) > 0:
            i, = args
            return cons[i]
        else:
            return cons

    def callback(self, x):
        print("Did an iteration at ", x)



