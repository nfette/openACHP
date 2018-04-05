"""
A useful compendium of executable wisdom.
"""

import numpy
import pandas
import CoolProp.CoolProp as CP
from ammonia_props import massFractionToMolar, AmmoniaProps
import ammonia1

amm=AmmoniaProps()

# Absorber: plotting T_abs,max vs x_rich to find maximum
def P_abs_water(x_rich,x_refrig,P_evap):
    x_molar_rich = massFractionToMolar(x_rich)
    x_molar_refrig = massFractionToMolar(x_refrig)
    alpha = (1 - x_molar_refrig) / (1 - x_molar_rich)
    result = alpha * P_evap
    return result

def P_abs_ammonia(x_rich,x_refrig,P_evap):
    x_molar_rich = massFractionToMolar(x_rich)
    x_molar_refrig = massFractionToMolar(x_refrig)
    #alpha = (1 - x_molar_refrig) / (1 - x_molar_rich)
    alpha = (x_molar_refrig) / (x_molar_rich)
    result = alpha * P_evap
    return result

def T_abs_max_func(x_rich,P_abs):
    return amm.props2(x=x_rich,P=P_abs,Qu=0,out='T')

n_constraints_aqua_ammonia = 10
spec_columns = ['T_abs','T_cond','T_evap','T_gen','T_rect','x_refrig']
constraint_columns = range(n_constraints_aqua_ammonia)
calc_columns = ['x_refrig_min','x_refrig_max','T_cond_max','T_evap_min',
      'P_cond','P_evap','T_abs_max1','x_rich_max','x_rich',
      'T_abs_max2','T_gen_min']
all_columns = pandas.MultiIndex.from_product([['spec'],spec_columns]) \
    .append(pandas.MultiIndex.from_product([['constraint'],constraint_columns])) \
    .append(pandas.MultiIndex.from_product([['calc'],calc_columns])) \
    .append(pandas.MultiIndex.from_product([['mapped'],spec_columns]))


def makeSpec(T_abs=300,T_cond=300,T_evap=278,T_gen=380,T_rect=310,x_refrig=0.9998):
    return pandas.Series(data=[T_abs, T_cond, T_evap, T_gen, T_rect, x_refrig],
                         index=spec_columns)


def makeChiller(spec):
    ch = ammonia1.AmmoniaChiller()
    ch.update(x_refrig=spec.x_refrig, T_evap=spec.T_evap, T_cond=spec.T_cond,
          T_abs_outlet=spec.T_abs, T_gen_outlet=spec.T_gen, T_rect=spec.T_rect)
    return ch


class AquaChillerSpec1:
    
    n_cons = n_constraints_aqua_ammonia
    Qu_evap = 0.998
    x_refrig_min = 0.9
    x_refrig_max = 0.999999
    x_rich_min = 0.2
    x_rich_max_max = 1.0
    P_cond_max = 80.85
    
    def __init__(self,spec,throws=True):
        self.spec = spec
        self.html_opts = dict(inputs=True,limits=True,mapped=True,constraints=True,warnings=True)
        # Constraints
        C = pandas.Series(index=constraint_columns)
        calc = pandas.Series(index=calc_columns)
        mapped = spec.copy()
        messages=[]
        
        self.C = C
        self.calc = calc
        self.mapped = mapped
        self.messages = messages
        
        try:
            # Todo: implement the constraints on T_rect.
            #  eg C[9] to enforce T_cond < T_rect,
            #  C[10] to enforce T_rect <= T(P=Pcond,x=xrefrig,Qu=1).
            # Todo: Do we need a lower bound for T_abs?
            
            # A trivial constraint that is always satisfied
            C[0] = 0
            
            # Basics
            C[1] = spec.x_refrig - self.x_refrig_min
            if C[1] < 0:
                messages.append("x_refrig ({:g}) should be greater than {:g} but is not."
                                     .format(spec.x_refrig, self.x_refrig_min))
                mapped.x_refrig = self.x_refrig_min
            
            C[2] = self.x_refrig_max - spec.x_refrig
            if C[2] < 0:
                messages.append("x_refrig ({:g}) should be less than {:g} but is not."
                                    .format(spec.x_refrig, self.x_refrig_max))
                mapped.x_refrig = self.x_refrig_max
            
            calc.T_evap_min = amm.props2(P=0.01, x=mapped.x_refrig, Qu=self.Qu_evap, out="T")
            C[3] = (spec.T_evap - calc.T_evap_min)
            if C[3] < 0:
                messages.append("T_evap ({:g}) should be greater than {:g} but is not."
                                     .format(self.T_evap, calc.T_evap_min))
                mapped.T_evap = calc.T_evap_min
            
            # convert mass fraction to molar only for Refprop
            x_molar_refrig = massFractionToMolar(mapped.x_refrig)
            T_lookup_max_bounds = []
            T_lookup_max_bounds.append(CP.PropsSI('Tcrit',
                'REFPROP::ammonia[{}]&water[{}]'.format(
                x_molar_refrig, 1.0 - x_molar_refrig)))
            #T_lookup_max_bounds.append(550)
            T_lookup_max_bounds.append(amm.props2(P=80, x=mapped.x_refrig, Qu=0, out='T'))
            calc.T_cond_max = numpy.min(T_lookup_max_bounds)
            C[4] = (calc.T_cond_max - spec.T_cond)
            if C[4] < 0:
                messages.append("T_cond ({:g}) should be less than {:g} but is not."
                                     .format(spec.T_cond, calc.T_cond_max))
                mapped.T_cond = calc.T_cond_max
                
            C[5] = (spec.T_cond - spec.T_evap)
            if C[5] < 0:
                messages.append("T_cond ({:g}) should be greater than T_evap ({:g}) but is not."
                                    .format(spec.T_cond, spec.T_evap))
                if mapped.T_evap < calc.T_cond_max:
                    mapped.T_cond = mapped.T_evap
                elif mapped.T_cond > calc.T_evap_min:
                    mapped.T_evap = mapped.T_cond
                else:
                    mapped.T_evap = calc.T_cond_max
                    mapped.T_cond = mapped.T_evap

            calc.P_cond = amm.props2(T=mapped.T_cond, x=mapped.x_refrig, Qu=0, out="P")
            calc.P_evap = amm.props2(T=mapped.T_evap, x=mapped.x_refrig, Qu=self.Qu_evap, out="P")

            # TODO: explain in write-up above
            # I had started off with constraint 'x_rich_max > 0',
            # but it appears that an equivalent statement is 'T_abs < T(x=0, ...)'.
            # This seems preferrable since following calculations depend on T_abs,
            # so we can adjust T_abs to the allowed limit,
            # whereas violation of limit on x_rich_max is not as easy to correct/override(?)
            calc.T_abs_max1 = amm.props2(x=self.x_rich_min, P=calc.P_evap, Qu=0, out="T")
            C[6] = (calc.T_abs_max1 - spec.T_abs)
            if C[6] < 0:
                messages.append("T_abs ({:}) should be less than {:g} but is not."
                                     .format(spec.T_abs, calc.T_abs_max1))
                mapped.T_abs = calc.T_abs_max1
            
            # Also, enforce 'x_rich_max < 1'
            # First determine maximum equilibrium pressure at this temperature,
            # that is, when x_rich is 1, so that we don't call State function with invalid inputs.
            # If evaporator pressure is higher, no problem; even pure ammonia can absorb
            # by simple condensation.
            P_max = amm.props2(T=mapped.T_abs, x=self.x_rich_max_max, Qu=0, out="P")
            if P_max < calc.P_evap:
                calc.x_rich_max = self.x_rich_max_max
                calc.x_rich = calc.x_rich_max
            else:
                # This must be guarded by the above if-statement.
                # For P_evap slightly in excess of max, the function
                # returns x slightly greater than 1. For larger values, it fails altogether.
                calc.x_rich_max = amm.props2(T=mapped.T_abs, P=calc.P_evap, Qu=0, out="x")
                calc.x_rich = calc.x_rich_max
                if calc.x_rich > self.x_rich_max_max:
                    calc.x_rich = self.x_rich_max_max
            
            calc.T_abs_max2 = T_abs_max_func(calc.x_rich,
                                            P_abs_ammonia(calc.x_rich, mapped.x_refrig, calc.P_evap))
            C[7] = (calc.T_abs_max2 - spec.T_abs)
            if C[7] < 0:
                messages.append("T_abs ({:g}) should be less than {:g} but is not.".format(
                    spec.T_abs, calc.T_abs_max2))
                # Should already have been fixed by previous constraint ... get rid of this one
                # self.T_abs_mapped = self.T_abs_max2
            
            # TODO: fix this ... should map T_cond sooner, above.
            # Sometimes, the given pressure exceeds what the function can handle (see benchmarks).
            # A guaranteed level of pressure is substantially lower, about 80 bar.
            # But to use that pressure, we should also need a lower T_cond.
            try:
                calc.T_gen_min = amm.props2(P=calc.P_cond, x=calc.x_rich, Qu=0).T
            except KeyError as e:
                calc.T_gen_min = amm.props2(P=self.P_cond_max, x=calc.x_rich, Qu=0).T
            C[8] = (spec.T_gen - calc.T_gen_min)
            if C[8] < 0:
                messages.append("T_gen ({:g}) should be greater than {:g} but is not.".format(
                    spec.T_gen, calc.T_gen_min))
                mapped.T_gen = calc.T_gen_min

            # TODO: reorder newly added constraint.
            # Enforces T_rect - T_cond > 0.
            C[9] = (spec.T_rect - spec.T_cond)
            if C[9] < 0:
                messages.append("T_rect ({:g}) should be greater than T_cond ({:g}) but is not.".format(
                    spec.T_rect, spec.T_cond))
                # Todo: decide what to do next.
                # Corrective Action, Case 1
                # mapped.T_rect = spec.T_cond + 1
                # Now we would have to go back and update variables depending on T_rect.
                # Corrective Action, Case 2
                # mapped.T_cond = spec.T_rect - 1
                # Now we would have to go back and update variables depending on T_cond.

        except KeyError as e:
            messages.append(e.__repr__())
            if throws:
                raise
        except ValueError as e:
            messages.append(e.__repr__())
            if throws:
                raise
        except:
            raise
    
    def _repr_html_(self):
        return self.html_helper(**self.html_opts)
    
    def html_helper(self,inputs=True,limits=True,mapped=True,constraints=True,warnings=True):
        result = ""
        if inputs:
            result += "<h5>Inputs</h5><pre>{}</pre>".format(self.spec)
        if limits:
            result += "<h5>Computed limits</h5><pre>{}</pre>".format(self.calc)
        if mapped:
            result += "<h5>Mapped inputs</h5><pre>{}</pre>".format(self.mapped)
        if constraints:
            result += "<h5>Constraint functions</h5><pre>{}</pre>".format(self.C)
        if warnings:
            result += """<h5>Warnings</h5>
                <ul>{}</ul>""".format("".join(
                map(lambda msg: "<li>{}</li>".format(msg),
                    ["None."] if len(self.messages) == 0 else self.messages)))
        return result
