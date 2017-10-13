
import numpy
import matplotlib.pyplot as plt
import ammonia_props
amm = ammonia_props.AmmoniaProps()

P0=1.0
x_range = numpy.linspace(0,1,101)
p_min = numpy.empty_like(x_range)
p_min.fill(numpy.nan)
t_out = p_min.copy()
for i, x in enumerate(x_range):
    try:
        P=P0
        for j in range(20):
            state=amm.props2(x=x,P=P,Qu=0)
            P *= 0.8
    except:
        pass
    p_min[i] = P
    t_out[i] = state.T

plt.figure()
plt.plot(x_range,p_min)
plt.figure()
plt.plot(x_range,t_out)

plt.show()
