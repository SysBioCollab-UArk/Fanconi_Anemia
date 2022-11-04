from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt

Model()

Monomer("A", ["b"])
Monomer("B", ["a"])

Parameter("A_0", 100)
Parameter("B_0", 80)

Initial(A(b=None), A_0)
Initial(B(a=None), B_0)

Parameter("Kf_AB", 1)
Parameter("Kr_AB", 10)

Rule("A_binds_B", A(b=None) + B(a=None) | A(b=1) % B(a=1), Kf_AB, Kr_AB)

Observable("A_free", A(b=None))
Observable("B_free", B(a=None))
Observable("AB_complex", A(b=1) % B(a=1))

# simulation commands

tspan = np.linspace(0,1,101)
sim = ScipyOdeSimulator(model, tspan, verbose=True)
result = sim.run()

for obs in model.observables:
    plt.plot(tspan, result.observables[obs.name], lw=2, label=obs.name)

plt.legend(loc="best")

plt.show()