from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt

Model()

Monomer("FANCA", ["fancg", "faap20"])
Monomer("FANCG", ["fanca"])
Monomer("FAAP20", ["fanca"])

Monomer("FANCB", ["fancl", "faap100"])
Monomer("FANCL", ["fancb", "faap100"])
Monomer("FAAP100", ["fancb", "fancl"])

Monomer("FANCF", ["fancc"])
Monomer("FANCC", ["fancf", "fance"])
Monomer("FANCE", ["fancc"])

Parameter("FANCA_0", 100)
Parameter("FANCG_0", 80)
Parameter("FAAP20_0", 120)

Initial(FANCA(fancg=None, faap20=None), FANCA_0)
Initial(FANCG(fanca=None), FANCG_0)
Initial(FAAP20(fanca=None), FAAP20_0)

# Formination of AG20 complex
Parameter('kf_AG', 1)
Parameter('kr_AG', 1)
Parameter('kf_A20', 1)
Parameter('kr_A20', 1)
Rule("FANCA_binds_FANCG", FANCA(fancg=None) + FANCG(fanca=None) | FANCA(fancg=1) % FANCG(fanca=1), kf_AG, kr_AG)
Rule("FANCA_binds_FAAP20", FANCA(faap20=None) + FAAP20(fanca=None) | FANCA(faap20=1) % FAAP20(fanca=1), kf_A20, kr_A20)

Observable("FANCA_free", FANCA(fancg=None, faap20=None))
Observable("FANCG_free", FANCG(fanca=None))
Observable("FAAP20_free", FAAP20(fanca=None))
Observable("AG20_complex", FANCA(fancg=1, faap20=2) % FANCG(fanca=1) % FAAP20(fanca=2))

# simulation commands

tspan = np.linspace(0,1,101)
sim = ScipyOdeSimulator(model, tspan, verbose=True)
result = sim.run()

for obs in model.observables:
    plt.plot(tspan, result.observables[obs.name], lw=2, label=obs.name)

plt.legend(loc="best")

plt.show()
