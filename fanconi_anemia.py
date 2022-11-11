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

Parameter("FANCB_0", 100)
Parameter("FANCL_0", 80)
Parameter("FAAP100_0", 120)

Initial(FANCB(fancl=None, faap100=None), FANCB_0)
Initial(FANCL(fancb=None, faap100=None), FANCL_0)
Initial(FAAP100(fancb=None, fancl=None), FAAP100_0)

# Formation of BL100 complex
Parameter('kf_BL', 1)
Parameter('kr_BL', 1)
Parameter('kf_B100', 1)
Parameter('kr_B100', 1)
Parameter('kf_L100', 1)
Parameter('kr_L100', 1)
Parameter('kf_BL_100', 1)
Parameter('kr_BL_100', 1)
Parameter('kf_B100_L', 1)
Parameter('kr_B100_L', 1)
Parameter('kf_L100_B', 1)
Parameter('kr_L100_B', 1)
Rule("FANCB_binds_FANCL", FANCB(fancl=None, faap100=None) + FANCL(fancb=None, faap100=None) |
     FANCB(fancl=1, faap100=None) % FANCL(fancb=1, faap100=None), kf_BL, kr_BL)
Rule("FANCB_FANCL_binds_FAAP100", FANCB(fancl=1, faap100=None) % FANCL(fancb=1, faap100=None) +
     FAAP100(fancb=None, fancl=None) | FANCB(fancl=1, faap100=2) % FANCL(fancb=1, faap100=3) %
     FAAP100(fancb=2, fancl=3), kf_BL_100, kr_BL_100)

Rule("FANCB_binds_FAAP100", FANCB(fancl=None, faap100=None) + FAAP100(fancb=None, fancl=None) |
     FANCB(fancl=None, faap100=1) % FAAP100(fancb=1, fancl=None), kf_B100, kr_B100)
Rule("FANCB_FAAP100_binds_FANCL", FANCB(fancl=None, faap100=1) % FAAP100(fancb=1, fancl=None) +
     FANCL(fancb=None, faap100=None) | FANCB(fancl=1, faap100=2) % FANCL(fancb=1, faap100=3) %
     FAAP100(fancb=2, fancl=3), kf_B100_L, kr_B100_L)

Rule("FANCL_binds_FAAP100", FANCL(fancb=None, faap100=None) + FAAP100(fancb=None, fancl=None) |
     FANCL(fancb=None, faap100=1) % FAAP100(fancb=None, fancl=1), kf_L100, kr_L100)
Rule("FANCL_FAAP100_binds_FANCB", FANCL(fancb=None, faap100=1) % FAAP100(fancb=None, fancl=1) +
     FANCB(fancl=None, faap100=None) | FANCB(fancl=1, faap100=2) % FANCL(fancb=1, faap100=3) %
     FAAP100(fancb=2, fancl=3), kf_L100_B, kr_L100_B)

Observable("FANCB_free", FANCB(fancl=None, faap100=None))
Observable("FANCL_free", FANCL(fancb=None, faap100=None))
Observable("FAAP100_free", FAAP100(fancb=None, fancl=None))
Observable("BL100_Complex", FANCB(fancl=1, faap100=2) % FANCL(fancb=1, faap100=3) %
     FAAP100(fancb=2, fancl=3))

# TODO: write rules for CEF complex

# simulation commands

tspan = np.linspace(0,1,101)
sim = ScipyOdeSimulator(model, tspan, verbose=True)
result = sim.run()

for obs in model.observables:
    plt.plot(tspan, result.observables[obs.name], lw=2, label=obs.name)

plt.legend(loc="best")

plt.show()
