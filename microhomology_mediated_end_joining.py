from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt

Model()

# STEPS:
# 1. Parp1 binds DSB
# 2. CtIP binds Parp1
# 3. MRE11 binds Parp1


Monomer("DSB", ["b"])
Monomer("Parp1", ["dsb", "ctip_mre11"])
Monomer("CtIP", ["parp1"])
Monomer("MRE11", ["parp1"])

Parameter("DSB_0", 100)
Parameter("Parp1_0", 100)
Parameter("CtIP_0", 100)
Parameter("MRE11_0", 100)

Initial(DSB(b=None), DSB_0)
Initial(Parp1(dsb=None, ctip_mre11=None), Parp1_0)
Initial(CtIP(parp1=None), CtIP_0)
Initial(MRE11(parp1=None), MRE11_0)

Parameter("kf_DSB_Parp1", 1)
Parameter("kr_DSB_Parp1", 10)
Parameter("kf_DSB_CtIP", 1)
Parameter("kr_DSB_CtIP", 10)


Rule("Parp1_binds_DSB",
     DSB(b=None) + Parp1(dsb=None, ctip_mre11=None) |
     DSB(b=1) % Parp1(dsb=1, ctip_mre11=None), kf_DSB_Parp1, kr_DSB_Parp1)
Rule("CtIP_binds_Parp1",
     CtIP(parp1=None) + DSB(b=1) % Parp1(dsb=1, ctip_mre11=None) |
     CtIP(parp1=2) % DSB(b=1) % Parp1(dsb=1, ctip_mre11=2), kf_DSB_CtIP, kr_DSB_CtIP)
Rule("MRE11_binds_Parp1",
     MRE11(parp1=None) + Parp1(dsb=None) |
     MRE11(parp1=1) % Parp1(dsb=1), kf_DSB_Parp1, kr_DSB_Parp1)

print(model)
print(model.monomers)
print(model.parameters)
print(model.rules)
