from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt

Model()

# STEPS:
# 1. Parp1 binds DSB
# 2. MRE11 binds Parp1
# 3. CtIP binds Parp1

# Monomer('A', ['a', 'a'], {'a': ['u', 'p']})
# Parameter('k1', 100)
# Parameter('A_0', 200)
# Rule('r1', None >> A(a=MultiState('u', 'p')), k1)
# Initial(A(a=MultiState(('u', 1), 'p')) %
#         A(a=MultiState(('u', 1), 'u')), A_0)

Monomer("DSB", ["b", "b"])
Monomer("Parp1", ["dsb", "ctip_mre11"])
Monomer("CtIP", ["parp1"])
Monomer("MRE11", ["parp1_rad50"])
# Monomer("MRN", ["dsb"])
Monomer("RAD50", ["mre11", "nbs1"])
Monomer("NBS1", ["rad50"])
Monomer("RPA", ["dsb"])

Parameter("DSB_0", 100)
Parameter("Parp1_0", 100)
Parameter("CtIP_0", 100)
Parameter("MRE11_0", 100)
# Parameter("MRN_0", 100)
Parameter("RAD50_0", 100)
Parameter("NBS1_0", 100)
Parameter("RPA_0", 100)

Initial(DSB(b=MultiState(None, None)), DSB_0)
Initial(Parp1(dsb=None, ctip_mre11=None), Parp1_0)
Initial(CtIP(parp1=None), CtIP_0)
Initial(MRE11(parp1_rad50=None), MRE11_0)
Initial(RAD50(mre11=None, nbs1=None), RAD50_0)
Initial(NBS1(rad50=None), NBS1_0)
Initial(RPA(dsb=None), RPA_0)


Parameter("kf_DSB_Parp1", 1)
Parameter("kr_DSB_Parp1", 10)
Rule("Parp1_binds_DSB",
     DSB(b=None) + Parp1(dsb=None, ctip_mre11=None) |
     DSB(b=1) % Parp1(dsb=1, ctip_mre11=None),
     kf_DSB_Parp1, kr_DSB_Parp1)

# TODO: Need to deal with unbinding correctly - worry about this later
Parameter("kf_DSB_CtIP", 1)
Parameter("kr_DSB_CtIP", 10)
Rule("CtIP_binds_Parp1",
     CtIP(parp1=None) + DSB(b=1) % Parp1(dsb=1, ctip_mre11=None) |
     CtIP(parp1=2) % DSB(b=1) % Parp1(dsb=1, ctip_mre11=2),
     kf_DSB_CtIP, kr_DSB_CtIP)

Parameter("kf_mre11_Parp1",1)
Parameter("kr_mre11_Parp1",10)
Rule("MRE11_binds_Parp1",
     MRE11(parp1_rad50=None) +
     DSB(b=MultiState(1,3)) % Parp1(dsb=1, ctip_mre11=None) % Parp1(dsb=3, ctip_mre11= 2) % CtIP(parp1=2) |
     MRE11(parp1_rad50=4) %
     DSB(b=MultiState(1,3)) % Parp1(dsb=1, ctip_mre11=4) % Parp1(dsb=3, ctip_mre11= 2) % CtIP(parp1=2),
     kf_mre11_Parp1, kr_mre11_Parp1)

Parameter("kf_mre11_binds_rad50", 1)
Parameter("kr_mre11_binds_rad50",1)
Rule("MRE11_binds_RAD50",
     MRE11(parp1_rad50=None) + RAD50(mre11=None, nbs1=None) |
     MRE11(parp1_rad50=1) % RAD50(mre11=1, nbs1=None),
     kf_mre11_binds_rad50, kr_mre11_binds_rad50)

Parameter("kf_nbs1_binds_mre11_rad50", 1)
Parameter("kr_nbs1_binds_mre11_rad50", 1)
Rule("NBS1_binds_MRE11_RAD50",
     MRE11(parp1_rad50=1) % RAD50(mre11=1, nbs1=None) + NBS1(rad50=None) |
     MRE11(parp1_rad50=1) % RAD50(mre11=1, nbs1=2) % NBS1(rad50=2),
     kf_nbs1_binds_mre11_rad50, kr_nbs1_binds_mre11_rad50)

# Parameter("kf_mrn_formation", 1)
# Parameter("kr_mrn_formation", 1)
# Rule("MRN_complex_formation",
#      MRE11(parp1_rad50=None) + RAD50(mre11=None, nbs1=None) + NBS1(rad50=None) |
#      MRE11(parp1_rad50=1) % RAD50(mre11=1, nbs1=2) % NBS1(rad50=2),
#      kf_mrn_formation, kr_mrn_formation)

Parameter("k_RPA_binds_MRE11_Parp1_DSB", 1)
Rule("RPA_binds_MRE11_Parp1_DSB",
     RPA(dsb=None) + MRE11(parp1_rad50=4) % DSB(b=MultiState(1, ANY)) % Parp1(dsb=1, ctip_mre11=4) >>
     RPA(dsb=1) % DSB(b=MultiState(1, ANY)) + MRE11(parp1_rad50=None) + Parp1(dsb=None, ctip_mre11=None),
     k_RPA_binds_MRE11_Parp1_DSB)

Parameter("k_RPA_binds_CtIP_Parp1_DSB", 1)
Rule("RPA_binds_CtIP_Parp1_DSB",
     RPA(dsb=None) + CtIP(parp1=4) % DSB(b=MultiState(1, ANY)) % Parp1(dsb=1, ctip_mre11=4) >>
     RPA(dsb=1) % DSB(b=MultiState(1, ANY)) + CtIP(parp1=None) + Parp1(dsb=None, ctip_mre11=None),
     k_RPA_binds_CtIP_Parp1_DSB)


print(model)
print(model.monomers)
print(model.parameters)
print(model.rules)
