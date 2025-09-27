from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt

Model()

# STEPS:
# 1. Parp1 binds DSB
# 2. MRE11 binds Parp1M
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
Monomer("MRE11", ["parp1_rad50"] )
Monomer("RAD50", ["mre11", "nbs1"])
Monomer("NBS1", ["rad50"])
Monomer("RPA", ["dsb", "parp1"])
Monomer("PolQ", ["dsb", "parp1", "polq"])
Monomer("LIG3", ["dsb", "dsb"])


Parameter("DSB_0", 100)
Parameter("Parp1_0", 1000)
Parameter("CtIP_0", 1000)
Parameter("MRE11_0", 1000)
Parameter("RAD50_0", 100)
Parameter("NBS1_0", 100)
Parameter("RPA_0", 1000)
Parameter("PolQ_0", 1000)
Parameter("LIG3_0", 1000)

Initial(DSB(b=MultiState(None, None)), DSB_0)
Initial(Parp1(dsb=None, ctip_mre11=None), Parp1_0)
Initial(CtIP(parp1=None), CtIP_0)
Initial(MRE11(parp1_rad50=None), MRE11_0)
Initial(RAD50(mre11=None, nbs1=None), RAD50_0)
Initial(NBS1(rad50=None), NBS1_0)
Initial(RPA(dsb=None, parp1=None), RPA_0)
Initial(PolQ(dsb=None, parp1=None, polq=None),PolQ_0)
Initial(LIG3(dsb=MultiState(None, None)), LIG3_0)

# STEP 0: Two-step binding process to form MRN from MRE11, RAD50, and NBS1:
# MRE11 is involved in MMEJ, MRN is involved in HR;
# Decision between pathways depends on relative amount of MRE11 vs. MRN.
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

# STEP 1a: Parp1 binds to DSB
Parameter("kf_DSB_Parp1", 1)
Parameter("kr_DSB_Parp1", 10)
Rule("Parp1_binds_DSB",
     DSB(b=None) + Parp1(dsb=None, ctip_mre11=None) >> DSB(b=1) % Parp1(dsb=1, ctip_mre11=None),
     kf_DSB_Parp1)

# STEP 1b: Parp1 unbinds from DSB
Rule("Parp1_unbinds_DSB",
     DSB(b=MultiState(1,None)) % Parp1(dsb=1, ctip_mre11=None) >>
     DSB(b=MultiState(None,None)) + Parp1(dsb=None, ctip_mre11=None),
     kr_DSB_Parp1)

# STEP 1c: Parp1 unbinds from DSB-Parp1
Rule("Parp1_unbinds_DSB_Parp1",
     DSB(b=MultiState(1,2)) % Parp1(dsb=1, ctip_mre11=None) % Parp1(dsb=2, ctip_mre11=None) >>
     DSB(b=MultiState(1,None)) % Parp1(dsb=1, ctip_mre11=None) + Parp1(dsb=None, ctip_mre11=None),
     kr_DSB_Parp1)

# STEP 2: Parp1 recruits CtIP
Parameter("kf_DSB_CtIP", 1)
Parameter("kr_DSB_CtIP", 10)
Rule("CtIP_unbinds_Parp1",
     CtIP(parp1=None) +
     DSB(b=MultiState(1,2)) % Parp1(dsb=1, ctip_mre11=None) % Parp1(dsb=2, ctip_mre11=None) |
     DSB(b=MultiState(1,2)) % Parp1(dsb=1, ctip_mre11=None) % Parp1(dsb=2, ctip_mre11=3) % CtIP(parp1=3),
     kf_DSB_CtIP, kr_DSB_CtIP)

# STEP 3: Parp1 recruits MRE11 after CtIP is bound
Parameter("kf_mre11_Parp1",1)
Parameter("kr_mre11_Parp1",10)
Rule("MRE11_binds_Parp1",
     MRE11(parp1_rad50=None) +
     DSB(b=MultiState(1,3)) % Parp1(dsb=1, ctip_mre11=None) % Parp1(dsb=3, ctip_mre11= 2) % CtIP(parp1=2) |
     MRE11(parp1_rad50=4) %
     DSB(b=MultiState(1,3)) % Parp1(dsb=1, ctip_mre11=4) % Parp1(dsb=3, ctip_mre11= 2) % CtIP(parp1=2),
     kf_mre11_Parp1, kr_mre11_Parp1)

# STEP 4a: RPA displaces Parp1 % CtIP
Parameter("k_RPA_binds_CtIP_Parp1_DSB", 1)
Rule("RPA_binds_CtIP_Parp1_DSB",
     RPA(dsb=None, parp1=None) + DSB(b=MultiState(1, ANY))
     % Parp1(dsb=1, ctip_mre11=3) % CtIP(parp1=3) % MRE11(parp1_rad50=ANY) >>
     RPA(dsb=1, parp1=5) % DSB(b=MultiState(1, ANY))
     % Parp1(dsb=5, ctip_mre11=3) % CtIP(parp1=3) % MRE11(parp1_rad50=ANY),
     k_RPA_binds_CtIP_Parp1_DSB)

# STEP 4b: RPA displaces Parp1 % MRE11
Parameter("k_RPA_binds_MRE11_Parp1_DSB", 1)
Rule("RPA_binds_MRE11_Parp1_DSB",
     RPA(dsb=None, parp1=None) + DSB(b=MultiState(1, ANY)) % Parp1(dsb=1, ctip_mre11=4) % MRE11(parp1_rad50=4) >>
     RPA(dsb=1,parp1=2) % DSB(b=MultiState(1, ANY)) % Parp1(dsb=2, ctip_mre11=4) % MRE11(parp1_rad50=4),
     k_RPA_binds_MRE11_Parp1_DSB)


# STEP 5: POLQ displaces RPA
Parameter("k_PolQ_displaces_RPA", 1)
Rule("POLQ_displaces_RPA_Parp1_Parp1",
     PolQ(dsb=None, parp1=None, polq=None) + RPA(dsb=3) %
     DSB(b=MultiState(1, 3)) % Parp1(dsb=2, ctip_mre11=ANY) % RPA(dsb=1, parp1=2) >>
     PolQ(dsb=1, parp1=2, polq=None) % RPA(dsb=3) %
     DSB(b=MultiState(1, 3)) % Parp1(dsb=2, ctip_mre11=ANY) + RPA(dsb=None, parp1=None),
     k_PolQ_displaces_RPA)

Rule("POLQ_displaces_RPA_Parp1_POLQ",
     PolQ(dsb=None, parp1=None, polq=None) + PolQ(dsb=3) %
     DSB(b=MultiState(1, 3)) % Parp1(dsb=2, ctip_mre11=ANY) % RPA(dsb=1, parp1=2) >>
     PolQ(dsb=1, parp1=2, polq=None) % PolQ(dsb=3) %
     DSB(b=MultiState(1, 3)) % Parp1(dsb=2, ctip_mre11=ANY) + RPA(dsb=None, parp1=None),
     k_PolQ_displaces_RPA)

#STEP 6: POLQ aligns microhomologies and extends DNA (POLQs bind to each other)
Parameter("k_PolQ_aligns_microhomolgies", 1)
Rule("PolQ_aligns_microhomologies",
     PolQ(dsb=ANY, parp1=ANY, polq=None) % PolQ(dsb=ANY, parp1=ANY, polq=None) >>
     PolQ(dsb=ANY, parp1=ANY, polq=1) % PolQ(dsb=ANY, parp1=ANY, polq=1),
     k_PolQ_aligns_microhomolgies)



#STEP 7: LIG3 binds DSB and seals nicks
Parameter("k_LIG3_binds_DSB", 1)
Rule("LIG3_binds_DSB",
     LIG3(dsb=MultiState(None,None)) + DSB(b=MultiState(1, 2)) %
     PolQ(dsb=1, parp1=3, polq=7) % CtIP(parp1=4) % Parp1(dsb=3, ctip_mre11=4) %
     PolQ(dsb=2, parp1=5, polq=7) % MRE11(parp1_rad50=6) % Parp1(dsb=5, ctip_mre11=6) >>
     LIG3(dsb=MultiState(1,2)) % DSB(b=MultiState(1, 2)) +
     PolQ(dsb=None, parp1=None, polq=None) + CtIP(parp1=None) + Parp1(dsb=None, ctip_mre11=None) +
     PolQ(dsb=None, parp1=None, polq=None) + MRE11(parp1_rad50=None) + Parp1(dsb=None, ctip_mre11=None),
     k_LIG3_binds_DSB)

#STEP 8: DNA is repaired
Parameter("k_LIG3_repairs_DSB", 1)
Rule("LIG3_repairs_DSB",
     LIG3(dsb=MultiState(1,2)) % DSB(b=MultiState(1, 2)) >>
     LIG3(dsb=MultiState(None,None)),
     k_LIG3_repairs_DSB)

Observable("DSB_tot",DSB())


# print(model)
# print(model.monomers)
# print(model.parameters)
# print(model.rules)

tspan=np.linspace(0,10,1001)
sim=ScipyOdeSimulator(model, tspan=tspan, verbose=True)
output=sim.run()

plt.figure(constrained_layout=True)
for obs in model.observables:
     plt.plot(tspan,output.observables[obs.name],lw=2,label=obs.name)
plt.xlabel("time")
plt.ylabel("concentration")
plt.legend(loc="best")

plt.show()

for sp in model.species:
     print(sp)
