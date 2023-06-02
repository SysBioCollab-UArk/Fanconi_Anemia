from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt

# This simple model is for the FA Digital Twin perspective article being written together with Carsten Carlberg,
# Eunike Velleuer, Alfredo Rodriguez, and Elisa Huttinger.

Model()

# 1. M + ICL <-> M:ICL
# 2. AG20 + BL100 <-> AG20:BL100
# 3. AG20:BL100 + CEF <-> AG20:BL100:CEF
# 4. AG20 + CEF <-> AG20:CEF
# 5. AG20:CEF + BL100 <-> AG20:BL100:CEF
# 6. BL100 + CEF <-> BL100:CEF
# 7. BL100:CEF + AG20 <-> AG20:BL100:CEF
# 8. AG20:BL100:CEF + M:ICL <-> AG20:BL100:CEF:M:ICL

Monomer('M', ['icl', 'ag20', 'bl100', 'cef'])
Monomer('ICL', ['m'])
Monomer('AG20', ['bl100', 'cef', 'm'])
Monomer('BL100', ['ag20', 'cef', 'm'])
Monomer('CEF', ['ag20', 'bl100', 'm'])

Observable('FANCM_free', M(icl=None, ag20=None, bl100=None, cef=None))
Observable('FANCM_ICL', M(icl=ANY, ag20=None, bl100=None, cef=None))
Observable('FA_complex', AG20(bl100=ANY, cef=ANY, m=None))
Observable('FAcpx_M_ICL', AG20(bl100=ANY, cef=ANY, m=ANY))

Initial(M(icl=None, ag20=None, bl100=None, cef=None), Parameter('M_0', 100))
Initial(ICL(m=None), Parameter('ICL_0', 100))
Initial(AG20(bl100=None, cef=None, m=None), Parameter('AG20_0', 100))
Initial(BL100(ag20=None, cef=None, m=None), Parameter('BL100_0', 100))
Initial(CEF(ag20=None, bl100=None, m=None), Parameter('CEF_0', 100))

# 1. M + ICL <-> M:ICL
Parameter('kf_M_binds_ICL', 1)
Parameter('kr_M_binds_ICL', 1)
Rule('M_binds_ICL',
     M(icl=None, ag20=None, bl100=None, cef=None) + ICL(m=None) |
     M(icl=1, ag20=None, bl100=None, cef=None) % ICL(m=1),
     kf_M_binds_ICL, kr_M_binds_ICL)

# 2. AG20 + BL100 <-> AG20:BL100
Parameter('kf_AG20_binds_BL100', 1)
Parameter('kr_AG20_binds_BL100', 1)
Rule('AG20_binds_BL100',
     AG20(bl100=None, cef=None, m=None) + BL100(ag20=None, cef=None, m=None) |
     AG20(bl100=1, cef=None, m=None) % BL100(ag20=1, cef=None, m=None),
     kf_AG20_binds_BL100, kr_AG20_binds_BL100)

# 3. AG20:BL100 + CEF <-> AG20:BL100:CEF
Parameter('kf_AG20_binds_CEF', 1)
Parameter('kr_AG20_binds_CEF', 1)
Rule('AG20_binds_CEF',
     AG20(cef=None, bl100=None, m=None) + CEF(ag20=None, bl100=None, m=None) |
     AG20(cef=1, bl100=None, m=None) % CEF(ag20=1, bl100=None, m=None),
     kf_AG20_binds_CEF, kr_AG20_binds_CEF)

# 4. AG20 + CEF <-> AG20:CEF
Parameter('kf_BL100_binds_CEF', 1)
Parameter('kr_BL100_binds_CEF', 1)
Rule('BL100_binds_CEF',
     BL100(cef=None, ag20=None, m=None) + CEF(bl100=None, ag20=None, m=None) |
     BL100(cef=1, ag20=None, m=None) % CEF(bl100=1, ag20=None, m=None),
     kf_BL100_binds_CEF, kr_BL100_binds_CEF)

# 5. AG20:CEF + BL100 <-> AG20:BL100:CEF
Parameter('kf_AG20_BL100_binds_CEF', 1)
Parameter('kr_AG20_BL100_binds_CEF', 1)
Rule('AG20_BL100_binds_CEF',
     AG20(bl100=1, cef=None, m=None) % BL100(ag20=1, cef=None, m=None) + CEF(ag20=None, bl100=None, m=None) |
     AG20(bl100=1, cef=2, m=None) % BL100(ag20=1, cef=3, m=None) % CEF(ag20=2, bl100=3, m=None),
     kf_AG20_BL100_binds_CEF, kr_AG20_BL100_binds_CEF)

# 6. BL100 + CEF <-> BL100:CEF
Parameter('kf_AG20_CEF_binds_BL100', 1)
Parameter('kr_AG20_CEF_binds_BL100', 1)
Rule('AG20_CEF_binds_BL100',
     AG20(cef=1, bl100=None, m=None) % CEF(ag20=1, bl100=None, m=None) + BL100(ag20=None, cef=None, m=None) |
     AG20(cef=1, bl100=2, m=None) % CEF(ag20=1, bl100=3, m=None) % BL100(ag20=2, cef=3, m=None),
     kf_AG20_CEF_binds_BL100, kr_AG20_CEF_binds_BL100)

# 7. BL100:CEF + AG20 <-> AG20:BL100:CEF
Parameter('kf_BL100_CEF_binds_AG20', 1)
Parameter('kr_BL100_CEF_binds_AG20', 1)
Rule('BL100_CEF_binds_AG20',
     BL100(cef=1, ag20=None, m=None) % CEF(bl100=1, ag20=None, m=None) + AG20(bl100=None, cef=None, m=None) |
     BL100(cef=1, ag20=2, m=None) % CEF(bl100=1, ag20=3, m=None) % AG20(bl100=2, cef=3, m=None),
     kf_BL100_CEF_binds_AG20, kr_BL100_CEF_binds_AG20)

# 8. AG20:BL100:CEF + M:ICL <-> AG20:BL100:CEF:M:ICL
Parameter('kf_AG20_BL100_CEF_binds_M_ICL', 1)
Parameter('kr_AG20_BL100_CEF_binds_M_ICL', 1)
Rule('AG20_BL100_CEF_binds_M_ICL',
     AG20(m=None) % BL100(m=None) % CEF(m=None) + M(icl=ANY, ag20=None, bl100=None, cef=None) |
     AG20(m=1) % BL100(m=2) % CEF(m=3) % M(icl=ANY, ag20=1, bl100=2, cef=3),
     kf_AG20_BL100_CEF_binds_M_ICL, kr_AG20_BL100_CEF_binds_M_ICL)

# simulation
tspan = np.linspace(0, 0.5, 501)
sim = ScipyOdeSimulator(model, tspan, verbose=True)
result = sim.run()

for obs in model.observables:
    plt.plot(tspan, result.observables[obs.name], lw=3, label=obs.name)
plt.xlabel('time', fontsize=18)
plt.ylabel('# of molecules', fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.legend(loc=0, fontsize=16, ncol=2)

plt.tight_layout()
plt.show()
