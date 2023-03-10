from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt

Model()

Monomer('ICL', ['b'])  # stalled replication fork
Monomer('DSB')
Monomer('Lesion', ['rev1'])

Parameter('ICL_0', 100)
Initial(ICL(b=0), ICL_0)

# Proteins

Monomer("FANCA", ["fancg", "faap20"])
Monomer("FANCG", ["fanca"])
Monomer("FAAP20", ["fanca"])

Monomer("FANCB", ["fancl", "faap100"])
Monomer("FANCL", ["fancb", "faap100"])
Monomer("FAAP100", ["fancb", "fancl"])

Monomer("FANCF", ["fancc"])
Monomer("FANCC", ["fancf", "fance"])
Monomer("FANCE", ["fancc"])

# Complexes

Monomer('AG20', ['bl100', 'cef'])
Monomer('BL100', ['ag20', 'cef'])
Monomer('CEF', ['ag20', 'bl100'])

# Formation of AG20 complex

Parameter("FANCA_0", 100)
Parameter("FANCG_0", 80)
Parameter("FAAP20_0", 120)

Initial(FANCA(fancg=None, faap20=None), FANCA_0)
Initial(FANCG(fanca=None), FANCG_0)
Initial(FAAP20(fanca=None), FAAP20_0)

Parameter('kf_AG', 1)
Parameter('kr_AG', 1)
Parameter('kf_A20', 1)
Parameter('kr_A20', 1)

Rule("FANCA_binds_FANCG", FANCA(fancg=None) + FANCG(fanca=None) | FANCA(fancg=1) % FANCG(fanca=1), kf_AG, kr_AG)
Rule("FANCA_binds_FAAP20", FANCA(faap20=None) + FAAP20(fanca=None) | FANCA(faap20=1) % FAAP20(fanca=1), kf_A20, kr_A20)

Parameter('k_AG20_lump', 1e9)
Rule('AG20_lump', FANCA(fancg=1, faap20=2) % FANCG(fanca=1) % FAAP20(fanca=2) >> AG20(bl100=None, cef=None),
     k_AG20_lump)

Rule('AG20_FAAP20_unlump', AG20(bl100=None, cef=None) >>
     FANCA(fancg=1, faap20=None) % FANCG(fanca=1) + FAAP20(fanca=None), kr_A20)
Rule('AG20_FANCG_unlump', AG20(bl100=None, cef=None) >>
     FANCA(fancg=None, faap20=1) % FAAP20(fanca=1) + FANCG(fanca=None), kr_AG)

# Observable("FANCA_free", FANCA(fancg=None, faap20=None))
# Observable("FANCG_free", FANCG(fanca=None))
# Observable("FAAP20_free", FAAP20(fanca=None))
# Observable("AG20_complex", FANCA(fancg=1, faap20=2) % FANCG(fanca=1) % FAAP20(fanca=2))

# Formation of BL100 complex

Parameter("FANCB_0", 100)
Parameter("FANCL_0", 80)
Parameter("FAAP100_0", 120)

Initial(FANCB(fancl=None, faap100=None), FANCB_0)
Initial(FANCL(fancb=None, faap100=None), FANCL_0)
Initial(FAAP100(fancb=None, fancl=None), FAAP100_0)

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

Parameter('k_BL100_lump', 1e9)
Rule('BL100_formation', FANCB(fancl=1, faap100=2) % FANCL(fancb=1, faap100=3) % FAAP100(fancb=2, fancl=3) >>
     BL100(ag20=None, cef=None), k_BL100_lump)

Rule('BL100_FAAP100_unlump', BL100(ag20=None, cef=None) >>
     FANCB(fancl=1, faap100=None) % FANCL(fancb=1, faap100=None) + FAAP100(fancb=None, fancl=None), kr_BL_100)
Rule('BL100_FANCL_unlump', BL100(ag20=None, cef=None) >>
     FANCB(fancl=None, faap100=1) % FAAP100(fancb=1, fancl=None) + FANCL(fancb=None, faap100=None), kr_B100_L)
Rule('BL100_FANCB_unlump', BL100(ag20=None, cef=None) >>
     FANCL(fancb=None, faap100=1) % FAAP100(fancb=None, fancl=1) + FANCB(fancl=None, faap100=None), kr_L100_B)

# Observable("FANCB_free", FANCB(fancl=None, faap100=None))
# Observable("FANCL_free", FANCL(fancb=None, faap100=None))
# Observable("FAAP100_free", FAAP100(fancb=None, fancl=None))
# Observable("BL100_complex", FANCB(fancl=1, faap100=2) % FANCL(fancb=1, faap100=3) % FAAP100(fancb=2, fancl=3))

# Formation of CEF complex

Parameter("FANCC_0", 150)
Parameter("FANCE_0", 100)
Parameter("FANCF_0", 80)

Initial(FANCC(fance=None, fancf=None), FANCC_0)
Initial(FANCE(fancc=None), FANCE_0)
Initial(FANCF(fancc=None), FANCF_0)

Parameter('kf_CE', 1)
Parameter('kr_CE', 1)
Parameter('kf_CF', 1)
Parameter('kr_CF', 1)

Rule("FANCC_binds_FANCE", FANCC(fance=None) + FANCE(fancc=None) | FANCC(fance=1) % FANCE(fancc=1), kf_CE, kr_CE)
Rule("FANCC_binds_FANCF", FANCC(fancf=None) + FANCF(fancc=None) | FANCC(fancf=1) % FANCF(fancc=1), kf_CF, kr_CF)

Parameter('k_CEF_lump', 1e9)
Rule('CEF_lump', FANCC(fance=1, fancf=2) % FANCE(fancc=1) % FANCF(fancc=2) >> CEF(ag20=None, bl100=None),
     k_CEF_lump)

Rule('CEF_FANCF_unlump', CEF(bl100=None, ag20=None) >>
     FANCC(fance=1, fancf=None) % FANCE(fancc=1) + FANCF(fancc=None), kr_CF)
Rule('CEF_FANCE_unlump', CEF(bl100=None, ag20=None) >>
     FANCC(fance=None, fancf=1) % FANCF(fancc=1) + FANCE(fancc=None), kr_CE)

# Observable("FANCC_free", FANCC(fance=None, fancf=None))
# Observable("FANCE_free", FANCE(fancc=None))
# Observable("FANCF_free", FANCF(fancc=None))
# Observable("CEF_complex", FANCC(fance=1, fancf=2) % FANCE(fancc=1) % FANCF(fancc=2))

# Formation of FA core complex (AG20 % BL100 % CEF)
Observable('AG20_total', AG20())
Observable('BL100_total', BL100())
Observable('CEF_total', CEF())

###########

Parameter('kf_AG20_BL100', 0.1)
Parameter('kr_AG20_BL100', 0.1)
Parameter('kf_AG20_BL100_CEF', 0.1)
Parameter('kr_AG20_BL100_CEF', 0.1)
Parameter('kf_AG20_CEF', 0.1)
Parameter('kr_AG20_CEF', 0.1)
Parameter('kf_AG20_CEF_BL100', 0.1)
Parameter('kr_AG20_CEF_BL100', 0.1)
Parameter('kf_BL100_CEF', 0.1)
Parameter('kr_BL100_CEF', 0.1)
Parameter('kf_BL100_CEF_AG20', 0.1)
Parameter('kr_BL100_CEF_AG20', 0.1)

Rule("AG20_binds_BL100", AG20(bl100=None, cef=None) + BL100(ag20=None, cef=None) |
     AG20(bl100=1, cef=None) % BL100(ag20=1, cef=None), kf_AG20_BL100, kr_AG20_BL100)
Rule("AG20_BL100_binds_CEF", AG20(bl100=1, cef=None) % BL100(ag20=1, cef=None) + CEF(ag20=None, bl100=None) |
     AG20(bl100=1, cef=2) % BL100(ag20=1, cef=3) % CEF(ag20=2, bl100=3), kf_AG20_BL100_CEF, kr_AG20_BL100_CEF)

Rule("AG20_binds_CEF", AG20(bl100=None, cef=None) + CEF(ag20=None, bl100=None) |
     AG20(bl100=None, cef=1) % CEF(ag20=1, bl100=None), kf_AG20_CEF, kr_AG20_CEF)
Rule("AG20_CEF_binds_BL100", AG20(bl100=None, cef=1) % CEF(ag20=1, bl100=None) + BL100(ag20=None, cef=None) |
     AG20(bl100=2, cef=1) % CEF(ag20=1, bl100=3) % BL100(ag20=2, cef=3), kf_AG20_CEF_BL100, kr_AG20_CEF_BL100)

Rule("BL100_binds_CEF", BL100(ag20=None, cef=None) + CEF(ag20=None, bl100=None) |
     BL100(ag20=None, cef=1) % CEF(ag20=None, bl100=1), kf_BL100_CEF, kr_BL100_CEF)
Rule("BL100_CEF_binds_AG20", BL100(ag20=None, cef=1) % CEF(ag20=None, bl100=1) + AG20(bl100=None, cef=None) |
     BL100(ag20=2, cef=1) % CEF(ag20=3, bl100=1) % AG20(bl100=2, cef=3), kf_BL100_CEF_AG20, kr_BL100_CEF_AG20)

Observable("AG20_BL100", AG20(bl100=ANY, cef=None))
Observable("AG20_CEF", AG20(bl100=None, cef=ANY))
Observable("BL100_CEF", BL100(ag20=None, cef=ANY))
# Observable("AG20_BL100_CEF", AG20(bl100=ANY, cef=ANY))

Monomer("FA_complex", ["fancm", "fanct", "fancd2", "rev1"])
Parameter('k_AG20_BL100_CEF_lump', 1e9)
Rule("AG20_BL100_CEF_lump", BL100(ag20=2, cef=1) % CEF(ag20=3, bl100=1) % AG20(bl100=2, cef=3) >>
     FA_complex(fancm=None, fanct=None, fancd2=None, rev1=None),
     k_AG20_BL100_CEF_lump)

Rule('FA_complex_BL100_unlump', FA_complex(fancm=None, fanct=None, fancd2=None, rev1=None) >>
     CEF(ag20=1, bl100=None) % AG20(cef=1, bl100=None) + BL100(ag20=None, cef=None), kr_AG20_CEF_BL100)
Rule('FA_complex_CEF_unlump', FA_complex(fancm=None, fanct=None, fancd2=None, rev1=None) >>
     AG20(cef=None, bl100=1) % BL100(ag20=1, cef=None) + CEF(ag20=None, bl100=None), kr_AG20_BL100_CEF)
Rule('FA_complex_AG20_unlump', FA_complex(fancm=None, fanct=None, fancd2=None, rev1=None) >>
     BL100(ag20=None, cef=1) % CEF(ag20=None, bl100=1) + AG20(bl100=None, cef=None), kr_BL100_CEF_AG20)

Observable("FA_complex_free", FA_complex(fancm=None, fanct=None, fancd2=None, rev1=None))

# 0. FANCM binds to DNA at damage site
Monomer("FANCM", ['icl', "facpx"])
Parameter("FANCM_0", 10)
Initial(FANCM(facpx=None), FANCM_0)
Parameter('kf_M_bind_icl', 1)
Parameter('kr_M_bind_icl', 0.01)
Rule('FANCM_binds_ICL', FANCM(icl=None, facpx=None) + ICL(b=None) | FANCM(icl=1, facpx=None) % ICL(b=1),
     kf_M_bind_icl, kr_M_bind_icl)

# 1. FA complex binds FANCM
Parameter("kf_FAcpx_M", 1)
Parameter("kr_FAcpx_M", 1)
Rule('FAcpx_binds_FANCM', FA_complex(fancm=None, fanct=None, fancd2=None, rev1=None) + FANCM(icl=ANY, facpx=None) |
     FA_complex(fancm=1, fanct=None, fancd2=None, rev1=None) % FANCM(icl=ANY, facpx=1), kf_FAcpx_M, kr_FAcpx_M)
Observable("FAcpx_FANCM", FA_complex(fancm=1, fanct=None, fancd2=None, rev1=None) % FANCM(icl=ANY, facpx=1))

# 2. FA complex % FANCM binds FANCT
Monomer("FANCT", ["facpx"])
Parameter("FANCT_0", 10)
Initial(FANCT(facpx=None), FANCT_0)
Parameter("kf_FAcpx_T", 1)
Parameter("kr_FAcpx_T", 1)
Rule('FAcpx_FANCM_binds_FANCT', FA_complex(fancm=1, fanct=None, fancd2=None, rev1=None) % FANCM(icl=ANY, facpx=1)
     + FANCT(facpx=None) | FA_complex(fancm=1, fanct=2, fancd2=None, rev1=None) % FANCM(icl=ANY, facpx=1)
     % FANCT(facpx=2), kf_FAcpx_T, kr_FAcpx_T)
Observable("FAcpx_FANCM_FANCT", FA_complex(fancm=ANY, fanct=ANY, fancd2=None, rev1=None))

# 3. I + D2 <> I % D2
Monomer("FANCI", ["fancd2", 'fancp', "state"], {"state": ["x", "ub"]})
Monomer("FANCD2", ["fanci", "facpx", 'fancp', 'icl', "state"], {"state": ["x", "ub"]})
Parameter("FANCI_0", 50)
Parameter("FANCD2_0", 50)
Initial(FANCI(fancd2=None, fancp=None, state="x"), FANCI_0)
Initial(FANCD2(fanci=None, facpx=None, fancp=None, icl=None, state="x"), FANCD2_0)
Parameter("kf_fanci_fancd2", 1)
Parameter("kr_fanci_fancd2", 1)
Rule('FANCI_binds_FANCD2',
     FANCI(fancd2=None, fancp=None, state="x") + FANCD2(fanci=None, facpx=None, fancp=None, icl=None, state="x") |
     FANCI(fancd2=1, fancp=None, state="x") % FANCD2(fanci=1, facpx=None, fancp=None, icl=None, state="x"),
     kf_fanci_fancd2, kr_fanci_fancd2)
Observable("FANCIx_FANCD2x",
           FANCI(fancd2=1, fancp=None, state="x") % FANCD2(fanci=1, facpx=None, fancp=None, icl=None, state="x"))

# 4. FA complex % FANCM % FANCT binds I % D2
Parameter("kf_FAcpxMT_binds_ID2", 1)
Parameter("kr_FAcpxMT_binds_ID2", 1)
Rule("FAcpx_M_T_binds_I_D2",
     FA_complex(fancm=ANY, fanct=ANY, fancd2=None, rev1=None) +
     FANCD2(fanci=ANY, facpx=None, fancp=None, icl=None, state="x") |
     FA_complex(fancm=ANY, fanct=ANY, fancd2=1, rev1=None) %
     FANCD2(fanci=ANY, facpx=1, fancp=None, icl=None, state="x"),
     kf_FAcpxMT_binds_ID2, kr_FAcpxMT_binds_ID2)
Observable("FAcpx_M_T_Ix_D2x",
           FA_complex(fancm=ANY, fanct=ANY, fancd2=1, rev1=None) %
           FANCD2(fanci=ANY, facpx=1, fancp=None, icl=None, state="x"))

# 5. ubiquitination of I % D2
Parameter("k_ID2_Ubiq", 1)
Rule("FAcpx_I_D2_Ubiq",
     FANCI(fancd2=1, fancp=None, state="x") % FANCD2(fanci=1, facpx=ANY, fancp=None, icl=None, state="x") >>
     FANCI(fancd2=1, fancp=None, state="ub") % FANCD2(fanci=1, facpx=ANY, fancp=None, icl=None, state="ub"), k_ID2_Ubiq)
Observable("FAcpx_FANCIub_FANCD2ub",
           FANCI(fancd2=1, fancp=None, state="ub") % FANCD2(fanci=1, facpx=ANY, fancp=None, icl=None, state="ub"))

# 6. release of ID2-Ub from FA complex % FANCM % FANCT
Parameter("k_FAcpxMT_release_ID2ub", 10)
Rule("FAcpx_M_T_release_ID2ub",
     FA_complex(fancm=3, fanct=ANY, fancd2=2, rev1=None) % FANCM(facpx=3, icl=4) % ICL(b=4) %
     FANCI(fancd2=1, fancp=None, state="ub") % FANCD2(fanci=1, facpx=2, fancp=None, state="ub") >>
     FA_complex(fancm=3, fanct=ANY, fancd2=None, rev1=None) % FANCM(facpx=3, icl=None) +
     FANCI(fancd2=1, fancp=None, state="ub") % FANCD2(fanci=1, facpx=None, fancp=None, icl=2, state="ub") % ICL(b=2),
     k_FAcpxMT_release_ID2ub)
Observable("ID2_Ub", FANCI(fancd2=1, fancp=None, state="ub") % FANCD2(fanci=1, facpx=None, fancp=None, state="ub"))

# 7. Reversible binding of UAF1 and USP1
Monomer('UAF1', ['usp1'])
Monomer('USP1', ['uaf1'])
Parameter('UAF1_0', 100)
Parameter('USP1_0', 100)
Initial(UAF1(usp1=None), UAF1_0)
Initial(USP1(uaf1=None), USP1_0)
Parameter('kf_uaf1_usp1', 0.01)
Parameter('kr_uaf1_usp1', 0.5)
Rule('UAF1_binds_USP1', UAF1(usp1=None) + USP1(uaf1=None) | UAF1(usp1=1) % USP1(uaf1=1), kf_uaf1_usp1, kr_uaf1_usp1)
Observable('UAF1_USP1', UAF1(usp1=1) % USP1(uaf1=1))

# 8. Deubiquitination of ID2-Ub by UAF1 and USP1
Parameter('k_ID2_deubiq', 0.01)
Rule('ID2_deubiqitination',
     FANCI(fancd2=1, fancp=None, state="ub") % FANCD2(fanci=1, facpx=None, fancp=None, icl=2, state="ub") % ICL(b=2) +
     UAF1(usp1=1) % USP1(uaf1=1) >>
     FANCI(fancd2=1, fancp=None, state="x") % FANCD2(fanci=1, facpx=None, fancp=None, icl=None, state="x") + ICL(b=None)
     + UAF1(usp1=1) % USP1(uaf1=1), k_ID2_deubiq)

# 9. FANCP binds to ID2-Ub
Monomer('FANCP', ['fanci', 'fancd2', 'fancq'])
Parameter('FANCP_0', 50)
Initial(FANCP(fanci=None, fancd2=None, fancq=None), FANCP_0)
Parameter('kf_fancp_ID2Ub', 1)
Parameter('kr_fancp_ID2Ub', 10)
Rule('FANCP_binds_ID2Ub',
     FANCP(fanci=None, fancd2=None) +
     FANCI(fancd2=1, fancp=None, state="ub") % FANCD2(fanci=1, facpx=None, fancp=None, icl=ANY, state="ub") |
     FANCP(fanci=2, fancd2=3) %
     FANCI(fancd2=1, fancp=2, state="ub") % FANCD2(fanci=1, facpx=None, fancp=3, icl=ANY, state="ub"),
     kf_fancp_ID2Ub, kr_fancp_ID2Ub)
Observable('FANCP_ID2Ub', FANCP(fanci=ANY, fancd2=ANY, fancq=None))

# 10. FANCQ binds to ID2-Ub % FANCP
Monomer('FANCQ', ['fancp'])
Parameter('FANCQ_0', 100)
Initial(FANCQ(fancp=None), FANCQ_0)
Parameter('kf_fancq_fancp', 1)
Parameter('kr_fancq_fancp', 100)
Rule('FANCQ_binds_FANCP_ID2Ub',
     FANCQ(fancp=None) + FANCP(fanci=ANY, fancd2=ANY, fancq=None) |
     FANCQ(fancp=1) % FANCP(fanci=ANY, fancd2=ANY, fancq=1), kf_fancq_fancp, kr_fancq_fancp)
Observable('FANCQ_FANCP_ID2Ub', FANCQ(fancp=ANY))

# 11. Unhooking machinery (FANCP, FANCQ) produce double-strand break and DNA lesion
Parameter('k_unhook', 10)
Rule('DSB_and_DNA_lesion_creation', FANCQ(fancp=1) % FANCP(fanci=ANY, fancd2=2, fancq=1) % FANCD2(fancp=2, icl=3)
     % ICL(b=3) >> FANCQ(fancp=1) % FANCP(fanci=ANY, fancd2=2, fancq=1) % FANCD2(fancp=2, icl=None)
     + DSB() + Lesion(rev1=None), k_unhook)

# TODO: Add rules for double strand break repair AND translesion synthesis
# TODO: Sabrina = DSB repair pathway
# TODO: Alyssa = DNA lesion repair pathway

# ...

# DSB >> None ==> DNA repair

# simulation commands

tspan = np.linspace(0, 10, 101)
sim = ScipyOdeSimulator(model, tspan, verbose=True)
result = sim.run()

for obs in model.observables:
    plt.plot(tspan, result.observables[obs.name], lw=2, label=obs.name)
plt.xlabel('time')
plt.ylabel('number of molecules')
plt.legend(loc="best", ncol=2, fontsize=6)

plt.show()
