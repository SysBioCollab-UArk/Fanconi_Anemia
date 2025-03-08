from pysb import *
from pysb.simulator import ScipyOdeSimulator
from homologous_recombination import create_hr_model_elements
from nucleotide_excision_repair import create_ner_model_elements
import numpy as np
import matplotlib.pyplot as plt

Model()

# Shared components

Monomer('ICL', ['b'])  # interstrand crosslinks
Monomer('DSB', ['b'])  # double strand breaks
Monomer('Lesion', ['fanc', 'ner'])  # DNA lesions
Monomer("Pol_Zeta", ["dna"])  # DNA polymerase zeta
Monomer("Ligase", ["dna"])  # DNA ligase

Parameter("ICL_0", 0)
Parameter('DSB_0', 0)
Parameter('Lesion_0', 0)
Parameter("Pol_Zeta_0", 100)
Parameter("Ligase_0", 100)

Initial(ICL(b=None), ICL_0)
Initial(DSB(b=None), DSB_0)
Initial(Lesion(fanc=None, ner=None), Lesion_0)
Initial(Pol_Zeta(dna=None), Pol_Zeta_0)
Initial(Ligase(dna=None), Ligase_0)

# ICL synthesis rule
Parameter("k_ICL_synth", 0)
Rule("ICL_synthesis", None >> ICL(b=None), k_ICL_synth)

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
Rule("FANCA_binds_FAAP20", FANCA(faap20=None) + FAAP20(fanca=None) | FANCA(faap20=1) % FAAP20(fanca=1),
     kf_A20, kr_A20)

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
Observable("AG20_BL100_CEF", AG20(bl100=ANY, cef=ANY))

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
Observable("FAcpx_free", FA_complex(fancm=None, fanct=None, fancd2=None, rev1=None))

# 0. FANCM binds to DNA at damage site
Monomer("FANCM", ['dna', "facpx"])
Parameter("FANCM_0", 100)
Initial(FANCM(dna=None, facpx=None), FANCM_0)
Parameter('kf_M_bind_icl', 1)
Parameter('kr_M_bind_icl', 100)
Rule('FANCM_binds_ICL', FANCM(dna=None, facpx=None) + ICL(b=None) | FANCM(dna=1, facpx=None) % ICL(b=1),
     kf_M_bind_icl, kr_M_bind_icl)
Observable('FANCM_tot', FANCM())
Observable('FANCM_free', FANCM(dna=None, facpx=None))
Observable('FANCM_icl', FANCM(dna=1, facpx=None) % ICL(b=1))

Parameter('kf_M_bind_lesion', 1)
Parameter('kr_M_bind_lesion', 100)
Rule('FANCM_binds_Lesion',
     FANCM(dna=None, facpx=None) + Lesion(fanc=None) | FANCM(dna=1, facpx=None) % Lesion(fanc=1),
     kf_M_bind_lesion, kr_M_bind_lesion)
Observable('FANCM_lesion', FANCM(dna=1, facpx=None) % Lesion(fanc=1))

# 1. FA complex binds FANCM
Parameter("kf_FAcpx_M", 1)
Parameter("kr_FAcpx_M", 10)
Rule('FAcpx_binds_FANCM',
     FA_complex(fancm=None, fanct=None, fancd2=None, rev1=None) + FANCM(dna=ANY, facpx=None) >>
     FA_complex(fancm=1, fanct=None, fancd2=None, rev1=None) % FANCM(dna=ANY, facpx=1), kf_FAcpx_M)
Rule('FAcpx_ubinds_FANCM',
     FA_complex(fancm=1, fanct=None, fancd2=None, rev1=None) % FANCM(facpx=1) >>
     FA_complex(fancm=None, fanct=None, fancd2=None, rev1=None) + FANCM(facpx=None), kr_FAcpx_M)
Observable("FAcpx_FANCM", FA_complex(fancm=ANY, fanct=None, fancd2=None, rev1=None))

# 2. FA complex % FANCM binds FANCT
Monomer("FANCT", ["facpx"])
Parameter("FANCT_0", 100)
Initial(FANCT(facpx=None), FANCT_0)
Parameter("kf_FAcpx_T", 1)
Parameter("kr_FAcpx_T", 10)
Rule('FAcpx_FANCM_binds_FANCT',
     FA_complex(fancm=1, fanct=None, fancd2=None, rev1=None) % FANCM(dna=ANY, facpx=1) + FANCT(facpx=None) >>
     FA_complex(fancm=1, fanct=2, fancd2=None, rev1=None) % FANCM(dna=ANY, facpx=1) % FANCT(facpx=2), kf_FAcpx_T)
Rule('FAcpx_FANCM_unbinds_FANCT',
     FA_complex(fancm=1, fanct=2, fancd2=None, rev1=None) % FANCM(facpx=1) % FANCT(facpx=2) >>
     FA_complex(fancm=1, fanct=None, fancd2=None, rev1=None) % FANCM(facpx=1) + FANCT(facpx=None), kr_FAcpx_T)
Observable('FANCT_tot', FANCT())
Observable('FANCT_free', FANCT(facpx=None))
Observable("FAcpx_FANCM_FANCT", FA_complex(fancm=ANY, fanct=ANY, fancd2=None, rev1=None))

# 3. I + D2 <> I % D2
Monomer("FANCI", ["fancd2", 'fancp', "state"], {"state": ["x", "ub"]})
Monomer("FANCD2", ["fanci", "facpx", 'fancp', 'dna', "state"], {"state": ["x", "ub"]})
Parameter("FANCI_0", 200)
Parameter("FANCD2_0", 200)
Initial(FANCI(fancd2=None, fancp=None, state="x"), FANCI_0)
Initial(FANCD2(fanci=None, facpx=None, fancp=None, dna=None, state="x"), FANCD2_0)
Parameter("kf_fanci_fancd2", 1)
Parameter("kr_fanci_fancd2", 1)
Rule('FANCI_binds_FANCD2',
     FANCI(fancd2=None, fancp=None, state="x") + FANCD2(fanci=None, facpx=None, fancp=None, dna=None, state="x") |
     FANCI(fancd2=1, fancp=None, state="x") % FANCD2(fanci=1, facpx=None, fancp=None, dna=None, state="x"),
     kf_fanci_fancd2, kr_fanci_fancd2)
Observable('FANCI_tot', FANCI())
Observable('FANCD2_tot', FANCD2())
Observable("FANCIx_FANCD2x",
           FANCI(fancd2=1, fancp=None, state="x") % FANCD2(fanci=1, facpx=None, fancp=None, dna=None, state="x"))

# 4. FA complex % FANCM % FANCT binds I % D2
Parameter("kf_FAcpxMT_binds_ID2", 1)
Parameter("kr_FAcpxMT_binds_ID2", 1)
Rule("FAcpx_M_T_binds_I_D2",
     FA_complex(fancm=1, fanct=ANY, fancd2=None, rev1=None) % FANCM(dna=ANY, facpx=1) +
     FANCD2(fanci=ANY, facpx=None, fancp=None, dna=None, state="x") |
     FA_complex(fancm=1, fanct=ANY, fancd2=2, rev1=None) % FANCM(dna=ANY, facpx=1) %
     FANCD2(fanci=ANY, facpx=2, fancp=None, dna=None, state="x"),
     kf_FAcpxMT_binds_ID2, kr_FAcpxMT_binds_ID2)
Observable("FAcpx_M_T_Ix_D2x",
           FA_complex(fancm=ANY, fanct=ANY, fancd2=1, rev1=None) %
           FANCD2(fanci=ANY, facpx=1, fancp=None, dna=None, state="x"))

# 5. ubiquitination of I % D2
Parameter("k_ID2_Ubiq", 1)
Rule("FAcpx_I_D2_Ubiq",
     FANCI(fancd2=1, fancp=None, state="x") % FANCD2(fanci=1, facpx=ANY, fancp=None, dna=None, state="x")
     % FANCM(dna=ANY) >>
     FANCI(fancd2=1, fancp=None, state="ub") % FANCD2(fanci=1, facpx=ANY, fancp=None, dna=None, state="ub")
     % FANCM(dna=ANY), k_ID2_Ubiq)
Observable("FAcpx_ID2ub",
           FANCI(fancd2=1, fancp=None, state="ub") % FANCD2(fanci=1, facpx=ANY, fancp=None, dna=None, state="ub"))
Observable('FANCIub_tot', FANCI(state='ub'))
Observable('FANCD2ub_tot', FANCD2(state='ub'))

# 6. release of ID2-Ub from FA complex % FANCM % FANCT
Parameter("k_FAcpxMT_ICL_release_ID2ub", 10)
Rule("FAcpx_M_T_ICL_release_ID2ub",
     FA_complex(fancm=3, fanct=ANY, fancd2=2, rev1=None) % FANCM(facpx=3, dna=4) % ICL(b=4) %
     FANCI(fancd2=1, fancp=None, state="ub") % FANCD2(fanci=1, facpx=2, fancp=None, dna=None, state="ub") >>
     FA_complex(fancm=3, fanct=ANY, fancd2=None, rev1=None) % FANCM(facpx=3, dna=None) +
     FANCI(fancd2=1, fancp=None, state="ub") % FANCD2(fanci=1, facpx=None, fancp=None, dna=2, state="ub") % ICL(b=2),
     k_FAcpxMT_ICL_release_ID2ub)

Parameter("k_FAcpxMT_Lesion_release_ID2ub", 1)
Rule("FAcpx_M_T_Lesion_release_ID2ub",
     FA_complex(fancm=3, fanct=ANY, fancd2=2, rev1=None) % FANCM(facpx=3, dna=4) % Lesion(fanc=4) %
     FANCI(fancd2=1, fancp=None, state="ub") % FANCD2(fanci=1, facpx=2, fancp=None, dna=None, state="ub") >>
     FA_complex(fancm=3, fanct=ANY, fancd2=None, rev1=None) % FANCM(facpx=3, dna=None) + Lesion(fanc=2) %
     FANCI(fancd2=1, fancp=None, state="ub") % FANCD2(fanci=1, facpx=None, fancp=None, dna=2, state="ub"),
     k_FAcpxMT_Lesion_release_ID2ub)

Observable("ID2_Ub",
           FANCI(fancd2=1, fancp=None, state="ub") % FANCD2(fanci=1, facpx=None, fancp=None, state="ub"))

# 7. Reversible binding of UAF1 and USP1
Monomer('UAF1', ['usp1'])
Monomer('USP1', ['uaf1'])
Parameter('UAF1_0', 100)
Parameter('USP1_0', 100)
Initial(UAF1(usp1=None), UAF1_0)
Initial(USP1(uaf1=None), USP1_0)
Parameter('kf_uaf1_usp1', 0.01)  # 0.01
Parameter('kr_uaf1_usp1', 0.5)  # 0.5
Rule('UAF1_binds_USP1',
     UAF1(usp1=None) + USP1(uaf1=None) | UAF1(usp1=1) % USP1(uaf1=1), kf_uaf1_usp1, kr_uaf1_usp1)
Observable('UAF1_USP1', UAF1(usp1=1) % USP1(uaf1=1))

# 8. Deubiquitination of ID2-Ub by UAF1 and USP1
Parameter('k_ID2_deubiq', 0.01)  # 0.01
Rule('ID2_ICL_deubiqitination',
     FANCI(fancd2=1, fancp=None, state="ub") % FANCD2(fanci=1, facpx=None, fancp=None, dna=2, state="ub") %
     ICL(b=2) + UAF1(usp1=1) % USP1(uaf1=1) >>
     FANCI(fancd2=1, fancp=None, state="x") % FANCD2(fanci=1, facpx=None, fancp=None, dna=None, state="x") +
     ICL(b=None) + UAF1(usp1=1) % USP1(uaf1=1), k_ID2_deubiq)

Rule('ID2_Lesion_deubiqitination',
     FANCI(fancd2=1, fancp=None, state="ub") % FANCD2(fanci=1, facpx=None, fancp=None, dna=2, state="ub") %
     Lesion(fanc=2) + UAF1(usp1=1) % USP1(uaf1=1) >>
     FANCI(fancd2=1, fancp=None, state="x") % FANCD2(fanci=1, facpx=None, fancp=None, dna=None, state="x") +
     Lesion(fanc=None) + UAF1(usp1=1) % USP1(uaf1=1), k_ID2_deubiq)

Rule("ID2_free_deubiqitination",
     FANCI(fancd2=1, fancp=None, state="ub") % FANCD2(fanci=1, facpx=None, fancp=None, dna=None, state="ub") +
     UAF1(usp1=1) % USP1(uaf1=1) >>
     FANCI(fancd2=1, fancp=None, state="x") % FANCD2(fanci=1, facpx=None, fancp=None, dna=None, state="x") +
     UAF1(usp1=1) % USP1(uaf1=1), k_ID2_deubiq)

# 9. FANCP binds to ID2-Ub
# TODO: Right now, this happens for ID2-Ub bound both to ICLs and Lesions. Not sure that should be the case.
Monomer('FANCP', ['fanci', 'fancd2', 'fancq'])
Parameter('FANCP_0', 50)
Initial(FANCP(fanci=None, fancd2=None, fancq=None), FANCP_0)
Parameter('kf_fancp_ID2Ub', 1)
Parameter('kr_fancp_ID2Ub', 10)
Rule('FANCP_binds_ID2Ub',
     FANCP(fanci=None, fancd2=None, fancq=None) +
     FANCI(fancd2=1, fancp=None, state="ub") % FANCD2(fanci=1, facpx=None, fancp=None, dna=ANY, state="ub") >>
     FANCP(fanci=2, fancd2=3, fancq=None) %
     FANCI(fancd2=1, fancp=2, state="ub") % FANCD2(fanci=1, facpx=None, fancp=3, dna=ANY, state="ub"),
     kf_fancp_ID2Ub)
Rule('FANCP_unbinds_ID2Ub',
     FANCP(fanci=2, fancd2=3, fancq=None) %
     FANCI(fancd2=1, fancp=2, state="ub") % FANCD2(fanci=1, facpx=None, fancp=3, state="ub") >>
     FANCP(fanci=None, fancd2=None, fancq=None) +
     FANCI(fancd2=1, fancp=None, state="ub") % FANCD2(fanci=1, facpx=None, fancp=None, state="ub"),
     kr_fancp_ID2Ub)
Observable('FANCP_ID2Ub', FANCP(fanci=ANY, fancd2=ANY, fancq=None))

# 9.5. FANCQ binds to ERCC1
Monomer("FANCQ", ["ercc1", 'b'])
Monomer("ERCC1", ["fancq"])
Parameter('FANCQ_0', 100)
Parameter('ERCC1_0', 100)
Initial(FANCQ(ercc1=None, b=None), FANCQ_0)
Initial(ERCC1(fancq=None), ERCC1_0)
Parameter('kf_FANCQ_ERCC1', 1)
Parameter('kr_FANCQ_ERCC1', 1)
Rule("FANCQ_binds_ERCC1",
         FANCQ(ercc1=None, b=None) + ERCC1(fancq=None) | FANCQ(ercc1=1, b=None) % ERCC1(fancq=1),
         kf_FANCQ_ERCC1, kr_FANCQ_ERCC1)
Observable('FANCQ_ERCC1_dimer', FANCQ(ercc1=1, b=None) % ERCC1(fancq=1))

# 10. FANCQ_ERCC1 binds to ID2-Ub % FANCP
Parameter('kf_fancq_fancp', 1)
Parameter('kr_fancq_fancp', 100)
Rule('FANCQ_ERCC1_binds_FANCP_ID2Ub',
     FANCQ(ercc1=ANY, b=None) + FANCP(fanci=ANY, fancd2=1, fancq=None) %
     FANCD2(fanci=ANY, facpx=None, fancp=1, dna=ANY, state="ub") >>
     FANCQ(ercc1=ANY, b=2) % FANCP(fanci=ANY, fancd2=1, fancq=2) %
     FANCD2(fanci=ANY, facpx=None, fancp=1, dna=ANY, state="ub"), kf_fancq_fancp)
Rule('FANCQ_unbinds_FANCP_ID2Ub',
     FANCQ(ercc1=ANY, b=2) % FANCP(fanci=ANY, fancd2=1, fancq=2) %
     FANCD2(fanci=ANY, facpx=None, fancp=1, state="ub") >>
     FANCQ(ercc1=ANY, b=None) + FANCP(fanci=ANY, fancd2=1, fancq=None) %
     FANCD2(fanci=ANY, facpx=None, fancp=1, state="ub"), kr_fancq_fancp)
Observable('FANCQ_FANCP_ID2Ub',
           FANCQ(ercc1=ANY, b=1) % FANCP(fancq=1, fancd2=2) % FANCD2(fancp=2, fanci=ANY, state='ub'))

# 11. Unhooking machinery (FANCP, FANCQ) produce double-strand break and DNA lesion
Parameter('k_unhook', 10)
Rule('DSB_and_DNA_lesion_creation',
     FANCQ(b=1) % FANCP(fanci=ANY, fancd2=2, fancq=1) % FANCD2(fancp=2, dna=3) % ICL(b=3) >>
     FANCQ(b=1) % FANCP(fanci=ANY, fancd2=2, fancq=1) % FANCD2(fancp=2, dna=None) + DSB(b=None) +
     Lesion(ner=None, fanc=None),
     k_unhook)
Observable('Interstrand_crosslinks', ICL())
Observable('Double_strand_breaks', DSB())
Observable('DNA_lesions', Lesion())

# additional Observables for debugging
# TODO: Probably want to remove these
Observable("Lesion_free", Lesion(ner=None, fanc=None))
Observable("Lesion_ner_ANY_fanc_None", Lesion(ner=ANY, fanc=None))
Observable("Lesion_ner_None_fanc_ANY", Lesion(ner=None, fanc=ANY))
Observable("Lesion_ner_ANY_fanc_ANY", Lesion(ner=ANY, fanc=ANY))

# Homologous recombination model elements
create_hr_model_elements()
Observable("Pol_Zeta_DSB", Pol_Zeta(dna=1) % DSB(b=1))
Observable("Ligase_DSB", Ligase(dna=1) % DSB(b=1))

# Nucleotide Excision Repair model elements
create_ner_model_elements()
Observable("Pol_Zeta_Lesion", Pol_Zeta(dna=1) % Lesion(ner=1))
Observable("Ligase_lesion", Ligase(dna=1) % Lesion(ner=1))

if __name__ == '__main__':

    ICL_0.value = 100

    # simulation commands
    tspan = np.linspace(0, 24, 241)
    sim = ScipyOdeSimulator(model, tspan, verbose=True)
    result = sim.run()

    complexes = ['AG20_total', 'BL100_total', 'CEF_total', 'FAcpx_free', 'FAcpx_ID2ub', 'FANCQ_FANCP_ID2Ub']

    mutations = ['Interstrand_crosslinks', 'Double_strand_breaks', 'DNA_lesions']  # , 'Pol_Zeta_DSB', 'Ligase_DSB',
                 # 'Pol_Zeta_Lesion', 'Ligase_lesion']

    plt.figure('complexes')
    plt.figure('mutations')
    for obs in complexes + mutations:
        if obs in mutations:
            plt.figure('mutations')
        elif obs in complexes:
            plt.figure('complexes')
        plt.plot(tspan, result.observables[obs], lw=2, label=obs)

    plt.figure('complexes')
    plt.xlabel('time (arbitrary units)', fontsize=16)
    plt.ylabel('# of molecules', fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(loc="best", ncol=2, fontsize=14)
    plt.tight_layout()

    plt.figure('mutations')
    plt.xlabel('time (arbitrary units)', fontsize=16)
    plt.ylabel('# of molecules', fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(loc="best", fontsize=14)
    plt.tight_layout()

    plt.figure()
    obs2plot=["DNA_lesions", "FANCM_lesion", 'FANCM_tot', 'FANCM_free', "Lesion_free", "Lesion_ner_ANY_fanc_None",
              "Lesion_ner_None_fanc_ANY", "Lesion_ner_ANY_fanc_ANY"]
    for obs in obs2plot:
        plt.plot(tspan, result.observables[obs], lw=2, label=obs)
    plt.plot(tspan, result.observable(Lesion(fanc=None)), lw=2, label='Lesion_fancm_None')
    plt.legend(loc="best", ncol=1)
    plt.tight_layout()

    plt.show()
