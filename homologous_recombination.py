from pysb import *
from pysb.util import alias_model_components


def create_model_elements(define_observables=True):

    # Monomer("MRN", ["dsb"])  # TODO
    # MRN Complex (Mre11-Rad50-Nbs1) comes in contact with damaged DNA at 5' ends
    #   - fxn is DNA end resection (pick 5' end and cut it back so we have a 3' overhang on either strand)
    # Monomer("RPA", ["dsb"])  # TODO
    # Replication protein A (RPA) cover single strand DNA to prevent DNA from coiling back until rad51 gets there
    #   - this complex of RPA and DNA is called ssDNA nucleofilament
    Monomer("Rad51_BRCA2", ["dsb"])
    # rad 51 and BRCA2 replace RPA (activated by ATP)
    # Strand Invasion: one strand of broken DNA attempts to bind to one sister homolog strand
    # DNA polymerase uses invaded strand as a template to remake broken strand
    # Helicase displaces now extended invading strand which re-base pairs with other damaged strand
    # second damaged strand anneals to complementary strand of DNA for another round of DNA synthesis
    # Sister strand dissociates
    # DNA ligase-1 restores knicks

    # Parameter("MRN_0", 100)
    # Parameter("RPA_0", 100)
    Parameter("Rad51_BRCA2_0", 100)

    Parameter("kf_MRN_DSB", 0.01)
    Parameter("kr_MRN_DSB", 1)
    Parameter("k_RPA_DSB", 1)
    Parameter("k_Rad51_BRCA2_DSB", 1)
    Parameter("k_Pol_Zeta_DSB", 1)
    Parameter("k_LIG1_DSB", 1)
    Parameter('k_LIG1_repairs_DNA', 1)

    alias_model_components()

    # Initial(MRN(dsb=None), MRN_0)
    # Initial(RPA(dsb=None), RPA_0)
    Initial(Rad51_BRCA2(dsb=None), Rad51_BRCA2_0)

    if define_observables:
        # Observable("MRN_free", MRN(dsb=None))
        # Observable("RPA_free", RPA(dsb=None))
        Observable("Rad51_BRCA2_free", Rad51_BRCA2(dsb=None))

    Rule('MRN_binds_DSB', MRN(dsb=None) + DSB(b=None) | MRN(dsb=1) % DSB(b=1), kf_MRN_DSB, kr_MRN_DSB)

    Rule('RPA_binds_DSB', RPA(dsb=None) + MRN(dsb=1) % DSB(b=1) >> RPA(dsb=1) % DSB(b=1) + MRN(dsb=None),
         k_RPA_DSB)

    Rule('Rad51_BRCA2_binds_DSB',
         Rad51_BRCA2(dsb=None) + RPA(dsb=1) % DSB(b=1) >> Rad51_BRCA2(dsb=1) % DSB(b=1) + RPA(dsb=None),
         k_Rad51_BRCA2_DSB)

    Rule('Pol_Zeta_binds_DSB',
         Pol_Zeta(dna=None) + Rad51_BRCA2(dsb=1) % DSB(b=1) >> Pol_Zeta(dna=1) % DSB(b=1) + Rad51_BRCA2(dsb=None),
         k_Pol_Zeta_DSB)

    Rule('LIG1_binds_DSB',
         LIG1(dna=None) + Pol_Zeta(dna=1) % DSB(b=1) >> LIG1(dna=1) % DSB(b=1) + Pol_Zeta(dna=None),
         k_LIG1_DSB)

    Rule('LIG1_repairs_DNA', LIG1(dna=1) % DSB(b=1) >> LIG1(dna=None), k_LIG1_repairs_DNA)

    if define_observables:
        Observable("Pol_Zeta_DSB", Pol_Zeta(dna=1) % DSB(b=1))
        Observable("LIG1_DSB", LIG1(dna=1) % DSB(b=1))


if __name__ == '__main__':
    from pysb.simulator import ScipyOdeSimulator
    import numpy as np
    import matplotlib.pyplot as plt

    Model()

    Monomer('DSB', ['b'])
    Parameter("DSB_0", 100)
    Initial(DSB(b=None), DSB_0)
    Observable("DSB_tot", DSB())

    Monomer("Pol_Zeta", ["dna"])
    Parameter("Pol_Zeta_0", 100)
    Initial(Pol_Zeta(dna=None), Pol_Zeta_0)
    Observable("Pol_Zeta_free", Pol_Zeta(dna=None))

    Monomer("LIG1", ["dna"])
    Parameter("LIG1_0", 100)
    Initial(LIG1(dna=None), LIG1_0)
    Observable("LIG1_free", LIG1(dna=None))

    Monomer("RPA", ["dsb"])
    Parameter("RPA_0", 100)
    Initial(RPA(dsb=None), RPA_0)
    Observable("RPA_free", RPA(dsb=None))

    Monomer("MRN", ["dsb"])
    Parameter("MRN_0", 100)
    Initial(MRN(dsb=None), MRN_0)
    Observable("MRN_free", MRN(dsb=None))

    create_model_elements()

    tspan = np.linspace(0,10,1001)
    sim = ScipyOdeSimulator(model,tspan,verbose=True)
    output = sim.run()

    plt.figure(constrained_layout=True)
    for obs in model.observables:
        plt.plot(tspan, output.observables[obs.name], lw=2, label=obs.name)
    plt.xlabel("time", fontsize=16)
    plt.ylabel("concentration", fontsize=16)
    plt.tick_params(labelsize=14)
    plt.legend(loc='best', frameon=False, fontsize=14)

    plt.show()
