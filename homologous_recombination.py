from pysb import *
from pysb.util import alias_model_components


def create_hr_model_elements(define_observables=True):

    Monomer("MRN", ["dsb"])
    # MRN Complex (Mre11-Rad50-Nbs1) comes in contact with damaged DNA at 5' ends
    #   - fxn is DNA end resection (pick 5' end and cut it back so we have a 3' overhang on either strand)
    Monomer("RPA", ["dsb"])
    # Replication protein A (RPA) cover single strand DNA to prevent DNA from coiling back until rad51 gets there
    #   - this complex of RPA and DNA is called ssDNA nucleofilament
    Monomer("Rad51_BRCA2", ["dsb"])
    # rad 51 and BRCA2 replace RPA (activated by ATP)
    # Strand Invasion: one strand of broken DNA attempts to bind to one sister homolog strand
    # DNA polymerase uses invaded strand as a template to remake broken strand
    # Helicase displaces now extended invading strand which re-base pairs with other damaged strand
    # second damaged strand anneals to complementary strand of DNA for another round of DNA synthesis
    # Sister strand dissociates
    # DNA ligase resotres knicks

    Parameter("MRN_0", 100)
    Parameter("RPA_0", 100)
    Parameter("Rad51_BRCA2_0", 100)

    Parameter("kf_MRN_DSB", 0.01)
    Parameter("kr_MRN_DSB", 1)
    Parameter("k_RPA_DSB", 1)
    Parameter("k_Rad51_BRCA2_DSB", 1)
    Parameter("k_Pol_Zeta_DSB", 1)
    Parameter("k_Ligase_DSB", 1)
    Parameter('k_Ligase_repairs_DNA', 1)

    alias_model_components()

    if define_observables:
        Observable("MRN_free", MRN(dsb=None))
        Observable("RPA_free", RPA(dsb=None))
        Observable("Rad51_BRCA2_free", Rad51_BRCA2(dsb=None))

    Initial(MRN(dsb=None), MRN_0)
    Initial(RPA(dsb=None), RPA_0)
    Initial(Rad51_BRCA2(dsb=None), Rad51_BRCA2_0)

    Rule('MRN_binds_DSB', MRN(dsb=None) + DSB(b=None) | MRN(dsb=1) % DSB(b=1), kf_MRN_DSB, kr_MRN_DSB)

    Rule('RPA_binds_DSB', RPA(dsb=None) + MRN(dsb=1) % DSB(b=1) >> RPA(dsb=1) % DSB(b=1) + MRN(dsb=None),
         k_RPA_DSB)

    Rule('Rad51_BRCA2_binds_DSB',
         Rad51_BRCA2(dsb=None) + RPA(dsb=1) % DSB(b=1) >> Rad51_BRCA2(dsb=1) % DSB(b=1) + RPA(dsb=None),
         k_Rad51_BRCA2_DSB)

    Rule('Pol_Zeta_binds_DSB',
         Pol_Zeta(dna=None) + Rad51_BRCA2(dsb=1) % DSB(b=1) >> Pol_Zeta(dna=1) % DSB(b=1) + Rad51_BRCA2(dsb=None),
         k_Pol_Zeta_DSB)

    Rule('Ligase_binds_DSB',
         Ligase(dna=None) + Pol_Zeta(dna=1) % DSB(b=1) >> Ligase(dna=1) % DSB(b=1) + Pol_Zeta(dna=None),
         k_Ligase_DSB)

    Rule('Ligase_repairs_DNA', Ligase(dna=1) % DSB(b=1) >> Ligase(dna=None), k_Ligase_repairs_DNA)


if __name__ == '__main__':
    from pysb.simulator import ScipyOdeSimulator
    import numpy as np
    import matplotlib.pyplot as plt

    Model()

    Monomer('DSB', ['b'])
    Monomer("Pol_Zeta", ["dna"])
    Monomer("Ligase", ["dna"])

    Parameter("DSB_0", 100)
    Parameter("Pol_Zeta_0", 100)
    Parameter("Ligase_0", 100)

    Initial(DSB(b=None), DSB_0)
    Initial(Pol_Zeta(dna=None), Pol_Zeta_0)
    Initial(Ligase(dna=None), Ligase_0)

    Observable("DSB_tot", DSB())
    Observable("Pol_Zeta_free", Pol_Zeta(dna=None))
    Observable("Ligase_free", Ligase(dna=None))

    create_hr_model_elements()

    tspan = np.linspace(0,10,1001)
    sim = ScipyOdeSimulator(model,tspan,verbose=True)
    output = sim.run()

    for obs in model.observables:
        plt.plot(tspan, output.observables[obs.name], lw=2, label=obs.name)
    plt.xlabel("time")
    plt.ylabel("concentration")
    plt.legend(loc=0)
    plt.tight_layout()

    plt.show()
