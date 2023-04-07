from pysb import *

def create_tls_model_elements(Lesion):

    # Monomers
    Monomer("RAD238", ["xpc"])
    Monomer("XPC", ["rad238", "lesion"])
    Monomer("TFIIH", ["lesion", "xpd"])
    Monomer("XPD", ["tfiih"])
    Monomer("ERCC1", ["xpf", "lesion"])
    Monomer("XPF", ["ercc1"])
    Monomer("Pol_Zeta", ["dna"])
    Monomer("Ligase", ["dna"])

    # Initials (Parameters)
    Parameter("RAD238_0", 100)
    Parameter("XPC_0", 100)
    Parameter("TFIIH_0", 100)
    Parameter("XPD_0", 100)
    Parameter("ERCC1_0", 100)
    Parameter("XPF_0", 100)
    Parameter("Pol_Zeta_0", 100)
    Parameter("Ligase_0", 100)

    Initial(RAD238(xpc=None), RAD238_0)
    Initial(XPC(rad238=None, lesion=None), XPC_0)
    Initial(TFIIH(lesion=None, xpd=None), TFIIH_0)
    Initial(XPD(tfiih=None), XPD_0)
    Initial(ERCC1(xpf=None, lesion=None), ERCC1_0)
    Initial(XPF(ercc1=None), XPF_0)
    Initial(Pol_Zeta(dna=None),Pol_Zeta_0)
    Initial(Ligase(dna=None), Ligase_0)

    # Observables
    Observable("XPC_RAD238_lesion", XPC(lesion=ANY, rad238=ANY))
    Observable("TFIIH_XPD_lesion", TFIIH(lesion=ANY, xpd=ANY))
    Observable("ERCC1_XPF_lesion", ERCC1(lesion=ANY, xpf=ANY))
    Observable("Pol_Zeta_Lesion", Pol_Zeta(dna=ANY))
    Observable("Ligase_lesion", Ligase(dna=ANY))



    # Rules (Parameters)
    # steps
    # 1. RAD238 and XPC detect lesion, leave
    # 2. TFIIH binds to DNA, XPD binds to DNA, unwinds
    # 3. ERCC1 and XPF cleave the damaged sections
    # 4. Polymerase zeta fills gap
    # 5. DNA ligase seals the strands

    Parameter('kf_XPC_lesion', 1)
    Parameter('kr_XPC_lesion', 1)
    Rule("XPC_binds_Lesion", XPC(lesion=None, rad238=None) + Lesion(b=None) | XPC(lesion=1, rad238=None) % Lesion(b=1),
         kf_XPC_lesion, kr_XPC_lesion)

    Parameter('kf_RAD238_XPC_lesion', 1)
    Parameter('kr_RAD238_XPC_lesion', 1)
    Rule("RAD238_binds_XPC_lesion",
         RAD238(xpc=None) + XPC(lesion=ANY, rad238=None) | RAD238(xpc=1) % XPC(lesion=ANY, rad238=1),
         kf_RAD238_XPC_lesion, kr_RAD238_XPC_lesion)

    Parameter('k_TFIIH_lesion', 1)
    Rule("TFIIH_binds_lesion", TFIIH(lesion=None, xpd=None) + Lesion(b=2) % RAD238(xpc=1) % XPC(lesion=2, rad238=1) >>
         TFIIH(lesion=3, xpd=None) % Lesion(b=3) + RAD238(xpc=None) + XPC(lesion=None, rad238=None),
         k_TFIIH_lesion)

    Parameter('kf_XPD_TFIIH_lesion', 1)
    Parameter('kr_XPD_TFIIH_lesion', 1)
    Rule("XPD_binds_TFIIH_lesion",
         XPD(tfiih=None) + TFIIH(lesion=ANY, xpd=None) | XPD(tfiih=1) % TFIIH(lesion=ANY, xpd=1),
         kf_XPD_TFIIH_lesion, kr_XPD_TFIIH_lesion)

    Parameter('kf_ERCC1_XPF', 1)
    Parameter('kr_ERCC1_XPF', 1)
    Rule("ERCC1_binds_XPF", ERCC1(xpf=None, lesion=None) + XPF(ercc1=None) | ERCC1(xpf=1, lesion=None) % XPF(ercc1=1),
         kf_ERCC1_XPF, kr_ERCC1_XPF)

    Parameter('k_ERCC1_XPF_lesion', 1)
    Rule("ERCC1_XPF_binds_lesion",
         ERCC1(xpf=1, lesion=None) % XPF(ercc1=1) + Lesion(b=3) % XPD(tfiih=2) % TFIIH(lesion=3, xpd=2) >>
         ERCC1(xpf=1, lesion=3) % XPF(ercc1=1) % Lesion(b=3) + XPD(tfiih=2) % TFIIH(lesion=None, xpd=2),
         k_ERCC1_XPF_lesion)

    Parameter('k_Pol_Zeta_lesion', 1)
    Rule("Pol_Zeta_binds_lesion", Pol_Zeta(dna=None) + ERCC1(xpf=1, lesion=2) % XPF(ercc1=1) % Lesion(b=2) >>
         Pol_Zeta(dna=3) % Lesion(b=3) + ERCC1(xpf=None, lesion=None) + XPF(ercc1=None),
         k_Pol_Zeta_lesion)

    Parameter('k_Ligase_lesion', 1)
    Rule("Ligase_binds_lesion", Ligase(dna=None) + Pol_Zeta(dna=1) % Lesion(b=1) >>
         Ligase(dna=1) % Lesion(b=1) + Pol_Zeta(dna=None), k_Ligase_lesion)

    Parameter('k_ligase_repairs_lesion', 1)
    Rule("Ligase_Repairs_Lesion", Ligase(dna=1) % Lesion(b=1) >> Ligase(dna=None), k_ligase_repairs_lesion)


if __name__ == "__main__":
    from pysb.simulator import ScipyOdeSimulator
    import numpy as np
    import matplotlib.pyplot as plt

    Model()
    Monomer('Lesion', ['b'])
    Parameter("Lesion_0", 100)
    Initial(Lesion(b=None), Lesion_0)
    Observable("Lesion_tot", Lesion())
    create_tls_model_elements(Lesion)

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




