from pysb import *
from pysb.util import alias_model_components


def create_ner_model_elements():

    # Monomers
    Monomer("RAD23B", ["xpc"])
    Monomer("XPC", ["rad23b", "lesion"])
    Monomer("TFIIH", ["lesion", "xpd"])
    Monomer("XPD", ["tfiih"])
    # Monomer("ERCC1", ["xpf", "lesion"])  # TODO: move binding of XPF, aka FANCQ, to ERCC1 to fanconi_anemia.py
    # Monomer("XPF", ["ercc1"])  # TODO: XPF is another name for FANCQ

    # Initials (Parameters)
    Parameter("RAD23B_0", 100)
    Parameter("XPC_0", 100)
    Parameter("TFIIH_0", 100)
    Parameter("XPD_0", 100)
    # Parameter("ERCC1_0", 100)  # TODO: move to fanconi_anemia.py
    # Parameter("XPF_0", 100)  # TODO: XPF is another name for FANCQ

    Parameter("k_bind_FANCD2_FANCI", 1)
    Parameter("k_unbind_FANCD2_FANCI", 1)
    Parameter('kf_XPC_lesion', 1)
    Parameter('kr_XPC_lesion', 1)
    Parameter('kf_RAD23B_XPC_lesion', 1)
    Parameter('kr_RAD23B_XPC_lesion', 1)
    Parameter('k_TFIIH_lesion', 1)
    Parameter('kf_XPD_TFIIH_lesion', 1)
    Parameter('kr_XPD_TFIIH_lesion', 1)
    # Parameter('kf_ERCC1_XPF', 1)  # TODO: XPF is another name for FANCQ
    # Parameter('kr_ERCC1_XPF', 1)  # TODO: XPF is another name for FANCQ
    Parameter('k_FANCQ_lesion', 1)  # TODO: XPF is another name for FANCQ
    Parameter('k_Pol_Zeta_lesion', 1)
    Parameter('k_Ligase_lesion', 1)
    Parameter('k_unbind_lesion', 1e3)
    Parameter('k_ligase_repairs_lesion', 1)

    alias_model_components()

    Initial(RAD23B(xpc=None), RAD23B_0)
    Initial(XPC(rad23b=None, lesion=None), XPC_0)
    Initial(TFIIH(lesion=None, xpd=None), TFIIH_0)
    Initial(XPD(tfiih=None), XPD_0)
    # Initial(ERCC1(xpf=None, lesion=None), ERCC1_0)  # TODO: move to fanconi_anemia.py
    # Initial(XPF(ercc1=None), XPF_0)  # TODO: XPF is another name for FANCQ

    # Observables
    Observable("XPC_RAD23B_lesion", XPC(lesion=ANY, rad23b=ANY))
    Observable("TFIIH_XPD_lesion", TFIIH(lesion=ANY, xpd=ANY))
    Observable("FANCQ_ERCC1_lesion", FANCQ(ercc1=ANY, b=1) % Lesion(ner=1))  # TODO: XPF is another name for FANCQ

    # Steps
    # 1. XPC binds DNA adduct
    # 2. RAD23B binds XPC
    # 3. TFIIH binds DNA adduct, displaces RAD23B and XPC
    # 4. XPD binds TFIIH

    # 5. FANCQ/ERCC1 binds DNA adduct, displaces XPD and TFIIH
    # 6. XPG binds to the **DNA adduct**

    # 6. Polymerase zeta binds DNA adduct, displaces FANCQ/ERCC1 + XPG
    # 7. DNA ligase binds DNA adduct, displaces polymerase zeta
    # 8. DNA ligase repairs the DNA adduct

    Rule("XPC_binds_Lesion",
         XPC(lesion=None, rad23b=None) + Lesion(ner=None, fanc=ANY) |
         XPC(lesion=2, rad23b=None) % Lesion(ner=2, fanc=ANY),
         kf_XPC_lesion, kr_XPC_lesion)

    Rule('XPC_unbinds_lesion',
         XPC(lesion=2, rad23b=None) % Lesion(ner=2, fanc=None) >>
         XPC(lesion=None, rad23b=None) + Lesion(ner=None, fanc=None),
         k_unbind_lesion)

    Rule("RAD23B_binds_XPC_lesion",
         RAD23B(xpc=None) + XPC(lesion=1, rad23b=None) % Lesion(ner=1, fanc=ANY) |
         RAD23B(xpc=2) % XPC(lesion=1, rad23b=2) % Lesion(ner=1, fanc=ANY),
         kf_RAD23B_XPC_lesion, kr_RAD23B_XPC_lesion)

    Rule('RAD23B_XPC_unbinds_lesion',
         RAD23B(xpc=2) % XPC(lesion=1, rad23b=2) % Lesion(ner=1, fanc=None) >>
         RAD23B(xpc=None) + XPC(lesion=None, rad23b=None) + Lesion(ner=None, fanc=None),
         k_unbind_lesion)

    Rule("TFIIH_binds_lesion",
         TFIIH(lesion=None, xpd=None) + Lesion(ner=2, fanc=ANY) % RAD23B(xpc=1) % XPC(lesion=2, rad23b=1) >>
         TFIIH(lesion=3, xpd=None) % Lesion(ner=3, fanc=ANY) + RAD23B(xpc=None) + XPC(lesion=None, rad23b=None),
         k_TFIIH_lesion)

    Rule('TFIIH_unbinds_lesion',
         TFIIH(lesion=3, xpd=None) % Lesion(ner=3, fanc=None) >>
         TFIIH(lesion=None, xpd=None) + Lesion(ner=None, fanc=None),
         k_unbind_lesion)

    Rule("XPD_binds_TFIIH_lesion",
         XPD(tfiih=None) + TFIIH(lesion=1, xpd=None) % Lesion(ner=1, fanc=ANY) |
         XPD(tfiih=2) % TFIIH(lesion=1, xpd=2) % Lesion(ner=1, fanc=ANY),
         kf_XPD_TFIIH_lesion, kr_XPD_TFIIH_lesion)

    Rule('XPD_TFIIH_unbinds_lesion',
         XPD(tfiih=2) % TFIIH(lesion=1, xpd=2) % Lesion(ner=1, fanc=None) >>
         XPD(tfiih=None) + TFIIH(lesion=None, xpd=None) + Lesion(ner=None, fanc=None),
         k_unbind_lesion)

    # TODO: XPF is another name for FANCQ. Move this rule to fanconi_anemia.py
    # Rule("ERCC1_binds_XPF",
    #      ERCC1(xpf=None, lesion=None) + XPF(ercc1=None) | ERCC1(xpf=1, lesion=None) % XPF(ercc1=1),
    #      kf_ERCC1_XPF, kr_ERCC1_XPF)

    # TODO: XPF is another name for FANCQ
    Rule("FANCQ_ERCC1_binds_Lesion",
         FANCQ(ercc1=ANY, b=None) + Lesion(ner=3, fanc=4) % FANCD2(dna=4) % XPD(tfiih=2) % TFIIH(lesion=3, xpd=2) >>
         FANCQ(ercc1=ANY, b=3) % Lesion(ner=3, fanc=4) % FANCD2(dna=4) + XPD(tfiih=None) + TFIIH(lesion=None, xpd=None),
         k_FANCQ_lesion)

    # TODO: XPF is another name for FANCQ
    Rule('FANCQ_ERCC1_unbinds_Lesion_fanc_unbound',
         FANCQ(ercc1=ANY, b=3) % Lesion(ner=3, fanc=None) >> FANCQ(ercc1=ANY, b=None) + Lesion(ner=None, fanc=None),
         k_unbind_lesion)

    # TODO: XPF is another name for FANCQ
    Rule('FANCQ_ERCC1_unbinds_Lesion_FANCM',
         FANCQ(ercc1=ANY, b=3) % Lesion(ner=3, fanc=2) % FANCM(dna=2) >>
         FANCQ(ercc1=ANY, b=None) + Lesion(ner=None, fanc=2) % FANCM(dna=2),
         k_unbind_lesion)

    # TODO: XPF is another name for FANCQ
    Rule("Pol_Zeta_binds_lesion",
         Pol_Zeta(dna=None) + Lesion(ner=2, fanc=4) % FANCD2(dna=4) % FANCQ(ercc1=ANY, b=2) >>
         Pol_Zeta(dna=3) % Lesion(ner=3, fanc=4) % FANCD2(dna=4) + FANCQ(ercc1=ANY, b=None),
         k_Pol_Zeta_lesion)

    Rule('Pol_Zeta_unbinds_Lesion_fanc_unbound',
         Pol_Zeta(dna=3) % Lesion(ner=3, fanc=None) >> Pol_Zeta(dna=None) + Lesion(ner=None, fanc=None),
         k_unbind_lesion)

    Rule('Pol_Zeta_unbinds_Lesion_FANCM',
         Pol_Zeta(dna=3) % Lesion(ner=3, fanc=1) % FANCM(dna=1) >>
         Pol_Zeta(dna=None) + Lesion(ner=None, fanc=1) % FANCM(dna=1),
         k_unbind_lesion)

    Rule("Ligase_binds_lesion",
         Ligase(dna=None) +  Lesion(ner=1, fanc=4) % FANCD2(dna=4) % Pol_Zeta(dna=1) >>
         Ligase(dna=1) % Lesion(ner=1, fanc=4) % FANCD2(dna=4) + Pol_Zeta(dna=None),
         k_Ligase_lesion)

    Rule('Ligase_unbinds_Lesion_fanc_unbound',
         Ligase(dna=1) % Lesion(ner=1, fanc=None) >> Ligase(dna=None) + Lesion(ner=None, fanc=None),
         k_unbind_lesion)

    Rule('Ligase_unbinds_Lesion_FANCM',
         Ligase(dna=1) % Lesion(ner=1, fanc=2) % FANCM(dna=2) >>
         Ligase(dna=None) + Lesion(ner=None, fanc=2) % FANCM(dna=2),
         k_unbind_lesion)

    Rule("Ligase_Repairs_Lesion",
         Ligase(dna=1) % Lesion(ner=1, fanc=2) % FANCD2(dna=2) >> Ligase(dna=None) + FANCD2(dna=None),
         k_ligase_repairs_lesion)


if __name__ == "__main__":
    from pysb.simulator import ScipyOdeSimulator
    import numpy as np
    import matplotlib.pyplot as plt

    Model()

    Monomer('Lesion', ['fanc', 'ner'])
    Monomer("Pol_Zeta", ["dna"])
    Monomer("Ligase", ["dna"])

    Parameter("Lesion_0", 100)
    Parameter("Pol_Zeta_0", 100)
    Parameter("Ligase_0", 100)

    Initial(Lesion(fanc=None, ner=None), Lesion_0)
    Initial(Pol_Zeta(dna=None), Pol_Zeta_0)
    Initial(Ligase(dna=None), Ligase_0)

    Observable("Lesion_tot", Lesion())
    Observable("Pol_Zeta_Lesion", Pol_Zeta(dna=ANY))
    Observable("Ligase_lesion", Ligase(dna=ANY))

    # FANCQ binds to ERCC1
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

    # FANCM and FANCD2 bind to the Lesion (simplified rules to make this example work)
    Monomer("FANCM", ['dna', "facpx"])
    Monomer("FANCD2", ["fanci", "facpx", 'fancp', 'dna', "state"], {"state": ["x", "ub"]})
    Parameter("FANCM_0", 100)
    Parameter("FANCD2_0", 200)
    Initial(FANCM(dna=None, facpx=None), FANCM_0)
    Initial(FANCD2(fanci=None, facpx=None, fancp=None, dna=None, state="x"), FANCD2_0)
    Parameter('kf_FANCM_Lesion', 1)
    Parameter('kr_FANCM_Lesion', 1)
    Rule('FANCM_binds_Lesion', FANCM(dna=None) + Lesion(fanc=None) | FANCM(dna=1) % Lesion(fanc=1),
         kf_FANCM_Lesion, kr_FANCM_Lesion)
    Parameter('kf_FANCD2_Lesion', 1)
    Parameter('kr_FANCD2_Lesion', 1)
    Rule('FANCD2_binds_Lesion',
         FANCD2(dna=None) + Lesion(fanc=1) % FANCM(dna=1) | FANCD2(dna=1) % Lesion(fanc=1) + FANCM(dna=None),
         kf_FANCD2_Lesion, kr_FANCD2_Lesion)

    create_ner_model_elements()

    tspan = np.linspace(0,100,101)
    sim = ScipyOdeSimulator(model,tspan,verbose=True)
    output = sim.run()

    for obs in model.observables:
        plt.plot(tspan, output.observables[obs.name], lw=2, label=obs.name)
    plt.xlabel("time")
    plt.ylabel("concentration")
    plt.legend(loc=0)
    plt.tight_layout()
    plt.show()




