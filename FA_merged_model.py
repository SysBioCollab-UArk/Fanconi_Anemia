from pysb import *

include_FACore = True
include_HR = True
include_NER = True
include_MMEJ = False

if include_FACore:
    from fanconi_anemia_core_pathway import create_model_elements as create_FACore_ME
if include_HR:
    from homologous_recombination import create_model_elements as create_HR_ME
if include_NER:
    from nucleotide_excision_repair import create_model_elements as create_NER_ME
if include_MMEJ:
    from microhomology_mediated_end_joining import create_model_elements as create_MMEJ_ME

DEFINE_OBSERVABLES = True

Model()

if include_MMEJ:
    raise Exception('MMEJ pathway not yet available')

# shared model elements
if include_FACore or include_HR or include_MMEJ:
    Monomer('DSB', ['b'])  # fanconi_anemia_pathway.py
    # Monomer('DSB', ['b'])  # homologous_recombination.py
    # Monomer("DSB", ["b", "b"])  # microhomology_mediated_end_joining.py  # TODO
    Parameter('DSB_0', 0)
    Initial(DSB(b=None), DSB_0)
    Observable('Double_strand_breaks', DSB())

if include_FACore or include_NER:
    Monomer('Lesion', ['fanc', 'ner'])  # fanconi_anemia_pathway.py
    # Monomer('Lesion', ['fanc', 'ner'])  # nucleotide_excision_repair.py
    Parameter('Lesion_0', 0)
    Initial(Lesion(fanc=None, ner=None), Lesion_0)
    Observable('DNA_lesions', Lesion())

if include_HR or include_NER:
    Monomer("Pol_Zeta", ["dna"])  # homologous_recombination.py
    # Monomer("Pol_Zeta", ["dna"])  # nucleotide_excision_repair.py
    Parameter("Pol_Zeta_0", 100)
    Initial(Pol_Zeta(dna=None), Pol_Zeta_0)
    Observable("Pol_Zeta_free", Pol_Zeta(dna=None))

    Monomer("LIG1", ["dna"])  # homologous_recombination.py
    # Monomer("LIG1", ["dna"])  # nucleotide_excision_repair.py
    Parameter("LIG1_0", 100)
    Initial(LIG1(dna=None), LIG1_0)
    Observable("LIG1_free", LIG1(dna=None))

if include_HR or include_MMEJ:
    Monomer("RPA", ["dsb"])  # homologous_recombination.py
    # Monomer("RPA", ["dsb", "parp1"])  # microhomology_mediated_end_joining.py  # TODO
    Parameter("RPA_0", 100)
    Initial(RPA(dsb=None), RPA_0)
    Observable("RPA_free", RPA(dsb=None))

    Monomer("MRN", ["dsb"])  # homologous_recombination.py
    # Monomer("MRN", ["dsb"]) # microhomology_mediated_end_joining.py
    Parameter("MRN_0", 100)
    Initial(MRN(dsb=None), MRN_0)
    Observable("MRN_free", MRN(dsb=None))

if include_FACore:
    create_FACore_ME(define_observables=DEFINE_OBSERVABLES)
if include_HR:
    create_HR_ME(define_observables=DEFINE_OBSERVABLES)
if include_NER:
    create_NER_ME(define_observables=DEFINE_OBSERVABLES)
if include_MMEJ:
    create_MMEJ_ME(define_observables=DEFINE_OBSERVABLES)


if __name__ == '__main__':
    from pysb.simulator import ScipyOdeSimulator
    import numpy as np
    import matplotlib.pyplot as plt
    import os
    from plot_expt_data import plot_expt_data

    # experimental data
    plot_expt_data(os.path.join('DATA', 'Alcon2024_Fig6B.csv'), obs_same_plot=False, show_plot=True)

    # simulations
    ICL_0.value = 1.9

    # simulation commands
    tspan = np.linspace(0, 250, 251)
    sim = ScipyOdeSimulator(model, tspan, verbose=True)
    result = sim.run()

    plt.figure(constrained_layout=True)
    plt.plot(tspan, result.observables['Interstrand_crosslinks'], lw=2, label='ICLs')
    plt.xlabel('time (min)')
    plt.ylabel('concentration (ng/ul')
    plt.legend(loc='best')

    plt.show()
