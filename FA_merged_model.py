from pysb import *
from fanconi_anemia_core_pathway import create_model_elements as create_FACore_ME
from homologous_recombination import create_model_elements as create_HR_ME
from nucleotide_excision_repair import create_model_elements as create_NER_ME
from microhomology_mediated_end_joining import create_model_elements as create_MMEJ_ME

DEFINE_OBSERVABLES = True

Model()

# FA core pathway model elements
create_FACore_ME(define_observables=DEFINE_OBSERVABLES)

Monomer("Pol_Zeta", ["dna"])  # DNA polymerase zeta
Monomer("Ligase", ["dna"])  # DNA ligase
Parameter("Pol_Zeta_0", 100)
Parameter("Ligase_0", 100)
Initial(Pol_Zeta(dna=None), Pol_Zeta_0)
Initial(Ligase(dna=None), Ligase_0)
if DEFINE_OBSERVABLES:
    Observable("Pol_Zeta_DSB", Pol_Zeta(dna=1) % DSB(b=1))
    Observable("Ligase_DSB", Ligase(dna=1) % DSB(b=1))
    Observable("Pol_Zeta_Lesion", Pol_Zeta(dna=1) % Lesion(ner=1))
    Observable("Ligase_lesion", Ligase(dna=1) % Lesion(ner=1))

# Homologous recombination model elements
# create_HR_ME(define_observables=DEFINE_OBSERVABLES)

# Nucleotide excision repair model elements
create_NER_ME(define_observables=DEFINE_OBSERVABLES)

# Microhomology-mediated end joining model elements
# create_MMEJ_ME(define_observables=DEFINE_OBSERVABLES)


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
