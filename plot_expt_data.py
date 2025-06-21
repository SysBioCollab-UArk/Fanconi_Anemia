from fanconi_anemia import model
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt
import os
import difflib

def find_common_and_unique(strings):
    """
    Find the common substring among multiple strings and extract unique parts.
    """
    # Step 1: Find the longest common substring across all strings
    common_part = strings[0]
    for s in strings[1:]:
        matcher = difflib.SequenceMatcher(None, common_part, s)
        common_part = ''.join(common_part[match.a: match.a + match.size] for match in matcher.get_matching_blocks()
                              if match.size > 0)

    # Step 2: Extract the unique parts
    unique_parts = [s.replace(common_part, "").strip() for s in strings]

    return common_part, unique_parts


####################################################
# Experiment A: 44 * 6.4 = 281.6 DNA adducts, 0 ICLs
# Experiment B: 0 DNA adducts, 2 * 6.4 = 12.8 ICLs
# Experiment C: 5.3 * 6.4 = 33.92 DNA adducts, 38.7 * 6.4 = 247.68 ICLs

expt_conds = {'A': {'lesion_val': 44 * 6.4, 'icl_val': 0, 'observables': ['DNA_lesions']},
              'B': {'lesion_val': 0, 'icl_val': 2 * 6.4, 'observables': ['Interstrand_crosslinks']},
              'C': {'lesion_val': 5.3 * 6.4, 'icl_val': 38.7 * 6.4, 'observables': ['Interstrand_crosslinks']}}

##### Run simulations #####
tspan = np.linspace(0, 24, 241)
sim = ScipyOdeSimulator(model, tspan, verbose=False)
sp_names = [str(sp) for sp in model.species]

# pre-equilibration simulation
output = sim.run(param_values={'DSB_0': 0, 'Lesion_0': 0, 'ICL_0': 0})
initials = output.species[-1]

# get indices for DNA lesions and ICLs
lesion_idx = sp_names.index('Lesion(fanc=None, ner=None)')
icl_idx = sp_names.index('ICL(b=None)')

# loop over experiments
for expt in expt_conds.keys():
    print('Experiment:', expt)

    # add mutations and simulate DNA damage repair
    initials[lesion_idx] = expt_conds[expt]['lesion_val']
    initials[icl_idx] = expt_conds[expt]['icl_val']
    output = sim.run(initials=initials)

    plt.figure(expt, constrained_layout=True)
    for obs in expt_conds[expt]['observables']:
        plt.plot(tspan, output.observables[obs], lw=2, label='%s (sim)' % obs)
    plt.xlabel('time (hr)')
    plt.ylabel('# of mutations')
    plt.legend(loc='best')

##### Read in experimental data #####
datafile_Normal = os.path.join('DATA', 'Averbeck1988_Normal.csv')
data_Normal = np.genfromtxt(datafile_Normal, dtype=None, delimiter=',', names=True, encoding="utf_8_sig")
print(data_Normal.dtype.names)

datafile_FA150 = os.path.join('DATA', 'Averbeck1988_FA150.csv')
data_FA150 = np.genfromtxt(datafile_FA150, dtype=None, delimiter=',', names=True, encoding="utf_8_sig")
print(data_FA150.dtype.names)

data_list = [data_Normal, data_FA150]
experiments = np.unique([d['expt_id'] for data in data_list for d in data])

for expt in expt_conds.keys():
    plt.figure(expt, constrained_layout=True)
    expt_list = [x for x in experiments if expt in x]  # [A, AA], [B, BB], etc.
    alt_expt_ids = np.array(
        [np.unique([d['alt_expt_id'] for data in data_list for d in data if d['expt_id'] == e])
         for e in expt_list]).flatten()
    common, unique = find_common_and_unique(alt_expt_ids)
    plt.title(common, fontweight='bold')
    observables = np.unique([d['observable'] for data in data_list for d in data if d['expt_id'] in expt_list])
    for obs in observables:
        for e, alt, cell_type in zip(expt_list, alt_expt_ids, unique):
            for data in data_list:
                xvals = [d['time'] for d in data if d['expt_id'] == e and d['observable'] == obs
                         and d['alt_expt_id'] == alt]
                yvals = [d['average'] for d in data if d['expt_id'] == e and d['observable'] == obs
                         and d['alt_expt_id'] == alt]
                yerrs = [d['stderr'] for d in data if d['expt_id'] == e and d['observable'] == obs
                         and d['alt_expt_id'] == alt]
                if len(xvals) > 0:
                    plt.errorbar(xvals, yvals, yerr=yerrs, fmt='o', ms=8, capsize=6, label='%s (%s)' % (obs, cell_type))
    plt.xlabel('time (hr)')
    plt.ylabel('# of mutations')
    plt.legend(loc='best')

plt.show()
