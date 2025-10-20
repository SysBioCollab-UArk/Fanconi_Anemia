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


def plot_expt_id_data(expt_data, expt_id, obs_same_plot=False, label_dict=None, legend_suffix=None):

    if label_dict is None:
        label_dict = {}

    legend_suffix = '' if legend_suffix is None \
        else ' (%s)' % label_dict.get(legend_suffix.strip(), legend_suffix.strip())


    observables = np.unique([d['observable'] for d in expt_data if d['expt_id'] == expt_id])
    for obs in observables:
        alt_expt_id = [d['alt_expt_id'] for d in expt_data if d['expt_id'] == expt_id and d['observable'] == obs]
        if len(np.unique(alt_expt_id)) == 1:
            alt_expt_id = alt_expt_id[0]
        else:
            raise Exception('Multiple alternate expt IDs for expt_id %s and observable %s' % (expt_id, obs))

        if obs_same_plot:
            figname = obs
            title = obs
            leg_label = alt_expt_id.strip() + legend_suffix
        else:
            figname = '%s_%s' % (expt_id, obs)
            title = alt_expt_id
            leg_label = obs.strip() + legend_suffix

        plt.figure(figname, constrained_layout=True)
        plt.title(title)

        xvals = [d['time'] for d in expt_data if d['expt_id'] == expt_id and d['observable'] == obs]
        yvals = [d['average'] for d in expt_data if d['expt_id'] == expt_id and d['observable'] == obs]
        yerrs = [d['stderr'] for d in expt_data if d['expt_id'] == expt_id and d['observable'] == obs]
        if len(xvals) > 0:
            plt.errorbar(xvals, yvals, yerr=yerrs, fmt='o', ms=8, capsize=6, label=leg_label)

        time_units = [d['time_units'] for d in expt_data if d['expt_id'] == expt_id and d['observable'] == obs]
        amount_units = [d['amount_units'] for d in expt_data if d['expt_id'] == expt_id and d['observable'] == obs]
        if len(np.unique(time_units)) > 1:
            raise Exception('Multiple time units for expt_id %s and observable %s' % (expt_id, obs))
        if len(np.unique(amount_units)) > 1:
            raise Exception('Multiple amount units for expt_id %s and observable %s' % (expt_id, obs))

        plt.xlabel('time (%s)' % time_units[0])
        plt.ylabel('amount (%s)' % amount_units[0])
        plt.legend(loc='best')


def plot_expt_data(expt_datafiles, expt_ids=None, **kwargs):

    if isinstance(expt_datafiles, str):
        expt_datafiles = [expt_datafiles]

    if isinstance(expt_ids, str):
        expt_ids = [expt_ids]

    if expt_ids is not None and len(expt_datafiles) != len(expt_ids):
        raise Exception("'expt_datafiles' (len=%d) and 'expt_ids' (len=%d) must have the same length." %
                        (len(expt_datafiles), len(expt_ids)))

    # Read in experimental data
    expt_data_list = []
    legend_suffix = []
    for expt_datafile in expt_datafiles:
        print(expt_datafile)
        legend_suffix.append(None if len(expt_datafiles) == 1 else os.path.splitext(os.path.basename(expt_datafile))[0])
        expt_data_list.append(np.genfromtxt(expt_datafile, dtype=None, delimiter=',', names=True, encoding="utf_8_sig"))
        print(expt_data_list[-1].dtype.names)

    # Loop over datasets
    for i, data in enumerate(expt_data_list):
        e_ids = np.unique([d['expt_id'] for d in data]) if expt_ids is None else expt_ids[i]
        for e_id in e_ids:
            print('expt_id:', e_id)
            plot_expt_id_data(data, e_id, legend_suffix=legend_suffix[i], **kwargs)


if __name__ == '__main__':

    ##### PLOT EXPERIMENTAL DATA #####
    '''
    label_dict = {
        'Averbeck1988_FA150': 'FA150',
        'Averbeck1988_Normal': 'Normal'
    }

    datafile_Normal = os.path.join('DATA', 'Averbeck1988_Normal.csv')
    datafile_FA150 = os.path.join('DATA', 'Averbeck1988_FA150.csv')

    plot_expt_data([datafile_Normal, datafile_FA150], obs_same_plot=False , label_dict=label_dict)
    '''

    plot_expt_data(os.path.join('DATA', 'Alcon2024_Fig6B.csv'), obs_same_plot=False)

    plt.show()

    quit()

    ##### PLOT SIMULATION RESULTS #####

    # Experiment A: 44 * 6.4 = 281.6 DNA adducts, 0 ICLs
    # Experiment B: 0 DNA adducts, 2 * 6.4 = 12.8 ICLs
    # Experiment C: 5.3 * 6.4 = 33.92 DNA adducts, 38.7 * 6.4 = 247.68 ICLs

    expt_conds = {'A': {'lesion_val': 44 * 6.4, 'icl_val': 0, 'observables': ['DNA_lesions']},
                  'B': {'lesion_val': 0, 'icl_val': 2 * 6.4, 'observables': ['Interstrand_crosslinks']},
                  'C': {'lesion_val': 5.3 * 6.4, 'icl_val': 38.7 * 6.4, 'observables': ['Interstrand_crosslinks']}}

    # Run simulations
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

    plt.show()
