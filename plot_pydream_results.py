import os
import importlib
import numpy as np
import pandas as pd
from util import plot_pydream_results

basepath = 'RESULTS'
directories = ['Averbeck1988_Normal', 'Averbeck1988_FA150']

for dirpath in [os.path.join(basepath, directory) for directory in directories]:

    # import everything from the run_fanconi_pydream.py file that's in the results file path
    run_pydream_file = os.path.join(dirpath, 'run_fanconi_pydream.py')
    import_string = run_pydream_file.replace('/', '.').replace('\\', '.').rstrip('.py')
    module = importlib.import_module(str(import_string))  # import the module

    calibrator = module.ParameterCalibration(module.model, module.exp_data_file, module.sim_protocols,
                                             priors=module.custom_priors, no_sample=module.no_sample)

    # get timepoints from experimental data to define `tspans`
    expt_data = pd.read_csv(module.exp_data_file)
    expt_ids = expt_data['expt_id'].unique()
    t_maxes = [expt_data.loc[expt_data['expt_id'] == expt_id, 'time'].max() for expt_id in expt_ids]
    tspans = [np.linspace(0, t_max, int(t_max) * 10 + 1) for t_max in t_maxes]

    plot_pydream_results(dirpath, calibrator, obs_labels=module.obs_labels, show_plots=True,
                         plot_ll_args={'cutoff': 2},
                         plot_pd_args={'sharex': 'all'},
                         plot_tc_args={'separate_plots': True, 'save_sim_data': True, 'tspans': tspans},
                         which_plots=3)
