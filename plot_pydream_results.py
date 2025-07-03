import sys
import os
import numpy as np
import pandas as pd
from util import plot_pydream_results

basepath = 'RESULTS'
directories = ['Averbeck1988_Normal', 'Averbeck1988_FA150']

for dirpath in [os.path.join(basepath, directory) for directory in directories]:

    sys.path.insert(0, dirpath)
    import run_fanconi_pydream as rfp

    this_exp_data_file = os.path.join(dirpath, rfp.exp_data_file)

    calibrator = rfp.ParameterCalibration(rfp.model, this_exp_data_file, rfp.sim_protocols,
                                          priors=rfp.custom_priors, no_sample=rfp.no_sample)

    # get timepoints from experimental data to define `tspans`
    expt_data = pd.read_csv(this_exp_data_file)
    expt_ids = expt_data['expt_id'].unique()
    t_maxes = [expt_data.loc[expt_data['expt_id'] == expt_id, 'time'].max() for expt_id in expt_ids]
    tspans = [np.linspace(0, t_max, int(t_max) * 10 + 1) for t_max in t_maxes]

    plot_pydream_results(dirpath, calibrator, obs_labels=rfp.obs_labels, show_plots=True,
                         plot_ll_args={'cutoff': 2},
                         plot_pd_args={'sharex': 'all'},
                         plot_tc_args={'separate_plots': True, 'save_sim_data': True, 'tspans': tspans},
                         which_plots=3)

    del sys.modules['run_fanconi_pydream']
