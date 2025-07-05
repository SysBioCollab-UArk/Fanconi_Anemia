import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

basepath = 'RESULTS'
directories = ['Averbeck1988_Normal', 'Averbeck1988_FA150']

# get the number of rows in the figure by finding the number of experiments (# of alt expt ids)
expt_datafiles = [os.path.join(basepath, directory, 'DATA', '%s.csv' % directory) for directory in directories]
all_alt_expt_ids = np.unique([pd.read_csv(expt_datafile)['alt_expt_id'].unique() for expt_datafile in expt_datafiles])
nrows = len(all_alt_expt_ids)
# get the number of columns in the figure by finding the number of observables for each experiment (ncols = max)
ncols = 0
for alt_expt_id in all_alt_expt_ids:
    observables = []
    for expt_datafile in expt_datafiles:
        data = pd.read_csv(expt_datafile)
        these_observables = data[data['alt_expt_id'] == alt_expt_id]['observable'].unique()
        observables += [obs for obs in these_observables if obs not in observables]
    if len(observables) > ncols:
        ncols = len(observables)
# create the figure
fig = plt.figure(constrained_layout=True, figsize=(6.4 * ncols, 4.8 * 0.8 * nrows))

# dictionary for changing strings displayed on the plots
labels_dict = {'Interstrand_crosslinks': 'Interstrand Crosslinks', 'DNA_lesions': 'Monoadducts', 'FA150': 'FA 150'}

# loop over datafiles and create the figure
fig_names = {}
for directory in directories:
    print('Directory:', directory)

    # read experimental data
    expt_data = pd.read_csv(os.path.join(basepath, directory, 'DATA', '%s.csv' % directory))

    # get experiment IDs
    expt_ids = expt_data['expt_id'].unique()
    print('expt_ids:', expt_ids)

    # get time units
    time_units = expt_data['time_units'].unique()
    if len(time_units) != 1:
        raise ValueError('There must be exactly one time unit but got %d:' % len(time_units), time_units)
    time_units = time_units[0]

    # get alt expt IDs for each experiment ID
    alt_expt_ids = {}
    for expt_id in expt_ids:
        alt_expt_id = expt_data[expt_data['expt_id'] == expt_id]['alt_expt_id'].unique()
        if len(alt_expt_id) != 1:
            raise ValueError("Expected exactly one 'alt_expt_id' but got %d:" % len(alt_expt_id), alt_expt_id)
        alt_expt_ids[expt_id] = alt_expt_id[0]

    # read simulation data
    sim_data = pd.read_csv(os.path.join(basepath, directory, 'SIM_DATA.csv'))

    # loop over experiments and plot each observable
    for expt_id in expt_ids:
        observables = sim_data[sim_data['sim_id'] == expt_id]['observable'].unique()
        for obs in observables:
            fig_name = '%s_%s' % (alt_expt_ids[expt_id], obs)
            # get subplot index
            if fig_name in fig_names:
                idx = fig_names[fig_name]
                ax = fig.get_axes()[idx]
            else:
                idx = len(fig_names)
                ax = fig.add_subplot(nrows, ncols, idx + 1)
                fig_names[fig_name] = idx

            # plot simulation data
            sim_data_subset = sim_data[(sim_data['sim_id'] == expt_id) & (sim_data['observable'] == obs)]
            tspan = sim_data_subset['time']
            yvals_min = sim_data_subset['yval_min']
            yvals_max = sim_data_subset['yval_max']
            label = expt_id.split('_')[-1]
            ax.fill_between(tspan, yvals_min, yvals_max, alpha=0.25, label=labels_dict.get(label, label))

            # plot experimental data
            expt_data_subset = expt_data[(expt_data['expt_id'] == expt_id) & (expt_data['observable'] == obs)]
            tdata = expt_data_subset['time']
            average = expt_data_subset['average']
            stderr = expt_data_subset['stderr']
            label = expt_id.split('_')[-1]
            ax.errorbar(tdata, average, yerr=stderr, fmt='o', ms=6, capsize=8, label=labels_dict.get(label, label))

            # get observable amount units
            y_units = expt_data_subset['amount_units'].unique()
            if len(y_units) != 1:
                raise ValueError('There must be exactly one amount unit for observable %s but got %d:' %
                                 (obs, len(y_units)), y_units)
            y_units = y_units[0]

            ax.set_title(' ', fontsize=16)  # reserve space for a single title across the row
            ax.set_xlabel('Time (%s)' % time_units, fontsize=14)
            ax.set_ylabel('%s (%s)' % (labels_dict.get(obs, obs), y_units), fontsize=14)
            ax.tick_params(axis='both', labelsize=14)
            ax.legend(loc='best')

# clean up the plots
DONE = []
for fig_name in fig_names:
    idx = fig_names[fig_name]
    ax = fig.get_axes()[idx]

    # Turn off x-axis and y-axis labels/ticks for all inner subplots
    row, col = divmod(idx, ncols)  # Convert flat index to (row, col)
    if row < nrows - 1:
        ax.set_xlabel('')
        ax.set_xticklabels([])
    if col > 0:
        ax.set_ylabel('')
        ax.set_yticklabels([])

    # Create one title for each row in the figure
    fig.canvas.draw()
    if row not in DONE:
        top_y = ax.get_position().y1  # y-position of the top of the row
        left_idx = row * ncols  # index of left-most subplot in this row
        left_x = fig.get_axes()[left_idx].get_position().x0  # x-position of left side of this row
        right_idx = (row + 1) * ncols - 1  # index of right-most subplot in this row
        right_x = fig.get_axes()[right_idx].get_position().x1  # x-position of right side of this row
        mid_x = left_x + (right_x - left_x) / 2
        title = None
        for alt_expt_id in all_alt_expt_ids:
            if alt_expt_id in fig_name:
                title = alt_expt_id
        fig.text(mid_x, top_y + 0.0025, title, ha='center', va='bottom', fontsize=14, fontweight='bold')
        DONE.append(row)

    # Adjust the legends of each subplot
    handles, labels = ax.get_legend_handles_labels()
    # figure out the indices where each unique label appears in 'labels'
    legend_dict = {}
    for i in range(len(labels)):
        if labels[i] not in legend_dict:
            legend_dict[labels[i]] = [i]
            for j in range(i + 1, len(labels)):
                if labels[j] == labels[i]:
                    legend_dict[labels[i]].append(j)
    # loop through the labels and create the new labels and handles
    new_labels = []
    new_handles = []
    for label in list(dict.fromkeys(labels)):
        new_labels.append(label)
        new_handles.append(tuple([handles[i] for i in legend_dict[label]]))
    # create the new legend
    ax.legend(new_handles, new_labels, loc='best', fontsize=12)

plt.show()
