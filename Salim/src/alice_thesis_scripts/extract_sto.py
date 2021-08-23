import numpy as np
import matplotlib.pyplot as plt


def extract_sto(sto_file):
    """
    Extracts variables from sto file:
    INPUTS: - sto_file: path to the sto file of interest
    OUTPUTS: - var_names: list of sto file variable names
             - var_tab: tab of previous variables during simulation
    """

    file = open(sto_file, "r")

    file.readline()
    file.readline()
    nrows = int(file.readline().split('=')[1])
    ncol = int(file.readline().split('=')[1])
    file.readline()
    file.readline()

    var_names = np.asarray(file.readline().split('\t'))

    var_tab = np.zeros((nrows, ncol))
    line = file.readline()
    i = 0
    while line:
        var_tab[i, :] = np.asarray(line.split('\t')).astype(np.float)
        line = file.readline()
        i += 1

    return var_names, var_tab


def stance_phase_indexes(var_names, var_tab):
    """
    Extracts left and right stance phase start and end indexes
    INPUTS: - var_names: list of variable names
            - var_tab: tab of previous variables during simulation
    OUTPUTS: - lstance_starts: list of left stance phase start indexes
             - lstance_ends: list of left stance phase end indexes
             - rstance_starts: list of right stance phase start indexes
             - rstance_ends: list of right stance phase end indexes
    """

    lstate_index = np.where(var_names == 'leg0_l.state')[0][0]
    rstate_index = np.where(var_names == 'leg1_r.state')[0][0]
    lstate = var_tab[:, lstate_index]
    rstate = var_tab[:, rstate_index]

    lstance_starts = np.where(lstate == 0)[0]
    lstance_starts = lstance_starts[np.concatenate(
        ([0], np.where(lstance_starts[1:] - lstance_starts[:len(lstance_starts) - 1] > 1)[0] + 1))]
    lstance_ends = np.where(lstate == 3)[0]
    lstance_ends = lstance_ends[np.concatenate(
        ([0], np.where(lstance_ends[1:] - lstance_ends[:len(lstance_ends) - 1] > 1)[0])) + 1] - 1

    rstance_starts = np.where(rstate == 0)[0]
    rstance_starts = rstance_starts[np.concatenate(
        ([0], np.where(rstance_starts[1:] - rstance_starts[:len(rstance_starts) - 1] > 1)[0] + 1))]
    rstance_ends = np.where(rstate == 3)[0]
    rstance_ends = rstance_ends[np.concatenate(
        ([0], np.where(rstance_ends[1:] - rstance_ends[:len(rstance_ends) - 1] > 1)[0] + 1))]

    return lstance_starts, lstance_ends, rstance_starts, rstance_ends


def mean_gait_phases(var_names, var_tab, var_name, side):
    """
    Interpolates and averages a variable during stance and swing phases over gait cycles
    INPUTS: - var_names: list of variable names
            - var_tab: tab of previous variables during simulation
            - var_name: name of the variable to average
            - side: 'r' or 'l' for right or left variable
    OUTPUTS: - av_stance_var: averaged variable during stance phase over gait cycles
             - av_swing_var: averaged variable during swing phase over gait cycles
             - av_gait_cycle: averaged gait cycles duration
             - av_stance_end: averaged stance phases duration
    """

    lstance_starts, lstance_ends, rstance_starts, rstance_ends = stance_phase_indexes(var_names, var_tab)

    var_index = np.where(var_names == var_name)[0][0]
    var = var_tab[:, var_index]

    if side == 'l':
        stance_starts = lstance_starts
        stance_ends = lstance_ends
    else:
        stance_starts = rstance_starts
        stance_ends = rstance_ends

    av_gait_cycle = int(np.mean(stance_starts[1:] - stance_starts[:len(stance_starts) - 1]))
    av_cycle_var = np.zeros((len(stance_starts) - 1, av_gait_cycle))

    # interpolation of each gait cycle on the averaged gait cycles duration
    for t in range(len(stance_starts) - 1):
        cycle_var = var[stance_starts[t]:stance_starts[t + 1]]
        av_cycle_var[t, :] = np.interp(np.arange(av_gait_cycle), np.arange(len(cycle_var)), cycle_var)
    av_cycle_var0 = np.mean(av_cycle_var, axis=0)

    if stance_starts[0] < stance_ends[0]:
        av_stance_end = int(np.mean(stance_ends[:min(len(stance_starts), len(stance_ends))] -
                                    stance_starts[:min(len(stance_starts), len(stance_ends))]))
    else:
        av_stance_end = int(np.mean(stance_ends[1:len(stance_ends)] -
                                    stance_starts[:len(stance_ends) - 1]))

    # average of the variable of interest during stance and swing phases over gait cycles
    av_stance_var = av_cycle_var0[0:av_stance_end]
    av_swing_var = av_cycle_var0[av_stance_end:av_gait_cycle]

    return av_stance_var, av_swing_var, av_gait_cycle, av_stance_end


def plot_mean_gc(sto_file, var_list, var_legend, side, ylabel, title, sto_file2=None):
    """
    Plots averaged variables during gait cycle, two simulations can be plotted together
    INPUTS: - sto_file: path to the sto file of interest
            - var_list: list of variables to plot
            - var_legend: list of legends for the previous variables
            - side: 'r' or 'l' for right or left variables
            - ylabel: plot y label
            - title: plot title
            - sto_file2: optional path to second sto file of interest
    OUTPUTS: - mean_gc_var: mean of the variable over gait cycle
             - std_gc_var: std of the variable over gait cycle
             - mean_stance_var: mean of the variable over stance phase
             - std_stance_var: std of the variable over stance phase
    """

    var_names, var_tab = extract_sto(sto_file)
    time = var_tab[:, 0]

    if sto_file2 is not None:
        var_names2, var_tab2 = extract_sto(sto_file2)

    fig, ax = plt.subplots()
    colors = ['b', 'g', 'c', 'r']
    c = 0

    for var in range(len(var_list)):
        av_stance_var, av_swing_var, av_gait_cycle, av_stance_end = mean_gait_phases(var_names, var_tab, var_list[var],
                                                                                     side)
        av_var = np.concatenate((av_stance_var, av_swing_var, av_stance_var))
        ax.plot(time[:len(av_var)]*100/time[av_gait_cycle], av_var, color=colors[c], label=var_legend[var])
        if sto_file2 is not None:
            av_stance_var2, av_swing_var2, av_gait_cycle2, av_stance_end2 = mean_gait_phases(var_names2, var_tab2,
                                                                                                 var_list[var], side)
            av_var2 = np.concatenate((av_stance_var2, av_swing_var2, av_stance_var2))
            ax.plot(time[:len(av_var2)]*100/time[av_gait_cycle2], av_var2, color=colors[c], linestyle='dashed',
                    label='Healthy '+var_legend[var])
        c += 1
    if sto_file2 is not None:
        ax.axvspan(time[av_stance_end2]*100/time[av_gait_cycle2], time[av_gait_cycle2]*100/time[av_gait_cycle2],
                   alpha=0.5, color='bisque', label='Healthy right swing')
    else:
        ax.axvspan(time[av_stance_end] * 100 / time[av_gait_cycle], time[av_gait_cycle] * 100 / time[av_gait_cycle],
                   alpha=0.5, color='bisque', label='Right swing')

    ax.legend()
    ax.set_xlabel('Gait cycle (%)')
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    plt.show()

    av_gc_var = np.concatenate((av_stance_var, av_swing_var))
    mean_gc_var = np.mean(av_gc_var)
    std_gc_var = np.std(av_gc_var)
    mean_stance_var = np.mean(av_stance_var)
    std_stance_var = np.std(av_stance_var)

    if sto_file2 is not None:
        av_gc_var2 = np.concatenate((av_stance_var2, av_swing_var2))
        mean_gc_var2 = np.mean(av_gc_var2)
        std_gc_var2 = np.std(av_gc_var2)
        mean_stance_var2 = np.mean(av_stance_var2)
        std_stance_var2 = np.std(av_stance_var2)
        return mean_gc_var, std_gc_var, mean_stance_var, std_stance_var, \
               mean_gc_var2, std_gc_var2, mean_stance_var2, std_stance_var2
    else:
        return mean_gc_var, std_gc_var, mean_stance_var, std_stance_var,


def get_mean_gc(sto_file, var_list, side, sto_file2=None):
    """
    Gets averaged variables during gait cycle, two simulations can be returned together
    INPUTS: - sto_file: path to the sto file of interest
            - var_list: list of variables to return
            - side: 'r' or 'l' for right or left variables
            - sto_file2: optional path to second sto file of interest
    OUTPUTS: - mean_gc_var: mean of the variable over gait cycle
             - std_gc_var: std of the variable over gait cycle
             - mean_stance_var: mean of the variable over stance phase
             - std_stance_var: std of the variable over stance phase
    """

    var_names, var_tab = extract_sto(sto_file)
    time = var_tab[:, 0]

    if sto_file2 is not None:
        var_names2, var_tab2 = extract_sto(sto_file2)
        
    c = 0

    for var in range(len(var_list)):
        av_stance_var, av_swing_var, av_gait_cycle, av_stance_end = mean_gait_phases(var_names, var_tab, var_list[var],
                                                                                     side)
        if sto_file2 is not None:
            av_stance_var2, av_swing_var2, av_gait_cycle2, av_stance_end2 = mean_gait_phases(var_names2, var_tab2,
                                                                                                 var_list[var], side)
            av_var2 = np.concatenate((av_stance_var2, av_swing_var2, av_stance_var2))
            
        c += 1

    av_gc_var = np.concatenate((av_stance_var, av_swing_var))
    mean_gc_var = np.mean(av_gc_var)
    std_gc_var = np.std(av_gc_var)
    mean_stance_var = np.mean(av_stance_var)
    std_stance_var = np.std(av_stance_var)
    mean_swing_var = np.mean(av_swing_var)
    std_swing_var = np.std(av_swing_var)

    if sto_file2 is not None:
        av_gc_var2 = np.concatenate((av_stance_var2, av_swing_var2))
        mean_gc_var2 = np.mean(av_gc_var2)
        std_gc_var2 = np.std(av_gc_var2)
        mean_stance_var2 = np.mean(av_stance_var2)
        std_stance_var2 = np.std(av_stance_var2)
        mean_swing_var2 = np.mean(av_swing_var2)
        std_swing_var2 = np.std(av_swing_var2)
        return mean_gc_var, std_gc_var, mean_stance_var, std_stance_var, mean_swing_var, std_swing_var, \
               mean_gc_var2, std_gc_var2, mean_stance_var2, std_stance_var2, mean_swing_var2, std_swing_var2
    else:
        return mean_gc_var, std_gc_var, mean_stance_var, std_stance_var, mean_swing_var, std_swing_var