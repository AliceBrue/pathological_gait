"""
This script contains various methods used to extract and plot metrics of SCONE
simulation results.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import extract_sto
import seaborn as sns

sns.set_theme()
sns.set_context("paper", font_scale=1.5)
sns.set_style("whitegrid")

results_folder = '../results'


def hex_to_rgb(value):
    '''Converts hex to rgb colours
    Parameters
    ----------
    value: (string) 6 characters representing a hex colour.

    Returns
    ----------
    list: length 3 of RGB values
    '''
    value = value.strip("#")  # removes hash symbol if present
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))


def rgb_to_dec(value):
    '''Converts rgb to decimal colours (i.e. divides each value by 256)
    Parameters
    ----------
    value: (list) list length 3 of RGB values

    Returns
    ----------
    list: length 3 of decimal values
    '''
    return [v / 256 for v in value]


def get_continuous_cmap(hex_list, float_list=None):
    '''Creates and returns a color map that can be used in heat map figures.
    If float_list is not provided, colour map graduates linearly between each color in hex_list.
    If float_list is provided, each color in hex_list is mapped to the respective location in float_list.

    Parameters
    ----------
    hex_list: (list) list of hex code strings
    float_list: (list) list of floats between 0 and 1, same length as hex_list. Must start with 0 and end with 1.

    Returns
    ----------
    colour map
    '''
    rgb_list = [rgb_to_dec(hex_to_rgb(i)) for i in hex_list]
    if float_list:
        pass
    else:
        float_list = list(np.linspace(0, 1, len(rgb_list)))

    cdict = dict()
    for num, col in enumerate(['red', 'green', 'blue']):
        col_list = [[float_list[i], rgb_list[i][num], rgb_list[i][num]] for i in range(len(float_list))]
        cdict[col] = col_list
    cmp = mcolors.LinearSegmentedColormap('my_cmp', segmentdata=cdict, N=256)
    return cmp


def plot_healthy_clinical_angles(ax):
    '''Plots the value ranges of ankle angle from healthy clinical observations.

    Parameters
    ----------
    ax : (Matplotlib <AxesSubplot:>) axes subplot object to plot on.

    Returns
    -------
    None.
    '''
    norm_min = [-7.7, -9.1, -10.3, -9.9, -8.2, -6.1, -4.2, -2.5,
                -1.0, 0.2, 1.4, 2.3, 3.1, 3.7, 4.3, 4.8, 5.3, 5.7, 6.1, 6.4,
                6.7, 7.0, 7.1, 6.9, 6.3, 4.8, 2.1, -2.2, -8.6, -16.2, -23.3,
                -27.6, -28.6, -26.7, -22.9, -18.7, -14.7, -11.1, -8.1, -5.6,
                -3.7, -2.4, -1.7, -1.7, -2.4, -3.4, -4.6, -5.7, -6.4, -7.0,
                -7.8]
    norm_mean = [-2.1, -3.9, -5.6, -5.6, -4.1, -2.1, -0.2, 1.5,
                 2.9, 4.1, 5.2, 6.1, 6.8, 7.4, 8.0, 8.6, 9.1, 9.7,
                 10.2, 10.7, 11.1, 11.5, 11.8, 11.8, 11.5, 10.4,
                 8.4, 4.9, -0.2, -6.9, -13.8, -18.5, -19.8, -18.0,
                 -14.7, -11.2, -7.9, -4.9, -2.4, -0.3, 1.1, 2.1,
                 2.5, 2.4, 1.8, 1.0, 0.1, -0.6, -1.1, -1.4, -2.3]
    norm_max = [3.5, 1.3, -0.8, -1.3, -0.1, 1.9, 3.8, 5.5, 6.9,
                8.0, 9.0, 9.8, 10.5, 11.1, 11.8, 12.4, 13.0, 13.7,
                14.3, 15.0, 15.5, 16.1, 16.5, 16.8, 16.7, 16.0,
                14.6, 12.1, 8.1, 2.4, -4.3, -9.5, -10.9, -9.3,
                -6.5, -3.7, -1.1, 1.3, 3.4, 4.9, 6.0, 6.6, 6.7,
                6.4, 5.9, 5.4, 4.8, 4.5, 4.3, 4.1, 3.1]
    std = np.array(norm_max) - np.array(norm_mean)
    norm_max = np.array(norm_mean) + 2 * std
    norm_min = np.array(norm_mean) - 2 * std

    time_vector = np.linspace(0, 100, len(norm_min))
    ax.fill_between(time_vector, norm_min, norm_max, facecolor='silver', alpha=0.5, label='clin.', zorder=5)


def plot_column(column_values, parameter_values, column_name, x_label, export_path, healthy_value, healthy_metric,
                inv=False, plot=False):
    '''Plots values of a column with parameter_values on x axis.
    Parameters
    ---------
    column_values: values to be plotted
    parameter_values: values on x axis
    column_name: name of the column, y axis legend
    x_label: name of the parameters
    export_path: folder path of export_file
    healthy value: value from healthy sto file of interest
    healthy metric: value of the metric for healthy optimisation

    Returns
    ---------
    None
    '''
    fig, ax = plt.subplots()
    for i in range(len(parameter_values)):
        if int(float(parameter_values[i])) == 100:
            parameter_values[i] = 100
        elif int(str(int(float(parameter_values[i])))[:3]) == 100:
            parameter_values[i] = int(str(int(parameter_values[i]))[3:])
        elif int(str(int(float(parameter_values[i])))[-3:]) == 100:
            parameter_values[i] = int(str(int(parameter_values[i]))[:-3])
    parameter_values_extended = parameter_values + [healthy_value]
    column_values_extended = column_values + [healthy_metric]
    sorted_idx = np.argsort(parameter_values_extended)
    column_values_extended = [column_values_extended[idx] for idx in sorted_idx]
    parameter_values_extended = [parameter_values_extended[idx] for idx in sorted_idx]

    ax.plot(parameter_values_extended, column_values_extended, color='blue')
    ax.plot(healthy_value, healthy_metric, color='black', marker='D')

    ax.set_xlabel(x_label + " [%]")
    ax.set_xticks(parameter_values_extended)
    if inv:
        ax.invert_xaxis()

    y_label = ""
    for s in range(len(column_name.split('_'))):
        y_label = y_label + " " + column_name.split('_')[s]
    ax.set_ylabel(y_label)
    if inv:
        ax.set_title(y_label + " for" + x_label + " from " + str(int(max(max(parameter_values), healthy_value)))
                     + " to " + str(int(min(min(parameter_values), healthy_value))) + "%")
    else:
        ax.set_title(y_label + " for" + x_label + " from " + str(int(min(min(parameter_values), healthy_value)))
                     + " to " + str(int(max(max(parameter_values), healthy_value))) + "%")

    plt.tight_layout()
    plt.savefig(os.path.join(export_path, column_name + '.png'))
    if not plot:
        plt.close(fig)


def plot_column_std_phases(column_values, std_values, st_column_values, st_std_values, sw_column_values, sw_std_values,
                           parameter_values, column_name, x_label, export_path, healthy_value, healthy_metric,
                           std_healthy, inv=False, plot=False):
    '''Plots values of a column during gait cycle, stance and swing phases with parameter_values on x axis.
    Parameters
    -------
    column_values: values to be plotted
    parameter_values: values on x axis
    column_name: name of the column, y axis legend
    x_label: name of the parameters
    export_path: folder path of export_file
    healthy value: value from healthy sto file of interest
    healthy metric: value of the metric for healthy optimisation

    Returns
    -------
    None
    '''
    fig, ax = plt.subplots()
    for i in range(len(parameter_values)):
        if int(float(parameter_values[i])) == 100:
            parameter_values[i] = 100
        elif int(str(int(float(parameter_values[i])))[:3]) == 100:
            parameter_values[i] = int(str(int(parameter_values[i]))[3:])
        elif int(str(int(parameter_values[i]))[-3:]) == 100:
            parameter_values[i] = int(str(int(float(parameter_values[i])))[:-3])
    parameter_values_extended = parameter_values + [healthy_value]
    st_column_values_extended = st_column_values + [healthy_metric]
    st_column_std_extended = st_std_values + [std_healthy]
    sw_column_values_extended = sw_column_values + [healthy_metric]
    sw_column_std_extended = sw_std_values + [std_healthy]
    column_values_extended = column_values + [healthy_metric]
    column_std_extended = std_values + [std_healthy]
    sorted_idx = np.argsort(parameter_values_extended)
    st_column_values_extended = [st_column_values_extended[idx] for idx in sorted_idx]
    sw_column_values_extended = [sw_column_values_extended[idx] for idx in sorted_idx]
    column_values_extended = [column_values_extended[idx] for idx in sorted_idx]
    parameter_values_extended = [parameter_values_extended[idx] for idx in sorted_idx]
    st_column_std_extended = [st_column_std_extended[idx] for idx in sorted_idx]
    sw_column_std_extended = [sw_column_std_extended[idx] for idx in sorted_idx]
    column_std_extended = [column_std_extended[idx] for idx in sorted_idx]

    delta = parameter_values_extended[1] - parameter_values_extended[0]
    if inv:
        ax.errorbar(np.array(parameter_values_extended) + delta / 10, column_values_extended, column_std_extended,
                    color="k", marker='D', linestyle="None", label="gait cycle")
        ax.errorbar(np.array(parameter_values_extended), st_column_values_extended, st_column_std_extended,
                    color="tab:cyan", marker='D', linestyle="None", label="stance")
        ax.errorbar(np.array(parameter_values_extended) - delta / 10, sw_column_values_extended, sw_column_std_extended,
                    color="tab:orange", marker='D', linestyle="None", label="swing")
    else:
        ax.errorbar(np.array(parameter_values_extended) - delta / 10, column_values_extended, column_std_extended,
                    color="k", marker='D', linestyle="None", label="gait cycle")
        ax.errorbar(np.array(parameter_values_extended), st_column_values_extended, st_column_std_extended,
                    color="tab:cyan", marker='D', linestyle="None", label="stance")
        ax.errorbar(np.array(parameter_values_extended) + delta / 10, sw_column_values_extended, sw_column_std_extended,
                    color="tab:orange", marker='D', linestyle="None", label="swing")

    ax.set_xlabel(x_label + " [%]")
    ax.set_xticks(parameter_values_extended)
    if inv:
        ax.invert_xaxis()

    if column_name.split("_")[0] == "ME":
        y_label = "ME " + "(" + column_name.split("_")[1] + ") [°]"
    else:
        y_label = ""
        for s in range(len(column_name.split('_')) - 1):
            y_label = y_label + " " + column_name.split('_')[s]
        y_label = y_label + " (" + column_name.split('_')[-1] + ")"
    ax.set_ylabel(y_label)
    if column_name.split("_")[0] == "ME":
        y_label = "ME " + "(" + column_name.split("_")[1] + ")"
    if inv:
        ax.set_title(y_label + " for" + x_label + " from " + str(int(max(max(parameter_values), healthy_value)))
                     + " to " + str(int(min(min(parameter_values), healthy_value))) + "%")
    else:
        ax.set_title(y_label + " for" + x_label + " from " + str(int(min(min(parameter_values), healthy_value)))
                     + " to " + str(int(max(max(parameter_values), healthy_value))) + "%")

    ax.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(export_path, column_name + '_std_phases.png'))
    if not plot:
        plt.close(fig)


def plot_column_std_two_muscles(column_values_1, std_values_1, column_values_2, std_values_2, parameter_values,
                                column_name, x_label, export_path, healthy_value, healthy_metric_1, healthy_metric_2,
                                std_healthy_1, std_healthy_2, side, inv=False, plot=False):
    '''Plots values of a column of both soleus and gastroc muscles with parameter_values on x axis.
    sto_file: (string) path to sto_file
    Parameters
    -------
    column_values: values to be plotted
    parameter_values: values on x axis
    column_name: name of the column, y axis legend
    x_label: name of the parameters
    export_path: folder path of export_file
    healthy value: value from healthy sto file of interest
    healthy metric: value of the metric for healthy optimisation

    Returns
    -------
    None
    '''
    fig, ax = plt.subplots()
    for i in range(len(parameter_values)):
        if int(float(parameter_values[i])) == 100:
            parameter_values[i] = 100
        elif int(str(int(float(parameter_values[i])))[:3]) == 100:
            parameter_values[i] = int(str(int(parameter_values[i]))[3:])
        elif int(str(int(parameter_values[i]))[-3:]) == 100:
            parameter_values[i] = int(str(int(float(parameter_values[i])))[:-3])
    parameter_values_extended = parameter_values + [healthy_value]
    column_values_extended_1 = column_values_1 + [healthy_metric_1]
    column_std_extended_1 = std_values_1 + [std_healthy_1]
    column_values_extended_2 = column_values_2 + [healthy_metric_2]
    column_std_extended_2 = std_values_2 + [std_healthy_2]
    sorted_idx = np.argsort(parameter_values_extended)
    st_column_values_extended = [column_values_extended_1[idx] for idx in sorted_idx]
    sw_column_values_extended = [column_values_extended_2[idx] for idx in sorted_idx]
    parameter_values_extended = [parameter_values_extended[idx] for idx in sorted_idx]
    st_column_std_extended = [column_std_extended_1[idx] for idx in sorted_idx]
    sw_column_std_extended = [column_std_extended_2[idx] for idx in sorted_idx]

    delta = parameter_values_extended[1] - parameter_values_extended[0]

    ax.errorbar(np.array(parameter_values_extended) - delta / 12, st_column_values_extended, st_column_std_extended,
                color="dimgrey", marker='D', linestyle="None", label="SOL")
    ax2 = ax.twinx()
    ax2.errorbar(np.array(parameter_values_extended) + delta / 12, sw_column_values_extended, sw_column_std_extended,
                 color="darkred", marker='D', linestyle="None", label="GAS")

    ax.set_xlabel(x_label + " [%]")
    ax.set_xticks(parameter_values_extended)
    if inv:
        ax.invert_xaxis()

    y_label = "spasticity index (" + side + ")"
    ax.set_ylabel("SOL " + y_label, color="dimgrey")
    ax2.set_ylabel("GAS " + y_label, color="darkred")
    if inv:
        ax.set_title(y_label + " for" + x_label + " from " + str(int(max(max(parameter_values), healthy_value)))
                     + " to " + str(int(min(min(parameter_values), healthy_value))) + "%")
    else:
        ax.set_title(y_label + " for" + x_label + " from " + str(int(min(min(parameter_values), healthy_value)))
                     + " to " + str(int(max(max(parameter_values), healthy_value))) + "%")

    lines, labels = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2)
    ax2.grid(None)
    plt.tight_layout()
    plt.savefig(os.path.join(export_path, column_name + '.png'))
    if not plot:
        plt.close(fig)


def plot_mean_gc(sto_files, parameter_values, var_list, side, ylabel, title, export_path, healthy_sto, healthy_value,
                 inv=False, plot=False):
    '''Plots averaged variables during gait cycle, for multiple sto files
    Parameters
    ----------
    sto_file: path to the sto file of interest
    var_list: list of variables to plot
    var_legend: list of legends for the previous variables
    side: 'r' or 'l' for right or left variables
    ylabel: plot y label
    title: plot title
    export_path: folder path of export file
    healthy_sto: path to healthy sto file of interest
    healthy_value: value to healthy parameter

    Returns
    -------
    None
    '''
    var_names_list = []
    var_tab_list = []
    time_list = []

    mme = []
    std_me = []
    st_mme = []
    st_std_me = []
    sw_mme = []
    sw_std_me = []

    fig, ax = plt.subplots()
    par_values = parameter_values
    for i in range(len(parameter_values)):
        if int(str(int(parameter_values[i]))) == 100:
            parameter_values[i] = 100
        elif int(str(int(parameter_values[i]))[:3]) == 100:
            parameter_values[i] = int(str(int(parameter_values[i]))[3:])
        elif int(str(int(parameter_values[i]))[-3:]) == 100:
            parameter_values[i] = int(str(int(parameter_values[i]))[:-3])

    min_value = np.min(parameter_values)
    max_value = np.max(parameter_values)
    color_offset = mcolors.Normalize(vmin=min_value, vmax=max_value)

    # start with healthy plot
    parameter_value = healthy_value
    sto_file = healthy_sto
    var_names, var_tab = extract_sto.extract_sto(sto_file)
    time = var_tab[:, 0]
    var_names_list.append(var_names)
    var_tab_list.append(var_tab)
    time_list.append(time)
    _, _, _, lend, lst = extract_sto.me_mean_gait_phases(var_names, var_tab, var_list, side)
    av_stance_var, av_swing_var, av_gait_cycle, av_stance_end = extract_sto.mean_gait_phases(var_names, var_tab,
                                                                                             var_list, side)
    av_var_h = np.concatenate((av_stance_var, av_swing_var))

    if "angle" in var_list:
        av_var_h = np.multiply(av_var_h, 180 / np.pi)
    label = str(parameter_value)
    ax.plot(time[:len(av_var_h)] * 100 / time[av_gait_cycle], av_var_h, label=label, color='black', linewidth=2,
            linestyle='dashed', zorder=10)
    ax.axvspan(time[av_stance_end] * 100 / time[av_gait_cycle], 100, alpha=0.5, color='blanchedalmond', zorder=0)
    if min_value > healthy_value:
        hex_list = ['#0000ff', '#ff0000']
    else:
        hex_list = ['#ff0000', '#0000ff']
    cmap = get_continuous_cmap(hex_list)

    for idx in range(0, len(sto_files)):
        sto_file = sto_files[idx]
        parameter_value = parameter_values[idx]
        var_names, var_tab = extract_sto.extract_sto(sto_file)
        time = var_tab[:, 0]
        var_names_list.append(var_names)
        var_tab_list.append(var_tab)
        time_list.append(time)
        av_stance_var, av_swing_var, av_gait_cycle, av_stance_end = extract_sto.mean_gait_phases(var_names, var_tab,
                                                                                                 var_list, side)
        av_var = np.concatenate((av_stance_var, av_swing_var))

        mm, std_m = extract_sto.me(sto_file, healthy_sto, var_list, side)
        st_mm, st_std_m, sw_mm, sw_std_m = extract_sto.me_phases(sto_file, healthy_sto, var_list, side)

        mme.append(mm)
        std_me.append(std_m)
        st_mme.append(st_mm)
        st_std_me.append(st_std_m)
        sw_mme.append(sw_mm)
        sw_std_me.append(sw_std_m)
        if "angle" in var_list:
            av_var = np.multiply(av_var, 180 / np.pi)
        label = str(parameter_value)
        ax.plot(time[:len(av_var)] * 100 / time[av_gait_cycle], av_var, label=label,
                color=cmap(color_offset(parameter_value)), zorder=10)

    ax.axhline(y=24, xmin=time[lst] / time[av_gait_cycle], xmax=(time[lend + 2]) / time[av_gait_cycle],
               color="tab:cyan", linewidth=4, label="MS and PS")
    ax.set_ylim((-45, 25))
    handles, labels = ax.get_legend_handles_labels()
    # sort both labels and handles by labels
    labls, handls = zip(*sorted(zip(labels[:-1], handles[:-1]), key=lambda t: float(t[0])))
    if float(labls[0]) < 100:
        labls = labls[::-1]
        handls = handls[::-1]
    Labels = [str(int(float(label))) + " %" for label in labls]

    for i in range(len(handles)):
        if i < len(handles) - 1:
            handles[i] = handls[i]
            labels[i] = Labels[i]
        else:
            handles[i] = handles[-1]
            labels[i] = labels[-1]
    ax.legend(handles, labels, loc=4)
    plot_healthy_clinical_angles(ax)

    ax.set_xlabel('Gait cycle [%]')
    ax.set_ylabel(ylabel + " [°]")
    parameter_values = np.array(labls).astype(float).astype(int)
    if inv:
        ax.set_title(ylabel + " for" + title + " from " + str(max(parameter_values))
                     + " to " + str(min(parameter_values)) + "%")
    else:
        ax.set_title(ylabel + " for" + title + " from " + str(min(parameter_values))
                     + " to " + str(max(parameter_values)) + "%")

    plt.tight_layout()
    plt.savefig(os.path.join(export_path, var_list + '.png'))
    if not plot:
        plt.close(fig)

    # ME
    if ylabel.split(" ")[1] == "ankle":
        plot_column_std_phases(mme, std_me, st_mme, st_std_me, sw_mme, sw_std_me, par_values, "ME_" + side, title,
                               export_path, healthy_value, 0, 0, inv=inv)


def get_sto_total_time(sto_file):
    '''Returns total simulation time from an sto file.

    Parameters
    ----------
    sto_file: (string) path to sto_file

    Returns
    -------
    float: simulation total time
    '''
    var_names, var_tab = extract_sto.extract_sto(sto_file)
    time = var_tab[:, 0]
    return time[-1]
