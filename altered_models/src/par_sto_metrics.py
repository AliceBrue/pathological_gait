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


def plot_healthy_clinical_angles(ax, var_name):
    '''Plots the value ranges of ankle angle from healthy clinical observations.

    Parameters
    ----------
    ax : (Matplotlib <AxesSubplot:>) axes subplot object to plot on.

    Returns
    -------
    None.
    '''
    if var_name in ['ankle_angle_r', 'ankle_angle_l']:
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
        norm_max = (np.array(norm_mean) + 2 * std)
        norm_min = (np.array(norm_mean) - 2 * std)

    elif var_name in ['knee_angle_r', 'knee_angle_l']:
        norm_min = [-0.1, 2.5, 5.3, 8.3, 10.8, 12.3, 13.0, 12.9, 12.4,
                    11.6, 10.6, 9.4, 8.2, 7.0, 5.8, 4.6, 3.5, 2.5,
                    1.6, 0.9, 0.5, 0.5, 1.1, 2.2, 3.8, 6.0, 8.7, 12.2,
                    16.4, 21.5, 27.3, 33.6, 39.9, 45.7, 50.0, 52.5,
                    53.3, 52.4, 50.1, 46.5, 42.0, 36.5, 30.3, 23.7,
                    16.8, 10.2, 4.4, 0.2, -2.0, -1.9, -0.2]
        norm_mean = [5.6, 7.9, 10.9, 14.1, 16.9, 18.6, 19.3, 19.2,
                     18.6, 17.6, 16.5, 15.2, 13.9, 12.6, 11.3, 10.1,
                     9.0, 8.0, 7.2, 6.6, 6.2, 6.3, 6.8, 7.9, 9.5,
                     11.6, 14.5, 18.1, 22.6, 27.9, 33.9, 40.2, 46.3,
                     51.7, 55.7, 58.2, 59.1, 58.7, 56.9, 53.9, 49.8,
                     44.7, 38.8, 32.4, 25.5, 18.7, 12.5, 7.6, 4.6,
                     4.0, 5.4]
        norm_max = [11.2, 13.4, 16.5, 20.0, 23.0, 24.8, 25.6, 25.5,
                    24.8, 23.7, 22.4, 21.0, 19.6, 18.2, 16.8, 15.6,
                    14.6, 13.6, 12.9, 12.3, 12.0, 12.0, 12.5, 13.5,
                    15.1, 17.3, 20.2, 24.0, 28.8, 34.3, 40.5, 46.8,
                    52.8, 57.7, 61.4, 63.8, 65.0, 65.0, 63.8, 61.3,
                    57.6, 52.9, 47.3, 41.0, 34.2, 27.2, 20.5, 14.9,
                    11.2, 10.0, 11.0]
        std = np.array(norm_max) - np.array(norm_mean)
        norm_max = -(np.array(norm_mean) + 2 * std)
        norm_min = -(np.array(norm_mean) - 2 * std)

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
    inv: bool to inverse decreasing parameter values
    plot: bool to show plot

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
    if column_name == "speed":
        ax.set_ylabel(y_label+' [m/s]')
    elif column_name == "total_time":
        ax.set_ylabel(y_label+' [s]')
    else:
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


def plot_columns_std(column_values1, std_values1, parameter_values, column_name, x_label, export_path,
                     healthy_value, healthy_metric, std_healthy, inv=False, column_values2=None, std_values2=None,
                     healthy_metric2=0, std_healthy2=0, plot=False):
    '''Plots values of a column during gait cycle, stance and swing phases with parameter_values on x axis.
    Parameters
    -------
    column_values1: values to be plotted
    std_values1: values std to be plotted
    parameter_values: values on x axis
    column_name: name of the column, y axis legend
    x_label: name of the parameters
    export_path: folder path of export_file
    healthy value: value from healthy sto file of interest
    healthy metric: value of the metric for healthy optimisation
    std_healthy: std value of the metric for healthy optimisation
    inv: bool to inverse decreasing parameter values
    column_values2: values to be plotted
    std_values2: values std to be plotted
    healthy metric2: value of the metric for healthy optimisation
    std_healthy2: std value of the metric for healthy optimisation
    plot: bool to show plot

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
    column_values_extended1 = column_values1 + [healthy_metric]
    column_std_extended1 = std_values1 + [std_healthy]
    if column_values2 is not None:
        column_values_extended2 = column_values2 + [healthy_metric2]
        column_std_extended2 = std_values2 + [std_healthy2]

    sorted_idx = np.argsort(parameter_values_extended)
    parameter_values_extended = [parameter_values_extended[idx] for idx in sorted_idx]
    column_values_extended1 = [column_values_extended1[idx] for idx in sorted_idx]
    column_std_extended1 = [column_std_extended1[idx] for idx in sorted_idx]
    if column_values2 is not None:
        column_std_extended2 = [column_std_extended2[idx] for idx in sorted_idx]
        column_std_extended2 = [column_std_extended2[idx] for idx in sorted_idx]

    if column_name.split("_")[0] == "ME" and column_values2 is None:
        label1 = 'ME ST'
        color1 = 'tab:cyan'
    elif column_name.split("_")[0] == "ME" and column_values2 is not None:
        label1 = 'min ES'
        color1 = 'darkorange'
        label2 = 'ME ST'
        color2 = 'tab:cyan'
    elif column_name.split("_")[0] == "mstance":
        label1 = 'stance T'
        color1 = 'tab:blue'
        label2 = 'step L'
        color2 = 'dimgrey'
    elif column_name.split("_")[0] == "spasticity":
        label1 = 'SOL'
        color1 = 'dimgrey'
        label2 = 'GAS'
        color2 = 'darkred'
    delta = parameter_values_extended[1] - parameter_values_extended[0]
    if inv:
        ax.errorbar(np.array(parameter_values_extended) + delta/12, column_values_extended1, column_std_extended1,
                    color=color1, marker='D', linestyle="None", label=label1)
        if column_values2 is not None:
            ax2 = ax.twinx()
            ax2.errorbar(np.array(parameter_values_extended) - delta/12, column_values_extended2, column_std_extended2,
                        color=color2, marker='D', linestyle="None", label=label2)
    else:
        ax.errorbar(np.array(parameter_values_extended) - delta/12, column_values_extended1, column_std_extended1,
                    color=color1, marker='D', linestyle="None", label=label1)
        if column_values2 is not None:
            ax2 = ax.twinx()
            ax2.errorbar(np.array(parameter_values_extended) + delta/12, column_values_extended2, column_std_extended2,
                        color=color2, marker='D', linestyle="None", label=label2)

    ax.set_xlabel(x_label + " [%]")
    ax.set_xticks(parameter_values_extended)
    if inv:
        ax.invert_xaxis()

    if column_name.split("_")[0] == "ME":
        title = 'Ankle comparison'
        if column_values2 is not None:
            y_label1 = "min angle ES " + "(" + column_name.split("_")[1] + ") [°]"
            y_label2 = "ME ST " + "(" + column_name.split("_")[1] + ") [°]"
        else:
            y_label1 = "ME ST " + "(" + column_name.split("_")[1] + ") [°]"
    elif column_name.split("_")[0] == "mstance":
        title = 'Gait features'
        y_label1 = "stance period (T) [s]"
        y_label2 = "step length (L) [m]"
    elif column_name.split("_")[0] == "spasticity":
        title = 'Spasticity indexes'
        y_label1 = "SOL spasticity index (" + column_name.split("_")[1] + ")"
        y_label2 = "GAS spasticity index (" + column_name.split("_")[1] + ")"
    ax.set_ylabel(y_label1, color=color1)
    if column_values2 is not None:
        ax2.set_ylabel(y_label2, color=color2)
        ax2.grid(None)
    if inv:
        ax.set_title(title + " for" + x_label + " from " + str(int(max(max(parameter_values), healthy_value)))
                     + " to " + str(int(min(min(parameter_values), healthy_value))) + "%")
    else:
        ax.set_title(title + " for" + x_label + " from " + str(int(min(min(parameter_values), healthy_value)))
                     + " to " + str(int(max(max(parameter_values), healthy_value))) + "%")

    if column_values2 is not None:
        lines, labels = ax.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax2.legend(lines + lines2, labels + labels2)
    plt.tight_layout()
    plt.savefig(os.path.join(export_path, column_name + '_std.png'))
    if not plot:
        plt.close(fig)


def plot_mean_gc(sto_files, parameter_values, var_name, side, title, export_path, healthy_sto, healthy_value,
                 inv=False, es='es', plot=False):
    '''Plots averaged variables during gait cycle, for multiple sto files.
    Parameters
    ----------
    sto_files: list of path to the sto file of interest
    parameter_values: values on x axis
    var_name: list of variables to plot
    side: 'r' or 'l' for right or left variables
    title: plot title
    export_path: folder path of export file
    healthy_sto: path to healthy sto file of interest
    healthy_value: value to healthy parameter
    inv: bool to inverse decreasing parameter values
    es: ankle measures to plot
    plot: bool to show plot

    Returns
    -------
    None
    '''
    var_names_list = []
    var_tab_list = []
    time_list = []

    st_mme = []
    st_std_me = []
    es_mmin = []
    es_std_min = []

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
    _, _, _, lend, lst = extract_sto.me_mean_gait_phases(var_names, var_tab, var_name, side)
    h_av_stance_var, h_av_swing_var, h_av_gait_cycle, h_av_stance_end = extract_sto.mean_gait_phases(var_names, var_tab,
                                                                                             var_name, side)
    av_var_h = np.concatenate((h_av_stance_var, h_av_swing_var))

    if "angle" in var_name:
        av_var_h = np.multiply(av_var_h, 180 / np.pi)
    label = str(parameter_value)
    ax.plot(time[:len(av_var_h)] * 100 / time[h_av_gait_cycle], av_var_h, label=label, color='black', linewidth=2,
            linestyle='dashed', zorder=10)
    ax.axvspan(time[h_av_stance_end] * 100 / time[h_av_gait_cycle], 100, alpha=0.5, color='blanchedalmond', zorder=0)

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
                                                                                                 var_name, side)
        av_var = np.concatenate((av_stance_var, av_swing_var))

        st_mm, st_std_m, sw_mm, sw_std_m = extract_sto.me_phases(sto_file, healthy_sto, var_name, side)
        st_mme.append(st_mm)
        st_std_me.append(st_std_m)

        es_min, es_smin = extract_sto.min_ankle_es(sto_file, var_name, side)
        es_mmin.append(es_min)
        es_std_min.append(es_smin)

        if "angle" in var_name:
            av_var = np.multiply(av_var, 180 / np.pi)
        label = str(parameter_value)
        ax.plot(time[:len(av_var)] * 100 / time[av_gait_cycle], av_var, label=label,
                color=cmap(color_offset(parameter_value)), zorder=10)

    if es == "es":
        ax.axhline(y=24, xmin=0, xmax=time[lst - 1] / time[h_av_gait_cycle], color="darkorange", linewidth=4,
                   label="ES")
    if "ankle" in var_name:
        xlabel = title
        ylabel = 'ankle angle'
        title = "Ankle angle"
        ax.axhline(y=24, xmin=time[lst] / time[h_av_gait_cycle], xmax=(time[h_av_stance_end]) / time[h_av_gait_cycle],
                   color="tab:cyan", linewidth=4, label="ST")
        ax.set_ylim((-45, 25))
    if "knee" in var_name:
        xlabel = title
        ylabel = 'knee angle'
        title = "Knee angle"
        ax.axhline(y=9, xmin=time[lst] / time[h_av_gait_cycle], xmax=(time[h_av_stance_end]) / time[h_av_gait_cycle],
                   color="tab:cyan", linewidth=4, label="ST")
        ax.set_ylim((-85, 10))

    handles, labels = ax.get_legend_handles_labels()
    # sort both labels and handles by labels
    if es == 'es':
        labls, handls = zip(*sorted(zip(labels[:-2], handles[:-2]), key=lambda t: float(t[0])))
        last = 2
    else:
        labls, handls = zip(*sorted(zip(labels[:-1], handles[:-1]), key=lambda t: float(t[0])))
        last = 1
    if float(labls[0]) < 100:
        labls = labls[::-1]
        handls = handls[::-1]
    Labels = [str(int(float(label))) + " %" for label in labls]

    for i in range(len(handles)):
        if i < len(handles) - last:
            handles[i] = handls[i]
            labels[i] = Labels[i]
        elif es=='es' and i == len(handles)-1:
            handles[i] = handles[-1]
            labels[i] = labels[-1]
        elif i == len(handles) - last:
            handles[i] = handles[-last]
            labels[i] = labels[-last]
    ax.legend(handles, labels)
    plot_healthy_clinical_angles(ax, var_name)

    ax.set_xlabel('gait cycle [%]')
    ax.set_xlim((0, 100))
    ax.set_ylabel(ylabel + " [°]")
    parameter_values = np.array(labls).astype(float).astype(int)
    if inv:
        ax.set_title(title + " for" + xlabel + " from " + str(max(parameter_values))
                     + " to " + str(min(parameter_values)) + "%")
    else:
        ax.set_title(title + " for" + xlabel + " from " + str(min(parameter_values))
                     + " to " + str(max(parameter_values)) + "%")

    plt.tight_layout()
    plt.savefig(os.path.join(export_path, var_name + '.png'))
    if not plot:
        plt.close(fig)

    # ME
    if ylabel.split(" ")[0] == "ankle":
        if es == 'es':
            h_es_min, h_es_smin = extract_sto.min_ankle_es(healthy_sto, var_name, side)
            plot_columns_std(es_mmin, es_std_min, par_values, "ME_" + side, xlabel, export_path, healthy_value, h_es_min,
                             h_es_smin, column_values2=st_mme, std_values2=st_std_me, inv=inv)
        else:
            plot_columns_std(st_mme, st_std_me, par_values, "ME_" + side, xlabel, export_path, healthy_value, 0, 0,
                            inv=inv)


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
