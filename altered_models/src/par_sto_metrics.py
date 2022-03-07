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
import opensim
from scipy.signal import find_peaks

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
    norm_mean = None
    norm_max = None
    if "hip_flexion" in var_name:
        norm_min = [30.4, 30.2, 29.9, 29.3, 28.3, 26.9, 25.1, 23.1,
                    20.9, 18.7, 16.5, 14.3, 12.1, 9.9, 7.6, 5.4, 3.2,
                    1.0, -1.1, -3.2, -5.1, -6.9, -8.5, -10.0, -11.2,
                    -12.1, -12.6, -12.5, -11.7, -10.1, -7.7, -4.4,
                    -0.5, 3.6, 7.8, 11.7, 15.4, 18.7, 21.7, 24.3,
                    26.5, 28.3, 29.7, 30.6, 31.1, 31.2, 31.0, 30.6,
                    30.3, 30.1, 29.9 ]
        norm_mean = [35.8, 35.7, 35.5, 35.1, 34.2, 32.8, 31.1, 29.1,
                     26.9, 24.7, 22.5, 20.3, 18.1, 15.8, 13.6, 11.4,
                     9.3, 7.3, 5.3, 3.4, 1.6, -0.1, -1.7, -3.2, -4.4,
                     -5.3, -5.8, -5.7, -4.9, -3.4, -0.9, 2.3, 6.1,
                     10.1, 14.1, 17.9, 21.4, 24.6, 27.4, 29.9, 32.0,
                     33.7, 35.1, 36.0, 36.6, 36.7, 36.5, 36.2, 35.8,
                     35.6, 35.5 ]
        norm_max = [41.3, 41.2, 41.1, 40.9, 40.2, 38.8, 37.1, 35.0,
                    32.9, 30.7, 28.5, 26.3, 24.0, 21.8, 19.6, 17.4,
                    15.4, 13.6, 11.8, 10.0, 8.3, 6.7, 5.1, 3.7, 2.5,
                    1.5, 1.0, 1.1, 1.8, 3.4, 5.8, 9.0, 12.7, 16.6,
                    20.4, 24.0, 27.4, 30.5, 33.2, 35.5, 37.5, 39.1,
                    40.5, 41.4, 42.0, 42.2, 42.1, 41.7, 41.3, 41.0,
                    41.0 ]
    elif var_name in ['ankle_angle_r', 'ankle_angle_l']:
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

    elif var_name in ['ankle_moment_l', 'ankle_moment_r']:
        norm_min = [-0.01, -0.14, -0.23, -0.22, -0.19, -0.14, -0.10,
                    -0.06, -0.02, 0.02, 0.08, 0.14, 0.20, 0.27, 0.34, 0.41, 0.48,
                    0.55, 0.63, 0.71, 0.79, 0.87, 0.93, 0.98, 0.99, 0.95, 0.83,
                    0.61, 0.31, 0.05, -0.06, -0.07, -0.06, -0.04, -0.03, -0.02,
                    -0.01, -0.01, -0.01, -0.01, -0.01, -0.02, -0.02, -0.02, -0.02,
                    -0.02, -0.02, -0.01, -0.01, 0.00, -0.01]
        norm_mean = [0.00, -0.08, -0.13, -0.11, -0.06, 0.00, 0.06,
                     0.12, 0.18, 0.25, 0.31, 0.38, 0.44, 0.50, 0.56,
                     0.62, 0.69, 0.76, 0.84, 0.92, 1.00, 1.08, 1.16,
                     1.21, 1.23, 1.19, 1.07, 0.86, 0.56, 0.25, 0.05,
                     -0.02, -0.03, -0.03, -0.02, -0.01, -0.01, -0.01,
                     -0.01, -0.01, -0.01, -0.01, -0.01, -0.02, -0.02,
                     -0.01, -0.01, 0.00, 0.00, 0.00, 0.00]
        norm_max = [0.01, -0.02, -0.03, 0.00, 0.07, 0.15, 0.23, 0.31,
                    0.39, 0.47, 0.55, 0.61, 0.67, 0.72, 0.78, 0.84,
                    0.90, 0.97, 1.04, 1.12, 1.21, 1.30, 1.39, 1.45,
                    1.47, 1.42, 1.31, 1.11, 0.81, 0.46, 0.17, 0.03,
                    -0.01, -0.01, -0.01, -0.01, -0.01, -0.01, -0.01,
                    -0.01, -0.01, -0.01, -0.01, -0.01, -0.01, -0.01,
                    0.00, 0.00, 0.01, 0.01, 0.01]
        std = np.array(norm_max) - np.array(norm_mean)
        print("Ankle std: ", np.mean(std), "° , during stance: ", np.mean(std[:int(0.6 * len(std))]), "°")

    if norm_mean is not None and norm_max is not None:
        std = np.array(norm_max) - np.array(norm_mean)
        norm_max = (np.array(norm_mean) + 2 * std)
        norm_min = (np.array(norm_mean) - 2 * std)
        time_vector = np.linspace(0, 100, len(norm_min))
        ax.fill_between(time_vector, norm_min, norm_max, facecolor='silver', alpha=0.5, label='clin.', zorder=5)


def plot_metrics_std3(metric_values1, std_values1, metric_values3, std_values3, parameter_values, metric_name, x_label,
                      export_path, healthy_value, healthy_metric1, std_healthy1, healthy_metric3, std_healthy3,
                      es=None, inv=False, metric_values2=None, std_values2=None, healthy_metric2=0, std_healthy2=0, plot=False):
    '''Plots 3 values with std for various parameter_values on x axis.
    Parameters
    -------
    metric_values1: values to be plotted
    std_values1: values std to be plotted
    parameter_values: values on x axis
    metric_name: name of the metric, y axis legend
    x_label: name of the parameters
    export_path: folder path of export_file
    healthy value: value from healthy sto file of interest
    healthy metric: value of the metric for healthy optimisation
    std_healthy: std value of the metric for healthy optimisation
    inv: bool to inverse decreasing parameter values
    metric_values2: values to be plotted
    std_values2: values std to be plotted
    healthy metric2: value of the metric for healthy optimisation
    std_healthy2: std value of the metric for healthy optimisation
    plot: bool to show plot

    Returns
    -------
    None
    '''
    fig, ax = plt.subplots()
    parameter_values_extended = parameter_values + [healthy_value]
    metric_values_extended1 = metric_values1 + [healthy_metric1]
    metric_std_extended1 = std_values1 + [std_healthy1]
    metric_values_extended3 = metric_values3 + [healthy_metric3]
    metric_std_extended3 = std_values3 + [std_healthy3]
    if metric_values2 is not None:
        metric_values_extended2 = metric_values2 + [healthy_metric2]
        metric_std_extended2 = std_values2 + [std_healthy2]
    sorted_idx = np.argsort(parameter_values_extended)
    parameter_values_extended = [parameter_values_extended[idx] for idx in sorted_idx]
    metric_values_extended1 = [metric_values_extended1[idx] for idx in sorted_idx]
    metric_std_extended1 = [metric_std_extended1[idx] for idx in sorted_idx]
    metric_values_extended3 = [metric_values_extended3[idx] for idx in sorted_idx]
    metric_std_extended3 = [metric_std_extended3[idx] for idx in sorted_idx]
    if metric_values2 is not None:
        metric_values_extended2 = [metric_values_extended2[idx] for idx in sorted_idx]
        metric_std_extended2 = [metric_std_extended2[idx] for idx in sorted_idx]
    if es != 'es':
        metric_values_extended3 = metric_values_extended3 - healthy_metric3

    if metric_name.split("_")[0] == "ME" and metric_values2 is None:
        #label1 = 'ME ST'
        label1 = 'mean ST'
        label2 = 'max ST'
        color1 = 'tab:cyan'
        color3 = 'darkorange'
    elif metric_name.split("_")[0] == "ME" and metric_values2 is not None:
        """label1 = 'min ES'
        color1 = 'darkorange'
        label2 = 'ME ST'
        color2 = 'tab:cyan'"""
        label1 = r'$\Delta_{h}$'+' mean ST' #r'$\Delta \bar{\alpha}^{ST}$'
        label2 = 'max ST'
        color1 = 'tab:cyan'
        color2 = 'tab:cyan'
        """if 'stance' in x_label:
            label3 = r'$\Delta$' + ' moment peaks'
        else:
            label3 = 'max moment'"""
        color3 = 'darkorange'
    elif metric_name.split("_")[0] == "mstance":
        label1 = 'stance T'
        color1 = 'tab:blue'
        label2 = 'step L'
        color2 = 'dimgrey'
    delta = parameter_values_extended[1] - parameter_values_extended[0]
    if inv:
        ax.errorbar(np.array(parameter_values_extended) + delta / 20, metric_values_extended1-healthy_metric1, metric_std_extended1,
                    color=color1, marker='D', linestyle="None", label=label1)
        if metric_values2 is not None:
            ax.errorbar(np.array(parameter_values_extended) + delta / 20, metric_values_extended2,
                     metric_std_extended2, color=color2, marker='X', linestyle="None", label=label2, ms=8)
        ax2 = ax.twinx()
        ax2.errorbar(np.array(parameter_values_extended) - delta / 20, metric_values_extended3, metric_std_extended3,
                    color=color3, marker='D', linestyle="None")

    else:
        ax.errorbar(np.array(parameter_values_extended) - delta / 20, metric_values_extended1-healthy_metric1, metric_std_extended1,
                    color=color1, marker='D', linestyle="None", label=label1)
        if metric_values2 is not None:
            ax.errorbar(np.array(parameter_values_extended) - delta / 20, metric_values_extended2,
                         metric_std_extended2, color=color2, marker='X', linestyle="None", label=label2, ms=8)
        ax2 = ax.twinx()
        ax2.errorbar(np.array(parameter_values_extended) + delta / 20, metric_values_extended3, metric_std_extended3,
                    color=color3, marker='D', linestyle="None")

    ax.set_xlabel(x_label + " [%]")
    ax.set_xticks(parameter_values_extended)
    if inv:
        ax.invert_xaxis()

    if metric_name.split("_")[0] == "ME":
        title = 'Ankle comparison'
        #if metric_values2 is not None:
        #y_label1 = "min angle ES " + "(" + metric_name.split("_")[-1] + ") [°]"
        #y_label2 = "mean error (ME) ST " + "(" + metric_name.split("_")[-1] + ") [°]"
        if es == 'es':
            #y_label1 = r'$\Delta \bar{\alpha}^{ST}$' + " (" + metric_name.split("_")[-1] + ") [°]"
            y_label1 = r'$\Delta_{h}$'+' mean angle ST' + " (" + metric_name.split("_")[-1] + ") [°]"
            y_label2 = r'$\Delta$' + " moment peaks " + "(" + metric_name.split("_")[-1] + ") [Nm/kg]"
            ax2.grid(None)
        else:
            y_label1 = "ankle plantarflexion " + "(" + metric_name.split("_")[-1] + ") [°]"
            y_label2 = r'$\Delta_{h}$' +" mean moment " + "(" + metric_name.split("_")[-1] + ") [Nm/kg]"

        ax2.set_ylabel(y_label2, color=color3)
        #ax2.set_yticks(ax2.get_yticks(), np.round(ax2.get_yticks(), 2))
    elif metric_name.split("_")[0] == "mstance":
        title = 'Gait features'
        y_label1 = "stance period (T) [s]"
        y_label2 = "step length (L) [m]"
    ax.set_ylabel(y_label1, color=color1)

    if metric_values2 is not None:
        ax2.set_ylabel(y_label2, color=color3)
        ax2.grid(None)
    if inv:
        ax.set_title(title + " for" + x_label + " from " + str(int(max(max(parameter_values), healthy_value)))
                     + " to " + str(int(min(min(parameter_values), healthy_value))) + "%")
    else:
        ax.set_title(title + " for" + x_label + " from " + str(int(min(min(parameter_values), healthy_value)))
                     + " to " + str(int(max(max(parameter_values), healthy_value))) + "%")

    if metric_values2 is not None:
        lines, labels = ax.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax2.legend(lines + lines2, labels + labels2, ncol=2, prop={'size': 12}, loc=2)
    plt.tight_layout()
    plt.savefig(os.path.join(export_path, metric_name + '_std.png'))
    if not plot:
        plt.close(fig)


def plot_metric(metric_values, parameter_values, metric_name, x_label, export_path, healthy_value, healthy_metric,
                inv=False, plot=False):
    '''Plots values for various parameter_values on x axis.
    Parameters
    ---------
    metric_values: values to be plotted
    parameter_values: values on x axis
    metric_name: name of the metric, y axis legend
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
    metric_values_extended = metric_values + [healthy_metric]
    sorted_idx = np.argsort(parameter_values_extended)
    metric_values_extended = [metric_values_extended[idx] for idx in sorted_idx]
    parameter_values_extended = [parameter_values_extended[idx] for idx in sorted_idx]

    ax.plot(parameter_values_extended, metric_values_extended, color='blue')
    ax.plot(healthy_value, healthy_metric, color='black', marker='D')

    ax.set_xlabel(x_label + " [%]")
    ax.set_xticks(parameter_values_extended)
    if inv:
        ax.invert_xaxis()

    y_label = ""
    for s in range(len(metric_name.split('_'))):
        y_label = y_label + " " + metric_name.split('_')[s]
    if metric_name == "speed":
        ax.set_ylabel(y_label+' [m/s]')
    elif metric_name == "total_time":
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
    plt.savefig(os.path.join(export_path, metric_name + '.png'))
    if not plot:
        plt.close(fig)


def plot_metrics_std(metric_values1, std_values1, parameter_values, metric_name, x_label, export_path,
                     healthy_value, healthy_metric, std_healthy, inv=False, metric_values2=None, std_values2=None,
                     healthy_metric2=0, std_healthy2=0, plot=False):
    '''Plots metrics with std for various parameter_values on x axis.
    Parameters
    -------
    metric_values1: values to be plotted
    std_values1: values std to be plotted
    parameter_values: values on x axis
    metric_name: name of the metric, y axis legend
    x_label: name of the parameters
    export_path: folder path of export_file
    healthy value: value from healthy sto file of interest
    healthy metric: value of the metric for healthy optimisation
    std_healthy: std value of the metric for healthy optimisation
    inv: bool to inverse decreasing parameter values
    metric_values2: values to be plotted
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
    metric_values_extended1 = metric_values1 + [healthy_metric]
    metric_std_extended1 = std_values1 + [std_healthy]
    if metric_values2 is not None:
        metric_values_extended2 = metric_values2 + [healthy_metric2]
        metric_std_extended2 = std_values2 + [std_healthy2]

    sorted_idx = np.argsort(parameter_values_extended)
    parameter_values_extended = [int(parameter_values_extended[idx]) for idx in sorted_idx]
    metric_values_extended1 = [metric_values_extended1[idx] for idx in sorted_idx]
    metric_std_extended1 = [metric_std_extended1[idx] for idx in sorted_idx]
    if metric_values2 is not None:
        metric_values_extended2 = [metric_values_extended2[idx] for idx in sorted_idx]
        metric_std_extended2 = [metric_std_extended2[idx] for idx in sorted_idx]

    if metric_name.split("_")[0] == "ME" and metric_values2 is None:
        label1 = 'ME ST'
        color1 = 'tab:cyan'
    elif metric_name.split("_")[0] == "ME" and metric_values2 is not None:
        label1 = 'min ES'
        color1 = 'darkorange'
        label2 = 'ME ST'
        color2 = 'tab:cyan'
    elif metric_name.split("_")[0] == "mstance":
        label1 = 'stance T'
        color1 = 'tab:blue'
        label2 = 'step L'
        color2 = 'dimgrey'
    elif metric_name == "ic":
        label1 = 'MAE ankle angle'
        color1 = 'dimgrey'
        label2 = 'MRE parameters'
        color2 = 'darkred'
    delta = parameter_values_extended[1] - parameter_values_extended[0]
    if inv:
        ax.errorbar(np.array(parameter_values_extended) + delta/12, metric_values_extended1, metric_std_extended1,
                    color=color1, marker='D', linestyle="None", label=label1)
        if metric_values2 is not None:
            ax2 = ax.twinx()
            ax2.errorbar(np.array(parameter_values_extended) - delta/12, metric_values_extended2, metric_std_extended2,
                        color=color2, marker='D', linestyle="None", label=label2)
    else:
        ax.errorbar(np.array(parameter_values_extended) - delta/12, metric_values_extended1, metric_std_extended1,
                    color=color1, marker='D', linestyle="None", label=label1)
        if metric_values2 is not None:
            ax2 = ax.twinx()
            ax2.errorbar(np.array(parameter_values_extended) + delta/12, metric_values_extended2, metric_std_extended2,
                        color=color2, marker='D', linestyle="None", label=label2)

    if metric_name == "ic":
        ax.set_xlabel('initial condition (IC) number')
    else:
        ax.set_xlabel(x_label + " [%]")
    ax.set_xticks(parameter_values_extended)
    if inv:
        ax.invert_xaxis()

    if metric_name.split("_")[0] == "ME":
        title = 'Ankle comparison'
        if metric_values2 is not None:
            y_label1 = "min angle ES " + "(" + metric_name.split("_")[-1] + ") [°]"
            y_label2 = "mean error (ME) ST " + "(" + metric_name.split("_")[-1] + ") [°]"
        else:
            y_label1 = "mean error (ME) ST " + "(" + metric_name.split("_")[-1] + ") [°]"
    elif metric_name.split("_")[0] == "mstance":
        title = 'Gait features'
        y_label1 = "stance period (T) [s]"
        y_label2 = "step length (L) [m]"
    elif metric_name == "ic":
        y_label1 = 'MAE ankle angle [°]'
        y_label2 = 'MRE parameters'
    ax.set_ylabel(y_label1, color=color1)
    if metric_values2 is not None:
        ax2.set_ylabel(y_label2, color=color2)
        ax2.grid(None)
    if metric_name == 'ic':
        ax.set_title('Mean errors for ' + x_label + " % \n for various initial conditions (IC)")
    else:
        if inv:
            ax.set_title(title + " for" + x_label + " from " + str(int(max(max(parameter_values), healthy_value)))
                         + " to " + str(int(min(min(parameter_values), healthy_value))) + "%")
        else:
            ax.set_title(title + " for" + x_label + " from " + str(int(min(min(parameter_values), healthy_value)))
                     + " to " + str(int(max(max(parameter_values), healthy_value))) + "%")

    if metric_values2 is not None:
        lines, labels = ax.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax2.legend(lines + lines2, labels + labels2, ncol=2, prop={'size': 12})
    plt.tight_layout()
    plt.savefig(os.path.join(export_path, metric_name + '_std.png'))
    if not plot:
        plt.close(fig)


def plot_mean_gc(sto_files, parameter_values, var_name, side, title, export_path, healthy_sto, healthy_value,
                 inv=False, es='es', plot=False, ic=False):
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
    es: 'es' or 'no' to compute toe or heel gait metrics
    plot: bool to show plot
    ic: bool to plot various initial conditions

    Returns
    -------
    None
    '''

    # store metrics
    if 'ankle' in var_name:
        st_mean = []
        st_std_mean = []
        st_maximum = []
        st_std_maximum = []

    fig, ax = plt.subplots()
    par_values = parameter_values

    min_value = np.min(parameter_values)
    max_value = np.max(parameter_values)
    color_offset = mcolors.Normalize(vmin=min_value, vmax=max_value)

    # plot healthy variable
    parameter_value = healthy_value
    sto_file = healthy_sto
    var_names, var_tab = extract_sto.extract_sto(sto_file)
    h_time = var_tab[:, 0]
    _, _, _, lend, lst = extract_sto.me_mean_gait_phases(var_names, var_tab, var_name, side)
    h_av_stance_var, h_av_swing_var, h_av_gait_cycle, h_av_stance_end = extract_sto.mean_gait_phases(var_names, var_tab,
                                                                                                     var_name, side)
    h_av_var = np.concatenate((h_av_stance_var, h_av_swing_var))

    if "angle" in var_name or 'flexion' in var_name:
        h_av_var = np.multiply(h_av_var, 180 / np.pi)
        if es=='es':
            # compute and plot 3rd ankle rocker time
            vel_name = var_name + '_u'
            av_stance_vel, av_swing_vel, av_gait_cycle, av_stance_end = extract_sto.mean_gait_phases(var_names, var_tab,
                                                                                                     vel_name, side)
            av_vel = np.concatenate((av_stance_vel, av_swing_vel))
            if av_vel[0] < 0:
                h_third_rocker = np.where(np.diff(np.sign(av_vel)) < 0)[0][0]
            else:
                h_third_rocker = np.where(np.diff(np.sign(av_vel)) < 0)[0][1]

    if 'knee' in var_name:
        h_av_var = -h_av_var
    label = str(parameter_value)

    ax.plot(h_time[:len(h_av_var)] * 100 / h_time[h_av_gait_cycle], h_av_var, label=label, color='black', linewidth=2,
                linestyle='dashed', zorder=10)
    ax.axvspan(h_time[h_av_stance_end] * 100 / h_time[h_av_gait_cycle], 100, alpha=0.5, color='blanchedalmond',
               zorder=0)

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
        try:
            av_stance_var, av_swing_var, av_gait_cycle, av_stance_end = extract_sto.mean_gait_phases(var_names, var_tab,
                                                                                                     var_name, side)
            av_var = np.concatenate((av_stance_var, av_swing_var))
        except:
            av_var = np.nan

        if "angle" in var_name or 'flexion' in var_name:
            av_var = np.multiply(av_var, 180 / np.pi)
        if 'knee' in var_name:
            av_var =-av_var
        label = str(parameter_value)
        if not np.isnan(av_var).any():
            ax.plot(time[:len(av_var)] * 100 / time[av_gait_cycle], av_var, label=label,
                    color=cmap(color_offset(parameter_value)), zorder=10)

        if "ankle" in var_name and not np.isnan(av_var).any():
            if es == 'es':
                # compute and plot 3rd ankle rocker time
                vel_name = var_name + '_u'
                av_stance_vel, av_swing_vel, av_gait_cycle, av_stance_end = extract_sto.mean_gait_phases(var_names, var_tab,
                                                                                                         vel_name, side)
                av_vel = np.concatenate((av_stance_vel, av_swing_vel))
                if av_vel[0] < 0:
                    third_rocker = np.where(np.diff(np.sign(av_vel)) < 0)[0][0]
                else:
                    third_rocker = np.where(np.diff(np.sign(av_vel)) < 0)[0][1]
                ax.plot(time[third_rocker] * 100 / time[av_gait_cycle], av_var[third_rocker], marker='x',
                        color=cmap(color_offset(parameter_value)), zorder=10)

            # compute mean and max during stance
            st_m, st_std_m, st_max, st_std_max = extract_sto.mean_stance(sto_file, var_name, side)
            st_mean.append(st_m)
            st_std_mean.append(st_std_m)

            st_maximum.append(st_max)
            st_std_maximum.append(st_std_max)

    if 'ankle' in var_name and es == 'es':
        ax.plot(h_time[h_third_rocker] * 100 / h_time[h_av_gait_cycle], h_av_var[h_third_rocker], marker='x',
                label='3rd rocker', linestyle='None', color='black', zorder=10)

    if "ankle" in var_name:
        xlabel = title
        ylabel = 'ankle plantarflexion (' + side + ')'
        title = "Ankle angle"
        if not ic:
            ax.axhline(y=24, xmin=time[lst-1] / time[h_av_gait_cycle], xmax=(time[h_av_stance_end]) / time[h_av_gait_cycle],
                       color="tab:cyan", linewidth=4, label="ST")
        ax.set_ylim((-45, 25))
    if "knee" in var_name:
        xlabel = title
        ylabel = 'knee flexion (' + side + ')'
        title = "Knee angle"
        ax.set_ylim((-10, 85))
    if "hip" in var_name:
        xlabel = title
        ylabel = 'hip flexion (' + side + ')'
        title = "Hip angle"
        ax.set_ylim((-45, 60))

    # sort both labels and handles by labels
    handles, labels = ax.get_legend_handles_labels()

    if "ankle" in var_name:
        if es == 'es':
            labls, handls = zip(*sorted(zip(labels[:-2], handles[:-2]), key=lambda t: float(t[0])))
            last = 2
        else:
            labls, handls = zip(*sorted(zip(labels[:-1], handles[:-1]), key=lambda t: float(t[0])))
            last = 1
    else:
        labls, handls = zip(*sorted(zip(labels[:], handles[:]), key=lambda t: float(t[0])))
        last = 0

    if float(labls[0]) < 100 and not ic:
        labls = labls[::-1]
        handls = handls[::-1]
    if not ic:  # various initial conditions
        Labels = [str(int(float(label))) + " %" for label in labls]
    if ic:
        Labels = ['IC '+ str(int(float(label))) for label in labls]

    # plot all altered files
    for i in range(len(handles)):
        if i < len(handles) - last:
            handles[i] = handls[i]
            labels[i] = Labels[i]
        elif i == len(handles)-1:
            handles[i] = handles[-1]
            labels[i] = labels[-1]
        elif i == len(handles) - last:
            handles[i] = handles[-last]
            labels[i] = labels[-last]
    ax.legend(handles, labels, prop={'size': 12})
    plot_healthy_clinical_angles(ax, var_name)

    if not ic:
        ax.set_xlabel('gait cycle [%]')
    else:
        ax.set_xlabel('initial condition (IC) number')
    ax.set_xlim((0, 100))
    ax.set_ylabel(ylabel + " [°]")
    parameter_values = np.array(labls).astype(float).astype(int)
    if not ic:
        if inv:
            ax.set_title(title + " for" + xlabel + " from " + str(max(parameter_values))
                         + " to " + str(min(parameter_values)) + "%")
        else:
            ax.set_title(title + " for" + xlabel + " from " + str(min(parameter_values))
                         + " to " + str(max(parameter_values)) + "%")
    if ic:
        ax.set_title(title + " for" + xlabel + " % \n for various initial conditions (IC)")
    plt.tight_layout()
    plt.savefig(os.path.join(export_path, var_name + '.png'))
    if not plot:
        plt.close(fig)

    # ME
    if ylabel.split(" ")[0] == "ankle" and not np.isnan(av_var).any() and not ic:
        # plot and compute ankle moment metrics
        mom_name = 'ankle_moment_' + side
        mdelta_mom_peaks, sdelta_mom_peaks, mmean_moment, smean_moment, h_mdelta, h_sdelta, h_mmean, h_smean = \
            plot_mean_moment_gc(sto_files, par_values, mom_name, side, xlabel, export_path,
                                healthy_sto, healthy_value, inv=inv, plot=False)
        if es == 'es':
            # healthy metrics
            h_mme, h_sme, h_es_min, h_es_smin = extract_sto.mean_stance(healthy_sto, var_name, side)
            plot_metrics_std3(st_mean, st_std_mean, mdelta_mom_peaks, sdelta_mom_peaks, par_values, "ME_" + side, xlabel,
                              export_path, healthy_value, h_mme, h_sme, h_mdelta, h_sdelta, es=es, inv=inv)
        else:
            h_st_m, h_st_std_m, h_st_max, h_st_std_max = extract_sto.mean_stance(healthy_sto, var_name, side)
            plot_metrics_std3(st_mean, st_std_mean, mmean_moment, smean_moment, par_values, "ME_" + side, xlabel,
                              export_path, healthy_value, h_st_m, h_st_std_m, h_mmean, h_smean, metric_values2=st_maximum,
                              std_values2=st_std_maximum, healthy_metric2=h_st_max, std_healthy2=h_st_std_max, es=es,
                              inv=inv)


def plot_mean_moment_gc(sto_files, parameter_values, var_name, side, title, export_path, healthy_sto,
                        healthy_value, inv=False, plot=False):
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
    plot: bool to show plot

    Returns
    -------
    None
    '''
    # store metrics
    mdelta_mom_peaks = []
    sdelta_mom_peaks = []
    mmean_moment = []
    smean_moment = []

    fig, ax = plt.subplots()
    fig2, ax2 = plt.subplots()

    min_value = np.min(parameter_values)
    max_value = np.max(parameter_values)
    color_offset = mcolors.Normalize(vmin=min_value, vmax=max_value)

    # start with healthy plot
    parameter_value = healthy_value
    sto_file = healthy_sto

    # moment
    moment_norm = extract_sto.norm_moment(sto_file, var_name, side)

    joint_name = var_name.split('_')[0]
    var_names, var_tab = extract_sto.extract_sto(sto_file)
    h_time = var_tab[:, 0]
    _, _, _, lend, lst = extract_sto.me_mean_gait_phases(var_names, var_tab, 'time', side)

    h_av_stance_mom, h_av_swing_mom, h_av_gait_cycle, h_av_stance_end = extract_sto.mean_gait_phases(var_names, var_tab,
                                                                                                     var_name, side,
                                                                                                     moment_norm)
    joint_vel_name = joint_name + '_angle_' + side + '_u'
    h_av_stance_vel, h_av_swing_vel, h_av_gait_cycle, h_av_stance_end = extract_sto.mean_gait_phases(var_names, var_tab,
                                                                                                     joint_vel_name, side)
    h_av_mom = np.concatenate((h_av_stance_mom, h_av_swing_mom))

    h_mom_peaks, _ = find_peaks(h_av_stance_mom, distance=len(h_av_stance_mom))
    h_mdelta, h_sdelta, h_mmean, h_smean = extract_sto.moment_metrics(var_names, var_tab, side, moment_norm)

    label = str(parameter_value)
    ax.plot(h_time[:len(h_av_mom)] * 100 / h_time[h_av_gait_cycle], h_av_mom, label=label, color='black', linewidth=2,
            linestyle='dashed', zorder=10)

    ax.axvspan(h_time[h_av_stance_end] * 100 / h_time[h_av_gait_cycle], 100, alpha=0.5, color='blanchedalmond', zorder=0)

    ax2.axvspan(h_time[h_av_stance_end] * 100 / h_time[h_av_gait_cycle], 100, alpha=0.5, color='blanchedalmond', zorder=0)

    if min_value > healthy_value:
        hex_list = ['#0000ff', '#ff0000']
    else:
        hex_list = ['#ff0000', '#0000ff']
    cmap = get_continuous_cmap(hex_list)

    # plot all altered files
    for idx in range(0, len(sto_files)):
        sto_file = sto_files[idx]
        parameter_value = parameter_values[idx]

        # moment
        moment_norm = extract_sto.norm_moment(sto_file, var_name, side)

        var_names, var_tab = extract_sto.extract_sto(sto_file)
        time = var_tab[:, 0]
        av_stance_mom, av_swing_mom, av_gait_cycle, av_stance_end = extract_sto.mean_gait_phases(var_names, var_tab,
                                                                                                 var_name, side,
                                                                                                 moment_norm)
        joint_vel_name = joint_name + '_angle_' + side + '_u'
        try:
            av_stance_vel, av_swing_vel, av_gait_cycle, av_stance_end = extract_sto.mean_gait_phases(var_names, var_tab,
                                                                                                     joint_vel_name, side)
            av_mom = np.concatenate((av_stance_mom, av_swing_mom))
        except:
            av_mom = np.nan

        # compute moment peaks
        if not np.isnan(av_mom).any():
            first_mom_peaks, _ = find_peaks(av_stance_mom[:int(0.33*len(av_stance_mom))],
                                            distance=int(0.33*len(av_stance_mom)))
            sec_mom_peaks, _ = find_peaks(av_stance_mom[int(0.33 * len(av_stance_mom)):],
                                            distance=int(0.66*len(av_stance_mom)))

            label = str(parameter_value)
            ax.plot(time[:len(av_mom)] * 100 / time[av_gait_cycle], av_mom, label=label,
                    color=cmap(color_offset(parameter_value)), zorder=10)
            if len(first_mom_peaks) > 0:
                ax.plot(time[first_mom_peaks] * 100 / time[av_gait_cycle], av_stance_mom[first_mom_peaks], marker='x',
                        color=cmap(color_offset(parameter_value)), linestyle='None', zorder=10)
            ax.plot(time[int(0.33 * len(av_stance_mom))+sec_mom_peaks] * 100 / time[av_gait_cycle],
                    av_stance_mom[int(0.33 * len(av_stance_mom))+sec_mom_peaks], marker='x',
                    color=cmap(color_offset(parameter_value)), linestyle='None', zorder=10)

            # compute moment metrics
            mdelta, sdelta, mmean, smean = extract_sto.moment_metrics(var_names, var_tab, side, moment_norm)
        else:
            mdelta, sdelta, mmean, smean = np.nan, np.nan, np.nan, np.nan
        mdelta_mom_peaks.append(mdelta)
        sdelta_mom_peaks.append(sdelta)
        mmean_moment.append(mmean)
        smean_moment.append(smean)

    # plot healthy moment peak
    ax.plot(h_time[h_mom_peaks] * 100 / h_time[h_av_gait_cycle], h_av_stance_mom[h_mom_peaks], marker='x', label='peaks',
            color='black', linestyle='None', zorder=10)

    xlabel = title
    ylabel = joint_name + ' moment'
    ax.set_ylabel(ylabel + ' (' + side + ') [Nm/kg]')
    title = ylabel.split(" ")
    title[0] = title[0].capitalize()
    title = " ".join(title)
    ax.set_ylim((-0.5, 2))

    # sort both labels and handles by labels
    handles, labels = ax.get_legend_handles_labels()
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
        elif i == len(handles) - last:
            handles[i] = handles[-last]
            labels[i] = labels[-last]
    ax.legend(handles, labels, prop={'size': 12})
    ax2.legend(handles, labels, prop={'size': 12})
    plot_healthy_clinical_angles(ax, var_name)

    ax.set_xlabel('gait cycle [%]')
    ax.set_xlim((0, 100))
    ax2.set_xlabel('gait cycle [%]')
    ax2.set_xlim((0, 100))
    parameter_values = np.array(labls).astype(float).astype(int)
    if inv:
        ax.set_title(title + " for" + xlabel + " from " + str(max(parameter_values))
                     + " to " + str(min(parameter_values)) + "%")
    else:
        ax.set_title(title + " for" + xlabel + " from " + str(min(parameter_values))
                     + " to " + str(max(parameter_values)) + "%")

    fig.tight_layout()
    fig.savefig(os.path.join(export_path, var_name + '.png'))
    if not plot:
        plt.close(fig)

    return mdelta_mom_peaks, sdelta_mom_peaks, mmean_moment, smean_moment, h_mdelta, h_sdelta, h_mmean, h_smean


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

