# Utility functions.
#
# author: Dimitar Stanev <jimstanev@gmail.com>
import re
import opensim
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import matplotlib.colors as mcolors
import seaborn as sns
import os
import opensim

sns.set_theme()
sns.set_context("paper", font_scale=1.5)
sns.set_style("whitegrid")


def osim_array_to_list(array):
    """Converts OpenSim::Array<T> to Python list.

    Parameters
    ----------
    array: (array)

    Returns
    -------
    list: (list)
    """
    temp = []
    for i in range(array.getSize()):
        temp.append(array.get(i))

    return temp


def to_gait_cycle(data_frame, t0, tf):
    """Converts gait cycle time to gait cycle pourcentage.

    Parameters
    ----------
    data_frame: (panda df)
    t0: (float) gait cycle initial time
    tf: (float) gait cycle final time

    Returns
    -------
    data_frame: (panda df)
    """
    temp = data_frame[(data_frame.time >= t0) & (data_frame.time <= tf)].copy()
    temp.time = data_frame['time'].transform(lambda x: 100.0 / (tf - t0) * (x - t0))
    temp.set_index('time', inplace=True)
    temp.index.names = ['gait cycle [%]']
    return temp


def read_from_storage(file_name, to_filter=False):
    """Reads OpenSim.Storage files.

    Parameters
    ----------
    file_name: (string) path to file

    Returns
    -------
    tuple: (labels, time, data)
    """
    sto = opensim.Storage(file_name)
    sto.resampleLinear(0.01)
    if to_filter:
        sto.lowpassFIR(4, 6)

    labels = osim_array_to_list(sto.getColumnLabels())
    time = opensim.ArrayDouble()
    sto.getTimeColumn(time)
    time = np.round(osim_array_to_list(time), 3)
    data = []
    for i in range(sto.getSize()):
        temp = osim_array_to_list(sto.getStateVector(i).getData())
        temp.insert(0, time[i])
        data.append(temp)

    df = pd.DataFrame(data, columns=labels)
    df.index = df.time
    return df


def add_emg_timings(ax, var_name):
    """Adds experimental EMG timings.

    Parameters
    ----------
    ax: (Matplotlib <AxesSubplot>) axes subplot object to plot on
    var_name: (string) name of the variable to plot

    Returns
    -------
    None
    """
    color = 'grey'
    if "glut_max" in var_name:
        ax.axhline(y=1.2, xmin=0.00, xmax=0.25, color=color, linewidth=10)
        ax.axhline(y=1.2, xmin=0.95, xmax=1.00, color=color, linewidth=10)
    if "psoas" in var_name:
        ax.axhline(y=1.2, xmin=0.65, xmax=0.75, color=color, linewidth=10)
    if "hamstrings" in var_name:
        ax.axhline(y=1.2, xmin=0.00, xmax=0.20, color=color, linewidth=10)
        ax.axhline(y=1.2, xmin=0.80, xmax=1.00, color=color, linewidth=10)
    if "bifemsh" in var_name:
        ax.axhline(y=1.2, xmin=0.65, xmax=0.85, color=color, linewidth=10)
    if "vast" in var_name:
        ax.axhline(y=1.2, xmin=0.00, xmax=0.20, color=color, linewidth=10)
        ax.axhline(y=1.2, xmin=0.85, xmax=1.00, color=color, linewidth=10)
    if "rect_fem" in var_name:
        ax.axhline(y=1.2, xmin=0.55, xmax=0.65, color=color, linewidth=10)
    if "tib_ant" in var_name:
        ax.axhline(y=1.2, xmin=0.00, xmax=0.15, color=color, linewidth=10)
        ax.axhline(y=1.2, xmin=0.60, xmax=1.00, color=color, linewidth=10)
    if "gas" in var_name:
        ax.axhline(y=1.2, xmin=0.05, xmax=0.50, color=color, linewidth=10)
    if "sol" in var_name:
        ax.axhline(y=1.2, xmin=0.05, xmax=0.55, color=color, linewidth=10, label='experimental data')


def add_healthy_range_schwartz(ax, var_name):
    """Plots the value ranges of 'var_name' from healthy clinical observations.

    Source (free gait): https://doi.org/10.1016/j.jbiomech.2008.03.015

    Parameters
    ----------
    ax: (Matplotlib <AxesSubplot>) axes subplot object to plot on
    var_name: (string) name of the variable to plot

    Returns
    -------
    None
    """
    norm_min = None
    norm_max = None
    if "pelvis_tilt" in var_name:
        norm_min = np.multiply([7.2, 7.1, 7.0, 6.7, 6.4, 6.1, 6.0,
                                5.9, 5.9, 6.0, 6.2, 6.4, 6.6, 6.8,
                                6.9, 7.0, 7.1, 7.2, 7.3, 7.4, 7.5,
                                7.5, 7.5, 7.4, 7.2, 7.1, 6.9, 6.7,
                                6.5, 6.2, 5.9, 5.7, 5.7, 5.7, 5.8,
                                5.9, 6.1, 6.3, 6.5, 6.6, 6.7, 6.8,
                                7.0, 7.1, 7.2, 7.3, 7.4, 7.3, 7.2,
                                7.1, 6.9 ], -1)
        norm_mean = np.multiply([12.2, 12.0, 11.8, 11.6, 11.2, 10.9,
                                 10.8, 10.7, 10.7, 10.8, 11.0, 11.2,
                                 11.4, 11.5, 11.6, 11.7, 11.8, 12.0,
                                 12.1, 12.3, 12.4, 12.5, 12.5, 12.4,
                                 12.2, 12.0, 11.9, 11.7, 11.4, 11.1,
                                 10.8, 10.6, 10.5, 10.6, 10.7, 10.9,
                                 11.0, 11.2, 11.4, 11.5, 11.6, 11.7,
                                 11.8, 12.0, 12.2, 12.3, 12.4, 12.4,
                                 12.3, 12.1, 11.9 ], -1)
        norm_max = np.multiply([17.2, 17.0, 16.7, 16.4, 16.0, 15.7,
                                15.5, 15.5, 15.5, 15.7, 15.9, 16.0,
                                16.2, 16.3, 16.4, 16.5, 16.6, 16.7,
                                16.9, 17.2, 17.3, 17.4, 17.4, 17.3,
                                17.2, 17.0, 16.8, 16.6, 16.3, 16.0,
                                15.7, 15.5, 15.4, 15.4, 15.6, 15.8,
                                16.0, 16.1, 16.2, 16.3, 16.4, 16.5,
                                16.7, 16.9, 17.1, 17.3, 17.4, 17.4,
                                17.3, 17.2, 17.0 ], -1)
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
    if "knee_angle" in var_name:
        norm_min = [-0.1, 2.5, 5.3, 8.3, 10.8, 12.3, 13.0, 12.9, 12.4,
                    11.6, 10.6, 9.4, 8.2, 7.0, 5.8, 4.6, 3.5, 2.5,
                    1.6, 0.9, 0.5, 0.5, 1.1, 2.2, 3.8, 6.0, 8.7, 12.2,
                    16.4, 21.5, 27.3, 33.6, 39.9, 45.7, 50.0, 52.5,
                    53.3, 52.4, 50.1, 46.5, 42.0, 36.5, 30.3, 23.7,
                    16.8, 10.2, 4.4, 0.2, -2.0, -1.9, -0.2 ]
        norm_mean = [5.6, 7.9, 10.9, 14.1, 16.9, 18.6, 19.3, 19.2,
                     18.6, 17.6, 16.5, 15.2, 13.9, 12.6, 11.3, 10.1,
                     9.0, 8.0, 7.2, 6.6, 6.2, 6.3, 6.8, 7.9, 9.5,
                     11.6, 14.5, 18.1, 22.6, 27.9, 33.9, 40.2, 46.3,
                     51.7, 55.7, 58.2, 59.1, 58.7, 56.9, 53.9, 49.8,
                     44.7, 38.8, 32.4, 25.5, 18.7, 12.5, 7.6, 4.6,
                     4.0, 5.4 ]
        norm_max = [11.2, 13.4, 16.5, 20.0, 23.0, 24.8, 25.6, 25.5,
                    24.8, 23.7, 22.4, 21.0, 19.6, 18.2, 16.8, 15.6,
                    14.6, 13.6, 12.9, 12.3, 12.0, 12.0, 12.5, 13.5,
                    15.1, 17.3, 20.2, 24.0, 28.8, 34.3, 40.5, 46.8,
                    52.8, 57.7, 61.4, 63.8, 65.0, 65.0, 63.8, 61.3,
                    57.6, 52.9, 47.3, 41.0, 34.2, 27.2, 20.5, 14.9,
                    11.2, 10.0, 11.0 ]
    if "ankle_angle" in var_name:
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
        print("Ankle std: ", np.mean(std), "° , during stance: ", np.mean(std[:int(0.6*len(std))]), "°")
    if "grf_norm_y" in var_name:
        norm_min = [0.05, 0.61, 0.76, 0.86, 1.00, 1.15, 1.22, 1.22,
                    1.17, 1.10, 1.03, 0.96, 0.92, 0.89, 0.89, 0.89,
                    0.91, 0.93, 0.96, 1.01, 1.06, 1.11, 1.16, 1.19,
                    1.19, 1.15, 1.06, 0.90, 0.65, 0.36, 0.15, 0.05,
                    0.02, 0.01, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
                    0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
                    0.00, 0.00, 0.00]
        norm_mean = [0.01, 0.47, 0.60, 0.69, 0.84, 0.99, 1.08, 1.10,
                     1.06, 1.00, 0.93, 0.87, 0.82, 0.80, 0.79, 0.80,
                     0.81, 0.84, 0.87, 0.91, 0.96, 1.01, 1.06, 1.08,
                     1.07, 1.02, 0.91, 0.71, 0.45, 0.21, 0.07, 0.02,
                     0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
                     0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
                     0.00, 0.00, 0.00]
        norm_max= [-0.03, 0.33, 0.43, 0.52, 0.68, 0.83, 0.93, 0.97,
                   0.95, 0.89, 0.82, 0.77, 0.73, 0.71, 0.70, 0.70,
                   0.72, 0.74, 0.78, 0.82, 0.86, 0.91, 0.95, 0.97,
                   0.96, 0.89, 0.75, 0.53, 0.25, 0.05, -0.01, -0.01,
                   -0.01, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
                   0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
                   0.00, 0.00, 0.00]
    if 'hip_flexion_moment' in var_name:
        norm_min = [-0.01, 0.37, 0.23, 0.33, 0.30, 0.28, 0.26, 0.21,
                    0.15, 0.08, 0.01, -0.02, -0.04, -0.07, -0.09,
                    -0.13, -0.17, -0.23, -0.27, -0.34, -0.42, -0.51,
                    -0.62, -0.72, -0.80, -0.85, -0.84, -0.76, -0.61,
                    -0.50, -0.46, -0.40, -0.33, -0.28, -0.22, -0.17,
                    -0.13, -0.10, -0.08, -0.07, -0.07, -0.08, -0.07,
                    -0.05, 0.00, 0.08, 0.15, 0.20, 0.20, 0.13, -0.04]
        norm_mean = [0.12, 0.73, 0.54, 0.52, 0.50, 0.48, 0.44, 0.39,
                     0.32, 0.24, 0.17, 0.12, 0.09, 0.06, 0.03, 0.00,
                     -0.03, -0.08, -0.12, -0.18, -0.25, -0.33, -0.42,
                     -0.50, -0.58, -0.62, -0.62, -0.55, -0.43, -0.34,
                     -0.32, -0.27, -0.21, -0.17, -0.14, -0.10, -0.08,
                     -0.05, -0.03, -0.03, -0.02, -0.02, 0.00, 0.04,
                     0.11, 0.20, 0.28, 0.33, 0.34, 0.28, 0.12]
        norm_max = [0.26, 1.08, 0.84, 0.72, 0.70, 0.67, 0.63, 0.57,
                    0.48, 0.40, 0.32, 0.27, 0.23, 0.19, 0.16, 0.13,
                    0.10, 0.08, 0.04, -0.02, -0.08, -0.14, -0.21,
                    -0.29, -0.36, -0.40, -0.40, -0.35, -0.25, -0.19,
                    -0.18, -0.14, -0.10, -0.07, -0.05, -0.03, -0.02,
                    0.00, 0.01, 0.02, 0.03, 0.04, 0.08, 0.14, 0.22,
                    0.32, 0.41, 0.47, 0.49, 0.43, 0.28]
    if 'knee_flexion_moment' in var_name:
        norm_min = [-0.21, -0.47, -0.25, -0.07, 0.03, 0.10, 0.13,
                    0.13, 0.09, 0.04, -0.02, -0.07, -0.13, -0.17,
                    -0.20, -0.24, -0.27, -0.30, -0.34, -0.36, -0.38,
                    -0.39, -0.37, -0.32, -0.25, -0.17, -0.08, 0.00,
                    0.04, 0.04, 0.04, 0.02, 0.01, 0.02, 0.02, 0.01,
                    0.00, -0.01, -0.01, -0.02, -0.03, -0.04, -0.07,
                    -0.10, -0.15, -0.21, -0.28, -0.34, -0.36, -0.32,
                    -0.22]
        norm_mean = [-0.13, -0.29, -0.07, 0.08, 0.23, 0.33, 0.37,
                     0.35, 0.30, 0.24, 0.17, 0.10, 0.03, -0.02, -0.06,
                     -0.10, -0.14, -0.17, -0.20, -0.23, -0.24, -0.24,
                     -0.22, -0.17, -0.10, -0.02, 0.06, 0.12, 0.14,
                     0.12, 0.10, 0.08, 0.08, 0.09, 0.08, 0.07, 0.05,
                     0.03, 0.02, 0.01, 0.00, -0.02, -0.03, -0.06,
                     -0.10, -0.15, -0.20, -0.24, -0.26, -0.22, -0.13]
        norm_max = [-0.05, -0.12, 0.11, 0.22, 0.42, 0.55, 0.60, 0.58,
                    0.52, 0.43, 0.35, 0.27, 0.19, 0.13, 0.08, 0.03,
                    -0.01, -0.04, -0.07, -0.09, -0.10, -0.10, -0.07,
                    -0.02, 0.05, 0.12, 0.20, 0.24, 0.23, 0.19, 0.17,
                    0.15, 0.14, 0.15, 0.14, 0.12, 0.10, 0.07, 0.05,
                    0.03, 0.02, 0.01, 0.00, -0.02, -0.04, -0.08,
                    -0.12, -0.15, -0.16, -0.13, -0.04]
    if 'ankle_angle_moment' in var_name:
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
    norm_max = np.array(norm_mean) + 2 * std
    norm_min = np.array(norm_mean) - 2 * std

    if norm_min is not None and norm_max is not None:
        time_vector = np.linspace(0, 100, len(norm_min))
        ax.fill_between(time_vector, norm_min, norm_max,
                        facecolor='silver', alpha=0.5, label='clin.', zorder=5)


def get_heel_strike_events(state, side):
    """Gets the time instances of heel strike.

    Parameters
    ----------
    state: (tuple labels, time, data)
    side: (string) 'r' or 'l' for right or left side

    Returns
    -------
    list: (list)
    """
    assert(side in ['r', 'l'])
    if side == 'r':
        leg_name = 'leg1_r.state'
    else:
        leg_name = 'leg0_l.state'

    leg = state[leg_name].round()
    heel_strik = leg.diff() < 0
    return leg[heel_strik].index.tolist()


def get_toe_strike_events(state, side):
    """Gets the time instances of heel strike.

    Parameters
    ----------
    state: (tuple labels, time, data)
    side: (string) 'r' or 'l' for right or left side

    Returns
    -------
    list: (list)
    """
    assert(side in ['r', 'l'])
    if side == 'r':
        leg_name = 'leg1_r.state'
    else:
        leg_name = 'leg0_l.state'

    leg = state[leg_name].round()
    val = leg.values
    toe_strik = val == 3
    swing = np.zeros(len(leg))
    for i in range(len(swing)):
        if toe_strik[i]:
            swing[i] = 1
    toe_strik = np.diff(swing) > 0
    toe_strik = np.concatenate(([False], toe_strik))
    return leg[toe_strik].index.tolist()


def plot_scone_kinematics(state, side, output_file=None):
    """Plots 2D kinematic angles.

    Parameters
    ----------
    state: (tuple labels, time, data)
    side: (string) 'r' or 'l' for right or left side
    output_file: (string) path to save the plot

    Returns
    -------
    None
    """
    assert(side in ['r', 'l'])

    # calculate heel strike events
    hs_time = get_heel_strike_events(state, side)
    ts_time = np.array(get_toe_strike_events(state, side))
    hs_t = np.array(hs_time)
    if hs_time[0] < ts_time[0]:
        sw = np.mean((ts_time[1:] - hs_t[1:]) * 100 / (hs_t[1:] - hs_t[:-1]))
    else:
        if len(ts_time[1:]) < len(hs_t[:]):
            sw = np.mean((ts_time[1:] - hs_t[:-1]) * 100 / (hs_t[1:] - hs_t[:-1]))  # hs_t[:-1]
        elif len(ts_time[1:]) > len(hs_t[:]):
            sw = np.mean((ts_time[1:-1] - hs_t[:]) * 100 / (hs_t[1:] - hs_t[:-1]))
        else:
            sw = np.mean((ts_time[1:-1] - hs_t[:-1]) * 100 / (hs_t[1:] - hs_t[:-1]))

    # get kinematic variables
    kinematics = state[['time',
                        'pelvis_tilt',
                        'hip_flexion_' + side,
                        'knee_angle_' + side,
                        'ankle_angle_' + side]]
    if side == 'r':
        leg_name = 'leg1_r.grf_norm_'

    else:
        leg_name = 'leg0_l.grf_norm_'

    forces = state[['time', leg_name + 'x', leg_name + 'y']]

    hex_list = ['#00e8ff', '#020291']
    cmap = get_continuous_cmap(hex_list)
    color_offset = mcolors.Normalize(vmin=0, vmax=len(hs_time) - 1)

    # plot kinematics
    fig, ax = plt.subplots(2, 3, figsize=(9, 6))
    ax = ax.flatten()
    ax[0].axvspan(sw, 100, alpha=0.5, color='blanchedalmond', zorder=0)
    ax[1].axvspan(sw, 100, alpha=0.5, color='blanchedalmond', zorder=0)
    ax[2].axvspan(sw, 100, alpha=0.5, color='blanchedalmond', zorder=0)
    ax[3].axvspan(sw, 100, alpha=0.5, color='blanchedalmond', zorder=0)
    ax[4].axvspan(sw, 100, alpha=0.5, color='blanchedalmond', zorder=0)
    for i in range(len(hs_time) - 1):
        normalized = to_gait_cycle(kinematics, hs_time[i], hs_time[i+1])
        normalized = normalized.apply(lambda x: np.rad2deg(x))
        normalized['pelvis_tilt'].plot(ax=ax[0], color=cmap(color_offset(i)), zorder=10)
        normalized['hip_flexion_' + side].plot(ax=ax[1], color=cmap(color_offset(i)), zorder=10)
        normalized['knee_angle_' + side].apply(lambda x: -x).plot(ax=ax[2],color=cmap(color_offset(i)), zorder=10)
        normalized['ankle_angle_' + side].plot(ax=ax[3],color=cmap(color_offset(i)), zorder=10)
        normalized = to_gait_cycle(forces, hs_time[i], hs_time[i + 1])
        normalized['leg1_'+side+'.grf_norm_y'].plot(ax=ax[4], color=cmap(color_offset(i)), zorder=10)

    ax[0].set_ylabel('pelvis tilt [°]')
    ax[0].set_title('Pelvis tilt', fontweight="bold")
    ax[0].set_xlabel('gait cycle [%]')
    ax[0].set_xlim([0, 100])
    ax[0].set_xticks(range(0, 120, 20))
    ax[1].set_ylabel('hip flexion [°]')
    ax[1].set_title('Hip angle', fontweight="bold")
    ax[1].set_xlabel('gait cycle [%]')
    ax[1].set_xlim([0, 100])
    ax[1].set_xticks(range(0, 120, 20))
    ax[2].set_ylabel('knee angle [°]')
    ax[2].set_title('Knee angle', fontweight="bold")
    ax[2].set_xlabel('gait cycle [%]')
    ax[2].set_xlim([0, 100])
    ax[2].set_xticks(range(0, 120, 20))
    ax[3].set_ylabel('ankle angle [°]')
    ax[3].set_title('Ankle angle', fontweight="bold")
    ax[3].set_xlabel('gait cycle [%]')
    ax[3].set_xlim([0, 100])
    ax[3].set_xticks(range(0, 120, 20))
    ax[4].set_ylabel('GRF [BW]')
    ax[4].set_title('Ground reaction force', fontweight="bold")
    ax[4].set_xlabel('gait cycle [%]')
    ax[4].set_xlim([0, 100])
    ax[4].set_xticks(range(0, 120, 20))
    add_healthy_range_schwartz(ax[0], 'pelvis_tilt')
    add_healthy_range_schwartz(ax[1], 'hip_flexion')
    add_healthy_range_schwartz(ax[2], 'knee_angle')
    add_healthy_range_schwartz(ax[3], 'ankle_angle')
    add_healthy_range_schwartz(ax[4], 'grf_norm_y')
    fig.tight_layout()
    plt.subplots_adjust()
    # remove unused plots
    fig.delaxes(ax[5])
    if output_file is not None:
        fig.savefig(output_file + '.kinematics.png', bbox_inches='tight')


def plot_scone_joint_kinematics(state, model_file, muscle_analysis_output_dir, side, output_file=None):
    """Computes the joint moments obtained from OpenSim's MuscleAnalysis and plot joint kinemtics
    tool.
    Parameters
    ----------
    state: (tuple labels, time, data)
    model_file: (string) path to osim model file
    muscle_analysis_output_dir: (string) path to muscle analysis directory
    side: (string) 'r' or 'l' for right or left side
    output_file: (string) path to save the plot

    Returns
    -------
    None
    """
    assert (side in ['r', 'l'])
    # calculate heel strike events
    hs_time = get_heel_strike_events(state, side)
    ts_time = np.array(get_toe_strike_events(state, side))
    hs_t = np.array(hs_time)
    if hs_time[0] < ts_time[0]:
        sw = np.mean((ts_time[1:min(len(ts_time), len(hs_time))] - hs_t[1:min(len(ts_time), len(hs_time))]) * 100 /
                     (hs_t[1:min(len(ts_time), len(hs_time))] - hs_t[:min(len(ts_time), len(hs_time))-1]))
    else:
        sw = np.mean((ts_time[1:min(len(ts_time), len(hs_time))] - hs_t[:min(len(ts_time), len(hs_time))-1]) * 100 /
                     (hs_t[1:min(len(ts_time), len(hs_time))] - hs_t[:min(len(ts_time), len(hs_time))-1]))

    # get all files generated by the analysis tool
    _, _, filenames = next(os.walk(muscle_analysis_output_dir))

    # find moment files of interest
    ankle_moment_file = ''
    knee_moment_file = ''
    hip_moment_file = ''
    for f in filenames:
        if '_Moment_' + 'ankle_angle_' + side in f:
            ankle_moment_file = f
        if '_Moment_' + 'knee_angle_' + side in f:
            knee_moment_file = f
        if '_Moment_' + 'hip_flexion_' + side in f:
            hip_moment_file = f

    # load moment storage
    ankle_moment = read_from_storage(os.path.join(muscle_analysis_output_dir,
                                                  ankle_moment_file))
    knee_moment = read_from_storage(os.path.join(muscle_analysis_output_dir,
                                                 knee_moment_file))
    hip_moment = read_from_storage(os.path.join(muscle_analysis_output_dir,
                                                hip_moment_file))

    # get total mass for normalisation
    model = opensim.Model(model_file)
    s = model.initSystem()
    mass = model.getTotalMass(s)

    # compute normalized joint moments and place into a DataFrame
    ankle_moment_norm = ankle_moment.sum(axis=1) / mass
    knee_moment_norm = knee_moment.sum(axis=1) / mass
    hip_moment_norm = hip_moment.sum(axis=1) / mass
    moments = pd.concat([hip_moment_norm, knee_moment_norm, ankle_moment_norm],
                        axis=1)
    moments['time'] = moments.index
    moments.columns = ['hip_flexion_moment_' + side, 'knee_angle_moment_' + side,
                       'ankle_angle_moment_' + side, 'time']

    # get kinematic variables
    kinematics = state[['time',
                        'pelvis_tilt',
                        'hip_flexion_' + side,
                        'knee_angle_' + side,
                        'ankle_angle_' + side]]
    if side == 'r':
        leg_name = 'leg1_r'
    else:
        leg_name = 'leg0_l'
    forces = state[['time', leg_name + '.grf_norm_y']]

    # plot kinematics
    hex_list = ['#00e8ff', '#020291']
    cmap = get_continuous_cmap(hex_list)
    color_offset = mcolors.Normalize(vmin=0, vmax=len(hs_time) - 1)
    fig, ax = plt.subplots(2, 4, figsize=(12, 6))
    ax = ax.flatten()
    for i in range(len(hs_time) - 1):
        normalized = to_gait_cycle(kinematics, hs_time[i], hs_time[i + 1])
        normalized = normalized.apply(lambda x: np.rad2deg(x))
        normalized['pelvis_tilt'].plot(ax=ax[0], color=cmap(color_offset(i)), zorder=10)
        normalized['hip_flexion_' + side].plot(ax=ax[1], color=cmap(color_offset(i)), zorder=10)
        normalized['knee_angle_' + side].apply(lambda x: -x).plot(ax=ax[2], color=cmap(color_offset(i)), zorder=10)
        normalized['ankle_angle_' + side].plot(ax=ax[3], color=cmap(color_offset(i)), zorder=10)
        normalized = to_gait_cycle(moments, hs_time[i], hs_time[i + 1])
        normalized['hip_flexion_moment_' + side].apply(lambda x: -x).plot(ax=ax[5], color=cmap(color_offset(i)),
                                                                          zorder=10)
        normalized['knee_angle_moment_' + side].plot(ax=ax[6], color=cmap(color_offset(i)), zorder=10)
        normalized['ankle_angle_moment_' + side].apply(lambda x: -x).plot(ax=ax[7], color=cmap(color_offset(i)),
                                                                          zorder=10)
        normalized = to_gait_cycle(forces, hs_time[i], hs_time[i + 1])
        normalized[leg_name + '.grf_norm_y'].plot(ax=ax[4], color=cmap(color_offset(i)), zorder=10)

    # plot healthy ranges
    add_healthy_range_schwartz(ax[0], 'pelvis_tilt')
    add_healthy_range_schwartz(ax[1], 'hip_flexion')
    add_healthy_range_schwartz(ax[2], 'knee_angle')
    add_healthy_range_schwartz(ax[3], 'ankle_angle')
    add_healthy_range_schwartz(ax[4], 'grf_norm_y')
    add_healthy_range_schwartz(ax[5], 'hip_flexion_moment')
    add_healthy_range_schwartz(ax[6], 'knee_flexion_moment')
    add_healthy_range_schwartz(ax[7], 'ankle_angle_moment')

    # add labels
    for p in range(8):
        ax[p].axvspan(sw, 100, alpha=0.5, color='blanchedalmond', zorder=0)
        ax[p].set_xlabel('gait cycle [%]')
        ax[p].set_xlim([0, 100])
        ax[p].set_xticks(range(0, 120, 20))
    ax[0].set_ylabel('pelvis tilt (' + side + ') [°]')
    ax[0].set_title('Pelvis tilt', fontweight="bold")
    ax[1].set_ylabel('hip flexion (' + side + ') [°]')
    ax[1].set_title('Hip angle', fontweight="bold")
    ax[2].set_ylabel('knee flexion (' + side + ') [°]')
    ax[2].set_title('Knee angle', fontweight="bold")
    ax[3].set_ylabel('ankle flexion (' + side + ') [°]')
    ax[3].set_title('Ankle angle', fontweight="bold")
    ax[4].set_ylabel('GRF (' + side + ') [BW]')
    ax[4].set_title('Ground reaction force', fontweight="bold")
    ax[5].set_title('Hip moment', fontweight="bold")
    ax[5].set_ylabel('hip moment (' + side + ') [Nm/kg]')
    ax[6].set_title('Knee moment', fontweight="bold")
    ax[6].set_ylabel('knee moment (' + side + ') [Nm/kg]')
    ax[7].set_title('Ankle moment', fontweight="bold")
    ax[7].set_ylabel('ankle moment (' + side + ') [Nm/kg]')

    fig.tight_layout()
    plt.subplots_adjust()

    if output_file is not None:
        fig.savefig(output_file + '.kinematics.png', bbox_inches='tight')


def hex_to_rgb(value):
    '''Converts hex to rgb colours

    Parameters
    ----------
    value: (string) 6 characters representing a hex colour.

    Returns
    -------
    list: (list) length 3 of RGB values
    '''
    value = value.strip("#")  # removes hash symbol if present
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))


def rgb_to_dec(value):
    '''Converts rgb to decimal colours (i.e. divides each value by 256)

    Parameters
    ----------
    value: (list) length 3 of RGB values

    Returns
    -------
    list: (list) length 3 of decimal values
    '''

    return [v/256 for v in value]


def get_continuous_cmap(hex_list, float_list=None):
    '''Creates and returns a color map that can be used in heat map figures.
    If float_list is not provided, colour map graduates linearly between each color in hex_list.
    If float_list is provided, each color in hex_list is mapped to the respective location in float_list.

    Parameters
    ----------
    hex_list: (list) hex code strings
    float_list: (list) floats between 0 and 1, same length as hex_list. Must start with 0 and end with 1.

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


def plot_scone_muscle_activations(state, muscles, side, col=4, output_file=None):
    """Plots muscle activations.

     Parameters
    ----------
    state: (tuple labels, time, data)
    muscles: (string array) muscle names
    side: (string) 'r' or 'l' for right or left side
    col: (int) number of subplot columns
    output_file: (string) path to save the plot

    Returns
    -------
    None
    """
    # get muscle activations
    mask = ['time']
    for muscle in muscles:
        mask.append(muscle + '_' + side + '.activation')

    activations = state[mask]

    # calculate heel strike events
    hs_time = get_heel_strike_events(state, side)
    ts_time = np.array(get_toe_strike_events(state, side))
    hs_t = np.array(hs_time)
    if hs_time[0] < ts_time[0]:
        sw = np.mean((ts_time[1:min(len(ts_time), len(hs_time))] - hs_t[1:min(len(ts_time), len(hs_time))]) * 100 /
                     (hs_t[1:min(len(ts_time), len(hs_time))] - hs_t[:min(len(ts_time), len(hs_time)) - 1]))
    else:
        sw = np.mean((ts_time[1:min(len(ts_time), len(hs_time))] - hs_t[:min(len(ts_time), len(hs_time)) - 1]) * 100 /
                     (hs_t[1:min(len(ts_time), len(hs_time))] - hs_t[:min(len(ts_time), len(hs_time)) - 1]))

    # visualize muscle activations
    M = len(muscles)
    if M % col == 0:
        N = int(np.floor(M / col))
    else:
        N = int(np.floor(M / col)) + 1

    fig, ax = plt.subplots(N, col, figsize=(3 * col, 3 * N))
    ax = ax.flatten()
    hex_list = ['#00e8ff', '#020291']
    cmap = get_continuous_cmap(hex_list)
    color_offset = mcolors.Normalize(vmin=0, vmax = len(hs_time)-1)

    # add gait cycle plots
    for i in range(len(hs_time) - 1):
        normalized = to_gait_cycle(activations, hs_time[i],
                                   hs_time[i+1])
        for j in range(len(muscles)):
            normalized[mask[j + 1]].plot(ax=ax[j], color=cmap(color_offset(i)))

    # configure subfigures and add range values
    muscle_names = ['ILPSO', 'GMAX', 'HAMS', 'VAS', 'TA', 'GAS', 'SOL']
    for i in range(len(muscles)):
        add_emg_timings(ax[i], mask[i + 1])
        ax[i].axvspan(sw, 100, alpha=0.5, color='blanchedalmond')
        ax[i].set_xlabel('gait cycle [%]')
        ax[i].set_title(muscle_names[i], fontweight="bold")
        ax[i].set_ylabel('muscle activation')
        ax[i].set_ylim([0, 1.2])
        ax[i].set_yticks(np.arange(0, 1.2, 0.2))
        ax[i].set_xlim([0, 100])
        ax[i].set_xticks(range(0, 120, 20))

    fig.tight_layout()

    # remove unused plots
    i = len(ax) - 1
    if M % col != 0:
        while i >= M:
            fig.delaxes(ax[i])
            i = i - 1
    fig.tight_layout()
    if output_file is not None:
        fig.savefig(output_file + '.activations.png', bbox_inches='tight')


def perform_muscle_analysis(model_file, state_file, output_dir):
    """Perform OpenSim MuscleAnalysis on SCONE state file generated
    through simulation.
    Parameters
    --------
    model_file: (string) path to osim model file
    state_file: (string) path to sto file
    output_dir: (string) path to directory to save analysis files

    Returns
    --------
    None
    """
    model = opensim.Model(model_file)

    # construct static optimization
    state_storage = opensim.Storage(state_file)
    muscle_analysis = opensim.MuscleAnalysis()
    muscle_analysis.setStartTime(state_storage.getFirstTime())
    muscle_analysis.setEndTime(state_storage.getLastTime())
    model.addAnalysis(muscle_analysis)

    # analysis
    analysis = opensim.AnalyzeTool(model)
    analysis.setName('muscle_analysis')
    analysis.setModel(model)
    analysis.setInitialTime(state_storage.getFirstTime())
    analysis.setFinalTime(state_storage.getLastTime())
    analysis.setStatesFileName(state_file)
    # analysis.setLowpassCutoffFrequency(6)
    analysis.setLoadModelAndInput(True)
    analysis.setResultsDir(os.path.abspath(output_dir))
    analysis.run()
