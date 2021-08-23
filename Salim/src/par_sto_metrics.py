"""
This file contains different methods used to extract and plot metrics of SCONE
simulation results.
"""

from pathlib import Path
import sys, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
sys.path.insert(0, 'alice_thesis_scripts')
import spasticity_index, extract_sto
import seaborn as sns

sns.set_theme()
sns.set_context("notebook", font_scale=1.2)

results_folder = '../results'

def hex_to_rgb(value):
    '''
    Converts hex to rgb colours
    value: string of 6 characters representing a hex colour.
    Returns: list length 3 of RGB values
    '''
    value = value.strip("#") # removes hash symbol if present
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))


def rgb_to_dec(value):
    '''
    Converts rgb to decimal colours (i.e. divides each value by 256)
    value: list (length 3) of RGB values
    Returns: list (length 3) of decimal values
    '''
    return [v/256 for v in value]


def get_continuous_cmap(hex_list, float_list=None):
    ''' 
    creates and returns a color map that can be used in heat map figures.
    If float_list is not provided, colour map graduates linearly between each color in hex_list.
    If float_list is provided, each color in hex_list is mapped to the respective location in float_list. 
    
    Parameters
    ----------
    hex_list: list of hex code strings
    float_list: list of floats between 0 and 1, same length as hex_list. Must start with 0 and end with 1.
    
    Returns
    ----------
    colour map
    '''
    rgb_list = [rgb_to_dec(hex_to_rgb(i)) for i in hex_list]
    if float_list:
        pass
    else:
        float_list = list(np.linspace(0,1,len(rgb_list)))
        
    cdict = dict()
    for num, col in enumerate(['red', 'green', 'blue']):
        col_list = [[float_list[i], rgb_list[i][num], rgb_list[i][num]] for i in range(len(float_list))]
        cdict[col] = col_list
    cmp = mcolors.LinearSegmentedColormap('my_cmp', segmentdata=cdict, N=256)
    return cmp

def plot_healthy_clinical_angles(ax, var_name):
    '''
    This method plots the value ranges of 'var_name' from healthy clinical
    observations.

    Parameters
    ----------
    ax : Matplotlib <AxesSubplot:>
        Axes subplot object to plot on.
    var_name : string
        Name of the variable to plot.

    Returns
    -------
    None.

    '''
    norm_min = None
    norm_mean = None
    norm_max = None
    if "pelvis_tilt" in var_name:
        norm_min = [-16.4,-16.4,-16.5,-16.7,-16.7,-16.7,-16.6,-16.5,-16.4,-16.3,-16.2,-16.2,-16.2,-16.3,-16.3,-16.4,-16.5,-16.5,-16.6,-16.6,-16.7,-16.7,-16.7,-16.7,-16.6,-16.6,-16.5,-16.5,-16.4,-16.4,-16.3,-16.3,-16.3,-16.2,-16.2,-16.2,-16.2,-16.1,-16.1,-16.1,-16.1,-16.1,-16.1,-16.1,-16.1,-16.2,-16.2,-16.2,-16.3,-16.3,-16.3,-16.3,-16.4,-16.4,-16.5,-16.5,-16.4,-16.3,-16.2,-16.1,-16.0,-16.0,-16.0,-16.0,-16.0,-16.1,-16.1,-16.2,-16.2,-16.2,-16.2,-16.2,-16.2,-16.1,-16.1,-16.1,-16.0,-16.0,-15.9,-15.9,-15.9,-15.8,-15.8,-15.8,-15.8,-15.7,-15.7,-15.7,-15.7,-15.7,-15.7,-15.6,-15.6,-15.6,-15.6,-15.7,-15.7,-15.8,-15.8,-15.9,-15.9]
        norm_mean = [-9.3,-9.4,-9.5,-9.6,-9.7,-9.7,-9.6,-9.4,-9.3,-9.2,-9.1,-9.0,-9.0,-9.0,-9.0,-9.1,-9.2,-9.2,-9.3,-9.3,-9.4,-9.4,-9.4,-9.4,-9.4,-9.3,-9.3,-9.2,-9.2,-9.1,-9.1,-9.1,-9.1,-9.0,-9.0,-9.0,-9.0,-9.0,-8.9,-8.9,-8.9,-8.9,-8.8,-8.9,-8.9,-8.9,-9.0,-9.0,-9.1,-9.1,-9.2,-9.2,-9.3,-9.4,-9.4,-9.4,-9.3,-9.2,-9.1,-8.9,-8.8,-8.8,-8.7,-8.7,-8.8,-8.8,-8.9,-8.9,-9.0,-9.0,-9.1,-9.1,-9.0,-9.0,-9.0,-9.0,-8.9,-8.9,-8.8,-8.8,-8.8,-8.7,-8.7,-8.7,-8.7,-8.7,-8.6,-8.6,-8.6,-8.5,-8.5,-8.5,-8.4,-8.4,-8.4,-8.5,-8.5,-8.6,-8.7,-8.7,-8.8]
        norm_max = [-2.3,-2.4,-2.5,-2.6,-2.6,-2.6,-2.5,-2.3,-2.2,-2.0,-1.9,-1.8,-1.8,-1.8,-1.8,-1.8,-1.8,-1.9,-2.0,-2.0,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.0,-2.0,-1.9,-1.9,-1.9,-1.9,-1.9,-1.9,-1.9,-1.8,-1.8,-1.8,-1.7,-1.7,-1.6,-1.6,-1.6,-1.6,-1.6,-1.7,-1.8,-1.8,-1.9,-2.0,-2.0,-2.1,-2.2,-2.3,-2.3,-2.3,-2.2,-2.1,-1.9,-1.8,-1.6,-1.5,-1.5,-1.5,-1.5,-1.5,-1.6,-1.7,-1.8,-1.9,-1.9,-1.9,-1.9,-1.9,-1.9,-1.9,-1.8,-1.8,-1.7,-1.7,-1.7,-1.6,-1.6,-1.6,-1.6,-1.6,-1.6,-1.5,-1.5,-1.4,-1.3,-1.3,-1.2,-1.2,-1.3,-1.3,-1.4,-1.5,-1.5,-1.6,-1.6]
    if "hip_flexion" in var_name:
        norm_min = [25.2,25.1,25.0,24.7,24.3,23.9,23.5,23.0,22.5,21.8,21.0,20.2,19.1,18.1,16.9,15.8,14.7,13.7,12.7,11.6,10.5,9.4,8.2,7.0,5.8,4.6,3.4,2.2,1.0,-0.2,-1.3,-2.4,-3.5,-4.6,-5.6,-6.6,-7.6,-8.6,-9.6,-10.5,-11.4,-12.2,-13.0,-13.6,-14.2,-14.7,-15.1,-15.4,-15.7,-16.0,-16.1,-16.2,-16.1,-15.9,-15.6,-15.1,-14.5,-13.7,-12.8,-11.7,-10.4,-9.0,-7.3,-5.5,-3.6,-1.7,0.3,2.4,4.4,6.4,8.3,10.1,11.8,13.5,15.1,16.6,18.0,19.3,20.6,21.7,22.6,23.5,24.2,24.7,25.1,25.3,25.4,25.3,25.1,24.8,24.4,24.1,23.8,23.6,23.5,23.5,23.6,23.7,23.8,24.0,24.1]
        norm_mean = [33.0,33.0,32.9,32.7,32.5,32.2,32.0,31.6,31.2,30.6,29.9,29.0,28.1,27.1,26.1,25.1,24.0,23.0,21.9,20.9,19.7,18.6,17.4,16.2,14.9,13.7,12.4,11.2,9.9,8.7,7.5,6.3,5.2,4.1,3.0,1.9,0.9,-0.1,-1.1,-2.0,-2.9,-3.7,-4.5,-5.2,-5.8,-6.4,-6.8,-7.3,-7.6,-7.9,-8.2,-8.3,-8.3,-8.1,-7.8,-7.3,-6.7,-5.9,-4.9,-3.7,-2.3,-0.7,1.1,3.0,5.0,7.0,9.0,10.9,12.9,14.9,16.7,18.5,20.2,21.9,23.4,24.9,26.3,27.7,28.9,30.0,30.9,31.8,32.5,33.0,33.4,33.7,33.7,33.7,33.5,33.1,32.8,32.4,32.0,31.7,31.5,31.4,31.4,31.5,31.7,31.9,32.0]
        norm_max = [40.8,40.9,40.8,40.8,40.7,40.6,40.5,40.3,39.9,39.3,38.7,37.9,37.1,36.2,35.3,34.3,33.3,32.3,31.2,30.1,29.0,27.8,26.6,25.4,24.1,22.8,21.5,20.2,18.9,17.6,16.4,15.1,13.9,12.7,11.6,10.5,9.4,8.4,7.5,6.5,5.6,4.8,4.0,3.2,2.5,1.9,1.4,0.9,0.5,0.1,-0.2,-0.4,-0.5,-0.4,0.0,0.4,1.1,1.9,3.0,4.3,5.8,7.6,9.5,11.5,13.5,15.6,17.6,19.5,21.5,23.3,25.2,26.9,28.6,30.2,31.8,33.3,34.7,36.0,37.2,38.3,39.2,40.1,40.8,41.3,41.7,42.0,42.0,42.0,41.8,41.5,41.2,40.7,40.3,39.9,39.6,39.4,39.3,39.4,39.6,39.8,40.0]
    if "knee_angle" in var_name:
        norm_min = np.multiply([4.4,5.9,6.9,7.5,8.2,9.2,10.6,12.0,13.1,13.8,14.3,14.6,14.7,14.5,14.1,13.5,12.9,12.3,11.7,11.2,10.6,10.1,9.6,9.0,8.5,8.0,7.5,6.9,6.4,5.9,5.3,4.8,4.4,4.0,3.6,3.4,3.1,3.0,2.9,2.9,2.9,3.1,3.4,3.8,4.4,5.1,5.9,6.8,7.8,8.8,9.8,11.0,12.3,13.7,15.3,17.1,19.2,21.6,24.2,27.1,30.2,33.5,36.9,40.3,43.6,46.7,49.5,52.0,54.2,55.9,57.2,58.1,58.5,58.6,58.2,57.5,56.4,54.9,53.0,50.8,48.3,45.5,42.3,38.8,35.0,31.1,27.0,22.7,18.5,14.5,10.6,7.2,4.3,2.0,0.2,-0.8,-1.3,-1.0,-0.1,1.3,2.9],-1)
        norm_mean = np.multiply([7.9,9.3,10.3,11.0,11.9,13.1,14.7,16.2,17.4,18.3,19.0,19.3,19.4,19.3,19.0,18.5,18.1,17.6,17.2,16.7,16.2,15.7,15.1,14.6,14.0,13.4,12.8,12.1,11.5,10.9,10.2,9.6,9.1,8.6,8.1,7.8,7.5,7.3,7.2,7.2,7.3,7.4,7.7,8.1,8.6,9.3,10.0,10.9,11.9,12.9,14.1,15.5,17.0,18.7,20.6,22.8,25.2,27.9,30.8,33.9,37.2,40.5,44.0,47.3,50.4,53.2,55.7,57.9,59.8,61.3,62.4,63.1,63.4,63.3,63.0,62.3,61.3,60.0,58.4,56.4,54.2,51.7,48.9,45.8,42.4,38.8,35.0,30.9,26.8,22.7,18.7,15.0,11.6,8.7,6.5,4.8,4.0,3.8,4.4,5.5,7.0],-1)
        norm_max = np.multiply([11.5,12.7,13.6,14.4,15.6,17.1,18.8,20.5,21.8,22.8,23.6,24.0,24.2,24.1,23.8,23.6,23.3,23.0,22.6,22.2,21.7,21.2,20.7,20.1,19.5,18.8,18.1,17.3,16.6,15.8,15.1,14.4,13.7,13.2,12.7,12.3,11.9,11.7,11.6,11.5,11.6,11.7,12.0,12.3,12.8,13.4,14.1,15.0,16.0,17.1,18.4,19.9,21.7,23.7,26.0,28.5,31.3,34.2,37.3,40.6,44.1,47.6,51.0,54.2,57.2,59.8,62.0,63.8,65.4,66.6,67.5,68.0,68.2,68.1,67.8,67.2,66.3,65.1,63.7,62.0,60.2,58.0,55.6,52.8,49.8,46.5,42.9,39.1,35.1,31.0,26.8,22.8,19.0,15.5,12.7,10.5,9.2,8.6,8.9,9.8,11.1],-1)
    if "ankle_angle" in var_name:
        norm_min = np.add([-25.1,-25.5,-26.7,-28.3,-29.3,-29.4,-28.7,-27.7,-26.7,-25.8,-25.0,-24.1,-23.2,-22.5,-21.8,-21.1,-20.4,-19.8,-19.3,-18.8,-18.3,-17.8,-17.3,-16.8,-16.4,-15.9,-15.5,-15.2,-14.8,-14.5,-14.2,-13.9,-13.6,-13.2,-12.8,-12.3,-11.9,-11.6,-11.2,-11.0,-10.7,-10.6,-10.5,-10.6,-10.7,-10.9,-11.3,-11.9,-12.6,-13.7,-14.9,-16.5,-18.5,-21.0,-23.9,-27.3,-31.1,-35.0,-38.7,-41.7,-44.2,-45.9,-47.1,-47.7,-47.7,-46.9,-45.2,-43.0,-40.6,-38.2,-36.1,-34.2,-32.4,-30.6,-28.9,-27.3,-26.0,-24.9,-24.0,-23.3,-22.5,-21.9,-21.3,-20.8,-20.5,-20.4,-20.6,-21.0,-21.7,-22.3,-22.9,-23.5,-24.0,-24.4,-24.8,-25.1,-25.3,-25.5,-25.4,-25.2,-25.1],20)
        norm_mean = np.add([-22.3,-22.8,-23.9,-25.2,-26.2,-26.5,-26.1,-25.5,-24.6,-23.6,-22.5,-21.5,-20.5,-19.7,-19.0,-18.3,-17.7,-17.1,-16.6,-16.0,-15.6,-15.1,-14.6,-14.1,-13.6,-13.2,-12.8,-12.4,-12.1,-11.7,-11.4,-11.1,-10.8,-10.4,-10.0,-9.6,-9.3,-8.9,-8.5,-8.2,-8.0,-7.7,-7.6,-7.5,-7.6,-7.7,-7.9,-8.3,-8.9,-9.7,-10.6,-11.9,-13.4,-15.3,-17.6,-20.3,-23.4,-26.6,-29.8,-32.8,-35.3,-37.4,-38.9,-39.8,-39.9,-39.3,-38.0,-36.4,-34.6,-32.6,-30.7,-28.9,-27.3,-25.8,-24.3,-23.0,-21.8,-20.8,-20.0,-19.3,-18.6,-18.1,-17.6,-17.4,-17.3,-17.3,-17.6,-18.1,-18.7,-19.4,-20.0,-20.6,-21.1,-21.5,-21.9,-22.2,-22.5,-22.7,-22.7,-22.5,-22.3],20)
        norm_max = np.add([-19.6,-20.1,-21.1,-22.2,-23.0,-23.5,-23.6,-23.2,-22.4,-21.3,-20.0,-18.8,-17.8,-16.9,-16.2,-15.5,-14.9,-14.3,-13.8,-13.3,-12.8,-12.3,-11.9,-11.4,-10.9,-10.5,-10.0,-9.7,-9.3,-8.9,-8.6,-8.3,-8.0,-7.6,-7.3,-6.9,-6.6,-6.2,-5.8,-5.5,-5.2,-4.9,-4.7,-4.5,-4.5,-4.5,-4.5,-4.8,-5.1,-5.7,-6.3,-7.2,-8.3,-9.7,-11.3,-13.3,-15.6,-18.2,-21.0,-23.8,-26.5,-28.9,-30.7,-31.8,-32.1,-31.6,-30.9,-29.8,-28.5,-27.0,-25.3,-23.6,-22.2,-20.9,-19.7,-18.6,-17.6,-16.7,-15.9,-15.3,-14.7,-14.3,-14.0,-13.9,-14.0,-14.3,-14.7,-15.2,-15.8,-16.5,-17.1,-17.7,-18.1,-18.6,-19.0,-19.3,-19.6,-19.9,-19.9,-19.8,-19.5],20)
    if "grf_norm_y" in var_name:
        norm_min = [0.00,0.31,0.41,0.45,0.47,0.51,0.56,0.62,0.69,0.75,0.81,0.86,0.91,0.94,0.96,0.97,0.97,0.96,0.94,0.92,0.89,0.85,0.82,0.80,0.77,0.75,0.74,0.73,0.72,0.72,0.72,0.72,0.73,0.74,0.76,0.78,0.81,0.83,0.86,0.89,0.91,0.94,0.96,0.99,1.01,1.03,1.05,1.06,1.07,1.06,1.05,1.01,0.95,0.87,0.76,0.62,0.47,0.32,0.19,0.10,0.02,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00]
        norm_mean = [0.00,0.40,0.52,0.56,0.58,0.61,0.67,0.73,0.80,0.88,0.94,0.99,1.03,1.05,1.05,1.05,1.04,1.03,1.01,0.98,0.95,0.93,0.90,0.87,0.85,0.83,0.82,0.81,0.80,0.80,0.80,0.80,0.81,0.82,0.83,0.85,0.87,0.89,0.91,0.94,0.97,1.00,1.03,1.06,1.08,1.11,1.13,1.14,1.14,1.14,1.12,1.09,1.05,0.99,0.91,0.80,0.68,0.54,0.41,0.30,0.19,0.12,0.07,0.03,0.01,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00]
        norm_max = [0.00,0.49,0.63,0.67,0.70,0.71,0.77,0.84,0.92,1.00,1.07,1.12,1.15,1.16,1.15,1.14,1.12,1.09,1.07,1.05,1.02,1.00,0.97,0.95,0.93,0.92,0.90,0.89,0.89,0.88,0.88,0.89,0.89,0.90,0.91,0.92,0.93,0.94,0.96,0.99,1.02,1.05,1.09,1.12,1.16,1.18,1.20,1.21,1.21,1.21,1.20,1.17,1.15,1.11,1.05,0.98,0.88,0.77,0.63,0.50,0.37,0.25,0.16,0.09,0.05,0.02,0.01,0.01,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00]
    if norm_min is not None and norm_max is not None:
        time_vector = np.linspace(0,100,len(norm_min))
        ax.fill_between(time_vector, norm_min, norm_max, facecolor='silver', alpha=0.5, label='clin.')
    
    

def plot_column(column_values, parameter_values, column_name, x_label, export_path, healthy_value, healthy_metric, plot = False):
    '''
    Plots values of a column with parameter_values on x axis.
    INPUTS:  - column_values: values to be plotted
             - parameter_values: values on x axis
             - column_name: name of the column, y axis legend
             - x_label: name of the parameters
             - export_path: folder path of export_file
             - Healthy value: value from healthy sto file of interest
             - Healthy metric: value of the metric for healthy optimization
    '''
    
    fig, ax = plt.subplots()
    parameter_values_extended = parameter_values + [healthy_value]
    column_values_extended = column_values + [healthy_metric]
    sorted_idx = np.argsort(parameter_values_extended)
    column_values_extended = [column_values_extended[idx] for idx in sorted_idx]
    parameter_values_extended = [parameter_values_extended[idx] for idx in sorted_idx]
    
    
    ax.plot(parameter_values_extended, column_values_extended, color='blue')
    ax.plot(healthy_value, healthy_metric, color='lime', marker='P', markersize=20)
    ax.set_xlabel(x_label + ', healthy = ' + str(healthy_value))
    ax.set_ylabel(column_name.replace('_',' '))
    plt.savefig(os.path.join(export_path,column_name+'.pdf'))
    if not plot:
        plt.close(fig)


def plot_mean_gc(sto_files, parameter_values, var_name, side, ylabel, title, export_path, healthy_sto, healthy_value, plot = False):
    '''
    Plots averaged variables during gait cycle, for multiple sto files
    INPUTS: - sto_file: path to the sto file of interest
            - var_list: list of variables to plot
            - var_legend: list of legends for the previous variables
            - side: 'r' or 'l' for right or left variables
            - ylabel: plot y label
            - title: plot title
            - export_path: folder path of export file
            - healthy_sto: path to healthy sto file of interest
            - healthy_value: value to healthy parameter
    '''
    var_names_list = []
    var_tab_list = []
    time_list = []
    
    fig, ax = plt.subplots()
    
    
    min_value = np.min(parameter_values)
    max_value = np.max(parameter_values)
    max_diff = np.max([max_value - healthy_value, healthy_value - min_value])
        
    try:
        color_offset = mcolors.TwoSlopeNorm(vmin=healthy_value - np.mean([max_diff, healthy_value - min_value]), vcenter=healthy_value, vmax=healthy_value + np.mean([max_diff, max_value - healthy_value]))
    except :
        color_offset = mcolors.TwoSlopeNorm(vcenter=healthy_value)
        
    color_offset = mcolors.Normalize(vmin=min_value,vmax = max_value)    
    # start with healthy plot
    parameter_value = healthy_value    
    sto_file = healthy_sto
    var_names, var_tab = extract_sto.extract_sto(sto_file)
    time = var_tab[:, 0]
    var_names_list.append(var_names)
    var_tab_list.append(var_tab)
    time_list.append(time)
    av_stance_var, av_swing_var, av_gait_cycle, av_stance_end = extract_sto.mean_gait_phases(var_names, var_tab, var_name, side)
    av_var = np.concatenate((av_stance_var, av_swing_var))
    if "angle" in var_name:
        av_var = np.multiply(av_var,180/np.pi)
    ax.plot(time[:len(av_var)]*100/time[av_gait_cycle], av_var, label=parameter_value, color= 'lime', linewidth=4, linestyle='dashed')
    
    hex_list = ['#0000ff', '#ff0000'] 
    cmap = get_continuous_cmap(hex_list)
    nb_plots = 5
    
    for idx in range(0,len(sto_files)):
        sto_file = sto_files[idx]
        parameter_value = parameter_values[idx]
        var_names, var_tab = extract_sto.extract_sto(sto_file)
        time = var_tab[:, 0]
        var_names_list.append(var_names)
        var_tab_list.append(var_tab)
        time_list.append(time)
        av_stance_var, av_swing_var, av_gait_cycle, av_stance_end = extract_sto.mean_gait_phases(var_names, var_tab, var_name, side)
        av_var = np.concatenate((av_stance_var, av_swing_var))
        if "angle" in var_name:
            av_var = np.multiply(av_var,180/np.pi)
        ax.plot(time[:len(av_var)]*100/time[av_gait_cycle], av_var, label=parameter_value,color= cmap(color_offset(parameter_value)))
        
    handles, labels = ax.get_legend_handles_labels()
    # sort both labels and handles by labels
    labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: float(t[0])))
    ax.legend(handles, labels)
    plot_healthy_clinical_angles(ax, var_name)

    ax.set_xlabel('Gait cycle (%)')
    ax.set_ylabel(ylabel)
    ax.set_title(title + ', healthy = ' + str(healthy_value))
    plt.savefig(os.path.join(export_path,var_name+'.pdf'))
    if not plot:
        plt.close(fig)
        
def get_sto_total_time(sto_file):
    '''
    Returns total simulation time from an sto file.

    '''
    var_names, var_tab = extract_sto.extract_sto(sto_file)
    time = var_tab[:, 0]
    return time[-1] # return max time

