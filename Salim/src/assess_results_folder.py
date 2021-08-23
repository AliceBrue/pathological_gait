"""
This script takes as input a path of a results folder with an organized folder structure. 
It evaluates different metrics for each of the parameters explored and 
generates a report folder containing visualization of these metrics.
To use this script, run the following terminal command:\
    python3  assess_results_folder.py *path/to/organized/results/folder*
"""

import sys, os
from pathlib import Path
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import par_sto_metrics
sys.path.insert(0, 'alice_thesis_scripts')
import spasticity_index, extract_sto

def get_experiment_title(parameter_folder):
    '''
    This method returns an experiment title based off the parameter folder and
    its parent folder.

    Parameters
    ----------
    parameter_folder : pathlib.Path
        Object defining parameter folder path.

    Returns
    -------
    string
        Title of the experiment.

    '''
    parameter_phase = parameter_folder.parent.name
    parameter_name = parameter_folder.name
    parameter_phase_name = parameter_phase + '_' + parameter_name
    
    return parameter_phase_name.replace('_', ' ')
    


def get_score_list(source_folder):
    '''
    Returns the child folder list of a source folder and sto scores corresponding 
    to each child folder.

    Parameters
    ----------
    source_folder : string
        path of source folder.

    Returns
    -------
    child_folders : list
        list containing child folder names.
    sto_scores : list
        scores corresponding to each child folder.

    '''
    child_folders = [child_folder for child_folder in os.listdir(source_folder) if os.path.isdir(os.path.join(source_folder, child_folder))]
    sto_files = [get_experiment_sto_file(os.path.join(source_folder, child_folder)) for child_folder in child_folders]
    
    sto_scores = [sto_file.rstrip('.par.sto').split('_')[-1] for sto_file in sto_files]
    
    return child_folders, sto_scores
    

def str_is_float(x):
    '''
    Returns whether a string can be converted to a float or not.

    Parameters
    ----------
    x : string
        potential float.

    Returns
    -------
    bool
        True if string can be converted to a float. False otherwise

    '''
    try:
        float(x)
        return True
    except ValueError:
        return False
    
def get_healthy_sto(parameter_folder):
    '''
    Returns path to healthy sto file based off parameter folder name.

    Parameters
    ----------
    parameter_folder : pathlib.Path
        Object defining parameter folder path.

    Returns
    -------
    str
        path to healthy sto file.

    '''
    parameter_folder_str = str(parameter_folder)
    folder_name = parameter_folder.name
    if 'simple' or 'solLV' in  parameter_folder_str:
        
        if ('KF' in folder_name) or ('KL' in folder_name):
            return '../models/healthy_initializations/init_geyer_swKF_005.38.par.sto'
        else:
            return '../models/healthy_initializations/init_geyer_good_signs_004.86.par.sto'
    else:
        return '../models/healthy_initializations/init_geyer_couplings_002.18.par.sto'
    
def get_healthy_gains_2d(parameter_folder):
    '''
    Returns the gains obtained from healthy optimization based off parameter
    folder path of 2d experiments.

    Parameters
    ----------
    parameter_folder : pathlib.Path
        Object defining parameter folder path.

    Returns
    -------
    output_1 : float
        Healthy stance gain.
    output_2 : float
        Healthy swing gain.

    '''
    parameter_folder_str = str(parameter_folder)
    if 'simple' or 'solLV' in  parameter_folder_str:
        healthy_gains = {'stance_KV': 0.12,
                          'swing_KV': 0.13,
                          'stance_KL': 1.75,
                          'swing_KL' : 0.62,
                          'stance_KF': 0.19, 
                          'swing_KF': -0.24}
    else:
        healthy_gains = {'stance_KV': -2.11,
                          'swing_KV': -0.21,
                          'stance_KL': -1.68,
                          'swing_KL' : -1.67}
    if 'stance_KV_swing_KV' in parameter_folder_str:
        output_1 = healthy_gains['stance_KV']
        output_2 = healthy_gains['swing_KV']
        
    if 'stance_KL_swing_KL' in parameter_folder_str:
        output_1 = healthy_gains['stance_KL']
        output_2 = healthy_gains['swing_KL']
        
    if 'stance_KF_swing_KF' in parameter_folder_str:
        output_1 = healthy_gains['stance_KF']
        output_2 = healthy_gains['swing_KF']
        
    return output_1, output_2


def get_healthy_gains_1d(parameter_folder):
    '''
    Returns the gains obtained from healthy optimization based off parameter
    folder path of 1d experiments.

    Parameters
    ----------
    parameter_folder : pathlib.Path
        Object defining parameter folder path.

    Returns
    -------
    float
        Healthy gain.

    '''
    parameter_folder_str = str(parameter_folder)
    simple_healthy_gains = {'KV' : {'stance' : 0.12, 'swing' : 0.13},
                            'KL' : {'stance' : 1.75, 'swing' : 0.62},
                            'KF' : {'stance' : 0.19, 'swing' : -0.24}
                           } 
    complex_healthy_gains = {'KV' : {'stance' : -2.11, 'swing' : -0.21},
                             'KL' : {'stance' : -1.68, 'swing' : -1.67}
                            } 
    
    physiological_gains = {'kshapeactive' : 0.5,
                           'max_isometric_force_ratio' : 1,
                           'optimal_fiber_length_ratio' : 1,
                           'tendon_slack_length_ratio' : 1}
    
    folder_name = parameter_folder.name
    parent_folder_name = parameter_folder.parent.name
    if parent_folder_name == 'physiological_gains':
        return physiological_gains[folder_name]

    if 'KV' in folder_name:
        parameter_name = 'KV'
    elif 'KL' in folder_name:
        parameter_name = 'KL'
    elif 'KF' in folder_name:
        parameter_name = 'KF'
        
    is_simple_model = ('solLV' in parameter_folder_str) or ('simple' in parameter_folder_str)
    
    try:
    
        if is_simple_model:
            return simple_healthy_gains[parameter_name][parent_folder_name]
        else:
            return complex_healthy_gains[parameter_name][parent_folder_name]

    except:
        print(parameter_folder)
        sys.exit(1)


def get_experiment_sto_file(experiment_folder):
    '''
    Returns the path to the sto file found under experiment_folder path

    Parameters
    ----------
    experiment_folder : string
        Experiment folder path.

    Returns
    -------
    sto_file_str : string
        Path to sto file found under experiment_folder.

    '''
    sto_files_str = [str(sto_file) for sto_file in Path(experiment_folder).rglob('*.par.sto')]
    assert len(sto_files_str) >= 1, str(experiment_folder)  + '\n' + str(sto_files_str)# check that at least one par.sto evaluation file is found
    sto_file_str =  sto_files_str[0]
    
    return sto_file_str

def assess_parameter_folder_1d(parameter_folder):
    '''
    This method assess the necessary metrics of a 1d experiment parameter folder 
    and analyzes performance in function of parameter values. It then exports
    results in a report folder.

    Parameters
    ----------
    experiment_folder : string
        Experiment folder path.


    Returns
    -------
    None.

    '''
    parameter_folder_str = str(parameter_folder)
    experiment_values = [parameter_value for parameter_value in os.listdir(parameter_folder_str) if os.path.isdir(os.path.join(parameter_folder_str, parameter_value)) and str_is_float(parameter_value) ]
    # sort experiments ascendingly
    experiment_values_idx_sorted = np.argsort([float(experiment_value) for experiment_value in experiment_values])
    experiment_values = [experiment_values[idx] for idx in experiment_values_idx_sorted]
    
    experiment_folders = [os.path.join(parameter_folder_str, str(experiment_value)) for experiment_value in experiment_values]
    experiment_sto_files = [get_experiment_sto_file(experiment_folder) for experiment_folder in experiment_folders]
    
    folder_name = parameter_folder.name # used to extract parameter name
    parent_folder_name = parameter_folder.parent.name # used to extract if stance or swing phase    
    
    experiment_scores = [experiment_sto.rstrip('.par.sto').split('_')[-1] for experiment_sto in experiment_sto_files]
    
    experiment_sto_success = []
    experiment_values_success = []
    experiment_values_float_success = []
    
    # scores experiment
    experiment_scores = []
    
    # spasticity index
    sto_files_str_error = []
    spasticity_idx_l = []
    spasticity_idx_r = []
    
    
    
    # mean ankle angles
    mean_ankle_angle_gc_l = []
    std_ankle_angle_gc_l = []
    
    mean_ankle_angle_st_l = []
    std_ankle_angle_st_l = []
    
    mean_ankle_angle_sw_l = []
    std_ankle_angle_sw_l = []
    
    side_l = 'l'
    side_name_l = 'Left'
    var_list_l = ['ankle_angle_'+side_l]
    
    
    mean_ankle_angle_gc_r = []
    std_ankle_angle_gc_r = []
    
    mean_ankle_angle_st_r = []
    std_ankle_angle_st_r = []
    
    mean_ankle_angle_sw_r = []
    std_ankle_angle_sw_r = []
    
    side_r = 'r'
    side_name_r = 'Right'
    var_list_r = ['ankle_angle_'+side_r]
    
    
    for idx in range(len(experiment_sto_files)):
        experiment_sto = experiment_sto_files[idx]
        try:
            spasticity_l = spasticity_index.all_cycles_spasticity_index(experiment_sto,'l',plot=False)
            spasticity_r = spasticity_index.all_cycles_spasticity_index(experiment_sto,'r',plot=False)
            
            mean_gc_l, std_gc_l, mean_stance_l, std_stance_l, mean_swing_l, std_swing_l = extract_sto.get_mean_gc(experiment_sto, var_list_l, side_l)
            mean_gc_r, std_gc_r, mean_stance_r, std_stance_r, mean_swing_r, std_swing_r = extract_sto.get_mean_gc(experiment_sto, var_list_r, side_r)
            
                       
            spasticity_idx_l.append(spasticity_l)
            spasticity_idx_r.append(spasticity_r)
            
            experiment_scores.append(float(experiment_sto.rstrip('.par.sto').split('_')[-1]))
            
            mean_ankle_angle_gc_l.append(mean_gc_l*180/np.pi)
            std_ankle_angle_gc_l.append(std_gc_l)
            mean_ankle_angle_st_l.append(mean_stance_l*180/np.pi)
            std_ankle_angle_st_l.append(std_stance_l)
            mean_ankle_angle_sw_l.append(mean_swing_l*180/np.pi)
            std_ankle_angle_sw_l.append(std_swing_l)
            

            mean_ankle_angle_gc_r.append(mean_gc_r*180/np.pi)
            std_ankle_angle_gc_r.append(std_gc_r)
            mean_ankle_angle_st_r.append(mean_stance_r*180/np.pi)
            std_ankle_angle_st_r.append(std_stance_r)
            mean_ankle_angle_sw_r.append(mean_swing_r*180/np.pi)
            std_ankle_angle_sw_r.append(std_swing_r)
            
            experiment_sto_success.append(experiment_sto)
            experiment_values_success.append(experiment_values[idx])
            experiment_values_float_success.append(float(experiment_values[idx]))
        except:
            sto_files_str_error.append(experiment_sto_files)
            pass
        
    experiments_dict = {'parameter_value' : experiment_values_success,
                        'score' : experiment_scores,
                        'spasticity_idx_l' : spasticity_idx_l, 
                        'mean_ankle_angle_gc_l' : mean_ankle_angle_gc_l,
                        'std_ankle_angle_gc_l' : std_ankle_angle_gc_l,
                        'mean_ankle_angle_st_l' : mean_ankle_angle_st_l,
                        'std_ankle_angle_st_l' : std_ankle_angle_st_l,
                        'mean_ankle_angle_sw_l' : mean_ankle_angle_sw_l,
                        'std_ankle_angle_sw_l' : std_ankle_angle_sw_l,
                        'spasticity_idx_r' : spasticity_idx_r,
                        'mean_ankle_angle_gc_r' : mean_ankle_angle_gc_r,
                        'std_ankle_angle_gc_r' : std_ankle_angle_gc_r,
                        'mean_ankle_angle_st_r' : mean_ankle_angle_st_r,
                        'std_ankle_angle_st_r' : std_ankle_angle_st_r,
                        'mean_ankle_angle_sw_r' : mean_ankle_angle_sw_r,
                        'std_ankle_angle_sw_r' : std_ankle_angle_sw_r}
    experiments_df = pd.DataFrame.from_dict(experiments_dict, orient='index').T
    
    
    report_folder =  os.path.join(parameter_folder_str, 'report')
    if not os.path.exists(report_folder):
        os.mkdir(report_folder)
    experiments_df.to_csv(os.path.join(report_folder,'experiment_results.csv'),sep='\t')
    
    
    # healthy metrics
    
    healthy_sto = get_healthy_sto(parameter_folder)
    score_healthy = float(healthy_sto.rstrip('.par.sto').split('_')[-1])
    spasticity_l_healthy = spasticity_index.all_cycles_spasticity_index(healthy_sto,'l',plot=False)
    spasticity_r_healthy = spasticity_index.all_cycles_spasticity_index(healthy_sto,'r',plot=False)
    mean_gc_l_healthy, std_gc_l_healthy, mean_stance_l_healthy, std_stance_l_healthy, mean_swing_l_healthy, std_swing_l_healthy = extract_sto.get_mean_gc(healthy_sto, var_list_l, side_l)
    mean_gc_r_healthy, std_gc_r_healthy, mean_stance_r_healthy, std_stance_r_healthy, mean_swing_r_healthy, std_swing_r_healthy = extract_sto.get_mean_gc(healthy_sto, var_list_r, side_r)
    total_time_healthy = par_sto_metrics.get_sto_total_time(healthy_sto)
    
    experiments_dict_healthy = {'score' : score_healthy,
                                'total_time' : total_time_healthy,
                                'spasticity_idx_l' : spasticity_l_healthy, 
                                'mean_ankle_angle_gc_l' : mean_gc_l_healthy*180/np.pi,
                                'std_ankle_angle_gc_l' : std_gc_l_healthy,
                                'mean_ankle_angle_st_l' : mean_stance_l_healthy*180/np.pi,
                                'std_ankle_angle_st_l' : std_stance_l_healthy,
                                'mean_ankle_angle_sw_l' : mean_swing_l_healthy*180/np.pi,
                                'std_ankle_angle_sw_l' : std_swing_l_healthy,
                                'spasticity_idx_r' : spasticity_r_healthy, 
                                'mean_ankle_angle_gc_r' : mean_gc_r_healthy*180/np.pi,
                                'std_ankle_angle_gc_r' : std_gc_r_healthy,
                                'mean_ankle_angle_st_r' : mean_stance_r_healthy*180/np.pi,
                                'std_ankle_angle_st_r' : std_stance_r_healthy,
                                'mean_ankle_angle_sw_r' : mean_swing_r_healthy*180/np.pi,
                                'std_ankle_angle_sw_r' : std_swing_r_healthy}
    
    experiment_title = get_experiment_title(parameter_folder)
    healthy_value = get_healthy_gains_1d(parameter_folder)
    

    # export ankle and knee angle plot
    side = 'l'
    var_names = ['ankle_angle_'+side, 'knee_angle_'+side]
    for var_name in var_names:
        ylabel = var_name.replace('_', ' ')
        par_sto_metrics.plot_mean_gc(experiment_sto_success, experiment_values_float_success,
                                     var_name, side, ylabel, experiment_title, report_folder,
                                     healthy_sto = healthy_sto, healthy_value = healthy_value)
    # export metrics plots
    
    for metric in experiments_dict:
        if metric == 'parameter_value':
            continue
        try:
            par_sto_metrics.plot_column(experiments_dict[metric], experiment_values_float_success, metric, experiment_title, report_folder, healthy_value = healthy_value, healthy_metric = experiments_dict_healthy[metric])
        except ValueError:
            print("Error in plot_column:")
            print("metric = " + metric)
            print("parameter folder = " + parameter_folder_str)
            print('Skipping...')
            pass
        
def assess_parameter_folder_2d(parameter_folder):
    '''
    This method assess the necessary metrics of a 2d experiment parameter folder 
    and analyzes performance in function of parameter values. It then exports
    results in a report folder.

    Parameters
    ----------
    experiment_folder : string
        Experiment folder path.


    Returns
    -------
    None.

    '''
    parameter_folder_str = str(parameter_folder)
    experiment_values = [parameter_value for parameter_value in os.listdir(parameter_folder_str) if os.path.isdir(os.path.join(parameter_folder_str, parameter_value)) and parameter_value != 'report' ]
    experiment_values.sort()
    
    
    
    experiment_folders = [os.path.join(parameter_folder_str, str(experiment_value)) for experiment_value in experiment_values]
    experiment_sto_files = [get_experiment_sto_file(experiment_folder) for experiment_folder in experiment_folders]
    
    experiment_scores = [experiment_sto.rstrip('.par.sto').split('_')[-1] for experiment_sto in experiment_sto_files]
    

    
    experiment_sto_success = []
    experiment_values_success = []
    experiment_scores_success = []
    
    experiment_values_stance = []
    experiment_values_swing = []
    
    experiment_total_time = []
    
    # scores experiment
    experiment_scores = []
    
    # spasticity index
    sto_files_str_error = []
    spasticity_idx_l = []
    spasticity_idx_r = []
    
    
    
    # mean ankle angles
    mean_ankle_angle_gc_l = []
    std_ankle_angle_gc_l = []
    
    mean_ankle_angle_st_l = []
    std_ankle_angle_st_l = []
    
    mean_ankle_angle_sw_l = []
    std_ankle_angle_sw_l = []
    
    side_l = 'l'
    side_name_l = 'Left'
    var_list_l = ['ankle_angle_'+side_l]
    
    
    mean_ankle_angle_gc_r = []
    std_ankle_angle_gc_r = []
    
    mean_ankle_angle_st_r = []
    std_ankle_angle_st_r = []
    
    mean_ankle_angle_sw_r = []
    std_ankle_angle_sw_r = []
    
    side_r = 'r'
    side_name_r = 'Right'
    var_list_r = ['ankle_angle_'+side_r]
    
    
    for idx in range(len(experiment_sto_files)):
        experiment_sto = experiment_sto_files[idx]
        try:
            spasticity_l = spasticity_index.all_cycles_spasticity_index(experiment_sto,'l',plot=False)
        except:
            spasticity_l = np.nan
        try:
            spasticity_r = spasticity_index.all_cycles_spasticity_index(experiment_sto,'r',plot=False)
        except:
            spasticity_r = np.nan
            
        try:
            mean_gc_l, std_gc_l, mean_stance_l, std_stance_l, mean_swing_l, std_swing_l = extract_sto.get_mean_gc(experiment_sto, var_list_l, side_l)
        except:
            mean_gc_l = np.nan
            std_gc_l = np.nan
            mean_stance_l = np.nan
            std_stance_l = np.nan
            mean_swing_l = np.nan
            std_swing_l = np.nan
        try:
            mean_gc_r, std_gc_r, mean_stance_r, std_stance_r, mean_swing_r, std_swing_r = extract_sto.get_mean_gc(experiment_sto, var_list_r, side_r)
        except:
            mean_gc_r = np.nan
            std_gc_r = np.nan
            mean_stance_r = np.nan
            std_stance_r = np.nan
            mean_swing_r = np.nan
            std_swing_r = np.nan
            
            
        spasticity_idx_l.append(spasticity_l)
        spasticity_idx_r.append(spasticity_r)
        
        experiment_scores.append(float(experiment_sto.rstrip('.par.sto').split('_')[-1]))
        
        mean_ankle_angle_gc_l.append(mean_gc_l*180/np.pi)
        std_ankle_angle_gc_l.append(std_gc_l)
        mean_ankle_angle_st_l.append(mean_stance_l*180/np.pi)
        std_ankle_angle_st_l.append(std_stance_l)
        mean_ankle_angle_sw_l.append(mean_swing_l*180/np.pi)
        std_ankle_angle_sw_l.append(std_swing_l)
        
        
        mean_ankle_angle_gc_r.append(mean_gc_r*180/np.pi)
        std_ankle_angle_gc_r.append(std_gc_r)
        mean_ankle_angle_st_r.append(mean_stance_r*180/np.pi)
        std_ankle_angle_st_r.append(std_stance_r)
        mean_ankle_angle_sw_r.append(mean_swing_r*180/np.pi)
        std_ankle_angle_sw_r.append(std_swing_r)
        
        experiment_sto_success.append(experiment_sto)
        experiment_values_success.append(experiment_values[idx])
        
        experiment_values_stance.append(float(experiment_values[idx].split('_')[0]))
        
        experiment_values_swing.append(float(experiment_values[idx].split('_')[1]))
        
        
        experiment_scores_success.append(experiment_scores[idx])
        try:
            experiment_total_time.append(par_sto_metrics.get_sto_total_time(experiment_sto))
        except:
            print(str(experiment_sto))
            experiment_total_time.append(np.nan)
        
    experiments_dict = {'score' : experiment_scores_success,
                        'total_time' : experiment_total_time,
                        'spasticity_idx_l' : spasticity_idx_l, 
                        'mean_ankle_angle_gc_l' : mean_ankle_angle_gc_l,
                        'std_ankle_angle_gc_l' : std_ankle_angle_gc_l,
                        'mean_ankle_angle_st_l' : mean_ankle_angle_st_l,
                        'std_ankle_angle_st_l' : std_ankle_angle_st_l,
                        'mean_ankle_angle_sw_l' : mean_ankle_angle_sw_l,
                        'std_ankle_angle_sw_l' : std_ankle_angle_sw_l,
                        'spasticity_idx_r' : spasticity_idx_r,
                        'mean_ankle_angle_gc_r' : mean_ankle_angle_gc_r,
                        'std_ankle_angle_gc_r' : std_ankle_angle_gc_r,
                        'mean_ankle_angle_st_r' : mean_ankle_angle_st_r,
                        'std_ankle_angle_st_r' : std_ankle_angle_st_r,
                        'mean_ankle_angle_sw_r' : mean_ankle_angle_sw_r,
                        'std_ankle_angle_sw_r' : std_ankle_angle_sw_r}
    
    
    
    experiment_values_tuple = list(zip(*[experiment_values_stance, experiment_values_swing]))
    df_index = pd.MultiIndex.from_tuples(experiment_values_tuple, names=['stance', 'swing'])
    experiments_df = pd.DataFrame(experiments_dict, index=df_index)
    
    
    # healthy metrics
    stance_gain_healthy, swing_gain_healthy = get_healthy_gains_2d(parameter_folder)
    healthy_sto = get_healthy_sto(parameter_folder)
    score_healthy = float(healthy_sto.rstrip('.par.sto').split('_')[-1])
    spasticity_l_healthy = spasticity_index.all_cycles_spasticity_index(healthy_sto,'l',plot=False)
    spasticity_r_healthy = spasticity_index.all_cycles_spasticity_index(healthy_sto,'r',plot=False)
    mean_gc_l_healthy, std_gc_l_healthy, mean_stance_l_healthy, std_stance_l_healthy, mean_swing_l_healthy, std_swing_l_healthy = extract_sto.get_mean_gc(healthy_sto, var_list_l, side_l)
    mean_gc_r_healthy, std_gc_r_healthy, mean_stance_r_healthy, std_stance_r_healthy, mean_swing_r_healthy, std_swing_r_healthy = extract_sto.get_mean_gc(healthy_sto, var_list_r, side_r)
    total_time_healthy = par_sto_metrics.get_sto_total_time(healthy_sto)
    
    experiments_dict_healthy = {'score' : score_healthy,
                                'total_time' : total_time_healthy,
                                'spasticity_idx_l' : spasticity_l_healthy, 
                                'mean_ankle_angle_gc_l' : mean_gc_l_healthy*180/np.pi,
                                'std_ankle_angle_gc_l' : std_gc_l_healthy,
                                'mean_ankle_angle_st_l' : mean_stance_l_healthy*180/np.pi,
                                'std_ankle_angle_st_l' : std_stance_l_healthy,
                                'mean_ankle_angle_sw_l' : mean_swing_l_healthy*180/np.pi,
                                'std_ankle_angle_sw_l' : std_swing_l_healthy,
                                'spasticity_idx_r' : spasticity_r_healthy, 
                                'mean_ankle_angle_gc_r' : mean_gc_r_healthy*180/np.pi,
                                'std_ankle_angle_gc_r' : std_gc_r_healthy,
                                'mean_ankle_angle_st_r' : mean_stance_r_healthy*180/np.pi,
                                'std_ankle_angle_st_r' : std_stance_r_healthy,
                                'mean_ankle_angle_sw_r' : mean_swing_r_healthy*180/np.pi,
                                'std_ankle_angle_sw_r' : std_swing_r_healthy}
    
    
    report_folder =  os.path.join(parameter_folder_str, 'report')
    if not os.path.exists(report_folder):
        os.mkdir(report_folder)
    experiments_df.to_csv(os.path.join(report_folder,'experiment_results.csv'),sep='\t')
    
    # export metrics plots
    
    folder_name = parameter_folder.name
    folder_name_split = folder_name.split('_')
    
    for metric in experiments_dict:
        df2 = experiments_df[metric].reset_index().pivot(columns='swing',index='stance', values=metric)
        df2.interpolate(method='akima', axis=1 ,inplace=True, limit_direction='both')
        fig, ax = plt.subplots()
        
        ax = sns.heatmap(df2, cmap="RdBu_r" ,center=experiments_dict_healthy[metric])
        ax.set(title=metric.replace('_',' ') + ', healthy in green.',
               xlabel=' '.join(folder_name_split[2:3]),
               ylabel=' '.join(folder_name_split[0:1]))
            
        plt.plot(swing_gain_healthy/0.4 + 1.5,stance_gain_healthy/0.2 + 0.5,label='healthy',marker='P', color ='lime',markersize=20)
        
        plt.title(metric)
        plt.savefig(os.path.join(report_folder,metric+'.pdf'))

# Main script
if __name__ == "__main__":
    # Set this to True to analyze 1d experiment results
    one_d = True
    # Set this to True to analyze 2d experiement results
    two_d = False
    
    
    results_folder = sys.argv[1] # path of folder containing optim_results
    results_folder_1d = os.path.join(results_folder,'1d')
    results_folder_2d = os.path.join(results_folder,'2d')
    
    # extract folders for each 1d optimization
    if one_d:
        sto_files_1d =  Path(results_folder_1d).rglob('*.par.sto')
        parameter_folders_1d = [sto_file.parent.parent.parent for sto_file in sto_files_1d]
        parameter_folders_1d = list(set(parameter_folders_1d)) # extract unique elements
        for parameter_folder in parameter_folders_1d:
            assess_parameter_folder_1d(parameter_folder)

    # extract folders for each 2d optimization
    if two_d:
        sto_files_2d =  Path(results_folder_2d).rglob('*.par.sto')
        parameter_folders_2d = [sto_file.parent.parent.parent for sto_file in sto_files_2d]
        parameter_folders_2d = list(set(parameter_folders_2d)) # extract unique elements
        for parameter_folder in parameter_folders_2d:
            assess_parameter_folder_2d(parameter_folder)
        
    