"""
This script contains multiple functions to evaluate various metrics for each of the parameters explored and
generates a report folder containing plots of these metrics.
"""
import sys, os
from pathlib import Path
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import par_sto_metrics as psm
import spasticity_index, extract_sto


def get_experiment_title(parameter_folder):
    '''Returns an experiment title based off the parameter folder and
    its parent folder.

    Parameters
    ----------
    parameter_folder: (pathlib.Path) object defining parameter folder path.

    Returns
    -------
    string: title of the experiment.
    '''
    parameter_phase = parameter_folder.parent.name
    parameter_name = parameter_folder.name
    parameter_phase_name = parameter_phase + '_' + parameter_name
    
    return parameter_phase_name.replace('_', ' ')


def get_score_list(source_folder):
    '''Returns the child folder list of a source folder and sto scores corresponding
    to each child folder.

    Parameters
    ----------
    source_folder: (string) path of source folder.

    Returns
    -------
    child_folders: (list) list containing child folder names.
    sto_scores: (list) list containing scores corresponding to each child folder.
    '''
    child_folders = [child_folder for child_folder in os.listdir(source_folder) if os.path.isdir(os.path.join(source_folder, child_folder))]
    sto_files = [get_experiment_sto_file(os.path.join(source_folder, child_folder)) for child_folder in child_folders]

    sto_scores = [sto_file.rstrip('.par.sto').split('_')[-1] for sto_file in sto_files]

    return child_folders, sto_scores
    

def str_is_float(x):
    '''Returns whether a string can be converted to a float or not.

    Parameters
    ----------
    x: (string) potential float.

    Returns
    -------
    bool: True if string can be converted to a float. False otherwise
    '''
    try:
        float(x)
        return True
    except ValueError:
        return False


def get_healthy_sto():
    '''Returns path to healthy sto file based off parameter folder name.

    Parameters
    ----------
    parameter_folder: (pathlib.Path) object defining parameter folder path.

    Returns
    -------
    string: path to healthy sto file.
    '''

    return '../models/optimisation/optimisation_0.6.sto'


def get_experiment_sto_file(experiment_folder):
    '''Returns the path to the sto file found under experiment_folder path

    Parameters
    ----------
    experiment_folder: (string) experiment folder path.

    Returns
    -------
    sto_file_str: (string) path to sto file found under experiment_folder.
    '''
    sto_files_str = [str(sto_file) for sto_file in Path(experiment_folder).rglob('*.par.sto')]
    assert len(sto_files_str) >= 1, str(experiment_folder) + '\n' + str(sto_files_str)
    sto_file_str = sto_files_str[0]
    
    return sto_file_str


def assess_parameter_folder_1d(parameter_folder, side):
    '''Assesses the necessary metrics of a 1d experiment parameter folder
    and analyzes performance in function of parameter values. It then exports
    results in a report folder.

    Parameters
    ----------
    experiment_folder: (string) experiment folder path.
    side: (sting) 'l' or 'r' for left or right leg

    Returns
    -------
    None.
    '''
    parameter_folder_str = str(parameter_folder)
    param = parameter_folder_str.split("\\")[-1]
    param = param.split('.')[0]

    experiment_values = [parameter_value.split('_')[-1] for parameter_value in os.listdir(parameter_folder_str)
                         if os.path.isdir(os.path.join(parameter_folder_str, parameter_value)) and
                         len(parameter_value.split('.')) > 1]
    experiment_values = [parameter_value.split('.')[0] for parameter_value in experiment_values]

    # sort experiments ascendingly
    experiment_values_idx_sorted = np.argsort([float(experiment_value) for experiment_value in
                                               experiment_values])
    experiment_values = [experiment_values[idx] for idx in experiment_values_idx_sorted]

    experiment_folders = [os.path.join(parameter_folder_str, param+'_'+str(experiment_value) +
                                       '.f0914m.GH2010v8.S10W.D15.I') for experiment_value in experiment_values]
    experiment_sto_files = [get_experiment_sto_file(experiment_folder) for experiment_folder in experiment_folders]

    experiment_sto_success = []
    experiment_values_success = []
    experiment_values_float_success = []
    
    # scores
    experiment_scores = []
    total_times = []
    
    # ankle angle
    var_name = 'ankle_angle_' + side

    # spasticity indexes
    spasticity_idx_sol = []
    spasticity_std_sol = []
    spasticity_idx_gas = []
    spasticity_std_gas = []

    # altered models metrics
    for idx in range(len(experiment_sto_files)):
        experiment_sto = experiment_sto_files[idx]
        print(experiment_sto)

        # score
        experiment_scores.append(float(experiment_sto.rstrip('.par.sto').split('_')[-1]))
        total_times.append(psm.get_sto_total_time(experiment_sto))

        # spasticity indexes
        spas_idx_sol, spas_std_sol = spasticity_index.all_cycles_spasticity_index(experiment_sto, "soleus", side,
                                                                                plot=False)
        spas_idx_gas, spas_std_gas = spasticity_index.all_cycles_spasticity_index(experiment_sto, "gastroc",
                                                                                        side, plot=False)
        spasticity_idx_sol.append(spas_idx_sol)
        spasticity_std_sol.append(spas_std_sol)
        spasticity_idx_gas.append(spas_idx_gas)
        spasticity_std_gas.append(spas_std_gas)

        experiment_sto_success.append(experiment_sto)
        experiment_values_success.append(experiment_values[idx])
        experiment_values_float_success.append(float(experiment_values[idx]))
    
    experiments_dict = {'parameter_value': experiment_values_success,
                        'score': experiment_scores,
                        'total_time': total_times,
                        'spasticity_index_sol': spasticity_idx_sol,
                        'spasticity_std_sol': spasticity_std_sol,
                        'spasticity_index_gas': spasticity_idx_gas,
                        'spasticity_std_gas': spasticity_std_gas}

    experiments_df = pd.DataFrame.from_dict(experiments_dict, orient='index').T
    
    report_folder = os.path.join(parameter_folder_str, 'report')
    if not os.path.exists(report_folder):
        os.mkdir(report_folder)
    experiments_df.to_csv(os.path.join(report_folder, 'experiment_results.csv'), sep='\t')
    
    # healthy metrics
    healthy_sto = get_healthy_sto()
    score_healthy = float(healthy_sto.rstrip('.par.sto').split('_')[-1])
    total_time_healthy = psm.get_sto_total_time(healthy_sto)

    spas_idx_sol_h, spas_std_sol_h = spasticity_index.all_cycles_spasticity_index(healthy_sto, "soleus",
                                                                                      side, plot=False)
    spas_idx_gas_h, spas_std_gas_h = spasticity_index.all_cycles_spasticity_index(healthy_sto, "gastroc",
                                                                                      side, plot=False)

    experiments_dict_healthy = {'score': score_healthy,
                                'total_time': total_time_healthy,
                                'spasticity_index_sol': spas_idx_sol_h,
                                'spasticity_std_sol': spas_std_sol_h,
                                'spasticity_index_gas': spas_idx_gas_h,
                                'spasticity_std_gas': spas_std_gas_h}
    
    experiment_title = get_experiment_title(parameter_folder)
    healthy_value = 100

    # export ankle plot
    ylabel = ""
    for s in range(len(var_name.split('_'))-1):
        ylabel = ylabel + " " + var_name.split('_')[s]
    ylabel = ylabel + " (" + var_name.split('_')[-1] + ")"
    title = ""
    inv = False
    if experiment_title.split(" ")[0] == "biomechanical":
        inv = True
        for n in range(len(experiment_title.split(" ")[2:])):
            title = title + " " + experiment_title.split(" ")[2+n]
    else:
        s = 1
        if "weak" in experiment_title.split(" "):
            inv = True
            for n in range(len(experiment_title.split(" "))-1):
                if experiment_title.split(" ")[s+n] != "weak":
                    title = title + " " + experiment_title.split(" ")[s+n]
        else:
            for n in range(len(experiment_title.split(" "))-1):
                title = title + " " + experiment_title.split(" ")[s+n]

    psm.plot_mean_gc(experiment_sto_success, experiment_values_float_success, var_name, side, ylabel, title,
                     report_folder, healthy_sto=healthy_sto, healthy_value=healthy_value, inv=inv)

    # export metrics plots
    for metric in experiments_dict:
        if metric.split("_")[0] in ['score', 'total']:
            title = ""
            inv = False
            if experiment_title.split(" ")[0] == "biomechanical":
                inv = True
                for n in range(len(experiment_title.split(" ")[2:])):
                    title = title + " " + experiment_title.split(" ")[2 + n]
            else:
                s = 1
                if "weak" in experiment_title.split(" "):
                    inv = True
                    for n in range(len(experiment_title.split(" ")) - 1):
                        if experiment_title.split(" ")[s + n] != "weak":
                            title = title + " " + experiment_title.split(" ")[s + n]
                else:
                    for n in range(len(experiment_title.split(" ")) - 1):
                        title = title + " " + experiment_title.split(" ")[s + n]
            psm.plot_column(experiments_dict[metric], experiment_values_float_success, metric, title,
                                report_folder, healthy_value, experiments_dict_healthy[metric], inv=inv)

        elif metric.split("_")[0] == "spasticity" and metric.split("_")[-1] == 'sol':
            title = ""
            inv = False
            if experiment_title.split(" ")[0] == "biomechanical":
                inv = True
                for n in range(len(experiment_title.split(" ")[2:])):
                    title = title + " " + experiment_title.split(" ")[2+n]
            else:
                s = 1
                if "weak" in experiment_title.split(" "):
                    inv = True
                    for n in range(len(experiment_title.split(" "))-1):
                        if experiment_title.split(" ")[s+n] != "weak":
                            title = title + " " + experiment_title.split(" ")[s+n]
                else:
                    for n in range(len(experiment_title.split(" "))-1):
                        title = title + " " + experiment_title.split(" ")[s+n]
            metric_1 = 'spasticity_index_sol'
            std_1 = "spasticity_std_sol"
            metric_2 = 'spasticity_index_gas'
            std_2 = "spasticity_std_gas"
            psm.plot_column_std_two_muscles(experiments_dict[metric_1], experiments_dict[std_1],
                                            experiments_dict[metric_2], experiments_dict[std_2],
                                            experiment_values_float_success, 'spasticity_indexes', title,
                                            report_folder, healthy_value, experiments_dict_healthy[metric_1],
                                            experiments_dict_healthy[metric_2], experiments_dict_healthy[std_1],
                                            experiments_dict_healthy[std_2], side, inv=inv)


def assess_parameter_folder_2d(parameter_folder, side):
    '''Assesses the necessary metrics of a 2d experiment parameter folder
    and analyzes performance in function of parameter values. It then exports
    results in a report folder.

    Parameters
    ----------
    experiment_folder: (string) experiment folder path.
    side: (string) 'l' or 'r' for left or right leg

    Returns
    -------
    None.
    '''
    parameter_folder_str = str(parameter_folder)
    param = parameter_folder_str.split("\\")[-1]
    param = param.split('.')[0]

    experiment_values = [parameter_value.split('_')[-2] + '_' + parameter_value.split('_')[-1] for parameter_value in
                         os.listdir(parameter_folder_str) if
                         os.path.isdir(os.path.join(parameter_folder_str, parameter_value)) and
                         len(parameter_value.split('.')) > 1]
    experiment_values = [parameter_value.split('.')[0] for parameter_value in experiment_values]

    experiment_folders = [os.path.join(parameter_folder_str, param + '_' + str(experiment_value.split('_')[0]) +
                                       '_' + str(experiment_value.split('_')[1]) + '.f0914m.GH2010v8.S10W.D15.I')
                          for experiment_value in experiment_values]
    experiment_sto_files = [get_experiment_sto_file(experiment_folder) for experiment_folder in experiment_folders]

    experiment_sto_success = []
    experiment_values_success = []
    experiment_scores_success = []
    
    experiment_values_stance = []
    experiment_values_swing = []
    
    # scores experiment
    experiment_scores = []
    experiment_total_time = []
    
    # ankle angle
    var_name = 'ankle_angle_' + side

    ME_gc = []
    ME_st = []
    ME_sw = []

    # spasticity index
    spasticity_idx_sol = []
    spasticity_idx_gas = []

    healthy_sto = get_healthy_sto()
    
    for idx in range(len(experiment_sto_files)):
        experiment_sto = experiment_sto_files[idx]
        print(experiment_sto)
        # scores
        experiment_scores.append(float(experiment_sto.rstrip('.par.sto').split('_')[-1]))

        # ME
        try:
            me_gc, me_std_gc = extract_sto.me(experiment_sto, healthy_sto, var_name, side)
            me_st, me_std_st, me_sw, me_std_sw = extract_sto.me_phases(experiment_sto, healthy_sto, var_name, side)
        except:
            me_gc, me_std_gc = np.nan, np.nan
            me_st, me_std_st = np.nan, np.nan
            me_sw, me_std_sw = np.nan, np.nan

        ME_gc.append(me_gc)
        ME_st.append(me_st)
        ME_sw.append(me_sw)

        # spasticity indexes
        try:
            spas_idx_sol, spas_std_sol = spasticity_index.all_cycles_spasticity_index(experiment_sto, "soleus",
                                                                                      side, plot=False)
        except:
            spas_idx_sol, spas_std_sol = np.nan, np.nan
        try:
            spas_idx_gas, spas_std_gas = spasticity_index.all_cycles_spasticity_index(experiment_sto, "gastroc",
                                                                                      side, plot=False)
        except:
            spas_idx_gas, spas_std_gas = np.nan, np.nan

        spasticity_idx_sol.append(spas_idx_sol)
        spasticity_idx_gas.append(spas_idx_gas)
        
        experiment_sto_success.append(experiment_sto)
        experiment_values_success.append(experiment_values[idx])
        experiment_values_stance.append(int(float(experiment_values[idx].split('_')[0])))
        
        experiment_values_swing.append(int(float(experiment_values[idx].split('_')[1])))
        
        experiment_scores_success.append(experiment_scores[idx])
        try:
            experiment_total_time.append(psm.get_sto_total_time(experiment_sto))
        except:
            experiment_total_time.append(np.nan)
        
    experiments_dict = {'score': experiment_scores_success,
                        'total_time': experiment_total_time,
                        'ME_gc': ME_gc,
                        'ME_st': ME_st,
                        'ME_sw': ME_sw,
                        'spasticity_index_sol': spasticity_idx_sol,
                        'spasticity_index_gas': spasticity_idx_gas}

    experiment_values_tuple = list(zip(*[experiment_values_stance, experiment_values_swing]))
    df_index = pd.MultiIndex.from_tuples(experiment_values_tuple, names=['stance', 'swing'])
    experiments_df = pd.DataFrame(experiments_dict, index=df_index)
    
    # healthy metrics
    score_healthy = float(healthy_sto.rstrip('.par.sto').split('_')[-1])
    total_time_healthy = psm.get_sto_total_time(healthy_sto)

    ME_gc_h = 0
    ME_st_h = 0
    ME_sw_h = 0

    spas_idx_sol_h, spas_std_sol_h = spasticity_index.all_cycles_spasticity_index(healthy_sto, "soleus",
                                                                                  side, plot=False)
    spas_idx_gas_h, spas_std_gas_h = spasticity_index.all_cycles_spasticity_index(healthy_sto, "gastroc",
                                                                                  side, plot=False)
    experiments_dict_healthy = {'score': score_healthy,
                                'total_time': total_time_healthy,
                                'ME_gc': ME_gc_h,
                                'ME_st': ME_st_h,
                                'ME_sw': ME_sw_h,
                                'spasticity_index_sol': spas_idx_sol_h,
                                'spasticity_index_gas': spas_idx_gas_h}
    
    report_folder = os.path.join(parameter_folder_str, 'report')
    if not os.path.exists(report_folder):
        os.mkdir(report_folder)
    experiments_df.to_csv(os.path.join(report_folder, 'experiment_results.csv'), sep='\t')
    
    # export metrics plots
    folder_name = parameter_folder.name
    folder_name_split = folder_name.split('_')
    
    for metric in experiments_dict:
        df2 = experiments_df[metric].reset_index().pivot(columns='swing', index='stance', values=metric)
        df2.interpolate(method='akima', axis=1, inplace=True, limit_direction='both')
        if folder_name in ["stance_KF_stance_KL", "stance_TA_KS_swing_TA_KS"]:
            df2 = df2[df2.columns[::-1]]
        else:
            df2 = df2[::-1]
        
        map = "RdBu"
        if metric in ["ME_gc", "ME_st", "ME_sw"]:
            if folder_name == "stance_KF_stance_KL":
                ax = sns.heatmap(df2, cmap=map, vmin=-4, vmax=8, center=experiments_dict_healthy[metric])
                ax.set(xlabel='stance KL [%]', ylabel='stance KF [%]')
            elif folder_name == "stance_TA_KS_swing_TA_KS":
                ax = sns.heatmap(df2, cmap=map, vmin=-8, vmax=4, center=experiments_dict_healthy[metric])
                ax.set(xlabel='swing TA KS [%]', ylabel='stance TA KS [%]')
            else:
                ax = sns.heatmap(df2, cmap=map, vmin=-8, vmax=4, center=experiments_dict_healthy[metric])
                ax.set(xlabel=folder_name_split[2]+" "+folder_name_split[3]+' [%]',
                       ylabel=folder_name_split[0]+" "+folder_name_split[1]+' [%]')
        else:
            if folder_name == "stance_KF_stance_KL":
                ax = sns.heatmap(df2, cmap=map, center=experiments_dict_healthy[metric])
                ax.set(xlabel='stance KL [%]', ylabel='stance KF [%]')
            elif folder_name == "stance_TA_KS_swing_TA_KS":
                ax = sns.heatmap(df2, cmap=map+"_r", center=experiments_dict_healthy[metric])
                ax.set(xlabel='swing TA KS [%]', ylabel='stance TA KS [%]')
            else:
                ax = sns.heatmap(df2, cmap=map+"_r", center=experiments_dict_healthy[metric])
                ax.set(xlabel=folder_name_split[2]+" "+folder_name_split[3]+' [%]',
                       ylabel=folder_name_split[0]+" "+folder_name_split[1]+' [%]')

        if folder_name == 'stance_KS_swing_KS':
            par = 'increased stance and swing KS'
        elif folder_name == 'stance_TA_KS_swing_TA_KS':
            par = 'increased stance and swing TA KS'
        elif folder_name == 'stance_KL_swing_KL':
            par = 'increased stance and swing KL'
        elif folder_name == 'stance_KF_stance_KL':
            par = 'decreased stance KF and KL'
        if metric == 'spasticity_index_sol':
            title = 'mean SOL spasticity index ('+side+') \n for '+par
        elif metric == 'spasticity_index_gas':
            title = 'mean GAS spasticity index ('+side+') \n for '+par
        elif metric == 'ME_gc':
            title = 'mean ME ('+side+') over the gait cycle \n for '+par
        elif metric == 'ME_st':
            title = 'mean ME ('+side+') over the stance phase \n for '+par
        elif metric == 'ME_sw':
            title = 'mean ME ('+side+') over the swing phase \n for '+par
        else:
            title = metric
        plt.title(title)
        plt.savefig(os.path.join(report_folder, metric+'.png'))
        plt.close()
