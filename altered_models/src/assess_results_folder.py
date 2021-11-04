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
    bool: True if string can be converted to a float. False otherwise.
    '''
    try:
        float(x)
        return True
    except ValueError:
        return False


def get_healthy_sto():
    '''Returns path to healthy sto file based off parameter folder name.

    Returns
    -------
    string: path to healthy sto file.
    '''

    return '../models/optimisation/optimisation_0.6.sto'


def get_experiment_sto_file(parameter_folder):
    '''Returns the path to the sto file found under experiment_folder path.
    Parameters
    ----------
    parameter_folder: (string) experiment folder path.

    Returns
    -------
    sto_file_str: (string) path to sto file found under experiment_folder.
    '''
    sto_files_str = [str(sto_file) for sto_file in Path(parameter_folder).rglob('*.par.sto')]
    assert len(sto_files_str) >= 1, str(parameter_folder) + '\n' + str(sto_files_str)
    sto_file_str = sto_files_str[0]
    
    return sto_file_str


def assess_parameter_folder_1d(parameter_folder, side):
    '''Assesses the necessary metrics of a 1d experiment parameter folder
    and analyzes performance in function of parameter values. It then exports
    results in a report folder.
    Parameters
    ----------
    parameter_folder: (string) experiment folder path.
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
    var_names = ['ankle_angle_' + side, 'knee_angle_' + side, 'hip_flexion_' + side]

    # spasticity indexes
    spasticity_idx_sol = []
    spasticity_std_sol = []
    spasticity_idx_gas = []
    spasticity_std_gas = []

    # gait features
    mstep_length = []
    sstep_length = []
    gspeed = []
    mstance_period = []
    sstance_period = []

    # altered models metrics
    for idx in range(len(experiment_sto_files)):
        experiment_sto = experiment_sto_files[idx]
        print(experiment_sto)

        # score
        experiment_scores.append(float(experiment_sto.rstrip('.par.sto').split('_')[-1]))
        total_times.append(psm.get_sto_total_time(experiment_sto))

        # gait features
        mstep_l, sstep_l, speed = extract_sto.gait_features(experiment_sto)
        mstep_length.append(mstep_l)
        sstep_length.append(sstep_l)
        gspeed.append(speed)

        mstance_p, sstance_p = extract_sto.stance_period(experiment_sto, 'l')
        mstance_period.append(mstance_p)
        sstance_period.append(sstance_p)

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
                        'spasticity_std_gas': spasticity_std_gas,
                        'mstep_length': mstep_length,
                        'sstep_length': sstep_length,
                        'speed': gspeed,
                        'mstance_period': mstance_period,
                        'sstance_period': sstance_period}

    experiments_df = pd.DataFrame.from_dict(experiments_dict, orient='index').T
    
    report_folder = os.path.join(parameter_folder_str, 'report')
    if not os.path.exists(report_folder):
        os.mkdir(report_folder)
    experiments_df.to_csv(os.path.join(report_folder, 'experiment_results.csv'), sep='\t')
    
    # healthy metrics
    healthy_sto = get_healthy_sto()
    score_healthy = float(healthy_sto.rstrip('.par.sto').split('_')[-1])
    total_time_healthy = psm.get_sto_total_time(healthy_sto)

    h_mstep_length, h_sstep_length, h_gspeed = extract_sto.gait_features(healthy_sto)
    h_mstance_period, h_sstance_period = extract_sto.stance_period(healthy_sto, 'l')

    spas_idx_sol_h, spas_std_sol_h = spasticity_index.all_cycles_spasticity_index(healthy_sto, "soleus",
                                                                                      side, plot=False)
    spas_idx_gas_h, spas_std_gas_h = spasticity_index.all_cycles_spasticity_index(healthy_sto, "gastroc",
                                                                                      side, plot=False)

    experiments_dict_healthy = {'score': score_healthy,
                                'total_time': total_time_healthy,
                                'spasticity_index_sol': spas_idx_sol_h,
                                'spasticity_std_sol': spas_std_sol_h,
                                'spasticity_index_gas': spas_idx_gas_h,
                                'spasticity_std_gas': spas_std_gas_h,
                                'mstep_length': h_mstep_length,
                                'sstep_length': h_sstep_length,
                                'speed': h_gspeed,
                                'mstance_period': h_mstance_period,
                                'sstance_period': h_sstance_period}
    
    experiment_title = get_experiment_title(parameter_folder)
    healthy_value = 100

    # export ankle plot
    title = ""
    inv = False
    es = 'es'
    if experiment_title.split(" ")[0] == "biomechanical":
        inv = True
        for n in range(len(experiment_title.split(" ")[1:])):
            title = title + " " + experiment_title.split(" ")[1+n]
        if "force" in experiment_title.split(" "):
            es = "no"
    else:
        s = 1
        if "weak" in experiment_title.split(" "):
            inv = True
            es = "no"
            for n in range(len(experiment_title.split(" "))-1):
                if experiment_title.split(" ")[s+n] != "weak":
                    title = title + " " + experiment_title.split(" ")[s+n]
        else:
            for n in range(len(experiment_title.split(" "))-1):
                title = title + " " + experiment_title.split(" ")[s+n]

    for var_name in var_names:
        psm.plot_mean_gc(experiment_sto_success, experiment_values_float_success, var_name, side, title,
                         report_folder, healthy_sto=healthy_sto, healthy_value=healthy_value, inv=inv, es=es)

    # export metrics plots
    for metric in experiments_dict:
        if metric.split("_")[0] in ['score', 'total', 'speed']:
            title = ""
            inv = False
            if experiment_title.split(" ")[0] == "biomechanical":
                inv = True
                for n in range(len(experiment_title.split(" ")[1:])):
                    title = title + " " + experiment_title.split(" ")[1 + n]
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

        elif metric.split("_")[0] == "mstance":
            title = ""
            inv = False
            if experiment_title.split(" ")[0] == "biomechanical":
                inv = True
                for n in range(len(experiment_title.split(" ")[1:])):
                    title = title + " " + experiment_title.split(" ")[1 + n]
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
            psm.plot_columns_std(experiments_dict[metric], experiments_dict['sstance_period'],
                                 experiment_values_float_success, metric, title, report_folder, healthy_value,
                                 experiments_dict_healthy[metric], experiments_dict_healthy['sstance_period'],
                                 column_values2=experiments_dict['mstep_length'],
                                 std_values2=experiments_dict['sstep_length'],
                                 healthy_metric2=experiments_dict_healthy['mstep_length'],
                                 std_healthy2=experiments_dict_healthy['sstep_length'], inv=inv)
        elif "biomechanical" not in experiment_title.split(" ") and "weak" not in experiment_title.split(" ") \
                and metric.split("_")[0] == "spasticity" and metric.split("_")[-1] == 'sol':
            title = ""
            inv = False
            for n in range(len(experiment_title.split(" "))-1):
                title = title + " " + experiment_title.split(" ")[s+n]
            psm.plot_columns_std(experiments_dict[metric], experiments_dict["spasticity_std_sol"],
                                 experiment_values_float_success, 'spasticity_indexes_l', title, report_folder,
                                 healthy_value, experiments_dict_healthy[metric],
                                 experiments_dict_healthy["spasticity_std_sol"],
                                 column_values2=experiments_dict['spasticity_index_gas'],
                                 std_values2=experiments_dict["spasticity_std_gas"],
                                 healthy_metric2=experiments_dict_healthy['spasticity_index_gas'],
                                 std_healthy2=experiments_dict_healthy['spasticity_std_gas'], inv=inv)


def assess_parameter_folder_2d(parameter_folder, side):
    '''Assesses the necessary metrics of a 2d experiment parameter folder
    and analyzes performance in function of parameter values. It then exports
    results in a report folder.
    Parameters
    ----------
    parameter_folder: (string) experiment folder path.
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
    ME_st = []
    min_es = []
    moment_peaks = []
    moment_max = []

    # gait features
    stance_period = []
    step_length = []
    speed = []

    # spasticity indexes
    spasticity_idx_sol = []
    spasticity_idx_gas = []

    healthy_sto = get_healthy_sto()
    # healthy metrics
    score_healthy = float(healthy_sto.rstrip('.par.sto').split('_')[-1])
    total_time_healthy = psm.get_sto_total_time(healthy_sto)

    me_h, _, min_es_h, min_std_es = extract_sto.mean_stance(healthy_sto, var_name,
                                                            side)  # extract_sto.min_ankle_es(healthy_sto, var_name, side)
    step_length_h, sstep_length, speed_h = extract_sto.gait_features(healthy_sto)
    stance_period_h, sstance_period = extract_sto.stance_period(healthy_sto, side)

    moment_norm = extract_sto.norm_moment(healthy_sto, var_name, side)
    var_names, var_tab = extract_sto.extract_sto(healthy_sto)
    mdelta_h, sdelta_h, mmax_h, smax_h = extract_sto.moment_peaks(var_names, var_tab, side, moment_norm)

    spas_idx_sol_h, spas_std_sol_h = spasticity_index.all_cycles_spasticity_index(healthy_sto, "soleus",
                                                                                  side, plot=False)
    spas_idx_gas_h, spas_std_gas_h = spasticity_index.all_cycles_spasticity_index(healthy_sto, "gastroc",
                                                                                  side, plot=False)
    experiments_dict_healthy = {'score': score_healthy,
                                'total_time': total_time_healthy,
                                'ME_st': me_h-me_h,
                                'min_es': min_es_h,
                                'stance_period': stance_period_h,
                                'step_length': step_length_h,
                                'speed': speed_h,
                                'moment_peaks': mdelta_h,
                                'moment_max': mmax_h-mmax_h,
                                'spasticity_index_sol': spas_idx_sol_h,
                                'spasticity_index_gas': spas_idx_gas_h}
    
    for idx in range(len(experiment_sto_files)):
        experiment_sto = experiment_sto_files[idx]
        print(experiment_sto)
        # scores
        experiment_scores.append(float(experiment_sto.rstrip('.par.sto').split('_')[-1]))

        # ankle ME
        try:
            #me_st, me_std_st, me_sw, me_std_sw = extract_sto.me_phases(experiment_sto, healthy_sto, var_name, side)
            me_st, me_std_st, mmin_es, min_std_es = extract_sto.mean_stance(experiment_sto, var_name, side) #extract_sto.min_ankle_es(experiment_sto, var_name, side)
        except:
            me_st, me_std_st, me_sw, me_std_sw = np.nan, np.nan, np.nan, np.nan
            mmin_es, min_std_es = np.nan, np.nan
        ME_st.append(me_st)
        min_es.append(mmin_es)

        # gait features
        try:
            mstep_length, sstep_length, gspeed = extract_sto.gait_features(experiment_sto)
            mstance_period, sstance_period = extract_sto.stance_period(experiment_sto, side)
        except:
            mstep_length, sstep_length, gspeed = np.nan, np.nan, np.nan
            mstance_period, sstance_period = np.nan, np.nan

        step_length.append(mstep_length)
        stance_period.append(mstance_period)
        speed.append(gspeed)

        # muscle analysis to compute joint moments
        moment_norm = extract_sto.norm_moment(experiment_sto, var_name, side)
        var_names, var_tab = extract_sto.extract_sto(experiment_sto)
        mdelta, sdelta, mmax, smax = extract_sto.moment_peaks(var_names, var_tab, side, moment_norm)
        moment_peaks.append(mdelta)
        moment_max.append(mmax)

        # spasticity indexes
        try:
            spas_idx_sol, spas_std_sol = spasticity_index.all_cycles_spasticity_index(experiment_sto, "soleus",
                                                                                      side, plot=False)
            spas_idx_gas, spas_std_gas = spasticity_index.all_cycles_spasticity_index(experiment_sto, "gastroc",
                                                                                      side, plot=False)
        except:
            spas_idx_sol, spas_std_sol = np.nan, np.nan
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
                        'ME_st': ME_st - me_h,
                        'min_es': min_es,
                        'stance_period': stance_period,
                        'step_length': step_length,
                        'speed': speed,
                        'moment_peaks': moment_peaks,
                        'moment_max': moment_max - mmax_h,
                        'spasticity_index_sol': spasticity_idx_sol,
                        'spasticity_index_gas': spasticity_idx_gas}

    experiment_values_tuple = list(zip(*[experiment_values_stance, experiment_values_swing]))
    df_index = pd.MultiIndex.from_tuples(experiment_values_tuple, names=['stance', 'swing'])
    experiments_df = pd.DataFrame(experiments_dict, index=df_index)
    
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

        es = 'es'
        if folder_name == 'stance_KS_swing_KS':
            par = 'increased stance and swing KS'
        elif folder_name == 'stance_TA_KS_swing_TA_KS':
            par = 'increased stance and swing TA KS'
        elif folder_name == 'stance_KL_swing_KL':
            par = 'increased stance and swing KL'
            es = 'no'
        elif folder_name == 'stance_KF_stance_KL':
            par = 'decreased stance KF and KL'
            es = 'no'

        map = "RdBu"
        if metric == 'ME_st':
            lab = r'$\Delta_{h}$'+' mean angle ST [°]' #'ME ST [°]'
            if es == 'es':
                vmin = -10
                vmax = 0
            else:
                vmin = 0
                vmax = 4
        elif metric == 'min_es':
            lab = 'max angle ST [°]' #'min ES [°]'
            vmin = 15
            vmax = 20
        elif metric == 'moment_peaks':
            lab = r'$\Delta$'+' moment peaks [Nm/kg]'
            vmin = -1
            vmax = 2
        elif metric == 'moment_max':
            lab = r'$\Delta_{h}$'+' mean moment [Nm/kg]'
            vmin = -0.2
            vmax = 0
        elif metric == 'step_length':
            lab = 'step length [m]'
            vmin = 0.5
            vmax = 1
        elif metric == 'stance_period':
            lab = 'stance period [s]'
            vmin = 50
            vmax = 64
        else:
            lab = metric
            vmin = None
            vmax = None
        ax = sns.heatmap(df2, cmap=map, center=experiments_dict_healthy[metric], cbar_kws={'label': lab},
                         vmax=vmax, vmin=vmin)
        if folder_name == "stance_KF_stance_KL":
            ax.set(xlabel='stance KL [%]', ylabel='stance KF [%]')
        elif folder_name == "stance_TA_KS_swing_TA_KS":
            ax.set(xlabel='swing TA KS [%]', ylabel='stance TA KS [%]')
        else:
            ax.set(xlabel=folder_name_split[2]+" "+folder_name_split[3]+' [%]',
                   ylabel=folder_name_split[0]+" "+folder_name_split[1]+' [%]')

        if metric == 'spasticity_index_sol':
            title = 'Mean SOL spasticity index ('+side+') \n for '+par
        elif metric == 'spasticity_index_gas':
            title = 'Mean GAS spasticity index ('+side+') \n for '+par
        elif metric == 'ME_st':
            #title = 'Mean ankle angle ME over the ST phase \n for '+par
            title = r'$\Delta_{h}$' + ' mean ankle angle over the ST phase \n for ' + par
        elif metric == 'min_es':
            #title = 'Mean ankle angle min during the ES phase \n for '+par
            title = 'Max ankle angle over the ST phase \n for ' + par
        elif metric == 'moment_peaks':
            title = r'$\Delta$'+' ankle moment peaks \n for ' + par
        elif metric == 'moment_max':
            title = r'$\Delta_{h}$' +' mean ankle moment \n for ' + par
        elif metric == 'step_length':
            title = 'Mean step length \n for '+par
        else:
            title = ''
            for m in metric.split():
                title += m + ' '
            title = 'mean ' + metric + '\n for ' + par
        plt.title(title)
        plt.savefig(os.path.join(report_folder, metric+'.png'))
        plt.close()
