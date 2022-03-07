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
import extract_sto


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


def assess_parameter_folder_1d(parameter_folder, side, scone_folder):
    '''Assesses the metrics for each of the explored parameters
    and generates a report folder containing the corresponding plots.
    Parameters
    ----------
    parameter_folder: (string) experiment folder path.
    side: (sting) 'l' or 'r' for left or right leg
    scone_folder: (string) scone results folders suffix depending on SCONE installation

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
                                       scone_folder) for experiment_value in experiment_values]
    experiment_sto_files = [get_experiment_sto_file(experiment_folder) for experiment_folder in experiment_folders]

    experiment_sto_success = []
    experiment_values_success = []
    experiment_values_float_success = []
    
    # scores and simulation times
    experiment_scores = []
    total_times = []
    
    # joint kinematics
    var_names = ['ankle_angle_' + side, 'knee_angle_' + side, 'hip_flexion_' + side]

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

        # score and simulation time
        experiment_scores.append(float(experiment_sto.rstrip('.par.sto').split('_')[-1]))
        total_times.append(psm.get_sto_total_time(experiment_sto))

        # gait features
        try:  # metric set to nan if there is not gait cycle over which to compute (ie the model falls)
            mstep_l, sstep_l, speed = extract_sto.gait_features(experiment_sto)
        except:
            mstep_l, sstep_l, speed = np.nan, np.nan, np.nan
        mstep_length.append(mstep_l)
        sstep_length.append(sstep_l)
        gspeed.append(speed)

        try:
            mstance_p, sstance_p = extract_sto.stance_period(experiment_sto, 'l')
        except:
            mstance_p, sstance_p = np.nan, np.nan
        mstance_period.append(mstance_p)
        sstance_period.append(sstance_p)

        experiment_sto_success.append(experiment_sto)
        experiment_values_success.append(experiment_values[idx])
        experiment_values_float_success.append(float(experiment_values[idx]))
    
    experiments_dict = {'parameter_value': experiment_values_success,
                        'score': experiment_scores,
                        'total_time': total_times,
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

    experiments_dict_healthy = {'score': score_healthy,
                                'total_time': total_time_healthy,
                                'mstep_length': h_mstep_length,
                                'sstep_length': h_sstep_length,
                                'speed': h_gspeed,
                                'mstance_period': h_mstance_period,
                                'sstance_period': h_sstance_period}
    
    experiment_title = get_experiment_title(parameter_folder)
    healthy_value = 100

    # export joints plots
    title = ""
    inv = False
    es = 'es'  # to compute toe gait metrics
    if experiment_title.split(" ")[0] == "biomechanical":
        inv = True
        for n in range(len(experiment_title.split(" ")[1:])):
            title = title + " " + experiment_title.split(" ")[1+n]
        if "force" in experiment_title.split(" "):
            es = "no"  # to compute heel gait metrics
    else:
        s = 1
        if "weak" in experiment_title.split(" "):
            inv = True
            es = "no"
            for n in range(len(experiment_title.split(" "))-1):
                if experiment_title.split(" ")[s+n] != "weak":
                    title = title + " " + experiment_title.split(" ")[s+n]
        else:
            if "TA" in experiment_title.split(" "):
                inv = True
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
            psm.plot_metric(experiments_dict[metric], experiment_values_float_success, metric, title,
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
            psm.plot_metrics_std(experiments_dict[metric], experiments_dict['sstance_period'],
                                 experiment_values_float_success, metric, title, report_folder, healthy_value,
                                 experiments_dict_healthy[metric], experiments_dict_healthy['sstance_period'],
                                 metric_values2=experiments_dict['mstep_length'],
                                 std_values2=experiments_dict['sstep_length'],
                                 healthy_metric2=experiments_dict_healthy['mstep_length'],
                                 std_healthy2=experiments_dict_healthy['sstep_length'], inv=inv)


def assess_parameter_folder_2d(parameter_folder, side, scone_folder):
    '''Assesses the metrics for each of the explored 2D parameter combinations
    and generates a report folder containing the corresponding plots.
    Parameters
    ----------
    parameter_folder: (string) experiment folder path.
    side: (string) 'l' or 'r' for left or right leg
    scone_folder: (string) scone results folders suffix depending on SCONE installation

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
                                       '_' + str(experiment_value.split('_')[1]) + scone_folder)
                          for experiment_value in experiment_values]
    experiment_sto_files = [get_experiment_sto_file(experiment_folder) for experiment_folder in experiment_folders]

    experiment_sto_success = []
    experiment_values_success = []
    experiment_scores_success = []
    
    experiment_values_stance = []
    experiment_values_swing = []
    
    # scores and simulation times
    experiment_scores = []
    experiment_total_time = []
    
    # joints kinematics
    var_name = 'ankle_angle_' + side
    st_me = []
    st_max = []
    moment_peaks = []
    moment_mean = []

    # gait features
    stance_period = []
    step_length = []
    speed = []

    healthy_sto = get_healthy_sto()
    # healthy metrics
    score_healthy = float(healthy_sto.rstrip('.par.sto').split('_')[-1])
    experiment_scores_success.append(score_healthy)
    total_time_healthy = psm.get_sto_total_time(healthy_sto)
    experiment_total_time.append(total_time_healthy)

    h_me, _, h_max, _ = extract_sto.mean_stance(healthy_sto, var_name, side)
    st_me.append(h_me)
    st_max.append(h_max)

    step_length_h, sstep_length, speed_h = extract_sto.gait_features(healthy_sto)
    stance_period_h, sstance_period = extract_sto.stance_period(healthy_sto, side)
    step_length.append(step_length_h)
    stance_period.append(stance_period_h)
    speed.append(speed_h)

    moment_norm = extract_sto.norm_moment(healthy_sto, var_name, side)
    var_names, var_tab = extract_sto.extract_sto(healthy_sto)
    mdelta_h, sdelta_h, mmean_h, smean_h = extract_sto.moment_metrics(var_names, var_tab, side, moment_norm)
    moment_peaks.append(mdelta_h)
    moment_mean.append(mmean_h)

    experiments_dict_healthy = {'score': score_healthy,
                                'total_time': total_time_healthy,
                                'ME_st': h_me-h_me,
                                'max_st': h_max,
                                'stance_period': stance_period_h,
                                'step_length': step_length_h,
                                'speed': speed_h,
                                'moment_peaks': mdelta_h,
                                'moment_mean': mmean_h-mmean_h
                                }

    experiment_sto_success.append(healthy_sto)
    experiment_values_success.append('100_100')
    experiment_values_stance.append(100)
    experiment_values_swing.append(100)

    for idx in range(len(experiment_sto_files)):
        experiment_sto = experiment_sto_files[idx]
        print(experiment_sto)
        
        # score and simulation time
        experiment_scores.append(float(experiment_sto.rstrip('.par.sto').split('_')[-1]))
        try:
            experiment_total_time.append(psm.get_sto_total_time(experiment_sto))
        except:
            experiment_total_time.append(np.nan)

        # ankle ME
        try:
            me, _, max, _ = extract_sto.mean_stance(experiment_sto, var_name, side)
        except:
            me, max = np.nan, np.nan
        st_me.append(me)
        st_max.append(max)

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
        try:
            mdelta, sdelta, mmean, smean = extract_sto.moment_metrics(var_names, var_tab, side, moment_norm)
        except:
            mdelta, sdelta, mmean, smean = np.nan, np.nan, np.nan, np.nan
        moment_peaks.append(mdelta)
        moment_mean.append(mmean)

        experiment_sto_success.append(experiment_sto)
        experiment_values_success.append(experiment_values[idx])
        experiment_values_stance.append(int(float(experiment_values[idx].split('_')[0])))
        experiment_values_swing.append(int(float(experiment_values[idx].split('_')[1])))
        experiment_scores_success.append(experiment_scores[idx])

    experiments_dict = {'score': experiment_scores_success,
                        'total_time': experiment_total_time,
                        'ME_st': st_me - h_me,
                        'max_st': st_max,
                        'stance_period': stance_period,
                        'step_length': step_length,
                        'speed': speed,
                        'moment_peaks': moment_peaks,
                        'moment_mean': moment_mean - mmean_h
                        }

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
        #df2.interpolate(method='akima', axis=1, inplace=True, limit_direction='both')
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
            lab = r'$\Delta_{h}$'+' mean angle ST [°]'
            if es == 'es':
                vmin = -10
                vmax = 0
            else:
                vmin = 0
                vmax = 4
        elif metric == 'max_st':
            lab = 'max angle ST [°]'
            vmin = 15
            vmax = 20
        elif metric == 'moment_peaks':
            lab = r'$\Delta$'+' moment peaks [Nm/kg]'
            vmin = -1
            vmax = 2
        elif metric == 'moment_mean':
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
        ax.patch.set(hatch='x', edgecolor='grey')
        if folder_name == "stance_KF_stance_KL":
            ax.set(xlabel='stance KL [%]', ylabel='stance KF [%]')
        elif folder_name == "stance_TA_KS_swing_TA_KS":
            ax.set(xlabel='swing TA KS [%]', ylabel='stance TA KS [%]')
        else:
            ax.set(xlabel=folder_name_split[2]+" "+folder_name_split[3]+' [%]',
                   ylabel=folder_name_split[0]+" "+folder_name_split[1]+' [%]')

        if metric == 'ME_st':
            title = r'$\Delta_{h}$' + ' mean ankle angle over the ST phase \n for ' + par
        elif metric == 'max_st':
            title = 'Max ankle angle over the ST phase \n for ' + par
        elif metric == 'moment_peaks':
            title = r'$\Delta$'+' ankle moment peaks \n for ' + par
        elif metric == 'moment_mean':
            title = r'$\Delta_{h}$' + ' mean ankle moment \n for ' + par
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


def assess_parameter_folder_rep(parameter_folder, side, scone_folder):
    '''Assesses the metrics for each of the repeated experiment (with various initial condition)
    and generates a report folder containing the corresponding plots.
    Parameters
    ----------
    parameter_folder: (string) experiment folder path.
    side: (string) 'l' or 'r' for left or right leg
    scone_folder: (string) scone results folders suffix depending on SCONE installation

    Returns
    -------
    None.
    '''
    parameter_folder_str = str(parameter_folder)
    param = parameter_folder_str.split("\\")[-1]
    param = param.split('.')[0]

    # various initial conditions (IC)
    experiment_values = [parameter_value.split('_')[-1] for parameter_value in os.listdir(parameter_folder_str)
                     if os.path.isdir(os.path.join(parameter_folder_str, parameter_value)) and
                     len(parameter_value.split('.')) > 1]
    experiment_values = [parameter_value.split('.')[0] for parameter_value in experiment_values]

    # sort IC ascendingly
    experiment_values_idx_sorted = np.argsort([float(experiment_value) for experiment_value in
                                               experiment_values])
    experiment_values = [experiment_values[idx] for idx in experiment_values_idx_sorted]

    experiment_folders = [os.path.join(parameter_folder_str, param + '_' + str(experiment_value) +
                                       scone_folder) for experiment_value in experiment_values]
    experiment_sto_files = [get_experiment_sto_file(experiment_folder) for experiment_folder in experiment_folders]

    experiment_sto_success = []
    experiment_values_success = []
    experiment_values_float_success = []

    # joints kinematics
    joint_names = ['ankle_angle_' + side, 'knee_angle_' + side, 'hip_flexion_' + side]

    # Mean absolute error (MAE) joints kinematics and optimised parameters
    mae_joints = []
    std_joints = []
    mae_param = []
    std_param = []

    # reference first initial condition
    ref_sto = experiment_sto_files[0]
    ref_par = ref_sto.split('.sto')[0]
    ref_var_names, ref_var_tab = extract_sto.extract_sto(ref_sto)
    for idx in range(len(experiment_sto_files)-1):
        experiment_sto = experiment_sto_files[1+idx]
        print(experiment_sto)

        # MAE joints kinematics
        sto_mae_joints = []
        sto_std_joints = []
        var_names, var_tab = extract_sto.extract_sto(experiment_sto)
        for joint_name in joint_names:
            try:
                _, _, av_gait_cycle, _ = extract_sto.mean_gait_phases(var_names, var_tab, joint_name, side)
                _, _, ref_av_gait_cycle, _ = extract_sto.mean_gait_phases(ref_var_names, ref_var_tab, joint_name, side)
                sto_mae_joints.append(np.mean(np.abs(av_gait_cycle-ref_av_gait_cycle)))
                sto_std_joints.append(np.std(np.abs(av_gait_cycle - ref_av_gait_cycle)))
            except:
                sto_mae_joints.append(np.nan)
                sto_mae_joints.append(np.nan)
        mae_joints.append(np.mean(sto_mae_joints))
        std_joints.append(np.mean(sto_std_joints))

        # MAE optimised parameters
        experiment_par = experiment_sto.split('.sto')[0]
        try:
            sto_mae_param, sto_std_param = extract_sto.mae_param(experiment_par, ref_par)
        except:
            sto_mae_param, sto_std_param = np.nan, np.nan
        mae_param.append(sto_mae_param)
        std_param.append(sto_std_param)

        experiment_sto_success.append(experiment_sto)
        experiment_values_success.append(experiment_values[idx])
        experiment_values_float_success.append(float(experiment_values[idx]))

    experiments_dict = {'parameter_value': experiment_values_success,
                        'mae_joints': mae_joints,
                        'std_joints': std_joints,
                        'mae_param': mae_param,
                        'std_param': std_param}

    experiments_df = pd.DataFrame.from_dict(experiments_dict, orient='index').T

    report_folder = os.path.join(parameter_folder_str, 'report')
    if not os.path.exists(report_folder):
        os.mkdir(report_folder)
    experiments_df.to_csv(os.path.join(report_folder, 'experiment_results.csv'), sep='\t')

    experiment_title = get_experiment_title(parameter_folder)

    # export joints plots
    title = ""
    es = False
    if experiment_title.split(" ")[0] == "biomechanical":
        for n in range(len(experiment_title.split(" ")[1:])):
            title = title + " " + experiment_title.split(" ")[1 + n]

    for joint_name in joint_names:
        psm.plot_mean_gc(experiment_sto_success, experiment_values_float_success, joint_name, side, title,
                         report_folder, healthy_sto=ref_sto, healthy_value=1, es=es, ic=True)

    # export MAE plots
    x_label = ""
    for n in range(len(experiment_title.split(" ")[1:])):
        x_label = x_label + " " + experiment_title.split(" ")[1 + n]
    metric = "ic"
    healthy_value = 1
    healthy_metric = 0
    healthy_std = 0
    psm.plot_metrics_std(experiments_dict['mae_joints'], experiments_dict['std_joints'], experiment_values_float_success,
                        metric, x_label, report_folder, healthy_value, healthy_metric, healthy_std,
                        metric_values2=experiments_dict['mae_param'], std_values2=experiments_dict['std_param'],
                        healthy_metric2=0, std_healthy2=0, plot=False)
