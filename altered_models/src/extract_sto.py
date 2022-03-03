"""
This script contains functions to extract variables from SCONE sto files and average them over gait phases.
"""
import numpy as np
import opensim
import os
from scipy.signal import find_peaks
import matplotlib.pyplot as plt


def extract_sto(sto_file):
    """Extracts variables from sto file:
    Parameters
    --------
    sto_file: (string) path to the sto file of interest

    Returns
    --------
    var_names: (list) list of sto file variable names
    var_tab: (array) previous variable values during simulation
    """
    file = open(sto_file, "r")

    file.readline()
    file.readline()
    nrows = int(file.readline().split('=')[1])
    ncol = int(file.readline().split('=')[1])

    line = file.readline()
    while len(line.split()) == 0 or line.split()[0] != 'endheader':
        line = file.readline()

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
    """Extracts left and right stance phase start and end indexes
    Parameters
    --------
    var_names: (list) list of variable names
    var_tab: (array) previous variable values during simulation

    Returns
    -------
    lstance_starts: (array) left stance phase start indexes
    lstance_ends: (array) left stance phase end indexes
    rstance_starts: (array) right stance phase start indexes
    rstance_ends: (array) right stance phase end indexes
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


def me_stance_phase_indexes(var_names, var_tab):
    """Extracts left and right late stance phase start and end indexes
    Parameters
    --------
    var_names: (list) list of variable names
    var_tab: (array) previous variable values during simulation

    Returns
    --------
    lstance_starts: (array) left stance phase start indexes
    lstance_ends: (array) left stance phase end indexes
    rstance_starts: (array) right stance phase start indexes
    rstance_ends: (array) right stance phase end indexes
    """
    lstate_index = np.where(var_names == 'leg0_l.state')[0][0]
    rstate_index = np.where(var_names == 'leg1_r.state')[0][0]
    lstate = var_tab[:, lstate_index]
    rstate = var_tab[:, rstate_index]

    lstance_starts = np.where(lstate == 1)[0]
    lstance_starts = lstance_starts[np.concatenate(
        ([0], np.where(lstance_starts[1:] - lstance_starts[:len(lstance_starts) - 1] > 1)[0] + 1))]
    lstance_ends = np.where(lstate == 3)[0]
    lstance_ends = lstance_ends[np.concatenate(
        ([0], np.where(lstance_ends[1:] - lstance_ends[:len(lstance_ends) - 1] > 1)[0])) + 1] - 1

    rstance_starts = np.where(rstate == 1)[0]
    rstance_starts = rstance_starts[np.concatenate(
        ([0], np.where(rstance_starts[1:] - rstance_starts[:len(rstance_starts) - 1] > 1)[0] + 1))]
    rstance_ends = np.where(rstate == 3)[0]
    rstance_ends = rstance_ends[np.concatenate(
        ([0], np.where(rstance_ends[1:] - rstance_ends[:len(rstance_ends) - 1] > 1)[0] + 1))]

    return lstance_starts, lstance_ends, rstance_starts, rstance_ends


def mean_gait_phases(var_names, var_tab, var_name, side, moment_norm=None):
    """Interpolates and averages a variable during stance and swing phases over gait cycles
    Parameters
    ---------
    var_names: (list) list of variable names
    var_tab: (array) previous variable valuess during simulation
    var_name: (string) name of the variable to average
    side: (string) 'r' or 'l' for right or left variable
    moment_norm: (array) normalised joint moment to average if not None

    Returns
    ----------
    av_stance_var: (array) averaged variable during stance phase over gait cycles
    av_swing_var: (array) averaged variable during swing phase over gait cycles
    av_gait_cycle: (float) averaged gait cycles duration
    av_stance_end: (float) averaged stance phases duration
    """
    lstance_starts, lstance_ends, rstance_starts, rstance_ends = stance_phase_indexes(var_names, var_tab)

    if moment_norm is None:
        var_index = np.where(var_names == var_name)[0][0]
        var = var_tab[:, var_index]
    else:
        var = moment_norm

    if side == 'l':
        stance_starts = lstance_starts
        stance_ends = lstance_ends
    else:
        stance_starts = rstance_starts
        stance_ends = rstance_ends
    if len(stance_starts):
        av_gait_cycle = stance_starts[1] - stance_starts[0] #int(np.mean(stance_starts[1:] - stance_starts[:len(stance_starts) - 1]))
        av_cycle_var = np.zeros((len(stance_starts) - 1, av_gait_cycle))

        # interpolation of each gait cycle on the averaged gait cycles duration
        for t in range(len(stance_starts) - 1):
            cycle_var = var[stance_starts[t]:stance_starts[t + 1]]
            av_cycle_var[t, :] = np.interp(np.arange(av_gait_cycle)/(av_gait_cycle-1)*100,
                                           np.arange(len(cycle_var))/(len(cycle_var)-1)*100, cycle_var)
        av_cycle_var0 = np.mean(av_cycle_var, axis=0)

        if stance_starts[0] < stance_ends[0]:   
            av_stance_end = int(np.mean(stance_ends[:min(len(stance_starts), len(stance_ends))] -
                                    stance_starts[:min(len(stance_starts), len(stance_ends))]))
        else:
            av_stance_end = int(np.mean(stance_ends[1:min(len(stance_starts)+1, len(stance_ends))] -
                                    stance_starts[:len(stan


            ce_ends)-1]))

        # average of the variable of interest during stance and swing phases over gait cycles
        av_stance_var = av_cycle_var0[0:av_stance_end]
        av_swing_var = av_cycle_var0[av_stance_end:av_gait_cycle]

        return av_stance_var, av_swing_var, av_gait_cycle, av_stance_end
    
    else:
        return None, None, None, None


def norm_moment(experiment_sto, var_name, side):

    sto_dir = os.path.dirname(experiment_sto)

    # perform muscle analysis
    if 'optimisation' in sto_dir:
        osim_model = '../models/scone_base_files/gait0914.osim'
    else:
        osim_model = sto_dir + '/' + 'gait0914.osim'

    if 'biomechanical' in sto_dir:
        osim_model = modify_model(osim_model, sto_dir)
    perform_muscle_analysis(osim_model, experiment_sto, sto_dir + '/muscle_analysis/')

    # find moment file of interest
    _, _, filenames = next(os.walk(sto_dir + '/muscle_analysis/'))
    moment_file = ''
    joint_name = var_name.split('_')[0]
    for f in filenames:
        if '_Moment_' + joint_name + '_angle_' + side in f:
            moment_file = f
    moment_names, moment_tab = extract_sto(sto_dir + '/muscle_analysis/' + moment_file)

    # get model total mass (for normalisation)
    model = opensim.Model(osim_model)
    s = model.initSystem()
    mass = model.getTotalMass(s)

    # normalise joint moment
    moment_norm = -np.sum(moment_tab[:, 1:], axis=1) / mass

    return moment_norm


def moment_peaks(var_names, var_tab, side, moment_norm):
    """Interpolates and averages a variable during stance and swing phases over gait cycles
    Parameters
    ---------
    var_names: (list) list of variable names
    var_tab: (array) previous variable valuess during simulation
    var_name: (string) name of the variable to average
    side: (string) 'r' or 'l' for right or left variable
    moment_norm: (array) normalised joint moment to average if not None

    Returns
    ----------
    av_stance_var: (array) averaged variable during stance phase over gait cycles
    av_swing_var: (array) averaged variable during swing phase over gait cycles
    av_gait_cycle: (float) averaged gait cycles duration
    av_stance_end: (float) averaged stance phases duration
    """
    delta_mom_peaks = []
    max_mom_peaks = []

    lstance_starts, lstance_ends, rstance_starts, rstance_ends = stance_phase_indexes(var_names, var_tab)
    var = moment_norm

    if side == 'l':
        stance_starts = lstance_starts
        stance_ends = lstance_ends
    else:
        stance_starts = rstance_starts
        stance_ends = rstance_ends

    if len(stance_starts):
        for t in range(len(stance_starts) - 1):
            if stance_starts[t] < stance_ends[t]:
                stance_mom = var[stance_starts[t]:stance_ends[t]]
            else:
                if t < len(stance_ends) - 1:
                    stance_mom = var[stance_starts[t]:stance_ends[t + 1]]

            max_mom_peaks.append(np.mean(stance_mom))
            # compute moment peaks
            first_mom_peaks, _ = find_peaks(stance_mom[:int(0.33 * len(stance_mom))],
                                            distance=int(0.33 * len(stance_mom)))
            sec_mom_peaks, _ = find_peaks(stance_mom[int(0.33 * len(stance_mom)):],
                                          distance=int(0.66 * len(stance_mom)))
            if len(sec_mom_peaks) == 0:
                if len(first_mom_peaks) > 0:
                    delta_mom_peaks.append(np.mean(stance_mom[int(0.33 * len(stance_mom)):]) -
                                           stance_mom[first_mom_peaks])
                else:
                    delta_mom_peaks.append(np.mean(stance_mom[int(0.33 * len(stance_mom)):]) -
                                           np.mean(stance_mom[:int(0.33 * len(stance_mom))]))
            else:
                if len(first_mom_peaks) > 0:
                    delta_mom_peaks.append(stance_mom[int(0.33 * len(stance_mom)) + sec_mom_peaks] -
                                           stance_mom[first_mom_peaks])
                else:
                    delta_mom_peaks.append(stance_mom[int(0.33 * len(stance_mom)) + sec_mom_peaks] -
                                           np.mean(stance_mom[:int(0.33 * len(stance_mom))]))

        return np.mean(delta_mom_peaks), np.std(delta_mom_peaks), np.mean(max_mom_peaks), np.std(max_mom_peaks)

    else:
        return None, None, None, None


def me_mean_gait_phases2(var_names, var_tab, var_name, side):
    """Interpolates and averages a variable during late stance and swing phases over gait cycles
        Parameters
    ---------
    var_names: (list) list of variable names
    var_tab: (array) previous variable values during simulation
    var_name: (string) name of the variable to average
    side: (string) 'r' or 'l' for right or left variable

    Returns
    ----------
    av_stance_var: (array) averaged variable during late stance phase over gait cycles
    av_swing_var: (array) averaged variable during swing phase over gait cycles
    av_gait_cycle: (float) averaged gait cycles duration
    av_stance_end: (float) averaged stance phases end
    av_stance_start: (float) averaged stance phases start
    """
    lstance_starts, lstance_ends, rstance_starts, rstance_ends = me_stance_phase_indexes(var_names, var_tab)
    l_starts, l_ends, r_starts, r_ends = stance_phase_indexes(var_names, var_tab)

    var_index = np.where(var_names == var_name)[0][0]
    var = var_tab[:, var_index]

    if side == 'l':
        stance_starts = lstance_starts
        stance_ends = lstance_ends
    else:
        stance_starts = rstance_starts
        stance_ends = rstance_ends
    if len(stance_starts):
        av_gait_cycle = int(np.mean(stance_starts[1:] - stance_starts[:len(stance_starts) - 1]))
        av_cycle_var = np.zeros((len(stance_starts) - 1, av_gait_cycle))

        # interpolation of each gait cycle on the averaged gait cycles duration
        for t in range(len(stance_starts) - 1):
            cycle_var = var[stance_starts[t]:stance_starts[t + 1]]
            av_cycle_var[t, :] = np.interp(np.arange(av_gait_cycle), np.arange(len(cycle_var)), cycle_var)
        av_cycle_var0 = np.mean(av_cycle_var, axis=0)

        if stance_starts[0] < stance_ends[0]:   
            av_stance_end = int(np.mean(stance_ends[:min(len(stance_starts), len(stance_ends))] -
                                    l_starts[:min(len(l_starts), len(l_ends))]))
            av_stance_start = int(np.mean(stance_starts[:min(len(stance_starts), len(stance_ends))] -
                            l_starts[:min(len(l_starts), len(l_ends))]))
        else:
            av_stance_end = int(np.mean(stance_ends[1:len(stance_ends)] -
                                    l_starts[:len(l_ends) - 1]))
            av_stance_start = int(np.mean(stance_starts[:len(stance_ends) - 1] -
                                    l_starts[:len(l_ends) - 1]))

        # average of the variable of interest during stance and swing phases over gait cycles
        av_stance_var = av_cycle_var0[0:av_stance_end-av_stance_start]
        av_swing_var = av_cycle_var0[av_stance_end-av_stance_start:av_gait_cycle-av_stance_start]

        return av_stance_var, av_swing_var, av_gait_cycle, av_stance_end, av_stance_start
    
    else:
        return None, None, None, None, None


def me_mean_gait_phases(var_names, var_tab, var_name, side):
    """Interpolates and averages a variable during late stance and swing phases over gait cycles
        Parameters
    ---------
    var_names: (list) list of variable names
    var_tab: (array) previous variable values during simulation
    var_name: (string) name of the variable to average
    side: (string) 'r' or 'l' for right or left variable

    Returns
    ----------
    av_stance_var: (array) averaged variable during late stance phase over gait cycles
    av_swing_var: (array) averaged variable during swing phase over gait cycles
    av_gait_cycle: (float) averaged gait cycles duration
    av_stance_end: (float) averaged stance phases end
    av_stance_start: (float) averaged stance phases start
    """
    lstance_starts, lstance_ends, rstance_starts, rstance_ends = me_stance_phase_indexes(var_names, var_tab)
    l_starts, l_ends, r_starts, r_ends = stance_phase_indexes(var_names, var_tab)

    var_index = np.where(var_names == var_name)[0][0]
    var = var_tab[:, var_index]

    if side == 'l':
        stance_starts = lstance_starts
        stance_ends = lstance_ends
    else:
        stance_starts = rstance_starts
        stance_ends = rstance_ends
    if len(stance_starts):
        av_gait_cycle = int(stance_starts[1] - stance_starts[0])
        av_cycle_var = np.zeros((len(stance_starts) - 1, av_gait_cycle))

        # interpolation of each gait cycle on the averaged gait cycles duration
        for t in range(len(stance_starts) - 1):
            cycle_var = var[stance_starts[t]:stance_starts[t + 1]]
            av_cycle_var[t, :] = np.interp(np.arange(av_gait_cycle)/(av_gait_cycle-1)*100,
                                           np.arange(len(cycle_var))/(len(cycle_var)-1)*100, cycle_var)
        av_cycle_var0 = np.mean(av_cycle_var, axis=0)

        if stance_starts[0] < stance_ends[0]:
            av_stance_end = int(np.mean(stance_ends[:min(len(stance_starts), len(stance_ends))] -
                                        l_starts[:min(len(l_starts), len(l_ends))]))
            av_stance_start = int(np.mean(stance_starts[:min(len(stance_starts), len(stance_ends))] -
                                          l_starts[:min(len(l_starts), len(l_ends))]))
        else:
            av_stance_end = int(np.mean(stance_ends[1:len(stance_ends)] -
                                        l_starts[:len(l_ends) - 1]))
            av_stance_start = int(np.mean(stance_starts[:len(stance_ends) - 1] -
                                          l_starts[:len(l_ends) - 1]))

        # average of the variable of interest during stance and swing phases over gait cycles
        av_stance_var = av_cycle_var0[0:av_stance_end - av_stance_start]
        av_swing_var = av_cycle_var0[av_stance_end - av_stance_start:av_gait_cycle - av_stance_start]

        return av_stance_var, av_swing_var, av_gait_cycle, av_stance_end, av_stance_start

    else:
        return None, None, None, None, None


def get_mean_gc(sto_file, var_list, side, sto_file2=None):
    """Gets averaged variables during the gait cycle, stance and swing phases for one or two simulations together.
    Parameters
    ---------
    sto_file: (string) path to the sto file of interest
    var_list: (list) list of variables to plot
    side: (string) 'r' or 'l' for right or left variables
    sto_file2: (string) optional, path to second sto file of interest. Default is None

    Returns
    ----------
    mean_gc_var: (float) mean of the variable over gait cycle
    std_gc_var: (float) std of the variable over gait cycle
    mean_stance_var: (float) mean of the variable over stance phase
    std_stance_var: (float) std of the variable over stance phase
    mean_swing_var: (float) mean of the variable over stance phase
    std_swing_var: (float) std of the variable over stance phase
    """
    var_names, var_tab = extract_sto(sto_file)

    if sto_file2 is not None:
        var_names2, var_tab2 = extract_sto(sto_file2)
        
    c = 0

    for var in range(len(var_list)):
        av_stance_var, av_swing_var, av_gait_cycle, av_stance_end = mean_gait_phases(var_names, var_tab, var_list[var],
                                                                                     side)
        if sto_file2 is not None:
            av_stance_var2, av_swing_var2, av_gait_cycle2, av_stance_end2 = mean_gait_phases(var_names2, var_tab2,
                                                                                                 var_list[var], side)
        c += 1

    if av_stance_var is not None:
        av_gc_var = np.concatenate((av_stance_var, av_swing_var))
        mean_gc_var = np.mean(av_gc_var)
        std_gc_var = np.std(av_gc_var)
        mean_stance_var = np.mean(av_stance_var)
        std_stance_var = np.std(av_stance_var)
        mean_swing_var = np.mean(av_swing_var)
        std_swing_var = np.std(av_swing_var)
    else:
        mean_gc_var = 0
        std_gc_var = 0
        mean_stance_var = 0
        std_stance_var = 0
        mean_swing_var = 0
        std_swing_var = 0

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


def me_stance(sto_file, sto_file_healthy, var_name, side):
    """Computes the mean ME between altered and healthy ankle angle during late stance period
    Parameters
    --------
    sto_file: (string) path to the sto file of interest
    sto_file_healthy: (string) path to the healthy sto file
    var_name: (string) name of the varaible of interest
    side: (string) 'l' or 'r' for left or right variable

    Returns
    --------
    mean_st_me: (float) mean ME during stance phase
    std_st_me: (float) ME std during stance phase
    mean_sw_me: (float) mean ME during swing phase
    std_sw_me: (float) ME std during swing phase
    """
    var_names, var_tab = extract_sto(sto_file_healthy)
    h_stance_var, h_swing_var, h_gait_cycle, h_stance_end = mean_gait_phases(var_names, var_tab, var_name, side)
    st_h_stance_var, st_h_swing_var, st_h_gait_cycle, st_h_stance_end, st_h_stance_start = \
        me_mean_gait_phases(var_names, var_tab, var_name, side)
    var_names, var_tab = extract_sto(sto_file)
        
    lstance_starts, lstance_ends, rstance_starts, rstance_ends = stance_phase_indexes(var_names, var_tab)
    st_lstance_starts, st_lstance_ends, st_rstance_starts, st_rstance_ends = me_stance_phase_indexes(var_names, var_tab)

    var_index = np.where(var_names == var_name)[0][0]
    var = var_tab[:, var_index]

    st_me = []

    if side == 'l':
        stance_starts = st_lstance_starts
        stance_ends = st_lstance_ends
    else:
        stance_starts = st_rstance_starts
        stance_ends = st_rstance_ends

    if len(stance_starts):
        av_stance_var = np.zeros((len(stance_starts) - 1, st_h_stance_end - st_h_stance_start))
        av_swing_var = np.zeros((len(stance_starts) - 1, st_h_stance_start))

        ## interpolation of each gait cycle on the h gait cycle duration
        for t in range(len(stance_starts) - 1):
            if stance_starts[t] < stance_ends[t]:
                cycle_var = var[stance_starts[t]:stance_ends[t]]
                """av_stance_var[t, :] = np.interp(np.arange(st_h_stance_end - st_h_stance_start), np.arange(len(cycle_var)),
                                      cycle_var)
                st_me.append(np.mean((av_stance_var[t, :] - st_h_stance_var)))"""
                st_me.append(np.mean(cycle_var) - np.mean(st_h_stance_var))
            else:
                if t < len(stance_ends) - 1:
                    cycle_var = var[stance_starts[t]:stance_ends[t + 1]]
                    """av_stance_var[t, :] = np.interp(np.arange(st_h_stance_end - st_h_stance_start), np.arange(len(cycle_var)),
                                          cycle_var)
                    st_me.append(np.mean((av_stance_var[t, :] - st_h_stance_var)))"""
                    st_me.append(np.mean(cycle_var) - np.mean(st_h_stance_var))

        return np.mean(st_me)*180/np.pi, np.std(st_me)*180/np.pi

    else:
        return None, None, None, None


def mean_stance(sto_file, var_name, side):
    """Computes the mean ankle angle during late stance period
    Parameters
    --------
    sto_file: (string) path to the sto file of interest
    var_name: (string) name of the varaible of interest
    side: (string) 'l' or 'r' for left or right variable

    Returns
    --------
    mean_st_me: (float) mean ME during stance phase
    std_st_me: (float) ME std during stance phase
    """
    var_names, var_tab = extract_sto(sto_file)
    st_lstance_starts, st_lstance_ends, st_rstance_starts, st_rstance_ends = me_stance_phase_indexes(var_names, var_tab)

    var_index = np.where(var_names == var_name)[0][0]
    var = var_tab[:, var_index]

    st_mean = []
    st_max = []

    if side == 'l':
        stance_starts = st_lstance_starts
        stance_ends = st_lstance_ends
    else:
        stance_starts = st_rstance_starts
        stance_ends = st_rstance_ends

    if len(stance_starts):
        for t in range(len(stance_starts) - 1):
            if stance_starts[t] < stance_ends[t]:
                cycle_var = var[stance_starts[t]:stance_ends[t]]
                st_mean.append(np.mean(cycle_var))
                st_max.append(np.max(cycle_var))
            else:
                if t < len(stance_ends) - 1:
                    cycle_var = var[stance_starts[t]:stance_ends[t + 1]]
                    st_mean.append(np.mean(cycle_var))
                    st_max.append(np.max(cycle_var))
        return np.mean(st_mean) * 180 / np.pi, np.std(st_mean) * 180 / np.pi, np.mean(st_max) * 180 / np.pi, \
               np.std(st_max) * 180 / np.pi

    else:
        return None, None, None, None


def min_ankle_es(sto_file, var_name, side):
    """Computes the mean min angle during ES period.
    Parameters
    --------
    sto_file: (string) path to the sto file of interest
    var_name: (string) name of the variable of interest
    side: (string) 'l' or 'r' for left or right variable

    Returns
    --------
    mean_min_es: (float) mean min angle during ES
    std_min_es: (float) min angle during ES std
    """
    var_names, var_tab = extract_sto(sto_file)

    lstance_starts, lstance_ends, rstance_starts, rstance_ends = stance_phase_indexes(var_names, var_tab)
    st_lstance_starts, st_lstance_ends, st_rstance_starts, st_rstance_ends = me_stance_phase_indexes(var_names, var_tab)

    var_index = np.where(var_names == var_name)[0][0]
    var = var_tab[:, var_index]

    min_ankle = []

    if side == 'l':
        stance_starts = st_lstance_starts
        stance_ends = st_lstance_ends
    else:
        stance_starts = st_rstance_starts
        stance_ends = st_rstance_ends

    if len(stance_starts):
        for t in range(len(stance_starts) - 1):
            if stance_starts[t] < stance_ends[t]:
                cycle_var = var[lstance_starts[t]:stance_starts[t]]
                min_ankle.append(np.min(cycle_var))
            else:
                if t < len(stance_ends) - 1:
                    cycle_var = var[lstance_starts[t]:stance_starts[t]]
                    min_ankle.append(np.min(cycle_var))

        return np.mean(min_ankle) * 180 / np.pi, np.std(min_ankle) * 180 / np.pi

    else:
        return None, None


def gait_features(sto_file):
    """Computes the mean step length and gait speed.
    Parameters
    --------
    sto_file: (string) path to the sto file of interest

    Returns
    --------
    mean_step_l: (float) mean step length
    std_step_l: (float) step length std
    speed: (float) mean gait speed
    """
    var_names, var_tab = extract_sto(sto_file)
    lcalc_x = np.where(var_names == 'calcn_l.com_pos_x')[0][0]
    rcalc_x = np.where(var_names == 'calcn_r.com_pos_x')[0][0]
    time = np.where(var_names == 'time')[0][0]
    lstance_starts, lstance_ends, rstance_starts, rstance_ends = stance_phase_indexes(var_names, var_tab)

    lcalc_x = var_tab[:, lcalc_x]
    rcalc_x = var_tab[:, rcalc_x]
    time = var_tab[:, time][-1]

    step_l = []
    if len(lstance_starts):
        if lstance_ends[0] > rstance_starts[0]:
            for s in range(len(lstance_ends)):
                step_l.append(lcalc_x[lstance_ends[s]] - rcalc_x[lstance_ends[s]])
        else:
            for s in range(len(lstance_ends) - 1):
                step_l.append(lcalc_x[lstance_ends[s + 1]] - rcalc_x[lstance_ends[s]])

        speed = (lcalc_x[-1] - lcalc_x[0]) / time

        return np.mean(step_l[1:]), np.std(step_l[1:]), speed

    else:
        return None, None, None


def stance_period(sto_file, side):
    """Computes the mean stance phase period.
    Parameters
    --------
    sto_file: (string) path to the sto file of interest
    side: (string) 'l' or 'r' for left or right variable

    Returns
    --------
    mean_stance: (float) mean stance phase
    std_stance: (float) stance period std
    """
    var_names, var_tab = extract_sto(sto_file)
    lstance_starts, lstance_ends, rstance_starts, rstance_ends = stance_phase_indexes(var_names, var_tab)

    t_index = np.where(var_names == "time")[0][0]
    time = var_tab[:, t_index]

    av_stance_end = []
    if side == 'l':
        stance_starts = lstance_starts
        stance_ends = lstance_ends
    else:
        stance_starts = rstance_starts
        stance_ends = rstance_ends

    if len(stance_starts):
        if stance_starts[0] < stance_ends[0]:
            for s in range(min(len(stance_ends), len(stance_starts) - 1)):
                av_stance_end.append((time[stance_ends[s]] -
                                      time[stance_starts[s]]) * 100 / (time[stance_starts[s + 1]] -
                                                                       time[stance_starts[s]]))
        else:
            for s in range(min(len(stance_starts), len(stance_ends) - 1)):
                av_stance_end.append((time[stance_ends[s + 1]] -
                                      time[stance_starts[s]]) * 100 / (time[stance_ends[s + 1]] - time[stance_ends[s]]))

        return np.mean(av_stance_end[1:]), np.std(av_stance_end[1:])

    else:
        return None, None


def perform_muscle_analysis(model_file, state_file, output_dir):
    """Perform OpenSim MuscleAnalysis on SCONE state file generated
    through simulation. This might be used to calculate joint moments
    induced by muscles.
    Parameters
    --------
    model_file: (string) path to osim model file
    state_file: (string) path to the sto file of interest
    output_dir: (string) path to sto file directory

    Returns
    --------
    muscle anaysis files
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


def modify_model(model_file, sto_dir):
    """
    Modify Millard muscles parameters
    Parameters
    --------
    model_file: (string) path to osim model file
    sto_dir: (string) path to sto file directory

    Returns
    --------
    osim_file: (string) path to modified osim model file
    """
    modified_muscles = ['soleus_r', 'soleus_l', 'gastroc_r', 'gastroc_l']
    file = open(model_file, 'r')
    lines = file.readlines()
    new_lines = lines
    muscle_name = ''
    for l in range(len(lines)):
        line = lines[l]
        if len(line.split()) > 0 and line.split()[0] == '<Millard2012EquilibriumMuscle':
            muscle_name = line.split('"')[-2]
        if muscle_name in modified_muscles and sto_dir.split()[0].split('\\')[-1].split('_')[0] == 'max' and \
                len(line.split()) > 0 and len(line.split()[0].split('>')) > 0 and \
                line.split()[0].split('>')[0] == '<max_isometric_force':
            new_value = str(float(sto_dir.split()[0].split('\\')[-1].split('_')[3].split('.')[0])/100 *
                        float(line.split()[0].split('>')[1].split('<')[0]))
            new_lines[l] = '\t\t\t\t\t<max_isometric_force>' + new_value + '</max_isometric_force>\n'
        if muscle_name in modified_muscles and sto_dir.split()[0].split('\\')[-1].split('_')[0] == 'opt' and \
                len(line.split()) > 0 and len(line.split()[0].split('>')) > 0 and\
                line.split()[0].split('>')[0] == '<optimal_fiber_length':
            new_value = str(float(sto_dir.split('\\')[-1].split('_')[3].split('.')[0])/100 *
                            float(line.split()[0].split('>')[1].split('<')[0]))
            new_lines[l] = '\t\t\t\t\t<optimal_fiber_length>' + new_value + '</optimal_fiber_length>\n'

    # write modified osim model file
    osim_file = sto_dir + 'modified_gait0914.osim'
    with open(osim_file, 'w') as file:
        file.writelines(new_lines)
    return osim_file
