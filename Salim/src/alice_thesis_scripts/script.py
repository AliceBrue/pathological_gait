"""
Script to compute and plot reflex contributions, ankle angles, spasticity indexes and hexagonal representations
"""

from spasticity_index import *
from hexagon import *
import matplotlib

# Path to results folder
results_path = 'C:/Users/Acer/Documents/Biorob/Master thesis/results/'
# Impaired evaluation of interest
model = 'cx'  # 'sp' or 'cx' for simple and complex models
param = 'swSOLKL'  # parameter of interest with 'st' or 'sw' for stance or swing phases, '', 'SOL' or 'TA' for direct
                   # pathways, coupling pathways from SOL or coupling pathways from TA, and 'KF', 'KV', 'KL' or 'C0' for
                   # force, velocity, length or constant pathways
value = '2'
impaired_sto = results_path+'evaluations/'+model+'_'+param+'_'+value+'.sto'
# Healthy optimisation to compare with
healthy_sto = results_path+model+'_healthy_optimisation/'+model+'_healthy_bestoptim.sto'
# Side of interest
side = 'r'  # 'l' or 'r' for left and right varaibles
# Plots of interest
reflex_contributions = False
spas_index = False
ankle_angle = True
hexagon = True
hexagon_file = 'sp_sw_hexagon'  # 'example', 'sp_st_hexagon', 'sp_sw_hexagon', 'cx_TA_hexagon' or 'cx_SOL_hexagon' for
                                # illustration, SOL direct gains during stance and swing, and TA to SOL and SOL to TA
                                # gains

if side == 'r':
    side_name = 'Right'
elif side == 'l':
    side_name = 'Left'

# Plot font
font = {'size': 22, 'weight': 'roman'}
matplotlib.rc('font', **font)
params = {'axes.labelsize': 22,
          'axes.labelweight': 'semibold',
          'axes.titlesize': 24,
          'axes.titleweight': 'semibold'}
plt.rcParams.update(params)

# Plot healthy reflex contributions
if reflex_contributions:
    var_list = ['soleus_'+side+'.F', 'soleus_'+side+'.V', 'soleus_'+side+'.L']
    var_legend = [side_name+' SOL F', side_name+' SOL V', side_name+' SOL L']
    ylab = 'SOL F / V / L ()'
    title = side_name+' SOL normalised F, V and L during gait'
    plot_mean_gc(healthy_sto, var_list, var_legend, side, ylab, title)
    if model == 'sp':
        var_list = ['soleus_'+side+'.RF', 'soleus_'+side+'.RV', 'soleus_'+side+'.RL']
        var_legend = [side_name+' SOL F excitation', side_name+' SOL V excitation', side_name+' SOL L excitation']
        ylab = 'SOL F / V / L excitation ()'
        title = side_name+' SOL F, V and L reflex excitation during gait'
        plot_mean_gc(healthy_sto, var_list, var_legend, side, ylab, title)
    elif model == 'cx':
        var_list = ['tib_ant_'+side+'-soleus_'+side+'.RF', 'tib_ant_'+side+'-soleus_'+side+'.RV',
                    'tib_ant_'+side+'-soleus_'+side+'.RL']
        var_legend = [side_name+' SOL F excitation', side_name+' SOL V excitation', side_name+' SOL L excitation']
        ylab = 'SOL F / V / L excitation ()'
        title = side_name+' SOL F, V and L reflex excitation during gait'
        plot_mean_gc(healthy_sto, var_list, var_legend, side, ylab, title)

# Spasticity index and related plot
if spas_index:
    mean_spas = mean_cycle_spasticity_index(impaired_sto, side)
    print(mean_spas)

# Mean and std ankle angle during all gc and stance and comparison plot
if ankle_angle:
    var_list = ['ankle_angle_'+side]
    ylab = 'Ankle angle (°)'
    title = side_name+' ankle angle comparison'
    var_legend = [side_name+' ankle angle (°)']
    mean_gc, std_gc, mean_stance, std_stance, mean_gc2, std_gc2, mean_stance2, std_stance2 = \
        plot_mean_gc(impaired_sto, var_list, var_legend, side, ylab, title, sto_file2=healthy_sto)
    print([mean_gc, std_gc, mean_stance, std_stance])
    print([mean_gc2, std_gc2, mean_stance2, std_stance2])

# Hexagonal representation
# Plot font
font = {'size': 16, 'weight': 'roman'}
matplotlib.rc('font', **font)
params = {'axes.titlesize': 16,
          'axes.titleweight': 'semibold'}
plt.rcParams.update(params)

if hexagon_file == 'example':
    param_names = ['A', 'B', 'C']
    plot_names = ['A', 'B', 'C']
    place_names = [-0.3, 0.2, 0.2, -0.3, -0.8, -0.8]
    title = 'Illustration of impairments ternary representation'
    hexagon_plot(hexagon_file, results_path, param_names, plot_names, place_names, title)
else:
    if hexagon_file == 'sp_st_hexagon':
        param_names = ['stKF', 'stKV', 'stC0']
        plot_names = ['$KF_{SOL}^{st}$', '$2*KV_{SOL}^{st}$', '$2*C0_{SOL}^{st}$']
        place_names = [-0.5, 0.2, 0.2, -0.5, -1.8, -1.8]
        param_factors = [1, 2, 2]
        axes_factor = 2
    elif hexagon_file == 'sp_sw_hexagon':
        param_names = ['swKF', 'swKV', 'swC0']
        plot_names = ['$KF_{SOL}^{sw}$', '10*' + '$KV_{SOL}^{sw}$', '10*' + '$C0_{SOL}^{sw}$']
        place_names = [-0.5, 0.05, 0.05, -0.5, -2, -2]
        param_factors = [1, 10, 10]
        axes_factor = 2
    elif hexagon_file == 'cx_SOL_hexagon':
        param_names = ['stSOLKF', 'stSOLKV', 'stSOLKL']
        plot_names = ['$KF_{SOL-TA}^{st}$', '$KV_{SOL-TA}^{st}$', '$KL_{SOL-TA}^{st}$']
        place_names = [-0.8, 0.1, 0.1, -0.8, -1.9, -1.9]
        param_factors = [1, 1, 1]
        axes_factor = 2
    elif hexagon_file == 'cx_TA_hexagon':
        param_names = ['stTAKL', 'swTAKV', 'swTAKL']
        plot_names = ['$KL_{TA-SOL}^{st}$', '$KV_{TA-SOL}^{sw}$', '$KL_{TA-SOL}^{sw}$']
        place_names = [-0.8, 0.1, 0.1, -0.8, -1.9, -1.9]
        param_factors = [1, 1, 1]
        axes_factor = 1
    title = 'Ternary representation of neural impairments \n leading to pathological gait'

    hexagon_plot(hexagon_file, results_path, param_names, plot_names, place_names, title, param_factors, axes_factor)
