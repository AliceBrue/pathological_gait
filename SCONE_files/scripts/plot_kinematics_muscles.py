# Plot kinematics and muscle activity of a SCONE simulation.
#
# author: Dimitar Stanev <dimitar.stanev@epfl.ch>
# %%
import os
from utils import *

state_file = os.path.abspath('C:/Users/Acer/Documents/BioRob/HBP/NR/SPD002/SCONE/SPD002_Ong_main.sto')  #MR012/MR012_Ong/healthy_gait_ong_thelen.Sto') #'../states/optimisation.sto')
state = read_from_storage(state_file)
side = 'r'

# perform muscle analysis to compute muscle induced moments
osim_model = 'C:/Users/Acer/Documents/BioRob/HBP/NR/SPD002/SPD002_gait9dof18musc_Thelen.osim'  #MR012/generic/ong_gait9dof18musc_Thelen.osim' #'../model/gait0914.osim'
analysis_dir = 'C:/Users/Acer/Documents/BioRob/HBP/NR/SPD002/SCONE/analysis'  #MR012/MR012_Ong/healthy/'  #'../states/muscle_analysis/'
#perform_muscle_analysis(osim_model, state_file, analysis_dir)

# plot joint kinematics
# plot_scone_kinematics(state, side, state_file)
plot_scone_joint_kinematics(state, osim_model, analysis_dir, side, output_file=state_file)

# plot musle activation
muscles = ['iliopsoas', 'glut_max', 'hamstrings', 'bifemsh', 'rect_fem', 'vasti', 'tib_ant',
           'gastroc', 'soleus']
plot_scone_muscle_activations(state, muscles, side, col=4, output_file=state_file)

plt.show()
