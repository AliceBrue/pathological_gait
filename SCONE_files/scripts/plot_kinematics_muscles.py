# Plot kinematics and muscle activity of a SCONE simulation.
#
# author: Dimitar Stanev <dimitar.stanev@epfl.ch>
# %%
import os
from utils import *

state_file = os.path.abspath('../states/optimisation.sto')
state = read_from_storage(state_file)
side = 'r'

# perform muscle analysis to compute muscle induced moments
osim_model = '../model/gait0914.osim'
analysis_dir = '../states/muscle_analysis/'
perform_muscle_analysis(osim_model, state_file, analysis_dir)

# plot joint kinematics
# plot_scone_kinematics(state, side, state_file)
plot_scone_joint_kinematics(state, osim_model, analysis_dir, side, output_file=state_file)

# plot musle activation
muscles = ['iliopsoas', 'glut_max', 'hamstrings', 'bifemsh', 'rect_fem', 'vasti', 'tib_ant',
           'gastroc', 'soleus']
plot_scone_muscle_activations(state, muscles, side, col=4, output_file=state_file)

plt.show()
