# Plot kinematics and muscle activity of a SCONE simulation.
#
# author: Dimitar Stanev <dimitar.stanev@epfl.ch>
# %%
import os
from utils import read_from_storage
from utils import plot_scone_kinematics
from utils import plot_scone_muscle_activations

# %%
# settings

state_file = os.path.abspath('../states/optimisation.sto')
state = read_from_storage(state_file)
side = 'r'

# %%
# plot joint kinematics

plot_scone_kinematics(state, side, state_file)


# %%
# perform muscle analysis to calculate muscle induced moments

muscles = ['iliopsoas', 'glut_max', 'hamstrings', 'vasti', 'tib_ant',
           'gastroc', 'soleus']
plot_scone_muscle_activations(state, muscles, side, col=4, output_file=state_file)

# %%
