# This script can be used to generate proper folder structure from SCONE results folder.
# First argument is SCONE results folder. Second argument is destination folder. Make sure that it ends with '1d'.
# Make sure to modify python_script according to where the scripts will be run.
results_dir=$1
destination_dir=$2
python_script='/Users/Selim/Documents/Etudes/Automne 2020 Master Robotique/BioRob Semester Project/switchdrive/Biorob_semester_project/SCONE-Automation-Scripts/src/rename_optim_folder_1d.py'
current_dir=$(pwd)
parent_dir="$(dirname "$current_dir")"
mkdir "$destination_dir"
cp -a "$results_dir"/. "$destination_dir"

cd "$destination_dir"
rm -rf *_strain_*
rm -rf *_stiffness_*

mkdir physiological_gains
mv *_simple_* physiological_gains
mkdir complex_model
mv *_couplings_* complex_model
mkdir simple_model
mv *_solLV_* simple_model
cd complex_model
mkdir swing
mv *_swing_* swing
mkdir stance
mv *_stance_* stance
cd swing
mkdir tib_ant-sol_KV
mv *_KV_* tib_ant-sol_KV
mkdir tib_ant-sol_KL
mv *_KL_* tib_ant-sol_KL
cd ../stance
mkdir tib_ant-sol_KV
mv *_KV_* tib_ant-sol_KV
mkdir tib_ant-sol_KL
mv *_KL_* tib_ant-sol_KL
cd ../../simple_model
mkdir swing
mv *_swing_* swing
mkdir stance
mv *_stance_* stance
cd swing
mkdir sol_KV
mv *_KV_* sol_KV
mkdir sol_KL
mv *_KL_* sol_KL
mkdir sol_KF
mv *_KF_* sol_KF
cd ../stance
mkdir sol_KV
mv *_KV_* sol_KV
mkdir sol_KL
mv *_KL_* sol_KL
mkdir sol_KF
mv *_KF_* sol_KF
cd ../../physiological_gains
mkdir optimal_fiber_length_ratio
mv *_optimal_fiber_length_* optimal_fiber_length_ratio
mkdir max_isometric_force_ratio
mv *_max_isometric_force_* max_isometric_force_ratio
mkdir tendon_slack_length_ratio
mv *_tendon_slack_length_* tendon_slack_length_ratio
mkdir kshapeactive
mv *_kshapeactive_* kshapeactive
cd ..


find "$destination_dir" -name 'geyer_*' -type d | while read dir_name; do python3 "$python_script" "$dir_name"; done 

cd "$current_dir"