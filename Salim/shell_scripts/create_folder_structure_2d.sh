# This script can be used to generate proper folder structure from SCONE results folder.
# First argument is SCONE results folder. Second argument is destination folder. Make sure that it ends with '2d'.
# Make sure to modify python_script according to where the scripts will be run.
results_dir=$1
destination_dir=$2
python_script='/Users/Selim/Documents/Etudes/Automne 2020 Master Robotique/BioRob Semester Project/switchdrive/Biorob_semester_project/SCONE-Automation-Scripts/src/rename_optim_folder_2d.py'
current_dir=$(pwd)
parent_dir="$(dirname "$current_dir")"
mkdir "$destination_dir"
cp -a "$results_dir"/. "$destination_dir"

cd "$destination_dir"
mkdir complex_model
mv *_couplings_* complex_model
mkdir simple_model
mv *_solLV_* simple_model
cd simple_model
mkdir stance_KV_swing_KV
mv *stance_KV_swing_KV* stance_KV_swing_KV
mkdir stance_KL_swing_KL
mv *stance_KL_swing_KL* stance_KL_swing_KL
mkdir stance_KF_swing_KF
mv *stance_KF_swing_KF* stance_KF_swing_KF
cd ../complex_model
mkdir stance_KV_swing_KV
mv *stance_KV_swing_KV* stance_KV_swing_KV
mkdir stance_KL_swing_KL
mv *stance_KL_swing_KL* stance_KL_swing_KL
mkdir stance_KF_swing_KF
mv *stance_KF_swing_KF* stance_KF_swing_KF

cd ..

find "$destination_dir" -name 'geyer_*' -type d | while read dir_name; do python3 "$python_script" "$dir_name"; done 

cd "$current_dir"