# This script can be used to evaluate with SCONE the parameter files with the best scores found under the simulation results folders under results_dir.
# results_dir is a folder path entered as first argument. It needs to contain folders which are obtained from SCONE optimizations.
# Make sure to update scone_install accordingly.
results_dir=$1
current_dir=$(pwd)
scone_install=/data/benghorb/SCONE
cd $scone_install
source ./tools/linux_config
cd $current_dir
python_script=/home/benghorb/SCONE/variable_controller_generation/batch_evaluate_par.py
python3 $python_script $1 | while read file_name; do "$scone_install"/build/bin/sconecmd -e "$file_name";  done
