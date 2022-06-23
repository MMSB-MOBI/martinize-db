#!/bin/bash
# GL 
# Following statements are left untouched, as in a centos-module/slurm
# execution context the correct python interpreter/venv is already set-up
# Place here commands to load the virtual env that contains ccmap
# Insert HERE the path to Python binary (do NOT include path to get_map.py)
source /data1/cecile/python_venv/martinize2venv/bin/activate

if [ ! -z "$venv" ]
then
  source $venv
fi

python_with_ccmap="python"

# ------------- #
##### USAGE #####
# ------------- #
# script_name.sh path_to_get_map.py pdb_filename distance_file

grep 'CA' "$2" > _backbones.pdb

echo run : python -m pcmap single _backbones.pdb --distance=10 --atomic
python -m pcmap single _backbones.pdb --distance=10 --atomic > distances.json
