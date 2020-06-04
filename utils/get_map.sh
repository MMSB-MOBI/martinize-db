#!/bin/bash
# GL 
# Following statements are left untouched, as in a centos-module/slurm
# execution context the correct python interpreter/venv is already set-up
# Place here commands to load the virtual env that contains ccmap
# Insert HERE the path to Python binary (do NOT include path to get_map.py)

if [ -z "$venv" ]
then
  source $venv/bin/activate
fi

python_with_ccmap="python"

# ------------- #
##### USAGE #####
# ------------- #
# script_name.sh path_to_get_map.py pdb_filename distance_file

grep 'CA' "$2" > _backbones.pdb

$python_with_ccmap $1 -f _backbones.pdb -o $3
