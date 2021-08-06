#!/bin/bash
# GL 
# Following statements are left untouched, as in a centos-module/slurm
# execution context the correct python interpreter/venv is already set-up
# Place here commands to load the virtual env that contains ccmap
# Insert HERE the path to Python binary (do NOT include path to get_map.py)


cd $WORKDIR

# ------------- #
##### USAGE #####
# ------------- #
# script_name.sh path_to_get_map.py pdb_filename distance_file

grep 'CA' "$INPUT_PDB" > _backbones.pdb
pwd
echo run : python -m pcmap single _backbones.pdb --distance=10 --atomic
python -m pcmap single _backbones.pdb --distance=10 --atomic > $DISTANCES
