#!/bin/bash

# GL
# Following statements are left untouched, as in a centos-module/slurm
# execution context the correct python interpreter/venv is already set-up
# Place here commands to load the virtual env that contains python3 with numpy (martinize2 venv)

cd $WORKDIR
pwd


# path

n_atoms=$(cut -f1 -d ' ' $INPUT_PDB | grep -c ATOM)

echo run python $GO_VIRT_SCRIPT -s $INPUT_PDB -f $MAP_FILE --moltype $MOLTYPE --Natoms $n_atoms
python $GO_VIRT_SCRIPT -s $INPUT_PDB -f $MAP_FILE --moltype $MOLTYPE --Natoms $n_atoms
