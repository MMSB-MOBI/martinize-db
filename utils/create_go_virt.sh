#!/bin/bash

# Place here commands to load the virtual env that contains python3 with numpy (martinize2 venv)


# path
python_path="python"

n_atoms=$(cut -f1 -d \' \' output.pdb | grep -c ATOM)

$python_path $@ --Natoms $n_atoms
