#!/bin/bash

# GL
# Following statements are left untouched, as in a centos-module/slurm
# execution context the correct python interpreter/venv is already set-up
# Place here commands to load the virtual env that contains python3 with numpy (martinize2 venv)

cd $WORKDIR
echo $ARGS

# path

echo run python $GO_VIRT_SCRIPT $GO_ARGS
python $GO_VIRT_SCRIPT $GO_ARGS