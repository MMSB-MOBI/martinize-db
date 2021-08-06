#!/bin/bash
# GL
# Following statements are left untouched, as in a centos-module/slurm
# execution context the correct python interpreter/venv is already set-up
# Place here commands to load the virtual env that contains insane


# insane path

source /home/freaky/Documents/stage/insanevenv/bin/activate

insane_path="insane"

$insane_path $@
