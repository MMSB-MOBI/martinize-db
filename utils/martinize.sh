#!/bin/bash
# GL
# Following statements are left untouched, as in a centos-module/slurm
# execution context the correct path to martinize2 is already set-up
# Place here commands to load the virtual env that contains martinize2
# path

source /home/freaky/Documents/stage/martinize2venv/bin/activate

martinize2_path="martinize2"

$martinize2_path $@ -maxwarn 100000
