#!/bin/bash
# GL
# Following statements are left untouched, as in a centos-module/slurm
# execution context the correct python interpreter/venv is already set-up
# Place here commands to load the virtual env that contains insane

echo "$(pwd)"
echo "$insaneArgs $basedir"
cd $basedir

# insane path
insane_path="insane"

$insane_path $insaneArgs 1>insane.log 2>insane.err
