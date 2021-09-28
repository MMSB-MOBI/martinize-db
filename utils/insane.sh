#!/bin/bash
# GL
# Following statements are left untouched, as in a centos-module/slurm
# execution context the correct python interpreter/venv is already set-up
# Place here commands to load the virtual env that contains insane
if [ ! -z "$venv" ]
then
  source $venv
fi

insane_path="insane"

insaneHackBefore=$1
insaneHackAfter=$2
inputFile=$3
shift 3

python3 $insaneHackBefore $inputFile output-insane-hack.pdb hacked-atoms.txt

$insane_path $@

python3 $insaneHackAfter system.gro system-insane-hack.gro hacked-atoms.txt
