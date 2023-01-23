#!/bin/bash
# GL
# Following statements are left untouched, as in a centos-module/slurm
# execution context the correct python interpreter/venv is already set-up
# Place here commands to load the virtual env that contains insane
if [ ! -z "$SLURM_SUBMIT_DIR" ]
then
    cd $SLURM_SUBMIT_DIR
fi


if [ ! -z "$venv" ]
then
  echo "we have venv"
  echo $venv
  source $venv
fi

insane_path="insane"

INPUT=input/input.pdb

if [[ ! -z $INPUT ]]; then 
    echo "input pdb so use insane hack"
    echo Launch : python3 $insaneHackBefore $INPUT output-insane-hack.pdb hacked-atoms.txt
    python3 $insaneHackBefore $INPUT output-insane-hack.pdb hacked-atoms.txt
    insaneArgs2=`echo $insaneArgs | perl -pe 's/^(.*\-f\s)([\S]+)(.*)$/${1}output-insane-hack.pdb${3}/'`
    echo Launch : $insane_path $insaneArgs2
    $insane_path $insaneArgs2 1> insane_redirect.stdout 2> insane_redirect.stderr
    echo Launch : python3 $insaneHackAfter system.gro system-insane-hack.gro hacked-atoms.txt
    python3 $insaneHackAfter system.gro system-insane-hack.gro hacked-atoms.txt
else
    echo "no input pdb"
    echo Launch : $insane_path $insaneArgs
    $insane_path $insaneArgs 1> insane_redirect.stdout 2> insane_redirect.stderr
fi
