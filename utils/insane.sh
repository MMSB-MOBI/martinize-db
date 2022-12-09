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

insaneArgs=$@

echo Launch : python3 $insaneHackBefore $inputFile output-insane-hack.pdb hacked-atoms.txt > cmd.txt
python3 $insaneHackBefore $inputFile output-insane-hack.pdb hacked-atoms.txt
insaneArgs2=`echo $insaneArgs | perl -pe 's/^(.*\-f\s)([\S]+)(.*)$/${1}output-insane-hack.pdb${3}/'`
echo Launch : $insane_path $insaneArgs2 >> cmd.txt
$insane_path $insaneArgs2 1> insane_redirect.stdout 2> insane_redirect.stderr
echo Launch : python3 $insaneHackAfter system.gro system-insane-hack.gro hacked-atoms.txt >> cmd.txt
python3 $insaneHackAfter system.gro system-insane-hack.gro hacked-atoms.txt
