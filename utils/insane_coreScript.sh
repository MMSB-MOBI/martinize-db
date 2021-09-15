#!/bin/bash
# GL
# Following statements are left untouched, as in a centos-module/slurm
# execution context the correct python interpreter/venv is already set-up
# Place here commands to load the virtual env that contains insane

echo "$(pwd)"
echo "$insaneArgs $basedir"
echo $inputFile; 
cd $basedir

# insane path
insane_path="insane"

python3 $insaneHackBefore $inputFile output-insane-hack.pdb hacked-atoms.txt

insaneArgs2=`echo $insaneArgs | perl -pe 's/^(.*\-f\s)([\S]+)(.*)$/${1}output-insane-hack.pdb${3}/'`

$insane_path $insaneArgs2 1> insane.stdout 2> insane.stderr

python3 $insaneHackAfter system.gro system-insane-hack.gro hacked-atoms.txt