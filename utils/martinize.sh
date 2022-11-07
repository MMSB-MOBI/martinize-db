#!/bin/bash
# GL
# Following statements are left untouched, as in a centos-module/slurm
# execution context the correct path to martinize2 is already set-up
# Place here commands to load the virtual env that contains martinize2
# path


if [ ! -z "$venv" ]
then
  source $venv
fi

martinize2_path="martinize2"
input="input/input.pdb"

cmd_line="$martinize2_path -f $input -x $OUTPUT_PDB -o $OUTPUT_TOP -ff $FORCE_FIELD -p $POSITION -dssp $DSSP_PATH -cys $CYSTEIN_BRIDGE -cter $C_TER -nter $N_TER -maxwarn 99999"

if [ ! -z "$DSSP_PATH" ] 
then 
  cmd_line="$cmd_line -dssp $DSSP_PATH"
fi

if [ ! -z $SIDE_CHAIN_FIX]
then  
  cmd_line="$cmd_line -scfix"
fi

echo $cmd_line

$cmd_line 2> martinize_redirect.stderr 

{ grep "WARNING" martinize_redirect.stderr > $MARTINIZE_WARN || true; }

