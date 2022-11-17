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

cmd_line="$martinize2_path -f $input $COMMAND_LINE -maxwarn 9999"

echo $cmd_line

$cmd_line 2> martinize_redirect.stderr 

{ grep "WARNING" martinize_redirect.stderr > $MARTINIZE_WARN || true; }

