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

#echo run : $martinize2_path $@ -maxwarn 100000

$martinize2_path $@ -maxwarn 100000 2> martinize_redirect.stderr
{ grep "WARNING" martinize_redirect.stderr > martinize_warnings.log || true; }

