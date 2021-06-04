#!/bin/bash
# GL
# Following statements are left untouched, as in a centos-module/slurm
# execution context the correct path to martinize2 is already set-up
# Place here commands to load the virtual env that contains martinize2
# path

echo "$(pwd)"
echo "$martinizeArgs $basedir"
cd $basedir
martinize2_path="martinize2"
$martinize2_path $martinizeArgs -maxwarn 100000000 1>martinize2.stdout 2>martinize2.stderr
