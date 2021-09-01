#!/bin/bash
# GL
# Following statements are left untouched, as in a centos-module/slurm
# execution context the correct path to martinize2 is already set-up
# Place here commands to load the virtual env that contains martinize2
# path

cd $basedir
martinize2_path="martinize2"
$martinize2_path $martinizeArgs -maxwarn 9999

# in case of missing all atoms coordinates in PDB files beads with nan coordinates are dumped
# which makes subsequent insane calls crash. (eq: 1PT4)
if [[ -e output.pdb ]]
then
        mv output.pdb output.pdb.may_nan
        grep -v 'nan     nan     nan' output.pdb.may_nan > output.pdb
fi
