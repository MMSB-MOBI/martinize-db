#!/bin/bash

#GROMACS_LOADER="source /usr/local/gromacs/bin/GMXRC"
GROMACS_LOADER="module load gromacs/2020.5"

# #########
# # USAGE #
# #########
# Create a PDB with -conect entries using GROMACS.
#
# ./script_name.sh <PDB_OR_GRO_FILE_PATH> <TOP_FILE_PATH> <MDP_FILE_PATH> [--remove-water]
#
# In order to remove water, a group "W" must exists.
#
# Write a file output-conect.pdb in current directory.
#


function to_stderr() {
  >&2 echo $@
}

if [ $# -lt 3 ]
then
  to_stderr "You need at least pdb/gro, top and mdp file specified in arguments."
  echo "Usage: ./<script>.sh <PDB_OR_GRO_FILE_PATH> <TOP_FILE_PATH> <MDP_FILE_PATH> [--remove-water]"
  exit 1
fi

# Load the gromacs module.
$GROMACS_LOADER

pdb="$1"
top="$2"
mdp="$3"
gro_box="__box__.gro"
tpr_run="__run__.tpr"
index_ndx="__index__.ndx"
tmp_stdin="tmp"
output_conect="output-conect.pdb"
output_conect_no_water="output-conect-no-w.pdb"

# Requires: pdb in argument $1, filled top in argument $2, in the right folder
# Requires: a .mdp file in $3


if [ $4 == "--remove-water" ]
then
  # do nothing, they're already a box
  gro_box="$pdb"
else
  # Create the box
  gmx editconf -f "$pdb" -o "$gro_box" -box 15 15 18 -noc
fi

# Create the computed topology .tpr
gmx grompp -f "$mdp" -c "$gro_box" -p "$top" -o "$tpr_run"

if [ $4 == "--remove-water" ]
then
  # File to give on stdin to make_ndx
  printf '!"W"\nq\n' > $tmp_stdin
  # Create index with a category without W
  gmx make_ndx -f "$gro_box" -o "$index_ndx" < $tmp_stdin

  # File to give on stdin to trjconv
  printf '!W\n' > $tmp_stdin
  # Create the PDB with conect entries without water 
  gmx trjconv -n "$index_ndx" -s "$tpr_run" -f "$gro_box" -o "$output_conect_no_water" -conect < $tmp_stdin

  echo "File $output_conect_no_water has been written."
fi

# File to give on stdin to trjconv
printf '0\n' > $tmp_stdin
# Create the PDB with conect entries with water 
gmx trjconv -s "$tpr_run" -f "$gro_box" -o "$output_conect" -conect < $tmp_stdin

echo "File $output_conect has been written."


