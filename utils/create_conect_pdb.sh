#!/bin/bash

#GROMACS_LOADER="source /usr/local/gromacs/bin/GMXRC"
GROMACS_LOADER="source /usr/local/gromacs/bin/GMXRC"

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

function index_creation(){
  echo "Index creation"
  gmx make_ndx -f "$gro_box" -o "$index_ndx" 1> makendx_dup.stdout 2> makendx_dup.stderr
  to_check=("W" "PW" "NA+" "CL-") #to pass through args
  to_del_cmd=""
  present_groups=""
  nbDel=0
  for g in ${to_check[@]};do
    nb=$(grep -w $g -c makendx_dup.stdout)
    if [[ $nb -gt 0 ]]; then
      present_groups+=$g" "
      if [[ $nb -gt 1 ]]; then
        toDel=$(grep -w $g -m 1 makendx_dup.stdout | sed 's/^ *//' | cut -f 1 -d " ")
        to_del_cmd+="del $(($toDel-$nbDel))\n"
        nbDel=$(($nbDel+1))
      fi
    fi
  done
  present_groups=$(echo "${present_groups%?}")
  if [[ $to_del_cmd ]]; then
    echo Groups duplicated
    printf "$to_del_cmd" > $index_cmd
  fi
  select_cmd=""
  sol_group_name=""
  for g in $present_groups; do
    select_cmd+='"'$g'"|'
    sol_group_name+=$g"_"
  done
  select_cmd=$(echo "${select_cmd%?}")
  sol_group_name=$(echo "${sol_group_name%?}")
  echo $select_cmd
  echo $sol_group_name
  printf "$select_cmd\n!\"$sol_group_name\"\nq\n" >> $index_cmd
  echo [make_ndx] run gmx make_ndx -f "$gro_box" -o "$index_ndx"
  gmx make_ndx -f "$gro_box" -o "$index_ndx" < $index_cmd > make_ndx.stdout 2>make_ndx.stderr
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
nGroup="$5"
gro_box="__box__.gro"
tpr_run="__run__.tpr"
index_ndx="__index__.ndx"
tmp_stdin="tmp"
output_conect="output-conect.pdb"
output_conect_no_water="output-conect-no-w.pdb"
index_cmd="index_cmd.txt"

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
  index_creation
  printf "!$sol_group_name" > $tmp_stdin
  # Create the PDB with conect entries without water 
  gmx trjconv -n "$index_ndx" -s "$tpr_run" -f "$gro_box" -o "$output_conect_no_water" -conect < $tmp_stdin >trjconv.stdout 2>trjconv.stderr
  echo "File $output_conect_no_water has been written."
fi

# File to give on stdin to trjconv
printf '0\n' > $tmp_stdin
# Create the PDB with conect entries with water 
gmx trjconv -s "$tpr_run" -f "$gro_box" -o "$output_conect" -conect < $tmp_stdin

echo "File $output_conect has been written."


