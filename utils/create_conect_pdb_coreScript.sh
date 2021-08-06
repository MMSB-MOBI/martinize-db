## TEST I/O FS status

# This a jobmanager wrappable version of create_conect_pdb.sh

function to_stderr() {
  >&2 echo $@
}

echo "cd $basedir"
cd $basedir

if [[ -z "$PDB_OR_GRO_FILE" || -z "$TOP_FILE" || -z "$MDP_FILE" || -z "$DEL_WATER_BOOL" ]]
then
  to_stderr "You need the following defined variables"
  echo "<PDB_OR_GRO_FILE_PATH> <TOP_FILE_PATH> <MDP_FILE_PATH> <DEL_WATER_BOOL>"
  exit 1
fi

pdb="$PDB_OR_GRO_FILE"
top="$TOP_FILE"
mdp="$MDP_FILE"
gro_box="__box__.gro"
tpr_run="__run__.tpr"
index_ndx="__index__.ndx"
tmp_stdin="tmp"
output_conect="output-conect.pdb"
output_conect_no_water="output-conect-no-w.pdb"

# Requires: pdb in argument $1, filled top in argument $2, in the right folder
# Requires: a .mdp file in $3

echo ">>$pdb $top $mdp<<"

if [ $DEL_WATER_BOOL == "YES" ]
then
  # do nothing, they're already a box
  gro_box="$pdb"
else
  # Create the box
  echo Create box : gmx editconf -f "$pdb" -o "$gro_box" -box 15 15 18 -noc
  gmx editconf -f "$pdb" -o "$gro_box" -box 15 15 18 -noc  > 1.editconf.stdout 2> 1.editconf.stderr
fi

# Create the computed topology .tpr
echo Create tpr : gmx grompp -f "$mdp" -c "$gro_box" -p "$top" -o "$tpr_run"
gmx grompp -f "$mdp" -c "$gro_box" -p "$top" -o "$tpr_run" > 2.grompp.stdout 2> 2.grompp.stderr

if [ $DEL_WATER_BOOL == "YES" ]
then
  echo "Delete water (1)" : gmx make_ndx -f "$gro_box" -o "$index_ndx"
  # File to give on stdin to make_ndx
  printf '!"W"\nq\n' > $tmp_stdin
  # Create index with a category without W
  gmx make_ndx -f "$gro_box" -o "$index_ndx" < $tmp_stdin > 3.make_ndx.stdout 2> 3.make_ndx.stderr

  echo "Delete water (2)" : gmx trjconv -n "$index_ndx" -s "$tpr_run" -f "$gro_box" -o "$output_conect_no_water" -conect
  # File to give on stdin to trjconv
  printf '!W\n' > $tmp_stdin
  # Create the PDB with conect entries without water 
  gmx trjconv -n "$index_ndx" -s "$tpr_run" -f "$gro_box" -o "$output_conect_no_water" -conect < $tmp_stdin > 4.trjconv.stdout 2> 4.trjconv.stderr

  echo "File $output_conect_no_water has been written."
fi

# File to give on stdin to trjconv
printf '0\n' > $tmp_stdin
# Create the PDB with conect entries with water 
echo Create pdb with conect entries : gmx trjconv -s "$tpr_run" -f "$gro_box" -o "$output_conect" -conect
gmx trjconv -s "$tpr_run" -f "$gro_box" -o "$output_conect" -conect < $tmp_stdin  > 5.trjconv.stdout 2> 5.trjconv.stderr

echo "File $output_conect has been written."
