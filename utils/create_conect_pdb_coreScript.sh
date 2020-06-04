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
  gmx editconf -f "$pdb" -o "$gro_box" -box 15 15 18 -noc  >editconf.out 2>editconf.err
fi

# Create the computed topology .tpr
gmx grompp -f "$mdp" -c "$gro_box" -p "$top" -o "$tpr_run" >grompp.out 2>grompp.err

if [ $DEL_WATER_BOOL == "YES" ]
then
  # File to give on stdin to make_ndx
  printf '!"W"\nq\n' > $tmp_stdin
  # Create index with a category without W
  gmx make_ndx -f "$gro_box" -o "$index_ndx" < $tmp_stdin >make_ndx.out 2>make_ndx.err

  # File to give on stdin to trjconv
  printf '!W\n' > $tmp_stdin
  # Create the PDB with conect entries without water 
  gmx trjconv -n "$index_ndx" -s "$tpr_run" -f "$gro_box" -o "$output_conect_no_water" -conect < $tmp_stdin >trjconv.out 2>trjconv.err

  echo "File $output_conect_no_water has been written."
fi

# File to give on stdin to trjconv
printf '0\n' > $tmp_stdin
# Create the PDB with conect entries with water 
gmx trjconv -s "$tpr_run" -f "$gro_box" -o "$output_conect" -conect < $tmp_stdin  >trjconv2.out 2>trjconv2.err

echo "File $output_conect has been written."
