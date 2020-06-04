

# This a jobmanager wrappable version of create_conect_pdb.sh

function to_stderr() {
  >&2 echo $@
}

if [[ -z "$PDB_OR_GRO_FILE_PATH" || -z "$TOP_FILE_PATH" || -z "$MDP_FILE_PATH" || -z "$DEL_WATER_BOOL" ]]
then
  to_stderr "You need the following defined variables"
  echo "<PDB_OR_GRO_FILE_PATH> <TOP_FILE_PATH> <MDP_FILE_PATH> <DEL_WATER_BOOL>"
  exit 1
fi


# Load the gromacs module.
#$GROMACS_LOADER
export GMXLIB="/data/databases/mobi/force_fields"
echo "...$GMXLIB..."
ls /data
echo "::$(whoami)"

echo "toto" > "/data/dev/mad/tmp/dummy.itp"

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
echo "cmd: cp input/PDB_OR_GRO_FILE_PATH.inp $pdb"

cp input/PDB_OR_GRO_FILE_PATH.inp $pdb
cp input/TOP_FILE_PATH.inp $top
cp input/MDP_FILE_PATH.inp $mdp

if [ $DEL_WATER_BOOL == "YES" ]
then
  # do nothing, they're already a box
  gro_box="$pdb"
else
  # Create the box
  gmx editconf -f "$pdb" -o "$gro_box" -box 15 15 18 -noc
fi

# Create the computed topology .tpr
gmx grompp -f "$mdp" -c "$gro_box" -p "$top" -o "$tpr_run"

if [ $DEL_WATER_BOOL == "YES" ]
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


