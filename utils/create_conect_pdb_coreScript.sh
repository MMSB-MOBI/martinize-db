## TEST I/O FS status

# This a jobmanager wrappable version of create_conect_pdb.sh

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
index_cmd="index_cmd.txt"

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
  gmx editconf -f "$pdb" -o "$gro_box" -box 15 15 18 -noc  >editconf.stdout 2>editconf.stderr
fi

# Create the computed topology .tpr
echo Create tpr : gmx grompp -f "$mdp" -c "$gro_box" -p "$top" -o "$tpr_run"
gmx grompp -f "$mdp" -c "$gro_box" -p "$top" -o "$tpr_run" >grompp.stdout 2>grompp.stderr

if [ $DEL_WATER_BOOL == "YES" ]
then
  index_creation
  echo "Delete water and ions (2)" : gmx trjconv -n "$index_ndx" -s "$tpr_run" -f "$gro_box" -o "$output_conect_no_water" -conect
  # File to give on stdin to trjconv
  printf "!$sol_group_name" > $tmp_stdin
  # Create the PDB with conect entries without water 
  gmx trjconv -n "$index_ndx" -s "$tpr_run" -f "$gro_box" -o "$output_conect_no_water" -conect < $tmp_stdin >trjconv.stdout 2>trjconv.stderr

  echo "File $output_conect_no_water has been written."
fi

# File to give on stdin to trjconv
printf '0\n' > $tmp_stdin
# Create the PDB with conect entries with water 
echo Create pdb with conect entries : gmx trjconv -s "$tpr_run" -f "$gro_box" -o "$output_conect" -conect
gmx trjconv -s "$tpr_run" -f "$gro_box" -o "$output_conect" -conect < $tmp_stdin  >trjconv2.stdout 2>trjconv2.stderr

echo "File $output_conect has been written."
