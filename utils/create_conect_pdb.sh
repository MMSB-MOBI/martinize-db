## TEST I/O FS status

# This a jobmanager wrappable version of create_conect_pdb.sh

if [ ! -z "$SLURM_SUBMIT_DIR" ]
then
    cd $SLURM_SUBMIT_DIR
fi

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

if [[ -z "$DEL_WATER_BOOL" ]]
then
  >&2 echo "You need the following defined variables : <DEL_WATER_BOOL>"
  return
fi

pdb="input/$INPUT_NAME"
top="input/input.top"
mdp="input/run.mdp"
gro_box="__box__.gro"
tpr_run="__run__.tpr"
index_ndx="__index__.ndx"
tmp_stdin="tmp"
output_conect="output-conect.pdb"
output_conect_no_water="output-conect-no-w.pdb"
index_cmd="index_cmd.txt"
output_no_water="output-no-w.pdb"
output="output.pdb"

# Requires: pdb in argument $1, filled top in argument $2, in the right folder
# Requires: a .mdp file in $3

cp input/input.top input.top #Qu'est-ce qui pourrait mal se passer ?



if [ $DEL_WATER_BOOL == "YES" ]
then
  # do nothing, they're already a box
  gro_box="$pdb"
else
  # Create the box
  echo Create box : gmx editconf -f "$pdb" -o "$gro_box" -box 15 15 18 -noc
  gmx editconf -f "$pdb" -o "$gro_box" -box 15 15 18 -noc  > 1.editconf.stdout 2> 1.editconf.stderr
fi

# Checking gro box
echo "#### Input is " $INPUT_NAME > gro_check.log
head -5 "$gro_box" >> gro_check.log

# Create the computed topology .tpr
echo Create tpr : gmx grompp -f "$mdp" -c "$gro_box" -p "$top" -o "$tpr_run"
gmx grompp -f "$mdp" -c "$gro_box" -p "$top" -o "$tpr_run" > 2.grompp.stdout 2> 2.grompp.stderr
[[ -f  "$tpr_run" ]] || { >&2 echo "[create_conect] Fatal: Create tpr failed"; return;}

if [ $DEL_WATER_BOOL == "YES" ]
then
  index_creation
  printf "!$sol_group_name" > $tmp_stdin
  echo "Delete water and ions (2)" : gmx trjconv -n "$index_ndx" -s "$tpr_run" -f "$gro_box" -o "$output_conect_no_water" -conect
  gmx trjconv -n "$index_ndx" -s "$tpr_run" -f "$gro_box" -o "$output_no_water" < $tmp_stdin >3.trjconv.stdout 2>3.trjconv.stderr
  [[ -f  "$output_no_water" ]] || { >&2 echo "[trjconv:output_no_water] Fatal: Create $output_no_water failed"; return;}
  echo "File $output_no_water has been written."
  nb_atoms=$(grep -c -w ATOM $output_no_water)
  if [[ $nb_atoms -gt 99999 ]]; then
    echo "[pdb without water] Too much atoms to create connection. We use unconnected pdb"
    ln -s $output_no_water $output_conect_no_water
  else
    echo "[pdb without water] gmx trjconv -n "$index_ndx" -s "$tpr_run" -f "$gro_box" -o "$output_conect_no_water" -conect"
    gmx trjconv -n "$index_ndx" -s "$tpr_run" -f "$gro_box" -o "$output_conect_no_water" -conect < $tmp_stdin >4.trjconv-connect.stdout 2>4.trjconv-connect.stderr
    [[ -f  "$output_conect_no_water" ]] || { >&2 echo "[trjconv:output_conect_no_water] Fatal: Create $output_conect_no_water failed"; return;}
    echo "File $output_conect_no_water has been written."
  fi
fi


printf '0\n' > $tmp_stdin
# Create the PDB with conect entries with water 
echo "[Create the PDB with conect entries with water] 1/2"
echo launch : gmx trjconv -s "$tpr_run" -f "$gro_box" -o "$output"
gmx trjconv -s "$tpr_run" -f "$gro_box" -o "$output" < $tmp_stdin >5.trjconv.stdout 2> 5.trjconv.stderr
[[ -f  "$output" ]] || { >&2 echo "[trjconv:5] Fatal: Create $output failed"; return;}
    
echo "File $output has been written."
nb_atoms=$(grep -c -w ATOM $output)
if [[ $nb_atoms -gt 99999 ]]; then
  echo "[Create the PDB with conect entries with water] 2/2"
  echo "[pdb with water] Too much atoms to create connection. We use unconnected pdb"
  ln -s $output $output_conect
else
  echo "[Create the PDB with conect entries with water] 2/2"
  echo launch : gmx trjconv -s "$tpr_run" -f "$gro_box" -o "$output_conect" -conect
  gmx trjconv -s "$tpr_run" -f "$gro_box" -o "$output_conect" -conect < $tmp_stdin >6.trjconv-conect.stdout 2> 6.trjconv-conect.stderr
  [[ -f  "$output_conect" ]] || { >&2 echo "[trjconv:output_conect] Fatal: Create $output_conect failed"; return;}
  echo "File $output_conect has been written."
fi

if [[ $INPUT_TYPE == "pdb" ]]; then
  echo "Keep gro box for $pdb to final_output.gro"
  cp $gro_box final_output.gro
fi

echo "create_conect_pdb done."