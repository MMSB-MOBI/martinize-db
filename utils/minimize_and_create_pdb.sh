#!/bin/bash

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

echo "BIP BIP, je suis l'ordinateur et je minimize and create pdb"

cp input/* .

gro="input/file.gro"
top="input/file.top"
mdp="$MDP_FILE"
gro_box="__box__.gro"
index_ndx="__index__.ndx"
tmp_stdin="tmp"
output_conect="output-conect.pdb"
output_conect_no_water="output-conect-no-w.pdb"
index_cmd="index_cmd.txt"
output_no_water="output-no-w.pdb"
output="output.pdb"
output_conect_temp="outputtemp.pdb"

# Requires: pdb in argument $1, filled top in argument $2, in the right folder
# Requires: a .mdp file in $3

echo ">>$gro $top $mdp<<"
 
#ajouter input em.mdp, water.gro et solvent

# Add solvent to run the minimization
gmx solvate -cp $gro -cs water.gro -o init.gro -radius 0.21 > 1solvate.stdout 2> 1solvate.stderr

#besoind le modifier le top
# #include "/home/rmarin/Downloads/martini_v3.0.0_solvents_v1.itp"
# W 81

# If we reload the script, in case of error
test=$(grep "W" $top)
if [ -z "$test" ]
then
    echo "On met W"
    echo "W " $(grep -c "W" init.gro) >> $top
    
else
    echo "RIEN"
fi

cp input/em.mdp .

gmx grompp -p $top -c init.gro -f em.mdp -o em.tpr -po em.mdout.mdp -maxwarn 10 > 2minimization.stdout 2> 2minimization.stderr
#### IF ERROR NEED TO RAISE AN ISSUE ABOUT THE BOX LENGHT

gmx mdrun -deffnm em -v > 3runminimization.stdout 2> 3runminimization.stderr

gro_box="init.gro"
#Center the molecule
echo "1" "0" | gmx trjconv -f em.gro -s em.tpr -pbc mol -conect -center -o "$output_conect_temp" > 4center.stdout 2> 4center.stderr

echo -e "del 0-29\n! a W\nq\n" | gmx make_ndx -f "$output_conect_temp" -o blah.ndx > 5make_ndx.stdout 2> 5make_ndx.stderr

gmx editconf -f em.gro -n blah.ndx -o noW.gro  > TRUC1make_ndx.stdout 2> TYRUC1make_ndx.stderr

gmx grompp -f em.mdp -c noW.gro -p file.top -n blah.ndx -o pdb.tpr -maxwarn 1 > TRUCmake_ndx.stdout 2> TYRUCmake_ndx.stderr

gmx editconf -f pdb.tpr -conect -n blah.ndx -o "$output_conect" > 6editconf.stdout 2> 6editconf.stderr