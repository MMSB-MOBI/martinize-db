#!/bin/bash
if [ ! -z "$SLURM_SUBMIT_DIR" ]
then
    cd $SLURM_SUBMIT_DIR
fi

if [ ! -z "$venv" ]
then
    source $venv
fi

echo "BIP BIP, je suis l'ordinateur et je minimize and create pdb"

pwd 
ls  
cp input/* .

gro="input/file.gro"
top="input/file.top"
output_conect="output-conect.pdb"
output_conect_temp="outputtemp.pdb"

# Requires: pdb in argument $1, filled top in argument $2, in the right folder
# Requires: a .mdp file in $3

echo ">>$gro $top  <<"

#ajouter input em.mdp, water.gro et solvent

cp input/em.mdp .

gmx grompp -p $top -c $gro -f em.mdp -o em.tpr -po em.mdout.mdp -maxwarn 100 > 1minimization.stdout 2> 1minimization.stderr
#### IF ERROR NEED TO RAISE AN ISSUE ABOUT THE BOX LENGHT

gmx mdrun -deffnm em -v > 2runminimization.stdout 2> 2runminimization.stderr

gro_box=$gro
#Center the molecule
echo "1" "0" | gmx trjconv -f em.gro -s em.tpr -pbc mol -conect -center -o "$output_conect_temp" > 3center.stdout 2> 3center.stderr

gmx grompp -f em.mdp -c $output_conect_temp -p file.top  -o pdb.tpr -maxwarn 1 > 4make_ndx.stdout 2> 4make_ndx.stderr
gmx editconf -f pdb.tpr -conect  -o "$output_conect" > 5editconf.stdout 2> 5editconf.stderr