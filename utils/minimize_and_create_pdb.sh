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
 
#cp input/* .

gro="input/file.gro"
top="input/file.top"
em_mdp="input/em.mdp"

output_conect="output-conect.pdb"
output_conect_temp="outputtemp.pdb"

# Requires: pdb in argument $1, filled top in argument $2, in the right folder
# Requires: a .mdp file in $3

echo ">>$gro $top  <<"

#ajouter input em.mdp, water.gro et solvent

#mv input/em.mdp input/martini_force_field.itp .
echo "[Running] gmx grompp -p $top -c $gro -f $em_mdp -o em.tpr -po em.mdout.mdp -maxwarn 100"
gmx grompp -p $top -c $gro -f $em_mdp -o em.tpr -po em.mdout.mdp -maxwarn 100 > 1.minimization.stdout 2> 1.minimization.stderr

#### IF ERROR NEED TO RAISE AN ISSUE ABOUT THE BOX LENGHT
echo "[Running] gmx mdrun -deffnm em -v"
gmx mdrun -deffnm em -v > 2.runminimization.stdout 2> 2.runminimization.stderr

gro_box=$gro

#Center the molecule
echo "[Running] gmx trjconv -f em.gro -s em.tpr -pbc mol -conect -center -o $output_conect_temp"
echo "1" "0" | gmx trjconv -f em.gro -s em.tpr -pbc mol -conect -center -o "$output_conect_temp" > 3.center.stdout 2> 3.center.stderr

echo "[Running] gmx grompp -f $em_mdp -c $output_conect_temp -p file.top  -o pdb.tpr -maxwarn 1"
gmx grompp -f $em_mdp -c $output_conect_temp -p $top  -o pdb.tpr -maxwarn 1 > 4.make_ndx.stdout 2> 4.make_ndx.stderr

echo "[Running] gmx editconf -f pdb.tpr -conect  -o $output_conect"
gmx editconf -f pdb.tpr -conect  -o "$output_conect" > 5editconf.stdout 2> 5editconf.stderr