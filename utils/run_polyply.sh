if [ ! -z "$SLURM_SUBMIT_DIR" ]
then
    cd $SLURM_SUBMIT_DIR
fi

if [ ! -z "$polyplyenv" ]
then
    source $polyplyenv
fi

polyply -V > version

ITPOUT="polymere.itp"
GROOUT="out.gro"

if [ $action == "itp" ]

then
    polyply gen_params -f input/monfichier.itp -lib $ff -o $ITPOUT -seqf input/polymer.json -name $name > polyply.out 2> polyply.err
    
    [ -f $ITPOUT ] && cat $ITPOUT
    echo "STOP"
    [ -f polyply.err ] && cat polyply.err
fi

if [ $action == "gro" ]

then
    cat input/polymere.itp > polymere.itp
    cat input/system.top > system.top
    
    polyply gen_coords -p system.top -o $GROOUT -name $name -dens $density > polyply2.out 2> polyply2.err
    

    [ -f $GROOUT ] && cat $GROOUT
    echo "STOP"
    [ -f polyply2.err ] && cat polyply2.err

fi