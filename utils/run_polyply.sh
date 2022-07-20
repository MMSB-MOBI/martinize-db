if [ ! -z "$SLURM_SUBMIT_DIR" ]
then
    cd $SLURM_SUBMIT_DIR
fi

if [ ! -z "$polyplyenv" ]
then
    source $polyplyenv
fi

polyply -V > version

if [ $action == "itp" ]

then
    ITPOUT="polymere.itp"
    GROOUT="out.gro"

    #cp $file monfichier.itp

    #polyply gen_params -f monfichier.itp -lib $ff -o $ITPOUT -seqf monjson.json -name $name > polyply.out 2> polyply.err
    polyply gen_params -lib $ff -o $ITPOUT -seqf input/*.json -name $name > polyply.out 2> polyply.err
    
    [ -f $ITPOUT ] && cat $ITPOUT
    echo "STOP"
    [ -f polyply.err ] && cat polyply.err
fi

if [ $action == "gro" ]

then
    
    ITPOUT="polymere.itp"
    GROOUT="out.gro"
    
    cat input/itp > polymere.itp
    cat input/top > system.top
    
    polyply gen_coords -p system.top -o $GROOUT -name $name -dens $density > polyply2.out 2> polyply2.err
    
    [ -f $GROOUT ] && cat $GROOUT
fi