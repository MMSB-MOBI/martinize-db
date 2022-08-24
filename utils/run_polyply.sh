if [ ! -z "$SLURM_SUBMIT_DIR" ]
then
    cd $SLURM_SUBMIT_DIR
fi

if [ ! -z "$venv" ]
then
    source $venv
fi

ITPOUT="polymere.itp"
GROOUT="out.gro"

if [ $action == "itp" ]

then
    
    if [ -s input/monfichier.itp ]; then
        
        cat input/monfichier.itp |
        awk -v RS=";NEWITP" '{ print $0 > "file_itp" NR"custom.itp" }'
        
        for i in $(grep -l 'connexion rule' *itp)
        
        do
            mv $i $i".ff"
        done
        
        polyply gen_params -f *custom.itp* -lib $ff -o $ITPOUT -seqf input/polymer.json -name $name > polyply.out 2> polyply.err
        
    else
        polyply gen_params -lib $ff -o $ITPOUT -seqf input/polymer.json -name $name > polyply.out 2> polyply.err
    fi
    
    
    [ -f $ITPOUT ] && cat $ITPOUT
    echo "STOP"
    [ -f polyply.err ] && cat polyply.err
fi

if [ $action == "gro" ]

then
    cat input/polymere.itp > polymere.itp
    cat input/system.top > system.top
    
    polyply gen_coords -p system.top -o $GROOUT -name $name -box $box $box $box > polyply2.out 2> polyply2.err
    
    [ -f $GROOUT ] && cat $GROOUT
    echo "STOP"
    [ -f polyply2.err ] && cat polyply2.err
    
fi