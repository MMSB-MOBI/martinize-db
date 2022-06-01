
if [ $action == "itp" ]

then
    if [ ! -z "$polyplyenv" ]
    then
        source $polyplyenv
    fi
    
    cp input/json.inp monjson.json
    ITPOUT="polymere.itp"
    GROOUT="out.gro"
    
    polyply gen_params -f $file -lib $ff -o $ITPOUT -seqf monjson.json -name $name > polyply.out 2> polyply.err
    
    #-f martini_v3.0.0_phospholipids_v1.itp
    
    [ -f $ITPOUT ] && cat $ITPOUT
    echo "STOP"
    [ -f polyply.err ] && cat polyply.err
fi

if [ $action == "gro" ]

then
    source $polyplyenv
    
    ITPOUT="polymere.itp"
    GROOUT="out.gro"
    
    cat $itp > polymere.itp
    echo $top > syst.top
    
    polyply gen_coords -p syst.top -o $GROOUT -name $name -dens $density > polyply2.out 2> polyply2.err
    
    [ -f $GROOUT ] && cat $GROOUT
fi