if [ ! -z "$SLURM_SUBMIT_DIR" ]
then
    cd $SLURM_SUBMIT_DIR
fi

if [ ! -z "$venv" ]
then
    source $venv
fi

polyply -list-lib > res


for ff in $(grep "^  [0-9]" res | cut -d'.' -f2);
do
    echo FORCEFIELD :$ff
    polyply -list-blocks $ff | grep -v "INFO" | grep -v "The following"
done;
