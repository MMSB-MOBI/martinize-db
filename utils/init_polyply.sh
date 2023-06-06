if [ ! -z "$SLURM_SUBMIT_DIR" ]
then
    cd $SLURM_SUBMIT_DIR
fi

if [ ! -z "$venv" ]
then
    source $venv
fi

echo VERSION :$( polyply -V | grep '^polyply')

for ff in $(echo $forcefields | sed -n 1'p' | tr ',' '\n');
do
    echo FORCEFIELD :$ff
    polyply -list-blocks $ff | grep -v "INFO" | grep -v "The following"
done;
