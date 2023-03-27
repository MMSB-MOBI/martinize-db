#!/bin/bash

if [ ! -z "$SLURM_SUBMIT_DIR" ]
then
    cd $SLURM_SUBMIT_DIR
fi


if [ ! -z "$RCSU_PATH" ]
then
  export PATH=$PATH:$RCSU_PATH #HACK
fi

PDB="input/input.pdb"
echo run contact_map $PDB redirected to $OUTPUT
contact_map $PDB > $OUTPUT

