#!/bin/bash

if [ ! -z "$RCSU_PATH" ]
then
  export PATH=$PATH:$RCSU_PATH #HACK
fi

PDB="input/input.pdb"
echo run contact_map $PDB redirected to $OUTPUT
contact_map $PDB > $OUTPUT

