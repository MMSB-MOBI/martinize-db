#!/bin/bash

cd $WORKDIR
echo run contact_map $PDB redirected to output.map
contact_map $PDB > $OUTPUT

