#!/bin/sh
infile=$1
awk -F"SMILES: " '/SMILES/{print $2}' $infile | sort | uniq
