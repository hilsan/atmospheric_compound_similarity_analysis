#!/bin/sh
infile=$1

awk -F"SMILES=" '/"SMILES/{print $2}' $infile | awk -F"\"" '/\"/{print $1}' | sort | uniq | sed '/\&gt\;\&gt\;/d'

