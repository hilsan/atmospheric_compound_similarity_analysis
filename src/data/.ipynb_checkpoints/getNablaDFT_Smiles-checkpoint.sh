#!/bin/bash
file=$1
awk -F"tar," '{print $2}' $file| awk -F "," '{print $1}' | sort | uniq 
