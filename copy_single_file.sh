#!/bin/bash

# bash script to copy any files not synched from SPARTAN

input="$1"

sshpass -p Feynman3142! scp mvirgato@spartan.hpc.unimelb.edu.au:/data/gpfs/projects/punim1360/WD_elec_capture/$1 ./$1

echo "File copied from Spartan"
