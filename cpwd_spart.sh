#!/bin/bash

#script to synch data with spartan server

echo "Copy TO or FROM Spartan?"

while : ; do
	read -p "Enter TO or FROM: " direction
	[[ $direction != 'TO' && $direction != 'FROM' ]] || break
done 

if [ $direction == 'TO' ]
then
	sshpass -p Feynman3142! rsync -avz --delete . mvirgato@spartan.hpc.unimelb.edu.au:/data/gpfs/projects/punim1360/WD_elec_capture/
elif [ $direction == 'FROM' ]
then
	sshpass -p Feynman3142! rsync -avz --update mvirgato@spartan.hpc.unimelb.edu.au:/data/gpfs/projects/punim1360/WD_elec_capture/* .
fi
	
