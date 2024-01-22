#!/bin/bash

FILE="$(pwd)/sequential-knn"

if [ $# -ne 4 ] 
then
	echo -e "The script should have 4 arguments:\n"
	echo -e "\t run.sh N K BS L	\n"
	echo -e "N: Number of points to generate \nK: Number of neighbours to account for\nBS: Size of block \nL: Upper boundary for point generation \n"
else
	gcc -o sequential-knn sequential-knn.c -lm
	if [ -a $FILE ]
	then
		./sequential-knn $1 $2 $3 $4
	fi
fi


