	#!/bin/bash

FILE="$(pwd)/sequential-knn"

N=(1024 2048 4096 8192 16384 32768)

# 65536


if [ $# -ne 2 ] 
then
	echo -e "The script should have 2 arguments:\n"
	echo -e "\t run.sh K L	\n"
	echo -e "K: Number of neighbours to account for\nL: Upper boundary for point generation \n"
else
	gcc -o sequential-knn sequential-knn.c -lm
	if [ -a $FILE ]
	then
		for n in ${N[@]}
		do
			echo -e "$n ========================================= \n" >> ./results.txt
			bs=64
			
			while [ $bs -le $n ]
			do
				./sequential-knn $n $1 $bs $2 >> ./results.txt
				bs=$(( $bs * 2 ))
			done
			
			echo -e "\n\n" >> ./results.txt
		done
	fi
fi


