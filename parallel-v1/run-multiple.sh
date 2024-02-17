	#!/bin/bash

FILE="$(pwd)/parallel-knn-v1"

N=(65536)

#2048 4096 8192 16384 32768

P=(2 4 8 16 32)


if [ $# -ne 1 ] 
then
	echo -e "The script should have 1 argument:\n"
	echo -e "\t run.sh K \n"
	echo -e "K: Number of neighbours to account for\n"
else
	mpicc -o parallel-knn-v1 parallel-knn-v1.c -lm
	for n in ${N[@]}
	do
		echo -e "$n ========================================= \n" >> ./results.txt
		
		for p in ${P[@]}
		do
			mpirun -np $p parallel-knn-v1 $n $1 100 >> ./results.txt
		done
		
		echo -e "\n\n" >> ./results.txt
	done
fi


