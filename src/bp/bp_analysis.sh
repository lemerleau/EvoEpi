#!/bin/bash

for u in `seq 0.05 0.05 0.40` ;
do 
	`mkdir ../../data/bp/mu/$u`; 
echo "Mutation rate:" $u
for N0 in `seq  1 1 100`;
do 
echo "N_0 = "$N0 ;
`python main.py -T 200 -sd 0.9 -mBN 5 -C 700 -N_0 $N0 -u $u -job 50`; 
done
done
