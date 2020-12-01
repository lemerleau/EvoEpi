#!/bin/bash

for ir in `seq 0.02 0.01 0.20` ;
do 
	`mkdir ../logs/analytic_seir/Ir/$ir`; 
	`mkdir /net/stdell/people7/nonosaha/Documents/PythonCodes/FragrobEVO/test1/analytic_seir/Ir/$ir/`;
for mu in `seq  0.02 0.02 0.40`;
do 
	`mkdir ../logs/analytic_seir/Ir/$ir/$mu`

`python ppseir_basic.py -C 100 -Ir $ir -Rr 0.2 -T 500 -N 5000 -lf $ir -t 1 -mu $mu`; 

`mv ../logs/analytic_seir/Ir/$ir/$mu /net/stdell/people7/nonosaha/Documents/PythonCodes/FragrobEVO/test1/analytic_seir/Ir/$ir/`
done
done
