#!/bin/bash

cd ~/piquelab/charvey/GxE/new_system/jointGenotyping/

myArray=( 1 2 4 5 6 7 8 9 11 12); 
for ii in ${myArray[@]}; do 
	find ../../QuASAR_results_DP${ii}/output -name "*readCounts.txt"; 
done | while read f; do 
	echo ${f/*output\//} | sed "s/_readCounts.txt//g" | tr "_" ":" | while read nm; do 
		less $f | grep -v "rsID" | awk -v nm=${nm} '{print nm":"$1,$2,$3}'; 
	done; 
done > MESH_QuASAR_master_readCounts_hets.txt
