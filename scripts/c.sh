#!/bin/sh
for i in 0 1 2 3 4
do
	head -n 58 ../24Apr08Datasets/tree${i}/1x-modeltree.nex  > tree${i}/modeltree.tre
	python ~/Documents/projects/dendropy/trunk/examples/tree_to_splits/tree_to_splits.py -nc `cat ../24Apr08Datasets/tree${i}/1x-modeltree.nex  | grep ModelTree | awk '{print $6}'` >> tree${i}/modeltree.tre 
	echo "; end ; " >> tree${i}/modeltree.tre 
done
