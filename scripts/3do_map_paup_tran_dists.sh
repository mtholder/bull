#!/bin/sh
set -x
source env_vars.sh
for i in 0 1 2 3 4 ; 
do
	for j in 0 1 2 3 4 ; 
	do
		for tag in "1x" "10x";
		do 
			psst="${topdirn}/tree$i/$tag/d${j}/map_transformed.tre"
			echo $i $j $tag
			cd  "${topdirn}/tree$i/$tag/d${j}"
			echo "#NEXUS" > paup_cmd.nex 
			echo "execute ${resultsdir}/correct/tree$i/modeltree.tre;" >> paup_cmd.nex 
			echo "gettrees  mode = 7 file = ${psst} ; " >> paup_cmd.nex 
			echo "treedist file = map_dists.txt replace ;" >> paup_cmd.nex 
			echo "quit;" >> paup_cmd.nex 
			paup -n paup_cmd.nex 
			cd -
		done
	done
done
