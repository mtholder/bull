#!/bin/sh
set -x
source env_vars.sh
if ! test -d "${resultsdir}/${analysistag}"
then
	mkdir  "${resultsdir}/${analysistag}"
fi
for tag in "1x" "10x" ;
	do 
	for i in 0 1 2 3 4 ; 
	do
		for j in 0 1 2 3 4; 
		do
			echo $i $j $tag
			psst="${topdirn}/tree$i/$tag/d${j}/map_dists.txt"
			echo $psst
			python ../parse_dists.py -n 47 ${psst} > ${topdirn}/tree$i/$tag/d${j}/map_dist_summary.txt  
			outfn="${resultsdir}/${analysistag}/${analysistag}-summary.txt"
			/bin/echo -n "$tag $i $j " >> "$outfn" 
			python ../parse_dists.py -n 47 -c ${psst} >> "$outfn" 
			echo >> "$outfn" 
		done
	done
done
