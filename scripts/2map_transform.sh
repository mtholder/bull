#!/bin/sh
#set -x
source env_vars.sh
for i in 0 1 2 3 4  ; 
do
	for j in 0 1 2 3 4 ; 
	do
		for tag in "1x" "10x" ;
		do 
			psst="${topdirn}/tree$i/$tag/d${j}/map_transformed.tre"
			echo $i $j $tag
			cat ${modeldir}/tree$i/1x-modeltree.nex | sed '/^;$/,$d' > $psst; 
			echo ";" >> ${psst}; 
			mfn="${topdirn}/tree$i/$tag/d$j/map.tre"
			if test -f $mfn ;
			then
				python ~/Documents/projects/dendropy/trunk/examples/tree_to_splits/tree_to_splits.py -n `cat $mfn | grep PAUP_1 | awk '{print $5}'` >> $psst  || exit
			else
				echo "${mfn} not found!"
			fi
			echo "end;" >> $psst ; 
		done
	done
done
