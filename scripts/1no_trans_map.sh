#!/bin/sh
set -x
source env_vars.sh
for i in 0 1 2 3 4 ; 
do
	for j in 0 1 2 3 4 ; 
	do
		for tag in "1x" "10x";
		do 
			cd  ${topdirn}/tree$i/$tag/d${j}
			echo "#NEXUS" > make_map_paup_cmd.nex 
			echo "execute ${resultsdir}/correct/tree$i/modeltree.tre;" >> make_map_paup_cmd.nex 
			echo "gettrees Unrooted from = 1 to = 1 mode=3 file= ${runn}.best.tre ; " >> make_map_paup_cmd.nex 
			echo "savetre from=1 to=1 file=map.tre format=nexus  replace ;" >> make_map_paup_cmd.nex 
			echo "quit;" >> make_map_paup_cmd.nex 
			paup -n make_map_paup_cmd.nex 
			cd -
		done
	done
done
