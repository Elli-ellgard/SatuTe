#!/usr/bin/env bash

set -e

#-----------------------------------------------------------------
#get options and input handling

input_dir=$1
if [ "$input_dir" == "" ]
then
        echo "You did not specify the input directory!"
        echo "exiting ..."
        exit
fi

#------------------------------------------------------
# input paths and files
fscripts=$(pwd)

#------------------------------------------------------
#Software and Scripts
annotateTree=$fscripts/script-annotate-tree-satute.R
iqtree2=/home/elgert/IQ-TREE/iqtree-2.2.2.4-Linux/bin/iqtree2
#------------------------------------------------------

number_rates=0
for res in $input_dir/resultsRate*
do
	number_rates=$(expr $number_rates + 1)
done

tree_file=$(basename $input_dir/*.treefile)
echo "---------------------------"
echo "Consider $treefile"
echo "---------------------------"
echo ""

for res_file in $input_dir/resultsRate*
do
	rate=$(basename ${res_file%.txt})
	echo "------------------------------------------------"
	echo " $rate from $number_rates rate categories"
	echo "------------------------------------------------"
	echo ""
	
	res_tree=$input_dir/tree_${rate}_table
	cat  $res_file | grep T2T | grep  Node | tr -d '*' > $res_tree

	outfile=$(echo ${tree_file}_extended.nex)
	Rscript $annotateTree $tree_file $res_tree $input_dir $outfile
	mv $new_dir/Tree_Rplot.pdf $new_dir/Tree_Rplot_${rate}.pdf
	
done
