#!/usr/bin/env bash

set -e

#-----------------------------------------------------------------
#get options and input handling
#example=/home/elgert/Desktop/Cassius/Satute2/test/visual_example
#input_dir=$example
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
visualTree=$fscripts/script-visualize-satute-tree.R
iqtree2=/home/elgert/IQ-TREE/iqtree-2.2.2.4-Linux/bin/iqtree2
#------------------------------------------------------

file=$(basename $input_dir/*.{phy,fasta})
name=$(echo ${file%.{phy,fasta}})
if [ -f $input_dir/${name}_resultsRate1.nwk ]
then
	rm -r  $input_dir/*.nwk $input_dir/*table $input_dir/*.nex $input_dir/*.pdf
fi

number_rates=0
for res in $input_dir/resultsRate*.txt
do
	number_rates=$(expr $number_rates + 1)
done

tree_file=$(basename $input_dir/*.treefile)
echo "---------------------------"
echo "Consider $treefile for $input_dir"
echo "---------------------------"
echo ""

for res_file in $input_dir/resultsRate*.txt
do
	rate=$(basename ${res_file%.txt})
	echo "------------------------------------------------"
	echo " $rate from $number_rates rate categories"
	echo "------------------------------------------------"
	echo ""
	
	res_tree=$input_dir/tree_${rate}_table
	cat  $res_file | grep T2T | grep  Node | tr -d '*' > $res_tree

	name=$(echo ${tree_file%.treefile})
	outfile=$(echo "${name%.phy}_${rate}")
	echo $outfile
	Rscript $annotateTree $input_dir/$tree_file $res_tree $input_dir $outfile
	Rscript $visualTree $input_dir/${outfile}.nwk $outfile $input_dir 0.01
	#mv $new_dir/Tree_Rplot.pdf $new_dir/Tree_Rplot_${rate}.pdf
	
done
