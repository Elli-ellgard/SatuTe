iqtree=iqtree2
#rm -r test/octo-kraken-msa-test/clades  
#rm -r test/octo-kraken-msa-test/subsequences      
#python3 remove_files_except.py -m  
#python3 satute_cli.py -dir /home/elgert/Desktop/Cassius/octo-kraken -iqtree iqtree 
#python3 satute_cli.py -iqtree $iqtree -dir /home/elgert/Desktop/Cassius/octo-kraken/test/octo-kraken-msa-test/ -model GTR+G -alpha 0.05
#python3 satute_cli.py -iqtree $iqtree -dir /home/elgert/Desktop/Neue_Daten_2023_05_12/results/pbmc3k/subsampled_alignment_pbmc3k_100hvg_genewisenormed/satute/

##  Test cases from Clemens
DIR=./Clemens/example_1
msa=ENSG00000087460_GNAS.fasta

#DIR=./Clemens/example_2
#msa=sim-JC+G-alpha1.2-taxa64-len1000bp-bla0.01-blb0.2-blc0.1-rep01.phy

DIR=./Clemens/example_3
msa=sim-JC+G-AC1-AG1-AT1-CG1-CT1-GT1-alpha1.2-taxa64-len1000bp-bla0.01-blb0.8-blc0.2-rep01.fasta

#DIR=./Clemens/example_4
#msa=example.txt

#DIR=./Clemens/example_sym_1
#msa=ENSG00000119574_ZBTB45.fasta

#DIR=./Clemens/example_sym_2
#msa=ENSG00000138316_ADAMTS14.fasta

#DIR=./Clemens/toy_example_GTR+G4
#msa=toy_example_ntaxa_7_run_1-alignment.phy

#DIR=./Clemens/toy_example_JC
#msa=toy_example_ntaxa_7_run_5-alignment.phy

PDIR=$(dirname $DIR)
if [ -e $DIR/${msa}.iqtree ]; then
    mv $DIR/$msa $PDIR
    rm -r $DIR
    mkdir $DIR
    mv $PDIR/$msa $DIR
fi
if [ -d "$DIR/clades" ]; then
    rm -r $DIR/clades
fi
if [ -d "$DIR/subsequences" ]; then
    rm -r $DIR/subsequences
fi

python3 satute_cli.py -iqtree $iqtree -msa $DIR/$msa -alpha 0.05

#python3 satute_cli.py -iqtree $iqtree -msa $alignment -tree random_generated_tree.tree -model $model -alpha 0.05


#cd visualisation
#./script-annotate-trees-from-satute.sh ../$DIR    
#cd ..
