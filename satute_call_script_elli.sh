iqtree=/home/elgert/IQ-TREE/iqtree-2.2.2.4-Linux/bin/iqtree2
rm -r test/octo-kraken-msa-test/clades  
rm -r test/octo-kraken-msa-test/subsequences      
python3 remove_files_except.py -m  
#python3 satute_cli.py -dir /home/elgert/Desktop/Cassius/octo-kraken -iqtree iqtree 
python3 satute_cli.py -iqtree $iqtree -msa /home/elgert/Desktop/Cassius/octo-kraken/test/octo-kraken-msa-test/example.phy -model GTR+G -alpha 0.05
# rm -r test/octo-kraken-msa-test/clades        
#python3 satute_cli.py -iqtree $iqtree -dir /home/elgert/Desktop/Neue_Daten_2023_05_12/results/pbmc3k/subsampled_alignment_pbmc3k_100hvg_genewisenormed/satute/

