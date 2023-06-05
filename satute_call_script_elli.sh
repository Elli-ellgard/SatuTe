iqtree=/home/elgert/IQ-TREE/iqtree-2.2.2.4-Linux/bin/iqtree2
rm -r test/octo-kraken-msa-test/clades        
python3 remove_files_except.py -m  
# python satute_cli.py -dir /Users/berksakalli/Projects/Satute/test/octo-kraken-msa-test -iqtree iqtree -model GTR+F
python3 satute_cli.py -iqtree $iqtree  -msa ./test/octo-kraken-msa-test/example.phy 
# rm -r test/octo-kraken-msa-test/clades        
