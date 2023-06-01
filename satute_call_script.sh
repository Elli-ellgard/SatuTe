rm -r test/octo-kraken-msa-test/clades        
python remove_files_except.py -m  
python satute_cli.py -dir ./test/octo-kraken-msa-test -iqtree iqtree -model GTR+F
rm -r test/octo-kraken-msa-test/clades        
