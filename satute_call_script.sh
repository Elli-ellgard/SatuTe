rm -r test/octo-kraken-msa-test/subtrees
rm -r test/octo-kraken-msa-test/subsequences        
python remove_files_except.py -m  
python satute_cli.py -dir /Users/berksakalli/Projects/Satute/test/octo-kraken-msa-test -iqtree iqtree -model GTR+F
# python satute_cli.py -msa ./test/octo-kraken-msa-test/example.phy -iqtree iqtree -model GTR+F+G4
# python satute_cli.py -tree ./test/octo-kraken-msa-test-without-gamma/example.phy.tree -msa ./test/octo-kraken-msa-test-without-gamma/example.phy -iqtree iqtree -model GTR+F

# rm -r test/octo-kraken-msa-test/clades  
# rm -r test/octo-kraken-msa-test/subsequences        
