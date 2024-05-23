import satute.cli

def test_simplet():
    satute.cli.main(['-msa', 'tests/data/data_dna/toy_example_JC/toy_example_ntaxa_7_run_5-alignment.phy', '-iqtree', 'iqtree'])