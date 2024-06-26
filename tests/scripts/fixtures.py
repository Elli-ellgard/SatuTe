import pytest


@pytest.fixture
def iqtree():
    return 'iqtree2'

@pytest.fixture
def python():
    return 'python3'

@pytest.fixture
def satute():
    return 'satute'

# AMINO ACID DATA
@pytest.fixture
def data_amino_acids_dir_path():
    return ("tests/data/data_aa/cbl3.5_A8+B8_100bp_LG", "ali-len100nt_cbl3.5.phy", "./simtree/simtree.tree")

@pytest.fixture
def dir_aa_path():
    return "tests/data/data_aa/cbl3.5_A8+B8_100bp_LG_iqtree_run", 

@pytest.fixture
def dir_category_aa_path():
    return "tests/data/data_aa/cbl3.5_A8+B8_100bp_LG+G_iqtree_run/", 

# Nucleotide DATA
@pytest.fixture
def data_dir_path():
    return ("tests/data/data_dna/toy_example_JC", "toy_example_ntaxa_7_run_5-alignment.phy", "tree_plain.treefile")

@pytest.fixture
def dir_path():
    return "tests/data/data_dna/toy_example_finished_iqtree_run", 


@pytest.fixture
def dir_path_categories():
    return "tests/data/data_dna/tree_cbl3.5_A8+B8_1500bp_JC+G-a0.8_iqtree_tree_run", 

@pytest.fixture
def dir_path_iqtree_files():
    return "tests/data/point_iqtree_files/", 

output_dir_path =  "./tests/"