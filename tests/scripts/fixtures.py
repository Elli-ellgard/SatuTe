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

@pytest.fixture
def data_dir_path():
    return ("tests/data/data_dna/toy_example_JC", "toy_example_ntaxa_7_run_5-alignment.phy", "tree_plain.treefile")


@pytest.fixture
def dir_path():
    return "tests/data/data_dna/toy_example_JC", 


@pytest.fixture
def dir_path_categories():
    return "tests/data/data_dna/tree_cbl3.5_A8+B8_1500bp_JC+G-a0.8_iqtree_tree_run", 


output_dir_path =  "./tests/"