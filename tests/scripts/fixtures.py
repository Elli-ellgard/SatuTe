import pytest

@pytest.fixture
def iqtree():
    return 'iqtree'

@pytest.fixture
def python():
    return 'python3'

@pytest.fixture
def satute():
    return 'satute'

@pytest.fixture
def data_dir_path():
    return ("tests/data/data_dna/toy_example_JC", "toy_example_ntaxa_7_run_5-alignment.phy", "tree_plain.treefile")


output_dir_path =  "./tests/"