from test_option_msa import test_option_msa
from test_option_model import test_option_model
from test_other_options import test_other_options


def main():

    # set paths to IQ-TREE and Python executable
    path_iqtree = "iqtree2"
    path_python = "python3"
    path_satute = "../../satute_cli.py"

    # smallest toy example
    data_dir_path = "../data/data_dna/toy_example_JC"
    msa = "toy_example_ntaxa_7_run_5-alignment.phy"

    output_dir_path =  "../test_results/"

    test_option_msa(path_iqtree, path_python, path_satute, data_dir_path, msa, output_dir_path)
    test_other_options(path_iqtree, path_python, path_satute, data_dir_path, msa, output_dir_path)
    test_option_model(path_iqtree, path_python, path_satute, data_dir_path, msa, output_dir_path)



if __name__ == "__main__":
    main()
