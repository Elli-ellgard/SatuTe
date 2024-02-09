from satute_test_utils import (
    run_external_command,
    print_test_name, 
    create_destination_dir,
    copy_files_to_dest_dir,
    check_iqtree_files_exist,
    check_satute_files,
    print_colored_message
    
)
import os
from pathlib import Path


def test_1(source_path, msa,treefile, iqtree, python, satute):
    suffix = "TEST 1: msa  tree  model"
    print_test_name(suffix)

    # create output directory
    dest_dir_path = create_destination_dir(source_path, suffix)

    # copy msa file to output directory
    files_to_copy = [msa, treefile]
    copy_files_to_dest_dir(source_path, dest_dir_path, files_to_copy)

    categories = []
    alpha = str(2.0)
    asr = True


    run_external_command(
        [
            python,
            satute,
            "-iqtree",
            iqtree,
            "-msa",
            os.path.join(dest_dir_path, msa),
            "-tree",
            os.path.join(dest_dir_path, treefile),
            "-model",
            "JC+G2",
            "-alpha",
            alpha,
            "-asr",
        ]
    )

    # check the files   
    if not  check_satute_files(msa, dest_dir_path, categories, alpha, asr):
        print_colored_message(f"{suffix} was successful", "32" )
    else: 
        print_colored_message(f"{suffix} failed", "31" )


def test_option_msa_model_tree(path_iqtree, path_python, path_satute, source_path, msa, treefile, results_path):

    print("")
    print_colored_message(" ============= Option MSA + MODEL + TREE .... ====================", "36")
    print("")


    test_1(source_path, msa, treefile, path_iqtree, path_python, path_satute)
    #test_3(source_path, msa, path_iqtree, path_python, path_satute)
    #test_4(source_path, msa, path_iqtree, path_python, path_satute)
    print("")


if __name__ == "__main__":
    # set paths to IQ-TREE and Python executable
    path_iqtree = "iqtree2"
    path_python = "python3"
    path_satute = "../../satute_cli.py"

    # smallest toy example
    data_dir_path = "../data/data_dna/toy_example_JC"
    msa = "toy_example_ntaxa_7_run_5-alignment.phy"
    treefile="tree_plain.treefile"

    output_dir_path =  "../test_results/"

    test_option_msa_model_tree(path_iqtree, path_python, path_satute, data_dir_path, msa, treefile, output_dir_path)
    
