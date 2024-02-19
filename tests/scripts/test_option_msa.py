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



def test_1(source_path, msa, iqtree, python, satute):
    suffix = "MSA TEST 1: only msa iqtree exists"
    print_test_name(suffix)

    # create output directory
    dest_dir_path = create_destination_dir(source_path, suffix)

    # copy msa file to output directory
    files_to_copy = [msa]
    copy_files_to_dest_dir(source_path, dest_dir_path, files_to_copy)

    categories = []
    alpha = str(0.05)
    asr = True

    # Satute run
    run_external_command(
        [
            python,
            satute,
            "-iqtree",
            iqtree,
            "-msa",
            os.path.join(dest_dir_path, msa),
            "-alpha",
            alpha,
            "-asr",
        ]
   )

    # check the files   
    if check_iqtree_files_exist(msa, dest_dir_path, ["m"] )  and check_satute_files(msa, dest_dir_path, categories, alpha, asr):
        print_colored_message(f"{suffix} was successful", "32" )
    else: 
        print_colored_message(f"{suffix} failed", "31" )


def test_2(source_path, msa, iqtree, python, satute):
    suffix = "MSA TEST 2: only msa ufboot"
    print_test_name(suffix)

    # create output directory
    dest_dir_path = create_destination_dir(source_path, suffix)

    # copy msa file to output directory
    files_to_copy = [msa]
    copy_files_to_dest_dir(source_path, dest_dir_path, files_to_copy)

    categories = []
    alpha = str(0.01)
    asr = False

    # Satute run
    run_external_command(
        [
            python,
            satute,
            "-iqtree",
            iqtree,
            "-msa",
            os.path.join(dest_dir_path, msa),
            "-ufboot",
            "1000",
            "-alpha",
            alpha,
        ]
    )

    
    # check the files   
    if check_iqtree_files_exist(msa, dest_dir_path, ["ufboot","m"]) and check_satute_files(msa, dest_dir_path, categories, alpha, asr):
        print_colored_message(f"{suffix} was successful", "32" )
    else: 
        print_colored_message(f"{suffix} failed", "31" )


def test_3(source_path, msa, iqtree, python, satute):
    suffix = "MSA TEST 3: only msa boot"
    print_test_name(suffix)

    # create output directory
    dest_dir_path = create_destination_dir(source_path, suffix)

    # copy msa file to output directory
    files_to_copy = [msa]
    copy_files_to_dest_dir(source_path, dest_dir_path, files_to_copy)

    categories = []
    alpha = str(0.02)
    asr = True


    run_external_command(
        [
            python,
            satute,
            "-iqtree",
            iqtree,
            "-msa",
            os.path.join(dest_dir_path, msa),
            "-boot",
            "100",
            "-alpha",
            alpha,
            "-asr",
        ]
    )

    # check the files   
    if check_iqtree_files_exist(msa, dest_dir_path, ["boot","m"])  and    check_satute_files(msa, dest_dir_path, categories, alpha, asr):
        print_colored_message(f"{suffix} was successful", "32" )
    else: 
        print_colored_message(f"{suffix} failed", "31" )


def test_option_msa(path_iqtree, path_python, path_satute, source_path, msa, results_path):

    print("")
    print_colored_message(" ============= MODI MSA....====================", "36")
    print("")


    test_1(source_path, msa, path_iqtree, path_python, path_satute)
    test_2(source_path, msa, path_iqtree, path_python, path_satute)
    test_3(source_path, msa, path_iqtree, path_python, path_satute)


if __name__ == "__main__":
    # set paths to IQ-TREE and Python executable
    path_iqtree = "iqtree"
    path_python = "python3"
    path_satute = "../../satute_cli.py"

    # smallest toy example
    data_dir_path = "../data/data_dna/toy_example_JC"
    msa = "toy_example_ntaxa_7_run_5-alignment.phy"

    output_dir_path =  "../test_results/"

    test_option_msa(path_iqtree, path_python, path_satute, data_dir_path, msa, output_dir_path)
    
