from tests.scripts.satute_test_utils import (
    run_external_command,
    print_test_name, 
    create_destination_dir,
    copy_files_to_dest_dir,
    check_iqtree_files_exist,
    check_satute_files,
    print_colored_message
    
)
import os
from tests.scripts.fixtures import *

def test_1(data_dir_path, iqtree, python, satute):
    source_path, msa, _ = data_dir_path 
    suffix = "OTHER TEST 1: only msa iqtree exists  path to executable"
    print_test_name(suffix)
    print(f"IQ-Tree executeable path: {iqtree}")

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


def test_2(data_dir_path, iqtree, python, satute):
    source_path, msa, _ = data_dir_path 
    suffix = "OTHER TEST 2: only msa iqtree not exists"
    print_test_name(suffix)
    print(f"IQ-Tree executeable path: {iqtree}")

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
    if not check_iqtree_files_exist(msa, dest_dir_path, ["m"]) and not check_satute_files(msa, dest_dir_path, categories, alpha, asr):
        print_colored_message(f"{suffix} was successful", "32" )
    else: 
        print_colored_message(f"{suffix} failed", "31" )


def test_3(data_dir_path, iqtree, python, satute):
    source_path, msa, _ = data_dir_path 
    suffix = "OTHER TEST 3: no msa no dir"
    print_test_name(suffix)

    # create output directory
    dest_dir_path = create_destination_dir(source_path, suffix)

    # copy msa file to output directory
    files_to_copy = [msa]
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


def test_4(data_dir_path, iqtree, python, satute):
    source_path, msa, _ = data_dir_path 
    suffix = "OTHER TEST 4: alpha greater 1"
    print_test_name(suffix)

    # create output directory
    dest_dir_path = create_destination_dir(source_path, suffix)

    # copy msa file to output directory
    files_to_copy = [msa]
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
            "-alpha",
            alpha,
            "-asr",
        ]
    )

    # check the files   
    if  not check_iqtree_files_exist(msa, dest_dir_path, ["m"])  and   not check_satute_files(msa, dest_dir_path, categories, alpha, asr):
        print_colored_message(f"{suffix} was successful", "32" )
    else: 
        print_colored_message(f"{suffix} failed", "31" )

def test_5(data_dir_path, iqtree, python, satute):
    source_path, msa, _ = data_dir_path 
    suffix = "CATEGORY TEST 1: msa model category exists"
    print_test_name(suffix)

    # create output directory
    dest_dir_path = create_destination_dir(source_path, suffix)

    # copy msa file to output directory
    files_to_copy = [msa]
    copy_files_to_dest_dir(source_path, dest_dir_path, files_to_copy)

    categories = [3]
    alpha = str(0.05)
    asr = True

    # Satute run
    run_external_command(
        [
            python,
            satute,
            "-msa",
            os.path.join(dest_dir_path, msa),
            "-model",
            "JC+G4",
            "-alpha",
            "0.05",
            "-category",
            "3",
            "-iqtree",
            iqtree,
            "-asr"
        ]
    )

    
    # check the files   
    if check_iqtree_files_exist(msa, dest_dir_path, []) and  check_satute_files(msa, dest_dir_path, categories, alpha, asr):
        print_colored_message(f"{suffix} was successful", "32" )
    else: 
        print_colored_message(f"{suffix} failed", "31" )

def test_6(data_dir_path, iqtree, python, satute):
    source_path, msa, _ = data_dir_path 
    suffix = "CATEGORY TEST 2: msa model category not exists"
    print_test_name(suffix)

    # create output directory
    dest_dir_path = create_destination_dir(source_path, suffix)

    # copy msa file to output directory
    files_to_copy = [msa]
    copy_files_to_dest_dir(source_path, dest_dir_path, files_to_copy)

    categories = [2]
    alpha = str(0.05)
    asr = True

    # Satute run
    run_external_command(
        [
            python,
            satute,
            "-msa",
            os.path.join(dest_dir_path, msa),
            "-model",
            "JC+G4",
            "-alpha",
            "0.05",
            "-category",
            "2",
            "-iqtree",
            iqtree,
            "-asr"
        ]
    )

    
    # check the files   
    if check_iqtree_files_exist(msa, dest_dir_path, []) and  check_satute_files(msa, dest_dir_path, categories, alpha, asr):
        print_colored_message(f"{suffix} was successful", "32" )
    else: 
        print_colored_message(f"{suffix} failed", "31" )