import os
from pathlib import Path
from tests.scripts.fixtures import *
from tests.scripts.satute_test_utils import (
    run_satute,
    create_destination_dir,
    copy_files_to_dest_dir,
    check_iqtree_files_exist,
    check_satute_files,
)


def test_1(data_amino_acids_dir_path, iqtree, python, satute):
    source_path, msa, _ = data_amino_acids_dir_path
    suffix = "Amino Acid MSA TEST 1: only msa IQ-TREE exists"

    # create output directory
    dest_dir_path = create_destination_dir(source_path, suffix)

    # copy msa file to output directory
    files_to_copy = [msa]
    copy_files_to_dest_dir(source_path, dest_dir_path, files_to_copy)

    categories = []
    alpha = str(0.05)
    asr = True
    
    run_satute(
        
             [
                 "-iqtree",
                 iqtree,
                 "-msa",
                 os.path.join(dest_dir_path, msa),
                 "-alpha",
                 alpha,
                 "-asr",
             ]
    )
    
        
    assert check_iqtree_files_exist(msa, Path(dest_dir_path).absolute(), []), "IQTree files check failed: Required files are missing or not created."
    assert check_satute_files(msa, Path(dest_dir_path).absolute(), categories, alpha, asr), "Satute files check failed: Required files are missing or not created."
    
    
def test_2(data_amino_acids_dir_path, iqtree, python, satute):
    source_path, msa, _ = data_amino_acids_dir_path 
    suffix = "Amino Acid MSA TEST 2: only msa ufboot"

    # create output directory
    dest_dir_path = create_destination_dir(source_path, suffix)

    # copy msa file to output directory
    files_to_copy = [msa]
    copy_files_to_dest_dir(source_path, dest_dir_path, files_to_copy)

    categories = []
    alpha = str(0.01)
    asr = False

    # Satute run
    run_satute(
        [
            "-iqtree",
            iqtree,
            "-msa",
            os.path.join(dest_dir_path, msa),
            "-ufboot",
            "1000",
            "-alpha",                            
            alpha,
            "-model",
            "LG"
        ]
    )

    
    # check the files   
    assert check_iqtree_files_exist(msa, dest_dir_path, ["ufboot","m"]) 
    assert check_satute_files(msa, dest_dir_path, categories, alpha, asr)

def test_3(data_amino_acids_dir_path, iqtree, python, satute):
    source_path, msa, _ = data_amino_acids_dir_path 
    suffix = "Amino Acid MSA TEST 3: only msa boot"

    # create output directory
    dest_dir_path = create_destination_dir(source_path, suffix)

    # copy msa file to output directory
    files_to_copy = [msa]
    copy_files_to_dest_dir(source_path, dest_dir_path, files_to_copy)

    categories = []
    alpha = str(0.02)
    asr = True


    run_satute(
        [
            "-iqtree",
            iqtree,            
            "-model",
            "DAYHOFF",
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
    assert check_iqtree_files_exist(msa, dest_dir_path, ["boot","m"])  
    assert check_satute_files(msa, dest_dir_path, categories, alpha, asr)
    
