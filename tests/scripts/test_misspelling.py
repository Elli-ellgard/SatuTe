import os
from pathlib import Path
from tests.scripts.fixtures import *
from tests.scripts.satute_test_utils import (
    run_satute,
    check_satute_files,    
    create_destination_dir,    
    copy_files_to_dest_dir,    
    check_iqtree_files_exist,
)


def test_1(data_dir_path, iqtree, python, satute):
    source_path, msa, _ = data_dir_path 
    suffix = "MISSPELL TEST 1: only msa IQ-TREE exists"

    # create output directory
    dest_dir_path = create_destination_dir(source_path, suffix)

    # copy msa file to output directory
    files_to_copy = [msa]
    copy_files_to_dest_dir(source_path, dest_dir_path, files_to_copy)

    categories = []
    alpha = str(0.05)
    asr = False

    # Satute run
    with pytest.raises(SystemExit) as excinfo:
        run_satute(
            [
                "-iqtree",
                iqtree,
                "-m",
                os.path.join(dest_dir_path, msa),
                "-alpha",
                alpha,
            ]
        )
    
    assert excinfo.value.code != 0  # assuming exit code 1 for failure

    # check the files   
    assert not check_iqtree_files_exist(msa, dest_dir_path, ["m"] )
    assert not check_satute_files(msa, dest_dir_path, categories, alpha, asr)

def test_2(data_dir_path, iqtree, python, satute):
    source_path, msa, _ = data_dir_path 
    suffix = "MISSPELL TEST 2: only msa"

    # create output directory
    dest_dir_path = create_destination_dir(source_path, suffix)

    # copy msa file to output directory
    files_to_copy = [msa]
    copy_files_to_dest_dir(source_path, dest_dir_path, files_to_copy)

    categories = []
    alpha = str(0.01)
    asr = False

    # Satute run
    with pytest.raises(SystemExit) as excinfo:
        run_satute(
            [
                "-iqtree2",
                iqtree,
                "-msa",
                os.path.join(dest_dir_path, msa),
                "-alpha",
                alpha,
            ]
        )

    assert excinfo.value.code != 0  # assuming exit code 1 for failure
    
    # check the files   
    assert not check_iqtree_files_exist(msa, dest_dir_path, ["m"]) 
    assert not check_satute_files(msa, dest_dir_path, categories, alpha, asr)

def test_3(data_dir_path, iqtree, python, satute):
    source_path, msa, _ = data_dir_path 
    suffix = "MSA TEST 3: only msa boot"

    # create output directory
    dest_dir_path = create_destination_dir(source_path, suffix)

    # copy msa file to output directory
    files_to_copy = [msa]
    copy_files_to_dest_dir(source_path, dest_dir_path, files_to_copy)

    categories = []
    alpha = str(0.02)
    asr = True

    with pytest.raises(SystemExit) as excinfo:
        run_satute(
            [
                "-iqtree",
                iqtree,
                "-msa",
                os.path.join(dest_dir_path, msa),
                "-booot",
                "100",
                "-alpha",
                alpha,
                "-asr",
            ]
        )

    assert excinfo.value.code != 0  # assuming exit code 1 for failure

    # check the files   
    assert not check_iqtree_files_exist(msa, dest_dir_path, ["boot","m"])  
    assert not check_satute_files(msa, dest_dir_path, categories, alpha, asr)
    