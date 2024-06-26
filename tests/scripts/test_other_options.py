import os

from tests.scripts.fixtures import *

from tests.scripts.satute_test_utils import (
    print_test_name, 
    create_destination_dir,
    copy_files_to_dest_dir,
    check_iqtree_files_exist,
    check_satute_files,
    run_satute    
)


def test_1(data_dir_path, iqtree, python, satute):
    source_path, msa, _ = data_dir_path 
    suffix = "OTHER TEST 1: only msa IQ-Tree exists path to executable"

    # create output directory
    dest_dir_path = create_destination_dir(source_path, suffix)

    # copy msa file to output directory
    files_to_copy = [msa]
    copy_files_to_dest_dir(source_path, dest_dir_path, files_to_copy)

    categories = []
    alpha = str(0.05)
    asr = True

    # Satute run
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

    # check the files   
    assert check_iqtree_files_exist(msa, dest_dir_path, ["m"] )
    assert check_satute_files(msa, dest_dir_path, categories, alpha, asr)
    
def test_2(data_dir_path, iqtree, python, satute):
    source_path, msa, _ = data_dir_path 
    suffix = "OTHER TEST 2: only msa iqtree not exists"

    # create output directory
    dest_dir_path = create_destination_dir(source_path, suffix)

    # copy msa file to output directory
    files_to_copy = [msa]
    copy_files_to_dest_dir(source_path, dest_dir_path, files_to_copy)

    categories = []
    alpha = str(0.05)
    asr = True

    # Satute run
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

    
    # check the files   
    assert check_iqtree_files_exist(msa, dest_dir_path, ["m"]), "IQTree files check failed: Required files are missing or not created."
    assert check_satute_files(msa, dest_dir_path, categories, alpha, asr), "Satute files check failed: Required files are missing or not created."

def test_3(data_dir_path, iqtree, python, satute):
    source_path, msa, _ = data_dir_path 
    suffix = "OTHER TEST 3: no msa no dir"

    # create output directory
    dest_dir_path = create_destination_dir(source_path, suffix)

    # copy msa file to output directory
    files_to_copy = [msa]
    copy_files_to_dest_dir(source_path, dest_dir_path, files_to_copy)

    categories = []
    alpha = str(0.05)
    asr = True

    # Capture sys.exit using pytest.raises
    with pytest.raises(SystemExit) as excinfo:
        run_satute(
            [
                "-iqtree",
                iqtree,
                "-alpha",
                alpha,
                "-asr",
            ]
        )
    # Verify the exit code if needed
    assert excinfo.value.code != 0  # assuming exit code 1 for failure

    # check the files   
    assert not check_satute_files(msa, dest_dir_path, categories, alpha, asr) , "Satute files check failed: Required files are missing or not created."

def test_4(data_dir_path, iqtree, python, satute):
    source_path, msa, _ = data_dir_path 
    suffix = "OTHER TEST 4: alpha greater 1"
 
    # create output directory
    dest_dir_path = create_destination_dir(source_path, suffix)

    # copy msa file to output directory
    files_to_copy = [msa]
    copy_files_to_dest_dir(source_path, dest_dir_path, files_to_copy)

    categories = []
    alpha = str(2.0)
    asr = True

   # Capture sys.exit using pytest.raises
    with pytest.raises(SystemExit) as excinfo:
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
    # Verify the exit code if needed
    assert excinfo.value.code != 0  # assuming exit code 1 for failure

    # check the files   
    assert  not check_iqtree_files_exist(msa, dest_dir_path, ["m"]) , "IQTree files check failed: Required files are missing or not created."
    assert not check_satute_files(msa, dest_dir_path, categories, alpha, asr), "Satute files check failed: Required files are missing or not created."

def test_5(data_dir_path, iqtree, python, satute):
    source_path, msa, _ = data_dir_path 
    suffix = "CATEGORY TEST 1: msa model category exists"

    # create output directory
    dest_dir_path = create_destination_dir(source_path, suffix)

    # copy msa file to output directory
    files_to_copy = [msa]
    copy_files_to_dest_dir(source_path, dest_dir_path, files_to_copy)

    categories = [4]
    alpha = str(0.05)
    asr = True

    # Satute run
    run_satute(
        [
            "-msa",
            os.path.join(dest_dir_path, msa),
            "-model",
            "JC+G4",
            "-alpha",
            alpha,
            "-category",
            "4",
            "-iqtree",
            iqtree,
            "-asr"
        ]
    )

    
    # check the files   
    assert check_iqtree_files_exist(msa, dest_dir_path, []), "IQTree files check failed: Required files are missing or not created."
    assert check_satute_files(msa, dest_dir_path, categories, alpha, asr), "Satute files check failed: Required files are missing or not created."

def test_6(data_dir_path, iqtree, python, satute):
    source_path, msa, _ = data_dir_path 
    suffix = "CATEGORY TEST 2: msa model category not exists"

    # create output directory
    dest_dir_path = create_destination_dir(source_path, suffix)

    # copy msa file to output directory
    files_to_copy = [msa]
    copy_files_to_dest_dir(source_path, dest_dir_path, files_to_copy)

    categories = [2]
    alpha = str(0.05)
    asr = True


    # Capture sys.exit using pytest.raises
    with pytest.raises(SystemExit) as excinfo:
        # Satute run
        run_satute(
            [
                "-msa",
                os.path.join(dest_dir_path, msa),
                "-model",
                "JC+G4",
                "-category",
                "2",
                "-iqtree",
                iqtree,
                "-asr"
            ]
        )

    # Verify the exit code if needed
    assert excinfo.value.code != 0  # assuming exit code 1 for failure
    # check the files    
    assert not check_satute_files(msa, dest_dir_path, categories, alpha, asr)
