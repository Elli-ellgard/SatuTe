from tests.scripts.satute_test_utils import (
    run_external_command,
    run_satute,
    print_test_name,
    create_destination_dir,
    copy_files_to_dest_dir,
    check_iqtree_files_exist,
    check_satute_files,
)
import os
from pathlib import Path
import sys
from tests.scripts.fixtures import *


def test_1(data_dir_path, iqtree, python, satute):
    source_path, msa, _ = data_dir_path
    suffix = "MODEL TEST 1: only model"
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

    # Capture sys.exit using pytest.raises
    with pytest.raises(SystemExit) as excinfo:
        # Satute run
        run_satute(
            [
                "-model",
                "JC",
                "-alpha",
                "0.05",
                "-asr",
                "-iqtree",
                iqtree,
            ]
        )
    # Verify the exit code if needed
    assert excinfo.value.code == 1  # assuming exit code 1 for failure
   # check the files
    assert not check_iqtree_files_exist(msa, dest_dir_path, []) , "IQTree files check failed: Required files are missing or not created."
    assert not check_satute_files(msa, dest_dir_path, categories, alpha, asr), "Satute files check failed: Required files are missing or not created."

def test_2(data_dir_path, iqtree, python, satute):

    source_path, msa, _ = data_dir_path
    
    suffix = "MODEL TEST 2: msa model no heterogeneity"
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
    run_satute(
        [
            "-msa",
            os.path.join(dest_dir_path, msa),
            "-model",
            "JC",
            "-alpha",
            "0.05",
            "-asr",
            "-iqtree",
            iqtree,
        ]
    )

    # check the files
    assert check_iqtree_files_exist(msa, dest_dir_path, []) , "IQTree files check failed: Required files are missing or not created."
    assert check_satute_files(msa, dest_dir_path, categories, alpha, asr), "Satute files check failed: Required files are missing or not created."

def test_3(data_dir_path, iqtree, python, satute):
    
    source_path, msa, _ = data_dir_path
    suffix = "MODEL TEST 3: msa model heterogeneity"
    print_test_name(suffix)

    # create output directory
    dest_dir_path = create_destination_dir(source_path, suffix)

    # copy msa file to output directory
    files_to_copy = [msa]
    copy_files_to_dest_dir(source_path, dest_dir_path, files_to_copy)

    categories = [1,3,4]
    alpha = str(0.05)
    asr = True

    # Satute run
    run_satute(
        [
            "-msa",
            os.path.join(dest_dir_path, msa),
            "-model",
            "TIM+G",
            "-alpha",
            alpha,
            "-asr",
            "-iqtree",
            iqtree,
        ]
    )
    
    # check the files
    assert check_iqtree_files_exist(msa, dest_dir_path, []) , "IQTree files check failed: Required files are missing or not created."
    assert check_satute_files(msa, dest_dir_path, categories, alpha, asr), "Satute files check failed: Required files are missing or not created."

def test_4(data_dir_path, iqtree, python, satute):
    source_path, msa, _ = data_dir_path
    suffix = "MODEL TEST 4: msa model  heterogeneity rate number set"
    print_test_name(suffix)

    # create output directory
    dest_dir_path = create_destination_dir(source_path, suffix)

    # copy msa file to output directory
    files_to_copy = [msa]
    copy_files_to_dest_dir(source_path, dest_dir_path, files_to_copy)

    categories = [1, 4]
    alpha = str(0.05)
    asr = False

    # Satute run
    run_satute(
        [
            "-msa",
            os.path.join(dest_dir_path, msa),
            "-model",
            "JC+G4",
            "-alpha",
            "0.05",
            "-iqtree",
            iqtree,
        ]
    )

    # check the files
    assert check_iqtree_files_exist(msa, dest_dir_path, []) , "IQTree files check failed: Required files are missing or not created."
    assert check_satute_files(msa, dest_dir_path, categories, alpha, asr), "Satute files check failed: Required files are missing or not created."

def test_5(data_dir_path, iqtree, python, satute):
    source_path, msa, _ = data_dir_path

    suffix = "MODEL TEST 5: msa model with parameter no heterogeneity"
    print_test_name(suffix)

    # create output directory
    dest_dir_path = create_destination_dir(source_path, suffix)

    # copy msa file to output directory
    files_to_copy = [msa]
    copy_files_to_dest_dir(source_path, dest_dir_path, files_to_copy)

    categories = []
    alpha = str(0.05)
    asr = False

    # Satute run
    run_satute(
        [
            "-msa",
            os.path.join(dest_dir_path, msa),
            "-model",
            "GTR{4.39/5.30/4.39/1.0/12.1}",
            "-alpha",
            "0.05",
            "-iqtree",
            iqtree,
        ]
    )

    # check the files
    assert check_iqtree_files_exist(msa, dest_dir_path, []) , "IQTree files check failed: Required files are missing or not created."
    assert check_satute_files(msa, dest_dir_path, categories, alpha, asr), "Satute files check failed: Required files are missing or not created."

def test_6(data_dir_path, iqtree, python, satute):
    source_path, msa, _ = data_dir_path
    suffix = "MODEL TEST 6: msa model with parameter gamma shape parameter"
    print_test_name(suffix)

    # create output directory
    dest_dir_path = create_destination_dir(source_path, suffix)

    # copy msa file to output directory
    files_to_copy = [msa]
    copy_files_to_dest_dir(source_path, dest_dir_path, files_to_copy)

    categories = [1, 2, 3, 4]
    alpha = str(0.05)
    asr = False

    # Satute run
    run_satute(
        [
            "-msa",
            os.path.join(dest_dir_path, msa),
            "-model",
            "TIM2{4.39/5.30/12.1}+G{0.7}",
            "-alpha",
            "0.05",
            "-iqtree",
            iqtree,
        ]
    )

    # check the files
    assert check_iqtree_files_exist(msa, dest_dir_path, []), "IQTree files check failed: Required files are missing or not created."
    assert check_satute_files(msa, dest_dir_path, categories, alpha, asr), "Satute files check failed: Required files are missing or not created."

def test_7(data_dir_path, iqtree, python, satute):
    source_path, msa, _ = data_dir_path
    
    suffix = "MODEL TEST 7: msa model with parameter free rate"
    print_test_name(suffix)

    # create output directory
    dest_dir_path = create_destination_dir(source_path, suffix)

    # copy msa file to output directory
    files_to_copy = [msa]
    copy_files_to_dest_dir(source_path, dest_dir_path, files_to_copy)

    categories = [1, 2]
    alpha = str(0.05)
    asr = False

    # Satute run
    run_satute(
        [
            "-msa",
            os.path.join(dest_dir_path, msa),
            "-model",
            "TIM2{4.39/5.30/12.1}+R2{0.7/0.1/0.3/0.9}",
            "-alpha",
            "0.05",
            "-iqtree",
            iqtree,
        ]
    )

    # check the files
    assert check_iqtree_files_exist(msa, dest_dir_path, []), "IQTree files check failed: Required files are missing or not created."
    assert check_satute_files(msa, dest_dir_path, categories, alpha, asr),  "Satute files check failed: Required files are missing or not created."

def test_8(data_dir_path, iqtree, python, satute):
    source_path, msa, _ = data_dir_path

    suffix = "MODEL TEST 8: msa model boot"
    print_test_name(suffix)

    # create output directory
    dest_dir_path = create_destination_dir(source_path, suffix)

    # copy msa file to output directory
    files_to_copy = [msa]
    copy_files_to_dest_dir(source_path, dest_dir_path, files_to_copy)

    categories = []
    alpha = str(0.05)
    asr = False

    # Satute run
    run_satute(
        [
            "-msa",
            os.path.join(dest_dir_path, msa),
            "-model",
            "GTR",
            "-alpha",
            "0.05",
            "-boot",
            "100",
            "-iqtree",
            iqtree,
        ]
    )

    assert check_iqtree_files_exist(msa, dest_dir_path, ["boot"]), "IQTree files check failed: Required files are missing or not created."
    assert check_satute_files(msa, dest_dir_path, categories, alpha, asr), "Satute files check failed: Required files are missing or not created."

def test_9(data_dir_path, iqtree, python, satute):
    source_path, msa, _ = data_dir_path
    suffix = "MODEL TEST 9: msa model ufboot"
    print_test_name(suffix)

    # create output directory
    dest_dir_path = create_destination_dir(source_path, suffix)

    # copy msa file to output directory
    files_to_copy = [msa]
    copy_files_to_dest_dir(source_path, dest_dir_path, files_to_copy)

    categories = [1, 2]
    alpha = str(0.05)
    asr = False


    # Satute run
    run_satute(
        [
            "-msa",
            os.path.join(dest_dir_path, msa),
            "-model",
            "GTR+G2",
            "-alpha",
            "0.05",
            "-ufboot",
            "1000",
            "-iqtree",
            iqtree,
        ]
    )

    # IQTree files check
    assert check_iqtree_files_exist(msa, dest_dir_path, ["ufboot"]), f"{suffix} was successful 32"
    # Satute files check
    assert check_satute_files(msa, dest_dir_path, categories, alpha, asr), f"{suffix} failed 31"

def test_10(data_dir_path, iqtree, python, satute):
    source_path, msa, _ = data_dir_path    
    suffix = "MODEL TEST 10: msa model category"
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
    run_satute(
        [
            "-msa",
            os.path.join(Path(dest_dir_path).absolute(), msa),
            "-model",
            "GTR+G4",
            "-alpha",
            "0.05",
            "-category",
            "3",
            "-iqtree",
            iqtree,
            "-asr",
        ]
    )

    assert check_iqtree_files_exist(msa, Path(dest_dir_path).absolute(), []), "IQTree files check failed: Required files are missing or not created."
    assert check_satute_files(msa, Path(dest_dir_path).absolute(), categories, alpha, asr), "Satute files check failed: Required files are missing or not created."

