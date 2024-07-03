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
    suffix = "MODEL TEST 1: only model"

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
    assert excinfo.value.code != 0  # assuming exit code 1 for failure
   # check the files
    assert not check_iqtree_files_exist(msa, dest_dir_path, []) , "IQ-Tree files check failed: Required files are missing or not created."
    assert not check_satute_files(msa, dest_dir_path, categories, alpha, asr), "Satute files check failed: Required files are missing or not created."



def test_2(data_dir_path, iqtree, python, satute):
    source_path, msa, _ = data_dir_path    
    suffix = "MODEL TEST 2: msa model no heterogeneity"

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

    # check the filesw
    assert check_iqtree_files_exist(msa, dest_dir_path, []) , "IQ-Tree files check failed: Required files are missing or not created."
    assert check_satute_files(msa, dest_dir_path, categories, alpha, asr), "Satute files check failed: Required files are missing or not created."

def test_3(data_dir_path, iqtree, python, satute):
    source_path, msa, _ = data_dir_path
    suffix = "MODEL TEST 3: msa model heterogeneity"

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
    assert check_iqtree_files_exist(msa, dest_dir_path, []) , "IQ-Tree files check failed: Required files are missing or not created."
    assert check_satute_files(msa, dest_dir_path, categories, alpha, asr), "Satute files check failed: Required files are missing or not created."

def test_4(data_dir_path, iqtree, python, satute):
    source_path, msa, _ = data_dir_path
    suffix = "MODEL TEST 4: msa model  heterogeneity rate number set"

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
    assert check_iqtree_files_exist(msa, dest_dir_path, []) , "IQ-Tree files check failed: Required files are missing or not created."
    assert check_satute_files(msa, dest_dir_path, categories, alpha, asr), "Satute files check failed: Required files are missing or not created."

def test_4_plus_i_g_four(data_dir_path, iqtree, python, satute):
    source_path, msa, _ = data_dir_path
    suffix = "MODEL TEST 4: msa model  heterogeneity rate number set Plus I Plus G4"

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
            "JC+I+G4",
            "-alpha",
            "0.05",
            "-iqtree",
            iqtree,
        ]
    )

    # check the files
    assert check_iqtree_files_exist(msa, dest_dir_path, []) , "IQ-Tree files check failed: Required files are missing or not created."
    assert check_satute_files(msa, dest_dir_path, categories, alpha, asr), "Satute files check failed: Required files are missing or not created."



def test_4_plus_i(data_dir_path, iqtree, python, satute):
    source_path, msa, _ = data_dir_path
    suffix = "MODEL TEST 4: msa model  heterogeneity rate number set plus I"
    
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
            "JC+I",
            "-alpha",
            "0.05",
            "-iqtree",
            iqtree,
        ]
    )

    # check the files
    assert check_iqtree_files_exist(msa, dest_dir_path, []) , "IQ-Tree files check failed: Required files are missing or not created."
    assert check_satute_files(msa, dest_dir_path, categories, alpha, asr), "Satute files check failed: Required files are missing or not created."

def test_5(data_dir_path, iqtree, python, satute):
    source_path, msa, _ = data_dir_path
    suffix = "MODEL TEST 5: msa model with parameter no heterogeneity"

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
    assert check_iqtree_files_exist(msa, dest_dir_path, []) , "IQ-Tree files check failed: Required files are missing or not created."
    assert check_satute_files(msa, dest_dir_path, categories, alpha, asr), "Satute files check failed: Required files are missing or not created."

def test_6(data_dir_path, iqtree, python, satute):
    source_path, msa, _ = data_dir_path
    suffix = "MODEL TEST 6: msa model with parameter gamma shape parameter"
    
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
    assert check_iqtree_files_exist(msa, dest_dir_path, []), "IQ-Tree files check failed: Required files are missing or not created."
    assert check_satute_files(msa, dest_dir_path, categories, alpha, asr), "Satute files check failed: Required files are missing or not created."

def test_7(data_dir_path, iqtree, python, satute):
    source_path, msa, _ = data_dir_path    
    suffix = "MODEL TEST 7: msa model with parameter free rate"

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
    assert check_iqtree_files_exist(msa, dest_dir_path, []), "IQ-Tree files check failed: Required files are missing or not created."
    assert check_satute_files(msa, dest_dir_path, categories, alpha, asr),  "Satute files check failed: Required files are missing or not created."

def test_8(data_dir_path, iqtree, python, satute):
    source_path, msa, _ = data_dir_path
    suffix = "MODEL TEST 8: msa model boot"
    
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

    assert check_iqtree_files_exist(msa, dest_dir_path, ["boot"]), "IQ-Tree files check failed: Required files are missing or not created."
    assert check_satute_files(msa, dest_dir_path, categories, alpha, asr), "Satute files check failed: Required files are missing or not created."

def test_9(data_dir_path, iqtree, python, satute):
    source_path, msa, _ = data_dir_path
    suffix = "MODEL TEST 9: msa model ufboot"

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

    # Create output directory
    dest_dir_path = create_destination_dir(source_path, suffix)

    # Copy msa file to output directory
    files_to_copy = [msa]
    copy_files_to_dest_dir(source_path, dest_dir_path, files_to_copy)

    categories = [3]
    alpha = str(0.05)
    asr = True

    try:
        # Satute run
        run_satute(
            [
                "-msa",
                os.path.join(Path(dest_dir_path).absolute(), msa),
                "-model",
                "JC+G4",
                "-alpha",
                alpha,
                "-category",
                "3",
                "-iqtree",
                iqtree,
                "-asr",
            ]
        )

        # Assertions for IQTree and Satute files
        assert check_iqtree_files_exist(msa, Path(dest_dir_path).absolute(), []), "IQ-Tree files check failed: Required files are missing or not created."
        assert check_satute_files(msa, Path(dest_dir_path).absolute(), categories, alpha, asr), "Satute files check failed: Required files are missing or not created."

    except SystemExit as e:
        # Check if the exit code matches the expected error scenario for empty category rate
        if e.code == 1:
            # Log and handle specific error message if needed
            print("Caught expected SystemExit due to empty category rate. Test will still pass.")
        else:
            # If it's a different SystemExit code, raise it again
            raise
    except Exception as e:
        # Catch other unexpected exceptions and fail the test
        pytest.fail(f"An unexpected exception occurred: {e}")


def test_11(data_dir_path, iqtree, python, satute):
    source_path, msa, _ = data_dir_path
    suffix = "MODEL TEST 11: msa model invariant sites"
    
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
            "GTR+I",
            "-alpha",
            "0.05",
            "-iqtree",
            iqtree,
        ]
    )

    # check the files
    assert check_iqtree_files_exist(msa, dest_dir_path, []), "IQ-Tree files check failed: Required files are missing or not created."
    assert check_satute_files(msa, dest_dir_path, categories, alpha, asr), "Satute files check failed: Required files are missing or not created."


def test_12(data_dir_path, iqtree, python, satute):
    source_path, msa, _ = data_dir_path
    suffix = "MODEL TEST 12: msa only one rate category"
    
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
            "GTR+G1",
            "-alpha",
            "0.05",
            "-iqtree",
            iqtree,
        ]
    )

    # check the files
    assert check_iqtree_files_exist(msa, dest_dir_path, []), "IQ-Tree files check failed: Required files are missing or not created."
    assert check_satute_files(msa, dest_dir_path, categories, alpha, asr), "Satute files check failed: Required files are missing or not created."
