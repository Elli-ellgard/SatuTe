import os
from tests.scripts.satute_test_utils import (
    run_satute,
    print_test_name,
    create_destination_dir,
    copy_files_to_dest_dir,
    check_satute_files,
)

from tests.scripts.fixtures import *

"""TODO: The following tests throw all errors, not checking the existence of output files but the standard output """

def test_1(data_dir_path, iqtree, python, satute):
    source_path, msa, tree_file = data_dir_path    
    suffix = "TEST 1: no msa no dir"
    print_test_name(suffix)

    # create output directory
    dest_dir_path = create_destination_dir(source_path, suffix)

    # copy msa file to output directory
    files_to_copy = []
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
    assert excinfo.value.code == 1  # assuming exit code 1 for failure

    # check the files    
    assert not check_satute_files(msa, dest_dir_path, categories, alpha, asr)
    
def test_2(data_dir_path, iqtree, python, satute):
    source_path, msa, tree_file = data_dir_path    
    suffix = "TEST 2: only IQ-Tree and -alpha value and -asr, should break!"
    print_test_name(suffix)

    # create output directory
    dest_dir_path = create_destination_dir(source_path, suffix)

    # copy msa file to output directory
    files_to_copy = [msa, tree_file]
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
    assert excinfo.value.code == 1  # assuming exit code 1 for failure
    # check the files
    assert not check_satute_files(msa, dest_dir_path, categories, alpha, asr)

def test_3(data_dir_path, iqtree, python, satute):
    source_path, msa, tree_file = data_dir_path    
    suffix = "TEST 3: msa  tree  model boot"
    print_test_name(suffix)
    # create output directory
    dest_dir_path = create_destination_dir(source_path, suffix)

    # copy msa file to output directory
    files_to_copy = [msa, tree_file]
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
                "-msa",
                os.path.join(dest_dir_path, msa),
                "-tree",
                os.path.join(dest_dir_path, tree_file),
                "-model",
                "JC+G2",
                "-boot",
                "100",
                "-alpha",
                alpha,
                "-asr",
            ]
        )

    # Verify the exit code if needed
    assert excinfo.value.code == 1  # assuming exit code 1 for failure
    # check the files
    assert not check_satute_files(msa, dest_dir_path, categories, alpha, asr)

def test_4(data_dir_path, iqtree, python, satute):
    source_path, msa, treefile = data_dir_path    
    suffix = "TEST 4: msa  tree model ufboot"
    print_test_name(suffix)

    # create output directory
    dest_dir_path = create_destination_dir(source_path, suffix)

    # copy msa file to output directory
    files_to_copy = [msa, treefile]
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
                "-msa",
                os.path.join(dest_dir_path, msa),
                "-tree",
                os.path.join(dest_dir_path, treefile),
                "-model",
                "JC+G2",
                "-ufboot",
                "1000",
                "-alpha",
                alpha,
                "-asr",
            ]
        )

    # Verify the exit code if needed
    assert excinfo.value.code == 1  # assuming exit code 1 for failure
    # check the files
    assert not check_satute_files(msa, dest_dir_path, categories, alpha, asr)



