import os
from pathlib import Path
from tests.scripts.fixtures import *
from tests.scripts.satute_test_utils import (
    create_destination_dir,
    copy_files_to_dest_dir,
    check_satute_files,
    run_satute,
    delete_taxon_from_msa,
    delete_taxon_from_tree,
)


def test_1(data_dir_path, iqtree, python, satute):
    source_path, msa, treefile = data_dir_path
    
    suffix = "TEST 1: msa + tree + model"

    # create output directory
    dest_dir_path = create_destination_dir(source_path, suffix)

    # copy msa file to output directory
    files_to_copy = [msa, treefile]
    copy_files_to_dest_dir(source_path, dest_dir_path, files_to_copy)

    categories = [1,2]
    alpha = str(0.05)
    asr = True

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
            "-alpha",
            alpha,
            "-asr",
        ]
    )

    assert check_satute_files(msa, dest_dir_path, categories, alpha, asr)    


def test_2(data_dir_path, iqtree, python, satute):
    source_path, msa, tree_file = data_dir_path
    
    suffix = "TEST 1: msa + tree + model"

    # create output directory
    dest_dir_path = create_destination_dir(source_path, suffix)

    # copy msa file to output directory
    files_to_copy = [msa, tree_file]
    copy_files_to_dest_dir(source_path, dest_dir_path, files_to_copy)

    alpha = str(0.05)

    # Capture sys.exit using pytest.raises
    with pytest.raises(SystemExit) as excinfo:
        run_satute(
            [
                "-iqtree",
                iqtree,
                "-msa",
                os.path.join(dest_dir_path, msa),
                "-tree",
                os.path.join(dest_dir_path),
                "-model",
                "JC+G2",
                "-alpha",
                alpha,
                "-asr",
            ]
        )
                
    # Verify the exit code if needed
    assert excinfo.value.code != 0  # assuming exit code 1 for failure

def test_3(data_dir_path, iqtree, python, satute):
    source_path, msa, tree_file = data_dir_path
    
    suffix = "TEST 3: msa + tree + model + deleted taxon"

    # create output directory
    dest_dir_path = create_destination_dir(source_path, suffix)

    # copy msa file to output directory
    files_to_copy = [msa, tree_file]
    copy_files_to_dest_dir(source_path, dest_dir_path, files_to_copy)

    categories = [1,2]
    alpha = str(0.05)
    asr = True

    delete_taxon_from_msa(os.path.join(dest_dir_path, msa), os.path.join(dest_dir_path, msa),"t7","phylip-sequential")

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
                "-alpha",
                alpha,
                "-asr",
            ]
        )

    # Verify the exit code if needed
    assert excinfo.value.code != 0  # assuming exit code 1 for failure
                
def test_4(data_dir_path, iqtree, python, satute):
    source_path, msa, tree_file = data_dir_path
    
    suffix = "TEST 4: msa + tree + model + deleted taxon in MSA"

    # create output directory
    dest_dir_path = create_destination_dir(source_path, suffix)

    # copy msa file to output directory
    files_to_copy = [msa, tree_file]
    copy_files_to_dest_dir(source_path, dest_dir_path, files_to_copy)

    categories = [1,2]
    alpha = str(0.05)
    asr = True

    delete_taxon_from_msa(os.path.join(dest_dir_path, msa), os.path.join(dest_dir_path, msa),"t7","phylip-sequential")

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
                "-alpha",
                alpha,
                "-asr",
            ]
        )
                    
        # Verify the exit code if needed
    assert excinfo.value.code != 0  # assuming exit code 1 for failure

def test_5(data_dir_path, iqtree, python, satute):
    source_path, msa, tree_file = data_dir_path
    
    suffix = "TEST 5: msa + tree + model + deleted taxon in TREE"

    # create output directory
    dest_dir_path = create_destination_dir(source_path, suffix)

    # copy msa file to output directory
    files_to_copy = [msa, tree_file]
    copy_files_to_dest_dir(source_path, dest_dir_path, files_to_copy)

    categories = [1,2]
    alpha = str(0.05)
    asr = True

    delete_taxon_from_tree(os.path.join(dest_dir_path, tree_file), "t7")

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
                "-alpha",
                alpha,
                "-asr",
            ]
        )                    
        # Verify the exit code if needed
    assert excinfo.value.code != 0  # assuming exit code 1 for failure
