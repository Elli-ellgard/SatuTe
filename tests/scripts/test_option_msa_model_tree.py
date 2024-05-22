import os
from pathlib import Path
from tests.scripts.satute_test_utils import (
    run_external_command,
    print_test_name,
    create_destination_dir,
    copy_files_to_dest_dir,
    check_iqtree_files_exist,
    check_satute_files,
    print_colored_message,
)

from tests.scripts.fixtures import *

def test_1(data_dir_path, iqtree, python, satute):
    source_path, msa, treefile = data_dir_path
    
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
    if not check_satute_files(msa, dest_dir_path, categories, alpha, asr):
        print_colored_message(f"{suffix} was successful", "32")
    else:
        print_colored_message(f"{suffix} failed", "31")
