import pytest
import os
from satute.dna_models import LIE_DNA_MODELS
from tests.scripts.satute_test_utils import (
    run_satute,
    check_satute_files,    
    create_destination_dir,
    copy_files_to_dest_dir,
    check_iqtree_files_exist,
)
from tests.scripts.fixtures import *


@pytest.mark.parametrize("model", LIE_DNA_MODELS)
def test_all_substitution_models(data_dir_path, iqtree, model):
    """
    Test function that iterates over all substitution models and verifies the Satute output.
    
    Args:
    - data_dir_path: Directory path containing dna data.
    - iqtree: Path to IQ-TREE executable.
    - model: Current substitution model to test.
    """
    source_path, msa, _ = data_dir_path
    suffix = f"DNA MODEL TEST: {model}"

    # Create output directory
    dest_dir_path = create_destination_dir(source_path, suffix)

    # Copy msa file to output directory
    files_to_copy = [msa]
    copy_files_to_dest_dir(source_path, dest_dir_path, files_to_copy)

    categories = []
    alpha = str(0.05)
    asr = False

    # Satute run
    run_satute(
        [
            "-msa", os.path.join(dest_dir_path, msa),
            "-model", model,
            "-alpha", alpha,
            "-iqtree", iqtree,
        ]
    )

    # Check the files
    assert check_iqtree_files_exist(msa, dest_dir_path, []), "IQ-Tree files check failed: Required files are missing or not created."
    assert check_satute_files(msa, dest_dir_path, categories, alpha, asr), "Satute files check failed: Required files are missing or not created."

