import os
from tests.scripts.fixtures import *
from satute.models.amino_acid_models import AA_STATE_FREQUENCIES, NOT_ACCEPTED_AA_MODELS

from tests.scripts.satute_test_utils import (
    run_satute,
    check_satute_files,
    create_destination_dir,
    copy_files_to_dest_dir,
    check_iqtree_files_exist,
)


@pytest.mark.parametrize("model", list(AA_STATE_FREQUENCIES.keys()))
def test_all_substitution_models(
    data_amino_acids_dir_path, iqtree, python, satute, model
):
    """
    Test function that iterates over all substitution models and verifies the Satute output.

    Args:
    - data_amino_acids_dir_path: Directory path containing amino acids data.
    - iqtree: Path to IQ-TREE executable.
    - python: Path to Python executable.
    - satute: Path to Satute executable.
    - model: Current substitution model to test.
    """
    source_path, msa, _ = data_amino_acids_dir_path
    suffix = f"Amino Acid MODEL TEST: {model}"

    # Create output directory
    dest_dir_path = create_destination_dir(source_path, suffix)

    # Copy msa file to output directory
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
            model,
            "-alpha",
            alpha,
            "-asr",
            "-iqtree",
            iqtree,
        ]
    )

    # Check the files
    assert check_iqtree_files_exist(
        msa, dest_dir_path, []
    ), "IQ-Tree files check failed: Required files are missing or not created."
    assert check_satute_files(
        msa, dest_dir_path, categories, alpha, asr
    ), "Satute files check failed: Required files are missing or not created."


@pytest.mark.parametrize("model", NOT_ACCEPTED_AA_MODELS)
def test_all_aa_not_substitution_models(
    data_amino_acids_dir_path, iqtree, python, satute, model
):
    """
    Test function that iterates over all substitution models and verifies the Satute output.

    Args:
    - data_amino_acids_dir_path: Directory path containing amino acids data.
    - iqtree: Path to IQ-TREE executable.
    - python: Path to Python executable.
    - satute: Path to Satute executable.
    - model: Current substitution model to test.
    """
    source_path, msa, _ = data_amino_acids_dir_path
    suffix = f"Amino Acid MODEL TEST: {model}"

    # Create output directory
    dest_dir_path = create_destination_dir(source_path, suffix)

    # Copy msa file to output directory
    files_to_copy = [msa]
    copy_files_to_dest_dir(source_path, dest_dir_path, files_to_copy)

    alpha = str(0.05)

    # Capture sys.exit using pytest.raises
    with pytest.raises(SystemExit) as excinfo:
        # Satute run
        run_satute(
            [
                "-msa",
                os.path.join(dest_dir_path, msa),
                "-model",
                model,
                "-alpha",
                alpha,
                "-asr",
                "-iqtree",
                iqtree,
            ]
        )
    # Verify the exit code if needed
    assert excinfo.value.code != 0  # assuming exit code 1 for failure
