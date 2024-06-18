from logging import Logger
from Bio.Align import MultipleSeqAlignment
from pathlib import Path
from typing import List
import os
import argparse


def validate_dir_conflicts(input_args: List[argparse.Namespace]):
    """
    Ensure '-dir' is not used with incompatible options.

    Args:
        input_args: Parsed command-line arguments.

    Raises:
        ValueError: If '-dir' is combined with '-msa', '-tree', '-model', '-ufboot', '-boot', or '-add_iqtree_options'.
    """
    if input_args.dir and any(
        getattr(input_args, opt)
        for opt in ["msa", "tree", "model", "ufboot", "boot", "add_iqtree_options"]
    ):
        raise ValueError(
            "Error: The '-dir' option cannot be used with '-msa', '-tree', '-model', '-ufboot', '-boot', or '-add_iqtree_options'. "
            "Choose either '-dir' or the other options."
        )

        
def validate_msa_presence(input_args: List[argparse.Namespace]):
    """
    Ensure either 'msa' or 'dir' is specified, not both.

    Args:
        input_args: Parsed command-line arguments.

    Raises:
        ValueError: If neither or both 'msa' and 'dir' are specified.
    """
    if not input_args.dir and not input_args.msa:
        raise ValueError("Error: An MSA file must be specified when '-dir' is not used.")
    if input_args.dir and input_args.msa:
        raise ValueError("Error: The 'msa' and 'dir' options cannot be used together.")
    
def validate_tree_and_model(input_args: List[argparse.Namespace]):
    """
    Ensure a model is specified when using a tree file.

    Args:
        input_args: Parsed command-line arguments.

    Raises:
        ValueError: If a tree file is used without specifying a model.
    """
    if input_args.tree and not input_args.model:
        raise ValueError("Error: A model must be specified when using a tree file.")


def validate_boot_options(input_args: List[argparse.Namespace]):
    """
    Ensure bootstrapping options are not used with a tree file.

    Args:
        input_args: Parsed command-line arguments.

    Raises:
        ValueError: If '-ufboot' or '-boot' is used with '-tree'.
    """
    if input_args.tree and (input_args.ufboot or input_args.boot):
        raise ValueError("Error: The '-ufboot' or '-boot' options cannot be used with '-tree'.")

def validate_category_range(input_args : List[argparse.Namespace], number_rates: int):
    """
    Validates the category value against the allowed range.
    
    Args:
    - input_args: Object containing model and category attributes.
    - number_rates: Maximum number of rates defining the category range.

    Raises:
    - ValueError: If the category is out of the valid range.
    """
    if input_args.model and input_args.category is not None:
        if not (1 <= input_args.category <= number_rates):
            raise ValueError(
                f"Invalid category: {input_args.category}. "
                f"The category must be between 1 and {number_rates}, inclusive. "
                "Please choose a valid category index."
            )

def validate_and_set_rate_category(input_category: int, number_rates: int, logger: Logger)->str:
    """
    Validate and set the rate category.

    Args:
        input_category: Category input from the user.
        number_rates: Number of rates from the model.
        logger: Logger for error messages.

    Raises:
        ValueError: If the category is out of the valid range.
    """
    if not 1 <= input_category <= number_rates:
        logger.error(f"Error: Chosen category '{input_category}' is out of the valid range of {number_rates}.")
        raise ValueError(f"Error: Chosen category '{input_category}' is out of the valid range of {number_rates}.")
    return str(input_category)

def validate_and_check_rate_categories(categorized_sites: dict[str, MultipleSeqAlignment], chosen_category: int, logger: Logger):
    """
    Validates rate categories and ensures that the chosen category is not empty.
    
    Args:
    - categorized_sites (dict): Dictionary where keys are rate categories and values are alignments (lists).
    - chosen_category (int): The category to validate.
    - logger (Logger): Logger for warning messages.

    Raises:
    - ValueError: If the chosen category is empty.
    """
    for rate, alignment in categorized_sites.items():
        if len(alignment) == 0:
            logger.warning(f"Skipping empty rate category '{rate}'.")
        if chosen_category and len(alignment) == 0 and str(chosen_category) in rate:
            raise ValueError(
                f"Chosen category rate '{chosen_category}' is empty. "
                "Choose a different category with assigned sites."
            )

#######  INPUT Validation ##########       
                
def valid_directory(path: Path)->Path:
    """
    Custom type function for argparse - checks if the provided path is a valid directory,
    is not empty, and contains a .iqtree file and at least one file with specified suffixes.

    Args:
    - path (str): Directory path to be validated.

    Returns:
    - pathlib.Path: Validated Path object.

    Raises:
    - argparse.ArgumentTypeError: If the provided path is not a directory, is empty, or does not contain the required files (.iqtree and one of the specified suffixes).
    """
    msa_file_types = {".fasta", ".nex", ".phy", ".txt"}

    if not os.path.isdir(path):
        raise argparse.ArgumentTypeError(f"{path} is not a valid directory")

    directory_files = os.listdir(path)
    if not directory_files:
        raise argparse.ArgumentTypeError(f"{path} directory is empty")

    # Check for the presence of a .iqtree file in the directory
    if not any(file.endswith(".iqtree") for file in directory_files):
        raise argparse.ArgumentTypeError(
            f"No .iqtree file found in the directory {path}"
        )
        
    # Check for the presence of at least one file with a specified suffix
    if not any(
        file.endswith(suffix) for suffix in msa_file_types for file in directory_files
    ):
        suffixes_str = ", ".join(msa_file_types)
        raise argparse.ArgumentTypeError(
            f"No file with suffixes {suffixes_str} found in the directory {path}"
        )

    return Path(path)


def valid_file(path: Path)->Path:
    """
    Custom type function for argparse - checks if the provided path is a valid file.

    Args:
    - path (str): File path to be validated.

    Returns:
    - pathlib.Path: Validated Path object.

    Raises:
    - argparse.ArgumentTypeError: If the provided path is not a file.
    """
    if not os.path.isfile(path):
        raise argparse.ArgumentTypeError(f"{path} is not a valid file")
    return Path(path)


def valid_alpha(alpha: float) -> float:
    """
    Validates that the provided alpha value is a valid significance level for a statistical test.

    Args:
        alpha (float): Alpha value to be validated.

    Returns:
        float: Validated alpha value.

    Raises:
        argparse.ArgumentTypeError: If the alpha value is not between 0 and 1 exclusive.
    """
    try:
        alpha = float(alpha)
        if not (0 <= alpha < 1):
            raise ValueError(
                f"Invalid alpha value '{alpha}'. The alpha value must be between 0 and 1, exclusive. "
                "Please provide a valid alpha value for the significance level.")
    except ValueError as e:
        # Raise an ArgumentTypeError for argparse compatibility
        raise argparse.ArgumentTypeError(e)
    return alpha