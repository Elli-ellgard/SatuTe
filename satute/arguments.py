# -*- coding: utf-8 -*-
import os
import argparse
import sys
from pathlib import Path


def valid_directory(path: Path):
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


def valid_file(path: Path):
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


def valid_alpha(alpha: float):
    """
    Custom type function for argparse - checks if the provided alpha value is valid.

    Args:
    - alpha (float): Alpha value to be validated.

    Returns:
    - float: Validated alpha value.

    Raises:
    - argparse.ArgumentTypeError: If the alpha value is not between 0 and 1 inclusive.
    """
    try:
        alpha = float(alpha)
        if not (0 <= alpha < 1):
            raise ValueError(f"The set alpha value of {alpha} is not valid")
    except ValueError as e:
        # Raise an ArgumentTypeError for argparse compatibility
        raise argparse.ArgumentTypeError(e)
    return alpha


ARGUMENT_LIST = [
    {
        "flag": "-dir",
        "help": (
            "Path to the input directory containing IQ-TREE output files. Use this option when you've already run IQ-TREE and want to avoid rerunning it. The directory should contain essential IQ-TREE output files including the .iqtree file, tree file(s), and possibly a .siteprob file."
        ),
        "metavar": "<directory_path>",
        "type": valid_directory,
    },
    {
        "flag": "-tree",
        "help": (
            "Path to the input tree file in Newick or Nexus format. This tree will be used as the basis for the saturation analysis."
        ),
        "metavar": "<tree_file_path>",
        "type": valid_file,
    },
    {
        "flag": "-msa",
        "help": (
            "Path to the Multiple Sequence Alignment (MSA) file you wish to analyze. The MSA can be in FASTA, NEXUS, PHYLIP, or TXT format."
        ),
        "metavar": "<msa_file_path>",
        "type": valid_file,
    },
    {
        "flag": "-iqtree",
        "help": (
            "Specifies the path to the IQ-TREE executable. If IQ-TREE is installed system-wide, just providing the executable name (`iqtree` or `iqtree2`) will suffice. Otherwise, give the complete path."
        ),
        "default": "iqtree2",
        "metavar": "<iqtree_path>",
        "type": Path,
    },
    {
        "flag": "-model",
        "help": (
            "Indicates the model of sequence evolution. Common models include `GTR`, `HKY`, etc. You can also specify rate heterogeneity and other model extensions, like `+G4` for gamma-distributed rates."
        ),
        "type": str,
        "metavar": "<evolution_model>",
    },
    {
        "flag": "-category",
        "help": (
            "Rate categories of interest. Relevant for models with gamma-distributed rate variations or FreeRate model. If the `-model` option includes rate variation (e.g., `+G4`), the `-category` should be a number between 1 and 4."
        ),
        "type": int,
        "metavar": "<rate_category>",
    },
    {
        "flag": "-ufboot",
        "help": (
            "Number of replicates for the ultrafast bootstrap analysis. Typically, a higher number like `1000` or `5000` is used. Ultrafast bootstrap provides rapid approximations to traditional bootstrap values."
        ),
        "type": int,
        "metavar": "<number_of_replicates>",
    },
    {
        "flag": "-boot",
        "help": (
            "Number of replicates for traditional bootstrap analysis. This also computes a Maximum Likelihood (ML) tree and a consensus tree. Common values are `1000` or `5000`."
        ),
        "type": int,
        "metavar": "<number_of_replicates>",
    },
    {
        "flag": "-alpha",
        "help": (
            "Significance level for the saturation test. A common threshold is `0.05`, indicating a 5% significance level. Lower values make the test more stringent."
        ),
        "type": valid_alpha,
        "default": 0.05,
        "metavar": "<significance_level>",
    },
    {
        "flag": "-edge",
        "help": (
            "Specify a branch or edge name to focus the analysis on. Useful when you want to check saturation on a specific branch."
        ),
        "type": str,
        "metavar": "<edge_name>",
    },
    {
        "flag": "-output_suffix",
        "help": "Specify a suffix for the output file.",
        "type": str,
        "metavar": "<output_suffix>",
        "default": "",
    },
    {
        "flag": "-add_iqtree_options",
        "help": "Specify additional options for the IQ-Tree run, if necessary.",
        "type": str,
        "metavar": "<additional_option>",
    },
    {
        "flag": "-asr",
        "help": "Write ancestral sequences (by empirical Bayesian method) for all nodes of the tree to a .asr.csv file.",
        "action": "store_true",
    },
    {
        "flag": "-category_assignment",
        "help": "Write assignment of the individual sites to the rate heterogeneity categories.",
        "action": "store_true",
    },
    {
        "flag": "-verbose",
        "help": "Enable verbose logging.",
        "action": "store_true",
    },
    {
        "flag": "-quiet",
        "help": "Even no warnings on the the terminal.",
        "action": "store_true",
    },
]
