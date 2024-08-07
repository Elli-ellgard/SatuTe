import logging
from pathlib import Path
from ete3 import Tree
from argparse import Namespace
from logging import Logger
from typing import List, Dict, Any
from satute.models.substitution_model import SubstitutionModel


def log_consider_iqtree_message(logger: Logger):
    logger.info("Running IQ-TREE")
    logger.warning(
        "Please consider for the analysis that IQ-Tree will be running with any advanced options."
    )
    logger.warning(
        "If specific options are required for the analysis, please run IQ-Tree separately."
    )


def construct_log_file_name(msa_file: Path, input_args: List[Namespace]):
    log_file = f"{msa_file.resolve()}_{input_args.alpha}.satute.log"
    if input_args.output_suffix:
        log_file = f"{msa_file.resolve()}_{input_args.alpha}_{input_args.output_suffix}.satute.log"
    return log_file


def setup_logging_configuration(
    logger: Logger, input_args: List[Namespace], msa_file: Path
):
    """
    Initializes the logging system for the application.
    Sets up two handlers:
    1. A file handler that always logs at the DEBUG level.
    2. A stream (console) handler that logs at the DEBUG level if verbose is true; otherwise, it logs at the WARNING level.
    The log file is named using the MSA file name, alpha value, and an output suffix.
    """
    # Logger level is set to DEBUG to capture all logs for the file handler
    logger.setLevel(logging.DEBUG)
    # File Handler - always active at DEBUG level
    log_file = construct_log_file_name(msa_file, input_args)
    file_handler = logging.FileHandler(log_file)
    # file_format = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    # file_handler.setFormatter(file_format)
    file_handler.setLevel(logging.DEBUG)  # Always log everything in file
    logger.addHandler(file_handler)

    # Set the default logging level
    if input_args.verbose:
        stream_level = logging.INFO
    elif input_args.quiet:
        stream_level = logging.CRITICAL
    else:
        stream_level = logging.WARNING

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(stream_level)  # Set level based on verbose flag
    logger.addHandler(stream_handler)


def log_iqtree_run_and_satute_info(
    input_args: List[Namespace],
    substitution_model: SubstitutionModel,
    active_directory: Path,
    rate_category: str,
    msa_file: Path,
    logger: logging.Logger,
) -> None:
    """
    Logs information about the initial IQ-TREE run and tests being performed.
    Args:
        iq_arguments_dict (dict): Dictionary containing IQ-TREE argument configurations.
        substitution_model: The substitution model used in the analysis.
        rate_category: The category of rates being considered.
        msa_file (Path): Path to the MSA file being used.
        multiplicity: Multiplicity value from spectral decomposition.
    """
    considered_rate_category_text = (
        f"\n\tConsidered rate category: {rate_category}"
        if substitution_model.number_rates > 1
        else ""
    )

    considered_edge = (
        f"Run test for saturation for the branch{input_args.edge}"
        if input_args.edge
        else "Run test for saturation for each branch"
    )

    logger.info(
        f"""
        Running tests and initial IQ-Tree with configurations:
        Significance Level Used for Test: {input_args.alpha}
        Run SatuTe on: {msa_file.resolve()}        
        Number of rate categories: {substitution_model.number_rates} {considered_rate_category_text}
        SatuTe will be applied on {considered_edge}        
        Results will be written to the directory: {active_directory.name}                
        """
    )


def log_iqtree_options(
    arguments_dict: Dict[str, List[str]],
    extra_arguments: List[str],
    logger: logging.Logger,
) -> None:
    """
    Logs the used IQ-TREE options.

    Args:
    - arguments_dict (dict): Dictionary containing IQ-TREE arguments.
    - extra_arguments (list): List of additional arguments used in the IQ-TREE run.
    - logger (logging.Logger): Logger instance to use for logging the options.
    """
    logger.info("Used IQ-TREE options:")
    logger.info(" ".join(arguments_dict["arguments"]))
    logger.info(" ".join(extra_arguments))


def log_tested_tree(logger: Logger, tree: Tree, option: str) -> None:
    if "tree" in option:
        logger.info(f"User defined Tree: {tree.write(format=1, format_root_node=True)}")
    else:
        logger.info(
            f"IQ-Tree inferred Tree: {tree.write(format=1, format_root_node=True)}"
        )


def log_rate_and_tree(
    logger: Logger, file_name: str, rate: str, results_set: Dict[str, Any]
) -> None:
    logger.info(f"Writing results for category rates to file: {file_name}")
    if "rescaled_tree" in results_set:
        logger.info(
            f"Tree for rate category {rate}: {results_set['rescaled_tree'].write(format=1, format_root_node=True)}"
        )
