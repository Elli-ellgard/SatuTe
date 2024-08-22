import logging
from pathlib import Path
from ete3 import Tree
from argparse import Namespace
from logging import Logger
from typing import List, Dict, Any
import os


def log_consider_iqtree_message(logger: Logger):
    logger.info("Running IQ-TREE")
    logger.warning(
        "Please consider for the analysis that IQ-Tree will be running with any advanced options."
    )
    logger.warning(
        "If specific options are required for the analysis, please run IQ-Tree separately."
    )


def construct_log_file_name(msa_file: Path, input_args: List[Namespace]) -> str:
    # Start with the base file name and mandatory alpha value
    parts = [msa_file.resolve(), input_args.alpha]

    # Append output suffix if it exists
    if input_args.output_suffix:
        parts.insert(1, input_args.output_suffix)

    # Append edge if it exists
    if input_args.edge:
        parts.append(input_args.edge)

    # Join all parts together with underscores and add the file extension
    log_file = f"{'_'.join(map(str, parts))}.satute.log"
    return log_file


def close_log_handlers(logger):
    for handler in logger.handlers:
        handler.close()
        logger.removeHandler(handler)


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
    close_log_handlers(logger)  # Close and remove all existing handlers

    log_file = construct_log_file_name(msa_file, input_args)

    if os.path.exists(log_file):
        os.remove(log_file)

    file_handler = logging.FileHandler(log_file, mode="w")
    file_format = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    file_handler.setFormatter(file_format)
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
    active_directory: Path,
    msa_file: Path,
    iq_tree_arguments: Dict,
    logger: logging.Logger,
) -> None:
    logs = [
        f"Run SatuTe on: {msa_file.resolve()}",
        f"Tree and parameters are read from: {msa_file.resolve()}.iqtree",
        f"Results will be written to the directory: {active_directory.resolve()}",
        f"Running tests and initial IQ-Tree Run with configurations: {iq_tree_arguments['option']}",
    ]

    for log in logs:
        logger.info(log)


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
        logger.info(f"User defined tree: {tree.write(format=1, format_root_node=True)}")
    else:
        logger.info(
            f"IQ-Tree inferred tree: {tree.write(format=1, format_root_node=True)}"
        )


def log_rate_and_tree(
    logger: Logger, file_name: str, rate: str, results_set: Dict[str, Any]
) -> None:
    logger.info(f"Writing results for category rates to file: {file_name}")
    if "rescaled_tree" in results_set:
        logger.info(
            f"Tree for rate category {rate}: {results_set['rescaled_tree'].write(format=1, format_root_node=True)}"
        )
