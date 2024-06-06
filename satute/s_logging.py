import logging
import numpy as np
from pathlib import Path
from typing import List
from satute.repository import SubstitutionModel

def format_matrix(matrix, precision: int = 4):
    """Format a matrix for pretty printing."""
    formatted_matrix = "\n".join(
        ["\t".join([f"{item:.{precision}f}" for item in row]) for row in matrix]
    )
    return formatted_matrix

def format_array(array, precision=4):
    """Format a 1D array for pretty printing."""
    formatted_array = "\t".join([f"{item:.{precision}f}" for item in array])
    return formatted_array


def construct_log_file_name(msa_file: Path, input_args):
    log_file = f"{msa_file.resolve()}_{input_args.alpha}.satute.log"
    if input_args.output_suffix:
        log_file = f"{msa_file.resolve()}_{input_args.alpha}_{input_args.output_suffix}.satute.log"
    return log_file

def setup_logging_configuration(logger, input_args, msa_file: Path):
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
    file_format = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    file_handler.setFormatter(file_format)
    file_handler.setLevel(logging.DEBUG)  # Always log everything in file
    logger.addHandler(file_handler)
        
    
    # Set the default logging level
    if input_args.verbose:
        stream_level = logging.DEBUG
    elif input_args.quiet:
        stream_level = logging.CRITICAL
    else:
        stream_level = logging.WARNING
        
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(file_format)
    stream_handler.setLevel(stream_level)  # Set level based on verbose flag
    
    logger.addHandler(stream_handler)
    
    
def log_iq_tree_run_and_satute_info(
    input_args,
    substitution_model: SubstitutionModel,
    active_directory,
    rate_category: str,
    msa_file: Path,
    multiplicity: int,
    logger: logging.Logger,
):
    """
    Logs information about the initial IQ-TREE run and tests being performed.
    Args:
        iq_arguments_dict (dict): Dictionary containing IQ-TREE argument configurations.
        substitution_model: The substitution model used in the analysis.
        rate_category: The category of rates being considered.
        msa_file (Path): Path to the MSA file being used.
        multiplicity: Multiplicity value from spectral decomposition.
    """
    logger.info(
        f"""
        Running tests and initial IQ-Tree with configurations:
        Model: {input_args.model}
        Alpha: {input_args.alpha}
        Running Saturation Test on file: {msa_file.resolve()}
        Number of rate categories: {substitution_model.number_rates}
        Considered rate category: {rate_category}
        Multiplicity: {multiplicity}
        Run test for saturation for each branch and category with {substitution_model.number_rates} rate categories
        Results will be written to the directory: {active_directory.name}
        """
    )
    
    
def log_substitution_model_info(
        logger: logging.Logger,
        input_args,
        substitution_model: SubstitutionModel,
        multiplicity: int,
        eigenvectors: List[np.array],
        eigenvalue: float,
    ):
        """
        Logs information about substitution model and its spectral decomposition

        Args:
            substitution_model: The substitution model used in the analysis.
            multiplicity: Multiplicity value from spectral decomposition.
            eigenvectors: eigenvector corresponding to eigenvalue from spectral decomposition
            eigenvalue:  dominant non-zero eigenvalue from spectral decomposition
        """
        # Formatting the rate matrix for logging
        rate_matrix_str = format_matrix(substitution_model.rate_matrix, precision=4)
        # Formatting the state frequencies for logging
        state_frequencies_str = format_array(
            np.array(list(substitution_model.state_frequencies)), precision=4
        )

        # Logging the formatted rate matrix and state frequencies
        logger.info(
            f"Substitution Model:\n\n"
            f"Model: {input_args.model}\n"
            f"Rate Matrix Q:\n{rate_matrix_str}\n"
            f"State Frequencies:\n{state_frequencies_str}\n"
        )

        eigenvector_str = ""
        for eigenvector in eigenvectors:
            eigenvector_str += f"\n{format_array(list(eigenvector))}"

        logger.info(
            f"Spectral Decomposition:\n\n"
            f"Eigenvalue: {eigenvalue}\n"
            f"Multiplicity: {multiplicity}\n"
            f"Eigenvectors: {eigenvector_str}\n"
        )
