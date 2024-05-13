# -*- coding: utf-8 -*-

import re
import pandas as pd
from logging import Logger
from Bio import AlignIO
from pathlib import Path
from ete3 import Tree
from pandas import DataFrame
from typing import Dict, Any, List
import numpy as np
from amino_acid_models import AMINO_ACIDS
from satute_result import TestStatisticComponentsContainer
from satute_statistic_posterior_distribution_components import calculate_posterior_probabilities_subtree_df


# New function to format float columns
def format_float_columns(data_frame: DataFrame):
    for col in data_frame.columns:
        if col != "branch_length":
            if data_frame[col].dtype == float:
                data_frame[col] = data_frame[col].apply(lambda x: round(x, 4))




def construct_file_name(
    msa_file: Path, output_suffix: str, rate: str, alpha: float, edge: str
) -> str:
    """
    Constructs a file name for saving analysis results, incorporating multiple parameters to generate a meaningful and unique file name.

    The function constructs the file name by starting with the path of the MSA (Multiple Sequence Alignment) file, followed by the specified evolutionary rate, and the alpha value. These components are concatenated using underscores. If an output suffix is provided, it is inserted between the MSA file path and the rate. If edge information is provided, it is appended at the end. The resulting file name ends with a ".satute" extension.
    Args:
        msa_file (Path): The path to the MSA file, which forms the base of the file name.
        output_suffix (str): An optional string to be included in the file name for additional context or differentiation.
        rate (str): A string representing the evolutionary rate category, included in the file name to indicate the analysis's focus.
        alpha (float): A float value representing the alpha parameter used in the analysis, included in the file name for specificity.
        edge (str): An optional string representing edge information in the phylogenetic analysis, appended to the file name if provided.

    Returns:
        str: The constructed file name, incorporating the provided parameters for clear identification of the analysis results.
    """
    file_name = f"{msa_file.resolve()}_{rate}_{alpha}"
    if output_suffix:
        file_name = f"{msa_file.resolve()}_{output_suffix}_{rate}_{alpha}"
    if edge:
        file_name += f"_{edge}"
    file_name += ".satute"
    return file_name


def write_components(
    components: TestStatisticComponentsContainer,
    msa_file: Path,
    output_suffix: str,
    rate: str,
    alpha: float,
    edge: str,
    site_indices: List[int],
):
    """
    Writes the components' DataFrame to a CSV file, including the same site indices for each identifier.

    Args:
        components (TestStatisticComponentsContainer): The container of test statistic components.
        msa_file (Path): Path to the MSA file, used for constructing the output file name.
        output_suffix (str): Suffix for the output file name.
        rate (str): Rate parameter for the file name.
        alpha (float): Alpha parameter for the file name.
        edge (str): Edge parameter for the file name.
        site_indices (List[int]): List of site indices to be included in the DataFrame for each identifier.
    """
    file_name = construct_file_name(msa_file, output_suffix, rate, alpha, edge)
    components_frame = components.to_dataframe()

    components_frame["rate"] = rate

    # Repeat the site_indices for each row in the DataFrame based on the identifier
    # Assuming every row/component should have an associated site index
    expanded_site_indices = []
    num_rows_per_identifier = len(site_indices)

    for identifier in components_frame["Edge"].unique():
        expanded_site_indices.extend(site_indices)

    # Check to ensure the expanded list matches the DataFrame's length
    # if len(expanded_site_indices) != len(components_frame):
    #    raise ValueError("Expanded site indices do not match the DataFrame length.")

    # Add the expanded site indices to the DataFrame
    components_frame["site"] = expanded_site_indices

    # Write the DataFrame to a CSV file
    components_frame.to_csv(f"{file_name}.components.csv", index=False)


def write_to_csv(data_frame: pd.DataFrame, file_name: str, logger: Logger) -> None:
    """
    Writes the content of a pandas DataFrame to a CSV file. This function constructs the CSV file name
    by appending ".csv" to the provided file_name parameter. It attempts to save the DataFrame into this
    CSV file and logs the outcome using the provided logger object.

    Extended error handling includes catching specific exceptions related to file writing operations
    and general exceptions, allowing for more granular logging and troubleshooting.

    Args:
        data_frame (pd.DataFrame): The DataFrame to be written to a CSV file.
        file_name (str): The base name for the output CSV file, without the ".csv" extension.
        logger (Logger): A logging.Logger object used for logging messages regarding the file writing operation.

    Returns:
        None

    Raises:
        IOError: If an I/O error occurs during file writing.
        Exception: For catching other unexpected exceptions that may occur during the process.
    """
    try:
        csv_file_name = f"{file_name}.csv"
        data_frame.to_csv(
            csv_file_name, index=False
        )  # Consider adding index=False to avoid writing row indices
        logger.info(f"Results saved to CSV file: {csv_file_name}")
    except PermissionError as e:
        logger.error(f"Permission denied when writing CSV file {csv_file_name}: {e}")
        raise IOError(f"Permission denied when writing CSV file {csv_file_name}: {e}")
    except IOError as e:
        logger.error(f"I/O error when writing CSV file {csv_file_name}: {e}")
        raise
    except Exception as e:
        logger.error(
            f"Unexpected error occurred when writing CSV file {csv_file_name}: {e}"
        )
        raise


def write_results_for_category_rates(
    results: Dict[str, Any],
    to_be_tested_tree: Tree,
    output_suffix: str,
    msa_file: Path,
    alpha: float,
    edge: str,
    categorized_sites: List[int],
    logger: Logger,
) -> None:
    """
    Iterates over a dictionary of results, organized by rate categories, and writes these results to files.
    Each rate category's results are processed and saved using a dedicated function. This process includes
    logging the original tree structure and handling results for each category rate by delegating to another
    function designed for processing and saving those results.

    This function serves as an orchestrator for writing the analysis results of different evolutionary rate
    categories into separate files, each named according to the provided parameters to reflect the analysis's specifics.

    Args:
        results (Dict[str, Any]): A dictionary where keys are rate categories (as strings) and values are the results
        data (of any type) associated with those categories.
        to_be_tested_tree (Tree): An ETE tree object that represents the phylogenetic tree to be tested or analyzed.
        output_suffix (str): A suffix to be appended to the output files' names, providing context or distinguishing
        between different analysis runs.
        msa_file (Path): The file path to the multiple sequence alignment (MSA) file that was analyzed. This is used
        as part of the basis for naming output files.
        alpha (float): A parameter value used in the analysis, included in the output file names for reference.
        edge (str): A string representing specific edge information in the phylogenetic tree, also included in the
        output file names.
        logger (Logger): A logging.Logger object used for logging informational messages and errors during the process.

    Returns:
        None

    Raises:
        Exception: If an error occurs during the logging of the original tree or the processing of category rates,
        exceptions are caught and logged. Specific handling or re-raising of exceptions should be implemented in the
        functions `log_original_tree` and `process_rate_category`.
    """
    try:
        log_original_tree(logger, to_be_tested_tree)
    except Exception as e:
        logger.error(f"Error logging original tree: {e}")
        # Consider whether to continue execution or raise an exception based on your application's needs.

    for rate, results_set in results.items():
        try:

            process_rate_category(
                rate,
                results_set,
                msa_file,
                output_suffix,
                alpha,
                edge,
                categorized_sites[rate],
                logger,
            )

        except Exception as e:
            logger.error(f"Error processing results for rate category '{rate}': {e}")
            # Decide on continuation or halting based on the severity of the error and your application's requirements.


def write_posterior_probabilities_for_rates(
    results: Dict[str, Any],
    state_frequencies: List[float],
    output_file: Path,
    output_suffix: str,
    alpha: float,
    edge: str,
    categorized_sites: Dict[str, List[int]],
) -> None:
    """
    Writes the posterior probabilities for different evolutionary rates to files. Each rate category in the results
    dictionary will have its posterior probabilities calculated and written to a separate file. The file names are
    constructed using the output_file path, output_suffix, the rate category, the alpha parameter, and the edge
    information, ensuring unique and descriptive file names for each rate category.

    Args:
        results (Dict[str, Any]): A dictionary where keys are rate categories (as strings) and values are dictionaries
            containing the results data for that rate, including "partial_likelihoods" which are necessary for
            calculating posterior probabilities.
        state_frequencies (List[float]): A list of state frequencies used in the calculation of posterior probabilities.
            Each element in the list represents the frequency of a particular state (e.g., nucleotide or amino acid)
            in the analyzed sequences.
        output_file (Path): The base path for output files. The actual output file names will be constructed by
            appending additional information to this base path.
        output_suffix (str): A string to be appended to the base file name to help differentiate files from different
            analyses or runs.
        alpha (float): The alpha parameter used in the analysis, included in the file name for specificity.
        edge (str): A string representing edge information in the phylogenetic tree, included in the file name if provided.
        categorized_sites (Dict[str, List[int]]): A dictionary mapping rate categories to lists of site indices. Each
            list contains integers that represent the positions (sites) in the sequence alignment associated with that rate.

    Returns:
        None
    """
    for rate, results_set in results.items():
        # Construct the base file name for this rate category
        base_file_name = construct_file_name(
            output_file, output_suffix, rate, alpha, edge
        )

        # Calculate and write posterior probabilities for this rate category
        calculate_and_write_posterior_probabilities(
            results_set["partial_likelihoods"],
            state_frequencies,
            base_file_name,
            categorized_sites[rate],
        )


def write_posterior_probabilities_single_rate(
    results: Dict[str, Any],
    state_frequencies: List[float],
    output_file: Path,
    output_suffix: str,
    alpha: float,
    edge: str,
    number_sites: int,
):
    """
    Writes the posterior probabilities to a file.

    Args:
    - results (dict): Dictionary containing the results data.
    - output_file (str): Path to the output file where results will be saved.
    """

    base_file_name = construct_file_name(
        output_file, output_suffix, "single_rate", alpha, edge
    )

    categorized_site = [i for i in range(0, number_sites + 1, 1)]

    calculate_and_write_posterior_probabilities(
        results["single_rate"]["partial_likelihoods"],
        state_frequencies,
        base_file_name,
        categorized_site,
    )


def process_rate_category(
    rate: str,
    results_set: Dict[str, Any],
    msa_file: Path,
    output_suffix: str,
    alpha: float,
    edge: str,
    for_categorized_rate_sites: List[int],
    logger: Logger,
):
    """
    Processes and writes results for a single rate category.

    Args:
    - rate (str): The rate category.
    - results_set (Dict[str, Any]): Results data for the rate category.
    - msa_file (Path): Path to the MSA file.
    - output_suffix (str): Suffix for output files.
    - alpha (float): Alpha value for calculations.
    - edge: Edge information for tree processing.
    - logger (Logger): Logger for logging events.
    """
    try:
        file_name = construct_file_name(msa_file, output_suffix, rate, alpha, edge)

        log_rate_info(logger, file_name, rate, results_set)

        results_data_frame = pd.DataFrame(results_set["result_list"].to_dataframe())
        # add sequence length for considered rate category and rate category
        results_data_frame["number_of_sites"] = len(for_categorized_rate_sites)
        results_data_frame["rate_category"] = rate

        format_float_columns(results_data_frame)

        write_results_to_files(
            results_set,
            file_name,
            results_data_frame,
            logger,
        )

    except Exception as e:
        logger.error(f"Error processing results for key '{rate}': {e}")


def log_original_tree(logger: Logger, tree: Tree):
    logger.info(f"Original Tree: {tree.write(format=1, format_root_node=True)}")


def log_rate_info(
    logger: Logger, file_name: str, rate: str, results_set: Dict[str, Any]
):
    logger.info(f"Writing results for category rates to file: {file_name}")
    if "rescaled_tree" in results_set:
        logger.info(
            f"Tree for rate category {rate}: {results_set['rescaled_tree'].write(format=1, format_root_node=True)}"
        )


def write_results_to_files(
    results_set: Dict[str, Any],
    file_name: str,
    results_data_frame: DataFrame,
    logger: Logger,
):
    if "rescaled_tree" in results_set:
        write_nexus_file(
            results_set["rescaled_tree"], file_name, results_data_frame, logger
        )
        logger.info(
            f"Saturation Test Results for rate category mapped to tree in Nexus file: {file_name}.nex"
        )
    # write_to_components(results_set['components'], file_name, [])
    write_to_csv(results_data_frame, file_name, logger)


def get_target_node(edge):
    return edge.split(",")[0].replace("(", "").strip()


def map_values_to_newick(newick: str, results_per_branch: DataFrame) -> str:
    """
    Maps values from the DataFrame onto the Newick string by dynamically modifying the metadata for each node.

    Args:
    - newick (str): The Newick string to be modified.
    - results_per_branch (DataFrame): DataFrame containing the data to map onto the Newick string.

    Returns:
    - str: Modified Newick string with updated metadata.
    """
    for index, row in results_per_branch.iterrows():
        newick = update_node_metadata(newick, row, results_per_branch.columns)
    return newick


def update_node_metadata(newick: str, row: DataFrame, columns: list) -> str:
    """
    Updates the metadata for a single node in the Newick string based on the DataFrame row.

    Args:
    - newick (str): Current Newick string.
    - row (DataFrame): A single row from the DataFrame.
    - columns (list): List of columns in the DataFrame.

    Returns:
    - str: Newick string with updated metadata for the node.
    """
    target_node = get_target_node(row["edge"])
    escaped_target_node = re.escape(target_node)
    meta_data = create_meta_data_string(row, columns)

    newick = insert_metadata_into_newick(newick, escaped_target_node, meta_data)
    return newick


def create_meta_data_string(row: DataFrame, columns: list) -> str:
    """
    Creates a metadata string from a DataFrame row.

    Args:
    - row (DataFrame): A single row from the DataFrame.
    - columns (list): List of columns in the DataFrame.

    Returns:
    - str: Metadata string.
    """
    meta_data_parts = [f"{col}={row[col]}" for col in columns if col != "edge"]
    return "&" + ",".join(meta_data_parts) if meta_data_parts else ""


def insert_metadata_into_newick(newick: str, target_node: str, meta_data: str) -> str:
    """
    Inserts metadata into the Newick string for a specific node.

    Args:
    - newick (str): The Newick string.
    - target_node (str): Target node for which metadata needs to be updated.
    - meta_data (str): Metadata string to be inserted.

    Returns:
    - str: Newick string with updated metadata.
    """
    pattern_with_brackets = re.compile(
        rf"({target_node}:\d+(\.\d+)?(e-?\d+)?)\[([^\]]+)\]"
    )
    pattern_without_brackets = re.compile(rf"({target_node}:\d+(\.\d+)?(e-?\d+)?)")

    if pattern_with_brackets.search(newick):
        return pattern_with_brackets.sub(rf"\1[\4,{meta_data}]", newick)
    else:
        return pattern_without_brackets.sub(rf"\1[{meta_data}]", newick)


def write_nexus_file(
    tree: Tree, file_name: str, results_data_frame: DataFrame, logger: Logger
):
    """
    Writes the given Newick string as a Nexus file after mapping values.
    Args:
    - newick_string (str): Newick formatted string representing the tree.
    - file_name (str): Name of the file to write.
    - results_data_frame (DataFrame): DataFrame with the data to map onto the Newick string.
    """

    try:
        newick_string = tree.write(format=1, format_root_node=True)
        # Validate input
        if not isinstance(newick_string, str):
            raise ValueError("Newick string must be a string.")
        if not file_name.endswith(".nex"):
            file_name += ".nex"

        # Map values to Newick string
        mapped_newick_string = map_values_to_newick(newick_string, results_data_frame)

        # Extract taxa names from the tree
        taxa_names = [leaf.name for leaf in tree.iter_leaves()]

        # Write the Nexus file
        with open(file_name, "w") as nexus_file:
            nexus_file.write("#NEXUS\n")
            nexus_file.write("BEGIN TAXA;\n")
            nexus_file.write(f"    DIMENSIONS NTAX={len(taxa_names)};\n")
            nexus_file.write("     TAXLABELS\n")
            for taxon in taxa_names:
                nexus_file.write(f"        {taxon}\n")
            nexus_file.write("    ;\n")
            nexus_file.write("END;\n\n")
            nexus_file.write("BEGIN TREES;\n")
            nexus_file.write(f"Tree tree1 = {mapped_newick_string}\n")
            nexus_file.write("END;\n")
        logger.info(f"Nexus file written successfully: {file_name}")
    except Exception as e:
        logger.error(f"Error writing Nexus file: {e}")


def write_alignment_and_indices(
    per_rate_category_alignment: Dict[str, AlignIO.MultipleSeqAlignment],
    categorized_sites: List[int],
    msa_file: Path,
):
    """
    Writes MultipleSeqAlignment objects and indices to files.
    Parameters:
    - per_rate_category_alignment (dict): A dictionary mapping rates to MultipleSeqAlignment objects.
    - categorized_sites (dict): A dictionary mapping rates to lists of integers (indices).
    Returns:
    - None
    """
    try:
        for rate in per_rate_category_alignment.keys():
            file_path = f"{msa_file.resolve()}.{rate}.phy.rate.indices"
            with open(file_path, "w") as file:
                if rate in per_rate_category_alignment and rate in categorized_sites:
                    if per_rate_category_alignment[rate].get_alignment_length() == 0:
                        continue
                    # Convert MultipleSeqAlignment to string in FASTA format
                    AlignIO.write(per_rate_category_alignment[rate], file, "phylip")
                    file.write(",".join([str(i) for i in categorized_sites[rate]]))
    except TypeError as e:
        print(f"TypeError: {e}")
    except IOError as e:
        print(f"IOError: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


def calculate_and_write_posterior_probabilities(
    partial_likelihood_per_site_storage: Dict,
    state_frequencies: List[float],
    output_file: str,
    categorized_sites: List[int],
):
    """
    Calculate and write posterior probabilities for left and right likelihoods for each edge.

    Args:
    - partial_likelihood_per_site_storage (dict): Storage containing left and right likelihoods for each edge.
    - state_frequencies (list or array): State frequencies used in probability calculations.
    - output_file (str): Path to the output file where results will be saved.
    """

    all_posterior_probabilities = pd.DataFrame()

    states = []  # Extend this list if dealing with amino acids
    dimension = len(state_frequencies)
    # Define state labels based on the dimension
    if len(state_frequencies) == 4:
        states = ["A", "C", "G", "T"]
    elif len(state_frequencies) == 20:
        states = AMINO_ACIDS
    else:
        raise ValueError(
            "Unsupported state frequency dimension. Expected 4 (nucleotides) or 20 (amino acids)."
        )

    for edge, likelihoods in partial_likelihood_per_site_storage.items():

        left_partial_likelihood = pd.DataFrame(likelihoods["left"]["likelihoods"])

        right_partial_likelihood = pd.DataFrame(likelihoods["right"]["likelihoods"])

        # Calculate left and right posterior probabilities
        left_posterior_probabilities = calculate_posterior_probabilities_subtree_df(
            dimension, state_frequencies, left_partial_likelihood
        )

        right_posterior_probabilities = calculate_posterior_probabilities_subtree_df(
            dimension, state_frequencies, right_partial_likelihood
        )

        left_posterior_probabilities["Node"] = left_partial_likelihood["Node"]

        right_posterior_probabilities["Site"] = [
            categorized_sites[i] for i in right_partial_likelihood["Site"]
        ]

        right_posterior_probabilities["Node"] = right_partial_likelihood["Node"]

        right_posterior_probabilities["Edge"] = edge

        # Map numerical indices to state names and add side suffix
        state_columns_left = {i: "p" + state for i, state in enumerate(states)}
        left_posterior_probabilities.rename(columns=state_columns_left, inplace=True)

        # Map numerical indices to state names and add side suffix
        state_columns_right = {i: "p" + state for i, state in enumerate(states)}
        right_posterior_probabilities.rename(columns=state_columns_right, inplace=True)

        # Rename columns with suffixes _left and _right
        left_posterior_probabilities = left_posterior_probabilities.add_suffix("_left")

        right_posterior_probabilities = right_posterior_probabilities.add_suffix(
            "_right"
        )

        # Ensure 'Site', 'Node', and 'Edge' columns are not suffixed
        for col in ["Site", "Edge"]:

            if col + "_left" in left_posterior_probabilities.columns:

                left_posterior_probabilities[col] = left_posterior_probabilities[
                    col + "_left"
                ]

                left_posterior_probabilities.drop(col + "_left", axis=1, inplace=True)

            if col + "_right" in right_posterior_probabilities.columns:

                right_posterior_probabilities[col] = right_posterior_probabilities[
                    col + "_right"
                ]

                right_posterior_probabilities.drop(col + "_right", axis=1, inplace=True)

        # Concatenate left and right posterior probabilities horizontally
        edge_posterior_probabilities = pd.concat(
            [left_posterior_probabilities, right_posterior_probabilities], axis=1
        )

        # Add the concatenated data to the master DataFrame
        all_posterior_probabilities = pd.concat(
            [all_posterior_probabilities, edge_posterior_probabilities],
            ignore_index=True,
        )

    # Write to the output file
    all_posterior_probabilities.to_csv(f"{output_file}.asr.csv", index=False)



