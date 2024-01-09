import re
import pandas as pd
from logging import Logger
from Bio import AlignIO
from pathlib import Path
from ete3 import Tree
from pandas import DataFrame
from typing import Dict, Any
import numpy as np


# New function to format float columns
def format_float_columns(data_frame: DataFrame):
    for col in data_frame.columns:
        if data_frame[col].dtype == float:
            data_frame[col] = data_frame[col].apply(lambda x: round(x, 4))


def construct_file_name(
    msa_file: Path, output_suffix: str, rate: str, alpha: float, edge
):
    """
    Constructs a file name based on the given parameters.
    """
    file_name = f"{msa_file.resolve()}_{rate}_{alpha}.satute"
    if output_suffix:
        file_name = f"{msa_file.resolve()}_{output_suffix}_{rate}_{alpha}.satute"
    if edge:
        file_name += f"_{edge}"
    return file_name


def write_to_csv(data_frame: DataFrame, file_name: str, logger: Logger):
    """
    Writes a DataFrame to a CSV file.
    """
    try:
        csv_file_name = f"{file_name}.csv"
        data_frame.to_csv(csv_file_name)
        logger.info(f"Results saved to CSV file: {csv_file_name}")
    except Exception as e:
        logger.error(f"Error writing CSV file {csv_file_name}: {e}")


def write_results_for_category_rates(
    results: Dict[str, Any],
    to_be_tested_tree: Tree,
    output_suffix: str,
    msa_file: Path,
    alpha: float,
    edge,
    logger: Logger,
):
    """
    Writes the results for category rates to appropriate files.

    Args:
    - results (dict): Dictionary containing the results data.
    - to_be_tested_tree (Tree): Tree object to be tested.
    - output_suffix (str): Suffix for output files.
    - msa_file (Path): Path to the MSA file.
    - alpha (float): Alpha value for calculations.
    - edge: Edge information for tree processing.
    - logger (Logger): Logger for logging events.
    """

    log_original_tree(logger, to_be_tested_tree)

    for rate, results_set in results.items():
        process_rate_category(
            rate, results_set, msa_file, output_suffix, alpha, edge, logger
        )


def write_posterior_probabilities_for_rates(
    results: Dict[str, Any],
    state_frequencies: list,
    output_file: Path,
    output_suffix: str,
    alpha: float,
    edge,
):
    """
    Writes the posterior probabilities to a file.

    Args:
    - results (dict): Dictionary containing the results data.
    - output_file (str): Path to the output file where results will be saved.
    """

    for key, results_set in results.items():
        base_file_name = construct_file_name(
            output_file, output_suffix, key, alpha, edge
        )

        calculate_and_write_posterior_probabilities(
            results_set["partial_likelihoods"], state_frequencies, base_file_name
        )


def process_rate_category(
    rate: str,
    results_set: Dict[str, Any],
    msa_file: Path,
    output_suffix: str,
    alpha: float,
    edge,
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
        results_data_frame = pd.DataFrame(results_set["result_list"])
        format_float_columns(results_data_frame)
        write_results_to_files(results_set, file_name, results_data_frame, logger)
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
    results_set: dict[str, Any],
    file_name: str,
    results_data_frame: pd.DataFrame,
    logger: Logger,
):
    if "rescaled_tree" in results_set:
        write_nexus_file(
            results_set["rescaled_tree"], file_name, results_data_frame, logger
        )
        logger.info(
            f"Saturation Test Results for rate category mapped to tree in Nexus file: {file_name}.nex"
        )

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

        # Write the Nexus file
        with open(file_name, "w") as nexus_file:
            nexus_file.write("#NEXUS\n")
            nexus_file.write("BEGIN TREES;\n")
            nexus_file.write(f"Tree tree1 = {mapped_newick_string}\n")
            nexus_file.write("END TREES;\n")
        logger.info(f"Nexus file written successfully: {file_name}")
    except Exception as e:
        logger.error(f"Error writing Nexus file: {e}")


def write_alignment_and_indices(
    per_rate_category_alignment: dict[str, AlignIO.MultipleSeqAlignment],
    categorized_sites: dict,
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


import pandas as pd


def calculate_and_write_posterior_probabilities(
    partial_likelihood_per_site_storage: dict, state_frequencies: list, output_file: str
):
    """
    Calculate and write posterior probabilities for left and right likelihoods for each edge.

    Args:
    - partial_likelihood_per_site_storage (dict): Storage containing left and right likelihoods for each edge.
    - state_frequencies (list or array): State frequencies used in probability calculations.
    - output_file (str): Path to the output file where results will be saved.
    """

    all_posterior_probabilities = pd.DataFrame()

    for edge, likelihoods in partial_likelihood_per_site_storage.items():
        left_partial_likelihood = pd.DataFrame(likelihoods["left"]["likelihoods"])
        right_partial_likelihood = pd.DataFrame(likelihoods["right"]["likelihoods"])
        
                # Calculate left and right posterior probabilities
        left_posterior_probabilities = calculate_posterior_probabilities_subtree_df(
            4, state_frequencies, left_partial_likelihood
        )
        right_posterior_probabilities = calculate_posterior_probabilities_subtree_df(
            4, state_frequencies, right_partial_likelihood
        )
        
        # left_posterior_probabilities["Site"] = left_partial_likelihood["Site"]
        left_posterior_probabilities["Node"] = left_partial_likelihood["Node"]
        # left_posterior_probabilities["Edge"] = edge
        right_posterior_probabilities["Site"] = right_partial_likelihood["Site"]
        right_posterior_probabilities["Node"] = right_partial_likelihood["Node"]
        right_posterior_probabilities["Edge"] = edge        

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


def calculate_posterior_probabilities_subtree_df(
    dimension: int, state_frequencies: list, partial_likelihood_df: pd.DataFrame
):
    diag = np.diag(list(state_frequencies))

    # Selecting the relevant columns for likelihoods
    likelihood_cols = partial_likelihood_df.iloc[:, 3 : (3 + dimension)]

    # Calculate the site likelihood for each site (row)
    site_likelihoods = likelihood_cols @ diag
    site_likelihoods_sum = site_likelihoods.sum(axis=1)

    # Calculate the posterior probabilities for each site
    posterior_probabilities = site_likelihoods.divide(site_likelihoods_sum, axis=0)

    return posterior_probabilities
