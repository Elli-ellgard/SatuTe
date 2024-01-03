import re
import pandas as pd
from logging import Logger
from Bio import AlignIO
from pathlib import Path
from ete3 import Tree
from pandas import DataFrame
from typing import Dict, Any


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


def write_results_for_single_rate(
    results: dict,
    output_suffix: str,
    msa_file: Path,
    to_be_tested_tree: Tree,
    alpha: float,
    edge: tuple,
    logger: Logger,
):
    """
    Writes the results for single rate to appropriate files.
    Args:
    - results (dict): Dictionary containing the results data.
    - input_args (object): Object containing input arguments.
    - to_be_tested_tree (Tree): The tree object containing the data to be written.
    """
    for key, results_set in results.items():
        file_name_base = construct_file_name(msa_file, output_suffix, key, alpha, edge)
        results_data_frame = pd.DataFrame(results_set)
        format_float_columns(results_data_frame)
        write_to_csv(results_data_frame, file_name_base, logger)
        write_nexus_file(to_be_tested_tree, file_name_base, results_data_frame, logger)


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
    categorized_sites,
    input_args,
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
            file_path = f"{input_args.msa.resolve()}.{rate}.phy.rate.indices"
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
