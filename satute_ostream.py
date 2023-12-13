import pandas as pd
from logging import Logger
from Bio import AlignIO
import re


# New function to format float columns
def format_float_columns(data_frame):
    for col in data_frame.columns:
        if data_frame[col].dtype == float:
            data_frame[col] = data_frame[col].apply(lambda x: round(x, 4))


def write_results_for_category_rates(results, msa_file, alpha, edge, logger: Logger):
    """
    Writes the results for category rates to appropriate files.

    Args:
    - results (dict): Dictionary containing the results data.
    - input_args (object): Object containing input arguments.
    """
    if not isinstance(results, dict):
        logger.error("Invalid input: results should be a dictionary.")
        return
    for key, results_set in results.items():
        if not isinstance(results_set, dict) or "result_list" not in results_set:
            logger.error(f"Invalid format for results_set under key '{key}'.")
            continue

        try:
            file_name = f"{msa_file}_{key}_{alpha}.satute"

            # Append the edge information to the file name if provided
            if edge:
                file_name = f"{file_name}_{edge}"
                
            results_data_frame = pd.DataFrame(results_set["result_list"])
            format_float_columns(results_data_frame)

            if "rescaled_tree" in results_set:
                tree_file_name = f"{file_name}.nex"
                newick_string = results_set["rescaled_tree"].write(format=1, format_root_node=True)
                write_nexus_file(
                    newick_string,
                    tree_file_name,
                    results_data_frame,
                    logger,
                )

            csv_file_name = f"{file_name}.csv"
            try:
                results_data_frame.to_csv(csv_file_name)
            except Exception as e:
                logger.error(f"Error writing CSV file {csv_file_name}: {e}")

        except Exception as e:
            logger.error(f"Error processing results for key '{key}': {e}")

    logger.info("Finished writing results for category rates to files")


def write_results_for_single_rate(
    results, msa_file, to_be_tested_tree, alpha, edge, logger
):
    """
    Writes the results for single rate to appropriate files.

    Args:
    - results (dict): Dictionary containing the results data.
    - input_args (object): Object containing input arguments.
    - to_be_tested_tree (Tree): The tree object containing the data to be written.
    """
    for key, results_set in results.items():
        file_name_base = f"{msa_file.resolve()}_{alpha}.satute"

        # Append the edge information to the file name if provided
        if edge:
            file_name_base += f"_{edge}"

        results_data_frame = pd.DataFrame(results_set)
        format_float_columns(results_data_frame)

        # Writing the .csv file
        csv_file_name = f"{file_name_base}.csv"
        try:
            results_data_frame.to_csv(csv_file_name)
            logger.info(f"Finished writing CSV file: {csv_file_name}")
        except Exception as e:
            logger.error(f"Error writing CSV file {csv_file_name}: {e}")

        # Writing the Nexus file
        if to_be_tested_tree:
            tree_file_name = file_name_base
            newick_string = to_be_tested_tree.write(format=1)
            write_nexus_file(
                newick_string,
                tree_file_name,
                results_data_frame,
                logger,
            )


def get_target_node(edge):
    return edge.split(",")[0].replace("(", "").strip()


def map_values_to_newick(newick, df):
    for index, row in df.iterrows():
        target_node = get_target_node(row["edge"])

        # Escape any special characters in target node for regex
        escaped_target_node = re.escape(target_node)

        # Adjusting the metadata as per the requirements
        meta_data = f"&delta={row['delta']},c_s={row['c_s']},p_value={row['p_value']},result_test={row['result_test']}"

        # Check for existing square brackets after the node
        # Use raw string notation for regex patterns
        pattern_with_brackets = re.compile(
            rf"({escaped_target_node}:\d+(\.\d+)?(e-?\d+)?)\[([^\]]+)\]"
        )
        pattern_without_brackets = re.compile(
            rf"({escaped_target_node}:\d+(\.\d+)?(e-?\d+)?)"
        )

        # If square brackets are present, append the metadata inside those brackets
        if pattern_with_brackets.search(newick):
            newick = pattern_with_brackets.sub(rf"\1[\4,{meta_data}]", newick)
        else:
            # If no square brackets, add them and insert the metadata
            newick = pattern_without_brackets.sub(rf"\1[{meta_data}]", newick)
    return newick


def write_nexus_file(newick_string, file_name, results_data_frame, logger):
    """
    Writes the given Newick string as a Nexus file after mapping values.
    Args:
    - newick_string (str): Newick formatted string representing the tree.
    - file_name (str): Name of the file to write.
    - results_data_frame (DataFrame): DataFrame with the data to map onto the Newick string.
    """
    try:
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
    per_rate_category_alignment, categorized_sites, input_args
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
                    # Convert MultipleSeqAlignment to string in FASTA format
                    AlignIO.write(per_rate_category_alignment[rate], file, "phylip")
                    file.write(",".join([str(i) for i in categorized_sites[rate]]))
    except TypeError as e:
        print(f"TypeError: {e}")
    except IOError as e:
        print(f"IOError: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
