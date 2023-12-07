import pandas as pd

def write_results_for_category_rates(results, input_args, map_values_to_newick, logger):
    """
    Writes the results for category rates to appropriate files.

    Args:
    - results (dict): Dictionary containing the results data.
    - input_args (object): Object containing input arguments.
    - map_values_to_newick (function): Function to map values from DataFrame to the Newick string.
    """
    if not isinstance(results, dict):
        logger.error("Invalid input: results should be a dictionary.")
        return
    for key, results_set in results.items():
        if not isinstance(results_set, dict) or "result_list" not in results_set:
            logger.error(f"Invalid format for results_set under key '{key}'.")
            continue

        try:
            file_name = f"{input_args.msa.resolve()}_{key}_{input_args.alpha}.satute"

            # Append the edge information to the file name if provided
            if hasattr(input_args, "edge") and input_args.edge:
                file_name = f"{file_name}_{input_args.edge}"

            results_data_frame = pd.DataFrame(
                results_set["result_list"]
            )

            if "rescaled_tree" in results_set:
                tree_file_name = f"{file_name}.nex"
                newick_string = results_set["rescaled_tree"].write(format=1)
                write_nexus_file(
                    newick_string,
                    tree_file_name,
                    results_data_frame,
                    map_values_to_newick,
                    logger,
                )

            csv_file_name = f"{file_name}.csv"
            try:
                results_data_frame.to_csv(csv_file_name, engine="python")
            except Exception as e:
                logger.error(f"Error writing CSV file {csv_file_name}: {e}")

        except Exception as e:
            logger.error(f"Error processing results for key '{key}': {e}")

    logger.info("Finished writing results for category rates to files")


def write_results_for_single_rate(
    results, input_args, to_be_tested_tree, map_values_to_newick, logger
):
    """
    Writes the results for single rate to appropriate files.

    Args:
    - results (dict): Dictionary containing the results data.
    - input_args (object): Object containing input arguments.
    - to_be_tested_tree (Tree): The tree object containing the data to be written.
    - map_values_to_newick (function): Function to map values from DataFrame to the Newick string.
    """
    for key, results_set in results.items():
        file_name_base = f"{input_args.msa.resolve()}_{input_args.alpha}.satute"

        # Append the edge information to the file name if provided
        if hasattr(input_args, "edge") and input_args.edge:
            file_name_base += f"_{input_args.edge}"

        results_data_frame = pd.DataFrame(results_set)

        # Writing the .csv file
        csv_file_name = f"{file_name_base}.csv"
        try:
            results_data_frame.to_csv(csv_file_name, engine="python")
            logger.info(f"Finished writing CSV file: {csv_file_name}")
        except Exception as e:
            logger.error(f"Error writing CSV file {csv_file_name}: {e}")

        # Writing the Nexus file
        if to_be_tested_tree:
            tree_file_name = file_name_base
            newick_string = to_be_tested_tree.write(format=1)
            write_nexus_file(
                newick_string, tree_file_name, results_data_frame, map_values_to_newick, logger
            )


def write_nexus_file(
    newick_string, file_name, results_data_frame, map_values_to_newick, logger
):
    """
    Writes the given Newick string as a Nexus file after mapping values.

    Args:
    - newick_string (str): Newick formatted string representing the tree.
    - file_name (str): Name of the file to write.
    - results_data_frame (DataFrame): DataFrame with the data to map onto the Newick string.
    - map_values_to_newick (function): Function to map values from DataFrame to the Newick string.
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
