from collections import defaultdict
import subprocess
import random
import string
import os
import shutil
from ete3 import Tree
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import shutil


iqtree = "iqtree"


def fetch_files_with_prefix_and_extension(directory, prefix, extension):
    """
    Fetch all files in a directory that start with a specific prefix and have a specific file type (extension).

    Parameters:
    - directory: The path to the directory to search.
    - prefix: The file prefix to match.
    - extension: The file extension to match (e.g., ".txt" for text files).

    Returns:
    - A list of filenames that match the prefix and extension.
    """

    # Ensure the extension starts with a dot
    if not extension.startswith("."):
        extension = "." + extension

    # List all files in the directory
    all_files = os.listdir(directory)

    # Filter files based on the prefix and extension
    matching_files = [
        f for f in all_files if f.startswith(prefix) and f.endswith(extension)
    ]

    return matching_files


import os


def fetch_files_with_prefix_and_extension_absolute(directory, prefix, extension):
    """
    Fetch all files in a directory that start with a specific prefix and have a specific file type (extension).

    Parameters:
      - directory: The path to the directory to search.
      - prefix: The file prefix to match.
      - extension: The file extension to match (e.g., ".txt" for text files).

    Returns:
      - A list of filenames that match the prefix and extension, absolute paths.
    """

    # Ensure the extension starts with a dot
    if not extension.startswith("."):
        extension = "." + extension

    # List all files in the directory
    all_files = os.listdir(directory)

    # Filter files based on the prefix and extension
    matching_files = [
        os.path.join(directory, f)
        for f in all_files
        if f.startswith(prefix) and f.endswith(extension)
    ]

    return matching_files


def execute_command(command):
    subprocess.run(command, shell=True)


def generate_tree(n, file_name="random_generated.tree"):
    cell_names = [f"{random.choice(string.ascii_uppercase)}_cell-{i}" for i in range(n)]

    t = Tree()

    t.populate(
        size=n,
        names_library=cell_names,
        reuse_names=False,
        branch_range=(0.1, 1.0),
        support_range=(0.5, 1.0),
    )

    file_name = "random_generated.tree"

    with open(f"{file_name}", "w") as f:
        f.write(t.write(format=1))


def delete_files_with_prefix(directory, prefix):
    """
    Delete all files in a directory that start with a specific prefix.

    Parameters:
    - directory: The path to the directory to search.
    - prefix: The file prefix to match.
    """
    # List all files in the directory
    all_files = os.listdir(directory)
    # Filter files based on the prefix and delete them
    for f in all_files:
        if f.startswith(prefix):
            file_path = os.path.join(directory, f)
            os.remove(file_path)
            print(f"Deleted: {file_path}")


def delete_directories_with_prefix(directory, prefix):
    """
    Delete all directories in a specified directory that start with a specific prefix.

    Parameters:
    - directory: The path to the directory to search.
    - prefix: The directory prefix to match.
    """
    # List all entries in the directory
    all_entries = os.listdir(directory)

    # Filter directories based on the prefix and delete them
    for entry in all_entries:
        entry_path = os.path.join(directory, entry)
        if os.path.isdir(entry_path) and entry.startswith(prefix):
            shutil.rmtree(entry_path)
            print(f"Deleted directory: {entry_path}")


def move_files_to_directory(prefix, alignment_num):
    # Assuming that the alignments and trees have specific extensions, e.g., ".phy" and ".tree"
    # Adjust these as needed
    alignment_extension = ".phy"

    for i in range(1, alignment_num + 1):
        alignment_file = f"{prefix}_{i}{alignment_extension}"

        # Create a new directory for the alignment and tree
        new_dir = f"{prefix}_{i}"
        os.makedirs(new_dir, exist_ok=True)

        if os.path.exists(os.path.join(new_dir, os.path.basename(alignment_file))):
            os.remove(os.path.join(new_dir, os.path.basename(alignment_file)))

        # Then line 96 remains as:
        shutil.move(alignment_file, new_dir)

    print(f"Moved {alignment_num} alignments and trees to individual directories.")


def execute_satute_cli_for_experiment_folders(main_directory, iqtree_executable):
    """
    Execute the satute_cli.py script for each experiment subdirectory in the main directory.

    Parameters:
    - main_directory (str): Path to the main directory containing experiment subdirectories.
    - iqtree_executable (str): Path to the IQ-TREE executable.
    - model (str): Substitution model to be used.

    Returns:
    - None
    """

    # List all subdirectories in the main directory
    subdirectories = [
        os.path.join(main_directory, d)
        for d in os.listdir(main_directory)
        if os.path.isdir(os.path.join(main_directory, d))
    ]

    for subdir in subdirectories:
        # Find the alignment file in the subdirectory (assuming they have .fasta extension)
        alignment_files = [f for f in os.listdir(subdir) if f.endswith(".phy")]

        # If there's an alignment file, process it
        if alignment_files:
            alignment_path = os.path.join(subdir, alignment_files[0])

            # Construct the command
            cmd = [
                "python",
                "./../satute_cli.py",
                "-iqtree",
                iqtree_executable,
                "-msa",
                alignment_path,
                "-tree",
                "./random_generated.tree",
            ]

            # Execute the command
            subprocess.run(cmd)


def move_phy_files_to_subdirectories(base_directory):
    """
    Moves every .phy file in the given directory into its own subdirectory.

    Parameters:
    - base_directory (str): The path to the directory containing .phy files.
    """
    # List all files in the base directory
    files_in_directory = os.listdir(base_directory)

    # Filter to get only .phy files
    phy_files = [f for f in files_in_directory if f.endswith(".phy")]

    # Move each .phy file into its own directory
    for phy_file in phy_files:
        # Extract the filename without the extension to use as directory name
        directory_name = os.path.splitext(phy_file)[0]

        # Create the directory if it doesn't exist
        new_directory_path = os.path.join(base_directory, directory_name)
        if not os.path.exists(new_directory_path):
            os.makedirs(new_directory_path)

        # Path to the destination file
        destination_path = os.path.join(new_directory_path, phy_file)

        # If the file already exists in the destination directory, remove it
        if os.path.exists(destination_path):
            os.remove(destination_path)

        # Move the .phy file to the new directory
        shutil.move(os.path.join(base_directory, phy_file), new_directory_path)


def move_tree_files_to_subdirectories(base_directory):
    """
    Moves every .phy file in the given directory into its own subdirectory.

    Parameters:
    - base_directory (str): The path to the directory containing .phy files.
    """
    # List all files in the base directory
    files_in_directory = os.listdir(base_directory)

    # Filter to get only .phy files
    phy_files = [f for f in files_in_directory if f.endswith(".tree")]

    # Move each .phy file into its own directory
    for phy_file in phy_files:
        # Extract the filename without the extension to use as directory name
        directory_name = os.path.splitext(phy_file)[0]

        # Create the directory if it doesn't exist
        new_directory_path = os.path.join(base_directory, directory_name)
        if not os.path.exists(new_directory_path):
            os.makedirs(new_directory_path)

        # Path to the destination file
        destination_path = os.path.join(new_directory_path, phy_file)

        # If the file already exists in the destination directory, remove it
        if os.path.exists(destination_path):
            os.remove(destination_path)

        # Move the .phy file to the new directory
        shutil.move(os.path.join(base_directory, phy_file), new_directory_path)


def move_files_by_suffix(source_folder, target_folder, suffix):
    """
    Move all files with a specific suffix from a source folder to a target folder.

    Parameters:
    - source_folder (str): The directory path where to look for files.
    - target_folder (str): The directory path where to move the files.
    - suffix (str): The file suffix to match against.

    Returns:
    - int: The number of files moved.
    """
    # Initialize counter for moved files
    moved_files_count = 0

    # Create target folder if it doesn't exist
    if not os.path.exists(target_folder):
        os.makedirs(target_folder)

    # Loop through all files in the source folder
    for filename in os.listdir(source_folder):
        # Construct the full path of the file
        source_file_path = os.path.join(source_folder, filename)

        # Check if the file has the specified suffix and is not a directory
        if filename.endswith(suffix) and os.path.isfile(source_file_path):
            # Construct the full path to move the file to
            target_file_path = os.path.join(target_folder, filename)

            # Move the file
            shutil.move(source_file_path, target_file_path)
            moved_files_count += 1

    return moved_files_count


# Updating the process function to include the move step
def power_test_one_process(model, sequence_length, alignment_num, taxa_count):
    prefix = f"{model}-{sequence_length}"
    # execute_command("source ./env/bin/activate")
    generate_tree(taxa_count)
    os.chdir("./power_test_one")
    execute_command(
        f"iqtree --alisim {model}-{sequence_length} -t random_generated.tree --length {sequence_length} --seqtype DNA --num-alignments {alignment_num} -m {model}"
    )
    # Move the generated files to individual directories
    move_files_to_directory(prefix, alignment_num)
    execute_satute_cli_for_experiment_folders("./", iqtree)
    os.chdir("./..")


def consolidate_data(main_directory, model="JC"):
    """
    Extract and consolidate the arbitrary .csv files from each experiment subdirectory.
    """
    # List all subdirectories in the main directory
    subdirectories = [
        os.path.join(main_directory, d)
        for d in os.listdir(main_directory)
        if os.path.isdir(os.path.join(main_directory, d))
    ]

    print(main_directory)
    # Collect the data from each subdirectory
    all_data_frames = []

    print(subdirectories)
    for subdir in subdirectories:
        result_file_paths = fetch_files_with_prefix_and_extension_absolute(
            subdir, model, ".csv"
        )
        for result_file_path in result_file_paths:
            if os.path.exists(result_file_path):
                df = pd.read_csv(result_file_path)
                # Extract sequence length from the subdirectory name and add it as a new column
                sequence_length = int(subdir.split("-")[1].split("_")[0])
                df["sequence_length"] = sequence_length
                all_data_frames.append(df)

    # Concatenate all individual dataframes to get the consolidated dataframe
    consolidated_df = pd.concat(all_data_frames, ignore_index=True)
    return consolidated_df


def save_consolidated_data(df, file_name="./consolidated_results.csv"):
    """
    Save the consolidated dataframe as a CSV file.
    """
    df.to_csv(file_name, index=False)


def plot_violin_graphs_to_directory(df, output_directory="saturation_plots"):
    """
    Plot violin plots for delta values for each edge/branch and save them to a directory.
    """
    # Ensure the output directory exists

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    branches = df["edge"].unique()

    for branch in branches:
        branch_data = df[df["edge"] == branch]

        plt.figure(figsize=(10, 6))
        sns.violinplot(x="sequence_length", y="delta", data=branch_data)
        plt.title(f"Delta values for branch: {branch}")

        # Save the plot to the output directory
        plt.savefig(os.path.join(output_directory, f"delta_values_{branch}.png"))
        plt.close()


def split_dataframe_by_column_value(df, column):
    """
    Split a DataFrame based on unique values in a column.

    Args:
    - df (pd.DataFrame): The DataFrame to be split.
    - column (str): The column name based on which the DataFrame will be split.

    Returns:
    - dict: A dictionary with unique column values as keys and sub-dataframes as values.
    """
    grouped = df.groupby(column)
    df_dict = {key: grouped.get_group(key) for key in grouped.groups}
    return df_dict


def power_test_two_process(
    model,
    sequence_length_begin,
    sequence_length_end,
    step,
    taxa_count,
    num_alignment=2,
):
    os.chdir("./power_test_two")
    generate_tree(taxa_count, file_name="./random_generated.tree")

    delete_directories_with_prefix("./", f"{model}")

    for sequence_length in range(sequence_length_begin, sequence_length_end, step):
        prefix = f"{model}-{sequence_length}"
        execute_command(
            f"iqtree --alisim {prefix} -t random_generated.tree --length {sequence_length} --num-alignments {num_alignment} --seqtype DNA  -m {model}"
        )
        move_phy_files_to_subdirectories("./")

    execute_satute_cli_for_experiment_folders("./", iqtree)
    # Define the main directory containing the experiment subdirectories

    main_directory = "./"
    consolidated_data = consolidate_data(main_directory)
    save_consolidated_data(consolidated_data)
    delete_directories_with_prefix("./", "saturation_plots")
    plot_violin_graphs_to_directory(consolidated_data)


def delete_files_by_suffix(folder_path, suffix):
    """
    Delete all files in a folder with a specific suffix (or file extension).

    Parameters:
    - folder_path (str): The directory path where to look for files.
    - suffix (str): The file suffix to match against.

    Returns:
    - int: The number of files deleted.
    """

    # Initialize counter for deleted files
    deleted_files_count = 0

    # Loop through all files in the folder
    for filename in os.listdir(folder_path):
        # Construct the full path of the file
        file_path = os.path.join(folder_path, filename)

        # Check if the file has the specified suffix and is not a directory
        if filename.endswith(suffix) and os.path.isfile(file_path):
            # Delete the file
            os.remove(file_path)
            deleted_files_count += 1

    return deleted_files_count


# Example usage
# folder_path = "/path/to/folder"
# print(f"Deleted {deleted_count} files with suffix '{suffix}' from folder '{folder_path}'")


def assign_branch_length(tree, factor=2.0, branch_idx=None):
    """
    Increase the length of a branch in the tree.

    Parameters:
    - tree (Tree): The ETE3 Tree object to modify.
    - factor (float): The factor by which to increase the branch length.
    - branch_idx (int, optional): The index of the branch to modify.
                                  If None, a random branch is selected.

    Returns:
    - Tree: The modified tree.
    """

    # Convert tree to a list of nodes (excluding the root)
    branches = [node for node in tree.traverse() if not node.is_root()]

    # Randomly choose a branch if branch_idx is not provided
    if branch_idx is None:
        branch_idx = random.randint(0, len(branches) - 1)

    # Select the branch to modify
    branch_to_modify = branches[branch_idx]

    # Modify the branch length
    if branch_to_modify.dist is not None:
        branch_to_modify.dist = factor

    return tree


def apply_satute_cli(tree_file, msa_folder, output_folder, model):
    """
    Apply the satute_cli method for a given tree and multiple MSAs.

    Parameters:
    - tree_file (str): The file path of the given tree in Newick format.
    - msa_folder (str): The directory where all the MSA files are located.
    - output_folder (str): The directory to store the results.

    Returns:
    - None
    """

    # Create output folder if it doesn't exist
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Loop through all files in the MSA folder
    for msa_file in os.listdir(msa_folder):
        # Construct the full path of the MSA file
        msa_file_path = os.path.join(msa_folder, msa_file)

        # Skip if not a file or doesn't end with .phy (assuming PHYLIP format)
        if not os.path.isfile(msa_file_path) or not msa_file.endswith(".phy"):
            continue

        print(f"Processing MSA file: {msa_file_path}")
        # Construct the command
        cmd = [
            "python",
            "./../satute_cli.py",
            "-tree",
            tree_file,
            "-iqtree",
            iqtree,
            "-msa",
            msa_file_path,
            "-model",
            model,
        ]

        print(" ".join(cmd))

        # Execute the command
        subprocess.run(cmd)


def summarize_saturation_results(base_folder, edge="(Node1*, s5)", branch_length=1.0):
    """
    Summarize saturation results from CSV files in subfolders.

    Parameters:
    - base_folder (str): The base directory containing the subfolders with CSV files.

    Returns:
    - summary_dict (dict): A dictionary summarizing the saturation results.
    """
    summary_dict = defaultdict(lambda: defaultdict(int))

    # List all relevant subfolders in the base folder
    all_subfolders = [
        folder
        for folder in os.listdir(base_folder)
        if folder.startswith("given_tree_branch") and not folder.endswith(".log")
    ]

    for subfolder in all_subfolders:
        length = subfolder.split("_")[4]
        # Define the path to this subfolder
        subfolder_path = os.path.join(base_folder, subfolder)

        # List all files in this subfolder
        all_files = os.listdir(subfolder_path)

        # Filter out .csv files
        csv_files = [file for file in all_files if file.endswith(".csv")]

        for csv_file in csv_files:
            # Define the path to the CSV file
            csv_file_path = os.path.join(subfolder_path, csv_file)

            # Read the CSV file into a DataFrame
            df = pd.read_csv(csv_file_path)

            # Loop through the DataFrame to summarize saturation results
            for index, row in df.iterrows():
                if edge in row["edge"]:
                    print(row["edge"], row["p_value"], row["result_test"], length)


def power_test_three_process(
    model, start, end, sequence_length, step, num_alignment=20
):
    os.chdir("./power_test_three")
    delete_directories_with_prefix("./", f"given_tree_branch")

    # Create a tree using the Newick format
    # Here, numbers in the format represent branch lengths
    newick_str = "(((s1:1,s2:2):0.5,(s3:1,s4:1):0.5),s5:1);"
    tree = Tree(newick_str)
    branch_index = 1

    # Make a specific branch longer (e.g., the 2nd branch)
    for increment_branch_length in np.arange(start, end + step, step):
        tree = assign_branch_length(
            tree, factor=increment_branch_length, branch_idx=branch_index
        )

        tree.write(
            outfile=f"./given_tree_branch_{branch_index}_{increment_branch_length}.tree",
            format=1,
        )

        move_tree_files_to_subdirectories("./")

        execute_command(
            f"iqtree --alisim given_tree_branch_{branch_index}_{increment_branch_length} --num-alignments {num_alignment} -t ./given_tree_branch_{branch_index}_{increment_branch_length}/given_tree_branch_{branch_index}_{increment_branch_length}.tree --length {sequence_length} --seqtype DNA  -m {model} --quiet"
        )

        move_files_by_suffix(
            "./",
            f"./given_tree_branch_{branch_index}_{increment_branch_length}",
            ".phy",
        )

        apply_satute_cli(
            f"./given_tree_branch_{branch_index}_{increment_branch_length}/given_tree_branch_{branch_index}_{increment_branch_length}.tree",
            f"./given_tree_branch_{branch_index}_{increment_branch_length}",
            f"./given_tree_branch_{branch_index}_{increment_branch_length}",
            model,
        )

        summarize_saturation_results("./", "(Node1*, s5)", increment_branch_length)


if __name__ == "__main__":
    # model_set = ["GTR", "GTR+G4", "GTR+G8", "JC", "JC+G4", "JC+G8"]
    # delete_files_with_prefix("./saturation_test", prefix)
    # power_test_one_process("JC", 100, 3)
    power_test_three_process(
        model="JC", start=1, end=10, step=1, sequence_length=500, num_alignment=2
    )

    # power_test_two_process(
    #     "JC",
    #     sequence_length_begin=100,
    #     sequence_length_end=300,
    #     step=100,
    #     taxa_count=5,
    #     num_alignment=2,
    # )
