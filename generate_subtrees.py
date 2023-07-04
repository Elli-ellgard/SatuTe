from ete3 import Tree
import os
import pandas as pd
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import numpy as np
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import shutil
from pathlib import Path
import subprocess
from multiprocessing import Pool
from functools import partial
from scipy.sparse.linalg import expm
from satute_repository import (
    parse_state_frequencies,
    parse_rate_matrices,
    parse_rate_parameters,
)
from satute_util import (
    node_type,
    branch_lengths,
)


def run_iqtree_for_each_clade(
    path_folder, number_rates, chosen_rate, iqtree_path, model
):
    """Prepares necessary information and runs the IQ-TREE."""
    path_new_folder = ""

    # Condition to define path and model
    if number_rates > 1:
        path_new_folder = os.path.join(
            path_folder, f"subsequence{chosen_rate}", "subtrees"
        )

    else:
        path_new_folder = os.path.join(path_folder, "subtrees")

    # Read model and frequency
    # Check if path_new_folder exists and is a directory
    # if not os.path.isdir(path_new_folder):
    #    raise NotADirectoryError(
    #        f"The path '{path_new_folder}' does not exist or is not a directory."
    #    )

    # Iterate over each clade directory
    for clade_dir in Path(path_new_folder).iterdir():
        if clade_dir.is_dir() and "external" not in str(clade_dir):
            cmd = [
                iqtree_path,
                "-s",
                "subtree.fasta",
                "-te",
                "subtree.tree",
                "-m",
                model,
                "-asr",
                "-blfix",
                "-o",
                "FOO",
                "-pre",
                "output",
                "-redo",
                "-quiet",
            ]

            # Check if sequence and tree files exist in the clade directory
            if not os.path.isfile(
                os.path.join(clade_dir, "subtree.treefile")
            ) or not os.path.isfile(os.path.join(clade_dir, "subtree.treefile")):
                raise FileNotFoundError(
                    "Either sequence.txt or tree.txt file does not exist in directory: "
                    + str(clade_dir)
                )

            print(clade_dir)

            # Run the command
            # result = subprocess.run(cmd, cwd=clade_dir)

            # Check if the command was successful
            # if result.returncode != 0:
            #    raise RuntimeError(
            #        f"The command '{' '.join(cmd)}' failed with return code: {result.returncode}"
            #    )


def name_nodes_by_level_order(tree):
    i = 1
    for node in tree.traverse("levelorder"):
        if not node.is_leaf():
            node.name = f"Node{i}*"
            i += 1
    return tree


def get_all_subtrees(tree):
    subtrees = []
    for node in tree.traverse("levelorder"):
        if not node.is_root():
            subtrees.append(node)
    return subtrees


def get_opposite_subtree(tree, subtree):
    # Create a copy of the tree so as not to modify the original
    tree_copy = Tree(tree.write())
    tree_copy = name_nodes_by_level_order(tree_copy)
    # Get the node in the copied tree that corresponds to the root of the subtree
    node = tree_copy.search_nodes(name=subtree.name)[0]
    # Detach this node (and its descendants) from the tree
    node.detach()
    return tree_copy


def get_leaves(tree):
    leaves = []
    for leaf in tree:
        if leaf.is_leaf():
            leaves.append(leaf.name)
    return leaves


def filter_alignment_by_ids(alignment, ids):
    """
    Filter a MultipleSeqAlignment object to include only sequences with specific IDs.

    Parameters:
    alignment (MultipleSeqAlignment): The alignment to filter.
    ids (list of str): The IDs of the sequences to include.

    Returns:
    MultipleSeqAlignment: The filtered alignment.
    """

    # Filter the alignment
    filtered_alignment = MultipleSeqAlignment(
        [record for record in alignment if record.id in ids]
    )

    return filtered_alignment


def generate_subtree_pair(subtrees, t):
    subtree_pairs = []
    for subtree in subtrees:
        subtree_copy = t.search_nodes(name=subtree.name)[0]
        opposite_subtree = get_opposite_subtree(t, subtree_copy)

        subtree_pair_entry = {
            "trees": (subtree_copy, opposite_subtree),
        }

        subtree_pairs.append(subtree_pair_entry)

    return subtree_pairs


""" TODO: 
    - check if the transition matrix is a real transition matrix
    - we should think about create a directory of the state_space at the beginning of the program and use it for dimension stuff
    - generalizes the functions below for other dimension (now only dimension 4)
    - the posterior distribution is not correct at the moment
"""


# get transition matrix using matrix exponential
def get_transition_matrix(rate_matrix, branch_length):
    return expm(rate_matrix * branch_length)


# calculate for the partial likelihoods the needed factors
def get_likelihood(transition_matrix, state, state_frequencies):
    likelihood_factors = []
    row = [0, 0, 0, 0]
    state_space = {"A": 0, "C": 1, "G": 2, "T": 3}
    if state == "-":
        row = list(state_frequencies.values())
    else:
        row[state_space[state]] = 1
    for i in range(len(state_space)):
        component = (float)(np.array(transition_matrix[i]) @ np.array(row))
        likelihood_factors.append(component)
    return likelihood_factors


def generate_output_state_file_for_cherry(
    alignment, file_path, state_frequencies, rate_matrix
):
    # get branch lengths
    tree_file = f"{file_path}subtree.treefile"
    T = parse_newick_file(tree_file)
    leaves, internal_nodes = node_type(T)
    vector_branches, vector_distances = branch_lengths(T)

    # get transitions matrix with the time given as branch lengths
    transition_matrices = []
    for idx in range(len(leaves)):
        transition_matrices.append(
            get_transition_matrix(rate_matrix, vector_distances[idx])
        )

        # get the sequences of the leaves in the correct order
        alignments = []
        for leaf in leaves:
            for record in alignment:
                if record.id == leaf:
                    alignments.append(record.seq)

    with open(f"{file_path}output.state", "w") as state_file_writer:
        header = "Node\tSite\tState\tp_A\tp_C\tp_G\tp_T"
        state_file_writer.write(header)
        site_number = len(alignments[0])
        for i in range(site_number):
            likelihood_left = get_likelihood(
                transition_matrices[0], alignments[0][i], state_frequencies
            )
            likelihood_right = get_likelihood(
                transition_matrices[1], alignments[1][i], state_frequencies
            )
            # calculate the partial likelihood vector
            likelihood = np.asarray(likelihood_left) * np.asarray(likelihood_right)

            scale_factor = np.asarray(likelihood) * np.asarray(
                list(state_frequencies.values())
            )

            # transform the partial likelihood vector into the posterior distribution
            distribution = (
                np.asarray(likelihood)
                * np.asarray(list(state_frequencies.values()))
                / scale_factor
            )

            values = "\t".join(
                "{:.5f}".format(distribution[value]) for value in range(4)
            )
            state_file_writer.write(f"\nNode1\t{i + 1}\t{alignments[0][i]}\t{values}")


def generate_output_state_file_for_external_branch(
    alignment, file_path, state_frequencies
):
    with open(f"{file_path}output.state", "w") as state_file_writer:
        header = "Node\tSite\tState\tp_A\tp_C\tp_G\tp_T"
        state_file_writer.write(header)
        for record in alignment:
            sequence = record.seq
            index = 0
            for character in sequence:
                index += 1
                row = {"A": 0, "C": 0, "G": 0, "T": 0}
                if character == "-":
                    state_freq = list(state_frequencies.values())
                    row["A"] = state_freq[0]
                    row["C"] = state_freq[1]
                    row["G"] = state_freq[2]
                    row["T"] = state_freq[3]
                else:
                    row[character] = 1
                values = "\t".join(
                    "{:.5f}".format(row[value]) for value in ["A", "C", "G", "T"]
                )
                state_file_writer.write(f"\nNode1\t{index}\t{character}\t{values}")


def write_subtree_and_sub_alignments(
    generated_subtree_pairs, alignment, state_frequencies, rate_matrix, path_prefix="./"
):
    for i, subtree_pair in enumerate(generated_subtree_pairs):
        first_subtree = subtree_pair["trees"][0]
        second_subtree = subtree_pair["trees"][1]

        first_subtree_leaves = get_leaves(first_subtree)
        first_sub_alignment = filter_alignment_by_ids(alignment, first_subtree_leaves)

        second_subtree_leaves = get_leaves(second_subtree)
        second_sub_alignment = filter_alignment_by_ids(alignment, second_subtree_leaves)

        if len(first_subtree.get_descendants()) + 1 == 1:
            first_subtree_dir = f"{path_prefix}subtrees/branch_{i}_subtree_one_leaf/"

        elif len(first_subtree.get_descendants()) + 1 == 3:
            first_subtree_dir = f"{path_prefix}subtrees/branch_{i}_subtree_one_cherry/"
        else:
            first_subtree_dir = f"{path_prefix}subtrees/branch_{i}_subtree_one/"

        os.makedirs(first_subtree_dir, exist_ok=True)

        # Write the first subtree
        subtree_writer = open(f"{first_subtree_dir}subtree.treefile", "w")
        subtree_writer.write(first_subtree.write(format=1))
        subtree_writer.close()

        # Write the alignment for the first subtree
        AlignIO.write(first_sub_alignment, f"{first_subtree_dir}subtree.fasta", "fasta")

        if len(second_subtree.get_descendants()) + 1 == 1:
            second_subtree_dir = f"{path_prefix}subtrees/branch_{i}_subtree_two_leaf/"

        elif len(second_subtree.get_descendants()) + 1 == 3:
            second_subtree_dir = f"{path_prefix}subtrees/branch_{i}_subtree_two_cherry/"

        else:
            second_subtree_dir = f"{path_prefix}subtrees/branch_{i}_subtree_two/"

        os.makedirs(second_subtree_dir, exist_ok=True)

        # Write the second subtree
        subtree_writer = open(f"{second_subtree_dir}subtree.treefile", "w")
        subtree_writer.write(second_subtree.write(format=1))
        subtree_writer.close()

        # Write the alignment for the second subtree
        AlignIO.write(
            second_sub_alignment, f"{second_subtree_dir}subtree.fasta", "fasta"
        )

        if len(first_subtree.get_descendants()) + 1 == 1:
            generate_output_state_file_for_external_branch(
                first_sub_alignment,
                first_subtree_dir,
                state_frequencies,
            )
        elif len(first_subtree.get_descendants()) + 1 == 3:
            generate_output_state_file_for_cherry(
                first_sub_alignment,
                first_subtree_dir,
                state_frequencies,
                rate_matrix,
            )

        if len(second_subtree.get_descendants()) + 1 == 1:
            generate_output_state_file_for_external_branch(
                second_sub_alignment,
                second_subtree_dir,
                state_frequencies,
            )
        elif len(second_subtree.get_descendants()) + 1 == 3:
            generate_output_state_file_for_cherry(
                second_sub_alignment,
                second_subtree_dir,
                state_frequencies,
                rate_matrix,
            )


def parse_file_to_dataframe(file_path):
    try:
        # Read the file into a dataframe
        df = pd.read_csv(file_path, delimiter="\t")

        return df

    except FileNotFoundError:
        raise Exception(f"File not found: {file_path}")


def parse_newick_file(file_path):
    try:
        # Open the file in read mode
        with open(file_path, "r") as f:
            # Read the content of the file
            newick_string = f.readlines()

        # Parse the newick string into a Tree object
        t = Tree(newick_string[0], format=1)

        # Return the parsed Tree object
        return t

    except FileNotFoundError:
        raise Exception("File not found: " + file_path)


def guess_alignment_format(file_name):
    with open(file_name, "r") as f:
        first_line = f.readline().strip()

    # Now we'll check for various signatures that might indicate the format.
    if first_line.startswith(">"):
        return "fasta"
    elif first_line.startswith("CLUSTAL"):
        return "clustal"
    elif first_line.startswith("# STOCKHOLM"):
        return "stockholm"
    elif first_line.startswith("#NEXUS"):
        return "nexus"
    elif first_line.startswith("PileUp"):
        return "pileup"
    elif first_line[0].isdigit():
        return "phylip"
    else:
        return "Unknown"


def rescale_branch_lengths(tree, rescale_factor):
    """
    Rescales the branch lengths of a phylogenetic tree by a given factor.

    Args:
        tree (ete3.Tree): The input phylogenetic tree object.
        rescale_factor (float): Scaling factor for the branch lengths.

    Returns:
        ete3.Tree: The phylogenetic tree object with rescaled branch lengths.
    """
    # Iterate over tree nodes and rescale branch lengths
    for node in tree.traverse():
        if node.up:
            node.dist *= rescale_factor

    # Return the tree with rescaled branch lengths
    return tree


def generate_write_subtree_pairs_and_msa(
    number_rates,
    t,
    file_path,
    msa_file_name,
    category_rate,
    state_frequencies,
    rate_matrix,
):
    # Write subtree pairs to files, assuming the function is defined elsewhere
    if number_rates == 1:
        try:
            # Call function to get all subtrees from t, assuming the function is defined elsewhere
            subtrees = get_all_subtrees(t)

            ## Discard the first subtree, assuming we don't need it
            # subtrees = subtrees[1:]

            # Generate subtree pairs, assuming the function is defined elsewhere
            generated_subtree_pairs = generate_subtree_pair(subtrees, t)
            alignment = read_alignment_file(msa_file_name)
            write_subtree_and_sub_alignments(
                generated_subtree_pairs,
                alignment,
                state_frequencies,
                rate_matrix,
                f"./{file_path}/",
            )
        except Exception as e:
            print(f"Error occurred during the first iteration: {e}")
    else:
        # Iterate from 0 to number_rates (inclusive)
        for i in range(1, number_rates + 1):
            rate = category_rate[i - 1]["Relative_rate"]
            # Rescale the branch lengths
            rescaled_tree = rescale_branch_lengths(t, rate)
            # Write subtree pairs to files in the new directory
            # Call function to get all subtrees from t, assuming the function is defined elsewhere
            subtrees = get_all_subtrees(rescaled_tree)

            # Discard the first subtree, assuming we don't need it
            # subtrees = subtrees[1:]
            # for tree in subtrees:
            #    print(tree.write())

            # Generate subtree pairs, assuming the function is defined elsewhere
            generated_subtree_pairs = generate_subtree_pair(subtrees, rescaled_tree)

            sub_alignment = read_alignment_file(
                f"./{file_path}/subsequence{i}/rate.fasta"
            )

            write_subtree_and_sub_alignments(
                generated_subtree_pairs,
                sub_alignment,
                state_frequencies,
                rate_matrix,
                path_prefix=f"./{file_path}/subsequence{i}/",
            )


def get_column_names_with_prefix(data_frame, prefix):
    # Filter the columns using the specified prefix
    columns_with_prefix = data_frame.columns[
        data_frame.columns.str.startswith(prefix)
    ].tolist()
    return columns_with_prefix


def build_categories_by_sub_tables(data_frame):
    rate_category_dictionary = {}

    # Assuming you already have a dataframe called 'dataframe'
    # Call the get_columns_with_prefix function to retrieve columns with a specific prefix
    prefix = "p"  # Specify the desired prefix

    columns_with_prefix = get_column_names_with_prefix(data_frame, prefix)

    # Create dictionaries using column names with the specified prefix as names
    rate_category_dictionary = {column: [] for column in columns_with_prefix}

    for index, row in data_frame.iterrows():
        p_row = row.filter(like="p")

        rate_category_dictionary[p_row.idxmax()].append(int(row["Site"]) - 1)

    return rate_category_dictionary


def read_alignment_file(file_name):
    # Guess the format of the file
    file_format = guess_alignment_format(file_name)

    # If the format could not be guessed, raise an error
    if file_format is None:
        raise ValueError("Could not guess the format of the file.")

    # Try to read the file in the guessed format
    try:
        alignment = AlignIO.read(file_name, file_format)
        return alignment
    except Exception as e:
        print(f"An error occurred while reading the file: {str(e)}")


def cut_alignment_columns(alignment, columns):
    # Convert the alignment to a NumPy array for easy column slicing
    alignment_array = np.array([list(rec) for rec in alignment], np.character)
    # Select the specified columns

    selected_columns = alignment_array[:, columns]
    selected_records = []
    for i, rec in enumerate(selected_columns):
        selected_records.append(
            SeqRecord(Seq(rec.tobytes().decode()), id=alignment[i].id)
        )
    selected_alignment = MultipleSeqAlignment(selected_records)

    return selected_alignment


def write_alignment(alignment, file_name, file_format):
    """
    Write a MultipleSeqAlignment object to a file.

    Parameters:
    alignment (MultipleSeqAlignment): The alignment to write.
    file_name (str): The name of the file to write to.
    file_format (str): The format of the alignment file. This must be a string specifying one
                       of the formats that Biopython supports, such as "fasta", "clustal", "phylip", etc.

    Returns:
    int: The number of alignments written (should be 1 for a MultipleSeqAlignment object).
    """

    # Write the alignment to the file
    num_alignments_written = AlignIO.write(alignment, file_name, file_format)

    return num_alignments_written


def delete_directory_contents(directory_path):
    for root, dirs, files in os.walk(directory_path, topdown=False):
        for file in files:
            file_path = os.path.join(root, file)
            os.remove(file_path)
        for dir in dirs:
            dir_path = os.path.join(root, dir)
            shutil.rmtree(dir_path)


def split_msa_into_rate_categories(site_probability, folder_path, msa_file_name):
    sub_category = build_categories_by_sub_tables(site_probability)
    alignment = read_alignment_file(msa_file_name)
    per_category_alignment_dict = {}

    for key, value in sub_category.items():
        per_category_alignment_dict[key] = cut_alignment_columns(alignment, value)

    for key, value in per_category_alignment_dict.items():
        # Make a new directory for this subsequence, if it doesn't exist yet
        os.makedirs(f"./{folder_path}/subsequence{key[1:]}/", exist_ok=True)

        write_alignment(
            value, f"./{folder_path}/subsequence{key[1:]}/rate.fasta", "fasta"
        )


def parse_table(log_file):
    f = open(log_file, "r")
    lines = f.readlines()
    f.close()
    table_data = []

    # Find the start and end indices of the table
    start_index = -1
    end_index = -1

    for i, line in enumerate(lines):
        if line.strip().startswith("Category"):
            start_index = i + 1
        elif line.strip().startswith("MAXIMUM LIKELIHOOD TREE"):
            end_index = i
            break

    if start_index == -1 or end_index == -1:
        raise ValueError("Table not found in the log file.")

    # Parse the table rows
    table_lines = lines[start_index : end_index - 2]

    for line in table_lines:
        line = line.strip()
        if line:
            row = line.split()
            category = row[0]
            relative_rate = float(row[1])
            proportion = float(row[2])
            if category != "0":
                table_data.append(
                    {
                        "Category": category,
                        "Relative_rate": relative_rate,
                        "Proportion": proportion,
                    }
                )
    return table_data


def execute_iqtree(sub_dir, iqtree_path, model_and_frequency):
    # Command to be executed for each clade
    cmd = [
        iqtree_path,
        "-s",
        "subtree.fasta",
        "-te",
        "subtree.treefile",
        "-m",
        model_and_frequency,
        "-asr",
        "-blfix",
        "-pre",
        "output",
        "-redo",
        "-quiet",
    ]

    # Check if sequence and tree files exist in the subdirectory
    if not os.path.isfile(os.path.join(sub_dir, "subtree.fasta")) or not os.path.isfile(
        os.path.join(sub_dir, "subtree.treefile")
    ):
        raise FileNotFoundError(
            f"Either sequence.txt or tree.txt file does not exist in directory: {sub_dir}"
        )

    # Run the IQ-TREE command
    result = subprocess.run(cmd, cwd=sub_dir)

    # Check if the command was successful
    if result.returncode != 0:
        raise RuntimeError(
            f"The command '{' '.join(cmd)}' failed with return code: {result.returncode}"
        )


def run_iqtree_for_each_clade_parallel(
    path_folder, number_rates, chosen_rate, iqtree_path
):
    """
    Run IQ-TREE for each clade directory in parallel.

    Args:
        path_folder (str): Root folder path.
        number_rates (int): Number of rates.
        chosen_rate (int): Chosen rate.
        iqtree_path (str): Path to the IQ-TREE executable.

    Raises:
        FileNotFoundError: If model and frequency file does not exist.
        NotADirectoryError: If the clade directory does not exist or is not a directory.
        FileNotFoundError: If sequence.txt or tree.txt file does not exist in the clade directory.
        RuntimeError: If IQ-TREE command fails.

    """
    model_and_frequency = ""
    clade_dir = ""

    # Determine path and model based on number of rates
    if number_rates > 1:
        subtrees_dir = os.path.join(
            path_folder, f"subsequence{chosen_rate}", "subtrees"
        )
        model_and_frequency_file = os.path.join(
            path_folder, f"subsequence{chosen_rate}", "model.txt"
        )
    else:
        subtrees_dir = os.path.join(path_folder, "subtrees")
        model_and_frequency_file = os.path.join(path_folder, "model.txt")

    # Check if model and frequency file exists
    if not os.path.isfile(model_and_frequency_file):
        raise FileNotFoundError(
            f"The model and frequency file '{model_and_frequency_file}' does not exist."
        )

    # Read model and frequency
    with open(model_and_frequency_file, "r") as toModel:
        model_and_frequency = toModel.readline().strip()

    # Check if clade directory exists and is a directory
    if not os.path.isdir(subtrees_dir):
        raise NotADirectoryError(
            f"The clade directory '{clade_dir}' does not exist or is not a directory."
        )

    # Create a list of subtree directories for internal branches
    subtrees_dirs = [
        sub_dir
        for sub_dir in Path(subtrees_dir).iterdir()
        if sub_dir.is_dir()
        and not os.path.isfile(os.path.join(sub_dir, "output.state"))
    ]

    # Create a partial function with fixed arguments for execute_iqtree
    execute_iqtree_partial = partial(
        execute_iqtree, iqtree_path=iqtree_path, model_and_frequency=model_and_frequency
    )

    # Create a multiprocessing pool and map the execute_iqtree_partial function to clade_dirs
    with Pool() as pool:
        pool.map(execute_iqtree_partial, subtrees_dirs)


def parse_rate_and_frequencies_and_create_model_files(
    path, number_rates, dimension, model="GTR"
):
    """
    Parse the rate parameter and state frequencies from the IQ-TREE log file and
    create model files based on these parameters.
    """

    def _write_model_file(path, content):
        """Helper function to write model and frequency tokens into a file."""
        with open(path, "w") as f:
            f.write(content)

    # Construct the model string with parsed rate parameters
    log_file_path = f"{path}.iqtree"
    model_final = parse_rate_parameters(log_file_path, dimension, model=model)

    # Parse state frequencies from the log content
    state_frequencies = parse_state_frequencies(log_file_path, dimension=dimension)

    # Create a string of state frequencies separated by a space
    concatenated_rates = " ".join(map(str, state_frequencies.values()))

    # Construct command line tokens for model and frequency
    model_and_frequency = f"{model_final}+FU{{{concatenated_rates}}}"

    # Get the directory of the path
    path_folder = os.path.dirname(path)

    # Write the model and frequency tokens into 'model.txt' files
    if number_rates == 1:
        _write_model_file(os.path.join(path_folder, "model.txt"), model_and_frequency)
    else:
        for i in range(1, number_rates + 1):
            subsequence_folder = os.path.join(path_folder, f"subsequence{i}")
            _write_model_file(
                os.path.join(subsequence_folder, "model.txt"), model_and_frequency
            )

    return state_frequencies.values(), model_and_frequency


##############################################################################################################
delete_directory_contents("./test_cladding_and_subsequence")
number_rates = 4
dimension = 4
site_probability = parse_file_to_dataframe("./Clemens/example_4/example.txt.siteprob")
folder_path = "test_cladding_and_subsequence"
msa_file_name = "./Clemens/example_4/example.txt"
t = name_nodes_by_level_order(
    parse_newick_file("./Clemens/example_4/example.txt.treefile")
)
category_rate = parse_table("./Clemens/example_4/example.txt.iqtree")
# Parse state frequencies from the log content
# state_frequencies = parse_rate_and_frequencies_and_create_model_files("./test_cladding_and_subsequence", number_rates, dimension, model="JC")
state_frequencies = parse_state_frequencies(
    "./Clemens/example_4/example.txt.iqtree", dimension=dimension
)
(rate_matrix, phi_matrix) = parse_rate_matrices(
    dimension, "./Clemens/example_4/example.txt"
)
if number_rates != 1:
    split_msa_into_rate_categories(site_probability, folder_path, msa_file_name)

generate_write_subtree_pairs_and_msa(
    number_rates,
    t,
    folder_path,
    msa_file_name,
    category_rate,
    state_frequencies,
    rate_matrix,
)

import shutil

src_path = r"./Clemens/example_4/model.txt"
dst_path = r"./test_cladding_and_subsequence/subsequence1/model.txt"
shutil.copy(src_path, dst_path)

run_iqtree_for_each_clade_parallel(
    "./test_cladding_and_subsequence",
    4,
    1,
    "iqtree2",
)

"""run_iqtree_for_each_clade(
    "./test_cladding_and_subsequence",
    4,
    1,
    "iqtree",
    "GTR{3.9907,5.5183,4.1388,0.4498,16.8174}+FU{0.3547 0.2282 0.1919 0.2252}",
)"""


""""summary TODO:
    - check if the transition matrix is a real transition matrix
    - we should think about create a directory of the state_space at the beginning of the programm and use it for dimension stuff
    - generalise the functions below for other dimension (now only dimension 4)
    - the posterior distribution is not correct at the moment
    
    subtrees: 
    - something is of with the iqtree-runs, it is not clear where the root for each subtree is (you see Node1 a lot)
    - for unrooted trees on n taxa we have to consider 2n-3 branches
    - distinction rooted and unrooted input trees (do you root the tree before using it?) 

"""
