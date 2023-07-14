from ete3 import Tree
import os
import numpy as np
from pathlib import Path
import subprocess
from multiprocessing import Pool
from functools import partial
from scipy.sparse.linalg import expm
from satute_trees_and_subtrees import (
    branch_lengths,
    get_leaves,
)



""" TODO: 
    - we should think about create a directory of the state_space at the beginning of the program and use it for dimension stuff
    - generalizes the functions below for other dimension (now only dimension 4)
"""


# get transition matrix using matrix exponential
def get_transition_matrix(rate_matrix, branch_length):
    return expm(rate_matrix * branch_length)


# calculate for the partial likelihoods the needed factors
def get_likelihood(transition_matrix, state, state_frequencies):
    likelihood_factors = []
    row = [0, 0, 0, 0]
    state_space = {"A": 0, "C": 1, "G": 2, "T": 3}
    if state == "-" or state == "N":
        row = list(state_frequencies.values())
    else:
        row[state_space[state]] = 1
    for i in range(len(state_space)):
        component = (float)(np.array(transition_matrix[i]) @ np.array(row))
        likelihood_factors.append(component)
    return likelihood_factors


# calculate partial likeklihoods and generate output.state file for cherries
def generate_output_state_file_for_cherry(
    alignment, tree, file_path, state_frequencies, rate_matrix
):
    # get branch lengths
    leaves = get_leaves(tree)
    vector_branches, vector_distances = branch_lengths(tree)

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
            scale_factor = np.asarray(likelihood) @ np.asarray(
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

# generate output.state file for leaves
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
                if character == "-" or character == "N":
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

    def remove_all_except(directory, ignore_list):
        for file in os.listdir(directory):
            if file not in ignore_list:
                os.remove(os.path.join(directory, file))

    # Check if sequence and tree files exist in the subdirectory
    if not os.path.isfile(os.path.join(sub_dir, "subtree.fasta")) or not os.path.isfile(
        os.path.join(sub_dir, "subtree.treefile")
    ):
        raise FileNotFoundError(
            f"Either sequence.txt or tree.txt file does not exist in directory: {sub_dir}"
        )

    # Run the IQ-TREE command
    result = subprocess.run(cmd, cwd=sub_dir)
    remove_all_except(sub_dir, ["subtree.fasta", "subtree.treefile", "output.state"])

    # Check if the command was successful
    if result.returncode != 0:
        raise RuntimeError(
            f"The command '{' '.join(cmd)}' failed with return code: {result.returncode}"
        )


def run_iqtree_for_each_subtree_parallel(
    path_folder, number_rates, chosen_rate, iqtree_path, model_and_frequency=""
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

    subtrees_dir = ""

    # Determine path and model based on number of rates
    if number_rates > 1:
        subtrees_dir = os.path.join(
            path_folder, f"subsequence{chosen_rate}", "subtrees"
        )

    else:
        subtrees_dir = os.path.join(path_folder, "subtrees")

    # Check if clade directory exists and is a directory
    if not os.path.isdir(subtrees_dir):
        raise NotADirectoryError(
            f"The subtrees directory '{subtrees_dir}' does not exist or is not a directory."
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
