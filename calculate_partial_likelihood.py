import numpy as np
from scipy.sparse.linalg import expm
from scipy.linalg import expm
from satute_rate_categories_and_alignments import (
    read_alignment_file,
)
import os
from pathlib import Path
import pandas as pd
from satute_trees_and_subtrees import parse_newick_file


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


# get transition matrix using matrix exponential
def get_transition_matrix(rate_matrix, branch_length):
    transition_matrix = expm(rate_matrix * branch_length)
    return transition_matrix


def get_initial_likelihood_vector(state, state_frequencies):
    nucleotide_code_vector = {
        "A": [1, 0, 0, 0],
        "C": [0, 1, 0, 0],
        "G": [0, 0, 1, 0],
        "T": [0, 0, 0, 1],
        "N": [1, 1, 1, 1],
        "-": [1, 1, 1, 1],
    }
    return nucleotide_code_vector[state]


def calculate_partial_likelihood_per_site(
    tree, pattern, rate_matrix, state_frequencies, dimension
):
    partial_likelihood_dictionary = {}

    for node in tree.traverse("postorder"):
        if node.is_leaf():
            node.add_feature(
                "partial_likelihood",
                get_initial_likelihood_vector(
                    get_sequence_by_taxon(pattern, node.name), state_frequencies
                ),
            )
        else:
            likelihood_factors = []

            for child_node in node.children:
                branch_length = child_node.dist

                transition_matrix = get_transition_matrix(rate_matrix, branch_length)
                factor = []

                for i in range(dimension):
                    component = (
                        np.array(transition_matrix[i]) @ child_node.partial_likelihood
                    )
                    factor.append(component)

                likelihood_factors.append(factor)

            partial_likelihood_vector = np.ones(dimension)

            for factor in likelihood_factors:
                partial_likelihood_vector = partial_likelihood_vector * factor

            partial_likelihood_dictionary[
                node.write(format=1)
            ] = partial_likelihood_vector

            node.add_feature("partial_likelihood", partial_likelihood_vector)


def get_sequence_by_taxon(alignment, taxon_name):
    """
    Retrieve a sequence from a multiple sequence alignment based on the taxon name.

    Parameters:
        - alignment_file (str): Path to the alignment file (in FASTA format).
        - taxon_name (str): The name of the taxon to search for.

    Returns:
        - str: The sequence corresponding to the taxon name if found, or None if not found.
    """
    # Iterate over the sequences in the alignment
    for record in alignment:
        if taxon_name in record.id:
            # Sequence with matching taxon name found
            return str(record.seq)

    # No sequence found with the specified taxon name
    return None


def run_calculation_partial_likelihood_for_directories(
    path_folder, number_rates, chosen_rate, rate_matrix, dimension, state_space
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
        sub_dir for sub_dir in Path(subtrees_dir).iterdir() if sub_dir.is_dir()
    ]

    for subtree_dir in subtrees_dirs:
        sub_tree = parse_newick_file(f"{subtree_dir}/subtree.treefile")

        sub_alignment = read_alignment_file(os.path.join(subtree_dir, "subtree.fasta"))

        partial_likelihoods = []

        for i in range(sub_alignment.get_alignment_length()):
            pattern = sub_alignment[:, (i + 1) - 1 : (i + 1)]

            calculate_partial_likelihood_per_site(
                sub_tree, pattern, rate_matrix, state_space, dimension
            )

            partial_likelihoods.append(
                {
                    "Site": i + 1,
                    "Node": sub_tree.name,
                    "State": pattern[0][0],
                    "p_A": sub_tree.partial_likelihood[0],
                    "p_C": sub_tree.partial_likelihood[1],
                    "p_G": sub_tree.partial_likelihood[2],
                    "p_T": sub_tree.partial_likelihood[3],
                }
            )

        partial_likelihoods = pd.DataFrame(partial_likelihoods)

        partial_likelihoods.to_csv(f"{subtree_dir}/output.state", sep="\t", index=False)


if __name__ == "__main__":
    # test_tree = Tree("(t7:0.0000029501,(((t3:1.1860038987,t5:0.3574240070)Node4:0.3519068150,t6:0.0000026009)Node3:0.9333701329,t1:0.0000020740)Node2:1.0874051066,(t4:1.7362743020,t2:3.4573102784)Node5:0.0372190293)Node1;", format=1)

    rate_matrix = np.asmatrix(
        [[-3.0, 1, 1, 1], [1, -3, 1, 1], [1, 1, -3, 1], [1, 1, 1, -3]], dtype=np.float64
    )
    state_space = {"A": 0, "C": 1, "G": 2, "T": 3}
    dimension = len(state_space)
    run_calculation_partial_likelihood_for_directories(
        "./test/octo-kraken-msa-test/", 4, 1, rate_matrix, dimension, state_space
    )
