from ete3 import Tree
import numpy as np
from scipy.sparse.linalg import expm
from scipy.linalg import expm
from satute_rate_categories_and_alignments import cut_alignment_columns, read_alignment_file
from Bio import AlignIO

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


def calculate_partial_likelihood_per_site(tree, pattern, rate_matrix, state_frequencies, dimension):
    for node in tree.traverse("postorder"):
        if node.is_leaf():
            node.add_feature(
                "partial_likelihood", get_initial_likelihood_vector(get_sequence_by_taxon(pattern,node.name), state_frequencies)
            )
        else:
            likelihood_factors = []
            for child_node in node.children:

                branch_length = child_node.dist
                transition_matrix = get_transition_matrix(rate_matrix, branch_length)
                factor = []

                for i in range(dimension):

                    component = (np.array(transition_matrix[i]) @ child_node.partial_likelihood)

                    factor.append(component)

                likelihood_factors.append(factor)

            partial_likelihood_vector = np.ones(dimension)

            for factor in likelihood_factors:
                partial_likelihood_vector = partial_likelihood_vector * factor

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


test_tree = Tree("(t7:0.0000029501,(((t3:1.1860038987,t5:0.3574240070)Node4:0.3519068150,t6:0.0000026009)Node3:0.9333701329,t1:0.0000020740)Node2:1.0874051066,(t4:1.7362743020,t2:3.4573102784)Node5:0.0372190293)Node1;", format=1)
rate_matrix = np.asmatrix([[-3.0, 1, 1, 1], [1, -3, 1, 1], [1, 1, -3, 1], [1, 1, 1, -3]], dtype=np.float64)
state_space = {"A": 0, "C": 1, "G": 2, "T": 3}
dimension = len(state_space)
alignment = read_alignment_file("./Clemens/toy_example_JC/toy_example_ntaxa_7_run_5-alignment.phy")

for i in range(alignment.get_alignment_length()): 
    pattern = cut_alignment_columns(alignment, [i])
    calculate_partial_likelihood_per_site(test_tree, pattern, rate_matrix, state_space, dimension)

# for node in test_tree.traverse("postorder"):
#    print(node, node.partial_likelihood)
#     print(node.dist)
