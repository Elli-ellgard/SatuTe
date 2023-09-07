import numpy as np
import pandas as pd
from typing import Optional
from functools import cache
from scipy.sparse.linalg import expm
import dataclasses
from ete3 import Tree
from satute_rate_categories_and_alignments import (
    parse_category_rates,
    read_alignment_file,
    split_msa_into_rate_categories_in_place,
)
from satute_trees_and_subtrees import rescale_branch_lengths, parse_newick_file
from satute_test_statistic_using_partial_likelihood import calculate_test_statistic
from satute_test_statistic_using_posterior_distribution import (
    calculate_test_statistic_posterior_distribution,
)
from satute_repository import parse_rate_matrices_from_file
from satute_util_new import spectral_decomposition

NUCLEOTIDE_CODE_VECTOR = {
    "A": [1, 0, 0, 0],
    "C": [0, 1, 0, 0],
    "G": [0, 0, 1, 0],
    "T": [0, 0, 0, 1],
    "U": [0, 0, 0, 1],
    "R": [1, 0, 1, 0],
    "Y": [0, 1, 0, 1],
    "K": [0, 0, 1, 1],
    "M": [1, 1, 0, 0],
    "S": [0, 1, 1, 0],
    "W": [1, 0, 0, 1],
    "B": [0, 1, 1, 1],
    "D": [1, 0, 1, 1],
    "H": [1, 1, 0, 1],
    "V": [1, 1, 1, 0],
    "N": [1, 1, 1, 1],
    "-": [1, 1, 1, 1],
}

RATE_MATRIX = np.array([[-3, 1, 1, 1], [1, -3, 1, 1], [1, 1, -3, 1], [1, 1, 1, -3]])


class RateMatrix:
    """Class representing a rate matrix for nucleotide substitutions."""

    def __init__(self, rate_matrix):
        """Initialize the RateMatrix with a given rate matrix."""
        self.rate_matrix = rate_matrix

    def __hash__(self):
        """Return the hash value of the RateMatrix."""
        return id(self)


class DirectedAcyclicGraph:
    """Class representing a directed acyclic graph (DAG) of interconnected nodes."""

    def __init__(self, edges):
        """Initialize the DAG with a list of directed edges."""
        for left, right, branch_length in edges:
            left.connect(right, branch_length)
            right.connect(left, branch_length)
        self.edges = edges

    def get_edges(self):
        """Return the list of edges in the DAG."""
        return self.edges

    def get_branch_length(self, edge):
        """Return the branch length of a given edge in the DAG."""
        return self.branch_lengths[edge]


class Node:
    """Class representing a node in the DAG."""

    name: str
    state: Optional[np.array] = dataclasses.field(default_factory=lambda: None)
    connected: dict = dataclasses.field(default_factory=dict)

    def __init__(self, name, state=None, connected=None):
        """Initialize a Node with a name, optional state vector, and optional connections."""
        self.name = name
        self.state = state
        self.connected = connected or {}

    def __hash__(self):
        """Return the hash value of the Node."""
        return id(self)

    def __eq__(self, other):
        """Check if two nodes are equal."""
        return self is other

    def connect(self, other_node, branch_length):
        """Connect the current node to another node with a given branch length."""
        self.connected[other_node] = branch_length

    def is_leaf(self):
        """Check if the node is a leaf (has a state vector)."""
        return self.state is not None

    def __repr__(self):
        """Return the string representation of the Node."""
        return f"Node({self.name})"


def get_initial_likelihood_vector(state):
    """Get the initial likelihood vector for a given nucleotide state."""
    return np.array(NUCLEOTIDE_CODE_VECTOR[state]).T


def name_nodes_by_level_order(tree):
    """Name nodes in a tree based on level-order traversal."""
    i = 1
    for node in tree.traverse("levelorder"):
        if not node.is_leaf():
            node.name = f"Node{i}"
            i += 1
    return tree


def convert_ete3_tree_to_directed_acyclic_graph(
    tree, msa_column, alignment_look_up_table
):
    """
    Convert an ETE3 tree to a DirectedAcyclicGraph representation.

    This function transforms an ETE3 tree into a directed acyclic graph (DAG) representation
    by creating Node objects for each node in the tree and establishing connections between
    them based on the tree structure. Leaf nodes are initialized with initial likelihood vectors
    derived from the sequence alignment data.

    Args:
        tree (Tree): An ETE3 tree object representing the phylogenetic tree.
        msa_column (np.array): A single column of the multiple sequence alignment for a site.
        alignment_look_up_table (dict): A lookup table mapping sequence IDs to alignment indices.

    Returns:
        DirectedAcyclicGraph: A directed acyclic graph representation of the phylogenetic tree.

    Note:
        - The ETE3 tree should be properly constructed and rooted.
        - The msa_column parameter provides data for a specific site in the alignment.
        - The alignment_look_up_table aids quick retrieval of alignment record indices.
        - The function creates Node and DirectedAcyclicGraph objects to represent the graph structure.
    """
    node_dictionary = {}

    # Create nodes for each node in the ETE3 tree using a level-order traversal.
    for node in tree.traverse("levelorder"):
        if node.is_leaf():
            # Initialize leaf nodes with initial likelihood vectors from the alignment.
            node_dictionary[node.name] = Node(
                node.name,
                get_initial_likelihood_vector(
                    msa_column[alignment_look_up_table[node.name]].seq
                ),
            )
        else:
            # Non-leaf nodes are initialized with empty states.
            node_dictionary[node.name] = Node(node.name)

    edge_list = []
    # Create edges between nodes based on the tree's hierarchical relationships.
    for node in tree.traverse("levelorder"):
        for child_node in node.children:
            edge_list.append(
                (
                    node_dictionary[node.name],
                    node_dictionary[child_node.name],
                    child_node.dist,
                )
            )

    # Return a DirectedAcyclicGraph representation of the tree using the constructed edges.
    return DirectedAcyclicGraph(edge_list)


def get_alignment_look_up_table(alignment):
    """
    Create a lookup table for quick retrieval of alignment record indices.

    This function generates a lookup table that maps sequence identifiers (record IDs)
    to their corresponding indices in the alignment. The lookup table facilitates quick
    access to alignment records during subsequent computations.

    Args:
        alignment (list): A list of sequence alignment records (e.g., Bio.SeqRecord objects).

    Returns:
        dict: A dictionary mapping record IDs to their corresponding indices.

    Note:
        - The alignment records should have unique identifiers.
        - The returned lookup table aids efficient indexing of records during calculations.
    """
    alignment_look_up_table = {}
    i = 0
    for record in alignment:
        alignment_look_up_table[record.id] = i
        i = i + 1
    return alignment_look_up_table


def calculate_partial_likelihoods_for_sites(tree, alignment, rate_matrix):
    """
    Calculate partial likelihoods for each site in the alignment given a phylogenetic tree.

    This function computes partial likelihoods for each site in the sequence alignment
    using the provided phylogenetic tree and rate matrix. Partial likelihoods represent
    the conditional probabilities of observing the data at each site, considering the
    evolutionary history and substitution rates.

    Args:
        tree (DirectedAcyclicGraph): The directed acyclic graph representing the phylogenetic tree.
        alignment (np.array): The sequence alignment data for multiple taxa and sites.
        rate_matrix (RateMatrix): The rate matrix representing nucleotide substitution rates.
    Returns:
        dict: A dictionary containing partial likelihood information for each site.
              The dictionary structure is as follows:
              {
                  "(Node1, Node2)": {
                      "left": [
                          {
                              "Node": node_name,
                              "Site": site_index,
                              "branch_length": branch_length,
                              "p_A": likelihood_A,
                              "p_C": likelihood_C,
                              "p_G": likelihood_G,
                              "p_T": likelihood_T,
                          },
                          ...
                      ],
                      "right": [
                          {
                              "Node": node_name,
                              "Site": site_index,
                              "branch_length": branch_length,
                              "p_A": likelihood_A,
                              "p_C": likelihood_C,
                              "p_G": likelihood_G,
                              "p_T": likelihood_T,
                          },
                          ...
                      ]
                  },
                  ...
              }

    Notes:
        - The function calculates partial likelihoods for each site in the alignment.
        - It constructs a dictionary to store partial likelihood information per site.
        - The partial likelihoods are calculated for both left and right descendants of each edge.
        - The resulting dictionary provides detailed information about likelihoods and branch lengths.
    Warning:
        Ensure that the provided tree structure and alignment data are accurate and compatible.
    """
    alignment_look_up_table = get_alignment_look_up_table(alignment)
    partial_likelihood_per_site_storage = {}

    for i in range(0, len(alignment[0].seq), 1):
        graph = convert_ete3_tree_to_directed_acyclic_graph(
            tree, alignment[:, (i + 1) - 1 : (i + 1)], alignment_look_up_table
        )

        for edge in graph.get_edges():
            right, left, branch_length = edge

            p1 = partial_likelihood(graph, left, right, rate_matrix)

            p2 = partial_likelihood(graph, right, left, rate_matrix)

            if (
                f"({left.name}, {right.name})"
                not in partial_likelihood_per_site_storage
            ):
                partial_likelihood_per_site_storage[f"({left.name}, {right.name})"] = {
                    "left": [
                        {
                            "Node": left.name,
                            "Site": i,
                            "branch_length": branch_length,
                            "p_A": p1[0],
                            "p_C": p1[1],
                            "p_G": p1[2],
                            "p_T": p1[3],
                        }
                    ],
                    "right": [
                        {
                            "Node": right.name,
                            "Site": i,
                            "branch_length": branch_length,
                            "p_A": p2[0],
                            "p_C": p2[1],
                            "p_G": p2[2],
                            "p_T": p2[3],
                        }
                    ],
                }
            else:
                partial_likelihood_per_site_storage[f"({left.name}, {right.name})"][
                    "left"
                ].append(
                    {
                        "Node": left.name,
                        "Site": i,
                        "branch_length": branch_length,
                        "p_A": p1[0],
                        "p_C": p1[1],
                        "p_G": p1[2],
                        "p_T": p1[3],
                    }
                )
                partial_likelihood_per_site_storage[f"({left.name}, {right.name})"][
                    "right"
                ].append(
                    {
                        "Node": right.name,
                        "Site": i,
                        "branch_length": branch_length,
                        "p_A": p2[0],
                        "p_C": p2[1],
                        "p_G": p2[2],
                        "p_T": p2[3],
                    }
                )
    return partial_likelihood_per_site_storage


@cache
def partial_likelihood(tree, node, coming_from, rate_matrix):
    """
    Compute the partial likelihood of a given node in a directed acyclic graph (DAG) tree.

    This function calculates the partial likelihood for a specific node in the DAG tree.
    The partial likelihood represents the conditional probability of observing the sequence
    data at a particular node, given the evolutionary history and rate matrix.

    Args:
        tree (DirectedAcyclicGraph): The directed acyclic graph representing the tree structure.
        node (Node): The current node for which the partial likelihood is being calculated.
        coming_from (Node): The previous node in the traversal path.
        rate_matrix (RateMatrix): The rate matrix representing nucleotide substitution rates.

    Returns:
        np.array: The calculated partial likelihood vector for the current node's state.

    Notes:
        - The function uses memorization with the @cache decorator to store and reuse computed results.
        - The DAG tree's nodes should have connections established using the "connect" method.
        - The coming_from parameter helps to avoid redundant calculations by skipping the reverse path.
    """
    results = 1
    # If the current node is a leaf, return its initial likelihood vector.
    if node.is_leaf():
        return node.state

    # Iterate through child nodes connected to the current node.
    for child in node.connected.keys():
        # Avoid traversing the path back to the parent node (coming_from).
        if not child.__eq__(coming_from):
            # Calculate the exponential matrix and partial likelihood for the child node.
            e = calculate_exponential_matrix(rate_matrix, node.connected[child])
            p = partial_likelihood(tree, child, node, rate_matrix)
            # Update the results by multiplying with the exponential matrix and partial likelihood.
            results = results * (e @ p)
    return results


@cache
def calculate_exponential_matrix(rate_matrix, branch_length):
    """
    Calculate the matrix exponential of a rate matrix scaled by a given branch length.

    This function computes the matrix exponential of a rate matrix, which represents
    the instantaneous rate of nucleotide substitutions between different states.
    The exponential matrix is calculated by exponentiation the rate matrix scaled by
    the given branch length, which represents the time interval of evolution along a branch.

    Args:
        rate_matrix (RateMatrix): The rate matrix representing nucleotide substitution rates.
        branch_length (float): The length of the evolutionary branch for which to calculate
                              the matrix exponential.

    Returns:
        np.array: The resulting matrix exponential of the scaled rate matrix.

    Notes:
        - The @cache decorator is used for memoization, storing previously computed results.
        - The resulting matrix describes the transition probabilities over the given time interval.
        - The matrix exponential plays a fundamental role in modeling evolutionary processes.

    Warning:
        Ensure that the rate matrix and branch length are appropriately defined for accurate results.
    """

    return expm(rate_matrix.rate_matrix * branch_length)


def test_one_partial_likelihood():
    a1, b2, u3, u4, c5, b6 = [
        Node("A", get_initial_likelihood_vector("A")),
        Node("B", get_initial_likelihood_vector("C")),
        Node("3"),
        Node("4"),
        Node("C", get_initial_likelihood_vector("A")),
        Node("D", get_initial_likelihood_vector("A")),
    ]

    tree = DirectedAcyclicGraph(
        [(a1, u3, 0.01), (u3, b2, 0.01), (u3, u4, 0.01), (u4, c5, 0.01), (u4, b6, 0.01)]
    )

    for edge in tree.get_edges():
        left, right, branch_length = edge
        p1 = partial_likelihood(tree, left, right)
        p2 = partial_likelihood(tree, right, left)
        test(p1, p2)


def single_rate_analysis(
    t,
    alignment,
    rate_matrix,
    state_frequencies,
    array_left_eigenvectors,
    array_right_eigenvectors,
    multiplicity,
):
    partial_likelihood_per_site_storage = calculate_partial_likelihoods_for_sites(
        t, alignment, rate_matrix
    )
    result_test_dictionary = {}
    result_list = []
    for edge in partial_likelihood_per_site_storage.keys():
        left_partial_likelihood = pd.DataFrame(
            partial_likelihood_per_site_storage[edge]["left"]
        )
        right_partial_likelihood = pd.DataFrame(
            partial_likelihood_per_site_storage[edge]["right"]
        )

        branch_type = "external"

        if "Node" in edge[1] and "Node" in edge[0]:
            branch_type = "internal"

        """
        (
            delta,
            c_s,
            c_sTwoSequence,
            p_value,
            result_test,
            result_test_tip2tip,
        ) = calculate_test_statistic(
            multiplicity,
            array_left_eigenvectors,
            state_frequencies,
            left_partial_likelihood,
            right_partial_likelihood,
            4,
            branch_type,
            alpha=0.05,
        )
        """

        (
            delta,
            c_s,
            c_sTwoSequence,
            p_value,
            result_test,
            result_test_tip2tip,
        ) = calculate_test_statistic_posterior_distribution(
            multiplicity,
            array_right_eigenvectors,
            state_frequencies,
            left_partial_likelihood,
            right_partial_likelihood,
            4,
            branch_type,
            alpha=0.05,
        )
        result_list.append(
            {
                "edge": edge,
                "delta": delta,
                "c_s": c_s,
                "c_sTwoSequence": c_sTwoSequence,
                "p_value": p_value,
                "result_test": result_test,
                "result_test_tip2tip": result_test_tip2tip,
                "branch_length": left_partial_likelihood["branch_length"][0],
            }
        )

    result_test_dictionary["single_rate"] = result_list
    return result_test_dictionary


def multiple_rate_analysis(
    t,
    category_rates_factors,
    rate_matrix,
    state_frequencies,
    array_left_eigenvectors,
    array_right_eigenvectors,
    multiplicity,
    per_rate_category_alignment,
):
    result_rate_dictionary = {}

    for rate, alignment in per_rate_category_alignment.items():
        relative_rate = category_rates_factors[rate]["Relative_rate"]
        result_list = []
        rescaled_tree = rescale_branch_lengths(t, relative_rate)
        partial_likelihood_per_site_storage = calculate_partial_likelihoods_for_sites(
            rescaled_tree, alignment, rate_matrix
        )

        for edge in partial_likelihood_per_site_storage.keys():
            left_partial_likelihood = pd.DataFrame(
                partial_likelihood_per_site_storage[edge]["left"]
            )
            right_partial_likelihood = pd.DataFrame(
                partial_likelihood_per_site_storage[edge]["right"]
            )
            branch_type = "external"

            if "Node" in edge[1] and "Node" in edge[0]:
                branch_type = "internal"

            """
            (
                delta,
                c_s,
                c_sTwoSequence,
                p_value,
                result_test,
                result_test_tip2tip,
            ) = calculate_test_statistic(
                multiplicity,
                array_left_eigenvectors,
                state_frequencies,
                left_partial_likelihood,
                right_partial_likelihood,
                4,
                branch_type,
                alpha=0.05,
            )
            """
            (
                delta,
                c_s,
                c_sTwoSequence,
                p_value,
                result_test,
                result_test_tip2tip,
            ) = calculate_test_statistic_posterior_distribution(
                multiplicity,
                array_right_eigenvectors,
                state_frequencies,
                left_partial_likelihood,
                right_partial_likelihood,
                4,
                branch_type,
                alpha=0.05,
            )

            result_list.append(
                {
                    "edge": edge,
                    "delta": delta,
                    "c_s": c_s,
                    "c_sTwoSequence": c_sTwoSequence,
                    "p_value": p_value,
                    "result_test": result_test,
                    "result_test_tip2tip": result_test_tip2tip,
                    "category_rate": rate,
                    "branch_length": left_partial_likelihood["branch_length"][0],
                }
            )

        result_rate_dictionary[rate] = result_list
    return result_rate_dictionary


def main_p():
    alignment_file = "./Clemens/example_3/sim-JC+G-AC1-AG1-AT1-CG1-CT1-GT1-alpha1.2-taxa64-len1000bp-bla0.01-blb0.8-blc0.2-rep01.fasta"
    newick_string = "(((S1:1,S2:2),(S3:1,S4:1),S5:1))"
    t = Tree(newick_string, format=1)
    t = name_nodes_by_level_order(t)
    alignment = read_alignment_file(alignment_file)
    rate_matrix = RateMatrix(RATE_MATRIX)
    partial_likelihood_per_site_storage = calculate_partial_likelihoods_for_sites(
        t, alignment, rate_matrix
    )


def main():
    alignment_file = "./Clemens/example_3/sim-JC+G-AC1-AG1-AT1-CG1-CT1-GT1-alpha1.2-taxa64-len1000bp-bla0.01-blb0.8-blc0.2-rep01.fasta"
    alignment = read_alignment_file(alignment_file)
    t = parse_newick_file(
        "./Clemens/example_3/sim-JC+G-AC1-AG1-AT1-CG1-CT1-GT1-alpha1.2-taxa64-len1000bp-bla0.01-blb0.8-blc0.2-rep01.fasta.treefile"
    )
    t = name_nodes_by_level_order(t)
    number_rate = 4

    RATE_MATRIX, psi_matrix = parse_rate_matrices_from_file(
        4,
        "./Clemens/example_3/sim-JC+G-AC1-AG1-AT1-CG1-CT1-GT1-alpha1.2-taxa64-len1000bp-bla0.01-blb0.8-blc0.2-rep01.fasta",
    )
    rate_matrix = RateMatrix(RATE_MATRIX)

    (
        array_left_eigenvectors,
        array_right_eigenvectors,
        multiplicity,
    ) = spectral_decomposition(RATE_MATRIX, psi_matrix)

    if number_rate == 1:
        single_rate_analysis(
            t, alignment, rate_matrix, array_left_eigenvectors, multiplicity
        )
    else:
        site_probability = parse_file_to_data_frame(
            "./Clemens/example_3/sim-JC+G-AC1-AG1-AT1-CG1-CT1-GT1-alpha1.2-taxa64-len1000bp-bla0.01-blb0.8-blc0.2-rep01.fasta.siteprob"
        )
        per_rate_category_alignment = split_msa_into_rate_categories_in_place(
            site_probability, alignment
        )
        category_rates_factors = parse_category_rates(
            "./Clemens/example_3/sim-JC+G-AC1-AG1-AT1-CG1-CT1-GT1-alpha1.2-taxa64-len1000bp-bla0.01-blb0.8-blc0.2-rep01.fasta.iqtree"
        )
        multiple_rate_analysis(
            t,
            category_rates_factors,
            rate_matrix,
            array_left_eigenvectors,
            array_right_eigenvectors,
            multiplicity,
            per_rate_category_alignment,
        )


from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def test_one():
    # Create SeqRecord objects for your sequences
    seq1 = SeqRecord(Seq("AGTATA"), id="A")
    seq2 = SeqRecord(Seq("CGTATG"), id="B")
    seq3 = SeqRecord(Seq("GGTATG"), id="C")
    seq4 = SeqRecord(Seq("GGTACG"), id="D")

    # Create a MultipleSeqAlignment object
    alignment = MultipleSeqAlignment([seq1, seq2, seq3, seq4])
    newick_string = "((A:0.2, B:0.4):0.3, (C:0.5,D:0.2):2);"
    t = Tree(newick_string, format=1)
    t = name_nodes_by_level_order(t)

    print(t.get_ascii(attributes=["name", "dist"]))

    rate_matrix = RateMatrix(RATE_MATRIX)

    partial_likelihood_per_site_storage = calculate_partial_likelihoods_for_sites(
        t, alignment, rate_matrix
    )

    # for key, values in partial_likelihood_per_site_storage.items():
    #    print(key)
    #    print("#### Left")
    #    for value in values["left"]:
    #        print(value)
    #        print("\n")


#
#
#    print("#### Right")
#    for value in values["right"]:
#        print(value)
#        print("\n")


if __name__ == "__main__":
    test_one()
