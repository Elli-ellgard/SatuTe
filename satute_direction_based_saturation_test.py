import numpy as np
import pandas as pd
from functools import cache
from scipy.sparse.linalg import expm
from ete3 import Tree
from satute_rate_categories_and_alignments import read_alignment_file
from satute_trees_and_subtrees import rescale_branch_lengths, parse_newick_file
from satute_test_statistic_using_posterior_distribution import (
    calculate_test_statistic_posterior_distribution,
)
from satute_repository import parse_rate_matrices_from_file
from satute_util_new import spectral_decomposition
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from concurrent.futures import ProcessPoolExecutor
from graph import Graph, Node
from rate_matrix import RateMatrix


NUCLEOTIDE_CODE_VECTOR = {
    "A": np.array([1, 0, 0, 0]),
    "C": np.array([0, 1, 0, 0]),
    "G": np.array([0, 0, 1, 0]),
    "T": np.array([0, 0, 0, 1]),
    "U": np.array([0, 0, 0, 1]),
    "R": np.array([1, 0, 1, 0]),
    "Y": np.array([0, 1, 0, 1]),
    "K": np.array([0, 0, 1, 1]),
    "M": np.array([1, 1, 0, 0]),
    "S": np.array([0, 1, 1, 0]),
    "W": np.array([1, 0, 0, 1]),
    "B": np.array([0, 1, 1, 1]),
    "D": np.array([1, 0, 1, 1]),
    "H": np.array([1, 1, 0, 1]),
    "V": np.array([1, 1, 1, 0]),
    #The following keys are treated as State_Unknown in IQ-Tree
    "N": np.array([1, 1, 1, 1]),
    "-": np.array([1, 1, 1, 1]),
    "?": np.array([1, 1, 1, 1]),
    ".": np.array([1, 1, 1, 1]),
    "~": np.array([1, 1, 1, 1]),
    "X": np.array([1, 1, 1, 1]),
    #Additional key from EvoNaps database
    "!": np.array([1, 1, 1, 1]),
}

RATE_MATRIX = np.array([[-3, 1, 1, 1], [1, -3, 1, 1], [1, 1, -3, 1], [1, 1, 1, -3]])


@cache
def partial_likelihood(tree, node, coming_from, rate_matrix, factor=0):
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
        return node.state, factor

    # Iterate through child nodes connected to the current node.
    for child in node.connected.keys():
        # Avoid traversing the path back to the parent node (coming_from).
        if not child == coming_from:
            # Calculate the exponential matrix and partial likelihood for the child node.
            e = calculate_exponential_matrix(rate_matrix, node.connected[child])
            p, factor = partial_likelihood(tree, child, node, rate_matrix, factor)
            # Update the results by multiplying with the exponential matrix and partial likelihood.
            results = results * (e @ p)

            # Check if scaling is needed
            if results.sum() < 1e-50:
                results *= 1e50
                factor += 50  # Keeping track of the total scaling

    return results, factor


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


def get_initial_likelihood_vector(state):
    """Get the initial likelihood vector for a given nucleotide state."""
    return np.array(NUCLEOTIDE_CODE_VECTOR[state]).T


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
    return Graph(edge_list)


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


def filter_graph_edges_by_focus(graph, focused_edge):
    """Filter the graph edges based on the focused edge.

    Args:
        graph: The directed acyclic graph.
        focused_edge (tuple): A tuple containing the names of the nodes representing the focused edge.

    Raises:
        ValueError: If the focused_edge does not exist in the graph.
    """
    # Check if the focused_edge exists in the graph's edges
    focused_edge_found = any(
        (node1.name in focused_edge and node2.name in focused_edge)
        for node1, node2, _ in graph.get_edges()
    )

    if not focused_edge_found:
        raise ValueError(f"Focused edge {focused_edge} not found in the graph!")

    # Filter the graph's edges to only include the focused_edge
    filtered_edges = [
        edge
        for edge in graph.get_edges()
        if edge[0].name in focused_edge and edge[1].name in focused_edge
    ]

    graph.set_edges(filtered_edges)


def get_partial_likelihood_dict(edge, p1, p2, i, branch_length):
    """Construct the partial likelihood dictionary for a specific edge and site."""
    left, right = edge
    return {
        "left": {
            "Node": left.name,
            "Site": i,
            "branch_length": branch_length,
            "p_A": p1[0],
            "p_C": p1[1],
            "p_G": p1[2],
            "p_T": p1[3],
        },
        "right": {
            "Node": right.name,
            "Site": i,
            "branch_length": branch_length,
            "p_A": p2[0],
            "p_C": p2[1],
            "p_G": p2[2],
            "p_T": p2[3],
        },
    }


def update_partial_likelihood_storage(storage, edge_name, likelihood_data):
    """Update the storage with new partial likelihood data."""
    if edge_name not in storage:
        storage[edge_name] = {"left": [], "right": []}
    storage[edge_name]["left"].append(likelihood_data["left"])
    storage[edge_name]["right"].append(likelihood_data["right"])


def calculate_partial_likelihoods_for_sites(
    tree, alignment, rate_matrix, focused_edge=None
):
    """... [Same docstring as before] ..."""
    alignment_look_up_table = get_alignment_look_up_table(alignment)
    partial_likelihood_per_site_storage = {}

    for i in range(0, len(alignment[0].seq), 1):
        graph = convert_ete3_tree_to_directed_acyclic_graph(
            tree, alignment[:, (i + 1) - 1 : (i + 1)], alignment_look_up_table
        )

        if focused_edge:
            filter_graph_edges_by_focus(graph, focused_edge)

        for edge in graph.get_edges():
            right, left, branch_length = edge

            p1, p1_factor = partial_likelihood(graph, left, right, rate_matrix)
            p2, p2_factor = partial_likelihood(graph, right, left, rate_matrix)

            likelihood_data = get_partial_likelihood_dict(
                (left, right), p1, p2, i, branch_length
            )
            edge_name = f"({left.name}, {right.name})"
            update_partial_likelihood_storage(
                partial_likelihood_per_site_storage, edge_name, likelihood_data
            )

    return partial_likelihood_per_site_storage


def store_test_results(edge, rate, left_partial_likelihood, results):
    """
    Store the test results in a structured dictionary format.

    Args:
        edge (str): The edge identifier.
        rate (float): The rate of the category.
        left_partial_likelihood (dict): Contains partial likelihoods.
        results (dict): Results from the test statistic calculation.

    Returns:
        dict: Structured result entry.
    """
    return {
        "edge": edge,
        "delta": results.get("delta"),
        "c_s": results.get("c_s"),
        "c_sTwoSequence": results.get("c_s_two_sequence"),
        "p_value": results.get("p_value"),
        "result_test": results.get("result_test"),
        "result_test_tip2tip": results.get("result_test_tip2tip"),
        "category_rate": rate,
        "branch_length": left_partial_likelihood.get("branch_length", [None])[0],
    }


def multiple_rate_analysis(
    initial_tree,
    category_rates_factors,
    rate_matrix,
    state_frequencies,
    array_left_eigenvectors,
    array_right_eigenvectors,
    multiplicity,
    per_rate_category_alignment,
    alpha=0.05,
    focused_edge=None,
):
    result_rate_dictionary = {}

    for rate, alignment in per_rate_category_alignment.items():
        relative_rate = category_rates_factors[rate]["Relative_rate"]
        result_list = []

        rescaled_copied_tree = initial_tree.copy("deepcopy")
        rescale_branch_lengths(rescaled_copied_tree, relative_rate)

        partial_likelihood_per_site_storage = calculate_partial_likelihoods_for_sites(
            rescaled_copied_tree, alignment, rate_matrix, focused_edge
        )

        for edge, likelihoods in partial_likelihood_per_site_storage.items():
            left_partial_likelihood = pd.DataFrame(likelihoods["left"])
            right_partial_likelihood = pd.DataFrame(likelihoods["right"])

            branch_type = determine_branch_type(edge)

            results = process_test_statistics_posterior(
                multiplicity,
                array_right_eigenvectors,
                state_frequencies,
                left_partial_likelihood,
                right_partial_likelihood,
                branch_type,
                alpha,
            )

            result_list.append(
                store_test_results(edge, rate, left_partial_likelihood, results)
            )

        result_rate_dictionary[rate] = {
            "result_list": result_list,
            "rescaled_tree": rescaled_copied_tree,
        }

    return result_rate_dictionary


def single_rate_analysis(
    initial_tree,
    alignment,
    rate_matrix,
    state_frequencies,
    array_left_eigenvectors,
    array_right_eigenvectors,
    multiplicity,
    alpha=0.05,
    focused_edge=None,
):
    """
    Performs a single rate analysis on a phylogenetic tree.

    Parameters:
    - initial_tree: The initial phylogenetic tree.
    - alignment: The alignment data.
    - rate_matrix: The matrix of rates.
    - state_frequencies: Frequencies of different states.
    - array_left_eigenvectors: Left eigenvectors array.
    - array_right_eigenvectors: Right eigenvectors array.
    - multiplicity: The multiplicity value.
    - alpha: The significance level for tests. Default is 0.05.
    - focused_edge: Specific edge to focus on (if any).

    Returns:
    A dictionary containing the results of the test for each edge.
    """

    # Calculate the partial likelihoods for all sites.
    partial_likelihood_per_site_storage = calculate_partial_likelihoods_for_sites(
        initial_tree, alignment, rate_matrix, focused_edge
    )

    # Dictionary to store the results.
    result_test_dictionary = {}

    # List to accumulate the results for each edge.
    result_list = []

    # Iterate over each edge and the associated likelihoods.
    for edge, likelihoods in partial_likelihood_per_site_storage.items():
        # Convert the left and right likelihoods to dataframes for easier processing.
        left_partial_likelihood = pd.DataFrame(likelihoods["left"])
        right_partial_likelihood = pd.DataFrame(likelihoods["right"])

        # Determine the type of branch (internal or external).
        branch_type = determine_branch_type(edge)

        # Calculate the test statistics using the posterior distribution.
        results = process_test_statistics_posterior(
            multiplicity,
            array_right_eigenvectors,
            state_frequencies,
            left_partial_likelihood,
            right_partial_likelihood,
            branch_type,
            alpha,
        )

        # Store the results of the test for the given edge.
        result_list.append(
            store_test_results(edge, "single_rate", left_partial_likelihood, results)
        )

    # Add the accumulated results to the main dictionary.
    result_test_dictionary["single_rate"] = result_list

    return result_test_dictionary


def determine_branch_type(edge):
    """
    Determine the type of the branch based on the edge.

    Args:
        edge (tuple): Edge details.

    Returns:
        str: 'internal' if the edge is internal, 'external' otherwise.
    """
    return "internal" if "Node" in edge[1] and "Node" in edge[0] else "external"


def process_test_statistics_posterior(
    multiplicity,
    array_right_eigenvectors,
    state_frequencies,
    left_partial_likelihood,
    right_partial_likelihood,
    branch_type,
    alpha,
):
    """
    Calculate test statistics for posterior distribution and store results in a dictionary.

    Args:
        ... [all the function arguments with their descriptions]

    Returns:
        dict: Dictionary containing calculated results.
    """
    results = calculate_test_statistic_posterior_distribution(
        multiplicity,
        array_right_eigenvectors,
        state_frequencies,
        left_partial_likelihood,
        right_partial_likelihood,
        4,
        branch_type,
        alpha,
    )

    result_keys = [
        "delta",
        "c_s",
        "c_s_two_sequence",
        "p_value",
        "result_test",
        "result_test_tip2tip",
    ]

    return {key: value for key, value in zip(result_keys, results)}


"""
    ======= Tests =======
"""


def name_nodes_by_level_order(tree):
    """Name nodes in a tree based on level-order traversal."""
    i = 1
    for node in tree.traverse("levelorder"):
        if not node.is_leaf():
            node.name = f"Node{i}"
            i += 1
    return tree


def test_one_partial_likelihood():
    a1, b2, u3, u4, c5, b6 = [
        Node("A", get_initial_likelihood_vector("A")),
        Node("B", get_initial_likelihood_vector("C")),
        Node("3"),
        Node("4"),
        Node("C", get_initial_likelihood_vector("A")),
        Node("D", get_initial_likelihood_vector("A")),
    ]

    tree = Graph(
        [(a1, u3, 0.01), (u3, b2, 0.01), (u3, u4, 0.01), (u4, c5, 0.01), (u4, b6, 0.01)]
    )

    for edge in tree.get_edges():
        left, right, branch_length = edge
        p1 = partial_likelihood(tree, left, right)
        p2 = partial_likelihood(tree, right, left)


def test_three():
    alignment_file = "./Clemens/example_3/sim-JC+G-AC1-AG1-AT1-CG1-CT1-GT1-alpha1.2-taxa64-len1000bp-bla0.01-blb0.8-blc0.2-rep01.fasta"
    newick_string = "(((S1:1,S2:2),(S3:1,S4:1),S5:1))"
    t = Tree(newick_string, format=1)
    t = name_nodes_by_level_order(t)
    alignment = read_alignment_file(alignment_file)
    rate_matrix = RateMatrix(RATE_MATRIX)
    partial_likelihood_per_site_storage = calculate_partial_likelihoods_for_sites(
        t, alignment, rate_matrix
    )


def test_two():
    alignment_file = "./Clemens/example_3/sim-JC+G-AC1-AG1-AT1-CG1-CT1-GT1-alpha1.2-taxa64-len1000bp-bla0.01-blb0.8-blc0.2-rep01.fasta"
    alignment = read_alignment_file(alignment_file)
    t = parse_newick_file(
        "./Clemens/example_3/sim-JC+G-AC1-AG1-AT1-CG1-CT1-GT1-alpha1.2-taxa64-len1000bp-bla0.01-blb0.8-blc0.2-rep01.fasta.treefile"
    )
    t = name_nodes_by_level_order(t)

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


if __name__ == "__main__":
    test_one()
