from functools import cache
from scipy.sparse.linalg import expm
from ete3 import Tree
from rate_matrix import RateMatrix
from Bio.Align import MultipleSeqAlignment
from graph import (
    filter_graph_edges_by_focus,
    get_alignment_look_up_table,
    convert_tree_to_state_graph,
    Node,
)
from amino_acid_models import AMINO_ACIDS


@cache
def partial_likelihood(node, coming_from: Node, rate_matrix: RateMatrix, factor=0):
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
        if child != coming_from:
            # Calculate the exponential matrix and partial likelihood for the child node.
            e = calculate_exponential_matrix(rate_matrix, node.connected[child])
            p, factor = partial_likelihood(child, node, rate_matrix, factor)
            # Update the results by multiplying with the exponential matrix and partial likelihood.
            results = results * (e @ p)

            # Check if scaling is needed
            if results.sum() < 1e-50:
                results *= 1e50
                factor += 50  # Keeping track of the total scaling
    return results, factor


@cache
def calculate_exponential_matrix(rate_matrix: RateMatrix, branch_length: float):
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


def get_partial_likelihood_dict(edge, p1, p2, i, branch_length, states):
    """Construct the partial likelihood dictionary for a specific edge and site with dynamic states."""
    left, right = edge

    # Helper function to create the likelihood dictionary for a node
    def create_likelihood_dict(node_name, probabilities):
        return {
            "Node": node_name,
            "Site": i,
            "branch_length": branch_length,
            **{f"p_{state}": prob for state, prob in zip(states, probabilities)},
        }

    return {
        "left": create_likelihood_dict(left.name, p1),
        "right": create_likelihood_dict(right.name, p2),
    }


def update_partial_likelihood_storage(storage: dict, edge_name, likelihood_data):
    """Update the storage with new partial likelihood data."""
    if edge_name not in storage:
        storage[edge_name] = {
            "left": {
                "likelihoods": [],
            },
            "right": {
                "likelihoods": [],
            },
        }
    storage[edge_name]["left"]["likelihoods"].append(likelihood_data["left"])
    storage[edge_name]["right"]["likelihoods"].append(likelihood_data["right"])


def calculate_partial_likelihoods_for_sites(
    tree: Tree,
    alignment: MultipleSeqAlignment,
    rate_matrix: RateMatrix,
    focused_edge=None,
):
    """Calculates the partial likelihoods of each site in a multiple sequence alignment across a given phylogeny.

    Args:
        tree: The phylogenetic tree represented as an ete3 Tree object.
        alignment: The multiple sequence alignment represented as a NumPy array.
        rate_matrix: The rate matrix governing nucleotide substitution probabilities.
        focused_edge: (Optional) A tuple indicating the specific edge to focus on for partial likelihood calculations.

    Returns:
        dict: A dictionary containing partial likelihood data for each edge and site.
    """

    state_type = "nucleotide"
    states = ["A", "C", "G", "T"]

    if rate_matrix.rate_matrix.shape != (4, 4):
        state_type = "amino_acid"
        states = AMINO_ACIDS

    # Create a lookup table for efficient sequence indexing
    alignment_look_up_table = get_alignment_look_up_table(alignment)

    # Initialize a storage dictionary for partial likelihood data
    partial_likelihood_per_site_storage = {}

    # Iterate over each site in the alignment
    for site_idx in range(0, alignment.get_alignment_length(), 1):
        # Convert the alignment column for the current site into a NumPy array
        msa_column = alignment[:, (site_idx + 1) - 1 : (site_idx + 1)]

        # Convert the ETE3 tree to a directed acyclic graph (DAG) for the current site
        graph = convert_tree_to_state_graph(
            tree=tree,
            msa_column=msa_column,
            alignment_look_up_table=alignment_look_up_table,
            state_type=state_type,
        )

        # If a specific edge is focused on, filter the graph to only include that edge
        if focused_edge:
            filter_graph_edges_by_focus(graph, focused_edge)

        # Calculate the partial likelihoods for each edge in the graph
        for edge in graph.get_edges():
            right, left, length = edge
            p1, p1_factor = partial_likelihood(left, right, rate_matrix)
            p2, p2_factor = partial_likelihood(right, left, rate_matrix)

            likelihood_data = get_partial_likelihood_dict(
                (left, right), p1, p2, site_idx, length, states
            )

            # Store partial likelihood data in the dictionary
            edge_name = f"({left.name}, {right.name})"

            update_partial_likelihood_storage(
                partial_likelihood_per_site_storage, edge_name, likelihood_data
            )

    return partial_likelihood_per_site_storage
