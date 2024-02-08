import pandas as pd
from pandas import DataFrame
from partial_likelihood import calculate_partial_likelihoods_for_sites
from graph import calculate_subtree_edge_metrics
from satute_sequences import dict_to_alignment
from ete3 import Tree
from Bio.Align import MultipleSeqAlignment
from rate_matrix import RateMatrix
from satute_statistic_posterior_distribution import (
    calculate_test_statistic_posterior_distribution,
)
from satute_trees import (
    rescale_branch_lengths,
    collapse_identical_leaf_sequences,
)


def single_rate_analysis(
    initial_tree: Tree,
    alignment: MultipleSeqAlignment,
    rate_matrix: RateMatrix,
    state_frequencies: list,
    array_right_eigenvectors: list,
    multiplicity: int,
    alpha: float = 0.05,
    focused_edge=None,
):
    """
    Performs a single rate analysis on a phylogenetic tree.

    Parameters:
    - initial_tree: The initial phylogenetic tree.
    - alignment: The alignment data.
    - rate_matrix: The matrix of rates.
    - state_frequencies: Frequencies of different states.
    - array_right_eigenvectors: Right eigenvectors array.
    - multiplicity: The multiplicity value.
    - alpha: The significance level for tests. Default is 0.05.
    - focused_edge: Specific edge to focus on (if any).

    Returns:
    A dictionary containing the results of the test for each edge.
    """
    # Calculate partial likelihoods for all sites
    partial_likelihood_per_site_storage = calculate_partial_likelihoods_for_sites(
        initial_tree, alignment, rate_matrix, focused_edge
    )

    # Count leaves and branches for subtrees
    edge_subtree_count_dict = calculate_subtree_edge_metrics(initial_tree, focused_edge)

    # Initialize a dictionary and a list to store results
    result_test_dictionary = {}
    result_list = []

    # Iterate over each edge and process likelihoods
    for edge, likelihoods in partial_likelihood_per_site_storage.items():
        # Convert left and right likelihoods to data frames
        left_partial_likelihood = pd.DataFrame(likelihoods["left"]["likelihoods"])
        right_partial_likelihood = pd.DataFrame(likelihoods["right"]["likelihoods"])

        # Count the number of leaves in the left and right subtree
        number_leaves_left_subtree = edge_subtree_count_dict[edge]["left"][
            "leave_count"
        ]
        number_leaves_right_subtree = edge_subtree_count_dict[edge]["right"][
            "leave_count"
        ]

        _type = edge_subtree_count_dict[edge]["type"]

        # Calculate test statistics using the posterior distribution
        results = process_test_statistics_posterior(
            multiplicity,
            array_right_eigenvectors,
            state_frequencies,
            left_partial_likelihood,
            right_partial_likelihood,
            number_leaves_left_subtree,
            number_leaves_right_subtree,
            _type,
            alpha,
        )

        # Store the results for the given edge
        result_list.append(
            store_test_results(edge, "single_rate", left_partial_likelihood, results)
        )

    # Add all results to the main result dictionary
    result_test_dictionary["single_rate"] = {
        "result_list": result_list,
        "rescaled_tree": initial_tree,
        "partial_likelihoods": partial_likelihood_per_site_storage,
    }

    return result_test_dictionary


def single_rate_analysis_collapsed_tree(
    initial_tree: Tree,
    alignment: MultipleSeqAlignment,
    rate_matrix: RateMatrix,
    state_frequencies: list,
    array_right_eigenvectors: list,
    multiplicity: int,
    alpha: float = 0.05,
    focused_edge=None,
):
    # Step 1: Map sequence IDs to their sequences from the alignment for easy access
    sequence_dict = {record.id: str(record.seq) for record in alignment}

    # Step 2: Create a deep copy of the initial tree to modify without altering the original
    collapsed_tree_one = initial_tree.copy("deepcopy")

    # Step 3: Collapse nodes in the tree with identical sequences to simplify the tree structure
    sequence_dict, collapsed_nodes = collapse_identical_leaf_sequences(
        collapsed_tree_one, sequence_dict
    )

    # Step 4: Convert the sequence dictionary back to a MultipleSeqAlignment object after collapsing nodes
    alignment = dict_to_alignment(sequence_dict)

    # Step 5: Carry out single rate analysis using the collapsed tree and updated alignment
    results = single_rate_analysis(
        collapsed_tree_one,
        alignment,
        rate_matrix,
        state_frequencies,
        array_right_eigenvectors,
        multiplicity,
        alpha,
        focused_edge,
    )

    # Step 6: Augment the results with additional data for the 'twin' nodes
    insert_nan_values_for_identical_taxa(
        results["single_rate"]["result_list"],
        collapsed_nodes,
        initial_tree,
        "single_rate",
        focused_edge,
    )

    # Step 7: Compile all the results into a dictionary and return it
    return results


def process_test_statistics_posterior(
    multiplicity: int,
    array_right_eigenvectors,
    state_frequencies: list,
    left_partial_likelihood: DataFrame,
    right_partial_likelihood: DataFrame,
    number_leaves_left_subtree: int,
    number_leaves_right_subtree: int,
    branch_type: str,
    alpha: type,
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
        number_leaves_left_subtree,
        number_leaves_right_subtree,
        branch_type,
        alpha,
    )

    result_keys = [
        "test_statistic",
        "delta",
        "variance",
        "p_value",
        "decision_test",
        "decision_corrected_test_tips",
        "decision_corrected_test_branches",
        "decision_test_tip2tip",
    ]

    return {key: value for key, value in zip(result_keys, results)}


def multiple_rate_analysis(
    initial_tree: Tree,
    category_rates_factors,
    rate_matrix: RateMatrix,
    state_frequencies: list,
    array_right_eigenvectors: list,
    multiplicity: int,
    per_rate_category_alignment: dict[str, MultipleSeqAlignment],
    alpha: float = 0.05,
    focused_edge=None,
):
    # Initialize a dictionary to store results for each rate category
    result_rate_dictionary = {}

    # Iterate over each rate category and its corresponding alignment
    for rate, sub_alignment in per_rate_category_alignment.items():
        if sub_alignment.get_alignment_length() == 0:
            continue

        # Step 1: Map sequence IDs to their sequences from the alignment for easy access
        sequence_dict = {record.id: str(record.seq) for record in sub_alignment}

        # Step 2: Retrieve the relative rate for the current category
        relative_rate = category_rates_factors[rate]["Relative_rate"]

        # Initialize a list to store results for the current rate category
        result_list = []

        # Step 3: Create a deep copy of the initial tree and collapse nodes with identical sequences
        collapsed_rescaled_tree_one = initial_tree.copy("deepcopy")

        sequence_dict, collapsed_nodes = collapse_identical_leaf_sequences(
            collapsed_rescaled_tree_one, sequence_dict
        )

        # Step 4: Rescale branch lengths according to the relative rate
        rescale_branch_lengths(collapsed_rescaled_tree_one, relative_rate)

        # Step 5: Convert the sequence dictionary back to a MultipleSeqAlignment object after collapsing nodes
        sub_alignment = dict_to_alignment(sequence_dict)

        if sub_alignment.get_alignment_length() != 0:
            # Step 6: Calculate partial likelihoods for all sites in the rescaled tree
            partial_likelihood_per_site_storage = (
                calculate_partial_likelihoods_for_sites(
                    collapsed_rescaled_tree_one,
                    sub_alignment,
                    rate_matrix,
                    focused_edge,
                )
            )

            # Step 7: Count leaves and branches for subtrees in the rescaled tree
            edge_subtree_metrics = calculate_subtree_edge_metrics(
                initial_tree, focused_edge
            )

            # Step 8: Process each edge and its associated likelihoods
            for edge, likelihoods in partial_likelihood_per_site_storage.items():

                # Convert left and right likelihoods to data frames
                left_partial_likelihood = pd.DataFrame(
                    likelihoods["left"]["likelihoods"]
                )

                right_partial_likelihood = pd.DataFrame(
                    likelihoods["right"]["likelihoods"]
                )

                # Get the number of leaves in the left and right subtree
                number_leaves_left_subtree = edge_subtree_metrics[edge]["left"][
                    "leave_count"
                ]

                number_leaves_right_subtree = edge_subtree_metrics[edge]["right"][
                    "leave_count"
                ]

                # Calculate test statistics using the posterior distribution
                results = process_test_statistics_posterior(
                    multiplicity,
                    array_right_eigenvectors,
                    state_frequencies,
                    left_partial_likelihood,
                    right_partial_likelihood,
                    number_leaves_left_subtree,
                    number_leaves_right_subtree,
                    edge_subtree_metrics[edge]["type"],
                    alpha,
                )

                # Store the results of the test for the given edge
                result_list.append(
                    store_test_results(edge, rate, left_partial_likelihood, results)
                )

            insert_nan_values_for_identical_taxa(
                result_list, collapsed_nodes, initial_tree, rate, focused_edge
            )

            # Step 10: Add the results for the current rate category to the main dictionary
            result_rate_dictionary[rate] = {
                "result_list": result_list,
                "rescaled_tree": collapsed_rescaled_tree_one,
                "partial_likelihoods": partial_likelihood_per_site_storage,
            }

    return result_rate_dictionary


def store_test_results(
    edge, rate: str, left_partial_likelihood: DataFrame, results: dict
):
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
        "test_statistic": results.get("test_statistic"),
        "edge": edge,
        "delta": results.get("delta"),
        "variance": results.get("variance"),
        "p_value": results.get("p_value"),
        "decision_test": results.get("decision_test"),
        "decision_corrected_test_tips": results.get("decision_corrected_test_tips"),
        "decision_corrected_test_branches": results.get(
            "decision_corrected_test_branches"
        ),
        "decision_test_tip2tip": results.get("decision_test_tip2tip"),
        "category_rate": rate,
        "branch_length": left_partial_likelihood.get("branch_length")[0],
    }


def insert_nan_values_for_identical_taxa(
    result_list: list,
    collapsed_nodes: dict,
    initial_tree: Tree,
    rate: str,
    focused_edge=None,
):
    """
    Insert NaN values into the result list for edges consisting of identical taxa.

    This function updates a result list by appending NaN values for specific metrics associated
    with edges between identical taxa. The function can focus on a specific edge if provided.

    Args:
        result_list (list): The list to which the NaN values are to be appended.
        collapsed_nodes (dict): A dictionary of parent nodes mapping to their child nodes.
        initial_tree (Tree): The initial phylogenetic tree used for distance calculation.
        rate (str): The category rate associated with the edge.
        focused_edge: A specific edge to focus on. If None,
            the function processes all edges in `collapsed_nodes`.

    """
    nan_values_dict = {
        "test_statistic": "Nan",
        "delta": "Nan",
        "variance": "Nan",
        "p_value": "Nan",
        "decision_test": "Nan",
        "decision_corrected_test_tips": "Nan",
        "decision_corrected_test_branches": "Nan",
        "decision_test_tip2tip": "Nan",
        "category_rate": rate,
    }

    for parent, children in collapsed_nodes.items():
        for child in children:
            if not focused_edge or (parent in focused_edge and child in focused_edge):
                edge_info = {
                    "edge": f"({parent},{child})",
                    "branch_length": initial_tree.get_distance(parent, child),
                    **nan_values_dict,
                }
                result_list.append(edge_info)
