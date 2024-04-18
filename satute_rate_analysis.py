import pandas as pd
from pandas import DataFrame
from partial_likelihood import calculate_partial_likelihoods_for_sites
from graph import calculate_subtree_edge_metrics
from satute_sequences import dict_to_alignment
from ete3 import Tree
from Bio.Align import MultipleSeqAlignment
from typing import List, Dict
from rate_matrix import RateMatrix
from satute_result import (
    TestResultsBranches,
    TestResultBranch,
    TestStatisticComponentsContainer,
)
from satute_statistic_posterior_distribution_components import (
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
    state_frequencies: List[float],
    array_right_eigenvectors: list,
    multiplicity: int,
    alpha: float = 0.05,
    focused_edge: str = None,
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
    results = TestResultsBranches()
    components_container = TestStatisticComponentsContainer()

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

        dimension = len(state_frequencies)

        components, result = calculate_test_statistic_posterior_distribution(
            multiplicity,
            array_right_eigenvectors,
            state_frequencies,
            left_partial_likelihood,
            right_partial_likelihood,
            dimension,
            number_leaves_left_subtree,
            number_leaves_right_subtree,
            _type,
            alpha,
        )

        result.add_result(
            "branch_length", left_partial_likelihood.get("branch_length")[0]
        )

        components_container.add_component("edge", components)
        # Store the results for the given edge
        results.add_branch(edge, result)

    # Add all results to the main result dictionary
    result_test_dictionary["single_rate"] = {
        "result_list": results,
        "rescaled_tree": initial_tree,
        "partial_likelihoods": partial_likelihood_per_site_storage,
        "components": components_container,
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
    focused_edge: str = None,
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
        focused_edge,
    )

    # Step 7: Compile all the results into a dictionary and return it
    return results


def multiple_rate_analysis(
    initial_tree: Tree,
    category_rates_factors,
    rate_matrix: RateMatrix,
    state_frequencies: List[float],
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
        results_list = TestResultsBranches()
        components_container = TestStatisticComponentsContainer()

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
                collapsed_rescaled_tree_one, focused_edge
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

                dimension = len(state_frequencies)

                test_components, test_result = (
                    calculate_test_statistic_posterior_distribution(
                        multiplicity,
                        array_right_eigenvectors,
                        state_frequencies,
                        left_partial_likelihood,
                        right_partial_likelihood,
                        dimension,
                        number_leaves_left_subtree,
                        number_leaves_right_subtree,
                        edge_subtree_metrics[edge]["type"],
                        alpha,
                    )
                )

                components_container.add_component(edge, test_components)

                # Store the results of the test for the given edge
                results_list.add_branch(edge, test_result)

                test_result.add_result(
                    "branch_length", left_partial_likelihood.get("branch_length")[0]
                )

            insert_nan_values_for_identical_taxa(
                results_list, collapsed_nodes, initial_tree, focused_edge
            )

            # Step 10: Add the results for the current rate category to the main dictionary
            result_rate_dictionary[rate] = {
                "result_list": results_list,
                "rescaled_tree": collapsed_rescaled_tree_one,
                "components": components_container,
                "partial_likelihoods": partial_likelihood_per_site_storage,
            }

    return result_rate_dictionary


def insert_nan_values_for_identical_taxa(
    test_results_branches: TestResultsBranches,
    collapsed_nodes: dict,
    initial_tree: Tree,
    focused_edge: str = None,
):
    """
    Insert NaN values for edges consisting of identical taxa by creating TestResultBranch instances
    with NaN values for specific metrics, and then add these instances to the TestResultsBranches container.

    Args:
        test_results_branches (TestResultsBranches): The container to which the TestResultBranch instances are added.
        collapsed_nodes (dict): A dictionary of parent nodes mapping to their child nodes.
        initial_tree (Tree): The initial phylogenetic tree used for distance calculation.
        focused_edge: A specific edge to focus on. If None, the function processes all edges in `collapsed_nodes`.
    """

    not_analyzed_constant = "duplicate sequence"

    for parent, children in collapsed_nodes.items():
        for child in children:
            if not focused_edge or (parent in focused_edge and child in focused_edge):
                # Prepare the nan_values_dict for this specific edge
                nan_values_dict = {
                    "test_statistic": not_analyzed_constant,
                    "p_value": not_analyzed_constant,
                    "decision_test": not_analyzed_constant,
                    "decision_corrected_test_tips": not_analyzed_constant,
                    "decision_corrected_test_branches": not_analyzed_constant,
                    "decision_test_tip2tip": not_analyzed_constant,
                    # Add the branch length directly to the dictionary
                    "branch_length": initial_tree.get_distance(parent, child),
                }

                # Create a TestResultBranch instance with NaN values for this specific edge
                edge_result = TestResultBranch(**nan_values_dict)
                # Generate a unique branch name for this edge
                branch_name = f"({child},{parent})"
                # Add this TestResultBranch instance to the TestResultsBranches container
                test_results_branches.add_branch(branch_name, edge_result)
