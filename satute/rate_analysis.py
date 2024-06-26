# -*- coding: utf-8 -*-
import pandas as pd
from ete3 import Tree
from Bio.Align import MultipleSeqAlignment
from typing import List

from satute.partial_likelihood import calculate_partial_likelihoods_for_sites
from satute.graph import calculate_subtree_edge_metrics
from satute.sequences import dict_to_alignment
from satute.rate_matrix import RateMatrix
from satute.ztest_posterior_distribution import (
    calculate_test_statistic_posterior_distribution,
)
from satute.trees import (
    rescale_branch_lengths,
)

from satute.result import (
    TestResultsBranches,
    TestStatisticComponentsContainer,
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
        rescaled_tree = initial_tree.copy("deepcopy")

        # Step 4: Rescale branch lengths according to the relative rate
        rescale_branch_lengths(rescaled_tree, relative_rate)

        # Step 5: Convert the sequence dictionary back to a MultipleSeqAlignment object after collapsing nodes
        sub_alignment = dict_to_alignment(sequence_dict)

        if sub_alignment.get_alignment_length() != 0:
            # Step 6: Calculate partial likelihoods for all sites in the rescaled tree
            partial_likelihood_per_site_storage = (
                calculate_partial_likelihoods_for_sites(
                    rescaled_tree,
                    sub_alignment,
                    rate_matrix,
                    focused_edge,
                )
            )

            # Step 7: Count leaves and branches for subtrees in the rescaled tree
            edge_subtree_metrics = calculate_subtree_edge_metrics(rescaled_tree, focused_edge)

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


            # Step 10: Add the results for the current rate category to the main dictionary
            result_rate_dictionary[rate] = {
                "result_list": results_list,
                "rescaled_tree": rescaled_tree,
                "components": components_container,
                "partial_likelihoods": partial_likelihood_per_site_storage,
            }

    return result_rate_dictionary


