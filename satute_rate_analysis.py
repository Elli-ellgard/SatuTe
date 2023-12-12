from satute_trees_and_subtrees import rescale_branch_lengths
import pandas as pd
from partial_likelihood import (
    calculate_partial_likelihoods_for_sites,
    determine_branch_type,
)
from satute_statistic_posterior_distribution import (
    calculate_test_statistic_posterior_distribution,
)


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
