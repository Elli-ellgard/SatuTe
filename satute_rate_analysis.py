import pandas as pd
from satute_trees_and_subtrees import rescale_branch_lengths, collapse_tree
from partial_likelihood import calculate_partial_likelihoods_for_sites
from graph import count_leaves_and_branches_for_subtrees
from satute_sequences import dict_to_alignment
from satute_statistic_posterior_distribution import (
    calculate_test_statistic_posterior_distribution,
)


def determine_branch_type(edge):
    """
    Determine the type of the branch based on the edge.

    Args:
        edge (tuple): Edge details.

    Returns:
        str: 'internal' if the edge is internal, 'external' otherwise.
    """
    return "internal" if "Node" in edge[1] and "Node" in edge[0] else "external"


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

    # Calculate partial likelihoods for all sites
    partial_likelihood_per_site_storage = calculate_partial_likelihoods_for_sites(
        initial_tree, alignment, rate_matrix, focused_edge
    )

    # Count leaves and branches for subtrees
    edge_subtree_count_dict = count_leaves_and_branches_for_subtrees(
        initial_tree, alignment, focused_edge
    )

    # Initialize a dictionary and a list to store results
    result_test_dictionary = {}
    result_list = []

    # Iterate over each edge and process likelihoods
    for edge, likelihoods in partial_likelihood_per_site_storage.items():
        # Convert left and right likelihoods to dataframes
        left_partial_likelihood = pd.DataFrame(likelihoods["left"]["likelihoods"])
        right_partial_likelihood = pd.DataFrame(likelihoods["right"]["likelihoods"])

        # Count the number of leaves in the left and right subtree
        number_leaves_left_subtree = edge_subtree_count_dict[edge]["left"][
            "leave_count"
        ]
        number_leaves_right_subtree = edge_subtree_count_dict[edge]["right"][
            "leave_count"
        ]

        # Determine branch type (internal or external)
        branch_type = determine_branch_type(edge)

        # Calculate test statistics using the posterior distribution
        results = process_test_statistics_posterior(
            multiplicity,
            array_right_eigenvectors,
            state_frequencies,
            left_partial_likelihood,
            right_partial_likelihood,
            number_leaves_left_subtree,
            number_leaves_right_subtree,
            branch_type,
            alpha,
        )

        # Store the results for the given edge
        result_list.append(
            store_test_results(edge, "single_rate", left_partial_likelihood, results)
        )

    # Add all results to the main result dictionary
    result_test_dictionary["single_rate"] = result_list
    return result_test_dictionary


def process_test_statistics_posterior(
    multiplicity,
    array_right_eigenvectors,
    state_frequencies,
    left_partial_likelihood,
    right_partial_likelihood,
    number_leaves_left_subtree,
    number_leaves_right_subtree,
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
        number_leaves_left_subtree,
        number_leaves_right_subtree,
        branch_type,
        alpha,
    )

    result_keys = [
        "delta",
        "p_value",
        "decision_corrected_test_tips",
        "decision_corrected_test_branches",
        "result_test_tip2tip",
        "decision_test_tip2tip",
    ]

    return {key: value for key, value in zip(result_keys, results)}


# TODO: Generate leave count, edge dictionary only once
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
    # Initialize a dictionary to store results for each rate category
    result_rate_dictionary = {}

    # Iterate over each rate category and its corresponding alignment
    for rate, sub_alignment in per_rate_category_alignment.items():
        # Step 1: Map sequence IDs to their sequences from the alignment for easy access
        sequence_dict = {record.id: str(record.seq) for record in sub_alignment}

        # Step 2: Retrieve the relative rate for the current category
        relative_rate = category_rates_factors[rate]["Relative_rate"]

        # Initialize a list to store results for the current rate category
        result_list = []

        # Step 3: Create a deep copy of the initial tree and collapse nodes with identical sequences
        collapsed_rescaled_tree_one = initial_tree.copy("deepcopy")
        sequence_dict, twin_dictionary = collapse_tree(
            collapsed_rescaled_tree_one, sequence_dict
        )
        # Step 4: Rescale branch lengths according to the relative rate
        rescale_branch_lengths(collapsed_rescaled_tree_one, relative_rate)

        # Step 5: Convert the sequence dictionary back to a MultipleSeqAlignment object after collapsing nodes
        sub_alignment = dict_to_alignment(sequence_dict)

        # Step 6: Calculate partial likelihoods for all sites in the rescaled tree
        partial_likelihood_per_site_storage = calculate_partial_likelihoods_for_sites(
            collapsed_rescaled_tree_one, sub_alignment, rate_matrix, focused_edge
        )

        # Step 7: Count leaves and branches for subtrees in the rescaled tree
        edge_subtree_count_dict = count_leaves_and_branches_for_subtrees(
            initial_tree, sub_alignment, focused_edge
        )

        # Step 8: Process each edge and its associated likelihoods
        for edge, likelihoods in partial_likelihood_per_site_storage.items():
            # Convert left and right likelihoods to dataframes
            left_partial_likelihood = pd.DataFrame(likelihoods["left"]["likelihoods"])
            right_partial_likelihood = pd.DataFrame(likelihoods["right"]["likelihoods"])

            # Get the number of leaves in the left and right subtree
            number_leaves_left_subtree = edge_subtree_count_dict[edge]["left"][
                "leave_count"
            ]
            number_leaves_right_subtree = edge_subtree_count_dict[edge]["right"][
                "leave_count"
            ]

            # Determine the branch type (internal or external)
            branch_type = determine_branch_type(edge)

            # Calculate test statistics using the posterior distribution
            results = process_test_statistics_posterior(
                multiplicity,
                array_right_eigenvectors,
                state_frequencies,
                left_partial_likelihood,
                right_partial_likelihood,
                number_leaves_left_subtree,
                number_leaves_right_subtree,
                branch_type,
                alpha,
            )

            # Store the results of the test for the given edge
            result_list.append(
                store_test_results(edge, rate, left_partial_likelihood, results)
            )

        for parent, value in twin_dictionary.items():
            for child in value:
                result_list.append(
                    {
                        "edge": f"({parent},{child})",
                        "delta": "Nan",
                        "p_value": "Nan",
                        "decision_corrected_test_tips": "Nan",
                        "decision_corrected_test_branches": "Nan",
                        "decision_test_tip2tip": "Nan",
                        "result_test": "Nan",
                        "result_test_tip2tip": "Nan",
                        "category_rate": rate,
                        "branch_length": initial_tree.get_distance(parent, child),
                    }
                )

        # Step 9: Augment the results with additional data for the 'twin' nodes
        # augmented_results = augment_results_with_twin_data(
        #    twin_dictionary, result_list, collapsed_rescaled_tree_one, rate=rate
        # )

        # Step 10: Add the results for the current rate category to the main dictionary
        result_rate_dictionary[rate] = {
            "result_list": result_list,
            "rescaled_tree": collapsed_rescaled_tree_one,
        }

    return result_rate_dictionary


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
        "p_value": results.get("p_value"),
        "decision_corrected_test_tips": results.get("decision_corrected_test_tips"),
        "decision_corrected_test_branches": results.get(
            "decision_corrected_test_branches"
        ),
        "decision_test_tip2tip": results.get("decision_test_tip2tip"),
        "category_rate": rate,
        "branch_length": left_partial_likelihood.get("branch_length", [None])[0],
    }


def single_rate_analysis_collapsed_tree(
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
    # Step 1: Map sequence IDs to their sequences from the alignment for easy access
    sequence_dict = {record.id: str(record.seq) for record in alignment}

    # Step 2: Create a deep copy of the initial tree to modify without altering the original
    collapsed_tree_one = initial_tree.copy("deepcopy")

    # Step 3: Collapse nodes in the tree with identical sequences to simplify the tree structure
    sequence_dict, twin_dictionary = collapse_tree(collapsed_tree_one, sequence_dict)

    # Step 4: Convert the sequence dictionary back to a MultipleSeqAlignment object after collapsing nodes
    alignment = dict_to_alignment(sequence_dict)

    # Step 5: Carry out single rate analysis using the collapsed tree and updated alignment
    results = single_rate_analysis(
        collapsed_tree_one,
        alignment,
        rate_matrix,
        state_frequencies,
        array_left_eigenvectors,
        array_right_eigenvectors,
        multiplicity,
        alpha,
        focused_edge,
    )

    # Step 6: Augment the results with additional data for the 'twin' nodes
    results = augment_results_with_twin_data(
        twin_dictionary, results, initial_tree, "single_rate"
    )

    # Step 7: Compile all the results into a dictionary and return it
    return results


def augment_results_with_twin_data(twin_dictionary, results, initial_tree, rate):
    """
    Augments the analysis results with additional data for the 'twin' nodes.

    Args:
    twin_dictionary (dict): A dictionary mapping parent nodes to their 'twin' children.
    results (dict): The results dictionary where augmented data will be appended.
    initial_tree (Tree): The initial tree used in the analysis for calculating distances.

    Returns:
    dict: The augmented results dictionary.
    """
    for parent, value in twin_dictionary.items():
        for child in value:
            results[rate].append(
                {
                    "edge": f"({parent},{child})",
                    "delta": "Nan",
                    "p_value": "Nan",
                    "decision_corrected_test_tips": "Nan",
                    "decision_corrected_test_branches": "Nan",
                    "decision_test_tip2tip": "Nan",
                    "result_test": "Nan",
                    "result_test_tip2tip": "Nan",
                    "category_rate": "single_rate_analysis",
                    "branch_length": initial_tree.get_distance(parent, child),
                }
            )
    return results
