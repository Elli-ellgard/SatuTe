import scipy.stats as st
import numpy as np
from satute_repository import (
    parse_state_frequencies,
    parse_rate_matrices_from_file,
)
from scipy.sparse.linalg import expm


# get transition matrix using matrix exponential
def get_transition_matrix(rate_matrix, branch_length):
    return expm(rate_matrix * branch_length)


"""## CALCULATION OF THE SAMPLE COHERENCE """


def calculate_sample_coherence(
    multiplicity, factors_left_subtree, factors_right_subtree, number_sites
):
    delta = 0
    for i in range(multiplicity):
        delta += np.asarray(factors_left_subtree[i]) @ np.asarray(
            factors_right_subtree[i]
        )

    delta = delta / number_sites
    return delta


"""## ESTIMATION OF THE SAMPLE VARIANCE"""


def calculate_sample_variance(
    multiplicity,
    factors_left_subtree,
    factors_right_subtree,
    number_sites,
    branch_type,
):
    variance = 0

    for i in range(multiplicity):
        for j in range(multiplicity):
            if branch_type == "internal":
                M_left = (
                    np.asarray(factors_left_subtree[i])
                    @ np.asarray(factors_left_subtree[i])
                    / number_sites
                )
            else:
                M_left = multiplicity
            M_right = (
                np.asarray(factors_right_subtree[j])
                @ np.asarray(factors_right_subtree[j])
                / number_sites
            )
            variance += M_right * M_left
    return variance


def calculate_alternative_variance(
    multiplicity, factors_left_subtree, factors_right_subtree, number_sites
):
    variance = 0
    for i in range(multiplicity):
        for j in range(multiplicity):
            square_factor_left = np.asarray(factors_left_subtree[i]) * np.asarray(
                factors_left_subtree[i]
            )
            square_factor_right = np.asarray(factors_right_subtree[j]) * np.asarray(
                factors_right_subtree[j]
            )
            variance += square_factor_left @ square_factor_right
    variance = variance / number_sites
    return variance


def calculate_sample_variance_optimized(
    multiplicity,
    factors_left_subtree,
    factors_right_subtree,
    number_sites,
    branch_type,
):
    # Convert the input lists to numpy arrays for efficient matrix operations
    factors_left_subtree = np.asarray(factors_left_subtree)
    factors_right_subtree = np.asarray(factors_right_subtree)

    # Compute M_left:
    # If the branch type is 'internal', compute M_left for each factor in the left subtree.
    # Otherwise, set M_left to a constant value of multiplicity for all factors.
    if branch_type == "internal":
        # '@' is the matrix multiplication operator in Python.
        # The expression computes the dot product of each factor with itself.
        # `.diagonal()` extracts the diagonal elements from the resulting matrix.
        # Each diagonal element is then divided by `number_sites`.
        M_left = (
            factors_left_subtree @ factors_left_subtree.T
        ).diagonal() / number_sites
    else:
        # If branch type is not 'internal', create an array filled with the value of `multiplicity`.
        M_left = np.full(multiplicity, multiplicity)

    # Compute M_right in a similar manner to M_left, but for the right subtree factors.
    M_right = (
        factors_right_subtree @ factors_right_subtree.T
    ).diagonal() / number_sites

    # Broadcasting:
    # The following line computes the outer product of M_left and M_right.
    # This is equivalent to multiplying each element of M_left with each element of M_right.
    # The resulting matrix contains all the products of pairs of elements from M_left and M_right.
    variance_matrix = M_left[:, None] * M_right

    # Sum all elements of the variance matrix to get the final variance value.
    return variance_matrix.sum()


# Optimized function with comments
def calculate_sample_coherence_optimized(
    multiplicity, factors_left_subtree, factors_right_subtree, number_sites
):
    # Convert the input lists of factors to numpy arrays.
    # This allows for efficient vectorized operations.
    factors_left_subtree = np.asarray(factors_left_subtree)
    factors_right_subtree = np.asarray(factors_right_subtree)

    # Multiply corresponding factors (elements) from the left and right subtrees element-wise.
    # Then, sum the results to get the total delta value.
    # This replaces the loop in the original function.
    delta = np.sum(factors_left_subtree * factors_right_subtree)

    # Normalize the delta value by dividing by the number of sites.
    delta = delta / number_sites

    return delta


"""## CALCULATION OF THE TEST STATISTIC FOR BRANCH SATURATION"""
# now the left eigenvectors h_i are necessary


def calculate_test_statistic(
    multiplicity,
    array_left_eigenvectors,
    partial_likelihood_left_subtree,
    partial_likelihood_right_subtree,
    dimension,
    branch_type="external",
    alpha=0.05,
):
    # quantiles of the standard normal distribution
    z_alpha = st.norm.ppf(1 - alpha)
    number_sites = len(partial_likelihood_left_subtree["Site"].unique())

    """ Calculation of the factors for the coefficient C_1 (correspond to the dominant non-zero eigenvalue)"""
    # list of vectors of the scalar products between right eigenvector and posterior probability per site
    factors_left_subtree = []  # list of vectors
    factors_right_subtree = []

    for i in range(multiplicity):
        h = array_left_eigenvectors[
            i
        ]  # left eigenvector h_i of the dominant non-zero eigenvalue

        a = (
            []
        )  # vector to store all scalar products v_i * site_partial_likelihood_left_subtree
        b = (
            []
        )  # vector to store all scalar products v_i * site_partial_likelihood_right_subtree

        for k in range(number_sites):
            a.append(
                h
                @ np.asarray(
                    partial_likelihood_left_subtree.iloc[k, 3 : (3 + dimension)]
                )
            )

            b.append(
                h
                @ np.asarray(
                    partial_likelihood_right_subtree.iloc[k, 3 : (3 + dimension)]
                )
            )

        factors_left_subtree.append(a)
        factors_right_subtree.append(b)

    """ calculation of the dominant sample coherence """
    delta = calculate_sample_coherence(
        multiplicity, factors_left_subtree, factors_right_subtree, number_sites
    )

    """ calculation of the sample variance """
    variance = calculate_sample_variance(
        multiplicity,
        factors_left_subtree,
        factors_right_subtree,
        number_sites,
        branch_type,
    )

    variance = variance / number_sites
    if variance < 0:
        print(
            "VARIANCE ESTIMATION IS NEGATIVE - CONSIDER INCREASING THE NUMBER OF STANDARD DEVIATIONS (number_standard_deviations) (CONFIDENCE INTERVAL)"
        )
        c_s = 999999999
        p_value = -1
    else:
        # computing the saturation coherence / critical value
        c_s = z_alpha * np.sqrt(variance)
        # computing the  p-value
        p_value = st.norm.sf(abs(delta / np.sqrt(variance)))

    if c_s > delta:
        result_test = "Saturated"
    else:
        result_test = "Informative"

    # computing the saturation coherence between two sequences
    c_s_two_sequence = multiplicity * z_alpha / np.sqrt(number_sites)

    if c_s_two_sequence > delta:
        result_test_tip2tip = "SatuT2T"
    else:
        result_test_tip2tip = "InfoT2T"

    return delta, c_s, c_s_two_sequence, p_value, result_test, result_test_tip2tip


"""## CALCULATION OF THE TEST STATISTIC FOR LIKELIHOOD RATIO TEST"""


def calculate_likelihood_ratio_test(
    input_directory,
    branch_length,
    partial_likelihood_left_subtree,
    partial_likelihood_right_subtree,
    dimension,
    alpha,
):
    state_frequencies = parse_state_frequencies(
        f"{input_directory}.iqtree", dimension=dimension
    )
    diag = np.diag(list(state_frequencies.values()))

    (rate_matrix, phi_matrix) = parse_rate_matrices_from_file(
        f"{input_directory}.iqtree"
    )
    transition_matrix = get_transition_matrix(rate_matrix, branch_length)

    number_sites = len(partial_likelihood_left_subtree["Site"].unique())
    number_nodes_1 = len(partial_likelihood_left_subtree["Node"].unique())
    number_nodes_2 = len(partial_likelihood_right_subtree["Node"].unique())

    """ Calculation of likelihood ratios per site"""
    likelihood_ratios = []
    # list of vectors
    factors_left = (
        []
    )  # partial likelihood of left subtree times diagonal matix of the stationary distribution
    factors_right = []  # transition matrix  times partial likelihood of right subtree
    site_likelihood_left = []
    site_likelihood_right = []
    vector = []

    for k in range(number_sites * (number_nodes_1 - 1), number_sites * number_nodes_1):
        vector = (
            np.asarray(partial_likelihood_left_subtree.iloc[k, 3 : (3 + dimension)])
            @ diag
        )
        factors_left.append(vector)
        site_likelihood_left.append(sum(vector))

    for k in range(number_sites * (number_nodes_2 - 1), number_sites * number_nodes_2):
        factors_right.append(
            transition_matrix
            @ np.asarray(
                partial_likelihood_right_subtree.iloc[k, 3 : (3 + dimension)]
            ).transpose()
        )
        vector = (
            np.asarray(partial_likelihood_right_subtree.iloc[k, 3 : (3 + dimension)])
            @ diag
        )
        site_likelihood_right.append(sum(vector))

    for k in range(number_sites):
        likelihood_ratios.append(
            np.dot(np.asarray(factors_left[k]), np.asarray(factors_right[k]))
        )
    test_statistic = 2 * np.sum(
        np.log(likelihood_ratios)
        - np.prod(site_likelihood_left)
        - np.prod(site_likelihood_right)
    )

    # critical value of the distribution for alpha = 0.05
    critical_value = 2.71
    if test_statistic > critical_value:
        result_lr_test = "Informative"
    else:
        result_lr_test = "Saturated"
    return test_statistic, result_lr_test
