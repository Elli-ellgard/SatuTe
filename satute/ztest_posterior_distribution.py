# -*- coding: utf-8 -*-
import scipy.stats as st
import numpy as np
import pandas as pd
from typing import List, Tuple
from satute.result import TestResultBranch, TestStatisticComponents


"""## CALCULATION OF POSTERIOR DISTRIBUTION """


def calculate_posterior_probabilities_subtree_df(
    dimension: int, state_frequencies: List[float], partial_likelihood_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Calculate the posterior probabilities for each site in a subtree based on the partial likelihoods of evolutionary states
    at each site and the state frequencies.

    Parameters:
    - dimension (int): The number of states considered.
    - state_frequencies (List[float]): A list of state frequencies. The length of this list should match the dimension.
    - partial_likelihood_df (pd.DataFrame): A DataFrame containing partial likelihoods. The relevant likelihood values should start from the fourth column.

    Returns:
    - pd.DataFrame: A DataFrame with the posterior probabilities for each site.

    """
    diag = np.diag(list(state_frequencies))

    # Selecting the relevant columns for likelihoods
    likelihood_cols = partial_likelihood_df.iloc[:, 3 : (3 + dimension)]

    # Calculate the site likelihood for each site (row)
    site_likelihoods = likelihood_cols @ diag
    site_likelihoods_sum = site_likelihoods.sum(axis=1)

    # Calculate the posterior probabilities for each site
    posterior_probabilities = site_likelihoods.divide(site_likelihoods_sum, axis=0)

    return posterior_probabilities


"""## CALCULATION OF FACTOR FOR C_1"""


def scalar_product_eigenvector_posterior_probability(
    multiplicity: int,
    array_eigenvectors: List[np.ndarray],
    posterior_probabilities_df: pd.DataFrame,
    number_sites: int,
) -> List[List[float]]:
    factors_subtree = []  # list of vectors

    # Convert DataFrame posterior probabilities to NumPy array
    posterior_probabilities = posterior_probabilities_df.values

    for i in range(multiplicity):
        v = array_eigenvectors[i]  # eigenvector v_i of the dominant non-zero eigenvalue
        a = []  # vector to store all scalar products v_i * site_posterior_probabilities

        for k in range(number_sites):
            a.append(v @ np.asarray(posterior_probabilities[k]))
        factors_subtree.append(a)
    return factors_subtree


"""## CALCULATION OF THE COEFFICIENTS C_1"""


def calculate_coefficients_of_test_statistic_for_each_site(
    multiplicity: int,
    factors_left_subtree: List[List[float]],
    factors_right_subtree: List[List[float]],
    number_sites: int,
) -> List[float]:
    delta = np.zeros(number_sites)

    for i in range(multiplicity):

        delta += np.asarray(factors_left_subtree[i]) * np.asarray(
            factors_right_subtree[i]
        )

    return delta


# def get_transition_matrix(rate_matrix: np.array, branch_length: float):
#     return expm(rate_matrix * branch_length)

# def calculate_sample_coherence(
#     multiplicity: int,
#     factors_left_subtree: List[float],
#     factors_right_subtree: List[float],
#     number_sites: int,
# ):
#     delta = 0
#     for i in range(multiplicity):
#         delta += np.asarray(factors_left_subtree[i]) @ np.asarray(
#             factors_right_subtree[i]
#         )

#     delta = delta / number_sites
#     return delta


"""## ESTIMATION OF THE SAMPLE VARIANCE """


def calculate_sample_variance(
    multiplicity: int,
    factors_left_subtree: List[float],
    factors_right_subtree: List[float],
    number_sites: int,
    branch_type: str,
) -> float:
    """
    Calculate the sample variance for a branch in a phylogenetic tree by considering contributions from both
    left and right subtrees. Variance computation differs based on the branch type: internal or external.

    Parameters:
    - multiplicity (int): Multiplicity of the considered dominant non-zero eigenvalue.
    - factors_left_subtree (List[float]): List of factors associated with the left subtree.
    - factors_right_subtree (List[float]): List of factors associated with the right subtree.
    - number_sites (int): The number of sites.
    - branch_type (str): Type of the branch. Expected values: 'internal' or 'external'.

    Returns:
    - float: The computed sample variance based on the factors and the branch type.
    """
    variance = 0.0

    for i in range(multiplicity):
        for j in range(multiplicity):
            if branch_type == "internal":
                m_left = (
                    np.asarray(factors_left_subtree[i])
                    @ np.asarray(factors_left_subtree[j])
                    / number_sites
                )
            else:
                m_left = i == j

            m_right = (
                np.asarray(factors_right_subtree[i])
                @ np.asarray(factors_right_subtree[j])
                / number_sites
            )
            variance += m_right * m_left
    return variance


""" CALCULATION OF THE TEST STATISTIC """


def calculate_test_statistic(
    coefficients: List[float], population_variance: float, number_sites: int
) -> Tuple[float, float, float]:
    """
    Calculate the test statistic for branch saturation.

    Parameters:
    - coefficients (List[float]): A list of coefficients C_1 (correspond to the dominant non-zero eigenvalue) per site.
    - population_variance (float): The variance of the population from which the coefficients are sampled.
    - number_sites (int): The number of sites (sample size), which corresponds to the number of coefficients.

    Returns:
    - Tuple[float, float, float]: A tuple containing the test statistic, sample mean and standard error of the mean.
      Returns (nan, nan, nan) if the calculation of the standard error is not possible (e.g., division by zero).
    """
    # Calculate the sample mean of coefficients
    sample_mean_sum = sum(coefficients)
    sample_mean = sample_mean_sum / number_sites

    # Calculate the variance of the sample mean
    sample_variance = population_variance / number_sites
    if sample_variance > 0:
        standard_error_of_mean = np.sqrt(sample_variance)
        test_statistic = sample_mean / standard_error_of_mean
    else:
        standard_error_of_mean = np.nan
        test_statistic = np.nan

    return test_statistic, sample_mean, standard_error_of_mean


# def calculate_test_statistic_exclude_zeros(
#     coefficients,
#     population_variance,
#     number_sites,
# ):
#     sample_mean_sum = 0
#     number_informative_sites = 0
#     for i in range(number_sites):
#         if abs(coefficients[i]) > 10**(-5):
#             sample_mean_sum += coefficients[i]
#             number_informative_sites += 1
#     if number_informative_sites > 0:
#         sample_mean = sample_mean_sum / number_informative_sites
#         sample_variance = population_variance/ number_informative_sites
#     else:
#         sample_mean = np.nan
#     if sample_variance> 0:
#         test_statistic = sample_mean / np.sqrt(sample_variance)
#     else:
#         test_statistic = np.nan
#     return test_statistic, sample_mean, number_informative_sites, sample_variance

""" ## DECISION OF STATISTICAL TEST """


def decision_z_test(test_statistic: float, alpha: float) -> Tuple[float, str]:
    """
    Determine the decision of the one-sided z-test for branch saturation given a test statistic and significance level (alpha).

    Parameters:
    - test_statistic (float): The test statistic value calculated from sample data.
    - alpha (float): The significance level used to determine the critical value from the z-distribution.

    Returns:
    - Tuple[float, str]: A tuple containing the critical value and the decision string.
      The decision is "Saturated" if the test statistic is below the critical value,
      otherwise "Informative".
    """
    # Calculate the critical value from the z-distribution for the given alpha
    z_alpha = st.norm.ppf(1 - alpha)

    # Decision based on the test statistic compared to the critical value
    if z_alpha > test_statistic:
        decision_test = "Saturated"
    else:
        decision_test = "Informative"

    return z_alpha, decision_test


def decision_tip2tip(
    coefficient_value: float, number_sites: int, multiplicity: int, alpha: float
) -> str:
    """
    Determines the classification of the saturation coherence between two sequences.

    Parameters:
    - coefficient_value (float): The sample mean of the coefficients C_1 (correspond to the dominant non-zero eigenvalue)
    - number_sites (int): The number of sites.
    - multiplicity (int): Multiplicity of the considered dominant non-zero eigenvalue.
    - alpha (float): The significance level.

    Returns:
    - str: Returns 'SatuT2T' if the test statistic is below the critical value, indicating saturation. Returns 'InfoT2T' if above, indicating informative results.
    """
    # Calculate the critical value from the z-distribution for the given alpha
    z_alpha = st.norm.ppf(1 - alpha)
    # Calculate the critical value adjusted for multiplicity and number of sites
    c_s_two_sequence = np.sqrt(multiplicity) * z_alpha / np.sqrt(number_sites)

    # Determine the decision based on the critical value
    if c_s_two_sequence > coefficient_value:
        decision_test_tip2tip = "SatuT2T"
    else:
        decision_test_tip2tip = "InfoT2T"
    return decision_test_tip2tip


""" ## BONFERRONI CORRECTION """


def bonferroni_test_correction_tips(
    p_value: float,
    number_tips_left_subtree: int,
    number_tips_right_subtree: int,
    alpha: float,
) -> Tuple[float, str]:
    """
    Apply Bonferroni correction to adjust the significance level for multiple comparisons and determine decision of the test based on the p-value and corrected alpha.

    Parameters:
    - p_value (float): The p-value obtained from the original statistical test.
    - number_tips_left_subtree (int): Number of tips in the left subtree.
    - number_tips_right_subtree (int): Number of tips in the right subtree.
    - alpha (float): Original significance level for the test.

    Returns:
    - Tuple containing:
        - corrected_critical_value (float): The critical value from the z-distribution based on the corrected alpha.
        - decision_corrected_test_tips (str): 'Saturated' if the test is considered statistically insignificant post-correction, 'Informative' otherwise.
    """
    # Calculate corrected significance level using Bonferroni adjustment
    corrected_alpha_tips = alpha / (
        number_tips_left_subtree * number_tips_right_subtree
    )

    # Calculate the critical value from the corrected significance level
    corrected_critical_value = st.norm.ppf(1 - corrected_alpha_tips)

    # Make a decision based on the corrected alpha and the observed p-value
    if p_value > corrected_alpha_tips:
        decision_corrected_test_tips = "Saturated"
    else:
        decision_corrected_test_tips = "Informative"
    return corrected_critical_value, decision_corrected_test_tips


# def get_number_of_branch_insertions(number_tips: int):
#     if number_tips == 1 or number_tips == 2:
#         number_branches = 1
#     elif number_tips == 3:
#         number_branches = 3
#     else:
#         number_branches = 2 * number_tips - 3
#     return number_branches

# def bonferroni_test_correction_branches(
#     p_value, number_tips_left_subtree, number_tips_right_subtree, alpha
# ):
#     # determine the number of possible branch insertion for the two subtrees
#     number_branch_insertion_left_subtree = get_number_of_branch_insertions(
#         number_tips_left_subtree
#     )
#     number_branch_insertion_right_subtree = get_number_of_branch_insertions(
#         number_tips_right_subtree
#     )
#     # calculate  corrected significance level
#     corrected_alpha_branches = alpha / (
#         number_branch_insertion_left_subtree * number_branch_insertion_right_subtree
#     )
#     # decision using p-value
#     decision_corrected_test_branches = ""
#     if p_value > corrected_alpha_branches:
#         decision_corrected_test_branches = "Saturated"
#     else:
#         decision_corrected_test_branches = "Informative"
#     return decision_corrected_test_branches

# """ ## SIDAK CORRECTION """

# def sidak_test_correction_tips(
#     test_statistic, number_tips_left_subtree, number_tips_right_subtree, alpha
# ):
#     # calculate corrected critical value
#     corrected_alpha_tips = 1 - (1 - alpha) ** (
#         1 / (number_tips_left_subtree * number_tips_right_subtree)
#     )
#     corrected_c_s_tips = st.norm.ppf(1 - corrected_alpha_tips)
#     # decision using critical value
#     decision_corrected_test_tips = ""
#     if corrected_c_s_tips > test_statistic:
#         decision_corrected_test_tips = "Saturated"
#     else:
#         decision_corrected_test_tips = "Informative"
#     return decision_corrected_test_tips

# def sidak_test_correction_branches(
#     test_statistic, number_tips_left_subtree, number_tips_right_subtree, alpha
# ):
#     # determine the number of possible branch insertion for the two subtrees
#     number_branch_insertion_left_subtree = get_number_of_branch_insertions(
#         number_tips_left_subtree
#     )
#     number_branch_insertion_right_subtree = get_number_of_branch_insertions(
#         number_tips_right_subtree
#     )
#     # calculate corrected critical value
#     corrected_alpha_branches = 1 - (1 - alpha) ** (
#         1
#         / (number_branch_insertion_left_subtree * number_branch_insertion_right_subtree)
#     )
#     corrected_c_s_branches = st.norm.ppf(1 - corrected_alpha_branches)
#     # decision using critical value
#     decision_corrected_test_branches = ""
#     if corrected_c_s_branches > test_statistic:
#         decision_corrected_test_branches = "Saturated"
#     else:
#         decision_corrected_test_branches = "Informative"
#     return decision_corrected_test_branches


"""## CALCULATION OF THE TEST STATISTIC FOR BRANCH SATURATION"""


def calculate_test_statistic_posterior_distribution(
    multiplicity: int,
    array_eigenvectors: List[np.ndarray],
    state_frequencies: List[float],
    partial_likelihood_left_subtree: pd.DataFrame,
    partial_likelihood_right_subtree: pd.DataFrame,
    dimension: int,
    number_tips_left_subtree: int,
    number_tips_right_subtree: int,
    branch_type: str = "external",
    alpha: float = 0.05,
) -> Tuple[TestStatisticComponents, TestResultBranch]:
    """Calculate various statistics and decisions for phylogenetic analysis of branch saturation.

    This function calculates coefficients, sample variances, test statistics, p-values, and decisions for hypothesis testing,
    incorporating Bonferroni corrections and tip-to-tip tests.

    Args:
        multiplicity (int): Multiplicity of the considered dominant non-zero eigenvalue.
        array_eigenvectors (List[np.ndarray]): List of eigenvectors.
        state_frequencies (List[float]): List of state frequencies.
        partial_likelihood_left_subtree (pd.DataFrame): DataFrame of partial likelihoods for the left subtree.
        partial_likelihood_right_subtree (pd.DataFrame): DataFrame of partial likelihoods for the right subtree.
        dimension (int): Number of states.
        number_tips_left_subtree (int): Number of tips in the left subtree.
        number_tips_right_subtree (int): Number of tips in the right subtree.
        branch_type (str, optional): Type of the branch (default 'external').
        alpha (float, optional): Significance level (default 0.05).

    Returns:
        Tuple[TestStatisticComponents, TestResultBranch]: Tuple containing components and results of the test.
    """
    number_sites = len(partial_likelihood_left_subtree["Site"].unique())

    """ Calculation of the posterior distributions """
    posterior_probabilities_left_subtree = calculate_posterior_probabilities_subtree_df(
        dimension, state_frequencies, partial_likelihood_left_subtree
    )
    posterior_probabilities_right_subtree = (
        calculate_posterior_probabilities_subtree_df(
            dimension, state_frequencies, partial_likelihood_right_subtree
        )
    )

    """ Calculation of the factors for the coefficient C_1 (correspond to the dominant non-zero eigenvalue)
    and the sample variance """
    # list of vectors of the scalar products between right eigenvector and posterior probability per site
    factors_left_subtree = scalar_product_eigenvector_posterior_probability(
        multiplicity,
        array_eigenvectors,
        posterior_probabilities_left_subtree,
        number_sites,
    )
    factors_right_subtree = scalar_product_eigenvector_posterior_probability(
        multiplicity,
        array_eigenvectors,
        posterior_probabilities_right_subtree,
        number_sites,
    )

    coefficients = calculate_coefficients_of_test_statistic_for_each_site(
        multiplicity,
        factors_left_subtree,
        factors_right_subtree,
        number_sites,
    )

    sample_variance = calculate_sample_variance(
        multiplicity,
        factors_left_subtree,
        factors_right_subtree,
        number_sites,
        branch_type,
    )

    components = TestStatisticComponents(coefficients, [sample_variance] * number_sites)

    """"Calculation of the test statistc and standard error of the mean"""
    (test_statistic, coefficient_value, standard_error_of_mean) = (
        calculate_test_statistic(coefficients, sample_variance, number_sites)
    )

    if test_statistic != np.nan:
        """Calculation of the p-value"""
        p_value = st.norm.sf(test_statistic)

        """ Decision of the statistcal tests"""
        # decision of the statistical test
        z_alpha, decision_test = decision_z_test(test_statistic, alpha)

        # decision of the test using Bonferroni correction
        # using number of tips of the considered subtrees
        z_alpha_corrected, decision_corrected_test_tips = (
            bonferroni_test_correction_tips(
                p_value, number_tips_left_subtree, number_tips_right_subtree, alpha
            )
        )

        # # using number of branch combinations
        # decision_corrected_test_branches = bonferroni_test_correction_branches(
        #     p_value, number_tips_left_subtree, number_tips_right_subtree, alpha
        # )
    else:
        p_value = np.nan
        z_alpha = np.nan
        decision_test = np.nan
        z_alpha_corrected = np.nan
        decision_corrected_test_tips = np.nan
        # decision_corrected_test_branches = np.nan

    """ Calculation of the saturation coherence between two sequences """
    decision_test_tip2tip = decision_tip2tip(
        coefficient_value, number_sites, multiplicity, alpha
    )

    result = TestResultBranch(
        coefficient_value=coefficient_value,
        standard_error_of_mean=standard_error_of_mean,
        z_score=test_statistic,
        p_value=p_value,
        z_alpha=z_alpha,
        decision_test=decision_test,
        z_alpha_corrected=z_alpha_corrected,
        decision_corrected_test_tips=decision_corrected_test_tips,
        decision_test_tip2tip=decision_test_tip2tip,
    )

    return components, result
