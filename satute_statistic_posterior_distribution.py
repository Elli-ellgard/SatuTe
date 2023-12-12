import scipy.stats as st
import numpy as np
from scipy.sparse.linalg import expm

"""## CALCULATION OF THE SAMPLE COHERENCE """


# get transition matrix using matrix exponential
def get_transition_matrix(rate_matrix, branch_length):
    return expm(rate_matrix * branch_length)


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


"""## CALCULATION OF THE TEST STATISTIC FOR BRANCH SATURATION"""


def calculate_test_statistic_posterior_distribution(
    multiplicity,
    array_eigenvectors,
    state_frequencies,
    partial_likelihood_left_subtree,
    partial_likelihood_right_subtree,
    dimension,
    branch_type="external",
    alpha=0.05,
):
    # quantiles of the standard normal distribution
    z_alpha = st.norm.ppf(1 - alpha)
    number_sites = len(partial_likelihood_left_subtree["Site"].unique())

    """Calculation of the posterior distributions """
    posterior_probabilities_left_subtree = []
    posterior_probabilities_right_subtree = []

    freq = np.array(list(state_frequencies))
    diag = np.diag(list(state_frequencies))

    site_likelihood_left_subtree = []
    site_likelihood_right_subtree = []

    for k in range(number_sites):
        sr = np.dot(
            np.asarray(partial_likelihood_right_subtree.iloc[k, 3 : (3 + dimension)]),
            freq,
        )
        sl = np.sum(
            np.asarray(partial_likelihood_left_subtree.iloc[k, 3 : (3 + dimension)])
            @ diag
        )
        site_likelihood_right_subtree.append(sr)
        site_likelihood_left_subtree.append(sl)
        posterior_probabilities_left_subtree.append(
            np.array(
                diag
                @ np.asarray(
                    partial_likelihood_left_subtree.iloc[k, 3 : (3 + dimension)]
                )
            )
            / sl
        )

        posterior_probabilities_right_subtree.append(
            np.array(
                diag
                @ np.asarray(
                    partial_likelihood_right_subtree.iloc[k, 3 : (3 + dimension)]
                )
            )
            / sr
        )

    """ Calculation of the factors for the coefficient C_1 (correspond to the dominant non-zero eigenvalue)"""
    # list of vectors of the scalar products between right eigenvector and posterior probability per site
    factors_left_subtree = []  # list of vectors
    factors_right_subtree = []

    for i in range(multiplicity):
        v = array_eigenvectors[i]  # eigenvector v_i of the dominant non-zero eigenvalue

        a = (
            []
        )  # vector to store all scalar products v_i * site_posterior_probabilities_left_subtree
        b = (
            []
        )  # vector to store all scalar products v_i * site_posterior_probabilities_right_subtree

        for k in range(number_sites):
            a.append(v @ np.asarray(posterior_probabilities_left_subtree[k]))
            b.append(v @ np.asarray(posterior_probabilities_right_subtree[k]))
        factors_right_subtree.append(b)
        factors_left_subtree.append(a)

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
        # computing the p-value
        p_value = st.norm.sf(delta / np.sqrt(variance))

    if c_s > delta:
        result_test = "Saturated"
    else:
        result_test = "Informative"

    # computing the saturation coherence between two sequences
    c_s_two_sequence = np.sqrt(multiplicity) * z_alpha / np.sqrt(number_sites)

    if c_s_two_sequence > delta:
        result_test_tip2tip = "SatuT2T"
    else:
        result_test_tip2tip = "InfoT2T"

    return delta, c_s, c_s_two_sequence, p_value, result_test, result_test_tip2tip
