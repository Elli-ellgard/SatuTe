import scipy.stats as st
import numpy as np
import pandas as pd

"""## CALCULATION OF THE SAMPLE COHERENCE """

def calculate_sample_coherence(
      multiplicity,
      factors_left_subtree,
      factors_right_subtree,
      number_sites    
):
    delta = 0
    for i in range(multiplicity):
        delta += np.asarray(factors_left_subtree[i]) @ np.asarray(factors_right_subtree[i])

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

    if multiplicity == 1:
            if branch_type == "internal":
                M_left = np.asarray(factors_left_subtree[0]) @ np.asarray(factors_left_subtree[0]) / number_sites 
                M_left = min(1, M_left)
            else:  # if clade A is a single leaf
                M_left = 1

            M_right = np.asarray(factors_right_subtree[0]) @ np.asarray(factors_right_subtree[0]) / number_sites
            M_right = min(1, M_right)
            variance = M_left * M_right
        
    else:
        for i in range(multiplicity):
            for j in range(multiplicity):
                    M_left = np.asarray(factors_left_subtree[i])@np.asarray(factors_left_subtree[i])/number_sites
                    M_left = min(1, M_left)
                    M_right = np.asarray(factors_right_subtree[j])@np.asarray(factors_right_subtree[j])/number_sites
                    M_right = min(1, M_right)
                    variance += M_right * M_left
    return variance

def calculate_alternative_variance(
    multiplicity,
    factors_left_subtree,
    factors_right_subtree,
    number_sites         
):
    variance = 0
    for i in range(multiplicity):
        for j in range(multiplicity):
            square_factor_left = np.asarray(factors_left_subtree[i]) * np.asarray(factors_left_subtree[i]) 
            square_factor_right = np.asarray(factors_right_subtree[j]) * np.asarray(factors_right_subtree[j])
            variance += square_factor_left @ square_factor_right
    variance = variance / number_sites
    return variance

"""## CALCULATION OF THE TEST STATISTIC"""

def calculate_test_statistic(
    multiplicity,
    array_eigenvectors,
    posterior_probabilities_left_subtree,
    posterior_probabilities_right_subtree,
    dimension,
    branch_type = "external",
    alpha = 0.05,
 ):  
    # quantile of the standard normal distribution
    z_alpha = st.norm.ppf(1-alpha)
    
    number_sites = len(posterior_probabilities_left_subtree["Site"].unique())
    number_nodes_1 = len(posterior_probabilities_left_subtree["Node"].unique())
    number_nodes_2 = len(posterior_probabilities_right_subtree["Node"].unique())

    """ Calculation of the factors for the coefficient C_1 (correspond to the dominant non-zero eigenvalue)"""
    # list of vetors of the scalar products between right eigenvector and posterior probability per site 
    factors_left_subtree = ([]) # list of vectors
    factors_right_subtree = ([])

    for i in range(multiplicity):
        v = array_eigenvectors[i] # eigenvector v_i of the dominant non-zero eigenvalue

        a = [] # vector to store all scalar products v_i * site_posterior_probabilities_left_subtree
        b = [] # vector to store all scalar products v_i * site_posterior_probabilities_right_subtree

        for k in range(
            number_sites * (number_nodes_1 - 1), number_sites * number_nodes_1
        ):
            a.append(
                v @ np.asarray(posterior_probabilities_left_subtree.iloc[k, 3:(3 + dimension)])
            )
        factors_left_subtree.append(a)

        for k in range(
            number_sites * (number_nodes_2 - 1), number_sites * number_nodes_2
        ):
            b.append(
                v @ np.asarray(posterior_probabilities_right_subtree.iloc[k, 3:(3 + dimension)])
            )
        factors_right_subtree.append(b)

    """ calculation of the dominant sample coherence """
    delta = calculate_sample_coherence(
        multiplicity,
        factors_left_subtree,
        factors_right_subtree,
        number_sites
    )
    print(delta)

    """ calculation of the sample variance """
    variance = calculate_sample_variance(
        multiplicity,
        factors_left_subtree,
        factors_right_subtree,
        number_sites,
        branch_type,
    )
    variance = variance / np.sqrt(number_sites)
    if variance < 0:
        print(
            "VARIANCE ESTIMATION IS NEGATIVE - CONSIDER INCREASING THE NUMBER OF STANDARD DEVIATIONS (number_standard_deviations) (CONFIDENCE INTERVAL)"
        )
        c_s = 999999999
        p_value = -1
    else:
        # computing the saturation coherence / critical value
        c_s = z_alpha * np.sqrt(variance)
        print(c_s)
        # computing the  p-value
        p_value = st.norm.sf(abs(delta / np.sqrt(variance)))

    if c_s > delta:
        result_test = "Saturated"
    else:
        result_test = "Informative"

    # computing the saturation coherence between two sequences
    c_sTwoSequence = (
                multiplicity * z_alpha / np.sqrt(number_sites)
            )  
    
    if c_sTwoSequence > delta:
        result_test_tip2tip = "SatuT2T"
    else:
            result_test_tip2tip = "InfoT2T"

    return delta, c_s, c_sTwoSequence, p_value, result_test, result_test_tip2tip 

def old_calculation(
    multiplicity,
    array_eigenvectors,
    posterior_probabilities_left_subtree,
    posterior_probabilities_right_subtree,
    dimension,
    branch_type = "external",
    alpha = 0.05,
):
        
    z_alpha = st.norm.ppf(1-alpha)
    number_sites = len(posterior_probabilities_left_subtree["Site"].unique())
    number_nodes_1 = len(posterior_probabilities_left_subtree["Node"].unique())
    number_nodes_2 = len(posterior_probabilities_right_subtree["Node"].unique())

    if multiplicity == 1:  # if D=1
        v1 = array_eigenvectors[0]

        a = (
            []
        )  # vector to store all products v1*rootsitesposteriorprobabilitiescladeA
        b = (
            []
        )  # vector to store all products v1*rootsitesposteriorprobabilitiescladeB

        for k in range(
            number_sites * (number_nodes_1 - 1), number_sites * number_nodes_1
        ):
            a.append(
                v1 @ np.asarray(posterior_probabilities_left_subtree.iloc[k, 3:7])
            )

        for k in range(
            number_sites * (number_nodes_2 - 1), number_sites * number_nodes_2
        ):
            b.append(
                v1 @ np.asarray(posterior_probabilities_right_subtree.iloc[k, 3:7])
            )

        delta = np.asarray(a) @ np.asarray(b) / number_sites
        # computing the dominant sample coherence

        if branch_type == "internal":
            M_a = np.asarray(a) @ np.asarray(a) / number_sites 
            M_a = min(1, M_a)
        else:  # if clade A is a single leaf
            M_a = 1

        M_b = np.asarray(b) @ np.asarray(b) / number_sites
        M_b = min(1, M_b)
        variance = M_a * M_b / np.sqrt(number_sites)
        c_s = z_alpha * np.sqrt(variance)  # computing the saturation coherence

        c_sTwoSequence = z_alpha / np.sqrt(
            number_sites
        )  # computing the saturation coherence between two sequences

        p_value = st.norm.sf(abs(delta / np.sqrt(variance)))
    else:
        c_sTwoSequence = (
            multiplicity * z_alpha / np.sqrt(number_sites)
        )  # computing the saturation coherence between two sequences
        delta = 0

        for j in range(multiplicity):
            a = []
            b = []
            v1 = array_eigenvectors[j]

            for k in range(
                number_sites * (number_nodes_1 - 1), number_sites * number_nodes_1
            ):
                a.append(
                    v1
                    @ np.asarray(posterior_probabilities_left_subtree.iloc[k, 3:7])
                )

            for k in range(
                number_sites * (number_nodes_2 - 1), number_sites * number_nodes_2
            ):
                b.append(
                    v1
                    @ np.asarray(posterior_probabilities_right_subtree.iloc[k, 3:7])
                )

            delta += np.asarray(a) @ np.asarray(b)
        delta = delta / number_sites

        variance = 0

        for j in range(multiplicity):
            for k in range(multiplicity):
                a = []
                b = []

                v_j = array_eigenvectors[j]
                v_k = array_eigenvectors[k]

                for l in range(
                    number_sites * (number_nodes_1 - 1),
                    number_sites * number_nodes_1,
                ):
                    a.append(
                        v_j
                        @ np.asarray(
                            posterior_probabilities_left_subtree.iloc[l, 3:7]
                        )
                    )

                for l in range(
                    number_sites * (number_nodes_2 - 1),
                    number_sites * number_nodes_2,
                ):
                    b.append(
                        v_k
                        @ np.asarray(
                            posterior_probabilities_right_subtree.iloc[l, 3:7]
                        )
                    )

                # variance = np.asarray(a)@np.asarray(b)
                variance += max(
                    np.asarray(a) @ np.asarray(b),
                    np.asarray(a) @ np.asarray(b),
                    np.asarray(a) @ np.asarray(b),
                    np.asarray(a) @ np.asarray(b),
                )

        variance = variance / (number_sites * number_sites)

        if variance < 0:
            print(
                "VARIANCE ESTIMATION IS NEGATIVE - CONSIDER INCREASING THE NUMBER OF STANDARD DEVIATIONS (number_standard_deviations) (CONFIDENCE INTERVAL)"
            )
            c_s = 999999999
            p_value = -1
        else:
            c_s = z_alpha * np.sqrt(variance)
            p_value = st.norm.sf(abs(delta / np.sqrt(variance)))

    if c_s > delta:
        result_test = "Saturated"
    else:
        result_test = "Informative"

    if c_sTwoSequence > delta:
        result_test_tip2tip = "SatuT2T"
    else:
        result_test_tip2tip = "InfoT2T"
    return delta, c_s, c_sTwoSequence, p_value, result_test, result_test_tip2tip