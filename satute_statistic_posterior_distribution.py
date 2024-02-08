import scipy.stats as st
import numpy as np
from scipy.sparse.linalg import expm




class TestStatisticComponents:
    def __init__(self, coefficients: list[float], variances: list[float]):
        self.coefficients = coefficients
        self.variances = variances

   
        

class TestResultBranch:
    """Class representing all results for the test for saturation."""

    def __init__(self, **result):
        self.result = result

    def add_result(self, result_name, score):
        self.results[result_name] = score

    def get_result(self):
        return self.result
        


"""## CALCULATION OF POSTERIOR DISTRIBUTION """

def calculate_posterior_probabilities_subtree(dimension, state_frequencies, partial_likelihood_subtree, number_sites):
    posterior_probabilities_subtree = []

    #freq = np.array(list(state_frequencies.values()))
    diag = np.diag(list(state_frequencies.values()))
    site_likelihood_subtree= []
    for k in range(number_sites):
        # sr = np.dot(
        #     np.asarray(partial_likelihood_subtree.iloc[k, 3 : (3 + dimension)]),
        #     freq,
        # )
        sl = np.sum(
            np.asarray(partial_likelihood_subtree.iloc[k, 3 : (3 + dimension)])
            @ diag
        )
        site_likelihood_subtree.append(sl)
        posterior_probabilities_subtree.append(
            np.array(
                diag
                @ np.asarray(
                    partial_likelihood_subtree.iloc[k, 3 : (3 + dimension)]
                )
            )
            / sl
        )
    return posterior_probabilities_subtree



"""## CALCULATION OF FACTOR FOR C_1"""

def scalar_product_eigenvector_posterior_probability(multiplicity, array_eigenvectors, posterior_probabilities, number_sites):
    factors_subtree = []  # list of vectors
    
    for i in range(multiplicity):
        v = array_eigenvectors[i]  # eigenvector v_i of the dominant non-zero eigenvalue
        a = (
            []
        )  # vector to store all scalar products v_i * site_posterior_probabilities

        for k in range(number_sites):
            a.append(v @ np.asarray(posterior_probabilities[k]))
        factors_subtree.append(a)
    return factors_subtree



"""## CALCULATION OF THE SAMPLE COHERENCE """

def get_number_of_branch_insertions(number_tips):
    if number_tips == 1 or number_tips == 2:
        number_branches = 1
    elif number_tips == 3:
        number_branches = 3
    else:
        number_branches =  2 * number_tips - 3
    return(number_branches)
    
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

"""## ESTIMATION OF THE SAMPLE VARIANCE """

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

""" New calculation of test-statistic by excluding zeros"""


def calculate_sample_coherence_for_each_site(
    multiplicity, factors_left_subtree, factors_right_subtree, number_sites
):
    delta = np.zeros(number_sites)
    for i in range(multiplicity):
        delta += np.asarray(factors_left_subtree[i]) * np.asarray(factors_right_subtree[i])

    return delta

def calculate_sample_coherence_without_zeros():

    result = 0
    return result


def calculate_components_of_test_statistic_for_each_site(
    multiplicity,
    factors_left_subtree,
    factors_right_subtree,
    number_sites,
    branch_type,
):
    variance = np.zeros(number_sites)
    delta = np.zeros(number_sites)
    for i in range(multiplicity):
        delta += np.asarray(factors_left_subtree[i]) * np.asarray(factors_right_subtree[i])
        for j in range(multiplicity):
            if branch_type == "internal":
                m_left = (
                    np.asarray(factors_left_subtree[i])
                    * np.asarray(factors_left_subtree[j])
                    
                )
            else:
                if i == j:
                    m_left = np.ones(number_sites)
                else:
                    m_left =np.zeros(number_sites)

            m_right = (
                np.asarray(factors_right_subtree[i])
                * np.asarray(factors_right_subtree[j])
            )
            variance += m_right * m_left

    return delta, variance


def calculate_test_statistic_exlude_zeros(
    coefficients,
    variances,
    number_sites,
    branch_type,
):
    sample_mean_sum = 0
    sample_variance_sum = 0
    number_informative_sites =0
    for i in range(number_sites):
        if(coefficients[i] != 0):
            sample_mean_sum += coefficients[i]
            sample_variance_sum += variances[i]
            number_informative_sites += 1
    if number_informative_sites > 0:
        sample_mean = sample_mean_sum/number_informative_sites
    else:
        sample_mean = np.nan
    if sample_variance_sum > 0:
        if branch_type == "internal":
            test_statistic = sample_mean_sum/np.sqrt(sample_variance_sum/number_informative_sites)
        else: 
            test_statistic = sample_mean_sum/np.sqrt(sample_variance_sum)
    else:
        test_statistic = np.nan
            
    #print("delta new:", sample_mean)
    # if branch_type == "internal":
    #     sample_variance= sample_variance/number_informative_sites/number_informative_sites/number_informative_sites
    # else: 
    #     sample_variance= sample_variance/number_informative_sites/number_informative_sites
    # #print("variance new:", sample_variance)
    # test_statistic= sample_mean / np.sqrt(sample_variance)
    return test_statistic, sample_mean, number_informative_sites


"""## DECISION OF STATTISTICAL TEST """

def decision_z_test(test_statistic, alpha):
    # quantile of the standard normal distribution
    z_alpha = st.norm.ppf(1 - alpha)
    # calculate the critical value
    c_s = z_alpha 
    # decision using critical value
    decision_test = ""   
    if c_s > test_statistic:
        decision_test = "Saturated"
    else:
        decision_test = "Informative"
    return decision_test

def decision_tip2tip(delta, number_sites, multiplicity, alpha):
    # quantile of the standard normal distribution
    z_alpha = st.norm.ppf(1 - alpha)
    # calculate the critical value
    c_s_two_sequence = np.sqrt(multiplicity) * z_alpha / np.sqrt(number_sites)
    # decision using critcial value
    decision_test_tip2tip = ""
    if c_s_two_sequence > delta:
        decision_test_tip2tip = "SatuT2T"
    else:
        decision_test_tip2tip = "InfoT2T"
    return decision_test_tip2tip



""" ## BONFERRONI CORRECTION """

def bonferroni_test_correction_tips(p_value, number_tips_left_subtree, number_tips_right_subtree, alpha):
    # calculate corrected significance level
    corrected_alpha_tips = alpha/(number_tips_left_subtree*number_tips_right_subtree)
    # decision using p-value
    decision_corrected_test_tips = ""
    if p_value > corrected_alpha_tips:
        decision_corrected_test_tips = "Saturated"
    else:
        decision_corrected_test_tips = "Informative" 
    return decision_corrected_test_tips

def bonferroni_test_correction_branches(p_value, number_tips_left_subtree, number_tips_right_subtree, alpha):
    # determine the number of possible branch insertion for the two subtrees
    number_branch_insertion_left_subtree = get_number_of_branch_insertions(number_tips_left_subtree)
    number_branch_insertion_right_subtree = get_number_of_branch_insertions(number_tips_right_subtree)
    # calculate  corrected significance level
    corrected_alpha_branches= alpha/(number_branch_insertion_left_subtree*number_branch_insertion_right_subtree)
    # decision using p-value
    decision_corrected_test_branches = ""
    if p_value > corrected_alpha_branches:
        decision_corrected_test_branches = "Saturated"
    else:
        decision_corrected_test_branches = "Informative" 
    return decision_corrected_test_branches



"""## SIDAK CORRECTION """

def sidak_test_correction_tips(test_statistic, number_tips_left_subtree, number_tips_right_subtree, alpha):
    # calculate corrected critical value
    corrected_alpha_tips = 1 - (1 - alpha)**(1/(number_tips_left_subtree*number_tips_right_subtree))
    corrected_c_s_tips = st.norm.ppf(1 - corrected_alpha_tips)
    # decision using critical value
    decision_corrected_test_tips = ""
    if corrected_c_s_tips > test_statistic:
        decision_corrected_test_tips = "Saturated"
    else:
        decision_corrected_test_tips = "Informative" 
    return decision_corrected_test_tips

def sidak_test_correction_tips(test_statistic, number_tips_left_subtree, number_tips_right_subtree, alpha):
    # determine the number of possible branch insertion for the two subtrees   
    number_branch_insertion_left_subtree = get_number_of_branch_insertions(number_tips_left_subtree)
    number_branch_insertion_right_subtree = get_number_of_branch_insertions(number_tips_right_subtree)
    # calculate corrected critical value
    corrected_alpha_branches = 1 - (1 - alpha)**(1/(number_branch_insertion_left_subtree*number_branch_insertion_right_subtree))
    corrected_c_s_branches = st.norm.ppf(1 - corrected_alpha_branches)
    # decision using critical value  
    decision_corrected_test_branches = ""
    if corrected_c_s_branches > test_statistic:
        decision_corrected_test_branches = "Saturated"
    else:
        decision_corrected_test_branches = "Informative" 
    return decision_corrected_test_branches



"""## CALCULATION OF THE TEST STATISTIC FOR BRANCH SATURATION"""

def calculate_test_statistic_posterior_distribution(
    multiplicity,
    array_eigenvectors,
    state_frequencies,
    partial_likelihood_left_subtree,
    partial_likelihood_right_subtree,
    dimension,
    number_tips_left_subtree,
    number_tips_right_subtree,
    branch_type="external",
    alpha=0.05,
): #-> tuple [TestStatisticComponents, TestResultBranch]:
    # quantiles of the standard normal distribution
    z_alpha = st.norm.ppf(1 - alpha)
    number_sites = len(partial_likelihood_left_subtree["Site"].unique())

    """ Calculation of the posterior distributions """
    posterior_probabilities_left_subtree = calculate_posterior_probabilities_subtree(dimension, state_frequencies, partial_likelihood_left_subtree, number_sites)
    posterior_probabilities_right_subtree = calculate_posterior_probabilities_subtree(dimension, state_frequencies, partial_likelihood_right_subtree, number_sites)

    """ Calculation of the factors for the coefficient C_1 (correspond to the dominant non-zero eigenvalue)"""
    # list of vectors of the scalar products between right eigenvector and posterior probability per site
    factors_left_subtree = scalar_product_eigenvector_posterior_probability(multiplicity, array_eigenvectors, posterior_probabilities_left_subtree, number_sites)
    factors_right_subtree = scalar_product_eigenvector_posterior_probability(multiplicity, array_eigenvectors, posterior_probabilities_right_subtree, number_sites)

    # """ Calculation of the sample mean of coefficient C_1 """
    # delta = calculate_sample_coherence(
    #     multiplicity, factors_left_subtree, factors_right_subtree, number_sites
    # )

    # """ calculation of the sample variance """
    # variance = calculate_sample_variance(
    #     multiplicity,
    #     factors_left_subtree,
    #     factors_right_subtree,
    #     number_sites,
    #     branch_type,
    # )
    # variance = variance / number_sites

    # """Calculation of the test-statistic and decision of the statistical tests """
    # if variance > 0:
    #     test_statistic = delta / np.sqrt(variance)
    
    #     """Calculation of the p-value"""
    #     p_value = st.norm.sf(test_statistic)

    #     """ Results of the statistcal tests"""
    #     # decision of the statistical test
    #     decision_test = decision_z_test(test_statistic, alpha)

    #     # decision of the test using Bonferroni correction 
    #     # using number of tips of the considered subtrees
    #     decision_corrected_test_tips = bonferroni_test_correction_tips(p_value, number_tips_left_subtree, number_tips_right_subtree, alpha)
        
    #     # using number of branch combinations
    #     decision_corrected_test_branches = bonferroni_test_correction_branches(p_value, number_tips_left_subtree, number_tips_right_subtree, alpha)
       
    # else: 
    #     test_statistic = np.nan
    #     p_value = np.nan
    #     decision_test = np.nan
    #     decision_corrected_test_tips = np.nan
    #     decision_corrected_test_branches = np.nan

    # """ Calculation of the saturation coherence between two sequences """
    # decision_test_tip2tip = decision_tip2tip(delta, number_sites, multiplicity, alpha)

    (
        coefficients,
        variances
    ) = calculate_components_of_test_statistic_for_each_site(
        multiplicity,
        factors_left_subtree,
        factors_right_subtree,
        number_sites,
        branch_type,
    )

    components = TestStatisticComponents(
        coefficients,
        variances
    )

    (
        test_statistic,
        coefficient_value, 
        number_informative_sites 
    )= calculate_test_statistic_exlude_zeros( 
        coefficients,
        variances,
        number_sites,
        branch_type
    )
    if test_statistic != np.nan:

        """Calculation of the p-value"""
        p_value = st.norm.sf(test_statistic)

        """ Results of the statistcal tests"""
        # decision of the statistical test
        decision_test = decision_z_test(test_statistic, alpha)

        # decision of the test using Bonferroni correction 
        # using number of tips of the considered subtrees
        decision_corrected_test_tips = bonferroni_test_correction_tips(p_value, number_tips_left_subtree, number_tips_right_subtree, alpha)
        
        # using number of branch combinations
        decision_corrected_test_branches = bonferroni_test_correction_branches(p_value, number_tips_left_subtree, number_tips_right_subtree, alpha)
       
    else: 
        p_value = np.nan
        decision_test = np.nan
        decision_corrected_test_tips = np.nan
        decision_corrected_test_branches = np.nan


    """ Calculation of the saturation coherence between two sequences """
    decision_test_tip2tip = decision_tip2tip(coefficient_value, number_informative_sites, multiplicity, alpha)
    
    result = TestResultBranch( 
        test_statistic = test_statistic,
        number_informative_sites = number_informative_sites,
        p_value = p_value,
        decision_test = decision_test,
        decision_corrected_test_tips= decision_corrected_test_tips,
        decision_corrected_test_branches = decision_corrected_test_branches,
        decision_test_tip2tip = decision_test_tip2tip
    )

    # return components, result
    return test_statistic, p_value, decision_test, decision_corrected_test_tips, decision_corrected_test_branches, decision_test_tip2tip


