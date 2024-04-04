import sys
sys.path.append("./../..")


from satute_statistic_posterior_distribution_components import(
    calculate_components_of_test_statistic_for_each_site,
    calculate_test_statistic_exclude_zeros
)

import unittest


#    (coefficients, variances) = calculate_components_of_test_statistic_for_each_site(
#         multiplicity,
#         factors_left_subtree,
#         factors_right_subtree,
#         number_sites,
#         branch_type,
#     )

#     (test_statistic, coefficient_value, number_informative_sites) = (
#         calculate_test_statistic_exclude_zeros(
#             coefficients, variances, number_sites, branch_type
#         )
#     )

def sum_positive_integers(n):
    return n*(n+1)/2

def sum_squares_positive_integers(n):
    return n*(n+1)*(2*n+1)/6


class TestCalculateComponents(unittest.TestCase):

    def run_test(self, number_sites, example_function, result_function):
        with self.subTest(number_sites=number_sites):
            multiplicity, factors_left_subtree, factors_right_subtree, number_sites = example_function(number_sites)
            sample_coherence , variances = calculate_components_of_test_statistic_for_each_site(
                multiplicity,
                factors_left_subtree,
                factors_right_subtree,
                number_sites,
                "external")
            result = result_function(number_sites)
            self.assertEqual(sample_coherence.sum(), result)

    def run_test_with_zeros(self, nzeros, n, example_function, result_function):
        with self.subTest(n=n):
            multiplicity, factors_left_subtree, factors_right_subtree, number_sites = example_function(nzeros, n)
            sample_coherence, variances = calculate_components_of_test_statistic_for_each_site(
                multiplicity,
                factors_left_subtree,
                factors_right_subtree,
                number_sites,
                "external")
            result = result_function(n)*multiplicity
            self.assertEqual(sample_coherence.sum(), result)

    def example1(self, n):
        # sum of the squares of the first n positive integers
        multiplicity = 1
        factors_left_subtree = [[i for i in range(1, n+1)]]
        factors_right_subtree = factors_left_subtree
        return multiplicity, factors_left_subtree, factors_right_subtree, n
    
    def example2(self, n):
        # sum of the first n positive integers
        multiplicity = 1
        factors_left_subtree = [[i for i in range(1, n+1)]]
        factors_right_subtree = [[1 for i in range(1, n+1)]]
        return multiplicity, factors_left_subtree, factors_right_subtree, n
    
    def example3(self, nzeros, n):
        # sum of the squares of the first n positive integers and nzeros zeros in front
        multiplicity = 1
        factors_left_subtree = [[0] * (nzeros) + list(range(1, n + 1))]
        factors_right_subtree = [[1]*(nzeros) + list(range(1, n + 1))]
        return multiplicity, factors_left_subtree, factors_right_subtree, nzeros + n
    
    def example4(self, nzeros, n):
        # sum of the squares of the first n positive integers and nzeros zeros in front
        multiplicity = 2
        factors_left_subtree = [[0] * (nzeros) + list(range(1, n + 1))] * multiplicity
        factors_right_subtree = [[1]*(nzeros) + list(range(1, n + 1))] * multiplicity
        return multiplicity, factors_left_subtree, factors_right_subtree, nzeros + n
    


    def test_calculate_sample_coherence(self):
        number_sites_list = [4, 10, 30]  
        for number_sites in number_sites_list:
            self.run_test(number_sites, self.example1, sum_squares_positive_integers)
            self.run_test(number_sites, self.example2, sum_positive_integers)

        for number_sites in number_sites_list:
            nzeros = 6
            self.run_test_with_zeros(nzeros, number_sites, self.example3, sum_squares_positive_integers)
            self.run_test_with_zeros(nzeros, number_sites, self.example4, sum_squares_positive_integers)


if __name__ == '__main__':
    unittest.main()
