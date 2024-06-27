import unittest
import pandas as pd
import numpy as np
from satute.ztest_posterior_distribution import (
    scalar_product_eigenvector_posterior_probability,
    calculate_coefficients_of_test_statistic_for_each_site,
    calculate_sample_variance,
    calculate_test_statistic,
)


def sum_positive_integers(n):
    return n*(n+1)/2

def sum_squares_positive_integers(n):
    return n*(n+1)*(2*n+1)/6



class TestScalarProduct(unittest.TestCase):

    def run_test(self, dimension, number_sites, example_function, result_function):
        with self.subTest(number_sites=number_sites):
            multiplicity, array_eigenvectors, posterior_probs_df, number_sites = (
                example_function(dimension, number_sites)
            )
            
            list_factors = scalar_product_eigenvector_posterior_probability(
                    multiplicity,
                    array_eigenvectors,
                    posterior_probs_df,
                    number_sites,
                )
        
            result = result_function(number_sites)
            self.assertEqual(list_factors, [[result]*number_sites])

    def run_test2(self, n, example_function, result_function):
        with self.subTest(n=n):
            multiplicity, array_eigenvectors, posterior_probs_df, number_sites = (
                example_function(n)
            )
            
            list_factors = scalar_product_eigenvector_posterior_probability(
                    multiplicity,
                    array_eigenvectors,
                    posterior_probs_df,
                    number_sites,
                )
            
            result =  [[i for i in range(1, n+1)]]* multiplicity
            self.assertEqual(list_factors, result)

    def example1(self, n, nsites):
        # sum of the squares of the first n positive integers
        multiplicity = 1
        array_eigenvectors = [np.array([i for i in range(1, n + 1)])]
        data = {
            f'State{i}': [i]*nsites for i in range(1, n + 1)
        }
        posterior_probs = pd.DataFrame(data)
        return multiplicity, array_eigenvectors, posterior_probs, n

    def example2(self, n):
        # sum of the first n positive integers
        multiplicity = 1
        array_eigenvectors = [np.array([i for i in range(1, n + 1)])]

        diag_matrix = np.diag(np.array([1]*n))
        posterior_probs = pd.DataFrame(diag_matrix, columns=[f'State{i}' for i in range(1,n+1)])
        return multiplicity, array_eigenvectors, posterior_probs, n
    
    def example3(self, n):
        # sum of the first n positive integers
        multiplicity = 2
        array_eigenvectors = [np.array([i for i in range(1, n + 1)])]*multiplicity

        diag_matrix = np.diag(np.array([1]*n))
        posterior_probs = pd.DataFrame(diag_matrix, columns=[f'State{i}' for i in range(1,n+1)])
        return multiplicity, array_eigenvectors, posterior_probs, n
    



    def test_calculate_sample_coherence(self):
        number_sites_list = [4, 10, 30]
        for number_sites in number_sites_list:
            self.run_test(4, number_sites, self.example1, sum_squares_positive_integers)
            self.run_test2(number_sites, self.example2, sum_positive_integers)
            self.run_test2(number_sites, self.example3, sum_positive_integers)


class TestCalculateComponents(unittest.TestCase):

    def run_test(self, number_sites, example_function, result_function):
        with self.subTest(number_sites=number_sites):
            multiplicity, factors_left_subtree, factors_right_subtree, number_sites = (
                example_function(number_sites)
            )

            sample_coherences = calculate_coefficients_of_test_statistic_for_each_site(
                    multiplicity,
                    factors_left_subtree,
                    factors_right_subtree,
                    number_sites,
                )
        
            result = result_function(number_sites)
            self.assertEqual(sample_coherences.sum(), result)

    def run_test_with_zeros(self, nzeros, n, example_function, result_function):
        with self.subTest(n=n):
            multiplicity, factors_left_subtree, factors_right_subtree, number_sites = (
                example_function(nzeros, n)
            )
            sample_coherences =  calculate_coefficients_of_test_statistic_for_each_site(
                    multiplicity,
                    factors_left_subtree,
                    factors_right_subtree,
                    number_sites,
                )
            
            result = result_function(n) * multiplicity
            self.assertEqual(sample_coherences.sum(), result)

    def example1(self, n):
        # sum of the squares of the first n positive integers
        multiplicity = 1
        factors_left_subtree = [[i for i in range(1, n + 1)]]
        factors_right_subtree = factors_left_subtree
        return multiplicity, factors_left_subtree, factors_right_subtree, n

    def example2(self, n):
        # sum of the first n positive integers
        multiplicity = 1
        factors_left_subtree = [[i for i in range(1, n + 1)]]
        factors_right_subtree = [[1 for i in range(1, n + 1)]]
        return multiplicity, factors_left_subtree, factors_right_subtree, n

    def example3(self, nzeros, n):
        # sum of the squares of the first n positive integers and nzeros zeros in front
        multiplicity = 1
        factors_left_subtree = [[0] * (nzeros) + list(range(1, n + 1))]
        factors_right_subtree = [[1] * (nzeros) + list(range(1, n + 1))]
        return multiplicity, factors_left_subtree, factors_right_subtree, nzeros + n

    def example4(self, nzeros, n):
        # sum of the squares of the first n positive integers and nzeros zeros in front
        multiplicity = 2
        factors_left_subtree = [[0] * (nzeros) + list(range(1, n + 1))] * multiplicity
        factors_right_subtree = [[1] * (nzeros) + list(range(1, n + 1))] * multiplicity
        return multiplicity, factors_left_subtree, factors_right_subtree, nzeros + n

    def test_calculate_sample_coherence(self):
        number_sites_list = [4, 10, 30]
        for number_sites in number_sites_list:
            self.run_test(number_sites, self.example1, sum_squares_positive_integers)
            self.run_test(number_sites, self.example2, sum_positive_integers)

        for number_sites in number_sites_list:
            nzeros = 6
            self.run_test_with_zeros(
                nzeros, number_sites, self.example3, sum_squares_positive_integers
            )
            self.run_test_with_zeros(
                nzeros, number_sites, self.example4, sum_squares_positive_integers
            )


class TestCalculateTestStatistic(unittest.TestCase):

    def run_test(self, number_sites, example_function):
        with self.subTest(number_sites=number_sites):
            coefficients, population_variance =  example_function(number_sites)
            test_statistic, sample_mean, se = calculate_test_statistic(coefficients,population_variance,number_sites)
            

            self.assertEqual(sample_mean, np.mean(coefficients))
            if population_variance > 0:
                self.assertEqual(test_statistic, np.mean(coefficients)/np.sqrt(population_variance/number_sites))
                self.assertEqual(se, np.sqrt(population_variance/number_sites) )
            else: 
                self.assertTrue(np.isnan(test_statistic))
                self.assertTrue(np.isnan(se))

    def example1(self, n):
        coefficients = np.array([i for i in range(1, n + 1)])
        variance = 1.0
        return coefficients, variance
    
    def example2(self, n):
        coefficients = np.array([np.random.rand() for i in range(1, n + 1)])
        variance = 2.2
        return coefficients, variance
    
    def example3(self, n):
        coefficients = np.array([np.random.rand() for i in range(1, n + 1)])
        variance = -1
        return coefficients, variance
    


    def test_calculate_sample_coherence(self):
        number_sites_list = [4, 10, 30]  
        for number_sites in number_sites_list:
            self.run_test(number_sites, self.example1)
            self.run_test(number_sites, self.example2)
            self.run_test(number_sites, self.example3)





if __name__ == '__main__':
    unittest.main()
