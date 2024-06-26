import unittest



def sum_positive_integers(n):
    return n*(n+1)/2

def sum_squares_positive_integers(n):
    return n*(n+1)*(2*n+1)/6


class TestCalculateSampleCoherence(unittest.TestCase):

    def run_test(self, number_sites, example_function, result_function):
        from satute.ztest_posterior_distribution import (calculate_sample_coherence)
        with self.subTest(number_sites=number_sites):
            multiplicity, factors_left_subtree, factors_right_subtree, number_sites = example_function(number_sites)
            sample_coherence = calculate_sample_coherence(multiplicity, factors_left_subtree, factors_right_subtree, number_sites)
            result = result_function(number_sites)
            self.assertEqual(sample_coherence, result/number_sites)

    def run_test_with_zeros(self, nzeros, n, example_function, result_function):
        from satute.ztest_posterior_distribution import (calculate_sample_coherence)
        with self.subTest(n=n):
            multiplicity, factors_left_subtree, factors_right_subtree, number_sites = example_function(nzeros, n)
            sample_coherence = calculate_sample_coherence(multiplicity, factors_left_subtree, factors_right_subtree, number_sites)
            result = result_function(n)/(n + nzeros)*multiplicity
            self.assertEqual(sample_coherence, result)

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
