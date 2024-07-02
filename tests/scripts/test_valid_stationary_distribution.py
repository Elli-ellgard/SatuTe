import unittest

class TestValidStationaryDistribution(unittest.TestCase):
    # def valid_stationary_distribution(self, frequencies):
    #     total = sum(frequencies.values())
    #     if total != 1.0:
    #         frequencies = {k: v / total for k, v in frequencies.items()}
    #     return frequencies

    def test_valid_stationary_distribution_case1(self):
        from satute.parser.iqtree_nucleotide_parser import (valid_stationary_distribution)
        frequencies = {'A': 0.2, 'B': 0.3, 'C': 0.5}
        frequencies = valid_stationary_distribution(frequencies)
        self.assertEqual(1.0, sum(frequencies.values()))

    def test_valid_stationary_distribution_case2(self):
        from satute.parser.iqtree_nucleotide_parser import (valid_stationary_distribution)
        frequencies = {'A': 0.2, 'B': 0.3, 'C': 0.24, 'D': 0.25}
        frequencies = valid_stationary_distribution(frequencies)
        self.assertEqual(1.0, sum(frequencies.values()))

if __name__ == "__main__":
    unittest.main()

