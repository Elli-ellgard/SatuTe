import sys

sys.path.append("./../../")

import unittest
from satute_repository import IqTreeParser

class TestIqTreeParser(unittest.TestCase):
    def setUp(self):
        # Mock data and instance setup for the tests
        self.parser = IqTreeParser(
            "../data/just_iq_tree_log_files/toy_example_ntaxa_7_run_1-alignment.phy.txt.iqtree"
        )
        self.parser.file_content = []

    def test_build_rate_matrix_basic(self):
        """Test the build_rate_matrix function with typical input."""
        only_rates = [0.5, 0.3, 0.2]
        list_state_freq = [0.25, 0.25, 0.25, 0.25]
        n = 2
        expected_matrix = np.array([[-0.5, 0.5], [0.5, -0.5]])

        result_matrix = self.parser.build_rate_matrix(only_rates, list_state_freq, n)

        np.testing.assert_array_almost_equal(result_matrix, expected_matrix)

    def test_build_rate_matrix_zero_rates(self):
        """Test the function with zero rates to check for no changes or errors."""
        only_rates = [0, 0]
        list_state_freq = [0.5, 0.5]
        n = 2
        expected_matrix = np.array([[0, 0], [0, 0]])
        result_matrix = self.parser.build_rate_matrix(only_rates, list_state_freq, n)
        np.testing.assert_array_almost_equal(result_matrix, expected_matrix)

    def test_build_rate_matrix_invalid_input(self):
        """Test the function with invalid (negative) inputs."""
        only_rates = [-0.5, -0.3]
        list_state_freq = [0.25, 0.75]
        n = 2
        with self.assertRaises(ValueError):
            self.parser.build_rate_matrix(only_rates, list_state_freq, n)

    def test_build_rate_matrix_non_square(self):
        """Test the function to ensure it handles non-square matrices correctly."""
        only_rates = [0.1, 0.2, 0.3]
        list_state_freq = [0.33, 0.33, 0.34]
        n = 3
        expected_matrix = np.zeros(
            (3, 3)
        )  # Should calculate based on the actual rates and frequencies
        result_matrix = self.parser.build_rate_matrix(only_rates, list_state_freq, n)
        self.assertEqual(
            result_matrix.shape, (3, 3)
        )  # Check matrix is square of correct size


if __name__ == "__main__":
    # unittest.main()
    parser = IqTreeParser(
        "../data/just_iq_tree_log_files/toy_example_ntaxa_7_run_1-alignment.phy.txt.iqtree"
    )
    current_substitution_model = parser.parse_substitution_model()

    parser.check_model(current_substitution_model)
    parser.model_type = parser.get_model_type(current_substitution_model)

    parser.load_iqtree_file_content()

    dict_state_frequencies, phi_matrix = parser.parse_state_frequencies()
    state_frequencies_dict = dict_state_frequencies

    rate_matrix = parser.construct_rate_matrix(state_frequencies_dict)
