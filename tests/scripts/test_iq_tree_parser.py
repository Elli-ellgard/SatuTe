import unittest
import numpy as np
from satute.repository import IqTreeParser
from satute.exceptions import ModelNotFoundError

class TestIqTreeParser(unittest.TestCase):
    def setUp(self):
        self.parser = IqTreeParser()
        # Simulating valid and invalid file contents for parsing tests
        self.valid_content = [
            "Best-fit model according to BIC: LG",
            "Rate matrix Q:",
            "A    -1    0.3333    0.3333    0.3333",
            "C    0.3333    -1    0.3333    0.3333",
            "G    0.3333    0.3333    -1    0.3333",
            "T    0.3333    0.3333    0.3333    -1"
        ]
        self.invalid_content = ["Random text", "Not the right content"]

    def test_get_aa_rate_matrix_valid_model(self):
        # Assuming 'LG' is a valid model known to be in the dictionary
        model = 'LG'
        expected_shape = (20, 20)  # Typically, the rate matrix for amino acids is 20x20
        rate_matrix = self.parser.get_aa_rate_matrix(model)
        self.assertEqual(rate_matrix.shape, expected_shape, "Rate matrix should be 20x20 for the LG model")

    def test_get_aa_rate_matrix_invalid_model(self):
        # Testing for an invalid model name
        model = 'InvalidModel'
        with self.assertRaises(ModelNotFoundError):
            self.parser.get_aa_rate_matrix(model)

    def test_parse_aa_rate_matrix_valid(self):
        self.parser.file_path("../data/just_iq_tree_log_files/ali-len1500nt_cbl3.5.fasta.iqtree") 
        self.parser.load_iqtree_file_content()
        matrix = self.parser.parse_aa_rate_matrix()
        expected_shape = (20, 20)  # Based on the mocked file content with a 4x4 matrix
        self.assertEqual(matrix.shape, expected_shape, "Parsed matrix should be 4x4")

    def test_parse_aa_rate_matrix_not_found(self):
        self.parser.file_content = self.invalid_content
        with self.assertRaises(ValueError):
            self.parser.parse_aa_rate_matrix()

    def test_parse_substitution_model_found(self):
        self.parser.file_content = ["Best-fit model according to BIC: LG"]
        expected_model = "LG"
        model = self.parser.parse_substitution_model() 
        self.assertEqual(model, expected_model, "Should return the correct model name")

    def test_parse_substitution_model_not_found(self):
        self.parser.file_content = self.invalid_content
        with self.assertRaises(ValueError):
            self.parser.parse_substitution_model()

    def test_load_substitution_model_integration(self):
        # This is a more integrated test to check model loading, assuming valid data
        self.parser.file_content = self.valid_content
        sub_model = self.parser.load_substitution_model()
        self.assertIsNotNone(sub_model, "Should successfully create a SubstitutionModel object")
        self.assertTrue(hasattr(sub_model, 'state_frequencies'), "SubstitutionModel should have state frequencies attribute")

if __name__ == '__main__':
    unittest.main()