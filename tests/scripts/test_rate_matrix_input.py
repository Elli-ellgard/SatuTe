from satute.amino_acid_models import create_rate_matrix_with_input, check_rate_matrix

import numpy as np
import unittest
from satute.amino_acid_models import AMINO_ACID_RATE_MATRIX, AA_STATE_FREQUENCIES

class TestRateMatrixCreation(unittest.TestCase):
    def test_create_rate_matrix_with_input(self):
        matrix_size = 4
        input_string = """1
                          1 1 
                          1 1 1"""
                            
        eq = np.array([0.25, 0.25, 0.25,0.25])
        
        # Expected result calculated manually or with a known correct algorithm
        expected_matrix = np.array([
            [-0.75, 0.25, 0.25, 0.25],
            [0.25, -0.75, 0.25, 0.25],
            [0.25, 0.25, -0.75, 0.25],
            [0.25, 0.25, 0.25, -0.75]
        ])
        
        result_matrix = create_rate_matrix_with_input(matrix_size, input_string, eq)

        # Use numpy's allclose function to compare floating point matrices
        self.assertTrue(np.allclose(result_matrix, expected_matrix),f"Expected matrix did not match result:\n{result_matrix}\nexpected:\n{expected_matrix}")
        
    def test_rate_matrix_properties(self):
        matrix_size = 4
        input_string = """1
                      1 1
                      1 1 1"""
        eq = np.array()
    
        result_matrix = create_rate_matrix_with_input(matrix_size, input_string, eq)
    
        # Check the properties of the rate matrix
        result, message = check_rate_matrix(result_matrix, eq)
        self.assertTrue(result, message)
        
    def test_rate_matrix_properties(self):
        matrix_size = 20
        input_string = AMINO_ACID_RATE_MATRIX['LG']
        eq = AA_STATE_FREQUENCIES["LG"]
    
        result_matrix = create_rate_matrix_with_input(matrix_size, input_string, eq)
    
        # Check the properties of the rate matrix
        result, message = check_rate_matrix(result_matrix, eq)
        self.assertTrue(result, message)
        
    def test_model_rate_matrices(self):
        for model_name, frequencies in AA_STATE_FREQUENCIES.items():
            with self.subTest(model=model_name):
                matrix_size = 20  # Assuming all matrices are of size 20x20
                rate_matrix_data = AMINO_ACID_RATE_MATRIX[model_name]
                
                # This function needs to be defined to create a rate matrix from the provided data
                result_matrix = create_rate_matrix_with_input(matrix_size, rate_matrix_data, frequencies)
                
                # Check the properties of the rate matrix
                result, message = check_rate_matrix(result_matrix, frequencies)
                self.assertTrue(result, f"{model_name} model failed: {message}")
        