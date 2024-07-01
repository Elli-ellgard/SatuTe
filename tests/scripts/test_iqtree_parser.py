import numpy as np
from satute.iqtree_parser import IqTreeParser


def test_load_substitution_gtr_f():
    parser = IqTreeParser("tests/data/point_iqtree_files/example_model_GTR_F.iqtree")
    # Simulating valid and invalid file contents for parsing tests        
    substitution_model = parser.load_substitution_model()   
    
    # Define the expected matrix
    expected_matrix = substitution_model.precomputed_q_matrix

    assert np.allclose(substitution_model.rate_matrix, expected_matrix, rtol=0.1, atol=0.1)
    assert substitution_model.model == "GTR+F"
    assert substitution_model.number_rates == 1        
    assert [0.267, 0.2783, 0.2962, 0.1585] == list(substitution_model.state_frequencies)

def test_load_substitution_gtr_f_g4():
    parser = IqTreeParser("tests/data/point_iqtree_files/example_model_GTR_F_G4.iqtree")
    # Simulating valid and invalid file contents for parsing tests        
    substitution_model = parser.load_substitution_model()        
    
    # Define the expected matrix
    expected_matrix = substitution_model.precomputed_q_matrix
    
    assert np.allclose(substitution_model.rate_matrix, expected_matrix, rtol=0.01, atol=0.01)
    assert substitution_model.number_rates == 4
    assert substitution_model.model == "GTR+F+G4"
    assert np.allclose([0.2515, 0.23, 0.3121, 0.2065], list(substitution_model.state_frequencies), atol=0.01)

def test_load_substitution_jc():
    parser = IqTreeParser("tests/data/point_iqtree_files/example_model_JC.iqtree")
    # Simulating valid and invalid file contents for parsing tests        
    substitution_model = parser.load_substitution_model()        
    substitution_model.state_frequencies   
    
    # Define the expected matrix
    expected_matrix = substitution_model.precomputed_q_matrix
    
    assert np.allclose(substitution_model.rate_matrix, expected_matrix, rtol=0.01, atol=0.01)
    assert substitution_model.model == "JC"
    assert substitution_model.number_rates == 1    
    assert [0.25,0.25,0.25,0.25] == list(substitution_model.state_frequencies)


def test_load_substitution_k2p_precomputed_rate_matrix():
    parser = IqTreeParser("tests/data/point_iqtree_files/example_model_K2P.iqtree")
    # Simulating valid and invalid file contents for parsing tests        
    substitution_model = parser.load_substitution_model()        
    substitution_model.state_frequencies   
    # Define the expected matrix
    precomputed_matrix = substitution_model.precomputed_q_matrix
        
    assert np.allclose(substitution_model.rate_matrix, precomputed_matrix, rtol=0.01, atol=0.01)
    assert substitution_model.model == "K2P"
    assert substitution_model.number_rates == 1        
    assert [0.25,0.25,0.25,0.25] == list(substitution_model.state_frequencies)

