import numpy as np
from satute.parser.iqtree_parser import IqTreeParser
from satute.models.amino_acid_models import check_rate_matrix

def test_load_substitution_gtr_f():
    parser = IqTreeParser("tests/data/point_iqtree_files/example_dna_model_GTR_F.iqtree")
    # Simulating valid and invalid file contents for parsing tests        
    substitution_model = parser.load_substitution_model()   
    
    # Define the expected matrix
    expected_matrix = substitution_model.precomputed_q_matrix

    assert np.allclose(substitution_model.rate_matrix, expected_matrix, rtol=0.1, atol=0.1)
    assert substitution_model.model == "GTR+F"
    assert substitution_model.number_rates == 1        
    assert [0.267, 0.2783, 0.2962, 0.1585] == list(substitution_model.state_frequencies)

def test_load_substitution_gtr_f_g4():
    parser = IqTreeParser("tests/data/point_iqtree_files/example_dna_model_GTR_F_G4.iqtree")
    # Simulating valid and invalid file contents for parsing tests        
    substitution_model = parser.load_substitution_model()        
    
    # Define the expected matrix
    expected_matrix = substitution_model.precomputed_q_matrix
    
    assert np.allclose(substitution_model.rate_matrix, expected_matrix, rtol=0.01, atol=0.01)
    assert substitution_model.number_rates == 4
    assert substitution_model.model == "GTR+F+G4"
    assert np.allclose([0.2515, 0.23, 0.3121, 0.2065], list(substitution_model.state_frequencies), atol=0.01)

def test_load_substitution_jc():
    parser = IqTreeParser("tests/data/point_iqtree_files/example_dna_model_JC.iqtree")
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
    parser = IqTreeParser("tests/data/point_iqtree_files/example_dna_model_K2P.iqtree")
    # Simulating valid and invalid file contents for parsing tests        
    substitution_model = parser.load_substitution_model()        
    substitution_model.state_frequencies   
    # Define the expected matrix
    precomputed_matrix = substitution_model.precomputed_q_matrix
        
    assert np.allclose(substitution_model.rate_matrix, precomputed_matrix, rtol=0.01, atol=0.01)
    assert substitution_model.model == "K2P"
    assert substitution_model.number_rates == 1        
    assert [0.25,0.25,0.25,0.25] == list(substitution_model.state_frequencies)

def test_load_substitution_k2p_precomputed_rate_matrix():
    parser = IqTreeParser("tests/data/point_iqtree_files/example_dna_model_K2Pb.iqtree")
    # Simulating valid and invalid file contents for parsing tests        
    substitution_model = parser.load_substitution_model()        
    substitution_model.state_frequencies   
    # Define the expected matrix
    precomputed_matrix = substitution_model.precomputed_q_matrix
        
    assert np.allclose(substitution_model.rate_matrix, precomputed_matrix, rtol=0.01, atol=0.01)
    assert substitution_model.model == "K2P"
    assert substitution_model.number_rates == 1        
    assert [0.25,0.25,0.25,0.25] == list(substitution_model.state_frequencies)


def test_load_substitution_ry3_4_precomputed_rate_matrix():
    parser = IqTreeParser("tests/data/point_iqtree_files/example_dna_model_RY3.4_G4.iqtree")
    # Simulating valid and invalid file contents for parsing tests        
    substitution_model = parser.load_substitution_model()        
    substitution_model.state_frequencies   
    # Define the expected matrix
    precomputed_matrix = substitution_model.precomputed_q_matrix
        
    assert np.allclose(substitution_model.rate_matrix, precomputed_matrix, rtol=0.1, atol=0.1)
    assert substitution_model.model == "RY3.4+G4"
    assert substitution_model.number_rates == 4        
    assert substitution_model.category_rates ==  {
                                                    'p1':
                                                    {
                                                        "Category": '1',
                                                        "Relative_rate": float(0.02862),
                                                        "Proportion": float(0.25)
                                                    },
                                                    'p2':{
                                                        "Category": '2',
                                                        "Relative_rate": float(0.2344),
                                                        "Proportion": float(0.25)                                                    
                                                        },
                                                    'p3': {
                                                        "Category": '3',
                                                        "Relative_rate": float(0.8001),
                                                        "Proportion": float(0.25)                                                    
                                                    },
                                                    'p4':{
                                                        "Category": '4',
                                                        "Relative_rate": float(2.937),
                                                        "Proportion": float(0.25)                                                    
                                                    }   
                                                }
    assert [0.2193,0.2807,0.2193,0.2807] == list(substitution_model.state_frequencies)


def test_load_substitution_ry3_4b_precomputed_rate_matrix():
    parser = IqTreeParser("tests/data/point_iqtree_files/example_dna_model_RY4.4b.iqtree")
    # Simulating valid and invalid file contents for parsing tests        
    substitution_model = parser.load_substitution_model()        
    substitution_model.state_frequencies   
    # Define the expected matrix
    precomputed_matrix = substitution_model.precomputed_q_matrix
        
    assert np.allclose(substitution_model.rate_matrix, precomputed_matrix, rtol=0.1, atol=0.1)
    assert substitution_model.model == "RY4.4b"
    assert substitution_model.number_rates == 1
    assert substitution_model.category_rates == None
    assert [0.2408,0.2592,0.2408,0.2592] == list(substitution_model.state_frequencies)


def test_load_substitution_gtr_20_precomputed_rate_matrix():
    parser = IqTreeParser("tests/data/point_iqtree_files/example_aa_model_GTR20.iqtree")
    # Simulating valid and invalid file contents for parsing tests        
    substitution_model = parser.load_substitution_model()        
    assert np.allclose(substitution_model.state_frequencies,[ 0.0790,0.0843, 0.0346, 0.0392, 0.0055, 0.0293, 0.0535, 0.0956, 0.0217, 0.0704, 0.0699, 0.1008, 0.0253, 0.0280, 0.0415, 0.0519, 0.0516, 0.0065, 0.0217, 0.0898, ], rtol=0.01)
    assert substitution_model.model == "GTR20+F"   
    assert substitution_model.number_rates == 1
    assert substitution_model.rate_matrix.shape == (20,20)
    outcome ,message = check_rate_matrix(substitution_model.rate_matrix, substitution_model.state_frequencies)
    assert outcome, f"{message}"
    
