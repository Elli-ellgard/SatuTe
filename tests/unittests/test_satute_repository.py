import numpy as np
import unittest
from satute.repository import IqTreeParser

class TestIqTreeParser(unittest.TestCase):
    def setUp(self):
        # Mock data and instance setup for the tests
        self.parser = IqTreeParser(
            "tests/data/just_iq_tree_log_files/toy_example_ntaxa_7_run_1-alignment.phy.txt.iqtree"
        )

    def test_build_rate_matrix_basic(self):
        """Test the build_rate_matrix function with typical input."""
        self.parser.load_iqtree_file_content()        
        substitution_model = self.parser.load_substitution_model()
        
                

if __name__ == "__main__":
    unittest.main()
