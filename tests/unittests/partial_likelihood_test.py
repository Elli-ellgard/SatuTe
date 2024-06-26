import numpy as np
from ete3 import Tree
import unittest
import time
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from satute.spectral_decomposition import spectral_decomposition
from satute.graph import Node
from satute.rate_matrix import RateMatrix
from satute.amino_acid_models import POISSON_RATE_MATRIX, AA_STATE_FREQUENCIES
from satute.rate_analysis import single_rate_analysis
from satute.sequences import dict_to_alignment
from satute.trees import (
    rename_internal_nodes_pre_order,
)
from satute.partial_likelihood import (
    partial_likelihood,
)

RATE_MATRIX = np.array([[-3, 1, 1, 1], [1, -3, 1, 1], [1, 1, -3, 1], [1, 1, 1, -3]])

""" ======= Tests ======= """


def calculate_stationary_distribution(rate_matrix) -> np.array:
    """
    Calculate the stationary distribution of a rate matrix.

    Args:
    rate_matrix (np.array): A square numpy array representing the rate matrix.

    Returns:
    np.array: The stationary distribution as a numpy array.
    """

    # Ensure the matrix is square
    if rate_matrix.shape[0] != rate_matrix.shape[1]:
        raise ValueError("Rate matrix must be square")

    print("Calculating stationary distribution for rate matrix:")

    # Find eigenvalues and eigenvectors
    eigenvalues, eigenvectors = np.linalg.eig(rate_matrix.T)

    # Find the eigenvector corresponding to the eigenvalue closest to zero
    stationary_vector = eigenvectors[:, np.isclose(eigenvalues, 0)].real

    # Normalize the stationary vector so its elements sum up to 1
    stationary_distribution = stationary_vector / stationary_vector.sum()

    return stationary_distribution.ravel()


def dict_to_alignment(sequence_dict):
    """
    Convert a dictionary of sequences to a MultipleSeqAlignment object.

    Args:
    sequence_dict (dict): A dictionary with sequence identifiers as keys and sequence strings as values.

    Returns:
    MultipleSeqAlignment: The corresponding MultipleSeqAlignment object.
    """
    alignment_list = []
    for id, sequence in sequence_dict.items():
        seq_record = SeqRecord(Seq(sequence), id=id)
        alignment_list.append(seq_record)

    return MultipleSeqAlignment(alignment_list)


class TestPartialLikelihoodCaching(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures before each test method."""
        # Create Nodes
        self.node_a = Node("A")
        self.node_b = Node("B")
        self.node_c = Node("C")

        # Connect nodes to form a simple graph
        self.node_a.connect(self.node_b, 0.1)
        self.node_b.connect(self.node_c, 0.2)

        # Create a simple RateMatrix for the test
        self.rate_matrix = RateMatrix(np.array([[-1, 1], [1, -1]]))

    def test_partial_likelihood_caching(self):
        """Test that partial_likelihood uses caching effectively."""
        # First call to partial_likelihood
        start_time = time.time()
        partial_likelihood(self.node_a, self.node_b, self.rate_matrix)
        first_call_duration = time.time() - start_time

        # Second call to partial_likelihood with the same parameters
        start_time = time.time()
        partial_likelihood(self.node_a, self.node_b, self.rate_matrix)
        second_call_duration = time.time() - start_time

        # Assert that the second call was faster, indicating caching
        self.assertLess(
            second_call_duration,
            first_call_duration,
            "Caching is not working as expected.",
        )


class TestPhylogeneticAnalysis(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures."""
        # Define your sequence records
        self.seq_records = [
            SeqRecord(Seq("SEKSQ"), id="t7"),
            SeqRecord(Seq("SEKSQ"), id="t3"),
            SeqRecord(Seq("SEKSQ"), id="t6"),
            SeqRecord(Seq("QQKTM"), id="t1"),
            SeqRecord(Seq("FSERL"), id="t2"),
            SeqRecord(Seq("SSRQQ"), id="t4"),
            SeqRecord(Seq("SSRQQ"), id="t5"),
        ]

        # Create a MultipleSeqAlignment object
        self.alignment = MultipleSeqAlignment(self.seq_records)

        # Define a newick string for your phylogenetic tree
        self.newick_string = "(t7:0.0000010000,t3:0.0000010000,(t6:0.0000010000,((t1:0.3068411287,t4:1.2358829187)Node4:0.4797066442,(t2:0.0000010000,t5:0.0000010000)Node5:0.5990502849)Node3:0.0291780183)Node2:0.0000010000)Node1;"

    def test_phylogenetic_analysis(self):
        """Perform a comprehensive test on phylogenetic analysis."""
        tree = Tree(self.newick_string, format=1)
        rename_internal_nodes_pre_order(tree)
        sequence_dict = {record.id: str(record.seq) for record in self.alignment}
        
        rate_matrix = RateMatrix(POISSON_RATE_MATRIX)
        psi_matrix = np.diag(AA_STATE_FREQUENCIES["POISSON"])
        array_left_eigenvectors, array_right_eigenvectors, multiplicity, eigenvalue = (
            spectral_decomposition(rate_matrix.rate_matrix, psi_matrix)
        )

        alignment = dict_to_alignment(sequence_dict)
        results = single_rate_analysis(
            tree,
            alignment,
            rate_matrix,
            np.array(AA_STATE_FREQUENCIES["POISSON"]),
            array_right_eigenvectors,
            multiplicity,
            0.05,
            None,
        )

        # Here, add assertions to validate your results as needed
        self.assertTrue(
            results
        )  # This is just a placeholder, customize it based on what `single_rate_analysis` returns


if __name__ == "__main__":
    unittest.main()
