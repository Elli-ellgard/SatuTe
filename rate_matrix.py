import numpy as np

class RateMatrix:
    """Class representing a rate matrix for nucleotide substitutions."""

    def __init__(self, rate_matrix):
        """Initialize the RateMatrix with a given rate matrix."""
        self.rate_matrix = rate_matrix

    def __hash__(self):
        """Return a hash value based on the content of the RateMatrix."""
        return hash(self.rate_matrix.tostring())

    def __eq__(self, other):
        """Check equality based on the content of the rate matrices."""
        return np.array_equal(self.rate_matrix, other.rate_matrix)
