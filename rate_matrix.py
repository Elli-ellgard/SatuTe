class RateMatrix:
    """Class representing a rate matrix for nucleotide substitutions."""

    def __init__(self, rate_matrix):
        """Initialize the RateMatrix with a given rate matrix."""
        self.rate_matrix = rate_matrix

    def __hash__(self):
        """Return the hash value of the RateMatrix."""
        return id(self)
