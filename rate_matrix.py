import numpy as np
from typing import Any


class RateMatrix:
    """
    A class representing a rate matrix for nucleotide substitutions.

    This class provides functionalities for handling rate matrices used in
    evolutionary biology and bioinformatics for modeling nucleotide substitutions.

    Attributes:
        rate_matrix (np.ndarray): A NumPy ndarray representing the substitution rate matrix.

    Methods:
        __init__(self, rate_matrix: np.ndarray): Initializes the RateMatrix with a given rate matrix.
        __hash__(self) -> int: Returns a hash value based on the content of the rate matrix.
        __eq__(self, other: Any) -> bool: Checks equality based on the content of the rate matrices.
    """

    def __init__(self, rate_matrix: np.ndarray) -> None:
        """
        Initializes the RateMatrix with a given rate matrix.

        Args:
            rate_matrix (np.ndarray): A NumPy ndarray representing the substitution rate matrix.
        """
        self.rate_matrix: np.ndarray = rate_matrix

    def __hash__(self) -> int:
        """
        Returns a hash value based on the content of the RateMatrix.

        The hash value is generated from the byte representation of the rate matrix,
        ensuring that identical matrices will produce the same hash value.

        Returns:
            int: The hash value of the rate matrix.
        """
        return hash(self.rate_matrix.tobytes())

    def __eq__(self, other: Any) -> bool:
        """
        Checks equality based on the content of the rate matrices.

        Determines if this RateMatrix is equal to another object. Equality is based
        on whether the other object is also a RateMatrix and has an identical rate matrix.

        Args:
            other (Any): The object to compare with this RateMatrix.

        Returns:
            bool: True if the other object is a RateMatrix and its rate matrix is equal to this RateMatrix's rate matrix; False otherwise.
        """
        if isinstance(other, RateMatrix):
            return np.array_equal(self.rate_matrix, other.rate_matrix)
        return False
