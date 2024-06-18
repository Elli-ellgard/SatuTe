import numpy as np
from satute.decomposition import spectral_decomposition
from tests.unittests.partial_likelihood_test import calculate_stationary_distribution
from satute.rate_matrix import RateMatrix
from satute.amino_acid_models import (
    POISSON_RATE_MATRIX,
    AA_STATE_FREQUENCIES,
    create_rate_matrix_with_input,
    AMINO_ACID_RATE_MATRIX,
    print_matrix,
)

def test_spectral_decomposition():
    rate_matrix = RateMatrix(POISSON_RATE_MATRIX)
    state_frequencies = calculate_stationary_distribution(rate_matrix.rate_matrix)
    psi_matrix = np.diag(AA_STATE_FREQUENCIES["POISSON"])

    (
        array_left_eigenvectors,
        array_right_eigenvectors,
        multiplicity,
        eigenvalue,
    ) = spectral_decomposition(rate_matrix.rate_matrix, psi_matrix)


def construct_rate_matrix(rate_matrix):
    model_string = """1
                      1 1 
                      1 1 1
                    """
    rate_matrix = create_rate_matrix_with_input(4, model_string)

    f = open("amino_acid_models.txt", "w")

    for key, string_matrix in AMINO_ACID_RATE_MATRIX.items():
        matrix = create_rate_matrix_with_input(20, string_matrix)
        f.write(key)
        f.write("\n")
        print_matrix(matrix, f)
        stationary_distribution = AA_STATE_FREQUENCIES[key]
        phi = np.diag(stationary_distribution)
        left, right, m, eigenvalue = spectral_decomposition(np.array(matrix), phi)

        # Custom printing function
        def print_eigenvectors(vectors):
            for vec in vectors:
                print(np.array_str(vec, precision=3, suppress_small=True))

        print(key)
        # Print left and right eigenvectors
        print("Left Eigenvectors:")
        print_eigenvectors(left)

        print("Right Eigenvectors:")
        print_eigenvectors(right)

        print("Multiplicity:", m)

        f.write("\n")
