import sys

sys.path.append("..")

import numpy as np
from satute_util import spectral_decomposition
from rate_matrix import RateMatrix
from amino_acid_models import POISSON_RATE_MATRIX, AA_STATE_FREQUENCIES
from tests.unittests.partial_likelihood_test import calculate_stationary_distribution
from partial_likelihood import (
    partial_likelihood,
    calculate_partial_likelihoods_for_sites,
    calculate_exponential_matrix,
)


def test_spectral_decomposition():
    rate_matrix = RateMatrix(POISSON_RATE_MATRIX)
    state_frequencies = calculate_stationary_distribution(rate_matrix.rate_matrix)
    psi_matrix = np.diag(AA_STATE_FREQUENCIES["POISSON"])

    (
        array_left_eigenvectors,
        array_right_eigenvectors,
        multiplicity,
    ) = spectral_decomposition(rate_matrix.rate_matrix, psi_matrix)

    print(array_left_eigenvectors, array_right_eigenvectors, multiplicity)


test_spectral_decomposition()
