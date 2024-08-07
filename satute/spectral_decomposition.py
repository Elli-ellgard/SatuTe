# -*- coding: utf-8 -*-
import scipy
import numpy as np


""" ## SPECTRAL DECOMPOSITION OF THE RATE MATRIX"""


def spectral_decomposition(
    rate_matrix, psi_matrix
) -> tuple[np.array, np.array, int, int]:
    """
    Then psi_matrix := Diag(pi). Recall that matrix Q is reversible iff M:= psi_matrix^1/2 x Q x psi_matrix^{-1/2} is symmetric.
    For a real symmetric matrix M, its eigenvectors can be chosen to be an orthonormal basis of R^n
    """
    M = scipy.linalg.fractional_matrix_power(psi_matrix, +1 / 2) @ rate_matrix
    M = M @ scipy.linalg.fractional_matrix_power(psi_matrix, -1 / 2)

    """ Schur decomposition of matrix M"""
    lamb_matrix, w = scipy.linalg.schur(M)  # Compute the eigenvalues and eigenvectors.
    lamb = np.diagonal(lamb_matrix)
    # lamb, w = np.linalg.eig(M)  # Compute the eigenvalues and eigenvectors.
    idx = lamb.argsort()[::-1]  # Order from small to large.
    lamb = lamb[idx]  # Order the eigenvalues according to idx.
    w = w[:, idx]  # Order the eigenvectors according to the eigenvalues"""
    # the first one should be the eigenvalue 0 in lamb, why are we doing the following?
    lamb_nozero = []  # list of eigenvalues without 0
    for i in lamb:
        if i > 0.00999 or i < -0.00999:
            lamb_nozero.append(i)

    max_lambda = max(lamb_nozero)  # dominant non-zero eigenvalue
    # get the indices of the dominant non-zero eigenvalue in lamb taking numerical inaccuracies into account and identical values
    index = []
    for i in range(len(lamb)):
        lambda_it = lamb[i]
        if abs(lambda_it - max_lambda) < 0.01:
            index.append(i)

    multiplicity = len(index)  # multiplicity of the dominant non-zero eigenvalue
    array_right_eigenvectors = []  # list of right eigenvectors for the dominant non-zero eigenvalue
    array_left_eigenvectors = []  # list of left eigenvectors for the dominant non-zero eigenvalue
    for i in range(multiplicity):
        # calculate the right eigenvectors for the dominant non-zero eigenvalue
        v1 = scipy.linalg.fractional_matrix_power(psi_matrix, -1 / 2) @ w[:, index[i]]
        array_right_eigenvectors.append(v1)
        # calculate the left eigenvectors for the dominant non-zero eigenvalue
        h1 = scipy.linalg.fractional_matrix_power(psi_matrix, +1 / 2) @ w[:, index[i]]
        array_left_eigenvectors.append(h1)

    return array_left_eigenvectors, array_right_eigenvectors, multiplicity, max_lambda
