from tests.scripts.fixtures import *
from tests.scripts.satute_test_utils import find_file_with_suffix

import numpy as np
from satute.spectral_decomposition import spectral_decomposition
from satute.partial_likelihood.rate_matrix import RateMatrix
from satute.parser.iqtree_parser import IqTreeParser


''' DNA MODELS '''

def test_JC_model(dir_path_iqtree_files):
    model="JC"
    source_path  = dir_path_iqtree_files[0]
    iqtree_file_path = find_file_with_suffix("example_dna_model", f"{model}.iqtree",source_path)

    satute_iq_tree_parser = IqTreeParser(iqtree_file_path)
    substitution_model = satute_iq_tree_parser.load_substitution_model()
    rate_matrix = RateMatrix(substitution_model.rate_matrix)

    (
        array_left_eigenvectors,
        array_right_eigenvectors,
        multiplicity,
        eigenvalue,
    ) = spectral_decomposition(substitution_model.rate_matrix, substitution_model.phi_matrix)

    assert multiplicity == 3, "The largest non-zero eigenvalue of JC should have multiplicity 3."
    assert len(array_right_eigenvectors) == multiplicity, "Dimension of right eigenspace is incorrect."
    assert len(array_left_eigenvectors) == multiplicity, "Dimension of left eigenspace is incorrect."
    

def test_K2P_model(dir_path_iqtree_files):
    model="K2P"
    source_path  = dir_path_iqtree_files[0]
    iqtree_file_path = find_file_with_suffix("example_dna_model", f"{model}.iqtree",source_path)

    satute_iq_tree_parser = IqTreeParser(iqtree_file_path)
    substitution_model = satute_iq_tree_parser.load_substitution_model()
    rate_matrix = RateMatrix(substitution_model.rate_matrix)

    (
        array_left_eigenvectors,
        array_right_eigenvectors,
        multiplicity,
        eigenvalue,
    ) = spectral_decomposition(substitution_model.rate_matrix, substitution_model.phi_matrix)

    assert multiplicity == 1, "The largest non-zero eigenvalue of this K2P model should have multiplicity 1."
    assert len(array_right_eigenvectors) == multiplicity, "Dimension of right eigenspace is incorrect."
    assert len(array_left_eigenvectors) == multiplicity, "Dimension of left eigenspace is incorrect."

def test_K2Pb_model(dir_path_iqtree_files):
    model="K2Pb"
    source_path  = dir_path_iqtree_files[0]
    iqtree_file_path = find_file_with_suffix("example_dna_model", f"{model}.iqtree",source_path)

    satute_iq_tree_parser = IqTreeParser(iqtree_file_path)
    substitution_model = satute_iq_tree_parser.load_substitution_model()
    rate_matrix = RateMatrix(substitution_model.rate_matrix)

    (
        array_left_eigenvectors,
        array_right_eigenvectors,
        multiplicity,
        eigenvalue,
    ) = spectral_decomposition(substitution_model.rate_matrix, substitution_model.phi_matrix)

    assert multiplicity == 2, "The largest non-zero eigenvalue of this K2P model should have multiplicity 2."
    assert len(array_right_eigenvectors) == multiplicity, "Dimension of right eigenspace is incorrect."
    assert len(array_left_eigenvectors) == multiplicity, "Dimension of left eigenspace is incorrect."




def test_GTR_model(dir_path_iqtree_files):
    model="GTR_F"
    source_path  = dir_path_iqtree_files[0]
    iqtree_file_path = find_file_with_suffix("example_dna_model", f"{model}.iqtree",source_path)

    satute_iq_tree_parser = IqTreeParser(iqtree_file_path)
    substitution_model = satute_iq_tree_parser.load_substitution_model()
    rate_matrix = RateMatrix(substitution_model.rate_matrix)

    (
        array_left_eigenvectors,
        array_right_eigenvectors,
        multiplicity,
        eigenvalue,
    ) = spectral_decomposition(substitution_model.rate_matrix, substitution_model.phi_matrix)

    assert multiplicity == 1, "The largest non-zero eigenvalue of GTR should have multiplicity 1."
    assert len(array_right_eigenvectors) == multiplicity, "Dimension of right eigenspace is incorrect."
    assert len(array_left_eigenvectors) == multiplicity, "Dimension of left eigenspace is incorrect."





''' PROTEIN MODELS '''


def test_LG_model(dir_path_iqtree_files):
    model="LG"
    source_path  = dir_path_iqtree_files[0]
    iqtree_file_path = find_file_with_suffix("example_aa_model", f"{model}.iqtree",source_path)

    satute_iq_tree_parser = IqTreeParser(iqtree_file_path)
    substitution_model = satute_iq_tree_parser.load_substitution_model()
    rate_matrix = RateMatrix(substitution_model.rate_matrix)

    (
        array_left_eigenvectors,
        array_right_eigenvectors,
        multiplicity,
        eigenvalue,
    ) = spectral_decomposition(substitution_model.rate_matrix, substitution_model.phi_matrix)

    assert multiplicity == 1, "The largest non-zero eigenvalue of LG should have multiplicity 1."
    assert len(array_right_eigenvectors) == multiplicity, "Dimension of right eigenspace is incorrect."
    assert len(array_left_eigenvectors) == multiplicity, "Dimension of left eigenspace is incorrect."



def test_poisson_model():
    from satute.models.amino_acid_models import POISSON_RATE_MATRIX, AA_STATE_FREQUENCIES
    rate_matrix = RateMatrix(POISSON_RATE_MATRIX)
    psi_matrix = np.diag(AA_STATE_FREQUENCIES["POISSON"])

    (
        array_left_eigenvectors,
        array_right_eigenvectors,
        multiplicity,
        eigenvalue,
    ) = spectral_decomposition(rate_matrix.rate_matrix, psi_matrix)

    assert multiplicity == 19, "The largest non-zero eigenvalue of POISSON should have multiplicity 19."
    assert len(array_right_eigenvectors) == multiplicity, "Dimension of right eigenspace is incorrect."
    assert len(array_left_eigenvectors) == multiplicity, "Dimension of left eigenspace is incorrect."
