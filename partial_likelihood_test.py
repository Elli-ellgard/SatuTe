import numpy as np
from ete3 import Tree
from satute_categories import read_alignment_file
from satute_repository import IqTreeParser
from satute_util import spectral_decomposition
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from graph import Graph, Node, get_initial_likelihood_vector
from rate_matrix import RateMatrix
from partial_likelihood import (
    partial_likelihood,
    calculate_partial_likelihoods_for_sites,
    calculate_exponential_matrix,
)
import pprint
from scipy.sparse.linalg import expm

RATE_MATRIX = np.array([[-3, 1, 1, 1], [1, -3, 1, 1], [1, 1, -3, 1], [1, 1, 1, -3]])


def parse_newick_file(file_path):
    try:
        # Open the file in read mode
        with open(file_path, "r") as f:
            # Read the content of the file
            newick_string = f.readlines()

        # Parse the newick string into a Tree object
        t = Tree(newick_string[0], format=1)

        # Return the parsed Tree object
        return t

    except FileNotFoundError:
        raise Exception("File not found: " + file_path)


"""
======= Tests =======
"""


def name_nodes_by_level_order(tree):
    """Name nodes in a tree based on level-order traversal."""
    i = 1
    for node in tree.traverse("levelorder"):
        if not node.is_leaf():
            node.name = f"Node{i}"
            i += 1
    return tree


def test_one_partial_likelihood():
    a1, b2, u3, u4, c5, b6 = [
        Node("A", get_initial_likelihood_vector("A")),
        Node("B", get_initial_likelihood_vector("C")),
        Node("3"),
        Node("4"),
        Node("C", get_initial_likelihood_vector("A")),
        Node("D", get_initial_likelihood_vector("A")),
    ]

    tree = Graph(
        [(a1, u3, 0.01), (u3, b2, 0.01), (u3, u4, 0.01), (u4, c5, 0.01), (u4, b6, 0.01)]
    )

    for edge in tree.get_edges():
        left, right, branch_length = edge
        p1 = partial_likelihood(tree, left, right)
        p2 = partial_likelihood(tree, right, left)


def test_three():
    alignment_file = "./Clemens/example_3/sim-JC+G-AC1-AG1-AT1-CG1-CT1-GT1-alpha1.2-taxa64-len1000bp-bla0.01-blb0.8-blc0.2-rep01.fasta"
    newick_string = "(((S1:1,S2:2),(S3:1,S4:1),S5:1))"
    t = Tree(newick_string, format=1)
    t = name_nodes_by_level_order(t)
    alignment = read_alignment_file(alignment_file)
    rate_matrix = RateMatrix(RATE_MATRIX)
    partial_likelihood_per_site_storage = calculate_partial_likelihoods_for_sites(
        t, alignment, rate_matrix
    )


def test_two():

    alignment = read_alignment_file(
        "./test/Clemens/toy_example_JC/toy_example_ntaxa_7_run_5-alignment.phy"
    )

    t = parse_newick_file(
        "./test/Clemens/toy_example_JC/toy_example_ntaxa_7_run_5-alignment.phy.treefile"
    )

    t = name_nodes_by_level_order(t)

    iqtree_parser = IqTreeParser(
        "./test/Clemens/toy_example_JC/toy_example_ntaxa_7_run_5-alignment.phy.iqtree"
    )

    substitution_model = iqtree_parser.load_substitution_model()

    rate_matrix = RateMatrix(RATE_MATRIX)

    (
        array_left_eigenvectors,
        array_right_eigenvectors,
        multiplicity,
    ) = spectral_decomposition(
        substitution_model.rate_matrix, substitution_model.phi_matrix
    )


def print_beautifully(dictionary):
    pprinter = pprint.PrettyPrinter(indent=4)
    pprinter.pprint(dictionary)


def test_one():
    # Create SeqRecord objects for your sequences
    seq1 = SeqRecord(Seq("AA"), id="A")
    seq2 = SeqRecord(Seq("CG"), id="B")
    seq3 = SeqRecord(Seq("GG"), id="C")
    seq4 = SeqRecord(Seq("GG"), id="D")

    # Create a MultipleSeqAlignment object
    alignment = MultipleSeqAlignment([seq1, seq2, seq3, seq4])

    newick_string = "((A:0.2, B:0.4):0.3, (C:0.5,D:0.2):2);"

    t = Tree(newick_string, format=1)

    t = name_nodes_by_level_order(t)

    rate_matrix = RateMatrix(RATE_MATRIX)

    partial_likelihood_per_site_storage = calculate_partial_likelihoods_for_sites(
        t, alignment, rate_matrix
    )


def test_verify_partial_likelihoods_comprehensively():
    # Setup with predefined sequences and phylogenetic tree
    seq1 = SeqRecord(Seq("AA"), id="A")
    seq2 = SeqRecord(Seq("CG"), id="B")
    seq3 = SeqRecord(Seq("GG"), id="C")
    seq4 = SeqRecord(Seq("GG"), id="D")
    alignment = MultipleSeqAlignment([seq1, seq2, seq3, seq4])
    newick_string = "((A:0.2, B:0.4):0.3, (C:0.5,D:0.2):2);"
    t = Tree(newick_string, format=1)
    t = name_nodes_by_level_order(t)
    rate_matrix = RateMatrix(RATE_MATRIX)

    # Calculate partial likelihoods for sites
    partial_likelihood_per_site_storage = calculate_partial_likelihoods_for_sites(
        t, alignment, rate_matrix
    )

    # Verify overall structure and correctness of the output
    assert isinstance(
        partial_likelihood_per_site_storage, dict
    ), "Expected output to be a dictionary."
    assert (
        len(partial_likelihood_per_site_storage) == alignment.get_alignment_length()
    ), "Mismatch in expected and actual number of sites."

    # Detailed verification of likelihood values and structure
    for edge_key, edge_data in partial_likelihood_per_site_storage.items():
        for direction, data in edge_data.items():
            likelihoods = data["likelihoods"]
            for likelihood in likelihoods:
                assert set(likelihood.keys()) >= {
                    "Node",
                    "Site",
                    "branch_length",
                    "p_A",
                    "p_C",
                    "p_G",
                    "p_T",
                }, "Likelihood entry missing required keys."
                assert all(
                    likelihood[nucleotide] >= 0
                    for nucleotide in ["p_A", "p_C", "p_G", "p_T"]
                ), "Found negative likelihood value."

    # Asserting consistency for identical sequences
    for site_index in range(alignment.get_alignment_length()):
        # Ensure this logic matches your data structure; adjust as necessary
        likelihoods = (
            partial_likelihood_per_site_storage.get(
                f"some_key_for_site_{site_index}", {}
            )
            .get("some_direction", {})
            .get("likelihoods", [])
        )
        if likelihoods:  # If likelihoods for this site and direction exist
            p_values = {
                nucleotide: sum(
                    likelihood.get(nucleotide, 0) for likelihood in likelihoods
                )
                for nucleotide in ["p_A", "p_C", "p_G", "p_T"]
            }
            total_p = sum(p_values.values())
            # Example check: total likelihoods across nodes for a site could be compared to a model-specific expectation
            assert (
                total_p > 0
            ), f"Total likelihood for site {site_index} should be positive."

    print("Comprehensive verification of partial likelihoods passed successfully.")


def test_calculate_exponential_matrix():
    # Test 1: Basic functionality
    rate_matrix = RateMatrix(np.array([[-1, 1], [1, -1]]))
    branch_length = 1.0
    expected = expm(rate_matrix.rate_matrix * branch_length)
    np.testing.assert_array_almost_equal(
        calculate_exponential_matrix(rate_matrix, branch_length), expected
    )

    # Test 2: Zero branch length (should return identity matrix)
    branch_length_zero = 0.0
    identity_matrix = np.eye(rate_matrix.rate_matrix.shape[0])
    np.testing.assert_array_almost_equal(
        calculate_exponential_matrix(rate_matrix, branch_length_zero), identity_matrix
    )

    # Test 3: Negative branch length (should still work theoretically, but check your model's assumptions)
    branch_length_negative = -1.0
    expected_negative = expm(rate_matrix.rate_matrix * branch_length_negative)
    np.testing.assert_array_almost_equal(
        calculate_exponential_matrix(rate_matrix, branch_length_negative),
        expected_negative,
    )

    # Test 4: Caching - repeat a calculation to ensure it's using the cache
    import time

    start_time = time.time()
    calculate_exponential_matrix(
        rate_matrix, branch_length
    )  # Assume this calculation is expensive
    first_duration = time.time() - start_time

    start_time = time.time()
    calculate_exponential_matrix(
        rate_matrix, branch_length
    )  # This should be quicker due to caching
    second_duration = time.time() - start_time

    assert (
        second_duration < first_duration
    ), "Caching does not seem to be working correctly."


def test_partial_likelihoods_output_structure_and_length():
    # Include setup code here (same as in the first test)
    seq1 = SeqRecord(Seq("AA"), id="A")
    seq2 = SeqRecord(Seq("CG"), id="B")
    seq3 = SeqRecord(Seq("GG"), id="C")
    seq4 = SeqRecord(Seq("GG"), id="D")
    alignment = MultipleSeqAlignment([seq1, seq2, seq3, seq4])
    newick_string = "((A:0.2, B:0.4):0.3, (C:0.5,D:0.2):2);"
    t = Tree(newick_string, format=1)
    t = name_nodes_by_level_order(t)
    rate_matrix = RateMatrix(RATE_MATRIX)

    # Calculate partial likelihoods for sites
    partial_likelihood_per_site_storage = calculate_partial_likelihoods_for_sites(
        t, alignment, rate_matrix
    )

    # Assertions
    assert isinstance(
        partial_likelihood_per_site_storage, dict
    ), "Expected output to be a dictionary."

    # Check if the structure is as expected
    for edge_key, edge_data in partial_likelihood_per_site_storage.items():
        assert isinstance(edge_key, str), "Edge key should be a string."
        assert isinstance(edge_data, dict), "Edge data should be a dictionary."
        assert "left" in edge_data, "Left child node data is missing."
        assert "right" in edge_data, "Right child node data is missing."

        left_child = edge_data["left"]
        right_child = edge_data["right"]

        assert isinstance(
            left_child, dict
        ), "Left child node data should be a dictionary."
        assert isinstance(
            right_child, dict
        ), "Right child node data should be a dictionary."

        assert "likelihoods" in left_child, "Left child node likelihoods are missing."
        assert "likelihoods" in right_child, "Right child node likelihoods are missing."

        left_likelihoods = left_child["likelihoods"]
        right_likelihoods = right_child["likelihoods"]

        assert isinstance(left_likelihoods, list), "Left likelihoods should be a list."
        assert isinstance(
            right_likelihoods, list
        ), "Right likelihoods should be a list."

        assert len(left_likelihoods) == len(
            right_likelihoods
        ), "Mismatch in left and right likelihoods."

        for likelihood in left_likelihoods + right_likelihoods:
            assert isinstance(likelihood, dict), "Likelihood should be a dictionary."
            assert set(likelihood.keys()) == {
                "Node",
                "Site",
                "branch_length",
                "p_A",
                "p_C",
                "p_G",
                "p_T",
            }, "Likelihood keys are incorrect."

    # Check if the number of sites matches the alignment length
    num_sites = len(
        next(iter(partial_likelihood_per_site_storage.values()))["left"]["likelihoods"]
    )
    assert (
        num_sites == alignment.get_alignment_length()
    ), "Mismatch in expected and actual number of sites."


def test_partial_likelihoods_values_integrity():
    # Include setup code here (same as in the first test)
    seq1 = SeqRecord(Seq("AA"), id="A")
    seq2 = SeqRecord(Seq("CG"), id="B")
    seq3 = SeqRecord(Seq("GG"), id="C")
    seq4 = SeqRecord(Seq("GG"), id="D")
    alignment = MultipleSeqAlignment([seq1, seq2, seq3, seq4])
    newick_string = "((A:0.2, B:0.4):0.3, (C:0.5,D:0.2):2);"
    t = Tree(newick_string, format=1)
    t = name_nodes_by_level_order(t)
    rate_matrix = RateMatrix(RATE_MATRIX)

    # Calculate partial likelihoods for sites
    partial_likelihood_per_site_storage = calculate_partial_likelihoods_for_sites(
        t, alignment, rate_matrix
    )

    # Detailed verification of likelihood values and structure
    for edge_data in partial_likelihood_per_site_storage.values():
        for child_data in edge_data.values():
            for likelihood in child_data["likelihoods"]:
                assert all(
                    likelihood[nucleotide] >= 0
                    for nucleotide in ["p_A", "p_C", "p_G", "p_T"]
                ), "Found negative likelihood value."


def test_partial_likelihoods_consistency_for_identical_sequences():
    # Include setup code here (same as in the first test)
    seq1 = SeqRecord(Seq("AA"), id="A")
    seq2 = SeqRecord(Seq("CG"), id="B")
    seq3 = SeqRecord(Seq("GG"), id="C")
    seq4 = SeqRecord(Seq("GG"), id="D")
    alignment = MultipleSeqAlignment([seq1, seq2, seq3, seq4])
    newick_string = "((A:0.2, B:0.4):0.3, (C:0.5,D:0.2):2);"
    t = Tree(newick_string, format=1)
    t = name_nodes_by_level_order(t)
    rate_matrix = RateMatrix(RATE_MATRIX)

    # Calculate partial likelihoods for sites
    partial_likelihood_per_site_storage = calculate_partial_likelihoods_for_sites(
        t, alignment, rate_matrix
    )

    # Asserting consistency for identical sequences
    for site_index in range(alignment.get_alignment_length()):
        p_values = {"p_A": 0, "p_C": 0, "p_G": 0, "p_T": 0}
        for edge_data in partial_likelihood_per_site_storage.values():
            for child_data in edge_data.values():
                for likelihood in child_data["likelihoods"]:
                    if likelihood["Site"] == site_index:
                        p_values["p_A"] += likelihood["p_A"]
                        p_values["p_C"] += likelihood["p_C"]
                        p_values["p_G"] += likelihood["p_G"]
                        p_values["p_T"] += likelihood["p_T"]

        total_p = sum(p_values.values())
        assert (
            total_p > 0
        ), f"Total likelihood for site {site_index} should be positive."


if __name__ == "__main__":
    test_one()
    test_calculate_exponential_matrix()
    test_partial_likelihoods_consistency_for_identical_sequences()
    test_partial_likelihoods_values_integrity()
    test_partial_likelihoods_output_structure_and_length()
