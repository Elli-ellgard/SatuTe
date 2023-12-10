import numpy as np
from ete3 import Tree
from satute_categories import read_alignment_file
from satute_trees_and_subtrees import parse_newick_file
from satute_repository import parse_rate_matrices_from_file
from satute_util import spectral_decomposition
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from graph import Graph, Node
from rate_matrix import RateMatrix

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
    alignment_file = "./Clemens/example_3/sim-JC+G-AC1-AG1-AT1-CG1-CT1-GT1-alpha1.2-taxa64-len1000bp-bla0.01-blb0.8-blc0.2-rep01.fasta"
    alignment = read_alignment_file(alignment_file)
    t = parse_newick_file(
        "./Clemens/example_3/sim-JC+G-AC1-AG1-AT1-CG1-CT1-GT1-alpha1.2-taxa64-len1000bp-bla0.01-blb0.8-blc0.2-rep01.fasta.treefile"
    )
    t = name_nodes_by_level_order(t)

    RATE_MATRIX, psi_matrix = parse_rate_matrices_from_file(
        4,
        "./Clemens/example_3/sim-JC+G-AC1-AG1-AT1-CG1-CT1-GT1-alpha1.2-taxa64-len1000bp-bla0.01-blb0.8-blc0.2-rep01.fasta",
    )
    rate_matrix = RateMatrix(RATE_MATRIX)

    (
        array_left_eigenvectors,
        array_right_eigenvectors,
        multiplicity,
    ) = spectral_decomposition(RATE_MATRIX, psi_matrix)


def test_one():
    # Create SeqRecord objects for your sequences
    seq1 = SeqRecord(Seq("AGTATA"), id="A")
    seq2 = SeqRecord(Seq("CGTATG"), id="B")
    seq3 = SeqRecord(Seq("GGTATG"), id="C")
    seq4 = SeqRecord(Seq("GGTACG"), id="D")

    # Create a MultipleSeqAlignment object
    alignment = MultipleSeqAlignment([seq1, seq2, seq3, seq4])
    newick_string = "((A:0.2, B:0.4):0.3, (C:0.5,D:0.2):2);"
    t = Tree(newick_string, format=1)
    t = name_nodes_by_level_order(t)

    print(t.get_ascii(attributes=["name", "dist"]))

    rate_matrix = RateMatrix(RATE_MATRIX)

    partial_likelihood_per_site_storage = calculate_partial_likelihoods_for_sites(
        t, alignment, rate_matrix
    )
