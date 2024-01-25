import numpy as np
from ete3 import Tree
import sys

sys.path.append("..")

from satute_categories import read_alignment_file
from satute_trees import (
    rename_internal_nodes_preorder,
    collapse_identical_leaf_sequences,
)
from file_handler import FileHandler
from satute_repository import IqTreeParser
from satute_util import spectral_decomposition
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from graph import Graph, Node
from rate_matrix import RateMatrix
from amino_acid_models import POISSON_RATE_MATRIX, AA_STATE_FREQUENCIES
from satute_rate_analysis import single_rate_analysis
import numpy as np
import pandas as pd
from graph import get_initial_likelihood_vector
from partial_likelihood import (
    partial_likelihood,
    calculate_partial_likelihoods_for_sites,
)


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
    file_handler = FileHandler()
    newick_string = file_handler.get_newick_string_from_file(
        "./Clemens/example_3/sim-JC+G-AC1-AG1-AT1-CG1-CT1-GT1-alpha1.2-taxa64-len1000bp-bla0.01-blb0.8-blc0.2-rep01.fasta.treefile"
    )
    t = Tree(newick_string, format=1)

    t = name_nodes_by_level_order(t)
    iq_tree_parser = IqTreeParser()

    RATE_MATRIX, psi_matrix = iq_tree_parser.parse_rate_matrices(
        4,
        "./Clemens/example_3/sim-JC+G-AC1-AG1-AT1-CG1-CT1-GT1-alpha1.2-taxa64-len1000bp-bla0.01-blb0.8-blc0.2-rep01.fasta.iqtree",
    )
    rate_matrix = RateMatrix(RATE_MATRIX)

    (
        array_left_eigenvectors,
        array_right_eigenvectors,
        multiplicity,
    ) = spectral_decomposition(RATE_MATRIX, psi_matrix)


def calculate_stationary_distribution(rate_matrix)->np.array:
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
    
    # print(eigenvalues)
    print(eigenvalues)
        
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


def test_one():
    # Step 1: Create SeqRecord objects for your sequences
    seq_records = [
        SeqRecord(Seq("ACGTAT"), id="ACGTAT_1"),
        SeqRecord(Seq("ACGTAT"), id="ACGTAT_2"),
        SeqRecord(Seq("ACGTAT"), id="ACGTAT_3"),
        SeqRecord(Seq("GGTATG"), id="C"),
        SeqRecord(Seq("GGTATG"), id="E"),
        SeqRecord(Seq("GGTACG"), id="D"),
    ]

    # Step 2: Create a MultipleSeqAlignment object from the SeqRecord objects
    alignment = MultipleSeqAlignment(seq_records)

    # Step 3: Create a phylogenetic tree from a Newick string
    newick_string = (
        "(((ACGTAT_1:0.2, ACGTAT_2:0.4):0.3,ACGTAT_3:1):1,(C:0.5,D:0.2, E:0.1):2);"
    )
    tree = Tree(newick_string, format=1)

    print(tree.write(format=1, format_root_node=True))

    # Step 4: Rename internal nodes in a preorder traversal
    rename_internal_nodes_preorder(tree)
    print(tree.get_ascii(show_internal=True))

    # Step 5: Create a dictionary to access sequences by their ID
    sequence_dict = {record.id: str(record.seq) for record in alignment}

    # Creating a deep copy of the tree for manipulation
    collapsed_tree_one = tree.copy("deepcopy")

    # Step 6: Process the tree to collapse nodes with identical sequences
    sequence_dict, twin_dictionary = collapse_identical_leaf_sequences(
        collapsed_tree_one, sequence_dict
    )

    print(collapsed_tree_one.write(format=1, format_root_node=True))

    # Step 7: Create a rate matrix and calculate stationary distribution
    rate_matrix = RateMatrix(RATE_MATRIX)
    
    
    state_frequencies = calculate_stationary_distribution(rate_matrix.rate_matrix)
    phi_matrix = np.diag(state_frequencies)

    # Step 8: Perform spectral decomposition
    (
        array_left_eigenvectors,
        array_right_eigenvectors,
        multiplicity,
    ) = spectral_decomposition(rate_matrix.rate_matrix, phi_matrix)

    print(collapsed_tree_one.get_ascii(show_internal=True))

    # Step 9: Convert the sequence dictionary back to a MultipleSeqAlignment object
    alignment = dict_to_alignment(sequence_dict)

    # Step 10: Print the rate matrix
    print("Rate Matrix:\n", rate_matrix.rate_matrix)

    # Step 11: Perform single rate analysis
    results = single_rate_analysis(
        collapsed_tree_one,
        alignment,
        rate_matrix,
        state_frequencies,
        array_left_eigenvectors,
        multiplicity,
        0.05,
        None,
    )

    # Step 12: Append additional data to results for twin nodes
    # for parent, value in twin_dictionary.items():
    #     for child in value:
    #         results["single_rate"].append(
    #            {
    #                 "edge": f"({parent}, {child})",
    #                 "delta": "Nan",
    #                 "c_s": "Nan",
    #                 "c_sTwoSequence": "Nan",
    #                 "p_value": "Nan",
    #                 "result_test": "Nan",
    #                 "result_test_tip2tip": "Nan",
    #                 "category_rate": 1,
    #                 "branch_length": tree.get_distance(parent, child),
    #             }
    #        )
    # Step 13: Convert results to a pandas DataFrame
    # pandas_data_frame = pd.DataFrame.from_dict(results["single_rate"])
    # pandas_data_frame.to_csv("satute_test_one.csv")


def test_two():
    # Step 1: Create SeqRecord objects for your sequences
    seq_records = [
        SeqRecord(Seq("SEKSQCLIGVAHVSVNATIDYLRQDPTAHLMGDLYQGWIMESKAKP"), id="t7"),
        SeqRecord(Seq("SEKSQCLIGVAHVSVNATIDYLRQDPTAHLMGDLYQGWIMESKAKP"), id="t3"),
        SeqRecord(Seq("SEKSQCLIGVAHVSVNATIDYLRQDPTAHLMGDLYQGWIMESKAKP"), id="t6"),
        SeqRecord(Seq("QQKTMALMAIPHVSANAAHEYLKSDVTKQRMGDPFQGQFVESKGKS"), id="t1"),
        SeqRecord(Seq("FSERLAMNGFIRITVNIALEYIRSYPIAQRVGDSFKGDLTESKSIF"), id="t2"),
        SeqRecord(Seq("SSRQQFLLAIAHVGINAATNYIRCDSTAYVTGDLYQGAVMEVQEKP"), id="t4"),
        SeqRecord(Seq("SSRQQFLLAIAHVGINAATNYIRCDSTAYVTGDLYQGAVMEVQEKP"), id="t5"),
    ]

    # Step 2: Create a MultipleSeqAlignment object from the SeqRecord objects
    alignment = MultipleSeqAlignment(seq_records)

    # Step 3: Create a phylogenetic tree from a Newick string
    newick_string = "(t7:0.0000010000,t3:0.0000010000,(t6:0.0000010000,((t1:0.3068411287,t4:1.2358829187)Node4:0.4797066442,(t2:0.0000010000,t5:0.0000010000)Node5:0.5990502849)Node3:0.0291780183)Node2:0.0000010000)Node1;"
    tree = Tree(newick_string, format=1)

    # Step 4: Rename internal nodes in a preorder traversal
    rename_internal_nodes_preorder(tree)
    print(tree.get_ascii(show_internal=True))

    # Step 5: Create a dictionary to access sequences by their ID
    sequence_dict = {record.id: str(record.seq) for record in alignment}

    # Creating a deep copy of the tree for manipulation
    collapsed_tree_one = tree.copy("deepcopy")

    # Step 6: Process the tree to collapse nodes with identical sequences
    sequence_dict, twin_dictionary = collapse_identical_leaf_sequences(
        collapsed_tree_one, sequence_dict
    )

    # Step 7: Create a rate matrix and calculate stationary distribution
    rate_matrix = RateMatrix(POISSON_RATE_MATRIX)
    
    state_frequencies = AA_STATE_FREQUENCIES['POISSON']
    phi_matrix = np.diag(state_frequencies)    
        
        
    # print(POISSON_RATE_MATRIX)        
    # state_frequencies = calculate_stationary_distribution(rate_matrix.rate_matrix)
    # phi_matrix = np.diag(state_frequencies)
    # 
    # Step 8: Perform spectral decomposition
    # stationary_distribution = calculate_stationary_distribution(rate_matrix.rate_matrix)    

    
    (
        array_left_eigenvectors,
        array_right_eigenvectors,
        multiplicity,
    ) = spectral_decomposition(rate_matrix.rate_matrix, phi_matrix)
    
    print(array_left_eigenvectors, array_right_eigenvectors, multiplicity)    


    # Step 9: Convert the sequence dictionary back to a MultipleSeqAlignment object
    alignment = dict_to_alignment(sequence_dict)

    """
    # Calculate partial likelihoods for all sites
    partial_likelihood_per_site_storage = calculate_partial_likelihoods_for_sites(
        tree=collapsed_tree_one,
        alignment=alignment,
        rate_matrix=rate_matrix,
        focused_edge=None,
    )
    """
    
    # print(partial_likelihood_per_site_storage)
    
    
    # Step 11: Perform single rate analysis
    # results = single_rate_analysis(
    #     collapsed_tree_one,
    #     alignment,
    #     rate_matrix,
    #     state_frequencies,
    #     array_left_eigenvectors,
    #     array_right_eigenvectors,
    #     multiplicity,
    #     0.05,
    #     None,
    # )
    # Step 12: Append additional data to results for twin nodes
    # for parent, value in twin_dictionary.items():
    #     for child in value:
    #         results["single_rate"].append(
    #             {
    #                 "edge": f"({parent}, {child})",
    #                 "delta": "Nan",
    #                 "c_s": "Nan",
    #                 "c_sTwoSequence": "Nan",
    #                 "p_value": "Nan",
    #                 "result_test": "Nan",
    #                 "result_test_tip2tip": "Nan",
    #                 "category_rate": 1,
    #                 "branch_length": tree.get_distance(parent, child),
    #             }
    #         )
    #
    # Step 13: Convert results to a pandas DataFrame
    # pandas_data_frame = pd.DataFrame.from_dict(results["single_rate"])
    # pandas_data_frame.to_csv("satute_test_one.csv")


if __name__ == "__main__":
    # test_one()
    test_two()
