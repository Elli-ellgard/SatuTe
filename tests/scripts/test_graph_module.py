import numpy as np
from ete3 import Tree
import unittest
from satute.graph import get_initial_likelihood_vector, convert_tree_to_graph

def test_get_initial_likelihood_vector_returns_numpy_array():
    # Test with a nucleotide state
    nucleotide_state = "A"
    nucleotide_likelihood_vector = get_initial_likelihood_vector(nucleotide_state, "nucleotide")
    assert isinstance(nucleotide_likelihood_vector, np.ndarray), "Expected a numpy array for nucleotide state"

    # Test with an amino acid state
    amino_acid_state = "A"
    amino_acid_likelihood_vector = get_initial_likelihood_vector(amino_acid_state, "amino_acid")
    assert isinstance(amino_acid_likelihood_vector, np.ndarray), "Expected a numpy array for amino acid state"

def test_get_initial_likelihood_vector_with_unknown_nucleotide_state():
    state = "-"
    state_type = "nucleotide"
    nucleotide_likelihood_vector = get_initial_likelihood_vector(state,state_type)
    assert isinstance(nucleotide_likelihood_vector, np.ndarray), "Expected a numpy array for amino acid state"

def test_get_initial_likelihood_vector_with_unknown_aa_state():
    state = "-"
    state_type = "amino_acid"
    nucleotide_likelihood_vector = get_initial_likelihood_vector(state,state_type)
    assert isinstance(nucleotide_likelihood_vector, np.ndarray), "Expected a numpy array for amino acid state"

class TestConvertTreeToGraph(unittest.TestCase):
    
    def test_single_node_tree(self):
        # Test conversion of a tree with a single node
        tree = Tree("A;")
        graph = convert_tree_to_graph(tree)
        self.assertEqual(len(graph.get_edges()), 0)
        self.assertEqual(len(tree.get_leaves()), 1)
        self.assertEqual(graph.get_edges(), [])

    def test_simple_tree(self):
        # Test conversion of a simple tree with three nodes
        newick = "((A:1,B:1)Node1:1,(C:1,D:1)Node2:1);"
        tree = Tree(newick, format=1)
        graph = convert_tree_to_graph(tree)
        # Ensure the graph has the correct number of edges
        self.assertEqual(len(graph.get_edges()), 6)
        # Check if the nodes are correctly connected and branch lengths are accurate
        expected_edges = set([("Node1", "A", 1), ("Node1", "B", 1), ("Node2", "C", 1), ("Node2", "D", 1), ("", "Node1", 1), ("", "Node2", 1)])
        result_edges = set([(edge[0].name, edge[1].name, edge[2]) for edge in graph.get_edges()])                
        self.assertEqual(result_edges, expected_edges)

    def test_six_taxa_tree(self):
        # Test conversion of a simple tree with three nodes
        newick = "((A:2,B:3)Node1:1,(C:1,D:1)Node2:1, (E:1,F:1)Node3);"
        tree = Tree(newick, format=1)
        graph = convert_tree_to_graph(tree)
        # Ensure the graph has the correct number of edges
        self.assertEqual(len(graph.get_edges()), 9)
        # Check if the nodes are correctly connected and branch lengths are accurate
        expected_edges = set([("Node1", "A", 2), ("Node1", "B", 3), ("Node2", "C", 1), ("Node2", "D", 1), ("", "Node1", 1), ("", "Node2", 1), ("Node3", "F", 1),("Node3", "E", 1), ("", "Node3", 1), ("", "Node3", 1)])
        result_edges = set([(edge[0].name, edge[1].name, edge[2]) for edge in graph.get_edges()])                
        self.assertEqual(result_edges, expected_edges)

    def test_branch_lengths(self):
        # Test tree with varying branch lengths
        newick = "((A:1,B:1)Node1:1,(C:1,D:1)Node2:1);"
        tree = Tree(newick, format=1)
        graph = convert_tree_to_graph(tree)
        
        # Validate branch lengths
        edges = graph.get_edges()
        print(edges)
        self.assertAlmostEqual(edges[0][2], 1)
        self.assertAlmostEqual(edges[1][2], 1)

    def test_complex_tree(self):
        # Test conversion of a more complex tree
        newick = "((A:1,B:1)Node1:2,(C:3,E:4)Node2:5)Root;"
        tree = Tree(newick, format=1)
        graph = convert_tree_to_graph(tree)
        
        # Ensure all nodes and edges are present
        self.assertEqual(len(graph.get_edges()), 6)
        
        # Optional: check more details like node connectivity and exact edge setup
        # This can be detailed based on specific connectivity checks among nodes

if __name__ == '__main__':
    unittest.main()
    