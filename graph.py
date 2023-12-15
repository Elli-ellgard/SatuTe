from typing import Optional
import dataclasses
import numpy as np
from ete3 import Tree
from nucleotide_code_vector import NUCLEOTIDE_CODE_VECTOR


class Graph:
    """Class representing a directed acyclic graph (DAG) of interconnected nodes."""

    def __init__(self, edges):
        """Initialize the DAG with a list of directed edges."""
        for left, right, branch_length in edges:
            left.connect(right, branch_length)
            right.connect(left, branch_length)
        self.edges = edges

    def get_edges(self):
        """Return the list of edges in the DAG."""
        return self.edges

    def set_edges(self, edges):
        self.edges = edges

    def get_branch_length(self, edge):
        """Return the branch length of a given edge in the DAG."""
        return self.branch_lengths[edge]


class Node:
    """Class representing a node in the DAG."""

    name: str
    state: Optional[np.array] = dataclasses.field(default_factory=lambda: None)
    connected: dict = dataclasses.field(default_factory=dict)

    def __init__(self, name, state=None, connected=None):
        """Initialize a Node with a name, optional state vector, and optional connections."""
        self.name = name
        self.state = state
        self.connected = connected or {}

    def __hash__(self):
        """Return the hash value of the Node."""
        return id(self)

    def __eq__(self, other):
        """Check if two nodes are equal."""
        return self is other

    def connect(self, other_node, branch_length):
        """Connect the current node to another node with a given branch length."""
        self.connected[other_node] = branch_length

    def is_leaf(self):
        """Check if the node is a leaf (has a state vector)."""
        return self.state is not None

    def __repr__(self):
        """Return the string representation of the Node."""
        return f"Node({self.name})"


def get_initial_likelihood_vector(state):
    """Get the initial likelihood vector for a given nucleotide state."""
    return np.array(NUCLEOTIDE_CODE_VECTOR[state]).T


def convert_ete3_tree_to_directed_acyclic_graph(
    tree, msa_column, alignment_look_up_table
):
    """
    Convert an ETE3 tree to a DirectedAcyclicGraph representation.

    This function transforms an ETE3 tree into a directed acyclic graph (DAG) representation
    by creating Node objects for each node in the tree and establishing connections between
    them based on the tree structure. Leaf nodes are initialized with initial likelihood vectors
    derived from the sequence alignment data.

    Args:
        tree (Tree): An ETE3 tree object representing the phylogenetic tree.
        msa_column (np.array): A single column of the multiple sequence alignment for a site.
        alignment_look_up_table (dict): A lookup table mapping sequence IDs to alignment indices.

    Returns:
        DirectedAcyclicGraph: A directed acyclic graph representation of the phylogenetic tree.

    Note:
        - The ETE3 tree should be properly constructed and rooted.
        - The msa_column parameter provides data for a specific site in the alignment.
        - The alignment_look_up_table aids quick retrieval of alignment record indices.
        - The function creates Node and DirectedAcyclicGraph objects to represent the graph structure.
    """
    node_dictionary = {}

    # Create nodes for each node in the ETE3 tree using a level-order traversal.
    for node in tree.traverse("levelorder"):
        if node.is_leaf():
            # Initialize leaf nodes with initial likelihood vectors from the alignment.
            node_dictionary[node.name] = Node(
                node.name,
                get_initial_likelihood_vector(
                    msa_column[alignment_look_up_table[node.name]].seq
                ),
            )
        else:
            # Non-leaf nodes are initialized with empty states.
            node_dictionary[node.name] = Node(node.name)

    edge_list = []
    # Create edges between nodes based on the tree's hierarchical relationships.
    for node in tree.traverse("levelorder"):
        for child_node in node.children:
            edge_list.append(
                (
                    node_dictionary[node.name],
                    node_dictionary[child_node.name],
                    child_node.dist,
                )
            )

    # Return a DirectedAcyclicGraph representation of the tree using the constructed edges.
    return Graph(edge_list)


def count_and_nodes_branches_nodes(tree, node, coming_from):
    # If the current node is a leaf
    if node.is_leaf():
        return 1, 0

    # Initialize leaf and branch counts
    leaf_count, branch_count = 0, 0

    # Iterate over children
    for child in node.connected.keys():
        if child != coming_from:
            # Recursively count leaves and branches for each child
            child_leaves, child_branches = count_and_nodes_branches_nodes(
                tree, child, node
            )
            leaf_count += child_leaves
            branch_count += child_branches + 1  # Include the branch to this child
    # Print the count for the current node
    return leaf_count, branch_count


def count_leaves_and_branches_for_subtrees(tree: Tree, alignment, focused_edges=None):
    # Create a lookup table for alignment
    alignment_look_up_table = get_alignment_look_up_table(alignment)

    # Convert the given tree to a directed acyclic graph
    count_graph = convert_ete3_tree_to_directed_acyclic_graph(
        tree, alignment[:, 1:2], alignment_look_up_table
    )

    # Filter the graph edges if focused edges are provided
    if focused_edges:
        filter_graph_edges_by_focus(count_graph, focused_edges)

    # Initialize a dictionary to store counts
    edge_subtree_count_dict = {}

    # Iterate over the edges of the graph
    for edge in count_graph.get_edges():
        right, left, branch_length = edge

        # Count nodes and branches on the left side
        (
            count_graph_nodes_left,
            count_graph_branches_left,
        ) = count_and_nodes_branches_nodes(tree, left, right)
        # Count nodes and branches on the right side
        (
            count_graph_nodes_right,
            count_graph_branches_right,
        ) = count_and_nodes_branches_nodes(tree, right, left)
        # Formulate the edge name and store the counts in the dictionary
        edge_name = f"({left.name}, {right.name})"
        edge_subtree_count_dict[edge_name] = {
            "left": {
                "leave_count": count_graph_nodes_left,
                "branch_count": count_graph_branches_left,
            },
            "right": {
                "leave_count": count_graph_nodes_right,
                "branch_count": count_graph_branches_right,
            },
        }

    return edge_subtree_count_dict


def filter_graph_edges_by_focus(graph, focused_edge):
    """Filter the graph edges based on the focused edge.

    Args:
        graph: The directed acyclic graph.
        focused_edge (tuple): A tuple containing the names of the nodes representing the focused edge.

    Raises:
        ValueError: If the focused_edge does not exist in the graph.
    """
    # Check if the focused_edge exists in the graph's edges
    focused_edge_found = any(
        (node1.name in focused_edge and node2.name in focused_edge)
        for node1, node2, _ in graph.get_edges()
    )

    if not focused_edge_found:
        raise ValueError(f"Focused edge {focused_edge} not found in the graph!")

    # Filter the graph's edges to only include the focused_edge
    filtered_edges = [
        edge
        for edge in graph.get_edges()
        if edge[0].name in focused_edge and edge[1].name in focused_edge
    ]

    graph.set_edges(filtered_edges)


def get_alignment_look_up_table(alignment):
    """
    Create a lookup table for quick retrieval of alignment record indices.

    This function generates a lookup table that maps sequence identifiers (record IDs)
    to their corresponding indices in the alignment. The lookup table facilitates quick
    access to alignment records during subsequent computations.

    Args:
        alignment (list): A list of sequence alignment records (e.g., Bio.SeqRecord objects).

    Returns:
        dict: A dictionary mapping record IDs to their corresponding indices.

    Note:
        - The alignment records should have unique identifiers.
        - The returned lookup table aids efficient indexing of records during calculations.
    """
    alignment_look_up_table = {}
    i = 0
    for record in alignment:
        alignment_look_up_table[record.id] = i
        i = i + 1
    return alignment_look_up_table
