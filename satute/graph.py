# -*- coding: utf-8 -*-
from typing import Optional
import numpy as np
from ete3 import Tree
from Bio.Align import MultipleSeqAlignment
from typing import List, Tuple, Optional, Dict
from functools import cache
from dataclasses import dataclass, field
from satute.dna_model import NUCLEOTIDE_CODE_VECTOR
from satute.amino_acid_models import AMINO_ACID_CODE_VECTOR


@dataclass
class Node:
    """Class representing a node in the DAG."""

    name: str
    state: Optional[np.array] = field(default_factory=lambda: None)
    connected: Dict["Node", float] = field(default_factory=dict)

    def __hash__(self):
        """Return the hash value of the Node."""
        return id(self)

    def __eq__(self, other):
        """Check if two nodes are equal."""
        return self is other

    def connect(self, other_node: "Node", branch_length: float):
        """Connect the current node to another node with a given branch length."""
        self.connected[other_node] = branch_length

    def is_leaf(self) -> bool:
        """Check if the node is a leaf (has a state vector)."""
        return len(self.connected) <= 1

    def __repr__(self):
        """Return the string representation of the Node."""
        return f"Node({self.name})"


class Graph:
    """Class representing a directed acyclic graph (DAG) of interconnected nodes."""

    def __init__(self, edges: List[Tuple[Node, Node, float]]):
        """Initialize the DAG with a list of directed edges."""
        for left, right, branch_length in edges:
            left.connect(right, branch_length)
            right.connect(left, branch_length)
        self.edges = edges

    def get_edges(self) -> List[Tuple[Node, Node, float]]:
        """Return the list of edges in the DAG."""
        return self.edges

    def set_edges(self, edges: List[Tuple[Node, Node, float]]):
        self.edges = edges

    def get_branch_length(self, edge: Tuple[Node, Node, float]) -> float:
        """Return the branch length of a given edge in the DAG."""
        return self.branch_lengths[edge]


def get_initial_likelihood_vector(state: str, state_type: str) -> np.array:
    if state_type == "nucleotide":
        """Get the initial likelihood vector for a given nucleotide state."""
        return np.array(NUCLEOTIDE_CODE_VECTOR[state]).T
    elif state_type == "amino_acid":
        """Get the initial likelihood vector for a given amino acid state."""
        return np.array(AMINO_ACID_CODE_VECTOR[state]).T


def convert_tree_to_graph(tree: Tree) -> Graph:
    """
    Convert an ETE3 tree to a Directed Acyclic Graph (DAG) representation.

    This function transforms an ETE3 tree into a DAG by creating Node objects for each node in the tree
    and establishing connections between them based on the tree's hierarchical relationships. Each edge in the
    DAG corresponds to a parent-child relationship in the tree, along with the branch length.

    Args:
        tree (Tree): An ETE3 tree object representing the phylogenetic tree.

    Returns:
        Graph: A DAG representation of the phylogenetic tree, with nodes corresponding to tree nodes and edges representing parent-child relationships.
    """

    node_dictionary = {}
    edge_list = []

    # Traverse the tree in level-order
    for ete_node in tree.traverse("levelorder"):
        # Retrieve or create the Node object for the current ETE node
        parent_node_obj = node_dictionary.setdefault(ete_node.name, Node(ete_node.name))

        # Create edges from the parent node to each of its children
        for child_node in ete_node.children:
            # Retrieve or create the Node object for the child
            child_node_obj = node_dictionary.setdefault(
                child_node.name, Node(child_node.name)
            )

            # Create an edge and add it to the edge list
            edge_list.append((parent_node_obj, child_node_obj, child_node.dist))

    # Create and return the Graph object using the constructed edges
    return Graph(edge_list)


def convert_tree_to_state_graph(
    tree: Tree,
    msa_column: MultipleSeqAlignment,
    alignment_look_up_table: Dict[str, str],
    state_type: str,
) -> Graph:
    """
    Convert an ETE3 tree to a Directed Acyclic Graph (DAG) with state information.

    Args:
        tree (Tree): An ETE3 tree object representing the phylogenetic tree.
        msa_column (MultipleSeqAlignment): A single column of the multiple sequence alignment.
        alignment_look_up_table (dict): Maps sequence IDs to alignment indices.

    Returns:
        Graph: A DAG representation of the phylogenetic tree. Nodes represent tree nodes,
        and edges represent parent-child relationships. Leaf nodes contain likelihood vectors.
    """

    node_dictionary = {}
    edge_list: List[Tuple[Node]] = []

    # Traverse the tree and create Node objects and edges
    for ete_node in tree.traverse("levelorder"):

        node = create_or_get_node(
            ete_node=ete_node,
            msa_column=msa_column,
            alignment_look_up_table=alignment_look_up_table,
            node_dictionary=node_dictionary,
            state_type=state_type,
        )

        add_child_edges(
            ete_node,
            node,
            msa_column,
            alignment_look_up_table,
            node_dictionary,
            edge_list,
            state_type,
        )

    return Graph(edge_list)


def create_or_get_node(
    ete_node: Tree,
    msa_column: MultipleSeqAlignment,
    alignment_look_up_table: Dict[str, str],
    node_dictionary: Dict[str, Node],
    state_type: str,
) -> Node:
    """
    Retrieve an existing Node object or create a new one for an ETE node.

    Args:
        ete_node: An individual node from an ETE tree.
        msa_column: The multiple sequence alignment column.
        alignment_look_up_table: Maps sequence IDs to alignment indices.
        node_dictionary: Dictionary holding created Node objects.

    Returns:
        Node: The retrieved or newly created Node object.
    """
    if ete_node.is_leaf():
        likelihood_vector = get_initial_likelihood_vector(
            state=msa_column[alignment_look_up_table[ete_node.name]].seq,
            state_type=state_type,
        )

        return node_dictionary.setdefault(
            ete_node.name, Node(ete_node.name, likelihood_vector)
        )
    else:
        return node_dictionary.setdefault(ete_node.name, Node(ete_node.name))


def add_child_edges(
    parent_ete_node: Tree,
    parent_node: Node,
    msa_column: MultipleSeqAlignment,
    alignment_look_up_table: dict,
    node_dictionary: Dict[str, Node],
    edges: List[Tuple[Node, Node, float]],
    state_type: str,
):
    """
    Create and add edges from a parent node to its children.

    Args:
        parent_ete_node: The current ETE node being processed.
        parent_node: The corresponding Node object for the parent ETE node.
        msa_column: The multiple sequence alignment column.
        alignment_look_up_table: Maps sequence IDs to alignment indices.
        node_dictionary: Dictionary holding created Node objects.
        edge_list: List to which new edges are added.
    """
    for child_ete_node in parent_ete_node.children:
        child_node = create_or_get_node(
            child_ete_node,
            msa_column,
            alignment_look_up_table,
            node_dictionary,
            state_type=state_type,
        )
        edges.append((parent_node, child_node, child_ete_node.dist))


def edge_type(left: Node, right: Node) -> str:
    """
    Determine whether an edge in a tree is internal or external.

    Args:
        left (Node): The left node of the edge.
        right (Node): The right node of the edge.

    Returns:
        str: "external" if either node is a leaf, otherwise "internal".
    """
    return "external" if left.is_leaf() or right.is_leaf() else "internal"


def calculate_subtree_edge_metrics(
    tree: Tree, focused_edges: str
) -> Dict[str, Dict[str, int]]:
    """
    Calculate the number of leaves and branches for each subtree edge within a given tree.

    This function transforms the tree into a directed acyclic graph (DAG) and calculates
    the number of leaves and branches for each subtree defined by the edges of the tree.
    It optionally focuses on specific edges if provided.

    Args:
        tree (Tree): The tree to analyze.
        focused_edges (Optional[List[Tuple[Node, Node]]]): Specific edges to focus on.
            If None, all edges in the tree are considered.

    Returns:
        Dict[str, Dict[str, int]]: A dictionary where each key is an edge represented as
            a string "(left.name, right.name)". Each value is another dictionary containing
            the counts of leaves and branches on both sides of the edge, along with the
            branch length and type.
    """

    count_graph = convert_tree_to_graph(tree)
    edge_metrics = {}

    # Get the number of taxa (leaf nodes)
    num_taxa = len(tree.get_leaves())

    if focused_edges:
        filter_graph_edges_by_focus(count_graph, focused_edges)

    for edge in count_graph.get_edges():
        right, left, branch_length = edge
        edge_name = f"({left.name}, {right.name})"
        _type = edge_type(right, left)

        left_leave_count, left_branch_count = count_leaves_and_branches_for_subtree(
            left, right
        )
        right_leave_count = num_taxa - left_leave_count
        right_branch_count = (2 * num_taxa - 3) - left_branch_count

        edge_metrics[edge_name] = {
            "left": {
                "leave_count": left_leave_count,
                "branch_count": left_branch_count,
            },
            "right": {
                "leave_count": right_leave_count,
                "branch_count": right_branch_count,
            },
            "length": branch_length,
            "type": _type,
        }

    return edge_metrics


@cache
def count_leaves_and_branches_for_subtree(node: Node, coming_from: Node):
    """
    Recursively counts the number of leaf nodes and branches in a subtree defined by a given node in a tree.

    This function is designed to work with a tree structure where each node maintains a list of connected nodes.
    It starts from a specified node ('node') and traverses the subtree, counting the number of leaf nodes and branches.
    The traversal avoids backtracking to the 'coming_from' node to prevent double counting.

    Args:
        node (Node): The node from which the subtree counting begins. This node acts as the root of the subtree.
        coming_from (Node): The node from which the traversal to the current 'node' occurred. This is used to
                            prevent the traversal from moving back up the tree, thereby avoiding the branch
                            that was used to reach the 'node'.

    Returns:
        tuple: A tuple (leaf_count, branch_count) representing the count of leaf nodes and branches within the subtree.
        leaf_count (int): The number of leaf nodes in the subtree rooted at 'node'.
        branch_count (int): The number of branches in the subtree rooted at 'node', excluding the branch
        from 'coming_from' to 'node'.

    Note:
        - The function assumes that each node has an attribute 'connected', which is a dictionary with keys
        representing connected child nodes.
        - It is primarily used for traversing phylogenetic trees or similar hierarchical structures where each node can have multiple children.
        - The function is recursive and may not be suitable for extremely large trees due to Python's recursion limit.
    """
    # If the current node is a leaf
    if node.is_leaf():
        return 1, 0

    # Initialize leaf and branch counts
    leaf_count, branch_count = 0, 0

    # Iterate over children
    for child in node.connected.keys():
        if child != coming_from:
            # Recursively count leaves and branches for each child
            child_leaves, child_branches = count_leaves_and_branches_for_subtree(
                child, node
            )
            leaf_count += child_leaves
            branch_count += child_branches + 1  # Include the branch to this child
    return leaf_count, branch_count


def filter_graph_edges_by_focus(graph: Graph, focused_edge: str):
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

    return filtered_edges


def get_alignment_look_up_table(alignment: MultipleSeqAlignment) -> dict:
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
