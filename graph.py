from typing import Optional
import dataclasses
import numpy as np


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