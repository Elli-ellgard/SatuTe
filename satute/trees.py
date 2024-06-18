# -*- coding: utf-8 -*-
from ete3 import Tree
from collections import Counter
from typing import Tuple, Set


###############              Printing                 ######################
def print_tree_with_inner_node_names(tree: Tree):
    print(tree.get_ascii(show_internal=True))


###############              Tree Annotation          ######################
def rescale_branch_lengths(tree: Tree, rescale_factor: float) -> None:
    """
    Rescales the branch lengths of a phylogenetic tree by a given factor.

    Args:
        tree (ete3.Tree): The input phylogenetic tree object.
        rescale_factor (float): Scaling factor for the branch lengths.

    Returns:
        ete3.Tree: The phylogenetic tree object with rescaled branch lengths.
    """
    # Iterate over tree nodes and rescale branch lengths
    for node in tree.traverse():
        if node.up:
            node.dist *= rescale_factor


def rename_internal_nodes_pre_order(tree: Tree) -> Tree:
    """
    Modifies the input tree in-place by naming its internal nodes using a pre-order traversal.
    Internal nodes are named as "NodeX*" where X is an incremental number. This function skips renaming
    for nodes that already start with "Node" to avoid redundancy and does not rename numeric node names,
    considering them as not annotated.
    Args:
        tree (Tree): The input tree to be modified.

    Returns:
        Tree: The same tree instance with updated internal node names. Note that this function modifies the tree in-place.
    """
    number_annotated_internal_nodes, number_unannotated_internal_nodes, seen_names = (
        get_set_and_number_annotated_internal_nodes(tree)
    )
    if number_unannotated_internal_nodes != 0:
        for node in tree.traverse("preorder"):
            if node.is_leaf() or node.name.startswith("Node"):
                continue  # Skip leaf nodes and nodes already properly named

            # Rename internal nodes with a placeholder name or numeric name
            if not node.name or node.name.isdigit():
                new_name = f"Node{number_annotated_internal_nodes + 1}*"
                while new_name in seen_names:
                    number_annotated_internal_nodes += 1
                    new_name = f"Node{number_annotated_internal_nodes + 1}*"

                node.name = new_name
                seen_names.add(new_name)
                number_annotated_internal_nodes += 1

    # Optional: Check for any duplicate node names which could cause issues
    check_and_raise_for_duplicate_nodes(tree)

    return tree


def get_set_and_number_annotated_internal_nodes(
    tree: Tree,
) -> Tuple[int, int, Set[str]]:
    """
    Analyzes the internal nodes of a given tree and categorizes them into annotated and unannotated.

    Annotated nodes are those that have names starting with "Node" and are neither empty nor purely numeric.
    Unannotated nodes are either unnamed, numeric named, or have a generic placeholder name.

    The function performs a pre-order traversal of the tree to ensure each node is visited in pre-order sequence.

    Args:
        tree (Tree): The tree to analyze for node annotation status.

    Returns:
        Tuple[int, int, Set[str]]:
            - First int represents the count of annotated internal nodes.
            - Second int represents the count of unannotated internal nodes.
            - Set[str] contains the unique names of the annotated internal nodes.

    This function does not modify the input tree.
    """
    seen_names = set()
    number_annotated = 0
    number_unannotated = 0
    for node in tree.traverse("preorder"):
        if node.is_leaf():
            continue  # Skip leaf nodes to focus on internal nodes

        if not node.name or node.name.isdigit() or not node.name.startswith("Node"):
            number_unannotated += 1
        else:
            seen_names.add(node.name)
            number_annotated += 1

    return number_annotated, number_unannotated, seen_names


def check_and_raise_for_duplicate_nodes(tree: Tree) -> None:
    """
    Checks the given tree for any duplicate node names and raises an exception if duplicates are found.
    This ensures the tree structure is uniquely identifiable and consistent for downstream analysis.

    Args:
        tree (Tree): The tree to check for duplicate node names.

    Raises:
        ValueError: If a duplicate node name is found in the tree.
    """
    seen_names = set()
    for node in tree.traverse("preorder"):
        if node.name in seen_names:
            raise ValueError(f"Duplicate node name found: '{node.name}'")
        seen_names.add(node.name)


def check_all_internal_nodes_annotated(tree: Tree) -> bool:
    """
    Checks if all internal nodes in the given tree are annotated, where an internal node is considered annotated
    if it has a non-empty, non-numeric name and does not start with a generic prefix like "Node".

    Args:
        tree (Tree): The tree to check for annotated internal nodes.

    Returns:
        bool: True if all internal nodes are annotated, False otherwise.
    """
    for node in tree.traverse("preorder"):
        if node.is_leaf():
            continue  # Focus on internal nodes
        # An internal node is unannotated if its name is empty, numeric, or a generic placeholder
        if not node.name or node.name.isdigit() or not node.name.startswith("Node"):
            return False
    return True


###############              Tree Check          ######################
def is_bifurcating(tree: Tree) -> bool:
    """
    Checks if a given tree is strictly bifurcating.

    Args:
    tree (ete3.Tree): The tree to be checked.

    Returns:
    bool: True if the tree is strictly bifurcating, False otherwise.
    """
    for node in tree.traverse():
        # Ignore leaf nodes
        if not node.is_leaf():
            # Check if the number of children is exactly two
            print(node.name, [child.name for child in node.children])
            if len(node.children) != 2:
                return False
    return True


def has_duplicate_leaf_sequences(node, multiple_sequence_alignment) -> bool:
    sequences = [multiple_sequence_alignment[leaf.name] for leaf in node.iter_leaves()]
    sequence_count = Counter(sequences)
    return any(count > 1 for count in sequence_count.values())
