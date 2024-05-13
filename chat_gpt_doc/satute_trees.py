# -*- coding: utf-8 -*-
from ete3 import Tree
from collections import Counter
from typing import Dict, List


def rescale_branch_lengths(tree: Tree, rescale_factor: float):
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
    if not check_all_internal_nodes_annotated(tree):
        idx = 1  # Start indexing for new names
        for node in tree.traverse("preorder"):
            if node.is_leaf() or node.name.startswith("Node"):
                continue  # Skip leaf nodes and nodes already properly named

            # Rename internal nodes with a placeholder name or numeric name
            if not node.name or node.name.isdigit():
                node.name = f"Node{idx}*"
                idx += 1

    # Optional: Check for any duplicate node names which could cause issues
    check_and_raise_for_duplicate_nodes(tree)

    return tree


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


def is_bifurcating(tree) -> bool:
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


def has_duplicate_leaf_sequences(node, multiple_sequence_alignment):
    sequences = [multiple_sequence_alignment[leaf.name] for leaf in node.iter_leaves()]
    sequence_count = Counter(sequences)
    return any(count > 1 for count in sequence_count.values())


def collapse_identical_leaf_sequences(
    tree: Tree, multiple_sequence_alignment: Dict[str, str]
) -> tuple[Dict, Dict]:
    """
    Collapses a tree by merging nodes with identical leaf sequences.

    Args:
    tree (Tree): The tree to be collapsed.
    multiple_sequence_alignment: The multiple sequence alignment associated with the tree's leaves.

    Returns:
    tuple: A tuple containing the updated multiple sequence alignment and a dictionary of collapsed nodes.
    """
    # Initialize a dictionary to track the collapsed nodes
    collapsed_nodes_dict = {}

    # Traverse the tree in postorder
    for node in tree.traverse("postorder"):
        # Check only internal (non-leaf) nodes
        if not node.is_leaf():
            # Get sequences of all leaf nodes under the current internal node
            leaf_sequences = [
                multiple_sequence_alignment[leaf.name] for leaf in node.iter_leaves()
            ]

            # Check if all leaf sequences under this node are identical
            if len(set(leaf_sequences)) == 1:
                # Delete the children nodes as they are identical and update the collapsed nodes dictionary
                delete_children_nodes(node, collapsed_nodes_dict)
                # Update the multiple sequence alignment for the current node
                multiple_sequence_alignment[node.name] = leaf_sequences[0]

    # Return the updated multiple sequence alignment and collapsed nodes dictionary
    return multiple_sequence_alignment, collapsed_nodes_dict


def has_only_leaf_children(node) -> bool:
    # Check if all children of the node are leaves
    return any(child.is_leaf() for child in node.children)


def delete_children_nodes(
    node: Tree, collapsed_nodes_info: Dict[str, List[str]]
) -> Dict:
    """
    Deletes the children nodes of a specified non-leaf node in a phylogenetic tree and updates a dictionary
    with the names of all leaves that were under these nodes.

    This function is primarily used in processes where it's necessary to simplify a tree by collapsing nodes
    with identical sequences, and there's a need to keep track of the original leaf distribution for further analysis.

    Args:
        node (Tree): A non-leaf node from which children will be detached. This node must be part of a larger Tree structure.
        collapsed_nodes_info (dict): A dictionary to be updated with information about the collapsed nodes. It maps node names
        to lists containing the names of all leaf nodes that were under each collapsed node.

    Returns:
        dict: The updated dictionary reflecting the changes made to the tree structure, including the names of leaves under the detached nodes.

    Raises:
        TypeError: If `node` is not an instance of `Tree` or `collapsed_nodes_info` is not a dictionary.
        ValueError: If `node` is a leaf node, indicating it has no children to delete.

    Example:
        >>> from ete3 import Tree
        >>> tree = Tree("((A:1,B:1)Node1:0.5,(C:1,D:1)Node2:0.5)Root;")
        >>> collapsed_nodes_info = {}
        >>> updated_info = delete_children_nodes(tree&"Node1", collapsed_nodes_info)
        >>> print(tree)
        ((C:1,D:1)Node2:0.5)Root;
        >>> print(updated_info)
        {'Node1': ['A', 'B']}

    Note:
        - The function modifies the input tree by detaching child nodes from the specified node. These child nodes are removed from the tree but not deleted from memory.
        - The function assumes the tree and the dictionary are correctly formatted and the node exists within the tree.
    """
    if not isinstance(node, Tree):
        raise TypeError("The node must be an instance of ete3.Tree.")
    if not isinstance(collapsed_nodes_info, dict):
        raise TypeError("collapsed_nodes_info must be a dictionary.")
    if node.is_leaf():
        raise ValueError(
            "The specified node is a leaf and does not have children to delete."
        )

    collapsed_nodes_info[node.name] = [leaf.name for leaf in node.iter_leaves()]
    for child in list(node.children):
        child.detach()

    return collapsed_nodes_info


def print_tree_with_inner_node_names(tree: Tree):
    print(tree.get_ascii(show_internal=True))


if __name__ == "__main__":
    tree = Tree("(((A:0.1,B:0.1):0.1,AB:1),(C:0.2,D:0.1,E:0.1):1);")

    multiple_sequence_alignment = {
        "A": "CATGC",
        "B": "CATGC",
        "AB": "CATGC",
        "AE": "CATGC",
        "C": "ATGT",
        "D": "ATGC",
        "E": "ATGC",
    }

    rename_internal_nodes_pre_order(tree)

    print_tree_with_inner_node_names(tree)

    collapse_identical_leaf_sequences(tree, multiple_sequence_alignment)

    print_tree_with_inner_node_names(tree)
