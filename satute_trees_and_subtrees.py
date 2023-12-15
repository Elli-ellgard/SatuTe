from ete3 import Tree
from collections import Counter


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


def rename_internal_nodes_preorder(tree: Tree):
    """
    Modify the input tree by naming its nodes using a preorder traversal.
    Nodes are named as "NodeX*" where X is an incremental number.
    If a node name is purely numeric, it is preserved as 'apriorism' feature of the node.
    Args:
        t (Tree): The input tree to be modified.
    Returns:
        Tree: The modified tree with updated node names.
    """
    idx = 1
    for node in tree.traverse("preorder"):
        if not node.is_leaf():
            # If a node name already starts with "Node", no further modification is required.
            if node.name.startswith("Node"):
                return tree
            # Preserve numeric node names as 'apriori' for reference.
            if node.name.isdigit():
                node.add_features(apriori=node.name)
            # Assign new node names based on the preorder traversal index.
            node.name = f"Node{idx}*"
            idx += 1
    return tree


def has_duplicate_leaf_sequences(node, multiple_sequence_alignment):
    sequences = [multiple_sequence_alignment[leaf.name] for leaf in node.iter_leaves()]
    sequence_count = Counter(sequences)
    return any(count > 1 for count in sequence_count.values())


def collapse_tree(tree: Tree, multiple_sequence_alignment):
    twin_dictionary = {}

    for node in tree.traverse("postorder"):
        if not node.is_leaf():
            leaf_sequences = [
                multiple_sequence_alignment[leaf.name] for leaf in node.iter_leaves()
            ]
            if len(set(leaf_sequences)) == 1:  # All leaf sequences are identical
                delete_children_nodes(node, twin_dictionary)
                multiple_sequence_alignment[node.name] = leaf_sequences[0]
    return multiple_sequence_alignment, twin_dictionary


def has_only_leaf_children(node):
    # Check if all children of the node are leaves
    return any(child.is_leaf() for child in node.children)


def summarize_similar_leaves(node, multiple_sequence_alignment, sibling_dictionary):
    seen_sequences = {}

    for child in list(node.children):
        if child.name in multiple_sequence_alignment:
            if multiple_sequence_alignment[child.name] in seen_sequences:
                sibling_dictionary[child.name] = seen_sequences[
                    multiple_sequence_alignment[child.name]
                ]

                seen_sequences[multiple_sequence_alignment[child.name]].append(
                    str(child.name)
                )

                child.delete()
            else:
                # Keep the first occurrence of each sequence
                seen_sequences[multiple_sequence_alignment[child.name]] = [
                    str(child.name)
                ]

    for child in list(node.children):
        if child.name in multiple_sequence_alignment:
            if multiple_sequence_alignment[child.name] in seen_sequences:
                sibling_dictionary[child.name] = seen_sequences


def delete_children_nodes(node: Tree, twin_dictionary: dict):
    twin_dictionary[node.name] = [leaf.name for leaf in node.iter_leaves()]
    # Gather names of all leaf children for record-keeping
    for child in list(node.children):
        # Delete the child node
        child.detach()


def print_tree_with_inner_node_names(tree):
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
    rename_internal_nodes_preorder(tree)
    print_tree_with_inner_node_names(tree)
    collapse_tree(tree, multiple_sequence_alignment)
    print_tree_with_inner_node_names(tree)
