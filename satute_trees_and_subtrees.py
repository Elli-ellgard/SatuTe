from ete3 import Tree
import pandas as pd
import re

def parse_file_to_data_frame(file_path):
    try:
        # Read the file into a dataframe
        df = pd.read_csv(file_path, delimiter="\t")

        return df

    except FileNotFoundError:
        raise Exception(f"File not found: {file_path}")


def parse_rate_matrix(file_content):
    rate_matrix_start = file_content.find("Rate matrix Q:")
    rate_matrix_text = file_content[rate_matrix_start:]
    rate_matrix_lines = rate_matrix_text.split("\n")

    rate_matrix = []
    row_ids = []

    for line in rate_matrix_lines:
        # Continue until line doesn't start with a letter (assumes row IDs are letters)
        if line.strip() == "" or not line.strip()[0].isalpha():
            break
        tokens = line.split()
        row_ids.append(tokens[0])
        rate_matrix.append([float(x) for x in tokens[1:]])

    return pd.DataFrame(rate_matrix, index=row_ids, columns=row_ids)


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


def get_leaves(tree):
    leaves = []
    for leaf in tree:
        if leaf.is_leaf():
            leaves.append(leaf.name)
    return leaves


def branch_lengths(T):
    vector_branches = []
    vector_distances = []
    internal_nodes = []

    for node in T.traverse("levelorder"):
        if (node not in T.get_leaves()) and (node in T.get_tree_root()):
            internal_nodes.append(node.name)

    for node in T.traverse("levelorder"):  # First internal branches.
        children = node.get_children()

        for child in children:
            if child.name in internal_nodes:
                vector_branches.append(node.name + "-" + child.name)
                vector_distances.append(T.get_distance(node, child))

    for node in T.traverse("levelorder"):  # Same for external branches.
        children = node.get_children()

        for child in children:
            if child.name not in internal_nodes:
                vector_branches.append(node.name + "-" + child.name)
                vector_distances.append(T.get_distance(node, child))

    return vector_branches, vector_distances


def name_nodes_by_level_order(tree):
    i = 1
    for node in tree.traverse("levelorder"):
        if not node.is_leaf():
            node.name = f"Node{i}*"
            i += 1
    return tree


def rescale_branch_lengths(tree, rescale_factor):
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


def get_target_node(edge):
    return edge.split(",")[0].replace("(", "").strip()


def map_values_to_newick(newick, df):
    for index, row in df.iterrows():
        target_node = get_target_node(row["edge"])

        # Escape any special characters in target node for regex
        escaped_target_node = re.escape(target_node)

        # Adjusting the metadata as per the requirements
        meta_data = f"&delta={row['delta']}, &c_s={row['c_s']}, &p_value={row['p_value']}, &result_test={row['result_test']}"

        # Check for existing square brackets after the node
        pattern_with_brackets = re.compile(
            f"({escaped_target_node}:\d+(\.\d+)?(e-?\d+)?)\[([^\]]+)\]"
        )
        pattern_without_brackets = re.compile(
            f"({escaped_target_node}:\d+(\.\d+)?(e-?\d+)?)"
        )

        # If square brackets are present, append the metadata inside those brackets
        if pattern_with_brackets.search(newick):
            newick = pattern_with_brackets.sub(f"\\1[\\4,{meta_data}]", newick)
        else:
            # If no square brackets, add them and insert the metadata
            newick = pattern_without_brackets.sub(f"\\1[{meta_data}]", newick)

    return newick


