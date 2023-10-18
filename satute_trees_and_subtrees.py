from ete3 import Tree
import re


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
