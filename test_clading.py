from ete3 import Tree
import os
import pandas as pd
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import numpy as np
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def name_nodes_by_level_order(tree):
    i = 1
    for node in tree.traverse("levelorder"):
        if not node.is_leaf():
            node.name = f"Node{i}*"
            i += 1
    return tree


def get_all_subtrees(tree):
    subtrees = []
    for node in tree.traverse("levelorder"):
        if not node.is_root():
            subtrees.append(node)
    return subtrees


def get_opposite_subtree(tree, subtree):
    # Create a copy of the tree so as not to modify the original
    tree_copy = Tree(tree.write())
    tree_copy = name_nodes_by_level_order(tree_copy)
    # Get the node in the copied tree that corresponds to the root of the subtree
    node = tree_copy.search_nodes(name=subtree.name)[0]
    # Detach this node (and its descendants) from the tree
    node.detach()
    return tree_copy


def generate_subtree_pairs(subtrees, t):
    subtree_pairs = []
    for subtree in subtrees:
        subtree_copy = t.search_nodes(name=subtree.name)[0]
        opposite_subtree = get_opposite_subtree(t, subtree_copy)

        subtree_pairs.append(
            (subtree.write(format=1), opposite_subtree.write(format=1))
        )

    return subtree_pairs


def write_subtree_pairs(generated_subtree_pairs, path_prefix="./"):
    for i in range(len(generated_subtree_pairs)):
        os.makedirs(f"{path_prefix}clades/Branch{i}_clade1/", exist_ok=True)
        os.makedirs(f"{path_prefix}clades/Branch{i}_clade2/", exist_ok=True)

        first_subtree = generated_subtree_pairs[i][0]
        second_subtree = generated_subtree_pairs[i][1]

        first_clade_writer = open(
            f"{path_prefix}clades/Branch{i}_clade1/subtree.treefile", "w"
        )
        second_clade_writer = open(
            f"{path_prefix}clades/Branch{i}_clade2/subtree.treefile", "w"
        )

        first_clade_writer.write(first_subtree)
        second_clade_writer.write(second_subtree)

        first_clade_writer.close()
        second_clade_writer.close()


def parse_file_to_dataframe(file_path):
    try:
        # Read the file into a dataframe
        print(file_path)
        df = pd.read_csv(file_path, delimiter="\t")

        return df

    except FileNotFoundError:
        raise Exception(f"File not found: {file_path}")


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


def guess_alignment_format(file_name):
    with open(file_name, "r") as f:
        first_line = f.readline().strip()

    # Now we'll check for various signatures that might indicate the format.
    if first_line.startswith(">"):
        return "fasta"
    elif first_line.startswith("CLUSTAL"):
        return "clustal"
    elif first_line.startswith("# STOCKHOLM"):
        return "stockholm"
    elif first_line.startswith("#NEXUS"):
        return "nexus"
    elif first_line.startswith("PileUp"):
        return "pileup"
    elif first_line[0].isdigit():
        return "phylip"
    else:
        return "Unknown"


def generate_and_write_subtree_pairs(number_rates, t, file_name):
    # Iterate from 0 to number_rates (inclusive)
    for i in range(number_rates + 1):
        # If this is the first iteration
        if i == 0:
            try:
                # Call function to get all subtrees from t, assuming the function is defined elsewhere
                subtrees = get_all_subtrees(t)
                # Discard the first subtree, assuming we don't need it
                subtrees = subtrees[1:]
                # Generate subtree pairs, assuming the function is defined elsewhere
                generated_subtree_pairs = generate_subtree_pairs(subtrees, t)
                # Write subtree pairs to files, assuming the function is defined elsewhere
                write_subtree_pairs(generated_subtree_pairs, f"./{file_name}/")
            except Exception as e:
                print(f"Error occurred during the first iteration: {e}")
        # For the rest of the iterations
        else:
            try:
                # Make a new directory for this subsequence, if it doesn't exist yet
                os.makedirs(f"./{file_name}/subsequence{i}/", exist_ok=True)
                # Write subtree pairs to files in the new directory
                write_subtree_pairs(
                    generated_subtree_pairs,
                    path_prefix=f"./{file_name}/subsequence{i}/",
                )
            except Exception as e:
                print(f"Error occurred during iteration {i}: {e}")


def get_column_names_with_prefix(dataframe, prefix):
    # Filter the columns using the specified prefix
    columns_with_prefix = dataframe.columns[
        dataframe.columns.str.startswith(prefix)
    ].tolist()
    return columns_with_prefix


def build_categories_by_subtables(dataframe):
    subtables = []

    rate_category_dictionary = {}

    # Assuming you already have a dataframe called 'dataframe'
    # Call the get_columns_with_prefix function to retrieve columns with a specific prefix
    prefix = "p"  # Specify the desired prefix

    columns_with_prefix = get_column_names_with_prefix(dataframe, prefix)

    # Create dictionaries using column names with the specified prefix as names
    rate_category_dictionary = {column: [] for column in columns_with_prefix}

    for index, row in dataframe.iterrows():
        p_row = row.filter(like="p")

        rate_category_dictionary[p_row.idxmax()].append(row["Site"])

    return rate_category_dictionary


def read_alignment_file(file_name):
    # Guess the format of the file
    file_format = guess_alignment_format(file_name)

    # If the format could not be guessed, raise an error
    if file_format is None:
        raise ValueError("Could not guess the format of the file.")

    # Try to read the file in the guessed format
    try:
        alignment = AlignIO.read(file_name, file_format)
        return alignment
    except Exception as e:
        print(f"An error occurred while reading the file: {str(e)}")


def cut_alignment_columns(alignment, columns):
    # Convert the alignment to a NumPy array for easy column slicing
    alignment_array = np.array([list(rec) for rec in alignment], np.character)
    # Select the specified columns
    selected_columns = alignment_array[:, columns]
    selected_records = []
    for i, rec in enumerate(selected_columns):
        selected_records.append(
            SeqRecord(Seq(rec.tobytes().decode()), id=alignment[i].id)
        )
    selected_alignment = MultipleSeqAlignment(selected_records)
    return selected_alignment


def write_alignment(alignment, file_name, file_format):
    """
    Write a MultipleSeqAlignment object to a file.

    Parameters:
    alignment (MultipleSeqAlignment): The alignment to write.
    file_name (str): The name of the file to write to.
    file_format (str): The format of the alignment file. This must be a string specifying one
                       of the formats that Biopython supports, such as "fasta", "clustal", "phylip", etc.

    Returns:
    int: The number of alignments written (should be 1 for a MultipleSeqAlignment object).
    """

    # Write the alignment to the file
    num_alignments_written = AlignIO.write(alignment, file_name, file_format)

    return num_alignments_written


##############################################################################################################
number_rates = 4
file_format = guess_alignment_format("./test/octo-kraken-msa-test/example.phy")
t = name_nodes_by_level_order(
    parse_newick_file("./test/octo-kraken-msa-test/example.phy.treefile")
)
file_name = "test_cladding_and_subsequence"
site_probability = parse_file_to_dataframe(
    "./test/octo-kraken-msa-test/example.phy.siteprob"
)
# Call the function
generate_and_write_subtree_pairs(number_rates, t, file_name)
# Assuming you already have a dataframe called 'dataframe'
# Call the split_dataframe_by_highest_value function to split the dataframe by highest values
sub_category = build_categories_by_subtables(site_probability)
alignment = read_alignment_file("./test/octo-kraken-msa-test/example.phy")
per_category_alignment_dict = {}

for key, value in sub_category.items():
    per_category_alignment_dict[key] = cut_alignment_columns(alignment, value)


for key, value in per_category_alignment_dict.items():
    write_alignment(value, f"./{file_name}/subsequence{key[1:]}/rate.fast", "fasta")
