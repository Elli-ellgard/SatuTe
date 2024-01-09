from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from pandas import DataFrame
import os

""" ## RATE CATEGORIES  """
def get_column_names_with_prefix(data_frame, prefix):
    # Filter the columns using the specified prefix
    columns_with_prefix = data_frame.columns[
        data_frame.columns.str.startswith(prefix)
    ].tolist()
    return columns_with_prefix


def build_categories_by_sub_tables(data_frame: DataFrame):
    rate_category_dictionary = {}

    # Assuming you already have a dataframe called 'dataframe'
    # Call the get_columns_with_prefix function to retrieve columns with a specific prefix
    prefix = "p"  # Specify the desired prefix

    columns_with_prefix = get_column_names_with_prefix(data_frame, prefix)

    # Create dictionaries using column names with the specified prefix as names
    rate_category_dictionary = {column: [] for column in columns_with_prefix}

    for index, row in data_frame.iterrows():
        p_row = row.filter(like="p")

        rate_category_dictionary[p_row.idxmax()].append(int(row["Site"]) - 1)

    return rate_category_dictionary


def parse_category_rates(log_file, number_rates):
    f = open(log_file, "r")
    lines = f.readlines()
    f.close()
    table_data = {}

    # Find the start and end indices of the table
    start_index = -1
    end_index = -1

    for i, line in enumerate(lines):
        if line.strip().startswith("Category"):
            start_index = i + 1
        elif line.strip().startswith(f"{number_rates}"):
            end_index = i + 1
            break

    if start_index == -1 or end_index == -1:
        raise ValueError("Table not found in the log file.")

    # Parse the table rows
    table_lines = lines[start_index:end_index]

    for line in table_lines:
        line = line.strip()
        if line:
            row = line.split()
            category = row[0]
            relative_rate = float(row[1])
            proportion = float(row[2])
            if category != "0":
                table_data[f"p{category}"] = {
                    "Category": category,
                    "Relative_rate": relative_rate,
                    "Proportion": proportion,
                }
    return table_data


""" ## HANDLE ALIGNMENTS  """
def filter_alignment_by_ids(alignment, ids):
    """
    Filter a MultipleSeqAlignment object to include only sequences with specific IDs.

    Parameters:
    alignment (MultipleSeqAlignment): The alignment to filter.
    ids (list of str): The IDs of the sequences to include.

    Returns:
    MultipleSeqAlignment: The filtered alignment.
    """

    # Filter the alignment
    filtered_alignment = MultipleSeqAlignment(
        [record for record in alignment if record.id in ids]
    )

    return filtered_alignment


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
        return None


def change_states_to_allowed(alignment):
    for record in alignment:
        record.seq = record.seq.upper()
        record.seq = record.seq.replace(".", "-")
        record.seq = record.seq.replace("!", "-")

    return alignment


def read_alignment_file(file_name) -> MultipleSeqAlignment:
    """
    Reads an alignment file and returns the alignment object.

    Args:
    - file_name (str): Path to the alignment file.

    Returns:
    - alignment: The alignment object.

    Raises:
    - FileNotFoundError: If the file does not exist or is not readable.
    - ValueError: If the file format could not be guessed or other issues with file content.
    """

    # Check if file exists and is readable
    if not os.path.exists(file_name) or not os.access(file_name, os.R_OK):
        raise FileNotFoundError(
            f"The file {file_name} does not exist or is not readable."
        )

    # Guess the format of the file
    file_format = guess_alignment_format(file_name)

    # If the format could not be guessed, raise an error
    if file_format is None:
        raise ValueError("Could not guess the format of MSA the file.")

    try:
        # Try to read the file in the guessed format
        alignment = AlignIO.read(file_name, file_format)
    except Exception as e:
        # Catch specific exceptions for better error handling
        raise ValueError(f"An error occurred while reading the file: {str(e)}")

    try:
        # Process the alignment
        alignment = change_states_to_allowed(alignment)
    except Exception as e:
        # Catch specific exceptions for better error handling
        raise ValueError(f"An error occurred while processing the alignment: {str(e)}")

    return alignment


def cut_alignment_columns_optimized(alignment, columns):
    """
    Extracts specified columns from a given multiple sequence alignment.

    Parameters:
    - alignment (MultipleSeqAlignment): The input alignment from which columns are to be extracted.
    - columns (list): A list of indices specifying the columns to be extracted.

    Returns:
    - MultipleSeqAlignment: A new alignment containing only the specified columns.
    """

    # Create a new MultipleSeqAlignment from the list of SeqRecord objects, using list comprehension
    selected_records = [
        SeqRecord(Seq("".join(rec.seq[column] for column in columns)), id=rec.id)
        for rec in alignment
    ]

    return MultipleSeqAlignment(selected_records)


def split_msa_into_rate_categories_in_place(
    site_probability, alignment, rate_category
) -> dict[str, MultipleSeqAlignment]:
    """
    Splits a multiple sequence alignment into sub-alignments based on rate categories.

    Parameters:
    - site_probability (dict): A dictionary mapping rate categories to lists of column indices.
    - alignment (MultipleSeqAlignment): The input alignment to be split.
    - rate_category (str): The specific rate category to extract, or "all" to extract all categories.

    Returns:
    - dict: A dictionary mapping rate categories to sub-alignments.
    """

    # Build a dictionary mapping rate categories to lists of column indices
    sub_category = build_categories_by_sub_tables(site_probability)

    # Initialize an empty dictionary to hold the sub-alignments
    per_category_alignment_dict = {}

    # Check if all rate categories should be extracted
    if rate_category == "all":
        # Iterate through each rate category and extract the corresponding columns
        for key, value in sub_category.items():
            per_category_alignment_dict[key] = cut_alignment_columns_optimized(
                alignment, value
            )
    else:
        # Extract only the specified rate category
        key = f"p{rate_category}"
        per_category_alignment_dict[key] = cut_alignment_columns_optimized(
            alignment, sub_category[key]
        )

    return per_category_alignment_dict
