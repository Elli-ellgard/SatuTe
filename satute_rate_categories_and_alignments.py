from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import numpy as np
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

""" ## RATE CATEGORIES  """


def get_column_names_with_prefix(data_frame, prefix):
    # Filter the columns using the specified prefix
    columns_with_prefix = data_frame.columns[
        data_frame.columns.str.startswith(prefix)
    ].tolist()
    return columns_with_prefix


def build_categories_by_sub_tables(data_frame):
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
        elif line.strip().startswith(
            f"{number_rates}"
        ):
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
        return "Unknown"


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


def split_msa_into_rate_categories_in_place(site_probability, alignment):
    sub_category = build_categories_by_sub_tables(site_probability)
    per_category_alignment_dict = {}
    for key, value in sub_category.items():
        per_category_alignment_dict[key] = cut_alignment_columns(alignment, value)
    return per_category_alignment_dict
