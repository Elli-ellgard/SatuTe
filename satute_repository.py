import os
import numpy as np
import pandas as pd
import re

# Define a function to parse rate parameters from a file
def parse_rate_parameters(file_path, dimension, model="GTR"):
    # Open the file in read mode
    with open(file_path, "r") as f:
        # Initialize variables
        found = 0
        number_lines = 0
        rates = []

        # Loop through each line in the file
        for line in f:
            # If we find the line that starts the rate parameters, set 'found' to true
            if "Rate parameter R:" in line:
                found = 1
            if found and number_lines < dimension * (dimension - 1) // 2 + 1:
                if number_lines > 1:
                    # Split the line at the ':' character and append the second part (the rate)
                    # to our list of rates, after stripping whitespace from its ends
                    line = line.split(":")
                    rates.append(line[1].strip())
                # Increment the number of lines read
                number_lines += 1

        # Split the model string at the '+' character
        splitted_model = model.split("+")
        # Remove duplicate rates and join them into a comma-separated string
        rates_list = ",".join(list(dict.fromkeys(rates)))

        # Format a new model string with the rates inside curly brackets
        model_with_rates_token = f"{splitted_model[0]}{{{rates_list}}}"

        # Return the new model string
        return model_with_rates_token


def parse_output_state_frequencies(file_path):
    df = pd.read_csv(file_path, comment="#", sep="\t", engine="c")
    return df


# Define a function to parse state frequencies from a log content
def parse_state_frequencies(log_file_path, dimension=4):
    # Read the content of the log file
    with open(log_file_path, "r") as f:
        log_content = f.read()
    # If the string "equal frequencies" is not in the log content, proceed to parse
    # Initialize an empty dictionary to hold the frequencies
    frequencies = {}
    if "equal frequencies" not in log_content:
        # Define the start and end markers for the section of the file we want to parse
        start_line = "State frequencies: (empirical counts from alignment)"
        end_line = "Rate matrix Q:"
        # Split the log content into lines
        lines = log_content.split("\n")
        # Find the indices of the start and end markers
        start_index = lines.index(start_line)
        end_index = lines.index(end_line)
        # Get the lines in between the start and end markers
        freq_lines = lines[start_index + 1 : end_index]

        # Loop through the lines we've extracted
        for line in freq_lines:
            # If the line is not empty after removing leading/trailing whitespace
            if line.strip():
                # Split the line on " = " into a key and a value, and add them to the frequencies dictionary
                key, value = line.strip().split(" = ")
                frequencies[key] = float(value)  # convert value to float before storing
    else:
        # If "equal frequencies" is in the log content, return a pseudo dictionary with equal frequencies
        for i in range(dimension):
            key = "key_" + str(i)
            frequencies[key] = 1 / dimension

    # Return the frequencies dictionary
    return frequencies


def extract_rate_matrix(file_path):
    file_content = ""
    # Load the file content into a string
    with open(file_path, "r") as file:
        file_content = file.read()

    rate_matrix_start = file_content.find("Rate matrix Q:")
    rate_matrix_text = file_content[rate_matrix_start:]
    rate_matrix_lines = rate_matrix_text.split("\n")

    rate_matrix = []
    row_ids = []

    for line in rate_matrix_lines[1:]:  # Skip the "Rate matrix Q:" line
        # Check if the line starts with a letter (assumes row IDs are letters)
        if not line.strip() or not line.strip()[0].isalpha():
            continue
        tokens = line.split()
        try:
            # Try to convert the tokens after the first one to floats
            row = [float(x) for x in tokens[1:]]
        except ValueError:
            # If the conversion fails, we've reached the end of the matrix
            break
        row_ids.append(tokens[0])
        rate_matrix.append(row)

    return pd.DataFrame(rate_matrix, index=row_ids, columns=row_ids)


def parse_rate_matrices_from_file(file_path):
    """
    Parse the rate matrix Q and stationary distribution pi from a given .iqtree file path.

    Parameters:
    - file_path (str): Path to the .iqtree file.

    Returns:
    - rate_matrix (np.array): The parsed rate matrix Q.
    - phi_matrix (np.array): The stationary distribution pi with values filled in the diagonal.
    """

    # Check if file exists
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"The file '{file_path}' does not exist.")

    with open(file_path, "r") as file:
        lines = file.readlines()

    # Dynamically determine the matrix dimension 'n'
    start_idx = next(
        (idx for idx, line in enumerate(lines) if "Rate matrix Q:" in line), None
    )
    if start_idx is None:
        raise ValueError("'Rate matrix Q:' not found in file.")

    # Detect the number of matrix rows based on numeric entries
    n = 0
    current_idx = start_idx + 2  # Adjusting to start from matrix values
    while current_idx < len(lines) and re.search(
        r"(\s*-?\d+\.\d+\s*)+", lines[current_idx]
    ):
        n += 1
        current_idx += 1

    rate_matrix = np.zeros((n, n))
    phi_matrix = np.zeros((n, n))

    # Parse the rate matrix Q
    for idx, line in enumerate(lines):
        if "Rate matrix Q:" in line:
            try:
                for j in range(n):
                    entries = lines[idx + j + 2].split()
                    for k in range(n):
                        rate_matrix[j, k] = (
                            0 if "e" in entries[k + 1] else float(entries[k + 1])
                        )
            except (IndexError, ValueError) as e:
                raise Exception(f"Error while parsing rate matrix. Exception: {e}")

        # Parse the state frequencies (empirical counts)
        elif "State frequencies: (empirical counts from alignment)" in line:
            try:
                for j in range(n):
                    phi_matrix[j, j] = float(lines[idx + j + 2].split()[2])
            except (IndexError, ValueError) as e:
                raise Exception(
                    f"Error while parsing empirical state frequencies. Exception: {e}"
                )

        # If the state frequencies are equal, fill the diagonal accordingly
        elif "State frequencies: (equal frequencies)" in line:
            np.fill_diagonal(phi_matrix, 1 / n)

    return rate_matrix, phi_matrix


# Return the function for review without executing it


# Sample test call (without an actual file for now)
# parse_rate_matrices("sample_path")


def parse_rate_and_frequencies_and_model(
    input_path, dimension, model="GTR"
):
    # Construct the model string with parsed rate parameters
    log_file_path = f"{input_path}.iqtree"
    model_final = parse_rate_parameters(log_file_path, dimension, model=model)

    # Parse state frequencies from the log content
    state_frequencies = parse_state_frequencies(log_file_path, dimension=dimension)

    # Create a string of state frequencies separated by a space
    concatenated_rates = " ".join(map(str, state_frequencies.values()))

    # Construct command line tokens for model and frequency
    model_and_frequency = f"{model_final}+FU{{{concatenated_rates}}}"

    return state_frequencies.values(), model_and_frequency
