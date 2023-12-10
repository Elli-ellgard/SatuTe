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

def parse_state_frequencies_from_file(file_path):
    """
    Parse the stationary distribution pi from a given .iqtree file path.

    Parameters:
    - file_path (str): Path to the .iqtree file.

    Returns
    - frequencies (directory): The stationary distribution.
    - phi_matrix (np.array): The stationary distribution pi with values filled in the diagonal.
    """

    # Check if file exists
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"The file '{file_path}' does not exist.")

    with open(file_path, "r") as file:
        lines = file.readlines()

    index = next(
        (idx for idx, line in enumerate(lines) if "State frequencies:" in line), None
    )
    if index is None:
        raise ValueError("'State frequencies:' not found in file.")

    # Dynamically determine the dimension 'n'
    start_idx = next(
        (idx for idx, line in enumerate(lines) if "Rate matrix Q:" in line), None
    )
    if start_idx is None:
        raise ValueError("'Rate matrix Q:' not found in file. Determination of dimension not possible.")

    # Detect the number of matrix rows based on numeric entries
    n = 0
    current_idx = start_idx + 2  # Adjusting to start from matrix values
    while current_idx < len(lines) and re.search(
        r"(\s*-?\d+\.\d+\s*)+", lines[current_idx]
    ):
        n += 1
        current_idx += 1
    
    # Initialize an empty dictionary to hold the frequencies
    frequencies = {}

    # Parse the state frequencies
    for idx, line in enumerate(lines):
        # Parse the state frequencies (empirical counts)
        if "State frequencies: (empirical counts from alignment)" in line:
            try:
                for i in range(n):
                   # Split the line on " = " into a key and a value, and add them to the frequencies dictionary
                   key, value = lines[idx + i + 2].split(" = ")
                   frequencies[key] = float(value)  # convert value to float before storing
            except (IndexError, ValueError) as e:
                raise Exception(
                    f"Error while parsing empirical state frequencies. Exception: {e}"
                )

        # If "equal frequencies" is in the log content, return a pseudo dictionary with equal frequencies
        elif "State frequencies: (equal frequencies)" in line:
            for i in range(n):
                key = "key_" + str(i)
                frequencies[key] = float(1 / n)

    phi_matrix = np.diag(list(frequencies.values()))

    return frequencies, phi_matrix

def parse_rate_matrices_from_file_new(file_path, state_frequencies):
    """
    Parse the rate parameters R  .iqtree file path and determine 
    the rate matrix Q using the rate parameters and stationary distribution

    Parameters:
    - file_path (str): Path to the .iqtree file.

    Returns:
    - rate_matrix (np.array): The parsed rate matrix Q.
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

    # Find the start and end indices of the rate parameter
    start_index = -1
    end_index = -1

    for i, line in enumerate(lines):
        if line.strip().startswith("Rate parameter R:"):
            start_index = i + 1
        elif line.strip().startswith("State frequencies"):
            end_index = i
            break

    if start_index == -1 or end_index == -1:
        raise ValueError("Rate parameter not found in the log file.")
    
    # Get the lines in between the start and end markers
    parameter = lines[start_index + 1 : end_index]
    rates = {}

    # Loop through the lines we've extracted
    for line in parameter:
        # If the line is not empty after removing leading/trailing whitespace
        if line.strip():
           # Split the line on " = " into a key and a value, and add them to the frequencies dictionary
            key, value = line.strip().split(":")
            rates[key] = float(value)  # convert value to float before storing
    
    only_rates=list(rates.values())
    list_state_freq=list(state_frequencies.values())
    idx = 0
    rate_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1,n):
            rate_matrix[i,j] = only_rates[idx] * list_state_freq[j]
            rate_matrix[j,i] = only_rates[idx] * list_state_freq[i]
            idx += 1
    average_rate = 0
    for i in range(n):
        rate_matrix[i,i] = - sum(rate_matrix[i,])
        average_rate = average_rate + rate_matrix[i,i] * list_state_freq[i]    
    rate_matrix = - rate_matrix / average_rate

    return rate_matrix


def parse_rate_matrices_from_file(file_path):
    """
    Parse the rate matrix Q  from a given .iqtree file path.

    Parameters:
    - file_path (str): Path to the .iqtree file.

    Returns:
    - rate_matrix (np.array): The parsed rate matrix Q.
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

    return rate_matrix

# Return the function for review without executing it


# Sample test call (without an actual file for now)


def parse_rate_and_frequencies_and_model(input_path, dimension, model="GTR"):
    # Construct the model string with parsed rate parameters
    log_file_path = f"{input_path}.iqtree"
    model_final = parse_rate_parameters(log_file_path, dimension, model=model)

    # Parse state frequencies from the log content
    state_frequencies = parse_state_frequencies_from_file(log_file_path)

    # Create a string of state frequencies separated by a space
    concatenated_rates = " ".join(map(str, state_frequencies.values()))

    # Construct command line tokens for model and frequency
    model_and_frequency = f"{model_final}+FU{{{concatenated_rates}}}"

    return state_frequencies.values(), model_and_frequency


def parse_substitution_model(file_path):
    try:
        with open(file_path, "r") as file:
            content = file.read()
            for line in content.splitlines():
                if "Best-fit model according to BIC:" in line:
                    model_string = line.split(":")[1].strip()
                    return model_string
                if "Model of substitution:" in line:
                    model_string = line.split(":")[1].strip()
                    return model_string
            raise ValueError("Could not parse the substitution model from the file.")
    except (IOError, ValueError):
        # If the file cannot be read or ':' is not found in the content
        # Return None or an appropriate value for error handling
        raise ValueError("Could not parse the substitution model from the file.")
    

def parse_rate_from_model(model):
    # Find the index of '+G' and '+R' in the model string
    plus_g_index = model.find("+G")
    plus_r_index = model.find("+R")

    if plus_g_index != -1 and plus_r_index != -1:
        raise ValueError("Cannot use +G and +R")

    if plus_g_index != -1:
        rate_start_index = plus_g_index + 2
    elif plus_r_index != -1:
        rate_start_index = plus_r_index + 2
    else:
        return 1  # default number_rates = 1 if no +G or +R model

    try:
        # Extract the substring after '+G'
        number = model[rate_start_index:]

        # Parse the extracted substring as an integer
        rate = int(number)

        return rate
    except ValueError:
        # If '+G' is not found or the number after '+G' is not a valid integer
        # Return None or an appropriate value for error handling
        raise ValueError("Could not parse the substitution model from the file.")
        

def parse_file_to_data_frame(file_path):
    try:
        # Read the file into a dataframe
        df = pd.read_csv(file_path, delimiter="\t")

        return df

    except FileNotFoundError:
<<<<<<< Updated upstream
        raise Exception(f"File not found: {file_path}")
=======
        raise Exception(f"File not found: {file_path}")


class SubstitutionModel:
    """
    A class representing a substitution model used in phylogenetics or related fields.

    Attributes:
    model (optional): An external model object that may be used for substitution computations.
    state_frequencies (optional): A data structure (e.g., list or array) containing the
                                  frequencies of different states in the model.
    phi_matrix (optional): A matrix containing phi values, which could represent transition probabilities or
                           other relevant parameters in the model.
    rate_matrix (optional): A matrix containing rate values which could represent substitution rates or
                            other relevant parameters in the model.
    number_rates (optional): A scalar indicating the number of different rates in the model.
    category_rates (optional): A data structure (e.g., list or array) containing the rates of different
                               categories in the model.
    """

    def __init__(
        self,
        model=None,  # An optional external model object
        state_frequencies=None,  # Optional state frequencies
        phi_matrix=None,  # Optional phi matrix
        rate_matrix=None,  # Optional rate matrix
        number_rates=None,  # Optional number of rates
        category_rates=None,  # Optional category rates
    ) -> None:
        """Initializes a new SubstitutionModel instance with the provided parameters."""
        self.model = model  # Set the model attribute
        self.number_rates = number_rates  # Set the number_rates attribute
        self.rate_matrix = rate_matrix  # Set the rate_matrix attribute
        self.state_frequencies = (
            state_frequencies  # Set the state_frequencies attribute
        )
        self.phi_matrix = phi_matrix  # Set the phi_matrix attribute
        self.category_rates = category_rates  # Set the category_rates attribute


class IqTreeParser:
    """
    A class to parse IQ-TREE output files and extract relevant information for
    constructing a SubstitutionModel object.
    """

    def __init__(self, file_path=None):
        """
        Initializes the IqTreeParser with the path to the IQ-TREE file.

        Parameters:
        - file_path (str, optional): The path to the IQ-TREE file to be parsed.
        """

        self.file_content = []
        self.file_path = file_path
        if self.file_path:
            self.load_iqtree_file_content()

    def load_iqtree_file_content(self):
        """
        Loads the content of the .iqtree file into the class attribute 'file_content'.

        This method checks if the file exists and then attempts to read its content.
        If any issues arise during this process, an appropriate exception will be raised.

        Raises:
        - FileNotFoundError: If the specified file does not exist.
        - IOError: If there is an issue reading the file.
        """

        # Check if the file exists
        if not os.path.exists(self.file_path):
            raise FileNotFoundError(f"File '{self.file_path}' not found.")

        try:
            # Attempt to open and read the file content
            with open(self.file_path, "r") as file:
                self.file_content = file.readlines()
        except IOError as e:
            # Raise an IOError if there's a problem reading the file
            raise IOError(
                f"An error occurred while reading the file '{self.file_path}': {str(e)}"
            )

    def load_substitution_model(self) -> SubstitutionModel:
        """
        Load and parse the content of the iqtree file to form the substitution model.

        This method parses various components of the iqtree file such as state frequencies,
        rate matrices, model details, number of rate categories, and category rates to
        construct a SubstitutionModel object.

        Returns:
        - SubstitutionModel: An object containing the parsed details of the substitution model.
        """

        # Ensure the iqtree file content is loaded
        self.load_iqtree_file_content()
        # Parse the required components from the file content
        state_frequencies, phi_matrix = self.parse_state_frequencies()
        rate_matrix = self.parse_rate_matrices(state_frequencies)
        model = self.parse_substitution_model()
        number_rates = self.parse_number_rate_categories()
        category_rates = self.parse_category_rates() if number_rates > 1 else None

        # Return the constructed SubstitutionModel object
        return SubstitutionModel(
            model=model,
            state_frequencies=state_frequencies,
            phi_matrix=phi_matrix,
            rate_matrix=rate_matrix,
            number_rates=number_rates,
            category_rates=category_rates,
        )

    def parse_rate_parameters(self, dimension, model="GTR"):
        """
        Parse rate parameters from the file content to format the model string.

        This method identifies lines in the file content that contain rate parameters.
        It then formats and returns a model string with the parsed rates.

        Parameters:
        - dimension (int): The dimension of the rate matrix.
        - model (str, optional): The substitution model to be used. Defaults to "GTR".

        Returns:
        - str: The formatted model string with parsed rate parameters.
        """

        if not isinstance(dimension, int) or dimension <= 0:
            raise ValueError("Invalid dimension value provided.")

        # Set to store unique rate values
        rates = set()

        # Flag to indicate when rate parameters are found in the file
        found = False

        # Loop through each line in the file content
        for line in self.file_content:
            # Check for the marker indicating the start of rate parameters
            if "Rate parameter R:" in line:
                found = True
                continue
            if found:
                # Extract rate values from the line and add them to the rates set
                line_parts = line.split(":")
                if len(line_parts) == 2:
                    rates.add(line_parts[1].strip())

        if not found:
            raise ValueError("Rate parameters not found in the file.")

        # Split the model string on '+' to extract the base model (e.g., "GTR")
        splitted_model = model.split("+")
        # Combine the base model and rates list to form the final model string
        model_with_rates_token = f"{splitted_model[0]}{{{' ,'.join(rates)}}}"

        return model_with_rates_token

    def parse_state_frequencies(self):
        """
        Parse the stationary distribution pi from a given .iqtree file path.

        Returns
        - frequencies (directory): The stationary distribution.
        - phi_matrix (np.array): The stationary distribution pi with values filled in the diagonal.
        """
        index = next(
            (
                idx
                for idx, line in enumerate(self.file_content)
                if "State frequencies:" in line
            ),
            None,
        )

        if index is None:
            raise ValueError("'State frequencies:' not found in file.")

        # Dynamically determine the dimension 'n'
        start_idx = next(
            (
                idx
                for idx, line in enumerate(self.file_content)
                if "Rate matrix Q:" in line
            ),
            None,
        )

        if start_idx is None:
            raise ValueError(
                "'Rate matrix Q:' not found in file. Determination of dimension not possible."
            )

        # Detect the number of matrix rows based on numeric entries
        n = 0
        current_idx = start_idx + 2  # Adjusting to start from matrix values
        while current_idx < len(self.file_content) and re.search(
            r"(\s*-?\d+\.\d+\s*)+", self.file_content[current_idx]
        ):
            n += 1
            current_idx += 1

        # Initialize an empty dictionary to hold the frequencies
        frequencies = {}

        # Parse the state frequencies
        for idx, line in enumerate(self.file_content):
            # Parse the state frequencies (empirical counts)
            if "State frequencies: (empirical counts from alignment)" in line:
                try:
                    for i in range(n):
                        # Split the line on " = " into a key and a value, and add them to the frequencies dictionary
                        key, value = self.file_content[idx + i + 2].split(" = ")
                        frequencies[key] = float(
                            value
                        )  # convert value to float before storing
                except (IndexError, ValueError) as e:
                    raise Exception(
                        f"Error while parsing empirical state frequencies. Exception: {e}"
                    )

            # If "equal frequencies" is in the log content, return a pseudo dictionary with equal frequencies
            elif "State frequencies: (equal frequencies)" in line:
                for i in range(n):
                    key = "key_" + str(i)
                    frequencies[key] = float(1 / n)

        phi_matrix = np.diag(list(frequencies.values()))
        return frequencies, phi_matrix

    def parse_rate_matrices(self, state_frequencies):
        """

        Parse the rate parameters R  .iqtree file path and determine
        the rate matrix Q using the rate parameters and stationary distribution

        Returns:
        - rate_matrix (np.array): The parsed rate matrix Q.
        """

        # Dynamically determine the matrix dimension 'n'
        start_idx = next(
            (
                idx
                for idx, line in enumerate(self.file_content)
                if "Rate matrix Q:" in line
            ),
            None,
        )

        if start_idx is None:
            raise ValueError("'Rate matrix Q:' not found in file.")

        # Detect the number of matrix rows based on numeric entries
        n = 0
        current_idx = start_idx + 2  # Adjusting to start from matrix values
        while current_idx < len(self.file_content) and re.search(
            r"(\s*-?\d+\.\d+\s*)+", self.file_content[current_idx]
        ):
            n += 1
            current_idx += 1

        # Find the start and end indices of the rate parameter
        start_index = -1
        end_index = -1

        for i, line in enumerate(self.file_content):
            if line.strip().startswith("Rate parameter R:"):
                start_index = i + 1
            elif line.strip().startswith("State frequencies"):
                end_index = i
                break

        if start_index == -1 or end_index == -1:
            raise ValueError("Rate parameter not found in the log file.")

        # Get the lines in between the start and end markers
        parameter = self.file_content[start_index + 1 : end_index]
        rates = {}

        # Loop through the lines we've extracted
        for line in parameter:
            # If the line is not empty after removing leading/trailing whitespace
            if line.strip():
                # Split the line on " = " into a key and a value, and add them to the frequencies dictionary
                key, value = line.strip().split(":")
                rates[key] = float(value)  # convert value to float before storing

        only_rates = list(rates.values())
        list_state_freq = list(state_frequencies.values())
        idx = 0
        rate_matrix = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1, n):
                rate_matrix[i, j] = only_rates[idx] * list_state_freq[j]
                rate_matrix[j, i] = only_rates[idx] * list_state_freq[i]
                idx += 1
        average_rate = 0
        for i in range(n):
            rate_matrix[i, i] = -sum(rate_matrix[i,])
            average_rate = average_rate + rate_matrix[i, i] * list_state_freq[i]
        rate_matrix = -rate_matrix / average_rate

        return rate_matrix

    def parse_number_rate_categories(self):
        """
        Parse the number of rate categories from the substitution model string.

        Parameters:
        - model (str): The substitution model string from the file.

        Returns:
        - rate (int): The number of rate categories parsed from the model string.
        """
        index = next(
            (
                idx
                for idx, line in enumerate(self.file_content)
                if "Model of rate heterogeneity:" in line
            ),
            None,
        )
        if index is None:
            raise ValueError("'Model of rate heterogeneity:' not found in file.")

        if "Uniform" in self.file_content[index]:
            number_rate_categories = 1
        else:
            # initialize substrings
            sub1 = "with "
            sub2 = " categories"
            # getting index of substrings

            idx1 = self.file_content[index].index(sub1)
            idx2 = self.file_content[index].index(sub2)

            res = ""
            # getting elements in between
            for idx3 in range(idx1 + len(sub1), idx2):
                res = res + self.file_content[index][idx3]

            number_rate_categories = int(res)
        return number_rate_categories

    def parse_rate_from_model(self, model):
        """
        Parse the number of rate categories from the substitution model string.

        Parameters:
        - model (str): The substitution model string from the file.

        Returns:
        - rate (int): The number of rate categories parsed from the model string.
        """
        plus_g_index = model.find("+G")
        plus_r_index = model.find("+R")

        if plus_g_index != -1 and plus_r_index != -1:
            raise ValueError("Cannot use +G and +R")

        if plus_g_index != -1:
            rate_start_index = plus_g_index + 2
        elif plus_r_index != -1:
            rate_start_index = plus_r_index + 2
        else:
            return 1  # default number_rates = 1 if no +G or +R model

        try:
            # Extract the substring after e.g.'+G'
            number = model[rate_start_index:]

            # Parse the extracted substring as an integer
            if "{" in number:
                # e.g. +G{0.9} will fix the Gamma shape parameter (alpha)to 0.9
                # discrete Gamma model: default 4 rate categories
                rate = 4
            else:
                # number of rate categories
                rate = int(number)

            return rate
        except ValueError:
            # If '+G' is not found or the number after '+G' is not a valid integer
            # Return None or an appropriate value for error handling
            raise ValueError("Could not parse the substitution model from the file.")

    def parse_substitution_model(self) -> str:
        """
        Parses the file content to extract the substitution model.

        The function searches for lines that contain specific keywords indicative of the substitution model.
        Once found, it extracts and returns the model as a string.

        Returns:
        - str: The extracted substitution model. If not found, raises a ValueError.
        """
        # Keywords indicating the presence of the substitution model in a line
        keywords = ["Best-fit model according to BIC:", "Model of substitution:"]
        for line in self.file_content:
            # Check if the line contains any of the keywords
            if any(keyword in line for keyword in keywords):
                model_string = line.split(":")[1].strip()

                # If a valid model string is found, return it
                if model_string:
                    return model_string

        # If the loop completes without returning, raise an error
        raise ValueError("Substitution model not found in the file content.")

    def parse_category_rates(self):
        """
        Parses the category rates from the file content and returns them in a structured format.

        The function identifies the table of category rates in the file content and extracts
        the category, relative rate, and proportion for each row in the table.

        Returns:
        - dict: A dictionary containing category rates. Each key represents a category and the
                associated value is another dictionary with details of that category.
        """

        # Get the number of rate categories from another method
        number_rates = self.parse_number_rate_categories()

        # Dictionary to store parsed table data
        table_data = {}

        # Variables to track the start and end of the table in the file content
        start_index = -1
        end_index = -1

        # Find the start and end indices of the table in a single pass
        for i, line in enumerate(self.file_content):
            stripped_line = line.strip()
            if stripped_line.startswith("Category"):
                start_index = i + 1
            elif start_index != -1 and stripped_line.startswith(f"{number_rates}"):
                end_index = i + 1
                break

        # Error handling in case the table isn't found
        if start_index == -1 or end_index == -1:
            raise ValueError("Table not found in the log file.")

        # Extract and parse the table rows
        for line in self.file_content[start_index:end_index]:
            parts = line.split()
            if len(parts) >= 3:
                category, relative_rate, proportion = parts[:3]
                if category != "0":
                    table_data[f"p{category}"] = {
                        "Category": category,
                        "Relative_rate": float(relative_rate),
                        "Proportion": float(proportion),
                    }

        return table_data


if __name__ == "__main__":
    iq_tree_parser = IqTreeParser(
        "./test/Clemens/example_3/sim-JC+G-AC1-AG1-AT1-CG1-CT1-GT1-alpha1.2-taxa64-len1000bp-bla0.01-blb0.8-blc0.2-rep01.fasta.iqtree"
    )
    substitution_model = iq_tree_parser.load_substitution_model_from_iqtree_file()

    print(substitution_model.model)
    print(substitution_model.state_frequencies)
    print(substitution_model.phi_matrix)
    print(substitution_model.rate_matrix)
>>>>>>> Stashed changes
