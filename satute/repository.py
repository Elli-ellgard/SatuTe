# -*- coding: utf-8 -*-
import os
import numpy as np
import pandas as pd
import re
from enum import Enum
from typing import List

from satute.amino_acid_models import get_aa_state_frequency_substitution_models
from satute.dna_model import NOT_ACCEPTED_DNA_MODELS
from satute.amino_acid_models import (
    AMINO_ACID_RATE_MATRIX,
    create_rate_matrix_with_input,
    AMINO_ACID_MODELS,
    AA_STATE_FREQUENCIES,
    NOT_ACCEPTED_AA_MODELS,
)


class ModelType(Enum):
    DNA = "DNA"
    PROTEIN = "Protein"


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

    def __init__(self, file_path: str = None):
        """
        Initializes the IqTreeParser with the path to the IQ-TREE file.

        Parameters:
        - file_path (str, optional): The path to the IQ-TREE file to be parsed.
        """

        self.file_content = []
        self.model_type = ModelType.DNA
        self.file_path = file_path
        if self.file_path:
            self.load_iqtree_file_content()

    def load_iqtree_file_content(self) -> list[str]:
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

        current_substitution_model = self.parse_substitution_model()
        self.check_model(current_substitution_model)
        self.model_type = self.get_model_type(current_substitution_model)
        self.load_iqtree_file_content()

        if self.model_type == ModelType.DNA:
            # Parse the rate matrix and stationary distribution for the DNA Substitution Model
            dict_state_frequencies, phi_matrix = self.parse_state_frequencies()
            state_frequencies = dict_state_frequencies.values()
            rate_matrix = self.construct_rate_matrix(dict_state_frequencies)
        else:
            # Parse the rate matrix and stationary distribution for the Protein Substitution Model
            state_frequencies, phi_matrix = get_aa_state_frequency_substitution_models(
                current_substitution_model
            )
            rate_matrix = self.get_aa_rate_matrix(current_substitution_model)

        number_rates = self.parse_number_rate_categories()
        category_rates = self.parse_category_rates() if number_rates > 1 else None

        return SubstitutionModel(
            model=current_substitution_model,
            state_frequencies=state_frequencies,
            phi_matrix=phi_matrix,
            rate_matrix=rate_matrix,
            number_rates=number_rates,
            category_rates=category_rates,
        )

    def check_model(self, model: str):
        # Check if the model is one of the not accepted DNA models
        for dna_model in NOT_ACCEPTED_DNA_MODELS:
            if dna_model in model:
                # Raise an exception signaling that the model is not accepted
                raise ValueError(
                    f"The model '{dna_model}' is not accepted for analysis, because it is non-reversible."
                )
        # Check if the model is one of the not accepted Protein models
        for aa_model in NOT_ACCEPTED_AA_MODELS:
            if aa_model in model:
                # Raise an exception signaling that the model is not accepted
                raise ValueError(
                    f"The model '{aa_model}' is not accepted for analysis."
                )

    def get_model_type(self, model: str) -> ModelType:
        # Check if it is a protein model
        for amino_acid_substitution_model in AMINO_ACID_MODELS:
            if amino_acid_substitution_model in model:
                return ModelType.PROTEIN

        # Default to DNA if no conditions above are met
        return ModelType.DNA

    def parse_rate_parameters(self, dimension: int, model: str = "GTR"):
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
                    frequencies = valid_stationary_distribution(frequencies)

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

    @staticmethod
    def find_line_index(lines: List[str], search_string: str) -> int:
        """
        Returns the index of the first line that contains the given search string.
        Raises an exception if the search string is not found.

        Args:
        - lines (List[str]): The list of lines to search through.
        - search_string (str): The string to search for.

        Returns:
        - int: The index of the first line containing the search string.

        Raises:
        - ValueError: If the search string is not found in any line.
        """
        for idx, line in enumerate(lines):
            if search_string in line:
                return idx
        raise ValueError(
            f"Search string '{search_string}' not found in the provided lines."
        )

    @staticmethod
    def find_dimension_by_rate_matrix_parsing(start_index, file_content) -> int:
        # Detect the number of matrix rows based on numeric entries
        n = 0
        current_idx = start_index + 2  # Adjusting to start from matrix values
        while current_idx < len(file_content) and re.search(
            r"(\s*-?\d+\.\d+\s*)+", file_content[current_idx]
        ):
            n += 1
            current_idx += 1
        return n

    def construct_rate_matrix(self, state_frequencies: List[float]):
        """

        Parse the rate parameters R  .iqtree file path and determine
        the rate matrix Q using the rate parameters and stationary distribution

        Returns:
        - rate_matrix (np.array): The parsed rate matrix Q.
        """
        start_index = IqTreeParser.find_line_index(self.file_content, "Rate matrix Q:")
        # Detect the number of matrix rows based on numeric entries
        n = IqTreeParser.find_dimension_by_rate_matrix_parsing(
            start_index=start_index, file_content=self.file_content
        )
        rates_dict = self.extract_rate_parameters()

        rates = list(rates_dict.values())
        list_state_freq = list(state_frequencies.values())
        return self.build_rate_matrix(n=n, rates=rates, list_state_freq=list_state_freq)

    def extract_rate_parameters(self):
        """
        Extracts rate parameters from the content of the log file. Assumes the rates are listed
        between the lines starting with 'Rate parameter R:' and 'State frequencies'.

        Returns:
        - dict: A dictionary of rate parameters with their corresponding values.

        Raises:
        - ValueError: If the section containing rate parameters is not properly defined.
        """
        start_index = -1
        end_index = -1
        # Identify the start and end indices for the rate parameters section
        for i, line in enumerate(self.file_content):
            stripped_line = line.strip()
            if stripped_line.startswith("Rate parameter R:"):
                start_index = i + 1  # Start reading rates from the next line
            elif stripped_line.startswith("State frequencies"):
                end_index = i
                break

        if start_index == -1 or end_index == -1:
            raise ValueError(
                "Rate parameter section is incomplete or missing in the log file."
            )

        # Extract the parameters and store them in a dictionary
        rates = {}
        for line in self.file_content[start_index:end_index]:
            line = line.strip()
            if line:
                key, value = line.split(":")
                rates[key.strip()] = float(
                    value.strip()
                )  # Ensure keys and values are cleanly extracted

        return rates

    def build_rate_matrix(
        self, rates: List[float], list_state_freq: List[float], n: int
    ):
        if len(list_state_freq) != n:
            raise ValueError("Length of state_frequencies must match the dimension n.")

        rate_matrix = IqTreeParser.assemble_rate_matrix_from_rates_and_frequencies(
            n,
            rates,
            list_state_freq,
        )

        normalized_rate_matrix = IqTreeParser.normalize_rate_matrix(
            rate_matrix, list_state_freq, n
        )

        return normalized_rate_matrix

    @staticmethod
    def normalize_rate_matrix(
        rate_matrix: np.ndarray, state_frequencies: List[float], n: int
    ) -> np.ndarray:
        """
        Normalizes the rate matrix by adjusting its scale based on the weighted average rate.

        This function ensures the sum of each row in the rate matrix is zero and normalizes
        the matrix by the average rate calculated as the weighted sum of diagonal elements,
        using state frequencies as weights.

        Args:
        - rate_matrix (np.ndarray): The rate matrix to be normalized, shape (n, n).
        - state_frequencies (List[float]): List of frequencies for each state, length n.
        - n (int): The dimension of the rate matrix.

        Returns:
        - np.ndarray: The normalized rate matrix with adjusted scale.

        Raises:
        - ValueError: If average rate results in zero, potentially causing division by zero.
        """
        average_rate = 0
        for i in range(n):
            rate_matrix[i, i] = -np.sum(rate_matrix[i, :]) + rate_matrix[i, i]
            average_rate += rate_matrix[i, i] * state_frequencies[i]

        if average_rate == 0:
            raise ValueError(
                "Division by zero in normalization due to zero average rate."
            )
        return -rate_matrix / average_rate

    @staticmethod
    def assemble_rate_matrix_from_rates_and_frequencies(
        n, rates: List[float], state_frequencies: List[float]
    ):
        """
        Fills the upper and lower triangles of a rate matrix based on the provided rates
        and state frequencies.

        Args:
        - n (int): The dimension of the rate matrix.
        - rates (list[float]): List of rate values for interactions between states. Should have n*(n-1)/2 elements.
        - state_frequencies (list[float]): List of state frequencies for each state.
        Returns:
        - np.ndarray: A symmetric rate matrix with both upper and lower triangles filled.

        Raises:
        - ValueError: If the length of rates does not match the expected number for a  symmetric matrix of dimension n.
        """
        expected_number_of_rates = n * (n - 1) // 2
        if len(rates) != expected_number_of_rates:
            raise ValueError(
                f"Expected {expected_number_of_rates} rates, got {len(rates)}."
            )

        # Initialize the rate matrix
        rate_matrix = np.zeros((n, n))
        idx = 0

        # Fill the upper and lower triangles of the matrix
        for i in range(n):
            for j in range(i + 1, n):
                rate_matrix[i, j] = rates[idx] * state_frequencies[j]  # Upper triangle
                rate_matrix[j, i] = rates[idx] * state_frequencies[i]  # Lower triangle
                idx += 1
        return rate_matrix

    def get_aa_rate_matrix(self, current_substitution_model: str) -> np.ndarray:
        """
        Retrieves and constructs the amino acid rate matrix for a given substitution model.

        This method sanitizes the model name to remove potential model extensions or parameters
        indicated by "+" or "{" symbols. It then looks up the core model in the
        AMINO_ACID_RATE_MATRIX dictionary to obtain the parameters required to construct the
        rate matrix. The matrix is constructed using these parameters and returned as a NumPy array.

        Parameters:
        - current_substitution_model (str): The substitution model string, which may include extensions or parameters.

        Returns:
        - np.ndarray: The amino acid rate matrix as a NumPy array.

        Raises:
        - KeyError: If the core model name is not found in the AMINO_ACID_RATE_MATRIX dictionary.
        """
        # Regular expression to extract the core model name
        core_model_match = re.match(r"^[^\+\{]+", current_substitution_model)
        if core_model_match:
            core_model = core_model_match.group()
        else:
            raise ValueError(
                f"Could not extract core model from: {current_substitution_model}"
            )

        if core_model not in AMINO_ACID_RATE_MATRIX:
            raise KeyError(f"Model '{core_model}' not found in AMINO_ACID_RATE_MATRIX.")

        # THIS WILL LEAD TO A BUG

        rate_matrix_params = AMINO_ACID_RATE_MATRIX[core_model]
        rate_matrix_eq = AA_STATE_FREQUENCIES[core_model]
        rate_matrix = create_rate_matrix_with_input(
            20, rate_matrix_params, rate_matrix_eq
        )

        return np.array(rate_matrix)

    def parse_aa_rate_matrix(self) -> np.matrix:
        """
        Parses the amino acid rate matrix from the file content.

        This function searches through the file content for a specific line
        that indicates the start of the rate matrix ("Rate matrix Q:").
        It then parses the following lines as the matrix, converting the
        values to floats and constructing a 2D numpy matrix. If the resulting
        matrix is not 20x20, a ValueError is raised.

        Raises:
            ValueError: If the "Rate matrix Q:" line is not found in the file.
            ValueError: If the resulting matrix is not 20x20.

        Returns:
            np.matrix: A 2D numpy matrix representing the parsed rate matrix.
        """

        # Find the line where the rate matrix starts
        start = None
        for i, line in enumerate(self.file_content):
            if line.startswith("Rate matrix Q:"):
                start = i + 2
                break

        if start is None:
            raise ValueError("Rate matrix Q not found in file.")

        # Parse the matrix using a list comprehension
        matrix = []
        for line in self.file_content[start:]:
            if line.strip() == "":
                break
            row = line.split(maxsplit=1)[1].split()
            matrix.append([float(x) for x in row])

        # Convert to numpy matrix and check its shape
        matrix = np.matrix(matrix, dtype=float)
        if matrix.shape != (20, 20):
            raise ValueError(
                "The rate matrix is not 20x20. Its shape is {}".format(matrix.shape)
            )

        # Return the numpy matrix
        return matrix

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

    def parse_rate_from_model(self, model: str):
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


def parse_substitution_model(file_path: str) -> str:
    """
    Parse the substitution model from an IQ-TREE log file.

    This function reads an IQ-TREE log file and extracts the substitution model
    based on specific lines containing "Best-fit model according to BIC:" or
    "Model of substitution:". The function returns the extracted model as a string.

    Parameters:
    - file_path (str): The path to the IQ-TREE log file.

    Returns:
    - str: The parsed substitution model.

    Raises:
    - ValueError: If the expected model strings are not found in the file.
    - ValueError: If there's an error reading the file.

    Example:
    >>> parse_substitution_model("path_to_iqtree_log.txt")
    'TEST_MODEL_1'

    Note:
    The function expects the IQ-TREE log file to contain specific lines indicating
    the substitution model. If these lines are not found or if there's an issue
    reading the file, a ValueError is raised.

    """
    try:
        with open(file_path, "r") as file:
            content = file.read()
            for line in content.splitlines():
                if ("Best-fit model according to BIC:" in line) or (
                    "Model of substitution:" in line
                ):
                    model_string = line.split(":")[1].strip()
                    if model_string:
                        return model_string
            raise ValueError("Expected model strings not found in the file.")
    except IOError:
        raise ValueError("Could not read the file.")


def parse_rate_from_cli_input(model: str) -> int:
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
        # Extract the substring after e.g.'+G'
        number = model[rate_start_index:]

        # Parse the extracted substring as an integer
        if "{" in number:
            # e.g. +G{0.9} will fix the Gamma shape parameter (alpha)to 0.9
            # discrete Gamma model: default 4 rate categories
            rate = 4
            return rate
        else:
            if number and str(number).isnumeric():
                # number of rate categories
                rate = int(number)
                return rate
            else:
                return "AMBIGUOUS"
    except ValueError:
        # If '+G' is not found or the number after '+G' is not a valid integer
        # Return None or an appropriate value for error handling
        raise ValueError("Could not parse the substitution model from the file.")


def parse_file_to_data_frame(file_path) -> pd.DataFrame:
    try:
        # Read the file into a dataframe
        df = pd.read_csv(file_path, delimiter="\t")

        return df

    except FileNotFoundError:
        raise Exception(f"File not found: {file_path}")


def valid_stationary_distribution(frequencies: List[float]) -> List[float]:
    sum_freqs = sum(frequencies.values())
    if sum_freqs == 1:
        # Valid stationary distribution
        return frequencies
    else:
        # Update frequencies dictionary with new values
        keys = list(frequencies.keys())
        for i in range(len(keys)):
            frequencies[keys[i]] /= sum_freqs
        return frequencies


if __name__ == "__main__":
    iq_tree_parser = IqTreeParser("./test/amino_acids/alignment_aa.phy.iqtree")
    substitution_model = iq_tree_parser.load_substitution_model()
