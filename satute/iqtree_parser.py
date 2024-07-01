# -*- coding: utf-8 -*-
import os
import re
import numpy as np
import pandas as pd
from enum import Enum
from typing import List, Dict

from satute.exceptions import InvalidModelNameError, ModelNotFoundError
from satute.dna_models import NOT_ACCEPTED_DNA_MODELS
from satute.substitution_model import SubstitutionModel
from satute.amino_acid_models import get_aa_state_frequency_substitution_models, normalize_stationary_distribution_aa
from satute.dna_models import LIE_DNA_MODELS

from satute.amino_acid_models import (
    AMINO_ACID_RATE_MATRIX,
    AMINO_ACID_MODELS,
    AA_STATE_FREQUENCIES,
    NOT_ACCEPTED_AA_MODELS,
    create_rate_matrix_with_input    
)

class ModelType(Enum):
    DNA = "DNA"
    PROTEIN = "Protein"

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

    def create_rate_matrix_for_lie_markov_models(self,substitution_rates: Dict[str, float], state_frequencies : List[float]) -> List[List[float]]:
        """
        Create a rate matrix from the given substitution rates.

        Parameters:
        - substitution_rates (Dict[str, float]): A dictionary of substitution rates.

        Returns:
        - rate_matrix (List[List[float]]): The resulting rate matrix.
        """
        rates = list(substitution_rates.values())
        num_rates = len(rates)
        matrix_size = int(num_rates**0.5) + 1  # Calculate the size of the matrix (square root of num_rates plus one)

        rate_matrix: List[List[float]] = []
        rate_index = 0
        for i in range(matrix_size):
            row: List[float] = []
            row_sum = 0
            for j in range(matrix_size):
                if i == j:
                    row.append(0)  # Placeholder for the diagonal
                else:
                    rate = rates[rate_index] * state_frequencies[j]
                    row.append(rate)
                    row_sum += rate
                    rate_index += 1

            row[i] = -row_sum  # Set the diagonal element
            rate_matrix.append(row)
        
        
        for i in range(matrix_size):
            row_sum = 0 
            for j in range(matrix_size):
                if i!= j:
                    row_sum += rate_matrix[i][j]

            for j in range(matrix_size):
                rate_matrix[i][j] = rate_matrix[i][j] / row_sum
                                    
        rate_matrix = np.array(rate_matrix)
        return rate_matrix

    def check_if_lie_model(self, current_substitution_model: str):        
        for lie_model in LIE_DNA_MODELS:
            if lie_model in current_substitution_model:
                return True        
        return False

    def load_substitution_model(self) -> SubstitutionModel:
        """
        Load and parse the content of the iqtree file to form the substitution model.

        This method parses various components of the iqtree file such as state frequencies,
        rate matrices, model details, number of rate categories, and category rates to
        construct a SubstitutionModel object.

        Returns:
        - SubstitutionModel: An object containing the parsed details of the substitution model.
        """

        self.load_iqtree_file_content()        
        current_substitution_model : str = self.parse_substitution_model()
        self.check_model(current_substitution_model)
        self.model_type : ModelType = self.get_model_type(current_substitution_model)
        
        state_frequencies = []
        phi_matrix = []
        rate_matrix = []
        number_rates = 0
        category_rates = []
        pre_computed_q_matrix = []
        
        if self.model_type == ModelType.DNA:
            
            pre_computed_q_matrix  = self.parse_nucleotide_q_matrix()            
            
            if self.check_if_lie_model(current_substitution_model):
                
                dict_state_frequencies, phi_matrix = self.parse_state_frequencies()
                substitution_rates = self.parse_substitution_rates()
                state_frequencies = list(dict_state_frequencies.values())
                rate_matrix = self.create_rate_matrix_for_lie_markov_models(substitution_rates, list(dict_state_frequencies.values()))
                
            else:                
                
                # Parse the rate matrix and stationary distribution for the DNA Substitution Model
                dict_state_frequencies, phi_matrix = self.parse_state_frequencies()
                state_frequencies = list(dict_state_frequencies.values())
                rate_matrix = self.construct_rate_matrix(dict_state_frequencies)
                
        else:
            
            current_substitution_model = current_substitution_model.upper()
            # Parse the rate matrix and stationary distribution for the Protein Substitution Model
            state_frequencies, phi_matrix = get_aa_state_frequency_substitution_models(
                current_substitution_model
            )
            state_frequencies = normalize_stationary_distribution_aa(state_frequencies)
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
            precomputed_q_matrix=pre_computed_q_matrix,
        )

    def check_model(self, model: str):
        """
        Check if the model is one of the not accepted DNA or Protein models.

        Parameters:
        - model (str): The model string to be checked.

        Raises:
        - ValueError: If the model is not accepted for analysis.
        """
        for dna_model in NOT_ACCEPTED_DNA_MODELS:
            if dna_model in model:
                raise ValueError(
                    f"The DNA model '{dna_model}' is not accepted for analysis because it is non-reversible."
                )

        for aa_model in NOT_ACCEPTED_AA_MODELS:
            if aa_model in model:
                raise ValueError(
                    f"The protein model '{aa_model}' is not accepted for analysis."
                )

    def get_model_type(self, model: str) -> ModelType:
    
        # Check if it is a protein model
        model_upper = model.upper()
        for amino_acid_substitution_model in AMINO_ACID_MODELS:
            if amino_acid_substitution_model in model_upper:
                return ModelType.PROTEIN

        # Default to DNA if no conditions above are met
        return ModelType.DNA

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
            if "State frequencies: (empirical counts from alignment)" in line or 'State frequencies: (estimated with maximum likelihood)' in line:
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
            if "State frequencies: (equal frequencies)" in line:
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
    def find_dimension_by_rate_matrix_parsing(start_index: int, file_content: str) -> int:
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
        n = IqTreeParser.find_dimension_by_rate_matrix_parsing(
            start_index=start_index, file_content=self.file_content
        )
        
        substitution_rates = self.parse_substitution_rates()
        substitution_rates = list(substitution_rates.values())
        list_state_freq = list(state_frequencies.values())
        
        
                
        return self.build_rate_matrix(n=n, rates=substitution_rates, list_state_freq=list_state_freq)

    def parse_nucleotide_q_matrix(self) -> np.array:
        # Flag to indicate if the next lines contain the Q matrix
        capture_matrix: bool = False
        # List to store the rows of the Q matrix
        string_based_q_matrix: List[str] = []

        for line in self.file_content:
            # Check for the Q matrix header in the file
            if 'Rate matrix Q:' in line:
                capture_matrix = True
                continue
                
            # Capture the matrix after the header is found
            if capture_matrix:
                if 'Model of rate heterogeneity:' in line:
                    break                
                if line.strip() != '':        
                    string_based_q_matrix.append(line.strip())
                    
        # Process captured lines to format them into a proper matrix
        
        # Extract the numeric values from each string
        numeric_values = []
        for line in string_based_q_matrix:
            parts = line.split()
            numeric_values.append([float(part) for part in parts[1:]])
        return np.array(numeric_values)

    def parse_substitution_rates(self):
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
        substitution_rates = {}
        for line in self.file_content[start_index:end_index]:
            line = line.strip()
            if line:
                key, value = line.split(":")
                substitution_rates[key.strip()] = float(
                    value.strip()
                )  # Ensure keys and values are cleanly extracted

        return substitution_rates

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
        n: int, rates: List[float], state_frequencies: List[float]
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

        Raises:
        - ModelNotFoundError: If the core model name is not found in the AMINO_ACID_RATE_MATRIX dictionary.
        - InvalidModelNameError: If the core model name cannot be extracted.

        Returns:
        - np.ndarray: The constructed amino acid rate matrix.
        """
        # Regular expression to extract the core model name
        core_model_match = re.match(r"^[^\+\{]+", current_substitution_model)
        if core_model_match:
            core_model = core_model_match.group()
        else:
            raise InvalidModelNameError(current_substitution_model)

        if core_model not in AMINO_ACID_RATE_MATRIX or core_model not in AA_STATE_FREQUENCIES:
            raise ModelNotFoundError(core_model)

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

    def parse_number_rate_categories(self) -> int:
        """
        Parse the number of rate categories from the substitution model string.

        Returns:
        - rate (int): The number of rate categories parsed from the model string.
        """
        index = next(
            (idx for idx, line in enumerate(self.file_content)
             if "Model of rate heterogeneity:" in line),
            None,
        )

        if index is None:
            raise ValueError("'Model of rate heterogeneity:' not found in file.")
        
        line = self.file_content[index]
        
        if "Uniform" in line:
            return 1
        if "with" in line and "categories" in line:
            start = line.index("with ") + len("with ")
            end = line.index(" categories")
            return int(line[start:end])
        if "Invar" in line:
            return 1
        raise ValueError("Unexpected format for 'Model of rate heterogeneity:' line.")

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

    def parse_category_rates(self)->dict:
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

def parse_file_to_data_frame(file_path: str) -> pd.DataFrame:
    try:
        # Read the file into a dataframe
        df = pd.read_csv(file_path, delimiter="\t")
        return df

    except FileNotFoundError:
        raise Exception(f"File not found: {file_path}")

def valid_stationary_distribution(frequencies: Dict[str, float]) -> Dict[str, float]:
    sum_frequencies = sum(frequencies.values())
    if sum_frequencies == 1:
        # Valid stationary distribution
        return frequencies
    else:
        # Normalize frequencies dictionary with new values
        for key in frequencies:
            frequencies[key] /= sum_frequencies
        return frequencies

