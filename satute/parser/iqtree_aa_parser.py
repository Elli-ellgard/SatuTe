import re
import numpy as np
from typing import List, Tuple
from satute.parser.base_parser import BaseParser
from satute.exceptions import ModelNotFoundError, InvalidModelNameError
from satute.models.amino_acid_models import (
    AMINO_ACID_RATE_MATRIX,
    AA_STATE_FREQUENCIES,
    get_core_model,
)


class AminoAcidParser(BaseParser):
    def parse_state_frequencies(self) -> Tuple[List[float], np.ndarray]:
        """
        Parses the state frequencies from the file content.

        Returns:
        - list: A list of state frequencies as floats.
        """
        capture_state_frequencies: bool = False
        state_frequencies: List[float] = []

        for line in self.file_content:
            # Check for the state frequencies header in the file
            if "State frequencies:" in line:
                capture_state_frequencies = True
                continue

            # Capture the frequencies after the header is found
            if capture_state_frequencies:
                if "Model of rate heterogeneity:" in line:
                    break
                if line.strip() != "":
                    frequency = line.split("=")[1].strip()
                    state_frequencies.append(float(frequency))
        return state_frequencies, np.diag(state_frequencies)

    @staticmethod
    def get_aa_rate_matrix(current_substitution_model: str) -> np.ndarray:
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

        if (
            core_model not in AMINO_ACID_RATE_MATRIX
            or core_model not in AA_STATE_FREQUENCIES
        ):
            raise ModelNotFoundError(core_model)

        rate_matrix_params = AMINO_ACID_RATE_MATRIX[core_model]
        rate_matrix_eq = AA_STATE_FREQUENCIES[core_model]
        rate_matrix = AminoAcidParser.create_rate_matrix_with_input(
            20, rate_matrix_params, rate_matrix_eq
        )

        return np.array(rate_matrix)

    def parse_substitution_rates(self) -> str:
        # Flag to indicate if the next lines contain the Q matrix
        capture_substitution_rate: bool = False
        # List to store the rows of the Q matrix
        string_based_substitution_rates: str = []

        for line in self.file_content:
            # Check for the Q matrix header in the file
            if "Substitution parameters" in line:
                capture_substitution_rate = True
                continue

            # Capture the matrix after the header is found
            if capture_substitution_rate:
                if "State frequencies:" in line:
                    break
                if line.strip() != "":
                    string_based_substitution_rates.append(line.strip())
        return "\n".join(string_based_substitution_rates[0:18])

    def get_aa_state_frequency_substitution_models(
        self, substitution_model: str
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Extracts amino acid state frequencies and creates a diagonal matrix for a substitution model.

        Args:
            substitution_model (str): The substitution model string, potentially with extensions.

        Returns:
            tuple: A tuple containing:
                - A numpy array of state frequencies.
                - A diagonal numpy matrix of these frequencies.

        Raises:
            ModelNotFoundError: If the core model is not found in AA_STATE_FREQUENCIES.
            InvalidModelNameError: If the core model cannot be extracted.
        """
        # Regular expression to extract the core model name
        core_model = get_core_model(substitution_model)

        if core_model not in AA_STATE_FREQUENCIES:
            raise ModelNotFoundError(core_model)

        frequencies = np.array(AA_STATE_FREQUENCIES[core_model])
        return frequencies, np.diag(frequencies)

    @staticmethod
    def create_rate_matrix_with_input(
        matrix_size: int, input_string: str, eq: List[float]
    ) -> np.array:
        # Split the string into lines
        # input_string = input_string.replace(" ", "")
        lines = input_string.split("\n")

        # Initialize a matrix with 1's for off-diagonal elements
        rate_matrix = [
            [0 if i != j else 0 for j in range(matrix_size)] for i in range(matrix_size)
        ]

        for i, row in enumerate(lines):
            for j, col in enumerate(row.split()):
                rate_matrix[i + 1][j] = float(col)

        # Mirror the lower triangle to the upper triangle
        for i in range(matrix_size):
            for j in range(i):
                rate_matrix[j][i] = rate_matrix[i][j]

        # Multiply the rates by the equilibrium distribution
        rate_matrix = rate_matrix @ np.diag(eq)

        for i in range(matrix_size):
            rate_matrix[i][i] = -np.sum(rate_matrix[i])

        return np.array(rate_matrix)

    @staticmethod
    def normalize_stationary_distribution_aa(frequencies: List[float]) -> List[float]:
        sum_freqs = sum(frequencies)
        if sum_freqs == 1:
            # Valid stationary distribution
            return frequencies
        else:
            # Normalize frequencies list with new values
            return [freq / sum_freqs for freq in frequencies]
