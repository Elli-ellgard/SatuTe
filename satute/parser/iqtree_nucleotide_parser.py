from satute.parser.base_parser import BaseParser
from typing import List, Dict
import numpy as np
import re


class NucleotideParser(BaseParser):
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
            if (
                "State frequencies: (empirical counts from alignment)" in line
                or "State frequencies: (estimated with maximum likelihood)" in line
            ):
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

    def construct_rate_matrix(self, state_frequencies: List[float]) -> np.array:
        start_index = self.find_line_index(self.file_content, "Rate matrix Q:")
        n = self.find_dimension_by_rate_matrix_parsing(
            start_index=start_index, file_content=self.file_content
        )
        substitution_rates = self.parse_substitution_rates()
        substitution_rates = list(substitution_rates.values())
        list_state_freq = list(state_frequencies.values())
        return self.build_rate_matrix(
            n=n, rates=substitution_rates, list_state_freq=list_state_freq
        )

    def parse_q_matrix(self) -> np.array:
        # Flag to indicate if the next lines contain the Q matrix
        capture_matrix: bool = False
        # List to store the rows of the Q matrix
        string_based_q_matrix: List[str] = []

        for line in self.file_content:
            # Check for the Q matrix header in the file
            if "Rate matrix Q:" in line:
                capture_matrix = True
                continue

            # Capture the matrix after the header is found
            if capture_matrix:
                if "Model of rate heterogeneity:" in line:
                    break
                if line.strip() != "":
                    string_based_q_matrix.append(line.strip())

        # Process captured lines to format them into a proper matrix

        # Extract the numeric values from each string
        numeric_values = []
        for line in string_based_q_matrix:
            parts = line.split()
            numeric_values.append([float(part) for part in parts[1:]])
        return np.array(numeric_values)

    def parse_substitution_rates(self) -> Dict[str, float]:
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
        substitution_rates: Dict[str, float] = {}
        for line in self.file_content[start_index:end_index]:
            line = line.strip()
            if line:
                key, value = line.split(":")
                substitution_rates[key.strip()] = float(
                    value.strip()
                )  # Ensure keys and values are cleanly extracted

        return substitution_rates

    def normalize_rate_matrix(
        self, rate_matrix: np.ndarray, state_frequencies: List[float], n: int
    ) -> np.ndarray:
        average_rate = 0
        for i in range(n):
            rate_matrix[i, i] = -np.sum(rate_matrix[i, :]) + rate_matrix[i, i]
            average_rate += rate_matrix[i, i] * state_frequencies[i]
        if average_rate == 0:
            raise ValueError(
                "Division by zero in normalization due to zero average rate."
            )
        return -rate_matrix / average_rate

    def assemble_rate_matrix_from_rates_and_frequencies(
        self, n: int, rates: List[float], state_frequencies: List[float]
    ) -> np.ndarray:
        expected_number_of_rates = n * (n - 1) // 2
        if len(rates) != expected_number_of_rates:
            raise ValueError(
                f"Expected {expected_number_of_rates} rates, got {len(rates)}."
            )
        rate_matrix = np.zeros((n, n))
        idx = 0
        for i in range(n):
            for j in range(i + 1, n):
                rate_matrix[i, j] = rates[idx] * state_frequencies[j]
                rate_matrix[j, i] = rates[idx] * state_frequencies[i]
                idx += 1
        return rate_matrix

    def build_rate_matrix(
        self, rates: List[float], list_state_freq: List[float], n: int
    ) -> np.ndarray:
        if len(list_state_freq) != n:
            raise ValueError("Length of state_frequencies must match the dimension n.")

        rate_matrix = self.assemble_rate_matrix_from_rates_and_frequencies(
            n,
            rates,
            list_state_freq,
        )

        normalized_rate_matrix = self.normalize_rate_matrix(
            rate_matrix, list_state_freq, n
        )

        return normalized_rate_matrix

    def create_rate_matrix_for_lie_markov_models(
        self, substitution_rates: Dict[str, float], state_frequencies: List[float]
    ) -> List[List[float]]:
        """
        Create a rate matrix from the given substitution rates.

        Parameters:
        - substitution_rates (Dict[str, float]): A dictionary of substitution rates.

        Returns:
        - rate_matrix (List[List[float]]): The resulting rate matrix.
        """
        rates = list(substitution_rates.values())
        num_rates = len(rates)
        matrix_size = (
            int(num_rates**0.5) + 1
        )  # Calculate the size of the matrix (square root of num_rates plus one)

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
                if i != j:
                    row_sum += rate_matrix[i][j]

            for j in range(matrix_size):
                rate_matrix[i][j] = rate_matrix[i][j] / row_sum

        rate_matrix = np.array(rate_matrix)
        return rate_matrix


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
