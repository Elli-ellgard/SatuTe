import re 
import numpy as np
import os
from satute_repository import (
    valid_stationary_distribution,
)

def read_file_content(file_path):
    try:
        with open(file_path, 'r') as file:
            content = file.read()
        return content
    except FileNotFoundError:
        print("File not found.")
        return None
    except Exception as e:
        print("An error occurred:", e)
        return None
    
   

def parse_state_frequencies(file_content):
        """
        Parse the stationary distribution pi from a given .iqtree file path.

        Returns
        - frequencies (directory): The stationary distribution.
        - phi_matrix (np.array): The stationary distribution pi with values filled in the diagonal.
        """
        index = next(
            (
                idx
                for idx, line in enumerate(file_content)
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
                for idx, line in enumerate(file_content)
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
        while current_idx < len(file_content) and re.search(
            r"(\s*-?\d+\.\d+\s*)+", file_content[current_idx]
        ):
            n += 1
            current_idx += 1

        # Initialize an empty dictionary to hold the frequencies
        frequencies = {}

        # Parse the state frequencies
        for idx, line in enumerate(file_content):
            # Parse the state frequencies (empirical counts)
            if "State frequencies: (empirical counts from alignment)" in line:
                try:
                    for i in range(n):
                        # Split the line on " = " into a key and a value, and add them to the frequencies dictionary
                        key, value = file_content[idx + i + 2].split(" = ")
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




def parse_rate_matrices(file_content, state_frequencies):
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
                for idx, line in enumerate(file_content)
                if "Rate matrix Q:" in line
            ),
            None,
        )

        if start_idx is None:
            raise ValueError("'Rate matrix Q:' not found in file.")

        # Detect the number of matrix rows based on numeric entries
        n = 0
        current_idx = start_idx + 2  # Adjusting to start from matrix values
        while current_idx < len(file_content) and re.search(
            r"(\s*-?\d+\.\d+\s*)+", file_content[current_idx]
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



def test_get_rate_matrix(file_path):
    content = read_file_content(file_path)
    if content:
        print(content)
        stationary_distribution, phi_matrix = parse_state_frequencies(content)
        print(stationary_distribution)

        


if __name__ == "__main__":
    # set paths to IQ-TREE and Python executable
    path_iqtree = "iqtree"
    path_python = "python3"
    
    # examples
    data_dir_path = "./tests/data/data_dna/model"
    output_dir_path =  "./tests/test_results/"
    iqtree_file1 = os.path.join(data_dir_path, "example1.iqtree")

    test_get_rate_matrix(iqtree_file1)
    
