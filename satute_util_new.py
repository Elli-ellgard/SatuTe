import os
import pandas as pd
import scipy
import scipy.linalg
import numpy as np
from rich import print
import logging


def initialize_logger(log_file):
    # Create a logger instance
    logger = logging.getLogger("my_logger")
    logger.setLevel(logging.DEBUG)

    # Create a file handler and set its level to DEBUG
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.DEBUG)

    # Create a console handler and set its level to INFO
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)

    # Create a formatter and set it for the handlers
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    file_handler.setFormatter(formatter)
    console_handler.setFormatter(formatter)

    # Add the handlers to the logger
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)

    return logger


""" ## SPECTRAL DECOMPOSITION OF THE RATE MATRIX"""

def spectral_decomposition(rate_matrix, psi_matrix):
    """ 
        Then psi_matrix := Diag(pi). Recall that matrix Q is reversible iff M:= psi_matrix^1/2 x Q x psi_matrix^{-1/2} is symmetric.
        For a real symmetric matrix M, its eigenvectors can be chosen to be an orthonormal basis of R^n 
    """
    M = scipy.linalg.fractional_matrix_power(psi_matrix, +1 / 2) @ rate_matrix
    M = M @ scipy.linalg.fractional_matrix_power(psi_matrix, -1 / 2)

    """ eigendecomposition of matrix M"""
    lamb, w = np.linalg.eig(M)  # Compute the eigenvalues and eigenvectors.
    idx = lamb.argsort()[::-1]  # Order from large to small.
    lamb = lamb[idx]  # Order the eigenvalues (large to small).
    w = w[:, idx]  # Order the eigenvectors according to the eigenvalues"""

    # the first one should be the eigenvalue 0 in lamb, why are we doing the following?
    lamb_nozero = []  # list of eigenvalues without 0
    for i in lamb:
        if i > 0.00999 or i < -0.00999:
            lamb_nozero.append(i)

    max_lambda = max(lamb_nozero)  # dominant non-zero eigenvalue
    # get the indices of the dominant non-zero eigenvalue in lamb taking numerical inaccuracies into account and identical values
    index = []
    for i in range(len(lamb)):
        lambda_it = lamb[i]
        if abs(lambda_it - max_lambda) < 0.01:
            index.append(i)

    multiplicity = len(index)  # multiplicity of the dominant non-zero eigenvalue
    array_right_eigenvectors = (
        []
    )  # list of right eigenvectors for the dominant non-zero eigenvalue
    array_left_eigenvectors = (
        []
    )  # list of left eigenvectors for the dominant non-zero eigenvalue
    for i in range(multiplicity):
        # calculate the right eigenvectors for the dominant non-zero eigenvalue
        v1 = scipy.linalg.fractional_matrix_power(psi_matrix, -1 / 2) @ w[:, index[i]]
        array_right_eigenvectors.append(v1)
        # calculate the left eigenvectors for the dominant non-zero eigenvalue
        h1 = scipy.linalg.fractional_matrix_power(psi_matrix, +1 / 2) @ w[:, index[i]]
        array_left_eigenvectors.append(h1)

    return array_left_eigenvectors, array_right_eigenvectors, multiplicity


""" ## GENERATE STRUCTURE FOR SUBTREES"""

def parse_file_to_data_frame(file_path):
    try:
        # Read the file into a dataframe
        df = pd.read_csv(file_path, delimiter="\t")

        return df

    except FileNotFoundError:
        raise Exception(f"File not found: {file_path}")


def write_results_and_newick_tree(
    results_list, newick_string, path_folder, chosen_rate, c_sTwoSequence, T
):
    """
    This function writes the saturation branches data to a file and then appends the saturation information
    newick string to the same file.

    :param results_list: List of results data to be written to file.
    :param newick_string: Newick formatted string representing the tree.
    :param path_folder: Path to the folder where results will be saved.
    :param chosen_rate: Chosen rate parameter to be included in the filename.
    :param c_sTwoSequence: Saturation coherence between two sequences.
    :param T: ETE Tree instance.
    """

    # Convert the results_list into a pandas dataframe
    saturation_branches_data_frame = pd.DataFrame(results_list)

    # Save the dataframe as a tab-separated CSV file
    saturation_branches_data_frame.to_csv(
        f"{path_folder}/resultsRate{chosen_rate}.satute.csv",
        header=True,
        index=None,
        sep="\t",
        mode="w",
    )

    print(saturation_branches_data_frame)

    # Generate a newick string with saturation information
    # saturation_information_newick_string = map_values_to_newick(
    #    results_list, newick_string
    # )

    # Open the results file in append mode
    with open(
        f"{path_folder}/resultsRate{chosen_rate}.satute.csv", "a"
    ) as satute_result_file:
        # Write additional information to the file
        satute_result_file.write(
            "\n\nThe T2T status uses as threshold the saturation coherence between two sequences, which is  {:6.4f}".format(
                c_sTwoSequence
            )
        )

        satute_result_file.write(
            "\n\nFor better reference, this is the reconstructed tree topology :\n\n"
        )

        # Write the ascii representation of the tree to the file
        satute_result_file.write(
            T.copy("newick").get_ascii(attributes=["name", "label", "distance"])
        )

        # Tree without saturation values
        satute_result_file.write(f"\n\n Tree with saturation values: {newick_string}")
        # Write the saturation information newick string to the file
        # satute_result_file.write(f"\n\n Tree with saturation values: {saturation_information_newick_string}")

