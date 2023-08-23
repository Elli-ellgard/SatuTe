import os
import pandas as pd
from Bio import AlignIO
import scipy
import scipy.linalg
import numpy as np
from rich import print
from multiprocessing import Pool
from satute_repository import (
    parse_state_frequencies,
    parse_rate_matrices,
)
import glob
from satute_repository import parse_output_state_frequencies
from satute_trees_and_subtrees import (
    name_nodes_by_level_order,
    get_leaves,
    branch_lengths,
    rescale_branch_lengths,
    generate_subtree_pair,
    get_all_subtrees,
    parse_newick_file,
)
from satute_calculate_likelihoods import (
    generate_output_state_file_for_cherry,
    generate_output_state_file_for_external_branch,
)

from satute_rate_categories_and_alignments import (
    filter_alignment_by_ids,
    read_alignment_file,
    split_msa_into_rate_categories,
    parse_category_rates,
)

from satute_test_statistic import (
    calculate_test_statistic,
    calculate_likelihood_ratio_test,
)


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


def spectral_decomposition(n, path):
    rate_matrix, psi_matrix = parse_rate_matrices(n, path)

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
    array_eigenvectors = (
        []
    )  # list of right eigenvectors for the dominant non-zero eigenvalue
    for i in range(multiplicity):
        # calculate the right eigenvectors for the dominant non-zero eigenvalue
        v1 = scipy.linalg.fractional_matrix_power(psi_matrix, -1 / 2) @ w[:, index[i]]
        array_eigenvectors.append(v1)
        """# calculate the left eigenvectors for the dominant non-zero eigenvalue
        h1 = (
            scipy.linalg.fractional_matrix_power(psi_matrix, +1 / 2)
            @ w[:, index[i]]
            ) """

    return array_eigenvectors, multiplicity


""" ## GENERATE STRUCTURE FOR SUBTREES"""


def write_subtree_and_sub_alignments(
    generated_subtree_pairs, alignment, state_frequencies, rate_matrix, path_prefix="./"
):
    for i, subtree_pair in enumerate(generated_subtree_pairs):
        first_subtree = subtree_pair["trees"][0]
        second_subtree = subtree_pair["trees"][1]

        first_subtree_leaves = get_leaves(first_subtree)
        first_sub_alignment = filter_alignment_by_ids(alignment, first_subtree_leaves)

        second_subtree_leaves = get_leaves(second_subtree)
        second_sub_alignment = filter_alignment_by_ids(alignment, second_subtree_leaves)

        if len(first_subtree.get_descendants()) + 1 == 1:
            first_subtree_dir = f"{path_prefix}subtrees/branch_{i}_subtree_one_leaf/"

        elif len(first_subtree.get_descendants()) + 1 == 3:
            first_subtree_dir = f"{path_prefix}subtrees/branch_{i}_subtree_one_cherry/"
        else:
            first_subtree_dir = f"{path_prefix}subtrees/branch_{i}_subtree_one/"

        os.makedirs(first_subtree_dir, exist_ok=True)

        # Write the first subtree
        subtree_writer = open(f"{first_subtree_dir}subtree.treefile", "w")
        subtree_writer.write(first_subtree.write(format=1))
        subtree_writer.close()

        # Write the alignment for the first subtree
        AlignIO.write(first_sub_alignment, f"{first_subtree_dir}subtree.fasta", "fasta")

        if len(second_subtree.get_descendants()) + 1 == 1:
            second_subtree_dir = f"{path_prefix}subtrees/branch_{i}_subtree_two_leaf/"

        elif len(second_subtree.get_descendants()) + 1 == 3:
            second_subtree_dir = f"{path_prefix}subtrees/branch_{i}_subtree_two_cherry/"

        else:
            second_subtree_dir = f"{path_prefix}subtrees/branch_{i}_subtree_two/"

        os.makedirs(second_subtree_dir, exist_ok=True)

        # Write the second subtree
        subtree_writer = open(f"{second_subtree_dir}subtree.treefile", "w")
        subtree_writer.write(second_subtree.write(format=1))
        subtree_writer.close()

        # Write the alignment for the second subtree
        AlignIO.write(
            second_sub_alignment, f"{second_subtree_dir}subtree.fasta", "fasta"
        )

        if len(first_subtree.get_descendants()) + 1 == 1:
            generate_output_state_file_for_external_branch(
                first_sub_alignment,
                first_subtree_dir,
                state_frequencies,
            )
        elif len(first_subtree.get_descendants()) + 1 == 3:
            generate_output_state_file_for_cherry(
                first_sub_alignment,
                first_subtree,
                first_subtree_dir,
                state_frequencies,
                rate_matrix,
            )

        if len(second_subtree.get_descendants()) + 1 == 1:
            generate_output_state_file_for_external_branch(
                second_sub_alignment,
                second_subtree_dir,
                state_frequencies,
            )
        elif len(second_subtree.get_descendants()) + 1 == 3:
            generate_output_state_file_for_cherry(
                second_sub_alignment,
                second_subtree,
                second_subtree_dir,
                state_frequencies,
                rate_matrix,
            )


def parse_file_to_data_frame(file_path):
    try:
        # Read the file into a dataframe
        df = pd.read_csv(file_path, delimiter="\t")

        return df

    except FileNotFoundError:
        raise Exception(f"File not found: {file_path}")


def generate_write_subtree_pairs_and_msa(
    number_rates,
    t,
    file_path,
    msa_file_name,
    category_rate,
    state_frequencies,
    rate_matrix,
):
    # Write subtree pairs to files, assuming the function is defined elsewhere
    if number_rates == 1:
        try:
            # Call function to get all subtrees from t, assuming the function is defined elsewhere
            subtrees = get_all_subtrees(t)

            ## Discard the first subtree, assuming we don't need it
            # subtrees = subtrees[1:]

            # Generate subtree pairs, assuming the function is defined elsewhere
            generated_subtree_pairs = generate_subtree_pair(subtrees, t)
            alignment = read_alignment_file(msa_file_name)

            write_subtree_and_sub_alignments(
                generated_subtree_pairs,
                alignment,
                state_frequencies,
                rate_matrix,
                f"./{file_path}/",
            )

        except Exception as e:
            print(f"Error occurred during the first iteration: {e}")
    else:
        # Iterate from 0 to number_rates (inclusive)
        for i in range(1, number_rates + 1):
            rate = category_rate[i - 1]["Relative_rate"]
            # Rescale the branch lengths
            rescaled_tree = rescale_branch_lengths(t, rate)
            # Write subtree pairs to files in the new directory
            # Call function to get all subtrees from t, assuming the function is defined elsewhere
            subtrees = get_all_subtrees(rescaled_tree)

            # Generate subtree pairs, assuming the function is defined elsewhere
            generated_subtree_pairs = generate_subtree_pair(subtrees, rescaled_tree)

            sub_alignment = read_alignment_file(
                f"./{file_path}/subsequence{i}/rate.fasta"
            )

            write_subtree_and_sub_alignments(
                generated_subtree_pairs,
                sub_alignment,
                state_frequencies,
                rate_matrix,
                path_prefix=f"./{file_path}/subsequence{i}/",
            )


def fetch_subdirectories(directory, prefix):
    pattern = os.path.join(directory, prefix + "*")
    subdirectories = glob.glob(pattern)
    return subdirectories


def build_tree_test_space(number_rates, msa_file_name, target_directory, dimension=4):
    site_probability = parse_file_to_data_frame(f"{msa_file_name}.siteprob")
    t = name_nodes_by_level_order(parse_newick_file(f"{msa_file_name}.treefile"))

    category_rate = ""
    state_frequencies = parse_state_frequencies(
        f"{msa_file_name}.iqtree", dimension=dimension
    )

    (
        rate_matrix,
        phi_matrix,
    ) = parse_rate_matrices(dimension, msa_file_name)

    valid_category_rates = []

    if number_rates != 1:
        valid_category_rates = split_msa_into_rate_categories(
            site_probability, target_directory, msa_file_name
        )

        category_rate = parse_category_rates(f"{msa_file_name}.iqtree")

    generate_write_subtree_pairs_and_msa(
        number_rates,
        t,
        target_directory,
        msa_file_name,
        category_rate,
        state_frequencies,
        rate_matrix,
    )

    return valid_category_rates


def run_saturation_test_for_branches_and_categories(
    input_directory,
    target_directory,
    number_rates,
    t,
    dimension,
    alpha=0.01,
    valid_category_rates=[],
):
    # number of branches
    branch_count = len(t.get_descendants())

    """ get the right eigenvector(s) of the dominate non-zero eigenvalue"""
    array_eigenvectors, multiplicity = spectral_decomposition(
        dimension, input_directory
    )

    vector_branches, vector_distances = branch_lengths(t)

    results_list = {}
    # special case: +G or +R model
    # test for branch saturation for each valid  rate category of the model
    if number_rates != 1:
        for rate in valid_category_rates:
            results_list[rate] = []

            for branch_index in range(0, branch_count):
                subtree_directories = fetch_subdirectories(
                    f"{target_directory}/subsequence{rate}/subtrees",
                    f"branch_{branch_index}",
                )
                """ determine branch type (internal or external) and 
                    set the subdirectory corresponding to the leaf as the left subtree
                """
                if subtree_directories[0].find("leaf") != -1:
                    left_subtree_dir = subtree_directories[0]
                    right_subtree_dir = subtree_directories[1]
                    branch_type = "external"
                elif subtree_directories[1].find("leaf") != -1:
                    left_subtree_dir = subtree_directories[1]
                    right_subtree_dir = subtree_directories[0]
                    branch_type = "external"
                else:
                    left_subtree_dir = subtree_directories[0]
                    right_subtree_dir = subtree_directories[1]
                    branch_type = "internal"

                # read the posterior probabilities of the left subtree from .state file
                posterior_probabilities_left_subtree = parse_output_state_frequencies(
                    f"{left_subtree_dir}/output.state"
                )
                # read the posterior probabilities of the right subtree from .state file
                posterior_probabilities_right_subtree = parse_output_state_frequencies(
                    f"{right_subtree_dir}/output.state"
                )

                # Calculation of the test for branch saturation
                results = run_saturation_test_for_branch(
                    multiplicity,
                    array_eigenvectors,
                    branch_type,
                    dimension,
                    alpha,
                    posterior_probabilities_left_subtree,
                    posterior_probabilities_right_subtree,
                    vector_distances[branch_index],
                    input_directory,
                )

                results["vector_branches"] = vector_branches[branch_index]
                results_list[rate].append(results)

            print(f"Results for rate {rate}")

            write_results_and_newick_tree(
                results_list[rate],
                t.write(format=1),
                target_directory,
                rate,
                results_list[rate][0]["c_s_two_sequence"],
                t,
            )
    # model of evolution does not include +G and +R
    else:
        results_list[number_rates] = []
        for branch_index in range(0, branch_count):
            subtree_directories = fetch_subdirectories(
                f"{target_directory}/subtrees",
                f"branch_{branch_index}",
            )

            """ determine branch type (internal or external) and 
                set the subdirectory corresponding to the leaf as the left subtree
            """
            if subtree_directories[0].find("leaf") != -1:
                left_subtree_dir = subtree_directories[0]
                right_subtree_dir = subtree_directories[1]
                branch_type = "external"
            elif subtree_directories[1].find("leaf") != -1:
                left_subtree_dir = subtree_directories[1]
                right_subtree_dir = subtree_directories[0]
                branch_type = "external"
            else:
                left_subtree_dir = subtree_directories[0]
                right_subtree_dir = subtree_directories[1]
                branch_type = "internal"

            # read the posterior probabilities of the left subtree from .state file
            posterior_probabilities_left_subtree = parse_output_state_frequencies(
                f"{left_subtree_dir}/output.state"
            )
            # # read the posterior probabilities of the right subtree from .state file
            posterior_probabilities_right_subtree = parse_output_state_frequencies(
                f"{right_subtree_dir}/output.state"
            )

            # Calculation of the  test for branch saturation
            results = run_saturation_test_for_branch(
                multiplicity,
                array_eigenvectors,
                branch_type,
                dimension,
                alpha,
                posterior_probabilities_left_subtree,
                posterior_probabilities_right_subtree,
                vector_distances[branch_index],
                input_directory,
            )

            results_list[number_rates].append(results)

            write_results_and_newick_tree(
                results_list[number_rates],
                t.write(format=1),
                target_directory,
                1,
                results_list[number_rates][0]["c_s_two_sequence"],
                t,
            )


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


def run_saturation_test_for_branch(
    multiplicity,
    array_eigenvectors,
    branch_type,
    dimension,
    alpha,
    posterior_probabilities_left_subtree,
    posterior_probabilities_right_subtree,
    branch_length,
    input_directory,
):
    (
        delta,
        c_s,
        c_s_two_sequence,
        p_value,
        result_test,
        result_test_tip2tip,
    ) = calculate_test_statistic(
        multiplicity,
        array_eigenvectors,
        posterior_probabilities_left_subtree,
        posterior_probabilities_right_subtree,
        dimension,
        branch_type,
        alpha,
    )

    (lr_test, result_lr_test) = calculate_likelihood_ratio_test(
        input_directory,
        branch_length,
        posterior_probabilities_left_subtree,
        posterior_probabilities_right_subtree,
        dimension,
        alpha,
    )

    return {
        "delta": delta,
        "c_s": c_s,
        "c_s_two_sequence": c_s_two_sequence,
        "p_value": p_value,
        "result_test": result_test,
        "result_test_tip2tip": result_test_tip2tip,
        "lr_test": lr_test,
        "result_lr_test": result_lr_test,
    }
