from satute_rate_categories_and_alignments import read_alignment_file
from satute_rate_categories_and_alignments import (
    split_msa_into_rate_categories_in_place,
)
from satute_util_new import parse_file_to_data_frame
from satute_trees_and_subtrees import rescale_branch_lengths, parse_newick_file
from satute_rate_categories_and_alignments import parse_category_rates
from direction_based_saturation_test import (
    name_nodes_by_level_order,
    calculate_partial_likelihoods_for_sites,
    RateMatrix,
)
from satute_rate_categories_and_alignments import (
    read_alignment_file,
    split_msa_into_rate_categories_in_place,
)
from satute_util_new import (
    parse_file_to_data_frame,
    spectral_decomposition_without_path,
)
from satute_trees_and_subtrees import rescale_branch_lengths, parse_newick_file
from satute_repository import parse_rate_matrices
import pandas as pd
from satute_test_statistic_using_partial_likelihood import calculate_test_statistic


def single_rate_analysis(t, alignment, rate_matrix, array_eigenvectors, multiplicity):
    partial_likelihood_per_site_storage = calculate_partial_likelihoods_for_sites(
        t, alignment, rate_matrix
    )
    for edge in partial_likelihood_per_site_storage.keys():
        left_partial_likelihood = pd.DataFrame(
            partial_likelihood_per_site_storage[edge]["left"]
        )
        right_partial_likelihood = pd.DataFrame(
            partial_likelihood_per_site_storage[edge]["right"]
        )

        branch_type = "external"
        if "Node" in edge[1]:
            branch_type = "internal"

        print(
            calculate_test_statistic(
                multiplicity,
                array_eigenvectors,
                left_partial_likelihood,
                right_partial_likelihood,
                4,
                branch_type,
                alpha=0.05,
            )
        )


def multiple_rate_analysis(
    t,
    category_rates_factors,
    rate_matrix,
    array_eigenvectors,
    multiplicity,
    per_rate_category_alignment,
):
    for rate, alignment in per_rate_category_alignment.items():
        relative_rate = category_rates_factors[rate]["Relative_rate"]

        if len(alignment) > 0:
            rescaled_tree = rescale_branch_lengths(t, relative_rate)
            partial_likelihood_per_site_storage = (
                calculate_partial_likelihoods_for_sites(
                    rescaled_tree, alignment, rate_matrix
                )
            )

            for edge in partial_likelihood_per_site_storage.keys():
                left_partial_likelihood = pd.DataFrame(
                    partial_likelihood_per_site_storage[edge]["left"]
                )
                right_partial_likelihood = pd.DataFrame(
                    partial_likelihood_per_site_storage[edge]["right"]
                )

                branch_type = "external"
                if "Node" in edge[1]:
                    branch_type = "internal"

                print(
                    calculate_test_statistic(
                        multiplicity,
                        array_eigenvectors,
                        left_partial_likelihood,
                        right_partial_likelihood,
                        4,
                        branch_type,
                        alpha=0.05,
                    )
                )


def main():
    alignment_file = "./Clemens/example_3/sim-JC+G-AC1-AG1-AT1-CG1-CT1-GT1-alpha1.2-taxa64-len1000bp-bla0.01-blb0.8-blc0.2-rep01.fasta"
    alignment = read_alignment_file(alignment_file)
    t = parse_newick_file(
        "./Clemens/example_3/sim-JC+G-AC1-AG1-AT1-CG1-CT1-GT1-alpha1.2-taxa64-len1000bp-bla0.01-blb0.8-blc0.2-rep01.fasta.treefile"
    )
    t = name_nodes_by_level_order(t)
    number_rate = 4

    RATE_MATRIX, psi_matrix = parse_rate_matrices(
        4,
        "./Clemens/example_3/sim-JC+G-AC1-AG1-AT1-CG1-CT1-GT1-alpha1.2-taxa64-len1000bp-bla0.01-blb0.8-blc0.2-rep01.fasta",
    )
    rate_matrix = RateMatrix(RATE_MATRIX)

    array_eigenvectors, multiplicity = spectral_decomposition_without_path(
        RATE_MATRIX, psi_matrix
    )

    if number_rate == 1:
        single_rate_analysis(
            t, alignment, rate_matrix, array_eigenvectors, multiplicity
        )
    else:
        site_probability = parse_file_to_data_frame(
            "./Clemens/example_3/sim-JC+G-AC1-AG1-AT1-CG1-CT1-GT1-alpha1.2-taxa64-len1000bp-bla0.01-blb0.8-blc0.2-rep01.fasta.siteprob"
        )
        per_rate_category_alignment = split_msa_into_rate_categories_in_place(
            site_probability, alignment
        )
        category_rates_factors = parse_category_rates(
            "./Clemens/example_3/sim-JC+G-AC1-AG1-AT1-CG1-CT1-GT1-alpha1.2-taxa64-len1000bp-bla0.01-blb0.8-blc0.2-rep01.fasta.iqtree"
        )
        multiple_rate_analysis(
            t,
            category_rates_factors,
            rate_matrix,
            array_eigenvectors,
            multiplicity,
            per_rate_category_alignment,
        )


if __name__ == "__main__":
    main()
