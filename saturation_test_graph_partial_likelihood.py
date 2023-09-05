from satute_rate_categories_and_alignments import read_alignment_file
from satute_rate_categories_and_alignments import (
    split_msa_into_rate_categories_in_place,
)
from satute_util_new import parse_file_to_data_frame
from satute_trees_and_subtrees import parse_newick_file
from satute_rate_categories_and_alignments import parse_category_rates
from satute_direction_based_saturation_test import (
    name_nodes_by_level_order,
    RateMatrix,
    single_rate_analysis,
    multiple_rate_analysis,
    
)

from satute_rate_categories_and_alignments import (
    read_alignment_file,
    split_msa_into_rate_categories_in_place,
)
from satute_util_new import parse_file_to_data_frame, spectral_decomposition
from satute_trees_and_subtrees import parse_newick_file
from satute_repository import parse_rate_matrices_from_file


def main():
    alignment_file = "./Clemens/example_3/sim-JC+G-AC1-AG1-AT1-CG1-CT1-GT1-alpha1.2-taxa64-len1000bp-bla0.01-blb0.8-blc0.2-rep01.fasta"
    alignment = read_alignment_file(alignment_file)
    t = parse_newick_file(
        "./Clemens/example_3/sim-JC+G-AC1-AG1-AT1-CG1-CT1-GT1-alpha1.2-taxa64-len1000bp-bla0.01-blb0.8-blc0.2-rep01.fasta.treefile"
    )
    t = name_nodes_by_level_order(t)
    number_rate = 4

    RATE_MATRIX, psi_matrix = parse_rate_matrices_from_file(
        4,
        "./Clemens/example_3/sim-JC+G-AC1-AG1-AT1-CG1-CT1-GT1-alpha1.2-taxa64-len1000bp-bla0.01-blb0.8-blc0.2-rep01.fasta",
    )
    rate_matrix = RateMatrix(RATE_MATRIX)

    state_frequencies = parse_state_frequencies(
        "./Clemens/example_3/sim-JC+G-AC1-AG1-AT1-CG1-CT1-GT1-alpha1.2-taxa64-len1000bp-bla0.01-blb0.8-blc0.2-rep01.fasta.iqtree",
        dimension=4,
    )

    (
        array_left_eigenvectors,
        array_right_eigenvectors,
        multiplicity,
    ) = spectral_decomposition(RATE_MATRIX, psi_matrix)

    if number_rate == 1:
        single_rate_analysis(
            t,
            alignment,
            rate_matrix,
            state_frequencies,
            array_left_eigenvectors,
            array_right_eigenvectors,
            multiplicity,
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

        


if __name__ == "__main__":
    main()
