from ete3 import Tree
from satute_direction_based_partial_likelihood import (
    name_nodes_by_level_order,
    calculate_partial_likelihoods_for_sites,
    RATE_MATRIX,
    NUCLEOTIDE_CODE_VECTOR,
    RateMatrix
)
from satute_rate_categories_and_alignments import read_alignment_file
import numpy as np
from satute_rate_categories_and_alignments import (
    split_msa_into_rate_categories_in_place,
)
from satute_util_new import parse_file_to_data_frame
from pathlib import Path
from satute_trees_and_subtrees import rescale_branch_lengths
from satute_rate_categories_and_alignments import parse_category_rates
import json
import pandas as pd


def dump_dict_to_json(dict_obj, filename):
    with open(filename, "w") as f:
        json.dump(dict_obj, f, indent=4)


if __name__ == "__main__":
    alignment_file = "./Clemens/example_3/sim-JC+G-AC1-AG1-AT1-CG1-CT1-GT1-alpha1.2-taxa64-len1000bp-bla0.01-blb0.8-blc0.2-rep01.fasta"
    newick_string = "(t000000:0.0067892257,t000001:0.0067787574,(((t000010:0.0059673375,((((((t010000:0.0030563760,t010001:0.0114093850):0.0094182376,(t010010:0.0131446398,t010011:0.0097395536):0.0133621566):0.0104108104,(((t011000:0.0124958339,t011001:0.0125011474):0.0075917776,(t011010:0.0051495665,t011011:0.0124730837):0.0069341413):0.0078912673,((t011100:0.0104118838,t011101:0.0072227008):0.0073307877,(t011110:0.0144514225,t011111:0.0055786880):0.0102245774):0.0059268453):0.0184390656):0.0061628417,(t010110:0.0082616199,t010111:0.0146332627):0.0101990740):0.0052544609,(t010100:0.0092991497,t010101:0.0114804795):0.0092491114):0.6992169228,(((((t100000:0.0114103154,t100001:0.0051613214):0.0061697358,(t100010:0.0136268730,t100011:0.0114221532):0.0072470758):0.0155632848,((t100100:0.0081218605,t100101:0.0073837399):0.0130036860,(t100110:0.0071915054,t100111:0.0125350507):0.0131619893):0.0090043100):0.0007243217,(((t101000:0.0106452037,t101001:0.0182293876):0.0045033981,(t101010:0.0121747946,t101011:0.0121664259):0.0105869819):0.0092355725,((t101100:0.0078255361,t101101:0.0176067822):0.0088049086,(t101110:0.0061914059,t101111:0.0114288359):0.0119823991):0.0077233513):0.0115530529):0.7027436378,(((((t110000:0.0152109080,t110001:0.0070318894):0.0081739309,(t110010:0.0101504845,t110011:0.0117088941):0.0091799218):0.0134467581,((t110100:0.0050989755,t110101:0.0103895458):0.0158105135,(t110110:0.0072259592,t110111:0.0124929711):0.0072567017):0.0115489293):0.0176840292,((t111000:0.0127615540,t111001:0.0059373002):0.0117927352,(t111010:0.0162749028,t111011:0.0080293329):0.0114795491):0.0181492686):0.0106254068,((t111100:0.0082902032,t111101:0.0114231618):0.0065691489,(t111110:0.0101567116,t111111:0.0058159250):0.0047846372):0.0005525044):0.8384907952):0.4960190049):0.7890571203):0.0055032299,t000011:0.0125092154):0.0125564734,(((t000100:0.0129052942,t000101:0.0081677003):0.0128058021,(t000110:0.0061739221,t000111:0.0124881756):0.0088603366):0.0105251504,(((t001000:0.0124981985,t001001:0.0103935537):0.0089686048,(t001010:0.0135912957,t001011:0.0093498130):0.0044762602):0.0161461211,((t001100:0.0093389059,t001101:0.0125044394):0.0064576918,(t001110:0.0104262512,t001111:0.0135413564):0.0135893360):0.0090709753):0.0287882859):0.0040229738):0.0188329058);"

    t = Tree(newick_string, format=1)

    t = name_nodes_by_level_order(t)

    alignment = read_alignment_file(alignment_file)

    print("Reading Rate Matrix")

    RATE_MATRIX = np.array([[-3, 1, 1, 1], [1, -3, 1, 1], [1, 1, -3, 1], [1, 1, 1, -3]])
    number_rate = 2
    rate_matrix = RateMatrix(RATE_MATRIX)

    if number_rate == 1:
        print("Number Rate: ", number_rate)
        partial_likelihood_per_site_storage = calculate_partial_likelihoods_for_sites(
            t, alignment
        )
    else:

        site_probability = parse_file_to_data_frame(
            "./Clemens/example_3/sim-JC+G-AC1-AG1-AT1-CG1-CT1-GT1-alpha1.2-taxa64-len1000bp-bla0.01-blb0.8-blc0.2-rep01.fasta.siteprob"
        )
        valid_category_rates = split_msa_into_rate_categories_in_place(
            site_probability, "./Clemens/example_3", alignment
        )
        category_rates_factors = parse_category_rates(
            "./Clemens/example_3/sim-JC+G-AC1-AG1-AT1-CG1-CT1-GT1-alpha1.2-taxa64-len1000bp-bla0.01-blb0.8-blc0.2-rep01.fasta.iqtree"
        )

        for valid_category_rate in valid_category_rates:
            print("Calculating for rate ", valid_category_rate)
            rate = category_rates_factors[valid_category_rate - 1]["Relative_rate"]
            rescaled_tree = rescale_branch_lengths(t, rate)
            vertical_sub_alignment = read_alignment_file(
                "./Clemens/example_3/subsequence"
                + str(valid_category_rate)
                + "/rate.fasta"
            )
            partial_likelihood_per_site_storage = (
                calculate_partial_likelihoods_for_sites(
                    rescaled_tree, vertical_sub_alignment, rate_matrix
                )
            )

            for edge, values in partial_likelihood_per_site_storage.items():

                left_partial_likelihood = pd.DataFrame(
                    partial_likelihood_per_site_storage[edge]["left"]
                )
                right_partial_likelihood = pd.DataFrame(
                    partial_likelihood_per_site_storage[edge]["right"]
                )

                print(edge)
                print("Left: \n", left_partial_likelihood)
                print("Right: \n", right_partial_likelihood)


#    print("Partial Likelihoods: ", partial_likelihood_per_site_storage)
