import numpy as np
import pandas as pd

def calculate_test_statistic(
        multiplicity,
        array_eigenvectors,
        posterior_probabilities_left_subtree,
        posterior_probabilities_right_subtree,
        dimension,
        alpha = 0.05,
 ):      
    number_sites = len(posterior_probabilities_left_subtree["Site"].unique())
    number_nodes_1 = len(posterior_probabilities_left_subtree["Node"].unique())
    number_nodes_2 = len(posterior_probabilities_right_subtree["Node"].unique())

        if multiplicity == 1:  # if D=1
            v1 = array_eigenvectors[0]

            a = (
                []
            )  # vector to store all products v1*rootsitesposteriorprobabilitiescladeA
            b = (
                []
            )  # vector to store all products v1*rootsitesposteriorprobabilitiescladeB

            for k in range(
                number_sites * (number_nodes_1 - 1), number_sites * number_nodes_1
            ):
                a.append(
                    v1 @ np.asarray(posterior_probabilities_left_subtree.iloc[k, 3:7])
                )

            for k in range(
                number_sites * (number_nodes_2 - 1), number_sites * number_nodes_2
            ):
                b.append(
                    v1 @ np.asarray(posterior_probabilities_right_subtree.iloc[k, 3:7])
                )

            delta = np.asarray(a) @ np.asarray(b) / number_sites
            # computing the dominant sample coherence

            if i < len(internal_nodes):
                M_a = np.asarray(a) @ np.asarray(a) / number_sites 
                M_a = min(1, M_a)
            else:  # if clade A is a single leaf
                M_a = 1

            M_b = np.asarray(b) @ np.asarray(b) / number_sites
            M_b = min(1, M_b)
            variance = M_a * M_b / np.sqrt(number_sites)
            c_s = z_alpha * np.sqrt(variance)  # computing the saturation coherence

            c_sTwoSequence = z_alpha / np.sqrt(
                number_sites
            )  # computing the saturation coherence between two sequences

            p_value = st.norm.sf(abs(delta / np.sqrt(variance)))
        else:
            c_sTwoSequence = (
                multiplicity * z_alpha / np.sqrt(number_sites)
            )  # computing the saturation coherence between two sequences
            delta = 0

            for j in range(multiplicity):
                a = []
                b = []

                v1 = array_eigenvectors[j]

                for k in range(
                    number_sites * (number_nodes_1 - 1), number_sites * number_nodes_1
                ):
                    a.append(
                        v1
                        @ np.asarray(posterior_probabilities_left_subtree.iloc[k, 3:7])
                    )

                for k in range(
                    number_sites * (number_nodes_2 - 1), number_sites * number_nodes_2
                ):
                    b.append(
                        v1
                        @ np.asarray(posterior_probabilities_right_subtree.iloc[k, 3:7])
                    )

                delta += np.asarray(a) @ np.asarray(b)

            delta = delta / number_sites

            variance = 0

            for j in range(multiplicity):
                for k in range(multiplicity):
                    a = []
                    b = []

                    v_j = array_eigenvectors[j]
                    v_k = array_eigenvectors[k]

                    for l in range(
                        number_sites * (number_nodes_1 - 1),
                        number_sites * number_nodes_1,
                    ):
                        a.append(
                            v_j
                            @ np.asarray(
                                posterior_probabilities_left_subtree.iloc[l, 3:7]
                            )
                        )

                    for l in range(
                        number_sites * (number_nodes_2 - 1),
                        number_sites * number_nodes_2,
                    ):
                        b.append(
                            v_k
                            @ np.asarray(
                                posterior_probabilities_right_subtree.iloc[l, 3:7]
                            )
                        )

                    # variance = np.asarray(a)@np.asarray(b)
                    variance += max(
                        np.asarray(a) @ np.asarray(b),
                        np.asarray(a) @ np.asarray(b),
                        np.asarray(a) @ np.asarray(b),
                        np.asarray(a) @ np.asarray(b),
                    )

            variance = variance / (number_sites * number_sites)

            if variance < 0:
                print(
                    "VARIANCE ESTIMATION IS NEGATIVE - CONSIDER INCREASING THE NUMBER OF STANDARD DEVIATIONS (number_standard_deviations) (CONFIDENCE INTERVAL)"
                )
                c_s = 999999999
                p_value = -1
            else:
                c_s = z_alpha * np.sqrt(variance)
                p_value = st.norm.sf(abs(delta / np.sqrt(variance)))

        if c_s > delta:
            result_test = "Saturated"
        else:
            result_test = "Informative"

        if c_sTwoSequence > delta:
            result_test_tip2tip = "SatuT2T"
        else:
            result_test_tip2tip = "InfoT2T"
return delta, c_s, p_value, 