import numpy as np
import pandas as pd
from ete3 import Tree
from pathlib import Path
from typing import List, Dict
from datetime import datetime

from satute.ostream import construct_file_name
from satute.models.substitution_model import SubstitutionModel
from satute.satute_file.file_writer import FileWriter
from satute.messages.messages import SATUTE_VERSION


class SatuteFileWriter(FileWriter):
    def write_substitution_model_info(
        self,
        substitution_model: SubstitutionModel,
        option: str,
    ) -> None:
        """
        Logs information about substitution model and its spectral decomposition

        Args:
            substitution_model: The substitution model used in the analysis.
            multiplicity: Multiplicity value from spectral decomposition.
            eigenvectors: eigenvector corresponding to eigenvalue from spectral decomposition
            eigenvalue:  dominant non-zero eigenvalue from spectral decomposition
        """
        # Formatting the rate matrix for loggings
        rate_matrix_str: str = format_matrix(
            substitution_model.rate_matrix, precision=4
        )
        # Formatting the state frequencies for logging
        state_distribution_str: str = format_array(
            np.array(list(substitution_model.state_frequencies)), precision=4
        )

        substitution_model_text: str = f"{substitution_model.model}\n\n"

        model_text: str = (
            "Best fitted substitution model:" + substitution_model_text
            if option == "msa"
            else "User chosen substitution model: " + substitution_model_text
        )

        # Logging the formatted rate matrix and state frequencies
        self.write_to_file(
            f"\n\n{model_text}"
            f"Rate Matrix Q:\n{rate_matrix_str}\n\n"
            f"Stationary Distribution:\n{state_distribution_str}\n\n"
        )

    def write_rate_categories(
        self,
        category_rates: Dict[str, Dict],
        categorized_sites: Dict[str, List],
        alignment_length: float,
        substitution_model: SubstitutionModel,
    ):
        if category_rates:
            if substitution_model.gamma_shape:
                self.write_to_file(
                    f"\nGamma shape parameter: {substitution_model.gamma_shape}\n"
                )

            for rate, value in category_rates.items():
                value["empirical"] = len(categorized_sites[rate]) / alignment_length
            # Convert the category rates dictionary to a pandas DataFrame for table formatting

            df = pd.DataFrame.from_dict(category_rates, orient="index")
            df.reset_index(drop=True, inplace=True)
            df.columns = [
                "Category",
                "Relative_rate",
                "Proportion",
                f"Empirical Proportion (alignment length {alignment_length})",
            ]

            self.write_to_file("\n" + df.to_string(index=False) + "\n\n")
        else:
            self.write_to_file("\n\nRate Category: None")
            self.write_to_file("\n\n")

    def write_which_tested_tree(self, tree: Tree, option: str) -> None:
        if "tree" in option:
            self.write_to_file(
                f"\nUser defined tree: {tree.write(format=1, format_root_node=True)}\n\n"
            )
        else:
            self.write_to_file(
                f"\nIQ-Tree inferred tree: {tree.write(format=1, format_root_node=True)}\n\n"
            )

    def write_intro(self):
        self.write_header("REFERENCE")
        self.write_to_file(SATUTE_VERSION)

    def write_alignment_info(self, msa_file: Path, option: str):
                
        if "dir" in option:
            self.write_to_file(
                f"\n\nUsed alignment file from directory: {msa_file.resolve()}\n"
            )
        else:
            self.write_to_file(f"\n\nUsed alignment file: {msa_file.resolve()}\n")

    def write_results_to_csv(self, msa_file, results, input_args):
        self.write_to_file(
            "\nThe satute.csv file provides a comprehensive overview of the saturation test results for a specific branch or all branches."
        )
        self.write_to_file(
            "\nContaining: branch, mean_coherence, standard_error_of_mean, z_score, p_value, z_alpha, decision_test,\nz_alpha_bonferroni_corrected,decision_bonferroni_corrected, branch_length, number_of_sites, rate_category\n\n"
        )

        for rate, results_set in results.items():
            replaced_rate_category = rate.replace("p", "c")
            self.write_to_file(
                f"Results for category {replaced_rate_category} to CSV File: {construct_file_name(msa_file, input_args.output_suffix, replaced_rate_category, input_args.alpha, input_args.edge)}.csv\n"
            )

    def write_results_to_nexus(self, msa_file, results, input_args):
        self.write_to_file(
            "\nThe file contains a block for the taxon labels and a block for the phylogenetic tree,\nwith the most important test results integrated into the NEWICK string as metadata."
        )
        self.write_to_file("\nContaining: z_score, p_value, decision_test, decision_bonferroni_corrected.\n\n")

        for rate, results_set in results.items():
            replaced_rate_category = rate.replace("p", "c")
            self.write_to_file(
                f"Results for category {replaced_rate_category} to Nexus File: {construct_file_name(msa_file, input_args.output_suffix, replaced_rate_category, input_args.alpha, input_args.edge)}.nex\n"
            )

    def write_components(self, msa_file, results, input_args):
        self.write_to_file(
            "\nThe components file provides the estimated (category) variance and the coherence coefficient for each site and branch in \nthe tree, enabling other analysis of the saturation status like sliding window analysis."
        )

        self.write_to_file(
            "\nContaining: branch, site, coherence, category_variance, rate_category.\n\n"
        )


        for rate, results_set in results.items():
            replaced_rate_category = rate.replace("p", "c")
            self.write_to_file(
                f"Results for category {replaced_rate_category} to CSV File: {construct_file_name(msa_file, input_args.output_suffix, replaced_rate_category, input_args.alpha, input_args.edge)}.components.csv\n"
            )

    def write_results_to_ancestral_states(self, msa_file, results, input_args):
        self.write_to_file(
            "The  file contains the posterior distributions of ancestral sequences for the left and right node of each branch in the tree.\n\n"
        )

        for rate, results_set in results.items():
            replaced_rate_category = rate.replace("p", "c")
            self.write_to_file(
                f"\nWrite Ancestral States for rate category {replaced_rate_category} to CSV: {construct_file_name(msa_file, input_args.output_suffix, replaced_rate_category, input_args.alpha, input_args.edge)}.asr.csv"
            )

    def write_results_to_rate_index(self, msa_file, results, input_args):
        self.write_to_file(
            "The file contains a multiple sequence alignment subdivided into categories which the sites are assigned to.\n"
        )

        for rate, results_set in results.items():
            replaced_rate_category = rate.replace("p", "c")
            self.write_to_file(
                f"\nWrite Ancestral States{replaced_rate_category} to CSV: {construct_file_name(msa_file, input_args.output_suffix, replaced_rate_category, input_args.alpha, input_args.edge)}.rateidx"
            )

    def write_considered_rate_category(
        self, rate_category: int, substitution_model: SubstitutionModel
    ):
        considered_rate_category_text = (
            f"\nTested rate categories: {rate_category}"
            if substitution_model.number_rates > 1
            else ""
        )
        self.write_to_file(f"{considered_rate_category_text}\n")

    def write_significance_level(self, alpha: float):
        considered_significance_level_text = f"\n\nSignificance level: {alpha}\n"
        self.write_to_file(f"{considered_significance_level_text}")

    def write_considered_branch(self, input_args):
        considered_edge = (
            f"Tested branch: {input_args.edge}\n"
            if input_args.edge
            else "Tested branches: all\n"
        )
        self.write_to_file(f"\n{considered_edge}")

    def write_input_source(self, iq_tree_file: str):
        self.write_to_file(
            f"\nTree and substitution model are read from: {iq_tree_file}\n"
        )

    def write_satute_file(
        self,
        msa_file: str,
        iq_tree_file: Path,
        test_tree: Tree,
        rate_category: int,
        substitution_model: SubstitutionModel,
        multiplicity: int,
        eigenvalue: float,
        array_right_eigenvectors: List[np.array],
        iqtree_arguments,
        input_args,
        categorized_sites,
        alignment_length,
        alpha,
        results,
    ):
        self.write_header("TIMESTAMP")

        # Get the current timestamp
        current_timestamp = datetime.now()

        # Format the timestamp
        formatted_timestamp = current_timestamp.strftime("%Y-%m-%d %H:%M:%S \n")

        # Print the formatted timestamp
        self.write_to_file(formatted_timestamp)

        self.write_intro()

        self.write_header("INPUT")

        self.write_alignment_info(
            msa_file=msa_file,
            option=iqtree_arguments["option"],
        )

        self.write_input_source(iq_tree_file=iq_tree_file)

        self.write_which_tested_tree(
            tree=test_tree,
            option=iqtree_arguments["option"],
        )

        self.write_header("SUBSTITUTION MODEL")

        self.write_substitution_model_info(
            substitution_model=substitution_model, option=iqtree_arguments["option"]
        )

        self.write_header("MODEL FOR RATE HETEROGENEITY")

        self.write_rate_categories(
            category_rates=substitution_model.category_rates,
            categorized_sites=categorized_sites,
            alignment_length=alignment_length,
            substitution_model=substitution_model,
        )

        self.write_spectral_decomposition(
            eigenvalue=eigenvalue,
            multiplicity=multiplicity,
            eigenvectors=array_right_eigenvectors,
        )

        self.write_header("SATURATION TEST")

        self.write_significance_level(alpha=alpha)

        self.write_considered_branch(input_args=input_args)

        self.write_considered_rate_category(
            rate_category=rate_category, substitution_model=substitution_model
        )

        self.write_header("OUTPUT")

        self.write_to_file("\n\nCSV FILES\n")

        self.write_results_to_csv(
            msa_file=msa_file,
            results=results,
            input_args=input_args,
        )

        self.write_to_file("\n\nCOHERENCE AND VARIANCE PER SITE FOR EACH CATEGORY\n")

        self.write_components(
            msa_file=msa_file,
            results=results,
            input_args=input_args,
        )

        self.write_to_file("\n\nNEXUS FILES\n")

        self.write_results_to_nexus(
            msa_file=msa_file,
            results=results,
            input_args=input_args,
        )

        if input_args.asr:
            self.write_to_file("\n\n\nPOSTERIOR DISTRIBUTIONS\n")

            self.write_results_to_ancestral_states(
                msa_file=msa_file,
                results=results,
                input_args=input_args,
            )

        if input_args.category_assignment:
            self.write_to_file("\n\nRATE CATEGORY ASSIGNMENTS\n")

            self.write_results_to_rate_index(
                msa_file=msa_file,
                results=results,
                input_args=input_args,
            )

    def write_spectral_decomposition(self, eigenvalue, multiplicity, eigenvectors):
        self.write_header("SPECTRAL DECOMPOSITION")

        eigenvector_str = ""
        for eigenvector in eigenvectors:
            eigenvector_str += f"\n{format_array(list(eigenvector))}"

        self.write_to_file(
            f"\n\n"
            f"Second Largest Eigenvalue: {eigenvalue}\n\n"
            f"Multiplicity: {multiplicity}\n\n"
            f"Eigenvectors: {eigenvector_str}\n\n"
        )


def format_matrix(matrix, precision: int = 4):
    """Format a matrix for pretty printing."""
    formatted_matrix = "\n".join(
        ["\t".join([f"{item:.{precision}f}" for item in row]) for row in matrix]
    )
    return formatted_matrix


def format_array(array, precision=4):
    """Format a 1D array for pretty printing."""
    formatted_array = "\t".join([f"{item:.{precision}f}" for item in array])
    return formatted_array
