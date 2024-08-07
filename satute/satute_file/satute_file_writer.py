import numpy as np
import pandas as pd
from ete3 import Tree
from pathlib import Path
from typing import List, Union, Dict

from satute.ostream import construct_file_name
from satute.models.substitution_model import SubstitutionModel
from satute.satute_file.file_writer import FileWriter


class SatuteFileWriter(FileWriter):
    def write_substitution_model_info(
        self,
        substitution_model: SubstitutionModel,
        multiplicity: int,
        eigenvectors: List[np.array],
        eigenvalue: float,
        option: str,
        category_rates: Dict[str, Dict[str, Union[int, float]]],
        categorized_sites: Dict[str, List[int]],
        alignment_length: int,
    ) -> None:
        """
        Logs information about substitution model and its spectral decomposition

        Args:
            substitution_model: The substitution model used in the analysis.
            multiplicity: Multiplicity value from spectral decomposition.
            eigenvectors: eigenvector corresponding to eigenvalue from spectral decomposition
            eigenvalue:  dominant non-zero eigenvalue from spectral decomposition
        """

        self.write_to_file("""\n\nSUBSTITUTION PROCESS\n\n""")
        # Formatting the rate matrix for loggings
        rate_matrix_str: str = format_matrix(
            substitution_model.rate_matrix, precision=4
        )
        # Formatting the state frequencies for logging
        state_distribution_str: str = format_array(
            np.array(list(substitution_model.state_frequencies)), precision=4
        )

        substitution_model_text: str = (
            f"Substitution Model: {substitution_model.model}\n\n"
        )

        model_text: str = (
            "Best fitted model: " + substitution_model_text
            if option == "msa"
            else "User chosen model: " + substitution_model_text
        )

        if substitution_model.gamma_shape:
            self.write_to_file(f"Gamma shape alpha: {substitution_model.gamma_shape}")

        # Logging the formatted rate matrix and state frequencies
        self.write_to_file(
            f"\n{model_text}"
            f"Rate Matrix Q:\n{rate_matrix_str}\n\n"
            f"Stationary Distribution:\n{state_distribution_str}\n"
        )

        self.write_rate_categories(
            category_rates=category_rates,
            categorized_sites=categorized_sites,
            alignment_length=alignment_length,
        )

        eigenvector_str = ""
        for eigenvector in eigenvectors:
            eigenvector_str += f"\n{format_array(list(eigenvector))}"

        self.write_to_file(
            f"\n\n"
            f"Spectral Decomposition:\n\n"
            f"Eigenvalue: {eigenvalue}\n\n"
            f"Multiplicity: {multiplicity}\n\n"
            f"Eigenvectors: {eigenvector_str}\n"
            f"\n\n"
        )

    def write_rate_categories(
        self,
        category_rates: Dict[str, Dict],
        categorized_sites: Dict[str, List],
        alignment_length: float,
    ):
        if category_rates:
            for rate, value in category_rates.items():
                value["empirical"] = len(categorized_sites[rate]) / alignment_length
            # Convert the category rates dictionary to a pandas DataFrame for table formatting
            df = pd.DataFrame.from_dict(category_rates, orient="index")
            df.reset_index(drop=True, inplace=True)
            df.columns = [
                "Category",
                "Relative_rate",
                "Proportion",
                "Empirical Proportion (Number of Sites)",
            ]
            self.write_to_file("\n" + df.to_string(index=False) + "\n")

    def write_which_tested_tree(self, tree: Tree, option: str) -> None:
        if "tree" in option:
            self.write_to_file(
                f"User defined Tree: {tree.write(format=1, format_root_node=True)}"
            )
        else:
            self.write_to_file(
                f"IQ-Tree inferred Tree: {tree.write(format=1, format_root_node=True)}"
            )

    def write_intro(self):
        self.write_to_file("SatuTe Version 0.0.1 \n\n")
        self.write_to_file(
            """To cite SatuTe please use: When the past fades: Detecting phylogenetic signal with SatuTe: Cassius Manuel, Christiane Elgert, Enes Sakalli, Heiko A. Schmidt1, Carme Viñas and Arndt von Haeseler\n\n"""
        )

    def write_alignment_info(self, msa_file, option):
        if "dir" in option:
            self.write_to_file(
                f"\nUsed Alignment File from directory {msa_file.resolve()}\n"
            )
        else:
            self.write_to_file(f"\nUsed Alignment File {msa_file.resolve()}\n")

    def write_results_to_csv(self, msa_file, results, input_args):
        self.write_to_file(
            "\nThe satute.csv file provides a comprehensive overview of the saturation test results for specific branches or all branches.\n\n"
        )

        for rate, results_set in results.items():
            replaced_rate_category = rate.replace("p", "c")
            self.write_to_file(
                f"Writing results for {replaced_rate_category} to CSV File: {construct_file_name(msa_file, input_args.output_suffix, replaced_rate_category, input_args.alpha, input_args.edge)}.csv\n"
            )

    def write_results_to_nexus(self, msa_file, results, input_args):
        self.write_to_file(
            f"""\nA satute.nex file consists of two blocks, each enclosed by BEGIN and END statements.\nThese blocks contain the taxon labels and the phylogenetic tree, with the most important test results integrated into the NEWICK string as metadata.\n"""
        )
        for rate, results_set in results.items():
            replaced_rate_category = rate.replace("p", "c")
            self.write_to_file(
                f"\nWriting results for {replaced_rate_category} to Nexus File: {construct_file_name(msa_file, input_args.output_suffix, replaced_rate_category, input_args.alpha, input_args.edge)}.nex"
            )

    def write_considered_rate_category(
        self, rate_category: int, substitution_model: SubstitutionModel
    ):
        considered_rate_category_text = (
            f"Test will be applied on the rate categories: {rate_category}"
            if substitution_model.number_rates > 1
            else ""
        )
        self.write_to_file(f"{considered_rate_category_text}")

    def write_significance_level(self, alpha: float):
        considered_significance_level_text = (
            f"\nTest will be applied with the significance level: {alpha}\n"
        )
        self.write_to_file(f"{considered_significance_level_text}")

    def write_considered_rate_category(
        self, rate_category: int, substitution_model: SubstitutionModel
    ):
        considered_rate_category_text = (
            f"Test will be applied on the rate categories: {rate_category}"
            if substitution_model.number_rates > 1
            else ""
        )
        self.write_to_file(f"{considered_rate_category_text}")

    def write_considered_branch(self, input_args):
        considered_edge = (
            f"Run test for saturation for the branch: {input_args.edge}"
            if input_args.edge
            else "Run test for saturation for each branch"
        )
        self.write_to_file(f"{considered_edge}")

    def write_input_source(self, iq_tree_file: str):
        self.write_to_file(f"\nTree and parameters are read from: {iq_tree_file}\n")

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
        self.write_intro()
        self.write_input_source(iq_tree_file=iq_tree_file)

        self.write_alignment_info(
            msa_file=msa_file,
            option=iqtree_arguments["option"],
        )

        self.write_considered_branch(input_args=input_args)

        self.write_which_tested_tree(
            tree=test_tree,
            option=iqtree_arguments["option"],
        )

        self.write_substitution_model_info(
            substitution_model=substitution_model,
            multiplicity=multiplicity,
            eigenvectors=array_right_eigenvectors,
            eigenvalue=eigenvalue,
            option=iqtree_arguments["option"],
            category_rates=substitution_model.category_rates,
            categorized_sites=categorized_sites,
            alignment_length=alignment_length,
        )

        self.write_significance_level(alpha=alpha)

        self.write_considered_rate_category(
            rate_category=rate_category,
            substitution_model=substitution_model,
        )

        self.write_to_file("\n\nCSV FILES\n\n")

        self.write_results_to_csv(
            msa_file=msa_file,
            results=results,
            input_args=input_args,
        )

        self.write_to_file("\nNEXUS FILES\n\n")

        self.write_results_to_nexus(
            msa_file=msa_file,
            results=results,
            input_args=input_args,
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
