# -*- coding: utf-8 -*-
from Bio.Align import MultipleSeqAlignment

import sys
import argparse
import logging
import numpy as np
from pathlib import Path
from ete3 import Tree
from typing import Dict, List

from satute.decomposition import spectral_decomposition
from satute.rate_matrix import RateMatrix
from satute.file_handler import FileHandler, IqTreeHandler
from satute.trees import rename_internal_nodes_pre_order
from satute.arguments import ARGUMENT_LIST
from satute.sequences import check_if_tree_has_same_taxa_as_msa
from satute.logging import format_matrix, format_array

from satute.rate_analysis import (
    multiple_rate_analysis,
    single_rate_analysis_collapsed_tree,
)

from satute.ostream import (
    write_results_for_category_rates,
    write_alignment_and_indices,
    write_posterior_probabilities_single_rate,
    write_posterior_probabilities_for_rates,
    write_components,
)

from satute.categories import (
    read_alignment_file,
    split_msa_into_rate_categories_in_place,
    build_categories_by_sub_tables,
)

from satute.repository import (
    parse_substitution_model,
    parse_rate_from_cli_input,
    parse_file_to_data_frame,
    SubstitutionModel,
    IqTreeParser,
)

class Satute:
    """Class representing Satute command-line tool for wrapping up functions of IQ-TREE."""

    def __init__(
        self,
        iqtree=None,
        logger=None,
    ):
        # IQ-TREE related attributes
        self.iqtree = iqtree
        self.iqtree_tree_file = None
        self.iqtree_handler = None

        # Directories and paths
        self.input_dir = None
        self.site_probabilities_file = None
        self.active_directory = None

        self.alpha = 0.05

        # Miscellaneous attributes
        self.output_prefix = None
        self.input_args = []
        self.number_rates = 1
        self.logger = logger

        self.iq_tree_arguments_dict = {}

    """
        BEGIN Input Arguments Validation functions   
    """

    def validate_satute_input_options(self):
        """
        Validates the combinations of input arguments for Satute analysis.

        This function ensures that:
        - The -dir option is not used with specific other options (like -msa, -tree).
        - The -msa option is provided if -dir is not used.
        - The -tree option, if used, must be accompanied by a -model option.
        - The -ufboot or -boot options are not used with a specific combination of options.
        - The chosen category (if provided) is within a valid range.

        Raises:
            ValueError: If an invalid combination of arguments is provided.
        """
        self.validate_dir_conflicts()
        self.validate_msa_presence()
        self.validate_tree_and_model()
        self.validate_boot_options()
        self.validate_category_range()

    def validate_dir_conflicts(self):
        if self.input_args.dir and any(
            getattr(self.input_args, opt)
            for opt in ["msa", "tree", "model", "ufboot", "boot", "add_iqtree_options"]
        ):
            raise ValueError(
                "Cannot run Satute with -dir option combined with -'add_iqtree_options', -msa, -tree, -model, -ufboot, or -boot."
            )

    def validate_msa_presence(self):
        if not self.input_args.dir and not self.input_args.msa:
            raise ValueError("MSA file is required when -dir is not used.")
        if self.input_args.dir and self.input_args.msa:
            raise ValueError("MSA and dir cannot be used together.")

    def validate_tree_and_model(self):
        if self.input_args.tree and not self.input_args.model:
            raise ValueError("Model must be specified when using a tree file.")

    def validate_boot_options(self):
        if self.input_args.tree and (self.input_args.ufboot or self.input_args.boot):
            raise ValueError("Cannot use -ufboot or -boot with -tree.")

    def validate_category_range(self):
        if self.input_args.model and self.input_args.category:
            self.number_rates = self.handle_number_rates()
            if not (1 <= self.input_args.category <= self.number_rates):
                raise ValueError("Chosen category of interest is out of range.")

    def validate_and_set_rate_category(self, input_category: int, number_rates: int):
        """
        Validates the input category against the number of rates and sets the rate category.

        Args:
        - input_category (int): The category input from the user.
        - number_rates (int): The number of rates from the substitution model.

        Returns:
        - str: The validated and set rate category.

        Raises:
        - ValueError: If the input category is out of the valid range.
        """
        if not 1 <= input_category <= number_rates:
            self.logger.error("Chosen category of interest is out of range.")
            raise ValueError("Chosen category of interest is out of range.")
        return str(input_category)

    """END Input Arguments Validation functions """

    """ BEGIN Logging Function"""

    def setup_logging_configuration(self):
        """
        Initializes the logging system for the application.

        Sets up two handlers:
        1. A file handler that always logs at the DEBUG level.
        2. A stream (console) handler that logs at the DEBUG level if verbose is true; otherwise, it logs at the WARNING level.

        The log file is named using the MSA file name, alpha value, and an output suffix.
        """
        # Logger level is set to DEBUG to capture all logs for the file handler
        self.logger.setLevel(logging.DEBUG)
        # File Handler - always active at DEBUG level
        log_file = self.construct_log_file_name(Path(self.file_handler.find_msa_file()))

        file_handler = logging.FileHandler(log_file)
        file_format = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
        file_handler.setFormatter(file_format)
        file_handler.setLevel(logging.DEBUG)  # Always log everything in file
        self.logger.addHandler(file_handler)

        # Stream Handler - level depends on verbose flag
        stream_level = logging.DEBUG if self.input_args.verbose else logging.WARNING
        stream_handler = logging.StreamHandler()
        stream_handler.setFormatter(file_format)
        stream_handler.setLevel(stream_level)  # Set level based on verbose flag
        self.logger.addHandler(stream_handler)

    def log_substitution_model_info(
        self,
        substitution_model: SubstitutionModel,
        multiplicity: int,
        eigenvectors: List[np.array],
        eigenvalue: float,
    ):
        """
        Logs information about substitution model and its spectral decomposition

        Args:
            substitution_model: The substitution model used in the analysis.
            multiplicity: Multiplicity value from spectral decomposition.
            eigenvectors: eigenvector corresponding to eigenvalue from spectral decomposition
            eigenvalue:  dominant non-zero eigenvalue from spectral decomposition
        """
        # Formatting the rate matrix for logging
        rate_matrix_str = format_matrix(substitution_model.rate_matrix, precision=4)
        # Formatting the state frequencies for logging
        state_frequencies_str = format_array(
            np.array(list(substitution_model.state_frequencies)), precision=4
        )

        # Logging the formatted rate matrix and state frequencies
        self.logger.info(
            f"Substitution Model:\n\n"
            f"Model: {self.input_args.model}\n"
            f"Rate Matrix Q:\n{rate_matrix_str}\n"
            f"State Frequencies:\n{state_frequencies_str}\n"
        )

        eigenvector_str = ""
        for eigenvector in eigenvectors:
            eigenvector_str += f"\n{format_array(list(eigenvector))}"

        self.logger.info(
            f"Spectral Decomposition:\n\n"
            f"Eigenvalue: {eigenvalue}\n"
            f"Multiplicity: {multiplicity}\n"
            f"Eigenvectors: {eigenvector_str}\n"
        )

    def log_iqtree_run_and_satute_info(
        self,
        substitution_model: SubstitutionModel,
        rate_category: str,
        msa_file: Path,
        multiplicity: int,
    ):
        """
        Logs information about the initial IQ-TREE run and tests being performed.

        Args:
            iq_arguments_dict (dict): Dictionary containing IQ-TREE argument configurations.
            substitution_model: The substitution model used in the analysis.
            rate_category: The category of rates being considered.
            msa_file (Path): Path to the MSA file being used.
            multiplicity: Multiplicity value from spectral decomposition.
        """
        self.logger.info(
            f"""
            Running tests and initial IQ-Tree with configurations:
            Model: {self.input_args.model}
            Alpha: {self.input_args.alpha}
            Running Saturation Test on file: {msa_file.resolve()}
            Number of rate categories: {substitution_model.number_rates}
            Considered rate category: {rate_category}
            Multiplicity: {multiplicity}
            Run test for saturation for each branch and category with {substitution_model.number_rates} rate categories
            Results will be written to the directory: {self.active_directory.name}
            """
        )

    def construct_log_file_name(self, msa_file: Path):
        log_file = f"{msa_file.resolve()}_{self.input_args.alpha}.satute.log"
        if self.input_args.output_suffix:
            log_file = f"{msa_file.resolve()}_{self.input_args.alpha}_{self.input_args.output_suffix}.satute.log"
        return log_file

    """ END Logging Function"""

    def parse_input(self, args=None):
        """
        Parse command-line arguments using the argparse module. It dynamically
        adds arguments to the parser based on a predefined list of argument
        configurations, ensuring flexibility and ease of updates.
        """
        parser = argparse.ArgumentParser(description="Satute")
        for argument in ARGUMENT_LIST:
            # Unpack the dictionary directly without modifying the original list
            parser.add_argument(
                argument["flag"], **{k: v for k, v in argument.items() if k != "flag"}
            )
        self.input_args = parser.parse_args(args)

    def run_iqtree_workflow(self, arguments_dict: Dict[str, List]):
        extra_arguments = []

        if self.input_args.add_iqtree_options:
            extra_arguments.append(self.input_args.add_iqtree_options)

        if arguments_dict["option"] == "dir":
            self.logger.info(
                "IQ-TREE will not to be needed the analysis will be done on the already existing iqtree files."
            )
        else:
            # For the other options IQ-Tree is necessary. Therefore, test if IQ-TREE exists
            self.iqtree_handler.check_iqtree_path(self.input_args.iqtree)

        if arguments_dict["option"] == "dir + site_probabilities":
            self.logger.info("Running Satute with site probabilities")
            self.logger.info(
                "IQ-TREE will be needed for the site probabilities for the corresponding rate categories."
            )

            extra_arguments = extra_arguments + [
                "-m",
                self.input_args.model,
                "--redo",
                # "-blfix",
                "-wspr",
                "--quiet",
                "--keep-ident",
            ]

            self.iqtree_handler.run_iqtree_with_arguments(
                arguments_dict["arguments"], extra_arguments
            )

        if arguments_dict["option"] == "msa":
            self.logger.info("Running IQ-TREE with constructed arguments")

            self.logger.info(
                "If no model is specified in input arguments, best-fit model will be extracted from log file."
            )

            self.logger.warning(
                "Please consider for the analysis that IQ-Tree will be running with default options."
            )

            self.logger.warning(
                "If specific options are required for the analysis, please run IQ-Tree separately."
            )

            modelfinder_arguments = [
                "-m MF",
                "--quiet",
            ]

            self.iqtree_handler.run_iqtree_with_arguments(
                arguments=arguments_dict["arguments"],
                extra_arguments=modelfinder_arguments,
            ),

            # Update model in input arguments and re-construct arguments

            iq_tree_parser = IqTreeParser(f"{arguments_dict['msa_file']}.iqtree")

            self.input_args.model = iq_tree_parser.parse_substitution_model()
            self.handle_number_rates()

            # Validate and append ufboot and boot parameters to extra_arguments
            bb_arguments = self.iqtree_handler.validate_and_append_boot_arguments(
                self.input_args.ufboot, self.input_args.boot
            )

            extra_arguments = bb_arguments + [
                "-m",
                self.input_args.model,
                "--redo",
                "--keep-ident",
                "--quiet",
            ]

            if self.number_rates > 1:
                extra_arguments = extra_arguments + ["-wspr"]

            self.iqtree_handler.run_iqtree_with_arguments(
                arguments=arguments_dict["arguments"], extra_arguments=extra_arguments
            )

        if arguments_dict["option"] == "msa + model":
            self.logger.info("Running IQ-TREE with constructed arguments")
            self.logger.warning(
                "Please consider for the analysis that IQ-Tree will be running with default options."
            )
            self.logger.warning(
                "If specific options are required for the analysis, please run IQ-Tree separately."
            )

            self.handle_number_rates()

            bb_arguments = self.iqtree_handler.validate_and_append_boot_arguments(
                self.input_args.ufboot, self.input_args.boot
            )

            extra_arguments = bb_arguments + [
                "-m",
                self.input_args.model,
                "--quiet",
                "--keep-ident",
                "--redo",
            ]

            if isinstance(self.number_rates, int):
                # Add the '-wspr' option if number_rates > 1
                if self.number_rates > 1:
                    extra_arguments.append("-wspr")
            if isinstance(self.number_rates, str):
                if self.number_rates == "AMBIGUOUS":
                    extra_arguments.append("-wspr")

            self.iqtree_handler.run_iqtree_with_arguments(
                arguments=arguments_dict["arguments"], extra_arguments=extra_arguments
            )

        if arguments_dict["option"] == "msa + tree + model":
            self.logger.info("Running IQ-TREE with constructed arguments")
            self.logger.warning(
                "Please consider for the analysis that IQ-Tree will be running with default options."
            )
            self.logger.warning(
                "If specific options are required for the analysis, please run IQ-Tree separately."
            )

            self.handle_number_rates()

            extra_arguments = extra_arguments + [
                "-m",
                self.input_args.model,
                "--quiet",
                # "-blfix",
                "--keep-ident",
            ]

            if isinstance(self.number_rates, int):
                # Add the '-wspr' option if number_rates > 1
                if self.number_rates > 1:
                    extra_arguments.append("-wspr")
            if isinstance(self.number_rates, str):
                if self.number_rates == "AMBIGUOUS":
                    extra_arguments.append("-wspr")

            # Call IQ-TREE with the constructed arguments
            self.iqtree_handler.run_iqtree_with_arguments(
                arguments_dict["arguments"], extra_arguments
            )

        self.logger.info("Used IQ-TREE options:")
        self.logger.info(" ".join(arguments_dict["arguments"]))
        self.logger.info(" ".join(extra_arguments))

    def extract_tree(self):
        """
        Extracts and modifies a tree based on the Newick representation from either a specified file or a directory.

        This method first checks if a specific tree file is provided via input arguments. If so, it extracts
        the Newick representation from that file. Otherwise, it attempts to find and extract the Newick string
        from a tree file within the specified directory. After obtaining the Newick string, the method modifies
        the tree to ensure unique naming of the nodes based on a preorder traversal strategy.

        Returns:
            Tree: A modified ETE Tree object based on the Newick representation.

        Raises:
            ValueError: If no suitable tree file is found in the specified directory or provided via arguments.
        """
        try:
            if self.input_args.tree:
                newick_string = self.file_handler.get_newick_string(
                    str(self.input_args.tree.resolve())
                )
                self.logger.info(f"Using the user-defined tree: {self.input_args.tree}")
            else:
                newick_string = self.file_handler.get_newick_string_from_args()
        except (FileNotFoundError, ValueError) as e:
            self.logger.error(str(e))
            raise ValueError(f"Error extracting tree: {e}")
        return rename_internal_nodes_pre_order(Tree(newick_string, format=1))

    def run(self):
        """
        Main entry point for running the Satute command-line tool.
        """
        # Parsing and checking input arguments and constructing IQ-TREE command-line arguments
        # ======== Arguments =================
        msa_file = self.input_args.msa
        # ======== Tree File Handling ========
        newick_string = self.file_handler.get_newick_string_from_iq_tree_file(
            msa_file.resolve()
        )
        test_tree = rename_internal_nodes_pre_order(Tree(newick_string, format=1))
        # ======== Model parameter ===========
        ## Get dictionary for stationary distribution and diagonal matrix of the stationary distribution
        iq_tree_file_path = f"{msa_file.resolve()}.iqtree"
        satute_iq_tree_parser = IqTreeParser(iq_tree_file_path)
        substitution_model = satute_iq_tree_parser.load_substitution_model()
        ## Convert representation of rate_matrix
        RATE_MATRIX = RateMatrix(substitution_model.rate_matrix)

        # Calculation of the spectral decomposition of the rate matrix
        (
            array_left_eigenvectors,
            array_right_eigenvectors,
            multiplicity,
            eigenvalue,
        ) = spectral_decomposition(
            substitution_model.rate_matrix, substitution_model.phi_matrix
        )

        # Get number of rate categories in case of a +G or +R model
        # Consider a specific rate category

        rate_category = "all"
        if self.input_args.category:
            rate_category = self.validate_and_set_rate_category(
                self.input_args.category, substitution_model.number_rates
            )

        if self.number_rates == "AMBIGUOUS":
            self.number_rates = substitution_model.number_rates

        # ======== Multiple Sequence Alignment
        alignment = read_alignment_file(msa_file.resolve())
        # ========  Test for Branch Saturation =========

        self.log_substitution_model_info(
            substitution_model,
            multiplicity,
            array_right_eigenvectors,
            eigenvalue,
        )

        self.log_iqtree_run_and_satute_info(
            substitution_model,
            rate_category,
            msa_file,
            multiplicity,
        )

        check_if_tree_has_same_taxa_as_msa(sequence_alignment=alignment, tree=test_tree)

        if substitution_model.number_rates == 1:
            self.run_single_rate_analysis(
                test_tree,
                alignment,
                RATE_MATRIX,
                substitution_model.state_frequencies,
                array_right_eigenvectors,
                multiplicity,
                msa_file,
                self.input_args.alpha,
                self.input_args.edge,
            )
        else:
            self.run_multiple_rate_analysis(
                test_tree,
                substitution_model.category_rates,
                RATE_MATRIX,
                substitution_model.state_frequencies,
                array_right_eigenvectors,
                multiplicity,
                alignment,
                f"{msa_file.resolve()}.siteprob",
                rate_category,
                msa_file,
                self.input_args.alpha,
                self.input_args.edge,
            )

    def run_single_rate_analysis(
        self,
        test_tree: Tree,
        alignment: MultipleSeqAlignment,
        rate_matrix: RateMatrix,
        state_frequencies: List[float],
        array_right_eigenvectors: List[np.array],
        multiplicity: int,
        msa_file: Path,
        alpha: float,
        focused_edge: str,
    ):
        results = single_rate_analysis_collapsed_tree(
            test_tree,
            alignment,
            rate_matrix,
            state_frequencies,
            array_right_eigenvectors,
            multiplicity,
            alpha,
            focused_edge,
        )

        singe_rate_indices = [
            i for i in range(1, alignment.get_alignment_length() + 1, 1)
        ]

        single_rate_category = {"single_rate": singe_rate_indices}

        write_results_for_category_rates(
            results,
            test_tree,
            self.input_args.output_suffix,
            msa_file,
            self.input_args.alpha,
            self.input_args.edge,
            single_rate_category,
            self.logger,
        )

        if self.input_args.asr:
            self.logger.info(
                f"Writing ancestral sequences to {self.active_directory.name}"
            )

            write_posterior_probabilities_single_rate(
                results,
                state_frequencies,
                msa_file,
                self.input_args.output_suffix,
                self.input_args.alpha,
                self.input_args.edge,
                alignment.get_alignment_length(),
            )

        self.logger.info(
            f"Rate single_rate has {alignment.get_alignment_length()} Sites"
        )

        write_components(
            results["single_rate"]["components"],
            msa_file,
            self.input_args.output_suffix,
            "single_rate",
            self.input_args.alpha,
            self.input_args.edge,
            singe_rate_indices,
        )

    def run_multiple_rate_analysis(
        self,
        test_tree: Tree,
        category_rates: Dict[str, Dict],
        rate_matrix: RateMatrix,
        state_frequencies: List[float],
        array_right_eigenvectors: List[np.array],
        multiplicity: int,
        alignment: MultipleSeqAlignment,
        site_probability_file: str,
        rate_category: str,
        msa_file: Path,
        alpha: float = 0.05,
        edge: str = None,
    ):
        """
        Run the multiple rate analysis.
        """
        site_probability = parse_file_to_data_frame(site_probability_file)

        per_rate_category_alignment = split_msa_into_rate_categories_in_place(
            site_probability, alignment, rate_category
        )

        categorized_sites = build_categories_by_sub_tables(site_probability)

        results = multiple_rate_analysis(
            test_tree,
            category_rates,
            rate_matrix,
            state_frequencies,
            array_right_eigenvectors,
            multiplicity,
            per_rate_category_alignment,
            alpha,
            edge,
        )

        write_results_for_category_rates(
            results,
            test_tree,
            self.input_args.output_suffix,
            msa_file,
            alpha,
            edge,
            categorized_sites,
            self.logger,
        )

        for rate, alignment in categorized_sites.items():
            if len(alignment) == 0:
                self.logger.warning(f"Will be skipping Rate category {rate}")

        self.logger.info(
            f"Writing alignment and indices to {self.active_directory.name}"
        )

        if self.input_args.asr:

            self.logger.info(
                f"Writing ancestral sequences to {self.active_directory.name}"
            )

            write_posterior_probabilities_for_rates(
                results,
                state_frequencies,
                msa_file,
                self.input_args.output_suffix,
                self.input_args.alpha,
                self.input_args.edge,
                categorized_sites,
            )

        if self.input_args.category_assignment:
            write_alignment_and_indices(
                per_rate_category_alignment, categorized_sites, msa_file
            )

        for rate, results_set in results.items():
            self.logger.info(f"Rate {rate} has {len(categorized_sites[rate])} Sites")

            write_components(
                results_set["components"],
                msa_file,
                self.input_args.output_suffix,
                rate,
                self.input_args.alpha,
                self.input_args.edge,
                categorized_sites[rate],
            )

    def handle_number_rates(self):
        self.number_rates = 1
        if self.input_args.model:
            self.number_rates = parse_rate_from_cli_input(self.input_args.model)
        return self.number_rates

    def initialize_handlers(self):
        """
        Initialize the FileHandler and IqTreeHandler.
        """
        self.file_handler = FileHandler(self.active_directory)
        self.iqtree_handler = IqTreeHandler(self.input_args.iqtree)

    def initialize_active_directory(self):
        if self.input_args.msa:
            self.active_directory = self.input_args.msa.parent
        elif self.input_args.dir:
            self.active_directory = self.input_args.dir

    def get_dir_argument_options(self, msa_file: Path, tree_file: str):

        argument_option = {
            "option": "dir",
            "msa_file": Path(msa_file),
            "arguments": ["-s", str(Path(msa_file).resolve())],
        }

        if tree_file:
            argument_option["arguments"].extend(["-te", str(tree_file)])
        return argument_option

    """BEGIN Input Argument Construction"""

    def construct_IQ_TREE_arguments(self):
        """
        Validate and process input arguments.

        Raises:
            InvalidDirectoryError: If the input directory does not exist.
            NoAlignmentFileError: If no multiple sequence alignment file is found.

        Returns:
            A dictionary with keys 'option', 'msa_file', and 'arguments' that represents the argument options for the process.
        """
        # Define the acceptable file types for sequence alignments and trees
        argument_option = {}
        if self.input_args.dir:
            self.active_directory = self.input_args.dir
            self.input_args.msa = Path(self.file_handler.find_msa_file())
            tree_file = self.file_handler.find_tree_file()
            self.iqtree_tree_file = self.file_handler.find_iqtree_file()
            substitution_model = parse_substitution_model(self.iqtree_tree_file)
            self.input_args.model = substitution_model
            self.handle_number_rates()

            # Check for site probabilities file
            if self.number_rates > 1:
                self.site_probabilities_file = self.file_handler.find_file_by_suffix(
                    {".siteprob"}
                )
                if not self.site_probabilities_file:
                    self.logger.error(
                        "For the inference of tree a model with rate categories has been specified. But no site probabilities file found in directory. Please rerun iqtree with the option -wspr"
                    )
                    exit(1)

            return self.get_dir_argument_options(self.input_args.msa, tree_file)
        else:
            argument_option = {}
            if self.input_args.msa:
                self.construct_arguments_for_msa(argument_option)
            if self.input_args.tree:
                self.construct_argument_for_tree(argument_option)
            # If a model was specified in the input arguments, add it to the argument options
            if self.input_args.model:
                self.construct_argument_for_model(argument_option)
        # Return the constructed argument options
        return argument_option

    def construct_arguments_for_msa(self, argument_option: Dict):
        """
        Constructs the arguments required for multiple sequence alignment (MSA).

        Args:
        argument_option (dict): The dictionary where MSA arguments are to be added.

        Returns:
        dict: Updated dictionary with MSA specific arguments.
        """
        # Specify the option type and MSA file path
        argument_option["option"] = "msa"
        argument_option["msa_file"] = self.input_args.msa
        # Add MSA specific command-line arguments
        argument_option["arguments"] = ["-s", str(self.input_args.msa.resolve())]
        return argument_option

    def construct_argument_for_tree(self, argument_option: Dict):
        """
        Appends tree-related arguments to the existing argument option.

        Args:
        argument_option (dict): The dictionary to which tree arguments are added.
        """
        # Specify that the option now includes tree data
        argument_option["option"] = "msa + tree"

        # Extend the existing arguments with tree specific command-line arguments
        argument_option["arguments"].extend(
            ["-te", str(self.input_args.tree.resolve())]
        )

    def construct_argument_for_model(self, argument_option: Dict):
        """
        Constructs the arguments required for the evolutionary model.

        Args:
        argument_option (dict): The dictionary where model arguments are to be added.
        """
        # Update the option to indicate inclusion of the model
        argument_option["option"] += " + model"

        # Add model specific command-line arguments
        argument_option["model_arguments"] = ["-m", self.input_args.model]

        # Add extra arguments if the model includes certain features
        if "+G" in self.input_args.model or "+R" in self.input_args.model:
            argument_option["model_arguments"].extend(["-wspr"])

    """END Input Argument Construction"""

"""
def main(args=None):
    # Instantiate the Satute class
    logger = logging.getLogger(__name__)
    satute = Satute(iqtree="iqtree", logger=logger)
    # Parse and validate input arguments
    satute.parse_input(args)
    satute.validate_satute_input_options()
    # Initialize file handler and logger
    satute.initialize_active_directory()
    satute.initialize_handlers()
    satute.setup_logging_configuration()
    # IQ-Tree run if necessary
    satute.iq_arguments_dict = satute.construct_IQ_TREE_arguments()
    satute.run_iqtree_workflow(satute.iq_arguments_dict)
    # Run the tool
    satute.run()


def cli():
    main()
"""


def main(args=None):
    # Instantiate the Satute class
    logger = logging.getLogger(__name__)
    satute = Satute(iqtree="iqtree", logger=logger)
    
    try:
        # Parse and validate input arguments
        satute.parse_input(args)
        satute.validate_satute_input_options()
        
        # Initialize file handler and logger
        satute.initialize_active_directory()
        satute.initialize_handlers()
        satute.setup_logging_configuration()
        
        # IQ-Tree run if necessary
        satute.iq_arguments_dict = satute.construct_IQ_TREE_arguments()
        satute.run_iqtree_workflow(satute.iq_arguments_dict)
        
        # Run the tool
        satute.run()
    except Exception as e:
        logger.error(f"An error occurred: {e}")
        sys.exit(1)
    
    # Optionally return an exit code or nothing
    return 0

def cli():
    # args = sys.argv[1:]
    exit_code = main()
    sys.exit(exit_code)

