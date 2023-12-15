#!/usr/bin/env python
import argparse
import logging
from pathlib import Path
from ete3 import Tree
from satute_util import spectral_decomposition
from rate_matrix import RateMatrix
from file_handler import FileHandler, IqTreeHandler
from satute_trees_and_subtrees import rename_internal_nodes_preorder
from satute_arguments import ARGUMENT_LIST
from satute_rate_analysis import (
    multiple_rate_analysis,
    single_rate_analysis_collapsed_tree,
)
from satute_ostream import (
    write_results_for_category_rates,
    write_results_for_single_rate,
    write_alignment_and_indices,
)
from satute_categories import (
    read_alignment_file,
    split_msa_into_rate_categories_in_place,
    build_categories_by_sub_tables,
)
from satute_repository import (
    parse_substitution_model,
    parse_rate_from_model,
    parse_file_to_data_frame,
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
            for opt in ["msa", "tree", "model", "ufboot", "boot"]
        ):
            raise ValueError(
                "Cannot run Satute with -dir option combined with -msa, -tree, -model, -ufboot, or -boot."
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

    def parse_input(self):
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
        self.input_args = parser.parse_args()

    def run_iqtree_workflow(self, arguments_dict):
        extra_arguments = []
        # if input.args.additional:
        #    extra_arguments.append(input.args.additional)
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

            extra_arguments = [
                "-m",
                self.input_args.model,
                "--redo",
                "-blfix",
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
            ]

            if self.number_rates > 1:
                extra_arguments = extra_arguments + ["-wspr"]

            self.iqtree_handler.run_iqtree_with_arguments(
                arguments=arguments_dict["arguments"], extra_arguments=extra_arguments
            )

        if arguments_dict["option"] == "msa + tree + model":
            logger.info("Running IQ-TREE with constructed arguments")
            logger.warning(
                "Please consider for the analysis that IQ-Tree will be running with default options."
            )
            logger.warning(
                "If specific options are required for the analysis, please run IQ-Tree separately."
            )

            self.handle_number_rates()

            extra_arguments = [
                "-m",
                self.input_args.model,
                "--quiet",
                "-blfix",
                "--keep-ident",
            ]

            # Add the '-wspr' option if number_rates > 1
            if self.number_rates > 1:
                extra_arguments.append("-wspr")

            # Call IQ-TREE with the constructed arguments
            self.iqtree_handler.run_iqtree_with_arguments(
                arguments_dict["arguments"], extra_arguments
            )

        logger.info("Used IQ-TREE options:")
        logger.info(" ".join(arguments_dict["arguments"]))
        logger.info(" ".join(extra_arguments))

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
        return rename_internal_nodes_preorder(Tree(newick_string, format=1))

    def initialise_logger(self):
        log_file = f"{self.input_args.msa.resolve()}_{self.input_args.alpha}.satute.log"
        file_handler = logging.FileHandler(log_file)
        # Add the FileHandler to the logger
        self.logger.addHandler(file_handler)

    def validate_and_set_rate_category(self, input_category, number_rates):
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
            logger.error("Chosen category of interest is out of range.")
            raise ValueError("Chosen category of interest is out of range.")
        return str(input_category)

    def run(self):
        """
        Main entry point for running the Satute command-line tool.
        """
        # Parsing and checking input arguments and constructing IQ-TREE command-line arguments
        # ======== Input =====================
        self.parse_input()
        # ======== Validation ================
        self.validate_satute_input_options()
        # ======== Arguments =================
        iq_arguments_dict = self.construct_iq_arguments()
        msa_file = iq_arguments_dict["msa_file"]
        # ======== Logger ====================
        self.initialise_logger()
        # ======== IQ-Tree ===================
        self.run_iqtree_workflow(iq_arguments_dict)
        # ======== Tree File Handling ========
        to_be_tested_tree = self.extract_tree()
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
        # ======== Multiple Sequence Alignment
        alignment = read_alignment_file(msa_file.resolve())
        # ========  Test for Branch Saturation =========
        self.logger.info(
            f"""
            Running tests and initial IQ-Tree with configurations:
            Mode {iq_arguments_dict['option']}
            Model: {self.input_args.model}
            Alpha: {self.input_args.alpha}
            Running Saturation Test on file: {msa_file.resolve()}
            Number of rate categories: {substitution_model.number_rates}
            Considered rate category: {rate_category}
            Options for Initial IQ-Tree run: {iq_arguments_dict['option']}
            Multiplicity: {multiplicity}
            Run test for saturation for each branch and category with {substitution_model.number_rates} rate categories
            Results will be written to the directory:{self.active_directory.name}
            """
        )

        if substitution_model.number_rates == 1:
            self.run_single_rate_analysis(
                to_be_tested_tree,
                alignment,
                RATE_MATRIX,
                substitution_model.state_frequencies.values(),
                array_left_eigenvectors,
                array_right_eigenvectors,
                multiplicity,
                msa_file,
                self.input_args.alpha,
                self.input_args.edge,
            )

        else:
            self.run_multiple_rate_analysis(
                to_be_tested_tree,
                substitution_model.category_rates,
                RATE_MATRIX,
                substitution_model.state_frequencies.values(),
                array_left_eigenvectors,
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
        to_be_tested_tree,
        alignment,
        rate_matrix,
        state_frequencies,
        array_left_eigenvectors,
        array_right_eigenvectors,
        multiplicity,
        msa_file,
        alpha,
        focused_edge,
    ):
        results = single_rate_analysis_collapsed_tree(
            to_be_tested_tree,
            alignment,
            rate_matrix,
            state_frequencies,
            array_left_eigenvectors,
            array_right_eigenvectors,
            multiplicity,
            alpha,
            focused_edge,
        )

        write_results_for_single_rate(
            results,
            msa_file,
            to_be_tested_tree,
            self.input_args.alpha,
            self.input_args.edge,
            logger,
        )

    def run_multiple_rate_analysis(
        self,
        to_be_tested_tree,
        category_rates,
        rate_matrix,
        state_frequencies,
        array_left_eigenvectors,
        array_right_eigenvectors,
        multiplicity,
        alignment,
        site_probability_file,
        rate_category,
        msa_file,
        alpha=0.05,
        edge=None,
    ):
        """
        Run the multiple rate analysis.
        """
        site_probability = parse_file_to_data_frame(site_probability_file)

        per_rate_category_alignment = split_msa_into_rate_categories_in_place(
            site_probability, alignment, rate_category
        )

        categorized_sites = build_categories_by_sub_tables(site_probability)

        write_alignment_and_indices(
            per_rate_category_alignment, categorized_sites, self.input_args
        )

        print(to_be_tested_tree)

        results = multiple_rate_analysis(
            to_be_tested_tree,
            category_rates,
            rate_matrix,
            state_frequencies,
            array_left_eigenvectors,
            array_right_eigenvectors,
            multiplicity,
            per_rate_category_alignment,
            alpha,
            edge,
        )
        write_results_for_category_rates(results, msa_file, alpha, edge, logger)

    def handle_number_rates(self):
        self.number_rates = 1
        if self.input_args.model:
            self.number_rates = parse_rate_from_model(self.input_args.model)
        return self.number_rates

    def initialize_handlers(self):
        """
        Initialize the FileHandler and IqTreeHandler.
        """
        self.file_handler = FileHandler(self.active_directory)
        self.iqtree_handler = IqTreeHandler(self.input_args.iqtree)

    def initialize_working_context(self):
        if self.input_args.msa:
            self.active_directory = self.input_args.msa.parent
        elif self.input_args.dir:
            self.active_directory = self.input_args.dir

    def construct_iq_arguments(self):
        """
        Validate and process input arguments.

        Raises:
            InvalidDirectoryError: If the input directory does not exist.
            NoAlignmentFileError: If no multiple sequence alignment file is found.

        Returns:
            A dictionary with keys 'option', 'msa_file', and 'arguments' that represents the argument options for the process.
        """
        # Define the acceptable file types for sequence alignments and trees
        msa_file_types = {".fasta", ".nex", ".phy", ".txt"}
        tree_file_types = {".treefile", ".nex", ".nwk"}
        self.initialize_working_context()
        self.initialize_handlers()
        argument_option = {}

        if self.input_args.dir:
            self.active_directory = self.input_args.dir

            # Check if the input directory exists

            # Find msa file in the directory
            self.input_args.msa_file = self.file_handler.find_file_by_suffix(
                msa_file_types
            )

            self.input_args.msa = Path(self.input_args.msa_file)

            argument_option = {
                "option": "dir",
                "msa_file": self.input_args.msa,
                "arguments": ["-s", str(self.input_args.msa.resolve())],
            }

            # Find the tree file in the directory
            tree_file = self.file_handler.find_file_by_suffix(tree_file_types)

            # Check if a tree file was found
            if tree_file:
                argument_option["arguments"].extend(["-te", str(tree_file)])
            else:
                raise FileNotFoundError("No tree file found in directory")

            # Find .iqtree file in the directory
            self.iqtree_tree_file = self.file_handler.find_file_by_suffix({".iqtree"})

            # Check if iqtree file was found
            if self.iqtree_tree_file:
                substitution_model = parse_substitution_model(self.iqtree_tree_file)
                self.input_args.model = substitution_model
            else:
                raise FileNotFoundError("No iqtree file found in directory")

            self.handle_number_rates()

            # How apply iqtree just to get the site probabilities file without --redo
            if self.number_rates > 1:
                self.site_probabilities_file = self.file_handler.find_file_by_suffix(
                    {".siteprob"}
                )
                if not self.site_probabilities_file:
                    argument_option["option"] += " + site_probabilities"

            return argument_option
        else:
            if self.input_args.msa:
                argument_option = {
                    "option": "msa",
                    "msa_file": self.input_args.msa,
                    "arguments": ["-s", str(self.input_args.msa.resolve())],
                }

            if self.input_args.tree:
                argument_option["option"] = "msa + tree"
                argument_option["arguments"].extend(
                    ["-te", str(self.input_args.tree.resolve())]
                )

            # If a model was specified in the input arguments, add it to the argument options
            if self.input_args.model:
                argument_option["option"] += " + model"
                argument_option["model_arguments"] = ["-m", self.input_args.model]
                # If the model includes a Gamma distribution, add the corresponding argument
                if "+G" in self.input_args.model or "+R" in self.input_args.model:
                    argument_option["model_arguments"].extend(["-wspr"])

        # Return the constructed argument options
        return argument_option


if __name__ == "__main__":
    # Configure the logging settings (optional)
    logging.basicConfig(
        level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s"
    )
    logger = logging.getLogger(__name__)
    # Instantiate the Satute class with your desired arguments
    iqtree_path = "iqtree2"
    satute = Satute(iqtree=iqtree_path, logger=logger)
    # Run the tool
    satute.run()
