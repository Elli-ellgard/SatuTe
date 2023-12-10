#!/usr/bin/env python
import argparse
import logging
from pathlib import Path
import pandas as pd
from ete3 import Tree
from satute_repository import (
    parse_rate_matrices_from_file_new,
    parse_state_frequencies_from_file,
    parse_substitution_model,
    parse_rate_from_model,
    parse_file_to_data_frame,
)
<<<<<<< Updated upstream
from satute_exception import InvalidDirectoryError, NoAlignmentFileError
from satute_util_new import spectral_decomposition
from satute_direction_based_saturation_test import (
=======
from satute_util import spectral_decomposition
from partial_likelihood import (
>>>>>>> Stashed changes
    single_rate_analysis,
    multiple_rate_analysis,
)
from rate_matrix import RateMatrix
from satute_rate_categories_and_alignments import (
    read_alignment_file,
    split_msa_into_rate_categories_in_place,
)
from satute_rate_categories_and_alignments import parse_category_rates
from FileHandler import FileHandler, IqTreeHandler
from satute_trees_and_subtrees import map_values_to_newick
<<<<<<< Updated upstream
=======
from satute_ostream import (
    write_results_for_category_rates,
    write_results_for_single_rate,
    write_alignment_and_indices,
)
from satute_arguments import ARGUMENT_LIST
from satute_trees_and_subtrees import modify_tree
>>>>>>> Stashed changes

# Configure the logging settings (optional)
logging.basicConfig(
    level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


class Satute:
    """Class representing Satute command-line tool for wrapping up functions of IQ-TREE."""

    def __init__(
        self,
        iqtree=None,
        model=None,
        nr=None,
        output_prefix=None,
        ufboot=None,
        boot=None,
    ):
        # IQ-TREE related attributes
        self.iqtree = iqtree
        self.iqtree_tree_file = None
        self.iqtree_handler = None

        # Directories and paths
        self.msa = None
        self.active_directory = None

        # Model and testing parameters
        self.model = model
        self.nr = nr
        self.ufboot = ufboot
        self.boot = boot
        self.alpha = 0.05

        # Miscellaneous attributes
        self.output_prefix = output_prefix
        self.input_args = []
        self.input_args_dict = {}
        self.number_rates = 1

<<<<<<< Updated upstream
        # Command-line argument configurations
        self.arguments = self._initialize_arguments()

    def _initialize_arguments(self):
        """Returns a list of dictionaries representing the command-line arguments."""
        return [
            {
                "flag": "-dir",
                "help": (
                    "Path to an existing directory containing IQ-TREE output files. "
                    "Use this option when you've already run IQ-TREE and want to avoid rerunning it. "
                    "The directory should contain essential IQ-TREE output files including the `.iqtree` file, tree file(s), and possibly a `.siteprob` file."
                ),
                "default": self.input_dir,
                "metavar": "<directory_path>",
            },
            {
                "flag": "-tree",
                "help": (
                    "Path to the input tree file in Newick or Nexus format. "
                    "This tree will be used as the basis for the saturation analysis."
                ),
                "default": self.tree,
                "metavar": "<tree_file_path>",
            },
            {
                "flag": "-msa",
                "help": (
                    "Path to the Multiple Sequence Alignment (MSA) file you wish to analyze."
                    "The MSA can be in FASTA, NEXUS, PHYLIP, or TXT format."
                ),
                "default": self.msa,
                "metavar": "<msa_file_path>",
            },
            {
                "flag": "-iqtree",
                "help": (
                    "Specifies the path to the IQ-TREE executable. If IQ-TREE is installed system-wide, "
                    "just providing the executable name (`iqtree` or `iqtree2`) will suffice. "
                    "Otherwise, give the complete path."
                ),
                "default": self.iqtree,
                "metavar": "<iqtree_path>",
            },
            {
                "flag": "-model",
                "help": (
                    "Indicates the model of sequence evolution. Common models include `GTR`, `HKY`, etc. "
                    "You can also specify rate heterogeneity and other model extensions, like `+G4` for gamma-distributed rates."
                ),
                "type": str,
                "default": self.model,
                "metavar": "<evolution_model>",
            },
            {
                "flag": "-nr",
                "help": (
                    "Number of rate categories for the model. Relevant for models with gamma-distributed rate variations. "
                    "If the `-model` option includes rate variation (e.g., `+G4`), the `-nr` should match the number in the model."
                ),
                "type": int,
                "default": self.nr,
                "metavar": "<number_of_rates>",
            },
            {
                "flag": "-ufboot",
                "help": (
                    "Number of replicates for the ultrafast bootstrap analysis. Typically, a higher number like `1000` or `5000` is used. "
                    "Ultrafast bootstrap provides rapid approximations to traditional bootstrap values."
                ),
                "type": int,
                "default": self.ufboot,
                "metavar": "<number_of_replicates>",
            },
            {
                "flag": "-boot",
                "help": (
                    "Number of replicates for traditional bootstrap analysis. This also computes a Maximum Likelihood (ML) tree and a consensus tree. "
                    "Common values are `1000` or `5000`."
                ),
                "type": int,
                "default": self.boot,
                "metavar": "<number_of_replicates>",
            },
            {
                "flag": "-alpha",
                "help": (
                    "Significance level for the saturation test. A common threshold is `0.05`, indicating a 5% significance level. "
                    "Lower values make the test more stringent."
                ),
                "type": float,
                "default": self.alpha,
                "metavar": "<significance_level>",
            },
            {
                "flag": "-edge",
                "help": (
                    "Specify a branch or edge name to focus the analysis on. Useful when you want to check saturation on a specific branch."
                ),
                "type": str,
                "default": None,
                "metavar": "<edge_name>",
            },
        ]

=======
>>>>>>> Stashed changes
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
        self.input_args_dict = vars(self.input_args)

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
        if self.input_args.dir and any(getattr(self.input_args, opt) for opt in ["msa", "tree", "model", "ufboot", "boot"]):
            raise ValueError("Cannot run Satute with -dir option combined with -msa, -tree, -model, -ufboot, or -boot.")

    def validate_msa_presence(self):
        if not self.input_args.dir and not self.input_args.msa:
            raise ValueError("MSA file is required when -dir is not used.")

<<<<<<< Updated upstream
            if self.input_args.nr is not None:
                logger.error(f"{error_str}")
                raise ValueError("Cannot run Satute with option -dir and -nr.")

            if self.input_args.ufboot is not None:
                logger.error(f"{error_str}")
                raise ValueError("Cannot run Satute with option -dir and -ufboot.")

            if self.input_args.boot is not None:
                logger.error(f"{error_str}")
                raise ValueError("Cannot run Satute with option -dir and -boot.")
        else:
            if self.input_args.msa is None:
                logger.error(
                    "Cannot run Satute without a given MSA file or a given directory path."
                )
                raise ValueError(
                    "Cannot run Satute without a given MSA file or a given directory path."
                )
            else:
                if self.input_args.tree is not None:
                    if self.input_args.model is None:
                        logger.error(
                            "Cannot run Satute with only a tree file and MSA file. The model must be specified."
                        )
                        raise ValueError(
                            "Cannot run Satute with only a tree file and MSA file. The model must be specified."
                        )
                    else:
                        if self.input_args.ufboot is not None:
                            logger.error(
                                "Cannot run the Satute modi msa+model+tree with otpion -ufboot."
                            )
                            raise ValueError(
                                "Cannot run the Satute modi msa+model+tree with otpion -ufboot."
                            )
                        if self.input_args.boot is not None:
                            logger.error(
                                "Cannot run the Satute modi msa+model+tree with otpion -boot."
                            )
                            raise ValueError(
                                "Cannot run the Satute modi msa+model+tree with otpion -boot."
                            )
=======
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
>>>>>>> Stashed changes

    def run_iqtree_workflow(self, arguments_dict):
        if arguments_dict["option"] == "dir":
            logger.info(
                "IQ-TREE will not to be needed the analysis will be done on the already existing iqtree files."
            )
        else:
            # For the other options IQ-Tree is necessary. Therefore, test if IQ-TREE exists
            self.iqtree_handler.check_iqtree_path(self.input_args.iqtree)

        if arguments_dict["option"] == "dir + site_probabilities":
            logger.info("Running Satute with site probabilities")
            logger.info(
                "IQ-TREE will be needed for the site probabilities for the corresponding rate categories."
            )
<<<<<<< Updated upstream
            number_rates = self.handle_number_rates()
=======
>>>>>>> Stashed changes

            iqtree_args = [
                "-m",
                self.input_args.model,
                "--redo",
<<<<<<< Updated upstream
                "--tree-fix",
                "-n 0",
                "-wspr",
                "--quiet",
=======
                "-blfix",
                "-wspr",
                "--quiet",
                "--keep-ident",
>>>>>>> Stashed changes
            ]

            logger.info("Used IQ-TREE options:")
            logger.info(" ".join(arguments_dict["arguments"]))
            logger.info(" ".join(iqtree_args))

            self.iqtree_handler.run_iqtree_with_arguments(
                arguments_dict["arguments"], iqtree_args
            )

        if arguments_dict["option"] == "msa":
            logger.info(
                "Running IQ-TREE with constructed arguments. Note: If no model is specified in input arguments, "
                "the best-fit model will be extracted from the log file."
            )
            logger.warning(
                "IQ-TREE will run with default options. For specific options required in the analysis, please run IQ-TREE separately."
            )

            extra_arguments = [
                "-m MF",
                "--quiet",
            ]

            logger.info("Used IQ-TREE options for Modelfinder run:")
            logger.info(" ".join(arguments_dict["arguments"]))
            logger.info(" ".join(extra_arguments))

            self.iqtree_handler.run_iqtree_with_arguments(
                arguments=arguments_dict["arguments"],
                extra_arguments=extra_arguments,
            ),

            # Update model in input arguments and re-construct arguments
            substitution_model = parse_substitution_model(
                f"{arguments_dict['msa_file']}.iqtree"
            )

<<<<<<< Updated upstream
            self.input_args.model = substitution_model
            number_rates = self.handle_number_rates()
=======
            iq_tree_parser = IqTreeParser(f"{arguments_dict['msa_file']}.iqtree")
            self.input_args.model = iq_tree_parser.parse_substitution_model()
            self.number_rates = self.handle_number_rates()
>>>>>>> Stashed changes

            # Validate and append ufboot and boot parameters to extra_arguments
            bb_arguments = self.iqtree_handler.validate_and_append_boot_arguments(
                self.input_args.ufboot, self.input_args.boot
            )

            extra_arguments = bb_arguments + [
                "-m",
                self.input_args.model,
                "--redo",
                "--quiet",
            ]

            if self.number_rates > 1:
                extra_arguments = extra_arguments + ["-wspr"]

            logger.info("Used IQ-TREE options for final run:")
            logger.info(" ".join(arguments_dict["arguments"]))
            logger.info(" ".join(extra_arguments))

            self.iqtree_handler.run_iqtree_with_arguments(
                arguments=arguments_dict["arguments"], extra_arguments=extra_arguments
            )

        if arguments_dict["option"] == "msa + model":
            logger.info("Running IQ-TREE with constructed arguments")
            logger.warning(
                "Please consider for the analysis that IQ-Tree will be running with default options."
            )
            logger.warning(
                "If specific options are required for the analysis, please run IQ-Tree separately."
            )

            self.number_rates = self.handle_number_rates()
            bb_arguments = self.iqtree_handler.validate_and_append_boot_arguments(
                self.input_args.ufboot, self.input_args.boot
            )

            extra_arguments = bb_arguments + [
                "-m",
                self.input_args.model,
                "--quiet",
            ]

            if self.number_rates > 1:
                extra_arguments = extra_arguments + ["-wspr"]

            logger.info("Used IQ-TREE options:")
            logger.info(" ".join(arguments_dict["arguments"]))
            logger.info(" ".join(extra_arguments))

            self.iqtree_handler.run_iqtree_with_arguments(
                arguments_dict["arguments"], extra_arguments
            )

        if arguments_dict["option"] == "msa + tree + model":
            logger.info("Running IQ-TREE with constructed arguments")
            logger.warning(
                """Please consider for the analysis that IQ-Tree will be running with default options.
                If specific options are required for the analysis, please run IQ-Tree separately."""
            )
            self.number_rates = self.handle_number_rates()

            iqtree_args = [
                "-m",
                self.input_args.model,
                "--quiet",
<<<<<<< Updated upstream
                "--tree-fix",
                "-n 0",
=======
                "-blfix",
                "--keep-ident",
>>>>>>> Stashed changes
            ]

            # Add the '-wspr' option if number_rates > 1
            if self.number_rates > 1:
                iqtree_args.append("-wspr")

            logger.info("Used IQ-TREE options:")
            logger.info(" ".join(arguments_dict["arguments"]))
            logger.info(" ".join(iqtree_args))

            # Call IQ-TREE with the constructed arguments
            self.iqtree_handler.run_iqtree_with_arguments(
                arguments_dict["arguments"], iqtree_args
            )

    def extract_tree(self):
        """
        Extracts the Newick representation of a tree from either a specified file or a directory.

        This method first checks if a specific tree file is provided via input arguments.
        If so, it extracts the Newick representation from that file. Otherwise, it attempts
        to find and extract the Newick string from a tree file within the specified directory.

        After obtaining the Newick string, the method then modifies the tree to ensure
        unique naming of the nodes based on a preorder traversal strategy.

        Returns:
            Tree: An ETE Tree object constructed and potentially modified from the Newick representation.

        Raises:
            ValueError: If no suitable tree file is found in the specified directory or provided via arguments.
        """
        try:
            if self.input_args.tree:
                logger.info(f"Using the user-defined tree: {self.input_args.tree}")
                newick_string = self.file_handler.get_newick_string(
                    self.input_args.tree
                )
            else:
                newick_string = self.file_handler.get_newick_string_from_args()
        except (FileNotFoundError, ValueError) as e:
            logger.error(str(e))
            raise ValueError(f"Error extracting tree: {e}")

<<<<<<< Updated upstream
    def modify_tree(self, t):
        """
        Modify the input tree by naming its nodes using a preorder traversal.

        Nodes are named as "NodeX*" where X is an incremental number.
        If a node name is purely numeric, it is preserved as 'apriori' feature of the node.

        Args:
            t (Tree): The input tree to be modified.

        Returns:
            Tree: The modified tree with updated node names.
        """
        idx = 1
        for node in t.traverse("preorder"):
            if not node.is_leaf():
                # If a node name already starts with "Node", no further modification is required.
                if node.name.startswith("Node"):
                    return t
                # Preserve numeric node names as 'apriori' for reference.
                if node.name.isdigit():
                    node.add_features(apriori=node.name)
                # Assign new node names based on the preorder traversal index.
                node.name = f"Node{idx}*"
                idx += 1
        return t

    def tree_handling(self):
        """
        Handle the extraction and subsequent modification of the tree.

        This function first extracts the tree from either the specified tree file
        or the directory. After extraction, the tree nodes are modified to
        have unique names based on a preorder traversal.

        Returns:
            Tree: The modified tree.
        """
        to_be_tested_tree = self.extract_tree()
        return self.modify_tree(to_be_tested_tree)

    def write_results_for_single_rate(self, results, to_be_tested_tree):
        """
        Writes the results for single rate to appropriate files.

        Args:
        - results (dict): Dictionary containing the results data.
        - to_be_tested_tree (Tree): The tree object containing the data to be written.
        """
        for key, results_set in results.items():
            file_name_base = (
                f"{self.input_args.msa.resolve()}_{self.input_args.alpha}.satute"
            )

            # Append the edge information to the file name if provided
            if self.input_args.edge:
                file_name_base += f"_{self.input_args.edge}"

            results_data_frame = pd.DataFrame(results_set)

            # Writing the .tree file
            with open(f"{file_name_base}.nex", "w") as tree_file:
                newick_string = to_be_tested_tree.write(format=1, features=["apriori"])
                newick_string = map_values_to_newick(newick_string, results_data_frame)
                tree_file.write("#NEXUS")
                tree_file.write("BEGIN TREES;")
                tree_file.write(f"Tree tree1 = {newick_string}")
                tree_file.write("END TREES;")

            # Writing the .csv file
            results_data_frame.to_csv(f"{file_name_base}.csv")
            logger.info("Finished writing results to files")

    def run(self):
        """
        Main entry point for running the Satute command-line tool.
        """
        # Parsing and checking input arguments and constructing IQ-TREE command-line arguments
        self.parse_input()
        self.check_input()

        arguments_dict = self.construct_arguments()

        self.run_iqtree_workflow(arguments_dict)

        # ======== Tree File Handling =========
        to_be_tested_tree = self.tree_handling()
=======
        return modify_tree(Tree(newick_string, format=1))

    def initialise_logger(self):
        log_file = f"{self.input_args.msa.resolve()}_{self.input_args.alpha}.satute.log"
        file_handler = logging.FileHandler(log_file)
        # Add the FileHandler to the logger
        logger.addHandler(file_handler)

    def validate_and_set_rate_category(self, input_category, number_rates):
        """
        Validates the input category against the number of rates and sets the rate category.
>>>>>>> Stashed changes

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
        # ======== Input =====================
        self.parse_input()
        # ======== Validation ================
        self.validate_satute_input_options()
        # ======== Arguments =================
        iq_parser_arguments_dict = self.construct_iq_tree_arguments()
        msa_file = iq_parser_arguments_dict["msa_file"]
        # ======== Logger ====================
        self.initialise_logger()
        # ======== IQ-Tree ===================
        self.run_iqtree_workflow(iq_parser_arguments_dict)
        # ======== Tree File Handling ========
        to_be_tested_tree = self.extract_tree()
        # ======== Model parameter ===========
<<<<<<< Updated upstream
        ## Get dictionary for stationary distribution and diagonal matrix of the stationary distribution
        state_frequencies, psi_matrix = parse_state_frequencies_from_file(
            f"{arguments_dict['msa_file'].resolve()}.iqtree"
        )

        ## Get rate matrix using rate parameters and stationay distribution
        rate_matrix = parse_rate_matrices_from_file_new(
            f"{arguments_dict['msa_file'].resolve()}.iqtree", state_frequencies
        )

        ## Convert representation of rate_matrix
        RATE_MATRIX = RateMatrix(rate_matrix)

        ## Calculation of the spectral decomposition of the rate matrix
=======

        # Get dictionary for stationary distribution and diagonal matrix of the stationary distribution
        iqtree_file = f"{msa_file.resolve()}.iqtree"
        satute_iq_tree_parser = IqTreeParser(iqtree_file)
        substitution_model = satute_iq_tree_parser.load_substitution_model()
        RATE_MATRIX = RateMatrix(substitution_model.rate_matrix)

        # Calculation of the spectral decomposition of the rate matrix
>>>>>>> Stashed changes
        (
            array_left_eigenvectors,
            array_right_eigenvectors,
            multiplicity,
<<<<<<< Updated upstream
        ) = spectral_decomposition(RATE_MATRIX.rate_matrix, psi_matrix)

        ## Get number of rate categories in case of a +G or +R model
        number_rates = self.handle_number_rates()

        # ======== Multiple Sequence Alignment
        alignment = read_alignment_file(arguments_dict["msa_file"].resolve())

        # ========  Test for Branch Saturation =========
=======
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
            
        # ======== Logger ====================            
>>>>>>> Stashed changes
        logger.info(
            f"""
            Running tests and initial IQ-Tree with configurations:
            Mode {iq_parser_arguments_dict['option']}
            Model: {self.input_args.model}
            Alpha: {self.input_args.alpha}
<<<<<<< Updated upstream
            Running Saturation Test on file: {arguments_dict['msa_file'].resolve()}
            Number of rate categories: {number_rates}
            Options for Initial IQ-Tree run: {arguments_dict['option']}
=======
            Running Saturation Test on file: {msa_file.resolve()}
            Number of rate categories: {substitution_model.number_rates}
            Considered rate category: {rate_category}
            Options for Initial IQ-Tree run: {iq_parser_arguments_dict['option']}
>>>>>>> Stashed changes
            Multiplicity: {multiplicity}
            Run test for saturation for each branch and category with {number_rates} rate categories
            Results will be written to the directory:{self.active_directory.name}
            """
        )

<<<<<<< Updated upstream
        if number_rates == 1:
=======
        # ======== Multiple Sequence Alignment =========
        alignment = read_alignment_file(msa_file)
        # ========  Test for Branch Saturation =========        
        if substitution_model.number_rates == 1:
>>>>>>> Stashed changes
            results = single_rate_analysis(
                to_be_tested_tree,
                alignment,
                RATE_MATRIX,
                state_frequencies,
                array_left_eigenvectors,
                array_right_eigenvectors,
                multiplicity,
                self.input_args.alpha,
                self.input_args.edge,
            )

            self.write_results_for_single_rate(results, to_be_tested_tree)
        else:
<<<<<<< Updated upstream
            site_probability = parse_file_to_data_frame(
                f"{self.input_args.msa.resolve()}.siteprob"
            )
            per_rate_category_alignment = split_msa_into_rate_categories_in_place(
                site_probability, alignment
=======
            site_probability_file = f"{msa_file.resolve()}.siteprob"
            site_probability = parse_file_to_data_frame(site_probability_file)
            per_rate_category_alignment = split_msa_into_rate_categories_in_place(
                site_probability, alignment, rate_category
            )
            categorized_sites = build_categories_by_sub_tables(site_probability)

            write_alignment_and_indices(
                per_rate_category_alignment, categorized_sites, self.input_args
>>>>>>> Stashed changes
            )

            for rate in per_rate_category_alignment.keys():
                logger.info(
                    f"""Category Rate {rate}, Site per category {len(per_rate_category_alignment[rate][0].seq)}"""
                )

            category_rates_factors = parse_category_rates(
                f"{self.input_args.msa.resolve()}.iqtree", number_rates
            )

            results = multiple_rate_analysis(
                to_be_tested_tree,
                category_rates_factors,
                RATE_MATRIX,
                state_frequencies,
                array_left_eigenvectors,
                array_right_eigenvectors,
                multiplicity,
                per_rate_category_alignment,
                self.input_args.alpha,
                self.input_args.edge,
            )

            self.write_results_for_category_rates(results)

<<<<<<< Updated upstream
    def rename_internal_nodes(self, t: Tree) -> Tree:
        """
        Rename the internal nodes of the tree to a standardized format.

        Args:
        - t (Tree): The tree to be modified.

        Returns:
        - Tree: The modified tree.
        """
        idx = 1  # Start index from 1 as per your requirement
        for node in t.traverse("preorder"):
            # Process only inner nodes
            if not node.is_leaf():
                # If node name starts with "Node", skip this node
                if node.name.startswith("Node"):
                    continue
                # Check if node name is just a number
                if node.name.isdigit():
                    node.add_features(apriori=node.name)
                # Set inner node names as "Node<index>*"
                node.name = f"Node{idx}*"
                idx += 1
        return t

    def write_results_for_category_rates(self, results):
        """
        Writes the results for category rates to appropriate files.

        Args:
        - results (dict): Dictionary containing the results data.
        """
        for key, results_set in results.items():
            file_name = (
                f"{self.input_args.msa.resolve()}_{key}_{self.input_args.alpha}.satute"
            )

            # Append the edge information to the file name if provided
            if self.input_args.edge:
                file_name = file_name + f"_{self.input_args.edge}"

            results_data_frame = pd.DataFrame(results_set["result_list"])

            if "rescaled_tree" in results_set and "result_list" in results_set:
                with open(f"{file_name}.nex", "w") as tree_file_writer:
                    newick_string = results_set["rescaled_tree"].write(format=1)
                    newick_string = map_values_to_newick(
                        newick_string, results_data_frame
                    )
                    # Writing the .tree file
                    tree_file_writer.write("#NEXUS\n")
                    tree_file_writer.write("BEGIN TREES;\n")
                    tree_file_writer.write(f"Tree tree1 = {newick_string}\n")
                    tree_file_writer.write("END TREES;\n")

                results_data_frame.to_csv(f"{file_name}.csv")
        logger.info("Finished writing results for category rates to files")

    def handle_number_rates(self):
        number_rates = 1
        if self.input_args.nr:
            if self.input_args.model:
                number_rates_model = parse_rate_from_model(self.input_args.model)
                if self.input_args.nr == number_rates_model:
                    number_rates = self.input_args.nr
                else:
                    raise ValueError(
                        "Input number of rates is unequal to number of rates of the model!"
                    )
            else:
                number_rates = self.input_args.nr
        else:
            if self.input_args.model:
                number_rates = parse_rate_from_model(self.input_args.model)

        return number_rates

    def convert_paths_to_objects(self):
        """
        Convert input paths to Path objects for easier handling.
        """
        paths_to_convert = ["dir", "iqtree", "msa", "tree"]
        for path_key in paths_to_convert:
            path_value = getattr(self.input_args, path_key)
            if path_value:
                setattr(self.input_args, path_key, Path(path_value))
=======
    def handle_number_rates(self):
        self.number_rates = 1
        if self.input_args.model:
            self.number_rates = parse_rate_from_model(self.input_args.model)
        return self.number_rates
>>>>>>> Stashed changes

    def initialize_working_context(self):
        if self.input_args.msa:
            self.active_directory = self.input_args.msa.parent
        elif self.input_args.dir:
            self.active_directory = self.input_args.dir

    def initialize_handlers(self):
        """
        Initialize the FileHandler and IqTreeHandler.
        """
        self.file_handler = FileHandler(self.active_directory)
        self.iqtree_handler = IqTreeHandler(self.input_args.iqtree)

    def construct_iq_tree_arguments(self):
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
            self.input_args.msa = Path(
                self.file_handler.find_file_by_suffix(msa_file_types)
            )

            argument_option = {
                "option": "dir",
                "msa_file": self.input_args.msa,
                "arguments": ["-s", str(self.input_args.msa.resolve())],
            }

            # Find the tree file in the directory
            tree_file = self.file_handler.find_file_by_suffix(tree_file_types)
            argument_option["arguments"].extend(["-te", str(tree_file)])

            # Find .iqtree file in the directory
            self.iqtree_tree_file = self.file_handler.find_file_by_suffix({".iqtree"})
            self.input_args.model = parse_substitution_model(self.iqtree_tree_file)
            self.number_rates = self.handle_number_rates()

            # TODO Test when there is no site probabilities file in the directory
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
    # Instantiate the Satute class with your desired arguments
<<<<<<< Updated upstream
    iqtree_path = "iqtree2"
    input_directory = None
    output_prefix = None
    model = None
    num_rate_categories = None
    ufboot_replicates = None
    boot_replicates = None
    alpha = None

    satute = Satute(
        iqtree=iqtree_path,
        input_dir=input_directory,
        model=model,
        nr=num_rate_categories,
        output_prefix=output_prefix,
        ufboot=ufboot_replicates,
        boot=boot_replicates,
    )
    # Run the tool
=======
    satute = Satute()
>>>>>>> Stashed changes
    satute.run()
