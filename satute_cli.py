#!/usr/bin/env python
import argparse
import os
import logging
from pathlib import Path
from ete3 import Tree
from satute_repository import (
    parse_substitution_model,
    parse_rate_from_model,
    parse_file_to_data_frame,
    IqTreeParser,
)
from satute_exception import InvalidDirectoryError, NoAlignmentFileError
from satute_util import spectral_decomposition
from partial_likelihood import (
    single_rate_analysis,
    multiple_rate_analysis,
)
from rate_matrix import RateMatrix
from satute_categories import (
    read_alignment_file,
    split_msa_into_rate_categories_in_place,
    build_categories_by_sub_tables,
)
from file_handler import FileHandler, IqTreeHandler
from satute_trees_and_subtrees import map_values_to_newick
from satute_ostream import (
    write_results_for_category_rates,
    write_results_for_single_rate,
)
from satute_ostream import write_alignment_and_indices


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
        input_dir=None,
        model=None,
        category=None,
        output_prefix=None,
        ufboot=None,
        boot=None,
    ):
        # IQ-TREE related attributes
        self.iqtree = iqtree
        self.iqtree_tree_file = None
        self.iqtree_handler = None

        # Directories and paths
        self.input_dir = None
        self.tree = None
        self.msa = None
        self.site_probabilities_file = None
        self.active_directory = None

        # Model and testing parameters
        self.model = model
        self.category = category
        self.ufboot = ufboot
        self.boot = boot
        self.alpha = 0.05

        # Miscellaneous attributes
        self.output_prefix = output_prefix
        self.input_args = []
        self.input_args_dict = {}
        self.number_rates = 1

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
                "flag": "-category",
                "help": (
                    "Rate categories of interest. Relevant for models with gamma-distributed rate variations or FreeRate model. "
                    "If the `-model` option includes rate variation (e.g., `+G4`), the `-category` should be a number between 1 and 4."
                ),
                "type": int,
                "default": self.category,
                "metavar": "<rate_category>",
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

    def parse_input(self):
        """
        Parse command-line arguments.
        """
        parser = argparse.ArgumentParser(description="Satute")
        for argument in self.arguments:
            flag = argument["flag"]
            argument.pop("flag")
            parser.add_argument(flag, **argument)
        self.input_args = parser.parse_args()
        self.input_args_dict = vars(self.input_args)

    def check_input(self):
        """
        Check allowed option combinations.
        """
        if self.input_args.dir is not None:
            error_str = "The option -dir is for specifying a path to an existing directory containing IQ-TREE output files and cannot be run with the options msa, tree, model, nr, ufboot, boot."
            if self.input_args.msa is not None:
                logger.error(f"{error_str}")
                raise ValueError("Cannot run Satute with option -dir and -msa.")

            if self.input_args.tree is not None:
                logger.error(f"{error_str}")
                raise ValueError("Cannot run Satute with option -dir and -tree.")

            if self.input_args.model is not None:
                logger.error(f"{error_str}")
                raise ValueError("Cannot run Satute with option -dir and -model.")

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
                                "Cannot run the Satute modi msa+model+tree with option -ufboot."
                            )
                            raise ValueError(
                                "Cannot run the Satute modi msa+model+tree with option -ufboot."
                            )
                        if self.input_args.boot is not None:
                            logger.error(
                                "Cannot run the Satute modi msa+model+tree with option -boot."
                            )
                            raise ValueError(
                                "Cannot run the Satute modi msa+model+tree with option -boot."
                            )
                if (
                    self.input_args.model is not None
                    and self.input_args.category is not None
                ):
                    number_rates = self.handle_number_rates()
                    if (
                        self.input_args.category < 1
                        or self.input_args.category > number_rates
                    ):
                        logger.error("Choosen category of interest is out of range.")
                        raise ValueError("Chosen category of interest is out of range.")

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
            # number_rates = self.handle_number_rates()

            iqtree_args = [
                "-m",
                self.input_args.model,
                "--redo",
                "-blfix"
                # "--tree-fix",
                # "-n 0",
                "-wspr",
                "--quiet" "--keep-ident",
            ]

            logger.info("Used IQ-TREE options:")
            logger.info(" ".join(arguments_dict["arguments"]))
            logger.info(" ".join(iqtree_args))

            self.iqtree_handler.run_iqtree_with_arguments(
                arguments_dict["arguments"], iqtree_args
            )

        elif arguments_dict["option"] == "msa":
            logger.info("Running IQ-TREE with constructed arguments")
            logger.info(
                "If no model is specified in input arguments, best-fit model will be extracted from log file."
            )
            logger.warning(
                "Please consider for the analysis that IQ-Tree will be running with default options."
            )
            logger.warning(
                "If specific options are required for the analysis, please run IQ-Tree separately."
            )

            extra_arguments = [
                "-m MF",
                "--quiet",
            ]

            logger.info("Used IQ-TREE options for Model Finder run:")
            logger.info(" ".join(arguments_dict["arguments"]))
            logger.info(" ".join(extra_arguments))

            self.iqtree_handler.run_iqtree_with_arguments(
                arguments=arguments_dict["arguments"],
                extra_arguments=extra_arguments,
            ),

            # Update model in input arguments and re-construct arguments

            iq_tree_parser = IqTreeParser(f"{arguments_dict['msa_file']}.iqtree")
            self.input_args.model = iq_tree_parser.parse_substitution_model()
            number_rates = self.handle_number_rates()

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

            if number_rates > 1:
                extra_arguments = extra_arguments + ["-wspr"]

            logger.info("Used IQ-TREE options for final run:")
            logger.info(" ".join(arguments_dict["arguments"]))
            logger.info(" ".join(extra_arguments))

            self.iqtree_handler.run_iqtree_with_arguments(
                arguments=arguments_dict["arguments"], extra_arguments=extra_arguments
            )

        elif arguments_dict["option"] == "msa + model":
            logger.info("Running IQ-TREE with constructed arguments")
            logger.warning(
                "Please consider for the analysis that IQ-Tree will be running with default options."
            )
            logger.warning(
                "If specific options are required for the analysis, please run IQ-Tree separately."
            )

            number_rates = self.handle_number_rates()
            bb_arguments = self.iqtree_handler.validate_and_append_boot_arguments(
                self.input_args.ufboot, self.input_args.boot
            )

            extra_arguments = bb_arguments + [
                "-m",
                self.input_args.model,
                "--quiet",
                "--keep-ident",
            ]

            if number_rates > 1:
                extra_arguments = extra_arguments + ["-wspr"]

            logger.info("Used IQ-TREE options:")
            logger.info(" ".join(arguments_dict["arguments"]))
            logger.info(" ".join(extra_arguments))

            self.iqtree_handler.run_iqtree_with_arguments(
                arguments=arguments_dict["arguments"], extra_arguments=extra_arguments
            )

        # elif arguments_dict["option"] == "msa + tree":
        #     logger.error(
        #         "Cannot run Satute with only a tree file and MSA file. The model must be specified."
        #     )
        #     raise ValueError(
        #         "Cannot run Satute with only a tree file and MSA file. The model must be specified."
        #     )

        if arguments_dict["option"] == "msa + tree + model":
            logger.info("Running IQ-TREE with constructed arguments")
            logger.warning(
                "Please consider for the analysis that IQ-Tree will be running with default options."
            )
            logger.warning(
                "If specific options are required for the analysis, please run IQ-Tree separately."
            )
            number_rates = self.handle_number_rates()

            iqtree_args = [
                "-m",
                self.input_args.model,
                "--quiet",
                "-blfix",
                "--keep-ident",
                # -blfix Fix branch lengths of user tree passed via -te
                # --tree-fix,Fix -t tree (no tree search performed)
                # "-n 0", Fix number of iterations to stop (default: OFF)
            ]

            # Add the '-wspr' option if number_rates > 1
            if number_rates > 1:
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
        Extract the Newick representation of a tree from a specified file or directory.

        If a specific tree file is provided via the input arguments, it uses that file.
        Otherwise, it attempts to retrieve the tree from the specified directory.

        Returns:
            Tree: An ETE Tree object constructed from the Newick representation.

        Raises:
            ValueError: If no tree file is found in the directory.
        """
        if self.input_args.tree:
            logger.info(f"Using the already defined tree: {self.input_args.tree}")
            with open(self.input_args.tree, "r") as file:
                tree_string = ""
                for line in file:
                    tree_string += line.strip()
                    if ";" in line:  # Check for the end of the Newick representation
                        break
                return Tree(tree_string, format=1)
        else:
            newick_string = self.file_handler.get_newick_string_from_args()
            if newick_string is None:
                raise ValueError("No tree file found in directory")
            return Tree(newick_string, format=1)

    def modify_tree(self, t):
        """
        Modify the input tree by naming its nodes using a preorder traversal.

        Nodes are named as "NodeX*" where X is an incremental number.
        If a node name is purely numeric, it is preserved as 'apriorism' feature of the node.

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

    def run(self):
        """
        Main entry point for running the Satute command-line tool.
        """
        # Parsing and checking input arguments and constructing IQ-TREE command-line arguments
        self.parse_input()
        self.check_input()
        arguments_dict = self.construct_arguments()

        log_file = f"{self.input_args.msa.resolve()}_{self.input_args.alpha}.satute.log"
        file_handler = logging.FileHandler(log_file)
        # Add the FileHandler to the logger
        logger.addHandler(file_handler)

        self.run_iqtree_workflow(arguments_dict)
        # ======== Tree File Handling =========
        to_be_tested_tree = self.tree_handling()

        # ======== Model parameter ===========
        ## Get dictionary for stationary distribution and diagonal matrix of the stationary distribution
        iq_tree_file_path = f"{arguments_dict['msa_file'].resolve()}.iqtree"
        satute_iq_tree_parser = IqTreeParser(iq_tree_file_path)
        substitution_model = satute_iq_tree_parser.load_substitution_model()
        ## Convert representation of rate_matrix
        RATE_MATRIX = RateMatrix(substitution_model.rate_matrix)
        ## Calculation of the spectral decomposition of the rate matrix
        (
            array_left_eigenvectors,
            array_right_eigenvectors,
            multiplicity,
        ) = spectral_decomposition(
            substitution_model.rate_matrix, substitution_model.phi_matrix
        )
        ## Get number of rate categories in case of a +G or +R model
        # Consider a specific rate category
        rate_category = "all"
        if self.input_args.category is not None:
            if (
                self.input_args.category < 1
                or self.input_args.category > substitution_model.number_rates
            ):
                logger.error("Chosen category of interest is out of range.")
                raise ValueError("Chosen category of interest is out of range.")
            else:
                rate_category = str(self.input_args.category)
        # ======== Multiple Sequence Alignment
        alignment = read_alignment_file(arguments_dict["msa_file"].resolve())
        # ========  Test for Branch Saturation =========
        logger.info(
            f"""
            Running tests and initial IQ-Tree with configurations:
            Mode {arguments_dict['option']}
            Model: {self.input_args.model}
            Alpha: {self.input_args.alpha}
            Running Saturation Test on file: {arguments_dict['msa_file'].resolve()}
            Number of rate categories: {substitution_model.number_rates}
            Considered rate category: {rate_category}
            Options for Initial IQ-Tree run: {arguments_dict['option']}
            Multiplicity: {multiplicity}
            Run test for saturation for each branch and category with {substitution_model.number_rates} rate categories
            Results will be written to the directory:{self.active_directory.name}
            """
        )

        if substitution_model.number_rates == 1:
            results = single_rate_analysis(
                to_be_tested_tree,
                alignment,
                RATE_MATRIX,
                substitution_model.state_frequencies,
                array_left_eigenvectors,
                array_right_eigenvectors,
                multiplicity,
                self.input_args.alpha,
                self.input_args.edge,
            )

            write_results_for_single_rate(
                results,
                self.input_args,
                to_be_tested_tree,
                map_values_to_newick,
                logger,
            )

        else:
            site_probability = parse_file_to_data_frame(
                f"{self.input_args.msa.resolve()}.siteprob"
            )

            per_rate_category_alignment = split_msa_into_rate_categories_in_place(
                site_probability, alignment, rate_category
            )

            categorized_sites = build_categories_by_sub_tables(site_probability)

            write_alignment_and_indices(
                per_rate_category_alignment, categorized_sites, self.input_args
            )

            for rate in per_rate_category_alignment.keys():
                logger.info(
                    f"""Category Rate {rate}, Site per category {len(per_rate_category_alignment[rate][0].seq)}"""
                )

            results = multiple_rate_analysis(
                to_be_tested_tree,
                substitution_model.category_rates,
                RATE_MATRIX,
                substitution_model.state_frequencies,
                array_left_eigenvectors,
                array_right_eigenvectors,
                multiplicity,
                per_rate_category_alignment,
                self.input_args.alpha,
                self.input_args.edge,
            )

            write_results_for_category_rates(
                results, self.input_args, map_values_to_newick, logger
            )

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

    def handle_number_rates(self):
        number_rates = 1
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

    def validate_directory(self):
        """
        Validate the input directory.
        """
        if not self.input_args.dir.is_dir():
            raise InvalidDirectoryError("Input directory does not exist")
        if not os.listdir(self.input_args.dir):
            raise InvalidDirectoryError("Input directory is empty")

    def construct_arguments(self):
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

        self.convert_paths_to_objects()
        self.initialize_handlers()

        argument_option = {}

        if self.input_args.dir:
            self.active_directory = self.input_args.dir

            # Check if the input directory exists
            self.validate_directory()

            # Find msa file in the directory
            self.input_args.msa_file = self.file_handler.find_file_by_suffix(
                msa_file_types
            )

            # Check if a sequence alignment file was found
            if self.input_args.msa_file:
                self.input_args.msa = Path(self.input_args.msa_file)
                argument_option = {
                    "option": "dir",
                    "msa_file": self.input_args.msa,
                    "arguments": ["-s", str(self.input_args.msa.resolve())],
                }
            else:
                raise NoAlignmentFileError("No multiple sequence alignment file found")

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

            self.number_rates = self.handle_number_rates()

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

        # If tree argument is provided, get newick string directly


if __name__ == "__main__":
    # Instantiate the Satute class with your desired arguments
    iqtree_path = "iqtree2"
    input_directory = None
    output_prefix = None
    model = None
    rate_category = None
    ufboot_replicates = None
    boot_replicates = None
    alpha = None

    satute = Satute(
        iqtree=iqtree_path,
        input_dir=input_directory,
        model=model,
        category=rate_category,
        output_prefix=output_prefix,
        ufboot=ufboot_replicates,
        boot=boot_replicates,
    )
    # Run the tool
    satute.run()
