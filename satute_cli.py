#!/usr/bin/env python
import argparse
import os
import logging
from pathlib import Path
import re
import pandas as pd
from ete3 import Tree
from satute_repository import parse_rate_matrices_from_file, parse_state_frequencies
from satute_exception import InvalidDirectoryError, NoAlignmentFileError
from satute_util_new import spectral_decomposition, parse_file_to_data_frame

from satute_direction_based_saturation_test import (
    single_rate_analysis,
    multiple_rate_analysis,
    RateMatrix,
)
from satute_rate_categories_and_alignments import (
    read_alignment_file,
    split_msa_into_rate_categories_in_place,
)
from satute_rate_categories_and_alignments import parse_category_rates
from FileHandler import FileHandler, IqTreeHandler
from satute_trees_and_subtrees import name_nodes_by_level_order

# Configure the logging settings (optional)
logging.basicConfig(
    level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def parse_rate_from_model(model):
    # Find the index of '+G' and '+R' in the model string
    plus_g_index = model.find("+G")
    plus_r_index = model.find("+R")

    if plus_g_index != -1 and plus_r_index != -1:
        raise ValueError("Cannot use +G and +R")

    if plus_g_index != -1:
        rate_start_index = plus_g_index + 2
    elif plus_r_index != -1:
        rate_start_index = plus_r_index + 2
    else:
        return 1  # default number_rates = 1 if no +G or +R model

    try:
        # Extract the substring after '+G'
        number = model[rate_start_index:]

        # Parse the extracted substring as an integer
        rate = int(number)

        return rate
    except ValueError:
        # If '+G' is not found or the number after '+G' is not a valid integer
        # Return None or an appropriate value for error handling
        logger.info("Could not find a Gamma distribution in the model.")
        logger.info("So not rates were found.")
        return 1


def parse_substitution_model(file_path):
    try:
        with open(file_path, "r") as file:
            content = file.read()
            for line in content.splitlines():
                if "Best-fit model according to BIC:" in line:
                    model_string = line.split(":")[1].strip()
                    return model_string
                if "Model of substitution: GTR+F+G4" in line:
                    model_string = line.split(":")[1].strip()
                    return model_string
            raise ValueError("Could not parse the substitution model from the file.")
    except (IOError, ValueError):
        # If the file cannot be read or ':' is not found in the content
        # Return None or an appropriate value for error handling
        raise ValueError("Could not parse the substitution model from the file.")


class Satute:
    """Class representing Satute command-line tool for wrapping up functions of IQ-TREE."""

    def __init__(
        self,
        iqtree=None,
        input_dir=None,
        model=None,
        nr="None",
        output_prefix=None,
        ufboot=None,
        boot=None,
    ):
        self.iqtree = iqtree
        self.iqtree_tree_file = None
        self.input_dir = None
        self.tree = None
        self.msa = None
        self.model = model
        self.nr = nr
        self.output_prefix = output_prefix
        self.ufboot = ufboot
        self.boot = boot
        self.input_args = []
        self.input_args_dict = {}
        self.site_probabilities_file = None
        self.alpha = 0.05
        self.number_rates = 1

        self.arguments = [
            {
                "flag": "-dir",
                "help": "Path to input directory",
                "default": self.input_dir,
                "metavar": "<file_name>",
            },
            {
                "flag": "-tree",
                "help": "Path to input tree file",
                "default": self.tree,
                "metavar": "<file_name>",
            },
            {
                "flag": "-msa",
                "help": "Path to MSA",
                "default": self.msa,
                "metavar": "<file_name>",
            },
            {
                "flag": "-iqtree",
                "help": "Path to IQ-TREE",
                "default": self.iqtree,
                "metavar": "<file_name>",
            },
            {
                "flag": "-model",
                "help": "Model of evolution",
                "type": str,
                "default": self.model,
                "metavar": "<str>",
            },
            {
                "flag": "-nr",
                "help": "Number of rate categories",
                "type": int,
                "default": self.nr,
                "metavar": "<num>",
            },
            {
                "flag": "-ufboot",
                "help": "Replicates for ultrafast bootstrap (>=1000)",
                "type": int,
                "default": self.ufboot,
                "metavar": "<num>",
            },
            {
                "flag": "-boot",
                "help": "Replicates for bootstrap + ML tree + consensus tree",
                "type": int,
                "default": self.boot,
                "metavar": "<num>",
            },
            {
                "flag": "-alpha",
                "help": "significance level of the test",
                "type": float,
                "default": self.alpha,
                "metavar": "<num>",
            },
            {
                "flag": "-edge",
                "help": "edge to be tested",
                "type": str,
                "default": None,
                "metavar": "str",
            },
        ]
        self.iqtree_handler = None
        self.active_directory = None

    def parse_input(self):
        """Parse command-line arguments."""
        parser = argparse.ArgumentParser(description="Satute")
        for argument in self.arguments:
            flag = argument["flag"]
            argument.pop("flag")
            parser.add_argument(flag, **argument)
        self.input_args = parser.parse_args()
        self.input_args_dict = vars(self.input_args)

    def run_iqtree_workflow(self, arguments_dict):
        if arguments_dict["option"] == "dir":
            logger.info("Running Satute without iqtree run")

        if arguments_dict["option"] == "dir + site_probabilities":
            logger.info("Running Satute with site probabilities")
            # TO call Iqtree only for .sibeprob file

        if arguments_dict["option"] == "msa + tree + model":
            number_rates = self.handle_number_rates()

            if number_rates > 1:
                self.iqtree_handler.run_iqtree_with_arguments(
                    arguments_dict["arguments"],
                    [
                        "-m",
                        self.input_args.model,
                        "--redo",
                        "--quiet",
                        "-blfix",
                        "-n 0",
                        "-wspr",
                    ],
                )
            else:
                self.iqtree_handler.run_iqtree_with_arguments(
                    arguments_dict["arguments"],
                    [
                        "-m",
                        self.input_args.model,
                        "--redo",
                        "--quiet",
                        "-blfix",
                        "-n 0",
                    ],
                )

        elif arguments_dict["option"] == "msa + tree":
            logger.error(
                "Cannot run Satute with only a tree file and MSA file. The model must be specified."
            )
            raise ValueError(
                "Cannot run Satute with only a tree file and MSA file. The model must be specified."
            )

        elif arguments_dict["option"] == "msa":
            logger.info("Running IQ-TREE with constructed arguments")
            logger.info(
                "If no model specified in input arguments, extract best model from log file"
            )
            # Validate and append ufboot and boot parameters to extra_arguments
            bb_arguments = self.iqtree_handler.validate_and_append_boot_arguments(
                self.input_args.ufboot, self.input_args.boot
            )

            logger.info(" ".join(arguments_dict["arguments"]))
            logger.info(" ".join(bb_arguments))

            self.iqtree_handler.run_iqtree_with_arguments(
                arguments=arguments_dict["arguments"],
                extra_arguments=bb_arguments
                + [
                    "-m MF",
                    "--redo",
                    "--quiet",
                ],
            ),

            # Update model in input arguments and re-construct arguments
            substitution_model = parse_substitution_model(
                f"{arguments_dict['msa_file']}.iqtree"
            )

            self.input_args.model = substitution_model
            number_rates = self.handle_number_rates()

            extra_arguments = bb_arguments + [
                "-m",
                self.input_args.model,
                "--redo",
                "--quiet",
            ]

            if number_rates > 1:
                self.iqtree_handler.run_iqtree_with_arguments(
                    arguments=arguments_dict["arguments"],
                    extra_arguments=extra_arguments + ["-wspr"],
                )
            else:
                self.iqtree_handler.run_iqtree_with_arguments(
                    arguments=arguments_dict["arguments"],
                )

        elif arguments_dict["option"] == "msa + model":
            number_rates = self.handle_number_rates()
            bb_arguments = self.iqtree_handler.validate_and_append_boot_arguments(
                self.input_args.ufboot, self.input_args.boot
            )

            extra_arguments = bb_arguments + [
                "-m",
                self.input_args.model,
                "--redo",
                "--quiet",
            ]

            if number_rates > 1:
                extra_arguments = extra_arguments + ["-wspr"]

            self.iqtree_handler.run_iqtree_with_arguments(
                arguments=arguments_dict["arguments"], extra_arguments=extra_arguments
            )

    def tree_handling(self):
        to_be_tested_tree = None
        if self.input_args.tree:
            logger.info(f"Using the already defined tree: {self.input_args.tree}")
            with open(self.input_args.tree, "r") as file:
                to_be_tested_tree = Tree(file.readlines()[0], format=1)
            to_be_tested_tree = self.modify_tree(to_be_tested_tree)
        else:
            newick_string = self.file_handler.get_newick_string_from_args()
            if newick_string is None:
                raise ValueError("No tree file found in directory")
            to_be_tested_tree = Tree(newick_string, format=1)
            to_be_tested_tree = self.modify_tree(to_be_tested_tree)

        return to_be_tested_tree
        # ======== End Tree File Handling =========

    def run(self):
        """Main entry point for running the Satute command-line tool."""
        # TODO Change number rated to None and test
        number_rates = 1
        dimension = 4
        # Parsing input arguments and constructing IQ-TREE command-line arguments
        self.parse_input()

        arguments_dict = self.construct_arguments()

        try:
            (
                path_to_exe,
                self.input_args.iqtree,
            ) = self.iqtree_handler.get_iqtree_version(self.input_args.iqtree)
            print(f"iqtree version: {self.input_args.iqtree}")
            print(f"Path to iqtree executable: {path_to_exe}")
        except Exception as e:
            print(f"Error: {e}")

        # TODO Check for model of substitution
        self.run_iqtree_workflow(arguments_dict)

        rate_matrix, psi_matrix = parse_rate_matrices_from_file(
            f"{arguments_dict['msa_file'].resolve()}.iqtree"
        )

        # TODO Change state frequency function
        state_frequencies = parse_state_frequencies(
            f"{arguments_dict['msa_file'].resolve()}.iqtree", dimension=dimension
        )

        # ======== Tree File Handling =========
        to_be_tested_tree = self.tree_handling()
        # ======== End Tree Handling =========
        RATE_MATRIX = RateMatrix(rate_matrix)

        (
            array_left_eigenvectors,
            array_right_eigenvectors,
            multiplicity,
        ) = spectral_decomposition(RATE_MATRIX.rate_matrix, psi_matrix)

        alignment = read_alignment_file(arguments_dict["msa_file"].resolve())
        number_rates = self.handle_number_rates()

        # ======== Saturation Test =========
        logger.info(
            f"""
            Running tests and initial IQ-Tree with configurations:
            Mode {arguments_dict['option']}
            Model: {self.input_args.model}
            Alpha: {self.input_args.alpha}
            Running Saturation Test on file: {arguments_dict['msa_file'].resolve()}
            Number of rate categories: {number_rates}
            Options for Initial IQ-Tree run: {arguments_dict['option']}
            Multiplicity: {multiplicity}
            Run test for saturation for each branch and category with {number_rates} rate categories
            Results will be written to the directory:{self.active_directory.name}
            """
        )

        # TODO Test when there is no site probabilities file in the directory
        to_be_tested_tree.write(
            format=1, outfile=f"{arguments_dict['msa_file'].resolve()}_satute_tree.tree"
        )

        if number_rates == 1:
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

            for key, results_set in results.items():
                pd.DataFrame(results_set).to_csv(
                    f"{self.input_args.msa.resolve()}_{self.input_args.alpha}_satute.csv"
                )

        else:
            site_probability = parse_file_to_data_frame(
                f"{self.input_args.msa.resolve()}.siteprob"
            )
            per_rate_category_alignment = split_msa_into_rate_categories_in_place(
                site_probability, alignment
            )

            for rate in per_rate_category_alignment.keys():
                logger.info(
                    f"""Category Rate {rate}, Site per category {len(per_rate_category_alignment[rate][0].seq)}"""
                )

            category_rates_factors = parse_category_rates(
                f"{self.input_args.msa.resolve()}.iqtree"
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

            for key, results_set in results.items():
                to_be_tested_tree.write(
                    "newick",
                    f"{self.input_args.msa.resolve()}_satute_rate_{key}_{self.input_args.alpha}_.tree",
                    format=1,
                )
                pd.DataFrame(results_set).to_csv(
                    f"{self.input_args.msa.resolve()}_satute_rate_{key}_{self.input_args.alpha}_.csv"
                )

        logger.info("Finished running Satute")

    def modify_tree(self, t):
        # Initialize counter for preorder traversal
        idx = 1  # Start index from 1 as per your requirement
        for node in t.traverse("preorder"):
            # Process only inner nodes
            if not node.is_leaf():
                # If node name starts with "Node", stop the modifications
                if node.name.startswith("Node"):
                    return t
                # Check if node name is just a number
                if node.name.isdigit():
                    node.add_features(apriori_knowledge=node.name)
                # Set inner node names as "Node<index>*"
                node.name = "Node" + str(idx) + "*"
                print(node.name)
                idx += 1
        return t

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
        """Convert input paths to Path objects for easier handling."""
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
        """Initialize the FileHandler and IqTreeHandler."""
        self.file_handler = FileHandler(self.active_directory)
        self.iqtree_handler = IqTreeHandler(self.input_args.iqtree)

    def validate_directory(self):
        """Validate the input directory."""
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
        tree_file_types = {".treefile", ".nex", ".nwk", ".tree"}

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
                    "arguments": ["-s", self.input_args.msa.resolve()],
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

            # Find iqtree file in the directory
            self.iqtree_tree_file = self.file_handler.find_file_by_suffix({".iqtree"})
            # Check if iqtree file was found
            if self.iqtree_tree_file:
                substitution_model = parse_substitution_model(self.iqtree_tree_file)
                self.input_args.model = substitution_model
            else:
                raise FileNotFoundError("No iqtree file found in directory")

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

        # If tree argument is provided, get newick string directly


if __name__ == "__main__":
    # Instantiate the Satute class with your desired arguments
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
    satute.run()
