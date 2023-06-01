#!/usr/bin/env python
import sys
import argparse
import os
import logging
import subprocess
from pathlib import Path
import re
import subprocess
import os
from satute_util import saturation_test_cli


# Configure the logging settings (optional)
logging.basicConfig(
    level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


class InvalidDirectoryError(Exception):
    """Exception raised when the input directory does not exist."""

    pass


class NoAlignmentFileError(Exception):
    """Exception raised when no multiple sequence alignment file is found."""

    pass


class Satute:
    """Class representing Satute command-line tool for wrapping up functions of IQ-TREE."""

    def __init__(self, iqtree, input_dir, model, nr, output_prefix, ufboot, boot):
        self.iqtree = iqtree
        self.input_dir = Path(
            input_dir
        )  # Convert to Path object for convenient file path manipulation
        self.model = model
        self.nr = nr
        self.output_prefix = output_prefix
        self.ufboot = ufboot
        self.boot = boot
        self.input_args = []
        self.input_args_dict = {}

        # Define the program header and command-line arguments
        self.header = f"{'=' * 100}\n\nSatute - asymptotic test for branch saturation\n\n{'=' * 100}\nAuthor:\nCitation:\n{'=' * 100}\n\n"
        self.arguments = [
            {
                "flag": "-o",
                "help": "Output files prefix",
                "default": self.output_prefix,
                "metavar": "<file_name>",
            },
            {
                "flag": "-tree",
                "help": "Path to input tree file",
                "default": self.input_dir,
                "metavar": "<file_name>",
            },
            {
                "flag": "-dir",
                "help": "Path to input directory",
                "default": self.input_dir,
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
                "flag": "-b",
                "help": "Replicates for bootstrap + ML tree + consensus tree",
                "type": int,
                "default": self.boot,
                "metavar": "<num>",
            },
        ]

    def parse_input(self):
        """Parse command-line arguments."""
        parser = argparse.ArgumentParser(description="Satute")
        for argument in self.arguments:
            flag = argument["flag"]
            argument.pop("flag")
            parser.add_argument(flag, **argument)
        self.input_args = parser.parse_args()
        self.input_args_dict = vars(self.input_args)

    def write_log(self):
        """Write the log file."""
        log_lines = []
        log_lines.append(self.header)
        log_lines.append("-" * 100)
        log_lines.append("Satute")
        log_lines.append("-" * 100)
        log_lines.append(f"Command line: {' '.join(sys.argv)}")
        log_lines.append("-" * 100)
        log_lines.append("Printing parameter values:")
        log_lines.append("-" * 100)
        for key, value in vars(self).items():
            log_lines.append(f"{key}:{value}")
        log_lines.append("-" * 100)

        with open(f"{self.input_args_dict['o']}.log", "w") as f:
            f.write("\n".join(log_lines))

    def run(self):
        """
        Main entry point for running the Satute command-line tool.
        """
        # Parsing input arguments and constructing IQ-TREE command-line arguments
        self.parse_input()
        arguments_dict = self.construct_arguments()

        # Running IQ-TREE with constructed arguments
        self.run_iqtree_with_arguments(arguments_dict["arguments"],["--quiet"])

        # If no model specified in input arguments, extract best model from log file
        if not self.input_args.model:
            best_model_log_path = f"{arguments_dict['msa_file']}.log"
            best_model = self.extract_best_model_from_log_file(best_model_log_path)

            if not best_model:
                raise ValueError("Could not find a model in the log file.")

            best_model_name = best_model[0]["Model"]

            logger.info(f"Best model: {best_model_name}")
            logger.info(f"Running a second time with the best model: {best_model_name}")

            # Update model in input arguments and re-construct arguments
            self.input_args.model = best_model_name
            arguments_dict = self.construct_arguments()
            # Running IQ-TREE a second time with updated arguments and --redo option
            self.run_iqtree_with_arguments(arguments_dict["arguments"], ["--quiet"])

        # Here then should come the code, where we should start the saturation test
        # I don't get idea how to do it, so I just call the function from main.py
        # I don't know if it is correct, but it works not
        # saturationTest(self.input_args.dir, self.input_args.o)
        # Something like his should be here:
        # saturationTest(pathDATA, pathIQTREE, runIQTREE = True, runBOOTSRAP = True, dimension = 4, number_rates = 4, chosen_rate = str(4), z_alpha = 2.33, newickformat = 1, epsilon = 0.01, rawMemory = True)

        number_rates = 1
        tree_file_path = self.find_file({".treefile", ".nex", ".nwk"})
        iqtree_file_path = self.find_file({".iqtree"})
        newick_string = ""

        if tree_file_path is not None:
            newick_string = self.get_newick_string(tree_file_path)
        else:
            newick_string = self.get_newick_string_from_iq_tree_file(iqtree_file_path)

        for i in range(number_rates):
            logger.info(f"Here comes the {i+1} th fastest evolving region: ")
            if i == 0:
                saturation_test_cli(
                    str(arguments_dict["msa_file"]),
                    newick_string,
                    self.input_args.iqtree,
                    4,
                    number_rates,
                    str(number_rates - i),
                    2.33,
                    1,
                    0.01,
                    True,
                )
            else:
                saturation_test_cli(
                    str(arguments_dict["msa_file"]),
                    newick_string,
                    self.input_args.iqtree,
                    4,
                    number_rates,
                    str(number_rates - i),
                    2.33,
                    1,
                    0.01,
                    True,
                )

        # End of the code, which should be here

        logger.info(f"Arguments: {arguments_dict}")

        # Writing log file
        self.write_log()

    def run_iqtree_with_arguments(self, arguments, extra_arguments=[]):
        """Run IQ-TREE with given arguments and extra arguments."""
        argument_list = []
        for key, value in self.input_args_dict.items():
            if key != "iqtree" and value and key != "dir":
                argument_list.append(f"-{key} {value}")

        logger.info(f"Running IQ-TREE with the following arguments: {argument_list}")
        extra_arguments_string = " ".join(extra_arguments)

        iq_tree_command = (
            f"{self.input_args.iqtree} {' '.join(arguments)} {extra_arguments_string}"
        )

        logger.info(f"Running IQ-TREE with the following command: {iq_tree_command}")

        try:
            subprocess.run(iq_tree_command, shell=True, check=True)
            logger.info("IQ-TREE execution completed successfully.")
        except subprocess.CalledProcessError as e:
            logger.error(f"IQ-TREE execution failed with the following error: {e}")

    def construct_arguments(self):
        """
        Validate and process input arguments.
        """

        # Convert input paths to Path objects for easier handling
        self.input_args.dir = Path(self.input_args.dir)
        self.input_args.iqtree = Path(self.input_args.iqtree)

        # Check if the input directory exists
        if not self.input_args.dir.is_dir():
            raise InvalidDirectoryError("Input directory does not exist")

        # Define the acceptable file types for sequence alignments and trees
        msa_file_types = {".fasta", ".nex", ".phy"}
        tree_file_types = {".treefile", ".nex", ".nwk"}

        # Find the tree and sequence alignment files in the directory
        tree_file = self.find_file(tree_file_types)
        msa_file = self.find_file(msa_file_types)

        # Check if a sequence alignment file was found
        if msa_file is None:
            raise NoAlignmentFileError("No multiple sequence alignment file found")

        # Initialize the dictionary that will store the argument options
        argument_option = {
            "option": "msa",
            "msa_file": msa_file,
            "arguments": ["-s", str(msa_file), "-asr"],
        }

        # Check if a tree file was found
        if tree_file:
            logger.info(f"{tree_file.suffix} file {tree_file} exists")
            argument_option["option"] = "msa + tree"
            argument_option["arguments"].extend(["-te", str(tree_file)])

        # If a model was specified in the input arguments, add it to the argument options
        if self.input_args.model:
            argument_option["option"] += " + model"
            argument_option["arguments"].extend(["-m", self.input_args.model])

            # If the model includes a Gamma distribution, add the corresponding argument
            if "Gamma" in self.input_args.model:
                argument_option["arguments"].extend(["-wspr"])

            # Todo Check if markov model is reversible if not throw error:

        # Return the constructed argument options
        return argument_option

    def find_file(self, suffixes):
        """Find file in input directory with given suffixes."""
        for file in self.input_args.dir.iterdir():
            if file.suffix in suffixes:
                return file
        return None

    def extract_best_model_from_log_file(self, file_path):
        """Extract best model from IQ-TREE log file."""
        if not self.file_exists(file_path):
            print("File does not exist.")
            return []

        with open(file_path, "r") as file:
            content = file.read()

        table_pattern = re.compile(
            r"No.\s+Model\s+-LnL\s+df\s+AIC\s+AICc\s+BIC\n(.*?\n)---", re.DOTALL
        )
        table_match = re.search(table_pattern, content)

        if table_match:
            table_content = table_match.group(1)
            rows = table_content.strip().split("\n")
            row_data = re.split(r"\s+", rows[0].strip())
            model_info = [
                {
                    "No.": row_data[0],
                    "Model": row_data[1],
                    "-LnL": row_data[2],
                    "df": row_data[3],
                    "AIC": row_data[4],
                    "AICc": row_data[5],
                    "BIC": row_data[6],
                }
            ]

            return model_info
        else:
            print("Table not found in the file.")
            return []

    def file_exists(self, file_path):
        """Check if a file exists."""
        return Path(file_path).exists()

    def get_newick_string_from_iq_tree_file(self, path):
        iq_tree_file = path + ".iqtree"
        newick_string = ""
        with open(iq_tree_file, "r+") as f:
            lines = f.readlines()
            for i in range(0, len(lines)):
                line = lines[i]
                if "Tree in newick format:" in line:
                    newick_string = lines[i + 2]
                    break
        return newick_string

    def get_newick_string(self, file_path):
        """Fetch the Newick string from a file."""
        from pathlib import Path

        # Check if file exists
        if not Path(file_path).is_file():
            raise FileNotFoundError(f"The file at path {file_path} does not exist.")

        with open(file_path, "r") as file:
            newick_string = file.read().strip()

        # Check if file is empty
        if not newick_string:
            raise ValueError(f"The file at path {file_path} is empty.")

        # Check if file contains a valid Newick string
        if not newick_string.endswith(";"):
            raise ValueError(
                f"The file at path {file_path} does not contain a valid Newick string."
            )

        return newick_string


if __name__ == "__main__":
    # Instantiate the Satute class with your desired arguments
    iqtree_path = "iqtree"
    input_directory = Path("./")
    output_prefix = None
    model = None
    num_rate_categories = None
    ufboot_replicates = None
    boot_replicates = None

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
