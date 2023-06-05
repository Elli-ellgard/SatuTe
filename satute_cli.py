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


def parse_rate_and_frequencies_alternative_and_create_model_files(
    path, number_rates, dimension, model="GTR"
):
    """
    Note: not my function
    Parse the rate parameter and state frequencies from the IQ-TREE log file.
    """

    with open(path + ".iqtree", "r") as f:
        found = 0

        number_lines = 0

        modelString = ""

        for line in f:
            if "Rate parameter R:" in line:
                found = 1
            if found and number_lines < dimension * (dimension - 1) // 2 + 1:
                if number_lines > 1:
                    modelString += line[7:]
                number_lines += 1
        separated = modelString.splitlines()

        modelFinal = model + "{"
        modelFinal += separated[0]

        for words in separated[1:]:
            modelFinal += "," + words

        modelFinal += "}"

    with open(path + ".iqtree", "r") as f:
        # FIX: Use split instead of character position, too delicate otherwise

        found = 0

        number_lines = 0

        frequencyString = ""

        for line in f:
            if "State frequencies:" in line:
                found = 1
            if found and number_lines < dimension + 2:
                if "equal frequencies" in line:
                    frequencyString = "0.25\n" * dimension
                    break

                if number_lines > 1:
                    frequencyString += line[10:]

                number_lines += 1

        separated = frequencyString.splitlines()

        frequencyFinal = "FU{"

        frequencyFinal += separated[0]

        for words in separated[1:]:
            frequencyFinal += "," + words
        frequencyFinal += "}"
        modelAndFrequency = modelFinal + "+" + frequencyFinal

    state_frequencies_vector = []

    # Convert string to float
    for state_frequency in separated:
        state_frequencies_vector.append(float(state_frequency))

    path_folder = remove_filename(path)

    if number_rates == 1:
        f1 = open(path_folder + "model.txt", "w")
        f1.write(modelAndFrequency)
    else:
        for i in range(number_rates):
            f1 = open(
                path_folder + "subsequences/subseq" + str(i + 1) + "/model.txt", "w"
            )
            f1.write(modelAndFrequency)

    return state_frequencies_vector


def parse_rate_from_model(model):
    try:
        # Find the index of '+G' in the model string
        plus_g_index = model.index("+G")

        # Extract the substring after '+G'
        number = model[plus_g_index + 2 :]

        # Parse the extracted substring as an integer
        rate = int(number)

        return rate
    except ValueError:
        # If '+G' is not found or the number after '+G' is not a valid integer
        # Return None or an appropriate value for error handling
        logger.info("Could not find a Gamma distribution in the model.")
        logger.info("So not rates were found.")
        return 1


def remove_filename(path):
    parts = path.split("/")
    parts.pop()
    path_folder = "/".join(parts) + "/"  # Path where the we create folders
    return path_folder


def parse_substitution_model(file_path):
    try:
        with open(file_path, "r") as file:
            content = file.read()
            # Find the index of ':' in the file content
            colon_index = content.index(":")

            # Extract the substring after ':'
            model_string = content[colon_index + 1 :].strip()

            return model_string
    except (IOError, ValueError):
        # If the file cannot be read or ':' is not found in the content
        # Return None or an appropriate value for error handling
        return None


class Satute:
    """Class representing Satute command-line tool for wrapping up functions of IQ-TREE."""

    def __init__(self, iqtree, input_dir, model, nr, output_prefix, ufboot, boot):
        self.iqtree = iqtree
        self.input_dir = None
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
                "flag": "-outdir",
                "help": "Output File Directory",
                "default": self.output_prefix,
                "metavar": "<file_name>",
            },
            {
                "flag": "-dir",
                "help": "Path to input directory",
                "default": self.input_dir,
                "metavar": "<file_name>",
            },
            {
                "flag": "-tree",
                "help": "Path to input tree file",
                "default": self.input_dir,
                "metavar": "<file_name>",
            },
            {
                "flag": "-msa",
                "help": "Path to MSA",
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

    def write_log(self, msa_file):
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

        if self.input_args.dir:
            with open(f"{str(msa_file)}-satute.log", "w") as f:
                f.write("\n".join(log_lines))

    def run(self):
        """Main entry point for running the Satute command-line tool."""
        # Parsing input arguments and constructing IQ-TREE command-line arguments
        self.parse_input()
        arguments_dict = self.construct_arguments()

        # Running IQ-TREE with constructed arguments
        # If no model specified in input arguments, extract best model from log file
        if not self.input_args.model:
            # best_model_log_path = f"{arguments_dict['msa_file']}.log"
            # best_model = # This is wrong # parse_rate_and_frequencies_alternative(best_model_log_path)
            # best_model_name = best_model[0]["Model"]
            # self.run_iqtree_with_arguments([ar,"-m TESTONLY"], ["--redo", "--quiet"])

            self.run_iqtree_with_arguments(
                arguments_dict["arguments"], ["-m", "TESTONLY", "--redo", "--quiet"]
            )
            substitution_model = parse_substitution_model(
                str(arguments_dict["msa_file"]) + ".iqtree"
            )
            logger.info(f"Best model: {substitution_model}")
            logger.info(
                f"Running a second time with the best model: {substitution_model}"
            )
            # Update model in input arguments and re-construct arguments
            self.input_args.model = substitution_model
            arguments_dict = self.construct_arguments()

        # =========  Number Rate Handling =========
        number_rates = 1
        # If no number of rate categories specified in input arguments, extract it from the model
        if self.input_args.nr:
            number_rates = self.input_args.nr
        else:
            if self.input_args.model:
                number_rates = parse_rate_from_model(self.input_args.model)
        # =========  End of Number Rate Handling =========

        # =========  IF Tree File  =========
        if self.input_args.dir:
            # CHECK IF Tree File IS DEFINED
            tree_file_path = self.find_file({".treefile", ".nex", ".nwk"})
            # CHECK IF NOT SEE IF IQ Tree File IS DEFINED
            iqtree_file_path = self.find_file({".iqtree"})

            # FIND Tree is file and get Newick String
            if tree_file_path is not None:
                newick_string = self.get_newick_string(tree_file_path)
            else:
                newick_string = self.get_newick_string_from_iq_tree_file(
                    iqtree_file_path
                )
        # =========  IF parameter -tree TREE IS DEFINED =========
        elif self.input_args.tree:
            # CHECK IF Tree File IS DEFINED
            tree_file_path = self.input_args.tree
            # GET Newick String
            newick_string = self.get_newick_string(tree_file_path)
        # ======== End Tree File Handling =========
        logger.info(f"Run Initial IQ-Tree with options {arguments_dict['option']}")

        extra_arguments = arguments_dict.get("model_arguments", []) + [
            "--quiet",
            "--redo",
        ]

        self.run_iqtree_with_arguments(
            arguments=arguments_dict["arguments"], extra_arguments=extra_arguments
        )

        state_frequencies_vector = (
            parse_rate_and_frequencies_alternative_and_create_model_files(
                str(arguments_dict["msa_file"]), number_rates, 4
            )
        )

        logger.info(f"Run Saturation Test with Model {self.input_args.model}")
        logger.info(f"Run Saturation Test with {number_rates} rate categories")
        logger.info(f"Run Saturation Test with Rates {state_frequencies_vector}")

        for i in range(number_rates):
            logger.info(f"Here comes the {i+1} th fastest evolving region: ")
            if i == 0:
                saturation_test_cli(
                    str(arguments_dict["msa_file"]),
                    newick_string,
                    str(self.input_args.iqtree),
                    4,
                    number_rates,
                    str(number_rates - i),
                    2.33,
                    1,
                    0.01,
                    True,
                    state_frequencies_vector
                )
            else:
                saturation_test_cli(
                    str(arguments_dict["msa_file"]),
                    newick_string,
                    str(self.input_args.iqtree),
                    4,
                    number_rates,
                    str(number_rates - i),
                    2.33,
                    1,
                    0.01,
                    True,
                    state_frequencies_vector,
                )

        # End of the code, which should be here
        logger.info(f"Arguments: {arguments_dict}")
        # Writing log file
        self.write_log(arguments_dict["msa_file"])

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

        Raises:
            InvalidDirectoryError: If the input directory does not exist.
            NoAlignmentFileError: If no multiple sequence alignment file is found.

        Returns:
            A dictionary with keys 'option', 'msa_file', and 'arguments' that represents the argument options for the process.
        """
        # Define the acceptable file types for sequence alignments and trees
        msa_file_types = {".fasta", ".nex", ".phy"}
        tree_file_types = {".treefile", ".nex", ".nwk"}

        # Convert input paths to Path objects for easier handling
        if self.input_args.dir:
            self.input_args.dir = Path(self.input_args.dir)

        if self.input_args.iqtree:
            self.input_args.iqtree = Path(self.input_args.iqtree)

        if self.input_args.msa:
            self.input_args.msa = Path(self.input_args.msa)
            self.input_args.dir = self.input_args.msa.parent

        if self.input_args.tree:
            self.input_args.tree = Path(self.input_args.tree)

        argument_option = {}
        if self.input_args.dir:
            # Check if the input directory exists
            if not self.input_args.dir.is_dir():
                raise InvalidDirectoryError("Input directory does not exist")

            # Find the tree and sequence alignment files in the directory
            tree_file = self.find_file(tree_file_types)
            msa_file = self.find_file(msa_file_types)

            # Check if a sequence alignment file was found
            if msa_file is None:
                raise NoAlignmentFileError("No multiple sequence alignment file found")

            argument_option = {
                "option": "msa",
                "msa_file": msa_file,
                "arguments": ["-s", str(msa_file), "-asr"],
            }

            # Check if a tree file was found
            if tree_file:
                print(f"{tree_file.suffix} file {tree_file} exists")
                argument_option["option"] = "msa + tree"
                argument_option["arguments"].extend(["-te", str(tree_file)])

        else:
            if self.input_args.msa:
                argument_option = {
                    "option": "msa",
                    "msa_file": self.input_args.msa,
                    "arguments": ["-s", str(self.input_args.msa), "-asr"],
                }

            if self.input_args.tree:
                argument_option["option"] = "msa + tree"
                argument_option["arguments"].extend(["-te", str(self.input_args.tree)])

        # If a model was specified in the input arguments, add it to the argument options
        if self.input_args.model:
            argument_option["option"] += " + model"
            argument_option["model_arguments"] = ["-m", self.input_args.model]
            # If the model includes a Gamma distribution, add the corresponding argument
            if "+G" in self.input_args.model:
                argument_option["model_arguments"].extend(["-wspr"])
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
    input_directory = None  # Path("./")
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
