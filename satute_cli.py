#!/usr/bin/env python
import sys
import argparse
import os
import logging
from pathlib import Path
import re
import subprocess
from satute_util_new import (
    build_tree_test_space,
    run_saturation_test_for_branches_and_categories,
)
from satute_calculate_likelihoods import run_iqtree_for_each_subtree_parallel
from satute_repository import parse_rate_and_frequencies_and_create_model_files
from ete3 import Tree
from satute_repository import parse_rate_matrices, extract_rate_matrix
from calculate_partial_likelihood import (
    run_calculation_partial_likelihood_for_directories,
)

# Configure the logging settings (optional)
logging.basicConfig(
    level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


from rich import print
from rich.console import Console
from rich.text import Text
from satute_exception import (
    InputArgumentsError,
    InvalidDirectoryError,
    NoAlignmentFileError,
)


def print_ascii_art(ascii_art):
    console = Console()
    console.print(Text(ascii_art, style="bold green"))


# Now you can call the function with your ASCII art
ASCII_ART = """
%%%%%%%%%%%%%%%%/   &%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%(  %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%                                           &%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%   @@@@    @                     @.@   %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%&   @   @,       @/       @   @@   &%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% @@ @      @@@@    @@@(     @  @@#%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%    &(@,@@@      @%@      /@@@@.@    &%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%&    (@@@@@@@@@@    /     @@@@@@@@(*     %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%      @@@@@@@@@@@@  @%  @@@@@@@@@@@@,@.   %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%     %@@@@@@@@@@@   &@@   @@@@@@@@@@@&@     %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%      @@@@@@@@@@ @   @   @ @@@@@@@@@ @      %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%     @(@ @@@ @      @@     @ @@@,  @      %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%  .       @@@. .@@@   ,@@    #@@@@@ @        @  %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%    @ @@      @    @@@@  @@ @@@  .@        @@ @  ,  %%%%%%%%%%%%%
%%%%%%%%%%%%%     @,  @ @ ( * @  @@@@@@@@@@@@@ @ @ @ @ @@ @@        &%%%%%%%%%%%
%%%%%%%%%%%& & @ .   &%@@   @@@ @@@ @*  @ (#*@    @@@/   /.@@@  @  . ,%%%%%%%%%%
%%%%%%%%%%   &   @  @@@@&@@@@@@@@@@@@@@@  (@@@@@%@@%@*&&,.&&@@    @ , %%%%%%%%%%
%%%%%%%%%% # &   ( @&,@@ @ @&(@% & @.@@@@@ %@@ @ @@@@@&@ @@@ ,@  %& &  %%%%%%%%%
%%%%%%%%%    @  @  .@@ @@ @& @ &.& @@ @/@# @ @@*/&&&*@ &*@&# @*     .  %%%%%%%%%
%%%%%%%%%  * @    @@ @#&& &&* @ &%& .@@@@** @,% &@@ @&@&&&&@@@ @       %%%%%%%%%
%%%%%%%%%  ./@,  (@@&%@/&&&/& @ #&& & &/%%%%&&%@ @ @ && && & @#@@      %%%%%%%%%
%%%%%%%%%% @ @  &@@%& @@%@& %&&&.&% &% #%%&%,&& &&@&& &@&@@@*&&@@@/    %%%%%%%%%
%%%%%%%%%%    @@@%#&@  &&%,&##%%%%%.,%.%%%&%.%%% &%%& &&&&  &&&(@@@    %%%%%%%%%
%####%%%%%%    @%%%(%%&%% %%%%%%%(%%%%%%%%%%%%%%&&%%%%%%&/&%%&&& &&    %%%%%%%%%
%%%%%%%%%%%   @&%&% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&@# &     %%%%%%%%%
%%%%%%%%%%%   .% &   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  @ @    %%%%%%%%%
%%%%%%%%%%%%   &.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% & @ @ #,       %%%%%%%%%
%%%%%%%%%%%%   #%/.  @ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  %%%%   @,&  %%%%%%%%%
%%%%%%%%%%%%    , (# *&## %@@&%@@%@%%%%%%%%%%%%%%@@@%@@@%&   %%%%%%   % %%%%%%%%
%%%%%%%%%%%%   *  .%  */@   *   @%%,  %%%%%%% @  @@@%  @@.  %%%%%%%%/ &%%%%%%%%%
%%%%%%%%%%%%%%   .(.  .%%&   *   &%%%  %%%%%%  @ %  %%%%&   .%%%%%%%%%%% %%%%%%%
%%%%%%%%%%%%%% ,##*& %%%%&  %%  %%%%%%%%%%%%&&%&   &%%%%* &%%%%%%%%%%%%%%%%%%%%%
        I will break your tree into parts check if they are saturated
=================================================================================
"""


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


def remove_filename(path):
    parts = path.split("/")
    parts.pop()
    path_folder = "/".join(parts) + "/"  # Path where the we create folders
    return path_folder


def parse_substitution_model(file_path):
    try:
        with open(file_path, "r") as file:
            content = file.read()
            for line in content.splitlines():
                if "Best-fit model according to BIC:" in line:
                    model_string = line.split(":")[1].strip()
                    return model_string
    except (IOError, ValueError):
        # If the file cannot be read or ':' is not found in the content
        # Return None or an appropriate value for error handling
        return None


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
        self.alpha = 0.01

        # Define the program header and command-line arguments
        self.header = f"{'=' * 100}\n\nSatute - asymptotic test for branch saturation\n\n{'=' * 100}\nAuthor:\nCitation:\n{'=' * 100}\n\n"
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
        ]

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
        print_ascii_art(ASCII_ART)

        """Main entry point for running the Satute command-line tool."""
        # Parsing input arguments and constructing IQ-TREE command-line arguments

        self.parse_input()
        arguments_dict = self.construct_arguments()

        # Running IQ-TREE with constructed arguments
        # If no model specified in input arguments, extract best model from log file

        if not self.input_args.model:
            self.run_iqtree_with_arguments(
                arguments_dict["arguments"], ["-m", "TESTONLY", "--redo", "--quiet"]
            )

            substitution_model = parse_substitution_model(
                str(arguments_dict["msa_file"]) + ".iqtree"
            )

            logger.info(
                f"""
                        Best model: {substitution_model}
                        Running a second time with the best model: {substitution_model}"""
            )
            # Update model in input arguments and re-construct arguments
            self.input_args.model = substitution_model
            arguments_dict = self.construct_arguments()

        # =========  Number Rate Handling =========
        number_rates = self.handle_number_rates()
        # =========  End of Number Rate Handling =========

        extra_arguments = arguments_dict.get("model_arguments", []) + [
            "--quiet",
            "--redo",
            # "-T AUTO",
        ]

        # Validate and append ufboot and boot parameters to extra_arguments
        extra_arguments = extra_arguments + self.validate_and_append_boot_arguments()
        # Check if the IQ-TREE state file exists and
        # check in the case of a +G or +R model if the IQ-Tree siteprob file exists
        # if necessary  we have to rerun IQ-TREE
        state_file = self.find_file({".state"})
        site_probability_file = self.find_file({".siteprob"})
        if state_file is None or (number_rates > 1 and site_probability_file is None):
            self.run_iqtree_with_arguments(
                arguments=arguments_dict["arguments"], extra_arguments=extra_arguments
            )
        # ======== Tree File Handling =========
        newick_string = self.get_newick_string_from_args()
        # ======== End Tree File Handling =========

        logger.info(
            f"""
            Running tests and initial IQ-Tree with configurations:
            Saturation Test Model: {self.input_args.model}
            Number of rate categories: {number_rates}
            Options for Initial IQ-Tree run: {arguments_dict['option']}
            Initial Arguments for IQ-Tree: {' '.join(arguments_dict['arguments'])}
            """
        )

        (
            state_frequencies,
            model_and_frequency,
        ) = parse_rate_and_frequencies_and_create_model_files(
            input_path=arguments_dict["msa_file"],
            dimension=4,
            model=self.input_args.model,
        )

        logger.info(
            f"""
                Run with model and frequencies: \n {' '.join(model_and_frequency)}
                Building tree test space with {number_rates} rate categories
            """
        )

        # Build the structure for the subtrees and msas
        valid_category_rates = build_tree_test_space(
            number_rates=number_rates,
            msa_file_name=str(arguments_dict["msa_file"]),
            target_directory=self.active_directory,
        )

        logger.info(
            f"""
                These are the valid {str(valid_category_rates)} rate categories")
                Run IQ-Tree for each subtree in parallel with {number_rates} rate categories
            """
        )

        # TODO make or generate state_space and dimension
        # Extract the rate matrix and row identifiers
        rate_matrix = extract_rate_matrix(str(arguments_dict["msa_file"]) + ".iqtree")
        state_space = rate_matrix.index.tolist()
        dimension = len(state_space)
        state_space = {id: index for index, id in enumerate(state_space)}
        rate_matrix, psi_matrix = parse_rate_matrices(
            dimension, str(arguments_dict["msa_file"])
        )
        
        logger.info(f"Rate Matrix: \n {rate_matrix}")

        # For each rate, run IQ-TREE for the calculation of the posterior distributions in parallel
        for category_rate in valid_category_rates:
            run_calculation_partial_likelihood_for_directories(
                self.active_directory,
                number_rates,
                category_rate,
                rate_matrix,
                dimension,
                state_space,
            )

        logger.info(
            f"""Run test for saturation for each branch and category with {number_rates} rate categories\n 
            Results will be written to the directory:{self.active_directory.name}
            """
        )

        # Run the test for branch saturation
        run_saturation_test_for_branches_and_categories(
            input_directory=str(arguments_dict["msa_file"]),
            target_directory=self.active_directory,
            dimension=dimension,
            number_rates=number_rates,
            t=Tree(newick_string, format=1),
            valid_category_rates=valid_category_rates,
        )

        # Writing log file
        self.write_log(arguments_dict["msa_file"])

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

    def run_iqtree_with_arguments(self, arguments, extra_arguments=[]):
        """Run IQ-TREE with given arguments and extra arguments."""
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
        msa_file_types = {".fasta", ".nex", ".phy", ".txt"}
        tree_file_types = {".treefile", ".nex", ".nwk", ".tree"}

        # Convert input paths to Path objects for easier handling
        if self.input_args.dir:
            self.input_args.dir = Path(self.input_args.dir)
            self.active_directory = self.input_args.dir

        if self.input_args.iqtree:
            self.input_args.iqtree = Path(self.input_args.iqtree)

        if self.input_args.msa:
            self.input_args.msa = Path(self.input_args.msa)
            self.active_directory = self.input_args.msa.parent

        if self.input_args.tree:
            self.input_args.tree = Path(self.input_args.tree)

        argument_option = {}

        if self.input_args.dir:
            # Check if the input directory exists
            if not self.input_args.dir.is_dir():
                raise InvalidDirectoryError("Input directory does not exist")

            if not os.listdir(self.input_args.dir):
                raise InvalidDirectoryError("Input directory is empty")

            # Find iqtree file and extract model
            iqtree_file = self.find_file({".iqtree"})

            if iqtree_file:
                substitution_model = parse_substitution_model(iqtree_file)
                self.input_args.model = substitution_model

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
            if "+G" in self.input_args.model or "+R" in self.input_args.model:
                argument_option["model_arguments"].extend(["-wspr"])

        # Return the constructed argument options

        return argument_option

    def find_file(self, suffixes):
        """Find file in input directory with given suffixes."""
        for file in self.active_directory.iterdir():
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

    def get_newick_string_from_args(self):
        """
        Get the newick string from the provided input arguments.

        Returns:
            newick_string (str): The newick string from the provided file.
        """
        # If directory is provided, check for tree file and iqtree file
        return self.get_newick_string(
            self.find_file({".treefile"})
        ) or self.get_newick_string_from_iq_tree_file(self.find_file({".iqtree"}))

        # If tree argument is provided, get newick string directly

    def validate_and_append_boot_arguments(self):
        """Validates the ufboot and boot parameters and appends them to extra_arguments if valid.

        Raises:
            ValueError: If both ufboot and boot parameters are defined, or if values are less than expected.
        """
        extra_arguments = []  # initialize an empty list for extra_arguments

        # Check if both ufboot and boot parameters are defined
        if self.input_args.ufboot and self.input_args.boot:
            # If both parameters are defined, raise a ValueError
            raise ValueError("Cannot run both ufboot and boot at the same time")
        else:
            # If only ufboot is defined, further check if its value is >= 1000
            if self.input_args.ufboot:
                if self.input_args.ufboot < 1000:
                    # If the value is less than 1000, raise a ValueError
                    raise ValueError("ufboot must be >= 1000")
                # If the value is correct, append it to the list of extra_arguments
                extra_arguments.append(f"--ufboot {self.input_args.ufboot}")
            # If only boot is defined, further check if its value is >= 100
            if self.input_args.boot:
                if self.input_args.boot < 100:
                    # If the value is less than 100, raise a ValueError
                    raise ValueError("boot must be >= 100")
                # If the value is correct, append it to the list of extra_arguments
                extra_arguments.append(f"--boot {self.input_args.boot}")

        return extra_arguments  # return the list of extra_arguments

    def print_args(self):
        print("")
        print("=" * 10)
        print(self.input_args)
        print("=" * 10)
        print("")


# ========= Mapping Code to Newick String ===========
def map_values_to_newick(value_rows, newick_string):
    # Parsing the file of values into a dictionary
    values_dict = {}
    """
     {
                "index": i + 1,
                "delta": delta,
                "c_s": c_s,
                "p_value": p_value,
                "result_test": result_test,
                "result_test_tip2tip": result_test_tip2tip,
                "vector_branches": vector_branches[i],
    }
    """
    for row in value_rows:
        node_one, node_two = row["vector_branches"].split("-")
        print(node_one, node_two)
        # node_two = node_two.split("*")[0]  # remove trailing '*'
        # values_dict[node_two] = {
        #     "delta": row["delta"],
        #     "c_s": row["c_s"],
        #     "p-value": row["p_value"],
        #     "result_test": row["result_test"],
        #     "status": row["result_test_tip2tip"],
        # }

    saturated_newick_string = map_values_to_newick_regex(values_dict, newick_string)

    return saturated_newick_string


def map_values_to_newick_regex(values_dict, newick_string):
    for node_name, values in values_dict.items():
        delta = values["delta"]
        c_s = values["c_s"]
        p_value = values["p-value"]
        branch_status = values["status"]
        newick_string = re.sub(
            rf"({node_name})",
            rf"\1[delta={delta}; c_s={c_s}; p_value={p_value}; branch_status={branch_status}]",
            newick_string,
        )

    return newick_string


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
