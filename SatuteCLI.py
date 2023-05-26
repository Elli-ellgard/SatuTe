#!/usr/bin/env python

import sys
import argparse
import os
import logging
import subprocess
from pathlib import Path
import re


# Configure the logging settings (optional)
logging.basicConfig(
    level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


class Satute:
    def __init__(self, iqtree, input_dir, model, nr, output_prefix, ufboot, boot):
        self.iqtree = iqtree
        self.input_dir = Path(input_dir)  # Convert to Path object
        self.model = model
        self.nr = nr
        self.output_prefix = output_prefix
        self.ufboot = ufboot
        self.boot = boot
        self.input_args = []
        self.input_args_dict = {}
        self.header = f"{'=' * 100}\n\nSatute - asymptotic test for branch saturation\n\n{'=' * 100}\nAuthor:\nCitation:\n{'=' * 100}\n\n"
        self.arguments = [
            {
                "flag": "-o",
                "help": "Output files prefix",
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
                "flag": "-iqtree",
                "help": "Path to IQ-TREE",
                "default": self.iqtree,
                "metavar": "<file_name>",
            },
            {
                "flag": "-m",
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
        ]

    def parse_input(self):
        parser = argparse.ArgumentParser(description="Satute")
        for argument in self.arguments:
            flag = argument["flag"]
            argument.pop("flag")
            parser.add_argument(flag, **argument)
        self.input_args = parser.parse_args()
        self.input_args_dict = vars(self.input_args)

    def write_log(self):
        with open(f"{self.input_args_dict['o']}.log", "w") as f:
            f.write(self.header)
            f.write("-" * 100 + "\n")
            f.write("Satute\n")
            f.write("-" * 100 + "\n")
            f.write(f"Command line: {' '.join(sys.argv)}\n")
            f.write("-" * 100 + "\n")
            f.write("Printing parameter values:\n")
            f.write("-" * 100 + "\n")
            for key, value in vars(self).items():
                f.write(f"{key}:{value}\n")
            f.write("-" * 100 + "\n")

    def run(self):
        # Run the Satute command-line tool.
        self.parse_input()
        file_arguments = self.check_input()
        self.run_iqtree_with_arguments(file_arguments)
        self.write_log()

    def run_iqtree_with_arguments(self, file_argument):
        argument_list = []
        for key, value in self.input_args_dict.items():
            if key != "iqtree" and value and key != "dir":
                argument_list.append(f"-{key} {value}")
        logger = logging.getLogger(__name__)
        logger.info(f"Running IQ-TREE with the following arguments: {argument_list}")
        iq_tree_command = f"{self.input_args.iqtree} {file_argument['argument']} {file_argument['file']} {' '.join(argument_list)}"
        logger.info(f"Running IQ-TREE with the following command: {iq_tree_command}")
        try:
            subprocess.run(iq_tree_command, shell=True, check=True)
            logger.info("IQ-TREE execution completed successfully.")
        except subprocess.CalledProcessError as e:
            logger.error(f"IQ-TREE execution failed with the following error: {e}")    

    def run_subprocess(self, command):
        try:
            process = subprocess.Popen(
                command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True
            )
            stdout, stderr = process.communicate()

            if process.returncode != 0:
                logging.error(f"Error occurred: {stderr.decode('utf-8')}")
            else:
                logging.info(f"Output: {stdout.decode('utf-8')}")
        except Exception as e:
            logging.error(f"An error occurred while running the subprocess: {str(e)}")

    def check_input(self):
        logger.info("Check available input")
        self.input_args.dir = Path(self.input_args.dir)
        self.input_args.iqtree = Path(self.input_args.iqtree)

        if not self.input_args.dir.is_dir():
            logger.error("INPUT_ERROR: input directory does not exist")
            sys.exit()

        msa_file_types = {".fasta", ".nex", ".phy"}
        tree_file_types = {".treefile", ".nex", ".nwk"}
            
        tree_file = self.find_file(tree_file_types)
        msa_file = self.find_file(msa_file_types)

        if msa_file is None:
            logger.error("INPUT_ERROR: no multiple sequence alignment is given")
            sys.exit()
        if tree_file:
            logger.info("Run iqtree with fixed tree")
            logger.info(f"{tree_file.suffix} file {tree_file} exists")
            return {"option": "fixed_tree", "file": tree_file, 'argument': '-t'}
        else:
            logger.info("Run iqtree")
            logger.info(f"{msa_file.suffix} file {msa_file} exists")
            return {"option": "msa", "file": msa_file, 'argument': "-s"}

    def find_file(self, suffixes):
        for file in self.input_args.dir.iterdir():
            if file.suffix in suffixes:
                return file
        return None

    def extract_best_model_from_log_file(self, file_path):
        if not self.file_exists(file_path):
            print("File does not exist.")
            return []

        with open(file_path, 'r') as file:
            content = file.read()

        table_pattern = re.compile(r'No.\s+Model\s+-LnL\s+df\s+AIC\s+AICc\s+BIC\n(.*?\n)---', re.DOTALL)
        table_match = re.search(table_pattern, content)

        if table_match:
            table_content = table_match.group(1)
            rows = table_content.strip().split('\n')
            row_data = re.split(r'\s+', rows[0].strip())
            model_info = []
            model_info.append({
                'No.': row_data[0],
                'Model': row_data[1],
                '-LnL': row_data[2],
                'df': row_data[3],
                'AIC': row_data[4],
                'AICc': row_data[5],
                'BIC': row_data[6]
            })

            return model_info
        else:
            print("Table not found in the file.")
            return []

    def file_exists(self, file_path):
        try:
            with open(file_path):
                pass
            return True
        except FileNotFoundError:
            return False

if __name__ == "__main__":
    # Instantiate the Satute class with your desired arguments
    iqtree_path = "iqtree2"
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
    # Example usage
    file_path = './test/octo-kraken-test/example.phy.log'
    model_info = satute.extract_best_model_from_log_file(file_path)



 

