#!/usr/bin/env python

import sys
import argparse
import os
import logging
import subprocess
from pathlib import Path

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
        self.check_input()
        self.run_subprocess(f"{self.input_args.iqtree} -s {self.input_args.dir}")
        # Run the desired functionality of your tool
        self.write_log()

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

        input_msa = self.find_file(msa_file_types)
        if input_msa is None:
            logger.error("INPUT_ERROR: no multiple sequence alignment is given")
            sys.exit()
        logger.info(f"{input_msa.suffix} file {input_msa} exists")

        tree_file = self.find_file(tree_file_types)
        if tree_file:
            logger.info(f"{tree_file.suffix} file {tree_file} exists")
            logger.info("run iqtree with fixed tree")
        else:
            logger.info("run iqtree")

    def find_file(self, suffixes):
        for file in self.input_args.dir.iterdir():
            if file.suffix in suffixes:
                return file
        return None

if __name__ == "__main__":
    # Instantiate the Satute class with your desired arguments
    iqtree_path = "iqtree2"
    input_directory = "no directory"
    model = ""
    num_rate_categories = 1
    output_prefix = "satute_output_file"
    ufboot_replicates = 0
    boot_replicates = 0

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
