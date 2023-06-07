# /usr/bin/python3
# -----------------------------------------------------
import sys
import argparse
import numpy as np
import regex as re
import pandas as pd
from ete3 import Tree
import os

# -----------------------------------------------------


def parse_input():
    """Parsing input arguments passed to the script"""
    #################################################
    print("=" * 100)
    print("")
    print("Satute - asymptotic test for branch saturation")
    print("")
    print("=" * 100)
    print("Author: ")
    print("Citation: ")
    print("=" * 100)
    print("")
    ##################################################
    parser = argparse.ArgumentParser(description="Satute")

    # parameters ---------------------------------------
    parser.add_argument(
        "-o",
        help="output files prefix",
        default="satute_output_file",
        metavar="<file_name>",
    )
    parser.add_argument(
        "-dir",
        help="path to input directory",
        default="no directory",
        metavar="<file_name>",
    )
    parser.add_argument(
        "-iqtree", help="path to IQ-TREE", default="iqtree2", metavar="<file_name>"
    )
    parser.add_argument(
        "-model", help="model of evolution", type=str, default="", metavar="<str>"
    )
    parser.add_argument(
        "-nr", help="number of rate categories", type=int, default=1, metavar="<num>"
    )

    parser.add_argument(
        "-ufboot",
        help="Replicates for ultrafast bootstrap (>=1000)",
        type=int,
        default=0,
        metavar="<num>",
    )
    parser.add_argument(
        "-boot",
        help="Replicates for bootstrap + ML tree + consensus tree",
        type=int,
        default=0,
        metavar="<num>",
    )

    # parsing arguments
    args = parser.parse_args(sys.argv[1:])
    # ------------------------------------
    d = vars(args)
    file = open(args.o + ".log", "w")
    file.write("-" * 100 + "\n")
    file.write("Satute" + "\n")
    file.write("-" * 100 + "\n")
    file.write("Command line: ")
    file.write(
        str(sys.argv)
        .replace("'", "")
        .replace(",", "")
        .replace("[", "")
        .replace("]", "")
        + "\n"
    )
    file.write("-" * 100 + "\n")
    file.write("Printing parameter values:" + "\n")
    file.write("-" * 100 + "\n")
    for i in d.keys():
        file.write(str(i) + ":" + str(d[i]) + "\n")
    file.write("-" * 100 + "\n")
    file.close()

    # -------------
    check_input(
        d["iqtree"], d["dir"], d["model"], d["nr"], d["o"], d["ufboot"], d["boot"]
    )


def check_input(iqtree, idir, m, nr, o, ufboot, boot):
    """
    Check if all parameters are correct and which files the input directory includes
    """
    print("Check available input")
    if os.path.isdir(idir):
        # check if it is a possible input msa for iqtree is present
        input_msa = ""
        lst = ["fasta", "nex", "phy"]
        for file in os.listdir(idir):
            for suffix in lst:
                if file.find(suffix) != -1:
                    input_msa = idir + file
                    print(suffix, " file", input_msa, " exists")
        if input_msa == "":
            sys.exit("INPUT_ERROR: no multiple sequence alignment is given")

        # check if treefile is present
        lst = ["treefile", "nex", "nwk"]
        tree_file = ""
        for file in os.listdir(idir):
            for suffix in lst:
                if file.find(suffix) != -1:
                    tree_file = idir + file
                    print(suffix, "file", tree_file, "exists")

        if tree_file != "":
            print("run iqtree with fixed tree")
        else:
            print("run iqtree")

    else:
        sys.exit("INPUT_ERROR: input directory does not exist")


if __name__ == "__main__":
    # parse_input()
    # check_input()
    # satute()

    import logging

    # Configure the logging settings (optional)
    logging.basicConfig(
        level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s"
    )

    # Create a logger
    logger = logging.getLogger(__name__)

    # Log messages with different severity levels
    logger.debug("This is a debug message")
    logger.info("This is an info message")
    logger.warning("This is a warning message")
    logger.error("This is an error message")
    logger.critical("This is a critical message")
