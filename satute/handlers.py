# -*- coding: utf-8 -*-
import os
import subprocess
import shutil
from typing import List, Optional
from pathlib import Path
from logging import Logger
from satute.exceptions import IqTreeNotFoundError

class FileHandler:
    """Handles file operations for the Satute project."""

    def __init__(self, base_directory:str):
        """
        Initialize with a base directory.

        Args:
        - base_directory (str): The directory where the handler should operate.
        """
        self.base_directory = base_directory

    def find_file_by_suffix(self, suffixes: List[str])-> Optional[str]:
        """
        Locate a file in the base directory based on its suffixes, ignoring files containing 'satute'.

        Args:
        - suffixes (List[str]): List of file extensions to search for.

        Returns:
        - Optional[str]: Full path to the file if found, None otherwise.
        """
        # Iterate through the base directory and its subdirectories
        for root, dirs, files in os.walk(self.base_directory):
            # Check each file to see if it ends with any of the specified suffixes
            for file in files:
                if "satute" in file:
                    # Skip files containing 'satute'
                    continue
                if any(file.endswith(suffix) for suffix in suffixes):
                    # Return the full path of the file if a match is found
                    return os.path.join(root, file)
        # Return None if no file with the specified suffixes is found
        return None

    def file_exists(self, filepath:str) -> bool:
        """
        Check if a file exists.

        Args:
        - filepath (str): Path to the file to check.

        Returns:
        - bool: True if the file exists, False otherwise.
        """
        return os.path.exists(filepath)

    def read_file(self, filepath)-> list[str]:
        """
        Read the contents of a file.

        Args:
        - filepath (str): Path to the file to read.

        Returns:
        - str: Contents of the file.

        Raises:
        - FileNotFoundError: If the file does not exist.
        - ValueError: If the file is empty.
        """
        # Check if the file exists
        if not self.file_exists(filepath):
            raise FileNotFoundError(f"{filepath} not found.")
        # Check if the file is empty
        if os.path.getsize(filepath) == 0:
            raise ValueError(f"The file at path {filepath} is empty.")
        # Open and read the file, then return its contents
        with open(filepath, "r") as file:
            return file.readlines()

    def get_newick_string_from_iq_tree_file(self, path: str) -> str:
        """
        Extracts Newick format string from an IQ-TREE file.

        Args:
        - path (str): Path to the IQ-TREE file without extension.

        Returns:
        - str: Newick format string.

        Raises:
        - FileNotFoundError: If the IQ-TREE file does not exist.
        - ValueError: If the IQ-TREE file does not contain a valid Newick format string.
        """
        iq_tree_file = Path(f"{path}.iqtree")
        if not iq_tree_file.is_file():
            raise FileNotFoundError(f"IQ-TREE file not found: {iq_tree_file}")

        lines = self.read_file(iq_tree_file)

        # Search for the Newick format header and extract the Newick string
        try:
            newick_line_index = next(
                i for i, line in enumerate(lines) if "Tree in newick format:" in line
            )
            newick_string = lines[newick_line_index + 2].strip()
        except (StopIteration, IndexError):
            raise ValueError(
                f"The IQ-TREE file at {iq_tree_file} does not contain a valid Newick format string."
            )

        if not newick_string.endswith(";"):
            raise ValueError(
                f"The IQ-TREE file at {iq_tree_file} contains an invalid Newick format string."
            )

        return newick_string

    def get_newick_string(self, file_path: str)->str:
        """
        Fetch the Newick string from a file.

        Args:
        - file_path (str): Path to the file containing the Newick string.

        Returns:
        - str: Newick format string.

        Raises:
        - FileNotFoundError: If the file does not exist.
        - ValueError: If the file is empty or does not contain a valid Newick string.
        """
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

    def find_msa_file(self)->str:
        msa_file_types = [".fasta", ".nex", ".phy", ".txt"]
        msa_file = self.find_file_by_suffix(msa_file_types)
        if not msa_file:
            raise FileNotFoundError("No MSA file found in directory")
        return msa_file

    def find_tree_file(self)-> str:
        tree_file_types = [".treefile", ".nex", ".nwk"]
        tree_file = self.find_file_by_suffix(tree_file_types)
        if not tree_file:
            raise FileNotFoundError("No tree file found in directory")
        return tree_file

    def find_iq_tree_file(self)->str:
        iq_tree_file = self.find_file_by_suffix([".iqtree"])
        if not iq_tree_file:
            raise FileNotFoundError("No .iqtree file found in directory")
        return iq_tree_file


class IqTreeHandler:
    """Handles IQ-TREE related operations for the Satute project."""

    def __init__(self, iqtree_path, logger : Logger=None):
        """
        Initialize with an optional base directory.

        Args:
        - base_directory (str, optional): The directory where the handler might operate. Defaults to None.
        """
        self.iq_tree_path = iqtree_path
        self.logger = logger

    def validate_and_append_boot_arguments(self, ufboot, bootstrap):
        """Validates the ufboot and boot parameters and appends them to extra_arguments if valid.

        Raises:
            ValueError: If both ufboot and boot parameters are defined, or if values are less than expected.
        """
        extra_arguments = []  # initialize an empty list for extra_arguments

        # Check if both ufboot and boot parameters are defined
        if ufboot and bootstrap:
            # If both parameters are defined, raise a ValueError
            raise ValueError("Cannot run both ufboot and boot at the same time")
        else:
            # If only ufboot is defined, further check if its value is >= 1000
            if ufboot:
                if ufboot < 1000:
                    # If the value is less than 1000, raise a ValueError
                    raise ValueError("ufboot must be >= 1000")
                # If the value is correct, append it to the list of extra_arguments
                extra_arguments.append(f"--ufboot {ufboot}")
            # If only boot is defined, further check if its value is >= 100
            if bootstrap:
                if bootstrap < 100:
                    # If the value is less than 100, raise a ValueError
                    raise ValueError("boot must be >= 100")
                # If the value is correct, append it to the list of extra_arguments
                extra_arguments.append(f"--boot {bootstrap}")

        return extra_arguments  # return the list of extra_arguments

    def run_iqtree_with_arguments(self, arguments: List[str], extra_arguments: List[str]=[]):
        """
        Run IQ-TREE with given arguments and extra arguments.

        Args:
            arguments (list): List of arguments for IQ-TREE.
            extra_arguments (list): List of additional arguments for IQ-TREE (optional).

        Raises:
            RuntimeError: If IQ-TREE execution fails.
        """
        extra_arguments_string = " ".join(extra_arguments)
        iq_tree_command = f"{self.iq_tree_path} {' '.join(arguments)} {extra_arguments_string}"

        self.logger.info(f"Running IQ-TREE command: {iq_tree_command}")

        try:
            result = subprocess.run(
                iq_tree_command,
                shell=True,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )

            # Log the output and error messages
            self.logger.info(f"IQ-TREE Output: {result.stdout}")
            if result.stderr:
                self.logger.warning(f"IQ-TREE Error: {result.stderr}")

        except subprocess.CalledProcessError as e:
            error_message = (
                f"IQ-TREE execution failed with error code {e.returncode}.\n"
                f"Command: {iq_tree_command}\n"
                f"Output: {e.stdout}\n"
                f"Error: {e.stderr}\n"
                "Please check the command and ensure that IQ-TREE is installed and accessible."
            )
            
            self.logger.error(error_message)
    
            # Raise a RuntimeError with details for further handling
            raise RuntimeError(error_message) from e

    def check_iqtree_path(self, iqtree_path):
        """Check if the given IQ-TREE path exists and raise an exception if it doesn't."""
        if os.path.exists(iqtree_path) and os.path.isfile(iqtree_path) or shutil.which(iqtree_path):
            return True
        else:
            raise IqTreeNotFoundError(f"IQ-TREE does not exist at {iqtree_path}")
