# -*- coding: utf-8 -*-
import os
import subprocess
import shutil


class FileHandler:
    """Handles file operations for the Satute project."""

    def __init__(self, base_directory):
        """
        Initialize with a base directory.

        Args:
        - base_directory (str): The directory where the handler should operate.
        """
        self.base_directory = base_directory

    def find_file_by_suffix(self, suffixes):
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

    def file_exists(self, filepath) -> bool:
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

    def write_to_file(self, filepath, content):
        """
        Write content to a file.

        Args:
        - filepath (str): Path to the file to write to.
        - content (str): Content to write to the file.
        """
        # Open and write the content to the file
        with open(filepath, "w") as file:
            file.write(content)

    def get_newick_string_from_iq_tree_file(self, path):
        """
        Extracts Newick format string from an IQ-TREE file.

        Args:
        - path (str): Path to the IQ-TREE file without extension.

        Returns:
        - str: Newick format string.

        Raises:
        - FileNotFoundError: If the IQ-TREE file does not exist.
        - ValueError: If the IQ-TREE file is empty or does not contain a valid Newick format string.
        """
        iq_tree_file = f"{path}.iqtree"
        lines = self.read_file(iq_tree_file)        

        # Instead of iterating through every line, we'll use a more Pythonic approach to find the desired line
        # This should help to slightly optimize the function.
        try:
            # Find the index of the line containing the Newick format header
            newick_line_index = next(
                i for i, line in enumerate(lines) if "Tree in newick format:" in line
            )
            # The Newick string is expected to be two lines after the header
            newick_string = lines[newick_line_index + 2].strip()
        except StopIteration:
            # This exception will be raised if the header line is not found in the file
            raise ValueError(
                f"The IQ-TREE file at path {iq_tree_file} does not contain a valid Newick format string."
            )

        # Additional validation to ensure the Newick string is properly formatted
        if not newick_string.endswith(";"):
            raise ValueError(
                f"The IQ-TREE file at path {iq_tree_file} does not contain a valid Newick format string."
            )

        return newick_string

    def get_newick_string(self, file_path):
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
        - str: Newick format string.
        """
        # If directory is provided, check for tree file and iqtree file
        return self.get_newick_string(
            self.find_file_by_suffix({".treefile", ".nwk"})
        ) or self.get_newick_string_from_iq_tree_file(
            self.find_file_by_suffix({".iqtree"})
        )

    def find_msa_file(self):
        msa_file_types = {".fasta", ".nex", ".phy", ".txt"}
        msa_file = self.find_file_by_suffix(msa_file_types)
        if not msa_file:
            raise FileNotFoundError("No MSA file found in directory")
        return msa_file

    def find_tree_file(self):
        tree_file_types = {".treefile", ".nex", ".nwk"}
        tree_file = self.find_file_by_suffix(tree_file_types)
        if not tree_file:
            raise FileNotFoundError("No tree file found in directory")
        return tree_file

    def find_iqtree_file(self):
        iqtree_file = self.find_file_by_suffix({".iqtree"})
        if not iqtree_file:
            raise FileNotFoundError("No .iqtree file found in directory")
        return iqtree_file


class IqTreeNotFoundError(Exception):
    """Exception raised when IQ-TREE is not found at the given path."""

    pass


class IqTreeHandler:
    """Handles IQ-TREE related operations for the Satute project."""

    def __init__(self, iqtree_path):
        """
        Initialize with an optional base directory.

        Args:
        - base_directory (str, optional): The directory where the handler might operate. Defaults to None.
        """
        self.iqtree_path = iqtree_path

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

    def run_iqtree_with_arguments(self, arguments, extra_arguments=[]):
        """Run IQ-TREE with given arguments and extra arguments."""
        extra_arguments_string = " ".join(extra_arguments)

        iq_tree_command = (
            f"{self.iqtree_path} {' '.join(arguments)} {extra_arguments_string}"
        )

        try:
            subprocess.run(iq_tree_command, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"IQ-TREE execution failed with the following error: {e}")


    def check_iqtree_path(self, iqtree_path):
        """Check if the given IQ-TREE path exists and raise an exception if it doesn't."""
        if os.path.exists(iqtree_path) and os.path.isfile(iqtree_path) or shutil.which(iqtree_path):
            return True
        else:
            raise IqTreeNotFoundError(f"IQ-TREE does not exist at {iqtree_path}")
