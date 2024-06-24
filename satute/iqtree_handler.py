# -*- coding: utf-8 -*-
import os
import shutil
import subprocess
import sys 
from pathlib import Path
from logging import Logger
from typing import List, Optional
from satute.exceptions import IqTreeNotFoundError

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

    def validate_and_append_boot_arguments(self, ufboot, bootstrap) -> List[str]:
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
                self.logger.error(f"IQ-TREE Error: {result.stderr}")

        except subprocess.CalledProcessError as e:
            error_message = (
                f"IQ-TREE execution failed with error code {e.returncode}.\n"
                f"Command: {iq_tree_command}\n"
                f"Output: {e.stdout}\n"
                f"Error: {e.stderr}\n"
                "Please check the command and ensure that IQ-TREE is installed and accessible."
            )
            
            self.logger.error(error_message)            
            raise RuntimeError(error_message) from e

    def check_iqtree_path(self, iqtree_path):
        """Check if the given IQ-TREE path exists and raise an exception if it doesn't."""
        if os.path.exists(iqtree_path) and os.path.isfile(iqtree_path) or shutil.which(iqtree_path):
            return True
        else:
            raise IqTreeNotFoundError(f"IQ-TREE does not exist at {iqtree_path}")
