# -*- coding: utf-8 -*-

class InputArgumentsError(Exception):
    """
    Exception raised for errors in the input arguments.

    Attributes:
        message -- explanation of the error
    """

    def __init__(
        self,
        message="Both 'msa' and 'dir' input arguments are defined. Please decide between 'msa' or 'dir' input.",
    ):
        self.message = message
        super().__init__(self.message)


class InvalidDirectoryError(Exception):
    """Exception raised when the input directory does not exist."""

    pass


class NoAlignmentFileError(Exception):
    """Exception raised when no multiple sequence alignment file is found."""

    pass
