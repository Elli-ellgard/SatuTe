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


class IqTreeNotFoundError(Exception):
    """Exception raised when IQ-TREE is not found at the given path."""

    pass


class ModelNotFoundError(Exception):
    """Exception raised when the amino acid substitution model is not found."""

    def __init__(self, model_name: str):
        self.model_name = model_name
        self.message = f"Model '{model_name}' not found in AMINO_ACID_RATE_MATRIX."
        super().__init__(self.message)


class InvalidModelNameError(Exception):
    """Exception raised when the core model name cannot be extracted."""

    def __init__(self, model_name: str):
        self.model_name = model_name
        self.message = f"Could not extract core model from: {model_name}"
        super().__init__(self.message)
