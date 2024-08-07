# -*- coding: utf-8 -*-
import pandas as pd
from enum import Enum
import os
import numpy as np
from satute.models.dna_models import NOT_ACCEPTED_DNA_MODELS
from satute.models.substitution_model import SubstitutionModel
from satute.models.dna_models import LIE_DNA_MODELS
from satute.parser.iqtree_aa_parser import AminoAcidParser
from satute.parser.iqtree_nucleotide_parser import NucleotideParser
from satute.parser.base_parser import BaseParser
from typing import List, Optional, Dict
from satute.models.amino_acid_models import (
    AMINO_ACID_MODELS,
    NOT_ACCEPTED_AA_MODELS,
    ESTIMATED_AA_MODELS,
    get_core_model,
)


class ModelType(Enum):
    DNA = "DNA"
    PROTEIN = "Protein"


class IqTreeParser:
    """
    A class to parse IQ-TREE output files and extract relevant information for
    constructing a SubstitutionModel object.
    """

    def __init__(self, file_path: str = None):
        """
        Initializes the IqTreeParser with the path to the IQ-TREE file.

        Parameters:
        - file_path (str, optional): The path to the IQ-TREE file to be parsed.
        """

        self.file_content = []
        self.model_type = ModelType.DNA
        self.file_path = file_path
        if self.file_path:
            self.load_iqtree_file_content()

    def load_iqtree_file_content(self) -> list[str]:
        """
        Loads the content of the .iqtree file into the class attribute 'file_content'.

        This method checks if the file exists and then attempts to read its content.
        If any issues arise during this process, an appropriate exception will be raised.

        Raises:
        - FileNotFoundError: If the specified file does not exist.
        - IOError: If there is an issue reading the file.
        """

        # Check if the file exists
        if not os.path.exists(self.file_path):
            raise FileNotFoundError(f"File '{self.file_path}' not found.")

        try:
            # Attempt to open and read the file content
            with open(self.file_path, "r") as file:
                self.file_content = file.readlines()
        except IOError as e:
            # Raise an IOError if there's a problem reading the file
            raise IOError(
                f"An error occurred while reading the file '{self.file_path}': {str(e)}"
            )

    def check_if_lie_model(self, current_substitution_model: str) -> bool:
        for lie_model in LIE_DNA_MODELS:
            if lie_model in current_substitution_model:
                return True
        return False

    def load_substitution_model(self) -> SubstitutionModel:
        """
        Load and parse the content of the iqtree file to form the substitution model.

        This method parses various components of the iqtree file such as state frequencies,
        rate matrices, model details, number of rate categories, and category rates to
        construct a SubstitutionModel object.

        Returns:
        - SubstitutionModel: An object containing the parsed details of the substitution model.
        """
        current_substitution_model: str = self.parse_substitution_model()
        self.check_model(current_substitution_model)
        self.model_type: ModelType = self.get_model_type(current_substitution_model)

        state_frequencies: List[float] = []
        phi_matrix: Optional[np.ndarray] = []
        rate_matrix: Optional[np.ndarray] = []
        number_rates: int = 1
        category_rates: Optional[Dict[str, Dict[str, float]]] = {}
        pre_computed_q_matrix: np.ndarray = []
        parser: Optional[BaseParser] = None

        if self.model_type == ModelType.DNA:
            parser = NucleotideParser(self.file_content)
            dict_state_frequencies, phi_matrix = parser.parse_state_frequencies()
            state_frequencies = list(dict_state_frequencies.values())
            pre_computed_q_matrix = parser.parse_q_matrix()

            if self.check_if_lie_model(current_substitution_model):
                substitution_rates = parser.parse_substitution_rates()
                rate_matrix = parser.create_rate_matrix_for_lie_markov_models(
                    substitution_rates, state_frequencies
                )
            else:
                # Parse the rate matrix and stationary distribution for the DNA Substitution Model
                state_frequencies = list(dict_state_frequencies.values())
                rate_matrix = parser.construct_rate_matrix(dict_state_frequencies)
        else:
            parser = AminoAcidParser(self.file_content)
            core_model = get_core_model(current_substitution_model)
            if core_model in ESTIMATED_AA_MODELS:
                # In the case when the substitution rates are not estimated so we get the static elements
                substitution_rates: str = parser.parse_substitution_rates()
                state_frequencies, phi_matrix = parser.parse_state_frequencies()
                state_frequencies = (
                    AminoAcidParser.normalize_stationary_distribution_aa(
                        state_frequencies
                    )
                )
                rate_matrix = AminoAcidParser.create_rate_matrix_with_input(
                    20, substitution_rates, state_frequencies
                )
            else:
                # In the case when the substitution rates are not estimated so we get the static elements
                current_substitution_model = current_substitution_model.upper()
                # Parse the rate matrix and stationary distribution for the Protein Substitution Model
                state_frequencies, phi_matrix = (
                    parser.get_aa_state_frequency_substitution_models(
                        current_substitution_model
                    )
                )
                state_frequencies = (
                    AminoAcidParser.normalize_stationary_distribution_aa(
                        state_frequencies
                    )
                )
                rate_matrix = AminoAcidParser.get_aa_rate_matrix(
                    current_substitution_model
                )

        number_rates = self.parse_number_rate_categories()
        category_rates = self.parse_category_rates() if number_rates > 1 else None
        gamma_shape = self.parse_gamma_shape()

        return SubstitutionModel(
            model=current_substitution_model,
            state_frequencies=state_frequencies,
            phi_matrix=phi_matrix,
            rate_matrix=rate_matrix,
            number_rates=number_rates,
            category_rates=category_rates,
            precomputed_q_matrix=pre_computed_q_matrix,
            gamma_shape=gamma_shape,
        )

    def check_model(self, model: str) -> None:
        """
        Check if the model is one of the not accepted DNA or Protein models.

        Parameters:
        - model (str): The model string to be checked.

        Raises:
        - ValueError: If the model is not accepted for analysis.
        """
        for dna_model in NOT_ACCEPTED_DNA_MODELS:
            if dna_model in model:
                raise ValueError(
                    f"The DNA model '{dna_model}' is not accepted for analysis because it is non-reversible."
                )

        for aa_model in NOT_ACCEPTED_AA_MODELS:
            if aa_model in model:
                raise ValueError(
                    f"The protein model '{aa_model}' is not accepted for analysis."
                )

    def get_model_type(self, model: str) -> ModelType:
        # Check if it is a protein model
        model_upper = model.upper()
        for amino_acid_substitution_model in AMINO_ACID_MODELS + ESTIMATED_AA_MODELS:
            if amino_acid_substitution_model in model_upper:
                return ModelType.PROTEIN

        # Default to DNA if no conditions above are met
        return ModelType.DNA

    def parse_number_rate_categories(self) -> int:
        """
        Parse the number of rate categories from the substitution model string.

        Returns:
        - rate (int): The number of rate categories parsed from the model string.
        """
        index = next(
            (
                idx
                for idx, line in enumerate(self.file_content)
                if "Model of rate heterogeneity:" in line
            ),
            None,
        )

        if index is None:
            raise ValueError("'Model of rate heterogeneity:' not found in file.")

        line = self.file_content[index]

        if "Uniform" in line:
            return 1
        if "with" in line and "categories" in line:
            start = line.index("with ") + len("with ")
            end = line.index(" categories")
            return int(line[start:end])
        if "Invar" in line:
            return 1
        raise ValueError("Unexpected format for 'Model of rate heterogeneity:' line.")

    def parse_category_rates(self) -> Dict[str, Dict[str, float]]:
        """
        Parses the category rates from the file content and returns them in a structured format.

        The function identifies the table of category rates in the file content and extracts
        the category, relative rate, and proportion for each row in the table.

        Returns:
        - dict: A dictionary containing category rates. Each key represents a category and the
                associated value is another dictionary with details of that category.
        """

        # Get the number of rate categories from another method
        number_rates = self.parse_number_rate_categories()

        # Dictionary to store parsed table data
        table_data = {}

        # Variables to track the start and end of the table in the file content
        start_index = -1
        end_index = -1

        # Find the start and end indices of the table in a single pass
        for i, line in enumerate(self.file_content):
            stripped_line = line.strip()
            if stripped_line.startswith("Category"):
                start_index = i + 1
            elif start_index != -1 and stripped_line.startswith(f"{number_rates}"):
                end_index = i + 1
                break

        # Error handling in case the table isn't found
        if start_index == -1 or end_index == -1:
            raise ValueError("Table not found in the log file.")

        # Extract and parse the table rows
        for line in self.file_content[start_index:end_index]:
            parts = line.split()
            if len(parts) >= 3:
                category, relative_rate, proportion = parts[:3]
                if category != "0":
                    table_data[f"p{category}"] = {
                        "Category": category,
                        "Relative_rate": float(relative_rate),
                        "Proportion": float(proportion),
                    }

        return table_data

    def parse_substitution_model(self) -> str:
        """
        Parses the file content to extract the substitution model.

        The function searches for lines that contain specific keywords indicative of the substitution model.
        Once found, it extracts and returns the model as a string.

        Returns:
        - str: The extracted substitution model. If not found, raises a ValueError.
        """
        # Keywords indicating the presence of the substitution model in a line
        keywords = ["Best-fit model according to BIC:", "Model of substitution:"]
        for line in self.file_content:
            # Check if the line contains any of the keywords
            if any(keyword in line for keyword in keywords):
                model_string = line.split(":")[1].strip()
                # If a valid model string is found, return it
                if model_string:
                    return model_string
        # If the loop completes without returning, raise an error
        raise ValueError("Substitution model not found in the file content.")

    def parse_gamma_shape(self) -> float:
        """
        Parses the file content to extract the substitution model.

        The function searches for lines that contain specific keywords indicative of the substitution model.
        Once found, it extracts and returns the model as a string.

        Returns:
        - str: The extracted substitution model. If not found, raises a ValueError.
        """
        # Keywords indicating the presence of the substitution model in a line
        keywords = ["Gamma shape alpha:"]

        for line in self.file_content:
            # Check if the line contains any of the keywords
            if any(keyword in line for keyword in keywords):
                alpha = line.split(":")[1].strip()
                # If a valid model string is found, return it
                if alpha:
                    return float(alpha)
        # If the loop completes without returning, raise an error
        return None


def parse_substitution_model(file_path: str) -> str:
    """
    Parse the substitution model from an IQ-TREE log file.

    This function reads an IQ-TREE log file and extracts the substitution model
    based on specific lines containing "Best-fit model according to BIC:" or
    "Model of substitution:". The function returns the extracted model as a string.
    """
    try:
        with open(file_path, "r") as file:
            content = file.read()
            for line in content.splitlines():
                if ("Best-fit model according to BIC:" in line) or (
                    "Model of substitution:" in line
                ):
                    model_string = line.split(":")[1].strip()
                    if model_string:
                        return model_string
            raise ValueError("Expected model strings not found in the file.")
    except IOError:
        raise ValueError("Could not read the file.")


def parse_rate_from_cli_input(model: str) -> int:
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
        # Extract the substring after e.g.'+G'
        number = model[rate_start_index:]

        # Parse the extracted substring as an integer
        if "{" in number:
            # e.g. +G{0.9} will fix the Gamma shape parameter (alpha)to 0.9
            # discrete Gamma model: default 4 rate categories
            rate = 4
            return rate
        else:
            if number and str(number).isnumeric():
                # number of rate categories
                rate = int(number)
                return rate
            else:
                return "AMBIGUOUS"
    except ValueError:
        # If '+G' is not found or the number after '+G' is not a valid integer
        # Return None or an appropriate value for error handling
        raise ValueError("Could not parse the substitution model from the file.")


def parse_file_to_data_frame(file_path: str) -> pd.DataFrame:
    try:
        # Read the file into a dataframe
        df = pd.read_csv(file_path, delimiter="\t")
        return df

    except FileNotFoundError:
        raise Exception(f"File not found: {file_path}")
