import numpy as np

class SubstitutionModel:
    """
    A class representing a substitution model used in phylogenetics or related fields.

    Attributes:
    model (optional): An external model object that may be used for substitution computations.
    state_frequencies (optional): A data structure (e.g., list or array) containing the frequencies of different states in the model.
    phi_matrix (optional): A matrix containing phi values, which could represent transition probabilities or other relevant parameters in the model.
    rate_matrix (optional): A matrix containing rate values which could represent substitution rates or other relevant parameters in the model.
    number_rates (optional): A scalar indicating the number of different rates in the model.
    category_rates (optional): A data structure (e.g., list or array) containing the rates of different categories in the model.
    """

    def __init__(
        self,
        model : str = None,  # An optional external model object
        state_frequencies: list[float] = None,  # Optional state frequencies
        phi_matrix : np.array = None,  # Optional phi matrix
        rate_matrix : np.array = None,  # Optional rate matrix
        number_rates: int = None,  # Optional number of rates
        category_rates: int = None,  # Optional category rates
        precomputed_q_matrix: np.array = None, # Optional
    ) -> None:
        """Initializes a new SubstitutionModel instance with the provided parameters."""
        self.model = model  # Set the model attribute
        self.number_rates = number_rates  # Set the number_rates attribute
        self.rate_matrix = rate_matrix  # Set the rate_matrix attribute
        self.state_frequencies = (
            state_frequencies  # Set the state_frequencies attribute
        )
        self.phi_matrix = phi_matrix  # Set the phi_matrix attribute
        self.category_rates = category_rates  # Set the category_rates attribute
        self.precomputed_q_matrix = precomputed_q_matrix # Set the precomputed