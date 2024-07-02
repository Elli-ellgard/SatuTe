import numpy as np
from typing import Optional, List
from dataclasses import dataclass

@dataclass
class SubstitutionModel:
    """
    A class representing a substitution model used in phylogenetic or related fields.

    Attributes:
        model (Optional[str]): An optional external model object that may be used for substitution computations.
        state_frequencies (Optional[List[float]]): A list containing the frequencies of different states in the model.
        phi_matrix (Optional[np.ndarray]): A matrix containing phi values, which could represent transition probabilities or other relevant parameters in the model.
        rate_matrix (Optional[np.ndarray]): A matrix containing rate values which could represent substitution rates or other relevant parameters in the model.
        number_rates (Optional[int]): A scalar indicating the number of different rates in the model.
        category_rates (Optional[List[float]]): A list containing the rates of different categories in the model.
        precomputed_q_matrix (Optional[np.ndarray]): An optional precomputed Q matrix, possibly used for faster computations.
    """
    model: Optional[str] = None
    state_frequencies: Optional[List[float]] = None
    phi_matrix: Optional[np.ndarray] = None
    rate_matrix: Optional[np.ndarray] = None
    number_rates: Optional[int] = None
    category_rates: Optional[List[float]] = None
    precomputed_q_matrix: Optional[np.ndarray] = None