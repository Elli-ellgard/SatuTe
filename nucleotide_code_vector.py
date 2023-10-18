import numpy as np

NUCLEOTIDE_CODE_VECTOR = {
    "A": np.array([1, 0, 0, 0]),
    "C": np.array([0, 1, 0, 0]),
    "G": np.array([0, 0, 1, 0]),
    "T": np.array([0, 0, 0, 1]),
    "U": np.array([0, 0, 0, 1]),
    "R": np.array([1, 0, 1, 0]),
    "Y": np.array([0, 1, 0, 1]),
    "K": np.array([0, 0, 1, 1]),
    "M": np.array([1, 1, 0, 0]),
    "S": np.array([0, 1, 1, 0]),
    "W": np.array([1, 0, 0, 1]),
    "B": np.array([0, 1, 1, 1]),
    "D": np.array([1, 0, 1, 1]),
    "H": np.array([1, 1, 0, 1]),
    "V": np.array([1, 1, 1, 0]),
    #The following keys are treated as State_Unknown in IQ-Tree
    "N": np.array([1, 1, 1, 1]),
    "-": np.array([1, 1, 1, 1]),
    "?": np.array([1, 1, 1, 1]),
    ".": np.array([1, 1, 1, 1]),
    "~": np.array([1, 1, 1, 1]),
    "X": np.array([1, 1, 1, 1]),
    #Additional key from EvoNaps database
    "!": np.array([1, 1, 1, 1]),
}