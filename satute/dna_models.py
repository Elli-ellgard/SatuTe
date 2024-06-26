# -*- coding: utf-8 -*-

import numpy as np



DNA_MODELS = [
    "JC", # Equal substitution rates and equal base frequencies (Jukes and Cantor, 1969).
    "JC69",
    "F81", # Equal rates but unequal base freq. (Felsenstein, 1981).	
    "K80", # Unequal transition/transversion rates and equal base freq. (Kimura, 1980).
    "K2P",
    "HKY", # Unequal transition/transversion rates and unequal base freq. (Hasegawa, Kishino and Yano, 1985).
    "HKY85",
    "TN", # Like HKY but unequal purine/pyrimidine rates (Tamura and Nei, 1993).
    "TN93",
    "TNe", # Like TN but equal base freq.
    "K81", # Three substitution types model and equal base freq. (Kimura, 1981).
    "K3P", 	
    "K81u", # Like K81 but unequal base freq.
    "TPM2", # AC=AT, AG=CT, CG=GT and equal base freq.	
    "TPM2u", # Like TPM2 but unequal base freq.
    "TPM3", # AC=CG, AG=CT, AT=GT and equal base freq.	
    "TPM3u", # Like TPM3 but unequal base freq.	
    "TIM", # Transition model, AC=GT, AT=CG and unequal base freq.	
    "TIMe", # Like TIM but equal base freq.	
    "TIM2", # AC=AT, CG=GT and unequal base freq.	
    "TIM2e", # Like TIM2 but equal base freq.	
    "TIM3", # AC=CG, AT=GT and unequal base freq.
    "TIM3e", #	Like TIM3 but equal base freq.
    "TVM", # Transversion model, AG=CT and unequal base freq.
    "TVMe", # Like TVM but equal base freq.
    "SYM", # Symmetric model with unequal rates but equal base freq. (Zharkikh, 1994).
    "GTR", # General time reversible model with unequal rates and unequal base freq. (Tavare, 1986).
    # Substitution Model Code
    "000000", # JC, JC69, F81
    "010010", # K80, K2P, HKY, HKY85 
    "010020", # TN, TN93, TNe
    "012210", # K81, K3P, K81u
    "010212", # TPM2, TPM2u
    "012012", # TPM3, TPM3u
    "012230", # TIM, TIMe
    "010232", # TIM2, TIM2e
    "012032", # TIM3, TIM3e
    "012314", # TVM, TVMe
    "012345", # SYM, GTR
    #Lie Markov Models
    "1.1",
    "2.2b",
    "3.3a", 
    "3.3c",
    "3.4",
    "4.4a",
    "4.4b",
]

NOT_ACCEPTED_DNA_MODELS = [
    #non-reversible Lie Markov Models
    "3.3b",
    "4.5a",
    "4.5b", 
    "5.6a",
    "5.6b", 
    "5.7a", 
    "5.7b",
    "5.7c",
    "5.11a",
    "5.11b",
    "5.11c",
    "5.16",
    "6.6",	#equiv. to STRSYM for strand-symmetric model (Bielawski and Gold, 2002)
    "6.7a", #F81+K3P
    "6.7b",
    "6.8a", 
    "6.8b",	
    "6.17a",
    "6.17b", 
    "8.8",
    "8.10a",
    "8.10b", 
    "8.16",
    "8.17",
    "8.18",
    "9.20a",
    "9.20b",
    "10.12",
    "10.34", 
    "12.12", #equiv. to UNREST (unrestricted model)
]

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