import sys
sys.path.append("./../..")
from amino_acid_models import AA_STATE_FREQUENCIES

def valid_stationary_distribution(frequencies):
    sum_freqs = sum(frequencies)
    if sum_freqs == 1:
        print("valid distribution")
    else:
        # Update frequencies dictionary with new values
        print(sum_freqs)
        print("NOT a valid distribution")


for model, freqs in AA_STATE_FREQUENCIES.items():
    print(model)
    valid_stationary_distribution(freqs)
    print("")