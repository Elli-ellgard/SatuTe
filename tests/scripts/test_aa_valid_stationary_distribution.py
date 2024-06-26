from tests.scripts.fixtures import *
from satute.amino_acid_models import AA_STATE_FREQUENCIES

@pytest.mark.parametrize("model", list(AA_STATE_FREQUENCIES.keys()))
def test_aa_substitution_models(model):
    from satute.amino_acid_models import normalize_stationary_distribution_aa
    frequencies = {model: AA_STATE_FREQUENCIES[model]}
    frequencies = normalize_stationary_distribution_aa(AA_STATE_FREQUENCIES[model])
    assert sum(frequencies) == pytest.approx(1.0), f"Sum of frequencies is not 1.0 but {sum(frequencies)}"




# def valid_stationary_distribution(frequencies):
#     sum_freqs = sum(frequencies)
#     if sum_freqs == 1:
#         print("valid distribution")
#     else:
#         # Update frequencies dictionary with new values
#         print(sum_freqs)
#         print("NOT a valid distribution")


# if __name__ == "__main__":
#     for model, freqs in AA_STATE_FREQUENCIES.items():
#         print(model)
#         valid_stationary_distribution(freqs)
#         print("")