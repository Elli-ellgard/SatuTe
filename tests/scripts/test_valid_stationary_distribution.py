import sys
sys.path.append("./../..")
from satute_repository import valid_stationary_distribution


# Example usage:
frequencies = {'A': 0.2, 'B': 0.3, 'C': 0.5}
print("Input frequencies: ", frequencies)
print("sum of frequnecies:", sum(frequencies.values()))
frequencies = valid_stationary_distribution(frequencies)
print("Updated frequencies:", frequencies)
print("sum of frequnecies:", sum(frequencies.values()))

frequencies = {'A': 0.2, 'B': 0.3, 'C': 0.24, 'D': 0.25}
print("Input frequencies: ", frequencies)
print("sum of frequnecies:", sum(frequencies.values()))
frequencies = valid_stationary_distribution(frequencies)
print("Updated frequencies:", frequencies)
print("sum of frequnecies:", sum(frequencies.values()))



