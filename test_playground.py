from scipy.sparse.linalg import expm
import numpy as np

RATE_MATRIX = np.array([[-3, 1, 1, 1], [1, -3, 1, 1], [1, 1, -3, 1], [1, 1, 1, -3]])

sp5_pl = [0, 0, 1, 0]
sp1_pl = [1, 0, 0, 0]


sp2_pl = [0, 1, 0, 0]
sp7_pl = [1, 0, 0, 0]

transition_matrix = expm(RATE_MATRIX * 0.1)
inner_node_1_partial_likelihood = (transition_matrix @ sp5_pl) * (
    transition_matrix @ sp1_pl
)

print("Inner Node Likelihood 1", inner_node_1_partial_likelihood)
print(inner_node_1_partial_likelihood)


inner_node_2_partial_likelihood = (transition_matrix @ sp2_pl) * (
    transition_matrix @ sp7_pl
)


print("Inner Node Likelihood 2", inner_node_1_partial_likelihood)
root_partial_likelihood = (transition_matrix @ inner_node_1_partial_likelihood) * (
    transition_matrix @ inner_node_1_partial_likelihood
)


print("Root Likelihood", root_partial_likelihood)
