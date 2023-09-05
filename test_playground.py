from scipy.sparse.linalg import expm
import math
from satute_test_statistic_using_partial_likelihood import (get_transition_matrix,)
import numpy as np

# special example: JC model
RATE_MATRIX = np.array([[-3, 1, 1, 1], [1, -3, 1, 1], [1, 1, -3, 1], [1, 1, 1, -3]])

def JC_transition_matrix(lambd, branch_length):
    p = (1 + 3* math.exp(-lambd*4*branch_length))/4
    q = (1 - math.exp(-lambd*4*branch_length))/4
    return np.array([[p, q, q, q], [q, p, q, q], [q, q, p, q], [q, q, q, p]])
           
           

# example tree 
sp5_pl = [0, 0, 1, 0]
sp5_bl = 0.1
sp1_pl = [1, 0, 0, 0]
sp1_bl = 0.1

sp2_pl = [0, 1, 0, 0]
sp2_bl= 0.1
sp7_pl = [1, 0, 0, 0]
sp7_bl = 0.1

N1_N3_bl= 0.1
N2_N3_bl= 0.1

sp3_pl = [1, 0, 0, 0]
sp3_bl = 0.1
sp4_pl = [1, 0, 0, 0]
sp4_bl = 0.1

N5_N4_bl= 0.1

sp6_pl = [0, 0, 1, 0]
sp6_bl = 0.1

branch_length = 0.1


# calculation

transition_matrix = get_transition_matrix(RATE_MATRIX, branch_length)
print(transition_matrix)
transition_matrix_2 = JC_transition_matrix(1, branch_length)
print(transition_matrix_2)


#caclualtion by hand
print("Calculate likelihood of N1")
factor_sp5=get_transition_matrix(RATE_MATRIX,sp5_bl)@np.array(sp5_pl)
factor_sp1=get_transition_matrix(RATE_MATRIX,sp1_bl)@np.array(sp1_pl)
N1_pl=factor_sp5*factor_sp1
print(N1_pl)

# inner_node_1_partial_likelihood = (transition_matrix @ sp5_pl) * (
#     transition_matrix @ sp1_pl
# )

# print("Inner Node Likelihood 1", inner_node_1_partial_likelihood)
# print(inner_node_1_partial_likelihood)

print("Calculate likelihood of N2")
factor_sp2=get_transition_matrix(RATE_MATRIX,sp2_bl)@np.array(sp2_pl)
factor_sp7=get_transition_matrix(RATE_MATRIX,sp7_bl)@np.array(sp7_pl)
N2_pl=factor_sp2*factor_sp7
print(N2_pl)

# inner_node_2_partial_likelihood = (transition_matrix @ sp2_pl) * (
#     transition_matrix @ sp7_pl
# )

# print("Inner Node Likelihood 2", inner_node_1_partial_likelihood)


print("Calculate likelihood of N3")
factor_N1=get_transition_matrix(RATE_MATRIX,N1_N3_bl) @ N1_pl
factor_N2=get_transition_matrix(RATE_MATRIX,N2_N3_bl) @ N2_pl
N3_pl=factor_N1*factor_N2


# root_partial_likelihood = (transition_matrix @ inner_node_1_partial_likelihood) * (
#     transition_matrix @ inner_node_1_partial_likelihood
# )


# print("Root Likelihood", root_partial_likelihood)
