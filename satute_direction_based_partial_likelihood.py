import numpy as np
from typing import Optional
from functools import cache
from scipy.sparse.linalg import expm
import dataclasses
from ete3 import Tree
from satute_rate_categories_and_alignments import read_alignment_file
from multiprocessing import Pool, cpu_count

NUCLEOTIDE_CODE_VECTOR = {
    "A": [1, 0, 0, 0],
    "C": [0, 1, 0, 0],
    "G": [0, 0, 1, 0],
    "T": [0, 0, 0, 1],
    "N": [1, 1, 1, 1],
    "-": [1, 1, 1, 1],
}

RATE_MATRIX = np.array([[-3, 1, 1, 1], [1, -3, 1, 1], [1, 1, -3, 1], [1, 1, 1, -3]])


class RateMatrix:
    def __init__(self, rate_matrix):
        self.rate_matrix = rate_matrix

    def __hash__(self):
        return id(self)


class DirectedAcyclicGraph:
    def __init__(self, edges):
        for left, right, branch_length in edges:
            left.connect(right, branch_length)
            right.connect(left, branch_length)

        self.edges = edges

    def get_edges(self):
        return self.edges

    def get_branch_length(self, edge):
        return self.branch_lengths[edge]


class Node:
    name: str
    state: Optional[np.array] = dataclasses.field(default_factory=lambda: None)
    connected: dict = dataclasses.field(default_factory=dict)

    def __init__(self, name, state=None, connected=None):
        self.name = name
        self.state = state
        self.connected = connected or {}

    def __hash__(self):
        return id(self)

    def __eq__(self, other):
        return self is other

    def connect(self, other_node, branch_length):
        self.connected[other_node] = branch_length

    def is_leaf(self):
        return self.state is not None


def get_initial_likelihood_vector(state):
    return np.array(NUCLEOTIDE_CODE_VECTOR[state]).T


def name_nodes_by_level_order(tree):
    i = 1
    for node in tree.traverse("levelorder"):
        if not node.is_leaf():
            node.name = f"Node{i}"
            i += 1
    return tree


def convert_ete3_tree_to_directed_acyclic_graph(
    tree, msa_column, alignment_look_up_table
):
    node_dictionary = {}

    for node in tree.traverse("levelorder"):
        if node.is_leaf():
            node_dictionary[node.name] = Node(
                node.name,
                get_initial_likelihood_vector(
                    msa_column[alignment_look_up_table[node.name]].seq
                ),
            )

        else:
            node_dictionary[node.name] = Node(node.name)

    edge_list = []
    for node in tree.traverse("levelorder"):
        for child_node in node.children:
            edge_list.append(
                (
                    node_dictionary[node.name],
                    node_dictionary[child_node.name],
                    child_node.dist,
                )
            )

    return DirectedAcyclicGraph(edge_list)


def get_alignment_look_up_table(alignment):
    alignment_look_up_table = {}
    i = 0
    for record in alignment:
        alignment_look_up_table[record.id] = i
        i = i + 1
    return alignment_look_up_table


#@profile
def calculate_partial_likelihoods_for_sites(tree, alignment, rate_matrix):
    alignment_look_up_table = get_alignment_look_up_table(alignment)
    partial_likelihood_per_site_storage = {}

    for i in range(1, len(alignment[0].seq), 1):
        graph = convert_ete3_tree_to_directed_acyclic_graph(
            tree, alignment[:, (i + 1) - 1 : (i + 1)], alignment_look_up_table
        )

        for edge in graph.get_edges():
            left, right, branch_length = edge

            p1 = partial_likelihood(graph, left, right, rate_matrix)
            p2 = partial_likelihood(graph, right, left, rate_matrix)
            
            if (
                f"({left.name}, {right.name})"
                not in partial_likelihood_per_site_storage
            ):
                partial_likelihood_per_site_storage[f"({left.name}, {right.name})"] = {
                    "left": [
                        {
                            "Node": left.name,
                            "Site": i,
                            "p_A": p1[0],
                            "p_C": p1[1],
                            "p_G": p1[2],
                            "p_T": p1[3],
                            "branch_length": branch_length,
                        }
                    ],
                    "right": [
                        {
                            "Node": right.name,
                            "Site": i,
                            "p_A": p2[0],
                            "p_C": p2[2],
                            "p_G": p2[2],
                            "p_T": p2[3],
                            "branch_length": branch_length,
                        }
                    ],
                }
            else:
                partial_likelihood_per_site_storage[f"({left.name}, {right.name})"][
                    "left"
                ].append(
                    {
                        "Node": left.name,
                        "Site": i,
                        "p_A": p1[0],
                        "p_C": p1[1],
                        "p_G": p1[2],
                        "p_T": p1[3],
                        "branch_length": branch_length,
                    }
                )
                partial_likelihood_per_site_storage[f"({left.name}, {right.name})"][
                    "right"
                ].append(
                    {
                        "Node": right.name,
                        "Site": i,
                        "p_A": p1[0],
                        "p_C": p1[1],
                        "p_G": p1[2],
                        "p_T": p1[3],
                        "branch_length": branch_length,
                    }
                )
    return partial_likelihood_per_site_storage


@cache
#@profile
def partial_likelihood(tree, node, coming_from, rate_matrix):
    results = 1
    if node.is_leaf():
        return node.state
    for child in node.connected.keys():
        if not child.__eq__(coming_from):

            e = calculate_exponential_matrix(rate_matrix, node.connected[child])
            p = partial_likelihood(tree, child, node, rate_matrix)

            results = results * (e @ p)

    return results


@cache
def calculate_exponential_matrix(rate_matrix, branch_length):
    return expm(rate_matrix.rate_matrix * branch_length)


def test_one_partial_likelihood():
    a1, b2, u3, u4, c5, b6 = [
        Node("A", get_initial_likelihood_vector("A")),
        Node("B", get_initial_likelihood_vector("C")),
        Node("3"),
        Node("4"),
        Node("C", get_initial_likelihood_vector("A")),
        Node("D", get_initial_likelihood_vector("A")),
    ]

    tree = DirectedAcyclicGraph(
        [(a1, u3, 0.01), (u3, b2, 0.01), (u3, u4, 0.01), (u4, c5, 0.01), (u4, b6, 0.01)]
    )

    for edge in tree.get_edges():
        left, right, branch_length = edge
        p1 = partial_likelihood(tree, left, right)
        p2 = partial_likelihood(tree, right, left)
        test(p1, p2)


#@profile
def main():
    alignment_file = "./Clemens/example_3/sim-JC+G-AC1-AG1-AT1-CG1-CT1-GT1-alpha1.2-taxa64-len1000bp-bla0.01-blb0.8-blc0.2-rep01.fasta"
    newick_string = "(t000000:0.0067892257,t000001:0.0067787574,(((t000010:0.0059673375,((((((t010000:0.0030563760,t010001:0.0114093850):0.0094182376,(t010010:0.0131446398,t010011:0.0097395536):0.0133621566):0.0104108104,(((t011000:0.0124958339,t011001:0.0125011474):0.0075917776,(t011010:0.0051495665,t011011:0.0124730837):0.0069341413):0.0078912673,((t011100:0.0104118838,t011101:0.0072227008):0.0073307877,(t011110:0.0144514225,t011111:0.0055786880):0.0102245774):0.0059268453):0.0184390656):0.0061628417,(t010110:0.0082616199,t010111:0.0146332627):0.0101990740):0.0052544609,(t010100:0.0092991497,t010101:0.0114804795):0.0092491114):0.6992169228,(((((t100000:0.0114103154,t100001:0.0051613214):0.0061697358,(t100010:0.0136268730,t100011:0.0114221532):0.0072470758):0.0155632848,((t100100:0.0081218605,t100101:0.0073837399):0.0130036860,(t100110:0.0071915054,t100111:0.0125350507):0.0131619893):0.0090043100):0.0007243217,(((t101000:0.0106452037,t101001:0.0182293876):0.0045033981,(t101010:0.0121747946,t101011:0.0121664259):0.0105869819):0.0092355725,((t101100:0.0078255361,t101101:0.0176067822):0.0088049086,(t101110:0.0061914059,t101111:0.0114288359):0.0119823991):0.0077233513):0.0115530529):0.7027436378,(((((t110000:0.0152109080,t110001:0.0070318894):0.0081739309,(t110010:0.0101504845,t110011:0.0117088941):0.0091799218):0.0134467581,((t110100:0.0050989755,t110101:0.0103895458):0.0158105135,(t110110:0.0072259592,t110111:0.0124929711):0.0072567017):0.0115489293):0.0176840292,((t111000:0.0127615540,t111001:0.0059373002):0.0117927352,(t111010:0.0162749028,t111011:0.0080293329):0.0114795491):0.0181492686):0.0106254068,((t111100:0.0082902032,t111101:0.0114231618):0.0065691489,(t111110:0.0101567116,t111111:0.0058159250):0.0047846372):0.0005525044):0.8384907952):0.4960190049):0.7890571203):0.0055032299,t000011:0.0125092154):0.0125564734,(((t000100:0.0129052942,t000101:0.0081677003):0.0128058021,(t000110:0.0061739221,t000111:0.0124881756):0.0088603366):0.0105251504,(((t001000:0.0124981985,t001001:0.0103935537):0.0089686048,(t001010:0.0135912957,t001011:0.0093498130):0.0044762602):0.0161461211,((t001100:0.0093389059,t001101:0.0125044394):0.0064576918,(t001110:0.0104262512,t001111:0.0135413564):0.0135893360):0.0090709753):0.0287882859):0.0040229738):0.0188329058);"
    t = Tree(newick_string, format=1)
    t = name_nodes_by_level_order(t)
    alignment = read_alignment_file(alignment_file)

    rate_matrix = RateMatrix(RATE_MATRIX)

    partial_likelihood_per_site_storage = calculate_partial_likelihoods_for_sites(
        t, alignment, rate_matrix#
    )


if __name__ == "__main__":
    main()
