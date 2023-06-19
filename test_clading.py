from ete3 import Tree
import os

def name_nodes_by_level_order(tree):
    i = 1
    for node in tree.traverse("levelorder"):
        if not node.is_leaf():
            node.name = f"node_{i}"
            i += 1
    return tree


def get_all_subtrees(tree):
    subtrees = []
    for node in tree.traverse("levelorder"):
        if not node.is_root():
            subtrees.append(node)
    return subtrees


def get_opposite_subtree(tree, subtree):
    # Create a copy of the tree so as not to modify the original
    tree_copy = Tree(tree.write())
    tree_copy = name_nodes_by_level_order(tree_copy)
    # Get the node in the copied tree that corresponds to the root of the subtree
    node = tree_copy.search_nodes(name=subtree.name)[0]
    # Detach this node (and its descendants) from the tree
    node.detach()
    return tree_copy


def generate_subtree_pairs(subtrees, t):
    subtree_pairs = []
    for subtree in subtrees:
        subtree_copy = t.search_nodes(name=subtree.name)[0]
        opposite_subtree = get_opposite_subtree(t, subtree_copy)  
        subtree_pairs.append((subtree.write(format=1), opposite_subtree.write(format=1)))

    return subtree_pairs


t = Tree("((A:1,(B:1,(E:1,D:1):0.5):0.5),(A1:1,(B1:1,(E1:1,D1:1):0.5):0.5));")
t = name_nodes_by_level_order(t)
subtrees = get_all_subtrees(t)
subtrees = subtrees[1:]
generated_subtree_pairs = generate_subtree_pairs(subtrees, t)

for i in range(len(generated_subtree_pairs)):
    os.makedirs(f"./clades/Branch/{i}_clade1/", exist_ok=True)
    os.makedirs(f"./clades/Branch/{i}_clade2/", exist_ok=True)  
    first_clade_writer = open(f"./clades/Branch/{i}_clade1/subtree.treefile", "w")
    second_clade_writer = open(f"./clades/Branch/{i}_clade2/subtree.treefile", "w")
    first_clade_writer.write(generated_subtree_pairs[i][0])
    second_clade_writer.write(generated_subtree_pairs[i][1])






