import numpy
from ete3 import Tree


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

def generate_subtree_pairs(subtrees,t ):
    subtree_pairs = []
    for subtree in subtrees:
        subtree_copy = t.search_nodes(name=subtree.name)[-1]
        opposite_subtree = get_opposite_subtree(t, subtree_copy)
        subtree_pairs.append((subtree.write(), opposite_subtree.write()))
        print(subtree, opposite_subtree)
    return subtree_pairs

def generate_clades(t):
    t = name_nodes_by_level_order(t)
    subtrees = get_all_subtrees(t) 
    subtrees = subtrees[1:] # exclude the first element, otherwise the root pair is present twice
    subtrees_pairs = generate_subtree_pairs(subtrees, t)
    for pair in subtrees_pairs:
        print("*"*20)
        print(pair[0])
        print("="*10)
        print(pair[1])
    return subtrees_pairs



#t = "(((((t110000:0.00269526,t110001:0.00124194):0.00145376,(t110010:0.00179943,t110011:0.00207133):0.00161743):0.00239711,((t110100:0.000897087,t110101:0.00183948):0.00279912,(t110110:0.00127609,t110111:0.00221269):0.0012904):0.00203528):0.00314861,((t111000:0.00226518,t111001:0.0010435):0.00208953,(t111010:0.00288125,t111011:0.00141789):0.00204724):0.00321358))Root:1;"
t = "(LngfishAu:0.1713680662,(LngfishSA:0.1887831949,LngfishAf:0.1651645456)Node2:0.1075382326,(Frog:0.2569055294,((((Turtle:0.2219745817,(Crocodile:0.3064867197,Bird:0.2316135443)Node8:0.0651886769)Node7:0.0365628305,Sphenodon:0.3455137640)Node6:0.0205129559,Lizard:0.3869325486)Node5:0.0741361926,(((Human:0.1854113942,(Seal:0.0945454145,(Cow:0.0824091856,Whale:0.1013725423)Node13:0.0404916927)Node12:0.0252768166)Node11:0.0341299963,(Mouse:0.0584577774,Rat:0.0906477493)Node14:0.1219957375)Node10:0.0608333563,(Platypus:0.1923158391,Opossum:0.1511978659)Node15:0.0373256833)Node9:0.1494087207)Node4:0.1277507551)Node3:0.0942698255)Node1;"
T= Tree(t, format=1)
T = name_nodes_by_level_order(T)
#print(t)
print(T.write())

subtrees_pairs = generate_clades(T)

