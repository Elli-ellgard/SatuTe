import copy


# TreeNode class
class TreeNode:
    def __init__(self, name=None, length=None):
        self.name = name
        self.length = length
        self.children = []

    def __repr__(self):
        return f"TreeNode({self.name}, length={self.length}, children={len(self.children)})"


# Function to parse Newick string
def parse_newick(s):
    stack = []
    current_node = None
    i = 0
    while i < len(s):
        if s[i] == "(":
            new_node = TreeNode()
            if current_node:
                current_node.children.append(new_node)
            stack.append(new_node)
            current_node = new_node
        elif s[i] == ")":
            current_node = stack.pop()
            i += 1
            name = ""
            while i < len(s) and s[i] not in [",", ";", ")", ":"]:
                name += s[i]
                i += 1
            if name:
                current_node.name = name.strip()
            if s[i] == ":":
                i += 1
                length = ""
                while i < len(s) and s[i] not in [",", ";", ")"]:
                    length += s[i]
                    i += 1
                current_node.length = float(length)
            if stack:
                current_node = stack[-1]
        elif s[i] == ",":
            i += 1
            continue
        else:
            name = ""
            while i < len(s) and s[i] not in [",", ";", ")", ":"]:
                name += s[i]
                i += 1
            length = None
            if s[i] == ":":
                i += 1
                length_str = ""
                while i < len(s) and s[i] not in [",", ";", ")"]:
                    length_str += s[i]
                    i += 1
                length = float(length_str)
            child = TreeNode(name.strip(), length)
            current_node.children.append(child)
        i += 1
    return current_node


# Function to convert tree to Newick string
def to_newick(node):
    if not node.children:
        return f"{node.name}:{node.length}" if node.length is not None else node.name
    children_newick = ",".join(to_newick(child) for child in node.children)
    if node.name and node.length is not None:
        return f"({children_newick}){node.name}:{node.length}"
    elif node.name:
        return f"({children_newick}){node.name}"
    elif node.length is not None:
        return f"({children_newick}):{node.length}"
    return f"({children_newick})"


# Function to attach rooted tree to all branches of unrooted tree
def attach_to_all_edges_with_lengths(unrooted, rooted):
    all_trees = []

    def recursive_attach(node, prev_node=None):
        if prev_node:
            index = prev_node.children.index(node)
            original_length = node.length if node.length is not None else 0.0
            new_internal_node = TreeNode(
                length=0.5 * original_length if original_length else None
            )
            node.length = 0.5 * original_length if original_length else None
            prev_node.children[index] = new_internal_node
            new_internal_node.children.extend([copy.deepcopy(rooted), node])
            all_trees.append(copy.deepcopy(unrooted))
            prev_node.children[index] = node
            node.length = original_length

        for child in node.children:
            recursive_attach(child, node)

    recursive_attach(unrooted)
    return all_trees


# Test examples
examples = [
    ("((A,B),(C,D));", "((E,F),(G,H));"),
    ("((A:1,B:1):1,(C:1,D:1):1);", "((E:1,F:1),(G:1,H:1):1);"),
    ("((A:1,B:1):1,(C:1,D:1):1);", "(((E:1,F:1),G:1):1,I:1);"),
    ("((A,B),(C,D));", "(((E,F),G),H);"),
]

for unrooted_str, rooted_str in examples:
    unrooted_tree = parse_newick(unrooted_str)
    rooted_tree = parse_newick(rooted_str)
    attached_trees = attach_to_all_edges_with_lengths(unrooted_tree, rooted_tree)
    newick_results = [to_newick(tree) + ";" for tree in attached_trees]
    print(
        f"\nFor Unrooted: {unrooted_str} and Rooted: {rooted_str} -> {len(newick_results)} Trees Generated:"
    )
    for tree_str in newick_results:
        print(tree_str)
