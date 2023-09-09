# Redefining the TreeNode class, functions, and test data

import copy

# TreeNode class updated to include branch lengths
class TreeNode:
    def __init__(self, name=None, length=None):
        self.name = name
        self.length = length
        self.children = []

    def __repr__(self):
        return f"TreeNode({self.name}, length={self.length}, children={len(self.children)})"

def parse_newick(s):
    stack = []
    current_node = None
    i = 0
    while i < len(s):
        if s[i] == '(':
            new_node = TreeNode()
            if current_node:
                current_node.children.append(new_node)
            stack.append(new_node)
            current_node = new_node
        elif s[i] == ')':
            current_node = stack.pop()
            i += 1
            name = ''
            while i < len(s) and s[i] not in [',', ';', ')', ':']:
                name += s[i]
                i += 1
            if name:
                current_node.name = name.strip()
            if s[i] == ':':  # Extract branch length
                i += 1
                length = ''
                while i < len(s) and s[i] not in [',', ';', ')']:
                    length += s[i]
                    i += 1
                current_node.length = float(length)
            if stack:
                current_node = stack[-1]
        elif s[i] == ',':
            i += 1
            continue
        else:
            name = ''
            while i < len(s) and s[i] not in [',', ';', ')', ':']:
                name += s[i]
                i += 1
            length = None
            if s[i] == ':':  # Extract branch length
                i += 1
                length_str = ''
                while i < len(s) and s[i] not in [',', ';', ')']:
                    length_str += s[i]
                    i += 1
                length = float(length_str)
            child = TreeNode(name.strip(), length)
            current_node.children.append(child)
        i += 1
    return current_node

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

def attach_to_all_edges_with_lengths(unrooted, rooted):
    all_trees = []

    def recursive_attach(node, prev_node=None):
        if prev_node:
            # Store the index of the node in its parent's children list
            index = prev_node.children.index(node)
            
            # Create a new internal node with half the branch length of the current node
            original_length = node.length if node.length is not None else 0.0
            new_internal_node = TreeNode(length=0.5 * original_length if original_length else None)
            
            # Adjust the length of the current node
            node.length = 0.5 * original_length if original_length else None
            
            # Reconnect the nodes
            prev_node.children[index] = new_internal_node
            new_internal_node.children.extend([copy.deepcopy(rooted), node])
            
            # Store the current structure
            all_trees.append(copy.deepcopy(unrooted))
            
            # Revert the changes to continue the recursion
            prev_node.children[index] = node
            node.length = original_length

        # Recurse to child nodes
        for child in node.children:
            recursive_attach(child, node)

    recursive_attach(unrooted)
    return all_trees

# Given Newick strings
unrooted_test_string_check = "((A:1,B:1):1,(C:1,D:1):1);"
rooted_test_string_check = "(((E:1,F:1),G:1):1,I:1);"

# Parse the Newick strings to tree structures
unrooted_tree_check = parse_newick(unrooted_test_string_check)
rooted_tree_check = parse_newick(rooted_test_string_check)

# Attach the rooted tree to every branch of the unrooted tree
attached_trees_check = attach_to_all_edges_with_lengths(unrooted_tree_check, rooted_tree_check)

# Convert the generated trees back to Newick strings
newick_results_check = [to_newick(tree) + ";" for tree in attached_trees_check]

newick_results_check, len(newick_results_check)

print(newick_results_check)
