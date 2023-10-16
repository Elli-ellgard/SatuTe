import copy


def generate_newick_from_graph(graph, node, parent=None):
    """
    Recursively generate a Newick string representation of a tree (graph).

    Parameters:
    graph: dict
        Adjacency list representation of a tree.
    node: str
        The current node being processed.
    parent: str, optional
        The parent of the current node being processed.

    Returns:
    str
        Newick string representation of the tree rooted at `node`.
    """
    # Get the children of the current node (excluding the parent to avoid cycles)
    children = [child for child in graph[node] if child != parent]

    # If the node has children, recurse on them and wrap their Newick strings in parentheses
    if children:
        return f"({','.join(generate_newick_from_graph(graph, child, node) for child in children)}){node}"

    # If the node has no children (i.e., it's a leaf), just return its name
    return node


def print_graph(graph):
    """
    Print the adjacency list representation of a graph in a readable format.

    Parameters:
    graph: dict
        The adjacency list representation of the graph.
    """
    for node, neighbors in graph.items():
        print(f"{node}: {', '.join(neighbors)}")
    print("\n")


def add_subtree(graph, edge, subtree):
    graph[edge[0]].remove(edge[1])
    graph[edge[1]].remove(edge[0])
    new_node = f"Node{len(graph) + 1}"
    graph[new_node] = [edge[0], edge[1]]
    graph[edge[0]].append(new_node)
    graph[edge[1]].append(new_node)
    subtree_root = max(subtree, key=lambda k: len(subtree[k]))
    subtree[subtree_root].append(new_node)
    graph[new_node].append(subtree_root)
    graph.update(subtree)
    return graph


def get_edges(graph):
    edges = set()
    for node, adjacency_list in graph.items():
        if len(adjacency_list) != 0:
            for adj_node in adjacency_list:
                # Ensure edges are stored in a consistent order
                edge = tuple([node, adj_node])
                edges.add(edge)
    return edges


def find_node_with_highest_postfix(graph):
    """
    Find the node with the highest numerical postfix in its name.

    Parameters:
        graph (dict): A dictionary representing the graph where keys are nodes and values
                      are lists of nodes that the key node points to.

    Returns:
        str: The name of the node with the highest numerical postfix.
    """
    # Initialize variables to keep track of the node with the highest postfix.
    max_postfix = -1
    max_node = None

    # Loop through each node in the graph.
    for node in graph.keys():
        # Try to extract the numerical postfix from the node name.
        try:
            postfix = int(
                node[4:]
            )  # Assuming the numerical part starts from the 5th character
            if postfix > max_postfix:
                max_postfix = postfix
                max_node = node
        except ValueError:
            # Handle the case where the postfix is not a number (e.g., for nodes 'A', 'B', etc.)
            pass

    return max_node


# Function to parse a node from a Newick string
def parse_node(pair_bracket_string):
    parent_count = 0
    tree = ""
    processed = ""
    index = 0

    # Iterate over each character in the string
    for char in pair_bracket_string:
        # If the character is an opening curved bracket
        if char == "(":
            # Increment the parent count
            parent_count += 1
            # If this is the first opening curved bracket, skip it
            if parent_count == 1:
                continue

        # If the character is a closing curved bracket
        elif char == ")":
            # Decrement the parent count
            parent_count -= 1

            # If there are no more open parents
            if parent_count == 0:
                # If the index is beyond the end of the string
                if index + 2 > len(pair_bracket_string):
                    # Break the loop
                    break
                else:
                    # Otherwise, set the tree to the rest of the string and break the loop
                    tree = pair_bracket_string[index + 2 :]
                    break

        # If the character is a comma
        if char == ",":
            # If there is more than one open parent
            if parent_count != 1:
                # Add a pipe to the processed string
                processed += "|"
            else:
                # Otherwise, add a comma to the processed string
                processed += ","
        else:
            # If the character is not a comma, add it to the processed string
            processed += char
        # Increment the index
        index += 1

    # Split the processed string by semicolons outside brackets
    data = split_on_semicolons_outside_brackets(processed)

    # Replace all pipes with commas in the data
    for i in range(len(data)):
        data[i] = data[i].replace("|", ",")

    # Parse the branch length and label
    label, dist = parse_branch_length_and_label(tree)
    # Return the label, the distance, the data, and the value set
    return (label, dist, data)


def split_on_semicolons_outside_brackets(string):
    result = []
    current = []
    bracket_count = 0

    # Iterate over each character in the string
    for char in string:
        if char == "(":
            bracket_count += 1
        elif char == ")":
            bracket_count -= 1
        elif char == "," and bracket_count == 0:
            result.append("".join(current))
            current = []
            continue
        current.append(char)

    result.append("".join(current))
    return result


# Function to parse the branch length and label from a Newick string
def parse_branch_length_and_label(pair_bracket_string):
    # If the string contains a colon
    if ":" in pair_bracket_string:
        # Split the string into a label and a branch length
        label, branch_length = pair_bracket_string.split(":", maxsplit=1)
        # Convert the branch length to a float
        branch_length = float(branch_length)
        # Return the label and the branch length
        return label, branch_length
    else:
        # If the string does not contain a colon, return the string as the label and an empty string as the branch length
        return pair_bracket_string, ""


def pair_bracket_to_adjacency(
    pair_bracket_string,
    node_count=1,
    adjacency_dictionary=None,
    adjacency_length_map=None,
    rooted=True,
):
    # Remove all semicolons from the string

    pair_bracket_string = pair_bracket_string.replace(";", "")

    # If the string does not contain any opening curved brackets
    if "(" not in pair_bracket_string:
        length = 1
        # If the string contains only one element when split by semicolons outside brackets
        if ":" not in pair_bracket_string:
            # The entire string is the label
            label = pair_bracket_string
            # There is no length
            length = None
        else:
            # The label is the part of the string before the first colon
            label = pair_bracket_string[: pair_bracket_string.find(":")]
            if label.find(":") != -1:
                length = float(pair_bracket_string[pair_bracket_string.find(":") + 1 :])
            # Return a dictionary with the label and the length

            adjacency_dictionary[label] = []
            return {"name": label, "length": length}
    else:
        # If the string contains an opening curved bracket, parse the node
        label, length, data = parse_node(pair_bracket_string)

        if label == "":
            label = f"Node{node_count}*"

        children = [
            pair_bracket_to_adjacency(
                item,
                node_count + 1,
                adjacency_dictionary=adjacency_dictionary,
                adjacency_length_map=adjacency_length_map,
            )
            for item in data
        ]
        adjacency_dictionary[label] = [child["name"] for child in children]

        for child in children:
            adjacency_length_map[(child["name"], label)] = child["length"]

        if not rooted:
            for child in children:
                if child["name"] not in adjacency_dictionary:
                    adjacency_dictionary[child["name"]] = [label]
                else:
                    adjacency_dictionary[child["name"]].append(label)

        # Return a dictionary with the label, the length, the children, and the values
        return {
            "name": label,
            "length": length,
            "children": children,
        }


def find_root(adjacency_dict):
    # Create a set of all nodes that are children
    all_children = set(
        child for children in adjacency_dict.values() for child in children
    )

    # Find a node that is not in the set of all children nodes (i.e., a root)
    for node in adjacency_dict:
        if node not in all_children:
            return node
    return None  # Return None if no root is found (should not happen in a tree)


def delete_and_rewire(node_to_delete, adjacency_dict):
    # Step 1: Identify nodes connected to node_to_delete
    children_of_deleted_node = adjacency_dict.get(node_to_delete, [])
    adjacency_dict.pop(node_to_delete)
    for child_i in children_of_deleted_node:
        adjacency_dict[child_i].remove(node_to_delete)
        for child_j in children_of_deleted_node:
            if child_i != child_j:
                if child_j not in adjacency_dict[child_i]:
                    adjacency_dict[child_i].append(child_j)
                if child_i not in adjacency_dict[child_j]:
                    adjacency_dict[child_j].append(child_i)


def graft_subtree_directed(graph, edge, subtree, subtree_root=None):
    parent, child = edge
    graph[parent].remove(child)  # Remove the old edge
    # Create a new internal node and update the graph    
    new_node = f"grafted_node_{len(graph) + len(subtree) + 1}"
    new_node = f"grafted_node_{id(new_node)}"
    graph[new_node] = [child]
    # Attach the subtree
    graph[new_node].append(subtree_root)  # Update the new node's children
    graph[edge[0]].append(new_node)  # Update the parent's children
    # Update the main graph with the subtree nodes and edges
    graph.update(subtree)
    return graph


def insert_root(newick, insert_text):
    # Find the index of the ")" character
    close_paren_index = newick.rfind(")")
    # Insert the desired text
    modified_newick = (
        newick[: close_paren_index + 1] + insert_text + newick[close_paren_index + 1 :]
    )
    return modified_newick


def to_newick(tree, root):
    def helper(node):
        if not tree[node]:  # If the node has no children (i.e., it's a leaf)
            return node
        else:
            children_str = ",".join(helper(child) for child in tree[node])
            return f"({children_str}){node}"

    newick_str = helper(root)
    return newick_str + ";"


newick_string_main_tree = "(A:1,B:2);"
main_root_node_prefix = "ROOT_NODE_X"
main_root_node = f"{main_root_node_prefix}{id(main_root_node_prefix)}"

newick_string_main_tree = insert_root(newick_string_main_tree, main_root_node)
adjacent_dictionary, adjacent_length_map = {}, {}
tree_structure = pair_bracket_to_adjacency(
    newick_string_main_tree,
    adjacency_dictionary=adjacent_dictionary,
    adjacency_length_map=adjacent_length_map,
    rooted=True,
)

subtree_root_node_prefix = "ROOT_NODE_Y"
subtree_root_node = f"{subtree_root_node_prefix}{id(main_root_node_prefix)}"
newick_subtree = "(E:4,F:5);"
newick_subtree = insert_root(newick_subtree, subtree_root_node)
subtree_adjacent_dictionary, subtree_adjacent_length_map = {}, {}
subtree_structure = pair_bracket_to_adjacency(
    newick_subtree,
    adjacency_dictionary=subtree_adjacent_dictionary,
    adjacency_length_map=subtree_adjacent_length_map,
    rooted=True,
)

edge_set = get_edges(adjacent_dictionary)
newick_list = []

for edge in edge_set:
    #     Make a deep copy of tree_graph and subtree_graph to work with in this iteration
    tree_graph_copy = copy.deepcopy(adjacent_dictionary)
    subtree_graph_copy = copy.deepcopy(subtree_adjacent_dictionary)
    #     Add the subtree to the copied graph, working with the copies in this iteration
    new_tree_graph = graft_subtree_directed(
        tree_graph_copy, edge, subtree_graph_copy, subtree_root_node
    )
    newick_list.append(to_newick(new_tree_graph, main_root_node))
