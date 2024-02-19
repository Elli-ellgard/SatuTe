import dendropy  # Import DendroPy for handling phylogenetic trees
import copy  # Import copy for deep copying complex objects like trees


def insert_taxon_next_to_each_node(tree, taxon_label, target_edge_length=0.000001):
    """
    Inserts a new taxon next to each node in the tree.

    :param tree: The DendroPy Tree object.
    :param taxon_label: The label of the new taxon to insert.
    :param target_edge_length: The length of the edge for the new taxon.
    :return: A list of trees with the new taxon inserted next to each node.
    """
    new_trees = []  # Initialize a list to hold the modified trees
    for (
        node
    ) in tree.levelorder_node_iter():  # Iterate over nodes in level order traversal
        if (
            node.parent_node is None
        ):  # Skip the root node since it doesn't have a parent
            continue
        # Clone the tree to avoid altering the original tree
        cloned_tree = copy.deepcopy(tree)
        # Define a filter to find the corresponding node in the cloned tree
        node_filter = lambda taxon: (
            taxon.label == node.taxon.label if node.taxon else False
        )
        cloned_node = cloned_tree.find_node_with_taxon(node_filter)
        taxon_namespace = (
            cloned_tree.taxon_namespace
        )  # Get the namespace to create new taxa

        if cloned_node:  # Check if the node is found in the cloned tree
            x_taxon = taxon_namespace.new_taxon(
                taxon_label
            )  # Create a new taxon with the given label

            parent_of_cloned_node = (
                cloned_node.parent_node
            )  # Get the parent of the cloned node

            new_node_taxon = taxon_namespace.new_taxon(
                "INSERTED_NODE"
            )  # Create a new taxon for the inserted node

            # Create a new node with the new taxon and specified edge length
            new_node = dendropy.Node(
                taxon=new_node_taxon, edge_length=target_edge_length
            )

            # Add the cloned node and the new taxon as children of the new node
            new_node.add_child(cloned_node)
            new_node.add_child(
                dendropy.Node(taxon=x_taxon, edge_length=target_edge_length)
            )

            # Add the new node to the parent of the cloned node
            parent_of_cloned_node.add_child(new_node)
            # Remove the cloned node from its original parent
            parent_of_cloned_node.remove_child(cloned_node)

            new_trees.append(
                cloned_tree
            )  # Add the modified tree to the list of new trees
    return new_trees  # Return the list of modified trees


if __name__ == "__main__":
    # Define the Newick string for the original tree
    newick_str = "((A:1,B:1)Node1:1,(C:1,D:1)Node2:1);"
    # Load the tree from the Newick string
    tree = dendropy.Tree.get(
        data=newick_str,
        schema="newick",
        rooting="default-rooted",
        suppress_internal_node_taxa=False,
        preserve_underscores=False,
    )
    tree.is_rooted = False  # Specify that the tree is not rooted

    # Insert 'SUBTREE_HERE' next to each node in the tree
    new_taxon_label = "SUBTREE_HERE"
    trees_with_x = insert_taxon_next_to_each_node(tree, new_taxon_label)

    # Create a TreeList object from the list of modified trees
    tree_list = dendropy.TreeList(trees_with_x)

    # Define the path and name of the Nexus file to write
    nexus_file_path = "modified_trees.nex"
    # Write the TreeList to a Nexus file
    tree_list.write(
        path=nexus_file_path,
        schema="nexus",
        suppress_rooting=True,
        unquoted_underscores=True,
    )
    print(f"Modified trees have been written to {nexus_file_path}.")
