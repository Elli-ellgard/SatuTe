import unittest
import sys

sys.path.append("./../../")

from ete3 import Tree
from typing import List

# Assuming rescale_branch_lengths is located in a module named mymodule
from satute_trees import (
    rename_internal_nodes_pre_order,
    rescale_branch_lengths,
    has_duplicate_leaf_sequences,
    collapse_identical_leaf_sequences,
    check_and_raise_for_duplicate_nodes,
    print_tree_with_inner_node_names,
    check_all_internal_nodes_annotated,
)

# Assuming all necessary functions (rescale_branch_lengths, has_duplicate_leaf_sequences, collapse_identical_leaf_sequences, and delete_children_nodes) are defined in the same module or imported correctly.
class TestPhylogeneticTreeFunctions(unittest.TestCase):

    def test_rescale_branch_lengths(self):
        test_tree = Tree("((A:1,B:2)Node2:0.5,C:3);", format=1)
        rescale_factor = 2
        expected_lengths = {"A": 2, "B": 4, "C": 6, "Node2": 1}

        rescale_branch_lengths(test_tree, rescale_factor)

        for node in test_tree.traverse("postorder"):
            if node.is_leaf() or node.name in expected_lengths:
                self.assertAlmostEqual(
                    node.dist,
                    expected_lengths.get(node.name, 0),
                    msg=f"Failed for node {node.name}",
                )

    def test_no_duplicates(self):
        tree = Tree("((A,B),(C,D));", format=1)
        msa = {"A": "ATCG", "B": "ATCGT", "C": "AGGG", "D": "AGGGC"}
        self.assertFalse(has_duplicate_leaf_sequences(tree, msa))

    def test_with_duplicates(self):
        tree = Tree("((A,B),(C,D));")
        msa = {"A": "ATCG", "B": "ATCG", "C": "AGGG", "D": "AGGG"}
        self.assertTrue(has_duplicate_leaf_sequences(tree, msa))

    def test_collapse_identical_leaf_sequences(self):
        tree = Tree("((A:1,B:1)Node1:0.5,(C:1,D:1)Node2:0.5)Root;", format=1)
        msa = {"A": "ATCG", "B": "ATCG", "C": "GGCC", "D": "GGCC"}

        updated_msa, collapsed_nodes = collapse_identical_leaf_sequences(tree, msa)

        expected_collapsed_nodes = {"Node1": ["A", "B"], "Node2": ["C", "D"]}
        expected_msa = {
            "A": "ATCG",
            "B": "ATCG",
            "C": "GGCC",
            "D": "GGCC",
            "Node1": "ATCG",
            "Node2": "GGCC",
        }

        self.assertEqual(updated_msa, expected_msa)
        self.assertEqual(collapsed_nodes, expected_collapsed_nodes)

        for node_name in expected_collapsed_nodes.keys():
            node = tree.search_nodes(name=node_name)[0]
            self.assertEqual(
                len(node.children),
                0,
                f"{node_name} should have no children after collapsing.",
            )

    def test_no_collapse_needed(self):
        tree = Tree("((A:1,B:1)Node1:0.5,(C:1,D:1)Node2:0.5)Root;", format=1)
        msa = {
            "A": "ATCG",
            "B": "ATGC",  # Different from A
            "C": "GGCC",
            "D": "GCCG",  # Different from C
        }
        updated_msa, collapsed_nodes = collapse_identical_leaf_sequences(tree, msa)

        expected_collapsed_nodes = {}
        expected_msa = (
            msa  # Unchanged because no sequences are identical under the same node
        )

        self.assertEqual(
            updated_msa,
            expected_msa,
            "MSA should remain unchanged when no collapsing is needed.",
        )
        self.assertEqual(
            collapsed_nodes,
            expected_collapsed_nodes,
            "No nodes should be collapsed when all leaf sequences are unique.",
        )
        self.assertTrue(
            all(
                len(node.children) == 2
                for node in tree.get_descendants()
                if not node.is_leaf()
            ),
            "Tree structure should remain unchanged when no collapsing is needed.",
        )

    def test_duplicate_node_names_raises_exception(self):
        # Tree with duplicate node names
        tree_with_duplicates = Tree("((A,B)C,(D,E)C);", format=1)

        # Expect check_and_raise_for_duplicate_nodes to raise ValueError due to duplicate node names
        with self.assertRaises(ValueError) as context:
            check_and_raise_for_duplicate_nodes(tree_with_duplicates)

        # Optionally, check if the message in ValueError is as expected
        self.assertIn("Duplicate node name found: 'C'", str(context.exception))

    def test_no_duplicate_node_names_passes(self):
        # Tree without duplicate node names
        tree_without_duplicates = Tree("((A,B)C,(D,E)F);", format=1)

        # Expect no exception to be raised for a tree without duplicate node names
        try:
            check_and_raise_for_duplicate_nodes(tree_without_duplicates)
        except ValueError:
            self.fail(
                "check_and_raise_for_duplicate_nodes raised ValueError unexpectedly!"
            )


class TestRenameInternalNodesPreOrder(unittest.TestCase):
    def setUp(self) -> None:
        """
        Set up test cases by initializing trees for testing.
        """
        self.base_tree = Tree("(A,(B,C),(D,E));", format=1)
        self.numeric_tree = Tree(
            "(A,(B,3),(D,E));", format=1
        )  # Includes a numeric internal node name
        self.prefixed_tree = Tree(
            "((A:1,B:1)Node2*:1,(D:1,E:1)Node3*:1)Node1*;", format=1
        )  # Already has 'Node' prefix

    def test_rename_internal_nodes_pre_order(self) -> None:
        """
        Test that internal nodes are renamed correctly using pre-order traversal.
        """
        tree_copy = self.base_tree.copy()
        print("Before renaming ")
        renamed_tree = rename_internal_nodes_pre_order(tree_copy)
        print("After renaming (test_rename_internal_nodes_pre_order):")
        print_tree_with_inner_node_names(renamed_tree)

        internal_node_names: List[str] = [
            node.name
            for node in renamed_tree.traverse("preorder")
            if not node.is_leaf()
        ]
        expected_names = ["Node1*", "Node2*", "Node3*"]

        self.assertListEqual(
            internal_node_names,
            expected_names,
            "Internal nodes were not renamed correctly.",
        )

    def test_no_modification_for_existing_node_prefix(self) -> None:
        """
        Ensure trees with existing 'Node' prefix on internal nodes are not modified.
        """
        tree_copy = self.prefixed_tree.copy()
        renamed_tree = rename_internal_nodes_pre_order(tree_copy)
        print(
            "Existing 'Node' prefix tree (test_no_modification_for_existing_node_prefix):"
        )
        print_tree_with_inner_node_names(renamed_tree)

        self.assertEqual(
            str(tree_copy),
            str(renamed_tree),
            "Tree was modified despite having 'Node' prefix.",
        )

    def test_all_internal_nodes_annotated(self) -> None:
        """
        Verify that all internal nodes are considered annotated after renaming.
        """
        tree_copy = self.base_tree.copy()
        renamed_tree = rename_internal_nodes_pre_order(tree_copy)
        print("After renaming (test_all_internal_nodes_annotated):")
        print_tree_with_inner_node_names(renamed_tree)

        self.assertTrue(
            check_all_internal_nodes_annotated(renamed_tree),
            "Not all internal nodes were considered annotated after renaming.",
        )

    def test_not_all_internal_nodes_annotated_before_renaming(self) -> None:
        """
        Confirm that not all internal nodes are considered annotated before the renaming operation.
        """
        print(
            "Before renaming (test_not_all_internal_nodes_annotated_before_renaming):"
        )
        print_tree_with_inner_node_names(self.base_tree)
        self.assertFalse(
            check_all_internal_nodes_annotated(self.base_tree),
            "Incorrectly identified unannotated nodes as annotated before renaming.",
        )


if __name__ == "__main__":
    unittest.main()
