import unittest
from ete3 import Tree
from collections import Counter
# Assuming rescale_branch_lengths is located in a module named mymodule
from satute_trees import (
    rescale_branch_lengths,
    has_duplicate_leaf_sequences,
    collapse_identical_leaf_sequences,
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


if __name__ == "__main__":
    unittest.main()
