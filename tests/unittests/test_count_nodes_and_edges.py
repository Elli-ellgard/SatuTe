import sys
sys.path.append("./../../")

import unittest
from ete3 import Tree
from typing import List
from satute.graph import calculate_subtree_edge_metrics


class TestCountNodeAndEdges(unittest.TestCase):
    def setUp(self):
        # Set up a simple tree and a more complex tree
        self.five_taxa_tree = Tree(
            "((A:1,B:1)N1,(C:1,(D:1,E:1)N4:1)N3:1)Root;", format=1
        )
        self.three_taxa_tree = Tree("((A:1,B:1)N1:1,C:1)Root;", format=1)
        self.two_taxa_tree = Tree("(A:1,B:1)R:1;", format=1)

    def test_three_taxa_tree_2(self):
        edge_metrics = calculate_subtree_edge_metrics(self.five_taxa_tree, None)

        expected_edge_metrics = {
            "(N1, Root)": {
                "left": {
                    "leave_count": 2,
                    "branch_count": 2,
                },
                "right": {
                    "leave_count": 3,
                    "branch_count": 5,
                },
                "length": 1,  # Assume branch length is 1 for simplicity
                "type": "internal",
            },
            "(A, N1)": {
                "left": {
                    "leave_count": 1,
                    "branch_count": 0,
                },
                "right": {
                    "leave_count": 4,
                    "branch_count": 7,
                },
                "length": 1,
                "type": "external",
            },
            "(B, N1)": {
                "left": {
                    "leave_count": 1,
                    "branch_count": 0,
                },
                "right": {
                    "leave_count": 4,
                    "branch_count": 7,
                },
                "length": 1,
                "type": "external",
            },
            "(N3, Root)": {
                "left": {
                    "leave_count": 3,
                    "branch_count": 4,
                },
                "right": {
                    "leave_count": 2,
                    "branch_count": 3,
                },
                "length": 1,
                "type": "internal",
            },
            "(C, N3)": {
                "left": {
                    "leave_count": 1,
                    "branch_count": 0,
                },
                "right": {
                    "leave_count": 4,
                    "branch_count": 7,
                },
                "length": 1,
                "type": "external",
            },
            "(D, N4)": {
                "left": {
                    "leave_count": 1,
                    "branch_count": 0,
                },
                "right": {
                    "leave_count": 4,
                    "branch_count": 7,
                },
                "length": 1,
                "type": "external",
            },
        }

        # Convert keys to a tuple format that matches the output from the function
        self.assertEdgeMetrics(edge_metrics, expected_edge_metrics)

    def test_three_taxa_tree(self):
        edge_metrics = calculate_subtree_edge_metrics(self.three_taxa_tree, None)

        expected_edge_metrics = {
            "(C, Root)": {
                "left": {
                    "leave_count": 1,
                    "branch_count": 0,
                },
                "right": {
                    "leave_count": 2,
                    "branch_count": 3,
                },
                "length": 1,  # Assume branch length is 1 for simplicity
                "type": "external",
            },
            "(A, N1)": {
                "left": {
                    "leave_count": 1,
                    "branch_count": 0,
                },
                "right": {
                    "leave_count": 2,
                    "branch_count": 3,
                },
                "length": 1,
                "type": "external",
            },
            "(B, N1)": {
                "left": {
                    "leave_count": 1,
                    "branch_count": 0,
                },
                "right": {
                    "leave_count": 2,
                    "branch_count": 3,
                },
                "length": 1,
                "type": "external",
            },
        }

        # Convert keys to a tuple format that matches the output from the function
        self.assertEdgeMetrics(edge_metrics, expected_edge_metrics)

    def assertEdgeMetrics(self, actual, expected):
        for key, metrics in expected.items():
            with self.subTest(edge=key):
                # Ensure the edge exists
                self.assertIn(key, actual, f"Missing edge: {key}")

                # Check metrics
                for metric_key in ["left", "right"]:
                    self.assertEqual(
                        actual[key][metric_key]["leave_count"],
                        metrics[metric_key]["leave_count"],
                        f"Incorrect leave_count for {key} on {metric_key}",
                    )
                    self.assertEqual(
                        actual[key][metric_key]["branch_count"],
                        metrics[metric_key]["branch_count"],
                        f"Incorrect branch_count for {key} on {metric_key}",
                    )

                # Check length and type
                self.assertEqual(
                    actual[key]["length"],
                    metrics["length"],
                    f"Incorrect length for {key}",
                )
                self.assertEqual(
                    actual[key]["type"], metrics["type"], f"Incorrect type for {key}"
                )


if __name__ == "__main__":
    unittest.main()
