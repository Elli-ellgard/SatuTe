import sys

sys.path.append("./../..")

import unittest
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from ete3 import Tree
from satute_sequences import check_if_alignment_has_same_taxa_as_msa, TaxaMismatchError


# Import the function and custom exception from the module where they are defined
# from your_module import check_if_alignment_has_same_taxa_as_msa, TaxaMismatchError


class TestCheckIfAlignmentHasSameTaxaAsMSA(unittest.TestCase):

    def setUp(self):
        # Create a MultipleSeqAlignment object with some sequences
        self.alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("ATGC"), id="taxon1"),
                SeqRecord(Seq("ATGC"), id="taxon2"),
                SeqRecord(Seq("ATGC"), id="taxon3"),
            ]
        )

        # Create a Tree object with matching taxa
        self.tree_matching = Tree(name="root")
        self.tree_matching.add_child(name="taxon1")
        self.tree_matching.add_child(name="taxon2")
        self.tree_matching.add_child(name="taxon3")

        # Create a Tree object with non-matching taxa
        self.tree_non_matching = Tree(name="root")
        self.tree_non_matching.add_child(name="taxon1")
        self.tree_non_matching.add_child(name="taxon2")
        self.tree_non_matching.add_child(name="taxon4")  # Different taxon
        
        # Create a Tree object with partially matching taxa
        self.tree_partial_matching = Tree(name="root")
        self.tree_partial_matching.add_child(name="taxon1")
        self.tree_partial_matching.add_child(name="taxon2")        

    def test_matching_taxa(self):
        try:
            check_if_alignment_has_same_taxa_as_msa(self.alignment, self.tree_matching)
        except TaxaMismatchError:
            self.fail("TaxaMismatchError raised unexpectedly with matching taxa.")

    def test_non_matching_taxa(self):
        with self.assertRaises(TaxaMismatchError) as cm:
            check_if_alignment_has_same_taxa_as_msa(
                self.alignment, self.tree_non_matching
            )

        # Check specific parts of the error message
        exception_message = str(cm.exception)
        self.assertIn(
            "The taxa sets in the alignment and the tree do not match.",
            exception_message,
        )

        self.assertIn(
            "Taxa missing in tree: taxon3.", exception_message
        )  # Missing in tree

        self.assertIn("Extra taxa in tree: taxon4.", exception_message)  # Extra in tree

    def test_partial_matching_taxa(self):
        with self.assertRaises(TaxaMismatchError) as cm:
            check_if_alignment_has_same_taxa_as_msa(self.alignment, self.tree_partial_matching)
                    
        # Check specific parts of the error message
        exception_message = str(cm.exception)
        self.assertIn("The taxa sets in the alignment and the tree do not match.", exception_message)
        self.assertIn("Taxa missing in tree: taxon3.", exception_message)  # Missing in tree
        self.assertNotIn("Extra taxa in tree:", exception_message)  # No extra taxa

if __name__ == "__main__":
    unittest.main()
