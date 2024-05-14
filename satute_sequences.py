# -*- coding: utf-8 -*-
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from ete3 import Tree


# Define these sets at a global level to avoid re-creation
DNA_CHARACTERS = set("ACGT")
PROTEIN_CHARACTERS = set("ARNDCQEGHILKMFPSTWYV")


def dict_to_alignment(sequence_dict: dict) -> MultipleSeqAlignment:
    """
    Convert a dictionary of sequences to a MultipleSeqAlignment object.

    Args:
    sequence_dict (dict): A dictionary with sequence identifiers as keys and sequence strings as values.

    Returns:
    MultipleSeqAlignment: The corresponding MultipleSeqAlignment object.
    """
    alignment_list = []
    for id, sequence in sequence_dict.items():
        seq_record = SeqRecord(Seq(sequence), id=id)
        alignment_list.append(seq_record)
    return MultipleSeqAlignment(alignment_list)


class TaxaMismatchError(Exception):
    """Custom exception to be raised when the taxa sets do not match."""
    pass

def check_if_alignment_has_same_taxa_as_msa(sequence_alignment: MultipleSeqAlignment, tree: Tree) -> None:
    """
    Check if the taxa set in the multiple sequence alignment is the same as the taxa set in the tree.

    Args:
        sequence_alignment (MultipleSeqAlignment): The multiple sequence alignment object.
        tree (Tree): The phylogenetic tree object.

    Raises:
        TaxaMismatchError: If the taxa sets in the alignment and the tree do not match, 
                           explaining which taxa are missing or extra.
    """
    # Extract taxa from the multiple sequence alignment
    alignment_taxa = set(record.id for record in sequence_alignment)
    
    # Extract taxa from the tree using ete3
    tree_taxa = set(leaf.name for leaf in tree.get_leaves())
    
    # Check if the sets of taxa are the same
    if alignment_taxa != tree_taxa:
        missing_in_tree = alignment_taxa - tree_taxa
        extra_in_tree = tree_taxa - alignment_taxa
        error_message = "The taxa sets in the alignment and the tree do not match. "
        
        if missing_in_tree:
            error_message += f"Taxa missing in tree: {', '.join(missing_in_tree)}. "
        
        if extra_in_tree:
            error_message += f"Extra taxa in tree: {', '.join(extra_in_tree)}."
        
        raise TaxaMismatchError(error_message)
