from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


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
