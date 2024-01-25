from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
sys.path.append('..')
from satute_sequences import classify_alignment_type, SequenceType

# Helper function to create a SeqRecord
def create_seq_record(sequence):
    return SeqRecord(Seq(sequence), id="")



# Test case 1: All DNA sequences
alignment_dna = MultipleSeqAlignment(
    [
        create_seq_record("AGTACGT"),
        create_seq_record("AGTACGT"),
        create_seq_record("AGTACGT"),
    ]
)
result_dna = classify_alignment_type(alignment_dna)
assert result_dna == SequenceType.DNA
print("Test case 1 passed!", result_dna)

# Test case 2: All protein sequences
alignment_protein = MultipleSeqAlignment(
    [
        create_seq_record("FPSTWYV"),
        create_seq_record("FPSTWYV"),
        create_seq_record("FPSTWYV"),
    ]
)
result_protein = classify_alignment_type(alignment_protein)
assert result_protein == SequenceType.PROTEIN
print("Test case 2 passed!", result_protein)

# Test case 3: Mixed DNA and protein sequences
alignment_mixed = MultipleSeqAlignment(
    [
        create_seq_record("ACGTTT"),
        create_seq_record("ARNDAA"),
        create_seq_record("AGTACG"),
    ]
)
result_mixed = classify_alignment_type(alignment_mixed)
assert result_mixed == SequenceType.PROTEIN


# Test case 4: Empty alignment
# alignment_empty = MultipleSeqAlignment([])
# result_empty = classify_alignment_type(alignment_empty)
# assert result_empty == SequenceType.UNKNOWN
