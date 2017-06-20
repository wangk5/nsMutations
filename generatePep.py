import sys

from Bio import SeqIO
from Bio.Seq import MutableSeq
from Bio.Alphabet import IUPAC

input_file = sys.argv[-1]
fasta_seq = next(SeqIO.parse(input_file, "fasta"))
list_9mers = []
mutable_seq = fasta_seq.tomutable()

# add mutations here

for x in len(mutable_seq):
	base = mutable_seq[x]
	sub_seq = ""
	if # mutation found:
		sub_seq = mutable_seq[x-8:x+8]
		for y in len(sub_seq):
			list_9mers.append(sub_seq[y:y+8])
