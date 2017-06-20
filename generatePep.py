#input: .maf, .fasta

import sys
import csv

from Bio import SeqIO
from Bio.Seq import MutableSeq
from Bio.Alphabet import IUPAC

input_file = sys.argv[1]
maf = sys.argv[2]
fasta_seq = next(SeqIO.parse(input_file, "fasta"))
#mutable_seq = fasta_seq.tomutable()

list_HGVSp = []
list_SWISS = []
list_9mers = []

#generates list of HGVSp and corrseponding SWISSPROT from .maf
maf_data = list(csv.reader(open('maf', 'r'), delimiter='\t'))
numrows = len(maf_data)
for i in numrows:
	if maf_data[8][i] == "Missense_Mutation":
		list_HGVSp.append(maf_data[35][i])
		list_SWISS.append(maf_data[67][i])

ref_seqs = []
for fasta_seq in SeqIO.parse(input_file, "fasta"):
	id_SWISS = fasta_seq.description[3:8]
	seq_SWISS = fasta_seq.seq
	ref_seqs.append(id_SWISS:seq_SWISS)

#adds mutation to SWISSPROT sequence
aminos = ["Gly":"G", "Ala":"A", "Leu":"L", "Met":"M", "Phe":"F", 
	  "Trp":"W", "Lys":"K", "Gln":"Q", "Glu":"E", "Ser":"S",
	  "Pro":"P", "Val":"V", "Ile":"I", "Cys":"C", "Tyr":"Y", 
	  "His":"H", "Arg":"R", "Asn":"N", "Asp":"D", "Thr":"T"]

mutated_seqs = []
mut_posit = []
for j in len(list_HGVSp):
	temp_seq = ref_seqs[list_SWISS[j]]
	mut_temp = temp_seq.tomutable()
	change_id = list_HGVSp[j][6:8]
	change_to = aminos[change_id]
	change_at = list_HGVSp[j][5:-4]
	mut_temp[change_at] = change_to #mutating the sequence

	mut_posit.append(change_at)
	mutated_seqs.append(mut_temp)

#generates 9-mers from mutated area
for x in len(mutated_seqs):
	if mut_posit[x] < 8:
		mutated_seqs[x] = mutated_seqs[x][0:mut_posit[x+8]]
	else:
		mutated_seqs[x] = mutated_seqs[x][mut_posit[x-8]:mut_posit[x+8]]

for y in len(mutated_seqs):
	for z in len(mutated_seqs[y])
		list_9mers.append(sub_seq[z:z+8])
