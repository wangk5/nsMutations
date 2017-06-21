#input: .maf, .fasta

import sys
import csv

from Bio import SeqIO
from Bio.Seq import MutableSeq
from Bio.Alphabet import IUPAC

input_file = sys.argv[-2]
maf = sys.argv[-1]
fasta_seq = next(SeqIO.parse(input_file, "fasta"))
#mutable_seq = fasta_seq.tomutable()

list_HGVSp = []
list_SWISS = []
list_9mers = []

#generates list of HGVSp and corrseponding SWISSPROT from .maf
maf_data = list(csv.reader(open(maf, 'r'), delimiter='\t'))
numrows = len(maf_data)
for i in range (1,numrows-1):
	if maf_data[i][8] == "Missense_Mutation":
		list_HGVSp.append(maf_data[i][35])
		list_SWISS.append(maf_data[i][67])

ref_seqs = {}
for fasta_seq in SeqIO.parse(input_file, "fasta"):
	id_SWISS = fasta_seq.description[3:9]
	seq_SWISS = fasta_seq.seq
	ref_seqs[id_SWISS] = seq_SWISS

#corresponding letter to three char amino acids
aminos = {"Gly":"G", "Ala":"A", "Leu":"L", "Met":"M", "Phe":"F", 
	      "Trp":"W", "Lys":"K", "Gln":"Q", "Glu":"E", "Ser":"S",
	      "Pro":"P", "Val":"V", "Ile":"I", "Cys":"C", "Tyr":"Y", 
	      "His":"H", "Arg":"R", "Asn":"N", "Asp":"D", "Thr":"T"}

#adds mutation to SWISSPROT sequences
mutated_seqs = []
mut_posit = []
for j in range (0, len(list_HGVSp)):
	if list_SWISS[j] in ref_seqs.keys():
		temp_seq = ref_seqs[list_SWISS[j]]
		mut_temp = temp_seq.tomutable()
		change_id = list_HGVSp[j][-3:]
		change_to = aminos[change_id]
		change_at = int(list_HGVSp[j][5:-3])
		mut_temp[change_at] = change_to #mutating the sequence
		mut_posit.append(change_at)
		mutated_seqs.append(mut_temp)

#generates 9-mers from mutated area
for x in range (0, len(mutated_seqs)):
	if mut_posit[x] < 8:
		mutated_seqs[x] = mutated_seqs[x][0:mut_posit[x]+9]
	if mut_posit[x] > len(mutated_seqs[x]) - 8:
		mutated_seqs[x] = mutated_seqs[x][lmut_posit[x]-8:len(mutated_seqs[x])]
	else:
		mutated_seqs[x] = mutated_seqs[x][mut_posit[x]-8:mut_posit[x]+9]

for y in range (0, len(mutated_seqs)):
	for z in range (0, len(mutated_seqs[y])-8):
		list_9mers.append(mutated_seqs[y][z:z+9])

#output to .txt file
with open("9mers.txt", "w") as txt_file:
	for k in range (0,len(list_9mers)):
		txt_file.write(str(list_9mers[k]) + "\n")
