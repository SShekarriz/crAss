import sys
import re
from Bio import SeqIO

#python in.fasta out.fasta

# returns seqname, number of Ns, total sequence, percentage of N's in a seq
# only if it finds a good crAss assembly saves it to the second argument
# makes the saved crAss's name shorter

FastaFile = open(sys.argv[1], 'r')
#CRASSpos= open(sys.argv[2], 'w')
for rec in SeqIO.parse(FastaFile, 'fasta'):
	name = rec.id
	seq = rec.seq
	seqLen = len(rec)
	Ns = seq.lower().count('n') 
	Nper= round(Ns / seqLen * 100, 2)
	if seqLen >= 80000 and Nper <= 20:
		rec.id= re.sub('_inter.consensus_threshold.+', '', rec.id)	
		rec.description=''
		CRASSpos= open(sys.argv[2], 'w')
		SeqIO.write(rec, CRASSpos, 'fasta')
	else:
		pass
print(name, Ns, seqLen, Nper)
FastaFile.close()
