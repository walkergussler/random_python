from Bio import SeqIO
from sys import argv
from itertools import izip
seqs=[]
with open(argv[1]) as f:
    for record in SeqIO.parse(f,'fasta'):
        seqs.append(record.seq)
        if len(seqs)==2:
            break
seq1=seqs[0]
seq2=seqs[1]
dist = sum(0 if a == b else 1 for a,b in izip(seq1, seq2))
print(dist)
