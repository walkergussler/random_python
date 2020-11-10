from Bio import SeqIO
import sys
from collections import defaultdict

seqs=[]
with open(sys.argv[1]) as f:
    for record in SeqIO.parse(f,'fasta'):
        seqs.append(record.seq)
with open(sys.argv[2]) as f:
    for record in SeqIO.parse(f,'fasta'):
        if record.seq in seqs:
            print(record.seq)
