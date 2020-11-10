from Bio import SeqIO
import sys
from collections import defaultdict

seqs=defaultdict(int)
with open(sys.argv[1]) as f:
    for record in SeqIO.parse(f,'fasta'):
        seqs[record.seq]+=1
a=seqs.values()
failco=0
for item in a:
    if item!=1:
        failco+=1
print(str(failco)+' sequences appear more than once')