from Bio import SeqIO
import os, sys

def read_and_count(file):
    c=0
    with open(file) as f:
        for record in SeqIO.parse(f,'fasta'):
            freq=int(record.id.split('_')[-1])
            c+=freq
    return c

print(read_and_count(sys.argv[1]))
