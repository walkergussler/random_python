from Bio import SeqIO
import os,sys

def see(file):
    seqs=[]
    with open(file) as f:
        for record in SeqIO.parse(f,'fasta'):
            return(str(len(record.seq)))

for file in os.listdir(os.getcwd()):
#    if file.startswith('WA'):
    if file.endswith('fasta'):
        print(file+' '+see(file))
