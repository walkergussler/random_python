from Bio import SeqIO
import os,sys

def see(file):
    seqs=[]
    with open(file) as f:
        for record in SeqIO.parse(f,'fasta'):
            seqs.append(record.seq)
    seqlen=len(seqs[0])
    for seq in seqs:
        if len(seq)!=seqlen:
            return(False)
    return(True)

for file in os.listdir(os.getcwd()):
    if file.endswith('fas'):
        if not see(file):
            print('mafft --quiet --auto --thread 20 --preservecase '+file+' > aligned/'+file)
            os.system('mafft --quiet --auto --thread 20 --preservecase '+file+' > aligned/'+file)
        else:
            print('cp '+file+' aligned/'+file)
            os.system('cp '+file+' aligned/'+file)
