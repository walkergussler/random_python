from Bio import SeqIO
import os,sys
from collections import defaultdict

def see(file):
    seqs=defaultdict(int)
    with open(file) as f:
        for record in SeqIO.parse(f,'fasta'):
            sl=len(record.seq)
            freq=int(record.id.split("_")[-1])
            seqs[sl]+=freq
    for item in seqs:
        print(item,seqs[item])

for file in os.listdir(os.getcwd()):
    if file.endswith('fasta') or file.endswith('fas'):
#    if '96' in file:
        print(file)
        see(file)
        
#################TO ALIGN IF THEY ARENT##########################
#        if not see(file):
#            print('mafft --quiet --auto --thread 20 --preservecase '+file+' > aligned/'+file)
#            os.system('mafft --quiet --auto --thread 20 --preservecase '+file+' > aligned/'+file)
#        else:
#            print('cp '+file+' aligned/'+file)
#            os.system('cp '+file+' aligned/'+file)
