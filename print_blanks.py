import os
from Bio import SeqIO

blanks=['_','-','N']
for f in os.listdir(os.getcwd()):
    if f.endswith('fasta'):
        for record in SeqIO.parse(f,'fasta'):
            for xd in blanks:
                if xd in record.seq:
                    print(f)
                    print(record.seq)
