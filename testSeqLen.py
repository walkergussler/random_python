
from Bio import SeqIO
from Bio.Seq import Seq
seqs=[]
input_handle = open('./indiana/files/IN126_1a.fas') # ('./indiana/files/IN136_1a.fas')
# input_handle2 = open(input2)
for record in SeqIO.parse(input_handle, "fasta"): # for FASTQ use "fastq", for fasta "fasta"
    if len(record.seq) > 0 and len(record.seq) < 50000:
        seqs.append(record.seq)
input_handle.close()

for seq in seqs:
    print(len(seq))