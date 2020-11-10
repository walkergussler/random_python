from Bio import SeqIO
import re
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# files=['AI002', 'AI005', 'AI142', 'AI201', 'AI004', 'AI283', 'AI240', 'AI007', 'AI102', 'AI107', 'AI003', 'AI105', 'AI052', 'AI202', 'AI001']
seqs={}

with open('NH_all.fas') as f:
    for record in SeqIO.parse(f, 'fasta'): # for FASTQ use 'fastq', for fasta 'fasta'
        name=re.findall('([^_]*)_.*',record.id)[0]
        freq=re.findall('_(\d*)$',record.id)[0]
        seqs[name]=[record.seq,freq]
for item in seqs:
    subseqs=[]
    for subitem in seqs[item]:
        [seq,freq]=seqs[item]
        record=SeqRecord(seq,id=item+'_'+freq,description="")
        subseqs.append(record)
    fname=item+'.fas'
    with open(fname,'w') as f:
        SeqIO.write(subseqs,fname,'fasta')

