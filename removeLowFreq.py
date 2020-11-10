import sys
from Bio import SeqIO

def read(file,freq):
    seqs=[]
    with open(file) as f:
        for record in SeqIO.parse(f,'fasta'):
            seqFreq=int(record.id.split('_')[-1])
            if seqFreq>=int(freq):
                seqs.append(record)
    return seqs


def main(input,output,freq):
    seqs=read(input,freq)
    print(len(seqs))
    with open(output,'w') as f:
        SeqIO.write(seqs,f,'fasta')

if __name__=='__main__':
    if len(sys.argv)==4:
        main(sys.argv[1],sys.argv[2],sys.argv[3])
    else:
        print('remove sequences with frequency below X')
        print('usage: python removeLowFreq.py <input file> <output file> <trim below this freq>')
