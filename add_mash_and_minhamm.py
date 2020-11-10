from __future__ import division
from Bio import SeqIO
from ghost.util.distance import hamming
import subprocess
import numpy as np
import os

def getseqs(file):
    seqs=[]
    with open(file) as f:
        for record in SeqIO.parse(f, 'fasta'): # for FASTQ use 'fastq', for fasta 'fasta'
            seqs.append(record.seq)
    return seqs

def trimfh(handle):
    return os.path.splitext(os.path.basename(handle))[0]

    
def calcDistanceMatrix(seqs1,seqs2): 
    l1=len(seqs1)
    l2=len(seqs2)
    arr=np.zeros([l1,l2])
    hdist=hamming(seqs1,seqs2,ignore_gaps=False)
    for id in range(len(hdist)):
        item=hdist[id]
        arr[:,id]=item[:,0]
    return arr

with open('combinations.csv') as f:
    lines=f.readlines()
with open('mash_hamm.csv','w') as f:
    f.write('f1,f2,mashdist,mindist\n')
    for line in lines:
        [f1,tmp]=line.split(',')
        f2=tmp.strip()
        f1fas='fas/'+f1
        f2fas='fas/'+f2
        f1msh='msh/'+f1+'.msh'
        f2msh='msh/'+f2+'.msh'
        seqs1=getseqs(f1fas)
        seqs2=getseqs(f2fas)
        seqlen=len(seqs1[0])
        mDist=np.amin(calcDistanceMatrix(seqs1,seqs2))/seqlen
        a=subprocess.check_output(['mash','dist',f1msh,f2msh]).split('\t')
        f.write(','.join([trimfh(a[0]),trimfh(a[1]),a[2],str(mDist)])+'\n')
        print(','.join([trimfh(a[0]),trimfh(a[1]),a[2],str(mDist)]))        