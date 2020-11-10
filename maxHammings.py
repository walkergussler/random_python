import os,re
import numpy as np
from Bio import SeqIO
from ghost.util.distance import hamming
from tempfile import NamedTemporaryFile
import subprocess
import itertools

def prealign(f1,f2):
    # print("prealign")
    uid=re.findall('([^_]*)_.*',os.path.basename(f1))[0]
    with NamedTemporaryFile(delete=False) as tmpcat:
        catname=tmpcat.name
        subprocess.check_call(['cat',f1,f2],stdout=tmpcat)
    with NamedTemporaryFile(delete=False) as aligned:
        alignname=aligned.name
        subprocess.check_call(['mafft', '--quiet', '--auto', '--thread', '20', '--preservecase', catname], stdout=aligned) 
    # os.system("grep -c '>' "+alignname)
    os.unlink(catname)
    seqs=SeqIO.parse(alignname,'fasta')
    current=uid
    fnameIt=0
    seenSecond=False
    seqs1=[]
    while not seenSecond:
        record=next(seqs)
        # print('first')
        thisseq=re.findall('([^_]*)_.*',record.id)[0]
        if thisseq!=uid:
            # print(thisseq,uid)
            seenSecond=True
        else:
            seqs1.append(record)
    with NamedTemporaryFile(delete=False) as align1:
        SeqIO.write(seqs1,align1,'fasta')
        f1name=align1.name
    seqs2=[]
    seqs2.append(record)#write record already accessed and ignored because it doesn't belong in first file
    co2=1
    done=False
    while not done:     
        try:
            record=next(seqs)
            # print('second')
            seqs2.append(record)
        except StopIteration:
            done=True
    with NamedTemporaryFile(delete=False) as align2:
        SeqIO.write(seqs2,align2,'fasta')
        f2name=align2.name
    os.unlink(alignname)
    # print(f1name)
    # print(f2name)
    return(f1name,f2name)


def getseqs(file):
    seqs=[]
    with open(file) as f:
        for record in SeqIO.parse(f, 'fasta'): # for FASTQ use 'fastq', for fasta 'fasta'
            seqs.append(record.seq)
    return seqs

def calcDistanceMatrix(seqs1,seqs2): 
    l1=len(seqs1)
    l2=len(seqs2)
    arr=np.zeros([l1,l2])
    hdist=hamming(seqs1,seqs2,ignore_gaps=False)
    for id in range(len(hdist)):
        item=hdist[id]
        arr[:,id]=item[:,0]
    return arr

def calcDistanceMatrix_slow(seqs1,seqs2):
    l1=len(seqs1)
    l2=len(seqs2)
    arr=np.zeros([l1,l2])
    for id1 in range(len(seqs1)):
        seq1=seqs1[id1]
        for id2 in range(len(seqs2)):
            seq2=seqs2[id2]
            dist = sum(0 if a == b else 1 for a,b in itertools.izip(seq1, seq2))
            arr[id1,id2]=dist
            arr[id2,id1]=dist
    return arr

    
f=open('../close5.csv')
lines=f.readlines()
f.close()
st='_freq5'
for line in lines:
    [p1,p2,hamm,gen1,tmp]=line.split(",")
    gen2=tmp.strip()
    f1=p1+"_"+gen1+st+".fas"
    f2=p2+"_"+gen2+st+".fas"
    dists=[]
    Tdists=[]
    if os.path.isfile(f1) and os.path.isfile(f2):
        # print(f1,f2)
        (f1name,f2name)=prealign(f1,f2)
        seqs1=getseqs(f1name)
        seqs2=getseqs(f2name)
        os.unlink(f1name)
        os.unlink(f2name)
################################################################################################
        # l1=len(seqs1)
        # l2=len(seqs2)
        # mins=[]
        # for seq1 in seqs1:
            # tmplist=[]
            # for seq2 in seqs2:
                # dist = sum(0 if a == b else 1 for a,b in itertools.izip(seq1, seq2))
                # tmplist.append(dist)
            # mins.append(min(tmplist))
        # val=max(mins)/float(len(seqs1[0]))
        # print(p1+','+p2+','+hamm+','+gen1+','+gen2+','+str(val))        
        # for seq1 in seqs2:
            # tmplist=[]
            # for seq2 in seqs1:
                # dist = sum(0 if a == b else 1 for a,b in itertools.izip(seq1, seq2))
                # tmplist.append(dist)
            # mins.append(min(tmplist))
        # val=max(mins)/float(len(seqs1[0]))
        # print(p2+','+p1+','+hamm+','+gen2+','+gen1+','+str(val))
################################################################################################
        mat=calcDistanceMatrix(seqs1,seqs2)
        for row in mat:
            dists.append(min(row))
        val=max(dists)/len(seqs1[0])
        print(p1+','+p2+','+hamm+','+gen1+','+gen2+','+str(val))
        for row in mat.T:
            Tdists.append(min(row))
        val=max(Tdists)/len(seqs1[0])
        print(p2+','+p1+','+hamm+','+gen2+','+gen1+','+str(val))
        