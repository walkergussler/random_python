import numpy as np
import sys
from ghost.util.distance import hamming
from Bio import SeqIO

def doHdistReturnProp(seqs1,seqs2): #calculate hamming proportions between two sets of sequences, return matrix
    keylen=len(seqs1[0])
    l1=len(seqs1)
    l2=len(seqs2)
    hdist=hamming(seqs1,seqs2,ignore_gaps=False)
    arr=np.zeros([l1,l2])
    for id in range(len(hdist)):
        item=hdist[id]
        arr[:,id]=item[:,0]
    return np.divide(arr,keylen,dtype=float)

def getseqs(input1):
    seqs=[]
    input_handle = open(input1) 
    for record in SeqIO.parse(input_handle, "fasta"): # for FASTQ use "fastq", for fasta "fasta"
        if len(record.seq) > 0 and len(record.seq) < 50000:
            seqs.append(record.seq)
    input_handle.close()
    return seqs

def printFiles(inmat,fileList): #takes in square matrices only
    # print(fileList)
    if len(fileList)==1:
        print(fileList[0])
        return 0
    maxes=np.zeros(len(fileList))
    for id in range(len(fileList)):
        # print("id=%i"%id)
        newarr=np.delete(np.delete(inmat,id,axis=1),id,axis=0)
        # print(newarr)
        (values,vectors)=np.linalg.eig(newarr)
        # print(values)
        maxes[id]=values.real.max()
        # print("-")
        # print(max(values),fileList[id])
        # print("-")
    # print(maxes)
    smallestSpec=maxes.real.argmax()
    newfiles=list(fileList)
    # print("its here")
    print(fileList[smallestSpec])
    # print("now its not")
    newfiles.remove(fileList[smallestSpec])
    newarr=np.delete(np.delete(inmat,smallestSpec,axis=1),smallestSpec,axis=0)
    printFiles(newarr,newfiles)
    
del sys.argv[0]
inputs=sys.argv
numFiles=len(inputs)
# print("---------running following inputs---------")
dsamp=np.zeros([numFiles,numFiles])
for i1 in range(numFiles):
    f1=inputs[i1]
    # print(f1)
    seqs1=getseqs(f1)
    l1=len(seqs1)
    for i2 in range(i1,numFiles):
        f2=inputs[i2]
        if i1!=i2:
            seqs2=getseqs(f2)
            l2=len(seqs2)
            tmp=doHdistReturnProp(seqs1,seqs2)
            asd=np.amin(tmp)
            if asd==0:
                dsamp[i1,i2]=.00001
                dsamp[i2,i1]=.00001
            else:
                dsamp[i1,i2]=asd
                dsamp[i2,i1]=asd
# print("------------------------------------------")
printFiles(dsamp,inputs)