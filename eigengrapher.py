#proposition on wiki page: potentially massive speed bonus
def whichIsSourceByEntropy(input1,input2): #get sequences from 2 files
    seqs1=[]
    seqs2=[]
    input_handle = open(input1) 
    input_handle2 = open(input2)
    for record in SeqIO.parse(input_handle, "fasta"): # for FASTQ use "fastq", for fasta "fasta"
        if len(record.seq) > 0 and len(record.seq) < 50000:
            seqs1.append(record.seq)
    input_handle.close()
    f1HaploNum=len(seqs1)
    for record in SeqIO.parse(input_handle2, "fasta"): # for FASTQ use "fastq", for fasta "fasta"
        if len(record.seq) > 0 and len(record.seq) < 50000:
            if record.seq not in seqs1:
                seqs2.append(record.seq)
            else:
                print("At least one of the sequences in f1 is identical to a sequence in f2; skipping")
    input_handle.close()
    input_handle2.close()
    print("-")
    haploSize1 = len(seqs1[0])
    print("-")
    haploSize2 = len(seqs2[0])
    haploNum1 = len(seqs1)
    haploNum2 = len(seqs2)
    entropy1=calcOrderedFrequencies(haploNum1,haploSize1,seqs1)
    entropy2=calcOrderedFrequencies(haploNum2,haploSize2,seqs2)
    e1=sum(entropy1)/len(entropy1)
    e2=sum(entropy2)/len(entropy2)
    if e1>e2:
        print("%s is the source"%input1)
        return 1
    else:
        print("%s is the source"%input2)
        return 2
    
def calcOrderedFrequencies(haploNum,haploSize,seqs): #change to order by eigenvalues?
    freqCount = np.zeros((haploSize, 5))
    productVector = np.zeros((haploSize, 5))
    for read in seqs:
        for pos in range(haploSize):
            if read[pos] == 'A':
                freqCount[pos, 0] = freqCount[pos, 0] + 1
            elif read[pos] == 'C':
                freqCount[pos, 1] = freqCount[pos, 1] + 1
            elif read[pos] == 'G':
                freqCount[pos, 2] = freqCount[pos, 2] + 1
            elif read[pos] == 'T':
                freqCount[pos, 3] = freqCount[pos, 3] + 1
            elif read[pos] == '-':
                freqCount[pos, 4] = freqCount[pos, 4] + 1
    freqRel = np.divide(freqCount, haploNum, dtype = float)
    for pos in range(haploSize):
        for i in range(5):
            freqPos = freqRel[pos, i]
            if freqPos > 0:
                logFreqRel = math.log(freqPos, 2)
                productVector[pos, i] = -1*(np.multiply(freqPos, logFreqRel, dtype = float))                
    return np.sum(productVector, axis = 1)

def makeEntropyGraph(inputs): #get sequences from 2 files
    l=len(inputs)
    pos=np.arange(l)
    entropies={}
    for i in range(l):
        seqs=[]
        file=inputs[i]
        input_handle = open(file) 
        for record in SeqIO.parse(input_handle, "fasta"): # for FASTQ use "fastq", for fasta "fasta"
            if len(record.seq) > 0 and len(record.seq) < 50000:
                seqs.append(record.seq)
        input_handle.close()
        haploNum=len(seqs)
        haploSize=len(seqs[0])
        entropy=calcOrderedFrequencies(haploNum,haploSize,seqs)
        ent=sum(entropy)/len(entropy)
        entropies[file]=ent
    return entropies,pos

def second_largest(numbers):
    count = 0
    m1 = m2 = float('-inf')
    for x in numbers:
        count += 1
        if x > m2:
            if x >= m1:
                m1, m2 = x, m1            
            else:
                m2 = x
    return m2 if count >= 2 else None
    
import sys
from Bio import SeqIO
import numpy as np
import math
import matplotlib.pyplot as plt
import re
from scipy import stats
del sys.argv[0]
a,pos=(makeEntropyGraph(sys.argv))
asdf=a.values()
second=float(second_largest(asdf))
m=float(max(asdf))
W=float(m/second)
print(stats.zscore(asdf))
print(max(stats.zscore(asdf)))
stor=max(asdf)
asdf.remove(max(asdf))
print(stats.zscore(asdf))
print(max(stats.zscore(asdf)))
print(stor/max(asdf))
objects=[]
s=a.values()
for v in a:
    name=re.findall('NH(.*)_unique.*',v)
    objects.append(name[0])
    print(name,a[v])
performance=a.values()
plt.bar(pos,performance,align='center',alpha=.5)
plt.xticks(pos,objects)
plt.ylabel('Average Entropy')
plt.xlabel('File Name')
plt.show()


# python eigengrapher.py ../Related_clipped/BB*
# python eigengrapher.py ../Related_clipped/BC*
# python eigengrapher.py ../Related_clipped/BJ*


def main(infils):
    l=len(inputs)
    