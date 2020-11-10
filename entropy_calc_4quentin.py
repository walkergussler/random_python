def parseInput(input1,input2): #get sequences from 2 files
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
    else:
        print("%s is the source"%input2)
    
    

def calcOrderedFrequencies(haploNum,haploSize,seqs): #calculate the ordered frequencies    
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

# def decideWhomsGiver(input1,input2):
    # [haploNum,haploSize,seqs,f1HaploNum]=parseInput(input1,input2)
    # f2HaploNum=haploNum-f1HaploNum
    # hVector=calcOrderedFrequencies(haploNum,haploSize,seqs)
    # print(hVector)
    # print(np.shape(hVector))

import sys
from Bio import SeqIO
import numpy as np
import math

parseInput(sys.argv[1],sys.argv[2])