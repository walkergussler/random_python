import numpy as np
import sys, os, re, time, math, glob, optparse, itertools
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from scipy.sparse.csgraph import connected_components, csgraph_from_dense, shortest_path

def parse_input(input): #get sequences from a file
    seqs={}
    with open(input,'r') as f:
        for record in SeqIO.parse(f,'fasta'):
            freq = int(record.id.split('_')[-1])
            if record.seq not in seqs:
                seqs[record.seq]=freq
            else:
                seqs[record.seq]+=freq
    return seqs    

def calc_ordered_frequencies(haploNum,haploSize,seqs,byFreq): 
    freqCount = np.zeros((haploSize, 5))
    productVector = np.zeros((haploSize, 5))
    order={'A':0,'C':1,'G':2,'T':3,'-':4}
    try:
        total_reads=0
        for read in seqs:
            if byFreq:
                freq=seqs[read]
            else:
                freq=1
            total_reads+=freq
            for pos in range(haploSize):
                num=order[read[pos]]
                freqCount[pos, num] = freqCount[pos, num] + freq
        freqRel = np.divide(freqCount, float(total_reads), dtype = float)
    
    except IndexError:
        print("Your files are not aligned and it caused an error! Try again with -a")
    
    for pos in range(haploSize):
        for i in range(5):
            freqPos = freqRel[pos, i]
            if freqPos > 0:
                logFreqRel = math.log(freqPos, 2)
                productVector[pos, i] = -1*(np.multiply(freqPos, logFreqRel, dtype = float))                
    return np.sum(productVector, axis = 1)

def order_positions(hVector,seqs,haploSize): #order positions for faster building of k-step network
    invH =  np.multiply(-1, hVector, dtype = float)
    ordH = np.argsort(invH)
    # reorder the sequences by entropy
    ordSeqs = []

    for i in seqs:
        newOne = ''
        for p in range(haploSize):
            newOne = ''.join([newOne, i[ordH[p]]])
        ordSeqs.append(newOne)
    return ordSeqs

def trimfh(file):
    return os.path.splitext(os.path.basename(file))[0]
    
def main(FILE,DISTANCE_LIMIT,OUTPUT):
    dict=parse_input(FILE)
    seqs=dict.keys()
    counts=dict.values()
    haploNum = len(seqs)
    haploSize=len(seqs[0])
    nReads = sum(counts)
    
    hVector = calc_ordered_frequencies(haploNum,haploSize,dict,False)
    ordSeqs=order_positions(hVector,seqs,haploSize)

    # Calculate distances
    compNum = haploSize
    compList = range(haploNum)
    adjMatrix = np.zeros((haploNum, haploNum))

    # check sequence pair below distLimit
    for r1 in range(haploNum-1):
        haplotype1 = ordSeqs[r1]
        for r2 in range(r1+1, haploNum):
            haplotype2 = ordSeqs[r2]
            tempDist = 0
            for a, b in itertools.izip(haplotype1, haplotype2):
                if a != b:
                    tempDist = tempDist + 1
                    if tempDist > DISTANCE_LIMIT:
                        break
            if tempDist <= DISTANCE_LIMIT: 
                adjMatrix[r1][r2] = 1
                adjMatrix[r2][r1] = 1
    # calculate components
    sparseMatrix = csgraph_from_dense(adjMatrix)
    connected = connected_components(sparseMatrix, directed=False, connection='weak', return_labels=True)
    compNum = connected[0]
    compList = connected[1]
    haploNumList = range(haploNum)
    for i in range(compNum):
        freqSub = []
        indices = []
        idx = 0
        for comp in compList:
            if comp == i:
                freqSub.append(counts[idx])
                indices.append(haploNumList[idx])
            idx = idx+1
        tempFreqCount = np.sum(freqSub)
        tempFreqFraction = np.divide(tempFreqCount, nReads, dtype = float)
        if tempFreqFraction >= .05:
            outputName2 = trimfh(FILE) + '_' + str(i) + '.fas'
            # compSize = len(indices)
            # finalDist = np.zeros((compSize, compSize))
            # adjMatrixComp = np.zeros((compSize, compSize))
            ## calculate distance between every pair of sequences	
            # for s1 in range(compSize-1):
                # haplotype1 = seqs[indices[s1]]
                # for s2 in range(s1+1, compSize):
                    # haplotype2 = seqs[indices[s2]]
                    # tempDist = 0
                    # for a, b in itertools.izip(haplotype1, haplotype2):
                        # if a != b:
                            # tempDist = tempDist + 1
                    # finalDist[s1][s2] = tempDist
                    # finalDist[s2][s1] = tempDist
                    # if tempDist <= DISTANCE_LIMIT:
                        # adjMatrixComp[s1][s2] = 1
                        # adjMatrixComp[s2][s1] = 1
            #save .fas file
            with open(os.path.join(OUTPUT, outputName2), "w") as f:
                n = 0
                for k in indices:
                    haplotype = str(seqs[k])
                    n += 1
                    seq_id = str(n) + '_' + str(counts[k])
                    seq_dna = Seq(haplotype,generic_dna)
                    seqFASTA = SeqRecord(seq_dna, id = seq_id)
                    SeqIO.write(seqFASTA, f, "fasta")

if __name__=="__main__":
    output='1step'
    for file in os.listdir(os.getcwd()):
        if file.endswith('fas'):
            print(file)
            main(file,1,output)