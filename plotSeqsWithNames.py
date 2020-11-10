#@profile
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
                productVector[pos, i] = -1 * (np.multiply(freqPos, logFreqRel, dtype = float))                
    return np.sum(productVector, axis = 1)

#@profile    
def orderPositions(hVector,seqs,haploSize): #order positions
    invH =  np.multiply(-1, hVector, dtype = float)
    ordH = np.argsort(invH)
    # reorder the sequences by entropy
    ordSeqs = {}
    for i in seqs:
        newOne = ''
        for p in range(haploSize):
            newOne = ''.join([newOne, i[ordH[p]]])
        ordSeqs[newOne]=seqs[i]
    return ordSeqs

#@profile
def calcKstep(haploNum,haploSize,seqs): #Calculate hamming distances between sequences
    ordSeqs=seqs.keys()
    colors=[]
    compNum = haploSize
    compList = range(haploNum)
    t = 0
    adjMatrix = np.zeros((haploNum, haploNum))
    kStepList = []
    while compNum > 1:
        t = t + 1
        # Check each query sequence 
        for r1 in range(haploNum-1):
            haplotype1 = ordSeqs[r1]
            for r2 in range(r1+1, haploNum):
                if compList[r1] != compList[r2]: 
                    haplotype2 = ordSeqs[r2]
                    tempDist = 0
                    for a, b in izip(haplotype1, haplotype2):
                        if a != b:
                            tempDist = tempDist + 1
                            if tempDist > t:
                                break
                    if tempDist == t: 
                        adjMatrix[r1][r2] = 1
                        # seqs[ordSeqs[r1]]
                        kStepList.append([seqs[ordSeqs[r1]], seqs[ordSeqs[r2]], t])
                        if seqs[ordSeqs[r1]] not in colors:
                            colors.append(seqs[ordSeqs[r1]])
                        if seqs[ordSeqs[r2]] not in colors:
                            colors.append(seqs[ordSeqs[r2]])
        # Calculate components
        sparseMatrix = csgraph_from_dense(adjMatrix)
        connected = connected_components(sparseMatrix, directed=False, connection='weak', return_labels=True)
        compNum = connected[0]
        compList = connected[1]
        if t/haploSize>.42:
            if sum(compList)==1:
                offender=ordSeqs[np.argmax(compList)]
                splitname=seqs[offender].split("_")
                if splitname[0]=="1":
                    f1counter-=1
                elif splitname[0]=="2":
                    f2counter-=1
                else:
                    both-=1
                del seqs[offender]
                altwrapper(seqs,colors,haploNum-1,haploSize,f1counter,f2counter,both,output,drawmode)
                sys.exit()
            else:
                sys.exit("FATAL ERROR: Your output will contain disconnected components at an unreconcilable distance from one another. Exiting")
    return kStepList,colors,t

#@profile
def kstepToNX(kStepList,colors): #go from a matrix to a weighted networkx graph
    g=nx.Graph()
    numNodes=4
    for seq in colors:
        g.add_node(seq,label=seq)
    for row in kStepList:
        print(row)
        if row[2]>=10:
            g.add_edge(row[0],row[1],len=10,color="red",label=str(row[2]))#,fontsize="24"
        else:
            g.add_edge(row[0],row[1],len=row[2],color="gray30",label=str(row[2]),fontsize="24")
    return g

from itertools import izip
from Bio import SeqIO
import numpy as np
import math
from scipy.sparse.csgraph import connected_components,csgraph_from_dense 
import networkx as nx
import os

seqs={}
with open('test.fas','rU') as f:
    for record in SeqIO.parse(f,"fasta"):
        seqs[record.seq]=record.id
haploNum=4
haploSize=len(seqs.keys()[0])
hVector=calcOrderedFrequencies(haploNum,haploSize,seqs)
ordSeqs=orderPositions(hVector,seqs,haploSize)
kStepList,colors,t=calcKstep(haploNum,haploSize,ordSeqs)
for row in kStepList:
    print(row)
g=kstepToNX(kStepList,colors)
nx.drawing.nx_agraph.write_dot(g,'tmp.dot')
os.system("neato -Tpng tmp.dot -O")
os.system("eog tmp.dot.png")
