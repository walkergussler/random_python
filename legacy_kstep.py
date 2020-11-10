from __future__ import absolute_import, division


import numpy as np
import math,sys,os,re
from scipy.sparse.csgraph import connected_components,csgraph_from_dense 
from itertools import izip
from Bio import SeqIO
from ghost.util.distance import hamming
import networkx as nx
import collections

__author__ = "CDC/OID/NCHHSTP/DVH bioinformatics team"

#REQUIREMENTS FROM INPUT
    #needs each read id to end with "_XXXX" where XXXX is the frequency of that read
    #needs file name to end with "XX.fas" where XX is a valid HCV genotype (1b,5a,etc.)

#better legend
    #number of dots for each color
    #genotype
    #min hamming
    
def parseInput(input1,input2): #get sequences from 2 files
    ali=True
    seqs={}
    input_handle = open(input1) 
    input_handle2 = open(input2)
    f1counter=0
    f2counter=0
    fbothcounter=0
    for record in SeqIO.parse(input_handle, "fasta"): # for FASTQ use "fastq", for fasta "fasta"
        if len(record.seq) > 0 and len(record.seq) < 50000:
            f1counter+=1
            freq=re.findall('_(\d*)$',record.id)
            seqs[record.seq]="1_"+str(f1counter)+"_"+freq[0]
    input_handle.close()
    for record in SeqIO.parse(input_handle2, "fasta"): # for FASTQ use "fastq", for fasta "fasta"
        if len(record.seq) > 0 and len(record.seq) < 50000:
            if record.seq not in seqs:
                f2counter+=1
                freq=re.findall('_(\d*)$',record.id)
                seqs[record.seq]="2_"+str(f2counter)+"_"+freq[0]
            else:
                fbothcounter+=1
                seqs[record.seq]="x_"+str(fbothcounter)       
    input_handle2.close()
    if ali:
        newseqs=align(seqs)
        try:
            haploSize = len(newseqs.keys()[0])
        except IndexError:
            sys.exit("Either mafft isn't loaded or both of your fasta files are empty. Exiting")
        haploNum = len(newseqs)
        return(newseqs,haploSize,haploNum,f1counter,f2counter,fbothcounter)
    else:
        haploSize = len(seqs.keys()[0])
        haploNum = len(seqs)
        return(seqs,haploSize,haploNum,f1counter,f2counter,fbothcounter)
    
def align(seqs):
    f=open("tmp.fas","w")
    for seq in seqs:
        f.write(">"+seqs[seq]+"\n")
        f.write(str(seq)+"\n")
    f.close()
    os.system('mafft --quiet --auto --thread 20 --preservecase tmp.fas > align.fas')
    outs={}
    input_handle = open('align.fas')
    for record in SeqIO.parse(input_handle, "fasta"): # for FASTQ use "fastq", for fasta "fasta"
        if len(record.seq) > 0 and len(record.seq) < 50000:
            if collections.Counter(record.seq).most_common(1)[0][0]!="-":
                outs[record.seq]=record.id
    input_handle.close()
    # rmstr="rm align.fas tmp.fas"
    # os.system(rmstr)
    # for item in outs:
        # print(item)
        # print(outs[item])
    return outs
    
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
                productVector[pos, i] = -1*(np.multiply(freqPos, logFreqRel, dtype = float))                
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
def calcDistances(haploNum,haploSize,asdf): #Calculate hamming distances between sequences
    ordSeqs=asdf.keys()
    # for seq in ordSeqs:
        # print(seq,asdf[seq])
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
                        asdf[ordSeqs[r1]]
                        kStepList.append([asdf[ordSeqs[r1]], asdf[ordSeqs[r2]], t])
                        if asdf[ordSeqs[r1]] not in colors:
                            colors.append(asdf[ordSeqs[r1]])
                        if asdf[ordSeqs[r2]] not in colors:
                            colors.append(asdf[ordSeqs[r2]])
        # Calculate components
        sparseMatrix = csgraph_from_dense(adjMatrix)
        connected = connected_components(sparseMatrix, directed=False, connection='weak', return_labels=True)
        compNum = connected[0]
        compList = connected[1]
        # print(compNum)
        # print(compList)
        # print(t)
        # print("=")
        if t/haploSize>.42:
            if sum(compList)==1:
                offender=ordSeqs[np.argmax(compList)]
                del asdf[offender]
                altwrapper(asdf,colors,haploNum-1,haploSize)
                sys.exit()
            else:
                sys.exit("Your output is going to look really messed up. Exiting")
    return kStepList,colors
  
def kstepToNX(kStepList,colors): #go from a matrix to a weighted networkx graph
    labels=False
    g=nx.Graph()
    for node in colors:
        splitnode=node.split("_")
        if splitnode[0]=="2":
            g.add_node(node,color='blue',shape='point',width=.1)
        if splitnode[0]=="1":
            g.add_node(node,color='red',shape='point',width=.1)
        if splitnode[0]=="x":
            g.add_node(node,color='green',shape='point',width=.1)
    for row in kStepList:
        if labels:
            g.add_edge(row[0],row[1],len=row[2],color='gray',label=str(row[2]))
        else:
            g.add_edge(row[0],row[1],len=row[2],color='gray')
                    # g.add_edge(row,col,len=matrix[row,col]*numRows,arrowsize=".01",color='gray',alpha='.5')
    return g

def splitseqs(seqs): #get sequences from 2 files
    firstFile=[]
    secondFile=[]
    # for key in seqs:
        # print(seqs[key])
    for key in seqs:
        if seqs[key].startswith("1"):
            firstFile.append(key)
        elif seqs[key].startswith("2"):
            secondFile.append(key)
        elif seqs[key].startswith("x"):
            firstFile.append(key)
            secondFile.append(key)
    return (firstFile,secondFile)

def calcDistanceMatrix(seqs1,seqs2): #calculate distance matrix from the 1-step list
    hdist=hamming(seqs1,seqs2,ignore_gaps=False)
    l=len(seqs1)
    w=len(seqs2)
    arr=np.zeros([l,w])
    for id in range(len(hdist)):
        item=hdist[id]
        arr[:,id]=item[:,0]
    return arr

def altwrapper(seqs,colors,haploNum,haploSize):
    # print("-------------altwrapper-------------")
    # print("-------------altwrapper-------------")
    # print("-------------altwrapper-------------")
    aligned=align(seqs)
    hVector=calcOrderedFrequencies(haploNum,haploSize,seqs)
    ordSeqs=orderPositions(hVector,seqs,haploSize)
    kStepList,colors=calcDistances(haploNum,haploSize,ordSeqs)
    g=kstepToNX(kStepList,colors)
    nx.drawing.nx_agraph.write_dot(g,'tmp.dot')
    f=open("tmp.dot","r")
    try:
        g=open(sys.argv[3],"w")
    except IndexError:
        g=open("tmper.dot","w")
    lines=f.readlines()
    f.close()
    lineco=0
    for line in lines:
        lineco+=1
        if lineco==2:
            g.write("\tgraph [outputorder=edgesfirst];\n")
        g.write(line)
    g.close()
    try:
        grapit="neato -Tpng "+sys.argv[3]+" -O"
    except IndexError:
        grapit="neato -Tpng tmper.dot -o tmp.png"
    os.system(grapit)
    # os.system("eog tmp.png")
    
if __name__=="__main__":
    checkMinHamm=False
    gen1=re.findall('_(\d[a-z]).fas',sys.argv[1])[0]
    gen2=re.findall('_(\d[a-z]).fas',sys.argv[2])[0]
    if gen1!=gen2:
        print("WARNING: your two files are not the same genotype! Your output will likely be bad")
    (seqs,haploSize,haploNum,f1counter,f2counter,both)=parseInput(sys.argv[1],sys.argv[2])
    (seqs1,seqs2)=splitseqs(seqs)
    array=calcDistanceMatrix(seqs1,seqs2)
    minHamm=np.amin(array)/len(seqs1[0])
    if checkMinHamm==True:
        if minHamm>.42:
            sys.exit("Your files do not appear to be the same genotype - manual curation likely required. Exiting")
    hVector=calcOrderedFrequencies(haploNum,haploSize,seqs)
    ordSeqs=orderPositions(hVector,seqs,haploSize)
    kStepList,colors=calcDistances(haploNum,haploSize,ordSeqs)
    g=kstepToNX(kStepList,colors)
    nx.drawing.nx_agraph.write_dot(g,'tmp.dot')
    f=open("tmp.dot","r")
    try:
        g=open(sys.argv[3],"w")
    except IndexError:
        g=open("tmper.dot","w")
    lines=f.readlines()
    f.close()
    lineco=0
    l=len(lines)
    for line in lines:
        lineco+=1
        if lineco==2:
            g.write("\tgraph [outputorder=edgesfirst];\n")
        if lineco!=l:
            g.write(line)
        else:
            g.write("{ rank = sink;\n")
            g.write("    Legend [shape=none, margin=0, label=<\n")
            g.write("""    <TABLE BORDER="0" CELLBORDER="1" CELLSPACING="0" CELLPADDING="4">\n""")
            g.write("     <TR>\n")
            g.write("""      <TD COLSPAN="2"><B>Genotype = %s, Min H-Dist = %f</B></TD>\n""" %(gen1,minHamm))
            g.write("     </TR>\n     <TR>\n")
            g.write("      <TD>%s</TD>\n" %sys.argv[1])
            g.write("""      <TD><FONT COLOR="red">%i</FONT></TD>\n"""%f1counter)
            g.write("     </TR>\n     <TR>\n")
            g.write("      <TD>%s</TD>\n" %sys.argv[2])
            g.write("""      <TD><FONT COLOR="blue">%i</FONT></TD>\n"""%f2counter)
            g.write("     </TR>\n     <TR>\n")
            g.write("      <TD>Both Files</TD> ")
            g.write("""      <TD><FONT COLOR="green">%i</FONT></TD>\n"""%both)
            g.write("     </TR>\n    </TABLE>\n   >];\n  }\n} ")
    g.close()
    try:
        grapit="neato -Tpng "+sys.argv[3]+" -O"
    except IndexError:
        grapit="neato -Tpng tmper.dot -o tmp.png"
    os.system(grapit)
    # os.system("eog tmp.png")