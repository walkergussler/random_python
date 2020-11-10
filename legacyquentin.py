from __future__ import absolute_import
import numpy as np
import networkx as nx
import math
import sys
import random
import fastcluster
import os
import re
import time
import six
import subprocess
import itertools
import lsqlin
from scipy.sparse.csgraph import connected_components,csgraph_from_dense,dijkstra 
from scipy.spatial.distance import pdist
from scipy.stats import zscore
from Bio import SeqIO
from ghost.util.distance import hamming
from shutil import copyfileobj
from tempfile import NamedTemporaryFile
from cvxopt.solvers import options

__author__ = "CDC/OID/NCHHSTP/DVH bioinformatics team"

"""QUENTIN, by J Walker Gussler. Translated from matlab version by Pavel Skums. (https://github.com/skumsp/quentin)
Given a number of .fas clinical HCV samples of HVR1, this program will calculate clusters (does an ok job)
and then it will calculate the directionality of infection (this is what the program is really for)
 
1. It requires as input the location of the aligned fasta files. All sequences are assumed to be different.
2. It calculates the k-step network using an iterative multi-index hashing approach.
3. It calculates all intermediate sequences so that we have a 1-step network (with some fabricated sequences) instead of a k-step network
4. From this 1-step network with intermediate sequences, it computes an all-against-all hamming distance matrix.
5. From this hamming distance matrix, it computes expected evolution times for file A to evolve into file B and vice versa 
6. Within an input cluster of n files, a nxn matrix of evolution simulation times is constructed
7. From these simulated evolution times, it iteratively constructs a phylogenetic tree and infers linkage data from the tree

Should be inside ghost virtual environment or else youll get errors, we use a legacy version of numpy (it will complain about dtype ill bet)
Must have lsqlin.py in same directory as this file or program will error out

Assumptions I make that pavel doesn't necessarily
timeInter=failtime-10
timeInter=2000
"""
################# TODO #################
#sequences should be pre-aligned and split
    #if sequences are not aligned, may produce murky results
    #if sequences are of different length from one another, the program will error out
    #sho
#figure out optimal container type for seqs

#the tree variable - integral to the program, and now documented!
#  1  2 | 3
#| 0  0 | 0 |
#| 0  0 | 1 |
#| 0  0 | 2 |
#| 0  0 | 3 |
#-------|----Above this line represents leaf nodes, below represents internal nodes (above the line never changes, below line is variable)
#| 0  1 | 4 |
#| 2  3 | 5 |
#| 4  5 | 6 |
#_______|___
#       |#this side dispalys the parent of the two child nodes
#this side idsplays the child nodes
#in this case, numBranches=5,numLeaves=4,numLabels=9
#this is how this matrix would be represented as a tree
#     6
#   /   \
#  4     5   
# / \   / \
#0   1 2   3 
#typically, the trees will not be in order like this
# AM -> adjacency matrix
# DM -> distance matrix

def prealign(f1,f2):
    if verbose:
        print('aligning: '+f1+', '+f2)
    uid=re.findall('([^_]*)_.*',os.path.basename(f1))[0]
    with NamedTemporaryFile(delete=False) as tmpcat:
        catname=tmpcat.name
        try:
            subprocess.check_call(['cat',f1,f2],stdout=tmpcat)
        except:
            print('Invalid file. Check -s argument if in csv mode')
            return 0
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
    return(f1name,f2name)

# def prealign(f1,f2):
    # uid1=re.findall('([^_]*)_.*',f1)[0]
    # uid2=re.findall('([^_]*)_.*',f2)[0]
    # print(uid1,uid2)
    # with NamedTemporaryFile(delete=False) as tmpcat:
        # catname=tmpcat.name
        # try:
            # subprocess.check_call(['cat',f1,f2],stdout=tmpcat)
        # except:
            # print("Invalid file. Check -s argument if in csv mode")
            # return 0
    # with NamedTemporaryFile(delete=False) as aligned:
        # alignname=aligned.name
        # subprocess.check_call(['mafft', '--quiet', '--auto', '--thread', '20', '--preservecase', catname], stdout=aligned) 
    # os.unlink(catname)
    # seqs=SeqIO.parse(alignname,'fasta')
    # current=uid1
    # fnameIt=0
    # seenSecond=False
    # seqs1=[]
    # while not seenSecond:
        # record=next(seqs)
        # print("first")
        # thisseq=re.findall('([^_]*)_.*',record.id)[0]
        # if thisseq!=uid1:
            # seenSecond=True
        # else:
            # seqs1.append(record)
    # with NamedTemporaryFile(delete=False) as align1:
        # SeqIO.write(seqs1,align1,'fasta')
        # f1name=align1.name
    # seqs2=[]
    # seqs2.append(record)#write record already accessed and ignored because it doesn't belong in first file
    # co2=1
    # done=False
    # while not done:     
        # try:
            # record=next(seqs)
            # print("second")
            # seqs2.append(record)
        # except StopIteration:
            # done=True
    # with NamedTemporaryFile(delete=False) as align2:
        # SeqIO.write(seqs2,align2,'fasta')
        # f2name=align2.name
    # os.unlink(alignname)    
    # if verbose:
        # print("%s have been aligned" %", and ".join([f1,f2]))
    # return(f1name,f2name)

#@profile
def parseInput(input1,input2): #get sequences from 2 files
    seqs=[]
    with open(input1) as input_handle:
        for record in SeqIO.parse(input_handle, "fasta"): # for FASTQ use "fastq", for fasta "fasta"
            if len(record.seq) > 0 and len(record.seq) < 50000:
                seqs.append(record.seq)
    f1HaploNum=len(seqs)
    skippedSeqs=0
    with open(input2,"rU") as input_handle2:
        for record in SeqIO.parse(input_handle2, "fasta"): # for FASTQ use "fastq", for fasta "fasta"
            if len(record.seq) > 0 and len(record.seq) < 50000:
                if record.seq not in seqs:
                    seqs.append(record.seq)
                else:
                    if skippedSeqs==0 and verbose:
                        print("At least one of the sequences in f1 is identical to a sequence in f2; skipping")
                    skippedSeqs+=1
    if skippedSeqs>0 and verbose:
        print("%i total sequences were skipped" % skippedSeqs)
    haploSize = len(seqs[0])
    for seq in seqs:
        if len(seq)!=haploSize:
            sys.exit("Your files are not aligned! Try again with -a")
    haploNum = len(seqs)
    return(haploNum,haploSize,seqs,f1HaploNum)

#@profile 
def calcMapVal(seqs): #calculates amount of heterogeneity in sample, used in evolution simulation
    #only called once
    numSeq=len(seqs)
    numPos=len(seqs[0])
    mapval=0;a=0;c=0;g=0;t=0;blank=0
    for pos in range(numPos,0,-1):
        for seq in seqs:
            try:
                if seq[pos-1].lower()=="a":
                    a+=1
                elif seq[pos-1].lower()=="c":
                    c+=1
                elif seq[pos-1].lower()=="t":
                    t+=1
                elif seq[pos-1].lower()=="g":
                    g+=1
                else:
                    blank+=1
            except IndexError:
                sys.exit("It looks like your data are not aligned - try again with -a!")
        consensus=max(a,c,t,g,blank)
        if a+c+t+g+blank<numSeq or consensus>=numSeq or blank>0:
            a=0;c=0;g=0;t=0;blank=0;
        else:
            mapval+=1;a=0;c=0;g=0;t=0;blank=0;
    return mapval
    
#@profile
def calcOrderedFrequencies(haploNum,haploSize,seqs,byFreq): #split into 2 programs? if byFreq==True, seqs must be dict[seq]=freq but if byFreq==False, seqs must be list
    freqCount = np.zeros((haploSize, 5))
    productVector = np.zeros((haploSize, 5))
    try:
        if byFreq:
            total_reads=0
            for read in seqs:
                freq=seqs[read]
                total_reads+=freq
                for pos in range(haploSize):
                    if read[pos] == 'A':
                        freqCount[pos, 0] = freqCount[pos, 0] + freq
                    elif read[pos] == 'C':
                        freqCount[pos, 1] = freqCount[pos, 1] + freq
                    elif read[pos] == 'G':
                        freqCount[pos, 2] = freqCount[pos, 2] + freq
                    elif read[pos] == 'T':
                        freqCount[pos, 3] = freqCount[pos, 3] + freq
                    elif read[pos] == '-':
                        freqCount[pos, 4] = freqCount[pos, 4] + freq
            freqRel = np.divide(freqCount, total_reads, dtype = float)
        else:
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
    except IndexError:
        print("Your files are not aligned and it caused an error! Try again with -a")
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
    ordSeqs = []
    for i in seqs:
        newOne = ''
        for p in range(haploSize):
            newOne = ''.join([newOne, i[ordH[p]]])
        ordSeqs.append(newOne)
    return ordSeqs

#@profile
def calcKstep(haploNum,haploSize,ordSeqs): #Calculate hamming distances between sequences
    compNum = haploSize
    compList = range(haploNum)
    t = 0
    adjMatrix = np.zeros((haploNum, haploNum))
    kStepList = []
    while compNum > 1:
        t = t + 1
        for r1 in range(haploNum-1):
            haplotype1 = ordSeqs[r1]
            for r2 in range(r1+1, haploNum):
                if compList[r1] != compList[r2]: 
                    haplotype2 = ordSeqs[r2]
                    tempDist = 0
                    for a, b in itertools.izip(haplotype1, haplotype2):
                        if a != b:
                            tempDist = tempDist + 1
                            if tempDist > t:
                                break
                    if tempDist == t: 
                        adjMatrix[r1][r2] = 1
                        kStepList.append([r1, r2, t])
        sparseMatrix = csgraph_from_dense(adjMatrix)
        connected = connected_components(sparseMatrix, directed=False, connection='weak', return_labels=True)
        compNum = connected[0]
        compList = connected[1]
    return kStepList

#@profile
def getIntermediateSequences(kStepList,seqs,haploSize,verbose): #make 1-step list of sequences based off k-step list, return list of sequences
    already=0
    fakeSeqs={}
    for element in kStepList:
        dist=element[2]
        if dist!=1:
            #element[1] -> id for seq1
            #seqs[element[1]] -> sequence for seq1
            s1=seqs[element[0]]
            s2=seqs[element[1]]
            #assume sequences are same length
            id=0
            while dist>1:
                id+=1
                for nucl in range(haploSize):
                    if s1[nucl]!=s2[nucl]:
                        takeseq1=s1[:nucl+1]+s2[nucl+1:]
                        takeseq2=s1[:nucl]+s2[nucl:]
                        if takeseq1!=s1:
                            newseq=takeseq1
                        else:
                            newseq=takeseq2
                dist-=1
                s1=newseq
                if newseq in seqs:
                    if verbose:
                        print("Intermediate sequence detected in input, skipping")
                elif newseq in fakeSeqs.values():
                    already+=1
                else:
                    newseqid=str(element[0])+"->"+str(element[1])+"_"+str(id)
                    fakeSeqs[newseqid]=newseq
    if verbose:
        print("%i intermediate sequences had already been added to set and were not added again" %already)
    #make sequence dictionary
    finalSeqs=[]
    for item in range(len(seqs)):
        finalSeqs.append(seqs[item])
    finalSeqs.extend(fakeSeqs.values())
    return finalSeqs

#@profile
def calcDistanceMatrix(finalSeqs): 
    l=len(finalSeqs)
    arr=np.zeros([l,l])
    hdist=hamming(finalSeqs,finalSeqs,ignore_gaps=False)
    for id in range(len(hdist)):
        item=hdist[id]
        arr[:,id]=item[:,0]
    return arr
    
#@profile
def cluster(Z): #constructs clusters from a hierarchical cluster tree
    #custom implementation of matlabs CLUSTER
    #flatten, flattenhoriz and clusterNum are all helper functions for this function
    maxclust=2
    m=len(Z)+1
    T=np.zeros(m)
    if m<=maxclust:
        T=np.arange(m)+1
    else:
        clsnum=1
        k=m-1
        i=int(Z[k-1,0])
        if i<=m-1:
            T[i]=clsnum
            clsnum+=1
        elif i<(2*m-maxclust+1):
            T=clusterNum(Z,T,int(i-m+1),clsnum)#assign leaves under cluster clsnum to clsnum
            clsnum+=1
        i=Z[k-1,1]
        if i<=m-1:
            T[i]=clsnum
            clsnum+=1
        elif i<(2*m-maxclust+1):
            T=clusterNum(Z,T,int(i-m+1),clsnum)#assign leaves under cluster clsnum to clsnum
            clsnum+=1
    return T

#@profile    
def flatten(inarr): #mimics behavior of out=in(:) as out=flatten(in)
    rows=len(inarr)
    try:
        cols=len(inarr[0])
    except TypeError:
        out=inarr
        return out
    out=np.zeros(rows*cols)
    it=-1
    for col in range(cols):
        for row in range(rows):
            it+=1
            out[it]=inarr[row][col]
    return out

#@profile
def flattenHoriz(inarr): #more robust than calling flatten on the transpose?
    rows=len(inarr)
    try:
        cols=len(inarr[0])
    except TypeError:
        out=inarr
        return out
    out=np.zeros(rows*cols)
    it=-1
    for row in range(rows):
        for col in range(cols):
            it+=1
            out[it]=inarr[row][col]
    return out

#@profile    
def clusterNum(X,T,k,c): #assign leaves under cluster c to c
    X=X.astype(int)
    T=T.astype(int)
    m=len(X)+1
    while type(k)==int or type(k)==np.float64 or len(k)>0:
        if isinstance(k,np.ndarray):
            k=k.astype(int)
        children=X[k-1,0:2]
        children=flatten(children)
        t=children<m
        if np.size(children[t])!=0:
            try:
                T[children[t]]=c
            except IndexError:
                ele=int(children[t][0])
                T[ele]=c
        k=children[~t]-m
        if len(k)<2 or type(k)==int or type(k)==float:
            if k==0 or k==[0]:
                return T
    return T

#@profile
def estimateHubs(TransNet): #estimates how many hubs we have in the network
    AM=(TransNet+np.transpose(TransNet))>0
    deg=sum(AM)
    deg=sorted(deg,reverse=True)
    transposed=[]
    for item in deg:
        transposed.append([item])
    dists=pdist(transposed)
    linkage_out=fastcluster.linkage(dists,method='single')
    clustered=cluster(linkage_out)
    c=clustered[0]
    k=len(np.where(clustered==c))
    return k

#@profile
def sCore(n,k): #no idea
    n=float(n)
    k=float(k)
    d1=(n-k)/k+k-1
    d2=(n-k)/k+1
    return (n-k)/k*d1+(k-1)*(n-k)/k*d2+(k-1)*d1*d2

#@profile
def objTransNetQuadDeg(DM,nhubs): #calculates prior probabilty of a possible transmission tree using an s-metric parameter
    n=len(DM)
    AM=np.zeros([n,n])
    for i in range(n):
        for j in range(n):
            if DM[i][j]+DM[j][i]>0:
                AM[i][j]=1
    degseq=np.zeros(n)
    for id in range(len(AM)):
        degseq[id]=np.sum(AM[id])
    val=0
    for i in range(n):
        for j in range(i+1,n):
            if AM[i][j]==1:
                val=val+degseq[i]*degseq[j]
    s=sCore(n,nhubs)
    val=math.exp(-abs(val-s)/s)
    return val

#@profile
def modifyTree(tree,steps,nSamp): #central to the whole altorithm: modifies nonzero internal node children stochastically (bottom left portion) 
    # hasZero=np.any(tree[nSamp:len(tree),0:2]==0)
    # if hasZero: #i think this is always true but im not 100% so i left it commented out
    tree[nSamp:len(tree),0:2]+=1
    for i in range(steps):
        u=random.randint(nSamp,nSamp*2-2)
        j=random.randint(0,1)
        jinv=1-j
        v=int(tree[u][j])-1
        newtree=tree
        k=random.randint(0,1)
        w=newtree[u][jinv]
        if newtree[v,k]==0:
            continue
        newtree[u][jinv]=newtree[v][k]
        newtree[v][k]=w
        tree=newtree
    # if hasZero:
    tree[nSamp:len(tree),0:2]-=1
    return tree

#@profile
def findDSampDir(DSamp): #returns similar DSamp in which DSamp(i,j)==0 if DSamp(i,j)>DSamp(j,i)
    bools=DSamp<=np.transpose(DSamp)
    AMSamp_dir=bools.astype(int)
    return DSamp*AMSamp_dir

#@profile
def phyTree(linkage_out,nSamp):#custom implementation of matlabs phytree
    #creates an ultrametric phylogenetic tree object
    numLeaves=nSamp
    numBranches=nSamp-1
    numLabels=numLeaves+numBranches
    B=linkage_out[:,0:2]
    B=B.astype(int)
    D=np.zeros([numLabels,1])
    names=[]
    for id in range(numLeaves):
        names.append("Leaf "+str(id))
    for id in range(numBranches):
        names.append('Branch '+str(id))
        D[id+numLeaves]=linkage_out[id,2]
    tmp1=np.squeeze(D[B])
    tmp2=np.squeeze(D[numLeaves+(np.arange(numBranches))*[[1],[1]]])-np.transpose(tmp1)
    t2=flatten(tmp2)
    b=flattenHoriz(B)
    if len(tmp2.shape)==1:
        for ind in range(numLabels-1):
            D[ind]=tmp2[ind]
    elif len(tmp2.shape)==2:
        t2=flatten(tmp2)
        b=flattenHoriz(B)
        for ind in range(numLabels-1):
            id=int(b[ind])
            D[ind]=t2[id]
    else:
        raise ValueError('tmp2 is neither 1 nor 2 dimensional, abort')
    D[-1]=0
    needsReorder=False
    mi=np.zeros(numLabels)
    ma=np.zeros(numLabels)
    sz=np.ones(numLabels)
    for ind in range(numLeaves):
        mi[ind]=ind
        ma[ind]=ind
    for ind in range(numBranches):
        mi[ind+numLeaves]=0
        ma[ind+numLeaves]=0 
    for i in range(numBranches):
        v=B[i,:]
        j=i+numLeaves
        mi[j]=min(mi[v])
        ma[j]=max(ma[v])
        sz[j]=sum(sz[v])
        if ma[j]-mi[j]+1!=sz[j]:
            needsReorder=True
    if needsReorder:
        [B,names]=prettyOrder(B,D,names,numBranches,numLeaves,numLabels)
    return(B,names)

#@profile
def calcLabelsPhyloTree(tree,DM,nSamp,failtime): #calculates values for bottom right portion of tree (internal nodes' parents)
    nVert=len(tree)
    s=np.zeros(nVert-1)
    t=np.zeros(nVert-1)
    e=0
    G=nx.Graph()
    for i in range(nVert): 
        G.add_node(i)
        if tree[i,0]>0:
            s[e]=i
            t[e]=tree[i,0]# +-1's were here
            e=e+1
        if tree[i,1]>0:
            s[e]=i
            t[e]=tree[i,1]# +-1's were here
            e=e+1
    for i in range(len(s)):
        G.add_edge(t[i],s[i])
    backOrder=[]
    for node in nx.dfs_preorder_nodes(G): #does this do the same depth first search as matlab's dfsearch? probably.
        backOrder.append(node)
    order=[]
    for ele in reversed(backOrder):
        order.append(ele)
    internal=[]
    for ele in range(nVert):
        zero=tree[ele,0]>0
        one=tree[ele,1]>0
        internal.append(any([zero,one]))
    w=1
    for i in range(nVert):
        v=int(order[i])
        if internal[v]==True:
            child1=int(tree[v,0])
            child2=int(tree[v,1])
            l1=int(tree[child1,2])
            l2=int(tree[child2,2])
            if DM[l1,l2]<DM[l2,l1]:
                tree[v,2]=l1
            else:
                tree[v,2]=l2
            if DM[l1,l2]==failtime and DM[l2,l1]==failtime:
                w=0
    return(tree,w)

#@profile
def prettyOrder(B,D,names,numBranches,numLeaves,numLabels): #reorders the leaf nodes to avoid branch crossings
    #custom implementation of matlabs PRETTYORDER
    L=np.zeros(numLabels)
    X=np.zeros(numLabels)
    for a in range(numLeaves):
        L[a]=1
    for ind in range(numBranches):
        v=B[ind,:]
        L[ind+numLeaves]=sum(L[v])
    for ind in range(numBranches-1,-1,-1):
        v=B[ind,:]
        X[v]=D[v]+X[ind+numLeaves]
    Li=np.zeros(numLabels)
    Ls=np.zeros(numLabels)
    Ls[-1]=numLeaves
    for ind in range(numBranches-1,-1,-1): 
        v=B[ind,:]
        Ls[v]=Ls[ind+numLeaves]
        Li[v]=Li[ind+numLeaves]
        if X[v[1]]-X[v[0]]>=0:
            Ls[B[ind,0]]=Li[B[ind,0]]+L[B[ind,0]]
            Li[B[ind,1]]=Ls[B[ind,1]]-L[B[ind,1]]
        else:
            Ls[B[ind,1]]=Li[B[ind,1]]+L[B[ind,1]]
            Li[B[ind,0]]=Ls[B[ind,0]]-L[B[ind,0]]
    Ls=Ls.astype(int)
   
    qw=np.arange(numLeaves,numLabels)
    for ind in range(numLeaves,numLabels):
        Ls[ind]=ind+1
    namesFin=[]
    for (a,b) in sorted(zip(Ls,names)):
        namesFin.append(b)
    B=Ls[B]
    B=B-1
    return[B,namesFin]

#@profile
def treeToTransNet(tree,DSamp,nSamp): #Turns tree into transmission network
    nTreeVert=2*nSamp-1
    TransNet=np.zeros([nSamp,nSamp])
    for i in range(nTreeVert):
        l=int(tree[i,2])
        child1=int(tree[i,0])
        child2=int(tree[i,1])
        if child1+child2!=0:
            l1=int(tree[child1,2])
            l2=int(tree[child2,2])
            if l==l1:
                TransNet[l,l2]=DSamp[l,l2]
            if l==l2:
                TransNet[l,l1]=DSamp[l,l1]
    return TransNet  

#@profile    
def treeAM(tree): #creates adjacency matrix from tree 
    #AMtree===numLabels,numLabels structure with {0,1} as entries which represents same tree as 'tree'
    numLabels=len(tree)
    n=len(tree)
    AMtree=np.zeros([n,n])
    for i in range(n):
        child1=int(tree[i,0])
        child2=int(tree[i,1])
        if child1+child2>0:
            AMtree[i,child1]=1
            AMtree[i,child2]=1
            AMtree[child1,i]=1
            AMtree[child2,i]=1
    return AMtree

#@profile    
def graphShortestPath(dijkstra_output,source,target): #walks through output of dijkstra algorithm to find shortest path in unweighted, undirected graph 
    numNodes=len(dijkstra_output)
    if source<0 or source>=numNodes or target<0 or target>=numNodes:
        raise ValueError("Source and target must be integers within [0,numNodes-1]")
    theList=[source]
    dist=int(dijkstra_output[source,target])
    for nodesVisited in range(dist):
        sourceRow=dijkstra_output[source,:]
        neighbors=np.where(sourceRow==1) 
        if len(neighbors[0])==1:
            source=int(neighbors[0])
            theList.append(source)
            if source==target:
                return theList
        else:
            possibilities=[]
            for item in neighbors[0]:
                if item==target:
                    theList.append(target)
                    return theList
                possibilities.append(dijkstra_output[target,item])
            tmp=np.argmin(possibilities)
            source=neighbors[0][tmp]
            if dijkstra_output[source,target]==1:
                theList.append(source)
                theList.append(target)
                return theList
            else:
                theList.append(source)

#@profile
def objTransNetPhyloFit(AMtree,tree,DM): #calculates likelihood of given genetic distance given a possible transmission tree using a least-square approach
    rec=0;n=len(DM);am=len(AMtree);E=[]; nPairs=0
    for row in range(am):
        for col in range(am):
            if row>=col and AMtree[row,col]!=0:
                E.append([col,row])
    nedges=len(E)
    for row in range(n):
        for col in range(n):
            if DM[row,col]!=0:
                nPairs+=1
    C=np.zeros([nPairs,nedges])
    d=np.zeros([nPairs]) 
    pathsPairs=[]; p=-1
    for i in range(n):
        for j in range(n):
            if DM[i,j]>0:
                d[p]=1
                p+=1
                sparse_amtree=csgraph_from_dense(AMtree)
                distances=dijkstra(sparse_amtree)
                path=graphShortestPath(distances,i,j)
                pathsPairs.append(path)
                for v in range(len(path)-1):
                    if path[v]<path[v+1]:
                        edge=[path[v],path[v+1]]
                    else:
                        edge=[path[v+1],path[v]]
                    for ind in range(len(E)):
                        ele=E[ind]
                        if edge==ele:
                            C[p,ind]=1/(float(DM[i,j]))
    deg=np.sum(AMtree,axis=0)
    root_tmp=np.where(deg==2)
    leaves_tmp=np.where(deg==1)
    root=int(root_tmp[0])
    leaves=leaves_tmp[0]
    nleaves=len(leaves)
    nbranches=nleaves-1
    treeLen=nleaves+nbranches
    charVectPaths=np.zeros([nleaves,nedges])
    for leaf in leaves:
        path=graphShortestPath(distances,root,leaf)
        for v in range(len(path)-1):
            if path[v]<path[v+1]:
                edge=[path[v],path[v+1]]
            else:
                edge=[path[v+1],path[v]]
            for ind in range(len(E)):
                ele=E[ind]
                if edge==ele:
                    charVectPaths[leaf,ind]=1
    Aeq=np.zeros([nbranches,nedges])
    beq=np.zeros([nbranches,1])
    for i in range(nbranches):
        Aeq[i]=charVectPaths[i]-charVectPaths[i+1]
    lsq=lsqlin.lsqlin(C,d,0,None,None,Aeq,beq,0)
    x=lsqlin.cvxopt_to_numpy_matrix(lsq['x'])
    DSamp_tree=np.zeros([n,n])
    pathID=-1; vec1id=-1
    vec1=np.zeros(n*n)
    vec2=np.zeros(n*n)
    for i in range(n):
        for j in range(n):
            vec1id+=1
            vec1[vec1id]=DM[j][i] #supposed to be analagous to vec1=DM(:);
            if DM[i][j]>0:
                pathID+=1
                pathEdges=np.where(C[pathID]>0)
                DSamp_tree[i][j]=np.sum(x[pathEdges])
    vec2id=-1
    for i in range(n):
        for j in range(n):
            vec2id+=1
            vec2[vec2id]=DSamp_tree[j][i]
    ind=np.where(vec1>0)
    vec1=vec1[ind]
    vec2=vec2[ind]
    val=np.corrcoef(vec1,vec2)
    if val.ndim==2:
        v=val[0][1]
        if v>rec:
            rec=v
    return rec

#@profile    
def simulEvol(D_all,nseq1,nseq2,len_eff,percent,failtime): #viral evolution simulator
    ev=nseq2*float(percent)/100.0
    maxPopl=10**12
    timeInter=failtime-10 # -10 is somewhat arbitrary
    mutprob=.01
    evolved=False
    sz=len(D_all);
    Q=np.zeros([sz,sz],float)
    for u in range(sz):
        Q[u,u]=(mutprob/3)**D_all[u,u]*(1-mutprob)**(len_eff-D_all[u,u])            #modifying this line? vvvv
        for v in range(u+1,sz):
            Q[u,v]=D_all[u,v]*(mutprob/3)**(D_all[u,v])*(1-mutprob)**(len_eff-D_all[u,v])
            Q[v,u]=Q[u,v]
    x=np.zeros([sz,timeInter],float)
    x[0:nseq1,0]=1                                                              #could we replace this line by ^^^^
    E=np.eye(sz)
    for t in range(1,timeInter):
        x[:,t]=(1-sum(x[:,t-1])/maxPopl) * np.dot((E+Q),x[:,t-1])
        where=np.where(x[nseq1:nseq1+nseq2,t]>=1)
        if len(where[0])>=ev:
            evolved=True
            break
    if evolved==True:
        time=t+1
    if evolved==False:
        time=failtime
    return time

#@profile
def reduceMat(DSamp,comp): #probably get rid of this program, only called once and is simple
    #replaces matlab's smallMat=originalMat(smallvector,smallvector)
    #with smallMat=reduceMat(originalMat,smallvector)
    DSsz=len(DSamp)
    cosz=len(comp)
    out=np.zeros([cosz,cosz])
    for i in range(cosz):
        for j in range(cosz):
            out[i,j]=DSamp[comp[i],comp[j]]
    return out

#@profile    
def checkDSamp(DSamp): #check to see if any values of DSamp are below failtime
    l=len(DSamp)
    for row in range(l):
        for col in range(l):
            if row!=col:
                if DSamp[row,col]!=failtime:
                    return True
    return False

#@profile        
def findTransNetMCMC(DSamp_mcmc,failtime): #this is the second 'wrapper' part of the program - iteratively constructs most likely transmission network
    nSamp=len(DSamp_mcmc)
    nVertTree=2*nSamp-1
    DSamp_dir=findDSampDir(DSamp_mcmc)
    DSamp_symm=np.zeros([nSamp,nSamp])
    for u in range (nSamp):
        for v in range(u+1,nSamp):
            if DSamp_mcmc[u,v]<=DSamp_mcmc[v,u]:
                DSamp_symm[u,v]=DSamp_mcmc[u,v]
                DSamp_symm[v,u]=DSamp_mcmc[u,v]
            else:
                DSamp_symm[u,v]=DSamp_mcmc[v,u]
                DSamp_symm[v,u]=DSamp_mcmc[v,u]
    l=fastcluster.linkage(DSamp_symm,method='average') #~10x faster than scipy's linkage
    #linkage function: builds hierarchical clustering on distance matrix 
    [branchLists,nodeNames]=phyTree(l,nSamp)
    leafID=np.zeros(nSamp)
    tree=np.zeros([nVertTree,3])
    for i in range(nSamp):
        name=nodeNames[i]
        leafID[i]=int(name[5:])
        tree[i,2]=i
    for i in range(nSamp-1):
        for j in range(2):
            if branchLists[i,j]<nSamp:
                tree[nSamp+i,j]=leafID[branchLists[i,j]] # +-1's were here
            else:
                tree[nSamp+i,j]=branchLists[i,j] # +-1's were here
    tree,w=calcLabelsPhyloTree(tree,DSamp_mcmc,nSamp,failtime)
    TransNet=treeToTransNet(tree,DSamp_mcmc,nSamp)
    k=estimateHubs(TransNet)
    AMtree=treeAM(tree)
    matrToFit=DSamp_dir
    obj1=objTransNetPhyloFit(AMtree,tree,matrToFit)
    obj2LinkNoPow=objTransNetQuadDeg(TransNet,k)
    obj2=obj2LinkNoPow
    record=w*obj1*obj2
    TransNetTrees=[]
    TransNetTrees.append(tree)
    delta=.0001;iEdgeMotif=1;nEdgeMotifCurr=2;nIter=250
    i=-1
    while i<=nIter*2:
        i+=1
        if iEdgeMotif==nIter+1:
            iEdgeMotif=1
            nEdgeMotifCurr-=1
        newtree=modifyTree(tree,nEdgeMotifCurr,nSamp)
        if np.all(newtree==0):
            continue
        AM_newtree=treeAM(newtree)
        sparse_amtree=csgraph_from_dense(AM_newtree)
        distances=dijkstra(sparse_amtree)
        newtree,w = calcLabelsPhyloTree(newtree,DSamp_mcmc,nSamp,failtime)
        TransNetNew=treeToTransNet(newtree,DSamp_mcmc,nSamp)
        matrToFit=DSamp_dir
        obj1=objTransNetPhyloFit(AM_newtree,newtree,matrToFit)
        obj2=objTransNetQuadDeg(TransNetNew,k)
        val=w*obj1*obj2
        try: #record is often zero
            prob=float(val)/float(record)
        except ZeroDivisionError:
            prob=0
        if 1<=prob:
            tree=newtree
        if val > record+delta: #I have never seen this happen
            record=val
            TransNetTrees.append(newtree)
        if abs(val-record)<delta:
            iscont=False
            nopt=len(TransNetTrees)
            for c in range(nopt):
                if np.all(np.equal(TransNetTrees[c],newtree)):
                    iscont=True
                    break
            if not iscont:
                TransNetTrees.append(newtree)
        iEdgeMotif+=1
    nopt=len(TransNetTrees)
    TransNets=[]
    for i in  range(nopt):
        try:
            tree=TransNetTrees[i]
        except UnboundLocalerror:
            tree=TransNetTrees
        AM_tree=treeToTransNet(tree,DSamp_mcmc,nSamp)
        isExist=False
        for j in range(len(TransNets)):
            if np.all(np.equal(TransNets[j],AM_tree)):
                isExist=True
                break
        if not isExist:
            TransNets.append(AM_tree)
    return TransNets

#@profile
def inferTransNetFromEvolTimes(DSamp,inputs):
    numFiles=len(DSamp)
    AMSamp=np.eye(numFiles)
    for u in range (numFiles):
        for v in range(u+1,numFiles):
            if DSamp[u,v]<=DSamp[v,u]:
                if DSamp[u,v]<failtime-20: 
                    AMSamp[u,v]=1
            else:
                if DSamp[v,u]<failtime-20: 
                    AMSamp[v,u]=1
    sparseMatrix = csgraph_from_dense(AMSamp)
    connected = connected_components(sparseMatrix, directed=False, connection='weak', return_labels=True)
    S = connected[0]
    C = connected[1]
    if S!=1:
        print("FATAL ERROR: Further analysis only makes sense for a connected component. Please rerun QUENTIN with the following arguments:")
        for a in range(S):
            runlist=[]
            for i in range(len(C)):
                if C[i]==a:
                    runlist.append(inputs[i])
            rerun="quentin -f"
            if len(runlist)>1:
                rerun="quentin -f"
                for item in runlist:
                    rerun=rerun+" "+item
                print(rerun)
            else:
                print("%s does not appear to be connected to any of the other files" % runlist[0])
        sys.exit("exiting QUENTIN...")
    transNets=[]
    for c in range(S):
        tmp=np.where(C==c)
        comp=tmp[0]
        if len(comp)>1:
            DSamp_comp=reduceMat(DSamp,comp)
            transNetsComp=findTransNetMCMC(DSamp_comp,failtime) 
            try:
                transNets.append(transNetsComp[0])
            except TypeError:
                transNets.append(transNetsComp)
    TransNet=transNets[0]
    return TransNet

def specialquentin2files(inputs,output,startTime):
    (forwards,backwards)=getEvolTimes(inputs[0],inputs[1])
    g=nx.DiGraph()
    if forwards>backwards:
        g.add_edge(inputs[1],inputs[0],label=backwards)
        print("If a transmission event occured between these two individuals, %s gave HCV to %s." %(inputs[1],inputs[0]))
    elif forwards<backwards:
        g.add_edge(inputs[0],inputs[1],label=forwards)
        print("If a transmission event occured between these two individuals, %s gave HCV to %s." %(inputs[0],inputs[1]))
    else:
        sys.exit("I cannot tell who the source is and who the target is! Exiting.")
    # with NamedTemporaryFile(delete=False) as tmpdot:
        # tmpname=tmpdot.name
        # nx.drawing.nx_agraph.write_dot(g,tmpname)
    # subprocess.check_call(['dot','-Tpng',tmpname,'-o',output+'.png'])
    # os.unlink(tmpname)
    endTime=time.time()
    workTime=endTime-startTime
    statement_time='Analysis took %.3f seconds' % workTime
    print(statement_time)
    sys.exit("Program finished!")

def trimfh(handle):
    return os.path.splitext(os.path.basename(handle))[0]
    
def matrixToWeightedNx(transNet,fileList): #go from a transNet to a weighted networkx graph
    g=nx.DiGraph()
    numRows=len(transNet)
    numCols=len(transNet[0])
    for row in range(numRows):
        for col in range(numCols):
            if row!=col:
                if transNet[row,col]!=0:
                    shortname=trimfh(fileList[row])
                    g.add_edge(shortname,fileList[col],label=transNet[row,col])
    return g

def calcEntropies(inputs): 
    l=len(inputs)
    entropies={}
    for i in range(l):
        seqs={}
        file=inputs[i]
        with open(file,"rU") as input_handle:
            for record in SeqIO.parse(input_handle, "fasta"):
                if len(record.seq) > 0 and len(record.seq) < 50000:
                    freq=re.findall('_(\d*)$',record.id)[0]
                    seqs[record.seq]=int(freq)
        haploNum=len(seqs)
        haploSize=len(six.next(six.iterkeys(seqs)))
        entropy=calcOrderedFrequencies(haploNum,haploSize,seqs,True)
        ent=sum(entropy)/len(entropy)
        entropies[file]=ent
    return entropies
    
def graphMeWithGV(TransNet,fileList,output,sourceNotPresent):
    try:
        shortNames=[]
        for item in fileList:
            shortNames.append(trimfh(item))
        gr=matrixToWeightedNx(TransNet,shortNames)
    except IndexError:
        gr=matrixToWeightedNx(TransNet,fileList)
    with NamedTemporaryFile(delete=False) as tmpdot:
        tmpdotname=tmpdot.name
        nx.drawing.nx_agraph.write_dot(gr,tmpdotname)
    if sourceNotPresent==1:
        lineco=0
        with open(tmpdotname,'r') as f:
            lines=f.readlines()
        l=len(lines)
        with NamedTemporaryFile(delete=False) as outdot:
            outdotname=outdot.name
            for line in lines:
                lineco+=1
                if lineco!=l:
                    outdot.write(line)
                else:
                    outdot.write("{\n")
                    outdot.write("    Legend [shape=none, margin=0, colorscheme=set23 label=<\n")
                    outdot.write("""    <TABLE BORDER="0" CELLBORDER="0" CELLSPACING="0" CELLPADDING="4">\n""")
                    outdot.write("     <TR>\n")
                    outdot.write("""      <TD>The node at the top of this chart does not have the highest nucleotide diversity</TD>\n""")
                    outdot.write("     </TR>\n     <TR>\n")
                    outdot.write("""      <TD>This means that the source is very likely not present among these samples</TD>\n""")
                    outdot.write("     </TR>\n    </TABLE>\n   >];\n  }\n} ")
        subprocess.check_call(['dot','-Tpng',outdotname,'-o',output+'.png'])
        os.unlink(outdotname)
    elif sourceNotPresent==2:
        lineco=0
        with open(tmpdotname,'r') as f:
            lines=f.readlines()
        l=len(lines)
        with NamedTemporaryFile(delete=False) as outdot:
            outdotname=outdot.name
            for line in lines:
                lineco+=1
                if lineco!=l:
                    outdot.write(line)
                else:
                    outdot.write("{\n")
                    outdot.write("    Legend [shape=none, margin=0, colorscheme=set23 label=<\n")
                    outdot.write("""    <TABLE BORDER="0" CELLBORDER="0" CELLSPACING="0" CELLPADDING="4">\n""")
                    outdot.write("     <TR>\n")
                    outdot.write("""      <TD>The nucleotide diversity from the node at the top of this chart is not very high </TD>\n""")
                    outdot.write("     </TR>\n     <TR>\n")
                    outdot.write("""      <TD>This means that the source may not be present among these samples</TD>\n""")
                    outdot.write("     </TR>\n    </TABLE>\n   >];\n  }\n} ")
        subprocess.check_call(['dot','-Tpng',outdotname,'-o',output+'.png'])
        os.unlink(outdotname)
    else:
        subprocess.check_call(['dot','-Tpng',tmpdotname,'-o',output+'.png'])
    os.unlink(tmpdotname)

def getEvolTimes(input1,input2,**kwargs): #get evolution times both ways from fasta files
    try:
        ali
    except NameError:
        global ali
        ali=False
    try:
        verbose
    except NameError:
        global verbose
        verbose=False
    try:
        failtime
    except NameError:
        global failtime
        failtime=3000
    try:
        percent
    except NameError:
        global percent
        percent=100
    fileName=trimfh(input1)
    fileName2=trimfh(input2)
    print("Running on: %s" %", ".join([fileName,fileName2]))
    if ali:
        try:
            (f1name,f2name)=prealign(input1,input2)
        except TypeError:
            print("I'll bet that you are passing in malformed file names")
    else:
        f1name=input1
        f2name=input2
    (haploNum,haploSize,seqs,f1HaploNum)=parseInput(f1name,f2name)
    if verbose:
        print("simulating evolution from %s to %s" % (fileName,fileName2))
    forwards=getEvolTime(haploNum,haploSize,seqs,f1HaploNum,verbose,percent,failtime)
    print('Forwards evolution time: %i' % forwards)
    if verbose:
        print("simulating evolution from %s to %s" % (fileName2,fileName))
    (haploNum,haploSize,seqs,f1HaploNum)=parseInput(f2name,f1name)
    backwards=getEvolTime(haploNum,haploSize,seqs,f1HaploNum,verbose,percent,failtime)
    print('Backwards evolution time: %i' % backwards)
    if ali:
        os.unlink(f1name)
        os.unlink(f2name)
    return forwards,backwards
        
#@profile
def getEvolTime(haploNum,haploSize,seqs,f1HaploNum,verbose,percent,failtime): #get evolution time one way from already parsed sequences
    f2HaploNum=haploNum-f1HaploNum
    len_eff=calcMapVal(seqs) 
    hVector=calcOrderedFrequencies(haploNum,haploSize,seqs,False)
    ordSeqs=orderPositions(hVector,seqs,haploSize)
    kStepList=calcKstep(haploNum,haploSize,ordSeqs)
    # for item in kStepList:
        # print(item)
    # print(len(kStepList))
    
    finalSeqs=getIntermediateSequences(kStepList,seqs,haploSize,verbose)
    # print('33',f1HaploNum)
    # print('23',f2HaploNum)
    # count=0
    # for seq in finalSeqs:
        # if count<f1HaploNum:
            # print('>1_12')
        # elif count<haploNum:
            # print('>2_12')
        # else:
            # print('>fake_12')
        # print(seq)
        # count+=1
    # sys.exit()
    D_all=calcDistanceMatrix(finalSeqs)
    evolTime=simulEvol(D_all,f1HaploNum,f2HaploNum,len_eff,percent,failtime)
    return evolTime
    
#@profile
def main(inputs,output): #overall wrapper - wraps the wrappers, checks some things
    if verbose:
        print("----------------------------------------------")
        print("Running QUENTIN software on the following inputs:")
        print("----------------------------------------------")
        for input in inputs:
            if not input.endswith("fas") and not input.endswith('fasta') and not input.endswith('fa'):
                print(input)
                sys.exit("Warning! One of your files may not be in fasta format, this will almost certainly cause an error later")
            print(input)
        print("----------------------------------------------")
    else:
        options['show_progress'] = False
        for input in inputs:
            if not input.endswith("fas") and not input.endswith('fasta') and not input.endswith('fa'):
                print(input)
                sys.exit("Warning! One of your files may not be in fasta format, this will almost certainly cause an error later")
    startTime = time.time()
    numFiles=len(inputs)
    if numFiles==2:
        specialquentin2files(inputs,output,startTime)
    
    DSamp=np.zeros([numFiles,numFiles],int)
    for i1,i2 in itertools.combinations(range(len(inputs)),2):
        f1=inputs[i1]
        f2=inputs[i2]
        (forwards,backwards)=getEvolTimes(f1,f2)
        DSamp[i1,i2]=forwards
        DSamp[i2,i1]=backwards
    print(DSamp)
    for item in inputs:
        print(item)
    isCluster=checkDSamp(DSamp)
    if isCluster==False:
        print("None of your samples appear to be related! Exiting")
        print(DSamp)
        endTime = time.time()
        workTime =  endTime - startTime
        statement_time='Analysis took %.3f seconds' % workTime
        sys.exit(statement_time)
    if verbose:
        print("Evolution simulations complete. Now building optimal tree")
    TransNet=inferTransNetFromEvolTimes(DSamp,inputs)
    entropies=calcEntropies(inputs)
    sourceFromTransNet=inputs[int(np.where(sum(TransNet)==0)[0])]
    
    print('-')
    print("source here")
    print(sourceFromTransNet)
    print(max(entropies))
    print(min(entropies.values()))
    print(max(entropies.values()))
    
    for item in entropies:
        print(item,entropies[item])
    print('-')
    sourceNotPresent=0
    if sourceFromTransNet!=max(entropies):
        sourceNotPresent=2
        print("ERROR: QUENTIN source prediction does not posses highest nucleotide diversity. Output is not very reliable in this case") #not meeting this condition means the tree is really not reliable
    else:
        entVals=entropies.values()
        oldMaxZ=max(zscore(entVals))
        entVals.remove(max(entVals))
        newMaxZ=max(zscore(entVals))
        if newMaxZ>oldMaxZ*1.1:
            sourceNotPresent=1
            print("Warning: source prediction may not be reliable") #this isn't quite so bad as the other condition, show figure with warning 
    graphMeWithGV(TransNet,inputs,output,sourceNotPresent)
    endTime = time.time()
    workTime =  endTime - startTime
    statement_time='Analysis took %.3f seconds. Thank you for using QUENTIN!' % workTime
    print(statement_time)
    
def csvmode(input,output,st):
    print("csvmode engage: writing to %s.csv" %output)
    with open(input,'r') as f:
        lines=f.readlines()
    ############################################################################ regular mode
    with open(output+'.csv',"w") as g:
        for line in lines:
            [p1,p2,dist,gen1,tmp]=line.split(",")
            gen2=tmp.strip()
            f1=p1+"_"+gen1+st+".fas"
            f2=p2+"_"+gen2+st+".fas"
            if os.path.isfile(f1) and os.path.isfile(f2):
                (forwards,backwards)=getEvolTimes(f1,f2)
                if forwards<backwards:
                    outline=",".join([p1,p2,dist,gen1,gen2,str(forwards),str(backwards)])+"\n"
                else:
                    outline=",".join([p2,p1,dist,gen2,gen1,str(backwards),str(forwards)])+"\n"
                g.write(outline)
                print(outline)
    ############################################################################ sparse csv file mode
    # with open(output+'.csv',"w") as g:
        # for line in lines:
            # [f1,tmp]=line.split(",")
            # f2=tmp.strip()
            # (forwards,backwards)=getEvolTimes(f1,f2)
            # if forwards<backwards:
                # outline=",".join([f1,f2,str(forwards),str(backwards)])+"\n"
            # else:
                # outline=",".join([f2,f1,str(backwards),str(forwards)])+"\n"
            # g.write(outline)
            # print(outline)            
    ############################################################################ special mode for time determination experiments
    # fileList=[]
    # days={}
    # for line in lines:
        # [pid,time,fid,trash]=line.split(",")
        # fileName=fid+'_asdf.fas'
        # fileList.append(fileName)
        # days[fileName]=time
    # with open(output+'.csv',"w") as g:
        # for f1,f2 in itertools.combinations(fileList,2):
            # (forwards,backwards)=getEvolTimes(f1,f2)
            # if forwards<backwards:
                # time=int(days[f1])-int(days[f2])
                # g.write(','.join([pid,f1,f2,str(forwards),str(backwards),str(time)])+'\n')
            # else:
                # time=int(days[f2])-int(days[f1])
                # g.write(','.join([pid,f2,f1,str(backwards),str(forwards),str(time)])+'\n')
    ############################################################################
            
def mashmode(inputs,output):
    fileList=[]
    for file in inputs:
        if file.endswith("fas") or file.endswith("fasta") or file.endswith("fas"):
            fileList.append(file+".msh")
            wr="mash sketch "+file
            os.system(wr)
    fileit=iter(fileList)
    first=next(fileit)
    second=next(fileit)
    f=open('mashdists.txt','a')
    dist=" ".join(['mash','dist',first,second,'>','mashdists.txt'])
    os.system(dist)
    paste=" ".join(['mash','paste','tmpmash',first,second])
    os.system(paste)
    finished=False
    while not finished:
        try:
            newfile=next(fileit)
            dist=" ".join(['mash','dist','tmpmash.msh',newfile,'>','mashdists.txt','>>','mashdists.txt'])
            os.system(dist)
            paste=" ".join(['mash','paste','pastedmash',newfile,'tmpmash.msh'])
            os.system(paste)
            move=" ".join(['mv','pastedmash.msh','tmpmash.msh'])
            os.system(move)
        except StopIteration:
            finished=True
    os.system('rm *.msh')
    # os.system("sort -k 3,3 mashdists.txt > tmp.txt")
    # os.system("mv tmp.txt mashdists.txt")
    with open("mashdists.txt","r") as f:
        lines=f.readlines()
    os.remove("mashdists.txt")
    g=nx.Graph()
    for line in lines:  
        [f1,f2,dist,trash1,trash2]=line.split("\t")
        if float(dist)<1:
            (forwards,backwards)=getEvolTimes(f1,f2)
            if forwards!=failtime or backwards!=failtime:
                g.add_edge(f1,f2,label=min(forwards,backwards))
                print("edge added!",f1,f2)
    nx.drawing.nx_agraph.write_dot(g,output+'.dot')
    print(nx.number_connected_components(g))
    for item in nx.connected_components(g):
        print(item)
    subprocess.check_call(['dot','-Tpng',output+'.dot','-o',output+'.png'])
    subprocess.check_call(['eog',output+'.png'])
            
if __name__ == '__main__':
    import argparse # possible arguments to add: delta, nIter
    parser = argparse.ArgumentParser(description='QUENTIN: predict directionality of HCV infection')
    parser.add_argument('files',
        nargs='+',
        help="List  of  files  to  be  analyzed,  order  of  files  does  not  matter")
    parser.add_argument('-o','--output', 
        type=str, required=False, default="quentin_links", 
        help="The name of the output .png file")
    parser.add_argument('-v','--verbose', 
        action='store_true', required=False, default=False, 
        help="Add this flag to get more information printed out")
    parser.add_argument('-c','--csvmode', 
        action='store_true', required=False, default=False, 
        help="Add this flag to use CSV mode - take as input threshold_links.csv and output a similar csv file")
    parser.add_argument('-m','--mash', 
        action='store_true', required=False, default=False, 
        help="Add this flag to use mash exploration mode - use mash to find close pairs and then analyze with quentin over current settings.")
    parser.add_argument('-a',  '--align', 
        action='store_true',  default=False,
        help="Pass  this  as  an  argument  to  align  your  files. Highly reccomended if they are not already aligned")
    parser.add_argument('-t',  '--time', 
        type=int, required=False,  default=3000,
        help="Time allowed for file1 to evolve into file2")
    parser.add_argument('-p',  '--percent', 
        type=int, required=False,  default=100,
        help="Percent evolution needed for you to say that file1 has evolved into file2")
    parser.add_argument('-s','--string', 
        type=str, required=False, default="", 
        help="String you need to append to your threshold_links.csv to correctly match your file names")
    
    args = parser.parse_args()
    global ali
    ali=args.align
    output=args.output
    files=args.files
    global verbose
    verbose=args.verbose
    global failtime
    failtime=args.time
    global percent
    percent=args.percent
    st=args.string
    if args.csvmode==True and args.mash==True:
        sys.exit('Please choose mash mode or csv mode, not both')
    elif args.csvmode==True:
        csvmode(files[0],output,st)
    elif args.mash==True:
        mashmode(files,output)
    else:
        main(files,output)
