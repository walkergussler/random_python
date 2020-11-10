from __future__ import absolute_import
import numpy as np
import math,sys,random,fastcluster,os,re,time
import networkx as nx
import matplotlib.pyplot as plt
from scipy.sparse.csgraph import connected_components,csgraph_from_dense,dijkstra 
from scipy.spatial.distance import pdist
from itertools import izip
from Bio import SeqIO
from ghost.util.distance import hamming

__author__ = "CDC/OID/NCHHSTP/DVH bioinformatics team"

"""QUENTIN, by J Walker Gussler. Translated from https://github.com/skumsp/quentin by Pavel Skums. 
Given a number of .fas clinical HCV samples of HVR1, this program will calculate clusters (does an ok job)
and then it will calculate the directionality of infection (this is what the program is really for)
 
1. It requires as input the location of the aligned fasta files. All sequences are assumed to be different.
2. It calculates the k-step network using an iterative multi-index hashing approach.
3. It calculates all intermediate sequences so that we have a 1-step network (with some fabricated sequences) instead of a k-step network
4. From this 1-step network with intermediate sequences, we compute an all-against-all hamming distance matrix.
5. From this hamming distance matrix, simulEvol computes expected evolution times for file A to evolve into file B and vice versa
6. Within an input cluster of n files, a nxn evolution simulation matrix is constructed
7. From these simulated evolution times, we iteratively construct a phylogenetic tree and infer linkage data from the tree

Should be inside ghost virtual environment or else youll get errors, we use a legacy version of numpy
This program is robust at comparing files with slightly different sequence lengths (i.e. 264 vs 266)
However, it has not been tested with files that have different sequence lengths within the individual files themselves
Must have lsqlin.py in same directory or program will error out

To run program:
python kquentin.py <file1> <file2> ... <fileN>
Or if you're a fancy person who names files so that their cluster matches some regex:
python kquentin.py ./Related_clipped/AA*

Assumptions I make that pavel doesn't necessarily
timeInter=failtime-20
timeInter=2000
"""
################# TODO #################
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

#@profile
def parseInput(input1,input2): #get sequences from 2 files
    seqs=[]
    input_handle = open(input1) 
    input_handle2 = open(input2)
    for record in SeqIO.parse(input_handle, "fasta"): # for FASTQ use "fastq", for fasta "fasta"
        if len(record.seq) > 0 and len(record.seq) < 50000:
            seqs.append(record.seq)
    input_handle.close()
    f1HaploNum=len(seqs)
    for record in SeqIO.parse(input_handle2, "fasta"): # for FASTQ use "fastq", for fasta "fasta"
        if len(record.seq) > 0 and len(record.seq) < 50000:
            if record.seq not in seqs:
                seqs.append(record.seq)
            else:
                print("At least one of the sequences in f1 is identical to a sequence in f2; skipping")
    input_handle2.close()
    haploSize = len(seqs[0])
    for seq in seqs:
        if len(seq)!=haploSize:
            haploSize,seqs=alignFixEarly(input1,input2)
    haploNum = len(seqs)
    return(haploNum,haploSize,seqs,f1HaploNum)

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
def alignFixEarly(input1,input2):
    catstr="cat "+input1+" "+input2+" > tmp.fas"
    os.system(catstr)
    os.system('mafft --quiet --auto --thread 20 --preservecase tmp.fas > align.fas')
    outs=[]
    input_handle = open('align.fas')
    for record in SeqIO.parse(input_handle, "fasta"): # for FASTQ use "fastq", for fasta "fasta"
        if len(record.seq) > 0 and len(record.seq) < 50000:
            outs.append(record.seq)
    input_handle.close()
    rmstr="rm align.fas tmp.fas"
    os.system(rmstr)
    return len(outs[0]),outs

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
def calcDistances(haploNum,haploSize,ordSeqs): #Calculate hamming distances between sequences
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
                        kStepList.append([r1, r2, t])
        # Calculate components
        sparseMatrix = csgraph_from_dense(adjMatrix)
        connected = connected_components(sparseMatrix, directed=False, connection='weak', return_labels=True)
        compNum = connected[0]
        compList = connected[1]
    return kStepList

#@profile
def makeNewEdgeList(kStepList,seqs,haploSize): #make 1-step list of sequences based off k-step list, return list of sequences
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
                    print("Intermediate sequence detected in input, skipping")
                elif newseq in fakeSeqs.values():
                    already+=1
                else:
                    newseqid=str(element[0])+"->"+str(element[1])+"_"+str(id)
                    fakeSeqs[newseqid]=newseq
    print("%i intermediate sequences had already been added to set and were not added again" %already)
    #make sequence dictionary
    finalSeqs=[]
    for item in range(len(seqs)):
        finalSeqs.append(seqs[item])
    finalSeqs.extend(fakeSeqs.values())
    return finalSeqs

#@profile
def calcDistanceMatrix(finalSeqs): #calculate distance matrix from the 1-step list
    l=len(finalSeqs)
    arr=np.zeros([l,l])
    hdist=hamming(finalSeqs,finalSeqs,ignore_gaps=False)
    for id in range(len(hdist)):
        item=hdist[id]
        arr[:,id]=item[:,0]
    return arr
    
#@profile
def alignFix(seqs,input1,input2): # this program may be unnecessary, but it probably catches niche errors so im afraid to delete it
    catstr="cat "+input1+" "+input2+" > tmp.fas"
    os.system(catstr)
    os.system('mafft --quiet --auto --thread 20 --preservecase tmp.fas > align.fas')
    outs=[]
    input_handle = open('align.fas')
    for record in SeqIO.parse(input_handle, "fasta"): # for FASTQ use "fastq", for fasta "fasta"
        if len(record.seq) > 0 and len(record.seq) < 50000:
            outs.append(record.seq)
    input_handle.close()
    rmstr="rm align.fas tmp.fas"
    calcMapVal(outs,input1,input2)
 
#@profile 
def calcMapVal(seqs,input1,input2): #calculates amount of heterogeneity in sample, used in evolution simulation
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
                alignFix(seqs,input1,input2)
        consensus=max(a,c,t,g,blank)
        if a+c+t+g+blank<numSeq or consensus>=numSeq or blank>0:
            a=0;c=0;g=0;t=0;blank=0;
        else:
            mapval+=1;a=0;c=0;g=0;t=0;blank=0;
    return mapval

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
    hasZero=np.any(tree[nSamp:len(tree),0:2]==0)
    if hasZero:
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
    if hasZero:
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
    #check input
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
    for node in nx.dfs_preorder_nodes(G): #does this do the same depth first search as matlab's dfsearch?
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
def treeAM(tree): #creates amtree from tree (what is amtree? - numLabels,numLabels structure with {0,1} as entries which represents same tree as 'tree')
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
def graphShortestPath(dijkstra_output,source,target): #walks through output of dijkstra algorithm to find shortest path in unweighted, undirected graph - i think it works
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
def objTransNetPhyloFit(AMtree,tree,DM): #calculates likelihood of given genetic distance given a possible transmission tree using a least-square approac
    import lsqlin
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
    if val.ndim>1:
        for row in range(len(val)):
            for col in range(len(val)):
                if abs(val[row,col]-1)>.00005:
                    if val[row,col]>rec:
                        rec=val[row,col]
    return rec

#@profile    
def simulEvol(D_all,nseq1,nseq2,len_eff,failtime): #viral evolution simulator, seems to work 
    maxPopl=10**12
    timeInter=failtime-20 # -20 is somewhat arbitrary
    mutprob=.01
    evolved=False
    sz=len(D_all);
    Q=np.zeros([sz,sz],float)
    for u in range(sz):
        Q[u,u]=(mutprob/3)**D_all[u,u]*(1-mutprob)**(len_eff-D_all[u,u])
        for v in range(u+1,sz):
            Q[u,v]=D_all[u,v]*(mutprob/3)**(D_all[u,v])*(1-mutprob)**(len_eff-D_all[u,v])
            Q[v,u]=Q[u,v]
    x=np.zeros([sz,timeInter],float)
    x[0:nseq1-1,0]=1
    E=np.eye(sz)
    for t in range(1,timeInter):
        x[:,t]=(1-sum(x[:,t-1])/maxPopl) * np.dot((E+Q),x[:,t-1])
        where=np.where(x[nseq1+1:nseq1+nseq2,t]>=1)
        if len(where[0])==len(x[nseq1+1:nseq1+nseq2,t]):
            evolved=True
            break
    if evolved==True:
        time=t+1
    if evolved==False:
        time=failtime
    return (time)

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
def graphme(TransNet,fileList,output): #plots nx graph
    sparse_graph=csgraph_from_dense(TransNet)
    graph_obj=nx.from_scipy_sparse_matrix(sparse_graph,create_using=nx.DiGraph())
    nodeLabels={}
    edgeLabels={}
    for node in graph_obj.node:
    #------------------------------------------------------------------------------------------
    #this part generates names using some regex grab from your file names for cleaner charts
    #for no regex, use 
        nodeLabels[node]=fileList[node]
        #or match the unique part of your filename and use that instead of the full name (below)
        # name=re.findall('.*_clipped/(.*)_unique.*',fileList[node])
        # nodeLabels[node]=name
        # name=re.findall('(.*).*',fileList[node])
        # nodeLabels[node]=name
    #------------------------------------------------------------------------------------------
    edgeListIterator=0
    for edge in sparse_graph.data:
        edgeLabels[graph_obj.edges()[edgeListIterator]]=str(edge)
        edgeListIterator+=1
    pos=nx.shell_layout(graph_obj,dim=2) 
    nx.draw_networkx_nodes(graph_obj,pos,node_shape="s")
    nx.draw_networkx_edges(graph_obj,pos,arrows=True)
    nx.draw_networkx_labels(graph_obj,pos,labels=nodeLabels)
    nx.draw_networkx_edge_labels(graph_obj,pos,edge_labels=edgeLabels)
    plt.axis('off')
    plt.savefig(output)

#@profile    
def checkDSamp(DSamp,failtime): #check to see if any values of DSamp are below failtime
    l=len(DSamp)
    for row in range(l):
        for col in range(l):
            if row!=col:
                if DSamp[row,col]!=failtime:
                    return True
    return False

#@profile    
def getEvolTime(input1,input2,failtime): #wrapper for 'pairwise' portion - gets data ready for evolution simulation 
    fileName, fileExtension = os.path.splitext(input1)
    fileName2, fileExtension2 = os.path.splitext(input2)
    statement_f1='File 1: %s' %input1;statement_f2='File 2: %s' %input2;
    print("Parsing the following fasta files...")
    print(statement_f1)
    print(statement_f2)
    [haploNum,haploSize,seqs,f1HaploNum]=parseInput(input1,input2)
    #f1HaploNum -> # of seqeunces in file 1
    #haploNum -> # of sequences in both files
    f2HaploNum=haploNum-f1HaploNum
    len_eff=calcMapVal(seqs,input1,input2) 
    hVector=calcOrderedFrequencies(haploNum,haploSize,seqs)
    ordSeqs=orderPositions(hVector,seqs,haploSize)
    kStepList=calcDistances(haploNum,haploSize,ordSeqs)
    finalSeqs=makeNewEdgeList(kStepList,seqs,haploSize)
    D_all=calcDistanceMatrix(finalSeqs)
    evolTime=simulEvol(D_all,f1HaploNum,f2HaploNum,len_eff,failtime)
    return evolTime

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
        except IndexError:
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
def inferTransNetFromEvolTimes(DSamp,failtime):
    # print(DSamp)
    # raw_input("press enter to continue")
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

def specialquentin2files(inputs):
    firstToSecond=getEvolTime(inputs[0],inputs[1])
    secondToFirst=getEvolTime(inputs[1],inputs[0])
    if firstToSecond>secondToFirst:
        return "I think that %s is the source, and %s is the target " %inputs[1],inputs[0]
    elif firstToSecond<secondToFirst:
        return "I think that %s is the source, and %s is the target " %inputs[0],inputs[1]
    else:
        return "I cannot tell who the source is and who the target is"

#@profile    
def main(inputs,output): #overall wrapper - wraps the wrappers, checks some things
    failtime=3000
    startTime = time.time() 
    numFiles=len(inputs)
    if numFiles==2:
        specialquentin2files(inputs)
    fileList=[]
    DSamp=np.zeros([numFiles,numFiles],int)
    for i1 in range(numFiles):
        f1=inputs[i1]
        fileList.append(f1)
        for i2 in range(numFiles):
            f2=inputs[i2]
            if i1!=i2:
                evolTime=getEvolTime(f1,f2,failtime)
                DSamp[i1,i2]=evolTime
    isCluster=checkDSamp(DSamp,failtime)
    if isCluster==False:
        print("None of your samples appear to be related! Exiting")
        print(DSamp)
        endTime = time.time()
        workTime =  endTime - startTime
        statement_time='Analysis took %.3f seconds' % workTime
        print(statement_time)
        return 0
    TransNet=inferTransNetFromEvolTimes(DSamp,failtime)
    graphme(TransNet,fileList,output)
    endTime = time.time()
    workTime =  endTime - startTime
    statement_time='Analysis took %.3f seconds. Thank you for using QUENTIN!' % workTime
    print(statement_time)
    showpic='eog '+output+'.png'
    os.system(showpic)
    
if __name__ == '__main__':
    print("Initializing...") 
    import argparse
    parser = argparse.ArgumentParser(description='QUENTIN: predict directionality of HCV infection')
    parser.add_argument('-o','--output', type=str,required=False, default="network", help="The name of the output .png file")
    parser.add_argument('-f','--files', nargs="+",required=True, help="List of HCV cluster to be analyzed, order of files does not matter")
    #possible arguments to add: failtime, (failtime-20), delta, nIter, smaller file names?, sequencing platform
    args = parser.parse_args()
    output=args.output
    files=args.files
    main(files,output)
