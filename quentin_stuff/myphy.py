import numpy as np


def cluster(Z):
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
            T=clusternum(Z,T,int(i-m+1),clsnum)
            clsnum+=1
        i=Z[k-1,1]
        if i<=m-1:
            T[i]=clsnum
            clsnum+=1
        elif i<(2*m-maxclust+1):
            T=clusternum(Z,T,int(i-m+1),clsnum)
            clsnum+=1
    return T

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

def flattenHoriz(inarr): #mimics behavior of out=in(:) as out=flatten(in)
    rows=len(inarr)
    try:
        cols=len(inarr[0])
    except TypeError:
        out=inarr
        return out
    out=np.zeros(rows*cols)
    # print(inarr)
    it=-1
    for row in range(rows):
        for col in range(cols):
            it+=1
            # print(col,row)
            out[it]=inarr[row][col]
    return out
    
    
def clusternum(X,T,k,c):
    X=X.astype(int)
    T=T.astype(int)
    m=len(X)+1
    while type(k)==int or type(k)==np.float64 or len(k)>0:
        if isinstance(k,np.ndarray):
            k=k.astype(int)
        children=X[k-1,0:2]
        children=flatten(children)
        t=children<m
        try:
            T[children[t]]=c
        except IndexError:
            ele=int(children[t][0])
            T[ele]=c
        k=children[~t]-m
    return T

def prettyorder(B,D,names,numBranches,numLeaves,numLabels):
    # print("---prettyorder called---")
    # print(B)
    # print(D)
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
        # print("----------")
        # print(ind)
        # print(v)
        # print(Ls)
        # print(Li)
        # print(X[v[0]])
        # print(X[v[1]])
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
    
def myphy(linkage_out,nSamp):
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
        [B,names]=prettyorder(B,D,names,numBranches,numLeaves,numLabels)
    return(B,names)


import math
linkage=np.array([[0,2,17,2],[1,4,36,3],[3,5,70,4]])
lists,names=myphy(linkage,4)
print("---prettyorder finished---")
print(lists)
print(names)








# def calcLabelsPhyloTree(tree,DM,nSamp): #update values for tree *might* be off in some cases
    # print("called")
    # print(tree,DM,nSamp)
    # nVert=len(tree)
    # s=np.zeros(nVert-1)
    # t=np.zeros(nVert-1)
    # e=0
    # G=nx.Graph()
    # for i in range(nVert): 
        # G.add_node(i)
        # if tree[i,0]>0:
            # s[e]=i
            # t[e]=tree[i,0]# +-1's were here
            # e=e+1
        # if tree[i,1]>0:
            # s[e]=i
            # t[e]=tree[i,1]# +-1's were here
            # e=e+1
    # for i in range(len(s)):
        # G.add_edge(t[i],s[i])
    # print(G.edges())
    # print(G.nodes())
    # backOrder=[]
    # for a in nx.dfs_preorder_nodes(G): #does this do the same depth first search as matlab's dfsearch?
        # backOrder.append(a)
    # order=[]
    # for ele in reversed(backOrder):
        # order.append(ele)
    # internal=[]
    # for ele in range(nVert):
        # zero=tree[ele,0]>0
        # one=tree[ele,1]>0
        # internal.append(any([zero,one]))
    # w=1
    # for i in range(nVert):
        # print(i)
        # v=int(order[i])
        # if internal[v]==True:
            # child1=int(tree[v,0])
            # child2=int(tree[v,1])
            # l1=int(tree[child1,2])
            # l2=int(tree[child2,2])
            # if DM[l1,l2]<DM[l2,l1]:
                # tree[v,2]=l1
            # else:
                # tree[v,2]=l2
            # if DM[l1,l2]==1200 and DM[l2,l1]==1200:
                # w=0
    # return(tree,w)
    
    
# import networkx as nx
# a=np.array([[0,0,0],[0,0,1],[0,0,2],[0,0,3],[2,4,0],[0,5,0],[1,6,0]])
# dm=np.array([[0,84,222,77],[15,0,200,54],[10,25,0,34],[36,269,401,0]])
# calcLabelsPhyloTree(a,dm,4)
# print(a)
# print(dm)

# matlabtree=[0,0,1
# 0,0,2
# 0,0,3
# 0,0,4
# 1,3,0
# 2,5,0
# 4,6,0]


# DM =[0,87,228,80
# 15,0,202,54
# 10,18,0,34
# 36,231,404,0]
