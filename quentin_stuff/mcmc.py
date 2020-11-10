import numpy as np
from scipy.sparse.csgraph import csgraph_from_dense, dijkstra
import fastcluster, math, random

#fix tmp

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
    
def clusternum(X,T,k,c):
    # print(X)
    # print(T)
    # print(k)
    # print(c)
    X=X.astype(int)
    T=T.astype(int)
    m=len(X)+1
    # print("m %i"%m)
    while type(k)==int or type(k)==np.float64 or len(k)>0:
        if isinstance(k,np.ndarray):
            k=k.astype(int)
        children=X[k-1,0:2]
        children=flatten(children)
        t=children<m
        # print("-")
        # print(k)
        # print(children)
        # print(t)
        # print(children[t])
        # print(type(children[t]))
        # print("-")
        if np.size(children[t])!=0:
            try:
                T[children[t]]=c
            except IndexError:
                ele=int(children[t][0])
                T[ele]=c
        k=children[~t]-m
        # print(k)
        # print(type(k))
        if len(k)<2 or type(k)==int or type(k)==float:
            if k==0 or k==[0]:
                return T
    return T
    
def estimateHubs(TransNet): 
    AM=(TransNet+np.transpose(TransNet))>0
    deg=sum(AM)
    deg=sorted(deg,reverse=True)
    transposed=[]
    for item in deg:
        transposed.append([item])
    dists=pdist(transposed)
    linkage_out=fastcluster.linkage(dists,method='single')
    # print(linkage_out)
    clustered=cluster(linkage_out)
    c=clustered[0]
    k=len(np.where(clustered==c))
    return k

def sCore(n,k):
    n=float(n)
    k=float(k)
    d1=(n-k)/k+k-1
    d2=(n-k)/k+1
    return (n-k)/k*d1+(k-1)*(n-k)/k*d2+(k-1)*d1*d2

def objTransNetQuadDeg2(DM,nhubs):
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

def modifyTree1(tree,steps,nSamp): 
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

def findDSampDir(DSamp): 
    AMSamp_dir=(DSamp<=np.transpose(DSamp))
    bools=DSamp<=np.transpose(DSamp)
    AMSamp_dir=bools.astype(int)
    return DSamp*AMSamp_dir

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

def calcLabelsPhyloTree(tree,DM,nSamp): #update values for tree *might* be off in some cases
    # print("called")
    # print(tree,DM,nSamp)
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
    for a in nx.dfs_preorder_nodes(G): #does this do the same depth first search as matlab's dfsearch?
        backOrder.append(a)
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
        # print(i)
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
            if DM[l1,l2]==1200 and DM[l2,l1]==1200:
                w=0
    return(tree,w)

def prettyorder(B,D,names,numBranches,numLeaves,numLabels):
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

def treeToTransNet(tree,DSamp,nSamp): #seems good
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
    
def treeAM(tree): 
    numLabels=len(tree)
    # tmp_tree=numpy2matlabTree(tree,numLabels) #old function i thought i needed
    n=len(tree)
    AMtree=np.zeros([n,n])
    for i in range(n):
        child1=int(tree[i,0])
        child2=int(tree[i,1])
        if child1+child2>0:
            # print(child1,child2,i)
            AMtree[i,child1]=1
            AMtree[i,child2]=1
            AMtree[child1,i]=1
            AMtree[child2,i]=1
    return AMtree
    
def graphshortestpath(dijkstra_output,source,target): #walks through output of dijkstra algorithm to find shortest path in unweighted, undirected graph - i think it works
    # print("graphshortestpath called")
    # print(dijkstra_output,source,target)
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

def objTransNetPhyloFitW(AMtree,tree,DM): 
    # print(AMtree)
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
                # print(AMtree)
                # print(sparse_amtree)
                # print(distances)
                path=graphshortestpath(distances,i,j)
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
    # print(tree)
    # print(AMtree)
    # print(root_tmp)
    root=int(root_tmp[0])
    leaves=leaves_tmp[0]
    nleaves=len(leaves)
    nbranches=nleaves-1
    treeLen=nleaves+nbranches
    charVectPaths=np.zeros([nleaves,nedges])
    for leaf in leaves:
        path=graphshortestpath(distances,root,leaf)
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
                
def findTransNetMCMC(DSamp_mcmc): 
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
    # print(l)
    [branchLists,nodeNames]=myphy(l,nSamp)
    # print("branchLists")
    # print(branchLists)
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
    # print(tree)
    tree,w=calcLabelsPhyloTree(tree,DSamp_mcmc,nSamp)
    TransNet=treeToTransNet(tree,DSamp_mcmc,nSamp)
    k=estimateHubs(TransNet)
    AMtree=treeAM(tree)
    matrToFit=DSamp_dir
    obj1=objTransNetPhyloFitW(AMtree,tree,matrToFit)
    obj2LinkNoPow=objTransNetQuadDeg2(TransNet,k)
    obj2=obj2LinkNoPow
    record=w*obj1*obj2
    TransNetTrees=[]
    TransNetTrees.append(tree)
    delta=.0001
    iEdgeMotif=1
    nEdgeMotifCurr=2
    nIter=250
    i=-1
    while i<=nIter*2:
        i+=1
        if iEdgeMotif==nIter+1:
            iEdgeMotif=1
            nEdgeMotifCurr-=1
        # print("trees")

        newtree=modifyTree1(tree,nEdgeMotifCurr,nSamp)
        if np.all(newtree==0):
            continue
        AM_newtree=treeAM(newtree)
        # for a in range(len(newtree)):
            # print(a,newtree[a])
        # print(newtree)
        # print(AM_newtree)
        # print(sum(AM_newtree))
        sparse_amtree=csgraph_from_dense(AM_newtree)
        distances=dijkstra(sparse_amtree)
        # print(sparse_amtree)
        # print(distances)        
        newtree,w = calcLabelsPhyloTree(newtree,DSamp_mcmc,nSamp)
        TransNetNew=treeToTransNet(newtree,DSamp_mcmc,nSamp)
        matrToFit=DSamp_dir
        obj1=objTransNetPhyloFitW(AM_newtree,newtree,matrToFit)
        obj2=objTransNetQuadDeg2(TransNetNew,k)
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
        tree=TransNetTrees[c]
        AM_tree=treeToTransNet(tree,DSamp_mcmc,nSamp)
        isExist=False
        for j in range(len(TransNets)):
            if np.all(np.equal(TransNets[j],AM_tree)):
                isExist=True
                break
        if not isExist:
            TransNets.append(AM_tree)
    return TransNets

    
import numpy as np  
import networkx as nx
from scipy.spatial.distance import pdist

# a=np.array([[0,908],[95,0]])
a=np.array([[0,910],[92,0]])
# a=np.array([[0,0,0],[55,0,368],[591,0,0]])

BB=np.array([[0,40,22,39,1200,499,1200],[57,0,35,42,1200,369,1200],[35,35,0,38,1200,354,1200],[33,27,21,0,1200,161,1200],[317,84,121,292,0,73,1200],[97,35,41,44,931,0,1200],[161,59,72,81,1200,590,0]])
AQ=np.array([[0,142,312,31,35,25,8,14,8],[48,0,401,96,269,62,38,84,49],[99,310,0,511,1200,213,42,100,100],[31,194,277,0,139,38,31,31,30],[14,468,654,9,0,33,14,15,14],[36,1200,599,36,125,0,14,41,41],[25,1043,518,31,70,23,0,35,36],[9,317,220,8,35,25,9,0,9],[2,146,276,2,35,25,9,9,0]])
# BA=np.array([[0,1200,1200,1200,1200,257],[1200,0,230,1200,620,1200],[1200,160,0,1200,268,1200],[671,20,36,0,186,982],[1200,315,91,1200,0,1200],[1200,1200,1200,1200,1200,0]])

print("------------------------BEGIN------------------------")

# BJ=np.array([[0.,84.,222.,77.],[15.,0.,200.,54.],[10.,25.,0.,34.],[36.,269.,401.,0.]])
# b=findTransNetMCMC(BJ)

# AQ=np.array([[0,142,312,31,35,25,8,14,8],[48,0,401,96,269,62,38,84,49],[99,310,0,511,1200,213,42,100,100],[31,194,277,0,139,38,31,31,30],[14,468,654,9,0,33,14,15,14],[36,1200,599,36,125,0,14,41,41],[25,1043,518,31,70,23,0,35,36],[9,317,220,8,35,25,9,0,9],[2,146,276,2,35,25,9,9,0]])
# b=findTransNetMCMC(AQ)

AW=np.array([[0,21,26,2,15,26,9,15,46,28,542,15,428,26,9,256,21,92,26],[32,0,30,14,20,26,14,20,116,37,484,21,415,32,15,287,26,136,30],[37,32,0,19,26,32,20,26,93,37,364,27,351,31,20,102,32,38,26],[26,19,25,0,14,25,8,14,44,30,423,14,280,25,8,204,19,64,25],[25,20,25,2,0,25,9,14,43,31,504,15,354,26,9,242,20,85,25],[38,31,36,19,25,0,20,26,283,45,667,26,646,38,20,304,31,144,31],[26,20,25,2,14,23,0,14,46,31,493,14,367,26,9,247,9,87,26],[26,20,26,2,14,26,9,0,40,31,447,9,373,26,9,249,20,88,26],[33,33,22,18,19,34,19,25,0,20,282,26,330,27,19,295,31,40,28],[40,36,52,23,29,38,24,30,57,0,1094,31,347,42,24,403,36,160,42],[104,43,43,41,53,99,42,50,560,194,0,55,1200,44,43,169,53,69,58],[27,21,27,2,15,26,9,15,52,32,594,0,427,27,9,282,21,99,27],[39,32,40,19,26,37,20,26,150,28,620,27,0,40,20,35,32,9,38],[27,26,25,14,21,26,15,21,94,37,357,21,500,0,15,263,26,44,31],[27,20,26,2,14,26,9,14,48,31,549,15,388,27,0,258,20,92,26],[52,26,40,26,33,34,27,33,273,55,650,35,36,57,27,0,39,15,33],[27,20,26,2,14,26,9,14,46,31,533,15,382,26,9,255,0,91,26],[215,210,83,60,139,214,82,141,690,294,562,145,367,140,84,149,208,0,127],[32,26,26,14,19,26,14,20,91,37,217,21,514,32,14,159,26,45,0]])
b=findTransNetMCMC(a)

print(b)


