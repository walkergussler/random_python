import numpy as np
import networkx as nx
from scipy.sparse.csgraph import connected_components,csgraph_from_dense,dijkstra
# from scipy.spatial.distance import pdist
# import fastcluster
# from scipy.optimize import lsq_linear

def graphshortestpath(dijkstra_output,source,target): #walks through output of dijkstra algorithm to find shortest path in unweighted, undirected graph - i think it works
    # print("graphshortestpath called")
    # print(dijkstra_output,source,target)
    numNodes=len(dijkstra_output)
    if source<0 or source>=numNodes or target<0 or target>=numNodes:
        raise ValueError("Source and target must be integers within [0,numNodes-1]")
    theList=[source]
    dist=int(dijkstra_output[source,target])
    for nodesVisited in range(dist):
        # print("nodes visited: %i" %nodesVisited)
        # print(theList)
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
            # print("customsource=%i"%source)
            if dijkstra_output[source,target]==1:
                theList.append(source)
                theList.append(target)
                return theList
            else:
                theList.append(source)

def objTransNetPhyloFitW(AMtree,tree,DM): #this program is stupid and mostly returns 0 as far as i can tell
    import lsqlin
    # print("-------------")
    # print(AMtree)
    # print(tree)
    # print(DM)
    # print("-------------")
    # modTree=matlabtree(tree)
    rec=0
    n=len(DM)
    am=len(AMtree)
    E=[]
    for row in range(am):
        for col in range(am):
            if row>=col and AMtree[row,col]!=0:
                E.append([col,row])
    nedges=len(E)
    nPairs=0
    for row in range(n):
        for col in range(n):
            if DM[row,col]!=0:
                nPairs+=1
    C=np.zeros([nPairs,nedges])
    d=np.zeros([nPairs]) 
    pathsPairs=[]
    p=-1
    for i in range(n):
        for j in range(n):
            if DM[i,j]>0:
                d[p]=1
                p+=1
                sparse_amtree=csgraph_from_dense(AMtree)
                distances=dijkstra(sparse_amtree)
                # print(distances)
                # print(distances)
                # print(i,j)
                path=graphshortestpath(distances,i,j)
                # print(path)
                # print("-------------")
                # print(path)
                pathsPairs.append(path)
                for v in range(len(path)-1):
                    if path[v]<path[v+1]:
                        # print("if",path[v],path[v+1])
                        edge=[path[v],path[v+1]]
                    else:
                        # print("else",path[v],path[v+1])
                        edge=[path[v+1],path[v]]
                    for ind in range(len(E)):
                        ele=E[ind]
                        if edge==ele:
                            C[p,ind]=1/(float(DM[i,j]))
    # print("p=%i"%p)
    deg=np.sum(AMtree,axis=0)
    # print(deg)
    root_tmp=np.where(deg==2)
    leaves_tmp=np.where(deg==1)
    # print(root_tmp)
    root=int(root_tmp[0])
    leaves=leaves_tmp[0]
    nleaves=len(leaves)
    nbranches=nleaves-1
    treeLen=nleaves+nbranches
    charVectPaths=np.zeros([nleaves,nedges])
    for leaf in leaves:
        path=graphshortestpath(distances,root,leaf)
        # print(path)
        # print("-------------")
        for v in range(len(path)-1):
            if path[v]<path[v+1]:
                edge=[path[v],path[v+1]]
            else:
                edge=[path[v+1],path[v]]
            for ind in range(len(E)):
                ele=E[ind]
                if edge==ele:
                    charVectPaths[leaf,ind]=1
    # print(charVectPaths)
    # for i in range(nleaves): #probably worthless loop @line55
        # for e in range(nedges):
            # edge=E[e]
            # if tree[edge[0],2] != tree[edge[1],2]: 
    Aeq=np.zeros([nbranches,nedges])
    beq=np.zeros([nbranches,1])
    # print(C)
    # print(d)
    for i in range(nbranches):
        Aeq[i]=charVectPaths[i]-charVectPaths[i+1]
    # print(Aeq)
    # lsq=lsq_linear(C,d,bounds=(0,10**9))
    lsq=lsqlin.lsqlin(C,d,0,None,None,Aeq,beq,0)
    x=lsqlin.cvxopt_to_numpy_matrix(lsq['x'])
    # print(type(x))
    # print(x)
    # x=lsq.x
    # print(x)
    # tree_weight=np.zeros([treeLen,5])
    # for e in range(nedges):
        # u=E[e][0]+1
        # v=E[e][1]+1
        # if modTree[u,0]==v or modTree[u,1]==v:
            # par=u
            # child=v
        # else:
            # par=v
            # child=u
        # cheatpar=par-1
        # if tree_weight[cheatpar][0]==0:
            # tree_weight[cheatpar][0]=child
            # tree_weight[cheatpar][2]=x[e]
        # else:
            # tree_weight[cheatpar][1]=child
            # tree_weight[cheatpar][3]=x[e]
    # for i in range(treeLen):
        # tree_weight[i][4]=modTree[i][2]
    DSamp_tree=np.zeros([n,n])
    p=-1
    pp=-1
    vec1=np.zeros(n*n)
    vec2=np.zeros(n*n)
    for i in range(n):
        for j in range(n):
            pp+=1
            vec1[pp]=DM[j][i] #am i building this incorrectly? supposed to be analagous to vec1=DM(:);
            if DM[i][j]>0:
                p+=1
                pathEdges=np.where(C[p]>0)
                # print(x[pathEdges])
                DSamp_tree[i][j]=np.sum(x[pathEdges])
    # print("tree here")
    # print(DSamp_tree)
    ppp=-1
    for i in range(n):
        for j in range(n):
            ppp+=1
            vec2[ppp]=DSamp_tree[j][i]
    ind=np.where(vec1>0)
    vec1=vec1[ind]
    vec2=vec2[ind]
    # try:
    # print(vec1)
    # print(vec2)
    val=np.corrcoef(vec1,vec2)
#    print(val)
    if val.ndim>1:
        for row in range(len(val)):
            for col in range(len(val)):
                if abs(val[row,col]-1)>.00005:
                    # print(val[row,col], "made it")
                    if val[row,col]>rec:
                        rec=val[row,col]
    return rec
    # except ValueError:
        # return 0


AMtree=np.array([[0,0,0,1,0],[0,0,0,1,0],[0,0,0,0,1],[1,1,0,0,1],[0,0,1,1,0]])

tree=np.array([[0,0,0],[0,0,1],[0,0,2],[0,1,1],[2,3,3]])
DM=np.array([[0,34,0],[0,0,0],[67,42,0]])
# DM=np.array([[0,0,0],[55,0,372],[595,0,0]])
# tree=np.array([[0,0,0],[0,0,1],[0,0,2],[2,0,0],[1,2,2]])    

# IdealTree=np.array([[0,0,1],[0,0,2],[0,0,3],[1,2,2],[3,4,2]])


newtree=objTransNetPhyloFitW(AMtree,tree,DM)
# newtree=calcLabelsPhyloTree(tree,DM,3)

print(newtree)
