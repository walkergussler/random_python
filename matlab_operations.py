def numpy2matlabTree(intree,numLabels): # i use this and it works
    outtree=np.zeros([numLabels,3])
    numNodes=(numLabels+1)/2
    for i in range(numLabels):
        if i<numNodes:
            outtree[i][2]+=1+intree[i][2]
        else:
            for j in range(3): 
                if j!=2:
                    outtree[i][j]+=intree[i][j]
                else:
                    outtree[i][j]+=1+intree[i][j]
    return outtree
    
    
def matlab2numpyTree(intree): # not useful (maybe), works
    numNodes=(len(intree)+1)/2
    for i in range(len(intree)):
        if i<numNodes:
            intree[i][2]-=1
        else:
            for j in range(3): #SHOULD THIS BE A 2?
                intree[i][j]-=1
    return intree  
    
#DSamp_flat code
DSamp_flat=np.zeros(nSamp*(nSamp-1)/2)
flatind=0
for u in range (nSamp):
    for v in range(u+1,nSamp):
        if DSamp_mcmc[u,v]<=DSamp_mcmc[v,u]:
            DSamp_flat[flatind]=DSamp_mcmc[u,v]
        else:
            DSamp_flat[flatind]=DSamp_mcmc[v,u]
