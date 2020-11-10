def calcDistanceMatrix(finalSeqs): #SLOW 
    l=len(finalSeqs)
    arr=np.zeros([l,l])
    for id1 in range(l):
        seq1=finalSeqs[id1]
        for id2 in range(l):
            seq2=finalSeqs[id2]
            dist=sum(0 if a==b else 1 for a,b in itertools.izip(seq1,seq2))
            arr[id1][id2]=dist
            arr[id2][id1]=dist
    return arr

