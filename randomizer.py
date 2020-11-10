def Randomizer(dic,seqnum,seqlen):
    tmp=[]
    seqs=dic.keys()
    freq=dic.values()
    for seq in seqs:
        tmp.append(list(seq))
    seqMat=np.array(tmp)
    tmp = np.zeros(np.shape(seqMat),dtype='S1')
    for i in range(seqlen):
        col = seqMat[:,i]
        tmp[:,i] = np.random.permutation(col)
    seqnew=[]
    for i in range(seqnum):
        seq='',join(seqMat[i])
        seqnew.append(seq)
    #find unique sequences among all shuffled sequences
    seqset=set(seqnew)
    #find frequences of unique sequences
    newFreq2=np.zeros(size(seqset,1),1)
    for k in range(len(seqset)):
        for l in range(len(seqnew)):
            if seqnew(l)==seqset(k):
            # if seqnew(l)==seqset(k):
                newFreq2[k]=newFreq2[k]+freq[l]
    return [newFreq2, seqset]
