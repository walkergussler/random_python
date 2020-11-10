def findIntermediateSequence(seq1,seq2):
    import random
    seqNew=""
    while len(seq1)<len(seq2):
        seq1=seq1+"_"
    while len(seq2)<len(seq1):
        seq2=seq2+"_"
    if len(seq1)!=len(seq2):
        print("+im dumb!")
    else:
        for pos in range(len(seq1)):
            if seq1[pos]==seq2[pos]:
                seqNew=seqNew+seq1[pos]ls
            else:
                seqNew=seqNew+random.sample([seq1[pos],seq2[pos]],1)[0]
    return seqNew
    
def findMedianSequence(seq1,seq2,seq3):
    import random
    all the same length?
    seqNew=""
    for pos in range(max(len(seq1),len(seq2),len(seq3))):
        if len(seq1)<=pos:
            seq1=seq1+"_"
        if len(seq2)<=pos:
            seq2=seq2+"_"
        if len(seq3)<=pos:
            seq3=seq3+"_"
        chars=seq1[pos]+seq2[pos]+seq3[pos]
        numblanks=chars.count("_")
        if numblanks==0:
            if seq1[pos]==seq2[pos]:
                seqNew=seqNew+seq1[pos]
            elif seq3[pos]==seq1[pos] or seq3[pos]==seq2[pos]:
                seqNew=seqNew+seq3[pos]
            else:
                seqNew=seqNew+random.sample([seq1[pos],seq2[pos],seq3[pos]],1)[0]
        elif numblanks==1:
            noblanks=chars.replace("_","")
            if noblanks.count(noblanks[0])==len(noblanks):
                seqNew=seqNew+noblanks[0]
            else:
                seqNew=seqNew+random.sample([noblanks[0],noblanks[1]],1)[0]
        elif numblanks==2:
            seqNew=seqNew+chars.replace("_","")
        elif numblanks==3:
            seqNew=seqNew+"_"
    return seqNew