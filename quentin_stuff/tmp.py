if __name__=='__main__':
    #setup
    # 
    printset=['AA','AC','AJ','AW','BA','BB','BC','BJ','AQ']
    sources=['AA45','AC124','AJ199','AW2','BA3','BB45','BC46','BJ28','AQ89']
    ###################################################################################
    # printset=['']
    files=[]
    
    # source='AA45_unique_1b_161.fas' # source='AA45_unique_1b_161.fas'    # source='AW2_unique_1b_49.fas'#################
    # for file in os.listdir(os.getcwd()):
        # if file.endswith('fas'): #################
            # files.append(file)
    ###################################################################################
    # source=sys.argv[1]
    # files=sys.argv[1:]
    ####################################################################################
    kmerCount={}
    kmerEntropy={}
    nucDivs={}
    entropies={}
    maxHammings={}
    files=['AA20_unique_1b_43.fas','AA8_unique_1b_88.fas','AA45_unique_1b_161.fas','AC122_unique_1a_48.fas','AC123_unique_1a_48.fas','AC124_unique_1a_47.fas','AC121_unique_1a_48.fas','AJ199_unique_1b_29.fas','AJ6_unique_1b_68.fas','AJ86_unique_1b_32.fas','AW28_unique_1b_48.fas','AW12_unique_1b_48.fas','AW21_unique_1b_48.fas','AW18_unique_1b_48.fas','AW10_unique_1b_48.fas','AW35_unique_1b_48.fas','AW34_unique_1b_48.fas','AW7_unique_1b_47.fas','AW22_unique_1b_47.fas','AW20_unique_1b_47.fas','AW16_unique_1b_52.fas','AW1_unique_1b_38.fas','AW9_unique_1b_48.fas','AW15_unique_1b_6.fas','AW27_unique_1b_47.fas','AW2_unique_1b_49.fas','AW4_unique_1b_47.fas','AW8_unique_1b_48.fas','AW13_unique_1b_48.fas','BA2_unique_1b_10.fas','BA3_unique_1b_52.fas','BA4_unique_1b_23.fas','BA8_unique_1b_33.fas','BA16_unique_1b_15.fas','BA1_unique_1b_9.fas','BB29_unique_1a_59.fas','BB42_unique_1a_84.fas','BB45_unique_1a_177.fas','BB41_unique_1a_30.fas','BB1_unique_1a_99.fas','BB44_unique_1a_12.fas','BB31_unique_1a_94.fas','BC46_unique_1a_69.fas','BC30_unique_1a_62.fas','BJ23_unique_1a_7.fas','BJ25_unique_1a_51.fas','BJ30_unique_1a_6.fas','BJ28_unique_1a_66.fas','AQ28_unique_1a_5.fas','AQ6_unique_1a_10.fas','AQ4_unique_1a_39.fas','AQ90_unique_1a_16.fas','AQ89_unique_1a_15.fas','AQ25_unique_1a_12.fas','AQ13_unique_1a_9.fas','AQ16_unique_1a_11.fas','AQ92_unique_1a_23.fas']
    legend=['filename','kmer entropy','nucleotide diversity','nucleotide entropy','kmer count','max hamming']
    print(", ".join(legend))
    for file in files:
        newkmers=[]
        seqs=getseqs(file)
        haploSize=len(six.next(six.iterkeys(seqs)))
        for seq in seqs:
            if len(seq)!=haploSize:
                sys.exit("go align")
        haploNum=len(seqs)
        kmers=getkmers(seqs,14)
        
        s=seqs.keys()
        array=calcDistanceMatrix(s,s)
        maxHammings[file]=np.amax(array)
        
        kmerCount[file]=len(getkmers(seqs,53))
        
        tot=float(sum(kmers.values()))
        for i in kmers.values():
            freq=i/tot
            newval=freq*math.log(freq)
            newkmers.append(newval)
        kmerEntropy[file]=sum(newkmers)*-1
        
        entropy=calcOrderedFrequencies(seqs,haploNum,haploSize)
        ent=sum(entropy)/len(entropy)
        entropies[file]=ent
        totalFreq=float(sum(seqs.values()))
        freqs=[x/totalFreq for x in seqs.values()]
        mat=calcDistanceMatrix(s,s)
        it=range(haploNum)
        nucDiv=0
        for a,b in combinations(it,2):
            nucDiv+=freqs[a]*freqs[b]*mat[a,b]*2/float(haploSize)
        nucDivs[file]=nucDiv
        plist=[file,kmerEntropy[file],nucDivs[file],entropies[file],kmerCount[file],maxHammings[file]]
        print(", ".join(map(str,plist)))
    print(kmerEntropy.values())
    print('kmer entropy vs nuc div')
    print(r2_score(kmerEntropy.values(),nucDivs.values()))
    print('kmer entropy vs nuc entropy')
    print(r2_score(kmerEntropy.values(),entropies.values()))
    print('kmer entropy vs kmer count')
    print(r2_score(kmerEntropy.values(),kmerCount.values()))
    print('kmer entropy vs max hamm')
    print(r2_score(kmerEntropy.values(),maxHammings.values()))
    print('nuc div vs nuc entoropy')
    print(r2_score(nucDivs.values(),entropies.values()))
    print('nuc div vs kmer count')
    print(r2_score(nucDivs.values(),kmerCount.values()))
    print('nuc div vs max hamming')
    print(r2_score(nucDivs.values(),maxHammings.values()))
    print('nuc entropy vs kmer count')
    print(r2_score(entropies.values(),kmerCount.values()))
    print('nuc entropy vs max hamming')
    print(r2_score(entropies.values(),maxHammings.values()))
    print('kmer count vs max hamming')
    print(r2_score(kmerCount.values(),maxHammings.values()))
    sys.exit()
    # AAset=[];ACset=[];AJset=[];AWset=[];BAset=[];BBset=[];BCset=[];BJset=[];AQset=[];
    
    # raw_input(len(kmerEntropy))
    # for letters in printset:
        # valDict={}
        # # for asdf in range(len(totalKmsers)):
            
        # for item in kmerEntropy:
            # if item.startswith(letters):
        # ######################################################################                
                # plist=[item,kmerEntropy[item],nucDivs[item],entropies[item],kmerCount[item],maxHammings[item]]
                # print(", ".join(map(str,plist)))
        # print('')
    ######################################################################
    sourceboyz=[]
    infecteds=[]
    for item1,item2 in combinations(kmerEntropy,2):
        print(item1,item2)
        v1=abs(kmerEntropy[item1]-kmerEntropy[item2])
        v2=abs(nucDivs[item1]-nucDivs[item2])
        v3=abs(entropies[item1]-entropies[item2])
        v4=abs(kmerCount[item1]-kmerCount[item2])
        v5=abs(maxHammings[item1]-maxHammings[item2])
        print([v1,v2,v3,v4,v5])
        if item1==source or item2==source:            
            sourceboyz.append([v1,v2,v3,v4,v5])
        else:
            infecteds.append([v1,v2,v3,v4,v5])
    print("source")
    print(sourceboyz)
    print("inf")
    print(infecteds)
    q=[]
    for a in range(5):
        roc_predict=[]
        roc_true=[]
        for item in sourceboyz:
            roc_predict.append(item[a])
            roc_true.append(1)
        for item in infecteds:
            roc_predict.append(item[a])
            roc_true.append(0)
        print(legend[a+1])
        auc=roc_auc_score(roc_true,roc_predict)
        print(auc)
        q.append(auc)
    print(q)
    print(q[0]*.2+q[1]*.4+q[2]*.3+q[3]*.05+q[4]*.05)