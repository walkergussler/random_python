def phaceliaAPI(INPUTS):
    #determine source via most chronic phacelia prediction
    #not used, but likely will be in the future
    seqsList=[]
    for input in INPUTS:
        with open(input) as f:
            for record in SeqIO.parse(f,'fasta'):
                splitid=record.id.split('_')                
                record.annotations['freq']=int(splitid[-1])
                record.annotations['sample']=input
                record.annotations['genotype']=splitid[1]
                seqsList.append(record)
    
    seqs=iter(seqsList)
    pdic={}
    for item in recency_bin([seqs]):
        pdic[item[2]]=item[0]
        # print(item[2],item[0])
    winner=max(pdic.keys())
    return pdic
