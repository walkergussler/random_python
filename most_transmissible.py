from Bio import SeqIO, SeqRecord, Seq
import sys, os
import numpy as np

def calc_distance_matrix(finalSeqs): #calculate distance matrix from the 1-step list
    from ghost.util.distance import hamming
    l=len(finalSeqs)
    arr=np.zeros([l,l])
    hdist=hamming(finalSeqs,finalSeqs,ignore_gaps=False)
    for id in range(len(hdist)):
        item=hdist[id]
        arr[:,id]=item[:,0]
    return arr

def main(input,output,target): #get sequences from a file
    records={}
    seqs=[]
    with open(input,'r') as f:
        for record in SeqIO.parse(f,'fasta'):
            records[record.seq]=record
            seqs.append(record.seq)
    if len(records)!=len(seqs):
        sys.exit('records different length from seqs! report')
    nseq = len(seqs)
    from collections import defaultdict
    transmissibility=defaultdict(list)
    DM=calc_distance_matrix(seqs)
    for id in range(len(DM)):
        seq=seqs[id]
        h1=DM[id]
        centralTemp = np.divide(sum(h1), nseq-1,   dtype = float)
        transmissibility[centralTemp].append(seq)
    outD=[]
    for transVal,subseqs in sorted(transmissibility.items(),reverse=True):
        for seq in subseqs:
            outD.append(SeqRecord.SeqRecord(seq,id='>seq_1'))
        print(len(outD),target)
        if len(outD)>target:
            with open(output,'w') as f:
                SeqIO.write(outD,f,'fasta')
            return(input,len(seqs),len(outD))
    
    # return outD
    
# def main(inputs):
    # if len(inputs)==2:
        # print('input file: '+inputs[0])
        # print('output file: '+inputs[1])
        # parse_input(inputs[0],inputs[1])
    # else:
        # sys.exit('only input and output please!')
        
# if __name__=='__main__':
    # import argparse # possible arguments to add: delta, nIter
    # parser = argparse.ArgumentParser(description='QUENTIN: predict directionality of HCV infection')
    # parser.add_argument('files',
        # nargs='+',
        # help="input file followed by output file: program will not run if argument is incorrect length")

    # args = parser.parse_args()
    # INPUTS=args.files
    # main(INPUTS)

s={'aa':39,'ac':22,'ba':24,'bb':26,'bc':37,'ai':11,'aj':16,'bj':16,'aq':9,'aw':11}

for file in s:
    print(file+'_all.fas','normalized/'+file+'._some.fas',s[file])
    main(file+'_all.fas','normalized/'+file+'._some.fas',s[file])