import itertools
import networkx as nx
from ghost.util.distance import hamming
import numpy as np

def calcDistanceMatrix_fast(finalSeqs): #BROKEN
    l=len(finalSeqs)
    arr=np.zeros([l,l])
    hdist=hamming(finalSeqs,finalSeqs,ignore_gaps=False)
    for id in range(len(hdist)):
        item=hdist[id]
        arr[:,id]=item[:,0]
    return arr

def calcDistanceMatrix_slow(finalSeqs): #SLOW AF
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
    
def tmp(s1,s2):
    seqs=[s1]
    if len(s1)==len(s2):
        dist=sum(0 if a==b else 1 for a,b in itertools.izip(s1,s2))
        if dist !=1:
            while dist>1:
                for nucl in range(len(s1)):
                    if s1[nucl]!=s2[nucl]:
                        takeseq1=s1[:nucl+1]+s2[nucl+1:]
                        takeseq2=s1[:nucl]+s2[nucl:]
                        if takeseq1!=s1:
                            newseq=takeseq1
                        else:
                            newseq=takeseq2
                print(takeseq1,takeseq2)
                dist-=1
                s1=newseq
                seqs.append(newseq)
    seqs.append(s2)
    for item in seqs:
        print(item)
    
s1='acg'
s2='tac'
    
    
# seqsin=['CAGAGACTTGTCGTTTTTATGTTACTACAACATATCATACCACGATGCGA','AAGTTAGGAGCCAATTCCTATGCTGATAGTGTTCTATGCTATTGAGACTG','TTGAATAGTTTATGAACAGTCGCTCTGAACCCCAGGGATTCTGGGGTCCT','TGAGTGATCCATGCGGCGGGGTTACCGCTATGATCAGCCCAGACGGGGCT','TTTAATACCCCCGCCTACTGTTGAGCAAATCCTGACTATGAAAACTCACC','CAAACGGACATGTGAGACGCAGGTCCATATCCCGACCGAACTATGTCGCC','CATCCAATCTGGCTCCCAACGTGTACAAGCGCTCTTACACTTTTGCATAC','AATGGGCGGGTTGTGCTTTCTATGATTGTCGCCCGTTCCTCTACATTTAA','TATCTGAAGTATCATTACGCTACAGATTGATTATCACTCATGGAGAACAG','TCTCCTCGACCGTCAGGCGTCGACTGCTTCATAGAAGCAGCGGGTAACTA']    
# print(calcDistanceMatrix_fast(seqsin))
# print(calcDistanceMatrix_slow(seqsin))