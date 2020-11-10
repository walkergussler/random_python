from __future__ import division
from itertools import combinations
from Bio import SeqIO
from ghost.util.distance import hamming
import sys
import copy
import community
import networkx as nx
import numpy as np
# from tempfile import NamedTemporaryFile
# import subprocess
# import os
# import re
# import time
# import collections
# import argparse
# import math

#program breaking conditions:
    #assumes same sequence will not appear twice in one file (not sure which error this will cause)
    #requires frequency for each read to follow last underscore. If anything following last underscore is not frequency,the program will error out

class UnionFind(object): # this could be replaced with an import statement
    # # # An implementation of union find data structure.
    # # # It uses weighted quick union by rank with path compression.

    def __init__(self,node_ids):
        """\
        Initialize an empty union find object.
        :param node_ids: The identifiers of the nodes.
        """
        self._sets={
                node_id : {
                    'rank' : 0,
                    'parent' : node_id
                }  for idx,node_id in enumerate(node_ids)
        }
        self.component_count=len(self._sets)

    def find(self,x):
        try:
            p_idx=self._sets[x]['parent']
            if p_idx != x:
                self._sets[x]['parent']=self.find(self._sets[p_idx]['parent'])
            return self._sets[x]['parent']
        except KeyError:
            raise KeyError('ID {0}is not a member of the union'.format(x))

    def join(self,p,q):
        # # # Combine sets containing p and q into a single set.
        p_id=self.find(p)
        q_id=self.find(q)
        pRoot=self._sets[p_id]
        qRoot=self._sets[q_id]
        if p_id != q_id:
            self.component_count -= 1
        if pRoot['rank']<qRoot['rank']:
            pRoot['parent']=q_id
        elif pRoot['rank'] > qRoot['rank']:
            qRoot['parent']=p_id
        else:
            qRoot['parent']=p_id
            pRoot['rank'] += 1

    def connected(self):
        it=iter(self._sets)
        f=self.find(next(it))
        for i in it:
            if self.find(i) != f:
                return False
        return True
        
def getseqs(file):
    seqs={}
    seqnames={}
    count=0
    with open(file,"rU") as input_handle:
        for record in SeqIO.parse(input_handle,"fasta"):
            if record.seq not in seqs:
                count+=1
                seqs[record.seq]=count
                seqnames[count]=record
    return seqs,seqnames
        
def make_allDistTuple_fast(seqs):
        dist_array=calcDistanceMatrix(seqs)
        ourStructure=[]
        if len(seqs)!=max(seqs.values()):
            print("len seqs != max seqs.values")
        else:
            print(len(seqs))
        for (seq1,seq2) in combinations(seqs,2):
            node1=int(seqs[seq1])
            node2=int(seqs[seq2])
            dist=dist_array[node1-1,node2-1]
            ourStructure.append((dist,(node1,node2)))
        retstruct=iter(sorted(ourStructure,key=lambda  t:t[0]))
        g=nx.Graph()
        for seq in seqs:
                g.add_node(seqs[seq])
        return retstruct,g
        
def calcDistanceMatrix(finalSeqs):  #calculate  distance  matrix  from  list  of  sequences
    l=len(finalSeqs)
    arr=np.zeros([l,l])
    hdist=hamming(finalSeqs,finalSeqs,ignore_gaps=False)
    
    for id in range(len(hdist)):
        item=hdist[id]
        arr[:,id]=item[:,0]
    return arr

def kstep(distances,g,threshold=float('inf')):
    """\
    Build a k-step network.
    :param d_iter: An iterable which yields distances in the form of a tuple of tuples (distance,(node name 1,node name 2))
    :param g: A NetworkX graph initialized with nodes using ids from the iterable.
    :param connect: Throw an exception if the graph is not connected by the k-step network.
    """
    d_iter=iter(distances)
    uf=UnionFind(g.nodes())
    current,(i,j)=next(d_iter)
    d_next,(n_i,n_j)=current,(i,j)
    
    try:
        while current<threshold and not uf.connected():
            next_uf=copy.deepcopy(uf)
            
            while d_next == current:
                if uf.find(n_i) != uf.find(n_j):
                    g.add_edge(n_i,n_j,len=d_next)
                    next_uf.join(n_i,n_j)
                d_next,(n_i,n_j)=next(d_iter)

            uf=next_uf
            current=d_next
            i=n_i
            j=n_j

    except StopIteration:
        if current<threshold and not uf.connected():
            raise ValueError('kstep network did not connect network.')
    return g
    
def dend(file):
    seqs,seqnames=getseqs(file)
    allDistTuple,g=make_allDistTuple_fast(seqs)
    g=kstep(allDistTuple,g)
    dend=community.generate_dendrogram(g,weight='len')
    # print('-')
    for item in dend:
        print(item)
    # print('-')
    return dend, seqnames
    
def best(file):
    seqs,seqnames=getseqs(file)
    allDistTuple,g=make_allDistTuple_fast(seqs)
    g=kstep(allDistTuple,g)
    best=community.best_partition(g,weight='len',randomize=True)
    for item in best:
        print(item)
    return best, seqnames
    
def p():
    print("You entered something weird! Welcome to this shoddy help menu!")
    print("this program has 2 modes currently, build a dendrogram and find best partition")
    print("best will compute the best parition of the graph nodes which maximizes modularity")
    print("dend will find communities in the graph and return the associated dendrogram [list of dictionaries]")
    print("you can run this script from the command line with the following options")
    print("./path/to/louvain.py d <input_fasta.fas>")
    print("for dend mode, or to do best mode,")
    print("./path/to/louvain.py b <input_fasta.fas>")
    print("Alternatively, you can import me as a module!")
    print("Just copy this script to whatever directory your thing is in, and say")
    print("import louvain")
    print("Then you can call functions as follows")
    print("bestPartition=louvain.best(<input_fasta.fas>")
    print("dendrogram=louvain.dend(<input_fasta.fas>")
    print("this module has other functions as well! I will print their short usage here, and can write similar functions to dend or best for these easily, just ask if you would like one!")
    print("-")
    print("Format: description of method goes here")
    print("and the example goes here")
    print("this line has return structure info")
    print("-")
    print('Induced graph: given a partition[PART] and graph[g], produce the graph where the nodes are communities')
    print("community.induced_graph(PART,g,weight='len')")
    print("nx graph where nodes are the parts")
    print("-")
    print('Modularity: compute the modularity of a partition of a graph')
    print("community.modularity(PART,g,weight='len')")
    print("floating point for modularity")
    print("-")
    print('partition_at_level: return the partition of the nodes at the given level')
    print("community.partition_at_level(dendrogram,level)")
    print("partition: dictionary where keys are the nodes and the values are the set it belongs to")
    print("-")
    
    
if __name__=="__main__":
    # try:d
    modeSelector=sys.argv[1]
    file=sys.argv[2]
    if modeSelector=='d':
        dend,seqnames=dend(file)
        #dend[0:2]
        betterdend={}
        for val in range(max(dend[0].values())+1):
            q=[]
            for item in dend[0]:
                if dend[0][item]==val:
                    q.append(item)
            betterdend[val]=q
        # for item in betterdend:
            # print(item,betterdend[item])
        # raw_input("ok?")
        besterdend={}
        for val in range(max(dend[1].values())+1):
            q=[]
            for item in dend[1]:
                if dend[1][item]==val:
                    for subitem in betterdend[item]:
                        q.append(subitem)
            besterdend[val]=q
        for item in besterdend:
            besterdend[item]=list(set(besterdend[item]))
        L=[]
        for item in besterdend:
            print(item,besterdend[item])
            for s in besterdend[item]:
                L.append(s)
            # print(item,len(besterdend[item]))
        print(len(L))
        print(len(besterdend))
        print(max(dend[1].values()))
        for item in besterdend:
            seqs=[]
            for sub in besterdend[item]:
                seqs.append(seqnames[sub])
            with open(str(item)+'.fas','w') as f:
                SeqIO.write(seqs,f,'fasta')
    elif modeSelector=='b':
        best(file)
    else:
        p()
    # except:
        # p()
    
# for level in range(len(dend) - 1) :
    # print("partition at level", level, "is", partition_at_level(dend, level))
    
    

#modularity
