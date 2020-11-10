#!/usr/bin/env python

__author__ = "CDC/OID/NCHHSTP/DVH bioinformatics team"

"""k_step_one file. Program for calculating a k-step network form a single file.
1. It requires as input the location of the aligned fasta files. All sequences are assumed to be different.
2. It calculates the k-step network using an iterative multi-index hashing approach.
3. It does quentin stuff 

module load Python/2.7.11

With default parameters:
python kquentin.py <file1> <file2> ... <fileN>
"""

#np.transpose -> faster than for loop?

def parseInput(input1,input2): # parse input & calculate frequences
	from Bio import SeqIO
	from Bio.Seq import Seq
	import re
	eqVert=0
	seqs=[]
	input_handle = open(input1) 
	input_handle2 = open(input2)
	for record in SeqIO.parse(input_handle, "fasta"): # for FASTQ use "fastq", for fasta "fasta"
		if len(record.seq) > 0 and len(record.seq) < 50000:
			seqs.append(record.seq)
	input_handle.close()
	
	f1HaploNum=len(seqs)
	for record in SeqIO.parse(input_handle2, "fasta"): # for FASTQ use "fastq", for fasta "fasta"
		if len(record.seq) > 0 and len(record.seq) < 50000:
			if record.seq not in seqs:
				seqs.append(record.seq)
			else:
				print("WARNING! At least one of the sequences in f1 is identical to a sequence in f2; skipping")
				eqVert+=1
	
	input_handle.close()
	input_handle2.close()
	haploNum = len(seqs)
	haploSize = len(seqs[0])
	return(haploNum,haploSize,seqs,f1HaploNum,eqVert)

def calcOrderedFrequencies(haploNum,haploSize,seqs): # calculate the ordered frequencies    
	import math
	freqCount = np.zeros((haploSize, 5))
	productVector = np.zeros((haploSize, 5))
	entropyVector = np.zeros(haploSize)
	for read in seqs:
		for pos in range(haploSize):
			if read[pos] == 'A':
				freqCount[pos, 0] = freqCount[pos, 0] + 1
			elif read[pos] == 'C':
				freqCount[pos, 1] = freqCount[pos, 1] + 1
			elif read[pos] == 'G':
				freqCount[pos, 2] = freqCount[pos, 2] + 1
			elif read[pos] == 'T':
				freqCount[pos, 3] = freqCount[pos, 3] + 1
			elif read[pos] == '-':
				freqCount[pos, 4] = freqCount[pos, 4] + 1
	freqRel = np.divide(freqCount, haploNum, dtype = float)
	for pos in range(haploSize):
		for i in range(5):
			freqPos = freqRel[pos, i]
			if freqPos > 0:
				logFreqRel = math.log(freqPos, 2)
				productVector[pos, i] = -1*(np.multiply(freqPos, logFreqRel, dtype = float))                
	#hVector = np.sum(productVector, axis = 1)
	return np.sum(productVector, axis = 1)

def orderPositions(hVector,seqs,haploSize): # order positions
	invH =  np.multiply(-1, hVector, dtype = float)
	ordH = np.argsort(invH)
	# reorder the sequences by entropy
	ordSeqs = []
	for i in seqs:
		newOne = ''
		for p in range(haploSize):
			newOne = ''.join([newOne, i[ordH[p]]])
		ordSeqs.append(newOne)
	return ordSeqs

def calcDistances(haploNum,haploSize,ordSeqs): # Calculate distances
	from itertools import izip

	chunkStart = range(0, haploSize, 1)
	compNum = haploSize
	compList = range(haploNum)
	t = 0
	m = haploSize
	adjMatrix = np.zeros((haploNum, haploNum))
	kStepList = []
	while compNum > 1:
		t = t + 1
		# Check each query sequence 
		for r1 in range(haploNum-1):
			haplotype1 = ordSeqs[r1]
			for r2 in range(r1+1, haploNum):
				if compList[r1] != compList[r2]: 
					haplotype2 = ordSeqs[r2]
					tempDist = 0
					for a, b in izip(haplotype1, haplotype2):
						if a != b:
							tempDist = tempDist + 1
							if tempDist > t:
								break
					if tempDist == t: 
						adjMatrix[r1][r2] = 1
						kStepList.append([r1, r2, t])
		# Recalculate components
		sparseMatrix = csgraph_from_dense(adjMatrix)
		connected = connected_components(sparseMatrix, directed=False, connection='weak', return_labels=True)
		compNum = connected[0]
		compList = connected[1]
	return [kStepList,t]

def makeNewEdgeList(kStepList,seqs): # made up edges need unique ID tags, now they are just all 0
	already=0
	edgeListNew=[]
	fakeEdges=[]
	fakeSeqs={}
	for element in kStepList:
		dist=element[2]
		if dist==1:
			edgeListNew.append(element)
		else:
			#element[1] -> id for seq1
			#seqs[element[1]] -> sequence for seq1
			s1=seqs[element[0]]
			s2=seqs[element[1]]
			#assume sequences are same length
			id=0
			while dist>1:
				id+=1
				for nucl in range(len(s1)):
					if s1[nucl]!=s2[nucl]:
						takeseq1=s1[:nucl+1]+s2[nucl+1:]
						takeseq2=s1[:nucl]+s2[nucl:]
						if takeseq1!=s1:
							newseq=takeseq1
						else:
							newseq=takeseq2
				dist-=1
				s1=newseq
				if newseq in seqs:
					print("WARNING: intermediate sequence detected in input, skipping")
				elif newseq in fakeSeqs.values():
					already+=1
				else:
					newseqid=str(element[0])+"->"+str(element[1])+"_"+str(id)
					fakeEdges.append([newseqid,element[1],1])
					fakeSeqs[newseqid]=newseq
	edgeListNew=edgeListNew.extend(fakeEdges)
	print("%i intermediate sequences had already been added to set and were not added again" %already)
	#make sequence dictionary
	finalSeqs={}
	for item in range(len(seqs)):
		finalSeqs[item]=seqs[item]
	finalSeqs.update(fakeSeqs)
	return(finalSeqs)

def calcDistanceMatrix(finalSeqs):
	from itertools import izip
	matSize=len(finalSeqs)
	indexList=finalSeqs.keys()
	D=np.zeros([matSize,matSize],int)
	for x in range(1,matSize):
		for y in range(0,x):
			seq1=finalSeqs[indexList[x]]
			seq2=finalSeqs[indexList[y]]
			pairDist=0
			for a,b in izip(seq1,seq2):
				if a!=b:
					pairDist+=1
			D[x][y]=pairDist
			D[y][x]=pairDist
	# outMat=outDirec+"/"+output+".data"
	# np.savetxt(outMat,D,fmt='%i') #output can be changed here
	return D

def calcMapVal(seqs):#calculate_mapval
	numSeq=len(seqs)
	numPos=len(seqs[0])
	mapval=0;a=0;c=0;g=0;t=0;blank=0
	for pos in range(numPos,0,-1):
		for seq in seqs:
			# print(pos)
			# print(len(seq))
			#print(seq)
			if seq[pos-1].lower()=="a":
				a+=1
			elif seq[pos-1].lower()=="c":
				c+=1
			elif seq[pos-1].lower()=="t":
				t+=1
			elif seq[pos-1].lower()=="g":
				g+=1
			else:
				blank+=1
				print("blank")
		consensus=max(a,c,t,g,blank)
		if a+c+t+g+blank<numSeq or consensus>=numSeq or blank>0:
			a=0;c=0;g=0;t=0;blank=0;
		else:
			mapval+=1;a=0;c=0;g=0;t=0;blank=0;
	return mapval

def simulEvol(D_all,nseq1,nseq2,len_eff): #viral evolution simulator, probably very broken
	maxPopl=10**12
	timeInter=1111
	mutprob=.01
	evolved=False
	sz=len(D_all);
	Q=np.zeros([sz,sz],float)
	for u in range(sz):
		Q[u,u]=(mutprob/3)**D_all[u,u]*(1-mutprob)**(len_eff-D_all[u,u])
		for v in range(u+1,sz):
			Q[u,v]=D_all[u,v]*(mutprob/3)**(D_all[u,v])*(1-mutprob)**(len_eff-D_all[u,v])
			Q[v,u]=Q[u,v]
	x=np.zeros([sz,timeInter],float)
	x[0:nseq1-1,0]=1
	E=np.eye(sz)
	for t in range(1,timeInter):
		x[:,t]=(1-sum(x[:,t-1])/maxPopl) * np.dot((E+Q),x[:,t-1])
		test=nseq1
		while x[test,t]>=1: #np.where
			test+=1
			if test==nseq1+nseq2:
				print("evolution complete")
				evolved=True
				break
		if evolved==True:
			time=t+1
			break
	if evolved==False:
		time=1200
	return (time)
			
def reduceMat(DSamp,comp):
    DSsz=len(DSamp)
    cosz=len(comp)
    out=np.zeros([cosz,cosz])
    for i in range(cosz):
        for j in range(cosz):
            out[i,j]=DSamp[comp[i],comp[j]]
	return out

def findDSampDir(DSamp):
	AMSamp_dir=(DSamp<=np.transpose(DSamp))
	return DSamp*AMSamp_dir
	
def findTransNetMCMC5(DSamp_comp): #incomplete
	nSamp=len(DSamp)
	nVertTree=2*nSamp-1
	DSamp_dir=findDSampDir(DSamp)
				
def analyze2FilesWrapper(input1,input2):
	import os  
	seqs = []
	seqNames = []
	fileName, fileExtension = os.path.splitext(input1)
	fileName2, fileExtension2 = os.path.splitext(input2)
	statement_f1='File 1: %s' %input1;statement_f2='File 2: %s' %input2;
	print("Parsing the following fasta files...")
	print(statement_f1)
	print(statement_f2)
	[haploNum,haploSize,seqs,f1HaploNum,eqVert]=parseInput(input1,input2)
	#f1HaploNum -> # of seqeunces in file 1
	#haploNum -> # of sequences in both files
	f2HaploNum=haploNum-f1HaploNum
	print("Calculating ordered frequencies...")
	len_eff=calcMapVal(seqs) 
	hVector=calcOrderedFrequencies(haploNum,haploSize,seqs)
	print("Ordering positions...")
	ordSeqs=orderPositions(hVector,seqs,haploSize)
	print("Calculating kSteplists...")
	[kStepList,t]=calcDistances(haploNum,haploSize,ordSeqs)
	print("Analyzing links...")
	finalSeqs=makeNewEdgeList(kStepList,seqs)
	print("Computing pairwise hamming distance...")
	D_all=calcDistanceMatrix(finalSeqs)
	statement_f1='File 1: %s' %input1;statement_f2='File 2: %s' %input2;
	return (D_all,f1HaploNum,f2HaploNum,eqVert,len_eff)
	

def main(inputs):
	import time

	startTime = time.time() 
	# f1=sys.argv[1]
	# f2=sys.argv[2]
	# (D_all,nseq1,nseq2,eqVert,len_eff)=analyze2FilesWrapper(f1,f2)
	# evolTime=simulEvol(D_all,nseq1,nseq2,len_eff)
	# print(evolTime)
	del sys.argv[0]
	numFiles=len(sys.argv)
	storage=np.zeros([numFiles,numFiles],int)
	for i1 in range(numFiles):
		f1=sys.argv[i1]
		for i2 in range(numFiles):
			f2=sys.argv[i2]
			if i1!=i2:
				(D_all,nseq1,nseq2,eqVert,len_eff)=analyze2FilesWrapper(f1,f2)
				evolTime=simulEvol(D_all,nseq1,nseq2,len_eff)
				storage[i1,i2]=evolTime
	AMSamp=np.eye(numFiles)
	for u in range (numFiles):
		for v in range(u+1,numFiles):
			if storage[u,v]<=storage[v,u]:
				if storage[u,v]<1100:
					AMSamp[u,v]=1
	sparseMatrix = csgraph_from_dense(AMSamp)
	connected = connected_components(sparseMatrix, directed=False, connection='weak', return_labels=True)
	S = connected[0]
	C = connected[1]
	for c in range(S):
		comp=np.where(C==c-1)
		if len(comp)>1:
			DSamp_comp=reduceMat(DSamp,comp)
			[transNetsComp,TransTreesComp]=findTransNetMCMC5(DSamp_comp)
	endTime = time.time()
	workTime =  endTime - startTime
	statement_time='Analysis took %.3f seconds' % workTime;
	print(statement_time)
	
if __name__ == '__main__':
	print("Initializing...") 
	import numpy as np
	import sys
	from scipy.sparse.csgraph import connected_components
	from scipy.sparse.csgraph import csgraph_from_dense
	main(sys.argv)
