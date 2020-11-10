


def updateVals1D(inVector,indexVector,updateValue): 
#replaces matlab's inVector(indexVector)=updateValue (returns outVector==inVector)
#with outVector=updateVals1D(inVector,indexVector,updateValue)
	outVector=inVector
	ind=0;
	for element in inVector:
		ind+=1
		if ind in indexVector:
			outVector[ind]=updateValue
	return outVector

def linkageold(z): #expects a square numpy matrix
	import numpy as np
	import math
	n=len(z[0])
	m=math.ceil((2*n)**.5)
	Z=np.zeros([m-1,3])
	N=np.zeros([1,2*m-1])
	N[0,0:m]=1
	n=m
	R=np.arange(n)+1
	# for sASDF in range(n-1)
		# s=sASDF+1
	p=np.arange(m-1,1,-1)
	I=np.zeros([m*(m-1)/2,1])
	cumVec=np.cumsum(np.append([1],p))
	outVector=updateVals1D(I,cumVec,1)
	print(outVector)	