#!/usr/bin/python
import sys

def lcs(S,T):
    m = len(S)
    n = len(T)
    counter = [[0]*(n+1) for x in range(m+1)]
    longest = 0
    lcs_set = set()
    for i in range(m):
        for j in range(n):
            if S[i] == T[j]:
                c = counter[i][j] + 1
                counter[i+1][j+1] = c
                if c > longest:
                    lcs_set = set()
                    longest = c
                    lcs_set.add(S[i-c+1:i+1])
                elif c == longest:
                    lcs_set.add(S[i-c+1:i+1])

 #   if len(list(lcs_set)[0])>10:
#        print(c1,c2)
 #       print(lcs_set)
    return max(lcs_set,key=len)

myarray=[]
theirarray=[]
f=open(sys.argv[1],"r")
lines=f.readlines()
f.close()
for line in lines:
    if line[0]!=">":
        myarray.append(line)

f=open(sys.argv[2],"r")
lines=f.readlines()
f.close()
for line in lines:
    if line[0]!=">":
        theirarray.append(line)

big=""
lenarray=[]
for my in myarray:
    for they in theirarray:
        seq=lcs(my,they)
        lenarray.append(len(seq))
        if len(seq)>len(big):
            big=seq
print(len(big),big)
print(sum(lenarray)/len(lenarray))
