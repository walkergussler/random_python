from scipy.io import savemat

f=open('mash_hamm.csv')
lines=f.readlines()
f.close()
mash=[]
mind=[]
for line in lines:
    s=line.split(',')
    mash.append(s[2])
    mind.append(s[3])
del mash[0]
del mind[0]

with open("tmp.mat","w") as f:
    savemat(f,{'mash':mash,'mind':mind})
