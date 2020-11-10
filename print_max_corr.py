
# from sklearn.model_selection import train_test_split
import numpy as np
from sys import argv
# from sklearn.feature_selection import SelectKBest, chi2, SelectFromModel
# from sklearn.metrics import accuracy_score, confusion_matrix
# from sklearn.cluster import DBSCAN
# from skrebate import ReliefF, MultiSURF
# from scipy.stats import ttest_ind
# from matplotlib import pyplot
# from sklearn.decomposition import PCA

def list2str(x):
    return ','.join(map(str,x))


X=[]
y=[]
c=0
with open(argv[1]) as f:
    for line in f.readlines():
        c+=1
        s=line.strip().split(',') 
        if c!=1:            
            x_items=list(map(float,s[1:-2]))
            X.append(x_items)
        else:
            names=s
X=np.transpose(np.array(X))
print(list2str(names))
corrmat=abs(np.corrcoef(X))
corrmat-=np.eye(len(corrmat))
max_corr=np.amax(corrmat)
print(max_corr)
q=np.where(corrmat==max_corr)
#print(names[int(q[0]+1)])
#print(names[int(q[1]+1)])
print(q)
# print(type(q))
# print(corrmat)
# print(names[-1*(q[0][0]+1)])
# print(names[-1*(q[1][0]+1)])

#print(names[np.argmax(corrmat)+1])
#print(names[np.argmax(corrmat[np.argmax(corrmat)])+1])
exit()
import seaborn
import matplotlib.pyplot as plt

plt.figure()
seaborn.clustermap(corrmat)
plt.show()
