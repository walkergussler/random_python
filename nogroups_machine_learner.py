from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier
from sklearn.model_selection import train_test_split, GroupKFold, cross_val_score
import numpy as np
from sklearn.feature_selection import SelectKBest, chi2, SelectFromModel, RFE, VarianceThreshold
from seaborn import heatmap
import matplotlib.pyplot as plt 
import pandas as pd
from sklearn.metrics import accuracy_score, confusion_matrix
from sklearn.naive_bayes import GaussianNB                                                                                                                                                                         
from sklearn.linear_model import LinearRegression, LogisticRegression
from copy import copy
data_end=-3
# sfdic={0:'regular',1:'1step'}
# sourcefiledic={0:'normalized_fulldata.csv',1:'normalized_1step.csv'}
def highest_x(vec,a):
    out=np.zeros(len(vec))
    toy_vec=[]
    for x in vec:
        toy_vec.append(x)
    for i in range(a):
        m_i=vec.index(max(toy_vec))
        toy_vec[m_i]=0
        out[m_i]=1
    return out

modeldic={
    0:'random forest',
    1:'svm',
    2:'naive_bayes',
    3:'logistic regression'
}
for modelid in range(4):
    print("=======")
    X=[]
    y=[]
    groups=[]
    names=[]
    c=0
    with open('74_data_fixed.csv') as f:
        for line in f.readlines():
            c+=1
            if c!=1:
                s=line.strip().split(',') 
                X.append(list(map(float,s[1:data_end])))
                y.append(int(s[data_end]))
                # groups.append(int(s[groups_col]))                    
            else:
                s=line.strip().split(',')
                names=s[1:data_end]
                # print(sfdic[i]+','+modeldic[modelid])
    ######################################## for groupkfold
    # group_kfold=GroupKFold(n_splits=3)
    # group_kfold.get_n_splits(X,y,groups)
    # for train_index,test_index in group_kfold.split(X,y,groups):
        # X_train=[]
        # X_test=[]
        # y_train=[]
        # y_test=[]
        # for id in train_index:
            # X_train.append(X[id])
        # for id in test_index:
            # X_test.append(X[id])
        # for id in train_index:
            # y_train.append(y[id])
        # for id in test_index:
            # y_test.append(y[id])
    ######################################## for random split
    clf=ExtraTreesClassifier(n_estimators=100).fit(X,y)
    feat_imp=clf.feature_importances_
    model=SelectFromModel(clf,prefit=True)
    X_new=model.transform(X)
    _,a=X_new.shape
    highest=highest_x(list(feat_imp),a)
    for i in range(len(names)):
        if highest[i]:
            print(names[i])
    # sel=VarianceThreshold().fit_transform(X)
    # print(sel)
    exit()
    if True:
        X_train,X_test,y_train,y_test=train_test_split(X,y)    
    ######################################## 
        if modelid==0:
            clf=RandomForestClassifier(n_estimators=100)
        elif modelid==1:
            clf=SVC(kernel='linear')
        elif modelid==2:
            clf=GaussianNB()
        elif modelid==3:
            clf=LogisticRegression(solver='liblinear')
        clf.fit(X_train,y_train)
        qwer=clf.score(X_test,y_test)
        print('accuracy='+','+str(qwer)+'\n')
    scores=cross_val_score(clf,X,y,cv=10)
    print('cv acc='+str(np.mean(scores)))
from sys import exit
exit()

