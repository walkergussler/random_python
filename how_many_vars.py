import numpy as np
from sklearn.model_selection import train_test_split, cross_val_score, GroupKFold
from sklearn.svm import SVC
# from sklearn.feature_selection import SelectKBest, chi2, SelectFromModel, RFE
from sklearn.metrics import accuracy_score, confusion_matrix
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier
from sklearn.linear_model import LogisticRegression
# from skrebate import ReliefF, MultiSURF
# from scipy.stats import ttest_ind
# from matplotlib import pyplot
from sklearn.naive_bayes import GaussianNB


def ml_data_parser(file):
    data_end=-3
    groups_col=-2
    c=0
    X=[]
    y=[]
    groups=[]
    with open(file) as f:
        for line in f.readlines():
            c+=1
            if c!=1:
                s=line.strip().split(',') 
                xx=s[1:data_end]
                # for i in xx:
                    # print(i,type(i))
                # print(xx)
                X.append(list(map(float,xx)))
                y.append(int(s[data_end]))
                groups.append(int(s[groups_col]))
            else:
                s=line.strip().split(',')
                names=s[1:data_end]
                groupname=s[groups_col]
                print(names)
                exit()
    return np.array(X),np.array(y),names,groups,groupname


X,y,names,groups,_=ml_data_parser('30_data.csv')
print('i,random_forest,svm,naiive_bayes,logistic_reg,extra_random_forest')
cv_num=10

for i in range(1,len(names)):
    new_x=X[:,:i]
    # f=[np.shape(new_x)[0]]
    # f.append(np.shape(new_x)[1])
    # for j in range(i):
        # f.append(names[j])
    # print(','.join(map(str,f)))
    # x_train,x_test,y_train,y_test=train_test_split(new_x,y,test_size=.33)
    group_kfold=GroupKFold(n_splits=cv_num)
    group_kfold.get_n_splits(new_x,y,groups)
    xlen=len(names)
    rfacc=[]
    svmacc=[]
    nbacc=[]
    lracc=[]
    extrarfacc=[]
    for train_index,test_index in group_kfold.split(new_x,y,groups):
        x_train=[]
        x_test=[]
        y_train=[]
        y_test=[]
        # print(train_index) 
        # print(test_index)
        for id in train_index:
            x_train.append(new_x[id])
        for id in test_index:
            x_test.append(new_x[id])
        for id in train_index:
            y_train.append(y[id])
        for id in test_index:
            y_test.append(y[id])    
        clf0=RandomForestClassifier(n_estimators=100)
        clf1=SVC(kernel='linear')
        clf2=GaussianNB()
        clf3=LogisticRegression(solver='liblinear')
        clf4=ExtraTreesClassifier()
        clf0.fit(x_train,y_train)
        clf1.fit(x_train,y_train)
        clf2.fit(x_train,y_train)
        clf3.fit(x_train,y_train) 
        clf4.fit(x_train,y_train) 
        s0=np.mean(cross_val_score(clf0,new_x,y,cv=cv_num))
        s1=np.mean(cross_val_score(clf1,new_x,y,cv=cv_num))
        s2=np.mean(cross_val_score(clf2,new_x,y,cv=cv_num))
        s3=np.mean(cross_val_score(clf3,new_x,y,cv=cv_num))
        s4=np.mean(cross_val_score(clf4,new_x,y,cv=cv_num))
        rfacc.append(s0)
        svmacc.append(s1)
        nbacc.append(s2)
        lracc.append(s3)
        extrarfacc.append(s4)
    
    # storage.append(s)
    s=[i,np.mean(rfacc),np.mean(svmacc),np.mean(nbacc),np.mean(lracc),np.mean(extrarfacc)]
    print(','.join(map(str,s)))