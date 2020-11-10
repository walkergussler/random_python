import numpy as np
from sklearn.model_selection import GroupKFold
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import ExtraTreesClassifier, RandomForestClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import SVC
# import tensorflow as tf
# import os
import pandas as pd
# from sklearn.preprocessing import MinMaxScaler

def ml_data_parser(file):
    # print(file)
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
                vars=s[1:data_end]
    return np.array(X),np.array(y),vars,groups 

def get_indices(vars,candidate):
    out=[]
    for i in range(len(vars)):
        name=vars[i]
        if name in candidate:
            out.append(i)
    return out

def list2str(x):
    return ','.join(map(str,x))

def nn(x_train,x_test,y_train,y_test):
    # constants
    learning_rate = 0.001
    training_epochs = 100
    number_of_inputs = 11
    number_of_outputs = 1
    layer_1_nodes = 50
    layer_2_nodes = 100
    layer_3_nodes = 50
    # Input Layer
    with tf.variable_scope('input'):
        X = tf.placeholder(tf.float32, shape=(None, number_of_inputs))
    # Layer 1
    with tf.variable_scope('layer_1',reuse=tf.AUTO_REUSE):
        weights = tf.get_variable("weights1", shape=[number_of_inputs, layer_1_nodes], initializer=tf.contrib.layers.xavier_initializer())
        biases = tf.get_variable(name="biases1", shape=[layer_1_nodes], initializer=tf.zeros_initializer())
        layer_1_output = tf.nn.relu(tf.matmul(X, weights) + biases)
    # Layer 2
    with tf.variable_scope('layer_2',reuse=tf.AUTO_REUSE):
        weights1 = tf.get_variable("weights2", shape=[layer_1_nodes, layer_2_nodes], initializer=tf.contrib.layers.xavier_initializer())
        biases = tf.get_variable(name="biases2", shape=[layer_2_nodes], initializer=tf.zeros_initializer())
        layer_2_output = tf.nn.relu(tf.matmul(layer_1_output, weights1) + biases)
    # Layer 3
    with tf.variable_scope('layer_3',reuse=tf.AUTO_REUSE):
        weights2 = tf.get_variable("weights3", shape=[layer_2_nodes, layer_3_nodes], initializer=tf.contrib.layers.xavier_initializer())
        biases = tf.get_variable(name="biases3", shape=[layer_3_nodes], initializer=tf.zeros_initializer())
        layer_3_output = tf.nn.relu(tf.matmul(layer_2_output, weights2) + biases)
    # Output Layer
    with tf.variable_scope('output',reuse=tf.AUTO_REUSE):
        weights3 = tf.get_variable("weights4", shape=[layer_3_nodes, number_of_outputs], initializer=tf.contrib.layers.xavier_initializer())
        biases = tf.get_variable(name="biases4", shape=[number_of_outputs], initializer=tf.zeros_initializer())
        prediction = tf.matmul(layer_3_output, weights3) + biases
    # Define the cost function of the neural network that will measure prediction accuracy during trainin
    with tf.variable_scope('cost',reuse=tf.AUTO_REUSE):
        Y = tf.placeholder(tf.float32, shape=(None, 1))
        cost = tf.reduce_mean(tf.squared_difference(prediction, Y))
    # Section Three: Define the optimizer function that will be run to optimize the neural network
    with tf.variable_scope('train',reuse=tf.AUTO_REUSE):
        optimizer = tf.train.AdamOptimizer(learning_rate).minimize(cost)
    # Initialize a session so that we can run TensorFlow operations
    with tf.Session() as session:
        # Run the global variable initializer to initialize all variables and layers of the neural network
        session.run(tf.global_variables_initializer())
        # Run the optimizer over and over to train the network.
        # One epoch is one full run through the training data set.
        for epoch in range(training_epochs):
            # Feed in the training data and do one step of neural network training
            y_train=np.array(y_train).reshape(len(y_train),1)
            session.run(optimizer, feed_dict={X: x_train, Y: y_train})
            # Every 5 training steps, log our progress
            # if epoch % 5==0:
                # train_cost=session.run(cost,feed_dict={X:x_train,Y:y_train})
                # test_cost=session.run(cost,feed_dict={X:x_test,Y:y_test})
                # print(epoch,train_cost,test_cost)
        # Training is now complete!
        # print("Training is complete!")
        y_test=np.array(y_test).reshape(len(y_test),1)
        train_cost=session.run(cost,feed_dict={X:x_train,Y:y_train})
        test_cost=session.run(cost,feed_dict={X:x_test,Y:y_test})
        # print(epoch,train_cost,test_cost)
    return test_cost

def group_test(X,y,model,groups):
    # X=pre_x[:,chosen_vars]
    # x_train,x_test,y_train,y_test=train_test_split(new_x,y,test_size=.3) 
    group_kfold=GroupKFold(n_splits=10)
    group_kfold.get_n_splits(X,y,groups)
    acc_arr=[]
    for train_index,test_index in group_kfold.split(X,y,groups):
        X_train=[]
        X_test=[]
        y_train=[]
        y_test=[]
        # print(train_index)
        # print(test_index)
        for id in train_index:
            X_train.append(X[id])
        for id in test_index:
            X_test.append(X[id])
        for id in train_index:
            y_train.append(y[id])
        for id in test_index:
            y_test.append(y[id])

        if model=='svm':
            clf=SVC(gamma='auto').fit(X_train,y_train)
            tmp_score=clf.score(X_test,y_test)
            acc_arr.append(tmp_score)
        elif model=='rf_extra':
            clf=ExtraTreesClassifier(n_estimators=100).fit(X_train,y_train)
            tmp_score=clf.score(X_test,y_test)
            acc_arr.append(tmp_score)
        elif model=='rf':
            clf=RandomForestClassifier(n_estimators=100).fit(X_train,y_train)
            tmp_score=clf.score(X_test,y_test)
            acc_arr.append(tmp_score)
        elif model=='nb':
            clf=GaussianNB().fit(X_train,y_train)
            tmp_score=clf.score(X_test,y_test)
            acc_arr.append(tmp_score)
        elif model=='lr':
            clf=LogisticRegression(solver='lbfgs').fit(X_train,y_train)
            tmp_score=clf.score(X_test,y_test)
            acc_arr.append(tmp_score)
        elif model=='nn':
            print('run nn')
        # score=np.mean(cross_val_score(svm,X,y,cv=10))
        # print(groupdic[groupid],modeldic[modelid],np.shape(X_test),np.shape(X_train))
        # print('accuracy='+','+str(qwer)+'\n')
    return np.mean(acc_arr)

def parse(file):
    data_end=-3
    data = pd.read_csv(file, index_col=0)
    X=data.iloc[:,0:data_end]
    y=data['status'].values
    groups=data['10_hamming'].values
    return X,y,groups

def scaler(X,y):
    x_scaler = MinMaxScaler(feature_range=(0, 1)).fit_transofrm(X)
    y_scaler = MinMaxScaler(feature_range=(0, 1)).fit_transform(y)
    x_scaled=x_scaler.fit_transform(X)
    y_scaled=x_scaler.fit_transform(y)
    return x_scaled,y_scaled

model_list=['rf_extra']
# model_list=['nn']
print(','+",".join(model_list))
files=['finalechlin_fullfile_normalized.csv','finalechlin_onestep_normalized.csv']
for file in files:
    print(file)
    # X,y,groups=parse('file_both.csv')
    X,y,vars_base,groups=ml_data_parser(file) #have to make a version of this file which also includes cherokee stuffs
    X_phacelia=X[:,0].reshape(-1,1)
    X_doris=X[:,1:]
    out_data=[]
    for a in ['whole','phacelia','doris']:
        if a=='whole':
            model_acc=group_test(X,y,'rf_extra',groups)
        elif a=='phacelia':
            model_acc=group_test(X_phacelia,y,'rf_extra',groups)
        elif a=='doris':
            model_acc=group_test(X_doris,y,'rf_extra',groups)
        else:
            model_acc='error'
        print(a,model_acc)
        # out_data.append(model_acc)
    print(list2str(out_data))
