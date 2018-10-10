#################################################
#created the 20/04/2018 12:57 by Alexis Blanchet#
#################################################
#-*- coding: utf-8 -*-
'''
Stack des différents algorithms
Automatique thresholding avec une logistic Regression
'''

'''
Améliorations possibles:

'''
import warnings
warnings.filterwarnings('ignore')
#################################################
###########        Imports      #################
#################################################
import sys
import os
from sklearn import linear_model
import pandas as pd
import numpy as np
import pickle
from sklearn.externals import joblib
import numpy as np
from itertools import cycle
from sklearn import svm, datasets
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import label_binarize
from sklearn.multiclass import OneVsRestClassifier
from scipy import interp
from ast import literal_eval
import time
#################################################
########### Global variables ####################
#################################################
PATH_IN = '/home/alexis/Bureau/finalproject/DatasIn/RTS/'
PATH_SCRIPT = '/home/alexis/Bureau/finalproject/scripts/'
PATH_OUT = '/home/alexis/Bureau/finalproject/Datas/'
LOG = "log.txt"
#################################################
########### Important functions #################
#################################################
def get_path():
    datas = pd.read_csv('path.csv')
    return datas['PathtoTempDatas'].values[0],datas['PathtoScripts'].values[0],datas['PathtoTempDatas'].values[0]
def Report(error):
    with open(LOG,'a+') as file:
        file.write(str(error)+' \n')
        print(str(error))

def plot_res(df,predict,y,h = [3,27],threshold=0.5):
    x = df.values
    x = x[(h[0]-3)*60:(h[1]-3)*60]
    t= [(i+3)/60+h[0] for i in range(len(x))]

    pred = list([1 if i[-1]>threshold else 0 for i in predict])
    pred = pred[(h[0]-3)*60:(h[1]-3)*60]
    y = y[(h[0]-3)*60:(h[1]-3)*60]
    tp = np.sum([z*x for z,x in zip(pred,y)])
    fp = np.sum([np.clip(z-x,0,1) for z,x in zip(pred,y)])
    fn = np.sum([np.clip(z-x,0,1) for z,x in zip(y,pred)])



    beta = 2
    p = tp/np.sum(pred)
    r = tp/np.sum(y)
    beta_squared = beta ** 2
    f = (beta_squared + 1) * (p * r) / (beta_squared * p + r)
    Report('---------------------------------<Stack>------------------------------------------------------------')
    Report("||precison: "+str(p)+"||recall: "+str(r)+"||fbeta: "+str(f))
    tp,fp,fn = mesure(pred,y)
    beta = 2
    p = tp/(tp+fp)
    r = tp/(tp+fn)
    beta_squared = beta ** 2
    f = (beta_squared + 1) * (p * r) / (beta_squared * p + r)

    Report("||precison: "+str(p)+"||recall: "+str(r)+"||fbeta: "+str(f))
    Report('----------------------------------------------------------------------------------------------------')

def mesure(y_pred,y_true):
    TP = 0
    FP = 0
    FN = 0
    for i in range(len(y_pred)-1):
        i = i+1
        if(y_pred[i] == 1):
            if(sum(y_true[i-1:i+1])>0):
                TP += 1
            else:
                FP += 1
    for i in range(len(y_true)-1):
        i = i+1
        if(y_true[i] == 1):
            if(sum(y_pred[i-1:i+1])>0):
                pass
            else:
                FN += 1
    return TP,FP,FN

def model_fit(X,y):
    clf = Classifier()
    clf.fit(X,y)
    return clf

def find_index(l,v):
    res = []
    for i, j in enumerate(l):
        if(j == v):
            res.append(i)
    return res

#################################################
########### main with options ###################
#################################################


def main(argv):
    global PATH_IN,PATH_SCRIPT,PATH_OUT
    PATH_IN,PATH_SCRIPT,PATH_OUT = get_path()
    #### get files names ###
    names = pd.read_csv('files.csv')
    fileX_train = literal_eval(names['fileX_train'][0])
    fileY_train = literal_eval(names['fileY_train'][0])

    fileX_valid =literal_eval(names['fileX_valid'][0])
    fileY_valid = literal_eval(names['fileY_valid'][0])
    fileX =literal_eval(names['fileX_test'][0])
    fileY = literal_eval(names['fileY_test'][0])
    y = pd.DataFrame()
    for filex,filey in zip(fileX,fileY  ):
        y_ = pd.read_csv(filey)
        y_train = y_['label'][3:]
        y = pd.concat([y,y_train])
    Y = y.values.reshape(-1, 1)

    y_valid = pd.DataFrame()
    for filex,filey in zip(fileX_valid,fileY_valid):
        y_ = pd.read_csv(filey)
        y_train = y_['label'][3:]
        y_valid = pd.concat([y_valid,y_train])
    Y_valid = y_valid.values.reshape(-1, 1)
    if(len(argv)==0):
        argv = [0.45]
    if(str(argv[0]) == 'trainclf'):
        Report('training models ...')
        Report("LGBM")
        #l1 = os.system("python "+PATH_SCRIPT+"LightGBM.py 0.5")
        Report("catboost")
        l2 = os.system("python "+PATH_SCRIPT+"CatBoost.py 0.5")
        Report("classic")
        l3 = os.system("python "+PATH_SCRIPT+"SVC.py 0.5")
        Report("NN")
        l4 = os.system("python "+PATH_SCRIPT+"NeuralNetwork.py 0.5")#TODO: add NN,LSTM back
        Report("XGB")
        l5 = os.system("python "+PATH_SCRIPT+"XgBoost.py 0.5")
        Report("KNN")
        l6 = os.system("python "+PATH_SCRIPT+"KNN.py 0.5")
        Report("LSTM")
        l7 = os.system("python "+PATH_SCRIPT+"LSTM.py 0.5")#TODO: add NN,LSTM back
        return 0
    if(str(argv[0]) == 'trainlogreg'):
        Report('training logistic regression...')
        l1 = pd.read_csv("lightGBM_valid.csv")
        l2 = pd.read_csv("catboost_valid.csv")
        l3 = pd.read_csv("SVC_valid.csv")
        l4 = pd.read_csv("NN_valid.csv")#TODO: add NN,LSTM back
        l5 = pd.read_csv("xgb_valid.csv")
        l6 = pd.read_csv("KNN_valid.csv")
        l7 = pd.read_csv("LSTM_valid.csv")#TODO: add NN,LSTM back


        X_valid = pd.concat([l2,l3,l4,l5,l6,l7], axis=1).values #TODO: add NN,LSTM back
        print(X_valid.shape)
        for i in [0.001]:
            Report("C="+str(i))
            np.random.seed(7)
            logistic = linear_model.LogisticRegression(C=i,class_weight='balanced',penalty='l2')
            logistic.fit(X_valid, Y_valid)
            pickle.dump(logistic, open('model/logistic_regression.sav', 'wb'))
            time.sleep(20)
            os.system("python "+PATH_SCRIPT+"Stack.py")
        return 0
    if(str(argv[0]) == 'trainall'):
        t = time.time()
        os.system("python "+PATH_SCRIPT+"Stack.py trainclf")
        Report("temps d'entraînement de l'ensemble des modèles: "+str(time.time()-t))
        os.system("python "+PATH_SCRIPT+"Stack.py trainlogreg")

    else:
        T = argv
        Report('Scoring...')
        #l1 = pd.read_csv("lightGBM.csv")
        l2 = pd.read_csv("catboost.csv")
        l3 = pd.read_csv("SVC.csv")
        l4 = pd.read_csv("NN.csv")#TODO: add NN,LSTM back
        l5 = pd.read_csv("xgb.csv")
        l6 = pd.read_csv("KNN.csv")
        l7 = pd.read_csv("LSTM.csv")#TODO: add NN,LSTM back

        X = pd.concat([l2,l3,l4,l5,l6,l7], axis=1).values #TODO: add NN,LSTM back
        np.random.seed(7)
        logistic = pickle.load(open('model/logistic_regression.sav', 'rb'))
        np.random.seed(7)
        Predict = logistic.predict_proba(X)
        for j in T:
            Report("Threshold="+str(j))
            #for h in [[3,27],[6,13],[13,20],[20,27],[6,24],[10,13],[12,15],[6,11],[13,16],[14,18],[16,19],[19,22],[20,23],[23,27],[10,18]]:
                #Report(h)
                #plot_res(pd.read_csv(fileX[0])['t'],Predict,Y,h,threshold = float(j))
                #pred = list([1 if i[-1]>float(j) else 0 for i in Predict])

            plot_res(pd.read_csv(fileX[0])['t'],Predict,Y,threshold = float(j))
    return ("process achevé sans erreures")


if __name__ == "__main__":
    # execute only if run as a script
    main(sys.argv[1:])
