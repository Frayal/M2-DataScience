#################################################
#created the 04/05/2018 09:52 by Alexis Blanchet#
#################################################
#-*- coding: utf-8 -*-
'''

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
import numpy as np
import pandas as pd
import scipy.stats
from sklearn.preprocessing import MinMaxScaler,StandardScaler
from sklearn.base import BaseEstimator
from sklearn.ensemble import RandomForestClassifier
import xgboost as xgb
from sklearn.model_selection import train_test_split
import pickle
from ast import literal_eval
#################################################
########### Global variables ####################
#################################################
### XGB modeling
params = {'eta': 0.001,
          'max_depth': 10,
          'subsample': 0.9,
          'colsample_bytree': 1,
          'colsample_bylevel':1,
          'min_child_weight':1,
          'alpha':1,
          'objective': 'binary:logistic',
          'eval_metric': 'logloss',
          'seed': 99,
          'silent': False,
         'scale_pos_weight':20,
         'silent' : 1}
params2 = {'eta': 0.001,
          'max_depth': 8,
          'subsample': 0.9,
          'colsample_bytree': 1,
          'colsample_bylevel':1,
          'min_child_weight':1,
          'alpha':1,
          'objective': 'binary:logistic',
          'eval_metric': 'logloss',
          'seed': 99,
          'silent': False,
         'scale_pos_weight':20,
          'silent' : 1}

######################################################
class Classifier(BaseEstimator):
    def __init__(self):
        pass

    def fit(self, x1, y1,x2,y2):
        watchlist = [(xgb.DMatrix(x1, y1), 'train'), (xgb.DMatrix(x2, y2), 'valid')]
        self.clf1 = xgb.train(params, xgb.DMatrix(x1, y1), 5000,  watchlist, maximize = False,verbose_eval=200, early_stopping_rounds=300)
        self.clf2 = xgb.train(params2, xgb.DMatrix(x1, y1), 5000,  watchlist, maximize = False,verbose_eval=200, early_stopping_rounds=300)


    def predict(self, X):
        return self.clf.predict(X)

    def predict_proba(self, X):
        res1 = self.clf1.predict(xgb.DMatrix(X), ntree_limit=self.clf1.best_ntree_limit)
        res2 = self.clf2.predict(xgb.DMatrix(X), ntree_limit=self.clf2.best_ntree_limit)
        res = [(r1+r2)*0.5 for r1,r2 in zip(res1,res2)]
        return np.array([ [1-c,c] for c in res])

PATH_IN = '/home/alexis/Bureau/finalproject/DatasIn/RTS/'
PATH_SCRIPT = '/home/alexis/Bureau/finalproject/scripts/'
PATH_OUT = '/home/alexis/Bureau/finalproject/Datas/'
LOG = "log.txt"

#################################################
########### Important functions #################
#################################################
def Report(error):
    with open(LOG,'a+') as file:
        file.write(str(error)+' \n')
        print(str(error))
def get_path():
    datas = pd.read_csv('path.csv')
    return datas['PathtoDatasIn'].values[0],datas['PathtoScripts'].values[0],datas['PathtoTempDatas'].values[0]

def load(fileX,fileY):
    X = pd.DataFrame()
    y = pd.DataFrame()
    for filex,filey in zip(fileX,fileY):
        df = pd.read_csv(filex)
        y_ = pd.read_csv(filey)
        df = df.replace([np.inf, -np.inf], np.nan)
        df = df.fillna(1)
        X_train = df
        y_train = y_['label'][3:]
        X = pd.concat([X,X_train])
        y = pd.concat([y,y_train])
    t = X.index.values

    scaler = StandardScaler()
    X = scaler.fit_transform(X.values)
    return  X,y.values.reshape(-1, 1),t




def model_fit(X1,y1,X2,y2):
    clf = Classifier()
    clf.fit(X1,[Y[0] for Y in y1],X2,[Y[0] for Y in y2])
    return clf

def find_index(l,v):
    res = []
    for i, j in enumerate(l):
        if(j == v):
            res.append(i)
    return res

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

def plot_res(df,pred,y):
    x = df
    t= [i/60 +3 for i in range(len(x))]
    tp = np.sum([z*x for z,x in zip(pred,y)])
    fp = np.sum([np.clip(z-x,0,1) for z,x in zip(pred,y)])
    fn = np.sum([np.clip(z-x,0,1) for z,x in zip(y,pred)])

    beta = 2
    p = tp/np.sum(pred)
    r = tp/np.sum(y)
    beta_squared = beta ** 2
    f = (beta_squared + 1) * (p * r) / (beta_squared * p + r)
    Report('-------------------<XGBoost>----------------------')
    Report("|| precison: "+str(p)+"|| recall: "+str(r)+"|| fbeta: "+str(f))

    tp,fp,fn = mesure(pred,y)
    beta = 2
    p = tp/(tp+fp)
    r = tp/(tp+fn)
    beta_squared = beta ** 2
    f = (beta_squared + 1) * (p * r) / (beta_squared * p + r)


    Report("|| precison: "+str(p)+"|| recall: "+str(r)+"|| fbeta: "+str(f))
    Report('--------------------------------------------------')


def save_model(model):
    pickle.dump(model.clf1, open("model/XGB1.pickle.dat", "wb"))
    pickle.dump(model.clf2, open("model/XGB2.pickle.dat", "wb"))

def load_(fileX):
    df = pd.read_csv(fileX)
    df = df.replace([np.inf, -np.inf], np.nan)
    df = df.fillna(1)
    X_train = df.values
    t = df.index.values

    scaler = MinMaxScaler(feature_range=(0, 1))#StandardScaler()
    s = StandardScaler()
    X_train_minmax = scaler.fit_transform(X_train)
    X_train_meanvar = s.fit_transform(X_train)
    return  X_train_minmax,X_train_meanvar,t


#################################################
########### main with options ###################
#################################################


def main(argv):
    global PATH_IN, PATH_SCRIPT, PATH_OUT
    PATH_IN, PATH_SCRIPT, PATH_OUT = get_path()
    if(len(argv)==2):
        fileX = argv[0]
        X_minmax, X_meanvar, t = load_(fileX)
        XGB = []
        XGB.append(pickle.load(open("model/XGB1.pickle.dat", "rb")))
        XGB.append(pickle.load(open("model/XGB2.pickle.dat", "rb")))
        res1 = XGB[0].predict(xgb.DMatrix(X_meanvar), ntree_limit=XGB[0].best_ntree_limit)
        res2 = XGB[1].predict(xgb.DMatrix(X_meanvar), ntree_limit=XGB[1].best_ntree_limit)
        res = pd.DataFrame([[1 - (r1 + r2) * 0.5, (r1 + r2) * 0.5] for r1, r2 in zip(res1, res2)])
        res.to_csv(str(fileX.split('.')[0])+'_temp_XGB.csv',index=False)
        return res
    else:
        if(len(argv)==0):
            argv = [0.04]
        THRESHOLD = float(argv)
        #### get files names ###
        names = pd.read_csv('files.csv')
        fileX_train = literal_eval(names['fileX_train'][0])
        fileY_train = literal_eval(names['fileY_train'][0])

        fileX_valid =literal_eval(names['fileX_valid'][0])
        fileY_valid = literal_eval(names['fileY_valid'][0])
        fileX_test =literal_eval(names['fileX_test'][0])
        fileY_test = literal_eval(names['fileY_test'][0])
        X_train,Y_train,_ = load(fileX_train,fileY_train)
        X_valid,Y_valid,_ = load(fileX_valid,fileY_valid)
        X_test,Y_test,t = load(fileX_test,fileY_test)

        model = model_fit(X_train,Y_train,X_valid,Y_valid)
        pred = model.predict_proba(X_test)
        testPredict = list([1 if i[1]>THRESHOLD else 0 for i in pred])
        # plot results
        plot_res(t,testPredict,Y_test)

        res = pd.DataFrame(pred)
        res.to_csv('xgb.csv',index=False)

        pred_valid = model.predict_proba(X_valid)
        res_valid = pd.DataFrame(pred_valid)
        res_valid.to_csv('xgb_valid.csv',index=False)
        save_model(model)
        return res

if __name__ == "__main__":
    # execute only if run as a script
    main(sys.argv[1:])
