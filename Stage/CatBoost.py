#-*- coding: utf-8 -*-
#################################################
#created the 04/05/2018 09:52 by Alexis Blanchet#
#################################################
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
import plotly
import plotly.graph_objs as go
import plotly.offline as offline
from plotly import tools
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler,StandardScaler
from sklearn.base import BaseEstimator
from catboost import CatBoostClassifier
from sklearn.model_selection import train_test_split
from ast import literal_eval

LOG = "log.txt"
#################################################
########### Global variables ####################
#################################################
'''
Définir ici une classe est beaucoup plus simple...
'''
class Classifier(BaseEstimator):
    def __init__(self):
        pass

    def fit(self, x1, y1,x2,y2):
        self.clf1 = CatBoostClassifier(iterations=2000,learning_rate=0.01, depth=10,metric_period = 50, loss_function='Logloss', eval_metric='Logloss', random_seed=99, od_type='Iter', od_wait=100,class_weights=[1,5])
        self.clf1.fit(x1,y1,verbose=True,eval_set=(x2,y2),use_best_model=True)
        self.clf2 = CatBoostClassifier(iterations=2000,learning_rate=0.001, depth=8,metric_period = 50, loss_function='Logloss', eval_metric='Logloss', random_seed=99, od_type='Iter', od_wait=100,class_weights=[1,5])
        self.clf2.fit(x1,y1,verbose=True,eval_set=(x2,y2),use_best_model=True)
    def predict(self, X):
        return self.clf.predict(X)

    def predict_proba(self, X):
        return np.array([[1-(v[1]+l[1])*0.5,(v[1]+l[1])*0.5] for v,l in zip(self.clf2.predict_proba(X),self.clf1.predict_proba(X))])


#################################################
########### Important functions #################
#################################################
#Duplicate of the def_context but as we use a different set of PATH and only the Report Functon
#redefining them here is a better option.
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
    '''
    Custom mesure to ensure that time is taken in account
    '''
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
    '''
    Used to be a wonderfull plot, thus the name.
    But the boss bullied me into getting ride of it.
    Now it's just a reporting function...what a waste....
    '''
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
    Report('-----------------<CatBoost>---------------------------')
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
    model.clf1.save_model("model/catboostmodel1")
    model.clf2.save_model("model/catboostmodel2")



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
    if (len(argv) == 2):
        fileX = argv[0]
        X_minmax, X_meanvar, t = load_(fileX)
        CatBoost = []
        CatBoost.append(CatBoostClassifier().load_model(fname="model/catboostmodel1"))
        CatBoost.append(CatBoostClassifier().load_model(fname="model/catboostmodel2"))
        res = pd.DataFrame([[1-(v[1]+l[1])*0.5,(v[1]+l[1])*0.5] for v,l in zip(CatBoost[0].predict_proba(X_meanvar),CatBoost[1].predict_proba(X_meanvar))])
        res.to_csv(str(fileX.split('.')[0])+'_temp_cat.csv',index=False)
        return res
    else:
        if(len(argv)==0):
            argv = [0.2]
        THRESHOLD = float(argv[0])
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

        pred_valid = model.predict_proba(X_valid)
        res_valid = pd.DataFrame(pred_valid)
        res_valid.to_csv('catboost_valid.csv',index=False)

        res = pd.DataFrame(pred)
        res.to_csv('catboost.csv',index=False)
        save_model(model)
        return res

if __name__ == "__main__":
    # execute only if run as a script
    main(sys.argv[1:])
