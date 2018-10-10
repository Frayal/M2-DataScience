#################################################
#created the 17/05/2018 14:14 by Alexis Blanchet#
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
import os
import numpy as np
import pandas as pd
import scipy.stats
from sklearn.preprocessing import MinMaxScaler
from sklearn.base import BaseEstimator
from sklearn.ensemble import RandomForestClassifier
import xgboost as xgb
from catboost import CatBoostClassifier
from sklearn.model_selection import train_test_split
from keras import backend as K
from keras.models import Sequential
from keras.layers import Dense, Dropout, Activation
from keras.wrappers.scikit_learn import KerasClassifier
from sklearn import preprocessing
import matplotlib.pyplot as plt
from itertools import cycle
from sklearn import svm, datasets
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import label_binarize
from sklearn.multiclass import OneVsRestClassifier
from scipy import interp
from keras.utils import np_utils
from sklearn.tree import DecisionTreeClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.ensemble import RandomForestClassifier,GradientBoostingClassifier
import pickle
from sklearn import linear_model
import def_context
#################################################
########### Global variables ####################
#################################################

#################################################
########### Important functions #################
#################################################
def get_temp_path():
    datas = pd.read_csv('path.csv')
    return datas['temp_path'].values[0]

def load(fileX):
    df = pd.read_csv(PATH_IN+'hop/'+fileX)
    if('labels' not in df.columns.values):
        def_context.Report('Pas de labels pour le fichier '+str(fileX))
    y = df['labels']
    y = y.fillna(0)
    return df.drop(['labels'],axis=1),y

def load_all(CHAINE):
    X = pd.DataFrame()
    Y = pd.DataFrame()
    files = os.listdir(PATH_IN+'hop/')
    for file in files:
        if(file.split('.')[-1] != 'csv'):
            pass
        elif(file.split('_')[-2] in ['2017-12-20'] or (file.split('_')[-2]).split('-')[0] == '2018'):
            pass

        else:
            def_context.Report(file.split('_')[-2])
            df,y = load(file)
            if(len(df)==1):
                continue
            df = df.replace([np.inf, -np.inf], np.nan)
            df = df.fillna(0)
            X_train = df
            y_train = y
            X = pd.concat([X,X_train])
            Y = pd.concat([Y,y_train])
    return def_context.process(X),Y
######################################################
######################################################
### XGB modeling
params = {'eta': 0.01,
          'max_depth': 10,
          'subsample': 0.9,
          'colsample_bytree': 0.9,
          'colsample_bylevel':0.9,
          'min_child_weight':0.9,
          'alpha':10,
          'objective': 'multi:softprob',
          'eval_metric': 'mlogloss',
          'seed': 99,
          'silent': 1,
         'num_class' : 3,
         }
params2 = {'eta': 0.01,
          'max_depth': 11,
          'subsample': 0.9,
          'colsample_bytree': 0.9,
          'colsample_bylevel':0.9,
          'min_child_weight':0.9,
          'alpha':5,
          'objective': 'multi:softprob',
          'eval_metric': 'mlogloss',
          'seed': 42,
          'silent': 1,
          'num_class' : 3,
         }
######################################################
class Classifier(BaseEstimator):
    def __init__(self):
        pass

    def fit(self,X,y):
        x1, x2, y1, y2 = train_test_split(X.values, y.values, test_size=0.2,random_state = 7)
        watchlist = [(xgb.DMatrix(x1, y1, weight = [int(y)*2+1 for y in y1]), 'train'), (xgb.DMatrix(x2, y2,weight = [int(y)*2+1 for y in y2]), 'valid')]
        self.clf1 = xgb.train(params, (xgb.DMatrix(x1, y1, weight = [int(y)*2+1 for y in y1])), 2000,  watchlist, maximize = False,verbose_eval=500, early_stopping_rounds=300)
        self.clf2 = xgb.train(params2, (xgb.DMatrix(x1, y1, weight = [int(y)*2+1 for y in y1])), 2000,  watchlist, maximize = False,verbose_eval=500, early_stopping_rounds=300)


    def predict(self, X):
        return self.clf.predict(X.values)

    def predict_proba(self, X):
        res1 = self.clf1.predict(xgb.DMatrix(X.values), ntree_limit=self.clf1.best_ntree_limit)
        res2 = self.clf2.predict(xgb.DMatrix(X.values), ntree_limit=self.clf2.best_ntree_limit)
        return np.array([[(a[0]+b[0])/2,(a[1]+b[1])/2,(a[2]+b[2])/2] for a,b in zip(res1,res2)])


#################################################
#################################################
class Classifier2(BaseEstimator):
    def __init__(self):
        pass

    def fit(self, X,y):
        x1, x2, y1, y2 = train_test_split(X, y, test_size=0.2,random_state = 7)
        self.clf1 = CatBoostClassifier(iterations=2000,learning_rate=0.1,l2_leaf_reg=1, depth=11,metric_period = 200, loss_function='MultiClass', eval_metric='MultiClass', random_seed=99, od_type='Iter', od_wait=300)
        self.clf1.fit(x1,y1,verbose=True,eval_set=(x2,y2),use_best_model=True)
        self.clf2 = CatBoostClassifier(iterations=2000,learning_rate=0.1, depth=10,metric_period = 200, loss_function='MultiClass', eval_metric='MultiClass', random_seed=99, od_type='Iter', od_wait=300)
        self.clf2.fit(x1,y1,verbose=True,eval_set=(x2,y2),use_best_model=True)
    def predict(self, X):
        return self.clf.predict(X)

    def predict_proba(self, X):
        return np.array(np.array([[(a[0]+b[0])/2,(a[1]+b[1])/2,(a[2]+b[2])/2] for a,b in zip(self.clf2.predict_proba(X),self.clf1.predict_proba(X))]))

#################################################
#################################################
def find_index(l,v):
    res = []
    for i, j in enumerate(l):
        if(j == v):
            res.append(i)
    return res

def get_label(y_score,p1=0.5,p2=0.5):
    return [np.argmax(y) for y in y_score]
    if(p1 == 0):
        return [np.argmax(y) for y in y_score]
    else:
        res = []
        for i in range(len(y_score)):
            #res.append(np.argmax(y_score[i]))
            if(y_score[i][1]>p1):
                res.append(1)
            elif(y_score[i][0]>p2):
                res.append(0)
            else:
                res.append(2)
        return res

def mesure_class(y_pred,y_true,j):
    TP = 0
    FP = 0
    FN = 0
    for i in range(len(y_pred)):
        if(y_pred[i] == j):
            if(y_true[i] == j):
                TP += 1
            else:
                FP += 1
    for i in range(len(y_true)):
        if(y_true[i] == j):
            if(y_pred[i] == j):
                pass
            else:
                FN += 1
    return TP,FP,FN

def score(tp,fp,fn,epsilon=10**-5):
    beta = 2
    p = tp/(tp+fp+epsilon)
    r = tp/(tp+fn+epsilon)
    beta_squared = beta ** 2
    f = (beta_squared + 1) * (p * r) / (beta_squared * p + r+epsilon)

    def_context.Report("|| precison: "+str(p)+"|| recall: "+str(r)+"|| fbeta: "+str(f))
    def_context.Report('--------------------------------------------------')

def mesure(y_score,y_test,p1=0.5,p2=0.5):
    y = get_label(y_score,p1,p2)
    TP1,FP1,FN1 = mesure_class(y,y_test,0)
    TP2,FP2,FN2 = mesure_class(y,y_test,1)
    TP3,FP3,FN3 = mesure_class(y,y_test,2)
    def_context.Report("pour la classe 0")
    score(TP1,FP1,FN1)
    def_context.Report("pour la classe 1")
    score(TP2,FP2,FN2)
    def_context.Report("pour la classe 2")
    score(TP3,FP3,FN3)

def mismatch(y_score,y_test,p1=0.5,p2=0.5):
    y = get_label(y_score,p1,p2)
    FP = 0
    FF = 0
    for i in range(len(y)):
        if(y[i]==1):
            if(y_test[i]==2):
                FP += 1
            else:
                pass
        if(y[i]==2):
            if(y_test[i]==1):
                FF += 1
            else:
                pass
        else:
            pass
    def_context.Report("fausses publicités")
    def_context.Report(FP)
    def_context.Report("fausses fins")
    def_context.Report(FF)
    return 0
def acc(y_score,y_test,p1=0.5,p2=0.5):
    res = 0
    y = get_label(y_score,p1,p2)
    for i in range(len(y)):
        if(y[i] == y_test[i]):
            res+=1
        else:
            pass
    def_context.Report("accuracy: "+str(res/len(y)))

def ROC_curve(y_score,Y_test):
    y_test = label_binarize(Y_test, classes=[0, 1, 2])
    n_classes = y_test.shape[1]
    fpr = dict()
    tpr = dict()
    roc_auc = dict()
    for i in range(n_classes):
        fpr[i], tpr[i], _ = roc_curve(y_test[:, i], y_score[:, i])
        roc_auc[i] = auc(fpr[i], tpr[i])

    # Compute micro-average ROC curve and ROC area
    fpr["micro"], tpr["micro"], _ = roc_curve(y_test.ravel(), y_score.ravel())
    roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])


    # Compute macro-average ROC curve and ROC area

    # First aggregate all false positive rates
    all_fpr = np.unique(np.concatenate([fpr[i] for i in range(n_classes)]))

    # Then interpolate all ROC curves at this points
    mean_tpr = np.zeros_like(all_fpr)
    for i in range(n_classes):
        mean_tpr += interp(all_fpr, fpr[i], tpr[i])

    # Finally average it and compute AUC
    mean_tpr /= n_classes

    fpr["macro"] = all_fpr
    tpr["macro"] = mean_tpr
    roc_auc["macro"] = auc(fpr["macro"], tpr["macro"])

    # Plot all ROC curves
    plt.figure(figsize=(15,15))
    plt.plot(fpr["micro"], tpr["micro"],
             label='micro-average ROC curve (area = {0:0.2f})'
                   ''.format(roc_auc["micro"]),
             color='deeppink', linestyle=':', linewidth=4)

    plt.plot(fpr["macro"], tpr["macro"],
             label='macro-average ROC curve (area = {0:0.2f})'
                   ''.format(roc_auc["macro"]),
             color='navy', linestyle=':', linewidth=4)

    colors = cycle(['aqua', 'darkorange', 'cornflowerblue'])
    for i, color in zip(range(n_classes), colors):
        plt.plot(fpr[i], tpr[i], color=color, lw=2,
                 label='ROC curve of class {0} (area = {1:0.2f})'
                 ''.format(i, roc_auc[i]))



    plt.plot([0, 1], [0, 1], 'k--', lw=2)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Some extension of Receiver operating characteristic to multi-class')
    plt.legend(loc="lower right")
    plt.savefig('auc.png')
    plt.show()


def save_model_cat(model):
    model.clf1.save_model(PATH_OUT+"model_PTV/catboostmodel1")
    model.clf2.save_model(PATH_OUT+"model_PTV/catboostmodel2")

def save_model_xgb(model):
    pickle.dump(model.clf1, open(PATH_OUT+"model_PTV/XGB1.pickle.dat", "wb"))
    pickle.dump(model.clf2, open(PATH_OUT+"model_PTV/XGB2.pickle.dat", "wb"))

def save_model(clf,name):
    pickle.dump(clf, open(PATH_OUT+"model_PTV/"+name+".pickle.dat", "wb"))

def use_logisticreg(l1,l2,l3,l4,l5,Y_train):
    X_valid = pd.concat([pd.DataFrame(l1),pd.DataFrame(l2),pd.DataFrame(l3),pd.DataFrame(l4),pd.DataFrame(l5)],axis = 1).values
    np.random.seed(7)
    logistic = linear_model.LogisticRegression(C=1,class_weight='balanced',penalty='l2')
    logistic.fit(X_valid, Y_train)
    pickle.dump(logistic, open(PATH_OUT+'model_PTV/logistic_regression.sav', 'wb'))
    return logistic

#################################################################
def main(argv):
    global PATH_IN,PATH_SCRIPT,PATH_OUT
    PATH_IN,PATH_SCRIPT,PATH_OUT = def_context.get_path()
    PATH_OUT = get_temp_path()
    if not os.path.exists(PATH_OUT+'model_PTV/'):
        os.makedirs(PATH_OUT+'model_PTV/')
    if(len(argv) == 0):
        argv = ['all']
    if(argv[0] == 'test'):
        Y_test = pd.read_csv('results.csv').values
        y_pred = pd.read_csv('y_pred.csv')
        y_pred2 = pd.read_csv('y_pred2.csv')
        y_pred3 = pd.read_csv('y_pred2.csv')
        y_pred4 = pd.read_csv('y_pred4.csv')
        y_pred5 = pd.read_csv('y_pred5.csv')

        logreg = use_logisticreg(y_pred,y_pred2,y_pred3,y_pred4,y_pred5,Y_test)
        res = pd.concat([y_pred,y_pred2,y_pred3,y_pred4,y_pred5],axis=1).values
        res = logreg.predict_proba(res)
        for p1 in [0]:
            for p2 in [0]:
                def_context.Report('################### '+str(p1)+' ### '+str(p2)+'###################')
                def_context.Report('############XGB##############')
                mesure(y_pred.values,Y_test,p1,p2)
                mismatch(y_pred.values,Y_test,p1,p2)
                acc(y_pred.values,Y_test,p1,p2)
                def_context.Report('############CatBoost##############')
                mesure(y_pred2.values,Y_test,p1,p2)
                mismatch(y_pred2.values,Y_test,p1,p2)
                acc(y_pred2.values,Y_test,p1,p2)
                def_context.Report('############GradientBoostingClassifier##############')
                mesure(y_pred4.values,Y_test,p1,p2)
                mismatch(y_pred4.values,Y_test,p1,p2)
                acc(y_pred4.values,Y_test,p1,p2)
                def_context.Report('############RandomForestClassifier##############')
                mesure(y_pred5.values,Y_test,p1,p2)
                mismatch(y_pred5.values,Y_test,p1,p2)
                acc(y_pred5.values,Y_test,p1,p2)
                def_context.Report('############Stack##############')
                mesure(res,Y_test,p1,p2)
                mismatch(res,Y_test,p1,p2)
                acc(res,Y_test,p1,p2)

    elif(len(argv) == 1):
        X,Y = load_all(argv[0])
        X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2)
        X_train = X_train.replace([np.inf, -np.inf], np.nan)
        X_train = X_train.fillna(0)
        X_test = X_test.replace([np.inf, -np.inf], np.nan)
        X_test = X_test.fillna(0)
        Y_test  = [Y[0] for Y in Y_test.values]
        ##########################################
        np.random.seed(42)
        clf = Classifier()
        clf.fit(X_train,Y_train)
        y_pred = clf.predict_proba(X_test)
        clf2 = Classifier2()
        clf2.fit(X_train,Y_train)
        y_pred2 = clf2.predict_proba(X_test)

        dtree_model = DecisionTreeClassifier(max_depth = 10).fit(X_train,Y_train)
        y_pred3 = dtree_model.predict_proba(X_test)

        tpot = GradientBoostingClassifier(learning_rate=0.05, max_depth=10, max_features=0.75, min_samples_leaf=7, min_samples_split=16, n_estimators=500, subsample=0.9)
        tpot.fit(X_train,Y_train)
        def_context.Report(tpot.score(X_test, Y_test))
        y_pred4 = tpot.predict_proba(X_test)

        RF_model = RandomForestClassifier(max_depth = 10).fit(X_train,Y_train)
        y_pred5 = RF_model.predict_proba(X_test)

        y_p = clf.predict_proba(X_train)
        y_p2 = clf2.predict_proba(X_train)
        y_p3 = dtree_model.predict_proba(X_train)
        y_p4 = tpot.predict_proba(X_train)
        y_p5 = RF_model.predict_proba(X_train)

        logreg = use_logisticreg(y_p,y_p2,y_p3,y_p4,y_p5,Y_train)

        ##########################################
        save_model_xgb(clf)
        save_model_cat(clf2)
        save_model(dtree_model,"DT")
        save_model(RF_model,"RF")
        pickle.dump(tpot, open(PATH_OUT+"model_PTV/GradientBoostingClassifier.pickle.dat", "wb"))
        pickle.dump(RF_model, open(PATH_OUT+"model_PTV/RandomForestClassifier.pickle.dat", "wb"))
        X = pd.concat([pd.DataFrame(y_pred),pd.DataFrame(y_pred2),pd.DataFrame(y_pred3),pd.DataFrame(y_pred4),pd.DataFrame(y_pred5)],axis = 1).values
        res = logreg.predict_proba(X)
        for p1,p2 in zip([0],[0]):
            def_context.Report('############XGB##############')
            mesure(y_pred,Y_test,p1,p2)
            mismatch(y_pred,Y_test,p1,p2)
            acc(y_pred,Y_test,p1,p2)
            def_context.Report('############CatBoost##############')
            mesure(y_pred2,Y_test,p1,p2)
            mismatch(y_pred2,Y_test,p1,p2)
            acc(y_pred2,Y_test,p1,p2)
            def_context.Report('############DecisionTreeClassifier##############')
            mesure(y_pred3,Y_test,p1,p2)
            mismatch(y_pred3,Y_test,p1,p2)
            acc(y_pred3,Y_test,p1,p2)
            def_context.Report('############GradientBoostingClassifier##############')
            mesure(y_pred4,Y_test,p1,p2)
            mismatch(y_pred4,Y_test,p1,p2)
            acc(y_pred4,Y_test,p1,p2)
            def_context.Report('############RandomForestClassifier##############')
            mesure(y_pred5,Y_test,p1,p2)
            mismatch(y_pred5,Y_test,p1,p2)
            acc(y_pred5,Y_test,p1,p2)
            def_context.Report('############Stack##############')
            mesure(res,Y_test,p1,p2)
            mismatch(res,Y_test,p1,p2)
            acc(res,Y_test,p1,p2)

        #ROC_curve(y_pred,Y_test)
        #ROC_curve(y_pred2,Y_test)
        pd.DataFrame(Y_test).to_csv('results.csv',index=False)
        pd.DataFrame(y_pred).to_csv('y_pred.csv',index=False)
        pd.DataFrame(y_pred2).to_csv('y_pred2.csv',index=False)
        pd.DataFrame(y_pred3).to_csv('y_pred3.csv',index=False)
        pd.DataFrame(y_pred4).to_csv('y_pred4.csv',index=False)
        pd.DataFrame(y_pred5).to_csv('y_pred5.csv',index=False)


    return ("process achevé sans erreures")

if __name__ == "__main__":
    # execute only if run as a script
    main(sys.argv[1:])
