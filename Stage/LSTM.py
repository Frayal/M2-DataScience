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
from sklearn.preprocessing import MinMaxScaler,StandardScaler
from sklearn.base import BaseEstimator
from sklearn.model_selection import train_test_split
from keras.layers import Dense, Dropout, Activation
from keras.layers import Dense, LSTM, SimpleRNN
from ast import literal_eval

#################################################
########### Global variables ####################
#################################################
THRESHOLD = 0.5
PATH_IN = '/home/alexis/Bureau/finalproject/DatasIn/RTS/'
PATH_SCRIPT = '/home/alexis/Bureau/finalproject/scripts/'
PATH_OUT = '/home/alexis/Bureau/finalproject/Datas/'
LOG = "log.txt"
######################################################
from keras import backend as K
from keras.models import Sequential
from keras.layers import Dense, Dropout, Activation


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


    scaler = MinMaxScaler(feature_range=(0, 1))#StandardScaler()
    X = scaler.fit_transform(X.values)
    X = np.reshape(X,(X.shape[0],1,X.shape[1]))
    return  X,y.values.reshape(-1, 1),t



def precision(y_true, y_pred, threshold_shift=0.5-THRESHOLD):
    beta = 1

    # just in case of hipster activation at the final layer
    y_pred = K.clip(y_pred, 0, 1)

    # shifting the prediction threshold from .5 if needed
    y_pred_bin = K.round(y_pred + threshold_shift)

    tp = K.sum(K.round(y_true * y_pred_bin)) + K.epsilon()
    fp = K.sum(K.round(K.clip(y_pred_bin - y_true, 0, 1)))
    fn = K.sum(K.round(K.clip(y_true - y_pred, 0, 1)))

    precision = tp / (tp + fp)
    return precision


def recall(y_true, y_pred, threshold_shift=0.5-THRESHOLD):
    beta = 1

    # just in case of hipster activation at the final layer
    y_pred = K.clip(y_pred, 0, 1)

    # shifting the prediction threshold from .5 if needed
    y_pred_bin = K.round(y_pred + threshold_shift)

    tp = K.sum(K.round(y_true * y_pred_bin)) + K.epsilon()
    fp = K.sum(K.round(K.clip(y_pred_bin - y_true, 0, 1)))
    fn = K.sum(K.round(K.clip(y_true - y_pred_bin, 0, 1)))

    recall = tp / (tp + fn)
    return recall


def fbeta(y_true, y_pred, threshold_shift=0.5-THRESHOLD):
    beta = 2

    # just in case of hipster activation at the final layer
    y_pred = K.clip(y_pred, 0, 1)

    # shifting the prediction threshold from .5 if needed
    y_pred_bin = K.round(y_pred + threshold_shift)

    tp = K.sum(K.round(y_true * y_pred_bin)) + K.epsilon()
    fp = K.sum(K.round(K.clip(y_pred_bin - y_true, 0, 1)))
    fn = K.sum(K.round(K.clip(y_true - y_pred, 0, 1)))

    precision = tp / (tp + fp)
    recall = tp / (tp + fn)

    beta_squared = beta ** 2
    return (beta_squared + 1) * (precision * recall) / (beta_squared * precision + recall)

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

def model_fit(X,y,X_t,y_t):
    class_weight={
    1: 1/(np.sum(y) / len(y)),
    0:1}
    # create and fit the LSTM network
    np.random.seed(7)
    model = Sequential()
    model.add(SimpleRNN(500,input_shape=(1, 29), return_sequences=True))
    model.add(LSTM(32, return_sequences=True))  # returns a sequence of vectors of dimension 32
    model.add(LSTM(32, return_sequences=True))  # returns a sequence of vectors of dimension 32
    model.add(LSTM(32))  # return a single vector of dimension 32
    model.add(Dense(1))
    model.add(Activation('sigmoid'))
    model.compile(loss='binary_crossentropy', optimizer='adamax',metrics=[fbeta,precision,recall])
    histories = Histories()
    model.fit(X, y,validation_data=(X_t,y_t), epochs=300, batch_size=50, verbose=2,class_weight = class_weight,callbacks=[histories])
    return model,histories

def find_index(l,v):
    res = []
    for i, j in enumerate(l):
        if(j == v):
            res.append(i)
    return res


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
    Report('------------------<LSTM>---------------------------')
    Report("|| precison: "+str(p)+"|| recall: "+str(r)+"|| fbeta: "+str(f))


    tp,fp,fn = mesure(pred,y)
    beta = 2
    p = tp/(tp+fp)
    r = tp/(tp+fn)
    beta_squared = beta ** 2
    f = (beta_squared + 1) * (p * r) / (beta_squared * p + r)


    Report("|| precison: "+str(p)+"|| recall: "+str(r)+"|| fbeta: "+str(f))
    Report('--------------------------------------------------')


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

from keras.callbacks import Callback
class Histories(Callback):
    def on_train_begin(self, logs={}):
        self.losses = []
        self.val_losses = []
        self.fbeta = []
        self.val_fbeta = []
        return
    def on_train_end(self, logs={}):
        return

    def on_epoch_begin(self, epoch, logs={}):
        return

    def on_epoch_end(self, epoch, logs={}):
        self.losses.append(logs.get('loss'))
        self.val_losses.append(logs.get('val_loss'))
        self.fbeta.append(logs.get('fbeta'))
        self.val_fbeta.append(logs.get('val_fbeta'))
        return
    def on_batch_begin(self, batch, logs={}):
        return
    def on_batch_end(self, batch, logs={}):
        return


#################################################
########### main with options ###################
#################################################


def main(argv):
    global PATH_IN, PATH_SCRIPT, PATH_OUT
    PATH_IN, PATH_SCRIPT, PATH_OUT = get_path()
    if (len(argv) == 2):
        from keras.models import model_from_json
        fileX = argv[0]
        X_minmax, X_meanvar, t = load_(fileX)
        json_file = open("model/LSTM.json", 'r')
        loaded_model_json = json_file.read()
        json_file.close()
        LSTM = model_from_json(loaded_model_json)
        # load weights into new model
        LSTM.compile(loss='binary_crossentropy', optimizer='adamax', metrics=['accuracy'])
        LSTM.load_weights("model/LSTM.h5")
       # Report("Loaded LSTM from disk")
        res = pd.DataFrame(LSTM.predict_proba(np.reshape(X_minmax,(X_minmax.shape[0],1,X_minmax.shape[1]))))
        res.to_csv(str(fileX.split('.')[0])+'_temp_LSTM.csv',index=False)
        return res
    else:
        if(len(argv)==0):
            argv = [0.35]
        THRESHOLD = float(argv[0])
        ##### get files names ###
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

        model,hist = model_fit(X_train,Y_train,X_valid,Y_valid)
        import matplotlib.pyplot as plt
        plt.plot(hist.losses,'b',hist.val_losses,'r')
        plt.savefig('plot3.png')
        plt.close()
        plt.plot(hist.fbeta,'b',hist.val_fbeta,'r')
        plt.savefig('plot4.png')
        pred = model.predict_proba(X_test)
        pred = np.clip(pred,0,1)
        testPredict = list([1 if i[0]>THRESHOLD else 0 for i in pred])


        # plot results
        plot_res(t,testPredict,Y_test)

        pred_valid = model.predict_proba(X_valid)
        res_valid = pd.DataFrame(pred_valid)
        res_valid.to_csv('LSTM_valid.csv',index=False)

        res = pd.DataFrame(pred)
        res.to_csv('LSTM.csv',index=False)

        model_json = model.to_json()
        with open("model/LSTM.json", "w") as json_file:
            json_file.write(model_json)
            # serialize weights to HDF5
        model.save_weights("model/LSTM.h5")
        return res

if __name__ == "__main__":
    # execute only if run as a script
    main(sys.argv[1:])
