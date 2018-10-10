#################################################
#created the 19/04/2018 12:57 by Alexis Blanchet#
#################################################
#-*- coding: utf-8 -*-
'''
exploration des données des chaînes
recherche de distribution
visualisation des données
récupération des événements notables dans l'historique
'''

'''
Améliorations possibles:
- meilleur crible des données pour trouver des sous catégories donnant
des fonnées plus utilisable (et utiliser ces spécificités)
- une fois l'historique des pubs trouvées, enlever le thresholding et le remplacer
par un xgboost afin d'optimiser les poids de chaque détécteur
- écrire le main pour que le fichier execute l'ensemble des actions automatiquement
(ce serait pas mal)
'''
import warnings
warnings.filterwarnings('ignore')
#################################################
###########        Imports      #################
#################################################
import numpy as np
import pandas as pd
from pandas import DataFrame
from pandas import concat
import random
import os
import sys
import scipy.stats
import time

import random

#################################################
########### Global variables ####################
#################################################
NOISE_LEVEL = 0
THRESHOLD = 2e7

PATH_IN = '/home/alexis/Bureau/Project/Datas/RTS/'
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
def load_timeserie(file,PATH):
    '''
    prend en entrée le nom d'un fichier se trouvant dans le dossier
    des times Series
    renvoie en sortie un DataFrame contenant les valeurs d'audition
    ATTENTION: le processus d'écriture peut se faire en même temps
    que la lecture ==> Tester l'option os.O_NONBLOCK qui doit régler
    le problème.
    '''
    temp = file.split('.')
    if(temp[-1] == 'txt'):
        data_file = open(PATH+str(file), 'r',os.O_NONBLOCK)
        data_file = data_file.read()
        datas = data_file.split('\n')
        data = pd.DataFrame([float(x) for x in datas[:-1]])
        return data
    if(temp[-1] == 'csv'):
        data = pd.read_csv(PATH+str(file))
        data = data.replace([np.inf, -np.inf], np.nan)
        data = data.fillna(1)
        return data
    else:
        print("mauvais format: veuillez fournir un .txt ou un .csv")
        return 0

def shift_preprocess(data,name):
    '''
    prend en entrée un dataframe contenant les valeurs d'audition
    renvoie un dataframe de 4 colonnes contenant les valeurs bruitées
    et un shift permettant d'utiliser l'historique (fixé a 3 ici)
    '''
    temps = data['values'].values
    temps = DataFrame(temps)
    dataframe = concat([temps.shift(3), temps.shift(2), temps.shift(1), temps], axis=1)
    dataframe.columns = ['t-3', 't-2', 't-1', 't']
    dataframe['minutes'] = dataframe.index
    dataframe = dataframe.drop(dataframe.index[[0,1,2]])
    '''
    plt.plot(dataframe['t'].values)
    plt.savefig('data/png/'+name+'-0.png')
    '''
    return dataframe



def KFDR(a,b):
    GAMMA = 10
    sigma = 0.5*a[1] + 0.5*b[1]
    return (((4**2/8)*(a[0]-b[0])**2)/(sigma + GAMMA))



def processing(dataframe,name):
    '''
    prend en entrée un dataframe contenant les shifts
    renvoie en sortie un dataframe de features
    permet le plot de certains features intéressants
    a compléter en fonction des features que l'on souhaite
    considérer ou non.

    '''
    SIGMA = 1
    L = 2000000

    THRESHOLD = max(dataframe['t'].values)/100
    dataframe["diff t t-1"]=dataframe.apply(lambda x:max(-THRESHOLD,min(THRESHOLD,x["t"]-x["t-1"])),axis = 1)
    dataframe["diff t t-2"]=dataframe.apply(lambda x:max(-THRESHOLD**2,min(THRESHOLD**2,x["t"]-x["t-2"])),axis = 1)
    dataframe["diff t t-3"]=dataframe.apply(lambda x:max(-THRESHOLD**2,min(THRESHOLD**2,x["t"]-x["t-3"])),axis = 1)
    dataframe["diff t-1 t-2"]=dataframe.apply(lambda x:max(-THRESHOLD,min(THRESHOLD**2,x["t-1"]-x["t-2"])),axis = 1)
    dataframe["diff t-1 t-3"]=dataframe.apply(lambda x:max(-THRESHOLD**2,min(THRESHOLD**2,x["t-1"]-x["t-3"])),axis = 1)
    dataframe["diff t-2 t-3"]=dataframe.apply(lambda x:max(-THRESHOLD,min(THRESHOLD,x["t-2"]-x["t-3"])),axis = 1)

    dataframe["mean"]=(dataframe["t"]+dataframe["t-1"]+dataframe["t-2"]+dataframe["t-3"])/4
    dataframe["distance to mean"] = dataframe.apply(lambda x:max(-THRESHOLD,min(10**7,(x["t"]-x["mean"])/x["t"])),axis = 1)
    dataframe["pente t t-2"] = dataframe.apply(lambda x:max(-THRESHOLD,min(THRESHOLD,x["diff t t-2"]/2)),axis = 1)
    dataframe["pente t t-3"] = dataframe.apply(lambda x:max(-THRESHOLD,min(THRESHOLD,x["diff t t-3"]/3)),axis = 1)
    dataframe["pente t-1 t-3"] = dataframe.apply(lambda x:max(-THRESHOLD,min(THRESHOLD,x["diff t-1 t-3"]/2)),axis = 1)


    dataframe["diff pente 1-2"] = (dataframe["diff t t-1"] - dataframe["diff t-1 t-2"])/dataframe["diff t t-1"]
    dataframe["diff pente 1-3"] = (dataframe["diff t t-1"] - dataframe["diff t-2 t-3"])/dataframe["diff t t-1"]
    dataframe['diff pentes'] = dataframe.apply(lambda x:max(-50,min(50,x["diff pente 1-2"]-x["diff pente 1-3"])),axis = 1)
    dataframe["diff pente 2-3"] = (dataframe["diff t-1 t-2"] - dataframe["diff t-2 t-3"])/dataframe["diff t-1 t-2"]

    dataframe['GP'] = (SIGMA**2)*(np.exp(-((dataframe['t']-dataframe['t-1'])**2)/dataframe['t']**2))
    dataframe['covariance'] = 0.25*((dataframe['t']-dataframe['mean'])**2 + (dataframe['t-1']-dataframe['mean'])**2 +(dataframe['t-2']-dataframe['mean'])**2 + (dataframe['t-3']-dataframe['mean'])**2)
    dataframe['pics'] = dataframe.apply(lambda row: 1 if((row['t-2']<row['t-1']) & (row['t']<row['t-1'])) else -1 if((row['t-2']>row['t-1']) & (row['t']>row['t-1'])) else 0,axis =1)

    h = [0,0]

    for index, row in dataframe.iterrows():
        dataframe.set_value(index, 'KDFR', KFDR([row['mean'],row['covariance']],h))
        h = [row['mean'],row['covariance']]
    dataframe.set_value(3, 'KDFR', 1.5)

    x = dataframe["diff t t-1"]
    m = np.mean(x)
    sd = sum([(y-m)**2 for y in x])/(len(x)-1)
    #print("moyenne: %s   standard deviation: %s" %(str(m),str(sd)))

    d = scipy.stats.norm(m, sd/1000)
    time = [i/((2*60)) for i in range(len(x))]
    prob = [ d.pdf(y) for y in x]
    dist = [ d.cdf(y) for y in x]
    dataframe["probability"] = prob
    dataframe["distribution"] = dist

    dataframe["skewness"]= ((3)**0.5)*(((dataframe['t']-dataframe["mean"])**3)/4+((dataframe['t-1']-dataframe["mean"])**3)/4
                            +((dataframe['t-2']-dataframe["mean"])**3)/4+((dataframe['t-3']-dataframe["mean"])**3)/4)/(sd**3)
    dataframe["kurtosis"]= (((dataframe['t']-dataframe["mean"])**4)/4+((dataframe['t-1']-dataframe["mean"])**4)/4
                            +((dataframe['t-2']-dataframe["mean"])**4)/4+((dataframe['t-3']-dataframe["mean"])**4)/4)/(sd**4)

    return dataframe



def annomalie_detection(df,clf,automatic_threshold = True):
    '''
    prend en entrée un DataFrame
    renvoie en sortie une liste de points ou il y a annomalie
    les threshold sont a déterminer manuellement pour le moment
    on peut cependant essayer de les fixer de manière automatique
    en fonction de l'historique de la chaîne/departement/regroupement(csp)
    '''
    res = []
    irrelevant = ('t-3', 't-2', 't-1', 't')
    allnames = df.columns
    names = list(set(allnames) - set(irrelevant))
    #print(names)
    temp_df = df[names]
    if(automatic_threshold):

        info = df.describe()
        temp_df['signe'] = np.sign(temp_df['pente t t-2']+temp_df['pente t-1 t-3']+temp_df['pente t t-3']+temp_df['diff t-1 t-2'])
        #['diff t-2 t-3', 'distribution', 'skewness', 'diff pente 1-3','covariance', 'KDFR', 'diff pente 1-2', 'pente t-1 t-3', 'diff t-1 t-2','pente t t-2', 'diff pente 2-3', 'distance to mean', 'pente t t-3','GP', 'probability']
        poids = [ 0.01635393,0.03110482,  0.01142569,  0.00263594,  0.02041714,  0.04325364,0.00982603,  0.01245762,  0.02220948,  0.00309737,0.00173761,0.04263638,0.03090817,  0.00556751,  0.01292914,  0.0071749,0.04850749,0.01424759,0.0060423,0.03402776,0.00756308,0.0317241,0.00425174,0.0046664 ]
        for i in range(len(names)-len(poids)):
            poids.append(random.randint(1,4))
        ptot = sum(poids)
        for name,p in zip(names,poids):
            m = info.loc[['mean']][name].values[0]
            sd = info.loc[['std']][name].values[0]
            d = scipy.stats.norm(m, sd)
            temp_df[name] = temp_df[name].apply(lambda x: p*d.cdf(x))


        temp_df['proba'] = temp_df[names].sum(axis = 1)
        temp_df['proba'] = temp_df['proba'].apply(lambda x: abs((x-ptot*0.5)/ptot))
        #temp_df['r'] = [ 1 if (v[1]+l[1])*0.5>0.05    else 0 for v,l in zip(clf[0].predict_proba(df.values),clf[1].predict_proba(df.values))]
        #temp_df['proba'] =  (temp_df['proba']*0.8+ temp_df['r']*0.2)

        temp_df['annomalies'] =  temp_df['proba'].apply(threshold)

        res = temp_df['annomalies']*temp_df['signe']



    else:
        df = df.replace([np.inf, -np.inf], np.nan)
        df = df.fillna(0)
        scaler = MinMaxScaler(feature_range=(0, 1))
        trainX = scaler.fit_transform(df)
        res1 = clf[0].predict_proba(trainX)
        res2 = clf[1].predict_proba(trainX)
        res = [(r1[1]+r2[1])*0.5 for r1,r2 in zip(res1,res2)]
        res = [1 if l>0.17 else 0 for l in res]


    return res

def threshold(x):
    if(x>0.17):
        return 1
    else:
        return 0

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
    if(len(argv) == 0):
        argv = ['test']
    if(len(argv)==1):
        T = time.time()
        if(argv[0] == "train"):
            files = os.listdir(PATH_IN+'train')
            for file in files:
                if((file.split('.'))[-1] == 'csv'):
                    if(((file.split('.'))[0].split('-'))[0] == 'label'):
                        pass
                    else:
                        print(str(file))
                        os.system("python "+PATH_SCRIPT+"processingdata.py "+str(file) +" train")
                else:
                    pass
            Report("Éxecution du scripte processingdata sur le train(%s fichiers) en %s" %(len(files),time.time()-T))
            return ("process achevé sans erreurs")
        else:
            files = os.listdir(PATH_IN+'RTS/')
            for file in files:
                if((file.split('.'))[-1] == 'csv' and file.split('_')[0] != 'pred'):
                    print(str(file))
                    os.system("python "+PATH_SCRIPT+"processingdata.py "+str(file)+" test")
                else:
                    pass
            Report("Éxecution du scripte processingdata sur le test(%s fichiers) en %s" %(len(files),time.time()-T))
            return ("process achevé sans erreurs")
    if(len(argv)==2):
        try:
            if(argv[1] == 'train'):
                df = load_timeserie(argv[0],PATH_IN+'train/')
                if(type(df) == int):
                    return("wrong imput file")
                if(df.shape[1]>5):
                    Report("already treated file: "+str(argv[0]))
                    return 0
                real_data = df['values'][3:]
                df = shift_preprocess(df,argv[0].split('.')[0])
                df = processing(df,argv[0].split('.')[0])
                clf = []
                annomalies = annomalie_detection(df,clf,True)
                df["label"] = annomalies
                df.to_csv(PATH_OUT+'RTS/'+argv[0].split('.')[0]+".csv",index=False)
                return ("process achevé sans erreurs")





            elif(argv[1] == 'test'):

                df = load_timeserie(argv[0],PATH_IN+'RTS/')
                if(type(df) == int):
                    return("wrong imput file")
                if(df.shape[1]>5):
                    Report("already treated file: "+str(argv[0]))
                    return 0
                real_data = df['values'][3:]
                df = shift_preprocess(df,argv[0].split('.')[0])
                df = processing(df,argv[0].split('.')[0])
                clf = []
                annomalies = annomalie_detection(df,clf,True)
                df["label"] = annomalies
                df.to_csv(PATH_OUT+'RTS/'+argv[0].split('.')[0]+".csv",index=False)
                return ("process achevé sans erreurs")
        except Exception as e:
            Report("Failed to process {0} {1}: {2}".format(str(argv[0]),str(argv[1]),str(e)))

    else:
        Report("Wrong imputs arguments, please see the doc and give us something good to wotk with")
        return 0


if __name__ == "__main__":
    # execute only if run as a script
    main(sys.argv[1:])
