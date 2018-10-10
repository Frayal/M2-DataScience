# DON T RUN THIS THING IT WILL BREAK MY MODEL 
# I REPEAT DON T RUN IT
# ACHTUNG CAREFULL DO NOT RUN THIS THING


import pandas as pd
import numpy as np
import sys
from random import shuffle


PATH_IN ='../DatasIn/'
PATH_SCRIPT = '../scripts/'
PATH_OUT = '../Datas/'
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

def make_names(date):
    d = str((date.split('-'))[1])+str((date.split('-'))[0])
    X = PATH_IN+'train/sfrdaily_2018'+d+'_0_192_0_cleandata.csv'
    y = PATH_IN+'train/label-'+date+'.csv'
    return X,y


def make_dataframe(dates):
    train = len(dates)-2
    valid=1
    test=1
    shuffle(dates)
    data = pd.DataFrame()
    fileX_test = []
    fileY_test = []
    fileX_train = []
    fileY_train = []
    fileX_valid = []
    fileY_valid = []

    j = 0
    for i in range(train):
        X,y = make_names(dates[j])
        fileX_train.append(X)
        fileY_train.append(y)
        j+=1
    for i in range(valid):
        X,y = make_names(dates[j])
        fileX_valid.append(X)
        fileY_valid.append(y)
        j+=1
    for i in range(test):
        X,y = make_names(dates[j])
        fileX_test.append(X)
        fileY_test.append(y)
        j+=1

    data['fileX_train'] = [fileX_train]
    data['fileY_train'] = [fileY_train]

    data['fileX_valid'] = [fileX_valid]
    data['fileY_valid'] = [fileY_valid]

    data['fileX_test'] = [fileX_test]
    data['fileY_test'] = [fileY_test]

    data.to_csv('files.csv',index=False)

    return data






def main(argv):
    global PATH_IN,PATH_SCRIPT,PATH_OUT
    PATH_IN,PATH_SCRIPT,PATH_OUT = get_path()
    Dates = ['30-04','07-05','09-05','18-05','23-05','28-05']
    make_dataframe(Dates)

    return "processe has end"

if __name__ == "__main__":
    # execute only if run as a script
    main(sys.argv[1:])
