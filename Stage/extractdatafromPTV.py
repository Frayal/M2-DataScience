#################################################
#created the 25/04/2018 15:41 by Alexis Blanchet#
#################################################
#-*- coding: utf-8 -*-
'''
NOTE A MOI MÊME: Qu'es que c'est chiant un fichier XML
PLEASE PLEASE supprimez moi ce truc!
'''

'''
Améliorations possibles:
ATTENTION dans le nom du fichier figure l'équivalence chaine-identification
il faut donc les deux pour lire le bon fichier mais une simple lecture des
noms peut donner l'équivalence des tables
'''
import warnings
warnings.filterwarnings('ignore')
#################################################
###########        Imports      #################
#################################################
import sys
import numpy as np
import pandas as pd
import xml.etree.ElementTree as ET
import xmltodict, json
from collections import OrderedDict
import os
import time
import datetime
#################################################
########### Global variables ####################
#################################################
PATH_IN ='/home/alexis/Bureau/finalproject/DatasIn/'
PATH_SCRIPT = '/home/alexis/Bureau/finalproject/scripts/'
PATH_OUT = '/home/alexis/Bureau/finalproject/Datas/'
LOG = "log.txt"


#################################################
########### Important functions #################
#################################################

def get_path():
    datas = pd.read_csv('path.csv')
    return datas['PathtoDatasIn'].values[0],datas['PathtoScripts'].values[0],datas['PathtoTempDatas'].values[0]

def Report(error):
    with open(LOG,'a+') as file:
        file.write(str(error)+' \n')
        print(str(error))


def processPTV(file):
    irrelevant = ['LESPERSONNES','DOCUMENTAIRE','LESTEXTES','PHOTO','MAGAZINEXY','SERIE','NATIONALITE','FILMTELEFILM','CS','DOLBY51','PT2','VM','RATIO','PT1','DIRECT','INEDIT_CRYPTE','REDIF','HD','CSA','STM','DOLBY','@CLEEMI','CLAIR','DERDIF'
                ,'OCCULTATIONMOBILE','PREMDIF','GENRE','RESUME','TEMPSFORT','VOST','INEDIT_EN_CLAIR','@DATEMODIF','STEREO','NOUVEAUTE','PT3','DIFFERE','SURTITRE','SOUSTITREDIF','TITREEMISSION','DATE']
    tree = ET.parse(file)
    root = tree.getroot()


    #Report(root.tag)
    l1 = [child for child in root]
    #Report(l1)
    l2 = [[child for child in root] for root in l1 ]
    #Report(l2[1])

    data = pd.DataFrame()
    for i,elem in enumerate(l2[1]):
        xmlstr = ET.tostring(elem)
        o = xmltodict.parse(xmlstr)
        df2 = pd.DataFrame.from_dict(o["DIFFUSION"],orient = 'index').rename(columns={0: i})
        data = [data,df2.T]
        data = pd.concat(data)


    temp_data = pd.DataFrame()
    for i in range(data.shape[0]):
        df = data['ATTRIBUTS']
        l = (df[[i]].values)[0]
        df2 = pd.DataFrame.from_dict(l,orient = 'index')
        df2 = df2.rename(columns = {0: i})
        #Report(df2.head())
        temp_data = [temp_data,df2]
        temp_data = pd.concat(temp_data,axis = 1)
    temp_data.head()
    data = data.join(temp_data.T)
    temp_data = pd.DataFrame()
    names = data.columns.values
    n = list(set(names) - set(irrelevant))
    data = data[n]
    for i in range(data.shape[0]):
        df = data['EMISSION']
        l = (df[[i]].values)[0]
        df2 = pd.DataFrame.from_dict(l,orient = 'index')
        df2 = df2.rename(columns = {0: i})
        #Report(df2.head())
        temp_data = [temp_data,df2]
        temp_data = pd.concat(temp_data,axis = 1)
        temp_data = temp_data.drop(['@DATEMODIF'])

    temp_data.head()

    data = data.join(temp_data.T)
    temp_data = pd.DataFrame()
    data = data.drop(['ATTRIBUTS','EMISSION'],axis = 1)

    for i in range(data.shape[0]):
        df = data['FORMAT']
        l = (df[[i]].values)[0]
        df2 = pd.DataFrame.from_dict(l,orient = 'index')
        df2 = df2.rename(columns = {0: i})
        #Report(df2.head())
        temp_data = [temp_data,df2]
        temp_data = pd.concat(temp_data,axis = 1)
    temp_data = temp_data.drop(['#text'])
    temp_data = temp_data.T.rename(columns={"@CLE": "CLE-FORMAT"})

    temp_data.head()
    data = data.join(temp_data)
    temp_data = pd.DataFrame()
    data = data.drop(['FORMAT'],axis=1)

    for i in range(data.shape[0]):
        df = data['GENRE']
        l = (df[[i]].values)[0]
        df2 = pd.DataFrame.from_dict(l,orient = 'index')
        df2 = df2.rename(columns = {0: i})
        #Report(df2.head())
        temp_data = [temp_data,df2]
        temp_data = pd.concat(temp_data,axis = 1)
        #temp_data = temp_data.drop(['@GENRESIMPLE'])

    temp_data = temp_data.T.rename(columns={"@CLE": "CLE-GENRE","#text":"description programme","@GENRESIMPLE":"GENRESIMPLE"})
    temp_data.head()
    data = data.join(temp_data)
    temp_data = pd.DataFrame()



    data['DATE'] = data['DATEHEURE'].apply(lambda x: x.split('T')[0])
    data['HEURE'] = data['DATEHEURE'].apply(lambda x: x.split('T')[1])
    data['debut'] = data['HEURE'].apply(lambda x: int(x.split(':')[0])*60 + int(x.split(':')[1]))
    names = data.columns.values
    n = list(set(names) - set(irrelevant))
    data = data[n]
    return data

def get_tuple(argv):
    df = pd.read_csv('Equivalence.csv',sep = ';')
    try:
        argv = int(argv)
        key = 'id_unique'
    except Exception:
        key = 'nom_chaine'
    try:
        return str(df[df[key] == argv]['id_unique'].values[0]),str(df[df[key] == argv]['nom_chaine'].values[0])
    except Exception as e:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        Report("Failed to process {0} at line {2} in {3}: {1}".format(str(argv), str(e),sys.exc_info()[-1].tb_lineno,fname))
        Report("Mauvais numéro/nom de chaîne")
        return 0,0

def get_previous_date(date):
    date = date.split('-')
    date = datetime.datetime(int(date[0]),int(date[1]),int(date[2])) + datetime.timedelta(days=-1)
    y = date.year
    m = date.month
    d = date.day
    m = "0"*(2-len(list(str(int(m)))))+str(int(m))
    d = "0"*(2-len(list(str(int(d)))))+str(int(d))
    past_date = str(y)+'-'+str(m)+'-'+str(d)
    return past_date



#################################################
########### main with options ###################
#################################################


def main(argv):
    global PATH_IN,PATH_SCRIPT,PATH_OUT
    PATH_IN,PATH_SCRIPT,PATH_OUT = get_path()
    Report("extractdatafromPTV called for chaine %s" %(str(argv)))
    t = time.time()
    chaine,name = get_tuple(argv)
    if name == 0:
        Report("erreur dans la récupération des informations de la chaîne")
        return 0
    c = list(chaine)
    n = 4-len(c)
    c = ['0']*n
    a = "".join(c)
    chaine = a+chaine
    files = os.listdir(PATH_IN+'PTV/')
    for file in files:
        try:
            date = file.split('_')[0]
            Report(date)
            f = str(PATH_IN+'PTV/'+file+'/IPTV_'+chaine+'_'+str(date)+'_'+name+'.xml')
            data = processPTV(f)
            nd = get_previous_date(date)
            try:
                g = str(PATH_IN+'PTV/'+nd+'_programmeTV/IPTV_'+chaine+'_'+str(nd)+'_'+name+'.xml')
                past_day = processPTV(g)
                entire_day = pd.concat([past_day,data],axis= 0).reset_index()
                l = entire_day.index[(pd.to_numeric(entire_day['debut']) <= 3*60) & (pd.to_numeric(entire_day['debut'])+pd.to_numeric(entire_day['DUREE']) > 3*60)].tolist()
                if(len(l) ==0):
                    Report("pas d'overlap entre journées. Ne devrait pas arriver.")
                if(l[0] == 0):
                    begining = l[1]
                else:
                    begining = l[0]
                end = l[-1]
                data = entire_day.iloc[begining:end+1].drop('index',axis=1)
                data.to_csv(PATH_OUT+'PTV/'+'IPTV_'+str(chaine)+'_'+str(date)+'_'+name+'.csv',index=False)
            except Exception as e:
                Report("Début du mois: le programme du jour précédent n'est pas disponible: %s"%(e))

        except Exception as e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            Report("Failed to process {0} at line {2} in {3}: {1}".format(str(file), str(e),sys.exc_info()[-1].tb_lineno,fname))


    Report("Processed %s files in %ss(%ss per file)" %(len(files),time.time()-t,(time.time()-t)/len(files)))
    return ("process achevé sans erreures")



if __name__ == "__main__":
    # execute only if run as a script
    main(sys.argv[1])
