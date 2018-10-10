#-*- coding: utf-8 -*-
#################################################
#created the 20/04/2018 12:57 by Alexis Blanchet#
#################################################

'''
pre-processing des datas
fabrication des Times Series
selection et agrégation de flux
prend en entrée:
    un fichier .csv
    un département(si 0 tous les départements, 99 étranger)
    une chaîne
    une catégorie csp (si 0 agrégation de toutes les csp)

renvoie:
    un fichier csv contenant la colonne d'index correspondant
    aux minutes de la journée et la colonne des valeurs d'audition
'''

'''
Améliorations possibles:
mettre en entrée des listes pour permettre l'agrégation de plusieures
csp ou départements (pas intéressant pour le moment)
'''
#################################################
###########        Imports      #################
#################################################
import warnings
warnings.filterwarnings('ignore')
import random
import os
import sys
from os import listdir
import pandas as pd
from os.path import isfile, join
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
    return datas['PathtoDatasIn'].values[0],datas['PathtoScripts'].values[0],datas['PathtoTempDatas'].values[0]

def Report(error):
    with open(LOG,'a+') as file:
        file.write(str(error)+' \n')
        print(str(error))

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
        Report("Mauvais numéro/nom de chaîne")
        return 0,0

#################################################
########### main with options ###################
#################################################


def main(argv):
    global PATH_IN,PATH_SCRIPT,PATH_OUT
    PATH_IN,PATH_SCRIPT,PATH_OUT = get_path()
    t = time.time()
    if(len(argv) == 3):
        try:
            Report("cleaningRTSfiles called for dep %s chaine %s csp %s" %(str(argv[0]),str(argv[1]),str(argv[2])))
            files = os.listdir(PATH_IN+"RTS/")
            for file in files:
                if(file.split('.')[-1]=='csv'):
                    print(file)
                    departement = argv[0]
                    csp = argv[2]
                    os.system("python "+PATH_SCRIPT+"cleaningRTSfiles.py "+str(file)+" " +str(departement)+' "'+str(argv[1])+'" '+str(csp))
                else:
                    pass
            Report("temps d'éxecution du script pour {0} fichiers: {1}({2} per file)".format(str(len(files)),str(time.time()-t),str((time.time()-t)/len(files))))
            return ("process fini")
        except Exception as e:
            Report("Failed to clean all files: {1}".format(str(argv[0]), str(e)))
    else:
        try:
            for chaine in argv[2].split('-'):
                chaine,name = get_tuple(chaine)
                if name == 0:
                    Report("erreur dans la récupération des informations de la chaîne")
                    return 0
                file = argv[0]
                departement = argv[1]
                csp = argv[3]

                f = file.split(".")[0]
                f = f.split('_')

                df = pd.read_csv(PATH_IN+'/RTS/'+str(file),sep=';')
                df = df.loc[df["IDCST"]== int(chaine)]
                df = df.loc[df["DPT"]== int(departement)]
                if(int(csp) == 0):
                    df = df.sum(axis = 0)
                else:
                    df = df.loc[df["IDCIBLE"] == int(csp)]
                df = df.drop(["IDCST","DPT","IDCIBLE"])
                c= list(str(chaine))
                chaine = "0"*(4-len(c))+chaine
                df.to_csv(PATH_OUT+"RTS/"+str(f[0])+"_"+str(f[1])+"_"+str(departement)+"_"+str(chaine)+"_"+str(csp)+"_cleandata.csv",header=['values'],index= False)
        except Exception as e:
            Report("Failed to process {0}: {1}".format(str(argv[0]), str(e)))



if __name__ == "__main__":
    # execute only if run as a script
    main(sys.argv[1:])
