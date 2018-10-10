#################################################
#created the 30/05/2018 11:56 by Alexis Blanchet#
#################################################
#-*- coding: utf-8 -*-
'''
Script pour l'établissement d'un premier historique
Peux appeller PTVM6 ou PTVF1
Dans le cas ou il s'agit d'une autre chaîne, appele d'un des deux fichiers,
mais avec des initialisations particulières.
(WIP)
Étude d'impact du bruitage volontaire de Dataset en lien avec l'étude de
l'inertie du projet. Recherche de Points d'équilibre stables à partir d'un
dataset bruité.



Prends en entrée un mois
Considère l'ensemble des fichiers correspondant à ce mois
Pour éviter la redondance, cherche la liste des jours dans le dossier des PTV
renvoie Deux fichiers sous deux formats à chaque fois
- Le nouveau PTV en html et en csv
- l'historique du context et la décision pour chaque Chnage Point en html et en csv
Écriture de ces fichiers dans le dossier Data (ne sert qu'à l'entraînement)

'''

'''
Améliorations possibles:
Réécriture nécéssaire si l'on veux utiliser cet algo en dernier recours
(quand les deux autres n'ont pas trouvé de bonne réponse)
'''
import warnings
warnings.filterwarnings('ignore')
#################################################
###########        Imports      #################
#################################################
import sys
import os
import pandas as pd
import numpy as np
import subprocess
import time
import def_context
#################################################
########### Global variables ####################
#################################################
PATH_IN = '/home/alexis/Bureau/finalproject/Datas/'
PATH_SCRIPT = '/home/alexis/Bureau/finalproject/scripts/'
PATH_OUT = '/home/alexis/Bureau/finalproject/Datas/'
LOG = "log.txt"

#################################################
########### Important functions #################
#################################################

def get_path():
    datas = pd.read_csv('path.csv')
    return datas['PathtoTempDatas'].values[0],datas['PathtoScripts'].values[0],datas['PathtoTempDatas'].values[0]



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
        def_context.Report("Failed to process {0} at line {2} in {3}: {1}".format(str(argv), str(e),sys.exc_info()[-1].tb_lineno,fname))
        Report("Mauvais numéro/nom de chaîne")
        return 0,0
#################################################
########### main with options ###################
#################################################



def main(argv):
    t = time.time()
    if(len(argv)==0):
        Report("Merci de renseigner une année et un mois (ex: 2017-12)")
    EPSILON = 1e-15
    err = 0
    m = 0
    err_TF1 = 0
    m_TF1 = 0
    err_M6 = 0
    m_M6 = 0
    err_F2 = 0
    m_F2 = 0
    err_F3 = 0
    m_F3 = 0

    err_type_1 = 0
    err_type_2 = 0
    err_type_3 = 0
    try:
        df = pd.read_csv('scores.csv')
    except Exception as e:
        df= pd.DataFrame()
        df['score TF1'] = [0]
        df['score M6'] = 0
        df['score France 2'] = 0
        df['score France 3'] = 0
        df['score Total'] = 0
        df['score sur la matinée'] = 0
        df["score sur l'après midi"] = 0
        df['score sur la soirée'] = 0
        df['part de relecture'] = 0
        df['temps de calcul'] = 0
        df['mois'] = '55-55'
        df.to_csv('scores.csv',index=False)

    files = os.listdir(PATH_IN+'PTV/')
    for file in files:
        def_context.Report('-------------------------------------')
        f = ((file.split('.'))[0].split('_'))[2]
        c = ((file.split('.'))[0].split('_'))[-1]
        if(f=='2017-12-20' or (f in ['2017-12-09','2017-12-06','2018-02-22'] and c=='TF1') or (f in ['2018-02-22'] and c=='M6') or (f.split('-')[0] != str(argv[0].split('-')[0])) or f.split('-')[1] != argv[0].split('-')[1]):
            def_context.Report(f)
        elif(c ==''):
            pass
        else:
            def_context.Report(c)
            if(c in ['M6','TF1']):
                chaine = c
            else:
                chaine = 'TF1'
            number,name = get_tuple(c)
            if(len(list(number))<4):
                number = "0"+number
            def_context.Report('Using PTV%s for %s'%(chaine,f))
            l = os.system('python '+PATH_SCRIPT+'PTV'+str(chaine)+'.py '+str(f)+' '+str(number))
            if(l/256 == 4):
                pass
            else:
                l = l/256
                err += int(l/100)+int((l%100)/10)+int((l%10))
                err_type_1 += int(l/100)
                err_type_2 += int((l%100)/10)
                err_type_3 += int((l%10))
                m+=3
                if(c == 'M6'):
                    err_M6 += int(l/100)+int((l%100)/10)+int((l%10))
                    m_M6 += 3
                if(c == 'TF1'):
                    err_TF1 += int(l/100)+int((l%100)/10)+int((l%10))
                    m_TF1 += 3
                if(c == 'France 3'):
                    err_F2 += int(l/100)+int((l%100)/10)+int((l%10))
                    m_F2 += 3
                if(c == 'France 3'):
                    err_F3 += int(l/100)+int((l%100)/10)+int((l%10))
                    m_F3 += 3


        def_context.Report(err)
    def_context.Report(m)
    if(m == 0):
        def_context("aucun fichier n'a été traité. Merci de vérifier la date et les données d'entrée.")
    def_context.Report("score Total:"+str(1-(err/(m+EPSILON))))
    def_context.Report("score TF1:"+str(1-(err_TF1/(m_TF1+EPSILON))))
    def_context.Report("score M6:"+str(1-(err_M6/(m_M6+EPSILON))))
    def_context.Report("score France 2:"+str(1-(err_F2/(m_F2+EPSILON))))
    def_context.Report("score France 3:"+str(1-(err_F3/(m_F3+EPSILON))))
    def_context.Report("score sur la matinée:"+str(1-((err_type_1*3)/(m+EPSILON))))
    def_context.Report("score sur l'après midi:"+str(1-((err_type_2*3)/(m+EPSILON))))
    def_context.Report("score sur la soirée:"+str(1-((err_type_3*3)/(m+EPSILON))))
    def_context.Report("temps de calcul:"+str(time.time()-t))
    try:
        df = pd.read_csv('scores.csv')
        df.loc[df.shape[0]] = [1-(err_TF1/(m_TF1+EPSILON)),1-(err_M6/(m_M6+EPSILON)),1-(err_F2/(m_F2+EPSILON)),1-(err_F3/(m_F3+EPSILON)),1-(err/(m+EPSILON)),1-((err_type_1*3)/(m+EPSILON)),1-((err_type_2*3)/(m+EPSILON)),1-((err_type_3*3)/(m+EPSILON)),0.5*m/(m+EPSILON),time.time()-t,argv[0]]
        df.to_csv('scores.csv',index=False)
    except Exception as e:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        def_context.Report("Failed to process {0} at line {2} in {3}: {1}".format('', str(e),sys.exc_info()[-1].tb_lineno,fname))
        Report("fichier non conforme ou non existant: %s" %(e))

if __name__ == "__main__":
    # execute only if run as a script
    main(sys.argv[1:])
