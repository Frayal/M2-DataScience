#-*- coding: utf-8 -*-
#################################################
#created the 27/07/2018 16:52 by Alexis Blanchet#
#################################################

'''
Objectif:
réunir et unifer les façcon d'extraire le context et de la traiter
Fichier contenat uniquement des focntions évitant une redondance et une disparité
 du traitement. N'existait pas dans la version antérieur du projet. Peut être sujet
 a modifications si détection de problèmes (WIP)

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
import pandas as pd
import numpy as np
import datetime
import pickle
from catboost import CatBoostClassifier
import xgboost as xgb
from tqdm import tqdm
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier
import time
from sklearn import linear_model
import re
#################################################
########### Global variables ####################
#################################################
PATH_IN = '/home/alexis/Bureau/finalproject/Datas/'
PATH_SCRIPT = '/home/alexis/Bureau/finalproject/scripts/'
PATH_OUT = '/home/alexis/Bureau/finalproject/Datas/'
LOG = "log.txt"
THRESHOLD = 0.4
#################################################
########### Important functions #################
#################################################

def get_tuple(argv):
    '''
    renvoie le nom de la chaîne et son numéro (écrit avec 4 chiffres ex: TF1 -> '0192')
    '''
    # Fichier conteant les paires
    df = pd.read_csv('Equivalence.csv',sep = ';')
    # On essaye de mettre la chaine sous la forme d'un numero pour savoir si l'on passe
    # en entrée le nom de la chaîne ou son numéro
    try:
        argv = int(argv)
        key = 'id_unique'
    except Exception:
        key = 'nom_chaine'
    # On va vérifier que la chaîne que l'on cherche existe vraiment
    try:
        number,name =  str(df[df[key] == argv]['id_unique'].values[0]),str(df[df[key] == argv]['nom_chaine'].values[0])
        number = "0"*(4-len(list(number)))+number
        return number,name
    except Exception as e:
        # eeeeeeet on oublie pas de dire tout au patron si on se rate et qu'on se prend un mur
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        Report("Failed to process {0} at line {2} in {3}: {1}".format(str(argv), str(e),sys.exc_info()[-1].tb_lineno,fname))
        Report("Mauvais numéro/nom de chaîne")
        return 0,0

def get_path():
    # Mouais
    datas = pd.read_csv('path.csv')
    return datas['PathtoTempDatas'].values[0],datas['PathtoScripts'].values[0],datas['PathtoDatasOut'].values[0]

def Report(error):
    #Obvious
    with open(LOG,'a+') as file:
        file.write(str(error)+' \n')
        print(str(error))

def part(i):
    #Useless subsubfunction....but it's still used
    if(i<12*60):
        return 1
    elif(i<20*60):
        return 2
    elif(i<27*60):
        return 3
    else:
        return 0


def load_file(date,c):
    #allllllller paf! one fait confiance a personne.
    # Retrouvons le tuple de la chaîne.
    numero,nom = get_tuple(c)
    try:
        JOINDATE = "".join(date.split('-'))
        c = list(numero)
        numero = "0"*(4-len(c))+numero
        PTV = pd.read_csv(PATH_IN+'PTV/IPTV_'+numero+'_'+date+'_'+nom+'.csv')
        PTV['fin'] = PTV['debut']+PTV['DUREE']
        #Report('IPTV_'+numero+'_'+date+'_'+nom+'.csv','pred_proba_'+str(JOINDATE)+'_'+numero+'.csv')
        proba = pd.read_csv(PATH_IN+'RTS/pred_proba_'+str(JOINDATE)+'_'+numero+'.csv').values
        return PTV,proba
    except:
        Report("Fichier Non existant: "+PATH_IN+'PTV/IPTV_'+numero+'_'+date+'_'+nom+'.csv')
        return [],[]

def init_newPTV(PTV,chaine):
    #Les deux premières options sont useless mais bon...C'est ca qu'on aime!
    if(chaine == 'TF1'):
        #Initialisation du NewPTV
        newPTV = pd.DataFrame()
        newPTV['minute'] = [180]
        newPTV['TITRE'] = 'Programmes de la nuit'
        newPTV['Change Point'] = 'Non'
        newPTV['pourcentage vu'] = 0
        newPTV['Évenement'] = 'Début de Détection'
        return newPTV
    elif(chaine == 'M6'):
        newPTV = pd.DataFrame()
        newPTV['minute'] = [180]
        newPTV['TITRE'] = 'M6 Music'
        newPTV['Change Point'] = 'Non'
        newPTV['pourcentage vu'] = 0
        newPTV['Évenement'] = 'Début de Détection'
        return newPTV
    else:
        newPTV = pd.DataFrame()
        newPTV['minute'] = [180]
        newPTV['TITRE'] = PTV[(PTV['debut'] < 3*60) & (PTV['debut']+PTV['DUREE'] > 3*60) ]['TITRE']
        newPTV['Change Point'] = 'Non'
        newPTV['pourcentage vu'] = 0
        newPTV['Évenement'] = 'Début de Détection'
        return newPTV
def init_history(chaine,PTV,lastend,currentduree):
    #Oui je sais! tout est a 0, pas besoin de faire une fonction pour ca...et bah si! moi j'en ai besoin!
    h = pd.DataFrame()
    h['minute'] = [179]
    h['partie de la journée'] = 'nuit'
    h['Change Point'] = 0
    h['pourcentage'] = 0
    h['partie du programme'] = 0
    h['programme'] = "programme de nuit"
    h['duree'] = 0
    h['nombre de pub potentiel'] =  0
    h['lastCP'] =  lastend
    h['lastPub'] = 500
    h['lastend'] =  lastend
    h['currentduree'] = currentduree
    h['Pubinhour'] =  0
    h['probability of CP'] = 0
    h['nb de pubs encore possible'] = 0
    h["chaine"]= chaine
    h['CLE-FORMAT'] = 0
    h['CLE-GENRE'] = 0
    h['day'] = 0
    h['part'] = 0
        #h['per'] = 1
    return h



def find_position(seen_percentage):
    if(seen_percentage<=0.25):
        return "début de programme"
    if(0.25<seen_percentage<=0.75):
        return "milieu de programme"
    if(0.75<seen_percentage<=0.95):
        return "début de fin de programme"
    if(0.95<seen_percentage<=1.05):
        return "fin de programme"
    if(seen_percentage>1.05):
        return "en dehors du programme"

def find_partofday(i):
    part = i/60
    if(3<=part<=7):
        return("fin de nuit")
    elif(7<part<=9):
        return("début de matinée")
    elif(9<part<=12):
        return("matinée")
    elif(12<part<=14):
        return("midi")
    elif(14<part<=17):
        return("après-midi")
    elif(17<part<=20):
        return("fin d'après-midi")
    elif(20<part<=24):
        return("prime time")
    else:
        return('nuit')

def find_threshold(cp,threshold):
    """
    Work In Progress (en gros je l'ai pas testé et j'ai aucune idée si c'est une bonne mauvaise idée)
    """
    global THRESHOLD
    num = 0
    if(threshold>=0.60):
        THRESHOLD = 0.6
        return 0.6
    if(threshold<=0.35):
        THRESHOLD = 0.35
        return 0.35
    for c in cp:
        if(c[0]>threshold):
            num+=1
        else:
            pass
    Report(num)
    if(60<=num<=100):
        THRESHOLD = threshold
        return threshold
    elif(num<60):
        return find_threshold(cp,threshold-0.01)
    elif(100<num):
        return find_threshold(cp,threshold+0.01)


def find_ifChangePoint(i,cp):
    if(cp[i-183]>THRESHOLD):
        return 1
    else:
        return 0

def categorize_duree(duree):
    if(duree<10):
        return("très court")
    elif(10<=duree<30):
        return("court")
    elif(30<=duree<60):
        return("moyen")
    elif(60<=duree<=100):
        return("long")
    elif(100<=duree<=180):
        return("très long")
    else:
        return('super long')

def categorize_type(description):
    mots = description.split(" ")
    if(mots[0] == "Série"):
        return "Série"
    if(mots[0] == "Téléfilm"):
        return "film"
    if(description == 'Magazine jeunesse'):
        return 'dessins animés'
    if(mots[0] == "Magazine"):
        return "magazine"
    if(mots[0]=="Feuilleton"):
        return "Feuilleton"
    else:
        return description

def categorize_pub(name,debut,duree,titre,chaine,PTV,index_PTV):
    '''
    Résumé de cette fonction: CANCER CANCER CANCER
    Non plus sérieusement, il y a ici sujet à discussion:
    la vérité importe peu, seul la vérité détectable est intéressante.
    Si on sait qu'il y a deux publicité mais qu'on n'en détecte qu'une seule à chaque fois,
    ne tentons pas le diable, soyons gentils, resposables, malhonnêtes et dions qu'il n'y en
    a qu'une...
    '''
    if("Journal" in titre.split(' ')):
        return 0
    elif(chaine == 'TF1' and (debut>20*60 or debut<180)):
        if(titre in ['Nos chers voisins'] and debut>20*60):
            return 0
        elif(titre in ['Nos chers voisins','Reportages découverte']):
            return 2
        elif(titre in ['Automoto']):
            return 1
        elif(titre in ['Journal','Téléshopping','Tirage du Loto']):
            return 0
        elif(titre in ['50mn Inside','Téléfoot']):
            if(duree in ['moyen']):
                return 1
            elif(duree in ["long"]):
                return 2
            elif(duree in ["court"]):
                return 1
            else:
                return 4
        elif(name in["Météo","Magazine","magazine"] and duree in ["court","très court"]):
            return 0
        elif(name in ['dessins animés'] and duree != "super long"):
            return 3
        elif(name in ['Jeu']):
            return 1
        elif(name in ['Feuilleton','film','Drame','Thriller']):
            return 2
        elif(name == 'Série' and (debut>21*60) and titre == PTV['TITRE'].loc[index_PTV-1]):
            return 2
        elif(name == 'Série' and (debut>21*60 or debut<180)):
            return 1
        elif(name == 'Série' and 180<debut<12*60):
            return 1
        elif(name == 'Série'and 12*60<debut<21*60):
            return 0
        elif(name == 'Série'):
            return 2
        elif(name in ['Téléréalité'] and debut < 20.5*60):
            return 2
        elif(name in ['Téléréalité'] and debut > 20.5*60):
            return 3
        elif(duree in ["très court","court"]):
            return 0
        else:
            return 4
    if(chaine == 'TF1' and 180<=debut<=20*60):
        if(titre in ['Nos chers voisins']):
            return 2
        elif(titre in ['Journal']):
            return 0
        elif(titre == '50mn Inside'):
            if(duree in ["moyen","long"]):
                return 2
            if(duree in ["court"]):
                return 1
            else:
                return 4
        elif(name in["Météo","Magazine","magazine"] and duree in ["court","très court"]):
            return 0
        elif(name in ['dessins animés'] and duree != "super long"):
            return 3
        elif(name in ['Jeu']):
            return 1
        elif(name in ['Feuilleton','film','Drame','Thriller']):
            return 2
        elif(name == 'Série' and 180<debut<12*60):
            return 1
        elif(name == 'Série'and 12*60<debut<21*60):
            return 0
        elif(name == 'Série'):
            return 2
        elif(name in ['Téléréalité'] and debut < 20.5*60):
            return 2
        elif(name in ['Téléréalité'] and debut > 20.5*60):
            return 3
        elif(duree in ["très court","court"]):
            return 0
        else:
            return 4

    if(chaine == 'M6' and (debut>20*60 or debut <180)):
        if(titre in ['Nos chers voisins','Les reines du shopping']):
            return 2
        elif(titre in ['Une superstar pour Noël']):
            return 3
        elif(titre in ['En famille']):
            if(duree in ["court","moyen",'très court','long']):
                return 0
            else:
                return 2
        elif(titre in ["Chasseurs d'appart'"]):
            if(duree in ["très long","super long"] ):
                return 12
            else:
                return 1
        elif(titre in ['Les aventures de Tintin','Absolument stars','Martine','Les Sisters','M6 boutique','Scènes de ménages','M6 Music']):
            return 0
        elif(titre in ['66 minutes : grand format']):
            return 1
        elif(titre =='Turbo' and debut<11*60):
            return 2
        elif(titre in ['66 minutes']):
            return 4
        elif(name in ['Dessin animé','dessins animés']):
            return 0
        elif(name in ['Journal']):
            return 0
        elif(titre == '50mn Inside'):
            if(duree in ["moyen","long"]):
                return 2
            else:
                return 4
        elif(name in ['magazine'] and duree not in ["court","très court"]):
            return 2
        elif(name in["Météo","Magazine","magazine"] and duree in ["court","très court"]):
            return 0
        elif(name in ['dessins animés'] and duree != "super long"):
            return 3
        elif(name in ['Jeu'] and duree in ['court','très court','moyen']):
            return 2
        elif(name in ['Feuilleton','film','Drame','Thriller']):
            return 2
        elif(name == 'Série' and (debut<13*60 or 25*60>debut>22*60+20)):
            return 1
        elif(name == 'Série'):
            return 2
        elif(name in ['Téléréalité'] and debut < 20.5*60):
            if( duree in ["court","moyen"]):
                return 1
            else:
                return 2
        elif(duree in ["très court","court"]):
            return 0
        else:
            return 10
    if(chaine == 'M6' and 180<=debut<=20*60):
        if(titre in ['Nos chers voisins','Les reines du shopping']):
            return 2
        elif(titre in ['Une superstar pour Noël']):
            return 3
        elif(titre in ['En famille']):
            if(duree in ["court","moyen",'très court','long']):
                return 0
            else:
                return 2
        elif(titre in ["Chasseurs d'appart'"]):
            if(duree in ["très long","super long"] ):
                return 12
            else:
                return 1
        elif(titre in ['Les aventures de Tintin','Absolument stars','Martine','Les Sisters','M6 boutique','Scènes de ménages']):
            return 0
        elif(titre in ['66 minutes : grand format']):
            return 1
        elif(titre =='Turbo' and debut<11*60):
            return 1
        elif(titre in ['66 minutes','M6 Music']):
            return 4
        elif(name in ['Dessin animé','dessins animés']):
            return 0
        elif(name in ['Journal']):
            return 0
        elif(titre == '50mn Inside'):
            if(duree in ["moyen","long"]):
                return 2
            else:
                return 4
        elif(name in ['magazine'] and duree not in ["court","très court"]):
            return 2
        elif(name in["Météo","Magazine","magazine"] and duree in ["court","très court"]):
            return 0
        elif(name in ['dessins animés'] and duree != "super long"):
            return 3
        elif(name in ['Jeu'] and duree in ['court','très court','moyen']):
            return 2
        elif(name in ['Feuilleton','film','Drame','Thriller']):
            return 2
        elif(name == 'Série' and 180<debut<13*60):
            return 1
        elif(name == 'Série'):
            return 2
        elif(name in ['Téléréalité'] and debut < 20.5*60):
            if( duree in ["court","moyen"]):
                return 1
            else:
                return 2
        elif(name in ['Téléréalité'] and debut > 20.5*60):
            return 3
        elif(duree in ["très court","court"]):
            return 0
        else:
            return 10
    elif(chaine in ['France 2','France 3'] and 6*60>debut>20*60 ):
        return 0
    elif(chaine in ['France 2','France 3']):
        return 1.5
    else:
        return(categorize_pub(name,debut,duree,titre,'TF1',PTV,index_PTV))

def encoding_partoftheday(x):
    if(x == 'nuit'):
        return 1
    if(x == 'fin de nuit'):
        return 2
    if(x == 'début de matinée'):
        return 3
    if(x == 'matinée'):
        return 4
    if(x == 'midi'):
        return 5
    if(x == 'après-midi'):
        return 6
    if(x == "fin d'après-midi"):
        return 7
    if(x == 'prime time'):
        return 8
    elif(type(x) == str):
        return 0
    else:
        return x

def encoding_partofprogramme(x):
    if(x == 'début de programme'):
        return 1
    if(x == 'milieu de programme'):
        return 2
    if(x == 'début de fin de programme'):
        return 3
    if(x == 'fin de programme'):
        return 4
    if(x == 'en dehors du programme'):
        return 5
    elif(type(x) == str):
        return 0
    else:
        return x


def encoding_duree(x):
    if(x == 'très court'):
        return 1
    if(x == 'court'):
        return 2
    if(x == 'moyen'):
        return 3
    if(x == 'long'):
        return 4
    if(x == 'très long'):
        return 5
    if(x == 'super long'):
        return 6
    elif(type(x) == str):
        return 0
    else:
        return x

def encoding_chaine(x):
    number,name = get_tuple(x)
    return(int(number))
def encoding_per(x):
    if(x>=1):
        return 1
    if(x<1):
        return 0

def process(df):
    df['Heure'] = df['minute'].apply(lambda x: str(int(x/60))+':'+str(x%60))
    df['partie de la journée'] = df['partie de la journée'].apply(lambda x: encoding_partoftheday(x))
    df['partie du programme'] = df['partie du programme'].apply(lambda x: encoding_partofprogramme(x))
    df["duree"] = df['duree'].apply(lambda x: encoding_duree(x))
    df['chaine'] = df['chaine'].apply(lambda x: encoding_chaine(x))
    if('per' in df.columns.values):
        df['per'] = df['per'].apply(lambda x: encoding_per(x))
    if('Heure' in df.columns.values):
        df['Time-h'] = df['Heure'].apply(lambda x: int((x.split(':'))[0]))
        df['Time-m'] = df['Heure'].apply(lambda x: int((x.split(':'))[1]))
        df = df.drop(['Heure'],axis = 1)
    if('programme' in df.columns.values):
        df = df.drop(['programme'],axis = 1)
    return df


def categorize_programme(programme,chaine,PTV,index_PTV):
    p = []
    p.append(categorize_type(programme['description programme']))
    p.append(categorize_duree(programme['DUREE']))
    p.append(categorize_pub(p[0],programme['debut'],p[-1],programme['TITRE'],chaine,PTV,index_PTV))
    return p

def get_important_points(chaine,PTV,index_PTV):
    num,chaine = get_tuple(chaine)
    index_PTV -= 1
    if(index_PTV<0):
        index_PTV = PTV.shape[0]
    if(chaine == 'TF1'):
        importantpts = [[13*60,"Journal"],[20*60,"Journal"],[(int(PTV['HEURE'].loc[index_PTV-1].split(':')[0])+24)*60+int(PTV['HEURE'].loc[index_PTV-1].split(':')[1]),PTV['TITRE'].loc[index_PTV-1]]]
    elif(chaine == 'M6'):
        importantpts = [[12*60+45,"Le 12.45"],[19*60+45,"Le 19.45"],[(int(PTV['HEURE'].loc[index_PTV-1].split(':')[0])+24)*60+int(PTV['HEURE'].loc[index_PTV-1].split(':')[1]),PTV['TITRE'].loc[index_PTV-1]]]
    elif(chaine == 'France 2'):
        importantpts = [[13*60,"Journal 13h00"],[20*60,"Journal 20h00"],[(int(PTV['HEURE'].loc[index_PTV-1].split(':')[0])+24)*60+int(PTV['HEURE'].loc[index_PTV-1].split(':')[1]),PTV['TITRE'].loc[index_PTV-1]]]
    elif(chaine == 'France 3'):
        importantpts = [[12*60+25,"12/13 : Journal national"],[19*60+30,"19/20 : Journal national"],[(int(PTV['HEURE'].loc[index_PTV-1].split(':')[0])+24)*60+int(PTV['HEURE'].loc[index_PTV-1].split(':')[1]),PTV['TITRE'].loc[index_PTV-1]]]

    else:
        Report('Pas de Points de contrôle pour cette chaîne. Merci de les reseigner dans le fichier def_context.py ligne 425')
    return importantpts

def get_help(chaine,PTV):
    num,chaine = get_tuple(chaine)
    help = []
    if(chaine == 'TF1'):
        if(PTV[PTV['TITRE'] =='TFou']['DUREE'].iloc[0]<=120):
            help = [8*60+25,8*60+30,10*60+55,12*60+50,19*60+50,11*60+52]
        else:
            help = [10*60+55,12*60+50,19*60+50,11*60+52]
    elif(chaine == 'M6'):
        help = [6*60,6*60+50,8*60+57,10*60,10*60+20]
    j = 0
    for t in PTV['GENRESIMPLE'].values:
        heure_in_text = re.search(r'\d{2}.\d{2}',t)
        if(t == 'Journal'):
            help.append(PTV['debut'].iloc[j])
        else:
            if heure_in_text is None:
                pass
            else:
                h = list(heure_in_text)
                h = int(h[0]+h[1])*60+int(h[3]+h[4])
                if(h == PTV['debut'].iloc[j]):
                    help.append(h)
        j+=1
    return help

def get_context(i,programme,lastCP,lastPub,lastend,currentduree,planifiedend,Pubinhour,probas,nbpub,chaine,per,PTV,index_PTV,date):
    #we create a list with different notes to understand the context
    # minute of the point and its situation in the day
    context = [i]
    context.append(find_partofday(i))
    # Is the Point a Change Point
    context.append(find_ifChangePoint(i,probas))
    # Where is the Point in the programme:
    seen_percentage = (i-lastend)/currentduree
    context.append(seen_percentage)
    context.append(find_position(seen_percentage))
    # which type of programme we are watching
    p = categorize_programme(programme,chaine,PTV,index_PTV)
    for j in range(len(p)):#3
        context.append(p[j])
    context.append(lastCP)
    context.append(lastPub)
    context.append(lastend)
    context.append(currentduree)
    context.append(Pubinhour)
    context.append(probas[i-183][0])
    context.append(context[7]-nbpub)
    context.append(chaine)
    context.append(programme['CLE-FORMAT'])
    context.append(programme['CLE-GENRE'])
    day = datetime.datetime(int(date.split('-')[0]), int(date.split('-')[1]), int(date.split('-')[2])).weekday()
    context.append(day)
    context.append(part(i))
    #context.append(per)
    return context

def get_temp_path():
    datas = pd.read_csv('path.csv')
    return datas['temp_path'].values[0]


def load_models(path = get_temp_path()):
    #Work of Art
    p = path.split('/')
    p[-2] = 'T'+str(int(''.join(list(p[-2][1:])))-1)
    path = '/'.join(p)
    XGB = []
    XGB.append(pickle.load(open(path+"model_PTV/XGB1.pickle.dat", "rb")))
    XGB.append(pickle.load(open(path+"model_PTV/XGB2.pickle.dat", "rb")))

    CatBoost = []
    CatBoost.append(CatBoostClassifier().load_model(fname=path+"model_PTV/catboostmodel1"))
    CatBoost.append(CatBoostClassifier().load_model(fname=path+"model_PTV/catboostmodel2"))

    rf = pickle.load(open(path+"model_PTV/RF.pickle.dat", "rb"))
    dt = pickle.load(open(path+"model_PTV/DT.pickle.dat", "rb"))
    gb = pickle.load(open(path+"model_PTV/GradientBoostingClassifier.pickle.dat", "rb"))
    logistic = pickle.load(open(path+'model_PTV/logistic_regression.sav', "rb"))
    return XGB,CatBoost,rf,dt,gb,logistic
