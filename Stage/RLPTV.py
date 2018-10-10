#################################################
#created the 30/05/2018 11:56 by Alexis Blanchet#
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
import pandas as pd
import numpy as np
import pickle
from catboost import CatBoostClassifier
import xgboost as xgb
from tqdm import tqdm
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier
import time
from sklearn import linear_model
import datetime
import def_context
#################################################
########### Global variables ####################
#################################################
THRESHOLD = 0.46
PATH_IN = '/home/alexis/Bureau/finalproject/Datas/'
PATH_SCRIPT = '/home/alexis/Bureau/finalproject/scripts/'
PATH_OUT = '/home/alexis/Bureau/finalproject/Datas/'
LOG = "log.txt"
#################################################
########### Important functions #################
#################################################

def make_predictedPTV(actu,PTV,chaine,j,index_PTV,lastCP,lastPub,lastend,currentduree,planifiedend,Pubinhour,proba,nbpub,XGB,CatBoost,rf,dt,gb,logistic,Predictiontimer,per,context,date,importantpts):
    #Initialisation des Variables
    starting_index = index_PTV
    verbose = False
    index_PTV = index_PTV
    ##########################
    a_0 = 200                                      #faute sur le point de contrôle
    a_1 = 5                                       #Coefficient pour le fait de manquer une publicité
    a_2 = 20                                       #Coeffeicient de la pénalité si le journal commence trop tôt
    a_3 = 1                                        #Coeffeicent de péanlité pour ne pas détecter une fin de programme
    a_4 = 20                                       #Paramètre de régularisation -- utilisation des points les plus probables.
    a_5 = 2                                        #Pénalitasion du la valeur absolue de l'éloignement du point prévu de fin.
    error_type_1 = 0                               # erreur sur les programmes de la demi journée
    error_type_2 = 0                               # erreur sur les Points de Contrôle
    error_type_3 = 0                               # Erreur de régularisation
    error_type_4 = 0                               # Erreur du noeud
    ##########################
    Predictiontimer = Predictiontimer
    Pubinhour = Pubinhour
    lastCP = lastCP
    lastPub= lastPub
    lastend = lastend
    currentduree = max(currentduree,1)
    planifiedend = planifiedend
    begin = True
    nbpub = nbpub
    Recall = 1
    wait = 4
    error = 0
    per = per
    context = context
    importantpts = importantpts
    help = def_context.get_help(chaine,PTV)
    if(j> importantpts[0][0]):
        end = importantpts[1][0]
        index_ipts = 1
    else:
        end = importantpts[0][0]
        index_ipts = 0
    #####################################
    if(actu == 0):
        error_type_4+= abs(((context[13]-THRESHOLD)*a_4)*context[2])
    elif(actu == 1):
        lastCP=0
        lastPub = 0
        Pubinhour+=4
        nbpub+=1
    elif(actu == 2):
        lastend = j
        lastCP=0
        index_PTV += 1
        index_PTV = index_PTV%(PTV.shape[0])
        currentduree = PTV['DUREE'].iloc[index_PTV]
        planifiedend = (lastend + currentduree)
        Predictiontimer = 200
        nbpub = 0
        if(context[14]-int(context[14])):
            error_type_4+= abs(context[14])*a_1
        if(context[14]<0):
            error_type_4+= abs(context[14])*50
        error_type_4 += abs((1-context[3])*context[11])*a_5
        per = context[3]
        if(PTV['TITRE'].iloc[index_PTV] == importantpts[index_ipts][1]):
            if(importantpts[index_ipts][0]-j>13):
                error_type_4+= (abs(importantpts[index_ipts][0]-j))*a_2
            else:
                error_type_4 -= abs(importantpts[index_ipts][0]-j)

    #####################################
    for i in range(j+1,end+1):
        ############Update loss################
        a_1 = 5*(0.8**(min(index_PTV-starting_index,1)))      #Coefficient pour le fait de manquer une publicité
        a_2 = 20*(0.8**(min(index_PTV-starting_index,1)))      #Coeffeicient de la pénalité si le journal commence trop tôt
        a_3 = 1 *(0.8**(min(index_PTV-starting_index,1)))      #Coeffeicent de péanlité pour ne pas détecter une fin de programme
        a_4 = 20*(0.8**(min(index_PTV-starting_index,1)))      #Paramètre de régularisation -- utilisation des points les plus probables.
        a_5 = 2 *(0.8**(min(index_PTV-starting_index,1)))     #Pénalitasion du la valeur absolue de l'éloignement du point prévu de fin.
        #Update time of commercials (Reset)
        if(i%60 == 0):
            Pubinhour = 0
        #Update timmers
        lastPub+=1
        lastCP+=1
        if(index_ipts==len(importantpts)):
            index_ipts-=1
        #let's get the context:
        context = def_context.get_context(i,PTV.iloc[index_PTV],lastCP,lastPub,lastend,currentduree,planifiedend,Pubinhour,proba,nbpub,chaine,per,PTV,index_PTV,date)
        #Sur M6 il y a 16 minutes de pub entre deux films!!!!!!!!!!!!.....!!!!!!!....!!.!.!.!.!....!.!...!..!.!.!.!
        if(PTV['GENRESIMPLE'].iloc[index_PTV].split(' ')[0] == PTV['GENRESIMPLE'].iloc[index_PTV-1].split(' ')[0] and PTV['GENRESIMPLE'].iloc[index_PTV].split(' ')[0] == 'Téléfilm'
            and (i-lastend)<2 and Recall > 0 and per<0.97 and chaine == 'M6'):

            lastend = i+5
            lastPub = -25
            Recall -= 0.5
        elif((i-lastend)<2 and Recall > 0 and per<0.95 and chaine == 'M6' and 15*60<i<16*60):

            lastend = i+5
            lastPub = -25
            Recall -= 0.5

        ###### Let's verify that the algo is not doing a crappy predicitions and if this the case, clean his historic #####
        elif(i == importantpts[index_ipts][0]):
            #### we are at an important point, let's now see what the algo has predict
            if(PTV['TITRE'].iloc[index_PTV] == importantpts[index_ipts][1]):
                #Well he doesn't have the programme wrong, that's a good start
                #let's now find out if we are at a logical point of the programme
                if(i-lastend>13):
                    #Wellllll, the programme began way too early...something went wrong before...Let's rest for now, we'll correct the algo later
                    error_type_2+= a_0

                else:
                    # OMG the ALGO IS RIGHT...here is a candy, let's rest a litle just in case...we never know....
                    nbpub = 0




            else:
                #maybe it's the next programme so calme the fuck down!
                if(PTV['TITRE'].iloc[(index_PTV+1)%(PTV.shape[0])] == importantpts[index_ipts][1]):
                    if(planifiedend-i<10):
                        error_type_2 += abs(planifiedend-i)
                        per = context[3]
                        if(context[14]-int(context[14])):
                            error_type_1+= abs(context[14])*a_1
                        error_type_1 += abs((1-context[3])*context[11])*a_5
                        if(context[14]<0):
                            error_type_1+= abs(context[14])*50

                        index_ipts+=1
                    else:

                        error_type_2+= a_0


                else:
                    #well the programme is wrong, and we are not even close to it, let's terminate this thing before it goes completly south. REBOOT The algo, erase the memory, just like in Westworld.
                    #BUT FIRST LET'S VERIFY THAT THERE IS INDEED AN IMPORTANT PROGRAMME THAT DAY...Don't go fuck everything up for no reason
                    l = PTV.index[(PTV['TITRE']==importantpts[index_ipts][1]) & (PTV['debut'] == i)].tolist()
                    if(len(l)>0):

                        error_type_2+= a_0

                    else:
                        print("PAS DE POINT DE CONTRÔLE POUR CETTE PARTIE DE LA JOURNÉE")
                        index_ipts+=1

        elif(context[2]):
            if(lastCP < min(6,currentduree/2)):
                if(lastCP == 1):
                    continue
                else:
                    error_type_3+= abs((context[13]-THRESHOLD)*a_4*(0.95*(4-lastCP)))

            else:
                X = def_context.process(pd.DataFrame([context],index=[0],columns=['minute','partie de la journée','Change Point','pourcentage','partie du programme','programme','duree','nombre de pub potentiel','lastCP','lastPub','lastend','currentduree','Pubinhour','probability of CP','nb de pubs encore possible','chaine','CLE-FORMAT','CLE-GENRE','day','part'])).values #,'per'
                res1 = CatBoost[0].predict_proba(X)
                res2 = CatBoost[1].predict_proba(X)
                res3 = XGB[0].predict(xgb.DMatrix(X), ntree_limit= XGB[0].best_ntree_limit)
                res4 = XGB[1].predict(xgb.DMatrix(X), ntree_limit= XGB[1].best_ntree_limit)
                res5 = rf.predict_proba(X)
                res6 = gb.predict_proba(X)
                res7 = dt.predict_proba(X)
                res = [(res1[0][0]+res2[0][0]+res3[0][0]+res4[0][0]+res5[0][0]+res6[0][0])/6,(res1[0][1]+res2[0][1]+res3[0][1]+res4[0][1]+res5[0][1]+res6[0][1])/6,(res1[0][2]+res2[0][2]+res3[0][2]+res4[0][2]+res5[0][2]+res6[0][2])/6]
                y_pred = [(res1[0][0]+res2[0][0])*0.5,(res1[0][1]+res2[0][1])*0.5,(res1[0][2]+res2[0][2])*0.5]
                y_pred2 = [(res3[0][0]+res4[0][0])*0.5,(res3[0][1]+res4[0][1])*0.5,(res3[0][2]+res4[0][2])*0.5]
                X = pd.concat([pd.DataFrame(y_pred).T,pd.DataFrame(y_pred2).T,pd.DataFrame(res7),pd.DataFrame(res6),pd.DataFrame(res5)],axis = 1)
                X = X.replace([np.inf, -np.inf], np.nan)
                X = X.fillna(1)
                X = X.values
                res = logistic.predict_proba(X)
                cla = np.argmax(res)

                if(cla == 1 and context[14]==0):
                    cla = np.argmax([res[0][0],res[0][2]])*2
                if(cla == 2 and context[3]<0.65 and context[11]>=20):
                    cla = 0
                if(cla == 2 and context[3]<0.9 and context[11]>=180):
                    cla = 0
                if(cla == 2 and context[3]<0.9 and context[0]<9*60 and context[11]>=20):
                    cla = 0
                if(context[3]>1):
                    cla = 2
                if(cla == 1 and context[9]<0):
                    cla = 0
                if(cla == 1):
                    lastCP=0
                    lastPub = 0
                    Pubinhour+=4
                    nbpub+=1
                    wait = 4
                elif(cla == 2):
                    lastend = i
                    lastCP=0
                    index_PTV += 1
                    index_PTV = index_PTV%(PTV.shape[0])
                    currentduree = PTV['DUREE'].iloc[index_PTV]
                    planifiedend = (lastend + currentduree)
                    Predictiontimer = 200
                    nbpub = 0
                    wait = 5
                    per = context[3]
                    if(context[14]-int(context[14])):
                        error_type_1+= abs(context[14])*a_1
                    error_type_1 += abs((1-context[3])*context[11])*a_5
                    if(context[14]<0):
                        error_type_1+= abs(context[14])*50
                    if(PTV['TITRE'].iloc[index_PTV] == importantpts[index_ipts][1]):
                        if(importantpts[index_ipts][0]-i>13):
                            error_type_2+= (abs(importantpts[index_ipts][0]-i))*a_2
                        else:

                            if(chaine == 'M6'):
                                error_type_2 -= importantpts[index_ipts][0]-i
                            else:
                                error_type_2 -= (importantpts[index_ipts][0]-i)*0.5

                else:
                    error_type_3+= abs((context[13]-THRESHOLD)*a_4)



        elif(i in help):
            if(lastCP < min(6,currentduree)):
                error_type_3-= 5
            else:
                X = def_context.process(pd.DataFrame([context],index=[0],columns=['minute','partie de la journée','Change Point','pourcentage','partie du programme','programme','duree','nombre de pub potentiel','lastCP','lastPub','lastend','currentduree','Pubinhour','probability of CP','nb de pubs encore possible','chaine','CLE-FORMAT','CLE-GENRE','day','part'])).values #,'per'
                res1 = CatBoost[0].predict_proba(X)
                res2 = CatBoost[1].predict_proba(X)
                res3 = XGB[0].predict(xgb.DMatrix(X), ntree_limit= XGB[0].best_ntree_limit)
                res4 = XGB[1].predict(xgb.DMatrix(X), ntree_limit= XGB[1].best_ntree_limit)
                res5 = rf.predict_proba(X)
                res6 = gb.predict_proba(X)
                res7 = dt.predict_proba(X)
                res = [(res1[0][0]+res2[0][0]+res3[0][0]+res4[0][0]+res5[0][0]+res6[0][0])/6,(res1[0][1]+res2[0][1]+res3[0][1]+res4[0][1]+res5[0][1]+res6[0][1])/6,(res1[0][2]+res2[0][2]+res3[0][2]+res4[0][2]+res5[0][2]+res6[0][2])/6]
                y_pred = [(res1[0][0]+res2[0][0])*0.5,(res1[0][1]+res2[0][1])*0.5,(res1[0][2]+res2[0][2])*0.5]
                y_pred2 = [(res3[0][0]+res4[0][0])*0.5,(res3[0][1]+res4[0][1])*0.5,(res3[0][2]+res4[0][2])*0.5]
                X = pd.concat([pd.DataFrame(y_pred).T,pd.DataFrame(y_pred2).T,pd.DataFrame(res7),pd.DataFrame(res6),pd.DataFrame(res5)],axis = 1)
                X = X.replace([np.inf, -np.inf], np.nan)
                X = X.fillna(1)
                X = X.values
                res = logistic.predict_proba(X)
                cla = np.argmax(res)

                if(cla == 1 and context[14]==0):
                    cla = np.argmax([res[0][0],res[0][2]])*2
                if(cla == 2 and context[3]<0.65 and context[11]>=20):
                    cla = 0
                if(cla == 2 and context[3]<0.9 and context[11]>=180):
                    cla = 0
                if(context[3]>1):
                    cla = 2
                if(cla == 2 and context[3]<0.9 and context[0]<9*60 and context[11]>=20):
                    cla = 0
                if(cla == 1 and context[9]<0):
                    cla = 0


                if(cla == 1):
                    lastCP=0
                    lastPub = 0
                    Pubinhour+=4
                    nbpub+=1
                    wait = 4
                elif(cla == 2):
                    lastend = i
                    lastCP=0
                    index_PTV += 1
                    index_PTV = index_PTV%(PTV.shape[0])
                    currentduree = PTV['DUREE'].iloc[index_PTV]
                    planifiedend = (lastend + currentduree)
                    Predictiontimer = 200
                    nbpub = 0
                    wait = 5
                    per = context[3]
                    if(context[14]-int(context[14])):
                        error_type_1+= abs(context[14])*a_1
                    error_type_1 += abs((1-context[3])*context[11])*a_5
                    if(context[14]<0):
                        error_type_1+= abs(context[14])*50
                    if(PTV['TITRE'].iloc[index_PTV] == importantpts[index_ipts][1]):
                        if(importantpts[index_ipts][0]-i>13):
                            error_type_2+= (abs(importantpts[index_ipts][0]-i))*a_2
                        else:
                            if(chaine == 'M6'):
                                error_type_2 -= importantpts[index_ipts][0]-i
                            else:
                                error_type_2 -= (importantpts[index_ipts][0]-i)*0.5

                else:
                    error_type_3-=5


        else:
            #labels.append(0)
            #Not a Change Point, we'll just check that nothing is wrong in the PTV at this time
            if(Predictiontimer <= 0):
                l = currentduree
                lastend = i
                lastCP=0
                index_PTV += 1
                index_PTV = index_PTV%(PTV.shape[0])
                currentduree = PTV['DUREE'].iloc[index_PTV]
                planifiedend = (lastend + currentduree)
                Predictiontimer = 200
                nbpub = 0
                if(context[14]-int(context[14])):
                    error_type_1+= abs(context[14])*a_1
                if(context[14]<0):
                    error_type_1+= abs(context[14])*50
                error_type_1 += abs((1-context[3])*context[11])*a_5
                error_type_1 += max(currentduree/5,4)*a_3
                per = context[3]
                if(PTV['TITRE'].iloc[index_PTV] == importantpts[index_ipts][1]):
                    if(importantpts[index_ipts][0]-i>13):
                        error_type_2+= (abs(importantpts[index_ipts][0]-i))*a_2
                    else:
                        if(chaine == 'M6'):
                            error_type_2 -= importantpts[index_ipts][0]-i
                        else:
                            error_type_2 -= (importantpts[index_ipts][0]-i)*0.5

            elif(context[3] == 1):
                #Dépassement autorisé: Modulable en fonction de la position dans la journée si besoin
                # C'est sur ces valeurs que l'on va jouer pour avoir le meilleur PTV possible
                # Plus les valeurs sont grandes, plus on fait confiance a l'algo
                # Il est important de bien découper la journée celon les périodes horaires que l'on qualifie
                # de "sous tension" si plusieurs programmes courts se succédent. Bien évidement une telle analyse sera
                #plus tard fait automatiquement.
                if(chaine == 'TF1'):
                    if(context[5] == 'Journal'):
                        if(i<20*60):
                            Predictiontimer = 10
                        else:
                            Predictiontimer = 0
                    elif(11.5*60<=i<=14*60 or 19.5*60<i<21*60):
                        Predictiontimer = 1
                    elif(context[6] == "très court"):
                        Predictiontimer = 0
                    elif(PTV['TITRE'].iloc[index_PTV] == 'Téléshopping'):
                        Predictiontimer = 5
                    elif(context[6] == "court"):
                        Predictiontimer = 5
                    elif(context[6] == "moyen"):
                        Predictiontimer = 5
                    elif(context[6] == "très long" or context[6] == "long"):
                        Predictiontimer = 15
                    else:
                        Predictiontimer = 5

                elif(chaine =='M6'):
                    #Dépassement autorisé: Modulable en fonction de la position dans la journée si besoin
                    #Dépassement autorisé: Modulable en fonction de la position dans la journée si besoin
                    if(i<9*60+56):
                        Predictiontimer = 2
                    elif(13*60<i<14*60):
                        Predictiontimer = 5
                    elif(PTV['TITRE'].iloc[index_PTV] in ['M6 boutique']):
                        Predictiontimer = 0
                    elif(context[6] == "très court"):
                        Predictiontimer = 0
                    elif(context[6] == "court"):
                        Predictiontimer = 5
                    elif(context[6] == "moyen"):
                        Predictiontimer = 10
                    elif(context[6] == "très long"):
                        Predictiontimer = 10
                    elif(context[6] == 'long'):
                        Predictiontimer = 15
                    else:
                        Predictiontimer = 5
                else:
                    if(context[5] == 'Journal'):
                        if(i<20*60):
                            Predictiontimer = 10
                        else:
                            Predictiontimer = 0
                    elif(context[6] == "très court"):
                        Predictiontimer = 4
                    elif(context[6] == "court"):
                        Predictiontimer = 5
                    elif(context[6] == "moyen"):
                        Predictiontimer = 5
                    elif(context[6] == "très long"):
                        Predictiontimer = 5
                    elif(context[6] == 'long'):
                        Predictiontimer = 15
                    else:
                        Predictiontimer = 5
            elif(context[3]>1):
                Predictiontimer -= 1
            else:
                pass
    error+= error_type_1 + error_type_2 + error_type_3 +error_type_4
    return error,[error_type_1, error_type_2,error_type_3,error_type_4]




def make_predictedPTV2(actu,PTV,chaine,j,index_PTV,lastCP,lastPub,lastend,currentduree,planifiedend,Pubinhour,proba,nbpub,XGB,CatBoost,rf,dt,gb,logistic,Predictiontimer,per,context,date,importantpts):
    #Initialisation des Variables
    starting_index = index_PTV
    verbose = False
    index_PTV = index_PTV
    ##########################
    a_0 = 150                                      #faute sur le point de contrôle
    a_1 = 5                                       #Coefficient pour le fait de manquer une publicité
    a_2 = 20                                       #Coeffeicient de la pénalité si le journal commence trop tôt
    a_3 = 1                                        #Coeffeicent de péanlité pour ne pas détecter une fin de programme
    a_4 = 20                                       #Paramètre de régularisation -- utilisation des points les plus probables.
    a_5 = 2                                        #Pénalitasion du la valeur absolue de l'éloignement du point prévu de fin.
    error_type_1 = 0                               # erreur sur les programmes de la demi journée
    error_type_2 = 0                               # erreur sur les Points de Contrôle
    error_type_3 = 0                               # Erreur de régularisation
    error_type_4 = 0                               # Erreur du noeud
    ##########################
    Predictiontimer = Predictiontimer
    Pubinhour = Pubinhour
    lastCP = lastCP
    lastPub= lastPub
    lastend = lastend
    currentduree = currentduree
    planifiedend = planifiedend
    begin = True
    nbpub = nbpub
    Recall = 1
    wait = 4
    error = 0
    per = per
    context = context
    help = def_context.get_help(chaine,PTV)

    importantpts = importantpts
    index_ipts = 2
    end = importantpts[2][0]
    #####################################
    if(actu == 0):
        error_type_4+= abs(((context[13]-THRESHOLD)*a_4)*context[2])
    elif(actu == 1):
        lastCP=0
        lastPub = 0
        Pubinhour+=4
        nbpub+=1
    elif(actu == 2):
        lastend = j
        lastCP=0
        index_PTV += 1
        index_PTV = index_PTV%(PTV.shape[0])
        currentduree = PTV['DUREE'].iloc[index_PTV]
        planifiedend = (lastend + currentduree)
        Predictiontimer = 200
        nbpub = 0
        if(context[14]-int(context[14])):
            error_type_4+= abs(context[14])*a_1
        if(context[14]<0):
            error+= abs(context[14])*a_4
        error_type_4 += abs((1-context[3])*context[11])*a_5
        per = context[3]
        if(PTV['TITRE'].iloc[index_PTV] == importantpts[index_ipts][1]):
            if(importantpts[index_ipts][0]-j>13):
                error_type_4+= (abs(importantpts[index_ipts][0]-j))*a_2
            else:
                error_type_4 -= abs(importantpts[index_ipts][0]-j)

    #####################################
    for i in range(j+1,min(end+13,1620)):
        ############Update loss################
        a_1 = 5*(0.9**(index_PTV-starting_index))     #Coefficient pour le fait de manquer une publicité
        a_2 = 20*(0.9**(index_PTV-starting_index))     #Coeffeicient de la pénalité si le journal commence trop tôt
        a_3 = 1 *(0.9**(index_PTV-starting_index))     #Coeffeicent de péanlité pour ne pas détecter une fin de programme
        a_4 = 20*(0.9**(index_PTV-starting_index))     #Paramètre de régularisation -- utilisation des points les plus probables.
        a_5 = 2 *(0.9**(index_PTV-starting_index))     #Pénalitasion du la valeur absolue de l'éloignement du point prévu de fin.
        #Update time of commercials (Reset)
        if(i%60 == 0):
            Pubinhour = 0
        #Update timmers
        lastPub+=1
        lastCP+=1
        if(index_ipts==len(importantpts)):
            index_ipts-=1
        #let's get the context:
        context = def_context.get_context(i,PTV.iloc[index_PTV],lastCP,lastPub,lastend,currentduree,planifiedend,Pubinhour,proba,nbpub,chaine,per,PTV,index_PTV,date)
        ###### Let's verify that the algo is not doing a crappy predicitions and if this the case, clean his historic #####
        if(i == importantpts[index_ipts][0]):
            #### we are at an important point, let's now see what the algo has predict
            if(PTV['TITRE'].iloc[index_PTV] == importantpts[index_ipts][1]):
                #Well he doesn't have the programme wrong, that's a good start
                #let's now find out if we are at a logical point of the programme
                if(i-lastend>13):
                    #Wellllll, the programme began way too early...something went wrong before...Let's rest for now, we'll correct the algo later
                    error_type_2+= a_0

                else:
                    # OMG the ALGO IS RIGHT...here is a candy, let's rest a litle just in case...we never know....
                    nbpub = 0




            else:
                #maybe it's the nexTirage du Lotot programme so calme the fuck down!
                if(PTV['TITRE'].iloc[(index_PTV+1)%PTV.shape[0]] == importantpts[index_ipts][1]):
                    if(planifiedend-i<13):
                        per = context[3]
                        if(context[14]-int(context[14])):
                            error_type_1+= abs(context[14])*a_1
                        error_type_1 += abs((1-context[3])*context[11])*a_5
                        if(context[14]<0):
                            error_type_1+= abs(context[14])*50
                        index_ipts+=1
                    else:

                        error_type_2+= a_0


                else:
                    #well the programme is wrong, and we are not even close to it, let's terminate this thing before it goes completly south. REBOOT The algo, erase the memory, just like in Westworld.
                    #BUT FIRST LET'S VERIFY THAT THERE IS INDEED AN IMPORTANT PROGRAMME THAT DAY...Don't go fuck everything up for no reason
                    l = PTV.index[(PTV['TITRE']==importantpts[index_ipts][1]) & (PTV['debut'] == i)].tolist()
                    if(len(l)>0):

                        error_type_2+= a_0

                    else:
                        print("PAS DE POINT DE CONTRÔLE POUR CETTE PARTIE DE LA JOURNÉE")
                        index_ipts+=1

        elif(context[2]):
            if(lastCP < min(6,currentduree/2)):
                if(lastCP == 1):
                    continue
                else:
                    error_type_3+= abs((context[13]-THRESHOLD)*a_4)

            else:
                X = def_context.process(pd.DataFrame([context],index=[0],columns=['minute','partie de la journée','Change Point','pourcentage','partie du programme','programme','duree','nombre de pub potentiel','lastCP','lastPub','lastend','currentduree','Pubinhour','probability of CP','nb de pubs encore possible','chaine','CLE-FORMAT','CLE-GENRE','day','part'])).values #,'per'
                res1 = CatBoost[0].predict_proba(X)
                res2 = CatBoost[1].predict_proba(X)
                res3 = XGB[0].predict(xgb.DMatrix(X), ntree_limit= XGB[0].best_ntree_limit)
                res4 = XGB[1].predict(xgb.DMatrix(X), ntree_limit= XGB[1].best_ntree_limit)
                res5 = rf.predict_proba(X)
                res6 = gb.predict_proba(X)
                res7 = dt.predict_proba(X)
                res = [(res1[0][0]+res2[0][0]+res3[0][0]+res4[0][0]+res5[0][0]+res6[0][0])/6,(res1[0][1]+res2[0][1]+res3[0][1]+res4[0][1]+res5[0][1]+res6[0][1])/6,(res1[0][2]+res2[0][2]+res3[0][2]+res4[0][2]+res5[0][2]+res6[0][2])/6]
                y_pred = [(res1[0][0]+res2[0][0])*0.5,(res1[0][1]+res2[0][1])*0.5,(res1[0][2]+res2[0][2])*0.5]
                y_pred2 = [(res3[0][0]+res4[0][0])*0.5,(res3[0][1]+res4[0][1])*0.5,(res3[0][2]+res4[0][2])*0.5]
                X = pd.concat([pd.DataFrame(y_pred).T,pd.DataFrame(y_pred2).T,pd.DataFrame(res7),pd.DataFrame(res6),pd.DataFrame(res5)],axis = 1)
                X = X.replace([np.inf, -np.inf], np.nan)
                X = X.fillna(1)
                X = X.values
                res = logistic.predict_proba(X)
                cla = np.argmax(res)

                if(cla == 1 and context[14]==0):
                    cla = np.argmax([res[0][0],res[0][2]])*2
                if(cla == 2 and context[3]<0.75 and context[11]>=20):
                    cla = 0
                if(cla == 2 and context[3]<0.9 and context[11]>=180):
                    cla = 0
                if(context[3]>1):
                    cla = 2
                if(cla == 1):
                    lastCP=0
                    lastPub = 0
                    Pubinhour+=4
                    nbpub+=1
                    wait = 4
                elif(cla == 2):
                    lastend = i
                    lastCP=0
                    index_PTV += 1
                    index_PTV = index_PTV%(PTV.shape[0])
                    currentduree = PTV['DUREE'].iloc[index_PTV]
                    planifiedend = (lastend + currentduree)
                    Predictiontimer = 200
                    nbpub = 0
                    wait = 5
                    per = context[3]
                    if(context[14]-int(context[14])):
                        error_type_1+= abs(context[14])*a_1
                    error_type_1 += abs((1-context[3])*context[11])*a_5
                    if(context[14]<0):
                        error_type_1+= abs(context[14])*a_4
                    if(PTV['TITRE'].iloc[index_PTV] == importantpts[index_ipts][1]):
                        if(importantpts[index_ipts][0]-i>13):
                            error_type_2+= (abs(importantpts[index_ipts][0]-i))*a_2
                        else:

                            if(chaine == 'M6'):
                                error_type_2 -= importantpts[index_ipts][0]-i
                            else:
                                pass

                else:
                    error_type_3+= abs((context[13]-THRESHOLD)*a_4)



        elif(i in help):
            if(lastCP < min(6,currentduree)):
                error_type_3-= 5
            else:
                X = def_context.process(pd.DataFrame([context],index=[0],columns=['minute','partie de la journée','Change Point','pourcentage','partie du programme','programme','duree','nombre de pub potentiel','lastCP','lastPub','lastend','currentduree','Pubinhour','probability of CP','nb de pubs encore possible','chaine','CLE-FORMAT','CLE-GENRE','day','part'])).values #,'per'
                res1 = CatBoost[0].predict_proba(X)
                res2 = CatBoost[1].predict_proba(X)
                res3 = XGB[0].predict(xgb.DMatrix(X), ntree_limit= XGB[0].best_ntree_limit)
                res4 = XGB[1].predict(xgb.DMatrix(X), ntree_limit= XGB[1].best_ntree_limit)
                res5 = rf.predict_proba(X)
                res6 = gb.predict_proba(X)
                res7 = dt.predict_proba(X)
                res = [(res1[0][0]+res2[0][0]+res3[0][0]+res4[0][0]+res5[0][0]+res6[0][0])/6,(res1[0][1]+res2[0][1]+res3[0][1]+res4[0][1]+res5[0][1]+res6[0][1])/6,(res1[0][2]+res2[0][2]+res3[0][2]+res4[0][2]+res5[0][2]+res6[0][2])/6]
                y_pred = [(res1[0][0]+res2[0][0])*0.5,(res1[0][1]+res2[0][1])*0.5,(res1[0][2]+res2[0][2])*0.5]
                y_pred2 = [(res3[0][0]+res4[0][0])*0.5,(res3[0][1]+res4[0][1])*0.5,(res3[0][2]+res4[0][2])*0.5]
                X = pd.concat([pd.DataFrame(y_pred).T,pd.DataFrame(y_pred2).T,pd.DataFrame(res7),pd.DataFrame(res6),pd.DataFrame(res5)],axis = 1)
                X = X.replace([np.inf, -np.inf], np.nan)
                X = X.fillna(1)
                X = X.values
                res = logistic.predict_proba(X)
                cla = np.argmax(res)

                if(cla == 1 and context[14]==0):
                    cla = np.argmax([res[0][0],res[0][2]])*2
                if(cla == 2 and context[3]<0.75 and context[11]>=20):
                    cla = 0
                if(cla == 2 and context[3]<0.9 and context[11]>=180):
                    cla = 0
                if(context[3]>1):
                    cla = 2


                if(cla == 1):
                    lastCP=0
                    lastPub = 0
                    Pubinhour+=4
                    nbpub+=1
                    wait = 4
                elif(cla == 2):
                    lastend = i
                    lastCP=0
                    index_PTV += 1
                    index_PTV = index_PTV%(PTV.shape[0])
                    currentduree = PTV['DUREE'].iloc[index_PTV]
                    planifiedend = (lastend + currentduree)
                    Predictiontimer = 200
                    nbpub = 0
                    wait = 5
                    per = context[3]
                    if(context[14]-int(context[14])):
                        error_type_1+= abs(context[14])*a_1
                    error_type_1 += abs((1-context[3])*context[11])*a_5
                    if(context[14]<0):
                        error_type_1+= abs(context[14])*a_4
                    if(PTV['TITRE'].iloc[index_PTV] == importantpts[index_ipts][1]):
                        if(importantpts[index_ipts][0]-i>13):
                            error_type_2+= (abs(importantpts[index_ipts][0]-i))*a_2
                        else:
                            if(chaine == 'M6'):
                                error_type_2 -= importantpts[index_ipts][0]-i
                            else:
                                pass

                else:
                    error_type_3-=5


        else:
            #labels.append(0)
            #Not a Change Point, we'll just check that nothing is wrong in the PTV at this time
            if(Predictiontimer <= 0):
                l = currentduree
                lastend = i
                lastCP=0
                index_PTV += 1
                index_PTV = index_PTV%(PTV.shape[0])
                currentduree = PTV['DUREE'].iloc[index_PTV]
                planifiedend = (lastend + currentduree)
                Predictiontimer = 200
                nbpub = 0
                if(context[14]-int(context[14])):
                    error_type_1+= abs(context[14])*a_1
                if(context[14]<0):
                    error_type_1+= abs(context[14])*a_4
                error_type_1 += abs((1-context[3])*context[11])*a_5
                error_type_1 += max(currentduree/5,1)*a_3
                per = context[3]
                if(PTV['TITRE'].iloc[index_PTV] == importantpts[index_ipts][1]):
                    if(importantpts[index_ipts][0]-i>13):
                        error_type_2+= (abs(importantpts[index_ipts][0]-i))*a_2
                    else:
                        if(chaine == 'M6'):
                            error_type_2 -= importantpts[index_ipts][0]-i
                        else:
                            pass

            elif(context[3] == 1 or (context[3]>1 and Predictiontimer>20)):
                #Dépassement autorisé: Modulable en fonction de la position dans la journée si besoin
                # C'est sur ces valeurs que l'on va jouer pour avoir le meilleur PTV possible
                # Plus les valeurs sont grandes, plus on fait confiance a l'algo
                # Il est important de bien découper la journée celon les périodes horaires que l'on qualifie
                # de "sous tension" si plusieurs programmes courts se succédent. Bien évidement une telle analyse sera
                #plus tard fait automatiquement.
                if(chaine == 'TF1'):
                    if(context[5] == 'Journal'):
                        if(i<20*60):
                            Predictiontimer = 10
                        else:
                            Predictiontimer = 0

                    elif(context[6] == "très court"):
                        Predictiontimer = 4
                    elif(PTV['TITRE'].iloc[index_PTV] == 'Téléshopping'):
                        Predictiontimer = 2
                    elif(context[6] == "court"):
                        Predictiontimer = 5
                    elif(context[6] == "moyen"):
                        Predictiontimer = 5
                    elif(context[6] == "très long" or context[6] == "long"):
                        Predictiontimer = 15
                    else:
                        Predictiontimer = 15
                    if(i>25*60+30):
                        Predictiontimer =2
                elif(chaine =='M6'):
                    #Dépassement autorisé: Modulable en fonction de la position dans la journée si besoin
                    #Dépassement autorisé: Modulable en fonction de la position dans la journée si besoin
                    if(i>25*60+30):
                        Predictiontimer = 2
                    elif(13*60<i<14*60):
                        Predictiontimer = 5
                    elif(PTV['TITRE'].iloc[index_PTV] in ['M6 boutique']):
                        Predictiontimer = 0
                    elif(context[6] == "très court"):
                        Predictiontimer = 15
                    elif(context[6] == "court"):
                        Predictiontimer = 15
                    elif(context[6] == "moyen"):
                        Predictiontimer = 15
                    elif(context[6] == "très long"):
                        Predictiontimer = 15
                    elif(context[6] == 'long'):
                        Predictiontimer = 15
                    else:
                        Predictiontimer = 5
                else:
                    if(context[5] == 'Journal'):
                        if(i<20*60):
                            Predictiontimer = 10
                        else:
                            Predictiontimer = 0
                    elif(context[6] == "très court"):
                        Predictiontimer = 4
                    elif(context[6] == "court"):
                        Predictiontimer = 5
                    elif(context[6] == "moyen"):
                        Predictiontimer = 5
                    elif(context[6] == "très long"):
                        Predictiontimer = 5
                    elif(context[6] == 'long'):
                        Predictiontimer = 15
                    else:
                        Predictiontimer = 5
            elif(context[3]>1):
                Predictiontimer -= 1
            else:
                pass
    error+= error_type_1 + error_type_2 + error_type_3+error_type_4
    return error,[error_type_1, error_type_2,error_type_3,error_type_4]

def most_pobable_path(error_0,error_1,error_2,XGB,CatBoost,rf,dt,gb,logistic,context):
    g = np.argmin([error_0,error_1,error_2])
    min = ([error_0,error_1,error_2])[g]
    res = [error_0,error_1,error_2]
    res.pop(g)
    if(error_0 != error_1 and error_0 != error_2 and error_1 != error_2):
        cla = g
    elif(min not in res):
        cla = g
    elif(g == 0 and context[13]<0.65):
        cla = 0
    else:
        possibles = []
        if(error_0 <= min):
            possibles.append(0)
        if(error_1 <= min):
            possibles.append(1)
        if(error_2 <= min):
            possibles.append(2)
        if(min >= 100):
            print("pas de bonne solution",str(int(context[0]/60))+':'+str(context[0]%60))
        if(len(possibles) == 0):
            print("pas de bonne solution",str(int(context[0]/60))+':'+str(context[0]%60))
            possibles = [0,1,2]

        X = def_context.process(pd.DataFrame([context],index=[0],columns=['minute','partie de la journée','Change Point','pourcentage','partie du programme','programme','duree','nombre de pub potentiel','lastCP','lastPub','lastend','currentduree','Pubinhour','probability of CP','nb de pubs encore possible','chaine','CLE-FORMAT','CLE-GENRE','day','part'])).values #,'per'
        res1 = CatBoost[0].predict_proba(X)
        res2 = CatBoost[1].predict_proba(X)
        res3 = XGB[0].predict(xgb.DMatrix(X), ntree_limit= XGB[0].best_ntree_limit)
        res4 = XGB[1].predict(xgb.DMatrix(X), ntree_limit= XGB[1].best_ntree_limit)
        res5 = rf.predict_proba(X)
        res6 = gb.predict_proba(X)
        res7 = dt.predict_proba(X)
        res = [(res1[0][0]+res2[0][0]+res3[0][0]+res4[0][0]+res5[0][0]+res6[0][0])/6,(res1[0][1]+res2[0][1]+res3[0][1]+res4[0][1]+res5[0][1]+res6[0][1])/6,(res1[0][2]+res2[0][2]+res3[0][2]+res4[0][2]+res5[0][2]+res6[0][2])/6]
        y_pred = [(res1[0][0]+res2[0][0])*0.5,(res1[0][1]+res2[0][1])*0.5,(res1[0][2]+res2[0][2])*0.5]
        y_pred2 = [(res3[0][0]+res4[0][0])*0.5,(res3[0][1]+res4[0][1])*0.5,(res3[0][2]+res4[0][2])*0.5]
        X = pd.concat([pd.DataFrame(y_pred).T,pd.DataFrame(y_pred2).T,pd.DataFrame(res7),pd.DataFrame(res6),pd.DataFrame(res5)],axis = 1)
        X = X.replace([np.inf, -np.inf], np.nan)
        X = X.fillna(1)
        X = X.values
        res = logistic.predict_proba(X)

        cl = [[res[0][i],i] for i in possibles]
        cla = cl[np.argmax([l[0] for l in cl])][1]

    if(cla == 2 and context[3]<0.65 and context[11]>=20):
        cla = most_pobable_path(error_0,error_1,np.Infinity,XGB,CatBoost,rf,dt,gb,logistic,context)
    if(cla == 2 and context[3]<0.9 and context[11]>=180):
        cla = most_pobable_path(error_0,error_1,np.Infinity,XGB,CatBoost,rf,dt,gb,logistic,context)
    if(cla == 1 and context[9]<0):
        cla = most_pobable_path(error_0,np.Infinity,error_2,XGB,CatBoost,rf,dt,gb,logistic,context)
    if(cla == 2 and context[3]<0.75 and context[11]>=20 and context[0]>19*60+50):
        cla = most_pobable_path(error_0,error_1,np.Infinity,XGB,CatBoost,rf,dt,gb,logistic,context)
    if(cla == 2 and context[3]<0.95 and 180<context[0]<9*60 and context[11]>=20):
        cla = most_pobable_path(error_0,error_1,np.Infinity,XGB,CatBoost,rf,dt,gb,logistic,context)
    if(cla == 1 and context[14]==0):
        cla = most_pobable_path(error_0,np.Infinity,error_2,XGB,CatBoost,rf,dt,gb,logistic,context)
    return cla





def make_newPTV(PTV,proba,chaine,index,lastPTV,lastcontext,index_PTV,importantpts,date,path):
    #Initialisation des Variables
    global THRESHOLD
    automatic_predtimer = True
    verbose = False
    index_PTV = index_PTV
    ##########################
    Predictiontimer = 200
    Pubinhour = lastcontext[12]
    lastCP = lastcontext[8]
    lastPub= lastcontext[9]
    lastend = lastcontext[10]
    currentduree = lastcontext[11]
    planifiedend = lastcontext[10]+lastcontext[11]
    begin = True
    nbpub = 0
    Recall = 1
    wait = 4
    error = 0
    per = 1
    importantpts = importantpts
    help = def_context.get_help(chaine,PTV)
    index_ipts = index
    newPTV = def_context.init_newPTV(PTV,chaine)
    historyofpoints = def_context.init_history(chaine,PTV,lastend,currentduree)
    ######################
    historyofpoints.loc[0] = lastcontext
    labels = [0]
    start = lastcontext[0]+1
    end = importantpts[index][0]
    print(str(start)+' '+str(end))
    #########init Classifier#############
    XGB,CatBoost,rf,dt,gb,logistic = def_context.load_models(path)
    ####################################
    for i in tqdm(range(start,end+5)):
        #Update time of commercials (Reset)
        if(i == end+5 and index == 2):
            newPTV.loc[newPTV.shape[0]] = [(i+currentduree)%1440,PTV['TITRE'].iloc[index_PTV],'non',1,"fin d'un programme"]
        if(i%60 == 0):
            Pubinhour = 0
        #Update timmers
        lastPub+=1
        lastCP+=1
        if(index_ipts==len(importantpts)):
            index_ipts-=1
        #let's get the context:
        context = def_context.get_context(i,PTV.iloc[index_PTV],lastCP,lastPub,lastend,currentduree,planifiedend,Pubinhour,proba,nbpub,chaine,per,PTV,index_PTV,date)
        #Sur M6 il y a 16 minutes de pub entre deux films!!!!!!!!!!!!.....!!!!!!!....!!.!.!.!.!....!.!...!..!.!.!.!
        if(PTV['GENRESIMPLE'].iloc[index_PTV].split(' ')[0] == PTV['GENRESIMPLE'].iloc[index_PTV-1].split(' ')[0] and PTV['GENRESIMPLE'].iloc[index_PTV].split(' ')[0] == 'Téléfilm'
            and (i-lastend)<2 and Recall > 0 and per<0.97 and chaine == 'M6'):
            print(i,PTV['GENRESIMPLE'].iloc[index_PTV])
            lastend = i+5
            lastPub = -25
            Recall -= 0.5
        elif((i-lastend)<2 and Recall > 0 and per<0.97 and chaine == 'M6' and 15*60<i<16*60):
            print(i,PTV['GENRESIMPLE'].iloc[index_PTV],PTV['GENRESIMPLE'].iloc[index_PTV-1])
            lastend = i+5
            lastPub = -25
            Recall -= 0.5

        ###### Let's verify that the algo is not doing a crappy predicitions and if this the case, clean his historic #####
        elif(i == importantpts[index_ipts][0]):
            if(3*60<i<22*60):
                #### we are at an important point, let's now see what the algo has predict
                if(PTV['TITRE'].iloc[index_PTV] == importantpts[index_ipts][1]):
                    #Well he doesn't have the programme wrong, that's a good start
                    #let's now find out if we are at a logical point of the programme
                    if(i-lastend>13):
                        #Wellllll, the programme began way too early...something went wrong before...Let's rest for now, we'll correct the algo later
                        Predictiontimer = 200
                        Pubinhour = 0
                        lastCP = 0
                        lastPub = 0
                        lastend = i
                        currentduree = PTV['DUREE'].iloc[index_PTV]
                        planifiedend = (lastend + currentduree)
                        nbpub = 0
                        if(index_ipts == 0):
                            print("erreur sur la matinée")
                        elif(index_ipts == 1):
                            print("erreur sur l'après midi")
                        else:
                            print("erreur sur la soirée")
                        error+=1
                        #we can now keep going throw the process like before
                        #we just add a line to the history to say that a reset occured
                        newPTV.loc[newPTV.shape[0]] = [i%1440,PTV['TITRE'].iloc[index_PTV],'non',context[3],"--HARD RESET OF ALGORITHM--(in programme)"]
                        index_ipts+=1
                    else:
                        # OMG the ALGO IS RIGHT...here is a candy, let's rest a litle just in case...we never know....
                        Predictiontimer = 200
                        Pubinhour = 0
                        lastCP = 0
                        lastPub= 0
                        lastend = i
                        currentduree = PTV['DUREE'].iloc[index_PTV]
                        planifiedend = (lastend + currentduree)
                        nbpub = 0
                        #we can now keep going throw the process like before
                        #we just add a line to the history to say that a reset occured
                        newPTV.loc[newPTV.shape[0]] = [i%1440,PTV['TITRE'].iloc[index_PTV],'non',context[3],"--soft reset to avoid any error--"]
                        index_ipts+=1



                else:
                    #maybe it's the next programme so calme the fuck down!
                    if(PTV['TITRE'].iloc[(index_PTV+1)%PTV.shape[0]] == importantpts[index_ipts][1]):
                        if(planifiedend-i<10):
                            #here you go, it's the next one...just terminate this one and we're good to go
                            newPTV.loc[newPTV.shape[0]] = [i%1440,PTV['TITRE'].iloc[index_PTV],'oui',context[3],"fin d'un programme"]
                            lastend = i
                            lastCP=0
                            index_PTV += 1
                            index_PTV = index_PTV%(PTV.shape[0])
                            currentduree = PTV['DUREE'].iloc[index_PTV]
                            planifiedend = (lastend + currentduree)
                            Predictiontimer = 200
                            nbpub = 0
                            index_ipts+=1
                        else:
                            #here you go, it's the next one...just terminate this one and we're good to go
                            newPTV.loc[newPTV.shape[0]] = [i%1440,PTV['TITRE'].iloc[index_PTV],'oui',context[3],"--HARD RESET OF ALGORITHM--(Oups)"]
                            lastend = i
                            lastCP=0
                            index_PTV += 1
                            index_PTV = index_PTV%(PTV.shape[0])
                            currentduree = PTV['DUREE'].iloc[index_PTV]
                            planifiedend = (lastend + currentduree)
                            Predictiontimer = 200
                            nbpub = 0
                            if(index_ipts == 0):
                                print("erreur sur la matinée")
                            elif(index_ipts == 1):
                                print("erreur sur l'après midi")
                            else:
                                print("erreur sur la soirée")
                            error+=1
                            index_ipts+=1



                    else:
                        #well the programme is wrong, and we are not even close to it, let's terminate this thing before it goes completly south. REBOOT The algo, erase the memory, just like in Westworld.
                        #BUT FIRST LET'S VERIFY THAT THERE IS INDEED AN IMPORTANT PROGRAMME THAT DAY...Don't go fuck everything up for no reason
                        l = PTV.index[(PTV['TITRE']==importantpts[index_ipts][1]) & (PTV['debut'] == i%1440)].tolist()
                        if(len(l)>0):
                            index_PTV = l[0]
                            ##########################
                            Predictiontimer = 200
                            Pubinhour = 0
                            lastCP = 0
                            lastPub= 0
                            lastend = i
                            currentduree = PTV['DUREE'].iloc[index_PTV]
                            planifiedend = (lastend + currentduree)
                            nbpub = 0
                            if(index_ipts == 0):
                                print("erreur sur la matinée")
                            elif(index_ipts == 1):
                                print("erreur sur l'après midi")
                            else:
                                print("erreur sur la soirée")
                            error+=1
                            #we can now keep going throw the process like before
                            #we just add a line to the history to say that a reset occured
                            newPTV.loc[newPTV.shape[0]] = [i%1440,PTV['TITRE'].iloc[index_PTV],'non',context[3],"--HARD RESET OF ALGORITHM--(out of programme)"]
                            index_ipts+=1
                        else:
                            index_ipts+=1
            else:
                #### we are at an important point, let's now see what the algo has predict
                if(PTV['TITRE'].iloc[index_PTV] == importantpts[index_ipts][1]):
                    #Well he doesn't have the programme wrong, that's a good start
                    #let's now find out if we are at a logical point of the programme
                    if(i-lastend>20):
                        #Wellllll, the programme began way too early...something went wrong before...Let's rest for now, we'll correct the algo later
                        Predictiontimer = 200
                        Pubinhour = 0
                        lastCP = 0
                        lastPub = 0
                        lastend = i
                        currentduree = PTV['DUREE'].iloc[index_PTV]
                        planifiedend = (lastend + currentduree)
                        nbpub = 0
                        if(index_ipts == 0):
                            print("erreur sur la matinée")
                        elif(index_ipts == 1):
                            print("erreur sur l'après midi")
                        else:
                            print("erreur sur la soirée")
                        error+=1
                        #we can now keep going throw the process like before
                        #we just add a line to the history to say that a reset occured
                        newPTV.loc[newPTV.shape[0]] = [i%1440,PTV['TITRE'].iloc[index_PTV],'non',context[3],"--HARD RESET OF ALGORITHM--(in programme)"]
                        index_ipts+=1
                    else:
                        # OMG the ALGO IS RIGHT...here is a candy, let's rest a litle just in case...we never know....
                        Predictiontimer = 200
                        Pubinhour = 0
                        lastCP = 0
                        lastPub= 0
                        lastend = i
                        currentduree = PTV['DUREE'].iloc[index_PTV]
                        planifiedend = (lastend + currentduree)
                        nbpub = 0
                        #we can now keep going throw the process like before
                        #we just add a line to the history to say that a reset occured
                        newPTV.loc[newPTV.shape[0]] = [i%1440,PTV['TITRE'].iloc[index_PTV],'non',context[3],"--soft reset to avoid any error--"]
                        index_ipts+=1



                else:
                    #maybe it's the next programme so calme the fuck down!
                    if(PTV['TITRE'].iloc[(index_PTV+1)%PTV.shape[0]] == importantpts[index_ipts][1]):
                        if(planifiedend-i<20):
                            #here you go, it's the next one...just terminate this one and we're good to go
                            newPTV.loc[newPTV.shape[0]] = [i%1440,PTV['TITRE'].iloc[index_PTV],'oui',context[3],"fin d'un programme"]
                            lastend = i
                            lastCP=0
                            index_PTV += 1
                            index_PTV = index_PTV%(PTV.shape[0])
                            currentduree = PTV['DUREE'].iloc[index_PTV]
                            planifiedend = (lastend + currentduree)
                            Predictiontimer = 200
                            nbpub = 0
                            index_ipts+=1
                        else:
                            #here you go, it's the next one...just terminate this one and we're good to go
                            newPTV.loc[newPTV.shape[0]] = [i%1440,PTV['TITRE'].iloc[index_PTV],'oui',context[3],"--HARD RESET OF ALGORITHM--(Oups)"]
                            lastend = i
                            lastCP=0
                            index_PTV += 1
                            index_PTV = index_PTV%(PTV.shape[0])
                            currentduree = PTV['DUREE'].iloc[index_PTV]
                            planifiedend = (lastend + currentduree)
                            Predictiontimer = 200
                            nbpub = 0
                            if(index_ipts == 0):
                                print("erreur sur la matinée")
                            elif(index_ipts == 1):
                                print("erreur sur l'après midi")
                            else:
                                print("erreur sur la soirée")
                            error+=1
                            index_ipts+=1



                    else:
                        #well the programme is wrong, and we are not even close to it, let's terminate this thing before it goes completly south. REBOOT The algo, erase the memory, just like in Westworld.
                        #BUT FIRST LET'S VERIFY THAT THERE IS INDEED AN IMPORTANT PROGRAMME THAT DAY...Don't go fuck everything up for no reason
                        l = PTV.index[(PTV['TITRE']==importantpts[index_ipts][1]) & (PTV['debut'] == i%1440)].tolist()
                        if(len(l)>0):
                            index_PTV = l[0]
                            ##########################
                            Predictiontimer = 200
                            Pubinhour = 0
                            lastCP = 0
                            lastPub= 0
                            lastend = i
                            currentduree = PTV['DUREE'].iloc[index_PTV]
                            planifiedend = (lastend + currentduree)
                            nbpub = 0
                            if(index_ipts == 0):
                                print("erreur sur la matinée")
                            elif(index_ipts == 1):
                                print("erreur sur l'après midi")
                            else:
                                print("erreur sur la soirée")
                            error+=1
                            #we can now keep going throw the process like before
                            #we just add a line to the history to say that a reset occured
                            newPTV.loc[newPTV.shape[0]] = [i%1440,PTV['TITRE'].iloc[index_PTV],'non',context[3],"--HARD RESET OF ALGORITHM--(out of programme)"]
                            index_ipts+=1
                        else:
                            index_ipts+=1
        elif(context[2]):
            historyofpoints.loc[historyofpoints.shape[0]] = context
            if(lastCP < min(6,currentduree/2)):
                labels.append(0)
                continue
            else:
                if(i<importantpts[1][0]):
                    # Ne fonctionne qu'avec un Point de Contrôle et donc qu'avec des points avant 20h
                    error_0,sub_error_0 = make_predictedPTV(0,PTV,chaine,i,index_PTV,lastCP,lastPub,lastend,currentduree,planifiedend,Pubinhour,proba,nbpub,XGB,CatBoost,rf,dt,gb,logistic,Predictiontimer,per,context,date,importantpts)
                    error_1,sub_error_1 = make_predictedPTV(1,PTV,chaine,i,index_PTV,lastCP,lastPub,lastend,currentduree,planifiedend,Pubinhour,proba,nbpub,XGB,CatBoost,rf,dt,gb,logistic,Predictiontimer,per,context,date,importantpts)
                    error_2,sub_error_2 = make_predictedPTV(2,PTV,chaine,i,index_PTV,lastCP,lastPub,lastend,currentduree,planifiedend,Pubinhour,proba,nbpub,XGB,CatBoost,rf,dt,gb,logistic,Predictiontimer,per,context,date,importantpts)
                    cla = most_pobable_path(error_0,error_1,error_2,XGB,CatBoost,rf,dt,gb,logistic,context)
                    print(context[3])
                    #print(error_0,sub_error_0,'|',error_1,sub_error_1,'|',error_2,sub_error_2,'|',str(int(context[0]/60))+':'+str(context[0]%60),cla)
                    print(error_0,sub_error_0,'|',error_1,sub_error_1,'|',error_2,sub_error_2,'|',str(int(context[0]/60))+':'+str(context[0]%60),cla)


                else:
                    error_0,sub_error_0 = make_predictedPTV2(0,PTV,chaine,i,index_PTV,lastCP,lastPub,lastend,currentduree,planifiedend,Pubinhour,proba,nbpub,XGB,CatBoost,rf,dt,gb,logistic,Predictiontimer,per,context,date,importantpts)
                    error_1,sub_error_1 = make_predictedPTV2(1,PTV,chaine,i,index_PTV,lastCP,lastPub,lastend,currentduree,planifiedend,Pubinhour,proba,nbpub,XGB,CatBoost,rf,dt,gb,logistic,Predictiontimer,per,context,date,importantpts)
                    error_2,sub_error_2 = make_predictedPTV2(2,PTV,chaine,i,index_PTV,lastCP,lastPub,lastend,currentduree,planifiedend,Pubinhour,proba,nbpub,XGB,CatBoost,rf,dt,gb,logistic,Predictiontimer,per,context,date,importantpts)
                    cla = most_pobable_path(error_0,error_1,error_2,XGB,CatBoost,rf,dt,gb,logistic,context)
                    print()
                    #print(error_0,sub_error_0,'|',error_1,sub_error_1,'|',error_2,sub_error_2,'|',str(int(context[0]/60))+':'+str(context[0]%60),cla)
                    print(error_0,sub_error_0,'|',error_1,sub_error_1,'|',error_2,sub_error_2,'|',str(int(context[0]/60))+':'+str(context[0]%60),cla)


                if(cla == 1):
                    newPTV.loc[newPTV.shape[0]] = [i%1440,"publicité",'oui',context[3],"publicité dans un programme"]
                    lastCP=0
                    lastPub = 0
                    Pubinhour+=4
                    nbpub+=1
                    wait = 4
                    labels.append(1)
                elif(cla == 2):
                    newPTV.loc[newPTV.shape[0]] = [i%1440,PTV['TITRE'].iloc[index_PTV],'oui',context[3],"fin d'un programme"]
                    lastend = i
                    lastCP=0
                    index_PTV += 1
                    index_PTV = index_PTV%(PTV.shape[0])
                    currentduree = PTV['DUREE'].iloc[index_PTV]
                    planifiedend = (lastend + currentduree)
                    Predictiontimer = 200
                    nbpub = 0
                    wait = 5
                    per = context[3]
                    labels.append(2)
                else:
                    labels.append(0)




        elif(i in help):
            historyofpoints.loc[historyofpoints.shape[0]] = context
            if(lastCP < min(6,currentduree)):
                labels.append(0)
                continue
            else:
                if(i<importantpts[1][0]):
                    print("utilisation de l'algo 1")
                    error_0,sub_error_0 = make_predictedPTV(0,PTV,chaine,i,index_PTV,lastCP,lastPub,lastend,currentduree,planifiedend,Pubinhour,proba,nbpub,XGB,CatBoost,rf,dt,gb,logistic,Predictiontimer,per,context,date,importantpts)
                    error_1,sub_error_1 = make_predictedPTV(1,PTV,chaine,i,index_PTV,lastCP,lastPub,lastend,currentduree,planifiedend,Pubinhour,proba,nbpub,XGB,CatBoost,rf,dt,gb,logistic,Predictiontimer,per,context,date,importantpts)
                    error_2,sub_error_2 = make_predictedPTV(2,PTV,chaine,i,index_PTV,lastCP,lastPub,lastend,currentduree,planifiedend,Pubinhour,proba,nbpub,XGB,CatBoost,rf,dt,gb,logistic,Predictiontimer,per,context,date,importantpts)
                    cla = most_pobable_path(error_0,error_1,error_2,XGB,CatBoost,rf,dt,gb,logistic,context)
                    print()
                    #print(error_0,sub_error_0,'|',error_1,sub_error_1,'|',error_2,sub_error_2,'|',str(int(context[0]/60))+':'+str(context[0]%60),cla)
                    print(error_0,sub_error_0,'|',error_1,sub_error_1,'|',error_2,sub_error_2,'|',str(int(context[0]/60))+':'+str(context[0]%60),cla)
                else:
                    print("utilisation de l'algo 2")
                    error_0,sub_error_0 = make_predictedPTV2(0,PTV,chaine,i,index_PTV,lastCP,lastPub,lastend,currentduree,planifiedend,Pubinhour,proba,nbpub,XGB,CatBoost,rf,dt,gb,logistic,Predictiontimer,per,context,date,importantpts)
                    error_1,sub_error_1 = make_predictedPTV2(1,PTV,chaine,i,index_PTV,lastCP,lastPub,lastend,currentduree,planifiedend,Pubinhour,proba,nbpub,XGB,CatBoost,rf,dt,gb,logistic,Predictiontimer,per,context,date,importantpts)
                    error_2,sub_error_2 = make_predictedPTV2(2,PTV,chaine,i,index_PTV,lastCP,lastPub,lastend,currentduree,planifiedend,Pubinhour,proba,nbpub,XGB,CatBoost,rf,dt,gb,logistic,Predictiontimer,per,context,date,importantpts)
                    cla = most_pobable_path(error_0,error_1,error_2,XGB,CatBoost,rf,dt,gb,logistic,context)
                    print()
                    #print(error_0,sub_error_0,'|',error_1,sub_error_1,'|',error_2,sub_error_2,'|',str(int(context[0]/60))+':'+str(context[0]%60),cla)
                    print(error_0,sub_error_0,'|',error_1,sub_error_1,'|',error_2,sub_error_2,'|',str(int(context[0]/60))+':'+str(context[0]%60),cla)


                if(cla == 1):
                    newPTV.loc[newPTV.shape[0]] = [i%1440,"publicité",'oui',context[3],"publicité dans un programme"]
                    lastCP=0
                    lastPub = 0
                    Pubinhour+=4
                    nbpub+=1
                    wait = 4
                    labels.append(1)
                elif(cla == 2):
                    newPTV.loc[newPTV.shape[0]] = [i%1440,PTV['TITRE'].iloc[index_PTV],'oui',context[3],"fin d'un programme"]
                    lastend = i
                    lastCP=0
                    index_PTV += 1
                    index_PTV = index_PTV%(PTV.shape[0])
                    currentduree = PTV['DUREE'].iloc[index_PTV]
                    planifiedend = (lastend + currentduree)
                    Predictiontimer = 200
                    nbpub = 0
                    wait = 5
                    per = context[3]
                    labels.append(2)
                else:
                    labels.append(0)



        else:
            #labels.append(0)
            #Not a Change Point, we'll just check that nothing is wrong in the PTV at this time
            if(Predictiontimer <= 0):
                labels.append(2)
                historyofpoints.loc[historyofpoints.shape[0]] = context
                newPTV.loc[newPTV.shape[0]] = [i%1440,PTV['TITRE'].iloc[index_PTV],'non',context[3],"fin non détectée d'un programme"]
                lastend = i
                lastCP=0
                index_PTV += 1
                index_PTV = index_PTV%(PTV.shape[0])
                currentduree = PTV['DUREE'].iloc[index_PTV]
                planifiedend = (lastend + currentduree)
                Predictiontimer = 200
                nbpub = 0
                per = context[3]
            elif(context[3] == 1 or (context[3]>1 and Predictiontimer>20)):
                if(automatic_predtimer):
                    try:
                        def_context.Report("Utilisation de la prediction de dépassement")
                        preds = []
                        if(i<importantpts[1][0]):
                            for p in range(20):
                                error_0,sub_error_0 = make_predictedPTV(0,PTV,chaine,i,index_PTV,lastCP,lastPub,lastend,currentduree,planifiedend,Pubinhour,proba,nbpub,XGB,CatBoost,rf,dt,gb,logistic,p,per,context,date,importantpts)
                                preds.append(error_0)


                        else:
                            for p in range(20):
                                error_0,sub_error_0 = make_predictedPTV2(0,PTV,chaine,i,index_PTV,lastCP,lastPub,lastend,currentduree,planifiedend,Pubinhour,proba,nbpub,XGB,CatBoost,rf,dt,gb,logistic,p,per,context,date,importantpts)
                                preds.append(error_0)

                        Predictiontimer = np.argmin(preds)
                        def_context.Report("best Predictiontimer found: "+str(Predictiontimer) +" "+str(preds))
                    except Exception as e:
                        def_context.Report("Failed to get prediction timmer:%s " %(e))
                        #Dépassement autorisé: Modulable en fonction de la position dans la journée si besoin
                        # C'est sur ces valeurs que l'on va jouer pour avoir le meilleur PTV possible
                        # Plus les valeurs sont grandes, plus on fait confiance a l'algo
                        # Il est important de bien découper	15	215	577	30	0	0.490790	1	TF1	3	1026	9:52	2 la journée celon les périodes horaires que l'on qualifie
                        # de "sous tension" si plusieurs programmes courts se succédent. Bien évidement une telle analyse sera
                        #plus tard fait automatiquement.
                        if(i<20*60+30):
                            if(chaine == 'TF1'):
                                if(11.5*60<=i<=14*60 or 19.5*60<i<21*60):
                                    Predictiontimer = 1
                                elif(context[6] == "très court"):
                                    Predictiontimer = 0
                                elif(PTV['TITRE'].iloc[index_PTV] == 'Téléshopping'):
                                    Predictiontimer = 5
                                elif(context[6] == "court"):
                                    Predictiontimer = 5
                                elif(context[6] == "moyen"):
                                    Predictiontimer = 5
                                elif(context[6] == "très long" or context[6] == "long"):
                                    Predictiontimer = 15
                                else:
                                    Predictiontimer = 5
                            elif(chaine =='M6'):
                                #Dépassement autorisé: Modulable en fonction de la position dans la journée si besoin
                                #Dépassement autorisé: Modulable en fonction de la position dans la journée si besoin
                                if(i<8*60+56):
                                    Predictiontimer = 0
                                elif(13*60<i<14*60):
                                    Predictiontimer = 5
                                elif(PTV['TITRE'].iloc[index_PTV] in ['M6 boutique']):
                                    Predictiontimer = 0
                                elif(context[6] == "très court"):
                                    Predictiontimer = 0
                                elif(context[6] == "court"):
                                    Predictiontimer = 2
                                elif(context[6] == "moyen"):
                                    Predictiontimer = 5
                                elif(context[6] == "très long"):
                                    Predictiontimer = 5
                                elif(context[6] == 'long'):
                                    Predictiontimer = 15
                                else:
                                    Predictiontimer = 5
                        else:
                            if(chaine == 'TF1'):
                                if(context[5] == 'Journal'):
                                    if(i<20*60):
                                        Predictiontimer = 10
                                    else:
                                        Predictiontimer = 0
                                elif(context[6] == "très court"):
                                    Predictiontimer = 4
                                elif(context[6] == "court"):
                                    Predictiontimer = 5
                                elif(context[6] == "moyen"):
                                    Predictiontimer = 5
                                elif(context[6] == "très long"):
                                    Predictiontimer = 5
                                elif(context[6] == 'long'):
                                    Predictiontimer = 15
                                else:
                                    Predictiontimer = 5
                            elif(chaine =='M6'):
                                #Dépassement autorisé: Modulable en fonction de la position dans la journée si besoin
                                #Dépassement autorisé: Modulable en fonction de la position dans la journée si besoin
                                if(context[6] == "très court"):
                                    Predictiontimer = 15
                                elif(context[6] == "court"):
                                    Predictiontimer = 15
                                elif(context[6] == "moyen"):
                                    Predictiontimer = 15
                                elif(context[6] == "très long"):
                                    Predictiontimer = 15
                                elif(context[6] == 'long'):
                                    Predictiontimer = 15
                                else:
                                    Predictiontimer = 5

                else:
                    #Dépassement autorisé: Modulable en fonction de la position dans la journée si besoin
                    # C'est sur ces valeurs que l'on va jouer pour avoir le meilleur PTV possible
                    # Plus les valeurs sont grandes, plus on fait confiance a l'algo
                    # Il est important de bien découper la journée celon les périodes horaires que l'on qualifie
                    # de "sous tension" si plusieurs programmes courts se succédent. Bien évidement une telle analyse sera
                    #plus tard fait automatiquement.
                    if(i<20*60+30):
                        if(chaine == 'TF1'):
                            if(11.5*60<=i<=14*60 or 19.5*60<i<21*60):
                                Predictiontimer = 1
                            elif(context[6] == "très court"):
                                Predictiontimer = 0
                            elif(PTV['TITRE'].iloc[index_PTV] == 'Téléshopping'):
                                Predictiontimer = 5
                            elif(context[6] == "court"):
                                Predictiontimer = 5
                            elif(context[6] == "moyen"):
                                Predictiontimer = 5
                            elif(context[6] == "très long" or context[6] == "long"):
                                Predictiontimer = 15
                            else:
                                Predictiontimer = 5
                        elif(chaine =='M6'):
                            #Dépassement autorisé: Modulable en fonction de la position dans la journée si besoin
                            #Dépassement autorisé: Modulable en fonction de la position dans la journée si besoin
                            if(i<8*60+56):
                                Predictiontimer = 0
                            elif(13*60<i<14*60):
                                Predictiontimer = 5
                            elif(PTV['TITRE'].iloc[index_PTV] in ['M6 boutique']):
                                Predictiontimer = 0
                            elif(context[6] == "très court"):
                                Predictiontimer = 0
                            elif(context[6] == "court"):
                                Predictiontimer = 2
                            elif(context[6] == "moyen"):
                                Predictiontimer = 5
                            elif(context[6] == "très long"):
                                Predictiontimer = 5
                            elif(context[6] == 'long'):
                                Predictiontimer = 15
                            else:
                                Predictiontimer = 5
                        else:
                            if(context[5] == 'Journal'):
                                if(i<20*60):
                                    Predictiontimer = 10
                                else:
                                    Predictiontimer = 0
                            elif(context[6] == "très court"):
                                Predictiontimer = 4
                            elif(context[6] == "court"):
                                Predictiontimer = 5
                            elif(context[6] == "moyen"):
                                Predictiontimer = 5
                            elif(context[6] == "très long"):
                                Predictiontimer = 5
                            elif(context[6] == 'long'):
                                Predictiontimer = 15
                            else:
                                Predictiontimer = 5
                    else:
                        if(chaine == 'TF1'):
                            if(context[5] == 'Journal'):
                                if(i<20*60):
                                    Predictiontimer = 10
                                else:
                                    Predictiontimer = 0
                            elif(context[6] == "très court"):
                                Predictiontimer = 4
                            elif(context[6] == "court"):
                                Predictiontimer = 5
                            elif(context[6] == "moyen"):
                                Predictiontimer = 5
                            elif(context[6] == "très long"):
                                Predictiontimer = 5
                            elif(context[6] == 'long'):
                                Predictiontimer = 15
                            else:
                                Predictiontimer = 5
                        elif(chaine =='M6'):
                            #Dépassement autorisé: Modulable en fonction de la position dans la journée si besoin
                            #Dépassement autorisé: Modulable en fonction de la position dans la journée si besoin
                            if(context[6] == "très court"):
                                Predictiontimer = 15
                            elif(context[6] == "court"):
                                Predictiontimer = 15
                            elif(context[6] == "moyen"):
                                Predictiontimer = 15
                            elif(context[6] == "très long"):
                                Predictiontimer = 15
                            elif(context[6] == 'long'):
                                Predictiontimer = 15
                            else:
                                Predictiontimer = 5
                        else:
                            if(context[5] == 'Journal'):
                                if(i<20*60):
                                    Predictiontimer = 10
                                else:
                                    Predictiontimer = 0
                            elif(context[6] == "très court"):
                                Predictiontimer = 4
                            elif(context[6] == "court"):
                                Predictiontimer = 5
                            elif(context[6] == "moyen"):
                                Predictiontimer = 5
                            elif(context[6] == "très long"):
                                Predictiontimer = 5
                            elif(context[6] == 'long'):
                                Predictiontimer = 15
                            else:
                                Predictiontimer = 5

            elif(context[3]>1):
                Predictiontimer -= 1
            else:
                pass
    return newPTV,historyofpoints,labels,error,index_PTV,context


#################################################
########### main with options ###################
#################################################

def main(argv):
    global PATH_IN,PATH_SCRIPT,PATH_OUT,THRESHOLD
    PATH_IN,PATH_SCRIPT,PATH_OUT = def_context.get_path()
    if(len(argv) == 2):
        import pandas as pd
        c = argv[0]
        f = argv[1]
        chaine = c
        date = f
        PTV,proba = def_context.load_file(str(f),str(c))
        if(len(PTV) == 0):
            return('Fichier Manquant')
        newPTV = def_context.init_newPTV(PTV,str(c))
        index_PTV = PTV.index[(PTV['debut'] <= 3*60) & (PTV['debut']+PTV['DUREE'] > 3*60+5)].tolist()[0]
        def_context.Report('Starting with: %s'%(PTV['TITRE'].iloc[index_PTV]))
        lastend = PTV['debut'].loc[index_PTV]
        currentduree = PTV['DUREE'].loc[index_PTV]
        historyofpoints = def_context.init_history(str(c),PTV,lastend,currentduree)
        importantpts = def_context.get_important_points(c,PTV,index_PTV)
        temp_context = historyofpoints.iloc[0]
        THRESHOLD = def_context.find_threshold(proba,0.46)
        def_context.Report(THRESHOLD)
        path = def_context.get_temp_path()
        for i in range(3):
            l,temp_newPTV,temp_history,index_PTV,temp_context = main([str(c),str(f),i,newPTV.iloc[newPTV.shape[0]-1],temp_context,index_PTV,importantpts,path])
            if(l == 4):
                pass
            else:
                if(l>0):
                    if(i == 0):
                        print("erreur sur la matinée")
                    elif(i == 1):
                        print("erreur sur l'après midi")
                    else:
                        print("erreur sur la soirée")
                newPTV = pd.concat([newPTV,temp_newPTV.iloc[1:]])
                historyofpoints = pd.concat([historyofpoints,temp_history])
                historyofpoints = historyofpoints[['minute','partie de la journée','Change Point','pourcentage','partie du programme','programme','duree','nombre de pub potentiel','lastCP','lastPub','lastend','currentduree','Pubinhour','probability of CP','nb de pubs encore possible','chaine','CLE-FORMAT','CLE-GENRE','day','part']]
        newPTV['Heure'] = newPTV['minute'].apply(lambda x: str(int(x/60))+':'+str(x%60))
        historyofpoints['Heure'] = historyofpoints['minute'].apply(lambda x: str(int(x/60))+':'+str(x%60))
        newPTV.to_html(PATH_IN+'new_ptv/new_PTV-'+date+'_'+chaine+'.html')
        newPTV.to_csv(PATH_IN+'new_ptv/new_PTV-'+date+'_'+chaine+'.csv',index=False)
        historyofpoints.to_html(PATH_IN+'hop/historyofpoints-'+date+'_'+chaine+'.html')
        historyofpoints.to_csv(PATH_IN+'hop/historyofpoints-'+date+'_'+chaine+'.csv',index=False)
        def_context.Report("%s,%s"%(date,historyofpoints.shape))
        sys.exit(l)
    else:
        chaine = argv[0]
        date = argv[1]
        index = argv[2]
        d = "".join(date.split('-'))
        PTV,proba = def_context.load_file(date,chaine)
        if(len(PTV) == 0):
            return 4,0,0,0,0,0
        new_PTV,historyofpoints,labels,error,index_PTV,temp_context = make_newPTV(PTV,proba,chaine,index,argv[3],argv[4],argv[5],argv[6],date,argv[7])
        print(len(labels),historyofpoints.shape)
        historyofpoints['labels'] = labels
        print(chaine,date,historyofpoints.shape,len(labels))
        return error,new_PTV,historyofpoints,index_PTV,temp_context


if __name__ == "__main__":
    # execute only if run as a script
    main(sys.argv[1:])
