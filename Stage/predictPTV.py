#######################################################################
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


#################################################
########### Important functions #################
#################################################

def make_newPTV(PTV,proba,chaine,index,lastPTV,lastcontext,index_PTV,importantpts,date,path):
    #Initialisation des Variables
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
    index_ipts = index
    importantpts = importantpts
    help = def_context.get_help(chaine,PTV)
    newPTV = def_context.init_newPTV(PTV,chaine)
    historyofpoints = def_context.init_history(chaine,PTV,lastend,currentduree)
    ####################################
    historyofpoints.loc[0] = lastcontext
    labels = [0]
    start = lastcontext[0]+1
    end = importantpts[index][0]
    #########init Classifier#############
    XGB,CatBoost,rf,dt,gb,logistic = def_context.load_models(path)
    ####################################
    for i in tqdm(range(start,min(end+5,1620))):
        if(i == end+5 and index == 2):
            newPTV.loc[newPTV.shape[0]] = [(i+currentduree)%1440,PTV['TITRE'].iloc[index_PTV],'non',1,"fin d'un programme"]
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
        elif((i-lastend)<2 and Recall > 0 and per<0.97 and chaine == 'M6' and 15*60<i<16*60):

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
                            def_context.Report("erreur sur la matinée")
                        elif(index_ipts == 1):
                            def_context.Report("erreur sur l'après midi")
                        else:
                            def_context.Report("erreur sur la soirée")
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
                            #here you go, it's the next one...But it's far far away
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
                                def_context.Report("erreur sur la matinée")
                            elif(index_ipts == 1):
                                def_context.Report("erreur sur l'après midi")
                            else:
                                def_context.Report("erreur sur la soirée")
                            error+=1
                            index_ipts+=1



                    else:
                        #well the programme is wrong, and we are not even close to it, let's terminate this thing before it goes completly south. REBOOT The algo, erase the memory, just like in Westworld.
                        #BUT FIRST LET'S VERIFY THAT THERE IS INDEED AN IMPORTANT PROGRAMME THAT DAY...Don't go fuck everything up for no reason
                        l = PTV.index[(PTV['TITRE']==importantpts[index_ipts][1]) & (PTV['debut'] == i)].tolist()
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
                            #we can now keep going throw the process like before
                            #we just add a line to the history to say that a reset occured
                            if(index_ipts == 0):
                                def_context.Report("erreur sur la matinée")
                            elif(index_ipts == 1):
                                def_context.Report("erreur sur l'après midi")
                            else:
                                def_context.Report("erreur sur la soirée")
                            error+=1
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
                            def_context.Report("erreur sur la matinée")
                        elif(index_ipts == 1):
                            def_context.Report("erreur sur l'après midi")
                        else:
                            def_context.Report("erreur sur la soirée")
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
                            #here you go, it's the next one...But it's far far away
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
                                def_context.Report("erreur sur la matinée")
                            elif(index_ipts == 1):
                                def_context.Report("erreur sur l'après midi")
                            else:
                                def_context.Report("erreur sur la soirée")
                            error+=1
                            index_ipts+=1




                    else:
                        #well the programme is wrong, and we are not even close to it, let's terminate this thing before it goes completly south. REBOOT The algo, erase the memory, just like in Westworld.
                        #BUT FIRST LET'S VERIFY THAT THERE IS INDEED AN IMPORTANT PROGRAMME THAT DAY...Don't go fuck everything up for no reason
                        l = PTV.index[(PTV['TITRE']==importantpts[index_ipts][1]) & (PTV['debut'] == i)].tolist()
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
                            #we can now keep going throw the process like before
                            #we just add a line to the history to say that a reset occured
                            if(index_ipts == 0):
                                def_context.Report("erreur sur la matinée")
                            elif(index_ipts == 1):
                                def_context.Report("erreur sur l'après midi")
                            else:
                                def_context.Report("erreur sur la soirée")
                            error+=1
                            newPTV.loc[newPTV.shape[0]] = [i%1440,PTV['TITRE'].iloc[index_PTV],'non',context[3],"--HARD RESET OF ALGORITHM--(out of programme)"]
                            index_ipts+=1
                        else:
                            index_ipts+=1
        if(context[2]):
            historyofpoints.loc[historyofpoints.shape[0]] = context
            if(lastCP < min(4,currentduree)):
                labels.append(0)
                continue
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
                    cla = 0
                if(cla == 2 and context[3]<0.5 and context[11]>30):
                    cla = 0
                if(cla == 2 and context[3]<0.9 and context[11]>=180):
                    cla = 0
                if( cla == 2 and PTV['TITRE'].loc[index_PTV] == 'Programmes de la nuit' and context[3]<1):
                    cla = 0

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
            if(lastCP < min(4,currentduree)):
                labels.append(0)
                continue
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
                    cla = 0
                if(cla == 2 and context[3]<0.5 and context[11]>30):
                    cla = 0
                if(cla == 2 and context[3]<0.9 and context[11]>=180):
                    cla = 0
                if( cla == 2 and PTV['TITRE'].loc[index_PTV] == 'Programmes de la nuit' and context[3]<1):
                    cla = 0

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
                historyofpoints.loc[historyofpoints.shape[0]] = context
                labels.append(2)
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
            elif(context[3] == 1):
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
    global PATH_IN,PATH_SCRIPT,PATH_OUT
    PATH_IN,PATH_SCRIPT,PATH_OUT = def_context.get_path()
    chaine = argv[0]
    date = argv[1]
    index = argv[2]
    d = "".join(date.split('-'))
    PTV,proba = def_context.load_file(date,chaine)
    if(len(PTV) == 0):
        sys.exit(4)
        return 4,0,0,0,0
    new_PTV,historyofpoints,labels,error,index_PTV,temp_context = make_newPTV(PTV,proba,chaine,index,argv[3],argv[4],argv[5],argv[6],date,argv[7])
    historyofpoints['labels'] = labels
    return error,new_PTV,historyofpoints,index_PTV,temp_context
if __name__ == "__main__":
    # execute only if run as a script
    main(sys.argv[1:])
