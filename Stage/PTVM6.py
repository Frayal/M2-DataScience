#################################################
#created the 01/07/2018 11:56 by Alexis Blanchet#
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
import datetime
from tqdm import tqdm
#################################################
########### Global variables ####################
#################################################
import def_context
PATH_IN = '/home/alexis/Bureau/finalproject/Datas/'
PATH_SCRIPT = '/home/alexis/Bureau/finalproject/scripts/'
PATH_OUT = '/home/alexis/Bureau/finalproject/Datas/'
LOG = "log.txt"
#################################################
########### Important functions #################
#################################################
def make_newPTV(PTV,proba,date,chaine):
    #Initialisation des Variables
    verbose = False
    index_CP = 0
    try:
        index_PTV = PTV.index[(PTV['debut'] <= 3*60) & (PTV['debut']+PTV['DUREE'] > 3*60)].tolist()[0]
    except Exception as e:
        def_context.Report("can't find a programe to make the link between the two days")
        index_PTV = PTV.shape[0]-1
    def_context.Report('Starting with: %s'%(PTV['TITRE'].iloc[index_PTV]))
    ##########################
    Predictiontimer = 200
    Pubinhour = 0
    lastCP = 3*60-PTV['debut'].loc[index_PTV]
    lastPub= 500
    lastend = PTV['debut'].loc[index_PTV]
    currentduree = PTV['DUREE'].loc[index_PTV]
    planifiedend = PTV['debut'].loc[index_PTV]+PTV['DUREE'].loc[index_PTV]
    begin = True
    nbpub = 0
    Recall = 1
    importantpts = def_context.get_important_points(chaine,PTV,index_PTV)
    help = def_context.get_help(chaine,PTV)
    index_ipts = 0
    error = 0
    chaine = chaine
    ######################
    newPTV = def_context.init_newPTV(PTV,chaine)
    historyofpoints = def_context.init_history(chaine,PTV,lastend,currentduree)
    labels = [0]
    per = 0
    ######################
    for i in tqdm(range(183,1620)):
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

        ###### Let's verify that 'M6 boutique',the algo is not doing a crappy predicitions and if this the case, clean his historic #####
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
                        error+=1*10**index_ipts
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
                            error+=1*10**index_ipts
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
                            error+=1*10**index_ipts
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
                        error+=1*10**index_ipts
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
                            error+=1*10**index_ipts
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
                            error+=1*10**index_ipts
                            #we can now keep going throw the process like before
                            #we just add a line to the history to say that a reset occured
                            newPTV.loc[newPTV.shape[0]] = [i%1440,PTV['TITRE'].iloc[index_PTV],'non',context[3],"--HARD RESET OF ALGORITHM--(out of programme)"]
                            index_ipts+=1
                        else:
                            index_ipts+=1

        elif(context[2]):
            historyofpoints.loc[historyofpoints.shape[0]] = context
            if(lastCP < min(max(currentduree/2,4),5)):
                labels.append(0)
                index_CP+=1
                continue
            elif(lastPub<=6):
                labels.append(0)
                index_CP+=1
                continue
            #Change Point ==> Decide what to do with it
            if(nbpub>=context[7] or Pubinhour >= 12):
                # La pub n'est plus possible dans le programme ==> Soit il s'agit de la fin dans le programme, Soit c'est un faux Change Points
                if(PTV['TITRE'].iloc[index_PTV] in ['M6 boutique','M6 Music'] and context[3]<0.98):
                    labels.append(0)
                elif(PTV['TITRE'].iloc[index_PTV] in ['Les reines du shopping','Absolument stars'] and context[3]<0.96):
                    labels.append(0)
                elif(PTV['TITRE'].iloc[index_PTV] in ['Desperate Housewives','Les Sisters'] and context[3]<0.65 and i<12*60):
                    labels.append(0)
                elif(PTV['TITRE'].iloc[index_PTV] in ['Desperate Housewives'] and context[3]<0.81 and i<12*60+30):
                    labels.append(0)
                elif(historyofpoints['programme'].loc[(historyofpoints.shape[0]-1)] == "Journal"):
                    if(context[3]>0.70 and i<15*60):
                        newPTV.loc[newPTV.shape[0]] = [i%1440,PTV['TITRE'].iloc[index_PTV],'oui',context[3],"fin d'un programme"]
                        lastend = i
                        lastCP=0
                        index_PTV += 1
                        index_PTV = index_PTV%(PTV.shape[0])
                        currentduree = PTV['DUREE'].iloc[index_PTV]
                        planifiedend = (lastend + currentduree)
                        Predictiontimer = 200
                        nbpub = 0
                        per = context[3]
                        labels.append(2)
                    elif(context[3]>0.90 and i>15*60):
                        newPTV.loc[newPTV.shape[0]] = [i%1440,PTV['TITRE'].iloc[index_PTV],'oui',context[3],"fin d'un programme"]
                        lastend = i
                        lastCP=0
                        index_PTV += 1
                        index_PTV = index_PTV%(PTV.shape[0])
                        currentduree = PTV['DUREE'].iloc[index_PTV]
                        planifiedend = (lastend + currentduree)
                        Predictiontimer = 200
                        nbpub = 0
                        per = context[3]
                        labels.append(2)
                    else:
                        labels.append(0)
                elif(historyofpoints['programme'].loc[(historyofpoints.shape[0]-1)] == "film"):
                    if(context[3]>0.75 and context[13]>=0.6):
                        newPTV.loc[newPTV.shape[0]] = [i%1440,PTV['TITRE'].iloc[index_PTV],'oui',context[3],"fin d'un programme"]
                        lastend = i
                        lastCP=0
                        index_PTV += 1
                        index_PTV = index_PTV%(PTV.shape[0])
                        currentduree = PTV['DUREE'].iloc[index_PTV]
                        planifiedend = (lastend + currentduree)
                        Predictiontimer = 200
                        nbpub = 0
                        per = context[3]
                        labels.append(2)
                    else:
                        labels.append(0)
                elif(historyofpoints['duree'].loc[(historyofpoints.shape[0]-1)] == "très court"):
                    newPTV.loc[newPTV.shape[0]] = [i%1440,PTV['TITRE'].iloc[index_PTV],'oui',context[3],"fin d'un programme"]
                    lastend = i
                    lastCP=0
                    index_PTV += 1
                    index_PTV = index_PTV%(PTV.shape[0])
                    currentduree = PTV['DUREE'].iloc[index_PTV]
                    planifiedend = (lastend + currentduree)
                    Predictiontimer = 200
                    nbpub = 0
                    per = context[3]
                    labels.append(2)

                elif(historyofpoints['duree'].loc[(historyofpoints.shape[0]-1)] == "court"):
                    if(context[3]>0.5):
                        newPTV.loc[newPTV.shape[0]] = [i%1440,PTV['TITRE'].iloc[index_PTV],'oui',context[3],"fin d'un programme"]
                        lastend = i
                        lastCP=0
                        index_PTV += 1
                        index_PTV = index_PTV%(PTV.shape[0])
                        currentduree = PTV['DUREE'].iloc[index_PTV]
                        planifiedend = (lastend + currentduree)
                        Predictiontimer = 200
                        nbpub = 0
                        per = context[3]
                        labels.append(2)
                    else:
                        labels.append(0)



                elif(historyofpoints['duree'].loc[(historyofpoints.shape[0]-1)] == "moyen"):
                    if(context[3]>0.80 ):
                        newPTV.loc[newPTV.shape[0]] = [i%1440,PTV['TITRE'].iloc[index_PTV],'oui',context[3],"fin d'un programme"]
                        lastend = i
                        lastCP=0
                        index_PTV += 1
                        index_PTV = index_PTV%(PTV.shape[0])
                        currentduree = PTV['DUREE'].iloc[index_PTV]
                        planifiedend = (lastend + currentduree)
                        Predictiontimer = 200
                        nbpub = 0
                        per = context[3]
                        labels.append(2)
                    elif(context[3]>0.5 and context[5] != 'magazine'):
                        newPTV.loc[newPTV.shape[0]] = [i%1440,PTV['TITRE'].iloc[index_PTV],'oui',context[3],"fin d'un programme"]
                        lastend = i
                        lastCP=0
                        index_PTV += 1
                        index_PTV = index_PTV%(PTV.shape[0])
                        currentduree = PTV['DUREE'].iloc[index_PTV]
                        planifiedend = (lastend + currentduree)
                        Predictiontimer = 200
                        nbpub = 0
                        per = context[3]
                        labels.append(2)
                    else:
                        labels.append(0)

                elif(historyofpoints['duree'].loc[(historyofpoints.shape[0]-1)] == "long"):
                    if(context[3]>0.84):
                        newPTV.loc[newPTV.shape[0]] = [i%1440,PTV['TITRE'].iloc[index_PTV],'oui',context[3],"fin d'un programme"]
                        lastend = i
                        lastCP=0
                        index_PTV += 1
                        index_PTV = index_PTV%(PTV.shape[0])
                        currentduree = PTV['DUREE'].iloc[index_PTV]
                        planifiedend = (lastend + currentduree)
                        Predictiontimer = 200
                        nbpub = 0
                        per = context[3]
                        labels.append(2)
                    else:
                        labels.append(0)

                elif(historyofpoints['duree'].loc[(historyofpoints.shape[0]-1)] == "super long"):
                    if(context[3]>0.90):
                        newPTV.loc[newPTV.shape[0]] = [i%1440,PTV['TITRE'].iloc[index_PTV],'oui',context[3],"fin d'un programme"]
                        lastend = i
                        lastCP=0
                        index_PTV += 1
                        index_PTV = index_PTV%(PTV.shape[0])
                        currentduree = PTV['DUREE'].iloc[index_PTV]
                        planifiedend = (lastend + currentduree)
                        Predictiontimer = 200
                        nbpub = 0
                        per = context[3]
                        labels.append(2)
                    else:
                        labels.append(0)


                elif(historyofpoints['duree'].loc[(historyofpoints.shape[0]-1)] == "très long"):
                    if(context[3]>0.80):
                        newPTV.loc[newPTV.shape[0]] = [i%1440,PTV['TITRE'].iloc[index_PTV],'oui',context[3],"fin d'un programme"]
                        lastend = i
                        lastCP=0
                        index_PTV += 1
                        index_PTV = index_PTV%(PTV.shape[0])
                        currentduree = PTV['DUREE'].iloc[index_PTV]
                        planifiedend = (lastend + currentduree)
                        Predictiontimer = 200
                        nbpub = 0
                        per = context[3]
                        labels.append(2)
                    else:
                        labels.append(0)

            else:
                #la pub est encore possible dans le programme mais pas certaine
                if(context[3]<0.95 and lastPub>=20 and (i-lastend)>15 and  context[5] in ["Feuilleton","Série"]):
                    newPTV.loc[newPTV.shape[0]] = [i%1440,"publicité",'oui',context[3],"publicité dans un programme"]
                    lastCP=0
                    lastPub = 0
                    Pubinhour+=4
                    nbpub+=1
                    labels.append(1)
                elif(context[3]<0.88 and lastPub>=20 and (i-lastend)>=10 and  context[5] in ["Feuilleton","Série"] and context[13]>0.6):
                    newPTV.loc[newPTV.shape[0]] = [i%1440,"publicité",'oui',context[3],"publicité dans un programme"]
                    lastCP=0
                    lastPub = 0
                    Pubinhour+=4
                    nbpub+=1
                    labels.append(1)
                elif(context[3]<0.95 and lastPub>=20 and (i-lastend)>=15 and  context[5] == "magazine"):
                    newPTV.loc[newPTV.shape[0]] = [i%1440,"publicité",'oui',context[3],"publicité dans un programme"]
                    lastCP=0
                    lastPub = 0
                    Pubinhour+=4
                    nbpub+=1
                    labels.append(1)
                elif(context[3]<0.90 and lastPub>=20 and (i-lastend)>=15 and historyofpoints['duree'].loc[(historyofpoints.shape[0]-1)] not in ["super long"] and context[5] not in ["film","Série"]):
                    newPTV.loc[newPTV.shape[0]] = [i%1440,"publicité",'oui',context[3],"publicité dans un programme"]
                    lastCP=0
                    lastPub = 0
                    Pubinhour+=4
                    nbpub+=1
                    labels.append(1)
                elif(context[3]<=0.95 and lastPub>22 and (i-lastend)>=15 and context[5] == "film"):
                    newPTV.loc[newPTV.shape[0]] = [i%1440,"publicité",'oui',context[3],"publicité dans un programme"]
                    lastCP=0
                    lastPub = 0
                    Pubinhour+=4
                    nbpub+=1
                    labels.append(1)
                elif(context[3]<=0.95 and lastPub>=20 and (i-lastend)>=15 and(historyofpoints['programme'].loc[(historyofpoints.shape[0]-1)] in ["Téléréalité","Feuilleton"]) ):
                    newPTV.loc[newPTV.shape[0]] = [i%1440,"publicité",'oui',context[3],"publicité dans un programme"]
                    lastCP=0
                    lastPub = 0
                    Pubinhour+=4
                    nbpub+=1
                    labels.append(1)

                elif(historyofpoints['programme'].loc[(historyofpoints.shape[0]-1)] == 'magazine' and context[3]<=0.95):
                    labels.append(0)

                elif(context[3]<=0.95 and lastPub>=20 and (i-lastend)>=15 and historyofpoints['duree'].loc[(historyofpoints.shape[0]-1)] == "très long" ):
                    newPTV.loc[newPTV.shape[0]] = [i%1440,"publicité",'oui',context[3],"publicité dans un programme"]
                    lastCP=0
                    lastPub = 0
                    Pubinhour+=4
                    nbpub+=1
                    labels.append(1)
                elif(context[3]<=0.95 and lastPub>=20 and (i-lastend)>=15 and historyofpoints['duree'].loc[(historyofpoints.shape[0]-1)] == "super long" ):
                    newPTV.loc[newPTV.shape[0]] = [i%1440,"publicité",'oui',context[3],"publicité dans un programme"]
                    lastCP=0
                    lastPub = 0
                    Pubinhour+=4
                    nbpub+=1
                    labels.append(1)

                elif(historyofpoints['duree'].loc[(historyofpoints.shape[0]-1)] == "très court"):
                    newPTV.loc[newPTV.shape[0]] = [i%1440,PTV['TITRE'].iloc[index_PTV],'oui',context[3],"fin d'un programme"]
                    lastend = i
                    lastCP=0
                    index_PTV += 1
                    index_PTV = index_PTV%(PTV.shape[0])
                    currentduree = PTV['DUREE'].iloc[index_PTV]
                    planifiedend = (lastend + currentduree)
                    Predictiontimer = 200
                    nbpub = 0
                    per = context[3]
                    labels.append(2)

                elif(historyofpoints['duree'].loc[(historyofpoints.shape[0]-1)] == "court"):
                    if(context[3]>0.85):
                        newPTV.loc[newPTV.shape[0]] = [i%1440,PTV['TITRE'].iloc[index_PTV],'oui',context[3],"fin d'un programme"]
                        lastend = i
                        lastCP=0
                        index_PTV += 1
                        index_PTV = index_PTV%(PTV.shape[0])
                        currentduree = PTV['DUREE'].iloc[index_PTV]
                        planifiedend = (lastend + currentduree)
                        Predictiontimer = 200
                        nbpub = 0
                        per = context[3]
                        labels.append(2)
                    else:
                        labels.append(0)



                elif(historyofpoints['duree'].loc[(historyofpoints.shape[0]-1)] == "moyen"):
                    if(context[3]>0.85):
                        newPTV.loc[newPTV.shape[0]] = [i%1440,PTV['TITRE'].iloc[index_PTV],'oui',context[3],"fin d'un programme"]
                        lastend = i
                        lastCP=0
                        index_PTV += 1
                        index_PTV = index_PTV%(PTV.shape[0])
                        currentduree = PTV['DUREE'].iloc[index_PTV]
                        planifiedend = (lastend + currentduree)
                        Predictiontimer = 200
                        nbpub = 0
                        per = context[3]
                        labels.append(2)
                    else:
                        labels.append(0)

                elif(historyofpoints['duree'].loc[(historyofpoints.shape[0]-1)] == "long"):
                    if(context[3]>0.90):
                        newPTV.loc[newPTV.shape[0]] = [i%1440,PTV['TITRE'].iloc[index_PTV],'oui',context[3],"fin d'un programme"]
                        lastend = i
                        lastCP=0
                        index_PTV += 1
                        index_PTV = index_PTV%(PTV.shape[0])
                        currentduree = PTV['DUREE'].iloc[index_PTV]
                        planifiedend = (lastend + currentduree)
                        Predictiontimer = 200
                        nbpub = 0
                        per = context[3]
                        labels.append(2)
                    else:
                        labels.append(0)
                elif(historyofpoints['duree'].loc[(historyofpoints.shape[0]-1)] == "très long"):
                    if(context[3]>0.90):
                        newPTV.loc[newPTV.shape[0]] = [i%1440,PTV['TITRE'].iloc[index_PTV],'oui',context[3],"fin d'un programme"]
                        lastend = i
                        lastCP=0
                        index_PTV += 1
                        index_PTV = index_PTV%(PTV.shape[0])
                        currentduree = PTV['DUREE'].iloc[index_PTV]
                        planifiedend = (lastend + currentduree)
                        Predictiontimer = 200
                        nbpub = 0
                        per = context[3]
                        labels.append(2)
                    else:
                        labels.append(0)
                elif(historyofpoints['duree'].loc[(historyofpoints.shape[0]-1)] == "super long"):
                    if(context[3]>0.99):
                        newPTV.loc[newPTV.shape[0]] = [i%1440,PTV['TITRE'].iloc[index_PTV],'oui',context[3],"fin d'un programme"]
                        lastend = i
                        lastCP=0
                        index_PTV += 1
                        index_PTV = index_PTV%(PTV.shape[0])
                        currentduree = PTV['DUREE'].iloc[index_PTV]
                        planifiedend = (lastend + currentduree)
                        Predictiontimer = 200
                        nbpub = 0
                        per = context[3]
                        labels.append(2)
                    else:
                        labels.append(0)
                else:
                    labels.append(0)
            index_CP+=1
        elif(i in help):
            historyofpoints.loc[historyofpoints.shape[0]] = context
            if(lastCP < min(max(currentduree/2,4),5)):
                labels.append(0)
                continue
            #Change Point ==> Decide what to do with it
            if(nbpub>=context[7] or Pubinhour >= 12):
                if(PTV['TITRE'].iloc[index_PTV] in ['M6 boutique','M6 Music'] and context[3]<0.98):
                    labels.append(0)
                elif(PTV['TITRE'].iloc[index_PTV] in ['Les reines du shopping','Absolument stars','Les Sisters'] and context[3]<0.95):
                    labels.append(0)
                elif(PTV['TITRE'].iloc[index_PTV] in ['Desperate Housewives'] and context[3]<0.65 and i<12*60):
                    labels.append(0)
                elif(PTV['TITRE'].iloc[index_PTV] in ['Desperate Housewives'] and context[3]<0.81 and i<12*60+30):
                    labels.append(0)
                # La pub n'est plus possible dans le programme ==> Soit il s'agit de la fin dans le programme, Soit c'est un faux Change Points
                elif(historyofpoints['duree'].loc[(historyofpoints.shape[0]-1)] == "très court"):
                    newPTV.loc[newPTV.shape[0]] = [i%1440,PTV['TITRE'].iloc[index_PTV],'oui',context[3],"fin d'un programme"]
                    lastend = i
                    lastCP=0
                    index_PTV += 1
                    index_PTV = index_PTV%(PTV.shape[0])
                    currentduree = PTV['DUREE'].iloc[index_PTV]
                    planifiedend = (lastend + currentduree)
                    Predictiontimer = 200
                    nbpub = 0
                    per = context[3]
                    labels.append(2)

                elif(historyofpoints['duree'].loc[(historyofpoints.shape[0]-1)] == "court"):
                    if(context[3]>0.5):
                        newPTV.loc[newPTV.shape[0]] = [i%1440,PTV['TITRE'].iloc[index_PTV],'oui',context[3],"fin d'un programme"]
                        lastend = i
                        lastCP=0
                        index_PTV += 1
                        index_PTV = index_PTV%(PTV.shape[0])
                        currentduree = PTV['DUREE'].iloc[index_PTV]
                        planifiedend = (lastend + currentduree)
                        Predictiontimer = 200
                        nbpub = 0
                        per = context[3]
                        labels.append(2)
                    else:
                        labels.append(0)



                elif(historyofpoints['duree'].loc[(historyofpoints.shape[0]-1)] == "moyen"):
                    if(context[3]>0.80 ):
                        newPTV.loc[newPTV.shape[0]] = [i%1440,PTV['TITRE'].iloc[index_PTV],'oui',context[3],"fin d'un programme"]
                        lastend = i
                        lastCP=0
                        index_PTV += 1
                        index_PTV = index_PTV%(PTV.shape[0])
                        currentduree = PTV['DUREE'].iloc[index_PTV]
                        planifiedend = (lastend + currentduree)
                        Predictiontimer = 200
                        nbpub = 0
                        per = context[3]
                        labels.append(2)
                    elif(context[3]>0.8 and context[5] != 'magazine'):
                        newPTV.loc[newPTV.shape[0]] = [i%1440,PTV['TITRE'].iloc[index_PTV],'oui',context[3],"fin d'un programme"]
                        lastend = i
                        lastCP=0
                        index_PTV += 1
                        index_PTV = index_PTV%(PTV.shape[0])
                        currentduree = PTV['DUREE'].iloc[index_PTV]
                        planifiedend = (lastend + currentduree)
                        Predictiontimer = 200
                        nbpub = 0
                        per = context[3]
                        labels.append(2)
                    else:
                        labels.append(0)

                elif(historyofpoints['duree'].loc[(historyofpoints.shape[0]-1)] == "long"):
                    if(context[3]>=0.85):
                        newPTV.loc[newPTV.shape[0]] = [i%1440,PTV['TITRE'].iloc[index_PTV],'oui',context[3],"fin d'un programme"]
                        lastend = i
                        lastCP=0
                        index_PTV += 1
                        index_PTV = index_PTV%(PTV.shape[0])
                        currentduree = PTV['DUREE'].iloc[index_PTV]
                        planifiedend = (lastend + currentduree)
                        Predictiontimer = 200
                        nbpub = 0
                        per = context[3]
                        labels.append(2)
                    else:
                        labels.append(0)

                elif(historyofpoints['duree'].loc[(historyofpoints.shape[0]-1)] == "super long"):
                    if(context[3]>0.90):
                        newPTV.loc[newPTV.shape[0]] = [i%1440,PTV['TITRE'].iloc[index_PTV],'oui',context[3],"fin d'un programme"]
                        lastend = i
                        lastCP=0
                        index_PTV += 1
                        index_PTV = index_PTV%(PTV.shape[0])
                        currentduree = PTV['DUREE'].iloc[index_PTV]
                        planifiedend = (lastend + currentduree)
                        Predictiontimer = 200
                        nbpub = 0
                        per = context[3]
                        labels.append(2)
                    else:
                        labels.append(0)


                elif(historyofpoints['duree'].loc[(historyofpoints.shape[0]-1)] == "très long"):
                    if(context[3]>0.80):
                        newPTV.loc[newPTV.shape[0]] = [i%1440,PTV['TITRE'].iloc[index_PTV],'oui',context[3],"fin d'un programme"]
                        lastend = i
                        lastCP=0
                        index_PTV += 1
                        index_PTV = index_PTV%(PTV.shape[0])
                        currentduree = PTV['DUREE'].iloc[index_PTV]
                        planifiedend = (lastend + currentduree)
                        Predictiontimer = 200
                        nbpub = 0
                        per = context[3]
                        labels.append(2)
                    else:
                        labels.append(0)

            else:
                #la pub est encore possible dans le programme mais pas certaine
                if(context[3]<=0.95 and lastPub>=20 and (i-lastend)>=max(10,0.2*currentduree) and historyofpoints['duree'].loc[(historyofpoints.shape[0]-1)] not in ["très long","super long"] and context[5] != "film"):
                    newPTV.loc[newPTV.shape[0]] = [i%1440,"publicité",'oui',context[3],"publicité dans un programme"]
                    lastCP=0
                    lastPub = 0
                    Pubinhour+=4
                    nbpub+=1
                    labels.append(1)

                elif(context[3]<=0.95 and lastPub>=20 and (i-lastend)>=20 and context[5] == "film"):
                    newPTV.loc[newPTV.shape[0]] = [i%1440,"publicité",'oui',context[3],"publicité dans un programme"]
                    lastCP=0
                    lastPub = 0
                    Pubinhour+=4
                    nbpub+=1
                    labels.append(1)
                elif(context[3]<=0.95  and lastPub>=20 and (historyofpoints['programme'].loc[(historyofpoints.shape[0]-1)] in ["Téléréalité","Feuilleton"]) ):
                    newPTV.loc[newPTV.shape[0]] = [i%1440,"publicité",'oui',context[3],"publicité dans un programme"]
                    lastCP=0
                    lastPub = 0
                    Pubinhour+=4
                    nbpub+=1
                    labels.append(1)
                elif(context[3]<=0.95 and lastPub>=20 and (i-lastend)>=10 and historyofpoints['duree'].loc[(historyofpoints.shape[0]-1)] == "très long" ):
                    newPTV.loc[newPTV.shape[0]] = [i%1440,"publicité",'oui',context[3],"publicité dans un programme"]
                    lastCP=0
                    lastPub = 0
                    Pubinhour+=4
                    nbpub+=1
                    labels.append(1)
                elif(historyofpoints['duree'].loc[(historyofpoints.shape[0]-1)] == "très court"):
                    newPTV.loc[newPTV.shape[0]] = [i%1440,PTV['TITRE'].iloc[index_PTV],'oui',context[3],"fin d'un programme"]
                    lastend = i
                    lastCP=0
                    index_PTV += 1
                    index_PTV = index_PTV%(PTV.shape[0])
                    currentduree = PTV['DUREE'].iloc[index_PTV]
                    planifiedend = (lastend + currentduree)
                    Predictiontimer = 200
                    nbpub = 0
                    per = context[3]
                    labels.append(2)

                elif(historyofpoints['duree'].loc[(historyofpoints.shape[0]-1)] == "court"):
                    if(context[3]>0.5):
                        newPTV.loc[newPTV.shape[0]] = [i%1440,PTV['TITRE'].iloc[index_PTV],'oui',context[3],"fin d'un programme"]
                        lastend = i
                        lastCP=0
                        index_PTV += 1
                        index_PTV = index_PTV%(PTV.shape[0])
                        currentduree = PTV['DUREE'].iloc[index_PTV]
                        planifiedend = (lastend + currentduree)
                        Predictiontimer = 200
                        nbpub = 0
                        per = context[3]
                        labels.append(2)
                    else:
                        labels.append(0)



                elif(historyofpoints['duree'].loc[(historyofpoints.shape[0]-1)] == "moyen"):
                    if(context[3]>0.5):
                        newPTV.loc[newPTV.shape[0]] = [i%1440,PTV['TITRE'].iloc[index_PTV],'oui',context[3],"fin d'un programme"]
                        lastend = i
                        lastCP=0
                        index_PTV += 1
                        index_PTV = index_PTV%(PTV.shape[0])
                        currentduree = PTV['DUREE'].iloc[index_PTV]
                        planifiedend = (lastend + currentduree)
                        Predictiontimer = 200
                        nbpub = 0
                        per = context[3]
                        labels.append(2)
                    else:
                        labels.append(0)

                elif(historyofpoints['duree'].loc[(historyofpoints.shape[0]-1)] == "long"):
                    if(context[3]>0.75):
                        newPTV.loc[newPTV.shape[0]] = [i%1440,PTV['TITRE'].iloc[index_PTV],'oui',context[3],"fin d'un programme"]
                        lastend = i
                        lastCP=0
                        index_PTV += 1
                        index_PTV = index_PTV%(PTV.shape[0])
                        currentduree = PTV['DUREE'].iloc[index_PTV]
                        planifiedend = (lastend + currentduree)
                        Predictiontimer = 200
                        nbpub = 0
                        per = context[3]
                        labels.append(2)
                    else:
                        labels.append(0)
                elif(historyofpoints['duree'].loc[(historyofpoints.shape[0]-1)] == "très long"):
                    if(context[3]>0.91):
                        newPTV.loc[newPTV.shape[0]] = [i%1440,PTV['TITRE'].iloc[index_PTV],'oui',context[3],"fin d'un programme"]
                        lastend = i
                        lastCP=0
                        index_PTV += 1
                        index_PTV = index_PTV%(PTV.shape[0])
                        currentduree = PTV['DUREE'].iloc[index_PTV]
                        planifiedend = (lastend + currentduree)
                        Predictiontimer = 200
                        nbpub = 0
                        per = context[3]
                        labels.append(2)
                    else:
                        labels.append(0)
                elif(historyofpoints['duree'].loc[(historyofpoints.shape[0]-1)] == "super long"):
                    if(context[3]>0.90):
                        newPTV.loc[newPTV.shape[0]] = [i%1440,PTV['TITRE'].iloc[index_PTV],'oui',context[3],"fin d'un programme"]
                        lastend = i
                        lastCP=0
                        index_PTV += 1
                        index_PTV = index_PTV%(PTV.shape[0])
                        currentduree = PTV['DUREE'].iloc[index_PTV]
                        planifiedend = (lastend + currentduree)
                        Predictiontimer = 200
                        nbpub = 0
                        per = context[3]
                        labels.append(2)
                    else:
                        labels.append(0)
                else:
                    labels.append(0)


        else:
            #labels.append(0)
            #Not a Change Point, we'll just check that nothing is wrong in the PTV at this time
            if(Predictiontimer == 0):
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
            elif(context[3] > 1 and Predictiontimer > 20):
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


    return newPTV,historyofpoints,labels,error




#################################################
########### main with options ###################
#################################################



def main(argv):
    if(len(argv) == 0):
        def_context.Report("Pas suffisement d'arguments pour le main de PTVM6.py")
    else:
        global PATH_IN,PATH_SCRIPT,PATH_OUT
        PATH_IN,PATH_SCRIPT,PATH_OUT = def_context.get_path()
        date = argv[0]
        def_context.Report(date)
        d = "".join(date.split('-'))
        if(len(argv) == 1):
            chaine = 'M6'
        else:
            chaine = argv[1]
        number,chaine = def_context.get_tuple(chaine)
        PTV,proba = def_context.load_file(date,chaine)
        if(len(PTV) == 0):
            sys.exit(4)
            return 0
        new_PTV,historyofpoints,labels,error = make_newPTV(PTV,proba,date,chaine)

        new_PTV['Heure'] = new_PTV['minute'].apply(lambda x: str(int(x/60))+':'+str(x%60))
        historyofpoints['Heure'] = historyofpoints['minute'].apply(lambda x: str(int(x/60))+':'+str(x%60))
        new_PTV.to_html(PATH_IN+'new_ptv/new_PTV_'+date+'_'+chaine+'.html')
        new_PTV.to_csv(PATH_IN+'new_ptv/new_PTV_'+date+'_'+chaine+'.csv',index=False)
        historyofpoints['labels'] = labels
        historyofpoints.to_html(PATH_IN+'hop/historyofpoints_'+date+'_'+chaine+'.html')
        historyofpoints.to_csv(PATH_IN+'hop/historyofpoints_'+date+'_'+chaine+'.csv',index=False)
        def_context.Report("%s,%s,%s"%(date,historyofpoints.shape,len(labels)))
        new_PTV.to_html(PATH_OUT+'T0/new_ptv/new_PTV_'+date+'_'+chaine+'.html')
        new_PTV.to_csv(PATH_OUT+'T0/new_ptv/new_PTV_'+date+'_'+chaine+'.csv')
        sys.exit(error)



if __name__ == "__main__":
    # execute only if run as a script
    main(sys.argv[1:])
