#################################################
#created the 27/07/2018 16:52 by Alexis Blanchet#
#################################################
#-*- coding: utf-8 -*-
'''
Objectif:
Contrôler l'ensemble de la conception des programmes TV
Prends entrée les données nétoyées et traitées et après divers interactions
renvoie le nouveau programme TV sur un mois pour toutes les chaînes dont on
posséde l'information.

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
import time
import def_context
import pickle
from subprocess import Popen
from multiprocessing import Process
#TODO
#################################################
########### Global variables ####################
#################################################
PATH_IN = '/home/alexis/Bureau/finalproject/Datas/'
PATH_SCRIPT = '/home/alexis/Bureau/finalproject/scripts/'
PATH_OUT = '/home/alexis/Bureau/finalproject/Datas/'
LOG = "log.txt"
MAX_PROCESSES = 15
#################################################
########### Important functions #################
#################################################


def get_path():
    datas = pd.read_csv('path.csv')
    return datas['PathtoTempDatas'].values[0],datas['PathtoScripts'].values[0],datas['PathtoDatasOut'].values[0]



def Report(error):
    with open(LOG,'a+') as file:
        file.write(str(error)+' \n')
        def_context.Report(str(error))



def get_temp_path():
    datas = pd.read_csv('path.csv')
    return datas['temp_path'].values[0]



def update_temp_path(i):
    datas = pd.read_csv('path.csv')
    datas['temp_path'] = datas['PathtoDatasOut']+'T'+str(i)+"/"
    def_context.Report('Updated Temp path to: '+datas['PathtoDatasOut'][0]+'T'+str(i)+"/")
    datas.to_csv('path.csv',index=False)

def pred(file):
    try:
        PATH_OUT = get_temp_path()
        relecture = True
        EPSILON = 1e-15
        f = ((file.split('.'))[0].split('_'))[2]
        c = ((file.split('.'))[0].split('_'))[-1]
        PTV,proba = def_context.load_file(str(f),str(c))
        if(len(PTV) == 0):
            return 0
        index_PTV = PTV.index[(PTV['debut'] <= 3*60) & (PTV['debut']+PTV['DUREE'] > 3*60+5)].tolist()[0]
        def_context.Report('Starting with: %s'%(PTV['TITRE'].iloc[index_PTV]))
        lastend = PTV['debut'].loc[index_PTV]
        currentduree = PTV['DUREE'].loc[index_PTV]
        newPTV = def_context.init_newPTV(PTV,str(c))
        historyofpoints = def_context.init_history(str(c),PTV,lastend,currentduree)
        temp_context = historyofpoints.iloc[0]
        importantpts = def_context.get_important_points(c,PTV,index_PTV)
        for i in range(3):
            def_context.Report(str(i)+' '+str(c)+' '+str(f))
            from predictPTV import main as pred
            l1,temp_newPTV1,temp_history1,index_PTV1,temp_context1 = pred([str(c),str(f),i,newPTV.iloc[newPTV.shape[0]-1],temp_context,index_PTV,importantpts])
            if(l1>0 and relecture):
                def_context.Report("Utilisation de la relecture "+str(i)+' '+str(c)+' '+str(f))
                from RLPTV import main as RL
                l2,temp_newPTV2,temp_history2,index_PTV2,temp_context2 = RL([str(c),str(f),i,newPTV.iloc[newPTV.shape[0]-1],temp_context,index_PTV,importantpts])
                if(l2>5):
                    def_context.Report("Utilisation de l'arbre de décision",f,c,i)
                    if(chaine == 'TF1'):
                        from PTVTF1 import main as arbre1
                        l3,temp_newPTV3,temp_history3,index_PTV3,temp_context3 = arbre1([str(c),str(f),i,newPTV.loc[newPTV.shape[0]-1],temp_context,index_PTV,importantpts])
                    elif(chaine == 'M6'):
                        from PTVM6 import main as arbre2
                        l3,temp_newPTV3,temp_history3,index_PTV3,temp_context3 = arbre2([str(c),str(f),i,newPTV.loc[newPTV.shape[0]-1],temp_context,index_PTV,importantpts])
                    else:
                        l3>5
                    if(l3>0):
                        def_context.Report("AUCUNE DÉCISION NE CONVIENT",f,c)
                        l,temp_newPTV,temp_history,index_PTV,temp_context = l2,temp_newPTV2,temp_history2,index_PTV2,temp_context2
                    else:
                        l,temp_newPTV,temp_history,index_PTV,temp_context = l3,temp_newPTV3,temp_history3,index_PTV3,temp_context3
                else:
                    l,temp_newPTV,temp_history,index_PTV,temp_context = l2,temp_newPTV2,temp_history2,index_PTV2,temp_context2
            else:
                l,temp_newPTV,temp_history,index_PTV,temp_context = l1,temp_newPTV1,temp_history1,index_PTV1,temp_context1
            if(l == 4):
                pass
            else:
                newPTV = pd.concat([newPTV,temp_newPTV.iloc[1:]])
                historyofpoints = pd.concat([historyofpoints,temp_history])

        newPTV['Heure'] = newPTV['minute'].apply(lambda x: str(int(x/60))+':'+str(x%60))
        historyofpoints['Heure'] = historyofpoints['minute'].apply(lambda x: str(int(x/60))+':'+str(x%60))
        newPTV.to_html(PATH_IN+'new_ptv/new_PTV_'+str(f)+'_'+str(c)+'.html')
        newPTV.to_csv(PATH_IN+'new_ptv/new_PTV_'+str(f)+'_'+str(c)+'.csv',index=False)
        historyofpoints.to_html(PATH_IN+'hop/historyofpoints_'+str(f)+'_'+str(c)+'.html')
        historyofpoints.to_csv(PATH_IN+'hop/historyofpoints_'+str(f)+'_'+str(c)+'.csv',index=False)
        newPTV.to_html(PATH_OUT+'new_ptv/new_PTV_'+str(f)+'_'+str(c)+'.html')
        newPTV.to_csv(PATH_OUT+'new_ptv/new_PTV_'+str(f)+'_'+str(c)+'.csv',index=False)
    except Exception as e:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        def_context.Report("Failed to process {0} at line {2} in {3}: {1}".format(str(file), str(e),sys.exc_info()[-1].tb_lineno,fname))

########################################
########## Main with some otpions ######
########################################

def main(argv):
    global PATH_IN,PATH_SCRIPT,PATH_OUT
    PATH_IN,PATH_SCRIPT,PATH_OUT = get_path()
    import pandas as pd
    import pickle
    end = 30

    if(len(argv) == 0):
        start = 0
        for i in range(start,end):
            t = time.time()
            update_temp_path(i)
            if(i == 0):
                p1 = Popen(['python',PATH_SCRIPT+'PTVall.py','2017-12'])
                p1.wait()
                time.sleep(60)
                os.system('python '+PATH_SCRIPT+'MLforPTV.py')
                def_context.Report("fin du tour "+str(i))
                def_context.Report(time.time()-t)
            else:
                p1 = Popen(['python',PATH_SCRIPT+'fasttrain.py','2017-12','TF1'])
                p2 = Popen(['python',PATH_SCRIPT+'fasttrain.py','2017-12','M6'])
                p3 = Popen(['python',PATH_SCRIPT+'fasttrain.py','2017-12','France 2'])
                p4 = Popen(['python',PATH_SCRIPT+'fasttrain.py','2017-12','France 3'])
                p1.wait()
                p2.wait()
                p3.wait()
                p4.wait()
                print('end of prediction')
                time.sleep(60)
                os.system('python '+PATH_SCRIPT+'MLforPTV.py')
                def_context.Report("fin du tour "+str(i))
                def_context.Report(time.time()-t)
                def_context.Report(time.time()-t)
    if(len(argv) == 1):
        try:
            start = int(argv[0])
            for i in range(start,end):
                update_temp_path(i)
                if(i == 0):
                    p1 = Popen(['python',PATH_SCRIPT+'PTVall.py','2017-12'])
                    p1.wait()
                    time.sleep(60)
                    os.system('python '+PATH_SCRIPT+'MLforPTV.py')
                    time.sleep(60)
                    def_context.Report("fin du tour "+str(i))
                else:
                    p1 = Popen(['python',PATH_SCRIPT+'fasttrain.py','2017-12','TF1'])
                    p2 = Popen(['python',PATH_SCRIPT+'fasttrain.py','2017-12','M6'])
                    p3 = Popen(['python',PATH_SCRIPT+'fasttrain.py','2017-12','France 2'])
                    p4 = Popen(['python',PATH_SCRIPT+'fasttrain.py','2017-12','France 3'])
                    p1.wait()
                    p2.wait()
                    p3.wait()
                    p4.wait()
                    print('end of prediction')
                    time.sleep(60)
                    os.system('python '+PATH_SCRIPT+'MLforPTV.py')
                    time.sleep(60)
                    def_context.Report("fin du tour "+str(i))
        except Exception as e:
            pred(argv[0])
    elif(len(argv) == 2):
        PATH_OUT = get_temp_path()
        import pandas as pd
        import random
        relecture = True
        EPSILON = 1e-15
        files = os.listdir(PATH_IN+'PTV/')
        nb_files = len(files)
        Processes = []
        for file in files:
            f = ((file.split('.'))[0].split('_'))[2]
            c = ((file.split('.'))[0].split('_'))[-1]
            if(f=='2017-12-20' or (f in ['2017-12-09','2017-12-06','2018-02-22'] and c=='TF1') or (f in ['2018-02-22'] and c=='M6') or (f.split('-')[0] != str(argv[0].split('-')[0])) or f.split('-')[1] != argv[0].split('-')[1]):
                #or (f in ['2018-02-22'] and c=='M6')
                pass
            elif(c ==''):
                pass
            elif('2018' in f):
                pass
            elif(c == argv[1]):
                def_context.Report(file)
                while(len(Processes)>= MAX_PROCESSES):
                    lenp = len(Processes)
                    for p in range(lenp):  # Check the processes in reverse order
                        if Processes[enp - 1 - p].poll() is not None:  # If the process hasn't finished will return None
                            del Processes[lenp - 1 - p]  # Remove from list - this is why we needed reverse order
                    time.sleep(5)
                Processes.append(Popen(['python',PATH_SCRIPT+'fasttrain.py',file]))
            else:
                pass
        while(len(Processes)):
            lenp = len(Processes)
            for p in range(lenp):  # Check the processes in reverse order
                if Processes[enp - 1 - p].poll() is not None:  # If the process hasn't finished will return None
                    del Processes[lenp - 1 - p]  # Remove from list - this is why we needed reverse order
            time.sleep(5)


if __name__ == "__main__":
    # execute only if run as a script
    main(sys.argv[1:])
