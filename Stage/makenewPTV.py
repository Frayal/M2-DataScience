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

########################################
########## Main with some otpions ######
########################################

def main(argv):
    global PATH_IN,PATH_SCRIPT,PATH_OUT
    PATH_IN,PATH_SCRIPT,PATH_OUT = get_path()
    import pandas as pd
    import pickle
    createfile = False
    end = 30
    t = time.time()
    if(len(argv) ==0):
        argv = ['2015']
    if(argv[0] == 'start'):
        if(len(argv) == 1):
            start = 0
            createfile = True
        else:
            start = int(argv[1])
            if(start == 0):
                createfile = True
        if(createfile):
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
            time.sleep(10)
        for i in range(start,end):
            update_temp_path(i)
            try:
                open(PATH_OUT+'res.txt', 'w').close()
                def_context.Report('file cleaned')
            except Exception as e:
                pass
            if(createfile and i == 0):
                p1 = Popen(['python',PATH_SCRIPT+'PTVall.py','2017-12'])
                p2 = Popen(['python',PATH_SCRIPT+'PTVall.py','2018-02'])
                p3 = Popen(['python',PATH_SCRIPT+'PTVall.py','2018-03'])
                p1.wait()
                p2.wait()
                p3.wait()
                """
                os.system('python '+PATH_SCRIPT+'PTVall.py 2017-12')
                os.system('python '+PATH_SCRIPT+'PTVall.py 2018-02')
                os.system('python '+PATH_SCRIPT+'PTVall.py 2018-03')
                """
                time.sleep(60)
                os.system('python '+PATH_SCRIPT+'MLforPTV.py')
                time.sleep(60)
                def_context.Report("fin du tour "+str(i))
            else:
                p1 = Popen(['python',PATH_SCRIPT+'makenewPTV.py','2017-12'])
                p2 = Popen(['python',PATH_SCRIPT+'makenewPTV.py','2018-02'])
                p3 = Popen(['python',PATH_SCRIPT+'makenewPTV.py','2018-03'])
                p1.wait()
                p2.wait()
                p3.wait()
                """
                os.system('python '+PATH_SCRIPT+'makenewPTV.py 2017-12')
                os.system('python '+PATH_SCRIPT+'makenewPTV.py 2018-02')
                os.system('python '+PATH_SCRIPT+'makenewPTV.py 2018-03')
                """
                time.sleep(60)
                os.system('python '+PATH_SCRIPT+'MLforPTV.py')
                time.sleep(60)
                def_context.Report("fin du tour "+str(i))


    elif(len(argv) == 1 and argv[0] != 'start'):
        PATH_OUT = get_temp_path()
        import pandas as pd
        import random
        relecture = True
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
        nb_rel= 0

        files = os.listdir(PATH_IN+'PTV/')
        nb_files = len(files)
        for file in files:
            def_context.Report('Il reste encore %s fichiers à traiter'%(nb_files))
            nb_files -=1
            try:
                f = ((file.split('.'))[0].split('_'))[2]
                c = ((file.split('.'))[0].split('_'))[-1]
                if(f=='2017-12-20' or (f in ['2017-12-09','2017-12-06','2018-02-22'] and c=='TF1') or (f in ['2018-02-22'] and c=='M6') or (f.split('-')[0] != str(argv[0].split('-')[0])) or f.split('-')[1] != argv[0].split('-')[1]):
                    #or (f in ['2018-02-22'] and c=='M6')
                    pass
                elif(c ==''):
                    pass
                else:
                    PTV,proba = def_context.load_file(str(f),str(c))
                    if(len(PTV) == 0):
                        continue
                    index_PTV = PTV.index[(PTV['debut'] <= 3*60) & (PTV['debut']+PTV['DUREE'] > 3*60+5)].tolist()[0]
                    def_context.Report('Starting with: %s'%(PTV['TITRE'].iloc[index_PTV]))
                    lastend = PTV['debut'].loc[index_PTV]
                    currentduree = PTV['DUREE'].loc[index_PTV]
                    newPTV = def_context.init_newPTV(PTV,str(c))
                    historyofpoints = def_context.init_history(str(c),PTV,lastend,currentduree)
                    temp_context = historyofpoints.iloc[0]
                    importantpts = def_context.get_important_points(c,PTV,index_PTV)
                    file_ = open(PATH_OUT+'res.txt', 'a+')
                    file_.write(str(f+' '+c+':').rstrip('\n'))
                    for i in range(3):
                        def_context.Report(str(i)+' '+str(c)+' '+str(f))
                        from predictPTV import main as pred
                        l1,temp_newPTV1,temp_history1,index_PTV1,temp_context1 = pred([str(c),str(f),i,newPTV.iloc[newPTV.shape[0]-1],temp_context,index_PTV,importantpts,PATH_OUT])
                        if(l1>0 and relecture):
                            nb_rel+=1
                            def_context.Report("Utilisation de la relecture "+str(i)+' '+str(c)+' '+str(f))
                            from RLPTV import main as RL
                            l2,temp_newPTV2,temp_history2,index_PTV2,temp_context2 = RL([str(c),str(f),i,newPTV.iloc[newPTV.shape[0]-1],temp_context,index_PTV,importantpts,PATH_OUT])
                            if(l2>5):
                                def_context.Report("Utilisation de l'arbre de décision",f,c,i)
                                if(chaine == 'TF1'):
                                    from PTVTF1 import main as arbre
                                elif(chaine == 'M6'):
                                    from PTVM6 import main as arbre
                                else:
                                    l3>5
                                l3,temp_newPTV3,temp_history3,index_PTV3,temp_context3 = arbre([str(c),str(f),i,newPTV.loc[newPTV.shape[0]-1],temp_context,index_PTV,importantpts])
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
                            err += l
                            if(i == 0):
                                err_type_1 += l
                            if(i == 1):
                                err_type_2 += l
                            if(i == 2):
                                err_type_3 += l
                            m+=1
                            if(c == 'M6'):
                                err_M6 += l
                                m_M6 += 1
                            if(c == 'TF1'):
                                err_TF1 += l
                                m_TF1 += 1
                            if(c == 'France 2'):
                                err_F2 += l
                                m_F2 += 1
                            if(c == 'France 3'):
                                err_F3 += l
                                m_F3 += 1
                            file_.write(str(l).rstrip('\n'))
                            file_.write(" ".rstrip('\n'))

                    newPTV['Heure'] = newPTV['minute'].apply(lambda x: str(int(x/60))+':'+str(x%60))
                    historyofpoints['Heure'] = historyofpoints['minute'].apply(lambda x: str(int(x/60))+':'+str(x%60))
                    newPTV.to_html(PATH_IN+'new_ptv/new_PTV_'+str(f)+'_'+str(c)+'.html')
                    newPTV.to_csv(PATH_IN+'new_ptv/new_PTV_'+str(f)+'_'+str(c)+'.csv',index=False)
                    historyofpoints.to_html(PATH_IN+'hop/historyofpoints_'+str(f)+'_'+str(c)+'.html')
                    historyofpoints.to_csv(PATH_IN+'hop/historyofpoints_'+str(f)+'_'+str(c)+'.csv',index=False)
                    newPTV.to_html(PATH_OUT+'new_ptv/new_PTV_'+str(f)+'_'+str(c)+'.html')
                    newPTV.to_csv(PATH_OUT+'new_ptv/new_PTV_'+str(f)+'_'+str(c)+'.csv',index=False)
                    #newPTV.to_csv(PATH_OUT+'new_ptv/new_PTV_'+str(f)+'_'+str(c)+'.csv',index=False)
                    #historyofpoints.to_html(PATH_OUT+'hop/historyofpoints_'+str(f)+'_'+str(c)+'.html')
                    #historyofpoints.to_csv(PATH_OUT+'hop/historyofpoints_'+str(f)+'_'+str(c)+'.csv',index=False)
                    file_.write("\n")
                    file_.close()

                def_context.Report(err)
            except Exception as e:
                exc_type, exc_obj, exc_tb = sys.exc_info()
                fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                def_context.Report("Failed to process {0} at line {2} in {3}: {1}".format(str(file), str(e),sys.exc_info()[-1].tb_lineno,fname))

        def_context.Report(m)
        def_context.Report("score Total:"+str(1-(err/(m+EPSILON))))
        def_context.Report("score TF1:"+str(1-(err_TF1/(m_TF1+EPSILON))))
        def_context.Report("score M6:"+str(1-(err_M6/(m_M6+EPSILON))))
        def_context.Report("score France 2:"+str(1-(err_F2/(m_F2+EPSILON))))
        def_context.Report("score France 3:"+str(1-(err_F3/(m_F3+EPSILON))))
        def_context.Report("score sur la matinée:"+str(1-((err_type_1*3)/(m+EPSILON))))
        def_context.Report("score sur l'après midi:"+str(1-((err_type_2*3)/(m+EPSILON))))
        def_context.Report("score sur la soirée:"+str(1-((err_type_3*3)/(m+EPSILON))))
        def_context.Report("temps de calcul:"+str(time.time()-t))
        df = pd.read_csv('scores.csv')
        df.loc[df.shape[0]] = [1-(err_TF1/(m_TF1+EPSILON)),1-(err_M6/(m_M6+EPSILON)),1-(err_F2/(m_F2+EPSILON)),1-(err_F3/(m_F3+EPSILON)),1-(err/(m+EPSILON)),1-((err_type_1*3)/(m+EPSILON)),1-((err_type_2*3)/(m+EPSILON)),1-((err_type_3*3)/(m+EPSILON)),nb_rel/(m+EPSILON),(time.time()-t)*3/(m+EPSILON),argv[0]]
        df.to_csv('scores.csv',index=False)

    elif(len(argv) == 2):
        PATH_OUT = get_temp_path()
        relecture = True
        import pandas as pd
        c = argv[0]
        f = argv[1]
        PTV,proba = def_context.load_file(str(f),str(c))
        index_PTV = PTV.index[(PTV['debut'] <= 3*60) & (PTV['debut']+PTV['DUREE'] > 3*60+5)].tolist()[0]
        def_context.Report('Starting with: %s'%(PTV['TITRE'].iloc[index_PTV]))
        lastend = PTV['debut'].loc[index_PTV]
        currentduree = PTV['DUREE'].loc[index_PTV]
        if(len(PTV) == 0):
            return("Fichier manquant")
        newPTV = def_context.init_newPTV(PTV,str(c))
        historyofpoints = def_context.init_history(str(c),PTV,lastend,currentduree)
        index_PTV = PTV.index[(PTV['debut'] <= 3*60) & (PTV['debut']+PTV['DUREE'] > 3*60)].tolist()[0]
        temp_context = historyofpoints.iloc[0]
        importantpts = def_context.get_important_points(c,PTV,index_PTV)
        help = def_context.get_help(c,PTV)
        print(help)
        for i in range(3):
            def_context.Report(str(i)+' '+str(c)+' '+str(f))
            from predictPTV import main as pred
            l1,temp_newPTV1,temp_history1,index_PTV1,temp_context1 = pred([str(c),str(f),i,newPTV.iloc[newPTV.shape[0]-1],temp_context,index_PTV,importantpts,PATH_OUT])
            if(l1>0 and relecture):
                def_context.Report("Utilisation de la relecture "+str(i)+' '+str(c)+' '+str(f))
                from RLPTV import main as RL
                l2,temp_newPTV2,temp_history2,index_PTV2,temp_context2 = RL([str(c),str(f),i,newPTV.iloc[newPTV.shape[0]-1],temp_context,index_PTV,importantpts,PATH_OUT])
                if(l2>5):
                    def_context.Report("Utilisation de l'arbre de décision",f,c,i)
                    if(chaine == 'TF1'):
                        from PTVTF1 import main as arbre
                    elif(chaine == 'M6'):
                        from PTVM6 import main as arbre
                    else:
                        l3>5
                    l3,temp_newPTV3,temp_history3,index_PTV3,temp_context3 = arbre([str(c),str(f),i,newPTV.loc[newPTV.shape[0]-1],temp_context,index_PTV,importantpts])
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

if __name__ == "__main__":
    # execute only if run as a script
    main(sys.argv[1:])
