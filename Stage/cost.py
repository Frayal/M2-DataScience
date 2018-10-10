#!/usr/bin/env python -W ignore::DeprecationWarning
#-*- coding: utf-8 -*-
#######################################################################
#created the 16/08/2018 14:48 by Alexis Blanchet#
#################################################
'''
Calcul de la fonction de coût sur l'ensemble des propositions de programmes TV
afin de renvoyer le plus probable. N'est appelé qu'après calcul de toutes les
possibilités.
/!| WIP /!| A NE PAS UTILISER. A MODIFIER. A TESTER.
En cours de test alors on se calme svp!
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
import def_context
import time
from subprocess import Popen
import filelock
#################################################
########### Global variables ####################
#################################################
MAX_PROCESSES = 5
LOCK = filelock.FileLock(file_name = 'cout.csv')
EXTENDED = True
#################################################
########## Important functions ##################
#################################################
def coef(x):
	if(180<x<=9*60+30):
		return 0.25
	if(9*60+30<x<=11*60):
		return 0.5
	if(11*60<x<=13*60):
		return 0.5
	if(13*60<x<=14*60):
		return 0.5
	if(13*60<x<=19*60):
		return 0.5
	if(19*60<x<=21*60):
		return 0.5
	if(21*60<x<=23*60):
		return 0.5
	else:
		return 0.25



def find_cost(date,numero,nom_chaine,i):
	'''
	création et remplissage d'un DataFrame pour calculer simplement le coûte
	enfin....simplement est un grand mot...certain PTV sont tellement WTF qu'on arrive
	même pas à calculer l'erreur. On devrait peut être directement les jeter. Enfin bref
	c'est pas simple d'être nul...

	'''
	####
	file = PATH_IN+'PTV/IPTV_'+numero+'_'+date+'_'+nom_chaine+'.csv'
	otherfile = PATH_OUT+'T'+str(i)+'/new_ptv/new_PTV_'+date+'_'+nom_chaine+'.csv'
	####
	try:
		ptv = pd.read_csv(file)
		new_ptv = pd.read_csv(otherfile)
		df = pd.DataFrame()
	except Exception as e:
		def_context.Report('petit problème: '+str(e))
		return [3+i/2000,1+i/2000,1+i/2000,1+i/2000]

	df['titre'] = ptv['TITRE']
	df['debut'] = ptv['debut']%1440
	df['duree'] = ptv['DUREE']%1440
	df['fin'] = (ptv['debut']+ptv['DUREE'])%1440
	df['vrai fin'] = 0
	df['coef'] = df['fin'].apply(lambda x:coef(x))
	df['ND'] = 0
	df['pourcentage vu'] = 0
	new_ptv_ = new_ptv[new_ptv['Évenement'].apply(lambda x: x.split(' ')[0]) == 'fin']
	current = 0
	for j in range(df.shape[0]):
	    for i in range(current,new_ptv_.shape[0]):
	        if(new_ptv_['TITRE'].iloc[i] == df['titre'].iloc[j]):
	            if(abs(df['fin'].iloc[j] - new_ptv_['minute'].iloc[i])<40 or df[df['titre'] == df['titre'].iloc[j]].shape[0] == 1):
	                df['vrai fin'].iloc[j] = new_ptv_['minute'].iloc[i]
	                df['pourcentage vu'].iloc[j] = new_ptv_['pourcentage vu'].iloc[i]
	                if(new_ptv_['Évenement'].iloc[i] == "fin d'un programme" ):
	                    df['ND'].iloc[j] = 0
	                else:
	                    df['ND'].iloc[j] = 1
	                current = i
	                break
	            else:
	                pass

	        else:
	            pass

	df['vrai debut'] = (df['vrai fin'] - df['duree']*df['pourcentage vu'])%1440
	df['vrai fin'].iloc[df.shape[0]-1] = df['vrai fin'].iloc[df.shape[0]-2] + df['duree'].iloc[df.shape[0]-1]
	df2 = df[df['pourcentage vu']==0]
	df = df[df['pourcentage vu']>0]
	df['cout'] = 0
	df2['cout'] = 1
	df = df.reset_index(drop = True)
	df2 = df2.reset_index(drop = True)
	for index, row in df.iterrows():
		df['cout'].iloc[index-1] = (min(abs(row['debut']-row['vrai debut'])%1440,abs(row['debut']%1440-row['vrai debut']%1440)) + min(abs(row['fin']-row['vrai fin'])%1440,abs(row['fin']%1440-row['vrai fin']%1440)))*row['coef']
	cout = np.sum(df['cout']) + 20*np.sum(df2['cout'])
	cout_matin = np.sum(df[(df['fin']<=13*60+40) & (df['fin']>180)]['cout']) + np.sum(df2[(df2['fin']<=13*60+40) & (df2['fin']>180)]['cout'])*20
	cout_aprem = np.sum(df[(df['fin']<=20*60+35) & (df['fin']>13*60+40)]['cout']) + 100*np.sum(df2[(df2['fin']<=20*60+35) & (df2['fin']>13*60+40)]['cout'])
	cout_soir = np.sum(df[df['fin']>20*60+35]['cout'])+np.sum(df[df['fin']<=180]['cout']) +100*(np.sum(df2[df2['fin']>20*60+35]['cout'])+np.sum(df2[df2['fin']<=180]['cout']))

	for index,x in new_ptv[['Évenement','minute']].iterrows():
		if 'HARD RESET OF ALGORITHM' in x['Évenement']:
			if(x['minute']<=13*60+40 and x['minute']>180):
				cout_matin += 2000+i
				cout += 2000+i
			elif(x['minute']<20*60+35 and x['minute']>13*60+40):
				cout_aprem += 2000+i
				cout += 2000+i
			else:
				cout_soir += 2000+i
				cout += 2000+i
	df['cout'] = df['cout']/2000
	df.to_html('test/'+date+'_'+nom_chaine+'_'+str(i)+'.html')
	return([cout/2000,cout_matin/2000,cout_aprem/2000,cout_soir/2000])


#################################################
########### Important functions #################
#################################################
def main(argv):
	global PATH_IN,PATH_SCRIPT,PATH_OUT
	PATH_IN,PATH_SCRIPT,PATH_OUT = def_context.get_path()
	if(len(argv) == 2):
		chaine = argv[0]
		date = argv[1]
		numero,nom_chaine = def_context.get_tuple(chaine)
		res = []
		files_ = os.listdir(PATH_OUT)
		for i in range(len(files_)):
			res.append(find_cost(date,numero,nom_chaine,i))
		def_context.Report(res)
		LOCK.acquire()
		try:

			couts = pd.read_csv('cout.csv')
		except Exception as e:
			couts = pd.DataFrame()
		try:
			couts[date+'_'+nom_chaine+'_tout'] = [i[0] for i in res]
			couts[date+'_'+nom_chaine+'_matinee'] = [i[1] for i in res]
			couts[date+'_'+nom_chaine+'_apresmidi'] =[i[2] for i in res]
			couts[date+'_'+nom_chaine+'_soiree'] = [i[3] for i in res]
		except Exception as e:
			def_context.Report('humm: '+str(e))
		couts.to_csv('cout.csv',index=False)
		LOCK.release()
	elif(len(argv) == 1):
		Processes = []
		files = os.listdir(PATH_IN+'PTV/')
		for file in files:
			date = file.split('_')[2]
			chaine = file.split('_')[-1].split('.')[0]
			print(type(date),type(chaine),date,chaine)
			if( str(argv[0]) in [str(date),str(chaine)]):
				print(argv[0])
				while(len(Processes)>= MAX_PROCESSES):
					lenp = len(Processes)
					for p in range(lenp):  # Check the processes in reverse order
						if Processes[enp - 1 - p].poll() is not None:  # If the process hasn't finished will return None
							del Processes[lenp - 1 - p]  # Remove from list - this is why we needed reverse order
					time.sleep(5)

				Processes.append(Popen(['python','cost.py',str(chaine),str(date)]))
				def_context.Report('calcul des coûts pour la journée du %s sur la chaîne %s'%(date,chaine))
			else:
				continue
	else:
		t = time.time()
		Processes = []
		files = os.listdir(PATH_IN+'PTV/')
		for file in files:

			print(len(Processes))
			date = file.split('_')[2]
			chaine = file.split('_')[-1].split('.')[0]
			while(len(Processes)>=5):
				time.sleep(5)
				for p in range(len(Processes)): # Check the processes in reverse order
					lenp = len(Processes)
					for p in range(lenp):  # Check the processes in reverse order
						if Processes[enp - 1 - p].poll() is not None:  # If the process hasn't finished will return None
							del Processes[lenp - 1 - p]  # Remove from list - this is why we needed reverse order
					time.sleep(5)
			Processes.append(Popen(['python','cost.py',str(chaine),str(date)]))
			def_context.Report('calcul des coûts pour la journée du %s sur la chaîne %s'%(date,chaine))
			time.sleep(2)
		while(len(Processes)):
			lenp = len(Processes)
			for p in range(lenp):  # Check the processes in reverse order
				if Processes[enp - 1 - p].poll() is not None:  # If the process hasn't finished will return None
					del Processes[lenp - 1 - p]  # Remove from list - this is why we needed reverse order
			time.sleep(5)
		def_context.Report(len(files))
		def_context.Report(time.time()-t)

if __name__ == "__main__":
	# execute only if run as a script
	main(sys.argv[1:])
