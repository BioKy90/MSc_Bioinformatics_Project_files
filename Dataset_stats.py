#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  6 13:25:24 2021

@author: Ky
"""
path = '/Users/kyvinguyen/Documents/Uni/Birkbeck/MSc_Bioinformatics/Research_Project/Project_python_files/'
import pandas as pd
import logomaker 
## Stats 
#df_interaction = pd.read_csv(path + 'Human_CDR3_MHCI_TRAandB_PDB_all_epitopes_filtered_resolution_interaction.csv')
#df = df_interaction

#uniqueCDR3 = df['CDR3'].unique()
#print('Unique CDR3s = ' + str(len(uniqueCDR3)))
#
#uniqueEpitope = df['Epitope'].unique()
#print('Unique epitopes = ' + str(len(uniqueEpitope)))
##
#uniquePDB = df['PDB'].unique()
#print('Unique PDB codes = '+ str(len(uniquePDB)))
#

logomaker.demo('fig1b')
