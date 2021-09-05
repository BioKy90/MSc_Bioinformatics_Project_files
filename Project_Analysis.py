#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  2 20:43:31 2021

@author: Ky
"""
import pandas as pd
import numpy as np
pd.set_option('display.max_columns', None)

path = '/Users/kyvinguyen/Documents/Uni/Birkbeck/MSc_Bioinformatics/Research_Project/Project_python_files/'
df = pd.read_csv(path + 'PDBSum_data/Interactions_df.csv')

## Duplicate residues occur where multiple distances between interacting receptor-peptide residues have been recorded. 
## drop_duplicates removes any repeating rows after the first occurance
## This will give a count of unique residue-residue interactions
df2 = df.drop_duplicates(['PDB', 'H bonded', 'Residue', 'Residue no.', 'Chain', 'Lig. residue', 'Lig. residue no.', 'Lig. chain'])
#df2.to_csv(path + 'PDBSum_data/Unique_interactions_df.csv', index = False)


df2_counts = df2.groupby(['PDB', 'H bonded']).size().to_frame(name = 'Count')

## Average no. of H bonded interactions 
df2_H = df2_counts.iloc[1::2, :]
df2_H = df2_H['Count']
print(df2_H.mean())
## Average no. of non H bonded interactions
df2_nonH = df2_counts.iloc[::2, :]
df2_nonH = df2_nonH['Count']
print(df2_nonH.mean())

#print(df2[(df2['PDB'] == '2f54') & (df2['Chain'] == 'E')])
#print(df2[df2['Chain'] == 'E'])
print(df2[(df2['Residue no.'] > 80) & (df2['Residue no.'] < 100)])

## Come back to this if 80-100 residue thing doesn't work
#df_group = df2.groupby(['PDB'])
#no_echain = []
#for name in df_group:
#    if 'E' not in df_group['Chain']:
#        no_echain += 
##    else:
##        print(name)
##        print(group[(df2['Chain'] == 'E')])
#
#print(no_echain)

