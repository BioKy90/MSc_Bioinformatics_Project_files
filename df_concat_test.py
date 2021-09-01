#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 28 16:33:32 2021

@author: Ky
"""

import pandas as pd
#import glob
import os 


path = '/Users/kyvinguyen/Documents/Uni/Birkbeck/MSc_Bioinformatics/Research_Project/Project_python_files/PDBSum_data/Interaction_df_files/'
#all_files = glob.glob(path + '/*.csv')
#

li = []

for filename in os.listdir(path):
    df = pd.read_csv(path + filename)
    name = filename.split('_')
    df['PDB'] = name[0]
    df = df[['PDB', 'H bonded', 'Residue', 'Residue no.', 'Chain', 'Lig. residue', 'Lig. residue no.', 'Lig. chain', 'Distance']]
    li.append(df)

frame = pd.concat(li, axis = 0, ignore_index = True)
print(frame)
frame.to_csv('Interactions_df.csv')