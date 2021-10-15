#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 19 16:34:57 2021

@author: Ky
"""

import pandas as pd
import itertools 
path = '/Users/kyvinguyen/Documents/Uni/Birkbeck/MSc_Bioinformatics/Research_Project/Project_python_files/'
pd.set_option('display.max_columns', None)

CDR3_df = pd.read_csv(path + 'PDBSum_data/CDR3_matched_df.csv')
Interaction_df = pd.read_csv(path + 'PDBSum_data/Unique_interactions_df.csv')

#CDR3_grouped_df = CDR3_df.groupby(['PDB'])

#for name, group in CDR3_grouped_df:
##    print(group['Chain'])
#    CDR3 = name
#
#    TRA_CDR3 = group.loc[group["Gene"] == "TRA"]
#    TRB_CDR3 = group.loc[group["Gene"] == "TRB"]
#    for group in TRA_CDR3:
#        TRA_chain = group['Chain']
#        print(TRA_chain)


common_df = pd.merge(CDR3_df, Interaction_df, on = ['PDB', 'Residue', 'Chain', 'Residue no.'], right_index=True)    

common_df.drop(common_df.columns[common_df.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
common_df = common_df.sort_values(by=['Epitope', 'PDB'])
common_df.reset_index(drop=True,inplace=True)

common_df.to_csv(path + '/PDBSum_data/CDR3_interactions.csv')
Selected_interactions_df = pd.read_csv(path + 'PDBSum_data/CDR3_interactions.csv')
Selected_epitopes = ['ELAGIGILTV', 'HMTEVVRHC', 'LLFGYPVYV']
Selected_interactions_df = Selected_interactions_df[Selected_interactions_df['Epitope'].isin(Selected_epitopes)].sort_values(by=['Epitope', 'PDB'])
#Selected_interactions_df.to_csv(path + '/PDBSum_data/Selected_CDR3_interactions.csv')

def Numerical_heatmap(df):
    epitope_grouped_df = df.groupby(['Epitope'])
    # for each EPITOPE
    Epitope_dict = {}
    for name, group in epitope_grouped_df:
        x = name

        # for each CDR3
        CDR3_groups = group
        CDR3_groups = CDR3_groups.groupby(['PDB'])
        CDR3_dict = {}

        for name, group in CDR3_groups:
#            print(group)
            y = name
            CDR3_resi_counts = group['Lig. residue no.'].value_counts()
            CDR3_resi_counts = CDR3_resi_counts.to_dict()
            CDR3_dict[y] = CDR3_resi_counts
        Epitope_dict[x] = CDR3_dict
#    print(Epitope_dict)
    
    heatmaps = []
    
    # Create heatmap_df     
    for epitope in Epitope_dict:
        epitope_name = epitope
#        print(epitope)
        epilen = len(epitope_name)
        epi_resi_numbers = list(range(1,1+epilen))
        epi_resi_letters = [letter for letter in epitope]
        heatmap_df = pd.DataFrame([epi_resi_letters], columns = epi_resi_numbers)
#        print(epi_resi_letters)
        CDR3_name_list = ['PDB']
        counter_list = []
#        print(Epitope_dict[epitope].values())
        for key, value in Epitope_dict[epitope].items():
            CDR3_name = key
            CDR3_name_list.append(CDR3_name)
            CDR3_resi_counts = value
#            print(CDR3_name)
#            print(CDR3_resi_counts)
            
        # Populate heatmap_df
            counter = []
            for col in heatmap_df:
                if col in CDR3_resi_counts:
                    counter.append(CDR3_resi_counts[col])
                else:
                    counter.append(0)
            counter_list.append(counter)
#            print(epitope_name)
#            print(counter)
        heatmap_df = pd.DataFrame([epi_resi_letters, *counter_list], columns = epi_resi_numbers)
#        heatmap_df.index= CDR3_name_list
        heatmap_df.index = CDR3_name_list
#        print(CDR3_name_list)
        heatmaps.append(heatmap_df)
    return heatmaps
#    for i in heatmap_dict:
#        print(i)

    
heatmaps_allbonds = Numerical_heatmap(Selected_interactions_df)
#print(heatmaps_allbonds)

H_bonded_interactions_df = Selected_interactions_df[Selected_interactions_df['H bonded'] == 'Y']
#print(H_bonded_interactions_df)

heatmaps_H_bonded = Numerical_heatmap(H_bonded_interactions_df)
print(heatmaps_H_bonded)

#for i in heatmaps_allbonds:
#    print(i)
#for i in heatmaps_H_bonded:
#    print(i)
#print(heatmaps_allbonds)
#print(heatmaps_H_bonded)

#ELAGIGILTV = heatmaps_allbonds[0]
#HMTEVVRHC = heatmaps_allbonds[1]
#LLFGYPVYV = heatmaps_allbonds[2]
#ELAGIGILTV.to_csv(path + '/PDBSum_data/ELAGIGILTV_heatmap_allbonds.csv')
#HMTEVVRHC.to_csv(path + '/PDBSum_data/HMTEVVRHC_heatmap_allbonds.csv')
#LLFGYPVYV.to_csv(path + '/PDBSum_data/LLFGYPVYV_heatmap_allbonds.csv')


