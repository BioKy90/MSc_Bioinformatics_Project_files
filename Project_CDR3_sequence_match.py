#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 15 17:55:42 2021

@author: Ky
"""
import pandas as pd
import os


PDBdf = pd.read_csv('PDBSum_data/Parsed_PDB_files/1ao7_df.csv')
#TRAdf = pd.read_csv('VDJDB_TRA.csv')
#TRBdf = pd.read_csv('VDJDB_TRB.csv')
#TRAandB = pd.read_csv('Human_CDR3_MHCI_TRAandB_PDB_all_epitopes_filtered_resolution_seqlen.csv')

CDR3seq = 'CASRPGLAGGRPEQYF'
# find where CDR3seq starts in the residue row,
seq = PDBdf.Residue.astype(str).str.cat()
indexstart = seq.index(CDR3seq)
indexend = indexstart + len(CDR3seq)

CDR3_df = PDBdf.iloc[indexstart:indexend]
print(CDR3_df)

uniqueCDR3 = TRBdf['CDR3']
#print('Unique CDR3s = ' + str(len(uniqueCDR3)))
#
uniqueEpitope = TRBdf['Epitope'].unique()
#print('Unique epitopes = ' + str(len(uniqueEpitope)))


#df_sort = TRBdf.groupby('Epitope')['Epitope'].count().reset_index(name='Count').sort_values(['Count'], ascending=False)
TRA_sorted = TRAdf.sort_values(by=['Epitope_len','Epitope'])
TRB_sorted = TRBdf.sort_values(by=['Epitope_len','Epitope'])
TRAandB_sorted = TRAandB.sort_values(by=['Epitope_len','Epitope'])
#df_sort.to_csv('TRB_epitope_sort.csv', index = False)
#TRA_sorted.to_csv('TRA_epitope_sort.csv', index = False)
#TRAandB_sorted.to_csv('TRAandB_epitope_sort.csv', index = False)

interaction_df = pd.read_csv('PDBSum_data/Interaction_df_files/1ao7_interaction_df.csv')

col_names = ['Unnamed: 0', 'H bonded', 'Residue', 'Residue no.', 'Chain', 'Lig. residue', 'Lig. residue no.', 'Lig. chain', 'Distance']
h_bond_df = pd.DataFrame(columns = col_names)
non_h_bond_df = pd.DataFrame(columns = col_names)
for filename in os.listdir('PDBSum_data/Interaction_df_files/'):  
    file = str(filename)
    df = pd.read_csv('PDBSum_data/Interaction_df_files/' + file)
    h_bond_df = h_bond_df.append(df[(df['H bonded'] == 'Y') & (df['Chain'] == 'E') & (df['Residue no.'].between(80, 110, inclusive=True))])
    non_h_bond_df = non_h_bond_df.append(df[(df['H bonded'] == 'N') & (df['Chain'] == 'E') & (df['Residue no.'].between(80, 110, inclusive=True))])


#print(h_bond_df['Distance'].mean())
#print(h_bond_df['Residue no.'].value_counts())
h_bond_df.to_csv('h_bonded_interactions.csv', index = False)
non_h_bond_df.to_csv('non_h_bonded_interactions.csv', index = False)

#print(non_h_bond_df['Distance'].mean())
#print(non_h_bond_df['Residue no.'].value_counts())
