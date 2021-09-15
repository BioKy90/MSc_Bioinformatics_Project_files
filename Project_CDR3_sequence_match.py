#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 15 17:55:42 2021

@author: Ky
"""
import pandas as pd
import os
path = '/Users/kyvinguyen/Documents/Uni/Birkbeck/MSc_Bioinformatics/Research_Project/Project_python_files/'
pd.set_option('display.max_columns', None)


PDBdf = pd.read_csv(path + 'PDBSum_data/Parsed_PDB_files/1ao7_df.csv')
#TRAdf = pd.read_csv(path +'VDJDB_TRA.csv')
#TRBdf = pd.read_csv(path +'VDJDB_TRB.csv')
#TRAandB = pd.read_csv(path +'Human_CDR3_MHCI_TRAandB_PDB_all_epitopes_filtered_resolution_seqlen.csv')
VDJdb = pd.read_csv(path + 'Human_CDR3_paired_MHCI_Res_multi_Interactions.csv')
AA_sequence_df = pd.read_csv(path + 'PDBSum_data/AA_sequence_df.csv')

bCDR3seq = 'CASRPGLAGGRPEQYF'
aCDR3seq = 'CAVTTDSWGKLQF'
# find where CDR3seq starts in the residue row,
#seq = PDBdf.Residue.astype(str).str.cat()
#indexstart = seq.index(aCDR3seq)
#indexend = indexstart + len(aCDR3seq)
#
#CDR3_df = PDBdf.iloc[indexstart:indexend]
#print(CDR3_df)

## Groupby VDJdb by epitope, then create a function to iterate each group that will:
## 1. match alpha and beta CDR3 seqeunce to the AA_seq, giving Residue no. and CHAIN data.
## save the Residue no. and chain in a dataframe: CDR3_df
## 2. cross-reference the dataframe to Unique_Interactions_df, to identify the LIGAND residue and residue no. 
## add the LIGAND residue and residue no to the dataframe: to_csv(CDR3_ligand_df.csv)
## Concatenate all dfs to a master_interactions_df.csv

VDJdb_group = VDJdb.groupby(['Epitope'])
#Epitopes = []
CDR3_matched_Epitopes_df = []
seq = AA_sequence_df.Residue.astype(str).str.cat()
for name, group in VDJdb_group:
#    Epitopes.append(str(name))

    a_CDR3 = group.loc[group["Gene"] == "TRA"]
    b_CDR3 = group.loc[group["Gene"] == "TRB"]
    CDR3_df = pd.concat([a_CDR3, b_CDR3])
#    print(name)
#    print('------------')
#    print(CDR3_df)
#    print(type(a_CDR3))
    
    Ligand_df = []
    for index, row in CDR3_df.iterrows():
        CDR3 = row['CDR3']
#        print(CDR3)
        try:
            indexstart = seq.index(CDR3)
            indexend = indexstart + len(CDR3)
            CDR3_df = AA_sequence_df.iloc[indexstart:indexend]
            ## Line below: Try using .loc[row_indexer,col_indexer] = value instead
            CDR3_df['Gene'] = row['Gene']
            Ligand_df.append(CDR3_df)
        except ValueError:
            pass
#        
    Ligand_df = pd.concat(Ligand_df)
    Ligand_df['Epitope'] = name
#    print(name)
#    print('------------')
#    print(Ligand_df)
    CDR3_matched_Epitopes_df.append(Ligand_df)

CDR3_matched_Epitopes_df = pd.concat(CDR3_matched_Epitopes_df)
#CDR3_matched_Epitopes_df.to_csv(path + 'PDBSum_data/CDR3_matched_Epitopes_df.csv')
    
##  Concat all dataframes into Ligand_df
#    Ligand_df = pd.DataFrame({'idx':[1,2,3], 'dfs':[df1, df2, df3]})

CDR = 'CAVNF'
#bCDR3seq = 'CASRPGLAGGRPEQYF'
#aCDR3seq = 'CAVTTDSWGKLQF'
# find where CDR3seq starts in the residue row,
#seq = AA_sequence_df.Residue.astype(str).str.cat()
indexstart = seq.index(CDR)
indexend = indexstart + len(CDR)
#
CDR3_df = AA_sequence_df.iloc[indexstart:indexend]
#
print(CDR3_df)






#uniqueCDR3 = TRBdf['CDR3']
#print('Unique CDR3s = ' + str(len(uniqueCDR3)))
#
#uniqueEpitope = TRBdf['Epitope'].unique()
#print('Unique epitopes = ' + str(len(uniqueEpitope)))


#df_sort = TRBdf.groupby('Epitope')['Epitope'].count().reset_index(name='Count').sort_values(['Count'], ascending=False)
#TRA_sorted = TRAdf.sort_values(by=['Epitope_len','Epitope'])
#TRB_sorted = TRBdf.sort_values(by=['Epitope_len','Epitope'])
#TRAandB_sorted = TRAandB.sort_values(by=['Epitope_len','Epitope'])
#df_sort.to_csv('TRB_epitope_sort.csv', index = False)
#TRA_sorted.to_csv('TRA_epitope_sort.csv', index = False)
#TRAandB_sorted.to_csv('TRAandB_epitope_sort.csv', index = False)

#interaction_df = pd.read_csv(path + 'PDBSum_data/Interaction_df_files/1ao7_interaction_df.csv')
#
#col_names = ['Unnamed: 0', 'H bonded', 'Residue', 'Residue no.', 'Chain', 'Lig. residue', 'Lig. residue no.', 'Lig. chain', 'Distance']
#h_bond_df = pd.DataFrame(columns = col_names)
#non_h_bond_df = pd.DataFrame(columns = col_names)
#for filename in os.listdir(path + 'PDBSum_data/Interaction_df_files/'):  
#    file = str(filename)
#    df = pd.read_csv(path + 'PDBSum_data/Interaction_df_files/' + file)
#    h_bond_df = h_bond_df.append(df[(df['H bonded'] == 'Y') & (df['Chain'] == 'E') & (df['Residue no.'].between(80, 110, inclusive=True))])
#    non_h_bond_df = non_h_bond_df.append(df[(df['H bonded'] == 'N') & (df['Chain'] == 'E') & (df['Residue no.'].between(80, 110, inclusive=True))])


#print(h_bond_df['Distance'].mean())
#print(h_bond_df['Residue no.'].value_counts())
#h_bond_df.to_csv(path + 'h_bonded_interactions.csv', index = False)
#non_h_bond_df.to_csv(path + 'non_h_bonded_interactions.csv', index = False)

#print(non_h_bond_df['Distance'].mean())
#print(non_h_bond_df['Residue no.'].value_counts())
