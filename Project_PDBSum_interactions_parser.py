#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 18:02:59 2021

@author: Ky
"""
import os
import urllib.request
import re
import pandas as pd
import time
path = '/Users/kyvinguyen/Documents/Uni/Birkbeck/MSc_Bioinformatics/Research_Project/Project_python_files/'


#df = pd.read_csv('Human_CDR3_MHCI_TRAandB_PDB_all_epitopes_filtered_Resolution.csv')
df = pd.read_csv(path + 'Human_CDR3_paired_MHCI_Res_multi.csv')

# Determine protein-ligand page URL for each PDB code
pdburl = 'https://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/'
geturl = 'GetLigInt.pl?pdb='
typeurl = '&ligtype=01&ligno='

## VDJdb > via PDBcode > PDBSum protein-ligand files
## Iterate through PDB codes in VDJdb df, dowload protein-ligand interaction.txt file 
## Output: .txt files in 'Protein_ligand' folder
## E.g. 'GetPage.pl?pdbcode=' + PDB CODE + '&template=align.html&l=' + CHAIN NUMBER

#for index, row in df.iterrows():
#    PDB = row['PDB'].lower()
#    url = str(pdburl+geturl+PDB+typeurl+'01')
#    with urllib.request.urlopen(url) as f:
#        html = f.read().decode('utf-8')
#        filename = str(path + 'PDBSum_data/Protein_ligand/' + PDB + '.txt')
#        f = open(filename, 'w')
#        f.write(html)
#        f.close()
#        time.sleep(2)

## Find PDB codes with no interaction data file, remove PDB entries from dataframe
## Output: dataframe 'df_interaction' 
#No_interaction_data = [] 
#for filename in os.listdir(path + 'PDBSum_data/Protein_ligand/'):  
#    with open(path + 'PDBSum_data/Protein_ligand/'+filename, errors='ignore') as file:
#        f = file.readlines()
#        if f[0].startswith('<!DOCTYPE'):
#            PDB = filename.split('.')
#            No_interaction_data.append(str(PDB[0]))

## Delete csv files with no data
#No_interaction_data_txt = [i +'.txt' for i in No_interaction_data]
#for i in No_interaction_data_txt: 
#    os.remove(path + 'PDBSum_data/Protein_ligand/' + i)

## Update df, removing rows matching PDBs with no interaction
#df = pd.read_csv(path + 'Human_CDR3_paired_MHCI_Res_multi.csv', index_col=[0])
#df['PDB'] = df['PDB'].str.lower()
#df_interaction = df[~df.PDB.isin(No_interaction_data)].reset_index(drop=True)
#df_interaction.to_csv(path + 'Human_CDR3_paired_MHCI_Res_multi_Interactions.csv')


## Stats 
#df_interaction = pd.read_csv('Human_CDR3_MHCI_TRAandB_PDB_all_epitopes_filtered_resolution_interaction.csv')
#df = df_interaction
#
#uniqueCDR3 = df['CDR3'].unique()
#print('Unique CDR3s = ' + str(len(uniqueCDR3)))
#
#uniqueEpitope = df['Epitope'].unique()
#print('Unique epitopes = ' + str(len(uniqueEpitope)))
##
#uniquePDB = df['PDB'].unique()
#print('Unique PDB codes = '+ str(len(uniquePDB)))


## Function for parsing interaction data into a df 
def PDB_interaction_parser(PDBcode):
    PDBdirectory = path + 'PDBSum_data/Protein_ligand/' + PDBcode +'.txt'
    
    with open(PDBdirectory, errors='ignore') as file:
        f = file.read()
        fsplit = f.split('Non-bonded contacts')
        HB = fsplit[0]
        nonHB = fsplit[1]
        
    #   Nested function for parsing ligand sequence 
        letters = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
        HBsplit = HB.splitlines()
        nonHBsplit = nonHB.splitlines()
        for line in HBsplit:
            line = line.strip()
            if line.startswith('PDB code: '):
                seq = line.split()[4]
                seq = seq.split('-')
        
        # Replace 3 letter with single letter aa using list comprehension      
#         Try without list comprehension
#        for key, value in letters.items():
#            if value not in seq:
#                continue
#            index = seq.index(value)
#            seq[index] = key
#            print(index)
#        aa = ''.join(seq)
        
#        aa = [letters[i] for i in seq] 
#        aa = ''.join(aa)           
#        print(aa)
#    
        # Nested function for parsing HB contact residues 
        HBcontacts = []
        for line in HBsplit:
            line = line.strip()
            if re.match('\d+\.', line):
                HBcontacts.append(line)
    #            print(line)
    #    print(HBcontacts)
        
        # Populate HB contact Dataframe 
        def HB_contact_df(HBcontacts):
            Residue = []
            Residue_no = []
            Chain = []
            Lig_residue = []
            Lig_residue_no = []
            Lig_chain = []
            Distance = []
            Bond_type = []
            for line in HBcontacts:
                Residue.append(line.split()[3])
                Residue_no.append(line.split()[4])
                Chain.append(line.split()[5])
                Lig_residue.append(line.split()[9])
                Lig_residue_no.append(line.split()[10])
                Lig_chain.append(line.split()[11])
                Distance.append(line.split()[12])                    
                Bond_type.append('HB')
#            Residue_single = [letters[i] for i in Residue]
            Residue_single = [x if x not in letters else letters[x] for x in Residue]
#            Lig_residue_single = [letters[i] for i in Lig_residue]
            Lig_residue_single = [x if x not in letters else letters[x] for x in Lig_residue]
        #    print(Residue_single)
        #    print(Residue_no)
        #    HBcontacts_df = pd.DataFrame(columns=['Residue', 'Residue no.', 'Chain', 'Lig. residue', 'Lig. residue no.', 'Lig. chain', 'Distance'])
            HBcontacts_df = pd.DataFrame({'Residue':Residue_single, 'Residue no.': Residue_no, 'Chain': Chain, 'Lig. residue': Lig_residue_single, 'Lig. residue no.': Lig_residue_no, 'Lig. chain': Lig_chain, 'Distance': Distance})
            return HBcontacts_df
    
        # Nested function for parsing non-HB contact residues 
        non_HBcontacts = []
        for line in nonHBsplit:
            line = line.strip()
            if re.match('\d+\.', line):
                non_HBcontacts.append(line)
    #            print(line)
    #    print(nonHBcontacts)
       
        # Nested function for parsing non-HB contact residues 
        def non_HB_contacts_df(non_HBcontacts):
            Residue = []
            Residue_no = []
            Chain = []
            Lig_residue = []
            Lig_residue_no = []
            Lig_chain = []
            Distance = []
            Bond_type = []
            for line in non_HBcontacts:
                Residue.append(line.split()[3])
                Residue_no.append(line.split()[4])
                Chain.append(line.split()[5])
                Lig_residue.append(line.split()[9])
                Lig_residue_no.append(line.split()[10])
                Lig_chain.append(line.split()[11])
                Distance.append(line.split()[12])                    
                Bond_type.append('Non-HB')
            Residue_single = [x if x not in letters else letters[x] for x in Residue]
            Lig_residue_single = [x if x not in letters else letters[x] for x in Lig_residue]
        #    print(Residue_single)
        #    print(Residue_no)
        #    HBcontacts_df = pd.DataFrame(columns=['Residue', 'Residue no.', 'Chain', 'Lig. residue', 'Lig. residue no.', 'Lig. chain', 'Distance'])
            non_HBcontacts_df = pd.DataFrame({'Residue':Residue_single, 'Residue no.': Residue_no, 'Chain': Chain, 'Lig. residue': Lig_residue_single, 'Lig. residue no.': Lig_residue_no, 'Lig. chain': Lig_chain, 'Distance': Distance})
            return non_HBcontacts_df
    
        df_HB = HB_contact_df(HBcontacts)
        df_HB.insert(0, 'H bonded', 'Y')
#        print(df_HB)
        df_nonHB = HB_contact_df(non_HBcontacts)
        df_nonHB.insert(0, 'H bonded', 'N')
#        print(df_nonHB)
#        df_join = pd.concat[df_HB, df_nonHB]
        frames = [df_HB, df_nonHB]
        result = pd.concat(frames)
        return result
        
#print(PDB_interaction_parser('1ao7'))
#a = PDB_interaction_parser('1ao7')
#print(a)
        
    
#letters = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
#seq = 'SER-LEU-LEU-MET-TRP-ILE-THR-GLN-CYS-XXX'
## SER-LEU-LEU-MET-TRP-ILE-THR-GLN-VAL
#seq = seq.split('-')
#print(seq)
##for key, value in letters.items():
##    if value not in seq:
##        continue
##    index = seq.index(value)
##    seq[index] = key
##    print(index)    
##[letters.get(item,item)  for item in seq]
#seq_short = [x if x not in letters else letters[x] for x in seq]
#print(seq_short)
    

#for filename in os.listdir(path + 'PDBSum_data/Protein_ligand/'):
#    filename_split = filename.split('.')
#    PDB = filename_split[0]
##    print(filename)
#    protein_ligand_df = PDB_interaction_parser(PDB)
##    print(protein_ligand_df.head())
#    protein_ligand_df.to_csv(path + 'PDBSum_data/Interaction_df_files/' + PDB + '_interaction_df' + '.csv')

## Concatenate individual interaction dfs into single df 
li = []
for filename in os.listdir(path + 'PDBSum_data/Interaction_df_files'):
    df = pd.read_csv(path + 'PDBSum_data/Interaction_df_files/' + filename, error_bad_lines=False)
    name = filename.split('_')
    df['PDB'] = name[0]
    df = df[['PDB', 'H bonded', 'Residue', 'Residue no.', 'Chain', 'Lig. residue', 'Lig. residue no.', 'Lig. chain', 'Distance']]
    li.append(df)

frame = pd.concat(li, axis = 0, ignore_index = True)
#print(frame)
frame.to_csv(path + 'PDBSum_data/Interaction_df_files/Interactions_df.csv')


# 
#for filename in os.listdir('PDBSum_data/Protein_ligand/'):
#    print(filename)
#    with open(filename) as file:
#        f = file.read()
#        fsplit = f.split('Non-bonded contacts')
#        hb = fsplit[0]
#        nonhb = fsplit[1]
#        
#        print(nonhb)
        


#for filename in os.listdir('PDBSum_data_protein_ligand/'):
#    with open(filename) as file:
#        f = file.split()
#
#    for line in f:
#        if pdb != line.split()[0]: 
#            continue
#        if pdb == line.split()[0]: 
#            Res_list.append([line.split()[4]])
#            break
#        else:
#            Res_list.append('None')
# convert list format of res values to string
#strRes_list = [str(l) for l in Res_list]
#Res_list = [i.strip('[\'\']') for i in strRes_list]
#Res_list = [float(i) for i in Res_list]
# add 
#df['Resolution'] = pd.Series(Res_list)
#return df
