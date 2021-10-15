#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  9 14:16:46 2021

@author: Ky
"""

import os
import urllib.request
import re
import time
import pandas as pd

path = '/Users/kyvinguyen/Documents/Uni/Birkbeck/MSc_Bioinformatics/Research_Project/Project_python_files/'


df = pd.read_csv(path + 'Human_CDR3_paired_MHCI_Res_multi_Interactions.csv')



# URL for PDB file page
pdburl = 'https://files.rcsb.org/view/'

# Iterate through PDB codes in db and find associated alignment sequence text file 
# E.g. 'GetPage.pl?pdbcode=' + PDB CODE + '&template=align.html&l=' + CHAIN NUMBER
# Save sequence alignment data as txt file
for index, row in df.iterrows():
    PDB = row['PDB'].lower()
    url = str(pdburl+PDB+'.pdb')
    with urllib.request.urlopen(url) as f:
        html = f.read().decode('utf-8')
        filename = str(path + 'PDBSum_data/PDB_files/' + PDB + '.txt')
        f = open(filename, 'w')
        f.write(html)
        f.close()
        time.sleep(2)
        
### This works
#letters = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
#for filename in os.listdir(path + 'PDBSum_data/PDB_files/'): 
#    PDB = filename.split('.')[0]
#    PDBdirectory = path + 'PDBSum_data/PDB_files/' + PDB + '.txt'
#    seq_residue_3 = []
#    seq_chain = []
#    seq_index = ['0']
#    with open(PDBdirectory, errors = 'ignore') as file:
#        f = file.readlines()
#        for line in f:
#            splitline = line.split()
#            if splitline[0] == 'ATOM' and splitline[5] != seq_index[-1]:
#                if len(splitline[2]) < 3:
##            if splitline[0] == 'ATOM' and splitline[5] != seq_index[-1] and len(splitline[2]) < 3:
#                    seq_residue_3.append(splitline[3][-3:])
#                    seq_chain.append(splitline[4])
#                    seq_index.append(splitline[5])
#                else: 
#                    seq_residue_3.append(splitline[2][-3:])
#                    seq_chain.append(splitline[3])
#                    seq_index.append(splitline[4])
##            elif splitline[0] == 'ATOM' and splitline[5] != seq_index[-1] and len(splitline[2]) > 2:
##                seq_residue_3.append(splitline[2][-3:])
##                seq_chain.append(splitline[3])
##                seq_index.append(splitline[4])
#        seq_index.pop(0)
#        seq_residue = [x if x not in letters else letters[x] for x in seq_residue_3]  
#        PDB_df = pd.DataFrame({'Residue': seq_residue, 'Chain': seq_chain, 'Residue no.': seq_index}).drop_duplicates(keep='first')
#        PDB_df.to_csv(path + 'PDBSum_data/Parsed_PDB_files/' + PDB + '_df.csv')

## Concatenate all dfs into single df
li = []
for filename in os.listdir(path + 'PDBSum_data/Parsed_PDB_files/'):
    df = pd.read_csv(path + 'PDBSum_data/Parsed_PDB_files/' + filename, error_bad_lines=False)
    name = filename.split('_')
    df['PDB'] = name[0]
    df = df[['PDB', 'Residue', 'Chain', 'Residue no.']]
    li.append(df)

frame = pd.concat(li, axis = 0, ignore_index = True)
#print(frame)
frame.to_csv(path + 'PDBSum_data/Parsed_PDB_files/AA_sequence_df.csv')






## BELOW TEST CODE ONLY
# Make below a function: def PDBfile_seq_extract(PDBfile): and then cycling through folder using function
    
#letters = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
#PDBdirectory = 'PDBSum_data/PDB_files/3qeq.txt'
#with open(PDBdirectory) as file:
#    seq_residue_3 = []
#    seq_chain = []
#    seq_index = ['0']
#    f = file.readlines()
#    for line in f:
#        splitline = line.split()
#        if splitline[0] == 'ATOM' and splitline[5] != seq_index[-1]:
#            seq_residue_3.append(splitline[3])
#            seq_chain.append(splitline[4])
#            seq_index.append(splitline[5])
#    seq_index.pop(0)
#    print(len(seq_residue_3))
#    seq_residue = [letters[i] for i in seq_residue_3]  
#    PDB_df = pd.DataFrame({'Residue': seq_residue, 'Chain': seq_chain, 'Residue no.': seq_index})
#    PDB_df.to_csv('PDBSum_data/PDB_files/Parsed_PDB_files/3qeq_df.csv')


# Function input: 1ao7.txt 
    
#def PDB_sequence_extract(filename):
#    letters = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
#    PDB_code = filename.split('.')[0]
#    PDBdirectory = path + 'PDBSum_data/PDB_files/' + filename
#    seq_residue_3 = []
#    seq_chain = []
#    seq_index = ['0']
#    with open(PDBdirectory, errors = 'ignore') as file:
#        f = file.readlines()
#        for line in f:
#            splitline = line.split()
#            if splitline[0] == 'ATOM' and splitline[5] != seq_index[-1]:
#                seq_residue_3.append(splitline[3])
#                seq_chain.append(splitline[4])
#                seq_index.append(splitline[5])
#        seq_index.pop(0)
#        seq_residue = [x if x not in letters else letters[x] for x in seq_residue_3]
#        PDB_df = pd.DataFrame({'Residue': seq_residue, 'Chain': seq_chain, 'Residue no.': seq_index})
#        Parsed_directory = path + 'PDBSum_data/Parsed_PDB_files/' + str(PDB_code) + '_df.csv'
#        PDB_df.to_csv(Parsed_directory)
#
#for filename in os.listdir(path + 'PDBSum_data/PDB_files/'): 
#    PDB_sequence_extract(filename)

# DEBUG UNICODE ERROR 



#with open('PDBSum_data/2bnq.txt') as file:
#    f = file.readlines()
## Parse PDB sequence text file for CDR3 and position data.
#def PDBparser(f):
#    Echain = []
#    for line in f:
#        if line.startswith('Hydrogen bonds'):
##            return seq.group(0)
##            if line.split()[6] == 'E':
#            Echain.append(line)
#            
#            seq = list(seq.group(0))
#            residue_no = list(range(3,len(seq)+3))
#            print(len(residue_no))
#            print(len(seq))
#            df_seq = pd.DataFrame({'Amino Acid': seq, 'Residue no': residue_no})
#            return df_seq
        
#df = PDBparser(f)
##print(df)
#print(df.loc[df['Residue no'] == 95])

