import re 
import pandas as pd
import Bio
pd.set_option('display.max_columns', None)
path = '/Users/kyvinguyen/Documents/Uni/Birkbeck/MSc_Bioinformatics/Research_Project/Project_python_files/'


## Import Raw VDJ database: Human_CDR3_MHCI_TRAandB_all_epitopes as a dataframe 
#df = pd.read_csv('/Users/kyvinguyen/Documents/Uni/Birkbeck/MSc_Bioinformatics/Research_Project/Project_python_files/Human_CDR3_MHCI_TRAandB_all_epitopes.tsv', sep='\t')
df2 = pd.read_csv(path + 'VDJdb_paired_MHCIandII.tsv', sep='\t')


## Use regular expression pattern to filter for entries with associated PDB code
pattern = r'\"structure.id\": \"....\"'
#df = df[df['Meta'].str.contains(pattern)]
df2 = df2[df2['Meta'].str.contains(pattern)]

## Group like epitope sequences and sort by descending order
#df_sort = df.groupby('Epitope')['Epitope'].count().reset_index(name='Count').sort_values(['Count'], ascending=False)
#print(df_sort)

## Extract PDB code and create new column 'PDB'
#df['PDB'] = df['Meta'].str.extract('\"structure.id\": \"(....)\"')
df2['PDB'] = df2['Meta'].str.extract('\"structure.id\": \"(....)\"')
#PDB_ID = df['PDB']
PDB_ID2 = df2['PDB']

## Change column order
#df = df[['complex.id', 'Gene', 'CDR3', 'PDB', 'V', 'J', 'Species', 'MHC A', 'MHC B',
#       'MHC class', 'Epitope', 'Epitope gene', 'Epitope species', 'Reference',
#       'Method', 'Meta', 'CDR3fix', 'Score']]
df2 = df2[['complex.id', 'Gene', 'CDR3', 'PDB', 'V', 'J', 'Species', 'MHC A', 'MHC B',
       'MHC class', 'Epitope', 'Epitope gene', 'Epitope species', 'Reference',
       'Method', 'Meta', 'CDR3fix', 'Score']]

## Reset the index 
## Old index is saved by default. To delete, use 
## Argument 'inplace=True' replaces the old df
df2 = df2.reset_index(drop=True)


## Create columns for CDR3 and Epitope lengths
#CDR3_len = []
#for x in df['CDR3']:
#    CDR3_len.append(len(x))
#df['CDR3_len'] = CDR3_len
#Epitope_len = []
#for x in df['Epitope']:
#    Epitope_len.append(len(x))
#df['Epitope_len'] = Epitope_len

CDR3_len2 = []
for x in df2['CDR3']:
    CDR3_len2.append(len(x))
df2['CDR3_len'] = CDR3_len2
Epitope_len2 = []
for x in df2['Epitope']:
    Epitope_len2.append(len(x))
df2['Epitope_len'] = Epitope_len2

df2 = df2[['complex.id', 'Gene', 'CDR3', 'CDR3_len', 'PDB', 'V', 'J', 'Species', 'MHC A', 'MHC B',
       'MHC class', 'Epitope', 'Epitope_len', 'Epitope gene', 'Epitope species', 'Reference',
       'Method', 'Meta', 'CDR3fix', 'Score']]

#print(df2)

## Write df to csv file
#df.to_csv('Human_CDR3_MHCI_TRAandB_PDB_all_epitopes.csv', index = True)
#df2.to_csv('Human_CDR3_paired_MHCIandII.csv', index = True)

## Import protnames.lst.txt as PDB_list 
with open(path + 'protnames.lst.txt') as file:
    PDB_list = file.readlines()
    
## Filter out PDB_list entries that do not contain enough indexes             
# PDB_list = [line for line in PDB_list if len(line.split()) > 4]

def df_Resolutionfinder(df, PDB_list):
    Res_list = []
    for index, row in df.iterrows():
        pdb = row['PDB'].lower()
#        pdb = lower(pdb)
        for line in PDB_list:
            if pdb != line.split()[0]: 
                continue
            if pdb == line.split()[0]: 
                Res_list.append([line.split()[4]])
                break
            else:
                Res_list.append('None')
    ## convert list format of res values to string
    strRes_list = [str(l) for l in Res_list]
    Res_list = [i.strip('[\'\']') for i in strRes_list]
    Res_list = [float(i) for i in Res_list]

    df['Resolution'] = pd.Series(Res_list)
    df_res = df[df['Resolution'] <= 3]
    df_res = df_res.reset_index(drop=True)
    
    return df_res

## Create new df from function and export to csv
# df_new = df_Resolutionfinder(df, PDB_list)
# df_new.to_csv('Human_CDR3_MHCI_TRAandB_PDB_all_epitopes_filtered_Resolution_seqlen.csv', index = True)

## Run df through function, remove MHC class II, remove ligands with single CDR3, sort by epitope length and export. 
df2_new = df_Resolutionfinder(df2, PDB_list)
df2_new = df2_new[df2_new['MHC class'] == 'MHCI']
df2_new = df2_new[df2_new.duplicated(subset=["Epitope"], keep=False)]
df2_new = df2_new.sort_values(by=['Epitope_len','Epitope'])
df2_final = df2_new.reset_index(drop=True)
#print(df2_new)
df2_final.to_csv(path + 'Human_CDR3_paired_MHCI_Res_multi.csv', index = True)


#VDJdb_df = pd.read_csv('Human_CDR3_MHCI_TRAandB_PDB_all_epitopes_filtered_resolution_seqlen.csv')
#VDJTRA = VDJdb_df[VDJdb_df['Gene'] == 'TRA']
#VDJTRB = VDJdb_df[VDJdb_df['Gene'] == 'TRB']
#print(VDJdb_df.shape)
#print(VDJTRA.shape)
#print(VDJTRB.shape)
#TRACDR3len = VDJTRA['CDR3_len'].tolist()
#print(TRACDR3len)
#TRBCDR3len = VDJTRB['CDR3_len'].tolist()
#print(TRBCDR3len)
    
#VDJTRA.to_csv('VDJDB_TRA.csv', index = False)
#VDJTRB.to_csv('VDJDB_TRB.csv', index = False)


