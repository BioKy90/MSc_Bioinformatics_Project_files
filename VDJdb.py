import re 
import pandas as pd
import Bio
pd.set_option('display.max_columns', None)

# Import Raw VDJ database: Human_CDR3_MHCI_TRAandB_all_epitopes as a dataframe 
df = pd.read_csv('/Users/kyvinguyen/Documents/Uni/Birkbeck/MSc_Bioinformatics/Research_Project/Project_python_files/Human_CDR3_MHCI_TRAandB_all_epitopes.tsv', sep='\t')

with open('/Users/kyvinguyen/Documents/Uni/Birkbeck/MSc_Bioinformatics/Research_Project/Project_python_files/protnames.lst.txt') as file:
    PDB_list = file.readlines()


#print(df[df['Meta'].str.contains(re.search(r'\"structure.id\": \".+\"', df))])

# 'Meta' cell data format: 
# {"cell.subset": "", "clone.id": "", "donor.MHC": "", "donor.MHC.method": "", "epitope.id": "", "replica.id": "", "samples.found": 1, "structure.id": "1zgl", 
# "studies.found": 1, "study.id": "", "subject.cohort": "", "subject.id": "", "tissue": ""}

# Use regular expression pattern to filter for entries with associated PDB code
pattern = r'\"structure.id\": \"....\"'
df = df[df['Meta'].str.contains(pattern)]

# Group like epitope sequences and sort by descending order
#df_sort = df.groupby('Epitope')['Epitope'].count().reset_index(name='Count').sort_values(['Count'], ascending=False)
#print(df_sort)

# Extract PDB code and create new column 'PDB'
df['PDB'] = df['Meta'].str.extract('\"structure.id\": \"(....)\"')
PDB_ID = df['PDB']
df = df[['complex.id', 'Gene', 'CDR3', 'PDB', 'V', 'J', 'Species', 'MHC A', 'MHC B',
       'MHC class', 'Epitope', 'Epitope gene', 'Epitope species', 'Reference',
       'Method', 'Meta', 'CDR3fix', 'Score']]

# Reset the index 
# Old index is saved by default. To delete, use 
# Argument 'inplace=True' replaces the old df
df = df.reset_index(drop=True)


CDR3_len = []
for x in df['CDR3']:
    CDR3_len.append(len(x))
df['CDR3_len'] = CDR3_len

Epitope_len = []
for x in df['Epitope']:
    Epitope_len.append(len(x))
df['Epitope_len'] = Epitope_len
#print(df.head())    

#Write new df to csv file
#df.to_csv('Human_CDR3_MHCI_TRAandB_PDB_all_epitopes.csv', index = True)

# Filter out PDB_list entries that do not contain enough indexes             
PDB_list = [line for line in PDB_list if len(line.split()) > 4]


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
    # convert list format of res values to string
    strRes_list = [str(l) for l in Res_list]
    Res_list = [i.strip('[\'\']') for i in strRes_list]
    Res_list = [float(i) for i in Res_list]
    # add 
    df['Resolution'] = pd.Series(Res_list)
    return df

#
#df_new = df_Resolutionfinder(df, PDB_list)
#df_mask = df_new['Resolution'] <= 3
#df_new = df[df_mask]
#df_new = df_new.reset_index(drop=True)
#print(df_new)
#print(df_new.Resolution.dtype)
#df_new.to_csv('Human_CDR3_MHCI_TRAandB_PDB_all_epitopes_filtered_Resolution_seqlen.csv', index = True)


#s = [['a'], ['b'], ['c']]
#ss = [str(l) for l in s]
#sss = [i.strip('[\'\']') for i in ss] 
#print(sss)

#
#df_interaction = pd.read_csv('Human_CDR3_MHCI_TRAandB_PDB_all_epitopes_filtered_resolution_interaction.csv')
#uniqueCDR3 = df_interaction['CDR3'].unique()
#print('Unique CDR3s = ' + str(len(uniqueCDR3)))
#
#uniqueEpitope = df_interaction['Epitope'].unique()
#print('Unique epitopes = ' + str(len(uniqueEpitope)))
##
#uniquePDB = df_interaction['PDB'].unique()
#print('Unique PDB codes = '+ str(len(uniquePDB)))

#PDB_list length = 182241 
#print(len(PDB_list))

#df.to_csv('Human_CDR3_MHCI_TRAandB_PDB_all_epitopes_filtered_Resolution_interaction_len.csv', index = True)

VDJdb_df = pd.read_csv('Human_CDR3_MHCI_TRAandB_PDB_all_epitopes_filtered_resolution_seqlen.csv')
VDJTRA = VDJdb_df[VDJdb_df['Gene'] == 'TRA']
VDJTRB = VDJdb_df[VDJdb_df['Gene'] == 'TRB']
print(VDJdb_df.shape)
print(VDJTRA.shape)
print(VDJTRB.shape)
TRACDR3len = VDJTRA['CDR3_len'].tolist()
print(TRACDR3len)
TRBCDR3len = VDJTRB['CDR3_len'].tolist()
print(TRBCDR3len)
    
VDJTRA.to_csv('VDJDB_TRA.csv', index = False)
VDJTRB.to_csv('VDJDB_TRB.csv', index = False)


