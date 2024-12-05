#!/usr/bin/env python3
import os
import sys
import re
import pandas as pd
import openpyxl
from openpyxl import Workbook
from openpyxl import load_workbook

# Get the command-line arguments
args = sys.argv

# Check if there are enough arguments
if len(args) < 1:
    print("python Aggregate_ABO_reports.py results_directory")
    sys.exit(1)

# Define input and output directories
input_dir = args[1]

# Function to extract exon 7 (position 422, 429) information from output tables
def parse_exon7(filename):    
    # Open the file for reading
    with open(filename, 'r') as f:
        lines = f.readlines()

        # Extract the required lines from the first section
        pos1 = int(lines[2].split(':')[1].strip())
        nuc1 = lines[3].split(':')[1].strip()
        blood1 = lines[4].split(':')[1].strip()
        count1 = int(lines[6].split(':')[1].strip())
        mat1, mis1, ins1, del1, a1, g1, c1, t1 = [int(x) for x in lines[8].split()]

        # Extract the required lines from the second section
        pos2 = int(lines[10].split(':')[1].strip())
        nuc2 = lines[11].split(':')[1].strip()
        blood2 = lines[12].split(':')[1].strip()
        count2 = int(lines[14].split(':')[1].strip())
        mat2, mis2, ins2, del2, a2, g2, c2, t2 = [int(x) for x in lines[16].split()]

    # Create a DataFrame from the extracted values
    df = pd.DataFrame({
        'Exon': ['7', '7'],
        'Position': [pos1, pos2],
        '#Reads': [count1, count2],
        'Mat': [mat1, mat2],
        'Mis': [mis1, mis2],
        'Ins': [ins1, ins2],
        'Del': [del1, del2],
        'A': [a1, a2],
        'G': [g1, g2],
        'C': [c1, c2],
        'T': [t1, t2]
    })

    # Remove the string "Exon 7 pos" from the Position column
    df['Position'] = df['Position'].apply(lambda x: x)

    # Add the 'Type' column based on values in the [Ins Del A G C T] columns
    def get_type(row):
        if row['Position'] == 422:
            if row['A'] >= 60:
                return 'B'
            elif row['C'] >= 60:
                return 'A or O'
            elif abs(row['A'] - row['C']) <= 40:
                return '(A or O) and B'
        elif row['Position'] == 429:
            if row['G'] >= 60:
                return 'A or O'
            elif row['C'] >= 60:
                return 'B'
            elif abs(row['G'] - row['C']) <= 40:
                return '(A or O) and B'
        return ''
    
    df['Type'] = df.apply(get_type, axis=1)

    # Reorder the columns
    df = df[['Exon', 'Position', '#Reads', 'Mat', 'Mis', 'Ins', 'Del', 'A', 'G', 'C', 'T', 'Type']]
    
    return df

# Function to extract exon 6 (position 22) information from output tables
def parse_exon6(filename):    
    # Open the file for reading
    with open(filename, 'r') as f:
        lines = f.readlines()

        # Extract the required lines from the first section
        pos1 = int(lines[2].split(':')[1].strip())
        count1 = int(lines[6].split(':')[1].strip())
        mat1, mis1, ins1, del1, a1, g1, c1, t1 = [int(x) for x in lines[8].split()]

    # Create a DataFrame from the extracted values
    df = pd.DataFrame({
        'Exon': ['6'],
        'Position': [pos1],
        '#Reads': [count1],
        'Mat': [mat1],
        'Mis': [mis1],
        'Ins': [ins1],
        'Del': [del1],
        'A': [a1],
        'G': [g1],
        'C': [c1],
        'T': [t1]
    })

    # Remove the string "Exon 6 pos" from the Position column
    df['Position'] = df['Position'].apply(lambda x: x)

    # Add the 'Type' column
    max_g = df['G'].max()
    max_del = df['Del'].max()
    if max_g > 60 and max_g > max_del:
        df['Type'] = 'A or B'
    elif max_del > 60 and max_del > max_g:
        df['Type'] = 'O'
    elif abs(max_g - max_del) <= 20:
        df['Type'] = 'O and (A or B)'
    else:
        df['Type'] = ''

    # Reorder the columns
    df = df[['Exon', 'Position', '#Reads', 'Mat', 'Mis', 'Ins', 'Del', 'A', 'G', 'C', 'T', 'Type']]

    return df

# Add Pheno, Geno, Expecetd columns to dataframe.|| This did not work TODO in different function
def add_phenotype_genotype(df):
    # Find the column numbers with "Type" in the column name
    type_cols = [i for i, col in enumerate(df.columns) if 'Type' in col]
    
    # Assign the column numbers to Type1, Type2, and Type3
    Type1 = df.columns[type_cols[0]]
    Type2 = df.columns[type_cols[1]]
    Type3 = df.columns[type_cols[2]]
    
    # print(Type1, Type2, Type3)
    # Add the Phenotype, Genotype, and Expected columns
    df['Phenotype'] = ''
    df['Genotype'] = ''
    df['Expected'] = ''

    for i in range(len(df)):
        if ((df.loc[i, Type1] == 'A or B').all() and (df.loc[i, Type2] == 'A or O').all() and (df.loc[i, Type3] == 'A or O')).any():
            df.at[i, 'Phenotype'] = 'A'
            df.at[i, 'Genotype'] = 'AA'
            df.at[i, 'Expected'] = 'AA'
        else:
            df.at[i, 'Phenotype'] = ''
            df.at[i, 'Genotype'] = ''
            df.at[i, 'Expected'] = ''
    return df

def assign_phenotype_genotype(df):
    exon6 = df[('Exon 6 pos 22', 'Type')]
    exon7_422 = df[('Exon 7 pos 422', 'Type')]
    exon7_429 = df[('Exon 7 pos 429', 'Type')]

    if (exon6 == 'A or B').all() and (exon7_422 == 'A or O').all() and (exon7_429 == 'A or O').all():
        Phenotype = 'A'
        Genotype = 'AA'
        Expected = 'Enter-manually'
    elif (exon6 == 'A or B').all() and (exon7_422 == 'B').all() and (exon7_429 == 'B').all():
        Phenotype = 'B'
        Genotype = 'BB'
        Expected = 'Enter-manually'
    elif (exon6 == 'A or B').all() and (exon7_422 == '(A or O) and B').all() & (exon7_429 == '(A or O) and B').all():
        Phenotype = 'AB'
        Genotype = 'AB'
        Expected = 'Enter-manually'
    elif (exon6 == 'O').all() and (exon7_422 == 'A or O').all() and (exon7_429 == 'A or O').all():
        Phenotype = 'O'
        Genotype = 'OO'
        Expected = 'Enter-manually'
    elif (exon6 == 'O and (A or B)').all() and (exon7_422 == 'A or O').all() and (exon7_429 == 'A or O').all():
        Phenotype = 'A'
        Genotype = 'AO'
        Expected = 'Enter-manually'
    elif (exon6 == 'O and (A or B)').all() and (exon7_422 == '(A or O) and B').all and (exon7_429 == '(A or O) and B').all():
        Phenotype = 'B'
        Genotype = 'BO'
        Expected = 'Enter-manually'
    else:
        Phenotype = 'Unknown'
        Genotype = 'Unknown'
        Expected = 'Enter-manually'

    df[('', 'Phenotype')] = Phenotype
    df[('', 'Genotype')] = Genotype
    df[('', 'Expected')] = Expected

    return df

# create an empty list to store the dataframes
results = []

# Define the column headers
exon6 = ['Exon 6 pos 22'] * 10
exon7_422 = ['Exon 7 pos 422'] * 10
exon7_429 = ['Exon 7 pos 429'] * 10
max_len = max(len(exon6), len(exon7_422), len(exon7_429))
exon6 += [''] * (max_len - len(exon6))
exon7_422 += [''] * (max_len - len(exon7_422))
exon7_429 += [''] * (max_len - len(exon7_429))
header_cols = ['', ''] + exon6 + exon7_422 + exon7_429 + ['', '', '']
header_rows = ['Barcode','Sample'] + ['#Reads', 'Mat', 'Mis', 'Ins', 'Del', 'A', 'G', 'C', 'T', 'Type'] * 3 + ['Phenotype','Genotype', 'Expected']
columns = pd.MultiIndex.from_arrays([header_cols, header_rows])

# Loop through all files in the input directory
for filename in os.listdir(input_dir):
    if os.path.isdir(os.path.join(input_dir, filename)):
        match = re.match(r"^IMM-[0-9]+-[0-9]+_barcode\d+", filename)
        if match:
            # process the file
            print("Processing file: " + filename)
            # Get the sample name and barcode from the folder name
            sample_name, barcode = filename.split("_")
            # Process exon 6 and 7 data using the files
            exon6_file = os.path.join(input_dir, filename, "exon6", "ABOPhenotype.txt")
            exon7_file = os.path.join(input_dir, filename, "exon7", "ABOPhenotype.txt")

            try:
                # Attempt to process the sample data
                if os.path.exists(exon6_file) and os.path.exists(exon7_file) \
                        and os.path.getsize(exon6_file) > 0 and os.path.getsize(exon7_file) > 0:
                    df_exon6 = parse_exon6(exon6_file)
                    df_exon7 = parse_exon7(exon7_file)

                    
                    # create two subset dataframes, one for each row of df for exon7
                    df_exon7_pos422 = df_exon7.iloc[[0]].reset_index(drop = True)
                    df_exon7_pos429 = df_exon7.iloc[[1]].reset_index(drop = True)
                    # print(df_exon7_pos422)
                    # print(df_exon7_pos429)
                    sample = pd.DataFrame({
                        'Barcode': [barcode],
                        'Sample': [sample_name]
                        })
                    
                    # concatenate the dataframes side by side
                    merged_df = pd.concat(
                        [sample, df_exon6, df_exon7_pos422, df_exon7_pos429],
                        axis = 1,
                        join = "inner").drop(['Exon', 'Position'],
                        axis = 1)
                    
                    # delete 'barcode' and leave only numbers in Barcode column
                    merged_df['Barcode'] = merged_df['Barcode'].str.replace('barcode', '', case = False)

                    # Convert the Barcode column to numeric
                    merged_df['Barcode'] = pd.to_numeric(merged_df['Barcode'], errors = 'coerce')

                    # Add pheno, gen, expected columns and data
                    merged_df = add_phenotype_genotype(merged_df)

                    # set the column headers to the MultiIndex
                    merged_df.columns = columns
                    merged_df = assign_phenotype_genotype(merged_df)

                    # Append the processed sample to the results list
                    results.append(merged_df)
                    print("Done adding Sample %s with barcode %s to merged data frame\n" % (barcode, sample_name))
                else:
                    print("Skipped processing Sample %s with barcode %s due to missing or empty files\n" % (barcode, sample_name))
            except Exception as e:
                print("Error processing Sample %s with barcode %s: %s\n" % (barcode, sample_name, str(e)))
                continue

# Concatenate all the dataframes in the results list into one dataframe
final_df = pd.concat(results)
final_df = final_df.sort_values(('', 'Barcode'))

# Write the concatenated dataframe to a text file
final_df.to_csv('./ABO_result.txt', sep='\t', index=False)

# Write the final dataframe to an Excel file
final_df.to_excel('./ABO_result.xlsx', sheet_name="ABO_Result", index=True)

# # Loop through all files in the input directory
# for filename in os.listdir(input_dir):
#     # Check if the file is a sample folder (for HSS Perth, these files start with IMM[0-9] with _barcode[0-9])
#     # if os.path.isdir(os.path.join(input_dir, filename)) and filename.startswith("IMM") and "barcode" in filename and "POS" not in filename:
#     if os.path.isdir(os.path.join(input_dir, filename)):
#         match = re.match(r"^IMM-[0-9]+-[0-9]+_barcode\d+", filename)
#         if match:
#             # process the file
#             print("Processing file: " + filename)
#             # Get the sample name and barcode from the folder name
#             sample_name, barcode = filename.split("_")
#             # print("\nProcessing Sample %s with barcode %s" %(barcode, sample_name))
#             # Process exon 6 and 7 data using the files
#             exon6_file = os.path.join(input_dir, filename, "exon6", "ABOPhenotype.txt")
#             exon7_file = os.path.join(input_dir, filename, "exon7", "ABOPhenotype.txt")
#             df_exon6 = parse_exon6(exon6_file)
#             df_exon7 = parse_exon7(exon7_file)
#             # create two subset dataframes, one for each row of df for exon7
#             df_exon7_pos422 = df_exon7.iloc[[0]].reset_index(drop = True)
#             df_exon7_pos429 = df_exon7.iloc[[1]].reset_index(drop = True)
#             # print(df_exon7_pos422)
#             # print(df_exon7_pos429)
#             sample = pd.DataFrame({
#                 'Barcode': [barcode],
#                 'Sample': [sample_name]
#                 })
#             # concatenate the dataframes side by side
#             merged_df = pd.concat(
#                 [sample, df_exon6, df_exon7_pos422, df_exon7_pos429],
#                 axis = 1,
#                 join = "inner").drop(['Exon', 'Position'],
#                 axis = 1)
#             # delete 'barcode' and leave only numbers in Barcode column
#             merged_df['Barcode'] = merged_df['Barcode'].str.replace('barcode', '', case = False)
#             # Convert the Barcode column to numeric
#             merged_df['Barcode'] = pd.to_numeric(merged_df['Barcode'], errors = 'coerce')
#             # Add pheno, gen, expected columns and data
#             merged_df = add_phenotype_genotype(merged_df)
#             # set the column headers to the MultiIndex
#             merged_df.columns = columns
#             # print(merged_df.shape)
#             # print(columns.shape)
#             merged_df = assign_phenotype_genotype(merged_df)
#             # append the merged_df to the results list
#             results.append(merged_df)
#             print("Done adding Sample %s with barcode %s to merged data frame\n" %(barcode, sample_name))

# # concatenate all the dataframes in the results list into one dataframe
# print("Creating merged dataframe for all samples.\n")
# final_df = pd.concat(results)
# final_df = final_df.sort_values(('', 'Barcode'))

# print(final_df.to_string(index = False), "\n")

# # Write the concatenated dataframe to a text file
# final_df.to_csv('./ABO_result.txt', sep = '\t', index = False)

# # write the final dataframe to an excel file
# final_df.to_excel('./ABO_result.xlsx', sheet_name = "ABO_Result", index = True)
