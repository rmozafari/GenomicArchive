# -*- coding: utf-8 -*-
"""
File: Map_Flow.py
Author: rezza.mozafari@gmail.com
Created: April 05, 2024

Description:
This script serves as a critical subprocess within the upload_file.py workflow, 
specifically handling the SNP map processing and validation.
It performs sophisticated data integrity checks and database operations to ensure
the accuracy and consistency of genomic data.

Key Features:
1. Data Validation:
   - Verifies the presence and integrity of the SNP_Name field
   - Implements rigorous checks for duplicate SNP identifiers
   - Ensures data completeness by identifying and handling empty records

2. Database Integration:
   - Interfaces with SQL Server to query existing map data
   - Performs intelligent map comparison and deduplication
   - Executes dynamic SQL operations for map insertion and updates

3. Error Handling and Logging:
   - Implements a multi-level logging system for comprehensive process tracking
   - Provides granular error reporting for troubleshooting and auditing

4. File Processing:
   - Handles compressed (zip) file reading with robust error management
   - Parses complex file structures with varying delimiters

5. Parametric Configuration:
   - Utilizes external configuration files for enhanced flexibility
   - Implements dynamic parameter updates for process flow control

This subprocess is crucial for maintaining data integrity in genomic databases,
ensuring that only valid and unique SNP maps are processed and stored.
"""

import sys  
import os 
import subprocess 
import pandas as pd
import warnings
import time 
import pyodbc  
from datetime import datetime
import shutil
import zipfile as zf  
import io
import json
import logging

# Load configuration
with open("config.json", encoding='utf-8-sig') as json_config_file:
    config = json.load(json_config_file)

# Append path to custom parameters module
sys.path.append(os.getcwd().replace("programmi", ""))

# Import custom parameters
import Parametri as P

# Define paths
pj = os.path.join 
path_programmi = config["path_programmi"]
pathTemplatesDir = config["path_parametri"] + 'templates'
templatesParametri = pj(pathTemplatesDir, 'Parametri.py')
mainParametri = pj(config["path_parametri"], 'Parametri.py')

# Set up logging
logLevel = int(config["log_level"])

logger = logging.getLogger('ProGenMap')
logger.setLevel(logging.DEBUG)

formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# Info handler
info_handler = logging.FileHandler('info.log')
info_handler.setLevel(logging.INFO)
info_handler.setFormatter(formatter)
if logLevel == 1:
    logger.addHandler(info_handler)

# Warning handler
warning_handler = logging.FileHandler('warning.log')
warning_handler.setLevel(logging.WARNING)
warning_handler.setFormatter(formatter)
if logLevel <= 2:
    logger.addHandler(warning_handler)

# Error handler
error_handler = logging.FileHandler('error.log')
error_handler.setLevel(logging.ERROR)
error_handler.setFormatter(formatter)
logger.addHandler(error_handler)

def DoLog(level, message):
    """Log messages based on the specified level."""
    if level == 3:
        logger.error(message)
    elif level == 2 and logLevel <= 2:
        logger.warning(message)
    elif level == 1 and logLevel == 1:
        logger.info(message)

def aggiorna_parametri(mainPar, pathTemplates):
    """Update parameters file."""
    if os.path.exists(mainPar) and os.path.exists(pathTemplates):
        shutil.copy2(mainPar, pathTemplates)
    elif not os.path.exists(mainPar):
        DoLog(3, f'Error updating parametri.py: {mainPar} does not exist')
        print("Error")
        exit()
    elif not os.path.exists(pathTemplates):
        DoLog(3, f'Error updating parametri.py: {pathTemplates} does not exist')  
        print("Error")
        exit()

def connect_to_database():
    """Establish a database connection."""
    try:
        conn = pyodbc.connect(config["connection_string"])
        cursor = conn.cursor()
        return conn, cursor, True
    except pyodbc.Error as e:
        DoLog(3, f'Database connection error: {e}')
        print("Error")
        exit()  

def aggiorna_Esiti_Caricamento(parameter_name, parameter_value, templatesParametri, mainPar, pathTemplates):
    """Updates the specified parameter in the 'parametri.py' file."""
    if not os.path.exists(templatesParametri):
        DoLog(3, f'Error updating parametri.py: {templatesParametri} does not exist')
        print("Error")
        exit()

    if not os.path.exists(mainPar):
        DoLog(3, f'Error updating parametri.py: {mainPar} does not exist')
        print("Error")
        exit()

    with open(mainParametri, 'w') as out, open(templatesParametri, 'r') as f:
        for line in f:
            if parameter_name in line:
                line = f"{parameter_name} = '{parameter_value}'\n"
            out.write(line)
    out.close()

# Main execution
conn, cursor, status = connect_to_database()

if not status:
    DoLog(3, 'Error connecting to database')
    print("Error")
    exit()

percorso_completo = []

# Process files
for root, dirs, files in os.walk(config["path_tmp"]):
    for file in files:
        if P.Nome_Map in file:
            try:
                percorso_completo = os.path.join(root, file)
                nome_file, ext_file = os.path.splitext(file)
                dimensione = os.path.getsize(percorso_completo)
                DoLog(1, f"OK, map file {percorso_completo} found!")
            except PermissionError as e:
                DoLog(3, f'Error: {e}')
                print("Error")
                exit()             
            except FileNotFoundError as e:
                DoLog(3, f'Error: {e}')
                print("Error")
                exit()
            except Exception as e:
                DoLog(3, f'Error: {e}')
                print("Error")
                exit()

# Read zipped data
n = 0
camp = []

for simbolo in config["lista_simbolo"]:
    if n == 1:
        break
    try:
        sep = simbolo
        with zf.ZipFile(percorso_completo, 'r') as zip_file:
            file_list = zip_file.namelist()
            if len(file_list) == 1:
                with zip_file.open(file_list[0]) as file:
                    file = io.TextIOWrapper(file, 'utf-8')
                    for en, line in enumerate(file):
                        if line.startswith('Index') and line.find('Name') != -1:
                            h = line.strip().split(sep)
                            if 'Index' in h and 'Name' in h:
                                index = h.index('Name')
                                n = 1
                        elif line.startswith('Index') and line.find('Name') == -1:
                            DoLog(3, f"Error reading zipped file: {file} Although Index is present in headers, the header 'Name' was not found")
                            print("B")
                            exit()
                        elif n == 1:
                            dati = line.strip().split(sep)
                            camp.append(dati[index])
            else:
                DoLog(3, 'Error reading zipped file: There are more than one file in the zip file')
                print("C")
                exit()
    except zf.BadZipFile as e:
        DoLog(3, f'Error reading zipped file: {e}')
        print("Error")
        exit()
    except ValueError as e:
        DoLog(3, f'Error reading zipped file (ValueError): {e}')
        print("Error")
        exit()
    except Exception as e:
        DoLog(3, f'Error reading zipped file (GeneralException): {e}')
        print("Error")
        exit()

if n == 0:
    DoLog(3, 'Error reading zipped file: Header is not defined properly')
    print("B")
    exit()

# Create dataframe with SNP names
snp_newmap = pd.DataFrame(camp, columns=['SNP_Name'])

# Check for duplicate SNP names
snp_newmap_controle1 = snp_newmap.drop_duplicates(['SNP_Name'])

if len(snp_newmap) != len(snp_newmap_controle1):
    DoLog(2, 'Map to be loaded contains duplicate SNP names')
    aggiorna_Esiti_Caricamento('Esito_caricamento_Mappa', 'A', templatesParametri, mainParametri, pathTemplatesDir)
    print("A")
    exit()

# Check for SNPs without names
snp_newmap_controle2 = snp_newmap.copy()
df_remove = snp_newmap.loc[(snp_newmap['SNP_Name']=='')]
snp_newmap = snp_newmap.drop(df_remove.index)

if len(snp_newmap) != len(snp_newmap_controle2):
    DoLog(2, 'Map to be loaded contains missing SNPs')
    aggiorna_Esiti_Caricamento('Esito_caricamento_Mappa', 'A', templatesParametri, mainParametri, pathTemplatesDir)
    print("A")
    exit()
else:
    DoLog(1, 'Map to be loaded does not contain duplicate Genotype Names')
    DoLog(1, 'Map to be loaded does not contain missing SNPs')                                            

snp_newmap.reset_index(inplace=True,drop=True)

# Retrieve map information from database
query = f'SELECT Map_Name, Number_snp, Map_Alias FROM GEN.[{config["Folder_Mappa"]}]'
table_Mappe = pd.read_sql(query,conn)
Table_Mappa_nsnp=table_Mappe.copy()

# Check if there is a map with the same number of SNPs
if not sys.warnoptions:
    warnings.simplefilter("ignore")

Table_Mappa_nsnp.pop('Map_Name') 
nsnp = Table_Mappa_nsnp.set_index('Number_snp').T.to_dict('list')

tof = False

for key in nsnp: 
    if (str(len(snp_newmap)) in str(key)):
        tof = True 
        Number_snp = key
        break

if tof:
    # Map with the same number of SNPs exists
    suffix = 'a'
    
    # Identify the name of the map
    DoLog(1, 'Map with the same number of SNPs as one already loaded')
    DoLog(1, f'Map_Alias: {nsnp[Number_snp][0]}')
    DoLog(1, f'Number_snp: {Number_snp}')

    # If the map name is '554_ICAR', keep it unchanged, otherwise add a suffix.
    if P.Nome_Map == '554_ICAR':
        mappa = P.Nome_Map 
    else:
        mappa = f'{Number_snp}_{suffix}'

    found = False
    
    # Check if the map name already exists in the database.
    # If it does, increment the suffix letter and check again until a unique map name is found.
    while mappa in table_Mappe['Map_Name'].values and found == False: 
        
        # Read the existing map from SQL
        query = f'SELECT SNP_Name FROM GEN.[{mappa}]'
        snpmap = pd.read_sql(query, conn)    
        controlle = pd.merge(snpmap, snp_newmap, on='SNP_Name', how="inner")
        DoLog(1, f'Checking {mappa}')
        
        # Check if the SNP names match by merging the two files
        if len(controlle) == len(snp_newmap) and len(controlle) == len(snpmap):
            found = True
            DoLog(1, 'Map matches one already loaded')
            
            aggiorna_Esiti_Caricamento('Esito_caricamento_Mappa', 'B', templatesParametri, mainParametri, pathTemplatesDir)
            
            print("D")  # It's a Flag 
            exit()
        
        # If the map name already exists, but the number of SNPs is different, increment the suffix letter and check again
        suffix = chr(ord(suffix) + 1)
        mappa = f'{Number_snp}_{suffix}'

    if found == False:
        DoLog(1, 'Map with matching number of SNPs but different SNP names')

        # Add a record to the GEN.Mappe table
        Map_Name = mappa
        valore = [str(Map_Name), len(snp_newmap), str(P.Tipo_Chip)]
        
        query = f"INSERT INTO GEN.[{config['Folder_Mappa']}] (Map_Name,Number_snp,Map_Alias) values(?,?,?)"
        cursor.execute(query, valore)
        conn.commit()

        # Create a new table for the new map
        query = f"CREATE TABLE GEN.[{Map_Name}] ([ID] [int] IDENTITY(1,1) NOT NULL,SNP_Name nvarchar(150))"
        cursor.execute(query)
        conn.commit() 

        # Add data to the new table
        query = f"INSERT INTO GEN.[{Map_Name}] (SNP_Name) values(?)"
        valore = snp_newmap['SNP_Name'].values.tolist()
        cursor.fast_executemany = True
        cursor.executemany(query, [(x,) for x in valore])
        conn.commit()

        DoLog(1, "------------> Inserting new Record into Mappe Tables")
        DoLog(1, f"--------------> Name of the new SQL server map: {Map_Name}")
        DoLog(1, f"----------------> Alias of the new SQL server map: {P.Tipo_Chip}")
        DoLog(1, f"------------------> Number of SNPs in the new SQL server map: {len(snp_newmap)}")
        DoLog(1, "--------------------> Loaded a new SQL server map")

        aggiorna_Esiti_Caricamento('Esito_caricamento_Mappa', 'E', templatesParametri, mainParametri, pathTemplatesDir)

        print("F")  # It's a Flag 
else:
    # No map with the same number of SNPs exists
    query = f'SELECT Nume_Cari, Data_Cari, User_Cari, Tipo_Cari, Nome_file FROM GEN.[{config["Folder_Caricamento"]}]'
    table_controlli = pd.read_sql(query,conn)
    DoLog(1, 'No Map with the same number of SNPs as one already loaded')

    # Add a record to the GEN.Mappe table
    Map_Name = f"{len(snp_newmap)}_a"
    valore = [str(Map_Name), len(snp_newmap), str(P.Tipo_Chip)]
    
    query = f"INSERT INTO GEN.[{config['Folder_Mappa']}] (Map_Name,Number_snp,Map_Alias) values(?,?,?)"
    cursor.execute(query, valore)
    conn.commit()

    # Create a new table for the new map
    query = f"CREATE TABLE GEN.[{Map_Name}] ([ID] [int] IDENTITY(1,1) NOT NULL,SNP_Name nvarchar(150))"
    cursor.execute(query)
    conn.commit() 

    # Add data to the new table
    query = f"INSERT INTO GEN.[{Map_Name}] (SNP_Name) values(?)"
    valore = snp_newmap['SNP_Name'].values.tolist()
    cursor.fast_executemany = True
    cursor.executemany(query, [(x,) for x in valore])
    conn.commit()

    DoLog(1, "------------> Inserting new Record into Mappe Tables")
    DoLog(1, f"--------------> Name of the new SQL server map: {Map_Name}")
    DoLog(1, f"----------------> Alias of the new SQL server map: {P.Tipo_Chip}")
    DoLog(1, f"------------------> Number of SNPs in the new SQL server map: {len(snp_newmap)}")
    DoLog(1, "--------------------> Loaded a new SQL server map")

    aggiorna_Esiti_Caricamento('Esito_caricamento_Mappa', 'C', templatesParametri, mainParametri, pathTemplatesDir)

    print("E")  # It's a Flag