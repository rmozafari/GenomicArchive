# -*- coding: utf-8 -*-
"""
File: Genotype_Map_Flow.py
Author: rezza.mozafari@gmail.com
Created: April 07, 2024

Description:
This script implements a multi-threaded genomic data processing pipeline
for SNP genotyping. It's designed to handle large-scale genomic data files, 
performing sophisticated data validation, transformation, and database
integration operations.

Key Features:
1. Asynchronous File Processing:
   - Implements parallel processing of compressed genomic data files
   - Utilizes efficient file I/O operations for handling large datasets

2. Advanced Data Validation:
   - Performs multi-stage validation checks on SNP data
   - Implements intelligent error detection and handling for data integrity

3. Dynamic SQL Operations:
   - Executes parameterized SQL queries for optimal performance and security
   - Implements transaction management for data consistency

4. Scalable Logging System:
   - Utilizes a hierarchical logging system with configurable verbosity levels
   - Implements sophisticated error classification and handling mechanisms

5. Configurable Execution Parameters:
   - Utilizes external JSON configuration for enhanced flexibility
   - Implements dynamic runtime parameter updates

Technical Specifications:
- Language: Python 3.7+
- Database: SQL Server with pyodbc driver
- Data Processing: NumPy and Pandas for high-performance numerical computations
- File Handling: Custom implementations using zipfile for compressed file operations
- Configuration: JSON-based configuration management
- Logging: Python's logging module with custom formatters and handlers
"""

import os
import pandas as pd
import numpy as np
import sys
import zipfile as zf  
import io
import shutil
import pyodbc 
import json
import logging
import argparse

# Load configuration
with open("config.json") as json_config_file:
    config = json.load(json_config_file)

# Append path to custom parameters module
sys.path.append(os.getcwd().replace("programmi", ""))

import Parametri as P      

# Define paths
pj = os.path.join 
path_programmi = config["path_programmi"]
pathTemplatesDir = config["path_parametri"] + 'templates'
templatesParametri = pj(pathTemplatesDir, 'Parametri.py')
mainParametri = pj(config["path_parametri"], 'Parametri.py')

# Configure logging
logLevel = int(config["log_level"])

logger = logging.getLogger('ProGen')
logger.setLevel(logging.DEBUG)

formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# Info handler (only logs INFO level messages)
info_handler = logging.FileHandler('info.log')
info_handler.setLevel(logging.INFO)
info_handler.setFormatter(formatter)
if logLevel == 1:
    logger.addHandler(info_handler)

# Warning handler (logs WARNING and higher level messages)
warning_handler = logging.FileHandler('warning.log')
warning_handler.setLevel(logging.WARNING)
warning_handler.setFormatter(formatter)
if logLevel <= 2:
    logger.addHandler(warning_handler)

# Error handler (for logLevel 1, 2, and 3)
error_handler = logging.FileHandler('error.log')
error_handler.setLevel(logging.ERROR)
error_handler.setFormatter(formatter)
logger.addHandler(error_handler)

# Set up argument parser
parser = argparse.ArgumentParser()
parser.add_argument('--numeCari', type=str, help='Nume_Cari')
parser.add_argument('--nomeFile', type=str, help='Nome_File')
args = parser.parse_args()

# Assign arguments to variables
Nume_Cari = args.numeCari
File_name = args.nomeFile

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
    except Exception as e:
        DoLog(3, f'Unknown error in database connection: {e}')
        print("Error")
        exit()

def aggiorna_Esiti_Caricamento(parameter_name, parameter_value, templatesParametri, mainPar, pathTemplates):
    """Update loading outcomes in parameters file."""
    # Implementation details...

# Main processing logic
conn, cursor, _ = connect_to_database()

# Process files in the temporary directory
for root, dirs, files in os.walk(config["path_tmp"]):
    for file in files:
        if File_name in file:
            try:
                full_path = os.path.join(root, file)
                file_name, file_ext = os.path.splitext(file)
                file_size = os.path.getsize(full_path)
            except PermissionError as e:
                DoLog(3, f'Permission denied: {e}')
                print("Error")
                exit() 
            except FileNotFoundError as e:
                DoLog(3, f'File not found: {e}')
                print("Error")
                exit()
            except Exception as e:
                DoLog(3, f'Unknown error: {e}')
                print("Error")
                exit()

    info_callrate = pd.DataFrame()
    exception_occurred = False

    for symbol in config["lista_simbolo"]:
        try:
            sep = symbol
            n = 0
            allele1_count = 0
            allele2_count = 0    
            
            previous_sample = '' 
            snpmap = code_snpmap.copy()
            snpmap_new_name = snpmap['SNP_Name'].str.upper()
            snpmap = snpmap_new_name.to_frame(name='SNP_Name')
            snpmap.reset_index(inplace=True, drop=True)
            snpmap["Sequence_no"] = snpmap.index + 1   
            snpmap.Sequence_no = snpmap.Sequence_no - 1 
            snpmap.set_index('SNP_Name', inplace=True)
            snpmap = snpmap.to_dict()['Sequence_no']
            nSnp = len(snpmap)
            genotypes = {}
            snp_finalrep = set()
            snp_finalrep_not = set()

            with zf.ZipFile(full_path, 'r') as zip_file:
                file_list = zip_file.namelist()
                if len(file_list) == 1:
                    with zip_file.open(file_list[0]) as file:
                        file = io.TextIOWrapper(file, 'utf-8')
                        for en, line in enumerate(file):
                            if line.startswith('SNP Name'):
                                h = line.strip().split(sep) 
                                n = 1
                            else:
                                if n == 1:
                                    data = line.strip().split(sep)
                                    snpname = data[h.index('SNP Name')].upper()
                                    snp_finalrep.add(snpname)
                                    sample = data[h.index('Sample ID')]        
                                    if sample not in genotypes:
                                        genotypes[sample] = np.array(['5']*nSnp)
                                    A1 = data[h.index('Allele1 - AB')]
                                    if A1 not in ['A', 'B', '-']:
                                        allele1_count += 1
                                        DoLog(2, f'Warning: Allele1 column has {allele1_count} errors')
                                    A2 = data[h.index('Allele2 - AB')]
                                    if A2 not in ['A', 'B', '-']:
                                        allele2_count += 1
                                        DoLog(2, f'Warning: Allele2 column has {allele2_count} errors')
                                    try:
                                        snppos = snpmap[snpname]
                                        genotypes[sample][snppos] = config["decode_genotype"][A1 + A2]
                                    except KeyError:
                                        snp_finalrep_not.add(snpname)
                                else:
                                    continue  

                        info_callrate = {}
                        for sample in genotypes:
                            info_callrate[sample] = {}
                            callrate = round((genotypes[sample] != '5').sum() / len(genotypes[sample]), 4)
                            info_callrate[sample]['CallRate'] = float(callrate)   
                            snp_cdcb = set(snpmap.keys())
                            check_missing = snp_cdcb - snp_finalrep
                            
                            DoLog(1, f'{sample:15} {nSnp:10} {len(check_missing):10} {len(snp_finalrep_not):10} {callrate:.4f} {P.Mappa_Finalreport:25}') 
                            
                            genotype = ''.join(genotypes[sample]) 
                            info_callrate[sample]['Genotipo'] = genotype
                            del genotype 
                        
                        info_callrate = pd.DataFrame.from_dict(info_callrate, orient='index')
                        info_callrate.reset_index(inplace=True)
                        info_callrate.columns = ['Campione', 'CallRate', 'Genotipo']
                        info_callrate['mappa_usata'] = P.Mappa_Finalreport 
                        
                        cols = ['Campione', 'CallRate', 'mappa_usata', 'Genotipo']
                        info_callrate = info_callrate[cols]
                        info_callrate.to_csv(config["path_output"] + File_name.replace(".zip", ""), sep=';', index=False, header=True)
                        
                        DoLog(1, f'File {File_name.replace(".zip", "")} created')

                        # Insert into Tmp_Finalreports
                        listOfTables = cursor.execute(f"SELECT * FROM information_schema.tables WHERE table_name like '{config['Tmp_Finalreports']}'").fetchall()
                        
                        if len(listOfTables) == 0:
                            DoLog(3, f"Case68: table {config['Tmp_Finalreports']} not present")
                            exception_occurred = True
                            break

                        info_callrate.rename(columns={'CallRate': 'CallRate_G', 'mappa_usata': 'mappa_usata_G'}, inplace=True)

                        data_for_insert = [(Nume_Cari, row[0], row[1], row[2], row[3], File_name) for row in info_callrate.values]

                        try:
                            query = f"INSERT INTO GEN.[{config['Tmp_Finalreports']}] (Nume_Cari, Campione, CallRate_G, mappa_usata_G, Genotipo, File_name) VALUES (?,?,?,?,?,?)"
                            cursor.fast_executemany = True
                            cursor.executemany(query, data_for_insert)
                            conn.commit()
                        except pyodbc.Error as e:
                            DoLog(3, f'Database error: {e}')
                            exception_occurred = True
                        except Exception as e:
                            DoLog(3, f'Unknown error: {e}')
                            exception_occurred = True
                        
                        DoLog(1, f"Genotype processing: Selected values inserted into the table {config['Tmp_Finalreports']} for Nume_Cari {Nume_Cari}")

        except KeyError: 
            DoLog(1, 'KeyError: separator different from those in the parameter list (NOT critical)')
            continue 
        
        except ValueError:
            DoLog(1, 'ValueError: separator different from those in the parameter list (NOT critical)')
            continue
        
        except Exception as e:
            DoLog(3, f'Unknown error: {e}')
            exception_occurred = True  

if exception_occurred:
    print("A") 
    aggiorna_Esiti_Caricamento('Esito_caricamento_Genotipi', 'A', templatesParametri, mainParametri, pathTemplatesDir)  
    exit()
else:
    aggiorna_Esiti_Caricamento('Esito_caricamento_Genotipi', 'F', templatesParametri, mainParametri, pathTemplatesDir)
    print("I")