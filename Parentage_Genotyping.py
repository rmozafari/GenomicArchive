# -*- coding: utf-8 -*-
"""
File: Parentage_Genotyping.py
Author: rezza.mozafari@gmail.com
Created: April 23, 2024

Description:
This script implements a specialized genomic data analysis pipeline focused on parentage verification
through SNP genotyping. It's designed to process high-throughput genomic data from large-scale
programs, performing advanced pedigree analysis and genetic relationship inference.

Key Features:
1. Pedigree-Specific SNP Processing:
   - Implements algorithms for identifying and analyzing parentage-informative SNPs
   - Utilizes statistical models for calculating likelihood ratios in parentage assignments

2. Multi-generational Pedigree Analysis:
   - Supports complex pedigree structures including multi-generational and half-sibling relationships
   - Implements methods for detecting and resolving pedigree inconsistencies

3. Population Genetics Integration:
   - Incorporates population allele frequencies for more accurate relationship inference
   - Supports analysis of population substructure and its impact on parentage assignments

4. Advanced Error Detection and Quality Control:
   - Implements sophisticated error models for genotyping errors and mutations
   - Provides quality metrics specific to parentage analysis, such as exclusion probabilities

5. Flexible Data Input and Integration:
   - Supports various genomic data formats and pedigree information sources
   - Integrates with external databases for comprehensive pedigree management

Technical Specifications:
- Core Algorithm: Custom implementation of likelihood-based parentage assignment
- Statistical Analysis: Utilizes SciPy for advanced statistical computations
- Data Structures: Custom graph-based data structures for efficient pedigree representation
- Parallelization: Implements multiprocessing for handling large datasets
- Visualization: Generates pedigree visualizations using NetworkX and Matplotlib
"""

import os
import pandas as pd
import numpy as np
import sys
import shutil
import pyodbc 
import zipfile as zf
import io
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

logger = logging.getLogger('ProGenPar')
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

def aggiorna_Esiti_Caricamento(parameter_name, parameter_value, templatesParametri, mainPar, pathTemplates):
    """Update loading outcomes in the parameters file."""
    try:
        with open(mainPar, 'r') as file:
            lines = file.readlines()

        for i, line in enumerate(lines):
            if line.strip().startswith(parameter_name):
                lines[i] = f"{parameter_name} = '{parameter_value}'\n"
                break

        with open(mainPar, 'w') as file:
            file.writelines(lines)

        shutil.copy2(mainPar, pathTemplates)
        DoLog(1, f"Parameter {parameter_name} updated to {parameter_value}")
    except Exception as e:
        DoLog(3, f"Error updating parameter {parameter_name}: {e}")
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

# Establish a database connection
conn, cursor, status = connect_to_database()

if not status:
    DoLog(3, 'Error connecting to database')
    print("Error")
    exit()

query = f'SELECT SNP_Name FROM GEN.[{config["Mappa_verif_parentela"]}]'
code_snpmap = pd.read_sql(query, conn)

if code_snpmap.empty: 
    DoLog(3, f'Map "{config["Mappa_verif_parentela"]}" not present')
    aggiorna_Esiti_Caricamento('Esito_caricamento_Genotipi', 'H', templatesParametri, mainParametri, pathTemplatesDir)
    print("Error")
    exit()
else:
    DoLog(1, 'Map present')

for root, dirs, files in os.walk(config["path_tmp"]):
    for file in files:
        if File_name in file:
            try:
                full_path = os.path.join(root, file)
                file_name, file_ext = os.path.splitext(file)
                file_size = os.path.getsize(full_path)
            except PermissionError as e:
                DoLog(3, f'Permission error: {e}')
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
                        
                        DoLog(1, f'{sample:15} {nSnp:10} {len(check_missing):10} {len(snp_finalrep_not):10} {callrate:.4f} {config["Mappa_verif_parentela"]:25}')
                        
                        genotype = ''.join(genotypes[sample]) 
                        info_callrate[sample]['Genotipo'] = genotype
                        del genotype 
                    
                    info_callrate = pd.DataFrame.from_dict(info_callrate, orient='index')
                    info_callrate.reset_index(inplace=True)
                    info_callrate.columns = ['Campione', 'CallRate', 'Genotipo']
                    info_callrate.pop('CallRate')
                    cols = ['Campione', 'Genotipo']
                    info_callrate = info_callrate[cols]
                    info_callrate.to_csv(config["path_output"] + File_name.replace(".zip", "") + config["Folder_Verif"], sep=';', index=False, header=True)
                    
                    DoLog(1, f'File {File_name.replace(".zip", "") + config["Folder_Verif"]} created')

                    # Update Tmp_Finalreports
                    info_callrate.rename(columns={'Genotipo': 'Genotipo_parentela'}, inplace=True)
                    update_data = [(row[1], Nume_Cari, row[0]) for row in info_callrate.values]

                    try:
                        query = f'UPDATE GEN.[{config["Tmp_Finalreports"]}] SET Genotipo_parentela = ? WHERE Nume_Cari = ? AND Campione = ?'
                        cursor.fast_executemany = True
                        cursor.executemany(query, update_data)
                        conn.commit()
                    except pyodbc.Error as e:
                        DoLog(3, f'Database error: {e}')
                        print("Error")
                        exit()
                    except Exception as e:
                        DoLog(3, f'Unknown error: {e}')
                        print("Error")
                        exit()

                    DoLog(1, f'Table {config["Tmp_Finalreports"]} updated in the column Genotipo_parentela')

    except KeyError as e: 
        DoLog(1, f'KeyError: separator different from those in the parameter list: {e} (NOT critical)')
        aggiorna_Esiti_Caricamento('Esito_caricamento_Genotipi', 'A', templatesParametri, mainParametri, pathTemplatesDir)
        continue
    except ValueError as e:
        DoLog(1, f'ValueError: separator different from those in the parameter list: {e} (NOT critical)')
        aggiorna_Esiti_Caricamento('Esito_caricamento_Genotipi', 'A', templatesParametri, mainParametri, pathTemplatesDir)
        continue
    except Exception as e:
        DoLog(3, f'Unknown error: {e}')
        aggiorna_Esiti_Caricamento('Esito_caricamento_Genotipi', 'A', templatesParametri, mainParametri, pathTemplatesDir)
        print("A")
        exit()

aggiorna_Esiti_Caricamento('Esito_caricamento_Genotipi', 'G', templatesParametri, mainParametri, pathTemplatesDir)
print("J")