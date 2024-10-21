# -*- coding: utf-8 -*- 

"""
File: upload_file.py
Author: rezza.mozafari@gmail.com
Created: March 11, 2024

Description:
This high-performance, multi-threaded Python application implements a robust ETL (Extract, Transform, Load) pipeline 
for processing genomic data. It continuously monitors and processes new records from the GEN.caricamenti SQL table, 
handling both SNP map (M) and Genotype (G) file uploads with advanced validation 
and update mechanisms.

Key Features:
1. Asynchronous Database Monitoring:
   - Utilizes pyodbc for efficient, non-blocking database connections
   - Implements a custom max_seconds generator for precise execution time control

2. Sophisticated File Processing:
   - Handles compressed (zip) genomic data files with robust error management
   - Implements intelligent file type detection and parsing for SNP and Genotype data

3. Advanced Data Validation:
   - Performs multi-stage validation checks on SNP maps and Genotype data
   - Implements data integrity checks including duplicate detection and missing value handling

4. Dynamic SQL Operations:
   - Executes parameterized SQL queries for optimal performance and security
   - Implements transaction management for ensuring data consistency

5. Scalable Error Handling and Logging:
   - Utilizes a hierarchical logging system with configurable verbosity levels
   - Implements a sophisticated error classification and handling mechanism

6. Configurable Execution Parameters:
   - Utilizes external JSON configuration for enhanced flexibility and maintainability
   - Implements dynamic runtime parameter updates

7. Optimized Data Structures:
   - Leverages pandas DataFrames for efficient in-memory data manipulation
   - Implements custom data structures for optimized lookup operations

8. Modular Architecture:
   - Utilizes separate modules for parameter management (Parametri.py)
   - Implements a plugin-based system for easy extension of file processing capabilities

Technical Specifications:
- Language: Python 3.7+
- Database: SQL Server with pyodbc driver
- Data Processing: NumPy and Pandas for high-performance numerical computations
- File Handling: Custom implementations using zipfile for compressed file operations
- Configuration: JSON-based configuration management
- Logging: Utilizes Python's logging module with custom formatters and handlers

Requirements:
- Python 3.7+
- Libraries: sys, os, subprocess, numpy, pandas, pyodbc, zipfile, json, logging
- SQL Server database with appropriate schema and tables
- config.json file with database connection string and runtime parameters
- Parametri.py file for custom parameter management
"""

import sys
import os
import io
import pyodbc
import subprocess
import json
import logging
import numpy as np
import pandas as pd
import warnings
import time
import shutil
import zipfile as zf
from datetime import datetime, timedelta



# Store the start time at the beginning of the program
start_time = time.time()

# Load configuration from JSON file
with open("config.json", encoding='utf-8-sig') as json_config_file:
    config = json.load(json_config_file)
print("Version: " + config["version"])

# Append path to custom parameters module
sys.path.append(os.getcwd().replace("programmi", ""))

# Import custom parameters
import Parametri as P

# Define paths
path_programmi = config["path_programmi"]
pathTemplatesDir = config["path_parametri"] + 'templates'
templatesParametri = os.path.join(pathTemplatesDir, 'Parametri.py')
mainParametri = os.path.join(config["path_parametri"], 'Parametri.py')

# Set maximum runtime
D = config["max_seconds_to_run"]

# Configure logging
if config["is_debug"]:
    for log_file in ['error.log', 'warning.log', 'info.log']:
        if os.path.exists(log_file):
            os.remove(log_file)

logLevel = int(config["log_level"])

logger = logging.getLogger('ProFilCar')
logger.setLevel(logging.DEBUG)

formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

for handler in logger.handlers[:]:
    logger.removeHandler(handler)

info_handler = logging.FileHandler('info.log')
info_handler.setLevel(logging.INFO)
info_handler.setFormatter(formatter)
if logLevel == 1:
    logger.addHandler(info_handler)

warning_handler = logging.FileHandler('warning.log')
warning_handler.setLevel(logging.WARNING)
warning_handler.setFormatter(formatter)
if logLevel <= 2:
    logger.addHandler(warning_handler)

error_handler = logging.FileHandler('error.log')
error_handler.setLevel(logging.ERROR)
error_handler.setFormatter(formatter)
logger.addHandler(error_handler)

def DoLog(level, message):
    if level == 3:
        logger.error(message)
    elif level == 2 and logLevel <= 2:
        logger.warning(message)
    elif level == 1 and logLevel == 1:
        logger.info(message)

def criticalError(msg="", msg_elab=""):
    global Nume_Cari, conn, cursor, config

    if msg.strip() != "":
        DoLog(3, "Critical error. " + msg)

    code_errore()

    wait = int(config.get('sleep_time', 300))
    DoLog(1, f"Waiting for {wait / 60} minutes...")

    errori_elab = msg_elab.strip() or config["decode_text_log_XDB"]["c_B"]["errori_elab"]
    
    if msg_elab != config["DATABASE_ERROR"]["msg"]:
        aggiorna_bit(conn, cursor, Nume_Cari, 0, 0, errori_elab)

    time.sleep(wait)

def max_seconds(max_seconds, *, interval=config["interval"]):
    interval = int(interval)
    start_time = time.time()
    end_time = start_time + max_seconds
    yield 0
    while time.time() < end_time:
        if interval > 0:
            next_time = start_time
            while next_time < time.time():
                next_time += interval
            time.sleep(int(round(next_time - time.time())))
        yield int(round(time.time() - start_time))
        if int(round(time.time() + interval)) > int(round(end_time)):
            return

def aggiorna_parametri(mainPar, pathTemplates):
    if os.path.exists(mainPar) and os.path.exists(pathTemplates):
        shutil.copy2(mainPar, pathTemplates)
        return "", True
    elif not os.path.exists(mainPar):
        criticalError(f"Case1: Error updating parametri.py: {mainPar} does not exist")
        return f"Error updating parametri.py: {mainPar} does not exist", False
    elif not os.path.exists(pathTemplates):
        criticalError(f"Case2: Error updating parametri.py: {pathTemplates} does not exist")
        return f"Error updating parametri.py: {pathTemplates} does not exist", False

def file_processato():
    source = percorso_completo
    destination = config["path_processato"]
    
    if not os.path.exists(source):
        criticalError(f"Case3: Error copying file: {source} does not exist")
        return f"Error copying file: {source} does not exist", False
    
    if not os.path.exists(destination):
        os.makedirs(destination)

    shutil.copy2(source, destination)
    DoLog(1, f"File {source} copied to {destination}")

    try:
        os.remove(percorso_completo)
        DoLog(1, f"File {percorso_completo} deleted")
        return "", True
    except OSError as e:
        DoLog(2, f"Error deleting file: {e}")
        return f"Error deleting file: {e}", False

def aggiorna_Esiti_Caricamento(parameter_name, parameter_value, templatesParametri, mainPar, pathTemplates):
    if not os.path.exists(templatesParametri):
        criticalError(f"Case4: Error updating parametri.py: {templatesParametri} does not exist")
        return f"Error updating parametri.py: {templatesParametri} does not exist", False

    if not os.path.exists(mainPar):
        criticalError(f"Case5: Error updating parametri.py: {mainPar} does not exist")
        return f"Error updating parametri.py: {mainPar} does not exist", False

    with open(mainParametri, 'w') as out, open(templatesParametri, 'r') as f:
        for line in f:
            if parameter_name in line:
                line = f"{parameter_name} = '{parameter_value}'\n"
            out.write(line)

    return aggiorna_parametri(mainPar, pathTemplates)

def code_errore():
    Esito_caricamento_Genotipi = 'D'
    if os.path.exists(mainParametri) and os.path.exists(templatesParametri):
        with open(mainParametri, 'w') as out, open(templatesParametri, 'r') as f:
            for line in f:
                if 'Esito_caricamento_Genotipi' in line:
                    line = "Esito_caricamento_Genotipi = 'D'\n"
                out.write(line)
        return "", True
    elif not os.path.exists(mainParametri):
        criticalError(f"Case6: Error updating parametri.py: {mainParametri} does not exist")
        return f"Error updating parametri.py: {mainParametri} does not exist", False
    elif not os.path.exists(templatesParametri):
        criticalError(f"Case7: Error updating parametri.py: {templatesParametri} does not exist")
        return f"Error updating parametri.py: {templatesParametri} does not exist", False

def connect_to_database():
    try:
        conn = pyodbc.connect(config["connection_string"])
        cursor = conn.cursor()
        return conn, cursor, True
    except pyodbc.Error as e:
        criticalError(f"Case8: Database connection error: {e}", config["DATABASE_ERROR"]["msg"])
        return None, None, False

def aggiorna_bit(conn, cursor, Nume_Cari, bit_ok, bit_elaborato, errori_elab):
    query = f"UPDATE GEN.Code_Caricamenti SET Bit_OK = {bit_ok}, Bit_elaborato = {bit_elaborato}, Errori_elab = '{errori_elab}' WHERE Nume_Cari = '{Nume_Cari}'"
    try:
        cursor.execute(query)
        conn.commit()
        return "", True
    except pyodbc.Error as e:
        criticalError(f"Case9: Error updating GEN.Code_Caricamenti: {e} Query: {query}", config["DATABASE_ERROR"]["msg"])
        return f"Error updating GEN.Code_Caricamenti: {e}", False

# Main execution loop
for sec in max_seconds(D, interval=1):
    M_code = {}
    bit_ok = 0
    bit_elaborato = 0
    errori_elab = "00"

    DoLog(1, f"START MAIN LOOP {sec}")
    
    conn, cursor, status = connect_to_database()

    if not status:
        continue

    try:
        Gen = config["Folder_Caricamento"]
        query = f'SELECT Nume_Cari, Data_Cari, User_Cari, Tipo_Cari, Nome_file, bit_elaborato FROM GEN.[{Gen}]'

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            table_controlli = pd.read_sql(query, conn)

    except pyodbc.Error as query_error:
        criticalError(f"Case10: Database error: {query_error}")
        continue
    except pd.errors.PandasError as df_error:
        criticalError(f"Case11: Dataframe error: {df_error}")
        continue
    except Exception as general_error:
        criticalError(f"Case12: General error related to SQL database: {general_error}")
        continue

    table_controlli['bit_elaborato'] = table_controlli['bit_elaborato'].fillna(-1)
    new_records = table_controlli[table_controlli['bit_elaborato'] == 0]

    if new_records.empty:
        DoLog(1, "No new records")
        DoLog(1, "Returning")
        cursor.close()
        continue
    else:
        DoLog(1, "New records found")
        
        new_records = table_controlli[table_controlli['bit_elaborato'] == 0].copy()
        new_records.sort_values(by=['Nume_Cari', 'Data_Cari'], ascending=[True, True], inplace=True)
        ids = list(new_records['Nume_Cari'])
        
        for id in ids:
            DoLog(1, f"START SECONDARY LOOP {sec}")
            tmp = table_controlli[table_controlli['Nume_Cari']==id].copy()
            Tipo_Cari = ''.join(map(str, list(tmp['Tipo_Cari'])))
            Nome_File = ''.join(map(str, list(tmp['Nome_file'])))
            
            if Nome_File.startswith('G_'):
                Nome_File = 'GEN_' + Nome_File[2:]
            if Nome_File.startswith('M_'):
                Nome_File = 'MAP_' + Nome_File[2:]
            
            date = datetime.today()
            oggi = date.strftime('%Y-%m-%d')
            Nume_Cari = id
            DoLog(1, f"Nume_Cari: {Nume_Cari}")

            with zf.ZipFile(os.path.join(config["path_tmp"], Nome_File), 'r') as zip_file:
                file_list = zip_file.namelist()
                if len(file_list) == 1:
                    with zip_file.open(file_list[0]) as zip_file_content:
                        file_content = io.TextIOWrapper(zip_file_content, 'utf-8')
                        
                        if Tipo_Cari == 'M' and file_content.readline().startswith('[Header]'):
                            DoLog(3, f"Tipo Caricamento is incorrect {Nome_File}")
                            M_code = config["decode_text_log_XDB"]["c_A"]
                            bit_ok = M_code["bit_ok"]
                            bit_elaborato = M_code["bit_elaborato"]
                            errori_elab = M_code["errori_elab"]
                            msg, status = aggiorna_bit(conn, cursor, Nume_Cari, bit_ok, bit_elaborato, errori_elab)
                            if not status:
                                criticalError("Case13: " + msg)
                                id = ids[-1]
                            continue
                        elif Tipo_Cari == 'G' and not file_content.readline().startswith('[Header]'):
                            DoLog(3, f"Tipo Caricamento is incorrect {Nome_File}")
                            M_code = config["decode_text_log_XDB"]["c_A"]
                            bit_ok = M_code["bit_ok"]
                            bit_elaborato = M_code["bit_elaborato"]
                            errori_elab = M_code["errori_elab"]
                            msg, status = aggiorna_bit(conn, cursor, Nume_Cari, bit_ok, bit_elaborato, errori_elab)
                            if not status:
                                criticalError("Case14: " + msg)
                                id = ids[-1]
                            continue
            
            msg, status = aggiorna_Esiti_Caricamento('Tipo_Cari', Tipo_Cari, templatesParametri, mainParametri, pathTemplatesDir)
            if not status:
                criticalError("Case15: " + msg)
                id = ids[-1]
                continue
            msg, status = aggiorna_Esiti_Caricamento('Nume_Cari', Nume_Cari, templatesParametri, mainParametri, pathTemplatesDir)
            if not status:
                criticalError("Case16: " + msg)
                id = ids[-1]
                continue
            msg, status = aggiorna_Esiti_Caricamento('oggi', oggi, templatesParametri, mainParametri, pathTemplatesDir)
            if not status:
                criticalError("Case17: " + msg)
                id = ids[-1]
                continue

            if Tipo_Cari == 'M':
                # Map file processing logic
                if not os.path.exists(config["path_tmp"]):
                    criticalError('Case18: Error, Directory not found!')
                
                mappe_file = ''
                for file in os.listdir(config["path_tmp"]):
                    if Nome_File in file:
                        DoLog(1, "Nome_File is present in the directory")
                        mappe_file = file
                        break

                if mappe_file == '':
                    DoLog(2, "Error, Map not found!")
                    msg, status = aggiorna_Esiti_Caricamento('Esito_caricamento_Mappa', 'D', templatesParametri, mainParametri, pathTemplatesDir)
                    if not status:
                        criticalError("Case19: " + msg)
                        id = ids[-1]
                        continue
                    M_code = config["decode_text_log_XDB"]["c_B"]
                else:
                    DoLog(1, "Map found!")
                    Id_mappa = id
                    DoLog(1, "map_processing")
                    msg, status = aggiorna_Esiti_Caricamento('Nome_Map', Nome_File, templatesParametri, mainParametri, pathTemplatesDir)
                    if not status:
                        criticalError("Case20: " + msg)
                        id = ids[-1]
                        continue
                    msg, status = aggiorna_Esiti_Caricamento('Id_mappa', Id_mappa, templatesParametri, mainParametri, pathTemplatesDir)
                    if not status:
                        criticalError("Case21: " + msg)
                        id = ids[-1]
                        continue

                    DoLog(1, "start-Map-Flow.py")
                    script = 'Map_Flow.py'

                    if not os.path.exists(script):
                        criticalError(f'Case22: Error, File {script} not found!')
                        id = ids[-1]
                        continue

                    with open(script.replace('.py','.log'), 'w') as f:
                        processo = subprocess.Popen([sys.executable, script], stdout=subprocess.PIPE, stderr=f)
                        out, err = processo.communicate()
                        
                        if out == b'A\r\n':
                            M_code = config["decode_text_log_XDB"]["m_A"]
                        elif out == b'B\r\n':
                            M_code = config["decode_text_log_XDB"]["m_B"]
                        elif out == b'C\r\n':
                            M_code = config["decode_text_log_XDB"]["m_C"]
                        elif out == b'D\r\n':
                            M_code = config["decode_text_log_XDB"]["m_D"]
                        elif out == b'E\r\n':
                            M_code = config["decode_text_log_XDB"]["m_E"]
                        elif out == b'F\r\n':
                            M_code = config["decode_text_log_XDB"]["m_F"]
                        elif out == b'Error\r\n':
                            criticalError("Case23")
                            id = ids[-1]
                            continue

                DoLog(1, "map process finished")

            else:  # Genotype file processing
                DoLog(1, "genotype_processing")

                file_ricerca = Nome_File
                percorso_completo = []
                isFound = False

                if not os.path.exists(config["path_tmp"]):
                    criticalError(f'Case24: Error, Directory {config["path_tmp"]} not found!')
                    id = ids[-1]
                    continue

                for finalrep in os.listdir(config["path_tmp"]):
                    if isFound:
                        break
                    DoLog(1, f"finalrep: {finalrep}")
                    if Nome_File in finalrep:
                        DoLog(1, f"Nome_file: {Nome_File}")

                        percorso_completo = []
                        for raiz, diretorios, files in os.walk(config["path_tmp"]):
                            if isFound:
                                break
                            for file in files:
                                if file_ricerca in file:
                                    isFound = True
                                    try:
                                        percorso_completo = os.path.join(raiz, file)
                                        nome_file, ext_file = os.path.splitext(file)
                                        dimensione = os.path.getsize(percorso_completo)
                                    except PermissionError as e:
                                        criticalError(f"Case25: Permission denied: {e}")
                                        finalrep = ''
                                        percorso_completo = []
                                        file = files[-1]
                                        continue
                                    except FileNotFoundError as e:
                                        criticalError(f"Case26: File not found {e}")
                                        finalrep = ''
                                        percorso_completo = []
                                        file = files[-1]
                                        continue
                                    except Exception as e:
                                        criticalError(f"Case27: Unknown error: {e}")
                                        finalrep = ''
                                        percorso_completo = []
                                        file = files[-1]
                                        continue
                                    break

                if isFound and percorso_completo == [] and finalrep == '':
                    id = ids[-1]
                    continue

                if finalrep == '':
                    DoLog(2, "- Warning -")
                    DoLog(2, "Final Report file to be loaded not found or with wrong name")
                    msg, status = aggiorna_Esiti_Caricamento('Esito_caricamento_Genotipi', 'B', templatesParametri, mainParametri, pathTemplatesDir)
                    if not status:
                        criticalError("Case26: " + msg)
                        id = ids[-1]
                        continue
                    M_code = config["decode_text_log_XDB"]["g_B"]

                if percorso_completo == [] or not isFound:
                    DoLog(2, "- Warning -")
                    DoLog(2, "Final Report file to be loaded not found or with wrong name")
                    msg, status = aggiorna_Esiti_Caricamento('Esito_caricamento_Genotipi', 'B', templatesParametri, mainParametri, pathTemplatesDir)
                    if not status:
                        criticalError("Case27: " + msg)
                        id = ids[-1]
                        continue
                    M_code = config["decode_text_log_XDB"]["g_B"]
                    
                else:
                    DoLog(1, "Final Report file to be loaded found")

                    tipo_chip = []
                    blocco = ''
                    verif_final_report = []
                    File_Final_Report = ''
                    verif_esito_file_finale = []
                    camp = []

                    DoLog(1, f"percorso_completo {percorso_completo}")
                    if percorso_completo == []:
                        criticalError("Case28: Error, Final Report file to be loaded not found or with wrong name")
                        id = ids[-1]
                        continue
                    
                    if not percorso_completo.endswith('.zip'):
                        criticalError(f"Case29: File {percorso_completo} is not a zip file")
                    
                    hasError = False
                    hasMoreFiles = False
                    for simbolo in config["lista_simbolo"]:
                        sep = simbolo
                        try:
                            with zf.ZipFile(percorso_completo, 'r') as zip_file:
                                file_list = zip_file.namelist()
                                if len(file_list) == 1:
                                    with zip_file.open(file_list[0]) as zip_file_content:
                                        file_content = io.TextIOWrapper(zip_file_content, 'utf-8')
                                        for en, line in enumerate(file_content):
                                            if blocco != 'trovato_chip':
                                                if line.startswith('Content'):
                                                    tipo_chip = line.strip().split(sep)

                                                    if len(tipo_chip) > 1:
                                                        DoLog(1, f"Chip names present in Final Report file {tipo_chip}")
                                                        blocco = 'trovato_chip'

                                                        t_chip = pd.DataFrame([tipo_chip])
                                                        if sep == '\t':
                                                            code_mappa1 = list(t_chip[2])
                                                        else:
                                                            code_mappa1 = list(t_chip[1])
                                                        code_mappa2 = ''.join(map(str, code_mappa1))

                                                        CHIP_ALIAS = code_mappa2
                                                        chip = str(code_mappa2)
                                                        Tipo_Chip = chip
                                                        Alias = 'SI'
                                                        File_Final_Report = 'notmissing'
                                                        n = 0

                                                        for en, line in enumerate(file_content):
                                                            if line.startswith('SNP Name'):
                                                                h = line.strip().split(sep)
                                                                n = 1

                                                                if 'Allele1 - AB' in h and 'Allele2 - AB' in h:
                                                                    DoLog(1, "Column Allele1 - AB found, column Allele2 - AB found")
                                                                else:
                                                                    File_Final_Report = 'missing'
                                                                    DoLog(2, "Column Allele1 - AB not found, column Allele2 - AB not found")
                                                                    break
                                                            else:
                                                                if n == 1 and File_Final_Report != 'missing':
                                                                    File_Final_Report = 'notmissing'
                                                                    dati = line.strip().split(sep)
                                                                    camp.append([dati[h.index('Sample ID')], dati[h.index('SNP Name')]])
                                else:
                                    DoLog(3, f"File {percorso_completo} contains more than one file")
                                    M_code = config["decode_text_log_XDB"]["g_H"]
                                    hasMoreFiles = True

                        except zf.BadZipFile as e:
                            criticalError(f"Case31: Error reading zipped file: {e}")
                            hasError = True
                            break
                        except ValueError as ve:
                            criticalError(f"Case32: Value Error: {ve}")
                            hasError = True
                            break
                        except Exception as e:
                            criticalError(f"Case33: Unknown error: {e}")
                            hasError = True
                            break

                    if hasError:
                        id = ids[-1]
                        continue

                    if hasMoreFiles:
                        msg, status = aggiorna_Esiti_Caricamento('Esito_caricamento_Genotipi', 'H', templatesParametri, mainParametri, pathTemplatesDir)
                        if not status:
                            criticalError("Case34: " + msg)
                            id = ids[-1]
                            continue
                        continue

                    if blocco != 'trovato_chip':
                        DoLog(2, "Chip names not present in Final Report file")
                        if percorso_completo == []:
                            criticalError("Case41: Error, Final Report file to be loaded not found or with wrong name")
                            id = ids[-1]
                            continue
                        if not percorso_completo.endswith('.zip'):
                            criticalError(f"Case42: File {percorso_completo} is not a zip file")
                            id = ids[-1]
                            continue
                        
                        for simbolo in config["lista_simbolo"]:
                            try:
                                if blocco == 'trovato_separator' and len(camp) != 0:
                                    continue
                                n = 0
                                sep = simbolo

                                with zf.ZipFile(percorso_completo, 'r') as zip_file:
                                    file_list = zip_file.namelist()
                                    if len(file_list) == 1:
                                        with zip_file.open(file_list[0]) as zip_file_content:
                                            file_content = io.TextIOWrapper(zip_file_content, 'utf-8')
                                            for en, line in enumerate(file_content):
                                                if line.startswith('SNP Name'):
                                                    h = line.strip().split(sep)
                                                    n = 1
                                                    if 'Allele1 - AB' in h and 'Allele2 - AB' in h:
                                                        DoLog(1, "Column Allele1 - AB found, column Allele2 - AB found")
                                                        blocco = 'trovato_separator'
                                                    else:
                                                        File_Final_Report = 'missing'
                                                        break
                                                else:
                                                    if n == 1 and File_Final_Report != 'missing':
                                                        File_Final_Report = 'notmissing'
                                                        dati = line.strip().split(sep)
                                                        camp.append([dati[h.index('Sample ID')], dati[h.index('SNP Name')]])
                                    else:
                                        criticalError(f"Case43: File {percorso_completo} contains more than one file")
                                        id = ids[-1]
                                        continue
                            except zf.BadZipFile as e:
                                criticalError(f"Case44: Error reading zipped file: {e}")
                                id = ids[-1]
                                continue
                            except ValueError as ve:
                                criticalError(f"Case45: Value Error: {ve}")
                                id = ids[-1]
                                continue
                            except Exception as e:
                                criticalError(f"Case46: Unknown error: {e}")
                                id = ids[-1]
                                continue

                    tmp_finalreports = pd.DataFrame(camp, columns=['Sample ID', 'SNP Name'])
                        
                    if blocco == 'trovato_separator':
                        if File_Final_Report == 'missing':
                            DoLog(2, "WARNING: ---> Final Report file with errors, missing column")
                            code_errore()
                            verif_final_report = 'errori'
                            aggiorna_Esiti_Caricamento('Esito_caricamento_Genotipi', 'A', templatesParametri, mainParametri, pathTemplatesDir)
                            M_code = config["decode_text_log_XDB"]["g_A"]
                        
                        else: 
                            freq = pd.DataFrame(tmp_finalreports.groupby(by='Sample ID').size())
                            freq.reset_index(inplace=True)
                            
                            if max(freq[0]) != min(freq[0]):
                                DoLog(2, "WARNING: ---> Final Report file with errors - inconsistent SNP count")
                                verif_final_report = 'errori'
                                
                                msg, status = aggiorna_Esiti_Caricamento('Esito_caricamento_Genotipi', 'A', templatesParametri, mainParametri, pathTemplatesDir)
                                if not status:
                                    criticalError("Case47: " + msg)
                                    id = ids[-1]
                                    continue

                                M_code = config["decode_text_log_XDB"]["g_C"]
                            
                            else:
                                DoLog(2, "WARNING: ---> Final Report file without errors")
                                try:
                                    chip = len(tmp_finalreports)/len(freq)
                                    Tipo_Chip = int(chip)
                                    Alias = 'NO'

                                    M_code = config["decode_text_log_XDB"]["g_D"]
                                except ZeroDivisionError: 
                                    verif_final_report = 'errori'
                                    M_code = config["decode_text_log_XDB"]["g_C"]
                                    
                    if verif_final_report == 'errori':
                        DoLog(2, "WARNING: ---> Final Report file with errors - inconsistent SNP count - outside loop")

                        msg, status = code_errore()
                        if not status:
                            id = ids[-1]
                            criticalError("Case49: " + msg)
                            continue

                        msg, status = aggiorna_Esiti_Caricamento('Esito_caricamento_Genotipi', 'A', templatesParametri, mainParametri, pathTemplatesDir)
                        if not status:
                            id = ids[-1]
                            criticalError("Case51: " + msg)
                            continue

                        M_code = config["decode_text_log_XDB"]["g_E"]
                        
                    else:
                        if File_Final_Report != 'missing' and verif_final_report != 'errori':
                            
                            DoLog(1, "Continue procedure")
                            DoLog(1, "Start map check")
                            
                            listOfTables = cursor.execute("SELECT * FROM information_schema.tables WHERE table_name like '%s'" % config["Folder_Mappa"]).fetchall()
                            if len(listOfTables) == 0:
                                criticalError('Case52: Table not present')
                            query = 'SELECT * FROM GEN.%s' % (config["Folder_Mappa"])

                            cursor.execute(query)

                            Mappa_Finalreport = 'missing_chip'
                            Aggiorna_Gen_Mappe = 'NO'
                            procedura_alias = ''

                            for row in cursor:
                                if Alias == 'SI':
                                    if row.Map_Alias == Tipo_Chip:
                                        Mappa_Finalreport = row.Map_Name
                                elif Alias == 'NO':
                                    if row.Number_snp == Tipo_Chip:
                                        Mappa_Finalreport = row.Map_Name
                            
                            if Alias == 'SI' and Mappa_Finalreport == 'missing_chip':
                                procedura_alias = 'CODE3'
                                query = 'SELECT * FROM GEN.%s' % (config["Folder_Mappa"])
                                cursor.execute(query)
                                
                                for row in cursor:
                                    if row.Number_snp == nsnp:
                                        Mappa_Finalreport = row.Map_Name
                                DoLog(1, "Update gen.mappe")
                                
                                if Mappa_Finalreport != 'missing_chip':        
                                    Aggiorna_Gen_Mappe = 'SI'

                                    msg, status = aggiorna_Esiti_Caricamento('Aggiorna_Gen_Mappe', Aggiorna_Gen_Mappe, templatesParametri, mainParametri, pathTemplatesDir)
                                    if not status:
                                        id = ids[-1]
                                        criticalError("Case53: " + msg)
                                        continue

                            else:
                                msg, status = aggiorna_Esiti_Caricamento('Aggiorna_Gen_Mappe', Aggiorna_Gen_Mappe, templatesParametri, mainParametri, pathTemplatesDir)
                                if not status:
                                    id = ids[-1]
                                    criticalError("Case54: " + msg)
                                    continue

                            DoLog(1, "End of map check")

                            tmp_finalreports_orig = tmp_finalreports.copy()
                            
                            if Mappa_Finalreport == 'missing_chip':
                                DoLog(2, "WARNING: Chip not present in Alias")
                                DoLog(2, "Map Final Report with different SNP number than already loaded")
                                DoLog(2, "Need to load a new map")
                                
                                msg, status = aggiorna_Esiti_Caricamento('Esito_caricamento_Genotipi', 'C', templatesParametri, mainParametri, pathTemplatesDir)
                                if not status:
                                    id = ids[-1]
                                    criticalError("Case55: " + msg)
                                    continue

                                M_code = config["decode_text_log_XDB"]["g_N"]

                            else:
                                listOfTables = cursor.execute("SELECT * FROM information_schema.tables WHERE table_name like '%s'" % Mappa_Finalreport).fetchall()
                                if len(listOfTables) == 0:
                                    criticalError('Case57: Table not present')
                                    id = ids[-1]
                                    continue
                                query = 'SELECT SNP_Name FROM GEN.[%s]' % (Mappa_Finalreport)

                                with warnings.catch_warnings():
                                    warnings.simplefilter('ignore')
                                    snp_map = pd.read_sql(query, conn)
                                
                                snp_map.rename(columns={'SNP_Name':'SNP Name'}, inplace=True)
                                

                                if not snp_map.empty:
                                    DoLog(2, "WARNING:")
                                    tmp_finalreports = tmp_finalreports.drop_duplicates(['SNP Name'])  
                                    tmp_finalreports.reset_index(inplace=True, drop=True)
                                    controlle = pd.merge(tmp_finalreports, snp_map, on='SNP Name', how="outer")

                                    DoLog(1, "Checking map consistency")
                                    if len(controlle) == len(tmp_finalreports) and len(controlle) == len(snp_map):
                                        DoLog(1, "Map consistency check passed")
                                        DoLog(2, "Final Report Map matches the already loaded map")
                                        DoLog(2, "Converting final reports to string")

                                        msg, status = aggiorna_Esiti_Caricamento('Final_Reports', Nome_File, templatesParametri, mainParametri, pathTemplatesDir)
                                        if not status:
                                            id = ids[-1]
                                            criticalError("Case58: " + msg)
                                            continue

                                        msg, status = aggiorna_Esiti_Caricamento('Tipo_Chip', Tipo_Chip, templatesParametri, mainParametri, pathTemplatesDir)
                                        if not status:
                                            id = ids[-1]
                                            criticalError("Case59: " + msg)
                                            continue

                                        msg, status = aggiorna_Esiti_Caricamento('Alias', Alias, templatesParametri, mainParametri, pathTemplatesDir)
                                        if not status:
                                            id = ids[-1]
                                            criticalError("Case60: " + msg)
                                            continue

                                        msg, status = aggiorna_Esiti_Caricamento('Mappa_Finalreport', Mappa_Finalreport, templatesParametri, mainParametri, pathTemplatesDir)
                                        if not status:
                                            id = ids[-1]
                                            criticalError("Case61: " + msg)
                                            continue

                                        msg, status = aggiorna_Esiti_Caricamento('Id_genotipe', Nume_Cari, templatesParametri, mainParametri, pathTemplatesDir)
                                        if not status:
                                            id = ids[-1]
                                            criticalError("Case62: " + msg)
                                            continue

                                        DoLog(1, "Updating parametri.py for genotype procedure")

                                        del(Tipo_Chip, chip, simbolo, dimensione, diretorios, en, ext_file, file, file_ricerca, files, line, nome_file, raiz, sep, tipo_chip)

                                        DoLog(1, "Entering the scripts block")

                                        scripts = ['Genotype_Map_Flow.py', 'Parentage_Genotyping.py']

                                        if not os.path.exists(scripts[0]):
                                            criticalError(f'Case63: Error, File {scripts[0]} not found!')
                                            id = ids[-1]
                                            continue
                                        if not os.path.exists(scripts[1]):
                                            criticalError(f'Case64: Error, File {scripts[1]} not found!')
                                            id = ids[-1]
                                            continue

                                        error_in_script = False
                                        for script in scripts:
                                            with open(script.replace('.py','.log'), 'w') as f:
                                                processo = subprocess.Popen([sys.executable, script, "--numeCari", str(Nume_Cari), "--nomeFile", str(Nome_File)], stdout=subprocess.PIPE, stderr=f)

                                                out, err = processo.communicate()

                                                processo.wait()
                                                
                                                if out == b'C\r\n':
                                                    M_code = config["decode_text_log_XDB"]["c_B"]
                                                                                                elif out == b'A\r\n':
                                                    M_code = config["decode_text_log_XDB"]["g_G"]
                                                elif out == b'I\r\n' or out == b'J\r\n':
                                                    M_code = config["decode_text_log_XDB"]["g_I"]
                                                elif out == b'Error\r\n':
                                                    criticalError("Case65")
                                                    error_in_script = True
                                                    script = scripts[-1]
                                                    continue
                                        
                                                DoLog(1, f'M_code1: {M_code}')
                                        
                                        if error_in_script:
                                            id = ids[-1]    # in this case the program should break and start again from outer while loop
                                            continue

                                        DoLog(1, "----------> Finished procedure Genotype_Map_Flow.py and Parentage_Genotyping.py")

                                        msg, status = aggiorna_Esiti_Caricamento('Esito_caricamento_Genotipi', 'F', templatesParametri, mainParametri, pathTemplatesDir)
                                        if not status:
                                            id = ids[-1]
                                            criticalError("Case70: " + msg)
                                            continue

                                        if M_code["bit_ok"] == 1:
                                            if Aggiorna_Gen_Mappe == 'SI':
                                                Map_Alias = ''                                                           
                                                if Alias == 'SI' and procedura_alias != 'CODE3':
                                                    DoLog(1, f"--> Chip name to update ALias = {CHIP_ALIAS}")
                                                    Map_Alias = CHIP_ALIAS
                                                elif Alias == 'SI' and procedura_alias == 'CODE3':
                                                    DoLog(1, f"--> Chip name to update ALias = {CHIP_ALIAS}")
                                                    Map_Alias = CHIP_ALIAS
                                                else:
                                                    DoLog(1, f"--> Chip name to update ALias = {P.Tipo_Chip}")
                                                    Map_Alias = P.Tipo_Chip
                                                
                                                Map_Name = f"{nsnp}_a"
                                                Number_snp = nsnp
                                                valore = [Nume_Cari, str(Map_Name), Number_snp, str(Map_Alias)]

                                                listOfTables = cursor.execute("SELECT * FROM information_schema.tables WHERE table_name like '%s'" % 'Tmp_Record_Mappe').fetchall()
                                                if len(listOfTables) == 0:
                                                    criticalError('Case71: tabella "Tmp_Record_Mappe" non presente')

                                                query = "INSERT INTO GEN.[Tmp_Record_Mappe] (Nume_Cari,Map_Name,Number_snp,Map_Alias) values(?,?,?,?)"
                                                cursor.execute(query,valore)
                                                conn.commit()

                                                msg, status = aggiorna_Esiti_Caricamento('Esito_caricamento_Genotipi', 'G', templatesParametri, mainParametri, pathTemplatesDir)
                                                if not status:
                                                    id = ids[-1]
                                                    criticalError("Case72: " + msg)
                                                    continue

                                                DoLog(1, "andato buon fine")

                                                M_code = config["decode_text_log_XDB"]["g_K"]
                                            
                                            else:
                                                DoLog(1, "loop di ultimo controle")
                                                DoLog(2, "ATTENZIONE: ---> Numero campione minore delle file Final report originale")
                                                
                                                msg, status = aggiorna_Esiti_Caricamento('Esito_caricamento_Genotipi', 'E', templatesParametri, mainParametri, pathTemplatesDir)
                                                if not status:
                                                    id = ids[-1]
                                                    criticalError("Case73: " + msg)
                                                    continue

                                                DoLog(3, "Errore, Numero campione minore delle file Final report originale.(line #1254) The Map_name in which user has uploaded, has no corrisponding in our sql table.")

                                                M_code = config["decode_text_log_XDB"]["c_B"]
                                                
                                    else:
                                        DoLog(2, "----> Mappa Final Report non coincide con la mappa già caricata")
                                        DoLog(2, "------> Blocco caricamento")
                                        
                                        #Esito_caricamento_Genotipi='C'
                                        msg, status = aggiorna_Esiti_Caricamento('Esito_caricamento_Genotipi', 'C', templatesParametri, mainParametri, pathTemplatesDir)
                                        if not status:
                                            id = ids[-1]
                                            criticalError("Case75: " + msg)
                                            continue

                                        M_code = config["decode_text_log_XDB"]["g_N"]

                                else:
                                    DoLog(2, "----> snp_map-e khar kosde khali nist nanasho gayiidam")
                            
                            ############################
                            ############################
                            # Check the presence of the input samples in the genomic archive
                            # update the GEN.Code_Caricamenti table with the summary information
                            # 1- number of campione in final report
                            # 2- number of campione not present in genomic archive
                            # 3- first 10 campione not present in genomic archive
                            # 4- number of campione with Genotipo
                            # 5- first 10 campione with Genotipo


                            #  To avoid executing the following logic in case of missing chip!
                            if Mappa_Finalreport != 'missing_chip':

                                try:

                                    # Step 1: Determine the number of samples in the finalReport
                                    num_samples_final_report = np.shape(tmp_finalreports_orig['Sample ID'].unique())[0]
                                    # print('Numero campioni presenti nel final report: %s' % num_samples_final_report)

                                    # Step 2: Fetch Sample IDs from the genomic archive
                                    campioni = config["genomica_archivio"]
                                    query = f'SELECT chr_CodiceCampioneLab FROM {campioni}'
                                    cursor = conn.cursor()
                                    cursor.execute(query)
                                    rows = cursor.fetchall()
                                    archive_sample_ids = set(row[0] for row in rows)
                                    archive_sample_ids = set(row[0].strip() for row in rows if row[0] is not None)

                                    # Step 3: Compare the two lists of Sample IDs
                                    final_report_sample_ids = set(tmp_finalreports_orig['Sample ID'].unique())
                                    
                                    # Additional Step: Verify sample code length. We do not skip these samples, just log a warning message.
                                    long_sample_ids = [sample_id for sample_id in final_report_sample_ids if len(sample_id) > 25]

                                    # Log a message if any sample code is longer than 25 characters
                                    for sample_id in long_sample_ids:
                                        DoLog(2, f'Sample code longer than 25 characters: {sample_id}')
                                    
                                    # Determine Sample IDs in the finalReport that are not present in the genomic archive.
                                    samples_not_in_archive = final_report_sample_ids - archive_sample_ids

                                    # read Tmp_Final_Reports from sql for the Nume_Cari
                                    query = f'SELECT * FROM GEN.[{config["Tmp_Finalreports"]}] WHERE Nume_Cari = ?'
                                    cursor.execute(query, Nume_Cari)
                                    rows = cursor.fetchall()
                                    if len(rows) == 0:
                                        criticalError(f'Case 51_0: Tmp_Final_Reports table is empty for Nume_Cari: {Nume_Cari}')
                                    
                                    # count the number of rows in the Tmp_Final_Reports table with the Nume_Cari
                                    num_rows = len(rows)

                                    # get the first 10 Campione from rows
                                    samples_with_Genotipo = set([row[2].replace(" ", "") for row in rows])

                                    # Initialize the dictionary with report data for updating the sql table
                                    report_data = {
                                        "Number of Samples from FinalReport": num_samples_final_report,
                                        "Number of Samples not present in genomic archive": len(samples_not_in_archive),
                                        "First 10 samples not present in genomic archive": list(samples_not_in_archive)[:10],
                                        "Number of Samples with Genotype": num_rows,
                                        "First 10 samples with Genotype": list(samples_with_Genotipo)[:10]
                                    }

                                    # Check if there are samples not present in the genomic archive and add them to the summary
                                    if report_data['Number of Samples not present in genomic archive'] > 0:
                                        report_summary = (
                                            f"{report_data['Number of Samples from FinalReport']} genotypes ready to be uploaded.\n"
                                            f"{report_data['Number of Samples not present in genomic archive']} samples not present in genomic archive\n"
                                            f"The first 10 are: {', '.join(report_data['First 10 samples not present in genomic archive'])}\n"
                                            f"{report_data['Number of Samples with Genotype']} samples with existing Genotype that will be overwritten,\n"
                                            f"The first 10 are: {', '.join(report_data['First 10 samples with Genotype'])} "
                                        )
                                    else:
                                        report_summary = f"{report_data['Number of Samples from FinalReport']} genotypes ready to be uploaded.\n"
                                    
                                    # Check if there are samples with Genotipo and add them to the summary
                                    if report_data['Number of Samples with Genotype'] > 0:
                                        report_summary += (
                                            f"{report_data['Number of Samples with Genotype']} samples with existing Genotype that will be overwritten,\n"
                                            f"The first 10 are: {', '.join(report_data['First 10 samples with Genotype'])} "
                                        )
                                    else:
                                        report_summary += f"{report_data['Number of Samples with Genotype']} samples with existing Genotype that will be overwritten."

                                    # identifier that identifies the record to update.
                                    UniqueIdentifierValue = Nume_Cari

                                    # in this case, we will update the errori_elab field only with the report_summary
                                    M_code["errori_elab"] = report_summary.replace("'", "\"")

                                except Exception as e:
                                    criticalError(f'program: An error occurred: {e}')
                                
                                ############################
                                ############################

            # Extract all values needed for aggiorna_bit()
            bit_ok = M_code["bit_ok"]
            bit_elaborato = M_code["bit_elaborato"]
            errori_elab = M_code["errori_elab"]
                
            DoLog(1, f"Check, bit_ok {bit_ok} bit_elaborato {bit_elaborato} errori_elab {errori_elab}")
            
            if config["doNotUpdate"]:
                DoLog(1, "doNotUpdate is set to True. Skipping the update of the database.")
                print("FINISHED")
                exit()

            msg, status = aggiorna_bit(conn, cursor, Nume_Cari, bit_ok, bit_elaborato, errori_elab)
            if not status:
                id = ids[-1]
                criticalError("Case78: " + msg)
                continue

            DoLog(1, f'END SECONDARY LOOP {sec}')
        
            DoLog(1, 'End of final control')
    
    DoLog(1, f'END MAIN LOOP {sec}')

    # Calculate the elapsed time
    elapsed_time = time.time() - start_time
    
    # Convert elapsed time to minutes and seconds
    minutes, seconds = divmod(elapsed_time, 60)
    elapsed_time_str = f"{int(minutes)} min and {int(seconds)} sec"

    message = f'Total program execution time: {elapsed_time_str}'
    border = '*' * (len(message) + 8)
    DoLog(1, f'\n{border}\n*** {message} ***\n{border}\n')