# Genomic Archive ETL Pipeline

## Overview

The **Genomic Archive** project is a high-performance, multi-threaded Python application designed to implement a robust ETL (Extract, Transform, Load) pipeline for processing and analyzing genomic data. This application continuously monitors and processes new records from a SQL Server database, handling both SNP map and Genotype file uploads with advanced validation, transformation, and analytical mechanisms.

The project not only uploads genomic data but also performs complex calculations and analyses, including pedigree verification and SNP genotyping. It is designed to support high-throughput genomic data processing, ensuring data integrity and providing insights into genetic relationships.

## Technical Description

The Genomic Archive project focuses on the efficient processing and management of genomic data, specifically Single Nucleotide Polymorphisms (SNPs) and genotype information. The ETL pipeline is designed to:

1. **Extract** genomic data from a SQL Server database, where it monitors a designated table for new records.
2. **Transform** the data by validating and processing SNP maps and genotype files. This includes:
   - Checking for data integrity, such as duplicate SNP identifiers and missing values.
   - Performing complex transformations to ensure that the data adheres to the expected formats and standards.
   - Utilizing multi-threading to handle large datasets efficiently, allowing for parallel processing of genomic files.
3. **Load** the processed data back into the database, updating the relevant tables with the new information while maintaining data consistency through transaction management.
4. **Analyze** the genomic data to derive meaningful insights, including:
   - Performing pedigree analysis to verify parentage through SNP genotyping.
   - Calculating allele frequencies and call rates for SNPs, providing quality metrics for genetic data.

The application employs a modular architecture, allowing for easy extension and maintenance. It utilizes external JSON configuration files for flexible parameter management, enabling dynamic updates to the execution parameters without modifying the codebase. The logging system provides detailed insights into the processing steps, facilitating troubleshooting and auditing.

## Key Features

- **Asynchronous Database Monitoring**: Utilizes `pyodbc` for efficient, non-blocking database connections and implements a custom max_seconds generator for precise execution time control.
  
- **Sophisticated File Processing**: Handles compressed genomic data files with robust error management and implements intelligent file type detection and parsing for SNP and Genotype data.

- **Advanced Data Validation**: Performs multi-stage validation checks on SNP maps and Genotype data, including data integrity checks such as duplicate detection and missing value handling.

- **Dynamic SQL Operations**: Executes parameterized SQL queries for optimal performance and security, with transaction management to ensure data consistency.

- **Scalable Error Handling and Logging**: Utilizes a hierarchical logging system with configurable verbosity levels and implements a sophisticated error classification and handling mechanism.

- **Configurable Execution Parameters**: Uses external JSON configuration for enhanced flexibility and maintainability, allowing dynamic runtime parameter updates.

- **Optimized Data Structures**: Leverages `pandas` DataFrames for efficient in-memory data manipulation and implements custom data structures for optimized lookup operations.

- **Modular Architecture**: Utilizes separate modules for parameter management and implements a plugin-based system for easy extension of file processing capabilities.

- **Genetic Analysis Capabilities**: Implements algorithms for parentage verification and SNP genotyping, providing insights into genetic relationships and data quality.

## Technical Specifications

- **Language**: Python 3.7+
- **Database**: SQL Server with `pyodbc` driver
- **Data Processing**: `NumPy` and `Pandas` for high-performance numerical computations
- **File Handling**: Custom implementations using `zipfile` for compressed file operations
- **Configuration**: JSON-based configuration management
- **Logging**: Utilizes Python's logging module with custom formatters and handlers

## Requirements

- Python 3.7+
- Libraries: `sys`, `os`, `subprocess`, `numpy`, `pandas`, `pyodbc`, `zipfile`, `json`, `logging`
- SQL Server database with appropriate schema and tables
- `config.json` file with database connection string and runtime parameters
- `Parametri.py` file for custom parameter management

## Contact

For any inquiries, please contact me at [rezza.mozafari@gmail.com](mailto:rezza.mozafari@gmail.com).