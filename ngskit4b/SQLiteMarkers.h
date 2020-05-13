/*
This toolkit is a source base clone of 'BioKanga' release 4.4.2 (https://github.com/csiro-crop-informatics/biokanga) and contains
significant source code changes enabling new functionality and resulting process parameterisation changes. These changes have resulted in
incompatibility with 'BioKanga'.

Because of the potential for confusion by users unaware of functionality and process parameterisation changes then the modified source base
and resultant compiled executables have been renamed to 'kit4b' - K-mer Informed Toolkit for Bioinformatics.
The renaming will force users of the 'BioKanga' toolkit to examine scripting which is dependent on existing 'BioKanga'
parameterisations so as to make appropriate changes if wishing to utilise 'kit4b' parameterisations and functionality.

'kit4b' is being released under the Opensource Software License Agreement (GPLv3)
'kit4b' is Copyright (c) 2019, 2020
Please contact Dr Stuart Stephen < stuartjs@g3web.com > if you have any questions regarding 'kit4b'.

Original 'BioKanga' copyright notice has been retained and immediately follows this notice..
*/
/*
 * CSIRO Open Source Software License Agreement (GPLv3)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENSE for the complete license information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */

#pragma once

const int cMaxIdntNameLen = 50;		// experiment identifiers and other names can be if this maximal length
const int cMaxSeqNameLen = 50;		// targeted sequence names can be if this maximal length
const int cMaxIdntDescrLen = 200;	// allow experiment conditions, and other descriptive fields, to be described by text of this maximal length
const int cMaxExprCultivars = 1000;	// can process markers from CSV file which were generated from at most this many cultivars

const int cMaxMRASeqs = 100;		// cache the last 100 sequence identifiers

typedef struct TAG_sCultivar {
	char szCultivarName[cMaxIdntNameLen+1];		// cultivar short name
	int CultIdx;					// when parsing CSV markers, cultivar starts at this CSV field
	int CultID;						// SQLite assigned cultivar identifier
	int SnpID;						// if SNP in this cultivar at current loci
	int MarkerID;					// if marker in this cultivar at current loci
	} tsCultivar;

typedef struct TAG_sMRASeq {
	char szSeqName[cMaxSeqNameLen+1];	// sequence name
	int SeqID;							// SQLite allocated sequence identifier
	int AccessCnt;						// number of times recently accessed
	} tsMRASeq;


class CSQLiteMarkers
{
	sqlite3 *m_pDB;						// pts to instance of SQLite
	int m_NumSeqMRA;					// number of entries in MRA sequence table
	static tsStmSQL m_StmSQL[7];		// SQLite table and index statements
	tsCultivar Cultivars[cMaxExprCultivars];	// can process upto this many cultivars in CSV marker file
	tsMRASeq m_MRASeqs[cMaxMRASeqs];	// MRA sequences

	bool m_bSafe;						// true if safe select required rather than simply getting last assigned ROWID
	int m_NumSeqs;						// number of seqs added to TblSeqs
	int m_NumSNPLoci;					// number of SNP loci added to TblLoci
	int m_NumSNPs;						// number of SNPs added to TblSnps
	int m_NumMarkers;					// number of markers add to TblMarkers

	sqlite3 *
		CreateDatabase(bool bSafe,		// true if sourcing from input CSV of unknown origin which may contain duplicates etc..
				char *pszDatabase);		// database to create (any existing database is deleted then clean created)

	int
		CloseDatabase(bool bNoIndexes = false);

	int												// errors if < eBSFSuccess, if positive then the ExprID
		CreateExperiment(int CSVtype,					// 0 if markers, 1 if SNPs
				char *pszInFile,				// parse from this input CSV file
				char *pszExprName,				// experiment identifier
				char *pszExprDescr,				// describes experiment
				char *pszAssembName);			// targeted assembly
	

	int												// errors if < eBSFSuccess
		CreateCultivars(int NumCultivars,				// number of cultivars to add
				tsCultivar *pCultivars);			// pts to array of cultivars

	int									// returned sequence identifier for sequence
		AddSeq(int ExprID,						// experiment identifier
		char *pszSeqName);				// target assembly sequence name

	int				// returned loci identifier 
		AddLoci(int ExprID,			// experiment identifier
				int SeqID,			// target assembly sequence identifier
				int Offset,			// offset on sequence			
				char Base);			// cannonical base at loci				

	int							// returned SNP indentifier
		AddSNP(	int ExprID,			// experiment identifier
				int CultID,			// SNP in this cultivar relative to target sequence
				int LociID,			// identifies target sequence loci at which SNP has been identified
                char SrcCnts,       // from where the counts were derived - 'S' SNP call, 'I' imputed from coverage or from SAM/BAM alignment sequences
				int Acnt,			// count of A's in reads covering loci
				int Ccnt,			// count of C's in reads covering loci
				int Gcnt,			// count of G's in reads covering loci
				int Tcnt,			// count of T's in reads covering loci
				int Ncnt,			// count of N's in reads covering loci
				int TotCovCnt,		// total count of reads covering loci
				int TotMMCnt);		// of which there were this many reads with mismatches
	

	int
		AddMarkerSnp(int ExprID,
				int MarkerID,
				int SnpID);
	
	int						// returned marker identifier
		AddMarker(int ExprID,		// experiment identifier
				int CultID,			// marker in this cultivar relative to other cultivars
				int LociID,			// identifies target sequence loci at which marker has been identified
				char MarkerBase,		// marker base
				int MarkerScore);	// score

	static char *RemoveQuotes(char *pszRawText);
				
	static int ExecCallbackCultID(void *pCallP1,	// callback function processing identifier (4th arg to sqlite3_exec())
					int NumCols,			// number of result columns 
					char **ppColValues,		// array of ptrs to column values 
					char **ppColName);		// array of ptrs to column names

	static int ExecCallbackID(void *pCallP1, // callback function processing identifier (4th arg to sqlite3_exec())
					int NumCols,			// number of result columns 
					char **ppColValues,		// array of ptrs to column values 
					char **ppColName);		// array of ptrs to column names

public:
	CSQLiteMarkers(void);
	~CSQLiteMarkers(void);
	int
		ProcessCSV2SQLite(int PMode,			// currently just the one mode...default is to parse from CSV and create/populate SQLite database
				  bool bSafe,					// if true then use indexing on all tables whilst inserting... much slower but perhaps safer if multiple threads ...
			      int CSVtype,					// input CSV file has this format (0: markers, 1: SNPs)
				  char *pszExprName,			// name by which this experiment is identified
				  char *pszExprDescr,			// describes experiment
				  char *pTargAssemb,			// assembly against which aligments for SNP discovery
				  int NumSpecies,				// number of species used in alignments
				  char *pszSpeciesNames[],		// names of species
				  char *pszInFile,				// parse from this input CSV file
				  char *pszDatabase);			// SQLite database file
};


