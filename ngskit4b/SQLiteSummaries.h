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

const int cMaxTextBuff = 4000;			// buffer upto this length results values text
const int cMaxNameLen = 50;				// maximum experiment/process/groupas name length
const int cMaxDescrText = 1000;			// max length descriptive text
const int cMaxSQLStatement = 1000;		// allow for this length SQL statements

typedef enum eSQLliteSummParamTypes {
	 ePTBool = 0,
	 ePTInt32,
	 ePTUint32,
	 ePTInt64,
	 ePTUint64,
	 ePTDouble,
	 ePTText
} teSQLliteSummParamTypes;

typedef struct TAG_sSummStmsSQL {
	char *pTblName;					// table name
	char *pszCreateTbl;				// SQL statement used to create the table
	char *pszInsert;				// SQL statement used to insert row into table
	sqlite3_stmt *pPrepInsert;		// prepared insert row statement; nullptr if not prepared
	char *pszOpenCreateIndexes;		// SQL statement used to create indexes on this table immediately after opening this database, essentially these are used whilst populating
	char *pszCreateIndexes;			// SQL statement used to create indexes on this table before closing the database
	char *pszDropIndexes;			// SQL statement used to drop indexes on this table
} tsSummStmSQL;

typedef struct TAG_sMonoSNP {
	int ExprID;				// SNPs identified within this experiment
	int ProcessingID;		// identifies processing instance
	int MonoSnpPID;			// SNP instance, processing instance unique
	char szElType[51];		// SNP type
	char szSpecies[81];		// SNP located for alignments againts this target/species assembly	
	char szChrom[81];		// SNP is on this chrom
	uint32_t StartLoci;		// offset (0..N) at which SNP located
	uint32_t EndLoci;			// offset (0..N) at which SNP located - allowing for future polymorphic varation covering multiple bases
	uint32_t Len;				// polymorphic variation is of this length
	char szStrand[2];		// SNP relative to this strand
	uint32_t Rank;			// ranking confidence in thisSNP - min 0, max 1000
	double PValue;			// probability of this SNP being a false positive
	uint32_t Bases;			// total number of bases aligning over the SNP loci
	uint32_t Mismatches;		// aligned bases were aligning with this many mismatches
	char szRefBase[2];			// target sequence base at the SNP locai
	uint32_t MMBaseA;			// this many mismatched bases were A
	uint32_t MMBaseC;			// this many mismatched bases were C
	uint32_t MMBaseG;			// this many mismatched bases were G
	uint32_t MMBaseT;			// this many mismatched bases were T
	uint32_t MMBaseN;			// this many mismatched bases were N
	double BackgroundSubRate; // background substitution rate within a window centered at SNP loci
	uint32_t TotWinBases;		// total number of bases within centeredwindow
	uint32_t TotWinMismatches;		// total number of mismatched bases within centered window
	uint32_t	MarkerID;			// marker identifier
	uint32_t NumPolymorphicSites;	// number of polymorphic sites within marker
} tsMonoSNP;

typedef struct TAG_sDiSNP {
	int ExprID;				// SNPs identified within this experiment
	int ProcessingID;		// identifies processing instance
	uint32_t DiSnpPID;		// DiSNP instance, processing instance unique
	char szElType[51];		// SNP type
	char szSpecies[81];		// DiSNP located for alignments againts this target/species assembly	
	char szChrom[81];		// DiSNP is on this chrom
	uint32_t SNP1Loci;		// offset (0..N) at which 1st SNP located
	char szSNP1RefBase[2];		// target sequence base at the 1st SNP loci
	uint32_t SNP1BaseAcnt;	// this many bases at SNP1 were A
	uint32_t SNP1BaseCcnt;	// this many bases at SNP1 were C
	uint32_t SNP1BaseGcnt;	// this many bases at SNP1 were G
	uint32_t SNP1BaseTcnt;	// this many bases at SNP1 were T
	uint32_t SNP1BaseNcnt;	// this many bases at SNP1 were N

	uint32_t SNP2Loci;		// offset (0..N) at which 2nd SNP located
	char szSNP2RefBase[2];		// target sequence base at the 2nd SNP loci
	uint32_t SNP2BaseAcnt;	// this many bases at SNP2 were A
	uint32_t SNP2BaseCcnt;	// this many bases at SNP2 were C
	uint32_t SNP2BaseGcnt;	// this many bases at SNP2 were G
	uint32_t SNP2BaseTcnt;	// this many bases at SNP2 were T
	uint32_t SNP2BaseNcnt;	// this many bases at SNP2 were N

	uint32_t Depth;			// coverage depth of reads covering both SNP1 and SNP2
	uint32_t Antisense;		// non-zero if reads were all antisense to target
	uint32_t Haplotypes;		// minimum number of haplotypes
	uint32_t HaplotypeCnts[16];	// number of instances of each (aa..tt) localised DiSNP haplotype
} tsDiSNP;

typedef struct TAG_sTriSNP {
	int ExprID;				// SNPs identified within this experiment
	int ProcessingID;		// identifies processing instance
	uint32_t TriSnpPID;		// DiSNP instance, processing instance unique
	char szElType[51];		// SNP type
	char szSpecies[81];		// DiSNP located for alignments againts this target/species assembly	
	char szChrom[81];		// DiSNP is on this chrom
	uint32_t SNP1Loci;		// offset (0..N) at which 1st SNP located
	char szSNP1RefBase[2];		// target sequence base at the 1st SNP loci
	uint32_t SNP1BaseAcnt;	// this many bases at SNP1 were A
	uint32_t SNP1BaseCcnt;	// this many bases at SNP1 were C
	uint32_t SNP1BaseGcnt;	// this many bases at SNP1 were G
	uint32_t SNP1BaseTcnt;	// this many bases at SNP1 were T
	uint32_t SNP1BaseNcnt;	// this many bases at SNP1 were N

	uint32_t SNP2Loci;		// offset (0..N) at which 2nd SNP located
	char szSNP2RefBase[2];		// target sequence base at the 2nd SNP loci
	uint32_t SNP2BaseAcnt;	// this many bases at SNP2 were A
	uint32_t SNP2BaseCcnt;	// this many bases at SNP2 were C
	uint32_t SNP2BaseGcnt;	// this many bases at SNP2 were G
	uint32_t SNP2BaseTcnt;	// this many bases at SNP2 were T
	uint32_t SNP2BaseNcnt;	// this many bases at SNP2 were N

	uint32_t SNP3Loci;		// offset (0..N) at which 3rd SNP located
	char szSNP3RefBase[2];		// target sequence base at the 3rd SNP loci
	uint32_t SNP3BaseAcnt;	// this many bases at SNP3 were A
	uint32_t SNP3BaseCcnt;	// this many bases at SNP3 were C
	uint32_t SNP3BaseGcnt;	// this many bases at SNP3 were G
	uint32_t SNP3BaseTcnt;	// this many bases at SNP3 were T
	uint32_t SNP3BaseNcnt;	// this many bases at SNP3 were N

	uint32_t Depth;			// coverage depth of reads covering both SNP1 and SNP2
	uint32_t Antisense;		// non-zero if reads were all antisense to target
	uint32_t Haplotypes;		// minimum number of haplotypes
	uint32_t HaplotypeCnts[64];	// number of instances of each (aaa..ttt) localised TriSNP haplotype
} tsTriSNP;


class CSQLiteSummaries
{
	bool m_bInitialised;				// set true after successful initialisation during class construction
	sqlite3 *m_pDB;						// pts to instance of SQLite
	static tsSummStmSQL m_StmSQL[10];	// SQLite table and index statements

	
	char *RemoveQuotes(char *pszRawText);	// remove quotes which may be a risk for aql injection attacks...

	int					// length - strlen() - of timestamp string 
		GetTimeStamp(char *pszTimeStamp);	// copy timestamp into this

	int Init(void);							// initialisation performed during class construction

	sqlite3 *OpenDatabase(char *pszDatabase,		// database to open, if not already existing then will be created
							bool bReplace = false);			// if database already exists then replace
	int CloseDatabase(void);

	static int ExecCallbackID(void *pCallP1, // callback function processing identifier (4th arg to sqlite3_exec())
					int NumCols,			// number of result columns 
					char **ppColValues,		// array of ptrs to column values 
					char **ppColName);		// array of ptrs to column names

	int							// eBSFSuccess if pszValueText contains textified value
		ValueToText(teSQLliteSummParamTypes ParamType,	// result value type
			int ValueSize,			// result value is of this byte length
			void *pValue,			// result value
			int MaxTextLen,			// truncate returned value text at this many chars - 1 (allows for terminating '\0') NOTE: error if insufficent length to hold max value for type
			char *pszValueText);		// user allocated to hold returned value as text
#ifdef _WIN32
	CRITICAL_SECTION m_hSCritSect;	// used to serialise access to SQLite functionality
#else
	pthread_spinlock_t m_hSpinLock;	// used to serialise access to SQLite functionality
#endif

	inline void SQLiteSerialise(void);	// serialise access to SQLite functionality
	inline void SQLiteRelease(void);

public:
	CSQLiteSummaries(void);
	~CSQLiteSummaries(void);

	int													// returned experiment identifier
			StartExperiment(char *pszDatabase,			// summary results to this SQLite database
							bool bReplace,				// if false then append to existing, if true then replace any existing database
							bool bContinue,				// if true then treat as if continuing pre-existing experiment
							char *pszExprimentName,		// experiment name
							char *pszExperimentTitle,	// experiment title
							char *pszExperimentDescr);  // describes experiment

	int					// returned process identifier 
			AddProcess(char *pszProcessName,			// process name
							char *pszProcessTitle,		// process title
							char *pszProcessDescr);		// describes process


	int													// uiniqely identifies this starting experiment process instance
		    StartProcessing(int ExprID,			// identifier returned by StartExperiment()
						 int ProcessID,					// identifier as returned by AddProcess()
						char *pszProcessVersion);		// process version

	int													// uiniqely identifies this added log text
			AddLog(int ExprID,			// identifier returned by StartExperiment()
					int ProcessingID,					// identifier returned by StartProcessing()
						 const char *pszFormat,...);	// printf style format
	int
			AddParameter(int ExprID,			// identifier returned by StartExperiment()
						int ProcessingID,		// identifier returned by StartProcessing()
						 teSQLliteSummParamTypes ParamType,	// parameter type
						 int ValueSize,					// parameter value is of this byte length
						 const char *pszParamName,		// parameter name
						 void *pParmValue);				// parameter value
	
	int
		AddMonoSNP(int ExprID,			// identifier returned by StartExperiment()
					int ProcessingID,			// identifier returned by StartProcessing()
					tsMonoSNP *pMonoSNP);		// add this MonoSNP to TblMonoSNPs table

	int
		AddDiSNP(int ExprID,			// identifier returned by StartExperiment()
				int ProcessingID,				// identifier returned by StartProcessing()
					tsDiSNP *pDiSNP);		// add this DiSNP to TblDiSNPs table

	int
		AddTriSNP(int ExprID,			// identifier returned by StartExperiment()
					int ProcessingID,				// identifier returned by StartProcessing()
					tsTriSNP *pTriSNP);		// add this TriSNP to TblTriSNPs table

	int
			AddResult(int ExprID,			// identifier returned by StartExperiment()
						int ProcessingID,		// identifier returned by StartProcessing()
						const char *GroupAs,			// result is part of this grouping
						 teSQLliteSummParamTypes ParamType,	// result value type
						 int ValueSize,					// result value is of this byte length
						 const char *pszResultName,		// result name
						 void *pResultValue);			// result value

	int
			AddResultXY(int ExprID,			// identifier returned by StartExperiment()
						int ProcessingID,			// identifier returned by StartProcessing()
						 const char *GroupAs,			// result is part of this grouping
						 teSQLliteSummParamTypes ParamXType,	// result value X type
						 int ValueXSize,				// result value is of this byte length
						 const char *pszResultXName,	// result X name
						 void *pResultXValue,			// result X value
 						 teSQLliteSummParamTypes ParamYType,	// result value Y type
						 int ValueYSize,				// result value is of this byte length
						 const char *pszResultYName,	// result Y name
						 void *pResultYValue);			// result Y value

	int					// returned process identifier
			EndProcessing(int ExprID,			// identifier returned by StartExperiment()
						  int ProcessingID,		// identifier returned by StartProcessing()
						  int ResultCode);		// processing result code

	int
			EndExperiment(int ExprID);			// identifier returned by StartExperiment()
};

