/*
This toolkit is a source base clone of 'BioKanga' release 4.4.2 (https://github.com/csiro-crop-informatics/biokanga) and contains
significant source code changes enabling new functionality and resulting process parameterisation changes. These changes have resulted in
incompatibility with 'BioKanga'.

Because of the potential for confusion by users unaware of functionality and process parameterisation changes then the modified source base
and resultant compiled executables have been renamed to 'kit4b' - K-mer Informed Toolkit for Bioinformatics.
The renaming will force users of the 'BioKanga' toolkit to examine scripting which is dependent on existing 'BioKanga'
parameterisations so as to make appropriate changes if wishing to utilise 'kit4b' parameterisations and functionality.

'kit4b' is being released under the Open-source Software License Agreement (GPLv3)
'kit4b' is Copyright (c) 2019
Please contact Dr Stuart Stephen < stuartjs@g3web.com > if you have any questions regarding 'kit4b'.

Original 'BioKanga' copyright notice has been retained and immediately follows this notice..
*/
/*
 * CSIRO Open Source Software License Agreement (GPLv3)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENSE for the complete license information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */
#include "stdafx.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#if _WIN32
#include <process.h>
#include "../libkit4b/commhdrs.h"
#else
#include <sys/mman.h>
#include <pthread.h>
#include "../libkit4b/commhdrs.h"
#endif

#include "ngskit4b.h"

#include "../libkit4b/bgzf.h"
#include "SQLitePSL.h"
#include "Blitz.h"

int
Process(etBLZPMode PMode,				// processing mode
		int SampleNthRawRead,		// sample every Nth raw read (or read pair) for processing (1..10000)
		char *pszExprName,				// experiment name
		char *pszExprDescr,				// experiment description
		char *pszParams,				// string containing blitz parameters
		bool KMerDist,					// true if K-mer counts distributions to be reported
		etBLZSensitivity Sensitivity,	// sensitivity 0 - standard, 1 - high, 2 - very high, 3 - low sensitivity
		eALStrand AlignStrand,			// align on to watson, crick or both strands of target
		int MismatchScore,				// decrease score by this for each mismatch bp
		int ExactMatchScore,			// increase score by this for each exactly matching bp
		int GapOpenScore,				// decrease score by this for each gap open
		int  CoreLen,					// use this core length as the exactly matching seed length to be 5' and 3' extended
		int  CoreDelta,					// offset cores by this many bp
		int MaxInsertLen,				// SAM output, accept observed insert sizes of at most this (default = 100000)
		int MaxOccKMerDepth,			// maximum depth to explore over-occurring core K-mers
		int  MinPathScore,				// only report alignment paths on any target sequence if the path score is >= this minimum score
		int QueryLenAlignedPct,			// only report alignment paths if the percentage of total aligned bases to the query sequence length is at least this percentage (1..100)
		int  MaxPathsToReport,			// report at most this many alignment paths for any query
		etBLZRsltsFomat RsltsFormat,	// output results format
		char *pszInputFile,				// name of input file containing query sequences (PE1 if PE processing)
		char *pszInputFilePE2,			// name of input file containing PE2 query sequences (only applies if output format is SAM)
		char *pszSfxFile,				// target as suffix array
		char *pszOutFile,				// where to write alignments
		int NumThreads);				// number of worker threads to use


#ifdef _WIN32
int Blitz(int argc, char* argv[])
{
// determine my process name
_splitpath(argv[0],NULL,NULL,gszProcName,NULL);
#else
int
Blitz(int argc, char** argv)
{
// determine my process name
CUtility::splitpath((char *)argv[0],NULL,gszProcName);
#endif

int iFileLogLevel;			// level of file diagnostics
int iScreenLogLevel;		// level of file diagnostics
char szLogFile[_MAX_PATH];	// write diagnostics to this file
int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure

int PMode;					// processing mode
int FMode;					// format output mode

bool KMerDist;				// true if K_mer counts distributions to be reported

int NumberOfProcessors;		// number of installed CPUs
int NumThreads;				// number of threads (0 defaults to number of CPUs)
int SampleNthRawRead;		// sample every Nth raw read (or read pair) for processing (1..10000)
int Sensitivity;			// sensitivity 0 - standard, 1 - high, 2 - very high, 3 - low sensitivity (default is 0)

int MismatchScore;			// decrease score by this for each mismatch bp
int ExactMatchScore;		// increase score by this for each exactly matching bp
int GapOpenScore;			// decrease score by this for each gap open

int CoreLen;				// use this core length as the exactly matching seed length to be 5' and 3' extended whilst no more than m_MaxSubRate
int CoreDelta;				// offset cores by this many bp
int MaxOccKMerDepth;		// maximum depth to explore over-occurring core K-mers
int MaxInsertLen;			// SAM output, accept observed insert sizes of at most this (default = 100000)

int  MinPathScore;			// only report alignment paths on any target sequence if the path score is >= this minimum score
int QueryLenAlignedPct;		// only report alignment paths if the percentage of total aligned bases to the query sequence length is at least this percentage (1..100)
int  MaxPathsToReport;		// report at most this many alignment paths for any query
int AlignStrand;			// align on to watson, crick or both strands of target
char szRsltsFile[_MAX_PATH];			// results to this file
char szTargFile[_MAX_PATH];				// align against this target suffix array genome file

char szInputFile[_MAX_PATH];	// input file containing sequences to be aligned (or PE1 if PE processing with SAM output format)
char szInputFilePE2[_MAX_PATH];	// name of input file containing PE2 query sequences (only applies if output format is SAM)

char szSQLiteDatabase[_MAX_PATH];	// results summaries to this SQLite file
char szExperimentName[cMaxDatasetSpeciesChrom+1];	// experiment name
char szExperimentDescr[1000];		// describes experiment
char szBlitzParams[2000];			// to hold Blitz parameters
//
struct arg_lit  *help    = arg_lit0("h","help",                 "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *pmode = arg_int0("m","mode","<int>",		    "alignment processing mode: 0 - standard");
struct arg_int* samplenthrawread = arg_int0("#", "samplenthrawread", "<int>", "sample every Nth raw read or read pair for processing (default 1, range 1..10000)");
struct arg_int *format = arg_int0("M","format","<int>",		    "output format: 0 - PSL, 1 - PSLX, 2 - MAF, 3 - BED, 4 - SQLite, 5 - SAM (default 0 - PSL)");
struct arg_file *inputfile = arg_file1("i","in","<file>",		"input sequences to align from this file (or PE1 if PE processing when SAM output format mode");
struct arg_file *inputfilepe2 = arg_file0("u", "inpe2", "<file>", "PE2 input sequences to align from this file if PE processing with SAM output format mode");

struct arg_lit  *kmerdist = arg_lit0("K", "kmerdist",			"report on K-mer counts distributions");

struct arg_int  *alignstrand = arg_int0("Q","alignstrand","<int>", "align to this strand: 0 either, 1 Watson '+', 2 Crick '-' (default is to align to either strand)");
struct arg_file *sfxfile = arg_file1("I","sfx","<file>",		"align against this suffix array ('kit4b index' generated) file");
struct arg_file *outfile = arg_file1("o","out","<file>",		"output alignments to this file");

struct arg_int *maxocckmerdepth = arg_int0("k","maxocckmerdepth","<int>",	"maximum depth to explore over-occurring seed K-mers (default is 0 for auto, range 100 to 20000)");
struct arg_int *sensitivity = arg_int0("s","sensitivity","<int>",	"sensitivity 0 - standard, 1 - high, 2 - very high, 3 - low sensitivity (default is 0)");

struct arg_int *mismatchscore = arg_int0("j","mismatchscore","<int>",	"penalise score for bp mismatches (default is 2, range 1..10)");
struct arg_int *exactmatchscore = arg_int0("J","exactmatchscore","<int>",	"score exact bp matching (default is 1, range 1..10)");
struct arg_int *gapopenscore = arg_int0("g","gapopenscore","<int>",	"penalise score for gap openings (default is 5, range 1..50)");

struct arg_int *coredelta = arg_int0("c","coredelta","<int>",	"core (seed) delta (default is 0 for auto, range 1..corelen)");
struct arg_int *corelen = arg_int0("C","corelen","<int>",		"core (seed) length (default is 0 for auto, range 8..100)");

struct arg_int  *maxinsertlen = arg_int0("D", "maxinsertlen", "<int>",		"SAM output, accept observed insert sizes of at most this (default = 100000, range 1000 to 1000000)");

struct arg_int *minpathscore = arg_int0("p","minpathscore","<int>",		"minimum alignment path score (default is 0 for auto, range 25..50000)");
struct arg_int *querylendpct = arg_int0("a","querylendpct","<int>",		"minimum required percentage of query sequence aligned (default is 75, range 20 to 100)");

struct arg_int *maxpathstoreport = arg_int0("P","maxpathstoreport","<int>",	"report at most this many highest scored alignment paths for each query (default is 1)");


struct arg_int *threads = arg_int0("T","threads","<int>",			"number of processing threads 0..128 (defaults to 0 which sets threads to number of CPU cores)");

struct arg_file *summrslts = arg_file0("q","sumrslts","<file>",		"Output results summary to this SQLite3 database file");
struct arg_str *experimentname = arg_str0("w","experimentname","<str>",		"experiment name SQLite3 database file");
struct arg_str *experimentdescr = arg_str0("W","experimentdescr","<str>",	"experiment description SQLite3 database file");

struct arg_end *end = arg_end(200);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					kmerdist,
					summrslts,experimentname,experimentdescr,
					pmode,samplenthrawread,sensitivity,alignstrand,mismatchscore,exactmatchscore,gapopenscore,coredelta,corelen,maxinsertlen,maxocckmerdepth,minpathscore,querylendpct,maxpathstoreport,format,
					inputfile,inputfilepe2,sfxfile,outfile,threads,
					end};

char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

/* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
        {
		printf("\n%s %s %s, Version %s\nOptions ---\n", gszProcName,gpszSubProcess->pszName,gpszSubProcess->pszFullDescr,kit4bversion);
        arg_print_syntax(stdout,argtable,"\n");
        arg_print_glossary(stdout,argtable,"  %-25s %s\n");
		printf("\nNote: Parameters can be entered into a parameter file, one parameter per line.");
		printf("\n      To invoke this parameter file then precede its name with '@'");
		printf("\n      e.g. %s %s @myparams.txt\n",gszProcName,gpszSubProcess->pszName);
		printf("\nPlease report any issues regarding usage of %s at https://github.com/kit4b/issues\n\n",gszProcName);
		return(1);
        }

    /* special case: '--version' takes precedence error reporting */
if (version->count > 0)
        {
		printf("\n%s %s Version %s\n",gszProcName,gpszSubProcess->pszName,kit4bversion);
		return(1);
        }

if (!argerrors)
	{
	if(FileLogLevel->count && !LogFile->count)
		{
		printf("\nError: FileLogLevel '-f%d' specified but no logfile '-F<logfile>\n'",FileLogLevel->ival[0]);
		exit(1);
		}

	iScreenLogLevel = iFileLogLevel = FileLogLevel->count ? FileLogLevel->ival[0] : eDLInfo;
	if(iFileLogLevel < eDLNone || iFileLogLevel > eDLDebug)
		{
		printf("\nError: FileLogLevel '-l%d' specified outside of range %d..%d\n",iFileLogLevel,eDLNone,eDLDebug);
		exit(1);
		}

	if(LogFile->count)
		{
		strncpy(szLogFile,LogFile->filename[0],_MAX_PATH);
		szLogFile[_MAX_PATH-1] = '\0';
		}
	else
		{
		iFileLogLevel = eDLNone;
		szLogFile[0] = '\0';
		}

	// now that log parameters have been parsed then initialise diagnostics log system
	if(!gDiagnostics.Open(szLogFile,(etDiagLevel)iScreenLogLevel,(etDiagLevel)iFileLogLevel,true))
		{
		printf("\nError: Unable to start diagnostics subsystem\n");
		if(szLogFile[0] != '\0')
			printf(" Most likely cause is that logfile '%s' can't be opened/created\n",szLogFile);
		exit(1);
		}

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Subprocess %s Version %s starting",gpszSubProcess->pszName,kit4bversion);
	gExperimentID = 0;
	gProcessID = 0;
	gProcessingID = 0;
	szSQLiteDatabase[0] = '\0';
	szExperimentName[0] = '\0';
	szExperimentDescr[0] = '\0';


	if(experimentname->count)
		{
		strncpy(szExperimentName,experimentname->sval[0],sizeof(szExperimentName));
		szExperimentName[sizeof(szExperimentName)-1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(szExperimentName);
		CUtility::ReduceWhitespace(szExperimentName);
		}

	if(strlen(szExperimentName) < 1)
		strcpy(szExperimentName,"N/A");

	if(experimentdescr->count)
		{
		strncpy(szExperimentDescr,experimentdescr->sval[0],sizeof(szExperimentDescr)-1);
		szExperimentDescr[sizeof(szExperimentDescr)-1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(szExperimentDescr);
		CUtility::ReduceWhitespace(szExperimentDescr);
		}
	if(strlen(szExperimentDescr) < 1)
		strcpy(szExperimentDescr,"N/A");

	gExperimentID = 0;
	gProcessID = 0;
	gProcessingID = 0;
	szSQLiteDatabase[0] = '\0';

	if(summrslts->count)
		{
		strncpy(szSQLiteDatabase,summrslts->filename[0],sizeof(szSQLiteDatabase)-1);
		szSQLiteDatabase[sizeof(szSQLiteDatabase)-1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(szSQLiteDatabase);
		if(strlen(szSQLiteDatabase) < 1)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no SQLite database specified with '-q<filespec>' option");
			return(1);
			}

		gExperimentID = gSQLiteSummaries.StartExperiment(szSQLiteDatabase,false,true,szExperimentName,szExperimentName,szExperimentDescr);
		if(gExperimentID < 1)
			return(1);
		gProcessID = gSQLiteSummaries.AddProcess((char *)gpszSubProcess->pszName,(char *)gpszSubProcess->pszName,(char *)gpszSubProcess->pszFullDescr);
		if(gProcessID < 1)
			return(1);
		gProcessingID = gSQLiteSummaries.StartProcessing(gExperimentID,gProcessID,(char *)kit4bversion);
		if(gProcessingID < 1)
			return(1);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Initialised SQLite database '%s' for results summary collection",szSQLiteDatabase);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"SQLite database experiment identifier for '%s' is %d",szExperimentName,gExperimentID);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"SQLite database process identifier for '%s' is %d",(char *)gpszSubProcess->pszName,gProcessID);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"SQLite database processing instance identifier is %d",gProcessingID);
		}
	else
		szSQLiteDatabase[0] = '\0';

	// ensure all filenames are initialised in case not user specified
	szRsltsFile[0] = '\0';
	szTargFile[0] = '\0';
	szInputFile[0] = '\0';
	szInputFilePE2[0] = '\0';

	PMode = pmode->count ? pmode->ival[0] : (int)eBLZPMdefault;
	if(PMode < eBLZPMdefault || PMode >= eBLZPMplaceholder)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Processing mode '-m%d' specified outside of range %d..%d\n",PMode,eBLZPMdefault,(int)eBLZPMplaceholder-1);
		exit(1);
		}

	SampleNthRawRead = samplenthrawread->count ? samplenthrawread->ival[0] : 1;
	if (SampleNthRawRead < 1)
		SampleNthRawRead = 1;
	else
		if (SampleNthRawRead > 10000)
			SampleNthRawRead = 10000;

	KMerDist = kmerdist->count > 0 ? true : false;

	Sensitivity = sensitivity->count ? sensitivity->ival[0] : (int)eBLZSdefault;
	if(Sensitivity < eBLZSdefault || Sensitivity >= eBLZSplaceholder)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Sensitivity mode '-s%d' specified outside of range %d..%d\n",Sensitivity,eBLZSdefault,(int)eBLZSplaceholder-1);
		exit(1);
		}

	QueryLenAlignedPct = querylendpct->count ? querylendpct->ival[0] : cDfltMinQueryLenAlignedPct;
	if(QueryLenAlignedPct < 1 || QueryLenAlignedPct > 100)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Percentage of query sequence aligned '-a%d' specified outside of range 1..100\n",QueryLenAlignedPct);
		exit(1);
		}

	CoreLen = 0;
	CoreDelta = 0;
	MinPathScore = 0;
	MaxInsertLen = 0;

	AlignStrand = (eALStrand)(alignstrand->count ? alignstrand->ival[0] : eALSboth);
	if(AlignStrand < eALSboth || AlignStrand >= eALSnone)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Aligned to strand '-Q%d' specified outside of range %d..%d\n",AlignStrand,eALSboth,(int)eALSnone-1);
		exit(1);
		}


	FMode = format->count ? format->ival[0] : (int)eBLZRsltsPSL;
	if(FMode < eBLZRsltsPSL || FMode >= eBLZRsltsplaceholder)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Output format mode '-M%d' specified outside of range %d..%d\n",FMode,eBLZRsltsPSL,(int)eBLZRsltsplaceholder-1);
		exit(1);
		}

	CoreLen = corelen->count ?  corelen->ival[0] : CoreLen;
	if(CoreLen != 0 && (CoreLen < cMinCoreLen) || CoreLen > cMaxCoreLen)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: core or seed length '-C%d' specified outside of range %d..%d\n",CoreLen,cMinCoreLen,cMaxCoreLen);
		exit(1);
		}

	CoreDelta = coredelta->count ?  coredelta->ival[0] : CoreDelta;
	if((CoreDelta != 0 && CoreDelta < cMinCoreDelta) || CoreDelta > CoreLen)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: core delta '-c%d' specified outside of range %d..%d\n",CoreDelta,cMinCoreDelta, cMaxCoreLen);
		exit(1);
		}

	if (FMode == eBLZRsltsSAM)
		{
		MaxInsertLen = maxinsertlen->count ? maxinsertlen->ival[0] : 100000;
		if(MaxInsertLen == 0)
			MaxInsertLen = 100000;
		if (MaxInsertLen < 1000 || MaxInsertLen > 1000000)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: maximum Insert length '-D%d' specified outside of range 1000..1000000\n", MaxInsertLen);
			exit(1);
			}
		}

	MismatchScore = mismatchscore->count ?  mismatchscore->ival[0] : cDfltMismatchScore;
	if(MismatchScore < cMinMismatchScore || MismatchScore > cMaxMismatchScore)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: mismatch penalty '-s%d' specified outside of range %d..%d\n",MismatchScore,cMinMismatchScore,cMaxMismatchScore);
		exit(1);
		}
	ExactMatchScore = exactmatchscore->count ?  exactmatchscore->ival[0] : cDfltExactMatchScore;
	if(ExactMatchScore < cMinExactMatchScore || ExactMatchScore > cMaxExactMatchScore)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: exact match score '-s%d' specified outside of range %d..%d\n",ExactMatchScore,cMinExactMatchScore,cMaxExactMatchScore);
		exit(1);
		}
	GapOpenScore = gapopenscore->count ?  gapopenscore->ival[0] : cDfltGapOpenScore;
	if(GapOpenScore < cMinGapOpenScore || GapOpenScore > cMaxGapOpenScore)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: gap open penalty '-s%d' specified outside of range %d..%d\n",GapOpenScore,cMinGapOpenScore,cMaxGapOpenScore);
		exit(1);
		}

	MaxOccKMerDepth = maxocckmerdepth->count ?  maxocckmerdepth->ival[0] : 0;
	if(MaxOccKMerDepth != 0 && (MaxOccKMerDepth < cMinOccKMerDepth) || MaxOccKMerDepth > cMaxOccKMerDepth)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: maximum depth to explore over-occurring seed K-mers '-k%d' specified outside of range %d..%d\n",MaxOccKMerDepth,cMinOccKMerDepth,cMaxOccKMerDepth);
		exit(1);
		}

	MinPathScore = minpathscore->count ?  minpathscore->ival[0] : MinPathScore;
	if(MinPathScore != 0 && (MinPathScore < cMinPathScore) || MinPathScore > cMaxPathScore)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: minimum path score '-p%d' specified outside of range %d..%d\n",MinPathScore,cMinPathScore,cMaxPathScore);
		exit(1);
		}

	MaxPathsToReport = maxpathstoreport->count ?  maxpathstoreport->ival[0] : cDfltMaxPathsToReport;
	if(MaxPathsToReport < 1 || MaxPathsToReport > 10)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: maximum number of highest scoring paths per query '-P%d' must be in range 1..10\n",MaxPathsToReport);
		exit(1);
		}

// show user current resource limits
#ifndef _WIN32
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Resources: %s",CUtility::ReportResourceLimits());
#endif

#ifdef _WIN32
	SYSTEM_INFO SystemInfo;
	GetSystemInfo(&SystemInfo);
	NumberOfProcessors = SystemInfo.dwNumberOfProcessors;
#else
	NumberOfProcessors = sysconf(_SC_NPROCESSORS_CONF);
#endif
	int MaxAllowedThreads = min(cMaxWorkerThreads,NumberOfProcessors);	// limit to be at most cMaxWorkerThreads
	if((NumThreads = threads->count ? threads->ival[0] : MaxAllowedThreads)==0)
		NumThreads = MaxAllowedThreads;
	if(NumThreads < 0 || NumThreads > MaxAllowedThreads)
		{
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"Warning: Number of threads '-T%d' specified was outside of range %d..%d",NumThreads,1,MaxAllowedThreads);
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"Warning: Defaulting number of threads to %d",MaxAllowedThreads);
		NumThreads = MaxAllowedThreads;
		}

	strcpy(szTargFile,sfxfile->filename[0]);
	strcpy(szInputFile,inputfile->filename[0]);
	if(FMode != eBLZRsltsSAM && inputfilepe2->count)
		{
		gDiagnostics.DiagOut(eDLWarn, gszProcName, "Warning: PE2 input file '-u%s' only processed for pairing when output format is SAM", inputfilepe2->filename[0]);
		gDiagnostics.DiagOut(eDLWarn, gszProcName, "Warning: Ignoring PE2 input file");
		}
	if (FMode == eBLZRsltsSAM && inputfilepe2->count)
		strcpy(szInputFilePE2, inputfilepe2->filename[0]);

	strcpy(szRsltsFile,outfile->filename[0]);

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing parameters:");

	const char *pszDescr;
	switch(PMode) {
		case eBLZPMdefault:
		default:
			pszDescr = "Standard alignment processing";
			break;
		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Alignment processing is : '%s'",pszDescr);
	gDiagnostics.DiagOutMsgOnly(eDLInfo, "Raw read or paired reads sampling is : every %u", SampleNthRawRead);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Experiment name : '%s'",szExperimentName);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Experiment description : '%s'",szExperimentDescr);
	gDiagnostics.DiagOutMsgOnly(eDLInfo, "Reporting K-mer counts distributions : '%s'", KMerDist ? "Yes" : "No");
	
	switch(Sensitivity) {
		case eBLZSdefault:
			pszDescr = "Standard alignment sensitivity";
			break;
		case eBLZSMoreSens:
			pszDescr = "High alignment sensitivity";
			break;
		case eBLZSUltraSens:
			pszDescr = "Very high alignment sensitivity - caution: very slow";
			break;
		default:
			pszDescr = "Less sensitive alignment (quicker)";
			break;
		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Sensitivity is : '%s'",pszDescr);
	switch(AlignStrand) {
		case eALSboth:
			pszDescr = "Watson '+' and Crick '-' strands";
			break;
		case eALSWatson:
			pszDescr = "Watson '+' strand only";
			break;
		case eALSCrick:
			pszDescr = "Crick '-' strand only";
			break;
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Mismatch score penalty : %d",MismatchScore);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Exact match score : %d",ExactMatchScore);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Gap open score penalty : %d",GapOpenScore);

	if(CoreLen == 0)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Core length : Auto");
	else
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Core length : %d",CoreLen);
	if(CoreDelta == 0)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Core delta : Auto");
	else
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Core delta : %d",CoreDelta);
	if(FMode == eBLZRsltsSAM)
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "Max insert Length : %d", MaxInsertLen);
	if(MaxOccKMerDepth == 0)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Maximum depth to explore over-occurring seed K-mers : Auto");
	else
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Maximum depth to explore over-occurring seed K-mers : %d",MaxOccKMerDepth);
	if(MinPathScore == 0)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum path score : Auto");
	else
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum path score : %d",MinPathScore);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum percentage of query sequence aligned : %d",QueryLenAlignedPct);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"maximum number of highest scoring paths per query : %d",MaxPathsToReport);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"alignments are to : %s",pszDescr);

	switch(FMode) {
		case eBLZRsltsPSL:
			pszDescr = "PSL";
			break;
		case eBLZRsltsPSLX:		
			pszDescr = "PSLX";
			break;
		case eBLZRsltsMAF:
			pszDescr = "MAF";
			break;
		case eBLZRsltsBED:
			pszDescr = "BED";
			break;
		case eBLZRsltsSQLite:
			pszDescr = "SQLite";
			break;
		case eBLZRsltsSAM:
			pszDescr = "SAM (all query sequences will be truncated to be a maximum length of 16000bp)";
			break;
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"output format is : '%s'",pszDescr);
	if(szInputFilePE2[0] == '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"input query sequences file: '%s'",szInputFile);
	else
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "input query sequences PE1 file: '%s'", szInputFile);
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "input query sequences PE2 file: '%s'", szInputFilePE2);
		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"input target sequence(s) suffix array file: '%s'",szTargFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"output results file: '%s'",szRsltsFile);


	if(szExperimentName[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"This processing reference: %s",szExperimentName);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"number of threads : %d",NumThreads);

	if(gExperimentID > 0)
		{
		int ParamID;
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTText,(int)strlen(szLogFile),"log",szLogFile);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,(int)sizeof(PMode),"mode",&PMode);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID, ePTInt32, (int)sizeof(SampleNthRawRead), "samplenthrawread", &SampleNthRawRead);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,(int)sizeof(KMerDist),"kmerdist",&KMerDist);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID, ePTInt32, (int)sizeof(Sensitivity), "sensitivity", &Sensitivity);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,(int)sizeof(MismatchScore),"mismatchscore",&MismatchScore);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,(int)sizeof(ExactMatchScore),"exactmatchscore",&ExactMatchScore);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,(int)sizeof(GapOpenScore),"gapopenscore",&GapOpenScore);

		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,(int)sizeof(FMode),"format",&FMode);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,(int)sizeof(AlignStrand),"alignstrand",&AlignStrand);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,(int)sizeof(CoreLen),"corelen",&CoreLen);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,(int)sizeof(CoreDelta),"coredelta",&CoreDelta);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID, ePTInt32, (int)sizeof(MaxInsertLen), "maxinsertlen", &MaxInsertLen);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,(int)sizeof(MaxOccKMerDepth),"maxocckmerdepth",&MaxOccKMerDepth);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,(int)sizeof(MinPathScore),"minpathscore",&MinPathScore);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,(int)sizeof(QueryLenAlignedPct),"querylendpct",&QueryLenAlignedPct);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,(int)sizeof(MaxPathsToReport),"maxpathstoreport",&MaxPathsToReport);

		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID, ePTText, (int)strlen(szInputFile), "in", szInputFile);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID, ePTText, (int)strlen(szInputFilePE2), "inpe2", szInputFilePE2);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTText,(int)strlen(szTargFile),"sfx",szTargFile);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTText,(int)strlen(szRsltsFile),"out",szRsltsFile);
		
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,(int)sizeof(NumThreads),"threads",&NumThreads);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,(int)sizeof(NumberOfProcessors),"cpus",&NumberOfProcessors);

		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTText,(int)strlen(szSQLiteDatabase),"sumrslts",szSQLiteDatabase);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTText,(int)strlen(szExperimentName),"experimentname",szExperimentName);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTText,(int)strlen(szExperimentDescr),"experimentdescr",szExperimentDescr);
		}

	sprintf(szBlitzParams,"mode: %d sensitivity: %d mismatchscore: %d exactmatchscore: %d gapopenscore: %d alignstrand: %d corelen: %d coredelta: %d maxinsertlen: %d maxocckmerdepth: %d minpathscore: %d querylendpct: %d maxpathstoreport: %d",
							PMode, Sensitivity, MismatchScore, ExactMatchScore, GapOpenScore,AlignStrand,CoreLen,CoreDelta, MaxInsertLen,MaxOccKMerDepth,MinPathScore,QueryLenAlignedPct,MaxPathsToReport);

#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = Process((etBLZPMode)PMode, SampleNthRawRead,szExperimentName,szExperimentDescr,szBlitzParams, KMerDist,(etBLZSensitivity)Sensitivity,(eALStrand)AlignStrand,MismatchScore,ExactMatchScore,GapOpenScore,CoreLen, CoreDelta, MaxInsertLen,
							MaxOccKMerDepth, MinPathScore,QueryLenAlignedPct,MaxPathsToReport,(etBLZRsltsFomat)FMode,szInputFile,szInputFilePE2,szTargFile,szRsltsFile,NumThreads);
	Rslt = Rslt >=0 ? 0 : 1;
	if(gExperimentID > 0)
		{
		if(gProcessingID)
			gSQLiteSummaries.EndProcessing(gExperimentID, gProcessingID,Rslt);
		gSQLiteSummaries.EndExperiment(gExperimentID);
		}
	gStopWatch.Stop();

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	exit(Rslt);
	}
else
	{
    printf("\n%s %s %s, Version %s\n", gszProcName,gpszSubProcess->pszName,gpszSubProcess->pszFullDescr,kit4bversion);
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nUse '-h' to view option and parameter usage\n");
	exit(1);
	}
return 0;
}


int
Process(etBLZPMode PMode,				// processing mode
		int SampleNthRawRead,		// sample every Nth raw read (or read pair) for processing (1..10000)
		char *pszExprName,				// experiment name
		char *pszExprDescr,				// experiment description
		char *pszParams,				// string containing blitz parameters
		bool KMerDist,				// true if K_mer counts distributions to be reported
		etBLZSensitivity Sensitivity,	// sensitivity 0 - standard, 1 - high, 2 - very high, 3 - low sensitivity
		eALStrand AlignStrand,			// align on to watson, crick or both strands of target
		int MismatchScore,				// decrease score by this for each mismatch bp
		int ExactMatchScore,			// increase score by this for each exactly matching bp
		int GapOpenScore,				// decrease score by this for each gap open
		int  CoreLen,					// use this core length as the exactly matching seed length to be 5' and 3' extended whilst no more than m_MaxSubRate
		int  CoreDelta,					// offset cores by this many bp
		int MaxInsertLen,				// SAM output, accept observed insert sizes of at most this (default = 50000)
		int MaxOccKMerDepth,			// maximum depth to explore over-occurring core K-mers
		int  MinPathScore,				// only report alignment paths on any target sequence if the path score is >= this minimum score
		int QueryLenAlignedPct,				// only report alignment paths if the percentage of total aligned bases to the query sequence length is at least this percentage (1..100)
		int  MaxPathsToReport,			// report at most this many alignment paths for any query
		etBLZRsltsFomat RsltsFormat,	// output results format
		char *pszInputFile,				// name of input file containing query sequences (PE1 if PE processing when output format is SAM)
		char *pszInputFilePE2,			// name of input file containing PE2 query sequences (only applies if output format is SAM)
		char *pszSfxFile,				// target as suffix array
		char *pszOutFile,				// where to write alignments
		int NumThreads)					// number of worker threads to use
{
int Rslt;
CBlitz *pBlitzer;

if((pBlitzer = new CBlitz)==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: Unable to instantiate CAligner");
	return(eBSFerrObj);
	}
Rslt = pBlitzer->Process(PMode, SampleNthRawRead,pszExprName,pszExprDescr,pszParams,KMerDist,Sensitivity,AlignStrand,MismatchScore,ExactMatchScore,GapOpenScore,CoreLen,CoreDelta,MaxInsertLen,MaxOccKMerDepth,
				MinPathScore,QueryLenAlignedPct,MaxPathsToReport,RsltsFormat,pszInputFile,pszInputFilePE2,pszSfxFile,pszOutFile,NumThreads);
delete pBlitzer;
return(Rslt);
}

CBlitz::CBlitz()
{
Init();
}

CBlitz::~CBlitz()
{
Reset(false);
}


void
CBlitz::Init(void)
{
m_hOutFile = -1;
m_hOutNonAlignedFilePE1 = -1;
m_hOutNonAlignedFilePE2 = -1;
m_pKmerOccsDist = NULL;
m_TotSeqIDs = 0;
m_NumQuerySeqs = 0;
m_NxtQuerySeqIdx = 0;
m_AllocdQuerySeqs = 0;
m_pQuerySeqs = NULL;
m_pSQLitePSL=NULL;
m_ProcMode = eBLZPMdefault;
m_SampleNthRawRead = 1;
m_Sensitivity = eBLZSdefault;				
m_AlignStrand = eALSboth;	
m_CoreLen = cDfltCoreLen;
m_CoreDelta = (cDfltCoreLen+1)/2;
m_MaxInsertLen = 0;
m_MaxQuerySeqLen = cMaxQuerySeqLen;
m_MaxOccKMerDepth = cDfltSensCoreIters;
m_QueryLenAlignedPct = cDfltMinQueryLenAlignedPct;
m_MinPathScore = cDfltPathScore;
m_MaxPathsToReport = cDfltMaxPathsToReport;
m_AlignPathID = 0;
m_RsltsFormat = eBLZRsltsPSL;	
m_pszInputFile = NULL;
m_pszInputFilePE2 = NULL;
m_pszSfxFile = NULL;		
m_pszOutFile = NULL;		
m_pSfxArray = NULL;
m_pszLineBuff = NULL;
m_szLineBuffIdx = 0;
m_ReportedPaths = 0;
m_QueriesPaths = 0;
m_NumQueriesProc = 0;
m_NumThreads = 0;		
m_bMutexesCreated = false;
m_TermBackgoundThreads = 0;
}

void
CBlitz::Reset(bool bSync)			// if bSync true then fsync before closing output file handles
{
m_TermBackgoundThreads = 0x01;	// need to require any background threads to self-terminate
if(m_hOutFile != -1)
	{
	if(bSync)
#ifdef _WIN32
		_commit(m_hOutFile);
#else
		fsync(m_hOutFile);
#endif
	close(m_hOutFile);
	m_hOutFile = -1;
	}

if (m_hOutNonAlignedFilePE1 != -1)
	{
	if (bSync)
#ifdef _WIN32
		_commit(m_hOutNonAlignedFilePE1);
#else
		fsync(m_hOutNonAlignedFilePE1);
#endif
	close(m_hOutNonAlignedFilePE1);
	m_hOutNonAlignedFilePE1 = -1;
	}

if (m_hOutNonAlignedFilePE2 != -1)
	{
	if (bSync)
#ifdef _WIN32
		_commit(m_hOutNonAlignedFilePE2);
#else
		fsync(m_hOutNonAlignedFilePE2);
#endif
	close(m_hOutNonAlignedFilePE2);
	m_hOutNonAlignedFilePE2 = -1;
	}

if(m_pKmerOccsDist != NULL)
	{
	delete m_pKmerOccsDist;
	m_pKmerOccsDist = NULL;
	}

if(m_pSQLitePSL != NULL)
	{
	delete m_pSQLitePSL;
	m_pSQLitePSL = NULL;
	}

if(m_pszLineBuff != NULL)
	{
	delete m_pszLineBuff;
	m_pszLineBuff = NULL;
	}

if(m_pSfxArray != NULL)
	{
	delete m_pSfxArray;
	m_pSfxArray = NULL;
	}

if(m_pQuerySeqs != NULL)
	{
	tsQuerySeq *pQuerySeq;
	if(m_NumQuerySeqs)
		{
		do {
			pQuerySeq = &m_pQuerySeqs[m_NxtQuerySeqIdx++];
			if(pQuerySeq->pQuerySeq != NULL)
				delete pQuerySeq->pQuerySeq;
			if(m_NxtQuerySeqIdx == m_AllocdQuerySeqs)
				m_NxtQuerySeqIdx = 0;
			}
		while(m_NumQuerySeqs--);
		}
	delete m_pQuerySeqs;
	m_pQuerySeqs = NULL;
	}

DeleteMutexes();

Init();
m_TermBackgoundThreads = 0x0;	// can startup any background thread processing
}



int
CBlitz::CreateMutexes(void)
{
if(m_bMutexesCreated)
	return(eBSFSuccess);

#ifdef _WIN32
InitializeSRWLock(&m_hRwLock);
#else
if(pthread_rwlock_init (&m_hRwLock,NULL)!=0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to create rwlock");
	return(eBSFerrInternal);
	}
#endif

#ifdef _WIN32
if((m_hMtxIterReads = CreateMutex(NULL,false,NULL))==NULL)
	{
#else
if(pthread_mutex_init (&m_hMtxIterReads,NULL)!=0)
	{
	pthread_rwlock_destroy(&m_hRwLock);
#endif
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to create mutex");
	return(eBSFerrInternal);
	}

#ifdef _WIN32
if((m_hMtxMHReads = CreateMutex(NULL,false,NULL))==NULL)
	{
#else
if(pthread_mutex_init (&m_hMtxMHReads,NULL)!=0)
	{
#endif
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to create mutex");
#ifdef _WIN32
	CloseHandle(m_hMtxIterReads);
#else
	pthread_rwlock_destroy(&m_hRwLock);
	pthread_mutex_destroy(&m_hMtxIterReads);
#endif
	return(eBSFerrInternal);
	}

m_bMutexesCreated = true;
return(eBSFSuccess);
}

void
CBlitz::DeleteMutexes(void)
{
if(!m_bMutexesCreated)
	return;
#ifdef _WIN32
CloseHandle(m_hMtxIterReads);
CloseHandle(m_hMtxMHReads);
#else
pthread_mutex_destroy(&m_hMtxIterReads);
pthread_mutex_destroy(&m_hMtxMHReads);
pthread_rwlock_destroy(&m_hRwLock);
#endif
m_bMutexesCreated = false;
}

void
CBlitz::AcquireSerialise(void)
{
#ifdef _WIN32
WaitForSingleObject(m_hMtxIterReads,INFINITE);
#else
pthread_mutex_lock(&m_hMtxIterReads);
#endif
}

void
CBlitz::ReleaseSerialise(void)
{
#ifdef _WIN32
ReleaseMutex(m_hMtxIterReads);
#else
pthread_mutex_unlock(&m_hMtxIterReads);
#endif
}

void
CBlitz::AcquireSerialiseMH(void)
{
#ifdef _WIN32
WaitForSingleObject(m_hMtxMHReads,INFINITE);
#else
pthread_mutex_lock(&m_hMtxMHReads);
#endif
}

void
CBlitz::ReleaseSerialiseMH(void)
{
#ifdef _WIN32
ReleaseMutex(m_hMtxMHReads);
#else
pthread_mutex_unlock(&m_hMtxMHReads);
#endif
}

void
CBlitz::AcquireLock(bool bExclusive)
{
#ifdef _WIN32
if(bExclusive)
	AcquireSRWLockExclusive(&m_hRwLock);
else
	AcquireSRWLockShared(&m_hRwLock);
#else
if(bExclusive)
	pthread_rwlock_wrlock(&m_hRwLock);
else
	pthread_rwlock_rdlock(&m_hRwLock);
#endif
}

void
CBlitz::ReleaseLock(bool bExclusive)
{

#ifdef _WIN32
if(bExclusive)
	ReleaseSRWLockExclusive(&m_hRwLock);
else
	ReleaseSRWLockShared(&m_hRwLock);
#else
pthread_rwlock_unlock(&m_hRwLock);
#endif
}

int
CBlitz::Process(etBLZPMode PMode,		// processing mode
		int SampleNthRawRead,			// sample every Nth raw read (or read pair) for processing (1..10000)
		char *pszExprName,				// experiment name
		char *pszExprDescr,				// experiment description
		char *pszParams,				// string containing blitz parameters
		bool KMerDist,					// true if K_mer counts distributions to be reported
		etBLZSensitivity Sensitivity,	// sensitivity 0 - standard, 1 - high, 2 - very high, 3 - low sensitivity
		eALStrand AlignStrand,			// align on to watson, crick or both strands of target
		int MismatchScore,				// decrease score by this for each mismatch bp
		int ExactMatchScore,			// increase score by this for each exactly matching bp
		int GapOpenScore,				// decrease score by this for each gap open
		int  CoreLen,					// use this core length (0 if determined from total target sequence length) as the exactly matching seed length to be 5' and 3' extended whilst no more than m_MaxSubRate
		int  CoreDelta,					// offset cores by this many bp
		int MaxInsertLen,				// SAM output, accept observed insert sizes of at most this (default = 50000)
		int MaxOccKMerDepth,			// maximum depth to explore over-occurring core K-mers
		int  MinPathScore,				// only report alignment paths on any target sequence if the path score is >= this minimum score
		int QueryLenAlignedPct,			// only report alignment paths if the percentage of total aligned bases to the query sequence length is at least this percentage (1..100)
		int  MaxPathsToReport,			// report at most this many alignment paths for any query
		etBLZRsltsFomat RsltsFormat,	// output results format
		char *pszInputFile,				// name of input file containing query sequences (PE1 if PE processing when output format is SAM)
		char *pszInputFilePE2,			// name of input file containing PE2 query sequences (only applies if output format is SAM)		
		char *pszSfxFile,				// target as suffix array
		char *pszOutFile,				// where to write alignments
		int NumThreads)					// number of worker threads to use
{
int Rslt;
Init();

m_ProcMode = PMode;
m_SampleNthRawRead = SampleNthRawRead;
m_Sensitivity = Sensitivity;				
m_AlignStrand = AlignStrand;	
m_MismatchScore = MismatchScore;
m_ExactMatchScore = ExactMatchScore;
m_GapOpenScore = GapOpenScore;
m_CoreLen = CoreLen;
m_CoreDelta = CoreDelta;
m_MaxInsertLen = MaxInsertLen;
m_MaxOccKMerDepth = MaxOccKMerDepth;		
m_MinPathScore = MinPathScore;	
m_QueryLenAlignedPct = QueryLenAlignedPct;
m_MaxPathsToReport = MaxPathsToReport;
m_RsltsFormat = RsltsFormat;	
m_pszInputFile = pszInputFile;
m_pszInputFilePE2 = pszInputFilePE2;
m_pszSfxFile = pszSfxFile;		
m_pszOutFile = pszOutFile;		
m_NumThreads = NumThreads;	

if((m_pszLineBuff = new char [cAlignRprtBufferSize]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Failed to allocate memory for alignment report buffering");
	Reset(false);
	return(eBSFerrMem);
	}

if((m_pQuerySeqs = new tsQuerySeq [cMaxReadAheadQuerySeqs]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Failed to allocate memory for query sequences");
	Reset(false);
	return(eBSFerrMem);
	}
m_AllocdQuerySeqs = cMaxReadAheadQuerySeqs;
m_NumQuerySeqs = 0;
m_NxtQuerySeqIdx = 0;
m_TotSeqIDs = 0;

if(CreateMutexes()!=eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Failed to create thread synchronisation mutexes");
	Reset(false);
	return(cBSFSyncObjErr);
	}

m_mtqsort.SetMaxThreads(NumThreads);	
// open bioseq file containing suffix array for targeted assembly to align reads against
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading suffix array file '%s'", pszSfxFile);
if((m_pSfxArray = new CSfxArray()) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CSfxArray");
	Reset(false);
	return(eBSFerrObj);
	}
if((Rslt=m_pSfxArray->Open(pszSfxFile,false,false,false))!=eBSFSuccess)
	{
	while(m_pSfxArray->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pSfxArray->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open input bioseq suffix array file '%s'",pszSfxFile);
	Reset(false);
	return(Rslt);
	}

// report to user some sfx array metadata as conformation the targeted assembly is correct
strcpy(m_szTargSpecies,m_pSfxArray->GetDatasetName());
tsSfxHeaderV3 SfxHeader;
m_pSfxArray->GetSfxHeader(&SfxHeader);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Genome Assembly Name: '%s' Descr: '%s' Title: '%s' Version: %d",
					 m_szTargSpecies,SfxHeader.szDescription,SfxHeader.szTitle,SfxHeader.Version);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Assembly has blocks: %d, max block size: %lu",SfxHeader.NumSfxBlocks,SfxHeader.SfxBlockSize);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading genome assembly suffix array...");
if((Rslt=m_pSfxArray->SetTargBlock(1))<0)
	{
	while(m_pSfxArray->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pSfxArray->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to load genome assembly suffix array");
	return(Rslt);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Genome assembly suffix array loaded");
// from the total sequence length then determine the core length to use
// the autodetermined core length is that length such that on average there would be expected less than one copy of the core sequence in a random sequence of same length as targeted genome
UINT64 TotSeqsLen = m_pSfxArray->GetTotSeqsLen();

if(CoreLen == 0)
	{
	int AutoCoreLen = 1;
	while (TotSeqsLen >>= 2)
		AutoCoreLen++;
	CoreLen = AutoCoreLen;

	switch(Sensitivity) {
		case eBLZSdefault:			// default sensitivity
			break;
		case eBLZSMoreSens:			// more sensitive - slower
			CoreLen -= 1;
			break;
		case eBLZSUltraSens:		// ultra sensitive - much slower
			CoreLen -= 2;
			break;
		case eBLZSLessSens:			// less sensitive - quicker
		default:
			CoreLen += 2;
		}
	if(CoreLen > cMaxCoreLen)
		CoreLen = cMaxCoreLen;
	m_CoreLen = CoreLen;
	}
gDiagnostics.DiagOut(eDLInfo, gszProcName, "Using core length : %d", m_CoreLen);

if(CoreDelta == 0)
	{
	switch(Sensitivity) {
		case eBLZSdefault:			// default sensitivity
			CoreDelta = max((CoreLen+1) / 2, 4);
			break;
		case eBLZSMoreSens:			// more sensitive - slower
			CoreDelta = max((CoreLen+2) / 3, 3);
			break;
		case eBLZSUltraSens:		// ultra sensitive - much slower
			CoreDelta = max((CoreLen+3) / 4, 2);
			break;
		case eBLZSLessSens:			// less sensitive - quicker
		default:
			CoreDelta = max((CoreLen + 1) / 2, 8);;
		}
	m_CoreDelta = CoreDelta;
	}
gDiagnostics.DiagOut(eDLInfo, gszProcName, "Using core delta : %d", m_CoreDelta);

m_MinPathScore = MinPathScore;
if(MinPathScore == 0)
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Using automatic query sequence length dependent minimum path score");
else
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Using minimum path score : %d", m_MinPathScore);
m_MaxIter = MaxOccKMerDepth;

// restrict the max core iterations and substitutions thresholding according to the requested sensitivity
switch(Sensitivity) {
	case eBLZSdefault:			// default sensitivity
		if(!MaxOccKMerDepth)
			m_MaxIter = cDfltSensCoreIters;
		break;
	case eBLZSMoreSens:			// more sensitive - slower
		if(!MaxOccKMerDepth)
			m_MaxIter = cMoreSensCoreIters;
		break;
	case eBLZSUltraSens:			// ultra sensitive - much slower
		if(!MaxOccKMerDepth)
			m_MaxIter = cUltraSensCoreIters;
		break;
	case eBLZSLessSens:			// less sensitive - quicker
	default:
		if(!MaxOccKMerDepth)
			m_MaxIter = cMinSensCoreIters;
		break;
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Using maximum depth to explore over-occurring seed K-mers : %d",m_MaxIter);

m_pSfxArray->SetMaxIter(m_MaxIter);

if(CoreLen <= cMaxKmerLen)
	{
	if((Rslt = m_pSfxArray->InitOverOccKMers(CoreLen,m_MaxIter))!=eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Failed to initialise for over occurring K-mers");
		Reset(false);
		return(Rslt);
		}
	}

if(KMerDist == true)
	{
	if ((m_pKmerOccsDist = new uint32_t[cMaxOccKMerDepth + 1]) == NULL)				// allocated to hold Kmer count distributions (up to cMaxOccKMerDepth counts)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Failed to allocate memory for KMer distribution counts");
		Reset(false);
		return(eBSFerrMem);
		}
	memset(m_pKmerOccsDist, 0, sizeof(int) * (cMaxOccKMerDepth + 1));
	m_pSfxArray->InitKMerCountDists(cMaxOccKMerDepth + 1, m_pKmerOccsDist);
	}

if(RsltsFormat == eBLZRsltsSAM)
	m_MaxQuerySeqLen = cSAMtruncSeqLen;		// truncating SAM query sequences at this length
else
	m_MaxQuerySeqLen = cMaxQuerySeqLen;

if(RsltsFormat == eBLZRsltsSQLite)
	{
	if((m_pSQLitePSL = new CSQLitePSL)==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate instance of CSQLitePSL");
		return(eBSFerrObj);
		}

	if(m_pSQLitePSL->CreateDatabase(pszOutFile,false)==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create SQLite database '%s'",pszOutFile);
		delete m_pSQLitePSL;
		m_pSQLitePSL = NULL;
		return(eBSFerrObj);
		}

	if((m_ExprID = m_pSQLitePSL->CreateExperiment(pszExprName,pszOutFile,pszInputFile,pszSfxFile,pszExprDescr,pszParams,0)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to initialise SQLite database '%s' with experiment details",pszOutFile);
		delete m_pSQLitePSL;
		m_pSQLitePSL = NULL;
		return(m_ExprID);
		}

	if((Rslt=m_pSQLitePSL->BeginPopulatingTables())!=eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to begin SQLite database '%s' table populating",pszOutFile);
		delete m_pSQLitePSL;
		m_pSQLitePSL = NULL;
		return(Rslt);
		}

	// add summary instances for all the target sequences
	char szSeqIdent[100];
	uint32_t SeqLen;
	int NumEntryIDs;
	int CurEntryID;
	NumEntryIDs = m_pSfxArray->GetNumEntries();
	for(CurEntryID = 1; CurEntryID <= NumEntryIDs; CurEntryID+=1)
		{
		m_pSfxArray->GetIdentName(CurEntryID,sizeof(szSeqIdent)-1,szSeqIdent);
		SeqLen = m_pSfxArray->GetSeqLen(CurEntryID);
		m_pSQLitePSL->AddAlignSummary(m_ExprID,NULL,0,szSeqIdent,SeqLen);		
		}
	}


// reads are loaded asynchronously to the alignment processing
if((Rslt=InitLoadQuerySeqs()) < eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Failed to load reads");
	Reset(false);
	return(Rslt);
	}

char szNonAlignedFilePE1[_MAX_PATH];
char szNonAlignedFilePE2[_MAX_PATH];

sprintf(szNonAlignedFilePE1,"%s.PE1.Unaligned.fasta", pszOutFile);
#ifdef _WIN32
m_hOutNonAlignedFilePE1 = open(szNonAlignedFilePE1, (O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC), (_S_IREAD | _S_IWRITE));
#else
if ((m_hOutNonAlignedFilePE1 = open(szNonAlignedFilePE1, O_WRONLY | O_CREAT, S_IREAD | S_IWRITE)) != -1)
if (ftruncate(m_hOutNonAlignedFilePE1, 0) != 0)
{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to truncate %s - %s", szNonAlignedFilePE1, strerror(errno));
	Reset(false);
	return(eBSFerrCreateFile);
}
#endif
if (m_hOutNonAlignedFilePE1 < 0)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: unable to create/truncate output file '%s'", szNonAlignedFilePE1);
	Reset(false);
	return(eBSFerrCreateFile);
	}

if(pszInputFilePE2 != NULL && pszInputFilePE2[0] != '\0')
	{
	sprintf(szNonAlignedFilePE2, "%s.PE2.Unaligned.fasta", pszOutFile);
#ifdef _WIN32
	m_hOutNonAlignedFilePE2 = open(szNonAlignedFilePE2, (O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC), (_S_IREAD | _S_IWRITE));
#else
	if ((m_hOutNonAlignedFilePE2 = open(szNonAlignedFilePE2, O_WRONLY | O_CREAT, S_IREAD | S_IWRITE)) != -1)
		if (ftruncate(m_hOutNonAlignedFilePE2, 0) != 0)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to truncate %s - %s", szNonAlignedFilePE2, strerror(errno));
		Reset(false);
		return(eBSFerrCreateFile);
		}
#endif
	if (m_hOutNonAlignedFilePE2 < 0)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: unable to create/truncate output file '%s'", szNonAlignedFilePE2);
		Reset(false);
		return(eBSFerrCreateFile);
		}
	}
else
	{
	szNonAlignedFilePE2[0] = '\0';
	m_hOutNonAlignedFilePE2 = -1;
	}

if(RsltsFormat != eBLZRsltsSQLite)
	{
#ifdef _WIN32
	m_hOutFile = open(pszOutFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE));
#else
if((m_hOutFile = open(pszOutFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE))!=-1)
	if(ftruncate(m_hOutFile,0)!=0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",pszOutFile,strerror(errno));
		Reset(false);
		return(eBSFerrCreateFile);
		}
#endif
	if(m_hOutFile < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output file '%s'",pszOutFile);
		Reset(false);
		return(eBSFerrCreateFile);
		}

	// write out format specific headers
	switch(RsltsFormat) {
		case eBLZRsltsPSL:
			m_szLineBuffIdx = sprintf(m_pszLineBuff,"psLayout version 3\nGenerated by %s %s, Version %s\n",gszProcName,gpszSubProcess->pszName,kit4bversion);
			m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"match	mis- 	rep. 	N's	Q gap	Q gap	T gap	T gap	strand	Q        	Q   	Q    	Q  	T        	T   	T    	T  	block	blockSizes 	qStarts	 tStarts\n");
			m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"     	match	match	   	count	bases	count	bases	      	name     	size	start	end	name     	size	start	end	count\n");
			m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"---------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
			break;
		case eBLZRsltsPSLX:
			m_szLineBuffIdx = sprintf(m_pszLineBuff,"psLayout version 3\nGenerated by %s %s, Version %s\n",gszProcName,gpszSubProcess->pszName,kit4bversion);
			m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"match	mis- 	rep. 	N's	Q gap	Q gap	T gap	T gap	strand	Q        	Q   	Q    	Q  	T        	T   	T    	T  	block	blockSizes 	qStarts	 tStarts\n");
			m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"     	match	match	   	count	bases	count	bases	      	name     	size	start	end	name     	size	start	end	count\n");
			m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"---------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
			break;
		case eBLZRsltsMAF:
			m_szLineBuffIdx = sprintf(m_pszLineBuff,"##maf version=1 scoring=blitz");
			break;
		case eBLZRsltsBED:
			m_szLineBuffIdx = sprintf(m_pszLineBuff,"track type=bed name=\"Blitz\" description=\"kit4b Blitz alignments\"\n");
			break;
		case eBLZRsltsSAM:			// output as SAM
			m_szLineBuffIdx = sprintf(m_pszLineBuff, "@HD\tVN:1.4\tSO:unsorted\n");
			// add the target sequences to SAM header
			char szSeqIdent[100];
			uint32_t SeqLen;
			int NumEntryIDs;
			int CurEntryID;
			NumEntryIDs = m_pSfxArray->GetNumEntries();
			for (CurEntryID = 1; CurEntryID <= NumEntryIDs; CurEntryID += 1)
				{
				m_pSfxArray->GetIdentName(CurEntryID, sizeof(szSeqIdent) - 1, szSeqIdent);
				SeqLen = m_pSfxArray->GetSeqLen(CurEntryID);
				m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx], "@SQ\tAS:%s\tSN:%s\tLN:%u\n", m_szTargSpecies,szSeqIdent, SeqLen);
				if (m_szLineBuffIdx > (cAlignRprtBufferSize * 9) / 10)
					{
					CUtility::SafeWrite(m_hOutFile, m_pszLineBuff, m_szLineBuffIdx);
					m_szLineBuffIdx = 0;
					}
				}
			break;
		}
	m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"@PG\tID:kit4b_blitz\tVN:%s\n", kit4bversion);
	if(m_szLineBuffIdx)
		CUtility::SafeWrite(m_hOutFile,m_pszLineBuff,m_szLineBuffIdx);
	m_szLineBuffIdx = 0;
	}

m_szLineBuffIdx = 0;

InitQuerySeqThreads(NumThreads,cNumAllocdAlignNodes);	

// pickup the query sequence loader thread, if the alignment processing threads all finished then the loader thread should also have finished
#ifdef _WIN32
if(m_hThreadLoadQuerySeqs != NULL)
	{
	while(WAIT_TIMEOUT == WaitForSingleObject(m_hThreadLoadQuerySeqs, 5000))
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: waiting on query sequence loader thread to complete");
		}
	CloseHandle(m_hThreadLoadQuerySeqs);
	m_hThreadLoadQuerySeqs = NULL;
	}
#else
if(m_ThreadLoadQuerySeqsID != 0)
	{
	struct timespec ts;
	int JoinRlt;
	clock_gettime(CLOCK_REALTIME, &ts);
	ts.tv_sec += 5;
	while((JoinRlt = pthread_timedjoin_np(m_ThreadLoadQuerySeqsID, NULL, &ts)) != 0)
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: waiting on query sequence loader thread to complete");
		ts.tv_sec += 60;
		}
	}
#endif

// Checking here that the reads were all loaded w/o any major dramas!
if(m_LoadQuerySeqsRslt < 0)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: problem loading");
	Reset(false);
	return(m_LoadQuerySeqsRslt);
	}

if(m_hOutFile != -1)
	{
	if(m_szLineBuffIdx > 0)
		CUtility::SafeWrite(m_hOutFile,m_pszLineBuff,m_szLineBuffIdx);
	m_szLineBuffIdx = 0;
#ifdef _WIN32
	_commit(m_hOutFile);
#else
	fsync(m_hOutFile);
#endif
	close(m_hOutFile);
	m_hOutFile = -1;
	}

// now for KMer distributions
if(m_pKmerOccsDist != NULL)
	{
	char szKMerDistFile[_MAX_PATH];
	sprintf(szKMerDistFile,"%s.KMer%dDist.csv", pszOutFile,m_CoreLen);
#ifdef _WIN32
	m_hOutFile = open(szKMerDistFile, (O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC), (_S_IREAD | _S_IWRITE));
#else
	if ((m_hOutFile = open(szKMerDistFile, O_WRONLY | O_CREAT, S_IREAD | S_IWRITE)) != -1)
		if (ftruncate(m_hOutFile, 0) != 0)
		{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to truncate %s - %s", szKMerDistFile, strerror(errno));
			Reset(false);
			return(eBSFerrCreateFile);
		}
#endif
	if (m_hOutFile < 0)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: unable to create/truncate output file '%s'", szKMerDistFile);
		Reset(false);
		return(eBSFerrCreateFile);
		}

	ReportKMerDist(m_CoreLen,cMaxOccKMerDepth, m_pKmerOccsDist);

#ifdef _WIN32
	_commit(m_hOutFile);
#else
	fsync(m_hOutFile);
#endif
	close(m_hOutFile);
	m_hOutFile = -1;
	}

if(RsltsFormat == eBLZRsltsSQLite)
	{
	if((Rslt=m_pSQLitePSL->EndPopulatingTables())!=eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to complete SQLite database '%s' table populating",pszOutFile);
		delete m_pSQLitePSL;
		m_pSQLitePSL = NULL;
		return(Rslt);
		}
	delete m_pSQLitePSL;
	m_pSQLitePSL = NULL;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"SQLite database ready for use");
	}

Reset(false);
return(Rslt);
}




#ifdef _WIN32
unsigned __stdcall AlignQuerySeqsThread(void * pThreadPars)
#else
void *AlignQuerySeqsThread(void * pThreadPars)
#endif
{
int Rslt;
tsThreadQuerySeqsPars *pPars = (tsThreadQuerySeqsPars *)pThreadPars;			// makes it easier not having to deal with casts!
CBlitz *pBlitzer = (CBlitz *)pPars->pThis;
if(pPars->bIsSAMOutput)
	{
	if(pPars->bIsSAMPE)
		Rslt = pBlitzer->ProcAlignSAMQuerySeqsPE(pPars);
	else
		Rslt = pBlitzer->ProcAlignSAMQuerySeqsSE(pPars);
	}
else
	Rslt = pBlitzer->ProcAlignQuerySeqs(pPars);
pPars->Rslt = Rslt;
#ifdef _WIN32
_endthreadex(0);
return(eBSFSuccess);
#else
pthread_exit(NULL);
#endif
}

int
CBlitz::InitQuerySeqThreads(int NumThreads,			// use this many threads
							int AlignNodes)			// each thread is allocatd this many subsequence alignment nodes
{
bool bIsSAMPE;
tsThreadQuerySeqsPars *pThreads;
tsThreadQuerySeqsPars *pThread;
m_QueriesPaths = 0;
m_ReportedPaths = 0;
bIsSAMPE = (m_pszInputFilePE2 == NULL || m_pszInputFilePE2[0] == 0) ? false : true;

if ((pThreads = new tsThreadQuerySeqsPars[NumThreads]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessReadsetDist: Memory allocation for thread context failed");
	return(eBSFerrMem);
	}
memset(pThreads,0,sizeof(tsThreadQuerySeqsPars) * NumThreads);
int ThreadIdx;
pThread = pThreads;
for (ThreadIdx = 1; ThreadIdx <= NumThreads; ThreadIdx++, pThread++)
	{
	pThread->ThreadIdx = ThreadIdx;
	pThread->pThis = this;
	pThread->pAllocdAlignNodes = new tsQueryAlignNodes [cNumAllocdAlignNodes];
	pThread->ppFirst2Rpts = new tsQueryAlignNodes * [cNumAllocdAlignNodes];
	pThread->NumAllocdAlignNodes = cNumAllocdAlignNodes;

	pThread->bIsSAMOutput = m_RsltsFormat == eBLZRsltsSAM ? true : false;
	pThread->bIsSAMPE = bIsSAMPE;
	if(pThread->bIsSAMPE)
		{
		pThread->pAllocdAlignNodesPE2 = new tsQueryAlignNodes[cNumAllocdAlignNodes];
		pThread->ppFirst2RptsPE2 = new tsQueryAlignNodes *[cNumAllocdAlignNodes];
		pThread->NumAllocdAlignNodesPE2 = cNumAllocdAlignNodes;
		}
	else
		{
		pThread->pAllocdAlignNodesPE2 = NULL;
		pThread->ppFirst2RptsPE2 = NULL;
		pThread->NumAllocdAlignNodesPE2 = 0;
		}

#ifdef _WIN32
	pThread->threadHandle = (HANDLE)_beginthreadex(NULL, 0x0fffff, AlignQuerySeqsThread, pThread, 0, &pThread->threadID);
#else
	pThread->threadRslt = pthread_create(&pThread->threadID, NULL, AlignQuerySeqsThread, pThread);
#endif
	}

// allow threads a few seconds to startup
#ifdef _WIN32
Sleep(2000);
#else
sleep(2);
#endif
uint32_t ReportedPaths;
uint32_t PrevReportedPaths = 0;
uint32_t QueriesPaths;
uint32_t PrevQueriesPaths = 0;
uint32_t NumQueriesProc;
uint32_t PrevNumQueriesProc = 0;


gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: Generated 0 alignment paths for 0 query %s sequences from 0 processed", bIsSAMPE ? "paired" : "single");
pThread = pThreads;
for (ThreadIdx = 0; ThreadIdx < NumThreads; ThreadIdx++, pThread++)
	{
#ifdef _WIN32
	while (WAIT_TIMEOUT == WaitForSingleObject(pThread->threadHandle, 60000))
		{
		AcquireSerialise();
		ReportedPaths = m_ReportedPaths;
		QueriesPaths = m_QueriesPaths;
		NumQueriesProc = m_NumQueriesProc;
		if(ReportedPaths > PrevReportedPaths || QueriesPaths > PrevQueriesPaths || NumQueriesProc > PrevNumQueriesProc)
			{
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: Generated %u alignment paths for %u query %s sequences from %d processed",ReportedPaths,QueriesPaths, bIsSAMPE ? "paired" : "single",NumQueriesProc);
			PrevReportedPaths = ReportedPaths;
			PrevQueriesPaths = QueriesPaths;
			PrevNumQueriesProc = NumQueriesProc;
			}
		ReleaseSerialise();
		};
	CloseHandle(pThread->threadHandle);
#else
	struct timespec ts;
	int JoinRlt;
	clock_gettime(CLOCK_REALTIME, &ts);
	ts.tv_sec += 60;
	while ((JoinRlt = pthread_timedjoin_np(pThread->threadID, NULL, &ts)) != 0)
		{
		AcquireSerialise();
		ReportedPaths = m_ReportedPaths;
		QueriesPaths = m_QueriesPaths;
		NumQueriesProc = m_NumQueriesProc;
		if(ReportedPaths > PrevReportedPaths || QueriesPaths > PrevQueriesPaths || NumQueriesProc > PrevNumQueriesProc)
			{
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: Generated %u alignment paths for %u query %s sequences from %d processed",ReportedPaths,QueriesPaths, bIsSAMPE ? "paired" : "single", NumQueriesProc);
			PrevReportedPaths = ReportedPaths;
			PrevQueriesPaths = QueriesPaths;
			PrevNumQueriesProc = NumQueriesProc;
			}
		ReleaseSerialise();
		ts.tv_sec += 60;
		}
#endif
	if(pThread->pAllocdAlignNodes != NULL)
		{
		delete pThread->pAllocdAlignNodes;
		pThread->pAllocdAlignNodes = NULL; 
		}
	if(pThread->ppFirst2Rpts != NULL)
		{
		delete pThread->ppFirst2Rpts;
		pThread->ppFirst2Rpts = NULL; 
		}
	pThread->NumAllocdAlignNodes = 0;
	if (pThread->pAllocdAlignNodesPE2 != NULL)
		{
		delete pThread->pAllocdAlignNodesPE2;
		pThread->pAllocdAlignNodesPE2 = NULL;
		}
	if (pThread->ppFirst2RptsPE2 != NULL)
		{
		delete pThread->ppFirst2RptsPE2;
		pThread->ppFirst2RptsPE2 = NULL;
		}
	pThread->NumAllocdAlignNodesPE2 = 0;
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Completed reporting %u alignment paths for %u query %s sequences from %d processed",m_ReportedPaths,m_QueriesPaths, bIsSAMPE ? "paired" : "single", m_NumQueriesProc);

if(pThreads != NULL)
	delete pThreads;
return(0);
}

int
CBlitz::ReportNonAligned(char *pszDescPE1, int LenSeqPE1, uint8_t *pSeqPE1, char *pszDescPE2, int LenSeqPE2, uint8_t *pSeqPE2)
{
char Base;
int LineLen;
int IdxPE1;
char szFastaBuffPE1[((cMaxQuerySeqIdentLen + 3) + cSAMtruncSeqLen + ((cSAMtruncSeqLen + 80) / 79)) * 2];
int FastaBuffIdxPE1;

int IdxPE2;
char szFastaBuffPE2[((cMaxQuerySeqIdentLen + 3) + cSAMtruncSeqLen + ((cSAMtruncSeqLen + 80) / 79)) * 2];
int FastaBuffIdxPE2;

szFastaBuffPE1[0] = '\0';
FastaBuffIdxPE1 = 0;
szFastaBuffPE2[0] = '\0';
FastaBuffIdxPE2 = 0;

if (pszDescPE1 != NULL && pszDescPE1[0] != '\0' && LenSeqPE1 > 0)
	{
	FastaBuffIdxPE1 = sprintf(szFastaBuffPE1,">%s\n", pszDescPE1);
	LineLen = 0;
	for(IdxPE1 = 0; IdxPE1 < LenSeqPE1; IdxPE1++)
		{
		switch (*pSeqPE1++) {
			case eBaseA:
				Base = 'A';
				break;
			case eBaseC:
				Base = 'C';
				break;
			case eBaseG:
				Base = 'G';
				break;
			case eBaseT:
				Base = 'T';
				break;
			default:
				Base = 'N';
			}
		szFastaBuffPE1[FastaBuffIdxPE1++] = Base;
		if(++LineLen >= 79)
			{
			szFastaBuffPE1[FastaBuffIdxPE1++] = '\n';
			LineLen = 0;
			}
		}
	if(szFastaBuffPE1[FastaBuffIdxPE1-1] != '\n')
		szFastaBuffPE1[FastaBuffIdxPE1++] = '\n';
	}

if(pszDescPE2 != NULL && pszDescPE2[0] != '\0' && LenSeqPE2 > 0)
	{
	FastaBuffIdxPE2 = sprintf(szFastaBuffPE2, ">%s\n", pszDescPE2);
	LineLen = 0;
	for (IdxPE2 = 0; IdxPE2 < LenSeqPE2; IdxPE2++)
	{
		switch (*pSeqPE2++) {
		case eBaseA:
			Base = 'A';
			break;
		case eBaseC:
			Base = 'C';
			break;
		case eBaseG:
			Base = 'G';
			break;
		case eBaseT:
			Base = 'T';
			break;
		default:
			Base = 'N';
		}
		szFastaBuffPE2[FastaBuffIdxPE2++] = Base;
		if (++LineLen >= 79)
		{
			szFastaBuffPE2[FastaBuffIdxPE2++] = '\n';
			LineLen = 0;
		}
	}
	if (szFastaBuffPE2[FastaBuffIdxPE2 - 1] != '\n')
		szFastaBuffPE2[FastaBuffIdxPE2++] = '\n';
	}

if(m_hOutNonAlignedFilePE1 != -1 && szFastaBuffPE1[0] != '\0' && FastaBuffIdxPE1 > 0)
	{
	AcquireSerialise();
	CUtility::SafeWrite(m_hOutNonAlignedFilePE1, szFastaBuffPE1, FastaBuffIdxPE1);

	if (m_hOutNonAlignedFilePE2 != -1 && szFastaBuffPE2[0] != '\0' && FastaBuffIdxPE2 > 0)
		CUtility::SafeWrite(m_hOutNonAlignedFilePE2, szFastaBuffPE2, FastaBuffIdxPE2);
	ReleaseSerialise();
	}
return(0);
}

int
CBlitz::ProcAlignSAMQuerySeqsSE(tsThreadQuerySeqsPars *pPars) // single ended processing only
{
int Rslt;
int MaxIter;
int NumQueryPathsRprtd;
int NumMatches;
int QuerySeqLen;
int SeqID;
uint8_t *pQuerySeq;
char szQuerySeqIdent[cMaxQuerySeqIdentLen + 1];
uint32_t NumHeadNodes;
uint32_t SAMFlags;

tsQueryAlignNodes *pHeadNode;
bool bStrand;
int PathHiScore;
int MinPathHiScore;
uint32_t SortedPathIdx;
uint32_t NumPathNodes;
uint32_t QueryPathStartOfs;
uint32_t QueryPathEndOfs;
uint32_t TargPathStartOfs;
uint32_t TargPathEndOfs;
uint32_t TargSeqLen;
char szTargName[100];
int MinPathScore;

while((Rslt = DequeueQuerySeq(cMaxQuerySeqIdentLen + 1, &SeqID, szQuerySeqIdent, &QuerySeqLen, &pQuerySeq)) == 1)
	{
	AcquireSerialise();
	m_NumQueriesProc += 1;
	ReleaseSerialise();

	NumQueryPathsRprtd = 0;
	MaxIter = m_MaxIter;

	NumMatches = m_pSfxArray->LocateQuerySeqs(SeqID, pQuerySeq, QuerySeqLen, m_CoreLen, m_CoreDelta, m_AlignStrand, pPars->NumAllocdAlignNodes, pPars->pAllocdAlignNodes, MaxIter, m_ExactMatchScore,m_MismatchScore);
	if(NumMatches >= 1)
		{
		if (NumMatches > 1)	// sorting by TargSeqID.QueryID.FlgStrand.TargStartOfs.QueryStartOfs
			qsort(pPars->pAllocdAlignNodes, NumMatches, sizeof(tsQueryAlignNodes), SortQueryAlignNodes);
		if(m_MinPathScore > 0)
			MinPathScore = m_MinPathScore;
		else
			MinPathScore = ((QuerySeqLen - 5) * m_ExactMatchScore) / 3;
		MinPathHiScore = MinPathScore;

		NumHeadNodes = IdentifyHighScorePaths(MinPathScore, m_MaxPathsToReport, QuerySeqLen, NumMatches, pPars->pAllocdAlignNodes, pPars->ppFirst2Rpts);
		if (NumHeadNodes)
			{
			for (SortedPathIdx = 0; SortedPathIdx < min(NumHeadNodes, (uint32_t)m_MaxPathsToReport); SortedPathIdx++)
				{
				pHeadNode = pPars->ppFirst2Rpts[SortedPathIdx];
				bStrand = pHeadNode->FlgStrand ? true : false;
				ConsolidateNodes(bStrand,QuerySeqLen, pQuerySeq, SortedPathIdx, pPars->pAllocdAlignNodes, pPars->ppFirst2Rpts);
				if ((NumPathNodes = CharacterisePath(QuerySeqLen,pQuerySeq,SortedPathIdx, pPars->pAllocdAlignNodes, pPars->ppFirst2Rpts, &QueryPathEndOfs, &TargPathEndOfs)) == 0)
					continue;
			
				PathHiScore = pHeadNode->HiScore;
				PathHiScore = (int)(((double)PathHiScore / ((int64_t)QuerySeqLen * m_ExactMatchScore)) * 255.0);
				if(PathHiScore >= MinPathHiScore)
					{
					bStrand = pHeadNode->FlgStrand ? true : false;
					TargPathStartOfs = pHeadNode->TargSeqLoci;
					QueryPathStartOfs = pHeadNode->QueryStartOfs;
			
					TargSeqLen = m_pSfxArray->GetSeqLen(pHeadNode->TargSeqID);
					m_pSfxArray->GetIdentName(pHeadNode->TargSeqID, sizeof(szTargName), szTargName);
					if (bStrand == true)
						SAMFlags = cSAMFlgAS;
					else
						SAMFlags = 0;
					if(NumQueryPathsRprtd > 0)
						SAMFlags |= cSAMFlgNotPrimary;
					ReportAsSAM(SAMFlags,				// use as reported SAM flags
						bStrand == true ? '-' : '+',	// query sequence strand, '+' or '-'
						PathHiScore,					// score for this path
						szQuerySeqIdent,				// identifies this query sequence
						QuerySeqLen,					// Query sequence length
						pQuerySeq,						// the query sequence as etSeqBase's
						QueryPathStartOfs,				// Alignment start offset in query, 0..QSize-1
						QueryPathEndOfs,				// Alignment end position in query
						szTargName,						// aligning to this target
						TargSeqLen,						// Target sequence size
						TargPathStartOfs+1,				// target alignment starts at this starting loci, 1..TargSeqLen
						TargPathEndOfs+1,				// ending at this loci (inclusive)
						'*',							// Reference sequence name of the primary alignment of the NEXT read in the template, '*' if unknown, '=' if PE and both ends align to same reference
						0,								// 1-based Position of the primary alignment of the NEXT read in the template. Set as 0 when the information is unavailable
						0,							    // signed template length, If all segments are mapped to the same reference, the signed observed template length equals the number of bases from the leftmost mapped base to the rightmost mapped base
						NumPathNodes,					// number of alignment nodes in alignment path
						SortedPathIdx,					// alignment path starts from this node - 0 based
						pPars->pAllocdAlignNodes,		// alignment nodes
						pPars->ppFirst2Rpts);			// allocated to hold ptrs to alignment nodes which are marked as being FlgFirst2tRpt
					NumQueryPathsRprtd += 1;
					}
				}
			}
		}
	if(NumQueryPathsRprtd)
		{
		AcquireSerialise();
		m_ReportedPaths += NumQueryPathsRprtd;
		m_QueriesPaths++;
		ReleaseSerialise();
		}
	else
		ReportNonAligned(szQuerySeqIdent, QuerySeqLen, pQuerySeq);
	delete pQuerySeq;
	}

gDiagnostics.DiagOut(eDLInfo, gszProcName, "Progress: Thread %d completed", pPars->ThreadIdx);
return(0);
}

int
CBlitz::ProcAlignSAMQuerySeqsPE(tsThreadQuerySeqsPars *pPars)
{
int Rslt;
int MaxIter;
int NumQueryPathsRprtdPE1;
int NumMatchesPE1;
int QuerySeqLenPE1;
int SeqIDPE1;
uint8_t *pQuerySeqPE1;
char szQuerySeqIdent[cMaxQuerySeqIdentLen + 1];
int NumQueryPathsRprtdPE2;
int NumMatchesPE2;
int QuerySeqLenPE2;
int SeqIDPE2;
uint8_t *pQuerySeqPE2;
char szQuerySeqIdentPE2[cMaxQuerySeqIdentLen + 1];

uint32_t NumHeadNodesPE1;
uint32_t NumHeadNodesPE2;

bool bQueryAtLeast1Path;
uint32_t SAMFlags;

tsQueryAlignNodes *pHeadNodePE1;
bool bStrandPE1;
int MinPathScorePE1;
int MinPathHiScorePE1;
int PathHiScorePE1;
uint32_t SortedPathIdxPE1;
uint32_t NumPathNodesPE1;
uint32_t QueryPathStartOfsPE1;
uint32_t QueryPathEndOfsPE1;
uint32_t TargPathStartOfsPE1;
uint32_t TargPathEndOfsPE1;
uint32_t TargSeqLenPE1;
tsQueryAlignNodes *pHeadNodePE2;
bool bStrandPE2;
int MinPathScorePE2;
int PathHiScorePE2;
int MinPathHiScorePE2;
uint32_t SortedPathIdxPE2;
uint32_t NumPathNodesPE2;
uint32_t QueryPathStartOfsPE2;
uint32_t QueryPathEndOfsPE2;
uint32_t TargPathStartOfsPE2;
uint32_t TargPathEndOfsPE2;
uint32_t TargSeqLenPE2;

INT64 DeltaPathStarts;
uint32_t CombinedScore;

uint32_t NoHeadNodes;
char szTargName[100];
char szTargNamePE2[100];

NumQueryPathsRprtdPE1 = 0;
NumQueryPathsRprtdPE2 = 0;
NoHeadNodes = 0;
bQueryAtLeast1Path = false;

MaxIter = m_MaxIter;
int CoreDelta = m_CoreDelta;

// getting 2 reads at a time, PE1 and PE2
while((Rslt = DequeueQuerySeq(cMaxQuerySeqIdentLen + 1, &SeqIDPE1, szQuerySeqIdent, &QuerySeqLenPE1, &pQuerySeqPE1, &SeqIDPE2, szQuerySeqIdentPE2, &QuerySeqLenPE2, &pQuerySeqPE2))==2)
	{
	AcquireSerialise();
	m_NumQueriesProc += 1;
	ReleaseSerialise();
	bQueryAtLeast1Path = false;

		// starting with none overlapping cores, expecting around 50% of reads to have at least a couple of non-overlapping cores
	int NumNoPE1Matches = 0;
	int NumNoPE2Matches = 0;

	do {
		NumQueryPathsRprtdPE2 = 0;
		NumQueryPathsRprtdPE1 = 0;
		NumHeadNodesPE1 = 0;
		NumHeadNodesPE2 = 0;
		NumMatchesPE1 = m_pSfxArray->LocateQuerySeqs(SeqIDPE1, pQuerySeqPE1, QuerySeqLenPE1, m_CoreLen, CoreDelta, m_AlignStrand, pPars->NumAllocdAlignNodes, pPars->pAllocdAlignNodes, MaxIter, m_ExactMatchScore, m_MismatchScore);
		if (NumMatchesPE1)
			{
			if (NumMatchesPE1 > 1)	// sorting by TargSeqID.QueryID.FlgStrand.TargStartOfs.QueryStartOfs
				qsort(pPars->pAllocdAlignNodes, NumMatchesPE1, sizeof(tsQueryAlignNodes), SortQueryAlignNodes);
			if (m_MinPathScore > 0)
				MinPathScorePE1 = m_MinPathScore;
			else
				MinPathScorePE1 = ((QuerySeqLenPE1 - 5) * m_ExactMatchScore) / 3;
			MinPathHiScorePE1 = MinPathScorePE1;
			NumHeadNodesPE1 = IdentifyHighScorePaths(MinPathScorePE1, m_MaxPathsToReport+2, QuerySeqLenPE1, NumMatchesPE1, pPars->pAllocdAlignNodes, pPars->ppFirst2Rpts);
			}
		else
			NumNoPE1Matches++;

		if(NumHeadNodesPE1 > 0)
			{
			NumMatchesPE2 = m_pSfxArray->LocateQuerySeqs(SeqIDPE2, pQuerySeqPE2, QuerySeqLenPE2, m_CoreLen, CoreDelta, m_AlignStrand, pPars->NumAllocdAlignNodesPE2, pPars->pAllocdAlignNodesPE2, MaxIter, m_ExactMatchScore, m_MismatchScore);
			if (NumMatchesPE2)
				{
				if (NumMatchesPE2 > 1)	// sorting by TargSeqID.QueryID.FlgStrand.TargStartOfs.QueryStartOfs
					qsort(pPars->pAllocdAlignNodesPE2, NumMatchesPE2, sizeof(tsQueryAlignNodes), SortQueryAlignNodes);
				if (m_MinPathScore > 0)
					MinPathScorePE2 = m_MinPathScore;
				else
					MinPathScorePE2 = ((QuerySeqLenPE2 - 5) * m_ExactMatchScore) / 3;
				MinPathHiScorePE2 = MinPathScorePE2;
				NumHeadNodesPE2 = IdentifyHighScorePaths(MinPathScorePE2, m_MaxPathsToReport+2, QuerySeqLenPE2, NumMatchesPE2, pPars->pAllocdAlignNodesPE2, pPars->ppFirst2RptsPE2);
				}
			else
				NumNoPE2Matches++;
			}

		// following is an attempt to actually pair the PEs!!!!
		if(NumHeadNodesPE1 && NumHeadNodesPE2)		// need to have both PE1 and PE2 alignments!
			{
			uint32_t BestCombinedScore = 0;
			uint32_t BestDeltaPathStart = m_MaxInsertLen;
			uint32_t BestSortedPathIdxPE1 = 0;
			uint32_t BestSortedPathIdxPE2 = 0;
			for (SortedPathIdxPE1 = 0; NumQueryPathsRprtdPE1 < m_MaxPathsToReport && SortedPathIdxPE1 < min(NumHeadNodesPE1, (uint32_t)m_MaxPathsToReport + 1); SortedPathIdxPE1++)
				{
				pHeadNodePE1 = pPars->ppFirst2Rpts[SortedPathIdxPE1];
				ConsolidateNodes(pHeadNodePE1->FlgStrand,QuerySeqLenPE1, pQuerySeqPE1, SortedPathIdxPE1, pPars->pAllocdAlignNodes, pPars->ppFirst2Rpts);
				if ((NumPathNodesPE1 = CharacterisePath(QuerySeqLenPE1,pQuerySeqPE1,SortedPathIdxPE1, pPars->pAllocdAlignNodes, pPars->ppFirst2Rpts, &QueryPathEndOfsPE1, &TargPathEndOfsPE1)) == 0)
					continue;
				
				PathHiScorePE1 = pHeadNodePE1->HiScore;
				PathHiScorePE1 = (int)(((double)PathHiScorePE1 / ((int64_t)QuerySeqLenPE1 * m_ExactMatchScore)) * 255.0);
				if (PathHiScorePE1 >= MinPathHiScorePE1)
					{
					bStrandPE1 = pHeadNodePE1->FlgStrand ? false : true;
					TargPathStartOfsPE1 = pHeadNodePE1->TargSeqLoci;
					QueryPathStartOfsPE1 = pHeadNodePE1->QueryStartOfs;
					for (SortedPathIdxPE2 = 0; SortedPathIdxPE2 < min(NumHeadNodesPE2, (uint32_t)m_MaxPathsToReport + 1); SortedPathIdxPE2++)
						{
						pHeadNodePE2 = pPars->ppFirst2RptsPE2[SortedPathIdxPE2];
						ConsolidateNodes(pHeadNodePE2->FlgStrand,QuerySeqLenPE2, pQuerySeqPE2, SortedPathIdxPE2, pPars->pAllocdAlignNodesPE2, pPars->ppFirst2RptsPE2);
						if ((NumPathNodesPE2 = CharacterisePath(QuerySeqLenPE2, pQuerySeqPE2, SortedPathIdxPE2, pPars->pAllocdAlignNodesPE2, pPars->ppFirst2RptsPE2, &QueryPathEndOfsPE2, &TargPathEndOfsPE2)) == 0)
							continue;
						if(pHeadNodePE1->TargSeqID != pHeadNodePE2->TargSeqID ||		// both PE1 and PE2 must align to same target and be antisense to each other 
							bStrandPE1 == (bStrandPE2 = (pHeadNodePE2->FlgStrand ? false : true)))
							continue;
						PathHiScorePE2 = pHeadNodePE2->HiScore;
						PathHiScorePE2 = (int)(((double)PathHiScorePE2 / ((int64_t)QuerySeqLenPE2 * m_ExactMatchScore)) * 255.0);
						if (PathHiScorePE2 >= MinPathHiScorePE2)
							{
							// select that pair which are shortest distance apart and having maximal combined scores
							// very rough pairing used because one or both of the ends may only be partially aligned and so true start/ends are unknown
							TargPathStartOfsPE2 = pHeadNodePE2->TargSeqLoci;
							DeltaPathStarts = abs((INT64)(UINT64)TargPathStartOfsPE2 - (INT64)(UINT64)TargPathStartOfsPE1);
							if(DeltaPathStarts > m_MaxInsertLen + 100)
								continue;
							if((uint32_t)DeltaPathStarts > BestDeltaPathStart)
								continue;
							CombinedScore = PathHiScorePE1 + PathHiScorePE2;
							if(BestCombinedScore > 0 && CombinedScore < BestCombinedScore)
								continue;
							if(CombinedScore >= BestCombinedScore && DeltaPathStarts <= BestDeltaPathStart)
								{
								BestCombinedScore = CombinedScore;
								BestDeltaPathStart = (uint32_t)DeltaPathStarts;
								BestSortedPathIdxPE1 = SortedPathIdxPE1;
								BestSortedPathIdxPE2 = SortedPathIdxPE2;
								}
							}
						}
					}
				}
			if (BestCombinedScore > 0)
				{
				NumPathNodesPE1 = CharacterisePath(QuerySeqLenPE1, pQuerySeqPE1, BestSortedPathIdxPE1, pPars->pAllocdAlignNodes, pPars->ppFirst2Rpts, &QueryPathEndOfsPE1, &TargPathEndOfsPE1);
				NumPathNodesPE2 = CharacterisePath(QuerySeqLenPE2, pQuerySeqPE2, BestSortedPathIdxPE2, pPars->pAllocdAlignNodesPE2, pPars->ppFirst2RptsPE2, &QueryPathEndOfsPE2, &TargPathEndOfsPE2);
				pHeadNodePE1 = pPars->ppFirst2Rpts[BestSortedPathIdxPE1];
				pHeadNodePE2 = pPars->ppFirst2RptsPE2[BestSortedPathIdxPE2];
				bStrandPE1 = pHeadNodePE1->FlgStrand ? true : false;
				bStrandPE2 = pHeadNodePE2->FlgStrand ? true : false;
				TargPathStartOfsPE1 = pHeadNodePE1->TargSeqLoci;
				TargPathStartOfsPE2 = pHeadNodePE2->TargSeqLoci;
				QueryPathStartOfsPE1 = pHeadNodePE1->QueryStartOfs;
				QueryPathStartOfsPE2 = pHeadNodePE2->QueryStartOfs;
				PathHiScorePE1 = pHeadNodePE1->HiScore;
				PathHiScorePE2 = pHeadNodePE2->HiScore;
				TargSeqLenPE1 = TargSeqLenPE2 = m_pSfxArray->GetSeqLen(pHeadNodePE1->TargSeqID);
				m_pSfxArray->GetIdentName(pHeadNodePE1->TargSeqID, sizeof(szTargName), szTargName);
				strcpy(szTargNamePE2, szTargName);
				PathHiScorePE1 = (int)(((double)PathHiScorePE1 / ((double)QuerySeqLenPE1 * m_ExactMatchScore)) * 255.0);
				PathHiScorePE2 = (int)(((double)PathHiScorePE2 / ((double)QuerySeqLenPE2 * m_ExactMatchScore)) * 255.0);
				SAMFlags = cSAMFlgReadPaired | cSAMFlgReadPairMap | cSAMFlgPE1;
				if(bStrandPE1 == true)
					SAMFlags |= cSAMFlgAS;
				else
					SAMFlags |= cSAMFlgMateAS;
				if (NumQueryPathsRprtdPE1 > 0)
					SAMFlags |= cSAMFlgNotPrimary;
				ReportAsSAM(SAMFlags,						// use as reported SAM flags
					bStrandPE1 == true ? '-' : '+',			// query sequence strand, '+' or '-'
					PathHiScorePE1,							// score for this path
					szQuerySeqIdent,						// identifies this query sequence
					QuerySeqLenPE1,							// Query sequence length
					pQuerySeqPE1,							// the query sequence as etSeqBase's
					QueryPathStartOfsPE1,					// Alignment start offset in query, 0..QSize-1
					QueryPathEndOfsPE1,						// Alignment end position in query
					szTargName,								// aligning to this target
					TargSeqLenPE1,							// Target sequence size
					TargPathStartOfsPE1 + 1,				// target alignment starts at this starting loci, 1..TargSeqLenPE1
					TargPathEndOfsPE1 + 1,					// ending at this loci (inclusive) 
					'=',									// Reference sequence name of the primary alignment of the NEXT read in the template, '*' if unknown, '=' if PE and both ends align to same reference
					TargPathStartOfsPE2 + 1,				// 1-based Position of the primary alignment of the NEXT read in the template. Set as 0 when the information is unavailable
					(TargPathEndOfsPE2 + 1) - TargPathStartOfsPE1,// signed template length, If all segments are mapped to the same reference, the signed observed template length equals the number of bases from the leftmost mapped base to the rightmost mapped base
					NumPathNodesPE1,						// number of alignment nodes in alignment path
					BestSortedPathIdxPE1,					// alignment path starts from this node - 0 based 
					pPars->pAllocdAlignNodes,				// alignment nodes
					pPars->ppFirst2Rpts);					// allocated to hold ptrs to alignment nodes which are marked as being FlgFirst2tRpt
				NumQueryPathsRprtdPE1 += 1;

				SAMFlags = cSAMFlgReadPaired | cSAMFlgReadPairMap | cSAMFlgPE2;
				if (bStrandPE2 == true)
					SAMFlags |= cSAMFlgAS;
				else
					SAMFlags |= cSAMFlgMateAS;
				if (NumQueryPathsRprtdPE2 > 0)
					SAMFlags |= cSAMFlgNotPrimary;
				ReportAsSAM(SAMFlags,						// use as reported SAM flags
					bStrandPE2 == true ? '-' : '+',			// query sequence strand, '+' or '-'
					PathHiScorePE2,							// score for this path
					szQuerySeqIdentPE2,						// identifies this query sequence
					QuerySeqLenPE2,							// Query sequence length
					pQuerySeqPE2,							// the query sequence as etSeqBase's
					QueryPathStartOfsPE2,					// Alignment start offset in query, 0..QSize-1
					QueryPathEndOfsPE2,						// Alignment end position in query
					szTargNamePE2,							// aligning to this target
					TargSeqLenPE2,							// Target sequence size
					TargPathStartOfsPE2 + 1,				// target alignment starts at this starting loci, 1..TargSeqLenPE2
					TargPathEndOfsPE2 + 1,					// target alignment ending at this offset (inclusive)
					'=',									// Reference sequence name of the primary alignment of the NEXT read in the template, '*' if unknown, '=' if PE and both ends align to same reference
					TargPathStartOfsPE1 + 1,				// 1-based Position of the primary alignment of the NEXT read in the template. Set as 0 when the information is unavailable
					(TargPathEndOfsPE2 + 1) - TargPathStartOfsPE1,// signed template length, If all segments are mapped to the same reference, the signed observed template length equals the number of bases from the leftmost mapped base to the rightmost mapped base
					NumPathNodesPE2,						// number of alignment nodes in alignment path
					BestSortedPathIdxPE2,					// alignment path starts from this node - 0 based
					pPars->pAllocdAlignNodesPE2,			// alignment nodes
					pPars->ppFirst2RptsPE2);				// allocated to hold ptrs to alignment nodes which are marked as being FlgFirst2tRpt
				NumQueryPathsRprtdPE2 += 1;
				bQueryAtLeast1Path = true;
				}
			}
		if (NumQueryPathsRprtdPE1 || NumQueryPathsRprtdPE2)
			{
			AcquireSerialise();
			m_ReportedPaths += (NumQueryPathsRprtdPE1 + NumQueryPathsRprtdPE2);
			m_QueriesPaths++;
			ReleaseSerialise();
			continue;
			}
		ReportNonAligned(szQuerySeqIdent, QuerySeqLenPE1, pQuerySeqPE2, szQuerySeqIdentPE2, QuerySeqLenPE2, pQuerySeqPE2);
		bQueryAtLeast1Path = true;
		}
	while (!bQueryAtLeast1Path);
	delete pQuerySeqPE1;
	delete pQuerySeqPE2;
	}
gDiagnostics.DiagOut(eDLInfo, gszProcName, "Progress: Thread %d completed", pPars->ThreadIdx);
return(0);
}


int
CBlitz::ProcAlignQuerySeqs(tsThreadQuerySeqsPars *pPars) 
{
int Rslt;
int MaxIter;
int NumQueryPathsRprtd;
int NumQueriesProc;
int NumMatches;
int QuerySeqLen;
int SeqID;
uint8_t* pQuerySeq;
char szQuerySeqIdent[cMaxQuerySeqIdentLen + 1];
uint32_t NumHeadNodes;
int MinPathScore;
NumQueriesProc = 0;
while((Rslt = DequeueQuerySeq(sizeof(szQuerySeqIdent),&SeqID,szQuerySeqIdent,&QuerySeqLen,&pQuerySeq))==1)
	{
	AcquireSerialise();
	m_NumQueriesProc += 1;
	ReleaseSerialise();
	NumQueriesProc += 1;
	NumQueryPathsRprtd = 0;
	MaxIter = m_MaxIter;
	NumMatches = m_pSfxArray->LocateQuerySeqs(SeqID,pQuerySeq,QuerySeqLen,m_CoreLen,m_CoreDelta,m_AlignStrand,pPars->NumAllocdAlignNodes,pPars->pAllocdAlignNodes,m_MaxIter, m_ExactMatchScore, m_MismatchScore);
	if(NumMatches)
		{
		if(NumMatches > 1)	// sorting by TargSeqID.QueryID.FlgStrand.TargStartOfs.QueryStartOfs
			qsort(pPars->pAllocdAlignNodes,NumMatches,sizeof(tsQueryAlignNodes),SortQueryAlignNodes);
		if (m_MinPathScore > 0)
			MinPathScore = m_MinPathScore;
		else
			MinPathScore = ((QuerySeqLen - 5) * m_ExactMatchScore) / 3;

		NumHeadNodes = IdentifyHighScorePaths(MinPathScore, m_MaxPathsToReport, QuerySeqLen, NumMatches, pPars->pAllocdAlignNodes, pPars->ppFirst2Rpts);

		NumQueryPathsRprtd = Report(MinPathScore,m_MaxPathsToReport,szQuerySeqIdent,QuerySeqLen,pQuerySeq,NumMatches,pPars->pAllocdAlignNodes,pPars->ppFirst2Rpts);
		AcquireSerialise();
		if(NumQueryPathsRprtd > 0)
			m_QueriesPaths += 1;
		ReleaseSerialise();
		}
	delete pQuerySeq;
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: Thread %d completed, processed %d query sequences",pPars->ThreadIdx,NumQueriesProc);
return(0);
}

// locates and identifies the highest scoring smith-waterman path
// returns the highest score for all paths explored starting at the current node
// NOTE: recursive!!!
uint32_t											// returned best score for paths starting at pAlignNodes[ExploreNodeIdx]
CBlitz::HighScoreSW(uint32_t QueryLen,			// query length
			uint32_t TargSeqLen,					// targeted sequence length
 			bool bStrand,						// scoring for series on this strand - false if sense, true if antisense
			uint32_t ExploreNodeIdx,				// node to be explored for maximally scored path
			uint32_t NumNodes,					// total number of alignment nodes 
			tsQueryAlignNodes *pAlignNodes)	// alignment nodes
{
uint32_t BestHighScore;
uint32_t NodeIdx;
uint32_t CurNodeScore;
uint32_t GapScore;
uint32_t GapLen;
uint32_t TargGapLen;
uint32_t QueryGapLen;
uint32_t PutHighScore;
tsQueryAlignNodes *pCurNode;
tsQueryAlignNodes *pExploreNode;
pCurNode = &pAlignNodes[ExploreNodeIdx - 1];

if(pCurNode->FlgScored)
	return(pCurNode->HiScore);

CurNodeScore = pCurNode->AlignLen * m_ExactMatchScore;	// score favors exact matches
if(pCurNode->NumMismatches)								// score heavily penalises for mismatches
	{
	if((pCurNode->NumMismatches * m_MismatchScore) >= CurNodeScore)
		CurNodeScore = 0;
	else
		CurNodeScore -= (pCurNode->NumMismatches * m_MismatchScore);
	}

pCurNode->HiScore = CurNodeScore;
pExploreNode = pAlignNodes;
BestHighScore = 0;
for(NodeIdx = 1; NodeIdx <= NumNodes; NodeIdx++,pExploreNode++)
	{
	if(ExploreNodeIdx == NodeIdx)		// skip self node
		continue;
	if(pExploreNode->TargSeqID != pCurNode->TargSeqID)  // skip if node to explore is for different target
		continue;
	if(pExploreNode->Flg2Rpt || pExploreNode->FlgStrand != (bStrand ? 1 : 0))		// skip if already path to be reported, or if not requested strand
		continue;

	if(pExploreNode->TargSeqLoci < (pCurNode->TargSeqLoci + pCurNode->AlignLen - cMaxOverlapFloat))	// allowing for possible overlaps on the target sequence
		continue;
	if(pExploreNode->TargSeqLoci > (pCurNode->TargSeqLoci + pCurNode->AlignLen + cGapMaxLength))		// if gap too large then assuming not on same path
		continue;
	if(pExploreNode->QueryStartOfs < (pCurNode->QueryStartOfs + pCurNode->AlignLen - cMaxOverlapFloat))   // allowing for possible overlaps on on the query sequence
		continue;
	QueryGapLen = abs((int)(pExploreNode->QueryStartOfs - (pCurNode->QueryStartOfs + pCurNode->AlignLen)));
	TargGapLen = abs((int)(pExploreNode->TargSeqLoci - (pCurNode->TargSeqLoci + pCurNode->AlignLen)));
	GapLen = (int)sqrt(((double)QueryGapLen * QueryGapLen) + ((double)TargGapLen * TargGapLen));
	GapScore = 1 + ((GapLen / 10) * cGapExtendCost);
	if(GapScore > cGapExtendCostLimit)
		GapScore = cGapExtendCostLimit;
	GapScore += m_GapOpenScore;

	if(pExploreNode->FlgScored)			// if node being explored already scored then accept that score
		PutHighScore = pExploreNode->HiScore;
	else
		PutHighScore = HighScoreSW(QueryLen,TargSeqLen,bStrand,NodeIdx,NumNodes,pAlignNodes);	// node being explored not scored so recursively score starting from that node
	PutHighScore += CurNodeScore;
	if(GapScore >= PutHighScore)
		PutHighScore = 0;
	else
		PutHighScore -= GapScore;
	if(PutHighScore > BestHighScore)
		{
		BestHighScore = PutHighScore;
		pCurNode->HiScore = BestHighScore; 
		pCurNode->HiScorePathNextIdx = NodeIdx;
		}
	}
pCurNode->FlgScored = 1;
return(pCurNode->HiScore);
}

// expectation is that nodes will have been sorted in TargSeqID.QueryID.FlgStrand.QueryStartOfs.TargStartOfs ascending order
// essentially is dynamic programming (a.k smith-waterman) using nodes instead of the sequences as the nodes already contain
// matching + mismatches along the diagonals
uint32_t		// returns count of paths meeting scoring threshold
CBlitz::IdentifyHighScorePaths(uint32_t QueryLen,	// query length
			uint32_t TargSeqLen,					// targeted sequence length
			bool bStrand,						// reporting for paths on this strand - false if sense, true if antisense
			uint32_t NumNodes,					// reporting on this number of nodes starting from StartNodeIdx
			uint32_t StartNodeIdx,				// report for nodes starting at this node index (1..NumNodes) which is expected to be the first alignment node of a new target sequence
			tsQueryAlignNodes *pAlignNodes,		// alignment nodes
			uint32_t MinPathScore,				// only interested in paths having at least this score
			uint32_t  MaxPathsToReport)			// report at most this many alignment paths for any query
{
uint32_t PutBestHighScore;
uint32_t BestHighScore;
uint32_t BestHighScoreNodeIdx;
uint32_t CurNodeIdx;
uint32_t NumPutPaths2Rpt;
tsQueryAlignNodes *pCurNode;
tsQueryAlignNodes *pAlignSeqNodes;

if(MinPathScore < 25)			// have to have some minimum!
	MinPathScore = 25;
pAlignSeqNodes =  &pAlignNodes[StartNodeIdx-1];
NumPutPaths2Rpt = 0;
do {
	pCurNode = pAlignSeqNodes;
	for(CurNodeIdx = 1; CurNodeIdx <= NumNodes; CurNodeIdx++,pCurNode++)
		{
		if(pCurNode->FlgStrand != (bStrand ? 1 : 0) || pCurNode->Flg2Rpt == 1) // path scoring is strand specific, and retain scores and flags if this node already in a putative path to be reported
			continue;
		pCurNode->FlgFirst2tRpt = 0;
		pCurNode->Flg2Rpt = 0;
		pCurNode->FlgScored = 0;
		pCurNode->HiScore = 0;
		pCurNode->HiScorePathNextIdx = 0;
		}
	BestHighScore = 0; // temp change back to: MinPathScore - 1;
	BestHighScoreNodeIdx = 0;
	pCurNode = pAlignSeqNodes;
	for(CurNodeIdx = 1; CurNodeIdx <= NumNodes; CurNodeIdx++,pCurNode++)
		{
		if(pCurNode->Flg2Rpt || pCurNode->FlgStrand != (bStrand ? 1 : 0))	// once a node marked for reporting in a given path then can't be reported in any other path 
			continue;
		// get best score for any path originating from pCurNode, and retain if higher than any other originating nodes highest path score
		if((PutBestHighScore=HighScoreSW(QueryLen,TargSeqLen,bStrand,CurNodeIdx,NumNodes,pAlignSeqNodes)) > BestHighScore)
			{
			BestHighScore = PutBestHighScore;		// best thus far
			BestHighScoreNodeIdx = CurNodeIdx;		// record which node is the originating node
			}
		} 

	// if best path score meets the minimum required then mark 1st node as being the first and all nodes on path to be putatively reported
	if(BestHighScore >= MinPathScore)
		{
		uint32_t CurPathNodeIdx = BestHighScoreNodeIdx;
		uint32_t PathAlignedLen = 0;

		do {
			pCurNode = &pAlignSeqNodes[CurPathNodeIdx-1];
			PathAlignedLen += pCurNode->AlignLen;
			CurPathNodeIdx = pCurNode->HiScorePathNextIdx;
			}
		while(CurPathNodeIdx != 0);
		if(((PathAlignedLen * 100) / QueryLen) >= (uint32_t)m_QueryLenAlignedPct)
			{
			CurPathNodeIdx = BestHighScoreNodeIdx;
			do {
				pCurNode = &pAlignSeqNodes[CurPathNodeIdx-1];
				if(CurPathNodeIdx == BestHighScoreNodeIdx)
					pCurNode->FlgFirst2tRpt=1;			// 1st node on path
				else
					pCurNode->FlgFirst2tRpt = 0;
				pCurNode->Flg2Rpt = 1;					// marking all nodes on path as being reportable
				CurPathNodeIdx = pCurNode->HiScorePathNextIdx;
				if(CurPathNodeIdx)
					pCurNode->HiScorePathNextIdx = CurPathNodeIdx + StartNodeIdx - 1;
				}
			while(CurPathNodeIdx != 0);
			NumPutPaths2Rpt += 1;
			}
		else
			BestHighScore = 0;
		}
	}
while(BestHighScore >= MinPathScore && NumPutPaths2Rpt < MaxPathsToReport);

return(NumPutPaths2Rpt);
}


int				// number of blocks processed
CBlitz::BlocksAlignStats(uint32_t *pMatches,			// returned number of bases that match that aren't repeats
					uint32_t *pmisMatches,			// returned number of bases that don't match
					uint32_t *prepMatches,			// returned number of bases that match but are part of repeats
					uint32_t *pnCount,				// returned number of 'N' bases
				char  Strand,						// query sequence strand, '+' or '-'
				uint8_t *pQuerySeq,					// the query sequence
				uint32_t qSize,						// Query sequence size
				uint32_t TargSeqID,					// CSfxArray sequence identifier
				uint32_t tSize,						// Target sequence size 
				uint32_t TargPathStartOfs,			// at this starting offset
				uint32_t TargPathEndOfs,				// ending at this offset (inclusive)
				uint32_t NumPathNodes,				// number of alignment nodes in alignment path
				int SortedPathIdx,
				tsQueryAlignNodes *pAlignNodes,		// alignment nodes
				tsQueryAlignNodes **ppFirst2Rpts) // allocated to hold ptrs to alignment nodes which are marked as being FlgFirst2tRpt
{
uint32_t Matches;
uint32_t misMatches;
uint32_t nCount;
uint32_t TotMatches;
uint32_t TotMisMatches;
uint32_t TotNCount;
uint32_t Idx;
int NumBlocks;
uint32_t MaxBlockSize;
tsQueryAlignNodes *pCurNode;
tsQueryAlignNodes *pHeadNode;
uint8_t *pTSeq;
uint8_t *pQSeq;
etSeqBase *pTargBase;
etSeqBase *pQueryBase;

if(pMatches != NULL)
	*pMatches = 0;
if(pmisMatches != NULL)
	*pmisMatches = 0;
if(prepMatches != NULL)
	*prepMatches = 0;
if(pnCount != NULL)
	*pnCount = 0;

// determine maximal sized block to allocate for
MaxBlockSize = 0;
pHeadNode = ppFirst2Rpts[SortedPathIdx];
// determine maximun sized block and allocate mem to hold both query and target seq of that maximum 
pCurNode = pHeadNode;
do {
	if(pCurNode->AlignLen > MaxBlockSize)
		MaxBlockSize = pCurNode->AlignLen;
	if(pCurNode->HiScorePathNextIdx > 0)
		pCurNode = &pAlignNodes[pCurNode->HiScorePathNextIdx-1];
	else
		pCurNode = NULL;
	}
while(pCurNode != NULL);

if((pTSeq = new uint8_t [MaxBlockSize])==NULL)
	return(eBSFerrMem);
if((pQSeq = new uint8_t [MaxBlockSize])==NULL)
	{
	delete pTSeq;
	return(eBSFerrMem);
	}

// iterate each block accruing matches, mismatches and number of indeterminates in that block to derive totals
NumBlocks = 0;
TotMatches = 0;
TotMisMatches = 0;
TotNCount =  0;
pCurNode = pHeadNode;
do {
	Matches = 0;
	misMatches = 0;
	nCount = 0;
	NumBlocks += 1;
	m_pSfxArray->GetSeq(TargSeqID,pCurNode->TargSeqLoci,pTSeq,pCurNode->AlignLen);
	if(Strand == '-')
		{
		memcpy(pQSeq,&pQuerySeq[qSize - (pCurNode->QueryStartOfs + pCurNode->AlignLen)],pCurNode->AlignLen);
		CSeqTrans::ReverseComplement(pCurNode->AlignLen,pQSeq);
		}
	else
		memcpy(pQSeq,&pQuerySeq[pCurNode->QueryStartOfs],pCurNode->AlignLen);

	pTargBase = pTSeq;
	pQueryBase = pQSeq;

	for(Idx = 0; Idx < pCurNode->AlignLen; Idx++,pTargBase++,pQueryBase++)
		{
		if((*pQueryBase & 0x07) > eBaseT ||  (*pTargBase & 0x07) > eBaseT)
			{
			nCount += 1;
			misMatches += 1;
			continue;
			}
		if(*pQueryBase == *pTargBase)
			Matches += 1;
		else
			misMatches += 1;
		}

	TotMatches += Matches;
	TotMisMatches += misMatches;
	TotNCount += nCount;
	if(pCurNode->HiScorePathNextIdx > 0)
		pCurNode = &pAlignNodes[pCurNode->HiScorePathNextIdx-1];
	else
		pCurNode = NULL;
	}
while(pCurNode != NULL);

delete pTSeq;
delete pQSeq;
if(pMatches != NULL)
	*pMatches = TotMatches;
if(pmisMatches != NULL)
	*pmisMatches = TotMisMatches;
if(pnCount != NULL)
	*pnCount = TotNCount;
return(NumBlocks);
}


// reporting alignment as SQLite PSL format
int 
CBlitz::ReportAsSQLitePSL(uint32_t Matches,				// Number of bases that match that aren't repeats
					uint32_t misMatches,			// Number of bases that don't match
					uint32_t repMatches,			// Number of bases that match but are part of repeats
					uint32_t nCount,				// Number of 'N' bases
					uint32_t	qNumInsert,			// Number of inserts in query
					uint32_t qBaseInsert,			// Number of bases inserted in query
					uint32_t tNumInsert,			// Number of inserts in target
					uint32_t tBaseInsert,			// Number of bases inserted in target
					char  Strand,				// query sequence strand, '+' or '-'
					char *pszQuerySeqIdent,     // this query sequence
					uint32_t qSize,				// Query sequence size
					uint32_t qStart,				// Alignment start position in query
					uint32_t qEnd,				// Alignment end position in query
					char *pszTargName,			// aligning to this target
					uint32_t tSize,				// Target sequence size 
					uint32_t TargPathStartOfs,	// at this starting offset
					uint32_t TargPathEndOfs,		// ending at this offset (inclusive)
					uint32_t NumPathNodes,		// number of alignment nodes in alignment path
					int SortedPathIdx,
					tsQueryAlignNodes *pAlignNodes,		// alignment nodes
					tsQueryAlignNodes **ppFirst2Rpts) // allocated to hold ptrs to alignment nodes which are marked as being FlgFirst2tRpt
{
int Score;				// Alignment score (using Blat pslScore() function)
int Identity;           // Alignment identity (using Blat 100.0 - pslCalcMilliBad(psl, TRUE) * 0.1)
int StrandQStart;
int StrandQEnd;
char szStrand[2];
uint32_t *pBlockLens;
uint32_t *pQueryBlockStarts;
uint32_t *pTargBlockStarts;
uint32_t BlockLens[5000];
uint32_t QueryBlockStarts[5000];
uint32_t TargBlockStarts[5000];

tsQueryAlignNodes *pCurNode;

StrandQStart = Strand == '+' ? qStart : qSize - (qEnd + 1),
StrandQEnd = Strand == '+' ? qEnd+1 : qSize - qStart,

szStrand[0] = Strand;
szStrand[1] = '\0';
// block sizes and starts
pCurNode = ppFirst2Rpts[SortedPathIdx];
pBlockLens = BlockLens;
pQueryBlockStarts = QueryBlockStarts;
pTargBlockStarts = TargBlockStarts;
do {
	*pBlockLens++ = pCurNode->AlignLen;
	*pQueryBlockStarts++ = pCurNode->QueryStartOfs;
	*pTargBlockStarts++ = pCurNode->TargSeqLoci;
	if(pCurNode->HiScorePathNextIdx > 0)
		pCurNode = &pAlignNodes[pCurNode->HiScorePathNextIdx-1];
	else
		pCurNode = NULL;
	}
while(pCurNode != NULL);

Score = m_pSQLitePSL->pslScore(Matches,misMatches,repMatches,qNumInsert,tNumInsert,szStrand,tSize,TargPathStartOfs,TargPathEndOfs+1,NumPathNodes,(int *)BlockLens,(int *)QueryBlockStarts,(int *)TargBlockStarts);
double pslIdent = (double)m_pSQLitePSL->pslCalcMilliBad(Matches,misMatches,repMatches,qNumInsert,tNumInsert,qSize,StrandQStart,StrandQEnd,szStrand,tSize,TargPathStartOfs,TargPathEndOfs,NumPathNodes,(int *)BlockLens,(int *)QueryBlockStarts,(int *)TargBlockStarts,true); 
Identity = (int)(100.0 - pslIdent * 0.1);
AcquireLock(true);

m_pSQLitePSL->AddAlignment(m_ExprID,Score,Identity,Matches,misMatches,repMatches,nCount,qNumInsert,qBaseInsert,tNumInsert,tBaseInsert,szStrand,pszQuerySeqIdent,qSize,StrandQStart,StrandQEnd,pszTargName,tSize,TargPathStartOfs,TargPathEndOfs+1,NumPathNodes,(int *)BlockLens,(int *)QueryBlockStarts,(int *)TargBlockStarts);

m_ReportedPaths += 1;				
ReleaseLock(true);
return(NumPathNodes);
}

// reporting alignment as PSL format
int 
CBlitz::ReportAsPSL(uint32_t Matches,				// Number of bases that match that aren't repeats
					uint32_t misMatches,			// Number of bases that don't match
					uint32_t repMatches,			// Number of bases that match but are part of repeats
					uint32_t nCount,				// Number of 'N' bases
					uint32_t	qNumInsert,			// Number of inserts in query
					uint32_t qBaseInsert,			// Number of bases inserted in query
					uint32_t tNumInsert,			// Number of inserts in target
					uint32_t tBaseInsert,			// Number of bases inserted in target
					char  Strand,				// query sequence strand, '+' or '-'
					char *pszQuerySeqIdent,     // this query sequence
					uint32_t qSize,				// Query sequence size
					uint32_t qStart,				// Alignment start position in query
					uint32_t qEnd,				// Alignment end position in query
					char *pszTargName,			// aligning to this target
					uint32_t tSize,				// Target sequence size 
					uint32_t TargPathStartOfs,	// at this starting offset
					uint32_t TargPathEndOfs,		// ending at this offset (inclusive)
					uint32_t NumPathNodes,		// number of alignment nodes in alignment path
					int SortedPathIdx,
					tsQueryAlignNodes *pAlignNodes,		// alignment nodes
					tsQueryAlignNodes **ppFirst2Rpts) // allocated to hold ptrs to alignment nodes which are marked as being FlgFirst2tRpt
{
tsQueryAlignNodes *pCurNode;

AcquireSerialise();
if(m_szLineBuffIdx > (cAlignRprtBufferSize * 9) / 10)
	{
	CUtility::SafeWrite(m_hOutFile,m_pszLineBuff,m_szLineBuffIdx);
	m_szLineBuffIdx = 0;
	}

m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%c\t%s\t%u\t%u\t%u\t%s\t%u\t%u\t%u\t%u\t",
					Matches,misMatches,repMatches,nCount,qNumInsert,qBaseInsert,tNumInsert,tBaseInsert,Strand,
					pszQuerySeqIdent,qSize,
					Strand == '+' ? qStart : qSize - (qEnd + 1),
					Strand == '+' ? qEnd+1 : qSize - qStart,
					pszTargName,tSize,TargPathStartOfs,TargPathEndOfs+1,NumPathNodes);

// block sizes
pCurNode = ppFirst2Rpts[SortedPathIdx];
do {
	m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"%u,",pCurNode->AlignLen);
	if(pCurNode->HiScorePathNextIdx > 0)
		pCurNode = &pAlignNodes[pCurNode->HiScorePathNextIdx-1];
	else
		pCurNode = NULL;
	}
while(pCurNode != NULL);
m_pszLineBuff[m_szLineBuffIdx++] = '\t';

// query starts
pCurNode = ppFirst2Rpts[SortedPathIdx];
do {
	m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"%u,", pCurNode->QueryStartOfs);
	if(pCurNode->HiScorePathNextIdx > 0)
		pCurNode = &pAlignNodes[pCurNode->HiScorePathNextIdx-1];
	else
		pCurNode = NULL;
	}
while(pCurNode != NULL);
m_pszLineBuff[m_szLineBuffIdx++] = '\t';

// target starts
pCurNode = ppFirst2Rpts[SortedPathIdx];
do {
	m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"%u,",pCurNode->TargSeqLoci);
	if(pCurNode->HiScorePathNextIdx > 0)
		pCurNode = &pAlignNodes[pCurNode->HiScorePathNextIdx-1];
	else
		pCurNode = NULL;
	}
while(pCurNode != NULL);
m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"\n");	
if(m_szLineBuffIdx > (cAlignRprtBufferSize * 9) / 10)
	{
	CUtility::SafeWrite(m_hOutFile,m_pszLineBuff,m_szLineBuffIdx);
	m_szLineBuffIdx = 0;
	}
m_ReportedPaths += 1;				
ReleaseSerialise();
return(NumPathNodes);
}

// reporting alignment as PSLX format
int 
CBlitz::ReportAsPSLX(uint32_t Matches,				// Number of bases that match that aren't repeats
					uint32_t misMatches,			// Number of bases that don't match
					uint32_t repMatches,			// Number of bases that match but are part of repeats
					uint32_t nCount,				// Number of 'N' bases
					uint32_t	qNumInsert,			// Number of inserts in query
					uint32_t qBaseInsert,			// Number of bases inserted in query
					uint32_t tNumInsert,			// Number of inserts in target
					uint32_t tBaseInsert,			// Number of bases inserted in target
					char  Strand,				// query sequence strand, '+' or '-'
					char *pszQuerySeqIdent,     // this query sequence
					uint32_t qSize,				// Query sequence size
					uint8_t *pQuerySeq,			// the query sequence
					uint32_t qStart,				// Alignment start position in query
					uint32_t qEnd,				// Alignment end position in query
					uint32_t TargSeqID,			// CSfxArray sequence identifier
					char *pszTargName,			// aligning to this target
					uint32_t tSize,				// Target sequence size 
					uint32_t TargPathStartOfs,	// at this starting offset
					uint32_t TargPathEndOfs,		// ending at this offset (inclusive)
					uint32_t NumPathNodes,		// number of alignment nodes in alignment path
					int SortedPathIdx,
					tsQueryAlignNodes *pAlignNodes,		// alignment nodes
					tsQueryAlignNodes **ppFirst2Rpts) // allocated to hold ptrs to alignment nodes which are marked as being FlgFirst2tRpt
{
uint32_t MaxBlockSize;
etSeqBase *pCurSeq;
tsQueryAlignNodes *pCurNode;

pCurSeq = NULL;
AcquireSerialise();
if(m_szLineBuffIdx > (cAlignRprtBufferSize * 9) / 10)
	{
	CUtility::SafeWrite(m_hOutFile,m_pszLineBuff,m_szLineBuffIdx);
	m_szLineBuffIdx = 0;
	}

m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%c\t%s\t%u\t%u\t%u\t%s\t%u\t%u\t%u\t%u\t",
					Matches,misMatches,repMatches,nCount,qNumInsert,qBaseInsert,tNumInsert,tBaseInsert,Strand,
					pszQuerySeqIdent,qSize,
					Strand == '+' ? qStart : qSize - (qEnd + 1),
					Strand == '+' ? qEnd+1 : qSize - qStart,
					pszTargName,tSize,TargPathStartOfs,TargPathEndOfs+1,NumPathNodes);

// block sizes
MaxBlockSize = 0;
pCurNode = ppFirst2Rpts[SortedPathIdx];
do {
	m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"%u,",pCurNode->AlignLen);
	if(pCurNode->AlignLen > MaxBlockSize)
		MaxBlockSize = pCurNode->AlignLen;
	if(pCurNode->HiScorePathNextIdx > 0)
		pCurNode = &pAlignNodes[pCurNode->HiScorePathNextIdx-1];
	else
		pCurNode = NULL;
	}
while(pCurNode != NULL);
m_pszLineBuff[m_szLineBuffIdx++] = '\t';

// query starts
pCurNode = ppFirst2Rpts[SortedPathIdx];
do {
	m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"%u,", pCurNode->QueryStartOfs);
	if(pCurNode->HiScorePathNextIdx > 0)
		pCurNode = &pAlignNodes[pCurNode->HiScorePathNextIdx-1];
	else
		pCurNode = NULL;
	}
while(pCurNode != NULL);
m_pszLineBuff[m_szLineBuffIdx++] = '\t';

// target starts
pCurNode = ppFirst2Rpts[SortedPathIdx];
do {
	m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"%u,",pCurNode->TargSeqLoci);
	if(pCurNode->HiScorePathNextIdx > 0)
		pCurNode = &pAlignNodes[pCurNode->HiScorePathNextIdx-1];
	else
		pCurNode = NULL;
	}
while(pCurNode != NULL);

// allocate to hold the maximally sized sequence block
if((pCurSeq = new uint8_t [10 + TargPathEndOfs - TargPathStartOfs])==NULL)	// inplace translation to ascii so allow a few extra for terminagting '\0'	
	{
	ReleaseSerialise();
	return(eBSFerrMem);
	}

// get query sequences for each block and report these
// query starts
pCurNode = ppFirst2Rpts[SortedPathIdx];
do {
	if((m_szLineBuffIdx + pCurNode->AlignLen) > (cAlignRprtBufferSize * 9) / 10)
		{
		CUtility::SafeWrite(m_hOutFile,m_pszLineBuff,m_szLineBuffIdx);
		m_szLineBuffIdx = 0;
		}
	if(Strand == '-')
		{
		memcpy(pCurSeq,&pQuerySeq[qSize - (pCurNode->QueryStartOfs + pCurNode->AlignLen)],pCurNode->AlignLen);
		CSeqTrans::ReverseComplement(pCurNode->AlignLen,pCurSeq);
		}
	else
		memcpy(pCurSeq,&pQuerySeq[pCurNode->QueryStartOfs],pCurNode->AlignLen);
	m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"%s,",CSeqTrans::MapSeq2Ascii(pCurSeq,pCurNode->AlignLen,(char *)pCurSeq));
	if(pCurNode->HiScorePathNextIdx > 0)
		pCurNode = &pAlignNodes[pCurNode->HiScorePathNextIdx-1];
	else
		pCurNode = NULL;
	}
while(pCurNode != NULL);
m_pszLineBuff[m_szLineBuffIdx++] = '\t';


// get target sequences for each block and report these
pCurNode = ppFirst2Rpts[SortedPathIdx];
do {
	if((m_szLineBuffIdx + pCurNode->AlignLen) > (cAlignRprtBufferSize * 9) / 10)
		{
		CUtility::SafeWrite(m_hOutFile,m_pszLineBuff,m_szLineBuffIdx);
		m_szLineBuffIdx = 0;
		}
	m_pSfxArray->GetSeq(TargSeqID,pCurNode->TargSeqLoci,pCurSeq,pCurNode->AlignLen);
	m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"%s,",CSeqTrans::MapSeq2Ascii(pCurSeq,pCurNode->AlignLen,(char *)pCurSeq));
	if(pCurNode->HiScorePathNextIdx > 0)
		pCurNode = &pAlignNodes[pCurNode->HiScorePathNextIdx-1];
	else
		pCurNode = NULL;
	}
while(pCurNode != NULL);

m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"\n");	
if(m_szLineBuffIdx > (cAlignRprtBufferSize * 9) / 10)
	{
	CUtility::SafeWrite(m_hOutFile,m_pszLineBuff,m_szLineBuffIdx);
	m_szLineBuffIdx = 0;
	}
m_ReportedPaths += 1;				
ReleaseSerialise();
if(pCurSeq!=NULL)
	delete pCurSeq;
return(NumPathNodes);
}

// Frequently there will be gaps between nodes and
// these gaps need to be closed through interpolation
uint32_t			// number of nodes after consolidation
CBlitz::ConsolidateNodes(bool bSense,				// if true then query sequence is aligning antisense to target
				uint32_t QueryLen,					// query length
				uint8_t* pQuerySeq,					// the query sequence as etSeqBase's
				uint32_t SortedPathIdx,				// alignment path starts from this node - 0 based 
				tsQueryAlignNodes* pAlignNodes,		// alignment nodes
				tsQueryAlignNodes** ppFirst2Rpts)	// allocated to hold ptrs to alignment nodes which are marked as being FlgFirst2tRpt
{
bool bCurBaseMatch;
bool bNxtBaseMatch;
int DeltaQuery;
int DeltaTarg;
uint32_t CurTargLoci;
uint32_t NxtTargLoci;
etSeqBase* pCurSeq;
etSeqBase CurTargBase;
etSeqBase* pNxtSeq;
etSeqBase NxtTargBase;
tsQueryAlignNodes* pCurNode;
tsQueryAlignNodes* pNxtNode;
uint32_t NumConsolidatedNodes;
NumConsolidatedNodes = 0;
pCurNode = ppFirst2Rpts[SortedPathIdx];			// starting from first node in path
if(pCurNode->HiScorePathNextIdx == 0)
	return(1);
if (bSense)
	CSeqTrans::ReverseComplement(QueryLen, (etSeqBase*)pQuerySeq);
while(pCurNode->HiScorePathNextIdx != 0)
	{
	NumConsolidatedNodes+=1;
	pNxtNode = &pAlignNodes[pCurNode->HiScorePathNextIdx - 1];
	DeltaQuery = pNxtNode->QueryStartOfs - (pCurNode->QueryStartOfs + pCurNode->AlignLen - 1);	// normally if at an InDel edge then difference should be 1 base
	DeltaTarg = pNxtNode->TargSeqLoci - (pCurNode->TargSeqLoci + pCurNode->AlignLen - 1);		// normally if just +1 then next node is a direct extension
	if (DeltaQuery < 0 || DeltaQuery > 1)
		{
		if (DeltaQuery > 1)			// mind {close} the gap ...
			{
			pCurSeq = &pQuerySeq[pCurNode->QueryStartOfs + pCurNode->AlignLen];
			pNxtSeq = &pQuerySeq[pNxtNode->QueryStartOfs - 1];
			CurTargLoci = pCurNode->TargSeqLoci + pCurNode->AlignLen;
			NxtTargLoci = pNxtNode->TargSeqLoci - 1;
			while (DeltaQuery-- > 1)
				{
				CurTargBase = m_pSfxArray->GetBase(pCurNode->TargSeqID, CurTargLoci);
				NxtTargBase = m_pSfxArray->GetBase(pNxtNode->TargSeqID, NxtTargLoci);
				bCurBaseMatch = *pCurSeq == CurTargBase;
				bNxtBaseMatch = *pNxtSeq == NxtTargBase;
				if(bCurBaseMatch)
					{
					pCurNode->AlignLen++;
					CurTargLoci++;
					pCurSeq++;
					}
				else
					{
					pNxtNode->AlignLen++;
					pNxtNode->QueryStartOfs--;
					pNxtNode->TargSeqLoci--;
					if (!bNxtBaseMatch)
						pNxtNode->NumMismatches++;
					NxtTargLoci--;
					pNxtSeq--;
					}
				}
			}
		else
			{
			pCurSeq = &pQuerySeq[pCurNode->QueryStartOfs + pCurNode->AlignLen - 1];
			pNxtSeq = &pQuerySeq[pNxtNode->QueryStartOfs];
			CurTargLoci = pCurNode->TargSeqLoci + pCurNode->AlignLen - 1;
			NxtTargLoci = pNxtNode->TargSeqLoci;
			while (DeltaQuery++ < 1)
				{
				CurTargBase = m_pSfxArray->GetBase(pCurNode->TargSeqID, CurTargLoci);
				NxtTargBase = m_pSfxArray->GetBase(pNxtNode->TargSeqID, NxtTargLoci);
				bCurBaseMatch = *pCurSeq == CurTargBase;
				bNxtBaseMatch = *pNxtSeq == NxtTargBase;
				if (bNxtBaseMatch)
					{
					pCurNode->AlignLen--;
					if (!bCurBaseMatch)
						pCurNode->NumMismatches--;
					CurTargLoci--;
					pCurSeq--;
					}
				else									
					{
					pNxtNode->AlignLen++;
					pNxtNode->QueryStartOfs++;
					pNxtNode->TargSeqLoci++;
					if (!bNxtBaseMatch)
						pNxtNode->NumMismatches--;
					NxtTargLoci++;
					pNxtSeq++;
					}
				}
			}
		}
	DeltaQuery = pNxtNode->QueryStartOfs - (pCurNode->QueryStartOfs + pCurNode->AlignLen - 1);	// normally if at an InDel edge then difference should be 1 base
	DeltaTarg = pNxtNode->TargSeqLoci - (pCurNode->TargSeqLoci + pCurNode->AlignLen - 1);		// normally if just +1 then next node is a direct extension
	if(DeltaQuery == 1 && DeltaTarg == 1)
		{
		pCurNode->HiScorePathNextIdx = pNxtNode->HiScorePathNextIdx;
		pCurNode->AlignLen += pNxtNode->AlignLen;
		pCurNode->NumMismatches += pNxtNode->NumMismatches;
		continue;
		}
	else
		pCurNode = pNxtNode;
	if (pCurNode->HiScorePathNextIdx == 0)
		NumConsolidatedNodes++;
	}
if (bSense)
	CSeqTrans::ReverseComplement(QueryLen, (etSeqBase *)pQuerySeq);
return(NumConsolidatedNodes);
}

int		// reporting alignment as SAM format
CBlitz::ReportAsSAM(uint32_t Flags,	// use as the reported SAM flags
	char  Strand,				// query sequence strand, '+' or '-'
	uint32_t PathScore,			// score for this path
	char *pszQuerySeqIdent,     // identifies this query sequence
	uint32_t qSize,				// Query sequence length
	uint8_t *pQuerySeq,			// the query sequence as etSeqBase's
	uint32_t qStart,				// Alignment start offset in query, 0..qSize-1
	uint32_t qEnd,				// Alignment end position in query
	char *pszTargName,			// aligning to this target
	uint32_t tSize,				// Target sequence size 
	uint32_t TargPathStartOfs,	// target alignment starts at this starting loci, 1..tSize
	uint32_t TargPathEndOfs,		// target alignment ending at this offset (inclusive)
	char RNEXT,					// Reference sequence name of the primary alignment of the NEXT read in the template, '*' if unknown, '=' if PE and both ends align to same reference
	uint32_t PNEXT,				// 1-based Position of the primary alignment of the NEXT read in the template. Set as 0 when the information is unavailable
	int TLEN,					// signed template length, If all segments are mapped to the same reference, the signed observed template length equals the number of bases from the leftmost mapped base to the rightmost mapped base
	uint32_t NumPathNodes,		// number of alignment nodes in alignment path
	int SortedPathIdx,			// alignment path starts from this node - 0 based 
	tsQueryAlignNodes *pAlignNodes,		// alignment nodes
	tsQueryAlignNodes **ppFirst2Rpts)  // allocated to hold ptrs to alignment nodes which are marked as being FlgFirst2tRpt
{
	etSeqBase* pCurSeq;
	tsQueryAlignNodes* pCurNode;
	etSeqBase CurSeq[cSAMtruncSeqLen + 1];
	char szLineBuff[cSAMtruncSeqLen + 2 * (cMaxQuerySeqIdentLen)+100];
	int BuffIdx;

	pCurSeq = NULL;

	BuffIdx = sprintf(szLineBuff, "%s\t%d\t%s\t%d\t%d\t", pszQuerySeqIdent, Flags, pszTargName, TargPathStartOfs, PathScore >= 254 ? 254 : PathScore);
	if (qStart > 0)
		BuffIdx += sprintf(&szLineBuff[BuffIdx], "%uS", qStart);

	uint32_t PrevQueryEndOfs = qStart;
	uint32_t PrevTargEndOfs = TargPathStartOfs;

	pCurNode = ppFirst2Rpts[SortedPathIdx];
	do {
		if (pCurNode->QueryStartOfs > PrevQueryEndOfs)
			BuffIdx += sprintf(&szLineBuff[BuffIdx], "%uI", pCurNode->QueryStartOfs - PrevQueryEndOfs);
		if (pCurNode->TargSeqLoci > PrevTargEndOfs)
		{
			if ((pCurNode->TargSeqLoci - PrevTargEndOfs) < 20)
				BuffIdx += sprintf(&szLineBuff[BuffIdx], "%uD", pCurNode->TargSeqLoci - PrevTargEndOfs);
			else
				BuffIdx += sprintf(&szLineBuff[BuffIdx], "%uN", pCurNode->TargSeqLoci - PrevTargEndOfs);	// larger insertion in target, treat as though insertion was due to intron spanning 
		}
		BuffIdx += sprintf(&szLineBuff[BuffIdx], "%uM", pCurNode->AlignLen);
		PrevQueryEndOfs = pCurNode->QueryStartOfs + pCurNode->AlignLen;
		PrevTargEndOfs = pCurNode->TargSeqLoci + pCurNode->AlignLen;
		if (pCurNode->HiScorePathNextIdx > 0)
			pCurNode = &pAlignNodes[pCurNode->HiScorePathNextIdx - 1];
		else
			pCurNode = NULL;
	} while (pCurNode != NULL);
	if ((qEnd + 1) < qSize)
		BuffIdx += sprintf(&szLineBuff[BuffIdx], "%uS", qSize - qEnd - 1);
	BuffIdx += sprintf(&szLineBuff[BuffIdx], "\t%c\t%u\t%d\t", RNEXT, PNEXT, TLEN);

	memcpy(CurSeq, pQuerySeq, qSize);
	if (Strand == '-')
		CSeqTrans::ReverseComplement(qSize, CurSeq);
	BuffIdx += sprintf(&szLineBuff[BuffIdx], "%s\t*\n", CSeqTrans::MapSeq2Ascii(CurSeq, qSize, (char*)CurSeq));

	AcquireSerialise();
	if (m_szLineBuffIdx >= (cAlignRprtBufferSize - (2 * BuffIdx)))
	{
		CUtility::SafeWrite(m_hOutFile, m_pszLineBuff, m_szLineBuffIdx);
		m_szLineBuffIdx = 0;
	}
	memcpy(&m_pszLineBuff[m_szLineBuffIdx], szLineBuff, BuffIdx);
	m_szLineBuffIdx += BuffIdx;
	ReleaseSerialise();
	return(NumPathNodes);
}

// reporting alignment as MAF format
int 
CBlitz::ReportAsMAF(int PathScore,				// score for this path
					uint32_t Matches,				// Number of bases that match that aren't repeats
					uint32_t misMatches,			// Number of bases that don't match
					uint32_t repMatches,			// Number of bases that match but are part of repeats
					uint32_t nCount,				// Number of 'N' bases
					uint32_t	qNumInsert,			// Number of inserts in query
					uint32_t qBaseInsert,			// Number of bases inserted in query
					uint32_t tNumInsert,			// Number of inserts in target
					uint32_t tBaseInsert,			// Number of bases inserted in target
					char  Strand,				// query sequence strand, '+' or '-'
					char *pszQuerySeqIdent,     // this query sequence
					uint32_t qSize,				// Query sequence size
					uint8_t *pQuerySeq,			// the query sequence
					uint32_t qStart,				// Alignment start position in query
					uint32_t qEnd,				// Alignment end position in query
					uint32_t TargSeqID,			// CSfxArray sequence identifier
					char *pszTargName,			// aligning to this target
					uint32_t tSize,				// Target sequence size 
					uint32_t TargPathStartOfs,	// at this starting offset
					uint32_t TargPathEndOfs,		// ending at this offset (inclusive)
					uint32_t NumPathNodes,		// number of alignment nodes in alignment path
					int SortedPathIdx,
					tsQueryAlignNodes *pAlignNodes,		// alignment nodes
					tsQueryAlignNodes **ppFirst2Rpts) // allocated to hold ptrs to alignment nodes which are marked as being FlgFirst2tRpt
{
uint32_t MaxBlockSize;
tsQueryAlignNodes *pCurNode;
uint8_t *pCurSeq;

// // allocate buffering for maximal sized alignment block
MaxBlockSize = 0;
pCurNode = ppFirst2Rpts[SortedPathIdx];
do {
	if(pCurNode->AlignLen > MaxBlockSize)
		MaxBlockSize = pCurNode->AlignLen;
	if(pCurNode->HiScorePathNextIdx > 0)
		pCurNode = &pAlignNodes[pCurNode->HiScorePathNextIdx-1];
	else
		pCurNode = NULL;
	}
while(pCurNode != NULL);

if((pCurSeq = new uint8_t [10 + MaxBlockSize])==NULL)	// inplace translation to ascii so allow a few extra for terminating '\0'	
	return(eBSFerrMem);

AcquireSerialise();
if(m_szLineBuffIdx > (cAlignRprtBufferSize * 9) / 10)
	{
	CUtility::SafeWrite(m_hOutFile,m_pszLineBuff,m_szLineBuffIdx);
	m_szLineBuffIdx = 0;
	}

pCurNode = ppFirst2Rpts[SortedPathIdx];
do {
	if((m_szLineBuffIdx + pCurNode->AlignLen) > (cAlignRprtBufferSize * 9) / 10)
		{
		CUtility::SafeWrite(m_hOutFile,m_pszLineBuff,m_szLineBuffIdx);
		m_szLineBuffIdx = 0;
		}	
	m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"\na\tscore=%d",pCurNode->HiScore);
	m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"\ns\t%s\t%u\t%u\t%c\t%u\t",pszTargName,pCurNode->TargSeqLoci,pCurNode->AlignLen,'+',tSize);
	m_pSfxArray->GetSeq(TargSeqID,pCurNode->TargSeqLoci,pCurSeq,pCurNode->AlignLen);
	m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"%s",CSeqTrans::MapSeq2Ascii(pCurSeq,pCurNode->AlignLen,(char *)pCurSeq));

	m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"\ns\t%s\t%u\t%u\t%c\t%u\t",pszQuerySeqIdent,pCurNode->TargSeqLoci,pCurNode->AlignLen,Strand,qSize);
	if(Strand == '-')
		{
		memcpy(pCurSeq,&pQuerySeq[qSize - (pCurNode->QueryStartOfs + pCurNode->AlignLen)],pCurNode->AlignLen);
		CSeqTrans::ReverseComplement(pCurNode->AlignLen,pCurSeq);
		}
	else
		memcpy(pCurSeq,&pQuerySeq[pCurNode->QueryStartOfs],pCurNode->AlignLen);

	m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"%s\n",CSeqTrans::MapSeq2Ascii(pCurSeq,pCurNode->AlignLen,(char *)pCurSeq));
	if(pCurNode->HiScorePathNextIdx > 0)
		pCurNode = &pAlignNodes[pCurNode->HiScorePathNextIdx-1];
	else
		pCurNode = NULL;
	}
while(pCurNode != NULL);
m_ReportedPaths += 1;
ReleaseSerialise();
if(pCurSeq != NULL)
	delete pCurSeq;

return(NumPathNodes);
}

// reporting KMer counts distributions as CSV
int	
CBlitz::ReportKMerDist(uint32_t KMerLen,// KMer count distributions are for this length K-mer 
	uint32_t MaxKMerCnts,				// max number of target occurrences which ranges from 0 to MaxKMerCnts inclusive
	uint32_t *pKmerOccsDist)			// array of length MaxKMerCnts + 1 which holds number of query K-mers having target occurrences
{
uint32_t DistIdx;
uint32_t *pOccsDist;
uint64_t CummulativeQKmers;
uint64_t CummulativeTKmers;
uint64_t TotalCummulativeQKmers;
uint64_t TotalCummulativeTKmers;
double PropCummulativeQKmers;
double PropCummulativeTKmers;

TotalCummulativeQKmers = 0;
TotalCummulativeTKmers = 0;
pOccsDist = pKmerOccsDist;
for (DistIdx = 0; DistIdx <= MaxKMerCnts; DistIdx++, pOccsDist++)
	{
	TotalCummulativeQKmers += (uint64_t)*pOccsDist;
	TotalCummulativeTKmers += DistIdx * (uint64_t)*pOccsDist;
	}


if (m_szLineBuffIdx > (cAlignRprtBufferSize * 9) / 10)
	{
	CUtility::SafeWrite(m_hOutFile, m_pszLineBuff, m_szLineBuffIdx);
	m_szLineBuffIdx = 0;
	}
CummulativeQKmers = 0;
CummulativeTKmers = 0;
for(DistIdx = 0; DistIdx <= MaxKMerCnts; DistIdx++, pKmerOccsDist++)
	{
	CummulativeQKmers += (uint64_t)*pKmerOccsDist;
	CummulativeTKmers += DistIdx * (uint64_t)*pKmerOccsDist;
	PropCummulativeQKmers = (double)CummulativeQKmers/(double)TotalCummulativeQKmers;
	PropCummulativeTKmers = (double)CummulativeTKmers / (double)TotalCummulativeTKmers;
#if _WIN32
	m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx], "%u,%u, %u, %Iu,%Iu,%f,%f\n",
#else
	m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx], "%u,%u, %u, %lu,%lu,%f,%f\n",
#endif
							DistIdx, *pKmerOccsDist, DistIdx == 0 ? *pKmerOccsDist : DistIdx * *pKmerOccsDist, CummulativeQKmers, CummulativeTKmers, PropCummulativeQKmers, PropCummulativeTKmers);
	if (m_szLineBuffIdx > (cAlignRprtBufferSize * 9) / 10)
		{
		CUtility::SafeWrite(m_hOutFile, m_pszLineBuff, m_szLineBuffIdx);
		m_szLineBuffIdx = 0;
		}
	}
if (m_szLineBuffIdx > 0)
	{
	CUtility::SafeWrite(m_hOutFile, m_pszLineBuff, m_szLineBuffIdx);
	m_szLineBuffIdx = 0;
	}
return(0);
}


// reporting alignment as BED format
int 
CBlitz::ReportAsBED(char *pszQuerySeqIdent,     // this query sequence
					char  Strand,				// query sequence strand, '+' or '-'
					uint32_t AlignScore,			// alignment has this score
					char *pszTargName,			// aligning to this target
					uint32_t NumPathNodes,		// number of alignment nodes in alignment path
					uint32_t TargPathStartOfs,	// at this starting offset
					uint32_t TargPathEndOfs,		// ending at this offset (inclusive)
					int SortedPathIdx,
					tsQueryAlignNodes *pAlignNodes,		// alignment nodes
					tsQueryAlignNodes **ppFirst2Rpts) // allocated to hold ptrs to alignment nodes which are marked as being FlgFirst2tRpt
{
tsQueryAlignNodes *pCurNode;
AlignScore = (uint32_t)(10 * sqrt(AlignScore));
if(AlignScore > 1000)
	AlignScore = 1000;
AcquireSerialise();
if(m_szLineBuffIdx > (cAlignRprtBufferSize * 9) / 10)
	{
	CUtility::SafeWrite(m_hOutFile,m_pszLineBuff,m_szLineBuffIdx);
	m_szLineBuffIdx = 0;
	}

m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"%s\t%d\t%d\t%s\t%d\t%c\t%d\t%d\t0\t%d\t",
					pszTargName,TargPathStartOfs,TargPathEndOfs+1,pszQuerySeqIdent,AlignScore,Strand,TargPathStartOfs,TargPathEndOfs+1,NumPathNodes);
pCurNode = ppFirst2Rpts[SortedPathIdx];
do {
	m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"%u,",pCurNode->AlignLen);
	if(pCurNode->HiScorePathNextIdx > 0)
		pCurNode = &pAlignNodes[pCurNode->HiScorePathNextIdx-1];
	else
		pCurNode = NULL;
	}
while(pCurNode != NULL);
m_pszLineBuff[m_szLineBuffIdx++] = '\t';
pCurNode = ppFirst2Rpts[SortedPathIdx];
do {
	m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"%u",pCurNode->TargSeqLoci - TargPathStartOfs);
	if(pCurNode->HiScorePathNextIdx > 0)
		{
		m_pszLineBuff[m_szLineBuffIdx++] = ',';
		pCurNode = &pAlignNodes[pCurNode->HiScorePathNextIdx-1];
		}
	else
		pCurNode = NULL;
	}
while(pCurNode != NULL);
m_szLineBuffIdx += sprintf(&m_pszLineBuff[m_szLineBuffIdx],"\n");
if(m_szLineBuffIdx > (cAlignRprtBufferSize * 9) / 10)
	{
	CUtility::SafeWrite(m_hOutFile,m_pszLineBuff,m_szLineBuffIdx);
	m_szLineBuffIdx = 0;
	}
m_ReportedPaths += 1;
ReleaseSerialise();
return(NumPathNodes);
}


uint32_t							// returned number of high scoring path head nodes
CBlitz::IdentifyHighScorePaths(int MinPathScore,			// only report paths having at least this minimum score
		int  MaxPathsToReport,		// report at most this many alignment paths for any query
		uint32_t QueryLen,			// query length
		uint32_t NumNodes,			// number of alignment nodes
		tsQueryAlignNodes *pAlignNodes,		// alignment nodes
		tsQueryAlignNodes **ppFirst2Rpts)	// allocated to hold ptrs to alignment nodes which are marked as being FlgFirst2tRpt
{
tsQueryAlignNodes *pCurrentNode;
uint32_t CurTargSeqID;
uint32_t CurTargMatchNodes;
uint32_t NodeIdx;
uint32_t CurTargSeqLen;
uint32_t CurNumHeadNodes;

pCurrentNode = pAlignNodes;
CurTargSeqID = 0;
CurTargMatchNodes = 0;
CurNumHeadNodes = 0;

// nodes will have been sorted in TargSeqID.QueryID.FlgStrand.QueryStartOfs.TargStartOfs ascending order
uint32_t StartTargNodeIdx;
uint32_t EndTargNodeIdx;
uint32_t NumPutPaths;

if(MaxPathsToReport < 1)
	MaxPathsToReport = 1;
NumPutPaths = 0;
CurTargSeqID = pCurrentNode->TargSeqID;
CurTargMatchNodes = 0;
for (StartTargNodeIdx = EndTargNodeIdx = 1; EndTargNodeIdx <= NumNodes; EndTargNodeIdx++, pCurrentNode++)
	{
	if (pCurrentNode->TargSeqID != CurTargSeqID || EndTargNodeIdx == NumNodes)				// onto a new target sequence or last node?
		{
		if (EndTargNodeIdx == NumNodes)
			CurTargMatchNodes += 1;
		CurTargSeqLen = m_pSfxArray->GetSeqLen(CurTargSeqID);
		if (m_AlignStrand != eALSCrick)
			NumPutPaths += IdentifyHighScorePaths(QueryLen, CurTargSeqLen, false, CurTargMatchNodes, StartTargNodeIdx, pAlignNodes, MinPathScore, MaxPathsToReport);	// sense/sense paths
		if (m_AlignStrand != eALSWatson)
			NumPutPaths += IdentifyHighScorePaths(QueryLen, CurTargSeqLen, true, CurTargMatchNodes, StartTargNodeIdx, pAlignNodes, MinPathScore, MaxPathsToReport);     // antisense/sense paths
		CurTargSeqID = pCurrentNode->TargSeqID;
		CurTargMatchNodes = 0;
		StartTargNodeIdx = EndTargNodeIdx;
		}
	CurTargMatchNodes += 1;
	}
if (NumPutPaths == 0)
	return(0);

pCurrentNode = pAlignNodes;
for (NodeIdx = 1; NodeIdx <= NumNodes; NodeIdx++, pCurrentNode++)
	{
	if (pCurrentNode->FlgFirst2tRpt)
		ppFirst2Rpts[CurNumHeadNodes++] = pCurrentNode;
	}

if (CurNumHeadNodes > 1)			// sort by highest scoring path descending
	qsort(ppFirst2Rpts, CurNumHeadNodes, sizeof(tsQueryAlignNodes *), SortHighScoreDescend);

return(min((uint32_t)MaxPathsToReport,CurNumHeadNodes));
}

uint32_t				// returns number of nodes in characterised path
CBlitz::CharacterisePath(uint32_t QueryLen,				// query length
					uint8_t* pQuerySeq,					// the query sequence as etSeqBase's
					uint32_t SortedPathIdx,				// index of paths head node
					tsQueryAlignNodes *pAlignNodes,		// alignment nodes
					tsQueryAlignNodes **ppFirst2Rpts,	// allocated to hold ptrs to alignment nodes which are marked as being FlgFirst2tRpt
					uint32_t *pQueryPathEndOfs,			// query path ends at this offset
					uint32_t *pTargPathEndOfs,			// target path ends at this offset
					uint32_t *pqNumInsert,				// qNumInsert, Number of inserts in query
					uint32_t *pqBaseInsert,				// qBaseInsert, Number of bases inserted in query
					uint32_t *ptNumInsert,				// tNumInsert, Number of inserts in target
					uint32_t *ptBaseInsert)				// tBaseInsert, Number of bases inserted in target
{
uint32_t CurPathNodeIdx;
uint32_t QueryPathStartOfs;
uint32_t TargPathStartOfs;
uint32_t QueryPathEndOfs;
uint32_t TargPathEndOfs;
uint32_t TotalMMs;
uint32_t NumPathNodes;
uint32_t	qNumInsert;			// Number of inserts in query
uint32_t qBaseInsert;			// Number of bases inserted in query
uint32_t tNumInsert;			// Number of inserts in target
uint32_t tBaseInsert;			// Number of bases inserted in target
int TargGap;
int QueryGap;
uint32_t TotalAlignLen;
bool bStrand;
int PathHiScore;
uint32_t TargSeqLen;
char szTargName[100];
tsQueryAlignNodes *pCurNode;

pCurNode = ppFirst2Rpts[SortedPathIdx];
bStrand = pCurNode->FlgStrand ? true : false;
PathHiScore = pCurNode->HiScore;
m_pSfxArray->GetIdentName(pCurNode->TargSeqID, sizeof(szTargName), szTargName);
TargSeqLen = m_pSfxArray->GetSeqLen(pCurNode->TargSeqID);
QueryPathStartOfs = pCurNode->QueryStartOfs;
TargPathStartOfs = pCurNode->TargSeqLoci;
TotalMMs = pCurNode->NumMismatches;
TotalAlignLen = pCurNode->AlignLen;
QueryPathEndOfs = pCurNode->QueryStartOfs + pCurNode->AlignLen - 1;
TargPathEndOfs = pCurNode->TargSeqLoci + pCurNode->AlignLen - 1;
NumPathNodes = 1;
qNumInsert = 0;
qBaseInsert = 0;
tNumInsert = 0;
tBaseInsert = 0;

tsQueryAlignNodes *pTT = pCurNode;
uint32_t QueryPathEndOfsTT = QueryPathEndOfs;
uint32_t TargPathEndOfsTT = TargPathEndOfs;
int DeltaGap;
while ((CurPathNodeIdx = pTT->HiScorePathNextIdx) != 0)
	{
	pTT = &pAlignNodes[CurPathNodeIdx - 1];
	QueryGap = pTT->QueryStartOfs - QueryPathEndOfsTT - 1;
	TargGap = pTT->TargSeqLoci - TargPathEndOfsTT - 1;
	if (QueryGap < 0 || TargGap < 0)
		{
		if (QueryGap >= 0)	// if query gap is >= 0 then targ gap must be < 0
			DeltaGap = abs(TargGap);
		else
			{
			if (TargGap >= 0)	// if targ gap >= 0 then query gap must be < 0
				DeltaGap = abs(QueryGap);
			else    // else both query and targ gap must have been < 0
				DeltaGap = abs(min(QueryGap, TargGap));
			}
		pTT->QueryStartOfs += DeltaGap;
		pTT->TargSeqLoci += DeltaGap;
		pTT->AlignLen -= DeltaGap;
		}

	QueryPathEndOfsTT = pTT->QueryStartOfs + pTT->AlignLen - 1;
	TargPathEndOfsTT = pTT->TargSeqLoci + pTT->AlignLen - 1;
	}

while ((CurPathNodeIdx = pCurNode->HiScorePathNextIdx) != 0)
	{
	pCurNode = &pAlignNodes[CurPathNodeIdx - 1];
	QueryGap = pCurNode->QueryStartOfs - QueryPathEndOfs - 1;
	TargGap = pCurNode->TargSeqLoci - TargPathEndOfs - 1;
	TotalMMs += pCurNode->NumMismatches;
	TotalAlignLen += pCurNode->AlignLen;
	QueryPathEndOfs = pCurNode->QueryStartOfs + pCurNode->AlignLen - 1;
	TargPathEndOfs = pCurNode->TargSeqLoci + pCurNode->AlignLen - 1;
	if (QueryGap > 0)
		{
		qNumInsert += 1;
		qBaseInsert += QueryGap;
		}
	if (TargGap > 0)
		{
		tNumInsert += 1;
		tBaseInsert += TargGap;
		}
	NumPathNodes += 1;
	}
if(pTargPathEndOfs != NULL)
	*pTargPathEndOfs = TargPathEndOfs;
if(pQueryPathEndOfs != NULL)
	*pQueryPathEndOfs = QueryPathEndOfs;
if(pqNumInsert + NULL)
	*pqNumInsert = qNumInsert;
if (pqBaseInsert + NULL)
	*pqBaseInsert = qBaseInsert;
if (ptNumInsert + NULL)
	*ptNumInsert = tNumInsert;
if (ptBaseInsert + NULL)
	*ptBaseInsert = tBaseInsert;
return(NumPathNodes);
}

int											// returns number of reported paths
CBlitz::Report(uint32_t MinPathScore,			// only report paths having at least this minimum score
				uint32_t  MaxPathsToReport,		// report at most this many alignment paths for any query
				char *pszQuerySeqIdent,		// query sequence 
				uint32_t QueryLen,			// query length
				uint8_t *pQuerySeq,			// the query sequence
				uint32_t NumNodes,			// number of alignment nodes
				tsQueryAlignNodes *pAlignNodes,		// alignment nodes
				tsQueryAlignNodes **ppFirst2Rpts)	// allocated to hold ptrs to alignment nodes which are marked as being FlgFirst2tRpt
{
bool bStrand;
int PathHiScore;
uint32_t TargSeqLen;
uint32_t SortedPathIdx;
uint32_t NumHeadNodes;
uint32_t NumPathNodes;
uint32_t QueryPathStartOfs;
uint32_t QueryPathEndOfs;
uint32_t TargPathStartOfs;
uint32_t TargPathEndOfs;
uint32_t	qNumInsert;			// Number of inserts in query
uint32_t qBaseInsert;			// Number of bases inserted in query
uint32_t tNumInsert;			// Number of inserts in target
uint32_t tBaseInsert;			// Number of bases inserted in target
uint32_t PathsReported;
char szTargName[100];
tsQueryAlignNodes *pHeadNode;

if((NumHeadNodes = IdentifyHighScorePaths(MinPathScore, MaxPathsToReport,QueryLen,NumNodes,pAlignNodes,ppFirst2Rpts))==0)
	return(0);

PathsReported = 0;
for(SortedPathIdx=0; SortedPathIdx < min((int)NumHeadNodes,MaxPathsToReport); SortedPathIdx++)
	{
	pHeadNode = ppFirst2Rpts[SortedPathIdx];
	bStrand = pHeadNode->FlgStrand ? true : false;
	ConsolidateNodes(pHeadNode->FlgStrand,QueryLen, pQuerySeq, SortedPathIdx, pAlignNodes, ppFirst2Rpts);
	if((NumPathNodes = CharacterisePath(QueryLen,pQuerySeq,SortedPathIdx, pAlignNodes, ppFirst2Rpts, &QueryPathEndOfs, &TargPathEndOfs, &qNumInsert, &qBaseInsert, &tNumInsert, &tBaseInsert)) == 0)
		continue;
	PathsReported += 1;
	
	TargPathStartOfs = pHeadNode->TargSeqLoci;
	QueryPathStartOfs = pHeadNode->QueryStartOfs;
	PathHiScore = pHeadNode->HiScore;
	TargSeqLen = m_pSfxArray->GetSeqLen(pHeadNode->TargSeqID);
	m_pSfxArray->GetIdentName(pHeadNode->TargSeqID, sizeof(szTargName), szTargName);

	// recalc total matches, mismatches, indeterminates as the block lengths and starting offsets may have been modified
	uint32_t Matches,MisMatches,NumbNs;
	BlocksAlignStats(&Matches,&MisMatches,NULL,&NumbNs,bStrand == true ? '-' : '+',pQuerySeq,QueryLen, pHeadNode->TargSeqID,TargSeqLen,TargPathStartOfs,TargPathEndOfs,NumPathNodes,SortedPathIdx,pAlignNodes,ppFirst2Rpts);

	switch(m_RsltsFormat) {
		case eBLZRsltsSQLite:
				ReportAsSQLitePSL(Matches,MisMatches,0,NumbNs,				// TotalAlignLen -TotalMMs,TotalMMs,0,0,
									qNumInsert,								// qNumInsert, Number of inserts in query
									qBaseInsert,							// qBaseInsert, Number of bases inserted in query
									tNumInsert,								// tNumInsert, Number of inserts in target
									tBaseInsert,							// tBaseInsert, Number of bases inserted in target
									bStrand == true ? '-' : '+',
									pszQuerySeqIdent,QueryLen,QueryPathStartOfs,QueryPathEndOfs,
									szTargName,TargSeqLen,TargPathStartOfs,TargPathEndOfs,
									NumPathNodes,SortedPathIdx,pAlignNodes,ppFirst2Rpts);
			break;



		case eBLZRsltsPSL:
			ReportAsPSL(Matches,MisMatches,0,NumbNs,						// TotalAlignLen -TotalMMs,TotalMMs,0,0,
									qNumInsert,								// qNumInsert, Number of inserts in query
									qBaseInsert,							// qBaseInsert, Number of bases inserted in query
									tNumInsert,								// tNumInsert, Number of inserts in target
									tBaseInsert,							// tBaseInsert, Number of bases inserted in target
									bStrand == true ? '-' : '+',
									pszQuerySeqIdent,QueryLen,QueryPathStartOfs,QueryPathEndOfs,
									szTargName,TargSeqLen,TargPathStartOfs,TargPathEndOfs,
									NumPathNodes,SortedPathIdx,pAlignNodes,ppFirst2Rpts);
			break;

		case eBLZRsltsPSLX:
			ReportAsPSLX(Matches,MisMatches,0,NumbNs,						// TotalAlignLen -TotalMMs,TotalMMs,0,0,
									qNumInsert,								// qNumInsert, Number of inserts in query
									qBaseInsert,							// qBaseInsert, Number of bases inserted in query
									tNumInsert,								// tNumInsert, Number of inserts in target
									tBaseInsert,							// tBaseInsert, Number of bases inserted in target
									bStrand == true ? '-' : '+',
									pszQuerySeqIdent,QueryLen,pQuerySeq,QueryPathStartOfs,QueryPathEndOfs,
									pHeadNode->TargSeqID,
									szTargName,TargSeqLen,TargPathStartOfs,TargPathEndOfs,
									NumPathNodes,SortedPathIdx,pAlignNodes,ppFirst2Rpts);
			break;

		case eBLZRsltsMAF:
			PathHiScore = (int)(((double)PathHiScore / ((double)QueryLen * m_ExactMatchScore)) * 999.0);
			ReportAsMAF(PathHiScore,Matches,MisMatches,0,NumbNs,			// TotalAlignLen -TotalMMs,TotalMMs,0,0,
									qNumInsert,								// qNumInsert, Number of inserts in query
									qBaseInsert,							// qBaseInsert, Number of bases inserted in query
									tNumInsert,								// tNumInsert, Number of inserts in target
									tBaseInsert,							// tBaseInsert, Number of bases inserted in target
									bStrand == true ? '-' : '+',
									pszQuerySeqIdent,QueryLen,pQuerySeq,QueryPathStartOfs,QueryPathEndOfs,
									pHeadNode->TargSeqID,
									szTargName,TargSeqLen,TargPathStartOfs,TargPathEndOfs,
									NumPathNodes,SortedPathIdx,pAlignNodes,ppFirst2Rpts);
			break;

		case eBLZRsltsBED:
			PathHiScore = (int)(((double)PathHiScore / ((double)QueryLen * m_ExactMatchScore)) * 999.0);
			ReportAsBED(pszQuerySeqIdent,bStrand == true ? '-' : '+',PathHiScore,szTargName,NumPathNodes,TargPathStartOfs,TargPathEndOfs,SortedPathIdx,pAlignNodes,ppFirst2Rpts);
			break;
		}
	}
return(PathsReported);
}


#ifdef _WIN32
unsigned __stdcall LoadQuerySeqsFileThread(void * pThreadPars)
#else
void *LoadQuerySeqsFileThread(void * pThreadPars)
#endif
{
int Rslt;
tsLoadQuerySeqsThreadPars *pPars = (tsLoadQuerySeqsThreadPars *)pThreadPars;			// makes it easier not having to deal with casts!
CBlitz *pBlitzer = (CBlitz *)pPars->pThis;

if(pPars->bIsSAMOutput)
	Rslt = pBlitzer->ProcLoadSAMQuerySeqsFile(pPars);
else
	Rslt = pBlitzer->ProcLoadQuerySeqsFile(pPars);
pPars->Rslt = Rslt;
*pPars->pRslt = Rslt;
#ifdef _WIN32
_endthreadex(0);
return(eBSFSuccess);
#else
pthread_exit(NULL);
#endif
}

int
CBlitz::InitLoadQuerySeqs(void)
{
static tsLoadQuerySeqsThreadPars ThreadPars;

// initiate loading the reads
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading query sequences from file...");
m_ThreadLoadQuerySeqsRslt = -1;

ThreadPars.pRslt = &m_ThreadLoadQuerySeqsRslt;
ThreadPars.pThis = this;
ThreadPars.ThreadIdx = 0;
ThreadPars.Rslt = 0;
ThreadPars.bIsSAMOutput = m_RsltsFormat == eBLZRsltsSAM ? true : false;
ThreadPars.bIsSAMPE = (m_pszInputFilePE2 == NULL || m_pszInputFilePE2[0] == '\0') ? false : true;

#ifdef _WIN32
m_hThreadLoadQuerySeqs = ThreadPars.threadHandle = (HANDLE)_beginthreadex(NULL,0x0fffff,LoadQuerySeqsFileThread,&ThreadPars,0,&m_ThreadLoadQuerySeqsID);
#else
int ThreadRslt = ThreadPars.threadRslt = pthread_create (&m_ThreadLoadQuerySeqsID , NULL , LoadQuerySeqsFileThread , &ThreadPars );
#endif

// wait a few seconds, if major problems with loading reads then should show very quickly
#ifdef _WIN32
if(WAIT_TIMEOUT != WaitForSingleObject(m_hThreadLoadQuerySeqs, 3000))
	{
	CloseHandle(m_hThreadLoadQuerySeqs);
	m_hThreadLoadQuerySeqs = NULL;
	return(m_ThreadLoadQuerySeqsRslt);
	}
#else
struct timespec ts;
int JoinRlt;
clock_gettime(CLOCK_REALTIME, &ts);
ts.tv_sec += 3;
if((JoinRlt = pthread_timedjoin_np(m_ThreadLoadQuerySeqsID, NULL, &ts)) == 0)
	{
	m_ThreadLoadQuerySeqsID = 0;
	return(m_ThreadLoadQuerySeqsRslt);
	}
#endif
return(eBSFSuccess);
}

int												// returned enqueued query identifier
CBlitz::EnqueueQuerySeq(char *pszQueryIdent,    // query identifier parsed from fasta descriptor
			int QuerySeqLen,					// query sequence length
			uint8_t *pQuerySeq,					// query sequence
			char *pszQueryIdentPE2,				// PE2 query identifier parsed from fasta descriptor
			int QuerySeqLenPE2,					// PE2 query sequence length
			uint8_t *pQuerySeqPE2)				// PE2 query sequence
{
int Idx;
int SeqID;
uint8_t *pSeq;
uint8_t *pSeqPE2;
tsQuerySeq *psQuery;
bool bTermBackgoundThreads;

if((pSeq = new uint8_t [QuerySeqLen]) == NULL)
	return(eBSFerrMem);
memcpy(pSeq,pQuerySeq,QuerySeqLen);

if(QuerySeqLenPE2 > 0 && pQuerySeqPE2 != NULL && pszQueryIdentPE2 != NULL && pszQueryIdentPE2[0] != '\0')
	{
	if ((pSeqPE2 = new uint8_t[QuerySeqLenPE2]) == NULL)
		{
		delete []pSeq;
		return(eBSFerrMem);
		}
	memcpy(pSeqPE2, pQuerySeqPE2, QuerySeqLenPE2);
	}
else
	pSeqPE2 = NULL;

// any room left in query sequence queue to accept this sequence(s)?
int BackoffMS = 1;	// backoff attempting to acquire lock progressively from 1 to 1024ms to give readers more opportunity to remove sequences from queue
while(1) {
	AcquireLock(true);
	if((m_NumQuerySeqs + 1) < m_AllocdQuerySeqs)  // + 1 so can push 2 reads as a minimum on to queue so can queue PE reads
		break;
	bTermBackgoundThreads = m_TermBackgoundThreads != 0 ? true : false;
	ReleaseLock(true);
	if(bTermBackgoundThreads)	// need to immediately self-terminate?
		{
		delete []pSeq;
		if(pSeqPE2 != NULL)
			delete []pSeqPE2;
		return(eBSFSuccess);
		}
#ifdef _WIN32
	Sleep((DWORD)BackoffMS);			// sleep for waiting for queue space to take additional sequences
#else
	usleep(BackoffMS * 1000);	// sleep for 5ms waiting for queue space to take additional sequences
#endif
	if(BackoffMS < 1024)
		BackoffMS *= 2;
	}
Idx = (m_NxtQuerySeqIdx + m_NumQuerySeqs) % m_AllocdQuerySeqs;
psQuery = &m_pQuerySeqs[Idx];
SeqID = ++m_TotSeqIDs;
psQuery->SeqID = SeqID;
psQuery->pQuerySeq = pSeq;
psQuery->QuerySeqLen = QuerySeqLen; 
strncpy(psQuery->szQueryIdent,pszQueryIdent,cMaxQuerySeqIdentLen);
psQuery->szQueryIdent[cMaxQuerySeqIdentLen] = '\0';
m_NumQuerySeqs += 1;

if(pSeqPE2 != NULL)
	{
	Idx = (m_NxtQuerySeqIdx + m_NumQuerySeqs) % m_AllocdQuerySeqs;
	psQuery = &m_pQuerySeqs[Idx];
	SeqID = ++m_TotSeqIDs;
	psQuery->SeqID = SeqID;
	psQuery->pQuerySeq = pSeqPE2;
	psQuery->QuerySeqLen = QuerySeqLenPE2;
	strncpy(psQuery->szQueryIdent, pszQueryIdentPE2, cMaxQuerySeqIdentLen);
	psQuery->szQueryIdent[cMaxQuerySeqIdentLen] = '\0';
	m_NumQuerySeqs += 1;
	}
ReleaseLock(true);
return(SeqID);
}

int										// number of sequences returned, 0 if none to be dequeued, < 0 if errors
CBlitz::DequeueQuerySeq(int MaxLenQueryIdent,			// maximum length query identifier
	int *pSeqID,					// returned sequence identifier
	char *pszQueryIdent,			// where to return query identifier
	int *pQuerySeqLen,				// where to return query sequence length
	uint8_t **pDequeuedQuerySeq,		// where to return ptr to dequeued query sequence, caller is responsible for deleting memory allocated to hold the returned sequence (delete *pDequeuedQuerySeq)
	int *pSeqIDPE2,			// if SAM PE: returned PE2 sequence identifier
	char *pszQueryIdentPE2,	// if SAM PE: where to return PE2 query identifier
	int *pQuerySeqLenPE2,	// if SAM PE: where to return PE2 query sequence length
	uint8_t **pDequeuedQuerySeqPE2)	// if SAM PE: where to return ptr to dequeued PE2 query sequence, caller is responsible for deleting memory allocated to hold the returned sequence (delete *pDequeuedQuerySeqPE2)
{
int Num2Dequeue;
bool bTermBackgoundThreads;
bool bAllQuerySeqsLoaded;
uint8_t *pSeq;
tsQuerySeq *psQuery;


*pSeqID = 0;
*pszQueryIdent = '\0';
*pQuerySeqLen = 0;
*pDequeuedQuerySeq = NULL;

if (pSeqIDPE2 == NULL || pszQueryIdentPE2 == NULL || pQuerySeqLenPE2 == NULL || pDequeuedQuerySeqPE2 == NULL)
	Num2Dequeue = 1;
else
	{
	*pSeqIDPE2 = 0;
	*pszQueryIdentPE2 = '\0';
	*pQuerySeqLenPE2 = 0;
	*pDequeuedQuerySeqPE2 = NULL;
	Num2Dequeue = 2;
	}

// any sequences available to be dequeued?
int BackoffMS = 1;	// backoff attempting to acquire lock progressively from 1 to 256ms to give writer more opportunity to add sequences to queue
while(1) {
	AcquireLock(true);
	if((bTermBackgoundThreads = m_TermBackgoundThreads != 0 ? true : false)==true)
		{
		ReleaseLock(true);
		return(eBSFSuccess);
		}
	if(m_NumQuerySeqs >= Num2Dequeue)
		break;
	bAllQuerySeqsLoaded = m_bAllQuerySeqsLoaded;
	ReleaseLock(true);
	if(bAllQuerySeqsLoaded)	
		return(eBSFSuccess);
#ifdef _WIN32
	Sleep((DWORD)BackoffMS);	// sleep waiting for at least Num2Dequeue sequences to be queued
#else
	usleep(BackoffMS *1000);	
#endif
	if(BackoffMS < 256)
		BackoffMS *= 2;
	}

psQuery = &m_pQuerySeqs[m_NxtQuerySeqIdx++];
if(m_NxtQuerySeqIdx == m_AllocdQuerySeqs)
	m_NxtQuerySeqIdx = 0;
*pSeqID = psQuery->SeqID;
pSeq = psQuery->pQuerySeq;
psQuery->pQuerySeq = NULL;
*pQuerySeqLen = psQuery->QuerySeqLen; 
strncpy(pszQueryIdent,psQuery->szQueryIdent,MaxLenQueryIdent);
pszQueryIdent[MaxLenQueryIdent-1] = '\0';
m_NumQuerySeqs -= 1;
if(m_RsltsFormat == eBLZRsltsSQLite)
	m_pSQLitePSL->AddAlignSummary(m_ExprID,pszQueryIdent,*pQuerySeqLen,NULL,0);
*pDequeuedQuerySeq = pSeq; 
if (Num2Dequeue == 2)
	{
	psQuery = &m_pQuerySeqs[m_NxtQuerySeqIdx++];
	if (m_NxtQuerySeqIdx == m_AllocdQuerySeqs)
		m_NxtQuerySeqIdx = 0;
	*pSeqIDPE2 = psQuery->SeqID;
	pSeq = psQuery->pQuerySeq;
	psQuery->pQuerySeq = NULL;
	*pQuerySeqLenPE2 = psQuery->QuerySeqLen;
	strncpy(pszQueryIdentPE2, psQuery->szQueryIdent, MaxLenQueryIdent);
	pszQueryIdentPE2[MaxLenQueryIdent - 1] = '\0';
	m_NumQuerySeqs -= 1;
	*pDequeuedQuerySeqPE2 = pSeq;
	}
ReleaseLock(true);
return(Num2Dequeue);
}


int 
CBlitz::ProcLoadSAMQuerySeqsFile(tsLoadQuerySeqsThreadPars *pPars)
{
return(LoadSAMRawReads((m_pszInputFilePE2 == NULL || m_pszInputFilePE2[0] == '\0') ? false : true, m_pszInputFile,m_pszInputFilePE2));
}

teBSFrsltCodes
CBlitz::LoadSAMRawReads(bool bIsPairReads,	// true if paired end processing - PE1 reads in pszPE1File and PE2 reads in pszPE2File
						char *pszPE1File,					// process PE1 reads from this file
						char *pszPE2File)					// optionally process PE2 reads from this file
{
teBSFrsltCodes Rslt;
int Idx;
uint8_t *pReadBuff;
int PE1NumDescrReads;
int PE1DescrLen;
uint8_t szPE1DescrBuff[cMaxDescrLen];
int PE1ReadLen;
uint8_t PE1ReadBuff[cSAMtruncSeqLen+1];
int PE1NumReadsAccepted;
int PE1NumUnderlength;
int PE1NumOverlength;

int PE2NumDescrReads;
int PE2DescrLen;
uint8_t szPE2DescrBuff[cMaxDescrLen];
int PE2ReadLen;
uint8_t PE2ReadBuff[cSAMtruncSeqLen+1];
int PE2NumReadsAccepted;
int PE2NumUnderlength;
int PE2NumOverlength;
int NxtToSample;

CFasta PE1Fasta;
CFasta PE2Fasta;

AcquireLock(true);
m_bAllQuerySeqsLoaded = false;
m_LoadQuerySeqsRslt = eBSFSuccess;		// presumed success, changed if any processing errors
ReleaseLock(true);

if ((Rslt = (teBSFrsltCodes)PE1Fasta.Open(pszPE1File, true)) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Load: Unable to open '%s' [%s] %s", pszPE1File, PE1Fasta.ErrText((teBSFrsltCodes)Rslt), PE1Fasta.GetErrMsg());
	AcquireLock(true);
	m_bAllQuerySeqsLoaded = true;
	m_LoadQuerySeqsRslt = (teBSFrsltCodes)Rslt;
	ReleaseLock(true);
	return(Rslt);
	}

if (bIsPairReads)
	{
	if ((Rslt = (teBSFrsltCodes)PE2Fasta.Open(pszPE2File, true)) != eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Load: Unable to open '%s' [%s] %s", pszPE2File, PE2Fasta.ErrText((teBSFrsltCodes)Rslt), PE2Fasta.GetErrMsg());
		PE1Fasta.Close();
		AcquireLock(true);
		m_bAllQuerySeqsLoaded = true;
		m_LoadQuerySeqsRslt = (teBSFrsltCodes)Rslt;
		ReleaseLock(true);
		return(Rslt);
		}
	}

PE1NumDescrReads = 0;
PE1NumReadsAccepted = 0;
PE1NumUnderlength = 0;
PE1NumOverlength = 0;
PE2NumDescrReads = 0;
PE2NumReadsAccepted = 0;
PE2NumUnderlength = 0;
PE2NumOverlength = 0;
NxtToSample = m_SampleNthRawRead;
while ((Rslt = (teBSFrsltCodes)(PE1ReadLen = PE1Fasta.ReadSequence(PE1ReadBuff, sizeof(PE1ReadBuff) - 1, true, false))) > eBSFSuccess)
	{
	if (m_TermBackgoundThreads != 0)	// need to immediately self-terminate?
		break;

	PE1NumDescrReads += 1;
	if (PE1ReadLen == eBSFFastaDescr)		// just read a descriptor line which is expected for multifasta or fastq
		{
		PE1DescrLen = PE1Fasta.ReadDescriptor((char *)szPE1DescrBuff, sizeof(szPE1DescrBuff) - 1);
		szPE1DescrBuff[cMaxDescrLen - 1] = '\0';
		for (Idx = 0; Idx < cMaxDescrIDLen - 1; Idx++)
			{
			if (szPE1DescrBuff[Idx] == '\0' || isspace(szPE1DescrBuff[Idx]))
				break;
			}
		szPE1DescrBuff[Idx] = '\0';
		PE1DescrLen = Idx;
		PE1ReadLen = PE1Fasta.ReadSequence(PE1ReadBuff, sizeof(PE1ReadBuff));
		if (PE1ReadLen < 0 || PE1ReadLen == eBSFFastaDescr)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Problem parsing sequence after %d reads parsed", PE1NumDescrReads);
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Last descriptor parsed: %s", szPE1DescrBuff);
			PE1Fasta.Close();
			if (bIsPairReads)
				PE2Fasta.Close();
			AcquireLock(true);
			m_bAllQuerySeqsLoaded = true;
			m_LoadQuerySeqsRslt = eBSFerrParse;
			ReleaseLock(true);
			return(eBSFerrParse);
			}

			// if paired end processing then also load PE2 read
		if (bIsPairReads)
			{
			if ((Rslt = (teBSFrsltCodes)(PE2ReadLen = PE2Fasta.ReadSequence(PE2ReadBuff, sizeof(PE2ReadBuff) - 1, true, false))) <= eBSFSuccess)
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Problem parsing sequence after %d reads parsed", PE2NumDescrReads);
				if (PE2NumDescrReads)
					gDiagnostics.DiagOut(eDLFatal, gszProcName, "Last descriptor parsed: %s", szPE2DescrBuff);
				PE1Fasta.Close();
				PE2Fasta.Close();
				AcquireLock(true);
				m_bAllQuerySeqsLoaded = true;
				m_LoadQuerySeqsRslt = eBSFerrParse;
				ReleaseLock(true);
				return(eBSFerrParse);
				}

			if (m_TermBackgoundThreads != 0)	// need to immediately self-terminate?
				break;

			PE2NumDescrReads += 1;
			if (PE2ReadLen == eBSFFastaDescr)		// just read a descriptor line which would be as expected for multifasta or fastq
				{
				PE2DescrLen = PE2Fasta.ReadDescriptor((char *)szPE2DescrBuff, sizeof(szPE2DescrBuff) - 1);
				szPE2DescrBuff[cMaxDescrLen - 1] = '\0';
				for (Idx = 0; Idx < cMaxDescrIDLen - 1; Idx++)
					{
					if (szPE2DescrBuff[Idx] == '\0' || isspace(szPE2DescrBuff[Idx]))
						break;
					}
				szPE2DescrBuff[Idx] = '\0';
				PE2DescrLen = Idx;
				PE2ReadLen = PE2Fasta.ReadSequence(PE2ReadBuff, sizeof(PE2ReadBuff));
				if (PE2ReadLen < 0 || PE2ReadLen == eBSFFastaDescr)
					{
					gDiagnostics.DiagOut(eDLFatal, gszProcName, "Problem parsing sequence after %d reads parsed", PE2NumDescrReads);
					gDiagnostics.DiagOut(eDLFatal, gszProcName, "Last descriptor parsed: %s", szPE2DescrBuff);
					PE1Fasta.Close();
					PE2Fasta.Close();
					AcquireLock(true);
					m_bAllQuerySeqsLoaded = true;
					m_LoadQuerySeqsRslt = eBSFerrParse;
					ReleaseLock(true);
					return(eBSFerrParse);
					}
				}
			else
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "raw sequence file '%s' processing error: %s ", pszPE2File,
						Rslt > 0 ? "over length read" : "not a multifasta short reads or fastq file");
				PE1Fasta.Close();
				if (bIsPairReads)
					PE2Fasta.Close();
				AcquireLock(true);
				m_bAllQuerySeqsLoaded = true;
				m_LoadQuerySeqsRslt = eBSFerrParse;
				ReleaseLock(true);
				return(eBSFerrParse);
				}

			}

		if (m_SampleNthRawRead > 1)
			{
			NxtToSample += 1;
			if (m_SampleNthRawRead > NxtToSample)
				continue;
			NxtToSample = 0;
			}

		// ensure sequence lengths are within acceptable range
		if (PE1ReadLen > cSAMtruncSeqLen)
			{
			PE1NumOverlength += 1;
			continue;
			}
		if (PE1ReadLen < m_CoreLen)
			{
			PE1NumUnderlength += 1;
			continue;
			}
		if(bIsPairReads)
			{
			if (PE2ReadLen > cSAMtruncSeqLen)
				{
				PE2NumOverlength += 1;
				continue;
				}
			if (PE2ReadLen < m_CoreLen)
				{
				PE2NumUnderlength += 1;
				continue;
				}
			}

		// remove any repeat masking flags
		pReadBuff = PE1ReadBuff;
		for (Idx = 0; Idx < PE1ReadLen; Idx++, pReadBuff++)
			*pReadBuff &= ~cRptMskFlg;

		if (bIsPairReads)
			{
			pReadBuff = PE2ReadBuff;
			for (Idx = 0; Idx < PE2ReadLen; Idx++, pReadBuff++)
				*pReadBuff &= ~cRptMskFlg;
			}

		if (!bIsPairReads)
			{
			if ((Rslt = (teBSFrsltCodes)EnqueueQuerySeq((char *)szPE1DescrBuff, PE1ReadLen, PE1ReadBuff)) < 0)
				break;
			PE1NumReadsAccepted += 1;
			}
		else
			{
			if((Rslt = (teBSFrsltCodes)EnqueueQuerySeq((char *)szPE1DescrBuff, PE1ReadLen, PE1ReadBuff, (char *)szPE2DescrBuff, PE2ReadLen, PE2ReadBuff)) < 0)
				break;
			PE1NumReadsAccepted += 1;
			PE2NumReadsAccepted += 1;
			}
		continue;
		}
	else
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "raw sequence file '%s' processing error: %s ", pszPE1File,
				Rslt > 0 ? "over length read" : "not a multifasta short reads or fastq file");
		PE1Fasta.Close();
		if (bIsPairReads)
			PE2Fasta.Close();
		AcquireLock(true);
		m_bAllQuerySeqsLoaded = true;
		m_LoadQuerySeqsRslt = eBSFerrParse;
		ReleaseLock(true);
		return(eBSFerrParse);
		}
	}
if (Rslt != eBSFSuccess)
	{
	if (m_TermBackgoundThreads == 0)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Errors processing file: %s ", pszPE1File);
		while (PE1Fasta.NumErrMsgs())
			gDiagnostics.DiagOut(eDLFatal, gszProcName, PE1Fasta.GetErrMsg());
		}
	PE1Fasta.Close();
	if (bIsPairReads)
		PE2Fasta.Close();
	AcquireLock(true);
	m_bAllQuerySeqsLoaded = true;
	m_LoadQuerySeqsRslt = Rslt;
	ReleaseLock(true);
	return(Rslt);
	}
PE1Fasta.Close();
if (bIsPairReads)
	PE2Fasta.Close();

if (m_TermBackgoundThreads != 0)	// need to immediately self-terminate?
	return(eBSErrSession);

gDiagnostics.DiagOut(eDLInfo, gszProcName, "LoadReads: Total of %1.9d reads parsed and loaded from %s", PE1NumReadsAccepted, pszPE1File);
if (PE1NumUnderlength > 0)
	gDiagnostics.DiagOut(eDLWarn, gszProcName, "Load: total of %d under length ( < %dbp ) sequences sloughed from file '%s'", PE1NumUnderlength, m_CoreLen, pszPE1File);
if (PE1NumOverlength > 0)
	gDiagnostics.DiagOut(eDLWarn, gszProcName, "Load: total of %d over length ( > %dbp ) sequences sloughed from file '%s'", PE1NumOverlength, cSAMtruncSeqLen, pszPE1File);

if (bIsPairReads)
	{
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "LoadReads: Total of %1.9d reads parsed and loaded from %s", PE2NumReadsAccepted, pszPE2File);
	if (PE2NumUnderlength > 0)
		gDiagnostics.DiagOut(eDLWarn, gszProcName, "Load: total of %d under length ( < %dbp ) sequences sloughed from file '%s'", PE2NumUnderlength, m_CoreLen, pszPE2File);
	if (PE2NumOverlength > 0)
		gDiagnostics.DiagOut(eDLWarn, gszProcName, "Load: total of %d over length sequences ( > %dbp ) sloughed from file '%s'", PE2NumOverlength, cSAMtruncSeqLen, pszPE2File);
	}
AcquireLock(true);
m_bAllQuerySeqsLoaded = true;
m_LoadQuerySeqsRslt = eBSFSuccess;
ReleaseLock(true);
return(eBSFSuccess);
return(eBSFSuccess);
}


int
CBlitz::ProcLoadQuerySeqsFile(tsLoadQuerySeqsThreadPars *pPars)
{
CFasta Fasta;
CFasta FastaPE2;
unsigned char *pSeqBuff;
unsigned char *pMskBase;
uint32_t MskIdx;
size_t BuffOfs;
size_t AllocdBuffSize;
size_t AvailBuffSize;
char szName[_MAX_PATH];
char szDescription[cMaxDescrLen+1];
uint32_t SeqLen;
int Descrlen;
bool bFirstEntry;
bool bEntryCreated;
bool bTruncSeq;
int Rslt;
int SeqID;
int NxtToSample;
int *pRslt = pPars->pRslt;
AcquireLock(true);
m_bAllQuerySeqsLoaded = false;
m_LoadQuerySeqsRslt = eBSFSuccess;		// presumed success, changed if any processing errors
ReleaseLock(true);

if((Rslt=Fasta.Open(m_pszInputFile,true))!=eBSFSuccess)
	{
	if(Rslt != eBSFerrNotFasta)
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcLoadQuerySeqsFile: Unable to open '%s' [%s] %s",m_pszInputFile,Fasta.ErrText((teBSFrsltCodes)Rslt),Fasta.GetErrMsg());
	AcquireLock(true);
	m_bAllQuerySeqsLoaded = true;
	m_LoadQuerySeqsRslt = (teBSFrsltCodes)Rslt;		
	ReleaseLock(true);
	return(Rslt);
	}

// note malloc is used as can then simply realloc to expand as may later be required
AllocdBuffSize = (size_t)cAllocQuerySeqLen;
if((pSeqBuff = (unsigned char *)malloc(AllocdBuffSize)) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcLoadQuerySeqsFile:- Unable to allocate memory (%u bytes) for sequence buffer",(uint32_t)cAllocQuerySeqLen);
	Fasta.Close();
	*pRslt = eBSFerrMem;
	AcquireLock(true);
	m_bAllQuerySeqsLoaded = true;
	m_LoadQuerySeqsRslt = (teBSFrsltCodes)Rslt;		
	ReleaseLock(true);
	return(eBSFerrMem);
	}
AvailBuffSize = cAllocQuerySeqLen;

bFirstEntry = true;
bEntryCreated = false;
bTruncSeq = false;
SeqID = 0;
BuffOfs = 0;
NxtToSample = m_SampleNthRawRead;
while((Rslt = SeqLen = Fasta.ReadSequence(&pSeqBuff[BuffOfs],(int)min(AvailBuffSize,(size_t)cMaxQuerySeqLen),true,false)) > eBSFSuccess)
	{
	if(m_TermBackgoundThreads != 0)	// requested to immediately self-terminate?
		{
		Rslt = eBSErrSession;
		break;
		}

	if(SeqLen == eBSFFastaDescr)		// just read a descriptor line
		{
		SeqID++;
		if(bEntryCreated)				// add any previous entry
			{
			if((Rslt=EnqueueQuerySeq(szName,(int)BuffOfs,pSeqBuff)) <= eBSFSuccess)
				break;
			}
		Descrlen = Fasta.ReadDescriptor(szDescription, cMaxDescrLen);
		// An assumption - will one day bite real hard - is that the
		// fasta descriptor line starts with some form of unique identifier.
		// Use this identifier as the entry name.
		if(sscanf(szDescription," %s[ ,]",szName)!=1)
			sprintf(szName,"%s.%d",m_pszInputFile,++SeqID);

		bFirstEntry = false;
		bEntryCreated = true;
		bTruncSeq = false;
		BuffOfs = 0;
		continue;
		}
	else
		if(bFirstEntry)	// if there was no descriptor then dummy up one...
			{
			SeqID++;
			sprintf(szName,"%s.%d",m_pszInputFile,SeqID);
			strcpy(szDescription,"No Description provided");
			bFirstEntry = false;
			bEntryCreated = true;
			bTruncSeq = false;
			}
	if(bTruncSeq)
		continue;

	if (m_SampleNthRawRead > 1)
		{
		NxtToSample += 1;
		if (m_SampleNthRawRead > NxtToSample)
			continue;
		NxtToSample = 0;
		}

	// remove any repeat masking flags
	pMskBase = &pSeqBuff[BuffOfs];
	for(MskIdx = 0; MskIdx < SeqLen; MskIdx++,pMskBase++)
		*pMskBase &= ~cRptMskFlg;

	BuffOfs += SeqLen;
	if (BuffOfs > m_MaxQuerySeqLen)
		{
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"ProcLoadQuerySeqsFile:- Truncating over length query sequence '%s' to %d",szName, m_MaxQuerySeqLen);
		BuffOfs = m_MaxQuerySeqLen;
		AvailBuffSize = AllocdBuffSize - m_MaxQuerySeqLen;
		bTruncSeq = true;
		continue;
		}
	AvailBuffSize -= SeqLen;

	if(AvailBuffSize < (size_t)(cAllocQuerySeqLen / 2))
		{
		size_t NewSize = (size_t)cAllocQuerySeqLen + AllocdBuffSize;
		unsigned char *pTmp;
		if((pTmp = (unsigned char *)realloc(pSeqBuff,NewSize))==NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcLoadQuerySeqsFile:- Unable to reallocate memory (%u bytes) for sequence buffer",(uint32_t)NewSize);
			Rslt = eBSFerrMem;
			break;
			}
		pSeqBuff = pTmp;
		AllocdBuffSize = NewSize;
		AvailBuffSize = AllocdBuffSize - BuffOfs;
		}
	}
if(Rslt < eBSFSuccess && Rslt != eBSErrSession)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcLoadQuerySeqsFile: Parsing errors");
	while(Fasta.NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,Fasta.GetErrMsg());
	}

if(Rslt >= eBSFSuccess && bEntryCreated && BuffOfs > 0)			// last entry
	Rslt=EnqueueQuerySeq(szName,(int)BuffOfs,pSeqBuff);
if(Rslt > eBSFSuccess)
	Rslt = eBSFSuccess;
if(pSeqBuff != NULL)
	free(pSeqBuff);
*pRslt = Rslt;
AcquireSerialise();
m_bAllQuerySeqsLoaded = true;
m_LoadQuerySeqsRslt = (teBSFrsltCodes)Rslt;
ReleaseSerialise();
return(Rslt);
}

// SortQueryAlignNodes
// Sort alignment nodes by TargSeqID.QueryID.FlgStrand.QueryStartOfs.TargStartOfs ascending
int
CBlitz::SortQueryAlignNodes(const void *arg1, const void *arg2)
{
tsQueryAlignNodes *pEl1 = (tsQueryAlignNodes *)arg1;
tsQueryAlignNodes *pEl2 = (tsQueryAlignNodes *)arg2;

if(pEl1->TargSeqID > pEl2->TargSeqID)
	return(1);
if(pEl1->TargSeqID < pEl2->TargSeqID)
	return(-1);
if(pEl1->QueryID > pEl2->QueryID)
	return(1);
if(pEl1->QueryID < pEl2->QueryID)
	return(-1);
if(pEl1->FlgStrand > pEl2->FlgStrand)
	return(1);
if(pEl1->FlgStrand < pEl2->FlgStrand)
	return(-1);
if(pEl1->QueryStartOfs > pEl2->QueryStartOfs)
	return(1);
if(pEl1->QueryStartOfs < pEl2->QueryStartOfs)
	return(-1);
if(pEl1->TargSeqLoci > pEl2->TargSeqLoci)
	return(1);
if(pEl1->TargSeqLoci < pEl2->TargSeqLoci)
	return(-1);
return(0);
}

// SortHighScoreDescend
// Sort alignment nodes which are the first in path by score descending
int
CBlitz::SortHighScoreDescend(const void *arg1, const void *arg2)
{
tsQueryAlignNodes *pEl1 = *(tsQueryAlignNodes **)arg1;
tsQueryAlignNodes *pEl2 = *(tsQueryAlignNodes **)arg2;
if(pEl1->HiScore < pEl2->HiScore)
	return(1);
if(pEl1->HiScore > pEl2->HiScore)
	return(-1);
if(pEl1->TargSeqID > pEl2->TargSeqID)
	return(1);
if(pEl1->TargSeqID < pEl2->TargSeqID)
	return(-1);
return(0);
}




