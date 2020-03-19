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

	QueryLenAlignedPct = querylendpct->count ? querylendpct->ival[0] : cDfltBlitzMinQueryLenAlignedPct;
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
	if(CoreLen != 0 && (CoreLen < cMinBlitzCoreLen) || CoreLen > cMaxBlitzCoreLen)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: core or seed length '-C%d' specified outside of range %d..%d\n",CoreLen,cMinBlitzCoreLen,cMaxBlitzCoreLen);
		exit(1);
		}

	CoreDelta = coredelta->count ?  coredelta->ival[0] : CoreDelta;
	if((CoreDelta != 0 && CoreDelta < cMinBlitzCoreDelta) || CoreDelta > CoreLen)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: core delta '-c%d' specified outside of range %d..%d\n",CoreDelta,cMinBlitzCoreDelta, cMaxBlitzCoreLen);
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

	MismatchScore = mismatchscore->count ?  mismatchscore->ival[0] : cDfltBlitzMismatchScore;
	if(MismatchScore < cMinBlitzMismatchScore || MismatchScore > cMaxBlitzMismatchScore)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: mismatch penalty '-s%d' specified outside of range %d..%d\n",MismatchScore,cMinBlitzMismatchScore,cMaxBlitzMismatchScore);
		exit(1);
		}
	ExactMatchScore = exactmatchscore->count ?  exactmatchscore->ival[0] : cDfltBlitzExactMatchScore;
	if(ExactMatchScore < cMinBlitzExactMatchScore || ExactMatchScore > cMaxBlitzExactMatchScore)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: exact match score '-s%d' specified outside of range %d..%d\n",ExactMatchScore,cMinBlitzExactMatchScore,cMaxBlitzExactMatchScore);
		exit(1);
		}
	GapOpenScore = gapopenscore->count ?  gapopenscore->ival[0] : cDfltBlitzGapOpenScore;
	if(GapOpenScore < cMinBlitzGapOpenScore || GapOpenScore > cMaxBlitzGapOpenScore)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: gap open penalty '-s%d' specified outside of range %d..%d\n",GapOpenScore,cMinBlitzGapOpenScore,cMaxBlitzGapOpenScore);
		exit(1);
		}

	MaxOccKMerDepth = maxocckmerdepth->count ?  maxocckmerdepth->ival[0] : 0;
	if(MaxOccKMerDepth != 0 && (MaxOccKMerDepth < cMinBlitzOccKMerDepth) || MaxOccKMerDepth > cMaxBlitzOccKMerDepth)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: maximum depth to explore over-occurring seed K-mers '-k%d' specified outside of range %d..%d\n",MaxOccKMerDepth,cMinBlitzOccKMerDepth,cMaxBlitzOccKMerDepth);
		exit(1);
		}

	MinPathScore = minpathscore->count ?  minpathscore->ival[0] : MinPathScore;
	if(MinPathScore != 0 && (MinPathScore < cMinBlitzPathScore) || MinPathScore > cMaxBlitzPathScore)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: minimum path score '-p%d' specified outside of range %d..%d\n",MinPathScore,cMinBlitzPathScore,cMaxBlitzPathScore);
		exit(1);
		}

	MaxPathsToReport = maxpathstoreport->count ?  maxpathstoreport->ival[0] : cDfltBlitzMaxPathsToReport;
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
	char* pszExprName,				// experiment name
	char* pszExprDescr,				// experiment description
	char* pszParams,				// string containing blitz parameters
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
	char* pszInputFile,				// name of input file containing query sequences (PE1 if PE processing when output format is SAM)
	char* pszInputFilePE2,			// name of input file containing PE2 query sequences (only applies if output format is SAM)
	char* pszSfxFile,				// target as suffix array
	char* pszOutFile,				// where to write alignments
	int NumThreads)					// number of worker threads to use
{
	int Rslt;
	CBlitz* pBlitzer; 

	if ((pBlitzer = new CBlitz) == NULL)
	{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Fatal: Unable to instantiate CAligner");
		return(eBSFerrObj);
	}
	Rslt = pBlitzer->Process((char *)gpszSubProcess->pszName,PMode, SampleNthRawRead, pszExprName, pszExprDescr, pszParams, KMerDist, Sensitivity, AlignStrand, MismatchScore, ExactMatchScore, GapOpenScore, CoreLen, CoreDelta, MaxInsertLen, MaxOccKMerDepth,
		MinPathScore, QueryLenAlignedPct, MaxPathsToReport, RsltsFormat, pszInputFile, pszInputFilePE2, pszSfxFile, pszOutFile, NumThreads);
	delete pBlitzer;
	return(Rslt);
}






