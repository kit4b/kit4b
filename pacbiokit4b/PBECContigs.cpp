/*
This toolkit is a source base clone of 'PacBioKanga' release 4.4.2 (https://github.com/csiro-crop-informatics/biokanga) and contains
significant source code changes enabling new functionality and resulting process parameterisation changes. These changes have resulted in
incompatibility with 'PacBioKanga'.

Because of the potential for confusion by users unaware of functionality and process parameterisation changes then the modified source base
and resultant compiled executables have been renamed to 'kit4pacbio' - PacBio K-mer Informed Toolkit for Bioinformatics.
The renaming will force users of the 'PacBioKanga' toolkit to examine scripting which is dependent on existing 'PacBioKanga'
parameterisations so as to make appropriate changes if wishing to utilise 'kit4b' parameterisations and functionality.

'kit4pacbio' is being released under the Opensource Software License Agreement (GPLv3)
'kit4pacbio' is Copyright (c) 2019, 2020, Dr Stuart Stephen
Please contact Dr Stuart Stephen < stuartjs@g3bio.com > if you have any questions regarding 'kit4b'.

Original 'BioKanga' copyright notice has been retained and is as follows.
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

#include "pacbiokit4b.h"

#include "../libkit4b/bgzf.h"
#include "SSW.h"
#include "pacbiocommon.h"
#include "SeqStore.h"
#include "AssembGraph.h"
#include "MAConsensus.h"
#include "PBECContigs.h"


int
ProcPacBioECContigs(etPBPMode PMode,	// processing mode
		int DeltaCoreOfs,			// offset by this many bp the core windows of coreSeqLen along the probe sequence when checking for overlaps
		int MaxSeedCoreDepth,		// only further process a seed core if there are no more than this number of matching cores in all targeted sequences
		int MinSeedCoreLen,			// use seed cores of this length when identifying putative overlapping scaffold sequences
		int MinNumSeedCores,        // require at least this many seed cores between overlapping scaffold sequences
		int SWMatchScore,			// score for matching bases (0..50)
		int SWMismatchPenalty,		// mismatch penalty (-50..0)
		int SWGapOpenPenalty,		// gap opening penalty (-50..0)
		int SWGapExtnPenalty,		// gap extension penalty (-50..0)
		int SWProgExtnPenaltyLen,	// progressive gap scoring then only apply gap extension score if gap at least this length (0..63) - use if aligning PacBio
		int MaxArtefactDev,			// classify overlaps as artefactual if sliding window of 1Kbp over any overlap deviates by more than this percentage from the overlap mean
		int MinContigLen,			// only accepting contigs of at least this length (defaults to 15Kbp)
		int MinHCSeqLen,			// only accepting hiconfidence reads of at least this length (defaults to 1Kbp)
		char *pszContigFile,		// input multifasta contig file
		char *pszHiConfFile,		// input hiconfidence file
	    char *pszErrCorFile,		// name of file into which write error corrected contigs
		int NumThreads);			// maximum number of worker threads to use


#ifdef _WIN32
int ProcECContigs(int argc, char* argv[])
{
// determine my process name
_splitpath(argv[0],NULL,NULL,gszProcName,NULL);
#else
int
ProcECContigs(int argc, char** argv)
{
// determine my process name
CUtility::splitpath((char *)argv[0],NULL,gszProcName);
#endif

int iFileLogLevel;			// level of file diagnostics
int iScreenLogLevel;		// level of file diagnostics
char szLogFile[_MAX_PATH];	// write diagnostics to this file
int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure

int PMode;					// processing mode
int MinSeedCoreLen;			// use seed cores of this length when identifying putative overlapping sequences
int MinNumSeedCores;        // require at least this many seed cores between overlapping sequences before attempting SW
int DeltaCoreOfs;			// offset by this many bp the core windows of coreSeqLen along the probe sequence when checking for overlaps
int MaxSeedCoreDepth;		// only further process a seed core if there are no more than this number of matching cores in all targeted sequences

int SWMatchScore;			// score for matching bases (0..50)
int SWMismatchPenalty;		// mismatch penalty (0..50)
int SWGapOpenPenalty;		// gap opening penalty (0..50)
int SWGapExtnPenalty;		// gap extension penalty (0..50)
int SWProgExtnPenaltyLen;	// progressive gap scoring then only apply gap extension score if gap at least this length (0..63) - use if aligning PacBio

int MaxArtefactDev;			// classify overlaps as artefactual if sliding window of 1000bp over any overlap deviates by more than this percentage from the overlap mean

int MinContigLen;			// only accepting assembled contigs of at least this length (defaults to 15Kbp)
int MinHCSeqLen;			// only accepting high confidence reads of at least this length (defaults to 1Kbp)

char szContigFile[_MAX_PATH];			// input multifasta contig file
char szHiConfFile[_MAX_PATH];			// input hiconfidence file

char szOutFile[_MAX_PATH];				// where to write error corrected contg sequences

int NumberOfProcessors;		// number of installed CPUs
int NumThreads;				// number of threads (0 defaults to number of CPUs)

char szSQLiteDatabase[_MAX_PATH];	// results summaries to this SQLite file
char szExperimentName[cMaxDatasetSpeciesChrom+1];			// experiment name
char szExperimentDescr[1000];		// describes experiment

//
struct arg_lit  *help    = arg_lit0("h","help",                 "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *pmode = arg_int0("m","pmode","<int>",				"processing mode - 0 error correct");
struct arg_int *mincontiglen = arg_int0("l","mincontiglen","<int>",	"minimum contig sequence lengths (default 15000, range 500 to 100000)");

struct arg_int *minhcseqlen = arg_int0("p","minhcseqlen","<int>",		"minimum individual high confidence sequence length (default 1000, range 100 to 100000)");

struct arg_int *minseedcorelen = arg_int0("c","seedcorelen","<int>",			"use seed cores of this length when identifying putative overlapping sequences (default 35, range 10 to 50)");
struct arg_int *minseedcores = arg_int0("C","minseedcores","<int>",				"require at least this many accepted seed cores between overlapping sequences to use SW (default 30, range 1 to 50)");

struct arg_int *deltacoreofs = arg_int0("d","deltacoreofs","<int>",				"offset cores (default 10, range 1 to 25)");
struct arg_int *maxcoredepth = arg_int0("D","maxcoredepth","<int>",				"explore cores of less than this maximum depth (default 10000, range 1000 to 20000)");

struct arg_int *matchscore = arg_int0("x","matchscore","<int>",					"SW score for matching bases (default 1, range 1 to 50)");
struct arg_int *mismatchpenalty = arg_int0("X","mismatchpenalty","<int>",		"SW mismatch penalty (default 10, range 1 to 50)");
struct arg_int *gapopenpenalty = arg_int0("y","gapopenpenalty","<int>",			"SW gap opening penalty (default 12, range 1 to 50)");
struct arg_int *gapextnpenalty = arg_int0("Y","gapextnpenalty","<int>",			"SW gap extension penalty (default 6, range 1 to 50)");
struct arg_int *progextnpenaltylen = arg_int0("z","progextnpenaltylen","<int>",	"SW gap extension penalty only applied for gaps of at least this number of bases (default 1, range 1 to 63)");

struct arg_int *maxartefactdev = arg_int0("a","artefactdev","<int>",			 "classify overlaps as artefactual if 1Kbp window score deviates by more than this percentage from complete overlap mean (0 to disable, range 1 to 25)");


struct arg_file *hiconffile = arg_file1("i","hiconfseqs","<file>",		"input file containing higher confidence reads or sequences to be used in error correcton of contigs");
struct arg_file *contigfile = arg_file1("I","contigs","<file>",		"input file containing assembled contig sequences to be error corrected");
struct arg_file *outfile = arg_file1("o","out","<file>",									"output error corrected PacBio reads to this file");
struct arg_int *threads = arg_int0("T","threads","<int>",						"number of processing threads 0..128 (defaults to 0 which sets threads to number of CPU cores)");

struct arg_file *summrslts = arg_file0("q","sumrslts","<file>",				"Output results summary to this SQLite3 database file");
struct arg_str *experimentname = arg_str0("w","experimentname","<str>",		"experiment name SQLite3 database file");
struct arg_str *experimentdescr = arg_str0("W","experimentdescr","<str>",	"experiment description SQLite3 database file");

struct arg_end *end = arg_end(200);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					pmode,minseedcorelen,minseedcores,deltacoreofs,maxcoredepth,
					matchscore,mismatchpenalty,gapopenpenalty,gapextnpenalty,progextnpenaltylen,
					mincontiglen,minhcseqlen,maxartefactdev,
					summrslts,contigfile,hiconffile,experimentname,experimentdescr,
					outfile,threads,
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
		printf("\nPlease report any issues regarding usage of %s to https://github.com/kit4b/kit4b/issues\n\n",gszProcName);
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
	else
		szExperimentName[0] = '\0';

	gExperimentID = 0;
	gProcessID = 0;
	gProcessingID = 0;
	szSQLiteDatabase[0] = '\0';
	szExperimentDescr[0] = '\0';

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

		if(strlen(szExperimentName) < 1)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no SQLite experiment name specified with '-w<str>' option");
			return(1);
			}
		if(experimentdescr->count)
			{
			strncpy(szExperimentDescr,experimentdescr->sval[0],sizeof(szExperimentDescr)-1);
			szExperimentDescr[sizeof(szExperimentDescr)-1] = '\0';
			CUtility::TrimQuotedWhitespcExtd(szExperimentDescr);
			}
		if(strlen(szExperimentDescr) < 1)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no SQLite experiment description specified with '-W<str>' option");
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
		{
		szSQLiteDatabase[0] = '\0';
		szExperimentDescr[0] = '\0';
		}

	PMode = (etPBPMode)(pmode->count ? pmode->ival[0] : (int)0);
	if(PMode < 0 || PMode > 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Processing mode '-m%d' must be 0",PMode);
		return(1);
		}

	MinSeedCoreLen = cDfltSeedCoreLen;
	MinNumSeedCores = cDfltNumSeedCores;
	DeltaCoreOfs = cDfltDeltaCoreOfs;
	MaxSeedCoreDepth = cDfltMaxSeedCoreDepth;
	SWMatchScore = cDfltSWMatchScore;
	SWMismatchPenalty = abs(cDfltSWMismatchPenalty);
	SWGapOpenPenalty = abs(cDfltSWGapOpenPenalty);
	SWGapExtnPenalty = abs(cDfltSWGapExtnPenalty);
	SWProgExtnPenaltyLen = cDfltSWProgExtnLen;
	MinContigLen = cDfltMinPBSeqLen;
	MinHCSeqLen = cDfltMinHCSeqLen;

	MinSeedCoreLen = minseedcorelen->count ? minseedcorelen->ival[0] : cDfltScaffSeedCoreLen;
	if(MinSeedCoreLen < cMinSeedCoreLen || MinSeedCoreLen > cMaxSeedCoreLen)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: seed core length '-c%d' must be in range %d..%d",MinSeedCoreLen,cMinSeedCoreLen,cMaxSeedCoreLen);
		return(1);
		}

	MinNumSeedCores = minseedcores->count ? minseedcores->ival[0] : 30;
	if(MinNumSeedCores < cMinNumSeedCores || MinNumSeedCores > cMaxNumSeedCores)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Minimum number of accepted seed cores '-C%d' must be in range %d..%d",MinNumSeedCores,cMinNumSeedCores,cMaxNumSeedCores);
		return(1);
		}

	DeltaCoreOfs = deltacoreofs->count ? deltacoreofs->ival[0] : 10;
	if(DeltaCoreOfs < 1 || DeltaCoreOfs > cMaxDeltaCoreOfs)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: offset seed cores '-d%d' must be in range 1..%d",DeltaCoreOfs,cMaxDeltaCoreOfs);
		return(1);
		}

	MaxSeedCoreDepth = maxcoredepth->count ? maxcoredepth->ival[0] : cDfltMaxSeedCoreDepth;
	if(MaxSeedCoreDepth < 1000 || MaxSeedCoreDepth > 25000)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: maximum depth to explore seed cores '-D%d' must be in range 1000..25000",MaxSeedCoreDepth);
		return(1);
		}

	SWMatchScore = matchscore->count ? matchscore->ival[0] : 1;
	if(SWMatchScore < 1 || SWMatchScore > cMaxAllowedSWScore)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: SW Match score '-x%d' must be in range 1..%d",SWMatchScore,cMaxAllowedSWScore);
		return(1);
		}
	SWMismatchPenalty = mismatchpenalty->count ? mismatchpenalty->ival[0] : 10;
	if(SWMismatchPenalty < 1 || SWMismatchPenalty > cMaxAllowedSWScore)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: SW Mismatch penalty '-X%d' must be in range 1..%d",SWMismatchPenalty,cMaxAllowedSWScore);
		return(1);
		}
	SWGapOpenPenalty = gapopenpenalty->count ? gapopenpenalty->ival[0] : 12;
	if(SWGapOpenPenalty < 1 || SWGapOpenPenalty > cMaxAllowedSWScore)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: SW Gap open penalty '-y%d' must be in range 1..%d",SWGapOpenPenalty,cMaxAllowedSWScore);
		return(1);
		}
	SWGapExtnPenalty = gapextnpenalty->count ? gapextnpenalty->ival[0] : 6;
	if(SWGapExtnPenalty < 1 || SWGapExtnPenalty > cMaxAllowedSWScore)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: SW Gap extension penalty '-Y%d' must be in range 1..%d",SWGapExtnPenalty,cMaxAllowedSWScore);
		return(1);
		}
	SWProgExtnPenaltyLen = progextnpenaltylen->count ? progextnpenaltylen->ival[0] : 1;
	if(SWProgExtnPenaltyLen < 1 || SWProgExtnPenaltyLen > 63)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: SW Apply gap extension progress penalty for gaps at least '-z%d' must be in range 1..%d",SWProgExtnPenaltyLen,63);
		return(1);
		}

	MaxArtefactDev = maxartefactdev->count ? maxartefactdev->ival[0] : cDfltScaffMaxArtefactDev;
	if(MaxArtefactDev <= 0 || MaxArtefactDev > cMaxMaxArtefactDev)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Max overlap artefactual deviation '-a%d' must be either 0 or in range %d..%d",MaxArtefactDev,cMinMaxArtefactDev,cMaxMaxArtefactDev);
		return(1);
		}
	if(MaxArtefactDev < 0)
		MaxArtefactDev = 0;

	MinHCSeqLen = minhcseqlen->count ? minhcseqlen->ival[0] : cDfltMinHCSeqLen;
	if(MinHCSeqLen < cMinHCSeqLen || MinHCSeqLen > cMaxMinPBSeqLen)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Minimum accepted high confidence length '-l%p' must be in range %d..%dbp",MinHCSeqLen, cMinHCSeqLen,cMaxMinPBSeqLen);
		return(1);
		}

	MinContigLen = mincontiglen->count ? mincontiglen->ival[0] : cDfltMinPBSeqLen;
	if(MinContigLen < cMinPBSeqLen || MinContigLen > cMaxMinPBSeqLen)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Minimum accepted error corrected sequence length '-l%d' must be in range %d..%dbp",MinContigLen,cMinPBSeqLen,cMaxMinPBSeqLen);
		return(1);
		}


	if(!hiconffile->count)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: No file containing higher confidence reads or sequences specified with '-i<filespec>' option)\n");
		exit(1);
		}

	strncpy(szHiConfFile,hiconffile->filename[0],sizeof(szHiConfFile));
	szOutFile[sizeof(szHiConfFile)-1] = '\0';
	CUtility::TrimQuotedWhitespcExtd(szHiConfFile);
	if(szHiConfFile[0] == '\0')
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no file containing higher confidence reads or sequences specified with '-i<filespec>' option)\n");
		exit(1);
		}

	if(!contigfile->count)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: No assembled contig sequence file specified with '-I<filespec>' option)\n");
		exit(1);
		}

	strncpy(szContigFile,contigfile->filename[0],sizeof(szContigFile));
	szOutFile[sizeof(szContigFile)-1] = '\0';
	CUtility::TrimQuotedWhitespcExtd(szContigFile);
	if(szContigFile[0] == '\0')
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no contig sequences file specified with '-I<filespec>' option)\n");
		exit(1);
		}

	if(!outfile->count)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: No output files specified with '-o<filespec>' option)\n");
		exit(1);
		}

	strncpy(szOutFile,outfile->filename[0],sizeof(szOutFile));
	szOutFile[sizeof(szOutFile)-1] = '\0';
	CUtility::TrimQuotedWhitespcExtd(szOutFile);
	if(szOutFile[0] == '\0')
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no output file specified with '-o<filespec>' option)\n");
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

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing parameters:");

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode: 'Error correct assembled contigs'");

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Use seed cores of this length when identifying putative overlapping sequences: %dbp",MinSeedCoreLen);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Require at least this many seed cores between overlapping sequences: %d",MinNumSeedCores);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Offset cores by this many bp: %d",DeltaCoreOfs);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Maximum seed core depth: %d",MaxSeedCoreDepth);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"SW score for matching bases: %d",SWMatchScore);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"SW mismatch penalty: %d",SWMismatchPenalty);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"SW gap opening penalty: %d",SWGapOpenPenalty);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"SW gap extension penalty: %d",SWGapExtnPenalty);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"SW gap extension penalty only applied for gaps of at least this size: %d",SWProgExtnPenaltyLen);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"classify overlaps as artefactual if sliding window of 1Kbp over any overlap deviates by more than this percentage: %d",MaxArtefactDev);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum contig sequence length: %dbp",MinContigLen);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum high confidence sequence length: %dbp",MinHCSeqLen);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input assembled contig sequences file: '%s'",szContigFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input high confidence sequences file: '%s'",szHiConfFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output error corrected sequences file: '%s'",szOutFile);

	if(szExperimentName[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"This processing reference: %s",szExperimentName);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"number of threads : %d",NumThreads);

	if(gExperimentID > 0)
		{
		int ParamID;
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szLogFile),"log",szLogFile);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(PMode),"pmode",&PMode);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MaxSeedCoreDepth),"maxcoredepth",&MaxSeedCoreDepth);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MinSeedCoreLen),"seedcorelen",&MinSeedCoreLen);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MinNumSeedCores),"minseedcores",&MinNumSeedCores);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(DeltaCoreOfs),"deltacoreofs",&DeltaCoreOfs);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(SWMatchScore),"matchscore",&SWMatchScore);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(SWMismatchPenalty),"mismatchpenalty",&SWMismatchPenalty);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(SWGapOpenPenalty),"gapopenpenalty",&SWGapOpenPenalty);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(SWGapExtnPenalty),"gapextnpenalty",&SWGapExtnPenalty);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MaxArtefactDev),"artefactdev",&MaxArtefactDev);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(SWProgExtnPenaltyLen),"progextnpenaltylen",&SWProgExtnPenaltyLen);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MinNumSeedCores),"minseedcores",&MinNumSeedCores);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MinContigLen),"mincontiglen",&MinContigLen);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MinHCSeqLen),"minhcseqlen",&MinHCSeqLen);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szContigFile),"contigfile",szContigFile);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szHiConfFile),"hiconffile",szHiConfFile);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szOutFile),"out",szOutFile);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(NumThreads),"threads",&NumThreads);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(NumberOfProcessors),"cpus",&NumberOfProcessors);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szSQLiteDatabase),"sumrslts",szSQLiteDatabase);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szExperimentName),"experimentname",szExperimentName);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szExperimentDescr),"experimentdescr",szExperimentDescr);
		}


#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = ProcPacBioECContigs((etPBPMode)PMode,DeltaCoreOfs,MaxSeedCoreDepth,MinSeedCoreLen,MinNumSeedCores,SWMatchScore,-1 * SWMismatchPenalty,-1 * SWGapOpenPenalty,-1 * SWGapExtnPenalty,SWProgExtnPenaltyLen,
								MaxArtefactDev,MinContigLen,MinHCSeqLen,szContigFile,szHiConfFile,szOutFile,NumThreads);
	Rslt = Rslt >=0 ? 0 : 1;
	if(gExperimentID > 0)
		{
		if(gProcessingID)
			gSQLiteSummaries.EndProcessing(gProcessingID,Rslt);
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
ProcPacBioECContigs(etPBPMode PMode,	// processing mode
		int DeltaCoreOfs,			// offset by this many bp the core windows of coreSeqLen along the probe sequence when checking for overlaps
		int MaxSeedCoreDepth,		// only further process a seed core if there are no more than this number of matching cores in all targeted sequences
		int MinSeedCoreLen,			// use seed cores of this length when identifying putative overlapping scaffold sequences
		int MinNumSeedCores,        // require at least this many seed cores between overlapping scaffold sequences
		int SWMatchScore,			// score for matching bases (0..50)
		int SWMismatchPenalty,		// mismatch penalty (-50..0)
		int SWGapOpenPenalty,		// gap opening penalty (-50..0)
		int SWGapExtnPenalty,		// gap extension penalty (-50..0)
		int SWProgExtnPenaltyLen,	// progressive gap scoring then only apply gap extension score if gap at least this length (0..63) - use if aligning PacBio
		int MaxArtefactDev,			// classify overlaps as artefactual if sliding window of 1Kbp over any overlap deviates by more than this percentage from the overlap mean
		int MinContigLen,			// only accepting contigs of at least this length (defaults to 15Kbp)
		int MinHCSeqLen,			// only accepting hiconfidence reads of at least this length (defaults to 1Kbp)
		char *pszContigFile,			// input multifasta contig file
		char *pszHiConfFile,			// input hiconfidence file
	    char *pszErrCorFile,		// name of file into which write error corrected contigs
		int NumThreads)			// maximum number of worker threads to use
{
int Rslt;
CPBECContigs *pPBContigs;

if((pPBContigs = new CPBECContigs)==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: Unable to instantiate CPBECContigs");
	return(eBSFerrObj);
	}

Rslt = pPBContigs->Process(PMode,DeltaCoreOfs,MaxSeedCoreDepth,MinSeedCoreLen,MinNumSeedCores,SWMatchScore,SWMismatchPenalty,SWGapOpenPenalty,SWGapExtnPenalty,SWProgExtnPenaltyLen,
								MaxArtefactDev,MinContigLen,MinHCSeqLen,pszContigFile,pszHiConfFile,pszErrCorFile,NumThreads);
delete pPBContigs;
return(Rslt);
}

CPBECContigs::CPBECContigs() // relies on base classes constructors
{
m_pSfxArray = NULL;
m_pSeqStore = NULL;
m_pPBScaffNodes = NULL;
m_pMapEntryID2NodeIDs = NULL;
m_pMAConsensus = NULL;
m_bMutexesCreated = false;
m_hErrCorFile = -1;
Init();
}

CPBECContigs::~CPBECContigs() // relies on base classes destructor's
{
Reset(false);
}


void
CPBECContigs::Init(void)
{
if(m_pSfxArray != NULL)
	{
	delete m_pSfxArray;
	m_pSfxArray = NULL;
	}

if(m_pSeqStore != NULL)
	{
	delete m_pSeqStore;
	m_pSeqStore = NULL;
	}

if(m_pPBScaffNodes != NULL)
	{
	delete m_pPBScaffNodes;
	m_pPBScaffNodes = NULL;
	}

if(m_pMAConsensus != NULL)
	{
	delete m_pMAConsensus;
	m_pMAConsensus = NULL;
	}

if(m_pMapEntryID2NodeIDs != NULL)
	{
	delete m_pMapEntryID2NodeIDs;
	m_pMapEntryID2NodeIDs = NULL;
	}

if(m_hErrCorFile != -1)
	{
#ifdef _WIN32
	_commit(m_hErrCorFile);
#else
	fsync(m_hErrCorFile);
#endif
	close(m_hErrCorFile);
	m_hErrCorFile = -1;
	}

m_NumPBScaffNodes = 0;
m_AllocdPBScaffNodes = 0;
m_NumHiConfSeqs = 0;
m_NumOverlapProcessed = 0;
m_ProvOverlapping = 0;
m_ProvContained = 0;
m_ProvArtefact = 0;
m_ProvSWchecked = 0;
m_MaxHiConfSeqLen = 0;

m_MaxTargSeqLen = 0;

m_PMode = ePBPMErrCorrect;
m_OverlapFloat = cDfltScaffMaxOverlapFloat;
m_MinContigLen = cDfltMinPBSeqLen;	
m_MinHCSeqLen = cDfltMinHCSeqLen;

m_DeltaCoreOfs = cDfltDeltaCoreOfs;
m_MaxSeedCoreDepth = cDfltMaxSeedCoreDepth;
m_MinSeedCoreLen = cDfltSeedCoreLen;
m_MinNumSeedCores = cDfltNumSeedCores;

m_SWMatchScore = cDfltSWMatchScore;
m_SWMismatchPenalty = cDfltSWMismatchPenalty;	
m_SWGapOpenPenalty = cDfltSWGapOpenPenalty;
m_SWGapExtnPenalty = cDfltSWGapExtnPenalty;
m_SWProgExtnPenaltyLen = cDfltSWProgExtnLen;	
m_MaxArtefactDev = cDfltScaffMaxArtefactDev;

memset(m_szContigFile,0,sizeof(m_szContigFile));
memset(m_szHiConfFile,0,sizeof(m_szHiConfFile));

m_NumThreads = 0;
if(m_bMutexesCreated)
	DeleteMutexes();
m_bMutexesCreated = false; 
}

void
CPBECContigs::Reset(bool bSync)			// if bSync true then fsync before closing output file handles
{

Init();
}


int
CPBECContigs::CreateMutexes(void)
{
if(m_bMutexesCreated)
	return(eBSFSuccess);

m_CASSerialise = 0;
m_CASLock = 0;

m_bMutexesCreated = true;
return(eBSFSuccess);
}

void
CPBECContigs::DeleteMutexes(void)
{
if(!m_bMutexesCreated)
	return;
m_bMutexesCreated = false;
}

void
CPBECContigs::AcquireCASSerialise(void)
{
int SpinCnt = 10;
int BackoffMS = 1;

#ifdef _WIN32
while(InterlockedCompareExchange(&m_CASSerialise,1,0)!=0)
	{
	if(SpinCnt -= 1)
		continue;
	CUtility::SleepMillisecs(BackoffMS);
	SpinCnt = 10;
	if(BackoffMS < 50)
		BackoffMS += 1;
	else
		BackoffMS = 1 + (rand() % 31);
	}
#else
while(__sync_val_compare_and_swap(&m_CASSerialise,0,1)!=0)
	{
	if(SpinCnt -= 1)
		continue;
	CUtility::SleepMillisecs(BackoffMS);
	SpinCnt = 10;
	if(BackoffMS < 50)
		BackoffMS += 1;
	else
		BackoffMS = 1 + (rand() % 31);
	}
#endif
}

void
CPBECContigs::ReleaseCASSerialise(void)
{
#ifdef _WIN32
InterlockedCompareExchange(&m_CASSerialise,0,1);
#else
__sync_val_compare_and_swap(&m_CASSerialise,1,0);
#endif
}


void
CPBECContigs::AcquireCASLock(void)
{
int SpinCnt = 10;
int BackoffMS = 1;

#ifdef _WIN32
while(InterlockedCompareExchange(&m_CASLock,1,0)!=0)
	{
	if(SpinCnt -= 1)
		continue;
	CUtility::SleepMillisecs(BackoffMS);
	SpinCnt = 10;
	if(BackoffMS < 50)
		BackoffMS += 1;
	else
		BackoffMS = 1 + (rand() % 31);
	}
#else
while(__sync_val_compare_and_swap(&m_CASLock,0,1)!=0)
	{
	if(SpinCnt -= 1)
		continue;
	CUtility::SleepMillisecs(BackoffMS);
	SpinCnt = 10;
	if(BackoffMS < 50)
		BackoffMS += 1;
	else
		BackoffMS = 1 + (rand() % 31);
	}
#endif
}

void
CPBECContigs::ReleaseCASLock(void)
{
#ifdef _WIN32
InterlockedCompareExchange(&m_CASLock,0,1);
#else
__sync_val_compare_and_swap(&m_CASLock,1,0);
#endif
}


// ProcessFastaFile
// Parse input fasta format file into either a suffix array or sequence store
int
CPBECContigs::ProcessFastaFile(int MinSeqLen,		// only accept for indexing sequences of at least this length
				char *pszFile,						// file containing sequences
				bool bSeqStore,						// false is to load into suffix array, true to load into sequence store
				int Flags)							// default is for flags = cFlgLCSeq used with PacBio read sequences
{
CFasta Fasta;
uint8_t *pSeqBuff;
uint8_t *pMskBase;
uint32_t MskIdx;
size_t BuffOfs;
size_t AllocdBuffSize;
size_t AvailBuffSize;
char szName[cBSFSourceSize];
char szDescription[cBSFDescriptionSize];
uint32_t SeqLen;
uint32_t MaxSeqLen;
int Descrlen;
bool bFirstEntry;
bool bEntryCreated;
int Rslt;
int SeqID;
int NumSeqsAccepted;
size_t TotAcceptedLen;
uint32_t NumSeqsUnderlength;

if((Rslt=Fasta.Open(pszFile,true))!=eBSFSuccess)
	{
	if(Rslt != eBSFerrNotFasta)
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile: Unable to open '%s' [%s] %s",pszFile,Fasta.ErrText((teBSFrsltCodes)Rslt),Fasta.GetErrMsg());
	return(Rslt);
	}

AllocdBuffSize = (size_t)cMaxAllocBuffChunk * 16;
// note malloc is used as can then simply realloc to expand as may later be required
if((pSeqBuff = (uint8_t *)malloc(AllocdBuffSize)) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile:- Unable to allocate memory (%u bytes) for sequence buffer",(uint32_t)AllocdBuffSize);
	Fasta.Close();
	return(eBSFerrMem);
	}
AvailBuffSize = AllocdBuffSize;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"ProcessFastaFile:- Adding %s..",pszFile);

bFirstEntry = true;
bEntryCreated = false;
SeqID = 0;
BuffOfs = 0;
NumSeqsUnderlength = 0;
NumSeqsAccepted = 0;
TotAcceptedLen = 0;
MaxSeqLen = 0;
while((Rslt = SeqLen = Fasta.ReadSequence(&pSeqBuff[BuffOfs],(int)min(AvailBuffSize,(size_t)cMaxAllocBuffChunk),true,false)) > eBSFSuccess)
	{
	if(SeqLen == eBSFFastaDescr)		// just read a descriptor line
		{
		SeqID++;
		if(bEntryCreated)				// add any previous entry
			{
			if(BuffOfs < (size_t)MinSeqLen)
				NumSeqsUnderlength += 1;
			else
				{
				if(bSeqStore)
					{
					if(m_pSeqStore->AddSeq(Flags,szName,(uint32_t)BuffOfs,pSeqBuff) == 0)
						Rslt = -1;
					else
						Rslt = eBSFSuccess;
					}
				else
					Rslt=m_pSfxArray->AddEntry(szName,pSeqBuff,(uint32_t)BuffOfs,Flags);

				if(Rslt < eBSFSuccess)
					{
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile - error %d %s",Rslt,bSeqStore ? "internal" : m_pSfxArray->GetErrMsg());
					break;
					}
				else
					{
					NumSeqsAccepted += 1;
					TotAcceptedLen += BuffOfs;
					if(MaxSeqLen < (uint32_t)BuffOfs)
						MaxSeqLen = (uint32_t)BuffOfs;
					}
				}
			Rslt = eBSFSuccess;
			}
		Descrlen = Fasta.ReadDescriptor(szDescription,cBSFDescriptionSize);
		// An assumption - will one day bite real hard - is that the
		// fasta descriptor line starts with some form of unique identifier.
		// Use this identifier as the entry name.
		if(sscanf(szDescription," %s[ ,]",szName)!=1)
			sprintf(szName,"%s.%d",pszFile,++SeqID);

		bFirstEntry = false;
		bEntryCreated = true;
		BuffOfs = 0;
		continue;
		}
	else
		if(bFirstEntry)	// if there was no descriptor then dummy up one...
			{
			SeqID++;
			sprintf(szName,"%s.%d",pszFile,SeqID);
			strcpy(szDescription,"No Description provided");
			bFirstEntry = false;
			bEntryCreated = true;
			}

	// remove any repeat masking flags so that sorts can actually sort
	// if run of more than 25 Ns and at least 5 Ns to end of buffer then randomly mutate
	// every 13th N
	//	e.g <25Ns>r<12Ns>r<12Ns> where r is a pseudorandom base
	pMskBase = &pSeqBuff[BuffOfs];
	int SeqNs = 0;
	for(MskIdx = 0; MskIdx < SeqLen; MskIdx++,pMskBase++)
		{
		*pMskBase &= ~cRptMskFlg;
		if(*pMskBase == eBaseN && (MskIdx+5) < SeqLen)
			{
			if(++SeqNs > 25 &&
				pMskBase[1] == eBaseN &&
				pMskBase[2] == eBaseN &&
				pMskBase[3] == eBaseN &&
				pMskBase[4] == eBaseN)
				{
				if(!(SeqNs % 13))	// mutate every 13th
					*pMskBase = rand() % 4;
				}
			}
		else
			SeqNs = 0;
		}

	BuffOfs += SeqLen;
	AvailBuffSize -= SeqLen;
	if(AvailBuffSize < (size_t)(cMaxAllocBuffChunk / 8))
		{
		size_t NewSize = (size_t)cMaxAllocBuffChunk + AllocdBuffSize;
		uint8_t *pTmp;
		if((pTmp = (uint8_t *)realloc(pSeqBuff,NewSize))==NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile:- Unable to reallocate memory (%u bytes) for sequence buffer",(uint32_t)NewSize);
			return(eBSFerrMem);
			}
		pSeqBuff = pTmp;
		AllocdBuffSize = NewSize;
		AvailBuffSize = AllocdBuffSize - BuffOfs;
		}
	}
if(Rslt < eBSFSuccess && Rslt != eBSErrSession)
	{
	while(Fasta.NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,Fasta.GetErrMsg());
	}


if(Rslt >= eBSFSuccess && bEntryCreated && BuffOfs > 0)			// close entry
	{
	if(BuffOfs < (size_t)MinSeqLen)
		{
		NumSeqsUnderlength += 1;
		Rslt = eBSFSuccess;
		}
	else
		{
		if(bSeqStore)
			{
			if(m_pSeqStore->AddSeq(Flags,szName,(uint32_t)BuffOfs,pSeqBuff) == 0)
				Rslt = -1;
			else
				Rslt = eBSFSuccess;
			}
		else
			Rslt=m_pSfxArray->AddEntry(szName,pSeqBuff,(uint32_t)BuffOfs,Flags);

		if(Rslt < eBSFSuccess)
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile - error %d %s",Rslt,bSeqStore ? "internal" : m_pSfxArray->GetErrMsg());
		else
			{
			Rslt = eBSFSuccess;
			NumSeqsAccepted += 1;
			TotAcceptedLen += BuffOfs;
			if(MaxSeqLen < (uint32_t)BuffOfs)
				MaxSeqLen = (uint32_t)BuffOfs;
			}
		}
	}
if(pSeqBuff != NULL)
	free(pSeqBuff);
if(bSeqStore)
	m_MaxHiConfSeqLen = MaxSeqLen;
gDiagnostics.DiagOut(eDLInfo,gszProcName,"ProcessFastaFile - %d parsed, %d accepted, %dbp mean length, %d sequences not accepted as length under %dbp ",
					SeqID,NumSeqsAccepted,NumSeqsAccepted == 0 ? 0 :(int)(TotAcceptedLen/NumSeqsAccepted),NumSeqsUnderlength,MinSeqLen);

return(Rslt);
}

int64_t
CPBECContigs::EstSumSeqLens(int NumTargFiles,char **pszTargFiles)		// guestimate total sequence length by simply summing the lengths of each file - likely to grossly over estimate
{
int64_t SumFileSizes;
int Idx;
int NumGlobs;

CSimpleGlob glob(SG_GLOB_FULLSORT);

	// determine crude estimate of sum of sequence lengths from the sum of file sizes
for(Idx = 0; Idx < NumTargFiles; Idx++)
	if (glob.Add(pszTargFiles[Idx]) < SG_SUCCESS)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to glob '%s",pszTargFiles[Idx]);
		Reset(false);
		return((int64_t)eBSFerrOpnFile);
   		}


SumFileSizes = 0;
for (NumGlobs = 0; NumGlobs < glob.FileCount(); NumGlobs += 1)
	{
#ifdef _WIN32
	struct _stat64 st;
	if(!_stat64(glob.File(NumGlobs),&st))
#else
	struct stat64 st;
	if(!stat64(glob.File(NumGlobs),&st))
#endif
		SumFileSizes += (int64_t)st.st_size;
	}
if(NumGlobs == 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to locate any input files");
	Reset(false);
	return(eBSFerrOpnFile);
	}

return(SumFileSizes);
}


int 
CPBECContigs::LoadSeqs(int MinSeqLen,int NumTargFiles,
							char **pszTargFiles,		// parse, and index sequences in these files into in memory suffix array; file expected to contain either fasta or fastq sequences
							bool bSeqStore,						// false is to load into suffix array, true to load into sequence store
							int Flags)			// which by default are low confidence PacBio read sequences
{
int Rslt;
int Idx;
int NumGlobs;
int64_t SumFileSizes;

CSimpleGlob glob(SG_GLOB_FULLSORT);

for(Idx = 0; Idx < NumTargFiles; Idx++)
	if (glob.Add(pszTargFiles[Idx]) < SG_SUCCESS)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to glob '%s",pszTargFiles[Idx]);
		m_pSfxArray->Close(false);
		Reset(false);
		return(eBSFerrOpnFile);
   		}


SumFileSizes = 0;
for (NumGlobs = 0; NumGlobs < glob.FileCount(); NumGlobs += 1)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Will process in this order: %d '%s", NumGlobs+1,glob.File(NumGlobs));
#ifdef _WIN32
	struct _stat64 st;
	if(!_stat64(glob.File(NumGlobs),&st))
#else
	struct stat64 st;
	if(!stat64(glob.File(NumGlobs),&st))
#endif
		SumFileSizes += (int64_t)st.st_size;
	}
if(NumGlobs == 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to locate any input files");
	m_pSfxArray->Close(false);
	Reset(false);
	return(eBSFerrOpnFile);
	}

Rslt = eBSFSuccess;

for (int n = 0; Rslt >= eBSFSuccess &&  n < glob.FileCount(); ++n)
	{
		// try opening as a fasta file
	Rslt = ProcessFastaFile(MinSeqLen,glob.File(n),bSeqStore,Flags);
	if(Rslt < eBSFSuccess)
		{
		m_pSfxArray->Close(false);
		Reset(false);
		return(Rslt);
		}
	}
if(bSeqStore)
	Rslt = m_pSeqStore->GenSeqDescrIdx();

return(Rslt);

}

int
CPBECContigs::Process(etPBPMode PMode,	// processing mode
		int DeltaCoreOfs,			// offset by this many bp the core windows of coreSeqLen along the probe sequence when checking for overlaps
		int MaxSeedCoreDepth,		// only further process a seed core if there are no more than this number of matching cores in all targeted sequences
		int MinSeedCoreLen,			// use seed cores of this length when identifying putative overlapping scaffold sequences
		int MinNumSeedCores,        // require at least this many seed cores between overlapping scaffold sequences
		int SWMatchScore,			// score for matching bases (0..50)
		int SWMismatchPenalty,		// mismatch penalty (-50..0)
		int SWGapOpenPenalty,		// gap opening penalty (-50..0)
		int SWGapExtnPenalty,		// gap extension penalty (-50..0)
		int SWProgExtnPenaltyLen,	// progressive gap scoring then only apply gap extension score if gap at least this length (0..63) - use if aligning PacBio
		int MaxArtefactDev,			// classify overlaps as artefactual if sliding window of 1Kbp over any overlap deviates by more than this percentage from the overlap mean
		int MinContigLen,			// only accepting contigs of at least this length (defaults to 10Kbp)
		int MinHCSeqLen,			// only accepting hiconfidence reads of at least this length (defaults to 1Kbp)
		char *pszContigFile,		// input multifasta contig file
		char *pszHiConfFile,		// input hiconfidence file
	    char *pszErrCorFile,		// name of file into which write error corrected contigs
		int NumThreads)				// maximum number of worker threads to use
{
int Rslt = eBSFSuccess;
uint32_t CurNodeID;
uint32_t MaxSeqLen;
tsPBECCScaffNode *pCurPBScaffNode;
uint32_t NumTargSeqs;
uint64_t TotTargSeqLen;
uint32_t MaxTargSeqLen;
uint32_t RefSeqID;
etSeqBase *pRefSeq;

Reset(false);

CreateMutexes();

m_SWMatchScore = SWMatchScore;
m_SWMismatchPenalty = SWMismatchPenalty;
m_SWGapOpenPenalty = SWGapOpenPenalty;
m_SWGapExtnPenalty = SWGapExtnPenalty;
m_SWProgExtnPenaltyLen = SWProgExtnPenaltyLen;
m_MaxArtefactDev = MaxArtefactDev;

m_PMode = PMode;
m_DeltaCoreOfs = DeltaCoreOfs;
m_MaxSeedCoreDepth = MaxSeedCoreDepth;
m_MinSeedCoreLen = MinSeedCoreLen;
m_MinNumSeedCores = MinNumSeedCores;

m_MinContigLen = MinContigLen;	
m_MinHCSeqLen = MinHCSeqLen;
m_OverlapFloat = cDfltScaffMaxOverlapFloat;
strncpy(m_szContigFile,pszContigFile,sizeof(m_szContigFile));
m_szContigFile[sizeof(m_szContigFile)-1] = '\0';
strncpy(m_szHiConfFile,pszHiConfFile,sizeof(m_szHiConfFile));
m_szHiConfFile[sizeof(m_szHiConfFile)-1] = '\0';
strncpy(m_szErrCorFile,pszErrCorFile,sizeof(m_szErrCorFile));
m_szErrCorFile[sizeof(m_szErrCorFile)-1] = '\0';

int64_t SumContigsSizes;

SumContigsSizes = EstSumSeqLens(1,&pszContigFile);	// guestimate sum of all sequences; will return -1 if errors
if(SumContigsSizes == -1)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Unable to access the contigs file");
	return(eBSFerrObj);
	}

int64_t HiConfFileSizes;
HiConfFileSizes = EstSumSeqLens(1,&pszHiConfFile);	// guestimate sum of all sequences; will return eBSFerrOpnFile if errors
if(HiConfFileSizes == (int64_t)eBSFerrOpnFile)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Unable to access the high confidence sequences file");
	return(eBSFerrObj);
	}

m_NumThreads = NumThreads;	
if(m_pSfxArray != NULL)
	delete m_pSfxArray;
if((m_pSfxArray = new CSfxArray) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadTargetSeqs: Unable to instantiate instance of CSfxArray");
	return(eBSFerrObj);
	}
m_pSfxArray->Reset(false);
m_pSfxArray->SetMaxQSortThreads(m_NumThreads);

Rslt=m_pSfxArray->Open(false,false);
if(Rslt !=eBSFSuccess)
	{
	while(m_pSfxArray->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pSfxArray->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"CreateBioseqSuffixFile, unable to create in-memory suffix array - error %d %s",Rslt,m_pSfxArray->GetErrMsg());
	Reset(false);
	return(Rslt);
	}

if((Rslt=m_pSfxArray->SetDescription((char *)"inmem")) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"CreateBioseqSuffixFile: Unable to set description 'inmem'");
	Reset(false);
	return(Rslt);
	}
if((Rslt=m_pSfxArray->SetTitle((char *)"inmem")) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"CreateBioseqSuffixFile: Unable to set title 'inmem'");
	Reset(false);
	return(Rslt);
	}

if((Rslt = m_pSfxArray->SetDatasetName((char *)"inmem")) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"CreateBioseqSuffixFile: Unable to set dataset name 'inmem'");
	Reset(false);
	return(Rslt);
	}
m_pSfxArray->SetInitalSfxAllocEls(SumContigsSizes);	// just a hint which is used for initial allocations by suffix processing

m_pSeqStore = new CSeqStore;

if((Rslt = LoadSeqs(m_MinContigLen,1,&pszContigFile,false,cFlgLCSeq)) < eBSFSuccess)
	{
	Reset(false);
	return(Rslt);
	}

	// get number of contig sequences accepted for error correction
if((int)(NumTargSeqs = m_pSfxArray->GetNumEntries()) < 1)
	{
	NumTargSeqs = 0;
	TotTargSeqLen = 0;
	}
else
	TotTargSeqLen = m_pSfxArray->GetTotSeqsLen();
if(NumTargSeqs > 0)
	MaxTargSeqLen = m_pSfxArray->GetMaxSeqLen();
else
	MaxTargSeqLen = 0;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loaded for error correcting %d assembled contig sequences totaling %zdbp",NumTargSeqs,TotTargSeqLen);

if(MaxTargSeqLen > cMaxRefSeqLen)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to continue, at least one contig sequence %u is longer than max %u allowed",MaxTargSeqLen,cMaxRefSeqLen);
	m_pSfxArray->Close(false);
	Reset(false);
	return(eBSFerrFastqSeq);
	}

if(NumTargSeqs < 1)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Need at least 1 contig sequences for error correction accepted from file");
	m_pSfxArray->Close(false);
	Reset(false);
	return(eBSFerrNoEntries);
	}

m_MaxTargSeqLen = MaxTargSeqLen;

int NumHCSeqs = 0;

if((Rslt = LoadSeqs(m_MinHCSeqLen,1,&pszHiConfFile,true,cFlgHCSeq)) < eBSFSuccess)
	{
	m_pSfxArray->Close(false);
	Reset(false);
	return(Rslt);
	}

	// get number of high confidence sequences accepted and the maximum length of any accepted sequence
NumHCSeqs = m_pSeqStore->GetNumSeqs();
if(NumHCSeqs < 1)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Need at least 1 high confidence sequence for errror correction accepted from file");
	m_pSfxArray->Close(false);
	Reset(false);
	return(eBSFerrNoEntries);
	}
m_NumHiConfSeqs = NumHCSeqs;
m_MaxHiConfSeqLen =  m_pSeqStore->GetMaxSeqLen();

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loaded and accepted for processing a total of %d high confidence sequences for error correcting contigs",m_NumHiConfSeqs);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"CreateBioseqSuffixFile: sorting suffix array...");
m_pSfxArray->Finalise();
gDiagnostics.DiagOut(eDLInfo,gszProcName,"CreateBioseqSuffixFile: sorting completed");

m_pSfxArray->SetMaxIter(m_MaxSeedCoreDepth);
if(MinSeedCoreLen <= cMaxKmerLen)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Initialising for over occurring ( > %d) K-mers of length %d",m_MaxSeedCoreDepth,MinSeedCoreLen);
	if((Rslt = m_pSfxArray->InitOverOccKMers((int)MinSeedCoreLen,m_MaxSeedCoreDepth))!=eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Failed to initialise for over occurring K-mers");
		Reset(false);
		return(Rslt);
		}
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Initialised for over occurring K-mers");
	}

if((m_pPBScaffNodes = new tsPBECCScaffNode [NumTargSeqs + 1]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate for %d nodes",NumTargSeqs);
	Reset(false);
	return(eBSFerrMem);
	}
memset(m_pPBScaffNodes,0,sizeof(tsPBECCScaffNode) * (NumTargSeqs+1));
m_AllocdPBScaffNodes = NumTargSeqs;
if((m_pMapEntryID2NodeIDs = new uint32_t [NumTargSeqs + 1]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate for %d mapping entries nodes",NumTargSeqs);
	Reset(false);
	return(eBSFerrMem);
	}
memset(m_pMapEntryID2NodeIDs,0,sizeof(uint32_t) * (NumTargSeqs+1));
m_NumPBScaffNodes = 0;

if((m_pMAConsensus = new CMAConsensus)==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate instance of CMAConsensus");
	return(eBSFerrObj);
	}
if((Rslt=m_pMAConsensus->Init(NumTargSeqs,TotTargSeqLen))!=eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"MAConsensus initialialisation for %u targeted sequences totaling %zubp failed",NumTargSeqs,TotTargSeqLen);
	return(eBSFerrObj);
	}


if((pRefSeq = new etSeqBase [MaxTargSeqLen+10])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate %u temp memory for processing target contigs",MaxTargSeqLen+10);
	return(eBSFerrObj);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Initialising for %d sequences",NumTargSeqs);	
MaxSeqLen = 0;
pCurPBScaffNode = m_pPBScaffNodes;
for(CurNodeID = 1; CurNodeID <= NumTargSeqs; CurNodeID++,pCurPBScaffNode++)
	{
	uint16_t SeqFlags;
	pCurPBScaffNode->SeqLen = m_pSfxArray->GetSeqLen(CurNodeID);
	if(pCurPBScaffNode->SeqLen > cMaxRefSeqLen)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadTargetSeqs: Exceeded allowed max length (%u) with target sequence length %u",cMaxRefSeqLen,pCurPBScaffNode->SeqLen);
		delete pRefSeq;
		return(eBSFerrObj);
		}

	pCurPBScaffNode->EntryID = CurNodeID;
	m_pSfxArray->GetSeq(CurNodeID,0,pRefSeq,pCurPBScaffNode->SeqLen);
	RefSeqID = m_pMAConsensus->AddRefSeq(pCurPBScaffNode->SeqLen,pRefSeq);
	if(RefSeqID == 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadTargetSeqs: Failed adding reference sequence of length %u to consensus class",pCurPBScaffNode->SeqLen);
		delete pRefSeq;
		return(eBSFerrObj);
		}

	pCurPBScaffNode->RefSeqID = RefSeqID;
	SeqFlags = m_pSfxArray->GetIdentFlags(CurNodeID);
	pCurPBScaffNode->flgHCseq = SeqFlags & cFlgHCSeq ? 1 : 0;
	pCurPBScaffNode->flgUnderlength = pCurPBScaffNode->SeqLen < m_MinContigLen ? 1 : 0;
	if(MaxSeqLen == 0 || pCurPBScaffNode->SeqLen > (uint32_t)MaxSeqLen)
		MaxSeqLen = pCurPBScaffNode->SeqLen;
	}
m_NumPBScaffNodes = NumTargSeqs;
delete pRefSeq;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Sorting %d contig sequences ...",NumTargSeqs);

if(m_NumPBScaffNodes > 1)
	{
	// sort scaffold nodes by sequence length descending
	m_mtqsort.SetMaxThreads(NumThreads);
	m_mtqsort.qsort(m_pPBScaffNodes,m_NumPBScaffNodes,sizeof(tsPBECCScaffNode),SortLenDescending);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Sorting %d sequences completed",NumTargSeqs);
pCurPBScaffNode = m_pPBScaffNodes;
for(CurNodeID = 1; CurNodeID <= NumTargSeqs; CurNodeID++,pCurPBScaffNode++)
	{
	pCurPBScaffNode->NodeID = CurNodeID;
	m_pMapEntryID2NodeIDs[pCurPBScaffNode->EntryID-1] = CurNodeID;
	}

if(m_szErrCorFile[0] != '\0')
	{
#ifdef _WIN32
	m_hErrCorFile = open(m_szErrCorFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE));
#else
	if((m_hErrCorFile = open(m_szErrCorFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE))!=-1)
		if(ftruncate(m_hErrCorFile,0)!=0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",m_szErrCorFile,strerror(errno));
			Reset(false);
			return(eBSFerrCreateFile);
			}
#endif
	if(m_hErrCorFile < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output file '%s'",m_szErrCorFile);
		Reset(false);
		return(eBSFerrCreateFile);
		}
	}
else
	m_hErrCorFile = -1;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Initiating error correction ...");
Rslt = InitiateECContigs(NumThreads);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Completed error correction");
m_pMAConsensus->GenMultialignConcensus();

if(m_hErrCorFile != -1)
	{
	int ConsensusLen;
	int AllocdConsensusSeqSize;
	etSeqBase *pConsensusSeq;
	int BaseIdx;
	int BasesLine;
	int LineBuffIdx;
	char szLineBuff[2048];
	char szContig[100];
	AllocdConsensusSeqSize = 0;
	LineBuffIdx = 0;
	pConsensusSeq = NULL;
	m_pMAConsensus->GenMultialignConcensus();
	pCurPBScaffNode = m_pPBScaffNodes;
	for(CurNodeID = 1; CurNodeID <= NumTargSeqs; CurNodeID++,pCurPBScaffNode++)
		{
		m_pSfxArray->GetIdentName(pCurPBScaffNode->EntryID,sizeof(szContig)-1,szContig);
		ConsensusLen = m_pMAConsensus->GetConsensus(pCurPBScaffNode->RefSeqID);
		if(pConsensusSeq == NULL || (AllocdConsensusSeqSize < (ConsensusLen + 10)))
			{
			if(pConsensusSeq != NULL)
				{
				delete pConsensusSeq;
				pConsensusSeq = NULL;
				}
			AllocdConsensusSeqSize = 10 + ((ConsensusLen * 4) / 3);
			pConsensusSeq = new uint8_t [AllocdConsensusSeqSize];
			}
		ConsensusLen = m_pMAConsensus->GetConsensus(pCurPBScaffNode->RefSeqID,1,ConsensusLen,pConsensusSeq);

		// asciify the consensus and report to m_hErrCorFile
		BaseIdx = 0;
		BasesLine = 0;
		LineBuffIdx += sprintf(&szLineBuff[LineBuffIdx],">%sec %d\n",szContig,ConsensusLen);
		while(BaseIdx < ConsensusLen)
			{
			BasesLine = ConsensusLen - BaseIdx > 80 ? 80 : ConsensusLen - BaseIdx;
			CSeqTrans::MapSeq2Ascii(&pConsensusSeq[BaseIdx],BasesLine,&szLineBuff[LineBuffIdx]);
			BaseIdx += BasesLine;
			LineBuffIdx += BasesLine;
			LineBuffIdx += sprintf(&szLineBuff[LineBuffIdx],"\n");
			if(LineBuffIdx + 100 > sizeof(szLineBuff))
				{
				CUtility::RetryWrites(m_hErrCorFile,szLineBuff,LineBuffIdx);
				LineBuffIdx = 0;
				}
			}

		}
	if(LineBuffIdx)
		CUtility::RetryWrites(m_hErrCorFile,szLineBuff,LineBuffIdx);	

	if(m_hErrCorFile != -1)
		{
#ifdef _WIN32
		_commit(m_hErrCorFile);
#else
		fsync(m_hErrCorFile);
#endif
		close(m_hErrCorFile);
		m_hErrCorFile = -1;
		}
	}


Reset(false);
return(Rslt);
}

#ifdef _WIN32
unsigned __stdcall PBECContigsThread(void * pThreadPars)
#else
void *PBECContigsThread(void * pThreadPars)
#endif
{
int Rslt;
tsThreadPBECContigs *pPars = (tsThreadPBECContigs *)pThreadPars;			// makes it easier not having to deal with casts!
CPBECContigs *pPBContigs = (CPBECContigs *)pPars->pThis;

Rslt = pPBContigs->ThreadPBECContigs(pPars);
pPars->Rslt = Rslt;
#ifdef _WIN32
_endthreadex(0);
return(eBSFSuccess);
#else
pthread_exit(&pPars->Rslt);
#endif
}

int
CPBECContigs::InitiateECContigs(int NumECThreads)	// initiate contig error correction using this many threads
{
tsThreadPBECContigs *pThreadPutOvlps;
int ThreadIdx;
tsThreadPBECContigs *pThreadPar;

pThreadPutOvlps = new tsThreadPBECContigs [NumECThreads];

pThreadPar = pThreadPutOvlps;
for(ThreadIdx = 0; ThreadIdx < NumECThreads; ThreadIdx++,pThreadPar++)
	{
	memset(pThreadPar,0,sizeof(tsThreadPBECContigs));
	pThreadPar->AllocdCoreHits = cAllocdNumCoreHits;
	pThreadPar->AllocdCoreHitsSize = cAllocdNumCoreHits * sizeof(tsPBECCCoreHit);
#ifdef _WIN32
	pThreadPar->pCoreHits = (tsPBECCCoreHit *)malloc(pThreadPar->AllocdCoreHitsSize);	
	if(pThreadPar->pCoreHits == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"IdentifySequenceOverlaps: Core hits memory allocation of %zu bytes - %s",pThreadPar->AllocdCoreHitsSize,strerror(errno));
		break;
		}
#else
	if((pThreadPar->pCoreHits = (tsPBECCCoreHit *)mmap(NULL,pThreadPar->AllocdCoreHitsSize, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0)) == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"IdentifySequenceOverlaps: Core hits memory allocation of %zu bytes - %s",pThreadPar->AllocdCoreHitsSize,strerror(errno));
		break;
		}
#endif

	pThreadPar->AllocdProbeSeqSize = m_MaxHiConfSeqLen + 100;
	if((pThreadPar->pProbeSeq = new etSeqBase [pThreadPar->AllocdProbeSeqSize])==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"IdentifySequenceOverlaps: Core hits probe sequence memory allocation of %d bytes - %s",pThreadPar->AllocdProbeSeqSize,strerror(errno));
		break;
		}

	pThreadPar->AllocdTargSeqSize = m_MaxTargSeqLen + 100;
	if((pThreadPar->pTargSeq = new etSeqBase [pThreadPar->AllocdTargSeqSize])==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"IdentifySequenceOverlaps: Core hits target sequence memory allocation of %d bytes - %s",pThreadPar->AllocdTargSeqSize,strerror(errno));
		break;
		}

	if((pThreadPar->pmtqsort = new CMTqsort) == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"IdentifySequenceOverlaps: Core hits instantiation of CMTqsort failed");
		break;
		}
	pThreadPar->pmtqsort->SetMaxThreads(4);
	pThreadPar->bRevCpl = false;
	pThreadPar->MaxSeedCoreDepth = m_MaxSeedCoreDepth;
	pThreadPar->DeltaCoreOfs = m_DeltaCoreOfs;
	pThreadPar->CoreSeqLen = m_MinSeedCoreLen;
	pThreadPar->MinNumCores = m_MinNumSeedCores;
	pThreadPar->MaxAcceptHitsPerSeedCore = cDfltMaxAcceptHitsPerSeedCore;
	pThreadPar->MinPBSeqLen = m_MinContigLen;
	}

if(ThreadIdx != NumECThreads)	// any errors whilst allocating memory for core hits?
	{
	do {
		if(pThreadPar->pmtqsort != NULL)
			delete pThreadPar->pmtqsort;

		if(pThreadPar->pCoreHits != NULL)
			{
#ifdef _WIN32
			free(pThreadPar->pCoreHits);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
			if(pThreadPar->pCoreHits != MAP_FAILED)
				munmap(pThreadPar->pCoreHits,pThreadPar->AllocdCoreHitsSize);
#endif	
			}
		if(pThreadPar->pProbeSeq != NULL)
			delete pThreadPar->pProbeSeq;
		if(pThreadPar->pSW != NULL)
			delete pThreadPar->pSW;
		pThreadPar -= 1;
		ThreadIdx -= 1;
		}
	while(ThreadIdx >= 0);
	delete pThreadPutOvlps;
	Reset(false);
	return((int64_t)eBSFerrMem);
	}

pThreadPar = pThreadPutOvlps;
for (ThreadIdx = 1; ThreadIdx <= NumECThreads; ThreadIdx++, pThreadPar++)
	{
	pThreadPar->ThreadIdx = ThreadIdx;
	pThreadPar->pThis = this;
#ifdef _WIN32
	pThreadPar->threadHandle = (HANDLE)_beginthreadex(NULL, 0x0fffff, PBECContigsThread, pThreadPar, 0, &pThreadPar->threadID);
#else
	pThreadPar->threadRslt = pthread_create(&pThreadPar->threadID, NULL, PBECContigsThread, pThreadPar);
#endif
	}

// allow threads a few seconds to startup
#ifdef _WIN32
Sleep(5000);
#else
sleep(5);
#endif
pThreadPar = pThreadPutOvlps;
for (ThreadIdx = 0; ThreadIdx < NumECThreads; ThreadIdx++, pThreadPar++)
	{
#ifdef _WIN32
	while (WAIT_TIMEOUT == WaitForSingleObject(pThreadPar->threadHandle, 60000))
		{
		AcquireCASSerialise();
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: %u processed, SW aligned: %u, Overlapping: %u, Contained: %u, Artefact: %u",
							m_NumOverlapProcessed,m_ProvSWchecked,m_ProvOverlapping,m_ProvContained,m_ProvArtefact);
		ReleaseCASSerialise();
		};
	CloseHandle(pThreadPar->threadHandle);
#else
	struct timespec ts;
	int JoinRlt;
	clock_gettime(CLOCK_REALTIME, &ts);
	ts.tv_sec += 60;
	while ((JoinRlt = pthread_timedjoin_np(pThreadPar->threadID, NULL, &ts)) != 0)
		{
		AcquireCASSerialise();
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: %u processed, SW aligned: %u, Overlapping: %u, Contained: %u, Artefact: %u",
							m_NumOverlapProcessed,m_ProvSWchecked,m_ProvOverlapping,m_ProvContained,m_ProvArtefact);
		ReleaseCASSerialise();
		ts.tv_sec += 60;
		}
#endif
	}

pThreadPar = pThreadPutOvlps;
for(ThreadIdx = 0; ThreadIdx < NumECThreads; ThreadIdx++,pThreadPar++)
	{
	if(pThreadPar->pmtqsort != NULL)
		delete pThreadPar->pmtqsort;

	if(pThreadPar->pCoreHits != NULL)
		{
#ifdef _WIN32
		free(pThreadPar->pCoreHits);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
		if(pThreadPar->pCoreHits != MAP_FAILED)
			munmap(pThreadPar->pCoreHits,pThreadPar->AllocdCoreHitsSize);
#endif	
		pThreadPar->pCoreHits = NULL;
		}
	if(pThreadPar->pProbeSeq != NULL)
		delete pThreadPar->pProbeSeq;
	if(pThreadPar->pTargSeq != NULL)
		delete pThreadPar->pTargSeq;
	}

delete pThreadPutOvlps;
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Completed: %u processed, SW aligned: %u, Overlapping: %u, Contained: %u, Artefact: %u",
							m_NumOverlapProcessed,m_ProvSWchecked,m_ProvOverlapping,m_ProvContained,m_ProvArtefact);

return(0);
}

int
CPBECContigs::ThreadPBECContigs(tsThreadPBECContigs *pThreadPar)
{
uint32_t AdjOverlapFloat;
uint32_t CurTargCoreHitCnts;
tsPBECCScaffNode *pTargNode;
uint32_t ProbeAlignLength;
uint32_t TargAlignLength;
uint32_t TargSeqLen;
uint32_t LongSAligns;
uint32_t LongAAligns;
tsSSWCell *pPeakMatchesCell;
tsSSWCell PeakMatchesCell;
#ifdef _PEAKSCOREACCEPT
tsSSWCell PeakScoreCell;
#endif
bool bProbeSense;
uint32_t Idx;
uint32_t CurSummaryHitCnts;
uint32_t LowestSummaryHitCnts;
sPBECCCoreHitCnts *pLowestSummaryHitCnts;
uint32_t HitIdx;
tsPBECCCoreHit *pCoreHit;
uint32_t CurTargHitOfs;
uint32_t CurProbeHitOfs;
uint32_t CurSEntryIDHits;
uint32_t CurAEntryIDHits;
uint32_t	CurSTargStartOfs;
uint32_t	CurSTargEndOfs;
uint32_t	CurATargStartOfs;
uint32_t	CurATargEndOfs;
uint32_t	CurSProbeStartOfs;
uint32_t	CurSProbeEndOfs;
uint32_t	CurAProbeStartOfs;
uint32_t	CurAProbeEndOfs;
uint32_t ProvOverlapping;
uint32_t ProvOverlapped;
uint32_t ProvContained;
uint32_t ProvArtefact;
uint32_t ProvSWchecked;
uint32_t MinOverlapLen;
sPBECCCoreHitCnts *pSummaryCnts;
uint32_t HiConfSeqID;
uint32_t HiConfSeqFlags;
uint32_t HiConfSeqLen;
int NumInMultiAlignment;
int Class;
char szProbeSeqName[cMaxDatasetSpeciesChrom];

pThreadPar->pSW = NULL;
NumInMultiAlignment = 0;
AdjOverlapFloat = m_OverlapFloat + pThreadPar->CoreSeqLen;
if (pThreadPar->CoreSeqLen < cMaxPacBioSeedExtn)
	AdjOverlapFloat += 120;

for(HiConfSeqID = 1; HiConfSeqID <= m_NumHiConfSeqs; HiConfSeqID++)
	{
	pThreadPar->NumCoreHits = 0;
	AcquireCASLock();
	HiConfSeqFlags = m_pSeqStore->GetFlags(HiConfSeqID);
	if(HiConfSeqFlags & cflgAlgnd)    // another thread already aligning or completed aligning this sequence? 
       	{
		ReleaseCASLock();
		continue;
		}
	m_pSeqStore->SetFlags(HiConfSeqID,cflgAlgnd);
	ReleaseCASLock();
	HiConfSeqLen = m_pSeqStore->GetLen(HiConfSeqID);

	pThreadPar->bRevCpl = false;
	IdentifyCoreHits(HiConfSeqID,pThreadPar->MinPBSeqLen,0,pThreadPar);

	pThreadPar->bRevCpl = true;
	IdentifyCoreHits(HiConfSeqID,pThreadPar->MinPBSeqLen,0,pThreadPar);


	ProvOverlapping = 0;
	ProvOverlapped = 0;
	ProvContained = 0;
	ProvArtefact = 0;
	ProvSWchecked = 0;

	pThreadPar->NumTargCoreHitCnts = 0;
	memset(pThreadPar->TargCoreHitCnts,0,sizeof(pThreadPar->TargCoreHitCnts));

	if(pThreadPar->NumCoreHits >= pThreadPar->MinNumCores)
		{
			// resort core hits by TargNodeID.TargOfs.ProbeNodeID.ProbeOfs ascending
		pThreadPar->pmtqsort->qsort(pThreadPar->pCoreHits,pThreadPar->NumCoreHits,sizeof(tsPBECCCoreHit),SortCoreHitsByTargProbeOfs);
		pCoreHit = pThreadPar->pCoreHits;
		for(HitIdx = 0; HitIdx < pThreadPar->NumCoreHits; HitIdx++, pCoreHit++)
			pCoreHit->flgMulti = 0;

	// with large target sequences then can have many artifactual core hits
	// process and mark these probable artifact hits by identifying the most spatially related cluster of hits; hits outside of
	// the most spatially related cluster are marked as being artefactual 

		uint32_t CurTargSeqID;
		uint32_t ProbeSeqLen;
		uint32_t MaxWinSize;

		uint32_t RelWinSize;
		bool bFirstHitNewTargSeq;
		uint32_t PrevAcceptedProbeOfs;
		uint32_t PrevAcceptedTargOfs;
		uint32_t MaxNoHitGapLen;
		uint32_t NumNoHitGaps;
		uint32_t SumNoHitGapLens;
		tsPBECCCoreHit *pFirstCoreHit;
		tsPBECCCoreHit *pNxtCoreHit;
		tsPBECCCoreHit *pMaxCoreHit;

		MaxNoHitGapLen = 1000;                 // expecting cores to be distributed such that there shouldn't be many separated by more than this intercore bp gap
		CurTargSeqID = 0;
		pFirstCoreHit = NULL;
		ProbeSeqLen = HiConfSeqLen;
		MaxWinSize = (ProbeSeqLen * 115) / 100;
		pCoreHit = pThreadPar->pCoreHits;
		pMaxCoreHit = NULL;
		for (HitIdx = 0; HitIdx < pThreadPar->NumCoreHits; HitIdx++, pCoreHit++)
			{
			if (CurTargSeqID == 0)    // 0 if 1st hit about to be processed for a new target sequence
				{
				pMaxCoreHit = NULL;
				pFirstCoreHit = pCoreHit;
				CurTargSeqID = pCoreHit->TargNodeID;
				TargSeqLen = m_pPBScaffNodes[CurTargSeqID-1].SeqLen;
				}

		// if just checked last core hit for the current target ...
			if (HitIdx + 1 == pThreadPar->NumCoreHits || pCoreHit[1].TargNodeID != CurTargSeqID)
				{
				while (pFirstCoreHit <= pCoreHit)
					{
					pFirstCoreHit->WinHits = 0;
					pFirstCoreHit->flgClustered = 0;
					pNxtCoreHit = pFirstCoreHit;
					RelWinSize = 0;
					PrevAcceptedProbeOfs = 0;
					PrevAcceptedTargOfs = 0;
					SumNoHitGapLens = 0;
					NumNoHitGaps = 0;

					if(pFirstCoreHit->ProbeOfs > MaxNoHitGapLen && pFirstCoreHit->TargOfs > MaxNoHitGapLen)
						{
						pFirstCoreHit += 1;
						continue;
						}
					
					while (pNxtCoreHit <= pCoreHit)
						{
						RelWinSize = pNxtCoreHit->TargOfs - pFirstCoreHit->TargOfs;
						if (RelWinSize > MaxWinSize)
							break;
						if (pNxtCoreHit->flgRevCpl == pFirstCoreHit->flgRevCpl)	// only interested in hits which are same sense as the first hit in window
							{
							if(PrevAcceptedProbeOfs == 0 || pNxtCoreHit->ProbeOfs >= PrevAcceptedProbeOfs + pThreadPar->CoreSeqLen) // a single matched seed core extension may have resulted in multiple hits, reduce counts by requiring a differential of at least m_SeedCoreLen
								{
								if((pNxtCoreHit->ProbeOfs - PrevAcceptedProbeOfs) >= MaxNoHitGapLen)
									{
									NumNoHitGaps += 1;
									SumNoHitGapLens += pNxtCoreHit->ProbeOfs - PrevAcceptedProbeOfs;
									}
								PrevAcceptedProbeOfs = pNxtCoreHit->ProbeOfs;
								PrevAcceptedTargOfs = pNxtCoreHit->TargOfs;
								pFirstCoreHit->WinHits += 1;
								}
							}
						pNxtCoreHit += 1;
						}

					if(pFirstCoreHit->WinHits >= pThreadPar->MinNumCores && ((PrevAcceptedProbeOfs + MaxNoHitGapLen) > ProbeSeqLen || (PrevAcceptedTargOfs + MaxNoHitGapLen) > TargSeqLen))
						{
						uint32_t PutativeOverlapLen = PrevAcceptedProbeOfs - pFirstCoreHit->ProbeOfs;
						if(((SumNoHitGapLens * 100) / PutativeOverlapLen) <= 20)		// only accepting if no more than 20% of alignment sums to gaps > MaxNoHitGapLen
							{
							if (pMaxCoreHit == NULL || pFirstCoreHit->WinHits > pMaxCoreHit->WinHits)
								pMaxCoreHit = pFirstCoreHit;
							}
						}
					pFirstCoreHit += 1;
					}

				if (pMaxCoreHit != NULL)
				{
					RelWinSize = 0;
					pNxtCoreHit = pMaxCoreHit;
					while (pNxtCoreHit <= pCoreHit)
					{
						RelWinSize = pNxtCoreHit->TargOfs - pMaxCoreHit->TargOfs;
						if (RelWinSize > MaxWinSize)
							break;
						if (pNxtCoreHit->flgRevCpl == pMaxCoreHit->flgRevCpl)	// only interested in hits which are same sense as the first hit in window
							pNxtCoreHit->flgClustered = 1;
						pNxtCoreHit->flgMulti = 0;
						pNxtCoreHit += 1;
					}
				}
				CurTargSeqID = 0;  // looking for cores on a new target
			}
		}

	// iterate and count hits for each TargNodeID whilst recording the loci of the first and last hit so can determine if overlap is a likely artefact
	CurSEntryIDHits = 0;
	CurAEntryIDHits = 0;
	CurTargSeqID = 0;
	CurTargHitOfs = 0;
	CurProbeHitOfs = 0;
	CurSTargStartOfs = 0;
	CurSTargEndOfs = 0;
	CurATargStartOfs = 0;
	CurATargEndOfs = 0;
	CurSProbeStartOfs = 0;
	CurSProbeEndOfs = 0;
	CurAProbeStartOfs = 0;
	CurAProbeEndOfs = 0;

	bFirstHitNewTargSeq = false;
	pCoreHit = pThreadPar->pCoreHits;
	for(HitIdx = 0; HitIdx < pThreadPar->NumCoreHits; HitIdx++,pCoreHit++)
		{
		if(CurTargSeqID == 0)    // 0 if 1st hit about to be processed for a new target sequence
			{
			bFirstHitNewTargSeq = true;
			CurTargSeqID = pCoreHit->TargNodeID;
			CurSEntryIDHits = 0;
			CurAEntryIDHits = 0;
			CurSTargStartOfs = 0;
			CurSTargEndOfs = 0;
			CurATargStartOfs = 0;
			CurATargEndOfs = 0;
			CurSProbeStartOfs = 0;
			CurSProbeEndOfs = 0;
			CurAProbeStartOfs = 0;
			CurAProbeEndOfs = 0;
			}

		if(pCoreHit->flgClustered && pCoreHit->TargNodeID == CurTargSeqID) // same target sequence so check for starting/ending offsets and accumulate hit counts 
			{
			CurTargHitOfs = pCoreHit->TargOfs;
			CurProbeHitOfs = pCoreHit->ProbeOfs;
			if(pCoreHit->flgRevCpl == 0)
				{
				if(bFirstHitNewTargSeq == true || CurTargHitOfs < CurSTargStartOfs)
					CurSTargStartOfs = CurTargHitOfs;
				if(CurTargHitOfs > CurSTargEndOfs)
					CurSTargEndOfs = CurTargHitOfs;	
				if(bFirstHitNewTargSeq == true || CurProbeHitOfs < CurSProbeStartOfs)
					CurSProbeStartOfs = CurProbeHitOfs;
				if(CurProbeHitOfs > CurSProbeEndOfs)
					CurSProbeEndOfs = CurProbeHitOfs;
				}
			else
				{
				if(bFirstHitNewTargSeq == true || CurTargHitOfs < CurATargStartOfs)
					CurATargStartOfs = CurTargHitOfs;
				if(CurTargHitOfs > CurATargEndOfs)
					CurATargEndOfs = CurTargHitOfs;	
				if(bFirstHitNewTargSeq == true || CurProbeHitOfs < CurAProbeStartOfs)
					CurAProbeStartOfs = CurProbeHitOfs;
				if(CurProbeHitOfs > CurAProbeEndOfs)
					CurAProbeEndOfs = CurProbeHitOfs;
				}
			bFirstHitNewTargSeq = false;
			if(pCoreHit->flgMulti != 1)
				{
				if(pCoreHit->flgRevCpl == 0)
					CurSEntryIDHits += 1;
				else
					CurAEntryIDHits += 1;
				}
			}

		// if just processed last core hit for the current target ...
		if(HitIdx + 1 == pThreadPar->NumCoreHits || pCoreHit[1].TargNodeID != CurTargSeqID)
			{
			// checking here that the first and last hit are consistent with either a completely contained or overlapped
			TargSeqLen = m_pPBScaffNodes[pCoreHit->TargNodeID-1].SeqLen;
			if(CurSEntryIDHits >= pThreadPar->MinNumCores) 
				{
				if((CurSProbeStartOfs >= m_OverlapFloat &&  CurSTargStartOfs >= m_OverlapFloat) ||
						((TargSeqLen - CurSTargEndOfs) >= AdjOverlapFloat && (ProbeSeqLen - CurSProbeEndOfs) >= AdjOverlapFloat))
					CurSEntryIDHits = 0;
				}
			else
				CurSEntryIDHits = 0;

			if(CurAEntryIDHits >= pThreadPar->MinNumCores)
				{
				if((CurAProbeStartOfs >= m_OverlapFloat && CurATargStartOfs >= m_OverlapFloat) ||
					((TargSeqLen - CurATargEndOfs) >= AdjOverlapFloat && (ProbeSeqLen - CurAProbeEndOfs) >= AdjOverlapFloat))
					CurAEntryIDHits = 0;
				}
			else
				CurAEntryIDHits = 0;

			if(CurSEntryIDHits >= pThreadPar->MinNumCores || CurAEntryIDHits >= pThreadPar->MinNumCores)
				{
				if(pThreadPar->NumTargCoreHitCnts == cSummaryTargCoreHitCnts)
					{
					LowestSummaryHitCnts = 0;
					pLowestSummaryHitCnts = NULL;
					pSummaryCnts = pThreadPar->TargCoreHitCnts; 
					for(Idx = 0; Idx < cSummaryTargCoreHitCnts; Idx++, pSummaryCnts++)
						{
						CurSummaryHitCnts = pSummaryCnts->NumSHits + pSummaryCnts->NumAHits;
						if(LowestSummaryHitCnts == 0 || CurSummaryHitCnts < LowestSummaryHitCnts)
							{
							LowestSummaryHitCnts = CurSummaryHitCnts;
							pLowestSummaryHitCnts = pSummaryCnts;
							}
						}

					if((CurSEntryIDHits + CurAEntryIDHits) <= LowestSummaryHitCnts)
						{
						CurTargSeqID = 0; 
						continue;
						}
					pSummaryCnts = pLowestSummaryHitCnts;
					}
				else
					pSummaryCnts = &pThreadPar->TargCoreHitCnts[pThreadPar->NumTargCoreHitCnts++];
				pSummaryCnts->TargNodeID = CurTargSeqID;
				pSummaryCnts->STargStartOfs = CurSTargStartOfs;
				pSummaryCnts->STargEndOfs = CurSTargEndOfs;
				pSummaryCnts->ATargStartOfs = CurATargStartOfs;
				pSummaryCnts->ATargEndOfs = CurATargEndOfs;
				pSummaryCnts->SProbeStartOfs = CurSProbeStartOfs;
				pSummaryCnts->SProbeEndOfs = CurSProbeEndOfs;
				pSummaryCnts->AProbeStartOfs = CurAProbeStartOfs;
				pSummaryCnts->AProbeEndOfs = CurAProbeEndOfs;
				pSummaryCnts->NumSHits = CurSEntryIDHits;
				pSummaryCnts->NumAHits = CurAEntryIDHits;
				}
			CurTargSeqID = 0; 
			}
		}
	}


	if(pThreadPar->NumTargCoreHitCnts > 1)
		pThreadPar->pmtqsort->qsort(pThreadPar->TargCoreHitCnts,pThreadPar->NumTargCoreHitCnts,sizeof(sPBECCCoreHitCnts),SortCoreHitsDescending);

	NumInMultiAlignment = 0;
	if(pThreadPar->NumTargCoreHitCnts > 0)
		{
		LongSAligns = 0;
		LongAAligns = 0;
		if(pThreadPar->pSW == NULL)
			{
			AcquireCASSerialise();
			pThreadPar->pSW = new CSSW;
			pThreadPar->pSW->SetScores(m_SWMatchScore,m_SWMismatchPenalty,m_SWGapOpenPenalty,m_SWGapExtnPenalty,m_SWProgExtnPenaltyLen,min(63,m_SWProgExtnPenaltyLen+3),cAnchorLen);
			pThreadPar->pSW->PreAllocMaxTargLen(m_MaxTargSeqLen + 100,(m_MaxHiConfSeqLen * 11) / 10);
			ReleaseCASSerialise();
			}

		pSummaryCnts = &pThreadPar->TargCoreHitCnts[0];
		NumInMultiAlignment = 0;
		for(CurTargCoreHitCnts = 0; CurTargCoreHitCnts < pThreadPar->NumTargCoreHitCnts; CurTargCoreHitCnts++,pSummaryCnts++)
			{
			if(pSummaryCnts->NumSHits < pThreadPar->MinNumCores && pSummaryCnts->NumAHits < pThreadPar->MinNumCores)
				continue;
			bProbeSense = pSummaryCnts->NumSHits >= pSummaryCnts->NumAHits ? true :  false;
			pTargNode = &m_pPBScaffNodes[pSummaryCnts->TargNodeID-1];
			MinOverlapLen = (HiConfSeqLen * 7) / 10;
			TargSeqLen = pTargNode->SeqLen; 
			if(TargSeqLen + 10 > (uint32_t)pThreadPar->AllocdTargSeqSize)
				{
				delete pThreadPar->pTargSeq;
				pThreadPar->AllocdTargSeqSize = (TargSeqLen * 150) / 100;
				pThreadPar->pTargSeq = new etSeqBase [pThreadPar->AllocdTargSeqSize];
				}
			m_pSfxArray->GetSeq(pTargNode->EntryID,0,pThreadPar->pTargSeq,TargSeqLen);

			pThreadPar->pTargSeq[TargSeqLen] = eBaseEOS;
			pThreadPar->pSW->SetTarg(TargSeqLen,pThreadPar->pTargSeq);

			if(!bProbeSense)
				CSeqTrans::ReverseComplement(HiConfSeqLen,pThreadPar->pProbeSeq);
			m_pSeqStore->GetDescr(HiConfSeqID,sizeof(szProbeSeqName),szProbeSeqName);
			pThreadPar->pSW->SetProbe(HiConfSeqLen,pThreadPar->pProbeSeq);

			
			// restrict the range over which the SW will be processed to that of the overlap +/- m_OverlapFloat
			int Rslt;
			if(bProbeSense)
				{
				if(pSummaryCnts->SProbeStartOfs < m_OverlapFloat)
					pSummaryCnts->SProbeStartOfs = 0;
				else
					pSummaryCnts->SProbeStartOfs -= m_OverlapFloat;
				if(pSummaryCnts->SProbeEndOfs + AdjOverlapFloat >= HiConfSeqLen)
					pSummaryCnts->SProbeEndOfs = HiConfSeqLen - 1;
				else
					pSummaryCnts->SProbeEndOfs += AdjOverlapFloat;
				if(pSummaryCnts->STargStartOfs < m_OverlapFloat)
					pSummaryCnts->STargStartOfs = 0;
				else
					pSummaryCnts->STargStartOfs -= m_OverlapFloat;
				if(pSummaryCnts->STargEndOfs + AdjOverlapFloat >= TargSeqLen)
					pSummaryCnts->STargEndOfs = TargSeqLen - 1;
				else
					pSummaryCnts->STargEndOfs += AdjOverlapFloat;
				Rslt = pThreadPar->pSW->SetAlignRange(pSummaryCnts->SProbeStartOfs,pSummaryCnts->STargStartOfs,
											pSummaryCnts->SProbeEndOfs + 1 - pSummaryCnts->SProbeStartOfs,pSummaryCnts->STargEndOfs + 1 - pSummaryCnts->STargStartOfs);
				}
			else
				{

				if(pSummaryCnts->AProbeStartOfs < m_OverlapFloat)
					pSummaryCnts->AProbeStartOfs = 0;
				else
					pSummaryCnts->AProbeStartOfs -= m_OverlapFloat;
				if(pSummaryCnts->AProbeEndOfs + AdjOverlapFloat >= HiConfSeqLen)
					pSummaryCnts->AProbeEndOfs = HiConfSeqLen - 1;
				else
					pSummaryCnts->AProbeEndOfs += AdjOverlapFloat;
				if(pSummaryCnts->ATargStartOfs < m_OverlapFloat)
					pSummaryCnts->ATargStartOfs = 0;
				else
					pSummaryCnts->ATargStartOfs -= m_OverlapFloat;
				if(pSummaryCnts->ATargEndOfs + AdjOverlapFloat >= TargSeqLen)
					pSummaryCnts->ATargEndOfs = TargSeqLen - 1;
				else
					pSummaryCnts->ATargEndOfs += AdjOverlapFloat;

				Rslt = pThreadPar->pSW->SetAlignRange(pSummaryCnts->AProbeStartOfs,pSummaryCnts->ATargStartOfs,
											pSummaryCnts->AProbeEndOfs + 1 - pSummaryCnts->AProbeStartOfs,pSummaryCnts->ATargEndOfs + 1 - pSummaryCnts->ATargStartOfs);
				}

			pPeakMatchesCell = pThreadPar->pSW->Align(NULL,m_MaxHiConfSeqLen);
			ProvSWchecked += 1;
			if(pPeakMatchesCell != NULL && pPeakMatchesCell->NumMatches >= MinOverlapLen)
				{
				ProvOverlapping += 1;
				PeakMatchesCell = *pPeakMatchesCell;
				ProbeAlignLength = PeakMatchesCell.EndPOfs - PeakMatchesCell.StartPOfs + 1;
				TargAlignLength = PeakMatchesCell.EndTOfs - PeakMatchesCell.StartTOfs + 1;

				// characterise the overlapped target
				// eOLCOverlapping if probe accepted as overlapping, either 5' or 3'
				// eOLCcontaining if both ends of target completely contained within probe
                // eOLCartefact if target is only partially contained
				int PathClass;
				Class = (int)eOLCOverlapping;		
				if((PeakMatchesCell.StartTOfs >= m_OverlapFloat &&  PeakMatchesCell.StartPOfs >= m_OverlapFloat) ||
					 ((TargSeqLen - PeakMatchesCell.EndTOfs) >= m_OverlapFloat && (HiConfSeqLen - PeakMatchesCell.EndPOfs) >= m_OverlapFloat))
					{
					Class = (int)eOLCartefact;
					ProvArtefact += 1;
					}

				if(Class == eOLCOverlapping && (PathClass = pThreadPar->pSW->ClassifyPath(m_MaxArtefactDev,
																					PeakMatchesCell.PFirstAnchorStartOfs,PeakMatchesCell.PLastAnchorEndOfs,
																					PeakMatchesCell.TFirstAnchorStartOfs,PeakMatchesCell.TLastAnchorEndOfs)) > 0)
					{
					Class = (int)eOLCartefact;
					ProvArtefact += 1;
					}

				if(Class == eOLCOverlapping)
					{
					if(PeakMatchesCell.StartPOfs < m_OverlapFloat && (HiConfSeqLen - PeakMatchesCell.EndPOfs) < m_OverlapFloat) // is probe completely contained within target?
						Class = (int)eOLCcontained;
					else
						if(PeakMatchesCell.StartTOfs < m_OverlapFloat && (TargSeqLen - PeakMatchesCell.EndTOfs) < m_OverlapFloat) // is target completely contained by probe?
							Class = (int)eOLCcontains;
					if(Class != eOLCOverlapping)
						{
						AcquireCASLock();
						if(Class == (int)eOLCcontained)
							pTargNode->flgContains = 1;
						else
							pTargNode->flgContained = 1;
						ReleaseCASLock();	
						ProvContained += 1;
						}
					}
				
				if(Class == eOLCcontained)
					{
					tMAOp *pAlignOps;
					int NumAlignOps;
					NumAlignOps = pThreadPar->pSW->TracebacksToAlignOps(PeakMatchesCell.PFirstAnchorStartOfs,PeakMatchesCell.PLastAnchorEndOfs,
															PeakMatchesCell.TFirstAnchorStartOfs,PeakMatchesCell.TLastAnchorEndOfs,&pAlignOps);

					AcquireCASLock();
					m_pMAConsensus->AddMultiAlignment(pTargNode->RefSeqID,
															PeakMatchesCell.TFirstAnchorStartOfs,PeakMatchesCell.TLastAnchorEndOfs,
															PeakMatchesCell.PFirstAnchorStartOfs,PeakMatchesCell.PLastAnchorEndOfs,
															pThreadPar->pProbeSeq,NumAlignOps,pAlignOps);	
					ReleaseCASLock();								
					NumInMultiAlignment += 1;
					}

				ProvOverlapped += 1;
				}
			if(!bProbeSense)
				{
				CSeqTrans::ReverseComplement(HiConfSeqLen,pThreadPar->pProbeSeq);
				bProbeSense = true;
				}
			}
		}

	AcquireCASSerialise();
	if(ProvOverlapping > 0)
		m_ProvOverlapping += 1;
	if(ProvContained > 0)
		m_ProvContained += ProvContained;
	if(ProvArtefact > 0)
		m_ProvArtefact += ProvArtefact;
	if(ProvSWchecked > 0)
		m_ProvSWchecked += ProvSWchecked;
	m_NumOverlapProcessed += 1;
	ReleaseCASSerialise();
	ProvOverlapping = 0;
	ProvOverlapped = 0;
	ProvContained = 0;
	ProvArtefact = 0;
	ProvSWchecked = 0;
	NumInMultiAlignment = 0;
	pThreadPar->NumTargCoreHitCnts = 0;
	pThreadPar->NumCoreHits = 0;
	}

if(pThreadPar->pSW != NULL)
	{
	delete pThreadPar->pSW;
	pThreadPar->pSW = NULL;
	}
return(0);
}


int					// returns 0 if core overlapped (uses a non-exhaustive search) a previously added core, index 1..N of just added core hit or -1 if errors
CPBECContigs::AddCoreHit(uint32_t ProbeNodeID,			// core hit was from this probe scaffold node 
			   bool bRevCpl,					// true if core sequence was revcpl'd before matching
			   uint32_t ProbeOfs,                 // hit started at this probe offset
			   uint32_t TargNodeID,               // probe core matched onto this target scaffold node
			   uint32_t TargOfs,                  // probe core matched starting at this target loci
			   uint32_t HitLen,					// hit was of this length
               tsThreadPBECContigs *pPars)			// thread specific
{
uint32_t ExpdHitLen;
uint32_t NumHits2Chk; 
uint32_t HitsChkd;
tsPBECCCoreHit *pCurCoreHit;

if((pPars->NumCoreHits + 5) > pPars->AllocdCoreHits)	// need to realloc memory to hold additional cores?
	{
		// realloc memory with a 25% increase over previous allocation 
	int coresreq;
	size_t memreq;
	void *pAllocd;
	coresreq = (int)(((int64_t)pPars->AllocdCoreHits * 125) / (int64_t)100);
	memreq = coresreq * sizeof(tsPBECCCoreHit);

#ifdef _WIN32
		pAllocd = realloc(pPars->pCoreHits,memreq);
#else
		pAllocd = mremap(pPars->pCoreHits,pPars->AllocdCoreHitsSize,memreq,MREMAP_MAYMOVE);
		if(pAllocd == MAP_FAILED)
			pAllocd = NULL;
#endif
		if(pAllocd == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"SavePartialSeqs: Memory re-allocation to %d bytes - %s",memreq,strerror(errno));
			Reset(false);
			return(eBSFerrMem);
			}

		pPars->pCoreHits = (tsPBECCCoreHit *)pAllocd;
		pPars->AllocdCoreHitsSize = memreq;
		pPars->AllocdCoreHits = coresreq; 
		}
		

// non-exhaustive check to see if core is overlapping a previously added core
// search back over at most 2000 recently added cores for overlaps
ExpdHitLen = min(50,HitLen * 3); // allowing for some float on the overlap loci and length
if(pPars->NumCoreHits > 0)
	{
	NumHits2Chk = min(2000,pPars->NumCoreHits);
	pCurCoreHit=&pPars->pCoreHits[pPars->NumCoreHits-1];
	for(HitsChkd = 0; HitsChkd < NumHits2Chk; HitsChkd+=1, pCurCoreHit-=1)
		{
		if(pCurCoreHit->TargNodeID == TargNodeID &&
		  pCurCoreHit->ProbeNodeID == ProbeNodeID &&
			pCurCoreHit->flgRevCpl == (bRevCpl ? 1 : 0) &&
			(pCurCoreHit->ProbeOfs >= (ProbeOfs < ExpdHitLen ? 0 : ProbeOfs - ExpdHitLen) && pCurCoreHit->ProbeOfs <= (ProbeOfs + ExpdHitLen)) &&
			(pCurCoreHit->TargOfs >= (TargOfs < ExpdHitLen ? 0 : TargOfs - ExpdHitLen) && pCurCoreHit->TargOfs <= (TargOfs + ExpdHitLen)))
			return(0);
		}	
	}
 
pCurCoreHit = &pPars->pCoreHits[pPars->NumCoreHits++];

pCurCoreHit->ProbeNodeID = ProbeNodeID;
pCurCoreHit->flgRevCpl = bRevCpl ? 1 : 0;
pCurCoreHit->flgMulti = 0;
pCurCoreHit->ProbeOfs = ProbeOfs;
pCurCoreHit->TargNodeID = TargNodeID;
pCurCoreHit->HitLen = HitLen;
pCurCoreHit->TargOfs = TargOfs;
memset(&pCurCoreHit[1],0,sizeof(tsPBECCCoreHit));	// ensuring that used cores are always terminated with a marker end of cores initialised to 0
return(pPars->NumCoreHits);
}



// MapEntryID2NodeID
// Given a suffix array entry identifier returns the corresponding node identifier
uint32_t										// returned tsPBScaffNode node identifier
CPBECContigs::MapEntryID2NodeID(uint32_t EntryID)			// suffix array entry identifier
{
if(EntryID == 0 || EntryID > m_NumPBScaffNodes || m_pMapEntryID2NodeIDs == NULL)
	return(0);
return(m_pMapEntryID2NodeIDs[EntryID-1]);
}



int
CPBECContigs::IdentifyCoreHits(uint32_t HiConfSeqID,	// identify all overlaps of this probe sequence HiConfSeqID onto target sequences
				uint32_t MinTargLen,				// accepted target hit sequences must be at least this length
				uint32_t MaxTargLen,				// and if > 0 then accepted targets no longer than this length
				tsThreadPBECContigs *pPars)		// thread specific
{
int64_t PrevHitIdx;
int64_t NextHitIdx;
uint32_t HitEntryID;
uint32_t HitLoci;
uint32_t HitsThisCore;
uint32_t HighHitsThisCore;
uint32_t TotHitsAllCores;
uint32_t HitSeqLen;
uint32_t ProbeOfs;
uint32_t LastProbeOfs;

uint32_t ChkOvrLapCoreProbeOfs;
uint32_t LastCoreProbeOfs;
int ChkOvrLapCoreStartIdx;

etSeqBase *pCoreSeq;
tsPBECCScaffNode *pTargNode;

if(HiConfSeqID < 1 || HiConfSeqID > m_NumHiConfSeqs)
	return(eBSFerrParams);
uint32_t ProbeSeqLen = m_pSeqStore->GetLen(HiConfSeqID);

if(pPars->pProbeSeq == NULL || pPars->AllocdProbeSeqSize < (ProbeSeqLen + 10))
	{
	if(pPars->pProbeSeq != NULL)
		delete pPars->pProbeSeq;
	pPars->AllocdProbeSeqSize = ProbeSeqLen;
	pPars->pProbeSeq = new etSeqBase [ProbeSeqLen + 10];
	}

m_pSeqStore->GetSeq(HiConfSeqID,0,ProbeSeqLen,pPars->pProbeSeq);
pPars->pProbeSeq[ProbeSeqLen] = eBaseEOS;

if(pPars->bRevCpl)
	CSeqTrans::ReverseComplement(ProbeSeqLen,pPars->pProbeSeq);
pCoreSeq = pPars->pProbeSeq;
ChkOvrLapCoreProbeOfs = 0;
ChkOvrLapCoreStartIdx = 0;
LastCoreProbeOfs = 0;
PrevHitIdx = 0;
HitsThisCore = 0;
HighHitsThisCore = 0;
TotHitsAllCores = 0;
LastProbeOfs = 1 + ProbeSeqLen - pPars->CoreSeqLen;
if(pPars->CoreSeqLen < cMaxPacBioSeedExtn)
	LastProbeOfs -= 120;

if(MinTargLen < 1)
	MinTargLen = 1;

for(ProbeOfs = 0; ProbeOfs < LastProbeOfs; ProbeOfs+=pPars->DeltaCoreOfs,pCoreSeq+=pPars->DeltaCoreOfs)
	{
	PrevHitIdx = 0;
	HitsThisCore = 0;

    while((NextHitIdx = m_pSfxArray->IteratePacBio(pCoreSeq,ProbeSeqLen - ProbeOfs,pPars->CoreSeqLen,0,pPars->MinPBSeqLen,MaxTargLen,PrevHitIdx,&HitEntryID,&HitLoci)) > 0)
		{
		PrevHitIdx = NextHitIdx;
		pTargNode = &m_pPBScaffNodes[MapEntryID2NodeID(HitEntryID)-1];
		HitSeqLen = pTargNode->SeqLen; 

		if( pTargNode->flgUnderlength == 1 ||	// not interested in under length targets
			 HitSeqLen < (uint32_t)pPars->MinPBSeqLen)		// not interested if target sequence length not in range of accepted target sequence length to be processed
			continue;

  		if(AddCoreHit(HiConfSeqID,pPars->bRevCpl,ProbeOfs,pTargNode->NodeID,HitLoci,pPars->CoreSeqLen,pPars)>0)
			{
			HitsThisCore += 1;
			if(HitsThisCore > pPars->MaxAcceptHitsPerSeedCore)
				break;
			}
		}
	if(HitsThisCore)	// if at least one hit from this core
		{
		if(HitsThisCore > HighHitsThisCore)
			HighHitsThisCore = HitsThisCore;
		TotHitsAllCores += HitsThisCore;
		}
	}
if(pPars->bRevCpl)
	CSeqTrans::ReverseComplement(ProbeSeqLen,pPars->pProbeSeq);

return(pPars->NumCoreHits);
}




// SortLenDescending
// Sort scaffolding nodes by length descending with entry identifiers as tie breaker
int
CPBECContigs::SortLenDescending(const void *arg1, const void *arg2)
{
tsPBECCScaffNode *pEl1 = (tsPBECCScaffNode *)arg1;
tsPBECCScaffNode *pEl2 = (tsPBECCScaffNode *)arg2;

if(pEl1->SeqLen < pEl2->SeqLen)
	return(1);
if(pEl1->SeqLen > pEl2->SeqLen)
	return(-1);
if(pEl1->EntryID < pEl2->EntryID)	
	return(-1);
if(pEl1->EntryID > pEl2->EntryID)
	return(1);
return(0);
}


// SortCoreHitsByTargNodeID
// Sort core hits by ProbeNodeID.TargNodeID.TargOfs.ProbeOfs.flgRevCpl ascending
int
CPBECContigs::SortCoreHitsByTargProbeOfs(const void *arg1, const void *arg2)
{
tsPBECCCoreHit *pEl1 = (tsPBECCCoreHit *)arg1;
tsPBECCCoreHit *pEl2 = (tsPBECCCoreHit *)arg2;

if(pEl1->ProbeNodeID < pEl2->ProbeNodeID)
	return(-1);
if(pEl1->ProbeNodeID > pEl2->ProbeNodeID)
	return(1);
if(pEl1->TargNodeID < pEl2->TargNodeID)
	return(-1);
if(pEl1->TargNodeID > pEl2->TargNodeID)
	return(1);
if(pEl1->TargOfs < pEl2->TargOfs)	
	return(-1);
if(pEl1->TargOfs > pEl2->TargOfs)
	return(1);
if(pEl1->ProbeOfs < pEl2->ProbeOfs)	
	return(-1);
if(pEl1->ProbeOfs > pEl2->ProbeOfs)
	return(1);
if(pEl1->flgRevCpl != pEl2->flgRevCpl)
	return(1);
return(0);
}

// SortCoreHitsByProbeNodeID
// Sort core hits by ProbeNodeID.TargNodeID.ProbeOfs.TargOfs.flgRevCpl ascending
int
CPBECContigs::SortCoreHitsByProbeTargOfs(const void *arg1, const void *arg2)
{
tsPBECCCoreHit *pEl1 = (tsPBECCCoreHit *)arg1;
tsPBECCCoreHit *pEl2 = (tsPBECCCoreHit *)arg2;

if(pEl1->ProbeNodeID < pEl2->ProbeNodeID)
	return(-1);
if(pEl1->ProbeNodeID > pEl2->ProbeNodeID)
	return(1);
if(pEl1->TargNodeID < pEl2->TargNodeID)
	return(-1);
if(pEl1->TargNodeID > pEl2->TargNodeID)
	return(1);
if(pEl1->ProbeOfs < pEl2->ProbeOfs)	
	return(-1);
if(pEl1->ProbeOfs > pEl2->ProbeOfs)
	return(1);
if(pEl1->TargOfs < pEl2->TargOfs)	
	return(-1);
if(pEl1->TargOfs > pEl2->TargOfs)
	return(1);
if(pEl1->flgRevCpl != pEl2->flgRevCpl)
	return(1);
return(0);
}

// SortCoreHitsDescending
// Sort target core hits by number of hits descending
int
CPBECContigs::SortCoreHitsDescending(const void *arg1, const void *arg2)
{
sPBECCCoreHitCnts *pEl1 = (sPBECCCoreHitCnts *)arg1;
sPBECCCoreHitCnts *pEl2 = (sPBECCCoreHitCnts *)arg2;

uint32_t El1NumHits;
uint32_t El2NumHits;
El1NumHits = max(pEl1->NumSHits,pEl1->NumAHits);
El2NumHits = max(pEl2->NumSHits,pEl2->NumAHits);
if(El1NumHits > El2NumHits)
	return(-1);
if(El1NumHits < El2NumHits)
	return(1);
return(0);
}





