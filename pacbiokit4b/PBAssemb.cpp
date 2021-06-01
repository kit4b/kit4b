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
#include "PBAssemb.h"

#define PBASCAFFSENSEANTI 1

int
ProcPacBioAssemb(etPBPMode PMode,			// processing mode
		int MinScaffSeqLen,					// individual scaffold sequences must be of at least this length (defaults to 5Kbp)
		int MinScaffOverlap,				// pairs of targeted scaffold sequences must have overlapped by at least this many bp to be considered for merging into a longer scaffold sequence (defaults to 5Kbp) 
	    int Min1kScore,                     // minimum normalised 1Kbp overlap score
		bool bAcceptOrphanSeqs,				// also report sequences which are not overlapped or overlapping any other sequence
		bool bSenseOnlyOvlps,				// if false then both sense/sense and sense/antisense overlaps will be accepted and processed, otherwise sense/sense only overlaps accepted and processed
		char *pszOverlapDetailFile,			// pregenerated sequence overlap loci details
		int NumErrCorrectedFiles,			// number of error corrected sequence specs
		char *pszErrCorrectedFiles[],		// input error corrected sequence files
	    char *pszOutFile,					// where to write merged scaffolded sequences
		int NumThreads);					// maximum number of worker threads to use


#ifdef _WIN32
int ProcAssemb(int argc, char* argv[])
{
// determine my process name
_splitpath(argv[0],NULL,NULL,gszProcName,NULL);
#else
int
ProcAssemb(int argc, char** argv)
{
// determine my process name
CUtility::splitpath((char *)argv[0],NULL,gszProcName);
#endif

int iFileLogLevel;			// level of file diagnostics
int iScreenLogLevel;		// level of file diagnostics
char szLogFile[_MAX_PATH];	// write diagnostics to this file
int Idx;                    // file index 
int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure

int PMode;					// processing mode

int MinScaffSeqLen;			// individual target scaffold sequences must be of at least this length (defaults to 5Kbp)
int MinScaffOverlap;		// pairs of targeted sequences must overlap by at least this many bp to be considered for merging into a longer scaffold sequence (defaults to 5Kbp) 
int Min1kScore;             // minimum normalised 1Kbp overlap score
bool bAcceptOrphanSeqs;		// also report sequences which are not overlapped or overlapping any other sequence
bool bSenseOnlyOvlps;	    // if false then both sense/sense and sense/antisense overlaps will be accepted and processed, otherwise sense/sense only overlaps accepted and processed

char szOverlapDetailFile[_MAX_PATH];	// pregenerated  overlap loci details
int NumPacBioFiles;						// number of input pacbio file specs
char *pszPacBioFiles[cMaxInFileSpecs];	// input pacbio files

char szOutFile[_MAX_PATH];				// where to write merged assembled sequences

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

struct arg_int *pmode = arg_int0("m","pmode","<int>",				"processing mode - 0 use pregenerated overlap loci details file and error corrected sequences, 1 generate GEXF output format, 2 generate GraphML output format");
struct arg_int *minscafflen = arg_int0("l","minscafflen","<int>",	"minimum individual sequence length (default 5000, range 500 to 100000)");
struct arg_int *minscaffovl = arg_int0("L","minscaffovl","<int>",	"minimum overlap required to merge reads into single contig (default 5000, range 250 to 100000)");

struct arg_int *min1kscore = arg_int0("s", "min1kscore", "<int>", "minimum normalised 1Kbp overlap score required (default 980, range 800 to 1000)");
struct arg_lit  *orphanseqs = arg_lit0("S", "orphanseqs", "accept orphan sequences not overlapped by any other sequence");

#ifdef PBASCAFFSENSEANTI
struct arg_lit  *senseonlyovlps = arg_lit0("a", "senseonlyovlps",		  "accepting sense only overlaps (default is to accept both sense/antisense and sense/sense)");
#else
struct arg_lit  *senseonlyovlps = arg_lit0("a", "senseonlyovlps",		  "currently accepting sense only overlaps if scaffolding");
#endif
struct arg_file *pacbiosovlps = arg_file1("I","pacbiosovlps","<file>",	"input file containing pregenerated error corrected overlap detail");

struct arg_file *pacbiofiles = arg_filen("i","pacbiofile","<file>",1,cMaxInFileSpecs,	"names of input files containing error corrected sequences to be used for contigs");
struct arg_file *outfile = arg_file1("o","out","<file>",					"output merged contig sequences (processing mode 0) or GEFX or GraphML format to this file");

struct arg_int *threads = arg_int0("T","threads","<int>",					"number of processing threads 0..128 (defaults to 0 which sets threads to number of CPU cores)");

struct arg_file *summrslts = arg_file0("q","sumrslts","<file>",				"Output results summary to this SQLite3 database file");
struct arg_str *experimentname = arg_str0("w","experimentname","<str>",		"experiment name SQLite3 database file");
struct arg_str *experimentdescr = arg_str0("W","experimentdescr","<str>",	"experiment description SQLite3 database file");

struct arg_end *end = arg_end(200);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					pmode,minscafflen,minscaffovl,min1kscore,orphanseqs,senseonlyovlps,
					summrslts,pacbiosovlps,pacbiofiles,experimentname,experimentdescr,
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

	PMode = (etPBPMode)(pmode->count ? pmode->ival[0] : (int)ePBPMScaffold);
	if(PMode < ePBPMScaffold || PMode > ePBPMToGraphML)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Processing mode '-m%d' must be in range 0..%d",PMode,ePBPMToGraphML);
		return(1);
		}

	if(PMode == ePBPMScaffold)
#ifdef PBASCAFFSENSEANTI
		bSenseOnlyOvlps = senseonlyovlps->count ? true : false;
#else
		bSenseOnlyOvlps = true;
#endif
	else
		bSenseOnlyOvlps = false;

	MinScaffSeqLen = minscafflen->count ? minscafflen->ival[0] : cDfltMinErrCorrectLen;
	if(MinScaffSeqLen < cMinPBSeqLen || MinScaffSeqLen > cMaxMinPBSeqLen)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Minimum accepted sequence length '-l%d' must be in range %d..%dbp",MinScaffSeqLen,cMinPBSeqLen,cMaxMinPBSeqLen);
		return(1);
		}

	MinScaffOverlap = minscaffovl->count ? minscaffovl->ival[0] : cDfltMinErrCorrectLen;
	if(MinScaffOverlap < cMinPBSeqLen || MinScaffOverlap > cMaxMinPBSeqLen)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Minimum sequence overlap required length '-L%d' must be in range %d..%dbp",MinScaffOverlap,cMinPBSeqLen,cMaxMinPBSeqLen);
		return(1);
		}

	Min1kScore = min1kscore->count ? min1kscore->ival[0] : cDfltMin1kScore;
	if (Min1kScore < cMinMin1kScore || Min1kScore > cMaxMin1kScore)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Minimum sequence 1Kbp normalised overlap score '-s%d' must be in range %d..%d", Min1kScore, cMinMin1kScore, cMaxMin1kScore);
		return(1);
		}

	bAcceptOrphanSeqs = orphanseqs->count ? true : false;

	for(NumPacBioFiles=Idx=0;NumPacBioFiles < cMaxInFileSpecs && Idx < pacbiofiles->count; Idx++)
			{
			pszPacBioFiles[Idx] = NULL;
			if(pszPacBioFiles[NumPacBioFiles] == NULL)
				pszPacBioFiles[NumPacBioFiles] = new char [_MAX_PATH];
			strncpy(pszPacBioFiles[NumPacBioFiles],pacbiofiles->filename[Idx],_MAX_PATH);
			pszPacBioFiles[NumPacBioFiles][_MAX_PATH-1] = '\0';
			CUtility::TrimQuotedWhitespcExtd(pszPacBioFiles[NumPacBioFiles]);
			if(pszPacBioFiles[NumPacBioFiles][0] != '\0')
				NumPacBioFiles++;
			}

	if(!NumPacBioFiles)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no input PacBio file(s) specified with '-i<filespec>' option)\n");
		exit(1);
		}

	strncpy(szOverlapDetailFile,pacbiosovlps->filename[0],sizeof(szOverlapDetailFile));
	szOverlapDetailFile[sizeof(szOverlapDetailFile)-1] = '\0';
	CUtility::TrimQuotedWhitespcExtd(szOverlapDetailFile);
	if(szOverlapDetailFile[0] == '\0')
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no input pregenerated sequence overlap loci detail file specified with '-I<filespec>' option)\n");
		exit(1);
		}

	strncpy(szOutFile,outfile->filename[0],sizeof(szOutFile));
	szOutFile[sizeof(szOutFile)-1] = '\0';
	CUtility::TrimQuotedWhitespcExtd(szOutFile);
	if(szOutFile[0] == '\0')
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no output assembly file specified with '-o<filespec>' option)\n");
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

	char *pszMode;
	switch(PMode) {
		case ePBPMScaffold:		// scaffolding
			pszMode = (char *)"Using pregenerated overlaps for scaffolding";
			break;
		case ePBPMToGEFX:		// converting input overlap detail into GEFX ready for graph visualisation
			pszMode = (char *)"Converting overlaps into GEXF ready for graph visualisation";
			break;
		case ePBPMToGraphML:		// converting input overlap detail into GraphML ready for graph visualisation
			pszMode = (char *)"Converting overlaps into GraphML ready for graph visualisation";
			break;
			}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode: '%s'",pszMode);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum individual input sequence length: %dbp",MinScaffSeqLen);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum sequence overlap required to merge into single config: %d",MinScaffOverlap);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum 1Kbp normalised overlap score: %d", Min1kScore);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Accepting orphan sequences: '%s'", bAcceptOrphanSeqs ? "Yes" : "No");
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Accepting overlaps : '%s'", bSenseOnlyOvlps ? "Sense/Sense" : "Sense/Sense and Sense/AntiSense");


	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input file containing pregenerated PacBio overlap loci detail: '%s'",szOverlapDetailFile);

	for(Idx=0; Idx < NumPacBioFiles; Idx++)
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input sequences file spec: '%s'",pszPacBioFiles[Idx]);

	if(PMode == ePBPMScaffold)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output contig sequences file: '%s'",szOutFile);
	else
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output GEFX file (sense overlap sense only reported): '%s'",szOutFile);

	if(szExperimentName[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"This processing reference: %s",szExperimentName);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"number of threads : %d",NumThreads);

	if(gExperimentID > 0)
		{
		int ParamID;
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szLogFile),"log",szLogFile);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(PMode),"pmode",&PMode);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MinScaffSeqLen),"minscafflen",&MinScaffSeqLen);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MinScaffOverlap),"minscaffovl",&MinScaffOverlap);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID, ePTInt32, (int)sizeof(Min1kScore), "min1kscore", &Min1kScore);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID, ePTBool, (int)sizeof(bAcceptOrphanSeqs), "orphanseqs", &bAcceptOrphanSeqs);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID, ePTBool, (int)sizeof(bSenseOnlyOvlps), "senseonlyovlps", &bSenseOnlyOvlps);

		for(Idx=0; Idx < NumPacBioFiles; Idx++)
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(pszPacBioFiles[Idx]),"pacbiofile",pszPacBioFiles[Idx]);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szOutFile),"out",szOutFile);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szOverlapDetailFile),"pacbiosovlps",szOverlapDetailFile);
		
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
	Rslt = ProcPacBioAssemb((etPBPMode)PMode,MinScaffSeqLen,MinScaffOverlap, Min1kScore, bAcceptOrphanSeqs,bSenseOnlyOvlps,szOverlapDetailFile,NumPacBioFiles,pszPacBioFiles,szOutFile,NumThreads);
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
ProcPacBioAssemb(etPBPMode PMode,		// processing mode
		int MinScaffSeqLen,			// individual scaffold sequences must be of at least this length (defaults to 5Kbp)
		int MinScaffOverlap,		// pairs of targeted scaffold sequences must overlap by at least this many bp to be considered for merging into a longer scaffold sequence (defaults to 5Kbp) 
	    int Min1kScore,             // minimum normalised 1Kbp overlap score
	    bool bAcceptOrphanSeqs,		// also report sequences which are not overlapped or overlapping any other sequence
		bool bSenseOnlyOvlps,				// if false then both sense/sense and sense/antisense overlaps will be accepted and processed, otherwise sense/sense only overlaps accepted and processed
		char *pszOverlapDetailFile,	// pregenerated sequence overlap loci details
		int NumErrCorrectedFiles,	// number of error corrected sequence specs
		char *pszErrCorrectedFiles[],		// input error corrected sequence files
	    char *pszOutFile,			// where to write merged scaffolded sequences
		int NumThreads)				// maximum number of worker threads to use
{
int Rslt;
CPBAssemb *pPacBioer;

if((pPacBioer = new CPBAssemb)==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: Unable to instantiate CPBAssemb");
	return(eBSFerrObj);
	}

Rslt = pPacBioer->Process(PMode,MinScaffSeqLen,MinScaffOverlap, Min1kScore, bAcceptOrphanSeqs,bSenseOnlyOvlps,
								pszOverlapDetailFile,NumErrCorrectedFiles,pszErrCorrectedFiles,pszOutFile,NumThreads);
delete pPacBioer;
return(Rslt);
}

CPBAssemb::CPBAssemb() // relies on base classes constructors
{
m_pSeqStore = NULL;
m_pAssembGraph = NULL;
m_pPBScaffNodes = NULL;
m_pMapEntryID2NodeIDs = NULL;
m_bMutexesCreated = false;
Init();
}

CPBAssemb::~CPBAssemb() // relies on base classes destructor
{
Reset(false);
}


void
CPBAssemb::Init(void)
{
if(m_pSeqStore != NULL)
	{
	delete m_pSeqStore;
	m_pSeqStore = NULL;
	}

if(m_pAssembGraph != NULL)
	{
	delete m_pAssembGraph;
	m_pAssembGraph = NULL;
	}

if(m_pPBScaffNodes != NULL)
	{
	delete m_pPBScaffNodes;
	m_pPBScaffNodes = NULL;
	}
if(m_pMapEntryID2NodeIDs != NULL)
	{
	delete m_pMapEntryID2NodeIDs;
	m_pMapEntryID2NodeIDs = NULL;
	}

m_NumPBScaffNodes = 0;
m_AllocdPBScaffNodes = 0;

m_NumOverlapProcessed = 0;
m_ProvOverlapping = 0;
m_ProvOverlapped = 0;
m_ProvContained = 0;
m_ProvArtefact = 0;
m_ProvSWchecked = 0;

m_PMode = ePBPMScaffold;

m_OverlapFloat = cDfltScaffMaxOverlapFloat;
m_MinScaffSeqLen = cDfltMinErrCorrectLen;	
m_MinScaffOverlap = cDfltMinErrCorrectLen;

m_ScaffScoreExact = cDfltScaffScoreExact;					
m_ScaffScoreMismatch = cDfltScaffScoreMismatch;				
m_ScaffScoreGapOpen = cDfltScaffScoreGapOpen;				
m_ScaffScoreGapExtn = cDfltScaffScoreGapExtn;				
m_MinScaffScoreThres = cDfltMin1kScore;
m_bAcceptOrphanSeqs = false;
m_bSenseOnlyOvlps = false;
m_NumRejectedMinSeqLen = 0;
m_NumRejectedMinScaffOverlap = 0;
m_NumRejectedScoreThres = 0;
m_NumRejectContained = 0;
m_NumRejectAntisense = 0;
m_NumRejectArtefact = 0;
m_NumAcceptedOverlaps = 0;		

m_NumErrCorrectedFiles = 0;
memset(m_szErrCorrectedFiles,0,sizeof(m_szErrCorrectedFiles));

m_szOutFile[0] = '\0';	

m_NumThreads = 0;
if(m_bMutexesCreated)
	DeleteMutexes();
m_bMutexesCreated = false; 
}

void
CPBAssemb::Reset(bool bSync)			// if bSync true then fsync before closing output file handles
{
Init();
}


int
CPBAssemb::CreateMutexes(void)
{
if(m_bMutexesCreated)
	return(eBSFSuccess);

m_CASSerialise = 0;
m_CASLock = 0;

m_bMutexesCreated = true;
return(eBSFSuccess);
}

void
CPBAssemb::DeleteMutexes(void)
{
if(!m_bMutexesCreated)
	return;
m_bMutexesCreated = false;
}


void
CPBAssemb::AcquireCASSerialise(void)
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
CPBAssemb::ReleaseCASSerialise(void)
{
#ifdef _WIN32
InterlockedCompareExchange(&m_CASSerialise,0,1);
#else
__sync_val_compare_and_swap(&m_CASSerialise,1,0);
#endif
}


void
CPBAssemb::AcquireCASLock(void)
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
CPBAssemb::ReleaseCASLock(void)
{
#ifdef _WIN32
InterlockedCompareExchange(&m_CASLock,0,1);
#else
__sync_val_compare_and_swap(&m_CASLock,1,0);
#endif
}



uint32_t												//  returns number of overlaps loaded and accepted, if > cMaxValidID then cast to teBSFrsltCodes for actual error 
CPBAssemb::LoadPacBioOvlps(char *pszPacBioOvlps,			// parse and load pregenerated PacBio sequence overlap loci CSV file 
						bool bValidateOnly,		// true if parse and validate only
						bool bSenseOnlyOvlps)				// if false then both sense/sense and sense/antisense overlaps will be accepted and processed, otherwise sense/sense only overlaps accepted and processed
{
int Rslt;
int NumFields;
uint32_t NumElsRead;
uint32_t NumEdges;
int NumNoProbeIdents;
int NumNoTargIdents;
int ProbeNameID;
int TargNameID;
int ScaffoldScore;
int ScoreAlignLen;
int Scaffold1KScore;	

char *pszTmp;
int Class;
int ProbeID;
char szProbeDescr[cMaxDescrIDLen+1];

int PrevProbeIdent;
char szPrevProbeDescr[cMaxDescrIDLen+1];

int TargID;
char szTargDescr[cMaxDescrIDLen+1];
int SeedHits;
char ProbeSense;
char TargSense;
int ProbeLen;
int TargLen;
int	ProbeAlignLength;
int	TargAlignLength;
int	PeakScore;
int FinalScore;
int	NumAlignedBases;
int	NumExactBases;
int	NumProbeInserts;
int	NumProbeInsertBases;
int	NumTargInserts;
int	NumTargInsertBases;
int	ProbeStartOfs;
int	TargStartOfs;
int	ProbeEndOfs;
int	TargEndOfs;
int	ProbeOfs5;
int	TargOfs5;
int	ProbeOfs3;
int	TargOfs3;

uint32_t ProbeSeqLen;
uint32_t TargSeqLen;

m_NumRejectedScoreThres = 0;
m_NumRejectedMinSeqLen = 0;
m_NumRejectedMinScaffOverlap = 0;
m_NumRejectContained = 0;
m_NumRejectAntisense = 0;
m_NumRejectArtefact = 0;
m_NumAcceptedOverlaps = 0;

if(!bValidateOnly && m_pAssembGraph == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected m_pAssembGraph to have been instantiated");
	return(eBSFerrObj);
	}

CCSVFile *pCSV = new CCSVFile;
if(pCSV == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CCSVfile");
	return((uint32_t)eBSFerrObj);
	}

if((Rslt=pCSV->Open(pszPacBioOvlps))!=eBSFSuccess)
	{
	while(pCSV->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pCSV->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open file: %s",pszPacBioOvlps);
	delete pCSV;
	return((uint32_t)Rslt);
	}
PrevProbeIdent = 0;
szPrevProbeDescr[0] = '\0';
NumNoProbeIdents = 0;
NumNoTargIdents = 0;
NumElsRead = 0;
NumEdges = 0;
while((Rslt=pCSV->NextLine()) > 0)			// onto next line containing fields
	{
	NumFields = pCSV->GetCurFields();		// PacBio scaffold overlap files contain 28 fields
	if(NumFields != 28)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected PacBio overlaps CSV file number of fields to be 28 in '%s', GetCurFields() returned '%d'",pszPacBioOvlps,NumFields);
		delete pCSV;
		return((uint32_t)eBSFerrFieldCnt);
		}
	if(!NumElsRead && pCSV->IsLikelyHeaderLine())
		continue;
	NumElsRead += 1;

	if((Rslt=pCSV->GetInt(1,&Class)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing Class at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((uint32_t)Rslt);
		}
	if((Rslt=pCSV->GetInt(2,&ProbeID)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing ProbeID at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((uint32_t)Rslt);
		};
	if((Rslt=pCSV->GetText(3,&pszTmp)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing ProbeName at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((uint32_t)Rslt);
		};
	strncpy(szProbeDescr,pszTmp,sizeof(szProbeDescr));
	szProbeDescr[sizeof(szProbeDescr)-1] = '\0';
	if((Rslt=pCSV->GetInt(4,&TargID)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing TargID at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((uint32_t)Rslt);
		};
	if((Rslt=pCSV->GetText(5,&pszTmp)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing TargName at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((uint32_t)Rslt);
		};
	strncpy(szTargDescr,pszTmp,sizeof(szTargDescr));
	szTargDescr[sizeof(szTargDescr)-1] = '\0';
	if((Rslt=pCSV->GetInt(6,&SeedHits)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing SeedHits at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((uint32_t)Rslt);
		};
	if((Rslt=pCSV->GetText(7,&pszTmp)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing ProbeSense at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((uint32_t)Rslt);
		};
	ProbeSense = pszTmp[0];
	switch(ProbeSense) {
		case 's': case 'S': break;
		case 'a': case 'A': break;
		default:
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing ProbeSense at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
			delete pCSV;
			return((uint32_t)eBSFerrParse);
		};

	if((Rslt=pCSV->GetText(8,&pszTmp)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing TargSense at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((uint32_t)Rslt);
		};
	TargSense = pszTmp[0];
	switch(TargSense) {
		case 's': case 'S': break;
		case 'a': case 'A': break;
		default:
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing TargSense at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
			delete pCSV;
			return((uint32_t)eBSFerrParse);
		};

	if((Rslt=pCSV->GetInt(9,&ProbeLen)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing ProbeLen at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((uint32_t)Rslt);
		};
	if((Rslt=pCSV->GetInt(10,&TargLen)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing TargLen at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((uint32_t)Rslt);
		};
	if((Rslt=pCSV->GetInt(11,&ProbeAlignLength)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing ProbeAlignLength at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((uint32_t)Rslt);
		};
	if((Rslt=pCSV->GetInt(12,&TargAlignLength)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing TargAlignLength at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((uint32_t)Rslt);
		};
	if((Rslt=pCSV->GetInt(13,&PeakScore)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing PeakScore at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((uint32_t)Rslt);
		};
	if((Rslt=pCSV->GetInt(14,&FinalScore)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing FinalScore at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((uint32_t)Rslt);
		};
	if((Rslt=pCSV->GetInt(15,&NumAlignedBases)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing NumAlignedBases at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((uint32_t)Rslt);
		};
	if((Rslt=pCSV->GetInt(16,&NumExactBases)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing NumExactBases at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((uint32_t)Rslt);
		};
	if((Rslt=pCSV->GetInt(17,&NumProbeInserts)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing NumProbeInserts at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((uint32_t)Rslt);
		};
	if((Rslt=pCSV->GetInt(18,&NumProbeInsertBases)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing NumProbeInsertBases at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((uint32_t)Rslt);
		};
	if((Rslt=pCSV->GetInt(19,&NumTargInserts)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing NumTargInserts at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((uint32_t)Rslt);
		};
	if((Rslt=pCSV->GetInt(20,&NumTargInsertBases)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing NumTargInsertBases at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((uint32_t)Rslt);
		};
	if((Rslt=pCSV->GetInt(21,&ProbeStartOfs)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing ProbeStartOfs at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((uint32_t)Rslt);
		};

	if((Rslt=pCSV->GetInt(22,&TargStartOfs)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing TargStartOfs at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((uint32_t)Rslt);
		};
	if((Rslt=pCSV->GetInt(23,&ProbeEndOfs)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing ProbeEndOfs at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((uint32_t)Rslt);
		};
	if((Rslt=pCSV->GetInt(24,&TargEndOfs)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing TargEndOfs at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((uint32_t)Rslt);
		};
	if((Rslt=pCSV->GetInt(25,&ProbeOfs5)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing ProbeOfs5 at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((uint32_t)Rslt);
		};
	if((Rslt=pCSV->GetInt(26,&TargOfs5)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing TargOfs5 at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((uint32_t)Rslt);
		};
	if((Rslt=pCSV->GetInt(27,&ProbeOfs3)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing ProbeOfs3 at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((uint32_t)Rslt);
		};
	if((Rslt=pCSV->GetInt(28,&TargOfs3)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error parsing TargOfs3 at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((uint32_t)Rslt);
		};

	// check for consistency in the PacBio read overlaps before edge processing
	// only accepted as consistent if sequence lengths are in the range of 1000..1000000, alignments are at least 100bp, starts+100 < ends, start/ends within read sequence lengths, and SW scores less than 100000000
	// note that all offsets are 1 based 
	if(ProbeID < 1 || TargID < 1 || ProbeID == TargID ||  // not expecting probe to be self overlapping!
	   PeakScore < 0 || PeakScore > 10000000 ||  FinalScore < 0 || FinalScore > PeakScore ||         // SW scores in the millions surely is an error????
	   Class < (int)eOLCOverlapping || Class > (int)eOLCartefact ||	// overlap classification expected to be in this range	
	   ProbeLen < 1000 || ProbeLen > 100000000 || TargLen < 1000 || TargLen > 100000000 ||    // allowing for PacBio reads of 1Kbp upto 100Mbp!	
	   ProbeStartOfs < 1 || (ProbeStartOfs + 100) > ProbeLen ||  // expecting alignments to be at least 100bp and alignment start/end within the read length
	   ProbeEndOfs < (ProbeStartOfs + 99) || ProbeEndOfs > ProbeLen || 
	   TargStartOfs < 1 || (TargStartOfs + 100) > TargLen ||  
	   TargEndOfs < (TargStartOfs + 99) || TargEndOfs > TargLen ||
	   !(ProbeSense == 'S' || ProbeSense == 's' || ProbeSense == 'A' || ProbeSense == 'a') ||
		!(TargSense == 'S' || TargSense == 's' || TargSense == 'A' || TargSense == 'a')) 
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error, inconsistencies in values parsed from line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((uint32_t)eBSFerrParse);
		}
	if(ProbeSense == 'A' || ProbeSense == 'a')
		{
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"Sorry, can not accept (treating as artefact) ProbeSense == 'A' at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		m_NumRejectArtefact += 1;
		continue;
		}

	switch(Class) {
		case eOLCartefact:		// slough any known previously classified as artefact
			m_NumRejectArtefact += 1;
			continue;

		case eOLCcontained:     // slough contained or containing edges; can't extend these
		case eOLCcontains:
			m_NumRejectContained += 1;
			continue;

		default:
			if(bSenseOnlyOvlps && (TargSense == 'A' || TargSense == 'a'))
				{
				m_NumRejectAntisense += 1;
				continue;
				}
			break;
		}


		// Score for scaffolding needs to be independent of the SW score as user may have changed the default scoring penalties when error correcting.
        // Also because scaffolding is normally with sequences which are error reduced relative to the original reads then the scoring used should be
        // more akin to that normally associated with relatively error free alignments where the InDel rate is greatly reduced compared with the original PacBio sequences.
        // Currently the default scoring is 1 for exacts, -1 for mismatches and -3 for InDels and -1 for InDel gap extensions; these may need to be exposed to the user as parameterised ....
	ScoreAlignLen = ((1 + ProbeAlignLength + TargAlignLength)/2);
	if(ScoreAlignLen < (int)m_MinScaffOverlap)
		{
		m_NumRejectedMinScaffOverlap += 1;
		continue;
		}

	ScaffoldScore = (NumExactBases * m_ScaffScoreExact) +							// score the exacts
				((NumAlignedBases - NumExactBases) * m_ScaffScoreMismatch) +		// score the mismatches
				((NumProbeInserts + NumTargInserts) *  m_ScaffScoreGapOpen) +     // score the gap (InDel) openings
				((NumProbeInsertBases - NumProbeInserts + NumTargInsertBases - NumTargInserts) *  m_ScaffScoreGapExtn); // score the gap extensions
	Scaffold1KScore = (int)((((int64_t)ScaffoldScore * 1000) + 500) / (1+ScoreAlignLen)); // rounding up

	if(Scaffold1KScore < m_MinScaffScoreThres)
		{
		m_NumRejectedScoreThres += 1;
		continue;
		}

	if(szPrevProbeDescr[0] == 0 || stricmp(szPrevProbeDescr,szProbeDescr) != 0)
		{
		ProbeNameID = m_pSeqStore->GetSeqID(szProbeDescr);
		if(ProbeNameID <= 0)
			{
			if(!NumNoProbeIdents)
				gDiagnostics.DiagOut(eDLWarn,gszProcName,"Unable to locate probe sequence named '%s' in loaded sequences at line %d in file '%s'",szProbeDescr,NumElsRead,pszPacBioOvlps);
			NumNoProbeIdents += 1;
			continue;
			}
		strcpy(szPrevProbeDescr,szProbeDescr);
		PrevProbeIdent = ProbeNameID;
		}
	else
		ProbeNameID = PrevProbeIdent;

	 ProbeSeqLen = m_pSeqStore->GetLen(ProbeNameID);
	 if(ProbeSeqLen != ProbeLen)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Actual probe fasta sequence length: %d, Overlap detail probe length: %d, for sequence named '%s' at line %d in file '%s'",
									ProbeSeqLen,ProbeLen,szProbeDescr,NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((uint32_t)eBSFerrParse);
		}


	TargNameID = m_pSeqStore->GetSeqID(szTargDescr);
	if(TargNameID <= 0)
		{
		if(!NumNoTargIdents)
			gDiagnostics.DiagOut(eDLWarn,gszProcName,"Unable to locate target sequence named '%s' in loaded sequences at line %d in file '%s'",szTargDescr,NumElsRead,pszPacBioOvlps);
		NumNoTargIdents += 1;
		continue;
		}

	 TargSeqLen = m_pSeqStore->GetLen(TargNameID);
	 if(TargSeqLen != TargLen)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Actual target fasta sequence length: %d, Overlap detail target length: %d, for sequence named '%s' at line %d in file '%s'",
									TargSeqLen,TargLen,szProbeDescr,NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return((uint32_t)eBSFerrParse);
		}

	if(ProbeLen < (int)m_MinScaffSeqLen || TargLen < (int)m_MinScaffSeqLen)
		{
		m_NumRejectedMinSeqLen += 1;
		continue;
		}

	m_NumAcceptedOverlaps += 1;

	if(bValidateOnly)
		continue;

	// note that AddEdge() is expecting offsets to be 0 based hence substracting 1
	if((NumEdges = m_pAssembGraph->AddEdge(ProbeNameID,TargNameID,ProbeLen,TargLen,(uint32_t)Scaffold1KScore,ScoreAlignLen,ProbeStartOfs-1,ProbeEndOfs-1,TargStartOfs-1,TargEndOfs-1,(eOverlapClass)Class,TargSense == 's' || TargSense == 'S' ? false : true)) > cMaxValidID)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddEdge failed at line %d in file '%s'",NumElsRead,pszPacBioOvlps);
		delete pCSV;
		return(NumEdges);
		}
	}

delete pCSV;


if(NumNoProbeIdents)
	gDiagnostics.DiagOut(eDLWarn,gszProcName,"Unable to locate %d probe sequences in loaded sequences in file '%s'",NumNoProbeIdents,pszPacBioOvlps);
if(NumNoTargIdents)
	gDiagnostics.DiagOut(eDLWarn,gszProcName,"Unable to locate %d target sequences in loaded sequences in file '%s'",NumNoTargIdents,pszPacBioOvlps);

if(!bValidateOnly)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Accepted or inferred %u overlaps for further processing (rejected %u min seq len, %u overlap len, %u score, %u contained, %u antisense, %u artefact) from file '%s'",
												m_NumAcceptedOverlaps,m_NumRejectedMinSeqLen,m_NumRejectedMinScaffOverlap,m_NumRejectedScoreThres,m_NumRejectContained,m_NumRejectAntisense,m_NumRejectArtefact,pszPacBioOvlps);

	m_NumAcceptedOverlaps = NumEdges;
	}
return(m_NumAcceptedOverlaps);
}


// ProcessFastaFile
// Parse input fasta format file into a CSeqStore
int
CPBAssemb::ProcessFastaFile(int MinSeqLen,		// only accept for indexing sequences of at least this length
				char *pszFile)						// file containing sequences
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
int Descrlen;
bool bFirstEntry;
bool bEntryCreated;
int Rslt;
int SeqID;
int NumSeqsAccepted;
size_t TotAcceptedLen;
uint32_t NumSeqsUnderlength;

if(MinSeqLen <= cMinSeqLen)			// set a floor of cMinSeqLen on the minimum accepted sequence length
	MinSeqLen = cMinSeqLen;

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
				if((Rslt=m_pSeqStore->AddSeq(0,szName,(uint32_t)BuffOfs,pSeqBuff)) == 0)
					{
					Rslt = -1;
					break;
					}
				else
					{
					NumSeqsAccepted += 1;
					TotAcceptedLen += BuffOfs;
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
		if((Rslt=m_pSeqStore->AddSeq(0,szName,(uint32_t)BuffOfs,pSeqBuff)) == 0)
			Rslt = -1;
		else
			{
			Rslt = eBSFSuccess;
			NumSeqsAccepted += 1;
			TotAcceptedLen += BuffOfs;
			}
		}
	}
if(pSeqBuff != NULL)
	free(pSeqBuff);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"ProcessFastaFile - %d parsed, %d accepted, %dbp mean length, %d sequences not accepted for indexing as length under %dbp ",
					SeqID,NumSeqsAccepted,NumSeqsAccepted == 0 ? 0 :(int)(TotAcceptedLen/NumSeqsAccepted),NumSeqsUnderlength,MinSeqLen);
return(Rslt);
}


int 
CPBAssemb::LoadTargetSeqs(int MinSeqLen,int NumTargFiles,char **pszTargFiles)		// parse, and index sequences in these files into in memory suffix array; file expected to contain either fasta or fastq sequences
{
int Rslt;
int Idx;
int NumGlobs;
int64_t SumFileSizes;

CSimpleGlob glob(SG_GLOB_FULLSORT);

	// determine crude estimate of total genome size from the sum of file sizes
for(Idx = 0; Idx < NumTargFiles; Idx++)
	if (glob.Add(pszTargFiles[Idx]) < SG_SUCCESS)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to glob '%s",pszTargFiles[Idx]);
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
	Reset(false);
	return(eBSFerrOpnFile);
	}

Rslt = eBSFSuccess;

for (int n = 0; Rslt >= eBSFSuccess &&  n < glob.FileCount(); ++n)
	{
		// try opening as a fasta file
	Rslt = ProcessFastaFile(cMinSeqLen,glob.File(n));
	if(Rslt < eBSFSuccess)
		{
		Reset(false);
		return(Rslt);
		}

	}
if(m_pSeqStore->GetNumSeqs() < 3)		// need at least 3 sequences if attempting to error correct
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"CreateBioseqSuffixFile: Insufficient ( < 3) sequences accepted for processing");
	Reset(false);
	return(eBSFerrNoEntries);
	}
return(Rslt);

}

int
CPBAssemb::Process(etPBPMode PMode,		// processing mode
		int MinScaffSeqLen,			// individual scaffold sequences must be of at least this length (defaults to 5Kbp)
		int MinScaffOverlap,		// pairs of targeted scaffold sequences must overlap by at least this many bp to be considered for merging into a longer scaffold sequence (defaults to 5Kbp) 
	    int Min1kScore,             // minimum normalised 1Kbp overlap score
	    bool bAcceptOrphanSeqs,		// also report sequences which are not overlapped or overlapping any other sequence
		bool bSenseOnlyOvlps,				// if false then both sense/sense and sense/antisense overlaps will be accepted and processed, otherwise sense/sense only overlaps accepted and processed
		char *pszMAFFile,			// pregenerated multialignment sequence overlap loci details
		int NumErrCorrectedFiles,	// number of error corrected sequence specs
		char *pszErrCorrectedFiles[],		// input error corrected sequence files
	    char *pszOutFile,			// where to write merged scaffolded sequences
		int NumThreads)				// maximum number of worker threads to use
{
int Rslt = eBSFSuccess;
int Idx;
uint32_t NumTargSeqs;
uint32_t CurNodeID;
uint32_t MaxSeqLen;
tsPBAScaffNode *pCurPBScaffNode;

Reset(false);
CreateMutexes();

m_PMode = PMode;

m_MinScaffSeqLen = MinScaffSeqLen;	
m_MinScaffOverlap = MinScaffOverlap; 
m_MinScaffScoreThres = Min1kScore;
m_bAcceptOrphanSeqs = bAcceptOrphanSeqs;
m_bSenseOnlyOvlps = bSenseOnlyOvlps;

memset(m_szErrCorrectedFiles,0,sizeof(m_szErrCorrectedFiles));
if(PMode == ePBPMScaffold)
	{
	m_NumErrCorrectedFiles = NumErrCorrectedFiles;
	for(Idx = 0; Idx < NumErrCorrectedFiles; Idx++)
		strcpy(m_szErrCorrectedFiles[Idx],pszErrCorrectedFiles[Idx]);
	}
else
	m_NumErrCorrectedFiles = 0;

strncpy(m_szOutFile,pszOutFile,sizeof(m_szOutFile));
m_szOutFile[sizeof(m_szOutFile)-1] = '\0';	

m_NumThreads = NumThreads;	
if(m_pSeqStore != NULL)
	{
	delete m_pSeqStore;
	m_pSeqStore = NULL;
	}

if((m_pSeqStore = new CSeqStore) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Unable to instantiate instance of CSeqStore");
	return(eBSFerrObj);
	}
m_pSeqStore->Reset();

if((Rslt = LoadTargetSeqs(MinScaffSeqLen,NumErrCorrectedFiles,pszErrCorrectedFiles)) < eBSFSuccess)
	{
	Reset(false);
	return(Rslt);
	}
// get number of sequences loaded as targets to be used for scaffolding
NumTargSeqs = m_pSeqStore->GetNumSeqs();
if(NumTargSeqs < 1)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"No targeted scaffold sequences in file(s)");
	Reset(false);
	return(eBSFerrNoEntries);
	}

if((Rslt=m_pSeqStore->GenSeqDescrIdx())!=eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to generate index over sequence descriptors in file(s)");
	Reset(false);
	return(Rslt);
	}

if((uint32_t)(Rslt = (int)LoadPacBioOvlps(pszMAFFile, true, bSenseOnlyOvlps)) > cMaxValidID)
	{
	Reset(false);
	return(Rslt);
	}
if(Rslt == 0)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadPacBioOvlps: Nothing to do, no accepted overlaps to process");
	if(m_NumRejectedMinScaffOverlap > 0 || m_NumRejectedScoreThres > 0 || m_NumRejectContained > 0 || m_NumRejectArtefact > 0 || m_NumRejectAntisense > 0)
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Rejected %u min seq len,  %u overlap len, %u score, %u contained, %u antisense, %u artifact) from file '%s'",
											m_NumRejectedMinSeqLen,m_NumRejectedMinScaffOverlap, m_NumRejectedScoreThres,m_NumRejectContained,m_NumRejectAntisense,m_NumRejectArtefact,pszMAFFile);
	Reset(false);
	return(Rslt);
	}

if(m_pAssembGraph != NULL)
	{
	delete m_pAssembGraph;
	m_pAssembGraph = NULL;
	}

if((m_pAssembGraph = new CAssembGraph) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Unable to instantiate instance of CAssembGraph");
	Reset(false);
	return(eBSFerrObj);
	}
if((Rslt=m_pAssembGraph->Init(m_bSenseOnlyOvlps, m_MinScaffScoreThres, m_bAcceptOrphanSeqs,min(4,NumThreads)))!=eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Unable to initialise CAssembGraph");
	Reset(false);
	return(eBSFerrObj);
	}

if((m_pPBScaffNodes = new tsPBAScaffNode [NumTargSeqs + 1]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate for %d scaffold nodes",NumTargSeqs);
	Reset(false);
	return(eBSFerrMem);
	}
memset(m_pPBScaffNodes,0,sizeof(tsPBAScaffNode) * (NumTargSeqs+1));
m_AllocdPBScaffNodes = NumTargSeqs;
if((m_pMapEntryID2NodeIDs = new uint32_t [NumTargSeqs + 1]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate for %d scaffold nodes",NumTargSeqs);
	Reset(false);
	return(eBSFerrMem);
	}
memset(m_pMapEntryID2NodeIDs,0,sizeof(uint32_t) * (NumTargSeqs+1));
m_NumPBScaffNodes = 0;

// initialise scaffold nodes
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Initialising %d scaffold nodes and graph vertices",NumTargSeqs);
MaxSeqLen = 0;
pCurPBScaffNode = m_pPBScaffNodes;
for(CurNodeID = 1; CurNodeID <= NumTargSeqs; CurNodeID++,pCurPBScaffNode++)
	{
	pCurPBScaffNode->SeqLen = m_pSeqStore->GetLen(CurNodeID);
	pCurPBScaffNode->EntryID = CurNodeID;
	pCurPBScaffNode->flgUnderlength = pCurPBScaffNode->SeqLen < m_MinScaffSeqLen ? 1 : 0;
	if(MaxSeqLen == 0 || pCurPBScaffNode->SeqLen > (uint32_t)MaxSeqLen)
		MaxSeqLen = pCurPBScaffNode->SeqLen;

	if((pCurPBScaffNode->VertexID = m_pAssembGraph->AddVertex(pCurPBScaffNode->SeqLen,CurNodeID)) > cMaxValidID)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddVertex failed");
		Reset(false);
		return((int)pCurPBScaffNode->VertexID);
		}
	}

m_NumPBScaffNodes = NumTargSeqs;
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Finalising %d graph vertices",NumTargSeqs);
m_pAssembGraph->FinaliseVertices();
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Finalised graph vertices");

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Sorting %d scaffold nodes",NumTargSeqs);

if(m_NumPBScaffNodes > 1)
	{
	// sort scaffold nodes by sequence length descending
	m_mtqsort.SetMaxThreads(NumThreads);
	m_mtqsort.qsort(m_pPBScaffNodes,m_NumPBScaffNodes,sizeof(tsPBAScaffNode),SortLenDescending);
	}

pCurPBScaffNode = m_pPBScaffNodes;
for(CurNodeID = 1; CurNodeID <= NumTargSeqs; CurNodeID++,pCurPBScaffNode++)
	{
	pCurPBScaffNode->NodeID = CurNodeID;
	m_pMapEntryID2NodeIDs[pCurPBScaffNode->EntryID-1] = CurNodeID;
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading previously generated overlap detail from file '%s' ...",pszMAFFile);
if((uint32_t)(Rslt = LoadPacBioOvlps(pszMAFFile, false,bSenseOnlyOvlps)) > cMaxValidID)
	return(Rslt);
if(m_NumAcceptedOverlaps == 0)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"No sequence overlaps loaded");
	return(Rslt);										
	}

if((uint32_t)(Rslt = m_pAssembGraph->FinaliseEdges())  > cMaxValidID)
	{
	Reset(false);
	return(Rslt);
	}
m_NumAcceptedOverlaps = Rslt;
if(m_NumAcceptedOverlaps == 0)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"FinaliseEdges: No sequence overlaps");
	return(Rslt);										
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Accepted %u sequence overlap details from file '%s' plus inferred overlaps ...",m_NumAcceptedOverlaps,pszMAFFile);


uint32_t NumDiscComponents = 0;
if((NumDiscComponents = m_pAssembGraph->IdentifyDiscComponents()) > cMaxValidID)
	{
	Reset(false);
	return((int)NumDiscComponents);
	}

switch(PMode) {
	case ePBPMScaffold:
		if((Rslt = m_pAssembGraph->FindHighestScoringPaths()) >= eBSFSuccess)
			Rslt = m_pAssembGraph->WriteContigSeqs(pszOutFile,m_pSeqStore);
		break;
	case ePBPMToGEFX:
		Rslt = m_pAssembGraph->ReportVerticesEdgesGEXF(pszOutFile,m_pSeqStore);
		break;
	case ePBPMToGraphML:
		Rslt = m_pAssembGraph->ReportVerticesEdgesGraphML(pszOutFile,m_pSeqStore);
		break;
	}

Reset(false);
return(Rslt);
}

// MapEntryID2NodeID
// Given a suffix array entry identifier returns the corresponding node identifier
uint32_t										// returned tsPBScaffNode node identifier
CPBAssemb::MapEntryID2NodeID(uint32_t EntryID)			// suffix array entry identifier
{
if(EntryID == 0 || EntryID > m_NumPBScaffNodes || m_pMapEntryID2NodeIDs == NULL)
	return(0);
return(m_pMapEntryID2NodeIDs[EntryID-1]);
}




// SortLenDescending
// Sort scaffolding nodes by length descending with entry identifiers as tie breaker
int
CPBAssemb::SortLenDescending(const void *arg1, const void *arg2)
{
tsPBAScaffNode *pEl1 = (tsPBAScaffNode *)arg1;
tsPBAScaffNode *pEl2 = (tsPBAScaffNode *)arg2;

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


