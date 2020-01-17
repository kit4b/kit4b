/*
This toolkit is a source base clone of 'PacBioKanga' release 4.4.2 (https://github.com/csiro-crop-informatics/biokanga) and contains
significant source code changes enabling new functionality and resulting process parameterisation changes. These changes have resulted in
incompatibilty with 'PacBioKanga'.

Because of the potentential for confusion by users unaware of functionality and process parameterisation changes then the modified source base
and resultant compiled executables have been renamed to 'kit4pacbio' - PacBio K-mer Informed Toolkit for Bioinformatics.
The renaming will force users of the 'PacBioKanga' toolkit to examine scripting which is dependent on existing 'PacBioKanga'
parameterisations so as to make appropriate changes if wishing to utilise 'kit4b' parameterisations and functionality.

'kit4pacbio' is being released under the Opensource Software License Agreement (GPLv3)
'kit4pacbio' is Copyright (c) 2019, Dr Stuart Stephen
Please contact Dr Stuart Stephen < stuartjs@g3bio.com > if you have any questions regarding 'kit4b'.

Orginal 'BioKanga' copyright notice has been retained and is as follows.
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
#include "SWAlign.h"
#include "../libkit4b/bgzf.h"
#include "SeqStore.h"
#include "SSW.h"
#include "PBFilter.h"


int
ProcPacBioFilter(etPBPMode PMode,	// processing mode
		int SMRTBellSensitivity,	// sensitivity of SMRTBell detection - 5: max, 1: min
		int Trim5,					// 5' trim accepted reads by this many bp
		int Trim3,					// 3' trim accepted reads by this many bp
		int MinReadLen,				// read sequences must be at least this length after any end timming
		int ContamARate,			// PacBio expected accuracy event rate, used when contaminate processing
		int ContamOvlpLen,			// minimum contaminate overlap length to check for
		char *pszContamFile,		// name of file containing contaminate sequences
		int NumInputFiles,			// number of input files
		char *pszInputFiles[],		// names (wildcards allowed) of input files containing reads to be filtered
		char *pszOutFile,			// name of file in which to write filter accepted and trimmed sequences
		int NumThreads);			// maximum number of worker threads to use


#ifdef _WIN32
int ProcFilter(int argc, char* argv[])
{
// determine my process name
_splitpath(argv[0],NULL,NULL,gszProcName,NULL);
#else
int
ProcFilter(int argc, char** argv)
{
// determine my process name
CUtility::splitpath((char *)argv[0],NULL,gszProcName);
#endif

int iFileLogLevel;			// level of file diagnostics
int iScreenLogLevel;		// level of file diagnostics
char szLogFile[_MAX_PATH];	// write diagnostics to this file
int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure

int Idx;
int PMode;					// processing mode
int SMRTBellSensitivity;	// sensitivity of SMRTBell detection - 5: max, 1: min
int Trim5;					// 5' trim accepted reads by this many bp
int Trim3;					// 3' trim accepted reads by this many bp
int MinReadLen;				// read sequences must be at least this length after any end trimming

int ContamARate;			// PacBio sequences accuracy rate, used when contaminate processing
int ContamOvlpLen;				// minimum contaminate overlap length to check for


int NumInputFiles;			// number of input files
char *pszInputFiles[cMaxInfiles];	// names (wildcards allowed) of input files containing reads to be filtered
char szOutFile[_MAX_PATH];	// name of file in which to write filter accepted and trimmed sequences
char szContamFile[_MAX_PATH];	// name of file containing contaminate sequences

int NumberOfProcessors;		// number of installed CPUs
int NumThreads;				// number of threads (0 defaults to number of CPUs)

char szSQLiteDatabase[_MAX_PATH];	// results summaries to this SQLite file
char szExperimentName[cMaxDatasetSpeciesChrom+1];			// experiment name
char szExperimentDescr[1000];		// describes experiment

struct arg_lit  *help    = arg_lit0("h","help",                 "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *pmode = arg_int0("m","pmode","<int>",			"processing mode - 0 SMRTBell trimming, 1 remove contaminate containing sequences");
struct arg_int *smrtbellsensitivity = arg_int0("s","smrtbellsensitivity","<int>",			"sensitivity of SMRTBell detection (default 3, range 1: min to 5: max)");

struct arg_int *trim5 = arg_int0("z","trim5","<int>",			"5' trim accepted reads by this many bp (default 0, range 0 to 10000)");
struct arg_int *trim3 = arg_int0("Z","trim3","<int>",			"3' trim accepted reads by this many bp (default 0, range 0 to 10000)");
struct arg_int *minreadlen = arg_int0("l","minreadlen","<int>",		"read sequences must be at least this length after any end trimming (default 1000, range 500 to 50000)");

struct arg_file *contamfile = arg_file0("I","contam","<file>",		"file containing contaminate sequences");
struct arg_int *contamarate = arg_int0("c","contamerate","<int>",	"PacBio sequences expected accuracy rate (default 80, range 70 to 99");
struct arg_int *contamovlplen = arg_int0("C","contamovlplen","<int>",	"Minimum contaminate overlap length (default 500, range 250 to 5000");

struct arg_file *inputfiles = arg_filen("i","in","<file>", 1, cMaxInfiles,		"input file(s) containing PacBio long reads to be filtered");
struct arg_file *outfile = arg_file1("o","out","<file>",			"output accepted filtered reads to this file");

struct arg_int *threads = arg_int0("T","threads","<int>",		"number of processing threads 0..128 (defaults to 0 which sets threads to number of CPU cores)");

struct arg_file *summrslts = arg_file0("q","sumrslts","<file>",				"Output results summary to this SQLite3 database file");
struct arg_str *experimentname = arg_str0("w","experimentname","<str>",		"experiment name SQLite3 database file");
struct arg_str *experimentdescr = arg_str0("W","experimentdescr","<str>",	"experiment description SQLite3 database file");

struct arg_end *end = arg_end(200);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					pmode,smrtbellsensitivity,trim5,trim3,minreadlen,summrslts,experimentname,experimentdescr,
					contamarate,contamovlplen,contamfile,	inputfiles,outfile,threads,
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
	szContamFile[0] = '\0';
	ContamARate = 0;
	ContamOvlpLen = 0;

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

	PMode = (etPBPMode)(pmode->count ? pmode->ival[0] : (int)ePBPMFilter);
	if(PMode < ePBPMFilter || PMode > ePBPMContam)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Processing mode '-m%d' must be in range 0..%d",PMode,ePBPMContam);
		return(1);
		}

	if(PMode == ePBPMFilter)
		{
		SMRTBellSensitivity = smrtbellsensitivity->count ? smrtbellsensitivity->ival[0] : 3;
		if(SMRTBellSensitivity < 1 || SMRTBellSensitivity > 5)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: SMRTBell detection sensitivity threshold '-s%d' must be in range 1..5",SMRTBellSensitivity);
			return(1);
			}
		}
	else
		SMRTBellSensitivity = 0;

	Trim5 = trim5->count ? trim5->ival[0] : cDfltTrim5;
	if(Trim5 < 0 || Trim5 > 10000)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: 5' end trim '-z%d' must be in range 0..10000",Trim5);
		return(1);
		}

	Trim3 = trim3->count ? trim3->ival[0] : cDfltTrim3;
	if(Trim3 < 0 || Trim3 > 10000)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: 3' end trim '-Z%d' must be in range 0..10000",Trim3);
		return(1);
		}

	MinReadLen = minreadlen->count ? minreadlen->ival[0] : cDfltMinReadLen;
	if(MinReadLen < 500 || MinReadLen > 50000)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: minimum read length after any trim '-l%d' must be in range 500..50000",MinReadLen);
		return(1);
		}

	if(PMode == ePBPMContam)
		{
		ContamARate = contamarate->count ? contamarate->ival[0] : 80;
		if(ContamARate < 70 || ContamARate > 100)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: PacBio sequences accuracy rate '-c%d' must be in range 70..99",ContamARate);
			return(1);
			}
		if(ContamARate == 100)	// whilst accepting user input of 100% accuracy internally allowing the odd error event!
			ContamARate = 99;
		ContamOvlpLen = contamovlplen->count ? contamovlplen->ival[0] : 500;
		if(ContamOvlpLen < 250 || ContamOvlpLen > 5000)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: minimum read length after any trim '-C%d' must be in range 250..5000",ContamOvlpLen);
			return(1);
			}

		if(contamfile->count)
			{
			strncpy(szContamFile,contamfile->filename[0],_MAX_PATH);
			szContamFile[_MAX_PATH-1] = '\0';
			CUtility::TrimQuotedWhitespcExtd(szContamFile);
			if(strlen(szContamFile) == 0)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace no contaminate sequences file specified with '-I<contamfile>'");
				return(1);
				}
			}
		else
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Contaminate filtering requested but no contaminate sequences file specified with '-I<contamfile>'");
			return(1);
			}
		}

	if (!inputfiles->count)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: No input file(s) specified with '-i<filespec>' option)");
		return(1);
		}

	for (NumInputFiles = Idx = 0; NumInputFiles < cMaxInfiles && Idx < inputfiles->count; Idx++)
		{
		pszInputFiles[Idx] = NULL;
		if (pszInputFiles[NumInputFiles] == NULL)
			pszInputFiles[NumInputFiles] = new char[_MAX_PATH];
		strncpy(pszInputFiles[NumInputFiles], inputfiles->filename[Idx], _MAX_PATH);
		pszInputFiles[NumInputFiles][_MAX_PATH - 1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(pszInputFiles[NumInputFiles]);
		if (pszInputFiles[NumInputFiles][0] != '\0')
			NumInputFiles++;
		}

	if (!NumInputFiles)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: After removal of whitespace, no input file(s) specified with '-i<filespec>' option");
		return(1);
		}

	strcpy(szOutFile,outfile->filename[0]);
	szOutFile[_MAX_PATH-1] = '\0';
	CUtility::TrimQuotedWhitespcExtd(szOutFile);
	if(strlen(szOutFile) == 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace no output file specified with '-o<outfile>'");
		return(1);
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
		case ePBPMFilter:									// filter PacBio reads which have retained hairpins
			pszMode = (char *)"Trim/split PacBio reads containing SMRTBell adapter sequences";
			break;
		case ePBPMContam:									// filter reads containing contaminate sequences
			pszMode = (char *)"Trim/split PacBio reads containing contaminate sequences";
			break;

	}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"processing mode: '%s'",pszMode);
	if(PMode == ePBPMFilter)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"SMRTBell detection sensitivity: %d",SMRTBellSensitivity);	
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"5' trim accepted reads by: %dbp",Trim5);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"3' trim accepted reads by: %dbp",Trim3);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"read sequences must be at least this length after any end trimming: %dbp",MinReadLen);

	
	if(PMode == ePBPMContam)
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"PacBio sequences minimum accuracy rate: %d",ContamARate);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum contaminate overlap length: %d",ContamOvlpLen);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"processing for contaminant sequences in file: '%s'",szContamFile);
		}

	for (Idx = 0; Idx < NumInputFiles; Idx++)
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "Input reads file (%d) : '%s'", Idx + 1, pszInputFiles[Idx]);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"name of file in which to write filter accepted and trimmed sequences: '%s'",szOutFile);

	if(szExperimentName[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"This processing reference: %s",szExperimentName);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"number of threads : %d",NumThreads);

	if(gExperimentID > 0)
		{
		int ParamID;
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szLogFile),"log",szLogFile);

		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(PMode),"pmode",&PMode);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(SMRTBellSensitivity),"smrtbellsensitivity",&SMRTBellSensitivity);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(Trim5),"trim5",&Trim5);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(Trim3),"trim3",&Trim3);
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(MinReadLen),"minreadlen",&MinReadLen);
		if(PMode == ePBPMContam)
			{
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTText,(int)strlen(szContamFile),"contamfile",szContamFile);
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(ContamARate),"contamarate",&ContamARate);
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID,ePTInt32,(int)sizeof(ContamOvlpLen),"contamovlplen",&ContamOvlpLen);
			}
		ParamID = gSQLiteSummaries.AddParameter(gProcessingID, ePTInt32, sizeof(NumInputFiles), "NumInputFiles", &NumInputFiles);
		for (Idx = 0; Idx < NumInputFiles; Idx++)
			ParamID = gSQLiteSummaries.AddParameter(gProcessingID, ePTText, (int)strlen(pszInputFiles[Idx]), "inpe1", pszInputFiles[Idx]);

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
	Rslt = ProcPacBioFilter((etPBPMode)PMode,SMRTBellSensitivity,Trim5,Trim3,MinReadLen,ContamARate,ContamOvlpLen,
												szContamFile,NumInputFiles,pszInputFiles,szOutFile,NumThreads);
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
ProcPacBioFilter(etPBPMode PMode,	// processing mode
		int SMRTBellSensitivity, // sensitivity of SMRTBell detection - 5: max, 1: min
		int Trim5,					// 5' trim accepted reads by this many bp
		int Trim3,					// 3' trim accepted reads by this many bp
		int MinReadLen,				// read sequences must be at least this length after any end timming
		int ContamARate,			// PacBio expected accuracy event rate, used when contaminate processing
		int ContamOvlpLen,			// minimum contaminate overlap length to check for
		char *pszContamFile,		// name of file containing contaminate sequences
		int NumInputFiles,			// number of input files
		char *pszInputFiles[],		// names (wildcards allowed) of input files containing reads to be filtered
		char *pszOutFile,			// name of file in which to write filter accepted and trimmed sequences
		int NumThreads)				// maximum number of worker threads to use
{
int Rslt;
CPBFilter *pPacBioer;

if((pPacBioer = new CPBFilter)==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: Unable to instantiate CPBFilter");
	return(eBSFerrObj);
	}
Rslt = pPacBioer->Process(PMode,SMRTBellSensitivity,Trim5,Trim3,MinReadLen,ContamARate,ContamOvlpLen,pszContamFile,NumInputFiles,pszInputFiles,pszOutFile,NumThreads);
delete pPacBioer;
return(Rslt);
}

CPBFilter::CPBFilter() // relies on base classes constructors
{
m_bMutexesCreated = false;
m_hOutFile = -1;
m_pszOutBuff = NULL;
m_pSWAlign = NULL;
Init();
}

CPBFilter::~CPBFilter() // relies on base classes destructors
{
Reset();
}


void
CPBFilter::Init(void)
{
if(m_hOutFile != -1)
	{
#ifdef _WIN32
	_commit(m_hOutFile);
#else
	fsync(m_hOutFile);
#endif
	close(m_hOutFile);
	m_hOutFile = -1;
	}
if(m_pszOutBuff != NULL)
	{
	delete m_pszOutBuff;
	m_pszOutBuff = NULL;
	}

if(m_pSWAlign != NULL)
	{
	delete m_pSWAlign;
	m_pSWAlign = NULL;
	}

m_PMode = ePBPMFilter;

m_Trim5 = cDfltTrim5;							
m_Trim3 = cDfltTrim3;							
m_MinReadLen = cDfltMinReadLen;		

m_TotProcessed = 0;
m_TotAccepted = 0;
m_TotSMRTBells = 0;
m_TotSMRTBellSubSeqs = 0;
m_TotContamTrimmed = 0;
m_TotRejected = 0;
m_TotContamRejected = 0;
m_TotUnderLen = 0;
m_OutBuffIdx = 0;
m_AllocOutBuffSize=0;
m_NumInputFiles = 0;
m_ppszInputFiles = NULL;			
m_szOutFile[0] = '\0';			
m_szContamFile[0] = '0';

m_NumThreads = 0;
if(m_bMutexesCreated)
	DeleteMutexes();
m_bMutexesCreated = false; 
m_PacBioUtility.Reset();
}

void
CPBFilter::Reset(void)
{
Init();
}


int
CPBFilter::CreateMutexes(void)
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

m_bMutexesCreated = true;
return(eBSFSuccess);
}

void
CPBFilter::DeleteMutexes(void)
{
if(!m_bMutexesCreated)
	return;
#ifdef _WIN32
CloseHandle(m_hMtxIterReads);
#else
pthread_mutex_destroy(&m_hMtxIterReads);
pthread_rwlock_destroy(&m_hRwLock);
#endif
m_bMutexesCreated = false;
}

void
CPBFilter::AcquireSerialise(void)
{
#ifdef _WIN32
WaitForSingleObject(m_hMtxIterReads,INFINITE);
#else
pthread_mutex_lock(&m_hMtxIterReads);
#endif
}

void
CPBFilter::ReleaseSerialise(void)
{
#ifdef _WIN32
ReleaseMutex(m_hMtxIterReads);
#else
pthread_mutex_unlock(&m_hMtxIterReads);
#endif
}

void
CPBFilter::AcquireLock(bool bExclusive)
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
CPBFilter::ReleaseLock(bool bExclusive)
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
CPBFilter::Process(etPBPMode PMode,	// processing mode
		int SMRTBellSensitivity, // sensitivity of SMRTBell detection - 5: max, 1: min
		int Trim5,					// 5' trim accepted reads by this many bp
		int Trim3,					// 3' trim accepted reads by this many bp
		int MinReadLen,				// read sequences must be at least this length after any end timming
		int ContamARate,			// PacBio expected accuracy event rate, used when contaminate processing
		int ContamOvlpLen,			// minimum contaminate overlap length to check for
		char *pszContamFile,		// name of file containing contaminate sequences
		int NumInputFiles,			// number of input files
		char *pszInputFiles[],		// names (wildcards allowed) of input files containing reads to be filtered
		char *pszOutFile,			// name of file in which to write filter accepted and trimmed sequences
		int NumThreads)			// maximum number of worker threads to use
{
int Rslt = eBSFSuccess;

Reset();
CreateMutexes();

m_PMode = PMode;
m_Trim5 = Trim5;					
m_Trim3 = Trim3;					
m_MinReadLen = MinReadLen;	
m_NumInputFiles = NumInputFiles;
m_ppszInputFiles = 	pszInputFiles;
strncpy(m_szOutFile,pszOutFile,sizeof(m_szOutFile));
m_szOutFile[sizeof(m_szOutFile)-1] = '\0';	

m_NumThreads = NumThreads;	
m_PacBioUtility.Reset();

if(pszContamFile == NULL || pszContamFile[0] == '\0')
	{
	m_ContamARate = 0;
	m_ContamOvlpLen = 0;
	m_szContamFile[0] = '\0';
	}
else
	{
	strcpy(m_szContamFile,pszContamFile);
	m_ContamARate = ContamARate;
	m_ContamOvlpLen = ContamOvlpLen;
	if((m_pSWAlign = new CSWAlign) == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate instance of CSWAlign");
		Reset();
		return(eBSFerrInternal);
		}

	UINT32 DeltaCoreOfs;				// ranges between 1 and 4
	UINT32 CoreSeqLen;                  // ranges between 10 and 25
	UINT32 MinNumCores;                 // ranges between 5 and 50 
	int MaxInitiatePathOfs;	// require SW paths to have started within this many bp (0 to disable) on either the probe or target - effectively anchoring the SW 
	int OverlapFloat;		// with this overlap float

	int MatchScore = 3;                 // PacBio has high proportion of InDels relative to simple base subs         
	int MismatchPenalty = -7;    
	int GapOpenPenalty = -4;     
	int GapExtnPenalty = -1;     
	int DlyGapExtn = 2;   

	if(ContamARate <= 80)		// 70..80
		{
		DeltaCoreOfs = 1;
		CoreSeqLen = 10;
		MinNumCores = 5;
		MaxInitiatePathOfs = 1000;
		OverlapFloat = 600;
		}
	else
		{
		if(ContamARate <= 85)    // 81..85 
			{
			DeltaCoreOfs = 2;
			CoreSeqLen = 10;
			MinNumCores = 5 + ContamARate - 81;
			MaxInitiatePathOfs = 800;
			OverlapFloat = 500;
			}
		else
			{
			if(ContamARate < 90) // 86..90
				{
				DeltaCoreOfs = 2;
				CoreSeqLen = 11;
				MinNumCores = 10 + 2 * (ContamARate - 86);
				MaxInitiatePathOfs = 600;
				OverlapFloat = 400;
				}
			else
				{
				if(ContamARate <= 95) // 91..95
					{
					MatchScore = 2;
					MismatchPenalty = -4;
					GapOpenPenalty = -3;
					GapExtnPenalty = -1;
					DlyGapExtn = 1;
					DeltaCoreOfs = 3;
					CoreSeqLen = 12 + (ContamARate - 91);
					MinNumCores = 20 + 2 * (ContamARate - 91);
					MaxInitiatePathOfs = 400;
					OverlapFloat = 300;
					}
				else                 // 96..100
					{
					MatchScore = 1;
					MismatchPenalty = -2;
					GapOpenPenalty = -3;
					GapExtnPenalty = -1;
					DlyGapExtn = 1;
					DeltaCoreOfs = 4;
					CoreSeqLen = 15 + (ContamARate - 96);
					MinNumCores = 36 + 4 * (ContamARate - 96);
					MaxInitiatePathOfs = 200;
					OverlapFloat = 150;
					}
				}
			}
		}      

	if((m_pSWAlign->Initialise(NumThreads,m_ContamOvlpLen,cMaxSeedCoreDepth,DeltaCoreOfs,CoreSeqLen,MinNumCores,cMaxAcceptHitsPerSeedCore,cDfltMaxProbeSeqLen)) != eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to initialise m_pSWAlign");
		Reset();
		return(eBSFerrOpnFile);
		}
	m_pSWAlign->SetScores(MatchScore, MismatchPenalty, GapOpenPenalty, GapExtnPenalty, DlyGapExtn);
	m_pSWAlign->SetCPScores();
	m_pSWAlign->SetMaxInitiatePathOfs(MaxInitiatePathOfs, OverlapFloat);

	if((m_pSWAlign->LoadTargetSeqs(cMinContamSeqLen, pszContamFile)) <= 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to load any contaminate sequences from file: '%s'",pszContamFile);
		Reset();
		return(eBSFerrOpnFile);
		}

	}
	

if((Rslt = m_PacBioUtility.StartAsyncLoadSeqs(NumInputFiles,pszInputFiles)) < eBSFSuccess)
	{
	Reset();
	return(Rslt);
	}

Rslt = ProcessFiltering(SMRTBellSensitivity,cDfltMaxPacBioSeqLen,NumThreads);

Reset();
return(Rslt);
}

#ifdef _WIN32
unsigned __stdcall PBFilterThread(void * pThreadPars)
#else
void *PBFilterThread(void * pThreadPars)
#endif
{
int Rslt;
tsThreadPBFilter *pPars = (tsThreadPBFilter *)pThreadPars;			// makes it easier not having to deal with casts!
CPBFilter *pPBFilter = (CPBFilter *)pPars->pThis;

Rslt = pPBFilter->PBFilterReads(pPars);
pPars->Rslt = Rslt;
#ifdef _WIN32
_endthreadex(0);
return(eBSFSuccess);
#else
pthread_exit(NULL);
#endif
}

int
CPBFilter::ProcessFiltering(int SMRTBellSensitivity, // sensitivity of SMRTBell detection - 5: max, 1: min
							int MaxSeqLen,			// max length sequence expected
							int NumOvlpThreads)	// filtering using at most this many threads
{
tsThreadPBFilter *pThreadPutOvlps;
int ThreadIdx;
tsThreadPBFilter *pThreadPar;

pThreadPutOvlps = new tsThreadPBFilter [NumOvlpThreads];

pThreadPar = pThreadPutOvlps;
for(ThreadIdx = 0; ThreadIdx < NumOvlpThreads; ThreadIdx++,pThreadPar++)
	memset(pThreadPar,0,sizeof(tsThreadPBFilter));

if(ThreadIdx != NumOvlpThreads)	// any errors whilst allocating memory?
	{
	do {
		if(pThreadPar->pSW != NULL)
			delete pThreadPar->pSW;
		pThreadPar -= 1;
		ThreadIdx -= 1;
		}
	while(ThreadIdx >= 0);
	delete pThreadPutOvlps;
	Reset();
	return((INT64)eBSFerrMem);
	}


if((m_pszOutBuff = new char [cAllocOutBuffSize]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate %d chars for buffering",cAllocOutBuffSize);
	Reset();
	return(eBSFerrMem);
	}
m_AllocOutBuffSize = cAllocOutBuffSize;

#ifdef _WIN32
m_hOutFile = open(m_szOutFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE));
#else
if((m_hOutFile = open(m_szOutFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE))!=-1)
	if(ftruncate(m_hOutFile,0)!=0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",m_szOutFile,strerror(errno));
		Reset();
		return(eBSFerrCreateFile);
		}
#endif
if(m_hOutFile < 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output file '%s'",m_szOutFile);
	Reset();
	return(eBSFerrCreateFile);
	}

pThreadPar = pThreadPutOvlps;
for (ThreadIdx = 1; ThreadIdx <= NumOvlpThreads; ThreadIdx++, pThreadPar++)
	{
	pThreadPar->ThreadIdx = ThreadIdx;
	pThreadPar->pThis = this;
	pThreadPar->SMRTBellSensitivity = SMRTBellSensitivity;
#ifdef _WIN32
	pThreadPar->threadHandle = (HANDLE)_beginthreadex(NULL, 0x0fffff, PBFilterThread, pThreadPar, 0, &pThreadPar->threadID);
#else
	pThreadPar->threadRslt = pthread_create(&pThreadPar->threadID, NULL, PBFilterThread, pThreadPar);
#endif
	}

// allow threads a few seconds to startup
#ifdef _WIN32
Sleep(10000);
#else
sleep(10);
#endif

AcquireLock(false);
if(m_PMode == ePBPMFilter)
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Processing: %d - accepted %d, rejected: %d, under length %d, SMRTBells identified: %d, SMRTBell subseqs: %d", m_TotProcessed, m_TotAccepted, m_TotRejected, m_TotUnderLen, m_TotSMRTBells, m_TotSMRTBellSubSeqs);
else
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Processing: %d - accepted: %d, rejected: %d, under length: %d, vector trimmed: %d", m_TotProcessed, m_TotAccepted, m_TotRejected, m_TotUnderLen, m_TotContamTrimmed);
ReleaseLock(false);

pThreadPar = pThreadPutOvlps;
for (ThreadIdx = 0; ThreadIdx < NumOvlpThreads; ThreadIdx++, pThreadPar++)
	{
#ifdef _WIN32
	while (WAIT_TIMEOUT == WaitForSingleObject(pThreadPar->threadHandle, 60000))
		{
		AcquireLock(false);
		if(m_PMode == ePBPMFilter)
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "Processing: %d - accepted %d, rejected: %d, under length %d, SMRTBells identified: %d, SMRTBell subseqs: %d", m_TotProcessed, m_TotAccepted, m_TotRejected, m_TotUnderLen, m_TotSMRTBells, m_TotSMRTBellSubSeqs);
		else
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "Processing: %d - accepted: %d, rejected: %d, under length: %d, vector trimmed: %d", m_TotProcessed, m_TotAccepted, m_TotRejected, m_TotUnderLen, m_TotContamTrimmed);
		ReleaseLock(false);
		};
	CloseHandle(pThreadPar->threadHandle);
#else
	struct timespec ts;
	int JoinRlt;
	clock_gettime(CLOCK_REALTIME, &ts);
	ts.tv_sec += 60;
	while ((JoinRlt = pthread_timedjoin_np(pThreadPar->threadID, NULL, &ts)) != 0)
		{
		AcquireLock(false);
		if(m_PMode == ePBPMFilter)
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "Processing: %d - accepted %d, rejected: %d, under length %d, SMRTBells identified: %d, SMRTBell subseqs: %d", m_TotProcessed, m_TotAccepted, m_TotRejected, m_TotUnderLen, m_TotSMRTBells, m_TotSMRTBellSubSeqs);
		else
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "Processing: %d - accepted: %d, rejected: %d, under length: %d, vector trimmed: %d", m_TotProcessed, m_TotAccepted, m_TotRejected, m_TotUnderLen, m_TotContamTrimmed);
		ReleaseLock(false);
		ts.tv_sec += 60;
		}
#endif
	}

if(m_PMode == ePBPMFilter)
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Processed: %d - accepted %d, rejected: %d, under length %d, SMRTBells identified: %d, SMRTBell subseqs: %d", m_TotProcessed, m_TotAccepted, m_TotRejected, m_TotUnderLen, m_TotSMRTBells, m_TotSMRTBellSubSeqs);
else
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Processed: %d - accepted: %d, rejected: %d, under length: %d, vector trimmed: %d", m_TotProcessed, m_TotAccepted, m_TotRejected, m_TotUnderLen, m_TotContamTrimmed);

pThreadPar = pThreadPutOvlps;
for(ThreadIdx = 0; ThreadIdx < NumOvlpThreads; ThreadIdx++,pThreadPar++)
	{
	if(pThreadPar->pSW != NULL)
		delete pThreadPar->pSW;
	}

delete pThreadPutOvlps;


if(m_hOutFile != -1)
	{
	if(m_OutBuffIdx > 0)
		CUtility::SafeWrite(m_hOutFile,m_pszOutBuff,m_OutBuffIdx);

#ifdef _WIN32
	_commit(m_hOutFile);
#else
	fsync(m_hOutFile);
#endif
	close(m_hOutFile);
	m_hOutFile = -1;
	delete m_pszOutBuff;
	m_pszOutBuff = NULL;
	}

return(0);
}


int
CPBFilter::PBFilterReads(tsThreadPBFilter *pThreadPar)
{
int Rslt;
int SeqID;
char szQuerySeqIdent[100];
UINT8 *pQuerySeq;
int QuerySeqLen;

int TrimStartOfs;
int TrimLen;
int Trim5AtOfs;
int Trim3AtOfs;
int Trim5AtContam;
int Trim3AtContam;

tsSMRTBellHit BellHits[cMaxSMRTBellHits+1];
tsSMRTBellHit *pBellHit;
int BellIdx;
int PrevBellTargOfs;
int TrimSeqLen;

// iterating over all sequences
while((pQuerySeq = m_PacBioUtility.DequeueQuerySeq(60,sizeof(szQuerySeqIdent),&SeqID,szQuerySeqIdent,&QuerySeqLen))!=NULL)				
	{
	// remove near homopolymer runs
	QuerySeqLen = m_PacBioUtility.ReduceHomopolymers(QuerySeqLen,pQuerySeq,16);

	AcquireLock(true);
	m_TotProcessed += 1;
	if(QuerySeqLen < m_MinReadLen + m_Trim5 + m_Trim3) // simply slough if would be under length after any end trimming
		{
		m_TotUnderLen += 1;
		m_TotRejected += 1;
		ReleaseLock(true);
		delete pQuerySeq;
		continue;
		}
	ReleaseLock(true);
	Trim5AtOfs = m_Trim5;
	Trim3AtOfs = QuerySeqLen - (m_Trim3 + 1);
	TrimLen = 1 + Trim3AtOfs - Trim5AtOfs;

	if(m_PMode == ePBPMFilter)
		{
		int NumBellHits;
		// Note: more chance of finding 5' and 3' end SMRTBells if no end trimming ... End trim later ...
		NumBellHits = m_PacBioUtility.IdentifySMRTBells(pThreadPar->SMRTBellSensitivity,cMaxSMRTBellHits,QuerySeqLen,pQuerySeq,BellHits);
		if(NumBellHits > 0)
			{
			AcquireLock(true);
			m_TotSMRTBells += NumBellHits;
			if(NumBellHits > 4)		// if more than 4 retained SMRTbells then don't attempt to trim
				{
				m_TotRejected += 1;
				ReleaseLock(true);
				delete pQuerySeq;
				continue;
				}
			ReleaseLock(true);
			}

		int NumSubSeqs = 0;
		PrevBellTargOfs = Trim5AtOfs;
		BellIdx = 0;
		if(NumBellHits > 0)
			{
			pBellHit = BellHits;
			for(BellIdx = 0; BellIdx < NumBellHits; BellIdx++,pBellHit+=1)
				{
				TrimSeqLen = pBellHit->TargOfs - PrevBellTargOfs - 25; // allowing 25bp float on SMRTBell loci
				if(TrimSeqLen < m_MinReadLen)
					{
					PrevBellTargOfs += 25;
					continue;
					}
				NumSubSeqs += 1;
				if((Rslt=WriteFastaFile(SeqID,NumSubSeqs,TrimSeqLen,&pQuerySeq[PrevBellTargOfs]))!=eBSFSuccess)
					{
					delete pQuerySeq;
					return(Rslt);
					}
				
				PrevBellTargOfs = pBellHit->TargOfs + 25;			// allowing 25bp float on SMRTBell loci
				}
			}
		TrimSeqLen = 1 + Trim3AtOfs - PrevBellTargOfs;
		if(TrimSeqLen < m_MinReadLen)
			{
			AcquireLock(true);
			if(NumSubSeqs > 0)
				{
				m_TotAccepted += 1;
				m_TotSMRTBellSubSeqs += NumSubSeqs;
				}
			else
				m_TotRejected += 1;
			ReleaseLock(true);
			delete pQuerySeq;
			continue;
			}
		
		if(NumBellHits > 0)
			NumSubSeqs += 1;

		if((Rslt=WriteFastaFile(SeqID,NumSubSeqs,TrimSeqLen,&pQuerySeq[PrevBellTargOfs]))!=eBSFSuccess)
			{
			delete pQuerySeq;
			return(Rslt);
			}

		AcquireLock(true);
		m_TotAccepted += 1;	
		m_TotSMRTBellSubSeqs += NumSubSeqs;
		ReleaseLock(true);
		delete pQuerySeq;
		continue;
		}

	Trim5AtContam = 0;
	Trim3AtContam = QuerySeqLen - 1;

	if(pThreadPar->SWAlignInstance == 0)
		pThreadPar->SWAlignInstance = m_pSWAlign->InitInstance();
	tsSSWCell Matched;
	if(m_pSWAlign->AlignProbeSeq(pThreadPar->SWAlignInstance,QuerySeqLen,(etSeqBase *)pQuerySeq,1,0,false,&Matched)!=eBSFSuccess)
		{
		int Trim5ContamLen;
		int Trim3ContamLen;
		bool bTrimContam = false;
		if(Matched.StartPOfs > 0 && Matched.EndPOfs > 0)
			{
			Matched.StartPOfs -= 1;
			Matched.EndPOfs -= 1;
			if((int)Matched.StartPOfs > Trim5AtOfs)
				Trim5ContamLen = (int)Matched.StartPOfs - Trim5AtOfs;
			else
				Trim5ContamLen = 0;
			if((int)Matched.EndPOfs < Trim3AtOfs)
				Trim3ContamLen = Trim3AtOfs - (int)Matched.EndPOfs;
			else
				Trim3ContamLen = 0;

			if(Trim5ContamLen >= m_MinReadLen || Trim3ContamLen >= m_MinReadLen)
				{
				if(Trim5ContamLen >= Trim3ContamLen)
					{
					Trim3AtContam = Matched.StartPOfs;
					Trim5AtContam = 0;
					}
				else
					{
					Trim3AtContam = QuerySeqLen - 1;
					Trim5AtContam = Matched.EndPOfs;
					}
				AcquireLock(true);
				m_TotContamTrimmed += 1;
				ReleaseLock(true);
				bTrimContam = true;
				}
			}
		
		if(!bTrimContam)
			{
			AcquireLock(true);
			m_TotContamRejected += 1;
			ReleaseLock(true);
			delete pQuerySeq;
			continue;
			}
		}

	TrimStartOfs = max(Trim5AtOfs, Trim5AtContam);
	TrimLen = 1 + min(Trim3AtOfs, Trim3AtContam) - TrimStartOfs;
	if(TrimLen < m_MinReadLen)
		{
		AcquireLock(true);
		m_TotUnderLen += 1;
		ReleaseLock(true);
		delete pQuerySeq;
		continue;
		}	

	if((Rslt=WriteFastaFile(SeqID,0,TrimLen,&pQuerySeq[TrimStartOfs]))!=eBSFSuccess)
		{
		delete pQuerySeq;
		return(Rslt);
		}

	AcquireLock(true);
	m_TotAccepted += 1;
	ReleaseLock(true);
	delete pQuerySeq;
	}
return(eBSFSuccess);
}


int										// eBSFSuccess if sequence written to file otherwise the error code
CPBFilter::WriteFastaFile(int SeqID,	// primary sequence identifier
			int SecSeqID,				// secondary sequence identifier
			int SeqLen,					// sequence to write is this length
			etSeqBase *pSeq)			// sequence
{
int SeqIdx;
int LineLen;
char Base;
char *pBuff;
etSeqBase *pBase;

if(SeqLen < 1)
	return(eBSFSuccess);

AcquireSerialise();

if(m_hOutFile == -1)
	{
	ReleaseSerialise();
	return(eBSFerrClosed);
	}

if(m_OutBuffIdx > (m_AllocOutBuffSize - (SeqLen * 2)))
	{
	if(!CUtility::SafeWrite(m_hOutFile,m_pszOutBuff,m_OutBuffIdx))
		{
		ReleaseSerialise();
		return(eBSFerrFileAccess);
		}
	m_OutBuffIdx = 0;
	}
m_OutBuffIdx += sprintf(&m_pszOutBuff[m_OutBuffIdx],">FiltSeq%d_%d %d\n",SeqID,SecSeqID,SeqLen);
LineLen = 0;
pBuff = &m_pszOutBuff[m_OutBuffIdx]; 
pBase = pSeq;
for(SeqIdx = 0; SeqIdx < SeqLen; SeqIdx++, pBase+=1)
	{
	switch(*pBase) {
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
			break;
			}
	*pBuff++ = Base;
	m_OutBuffIdx += 1;
	LineLen += 1;
	if(LineLen > 80)
		{
		*pBuff++ = '\n';
		m_OutBuffIdx += 1;
		LineLen = 0;
		}
	}	
*pBuff++ = '\n';
m_OutBuffIdx += 1;
ReleaseSerialise();
return(eBSFSuccess);
}



// sort putative SMRTBell adapters by StartTOfs ascending with NumExacts descending
int
CPBFilter::SortSSWCells(const void *arg1, const void *arg2)
{
tsSSWCell *pEl1 = (tsSSWCell *)arg1;
tsSSWCell *pEl2 = (tsSSWCell *)arg2;

if(pEl1->StartTOfs < pEl2->StartTOfs)	
	return(-1);
if(pEl1->StartTOfs > pEl2->StartTOfs)
	return(1);
if(pEl1->NumExacts > pEl2->NumExacts)	
	return(-1);
if(pEl1->NumExacts < pEl2->NumExacts)
	return(1);

return(0);
}


