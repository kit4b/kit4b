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

#include "./ngskit4b.h"
#include "CallHaplotypes.h"


int CallHaplotypes(eModeCSH PMode,	// processing mode 0: call observed haplotype counts, mode 1: GBS processing
			char *pszTrackName,		// track name
			char *pszTrackDescr,	// track descriptor
			uint32_t FndrTrim5,			// trim this many aligned PBAs from 5' end of founder aligned segments - reduces false alleles due to sequencing errors
			uint32_t FndrTrim3,			// trim this many aligned PBAs from 3' end of founder aligned segments - reduces false alleles due to sequencing errors
			uint32_t ProgTrim5,			// trim this many aligned PBAs from 5' end of progeny aligned segments - reduces false alleles due to sequencing errors
			uint32_t ProgTrim3,			// trim this many aligned PBAs from 3' end of progeny aligned segments - reduces false alleles due to sequencing errors
			char *pszMaskBPAFile,	// optional masking input BPA file, only process BPAs which are intersect of these BPAs and progeny plus founder BPAs
			int NumFounderInputFiles,	// number of input founder file specs
			char *pszFounderInputFiles[],	// names of input founder PBA files (wildcards allowed)
			int NumProgenyInputFiles,	// number of input progeny file specs
			char* pszProgenyInputFiles[],		// names of input progeny PBA files (wildcards allowed)
			char* pszOutFile,		// Windowed haplotype calls output file (CSV format)
			int NumThreads);		// number of worker threads to use

#ifdef _WIN32
int callhaplotypes(int argc, char *argv[])
{
	// determine my process name
	_splitpath (argv[0], NULL, NULL, gszProcName, NULL);
#else
int
callhaplotypes(int argc, char **argv)
{
	// determine my process name
	CUtility::splitpath ((char *)argv[0], NULL, gszProcName);
#endif
	int iFileLogLevel;			// level of file diagnostics
	int iScreenLogLevel;		// level of file diagnostics
	char szLogFile[_MAX_PATH];	// write diagnostics to this file
	int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure
	int NumberOfProcessors;		// number of installed CPUs
	int NumThreads;				// number of threads (0 defaults to number of CPUs or a maximum of cMaxPBAWorkerThreads)

	int Idx;

	eModeCSH PMode;				// processing mode
	 uint32_t FndrTrim5;		// trim this many aligned PBAs from 5' end of founder aligned segments - reduces false alleles due to sequencing errors
	 uint32_t FndrTrim3;		// trim this many aligned PBAs from 3' end of founder aligned segments - reduces false alleles due to sequencing errors
	 uint32_t ProgTrim5;		// trim this many aligned PBAs from 5' end of progeny aligned segments - reduces false alleles due to sequencing errors
	 uint32_t ProgTrim3;		// trim this many aligned PBAs from 3' end of progeny aligned segments - reduces false alleles due to sequencing errors
	 char szTrackName[_MAX_PATH];	// track name
	 char szTrackDescr[_MAX_PATH];	// track description
	 int NumFounderInputFiles;		// number of input founder BPA files
	 char *pszFounderInputFiles[cMaxFounderFileSpecs];		// names of input founder BPA files (wildcards allowed)
	 int NumProgenyInputFiles;		// number of input progeny BPA files
	 char *pszProgenyInputFiles[cMaxProgenyFileSpecs];		// names of input progeny BPA files (wildcards allowed)

	 char szOutFile[_MAX_PATH];		// Windowed haplotype calls output file base name, progeny readset identifier is appended to this base name
	 char szMaskBPAFile[_MAX_PATH];	// optional masking input BPA file, only process BPAs which are intersect of these BPAs and progeny plus founder BPAs

	struct arg_lit *help = arg_lit0 ("h", "help", "print this help and exit");
	struct arg_lit *version = arg_lit0 ("v", "version,ver", "print version information and exit");
	struct arg_int *FileLogLevel = arg_int0 ("f", "FileLogLevel", "<int>", "Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
	struct arg_file *LogFile = arg_file0 ("F", "log", "<file>", "diagnostics log file");

	struct arg_int *pmode = arg_int0 ("m", "mode", "<int>", "processing mode: 0 default");

	struct arg_int* fndrtrim5 = arg_int0("y", "fndrtrim5", "<int>", "trim this many aligned PBAs from 5' end of founder aligned segments (default 10)");
	struct arg_int* fndrtrim3 = arg_int0("Y", "fndrtrim3", "<int>", "trim this many aligned PBAs from 3' end of founder aligned segments (default 10)");
	struct arg_int* progtrim5 = arg_int0("w", "progtrim5", "<int>", "trim this many aligned PBAs from 5' end of progeny aligned segments (default 5)");
	struct arg_int* progtrim3 = arg_int0("W", "progtrim3", "<int>", "trim this many aligned PBAs from 3' end of progeny aligned segments (default 5)");

	struct arg_str *trackname = arg_str0("t","trackname","<str>","BED Track name");
	struct arg_str *trackdescr = arg_str0("D","trackdescr","<str>","BED Track description");
	struct arg_file *founderfiles = arg_filen("I", "founderfiles", "<file>", 0,cMaxFounderFileSpecs,"founder input BPA file(s), wildcards allowed, limit of 100 founder filespecs supported");
	struct arg_file *progenyfiles = arg_filen("i", "inprogenyfile", "<file>",0, cMaxProgenyFileSpecs, "progeny input BPA file(s), wildcards allowed, limit of 500 progeny filespecs supported");
	struct arg_file *maskbpafile = arg_file0("c", "inmaskbpa", "<file>", "optional masking input BPA file, only process BPAs which are intersect of these BPAs and progeny plus founder BPAs");
	struct arg_file *outfile = arg_file1("o", "out", "<file>", "Windowed haplotype calls output prefix (outputs CSV, BED and WIG format)");
	struct arg_int *threads = arg_int0("T","threads","<int>","number of processing threads 0..64 (defaults to 0 which limits threads to maximum of 64 CPU cores)");
	struct arg_end *end = arg_end (200);

	void *argtable[] = { help,version,FileLogLevel,LogFile,
						pmode,fndrtrim5,fndrtrim3,progtrim5,progtrim3,
						trackname,trackdescr,maskbpafile,progenyfiles,founderfiles, outfile,threads,end };

	char **pAllArgs;
	int argerrors;
	argerrors = CUtility::arg_parsefromfile (argc, (char **)argv, &pAllArgs);
	if (argerrors >= 0)
		argerrors = arg_parse (argerrors, pAllArgs, argtable);

	/* special case: '--help' takes precedence over error reporting */
	if (help->count > 0)
	{
		printf ("\n%s %s %s, Version %s\nOptions ---\n", gszProcName, gpszSubProcess->pszName, gpszSubProcess->pszFullDescr, kit4bversion);
		arg_print_syntax (stdout, argtable, "\n");
		arg_print_glossary (stdout, argtable, "  %-25s %s\n");
		printf ("\nNote: Parameters can be entered into a parameter file, one parameter per line.");
		printf ("\n      To invoke this parameter file then precede its name with '@'");
		printf ("\n      e.g. %s %s @myparams.txt\n", gszProcName, gpszSubProcess->pszName);
		printf ("\nPlease report any issues regarding usage of %s at https://github.com/kit4b/issues\n\n", gszProcName);
		return(1);
	}

	/* special case: '--version' takes precedence error reporting */
	if (version->count > 0)
		{
		printf ("\n%s %s Version %s\n", gszProcName, gpszSubProcess->pszName, kit4bversion);
		return(1);
		}

	if (!argerrors)
	{
		if (FileLogLevel->count && !LogFile->count)
		{
			printf ("\nError: FileLogLevel '-f%d' specified but no logfile '-F<logfile>\n'", FileLogLevel->ival[0]);
			exit (1);
		}

		iScreenLogLevel = iFileLogLevel = FileLogLevel->count ? FileLogLevel->ival[0] : eDLInfo;
		if (iFileLogLevel < eDLNone || iFileLogLevel > eDLDebug)
		{
			printf ("\nError: FileLogLevel '-l%d' specified outside of range %d..%d\n", iFileLogLevel, eDLNone, eDLDebug);
			exit (1);
		}

		if (LogFile->count)
		{
			strncpy (szLogFile, LogFile->filename[0], _MAX_PATH);
			szLogFile[_MAX_PATH - 1] = '\0';
		}
		else
		{
		iFileLogLevel = eDLNone;
		szLogFile[0] = '\0';
		}

		// now that log parameters have been parsed then initialise diagnostics log system
		if(!gDiagnostics.Open(szLogFile, (etDiagLevel)iScreenLogLevel, (etDiagLevel)iFileLogLevel, true))
		{
			printf("\nError: Unable to start diagnostics subsystem\n");
			if(szLogFile[0] != '\0')
				printf(" Most likely cause is that logfile '%s' can't be opened/created\n", szLogFile);
			exit(1);
		}

		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Subprocess %s Version %s starting", gpszSubProcess->pszName, kit4bversion);
		gExperimentID = 0;
		gProcessID = 0;
		gProcessingID = 0;

		PMode = pmode->count ? (eModeCSH)pmode->ival[0] : eMCSHDefault;
		if(PMode < eMCSHDefault || PMode >= eMCSHPlaceHolder)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Unsupported processing mode '-m%d'\n",PMode);
			exit(1);
			}

		FndrTrim5 = fndrtrim5->count ? fndrtrim5->ival[0] : cDfltFndrTrim5;
		if(FndrTrim5 < 0)
			FndrTrim5 = 0;
		else
			if(FndrTrim5 > 1000)
				FndrTrim5 = 1000;

		FndrTrim3 = fndrtrim3->count ? fndrtrim3->ival[0] : cDfltFndrTrim3;
		if(FndrTrim3 < 0)
			FndrTrim3 = 0;
		else
			if(FndrTrim3 > 1000)
				FndrTrim3 = 1000;

		ProgTrim5 = progtrim5->count ? progtrim5->ival[0] : cDfltProgTrim5;
		if(ProgTrim5 < 0)
			ProgTrim5 = 0;
		else
			if(ProgTrim5 > 1000)
				ProgTrim5 = 1000;

		ProgTrim3 = progtrim3->count ? progtrim3->ival[0] : cDfltProgTrim3;
		if(ProgTrim3 < 0)
			ProgTrim3 = 0;
		else
			if(ProgTrim3 > 1000)
				ProgTrim3 = 1000;

		if(trackname->count)
		{
			strncpy(szTrackName, trackname->sval[0], 80);
			szTrackName[80] = '\0';
			CUtility::TrimQuotedWhitespcExtd(szTrackName);
			CUtility::CleanText(szTrackName);
		}
		else
			szTrackName[0] = '\0';
		if(szTrackName[0] == '\0')
			strcpy(szTrackName, "SH");


		if(trackdescr->count)
		{
			strncpy(szTrackDescr, trackdescr->sval[0], 80);
			szTrackDescr[80] = '\0';
			CUtility::TrimQuotedWhitespcExtd(szTrackDescr);
			CUtility::CleanText(szTrackDescr);
		}
		else
			szTrackDescr[0] = '\0';
		if(szTrackDescr[0] == '\0')
			strcpy(szTrackDescr, szTrackName);
		szTrackName[40] = '\0';

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

	int MaxAllowedThreads = min(cMaxPBAWorkerThreads,NumberOfProcessors);	// limit to be at most cMaxWorkerThreads
	if((NumThreads = threads->count ? threads->ival[0] : MaxAllowedThreads)==0)
		NumThreads = MaxAllowedThreads;
	if(NumThreads < 0 || NumThreads > MaxAllowedThreads)
		{
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"Warning: Number of threads '-T%d' specified was outside of range %d..%d",NumThreads,1,MaxAllowedThreads);
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"Warning: Defaulting number of threads to %d",MaxAllowedThreads);
		NumThreads = MaxAllowedThreads;
		}

	NumProgenyInputFiles = 0;

	if(progenyfiles->count)
		{
		for(Idx = 0; NumProgenyInputFiles < cMaxProgenyReadsets && Idx < progenyfiles->count; Idx++)
			{
			pszProgenyInputFiles[Idx] = NULL;
			if(pszProgenyInputFiles[NumProgenyInputFiles] == NULL)
				pszProgenyInputFiles[NumProgenyInputFiles] = new char[_MAX_PATH];
			strncpy(pszProgenyInputFiles[NumProgenyInputFiles], progenyfiles->filename[Idx], _MAX_PATH);
			pszProgenyInputFiles[NumProgenyInputFiles][_MAX_PATH - 1] = '\0';
			CUtility::TrimQuotedWhitespcExtd(pszProgenyInputFiles[NumProgenyInputFiles]);
			if(pszProgenyInputFiles[NumProgenyInputFiles][0] != '\0')
				NumProgenyInputFiles++;
			}
		}


	if(!NumProgenyInputFiles)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: After removal of whitespace, there are no input progeny file(s) specified with '-i<filespec>' option)\n");
		exit(1);
		}

	NumFounderInputFiles = 0;
	if(founderfiles->count)
		{
		for(NumFounderInputFiles = Idx = 0; NumFounderInputFiles < cMaxFounderFileSpecs && Idx < founderfiles->count; Idx++)
			{
			pszFounderInputFiles[Idx] = NULL;
			if(pszFounderInputFiles[NumFounderInputFiles] == NULL)
				pszFounderInputFiles[NumFounderInputFiles] = new char[_MAX_PATH];
			strncpy(pszFounderInputFiles[NumFounderInputFiles], founderfiles->filename[Idx], _MAX_PATH);
			pszFounderInputFiles[NumFounderInputFiles][_MAX_PATH - 1] = '\0';
			CUtility::TrimQuotedWhitespcExtd(pszFounderInputFiles[NumFounderInputFiles]);
			if(pszFounderInputFiles[NumFounderInputFiles][0] != '\0')
				NumFounderInputFiles++;
			}
		}

	if(!NumFounderInputFiles)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: After removal of whitespace, there are no input founder file(s) specified with '-I<filespec>' option)\n");
		exit(1);
		}

	if(PMode == eMCSHDefault && maskbpafile->count)
		{
		strcpy (szMaskBPAFile, maskbpafile->filename[0]);
		CUtility::TrimQuotedWhitespcExtd (szMaskBPAFile);
		if(szMaskBPAFile[0] == '\0')
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: After removal of whitespace, no input masking BPA file specified with '-c<filespec>' option)\n");
			exit(1);
			}
		}
	else
		szMaskBPAFile[0] = '\0';

	strcpy (szOutFile, outfile->filename[0]);
	CUtility::TrimQuotedWhitespcExtd (szOutFile);
	if (szOutFile[0] == '\0')
		{
		gDiagnostics.DiagOut (eDLFatal, gszProcName, "No output file specified");
		exit (1);
		}

	gDiagnostics.DiagOut (eDLInfo, gszProcName, "Processing parameters:");
	const char *pszDescr;
	switch (PMode) {
		case eMCSHDefault:
			pszDescr = "Generate founder unique allele stacks and process progeny alleles against founder allele stacks";
			break;

		}

	gDiagnostics.DiagOutMsgOnly (eDLInfo, "Calling haplotypes : '%s'", pszDescr);
	gDiagnostics.DiagOutMsgOnly(eDLInfo, "Trim this many aligned PBAs from 5' end of founder aligned segments: %d", FndrTrim5);
	gDiagnostics.DiagOutMsgOnly(eDLInfo, "Trim this many aligned PBAs from 3' end of founder aligned segments: %d", FndrTrim3);
	gDiagnostics.DiagOutMsgOnly(eDLInfo, "Trim this many aligned PBAs from 5' end of progeny aligned segments: %d", ProgTrim5);
	gDiagnostics.DiagOutMsgOnly(eDLInfo, "Trim this many aligned PBAs from 3' end of progeny aligned segments: %d", ProgTrim3);

	if(szMaskBPAFile[0] != '\0')
		gDiagnostics.DiagOutMsgOnly (eDLInfo, "Masking BPAs loaded from : '%s'", szMaskBPAFile);
		

	for(Idx = 0; Idx < NumFounderInputFiles; Idx++)
		gDiagnostics.DiagOutMsgOnly (eDLInfo, "Founder file : '%s'", pszFounderInputFiles[Idx]);

	for(Idx = 0; Idx < NumProgenyInputFiles; Idx++)
		gDiagnostics.DiagOutMsgOnly (eDLInfo, "Progeny file : '%s'", pszProgenyInputFiles[Idx]);
	gDiagnostics.DiagOutMsgOnly (eDLInfo, "Output file : '%s'", szOutFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"number of threads : %d",NumThreads);

#ifdef _WIN32
	SetPriorityClass (GetCurrentProcess (), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start ();
	Rslt = 0;
	Rslt = CallHaplotypes(PMode,			// processing mode
					szTrackName,			// track name
					szTrackDescr,			// track descriptor
					FndrTrim5,			// trim this many aligned PBAs from 5' end of founder aligned segments - reduces false alleles due to sequencing errors
					FndrTrim3,			// trim this many aligned PBAs from 3' end of founder aligned segments - reduces false alleles due to sequencing errors
					ProgTrim5,			// trim this many aligned PBAs from 5' end of progeny aligned segments - reduces false alleles due to sequencing errors
					ProgTrim3,			// trim this many aligned PBAs from 3' end of progeny aligned segments - reduces false alleles due to sequencing errors
					szMaskBPAFile,			// optional input masking BPA file, only process BPAs which are intersect of these BPAs and progeny plus founder BPAs
					NumFounderInputFiles,	// number of input founder file specs
					pszFounderInputFiles,	// names of input founder PBA files (wildcards allowed)
					NumProgenyInputFiles,		// number of input progeny file specs
					pszProgenyInputFiles,		// names of input progeny PBA files (wildcards allowed)
					szOutFile,				// output to this file
					NumThreads);			// number of worker threads to use
	Rslt = Rslt >= 0 ? 0 : 1;
	gStopWatch.Stop ();

	gDiagnostics.DiagOut (eDLInfo, gszProcName, "Exit code: %d Total processing time: %s", Rslt, gStopWatch.Read ());
	exit (Rslt);
	}
else
	{
	printf ("\n%s %s %s, Version %s\n", gszProcName, gpszSubProcess->pszName, gpszSubProcess->pszFullDescr, kit4bversion);
	arg_print_errors (stdout, end, gszProcName);
	arg_print_syntax (stdout, argtable, "\nUse '-h' to view option and parameter usage\n");
	exit (1);
	}
}

int CallHaplotypes(eModeCSH PMode,	// processing mode 0: call observed haplotype counts, mode 1: GBS processing
				   char* pszTrackName,		// track name
				   char* pszTrackDescr,	// track descriptor
				   uint32_t FndrTrim5,			// trim this many aligned PBAs from 5' end of founder aligned segments - reduces false alleles due to sequencing errors
				   uint32_t FndrTrim3,			// trim this many aligned PBAs from 3' end of founder aligned segments - reduces false alleles due to sequencing errors
				   uint32_t ProgTrim5,			// trim this many aligned PBAs from 5' end of progeny aligned segments - reduces false alleles due to sequencing errors
				   uint32_t ProgTrim3,			// trim this many aligned PBAs from 3' end of progeny aligned segments - reduces false alleles due to sequencing errors
				   char* pszMaskBPAFile,	// optional masking input BPA file, only process BPAs which are intersect of these BPAs and progeny plus founder BPAs
				   int NumFounderInputFiles,	// number of input founder file specs
				   char* pszFounderInputFiles[],	// names of input founder PBA files (wildcards allowed)
				   int NumProgenyInputFiles,	// number of input progeny file specs
				   char* pszProgenyInputFiles[],		// names of input progeny PBA files (wildcards allowed)
				   char* pszOutFile,		// Windowed haplotype calls output file (CSV format)
				   int NumThreads)		// number of worker threads to use
{
int Rslt;
CCallHaplotypes* pCallHaplotypes;

if((pCallHaplotypes = new CCallHaplotypes) == NULL)
	{
	gDiagnostics.DiagOut (eDLFatal, gszProcName, "Unable to instantiate instance of CCallHaplotypes");
	return(eBSFerrObj);
	}
Rslt = pCallHaplotypes->Process(PMode,pszTrackName,pszTrackDescr, FndrTrim5, FndrTrim3, ProgTrim5, ProgTrim3, 
				pszMaskBPAFile,NumFounderInputFiles,
				pszFounderInputFiles,NumProgenyInputFiles,pszProgenyInputFiles,pszOutFile,NumThreads);
delete pCallHaplotypes;
return(Rslt);
}


CCallHaplotypes::CCallHaplotypes()
{
m_pChromMetadata = NULL;
m_pProgenyFndrAligns = NULL;
m_pWorkQueueEls = NULL;
m_pInBuffer = NULL;
m_pszOutBuffer = NULL;
m_pAlleleStacks = NULL;
m_hOutFile = -1;
m_hInFile = -1;
m_bMutexesCreated = false;
Reset();
}

CCallHaplotypes::~CCallHaplotypes()
{
if(m_hInFile != -1)
	close(m_hInFile);
if(m_hOutFile != -1)
	close(m_hOutFile);
if(m_pWorkQueueEls != NULL)
	delete []m_pWorkQueueEls;
if(m_pszOutBuffer != NULL)
	delete []m_pszOutBuffer;
if(m_pInBuffer != NULL)
	delete []m_pInBuffer;

if(m_pChromMetadata != NULL)
	{
	tsChromMetadata *pChromMetadata = m_pChromMetadata;
	for(uint32_t Idx = 0; Idx < m_UsedNumChromMetadata; Idx++, pChromMetadata++)
		{
		if(pChromMetadata->pPBAs != NULL)
#ifdef _WIN32
			free(pChromMetadata->pPBAs);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
		if(pChromMetadata->pPBAs != MAP_FAILED)
			munmap(pChromMetadata->pPBAs, pChromMetadata->ChromLen);
#endif
		}

#ifdef _WIN32
	free(m_pChromMetadata);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pChromMetadata != MAP_FAILED)
		munmap(m_pChromMetadata, m_AllocdChromMetadataMem);
#endif
	}

if(m_pProgenyFndrAligns != NULL)
	{
#ifdef _WIN32
	free(m_pProgenyFndrAligns);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pProgenyFndrAligns != MAP_FAILED)
		munmap(m_pProgenyFndrAligns, m_AllocdProgenyFndrAlignsMem);
#endif
	}

if(m_pAlleleStacks != NULL)
{
#ifdef _WIN32
	free(m_pAlleleStacks);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pAlleleStacks != MAP_FAILED)
		munmap(m_pAlleleStacks, m_AllocdAlleleStacksMem);
#endif
	}
DeleteMutexes();
}

void
CCallHaplotypes::Reset(void)	// resets class instance state back to that immediately following instantiation
{
if(m_hInFile != -1)
	{
	close(m_hInFile);
	m_hInFile = -1;
	}
if(m_hOutFile != -1)
	{
	close(m_hOutFile);
	m_hOutFile = -1;
	}

if(m_pWorkQueueEls != NULL)
	{
	delete []m_pWorkQueueEls;
	m_pWorkQueueEls = NULL;
	}
m_AllocWorkQueueEls = 0;
m_TotWorkQueueEls = 0;
m_NumQueueElsProcessed = 0;

if(m_pInBuffer != NULL)
	{
	delete []m_pInBuffer;
	m_pInBuffer = NULL;
	}

if(m_pszOutBuffer != NULL)
	{
	delete []m_pszOutBuffer;
	m_pszOutBuffer = NULL;
	}

if(m_pChromMetadata != NULL)
	{
	tsChromMetadata *pChromMetadata = m_pChromMetadata;
	for(uint32_t Idx = 0; Idx < m_UsedNumChromMetadata; Idx++, pChromMetadata++)
		{
		if(pChromMetadata->pPBAs != NULL)
#ifdef _WIN32
			free(pChromMetadata->pPBAs);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
		if(pChromMetadata->pPBAs != MAP_FAILED)
			munmap(pChromMetadata->pPBAs, pChromMetadata->ChromLen);
#endif
		}

#ifdef _WIN32
	free(m_pChromMetadata);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pChromMetadata != MAP_FAILED)
		munmap(m_pChromMetadata, m_AllocdChromMetadataMem);
#endif
	m_pChromMetadata = NULL;
	}
m_UsedNumChromMetadata = 0;
m_AllocdChromMetadata = 0;
m_AllocdChromMetadataMem = 0;

if(m_pAlleleStacks != NULL)
	{
#ifdef _WIN32
	free(m_pAlleleStacks);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pAlleleStacks != MAP_FAILED)
		munmap(m_pAlleleStacks, m_AllocdAlleleStacksMem);
#endif
	m_pAlleleStacks = NULL;
	}
m_AllocdAlleleStacksMem = 0;
m_UsedAlleleStacks = 0;
m_AllocdAlleleStacks = 0;


DeleteMutexes();
if(m_pProgenyFndrAligns != NULL)
	{
#ifdef _WIN32
	free(m_pProgenyFndrAligns);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pProgenyFndrAligns != MAP_FAILED)
		munmap(m_pProgenyFndrAligns, m_AllocdProgenyFndrAlignsMem);
#endif
	m_pProgenyFndrAligns = NULL;
	}
m_UsedProgenyFndrAligns = 0;
m_AllocdProgenyFndrAligns = 0;
m_AllocdProgenyFndrAlignsMem = 0;

m_NumFounders = 0;
m_NumProgenies = 0;

m_FndrTrim5 = cDfltFndrTrim5;
m_FndrTrim3 = cDfltFndrTrim3;
m_ProgTrim5 = cDfltProgTrim5;
m_ProgTrim3 = cDfltProgTrim3;

m_InNumBuffered = 0;
m_AllocInBuff = 0;

m_OutBuffIdx = 0;	
m_AllocOutBuff = 0;

m_LAReadsetNameID = 0;
m_NumReadsetNames = 0;
m_NxtszReadsetIdx = 0;
m_szReadsetNames[0] = '\0';

m_LAChromNameID = 0;
m_NumChromNames = 0;
m_NxtszChromIdx = 0;
m_szChromNames[0] = '\0';

m_NumFounders = 0;
m_FndrsProcMap = (uint64_t)-1;

if(m_bMutexesCreated)
	DeleteMutexes();
m_bMutexesCreated = false;
}


int
CCallHaplotypes::CreateMutexes(void)
{
if(m_bMutexesCreated)
	return(eBSFSuccess);

#ifdef _WIN32
if((m_hSerialiseAccess = CreateMutex(NULL,false,NULL))==NULL)
	{
#else
if(pthread_mutex_init (&m_hSerialiseAccess,NULL)!=0)
	{
#endif
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to create mutex");
	return(eBSFerrInternal);
	}
m_bMutexesCreated = true;
return(eBSFSuccess);
}

void
CCallHaplotypes::DeleteMutexes(void)
{
if(!m_bMutexesCreated)
	return;
#ifdef _WIN32
CloseHandle(m_hSerialiseAccess);
#else
pthread_mutex_destroy(&m_hSerialiseAccess);
#endif
m_bMutexesCreated = false;
}

void
CCallHaplotypes::AcquireSerialise(void)
{
uint32_t WaitRslt;
#ifdef _WIN32
WaitRslt = WaitForSingleObject(m_hSerialiseAccess,INFINITE);
if(WaitRslt !=  WAIT_OBJECT_0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal:WaitForSingleObject() returned error %u",WaitRslt);
	Reset();
	exit(1);
	}
#else
pthread_mutex_lock(&m_hSerialiseAccess);
#endif
}

void
CCallHaplotypes::ReleaseSerialise(void)
{
#ifdef _WIN32
ReleaseMutex(m_hSerialiseAccess);
#else
pthread_mutex_unlock(&m_hSerialiseAccess);
#endif
}


int
CCallHaplotypes::Process(eModeCSH PMode,	// processing mode 0: call observed haplotype counts, mode 1: GBS processing
				char* pszTrackName,		// track name
				char* pszTrackDescr,	// track descriptor
				uint32_t FndrTrim5,			// trim this many aligned PBAs from 5' end of founder aligned segments - reduces false alleles due to sequencing errors
				uint32_t FndrTrim3,			// trim this many aligned PBAs from 3' end of founder aligned segments - reduces false alleles due to sequencing errors
				uint32_t ProgTrim5,			// trim this many aligned PBAs from 5' end of progeny aligned segments - reduces false alleles due to sequencing errors
				uint32_t ProgTrim3,			// trim this many aligned PBAs from 3' end of progeny aligned segments - reduces false alleles due to sequencing errors
				char* pszMaskBPAFile,	// optional masking input BPA file, only process BPAs which are intersect of these BPAs and progeny plus founder BPAs
				int NumFounderInputFiles,	// number of input founder file specs
				char* pszFounderInputFiles[],	// names of input founder PBA files (wildcards allowed)
				int NumProgenyInputFiles,	// number of input progeny file specs
				char* pszProgenyInputFiles[],		// names of input progeny PBA files (wildcards allowed)
				char* pszOutFile,		// Windowed haplotype calls output file (CSV format)
				int NumThreads)		// number of worker threads to use
{
int Rslt;
int NumFiles;
int TotNumFiles;
size_t memreq;

Reset();

CreateMutexes();
m_FndrTrim5 = FndrTrim5;
m_FndrTrim3 = FndrTrim3;
m_ProgTrim5 = ProgTrim5;
m_ProgTrim3 = ProgTrim3;

m_NumThreads = NumThreads;

if(NumFounderInputFiles > cMaxFounderFileSpecs)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Can accept at most %d founder wildcarded file specs for processing, %d requested", cMaxFounderFileSpecs, NumFounderInputFiles);
	Reset();
	return(eBSFerrMem);
	}

if(NumProgenyInputFiles > cMaxProgenyFileSpecs)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Can accept at most %d progeny wildcarded file specs for processing, %d requested", cMaxFounderFileSpecs, NumProgenyInputFiles);
	Reset();
	return(eBSFerrMem);
	}

if(pszMaskBPAFile == NULL || pszMaskBPAFile[0] == '\0')
	m_pszMaskBPAFile = NULL;
else
	m_pszMaskBPAFile = pszMaskBPAFile;


if((m_pInBuffer = new uint8_t[cInBuffSize]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Memory allocation of %u bytes for input buffering failed - %s", cInBuffSize, strerror(errno));
	Reset();
	return(eBSFerrMem);
	}
m_AllocInBuff = cInBuffSize;
m_InNumBuffered = 0;

// initial allocation, will be realloc'd as required if more memory required
memreq = (size_t)cAllocChromMetadata * sizeof(tsChromMetadata);
#ifdef _WIN32
m_pChromMetadata = (tsChromMetadata *)malloc(memreq);	// initial and perhaps the only allocation
if(m_pChromMetadata == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Initial memory allocation of %lld bytes for chromosome metadata failed - %s", (int64_t)memreq, strerror(errno));
	Reset();
	return(eBSFerrMem);
	}
#else
m_pChromMetadata = (tsChromMetadata *)mmap(NULL, memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
if(m_pChromMetadata == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Memory allocation of %lld bytes through mmap() for chromosome metadata failed - %s", (int64_t)memreq, strerror(errno));
	m_pChromMetadata = NULL;
	Reset();
	return(eBSFerrMem);
	}
#endif
m_AllocdChromMetadataMem = memreq;
m_UsedNumChromMetadata = 0;
m_AllocdChromMetadata = cAllocChromMetadata;


memreq = (size_t)cAllocProgenyFndrAligns * sizeof(tsProgenyFndrAligns);
#ifdef _WIN32
m_pProgenyFndrAligns = (tsProgenyFndrAligns*)malloc(memreq);	// initial and perhaps the only allocation
if(m_pProgenyFndrAligns == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Initial memory allocation of %lld bytes for chromosome windowed founder counts failed - %s", (int64_t)memreq, strerror(errno));
	Reset();
	return(eBSFerrMem);
	}
#else
m_pProgenyFndrAligns = (tsProgenyFndrAligns*)mmap(NULL, memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
if(m_pProgenyFndrAligns == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Memory allocation of %lld bytes through mmap() for chromosome windowed founder counts failed - %s", (int64_t)memreq, strerror(errno));
	m_pProgenyFndrAligns = NULL;
	Reset();
	return(eBSFerrMem);
	}
#endif
m_AllocdProgenyFndrAlignsMem = memreq;
m_UsedProgenyFndrAligns = 0;
m_AllocdProgenyFndrAligns = cAllocProgenyFndrAligns;


memreq = (size_t)cAllocAlleleStacks * sizeof(tsAlleleStack);
#ifdef _WIN32
m_pAlleleStacks = (tsAlleleStack*)malloc(memreq);	// initial and perhaps the only allocation
if(m_pAlleleStacks == NULL)
{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Initial memory allocation of %lld bytes for allele stacks failed - %s", (int64_t)memreq, strerror(errno));
	Reset();
	return(eBSFerrMem);
}
#else
m_pAlleleStacks = (tsAlleleStack*)mmap(NULL, memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
if(m_pAlleleStacks == MAP_FAILED)
{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Memory allocation of %lld bytes through mmap() for allele stacks failed - %s", (int64_t)memreq, strerror(errno));
	m_pAlleleStacks = NULL;
	Reset();
	return(eBSFerrMem);
}
#endif
m_AllocdAlleleStacksMem = memreq;
m_UsedAlleleStacks = 0;
m_AllocdAlleleStacks = cAllocAlleleStacks;

gDiagnostics.DiagOut(eDLInfo, gszProcName, "Process: Starting to load pool of founders");
Rslt = eBSFSuccess;		// assume success!
CSimpleGlob glob(SG_GLOB_FULLSORT);
int Idx;
char* pszInFile;
int32_t ReadsetID = 0;
TotNumFiles = 0;
for(Idx = 0; Idx < NumFounderInputFiles; Idx++)	
	{
	glob.Init();
	if(glob.Add(pszFounderInputFiles[Idx]) < SG_SUCCESS)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to glob '%s' founder pool files", pszFounderInputFiles[Idx]);
		Reset();
		return(eBSFerrOpnFile);	// treat as though unable to open file
		}
	if((NumFiles = glob.FileCount()) <= 0)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to locate any founder pool file matching '%s", pszFounderInputFiles[Idx]);
		Reset();
		return(eBSFerrFileName);
		}
	TotNumFiles += NumFiles;
	if(TotNumFiles > cMaxFounderReadsets)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Can accept at most %d founder pool readsets for processing, after wildcard file name expansions there are %d requested", cMaxFounderReadsets, TotNumFiles);
		Reset();
		return(eBSFerrMem);
		}

	Rslt = eBSFSuccess;
	for(int FileID = 0; Rslt >= eBSFSuccess && FileID < NumFiles; ++FileID)
		{
		pszInFile = glob.File(FileID);
		if((ReadsetID = LoadPBAFile(pszInFile,0)) <= 0)
			{
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "Process: Errors loading pool file '%s'",pszInFile);
			Reset();
			return(ReadsetID);
			}
		m_Fndrs2Proc[ReadsetID-1] = 0x01;	// if loaded then assumption is that this founder will be processed
		}
	}
m_NumFounders = (uint32_t)ReadsetID;
gDiagnostics.DiagOut(eDLInfo, gszProcName, "Process: Completed loading founders (%d) pool", m_NumFounders);

// load the control PBA file if it has been specified
m_MaskReadsetID = 0;
if(m_pszMaskBPAFile != NULL)
	{
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Process: Loading control PBA file '%s'",m_pszMaskBPAFile);
	if((ReadsetID = LoadPBAFile(m_pszMaskBPAFile,2)) <= 0)
		{
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Process: Errors loading control PBA file");
		Reset();
		return(ReadsetID);
		}
	m_MaskReadsetID = (uint32_t)ReadsetID;
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Completed loading control PBA file");
	}

if((Rslt = GenAlleleStacks(m_NumFounders)) < 1)
	{
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Process: No allele stacks generated");
	Reset();
	return(Rslt);
	}

// next is to individually process each progeny PBAs against the founder panel allele stacks
TotNumFiles = 0;
m_CurProgenyReadsetID = 0;
for(Idx = 0; Idx < NumProgenyInputFiles; Idx++)	
	{
	glob.Init();
	if(glob.Add(pszProgenyInputFiles[Idx]) < SG_SUCCESS)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to glob progeny '%s", pszProgenyInputFiles[Idx]);
		Reset();
		return(eBSFerrOpnFile);	// treat as though unable to open file
		}
	if((NumFiles = glob.FileCount()) <= 0)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to locate any progeny PBA file matching '%s", pszProgenyInputFiles[Idx]);
		Reset();
		return(eBSFerrOpnFile);	// treat as though unable to open file
		}

	TotNumFiles += NumFiles;
	if(TotNumFiles > cMaxProgenyReadsets)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Can accept at most %d progeny readsets for processing, after wildcard file name expansions there are %d requested", cMaxProgenyReadsets, TotNumFiles);
		Reset();
		return(eBSFerrMem);
		}

	Rslt = eBSFSuccess;
	for(int FileID = 0; Rslt >= eBSFSuccess && FileID < NumFiles; ++FileID)
		{
		pszInFile = glob.File(FileID);
		if((Rslt = ProcessProgenyPBAFile(pszInFile,pszOutFile)) < eBSFSuccess)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Failed processing progeny PBA file '%s", pszInFile);
			Reset();
			return(Rslt);
			}
		}
	}
if(m_UsedProgenyFndrAligns == 0)
	{
	gDiagnostics.DiagOut(eDLWarn, gszProcName, "Process: No progeny to  founder allele stack alignments");
	Reset();
	return(eBSFSuccess);
	}

// sort progeny allele stack overlaps by readset.chrom.loci ascending
m_mtqsort.SetMaxThreads(m_NumThreads);
m_mtqsort.qsort(m_pProgenyFndrAligns, (int64_t)m_UsedProgenyFndrAligns, sizeof(tsProgenyFndrAligns), SortProgenyFndrAligns);

// report haplotypes present for individual progenies
for(uint32_t ProgIdx = 0; ProgIdx < m_NumProgenies; ProgIdx++)
	ReportHaplotypesByProgeny(pszOutFile, m_ProgenyIDs[ProgIdx]);

// generate a matrix with chrom.loci (Y) by progeny (X)
m_mtqsort.qsort(m_pProgenyFndrAligns, (int64_t)m_UsedProgenyFndrAligns, sizeof(tsProgenyFndrAligns), SortProgenyChromLociReadset);

ReportMatrix(pszOutFile);

Reset();
return(Rslt);
}


int 
CCallHaplotypes::ReportHaplotypesByProgeny(char* pszRsltsFileBaseName,		// haplotype results are written to this file base name with '.haplotypes.csv' appended
											uint32_t ReadsetID)				// report on this progeny readset only, or if 0 then report on all progeny readsets
{
char szOutFile[_MAX_PATH];
tsProgenyFndrAligns *pCurPFA;
uint32_t PFAIdx;
uint32_t FndrIdx;

int32_t PrevReadsetID;
int32_t PrevChromID;

char *pszProgenyReadset;
char *pszChrom;

if(m_pszOutBuffer == NULL)
	{
	if((m_pszOutBuffer = new uint8_t[cOutBuffSize]) == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "ReportHaplotypesByProgeny: Memory allocation of %u bytes for output buffering failed - %s", cOutBuffSize, strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
	m_AllocOutBuff = cOutBuffSize;
	}
m_OutBuffIdx = 0;
if(ReadsetID == 0)
	sprintf(szOutFile, "%s.progeny.all.csv", pszRsltsFileBaseName);
else
	sprintf(szOutFile, "%s.progeny.%s.csv", pszRsltsFileBaseName, LocateReadset(ReadsetID));
#ifdef _WIN32
m_hOutFile = open(szOutFile, (O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC), (_S_IREAD | _S_IWRITE));
#else
if((m_hOutFile = open64(szOutFile, O_WRONLY | O_CREAT, S_IREAD | S_IWRITE)) != -1)
if(ftruncate(m_hOutFile, 0) != 0)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to create/truncate %s - %s", szOutFile, strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}
#endif
if(m_hOutFile < 0)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "ReportHaplotypesByProgeny: Unable to create/truncate %s - %s", szOutFile, strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}

// header line
m_OutBuffIdx = sprintf((char*)m_pszOutBuffer, "\"Progeny\",\"Chrom\",\"Loci\"");
for(FndrIdx = 0; FndrIdx < m_NumFounders; FndrIdx++)
	{
	if(m_Fndrs2Proc[FndrIdx] & 0x01)
		m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], ",\"Fndr:%s\"", LocateReadset(FndrIdx+1));
	}
m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], "\n");

PrevReadsetID = 0;
PrevChromID = 0;
pCurPFA = &m_pProgenyFndrAligns[0];
for(PFAIdx = 0; PFAIdx < m_UsedProgenyFndrAligns; PFAIdx++, pCurPFA++)
	{
	if(ReadsetID != 0 && pCurPFA->ReadsetID != ReadsetID)
		continue;
	if(pCurPFA->ReadsetID != PrevReadsetID || pCurPFA->ChromID != PrevChromID)
		{
		pszProgenyReadset = LocateReadset(pCurPFA->ReadsetID);
		pszChrom = LocateChrom(pCurPFA->ChromID);
		PrevReadsetID = pCurPFA->ReadsetID;
		PrevChromID = pCurPFA->ChromID;
		}

	m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx],"\"%s\",\"%s\",%u", pszProgenyReadset, pszChrom, pCurPFA->Loci);
	for(FndrIdx=0; FndrIdx < m_NumFounders; FndrIdx++)
		{
		if(m_Fndrs2Proc[FndrIdx] & 0x01)
			m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], ",%d", Bits256Test(FndrIdx, pCurPFA->ProgenyFounders) ? 1 : 0);
		}
	m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], "\n");
	if((m_OutBuffIdx + 1000) > m_AllocOutBuff)
		{
		if(!CUtility::RetryWrites(m_hOutFile, m_pszOutBuffer, m_OutBuffIdx))
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "ReportMatrix: Fatal error in RetryWrites()");
			Reset();
			return(eBSFerrFileAccess);
			}
		m_OutBuffIdx = 0;
		}
	}
if(m_OutBuffIdx && m_hOutFile != -1)
	{
	if(!CUtility::RetryWrites(m_hOutFile, m_pszOutBuffer, m_OutBuffIdx))
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "ReportMatrix: Fatal error in RetryWrites()");
		Reset();
		return(eBSFerrFileAccess);
		}
	m_OutBuffIdx = 0;
	}

if(m_hOutFile != -1)
	{
// commit output file
#ifdef _WIN32
	_commit(m_hOutFile);
#else
	fsync(m_hOutFile);
#endif
	close(m_hOutFile);
	m_hOutFile = -1;
	}

return(eBSFSuccess);
}

int
CCallHaplotypes::ReportMatrix(char* pszRsltsFileBaseName)		// matrix results are written to this file base name with 'matrix.csv' appended
{
char szOutFile[_MAX_PATH];
tsProgenyFndrAligns* pCurPFA;
uint32_t PFAIdx;
uint32_t FndrIdx;
uint32_t ProgIdx;
uint8_t ProgenyHaplotypes[cMaxProgenyReadsets];
uint8_t ProgHaplotype;
int32_t CurChromID;
uint32_t CurLoci;
bool bNewLoci;

char* pszChrom;

if(m_pszOutBuffer == NULL)
{
	if((m_pszOutBuffer = new uint8_t[cOutBuffSize]) == NULL)
	{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "ReportMatrix: Memory allocation of %u bytes for output buffering failed - %s", cOutBuffSize, strerror(errno));
		Reset();
		return(eBSFerrMem);
	}
	m_AllocOutBuff = cOutBuffSize;
}
m_OutBuffIdx = 0;
sprintf(szOutFile, "%s.matrix.csv", pszRsltsFileBaseName);
#ifdef _WIN32
	m_hOutFile = open(szOutFile, (O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC), (_S_IREAD | _S_IWRITE));
#else
if((m_hOutFile = open64(szOutFile, O_WRONLY | O_CREAT, S_IREAD | S_IWRITE)) != -1)
	if(ftruncate(m_hOutFile, 0) != 0)
	{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to create/truncate %s - %s", szOutFile, strerror(errno));
		Reset();
		return(eBSFerrCreateFile);
	}
#endif
if(m_hOutFile < 0)
{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "ReportMatrix: Unable to create/truncate %s - %s", szOutFile, strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
}

// header line
m_OutBuffIdx = sprintf((char*)m_pszOutBuffer, "\"Chrom\",\"Loci\"");

for(ProgIdx = 0; ProgIdx < m_NumProgenies; ProgIdx++)
	m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], ",\"Progeny:%s\"", LocateReadset(m_ProgenyIDs[ProgIdx]));

bNewLoci = false;
CurLoci = 0;
CurChromID = 0;
pCurPFA = &m_pProgenyFndrAligns[0];
for(PFAIdx = 0; PFAIdx < m_UsedProgenyFndrAligns; PFAIdx++, pCurPFA++)
	{
	if(pCurPFA->ChromID != CurChromID || pCurPFA->Loci != CurLoci)
		{
		if(bNewLoci)
			{
			m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], "\n\"%s\",%u", pszChrom, CurLoci);
			for(ProgIdx = 0; ProgIdx < m_NumProgenies; ProgIdx++)
				m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], ",%d", ProgenyHaplotypes[ProgIdx]);
			if((m_OutBuffIdx + 10000) > m_AllocOutBuff)
			{
				if(!CUtility::RetryWrites(m_hOutFile, m_pszOutBuffer, m_OutBuffIdx))
				{
					gDiagnostics.DiagOut(eDLFatal, gszProcName, "ReportMatrix: Fatal error in RetryWrites()");
					Reset();
					return(eBSFerrFileAccess);
				}
				m_OutBuffIdx = 0;
			}
			bNewLoci = false;
			}
		CurLoci = pCurPFA->Loci;
		CurChromID = pCurPFA->ChromID;
		pszChrom = LocateChrom(CurChromID);
		memset(ProgenyHaplotypes,0, m_NumProgenies);
		}
	bNewLoci = true;
	// locate corresponding offset into ProgenyHaplotypes[] for the progeny readset
	for(ProgIdx = 0; ProgIdx < m_NumProgenies; ProgIdx++)
		{
		if(pCurPFA->ReadsetID == m_ProgenyIDs[ProgIdx])
			break;
		}
	if(ProgIdx >= m_NumProgenies)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error locating progeny identifier");
		Reset();
		return(eBSFerrInternal);
		}

	ProgHaplotype = 0;
	for(FndrIdx = m_NumFounders; FndrIdx >= 1; FndrIdx--)
		{
		ProgHaplotype <<= 1;
		if(m_Fndrs2Proc[FndrIdx-1] & 0x01)
			ProgHaplotype |= Bits256Test(FndrIdx-1, pCurPFA->ProgenyFounders) ? 1 : 0;
		}
	ProgenyHaplotypes[ProgIdx] = ProgHaplotype;
	}
if(bNewLoci)
	{
	m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], "\n\"%s\",%u", pszChrom, CurLoci);
	for(ProgIdx = 0; ProgIdx < m_NumProgenies; ProgIdx++)
		m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], ",%d", ProgenyHaplotypes[ProgIdx]);
	}
if(m_OutBuffIdx && m_hOutFile != -1)
	{
	if(!CUtility::RetryWrites(m_hOutFile, m_pszOutBuffer, m_OutBuffIdx))
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "ReportMatrix: Fatal error in RetryWrites()");
		Reset();
		return(eBSFerrFileAccess);
		}
	m_OutBuffIdx = 0;
	}

if(m_hOutFile != -1)
	{
// commit output file
#ifdef _WIN32
		_commit(m_hOutFile);
#else
		fsync(m_hOutFile);
#endif
		close(m_hOutFile);
		m_hOutFile = -1;
	}

return(eBSFSuccess);
}


int
CCallHaplotypes::ProcessProgenyPBAFile(char* pszProgenyPBAFile,	// load, process and call haplotypes for this progeny PBA file against previously loaded panel founder PBAs
						char *pszRsltsFileBaseName)			// results are written to this file base name with progeny readset identifier and type appended
{
m_CurProgenyReadsetID = 0;
m_LAChromNameID = 0;
m_LAReadsetNameID = 0;
if((m_CurProgenyReadsetID = LoadPBAFile(pszProgenyPBAFile,1)) <= 0)
	{
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessProgenyPBAFile: Errors loading progeny PBA file '%s'",pszProgenyPBAFile);
	Reset();
	return(-1);
	}
gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessProgenyPBAFile: Completed loading progeny PBA file");
m_ProgenyIDs[m_NumProgenies++] = m_CurProgenyReadsetID;
if(m_pszOutBuffer == NULL)
	{
	if((m_pszOutBuffer = new uint8_t[cOutBuffSize]) == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessProgenyPBAFile: Memory allocation of %u bytes for output buffering failed - %s", cOutBuffSize, strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
	m_AllocOutBuff = cOutBuffSize;
	}
m_OutBuffIdx = 0;

char *pszProgenyReadset = LocateReadset(m_CurProgenyReadsetID);

tsReadsetMetadata *pProgeny = &m_Readsets[m_CurProgenyReadsetID -1];
tsProgenyFndrAligns ProgenyFndrAligns;
uint32_t NumOverlapping;
uint8_t ProgenyAlleles;
uint8_t Allele;
bool bAlleleAccept;
uint32_t AlleleStackIdx;
tsAlleleStack *pCurAlleleStack;
uint32_t AlleleIdx;
uint8_t AlleleMsk;

uint32_t NumBi = 0;
uint32_t NumMono = 0;
uint32_t NumFa = 0;
uint32_t NumFb = 0;
NumOverlapping = 0;
pCurAlleleStack = m_pAlleleStacks;
for(AlleleStackIdx = 0; AlleleStackIdx < m_UsedAlleleStacks; AlleleStackIdx++, pCurAlleleStack++)
	{
	ProgenyAlleles = LocateReadsetChromLociAlleles(m_CurProgenyReadsetID, pCurAlleleStack->ChromID,pCurAlleleStack->Loci);
	if(ProgenyAlleles)
		{
		memset(&ProgenyFndrAligns, 0, sizeof(tsProgenyFndrAligns));
		bAlleleAccept = false;
		AlleleMsk = 0x03;
		for(AlleleIdx = 0; AlleleIdx < 4; AlleleIdx++, AlleleMsk <<= 2)
			{
			Allele = (ProgenyAlleles & AlleleMsk) >> (AlleleIdx * 2);
			if(Allele < 1)		// treating non-significant alleles as noise and discarding these
				continue;
			if(pCurAlleleStack->NumAlleleFndrs[AlleleIdx]==0) // if no founders having this allele then treat as introgression and treat as if no progeny allele
				{
				bAlleleAccept = false;
				break;
				}
			ProgenyFndrAligns.NumProgenyFounders = Bits256Combine(ProgenyFndrAligns.ProgenyFounders, pCurAlleleStack->Alleles[AlleleIdx]);
			bAlleleAccept = true;
			}
		if(!bAlleleAccept)
			continue;

		if(ProgenyFndrAligns.NumProgenyFounders == 2)
			NumBi++;
		else
			NumMono++;

		if(Bits256Test(0,ProgenyFndrAligns.ProgenyFounders))
			NumFa++;
		if(Bits256Test(1,ProgenyFndrAligns.ProgenyFounders))
			NumFb++;

		ProgenyFndrAligns.Alleles = ProgenyAlleles;
		ProgenyFndrAligns.ChromID = pCurAlleleStack->ChromID;
		ProgenyFndrAligns.Loci = pCurAlleleStack->Loci;
		ProgenyFndrAligns.ReadsetID = m_CurProgenyReadsetID;
		AddProgenyFndrAligns(&ProgenyFndrAligns);
		NumOverlapping++;
		}
	}

gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessProgenyPBAFile: Located %d progeny loci overlapping %d founder AlleleStacks loci", NumOverlapping, m_UsedAlleleStacks);
return(eBSFSuccess);
}


uint32_t				// returns number of unprocessed bytes in buffer
CCallHaplotypes::FillInBuffer(uint32_t MinRequired) // try and fill input buffer with at least this many bytes reading from currently opened input file handle
{
int NumRead;
uint32_t UnProcessed;
UnProcessed = m_InNumBuffered - m_InNumProcessed;

// if already filled to MinRequired then no further action required
if(MinRequired <= UnProcessed)
	return(UnProcessed);

// copy down any bytes which are yet to be processed
if(UnProcessed)	{
	memmove(m_pInBuffer, &m_pInBuffer[m_InNumProcessed], UnProcessed);
	m_InNumBuffered = UnProcessed;
	m_InNumProcessed = 0;
}

// attempt to fill input buffer to it's allocation maximum
do{
	if((NumRead = read(m_hInFile, &m_pInBuffer[m_InNumProcessed], m_AllocInBuff - m_InNumBuffered)) < 0)		
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error %s attempting to read from file", strerror(errno));
		return(0);
		}
	m_InNumBuffered += (uint32_t)NumRead;
}
while(NumRead > 0 && m_InNumBuffered < MinRequired);
return(m_InNumBuffered);
}

int32_t					// returned readset identifier (1..n) or < 0 if errors
CCallHaplotypes::LoadPBAFile(char* pszFile,	// load chromosome metadata and PBA data from this file
							uint8_t ReadsetType)	// 0: founder, 1: progeny, 2: control
{
int Rslt;
int Version;
char szExperimentID[100];
char szRefAssemblyID[100];
char szReadsetID[100];
uint32_t PrevChromMetadataIdx;

uint32_t ChromID;
uint32_t ReadsetID;
tsReadsetMetadata *pReadsetMetadata;
tsChromMetadata *pChromMetadata;
tsChromMetadata *pPrevChromMetadata;
uint8_t* pBuff;
int scanlen;
int NumTags;

gDiagnostics.DiagOut(eDLInfo, gszProcName, "Loading '%s' file",pszFile);
#ifdef _WIN32
m_hInFile = open(pszFile, O_READSEQ);		// file access is normally sequential..
#else
m_hInFile = open64(pszFile, O_READSEQ);		// file access is normally sequential..
#endif
if(m_hInFile == -1)							// check if file open succeeded
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to open input file '%s' : %s", pszFile, strerror(errno));
	return(eBSFerrOpnFile);
	}

// attempt to maximally fill input buffer
m_InNumBuffered = 0;
m_InNumProcessed = 0;
if(FillInBuffer(m_AllocInBuff) == 0 || m_InNumBuffered < 9)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to read at least a partial header from input file '%s'", pszFile);
	return(eBSFerrOpnFile);
	}

// check file type is PbA
if(strncmp((char *)m_pInBuffer,"Type:PbA\n",9))
	{
	gDiagnostics.DiagOut (eDLFatal, gszProcName, "Input file '%s' exists, unable to parse file type header tag, is it a packed base allele file",pszFile);
	return(eBSFerrOpnFile);
	}
m_InNumProcessed = 9;

if(m_InNumBuffered < 500)
	{
	gDiagnostics.DiagOut (eDLFatal, gszProcName, "Input file '%s' exists, unable to read complete header, is it a packed base allele file",pszFile);
	return(eBSFerrOpnFile);
	}

// parse out tagnames and associated values
pBuff = &m_pInBuffer[m_InNumProcessed];
NumTags=sscanf((char *)pBuff,"Version:%d\nExperimentID:%99[^\n]\nReferenceID:%99[^\n]\nReadsetID:%99[^\n]%n",&Version,szExperimentID,szRefAssemblyID,szReadsetID,&scanlen);
if(NumTags != 4 || Version != 1 || szExperimentID[0] == '\0' || szRefAssemblyID[0] == '\0' || szReadsetID[0] == '\0')
	{
	gDiagnostics.DiagOut (eDLFatal, gszProcName, "Input file '%s' exists as a packed base allele file but inconsistencies in header tag values",pszFile);
	return(eBSFerrOpnFile);
	}
gDiagnostics.DiagOut(eDLInfo, gszProcName, "Metadata ExperimentID:%s, ReferenceID:%s, ReadsetID:%s",szExperimentID,szRefAssemblyID,szReadsetID);

m_InNumProcessed += scanlen + 1;		// header tags were '\0' terminated 
if((ReadsetID = AddReadset(szReadsetID, ReadsetType))==0)
	{
	gDiagnostics.DiagOut (eDLFatal, gszProcName, "Input file '%s' duplicates the ReadsetID '%s' of a previously loaded readset",pszFile,szReadsetID);
	return(eBSFerrOpnFile);
	}

pReadsetMetadata = &m_Readsets[ReadsetID-1];
memset(pReadsetMetadata,0,sizeof(*pReadsetMetadata));
pReadsetMetadata->ReadsetType = ReadsetType;
pReadsetMetadata->NumChroms = 0;
strcpy(pReadsetMetadata->szExperimentID,szExperimentID);
strcpy(pReadsetMetadata->szRefAssemblyID,szRefAssemblyID);
pReadsetMetadata->ReadsetID = ReadsetID;
pReadsetMetadata->StartChromID = 0;
pReadsetMetadata->StartChromMetadataIdx = 0;

int ChromNameLen;
char *pszChromName;
uint32_t ChromLen;

// iterate over all chromosomes

PrevChromMetadataIdx = 0;
while((m_InNumProcessed + 110) <= m_InNumBuffered)
	{
	pBuff = &m_pInBuffer[m_InNumProcessed++];
	ChromNameLen = (int)*pBuff++;
	pszChromName = (char *)pBuff;
	pBuff += 1+ChromNameLen;
	ChromLen = *(uint32_t *)pBuff;
	pBuff+=4;
	m_InNumProcessed += ChromNameLen + 1 + sizeof(uint32_t);
	if((m_InNumBuffered - m_InNumProcessed) < ChromLen && FillInBuffer(ChromLen) == 0)
		break;

	ChromID = AddChrom(pszChromName);

	if((Rslt = AllocChromMetadata()) != eBSFSuccess)
		return(Rslt);

	pChromMetadata = &m_pChromMetadata[m_UsedNumChromMetadata++];
	if(PrevChromMetadataIdx != 0)
		{
		pPrevChromMetadata = &m_pChromMetadata[PrevChromMetadataIdx-1];
		pPrevChromMetadata->NxtChromMetadataIdx = m_UsedNumChromMetadata;
		}
	else
		{
		pReadsetMetadata->StartChromID = ChromID;
		pReadsetMetadata->StartChromMetadataIdx = m_UsedNumChromMetadata;
		}
	pReadsetMetadata->NumChroms++;
	PrevChromMetadataIdx = m_UsedNumChromMetadata;
	pChromMetadata->ChromID = ChromID;
	pChromMetadata->ChromMetadataIdx = m_UsedNumChromMetadata;
	pChromMetadata->NxtChromMetadataIdx = 0;
	pChromMetadata->ChromLen = ChromLen;
	pChromMetadata->ReadsetID = ReadsetID;
	pChromMetadata->pPBAs = AllocPBAs(ChromLen);
	memcpy(pChromMetadata->pPBAs, pBuff, ChromLen);
	if(ReadsetType == 0)
		TrimPBAs(m_FndrTrim5, m_FndrTrim3, ChromLen,pChromMetadata->pPBAs);
	else // treating controls as if progeny when trimming
		TrimPBAs(m_ProgTrim5, m_ProgTrim3, ChromLen, pChromMetadata->pPBAs);
	m_InNumProcessed += ChromLen;
	if((m_InNumBuffered - m_InNumProcessed) < 110 && FillInBuffer(ChromLen) == 0)
		break;
	}

close(m_hInFile);
m_hInFile = -1;
return((int32_t)ReadsetID);
}

// trim 5' and 3' aligned segments within the PBAs attempting to reduce sequencing error induced false alleles
int			// returns number of non-trimmed loci in the pPBAs
CCallHaplotypes::TrimPBAs(uint32_t Trim5,	// trim 5' this many aligned PBA bases from each aligned segment
		 uint32_t Trim3,	// trim 3' this many aligned PBA bases from each aligned segment
		 uint32_t PBALen,	// pPBAs contains this many packed base alleles
		uint8_t* pPBAs)		// packed base alleles to be processed
{
uint32_t NonTrimmed;
uint32_t ToTrim;
uint32_t Loci;
uint8_t* pPBA;
if(pPBAs == NULL || PBALen == 0)
	return(0);
if(Trim5 == 0 && Trim3 == 0)
	return(PBALen);
NonTrimmed = PBALen;

// first the 5' trim
if(Trim5)
	{
	ToTrim = Trim5;
	pPBA = pPBAs;
	for(Loci=0;Loci < PBALen;Loci++, pPBA++)
		{
		if(*pPBA == 0)
			{
			ToTrim = Trim5;
			continue;
			}
		if(ToTrim)
			{
			ToTrim-=1;
			NonTrimmed-=1;
			*pPBA = 0;
			}
		}
	}
// then the 3' trim
if(Trim3)
	{
	ToTrim = Trim3;
	pPBA = &pPBAs[PBALen-1];
	for(Loci = 0; Loci < PBALen; Loci++, pPBA--)
		{
		if(*pPBA == 0)
			{
			ToTrim = Trim3;
			continue;
			}
		if(ToTrim)
			{
			NonTrimmed-=1;
			ToTrim -= 1;
			*pPBA = 0;
			}
		}
	}
return(NonTrimmed);
}


int
CCallHaplotypes::AllocChromMetadata(void)
{
uint32_t ToAllocdChromMetadata;
tsChromMetadata *pChromMetadata;
size_t memreq;
if (m_pChromMetadata == NULL)					// may be NULL first time in
	{
	memreq = cAllocChromMetadata * sizeof(tsChromMetadata);
#ifdef _WIN32
	m_pChromMetadata = (tsChromMetadata *)malloc((size_t)memreq);
	if (m_pChromMetadata == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AllocChromMetadata: Memory allocation of %lld bytes failed", (int64_t)memreq);
		return(eBSFerrMem);
		}
#else
	m_pChromMetadata = (tsChromMetadata *)mmap(NULL, (size_t)memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	if (m_pChromMetadata == MAP_FAILED)
		{
		m_pChromMetadata = NULL;
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AllocChromMetadata: Memory allocation of %lld bytes through mmap()  failed", (int64_t)memreq, strerror(errno));
		return(eBSFerrMem);
		}
#endif
	m_AllocdChromMetadataMem = memreq;
	m_AllocdChromMetadata = cAllocChromMetadata;
	m_UsedNumChromMetadata = 0;
	}
else
		// needing to allocate more memory?
	if (m_UsedNumChromMetadata == m_AllocdChromMetadata)
		{
		ToAllocdChromMetadata = m_UsedNumChromMetadata + cAllocChromMetadata;
		size_t memreq = ToAllocdChromMetadata * sizeof(tsChromMetadata);
#ifdef _WIN32
		pChromMetadata = (tsChromMetadata*)realloc(m_pChromMetadata, memreq);
		if (pChromMetadata == NULL)
			{
#else
			pChromMetadata = (tsChromMetadata*)mremap(m_pChromMetadata, m_AllocdChromMetadataMem, memreq, MREMAP_MAYMOVE);
			if (pChromMetadata == MAP_FAILED)
			{
#endif
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "AllocChromMetadata: Memory reallocation to %lld bytes failed - %s", memreq, strerror(errno));
			return(eBSFerrMem);
			}
		m_pChromMetadata = pChromMetadata;
		m_AllocdChromMetadataMem = memreq;
		m_AllocdChromMetadata=ToAllocdChromMetadata;
		}
return(eBSFSuccess);
}

uint8_t *
CCallHaplotypes::AllocPBAs(uint32_t ChromLen)	// allocate memory to hold at least this many packed base alleles
{
uint8_t *pPBAs;
size_t memreq;
memreq = (size_t)ChromLen;	// no safety margin!
#ifdef _WIN32
pPBAs = (uint8_t*)malloc((size_t)memreq);
if (pPBAs == NULL)
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "AllocPBAs Memory allocation of %lld bytes failed", (int64_t)memreq);
#else
pPBAs = (uint8_t*)mmap(NULL, (size_t)memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
if (pPBAs == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "AllocPBAs: Memory allocation of %lld bytes through mmap()  failed", (int64_t)memreq, strerror(errno));
	pPBAs = NULL;
	}
#endif
return(pPBAs);
}


uint8_t *								// returned pointer to start of PBA
CCallHaplotypes::LocatePBAfor(int32_t ReadSetID,		// readset identifier 
			 int32_t ChromID)			// chrom identifier
{
tsReadsetMetadata *pReadsetMetadata;
tsChromMetadata *pChromMetadata;
uint32_t CurChromMetadataIdx;

if(ReadSetID > m_NumReadsetNames)
	return(NULL);
pReadsetMetadata = &m_Readsets[ReadSetID-1];
CurChromMetadataIdx = pReadsetMetadata->StartChromMetadataIdx;
for(uint32_t ChromIdx = 0; ChromIdx < pReadsetMetadata->NumChroms && CurChromMetadataIdx != 0; ChromIdx++)
	{
	pChromMetadata = &m_pChromMetadata[CurChromMetadataIdx - 1];
	if(pChromMetadata->ChromID == ChromID)
		return(pChromMetadata->pPBAs);
	CurChromMetadataIdx = pChromMetadata->NxtChromMetadataIdx;
	}
return(NULL);
}



uint32_t									// returned index+1 into  m_pProgenyFndrAligns[] to allocated and initialised ProgenyFndrAligns, 0 if errors
CCallHaplotypes::AddProgenyFndrAligns(tsProgenyFndrAligns *pInitProgenyFndrAligns)	// allocated tsProgenyFndrAligns to be initialised with a copy of pInitProgenyFndrAligns
{
uint32_t ToAllocdProgenyFndrAligns;
tsProgenyFndrAligns *pProgenyFndrAligns;
size_t memreq;
AcquireSerialise();
if (m_pProgenyFndrAligns == NULL)					// may be NULL first time in
	{
	memreq = cAllocProgenyFndrAligns * sizeof(tsProgenyFndrAligns);
#ifdef _WIN32
	m_pProgenyFndrAligns = (tsProgenyFndrAligns*)malloc((size_t)memreq);
	if (m_pProgenyFndrAligns == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddWinBinCnts: Memory allocation of %lld bytes failed", (int64_t)memreq);
		return(0);
		}
#else
	m_pProgenyFndrAligns = (tsProgenyFndrAligns*)mmap(NULL, (size_t)memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	if (m_pProgenyFndrAligns == MAP_FAILED)
		{
		m_pProgenyFndrAligns = NULL;
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddWinBinCnts: Memory allocation of %lld bytes through mmap()  failed", (int64_t)memreq, strerror(errno));
		return(0);
		}
#endif
	m_AllocdProgenyFndrAlignsMem = memreq;
	m_AllocdProgenyFndrAligns = cAllocProgenyFndrAligns;
	m_UsedProgenyFndrAligns = 0;
	}
else
		// needing to allocate more memory?
	if ((m_UsedProgenyFndrAligns) >= m_AllocdProgenyFndrAligns)
		{
		ToAllocdProgenyFndrAligns = m_UsedProgenyFndrAligns + cAllocProgenyFndrAligns;
		size_t memreq = ToAllocdProgenyFndrAligns * sizeof(tsProgenyFndrAligns);
#ifdef _WIN32
		pProgenyFndrAligns = (tsProgenyFndrAligns*)realloc(m_pProgenyFndrAligns, memreq);
		if (pProgenyFndrAligns == NULL)
			{
#else
		pProgenyFndrAligns = (tsProgenyFndrAligns*)mremap(m_pProgenyFndrAligns, m_AllocdProgenyFndrAlignsMem, memreq, MREMAP_MAYMOVE);
		if (pProgenyFndrAligns == MAP_FAILED)
			{
#endif
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddProgenyFndrAligns: Memory reallocation to %lld bytes failed - %s", memreq, strerror(errno));
			return(0);
			}
		m_pProgenyFndrAligns = pProgenyFndrAligns;
		m_AllocdProgenyFndrAlignsMem = memreq;
		m_AllocdProgenyFndrAligns =ToAllocdProgenyFndrAligns;
		}
pProgenyFndrAligns = &m_pProgenyFndrAligns[m_UsedProgenyFndrAligns++];
if(pInitProgenyFndrAligns != NULL)
	*pProgenyFndrAligns = *pInitProgenyFndrAligns;
else
	memset(pProgenyFndrAligns,0,sizeof(tsProgenyFndrAligns));
ReleaseSerialise();
return(m_UsedProgenyFndrAligns);
}


uint32_t									// returned index+1 into m_pAlleleStacks[] to allocated and initialised allele stack, 0 if errors
CCallHaplotypes::AddAlleleStack(tsAlleleStack* pInitAlleleStack)	// allocated tsAlleleStack to be initialised with a copy of pInitAlleleStack
{
uint32_t ToAllocdAlleleStacks;
tsAlleleStack* pAlleleStack;
size_t memreq;
AcquireSerialise();
if(m_pAlleleStacks == NULL)					// may be NULL first time in
	{
	memreq = (size_t)cAllocAlleleStacks * sizeof(tsAlleleStack);
#ifdef _WIN32
	m_pAlleleStacks = (tsAlleleStack*)malloc((size_t)memreq);
	if(m_pAlleleStacks == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddAlleleStack: Memory allocation of %lld bytes failed", (int64_t)memreq);
		return(0);
		}
#else
	m_pAlleleStacks = (tsAlleleStack*)mmap(NULL, (size_t)memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	if(m_pAlleleStacks == MAP_FAILED)
		{
		m_pAlleleStacks = NULL;
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddAlleleStack: Memory allocation of %lld bytes through mmap()  failed", (int64_t)memreq, strerror(errno));
		return(0);
		}
#endif
	m_AllocdAlleleStacksMem = memreq;
	m_AllocdAlleleStacks = cAllocAlleleStacks;
	m_UsedAlleleStacks = 0;
	}
else
		// needing to allocate more memory?
	if((m_UsedAlleleStacks) >= m_AllocdAlleleStacks)
		{
		ToAllocdAlleleStacks = m_UsedAlleleStacks + cReallocAlleleStacks;
		size_t memreq = (size_t)ToAllocdAlleleStacks * sizeof(tsAlleleStack);
#ifdef _WIN32
		pAlleleStack = (tsAlleleStack*)realloc(m_pAlleleStacks, memreq);
		if(pAlleleStack == NULL)
			{
#else
		pAlleleStack = (tsAlleleStack*)mremap(m_pAlleleStacks, m_AllocdAlleleStacksMem, memreq, MREMAP_MAYMOVE);
		if(pAlleleStack == MAP_FAILED)
			{
#endif
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddAlleleStack: Memory reallocation to %lld bytes failed - %s", memreq, strerror(errno));
			return(0);
			}
		m_pAlleleStacks = pAlleleStack;
		m_AllocdAlleleStacksMem = memreq;
		m_AllocdAlleleStacks = ToAllocdAlleleStacks;
		}
pAlleleStack = &m_pAlleleStacks[m_UsedAlleleStacks++];
if(pInitAlleleStack != NULL)
	*pAlleleStack = *pInitAlleleStack;
else
	memset(pAlleleStack, 0, sizeof(tsAlleleStack));
ReleaseSerialise();
return(m_UsedAlleleStacks);
}


#ifdef _WIN32
unsigned __stdcall WorkerPBAInstance(void * pThreadPars)
#else
void *WorkerPBAInstance(void * pThreadPars)
#endif
{
int Rslt;
tsWorkerInstance *pPars = (tsWorkerInstance *)pThreadPars;			// makes it easier not having to deal with casts!
CCallHaplotypes *pWorkerInstance = (CCallHaplotypes *)pPars->pThis;

Rslt = pWorkerInstance->ProcWorkerThread(pPars);
pPars->Rslt = Rslt;
#ifdef _WIN32
_endthreadex(0);
return(eBSFSuccess);
#else
pthread_exit(NULL);
#endif
}

// initialise and start pool of worker threads
int
CCallHaplotypes::StartWorkerThreads(uint32_t NumThreads,		// there are this many threads in pool
									uint32_t NumChroms)			// processing this number of chromosomes
{
uint32_t MaxWait;
uint32_t StartQuerySeqIdx;
uint32_t ThreadIdx;
uint32_t StartedInstances;
tsWorkerInstance *pThreadPar;

if(NumThreads > m_TotWorkQueueEls)
	NumThreads = m_TotWorkQueueEls;

StartQuerySeqIdx = 0;
pThreadPar = m_WorkerInstances;
for(ThreadIdx = 0; ThreadIdx < NumThreads; ThreadIdx++,pThreadPar++)
	{
	memset(pThreadPar,0,sizeof(tsWorkerInstance));
#ifdef _WIN32
	pThreadPar->threadHandle = NULL;
#else
	pThreadPar->threadID = 0;
#endif
	}


m_NumWorkerInsts = 0;
m_CompletedWorkerInsts = 0;
m_ReqTerminate = 0;
pThreadPar = m_WorkerInstances;
for (ThreadIdx = 1; ThreadIdx <= NumThreads; ThreadIdx++, pThreadPar++)
	{
	memset(pThreadPar,0,sizeof(tsWorkerInstance));
	pThreadPar->ThreadIdx = ThreadIdx;
	pThreadPar->pThis = this;
#ifdef _WIN32
	pThreadPar->threadHandle = (HANDLE)_beginthreadex(NULL, 0x0fffff, WorkerPBAInstance, pThreadPar, 0, &pThreadPar->threadID);
#else
	pThreadPar->threadRslt = pthread_create(&pThreadPar->threadID, NULL, WorkerPBAInstance, pThreadPar);
#endif
	}

// allow threads a few seconds to startup
// check if all threads did actually startup; if they did so then m_NumWorkerInsts will have been incremented to NumInstances
MaxWait = 120;		// allowing at most 120 secs for threads to startup
do {
#ifdef WIN32
	Sleep(5000);
#else
	sleep(5);
#endif
	AcquireSerialise();
	StartedInstances = m_NumWorkerInsts;
	ReleaseSerialise();
	MaxWait -= 1;
	}
while(StartedInstances != NumThreads && MaxWait > 0);
if(StartedInstances != NumThreads)
	{
	TerminateWorkerThreads();
	StartedInstances = 0;
	}
return(StartedInstances);
}

// stop all threads in worker pool
int
CCallHaplotypes::TerminateWorkerThreads(int WaitSecs)				// allow at most this many seconds before force terminating worker pool threads
{
int NumForceTerminated;
uint32_t Idx;
uint32_t StartedInstances; 
tsWorkerInstance *pThreadPar;
time_t Then;
time_t Now;

AcquireSerialise();
StartedInstances = m_NumWorkerInsts;
ReleaseSerialise();
if(StartedInstances == 0)
	{
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "TerminateWorkerThreads: No worker threads to terminate");
	return(0);
	}

// request all worker threads to self terminate
AcquireSerialise();
m_ReqTerminate = 1;
ReleaseSerialise();
gDiagnostics.DiagOut(eDLInfo, gszProcName, "TerminateWorkerThreads: Requesting %u worker threads to terminate",StartedInstances);
Then = time(NULL) + WaitSecs;
NumForceTerminated = 0;
pThreadPar = m_WorkerInstances;
for(Idx = 0; Idx < StartedInstances; Idx++, pThreadPar += 1)
	{
	Now = time(NULL);
	if(Now >= Then)
		Now = 1;
	else
		Now = Then - Now;

#ifdef WIN32
	if(pThreadPar->threadHandle != NULL)
		{
		if(WAIT_TIMEOUT == WaitForSingleObject(pThreadPar->threadHandle, (uint32_t)Now * 1000))
			{
			NumForceTerminated += 1;
			TerminateThread(pThreadPar->threadHandle,0);
			}
		pThreadPar->threadHandle = NULL;
		}
#else
	if(pThreadPar->threadID != 0)
		{
		struct timespec ts;
		int JoinRlt;
		void *pExitRslt;
		clock_gettime(CLOCK_REALTIME, &ts);
		ts.tv_sec += Now;
		if ((JoinRlt = pthread_timedjoin_np(pThreadPar->threadID, &pExitRslt, &ts)) != 0)
			{
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "TerminateWorkerThreads: Force terminating thread %u, pthread_timedjoin_np() returned %d",pThreadPar->ThreadIdx,JoinRlt);
			NumForceTerminated += 1;
			pthread_cancel(pThreadPar->threadID);	
			pthread_join(pThreadPar->threadID, NULL);
			}
		pThreadPar->threadID = 0;
		}
#endif
	}

pThreadPar = m_WorkerInstances;
for(Idx = 0; Idx < StartedInstances; Idx++, pThreadPar += 1)
	memset(pThreadPar,0,sizeof(tsWorkerInstance));

m_ReqTerminate = 0;	
return(NumForceTerminated);
}


int
CCallHaplotypes::ProcWorkerThread(tsWorkerInstance *pThreadPar)	// worker thread parameters
{
int Rslt;
uint32_t NumQueueElsProcessed;
tsWorkQueueEl *pWorkQueueEl;
uint32_t ReqTerminate;
// this thread has started, one more worker thread
AcquireSerialise();
m_NumWorkerInsts++;
ReleaseSerialise();
Rslt = 0;
while(Rslt >= 0) {
	// check if requested to terminate
	AcquireSerialise();
	ReqTerminate = m_ReqTerminate;
	ReleaseSerialise();
	if(ReqTerminate)
		break;

	AcquireSerialise();
	NumQueueElsProcessed = m_NumQueueElsProcessed;
	if(NumQueueElsProcessed == m_TotWorkQueueEls)
		{
		ReleaseSerialise();
		break;
		}
	m_NumQueueElsProcessed++;
	ReleaseSerialise();
	pWorkQueueEl = &m_pWorkQueueEls[NumQueueElsProcessed];
	Rslt = GenChromAlleleStacks(pWorkQueueEl->ChromID, pWorkQueueEl->ChromLen, pWorkQueueEl->StartLoci, pWorkQueueEl->MaxNumLoci, pWorkQueueEl->NumFndrs, pWorkQueueEl->pFounderPBAs, pWorkQueueEl->pMskPBA);
	}

AcquireSerialise();
m_CompletedWorkerInsts++;
ReleaseSerialise();
return(Rslt);
}

bool
CCallHaplotypes::WaitAlignments(int WaitSecs)	// allow at most this many seconds for pool of worker threads to complete PBA scoring
{
int64_t WaitMS;
WaitMS = (int64_t)WaitSecs * 1000;

uint32_t CompletedInstances;
uint32_t NumWorkerInsts;
do {
	CUtility::SleepMillisecs(100);
	AcquireSerialise();
	CompletedInstances = m_CompletedWorkerInsts;
	NumWorkerInsts = m_NumWorkerInsts;
	ReleaseSerialise();
	WaitMS -= 100;
	}
while(CompletedInstances < NumWorkerInsts && WaitMS > 0);

return(CompletedInstances == NumWorkerInsts ? true : false);
}


int				// < 0 if errors, otherwise success
CCallHaplotypes::GenChromAlleleStacks(uint32_t ChromID,	// processing is for this chromosome
										  uint32_t ChromSize,		// chromosome is this size
										  uint32_t Loci,			// processing for allele stacks starting from this loci
										  uint32_t MaxNumLoci,		// processing for this maximum number of loci
										  uint32_t NumFndrs,		// number of founders to be processed 1st PBA at *pPBAs[0]
										  uint8_t* pFounderPBAs[],	// pPBAs[] pts to chromosome PBAs, for each of the chromosome founder PBAs
										  uint8_t* pMskPBA)			// pts to optional masking PBA, scoring only for segments contained in mask which are non-zero and where progeny and founder PBAs also have an allele.
																	// enables a founder to be processed as if a progeny, restricting scoring Kmers to be same as if a progeny

{
	uint32_t FounderIdx;
	tsAlleleStack AlleleStack;
	uint32_t NumAlleleStacks;
	bool bFndrAllele;
	uint32_t PrevAlleleIdx;
	uint32_t NumFndrs2Proc;
	uint8_t* pFndrLoci;
	uint32_t AlleleIdx;
	uint8_t AlleleMsk;
	uint32_t AlleleLoci;
	uint32_t FndrIdx;
	uint32_t DivAlleles;

	uint32_t LociFounders;
	uint8_t FndrBaseAllele;
	ts256Bits Fndrs2Proc;

	if(NumFndrs < 1 || NumFndrs > 256)
		return(-1);

	// check that at least 1 founders actually has a PBAs
	memset(&Fndrs2Proc, 0, sizeof(Fndrs2Proc));
	NumFndrs2Proc = 0;
	for(FounderIdx = 0; FounderIdx < NumFndrs; FounderIdx++)
	{
		if(!(m_Fndrs2Proc[FounderIdx] & 0x01))		// not interested as founder not marked for scoring?
			continue;
		if((pFndrLoci = pFounderPBAs[FounderIdx]) == NULL)	// can't score if no PBA for this founder
			continue;
		Bits256Set(NumFndrs2Proc, Fndrs2Proc);
		NumFndrs2Proc++;							// can use this founder for generating an allelestack
	}
	if(NumFndrs2Proc < 1)							// must have at least 1 founder
		return(-2);

	uint32_t EndLoci = min(Loci + MaxNumLoci, ChromSize) - 1;

	memset(&AlleleStack, 0, sizeof(AlleleStack));
	AlleleStack.ChromID = ChromID;
	AlleleStack.NumFndrs = NumFndrs;
	AlleleStack.NumProcFndrs = NumFndrs2Proc;
	memcpy(&AlleleStack.ProcFndrs, &Fndrs2Proc, sizeof(ts256Bits));

	bFndrAllele = false;
	NumAlleleStacks = 0;
	for(AlleleLoci = Loci; AlleleLoci <= EndLoci; AlleleLoci++)
	{
		if(pMskPBA != NULL && pMskPBA[AlleleLoci] == 0)	// skipping over any loci in control which has no PBA
			continue;

	// build stack of alleles containing founder alleles at AlleleLoci which are at least level 2
	// to be accepted as an allele stack there must be at least one founder having a level 2 allele unique to this founder only
		AlleleStack.Loci = AlleleLoci;
		memset(AlleleStack.NumAlleleFndrs, 0, sizeof(AlleleStack.NumAlleleFndrs));
		memset(AlleleStack.Alleles, 0, sizeof(AlleleStack.Alleles));
		LociFounders = 0;
		PrevAlleleIdx = 0;
		for(FndrIdx = 0; FndrIdx < NumFndrs; FndrIdx++)
		{
			if(!(m_Fndrs2Proc[FndrIdx] & 0x01))	// skip this founder?
				continue;
			bFndrAllele = false;
			pFndrLoci = pFounderPBAs[FndrIdx] + AlleleLoci;
			if(*pFndrLoci == 0)					// all founders must have alignment alleles at the AlleleLoci
				break;							// founder has no allele - can't have had an alignment to reference at AlleleLoci

			pFndrLoci = pFounderPBAs[FndrIdx] + AlleleLoci;
			AlleleMsk = 0x03;					// differentiation is at the allele level
			for(AlleleIdx = 0; AlleleIdx < 4; AlleleIdx++, AlleleMsk <<= 2)
			{
				FndrBaseAllele = (*pFndrLoci & AlleleMsk) >> (AlleleIdx * 2);
				if(FndrBaseAllele == 0)		// no allele
					continue;				// try next allele
				if(bFndrAllele)				// founders must have a single unique allele
				{
					bFndrAllele = false;
					break;
				}
				bFndrAllele = true;			// accepting this founder as having at least one allele, if more alleles then founder loci treated as having no alleles
				AlleleStack.NumAlleleFndrs[AlleleIdx]++;
				Bits256Set(FndrIdx, AlleleStack.Alleles[AlleleIdx]);
			}
			if(!bFndrAllele)				// if no alleles of required level for this founder then don't bother with other founders
				break;
		}

		if(!bFndrAllele)				// if no alleles then process next loci
			continue;

			// is there diversity in alleles between all founders - need diversity in order to differentiate founder alleles in progeny
		for(DivAlleles = AlleleIdx = 0; AlleleIdx < 4; AlleleIdx++)
		{
			if(AlleleStack.NumAlleleFndrs[AlleleIdx] > 0)
				DivAlleles++;

		}
		if(DivAlleles <= 1)						// if no uniques then non-informative...
			continue;

		AddAlleleStack(&AlleleStack);
		NumAlleleStacks++;
	}
	return(NumAlleleStacks);
}


int				// returns number of allele stacks generated
CCallHaplotypes::AlignAlleleStacks(uint32_t NumFndrs,			// number of founders to be processed against
									uint32_t MaxNumLoci)		// work queue items specify this number of loci for processing by threads
{
uint32_t FounderID;
tsReadsetMetadata *pReadsetMetadata;
tsChromMetadata *pChromMetadata;
tsWorkQueueEl *pWorkQueueEl;
uint32_t CurChromMetadataIdx;
uint8_t *pPBAs[cMaxFounderReadsets+1];							// additional is to allow for the progeny readset PBAs
uint8_t *pMskPBA;
uint32_t ChromIdx;
uint32_t NumUnits;			// NumUnits is the total number of work units to be distributed over available threads for processing chromosomes
uint32_t UnitSize;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"AlignAlleleStacks: Starting to generate allele stacks over %d founders ",NumFndrs);
pReadsetMetadata = &m_Readsets[0];		// founders were loaded first so 1st founder will be here
CurChromMetadataIdx = pReadsetMetadata->StartChromMetadataIdx;
if(m_pWorkQueueEls != NULL)				// ensuring only one work queue exists
	delete []m_pWorkQueueEls;
m_pWorkQueueEls = NULL;
m_AllocWorkQueueEls = 0;
pMskPBA = NULL;

// determine maximal number of work queue elements required
CurChromMetadataIdx = pReadsetMetadata->StartChromMetadataIdx;
for(NumUnits = ChromIdx = 0; ChromIdx < pReadsetMetadata->NumChroms && CurChromMetadataIdx != 0; ChromIdx++)
	{
	pChromMetadata = &m_pChromMetadata[CurChromMetadataIdx - 1];
	if(pChromMetadata->ChromLen >= MaxNumLoci)
		NumUnits += (pChromMetadata->ChromLen + MaxNumLoci - 1) / MaxNumLoci;
	NumUnits += 1; // safety - accounting for potential partial AccumBinSizes
	CurChromMetadataIdx = pChromMetadata->NxtChromMetadataIdx;
	}
if(NumUnits == 0)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "AlignAlleleStacks: Number of thread work units is 0");
	return(0);
	}
// sizing work queue to contain a maximum of NumUnits elements
m_AllocWorkQueueEls = NumUnits;
m_pWorkQueueEls = new tsWorkQueueEl[m_AllocWorkQueueEls];
memset(m_pWorkQueueEls, 0, sizeof(tsWorkQueueEl) * m_AllocWorkQueueEls);
m_TotWorkQueueEls = 0;
m_NumQueueElsProcessed = 0;
pWorkQueueEl = m_pWorkQueueEls;

// iterate along PBAs which are common to all founders, plus mask readset if present, and then initialise elements of work to be undertaken by each thread
CurChromMetadataIdx = pReadsetMetadata->StartChromMetadataIdx;
for(ChromIdx = 0; ChromIdx < pReadsetMetadata->NumChroms && CurChromMetadataIdx != 0; ChromIdx++)
	{
	pChromMetadata = &m_pChromMetadata[CurChromMetadataIdx - 1];

	// all founders, and mask readset if present, must have PBAs for same chromosome
	for(FounderID = 1; FounderID <= NumFndrs; FounderID++)
		{
		if((pPBAs[FounderID-1] = LocatePBAfor(FounderID, pChromMetadata->ChromID)) == NULL)
			{
			char *pszChrom = LocateChrom(pChromMetadata->ChromID);
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "AlignAlleleStacks: No PBA for chromosome '%s' in at least one founder '%s', skipping this chromosome", pszChrom, LocateReadset(FounderID));
			break;
			}
		}

	CurChromMetadataIdx = pChromMetadata->NxtChromMetadataIdx;		// ready for next iteration over the chromosomes
	if(FounderID <= NumFndrs)	// true only if PBAs not common to progeny and all founders, including mask readset
		continue;

	if(m_MaskReadsetID == 0)				// if masking then chromosome PBA for mask also required
		pMskPBA = NULL;
	else
		if((pMskPBA = LocatePBAfor(m_MaskReadsetID, pChromMetadata->ChromID)) == NULL)
			continue;

	// what is the maximal sized unit that a thread would be expected to process in this chromosome
	if(pChromMetadata->ChromLen < MaxNumLoci)	// too much overhead if too many threads processing shorter chromosomes
		UnitSize = pChromMetadata->ChromLen; // a single thread to process this small chromosome
	else
		UnitSize = MaxNumLoci;

	uint32_t StartLoci = 0;
	while(StartLoci < pChromMetadata->ChromLen)
		{
		if(m_TotWorkQueueEls == m_AllocWorkQueueEls)
			break;
		pWorkQueueEl->NumFndrs = NumFndrs;
		pWorkQueueEl->ChromID = pChromMetadata->ChromID;
		pWorkQueueEl->StartLoci = StartLoci;
		pWorkQueueEl->ChromLen = pChromMetadata->ChromLen;
		pWorkQueueEl->MaxNumLoci = min(MaxNumLoci, pChromMetadata->ChromLen - StartLoci);
		pWorkQueueEl->pMskPBA = pMskPBA;
		memcpy(pWorkQueueEl->pFounderPBAs,pPBAs,((size_t)NumFndrs * sizeof(uint8_t *)));
		m_TotWorkQueueEls++;
		pWorkQueueEl++;
		StartLoci += UnitSize;
		}
	}

// startup threads
StartWorkerThreads(m_NumThreads,		// there are this many threads in pool
					pReadsetMetadata->NumChroms);		// processing this number of chromosomes
while(!WaitAlignments(60))
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"AlignAlleleStacks: Generating allele stacks over %d founders",NumFndrs);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"AlignAlleleStacks: Completed generating %u allele stacks over %d founders", m_UsedAlleleStacks,NumFndrs);
if(m_pWorkQueueEls != NULL)
	{
	delete []m_pWorkQueueEls;
	m_pWorkQueueEls = NULL;	
	m_NumQueueElsProcessed = 0;
	}
return(m_UsedAlleleStacks);
}


int
CCallHaplotypes::GenAlleleStacks(uint32_t NumFndrs,			// number of founders to be processed against
								 uint32_t MaxNumLoci)		// work queue items specify this number of loci for processing by threads
{
int Rslt;
if((Rslt = AlignAlleleStacks( NumFndrs,		// number of founders to be processed against
		  MaxNumLoci)) < 1)	// comparing ReadsetID's PBA against founder PBAs then work queue items specify this number of loci for processing by threads
	return(Rslt);
m_mtqsort.SetMaxThreads(m_NumThreads);
m_mtqsort.qsort(m_pAlleleStacks,(int64_t)m_UsedAlleleStacks,sizeof(tsAlleleStack), SortAlleleStacks);
return(m_UsedAlleleStacks);
}


inline void
CCallHaplotypes::Bits256Shl1(ts256Bits &Bits256)	// SHL 
{
uint64_t MSB;
MSB = Bits256.Bits[2] >> 63; // MSB of word 2 into LSB of word 3 post SHL of word 3
Bits256.Bits[3] <<= 1;
Bits256.Bits[3] |= MSB & 0x01;
MSB = Bits256.Bits[1] >> 63; // MSB of word 1 into LSB of word 2 post SHL of word 2
Bits256.Bits[2] <<= 1;
Bits256.Bits[2] |= MSB & 0x01;
MSB = Bits256.Bits[0] >> 63; // MSB of word 0 into LSB of word 1 post SHL of word 1
Bits256.Bits[1] <<= 1;
Bits256.Bits[1] |= MSB & 0x01;
Bits256.Bits[0] <<= 1;		  // LSB of word 0 will be set to 0 post the SHL
}

inline void
CCallHaplotypes::Bits256Set(uint8_t Bit,		// bit to set, range 0..255
			ts256Bits & Bits256)
{
	Bits256.Bits[Bit / 64] |= ((uint64_t)0x01 << (Bit % 64));
}

inline void
CCallHaplotypes::Bits256Reset(uint8_t Bit,		// bit to reset, range 0..255
		   ts256Bits & Bits256)
{
	Bits256.Bits[Bit / 64] &= ~((uint64_t)0x01 << (Bit % 64));
}


inline bool
CCallHaplotypes::Bits256Test(uint8_t Bit,		// bit to test, range 0..255
			ts256Bits & Bits256)
{
return(Bits256.Bits[Bit / 64] & ((uint64_t)0x01 << (Bit % 64)) ? true : false);
}

inline void
CCallHaplotypes::Bits256Initialise(bool Set,			// if true then initialse all bits as set, otherwise initialise all bits as reset
				ts256Bits & Bits256)
{
uint64_t InitWith;
if(Set)
	InitWith = (uint64_t)(-1);
else
	InitWith = 0;
Bits256.Bits[0] = InitWith;
Bits256.Bits[1] = InitWith;
Bits256.Bits[2] = InitWith;
Bits256.Bits[3] = InitWith;
}

inline uint32_t
CCallHaplotypes::Bits256Count(ts256Bits& Bits256)		// count number of set bits
{
uint32_t Count = 0;
uint64_t *pWord = Bits256.Bits;
uint64_t Bit;
for(uint32_t WordIdx = 0; WordIdx < 4; WordIdx++, pWord++)
	{
	for(Bit = 0x01; Bit != 0; Bit <<= 1)
		if(*pWord & Bit)
			Count++;
	}
return(Count);
}

inline uint32_t
CCallHaplotypes::Bits256Combine(ts256Bits &Bits256A, ts256Bits& Bits256B)		// combine (effective |= ) bits in Bits256A with Bits256B with Bits256A updated 
{
uint32_t Count = 0;
uint64_t* pWordA = Bits256A.Bits;
uint64_t* pWordB = Bits256B.Bits;
uint64_t Bit;
for(uint32_t WordIdx = 0; WordIdx < 4; WordIdx++, pWordA++, pWordB++)
	{
	*pWordA |= *pWordB;
	for(Bit = 0x01; Bit != 0; Bit <<= 1)
		if(*pWordA & Bit)
			Count++;
	}
return(Count);
}


double
CCallHaplotypes::ChiSqr(int NumRows, int32_t* pExp,int32_t *pObs)
{
int64_t Diff;
double Sum = 0;
if(*pExp == 0)
	return(0);
// iterate and sum over all rows
for(int Row = 0; Row < NumRows; Row++, pExp++, pObs++)
	{
	Diff = (int64_t)*pExp - *pObs;
	Diff *= Diff;
	Sum += (double)Diff / *pExp == 0 ? 0.0001 : (double)*pExp;	// avoiding divide by 0 issues!
	}
return(Sum);
}


int32_t		// returned chrom identifier, < 1 if unable to accept this chromosome name
CCallHaplotypes::AddChrom(char* pszChrom) // associate unique identifier with this chromosome name
{
int32_t ChromNameIdx;
int ChromNameLen;
char *pszLAname;

// with any luck the sequence name will be same as the last accessed
if((pszLAname = LocateChrom(m_LAChromNameID)) != NULL)
	if(!stricmp(pszChrom,pszLAname))
		return(m_LAChromNameID);

// iterate over all known chroms in case this chrom to add is a duplicate
for(ChromNameIdx = 0; ChromNameIdx < m_NumChromNames; ChromNameIdx++)
	if(!stricmp(pszChrom, &m_szChromNames[m_szChromIdx[ChromNameIdx]]))
		{
		m_LAChromNameID = ChromNameIdx + 1;
		return(m_LAChromNameID);
		}

// chrom is not a duplicate
ChromNameLen = (int)strlen(pszChrom);
if((m_NxtszChromIdx + ChromNameLen + 1) > (int)sizeof(m_szChromNames))
	return(0);
if(m_NumChromNames == cMaxChromNames)
	return(0);

m_szChromIdx[m_NumChromNames++] = m_NxtszChromIdx;
strcpy(&m_szChromNames[m_NxtszChromIdx], pszChrom);
m_NxtszChromIdx += ChromNameLen + 1;
m_LAChromNameID = m_NumChromNames;
return(m_LAChromNameID);
}


int32_t		// returned chrom identifier, < 1 if unable to locate this chromosome name
CCallHaplotypes::LocateChrom(char* pszChrom) // return unique identifier associated with this chromosome name
{
int32_t ChromNameIdx;
char *pszLAname;

// with any luck the sequence name will be same as the last accessed
if((pszLAname = LocateChrom(m_LAChromNameID)) != NULL)
	if(!stricmp(pszChrom,pszLAname))
		return(m_LAChromNameID);

// iterate over all known chroms
for(ChromNameIdx = 0; ChromNameIdx < m_NumChromNames; ChromNameIdx++)
	if(!stricmp(pszChrom, &m_szChromNames[m_szChromIdx[ChromNameIdx]]))
		{
		m_LAChromNameID = ChromNameIdx + 1;
		return(m_LAChromNameID);
		}
return(0);
}

char* 
CCallHaplotypes::LocateChrom(int32_t ChromID)
{
if(ChromID < 1 || ChromID > m_NumChromNames)
	return(NULL);
return(&m_szChromNames[m_szChromIdx[ChromID-1]]);
}

uint8_t 
CCallHaplotypes::LocateReadsetChromLociAlleles(int32_t ReadsetID,	// return alleles for this readset 
									  int32_t ChromID,		// on this chromosome
									  uint32_t Loci)		// at this loci
{
uint8_t *pPBA;
if((pPBA=LocatePBAfor(ReadsetID,ChromID))==NULL)
	return(0);
return(pPBA[Loci]);
}

// NOTE: Readsets are checked for uniqueness as readsets must be unique within a given readetset type
int32_t		// returned readset identifier, < 1 if unable to accept this readset name
CCallHaplotypes::AddReadset(char* pszReadset, // associate unique identifier with this readset name
							uint8_t ReadsetType)	// 0: founder, 1: progeny, 2: control
{
int32_t ReadsetNameIdx;
int ReadsetNameLen;
char Type;
char *pszLAname;
Type = '0' + (char)ReadsetType;

// with any luck the sequence name will be same as the last accessed
if((pszLAname = LocateReadset(m_LAReadsetNameID)) != NULL)
	{
	if(pszLAname[-1] == Type &&	// prefixed with ReadsetType
		!stricmp(pszReadset,pszLAname))
			return(0); // non-unique within ReadsetType!
	}
	
// iterate over all known readsets in case this readset to add is a duplicate
for(ReadsetNameIdx = 0; ReadsetNameIdx < m_NumReadsetNames; ReadsetNameIdx++)
	{
	pszLAname = &m_szReadsetNames[m_szReadsetIdx[ReadsetNameIdx]];
	if(*pszLAname++ != Type)
		continue;
	if(!stricmp(pszReadset, pszLAname))
		{
		m_LAReadsetNameID = ReadsetNameIdx + 1;
		return(0); // non-unique within ReadsetType!
		}
	}

// not a duplicate
ReadsetNameLen = (int)strlen(pszReadset);
if((m_NxtszReadsetIdx + ReadsetNameLen + 2) > (int)sizeof(m_szReadsetNames))
	return(0);		// no more space to hold reeadset name, treating as if a non-unique!
if(m_NumReadsetNames == (cMaxProgenyReadsets + cMaxFounderReadsets))
	return(0);		// unable to hold any more readsets, treating as if a non-unique!

m_szReadsetIdx[m_NumReadsetNames++] = m_NxtszReadsetIdx;
pszLAname = &m_szReadsetNames[m_NxtszReadsetIdx];
*pszLAname++ = Type;
strcpy(pszLAname, pszReadset);
m_NxtszReadsetIdx += ReadsetNameLen + 2;
m_LAReadsetNameID = m_NumReadsetNames;
return(m_LAReadsetNameID);
}


int32_t		// returned Readset identifier, < 1 if unable to locate this Readset name
CCallHaplotypes::LocateReadset(char* pszReadset, // return unique identifier associated with this Readset name
							   uint8_t ReadsetType)	// 0: founder, 1: progeny, 2: control
{
int32_t ReadsetNameIdx;
char Type;
char *pszLAReadset;

Type = '0'+ (char)ReadsetType;
// with any luck the Readset name will be same as the last accessed
if((pszLAReadset = LocateReadset(m_LAReadsetNameID)) != NULL)
	{
	if(pszLAReadset[-1] == Type &&	// prefixed with ReadsetType
		!stricmp(pszReadset,pszLAReadset))
			return(m_LAReadsetNameID);
	}

// iterate over all known Readsets
for(ReadsetNameIdx = 0; ReadsetNameIdx < m_NumReadsetNames; ReadsetNameIdx++)
	{
	pszLAReadset = &m_szReadsetNames[m_szReadsetIdx[ReadsetNameIdx]];
	if(*pszLAReadset++ != Type)		// readset must be that of ReadsetType
		continue;
	if(!stricmp(pszReadset, pszLAReadset))
		{
		m_LAReadsetNameID = ReadsetNameIdx + 1;
		return(m_LAReadsetNameID);
		}
	}
return(0);
}

char* 
CCallHaplotypes::LocateReadset(int32_t ReadsetID)
{
int Idx;
Idx = ReadsetID & 0x0fffffff;			// mask out any potential ReadsetType
if(Idx < 1 || Idx > m_NumReadsetNames)
	return(NULL);
return(&(m_szReadsetNames[m_szReadsetIdx[Idx -1]+1])); // skipping lead char which is the ReadsetType
}

// sorting by ChromID.StartLoci ascending
int
CCallHaplotypes::SortAlleleStacks(const void* arg1, const void* arg2)
{
tsAlleleStack* pEl1 = (tsAlleleStack*)arg1;
tsAlleleStack* pEl2 = (tsAlleleStack*)arg2;
if(pEl1->ChromID < pEl2->ChromID)
	return(-1);
if(pEl1->ChromID > pEl2->ChromID)
	return(1);
if(pEl1->Loci < pEl2->Loci)
	return(-1);
if(pEl1->Loci > pEl2->Loci)
	return(1);
return(0);
}

// sorting by ReadsetID.ChromID.StartLoci ascending
int
CCallHaplotypes::SortProgenyFndrAligns(const void* arg1, const void* arg2)
{
tsProgenyFndrAligns* pEl1 = (tsProgenyFndrAligns*)arg1;
tsProgenyFndrAligns* pEl2 = (tsProgenyFndrAligns*)arg2;
if(pEl1->ReadsetID < pEl2->ReadsetID)
	return(-1);
if(pEl1->ReadsetID > pEl2->ReadsetID)
return(1);
if(pEl1->ChromID < pEl2->ChromID)
	return(-1);
if(pEl1->ChromID > pEl2->ChromID)
	return(1);
if(pEl1->Loci < pEl2->Loci)
	return(-1);
if(pEl1->Loci > pEl2->Loci)
	return(1);
return(0);
}

// sorting by ChromID.StartLoci.Readset ascending
int
CCallHaplotypes::SortProgenyChromLociReadset(const void* arg1, const void* arg2)
{
tsProgenyFndrAligns* pEl1 = (tsProgenyFndrAligns*)arg1;
tsProgenyFndrAligns* pEl2 = (tsProgenyFndrAligns*)arg2;

if(pEl1->ChromID < pEl2->ChromID)
	return(-1);
if(pEl1->ChromID > pEl2->ChromID)
	return(1);
if(pEl1->Loci < pEl2->Loci)
	return(-1);
if(pEl1->Loci > pEl2->Loci)
	return(1);
if(pEl1->ReadsetID < pEl2->ReadsetID)
	return(-1);
if(pEl1->ReadsetID > pEl2->ReadsetID)
	return(1);
return(0);
}


