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
#include "CWIGutils.h"

int
Process(eWIGuMode PMode,	// processing mode: eWIGu2DEseq WIG to DEseq, eWIGu2TODEseq WIG to turnover DEseq
		int32_t NumInputFiles,		// number of input WIG file specs
		char* pszInputFiles[],		// names of input WIG files (wildcards allowed)
		char* pszOutFile,			// output to this file
		char* pszChromFile,			// BED file containing reference assembly chromosome names and sizes
		char* pszROIFile,			// BED file containing regions/genes of interest
		int NumIncludeChroms,		// number of chromosome regular expressions to include
		char *pszIncludeChroms[],	// array of include chromosome regular expressions
		int NumExcludeChroms,		// number of chromosome expressions to exclude
		char *pszExcludeChroms[],	// array of exclude chromosome regular expressions
		int32_t NumThreads);		// number of worker threads to use

#ifdef _WIN32
int wigutils(int argc, char* argv[])
{
	// determine my process name
_splitpath(argv[0], nullptr, nullptr, gszProcName, nullptr);
#else
int
wigutils(int argc, char** argv)
{
	// determine my process name
CUtility::splitpath((char*)argv[0], nullptr, gszProcName);
#endif
int iFileLogLevel;			// level of file diagnostics
int iScreenLogLevel;		// level of file diagnostics
char szLogFile[_MAX_PATH];	// write diagnostics to this file

int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure
int NumberOfProcessors;		// number of installed CPUs
int NumThreads;				// number of threads (0 defaults to number of CPUs or a maximum of cMaxPBAWorkerThreads)

int Idx;

eWIGuMode PMode;			// processing mode: eWIGu2DEseq WIG to DEseq, eWIGu2TODEseq WIG to turnover DEseq
char szChromFile[_MAX_PATH];	// BED file containing chromosome names and sizes
char szROIFile[_MAX_PATH];	// BED file containing regions of interest
char szOutFile[_MAX_PATH];	// output DEseq written to this file
int NumInputFiles;			// number of input files
char* pszInputFiles[cMaxWildCardFileSpecs];		// names of input files (wildcards allowed)
int NumIncludeChroms;
char* pszIncludeChroms[cMaxIncludeChroms];
int NumExcludeChroms;
char* pszExcludeChroms[cMaxExcludeChroms];

struct arg_lit* help = arg_lit0("h", "help", "print this help and exit");
struct arg_lit* version = arg_lit0("v", "version,ver", "print version information and exit");
struct arg_int* FileLogLevel = arg_int0("f", "FileLogLevel", "<int>", "Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file* LogFile = arg_file0("F", "log", "<file>", "diagnostics log file");
struct arg_int* pmode = arg_int0("m", "mode", "<int>", "processing mode: 0 WIG to DEseq, 2 WIG to turnover DEseq");
struct arg_file* chromfile = arg_file1("c", "chromfile", "<file>", "input BED file containing chromosome names and sizes");
struct arg_file* roifile = arg_file1("C", "roifile", "<file>", "BED file containing regions of interest/genes");

struct arg_str* ExcludeChroms = arg_strn("Z", "chromexclude", "<string>", 0, cMaxExcludeChroms, "high priority - regular expressions defining chromosomes to exclude");
struct arg_str* IncludeChroms = arg_strn("z", "chromeinclude", "<string>", 0, cMaxIncludeChroms, "low priority - regular expressions defining chromosomes to include");
struct arg_file* infiles = arg_filen("i", "infiles", "<file>", 1, cMaxWildCardFileSpecs, "input WIG file(s), wildcards allowed, limit of 200 filespecs supported");

struct arg_file* outfile = arg_file1("o", "out", "<file>", "output to this file");

struct arg_int* threads = arg_int0("T", "threads", "<int>", "number of processing threads 0..4 (defaults to 0 which limits threads to maximum of 2 CPU cores)");

struct arg_end* end = arg_end(200);
void* argtable[] = {help,version,FileLogLevel,LogFile,
					pmode,chromfile,roifile,ExcludeChroms,IncludeChroms,infiles,outfile,
					threads,end};
char** pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc, (char**)argv, &pAllArgs);
if(argerrors >= 0)
argerrors = arg_parse(argerrors, pAllArgs, argtable);

/* special case: '--help' takes precedence over error reporting */
if(help->count > 0)
	{
	printf("\n%s %s %s, Version %s\nOptions ---\n", gszProcName, gpszSubProcess->pszName, gpszSubProcess->pszFullDescr, kit4bversion);
	arg_print_syntax(stdout, argtable, "\n");
	arg_print_glossary(stdout, argtable, "  %-25s %s\n");
	printf("\nNote: Parameters can be entered into a parameter file, one parameter per line.");
	printf("\n      To invoke this parameter file then precede its name with '@'");
	printf("\n      e.g. %s %s @myparams.txt\n", gszProcName, gpszSubProcess->pszName);
	printf("\nPlease report any issues regarding usage of %s at https://github.com/kit4b/issues\n\n", gszProcName);
	return(1);
	}

/* special case: '--version' takes precedence error reporting */
if(version->count > 0)
	{
	printf("\n%s %s Version %s\n", gszProcName, gpszSubProcess->pszName, kit4bversion);
	return(1);
	}

if(!argerrors)
	{
	if(FileLogLevel->count && !LogFile->count)
		{
		printf("\nError: FileLogLevel '-f%d' specified but no logfile '-F<logfile>\n'", FileLogLevel->ival[0]);
		exit(1);
		}

	iScreenLogLevel = iFileLogLevel = FileLogLevel->count ? FileLogLevel->ival[0] : eDLInfo;
	if(iFileLogLevel < eDLNone || iFileLogLevel > eDLDebug)
		{
		printf("\nError: FileLogLevel '-l%d' specified outside of range %d..%d\n", iFileLogLevel, eDLNone, eDLDebug);
		exit(1);
		}

	if(LogFile->count)
		{
		strncpy(szLogFile, LogFile->filename[0], _MAX_PATH);
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

	PMode = pmode->count ? (eWIGuMode)pmode->ival[0] : eWIGu2DEseq;
	if(PMode < eWIGu2DEseq || PMode >= eWIGuPlaceholder)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Unsupported processing mode '-m%d'\n", PMode);
		exit(1);
		}

	// show user current resource limits
#ifndef _WIN32
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Resources: %s", CUtility::ReportResourceLimits());
#endif

#ifdef _WIN32
	SYSTEM_INFO SystemInfo;
	GetSystemInfo(&SystemInfo);
	NumberOfProcessors = SystemInfo.dwNumberOfProcessors;
#else
	NumberOfProcessors = sysconf(_SC_NPROCESSORS_CONF);
#endif
	int MaxAllowedThreads = min((int)cMaxWIGutilityThreads, NumberOfProcessors);	// limit to be at most cMaxWorkerThreads
	if ((NumThreads = threads->count ? threads->ival[0] : MaxAllowedThreads) == 0)
		NumThreads = MaxAllowedThreads;
	if (NumThreads < 0 || NumThreads > MaxAllowedThreads)
		{
		gDiagnostics.DiagOut(eDLWarn, gszProcName, "Warning: Number of threads '-T%d' specified was outside of range %d..%d", NumThreads, 1, MaxAllowedThreads);
		gDiagnostics.DiagOut(eDLWarn, gszProcName, "Warning: Defaulting number of threads to %d", MaxAllowedThreads);
		NumThreads = MaxAllowedThreads;
		}
	if (chromfile->count)
		{
		strncpy(szChromFile, chromfile->filename[0], _MAX_PATH);
		szChromFile[_MAX_PATH-1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(szChromFile);
		}
	else
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "No BED file containing chromosome names and sizes specified");
		exit(1);
		}

	NumIncludeChroms = 0;
	if (IncludeChroms->count)
		{
		for (NumIncludeChroms = Idx = 0; NumIncludeChroms < cMaxIncludeChroms && Idx < IncludeChroms->count; Idx++)
			{
			pszIncludeChroms[Idx] = nullptr;
			if (pszIncludeChroms[NumIncludeChroms] == nullptr)
				pszIncludeChroms[NumIncludeChroms] = new char[_MAX_PATH];
			strncpy(pszIncludeChroms[NumIncludeChroms], IncludeChroms->sval[Idx], _MAX_PATH);
			pszIncludeChroms[NumIncludeChroms][_MAX_PATH - 1] = '\0';
			CUtility::TrimQuotedWhitespcExtd(pszIncludeChroms[NumIncludeChroms]);
			if (pszIncludeChroms[NumIncludeChroms][0] != '\0')
				NumIncludeChroms++;
			}
		}
	

	NumExcludeChroms = 0;
	if (ExcludeChroms->count)
		{
		for (NumExcludeChroms = Idx = 0; NumExcludeChroms < cMaxExcludeChroms && Idx < ExcludeChroms->count; Idx++)
			{
			pszExcludeChroms[Idx] = nullptr;
			if (pszExcludeChroms[NumExcludeChroms] == nullptr)
				pszExcludeChroms[NumExcludeChroms] = new char[_MAX_PATH];
			strncpy(pszExcludeChroms[NumExcludeChroms], ExcludeChroms->sval[Idx], _MAX_PATH);
			pszExcludeChroms[NumExcludeChroms][_MAX_PATH - 1] = '\0';
			CUtility::TrimQuotedWhitespcExtd(pszExcludeChroms[NumExcludeChroms]);
			if (pszExcludeChroms[NumExcludeChroms][0] != '\0')
				NumExcludeChroms++;
			}
		}
	NumInputFiles = 0;
	if (infiles->count)
		{
		for (Idx = 0; NumInputFiles < cMaxWildCardFileSpecs && Idx < infiles->count; Idx++)
			{
			pszInputFiles[Idx] = nullptr;
			if (pszInputFiles[NumInputFiles] == nullptr)
				pszInputFiles[NumInputFiles] = new char[_MAX_PATH];
			strncpy(pszInputFiles[NumInputFiles], infiles->filename[Idx], _MAX_PATH);
			pszInputFiles[NumInputFiles][_MAX_PATH - 1] = '\0';
			CUtility::TrimQuotedWhitespcExtd(pszInputFiles[NumInputFiles]);
			if (pszInputFiles[NumInputFiles][0] != '\0')
				NumInputFiles++;
			}
		}

	if (!NumInputFiles)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: After removal of whitespace, there are no input file(s) specified with '-i<filespec>' option)\n");
		exit(1);
		}

	strcpy(szOutFile, outfile->filename[0]);
	CUtility::TrimQuotedWhitespcExtd(szOutFile);
	if(szOutFile[0] == '\0')
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "No output file specified");
		exit(1);
		}

	szROIFile[0] = '\0';
	if(roifile->count)
		{
		strncpy(szROIFile, roifile->filename[0], _MAX_PATH);
		szROIFile[_MAX_PATH-1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(szROIFile);
		}

	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Processing parameters:");
	const char* pszDescr;
	switch(PMode) {
		case eWIGu2DEseq:
			pszDescr = "generating DEseq file from WIG file";
			break;
		case eWIGu2TODEseq:
			pszDescr = "generating turnover DEseq file from WIG file";
			break;

	}

	gDiagnostics.DiagOutMsgOnly(eDLInfo, "WIG utilities : '%s'", pszDescr);

	if(szROIFile[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "Regions of interest file : '%s'", szROIFile);

	gDiagnostics.DiagOutMsgOnly(eDLInfo, "BED containing chromosome names and sizes : '%s'", szChromFile);

	for (Idx = 0; Idx < NumIncludeChroms; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "reg expressions defining chroms to include: '%s'", pszIncludeChroms[Idx]);
	for (Idx = 0; Idx < NumExcludeChroms; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "reg expressions defining chroms to exclude: '%s'", pszExcludeChroms[Idx]);

	for (Idx = 0; Idx < NumInputFiles; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "Input file spec (%d) : '%s'", Idx+1, pszInputFiles[Idx]);

	gDiagnostics.DiagOutMsgOnly(eDLInfo, "Output to file : '%s'", szOutFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo, "number of threads : %d", NumThreads);
#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = 0;
	Rslt = Process(PMode,	// processing mode: eWIGu2DEseq WIG to DEseq, eWIGu2TODEseq WIG to turnover DEseq
					NumInputFiles,		// number of input WIG file specs
					pszInputFiles,		// names of input WIG files (wildcards allowed)
					szOutFile,			// output to this file
					szChromFile,			// BED file containing reference assembly chromosome names and sizes
					szROIFile,			// BED file containing regions/genes of interest
					NumIncludeChroms,		// number of chromosome regular expressions to include
					pszIncludeChroms,	// array of include chromosome regular expressions
					NumExcludeChroms,		// number of chromosome expressions to exclude
					pszExcludeChroms,	// array of exclude chromosome regular expressions
					NumThreads);		// number of worker threads to use

	Rslt = Rslt >= 0 ? 0 : 1;
	gStopWatch.Stop();

	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Exit code: %d Total processing time: %s", Rslt, gStopWatch.Read());
	exit(Rslt);
	}
else
	{
	printf("\n%s %s %s, Version %s\n", gszProcName, gpszSubProcess->pszName, gpszSubProcess->pszFullDescr, kit4bversion);
	arg_print_errors(stdout, end, gszProcName);
	arg_print_syntax(stdout, argtable, "\nUse '-h' to view option and parameter usage\n");
	exit(1);
	}
}

int
Process(eWIGuMode PMode,	// processing mode: eWIGu2DEseq WIG to DEseq, eWIGu2TODEseq WIG to turnover DEseq
		int32_t NumInputFiles,		// number of input WIG file specs
		char* pszInputFiles[],		// names of input WIG files (wildcards allowed)
		char* pszOutFile,			// output to this file
		char* pszChromFile,			// BED file containing reference assembly chromosome names and sizes
		char* pszROIFile,			// BED file containing regions/genes of interest
		int NumIncludeChroms,		// number of chromosome regular expressions to include
		char *pszIncludeChroms[],	// array of include chromosome regular expressions
		int NumExcludeChroms,		// number of chromosome expressions to exclude
		char *pszExcludeChroms[],	// array of exclude chromosome regular expressions
		int32_t NumThreads)		// number of worker threads to use
{
int Rslt;
CWIGutils *pWIGutils;
if((pWIGutils = new CWIGutils) == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to instantiate instance of CWIGutils");
	return(eBSFerrObj);
	}
Rslt = pWIGutils->Process(PMode,NumInputFiles,pszInputFiles,pszOutFile,pszChromFile,pszROIFile,NumIncludeChroms,pszIncludeChroms,NumExcludeChroms,pszExcludeChroms,NumThreads);
delete pWIGutils;
return(Rslt);
}

CWIGutils::CWIGutils()	// constructor
{
m_pChromMetadata = nullptr;
m_pROIFile = nullptr;
m_pBedFile = nullptr;
m_bMutexesCreated = false;
Reset();
}

CWIGutils::~CWIGutils()	// destructor
{
if(m_pROIFile != nullptr)
	delete m_pROIFile;

if(m_pBedFile != nullptr)
	delete m_pBedFile;

if(m_pChromMetadata != nullptr)
	{
	tsWUChromMetadata* pChromMetadata = m_pChromMetadata;
	for (uint32_t Idx = 0; Idx < m_UsedNumChromMetadata; Idx++, pChromMetadata++)
		{
		if (pChromMetadata->pCnts != nullptr)
#ifdef _WIN32
			free(pChromMetadata->pCnts);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
			if (pChromMetadata->pCnts != MAP_FAILED)
				munmap(pChromMetadata->pCnts, pChromMetadata->ChromLen * sizeof(uint32_t));
#endif
		}

#ifdef _WIN32
	free(m_pChromMetadata);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (m_pChromMetadata != MAP_FAILED)
		munmap(m_pChromMetadata, m_AllocdChromMetadataMem);
#endif
	}
if(m_bMutexesCreated)
	DeleteMutexes();
}

void
CWIGutils::Reset(void)
{
if(m_pBedFile != nullptr)
	{
	delete m_pBedFile;
	m_pBedFile = nullptr;
	}

if(m_pROIFile != nullptr)
	{
	delete m_pROIFile;
	m_pROIFile = nullptr;
	}

if (m_pChromMetadata != nullptr)
	{
	tsWUChromMetadata* pChromMetadata = m_pChromMetadata;
	for (uint32_t Idx = 0; Idx < m_UsedNumChromMetadata; Idx++, pChromMetadata++)
		{
		if (pChromMetadata->pCnts != nullptr)
#ifdef _WIN32
			free(pChromMetadata->pCnts);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
		if (pChromMetadata->pCnts != MAP_FAILED)
				munmap(pChromMetadata->pCnts, pChromMetadata->ChromLen * sizeof(uint32_t));
#endif
		}

#ifdef _WIN32
	free(m_pChromMetadata);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (m_pChromMetadata != MAP_FAILED)
		munmap(m_pChromMetadata, m_AllocdChromMetadataMem);
#endif
	m_pChromMetadata = nullptr;
	}

m_UsedNumChromMetadata = 0;
m_AllocdChromMetadata = 0;
m_AllocdChromMetadataMem = 0;

m_NumReadsetIDs = 0;

m_LAReadsetNameID = 0;
m_NumReadsetNames = 0;
m_NxtszReadsetIdx = 0;
m_szReadsetNames[0] = '\0';

m_LAChromNameID = 0;
m_NumChromNames = 0;
m_NxtszChromIdx = 0;
m_szChromNames[0] = '\0';

m_NumChromSizes = 0;

m_szOutFile[0] = '\0';
if (m_bMutexesCreated)
	DeleteMutexes();
m_bMutexesCreated = false;
m_bMutexesCreated = false;
m_NumThreads = cMaxWIGutilityThreads;
}

int
CWIGutils::CreateMutexes(void)
{
if (m_bMutexesCreated)
	return(eBSFSuccess);

#ifdef _WIN32
if ((m_hSerialiseAccess = CreateMutex(nullptr, false, nullptr)) == nullptr)
	{
#else
if (pthread_mutex_init(&m_hSerialiseAccess, nullptr) != 0)
	{
#endif
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Fatal: unable to create mutex");
	return(eBSFerrInternal);

	}
m_FastSerialise = 0;
m_bMutexesCreated = true;
return(eBSFSuccess);
}

void
CWIGutils::DeleteMutexes(void)
{
	if (!m_bMutexesCreated)
		return;
#ifdef _WIN32
	CloseHandle(m_hSerialiseAccess);
#else
	pthread_mutex_destroy(&m_hSerialiseAccess);
#endif
	m_bMutexesCreated = false;
}

void
CWIGutils::AcquireSerialise(void)
{
	uint32_t WaitRslt;
#ifdef _WIN32
	WaitRslt = WaitForSingleObject(m_hSerialiseAccess, INFINITE);
	if (WaitRslt != WAIT_OBJECT_0)
	{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Fatal:WaitForSingleObject() returned error %u", WaitRslt);
		Reset();
		exit(1);
	}
#else
	pthread_mutex_lock(&m_hSerialiseAccess);
#endif
}

void
CWIGutils::ReleaseSerialise(void)
{
#ifdef _WIN32
	ReleaseMutex(m_hSerialiseAccess);
#else
	pthread_mutex_unlock(&m_hSerialiseAccess);
#endif
}

void
CWIGutils::AcquireFastSerialise(void)
{
	int SpinCnt = 500;
	int BackoffMS = 5;

#ifdef _WIN32
	while (InterlockedCompareExchange(&m_FastSerialise, 1, 0) != 0)
	{
		if (SpinCnt -= 1)
			continue;
		CUtility::SleepMillisecs(BackoffMS);
		SpinCnt = 100;
		if (BackoffMS < 500)
			BackoffMS += 2;
	}
#else
	while (__sync_val_compare_and_swap(&m_FastSerialise, 0, 1) != 0)
	{
		if (SpinCnt -= 1)
			continue;
		CUtility::SleepMillisecs(BackoffMS);
		SpinCnt = 100;
		if (BackoffMS < 500)
			BackoffMS += 2;
	}
#endif
}

void
CWIGutils::ReleaseFastSerialise(void)
{
#ifdef _WIN32
	InterlockedCompareExchange(&m_FastSerialise, 0, 1);
#else
	__sync_val_compare_and_swap(&m_FastSerialise, 1, 0);
#endif
}


uint32_t*
CWIGutils::AllocCnts(uint32_t ChromLen)	// allocate memory to hold at least this many coverage counts
{
	uint32_t* pCnts;
	size_t memreq;
	memreq = (size_t)ChromLen * sizeof(uint32_t);	// no safety margin!
#ifdef _WIN32
	pCnts = (uint32_t*)malloc((size_t)memreq);
	if (pCnts == nullptr)
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AllocCnts Memory allocation of %zd bytes failed", (int64_t)memreq);
#else
	pCnts = (uint32_t*)mmap(nullptr, (size_t)memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	if (pCnts == MAP_FAILED)
	{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AllocCnts: Memory allocation of %zd bytes through mmap()  failed", (int64_t)memreq, strerror(errno));
		pCnts = nullptr;
	}
#endif
	return(pCnts);
}

void
CWIGutils::DeleteAllChromCnts(void) // delete all currently loaded PBAs - all sample chroms
{
if (m_pChromMetadata != nullptr)
	{
	tsWUChromMetadata* pChromMetadata = m_pChromMetadata;
	for (uint32_t Idx = 0; Idx < m_UsedNumChromMetadata; Idx++, pChromMetadata++)
		{
		if (pChromMetadata->pCnts != nullptr)
#ifdef _WIN32
			free(pChromMetadata->pCnts);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
		if (pChromMetadata->pCnts != MAP_FAILED)
				munmap(pChromMetadata->pCnts, pChromMetadata->ChromLen * sizeof(uint32_t));
#endif
		pChromMetadata->pCnts = nullptr; 
		}
	}
}

bool  // false: no previous memory allocated for containing the chromosome PBAs, true: allocation was present and has been deleted
CWIGutils::DeleteSampleChromCnts(uint32_t SampleID,   // Sample identifier
	uint32_t ChromID)    // chrom identifier
{
	tsWUChromMetadata* pChromMetadata;

	// returned pointer to chromosome metadata
	if ((pChromMetadata = LocateChromMetadataFor(SampleID, ChromID)) == nullptr)
		return(false);
	if (pChromMetadata->pCnts == nullptr)
		return(false);
#ifdef _WIN32
	free(pChromMetadata->pCnts);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (pChromMetadata->pCnts != MAP_FAILED)
		munmap(pChromMetadata->pCnts, pChromMetadata->ChromLen * sizeof(uint32_t));
#endif
	pChromMetadata->pCnts = nullptr;
	return(true);
}

// NOTE: Readsets are checked for uniqueness as readsets must be unique within a given readetset type
uint32_t		// returned readset identifier, 0 if unable to accept this readset name
CWIGutils::AddReadset(char* pszReadset, // associate unique identifier with this readset name
	uint8_t ReadsetType)	// 0: founder, 1: progeny, 2: control
{
	uint32_t ReadsetNameIdx;
	int ReadsetNameLen;
	char Type;
	char* pszLAname;
	Type = '0' + (char)ReadsetType;

	// with any luck the sequence name will be same as the last accessed
	if ((pszLAname = LocateReadset(m_LAReadsetNameID)) != nullptr)
	{
		if (pszLAname[-1] == Type &&	// prefixed with ReadsetType
			!stricmp(pszReadset, pszLAname))
			return(0); // non-unique within ReadsetType!
	}

	// iterate over all known readsets in case this readset to add is a duplicate
	for (ReadsetNameIdx = 0; ReadsetNameIdx < m_NumReadsetNames; ReadsetNameIdx++)
	{
		pszLAname = &m_szReadsetNames[m_szReadsetIdx[ReadsetNameIdx]];
		if (*pszLAname++ != Type)
			continue;
		if (!stricmp(pszReadset, pszLAname))
		{
			m_LAReadsetNameID = ReadsetNameIdx + 1;
			return(0); // non-unique within ReadsetType!
		}
	}

	// not a duplicate
	ReadsetNameLen = (int)strlen(pszReadset);
	if ((m_NxtszReadsetIdx + ReadsetNameLen + 2) > (int)sizeof(m_szReadsetNames))
		return(0);		// no more space to hold reeadset name, treating as if a non-unique!
	if (m_NumReadsetNames == cMaxWIGReadsets)
		return(0);		// unable to hold any more readsets, treating as if a non-unique!

	m_szReadsetIdx[m_NumReadsetNames++] = m_NxtszReadsetIdx;
	pszLAname = &m_szReadsetNames[m_NxtszReadsetIdx];
	*pszLAname++ = Type;
	strcpy(pszLAname, pszReadset);
	m_NxtszReadsetIdx += ReadsetNameLen + 2;
	m_LAReadsetNameID = m_NumReadsetNames;
	return(m_LAReadsetNameID);
}

uint32_t		// returned Readset identifier, 0 if unable to locate this Readset name
CWIGutils::LocateReadset(char* pszReadset, // return unique identifier associated with this Readset name
	uint8_t ReadsetType)	// 0: founder, 1: progeny, 2: control
{
	uint32_t ReadsetNameIdx;
	char Type;
	char* pszLAReadset;

	Type = '0' + (char)ReadsetType;
	// with any luck the Readset name will be same as the last accessed
	if ((pszLAReadset = LocateReadset(m_LAReadsetNameID)) != nullptr)
	{
		if (pszLAReadset[-1] == Type &&	// prefixed with ReadsetType
			!stricmp(pszReadset, pszLAReadset))
			return(m_LAReadsetNameID);
	}

	// iterate over all known Readsets
	for (ReadsetNameIdx = 0; ReadsetNameIdx < m_NumReadsetNames; ReadsetNameIdx++)
	{
		pszLAReadset = &m_szReadsetNames[m_szReadsetIdx[ReadsetNameIdx]];
		if (*pszLAReadset++ != Type)		// readset must be that of ReadsetType
			continue;
		if (!stricmp(pszReadset, pszLAReadset))
		{
			m_LAReadsetNameID = ReadsetNameIdx + 1;
			return(m_LAReadsetNameID);
		}
	}
	return(0);
}

char*
CWIGutils::LocateReadset(uint32_t ReadsetID)
{
	uint32_t Idx;
	Idx = ReadsetID & 0x0fffffff;			// mask out any potential ReadsetType
	if (Idx < 1 || Idx > m_NumReadsetNames)
		return(nullptr);
	return(&(m_szReadsetNames[m_szReadsetIdx[Idx - 1] + 1])); // skipping lead char which is the ReadsetType
}




uint32_t		// returned chrom identifier, 0 if unable to accept this chromosome name
CWIGutils::AddChrom(char* pszChrom) // associate unique identifier with this chromosome name
{
int32_t ChromNameIdx;
int ChromNameLen;
char* pszLAname;

// with any luck the sequence name will be same as the last accessed
if ((pszLAname = LocateChrom(m_LAChromNameID)) != nullptr)
	if (!stricmp(pszChrom, pszLAname))
		return(m_LAChromNameID);

// iterate over all known chroms in case this chrom to add is a duplicate
for (ChromNameIdx = 0; ChromNameIdx < m_NumChromNames; ChromNameIdx++)
	if (!stricmp(pszChrom, &m_szChromNames[m_szChromIdx[ChromNameIdx]]))
	{
	m_LAChromNameID = ChromNameIdx + 1;
	return(m_LAChromNameID);
	}

// chrom is not a duplicate
ChromNameLen = (int)strlen(pszChrom);
if ((m_NxtszChromIdx + ChromNameLen + 1) > (int)sizeof(m_szChromNames))
	return(0);
if (m_NumChromNames == cMaxChromNames)
	return(0);

m_szChromIdx[m_NumChromNames++] = m_NxtszChromIdx;
strcpy(&m_szChromNames[m_NxtszChromIdx], pszChrom);
m_NxtszChromIdx += ChromNameLen + 1;
m_LAChromNameID = m_NumChromNames;
return(m_LAChromNameID);
}


uint32_t		// returned chrom identifier, 0 if unable to locate this chromosome name
CWIGutils::LocateChrom(char* pszChrom) // return unique identifier associated with this chromosome name
{
	int32_t ChromNameIdx;
	char* pszLAname;

	// with any luck the sequence name will be same as the last accessed
	if ((pszLAname = LocateChrom(m_LAChromNameID)) != nullptr)
		if (!stricmp(pszChrom, pszLAname))
			return(m_LAChromNameID);

	// iterate over all known chroms
	for (ChromNameIdx = 0; ChromNameIdx < m_NumChromNames; ChromNameIdx++)
		if (!stricmp(pszChrom, &m_szChromNames[m_szChromIdx[ChromNameIdx]]))
		{
			m_LAChromNameID = ChromNameIdx + 1;
			return(m_LAChromNameID);
		}
	return(0);
}

char*
CWIGutils::LocateChrom(uint32_t ChromID)
{
	if (ChromID < 1 || (int32_t)ChromID > m_NumChromNames)
		return(nullptr);
	return(&m_szChromNames[m_szChromIdx[ChromID - 1]]);
}

bool					// true if chrom is accepted, false if chrom not accepted
CWIGutils::AcceptThisChromName(char* pszChrom,   // chromosome name
									 bool bKnown)	// if true then chromosome must have been previously processed and accepted by LoadChromSizes() processing
{
bool bMatch;
if (bKnown)
	return(LocateChrom(pszChrom) < 1 ? false : true);
AcquireSerialise();
bMatch = !m_RegExprs.MatchExcludeRegExpr(pszChrom);
if(bMatch)
    bMatch = m_RegExprs.MatchIncludeRegExpr(pszChrom);
ReleaseSerialise();
return(bMatch);
}

int
CWIGutils::AllocChromMetadata(void)
{
uint32_t ToAllocdChromMetadata;
tsWUChromMetadata* pChromMetadata;
size_t memreq;
if (m_pChromMetadata == nullptr)					// may be nullptr first time in
{
	memreq = cAllocChromMetadata * sizeof(tsWUChromMetadata);
#ifdef _WIN32
	m_pChromMetadata = (tsWUChromMetadata*)malloc((size_t)memreq);
	if (m_pChromMetadata == nullptr)
	{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AllocChromMetadata: Memory allocation of %zd bytes failed", (int64_t)memreq);
		return(eBSFerrMem);
	}
#else
	m_pChromMetadata = (tsWUChromMetadata*)mmap(nullptr, (size_t)memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	if (m_pChromMetadata == MAP_FAILED)
	{
		m_pChromMetadata = nullptr;
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AllocChromMetadata: Memory allocation of %zd bytes through mmap()  failed", (int64_t)memreq, strerror(errno));
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
		size_t memreq = ToAllocdChromMetadata * sizeof(tsWUChromMetadata);
#ifdef _WIN32
		pChromMetadata = (tsWUChromMetadata*)realloc(m_pChromMetadata, memreq);
		if (pChromMetadata == nullptr)
		{
#else
		pChromMetadata = (tsWUChromMetadata*)mremap(m_pChromMetadata, m_AllocdChromMetadataMem, memreq, MREMAP_MAYMOVE);
		if (pChromMetadata == MAP_FAILED)
		{
#endif
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "AllocChromMetadata: Memory reallocation to %zd bytes failed - %s", (int64_t)memreq, strerror(errno));
			return(eBSFerrMem);
		}
		m_pChromMetadata = pChromMetadata;
		m_AllocdChromMetadataMem = memreq;
		m_AllocdChromMetadata = ToAllocdChromMetadata;
		}
return(eBSFSuccess);
	}


tsWUChromMetadata*								// returned pointer to chromosome metadata
CWIGutils::LocateChromMetadataFor(uint32_t ReadSetID,		// readset identifier 
		uint32_t ChromID)			// chrom identifier
{
tsWUReadsetMetadata* pReadsetMetadata;
tsWUChromMetadata* pChromMetadata;
uint32_t CurChromMetadataIdx;

if (ReadSetID > m_NumReadsetNames || ReadSetID == 0)
	return(nullptr);
pReadsetMetadata = &m_Readsets[ReadSetID - 1];
CurChromMetadataIdx = pReadsetMetadata->StartChromMetadataIdx;
for (uint32_t ChromIdx = 0; ChromIdx < pReadsetMetadata->NumChroms && CurChromMetadataIdx != 0; ChromIdx++)
	{
	pChromMetadata = &m_pChromMetadata[CurChromMetadataIdx - 1];
	if (pChromMetadata->ChromID == ChromID)
		return(pChromMetadata);
	CurChromMetadataIdx = pChromMetadata->NxtChromMetadataIdx;
	}
return(nullptr);
}

// Note: an assumption is that within the WIG input file ordering is by chromosome ascending
uint32_t* // loaded PBAs for requested chromosome or nullptr if unable to load
CWIGutils::LoadChromCoverage(uint32_t ReadsetID, // loading chromosome coverage for this readset and mapping coverage as though PBAs
	uint32_t ChromID)    // coverage is for this chromosome
{
int Rslt;
char* pszInWIG;
char szLineBuff[1000];
tsWUReadsetMetadata* pReadsetMetadata;
tsWUChromMetadata* pChromMetadata;
char* pszReadset;
FILE* pInStream;
uint32_t LineNumb;
uint32_t WIGtype;    // 0: unknown, 1: fixed steps, 2 variable steps
uint32_t CurChromID;
uint32_t Coverage;
uint32_t StartLoci;
uint32_t Span;
uint32_t* pCnt;
uint32_t StartLociOfs;
char* pszChrom;
char szChrom[100];
char* pszInBuff;
char Chr;
uint32_t NumElsParsed;

pReadsetMetadata = &m_Readsets[ReadsetID - 1];
pszInWIG = pReadsetMetadata->szFileName;

if ((pszReadset = LocateReadset(ReadsetID)) == nullptr)
	return(nullptr);

if ((pszChrom = LocateChrom(ChromID)) == nullptr)
	return(nullptr);

pChromMetadata = LocateChromMetadataFor(ReadsetID, ChromID);
if (pChromMetadata->pCnts == nullptr)
	{
	if ((pChromMetadata->pCnts = AllocCnts(pChromMetadata->ChromLen)) == nullptr)
		return(nullptr);
	}
memset(pChromMetadata->pCnts, 0, pChromMetadata->ChromLen * sizeof(uint32_t));

// WIG file can now be opened and coverage for requested chromosome parsed and mapped as though PBAs
if ((pInStream = fopen(pszInWIG, "r")) == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadChromCoverage: Unable to open WIG file %s for reading, error: %s", pszInWIG, strerror(errno));
	return(nullptr);
	}
#if _WIN32
_fseeki64(pInStream, pChromMetadata->FileOfs, SEEK_SET);
#else
fseeko64(pInStream, pChromMetadata->FileOfs, SEEK_SET);
#endif
	// load the WIG and where there is coverage then set the corresponding loci in pCnts
	// Note: an assumption is that within the WIG input file ordering is by chromosome ascending
LineNumb = 0;
WIGtype = 0;
CurChromID = 0;
Rslt = eBSFSuccess;

while (fgets(szLineBuff, sizeof(szLineBuff) - 1, pInStream) != nullptr)
	{
	LineNumb++;
	pszInBuff = szLineBuff;
	while((Chr = *pszInBuff) && isspace(Chr))
		pszInBuff++;
	if(Chr == '\0' || Chr == '\n')
		continue;

	if (pszInBuff[0] == 'v' && !strncmp(pszInBuff, "variableStep ", 13))
		{
		NumElsParsed = sscanf(pszInBuff, "variableStep chrom=%s span=%d", szChrom, &Span);
		if (NumElsParsed != 2)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadChromCoverage: Errors parsing WIG line %d - \"%s\" - in WIG file '%s'", LineNumb, pszInBuff, pszInWIG);
			Rslt = -1;
			break;
			}

		// check for valid chrom, a missing chrom relative to PBA could be that the user has specified a limit on number of chromosomes to process so can't treat as a fatal error
		if ((CurChromID = LocateChrom(szChrom)) < 1)
			{
			CurChromID = 0;
			continue;
			}

		if (CurChromID < ChromID) // not yet onto targeted chrom coverage?
			{
			CurChromID = 0;
			continue;
			}

		if (CurChromID > ChromID) // onto coverage for next chromosome?
			break;

		if (Span < 1 || Span > pChromMetadata->ChromLen)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadChromCoverage: Span must be between 1..%d \"%s\" at line %d - \"%s\" - for chromosome '%s' in WIG file '%s'", pChromMetadata->ChromLen, LineNumb, pszInBuff, pszChrom, pszInWIG);
			Rslt = -1;
			break;
			}
		WIGtype = 2;
		continue;
		}
	else
		if (pszInBuff[0] == 'f' && !strncmp(pszInBuff, "fixedStep ", 11))
			{
			NumElsParsed = sscanf(pszInBuff, "fixedStep chrom=%s start=%d step=%d", szChrom, &StartLoci, &Span);
			if (NumElsParsed != 3)
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadChromCoverage: Errors parsing line %d - \"%s\" - in WIG file '%s'", LineNumb, pszInBuff, pszInWIG);
				Rslt = -1;
				break;
				}

			// check for valid chrom, a missing chrom relative to PBA could be that the user has specified a limit on number of chromosomes to process so can't treat as a fatal error
			if ((CurChromID = LocateChrom(szChrom)) < 1)
				{
				CurChromID = 0;
				continue;
				}

			if (CurChromID < ChromID) // not yet onto targeted chrom coverage?
				{
				CurChromID = 0;
				continue;
				}

			if (CurChromID > ChromID) // onto coverage for next chromosome?
				break;

			if (StartLoci < 1 || (StartLoci + Span - 1) >  pChromMetadata->ChromLen)
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadChromCoverage: StartLoci must be between 1..%d, Span %d extending past end of chromosome \"%s\" at line %d - \"%s\" - in WIG file '%s'", pChromMetadata->ChromLen, Span, pszChrom, LineNumb, pszInBuff, pszInWIG);
				Rslt = -1;
				break;
				}
			WIGtype = 1;
			continue;
			}
	if (CurChromID == 0) // 0 if yet to parse out coverage for targeted chromosome
		continue;

	switch (WIGtype) {
		case 0: // unrecognised format, assume some kind of comment or directive, slough
			CurChromID = 0;
			continue;

		case 1: // fixed step
			NumElsParsed = sscanf(pszInBuff, "%d", &Coverage);
			break;

		case 2: // variable step
			NumElsParsed = sscanf(pszInBuff, "%d %d", &StartLoci, &Coverage);
			break;
		}
	if (StartLoci < 1 || StartLoci > pChromMetadata->ChromLen) // ensure not writing past end of chromosome
		continue;
	pCnt = &pChromMetadata->pCnts[StartLoci - 1];
	for (StartLociOfs = 0; StartLociOfs < Span; StartLociOfs++, StartLoci++, pCnt++)
		{
		if (StartLoci > pChromMetadata->ChromLen) // ensure not writing past end of chromosome
			break;
		*pCnt = Coverage;
		}
}
fclose(pInStream);
if (Rslt < 0)
	{
	DeleteSampleChromCnts(ReadsetID, CurChromID);
	pChromMetadata->pCnts = nullptr;
	}
return(pChromMetadata->pCnts);
}

// loading BED which specifies chrom names and sizes
int		// returning number of chromosomes parsed from BED file and accepted after filtering for wildcards
CWIGutils::LoadChromSizes(char* pszBEDFile) // BED file containing chromosome names and sizes
{
int Rslt;
uint32_t ChromID;
int CurFeatureID;
char szFeatName[cMaxDatasetSpeciesChrom];
char szChromName[cMaxDatasetSpeciesChrom];
int32_t StartLoci;
int32_t EndLoci;
int NumChroms;

if ((m_pBedFile = new CBEDfile) == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to instantiate CBEDfile");
	return(eBSFerrObj);
	}

if ((Rslt = m_pBedFile->Open(pszBEDFile, eBTAnyBed)) != eBSFSuccess)
	{
	while (m_pBedFile->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal, gszProcName, m_pBedFile->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadChromSizes: Unable to process '%s' as BED", pszBEDFile);
	return(eBSFerrOpnFile);
	}

NumChroms = 0;
CurFeatureID = 0;
while (Rslt == eBSFSuccess && (CurFeatureID = m_pBedFile->GetNextFeatureID(CurFeatureID)) > 0)
	{
	m_pBedFile->GetFeature(CurFeatureID,	// feature instance identifier
		szFeatName,				// where to return feature name
		szChromName,			// where to return chromosome name
		&StartLoci,				// where to return feature start on chromosome (0..n) 
		&EndLoci);				// where to return feature end on chromosome

	// is this chromosome to be be filtered out?
	if (!AcceptThisChromName(szChromName,false))
		{
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "LoadChromSizes: filtered out chromosome '%s'", szChromName);
		continue;
		}
 
	if (StartLoci != 0)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadChromSizes: unable to accept chromosome '%s' because start loci not 0", szChromName);
		return(eBSFerrChrom);
		}

	if (EndLoci < 0 || (size_t)EndLoci >= cAllocPackedBaseAlleles)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadChromSizes: unable to accept chromosome '%s' because end loci not in acceptable range, must be in range 0..%zd", szChromName, cAllocPackedBaseAlleles);
		return(eBSFerrChrom);
		}

	// expecting a single instance of each chromosome 
	if (LocateChrom(szChromName) != 0)	// must be a unique chromosome name
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadChromSizes: multiple instances of chromosome '%s', chromosome names must be unique", szChromName);
		return(eBSFerrChrom);
		}

	ChromID = AddChrom(szChromName);		// new chromosome
	m_ChromSizes[ChromID-1] = EndLoci+1;	// record it's size
	NumChroms += 1;
	}
m_NumChromSizes = NumChroms;
gDiagnostics.DiagOut(eDLInfo, gszProcName, "LoadChromSizes: accepted %d chromosome sizes from '%s'", NumChroms, pszBEDFile);
delete m_pBedFile;
m_pBedFile = nullptr;
return(NumChroms);
}

int32_t				// returned readset identifier (1..n) or < 0 if errors
CWIGutils::InitialiseMetadata(char *pszReadset,	  // readset name
								char* pszInWIGFile)   // initialise chromosome metadata from file containing WIG coverage
{
int32_t ReadsetID;
FILE* pInStream;
char* pszInWIG;
tsWUReadsetMetadata* pReadsetMetadata;
tsWUChromMetadata* pChromMetadata;
tsWUChromMetadata* pPrevChromMetadata;

if ((ReadsetID = AddReadset(pszReadset, 0)) == 0)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "InitialiseMetadata: Input file '%s' duplicates the ReadsetID '%s' of a previously loaded readset", pszInWIGFile, pszReadset);
	return(eBSFerrOpnFile);
	}

pszInWIG = pszInWIGFile;
if ((pInStream = fopen(pszInWIG, "r")) == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "InitialiseMetadata: Unable to open WIG file %s for reading, error: %s", pszInWIG, strerror(errno));
	Reset();
	return(eBSFerrOpnFile);
	}

pReadsetMetadata = &m_Readsets[ReadsetID - 1];
memset(pReadsetMetadata, 0, sizeof(*pReadsetMetadata));
pReadsetMetadata->NumChroms = 0;
strcpy(pReadsetMetadata->szFileName, pszInWIGFile);
pReadsetMetadata->ReadsetID = ReadsetID;
pReadsetMetadata->StartChromID = 0;
pReadsetMetadata->StartChromMetadataIdx = 0;

int32_t PrevChromMetadataIdx = 0;
int32_t Span;
int32_t StartLoci;
uint32_t LineNumb = 0;
uint32_t CurChromID = 0;
int32_t Rslt = eBSFSuccess;
char Chr;
char szLineBuff[1000];
char szChrom[100];
char* pszInBuff;
uint32_t NumElsParsed;
uint32_t ChromID = 0;
size_t CurFileOfs = 0;
while (fgets(szLineBuff, sizeof(szLineBuff) - 1, pInStream) != nullptr)
	{
	LineNumb++;
	pszInBuff = szLineBuff;
	while((Chr = *pszInBuff) && isspace(Chr))
		pszInBuff++;
	if(Chr == '\0' || Chr == '\n' || !(Chr == 'v' || Chr == 'f'))
		continue;

	if (pszInBuff[0] == 'v' && !strncmp(pszInBuff, "variableStep ", 13))
		{
		NumElsParsed = sscanf(pszInBuff, "variableStep chrom=%s span=%d", szChrom, &Span);
		if (NumElsParsed != 2)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadChromCoverage: Errors parsing WIG line %d - \"%s\" - in WIG file '%s'", LineNumb, pszInBuff, pszInWIG);
			Rslt = -1;
			break;
			}
		}
	else
		if (pszInBuff[0] == 'f' && !strncmp(pszInBuff, "fixedStep ", 11))
			{
			NumElsParsed = sscanf(pszInBuff, "fixedStep chrom=%s start=%d step=%d", szChrom, &StartLoci, &Span);
			if (NumElsParsed != 3)
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadChromCoverage: Errors parsing line %d - \"%s\" - in WIG file '%s'", LineNumb, pszInBuff, pszInWIG);
				Rslt = -1;
				break;
				}
			}
		else
			continue;

	// check for chrom being of interest 
	if ((CurChromID = LocateChrom(szChrom)) < 1)
		{
		CurChromID = 0;		// of no interest
		continue;
		}

	if (CurChromID == ChromID) // same chromosome as previously encountered
		continue;

	ChromID = CurChromID;
	if ((Rslt = AllocChromMetadata()) != eBSFSuccess)
		return(Rslt);

	pChromMetadata = &m_pChromMetadata[m_UsedNumChromMetadata++];
	if (PrevChromMetadataIdx != 0)
		{
		pPrevChromMetadata = &m_pChromMetadata[PrevChromMetadataIdx - 1];
		pPrevChromMetadata->NxtChromMetadataIdx = m_UsedNumChromMetadata;
		}
	else
		{
		pReadsetMetadata->StartChromID = ChromID;
		pReadsetMetadata->StartChromMetadataIdx = m_UsedNumChromMetadata;
		}
#if _WIN32
	CurFileOfs = _ftelli64(pInStream);
#else
	CurFileOfs = ftello64(pInStream);
#endif

	CurFileOfs -= 3*strlen(szLineBuff); // allows for source file to have been in 16-bit unicode character format
	if((int64_t)CurFileOfs < 0)
		CurFileOfs = 0;
	pReadsetMetadata->NumChroms++;
	PrevChromMetadataIdx = m_UsedNumChromMetadata;
	pChromMetadata->ChromID = ChromID;
	pChromMetadata->ChromMetadataIdx = m_UsedNumChromMetadata;
	pChromMetadata->NxtChromMetadataIdx = 0;
	pChromMetadata->ChromLen = m_ChromSizes[CurChromID-1];
	pChromMetadata->ReadsetID = ReadsetID;
	pChromMetadata->FileOfs = CurFileOfs;
	pChromMetadata->pCnts = nullptr;
	}
fclose(pInStream);

return(ReadsetID);
}


int 
CWIGutils::Process(eWIGuMode PMode,	// processing mode: eWIGu2DEseq WIG to DEseq, eWIGu2TODEseq WIG to turnover DEseq
			int32_t NumInputFiles,		// number of input WIG file specs
			char* pszInputFiles[],		// names of input WIG files (wildcards allowed)
			char* pszOutFile,			// output to this file
			char* pszChromFile,			// BED file containing reference assembly chromosome names and sizes
			char* pszROIFile,			// BED file containing regions/genes of interest
			int NumIncludeChroms,		// number of chromosome regular expressions to include
			char *pszIncludeChroms[],	// array of include chromosome regular expressions
			int NumExcludeChroms,		// number of chromosome expressions to exclude
			char *pszExcludeChroms[],	// array of exclude chromosome regular expressions
			int32_t NumThreads)			// number of worker threads to use
{
int Rslt;
Reset();
CreateMutexes();
m_PMode = PMode;
strcpy(m_szOutFile, pszOutFile);

// compile include/exclude chromosome regexpr if user has specified alignments to be filtered by chrom
if(Rslt = (m_RegExprs.CompileREs(NumIncludeChroms, pszIncludeChroms,NumExcludeChroms, pszExcludeChroms)) < eBSFSuccess)
	{
	Reset();
	return(Rslt);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing chromosome sizes from file '%s'", pszChromFile);
if ((Rslt = LoadChromSizes(pszChromFile)) < 1) // BED file containing chromosome names and sizes - NOTE chromosomes will be filtered by include/exclude wildcards
	{
	Reset();
	return(Rslt);
	}

// load regions of interest (genes/exons) from file
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading regions of interest file '%s'", pszROIFile);
if((m_pROIFile = new CBEDfile()) == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBEDfile");
	Reset();
	return(eBSFerrObj);
	}
if((Rslt=m_pROIFile->Open(pszROIFile))!=eBSFSuccess)
	{
	while(m_pROIFile->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pROIFile->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open regions of interest file '%s'",pszROIFile);
	Reset();
	return(Rslt);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Regions of interest file '%s' loaded", pszROIFile);


m_NumThreads = NumThreads;

	// initial allocation, will be realloc'd if more memory required
size_t memreq = (size_t)cAllocChromMetadata * sizeof(tsWUChromMetadata);
#ifdef _WIN32
m_pChromMetadata = (tsWUChromMetadata*)malloc(memreq);	// initial and perhaps the only allocation
if (m_pChromMetadata == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Initial memory allocation of %zd bytes for chromosome metadata failed - %s", (int64_t)memreq, strerror(errno));
	Reset();
	return(eBSFerrMem);
	}
#else
m_pChromMetadata = (tsWUChromMetadata*)mmap(nullptr, memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
if (m_pChromMetadata == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Memory allocation of %zd bytes through mmap() for chromosome metadata failed - %s", (int64_t)memreq, strerror(errno));
	m_pChromMetadata = nullptr;
	Reset();
	return(eBSFerrMem);
	}
#endif
m_AllocdChromMetadataMem = memreq;
m_UsedNumChromMetadata = 0;
m_AllocdChromMetadata = cAllocChromMetadata;

CSimpleGlob glob(SG_GLOB_FULLSORT);
int32_t ReadsetID;
int32_t NumFiles;
int32_t TotNumFiles;
int32_t Idx;
char* pszInFile;
gDiagnostics.DiagOut(eDLInfo, gszProcName, "Process: Starting to load input PBA files");
Rslt = eBSFSuccess;		// assume success!
TotNumFiles = 0;
for (Idx = 0; Idx < NumInputFiles; Idx++)
	{
	glob.Init();
	if (glob.Add(pszInputFiles[Idx]) < SG_SUCCESS)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to glob '%s' input file spec", pszInputFiles[Idx]);
		Reset();
		return(eBSFerrOpnFile);	// treat as though unable to open file
		}
	if ((NumFiles = glob.FileCount()) <= 0)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to locate any input file matching '%s'", pszInputFiles[Idx]);
		Reset();
		return(eBSFerrFileName);
		}

	TotNumFiles += NumFiles;

	Rslt = eBSFSuccess;
	for (int FileID = 0; Rslt >= eBSFSuccess && FileID < NumFiles; ++FileID)
		{
		pszInFile = glob.File(FileID);
		char szReadsetFileName[_MAX_FNAME];
		char szReadsetID[_MAX_FNAME];
		CUtility::splitpath(pszInFile,nullptr,szReadsetFileName);
		// Ugh, Ugh, and Ugh again
		// An assumption is that all WIG file names consist of the readset identifier suffixed by '_vs_' followed by the targeted assembly name
		// Thus the readset identifer is parsed from the file name using the above assumption!
		char Chr;
		char *pFN;
		char *pSID;
		pFN = szReadsetFileName;
		pSID = szReadsetID;
		*pSID = '\0';
		while (Chr = *pFN++)
			{
			if(*pFN == '\0'|| pFN[1] == '\0' || pFN[2] == '\0')
				break;
			if (Chr == '_')
				{
				if(*pFN == 'v' && pFN[1] == 's' && pFN[2] == '_')
					break;
				}
			*pSID++ = Chr;
			*pSID = '\0';
			}
		if(szReadsetID[0] == '\0')
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Error parsing readset identifier from file name '%s'", szReadsetFileName);
			Reset();
			return(ReadsetID);
			}
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Process: Initialising metadata for input file '%s', readset identifier: '%s'", pszInFile,szReadsetID);
		ReadsetID = InitialiseMetadata(szReadsetID, pszInFile);		 // will initialise readset + chrom metadata
		if (ReadsetID <= 0)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Errors loading file '%s'", pszInFile);
			Reset();
			return(ReadsetID);
			}
		}
	}
m_NumReadsetIDs = ReadsetID;
gDiagnostics.DiagOut(eDLInfo, gszProcName, "Process: Completed initialising metadata for input files");

// try loading coverage counts for all readsets
char *pszLineBuff;
int LineOfs;
tsChromROI *pRSChromROI;
tsChromROI *pCurRS;
uint32_t CurChromMetadataIdx;
tsWUChromMetadata* pChromMetadata;
tsWUReadsetMetadata* pReadsetMetadata;
uint32_t Cnts;
uint32_t ChromID;
uint32_t ChromIdx;

char szFileFeatExonTurnOver[_MAX_FNAME];
FILE *pOutFeatExonTurnOver;
char szFileFeatExonCnts[_MAX_FNAME];
FILE *pOutFeatExonCnts;
char szFileFeatTurnOver[_MAX_FNAME];
FILE *pOutFeatTurnOver;
char szFileFeatCnts[_MAX_FNAME];
FILE *pOutFeatCnts;

CUtility::AppendFileNameSuffix(szFileFeatExonTurnOver,m_szOutFile,(char *)".exonturnover.csv",'.');
CUtility::AppendFileNameSuffix(szFileFeatTurnOver,m_szOutFile,(char *)".featturnover.csv",'.');
CUtility::AppendFileNameSuffix(szFileFeatExonCnts,m_szOutFile,(char *)".exoncnts.csv",'.');
CUtility::AppendFileNameSuffix(szFileFeatCnts,m_szOutFile,(char *)".featcnts.csv",'.');

pszLineBuff = new char [500000];
pRSChromROI = new tsChromROI[m_NumReadsetIDs];

if((pOutFeatExonTurnOver = fopen(szFileFeatExonTurnOver,"w"))==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create/truncate file %s for writing error: %s",szFileFeatExonTurnOver,strerror(errno));
	Reset();
	return(eBSFerrOpnFile);
	}
if((pOutFeatTurnOver = fopen(szFileFeatTurnOver,"w"))==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create/truncate file %s for writing error: %s",szFileFeatTurnOver,strerror(errno));
	fclose(pOutFeatExonTurnOver);
	Reset();
	return(eBSFerrOpnFile);
	}
if((pOutFeatCnts = fopen(szFileFeatCnts,"w"))==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create/truncate file %s for writing error: %s",szFileFeatCnts,strerror(errno));
	fclose(pOutFeatExonTurnOver);
	fclose(pOutFeatTurnOver);
	Reset();
	return(eBSFerrOpnFile);
	}
if((pOutFeatExonCnts = fopen(szFileFeatExonCnts,"w"))==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create/truncate file %s for writing error: %s",szFileFeatCnts,strerror(errno));
	fclose(pOutFeatExonTurnOver);
	fclose(pOutFeatTurnOver);
	fclose(pOutFeatCnts);
	Reset();
	return(eBSFerrOpnFile);
	}
LineOfs = sprintf(pszLineBuff,"\"Feature\",\"Chrom\",\"Start\",\"End\",\"Strand\",\"NumExons\",\"FeatLength\",\"TransLength\"");
for (ReadsetID = 1; ReadsetID <= m_NumReadsetIDs; ReadsetID++)
	LineOfs += sprintf(&pszLineBuff[LineOfs],",\"%s\"",LocateReadset(ReadsetID));
LineOfs+=sprintf(&pszLineBuff[LineOfs],"\n");
fwrite(pszLineBuff,1,LineOfs,pOutFeatExonTurnOver);
fwrite(pszLineBuff,1,LineOfs,pOutFeatTurnOver);
fwrite(pszLineBuff,1,LineOfs,pOutFeatCnts);
fwrite(pszLineBuff,1,LineOfs,pOutFeatExonCnts);
LineOfs = 0;
pReadsetMetadata = &m_Readsets[0];
CurChromMetadataIdx = pReadsetMetadata->StartChromMetadataIdx;
for (ChromIdx = 0; ChromIdx < pReadsetMetadata->NumChroms && CurChromMetadataIdx != 0; ChromIdx++)
	{
	pChromMetadata = &m_pChromMetadata[CurChromMetadataIdx - 1];
	ChromID = pChromMetadata->ChromID;
	int ROIChromID;
	int ROILoci;
	int ROIFeatID;
	int ROIFeatStart;
	int ROIFeatEnd;
	int ROINumExons;
	char ROIFeatStrand;
	bool bInExon;
	int TranscribedLength;
	int FeatLength;
	bool bTot50Cnts;
	bool bTot50ExonCnts;
	int CurTranscriptOfs;
	int CurFeatOfs;
	char* pszChrom;
	char szROIFeatName[200];
	char szROIFeatChrom[200];
	pszChrom = LocateChrom(pChromMetadata->ChromID);
	ROIChromID = m_pROIFile->LocateChromIDbyName(pszChrom);
	if (ROIChromID < 0)
		continue;
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Process: Loading coverage for: '%s'", LocateChrom(ChromID));
	pCurRS = pRSChromROI;
	for (ReadsetID = 0; ReadsetID < m_NumReadsetIDs; ReadsetID++,pCurRS++)
		{
		memset(pCurRS,0,sizeof(tsChromROI));
		pCurRS->ReadsetID = ReadsetID;
		if ((pCurRS->pCnts = LoadChromCoverage(ReadsetID+1, pChromMetadata->ChromID)) == nullptr)
			break;
		}
	if(ReadsetID < m_NumReadsetIDs)	// all readsets must have read coverage on the current chromosome 
		continue;

		// iterate over features on this chromosome
	ROIFeatID = 0;
	while ((ROIFeatID = m_pROIFile->GetNextChromFeatureID(ROIChromID, ROIFeatID)) > 0)
		{
		pCurRS = pRSChromROI;
		for (ReadsetID = 0; ReadsetID < m_NumReadsetIDs; ReadsetID++,pCurRS++)
			memset(&pCurRS->ROICnts,0,sizeof(pCurRS->ROICnts));
		m_pROIFile->GetFeature(ROIFeatID,szROIFeatName,szROIFeatChrom,&ROIFeatStart,&ROIFeatEnd,nullptr,&ROIFeatStrand);
		TranscribedLength = m_pROIFile->GetTranscribedLen(ROIFeatID); // sum of all exon lengths
		FeatLength = m_pROIFile->GetFeatLen(ROIFeatID) + 1;
		ROINumExons = m_pROIFile->GetNumExons(ROIFeatID);

		LineOfs = sprintf(pszLineBuff,"\"%s\",\"%s\",%d,%d,\"%c\",%d,%d,%d",szROIFeatName,szROIFeatChrom,ROIFeatStart,ROIFeatEnd,ROIFeatStrand,ROINumExons,FeatLength,TranscribedLength);
		int MarkVarOfs = LineOfs;
		CurTranscriptOfs = 0;
		CurFeatOfs = 0;
		for (ROILoci = ROIFeatStart; ROILoci <= ROIFeatEnd; ROILoci++)
			{
			CurFeatOfs += 1;
			bTot50Cnts = CurFeatOfs <= FeatLength/2 ? true : false; 
			if(bInExon = m_pROIFile->InExon(ROIFeatID,ROILoci,ROILoci))
				CurTranscriptOfs += 1;
			bTot50ExonCnts = CurTranscriptOfs <= TranscribedLength/2 ? true:false;
			pCurRS = pRSChromROI;
			for (ReadsetID = 0; ReadsetID < m_NumReadsetIDs; ReadsetID++,pCurRS++)
				{
				Cnts = pCurRS->pCnts[ROILoci];
				if(Cnts == 0)
					continue;
				pCurRS->ROICnts.NumLociCoverage++;
				pCurRS->ROICnts.TotCnts += Cnts;
				if(bTot50Cnts)
					pCurRS->ROICnts.Tot50Cnts += Cnts;
				if(bInExon)
					{
					pCurRS->ROICnts.NumExonLociCoverage++;
					pCurRS->ROICnts.TotExonCnts += Cnts;
					if(bTot50ExonCnts)
						pCurRS->ROICnts.Tot50ExonCnts += Cnts;	
					}
				}
			}
		pCurRS = pRSChromROI;
		for (ReadsetID = 0; ReadsetID < m_NumReadsetIDs; ReadsetID++,pCurRS++)
			{
			if(pCurRS->ROICnts.TotExonCnts >= 5000)	// have to have a threshold else may be too noisy!
				pCurRS->ROICnts.ExonTurnOver = max(0.05,(double)pCurRS->ROICnts.Tot50ExonCnts / pCurRS->ROICnts.TotExonCnts);
			else
				pCurRS->ROICnts.ExonTurnOver = 0.0;

			if(pCurRS->ROICnts.TotCnts >= 5000)	// have to have a threshold else may be too noisy!
				pCurRS->ROICnts.TurnOver = max(0.05,(double)pCurRS->ROICnts.Tot50Cnts / pCurRS->ROICnts.TotCnts);
			else
				pCurRS->ROICnts.TurnOver = 0.0;
			}
		pCurRS = pRSChromROI;
		LineOfs = MarkVarOfs;
		for (ReadsetID = 0; ReadsetID < m_NumReadsetIDs; ReadsetID++,pCurRS++)
			LineOfs+=sprintf(&pszLineBuff[LineOfs],",%.3f",pCurRS->ROICnts.ExonTurnOver);
		LineOfs+=sprintf(&pszLineBuff[LineOfs],"\n");
		fwrite(pszLineBuff,1,LineOfs,pOutFeatExonTurnOver);

		pCurRS = pRSChromROI;
		LineOfs = MarkVarOfs;
		for (ReadsetID = 0; ReadsetID < m_NumReadsetIDs; ReadsetID++,pCurRS++)
			LineOfs+=sprintf(&pszLineBuff[LineOfs],",%.3f",pCurRS->ROICnts.TurnOver);
		LineOfs+=sprintf(&pszLineBuff[LineOfs],"\n");
		fwrite(pszLineBuff,1,LineOfs,pOutFeatTurnOver);

		pCurRS = pRSChromROI;
		LineOfs = MarkVarOfs;
		for (ReadsetID = 0; ReadsetID < m_NumReadsetIDs; ReadsetID++,pCurRS++)
			LineOfs+=sprintf(&pszLineBuff[LineOfs],",%zd",pCurRS->ROICnts.TotCnts);
		LineOfs+=sprintf(&pszLineBuff[LineOfs],"\n");
		fwrite(pszLineBuff,1,LineOfs,pOutFeatCnts);

		pCurRS = pRSChromROI;
		LineOfs = MarkVarOfs;
		for (ReadsetID = 0; ReadsetID < m_NumReadsetIDs; ReadsetID++,pCurRS++)
			LineOfs+=sprintf(&pszLineBuff[LineOfs],",%zd",pCurRS->ROICnts.TotExonCnts);
		LineOfs+=sprintf(&pszLineBuff[LineOfs],"\n");
		fwrite(pszLineBuff,1,LineOfs,pOutFeatExonCnts);
		LineOfs = 0;
		}
	CurChromMetadataIdx = pChromMetadata->NxtChromMetadataIdx;
	DeleteAllChromCnts();	// ready for next chromosome
	}
DeleteAllChromCnts();
fclose(pOutFeatExonTurnOver);
fclose(pOutFeatTurnOver);
fclose(pOutFeatExonCnts);
fclose(pOutFeatCnts);
delete []pRSChromROI;
delete []pszLineBuff;
Reset();
return(eBSFSuccess);
}