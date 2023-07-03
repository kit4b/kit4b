/*
'kit4b' is being released under the Opensource Software License Agreement (GPLv3)
'kit4b' is Copyright (c) 2019
Please contact Dr Stuart Stephen < stuartjs@g3web.com > if you have any questions regarding 'kit4b'.
*/


#include "stdafx.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#if _WIN32
//#define _CRTDBG_MAP_ALLOC
//#include <crtdbg.h>
#include <process.h>
#include "../libkit4b/commhdrs.h"
#else
#include <sys/mman.h>
#include <pthread.h>
#include "../libkit4b/commhdrs.h"
#endif

#include "./ngskit4b.h"
#include "CDGTvQTLs.h"

#include "../libkit4b/bgzf.h"


int Process(eModeDGTA PMode,				// processing mode
			double MinCoverage,				// if coverage < this threshold then class as being low coverage
			double HomozPropThres,			// if proportion of samples in Grp1Prop is >= this proportion then characterise as homozygous
			char* pszAssembRefFile,			// contains original reference assembly, PBA, against which samples were aligned
			char *pszChromFile,				// BED file, contains chromosome names and sizes
			char *pszDGTsFile,				// file containing DGT loci and allele sample groupings
			char *pszQTLsFile,				// file containing QTL SNP loci
			int NumPBAFiles,				// number of input PBA file specs, these are the source PBA files used during the grouping generation
			char **ppszPBAFiles,			// names of input PBA files (wildcards allowed)
			int	NumIncludeChroms,			// number of chromosome regular expressions to include
			char** ppszIncludeChroms,		// array of include chromosome regular expressions
			int	NumExcludeChroms,			// number of chromosome expressions to exclude
			char** ppszExcludeChroms,		// array of exclude chromosome regular expressions
			char* pszRsltsFileBaseName,		// write results to this file base name - will be suffixed with result types
			int NumThreads);				// number of worker threads to use

#ifdef _WIN32
int DGTs(int argc, char* argv[])
{
	// determine my process name
	_splitpath(argv[0], nullptr, nullptr, gszProcName, nullptr);
#else
int
DGTs(int argc, char** argv)
{
	// determine my process name
	CUtility::splitpath((char*)argv[0], nullptr, gszProcName);
#endif
int Idx;
int iFileLogLevel;			// level of file diagnostics
int iScreenLogLevel;		// level of file diagnostics
char szLogFile[_MAX_PATH];	// write diagnostics to this file
int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure
int NumberOfProcessors;		// number of installed CPUs
int NumThreads;				// number of threads (0 defaults to number of CPUs or a maximum of cMaxPBAWorkerThreads)


eModeDGTA PMode;				// processing mode 

double MinCoverage;				// if coverage < this threshold then class as being low coverage
double HomozPropThres;			// if proportion of samples in Grp1Prop is >= this proportion then characterise as homozygous

int NumIncludeChroms;
char* pszIncludeChroms[cMaxIncludeChroms];
int NumExcludeChroms;
char* pszExcludeChroms[cMaxExcludeChroms];

int NumPBAFiles;		// number of input BPA files
char* pszPBAFiles[cMaxPBAFileSpecs];		// names of input BPA files (wildcards allowed)
char szDGTsFile[_MAX_PATH];				// input QTLs file
char szQTLsFile[_MAX_PATH];				// input QTLs file
char szAssembRefFile[_MAX_PATH];		// contains orginal reference assembly, diplotype PBA, against which samples were aligned
char szChromFile[_MAX_PATH];	// BED file containing chromosome names and sizes
char szRsltsFileBaseName[_MAX_PATH];	// results written to this file base name - will be suffixed with results type

struct arg_lit* help = arg_lit0("h", "help", "print this help and exit");
struct arg_lit* version = arg_lit0("v", "version,ver", "print version information and exit");
struct arg_int* FileLogLevel = arg_int0("f", "FileLogLevel", "<int>", "Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file* LogFile = arg_file0("F", "log", "<file>", "diagnostics log file");

struct arg_int* pmode = arg_int0("m", "mode", "<int>", "processing mode 0: QTL only processing, 1: DGT and QTL processing (default 0)");

struct arg_dbl* mincovp = arg_dbl0("k", "mincovp", "<dbl>", "if coverage < this threshold then characterise loci as being low coverage (default 0.8)");
struct arg_dbl* homozp = arg_dbl0("p", "homozp", "<dbl>", "minimum proportion of samples in Grp1Prop to characterise loci as being homogenous (default 0.95)");

struct arg_str* chromexclude = arg_strn("Z", "chromexclude", "<string>", 0, cMaxExcludeChroms, "high priority - regular expressions defining chromosomes to exclude");
struct arg_str* chromeinclude = arg_strn("z", "chromeinclude", "<string>", 0, cMaxIncludeChroms, "low priority - regular expressions defining chromosomes to include");

struct arg_file* chromfile = arg_file1("c", "chromfile", "<file>", "input BED file containing chromosome names and sizes");
struct arg_file* qtlsfile = arg_file1("Q", "qtlsfile", "<file>", "input CSV QTLs alleles file");
struct arg_file* dgtsfile = arg_file1("q", "dgtsfile", "<file>", "input CSV DGTs alleles file");
struct arg_file* samplefiles = arg_filen("i", "samplefiles", "<file>", 1, cMaxPBAFileSpecs, "input sample BPA file(s), wildcards allowed, limit of 500 filespecs supported");
struct arg_file* assembrefile = arg_file1("I", "inassembfasta", "<file>", "reference assembly (diplotype PBA) 2file against which alignments were made");
struct arg_file* outfile = arg_file1("o", "out", "<file>", "generated analysis output to this file base name - willl be suffixed with file type");
struct arg_int* threads = arg_int0("T", "threads", "<int>", "number of processing threads 0..64 (defaults to 0 which limits threads to maximum of 64 CPU cores)");
struct arg_end* end = arg_end(200);

void* argtable[] = { help,version,FileLogLevel,LogFile,
						pmode, mincovp,homozp,assembrefile,chromfile,chromeinclude,chromexclude,samplefiles, qtlsfile,dgtsfile, outfile,threads,end };

char** pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc, (char**)argv, &pAllArgs);
if (argerrors >= 0)
	argerrors = arg_parse(argerrors, pAllArgs, argtable);

/* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
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
if (version->count > 0)
	{
	printf("\n%s %s Version %s\n", gszProcName, gpszSubProcess->pszName, kit4bversion);
	return(1);
	}

if (!argerrors)
	{
	if (FileLogLevel->count && !LogFile->count)
		{
		printf("\nError: FileLogLevel '-f%d' specified but no logfile '-F<logfile>\n'", FileLogLevel->ival[0]);
		exit(1);
		}

	iScreenLogLevel = iFileLogLevel = FileLogLevel->count ? FileLogLevel->ival[0] : eDLInfo;
	if (iFileLogLevel < eDLNone || iFileLogLevel > eDLDebug)
		{
		printf("\nError: FileLogLevel '-l%d' specified outside of range %d..%d\n", iFileLogLevel, eDLNone, eDLDebug);
		exit(1);
		}

	if (LogFile->count)
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
	if (!gDiagnostics.Open(szLogFile, (etDiagLevel)iScreenLogLevel, (etDiagLevel)iFileLogLevel, true))
		{
		printf("\nError: Unable to start diagnostics subsystem\n");
		if (szLogFile[0] != '\0')
			printf(" Most likely cause is that logfile '%s' can't be opened/created\n", szLogFile);
		exit(1);
		}

	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Subprocess %s Version %s starting", gpszSubProcess->pszName, kit4bversion);

	PMode = pmode->count ? (eModeDGTA)pmode->ival[0] : eDGTADefault;
	if (PMode < eDGTADefault || PMode >= eDGTAPlaceHolder)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Unsupported processing mode '-m%d'\n", PMode);
		exit(1);
		}

	MinCoverage = mincovp->count ? mincovp->dval[0] : cDfltMinCoverage; // silently clamp to a reasonable value
	if (MinCoverage < 0.75)
		MinCoverage = 0.75;
	else
		if(MinCoverage > 0.95)
			MinCoverage = 0.95;

	HomozPropThres = homozp->count ? homozp->dval[0] : cDfltHomozPropThres; // silently clamp to a reasonable value
	if (HomozPropThres < MinCoverage)
		HomozPropThres = MinCoverage + 0.025;
	else
		if (HomozPropThres > 0.99)
			HomozPropThres = 0.99;

	if (assembrefile->count)
		{
		strncpy(szAssembRefFile, assembrefile->filename[0], _MAX_PATH);
		szAssembRefFile[_MAX_PATH - 1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(szAssembRefFile);
		}
	else
		szAssembRefFile[0] = 0;
	if (szAssembRefFile[0] == 0)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "No input reference assembly file specified");
		exit(1);
		}

	if (chromfile->count)
		{
		strncpy(szChromFile, chromfile->filename[0], _MAX_PATH);
		szChromFile[_MAX_PATH - 1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(szChromFile);
		}
	else
		szChromFile[0] = 0;
	if (szChromFile[0] == 0)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "No BED file containing chromosome names and sizes specified");
		exit(1);
		}

	if (chromeinclude->count)
		{
		for (NumIncludeChroms = Idx = 0; NumIncludeChroms < cMaxIncludeChroms && Idx < chromeinclude->count; Idx++)
			{
			pszIncludeChroms[Idx] = nullptr;
			if (pszIncludeChroms[NumIncludeChroms] == nullptr)
				pszIncludeChroms[NumIncludeChroms] = new char[_MAX_PATH];
			strncpy(pszIncludeChroms[NumIncludeChroms], chromeinclude->sval[Idx], _MAX_PATH);
			pszIncludeChroms[NumIncludeChroms][_MAX_PATH - 1] = '\0';
			CUtility::TrimQuotedWhitespcExtd(pszIncludeChroms[NumIncludeChroms]);
			if (pszIncludeChroms[NumIncludeChroms][0] != '\0')
				NumIncludeChroms++;
			}
		}
	else
		NumIncludeChroms = 0;


	if (chromexclude->count)
		{
		for (NumExcludeChroms = Idx = 0; NumExcludeChroms < cMaxExcludeChroms && Idx < chromexclude->count; Idx++)
			{
			pszExcludeChroms[Idx] = nullptr;
			if (pszExcludeChroms[NumExcludeChroms] == nullptr)
				pszExcludeChroms[NumExcludeChroms] = new char[_MAX_PATH];
			strncpy(pszExcludeChroms[NumExcludeChroms], chromexclude->sval[Idx], _MAX_PATH);
			pszExcludeChroms[NumExcludeChroms][_MAX_PATH - 1] = '\0';
			CUtility::TrimQuotedWhitespcExtd(pszExcludeChroms[NumExcludeChroms]);
			if (pszExcludeChroms[NumExcludeChroms][0] != '\0')
				NumExcludeChroms++;
			}
		}
	else
		NumExcludeChroms = 0;


	if (dgtsfile->count)
		{
		strncpy(szDGTsFile, dgtsfile->filename[0], _MAX_PATH);
		szDGTsFile[_MAX_PATH - 1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(szDGTsFile);
		}
	else
		szDGTsFile[0] = 0;
	if (szDGTsFile[0] == 0)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "No input DGTs file specified");
		exit(1);
		}


	if (qtlsfile->count)
		{
		strncpy(szQTLsFile, qtlsfile->filename[0], _MAX_PATH);
		szQTLsFile[_MAX_PATH - 1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(szQTLsFile);
		}
	else
		szQTLsFile[0] = 0;
	if (szQTLsFile[0] == 0)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "No input QTLs file specified");
		exit(1);
		}


	NumPBAFiles = 0;
	if (samplefiles->count)
		{
		for (NumPBAFiles = Idx = 0; NumPBAFiles < cMaxPBAFileSpecs && Idx < samplefiles->count; Idx++)
			{
			pszPBAFiles[Idx] = nullptr;
			if (pszPBAFiles[NumPBAFiles] == nullptr)
				pszPBAFiles[NumPBAFiles] = new char[_MAX_PATH];
			strncpy(pszPBAFiles[NumPBAFiles], samplefiles->filename[Idx], _MAX_PATH);
			pszPBAFiles[NumPBAFiles][_MAX_PATH - 1] = '\0';
			CUtility::TrimQuotedWhitespcExtd(pszPBAFiles[NumPBAFiles]);
			if (pszPBAFiles[NumPBAFiles][0] != '\0')
				NumPBAFiles++;
			}
		}

	if (!NumPBAFiles)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: After removal of whitespace, there are no input sample PBA file(s) specified with '-i<filespec>' option)\n");
		exit(1);
		}



	strcpy(szRsltsFileBaseName, outfile->filename[0]);
	CUtility::TrimQuotedWhitespcExtd(szRsltsFileBaseName);
	if (szRsltsFileBaseName[0] == '\0')
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "No output file base name specified");
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

	int MaxAllowedThreads = min(cMaxPBAWorkerThreads, NumberOfProcessors);	// limit to be at most cMaxWorkerThreads
	if ((NumThreads = threads->count ? threads->ival[0] : MaxAllowedThreads) == 0)
		NumThreads = MaxAllowedThreads;
	if (NumThreads < 0 || NumThreads > MaxAllowedThreads)
		{
		gDiagnostics.DiagOut(eDLWarn, gszProcName, "Warning: Number of threads '-T%d' specified was outside of range %d..%d", NumThreads, 1, MaxAllowedThreads);
		gDiagnostics.DiagOut(eDLWarn, gszProcName, "Warning: Defaulting number of threads to %d", MaxAllowedThreads);
		NumThreads = MaxAllowedThreads;
		}


	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Processing parameters:");
	const char* pszDescr;
	switch (PMode) {
	
		case eDGTADGT:
			pszDescr = "DGT and QTL loci processing";
			break;
		default:
			PMode = eDGTADefault;
			pszDescr = "QTL loci only processing";
			break;
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo, "DGT vs. QTL: '%s'", pszDescr);

	gDiagnostics.DiagOutMsgOnly(eDLInfo, "Minimum coverage : '%.4f'", MinCoverage);
	gDiagnostics.DiagOutMsgOnly(eDLInfo, "Minimum Grp1Prop required before characterising as homozygous : '%.4f'", HomozPropThres);
	
	gDiagnostics.DiagOutMsgOnly(eDLInfo, "Reference assembly : '%s'", szAssembRefFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo, "BED containing chromosome names and sizes : '%s'", szChromFile);
	for (Idx = 0; Idx < NumIncludeChroms; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "reg expressions defining chroms to include: '%s'", pszIncludeChroms[Idx]);
	for (Idx = 0; Idx < NumExcludeChroms; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "reg expressions defining chroms to exclude: '%s'", pszExcludeChroms[Idx]);
	gDiagnostics.DiagOutMsgOnly(eDLInfo, "DGTs file : '%s'", szDGTsFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo, "QTLs file : '%s'", szQTLsFile);
	for (Idx = 0; Idx < NumPBAFiles; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "Sample PBA file : '%s'", pszPBAFiles[Idx]);
	gDiagnostics.DiagOutMsgOnly(eDLInfo, "Output file base name : '%s'", szRsltsFileBaseName);

gDiagnostics.DiagOutMsgOnly(eDLInfo, "number of threads : %d", NumThreads);

#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = Process(PMode,		// processing mode
					MinCoverage,				// if coverage < this threshold then class as being low coverage
					HomozPropThres,		// if proportion of samples in Grp1Prop is >= this proportion then characterise as homozygous
					szAssembRefFile,	// contains original reference assembly, fasta or PBA, against which samples were aligned
					szChromFile,		// BED file, contains chromosome names and sizes
					szDGTsFile,			// file containing DGT loci and allele sample groupings
					szQTLsFile,			// file containing QTL allele loci
					NumPBAFiles,		// number of input PBA file specs, these are the source PBA files used during the grouping generation
					pszPBAFiles,		// names of input PBA files (wildcards allowed)
					NumIncludeChroms,	// number of chromosome regular expressions to include
					pszIncludeChroms,	// array of include chromosome regular expressions
					NumExcludeChroms,	// number of chromosome expressions to exclude
					pszExcludeChroms,	// array of exclude chromosome regular expressions
					szRsltsFileBaseName,		// write results to this file base name - will be suffixed with result types
					NumThreads);				// number of worker threads to use
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

int Process(eModeDGTA PMode,		// processing mode
	double MinCoverage,				// if coverage < this threshold then class as being low coverage
	double HomozPropThres,			// if proportion of samples in Grp1Prop is >= this proportion then characterise as homozygous	char* pszAssembRefFile,	// contains orginal reference assembly, diplotype PBA, against which samples were aligned
	char *pszAssembRefFile,			// contains original reference assembly, PBA, against which samples were aligned
	char* pszChromFile,				// BED file, contains chromosome names and sizes
	char* pszDGTsFile,				// file containing DGT loci and allele sample groupings
	char* pszQTLsFile,				// file containing QTL allele loci
	int NumPBAFiles,				// number of input PBA file specs, these are the source PBA files used during the grouping generation
	char** ppszPBAFiles,			// names of input PBA files (wildcards allowed)
	int	NumIncludeChroms,			// number of chromosome regular expressions to include
	char** ppszIncludeChroms,		// array of include chromosome regular expressions
	int	NumExcludeChroms,			// number of chromosome expressions to exclude
	char** ppszExcludeChroms,		// array of exclude chromosome regular expressions
	char* pszRsltsFileBaseName,		// write results to this file base name - will be suffixed with result types
	int NumThreads)					// number of worker threads to use
{
int Rslt;
CDGTvQTLs* pDGTvQTLs;

if ((pDGTvQTLs = new CDGTvQTLs) == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to instantiate instance of CDGTvQTLs");
	return(eBSFerrObj);
	}
Rslt = pDGTvQTLs->Process(PMode, MinCoverage, HomozPropThres,pszAssembRefFile, pszChromFile, pszDGTsFile, pszQTLsFile, NumPBAFiles, ppszPBAFiles, NumIncludeChroms, ppszIncludeChroms, NumExcludeChroms, ppszExcludeChroms, pszRsltsFileBaseName, NumThreads);
delete pDGTvQTLs;
return(Rslt);
}

CDGTvQTLs::CDGTvQTLs()
{
m_hOutFile = -1;
m_hInFile = -1;
m_pChromMetadata = nullptr;
m_pDGTQTLAlleles = nullptr;
m_pInBuffer = nullptr;
m_pszOutBuffer = nullptr;
m_pBedFile = nullptr;
m_pInCSVFile = nullptr;
m_bMutexesCreated = false;
Reset();
}

CDGTvQTLs::~CDGTvQTLs()
{
if (m_hInFile != -1)
	close(m_hInFile);
if (m_hOutFile != -1)
	close(m_hOutFile);

if (m_pszOutBuffer != nullptr)
	delete[]m_pszOutBuffer;
if (m_pInBuffer != nullptr)
	delete[]m_pInBuffer;

if (m_pBedFile != nullptr)
	delete m_pBedFile;

if(m_pInCSVFile != nullptr)
	delete m_pInCSVFile;

if (m_pChromMetadata != nullptr)
	{
	tsCHChromMetadata* pChromMetadata = m_pChromMetadata;
	for (int32_t Idx = 0; Idx < m_UsedNumChromMetadata; Idx++, pChromMetadata++)
		{
		if (pChromMetadata->pPBAs != nullptr)
#ifdef _WIN32
			free(pChromMetadata->pPBAs);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
		if (pChromMetadata->pPBAs != MAP_FAILED)
			munmap(pChromMetadata->pPBAs, pChromMetadata->ChromLen);
#endif
		}

#ifdef _WIN32
	free(m_pChromMetadata);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (m_pChromMetadata != MAP_FAILED)
		munmap(m_pChromMetadata, m_AllocdChromMetadataMem);
#endif
	}

if (m_pDGTQTLAlleles != nullptr)
	{
#ifdef _WIN32
	free(m_pDGTQTLAlleles);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (m_pDGTQTLAlleles != MAP_FAILED)
		munmap(m_pDGTQTLAlleles, m_AllocDGTQTLAllelesMem);
#endif
	}
}

void
CDGTvQTLs::Reset(void)						// resets class instance state back to that immediately following instantiation
{
if (m_hInFile != -1)
	{
	close(m_hInFile);
	m_hInFile = -1;
	}
if (m_hOutFile != -1)
	{
	close(m_hOutFile);
	m_hOutFile = -1;
	}

if (m_pInBuffer != nullptr)
	{
	delete[]m_pInBuffer;
	m_pInBuffer = nullptr;
	}

if (m_pszOutBuffer != nullptr)
	{
	delete[]m_pszOutBuffer;
	m_pszOutBuffer = nullptr;
	}

if (m_pBedFile != nullptr)
	{
	delete m_pBedFile;
	m_pBedFile = nullptr;
	}

if (m_pInCSVFile != nullptr)
	{
	delete m_pInCSVFile;
	m_pInCSVFile = nullptr;
	}

if (m_pChromMetadata != nullptr)
	{
	tsCHChromMetadata* pChromMetadata = m_pChromMetadata;
	for (int32_t Idx = 0; Idx < m_UsedNumChromMetadata; Idx++, pChromMetadata++)
		{
		if (pChromMetadata->pPBAs != nullptr)
#ifdef _WIN32
			free(pChromMetadata->pPBAs);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
		if (pChromMetadata->pPBAs != MAP_FAILED)
			munmap(pChromMetadata->pPBAs, pChromMetadata->ChromLen);
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

if(m_pDGTQTLAlleles != nullptr)
	{
#ifdef _WIN32
	free(m_pDGTQTLAlleles);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (m_pDGTQTLAlleles != MAP_FAILED)
		munmap(m_pDGTQTLAlleles, m_AllocDGTQTLAllelesMem);
#endif
	m_pDGTQTLAlleles = nullptr;
	}
m_UsedDGTQTLAlleles = 0;
m_AllocDGTQTLAlleles = 0;
m_AllocDGTQTLAllelesMem = 0;
memset(m_DGTQTLAllelesHashes,0,sizeof(m_DGTQTLAllelesHashes));


m_UsedNumChromMetadata = 0;
m_AllocdChromMetadata = 0;
m_AllocdChromMetadataMem = 0;
m_PMode = eDGTADefault;
m_FndrTrim5 = 0;
m_FndrTrim3 = 0;
m_ProgTrim5 = 0;
m_ProgTrim3 = 0;
m_LAReadsetNameID = 0;
m_NumReadsetNames = 0;
m_NxtszReadsetIdx = 0;
m_szReadsetNames[0] = '\0';
m_LAChromNameID = 0;
m_NumChromNames = 0;
m_NxtszChromIdx = 0;
m_szChromNames[0] = '\0';
m_InNumProcessed = 0;
m_InNumBuffered = 0;
m_AllocInBuff = 0;
m_OutBuffIdx = 0;
m_AllocOutBuff = 0;
m_InFileOfs = 0;
memset(m_ChromSizes, 0, sizeof(m_ChromSizes));
m_TotChromSizes = 0;
m_NumChromSizes = 0;

m_MinCoverage = cDfltMinCoverage;
m_HomozPropThres = cDfltHomozPropThres;

memset(m_Fndrs2Proc, 0, sizeof(m_Fndrs2Proc));
m_NumFounders = 0;

if (m_bMutexesCreated)
	DeleteMutexes();
m_bMutexesCreated = false;
}

int
CDGTvQTLs::CreateMutexes(void)
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
CDGTvQTLs::DeleteMutexes(void)
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
CDGTvQTLs::AcquireSerialise(void)
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
CDGTvQTLs::ReleaseSerialise(void)
{
#ifdef _WIN32
ReleaseMutex(m_hSerialiseAccess);
#else
pthread_mutex_unlock(&m_hSerialiseAccess);
#endif
}


void
CDGTvQTLs::AcquireFastSerialise(void)
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
CDGTvQTLs::ReleaseFastSerialise(void)
{
#ifdef _WIN32
InterlockedCompareExchange(&m_FastSerialise, 0, 1);
#else
__sync_val_compare_and_swap(&m_FastSerialise, 1, 0);
#endif
}


int32_t		// returned chrom identifier, 0 if unable to accept this chromosome name
CDGTvQTLs::AddChrom(char* pszChrom) // associate unique identifier with this chromosome name
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


int32_t		// returned chrom identifier, 0 if unable to locate this chromosome name
CDGTvQTLs::LocateChrom(char* pszChrom) // return unique identifier associated with this chromosome name
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

bool					// true if chrom is accepted, false if chrom not accepted
CDGTvQTLs::AcceptThisChromID(uint32_t ChromID)
{
	char* pzChrom;
	if ((pzChrom = LocateChrom(ChromID)) == nullptr)
		return(false);
	return(AcceptThisChromName(pzChrom));
}

bool					// true if chrom is accepted, false if chrom not accepted
CDGTvQTLs::AcceptThisChromName(char* pszChrom,   // chromosome name
	bool bKnown)	// if true then chromosome must have been previously processed and accepted by LoadChromSizes() processing
{
bool bMatch;
if (bKnown)
	return(LocateChrom(pszChrom) < 1 ? false : true);
AcquireSerialise();
bMatch = !m_RegExprs.MatchExcludeRegExpr(pszChrom);
if (bMatch)
	bMatch = m_RegExprs.MatchIncludeRegExpr(pszChrom);
ReleaseSerialise();
return(bMatch);
}

uint8_t
CDGTvQTLs::LocateReadsetChromLociAlleles(int32_t ReadsetID,	// return alleles for this readset 
	int32_t ChromID,		// on this chromosome
	int32_t Loci)		// at this loci
{
	uint8_t* pPBA;
	if ((pPBA = LocatePBAfor(ReadsetID, ChromID)) == nullptr)
		return(0);

	return(pPBA[Loci]);
}

uint8_t*								// returned pointer to start of PBA
CDGTvQTLs::LocatePBAfor(int32_t ReadSetID,		// readset identifier 
	int32_t ChromID)			// chrom identifier
{
	tsCHReadsetMetadata* pReadsetMetadata;
	tsCHChromMetadata* pChromMetadata;
	int32_t CurChromMetadataIdx;

	if (ReadSetID > m_NumReadsetNames || ReadSetID == 0)
		return(nullptr);
	pReadsetMetadata = &m_Readsets[ReadSetID - 1];
	CurChromMetadataIdx = pReadsetMetadata->StartChromMetadataIdx;
	for (int32_t ChromIdx = 0; ChromIdx < pReadsetMetadata->NumChroms && CurChromMetadataIdx != 0; ChromIdx++)
	{
		pChromMetadata = &m_pChromMetadata[CurChromMetadataIdx - 1];
		if (pChromMetadata->ChromID == ChromID)
			return(pChromMetadata->pPBAs);
		CurChromMetadataIdx = pChromMetadata->NxtChromMetadataIdx;
	}
	return(nullptr);
}

tsCHChromMetadata*								// returned pointer to chromosome metadata
CDGTvQTLs::LocateChromMetadataFor(int32_t ReadSetID,		// readset identifier 
	int32_t ChromID)			// chrom identifier
{
	tsCHReadsetMetadata* pReadsetMetadata;
	tsCHChromMetadata* pChromMetadata;
	int32_t CurChromMetadataIdx;

	if (ReadSetID > m_NumReadsetNames || ReadSetID == 0)
		return(nullptr);
	pReadsetMetadata = &m_Readsets[ReadSetID - 1];
	CurChromMetadataIdx = pReadsetMetadata->StartChromMetadataIdx;
	for (int32_t ChromIdx = 0; ChromIdx < pReadsetMetadata->NumChroms && CurChromMetadataIdx != 0; ChromIdx++)
	{
		pChromMetadata = &m_pChromMetadata[CurChromMetadataIdx - 1];
		if (pChromMetadata->ChromID == ChromID)
			return(pChromMetadata);
		CurChromMetadataIdx = pChromMetadata->NxtChromMetadataIdx;
	}
	return(nullptr);
}

// NOTE: Readsets are checked for uniqueness as readsets must be unique within a given readetset type
int32_t		// returned readset identifier, 0 if unable to accept this readset name
CDGTvQTLs::AddReadset(char* pszReadset, // associate unique identifier with this readset name
	uint8_t ReadsetType)	// 0: founder, 1: progeny, 2: control
{
	int32_t ReadsetNameIdx;
	int ReadsetNameLen;
	char Type;
	char* pszLAname;
	Type = '0' + (char)ReadsetType;
	if (m_NumReadsetNames == 0)
	{
		memset(m_NumReadsetTypes, 0, sizeof(m_NumReadsetTypes));
		memset(m_FndrIDMappings, 0, sizeof(m_FndrIDMappings));
	}

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
	if (m_NumReadsetNames == (cMaxPBAFiles + 1))
		return(0);		// unable to hold any more readsets, treating as if a non-unique!

	m_szReadsetIdx[m_NumReadsetNames++] = m_NxtszReadsetIdx;
	pszLAname = &m_szReadsetNames[m_NxtszReadsetIdx];
	*pszLAname++ = Type;
	strcpy(pszLAname, pszReadset);
	m_NxtszReadsetIdx += ReadsetNameLen + 2;
	m_LAReadsetNameID = m_NumReadsetNames;
	m_NumReadsetTypes[ReadsetType]++;
	if (ReadsetType == 0)			// founder types are special, need to maintain mappings of readset identifiers
		m_FndrIDMappings[m_NumReadsetTypes[ReadsetType] - 1] = m_LAReadsetNameID;
	return(m_LAReadsetNameID);
}

int32_t	// returned readset identifier index for founder type readsets, -1 if unable to locate readset identifier
CDGTvQTLs::LocateReadsetIDIdx(int32_t ReadsetID,		// requiring idx for this readset identifier
	uint8_t ReadsetType)	// 0: founder, 1: progeny, 2: control
{
	int32_t ReadsetIDIdx;
	if (ReadsetType != 0 || m_NumReadsetNames == 0 || m_NumReadsetTypes[ReadsetType] == 0)
		return(-1);
	// currently a linear search but expecting only a few hundred readset identifers to be founder types so shouldn't be nuch of an overhead
	for (ReadsetIDIdx = 0; ReadsetIDIdx < m_NumReadsetTypes[ReadsetType]; ReadsetIDIdx++)
		if (ReadsetID == m_FndrIDMappings[ReadsetIDIdx])
			return(ReadsetIDIdx);
	return(-1);
}

int32_t		// returned Readset identifier, 0 if unable to locate this Readset name
CDGTvQTLs::LocateReadset(char* pszReadset, // return unique identifier associated with this Readset name
	uint8_t ReadsetType)	// 0: founder, 1: progeny, 2: control
{
	int32_t ReadsetNameIdx;
	char Type;
	char* pszLAReadset;

	if (m_NumReadsetNames > 0 && m_NumReadsetTypes[ReadsetType] == 0)
		return(0);

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
CDGTvQTLs::LocateReadset(int32_t ReadsetID)
{
	int32_t Idx;
	Idx = ReadsetID & 0x00ffffff;			// mask out any potential ReadsetType
	if (Idx < 1 || Idx > m_NumReadsetNames)
		return(nullptr);
	return(&(m_szReadsetNames[m_szReadsetIdx[Idx - 1] + 1])); // skipping lead char which is the ReadsetType
}


char*
CDGTvQTLs::LocateChrom(int32_t ChromID)
{
if (ChromID < 1 || ChromID > m_NumChromNames)
	return(nullptr);
return(&m_szChromNames[m_szChromIdx[ChromID - 1]]);
}

// loading BED which specifies chrom names and sizes
int		// returning number of chromosomes parsed from BED file and accepted after filtering for wildcards
CDGTvQTLs::LoadChromSizes(char* pszBEDFile) // BED file containing chromosome names and sizes
{
int Rslt;
uint32_t ChromID;
int CurFeatureID;
char szFeatName[cMaxDatasetSpeciesChrom];
char szChromName[cMaxDatasetSpeciesChrom];
int32_t StartLoci;
int32_t EndLoci;
int NumChroms;
int64_t TotChromSizes;


if ((m_pBedFile = new CBEDfile) == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to instantiate CBEDfile");
	return(eBSFerrObj);
	}

if ((Rslt = m_pBedFile->Open(pszBEDFile, eBTAnyBed)) != eBSFSuccess)
	{
	while (m_pBedFile->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal, gszProcName, m_pBedFile->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to open '%s' for processing", pszBEDFile);
	return(eBSFerrOpnFile);
	}

m_TotChromSizes = 0;
m_NumChromSizes = 0;
NumChroms = 0;
TotChromSizes = 0;
CurFeatureID = 0;
while (Rslt == eBSFSuccess && (CurFeatureID = m_pBedFile->GetNextFeatureID(CurFeatureID)) > 0)
	{
	m_pBedFile->GetFeature(CurFeatureID,	// feature instance identifier
		szFeatName,				// where to return feature name
		szChromName,			// where to return chromosome name
		&StartLoci,				// where to return feature start on chromosome (0..n) 
		&EndLoci);				// where to return feature end on chromosome

	// is this chromosome to be be filtered out?
	if (!AcceptThisChromName(szChromName, false))
		{
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "LoadChromSizes: filtered out chromosome '%s'", szChromName);
		continue;
		}

	if (StartLoci != 0)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadChromSizes: unable to accept chromosome '%s' because start loci not 0", szChromName);
		delete m_pBedFile;
		m_pBedFile = nullptr;
		return(eBSFerrChrom);
		}

	if (EndLoci < 0 || (size_t)EndLoci >= cAllocPackedBaseAlleles)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadChromSizes: unable to accept chromosome '%s' because end loci not in acceptable range, must be in range 0..%zd", szChromName, cAllocPackedBaseAlleles);
		delete m_pBedFile;
		m_pBedFile = nullptr;
		return(eBSFerrChrom);
		}

	// expecting a single instance of each chromosome 
	if (LocateChrom(szChromName) != 0)	// must be a unique chromosome name
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadChromSizes: multiple instances of chromosome '%s', chromosome names must be unique", szChromName);
		delete m_pBedFile;
		m_pBedFile = nullptr;
		return(eBSFerrChrom);
		}
	// ensure not about to overflow m_ChromSizes[] capacity
	if (NumChroms == cMaxChromNames)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadChromSizes: chromosome '%s' - can only accept a maximum of %d from '%s'", szChromName, NumChroms, pszBEDFile);
		delete m_pBedFile;
		m_pBedFile = nullptr;
		return(eBSFerrChrom);
		}
	ChromID = AddChrom(szChromName);		// new chromosome
	m_ChromSizes[ChromID - 1] = EndLoci + 1;	// record it's size (note that EndLoci is actually inclusive, so need to add 1 to obtain length!)
	TotChromSizes += (int64_t)EndLoci + 1;
	NumChroms += 1;
	}
m_NumChromSizes = NumChroms;
m_TotChromSizes = TotChromSizes;
gDiagnostics.DiagOut(eDLInfo, gszProcName, "LoadChromSizes: accepted %d chromosome sizes from '%s'", NumChroms, pszBEDFile);
delete m_pBedFile;
m_pBedFile = nullptr;
return(NumChroms);
}

uint32_t				// returns number of unprocessed bytes in buffer
CDGTvQTLs::FillInBuffer(uint32_t MinRequired, uint32_t MaxRequired) // try and fill input buffer with at least MinRequired or refill (if < MinRequired) up to MaxRequired (if 0 then max buffer allocation)
{
	int NumRead;
	uint32_t UnProcessed;
	if (MaxRequired == 0)
		MaxRequired = m_AllocInBuff;
	else
		if (MaxRequired < MinRequired)
			MaxRequired = MinRequired;
	UnProcessed = m_InNumBuffered - m_InNumProcessed;

	// if already filled to MinRequired then no further action required
	if (MinRequired <= UnProcessed)
		return(UnProcessed);

	// copy down any bytes which are yet to be processed
	if (UnProcessed) {
		memmove(m_pInBuffer, &m_pInBuffer[m_InNumProcessed], UnProcessed);
		m_InNumBuffered = UnProcessed;
		m_InNumProcessed = 0;
	}
	else
	{
		m_InNumBuffered = 0;
		m_InNumProcessed = 0;
	}

	// attempt to fill input buffer to it's MaxRequired
	do {
		if ((NumRead = read(m_hInFile, &m_pInBuffer[m_InNumProcessed], MaxRequired - m_InNumBuffered)) < 0)
		{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error %s attempting to read from file", strerror(errno));
			return(0);
		}
		m_InNumBuffered += (uint32_t)NumRead;
	} while (NumRead > 0 && m_InNumBuffered < MinRequired);
	m_InFileOfs = _lseeki64(m_hInFile, 0, SEEK_CUR);             // record file offset at which next file read will start from
	return(m_InNumBuffered);
}

int
CDGTvQTLs::AllocChromMetadata(void)
{
uint32_t ToAllocdChromMetadata;
tsCHChromMetadata* pChromMetadata;
size_t memreq;
if (m_pChromMetadata == nullptr)					// may be nullptr first time in
{
	memreq = cAllocChromMetadata * sizeof(tsCHChromMetadata);
#ifdef _WIN32
	m_pChromMetadata = (tsCHChromMetadata*)malloc((size_t)memreq);
	if (m_pChromMetadata == nullptr)
	{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AllocChromMetadata: Memory allocation of %zd bytes failed", (int64_t)memreq);
		return(eBSFerrMem);
	}
#else
	m_pChromMetadata = (tsCHChromMetadata*)mmap(nullptr, (size_t)memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
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
		size_t memreq = ToAllocdChromMetadata * sizeof(tsCHChromMetadata);
#ifdef _WIN32
		pChromMetadata = (tsCHChromMetadata*)realloc(m_pChromMetadata, memreq);
		if (pChromMetadata == nullptr)
		{
#else
		pChromMetadata = (tsCHChromMetadata*)mremap(m_pChromMetadata, m_AllocdChromMetadataMem, memreq, MREMAP_MAYMOVE);
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

uint8_t*
CDGTvQTLs::AllocPBAs(int32_t ChromLen)	// allocate memory to hold at least this many packed base alleles
{
uint8_t* pPBAs;
size_t memreq;
memreq = (size_t)ChromLen;	// no safety margin!
#ifdef _WIN32
pPBAs = (uint8_t*)malloc((size_t)memreq);
if (pPBAs == nullptr)
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "AllocPBAs Memory allocation of %zd bytes failed", (int64_t)memreq);
#else
pPBAs = (uint8_t*)mmap(nullptr, (size_t)memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
if (pPBAs == MAP_FAILED)
{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "AllocPBAs: Memory allocation of %zd bytes through mmap()  failed", (int64_t)memreq, strerror(errno));
	pPBAs = nullptr;
}
#endif
return(pPBAs);
}


int32_t					// returned readset identifier (1..n) or < 0 if errors
CDGTvQTLs::LoadPBAFile(char* pszFile,	// load chromosome metadata and PBA data from this file, filters out chroms - AcceptThisChromName()
	uint8_t ReadsetType,	// 0: founder, 1: progeny, 2: control
	bool bChromMetaOnly)  // load chrom metadata (chrom name,length, file offset at which chrom PBAs start) but don't actually load the chromosome PBAs
{
	int Rslt;
	int Version;
	char szExperimentID[100];
	char szRefAssemblyID[100];
	char szReadsetID[100];
	uint32_t PrevChromMetadataIdx;
	uint32_t ChromID;
	uint32_t ReadsetID;
	tsCHReadsetMetadata* pReadsetMetadata;
	tsCHChromMetadata* pChromMetadata;
	tsCHChromMetadata* pPrevChromMetadata;
	int64_t FileOfsPBA;
	uint8_t* pBuff;
	int scanlen;
	int NumTags;

	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Loading '%s' file", pszFile);
#ifdef _WIN32
	m_hInFile = open(pszFile, O_READSEQ);		// file access is normally sequential..
#else
	m_hInFile = open64(pszFile, O_READSEQ);		// file access is normally sequential..
#endif
	if (m_hInFile == -1)							// check if file open succeeded
	{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to open input file '%s' : %s", pszFile, strerror(errno));
		return(eBSFerrOpnFile);
	}

	// attempt to load the readset metadata 
	m_InNumBuffered = 0;
	m_InNumProcessed = 0;
	if ((FillInBuffer(500, min(m_AllocInBuff, 500u)) == 0) || m_InNumBuffered < 9) // 500 will cover the maximally sized header
	{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to read at least a partial header from input file '%s'", pszFile);
		return(eBSFerrOpnFile);
	}


	// check file type is PBA
	if (strncmp((char*)m_pInBuffer, "Type:PbA\n", 9))
	{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Input file '%s' exists, unable to parse file type header tag, is it a packed base allele file", pszFile);
		return(eBSFerrOpnFile);
	}
	m_InNumProcessed = 9;

	if (m_InNumBuffered < 500)
	{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Input file '%s' exists, unable to read complete header and initial chromosome metadata, expecting file size >= 500, is it a packed base allele file", pszFile);
		return(eBSFerrOpnFile);
	}

	// parse out tagnames and associated values
	pBuff = &m_pInBuffer[m_InNumProcessed];
	NumTags = sscanf((char*)pBuff, "Version:%d\nExperimentID:%99[^\n]\nReferenceID:%99[^\n]\nReadsetID:%99[^\n]%n", &Version, szExperimentID, szRefAssemblyID, szReadsetID, &scanlen);
	if (NumTags != 4 || Version != 1 || szExperimentID[0] == '\0' || szRefAssemblyID[0] == '\0' || szReadsetID[0] == '\0')
	{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Input file '%s' exists as a packed base allele file but inconsistencies in header tag values", pszFile);
		return(eBSFerrOpnFile);
	}
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Metadata ExperimentID:%s, ReferenceID:%s, ReadsetID:%s", szExperimentID, szRefAssemblyID, szReadsetID);

	m_InNumProcessed += scanlen + 1;		// header tags were '\n' terminated except for final which was '\0' terminated

	m_InFileOfs = m_InNumProcessed;         // 1st chromosome starts immediately following the header metadata
	if (_lseeki64(m_hInFile, m_InFileOfs, SEEK_SET) != m_InFileOfs)
	{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Input file '%s' exists as a packed base allele file but unable to seek past header metadata", pszFile);
		return(eBSFerrOpnFile);
	}

	if ((ReadsetID = AddReadset(szReadsetID, ReadsetType)) == 0)
	{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Input file '%s' duplicates the ReadsetID '%s' of a previously loaded readset", pszFile, szReadsetID);
		return(eBSFerrOpnFile);
	}

	pReadsetMetadata = &m_Readsets[ReadsetID - 1];
	memset(pReadsetMetadata, 0, sizeof(*pReadsetMetadata));
	pReadsetMetadata->ReadsetType = ReadsetType;
	pReadsetMetadata->NumChroms = 0;
	strcpy(pReadsetMetadata->szFileName, pszFile);
	strcpy(pReadsetMetadata->szExperimentID, szExperimentID);
	strcpy(pReadsetMetadata->szRefAssemblyID, szRefAssemblyID);
	pReadsetMetadata->ReadsetID = ReadsetID;
	pReadsetMetadata->StartChromID = 0;
	pReadsetMetadata->StartChromMetadataIdx = 0;
	pReadsetMetadata->NxtFileChromOfs = m_InFileOfs;

	int ChromNameLen;
	char* pszChromName;
	int32_t ChromLen;

	// iterate over all chromosomes
	m_InNumProcessed = 0;
	m_InNumBuffered = 0;
	PrevChromMetadataIdx = 0;
	while (FillInBuffer((uint32_t)110, 110) == 110) // reading chromosome metadata
	{
		pBuff = &m_pInBuffer[m_InNumProcessed++];
		ChromNameLen = (int)*pBuff++;
		pszChromName = (char*)pBuff;
		pBuff += (int64_t)1 + ChromNameLen;
		ChromLen = *(int32_t*)pBuff;
		FileOfsPBA = pReadsetMetadata->NxtFileChromOfs + ChromNameLen + 6;
		pReadsetMetadata->NxtFileChromOfs += (int64_t)ChromNameLen + 6 + ChromLen;
		// check if this chromosome is to be retained for further processing
		if (!AcceptThisChromName(pszChromName))
		{
			// not accepting this chromosome
			m_InFileOfs = pReadsetMetadata->NxtFileChromOfs;         // skip over curent chrom's PBAs to start of next chromosome
			if ((pReadsetMetadata->NxtFileChromOfs = _lseeki64(m_hInFile, m_InFileOfs, SEEK_SET)) != m_InFileOfs)
				break;
			m_InNumProcessed = 0;
			m_InNumBuffered = 0;
			continue;
		}

		// before accepting chrom then ensure that it's PBA length matches the BED chromosome sizes
		ChromID = AddChrom(pszChromName);
		if (ChromLen != m_ChromSizes[ChromID - 1])
		{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Input file '%s' with readset '%s' has chromosome '%s' size mismatch - expected size was %d, actual size %d", pszFile, szReadsetID, pszChromName, m_ChromSizes[ChromID - 1], ChromLen);
			return(eBSFerrChrom);
		}

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
		pReadsetMetadata->NumChroms++;
		PrevChromMetadataIdx = m_UsedNumChromMetadata;
		pChromMetadata->ChromID = ChromID;
		pChromMetadata->ChromMetadataIdx = m_UsedNumChromMetadata;
		pChromMetadata->NxtChromMetadataIdx = 0;
		pChromMetadata->ChromLen = ChromLen;
		pChromMetadata->ReadsetID = ReadsetID;
		pChromMetadata->HGBinID = 0;
		pChromMetadata->FileOfsPBA = FileOfsPBA;
		pChromMetadata->pPBAs = nullptr;
		if (bChromMetaOnly) // loading chrom metadata only without actually loading the chromosome PBAs?
		{
			m_InFileOfs = pReadsetMetadata->NxtFileChromOfs;         // skip over curent chrom's PBAs to start of next chromosome
			if ((pReadsetMetadata->NxtFileChromOfs = _lseeki64(m_hInFile, m_InFileOfs, SEEK_SET)) != m_InFileOfs)
				break;
			m_InNumProcessed = 0;
			m_InNumBuffered = 0;
			continue;
		}
		// loading PBAs
		pChromMetadata->pPBAs = AllocPBAs(ChromLen);
		m_InNumProcessed += ChromNameLen + 5;
		FillInBuffer((uint32_t)ChromLen, ChromLen);
		pBuff = &m_pInBuffer[m_InNumProcessed];


		memcpy(pChromMetadata->pPBAs, pBuff, ChromLen);
		m_InNumProcessed += ChromLen;
		//	validate PBA allele composition, earlier releases were retaining very low allele proportions so validating and in-place replacing these as being non-alignments
		int NumErrs = ValidatePBAs(ChromLen, pChromMetadata->pPBAs, true, true);
		//
		if (ReadsetType == 0)
			TrimPBAs(m_FndrTrim5, m_FndrTrim3, ChromLen, pChromMetadata->pPBAs);
		else // treating controls as if progeny when trimming
			TrimPBAs(m_ProgTrim5, m_ProgTrim3, ChromLen, pChromMetadata->pPBAs);
	}

	close(m_hInFile);
	m_hInFile = -1;
	return((int32_t)ReadsetID);
}

bool  // false: no previous memory allocated for containing the chromosome PBAs, true: allocation was present and has been deleted
CDGTvQTLs::DeleteSampleChromPBAs(int32_t SampleID,   // Sample identifier
	int32_t ChromID)    // chrom identifier
{
	tsCHChromMetadata* pChromMetadata;

	// returned pointer to chromosome metadata
	if ((pChromMetadata = LocateChromMetadataFor(SampleID, ChromID)) == nullptr)
		return(0);
	if (pChromMetadata->pPBAs == nullptr)
		return(false);
#ifdef _WIN32
	free(pChromMetadata->pPBAs);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (pChromMetadata->pPBAs != MAP_FAILED)
		munmap(pChromMetadata->pPBAs, pChromMetadata->ChromLen);
#endif
	pChromMetadata->pPBAs = nullptr;
	return(true);
}

uint8_t*  // returned PBA or nullptr if unable to load PBAs for requested sample.chrom 
CDGTvQTLs::LoadSampleChromPBAs(int32_t SampleID,   // Sample identifier
	int32_t ChromID)    // chrom identifier specifying which PBAs is to be loaded from SampleID file
{
	int hInFile;
	char* pszChrom;
	int64_t ChromSeekOfs;
	tsCHChromMetadata* pChromMetadata;
	tsCHReadsetMetadata* pReadsetMetadata;
	int32_t NumRead;
	int32_t NumLoaded;

	// returned pointer to chromosome metadata
	if ((pChromMetadata = LocateChromMetadataFor(SampleID, ChromID)) == nullptr)
		return(nullptr);

	pszChrom = LocateChrom(ChromID);

	if (pChromMetadata->pPBAs != nullptr) // PBAs may already have been loaded for this chromosome, delete as may be stale
	{
#ifdef _WIN32
		free(pChromMetadata->pPBAs);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
		if (pChromMetadata->pPBAs != MAP_FAILED)
			munmap(pChromMetadata->pPBAs, pChromMetadata->ChromLen);
#endif
		pChromMetadata->pPBAs = nullptr;
	}

	// need to actually load from file
	pReadsetMetadata = &m_Readsets[SampleID - 1];
#ifdef _WIN32
	hInFile = open(pReadsetMetadata->szFileName, O_READSEQ);		// file access is normally sequential..
#else
	hInFile = open64(pReadsetMetadata->szFileName, O_READSEQ);		// file access is normally sequential..
#endif
	if (hInFile == -1)							// check if file open succeeded
	{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadSampleChromPBA: Unable to open input file '%s' : %s", pReadsetMetadata->szFileName, strerror(errno));
		return(nullptr);
	}
	if ((ChromSeekOfs = _lseeki64(hInFile, pChromMetadata->FileOfsPBA, SEEK_SET)) != pChromMetadata->FileOfsPBA)
	{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Input file '%s' opened, _lseek() to file offset %zd failed for chrom '%s', returned %zd, error: '%s'",
			pReadsetMetadata->szFileName, pChromMetadata->FileOfsPBA, ChromSeekOfs, pszChrom, strerror(errno));
		close(hInFile);
		hInFile = -1;
		return(nullptr);
	}
	pChromMetadata->pPBAs = AllocPBAs(pChromMetadata->ChromLen);
	// attempt to load Chromosomes PBAs from file
	NumLoaded = 0;
	do {
		if ((NumRead = read(hInFile, &pChromMetadata->pPBAs[NumLoaded], pChromMetadata->ChromLen - NumLoaded)) < 0)
		{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error %s attempting to read chromosome '%s' from file '%s'", strerror(errno), pszChrom, pReadsetMetadata->szFileName);
			close(hInFile);
			hInFile = -1;
			if (pChromMetadata->pPBAs == nullptr)
				return(nullptr);
#ifdef _WIN32
			free(pChromMetadata->pPBAs);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
			if (pChromMetadata->pPBAs != MAP_FAILED)
				munmap(pChromMetadata->pPBAs, pChromMetadata->ChromLen);
#endif
			pChromMetadata->pPBAs = nullptr;
			return(nullptr);
		}
		NumLoaded += NumRead;
	} while (NumRead > 0 && NumLoaded < pChromMetadata->ChromLen);
	close(hInFile);
	hInFile = -1;
	//	validate PBA allele composition, earlier releases of 'kalign' were retaining very low allele proportions so validating and in-place replacing these as being non-alignments
	int NumErrs = ValidatePBAs(pChromMetadata->ChromLen, pChromMetadata->pPBAs, true, true);
	//
	if (pReadsetMetadata->ReadsetType == 0)
		TrimPBAs(m_FndrTrim5, m_FndrTrim3, pChromMetadata->ChromLen, pChromMetadata->pPBAs);
	else // treating controls as if progeny when trimming
		TrimPBAs(m_ProgTrim5, m_ProgTrim3, pChromMetadata->ChromLen, pChromMetadata->pPBAs);
	return(pChromMetadata->pPBAs);
}

// trim 5' and 3' aligned segments within the PBAs attempting to reduce sequencing error induced false alleles
int			// returns number of non-trimmed loci in the pPBAs
CDGTvQTLs::TrimPBAs(uint32_t Trim5,	// trim 5' this many aligned PBA bases from each aligned segment
	uint32_t Trim3,	// trim 3' this many aligned PBA bases from each aligned segment
	uint32_t PBALen,	// pPBAs contains this many packed base alleles
	uint8_t* pPBAs)		// packed base alleles to be processed
{
	uint32_t NonTrimmed;
	uint32_t ToTrim;
	uint32_t Loci;
	uint8_t* pPBA;
	if (pPBAs == nullptr || PBALen == 0)
		return(0);
	if (Trim5 == 0 && Trim3 == 0)
		return(PBALen);
	NonTrimmed = PBALen;

	// first the 5' trim
	if (Trim5)
	{
		ToTrim = Trim5;
		pPBA = pPBAs;
		for (Loci = 0; Loci < PBALen; Loci++, pPBA++)
		{
			if (*pPBA == 0)
			{
				ToTrim = Trim5;
				continue;
			}
			if (ToTrim)
			{
				ToTrim -= 1;
				NonTrimmed -= 1;
				*pPBA = 0;
			}
		}
	}
	// then the 3' trim
	if (Trim3)
	{
		ToTrim = Trim3;
		pPBA = &pPBAs[PBALen - 1];
		for (Loci = 0; Loci < PBALen; Loci++, pPBA--)
		{
			if (*pPBA == 0)
			{
				ToTrim = Trim3;
				continue;
			}
			if (ToTrim)
			{
				NonTrimmed -= 1;
				ToTrim -= 1;
				*pPBA = 0;
			}
		}
	}
	return(NonTrimmed);
}

int			// returned number of PBAs which are non-conformant
CDGTvQTLs::ValidatePBAs(int32_t Length,
	uint8_t* pPBAs,
	bool bSetNoAlleles, // if non-conformant then overwrite *pAlleles to be no alleles present
	bool bNormalise)    // normalise for coverage, low coverage major alleles (2,0,0,0) are normalised to high coverage (3,0,0,0)
{
	int32_t Idx;
	int32_t NumErrors = 0;
	for (Idx = 0; Idx < Length; Idx++, pPBAs++)
	{
		if (!ValidatePBA(pPBAs, bSetNoAlleles, bNormalise))
			NumErrors++;
	}
	return(NumErrors);
}

bool
CDGTvQTLs::ValidatePBA(uint8_t* pAlleles,	// validate that PBA alleles are properly conformant
	bool bSetNoAlleles, // if non-conformant then overwrite *pAlleles to be no alleles present
	bool bNormalise)    // normalise for coverage, low coverage major alleles (2,0,0,0) are normalised to high coverage (3,0,0,0)
{
	uint32_t AlleleIdx;
	uint8_t Allele;
	uint8_t AlleleMsk;
	uint8_t LoMinorAllele;
	int32_t Num2s;
	int32_t Num1s;

	if (*pAlleles == 0)
		return(true);

	Num2s = Num1s = 0;
	AlleleMsk = 0x03;
	for (AlleleIdx = 0; AlleleIdx < 4; AlleleIdx++, AlleleMsk <<= 2)
	{
		Allele = (*pAlleles & AlleleMsk) >> (AlleleIdx * 2);
		switch (Allele) {
		case 0x03:
			if ((LoMinorAllele = (*pAlleles & ~AlleleMsk)) != 0) // can only be one major, but could be a very low cover minor ...
				{
				if(LoMinorAllele == 0x01 || LoMinorAllele == 0x04 || LoMinorAllele == 0x10 || LoMinorAllele == 0x40) // when low cover minor then remove that minor
					{
					*pAlleles = AlleleMsk;
					return(true);
					}

				if (bSetNoAlleles)
					*pAlleles = 0;
				return(false);
				}
			else
				return(true);
		case 0x02:
			if (bNormalise && (*pAlleles & ~AlleleMsk) == 0) // if low coverage, a major allele is represented as 2,0,0,0 whereas if high coverage then that same allele would have been represented as 3,0,0,0
			{
				*pAlleles = AlleleMsk; // normalisation transforms low coverage major alleles into same representation as high coverage
				return(true);
			}
			Num2s++;
			continue;
		case 0x01:
			Num1s++;
			continue;
		}
	}
	if (Num2s > 2 || Num1s > 2 || (Num2s + Num1s > 2))
	{
		if (bSetNoAlleles)
			*pAlleles = 0;
		return(false);
	}
	return(true);
}


int32_t				// error or success (>=0) code 
CDGTvQTLs::ProcessDGTs(char* pszHapGrpDGTFile,				// input, previously generated by 'callhaplotypes', haplotype group DGT file (CSV format) 
					int32_t NumSamplePBAsInputFiles,	     // number of input sample PBAs file specs
					char* pszSamplePBAsInputFiles[],		// names of input sample PBAs files (wildcards allowed)
					char* pszOutFile)						// write results to this output file (CSV format)
{
int Rslt;
gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessDGTs: Initialising for processing of correspondence between DGTs and QTLs");
Rslt = eBSFSuccess;		// assume success!
CSimpleGlob glob(SG_GLOB_FULLSORT);
int32_t Idx;
char* pszInFile;
int32_t ReadsetID = 0;
int32_t NumFiles;
int32_t TotNumFiles = 0;

if (m_pInBuffer == nullptr)
	{
	if ((m_pInBuffer = new uint8_t[cInBuffSize]) == nullptr)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessGrpLociDGTs: Memory allocation of %u bytes for input buffering failed - %s", cInBuffSize, strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
	m_AllocInBuff = cInBuffSize;
	}
m_InNumBuffered = 0;

	// load all samples with out actually allocating memory for each individual chromosome PBAs, but the file offset at which the chromosome PBA starts will be recorded  
for (Idx = 0; Idx < NumSamplePBAsInputFiles; Idx++)
{
	glob.Init();
	if (glob.Add(pszSamplePBAsInputFiles[Idx]) < SG_SUCCESS)
	{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessGrpLociDGTs: Unable to glob '%s' sample PBA files", pszSamplePBAsInputFiles[Idx]);
		Reset();
		return(eBSFerrOpnFile);	// treat as though unable to open file
	}
	if ((NumFiles = glob.FileCount()) <= 0)
	{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessGrpLociDGTs: Unable to locate any sample PBA file matching '%s", pszSamplePBAsInputFiles[Idx]);
		Reset();
		return(eBSFerrFileName);
	}

	TotNumFiles += NumFiles;

	if (TotNumFiles > cMaxPBAFiles)
	{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessGrpLociDGTs: Can accept at most %d sample PBA files for processing, after wildcard file name expansions there are %d requested", cMaxPBAFiles, TotNumFiles);
		Reset();
		return(eBSFerrMem);
	}

	Rslt = eBSFSuccess;
	for (int32_t FileID = 0; Rslt >= eBSFSuccess && FileID < NumFiles; ++FileID)
	{
		pszInFile = glob.File(FileID);
		ReadsetID = LoadPBAFile(pszInFile, 0, true); // loading without allocation for chromosome PBAs
		if (ReadsetID <= 0)
		{
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessGrpLociDGTs: Errors loading file '%s'", pszInFile);
			Reset();
			return(ReadsetID);
		}
		m_Fndrs2Proc[ReadsetID - 1] = 0x01;	// if loaded then assumption is that this founder will be processed
	}
}
m_NumFounders = TotNumFiles;


if (m_hOutFile != -1)
{
	if (m_OutBuffIdx)
	{
		if (!CUtility::RetryWrites(m_hOutFile, m_pszOutBuffer, m_OutBuffIdx))
		{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessGrpLociDGTs: Fatal error in RetryWrites()");
			Reset();
			return(eBSFerrFileAccess);
		}
		m_OutBuffIdx = 0;
	}
	// commit output file
#ifdef _WIN32
	_commit(m_hOutFile);
#else
	fsync(m_hOutFile);
#endif
	close(m_hOutFile);
	m_hOutFile = -1;
}
Reset();
gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessGrpLociDGTs: Haplotype grouping DGT reporting completed");
return(0);
}




tsDGTQTLAlleles*							// returned instance or nullptr if unable to locate
CDGTvQTLs::LocateDGTQTLAlleles(int32_t ChromID,	// instance on this chrom
			int32_t Loci)						// and at this loci
{
uint32_t Nxt;
uint32_t Hash;
tsDGTQTLAlleles* pDGTQTLAlleles;
if(m_pDGTQTLAlleles == nullptr || ChromID < 1 || Loci < 0)
	return(nullptr);
	
Hash = (uint32_t)((ChromID | (Loci << 3)) % cMaxDGTQTLAllelesHashes); // should be reasonably normally distributed
if (m_DGTQTLAllelesHashes[Hash] == 0)
	return(nullptr);
Nxt = m_DGTQTLAllelesHashes[Hash];		// following linked by hash
do {
	pDGTQTLAlleles = &m_pDGTQTLAlleles[Nxt-1];
	if(pDGTQTLAlleles->ChromID == ChromID && pDGTQTLAlleles->Loci == Loci)
		return(pDGTQTLAlleles);
	}
while((Nxt = pDGTQTLAlleles->Nxt) != 0);
return(nullptr);
}

// Expectation is that m_pDGTQTLAlleles has been sorted ChromID.Loci ascending!!!!
tsDGTQTLAlleles*									// returned instance or nullptr if unable to locate
CDGTvQTLs::LocateDGTQTLAlleles(int32_t ChromID)		// returned instance required on this chrom and at the lowest loci (5')
{
tsDGTQTLAlleles* pDGTQTLAlleles;
tsDGTQTLAlleles* pHitDGTQTLAlleles;
int32_t IdxLo;
int32_t MidIdx;
int32_t IdxHi;

if (m_pDGTQTLAlleles == nullptr || m_UsedDGTQTLAlleles == 0 || ChromID < 1)
	return(nullptr);

// linear search if just a few
if(m_UsedDGTQTLAlleles < 50)
	{
	uint32_t Idx;
	pDGTQTLAlleles = m_pDGTQTLAlleles;
	for(Idx = 0; Idx < m_UsedDGTQTLAlleles; Idx++, pDGTQTLAlleles++)
		if (pDGTQTLAlleles->ChromID == ChromID)
			return(pDGTQTLAlleles);
	return(nullptr);
	}
pHitDGTQTLAlleles = nullptr;
IdxLo = 0;
IdxHi = (int32_t)m_UsedDGTQTLAlleles - 1;
do {
	MidIdx = (IdxHi + IdxLo) / 2;
	pDGTQTLAlleles = &m_pDGTQTLAlleles[MidIdx];
	if (pDGTQTLAlleles->ChromID == ChromID)		// matching on chrom but is it the lowest loci?
		{
		if(MidIdx == 0)
			return(pDGTQTLAlleles);
		if (pDGTQTLAlleles[-1].ChromID != ChromID)
			return(pDGTQTLAlleles);
		IdxHi = MidIdx - 1;
		}
	else		// no match on ChromID
		{
		if (pDGTQTLAlleles->ChromID > ChromID)
			IdxHi = MidIdx - 1;
		else
			IdxLo = MidIdx + 1;
		}
	} 
while (IdxHi >= IdxLo);
return(nullptr);
}

int32_t									// eBSFSuccess or error
CDGTvQTLs::AddDGTQTLAlleles(tsDGTQTLAlleles* pInitAlleles, bool bIsQTL) // add instance of tsDGTQTLAlleles
{
uint32_t ToAllocDGTQTLAlleles;
tsDGTQTLAlleles* pDGTQTLAlleles;
uint32_t Hash;
size_t memreq;

if (pInitAlleles == nullptr || pInitAlleles->ChromID < 1 || pInitAlleles->Loci < 0)	// must be initialised!
	return(eBSFerrParams);

AcquireFastSerialise();
if((pDGTQTLAlleles = LocateDGTQTLAlleles(pInitAlleles->ChromID, pInitAlleles->Loci)) != nullptr) // chrom.loci instance already exists?
	{
	if(bIsQTL)
		{
		memcpy(pDGTQTLAlleles->QTLAlleles, pInitAlleles->QTLAlleles,sizeof(pDGTQTLAlleles->QTLAlleles));
		pDGTQTLAlleles->flgQTL = true;
		}
	else
		{
		memcpy(pDGTQTLAlleles->DGTAlleles, pInitAlleles->DGTAlleles, sizeof(pDGTQTLAlleles->DGTAlleles));
		pDGTQTLAlleles->flgDGT = true;
		}
	ReleaseFastSerialise();
	return(eBSFSuccess);
	}

if (m_pDGTQTLAlleles == nullptr)					// may be nullptr first time in
	{
	memset(m_DGTQTLAllelesHashes,0,sizeof(m_DGTQTLAllelesHashes));
	memreq = (size_t)cAllocDGTQTLAlleles * sizeof(tsDGTQTLAlleles);
#ifdef _WIN32
	m_pDGTQTLAlleles = (tsDGTQTLAlleles*)malloc((size_t)memreq);
	if (m_pDGTQTLAlleles == nullptr)
		{
		ReleaseFastSerialise();
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddDGTQTLAlleles: Memory allocation of %zd bytes failed", (int64_t)memreq);
		return(0);
		}
#else
	m_pDGTQTLAlleles = (tsDGTQTLAlleles*)mmap(nullptr, (size_t)memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	if (m_pDGTQTLAlleles == MAP_FAILED)
		{
		m_pDGTQTLAlleles = nullptr;
		ReleaseFastSerialise();
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddDGTQTLAlleles: Memory allocation of %zd bytes through mmap()  failed", (int64_t)memreq, strerror(errno));
		return(0);
		}
#endif
	m_AllocDGTQTLAllelesMem = memreq;
	m_AllocDGTQTLAlleles = cAllocDGTQTLAlleles;
	m_UsedDGTQTLAlleles = 0;
	}
else
		// needing to allocate more memory?
	if ((m_UsedDGTQTLAlleles) >= m_AllocDGTQTLAlleles)
		{
		ToAllocDGTQTLAlleles = m_UsedDGTQTLAlleles + cAllocDGTQTLAlleles;
		size_t memreq = (size_t)ToAllocDGTQTLAlleles * sizeof(tsDGTQTLAlleles);
#ifdef _WIN32
		pDGTQTLAlleles = (tsDGTQTLAlleles*)realloc(m_pDGTQTLAlleles, memreq);
		if (pDGTQTLAlleles == nullptr)
			{
#else
		pDGTQTLAlleles = (tsDGTQTLAlleles*)mremap(m_pDGTQTLAlleles, m_AllocDGTQTLAllelesMem, memreq, MREMAP_MAYMOVE);
		if (pDGTQTLAlleles == MAP_FAILED)
			{
#endif
			ReleaseFastSerialise();
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddDGTQTLAlleles: Memory reallocation to %zd bytes failed - %s", (int64_t)memreq, strerror(errno));
			return(0);
			}
		m_pDGTQTLAlleles = pDGTQTLAlleles;
		m_AllocDGTQTLAllelesMem = memreq;
		m_AllocDGTQTLAlleles = ToAllocDGTQTLAlleles;
		}
pDGTQTLAlleles = &m_pDGTQTLAlleles[m_UsedDGTQTLAlleles++];
*pDGTQTLAlleles = *pInitAlleles;
if(bIsQTL)
	pDGTQTLAlleles->flgQTL = true;
else
	pDGTQTLAlleles->flgDGT = true;
Hash = (uint32_t)((pInitAlleles->ChromID | (pInitAlleles->Loci << 3)) % cMaxDGTQTLAllelesHashes); // should be reasonably normally distributed
if(m_DGTQTLAllelesHashes[Hash] != 0)
	pDGTQTLAlleles->Nxt = m_DGTQTLAllelesHashes[Hash];
else
	pDGTQTLAlleles->Nxt = 0;
pDGTQTLAlleles->Idx = m_UsedDGTQTLAlleles-1;
m_DGTQTLAllelesHashes[Hash] = m_UsedDGTQTLAlleles;

ReleaseFastSerialise();
return(m_UsedDGTQTLAlleles);
}


int 
CDGTvQTLs::LoadDGTs(char* pszDGTsFile)		// CSV DGTs file containing chrom.loci and alleles for each group tag
{
int Rslt = 0;
int32_t CurLineNumber;
int32_t NumFields;
char* pszChrom;
int32_t CurLoci;
int32_t CurChromID;
int32_t AlleleAGrp;
int32_t AlleleCGrp;
int32_t AlleleGGrp;
int32_t AlleleTGrp;
tsDGTQTLAlleles DGTQTLAlleles;
if (m_pInCSVFile != nullptr) // should be, but better to be sure!
	{
	delete m_pInCSVFile;
	m_pInCSVFile = nullptr;
	}

if ((m_pInCSVFile = new CCSVFile) == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadDGTs: Unable to instantiate CCSVfile");
	Reset();
	return(eBSFerrObj);
	}

m_pInCSVFile->SetMaxFields(52); // only processing 1st 51 fields, allowing 1 extra so can detect if more thasn 51

if ((Rslt = m_pInCSVFile->Open(pszDGTsFile)) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadDGTs: Unable to open file: '%s'", pszDGTsFile);
	Reset();
	return(Rslt);
	}

m_NumDGTsLoaded = 0;
CurLineNumber = 0;
NumFields = 0;
while ((Rslt = m_pInCSVFile->NextLine()) > 0)		// onto next line containing fields
	{
	CurLineNumber++;
	if ((NumFields = m_pInCSVFile->GetCurFields()) < 51)	// must contain at least 51 fields
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadDGTs: Allele association score file '%s' expected to contain a minimum of 51 fields, it contains %d at line %d", pszDGTsFile, NumFields, CurLineNumber);
		Reset();
		return(eBSFerrParse);
		}

	// 1st row may be a header row, if so then skip this row
	if (CurLineNumber == 1 && m_pInCSVFile->IsLikelyHeaderLine())
		continue;

	Rslt = m_pInCSVFile->GetText(5, &pszChrom);
	if (Rslt < 0 || pszChrom == nullptr || pszChrom[0] == '\0')
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadDGTs: Unable to parse chromosome name at line %d in file '%s'", CurLineNumber, pszDGTsFile);
		Reset();
		return(eBSFerrParse);
		}

	if (!AcceptThisChromName(pszChrom))
		continue;
	CurChromID = LocateChrom(pszChrom);
	Rslt = m_pInCSVFile->GetInt(6, &CurLoci);    // note 0 based loci
	Rslt = m_pInCSVFile->GetInt(7, &AlleleAGrp);
	Rslt = m_pInCSVFile->GetInt(8, &AlleleCGrp);
	Rslt = m_pInCSVFile->GetInt(9, &AlleleGGrp);
	Rslt = m_pInCSVFile->GetInt(10, &AlleleTGrp);
	memset(&DGTQTLAlleles,0,sizeof(DGTQTLAlleles));
	DGTQTLAlleles.ChromID = CurChromID;
	DGTQTLAlleles.Loci = CurLoci;
	DGTQTLAlleles.DGTAlleles[0] = AlleleAGrp;
	DGTQTLAlleles.DGTAlleles[1] = AlleleCGrp;
	DGTQTLAlleles.DGTAlleles[2] = AlleleGGrp;
	DGTQTLAlleles.DGTAlleles[3] = AlleleTGrp;
	AddDGTQTLAlleles(&DGTQTLAlleles,false);
	m_NumDGTsLoaded++;
	}

if (m_pInCSVFile != nullptr)
	{
	delete m_pInCSVFile;
	m_pInCSVFile = nullptr;
	}
return(m_NumDGTsLoaded);
}

int
CDGTvQTLs::LoadQTLs(char* pszQTLsFile)		// CSV QTLs file containing chrom.loci and alleles for each QTL
{
int Rslt = 0;
int32_t CurLineNumber;
int32_t NumFields;
char* pszChrom;
char* pszAllele;
eSeqBase RefAllele;
eSeqBase AltAllele;
int32_t CurLoci;
int32_t CurChromID;
tsDGTQTLAlleles DGTQTLAlleles;

if (m_pInCSVFile != nullptr) // should be, but better to be sure!
	{
	delete m_pInCSVFile;
	m_pInCSVFile = nullptr;
	}

if ((m_pInCSVFile = new CCSVFile) == nullptr)
{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadQTLs: Unable to instantiate CCSVfile");
	Reset();
	return(eBSFerrObj);
}

m_pInCSVFile->SetMaxFields(9); // only processing 1st 8 fields, allowing 1 extra so can detect if more thasn 8

if ((Rslt = m_pInCSVFile->Open(pszQTLsFile)) != eBSFSuccess)
{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadQTLs: Unable to open file: '%s'", pszQTLsFile);
	Reset();
	return(Rslt);
}

m_NumQTLsLoaded = 0;
CurLineNumber = 0;
NumFields = 0;
while ((Rslt = m_pInCSVFile->NextLine()) > 0)		// onto next line containing fields
{
	CurLineNumber++;
	if ((NumFields = m_pInCSVFile->GetCurFields()) < 8)	// must contain at least 8 fields
	{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadQTLs: QTLs file '%s' expected to contain a minimum of 8 fields, it contains %d at line %d", pszQTLsFile, NumFields, CurLineNumber);
		Reset();
		return(eBSFerrParse);
	}

	// 1st row may be a header row, if so then skip this row
	if (CurLineNumber == 1 && m_pInCSVFile->IsLikelyHeaderLine())
		continue;

	Rslt = m_pInCSVFile->GetText(2, &pszChrom);
	if (Rslt < 0 || pszChrom == nullptr || pszChrom[0] == '\0')
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadQTLs: Unable to parse chromosome name at line %d in file '%s'", CurLineNumber, pszQTLsFile);
		Reset();
		return(eBSFerrParse);
		}



	if (!AcceptThisChromName(pszChrom))
		continue;
	CurChromID = LocateChrom(pszChrom);
	Rslt = m_pInCSVFile->GetInt(3, &CurLoci);

	Rslt = m_pInCSVFile->GetText(4, &pszAllele);
	switch(*pszAllele) {
		case 'a': case 'A':
			RefAllele = eBaseA;
			break;
		case 'c': case 'C':
			RefAllele = eBaseC;
			break;
		case 'g': case 'G':
			RefAllele = eBaseG;
			break;
		case 't': case 'T':
			RefAllele = eBaseT;
			break;
		default:
			RefAllele = eBaseN;
		}
	Rslt = m_pInCSVFile->GetText(5, &pszAllele);
	switch (*pszAllele) {
		case 'a': case 'A':
			AltAllele = eBaseA;
			break;
		case 'c': case 'C':
			AltAllele = eBaseC;
			break;
		case 'g': case 'G':
			AltAllele = eBaseG;
			break;
		case 't': case 'T':
			AltAllele = eBaseT;
			break;
		default:
			AltAllele = eBaseN;
		}
	memset(&DGTQTLAlleles, 0, sizeof(DGTQTLAlleles));
	DGTQTLAlleles.ChromID = CurChromID;
	DGTQTLAlleles.Loci = CurLoci - 1;  // DGTs are using 0 based loci, QTLs are 1 based
	DGTQTLAlleles.QTLAlleles[0] = RefAllele;
	DGTQTLAlleles.QTLAlleles[1] = AltAllele;
	AddDGTQTLAlleles(&DGTQTLAlleles, true);
	m_NumQTLsLoaded++;
	}

if (m_pInCSVFile != nullptr)
	{
	delete m_pInCSVFile;
	m_pInCSVFile = nullptr;
	}
return(m_NumQTLsLoaded);
}

int				// returns 0 if all threads started and completed processing, 1 if all threads started but some yet to complete within WaitSecs, 2 if not all threads have started within WaitSecs
CDGTvQTLs::WaitWorkerThreadStatus(int WaitSecs)	// // wait at most this many seconds for all threads to start and complete processing before reporting on worker thread status
{
	int64_t WaitMS;
	if (WaitSecs < 1)
		WaitSecs = 1;
	WaitMS = (int64_t)WaitSecs * 1000;

	uint32_t CompletedInstances;
	uint32_t NumWorkerInsts;
	uint32_t ExpNumWorkerInsts;
	do {
		CUtility::SleepMillisecs(500);		// poll every 0.5 secs, restricting poll frequency
		AcquireSerialise();
		CompletedInstances = m_CompletedWorkerInsts;
		NumWorkerInsts = m_NumWorkerInsts;
		ExpNumWorkerInsts = m_ExpNumWorkerInsts;
		ReleaseSerialise();
		WaitMS -= 500;
	} while (ExpNumWorkerInsts > 0 &&		// surely must be some workers expected!
		CompletedInstances < ExpNumWorkerInsts && // have all workers completed 
		WaitMS > 0);					// within wait period (WaitSecs)

	return((NumWorkerInsts == ExpNumWorkerInsts && CompletedInstances == ExpNumWorkerInsts) ? 0 : (NumWorkerInsts == ExpNumWorkerInsts ? 1 : 2));
}


#ifdef _WIN32
unsigned __stdcall WorkerLoadChrPBAInstance(void* pThreadPars)
#else
void* WorkerLoadChrPBAInstance(void* pThreadPars)
#endif
{
	int Rslt;
	tsCHWorkerLoadChromPBAsInstance* pPars = (tsCHWorkerLoadChromPBAsInstance*)pThreadPars;			// makes it easier not having to deal with casts!
	CDGTvQTLs* pWorkerInstance = (CDGTvQTLs*)pPars->pThis;

	Rslt = pWorkerInstance->ProcWorkerLoadChromPBAThread(pPars);
	pPars->Rslt = Rslt;
#ifdef _WIN32
	_endthreadex(0);
	return(eBSFSuccess);
#else
	pthread_exit(&pPars->Rslt);
#endif
}

// initialise and start pool of worker threads to concurrently load chromosome PBAs
int				// returns 0 if all threads have started processing, 1 if all threads started and some have completed processing, 2 if not all have started
CDGTvQTLs::StartWorkerLoadChromPBAThreads(int32_t NumThreads,		// there are this many threads in pool
	int32_t StartSampleID,				// processing to start from this sample identifer
	int32_t EndSampleID,				// ending with this sample identifier inclusive
	int32_t ChromID,					// loading PBAs for this chromosome
	bool bNormAlleles)					// normalise alleles such that individual alleles can be compared without regard to the proportional coverage (0x22 -> 0x33 as an example)
{
	int Rslt = eBSFSuccess;
	int32_t FromSampleID;
	int32_t NumSamples;
	int32_t SamplesThisThread;
	int32_t ThreadIdx;
	tsCHWorkerLoadChromPBAsInstance* pThreadPar;

#ifndef _WIN32
	// increase the default thread stack to at least cWorkThreadStackSize
	size_t defaultStackSize;
	pthread_attr_t threadattr;
	pthread_attr_init(&threadattr);
	pthread_attr_getstacksize(&threadattr, &defaultStackSize);
	if (defaultStackSize < (size_t)cWorkThreadStackSize)
	{
		if ((Rslt = pthread_attr_setstacksize(&threadattr, (size_t)cWorkThreadStackSize)) != 0)
		{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "StartWorkerLoadChromPBAThreads: pthread_attr_setstacksize(%d) failed, default was %zd", cWorkThreadStackSize, defaultStackSize);
			return(eBSFerrInternal);
		}
	}
#endif


	NumSamples = 1 + EndSampleID - StartSampleID;
	NumThreads = min(NumSamples, NumThreads);
	FromSampleID = StartSampleID;
	m_NumWorkerInsts = 0;
	m_ExpNumWorkerInsts = NumThreads;
	m_CompletedWorkerInsts = 0;
	m_ReqTerminate = 0;
	pThreadPar = m_WorkerLoadChromPBAInstances;
	for (ThreadIdx = 1; ThreadIdx <= NumThreads; ThreadIdx++, pThreadPar++)
	{
		memset(pThreadPar, 0, sizeof(tsCHWorkerLoadChromPBAsInstance));
		pThreadPar->bNormAlleles = bNormAlleles;
		pThreadPar->ChromID = ChromID;
		SamplesThisThread = NumSamples / (1 + NumThreads - ThreadIdx);
		NumSamples -= SamplesThisThread;
		pThreadPar->StartSampleID = FromSampleID;
		pThreadPar->EndSampleID = SamplesThisThread - 1 + FromSampleID;
		FromSampleID = pThreadPar->EndSampleID + 1;

		pThreadPar->ThreadIdx = ThreadIdx;
		pThreadPar->pThis = this;
#ifdef _WIN32
		pThreadPar->threadHandle = (HANDLE)_beginthreadex(nullptr, cWorkThreadStackSize, WorkerLoadChrPBAInstance, pThreadPar, 0, &pThreadPar->threadID);
#else
		pThreadPar->threadRslt = pthread_create(&pThreadPar->threadID, &threadattr, WorkerLoadChrPBAInstance, pThreadPar);
#endif
	}

	// allow threads time to all startup
	Rslt = WaitWorkerThreadStatus(cMaxWaitThreadsStartup);	// waiting until all threads have started (some may have completed!)
	if (Rslt == 2)
	{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "StartWorkerThreads: All %d expected threads did not register as started within allowed %d seconds", NumThreads, cMaxWaitThreadsStartup);
		Rslt = TerminateWorkerThreads();
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "StartWorkerThreads: TerminateWorkThreads(60) returned %d", Rslt);
		return(eBSFerrInternal);
	}
	return(Rslt);
}

// stop all threads in worker pool
int	// returns number of threads forced to terminate
CDGTvQTLs::TerminateWorkerThreads(int WaitSecs)				// allow at most this many seconds before force terminating worker pool threads
{
	int NumForceTerminated;
	uint32_t Idx;
	uint32_t NumWorkerInsts;
	uint32_t ExpNumWorkerInsts;
	uint32_t CompletedWorkerInsts;
	tsCHWorkerInstance* pThreadPar;
	time_t Then;
	time_t Now;

	AcquireSerialise();
	NumWorkerInsts = m_NumWorkerInsts;
	ExpNumWorkerInsts = m_ExpNumWorkerInsts;
	CompletedWorkerInsts = m_CompletedWorkerInsts;
	ReleaseSerialise();
	if (ExpNumWorkerInsts == 0)
	{
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "TerminateWorkerThreads: No worker threads to terminate");
		pThreadPar = m_WorkerInstances;
		AcquireSerialise();
		for (Idx = 0; Idx < ExpNumWorkerInsts; Idx++, pThreadPar += 1)
			memset(pThreadPar, 0, sizeof(tsCHWorkerInstance));
		m_NumWorkerInsts = 0;
		m_ExpNumWorkerInsts = 0;
		m_CompletedWorkerInsts = 0;
		m_ReqTerminate = 0;
		ReleaseSerialise();
		return(0);
	}

	// request all worker threads to self terminate
	AcquireSerialise();
	m_ReqTerminate = 1;
	ReleaseSerialise();
	if (CompletedWorkerInsts < ExpNumWorkerInsts)
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "TerminateWorkerThreads: Requesting %u worker threads to terminate", ExpNumWorkerInsts - CompletedWorkerInsts);
	Then = time(nullptr) + WaitSecs;
	NumForceTerminated = 0;
	pThreadPar = m_WorkerInstances;
	for (Idx = 0; Idx < ExpNumWorkerInsts; Idx++, pThreadPar += 1)
	{
		Now = time(nullptr);
		if (Now >= Then)
			Now = 1;
		else
			Now = Then - Now;

#ifdef WIN32
		if (pThreadPar->threadHandle != nullptr)
		{
			if (WAIT_TIMEOUT == WaitForSingleObject(pThreadPar->threadHandle, (uint32_t)Now * 1000))
			{
				NumForceTerminated += 1;
				gDiagnostics.DiagOut(eDLInfo, gszProcName, "TerminateWorkerThreads: Force terminating thread %u, WaitForSingleObject() returned WAIT_TIMEOUT", pThreadPar->ThreadIdx);
				TerminateThread(pThreadPar->threadHandle, 0);
			}
			pThreadPar->threadHandle = nullptr;
		}
#else
		if (pThreadPar->threadID != 0)
		{
			struct timespec ts;
			int JoinRlt;
			void* pExitRslt;
			clock_gettime(CLOCK_REALTIME, &ts);
			ts.tv_sec += Now;
			if ((JoinRlt = pthread_timedjoin_np(pThreadPar->threadID, &pExitRslt, &ts)) != 0)
			{
				gDiagnostics.DiagOut(eDLInfo, gszProcName, "TerminateWorkerThreads: Force terminating thread %u, pthread_timedjoin_np() returned %d", pThreadPar->ThreadIdx, JoinRlt);
				NumForceTerminated += 1;
				pthread_cancel(pThreadPar->threadID);
				pthread_join(pThreadPar->threadID, nullptr);
			}
			pThreadPar->threadID = 0;
		}
#endif
	}

	pThreadPar = m_WorkerInstances;
	AcquireSerialise();
	for (Idx = 0; Idx < ExpNumWorkerInsts; Idx++, pThreadPar += 1)
		memset(pThreadPar, 0, sizeof(tsCHWorkerInstance));
	m_NumWorkerInsts = 0;
	m_ExpNumWorkerInsts = 0;
	m_CompletedWorkerInsts = 0;
	m_ReqTerminate = 0;
	ReleaseSerialise();
	return(NumForceTerminated);
}

// stop all threads in worker load chromosome PBAs thread pool
int
CDGTvQTLs::TerminateLoadChromPBAsThreads(int WaitSecs)				// allow at most this many seconds before force terminating worker pool threads
{
	int NumForceTerminated;
	uint32_t Idx;
	uint32_t NumWorkerInsts;
	uint32_t ExpNumWorkerInsts;
	uint32_t CompletedWorkerInsts;
	tsCHWorkerLoadChromPBAsInstance* pThreadPar;
	time_t Then;
	time_t Now;

	AcquireSerialise();
	NumWorkerInsts = m_NumWorkerInsts;
	ExpNumWorkerInsts = m_ExpNumWorkerInsts;
	CompletedWorkerInsts = m_CompletedWorkerInsts;
	ReleaseSerialise();
	if (ExpNumWorkerInsts == 0)
	{
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "TerminateLoadChromPBAsThreads: No worker threads to terminate");
		pThreadPar = m_WorkerLoadChromPBAInstances;
		AcquireSerialise();
		for (Idx = 0; Idx < ExpNumWorkerInsts; Idx++, pThreadPar += 1)
			memset(pThreadPar, 0, sizeof(tsCHWorkerLoadChromPBAsInstance));
		m_NumWorkerInsts = 0;
		m_ExpNumWorkerInsts = 0;
		m_CompletedWorkerInsts = 0;
		m_ReqTerminate = 0;
		ReleaseSerialise();
		return(0);
	}

	// request all worker threads to self terminate
	AcquireSerialise();
	m_ReqTerminate = 1;
	ReleaseSerialise();
	if (NumWorkerInsts != CompletedWorkerInsts)
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "TerminateLoadChromPBAsThreads: Requesting %u worker threads to terminate", NumWorkerInsts);
	Then = time(nullptr) + WaitSecs;
	NumForceTerminated = 0;
	pThreadPar = m_WorkerLoadChromPBAInstances;
	for (Idx = 0; Idx < ExpNumWorkerInsts; Idx++, pThreadPar += 1)
	{
		Now = time(nullptr);
		if (Now >= Then)
			Now = 1;
		else
			Now = Then - Now;

#ifdef WIN32
		if (pThreadPar->threadHandle != nullptr)
		{
			if (WAIT_TIMEOUT == WaitForSingleObject(pThreadPar->threadHandle, (uint32_t)Now * 1000))
			{
				NumForceTerminated += 1;
				TerminateThread(pThreadPar->threadHandle, 0);
			}
			pThreadPar->threadHandle = nullptr;
		}
#else
		if (pThreadPar->threadID != 0)
		{
			struct timespec ts;
			int JoinRlt;
			void* pExitRslt;
			clock_gettime(CLOCK_REALTIME, &ts);
			ts.tv_sec += Now;
			if ((JoinRlt = pthread_timedjoin_np(pThreadPar->threadID, &pExitRslt, &ts)) != 0)
			{
				gDiagnostics.DiagOut(eDLInfo, gszProcName, "TerminateLoadChromPBAsThreads: Force terminating thread %u, pthread_timedjoin_np() returned %d", pThreadPar->ThreadIdx, JoinRlt);
				NumForceTerminated += 1;
				pthread_cancel(pThreadPar->threadID);
				pthread_join(pThreadPar->threadID, nullptr);
			}
			pThreadPar->threadID = 0;
		}
#endif
	}

	pThreadPar = m_WorkerLoadChromPBAInstances;
	AcquireSerialise();
	for (Idx = 0; Idx < ExpNumWorkerInsts; Idx++, pThreadPar += 1)
		memset(pThreadPar, 0, sizeof(tsCHWorkerLoadChromPBAsInstance));
	m_NumWorkerInsts = 0;
	m_ExpNumWorkerInsts = 0;
	m_CompletedWorkerInsts = 0;
	m_ReqTerminate = 0;
	ReleaseSerialise();
	return(NumForceTerminated);
}

// Using threads to concurrently load chromosome PBAs
// ChromID is specified and range of sample identifiers
int
CDGTvQTLs::ProcWorkerLoadChromPBAThread(tsCHWorkerLoadChromPBAsInstance* pThreadPar)	// worker thread parameters
{
	int Rslt;
	uint8_t PBA;
	uint8_t* pPBAs;
	int32_t Ofs;
	int32_t CurSampleID;
	uint32_t ReqTerminate;
	// one more thread instance has started
	AcquireSerialise();
	m_NumWorkerInsts++;
	ReleaseSerialise();

	for (CurSampleID = pThreadPar->StartSampleID; CurSampleID <= pThreadPar->EndSampleID; CurSampleID++)
	{
		// check if requested to terminate
		AcquireSerialise();
		ReqTerminate = m_ReqTerminate;
		ReleaseSerialise();
		if (ReqTerminate)
			break;

		if ((pPBAs = LoadSampleChromPBAs(CurSampleID, pThreadPar->ChromID)) == nullptr)
			break;

		if (pThreadPar->bNormAlleles)
		{
			tsCHChromMetadata* pChromMetadata = LocateChromMetadataFor(CurSampleID, pThreadPar->ChromID);
			for (Ofs = 0; Ofs < pChromMetadata->ChromLen; Ofs++, pPBAs++)
			{
				if ((PBA = *pPBAs) == 0)
					continue;

				if (PBA & 0xc0)
					PBA |= 0xc0;
				if (PBA & 0x30)
					PBA |= 0x30;
				if (PBA & 0x0c)
					PBA |= 0x0c;
				if (PBA & 0x03)
					PBA |= 0x03;
				*pPBAs = PBA;
			}
		}
	}
	if (CurSampleID <= pThreadPar->EndSampleID)
		Rslt = -1;
	else
		Rslt = eBSFSuccess;
	AcquireSerialise();
	m_CompletedWorkerInsts++;
	ReleaseSerialise();
	return(Rslt);
}

int
CDGTvQTLs::ProcessDGTQTLsPBAs(eModeDGTA PMode)	// actual processing of DGTs/QTLs against reference and sample PBAs
{
int32_t ChromID;
int32_t SampleID;
int32_t DelSampleID;
uint8_t** ppPBAs = nullptr;

m_NumQTLInstances = 0;					// number of QTLs characterised
m_NumDGTInstances = 0;
m_NumDGTQTLbothInstances = 0;
m_NumQTLRefMismatched = 0;
m_NumSamplesLoCov = 0;
m_NumSamplesHeteroz = 0;
m_NumSamplesHomoz = 0;
m_NumSamplesRefMismatched = 0;
m_NumSamplesMonoAllelic = 0;
m_NumSamplesPolyAllelic = 0;

ppPBAs = new uint8_t * [cMaxPBAFiles + 1];
memset(ppPBAs, 0, sizeof(uint8_t*) * (cMaxPBAFiles + 1));

for (ChromID = 1; ChromID <= m_NumChromNames; ChromID++)
	{
	uint32_t ChromSize;
	char* pszChrom = LocateChrom(ChromID);

	gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessDGTQTLsPBAs: Beginning to load PBAs for %s chromosome for DGT processing .... ", pszChrom);
	int WorkerThreadStatus = StartWorkerLoadChromPBAThreads(m_NumThreads, 1, m_NumFounders, ChromID);
	if (WorkerThreadStatus < 0)
		{
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessGrpLociDGTs: Errors loading PBAs for %s chromosome for DGT processing .... ", pszChrom);
		delete[]ppPBAs;
		return(WorkerThreadStatus);
		}
	while (WorkerThreadStatus > 0)
		{
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessGrpLociDGTs: Continuing to load chromosome PBAs for %s .... ", pszChrom);
		WorkerThreadStatus = WaitWorkerThreadStatus(60);
		}
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessGrpLociDGTs: Completed loading chromosome PBAs for %s for DGT processing", pszChrom);
	TerminateLoadChromPBAsThreads(60);

	for (SampleID = 1; SampleID <= m_NumFounders; SampleID++)
		{
		tsCHChromMetadata* pChromMetadata;
		// returned pointer to chromosome metadata
		if ((pChromMetadata = LocateChromMetadataFor(SampleID, ChromID)) == nullptr)
			{
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessGrpLociDGTs: No metadata for chromosome '%s' in at least one founder '%s', skipping this chromosome", pszChrom, LocateReadset(SampleID));
			break;;
			}
		ChromSize = pChromMetadata->ChromLen;
		if ((ppPBAs[SampleID - 1] = pChromMetadata->pPBAs) == nullptr)
			{
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessGrpLociDGTs: No PBA for chromosome '%s' in at least one founder '%s', skipping this chromosome", pszChrom, LocateReadset(SampleID));
			break;
			}
		}

	if (SampleID <= m_NumFounders)
		{
		for (DelSampleID = 1; DelSampleID < SampleID; DelSampleID++)
			DeleteSampleChromPBAs(DelSampleID, ChromID);
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessGrpLociDGTs: Unloaded PBAs for chromosome %s completed, ready for next iteration of chromosome DGT processing", pszChrom);
		continue; // slough this chromosome and try next chromosome
		}
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessDGTQTLsPBAs: Loading all PBAs for chromosome %s completed, now processing ...", pszChrom);

	// iterate over the DGT/QTLs and process against the PBAs
	tsDGTQTLAlleles *pDGTQTLAlleles;

	int32_t ProcessedThisChrom = 0;
	int32_t Characterised = 0;
	if((pDGTQTLAlleles = LocateDGTQTLAlleles(ChromID))==nullptr)
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessDGTQTLsPBAs: No DGTQTLAlleles for chromosome %s located ...", pszChrom);
	else
		{
		for(ProcessedThisChrom = 0; pDGTQTLAlleles->ChromID == ChromID && pDGTQTLAlleles->Idx < m_UsedDGTQTLAlleles; ProcessedThisChrom++,pDGTQTLAlleles++)
			{
			// hand over this instance to generative function for actual processing
			if(AnalyseInstance(PMode,pDGTQTLAlleles, m_NumFounders, ppPBAs)==0)
				Characterised++;
			}
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessDGTQTLsPBAs: Processing of %d DGTQTLAlleles of which %d accepted for characterisation on chromosome %s completed", ProcessedThisChrom, Characterised, pszChrom);
		}

	for (DelSampleID = 1; DelSampleID <= SampleID; DelSampleID++)
		DeleteSampleChromPBAs(DelSampleID, ChromID);
	}
if (ppPBAs != nullptr)
	delete[]ppPBAs;


gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessDGTQTLsPBAs:  Number of DGT instance loci characterised: %d", m_NumDGTInstances);
gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessDGTQTLsPBAs:  Number of QTL instance loci characterised: %d", m_NumQTLInstances);
gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessDGTQTLsPBAs:  Number of DGT and QTL shared instance loci characterised: %d", m_NumDGTQTLbothInstances);
gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessDGTQTLsPBAs:  Number of QTL instance loci claimed reference base not matching reference assembly base: %d", m_NumQTLRefMismatched);
gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessDGTQTLsPBAs:  Number of DGT/QTL instance loci with low sample stacking coverage: %d", m_NumSamplesLoCov);
gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessDGTQTLsPBAs:  Number of DGT/QTL instance loci with heterozygous sample stacking : %d", m_NumSamplesHeteroz);
gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessDGTQTLsPBAs:  Number of DGT/QTL instance loci with homozygous sample stacking: %d", m_NumSamplesHomoz);
gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessDGTQTLsPBAs:  Number of DGT/QTL instance loci with Grp1 homozygous sample stacking and alleles not matching: reference assembly base: %d", m_NumSamplesRefMismatched);
gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessDGTQTLsPBAs:  Number of DGT/QTL instance loci with Grp1 homozygous sample stacking which are monoallelic: %d", m_NumSamplesMonoAllelic);
gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessDGTQTLsPBAs:  Number of DGT/QTL instance loci with Grp1 homozygous sample stacking which are polyallelic: %d", m_NumSamplesPolyAllelic);

return(0);
}

const char *
CDGTvQTLs::Diplotype2Txt(uint8_t Alleles)
{
	switch (Alleles) {
		case 0xc0:
			return((const char *)"A/A");
		case 0x30:
			return((const char*)"C/C");
		case 0x0c:
			return((const char*)"G/G");
		case 0x03:
			return((const char*)"T/T");
		case 0xf0:
			return((const char*)"A/C");
		case 0xcc:
			return((const char*)"A/G");
		case 0xc3:
			return((const char*)"A/T");
		case 0x3c:
			return((const char*)"C/G");
		case 0x33:
			return((const char*)"C/T");
		case 0x0f:
			return((const char*)"G/T");
		case 0x00:
			break;
		default:
			return((const char*)"-/-");
		}
return((const char *)"N/N");
}

int				// 0 if accepted, 1 if unprocessed
CDGTvQTLs::AnalyseInstance(eModeDGTA PMode,				// processing mode
				tsDGTQTLAlleles* pDGTQTLAlleles,		// processing this instance
				int32_t NumPBAs,						// against this number of PBAs 
				uint8_t** ppPBAs)						// ptrs to loaded PBAs
{
int32_t PBADist[256];
uint8_t *pPBA;
uint8_t PBA;
uint8_t RefPBA;
uint8_t HiFreqPBA;
uint8_t NxtHiFreqPBA;
int32_t SampleIdx;
double Coverage;
double Grp1Prop;
double Grp2Prop;
double GrpNProp;

switch(PMode) {
	case eDGTADefault:							// instance must contain QTL
		if(!pDGTQTLAlleles->flgQTL)
			return(1);
		break;
	case eDGTADGT:								// reporting DGT or QTL loci instances
		break;
	default:
		return(1);					// treating as unprocessed
}
// tracking down as to why there are discrepancies between the DGTs and QTLs
pDGTQTLAlleles->flgQTLRefMismatch = false;
pDGTQTLAlleles->flgSamplesRefMismatch = false;
pDGTQTLAlleles->flgSamplesHomoz = false;
pDGTQTLAlleles->flgSamplesPolyAllelic = false;
pDGTQTLAlleles->flgSamplesMonoAllelic = false;

// what are the reference assembly alleles at the current loci?
pPBA = ppPBAs[0];
RefPBA = pPBA[pDGTQTLAlleles->Loci];

if (pDGTQTLAlleles->flgDGT && pDGTQTLAlleles->flgQTL)				// contains both a DGT and QTL?
	m_NumDGTQTLbothInstances++;

if (pDGTQTLAlleles->flgQTL)				// contains a QTL?
	{
	m_NumQTLInstances++;

		// Is assembly reference base same as claimed by QTL?
	pPBA = ppPBAs[0];
	if (RefPBA != (0xc0 >> (pDGTQTLAlleles->QTLAlleles[0]*2)))
		{
		pDGTQTLAlleles->flgQTLRefMismatch = true;
		m_NumQTLRefMismatched++;
		}
	}
if (pDGTQTLAlleles->flgDGT)				// contains a DGT?
	m_NumDGTInstances++;

// determine proportions of allelic combinations over all samples
memset(PBADist, 0, sizeof(PBADist));
HiFreqPBA = 0;
for (SampleIdx = 1; SampleIdx < NumPBAs; SampleIdx++)
	{
	pPBA = ppPBAs[SampleIdx];		// reference assembly was at ppPBAs[0]
	PBA = pPBA[pDGTQTLAlleles->Loci];
	if(PBA & 0xc0)		// not interested in major/minor, just allele presence
		PBA |= 0xc0;
	if (PBA & 0x30)
		PBA |= 0x30;
	if (PBA & 0x0c)
		PBA |= 0x0c;
	if (PBA & 0x03)
		PBA |= 0x03;

	PBADist[PBA] += 1;
	if (PBADist[PBA] > PBADist[HiFreqPBA]) // new highest frequency?
		HiFreqPBA = PBA;
	}


Coverage = 1.0 - ((double)PBADist[0] / (NumPBAs - 1));		// PBADist[0] contains counts of samples with no coverage
Grp1Prop = Grp2Prop = GrpNProp = 0;

// if sufficent coverage and highest frequency known then find 2nd highest
if(Coverage >= m_MinCoverage)
	{
	NxtHiFreqPBA = 1;
	for (PBA = 1; PBA < 255; PBA++)
		{
		if(PBA == HiFreqPBA)
			continue;
		if (PBADist[PBA] > PBADist[NxtHiFreqPBA]) // new highest frequency?
			NxtHiFreqPBA = PBA;
		}
	Diplotype2Txt(NxtHiFreqPBA);
	Grp1Prop = (double)PBADist[HiFreqPBA] / (NumPBAs - PBADist[0] - 1);
	Grp2Prop = (double)PBADist[NxtHiFreqPBA] / (NumPBAs - PBADist[0] - 1);
	GrpNProp = max(1.0 - (Grp1Prop + Grp2Prop),0.0);
	}
else
	{
	NxtHiFreqPBA = 0;
	HiFreqPBA = 0;
	}

if (Coverage < m_MinCoverage)			// was defaulted as 0.8, if coverage < this threshold then don't further characterise
	{
	pDGTQTLAlleles->flgSamplesLoCov = true;
	m_NumSamplesLoCov++;
	}
else	// else have sufficient coverage to further characterise
	{   // note that the non-coverage counts are excluded!
	switch (HiFreqPBA) {
		case 0xc0: case 0x30: case 0x0c: case 0x03:
			pDGTQTLAlleles->flgSamplesMonoAllelic = true;
			m_NumSamplesMonoAllelic++;
			break;
		default:
			pDGTQTLAlleles->flgSamplesPolyAllelic = true;
			m_NumSamplesPolyAllelic++;
			break;
		}

	if (HiFreqPBA != RefPBA)
		{
		pDGTQTLAlleles->flgSamplesRefMismatch = true;
		m_NumSamplesRefMismatched++;
		}

	if (Grp1Prop >= m_HomozPropThres)  // if Grp1Prop is more than this proportion then classify as homozygous at that loci
		{
		pDGTQTLAlleles->flgSamplesHomoz = true;
		m_NumSamplesHomoz++;
		}
	else
		{
		pDGTQTLAlleles->flgSamplesHeteroz = true;
		m_NumSamplesHeteroz++;
		}
	}


if (m_hOutFile != -1 && m_pszOutBuffer != nullptr)
	{
	int LociType;
	if(pDGTQTLAlleles->flgDGT && pDGTQTLAlleles->flgQTL)
		LociType = 0x03;
	else
		if(pDGTQTLAlleles->flgQTL)
			LociType = 0x02;
		else
			LociType = 0x01;

	m_OutBuffIdx += sprintf((char *) & m_pszOutBuffer[m_OutBuffIdx], "\"%s\",%d,\"%s\",%d,\"%s\",\"%s\",%d,%.3f,%.3f,\"%s\",%.3f,\"%s\",%.3f,%d,%d,%d,%d,%d,%d\n", (char*)LocateChrom(pDGTQTLAlleles->ChromID), pDGTQTLAlleles->Loci, Diplotype2Txt(RefPBA),LociType, 
				pDGTQTLAlleles->flgQTL ? Diplotype2Txt((0xc0 >> (pDGTQTLAlleles->QTLAlleles[0] * 2))) : "-/-", pDGTQTLAlleles->flgQTL ? Diplotype2Txt((0xc0 >> (pDGTQTLAlleles->QTLAlleles[1] * 2))) : "-/-",
				pDGTQTLAlleles->flgSamplesLoCov ? 0 : 1, Coverage, Grp1Prop, Diplotype2Txt(HiFreqPBA), Grp2Prop, Diplotype2Txt(NxtHiFreqPBA), GrpNProp,
				pDGTQTLAlleles->flgQTLRefMismatch ? 1 : 0,
				pDGTQTLAlleles->flgSamplesRefMismatch ? 1 : 0,
				pDGTQTLAlleles->flgSamplesHomoz ? 1 : 0,
				pDGTQTLAlleles->flgSamplesHeteroz ? 1 : 0,
				pDGTQTLAlleles->flgSamplesPolyAllelic ? 1 : 0,
				pDGTQTLAlleles->flgSamplesMonoAllelic ? 1 : 0);

	if ((m_OutBuffIdx + 1000) > m_AllocOutBuff)
		{
		if (!CUtility::RetryWrites(m_hOutFile, m_pszOutBuffer, m_OutBuffIdx))
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "AnalyseInstance: Fatal error in RetryWrites()");
			Reset();
			return(eBSFerrFileAccess);
			}
		m_OutBuffIdx = 0;
		}
	}
return(eBSFSuccess);
}

int 
CDGTvQTLs::Process(eModeDGTA PMode,		// processing mode
	double MinCoverage,				// if coverage < this threshold then class as being low coverage
	double HomozPropThres,			// if proportion of samples in Grp1Prop is >= this proportion then characterise as homozygous	char* pszAssembRefFile,	// contains orginal reference assembly, diplotype PBA, against which samples were aligned
	char* pszAssembRefFile,			// contains orginal reference assembly, diplotype PBA, against which samples were aligned
	char* pszChromFile,				// BED file, contains chromosome names and sizes
	char* pszDGTsFile,				// file containing DGT loci and allele sample groupings
	char* pszQTLsFile,				// file containing QTL SNP loci
	int NumPBAFiles,				// number of input PBA file specs, these are the source PBA files used during the grouping generation
	char** ppszPBAFiles,			// names of input PBA files (wildcards allowed)
	int	NumIncludeChroms,			// number of chromosome regular expressions to include
	char** ppszIncludeChroms,		// array of include chromosome regular expressions
	int	NumExcludeChroms,			// number of chromosome expressions to exclude
	char** ppszExcludeChroms,		// array of exclude chromosome regular expressions
	char* pszRsltsFileBaseName,		// write results to this file base name - will be suffixed with result types
	int NumThreads)					// number of worker threads to use
{
int Rslt;
size_t MemReq;
Reset();

CreateMutexes();
m_PMode = PMode;
m_MinCoverage = MinCoverage;
m_HomozPropThres = HomozPropThres;
m_pszRsltsFileBaseName = pszRsltsFileBaseName;
m_NumThreads = NumThreads;

if(NumIncludeChroms)
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Process: Loading inclusion chromosomes'");
if (NumExcludeChroms)
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Process: Loading exclusion chromosomes'");
	// compile include/exclude chromosome regexpr if user has specified alignments to be filtered by chrom
if (Rslt = (m_RegExprs.CompileREs(NumIncludeChroms, ppszIncludeChroms, NumExcludeChroms, ppszExcludeChroms)) < eBSFSuccess)
	{
	Reset();
	return(Rslt);
	}

if (m_pInBuffer == nullptr)
	{
	if ((m_pInBuffer = new uint8_t[cInBuffSize]) == nullptr)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Memory allocation of %u bytes for input buffering failed - %s", cInBuffSize, strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
	m_AllocInBuff = cInBuffSize;
	}
m_InNumBuffered = 0;

gDiagnostics.DiagOut(eDLInfo, gszProcName, "Process: Loading chromosome sizes from file: '%s'", pszChromFile);
if ((Rslt = LoadChromSizes(pszChromFile)) < 1) // BED file containing chromosome names and sizes - NOTE chromosomes will be filtered by include/exclude wildcards
	{
	Reset();
	return(Rslt);
	}

if (NumPBAFiles > cMaxPBAFileSpecs)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Can accept at most %d source PBA wildcarded file specs for processing, %d requested", cMaxPBAFileSpecs, NumPBAFiles);
	Reset();
	return(eBSFerrMem);
	}

// initial allocation, will be realloc'd as required if more memory required
MemReq = (size_t)cAllocChromMetadata * sizeof(tsCHChromMetadata);
#ifdef _WIN32
m_pChromMetadata = (tsCHChromMetadata*)malloc(MemReq);	// initial and perhaps the only allocation
if (m_pChromMetadata == nullptr)
{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Initial memory allocation of %zd bytes for chromosome metadata failed - %s", (int64_t)MemReq, strerror(errno));
	Reset();
	return(eBSFerrMem);
}
#else
m_pChromMetadata = (tsCHChromMetadata*)mmap(nullptr, MemReq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
if (m_pChromMetadata == MAP_FAILED)
{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Memory allocation of %zd bytes through mmap() for chromosome metadata failed - %s", (int64_t)MemReq, strerror(errno));
	m_pChromMetadata = nullptr;
	Reset();
	return(eBSFerrMem);
}
#endif
m_AllocdChromMetadataMem = MemReq;
m_UsedNumChromMetadata = 0;
m_AllocdChromMetadata = cAllocChromMetadata;

gDiagnostics.DiagOut(eDLInfo, gszProcName, "Process: Loading DGTs from file: '%s'", pszDGTsFile);
if ((Rslt = LoadDGTs(pszDGTsFile)) < 0)
	{
	Reset();
	return(Rslt);
	}
gDiagnostics.DiagOut(eDLInfo, gszProcName, "Process: Loading QTLs from file: '%s'", pszQTLsFile);
if ((Rslt = LoadQTLs(pszQTLsFile)) < 0)
	{
	Reset();
	return(Rslt);
	}

// sort by ChromID.StartLoci ascending so can easily iterate
gDiagnostics.DiagOut(eDLInfo, gszProcName, "Process: Sorting/rehashing loaded DGTs/QTLs");
m_mtqsort.SetMaxThreads(m_NumThreads);
m_mtqsort.qsort(m_pDGTQTLAlleles, (int64_t)m_UsedDGTQTLAlleles, sizeof(tsDGTQTLAlleles), SortDGTQTLAlleles);
// recreate hashing
uint32_t Hash;
tsDGTQTLAlleles *pDGTQTLAlleles;
memset(m_DGTQTLAllelesHashes, 0, sizeof(m_DGTQTLAllelesHashes));
pDGTQTLAlleles = m_pDGTQTLAlleles;
for (uint32_t Idx = 0; Idx < m_UsedDGTQTLAlleles; Idx++, pDGTQTLAlleles++)
	{
	pDGTQTLAlleles->Nxt = 0;
	Hash = (uint32_t)((pDGTQTLAlleles->ChromID | (pDGTQTLAlleles->Loci << 3)) % cMaxDGTQTLAllelesHashes); // should be reasonably normally distributed
	if (m_DGTQTLAllelesHashes[Hash] != 0)
		pDGTQTLAlleles->Nxt = m_DGTQTLAllelesHashes[Hash];
	else
		pDGTQTLAlleles->Nxt = 0;
	pDGTQTLAlleles->Idx = Idx;
	m_DGTQTLAllelesHashes[Hash] = Idx+1;
	}

// load the reference assembly metadata, assembly must be as a diplotype PBA
gDiagnostics.DiagOut(eDLInfo, gszProcName, "Process: Loading reference diplotype assembly metadata: '%s'", pszAssembRefFile);
int ReadsetID =				// returned readset identifier (1..n) or < 0 if errors
	LoadPBAFile(pszAssembRefFile,	// load chromosome metadata and PBA data from this file, filters out chroms - AcceptThisChromName()
								0,	// 0: founder, 1: progeny, 2: control
				true);  // load chrom metadata (chrom name,length, file offset at which chrom PBAs start) but don't actually load the chromosome PBAs
if (ReadsetID < 1)	// errors loading
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Failed loading reference assembly PBA - %s", pszAssembRefFile);
	Reset();
	return(ReadsetID);
	}

char szOutFile[_MAX_PATH];
CUtility::AppendFileNameSuffix(szOutFile, pszRsltsFileBaseName, (char *)".DGTQTL.csv", '.');

m_Fndrs2Proc[ReadsetID-1] = 0x01;	// if loaded then assumption is that this sample will be processed 
// now load all the sample PBAs
gDiagnostics.DiagOut(eDLInfo, gszProcName, "Processs: Loading all sample PBAs metadata");
Rslt = eBSFSuccess;		// assume success!
CSimpleGlob glob(SG_GLOB_FULLSORT);
int32_t Idx;
char* pszInFile;
int32_t NumFiles;
int32_t TotNumFiles = 0;

// load all samples with out actually allocating memory for each individual chromosome PBAs, but the file offset at which the chromosome PBA starts will be recorded  
for (Idx = 0; Idx < NumPBAFiles; Idx++)
	{
	glob.Init();
	if (glob.Add(ppszPBAFiles[Idx]) < SG_SUCCESS)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessGrpLociDGTs: Unable to glob '%s' sample PBA files", ppszPBAFiles[Idx]);
		Reset();
		return(eBSFerrOpnFile);	// treat as though unable to open file
		}
	if ((NumFiles = glob.FileCount()) <= 0)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessGrpLociDGTs: Unable to locate any sample PBA file matching '%s", ppszPBAFiles[Idx]);
		Reset();
		return(eBSFerrFileName);
		}

	TotNumFiles += NumFiles;

	if (TotNumFiles > cMaxPBAFiles)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Can accept at most %d sample PBA files for processing, after wildcard file name expansions there are %d requested", cMaxPBAFiles, TotNumFiles);
		Reset();
		return(eBSFerrMem);
		}

	Rslt = eBSFSuccess;
	for (int32_t FileID = 0; Rslt >= eBSFSuccess && FileID < NumFiles; ++FileID)
	{
		pszInFile = glob.File(FileID);
		ReadsetID = LoadPBAFile(pszInFile, 1, true); // loading without allocation for chromosome PBAs
		if (ReadsetID <= 0)
			{
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "Process: Errors loading file '%s'", pszInFile);
			Reset();
			return(ReadsetID);
			}
		m_Fndrs2Proc[ReadsetID - 1] = 0x01;	// if loaded then assumption is that this sample will be processed
		}
	}
m_NumFounders = TotNumFiles;
gDiagnostics.DiagOut(eDLInfo, gszProcName, "Process: Completed loading of %d sample PBAs metadata", TotNumFiles);

gDiagnostics.DiagOut(eDLInfo, gszProcName, "Process: Processing DGTs/QTLs against reference and sample PBAs");

if (m_pszOutBuffer == nullptr)
{
	if ((m_pszOutBuffer = new uint8_t[cOutBuffSize]) == nullptr)
	{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Memory allocation of %u bytes for output buffering failed - %s", cOutBuffSize, strerror(errno));
		Reset();
		return(eBSFerrMem);
	}
	m_AllocOutBuff = cOutBuffSize;
}
m_OutBuffIdx = 0;

;
#ifdef _WIN32
m_hOutFile = open(szOutFile, (O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC), (_S_IREAD | _S_IWRITE));
#else
if ((m_hOutFile = open64(szOutFile, O_WRONLY | O_CREAT, S_IREAD | S_IWRITE)) != -1)
if (ftruncate(m_hOutFile, 0) != 0)
{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to create/truncate %s - %s", szOutFile, strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
}
#endif
if (m_hOutFile < 0)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to create/truncate %s - %s", szOutFile, strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}

// header line
m_OutBuffIdx = sprintf((char*)m_pszOutBuffer, "\"Chrom\",\"Loci\",\"RefAlleles\",\"DGT.QTL\",\"QTLRefAlleles\",\"QTLAltAlleles\",\"AcceptCoverage\",\"PropCoverage\",\"PropGrp1Cov\",\"Grp1Alleles\",\"PropGrp2Cov\",\"Grp2Alleles\",\"PropGrpNCov\",\"QTLRefMismatch\",\"SamplesGrpRefMismatch\",\"SamplesHomo\",\"SamplesHetero\",\"SamplesPolyAllelic\",\"SamplesMonoAllelic\"\n");

ProcessDGTQTLsPBAs(PMode);

if (m_OutBuffIdx && m_hOutFile != -1)
	{
	if (!CUtility::RetryWrites(m_hOutFile, m_pszOutBuffer, m_OutBuffIdx))
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Fatal error in RetryWrites()");
		Reset();
		return(eBSFerrFileAccess);
		}
	m_OutBuffIdx = 0;
	}

if (m_hOutFile != -1)
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

Reset();
return(0);
}

// sorting by ChromID.StartLoci ascending
int
CDGTvQTLs::SortDGTQTLAlleles(const void* arg1, const void* arg2)
{
tsDGTQTLAlleles* pEl1 = (tsDGTQTLAlleles*)arg1;
tsDGTQTLAlleles* pEl2 = (tsDGTQTLAlleles*)arg2;
if (pEl1->ChromID < pEl2->ChromID)
	return(-1);
if (pEl1->ChromID > pEl2->ChromID)
	return(1);
if (pEl1->Loci < pEl2->Loci)
	return(-1);
if (pEl1->Loci > pEl2->Loci)
	return(1);
return(0);
}
