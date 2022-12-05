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
#include <process.h>
#include "../libkit4b/commhdrs.h"
#else
#include <sys/mman.h>
#include <pthread.h>
#include "../libkit4b/commhdrs.h"
#endif

#include "./ngskit4b.h"
#include "CGenMLdatasets.h"

#include "../libkit4b/bgzf.h"


int32_t
Process(teModeGMLD PMode,                  // processing mode
				teTypeGMLD FType,                  // input sample feature file format type
				teReductGMLD RMode,                  // feature reduction mode
				char* pszInSampleFeats,     // input sample feature file
				char* pszInSampleLabels,    // input sample labels (classes) file
				char* pszOutMLdataset,     // output ML dataset - row per sample, column per feature, last column containing sample labels
				int32_t	NumThreads);		// maximum number of worker threads to use

#ifdef _WIN32
int genmldatasets(int argc, char *argv[])
{
// determine my process name
_splitpath (argv[0], NULL, NULL, gszProcName, NULL);
#else
int
genmldatasets(int argc, char **argv)
{
// determine my process name
CUtility::splitpath ((char *)argv[0], NULL, gszProcName);
#endif
int32_t iFileLogLevel;			// level of file diagnostics
int32_t iScreenLogLevel;		// level of file diagnostics
char szLogFile[_MAX_PATH];	    // write diagnostics to this file
int32_t Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure
int32_t NumberOfProcessors;		// number of installed CPUs
int32_t NumThreads;				// number of threads (0 defaults to number of CPUs or a maximum of cMaxPBAWorkerThreads)

teModeGMLD PMode;	// processing mode
teTypeGMLD FType;	// input sample feature file format type
teReductGMLD RMode;	// feature reduction mode

char szInSampleFeats[_MAX_PATH]; // input sample feature file
char szInSampleLabels[_MAX_PATH];// input sample labels (classes) file

char szOutMLdataset[_MAX_PATH]; // output ML dataset - row per sample, column per feature, last column containing sample labels

struct arg_lit *help = arg_lit0 ("h", "help", "print this help and exit");
struct arg_lit *version = arg_lit0 ("v", "version,ver", "print version information and exit");
struct arg_int *FileLogLevel = arg_int0 ("f", "FileLogLevel", "<int>", "Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0 ("F", "log", "<file>", "diagnostics log file");

struct arg_int *pmode = arg_int0 ("m", "pmode", "<int>", "processing mode 0: (default 0)");
struct arg_int *ftype = arg_int0 ("M", "ftype", "<int>", "expected sample features input type (default 0)");
struct arg_int *rmode = arg_int0 ("r", "rmode", "<int>", "feature reduction mode 0: (default 0)");

struct arg_file *insamplefeats = arg_file1("i", "infeats", "<file>", "input sample feats file");
struct arg_file *insamplelabels = arg_file0("I", "inlabels", "<file>", "input sample labels file");
struct arg_file *outdataset =  arg_file1("o", "out", "<file>", "output ML dataset file");
struct arg_int *threads = arg_int0("T","threads","<int>","number of processing threads 0..64 (defaults to 0 which limits threads to maximum of 64 CPU cores)");
struct arg_end *end = arg_end (200);

void *argtable[] = { help,version,FileLogLevel,LogFile,
					pmode,ftype,rmode,insamplefeats,insamplelabels,outdataset,end };

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

	PMode = (teModeGMLD)(pmode->count ? pmode->ival[0] : eGMLDDefault);
	if(PMode < 0 || PMode >= eGMLDPlaceHolder)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Unsupported processing mode '-m%d'\n",PMode);
		exit(1);
		}


	FType = (teTypeGMLD)(ftype->count ? ftype->ival[0] : eTGMLDDefault);
	if(FType < 0 || FType >= eTGMLDPlaceHolder)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Unsupported sample feature file format type '-m%d'\n",FType);
		exit(1);
		}

	RMode = (teReductGMLD)(rmode->count ? rmode->ival[0] : eRGMLDDefault);
	if(RMode < 0 || RMode >= eRGMLDPlaceHolder)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Unsupported feature reduction mode '-r%d'\n",RMode);
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

	int MaxAllowedThreads = min(cMaxGMLDWorkerThreads,NumberOfProcessors);	// limit to be at most cMaxWorkerThreads
	if((NumThreads = threads->count ? threads->ival[0] : MaxAllowedThreads)==0)
		NumThreads = MaxAllowedThreads;
	if(NumThreads < 0 || NumThreads > MaxAllowedThreads)
		{
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"Warning: Number of threads '-T%d' specified was outside of range %d..%d",NumThreads,1,MaxAllowedThreads);
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"Warning: Defaulting number of threads to %d",MaxAllowedThreads);
		NumThreads = MaxAllowedThreads;
		}

	if(insamplefeats->count)
		{
		strcpy (szInSampleFeats, insamplefeats->filename[0]);
		CUtility::TrimQuotedWhitespcExtd (szInSampleFeats);
		if(szInSampleFeats[0] == '\0')
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: After removal of whitespace, no input sample features file specified with '-i<filespec>' option)\n");
			exit(1);
			}
		}
	else
		szInSampleFeats[0] = '\0';


	if(insamplelabels->count)
		{
		strcpy (szInSampleLabels, insamplefeats->filename[0]);
		CUtility::TrimQuotedWhitespcExtd (szInSampleLabels);
		if(szInSampleLabels[0] == '\0')
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: After removal of whitespace, no input sample labels file specified with '-I<filespec>' option)\n");
			exit(1);
			}
		}
	else
		szInSampleLabels[0] = '\0';


	if(outdataset->count)
		{
		strcpy (szOutMLdataset, outdataset->filename[0]);
		CUtility::TrimQuotedWhitespcExtd (szOutMLdataset);
		if(szOutMLdataset[0] == '\0')
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: After removal of whitespace, no output file specified with '-o<filespec>' option)\n");
			exit(1);
			}
		}
	else
		szOutMLdataset[0] = '\0';

	gDiagnostics.DiagOut (eDLInfo, gszProcName, "Processing parameters:");
	const char *pszDescr;
	switch (PMode) {
		case eGMLDDefault:
			pszDescr = "Transposition";
			break;
		default:
			pszDescr = "Defaulting to Transposition";
			PMode = eGMLDDefault;
			break;
		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo, "Generating ML dataset : '%s'", pszDescr);
	
	switch (FType) {
		case eTGMLDDefault:
			pszDescr = "Default";
			break;
		default:
			pszDescr = "Defaulting to default";
			FType = eTGMLDDefault;
			break;
		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo, "Input sample feature file type : '%s'", pszDescr);

	switch (RMode) {
		case eRGMLDDefault:
			pszDescr = "Default";
			break;
		default:
			pszDescr = "Defaulting to default";
			RMode = eRGMLDDefault;
			break;
		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo, "Feature reduction : '%s'", pszDescr);

	gDiagnostics.DiagOutMsgOnly(eDLInfo, "Input from sample feature file : '%s'", szInSampleFeats);
	if(szInSampleLabels[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "Input from sample labels file : '%s'", szInSampleLabels);
	gDiagnostics.DiagOutMsgOnly(eDLInfo, "Output to file : '%s'", szOutMLdataset);

	#ifdef _WIN32
	SetPriorityClass (GetCurrentProcess (), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start ();
	Rslt = Process(PMode,                  // processing mode
				   FType,                  // input sample feature file format type
				   RMode,                  // feature reduction mode
				   szInSampleFeats,        // input sample feature file
				   szInSampleLabels,       // input sample labels (classes) file
				   szOutMLdataset,         // output ML dataset - row per sample, column per feature, last column containing sample labels
				   NumThreads);		       // maximum number of worker threads to use
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

int32_t
Process(teModeGMLD PMode,                  // processing mode
			teTypeGMLD FType,              // input sample feature file format type
			teReductGMLD RMode,              // feature reduction mode
			char* pszInSampleFeats,     // input sample feature file
			char* pszInSampleLabels,    // input sample labels (classes) file
			char* pszOutMLdataset,      // output ML dataset - row per sample, column per feature, last column containing sample labels
			int32_t	NumThreads)		    // maximum number of worker threads to use
{
int Rslt;
CGenMLdatasets *pML;
if((pML = new CGenMLdatasets)==NULL)
	{
	gDiagnostics.DiagOut (eDLFatal, gszProcName, "Unable to instantiate instance of CGenMLdatasets");
	return(eBSFerrObj);
	}
Rslt = pML->Process(PMode,FType,RMode,pszInSampleFeats,pszInSampleLabels,pszOutMLdataset,NumThreads);
delete pML;
return(Rslt);
}

CGenMLdatasets::CGenMLdatasets()       // constructor
{
m_pSFFile = nullptr;
m_pszSampleRefsMem = nullptr;
m_pszFeatRefsMem = nullptr;
m_pSampleFeatsMem = nullptr;
Reset();
}

CGenMLdatasets::~CGenMLdatasets()      // destructor
{
if (m_pSFFile != nullptr)
	delete m_pSFFile;

if(m_pszSampleRefsMem != nullptr)
	{
#ifdef _WIN32
	free(m_pszSampleRefsMem);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (m_pszSampleRefsMem != MAP_FAILED)
		munmap(m_pszSampleRefsMem, m_AllocdSampleRefsMem);
#endif
	}

if(m_pszFeatRefsMem != nullptr)
	{
#ifdef _WIN32
	free(m_pszFeatRefsMem);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (m_pszFeatRefsMem != MAP_FAILED)
		munmap(m_pszFeatRefsMem, m_AllocdFeatRefsMem);
#endif
	}


if(m_pSampleFeatsMem != nullptr)
	{
#ifdef _WIN32
	free(m_pSampleFeatsMem);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (m_pSampleFeatsMem != MAP_FAILED)
		munmap(m_pSampleFeatsMem, m_AllocdSampleFeatsMem);
#endif
	}
}

void
CGenMLdatasets::Reset(void)       // reset back to instantiation
{
if (m_pSFFile != nullptr)
	{
	delete m_pSFFile;
	m_pSFFile = nullptr;
	}
if(m_pszSampleRefsMem != nullptr)
	{
#ifdef _WIN32
	free(m_pszSampleRefsMem);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (m_pszSampleRefsMem != MAP_FAILED)
		munmap(m_pszSampleRefsMem, m_AllocdSampleRefsMem);
#endif
	m_pszSampleRefsMem = nullptr;
	}
m_AllocdSampleRefsMem = 0;
m_UsedSampleRefsMem = 0;
m_NumSampleRefs = 0;

if(m_pszFeatRefsMem != nullptr)
	{
#ifdef _WIN32
	free(m_pszFeatRefsMem);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (m_pszFeatRefsMem != MAP_FAILED)
		munmap(m_pszFeatRefsMem, m_AllocdFeatRefsMem);
#endif
	m_pszFeatRefsMem = nullptr;
	}
m_AllocdFeatRefsMem = 0;
m_UsedFeatRefsMem = 0;
m_NumFeatRefs = 0;

if(m_pSampleFeatsMem != nullptr)
	{
#ifdef _WIN32
	free(m_pSampleFeatsMem);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (m_pSampleFeatsMem != MAP_FAILED)
		munmap(m_pSampleFeatsMem, m_AllocdSampleFeatsMem);
#endif
	m_pSampleFeatsMem = nullptr;
	}
m_AllocdSampleFeatsMem = 0;
m_UsedSampleFeatsMem = 0;
m_SampleFeatsSize = sizeof(tsSampleFeats);
m_PMode = eGMLDDefault;
m_FType = eTGMLDDefault;
m_RMode = eRGMLDDefault;
}

int32_t 
CGenMLdatasets::Process(teModeGMLD PMode,  // processing mode
				teTypeGMLD FType,          // input sample feature file format type
				teReductGMLD RMode,          // feature reduction mode
				char* pszInSampleFeats, // input sample feature file
				char* pszInSampleLabels,    // input sample labels (classes) file
				char* pszOutMLdataset,  // output ML dataset - row per sample, column per feature, last column containing sample labels
				int32_t	NumThreads)		// maximum number of worker threads to use
{
int32_t Rslt = eBSFSuccess;
Reset();
m_PMode = PMode;
m_FType = FType;
m_RMode = RMode;

LoadSampleFeatures(FType,pszInSampleFeats);
       ReduceSampleFeatures(RMode);
if(pszInSampleLabels != nullptr && pszInSampleLabels[0] != '\0')
      AssociateSampleLabels(pszInSampleLabels);

Reset();
return(Rslt);
}

int
CGenMLdatasets::AllocateMemSamples(int NumSamples,  // initially allocate for this number of samples, memory will be reallocd to hold more if subsequently required
					   int NumFeatures) // initially allocate for this number of features
	{
	// initial allocations, will be realloc'd to larger sizes if later required
size_t memreq = (size_t)(NumSamples * (cMaxSampleRefNameLen+1));
#ifdef _WIN32
m_pszSampleRefsMem = (uint8_t*)malloc(memreq);	// initial and perhaps the only allocation
if (m_pszSampleRefsMem == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "AllocateSamples: memory allocation of %zd bytes for sample reference names failed - %s", (int64_t)memreq, strerror(errno));
	Reset();
	return(eBSFerrMem);
	}
#else
m_pszSampleRefsMem = (uint8_t *)mmap(NULL, memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
if (m_pszSampleRefsMem == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "AllocateSamples: Memory allocation of %zd bytes through mmap() for sample reference names failed - %s", (int64_t)memreq, strerror(errno));
	m_pszSampleRefsMem = NULL;
	Reset();
	return(eBSFerrMem);
	}
#endif
m_AllocdSampleRefsMem = memreq;
m_NumSampleRefs = 0;
m_UsedSampleRefsMem = 0;

memreq = (size_t)(NumFeatures * (cMaxFeatRefNameLen+1));
#ifdef _WIN32
m_pszFeatRefsMem = (uint8_t*)malloc(memreq);	// initial and perhaps the only allocation
if (m_pszFeatRefsMem == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "AllocateSamples: Memory allocation of %zd bytes for feature reference names failed - %s", (int64_t)memreq, strerror(errno));
	Reset();
	return(eBSFerrMem);
	}
#else
m_pszFeatRefsMem = (uint8_t *)mmap(NULL, memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
if (m_pszFeatRefsMem == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "AllocateSamples: Memory allocation of %zd bytes through mmap() for feature reference names failed - %s", (int64_t)memreq, strerror(errno));
	m_pszFeatRefsMem = NULL;
	Reset();
	return(eBSFerrMem);
	}
#endif
m_AllocdFeatRefsMem = memreq;
m_NumFeatRefs = 0;
m_UsedFeatRefsMem = 0;

m_SampleFeatsSize = sizeof(tsSampleFeats);
memreq = (size_t)NumSamples * ((size_t)NumFeatures + m_SampleFeatsSize);
#ifdef _WIN32
m_pSampleFeatsMem = (uint8_t*)malloc(memreq);	// initial and perhaps the only allocation
if (m_pSampleFeatsMem == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "AllocateSamples: Memory allocation of %zd bytes for sample features failed - %s", (int64_t)memreq, strerror(errno));
	Reset();
	return(eBSFerrMem);
	}
#else
m_pSampleFeatsMem = (uint8_t *)mmap(NULL, memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
if (m_pSampleFeatsMem == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "AllocateSamples: Memory allocation of %zd bytes through mmap() for sample features failed - %s", (int64_t)memreq, strerror(errno));
	m_pSampleFeatsMem = NULL;
	Reset();
	return(eBSFerrMem);
	}
#endif
m_AllocdSampleFeatsMem = memreq;
m_UsedSampleFeatsMem = 0;

return(eBSFSuccess);
}


int
CGenMLdatasets::ReallocSamples(void)      // realloc as may be required to have room for at least one more sample reference name and associated sample features
{

uint8_t *pReallocd;
		// needing to allocate more memory?
if ((m_UsedSampleRefsMem + (size_t)cMaxSampleRefNameLen + 1) >= m_AllocdSampleRefsMem)
	{
	size_t memreq = m_AllocdSampleRefsMem + (size_t)(1000 * cMaxSampleRefNameLen);
#ifdef _WIN32
	pReallocd = (uint8_t *)realloc(m_pszSampleRefsMem, memreq);
	if (pReallocd == NULL)
		{
#else
	pReallocd = (uint8_t *)mremap(m_pszSampleRefsMem, m_AllocdSampleRefsMem, memreq, MREMAP_MAYMOVE);
	if (pReallocd == MAP_FAILED)
		{
#endif
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "ReallocSamples: Memory reallocation to %zd bytes failed - %s", (int64_t)memreq, strerror(errno));
		return(eBSFerrMem);
		}
	m_pszSampleRefsMem = pReallocd;
	m_AllocdSampleRefsMem = memreq;
	}

if ((m_UsedSampleFeatsMem + (size_t)m_SampleFeatsSize) >= m_AllocdSampleFeatsMem)
	{
	size_t memreq = m_AllocdSampleFeatsMem + (size_t)(1000 * m_SampleFeatsSize);
#ifdef _WIN32
	pReallocd = (uint8_t *)realloc(m_pSampleFeatsMem, memreq);
	if (pReallocd == NULL)
		{
#else
	pReallocd = (uint8_t *)mremap(m_pSampleFeatsMem, m_AllocdSampleFeatsMem, memreq, MREMAP_MAYMOVE);
	if (pReallocd == MAP_FAILED)
		{
#endif
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "ReallocSamples: Memory reallocation to %zd bytes failed - %s", (int64_t)memreq, strerror(errno));
		return(eBSFerrMem);
		}
	m_pSampleFeatsMem = pReallocd;
	m_AllocdSampleFeatsMem = memreq;
	}
return(eBSFSuccess);
}

int32_t 
CGenMLdatasets::AddSampleRefName(char *pszSampleRef) // returns the sample name identifier, 1..N
{
char *pszName;
int32_t SampleRefsOfs;
SampleRefsOfs = (int32_t)m_UsedSampleRefsMem;
pszName = (char *) & m_pszSampleRefsMem[m_UsedSampleRefsMem];
strcpy(pszName,pszSampleRef);
m_UsedSampleRefsMem += strlen(pszSampleRef) + 1;
m_NumSampleRefs += 1;
return(m_NumSampleRefs);
}

int32_t 
CGenMLdatasets::AddFeatRefName(char *pszFeatRef)  // returns the feature name identifier, 1..N
{
char *pszName;
int32_t FeatsRefsOfs;
FeatsRefsOfs = (int32_t)m_UsedFeatRefsMem;
pszName = (char *) & m_pszFeatRefsMem[m_UsedFeatRefsMem];
strcpy(pszName,pszFeatRef);
m_UsedFeatRefsMem += strlen(pszFeatRef) + 1;
m_NumFeatRefs += 1;
return(m_NumFeatRefs);
}

int32_t
CGenMLdatasets::AddSampleFeatValue(int32_t SampleID,  // sample
				   int32_t FeatID,    // feature
				   int32_t Value) // has this value
{
tsSampleFeatValue SampleFeatValue;
SampleFeatValue.SampleID = SampleID;
SampleFeatValue.FeatID = FeatID;
SampleFeatValue.Value = Value;
return(0);
}

int
CGenMLdatasets::LoadSampleFeatures(teTypeGMLD FType,          // input sample feature file format type
	char* pszInSampleFeats) // load sample features from this file
{
uint32_t EstNumRows;
int64_t FileSize;
int32_t MaxFields;
int32_t MeanNumFields;
int FieldIdx;
int32_t Rslt = 0;
int32_t CurLineNumber;
int32_t NumFields;
int32_t LastSampleRefIdx;
int32_t SampleRefNameOfs;
char *pszSampleRef;
char szFeatRef[200];
char *pszFeatChrom;
char *pszFeatLoci;
int32_t FeatValue;
int32_t FeatID;

if(m_pSFFile != nullptr) // shouldn't have been instantiated, but better to be sure!
	{
	delete m_pSFFile;
	m_pSFFile = nullptr;
	}
if((m_pSFFile = new CCSVFile) == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to instantiate CCSVfile");
	Reset();
	return(eBSFerrObj);
	}

	// get an estimate of number of rows and fields
if ((EstNumRows = m_pSFFile->CSVEstSizes(pszInSampleFeats, &FileSize, &MaxFields, &MeanNumFields)) < 2)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to estimate number of rows in file: '%s'", pszInSampleFeats);
	Reset();
	return(eBSFerrFieldCnt);
	}

m_pSFFile->SetMaxFields(MaxFields);


AllocateMemSamples(MaxFields,EstNumRows); // CSV fields are per sample, rows are per feature

if((Rslt = m_pSFFile->Open(pszInSampleFeats)) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to open file: '%s'", pszInSampleFeats);
	Reset();
	return(Rslt);
	}

// header line contains sample identifiers
// rows contain features associated with samples
CurLineNumber = 0;
LastSampleRefIdx = 0;
FeatID = 0;
while((Rslt = m_pSFFile->NextLine()) > 0)		// onto next line containing fields
	{
	CurLineNumber++;
	if((NumFields = m_pSFFile->GetCurFields()) < 11)	// must contain at least 11 fields
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Haplotype group bin specification file '%s' expected to contain a minimum of 11 fields, it contains %d at line %d", pszInSampleFeats, NumFields, CurLineNumber);
		Reset();
		return(eBSFerrParse);
		}
	if (CurLineNumber == 1) // 1 if header containing sample identifiers
		{
		// parse in sample identifiers until last column
		for(FieldIdx = 11; FieldIdx <= NumFields; FieldIdx++)
			{
			m_pSFFile->GetText(FieldIdx, &pszSampleRef);
			if(!strncmp((char *)"GrpMembers:1", pszSampleRef,37))
				{
				LastSampleRefIdx = FieldIdx;
				break;
				}
			// sample reference name parsed out
			SampleRefNameOfs = AddSampleRefName(pszSampleRef);
			}
		continue;
		}

	ReallocSamples();	// ensure sufficient memory preallocd to hold at least one other sample with feature values
	// each row contains the features
	// 
	// parse in features to associate with sample references
	m_pSFFile->GetText(3, &pszFeatChrom);
	m_pSFFile->GetText(4, &pszFeatLoci);
	sprintf(szFeatRef,"%s:%s",pszFeatChrom,pszFeatLoci);
	FeatID = AddFeatRefName(szFeatRef);
	int32_t SampleID = 1;
	for (FieldIdx = 11; FieldIdx < LastSampleRefIdx; FieldIdx++,SampleID++)
		{
		m_pSFFile->GetInt(FieldIdx, &FeatValue);
		AddSampleFeatValue(SampleID,FeatID,FeatValue);
		}
	}
delete m_pSFFile;
m_pSFFile = nullptr;
return(Rslt);
}

int
CGenMLdatasets::ReduceSampleFeatures(teReductGMLD RMode)          // feature reduction mode
{
int32_t Rslt = 0;

return(Rslt);
}

int
CGenMLdatasets::AssociateSampleLabels(char* pszInSampleLabels) // associate samples with these labels
{
int32_t Rslt = 0;

return(Rslt);
}