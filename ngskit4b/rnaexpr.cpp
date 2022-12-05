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
#include "rnaexpr.h"

// forward declaration
int procrnaexpr(eModeRNAE PMode,			// processing mode
				char* pszMaterialRNAGBSWGSFile, // load RNA/GBS/WIG sample names and RNA associated metadata
				char* pszInCntsFile,		// load coverage counts from this file
				char* pszInWGSvWGSFile,		// load WGS vs. WGS homozygosity scores from this file
				char* pszInRNAvWGSFile,		// load RNA vs. WGS homozygosity scores from this file
				char* pszOutRslts);			// write results to this file, will be suffixed appropriately

#ifdef _WIN32
int rnaexpr(int argc, char *argv[])
{
// determine my process name
_splitpath (argv[0], NULL, NULL, gszProcName, NULL);
#else
int
rnaexpr(int argc, char **argv)
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

eModeRNAE PMode;			// processing mode
char szMaterialRNAGBSWGSFile[_MAX_PATH]; // load RNA/GBS/WIG sample names and RNA associated metadata
char szInCntsFile[_MAX_PATH];	// load expression  level counts from this file
char szInWGSvWGSFile[_MAX_PATH];	// load WGS vs. WGS homozygosity scores from this file
char szInRNAvWGSFile[_MAX_PATH];	// load RNA vs. WGS homozygosity scores from this file

char szOutRsltsFile[_MAX_PATH];		// write results to this file

struct arg_lit *help = arg_lit0 ("h", "help", "print this help and exit");
struct arg_lit *version = arg_lit0 ("v", "version,ver", "print version information and exit");
struct arg_int *FileLogLevel = arg_int0 ("f", "FileLogLevel", "<int>", "Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0 ("F", "log", "<file>", "diagnostics log file");

struct arg_int *pmode = arg_int0 ("m", "mode", "<int>", "processing mode 0: RNA replicates inconsistencies 1: RNA replicate WGS homozygosity score matching (default 0)");
struct arg_file* cntsfile = arg_file0("i", "cntsfile", "<file>", "input RNA expression level counts file");
struct arg_file* samplesfile = arg_file1("c", "samplesfile", "<file>", "input RNA/GBS/WIG sample names and RNA associated metadata file");
struct arg_file* wgsscorefile = arg_file0("w", "wgsscorefile", "<file>", "input WGS vs. WGS homozygosity scores file");
struct arg_file* rnascorefile = arg_file0("r", "rnascorefile", "<file>", "input RNA vs. WGS homozygosity scores file");
struct arg_file* rsltsfile = arg_file1("o", "rsltfile", "<file>", "output results base file name, will be suffixed with processing mode specific extension");
struct arg_int *threads = arg_int0("T","threads","<int>","number of processing threads 0..64 (defaults to 0 which limits threads to maximum of 64 CPU cores)");
struct arg_end *end = arg_end (200);

void *argtable[] = { help,version,FileLogLevel,LogFile,
					pmode,cntsfile,wgsscorefile,rnascorefile,samplesfile,rsltsfile,
					threads,end };

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

	PMode = pmode->count ? (eModeRNAE)pmode->ival[0] : eMRNAEDefault;
	if(PMode < eMRNAEDefault || PMode >= eMRNAEPlaceHolder)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Unsupported processing mode '-m%d'\n",PMode);
		exit(1);
		}

	szInCntsFile[0] = '\0';
	szInWGSvWGSFile[0] = '\0';
	szInRNAvWGSFile[0] = '\0';
	if (cntsfile->count)
		{
		strncpy(szInCntsFile, cntsfile->filename[0], _MAX_PATH);
		szInCntsFile[_MAX_PATH-1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(szInCntsFile);
		}
	else
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "No input RNA expression level counts file specified");
		exit(1);
		}
	if(PMode == eMRNAEHomozScores)
		{
		if (wgsscorefile->count)
			{
			strncpy(szInWGSvWGSFile, wgsscorefile->filename[0], _MAX_PATH);
			szInWGSvWGSFile[_MAX_PATH-1] = '\0';
			CUtility::TrimQuotedWhitespcExtd(szInWGSvWGSFile);
			}
		else
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "No input WGS vs. WGS homozygosity scores file specified");
			exit(1);
			}		

		if (rnascorefile->count)
			{
			strncpy(szInRNAvWGSFile, rnascorefile->filename[0], _MAX_PATH);
			szInRNAvWGSFile[_MAX_PATH-1] = '\0';
			CUtility::TrimQuotedWhitespcExtd(szInRNAvWGSFile);
			}
		else
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "No input RNA vs. WGS homozygosity scores file specified");
			exit(1);
			}
		}



	if (samplesfile->count)
		{
		strncpy(szMaterialRNAGBSWGSFile, samplesfile->filename[0], _MAX_PATH);
		szMaterialRNAGBSWGSFile[_MAX_PATH-1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(szMaterialRNAGBSWGSFile);
		}
	else
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "No input RNA/GBS/WIG sample names and RNA associated metadata file specified");
		exit(1);
		}

	if (rsltsfile->count)
		{
		strncpy(szOutRsltsFile, rsltsfile->filename[0], _MAX_PATH);
		szOutRsltsFile[_MAX_PATH-1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(szOutRsltsFile);
		}
	else
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "No output results file specified");
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

	int MaxAllowedThreads = min(cMaxRNAEWorkerThreads,NumberOfProcessors);	// limit to be at most cMaxWorkerThreads
	if((NumThreads = threads->count ? threads->ival[0] : MaxAllowedThreads)==0)
		NumThreads = MaxAllowedThreads;
	if(NumThreads < 0 || NumThreads > MaxAllowedThreads)
		{
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"Warning: Number of threads '-T%d' specified was outside of range %d..%d",NumThreads,1,MaxAllowedThreads);
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"Warning: Defaulting number of threads to %d",MaxAllowedThreads);
		NumThreads = MaxAllowedThreads;
		}

	gDiagnostics.DiagOut (eDLInfo, gszProcName, "Processing parameters:");
	const char *pszDescr;
	switch (PMode) {
			case eMRNAEDefault:	
				pszDescr = "check biological sample replicates for expression level counts inconsistencies";
				break;
			case eMRNAEHomozScores:	
				pszDescr = "check biological sample replicates for allelic homozygosity score inconsistencies";
				break;
			}

	gDiagnostics.DiagOutMsgOnly (eDLInfo, "RNA expression analysis mode: '%s'", pszDescr);
	gDiagnostics.DiagOutMsgOnly(eDLInfo, "Sample identifiers and RNA metadata : '%s'", szMaterialRNAGBSWGSFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo, "Input RNA expression levels counts file : '%s'", szInCntsFile);
	if(PMode == eMRNAEHomozScores)
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "input WGS vs. WGS homozygosity scores file : '%s'", szInWGSvWGSFile);
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "input RNA vs. WGS homozygosity scores file : '%s'", szInRNAvWGSFile);
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo, "Output base file name : '%s'", szOutRsltsFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"number of threads : %d",NumThreads);

#ifdef _WIN32
	SetPriorityClass (GetCurrentProcess (), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start ();
	Rslt = 0;
	Rslt = procrnaexpr(PMode,					// processing mode
						szMaterialRNAGBSWGSFile, // load RNA/GBS/WIG sample names and RNA associated metadata
						szInCntsFile,			// load coverage counts from this file
						szInWGSvWGSFile,		// load WGS vs. WGS homozygosity scores from this file
						szInRNAvWGSFile,		// load RNA vs. WGS homozygosity scores from this file
						szOutRsltsFile);		// write results to this file
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
return(0);
}

int
procrnaexpr(eModeRNAE PMode,			// processing mode
				char* pszMaterialRNAGBSWGSFile, // load RNA/GBS/WIG sample names and RNA associated metadata
				char* pszInCntsFile,		// load coverage counts from this file
				char* pszInWGSvWGSFile,		// load WGS vs. WGS homozygosity scores from this file
				char* pszInRNAvWGSFile,		// load RNA vs. WGS homozygosity scores from this file
				char* pszOutRsltsFile)			// write results to this file, will be suffixed appropriately
{
int Rslt;
CRNAExpr *pCRNAExpr;
if((pCRNAExpr = new CRNAExpr) == NULL)
	{
	gDiagnostics.DiagOut (eDLFatal, gszProcName, "Unable to instantiate instance of CRNAExpr");
	return(eBSFerrObj);
	}
Rslt = pCRNAExpr->Process(PMode,			// processing mode
					pszMaterialRNAGBSWGSFile, // load RNA/GBS/WIG sample names and RNA associated metadata
					pszInCntsFile,	// load coverage counts from this file
					pszInWGSvWGSFile,		// load WGS vs. WGS homozygosity scores from this file
					pszInRNAvWGSFile,		// load RNA vs. WGS homozygosity scores from this file
					pszOutRsltsFile);		// write results to this file);
delete pCRNAExpr;

return(Rslt);
}


CRNAExpr::CRNAExpr()	// constructor
{
m_pRNAGBSWGSFile = nullptr;
m_pInCntsFile = nullptr;
m_pInWGSvWGSFile = nullptr;
m_pInRNAvWGSFile = nullptr;
m_pRNACntsMem = nullptr;
m_pRNAhomozScoresMem = nullptr;
m_pWGShomozScoresMem = nullptr;
Reset();
}

CRNAExpr::~CRNAExpr()	// destructor
{
if(m_pRNAGBSWGSFile != nullptr)
	delete m_pRNAGBSWGSFile;
if(m_pInCntsFile != nullptr)
	delete m_pInCntsFile;
if(m_pInWGSvWGSFile != nullptr)
	delete m_pInWGSvWGSFile;
if(m_pInRNAvWGSFile != nullptr)
	delete m_pInRNAvWGSFile;

if(m_pRNACntsMem != nullptr)
	{
#ifdef _WIN32
	free(m_pRNACntsMem);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (m_pRNACntsMem != MAP_FAILED)
		munmap(m_pRNACntsMem, m_AllocdRNACntsMem);
#endif
	}

if(m_pRNAhomozScoresMem != nullptr)
	{
#ifdef _WIN32
	free(m_pRNAhomozScoresMem);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (m_pRNAhomozScoresMem != MAP_FAILED)
		munmap(m_pRNAhomozScoresMem, m_AllocdRNAhomozScoresMem);
#endif
	}

if(m_pWGShomozScoresMem != nullptr)
	{
#ifdef _WIN32
	free(m_pWGShomozScoresMem);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (m_pWGShomozScoresMem != MAP_FAILED)
		munmap(m_pWGShomozScoresMem, m_AllocdWGShomozScoresMem);
#endif
	}

}

void 
CRNAExpr::Reset(void)					// reset instance state back to that immediately following instantiation
{
if (m_pInCntsFile != nullptr)
	{
	delete m_pInCntsFile;
	m_pInCntsFile = nullptr;
	}

if (m_pInWGSvWGSFile != nullptr)
	{
	delete m_pInWGSvWGSFile;
	m_pInWGSvWGSFile = nullptr;
	}

if (m_pInRNAvWGSFile != nullptr)
	{
	delete m_pInRNAvWGSFile;
	m_pInRNAvWGSFile = nullptr;
	}

if (m_pRNAGBSWGSFile != nullptr)
	{
	delete m_pRNAGBSWGSFile;
	m_pRNAGBSWGSFile = nullptr;
	}

if(m_pRNACntsMem != nullptr)
	{
#ifdef _WIN32
	free(m_pRNACntsMem);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (m_pRNACntsMem != MAP_FAILED)
		munmap(m_pRNACntsMem, m_AllocdRNACntsMem);
#endif
	m_pRNACntsMem = nullptr;
	}
m_AllocdRNACntsMem = 0;
m_UsedRNACntsMem = 0;

if(m_pRNAhomozScoresMem != nullptr)
	{
#ifdef _WIN32
	free(m_pRNAhomozScoresMem);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (m_pRNAhomozScoresMem != MAP_FAILED)
		munmap(m_pRNAhomozScoresMem, m_AllocdRNAhomozScoresMem);
#endif
	m_pRNAhomozScoresMem = nullptr;
	}
m_UsedRNAhomozScoresMem = 0;
m_AllocdRNAhomozScoresMem = 0;

if(m_pWGShomozScoresMem != nullptr)
	{
#ifdef _WIN32
	free(m_pWGShomozScoresMem);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (m_pWGShomozScoresMem != MAP_FAILED)
		munmap(m_pWGShomozScoresMem, m_AllocdWGShomozScoresMem);
#endif
	m_pWGShomozScoresMem = nullptr;
	}
m_UsedWGShomozScoresMem = 0;
m_AllocdWGShomozScoresMem = 0;

m_NumbHdrRNAMappings = 0;
memset(m_HdrRNAMappings,0,sizeof(m_HdrRNAMappings));

m_NumbHdrWGSMappings = 0;
memset(m_HdrWGSMappings,0,sizeof(m_HdrWGSMappings));

m_LARNASampleNameID = 0;
m_NumRNASampleNames =0;
m_NxtszRNASampleNameIdx= 0;
m_szRNASampleNames[0] = 0;
memset(m_szRNASampleNames,0,sizeof(m_szRNASampleNames));
memset(m_szRNASampleNamesIdx,0,sizeof(m_szRNASampleNamesIdx));

m_LAWGSSampleNameID = 0;
m_NumWGSSampleNames =0;
m_NxtszWGSSampleNameIdx= 0;
m_szWGSSampleNames[0] = 0;
memset(m_szWGSSampleNames,0,sizeof(m_szWGSSampleNames));
memset(m_szWGSSampleNamesIdx,0,sizeof(m_szWGSSampleNamesIdx));

m_LAMaterialNameID = 0;
m_NumMaterialNames =0;
m_NxtszMaterialNameIdx= 0;
m_szMaterialNames[0] = 0;
memset(m_szMaterialNames,0,sizeof(m_szMaterialNames));
memset(m_szMaterialNamesIdx,0,sizeof(m_szMaterialNamesIdx));


m_LAFeatureNameID = 0;
m_NumFeatureNames = 0;
m_NxtszFeatureNameIdx = 0;
memset(m_szFeatureNames,0,sizeof(m_szFeatureNames));
memset(m_szFeatureNamesIdx,0,sizeof(m_szFeatureNamesIdx));

m_LAChromNameID = 0;
m_NumChromNames = 0;
m_NxtszChromNameIdx = 0;
memset(m_szChromNames,0,sizeof(m_szChromNames));
memset(m_szChromNamesIdx,0,sizeof(m_szChromNamesIdx));
}


int
CRNAExpr::LoadRNACntsFile(char* pszInCntsFile,	// load coverage counts from this file
						 bool bNormaliseCnts)	// normalise individual RNA feature counts, aka library size, to maximal total counts of any replicate
{
int32_t Rslt;
uint32_t EstNumRows;
int64_t FileSize;
int32_t MaxFields;
int32_t MeanNumFields;
int FieldIdx;
int32_t CurLineNumber;
int32_t NumFields;
int32_t ExpNumFields;
int32_t SampleNameID;
char *pszSampleRef;
char *pszFeatName;
char *pszFeatChrom;
char FeatStrand;
int32_t *pSampleValue;
int32_t FeatStartLoci;
int32_t FeatEndLoci;
int32_t NumExons;
int32_t FeatLength;
int32_t TransLength;
int32_t FeatValue;
int32_t FeatureNameID;

if(m_pInCntsFile != nullptr) // shouldn't have been instantiated, but better to be sure!
	{
	delete m_pInCntsFile;
	m_pInCntsFile = nullptr;
	}
if((m_pInCntsFile = new CCSVFile) == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to instantiate CCSVfile");
	Reset();
	return(eBSFerrObj);
	}

	// get an estimate of number of rows and fields
if ((EstNumRows = m_pInCntsFile->CSVEstSizes(pszInCntsFile, &FileSize, &MaxFields, &MeanNumFields)) < 2)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to estimate number of rows in file: '%s'", pszInCntsFile);
	Reset();
	return(eBSFerrFieldCnt);
	}

if((Rslt=AllocateRNAcnts(EstNumRows*MaxFields)) != eBSFSuccess)
	{
	Reset();
	return(Rslt);
	}

m_pInCntsFile->SetMaxFields(MaxFields);

if((Rslt = m_pInCntsFile->Open(pszInCntsFile)) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to open file: '%s'", pszInCntsFile);
	Reset();
	return(Rslt);
	}

// header line contains sample identifiers
// rows contain features associated with samples
CurLineNumber = 0;
ExpNumFields = 0;
FeatureNameID = 0;
m_NumbHdrRNAMappings = 0;
while((Rslt = m_pInCntsFile->NextLine()) > 0)		// onto next line containing fields
	{
	CurLineNumber++;
	if((NumFields = m_pInCntsFile->GetCurFields()) < 9)	// must contain at least 9 fields (assumes counts for at least 1 sample!)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Input counts file '%s' expected to contain a minimum of 9 fields, it contains %d at line %d", pszInCntsFile, NumFields, CurLineNumber);
		Reset();
		return(eBSFerrParse);
		}

	if (CurLineNumber == 1) // 1 if header containing sample identifiers
		{
		if (!m_pInCntsFile->IsLikelyHeaderLine())
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Input counts file '%s' line 1 does not parse as a header line, all fields are expected to be double quoted as containing string values", pszInCntsFile);
			Reset();
			return(eBSFerrParse);
			}

		ExpNumFields = NumFields;
		// parse out sample identifiers until last column
		for(FieldIdx = 9; FieldIdx <= NumFields; FieldIdx++)
			{
			m_pInCntsFile->GetText(FieldIdx, &pszSampleRef);
			if ((SampleNameID = AddRNASampleName(pszSampleRef)) != 0) // RNA sample name must already be known from the material mapping file
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Input counts file '%s' contains RNA sample name '%s' not in materials file at field %d at line %d", pszInCntsFile, pszSampleRef,FieldIdx, CurLineNumber);
				Reset();
				return(eBSFerrParse);
				}
			SampleNameID = LocateRNASampleNameID(pszSampleRef);
			m_HdrRNAMappings[m_NumbHdrRNAMappings++] = SampleNameID;
			}
		continue;
		}
	if((NumFields = m_pInCntsFile->GetCurFields()) != ExpNumFields)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Input counts file '%s' expected to contain same number of fields as header line (%d), it contains %d at line %d", pszInCntsFile,ExpNumFields, NumFields, CurLineNumber);
		Reset();
		return(eBSFerrParse);
		}


	if((Rslt = ReallocRNAcnts(m_NumbHdrRNAMappings))!=eBSFSuccess)	// ensure sufficient memory preallocd to hold at least one additional row of feature values
		{
		Reset();
		return(Rslt);
		}

	// each row contains the features
	// parse in features to associate with sample references
	if (m_NumFeatureNames == cMaxRNAEFeatures)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Input counts file '%s' expected to contain a maximum of %d features at line %d", pszInCntsFile,cMaxRNAEFeatures, CurLineNumber);
		Reset();
		return(eBSFerrParse);
		}
	m_pInCntsFile->GetText(1, &pszFeatName);
	if ((FeatureNameID = AddFeatureName(pszFeatName)) == 0)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Input counts file '%s' expected to contain unique features, duplicate feature '%s' at line %d", pszInCntsFile,pszFeatName, CurLineNumber);
		Reset();
		return(eBSFerrParse);
		}

	m_pInCntsFile->GetText(2, &pszFeatChrom);
	m_pInCntsFile->GetInt(3, &FeatStartLoci);
	m_pInCntsFile->GetInt(4, &FeatEndLoci);
	m_pInCntsFile->GetChar(5, &FeatStrand);
	m_pInCntsFile->GetInt(6,&NumExons);
	m_pInCntsFile->GetInt(7,&FeatLength);
	m_pInCntsFile->GetInt(8,&TransLength);

	pSampleValue = (int32_t *)&m_pRNACntsMem[((size_t)(FeatureNameID-1) * m_NumbHdrRNAMappings) * sizeof(int32_t)];
	for (FieldIdx = 9; FieldIdx <= ExpNumFields; FieldIdx++,pSampleValue++)
		{
		m_pInCntsFile->GetInt(FieldIdx, &FeatValue);
		*pSampleValue = FeatValue;
		}
	}
if (bNormaliseCnts && m_NumbHdrRNAMappings > 0 && m_NumFeatureNames > 0)
	{
	int32_t RNAIdx;
	int32_t FeatIdx;
	int32_t MaxRepTotalsRNAIdx;
	int64_t RepTotals[cMaxRNAESamples];

	// firstly total up counts for each replicate over all it's features into RepTotals
	// also noting which replicate has the maximal number of total counts
	MaxRepTotalsRNAIdx = 0;
	for(RNAIdx = 0; RNAIdx < m_NumbHdrRNAMappings;RNAIdx++)
		{
		RepTotals[RNAIdx] = 0;
		pSampleValue = (int32_t *)&m_pRNACntsMem[((size_t)RNAIdx * sizeof(int32_t))];
		for(FeatIdx = 0; FeatIdx < m_NumFeatureNames; FeatIdx++,pSampleValue += m_NumbHdrRNAMappings)
			RepTotals[RNAIdx] += *pSampleValue;
		if(RepTotals[RNAIdx] > RepTotals[MaxRepTotalsRNAIdx])
			MaxRepTotalsRNAIdx = RNAIdx;
		}
	// totals for each replicate are now known, normalise by scaling each replicate feature counts such that replicate counts will now sum to approximately RepTotals[MaxRepTotalsRNAIdx]
	for(RNAIdx = 0; RNAIdx < m_NumbHdrRNAMappings;RNAIdx++)
		{
		double ScaleFact = (double)RepTotals[MaxRepTotalsRNAIdx] / RepTotals[RNAIdx];
		pSampleValue = (int32_t *)&m_pRNACntsMem[((size_t)RNAIdx * sizeof(int32_t))];
		for(FeatIdx = 0; FeatIdx < m_NumFeatureNames; FeatIdx++,pSampleValue += m_NumbHdrRNAMappings)
			*pSampleValue = (int32_t)(*pSampleValue * ScaleFact);
		}
	}

delete m_pInCntsFile;
m_pInCntsFile = nullptr;

return((m_NumbHdrRNAMappings > 0 && m_NumFeatureNames > 0) ? eBSFSuccess : eBSFerrParse);
}


// Load sample material mappings
int
CRNAExpr::LoadMaterialRNAGBSWGSFile(char* pszMaterialRNAGBSWGSFile)
{
int32_t Rslt;
uint32_t EstNumRows;
int64_t FileSize;
int32_t MaxFields;
int32_t MeanNumFields;
int32_t CurLineNumber;
int32_t NumFields;
int32_t ExpNumFields;
char *pszHdrField;
char *pszValue;
int NumbRNASamples;
uint32_t RNASampleID;
uint32_t WGSSampleID;
uint32_t MaterialID;

if(m_pRNAGBSWGSFile != nullptr) // shouldn't have been instantiated, but better to be sure!
	{
	delete m_pRNAGBSWGSFile;
	m_pRNAGBSWGSFile = nullptr;
	}
if((m_pRNAGBSWGSFile = new CCSVFile) == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to instantiate CCSVfile");
	Reset();
	return(eBSFerrObj);
	}

	// get an estimate of number of rows and fields
if ((EstNumRows = m_pRNAGBSWGSFile->CSVEstSizes(pszMaterialRNAGBSWGSFile, &FileSize, &MaxFields, &MeanNumFields)) < 2)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to estimate number of rows in file: '%s'", pszMaterialRNAGBSWGSFile);
	Reset();
	return(eBSFerrFieldCnt);
	}


m_pRNAGBSWGSFile->SetMaxFields(MaxFields);

if((Rslt = m_pRNAGBSWGSFile->Open(pszMaterialRNAGBSWGSFile)) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to open file: '%s'", pszMaterialRNAGBSWGSFile);
	Reset();
	return(Rslt);
	}

CurLineNumber = 0;
ExpNumFields = 0;
NumbRNASamples = 0;
while((Rslt = m_pRNAGBSWGSFile->NextLine()) > 0)		// onto next line containing fields
	{
	CurLineNumber++;
	if((NumFields = m_pRNAGBSWGSFile->GetCurFields()) < 11)	// must contain at least 11 fields
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Input materials mapping file '%s' expected to contain a minimum of 11 fields, it contains %d at line %d", pszMaterialRNAGBSWGSFile, NumFields, CurLineNumber);
		Reset();
		return(eBSFerrParse);
		}

	if (CurLineNumber == 1) // 1 if header containing sample identifiers
		{
		if (!m_pRNAGBSWGSFile->IsLikelyHeaderLine())
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Input materials mapping file '%s' line 1 does not parse as a header line, all fields are expected to be double quoted as containing string values", pszMaterialRNAGBSWGSFile);
			Reset();
			return(eBSFerrParse);
			}

		ExpNumFields = NumFields;
		m_pRNAGBSWGSFile->GetText(1, &pszHdrField);		// expected to be RNASampleName
		if(stricmp(pszHdrField,"RNASampleName"))
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Input materials mapping file '%s' header does not have expected field %d named as '%', instead '%s' was parsed",pszMaterialRNAGBSWGSFile,1,"RNASampleName",pszHdrField);
			Reset();
			return(eBSFerrParse);
			}

		m_pRNAGBSWGSFile->GetText(3, &pszHdrField);		// expected to be WGSSampleName
		if(stricmp(pszHdrField,"WGSSampleName"))
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Input materials mapping file '%s' header does not have expected field %d named as '%', instead '%s' was parsed",pszMaterialRNAGBSWGSFile,3,"WGSSampleName",pszHdrField);
			Reset();
			return(eBSFerrParse);
			}

		m_pRNAGBSWGSFile->GetText(4, &pszHdrField);		// expected to be MaterialName
		if(stricmp(pszHdrField,"MaterialName"))
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Input materials mapping file '%s' header does not have expected field %d named as '%', instead '%s' was parsed",pszMaterialRNAGBSWGSFile,1,"MaterialName",pszHdrField);
			Reset();
			return(eBSFerrParse);
			}
		continue;
		}

	if((NumFields = m_pRNAGBSWGSFile->GetCurFields()) != ExpNumFields)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Input materials mapping file '%s' expected to contain same number of fields as header line (%d), it contains %d at line %d", pszMaterialRNAGBSWGSFile,ExpNumFields, NumFields, CurLineNumber);
		Reset();
		return(eBSFerrParse);
		}

	if(m_NumRNASampleNames == cMaxRNAESamples)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Input materials mapping file '%s' contains more than the maximum number (%d)of unique RNA sample names which can be accepted at line %d", pszMaterialRNAGBSWGSFile,cMaxRNAESamples,CurLineNumber);
		Reset();
		return(eBSFerrParse);
		}
	m_pRNAGBSWGSFile->GetText(1, &pszValue);		// expected to be RNASampleName, these must be unique
	if ((RNASampleID = AddRNASampleName(pszValue)) == 0)
		{
		if(LocateRNASampleNameID(pszValue) != 0)
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Input materials mapping file '%s' unable to accept duplicate RNA sample name '%s' at line %d", pszMaterialRNAGBSWGSFile,pszValue,CurLineNumber);
		else
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Input materials mapping file '%s' unable to accept RNA sample name '%s' at line %d", pszMaterialRNAGBSWGSFile,pszValue,CurLineNumber);
		Reset();
		return(eBSFerrParse);
		}

	m_pRNAGBSWGSFile->GetText(3, &pszValue);			// expected to be WGSSampleName
	if((WGSSampleID = AddWGSSampleName(pszValue)) == 0) // can be multiple instances of WGS sample names
		WGSSampleID = LocateWGSSampleNameID(pszValue);
	if(WGSSampleID == 0)
		{
		if(m_NumWGSSampleNames == cMaxWGSSamples)
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Input materials mapping file '%s' contains more than the maximum number (%d)of unique WGS sample names which can be accepted at line %d", pszMaterialRNAGBSWGSFile,cMaxWGSSamples,CurLineNumber);
		else
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Input materials mapping file '%s' unable to accept WGS sample name '%s' at line %d", pszMaterialRNAGBSWGSFile,pszValue,CurLineNumber);
		Reset();
		return(eBSFerrParse);
		}
	
	m_pRNAGBSWGSFile->GetText(4, &pszValue);			// expected to be MaterialName
	if((MaterialID = AddMaterialName(pszValue)) == 0)	// can be multiple instances of same material
		MaterialID = LocateMaterialNameID(pszValue);
	if(MaterialID == 0)
		{
		if(m_NumMaterialNames == cMaxMaterials)
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Input materials mapping file '%s' contains more than the maximum number (%d)of unique Material names which can be accepted at line %d", pszMaterialRNAGBSWGSFile,cMaxMaterials,CurLineNumber);
		else
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Input materials mapping file '%s' unable to accept Material name '%s' at line %d", pszMaterialRNAGBSWGSFile,pszValue,CurLineNumber);
		Reset();
		return(eBSFerrParse);
		}
	m_RNAGBSWGSMaterials[RNASampleID-1].RNASampleID = RNASampleID;
	m_RNAGBSWGSMaterials[RNASampleID-1].WGSSampleID = WGSSampleID;
	m_RNAGBSWGSMaterials[RNASampleID-1].MaterialID = MaterialID;
	}

delete m_pRNAGBSWGSFile;
m_pRNAGBSWGSFile = nullptr;
return(NumbRNASamples);
}

double	// returned correlation coefficient
CRNAExpr::PearsonCorrelation(int32_t StartFeatureID,	// starting from this feature (row)
			int32_t EndFeatureID,	// ending with this inclusive feature (row)
		  int32_t xSampleID,		// x sample (column)
		  int32_t ySampleID)		// y sample (column)
{
double xMean;
double yMean;
double SumDiffs;
double xSumDiffs2;
double ySumDiffs2;
double xySumDiffs2;
double r;
int32_t *pxValue;
int32_t *pyValue;

int32_t RowIdx;

pxValue = (int32_t *)& m_pRNACntsMem[(((size_t)(StartFeatureID-1) * m_NumbHdrRNAMappings) + xSampleID-1) * sizeof(int32_t)];
pyValue = (int32_t *)& m_pRNACntsMem[(((size_t)(StartFeatureID-1) * m_NumbHdrRNAMappings) + ySampleID-1) * sizeof(int32_t)];
xMean = 0.0;
yMean = 0.0;
for(RowIdx = StartFeatureID-1; RowIdx < EndFeatureID; RowIdx++,pxValue+=m_NumbHdrRNAMappings,pyValue+=m_NumbHdrRNAMappings)
	{
	xMean += *pxValue;
	yMean += *pyValue;
	}
xMean /= 1 + EndFeatureID - StartFeatureID;
yMean /= 1 + EndFeatureID - StartFeatureID;

pxValue = (int32_t *)& m_pRNACntsMem[(((size_t)(StartFeatureID-1) * m_NumbHdrRNAMappings) + xSampleID-1) * sizeof(int32_t)];
pyValue = (int32_t *)& m_pRNACntsMem[(((size_t)(StartFeatureID-1) * m_NumbHdrRNAMappings) + ySampleID-1) * sizeof(int32_t)];
SumDiffs = 0.0;
xSumDiffs2 = 0.0;
ySumDiffs2 = 0.0;
for(RowIdx = StartFeatureID-1; RowIdx < EndFeatureID; RowIdx++,pxValue+=m_NumbHdrRNAMappings,pyValue+=m_NumbHdrRNAMappings)
	{
	SumDiffs += (*pxValue - xMean) * (*pyValue - yMean);
	xSumDiffs2 +=  (*pxValue - xMean) * (*pxValue - xMean);
	ySumDiffs2 +=  (*pyValue - yMean) * (*pyValue - yMean);
	}
xySumDiffs2 = sqrt(xSumDiffs2 * ySumDiffs2);
r = SumDiffs/xySumDiffs2;
return(r);
}

double	// returned correlation coefficient
CRNAExpr::PearsonCorrelation(int32_t RNArowIdx,	// correlate this RNA row
							 int32_t WGSrowIdx)	// with this WGS row
{
double xMean;
double yMean;
double SumDiffs;
double xSumDiffs2;
double ySumDiffs2;
double xySumDiffs2;
double r;
double *pRNA;
double *pWGS;

int32_t ScoreIdx;

pRNA = (double *)& m_pRNAhomozScoresMem[((size_t)RNArowIdx * m_NumbHdrWGSMappings * sizeof(double))];
pWGS = (double *)& m_pWGShomozScoresMem[((size_t)WGSrowIdx * m_NumbHdrWGSMappings * sizeof(double))];
xMean = 0.0;
yMean = 0.0;
for(ScoreIdx = 0; ScoreIdx < m_NumbHdrWGSMappings; ScoreIdx++,pRNA++,pWGS++)
	{
	xMean += *pRNA;
	yMean += *pWGS;
	}
xMean /= m_NumbHdrWGSMappings;
yMean /= m_NumbHdrWGSMappings;

pRNA = (double *)& m_pRNAhomozScoresMem[((size_t)RNArowIdx * m_NumbHdrWGSMappings * sizeof(double))];
pWGS = (double *)& m_pWGShomozScoresMem[((size_t)WGSrowIdx * m_NumbHdrWGSMappings * sizeof(double))];
SumDiffs = 0.0;
xSumDiffs2 = 0.0;
ySumDiffs2 = 0.0;
for(ScoreIdx = 0; ScoreIdx < m_NumbHdrWGSMappings; ScoreIdx++,pRNA++,pWGS++)
	{
	SumDiffs += (*pRNA - xMean) * (*pWGS - yMean);
	xSumDiffs2 +=  (*pRNA - xMean) * (*pRNA - xMean);
	ySumDiffs2 +=  (*pWGS - yMean) * (*pWGS - yMean);
	}
xySumDiffs2 = sqrt(xSumDiffs2 * ySumDiffs2);
r = SumDiffs/xySumDiffs2;
return(r);
}

int	// checking for expression level counts inconsistencies - if replicate 1 and replicate 2 from same biological sample then expecting expression level counts to track over the features
CRNAExpr::GenExprCntsPearsons(eModeRNAE PMode,			// processing mode
					char* pszMaterialRNAGBSWGSFile,		// load RNA/GBS/WIG sample names and RNA associated metadata
					char* pszInCntsFile,					// load coverage counts from this file
					char* pszOutRsltsFile)				// write results to this file, will be suffixed appropriately
{
int32_t SampleID;
int32_t BuffIdx;
char *pszSample1;
char *pszSample2;
char *pszMaterial1;
char *pszMaterial2;
char szOutBuff[4000];
char szOutFile[_MAX_PATH];
FILE *pOutStream;

double r;
sprintf(szOutFile,"%s%s",pszOutRsltsFile,".csv");
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Generating Pearson correlation coefficient 'r' for labeled biological replicates over all features into file '%s'",szOutFile);

if((pOutStream = fopen(szOutFile,"w"))==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create/truncate file %s for writing error: %s",szOutFile,strerror(errno));
	Reset();
	return(eBSFerrOpnFile);
	}

BuffIdx = sprintf(szOutBuff,"\"SampleRep1\",\"SampleRep2\",\"Pearson\",\"Material1\",\"Material2\"");
for(SampleID = 1; SampleID <= (int32_t)m_NumbHdrRNAMappings; SampleID+=2)
	{
	r = PearsonCorrelation(1,m_NumFeatureNames,SampleID,SampleID+1);
	pszSample1 = LocateRNASampleName(m_HdrRNAMappings[SampleID-1]);
	pszSample2 = LocateRNASampleName(m_HdrRNAMappings[SampleID]);
	pszMaterial1 = LocateMaterialName(m_RNAGBSWGSMaterials[m_HdrRNAMappings[SampleID-1]-1].MaterialID);
	pszMaterial2 = LocateMaterialName(m_RNAGBSWGSMaterials[m_HdrRNAMappings[SampleID]-1].MaterialID);
	BuffIdx += sprintf(&szOutBuff[BuffIdx],"\n\"%s\",\"%s\",%.5f,\"%s\",\"%s\"",pszSample1,pszSample2,r,pszMaterial1,pszMaterial2);
	if(BuffIdx + 200 > sizeof(szOutBuff))
		{
		fwrite(szOutBuff,1,BuffIdx,pOutStream);
		BuffIdx = 0;
		}
	}
if(BuffIdx)
	fwrite(szOutBuff,1,BuffIdx,pOutStream);
fflush(pOutStream);
fclose(pOutStream);

// let's try finding the best pairings assuming that the low r's in some pairs are due to a mixup in replicate labeling
sprintf(szOutFile,"%s%s",pszOutRsltsFile,"recoverreps.csv");
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Generating Pearson correlation coefficient 'r' targeting discovery of potentially mislabeled replicates over all features into file '%s'",szOutFile);
if ((pOutStream = fopen(szOutFile, "w")) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create/truncate file %s for writing error: %s",szOutFile,strerror(errno));
	Reset();
	return(eBSFerrOpnFile);
	}

double Maximalr;
int32_t ChkSampleID;
int32_t ClosestSampleID;
BuffIdx = sprintf(szOutBuff,"\"SampleRep1\",\"ClosestRep2\",\"Pearson\",\"Material1\",\"Material2\"");
for(SampleID = 1; SampleID <= (int32_t)m_NumbHdrRNAMappings; SampleID+=2)
	{
	Maximalr = 0.0;
	ClosestSampleID = 0;
	for(ChkSampleID = 1; ChkSampleID <= (int32_t)m_NumbHdrRNAMappings; ChkSampleID++)
		{
		if(ChkSampleID == SampleID)
			continue;
		r = PearsonCorrelation(1,m_NumFeatureNames,SampleID,ChkSampleID);
		if(r > Maximalr)
			{
			ClosestSampleID = ChkSampleID;
			Maximalr = r;
			}
		}
	pszSample1 = LocateRNASampleName(m_HdrRNAMappings[SampleID-1]);
	pszSample2 = LocateRNASampleName(m_HdrRNAMappings[ClosestSampleID-1]);
	pszMaterial1 = LocateMaterialName(m_RNAGBSWGSMaterials[m_HdrRNAMappings[SampleID-1]-1].MaterialID);
	pszMaterial2 = LocateMaterialName(m_RNAGBSWGSMaterials[m_HdrRNAMappings[ClosestSampleID-1]-1].MaterialID);
	BuffIdx += sprintf(&szOutBuff[BuffIdx],"\n\"%s\",\"%s\",%.5f,\"%s\",\"%s\"",pszSample1,pszSample2,Maximalr,pszMaterial1,pszMaterial2);
	if(BuffIdx + 200 > sizeof(szOutBuff))
		{
		fwrite(szOutBuff,1,BuffIdx,pOutStream);
		BuffIdx = 0;
		}
	}
if(BuffIdx)
	fwrite(szOutBuff,1,BuffIdx,pOutStream);
fflush(pOutStream);
fclose(pOutStream);
return(eBSFSuccess);
}

int
CRNAExpr::LoadWGSScoresFile(char* pszInWGSvWGSFile) // parse and load WGS vs. WGS homozygosity scores
{
int32_t Rslt;
uint32_t EstNumRows;
int64_t FileSize;
int32_t MaxFields;
int32_t MeanNumFields;
int FieldIdx;
int32_t CurLineNumber;
int32_t NumFields;
int32_t ExpNumFields;
int32_t SampleNameID;
char *pszWGSChrom;
char *pszWGSSample;
double ScoreValue;

int32_t FeatureNameID;

if(m_pInWGSvWGSFile != nullptr) // shouldn't have been instantiated, but better to be sure!
	{
	delete m_pInWGSvWGSFile;
	m_pInWGSvWGSFile = nullptr;
	}
if((m_pInWGSvWGSFile = new CCSVFile) == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to instantiate CCSVfile");
	Reset();
	return(eBSFerrObj);
	}

	// get an estimate of number of rows and fields
if ((EstNumRows = m_pInWGSvWGSFile->CSVEstSizes(pszInWGSvWGSFile, &FileSize, &MaxFields, &MeanNumFields)) < 2)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to estimate number of rows in file: '%s'", pszInWGSvWGSFile);
	Reset();
	return(eBSFerrFieldCnt);
	}

m_pInWGSvWGSFile->SetMaxFields(MaxFields);

if((Rslt = m_pInWGSvWGSFile->Open(pszInWGSvWGSFile)) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to open WGS vs. WGS scores file: '%s'", pszInWGSvWGSFile);
	Reset();
	return(Rslt);
	}

// header line contains sample identifiers
// rows contain features associated with samples
CurLineNumber = 0;
ExpNumFields = 0;
FeatureNameID = 0;
m_NumbRowWGSMappings = 0;
m_NumbRowWGSMappings = 0;
while((Rslt = m_pInWGSvWGSFile->NextLine()) > 0)		// onto next line containing fields
	{
	CurLineNumber++;
	if((NumFields = m_pInWGSvWGSFile->GetCurFields()) < 3)	// must contain at least 3 fields (assumes scores for at least 1 sample!)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Input WGS vs. WGS scores file '%s' expected to contain a minimum of 3 fields, it contains %d at line %d", pszInWGSvWGSFile, NumFields, CurLineNumber);
		Reset();
		return(eBSFerrParse);
		}

	if (CurLineNumber == 1) // 1 if header containing sample identifiers
		{
		if (!m_pInWGSvWGSFile->IsLikelyHeaderLine())
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Input WGS vs. WGS scores file '%s' line 1 does not parse as a header line, all fields are expected to be double quoted as containing string values", pszInWGSvWGSFile);
			Reset();
			return(eBSFerrParse);
			}

		ExpNumFields = NumFields;
		// parse out WGS sample identifiers until last column
		for(FieldIdx = 3; FieldIdx <= NumFields; FieldIdx++)
			{
			m_pInWGSvWGSFile->GetText(FieldIdx, &pszWGSSample);
			if ((SampleNameID = AddWGSSampleName(pszWGSSample)) != 0) // WGS sample name should already be known from the material mapping file, but this is not assured!
				gDiagnostics.DiagOut(eDLWarn, gszProcName, "Input WGS vs. WGS scores file '%s' contains a WGS sample name '%s' not in materials file at field %d at line %d", pszInWGSvWGSFile, pszWGSSample, FieldIdx, CurLineNumber);
			SampleNameID = LocateWGSSampleNameID(pszWGSSample);
			m_HdrWGSMappings[m_NumbHdrWGSMappings++] = SampleNameID;
			}
		
		if ((Rslt = AllocateWGShomozScores(m_NumbHdrWGSMappings * m_NumbHdrWGSMappings)) != eBSFSuccess) // score matrix is expected to be orthogonal
			{
			Reset();
			return(Rslt);
			}
		continue;
		}
	if((NumFields = m_pInWGSvWGSFile->GetCurFields()) != ExpNumFields)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Input WGS vs. WGS scores file '%s' expected to contain same number of fields as header line (%d), it contains %d at line %d", pszInWGSvWGSFile,ExpNumFields, NumFields, CurLineNumber);
		Reset();
		return(eBSFerrParse);
		}

	// each row contains the WGS which was scored against all other WGSs
	m_pInWGSvWGSFile->GetText(1, &pszWGSChrom);	// currently individual chromosome scores are skipped, interest is only in the complete assembly scores

	if(stricmp(pszWGSChrom,"AllChroms"))
		if(stricmp(pszWGSChrom,"AllChroms\"")) 	// early releases of callhaplotypes had a double quote bug when generating the AllChroms summary counts/scores
		continue;

	m_pInWGSvWGSFile->GetText(2, &pszWGSSample);		// currently individual chromosome scores are skipped, interest is only in the complete assembly scores
	if ((SampleNameID = LocateWGSSampleNameID(pszWGSSample)) == 0) // WGS sample name must already be known from the headings!
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Input WGS vs. WGS scores file '%s' contains a WGS sample name '%s' not already known from header at field 2 at line %d", pszInWGSvWGSFile, pszWGSSample, CurLineNumber);
		Reset();
		return(eBSFerrParse);
		}
	if(m_HdrWGSMappings[m_NumbRowWGSMappings] != SampleNameID)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Input WGS vs. WGS scores file '%s' contains a row WGS sample name '%s' not orthogonal to header WGS sample name at field 2 at line %d", pszInWGSvWGSFile, pszWGSSample, CurLineNumber);
		Reset();
		return(eBSFerrParse);
		}
	m_RowWGSMappings[m_NumbRowWGSMappings++] = SampleNameID;
	if((Rslt=ReallocWGShomozScores(m_NumbHdrWGSMappings))!=eBSFSuccess)	// ensure sufficient memory preallocd to hold at least one additional row of WGS scores
		{
		Reset();
		return(Rslt);
		}
	double *pScore;
	pScore = (double *)&m_pWGShomozScoresMem[((size_t)(m_NumbRowWGSMappings-1) * m_NumbHdrWGSMappings) * sizeof(double)];
	for (FieldIdx = 3; FieldIdx <= ExpNumFields; FieldIdx++,pScore++)
		{
		m_pInWGSvWGSFile->GetDouble(FieldIdx, &ScoreValue);
		*pScore = ScoreValue;
		}
	}
delete m_pInWGSvWGSFile;
m_pInWGSvWGSFile = nullptr;

if (m_NumbRowWGSMappings != m_NumbHdrWGSMappings)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Input WGS vs. WGS scores file '%s' must contain orthogonal matrix scores", pszInWGSvWGSFile);
	Reset();
	return(eBSFerrParse);
	}

return((m_NumbHdrWGSMappings > 0 && m_NumbRowWGSMappings == m_NumbHdrWGSMappings) ? eBSFSuccess : eBSFerrParse);
}

int
CRNAExpr::LoadRNAScoresFile(char* pszInRNAvWGSFile)		// parse and load RNA vs. WGS homozygosity scores
{
int32_t Rslt;
uint32_t EstNumRows;
int64_t FileSize;
int32_t MaxFields;
int32_t MeanNumFields;
int FieldIdx;
int32_t CurLineNumber;
int32_t NumFields;
int32_t ExpNumFields;
int32_t SampleNameID;
char *pszWGSSample;
char *pszRNAChrom;
char *pszRNASample;
double ScoreValue;

if(m_pInRNAvWGSFile != nullptr) // shouldn't have been instantiated, but better to be sure!
	{
	delete m_pInRNAvWGSFile;
	m_pInRNAvWGSFile = nullptr;
	}
if((m_pInRNAvWGSFile = new CCSVFile) == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to instantiate CCSVfile");
	Reset();
	return(eBSFerrObj);
	}

	// get an estimate of number of rows and fields
if ((EstNumRows = m_pInRNAvWGSFile->CSVEstSizes(pszInRNAvWGSFile, &FileSize, &MaxFields, &MeanNumFields)) < 2)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to estimate number of rows in file: '%s'", pszInRNAvWGSFile);
	Reset();
	return(eBSFerrFieldCnt);
	}

m_pInRNAvWGSFile->SetMaxFields(MaxFields);

if((Rslt = m_pInRNAvWGSFile->Open(pszInRNAvWGSFile)) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to open RNA vs. WGS scores file: '%s'", pszInRNAvWGSFile);
	Reset();
	return(Rslt);
	}

// header line contains WGS sample names
// the header WGS sample names must be in exactly the same order as those in the WGS vs. WGS score file
// rows contain RNA samples and scores for each WGS
CurLineNumber = 0;
ExpNumFields = 0;
m_NumbRowRNAMappings = 0;
while((Rslt = m_pInRNAvWGSFile->NextLine()) > 0)		// onto next line containing fields
	{
	CurLineNumber++;
	if((NumFields = m_pInRNAvWGSFile->GetCurFields()) < 3)	// must contain at least 3 fields (assumes scores for at least 1 sample!)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Input RNA vs. WGS scores file '%s' expected to contain a minimum of 3 fields, it contains %d at line %d", pszInRNAvWGSFile, NumFields, CurLineNumber);
		Reset();
		return(eBSFerrParse);
		}

	if (CurLineNumber == 1) // 1 if header containing sample identifiers
		{
		if (!m_pInRNAvWGSFile->IsLikelyHeaderLine())
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Input RNA vs. WGS scores file '%s' line 1 does not parse as a header line, all fields are expected to be double quoted as containing string values", pszInRNAvWGSFile);
			Reset();
			return(eBSFerrParse);
			}

		if (NumFields - 2 != m_NumbHdrWGSMappings)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Input RNA vs. WGS scores file '%s' header line must contain same number of WGS vs. WGS sample names and in same order", pszInRNAvWGSFile);
			Reset();
			return(eBSFerrParse);
			}

		ExpNumFields = NumFields;
		// parse out WGS sample identifiers until last column
		for(FieldIdx = 3; FieldIdx <= NumFields; FieldIdx++)
			{
			m_pInRNAvWGSFile->GetText(FieldIdx, &pszWGSSample);
			if ((SampleNameID = AddWGSSampleName(pszWGSSample)) != 0) // WGS sample name must already be known from the WGS vs. WGS file
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Input RNA vs. WGS scores file '%s' contains a header WGS sample name '%s' not in WGS vs. WGS file at field %d at line %d", pszInRNAvWGSFile, pszWGSSample, FieldIdx, CurLineNumber);
				Reset();
				return(eBSFerrParse);
				}
			SampleNameID = LocateWGSSampleNameID(pszWGSSample);
			if (m_HdrWGSMappings[FieldIdx-3] != SampleNameID)
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Input RNA vs. WGS scores file '%s' header line must contain same ordering as WGS vs. WGS sample name '%s'", pszInRNAvWGSFile,pszWGSSample);
				Reset();
				return(eBSFerrParse);
				}
			}
		if ((Rslt = AllocateRNAhomozScores(m_NumbHdrWGSMappings * m_NumbHdrRNAMappings)) != eBSFSuccess)
			{
			Reset();
			return(Rslt);
			}
		continue;
		}
	if((NumFields = m_pInRNAvWGSFile->GetCurFields()) != ExpNumFields)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Input RNA vs. WGS scores file '%s' expected to contain same number of fields as header line (%d), it contains %d at line %d", pszInRNAvWGSFile,ExpNumFields, NumFields, CurLineNumber);
		Reset();
		return(eBSFerrParse);
		}

	// each row contains the RNA sample which was scored against all WGSs
	m_pInRNAvWGSFile->GetText(1, &pszRNAChrom);	// currently individual chromosome scores are skipped, interest is only in the complete assembly scores

	if(stricmp(pszRNAChrom,"AllChroms"))
		if(stricmp(pszRNAChrom,"AllChroms\"")) 	// early releases of callhaplotypes had a double quote bug when generating the AllChroms summary counts/scores
		continue;

	m_pInRNAvWGSFile->GetText(2, &pszRNASample);		// currently individual chromosome scores are skipped, interest is only in the complete assembly scores
	if ((SampleNameID = LocateRNASampleNameID(pszRNASample)) == 0) // RNA sample name must already be known from the material mappings
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Input RNA vs. WGS scores file '%s' contains a RNA sample name '%s' not already known from material mappings at field 2 at line %d", pszInRNAvWGSFile, pszRNASample, CurLineNumber);
		Reset();
		return(eBSFerrParse);
		}
	m_RowRNAMappings[m_NumbRowRNAMappings++] = SampleNameID;
	if((Rslt=ReallocRNAhomozScores(m_NumbHdrWGSMappings))!=eBSFSuccess)	// ensure sufficient memory preallocd to hold at least one additional row of RNA vs.WGS scores
		{
		Reset();
		return(Rslt);
		}
		double *pScore;
	pScore = (double *)&m_pRNAhomozScoresMem[((size_t)(m_NumbRowRNAMappings-1) * m_NumbHdrWGSMappings) * sizeof(double)];
	for (FieldIdx = 3; FieldIdx <= ExpNumFields; FieldIdx++,pScore++)
		{
		m_pInRNAvWGSFile->GetDouble(FieldIdx, &ScoreValue);
		*pScore = ScoreValue;
		}
	}
delete m_pInRNAvWGSFile;
m_pInRNAvWGSFile = nullptr;

return((m_NumbHdrWGSMappings > 0 && m_NumbRowRNAMappings > 0) ? eBSFSuccess : eBSFerrParse);
}

int
CRNAExpr::Process(eModeRNAE PMode,			// processing mode
				char* pszMaterialRNAGBSWGSFile, // load RNA/GBS/WIG sample names and RNA associated metadata
				char* pszInCntsFile,		// load coverage counts from this file
				char* pszInWGSvWGSFile,		// load WGS vs. WGS homozygosity scores from this file
				char* pszInRNAvWGSFile,		// load RNA vs. WGS homozygosity scores from this file
				char* pszOutRsltsFile)			// write results to this file, will be suffixed appropriately
{
int Rslt;
int32_t BuffIdx;
char szOutBuff[4000];
char szOutFile[_MAX_PATH];
FILE *pOutStream;
Reset();

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading RNA material mappings from '%s'",pszMaterialRNAGBSWGSFile);
Rslt = LoadMaterialRNAGBSWGSFile(pszMaterialRNAGBSWGSFile);
if(Rslt < eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to load material mappings from '%s'",pszMaterialRNAGBSWGSFile);
	Reset();
	return(Rslt);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading RNA raw counts from '%s'",pszInCntsFile);
Rslt = LoadRNACntsFile(pszInCntsFile,true);
if(Rslt < eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to load RNA raw counts from %s",pszInCntsFile);
	Reset();
	return(Rslt);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"RNA raw counts loaded from '%s'",pszInCntsFile);

// checking for expression level counts inconsistencies - if replicate 1 and replicate 2 from same biological sample then expecting expression level counts to track over the features
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Raw counts normalised and now checking for replicate expression inconsistencies");
Rslt = GenExprCntsPearsons(PMode,pszMaterialRNAGBSWGSFile,pszInCntsFile,pszOutRsltsFile);
if(Rslt != eBSFSuccess || PMode == eMRNAEDefault)
	return(Rslt);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading WGS vs. WGS homozygosity scores from '%s'",pszInWGSvWGSFile);
Rslt = LoadWGSScoresFile(pszInWGSvWGSFile);
if(Rslt < eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to load WGS homozygosity scores from %s",pszInWGSvWGSFile);
	Reset();
	return(Rslt);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"WGS vs. WGS homozygosity scores loaded from %s",pszInWGSvWGSFile);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading RNA vs. WGS homozygosity scores from '%s'",pszInRNAvWGSFile);
Rslt = LoadRNAScoresFile(pszInRNAvWGSFile);
if(Rslt < eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to load RNA vs. WGS homozygosity scores from %s",pszInRNAvWGSFile);
	Reset();
	return(Rslt);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"RNA vs. WGS homozygosity scores loaded from %s",pszInRNAvWGSFile);

double r;
sprintf(szOutFile,"%s%s",pszOutRsltsFile,".RNAvWGShomozygosity.csv");
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Generating Pearson correlation coefficient 'r' for WGS vs. WGS homozygosity into file '%s'",szOutFile);

if((pOutStream = fopen(szOutFile,"w"))==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create/truncate file %s for writing error: %s",szOutFile,strerror(errno));
	Reset();
	return(eBSFerrOpnFile);
	}
	

BuffIdx = sprintf(szOutBuff,"\"RNA Rep\",\"RNA Rep Material\",\"Match WGS\",\"Match WGS Material\",\"Match Pearson\"");
int32_t RNArowIdx;
int32_t WGSrowIdx;
int32_t MaxWGSrowIdx = 0;
double MaxPearson = 0.0;
for (RNArowIdx = 0; RNArowIdx < (int32_t)m_NumbRowRNAMappings; RNArowIdx++)
	{
	MaxWGSrowIdx = 0;
	MaxPearson = 0.0;
	for(WGSrowIdx = 0; WGSrowIdx < (int32_t)m_NumbRowWGSMappings; WGSrowIdx++)
		{
		r = PearsonCorrelation(RNArowIdx,WGSrowIdx);
		if(r > MaxPearson)
			{
			MaxPearson = r;
			MaxWGSrowIdx = WGSrowIdx;
			}
		}
	char *pWGSname = LocateWGSSampleName(m_RowWGSMappings[MaxWGSrowIdx]);
	char *pRNAname = LocateRNASampleName(m_RowRNAMappings[RNArowIdx]);
	char *pszWGSmaterial;
	char *pszRNAmaterial;
	int RepRNAWGSMaterialSame = 0;
	int PairRNAWGSMaterialSame = 0;
	int PairWGSMaterialMatch = 0;

	pszRNAmaterial = LocateMaterialName(m_RNAGBSWGSMaterials[m_HdrRNAMappings[RNArowIdx]-1].MaterialID);
	pszWGSmaterial = LocateWGSMaterialName(m_RowWGSMappings[MaxWGSrowIdx]);
	if(pszWGSmaterial == nullptr)
		pszWGSmaterial = (char *)"#N/A";

	BuffIdx += sprintf(&szOutBuff[BuffIdx],"\n\"%s\",\"%s\",\"%s\",\"%s\",%.6f",
						pRNAname,pszRNAmaterial,pWGSname,pszWGSmaterial,MaxPearson);
	if((BuffIdx + 500) > sizeof(szOutBuff))
		{
		fwrite(szOutBuff,1,BuffIdx,pOutStream);
		BuffIdx = 0;
		}
	}
if(BuffIdx)
	fwrite(szOutBuff,1,BuffIdx,pOutStream);
fflush(pOutStream);
fclose(pOutStream);

Reset();
return(Rslt);
}





// NOTE: SampleNames are checked for uniqueness as SampleNames must be unique
uint32_t		// returned SampleName identifier, 0 if unable to accept this SampleName name
CRNAExpr::AddRNASampleName(char* pszSampleName) // associate unique identifier with this SampleName name
{
uint32_t SampleNameIdx;
int SampleNameLen;
char* pszLAname;

// with any luck the sequence name will be same as the last accessed
if ((pszLAname = LocateRNASampleName(m_LARNASampleNameID)) != nullptr)
	{
	if (!stricmp(pszSampleName, pszLAname))
		return(0); // non-unique!
	}

// iterate over all known SampleNames in case this SampleName to add is a duplicate
for (SampleNameIdx = 0; SampleNameIdx < m_NumRNASampleNames; SampleNameIdx++)
	{
	pszLAname = &m_szRNASampleNames[m_szRNASampleNamesIdx[SampleNameIdx]];
	if (!stricmp(pszSampleName, pszLAname))
		{
		m_LARNASampleNameID = SampleNameIdx + 1;
		return(0); // non-unique!
		}
	}

// not a duplicate
SampleNameLen = (int)strlen(pszSampleName);
if ((m_NxtszRNASampleNameIdx + SampleNameLen + 1) > (int)sizeof(m_szRNASampleNames))
	return(0);		// no more space to hold reeadset name, treating as if a non-unique!
if (m_NumRNASampleNames == cMaxRNAESamples)
	return(0);		// unable to hold any more SampleNames, treating as if a non-unique!

m_szRNASampleNamesIdx[m_NumRNASampleNames++] = m_NxtszRNASampleNameIdx;
strcpy(&m_szRNASampleNames[m_NxtszRNASampleNameIdx], pszSampleName);
m_NxtszRNASampleNameIdx += SampleNameLen + 1;
m_LARNASampleNameID = m_NumRNASampleNames;
return(m_LARNASampleNameID);
}

uint32_t		// returned SampleName identifier, 0 if unable to locate this SampleName name
CRNAExpr::LocateRNASampleNameID(char* pszSampleName) // return unique identifier associated with this SampleName name
{
uint32_t SampleNameIdx;
char* pszLASampleName;

// with any luck the SampleName name will be same as the last accessed
if ((pszLASampleName = LocateRNASampleName(m_LARNASampleNameID)) != nullptr)
	if (!stricmp(pszSampleName, pszLASampleName))
		return(m_LARNASampleNameID);

	// iterate over all known SampleNames
for (SampleNameIdx = 0; SampleNameIdx < m_NumRNASampleNames; SampleNameIdx++)
	{
	pszLASampleName = &m_szRNASampleNames[m_szRNASampleNamesIdx[SampleNameIdx]];
	if (!stricmp(pszSampleName, pszLASampleName))
		{
		m_LARNASampleNameID = SampleNameIdx + 1;
		return(m_LARNASampleNameID);
		}
	}
return(0);
}

char*
CRNAExpr::LocateRNASampleName(uint32_t SampleNameID)
{
uint32_t Idx;
Idx = SampleNameID & 0x0fffffff;			// mask out any potential flag bits
if (Idx < 1 || Idx > m_NumRNASampleNames)
	return(nullptr);
return(&(m_szRNASampleNames[m_szRNASampleNamesIdx[Idx - 1]]));
}

// NOTE: WGSSampleNames are checked for uniqueness as SampleNames must be unique
uint32_t		// returned SampleName identifier, 0 if unable to accept this SampleName name
CRNAExpr::AddWGSSampleName(char* pszSampleName) // associate unique identifier with this WGS SampleName name
{
uint32_t SampleNameIdx;
int SampleNameLen;
char* pszLAname;

// with any luck the sequence name will be same as the last accessed
if ((pszLAname = LocateWGSSampleName(m_LAWGSSampleNameID)) != nullptr)
	{
	if (!stricmp(pszSampleName, pszLAname))
		return(0); // non-unique!
	}

// iterate over all known SampleNames in case this SampleName to add is a duplicate
for (SampleNameIdx = 0; SampleNameIdx < m_NumWGSSampleNames; SampleNameIdx++)
	{
	pszLAname = &m_szWGSSampleNames[m_szWGSSampleNamesIdx[SampleNameIdx]];
	if (!stricmp(pszSampleName, pszLAname))
		{
		m_LAWGSSampleNameID = SampleNameIdx + 1;
		return(0); // non-unique!
		}
	}

// not a duplicate
SampleNameLen = (int)strlen(pszSampleName);
if ((m_NxtszWGSSampleNameIdx + SampleNameLen + 1) > (int)sizeof(m_szWGSSampleNames))
	return(0);		// no more space to hold readset name, treating as if a non-unique!
if (m_NumWGSSampleNames == cMaxWGSSamples)
	return(0);		// unable to hold any more SampleNames, treating as if a non-unique!

m_szWGSSampleNamesIdx[m_NumWGSSampleNames++] = m_NxtszWGSSampleNameIdx;
strcpy(&m_szWGSSampleNames[m_NxtszWGSSampleNameIdx], pszSampleName);
m_NxtszWGSSampleNameIdx += SampleNameLen + 1;
m_LAWGSSampleNameID = m_NumWGSSampleNames;
return(m_LAWGSSampleNameID);
}

uint32_t		// returned SampleName identifier, 0 if unable to locate this SampleName name
CRNAExpr::LocateWGSSampleNameID(char* pszSampleName) // return unique identifier associated with this SampleName name
{
uint32_t SampleNameIdx;
char* pszLASampleName;

// with any luck the SampleName name will be same as the last accessed
if ((pszLASampleName = LocateWGSSampleName(m_LAWGSSampleNameID)) != nullptr)
	if (!stricmp(pszSampleName, pszLASampleName))
		return(m_LAWGSSampleNameID);

	// iterate over all known SampleNames
for (SampleNameIdx = 0; SampleNameIdx < m_NumWGSSampleNames; SampleNameIdx++)
	{
	pszLASampleName = &m_szWGSSampleNames[m_szWGSSampleNamesIdx[SampleNameIdx]];
	if (!stricmp(pszSampleName, pszLASampleName))
		{
		m_LAWGSSampleNameID = SampleNameIdx + 1;
		return(m_LAWGSSampleNameID);
		}
	}
return(0);
}

char*
CRNAExpr::LocateWGSSampleName(uint32_t SampleNameID)
{
uint32_t Idx;
Idx = SampleNameID & 0x0fffffff;			// mask out any potential flag bits
if (Idx < 1 || Idx > m_NumWGSSampleNames)
	return(nullptr);
return(&(m_szWGSSampleNames[m_szWGSSampleNamesIdx[Idx - 1]]));
}


// NOTE: MaterialNames are checked for uniqueness as MaterialNames must be unique
uint32_t		// returned MaterialName identifier, 0 if unable to accept this MaterialName name
CRNAExpr::AddMaterialName(char* pszMaterialName) // associate unique identifier with this MaterialName name
{
uint32_t MaterialNameIdx;
int MaterialNameLen;
char* pszLAname;

// with any luck the material name will be same as the last accessed
if ((pszLAname = LocateMaterialName(m_LAMaterialNameID)) != nullptr)
	{
	if (!stricmp(pszMaterialName, pszLAname))
		return(0); // non-unique!
	}

// iterate over all known material names in case this MaterialName to add is a duplicate
for (MaterialNameIdx = 0; MaterialNameIdx < m_NumMaterialNames; MaterialNameIdx++)
	{
	pszLAname = &m_szMaterialNames[m_szMaterialNamesIdx[MaterialNameIdx]];
	if (!stricmp(pszMaterialName, pszLAname))
		{
		m_LAMaterialNameID = MaterialNameIdx + 1;
		return(0); // non-unique!
		}
	}

// not a duplicate
MaterialNameLen = (int)strlen(pszMaterialName);
if ((m_NxtszMaterialNameIdx + MaterialNameLen + 1) > (int)sizeof(m_szMaterialNames))
	return(0);		// no more space to hold reeadset name, treating as if a non-unique!
if (m_NumMaterialNames == cMaxMaterials)
	return(0);		// unable to hold any more MaterialNames, treating as if a non-unique!

m_szMaterialNamesIdx[m_NumMaterialNames++] = m_NxtszMaterialNameIdx;
strcpy(&m_szMaterialNames[m_NxtszMaterialNameIdx], pszMaterialName);
m_NxtszMaterialNameIdx += MaterialNameLen + 1;
m_LAMaterialNameID = m_NumMaterialNames;
return(m_LAMaterialNameID);
}

uint32_t		// returned MaterialName identifier, 0 if unable to locate this MaterialName name
CRNAExpr::LocateMaterialNameID(char* pszMaterialName) // return unique identifier associated with this MaterialName name
{
uint32_t MaterialNameIdx;
char* pszLAMaterialName;

// with any luck the MaterialName name will be same as the last accessed
if ((pszLAMaterialName = LocateMaterialName(m_LAMaterialNameID)) != nullptr)
	if (!stricmp(pszMaterialName, pszLAMaterialName))
		return(m_LAMaterialNameID);

	// iterate over all known MaterialNames
for (MaterialNameIdx = 0; MaterialNameIdx < m_NumMaterialNames; MaterialNameIdx++)
	{
	pszLAMaterialName = &m_szMaterialNames[m_szMaterialNamesIdx[MaterialNameIdx]];
	if (!stricmp(pszMaterialName, pszLAMaterialName))
		{
		m_LAMaterialNameID = MaterialNameIdx + 1;
		return(m_LAMaterialNameID);
		}
	}
return(0);
}

char*
CRNAExpr::LocateMaterialName(uint32_t MaterialNameID)
{
uint32_t Idx;
Idx = MaterialNameID & 0x0fffffff;			// mask out any potential flag bits
if (Idx < 1 || Idx > m_NumMaterialNames)
	return(nullptr);
return(&(m_szMaterialNames[m_szMaterialNamesIdx[Idx - 1]]));
}

char*	// returns material name corresponding to WGS identifier
CRNAExpr::LocateWGSMaterialName(uint32_t WGSNameID) // WGS identifier
{
tsRNAGBSWGSMaterial *pMaterial;
// need to iterate over all materials until a match on WGS identifier
pMaterial = m_RNAGBSWGSMaterials;
for(uint32_t Idx = 0; Idx < m_NumRNASampleNames; Idx++,pMaterial++)
	if(pMaterial->WGSSampleID == WGSNameID)
		return(LocateMaterialName(pMaterial->MaterialID));
return(nullptr);
}

// NOTE: FeatureNames are checked for uniqueness as FeatureNames must be unique
uint32_t		// returned FeatureName identifier, 0 if unable to accept this FeatureName
CRNAExpr::AddFeatureName(char* pszFeatureName) // associate unique identifier with this FeatureName
{
int32_t FeatureNameIdx;
int FeatureNameLen;
char* pszLAname;

// with any luck the sequence name will be same as the last accessed
if ((pszLAname = LocateFeatureName(m_LAFeatureNameID)) != nullptr)
	{
	if (!stricmp(pszFeatureName, pszLAname))
		return(0); // non-unique!
	}

// iterate over all known FeatureNames in case this FeatureName to add is a duplicate
for (FeatureNameIdx = 0; FeatureNameIdx < m_NumFeatureNames; FeatureNameIdx++)
	{
	pszLAname = &m_szFeatureNames[m_szFeatureNamesIdx[FeatureNameIdx]];
	if (!stricmp(pszFeatureName, pszLAname))
		{
		m_LAFeatureNameID = FeatureNameIdx + 1;
		return(0); // non-unique!
		}
	}

// not a duplicate
FeatureNameLen = (int)strlen(pszFeatureName);
if ((m_NxtszFeatureNameIdx + FeatureNameLen + 1) > (int)sizeof(m_szFeatureNames))
	return(0);		// no more space to hold feature name, treating as if a non-unique!
if (m_NumFeatureNames == cMaxRNAEFeatures)
	return(0);		// unable to hold any more FeatureNames, treating as if a non-unique!

m_szFeatureNamesIdx[m_NumFeatureNames++] = m_NxtszFeatureNameIdx;
m_NxtszFeatureNameIdx += FeatureNameLen + 1;
m_LAFeatureNameID = m_NumFeatureNames;
return(m_LAFeatureNameID);
}

uint32_t		// returned FeatureName identifier, 0 if unable to locate this FeatureName name
CRNAExpr::LocateFeatureNameID(char* pszFeatureName) // return unique identifier associated with this FeatureName
{
int32_t FeatureNameIdx;
char* pszLAFeatureName;

// with any luck the FeatureName will be same as the last accessed
if ((pszLAFeatureName = LocateFeatureName(m_LAFeatureNameID)) != nullptr)
	if (!stricmp(pszFeatureName, pszLAFeatureName))
		return(m_LAFeatureNameID);

	// iterate over all known FeatureNames
for (FeatureNameIdx = 0; FeatureNameIdx < m_NumFeatureNames; FeatureNameIdx++)
	{
	pszLAFeatureName = &m_szFeatureNames[m_szFeatureNamesIdx[FeatureNameIdx]];
	if (!stricmp(pszFeatureName, pszLAFeatureName))
		{
		m_LAFeatureNameID = FeatureNameIdx + 1;
		return(m_LAFeatureNameID);
		}
	}
return(0);
}

char*
CRNAExpr::LocateFeatureName(uint32_t FeatureNameID)
{
int32_t Idx;
Idx = FeatureNameID & 0x0fffffff;			// mask out any potential flag bits
if (Idx < 1 || Idx > m_NumFeatureNames)
	return(nullptr);
return(&(m_szFeatureNames[m_szFeatureNamesIdx[Idx - 1]])); 
}



uint32_t		// returned chrom identifier, 0 if unable to accept this chromosome name
CRNAExpr::AddChromName(char* pszChromName) // associate unique identifier with this chromosome name
{
int32_t ChromNameIdx;
int ChromNameLen;
char* pszLAname;

// with any luck the sequence name will be same as the last accessed
if ((pszLAname = LocateChromName(m_LAChromNameID)) != nullptr)
	if (!stricmp(pszChromName, pszLAname))
		return(m_LAChromNameID);

// iterate over all known chroms in case this chrom to add is a duplicate
for (ChromNameIdx = 0; ChromNameIdx < m_NumChromNames; ChromNameIdx++)
	if (!stricmp(pszChromName, &m_szChromNames[m_szChromNamesIdx[ChromNameIdx]]))
	{
	m_LAChromNameID = ChromNameIdx + 1;
	return(m_LAChromNameID);
	}

// chrom is not a duplicate
ChromNameLen = (int)strlen(pszChromName);
if ((m_NxtszChromNameIdx + ChromNameLen + 1) > (int)sizeof(m_szChromNames))
	return(0);
if (m_NumChromNames == cMaxChromNames)
	return(0);

m_szChromNamesIdx[m_NumChromNames++] = m_NxtszChromNameIdx;
strcpy(&m_szChromNames[m_NxtszChromNameIdx], pszChromName);
m_NxtszChromNameIdx += ChromNameLen + 1;
m_LAChromNameID = m_NumChromNames;
return(m_LAChromNameID);
}


uint32_t		// returned chrom identifier, 0 if unable to locate this chromosome name
CRNAExpr::LocateChromID(char* pszChrom) // return unique identifier associated with this chromosome name
{
int32_t ChromNameIdx;
char* pszLAname;

// with any luck the sequence name will be same as the last accessed
if ((pszLAname = LocateChromName(m_LAChromNameID)) != nullptr)
	if (!stricmp(pszChrom, pszLAname))
		return(m_LAChromNameID);

// iterate over all known chroms
for (ChromNameIdx = 0; ChromNameIdx < m_NumChromNames; ChromNameIdx++)
	if (!stricmp(pszChrom, &m_szChromNames[m_szChromNamesIdx[ChromNameIdx]]))
	{
		m_LAChromNameID = ChromNameIdx + 1;
		return(m_LAChromNameID);
	}
return(0);
}

char*
CRNAExpr::LocateChromName(uint32_t ChromID)
{
int32_t Idx = ChromID & 0x0fffffff;			// mask out any potential flag bits
if (Idx < 1 || (int32_t)Idx > m_NumChromNames)
	return(nullptr);
return(&m_szChromNames[m_szChromNamesIdx[Idx - 1]]);
}

int
CRNAExpr::AllocateRNAcnts(int EstNumRNACnts)	// initially allocate for this estimated total number of RNA feature counts, memory will be reallocd to hold more if subsequently required
{
size_t memreq;
if(m_pRNACntsMem != nullptr)
	{
#ifdef _WIN32
	free(m_pRNACntsMem);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (m_pRNACntsMem != MAP_FAILED)
		munmap(m_pRNACntsMem, m_AllocdRNACntsMem);
#endif
	m_pRNACntsMem = nullptr;
	}
m_AllocdRNACntsMem = 0;
m_UsedRNACntsMem = 0;

	// initial allocations, will be realloc'd to larger size if later required
memreq = (size_t)EstNumRNACnts * sizeof(int32_t);
#ifdef _WIN32
m_pRNACntsMem = (uint8_t*)malloc(memreq);	// initial and perhaps the only allocation
if (m_pRNACntsMem == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "AllocateSamples: Memory allocation of %zd bytes for sample features failed - %s", (int64_t)memreq, strerror(errno));
	Reset();
	return(eBSFerrMem);
	}
#else
m_pRNACntsMem = (uint8_t *)mmap(NULL, memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
if (m_pRNACntsMem == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "AllocateSamples: Memory allocation of %zd bytes through mmap() for sample features failed - %s", (int64_t)memreq, strerror(errno));
	m_pRNACntsMem = NULL;
	Reset();
	return(eBSFerrMem);
	}
#endif
m_AllocdRNACntsMem = memreq;
m_UsedRNACntsMem = 0;
return(eBSFSuccess);
}


int
CRNAExpr::ReallocRNAcnts(int32_t NumFeatCntsPerRow)	// realloc as may be required to have room for at least one more row of RNA feature counts
{
uint8_t *pReallocd;

if(m_pRNACntsMem == nullptr)
	return(AllocateRNAcnts(NumFeatCntsPerRow * 2000));

		// needing to allocate more memory?
if ((m_UsedRNACntsMem + (sizeof(int32_t) * NumFeatCntsPerRow)) >= m_AllocdRNACntsMem)
	{
	size_t memreq = m_AllocdRNACntsMem + ((size_t)NumFeatCntsPerRow * sizeof(int32_t) * 100);	// alloc extra to reduce number of realloc's subsequently required
#ifdef _WIN32
	pReallocd = (uint8_t *)realloc(m_pRNACntsMem, memreq);
	if (pReallocd == NULL)
		{
#else
	pReallocd = (uint8_t *)mremap(m_pRNACntsMem, m_AllocdRNACntsMem, memreq, MREMAP_MAYMOVE);
	if (pReallocd == MAP_FAILED)
		{
#endif
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "ReallocSamples: Memory reallocation to %zd bytes failed - %s", (int64_t)memreq, strerror(errno));
		return(eBSFerrMem);
		}
	m_pRNACntsMem = pReallocd;
	m_AllocdRNACntsMem = memreq;
	}
return(eBSFSuccess);
}

int
CRNAExpr::AllocateRNAhomozScores(int EstRNAhomozScores)	// initially allocate for this estimated total number of RNA vs. WGS homozygosity scores, memory will be reallocd to hold more if subsequently required
{
size_t memreq;
if(m_pRNAhomozScoresMem != nullptr)
	{
#ifdef _WIN32
	free(m_pRNAhomozScoresMem);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (m_pRNAhomozScoresMem != MAP_FAILED)
		munmap(m_pRNAhomozScoresMem, m_AllocdRNAhomozScoresMem);
#endif
	m_pRNAhomozScoresMem = nullptr;
	}
m_UsedRNAhomozScoresMem = 0;
m_AllocdRNAhomozScoresMem = 0;

	// initial allocations, will be realloc'd to larger sizes if later required
memreq = (size_t)EstRNAhomozScores * sizeof(double);
#ifdef _WIN32
m_pRNAhomozScoresMem = (uint8_t*)malloc(memreq);	// initial and perhaps the only allocation
if (m_pRNAhomozScoresMem == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "AllocateRNAhomozScores: Memory allocation of %zd bytes for scores failed - %s", (int64_t)memreq, strerror(errno));
	Reset();
	return(eBSFerrMem);
	}
#else
m_pRNAhomozScoresMem = (uint8_t *)mmap(NULL, memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
if (m_pRNAhomozScoresMem == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "AllocateRNAhomozScores: Memory allocation of %zd bytes through mmap() for sample features failed - %s", (int64_t)memreq, strerror(errno));
	m_pRNAhomozScoresMem = nullptr;
	Reset();
	return(eBSFerrMem);
	}
#endif
m_AllocdRNAhomozScoresMem = memreq;
m_UsedRNAhomozScoresMem = 0;
return(eBSFSuccess);
}


int
CRNAExpr::ReallocRNAhomozScores(int32_t NumScoresPerRow)	// realloc as may be required to have room (columns) for at least one more row of RNA homozygosity scores
{
double *pReallocd;
if(m_pRNAhomozScoresMem == nullptr)
	return(AllocateRNAhomozScores(NumScoresPerRow * 2000));
if ((m_UsedRNAhomozScoresMem + (sizeof(double) * NumScoresPerRow)) >= m_AllocdRNAhomozScoresMem)
	{
	size_t memreq = m_AllocdRNAhomozScoresMem + ((size_t)NumScoresPerRow * sizeof(double) * 100);	// alloc extra to reduce number of realloc's subsequently required
#ifdef _WIN32
	pReallocd = (double *)realloc(m_pRNAhomozScoresMem, memreq);
	if (pReallocd == NULL)
		{
#else
	pReallocd = (double *)mremap(m_pRNAhomozScoresMem, m_AllocdRNAhomozScoresMem, memreq, MREMAP_MAYMOVE);
	if (pReallocd == MAP_FAILED)
		{
#endif
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "ReallocRNAhomozScores: Memory reallocation to %zd bytes failed - %s", (int64_t)memreq, strerror(errno));
		return(eBSFerrMem);
		}
	m_pRNAhomozScoresMem = (uint8_t *)pReallocd;
	m_AllocdRNAhomozScoresMem = memreq;
	}
return(eBSFSuccess);
}


int
CRNAExpr::AllocateWGShomozScores(int EstWGShomozScores)	// initially allocate for this estimated total number of WGS vs. WGS homozygosity scores, memory will be reallocd to hold more if subsequently required
{
size_t memreq;
if(m_pWGShomozScoresMem != nullptr)
	{
#ifdef _WIN32
	free(m_pWGShomozScoresMem);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (m_pWGShomozScoresMem != MAP_FAILED)
		munmap(m_pWGShomozScoresMem, m_AllocdWGShomozScoresMem);
#endif
	m_pWGShomozScoresMem = nullptr;
	}
m_UsedWGShomozScoresMem = 0;
m_AllocdWGShomozScoresMem = 0;

	// initial allocations, will be realloc'd to larger sizes if later required
memreq = (size_t)EstWGShomozScores * sizeof(double);
#ifdef _WIN32
m_pWGShomozScoresMem = (uint8_t*)malloc(memreq);	// initial and perhaps the only allocation
if (m_pWGShomozScoresMem == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "AllocateGWShomozScores: Memory allocation of %zd bytes for scores failed - %s", (int64_t)memreq, strerror(errno));
	Reset();
	return(eBSFerrMem);
	}
#else
m_pWGShomozScoresMem = (uint8_t *)mmap(NULL, memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
if (m_pWGShomozScoresMem == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "AllocateGWShomozScores: Memory allocation of %zd bytes through mmap() for sample features failed - %s", (int64_t)memreq, strerror(errno));
	m_pWGShomozScoresMem = nullptr;
	Reset();
	return(eBSFerrMem);
	}
#endif
m_AllocdWGShomozScoresMem = memreq;
m_UsedWGShomozScoresMem = 0;
return(eBSFSuccess);
}

int
CRNAExpr::ReallocWGShomozScores(int32_t NumScoresPerRow)	// realloc as may be required to have room (columns) for at least one more row of WGS homozygosity scores
{
double *pReallocd;
if(m_pWGShomozScoresMem == nullptr)
	return(AllocateWGShomozScores(NumScoresPerRow * 1000));
if ((m_UsedWGShomozScoresMem + (sizeof(double) * NumScoresPerRow)) >= m_AllocdWGShomozScoresMem)
	{
	size_t memreq = m_AllocdWGShomozScoresMem + ((size_t)NumScoresPerRow * sizeof(double) * 100);	// alloc extra to reduce number of realloc's subsequently required
#ifdef _WIN32
	pReallocd = (double *)realloc(m_pWGShomozScoresMem, memreq);
	if (pReallocd == NULL)
		{
#else
	pReallocd = (double *)mremap(m_pWGShomozScoresMem, m_AllocdWGShomozScoresMem, memreq, MREMAP_MAYMOVE);
	if (pReallocd == MAP_FAILED)
		{
#endif
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "ReallocWGShomozScores: Memory reallocation to %zd bytes failed - %s", (int64_t)memreq, strerror(errno));
		return(eBSFerrMem);
		}
	m_pWGShomozScoresMem = (uint8_t *)pReallocd;
	m_AllocdWGShomozScoresMem = memreq;
	}
return(eBSFSuccess);
}


