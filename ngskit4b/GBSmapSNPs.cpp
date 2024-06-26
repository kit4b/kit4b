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
#include "GBSmapSNPs.h"

// forward declaration
int
Process(eModeGBSMapSNPs PMode,							// processing mode
						uint32_t ExprID,			// assign this experiment identifier for this SNP to haplotype analysis
						char *pszInNMFile,				// processing chromosome name mapping from this file
						char *pszInGBSFile,				// processing GBS SNP calls from this file
						char *pszOutFile);				// GBS haplotype calls output file

#ifdef _WIN32
int gbsmapsnps(int argc, char *argv[])
{
	// determine my process name
	_splitpath (argv[0], nullptr, nullptr, gszProcName, nullptr);
#else
int
gbsmapsnps(int argc, char **argv)
{
	// determine my process name
	CUtility::splitpath ((char *)argv[0], nullptr, gszProcName);
#endif
	int iFileLogLevel;			// level of file diagnostics
	int iScreenLogLevel;		// level of file diagnostics
	char szLogFile[_MAX_PATH];	// write diagnostics to this file
	int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure
	int NumberOfProcessors;		// number of installed CPUs

	 eModeGBSMapSNPs PMode;			// processing mode
	 uint32_t ExprID;				// assign this experiment identifier for this PBA analysis
	char szInNMFile[_MAX_PATH];		// processing chromosome name mapping from this file
	char szInGBSFile[_MAX_PATH];	// processing GBS SNP calls from this file
	char szOutFile[_MAX_PATH];		// GBS haplotype calls output file

	struct arg_lit *help = arg_lit0 ("h", "help", "print this help and exit");
	struct arg_lit *version = arg_lit0 ("v", "version,ver", "print version information and exit");
	struct arg_int *FileLogLevel = arg_int0 ("f", "FileLogLevel", "<int>", "Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
	struct arg_file *LogFile = arg_file0 ("F", "log", "<file>", "diagnostics log file");

	struct arg_int *pmode = arg_int0 ("m", "mode", "<int>", "processing mode: 0 Map SNP GBS to PBA GBS haplotypes generating a matrix, 1: combine two haplotype matrices");
	struct arg_int* exprid = arg_int0("e","exprid","<int>","assign this experiment identifier for haplotypes called (default 1)");
	struct arg_file *incnmapfile = arg_file1("I", "cnmap", "<file>", "-m0 mode: chromosome name mappings from this file, -m1 mode: processing haplotype matrix M2 (CSV format)");
	struct arg_file *ingbsfile = arg_file1("i", "in", "<file>", "-m0 mode: processing GBS SNP calls from this file, -m1 mode: processing haplotype matrix M1  (CSV format)");
	struct arg_file *outfile = arg_file1("o", "out", "<file>", "GBS haplotype calls output file (CSV format)");
	struct arg_end *end = arg_end (200);

	void *argtable[] = { help,version,FileLogLevel,LogFile,
						pmode,exprid,incnmapfile,ingbsfile,outfile,end };

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

		PMode = pmode->count ? (eModeGBSMapSNPs)pmode->ival[0] : eMGBSMDefault;
		if(PMode < eMGBSMDefault || PMode >= eMGBSMPlaceholder)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Processing mode is not supported '-m'\n");
			exit(1);
			}
		if(PMode == eMGBSMCombine && !incnmapfile->count)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Second matrix file must be specified '-I<filespec>' when combining two matrices\n");
			exit(1);
			}

		ExprID = exprid->count ? exprid->ival[0] : 1;
		if(ExprID < 1)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Experiment identifiers '-e%d' must be in range 1..10000\n",ExprID);
			exit(1);
			}
		else
			if(ExprID > 10000)
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Experiment identifiers '-e%d' must be in range 1..10000\n",ExprID);
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

	if(incnmapfile->count)
		{
		strcpy(szInNMFile, incnmapfile->filename[0]);
		CUtility::TrimQuotedWhitespcExtd(szInNMFile);
		if(szInNMFile[0] == '\0')
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: After removal of whitespace, no file specified with '-I<filespec>' option)\n");
			exit(1);
			}
		}
	else
		szInNMFile[0] = '\0';

	strcpy (szInGBSFile, ingbsfile->filename[0]);
	CUtility::TrimQuotedWhitespcExtd (szInGBSFile);
	if(szInGBSFile[0] == '\0')
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: After removal of whitespace, no input GBS SNP loci file specified with '-i<filespec>' option)\n");
		exit(1);
		}

	strcpy (szOutFile, outfile->filename[0]);
	CUtility::TrimQuotedWhitespcExtd (szOutFile);
	if(szOutFile[0] == '\0')
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: After removal of whitespace, no output GBS haplotype output file specified with '-o<filespec>' option)\n");
		exit(1);
		}


	gDiagnostics.DiagOut (eDLInfo, gszProcName, "Processing parameters:");
	const char *pszDescr;
	switch (PMode) {
		case eMGBSMDefault:
			pszDescr = "Map SNP GBS to PBA GBS haplotypes";
			break;
		case eMGBSMCombine:
			pszDescr = "Combine two haplotype matrices";
		}
	gDiagnostics.DiagOutMsgOnly (eDLInfo, "GBS SNP call processing : '%s'", pszDescr);
	gDiagnostics.DiagOutMsgOnly(eDLInfo, "Experiment identifier: %d", ExprID);
	if(PMode == eMGBSMCombine)
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "First matrix loaded from this input file : '%s'", szInNMFile);
	else
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "Processing GBS SNP calls from this input file : '%s'", szInGBSFile);
	if(PMode == eMGBSMCombine)
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "Second matrix loaded from this input file : '%s'", szInNMFile);
	else
		if(PMode == eMGBSMDefault && szInNMFile[0] != '\0')
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "Using chromosome name mappings from this input file : '%s'", szInNMFile);
	
	gDiagnostics.DiagOutMsgOnly (eDLInfo, "Writing GBS haplotypes to this output file : '%s'", szOutFile);

#ifdef _WIN32
	SetPriorityClass (GetCurrentProcess (), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start ();
	Rslt = 0;
	Rslt = Process(PMode,			// processing mode
						ExprID,			// experiment identifier
						szInNMFile,		// processing chromosome name mapping from this file
						szInGBSFile,	// processing GBS SNP calls from this file
						szOutFile);		// GBS haplotypes calls output file

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

int
Process(eModeGBSMapSNPs PMode,			// processing mode
						uint32_t ExprID,			// assign this experiment identifier for this SNP to haplotype analysis
						char *pszInNMFile,		// processing chromosome name mapping from this file
						char *pszInGBSFile,	// processing GBS SNP calls from this file
						char *pszOutFile)		// GBS haplotypes output file
{
int Rslt;
CGBSmapSNPs* pGBSmapSNPs;

if((pGBSmapSNPs = new CGBSmapSNPs) == nullptr)
	{
	gDiagnostics.DiagOut (eDLFatal, gszProcName, "Unable to instantiate instance of CGBSmapSNPs");
	return(eBSFerrObj);
	}
Rslt = pGBSmapSNPs->Process(PMode,ExprID,pszInNMFile,pszInGBSFile,pszOutFile);
delete pGBSmapSNPs;
return(Rslt);
}


CGBSmapSNPs::CGBSmapSNPs()
{
m_hOutFile = -1;
m_pInNMFile = nullptr;
m_pInGBSFile = nullptr;
m_pChromMappings = nullptr;
m_pProgenyFndrAligns = nullptr;
m_pszOutBuffer = nullptr;
Reset();
}

CGBSmapSNPs::~CGBSmapSNPs()
{
if(m_hOutFile != -1)
	close(m_hOutFile);
if(m_pInNMFile != nullptr)
	delete m_pInNMFile;
if(m_pInGBSFile != nullptr)
	delete m_pInGBSFile;
if(m_pChromMappings != nullptr)
	delete m_pChromMappings;
if(m_pszOutBuffer != nullptr)
	delete []m_pszOutBuffer;
if(m_pProgenyFndrAligns != nullptr)
	{
#ifdef _WIN32
	free(m_pProgenyFndrAligns);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pProgenyFndrAligns != MAP_FAILED)
		munmap(m_pProgenyFndrAligns, m_AllocdProgenyFndrAlignsMem);
#endif
	}
}

void 
CGBSmapSNPs::Reset(void)			// reset class states back to that immediately following class instantiation
{
if(m_hOutFile != -1)
	{
	close(m_hOutFile);
	m_hOutFile = -1;
	}
if(m_pInNMFile != nullptr)
	{
	delete m_pInNMFile;
	m_pInNMFile = nullptr;
	}
if(m_pInGBSFile != nullptr)
	{
	delete m_pInGBSFile;
	m_pInGBSFile = nullptr;
	}
if(m_pChromMappings != nullptr)
	{
	delete m_pChromMappings;
	m_pChromMappings = nullptr;
	}
m_NumChromMappings = 0;

if(m_pszOutBuffer != nullptr)
	{
	delete []m_pszOutBuffer;
	m_pszOutBuffer = nullptr;
	}

if(m_pProgenyFndrAligns != nullptr)
	{
#ifdef _WIN32
	free(m_pProgenyFndrAligns);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pProgenyFndrAligns != MAP_FAILED)
		munmap(m_pProgenyFndrAligns, m_AllocdProgenyFndrAlignsMem);
#endif
	m_pProgenyFndrAligns = nullptr;
	}

m_UsedProgenyFndrAligns = 0;
m_AllocdProgenyFndrAligns = 0;
m_AllocdProgenyFndrAlignsMem = 0;

m_OutBuffIdx = 0;
m_AllocOutBuff = 0;
m_PMode = eMGBSMDefault;
m_pszInNMFile = nullptr;
m_pszInGBSFile = nullptr;
m_pszOutFile = nullptr;

m_NumFounders = 0;
memset(m_FndrIDs,0,sizeof(m_FndrIDs));

m_NumProgenies = 0;
memset(m_ProgenyIDs,0,sizeof(m_ProgenyIDs));

m_LAReadsetNameID = 0;
m_NumReadsetNames = 0;
m_NxtszReadsetIdx = 0;
m_szReadsetNames[0] = '\0';
m_ExprID = 0;
m_LAChromNameID=0;
m_NumChromNames=0;
m_NxtszChromIdx=0;
memset(m_szChromNames,0,sizeof(m_szChromNames));
memset(m_szChromIdx,0,sizeof(m_szChromIdx));
}

uint32_t		// returned chrom identifier, 0 if unable to accept this chromosome name
CGBSmapSNPs::AddChrom(char* pszChrom) // associate unique identifier with this chromosome name
{
uint32_t ChromNameIdx;
int ChromNameLen;
char *pszLAname;

// with any luck the sequence name will be same as the last accessed
if((pszLAname = LocateChrom(m_LAChromNameID)) != nullptr)
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


uint32_t		// returned chrom identifier, 0 if unable to locate this chromosome name
CGBSmapSNPs::LocateChrom(char* pszChrom) // return unique identifier associated with this chromosome name
{
uint32_t ChromNameIdx;
char *pszLAname;

// with any luck the sequence name will be same as the last accessed
if((pszLAname = LocateChrom(m_LAChromNameID)) != nullptr)
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
CGBSmapSNPs::LocateChrom(uint32_t ChromID)
{
if(ChromID < 1 || ChromID > m_NumChromNames)
	return(nullptr);
return(&m_szChromNames[m_szChromIdx[ChromID-1]]);
}

int
CGBSmapSNPs::LoadNM(char* pszInNMFile)		// processing chromosome name mapping from this file)
{
int Rslt;
int32_t NumFields;
int32_t CurLineNumber;
uint32_t MapFromChromID;
uint32_t MapToChromID;
uint32_t ChromSize;
char *pszFromChrom;
char *pszToChrom;

tsChromMapping *pChromMapping;

if((m_pInNMFile = new CCSVFile) == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to instantiate CCSVfile");
	Reset();
	return(eBSFerrObj);
	}

if((Rslt = m_pInNMFile->Open(pszInNMFile)) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to open file: '%s'", pszInNMFile);
	Reset();
	return(Rslt);
	}

if(m_pChromMappings == nullptr)
	{
	if((m_pChromMappings = new tsChromMapping [cMaxChromMappings]) == nullptr)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to allocate memory for holding chromosome sizes and names");
		Reset();
		return(eBSFerrMem);
		}
	}

CurLineNumber = 0;
while ((Rslt = m_pInNMFile->NextLine()) > 0)				// onto next line containing fields
	{
	CurLineNumber++;
	if((NumFields = m_pInNMFile->GetCurFields()) != 3)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Chromosome name mapping file '%s' expected to contain 3 fields, it contains %d at line %d",pszInNMFile,NumFields,CurLineNumber);
		Reset();
		return(eBSFerrParse);
		}
	// parse out mapping from chromosome name, followed by mapping to chromosome name, followed by the mapping to chromosome size (bp)
	// 1st row is expected to be a header row
	if(CurLineNumber == 1 && m_pInNMFile->IsLikelyHeaderLine())
		continue;

	if(m_NumChromMappings == cMaxChromMappings)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Chromosome name mapping file '%s' contains more than %d chromosome mappings at line %d",pszInNMFile,cMaxChromMappings,CurLineNumber);
		Reset();
		return(eBSFerrParse);
		}

	m_pInNMFile->GetText(1, &pszFromChrom);			
	m_pInNMFile->GetText(2, &pszToChrom);			
	m_pInNMFile->GetInt(3, (int*)&ChromSize);
	MapFromChromID = AddChrom(pszFromChrom);
	MapToChromID = AddChrom(pszToChrom);
	pChromMapping = &m_pChromMappings[m_NumChromMappings++];
	pChromMapping->AliasChromID = MapFromChromID;
	pChromMapping->RefChromID = MapToChromID;
	pChromMapping->Size = ChromSize;
	}
delete m_pInNMFile;
m_pInNMFile = nullptr;

pChromMapping = m_pChromMappings;

return(eBSFSuccess);
}


// NOTE: Readsets are checked for uniqueness as readsets must be unique within a given readetset type
int32_t		// returned readset identifier, < 1 if unable to accept this readset name
CGBSmapSNPs::AddReadset(char* pszReadset, // associate unique identifier with this readset name
							uint8_t ReadsetType)	// 0: founder, 1: progeny, 2: control
{
int32_t ReadsetNameIdx;
int ReadsetNameLen;
char Type;
char *pszLAname;
Type = '0' + (char)ReadsetType;

// with any luck the sequence name will be same as the last accessed
if((pszLAname = LocateReadset(m_LAReadsetNameID)) != nullptr)
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

char* 
CGBSmapSNPs::LocateReadset(int32_t ReadsetID)
{
int Idx;
Idx = ReadsetID & 0x0fffffff;			// mask out any potential ReadsetType
if(Idx < 1 || Idx > m_NumReadsetNames)
	return(nullptr);
return(&(m_szReadsetNames[m_szReadsetIdx[Idx -1]+1])); // skipping lead char which is the ReadsetType
}


int32_t		// returned Readset identifier, < 1 if unable to locate this Readset name
CGBSmapSNPs::LocateReadset(char* pszReadset, // return unique identifier associated with this Readset name
							   uint8_t ReadsetType)	// 0: founder, 1: progeny, 2: control
{
int32_t ReadsetNameIdx;
char Type;
char *pszLAReadset;

Type = '0'+ (char)ReadsetType;
// with any luck the Readset name will be same as the last accessed
if((pszLAReadset = LocateReadset(m_LAReadsetNameID)) != nullptr)
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


// linear search through alias names checking for longest match on the prefix of pszChrom
// can't do a simple exact match as sometimes there may be a suffix appended - could be SNP loci etc. - by the GBS SNP loci file generator. UGH, UGH ....
tsChromMapping *
CGBSmapSNPs::LocateAliasChromMapping(char* pszChrom)	// locates and returns chromosome mapping for a alias name which matches the prefix of this chromosome 
{
uint32_t MapIdx;
char *pszAlias;
tsChromMapping *pChromMapping;

pChromMapping = m_pChromMappings;
for(MapIdx = 0; MapIdx < m_NumChromMappings; MapIdx++,pChromMapping++)
	{
	pszAlias = LocateChrom(pChromMapping->AliasChromID);
	if(!strnicmp(pszAlias,pszChrom,strlen(pszAlias)))
		return(pChromMapping);
	}
return(nullptr);
}

// inplace modification of sample identifier ensuring that sample identifiers have a prefix of 'S0' where original identifier starts with 'S[1-9]'
// additionally, if sample prefixed by 'Progeny:' then remove that prefix so that only the sample identifier is retained
int32_t		// updated strlen()
CGBSmapSNPs::MakeConsistentSampleName(char* pszSampleID)
{
char szUpdated[200];
char CurChr;
char *pszUpdated = szUpdated;
char *pszOrig = pszSampleID;
if(!strncmp("Progeny:", pszSampleID,8))
    pszOrig += 8;
while((CurChr = *pszUpdated++ = *pszOrig++) != '\0')
	{
	if(CurChr == 'S')
		{
		if(*pszOrig >= '1' && *pszOrig <= '9')
			*pszUpdated++ = '0';
		}
	}
strcpy(pszSampleID,szUpdated);
return((int32_t)strlen(szUpdated));
}

int
CGBSmapSNPs::LoadHaplotypeMatrix(char* pszInGBSFile)							// processing chromosome name mapping from this file)
{
int Rslt;
int F4Idx;
int ExptdNumFields;
int CurNumFields;
int CurLineNumber;
int ExprID;
char *pszChrom;
uint32_t RefChromID;
int Haplotypes;
int F4Field;
uint32_t SNPLoci;

tsProgenyFndrAligns ProgenyFndrAligns;
char *pszGBSReadset;

if((m_pInGBSFile = new CCSVFile) == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to instantiate CCSVfile");
	Reset();
	return(eBSFerrObj);
	}

if((Rslt = m_pInGBSFile->Open(pszInGBSFile)) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to open file: '%s'", pszInGBSFile);
	Reset();
	return(Rslt);
	}
m_pInGBSFile->SetMaxFields(cMaxProgenyReadsets + 5);		// could be quite a few F4s!

ExprID = 0;
SNPLoci = 0;

ExptdNumFields = 0;
CurNumFields = 0;
CurLineNumber = 0;
F4Field = 6;	// defaulting as the first F4

while ((Rslt = m_pInGBSFile->NextLine()) > 0)				// onto next line containing fields
	{
	CurLineNumber++;
	if(CurLineNumber == 1 && m_pInGBSFile->IsLikelyHeaderLine())
		{
		if((ExptdNumFields = m_pInGBSFile->GetCurFields()) < 4)		// expected to be at least 3 + 1 F4!
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "GBS SNP Genotyping matrix file '%s' expected to contain at least 4 fields at line %d",pszInGBSFile,CurLineNumber);
			Reset();
			return(eBSFerrParse);
			}

		// inside knowledge - GBS matrix file currently only has 2 founders!
		// dummy founders used
		if(!LocateReadset((char *)"Fa",0))
			{
			m_FndrIDs[m_NumFounders] = AddReadset((char *)"Fa",0);
			m_Fndrs2Proc[m_NumFounders++] = 0x01;
			}
		if(!LocateReadset((char *)"Fb",0))
			{
			m_FndrIDs[m_NumFounders] = AddReadset((char *)"Fb",0);
			m_Fndrs2Proc[m_NumFounders++] = 0x01;
			}
		// iterate over header fields and locate the F4 field, record it's index
		int32_t CurReadset;
		for(F4Field = 4,F4Idx=0; F4Field <= ExptdNumFields; F4Field++,F4Idx++)
			{
			m_pInGBSFile->GetText(F4Field, &pszGBSReadset);
			// could be of the form 'S99693' or 'S099693' or 'Progeny:S99693' or 'Progeny:S099693' dependent on the source generator
			// prefix sample identifiers with a leading 0, if missing, so as to maintain naming consistency
			MakeConsistentSampleName(pszGBSReadset);
			if((CurReadset = LocateReadset(pszGBSReadset,1)) == 0)
				CurReadset = AddReadset(pszGBSReadset,1);
			if(CurReadset != m_ProgenyIDs[F4Idx])
				m_ProgenyIDs[F4Idx] = CurReadset;
			}
		if(m_NumProgenies > 0 && m_NumProgenies != F4Idx)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "GBS Haplotype matrix file '%s' contains different progeny set",pszInGBSFile,pszGBSReadset);
			Reset();
			return(eBSFerrEntry);
			}
		m_NumProgenies = F4Idx;
		continue;
		}

	if((CurNumFields = m_pInGBSFile->GetCurFields()) != ExptdNumFields)
		{
		if(ExptdNumFields == 0)
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "GBS haplotype matrix file '%s' expected to contain header line at line %",pszInGBSFile,CurLineNumber);
		else
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "GBS haplotype matrix file '%s' expected to contain at %d fields (same number as header) at line %",pszInGBSFile,ExptdNumFields,CurLineNumber);
		Reset();
		return(eBSFerrParse);
		}
	memset(&ProgenyFndrAligns,0,sizeof(ProgenyFndrAligns));
	// expectation is that -
	// 	field 1 contains experiment identifier
	// 	field 2 contains chrom name
	//	field 3 contains loci
	//	field 4..n contains F4n haplotype call,-1: no coverage, 0: no haplotype call, 1: Fa only, 2: Fb only, 3: Fa and Fb
	m_pInGBSFile->GetInt(1, &ExprID);
	m_pInGBSFile->GetText(2, &pszChrom);
	RefChromID = AddChrom(pszChrom);
	m_pInGBSFile->GetInt(3, (int*)&SNPLoci);

	for(F4Field = 4,F4Idx=0; F4Field <= ExptdNumFields; F4Field++,F4Idx++)
		{
		memset(&ProgenyFndrAligns,0,sizeof(ProgenyFndrAligns));
		ProgenyFndrAligns.ReadsetID = m_ProgenyIDs[F4Idx];
		ProgenyFndrAligns.ExprID = ExprID;
		ProgenyFndrAligns.Loci = SNPLoci;
		ProgenyFndrAligns.Source=1;
		ProgenyFndrAligns.ChromID = RefChromID;
		bool bFa = false;
		bool bFb = false;
		m_pInGBSFile->GetInt(F4Field, (int*)&Haplotypes);
		if(Haplotypes < 0)		// -1 if progeny was unaligned at this loci
			continue;
		ProgenyFndrAligns.Alleles = 1;	// dummy allele to flag that progeny was aligned even though there may have been no haplotype called
		if(Haplotypes & 0x01)
			bFa = true;
		if(Haplotypes & 0x02)
			bFb = true;

		if(bFa)
			{
			BitsVectSet(0,ProgenyFndrAligns.ProgenyFounders);
			ProgenyFndrAligns.NumProgenyFounders++;
			}
		if(bFb)
			{
			BitsVectSet(1,ProgenyFndrAligns.ProgenyFounders);
			ProgenyFndrAligns.NumProgenyFounders++;
			}
		AddProgenyFndrAligns(&ProgenyFndrAligns);
		}
	}

delete m_pInGBSFile;
m_pInGBSFile = nullptr;
return(Rslt);
}


int
CGBSmapSNPs::LoadGBSSNPs(uint32_t MatrixID,			// identifies matrix source
				char* pszInGBSFile)							// processing chromosome name mapping from this file)
{
int Rslt;
int F4Idx;
int ExptdNumFields;
int CurNumFields;
int CurLineNumber;
char *pszAliasChrom;
char *pszFa;
char *pszFb;
char *pszFaSNPs;
char *pszFbSNPs;
int F4Field;
char *pszF4SNPs;
uint32_t SNPLoci;
tsProgenyFndrAligns ProgenyFndrAligns;
tsChromMapping *pChromMapping;
char *pszGBSReadset;

if((m_pInGBSFile = new CCSVFile) == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to instantiate CCSVfile");
	Reset();
	return(eBSFerrObj);
	}

if((Rslt = m_pInGBSFile->Open(pszInGBSFile)) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to open file: '%s'", pszInGBSFile);
	Reset();
	return(Rslt);
	}
m_pInGBSFile->SetMaxFields(cMaxProgenyReadsets + 5);		// could be quite a few F4s!

pszAliasChrom = nullptr;
SNPLoci = 0;
ExptdNumFields = 0;
CurNumFields = 0;
CurLineNumber = 0;
F4Field = 6;	// defaulting as the first F4

while ((Rslt = m_pInGBSFile->NextLine()) > 0)				// onto next line containing fields
	{
	CurLineNumber++;
	if(CurLineNumber == 1 && m_pInGBSFile->IsLikelyHeaderLine())
		{
		if((ExptdNumFields = m_pInGBSFile->GetCurFields()) < 6)		// expected to be at least 1 F4!
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "GBS SNP file '%s' expected to contain at least 6 fields at line %d",pszInGBSFile,CurLineNumber);
			Reset();
			return(eBSFerrParse);
			}

		m_NumFounders = 0;
		// get Fa name
		m_pInGBSFile->GetText(4, &pszFa);	// inside knowledge - 150QTL GBS snps formated file only has 2 founders!
		m_FndrIDs[m_NumFounders] = AddReadset(pszFa,0);
		m_Fndrs2Proc[m_NumFounders++] = 0x01;
		// get Fb name
		m_pInGBSFile->GetText(5, &pszFb);
		m_FndrIDs[m_NumFounders] = AddReadset(pszFb,0);
		m_Fndrs2Proc[m_NumFounders++] = 0x01;
		
		// iterate over header fields and locate the F4 field, record it's index
		for(F4Field = 6,F4Idx=0; F4Field <= ExptdNumFields; F4Field++,F4Idx++)
			{
			m_pInGBSFile->GetText(F4Field, &pszGBSReadset);
			m_ProgenyIDs[F4Idx] = AddReadset(pszGBSReadset,1);
			}
		m_NumProgenies = F4Idx;
		continue;
		}

	if((CurNumFields = m_pInGBSFile->GetCurFields()) != ExptdNumFields)
		{
		if(ExptdNumFields == 0)
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "GBS SNP Genotyping file '%s' expected to contain header line at line %",pszInGBSFile,CurLineNumber);
		else
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "GBS SNP Genotyping file '%s' expected to contain at %d fields (same number as header) at line %",pszInGBSFile,ExptdNumFields,CurLineNumber);
		Reset();
		return(eBSFerrParse);
		}
	memset(&ProgenyFndrAligns,0,sizeof(ProgenyFndrAligns));

	// expectation is that -
	//	field 1 contains chrom name appended with loci
	// 	field 2 contains chrom name
	//	field 3 contains loci
	//	field 4 contains Fa diploid base pair
	//	field 5 contains Fb diploid base pair
	//	field 6..n contains F4n diploid base pair

	// parse out each GBS SNP loci and accumulate bin counts
	// may have to strip off alias chrom suffixes
	m_pInGBSFile->GetText(2, &pszAliasChrom);			
	m_pInGBSFile->GetInt(3, (int*)&SNPLoci);

	if((pChromMapping = LocateAliasChromMapping(pszAliasChrom))==nullptr)	// currently a linear search
		continue;

	if(SNPLoci > pChromMapping->Size)		// clamp
		SNPLoci = pChromMapping->Size;

	uint8_t FaAlleles;
	uint8_t FbAlleles;
	
	// requiring that Fa and Fb are single allele - major - only so progeny parents can be differentiated
	m_pInGBSFile->GetText(4, &pszFaSNPs);
	if((FaAlleles = SNPs2Alleles(pszFaSNPs,true))==0)
		continue;
	m_pInGBSFile->GetText(5, &pszFbSNPs);
	if((FbAlleles = SNPs2Alleles(pszFbSNPs,true))==0)
		continue;
	if(FaAlleles == FbAlleles)	// if same alleles for Fa and Fb then can't differentiate
		continue;

	for(F4Field = 6,F4Idx=0; F4Field <= ExptdNumFields; F4Field++,F4Idx++)
		{
		memset(&ProgenyFndrAligns,0,sizeof(ProgenyFndrAligns));
		ProgenyFndrAligns.ExprID = m_ExprID;
		ProgenyFndrAligns.ReadsetID = m_ProgenyIDs[F4Idx];
		ProgenyFndrAligns.Loci = SNPLoci;
		ProgenyFndrAligns.Source=MatrixID;
		ProgenyFndrAligns.ChromID = pChromMapping->RefChromID;
		bool bFa = false;
		bool bFb = false;
		m_pInGBSFile->GetText(F4Field, &pszF4SNPs);
		ProgenyFndrAligns.Alleles = SNPs2Alleles(pszF4SNPs);
		if(ProgenyFndrAligns.Alleles == 0)	// if non-canonical alleles then assuming there was no alignment, these will be reported as having a '-1' haplotype
			continue;

		if(ProgenyFndrAligns.Alleles == FaAlleles)	// F4 is exactly matching Fa as a dirac?
			bFa = true;
		else
			if(ProgenyFndrAligns.Alleles == FbAlleles) // F4 is exactly matching Fb as a dirac?
				bFb = true;
			else
				if(ProgenyFndrAligns.Alleles == ((FaAlleles | FbAlleles) & 0x0aa)) // if 2 alleles present then these would have been returned as non-diracs
					bFa = bFb = true;
				else
					bFa = bFb = false;

		if(bFa)
			{
			BitsVectSet(0,ProgenyFndrAligns.ProgenyFounders);
			ProgenyFndrAligns.NumProgenyFounders++;
			}
		if(bFb)
			{
			BitsVectSet(1,ProgenyFndrAligns.ProgenyFounders);
			ProgenyFndrAligns.NumProgenyFounders++;
			}
		AddProgenyFndrAligns(&ProgenyFndrAligns);
		}
	}

delete m_pInGBSFile;
m_pInGBSFile = nullptr;
return(Rslt);
}

uint8_t			// returned PBA alleles
CGBSmapSNPs::SNPs2Alleles(char* pszSNPs,	// translate char representation of major/minor SNPs into it's packed byte allelic representation
				bool bMajorOnly)	// true if only major/major single allele to be returned as PBA

{
uint8_t Alleles;
char szlcSNPs[80];
char *pDst;
if(pszSNPs == nullptr || pszSNPs[0] == '\0')
	return(0);
pDst = szlcSNPs;
while((*pDst++ = tolower(*pszSNPs++))!= '\0');
if(szlcSNPs[0] == 'n' && szlcSNPs[1] == 'a')
	return(0);
if(bMajorOnly && szlcSNPs[1] != '\0')
	{
	if(szlcSNPs[0] != szlcSNPs[1])
		return(0);
	}
Alleles = 0;
switch(szlcSNPs[0]) {
		case 'a':
			if(szlcSNPs[1] == 'a' || szlcSNPs[1] == '\0')
				Alleles |= 0x03;
			else
				Alleles |= 0x02;
			break;
		case 'c':
			if(szlcSNPs[1] == 'c' || szlcSNPs[1] == '\0')
				Alleles |= 0x0c;
			else
				Alleles |= 0x08;
			break;
		case 'g':
			if(szlcSNPs[1] == 'g' || szlcSNPs[1] == '\0')
				Alleles |= 0x030;
			else
				Alleles |= 0x020;
			break;
		case 't':
			if(szlcSNPs[1] == 't' || szlcSNPs[1] == '\0')
				Alleles |= 0x0c0;
			else
				Alleles |= 0x080;
			break;
		default:
			return(0);
		}

	switch(szlcSNPs[1]) {
		case 'a':
			if(szlcSNPs[0] == 'a')
				Alleles |= 0x03;
			else
				Alleles |= 0x02;
			break;
		case 'c':
			if(szlcSNPs[0] == 'c')
				Alleles |= 0x0c;
			else
				Alleles |= 0x08;
			break;
		case 'g':
			if(szlcSNPs[0] == 'g')
				Alleles |= 0x030;
			else
				Alleles |= 0x020;
			break;
		case 't':
			if(szlcSNPs[0] == 't')
				Alleles |= 0x0c0;
			else
				Alleles |= 0x080;
			break;
		case '\0':
			break;
		default:
			return(0);
		}
return(Alleles);
}



int			// success or otherwise ( >= 0 success, < 0 if processing failed
CGBSmapSNPs::Process(eModeGBSMapSNPs PMode,			// processing mode
					uint32_t ExprID,			// assign this experiment identifier for this SNP to haplotype analysis
					 char* pszInNMFile,		// processing chromosome name mapping from this file
					 char* pszInGBSFile,	// processing GBS SNP calls from this file
					 char* pszOutFile)		// GBS haplotype matrix output file
{
int Rslt;
Reset();
m_PMode = PMode;
m_pszInNMFile = pszInNMFile;
m_pszInGBSFile = pszInGBSFile;
m_pszOutFile = pszOutFile;
m_ExprID = ExprID;

if(PMode == eMGBSMDefault)
	{
	// parse and load chromosome sizes and name mappings
	if(pszInNMFile != nullptr && pszInNMFile[0] != '\0' && (Rslt = LoadNM(pszInNMFile)) < eBSFSuccess)
		{
		Reset();
		return(Rslt);
		}

	// parse and load the GBS SNP loci file
	if((Rslt = LoadGBSSNPs(ExprID,pszInGBSFile)) < eBSFSuccess)
		{
		Reset();
		return(Rslt);
		}
	// sort progeny allele stack overlaps by readset.chrom.loci ascending
	m_mtqsort.SetMaxThreads(4);	// not expecting too many progeny loci so just use a few threads
	m_mtqsort.qsort(m_pProgenyFndrAligns, (int64_t)m_UsedProgenyFndrAligns, sizeof(tsProgenyFndrAligns), SortProgenyFndrAligns);

	// report haplotypes present for individual progenies
	for(uint32_t ProgIdx = 0; ProgIdx < (uint32_t)m_NumProgenies; ProgIdx++)
		ReportHaplotypesByProgeny(pszOutFile, m_ProgenyIDs[ProgIdx]);
	}
else // else must be combining two matrices
	{
	if((Rslt = LoadHaplotypeMatrix(pszInGBSFile)) < eBSFSuccess)
		{
		Reset();
		return(Rslt);
		}
	if((Rslt = LoadHaplotypeMatrix(pszInNMFile)) < eBSFSuccess)
		{
		Reset();
		return(Rslt);
		}
	}

// generate a matrix with chrom.loci (Y) by progeny (X)
m_mtqsort.qsort(m_pProgenyFndrAligns, (int64_t)m_UsedProgenyFndrAligns, sizeof(tsProgenyFndrAligns), SortProgenyChromLociReadset);

ReportMatrix(pszOutFile);

Reset();
return(Rslt);
}


int 
CGBSmapSNPs::ReportHaplotypesByProgeny(char* pszRsltsFileBaseName,		// haplotype results are written to this file base name with '.progeny,m_ExprID.csv' appended
											uint32_t ReadsetID)			// report on this progeny readset only, or if 0 then report on all progeny readsets
{
char szOutFile[_MAX_PATH];
tsProgenyFndrAligns *pCurPFA;
size_t PFAIdx;
int32_t FndrIdx;

int32_t PrevReadsetID;
int32_t PrevChromID;

char *pszProgenyReadset;
char *pszChrom;

if(m_pszOutBuffer == nullptr)
	{
	if((m_pszOutBuffer = new uint8_t[cOutBuffSize]) == nullptr)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "ReportHaplotypesByProgeny: Memory allocation of %u bytes for output buffering failed - %s", cOutBuffSize, strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
	m_AllocOutBuff = cOutBuffSize;
	}
m_OutBuffIdx = 0;
if(ReadsetID == 0)
	sprintf(szOutFile, "%s.progeny.%d.all.csv", pszRsltsFileBaseName, m_ExprID);
else
	sprintf(szOutFile, "%s.progeny.%d.%s.csv", pszRsltsFileBaseName, m_ExprID,LocateReadset(ReadsetID));
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
m_OutBuffIdx = sprintf((char*)m_pszOutBuffer, "\"ExprID\",\"Progeny\",\"Chrom\",\"Loci\"");
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
	if(!pCurPFA->NumProgenyFounders)
		continue;

	m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx],"%d,\"%s\",\"%s\",%u",m_ExprID, pszProgenyReadset, pszChrom, pCurPFA->Loci);
	for(FndrIdx=0; FndrIdx < m_NumFounders; FndrIdx++)
		{
		if(m_Fndrs2Proc[FndrIdx] & 0x01)
			m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], ",%d", BitsVectTest(FndrIdx, pCurPFA->ProgenyFounders) ? 1 : 0);
		}
	m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], "\n");
	if((m_OutBuffIdx + 1000) > m_AllocOutBuff)
		{
		if(!CUtility::RetryWrites(m_hOutFile, m_pszOutBuffer, m_OutBuffIdx))
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "ReportHaplotypesByProgeny: Fatal error in RetryWrites()");
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
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "ReportHaplotypesByProgeny: Fatal error in RetryWrites()");
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
CGBSmapSNPs::ReportMatrix(char* pszRsltsFileBaseName)		// matrix results are written to this file base name with 'matrix.csv' appended
{
char szOutFile[_MAX_PATH];
tsProgenyFndrAligns* pCurPFA;
size_t PFAIdx;
uint32_t FndrIdx;
uint32_t ProgIdx;
int32_t ProgenyHaplotypes[cMaxProgenyReadsets];
uint8_t ProgHaplotype;
int32_t CurChromID;
uint32_t CurLoci;
uint32_t CurExprID;
bool bNewLoci;
bool bProgFndrs;

char* pszChrom;

if(m_pszOutBuffer == nullptr)
{
	if((m_pszOutBuffer = new uint8_t[cOutBuffSize]) == nullptr)
	{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "ReportMatrix: Memory allocation of %u bytes for output buffering failed - %s", cOutBuffSize, strerror(errno));
		Reset();
		return(eBSFerrMem);
	}
	m_AllocOutBuff = cOutBuffSize;
}
m_OutBuffIdx = 0;
sprintf(szOutFile, "%s.%d.matrix.csv", pszRsltsFileBaseName,m_ExprID);
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
m_OutBuffIdx = sprintf((char*)m_pszOutBuffer, "\"ExprID\",\"Chrom\",\"Loci\"");

for(ProgIdx = 0; ProgIdx < m_NumProgenies; ProgIdx++)
	m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], ",\"Progeny:%s\"", LocateReadset(m_ProgenyIDs[ProgIdx]));

bNewLoci = false;
bProgFndrs = false;
memset(ProgenyHaplotypes,-1, sizeof(ProgenyHaplotypes));
CurLoci = 0;
CurChromID = 0;
CurExprID = 0;
pCurPFA = &m_pProgenyFndrAligns[0];
for(PFAIdx = 0; PFAIdx < m_UsedProgenyFndrAligns; PFAIdx++, pCurPFA++)
	{
	if(pCurPFA->ChromID != CurChromID || pCurPFA->Loci != CurLoci ||  pCurPFA->ExprID != CurExprID)
		{
		if(bNewLoci && bProgFndrs)
			{
			m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], "\n%d,\"%s\",%u", CurExprID,pszChrom, CurLoci);
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
			bProgFndrs = false;
			}
		CurExprID = pCurPFA->ExprID;
		CurLoci = pCurPFA->Loci;
		CurChromID = pCurPFA->ChromID;
		pszChrom = LocateChrom(CurChromID);
		memset(ProgenyHaplotypes,-1, sizeof(ProgenyHaplotypes));
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

	if(pCurPFA->Alleles != 0)
		{
		ProgHaplotype = 0;
		if(pCurPFA->NumProgenyFounders)
			{
			bProgFndrs=true;
			for(FndrIdx = m_NumFounders; FndrIdx >= 1; FndrIdx--)
				{
				ProgHaplotype <<= 1;
				if(m_Fndrs2Proc[FndrIdx-1] & 0x01)
					ProgHaplotype |= BitsVectTest(FndrIdx-1, pCurPFA->ProgenyFounders) ? 1 : 0;
				}
			}
		ProgenyHaplotypes[ProgIdx] = ProgHaplotype;
		}
	else
		ProgenyHaplotypes[ProgIdx] = -1;
	}
if(bNewLoci && bProgFndrs)
	{
	m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], "\n%d,\"%s\",%u", CurExprID,pszChrom, CurLoci);
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

uint32_t									// returned index+1 into  m_pProgenyFndrAligns[] to allocated and initialised ProgenyFndrAligns, 0 if errors
CGBSmapSNPs::AddProgenyFndrAligns(tsProgenyFndrAligns *pInitProgenyFndrAligns)	// allocated tsProgenyFndrAligns to be initialised with a copy of pInitProgenyFndrAligns
{
uint32_t ToAllocdProgenyFndrAligns;
tsProgenyFndrAligns *pProgenyFndrAligns;
size_t memreq;

if (m_pProgenyFndrAligns == nullptr)					// may be nullptr first time in
	{
	memreq = cAllocProgenyFndrAligns * sizeof(tsProgenyFndrAligns);
#ifdef _WIN32
	m_pProgenyFndrAligns = (tsProgenyFndrAligns*)malloc((size_t)memreq);
	if (m_pProgenyFndrAligns == nullptr)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddWinBinCnts: Memory allocation of %zd bytes failed", (int64_t)memreq);
		return(0);
		}
#else
	m_pProgenyFndrAligns = (tsProgenyFndrAligns*)mmap(nullptr, (size_t)memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	if (m_pProgenyFndrAligns == MAP_FAILED)
		{
		m_pProgenyFndrAligns = nullptr;
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddWinBinCnts: Memory allocation of %zd bytes through mmap()  failed", (int64_t)memreq, strerror(errno));
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
		if (pProgenyFndrAligns == nullptr)
			{
#else
		pProgenyFndrAligns = (tsProgenyFndrAligns*)mremap(m_pProgenyFndrAligns, m_AllocdProgenyFndrAlignsMem, memreq, MREMAP_MAYMOVE);
		if (pProgenyFndrAligns == MAP_FAILED)
			{
#endif
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddProgenyFndrAligns: Memory reallocation to %zd bytes failed - %s", memreq, strerror(errno));
			return(0);
			}
		m_pProgenyFndrAligns = pProgenyFndrAligns;
		m_AllocdProgenyFndrAlignsMem = memreq;
		m_AllocdProgenyFndrAligns =ToAllocdProgenyFndrAligns;
		}
pProgenyFndrAligns = &m_pProgenyFndrAligns[m_UsedProgenyFndrAligns++];
if(pInitProgenyFndrAligns != nullptr)
	*pProgenyFndrAligns = *pInitProgenyFndrAligns;
else
	memset(pProgenyFndrAligns,0,sizeof(tsProgenyFndrAligns));
return(m_UsedProgenyFndrAligns);
}

inline void
CGBSmapSNPs::BitsVectSet(uint16_t Bit,		// bit to set, range 0..cMaxBitVectBits-1
	tsBitsVect& BitsVect)
{
	BitsVect.Bits[Bit / 64] |= ((uint64_t)0x01 << (Bit % 64)); // 64bit words!
}

inline void
CGBSmapSNPs::BitsVectReset(uint16_t Bit,		// bit to reset, range 0..cMaxBitVectBits-1
	tsBitsVect& BitsVect)
{
	BitsVect.Bits[Bit / 64] &= ~((uint64_t)0x01 << (Bit % 64)); // 64bit words!
}


inline bool
CGBSmapSNPs::BitsVectTest(uint16_t Bit,		// bit to test, range 0..cMaxBitVectBits-1
	tsBitsVect& BitsVect)
{
	return(BitsVect.Bits[Bit / 64] & ((uint64_t)0x01 << (Bit % 64)) ? true : false);
}

inline bool
CGBSmapSNPs::BitsVectEqual(tsBitsVect& BitsVectA,	// compare for equality
	tsBitsVect& BitsVectB)
{
	if (!memcmp(&BitsVectA, &BitsVectB, sizeof(tsBitsVect)))
		return(true);
	return(false);
}

inline void
CGBSmapSNPs::BitsVectInitialise(bool Set,			// if true then initialise all bits as set, otherwise initialise all bits as reset
	tsBitsVect& BitsVect)
{
	uint8_t InitWith;
	if (Set)
		InitWith = 0xff;
	else
		InitWith = 0;
	memset(&BitsVect, InitWith, sizeof(tsBitsVect));
}

inline uint32_t
CGBSmapSNPs::BitsVectCount(tsBitsVect& BitsVect)		// count number of set bits
{
	uint32_t Count = 0;
	uint64_t* pWord = BitsVect.Bits;
	uint64_t Bit;
	for (uint32_t WordIdx = 0; WordIdx < cBitVectWords; WordIdx++, pWord++)
	{
		if (*pWord == 0)
			continue;
		for (Bit = 0x01; Bit != 0; Bit <<= 1)
			if (*pWord & Bit)
				Count++;
	}
	return(Count);
}

inline uint32_t
CGBSmapSNPs::BitsVectUnion(tsBitsVect& BitsVectA, tsBitsVect& BitsVectB)		// union (effective BitsVectA |= BitsVectB) bits in BitsVectA with BitsVectB with BitsVectA updated, returns number of bits set in BitsVectA 
{
	uint32_t Count = 0;
	uint64_t* pWordA = BitsVectA.Bits;
	uint64_t* pWordB = BitsVectB.Bits;
	uint64_t Bit;
	for (uint32_t WordIdx = 0; WordIdx < cBitVectWords; WordIdx++, pWordA++, pWordB++)
	{
		*pWordA |= *pWordB;
		if (*pWordA == 0)
			continue;
		for (Bit = 0x01; Bit != 0; Bit <<= 1)
			if (*pWordA & Bit)
				Count++;
	}
	return(Count);
}

inline uint32_t
CGBSmapSNPs::BitsVectIntersect(tsBitsVect& BitsVectA, tsBitsVect& BitsVectB)	// intersect (effective BitsVectA &= BitsVectB) of bits in BitsVectA with BitsVectB with BitsVectA updated, returns number of set bits in BitsVectA
{
	uint32_t Count = 0;
	uint64_t* pWordA = BitsVectA.Bits;
	uint64_t* pWordB = BitsVectB.Bits;
	uint64_t Bit;
	for (uint32_t WordIdx = 0; WordIdx < cBitVectWords; WordIdx++, pWordA++, pWordB++)
	{
		*pWordA &= *pWordB;
		if (*pWordA == 0)
			continue;
		for (Bit = 0x01; Bit != 0; Bit <<= 1)
			if (*pWordA & Bit)
				Count++;
	}
	return(Count);
}

uint32_t
CGBSmapSNPs::BitsVectClear(tsBitsVect& BitsVectA, tsBitsVect& BitsVectB)	    // clear bits in BitsVectA which are set in BitsVectB with BitsVectA updated, returns number of set bits in BitsVectA  
{
	uint32_t Count = 0;
	uint64_t* pWordA = BitsVectA.Bits;
	uint64_t* pWordB = BitsVectB.Bits;
	uint64_t Bit;
	for (uint32_t WordIdx = 0; WordIdx < cBitVectWords; WordIdx++, pWordA++, pWordB++)
	{
		*pWordA &= ~(*pWordB);
		if (*pWordA == 0)
			continue;
		for (Bit = 0x01; Bit != 0; Bit <<= 1)
			if (*pWordA & Bit)
				Count++;
	}
	return(Count);
}

// sorting by ReadsetID.ChromID.StartLoci.ExprID.Source ascending
int
CGBSmapSNPs::SortProgenyFndrAligns(const void* arg1, const void* arg2)
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
if(pEl1->ExprID < pEl2->ExprID)
	return(-1);
if(pEl1->ExprID > pEl2->ExprID)
	return(1);
if(pEl1->Source < pEl2->Source)
	return(-1);
if(pEl1->Source > pEl2->Source)
	return(1);

return(0);
}

// sorting by ChromID.StartLoci.ExprID.Readset.Source ascending
int
CGBSmapSNPs::SortProgenyChromLociReadset(const void* arg1, const void* arg2)
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
if(pEl1->ExprID < pEl2->ExprID)
	return(-1);
if(pEl1->ExprID > pEl2->ExprID)
	return(1);
if(pEl1->ReadsetID < pEl2->ReadsetID)
	return(-1);
if(pEl1->ReadsetID > pEl2->ReadsetID)
	return(1);
if(pEl1->Source < pEl2->Source)
	return(-1);
if(pEl1->Source > pEl2->Source)
	return(1);
return(0);
}


