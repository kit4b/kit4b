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
Process(eModeGBSMapSNPs PMode,			// processing mode
						int BinSize,					// SNP loci counts are accumulated into this sized (Kbp) non-overlapping bins
						char *pszGBSName,				// column name referencing the GBS readset which is to be processed
						char *pszInNMFile,		// processing chromosome name mapping from this file
						char *pszInGBSFile,	// processing GBS SNP calls from this file
						char *pszOutFile);		// windowed GBS SNP calls output file

#ifdef _WIN32
int gbsmapsnps(int argc, char *argv[])
{
	// determine my process name
	_splitpath (argv[0], NULL, NULL, gszProcName, NULL);
#else
int
gbsmapsnps(int argc, char **argv)
{
	// determine my process name
	CUtility::splitpath ((char *)argv[0], NULL, gszProcName);
#endif
	int iFileLogLevel;			// level of file diagnostics
	int iScreenLogLevel;		// level of file diagnostics
	char szLogFile[_MAX_PATH];	// write diagnostics to this file
	int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure
	int NumberOfProcessors;		// number of installed CPUs

	 eModeGBSMapSNPs PMode;			// processing mode
	int BinSize;					// SNP loci counts are accumulated into this sized (Kbp) non-overlapping bins
	char szGBSName[80];				// column name referencing the GBS readset which is to be processed
	char szInNMFile[_MAX_PATH];		// processing chromosome name mapping from this file
	char szInGBSFile[_MAX_PATH];	// processing GBS SNP calls from this file
	char szOutFile[_MAX_PATH];		// windowed GBS SNP calls output file

	struct arg_lit *help = arg_lit0 ("h", "help", "print this help and exit");
	struct arg_lit *version = arg_lit0 ("v", "version,ver", "print version information and exit");
	struct arg_int *FileLogLevel = arg_int0 ("f", "FileLogLevel", "<int>", "Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
	struct arg_file *LogFile = arg_file0 ("F", "log", "<file>", "diagnostics log file");

	struct arg_int *pmode = arg_int0 ("m", "mode", "<int>", "processing mode: 0 default");

	struct arg_int *accumbinsize = arg_int0 ("b", "binsize", "<int>", "SNP loci counts accumulated into non-overlapping bins of this size in Kbp (defaults to 100Kbp)");
	struct arg_str *gbsname = arg_str1("r","gbsname","<str>","header field name referencing the GBS readset which is to be processed");

	struct arg_file *incnmapfile = arg_file1("I", "cnmap", "<file>", "chromosome name mappings from this file (CSV format)");
	struct arg_file *ingbsfile = arg_file1("i", "in", "<file>", "processing GBS SNP calls from this file (CSV format)");
	struct arg_file *outfile = arg_file1("o", "out", "<file>", "windowed GBS SNP loci counts output file (CSV format)");
	struct arg_end *end = arg_end (200);

	void *argtable[] = { help,version,FileLogLevel,LogFile,
						pmode,accumbinsize,gbsname,incnmapfile,ingbsfile,outfile,end };

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

		PMode = pmode->count ? (eModeGBSMapSNPs)pmode->ival[0] : eMGBSMdefault;
		if(PMode < eMGBSMdefault || PMode > eMGBSMdefault)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Currently only the one processing mode is supported '-m0'\n");
			exit(1);
			}


		BinSize = accumbinsize->count ? accumbinsize->ival[0] : 100;
		if(BinSize < 10 || BinSize > 2000)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Accumulating non-overlapping bin size '-b%dKbp' specified outside of range 10..2000 Kbp\n", BinSize);
			exit(1);
			}

		szGBSName[0] = '\0';
		if(gbsname->count)
			{
			strncpy(szGBSName, gbsname->sval[0], 80);
			szGBSName[79] = '\0';
			CUtility::TrimQuotedWhitespcExtd(szGBSName);
			CUtility::CleanText(szGBSName);
			}
		if(szGBSName[0] == '\0')
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: GBS readset name is empty after whitespace trimming, must be specified with '-r<readset>'\n", BinSize);
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


	strcpy (szInNMFile, incnmapfile->filename[0]);
	CUtility::TrimQuotedWhitespcExtd (szInNMFile);
	if(szInNMFile[0] == '\0')
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: After removal of whitespace, no chromosome name mapping file specified with '-c<filespec>' option)\n");
		exit(1);
		}

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
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: After removal of whitespace, no output windowed GBS SNP loci output file specified with '-o<filespec>' option)\n");
		exit(1);
		}


	gDiagnostics.DiagOut (eDLInfo, gszProcName, "Processing parameters:");
	const char *pszDescr;
	switch (PMode) {
		case eMGBSMdefault:
			pszDescr = "Binning input GBS SNP loci into fixed size windowed bins as loci counts";
			break;
		}

	gDiagnostics.DiagOutMsgOnly (eDLInfo, "Accumulating counts in non-overlapping windowed bins of size : %dKbp", BinSize);
	gDiagnostics.DiagOutMsgOnly (eDLInfo, "Filtering for GBS readsets : '%s'", szGBSName);
	
	gDiagnostics.DiagOutMsgOnly (eDLInfo, "Using chromosome name mappings from this input file : '%s'", szInNMFile);
	gDiagnostics.DiagOutMsgOnly (eDLInfo, "Processing GBS SNP calls from this input file : '%s'", szInGBSFile);
	gDiagnostics.DiagOutMsgOnly (eDLInfo, "Writing windowed GBS SNP loci counts to this output file : '%s'", szOutFile);

#ifdef _WIN32
	SetPriorityClass (GetCurrentProcess (), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start ();
	Rslt = 0;
	Rslt = Process(PMode,			// processing mode
						BinSize,					// SNP loci counts are accumulated into this sized (Kbp) non-overlapping bins
						szGBSName,				// column name referencing the GBS readset which is to be processed
						szInNMFile,		// processing chromosome name mapping from this file
						szInGBSFile,	// processing GBS SNP calls from this file
						szOutFile);		// windowed GBS SNP calls output file

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
						int BinSize,					// SNP loci counts are accumulated into this sized (Kbp) non-overlapping bins
						char *pszGBSName,				// column name referencing the GBS readset which is to be processed
						char *pszInNMFile,		// processing chromosome name mapping from this file
						char *pszInGBSFile,	// processing GBS SNP calls from this file
						char *pszOutFile)		// windowed GBS SNP calls output file
{
int Rslt;
CGBSmapSNPs* pGBSmapSNPs;

if((pGBSmapSNPs = new CGBSmapSNPs) == NULL)
	{
	gDiagnostics.DiagOut (eDLFatal, gszProcName, "Unable to instantiate instance of CGBSmapSNPs");
	return(eBSFerrObj);
	}
Rslt = pGBSmapSNPs->Process(PMode,BinSize,pszGBSName,pszInNMFile,pszInGBSFile,pszOutFile);
delete pGBSmapSNPs;
return(Rslt);
}


CGBSmapSNPs::CGBSmapSNPs()
{
m_hOutFile = -1;
m_pInNMFile = NULL;
m_pInGBSFile = NULL;
m_pChromMappings = NULL;
m_pWinBins = NULL;
m_pszOutBuffer = NULL;
Reset();
}

CGBSmapSNPs::~CGBSmapSNPs()
{
if(m_hOutFile != -1)
	close(m_hOutFile);
if(m_pInNMFile != NULL)
	delete m_pInNMFile;
if(m_pInGBSFile != NULL)
	delete m_pInGBSFile;
if(m_pChromMappings != NULL)
	delete m_pChromMappings;
if(m_pWinBins != NULL)
	delete m_pWinBins;
if(m_pszOutBuffer != NULL)
	delete []m_pszOutBuffer;
}

void 
CGBSmapSNPs::Reset(void)			// reset class states back to that immediately following class instantiation
{
if(m_hOutFile != -1)
	{
	close(m_hOutFile);
	m_hOutFile = -1;
	}
if(m_pInNMFile != NULL)
	{
	delete m_pInNMFile;
	m_pInNMFile = NULL;
	}
if(m_pInGBSFile != NULL)
	{
	delete m_pInGBSFile;
	m_pInGBSFile = NULL;
	}
if(m_pChromMappings != NULL)
	{
	delete m_pChromMappings;
	m_pChromMappings = NULL;
	}
m_NumChromMappings = 0;

if(m_pWinBins != NULL)
	{
	delete m_pWinBins;
	m_pWinBins = NULL;
	}
if(m_pszOutBuffer != NULL)
	{
	delete []m_pszOutBuffer;
	m_pszOutBuffer = NULL;
	}


m_OutBuffIdx = 0;
m_AllocOutBuff = 0;

m_NumWinBins = 0;
m_BinSize = 0;
m_PMode = eMGBSMdefault;
m_BinSize = 0;
m_pszGBSName = NULL;
m_pszInNMFile = NULL;
m_pszInGBSFile = NULL;
m_pszOutFile = NULL;

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


uint32_t		// returned chrom identifier, 0 if unable to locate this chromosome name
CGBSmapSNPs::LocateChrom(char* pszChrom) // return unique identifier associated with this chromosome name
{
uint32_t ChromNameIdx;
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
CGBSmapSNPs::LocateChrom(uint32_t ChromID)
{
if(ChromID < 1 || ChromID > m_NumChromNames)
	return(NULL);
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
uint32_t StartLoci;
uint32_t ReqBins;
uint32_t CMIdx;
tsWinBin *pBin;
char *pszFromChrom;
char *pszToChrom;

tsChromMapping *pChromMapping;

if((m_pInNMFile = new CCSVFile) == NULL)
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

if(m_pChromMappings == NULL)
	{
	if((m_pChromMappings = new tsChromMapping [cMaxChromMappings]) == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to allocate memory for holding chromosome sizes and names");
		Reset();
		return(eBSFerrMem);
		}
	}

ReqBins = 0;
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
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Chromosome name mapping file '%s' contains in excess of %d chromosome mappings at line %d",pszInNMFile,cMaxChromMappings,CurLineNumber);
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
	ReqBins += (pChromMapping->Size+m_BinSize-1) / m_BinSize;
	}
delete m_pInNMFile;
m_pInNMFile = NULL;

m_pWinBins = new tsWinBin [ReqBins];
memset(m_pWinBins,0,sizeof(tsWinBin)*ReqBins);
pChromMapping = m_pChromMappings;
pBin = m_pWinBins;
m_NumWinBins = 0;
for(CMIdx = 0; CMIdx < m_NumChromMappings; CMIdx++,pChromMapping++)
	{
	pChromMapping->StartBinID = m_NumWinBins + 1;
	StartLoci = 0;
	while(StartLoci < pChromMapping->Size) {
		pBin->BinID = ++m_NumWinBins;
		pChromMapping->EndBinID = m_NumWinBins;
		pBin->RefChromID = pChromMapping->RefChromID;
		pBin->StartLoci = StartLoci;
		StartLoci += m_BinSize-1;
		if(StartLoci > pChromMapping->Size)
			StartLoci = pChromMapping->Size - 1;
		pBin->EndLoci = StartLoci;
		StartLoci += 1;
		pBin++;
		}
	}
return(eBSFSuccess);
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
return(NULL);
}

int
CGBSmapSNPs::LoadGBSSNPs(char* pszGBSName,				// column name referencing the GBS readset which is to be processed
				char* pszInGBSFile)							// processing chromosome name mapping from this file)
{
int Rslt;
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
tsChromMapping *pChromMapping;
tsWinBin *pWinBin;
char *pszGBSReadset;

if((m_pInGBSFile = new CCSVFile) == NULL)
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
m_pInGBSFile->SetMaxFields(cMaxGBSF4s + 5);		// could be quite a few F4s!

pszAliasChrom = NULL;
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
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "GBS SNP Genotyping file '%s' expected to contain at least 6 fields at line %",pszInGBSFile,CurLineNumber);
			Reset();
			return(eBSFerrParse);
			}

		// get Fa name
		m_pInGBSFile->GetText(4, &pszFa);
		strncpy(m_szFndrA,pszFa,sizeof(m_szFndrA)-1);
		m_szFndrA[sizeof(m_szFndrA)-1] = '\0';
		// get Fb name
		m_pInGBSFile->GetText(5, &pszFb);
		strncpy(m_szFndrB,pszFb,sizeof(m_szFndrB)-1);
		m_szFndrB[sizeof(m_szFndrB)-1] = '\0';

		// iterate over header fields and locate the F4 field, record it's index
		for(F4Field = 6; F4Field <= ExptdNumFields; F4Field++)
			{
			m_pInGBSFile->GetText(F4Field, &pszGBSReadset);
			if(!stricmp(pszGBSName,pszGBSReadset))
				break;
			}
		if(F4Field > ExptdNumFields)	// if can't locate field matching readset to report on...
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "GBS SNP Genotyping file '%s' has no header field matching requested GBS readset \"%s\" at line %",pszInGBSFile,pszInGBSFile,CurLineNumber);
			Reset();
			return(eBSFerrParse);
			}
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

	if((pChromMapping = LocateAliasChromMapping(pszAliasChrom))==NULL)	// currently a linear search
		continue;
	if(SNPLoci > pChromMapping->Size)		// clamp
		SNPLoci = pChromMapping->Size;

	pWinBin = &m_pWinBins[pChromMapping->StartBinID + (SNPLoci / m_BinSize)-1];
	pWinBin->LociCounts++;

	m_pInGBSFile->GetText(4, &pszFaSNPs);
	m_pInGBSFile->GetText(5, &pszFbSNPs);
	m_pInGBSFile->GetText(F4Field, &pszF4SNPs);
	bool bFa = false;
	bool bFb = false;
	bool bFaFb = false;
	bool bCntNA = true;
	if(stricmp(pszF4SNPs,"NA"))	// ensure there was a F4 SNP call
		{
		if(!stricmp(pszF4SNPs,pszFaSNPs))	// F4 is exactly matching Fa and/or Fb?
			bFa = true;
		if(!stricmp(pszF4SNPs,pszFbSNPs))
			bFb = true;
		if(!(bFa || bFb))		// if no exact match to either founder then check for match on allele
			{
			if(stricmp(pszFaSNPs,"NA"))
				{
				if(tolower(pszF4SNPs[0]) == tolower(pszFaSNPs[0]) ||	// F4 has an allele matching a Fa allele?
					tolower(pszF4SNPs[1]) == tolower(pszFaSNPs[0]) ||
					tolower(pszF4SNPs[0]) == tolower(pszFaSNPs[1]) ||
					tolower(pszF4SNPs[1]) == tolower(pszFaSNPs[1]))
					bFa = true;
				}
			if(stricmp(pszFbSNPs,"NA"))
				{
				if(tolower(pszF4SNPs[0]) == tolower(pszFbSNPs[0]) ||	// F4 has an allele matching a Fb allele?
					tolower(pszF4SNPs[1]) == tolower(pszFbSNPs[0]) ||
					tolower(pszF4SNPs[0]) == tolower(pszFbSNPs[1]) ||
					tolower(pszF4SNPs[1]) == tolower(pszFbSNPs[1]))
					bFb = true;
				}
			}
		if(bFa || bFb)
			bCntNA = false;
		if(bFa && bFb)
			{
			bFaFb = true;
			bFa = false;
			bFb = false;
			}
		}
	if(bCntNA)
		pWinBin->CntNA++;
	else
		{
		if(bFa)
			pWinBin->CntExclFa++;
		if(bFb)
			pWinBin->CntExclFb++;
		if(bFaFb)
			pWinBin->CntFaFb++;
		}
	}

delete m_pInGBSFile;
m_pInGBSFile = NULL;
return(Rslt);
}

int
CGBSmapSNPs::ReportBinCounts(char* pszOutFile,
							char* pszGBSName)				// column name referencing the GBS readset which is to be processed)
{
uint32_t BinIdx;
tsWinBin *pBin;

if(m_pszOutBuffer == NULL)
	{
	m_pszOutBuffer = new uint8_t [cMaxAllocRsltsBuff];
	m_AllocOutBuff = cMaxAllocRsltsBuff;
	}
m_OutBuffIdx = 0;

#ifdef _WIN32
m_hOutFile = open(pszOutFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE));
#else
	if((m_hOutFile = open64(pszOutFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE)) != -1)
		if(ftruncate(m_hOutFile,0)!=0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create/truncate %s - %s",pszOutFile,strerror(errno));
			Reset();
			return(eBSFerrCreateFile);
			}
#endif
if(m_hOutFile < 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create/truncate %s - %s",pszOutFile,strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}
m_OutBuffIdx = sprintf((char *)m_pszOutBuffer,"\"Chrom\",\"Loci\",\"%s:Counts\",\"%s:NA\",\"%s\",\"%s\",\"%s:%s\"\n",pszGBSName,pszGBSName,m_szFndrA,m_szFndrB,m_szFndrA,m_szFndrB);
pBin = m_pWinBins;
for(BinIdx = 0; BinIdx < m_NumWinBins; BinIdx++,pBin++)
	{
	m_OutBuffIdx+=sprintf((char *)&m_pszOutBuffer[m_OutBuffIdx],"%s,%d,%d,%d,%d,%d,%d\n",LocateChrom(pBin->RefChromID),pBin->StartLoci,pBin->LociCounts,pBin->CntNA,pBin->CntExclFa,pBin->CntExclFb,pBin->CntFaFb);
	if((m_OutBuffIdx + 1000) > m_AllocOutBuff)
		{
		if (!CUtility::RetryWrites(m_hOutFile,m_pszOutBuffer,m_OutBuffIdx))
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "ReportBinCounts: Fatal error in RetryWrites()");
			Reset();
			return(eBSFerrFileAccess);
			}
		m_OutBuffIdx = 0;
		}
	}
if(m_OutBuffIdx != 0)
	CUtility::RetryWrites(m_hOutFile,m_pszOutBuffer,m_OutBuffIdx);
m_OutBuffIdx = 0;
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


void
CGBSmapSNPs::InitialiseWIGSpan(void) // initialise WIG span vars to values corresponding to no spans having been previously reported
{
m_WIGChromID = 0;
m_WIGRptdChromID = 0;
m_WIGSpanLoci = 0;
m_WIGSpanLen = 0;
m_WIGRptdSpanLen = 0;
m_WIGSpanCnts = 0;
}

int
CGBSmapSNPs::CompleteWIGSpan(bool bWrite)				// close off any current WIG span ready to start any subsequent span
{
char *pszCurChrom;

	// if existing span then write that span out
if(m_WIGChromID != 0 && m_WIGSpanLen > 0 && m_WIGSpanLoci > 0 && m_WIGSpanCnts > 0)
	{
	// has chrom and/or span changed since previously writing out a span?
	if(m_WIGChromID != m_WIGRptdChromID || m_WIGSpanLen != m_WIGRptdSpanLen)
		{
		pszCurChrom = LocateChrom(m_WIGChromID);
		m_OutBuffIdx += sprintf((char *)&m_pszOutBuffer[m_OutBuffIdx],"variableStep chrom=%s span=%d\n",pszCurChrom,m_WIGSpanLen);
		m_WIGRptdChromID = m_WIGChromID;
		m_WIGRptdSpanLen = m_WIGSpanLen;
		}
	m_OutBuffIdx += sprintf((char *)&m_pszOutBuffer[m_OutBuffIdx],"%d %d\n",m_WIGSpanLoci,(uint32_t)((m_WIGSpanCnts + m_WIGSpanLen-1)/m_WIGSpanLen));
	}
if((bWrite && m_OutBuffIdx) || (m_OutBuffIdx + 1000) >  m_AllocOutBuff)
	{
	if(!CUtility::RetryWrites(m_hOutFile, m_pszOutBuffer, m_OutBuffIdx))
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"RetryWrites() error");
		return(eBSFerrWrite);
		}
	m_OutBuffIdx=0;
	}
m_WIGSpanLoci = 0;
m_WIGSpanLen = 0;
m_WIGSpanCnts = 0;
return(eBSFSuccess);
}

int
CGBSmapSNPs::AccumWIGBinCnts(uint32_t ChromID,	// accumulate wiggle counts into variableStep spans over this chromosome
			uint32_t Loci,		// bin counts starting from this loci - WIG uses 1 based loci instead of the usual UCSC BED 0 based loci
			uint32_t Cnts,		// bin has this many counts attributed
			uint32_t BinLen)	// bin is this length
{
	int Rslt;
	uint32_t Meanx100;
	if(ChromID != m_WIGChromID || Cnts == 0)		// onto a different chromosome, or current span is at maximal length?
	{
		if(m_WIGChromID != 0)
		{
			if((Rslt = CompleteWIGSpan()) < 0)
				return(Rslt);
		}
		if(Cnts > 0 && BinLen > 0)
		{
			m_WIGChromID = ChromID;
			m_WIGSpanLoci = Loci;
			m_WIGSpanLen = BinLen;
			m_WIGSpanCnts = (uint64_t)Cnts * BinLen;
		}
		return(eBSFSuccess);
	}

	if(m_WIGSpanLen == 0 || m_WIGSpanCnts == 0)
	{
		m_WIGSpanLoci = Loci;
		m_WIGSpanLen = BinLen;
		m_WIGSpanCnts = (uint64_t)Cnts * BinLen;
		return(eBSFSuccess);
	}

	Meanx100 = 100 * (uint32_t)(m_WIGSpanCnts / (uint64_t)m_WIGSpanLen);
	if((Cnts <= 5 && (Cnts * 100) != Meanx100) || (Meanx100 < (Cnts * 75) || Meanx100 >= (Cnts * 125)))
	{
	// write current span out
		if((Rslt = CompleteWIGSpan()) < 0)
			return(Rslt);
		m_WIGSpanLoci = Loci;
		m_WIGSpanLen = BinLen;
		m_WIGSpanCnts = (uint64_t)Cnts * BinLen;
		return(eBSFSuccess);
	}
	m_WIGSpanCnts += (uint64_t)Cnts * BinLen;
	m_WIGSpanLen = Loci - m_WIGSpanLoci + BinLen;
	return(eBSFSuccess);
}


int
CGBSmapSNPs::AccumWIGCnts(uint32_t ChromID,	// accumulate wiggle counts into variableStep spans over this chromosome
			uint32_t Loci,		// this loci  - WIG uses 1 based loci instead of the usual UCSC BED 0 based loci
			uint32_t Cnts,		 // has this many counts attributed
			uint32_t MaxSpanLen) // allow WIG spans to be this maximal length
{
	int Rslt;
	uint32_t Meanx100;
	if(ChromID != m_WIGChromID || m_WIGSpanLen >= MaxSpanLen || Cnts == 0)		// onto a different chromosome, or current span is at maximal length?
	{
		if(m_WIGChromID != 0)
		{
			if((Rslt = CompleteWIGSpan()) < 0)
				return(Rslt);
		}
		if(Cnts > 0)
		{
			m_WIGChromID = ChromID;
			m_WIGSpanLoci = Loci;
			m_WIGSpanLen = 1;
			m_WIGSpanCnts = Cnts;
		}
		return(eBSFSuccess);
	}

	if(m_WIGSpanLen == 0 || m_WIGSpanCnts == 0)
	{
		m_WIGSpanLoci = Loci;
		m_WIGSpanLen = 1;
		m_WIGSpanCnts = (uint64_t)Cnts;
		return(eBSFSuccess);
	}

	Meanx100 = 100 * (uint32_t)(m_WIGSpanCnts / (uint64_t)m_WIGSpanLen);
	if((Cnts <= 5 && (Cnts * 100) != Meanx100) || (Meanx100 < (Cnts * 75) || Meanx100 >= (Cnts * 125)))
	{
		// write current span out
		if((Rslt = CompleteWIGSpan()) < 0)
			return(Rslt);
		m_WIGSpanLoci = Loci;
		m_WIGSpanLen = 1;
		m_WIGSpanCnts = Cnts;
		return(eBSFSuccess);
	}
	m_WIGSpanCnts += Cnts;
	m_WIGSpanLen = Loci - m_WIGSpanLoci + 1;
	return(eBSFSuccess);
}



int
CGBSmapSNPs::ReportCountsWIG(char* pszOutFile,
							 bool bFb,						// false if founder Fa, true if founder Fb
							 char* pszGBSName,				// column name referencing the GBS readset which is to be processed
							 char *pszFounder)				// founder
{
uint32_t BinIdx;
tsWinBin *pBin;

if(m_pszOutBuffer == NULL)
	{
	m_pszOutBuffer = new uint8_t [cMaxAllocRsltsBuff];
	m_AllocOutBuff = cMaxAllocRsltsBuff;
	}
m_OutBuffIdx = 0;

#ifdef _WIN32
m_hOutFile = open(pszOutFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE));
#else
	if((m_hOutFile = open64(pszOutFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE)) != -1)
		if(ftruncate(m_hOutFile,0)!=0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create/truncate %s - %s",pszOutFile,strerror(errno));
			Reset();
			return(eBSFerrCreateFile);
			}
#endif
if(m_hOutFile < 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create/truncate %s - %s",pszOutFile,strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}

m_OutBuffIdx = sprintf((char *)m_pszOutBuffer,"track name=\"SNP GBS %s->Fndr %s\" description=\"SNP GBS readset %s aligned against Founder %s\" useScore=1\n", pszGBSName,pszFounder, pszGBSName,pszFounder);
uint32_t CurChromID = 0;
uint32_t CurBinSpan = 0;
uint32_t PrevBinSpan = 0;
InitialiseWIGSpan();
pBin = m_pWinBins;
for(BinIdx = 0; BinIdx < m_NumWinBins; BinIdx++,pBin++)
	{
	CurBinSpan = pBin->EndLoci-pBin->StartLoci+1;
	AccumWIGBinCnts(pBin->RefChromID,pBin->StartLoci + 1,bFb ? pBin->CntExclFb :  pBin->CntExclFa,CurBinSpan); // WIG uses 1 based loci instead of the usual UCSC BED 0 based loci
	}

CompleteWIGSpan(true);

if(m_OutBuffIdx != 0)
	CUtility::RetryWrites(m_hOutFile,m_pszOutBuffer,m_OutBuffIdx);
m_OutBuffIdx = 0;
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

int			// success or otherwise ( >= 0 success, < 0 if processing failed)
CGBSmapSNPs::Process(eModeGBSMapSNPs PMode,			// processing mode
					 int BinSize,					// SNP loci counts are accumulated into this sized (Kbp) non-overlapping bins
					 char* pszGBSName,				// column name referencing the GBS readset which is to be processed
					 char* pszInNMFile,		// processing chromosome name mapping from this file
					 char* pszInGBSFile,	// processing GBS SNP calls from this file
					 char* pszOutFile)		// windowed GBS SNP calls output file
{
int Rslt;
Reset();
m_PMode = PMode;
m_BinSize = BinSize * 1000;	// parameter was in Kbp
m_pszGBSName = pszGBSName;
m_pszInNMFile = pszInNMFile;
m_pszInGBSFile = pszInGBSFile;
m_pszOutFile = pszOutFile;

// parse and load chromosome sizes and name mappings
if((Rslt = LoadNM(pszInNMFile)) < eBSFSuccess)
	{
	Reset();
	return(Rslt);
	}

// parse and load the GBS SNP loci file
if((Rslt = LoadGBSSNPs(pszGBSName,pszInGBSFile)) < eBSFSuccess)
	{
	Reset();
	return(Rslt);
	}

Rslt = ReportBinCounts(pszOutFile,pszGBSName);
char szSuffixName[_MAX_FNAME];
char szWIGFileName[_MAX_FNAME];
sprintf(szSuffixName,".%s_%s.wig",pszGBSName,m_szFndrA);
CUtility::AppendFileNameSuffix(szWIGFileName, m_pszOutFile, (char*)szSuffixName,'.');
Rslt = ReportCountsWIG(szWIGFileName,false,pszGBSName,m_szFndrA);
sprintf(szSuffixName,".%s_%s.wig",pszGBSName,m_szFndrB);
CUtility::AppendFileNameSuffix(szWIGFileName, m_pszOutFile, (char*)szSuffixName,'.');
Rslt = ReportCountsWIG(szWIGFileName,true,pszGBSName,m_szFndrB);

Reset();
return(Rslt);
}

