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

#include "ngskit4b.h"

#include "repassemb.h"



int
Process(etRPPMode PMode,				// processing mode
	char *pszSNPsFile,					// input SNPs file
	char *pszAssembFile,				// assembly to repurpose
	char *pszRepAssembFile);			// write out repurposed assembly to this file

#ifdef _WIN32
int repassemb(int argc, char *argv[])
{
	// determine my process name
	_splitpath (argv[0], NULL, NULL, gszProcName, NULL);
#else
int
repassemb(int argc, char **argv)
{
	// determine my process name
	CUtility::splitpath ((char *)argv[0], NULL, gszProcName);
#endif
	int iFileLogLevel;			// level of file diagnostics
	int iScreenLogLevel;		// level of file diagnostics
	char szLogFile[_MAX_PATH];	// write diagnostics to this file
	int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure
	int NumberOfProcessors;		// number of installed CPUs

	 etRPPMode PMode;				// processing mode
	 char szInFile[_MAX_PATH];		// input fasta assembly file
	 char szSNPsFile[_MAX_PATH];		// file containing kalign called SNPs
	 char szOutFile[_MAX_PATH];		 // output fasta assembly file

	struct arg_lit *help = arg_lit0 ("h", "help", "print this help and exit");
	struct arg_lit *version = arg_lit0 ("v", "version,ver", "print version information and exit");
	struct arg_int *FileLogLevel = arg_int0 ("f", "FileLogLevel", "<int>", "Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
	struct arg_file *LogFile = arg_file0 ("F", "log", "<file>", "diagnostics log file");

	struct arg_int *pmode = arg_int0 ("m", "mode", "<int>", "processing mode: 0 only currently supported processing mode");
	struct arg_file *snpsfile = arg_file1("s","snps","<file>", "input file containing kalign called SNPs");
	struct arg_file *infile = arg_file1("i", "in", "<file>", "input file containing fasta assembly sequences to be repurposed");
	struct arg_file *outfile = arg_file1 ("o", "out", "<file>", "output file into which to write fasta assembly sequences with SNP loci bases replaced by SNP call major allele bases");
	struct arg_end *end = arg_end (200);

	void *argtable[] = { help,version,FileLogLevel,LogFile,
						pmode,snpsfile,infile,outfile,end };

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
		if (!gDiagnostics.Open (szLogFile, (etDiagLevel)iScreenLogLevel, (etDiagLevel)iFileLogLevel, true))
		{
			printf ("\nError: Unable to start diagnostics subsystem\n");
			if (szLogFile[0] != '\0')
				printf (" Most likely cause is that logfile '%s' can't be opened/created\n", szLogFile);
			exit (1);
		}

		gDiagnostics.DiagOut (eDLInfo, gszProcName, "Subprocess %s Version %s starting", gpszSubProcess->pszName, kit4bversion);
		gExperimentID = 0;
		gProcessID = 0;
		gProcessingID = 0;

		PMode = pmode->count ? (etRPPMode)pmode->ival[0] : eRPMdefault;
		if (PMode < eRPMdefault || PMode >= eRPMplaceholder)
		{
			gDiagnostics.DiagOut (eDLFatal, gszProcName, "Error: Processing mode '-m%d' specified outside of range %d..%d\n", PMode, eRPMdefault, eRPMplaceholder-1);
			exit (1);
		}


		// show user current resource limits
#ifndef _WIN32
		gDiagnostics.DiagOut (eDLInfo, gszProcName, "Resources: %s", CUtility::ReportResourceLimits ());
#endif

#ifdef _WIN32
		SYSTEM_INFO SystemInfo;
		GetSystemInfo (&SystemInfo);
		NumberOfProcessors = SystemInfo.dwNumberOfProcessors;
#else
		NumberOfProcessors = sysconf (_SC_NPROCESSORS_CONF);
#endif

		strcpy (szSNPsFile, snpsfile->filename[0]);
		CUtility::TrimQuotedWhitespcExtd (szSNPsFile);
		if (szSNPsFile[0] == '\0')
			{
			gDiagnostics.DiagOut (eDLFatal, gszProcName, "No input file containing kalign generated SNP calls specified");
			exit (1);
			}


		strcpy (szInFile, infile->filename[0]);
		CUtility::TrimQuotedWhitespcExtd (szInFile);
		if (szInFile[0] == '\0')
			{
			gDiagnostics.DiagOut (eDLFatal, gszProcName, "No input file containing fasta assembly sequences specified");
			exit (1);
			}

		strcpy (szOutFile, outfile->filename[0]);
		CUtility::TrimQuotedWhitespcExtd (szOutFile);
		if (szOutFile[0] == '\0')
			{
			gDiagnostics.DiagOut (eDLFatal, gszProcName, "No output file into which to write fasta assembly sequences with SNP loci bases replaced by SNP call major allele bases");
			exit (1);
			}

		gDiagnostics.DiagOut (eDLInfo, gszProcName, "Processing parameters:");
		const char *pszDescr;
		switch (PMode) {
			case eRPMdefault:
				pszDescr = "Default SNP loci base replacements";
				break;
		}



		gDiagnostics.DiagOutMsgOnly (eDLInfo, "Processing : '%s'", pszDescr);
		gDiagnostics.DiagOutMsgOnly (eDLInfo, "SNP file : '%s'", szSNPsFile);
		gDiagnostics.DiagOutMsgOnly (eDLInfo, "Input file : '%s'", szInFile);
		gDiagnostics.DiagOutMsgOnly (eDLInfo, "Output file : '%s'", szOutFile);


#ifdef _WIN32
		SetPriorityClass (GetCurrentProcess (), BELOW_NORMAL_PRIORITY_CLASS);
#endif
		gStopWatch.Start ();
		Rslt = 0;
		Rslt = Process (PMode,			// processing mode
						szSNPsFile,		// kalign called SNPs
						szInFile,		// input fasta or SAM file
						szOutFile);		// output to this file
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
	return 0;
}

int
Process(etRPPMode PMode,				// processing mode
		char* pszSNPsFile,					// input SNPs file
		char* pszAssembFile,				// assembly to repurpose
		char* pszRepAssembFile)				// write out repurposed assembly to this file
{
int Rslt;
CRepAssemb *pRepAssemb;
if((pRepAssemb = new CRepAssemb) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to instantiate instance of CRepAssemb");
	return(eBSFerrInternal);
	}

Rslt = pRepAssemb->ProccessAssembly(pszSNPsFile,			// input kalign generated SNPs file
									pszAssembFile,			// input fasta assembly file
									pszRepAssembFile);		// output to this fasta assembly file

if(pRepAssemb != NULL)
	delete pRepAssemb;
return(Rslt);
}

CRepAssemb::CRepAssemb(void)			// constructor
{
m_pSNPSites = NULL;
m_pCSV = NULL;
m_pSeqBuffer = NULL;
m_hInFile = -1;
m_hOutFile = -1;
Reset();
}

CRepAssemb::~CRepAssemb(void)			// destructor
{
if(m_pCSV != NULL)
	delete m_pCSV;
if(m_pSNPSites != NULL)
	{
#ifdef _WIN32
	free(m_pSNPSites);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pSNPSites != MAP_FAILED)
		munmap(m_pSNPSites, m_AllocdSNPSitesMem);
#endif
	}
if(m_pSeqBuffer != NULL)
	{
#ifdef _WIN32
	free(m_pSeqBuffer);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pSeqBuffer != MAP_FAILED)
		munmap(m_pSeqBuffer, m_AllocSeqBuffMem);
#endif
	}

if(m_hInFile != -1)
	close(m_hInFile);
if(m_hOutFile != -1)
	close(m_hOutFile);
}

void
CRepAssemb::Reset(void)					// resets instance state back to that immediately following instantiation
{
if(m_pCSV != NULL)
	{
	delete m_pCSV;
	m_pCSV = NULL;
	}
if(m_pSNPSites != NULL)
	{
#ifdef _WIN32
	free(m_pSNPSites);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pSNPSites != MAP_FAILED)
		munmap(m_pSNPSites, m_AllocdSNPSitesMem);
#endif
	m_pSNPSites = NULL;
	}

if(m_pSeqBuffer != NULL)
	{
#ifdef _WIN32
	free(m_pSeqBuffer);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pSeqBuffer != MAP_FAILED)
		munmap(m_pSeqBuffer, m_AllocSeqBuffMem);
#endif
	m_pSeqBuffer = NULL;
	}

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

m_NumSNPSites = 0;
m_NumAllocdSNPSites = 0;
m_AllocdSNPSitesMem = 0;

m_SeqBuffIdx = 0;
m_AllocSeqBuffMem = 0;

m_szTargAssemblyName[0] = '\0';
m_LAChromNameID = 0;
m_NumChromNames = 0;
m_NxtszChromIdx=0;
m_szChromNames[0] = '\0';
m_szChromIdx[0] = 0;

m_MinCoverage = cDfltMinSiteCoverage;
m_MinAlleleProp = cDfltMinAlleleProp;
m_PValueThres = cDfltPValueThres;
}

int						// eBSFSuccess or error
CRepAssemb::LoadKalignSNPs(char *pszSNPFile)	// load and parse kalign generated SNPs
{
int Rslt;
int NumFields;
UINT32 EstNumSNPs;
UINT32 NumSNPsAccepted;
UINT32 NumSNPsParsed;
UINT32 RowNumber;
uint32_t ReadsetSiteId;

int ExpNumFields;
uint32_t ChromID;
uint32_t Loci;
etSeqBase RefBase;
etSeqBase AlleleBase;
uint32_t CoveringBases;
uint32_t TotMismatches;
uint32_t BaseCnts[5];
char szChrom[cMaxDatasetSpeciesChrom+1];

int CntRef;
double PValue;
double LocalSeqErrRate;

tsSNPSite *pSNPSite;

if ((m_pCSV = new CCSVFile) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to instantiate CCSVfile");
	Reset();
	return(eBSFerrObj);
	}
m_pCSV->SetMaxFields(cAlignNumSSNPXfields);

if((Rslt = m_pCSV->Open(pszSNPFile)) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to open file: '%s'", pszSNPFile);
	Reset();
	return(Rslt);
	}
gDiagnostics.DiagOut(eDLInfo, gszProcName, "Loading SNP calls from file: '%s'", pszSNPFile);
EstNumSNPs = m_pCSV->EstNumRows();
NumSNPsParsed = 0;
NumSNPsAccepted = 0;
RowNumber = 0;
ExpNumFields = 0;
szChrom[0] = '\0';
while ((Rslt = m_pCSV->NextLine()) > 0)				// onto next line containing fields
	{
	RowNumber += 1;
	NumFields = m_pCSV->GetCurFields();
	if (ExpNumFields && ExpNumFields != NumFields)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Inconsistency in number of fields, previously %d but %d fields parsed from '%s' near line %d", ExpNumFields, NumFields, pszSNPFile, RowNumber);
		Reset();
		return(eBSFerrFieldCnt);
		}
	if (!ExpNumFields)
		{
		if (!(NumFields == cAlignNumSSNPfields || NumFields == cAlignNumSSNPXfields))								// must be exactly this many if 'ngskit4b kalign' format
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "%d fields parsed from '%s' near line %d, expected %d, was file generated by 'kit4b kalign'?", NumFields, pszSNPFile, RowNumber, cAlignNumSSNPfields);
			Reset();
			return(eBSFerrFieldCnt);
			}
		ExpNumFields = NumFields;
		}

	if (RowNumber == 1)
		{
		if (m_pCSV->IsLikelyHeaderLine())
			continue;
		else
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Expected CSV file '%s' first line to be a header line with  fully quoted field names", pszSNPFile);
			Reset();
			return(eBSFerrFieldCnt);
			}
		}

	if (!NumSNPsParsed)
		{
		// parse out target assembly name against which alignments were made and the chromosome name
		char* pSrc;
		char* pDst;
		char Chr;
		int Len;

		m_pCSV->GetText(3, &pSrc);
		pDst = m_szTargAssemblyName;
		Len = 0;
		while (Len < sizeof(m_szTargAssemblyName) - 1 && (Chr = *pSrc++))
			{
			*pDst++ = Chr;
			Len++;
			*pDst = '\0';
			}
		}

	char* pszTxt;
	int SNPlen;
	m_pCSV->GetInt(7, &SNPlen);			// check that we are processing SNPs!
	if (SNPlen != 1)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Expected SNP CSV file '%s' to only contain 1 base SNPs, 'Len' = %d near line %d", pszSNPFile, SNPlen, RowNumber);
		Reset();
		return(eBSFerrFieldCnt);
		}
	m_pCSV->GetUint(11, &CoveringBases);		// get "Bases"
	if(CoveringBases < m_MinCoverage)
		continue;

	NumSNPsParsed++;
	m_pCSV->GetText(4, &pszTxt);			// get "Chrom"
	strncpy(szChrom, pszTxt, sizeof(szChrom));
	szChrom[sizeof(szChrom) - 1] = '\0';
	if(m_szCurChrom[0] == 0 || m_CurChromID == 0 || stricmp(szChrom, m_szCurChrom))
		{
		if((Rslt = ChromID = AddChrom(szChrom)) < 1)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Failed generating a unique chromosome name identifier '%s'", szChrom);
			Reset();
			return(Rslt);
			}
		strncpy(m_szCurChrom, szChrom, sizeof(m_szCurChrom));
		m_CurChromID = ChromID;
		}
	m_pCSV->GetUint(1, &ReadsetSiteId);
	m_pCSV->GetUint(5, &Loci);			// get "StartLoci"
	m_pCSV->GetUint(12, &TotMismatches);		// get "Mismatches"
	CntRef = CoveringBases - TotMismatches;
	m_pCSV->GetUint(14, &BaseCnts[0]);	// get "MMBaseA"
	m_pCSV->GetUint(15, &BaseCnts[1]);	// get "MMBaseC"
	m_pCSV->GetUint(16, &BaseCnts[2]);	// get "MMBaseG"
	m_pCSV->GetUint(17, &BaseCnts[3]);	// get "MMBaseT"
	m_pCSV->GetUint(18, &BaseCnts[4]);	// get "MMBaseN"
	m_pCSV->GetText(13, &pszTxt);		// get "RefBase"
	switch (*pszTxt) {
		case 'a': case 'A':
			RefBase = eBaseA;
			break;
		case 'c': case 'C':
			RefBase = eBaseC;
			break;
		case 'g': case 'G':
			RefBase = eBaseG;
			break;
		case 't': case 'T': case 'u': case 'U':	// U in case RNA alignments..
			RefBase = eBaseT;
			break;
		case 'n': case 'N':				// unlikely to have a SNP against an indeterminate base but you never know...
			RefBase = eBaseN;
			break;
		default:
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Expected SNP RefBase ('%s') to be one of 'ACGTN' in CSV file '%s' near line %d", pszTxt, pszSNPFile, RowNumber);
			Reset();
			return(eBSFerrFieldCnt);
		}
	m_pCSV->GetDouble(19, &LocalSeqErrRate);	// get "BackgroundSubRate"

	if(m_MinAlleleProp > 0.0 && (CoveringBases == 0 || (TotMismatches / (double)CoveringBases) < m_MinAlleleProp))
		continue;

	PValue = 1.0 - m_Stats.Binomial(CoveringBases, TotMismatches, LocalSeqErrRate);
	if(PValue > m_PValueThres)
		continue;

	AlleleBase = eBaseN;
	for (int Idx = 0; Idx < 4; Idx++)
		{
		if(Idx == RefBase)
			continue;
		if(m_MinAlleleProp > 0.0 && (BaseCnts[Idx] == 0 || (((double)BaseCnts[Idx] * 4) / (double)CoveringBases) < m_MinAlleleProp))
			continue;

		PValue = 1.0 - m_Stats.Binomial(CoveringBases, BaseCnts[Idx], LocalSeqErrRate);
		if(PValue > m_PValueThres)
			continue;

		AlleleBase = Idx;
		}
	if(AlleleBase == eBaseN)
		continue;
	NumSNPsAccepted++;

	if(m_pSNPSites == NULL)
		{
		size_t memreq = (size_t)(sizeof(tsSNPSite) * cAllocNumSNPs);

#ifdef _WIN32
		m_pSNPSites = (tsSNPSite*)malloc(memreq);	// initial and perhaps the only allocation
		if(m_pSNPSites == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "SNPs: Initial memory allocation of %lld bytes - %s", (INT64)memreq, strerror(errno));
			Reset();
			return(eBSFerrMem);
			}
#else
		// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
		m_pSNPSites = (tsSNPSite *)mmap(NULL, memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
		if(m_pSNPSites == MAP_FAILED)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "SNPs: Memory allocation of %lld bytes through mmap()  failed - %s", (INT64)memreq, strerror(errno));
			m_pSNPSites = NULL;
			Reset();
			return(eBSFerrMem);
			}
#endif
		m_AllocdSNPSitesMem = memreq;
		m_NumSNPSites = 0;
		m_NumAllocdSNPSites = cAllocNumSNPs;
		}
	else
		{
		if((m_NumSNPSites + 100) > m_NumAllocdSNPSites)
			{
			size_t memreq = m_AllocdSNPSitesMem + (cAllocNumSNPs * sizeof(tsSNPSite));
				
#ifdef _WIN32
			pSNPSite = (tsSNPSite*)realloc(m_pSNPSites, memreq);
#else
			pSNPSite = (tsSNPSite*)mremap(m_pSNPSites, m_AllocdSNPSitesMem, memreq, MREMAP_MAYMOVE);
			if(pSNPSite == MAP_FAILED)
				pSNPSite = NULL;
#endif
			if(pSNPSite == NULL)
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "SNPs: Memory re-allocation to %lld bytes - %s", memreq, strerror(errno));
				return(eBSFerrMem);
				}
			m_pSNPSites = pSNPSite;
			m_NumAllocdSNPSites += cAllocNumSNPs;
			m_AllocdSNPSitesMem = memreq;
			}
		}
	pSNPSite = &m_pSNPSites[m_NumSNPSites++];
	pSNPSite->AlleleBase = AlleleBase;
	pSNPSite->ChromID = ChromID;
	pSNPSite->Loci = Loci;
	pSNPSite->RefBase = RefBase;
	}
if(m_NumSNPSites > 1)
	qsort(m_pSNPSites, m_NumSNPSites, sizeof(tsSNPSite), SortSNPChromLoci);

tsChromSites *pChromSites;
uint32_t SNPIdx;
uint32_t PrevChromIdx;

m_NumChromSites = 0;
if(m_NumSNPSites)
	{
	pSNPSite = m_pSNPSites;
	PrevChromIdx = 0;
	pChromSites = &m_ChromSites[0];
	pChromSites->ChromID = pSNPSite->ChromID;
	pChromSites->NumSNPSites=1;
	pChromSites->SNPSiteIdx=0;
	m_NumChromSites = 1;
	pSNPSite+=1;
	for(SNPIdx = 1; SNPIdx < m_NumSNPSites; SNPIdx++,pSNPSite++)
		{
		if(pChromSites->ChromID != pSNPSite->ChromID)
			{
			pChromSites+=1;
			pChromSites->ChromID = pSNPSite->ChromID;
			pChromSites->SNPSiteIdx=SNPIdx;
			pChromSites->NumSNPSites=1;
			m_NumChromSites+=1;
			}
		else
			pChromSites->NumSNPSites+=1;
		}
	}
gDiagnostics.DiagOut(eDLInfo, gszProcName, "SNPs: Parsed %u called SNPs, accepted %u as major allele SNPs on %u chromosomes", NumSNPsParsed,m_NumSNPSites,m_NumChromSites);
return(eBSFSuccess);
}

int
CRepAssemb::ProccessAssembly(char *pszSNPFile,	// input kalign generated SNPs file
				char* pszInFile,		// input fasta assembly file
				char* pszOutFile)		// output to this fasta assembly file
{
int Rslt;
bool bNL;
bool bDescr;
int NumRead;
uint8_t InChr;
uint8_t *pInChr;
size_t FastaIdx;
uint8_t FastaChr;
uint8_t *pszChromName;
uint8_t szChromName[cMaxDatasetSpeciesChrom+1];
uint8_t *pChromNameChr;
size_t memreq;
int BasesReplaced;
Reset();
if((Rslt = LoadKalignSNPs(pszSNPFile)) < 0)
	{
	Reset();
	return(Rslt);
	}

#ifdef _WIN32
m_hInFile = open(pszInFile, O_READSEQ );		// file access is normally sequential..
#else
m_hInFile = open64(pszInFile, O_READSEQ );		// file access is normally sequential..
#endif
if(m_hInFile == -1)							// check if file open succeeded
	{
	gDiagnostics.DiagOut (eDLFatal, gszProcName, "Unable to open input file '%s' : %s",pszInFile,strerror(errno));
	return(eBSFerrOpnFile);
	}

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

// an initial allocation for buffering input assembly, will be realloc'd if required so as to contain complete assembly
memreq = (size_t)cAllocAssembFastaMem;

#ifdef _WIN32
m_pSeqBuffer = (uint8_t *)malloc(memreq);	// initial and perhaps the only allocation
if(m_pSeqBuffer == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Initial memory allocation of %lld bytes - %s", (INT64)memreq, strerror(errno));
	Reset();
	return(eBSFerrMem);
	}
#else
// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
m_pSeqBuffer = (uint8_t *)mmap(NULL, memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
if(m_pSeqBuffer == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Assembly: Memory allocation of %lld bytes through mmap()  failed - %s", (INT64)memreq, strerror(errno));
	m_pSeqBuffer = NULL;
	Reset();
	return(eBSFerrMem);
	}
#endif

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading fasta assembly from : '%s'",pszOutFile);
m_AllocSeqBuffMem = memreq;
m_SeqBuffIdx = 0;

while((NumRead = (int)read(m_hInFile, &m_pSeqBuffer[m_SeqBuffIdx], (int)(m_AllocSeqBuffMem - m_SeqBuffIdx))) > 0)
	{
	if((NumRead + m_SeqBuffIdx + (cAllocAssembFastaMem/100)) >= m_AllocSeqBuffMem)
		{
		memreq = m_AllocSeqBuffMem + cAllocAssembFastaMem;
#ifdef _WIN32
		pInChr = (uint8_t *)realloc(m_pSeqBuffer, memreq);
#else
		pInChr = (uint8_t*)mremap(m_pSeqBuffer, m_AllocSeqBuffMem, memreq, MREMAP_MAYMOVE);
		if(pInChr == MAP_FAILED)
			pInChr = NULL;
#endif
		if(pInChr == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Memory re-allocation to %lld bytes - %s", memreq, strerror(errno));
			return(eBSFerrMem);
			}
		m_pSeqBuffer = pInChr;
		m_AllocSeqBuffMem = memreq;
		}
	m_SeqBuffIdx += NumRead;
	}

if(NumRead < 0)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Loading assembly '%s' failed - %s", pszOutFile,strerror(errno));
	Reset();
	return(eBSFerrFileAccess);
	}
BasesReplaced = 0;
gDiagnostics.DiagOut(eDLFatal, gszProcName, "Replacing assembly sequence major allele sites ...");
if(m_NumChromSites)
	{
	bNL = true; // 1st line in fasta may be a descriptor with no preceding NL so pretend there was a preceding NL
	bDescr = false;	// only true if sequence name parsed out from descriptor line
	pInChr = m_pSeqBuffer;
	for(FastaIdx = 0; FastaIdx < m_SeqBuffIdx; FastaIdx++,pInChr++)
		{
		if((InChr = *pInChr) == '\0')
			break;
		switch(InChr) {
			case '\t': case '\n': case '\r':	// whitespace or line endings are accepted
				break;
			default:
				if(InChr >= 0x20 && InChr <= 0x7f)
					break;
				gDiagnostics.DiagOut (eDLFatal, gszProcName, "Expected Fasta '%s' to only contain ascii chars",pszInFile);
				Reset();
				return(eBSFerrFastqChr);
			}

		if(InChr == '\n')
			{
			if(bDescr)
				BasesReplaced += ReplaceBases((char *)szChromName,(uint32_t)min(0x0ffffffff,m_SeqBuffIdx-FastaIdx),(char *)pInChr+1);
			bDescr = false;
			bNL = true;
			}
		else
			{
			if(bNL == true && InChr == '>')		// start of descriptor, sequence name is terminated by whitespace
				{
				// extract sequence name
				pszChromName = szChromName;
				pChromNameChr = pInChr+1;
				while((FastaChr = *pChromNameChr++) != ' ' && FastaChr != '\t' && FastaChr != '\r' && FastaChr != '\n')
					*pszChromName++ = FastaChr;
				*pszChromName = '\0';
				bDescr = true;
				}
			bNL = false;
			}
		}
	}

gDiagnostics.DiagOut(eDLInfo, gszProcName, "%d replacements completed, writing repurposed assembly to '%s'",BasesReplaced,pszOutFile);

if(m_SeqBuffIdx)
	CUtility::RetryWrites(m_hOutFile, m_pSeqBuffer, m_SeqBuffIdx);	

// close input file
close(m_hInFile);
m_hInFile = -1;
// commit and close output file
#ifdef _WIN32
_commit(m_hOutFile);
#else
fsync(m_hOutFile);
#endif
close(m_hOutFile);
m_hOutFile = -1;
Reset();
return(eBSFSuccess);
}

int		// returned chrom identifier, < 1 if unable to accept this chromosome name
CRepAssemb::AddChrom(char* pszChrom) // associate unique identifier with this chromosome name
{
int ChromNameIdx;
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
	return(eBSFerrMaxEntries);
if(m_NumChromNames == cMaxSNPChromNames)
	return(eBSFerrMaxEntries);

m_szChromIdx[m_NumChromNames++] = m_NxtszChromIdx;
strcpy(&m_szChromNames[m_NxtszChromIdx], pszChrom);
m_NxtszChromIdx += ChromNameLen + 1;
m_LAChromNameID = m_NumChromNames;
return(m_LAChromNameID);
}

char* 
CRepAssemb::LocateChrom(int ChromID)
{
if(ChromID < 1 || ChromID > m_NumChromNames)
	return(NULL);
return(&m_szChromNames[m_szChromIdx[ChromID-1]]);
}

tsChromSites * // returned ptr
CRepAssemb::LocateChromSites(char* pszChrom)	// sites are to be those on this chromosome
{
tsChromSites *pChromSite;
char *pszChromSite;
int ChromIdx;
pChromSite = m_ChromSites;
for(ChromIdx = 0; ChromIdx < m_NumChromSites; ChromIdx++,pChromSite++)
	{
	if((pszChromSite = LocateChrom(pChromSite->ChromID))==NULL)
		break;
	if(!strcmp(pszChrom,pszChromSite))
		return(pChromSite);
	}
return(NULL);
}

int									// number of bases replaced
CRepAssemb::ReplaceBases(char *pszChrom, // chromosome name
			 uint32_t SeqLen,			// max of this number of bases in fasta sequence length - includes any whitespace
			 char *pszSeq)			// sequence
{
uint32_t RefLoci;
tsChromSites *pChromSites;
tsSNPSite *pSNPSite;
uint32_t NumSNPBases;

if(SeqLen == 0 || pszChrom == NULL || pszChrom[0] == '\0' || pszSeq == NULL || pszSeq[0] == '\0')
	return(0);

// locate ChromSites for name
if((pChromSites = LocateChromSites(pszChrom))==NULL)
	return(0);

NumSNPBases = 0;
RefLoci = 0;
pSNPSite = &m_pSNPSites[pChromSites->SNPSiteIdx]; 
while(SeqLen-- && NumSNPBases != pChromSites->NumSNPSites)
	{
	switch(*pszSeq) {
		case 'a': case 'A':
		case 'c': case 'C':
		case 'g': case 'G':
		case 't': case 'T':
		case 'n': case 'N':
			break;
		case '>':
		case '\0':
			return(NumSNPBases);
		default:						// skipping over any whitespace or non-base character
			pszSeq+=1;
			continue;
		}
	
	if(RefLoci++ == pSNPSite->Loci)
		{
		switch(pSNPSite->AlleleBase){
			case eBaseA:
				*pszSeq = 'a';
				break;
			case eBaseC:
				*pszSeq = 'c';
				break;
			case eBaseG:
				*pszSeq = 'g';
				break;
			case eBaseT:
				*pszSeq = 't';
				break;
			default:
				*pszSeq = 'N';
				break;
			}
		NumSNPBases++;
		pSNPSite++;
		}
	pszSeq+=1;
	}

return(NumSNPBases);
}

// SortSNPChromIDLoci
// Sort m_pSNPSites by ascending ChromID.Loci
int
CRepAssemb::SortSNPChromLoci(const void* arg1, const void* arg2)
{
tsSNPSite* pEl1 = (tsSNPSite*)arg1;
tsSNPSite* pEl2 = (tsSNPSite*)arg2;

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

