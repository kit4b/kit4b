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
Process(etRPPMode PMode,	// processing mode
	int MinSNPreads,				// must be at least this number of reads covering any loci before processing for SNPs at this loci
	double SNPMajorAllele,			// SNP major allele must be at least this proportion of allele counts at SNP site
	double PValueThres,				// SNP PValue must be <= this threshold
	char *pszProcFile,		// input kalign generated SNPs or alignment segments file name, filename can be wildcardeds
	char *pszAssembFile,	// input assembly to repurpose
	char *pszRepAssembFile);	// write out repurposed assembly to this file

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
	int MinSNPreads;				// must be at least this number of reads covering any loci before processing for SNPs at this loci
	double SNPMajorAllele;			// SNP major allele must be at least this proportion of allele counts at SNP site
	double PValueThres;				// SNP PValue must be <= this threshold
	char szInFile[_MAX_PATH];		// input fasta assembly file
	char szSNPsFile[_MAX_PATH];		// file containing kalign called SNPs
	char szOutFile[_MAX_PATH];		 // output fasta assembly file

	struct arg_lit *help = arg_lit0 ("h", "help", "print this help and exit");
	struct arg_lit *version = arg_lit0 ("v", "version,ver", "print version information and exit");
	struct arg_int *FileLogLevel = arg_int0 ("f", "FileLogLevel", "<int>", "Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
	struct arg_file *LogFile = arg_file0 ("F", "log", "<file>", "diagnostics log file");

	struct arg_int *pmode = arg_int0 ("m", "mode", "<int>", "processing mode: 0 SNP call major allele base replacement, 1 non-alignment segment bases replaced by 'N'");
	struct arg_int *minsnpreads = arg_int0("d","snpreadsmin","<int>","filter out SNP loci having less than this read coverage (default is 20)");
	struct arg_dbl *snpmajorallele = arg_dbl0("P","snpmajorallele","<dbl>", "SNP major allele must be at least this proportion of all counts at SNP site (defaults to 0.5, range 0.25 to 1.0)");
	struct arg_dbl *pvaluethres = arg_dbl0("p","pvaluethres","<dbl>", "SNP PValue must be <= this threshold (defaults to 0.05, range 0.0 to 0.05)");
	struct arg_file *snpsfile = arg_file1("i","ref","<file>", "input file containing kalign called SNPs - use wild cards if multiple SNP files to be processed, or aligned segments BED");
	struct arg_file *infile = arg_file1("I", "in", "<file>", "input file containing fasta assembly sequences to be repurposed");
	struct arg_file *outfile = arg_file1 ("o", "out", "<file>", "output file into which to write repurposed fasta assembly sequences");
	struct arg_end *end = arg_end (200);

	void *argtable[] = { help,version,FileLogLevel,LogFile,
						pmode,minsnpreads,snpmajorallele,pvaluethres,snpsfile,infile,outfile,end };

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

	if(PMode == eRPMdefault)
		{
		MinSNPreads = minsnpreads->count ? minsnpreads->ival[0] : cDfltMinSiteCoverage;
		if(MinSNPreads < 5 || MinSNPreads > 100)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Minimum read coverage at any loci '-d%d' must be in range %d..%d\n",MinSNPreads,5,100);
			exit(1);
			}

		SNPMajorAllele = snpmajorallele->count ? snpmajorallele->dval[0] : cDfltMinAlleleProp;
		if(SNPMajorAllele < 0.25 || SNPMajorAllele > 1.0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: SNP major allele proportion  '-P%f' must be in range %f..%f\n",SNPMajorAllele,0.25,1.0);
			exit(1);
			}

		PValueThres = pvaluethres->count ? pvaluethres->dval[0] : cDfltPValueThres;
		if(PValueThres < 0.0 || PValueThres > 0.10)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: PValue threshold  '-p%f' must be in range %f..%f\n",SNPMajorAllele,0.0,0.05);
			exit(1);
			}
		}
	else
		{
		MinSNPreads = 0;
		SNPMajorAllele = 0.0;
		PValueThres = 0.0;
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
			gDiagnostics.DiagOut (eDLFatal, gszProcName, "No input file containing kalign generated calls specified");
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
			gDiagnostics.DiagOut (eDLFatal, gszProcName, "No output file into which to write repurposed fasta assembly sequences");
			exit (1);
			}

		gDiagnostics.DiagOut (eDLInfo, gszProcName, "Processing parameters:");
		const char *pszDescr;
		switch (PMode) {
			case eRPMdefault:
				pszDescr = "SNP loci major allele base replacements";
				break;
			case eRPMNonalignSegs:
				pszDescr = "Non-aligned segment base replacements";
				break;
		}



		gDiagnostics.DiagOutMsgOnly (eDLInfo, "Processing : '%s'", pszDescr);
		
		switch (PMode) {
			case eRPMdefault:
				gDiagnostics.DiagOutMsgOnly (eDLInfo, "SNP minimum read coverage : %d", MinSNPreads);
				gDiagnostics.DiagOutMsgOnly (eDLInfo, "SNP major allele proportion : %f", SNPMajorAllele);
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"SNP PValue threshold : %f",PValueThres);
				gDiagnostics.DiagOutMsgOnly (eDLInfo, "SNPs file : '%s'", szSNPsFile);
				break;

			case eRPMNonalignSegs:
				gDiagnostics.DiagOutMsgOnly (eDLInfo, "BED file : '%s'", szSNPsFile);
				break;
			}

		gDiagnostics.DiagOutMsgOnly (eDLInfo, "Input file : '%s'", szInFile);
		gDiagnostics.DiagOutMsgOnly (eDLInfo, "Output file : '%s'", szOutFile);


#ifdef _WIN32
		SetPriorityClass (GetCurrentProcess (), BELOW_NORMAL_PRIORITY_CLASS);
#endif
		gStopWatch.Start ();
		Rslt = 0;
		Rslt = Process(PMode,				// processing mode
						MinSNPreads,		// must be at least this number of reads covering any loci before processing for SNPs at this loci
						SNPMajorAllele,		// SNP major allele must be at least this proportion of allele counts at SNP site
						PValueThres,		// SNP PValue must be <= this threshold
						szSNPsFile,			// kalign called SNPs or non-aligned BED file
						szInFile,			// input fasta or SAM file
						szOutFile);			// output to this file
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
		int MinSNPreads,				// must be at least this number of reads covering any loci before processing for SNPs at this loci
		double SNPMajorAllele,			// SNP major allele must be at least this proportion of allele counts at SNP site
		double PValueThres,				// SNP PValue must be <= this threshold
		char *pszProcFile,		// input kalign generated SNPs or alignment segments file name, filename can be wildcarded
		char *pszAssembFile,	// input assembly to repurpose
		char *pszRepAssembFile)	// write out repurposed assembly to this file
{
int Rslt;
CRepAssemb *pRepAssemb;
if((pRepAssemb = new CRepAssemb) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to instantiate instance of CRepAssemb");
	return(eBSFerrInternal);
	}

Rslt = pRepAssemb->ProccessAssembly(PMode,	// processing mode
						MinSNPreads,			// must be at least this number of reads covering any loci before processing for SNPs at this loci
						SNPMajorAllele,			// SNP major allele must be at least this proportion of allele counts at SNP site
						PValueThres,			// SNP PValue must be <= this threshold
						pszProcFile,			// input kalign generated SNPs or alignment segments file name, filename can be wildcarded
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
m_pASegments = NULL;
m_hInFile = -1;
m_hOutFile = -1;
m_hOutRepSNPsFile = -1;
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
if(m_pASegments != NULL)
	{
#ifdef _WIN32
	free(m_pASegments);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pASegments != MAP_FAILED)
		munmap(m_pASegments, m_AllocdASegmentsMem);
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
if(m_hOutRepSNPsFile != -1)
	close(m_hOutRepSNPsFile);
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

if(m_pASegments != NULL)
	{
#ifdef _WIN32
	free(m_pASegments);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pASegments != MAP_FAILED)
		munmap(m_pASegments, m_AllocdASegmentsMem);
#endif
	m_pASegments = NULL;
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
if(m_hOutRepSNPsFile != -1)
	{
	close(m_hOutRepSNPsFile);
	m_hOutRepSNPsFile = -1;
	}

m_NumSortedSites = 0;
m_NumSNPSites = 0;
m_NumAllocdSNPSites = 0;
m_AllocdSNPSitesMem = 0;

m_MaxASegmentLen = 0;
m_NumASegments = 0;
m_NumAllocdASegments = 0;
m_AllocdASegmentsMem = 0;

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
CRepAssemb::LoadKalignSNPs(char *pszSNPFiles)	// load and parse kalign generated SNPs from this - can be wild-carded
{
int Rslt;
char *pszSNPFile;
int NumFields;
uint32_t EstNumSNPs;
uint32_t NumSNPsAccepted;
uint32_t NumSNPsParsed;
uint32_t RowNumber;
uint32_t ReadsetSiteId;

int ExpNumFields;
uint32_t ChromID;
uint32_t Loci;
etSeqBase RefBase;
double PValue;
etSeqBase AlleleBase;
uint32_t CoveringBases;
uint32_t TotMismatches;
uint32_t BaseCnts[5];
char szChrom[cMaxDatasetSpeciesChrom+1];

int CntRef;

CSimpleGlob glob(SG_GLOB_FULLSORT);
glob.Init();
if(glob.Add(pszSNPFiles) < SG_SUCCESS)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to glob '%s",pszSNPFiles);
	return(eBSFerrOpnFile);	// treat as though unable to open file
	}

tsSNPSite *pSNPSite;
m_NumSortedSites = 0;
m_NumSNPFiles = glob.FileCount();

for(int FileIdx = 0;FileIdx < m_NumSNPFiles; FileIdx++)
	{
	pszSNPFile = glob.File(FileIdx);;
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
		m_pCSV->GetDouble(10, &PValue);			// get "PValue"
		m_pCSV->GetUint(12, &TotMismatches);	// get "Mismatches"
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

		if((uint32_t)CntRef > TotMismatches)	// can't have a major allele if reference count >= mismatches
			continue;

		if(PValue > m_PValueThres || CoveringBases < m_MinCoverage)
			continue;

		if((((double)TotMismatches+1)/CoveringBases) < m_MinAlleleProp)
			continue;

		AlleleBase = eBaseN;
		uint32_t AlleleBaseCnts = 0;
		for (int Idx = 0; Idx < 4; Idx++)
			{
			if(BaseCnts[Idx] == 0 || (((double)BaseCnts[Idx]+1) / (double)CoveringBases) < m_MinAlleleProp)
				continue;

			if(BaseCnts[Idx] > AlleleBaseCnts)	// major allele will be that with highest counts!
				{
				AlleleBaseCnts = BaseCnts[Idx];
				AlleleBase = Idx;
				}
			}
		if(AlleleBase == eBaseN || AlleleBase == RefBase)
			continue;

		NumSNPsAccepted++;

		if((pSNPSite = LocateSite(ChromID, Loci)) != NULL)
			{
			if(pSNPSite->ChromID != ChromID || pSNPSite->Loci != Loci)
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "SNPs: Inconsistency in LocateSite()");
				
			if(pSNPSite->RefBase != RefBase)	// double check for consistency in reference bases - same chrom and loci expected to have same reference base!
				{
				char SiteRefBaseChr = CFasta::Base2Chr(pSNPSite->RefBase);
				char RefBaseChr = CFasta::Base2Chr(RefBase);
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "SNPs: Inconsistency in reference base, previously was '%c', in this readset called as being '%c' in CSV file '%s' near line %d", SiteRefBaseChr,RefBaseChr, pszSNPFile, RowNumber);
				
				Reset();
				return(eBSFErrBase);
				}
			pSNPSite->AlleleBaseCnts[AlleleBase]++;
			continue;
			}

		if(m_pSNPSites == NULL)
			{
			size_t memreq = (size_t)(sizeof(tsSNPSite) * cAllocNumSNPs);

	#ifdef _WIN32
			m_pSNPSites = (tsSNPSite*)malloc(memreq);	// initial and perhaps the only allocation
			if(m_pSNPSites == NULL)
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "SNPs: Initial memory allocation of %I64d bytes - %s", (int64_t)memreq, strerror(errno));
				Reset();
				return(eBSFerrMem);
				}
	#else
			// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
			m_pSNPSites = (tsSNPSite *)mmap(NULL, memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
			if(m_pSNPSites == MAP_FAILED)
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "SNPs: Memory allocation of %I64d bytes through mmap()  failed - %s", (int64_t)memreq, strerror(errno));
				m_pSNPSites = NULL;
				Reset();
				return(eBSFerrMem);
				}
	#endif
			m_AllocdSNPSitesMem = memreq;
			m_NumSNPSites = 0;
			m_NumSortedSites = 0;
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
					gDiagnostics.DiagOut(eDLFatal, gszProcName, "SNPs: Memory re-allocation to %I64d bytes - %s", memreq, strerror(errno));
					return(eBSFerrMem);
					}
				m_pSNPSites = pSNPSite;
				m_NumAllocdSNPSites += cAllocNumSNPs;
				m_AllocdSNPSitesMem = memreq;
				}
			}

		pSNPSite = &m_pSNPSites[m_NumSNPSites++];
		memset(pSNPSite,0,sizeof(tsSNPSite));
		pSNPSite->AlleleBaseCnts[AlleleBase] = 1;
		pSNPSite->ChromID = ChromID;
		pSNPSite->Loci = Loci;
		pSNPSite->RefBase = RefBase;
		}
	delete m_pCSV;
	m_pCSV = NULL;
	if(m_NumSNPSites > 1 && m_NumSNPSites != m_NumSortedSites)
		qsort(m_pSNPSites, m_NumSNPSites, sizeof(tsSNPSite), SortSNPChromLoci);
	m_NumSortedSites = m_NumSNPSites;
	}

uint32_t NumAcceptedSites;
NumAcceptedSites = 0;
pSNPSite = m_pSNPSites;
for(uint32_t SNPIdx = 0; SNPIdx < m_NumSNPSites; SNPIdx++,pSNPSite++)
	{
	pSNPSite->RepBase = eBaseN;
	for(int AlleleIdx = 0; AlleleIdx < 4; AlleleIdx++)
		{
		if(pSNPSite->AlleleBaseCnts[AlleleIdx] > (uint32_t)(m_NumSNPFiles/2))
			{
			NumAcceptedSites+=1;
			pSNPSite->RepBase = AlleleIdx;
			break;
			}
		}
	}

gDiagnostics.DiagOut(eDLInfo, gszProcName, "SNPs: Processed %u files containing %u unique major allele SNP loci of which %u were accepted as being consensus sites", m_NumSNPFiles, m_NumSNPSites, NumAcceptedSites);


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

return(eBSFSuccess);
}


int
CRepAssemb::MergeASegs(char* pszChrom,		// segment is on this chromosome
			uint32_t Loci,		// starting at this loci inclusive
			uint32_t Len,		// segment is this length
			uint32_t Score)		// 0 if segment is an unaligned segment
{
int Rslt;
int ChromID;

tsASegment *pASeg;
static tsASegment *pMRAASeg = NULL;

if(m_pASegments == NULL)
	{
	size_t memreq = (size_t)(sizeof(tsASegment) * cAllocNumASegments);

#ifdef _WIN32
	m_pASegments = (tsASegment*)malloc(memreq);	// initial and perhaps the only allocation
	if(m_pASegments == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "MergeASegs: Initial memory allocation of %I64d bytes - %s", (int64_t)memreq, strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	m_pASegments = (tsASegment *)mmap(NULL, memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	if(m_pASegments == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "MergeASegs: Memory allocation of %I64d bytes through mmap()  failed - %s", (int64_t)memreq, strerror(errno));
		m_pASegments = NULL;
		Reset();
		return(eBSFerrMem);
		}
#endif
	m_AllocdASegmentsMem = memreq;
	m_NumASegments = 0;
	m_NumAllocdASegments = cAllocNumASegments;
	pMRAASeg = NULL;
	}
else
	{
	if((m_NumASegments + 100) > m_NumAllocdASegments)
		{
		size_t memreq = m_AllocdASegmentsMem + (cAllocNumASegments * sizeof(tsASegment));
				
#ifdef _WIN32
		pASeg = (tsASegment*)realloc(m_pASegments, memreq);
#else
		pASeg = (tsASegment*)mremap(m_pASegments, m_AllocdASegmentsMem, memreq, MREMAP_MAYMOVE);
		if(pASeg == MAP_FAILED)
			pASeg = NULL;
#endif
		if(pASeg == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "MergeASegs: Memory re-allocation to %I64d bytes - %s", memreq, strerror(errno));
			return(eBSFerrMem);
			}
		m_pASegments = pASeg;
		m_NumAllocdASegments += cAllocNumASegments;
		m_AllocdASegmentsMem = memreq;
		pMRAASeg = NULL;
		}
	}

if((Rslt = ChromID = AddChrom(pszChrom)) < 1)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Failed generating a unique chromosome name identifier '%s'", pszChrom);
	Reset();
	return(Rslt);
	}

if(Score == 0)			// only accepting segments which were aligned
	return(eBSFSuccess);

// It is highly probable that the previously processed alignment segment is adjacent to this new segment because
// the BED file processing will have resulted in ordering of BED features by chrom.loci
// So optimisation is to check against last accessed segment
if(pMRAASeg != NULL && pMRAASeg->ChromID == ChromID)
	{
	if(((Loci + Len) >= pMRAASeg->Loci) && (Loci <= (pMRAASeg->Loci + pMRAASeg->Len)))
		{
		if(pMRAASeg->Loci <= Loci)
			pMRAASeg->Len = max(pMRAASeg->Len,Loci - pMRAASeg->Loci + Len);
		else
			{
			pMRAASeg->Loci = Loci;
			pMRAASeg->Len = pMRAASeg->Loci - max(pMRAASeg->Len,Loci+Len);
			}
		return(eBSFSuccess);
		}
	}


// iterate over all previously merged segments; if an overlapping or adjoining then merge with that previous segment
// objective is to reduce the number of segments needing to be stored
if(m_NumASegments)
	{
	pASeg = m_pASegments;
	for(uint32_t Idx = 0; Idx < m_NumASegments; Idx++,pASeg++) {
		if(pASeg->ChromID != ChromID)
			continue;
		if((pASeg->Loci + pASeg->Len) < Loci)
			continue;
		if(Loci > pASeg->Loci + pASeg->Len)
			continue;
		// merge
		if(pASeg->Loci <= Loci)
			pASeg->Len = max(pASeg->Len,Loci - pASeg->Loci + Len);
		else
			{
			pASeg->Loci = Loci;
			pASeg->Len = pASeg->Loci - max(pASeg->Len,Loci+Len);
			}
		pMRAASeg = pASeg;
		return(eBSFSuccess);
		}
	}

// unable to merge, use a new aligned segment
pASeg = &m_pASegments[m_NumASegments++];
pASeg->ChromID = ChromID;
pASeg->Loci = Loci;
pASeg->Len = Len;
pMRAASeg = pASeg;
return(eBSFSuccess);
}

int						// eBSFSuccess or error
CRepAssemb::LoadKalignASegs(char *pszASegsFiles)	// load and parse kalign generated aligned segments from this BED file- can be wild-carded
{
int Rslt;
char *pszASegFile;
CBEDfile *pBedFile;
tsASegment *pASegment;

int StartLoci;
int EndLoci;
int Score;

char szChrom[cMaxDatasetSpeciesChrom+1];
char szFeatName[cMaxDatasetSpeciesChrom+1];

if((pBedFile = new CBEDfile) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to instantiate instance of CBEDfile");
	Reset();
	return(eBSFerrInternal);
	}
CSimpleGlob glob(SG_GLOB_FULLSORT);
glob.Init();
if(glob.Add(pszASegsFiles) < SG_SUCCESS)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to glob '%s",pszASegsFiles);
	return(eBSFerrOpnFile);	// treat as though unable to open file
	}

m_MaxASegmentLen = 0;
m_NumSortedSites = 0;
if((m_NumSNPFiles = glob.FileCount()) <= 0)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to locate file(s): '%s'", pszASegsFiles);
	delete pBedFile;
	Reset();
	return(eBSFerrOpnFile);
	}

for(int FileIdx = 0;FileIdx < m_NumSNPFiles; FileIdx++)
	{
	pszASegFile = glob.File(FileIdx);
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Loading alignment regions from file: '%s'", pszASegFile);
	if((Rslt = pBedFile->Open(pszASegFile)) != eBSFSuccess) // loads and processes complete BED file
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to open file: '%s'", pszASegFile);
		delete pBedFile;
		Reset();
		return(Rslt);
		}

	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Completed loading %d alignment segments from file: '%s'", pBedFile->GetNumFeatures(),pszASegFile);
	int CurFeatureID = 0;
	while((CurFeatureID = pBedFile->GetNextFeatureID(CurFeatureID)) > 0)
		{
		pBedFile->GetFeature(CurFeatureID,	// feature instance identifier
					szFeatName,				// where to return feature name
					szChrom,				// where to return chromosome name
					&StartLoci,				// where to return feature start on chromosome (0..n) 
					&EndLoci,				// where to return feature end on chromosome
					&Score);				// where to return score
		if(Score > 0)
			MergeASegs(szChrom,StartLoci,1+EndLoci-StartLoci,Score);
		}
	if(m_NumASegments)
		qsort(m_pASegments,m_NumASegments,sizeof(tsASegment),SortASegmentChromIDLoci);
	pASegment = m_pASegments;
	for(uint32_t Idx = 0; Idx < m_NumASegments; Idx++, pASegment++)
		if(pASegment->Len > m_MaxASegmentLen)
			m_MaxASegmentLen = pASegment->Len;
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "LoadKalignASegs: Processed '%s' resulting in %u merged aligned segments", pszASegFile,m_NumASegments);
	}

gDiagnostics.DiagOut(eDLInfo, gszProcName, "LoadKalignASegs: Processed %u file(s) resulting in %u merged segments", m_NumSNPFiles,m_NumASegments);

return(eBSFSuccess);
}

int
CRepAssemb::LoadAssembly(char* pszAssembly)	// load assembly sequences in this fasta file into memory (m_pSeqBuffer will be allocated to hold raw fasta sequences)
{
int NumRead;
size_t memreq;
uint8_t *pInChr;

#ifdef _WIN32
m_hInFile = open(pszAssembly, O_READSEQ );		// file access is normally sequential..
#else
m_hInFile = open64(pszAssembly, O_READSEQ );		// file access is normally sequential..
#endif
if(m_hInFile == -1)							// check if file open succeeded
	{
	gDiagnostics.DiagOut (eDLFatal, gszProcName, "Unable to open input file '%s' : %s",pszAssembly,strerror(errno));
	return(eBSFerrOpnFile);
	}

// an initial allocation for buffering input assembly, will be realloc'd if required so as to contain complete assembly
memreq = (size_t)cAllocAssembFastaMem;

#ifdef _WIN32
m_pSeqBuffer = (uint8_t *)malloc(memreq);	// initial and perhaps the only allocation
if(m_pSeqBuffer == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Initial memory allocation of %I64d bytes - %s", (int64_t)memreq, strerror(errno));
	Reset();
	return(eBSFerrMem);
	}
#else
// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
m_pSeqBuffer = (uint8_t *)mmap(NULL, memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
if(m_pSeqBuffer == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Assembly: Memory allocation of %I64d bytes through mmap()  failed - %s", (int64_t)memreq, strerror(errno));
	m_pSeqBuffer = NULL;
	Reset();
	return(eBSFerrMem);
	}
#endif

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading fasta assembly from : '%s'",pszAssembly);
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
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Memory re-allocation to %I64d bytes - %s", memreq, strerror(errno));
			return(eBSFerrMem);
			}
		m_pSeqBuffer = pInChr;
		m_AllocSeqBuffMem = memreq;
		}
	m_SeqBuffIdx += NumRead;
	}

if(NumRead < 0)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Loading assembly '%s' failed - %s", pszAssembly,strerror(errno));
	Reset();
	return(eBSFerrFileAccess);
	}

return(eBSFSuccess);
}


int
CRepAssemb::ReportSNPdist(char *pszRepAssembFile)	// base file name used for reporting major allele SNPs used when purposing assembly file
{
tsSNPSite *pSNPSite;
int RepSNPs;
char szOutRepSNPsFile[_MAX_PATH];
char szBuff[4096];
int BuffIdx;
strcpy(szOutRepSNPsFile,pszRepAssembFile);
strcat(szOutRepSNPsFile,".snps.csv");

#ifdef _WIN32
m_hOutRepSNPsFile = open(szOutRepSNPsFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE));
#else
if((m_hOutRepSNPsFile = open64(szOutRepSNPsFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE)) != -1)
	if(ftruncate(m_hOutRepSNPsFile,0)!=0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create/truncate %s - %s",szOutRepSNPsFile,strerror(errno));
		Reset();
		return(eBSFerrCreateFile);
		}
#endif
if(m_hOutRepSNPsFile < 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create/truncate %s - %s",szOutRepSNPsFile,strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}

// write out header line
BuffIdx = sprintf(szBuff,"\"ID\",\"Chrom\",\"Loci\",\"RefBase\",\"RepBase\",\"Instances\"\n");

pSNPSite = m_pSNPSites;
RepSNPs = 0;
for(uint32_t SNPIdx = 0; SNPIdx < m_NumSNPSites; SNPIdx++,pSNPSite++)
	{
	if(pSNPSite->RepBase == eBaseN)
		continue;
	RepSNPs+=1;
	BuffIdx += sprintf(&szBuff[BuffIdx],"%u,\"%s\",%u,%c,%c,%u\n",RepSNPs,LocateChrom(pSNPSite->ChromID),pSNPSite->Loci,CFasta::Base2Chr(pSNPSite->RefBase),CFasta::Base2Chr(pSNPSite->RepBase),pSNPSite->AlleleBaseCnts[pSNPSite->RepBase]);
	if((BuffIdx + 500) > (int)sizeof(szBuff))
		{
		CUtility::RetryWrites(m_hOutRepSNPsFile,szBuff,BuffIdx);
		BuffIdx = 0;
		}
	}
if(BuffIdx)
	CUtility::RetryWrites(m_hOutRepSNPsFile,szBuff,BuffIdx);
// commit and close output file
#ifdef _WIN32
_commit(m_hOutRepSNPsFile);
#else
fsync(m_hOutRepSNPsFile);
#endif
close(m_hOutRepSNPsFile);
m_hOutRepSNPsFile = -1;
return(eBSFSuccess);
}

int
CRepAssemb::ProccessAssembly(etRPPMode PMode,	// processing mode
				uint32_t MinCoverage,			// must be at least this number of covering bases at SNP site to be accepted
				double MinAlleleProp,			// major allele must be at least this proportion of all alleles (incl ref base) at SNP site to be accepted
				double PValueThres,				// SNP PValue must be <= this threshold
				char *pszProcFile,				// input kalign generated SNPs or alignment segments file name, filename can be wildcarded
				char *pszAssembFile,			// input assembly to repurpose
				char *pszRepAssembFile)			// write out repurposed assembly to this file
{
int Rslt;
bool bNL;
bool bDescr;
uint8_t InChr;
uint8_t *pInChr;
size_t FastaIdx;
uint8_t FastaChr;
uint8_t *pszChromName;
uint8_t szChromName[cMaxDatasetSpeciesChrom+1];
uint8_t *pChromNameChr;
int BasesReplaced;
Reset();

if((Rslt = LoadAssembly(pszAssembFile)) < 0)
	{
	Reset();
	return(Rslt);
	}

if(PMode == eRPMdefault)
	{
	m_MinCoverage = MinCoverage;
	m_MinAlleleProp = MinAlleleProp;
	m_PValueThres = PValueThres;
	if((Rslt = LoadKalignSNPs(pszProcFile)) < 0)
		{
		Reset();
		return(Rslt);
		}
	}
else
	{
	if((Rslt = LoadKalignASegs(pszProcFile)) < 0)
		{
		Reset();
		return(Rslt);
		}
	}

#ifdef _WIN32
m_hOutFile = open(pszRepAssembFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE));
#else
if((m_hOutFile = open64(pszRepAssembFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE)) != -1)
	if(ftruncate(m_hOutFile,0)!=0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create/truncate %s - %s",pszRepAssembFile,strerror(errno));
		Reset();
		return(eBSFerrCreateFile);
		}
#endif
if(m_hOutFile < 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create/truncate %s - %s",pszRepAssembFile,strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}

if(PMode == eRPMdefault && ((Rslt = ReportSNPdist(pszRepAssembFile))<eBSFSuccess))
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to report major allele SNPs used");
	Reset();
	return(Rslt);
	}

BasesReplaced = 0;
gDiagnostics.DiagOut(eDLFatal, gszProcName, "Updating assembly sequence with replacement bases ...");

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
			gDiagnostics.DiagOut (eDLFatal, gszProcName, "Expected Fasta '%s' to only contain ascii chars",pszAssembFile);
			Reset();
			return(eBSFerrFastqChr);
		}

	if(InChr == '\n')
		{
		if(bDescr)
			{
			if(PMode == eRPMdefault)
				BasesReplaced += ReplaceSNPBases((char *)szChromName,(uint32_t)min(0x0ffffffff,m_SeqBuffIdx-FastaIdx),(char *)pInChr+1);
			else
				BasesReplaced += ReplaceASegmentBases((char *)szChromName,(uint32_t)min(0x0ffffffff,m_SeqBuffIdx-FastaIdx),(char *)pInChr+1);
			}
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

gDiagnostics.DiagOut(eDLInfo, gszProcName, "%d base replacements completed, writing updated assembly to '%s'",BasesReplaced,pszRepAssembFile);

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


int		// returned chrom identifier, < 1 if unable to locate this chromosome name
CRepAssemb::LocateChrom(char* pszChrom) // return unique identifier associated with this chromosome name
{
int ChromNameIdx;
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
return(eBSFerrChrom);
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
CRepAssemb::ReplaceSNPBases(char *pszChrom, // chromosome name
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
		switch(pSNPSite->RepBase) {
			case eBaseA:
				NumSNPBases++;
				*pszSeq = 'a';
				break;
			case eBaseC:
				NumSNPBases++;
				*pszSeq = 'c';
				break;
			case eBaseG:
				NumSNPBases++;
				*pszSeq = 'g';
				break;
			case eBaseT:
				NumSNPBases++;
				*pszSeq = 't';
				break;
			default:
				break;
			}
		pSNPSite++;
		}
	pszSeq+=1;
	}

return(NumSNPBases);
}

int									// number of bases replaced
CRepAssemb::ReplaceASegmentBases(char *pszChrom, // chromosome name
			 uint32_t SeqLen,			// max of this number of bases in fasta sequence length - includes any whitespace
			 char *pszSeq)			// sequence
{
uint32_t RefLoci;
int ChromID;
tsASegment *pASegment;
uint32_t NumRepBases;

if(SeqLen == 0 || pszChrom == NULL || pszChrom[0] == '\0' || pszSeq == NULL || pszSeq[0] == '\0')
	return(0);

// locate ChromSites for name
if((ChromID = LocateChrom(pszChrom))< 1)
	return(0);

NumRepBases = 0;
RefLoci = 0;

while(SeqLen--)
	{
	switch(*pszSeq) {
		case 'a': case 'A':
		case 'c': case 'C':
		case 'g': case 'G':
		case 't': case 'T':
			break;
		case '>':
		case '\0':
			return(NumRepBases);
		case 'n': case 'N':				// existing indeterminates occupy a loci
			RefLoci++;;
			pszSeq+=1;
			continue;
		default:						// skipping over any whitespace or non-base character
			pszSeq+=1;					// but loci is unchanged
			continue;
		}

	if((pASegment = LocateASegment(ChromID, RefLoci)) == NULL)
		{
		*pszSeq = 'N';
		NumRepBases++;
		}

	RefLoci++;
	pszSeq+=1;
	}

return(NumRepBases);
}

tsSNPSite *
CRepAssemb::LocateSite(uint32_t ChromID, uint32_t Loci)	// locate existing site using a binary search over sorted sites
{
tsSNPSite* pEl1;
int32_t IdxLo;
int32_t MidIdx;
int32_t IdxHi;

if(m_NumSortedSites == 0)
	return(NULL);

if(m_NumSortedSites < 10)				// if just a few then do a simple linear search
	{
	pEl1 = m_pSNPSites;
	for(IdxLo = 0; IdxLo < (int32_t)m_NumSortedSites; IdxLo++,pEl1++)
		if(pEl1->ChromID == ChromID && pEl1->Loci == Loci)
			return(pEl1);
	return(NULL);
	}

IdxLo = 0;
IdxHi = m_NumSortedSites-1;
do {
	MidIdx = (IdxHi + IdxLo)/2;
	pEl1 = &m_pSNPSites[MidIdx];
	if(pEl1->ChromID == ChromID)
		{
		if(pEl1->Loci == Loci)
			return(pEl1);
		if(pEl1->Loci > Loci)
			IdxHi = MidIdx - 1;
		else
			IdxLo = MidIdx + 1;
		}
	else
		{
		if(pEl1->ChromID > ChromID)
			IdxHi = MidIdx - 1;
		else
			IdxLo = MidIdx + 1;
		}
	}
while(IdxHi >= IdxLo);
return(NULL);
}

tsASegment *
CRepAssemb::LocateASegment(uint32_t ChromID, uint32_t Loci)	// locate a tsASegment containing the requested ChromID.Loci
{
static tsASegment *pLastLocatedASegment = NULL;
tsASegment* pEl1;
int32_t MidIdx;
int32_t HiIdx = m_NumASegments - 1;
int32_t LoIdx = 0;
if(m_NumASegments == 0)
	return(NULL);

if(pLastLocatedASegment != NULL && 
		ChromID == pLastLocatedASegment->ChromID && pLastLocatedASegment->Loci <= Loci && (pLastLocatedASegment->Loci+pLastLocatedASegment->Len >  Loci))
	return(pLastLocatedASegment);
pLastLocatedASegment = NULL;

if(m_NumASegments < 10)				// if just a few then do a simple linear search
	{
	pEl1 = m_pASegments;
	for(LoIdx = 0; LoIdx < (int32_t)m_NumASegments; LoIdx++,pEl1++)
		if(pEl1->ChromID == ChromID && pEl1->Loci <= Loci && (pEl1->Loci+pEl1->Len >  Loci))
			{
			pLastLocatedASegment = pEl1;
			return(pEl1);
			}
	return(NULL);
	}

tsASegment *pCurASegment;

int CurStart;
int CurEnd;

while(HiIdx >= LoIdx) {
	MidIdx = (HiIdx + LoIdx)/2;
	pCurASegment = &m_pASegments[MidIdx];
	if(pCurASegment->ChromID == ChromID)
		{
		// read aligned to requested chrom
		CurStart = pCurASegment->Loci;
		CurEnd = pCurASegment->Loci+pCurASegment->Len-1;
		if(Loci >= (uint32_t)CurStart && Loci <= (uint32_t)CurEnd)	// have an overlap?
			{
			pLastLocatedASegment = pCurASegment;
			return(pCurASegment);
			}

		if(CurStart > (int32_t)Loci)		// CurStart after loci then reads sorted higher can't overlap - places an upper limit on search bounds
			HiIdx = MidIdx - 1;
		else                    // CurStart is <= loci but CurEnd is also < target loci
			{
			if((Loci - CurStart) > (int)m_MaxASegmentLen)
				LoIdx = MidIdx + 1;
			else
				LoIdx += 1;
			}
		continue;
		}

	// need to locate chromosome before can start looking for read overlaps onto loci
	if(pCurASegment->ChromID < ChromID)
		LoIdx = MidIdx+1;
	else
		HiIdx = MidIdx - 1;
	};

return(NULL);
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

// SortASegmentChromIDLoci
// Sort m_pASegments by ascending ChromID.Loci
int
CRepAssemb::SortASegmentChromIDLoci(const void* arg1, const void* arg2)
{
tsASegment* pEl1 = (tsASegment*)arg1;
tsASegment* pEl2 = (tsASegment*)arg2;

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


