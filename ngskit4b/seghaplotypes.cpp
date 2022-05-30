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
#include "seghaplotypes.h"


int SegHaplotypes(eModeSH PMode,	// processing mode
			bool bNoSplit,			// true if not to split output files by haplotype tag
			int MinBinScore,		// founder bin must be at least this absolute minimum score before being counted as founder segment presence (default 10)
			double MinBinProp,		// founder bin must be at least this minimum proportion before counted as founder segment presence (default 0.3)
			int SnpMarkerMult,		// boost alignments overlaying SNP markers confidence by this confidence multiplier (default 25)
			char *pszTrackName,		// track name
			char *pszTrackDescr,	// track descriptor
			bool bDontScore,		// don't score haplotype bin segments
			int BinSizeKbp,			// segmentation sized bins
			char *pszSNPMarkers,	// SNP marker loci association file
			char* pszInFile,		// input SAM file
			char* pszOutFile);		// output to this file

#ifdef _WIN32
int seghaplotypes(int argc, char *argv[])
{
	// determine my process name
	_splitpath (argv[0], NULL, NULL, gszProcName, NULL);
#else
int
seghaplotypes(int argc, char **argv)
{
	// determine my process name
	CUtility::splitpath ((char *)argv[0], NULL, gszProcName);
#endif
	int iFileLogLevel;			// level of file diagnostics
	int iScreenLogLevel;		// level of file diagnostics
	char szLogFile[_MAX_PATH];	// write diagnostics to this file
	int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure
	int NumberOfProcessors;		// number of installed CPUs

	 eModeSH PMode;					// processing mode
	 bool bNoSplit;					// true if not splitting haplotype founders
	 bool bDontScore;				// don't score haplotype bin segments

	 int MinBinScore;				// founder bin must be at least this absolute minimum score before being counted as founder segment presence (default 10)
	 double MinBinProp;				// founder bin must be at least this minimum proportion before counted as founder segment presence (default 0.3)
	 int SnpMarkerMult;				// boost alignments overlaying SNP markers confidence by this confidence multiplier (default 25)

	 char szTrackName[_MAX_PATH];	// track name
	 char szTrackDescr[_MAX_PATH];	// track description
	 int BinSizeKbp;				// number of alignments over these sized bins
	 char szInFile[_MAX_PATH];		// input SAM file
	 char szSNPMarkers[_MAX_PATH];	// input SNP marker loci association file
	 char szOutFile[_MAX_PATH];		 // output file

	struct arg_lit *help = arg_lit0 ("h", "help", "print this help and exit");
	struct arg_lit *version = arg_lit0 ("v", "version,ver", "print version information and exit");
	struct arg_int *FileLogLevel = arg_int0 ("f", "FileLogLevel", "<int>", "Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
	struct arg_file *LogFile = arg_file0 ("F", "log", "<file>", "diagnostics log file");

	struct arg_int *pmode = arg_int0 ("m", "mode", "<int>", "processing mode: 0 default is to generate segmentations using bin counts of unique loci only, 1 use all alignments");
	struct arg_lit *nosplit = arg_lit0 ("s", "split", "don't split output files by haplotype tag");
	struct arg_lit *dontscore = arg_lit0 ("n", "noscore", "don't score haplotype segment bins");
	struct arg_int *binsizekbp = arg_int0 ("b", "binsizekbp", "<int>", "Maximum Kbp bin size (default 10, range 1..1000");

	struct arg_int *minbinscore = arg_int0 ("m", "minbinscore", "<int>", "founder bin must be at least this absolute minimum score before being counted as founder segment presence (default 10)");
	struct arg_dbl *minbinprop = arg_dbl0 ("M", "minbinprop", "<int>", "founder bin must be at least this minimum proportion of all founders before accepted as founder segment presence (default 0.3)");
	struct arg_int *snpmarkermult = arg_int0 ("c", "snpmarkermult", "<int>", " boost alignments overlaying SNP markers confidence by this confidence multiplier (default 25)");

	struct arg_str *trackname = arg_str0("t","trackname","<str>","BED Track name");
	struct arg_str *trackdescr = arg_str0("d","trackdescr","<str>","BED Track description");
	struct arg_file *infile = arg_file1("i", "in", "<file>", "SAM/BAM input file(s) (can contain wildcards for multiple SAM file processing)");
	struct arg_file *insnpmarkers = arg_file0("I", "snpmarkers", "<file>", "SNP marker loci association file");
	struct arg_file *outfile = arg_file1 ("o", "out", "<file>", "BED output file");
	struct arg_end *end = arg_end (200);

	void *argtable[] = { help,version,FileLogLevel,LogFile,
						pmode,nosplit,minbinscore,minbinprop,snpmarkermult,trackname,trackdescr,dontscore,binsizekbp,insnpmarkers,infile,outfile,end };

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

		PMode = pmode->count ? (eModeSH)pmode->ival[0] : eMSHDefault;
		if (PMode < eMSHDefault || PMode >= eMSHPlaceHolder)
			{
			gDiagnostics.DiagOut (eDLFatal, gszProcName, "Error: Processing mode '-m%d' specified outside of range %d..%d\n", PMode, eMSHDefault, eMSHPlaceHolder-1);
			exit (1);
			}

		bDontScore = dontscore->count ? true : false;		// don't score haplotype bin segments
		if(trackname->count)
			{
			strncpy (szTrackName, trackname->sval[0],80);
			szTrackName[80] = '\0';
			CUtility::TrimQuotedWhitespcExtd (szTrackName);
			CUtility::CleanText(szTrackName);
			}
		else
			szTrackName[0] = '\0';
		if (szTrackName[0] == '\0')
			strcpy(szTrackName,"SH");

		MinBinScore = minbinscore->count ? minbinscore->ival[0] : 10;
		if(MinBinScore < 1)			// force to be within a reasonable range
			MinBinScore = 1;
		else
			if(MinBinScore > 10000)
				MinBinScore = 10000;

		MinBinProp = minbinprop->count ? minbinprop->dval[0] : 0.3;
		if(MinBinProp < 0.05)			// force to be within a reasonable range
			MinBinProp = 0.05;
		else
			if(MinBinProp > 0.5)
				MinBinProp = 0.5;

		SnpMarkerMult = snpmarkermult->count ? snpmarkermult->ival[0] : 25;
		if(SnpMarkerMult < 1)			// force to be within a reasonable range
			SnpMarkerMult = 1;
		else
			if(SnpMarkerMult > 100)
				SnpMarkerMult = 100;

		bNoSplit = nosplit->count ? true : false;

		if(trackdescr->count)
			{
			strncpy(szTrackDescr, trackdescr->sval[0],80);
			szTrackDescr[80] = '\0';
			CUtility::TrimQuotedWhitespcExtd (szTrackDescr);
			CUtility::CleanText(szTrackDescr);
			}
		else
			szTrackDescr[0] = '\0';
		if (szTrackDescr[0] == '\0')
			strcpy(szTrackDescr,szTrackName);
		szTrackName[40] = '\0';

		BinSizeKbp = binsizekbp->count ? binsizekbp->ival[0] : cDfltSHBinSize;
		if(BinSizeKbp < cMinSHBinSize || BinSizeKbp > cMaxSHBinSize)
			{
			gDiagnostics.DiagOut (eDLFatal, gszProcName, "Error: Windowed maximum bin size specified as '-b%d'Kbp , outside of range %d..%d\n", BinSizeKbp, cMinSHBinSize,cMaxSHBinSize);
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
		if(insnpmarkers->count)
			{
			strcpy (szSNPMarkers, insnpmarkers->filename[0]);
			CUtility::TrimQuotedWhitespcExtd (szSNPMarkers);
			}
		else
			szSNPMarkers[0] = '\0';

		strcpy (szInFile, infile->filename[0]);
		CUtility::TrimQuotedWhitespcExtd (szInFile);
		if (szInFile[0] == '\0')
			{
			gDiagnostics.DiagOut (eDLFatal, gszProcName, "No input file specified");
			exit (1);
			}

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
			case eMSHDefault:
				pszDescr = "generate segmentations using bin counts of unique loci only";
				break;

			case eMSHSegAll:
				pszDescr = "generate segmentations using bin counts of all alignment loci incl non-unique";
				break;
		}


		gDiagnostics.DiagOutMsgOnly (eDLInfo, "Split output files by haplotype segment tag", bNoSplit ? "No" : "Yes");
		gDiagnostics.DiagOutMsgOnly (eDLInfo, "founder bin must be at least this absolute minimum score before being counted as founder segment presence : %d",MinBinScore);
		gDiagnostics.DiagOutMsgOnly (eDLInfo, "founder bin must be at least this minimum proportion before counted as founder segment presence : %f",MinBinProp);
		gDiagnostics.DiagOutMsgOnly (eDLInfo, "boost alignments overlaying SNP markers confidence by this confidence multiplier : %d",SnpMarkerMult);
		
		gDiagnostics.DiagOutMsgOnly (eDLInfo, "Segment haplotypes : '%s'", pszDescr);
		gDiagnostics.DiagOutMsgOnly (eDLInfo, "Track name : '%s'", szTrackName);
		gDiagnostics.DiagOutMsgOnly (eDLInfo, "Track name : '%s'", szTrackDescr);
		gDiagnostics.DiagOutMsgOnly (eDLInfo, "Score haplotype segment bins : '%s'", bDontScore ? "No" : "Yes");
		gDiagnostics.DiagOutMsgOnly (eDLInfo, "Counts accumulated into maximal sized bins of : %dKbp'", BinSizeKbp);
		if(szSNPMarkers[0] != '\0')
			gDiagnostics.DiagOutMsgOnly (eDLInfo, "SNP marker loci file : '%s'", szSNPMarkers);
		gDiagnostics.DiagOutMsgOnly (eDLInfo, "Output file : '%s'", szOutFile);
		gDiagnostics.DiagOutMsgOnly (eDLInfo, "Input file(s) : '%s'", szInFile);
		gDiagnostics.DiagOutMsgOnly (eDLInfo, "Output file : '%s'", szOutFile);


#ifdef _WIN32
		SetPriorityClass (GetCurrentProcess (), BELOW_NORMAL_PRIORITY_CLASS);
#endif
		gStopWatch.Start ();
		Rslt = 0;
		Rslt = SegHaplotypes(PMode,			// processing mode
						bNoSplit,			// true if not to split output files by haplotype tag
						MinBinScore,		// founder bin must be at least this absolute minimum score before being counted as founder segment presence (default 10)
						MinBinProp,			// founder bin must be at least this minimum proportion before counted as founder segment presence (default 0.2)
						SnpMarkerMult,		// boost alignments overlaying SNP markers confidence by this confidence multiplier (default 5)
						szTrackName,		// track name
						szTrackDescr,		// track descriptor
						bDontScore,			// don't score haplotype segment bins
						BinSizeKbp,			// bin size in Kbp
						szSNPMarkers,		// SNP marker loci association file
						szInFile,			// input SAM file
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


CSegHaplotypes::CSegHaplotypes()
{
m_pInBuffer = NULL;
m_pOutBuffer = NULL;
m_pSAMfile=NULL;
m_pBAMalignment =NULL;
m_pSAMloci = NULL;
m_pBins = NULL;
m_pAllocSNPSites = NULL;
m_pSNPMarkerCSV = NULL;
m_hOutFile = -1;			// output file handle
Reset();
}

CSegHaplotypes::~CSegHaplotypes()
{
if(m_hOutFile != -1)
	close(m_hOutFile);
if(m_pInBuffer != NULL)
	delete []m_pInBuffer;
if(m_pOutBuffer != NULL)
	delete []m_pOutBuffer;
if(m_pSAMfile!=NULL)
	delete m_pSAMfile;
if(m_pSNPMarkerCSV != NULL)
	delete m_pSNPMarkerCSV;

if(m_pBAMalignment !=NULL)
	delete m_pBAMalignment;


if(m_pBins != NULL)
	delete []m_pBins;

if (m_pSAMloci != NULL)
	{
#ifdef _WIN32
	free(m_pSAMloci);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (m_pSAMloci != MAP_FAILED)
		munmap(m_pSAMloci, m_AllocdSAMlociMem);
#endif
	}

if(m_pAllocSNPSites != NULL)
	{
#ifdef _WIN32
	free(m_pAllocSNPSites);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (m_pAllocSNPSites != MAP_FAILED)
		munmap(m_pAllocSNPSites, m_AllocMemSNPSites);
#endif
	}
}

void
CSegHaplotypes::Reset(void)	// resets class instance state back to that immediately following instantiation
{
if(m_hOutFile != -1)
	{
	close(m_hOutFile);
	m_hOutFile = -1;
	}

if(m_pInBuffer != NULL)
	{
	delete []m_pInBuffer;
	m_pInBuffer = NULL;
	}

if(m_pOutBuffer != NULL)
	{
	delete []m_pOutBuffer;
	m_pOutBuffer = NULL;
	}

if (m_pSAMloci != NULL)
	{
#ifdef _WIN32
	free(m_pSAMloci);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (m_pSAMloci != MAP_FAILED)
		munmap(m_pSAMloci, m_AllocdSAMlociMem);
#endif
	m_pSAMloci = NULL;
	}

if(m_pAllocSNPSites != NULL)
	{
#ifdef _WIN32
	free(m_pAllocSNPSites);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (m_pAllocSNPSites != MAP_FAILED)
		munmap(m_pAllocSNPSites, m_AllocMemSNPSites);
#endif
	m_pAllocSNPSites = NULL;
	}

if(m_pSNPMarkerCSV != NULL)
	{
	delete m_pSNPMarkerCSV;
	m_pSNPMarkerCSV = NULL;
	}
if(m_pBins != NULL)
	{
	delete []m_pBins;
	m_pBins = NULL;
	}

m_BinSizeKbp = 0;

m_InBuffIdx = 0;
m_AllocInBuff = 0;

m_OutBuffIdx = 0;	
m_AllocOutBuff = 0;

m_CurNumSAMloci = 0;
m_AllocdSAMloci = 0;
m_AllocdSAMlociMem = 0;

m_NumAlignedTargSeqs = 0;
m_LASeqNameID = 0;
m_NumSeqNames = 0;
m_NxtszSeqNameIdx = 0;
m_szSeqNames[0] = '\0';
memset(m_TargSeqs,0,sizeof(m_TargSeqs));

m_LAFounderID = 0;
m_NumFounders = 0;
m_NxtszFounderIdx = 0;
m_szFounders[0] = '\0';
m_szFounderIdx[0] = 0;

m_AllocdBins = 0;

m_MinBinScore = 10;
m_MinBinProp = 0.2;
m_SnpMarkerMult = 5;

m_ParentsCommonAncestor[0] = '\0';
m_NumParents = 0;
memset(m_Parents,0,sizeof(m_Parents));
m_ReadsetIdentifier[0] = '\0';

m_UsedSNPSites = 0;
m_AllocSNPSites = 0;
m_AllocMemSNPSites = 0;

m_NumReadsetNames=0;
m_NxtszReadsetIdx=0; 
m_szReadsetNames[0]='\0';
memset(m_szReadsetIdx,0,sizeof(m_szReadsetIdx));
}


int SegHaplotypes(eModeSH PMode,	// processing mode
			bool bNoSplit,			// true if not to split output files by haplotype tag
			int MinBinScore,		// founder bin must be at least this absolute minimum score before being counted as founder segment presence (default 10)
			double MinBinProp,		// founder bin must be at least this minimum proportion before counted as founder segment presence (default 0.2)
			int SnpMarkerMult,		// boost alignments overlaying SNP markers confidence by this confidence multiplier (default 5)
			char *pszTrackName,		// track name
			char *pszTrackDescr,	// track descriptor
			bool bDontScore,		// don't score haplotype bin segments
			int BinSizeKbp,			// segmentation sized bins
			char *pszSNPMarkers,		// SNP marker loci association file
			char* pszInFile,		// input SAM file
			char* pszOutFile)		// output to this file
{
int Rslt;
CSegHaplotypes* pSegHaplotypes;

if((pSegHaplotypes = new CSegHaplotypes) == NULL)
	{
	gDiagnostics.DiagOut (eDLFatal, gszProcName, "Unable to instantiate instance of CSegHaplotypes");
	return(eBSFerrObj);
	}
Rslt = pSegHaplotypes->Process(PMode,bNoSplit,MinBinScore,MinBinProp,SnpMarkerMult,pszTrackName,pszTrackDescr,bDontScore,BinSizeKbp,pszSNPMarkers,pszInFile,pszOutFile);
delete pSegHaplotypes;
return(Rslt);
}

int 
CSegHaplotypes::Process(eModeSH PMode,			// processing mode
			bool bNoSplit,			// true if not to split output files by haplotype tag
			int MinBinScore,		// founder bin must be at least this absolute minimum score before being counted as founder segment presence (default 10)
			double MinBinProp,		// founder bin must be at least this minimum proportion before counted as founder segment presence (default 0.2)
			int SnpMarkerMult,		// boost alignments overlaying SNP markers confidence by this confidence multiplier (default 5)
			char *pszTrackName,		// track name
			char *pszTrackDescr,	// track descriptor
			bool bDontScore,		// don't score haplotype bin segments
			int BinSizeKbp,			// bin size in Kbp
			char *pszSNPMarkers,		// SNP marker loci association file
			char* pszInFile,		// input SAM file
			char* pszOutFile)		// output to this file
{
int Rslt;
Rslt = GenBinnedSegments(PMode,			// processing mode
			bNoSplit,					// true if not to split output files by haplotype tag (default is to split into separate files)
			MinBinScore,				// founder bin must be at least this absolute minimum score before being counted as founder segment presence (default 10)
			MinBinProp,					// founder bin must be at least this minimum proportion before counted as founder segment presence (default 0.2)
			SnpMarkerMult,				// boost alignments overlaying SNP markers confidence by this confidence multiplier (default 5)
			pszTrackName,				// track name
			pszTrackDescr,				// track descriptor
			bDontScore,					// don't score haplotype bin segments
			BinSizeKbp,					// Wiggle score is number of alignments over this sized bin
			pszSNPMarkers,				// SNP marker loci association file
			pszInFile,					// alignments are in this SAM/BAM file 
			pszOutFile);				// write out Wiggle to this file
return(Rslt);
}

int			// returned len of parsed founder name, 0 if unable to parse out a founder (add 2 to len to obtain start of non-prefixed original string 
CSegHaplotypes::ParseFounder(char* pszIn)			// input null terminated string assumed to contain founder name within first cMaxSHLenPrefix chars
{
int TagLen;
char Chr;
TagLen = 0;
while((Chr = *pszIn++) && Chr != cTagSHTerm1)
	{
	if(!isalnum(Chr))					// tags can only contain alpha-numeric
		return(0);
	if(++TagLen > cMaxSHLenPrefix)
		return(0);
	}
if(Chr == cTagSHTerm1 && *pszIn == cTagSHTerm2 && TagLen >= 1)
	return(TagLen);
return(0);
}


int			// returned number of alignments accepted for the SAM file processed 
CSegHaplotypes::ParseSAMAlignments(char *pszSAMFile)	// SAM file to be processed
{
teBSFrsltCodes Rslt;
int LineLen;
char *pszTxt;
int NumAcceptedAlignments;
int NumAlignmentsProc;
uint32_t NumMissingFeatures = 0;
uint32_t NumUnmapped = 0;

uint32_t NumSNPsOverlaid;
uint32_t AlignmentsWithSNPSites;
char szGenome[cMaxDatasetSpeciesChrom];
char szContig[cMaxDatasetSpeciesChrom+cMaxSHLenPrefix+1];
int ContigLen;
int NumMappedChroms;

tsTargSeq *pTargSeq;
char *pszRefSeqName;
int FounderNameLen;
char szFounder[cMaxSHLenPrefix+3];
int FounderID;
int NumBinsRequired;
tsSHSAMloci *pSAMloci;
int TargID;

NumMappedChroms = 0;
NumBinsRequired = 0;
NumSNPsOverlaid = 0;
AlignmentsWithSNPSites = 0;

// open SAM for reading
if(pszSAMFile == NULL || *pszSAMFile == '\0')
	return(eBSFerrParams);
if(m_pSAMfile != NULL)
	{
	delete m_pSAMfile;
	m_pSAMfile = NULL;
	}

if((m_pSAMfile = new CSAMfile) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenBinnedSegments: Unable to instantiate class CSAMfile");
	return(eBSFerrInternal);
	}

if((Rslt = (teBSFrsltCodes)m_pSAMfile->Open(pszSAMFile)) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenBinnedSegments: Unable to open SAM/BAM format file %s",pszSAMFile);
	Reset();
	return((teBSFrsltCodes)Rslt);
	}

NumAcceptedAlignments = 0;
NumAlignmentsProc = 0;
int NumFounder1IDs = 0;
int NumFounder2IDs = 0;
time_t Then = time(NULL);
time_t Now;
while(Rslt >= eBSFSuccess && (LineLen = m_pSAMfile->GetNxtSAMline((char *)m_pInBuffer)) > 0)
	{
	m_pInBuffer[LineLen] = '\0';
	pszTxt = CUtility::TrimWhitespc((char *)m_pInBuffer);
	if (*pszTxt == '\0')			// simply slough lines which are just whitespace
		continue;

	if(*pszTxt == '@')			// reference sequence dictionary entry?
		{
		if(pszTxt[1] != 'S' && pszTxt[2] != 'Q')
			continue;

		if(3 != sscanf(&pszTxt[3]," AS:%s SN:%s LN:%d",szGenome,szContig,&ContigLen))
			continue;
		
			// parse out the founder name tag - if present - must be at most cMaxSHLenPrefix chars long with separator cTagSHTerm1 and cTagSHTerm1 (currently "|#") between it and actual target sequence name!
		FounderNameLen = ParseFounder(szContig);
		if(FounderNameLen > 0)
			{
			strncpy(szFounder,szContig,FounderNameLen);
			szFounder[FounderNameLen] = '\0';
			pszRefSeqName = &szContig[FounderNameLen + 2];
			if((FounderID = AddFounder(szFounder))==0) // treating founder errors as if chrom name errors
				{
				Reset();
				return(eBSFerrChrom);
				}
			}
		else
			{
			pszRefSeqName = szContig;
			FounderID = AddFounder((char *)"NA"); // using this as the default if founder not specified
			}

		if((TargID = AddTargSeqName(pszRefSeqName, ContigLen)) == 0) // treating founder errors as if chrom name errors
			{
			Reset();
			return(eBSFerrChrom);
			}
		NumMappedChroms += 1;
		continue;
		}

// number of founders is now known
// create bed files to hold alignments for each of the founders

//
	// now processing the alignments
	NumAlignmentsProc += 1;
	if (!(NumAlignmentsProc % 100000) || NumAlignmentsProc == 1)
		{
		Now = time(NULL);
		if ((Now - Then) >= 60)
			{
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "Accepted %zd SAM/BAM alignments", NumAcceptedAlignments);
			Then += 60;
			}
		}

	// only interest is in the reference chrom name, startloci, length
	if ((Rslt = (teBSFrsltCodes)m_pSAMfile->ParseSAM2BAMalign(pszTxt, m_pBAMalignment, NULL,true)) < eBSFSuccess)
		{
		if (Rslt == eBSFerrFeature)	// simply ignoring SAM alignments which are incomplete - I have observed some aligners with inconsistent/incomplete alignment detail
			{
			NumMissingFeatures++;
			Rslt = eBSFSuccess;
			continue;
			}
		break;
		}

	// check if read has been mapped, if not then slough ...
	if (m_pBAMalignment->refID == -1 || (m_pBAMalignment->flag_nc >> 16) & 0x04)
		{
		NumUnmapped++;
		Rslt = eBSFSuccess;
		continue;
		}

	// parse out the founder name tag - if present!
	// if not present then treat as if founder "NA"
	pszRefSeqName = m_pBAMalignment->szRefSeqName;
	FounderNameLen = ParseFounder(pszRefSeqName);
	if(FounderNameLen > 0)
		{
		strncpy(szFounder,pszRefSeqName,FounderNameLen);
		szFounder[FounderNameLen] = '\0';
		pszRefSeqName += (size_t)FounderNameLen + 2;
		FounderID = LocateFounder(szFounder); 
		}
	else
		FounderID = LocateFounder((char *)"NA");

	if(FounderID ==0) // treating founder errors as if chrom name errors
			{
			Reset();
			return(eBSFerrChrom);
			}

	if((pTargSeq = LocateTargSeq(pszRefSeqName)) == NULL) // target unknown then simply slough
		{
		NumMissingFeatures++;
		Rslt = eBSFSuccess;
		continue;
		}
	if(!pTargSeq->fAligned)  // first alignment to this targ seq?
		{
		pTargSeq->fAligned = true;												// at least one alignment to this target
		m_NumAlignedTargSeqs++;
		pTargSeq->NumBins = (pTargSeq->TargSeqLen / (m_BinSizeKbp * 1000)) + 1;
		pTargSeq->BinsOfs = m_AllocdBins;
		m_AllocdBins += pTargSeq->NumBins;
		}

	// can now access the alignment loci info
	// hard or soft clipping is currently of no interest as using a large sliding window bin and counting loci in that bin
	if(m_CurNumSAMloci > m_AllocdSAMloci) // needing to realloc?
		{
		size_t memreq = m_AllocdSAMlociMem + ((size_t)cAllocSHNumSAMloci * sizeof(tsSHSAMloci));
		uint8_t* pTmp;
#ifdef _WIN32
		pTmp = (uint8_t*)realloc(m_pSAMloci, memreq);
#else
		pTmp = (uint8_t*)mremap(m_pSAMloci, m_AllocdSAMlociMem, memreq, MREMAP_MAYMOVE);
		if (pTmp == MAP_FAILED)
			pTmp = NULL;
#endif
		if (pTmp == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenBinnedSegments: Memory re-allocation to %zd bytes - %s", (int64_t)(memreq), strerror(errno));
			Reset();
			return(eBSFerrMem);
			}
		m_pSAMloci = (tsSHSAMloci *)pTmp;
		m_AllocdSAMlociMem = memreq;
		m_AllocdSAMloci += (size_t)cAllocSHNumSAMloci;
		}
	pSAMloci = &m_pSAMloci[m_CurNumSAMloci++];
	pSAMloci->TargID = pTargSeq->TargSeqID;
	pSAMloci->FounderID = FounderID;
	pSAMloci->TargLoci = m_pBAMalignment->pos;
	pSAMloci->AlignLen = m_pBAMalignment->end + 1 - m_pBAMalignment->pos;
	pSAMloci->bASense = ((m_pBAMalignment->flag_nc >> 16) & cSAMFlgAS) ? 1 : 0;
	pSAMloci->NumMarkerSNPs = 0;
	pSAMloci->Cnt = 1;

// following code determining if read alignment overlapped a marker SNP  will likely need much optimization!
// currently just a linear scan along the alignment checking if a loci on the alignment same as a marker SNP loci
// not checking if alignment base is the expected, as assuming no substitutions were allowed in the alignment, and
// even if a few substitution were allowed then probability of a substitution at the exact marker loci is relatively low
	if(m_UsedSNPSites)
		{
		tsSHSNPSSite *pSNPSite;
		int SitesInAlignment = 0;
		for(uint32_t ScanIdx = pSAMloci->TargLoci; ScanIdx < pSAMloci->TargLoci+pSAMloci->AlignLen; ScanIdx++)
			if((pSNPSite=LocateSNPSite(pSAMloci->TargID,ScanIdx))!=NULL)
				{
				if(pSAMloci->NumMarkerSNPs < 127)
					pSAMloci->NumMarkerSNPs += 1;
				NumSNPsOverlaid++;
				SitesInAlignment++;
				}
		if(SitesInAlignment)
			AlignmentsWithSNPSites++;
		}
	NumAcceptedAlignments++;
	}

if(m_pSAMfile != NULL)
	{
	delete m_pSAMfile;
	m_pSAMfile = NULL;
	}

if(NumAcceptedAlignments == 0)		// ugh, no alignments!
	{
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "ParseSAMAlignments: No alignments accepted for processing!");
	Reset();
	return(eBSFSuccess);			// not an error!
	}
else
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "ParseSAMAlignments: Accepted %d total alignments overlaying %d SNP markers onto %d target seqs for binning from  %d processed",NumAcceptedAlignments, AlignmentsWithSNPSites, m_NumAlignedTargSeqs, NumAlignmentsProc);

return(NumAcceptedAlignments);
}


// generate BEDs for each founder skim read alignments
// objective is to be able to view just what the distribution of these skim alignments along the chromosomes looks like - are they IID distributed or clustering?
int
CSegHaplotypes::GenerateAlignmentBEDs(char *pszSAMfile)	// generate BEDs for each founder skim read alignments which have been parsed from this SAM file
{
uint32_t FounderID;
uint32_t LociIdx;
tsSHSAMloci *pSAMloci;
char szFounderOutFile[_MAX_PATH];
char szFileName[_MAX_PATH];

char *pszFounder;

if(m_hOutFile != -1)
	{
	close(m_hOutFile);
	m_hOutFile = -1;
	}
m_OutBuffIdx = 0;

// <chrom> <chromStart (0..n)> <chromEnd (chromStart + Len)> <name> <score (0..999)
for(FounderID = 1; FounderID <= m_NumFounders; FounderID++)
	{
	pszFounder = LocateFounder(FounderID);
	sprintf(szFounderOutFile,"%s.%s.bed",pszSAMfile,pszFounder);
#ifdef _WIN32
	m_hOutFile = open(szFounderOutFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE));
#else
	if((m_hOutFile = open64(szFounderOutFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE)) != -1)
		if(ftruncate(m_hOutFile,0)!=0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create/truncate %s - %s",szFounderOutFile,strerror(errno));
			Reset();
			return(eBSFerrCreateFile);
			}
#endif
	if(m_hOutFile < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create/truncate %s - %s",szFounderOutFile,strerror(errno));
		Reset();
		return(eBSFerrCreateFile);
		}
	CUtility::splitpath(pszSAMfile,NULL,szFileName);
	szFileName[strlen(szFileName)-4] = '\0';
	m_OutBuffIdx = sprintf((char *)m_pOutBuffer,"track name=\"FAL %s.%s\" description=\"Founder Alignment Loci %s.%s\"\n",szFileName,pszFounder,szFileName,pszFounder);
	
	pSAMloci = m_pSAMloci;
	for(LociIdx = 0; LociIdx < m_CurNumSAMloci; LociIdx++,pSAMloci++)
		{
		if(pSAMloci->FounderID != FounderID)
			continue;

		m_OutBuffIdx += sprintf((char *)&m_pOutBuffer[m_OutBuffIdx],"%s\t%u\t%u\t%s\n",LocateTargSeqName(pSAMloci->TargID),pSAMloci->TargLoci,pSAMloci->TargLoci+pSAMloci->AlignLen,pszFounder);
		if((m_OutBuffIdx + 5000) > m_AllocOutBuff)
			{
			if (!CUtility::RetryWrites(m_hOutFile,m_pOutBuffer,m_OutBuffIdx))
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenerateAlignmentBEDs: Fatal error in RetryWrites()");
				Reset();
				return(eBSFerrFileAccess);
				}
			m_OutBuffIdx = 0;
			}
		}

	if(m_OutBuffIdx)
		{
		if (!CUtility::RetryWrites(m_hOutFile,m_pOutBuffer,m_OutBuffIdx))
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenerateAlignmentBEDs: Fatal error in RetryWrites()");
			Reset();
			return(eBSFerrFileAccess);
			}
		m_OutBuffIdx = 0;
		}
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

// 
// generate segmented haplotypes file from SAM/BAM alignment input file
int	
CSegHaplotypes::GenBinnedSegments(eModeSH PMode,			// processing mode
			bool bNoSplit,			// true if not to split output files by haplotype tag (default is to split into separate files)
			int MinBinScore,		// founder bin must be at least this absolute minimum score before being counted as founder segment presence (default 10)
			double MinBinProp,		// founder bin must be at least this minimum proportion before counted as founder segment presence (default 0.2)
			int SnpMarkerMult,		// boost alignments overlaying SNP markers confidence by this confidence multiplier (default 5)
			char *pszTrackName,		// track name
			char *pszTrackDescr,	// track descriptor
			bool bDontScore,		// don't score haplotype bin segments
			uint32_t BinSizeKbp,	// Wiggle score is number of alignments over this sized bins
			char *pszSNPMarkers,		// SNP marker loci association file
			char* pszInFile,		// alignments are in this SAM/BAM file 
			char* pszOutFile)		// write out segments to this file
{
teBSFrsltCodes Rslt;
uint32_t TargSeqIdx;
uint32_t BinStartLoci;
tsSHSAMloci *pSAMloci;
tsSHSAMloci *pUniqueSamLoci;
pSAMloci = m_pSAMloci;
uint32_t NumUniqueLoci;
size_t LociIdx;
tsSHBin *pBin;
tsTargSeq *pTargSeq;
size_t NumAcceptedAlignments;
int WindowSize = BinSizeKbp * 1000;
m_BinSizeKbp = BinSizeKbp;
m_bNoSplit = bNoSplit;
m_MinBinScore = MinBinScore;
m_MinBinProp = MinBinProp;
m_SnpMarkerMult = SnpMarkerMult;

if((m_pBAMalignment = new tsBAMalign) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenBinnedSegments: Unable to instantiate tsBAMalign");
	Reset();
	return(eBSFerrInternal);
	}

if((m_pInBuffer = new uint8_t[cMaxReadLen * 3]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenBinnedSegments: Unable to allocate buffers");
	Reset();
	return(eBSFerrMem);
	}
m_AllocInBuff = cMaxReadLen * 3;

if((m_pOutBuffer = new uint8_t[cAllocSHBuffOutSize]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenBinnedSegments: Unable to allocate buffers");
	Reset();
	return(eBSFerrMem);
	}
m_AllocOutBuff = cAllocSHBuffOutSize;



m_AllocdSAMlociMem = (size_t)cAllocSHNumSAMloci * sizeof(tsSHSAMloci);
#ifdef _WIN32
m_pSAMloci = (tsSHSAMloci*)malloc(m_AllocdSAMlociMem);
if (m_pSAMloci == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenBinnedSegments: Memory allocation of %zd bytes failed", (int64_t)m_AllocdSAMlociMem);
	m_AllocdSAMlociMem = 0;
	Reset();
	return(eBSFerrMem);
	}
#else
// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
m_pSAMloci = (tsSHSAMloci*)mmap(NULL, m_AllocdSAMlociMem, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
if (m_pSAMloci == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenBinnedSegments: Memory allocation of %zd bytes through mmap()  failed", (int64_t)m_AllocdSAMlociMem, strerror(errno));
	m_pSAMloci = NULL;
	m_AllocdSAMlociMem = 0;
	Reset();
	return(eBSFerrMem);
	}
#endif
m_AllocdSAMloci = cAllocSHNumSAMloci;

NumAcceptedAlignments = 0;
m_AllocdBins = 0;

// load SNP markers if file has been specified
if(!(pszSNPMarkers == NULL || pszSNPMarkers[0] == '\0'))
	{
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Processing SNPs from snpmarkers file '%s", pszSNPMarkers);
	if((m_pSNPMarkerCSV = new CCSVFile)==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenBinnedSegments: Unable to instantiate class CCSVFile for '%s'",pszSNPMarkers);
		Reset();
		return(eBSFerrInternal);
		}
	if((m_pSNPMarkerCSV->Open(pszSNPMarkers))!=eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenBinnedSegments: Unable to open SNP markers file '%s'",pszSNPMarkers);
		Reset();
		return(eBSFerrOpnFile);
		}
	if((Rslt = (teBSFrsltCodes)ProcessSnpmarkersSNPs(pszSNPMarkers)) < eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenBinnedSegments: Errors processing SNP markers file '%s'",pszSNPMarkers);
		Reset();
		return(Rslt);
		}
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Completed processing SNPs (%u accepted) from snpmarkers file '%s", (uint32_t)(m_UsedSNPSites & 0x0ffffffff), pszSNPMarkers);
	}
	
CSimpleGlob glob(SG_GLOB_FULLSORT);
glob.Init();
if(glob.Add(pszInFile) < SG_SUCCESS)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to glob '%s", pszInFile);
	Reset();
	return(eBSFerrOpnFile);	// treat as though unable to open file
	}

if(glob.FileCount() <= 0)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to locate any SAM file matching '%s", pszInFile);
	Reset();
	return(eBSFerrOpnFile);	// treat as though unable to open file
	}

Rslt = eBSFSuccess;
m_NumAlignedTargSeqs = 0;
m_CurNumSAMloci = 0;
char *pszSAMFile;
for(int FileID = 0; Rslt >= eBSFSuccess && FileID < glob.FileCount(); ++FileID)
	{
	pszSAMFile = glob.File(FileID);
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Processing SAM file '%s'\n", pszSAMFile);
	if((Rslt = (teBSFrsltCodes)ParseSAMAlignments(pszSAMFile)) < eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenBinnedSegments: Fatal error processing '%s'",pszSAMFile);
		return(Rslt);
		}
	NumAcceptedAlignments += (size_t)Rslt;
	}

if(m_CurNumSAMloci == 0)		// ugh, no alignments!
	{
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenBinnedSegments: No alignments accepted for processing!");
	Reset();
	return(eBSFSuccess);			// not an error!
	}
else
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenBinnedSegments: Accepted %zd total alignments onto %d targets for binning",m_CurNumSAMloci, m_NumAlignedTargSeqs);

if((m_pBins = new tsSHBin[m_AllocdBins]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenBinnedSegments: Unable to allocate memory for %u alignment bins",m_AllocdBins);
	Reset();
	return(eBSFerrMem);
	}
memset(m_pBins,0,m_AllocdBins * sizeof(tsSHBin));

pTargSeq = m_TargSeqs;
for(TargSeqIdx = 0; TargSeqIdx < m_NumSeqNames; TargSeqIdx++,pTargSeq++)
	{
	if(!pTargSeq->fAligned)
		continue;

	BinStartLoci = 0;
	pBin = &m_pBins[pTargSeq->BinsOfs];

	while(BinStartLoci < pTargSeq->TargSeqLen)
		{
		pBin->TargID = pTargSeq->TargSeqID;
		pBin->TargLoci = BinStartLoci;
		pBin->BinLen = min(m_BinSizeKbp * 1000,pTargSeq->TargSeqLen - BinStartLoci);
		BinStartLoci += pBin->BinLen;
		pBin++;
		}
	}

// sort SAMloci by FounderID.TargID.TargLoci ascending
if(m_CurNumSAMloci > 1)
	qsort(m_pSAMloci, m_CurNumSAMloci, sizeof(tsSHSAMloci), SortSAMFounderTargLoci);

if((Rslt =(teBSFrsltCodes)GenerateAlignmentBEDs(pszSAMFile)) < eBSFSuccess)	
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenBinnedSegments: Fatal error processing '%s' for reporting of alignments",pszSAMFile);
	return(Rslt);
	}

// could be multiple alignments to same loci so reduce such that each loci is unique with a count of alignments to that loci unless unique loci in which case only 1 count attributed
pSAMloci = m_pSAMloci;
pUniqueSamLoci = pSAMloci++;
pUniqueSamLoci->Cnt = 1;
NumUniqueLoci = 1;
int NumFounder1IDs = 0;
int NumFounder2IDs = 0;
for(LociIdx = 1; LociIdx < m_CurNumSAMloci; LociIdx++, pSAMloci++)
	{
	if(pSAMloci->FounderID == 1)
		NumFounder1IDs++;
	else
		NumFounder2IDs++;	
	if(pSAMloci->FounderID == pUniqueSamLoci->FounderID && pSAMloci->TargID == pUniqueSamLoci->TargID && pSAMloci->TargLoci == pUniqueSamLoci->TargLoci)
		{
		if(PMode == eMSHSegAll) // if all to be counted then counting multiple alignments to this loci
			pUniqueSamLoci->Cnt++; 
		}
	else // onto another loci
		{
		pUniqueSamLoci++;
		*pUniqueSamLoci = *pSAMloci;
		pUniqueSamLoci->Cnt = 1;
		NumUniqueLoci++;
		}
	}
m_CurNumSAMloci = NumUniqueLoci;
gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenBinnedSegments: Alignments were to %zd unique loci",m_CurNumSAMloci);

pSAMloci = m_pSAMloci;
NumFounder1IDs = 0;
NumFounder2IDs = 0;
for(LociIdx = 0; LociIdx < m_CurNumSAMloci; LociIdx++,pSAMloci++)
	{
	pTargSeq = &m_TargSeqs[pSAMloci->TargID-1];
	pBin = &m_pBins[pTargSeq->BinsOfs];
	pBin += pSAMloci->TargLoci / (m_BinSizeKbp * 1000);
	pBin->RawCnts[pSAMloci->FounderID-1]++;
	if(pSAMloci->NumMarkerSNPs)
		pBin->RawCnts[pSAMloci->FounderID-1]+=((m_SnpMarkerMult-1) * pSAMloci->NumMarkerSNPs * pSAMloci->Cnt);
	if(pSAMloci->FounderID == 1)
		NumFounder1IDs++;
	else
		NumFounder2IDs++;	
	}

// apply weighted smoothing. 50% of adjacent raw bin counts contribute to each bin smoothed counts
ApplySmoothing();
// now attempt to call the haplotype segments!!!!!!!!!!!!!!!!!!!!!!!!!!
// there can be multiple haplotypes in the same bin
int TotCalledBins=0;
int CalledBins;

TotCalledBins = IdentifySegments(false,bDontScore);	// identify higher confidence bins
while((CalledBins = IdentifySegments(true,bDontScore)) > 0) // these will 'fill in the gaps' ...
	TotCalledBins += CalledBins;

uint32_t FounderIdx;
uint32_t FounderID;
char *pszFounder;
char szFounderOutFile[_MAX_PATH];
char szFounderTrackName[_MAX_PATH];
char szFounderTrackDescr[_MAX_PATH];

for(FounderIdx = 0; FounderIdx < m_NumFounders; FounderIdx++)
	{
	if(bNoSplit)
		{
		FounderID = 0;
		strcpy(szFounderOutFile,pszOutFile);
		strcpy(szFounderTrackName,pszTrackName);
		strcpy(szFounderTrackDescr,pszTrackDescr);
		}
	else
		{
		FounderID = FounderIdx+1;
		pszFounder=LocateFounder(FounderIdx+1);
		sprintf(szFounderOutFile,"%s.%s.bed",pszOutFile,pszFounder);
		sprintf(szFounderTrackName,"%s %s:%s",pszTrackName,pszOutFile,pszFounder);
		sprintf(szFounderTrackDescr,"%s %s:%s",pszTrackDescr,pszOutFile,pszFounder);
		}
#ifdef _WIN32
	m_hOutFile = open(szFounderOutFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE));
#else
	if((m_hOutFile = open64(szFounderOutFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE)) != -1)
		if(ftruncate(m_hOutFile,0)!=0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create/truncate %s - %s",szFounderOutFile,strerror(errno));
			Reset();
			return(eBSFerrCreateFile);
			}
#endif
	if(m_hOutFile < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create/truncate %s - %s",szFounderOutFile,strerror(errno));
		Reset();
		return(eBSFerrCreateFile);
		}
	
	// output as a BED file with scoring
	genBED(FounderID,szFounderTrackName,szFounderTrackDescr);

	if(m_OutBuffIdx)
		{
		if (!CUtility::RetryWrites(m_hOutFile,m_pOutBuffer,m_OutBuffIdx))
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenBinnedSegments: Fatal error in RetryWrites()");
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
	if(bNoSplit)
		break;
	}

Reset();
return(eBSFSuccess);
}


// Smoothing
int
CSegHaplotypes::ApplySmoothing(void)	// Cnts in immediately adjacent bin counts are weighted at 0.5 and these weighted counts added to the current bin cnts as smoothed bin counts
{
tsTargSeq *pTargSeq;
tsSHBin *pBin;
tsSHBin *pWeightingBin;
uint32_t BinIdx;
uint32_t FounderIdx;
uint32_t TargSeqIdx;
uint32_t WeightedCnts;
for(FounderIdx = 0; FounderIdx < m_NumFounders; FounderIdx++)
	{
	pTargSeq = m_TargSeqs;
	for(TargSeqIdx = 0; TargSeqIdx < m_NumSeqNames; TargSeqIdx++, pTargSeq++)
		{
		pBin = &m_pBins[pTargSeq->BinsOfs];
		for(BinIdx = 0; BinIdx < pTargSeq->NumBins; BinIdx++, pBin++)
			{
			WeightedCnts = pBin->RawCnts[FounderIdx];
			if(BinIdx > 0)
				{
				pWeightingBin = pBin - 1;
				WeightedCnts += pWeightingBin->RawCnts[FounderIdx]/2;
				}
			if(BinIdx < pTargSeq->NumBins-1)
				{
				pWeightingBin = pBin + 1;
				WeightedCnts += pWeightingBin->RawCnts[FounderIdx]/2;
				}
			pBin->SmoothedCnts[FounderIdx] = WeightedCnts;
			}
		}
	}
return(0);
}

// now attempt to call the haplotype segments!!!!!!!!!!!!!!!!!!!!!!!!!!
// when bin calling, counts in immediately adjacent bins contribute to the haplotype(s) call for any individual bin  
int
CSegHaplotypes::IdentifySegments(bool bInterpolate,	// if true - if bin counts less than minimum coverage then interpolate from adjacent bins which were called
						bool bDontScore)	// if true then mark bin haplotype with score of cBEDNoScore without actually scoring according to bin counts
{
tsTargSeq *pTargSeq;
tsSHBin *pBin;
tsSHBin *pCpyBin;

tsSHBin *pPrevBin;
tsSHBin *pNextBin;

uint32_t BinIdx;
uint32_t TotBinCnts;
uint32_t FounderIdx;
uint32_t TargSeqIdx;
CStats Stats;
int NumCalledBins;

NumCalledBins = 0;
pTargSeq = m_TargSeqs;
for(TargSeqIdx = 0; TargSeqIdx < m_NumSeqNames; TargSeqIdx++,pTargSeq++)
	{
	if(!pTargSeq->fAligned)
		continue;
	pBin = &m_pBins[pTargSeq->BinsOfs];
	for(BinIdx = 0; BinIdx < pTargSeq->NumBins; BinIdx++, pBin++)
		{
		if(pBin->fCalled)	// if already called then don't call again
			continue;

		// uncalled bin
		TotBinCnts = 0;
		for(FounderIdx = 0; FounderIdx < m_NumFounders; FounderIdx++)
			TotBinCnts += pBin->SmoothedCnts[FounderIdx];

		if(TotBinCnts < (uint32_t)m_MinBinScore) // need at least this many counts to have any degree of confidence for this bin to be a seed bin
			{
			if(bInterpolate)	
				{
				// check neighbors as a last resort, if they were called then assume this bin is an extension of that haplotype
				if(BinIdx > 0)
					pPrevBin = pBin - 1;
				else
					pPrevBin = NULL;
				if((BinIdx + 1) < pTargSeq->NumBins)
					pNextBin = pBin + 1;
				else
					pNextBin = NULL;
				pCpyBin = NULL;

				// choose which bin to copy, nearly pseudo-random!
				if(BinIdx & 0x01)
					{
					if(pPrevBin != NULL && pPrevBin->fCalled)
						pCpyBin = pPrevBin;
					else
						if(pNextBin != NULL && pNextBin->fCalled)
							pCpyBin = pNextBin;
					}
				else
					{
					if(pNextBin != NULL && pNextBin->fCalled)
						pCpyBin = pNextBin;
					else
						if(pPrevBin != NULL && pPrevBin->fCalled)
							pCpyBin = pPrevBin;
					}

				if(pCpyBin != NULL && pCpyBin->fCalled)
					{
					for(FounderIdx = 0; FounderIdx < m_NumFounders; FounderIdx++)
						{
						if(pCpyBin->CalledHaplotype[FounderIdx] > 0)
							pBin->CalledHaplotype[FounderIdx]=bDontScore ? cBEDNoScore : 1;
						else
							pBin->CalledHaplotype[FounderIdx]=0;
						}
					pBin->fCalled = true;
					pBin->fInfer = true;
					NumCalledBins++;
					}
				}
			continue;
			}

		// look at proportions of founder counts in this bin
		// only calling a founder haplotype if that founder is a least MinBinProp as a proportion of total of all founders in that bin
		// otherwise call as having multiple haplotypes in that bin
		for(FounderIdx = 0; FounderIdx < m_NumFounders; FounderIdx++)
			{
			if(pBin->SmoothedCnts[FounderIdx] >= (uint32_t)m_MinBinScore && ((double)pBin->SmoothedCnts[FounderIdx] / TotBinCnts) >= m_MinBinProp)
				{
				if(bDontScore)
					pBin->CalledHaplotype[FounderIdx] = cBEDNoScore;
				else
					{
					int NormalisedCnts;
					NormalisedCnts = max(2,(pBin->SmoothedCnts[FounderIdx] * 100000) / pBin->BinLen);		// counts per 100Kbp
					pBin->CalledHaplotype[FounderIdx] = min(NormalisedCnts, 999);
					}
				}
			}
		pBin->fCalled = true;
		NumCalledBins++;
		}
	}
return(NumCalledBins);
}

int
CSegHaplotypes::genBED(uint32_t SplitByFounderID, // if 0 then combining all founders in same output BED, otherwise BED to only contain haplotype segments specific to the founder specified
						char *pszTrackName,		// track name
						char *pszTrackDescr)	// track descriptor
{
uint32_t FounderIdx;
uint32_t TargSeqIdx;
uint32_t StartLoci;
uint32_t EndLoci;
uint32_t BinIdx;
uint32_t Score;
tsTargSeq *pTargSeq;
tsSHBin *pBin;

if((m_OutBuffIdx + 5000) > m_AllocOutBuff)
	{
	if (!CUtility::RetryWrites(m_hOutFile,m_pOutBuffer,m_OutBuffIdx))
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "genBED: Fatal error in RetryWrites()");
		Reset();
		return(eBSFerrFileAccess);
		}
	m_OutBuffIdx = 0;
	}

m_OutBuffIdx = sprintf((char *)m_pOutBuffer,"track name=\"%s\" description=\"%s\" useScore=1\n",pszTrackName,pszTrackDescr);

// <chrom> <chromStart (0..n)> <chromEnd (chromStart + Len)> <name> <score (0..999)
for(FounderIdx = 0; FounderIdx < m_NumFounders; FounderIdx++)
	{
	if(SplitByFounderID != 0 && FounderIdx != SplitByFounderID-1)
		continue;
	pTargSeq = m_TargSeqs;
	for(TargSeqIdx = 0; TargSeqIdx < m_NumSeqNames; TargSeqIdx++,pTargSeq++)
		{
		if(!pTargSeq->fAligned)
			continue;
		Score = 0;
		pBin = &m_pBins[pTargSeq->BinsOfs]; 
		for(BinIdx = 0; BinIdx < pTargSeq->NumBins; BinIdx++, pBin++)
			{
			if(!pBin->fCalled)	// if not called then skip over
				continue;
			if(pBin->CalledHaplotype[FounderIdx]==0)
				continue;
			Score = pBin->CalledHaplotype[FounderIdx];
			StartLoci = EndLoci = pBin->TargLoci;
			// look ahead and extend segment if adjacent and having same score
			for(; BinIdx < pTargSeq->NumBins; BinIdx++, pBin++)
				{
				if(pBin->CalledHaplotype[FounderIdx] != Score) // must be same score (will be 0 if not extension of same haplotype segment)
					{
					if(pBin->CalledHaplotype[FounderIdx] != 0)	// same haplotype segment but different score?
						{
						BinIdx--;		// ensure segment will be extended next iteration
						pBin--;
						}
					break;
					}
				EndLoci += pBin->BinLen;					 // bin is being extended as is same score or not reporting scores
				}

			m_OutBuffIdx += sprintf((char *)&m_pOutBuffer[m_OutBuffIdx],"%s\t%u\t%u\t%s\t%d\n",LocateTargSeqName(pTargSeq->TargSeqID),StartLoci,EndLoci,LocateFounder(FounderIdx+1),Score);
			if((m_OutBuffIdx + 5000) > m_AllocOutBuff)
				{
				if (!CUtility::RetryWrites(m_hOutFile,m_pOutBuffer,m_OutBuffIdx))
					{
					gDiagnostics.DiagOut(eDLFatal, gszProcName, "genBED: Fatal error in RetryWrites()");
					Reset();
					return(eBSFerrFileAccess);
					}
				m_OutBuffIdx = 0;
				}
			}
		}
	}
if(m_OutBuffIdx)
	{
	if (!CUtility::RetryWrites(m_hOutFile,m_pOutBuffer,m_OutBuffIdx))
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "genBED: Fatal error in RetryWrites()");
		Reset();
		return(eBSFerrFileAccess);
		}
	m_OutBuffIdx = 0;
	}

return(eBSFSuccess);
}

int		// bin score assigned
CSegHaplotypes::GenBinnedCoverage(int TargID,	// coverage is on this targeted chrom/seq
				uint32_t BinStart,		// coverage starts at this loci inclusive
				uint32_t BinEnd,		// coverage ends at this loci inclusive
				uint32_t BinCnts)		// total number of coverage counts
{
uint32_t BinLen;
BinLen = 1 + BinEnd - BinStart; // inclusive!
if((m_OutBuffIdx + 1000) > m_AllocOutBuff)
	{
	if (!CUtility::RetryWrites(m_hOutFile,m_pOutBuffer,m_OutBuffIdx))
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenBinnedCoverage: Fatal error in RetryWrites()");
		Reset();
		return(eBSFerrFileAccess);
		}
	m_OutBuffIdx = 0;
	}
m_OutBuffIdx += sprintf((char *)&m_pOutBuffer[m_OutBuffIdx],"variableStep chrom=%s span=%d\n%d %d\n",LocateTargSeqName(TargID),BinLen,BinStart + 1,BinCnts); // Wiggle uses 1-start coordinate system

return((int)BinCnts);
}

uint32_t		// returned sequence name identifier, 0 if unable to accept this chromosome name
CSegHaplotypes::AddTargSeqName(char* pszSeqName,	// associate unique identifier with this sequence name
								uint32_t SeqLen)			// sequence is this length - SeqLen associated to sequence identifier will be maximum of all SeqLens specified for this pszSeqName
{
uint32_t SeqNameIdx;
int SeqNameLen;
tsTargSeq *pTargSeq;
char *pszLAname;

if(pszSeqName == NULL || pszSeqName[0] == '\0')
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddTargSeqName: Invalid parameters");
	return(0);
	}

// with any luck the sequence name will be same as the last accessed
if((pTargSeq = LocateTargSeq(m_LASeqNameID)) != NULL)
	{
	pszLAname = &m_szSeqNames[pTargSeq->TargSeqNameOfs];
	if(!stricmp(pszSeqName,pszLAname))
		{
		if(SeqLen > pTargSeq->TargSeqLen)
			pTargSeq->TargSeqLen = SeqLen;
		return(m_LASeqNameID);
		}
	}

// iterate over all known sequence names in case this name to add is a duplicate
for(SeqNameIdx = 0; SeqNameIdx < m_NumSeqNames; SeqNameIdx++)
	{
	pszLAname = &m_szSeqNames[m_TargSeqs[SeqNameIdx].TargSeqNameOfs];
	if(!stricmp(pszSeqName, pszLAname))
		{
		if(SeqLen > m_TargSeqs[SeqNameIdx].TargSeqLen)
			m_TargSeqs[SeqNameIdx].TargSeqLen = SeqLen;
		m_LASeqNameID = m_TargSeqs[SeqNameIdx].TargSeqID;
		return(m_LASeqNameID);
		}
	}

// sequence name is not a duplicate
SeqNameLen = (int)strlen(pszSeqName);
if((m_NxtszSeqNameIdx + SeqNameLen + 1) > (int)sizeof(m_szSeqNames))
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddTargSeqName: Can't accept name '%s', would exceed name space",pszSeqName);
	return(0);
	}
if(m_NumSeqNames == cMaxSHSeqNames)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddTargSeqName: Can't accept name '%s', would exceed limit of %d names",pszSeqName, cMaxSHSeqNames);
	return(0);
	}
pTargSeq = &m_TargSeqs[m_NumSeqNames++];
memset(pTargSeq,0,sizeof(tsTargSeq));
pTargSeq->TargSeqLen = SeqLen;
pTargSeq->TargSeqID = m_NumSeqNames;
pTargSeq->TargSeqNameOfs = m_NxtszSeqNameIdx;
strcpy(&m_szSeqNames[m_NxtszSeqNameIdx], pszSeqName);
m_NxtszSeqNameIdx += SeqNameLen + 1;
m_LASeqNameID = m_NumSeqNames;
return(m_NumSeqNames);
}

tsTargSeq *							// returned sequence detail
CSegHaplotypes::LocateTargSeq(char* pszSeqName)	// for this sequence name
{
uint32_t SeqNameIdx;
tsTargSeq *pTargSeq;
char *pszLAname;

if(pszSeqName == NULL || pszSeqName[0] == '\0')
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "LocateTargSeq: Invalid parameters");
	return(NULL);
	}

// with any luck the sequence name will be same as the last accessed
if((pTargSeq = LocateTargSeq(m_LASeqNameID)) != NULL)
	{
	pszLAname = &m_szSeqNames[pTargSeq->TargSeqNameOfs];
	if(!stricmp(pszSeqName,pszLAname))
		return(pTargSeq);
	}

// iterate over all known sequence names
for(SeqNameIdx = 0; SeqNameIdx < m_NumSeqNames; SeqNameIdx++)
	{
	pszLAname = &m_szSeqNames[m_TargSeqs[SeqNameIdx].TargSeqNameOfs];
	if(!stricmp(pszSeqName, pszLAname))
		{
		m_LASeqNameID = m_TargSeqs[SeqNameIdx].TargSeqID;
		return(&m_TargSeqs[SeqNameIdx]);
		}
	}
return(NULL);
}


char *												// returned sequence name
CSegHaplotypes::LocateTargSeqName(uint32_t SeqNameID)	// identifier returned by call to AddTargSeqName
{
if(SeqNameID < 1 || SeqNameID > m_NumSeqNames)
	return(NULL);
return(&m_szSeqNames[m_TargSeqs[SeqNameID-1].TargSeqNameOfs]);
}

tsTargSeq *										// returned sequence
CSegHaplotypes::LocateTargSeq(uint32_t SeqNameID)	// identifier returned by call to AddTargSeqName
{
if(SeqNameID < 1 || SeqNameID > m_NumSeqNames)
	return(NULL);
return(&m_TargSeqs[SeqNameID-1]);
}


uint32_t		// returned founder name identifier, 0 if unable to accept this founder name
CSegHaplotypes::AddFounder(char* pszFounder) // associate unique identifier with this founder name
{
uint32_t FounderNameIdx;
int FounderNameLen;
char *pszLAFounder;

// with any luck the founder name will be same as the last accessed
if((pszLAFounder = LocateFounder(m_LAFounderID)) != NULL)
	if(!stricmp(pszFounder,pszLAFounder))
		return(m_LAFounderID);

// iterate over all known found names in case this name to add is a duplicate
for(FounderNameIdx = 0; FounderNameIdx < m_NumFounders; FounderNameIdx++)
	if(!stricmp(pszFounder, &m_szFounders[m_szFounderIdx[FounderNameIdx]]))
		{
		m_LAFounderID = FounderNameIdx + 1;
		return(m_LAFounderID);
		}

// founder name is not a duplicate
FounderNameLen = (int)strlen(pszFounder);
if((m_NxtszFounderIdx + FounderNameLen + 1) > (int)sizeof(m_szFounders))
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddFounder: Can't accept name '%s', would exceed name space",pszFounder);
	return(0);
	}

if(m_NumFounders == cMaxSHFounders)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddFounder: Can't accept name '%s', would exceed limit of %d names",pszFounder, cMaxSHFounders);
	return(0);
	}

m_szFounderIdx[m_NumFounders++] = m_NxtszFounderIdx;
strcpy(&m_szFounders[m_NxtszFounderIdx], pszFounder);
m_NxtszFounderIdx += FounderNameLen + 1;
m_LAFounderID = m_NumFounders;
return(m_NumFounders);
}

int		// returned readset identifier, < 1 if unable to accept this readset name
CSegHaplotypes::AddReadset(char* pszReadset) // associate unique identifier with this readset name
{
int ReadsetIdx;
int ReadsetNameLen;

// iterate over all known readset in case this readset to add is a duplicate
for(ReadsetIdx = 0; ReadsetIdx < m_NumReadsetNames; ReadsetIdx++)
	if(!stricmp(pszReadset, &m_szReadsetNames[m_szReadsetIdx[ReadsetIdx]]))
		return(ReadsetIdx + 1);

// readset is not a duplicate
ReadsetNameLen = (int)strlen(pszReadset);
if((m_NxtszReadsetIdx + ReadsetNameLen + 1) > sizeof(m_szReadsetNames))
	return(eBSFerrMaxEntries);
if(m_NumReadsetNames == cMaxSHParents)
	return(eBSFerrMaxEntries);

m_szReadsetIdx[m_NumReadsetNames] = m_NxtszReadsetIdx;
strcpy(&m_szReadsetNames[m_NxtszReadsetIdx], pszReadset);
m_NxtszReadsetIdx += ReadsetNameLen + 1;
return(++m_NumReadsetNames);
}

int 
CSegHaplotypes::ProcessSnpmarkersSNPs(char *pszSNPMarkersFile)	// file containing SNP marker loci
{
	int Rslt;
	int NumFields;
	uint32_t ParentIdx;
	uint32_t EstNumSNPs;
	uint32_t NumSNPsParsed;
	uint32_t RowNumber;
	int ExpNumFields;

	tsSHParent* pParent;

	int NumParentsWithCnts;
	uint32_t SiteSeqID;
	uint32_t SiteLoci;
	etSeqBase SiteBase[cMaxSHParents];

	EstNumSNPs = m_pSNPMarkerCSV->EstNumRows();
	NumSNPsParsed = 0;
	RowNumber = 0;
	ExpNumFields = 0;
	m_NumParents = 0;

	memset(SiteBase,eBaseN,sizeof(SiteBase));
	while((Rslt = m_pSNPMarkerCSV->NextLine()) > 0)				// onto next line containing fields
		{
		RowNumber += 1;
		NumFields = m_pSNPMarkerCSV->GetCurFields();
		if(NumFields < cSSNPMarkerNumFields)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Number of fields (%d) at row %d does not match expected fields for a snpmarkers file '%s'", NumFields, RowNumber, pszSNPMarkersFile, RowNumber);
			Reset();
			return(eBSFerrFieldCnt);
			}
		if(ExpNumFields && ExpNumFields != NumFields)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Inconsistency in number of fields, previously %d but %d fields parsed from '%s' near line %d", ExpNumFields, NumFields, pszSNPMarkersFile, RowNumber);
			Reset();
			return(eBSFerrFieldCnt);
			}
		if(!ExpNumFields)
			{
			m_NumParents = (NumFields - 4) / 9;
			if(((m_NumParents * 9) + 4) != NumFields)
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Invalid number (%d) of fields parsed from '%s' near line %d, was this file generated by 'ngskit4b snpmarkers'?", NumFields, pszSNPMarkersFile, RowNumber);
				Reset();
				return(eBSFerrFieldCnt);
				}
			if(m_NumParents > cMaxSHParents)
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Too many parents (%d) in '%s', max allowed is %d", m_NumParents, pszSNPMarkersFile, cMaxSHParents);
				Reset();
				return(eBSFerrFieldCnt);
				}
			ExpNumFields = NumFields;
			}

		if(RowNumber == 1)
			{
			if(m_pSNPMarkerCSV->IsLikelyHeaderLine())
				{
				if(!NumSNPsParsed)
				{
				// parse out target assembly name against which alignments were made and the parent names
					char* pSrc;
					char* pDst;
					char Chr;
					int Len;

					m_pSNPMarkerCSV->GetText(1, &pSrc);
					pDst = m_ParentsCommonAncestor;
					Len = 0;
					while(Len < sizeof(m_ParentsCommonAncestor) - 1 && (Chr = *pSrc++) && Chr != ':')
						{
						*pDst++ = Chr;
						Len++;
						*pDst = '\0';
						}
					memset(m_Parents, 0, sizeof(m_Parents));
					for(ParentIdx = 0; ParentIdx < m_NumParents; ParentIdx++)
						{
						m_pSNPMarkerCSV->GetText(5 + (ParentIdx * 9), &pSrc);
						pDst = m_Parents[ParentIdx].szName;
						Len = 0;
						while(Len < sizeof(m_Parents[ParentIdx].szName) - 1 && (Chr = *pSrc++) && Chr != ':')
							{
							*pDst++ = Chr;
							Len++;
							*pDst = '\0';
							}
						if((m_Parents[ParentIdx].FounderID = AddFounder(m_Parents[ParentIdx].szName))==0)
							{
							gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to add founder '%s'", m_Parents[ParentIdx].szName);
							Reset();
							return(eBSFerrSpecies);
							}
						}
					}
				}
			else
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Expected CSV file '%s' first line to be a header line with  fully quoted field names", pszSNPMarkersFile);
				Reset();
				return(eBSFerrFieldCnt);
				}
			continue;
			}

		char* pszTxt;
		m_pSNPMarkerCSV->GetText(1, &pszTxt);
		if((SiteSeqID = AddTargSeqName(pszTxt)) == 0)
			{
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "Failed to add identifier for sequence '%s'", pszTxt);
			Reset();
			return(eBSFerrMaxChroms);
			}
		m_pSNPMarkerCSV->GetInt(2, (int*)&SiteLoci);			// get "Loci"
		m_pSNPMarkerCSV->GetInt(4, (int *)&NumParentsWithCnts);	// how many parents do have a representative base - must be same as number of parents in snpmarkers generated file
		if(NumParentsWithCnts != m_NumParents)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Expected number of parents with called bases (%d) to be %d in CSV file '%s' near line %d", NumParentsWithCnts, m_NumParents, pszSNPMarkersFile, RowNumber);
			Reset();
			return(eBSFerrNumRange);
			}

		int FieldIdx = 6;
		pParent = m_Parents;
		for(ParentIdx = 0; ParentIdx < m_NumParents; ParentIdx++, pParent++)
			{
			m_pSNPMarkerCSV->GetText(FieldIdx++, &pszTxt);			// get parent ":Base"
			switch(*pszTxt) {
				case 'a': case 'A':
					SiteBase[ParentIdx] = eBaseA;
					break;
				case 'c': case 'C':
					SiteBase[ParentIdx] = eBaseC;
					break;
				case 'g': case 'G':
					SiteBase[ParentIdx] = eBaseG;
					break;
				case 't': case 'T': case 'u': case 'U':	// U in case RNA alignments..
					SiteBase[ParentIdx] = eBaseT;
					break;
				default:
					gDiagnostics.DiagOut(eDLFatal, gszProcName, "Expected SNP TargBase ('%s') to be one of 'ACGT' in CSV file '%s' near line %d", pszTxt, pszSNPMarkersFile, RowNumber);
					Reset();
					return(eBSFerrFieldCnt);
				}
			FieldIdx+=8;
			}

		// have all parents and their counts at the SNP loci
		// iterate the parents and for every parent determine which bases are allelic and count occurrences over all parents
	if(!AddSNPSite(SiteSeqID,SiteLoci,SiteBase[0],SiteBase[1],SiteBase[2],SiteBase[3]))
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Failed adding SNP site from CSV file '%s' near line %d", pszSNPMarkersFile, RowNumber);
		Reset();
		return(eBSFerrEntry);
		}
	}

m_pSNPMarkerCSV->Close();
delete m_pSNPMarkerCSV;
m_pSNPMarkerCSV = NULL;

if(m_UsedSNPSites > 1)
	qsort(m_pAllocSNPSites, m_UsedSNPSites, sizeof(tsSHSNPSSite), SortSNPSites);

return(eBSFSuccess);
}

size_t 
CSegHaplotypes::AddSNPSite( uint32_t SiteSeqID,		// alignments were to this sequence - could be a chrom/contig/transcript
			  uint32_t SiteLoci,			// loci within target sequence at which SNP for parents were observed
			  etSeqBase P1Base,			// parent 1 base
			  etSeqBase P2Base,			// parent 2 base
			  etSeqBase P3Base,			// parent 3 base
			  etSeqBase P4Base)			// parent 4 base
{
tsSHSNPSSite *pSNPSite;

if(m_pAllocSNPSites == NULL)		// initial allocation?
	{
	size_t memreq = (size_t)cAllocSNPSites * sizeof(tsSHSNPSSite);

#ifdef _WIN32
	m_pAllocSNPSites = (tsSHSNPSSite *) malloc(memreq);	// initial and perhaps the only allocation
	if(m_pAllocSNPSites == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddLoci: Memory allocation of %zd bytes - %s",(int64_t)memreq,strerror(errno));
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	m_pAllocSNPSites = (tsSHSNPSSite *)mmap(NULL,memreq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pAllocSNPSites == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddLoci: Memory allocation of %zd bytes through mmap()  failed - %s",(int64_t)memreq,strerror(errno));
		m_pAllocSNPSites = NULL;
		return(eBSFerrMem);
		}
#endif
	m_AllocMemSNPSites = memreq;
	m_UsedSNPSites = 0;
	m_AllocSNPSites = cAllocSNPSites;
	}
else
	{
	if(m_AllocSNPSites <= (m_UsedSNPSites + 100))	// play safe when increasing allocation
		{
		size_t memreq;
		size_t AllocTo;
		AllocTo = (size_t)cAllocSNPSites + m_AllocMemSNPSites;
		memreq = AllocTo * sizeof(tsSHSNPSSite);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"AddLoci: memory re-allocation to %zd from %zd bytes",(int64_t)memreq,(int64_t)m_AllocMemSNPSites);

#ifdef _WIN32
		pSNPSite = (tsSHSNPSSite *) realloc(m_pAllocSNPSites,memreq);
#else
		pSNPSite = (tsSHSNPSSite *)mremap(m_pAllocSNPSites,m_AllocMemSNPSites,memreq,MREMAP_MAYMOVE);
		if(pSNPSite == MAP_FAILED)
			pSNPSite = NULL;
#endif
		if(pSNPSite == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddLoci: Memory re-allocation to %zd bytes - %s",(int64_t)memreq,strerror(errno));
			return(eBSFerrMem);
			}
		m_pAllocSNPSites = pSNPSite;
		m_AllocMemSNPSites = memreq;
		m_AllocSNPSites = AllocTo; 
		}
	}

pSNPSite = &m_pAllocSNPSites[m_UsedSNPSites++];

pSNPSite->SiteId = m_UsedSNPSites;
pSNPSite->SiteSeqID = SiteSeqID;
pSNPSite->SiteLoci = SiteLoci;
pSNPSite->ParentBase[0] = P1Base;
pSNPSite->ParentBase[1] = P2Base;
pSNPSite->ParentBase[2] = P3Base;
pSNPSite->ParentBase[3] = P4Base;
return(m_UsedSNPSites);
}

uint32_t		// returned founder name identifier, 0 if unable to accept this founder name
CSegHaplotypes::LocateFounder(char* pszFounder) // associate unique identifier with this founder name
{
uint32_t FounderNameIdx;
char *pszLAFounder;

// with any luck the founder name will be same as the last accessed
if((pszLAFounder = LocateFounder(m_LAFounderID)) != NULL)
	if(!stricmp(pszFounder,pszLAFounder))
		return(m_LAFounderID);

// iterate over all known founder names
for(FounderNameIdx = 0; FounderNameIdx < m_NumFounders; FounderNameIdx++)
	if(!stricmp(pszFounder, &m_szFounders[m_szFounderIdx[FounderNameIdx]]))
		{
		m_LAFounderID = FounderNameIdx + 1;
		return(m_LAFounderID);
		}
return(0);		// founder is unknown
}


char*							// returned founder name
CSegHaplotypes::LocateFounder(uint32_t FounderID)	// identifier returned by call to AddFounderName
{
if(FounderID < 1 || FounderID > m_NumFounders)
	return(NULL);
return(&m_szFounders[m_szFounderIdx[FounderID-1]]);
}

tsSHSNPSSite *
CSegHaplotypes::LocateSNPSite(uint32_t ChromID, uint32_t Loci)	// locate SNP site using a binary search over sorted sites
{
tsSHSNPSSite* pEl1;
int32_t IdxLo;
int32_t MidIdx;
int32_t IdxHi;

if(m_UsedSNPSites == 0)
	return(NULL);

if(m_UsedSNPSites < 10)				// if just a few then do a simple linear search
	{
	pEl1 = m_pAllocSNPSites;
	for(IdxLo = 0; IdxLo < (int32_t)m_UsedSNPSites; IdxLo++,pEl1++)
		if(pEl1->SiteSeqID == ChromID && pEl1->SiteLoci == Loci)
			return(pEl1);
	return(NULL);
	}

IdxLo = 0;
IdxHi = (int32_t)m_UsedSNPSites-1;
do {
	MidIdx = (IdxHi + IdxLo)/2;
	pEl1 = &m_pAllocSNPSites[MidIdx];
	if(pEl1->SiteSeqID == ChromID)
		{
		if(pEl1->SiteLoci == Loci)
			return(pEl1);
		if(pEl1->SiteLoci > Loci)
			IdxHi = MidIdx - 1;
		else
			IdxLo = MidIdx + 1;
		}
	else
		{
		if(pEl1->SiteSeqID > ChromID)
			IdxHi = MidIdx - 1;
		else
			IdxLo = MidIdx + 1;
		}
	}
while(IdxHi >= IdxLo);
return(NULL);
}

// SortSAMFounderTargLoci
// Sort m_pSAMloci by ascending Founder.Targ.Loci
int
CSegHaplotypes::SortSAMFounderTargLoci(const void* arg1, const void* arg2)
{
tsSHSAMloci* pEl1 = (tsSHSAMloci*)arg1;
tsSHSAMloci* pEl2 = (tsSHSAMloci*)arg2;

if(pEl1->FounderID < pEl2->FounderID)
	return(-1);
if(pEl1->FounderID > pEl2->FounderID)
	return(1);

if(pEl1->TargID < pEl2->TargID)
	return(-1);
if(pEl1->TargID > pEl2->TargID)
	return(1);

if(pEl1->TargLoci < pEl2->TargLoci)
	return(-1);
if(pEl1->TargLoci > pEl2->TargLoci)
	return(1);
return(0);
}

// SortSAMTargLoci
// Sort m_pSAMloci by ascending Targ.Loci
int
CSegHaplotypes::SortSAMTargLoci(const void* arg1, const void* arg2)
{
tsSHSAMloci* pEl1 = (tsSHSAMloci*)arg1;
tsSHSAMloci* pEl2 = (tsSHSAMloci*)arg2;

if(pEl1->TargID < pEl2->TargID)
	return(-1);
if(pEl1->TargID > pEl2->TargID)
	return(1);

if(pEl1->TargLoci < pEl2->TargLoci)
	return(-1);
if(pEl1->TargLoci > pEl2->TargLoci)
	return(1);
return(0);
}

// SortSNPSites
// 	// Sort SNPSites by ascending SiteSeqID.SiteLoci
int CSegHaplotypes::SortSNPSites(const void* arg1, const void* arg2)
{
tsSHSNPSSite* pEl1 = (tsSHSNPSSite*)arg1;
tsSHSNPSSite* pEl2 = (tsSHSNPSSite*)arg2;

if(pEl1->SiteSeqID < pEl2->SiteSeqID)
	return(-1);
if(pEl1->SiteSeqID > pEl2->SiteSeqID)
	return(1);

if(pEl1->SiteLoci < pEl2->SiteLoci)
	return(-1);
if(pEl1->SiteLoci > pEl2->SiteLoci)
	return(1);
return(0);
}


