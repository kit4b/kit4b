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


int CallHaplotypes(eModeSH PMode,	// processing mode 0: call observed haplotype counts, mode 1: call Skim haplotypes
			int NChroms,			// debugging - limit processing: 0 All chromosomes, 1..n 1st to nTh chromosome
			char *pszTrackName,		// track name
			char *pszTrackDescr,	// track descriptor
			int KmerSize,			// aligning founders using this Kmer size
			int AcceptKmerMatchPerc, // only accept as Kmer match if this percentage of base alleles match 
			int AccumBinSize,		// counts are accumulated into this sized (bp) non-overlapping bins
			int WinMeanBins,		// derive haplotype calls from mean of counts in a sliding window containing this many bins
			bool bTrimmedMean,		// mean is a trimmed mean with highest and lowest bin counts in sliding window removed
			bool bIsMonoploid,		// false: diploid, true: monoploid
			int NumFounderInputFiles,	// number of input founder file specs
			char *pszFounderInputFiles[],	// names of input founder PBA files (wildcards allowed)
			int NumSkimInputFiles,	// number of input skim file specs
			char* pszSkimInputFiles[],		// names of input skim PBA files (wildcards allowed)
			char *pszROIInFile,		// defining regions of input file
			char *pszROIOutFile,	// Regions of interest haplotype calls output file (CSV format)
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

	 eModeSH PMode;					// processing mode
	 int NChroms;					// Limit processing: 0 All chromosomes, 1..n 1st to nTh chromosome
	 int KmerSize;					// aligning founders using this Kmer size
	 int AcceptKmerMatchPerc;		// only accept as Kmer match if this percentage of base alleles match
	 int AccumBinSize;				// counts are accumulated into this sized (Kbp) non-overlapping bins
	 int WinMeanBins;				// derive haplotype calls from mean of counts in a sliding window containing this many bins
	bool bTrimmedMean;				// mean is a trimmed mean with highest and lowest bin counts in sliding window removed
	bool bIsMonoploid;		// false: diploid, true: monoploid

	 char szTrackName[_MAX_PATH];	// track name
	 char szTrackDescr[_MAX_PATH];	// track description
	 int NumFounderInputFiles;		// number of input founder BPA files
	 char *pszFounderInputFiles[cMaxFounderReadsets];		// names of input founder BPA files (wildcards allowed)
	 int NumSkimInputFiles;		// number of input skim BPA files
	 char *pszSkimInputFiles[cMaxSkimReadsets];		// names of input skim BPA files (wildcards allowed)

	 char szOutFile[_MAX_PATH];		// Windowed haplotype calls output file base name, skim readset identifier is appended to this base name
	 char szROIInFile[_MAX_PATH];	// Regions of interest haplotype calls input file (BED format)
	 char szROIOutFile[_MAX_PATH];	// Regions of interest haplotype calls output file (CSV format)

	struct arg_lit *help = arg_lit0 ("h", "help", "print this help and exit");
	struct arg_lit *version = arg_lit0 ("v", "version,ver", "print version information and exit");
	struct arg_int *FileLogLevel = arg_int0 ("f", "FileLogLevel", "<int>", "Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
	struct arg_file *LogFile = arg_file0 ("F", "log", "<file>", "diagnostics log file");

	struct arg_int *pmode = arg_int0 ("m", "mode", "<int>", "processing mode: 0 generate haplotype bin counts over all chromosomes using skim BPAs compared with founder BPAs");
	struct arg_int *nchroms = arg_int0 ("n", "nchroms", "<int>","Limit processing: 0 All chromosomes, 1..n 1st to nTh chromosome");
	struct arg_int *kmersize = arg_int0 ("k", "kmersize", "<int>", "allele alignments using this sized Kmer (defaults to 100bp)");
	struct arg_int *kmermatchperc = arg_int0 ("K", "kmermatchperc", "<int>", "only accept as Kmer match if at least this percentage of base alleles match (defaults to 97)");
	struct arg_int *accumbinsize = arg_int0 ("b", "binsize", "<int>", "counts accumulated into non-overlapping bins of this size in Kbp (defaults to 100Kbp)");
	struct arg_int *winmeanbins = arg_int0 ("B", "binsmean", "<int>", "derive haplotype calls from mean of counts in a sliding window containing this many bins (defaults to 10)");
	struct arg_lit *trimmedmeans = arg_lit0 ("t", "trimmedmeans", "mean is a trimmed mean with highest and lowest bin counts in sliding window removed");
	struct arg_lit *monoploid = arg_lit0 ("p", "ploidy", "haplotypes are for a monoploid, default is for a diploid");
	struct arg_str *trackname = arg_str0("t","trackname","<str>","BED Track name");
	struct arg_str *trackdescr = arg_str0("d","trackdescr","<str>","BED Track description");
	struct arg_file *founderfiles = arg_filen("I", "founderfiles", "<file>", 1,cMaxFounderReadsets,"founder input BPA file(s), wildcards allowed, limit of 32 founders supported");
	struct arg_file *skimfiles = arg_filen("i", "inskimfile", "<file>",1, cMaxSkimReadsets, "skim input BPA file(s), wildcards allowed, limit of 100 skim filespecs supported");
	struct arg_file *roiinfile = arg_file0("r", "inroi", "<file>", "Regions of interest input file (BED format)");
	struct arg_file *roioutfile = arg_file0("R", "outroi", "<file>", "Regions of interest haplotype calls output file (CSV format)");
	struct arg_file *outfile = arg_file1("o", "out", "<file>", "Windowed haplotype calls output file (CSV format)");
	struct arg_int *threads = arg_int0("T","threads","<int>","number of processing threads 0..32 (defaults to 0 which limits threads to maximum of 32 CPU cores)");
	struct arg_end *end = arg_end (200);

	void *argtable[] = { help,version,FileLogLevel,LogFile,
						pmode,monoploid,nchroms,kmersize,kmermatchperc,accumbinsize,winmeanbins,trimmedmeans,trackname,trackdescr,skimfiles,founderfiles,roiinfile,roioutfile,outfile,threads,end };

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

		PMode = pmode->count ? (eModeSH)pmode->ival[0] : eMSHDefault;
		if(PMode < eMSHDefault || PMode > eMSHDefault)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Currently only the one processing mode is supported '-m0'\n");
			exit(1);
			}
		NChroms = nchroms->count ? nchroms->ival[0] : 0;
		if(NChroms < 0)
			NChroms = 0;
		else
			if(NChroms > 1000)
				NChroms = 1000;

		KmerSize = kmersize->count ? kmersize->ival[0] : cDfltKmerSize;
		if(KmerSize < 50 || KmerSize >= 1000)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Kmer size '-k%d' specified outside of range 50.1000bp\n", KmerSize);
			exit(1);
			}

		AcceptKmerMatchPerc = kmermatchperc->count ? kmermatchperc->ival[0] : cDfltAcceptKmerMatchPerc;
		if(AcceptKmerMatchPerc < 90 || AcceptKmerMatchPerc > 100)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Minimum Kmer allele match percentage '-K%d' specified outside of range 90.100\n", AcceptKmerMatchPerc);
			exit(1);
			}

		AccumBinSize = accumbinsize->count ? accumbinsize->ival[0] : cDfltAccumBinSize;
		if(AccumBinSize < 10 || AccumBinSize >= 1000)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Accumulating non-overlapping bin size '-b%dKbp' specified outside of range 10..1000 Kbp\n", AccumBinSize);
			exit(1);
			}

		WinMeanBins = winmeanbins->count ? winmeanbins->ival[0] : cDfltNumMeanBins;
		if(WinMeanBins < 1 || WinMeanBins >= 100)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Mean of counts in '-B%d' bins specified outside of range 1..100\n", WinMeanBins);
			exit(1);
			}

		bTrimmedMean = trimmedmeans->count ? true : false;
		if(bTrimmedMean && WinMeanBins < 5)
			{
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "Warning: Mean of counts over '-t%d' bins must be at least 10 to use trimmed means, defaulting to not using trimmed means\n", WinMeanBins);
			bTrimmedMean = false;
			}

		bIsMonoploid = monoploid->count ? true : false;

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

	for(NumSkimInputFiles = Idx = 0; NumSkimInputFiles < cMaxSkimReadsets && Idx < skimfiles->count; Idx++)
		{
		pszSkimInputFiles[Idx] = NULL;
		if(pszSkimInputFiles[NumSkimInputFiles] == NULL)
			pszSkimInputFiles[NumSkimInputFiles] = new char[_MAX_PATH];
		strncpy(pszSkimInputFiles[NumSkimInputFiles], skimfiles->filename[Idx], _MAX_PATH);
		pszSkimInputFiles[NumSkimInputFiles][_MAX_PATH - 1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(pszSkimInputFiles[NumSkimInputFiles]);
		if(pszSkimInputFiles[NumSkimInputFiles][0] != '\0')
			NumSkimInputFiles++;
		}

	if(!NumSkimInputFiles)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: After removal of whitespace, there are no input skim file(s) specified with '-i<filespec>' option)\n");
		exit(1);
		}


	for(NumFounderInputFiles = Idx = 0; NumFounderInputFiles < cMaxFounderReadsets && Idx < founderfiles->count; Idx++)
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

	if(!NumFounderInputFiles)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: After removal of whitespace, there are no input founder file(s) specified with '-I<filespec>' option)\n");
		exit(1);
		}

	if(roiinfile->count)
		{
		strcpy (szROIInFile, roiinfile->filename[0]);
		CUtility::TrimQuotedWhitespcExtd (szROIInFile);
		if(szROIInFile[0] == '\0')
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: After removal of whitespace, no input ROI file specified with '-r<filespec>' option)\n");
			exit(1);
			}
		if(!roioutfile->count)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: ROI input file specified but no output file specified with '-R<filespec>' option)\n");
			exit(1);
			}
		}
	else
		szROIInFile[0] = '\0';

	if(roioutfile->count)
		{
		strcpy (szROIOutFile, roioutfile->filename[0]);
		CUtility::TrimQuotedWhitespcExtd (szROIOutFile);
		if(szROIOutFile[0] == '\0')
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: After removal of whitespace, no output file specified with '-R<filespec>' option)\n");
			exit(1);
			}
		}
	else
		szROIOutFile[0] = '\0';

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
			pszDescr = "Calling haplotypes in skim PBAs comparing with founder PBAs";
			break;
		}

	if(NChroms > 0)
		gDiagnostics.DiagOutMsgOnly (eDLInfo, "Debugging - only processing 1st %d chromosomes",NChroms); 
	gDiagnostics.DiagOutMsgOnly (eDLInfo, "Kmer size : %dbp", KmerSize);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Only accepting as Kmer match if at least this percentage of base alleles match: %d",AcceptKmerMatchPerc);
	gDiagnostics.DiagOutMsgOnly (eDLInfo, "Accumulating counts in non-overlapping bins of size : %dKbp", AccumBinSize);
	gDiagnostics.DiagOutMsgOnly (eDLInfo, "Using mean of counts in sliding window containing this many bins : %d", WinMeanBins);
	gDiagnostics.DiagOutMsgOnly (eDLInfo, "Using trimmed means : %s", bTrimmedMean ? "Yes" : "No");
	gDiagnostics.DiagOutMsgOnly (eDLInfo, "Haploids for : %s", bIsMonoploid ? "Monoploid" : "Diploid");
	gDiagnostics.DiagOutMsgOnly (eDLInfo, "Segment haplotypes : '%s'", pszDescr);
	gDiagnostics.DiagOutMsgOnly (eDLInfo, "Track name : '%s'", szTrackName);
	gDiagnostics.DiagOutMsgOnly (eDLInfo, "Track name : '%s'", szTrackDescr);
	
	for(Idx = 0; Idx < NumFounderInputFiles; Idx++)
		gDiagnostics.DiagOutMsgOnly (eDLInfo, "Founder file : '%s'", pszFounderInputFiles[Idx]);

	for(Idx = 0; Idx < NumSkimInputFiles; Idx++)
		gDiagnostics.DiagOutMsgOnly (eDLInfo, "Skim file : '%s'", pszSkimInputFiles[Idx]);

	if(szROIInFile[0] != '\0')
		gDiagnostics.DiagOutMsgOnly (eDLInfo, "Regions of interest loaded from : '%s'", szROIInFile);
	if(szROIOutFile[0] != '\0')
		gDiagnostics.DiagOutMsgOnly (eDLInfo, "Haplotypes in regions of interest written to : '%s'", szROIOutFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"number of threads : %d",NumThreads);

#ifdef _WIN32
	SetPriorityClass (GetCurrentProcess (), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start ();
	Rslt = 0;
	Rslt = CallHaplotypes(PMode,			// processing mode
					NChroms,				// debugging - limit processing: 0 All chromosomes, 1..n 1st to nTh chromosome
					szTrackName,			// track name
					szTrackDescr,			// track descriptor
					KmerSize,				// aligning founders using this Kmer size
					AcceptKmerMatchPerc,	// only accept as Kmer match if this percentage of base alleles match
					AccumBinSize*1000,		// now use bp rather than Kbp! Counts are accumulated into this sized non-overlapping bins
					WinMeanBins,			// derive haplotype calls from mean of counts in a sliding window containing this many bins
					bTrimmedMean,			// mean is a trimmed mean with highest and lowest bin counts in sliding window removed
					bIsMonoploid,			// false: diploid, true: monoploid
					NumFounderInputFiles,	// number of input founder file specs
					pszFounderInputFiles,	// names of input founder PBA files (wildcards allowed)
					NumSkimInputFiles,		// number of input skim file specs
					pszSkimInputFiles,		// names of input skim PBA files (wildcards allowed)
					szROIInFile,			// defining regions of input file
					szROIOutFile,			// Regions of interest haplotype calls output file (CSV format)
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

int CallHaplotypes(eModeSH PMode,	// processing mode
			int NChroms,			// debugging - limit processing: 0 All chromosomes, 1..n 1st to nTh chromosome
			char *pszTrackName,		// track name
			char *pszTrackDescr,	// track descriptor
			int KmerSize,			// aligning founders using this Kmer size
			int AcceptKmerMatchPerc, // only accept as Kmer match if this percentage of base alleles match
			int AccumBinSize,		// accumulating counts into this sized (bp) non-overlapping bins
			int WinMeanBins,		// derive haplotype calls from mean of counts in a sliding window containing this many bins
			bool bTrimmedMean,		// mean is a trimmed mean with highest and lowest bin counts in sliding window removed
			bool bIsMonoploid,		// false: diploid, true: monoploid
			int NumFounderInputFiles,	// number of input founder file specs
			char *pszFounderInputFiles[],	// names of input founder PBA files (wildcards allowed)
			int NumSkimInputFiles,	// number of input skim file specs
			char* pszSkimInputFiles[],		// names of input skim PBA files (wildcards allowed)
			char *pszROIInFile,		// defining regions of input file
			char *pszROIOutFile,	// Regions of interest haplotype calls output file (CSV format)
			char* pszOutFile,		// Windowed haplotype calls output file (CSV format)
			int NumThreads)			// number of worker threads to use
{
int Rslt;
CCallHaplotypes* pCallHaplotypes;

if((pCallHaplotypes = new CCallHaplotypes) == NULL)
	{
	gDiagnostics.DiagOut (eDLFatal, gszProcName, "Unable to instantiate instance of CCallHaplotypes");
	return(eBSFerrObj);
	}
Rslt = pCallHaplotypes->Process(PMode,NChroms,pszTrackName,pszTrackDescr,KmerSize,AcceptKmerMatchPerc,AccumBinSize,WinMeanBins,bTrimmedMean,bIsMonoploid,NumFounderInputFiles,pszFounderInputFiles,NumSkimInputFiles,pszSkimInputFiles,pszROIInFile,pszROIOutFile,pszOutFile,NumThreads);
delete pCallHaplotypes;
return(Rslt);
}


CCallHaplotypes::CCallHaplotypes()
{
m_pChromMetadata = NULL;
m_pWinBinCnts = NULL;
m_pWorkQueueEls = NULL;
m_pInBuffer = NULL;
m_pOutBuffer = NULL;
m_hOutFile = -1;
m_hInFile = -1;
m_hROIInFile = -1;
m_hROIOutFile = -1;
m_bMutexesCreated = false;
Reset();
}

CCallHaplotypes::~CCallHaplotypes()
{
if(m_hInFile != -1)
	close(m_hInFile);
if(m_hOutFile != -1)
	close(m_hOutFile);

if(m_hROIInFile != -1)
	close(m_hROIInFile);
if(m_hROIOutFile != -1)
	close(m_hROIOutFile);
if(m_pWorkQueueEls != NULL)
	delete []m_pWorkQueueEls;
if(m_pOutBuffer != NULL)
	delete []m_pOutBuffer;
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

if(m_pWinBinCnts != NULL)
	{
#ifdef _WIN32
	free(m_pWinBinCnts);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pWinBinCnts != MAP_FAILED)
		munmap(m_pWinBinCnts, m_AllocdWinBinCntsMem);
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

if(m_hROIInFile != -1)
	{
	close(m_hROIInFile);
	m_hROIInFile = -1;
	}
if(m_hROIOutFile != -1)
	{
	close(m_hROIOutFile);
	m_hROIOutFile = -1;
	}

if(m_pWorkQueueEls != NULL)
	{
	delete []m_pWorkQueueEls;
	m_pWorkQueueEls = NULL;
	}
m_TotWorkQueueEls = 0;
m_NumQueueElsProcessed = 0;

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

if(m_pWinBinCnts != NULL)
	{
#ifdef _WIN32
	free(m_pWinBinCnts);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pWinBinCnts != MAP_FAILED)
		munmap(m_pWinBinCnts, m_AllocdWinBinCntsMem);
#endif
	m_pWinBinCnts = NULL;
	}
m_UsedWinBinCnts = 0;
m_AllocdWinBinCnts = 0;
m_AllocdWinBinCntsMem = 0;

m_NumFndrs = 0;
m_NChroms = 0;
m_ScoreKmerSize = cDfltKmerSize;
m_AcceptKmerMatchPerc = cDfltAcceptKmerMatchPerc;
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

m_NumFndrs = 0;
m_WinMeanBins = cDfltNumMeanBins;
m_bTrimmedMean = false;
m_bIsMonoploid = false;
m_AccumBinSize = cDfltAccumBinSize * 1000;
m_ScoreKmerSize = cDfltKmerSize;

m_pszROIInFile = NULL;
m_pszROIOutFile = NULL;

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
CCallHaplotypes::Process(eModeSH PMode,	// processing mode
			int NChroms,			// debugging - limit processing: 0 All chromosomes, 1..n 1st to nTh chromosome
			char *pszTrackName,		// track name
			char *pszTrackDescr,	// track descriptor
			int KmerSize,			// aligning founders using this Kmer size
			int AcceptKmerMatchPerc, // only accept as Kmer match if this percentage of base alleles match
			uint32_t AccumBinSize,	// accumulate counts into this sized (bp) non-overlapping bins
			int WinMeanBins,		// derive haplotype calls from mean of counts in a sliding window containing this many bins
			bool bTrimmedMean,		// mean is a trimmed mean with highest and lowest bin counts in sliding window removed
			bool bIsMonoploid,		// false: diploid, true: monoploid
			int NumFounderInputFiles,	// number of input founder file specs
			char *pszFounderInputFiles[],	// names of input founder PBA files (wildcards allowed)
			int NumSkimInputFiles,	// number of input skim file specs
			char* pszSkimInputFiles[],		// names of input skim PBA files (wildcards allowed)
			char *pszROIInFile,		// defining regions of input file
			char *pszROIOutFile,	// Regions of interest haplotype calls output file (CSV format)
			char* pszOutFile,		// output to this file
			int NumThreads)			// number of worker threads to use
{
int Rslt;
int NumFiles;
int TotNumFiles;
size_t memreq;

Reset();

CreateMutexes();
m_bIsMonoploid = bIsMonoploid;
m_AccumBinSize = AccumBinSize;
m_WinMeanBins = WinMeanBins;
m_bTrimmedMean = bTrimmedMean;
m_ScoreKmerSize = KmerSize;
m_AcceptKmerMatchPerc = AcceptKmerMatchPerc;
m_NumThreads = NumThreads;
m_NChroms = NChroms;

if(NumFounderInputFiles > cMaxFounderReadsets)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Can accept at most %d founder readsets for processing, %d requested", cMaxFounderReadsets, NumFounderInputFiles);
	Reset();
	return(eBSFerrMem);
	}

if((m_pInBuffer = new uint8_t[cInBuffSize]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Memory allocation of %u bytes for input buffering failed - %s", cInBuffSize, strerror(errno));
	Reset();
	return(eBSFerrMem);
	}
m_AllocInBuff = cInBuffSize;
m_InNumBuffered = 0;

if((m_pOutBuffer = new uint8_t[cOutBuffSize]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Memory allocation of %u bytes for output buffering failed - %s", cOutBuffSize, strerror(errno));
	Reset();
	return(eBSFerrMem);
	}
m_AllocOutBuff = cOutBuffSize;
m_OutBuffIdx = 0;

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
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Assembly: Memory allocation of %lld bytes through mmap() for chromosome metadata failed - %s", (int64_t)memreq, strerror(errno));
	m_pChromMetadata = NULL;
	Reset();
	return(eBSFerrMem);
	}
#endif
m_AllocdChromMetadataMem = memreq;
m_UsedNumChromMetadata = 0;
m_AllocdChromMetadata = cAllocChromMetadata;

memreq = (size_t)cAllocWinBinCnts * sizeof(tsBinCnts);
#ifdef _WIN32
m_pWinBinCnts = (tsBinCnts *)malloc(memreq);	// initial and perhaps the only allocation
if(m_pWinBinCnts == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Initial memory allocation of %lld bytes for for chromosome windowed founder counts failed - %s", (int64_t)memreq, strerror(errno));
	Reset();
	return(eBSFerrMem);
	}
#else
m_pWinBinCnts = (tsBinCnts *)mmap(NULL, memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
if(m_pWinBinCnts == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Assembly: Memory allocation of %lld bytes through mmap() for chromosome windowed founder counts failed - %s", (int64_t)memreq, strerror(errno));
	m_pWinBinCnts = NULL;
	Reset();
	return(eBSFerrMem);
	}
#endif
m_AllocdWinBinCntsMem = memreq;
m_UsedWinBinCnts = 0;
m_AllocdWinBinCnts = cAllocWinBinCnts;

if(pszROIInFile!= NULL && pszROIInFile[0] != '\0' && pszROIOutFile!= NULL && pszROIOutFile[0] != '\0')
	{
	m_pszROIInFile = pszROIInFile;
	m_pszROIOutFile = pszROIOutFile;
	// try parsing the regions of interest which could be either BED or CSV
	etClassifyFileType FileType = CUtility::ClassifyFileType(m_pszROIInFile);
	switch(FileType) {
		case eCFTopenerr:		// unable to open file for reading
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open file: '%s'",pszROIInFile);
			Reset();
			return(eBSFerrOpnFile);

		case eCFTlenerr:		// file length is insufficient to classify type
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Unable to classify file type (insuffient data points): '%s'",pszROIInFile);
			Reset();
			return(eBSFerrFileAccess);

		case eCFTunknown:		// unable to reliably classify
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Unable to reliably classify file type: '%s'",pszROIInFile);
			Reset();
			return(eBSFerrFileType);

		case eCFTCSV:		// CSV loci processing
			Rslt = 0;
			break;

		case eCFTBED:			// UCSC BED processing
			Rslt = 0;
			break;

		default:
			break;
		}
	}


Rslt = eBSFSuccess;		// assume success!
CSimpleGlob glob(SG_GLOB_FULLSORT);
int Idx;
char* pszInFile;
TotNumFiles = 0;
for(Idx = 0; Idx < NumFounderInputFiles; Idx++)	
	{
	glob.Init();
	if(glob.Add(pszFounderInputFiles[Idx]) < SG_SUCCESS)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to glob '%s", pszFounderInputFiles[Idx]);
		Reset();
		return(eBSFerrOpnFile);	// treat as though unable to open file
		}
	if((NumFiles = glob.FileCount()) <= 0)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to locate any founder file matching '%s", pszFounderInputFiles[Idx]);
		Reset();
		return(eBSFerrFileName);
		}
	TotNumFiles += NumFiles;
	if(TotNumFiles > cMaxFounderReadsets)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Can accept at most %d founder readsets for processing, after wildcard file name expansions there are %d requested", cMaxFounderReadsets, TotNumFiles);
		Reset();
		return(eBSFerrMem);
		}

	Rslt = eBSFSuccess;
	for(int FileID = 0; Rslt >= eBSFSuccess && FileID < NumFiles; ++FileID)
		{
		pszInFile = glob.File(FileID);
		if((Rslt = LoadPBAFile(pszInFile,false)) < eBSFSuccess)			
			{
			Reset();
			return(Rslt);
			}
		}
	}
m_NumFndrs = m_NumReadsetNames;
gDiagnostics.DiagOut(eDLInfo, gszProcName, "Completed loading founder %d PBA files",m_NumFndrs);

// next is to individually process each skim PBAs against the founder panel PBAs
for(Idx = 0; Idx < NumSkimInputFiles; Idx++)	
	{
	glob.Init();
	if(glob.Add(pszSkimInputFiles[Idx]) < SG_SUCCESS)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to glob skim '%s", pszSkimInputFiles[Idx]);
		Reset();
		return(eBSFerrOpnFile);	// treat as though unable to open file
		}
	if((NumFiles = glob.FileCount()) <= 0)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to locate any skim PBA file matching '%s", pszSkimInputFiles[Idx]);
		Reset();
		return(eBSFerrOpnFile);	// treat as though unable to open file
		}

	Rslt = eBSFSuccess;
	for(int FileID = 0; Rslt >= eBSFSuccess && FileID < NumFiles; ++FileID)
		{
		pszInFile = glob.File(FileID);
		if((Rslt = ProcessSkimPBAFile(pszInFile,pszOutFile)) < eBSFSuccess)			
			{
			Reset();
			return(Rslt);
			}
		}
	}
Reset();
return(Rslt);
}



int
CCallHaplotypes::ProcessSkimPBAFile(char* pszSkimPBAFile,	// load, process and call haplotypes for this skim PBA file against previously loaded panel founder PBAs
						char *pszRsltsFileBaseName)			// results are written to this file base name with skim readset identifier and type appended
{
int Rslt;
char szOutFile[_MAX_PATH];

// mark state so can be restored ready for next skim PBA to be loaded and processed
uint32_t NumReadsetNames = m_NumReadsetNames;
uint32_t NxtszReadsetIdx = m_NxtszReadsetIdx;
uint32_t NumChromNames = m_NumChromNames;
uint32_t NxtszChromIdx = m_NxtszChromIdx;
uint32_t UsedNumChromMetadata = m_UsedNumChromMetadata;
uint32_t UsedWinBinCnts = m_UsedWinBinCnts;

m_LAChromNameID = 0;
m_LAReadsetNameID = 0;
if((Rslt = LoadPBAFile(pszSkimPBAFile,true)) < eBSFSuccess)
	{
	Reset();
	return(Rslt);
	}
gDiagnostics.DiagOut(eDLInfo, gszProcName, "Completed loading skim PBA file");

char *pszSkimReadset = LocateReadset(m_NumReadsetNames);
sprintf(szOutFile,"%s_Skim.%s.csv",pszRsltsFileBaseName,pszSkimReadset);

#ifdef _WIN32
m_hOutFile = open(szOutFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE));
#else
	if((m_hOutFile = open64(szOutFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE)) != -1)
		if(ftruncate(m_hOutFile,0)!=0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create/truncate %s - %s",szOutFile,strerror(errno));
			Reset();
			return(eBSFerrCreateFile);
			}
#endif
if(m_hOutFile < 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create/truncate %s - %s",szOutFile,strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}

m_Readsets[m_NumReadsetNames-1].StartWinBinCntsIdx = m_UsedWinBinCnts;
m_Readsets[m_NumReadsetNames-1].NumWinBins = 0;
m_Readsets[m_NumReadsetNames-1].LRAWinBinCntIdx = m_UsedWinBinCnts;
Rslt = CountHaplotypes(m_NumReadsetNames,m_NumReadsetNames-1,m_ScoreKmerSize,m_AccumBinSize);
m_Readsets[m_NumReadsetNames-1].NumWinBins = m_UsedWinBinCnts - m_Readsets[m_NumReadsetNames-1].StartWinBinCntsIdx;
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Process: Calling haplotype from bin counts ...");
ChooseHaplotypes();  // choose which set of haplotypes are most probable for each bin
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Process: Haplotype calls completed");
// report
m_OutBuffIdx=sprintf((char *)m_pOutBuffer,"\"Chrom\",\"Loci\",\"Skim-Haplotype:%s\",\"Skim-RawHaplotype:%s\",\"Skim-ConfHaplotype:%s\",\"Skim-NonUniques:%s\",\"Skim-Unaligned:%s\"",pszSkimReadset,pszSkimReadset,pszSkimReadset,pszSkimReadset,pszSkimReadset);
for(uint32_t FounderID = 1; FounderID <= m_NumFndrs; FounderID++)
	m_OutBuffIdx+=sprintf((char *)&m_pOutBuffer[m_OutBuffIdx],",\"FndrUniques-%s:%s\"",pszSkimReadset,LocateReadset(FounderID));
for(uint32_t FounderID = 1; FounderID <= m_NumFndrs; FounderID++)
	m_OutBuffIdx+=sprintf((char *)&m_pOutBuffer[m_OutBuffIdx],",\"FndrNonUniques-%s:%s\"",pszSkimReadset,LocateReadset(FounderID));


tsBinCnts *pSkimWinBinCnts;
uint32_t WndFndrCntsIdx;
uint32_t FounderIdx;
uint32_t NumWinBins = m_Readsets[m_NumReadsetNames-1].NumWinBins;
pSkimWinBinCnts = &m_pWinBinCnts[m_Readsets[m_NumReadsetNames-1].StartWinBinCntsIdx];
for(WndFndrCntsIdx = 0; WndFndrCntsIdx < NumWinBins; WndFndrCntsIdx++, pSkimWinBinCnts++)
	{
	m_OutBuffIdx+=sprintf((char *)&m_pOutBuffer[m_OutBuffIdx],"\n\"%s\",%u,%u,%u,%u,%u,%u",
							LocateChrom(pSkimWinBinCnts->ChromID),pSkimWinBinCnts->StartLoci,pSkimWinBinCnts->HaplotypeClass,pSkimWinBinCnts->RawHaplotypeClass,pSkimWinBinCnts->HaplotypeConf,pSkimWinBinCnts->MultiFounder,pSkimWinBinCnts->NoFounder);
	for(FounderIdx = 1; FounderIdx <= m_NumFndrs; FounderIdx++)
		m_OutBuffIdx+=sprintf((char *)&m_pOutBuffer[m_OutBuffIdx],",%u",pSkimWinBinCnts->Uniques[FounderIdx-1]);
	for(FounderIdx = 1; FounderIdx <= m_NumFndrs; FounderIdx++)
		m_OutBuffIdx+=sprintf((char *)&m_pOutBuffer[m_OutBuffIdx],",%u",pSkimWinBinCnts->NonUniques[FounderIdx-1]);

	if((m_OutBuffIdx + 1000) > m_AllocOutBuff)
		{
		if (!CUtility::RetryWrites(m_hOutFile,m_pOutBuffer,m_OutBuffIdx))
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "ClassifyHaplotypes: Fatal error in RetryWrites()");
			Reset();
			return(eBSFerrFileAccess);
			}
		m_OutBuffIdx = 0;
		}
	}
if(m_OutBuffIdx && m_hOutFile != -1)
	{
	if (!CUtility::RetryWrites(m_hOutFile,m_pOutBuffer,m_OutBuffIdx))
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "ClassifyHaplotypes: Fatal error in RetryWrites()");
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

// now for the BEDs, one for each founder
for(FounderIdx = 1; FounderIdx <= m_NumFndrs; FounderIdx++)
	{
	sprintf(szOutFile,"%s_Skim.%s.F%c.bed",pszRsltsFileBaseName,pszSkimReadset,(char)((uint32_t)'a'+ FounderIdx-1));
	ReportHaplotypesBED(szOutFile,LocateReadset(FounderIdx),pszSkimReadset,0x01 << (FounderIdx-1));
	}

// 
// restore state ready for next skim PBA to be loaded and processed
m_NumReadsetNames = NumReadsetNames;
m_NxtszReadsetIdx = NxtszReadsetIdx;
m_NumReadsetNames = NumReadsetNames;
m_NumChromNames = NumChromNames;
m_NxtszChromIdx = NxtszChromIdx;
m_UsedWinBinCnts = UsedWinBinCnts;
if(m_UsedNumChromMetadata > UsedNumChromMetadata)
	{
	tsChromMetadata *pChromMetadata = &m_pChromMetadata[UsedNumChromMetadata];
	for(uint32_t MetaIdx = UsedNumChromMetadata; MetaIdx < m_UsedNumChromMetadata; MetaIdx++,pChromMetadata++)
		{
		if(pChromMetadata->pPBAs != NULL)
#ifdef _WIN32
			free(pChromMetadata->pPBAs);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
		if(pChromMetadata->pPBAs != MAP_FAILED)
			munmap(pChromMetadata->pPBAs, pChromMetadata->ChromLen);
#endif
		pChromMetadata->pPBAs = NULL;
		}
	}
m_UsedNumChromMetadata = UsedNumChromMetadata;
m_LAChromNameID = 0;
m_LAReadsetNameID = 0;
return(eBSFSuccess);
}

int
CCallHaplotypes::ChooseHaplotypes(void) // iterate over all the window bins and choose which haplotype to call for that bin interpolating those bins which were initially called as being indeterminate
{
uint32_t CurChromID;
uint64_t HMMBinHaplotypes;		// 16 state HMM, each state (called haplotype) packed into 4 bits
tsReadsetMetadata *pSkim;
tsBinCnts *pSkimWinBinCnts;
tsBinCnts *pCurInterpWinBinCnts;
tsBinCnts InterpWinBinCnts;
tsBinCnts InterpWinBinCntsMax;
tsBinCnts InterpWinBinCntsMin;
uint32_t WinIdx;
pSkim = &m_Readsets[m_NumFndrs];	// skim always is last loaded PBA following the founders
pSkimWinBinCnts = &m_pWinBinCnts[pSkim->StartWinBinCntsIdx];

// calculate trimmed means for counts
HMMBinHaplotypes = 0;
CurChromID = 0;
memset(&InterpWinBinCnts,0,sizeof(InterpWinBinCnts));
for(uint32_t SkimBinIdx = 0; SkimBinIdx < pSkim->NumWinBins; SkimBinIdx++,pSkimWinBinCnts++)
	{
	if(CurChromID != pSkimWinBinCnts->ChromID)	// onto a new chrom?
		{
		CurChromID = pSkimWinBinCnts->ChromID;
		HMMBinHaplotypes = 0;
		}
	if(m_WinMeanBins > 1 && (pSkim->NumWinBins - SkimBinIdx) >= 2)
		{
		pCurInterpWinBinCnts = pSkimWinBinCnts;
		memset(&InterpWinBinCnts,0,sizeof(InterpWinBinCnts));
		memset(&InterpWinBinCntsMax,0,sizeof(InterpWinBinCnts));
		memset(&InterpWinBinCntsMin,0xff,sizeof(InterpWinBinCnts));
		for(WinIdx = 0; WinIdx < m_WinMeanBins && CurChromID == pCurInterpWinBinCnts->ChromID; WinIdx++,pCurInterpWinBinCnts++)
			{
			for(uint32_t CntIdx = 0; CntIdx < m_NumFndrs; CntIdx++)
				{
				if(pCurInterpWinBinCnts->Uniques[CntIdx] > InterpWinBinCntsMax.Uniques[CntIdx])
					InterpWinBinCntsMax.Uniques[CntIdx] = pCurInterpWinBinCnts->Uniques[CntIdx];
				if(pCurInterpWinBinCnts->Uniques[CntIdx] < InterpWinBinCntsMin.Uniques[CntIdx])
					InterpWinBinCntsMin.Uniques[CntIdx] = pCurInterpWinBinCnts->Uniques[CntIdx];
				InterpWinBinCnts.Uniques[CntIdx] += pCurInterpWinBinCnts->Uniques[CntIdx];
				}
			if(pCurInterpWinBinCnts->MultiFounder > InterpWinBinCntsMax.MultiFounder)
				InterpWinBinCntsMax.MultiFounder = pCurInterpWinBinCnts->MultiFounder;
			if(pCurInterpWinBinCnts->MultiFounder < InterpWinBinCntsMin.MultiFounder)
				InterpWinBinCntsMin.MultiFounder = pCurInterpWinBinCnts->MultiFounder;
			if(pCurInterpWinBinCnts->NoFounder > InterpWinBinCntsMax.NoFounder)
				InterpWinBinCntsMax.NoFounder = pCurInterpWinBinCnts->NoFounder;
			if(pCurInterpWinBinCnts->NoFounder < InterpWinBinCntsMin.NoFounder)
				InterpWinBinCntsMin.NoFounder = pCurInterpWinBinCnts->NoFounder;
			InterpWinBinCnts.MultiFounder += pCurInterpWinBinCnts->MultiFounder;
			InterpWinBinCnts.NoFounder += pCurInterpWinBinCnts->NoFounder;
			}

		if(m_bTrimmedMean && WinIdx >= m_WinMeanBins) // trimmed accumulated counts (top and bottom removed) only applied if at least this window size, otherwise no trimming applied
			{
			for(uint32_t CntIdx = 0; CntIdx < m_NumFndrs; CntIdx++) 
				{
				InterpWinBinCnts.Uniques[CntIdx] -= InterpWinBinCntsMax.Uniques[CntIdx];
				InterpWinBinCnts.Uniques[CntIdx] -= InterpWinBinCntsMin.Uniques[CntIdx];
				}
			InterpWinBinCnts.MultiFounder -= InterpWinBinCntsMax.MultiFounder;
			InterpWinBinCnts.MultiFounder -= InterpWinBinCntsMin.MultiFounder;
			InterpWinBinCnts.NoFounder -= InterpWinBinCntsMax.NoFounder;
			InterpWinBinCnts.NoFounder -= InterpWinBinCntsMin.NoFounder;
			}
		HMMBinHaplotypes = ChooseHaplotype(HMMBinHaplotypes,&InterpWinBinCnts); // choose haplotype using accumulated counts
		}
	else
		HMMBinHaplotypes = ChooseHaplotype(HMMBinHaplotypes,pSkimWinBinCnts);	// choose haplotype using non-accumulated counts
	
	// backtrack for up to last 16 skim haplotype calls and update
	pCurInterpWinBinCnts = pSkimWinBinCnts;
	uint64_t Haplotypes = HMMBinHaplotypes;
	do {
		pCurInterpWinBinCnts->HaplotypeClass = Haplotypes & 0x0f;
		Haplotypes >>= 4;
		Haplotypes &= 0x0fffffffffffffff;
		pCurInterpWinBinCnts--;
		}
	while(Haplotypes != 0);
	}
return(0);
}

int
CCallHaplotypes::ReportHaplotypesBED(char *pszOutFile,		// BED file
							char *pszFounder,	// this was the founder haplotype
							char *pszSkim,		// against which this skim readset was aligned
					uint32_t HaplotypeMsk)	// reporting on haplotypes matching this mask, 0x01 for Fa, 0x02 for Fb, 0x04 for Fc, 0x010 for Fd
{
uint32_t NumWinBins;
tsBinCnts *pWinBinCnts;

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
m_OutBuffIdx = sprintf((char *)m_pOutBuffer,"track name=\"%s->%s\" description=\"Skim readset %s PBA aligned against Founder %s\" useScore=0\n", pszSkim,pszFounder, pszSkim,pszFounder);

NumWinBins = m_Readsets[m_NumReadsetNames-1].NumWinBins;
pWinBinCnts = &m_pWinBinCnts[m_Readsets[m_NumReadsetNames-1].StartWinBinCntsIdx];
char *pszCurChrom;

uint32_t CurStartLoci;
uint32_t WndFndrCntsIdx;
uint32_t CurChromID = 0;
uint32_t CurEndLoci = 0;
for(WndFndrCntsIdx = 0; WndFndrCntsIdx < NumWinBins; WndFndrCntsIdx++, pWinBinCnts++)
	{
	// simple chrom\t\loci\t\endloci is all that is required
	if(pWinBinCnts->ChromID != CurChromID) // starting a new chrom?
		{
		if(CurEndLoci != 0)
			m_OutBuffIdx+=sprintf((char *)&m_pOutBuffer[m_OutBuffIdx],"%s\t%u\t%u\n",pszCurChrom,CurStartLoci,CurEndLoci);
		CurEndLoci = 0;
		CurChromID = pWinBinCnts->ChromID;
		pszCurChrom = LocateChrom(CurChromID);
		}

	if(pWinBinCnts->HaplotypeClass & HaplotypeMsk)
		{
		if(CurEndLoci == 0)
			CurStartLoci = pWinBinCnts->StartLoci;
		CurEndLoci = pWinBinCnts->EndLoci + 1;
		}
	else // else haplotype not called
		{
		if(CurEndLoci != 0)
			{
			m_OutBuffIdx+=sprintf((char *)&m_pOutBuffer[m_OutBuffIdx],"%s\t%u\t%u\n",pszCurChrom,CurStartLoci,CurEndLoci);
			CurEndLoci = 0;
			}
		}
	if((m_OutBuffIdx + 1000) > m_AllocOutBuff)
		{
		if (!CUtility::RetryWrites(m_hOutFile,m_pOutBuffer,m_OutBuffIdx))
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "ClassifyHaplotypes: Fatal error in RetryWrites()");
			Reset();
			return(eBSFerrFileAccess);
			}
		m_OutBuffIdx = 0;
		}
	}
if(CurEndLoci != 0)
	{
	m_OutBuffIdx+=sprintf((char *)&m_pOutBuffer[m_OutBuffIdx],"%s\t%u\t%u\n",pszCurChrom,CurStartLoci,CurEndLoci);
	CurEndLoci = 0;
	}

if(m_OutBuffIdx && m_hOutFile != -1)
	{
	if (!CUtility::RetryWrites(m_hOutFile,m_pOutBuffer,m_OutBuffIdx))
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "ClassifyHaplotypes: Fatal error in RetryWrites()");
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

tsBinCnts *			// returned WinBinCnts if located in pReadset
CCallHaplotypes::LocateWinBinCnts(tsBinCnts *pWinBinCnts,	// chrom and start loci in this bin are to be located
								tsReadsetMetadata *pReadset)	// in this readset
{
tsBinCnts *pCnts;
uint32_t LRAWinBinCntIdx;
uint32_t WinBinCntIdx;
uint32_t NumWinBins;
if(pReadset->LRAWinBinCntIdx >= m_UsedWinBinCnts)
	pReadset->LRAWinBinCntIdx = pReadset->StartWinBinCntsIdx;
pCnts = &m_pWinBinCnts[pReadset->LRAWinBinCntIdx];
if(pCnts->ReadsetID != pReadset->ReadsetID ||			// check cnts are for correct readset
	pCnts->ChromID > pWinBinCnts->ChromID ||			// WinBinCnts expected to have been sorted in ReadsetID.ChromID.Loci ascending order
	(pCnts->ChromID == pWinBinCnts->ChromID && pCnts->StartLoci > pWinBinCnts->StartLoci))		// so if LRAWinBinCntIdx is ahead of pWinBinCnts then do a linear search from the start of bins for the readset
	{
	pReadset->LRAWinBinCntIdx = pReadset->StartWinBinCntsIdx;
	pCnts = &m_pWinBinCnts[pReadset->StartWinBinCntsIdx];
	}
LRAWinBinCntIdx = pReadset->LRAWinBinCntIdx;
NumWinBins = pReadset->StartWinBinCntsIdx + pReadset->NumWinBins - LRAWinBinCntIdx;
for(WinBinCntIdx = 0;WinBinCntIdx < NumWinBins;WinBinCntIdx++,LRAWinBinCntIdx++,pCnts++)
	{
	if(pCnts->ReadsetID != pReadset->ReadsetID ||
		pCnts->ChromID > pWinBinCnts->ChromID ||
		(pCnts->ChromID == pWinBinCnts->ChromID && pCnts->StartLoci > pWinBinCnts->StartLoci))
			break;
	if(pCnts->ChromID == pWinBinCnts->ChromID && pCnts->StartLoci == pWinBinCnts->StartLoci)
		{
		pReadset->LRAWinBinCntIdx = LRAWinBinCntIdx;
		return(pCnts);
		}
	}
pReadset->LRAWinBinCntIdx = pReadset->StartWinBinCntsIdx;
return(NULL);	// unable to locate
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

int
CCallHaplotypes::LoadPBAFile(char* pszFile,bool bIsSkim)	// load chromosome metadata and PBA data from this file
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
if((ReadsetID = AddReadset(szReadsetID))==0)
	{
	gDiagnostics.DiagOut (eDLFatal, gszProcName, "Input file '%s' duplicates the ReadsetID '%s' of a previously loaded readset",pszFile,szReadsetID);
	return(eBSFerrOpnFile);
	}

pReadsetMetadata = &m_Readsets[ReadsetID-1];
pReadsetMetadata->bIsSkim = bIsSkim;
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
// if m_NChroms > 0 then limit number of chroms accepted to m_NChroms

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
	m_InNumProcessed += ChromLen;
	if((m_InNumBuffered - m_InNumProcessed) < 110 && FillInBuffer(ChromLen) == 0)
		break;
	if(m_NChroms > 0 && ChromID == m_NChroms)
		break;
	}

close(m_hInFile);
m_hInFile = -1;

return(eBSFSuccess);
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
CCallHaplotypes::LocatePBAfor(uint32_t ReadSetID,		// readset identifier 
			 uint32_t ChromID)			// chrom identifier
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



tsBinCnts *							// returned ptr to allocated and initialised windowed founder counts
CCallHaplotypes::AddWinBinCnts(tsBinCnts *pInitWinBinCnts)	// allocated tsWinFnderCnts to be initialised with WinBinCnts
{
uint32_t ToAllocdWinBinCnts;
tsBinCnts *pWinBinCnts;
size_t memreq;
AcquireSerialise();
if (m_pWinBinCnts == NULL)					// may be NULL first time in
	{
	memreq = cAllocWinBinCnts * sizeof(tsBinCnts);
#ifdef _WIN32
	m_pWinBinCnts = (tsBinCnts *)malloc((size_t)memreq);
	if (m_pWinBinCnts == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AllocWinBinCnts: Memory allocation of %lld bytes failed", (int64_t)memreq);
		return(NULL);
		}
#else
	m_pWinBinCnts = (tsBinCnts *)mmap(NULL, (size_t)memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	if (m_pWinBinCnts == MAP_FAILED)
		{
		m_pWinBinCnts = NULL;
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddWinBinCnts: Memory allocation of %lld bytes through mmap()  failed", (int64_t)memreq, strerror(errno));
		return(NULL);
		}
#endif
	m_AllocdWinBinCntsMem = memreq;
	m_AllocdWinBinCnts = cAllocWinBinCnts;
	m_UsedWinBinCnts = 0;
	}
else
		// needing to allocate more memory?
	if ((m_UsedWinBinCnts) >= m_AllocdWinBinCnts)
		{
		ToAllocdWinBinCnts = m_UsedWinBinCnts + cAllocWinBinCnts;
		size_t memreq = ToAllocdWinBinCnts * sizeof(tsBinCnts);
#ifdef _WIN32
		pWinBinCnts = (tsBinCnts*)realloc(m_pWinBinCnts, memreq);
		if (pWinBinCnts == NULL)
			{
#else
		pWinBinCnts = (tsBinCnts*)mremap(m_pWinBinCnts, m_AllocdWinBinCntsMem, memreq, MREMAP_MAYMOVE);
		if (pWinBinCnts == MAP_FAILED)
			{
#endif
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "AllocWinBinCnts: Memory reallocation to %lld bytes failed - %s", memreq, strerror(errno));
			return(NULL);
			}
		m_pWinBinCnts = pWinBinCnts;
		m_AllocdWinBinCntsMem = memreq;
		m_AllocdWinBinCnts=ToAllocdWinBinCnts;
		}
pWinBinCnts = &m_pWinBinCnts[m_UsedWinBinCnts++];
if(pInitWinBinCnts != NULL)
	*pWinBinCnts = *pInitWinBinCnts;
else
	memset(pWinBinCnts,0,sizeof(tsBinCnts));
ReleaseSerialise();
return(pWinBinCnts);
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
	Sleep(1000);
#else
	sleep(1);
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
	Rslt = ClassifyChromHaplotypes(pWorkQueueEl->ReadsetID,pWorkQueueEl->ChromID,pWorkQueueEl->StartLoci,pWorkQueueEl->EndLoci,pWorkQueueEl->ChromLen,pWorkQueueEl->AccumBinSize,pWorkQueueEl->NumFndrs,pWorkQueueEl->pPBAs);
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


int
CCallHaplotypes::CountHaplotypes(uint32_t ReadsetID,			// this is the readset being aligned to founders - could be a skim or founder
									uint32_t NumFndrs,			// number of founders to be processed against
									int KmerSize,				// aligning founders using this Kmer size
									uint32_t AccumBinSize)		// comparing ReadsetID's PBA against founder PBAs then counting using this bin size
{
uint32_t FounderID;
tsReadsetMetadata *pReadsetMetadata;
tsChromMetadata *pChromMetadata;
tsWorkQueueEl *pWorkQueueEl;
uint32_t CurChromMetadataIdx;
uint8_t *pPBAs[cMaxFounderReadsets+1];							// additional is to allow for the skim readset PBAs

char *pszReadset = LocateReadset(ReadsetID);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Process: PBA Aligning '%s' ...",pszReadset);

pReadsetMetadata = &m_Readsets[ReadsetID-1]; 
CurChromMetadataIdx = pReadsetMetadata->StartChromMetadataIdx;
m_ScoreKmerSize = KmerSize;
m_pWorkQueueEls = new tsWorkQueueEl [pReadsetMetadata->NumChroms * m_NumThreads * 2]; 
m_TotWorkQueueEls = 0;
m_NumQueueElsProcessed = 0;
pWorkQueueEl = m_pWorkQueueEls;
// iterate along PBAs for skim readset and align Kmers of length KmerSize to all founder PBAs at same loci
for(uint32_t ChromIdx = 0; ChromIdx < pReadsetMetadata->NumChroms && CurChromMetadataIdx != 0; ChromIdx++)
	{
	pChromMetadata = &m_pChromMetadata[CurChromMetadataIdx - 1];
	pPBAs[0] = pChromMetadata->pPBAs;

	// all founders must have PBAs for same chromosome as the skim
	for(FounderID = 1; FounderID <= NumFndrs; FounderID++)
		{
		if((pPBAs[FounderID] = LocatePBAfor(FounderID, pChromMetadata->ChromID)) == NULL)
			{
			char *pszChrom = LocateChrom(pChromMetadata->ChromID);
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "CountHaplotypes: No PBA for chromosome '%s' in founder '%s', skipping this chromosome", pszChrom, LocateReadset(FounderID));
			break;
			}
		}
	CurChromMetadataIdx = pChromMetadata->NxtChromMetadataIdx;
	if(FounderID <= NumFndrs)
		continue;

// divide chroms into units of work, such that maximal number of threads can be utilised
	uint32_t NumUnits;
	uint32_t UnitSize;
	// what is the maximal sized unit that a thread would be expected to process in this chromosome
	if(pChromMetadata->ChromLen < 1000000)	// too much overhead if too many threads processing shorter chromosomes
		{
		UnitSize = pChromMetadata->ChromLen;
		NumUnits = 1;
		}
	else
		{
		// unit size has to be a multiple of WindowSize, except for last unit which can be shorter
		NumUnits = (pChromMetadata->ChromLen + AccumBinSize - 1) / (AccumBinSize * m_NumThreads);
		if(NumUnits < 2)
			NumUnits = 2;
		UnitSize = NumUnits * AccumBinSize;
		}

	uint32_t StartLoci = 0;
	uint32_t EndLoci = UnitSize - 1;
	while(StartLoci < (pChromMetadata->ChromLen - KmerSize))
		{
		pWorkQueueEl->ReadsetID = ReadsetID;
		pWorkQueueEl->NumFndrs = NumFndrs;
		pWorkQueueEl->ChromID = pChromMetadata->ChromID;
		pWorkQueueEl->StartLoci = StartLoci;
		pWorkQueueEl->EndLoci = EndLoci;
		pWorkQueueEl->ChromLen = pChromMetadata->ChromLen;
		pWorkQueueEl->AccumBinSize = AccumBinSize;
		memcpy(pWorkQueueEl->pPBAs,pPBAs,(NumFndrs + 1) * sizeof(uint8_t *));
		m_TotWorkQueueEls++;
		pWorkQueueEl++;
		StartLoci += UnitSize;
		EndLoci = min(StartLoci + UnitSize - 1,pChromMetadata->ChromLen);
		}
	}

// startup threads
StartWorkerThreads(m_NumThreads,		// there are this many threads in pool
					pReadsetMetadata->NumChroms);		// processing this number of chromosomes

while(!WaitAlignments(60))
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Process: PBA Aligning '%s' ...",pszReadset);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Process: PBA Alignment for '%s' completed",pszReadset);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Process: Sorting bin counts ...");
m_mtqsort.SetMaxThreads(m_NumThreads);
m_mtqsort.qsort(m_pWinBinCnts,(int64_t)m_UsedWinBinCnts,sizeof(tsBinCnts),SortWinBinCnts);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Process: Sorting bin counts completed");
// ready for next founder or skim PBA processing
if(m_pWorkQueueEls != NULL)
	{
	delete []m_pWorkQueueEls;
	m_pWorkQueueEls = NULL;	
	m_NumQueueElsProcessed = 0;
	}

return(eBSFSuccess);
}


int
CCallHaplotypes::ClassifyChromHaplotypes(uint32_t ReadsetID,			// this readset is used for identifying the readset being aligned to founders - could be a skim or founder
										uint32_t ChromID,				// processing is for this chromosome
										uint32_t StartLoci,				// starting from this loci inclusive
										uint32_t EndLoci,				// and ending at this loci inclusive
										uint32_t ChromLen,				// chromosome length
										uint32_t AccumBinSize,		// accumulating counts into this sized non-overlapping bins
										uint32_t NumFndrs,				// number of founders to be processed against the skim (or founder) PBA at *pPBAs[0]
										uint8_t *pPBAs[])				// pPBAs[0] pts to chromosome skim PBA, followed by ptrs to each of the chromosome founder PBAs 
{
tsBinCnts WinBinCnts;
uint32_t Loci;
uint32_t NumWinBinCnts = 0;
int64_t FndrsScoreClassifications;

memset(&WinBinCnts,0,sizeof(WinBinCnts));
WinBinCnts.ReadsetID = ReadsetID;
WinBinCnts.ChromID = ChromID;
WinBinCnts.StartLoci = StartLoci; 
for(Loci = StartLoci; Loci <=  EndLoci && Loci <= (ChromLen - m_ScoreKmerSize); Loci++)
	{
	if((Loci % AccumBinSize) == 0)	// starting a new window to integrate counts over
		{
		if(Loci != StartLoci)			// reporting counts for previous window
			{
			WinBinCnts.EndLoci = Loci-1;
			AddWinBinCnts(&WinBinCnts);
			NumWinBinCnts++;
			memset(&WinBinCnts,0,sizeof(WinBinCnts));
			WinBinCnts.ReadsetID = ReadsetID;
			WinBinCnts.ChromID = ChromID;
			WinBinCnts.StartLoci = Loci;
			}
		}

	if((FndrsScoreClassifications = ScoreSkimLoci(ReadsetID, Loci, ChromLen, m_ScoreKmerSize, NumFndrs, pPBAs)) < 0)
		return((int)FndrsScoreClassifications);

	// classify founders from alignment scores
	// 2bits per founder, founder 1 in bits 0 to 1, through to founder 16 in bits 30..31
	// for each founder -
	//	if highest scorer then will be classified as 3
	//	if equal highest then will be classified as 2
	//	if not highest or equal highest, but still accepted as aligned, then will be classified as 1
	//	if not accepted as aligned then will be classified as 0
	if(FndrsScoreClassifications == 0)	// only 0 if no PBA alignments were accepted for skim against any founder
		WinBinCnts.NoFounder++;
	else
		{
		uint32_t FounderID;
		bool bUnique = false;
		for(FounderID = 1; FounderID <= NumFndrs && FndrsScoreClassifications != 0; FounderID++, FndrsScoreClassifications >>= 2)
			{
			switch(FndrsScoreClassifications & 0x03)
				{
				case 0x03:		// uniquely aligned
					WinBinCnts.Uniques[FounderID-1]++;
					bUnique = true;
					break;
				case 0x02:		// non-uniquely aligned - currently no interest if this was next best aligned
				case 0x01:
					WinBinCnts.NonUniques[FounderID-1]++;
					break;
				case 0x00:		// no alignments to this founder
					break;
				}
			}
		if(!bUnique)	// if none exclusively aligned then must have been alignments to multiple founders
			WinBinCnts.MultiFounder++;
		}
	}
// report here if partial AccumBinSize 
if(WinBinCnts.StartLoci < Loci-1)
	{
	WinBinCnts.EndLoci = Loci-1;
	AddWinBinCnts(&WinBinCnts);
	NumWinBinCnts++;
	}
return(eBSFSuccess);
}


int32_t				// returned Kmer score resulting from PBA alignment of pFndrA against pFndrB
CCallHaplotypes::ScoreFounderLoci(uint32_t Loci,// processing for Kmer haplotypes starting from this loci
						uint32_t SeqLen,		// sequence or chromosome is this length
						uint32_t KmerSize,		// Kmer to score on is of this size
						uint8_t *pAFndrPBA,		// PBA for FndrA
						uint8_t *pBFndrPBA)		// PBA for FndrB
{
uint8_t *pFndrALoci;
uint8_t *pFndrBLoci;
uint32_t AlleleIdx;
uint8_t AlleleMsk;
uint32_t KmerOfs;
uint32_t KmerLoci;

// quick check for skim and all founders having called PBA at specified loci
if(pAFndrPBA == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ScoreFounderLoci: PBA NULL pFndrA");
	return(-1);
	}
if(pBFndrPBA == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ScoreFounderLoci: PBA NULL pFndrB");
	return(-1);
	}

if((Loci + KmerSize) > SeqLen)
	return(0);

if(*(pFndrALoci = pAFndrPBA+Loci)==0 || *(pFndrBLoci = pBFndrPBA+Loci)==0)
	return(0);

uint32_t FounderScore;
FounderScore = 0;
KmerOfs = 0;
for(KmerLoci = Loci, KmerOfs = 0; KmerLoci <= (SeqLen - KmerSize) && KmerOfs < KmerSize; KmerOfs++,KmerLoci++,pFndrALoci++,pFndrBLoci++)
	{
	if(*pFndrALoci == 0 || *pFndrBLoci == 0)
		return(0);

	AlleleMsk = 0x03;	
	uint32_t FndrABaseAllele;
	uint32_t FndrBBaseAllele;
	uint32_t MaxAlleleScore = 0;
	uint32_t AlleleScore = 0;
	for(AlleleIdx = 0; AlleleIdx < 4; AlleleIdx++, AlleleMsk <<= 2)
		{
		FndrABaseAllele = (*pFndrALoci & AlleleMsk) >> (AlleleIdx * 2);
		FndrBBaseAllele = (*pFndrBLoci & AlleleMsk) >> (AlleleIdx * 2);
		AlleleScore = FndrABaseAllele * FndrBBaseAllele;
		if(AlleleScore)
			MaxAlleleScore = 10;
		}
	FounderScore += MaxAlleleScore;
	}
return(FounderScore);
}

int64_t				// < 0 if errors, returned packed Kmer scores, (2 bits per founder, founder 1 in bits 0..1, founder 31 in bits 60..61) scored as being as being present at this skim loci
CCallHaplotypes::ScoreSkimLoci(uint32_t ReadsetID,	// this readset is used for identifying the readset being aligned to founders - could be a skim or founder
						uint32_t Loci,			// processing for Kmer haplotypes starting from this loci
						uint32_t SeqLen,		// sequence or chromosome is this length
						uint32_t KmerSize,		// Kmer to score on is of this size
						uint32_t NumFndrs,		// number of founders to be processed against the skim (or founder) PBA at *pPBAs[0]
						uint8_t *pPBA[])		// pts to loci 0 of PBAs in readset order
{
uint32_t FounderID;
uint32_t NumFndrScoringAlleles;	// keeping count of number of scoring alleles for this founder so can determine if this founders proportion of scoring alleles reaches the acceptance threshold
int64_t InferencedFounders;
uint8_t *pSkimLoci;
uint8_t *pFndrLoci;
uint32_t AlleleIdx;
uint8_t AlleleMsk;
uint32_t KmerOfs;
uint32_t KmerLoci;

// quick check for skim and all founders having called PBA at specified loci
if( pPBA[0] == NULL)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"ScoreSkimLoci: PBA NULL skim ...");
	return(-1);
	}
pSkimLoci = pPBA[0] + Loci;
if(*pSkimLoci == 0)	// will be 0 if no skim allele aligned to starting Kmer Loci
	return(0);		// can't score ...
// all founders must have a PBA and alleles at initial Loci
for(FounderID = 1; FounderID <= NumFndrs; FounderID++)
	{
	if((pFndrLoci=pPBA[FounderID]) == NULL)
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"ScoreSkimLoci: PBA NULL founder at loci %u...",Loci);
		return(-1);
		}
	if(*(pFndrLoci+Loci)==0)
		return(0);
	}

// use KmerSize as being a minimum, enabling extension out to a maximum of 10x KmerSize
uint32_t ExtdKmerSize = KmerSize * 10;
for(FounderID = 0; FounderID <= NumFndrs; FounderID++) // including extension of skim as if a founder when determining maximal Kmer
	{
	pFndrLoci = pPBA[FounderID] + Loci;
	KmerOfs = 0;
	for(KmerLoci = Loci, KmerOfs = 0; KmerLoci <= (SeqLen - KmerSize) && KmerOfs < ExtdKmerSize; KmerOfs++,KmerLoci++,pFndrLoci++)
		if(*pFndrLoci == 0)	// will be 0 if no skim/fndr alleles aligned to this Loci within the Kmer
			break;		
	if(KmerOfs < ExtdKmerSize)		// reducing down to minimum extension of any founder or skim
		ExtdKmerSize = KmerOfs;
	}
if(ExtdKmerSize < KmerSize)			// extension has to be at least original minimal KmerSize
	return(0);
KmerSize = ExtdKmerSize; 

int32_t FounderScores[cMaxFounderReadsets];
int32_t *pFounderScore;
memset(FounderScores,0,sizeof(FounderScores));

InferencedFounders = 0;

pFounderScore = FounderScores;
for(FounderID = 1; FounderID <= NumFndrs; FounderID++,pFounderScore++)
	{
	NumFndrScoringAlleles = 0;		
	pFndrLoci = pPBA[FounderID] + Loci;
	pSkimLoci = pPBA[0] + Loci;
	KmerOfs = 0;
	for(KmerLoci = Loci, KmerOfs = 0; KmerLoci <= (SeqLen - KmerSize) && KmerOfs < KmerSize; KmerOfs++,KmerLoci++,pSkimLoci++,pFndrLoci++)
		{
		if(*pSkimLoci == 0 || *pFndrLoci == 0)	// will be 0 if no skim or founder alleles aligned to this Loci within the KmerSize
			break;		

		AlleleMsk = 0x03;	
		int32_t FndrBaseAllele;
		int32_t SkimBaseAllele;
		int32_t AlleleScore = 0;
		for(AlleleIdx = 0; AlleleIdx < 4; AlleleIdx++, AlleleMsk <<= 2)
			{
			SkimBaseAllele =  (*pSkimLoci & AlleleMsk) >> (AlleleIdx * 2);
			if(SkimBaseAllele < 2)
				continue;
			FndrBaseAllele = (*pFndrLoci & AlleleMsk) >> (AlleleIdx * 2);
			if(FndrBaseAllele < 2)
				continue;
			AlleleScore = 6;					// fixed score irrespective of the confidence in the allele match
			break;
			}
		if(AlleleScore)
			NumFndrScoringAlleles++;
		else
			AlleleScore = -18;				// penalise (3x a match) if no allele match at current loci within the Kmer
		*pFounderScore += AlleleScore;
		}
	// did this founder reach minimum proportion of required allele matches for score acceptance?
	if((int32_t)((NumFndrScoringAlleles * 100) / KmerSize) < m_AcceptKmerMatchPerc)
		*pFounderScore = 0;
	}

// classify from alignment scores
// if highest scorer then classify as 3
// if equal highest then classify as 2
// if not highest or equal highest, and accepted as aligned, then classify as 1
// if not accepted as aligned then classify as 0
uint64_t FndrHighest;
uint64_t FndrNxtHighest;
uint64_t FndrAccepted;
int32_t FndrHighestScore = 0;
int32_t FndrNxtHighestScore = 0;

// locate highest and second highest founder scores
pFounderScore = FounderScores;
for(FounderID = 1; FounderID <= NumFndrs; FounderID++,pFounderScore++)
	{
	if(*pFounderScore <= 0)			
		continue;					

	if(*pFounderScore > FndrHighestScore)
		{
		FndrNxtHighestScore = FndrHighestScore;		// next highest becomes the previous highest
		FndrHighestScore = *pFounderScore;
		}
	else
		if(*pFounderScore <= FndrHighestScore && *pFounderScore > FndrNxtHighestScore)
			FndrNxtHighestScore = *pFounderScore;
	}

FndrHighest = 0x03;
FndrNxtHighest = 0x02;
FndrAccepted = 0x01;
pFounderScore = FounderScores;
for(FounderID = 1; FounderID <= NumFndrs;FounderID++,pFounderScore++,FndrAccepted <<= 2,FndrHighest<<=2,FndrNxtHighest<<=2)
	{
	if(*pFounderScore == 0)
		continue;
	if(FndrHighestScore != 0 && *pFounderScore >= FndrHighestScore && FndrHighestScore > FndrNxtHighestScore)
		InferencedFounders |= FndrHighest;
	else
		if(*pFounderScore >= FndrNxtHighestScore)
			InferencedFounders |= FndrNxtHighest;
		else
			InferencedFounders |= FndrAccepted;
	}

return(InferencedFounders);
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

	// each returned packed haplotype occupies 4 bits, allowing for up to a total of 16 haplotype states
	//			0 - unable to call as either there are no founder and/or skim counts
	//			1 - Fa as likely the only haplotype
	//			2 - Fb as likely the only haplotype
	//			3 - Fa + Fb as likely both haplotypes present
uint64_t						
CCallHaplotypes::ChooseHaplotype(uint64_t HMMBinHaplotypes,		// 16 packed previously called haplotypes for use in as HMM transition probabilities, bits 0..3 contain previously called, through to bits 60..63 as oldest called (32 state)
						tsBinCnts* pSkimWinBinCnts)	// bin counts for skim
{
uint32_t SFaCnts;	// observed counts for Skim exclusively aligning to Fa when Skim aligned to Fa+Fb
uint32_t SFbCnts;	// observed counts for Skim exclusively aligning to Fb when Skim aligned to Fa+Fb
uint32_t SFaFbCnts;	// observed counts for Skim multifounders when Skim aligned to Fa+Fb

// because of noise there may have been an inconsistent haplotype call made earlier
// treat as inconsistent a call which is bracketed by calls which are same but different to the potentially incorrect call
// inconsistent calls are set to be same as bracketing calls
uint32_t HapIdx;
uint64_t Hap1Msk = 0x000f;
uint64_t CentralHapMsk = 0x00f0;
uint64_t Hap2Msk = 0x0f00;
uint64_t Haplotypes = HMMBinHaplotypes;
for(HapIdx = 1; HapIdx < 15; HapIdx++, Hap1Msk <<= 4,CentralHapMsk <<= 4, Hap2Msk <<= 4)
	{
	if((Haplotypes & Hap2Msk) == ((Haplotypes & Hap1Msk) << 8))
		{
		if((Haplotypes & CentralHapMsk) != ((Haplotypes & Hap1Msk) << 4))
			{
			HMMBinHaplotypes &= ~CentralHapMsk;
			HMMBinHaplotypes |= ((Haplotypes & Hap1Msk) << 4);
			}
		}
	}

SFaCnts = pSkimWinBinCnts->Uniques[0];
SFbCnts = pSkimWinBinCnts->Uniques[1];
SFaFbCnts = pSkimWinBinCnts->MultiFounder;

if((SFaCnts + SFbCnts + SFaFbCnts) < 5)	// can't call haplotype if there are insufficient aligning skim counts
	return(HMMBinHaplotypes << 4);		// if all counts are near 0 then must assume that there is a deletion in all haplotypes at the windowed bin loci
if((SFaCnts + SFbCnts) == 0)	// if no uniques and all multialigners then use previously called haplotype 
	return((HMMBinHaplotypes << 4) | HMMBinHaplotypes & 0x0f);
double FaFbProp = (double)(SFaCnts + SFbCnts)/(double)(SFaCnts + SFbCnts + SFaFbCnts);
double FaProp = (double)SFaCnts/ (SFaCnts + SFbCnts);
double FbProp = (double)SFbCnts/ (SFaCnts + SFbCnts);
double FaPropThres;
double FbPropThres;

// thresholds are dependent on what the previous call was ...
switch(HMMBinHaplotypes & 0x0f) {
	case 0:		// no call, or was a deletion, biasing will be for both haplotypes being present
		FaPropThres = 0.950;
		FbPropThres = 0.950;
		break;

	case 1:		// Fa was previously called so bias thresholds towards Fa and against Fb
		FaPropThres = 0.800;
		FbPropThres = 0.900;
		break;

	case 2:		// Fb was previously called so bias thresholds towards Fb and against Fa
		FaPropThres = 0.900;
		FbPropThres = 0.800;
		break;

	case 3:		// Fa+Fb previously called so bias thresholds towards Fa+Fb
		FaPropThres = 0.970;
		FbPropThres = 0.970;
		break;
	}

if(SFaCnts == 0)		// if only uniques for a single haplotype then call that as being the haplotype
	return((HMMBinHaplotypes << 4) | 2);

if(SFbCnts == 0)
	return((HMMBinHaplotypes << 4) | 1);

// have a mixture of haplotypes
if(FaProp >= FaPropThres)
	return((HMMBinHaplotypes << 4) | 1);
if(FbProp >= FbPropThres)
	return((HMMBinHaplotypes << 4) | 2);
if(SFaCnts >= 5 || SFbCnts >= 5)
	return((HMMBinHaplotypes << 4) | 3);
return((HMMBinHaplotypes << 4) | HMMBinHaplotypes & 0x0f);	// too close to call, use previously called haplotype
}



uint32_t		// returned chrom identifier, 0 if unable to accept this chromosome name
CCallHaplotypes::AddChrom(char* pszChrom) // associate unique identifier with this chromosome name
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
CCallHaplotypes::LocateChrom(char* pszChrom) // return unique identifier associated with this chromosome name
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
CCallHaplotypes::LocateChrom(uint32_t ChromID)
{
if(ChromID < 1 || ChromID > m_NumChromNames)
	return(NULL);
return(&m_szChromNames[m_szChromIdx[ChromID-1]]);
}


// NOTE: Readsets are checked for uniqueness as readsets must be unique
uint32_t		// returned readset identifier, 0 if unable to accept this readset name
CCallHaplotypes::AddReadset(char* pszReadset) // associate unique identifier with this readset name
{
uint32_t ReadsetNameIdx;
int ReadsetNameLen;
char *pszLAname;

// with any luck the sequence name will be same as the last accessed
if((pszLAname = LocateReadset(m_LAReadsetNameID)) != NULL)
	if(!stricmp(pszReadset,pszLAname))
		return(0);
		
// iterate over all known readsets in case this readset to add is a duplicate
for(ReadsetNameIdx = 0; ReadsetNameIdx < m_NumReadsetNames; ReadsetNameIdx++)
	if(!stricmp(pszReadset, &m_szReadsetNames[m_szReadsetIdx[ReadsetNameIdx]]))
		{
		m_LAReadsetNameID = ReadsetNameIdx + 1;
		return(0);
		}

// not a duplicate
ReadsetNameLen = (int)strlen(pszReadset);
if((m_NxtszReadsetIdx + ReadsetNameLen + 1) > (int)sizeof(m_szReadsetNames))
	return(0);
if(m_NumReadsetNames == (cMaxFounderReadsets + 1))
	return(0);

m_szReadsetIdx[m_NumReadsetNames++] = m_NxtszReadsetIdx;
strcpy(&m_szReadsetNames[m_NxtszReadsetIdx], pszReadset);
m_NxtszReadsetIdx += ReadsetNameLen + 1;
m_LAReadsetNameID = m_NumReadsetNames;
return(m_LAReadsetNameID);
}


uint32_t		// returned Readset identifier, 0 if unable to locate this Readset name
CCallHaplotypes::LocateReadset(char* pszReadset) // return unique identifier associated with this Readset name
{
uint32_t ReadsetNameIdx;
char *pszLAReadset;

// with any luck the Readset name will be same as the last accessed
if((pszLAReadset = LocateReadset(m_LAReadsetNameID)) != NULL)
	if(!stricmp(pszReadset,pszLAReadset))
		return(m_LAReadsetNameID);

// iterate over all known Readsets
for(ReadsetNameIdx = 0; ReadsetNameIdx < m_NumReadsetNames; ReadsetNameIdx++)
	if(!stricmp(pszReadset, &m_szReadsetNames[m_szReadsetIdx[ReadsetNameIdx]]))
		{
		m_LAReadsetNameID = ReadsetNameIdx + 1;
		return(m_LAReadsetNameID);
		}
return(0);
}

char* 
CCallHaplotypes::LocateReadset(uint32_t ReadsetID)
{
if(ReadsetID < 1 || ReadsetID > m_NumReadsetNames)
	return(NULL);
return(&m_szReadsetNames[m_szReadsetIdx[ReadsetID-1]]);
}

// sorting by ReadsetID.ChromID.StartLoci ascending
int
CCallHaplotypes::SortWinBinCnts(const void* arg1, const void* arg2)
{
tsBinCnts* pEl1 = (tsBinCnts*)arg1;
tsBinCnts* pEl2 = (tsBinCnts*)arg2;
if(pEl1->ReadsetID < pEl2->ReadsetID)
	return(-1);
if(pEl1->ReadsetID > pEl2->ReadsetID)
	return(1);
if(pEl1->ChromID < pEl2->ChromID)
	return(-1);
if(pEl1->ChromID > pEl2->ChromID)
	return(1);
if(pEl1->StartLoci < pEl2->StartLoci)
	return(-1);
if(pEl1->StartLoci > pEl2->StartLoci)
	return(1);
return(0);
}

