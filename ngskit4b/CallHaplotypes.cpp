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
			 int MinKmerSize,		// aligning founders using this minimum Kmer size
			 int MaxKmerSize,				// aligning founders using this maximum Kmer size
			int32_t MinAlignedScoringPBA,		// to be accepted as aligned and scored then PBAs in both founder and skim must be at least this PBA level
			int32_t MinAlignedPBA,		// to be accepted as aligned ( < m_MinScoringPBA) with no penalty applied then PBAs in both founder and skim must be at least this PBA level
			int32_t AlignedPBAScore,		// if PBA accepted for scoring (>= m_MinAlignedScoringPBA) then add this score to the summed Kmer alignment score
			int32_t PenaltyPBAScore,		// if PBA not accepted as aligned ( < m_MinAlignedPBA) then apply this penalty to the summed Kmer alignment score
			int32_t AcceptKmerMatchPerc,	// only accept Kmer alignment if at least this percentage of base alleles matched at >= m_MinAlignedPBA
			int32_t MinKmerScoreDiff,		// must be at least this Kmer percentage score differential between highest scoring founders before accepting a founder as being a uniquely highest scoring
			int AccumBinSize,		// counts are accumulated into this sized (bp) non-overlapping bins
			int WinMeanBins,		// derive haplotype calls from mean of raw haplotypes in a sliding window containing this many bins
			int MaxFndrHaps,		// process for reduction down to this maximum number of called haplotypes
			int Ploidy,				// smoothed haplotypes for organism having this ploidy
			char *pszMaskBPAFile,	// optional masking input BPA file, only process BPAs which are intersect of these BPAs and skim plus founder BPAs
			int NumFounderInputFiles,	// number of input founder file specs
			char *pszFounderInputFiles[],	// names of input founder PBA files (wildcards allowed)
			int NumSkimInputFiles,	// number of input skim file specs
			char* pszSkimInputFiles[],		// names of input skim PBA files (wildcards allowed)
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
	 int MinKmerSize;				// aligning founders using this minimum Kmer size
	 int MaxKmerSize;				// aligning founders using this maximum Kmer size
	int32_t MinAlignedScoringPBA;	// to be accepted as aligned and scored then PBAs in both founder and skim must be at least this PBA level
	int32_t MinAlignedPBA;			// to be accepted as aligned ( < m_MinScoringPBA) with no penalty applied then PBAs in both founder and skim must be at least this PBA level
	int32_t AlignedPBAScore;		// if PBA accepted for scoring (>= m_MinAlignedScoringPBA) then add this score to the summed Kmer alignment score
	int32_t PenaltyPBAScore;		// if PBA not accepted as aligned ( < m_MinAlignedPBA) then apply this penalty to the summed Kmer alignment score
	int32_t AcceptKmerMatchPerc;	// only accept Kmer alignment if at least this percentage of base alleles matched at >= m_MinAlignedPBA
	int32_t MinKmerScoreDiff;		// must be at least this Kmer percentage score differential between highest scoring founders before accepting a founder as being a uniquely highest scoring
	int AccumBinSize;				// counts are accumulated into this sized (Kbp) non-overlapping bins
	int WinMeanBins;				// derive haplotype calls from mean of counts in a sliding window containing this many bins
	int MaxFndrHaps;				// process for reduction down to this maximum number of called founder haplotypes
	int Ploidy;						// smoothed haplotypes for organism having this ploidy
	 char szTrackName[_MAX_PATH];	// track name
	 char szTrackDescr[_MAX_PATH];	// track description
	 int NumFounderInputFiles;		// number of input founder BPA files
	 char *pszFounderInputFiles[cMaxFounderReadsets];		// names of input founder BPA files (wildcards allowed)
	 int NumSkimInputFiles;		// number of input skim BPA files
	 char *pszSkimInputFiles[cMaxSkimReadsets];		// names of input skim BPA files (wildcards allowed)

	 char szOutFile[_MAX_PATH];		// Windowed haplotype calls output file base name, skim readset identifier is appended to this base name
	 char szMaskBPAFile[_MAX_PATH];	// optional masking input BPA file, only process BPAs which are intersect of these BPAs and skim plus founder BPAs

	struct arg_lit *help = arg_lit0 ("h", "help", "print this help and exit");
	struct arg_lit *version = arg_lit0 ("v", "version,ver", "print version information and exit");
	struct arg_int *FileLogLevel = arg_int0 ("f", "FileLogLevel", "<int>", "Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
	struct arg_file *LogFile = arg_file0 ("F", "log", "<file>", "diagnostics log file");

	struct arg_int *pmode = arg_int0 ("m", "mode", "<int>", "processing mode: 0 generate haplotype bin counts over all chromosomes using skim BPAs compared with founder BPAs, 1: generate IGV for each input PBA");
	struct arg_int *nchroms = arg_int0 ("n", "nchroms", "<int>","Limit processing: 0 All chromosomes, 1..n 1st to nTh chromosome");
	struct arg_int *minkmersize = arg_int0 ("k", "minkmersize", "<int>", "allele alignments using this minimum sized Kmer (defaults to 50bp)");
	struct arg_int *maxkmersize = arg_int0 ("K", "maxkmersize", "<int>", "allele alignments using this maximum sized Kmer (defaults to 10x minimum Kmer)");

	struct arg_int *kmermatchperc = arg_int0 ("x", "kmermatchperc", "<int>", "only accept as extended Kmer alignment if at least this percentage of base alleles match over that Kmer(default 97)");
	struct arg_int *minalignedscoringpba = arg_int0 ("a", "minalignedscoringpba", "<int>", "only aligned individual base PBAs scored if both founder and skim are at least this PBA level (default 2)");
	struct arg_int *minalignedpba = arg_int0 ("A", "minalignedpba", "<int>", "individual base PBAs accepted as aligned if PBAs in both founder and skim are at least this PBA level (default 1)");
	struct arg_int *alignedpbascore = arg_int0 ("s", "alignedpbascore", "<int>", "scored base PBAs alignment score (default 10)");
	struct arg_int *penaltypbascore = arg_int0 ("S", "penaltypbascore", "<int>", "unaligned base PBAs penalty score (default -30)");
	struct arg_int *minkmerscorediff = arg_int0 ("d", "minkmerscorediff", "<int>", "accumulated score differential between founder Kmers must be at least this percentage to accept as uniquely aligning (default 2)");

	struct arg_int *accumbinsize = arg_int0 ("b", "binsize", "<int>", "Kmer alignment counts accumulated into non-overlapping bins of this Kbp size (default 100Kbp)");
	struct arg_int *winmeanbins = arg_int0 ("B", "binsmean", "<int>", "derive smoothed haplotype calls from raw haplotype calls in sliding window containing this many bins (defaults to 20)");
	struct arg_int *maxfndrhaps = arg_int0 ("p", "maxfndrhaps", "<int>","process for reduction down to this maximum number of called founder haplotypes in any bin (default 2)");
	struct arg_int *ploidy = arg_int0 ("P", "ploidy", "<int>","smoothed haplotypes for organism having this ploidy (default 2)");
	struct arg_str *trackname = arg_str0("t","trackname","<str>","BED Track name");
	struct arg_str *trackdescr = arg_str0("d","trackdescr","<str>","BED Track description");
	struct arg_file *founderfiles = arg_filen("I", "founderfiles", "<file>", 0,cMaxFounderReadsets,"founder input BPA file(s), wildcards allowed, limit of 500 founders supported");
	struct arg_file *skimfiles = arg_filen("i", "inskimfile", "<file>",0, cMaxSkimReadsets, "skim input BPA file(s), wildcards allowed, limit of 1000 skim filespecs supported");
	struct arg_file *maskbpafile = arg_file0("c", "inmaskbpa", "<file>", "optional masking input BPA file, only process BPAs which are intersect of these BPAs and skim plus founder BPAs");
	struct arg_file *outfile = arg_file1("o", "out", "<file>", "Windowed haplotype calls output prefix (outputs CSV, BED and WIG format)");
	struct arg_int *threads = arg_int0("T","threads","<int>","number of processing threads 0..64 (defaults to 0 which limits threads to maximum of 64 CPU cores)");
	struct arg_end *end = arg_end (200);

	void *argtable[] = { help,version,FileLogLevel,LogFile,
						pmode,maxfndrhaps,ploidy,nchroms,minkmersize,maxkmersize,kmermatchperc,minalignedscoringpba,minalignedpba,alignedpbascore,penaltypbascore,minkmerscorediff,
						accumbinsize,winmeanbins,trackname,trackdescr,maskbpafile,skimfiles,founderfiles, outfile,threads,end };

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
		if(PMode < eMSHDefault || PMode >= eMSHPlaceHolder)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Unsupported processing mode '-m%d'\n",PMode);
			exit(1);
			}
		NChroms = nchroms->count ? nchroms->ival[0] : 0;
		if(NChroms < 0)
			NChroms = 0;
		else
			if(NChroms > 1000)
				NChroms = 1000;

		if(PMode == eMSHDefault)
			{
			MinKmerSize = minkmersize->count ? minkmersize->ival[0] : cDfltMinKmerSize;
			if(MinKmerSize < cMinKmerSize || MinKmerSize > cMaxKmerSize)
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: minimum Kmer size '-k%d' specified outside of range %d .. %dbp\n", MinKmerSize,cMinKmerSize,cMaxKmerSize);
				exit(1);
				}


			MaxKmerSize = maxkmersize->count ? maxkmersize->ival[0] : min(MinKmerSize*10,cMaxKmerSize);
			if(MaxKmerSize < MinKmerSize || MaxKmerSize > cMaxKmerSize)
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: maximum Kmer size '-K%d' specified outside of range %d .. %dbp\n", MaxKmerSize,MinKmerSize,cMaxKmerSize);
				exit(1);
				}

			AcceptKmerMatchPerc = kmermatchperc->count ? kmermatchperc->ival[0] : cDfltAcceptKmerMatchPerc;
			if(AcceptKmerMatchPerc < 90 || AcceptKmerMatchPerc > 100)
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Minimum Kmer allele match percentage '-K%d' specified outside of range 90.100\n", AcceptKmerMatchPerc);
				exit(1);
				}

			MinAlignedScoringPBA = minalignedscoringpba->count ? minalignedscoringpba->ival[0] : cDfltMinAlignedScoringPBA;
			if(MinAlignedScoringPBA < 1)
				MinAlignedScoringPBA = 1;
			else
				if(MinAlignedScoringPBA > 3)
					MinAlignedScoringPBA = 3;

			MinAlignedPBA = minalignedpba->count ? minalignedpba->ival[0] : cDfltMinAlignedPBA;
			if(MinAlignedPBA < 1)
				MinAlignedPBA = 1;
			else
				if(MinAlignedPBA > MinAlignedScoringPBA)
					MinAlignedPBA = MinAlignedScoringPBA;

			AlignedPBAScore = alignedpbascore->count ? alignedpbascore->ival[0] : cDfltAlignedPBAScore;
			if(AlignedPBAScore < 1)
				AlignedPBAScore = 1;
			else
				if(AlignedPBAScore > 100)
					AlignedPBAScore = 100;

			PenaltyPBAScore = penaltypbascore->count ? penaltypbascore->ival[0] : cDfltPenaltyPBAScore;
			if(PenaltyPBAScore > 0)
				PenaltyPBAScore *= -1;	
			if(PenaltyPBAScore < -1000)
				PenaltyPBAScore = -1000;

			MinKmerScoreDiff = minkmerscorediff->count ? minkmerscorediff->ival[0] : cDfltMinKmerScoreDiff;
			if(MinKmerScoreDiff > 50)
				MinKmerScoreDiff = 50;
			else
				if(MinKmerScoreDiff < 1)
					MinKmerScoreDiff = 1;

			AccumBinSize = accumbinsize->count ? accumbinsize->ival[0] : cDfltAccumBinSize;
			if(AccumBinSize < 1 || AccumBinSize > 2000)
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Accumulating non-overlapping bin size '-b%dKbp' specified outside of range 1..2000 Kbp\n", AccumBinSize);
				exit(1);
				}

			WinMeanBins = winmeanbins->count ? winmeanbins->ival[0] : cDfltNumMeanBins;
			if(WinMeanBins < 1 || WinMeanBins >= 100)
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Mean of counts in '-B%d' bins specified outside of range 1..100\n", WinMeanBins);
				exit(1);
				}

			MaxFndrHaps = maxfndrhaps->count ? maxfndrhaps->ival[0] : 2;
			if(MaxFndrHaps < 1 || MaxFndrHaps > 5)
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Process for reduction down to '-p%d' founder haplotypes specified outside of range 1..5\n", MaxFndrHaps);
				exit(1);
				}

			Ploidy = ploidy->count ? ploidy->ival[0] : 2;
			if(Ploidy < 1 || Ploidy > MaxFndrHaps)
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Process for organism poidy  '-P%d' must not be more than maximum haplotypes 1..%d\n", MaxFndrHaps);
				exit(1);
				}
			}
		else
			{
			MinKmerSize = cDfltMinKmerSize;
			MaxKmerSize = min(cDfltMinKmerSize * 10,cMaxKmerSize);
			AcceptKmerMatchPerc = cDfltAcceptKmerMatchPerc;
			MinAlignedScoringPBA = cDfltMinAlignedScoringPBA;
			MinAlignedPBA = cDfltMinAlignedPBA;
			AlignedPBAScore = cDfltAlignedPBAScore;
			PenaltyPBAScore = cDfltPenaltyPBAScore;
			MinKmerScoreDiff = cDfltMinKmerScoreDiff;
			AccumBinSize = cDfltAccumBinSize;
			WinMeanBins = cDfltNumMeanBins;
			MaxFndrHaps = 2;
			Ploidy = 2;
			}

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

	NumSkimInputFiles = 0;
	if(PMode == eMSHDefault)
		{
		if(skimfiles->count)
			{
			for(Idx = 0; NumSkimInputFiles < cMaxSkimReadsets && Idx < skimfiles->count; Idx++)
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
			}
		}

	if(PMode == eMSHDefault && !NumSkimInputFiles)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: After removal of whitespace, there are no input skim file(s) specified with '-i<filespec>' option)\n");
		exit(1);
		}

	NumFounderInputFiles = 0;
	if(founderfiles->count)
		{
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
		}

	if(!NumFounderInputFiles)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: After removal of whitespace, there are no input founder file(s) specified with '-I<filespec>' option)\n");
		exit(1);
		}

	if(PMode == eMSHDefault && maskbpafile->count)
		{
		strcpy (szMaskBPAFile, maskbpafile->filename[0]);
		CUtility::TrimQuotedWhitespcExtd (szMaskBPAFile);
		if(szMaskBPAFile[0] == '\0')
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: After removal of whitespace, no input masking BPA file specified with '-c<filespec>' option)\n");
			exit(1);
			}
		}
	else
		szMaskBPAFile[0] = '\0';

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
			pszDescr = "Calling haplotypes in skim PBAs through alignments with founder PBAs";
			break;
		case eMSHAlleles:
			pszDescr = "Reporting PBA allele detail in founder PBAs";
			break;
		}

	gDiagnostics.DiagOutMsgOnly (eDLInfo, "Segment haplotypes : '%s'", pszDescr);
	if(NChroms > 0)
		gDiagnostics.DiagOutMsgOnly (eDLInfo, "Debugging - only processing 1st %d chromosomes",NChroms);

	if(PMode == eMSHDefault)
		{
		gDiagnostics.DiagOutMsgOnly (eDLInfo, "Minimum Kmer size : %dbp", MinKmerSize);
				gDiagnostics.DiagOutMsgOnly (eDLInfo, "Maximum Kmer size : %dbp", MaxKmerSize);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Only accepting as aligned if at least this percentage of base alleles over Kmer match: %d",AcceptKmerMatchPerc);

		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Bases scored if PBAs in both founder and skim are at least this PBA level: %d",MinAlignedScoringPBA);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Bases aligned if PBAs in both founder and skim are be at least this PBA level: %d",MinAlignedPBA);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Scoring base PBAs accumulate this score to the summed Kmer alignment score: %d",AlignedPBAScore);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Unaligned base PBAs have this penalty accumulated to the summed Kmer alignment score: %d",PenaltyPBAScore);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Kmer percentage score differential between highest scoring founders required before accepting as uniquely highest scoring: %d",MinKmerScoreDiff);

		gDiagnostics.DiagOutMsgOnly (eDLInfo, "Accumulating Kmer alignments in non-overlapping bins of size : %dKbp", AccumBinSize);
		gDiagnostics.DiagOutMsgOnly (eDLInfo, "Using mean of raw haplotypes in sliding window containing this many bins to derived smoothed haplotypes : %d", WinMeanBins);
		gDiagnostics.DiagOutMsgOnly (eDLInfo, "Processing for reduction of founder haplotypes down to : %d",MaxFndrHaps);
		gDiagnostics.DiagOutMsgOnly (eDLInfo, "Processing for organism having ploidy : %d",Ploidy);

		if(szMaskBPAFile[0] != '\0')
			gDiagnostics.DiagOutMsgOnly (eDLInfo, "Masking BPAs loaded from : '%s'", szMaskBPAFile);
		}

	for(Idx = 0; Idx < NumFounderInputFiles; Idx++)
		gDiagnostics.DiagOutMsgOnly (eDLInfo, "Founder file : '%s'", pszFounderInputFiles[Idx]);

	if(PMode == eMSHDefault)
		for(Idx = 0; Idx < NumSkimInputFiles; Idx++)
			gDiagnostics.DiagOutMsgOnly (eDLInfo, "Skim file : '%s'", pszSkimInputFiles[Idx]);
	gDiagnostics.DiagOutMsgOnly (eDLInfo, "Output file : '%s'", szOutFile);
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
					MinKmerSize,			// aligning founders using this minimum Kmer size
					MaxKmerSize,			// aligning founders using this maximum Kmer size
					MinAlignedScoringPBA,		// to be accepted as aligned and scored then PBAs in both founder and skim must be at least this PBA level
					MinAlignedPBA,			// to be accepted as aligned ( < m_MinScoringPBA) with no penalty applied then PBAs in both founder and skim must be at least this PBA level
					AlignedPBAScore,		// if PBA accepted for scoring (>= m_MinAlignedScoringPBA) then add this score to the summed Kmer alignment score
					PenaltyPBAScore,		// if PBA not accepted as aligned ( < m_MinAlignedPBA) then apply this penalty to the summed Kmer alignment score
					AcceptKmerMatchPerc,	// only accept Kmer alignment if at least this percentage of base alleles matched at >= m_MinAlignedPBA
					MinKmerScoreDiff,		// must be at least this Kmer percentage score differential between highest scoring founders before accepting a founder as being a uniquely highest scoring
					AccumBinSize*1000,		// now use bp rather than Kbp! Counts are accumulated into this sized non-overlapping bins
					WinMeanBins,			// derive haplotype calls from mean of counts in a sliding window containing this many bins
					MaxFndrHaps,			// process for reduction down to this maximum number of called haplotypes
					Ploidy,					// smoothed haplotypes for organism having this ploidy
					szMaskBPAFile,			// optional input masking BPA file, only process BPAs which are intersect of these BPAs and skim plus founder BPAs
					NumFounderInputFiles,	// number of input founder file specs
					pszFounderInputFiles,	// names of input founder PBA files (wildcards allowed)
					NumSkimInputFiles,		// number of input skim file specs
					pszSkimInputFiles,		// names of input skim PBA files (wildcards allowed)
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
			 int MinKmerSize,				// aligning founders using this minimum Kmer size
			 int MaxKmerSize,				// aligning founders using this maximum Kmer size
			int32_t MinAlignedScoringPBA,		// to be accepted as aligned and scored then PBAs in both founder and skim must be at least this PBA level
			int32_t MinAlignedPBA,		// to be accepted as aligned ( < m_MinScoringPBA) with no penalty applied then PBAs in both founder and skim must be at least this PBA level
			int32_t AlignedPBAScore,		// if PBA accepted for scoring (>= m_MinAlignedScoringPBA) then add this score to the summed Kmer alignment score
			int32_t PenaltyPBAScore,		// if PBA not accepted as aligned ( < m_MinAlignedPBA) then apply this penalty to the summed Kmer alignment score
			int32_t AcceptKmerMatchPerc,	// only accept Kmer alignment if at least this percentage of base alleles matched at >= m_MinAlignedPBA
			int32_t MinKmerScoreDiff,		// must be at least this Kmer percentage score differential between highest scoring founders before accepting a founder as being a uniquely highest scoring
			int AccumBinSize,		// accumulating counts into this sized (bp) non-overlapping bins
			int WinMeanBins,		// derive haplotype calls from mean of counts in a sliding window containing this many bins
			int MaxFndrHaps,		// process for reduction down to this maximum number of called haplotypes
			int Ploidy,				// smoothed haplotypes for organism having this ploidy
			char *pszMaskBPAFile,	// optional input masking BPA file, only process BPAs which are intersect of these BPAs and skim plus founder BPAs
			int NumFounderInputFiles,	// number of input founder file specs
			char *pszFounderInputFiles[],	// names of input founder PBA files (wildcards allowed)
			int NumSkimInputFiles,	// number of input skim file specs
			char* pszSkimInputFiles[],		// names of input skim PBA files (wildcards allowed)
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
Rslt = pCallHaplotypes->Process(PMode,NChroms,pszTrackName,pszTrackDescr,MinKmerSize,MaxKmerSize,
				MinAlignedScoringPBA,MinAlignedPBA,AlignedPBAScore,PenaltyPBAScore,AcceptKmerMatchPerc,MinKmerScoreDiff,
				AccumBinSize,WinMeanBins,MaxFndrHaps,Ploidy,pszMaskBPAFile,NumFounderInputFiles,
				pszFounderInputFiles,NumSkimInputFiles,pszSkimInputFiles,pszOutFile,NumThreads);
delete pCallHaplotypes;
return(Rslt);
}


CCallHaplotypes::CCallHaplotypes()
{
m_pChromMetadata = NULL;
m_pWinBinCnts = NULL;
m_pWorkQueueEls = NULL;
m_pInBuffer = NULL;
m_pszOutBuffer = NULL;
m_hOutFile = -1;
m_hInFile = -1;
m_bMutexesCreated = false;
Reset();
}

CCallHaplotypes::~CCallHaplotypes()
{
if(m_hInFile != -1)
	close(m_hInFile);
if(m_hOutFile != -1)
	close(m_hOutFile);
if(m_pWorkQueueEls != NULL)
	delete []m_pWorkQueueEls;
if(m_pszOutBuffer != NULL)
	delete []m_pszOutBuffer;
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

if(m_pszOutBuffer != NULL)
	{
	delete []m_pszOutBuffer;
	m_pszOutBuffer = NULL;
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
m_Ploidy = 2;
m_MinScoreKmerSize = cDfltMinKmerSize;
m_MaxScoreKmerSize = min(cDfltMinKmerSize * 10,cMaxKmerSize);

m_MinAlignedScoringPBA = cDfltMinAlignedScoringPBA;
m_MinAlignedPBA = cDfltMinAlignedPBA;
m_AlignedPBAScore = cDfltAlignedPBAScore;
m_PenaltyPBAScore = cDfltPenaltyPBAScore;
m_AcceptKmerMatchPerc = cDfltAcceptKmerMatchPerc;
m_MinKmerScoreDiff = cDfltMinKmerScoreDiff;
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
m_MaxFndrHaps = 2;
m_AccumBinSize = cDfltAccumBinSize * 1000;

m_FndrsProcMap = (uint64_t)-1;

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
			int MinKmerSize,				// aligning founders using this minimum Kmer size
			int MaxKmerSize,				// aligning founders using this maximum Kmer size
			int32_t MinAlignedScoringPBA,		// to be accepted as aligned and scored then PBAs in both founder and skim must be at least this PBA level
			int32_t MinAlignedPBA,		// to be accepted as aligned ( < m_MinScoringPBA) with no penalty applied then PBAs in both founder and skim must be at least this PBA level
			int32_t AlignedPBAScore,		// if PBA accepted for scoring (>= m_MinAlignedScoringPBA) then add this score to the summed Kmer alignment score
			int32_t PenaltyPBAScore,		// if PBA not accepted as aligned ( < m_MinAlignedPBA) then apply this penalty to the summed Kmer alignment score
			int32_t AcceptKmerMatchPerc,	// only accept Kmer alignment if at least this percentage of base alleles matched at >= m_MinAlignedPBA
			int32_t MinKmerScoreDiff,		// must be at least this Kmer percentage score differential between highest scoring founders before accepting a founder as being a uniquely highest scoring
			uint32_t AccumBinSize,	// accumulate counts into this sized (bp) non-overlapping bins
			int WinMeanBins,		// derive haplotype calls from mean of counts in a sliding window containing this many bins
			int MaxFndrHaps,		// process for reduction down to this maximum number of called haplotypes
			int Ploidy,				// smoothed haplotypes for organism having this ploidy
			char *pszMaskBPAFile,	// optional input masking BPA file, only process BPAs which are intersect of these BPAs and skim plus founder BPAs
			int NumFounderInputFiles,	// number of input founder file specs
			char *pszFounderInputFiles[],	// names of input founder PBA files (wildcards allowed)
			int NumSkimInputFiles,	// number of input skim file specs
			char* pszSkimInputFiles[],		// names of input skim PBA files (wildcards allowed)
			char* pszOutFile,		// output to this files with this prefix
			int NumThreads)			// number of worker threads to use
{
int Rslt;
int NumFiles;
int TotNumFiles;
size_t memreq;

Reset();

CreateMutexes();
m_MaxFndrHaps = MaxFndrHaps;
m_Ploidy = Ploidy;
m_AccumBinSize = AccumBinSize;
m_WinMeanBins = WinMeanBins;
m_MinScoreKmerSize = MinKmerSize;
m_MaxScoreKmerSize = MaxKmerSize;
m_MinAlignedScoringPBA = MinAlignedScoringPBA;		// to be accepted as aligned and scored then PBAs in both founder and skim must be at least this PBA level
m_MinAlignedPBA = MinAlignedPBA;		// to be accepted as aligned ( < m_MinScoringPBA) with no penalty applied then PBAs in both founder and skim must be at least this PBA level
m_AlignedPBAScore = AlignedPBAScore;		// if PBA accepted for scoring (>= m_MinAlignedScoringPBA) then add this score to the summed Kmer alignment score
m_PenaltyPBAScore = PenaltyPBAScore;		// if PBA not accepted as aligned ( < m_MinAlignedPBA) then apply this penalty to the summed Kmer alignment score
m_AcceptKmerMatchPerc = AcceptKmerMatchPerc;	// only accept Kmer alignment if at least this percentage of base alleles matched at >= m_MinAlignedPBA
m_MinKmerScoreDiff = MinKmerScoreDiff;		// must be at least this Kmer percentage score differential between highest scoring founders before accepting a founder as being a uniquely highest scoring

m_NumThreads = NumThreads;
m_NChroms = NChroms;

if(PMode == eMSHAlleles)
	{
	Rslt = eBSFSuccess;		// assume success!
	CSimpleGlob glob(SG_GLOB_FULLSORT);
	int Idx;
	char* pszInFile;
	uint32_t ReadsetID = 0;
	TotNumFiles = 0;
	for(Idx = 0; Idx < NumFounderInputFiles; Idx++)	
		{
		glob.Init();
		if(glob.Add(pszFounderInputFiles[Idx]) < SG_SUCCESS)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to glob '%s' PBA Allele reporting files", pszFounderInputFiles[Idx]);
			Reset();
			return(eBSFerrOpnFile);	// treat as though unable to open file
			}
		if((NumFiles = glob.FileCount()) <= 0)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to locate any PBA Allele file matching '%s", pszFounderInputFiles[Idx]);
			Reset();
			return(eBSFerrFileName);
			}

		Rslt = eBSFSuccess;
		for(int FileID = 0; Rslt >= eBSFSuccess && FileID < NumFiles; ++FileID)
			{
			pszInFile = glob.File(FileID);
			if((ReadsetID = PBAReportAlleles(pszInFile,pszOutFile)) < 0)
				{
				gDiagnostics.DiagOut(eDLInfo, gszProcName, "Process: Errors loading Allele reporting file '%s'",pszInFile);
				Reset();
				return(-1);
				}
			}
		}
	Reset();
	return(0);
	}

if(pszMaskBPAFile == NULL || pszMaskBPAFile[0] == '\0')
	m_pszMaskBPAFile = NULL;
else
	m_pszMaskBPAFile = pszMaskBPAFile;

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

gDiagnostics.DiagOut(eDLInfo, gszProcName, "Process: Starting to load pool of founders");
Rslt = eBSFSuccess;		// assume success!
CSimpleGlob glob(SG_GLOB_FULLSORT);
int Idx;
char* pszInFile;
uint32_t ReadsetID = 0;
TotNumFiles = 0;
for(Idx = 0; Idx < NumFounderInputFiles; Idx++)	
	{
	glob.Init();
	if(glob.Add(pszFounderInputFiles[Idx]) < SG_SUCCESS)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to glob '%s' founder pool files", pszFounderInputFiles[Idx]);
		Reset();
		return(eBSFerrOpnFile);	// treat as though unable to open file
		}
	if((NumFiles = glob.FileCount()) <= 0)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to locate any founder pool file matching '%s", pszFounderInputFiles[Idx]);
		Reset();
		return(eBSFerrFileName);
		}
	TotNumFiles += NumFiles;
	if(TotNumFiles > cMaxFounderReadsets)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Can accept at most %d founder pool readsets for processing, after wildcard file name expansions there are %d requested", cMaxFounderReadsets, TotNumFiles);
		Reset();
		return(eBSFerrMem);
		}

	Rslt = eBSFSuccess;
	for(int FileID = 0; Rslt >= eBSFSuccess && FileID < NumFiles; ++FileID)
		{
		pszInFile = glob.File(FileID);
		if((ReadsetID = LoadPBAFile(pszInFile,0)) == 0)
			{
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "Process: Errors loading pool file '%s'",pszInFile);
			Reset();
			return(-1);
			}
		}
	}
m_NumFndrs = ReadsetID;
gDiagnostics.DiagOut(eDLInfo, gszProcName, "Process: Completed loading founders (%d) pool",m_NumFndrs);

// load the control PBA file if it has been specified
m_MaskReadsetID = 0;
if(m_pszMaskBPAFile != NULL)
	{
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Process: Loading control PBA file '%s'",m_pszMaskBPAFile);
	if((ReadsetID = LoadPBAFile(m_pszMaskBPAFile,2)) == 0)
		{
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Process: Errors loading control PBA file");
		Reset();
		return(-1);
		}
	m_MaskReadsetID = ReadsetID;
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Completed loading control PBA file");
	}

// next is to individually process each skim PBAs against the founder panel PBAs
m_SkimReadsetID = 0;
for(Idx = 0; Idx < NumSkimInputFiles; Idx++)	
	{
	glob.Init();
	if(glob.Add(pszSkimInputFiles[Idx]) < SG_SUCCESS)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to glob skim '%s", pszSkimInputFiles[Idx]);
		Reset();
		return(eBSFerrOpnFile);	// treat as though unable to open file
		}
	if((NumFiles = glob.FileCount()) <= 0)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to locate any skim PBA file matching '%s", pszSkimInputFiles[Idx]);
		Reset();
		return(eBSFerrOpnFile);	// treat as though unable to open file
		}

	Rslt = eBSFSuccess;
	for(int FileID = 0; Rslt >= eBSFSuccess && FileID < NumFiles; ++FileID)
		{
		pszInFile = glob.File(FileID);
		if((Rslt = ProcessSkimPBAFile(pszInFile,pszOutFile)) < eBSFSuccess)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Failed processing skim PBA file '%s", pszInFile);
			Reset();
			return(Rslt);
			}
		}
	}

Reset();
return(Rslt);
}

int
CCallHaplotypes::PBAReportAlleles(char* pszPBAFile,		// load and report allele details present in this PBA file
						char *pszRsltsFileBaseName)		// results are written to this file base name with skim readset identifier and type appended
{
uint32_t ReadsetID;

if(m_pInBuffer == NULL)
	{
	if((m_pInBuffer = new uint8_t[cInBuffSize]) == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Memory allocation of %u bytes for input buffering failed - %s", cInBuffSize, strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
	m_AllocInBuff = cInBuffSize;
	}

if(m_pChromMetadata == NULL)
	{
	// initial allocation, will be realloc'd as required if more memory required
	size_t memreq = (size_t)cAllocChromMetadata * sizeof(tsChromMetadata);
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
	m_AllocdChromMetadata = cAllocChromMetadata;
	}
m_UsedNumChromMetadata = 0;
m_InNumBuffered = 0;
m_InNumProcessed = 0;
m_LAReadsetNameID=0;
m_NumReadsetNames=0;
m_NxtszReadsetIdx=0;
m_LAChromNameID=0;
m_NumChromNames=0;
m_NxtszChromIdx=0;

if((ReadsetID = LoadPBAFile(pszPBAFile,0)) == 0)
	{
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Process: Errors loading PBA file '%s'",pszPBAFile);
	Reset();
	return(-1);
	}
m_NumFndrs = ReadsetID;

// iterate over each chromosome
uint32_t CurChromMetadataIdx;
tsReadsetMetadata *pReadsetMetadata;
tsChromMetadata *pChromMetadata;
uint32_t LociOfs;
uint8_t *pPBAs;
uint8_t *pPBA;
int32_t SkimBaseAllele;
uint32_t AlleleIdx;
uint8_t AlleleMsk;
uint32_t NumNonAligned;
pReadsetMetadata = &m_Readsets[ReadsetID-1];
CurChromMetadataIdx = pReadsetMetadata->StartChromMetadataIdx;

for(uint32_t ChromIdx = 0; ChromIdx < pReadsetMetadata->NumChroms && CurChromMetadataIdx != 0; ChromIdx++)
	{
	pChromMetadata = &m_pChromMetadata[CurChromMetadataIdx - 1];
	pPBA = pPBAs = LocatePBAfor(ReadsetID,pChromMetadata->ChromID);
	NumNonAligned = 0;
	for(LociOfs = 0; LociOfs < pChromMetadata->ChromLen; LociOfs++,pPBA++)
		{
		if(*pPBA == 0)
			NumNonAligned++;
		AlleleMsk = 0x03;
		for(AlleleIdx = 0; AlleleIdx < 4; AlleleIdx++, AlleleMsk <<= 2)
			{
			SkimBaseAllele =  (*pPBA & AlleleMsk) >> (AlleleIdx * 2);
			
			}
		}
	CurChromMetadataIdx = pChromMetadata->NxtChromMetadataIdx;
	}

 

if(m_hInFile != -1)
	{
	close(m_hInFile);
	m_hInFile = -1;
	}
return(eBSFSuccess);
}

int
CCallHaplotypes::ProcessSkimPBAFile(char* pszSkimPBAFile,	// load, process and call haplotypes for this skim PBA file against previously loaded panel founder PBAs
						char *pszRsltsFileBaseName)			// results are written to this file base name with skim readset identifier and type appended
{
int Rslt;
char szOutFile[_MAX_PATH];
char *pszFndrReadset;
uint32_t FounderID;
// mark state so can be restored ready for next skim PBA to be loaded and processed
uint32_t NumReadsetNames = m_NumReadsetNames;
uint32_t NxtszReadsetIdx = m_NxtszReadsetIdx;
uint32_t NumChromNames = m_NumChromNames;
uint32_t NxtszChromIdx = m_NxtszChromIdx;
uint32_t UsedNumChromMetadata = m_UsedNumChromMetadata;
uint32_t UsedWinBinCnts = 0;
m_SkimReadsetID = 0;
m_LAChromNameID = 0;
m_LAReadsetNameID = 0;
if((m_SkimReadsetID = LoadPBAFile(pszSkimPBAFile,1)) <= 0)
	{
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessSkimPBAFile: Errors loading skim PBA file '%s'",pszSkimPBAFile);
	Reset();
	return(-1);
	}
gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessSkimPBAFile: Completed loading skim PBA file");

if((m_pszOutBuffer = new uint8_t[cOutBuffSize]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessSkimPBAFile: Memory allocation of %u bytes for output buffering failed - %s", cOutBuffSize, strerror(errno));
	Reset();
	return(eBSFerrMem);
	}
m_AllocOutBuff = cOutBuffSize;
m_OutBuffIdx = 0;

char *pszSkimReadset = LocateReadset(m_SkimReadsetID);
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
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessSkimPBAFile: Unable to create/truncate %s - %s",szOutFile,strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}

tsReadsetMetadata *pSkim = &m_Readsets[m_SkimReadsetID-1];
pSkim->StartWinBinCntsIdx = m_UsedWinBinCnts;
pSkim->NumWinBins = 0;
pSkim->LRAWinBinCntIdx = m_UsedWinBinCnts;
Rslt = CountHaplotypes(m_SkimReadsetID,m_NumFndrs,m_MinScoreKmerSize,m_MaxScoreKmerSize,m_AccumBinSize);
pSkim->NumWinBins = m_UsedWinBinCnts - pSkim->StartWinBinCntsIdx;
gDiagnostics.DiagOut(eDLInfo,gszProcName,"ProcessSkimPBAFile: Calling haplotypes ...");
ChooseHaplotypes();  // choose which set of haplotypes are most probable for each bin
gDiagnostics.DiagOut(eDLInfo,gszProcName,"ProcessSkimPBAFile: Completed haplotype calling");
// report

m_OutBuffIdx=sprintf((char *)m_pszOutBuffer,"\"Chrom\",\"Loci\",\"Skim-NonUniques:%s\",\"Skim-Unaligned:%s\"",pszSkimReadset,pszSkimReadset);
for(FounderID = 1; FounderID <= m_NumFndrs; FounderID++)
	{
	if(m_Fndrs2Proc[FounderID-1] &  0x01)
		{
		pszFndrReadset = LocateReadset(FounderID);
		m_OutBuffIdx+=sprintf((char *)&m_pszOutBuffer[m_OutBuffIdx],",\"FndrUniques-%s:%s\"",pszSkimReadset,pszFndrReadset);
		}
	}

for(FounderID = 1; FounderID <= m_NumFndrs; FounderID++)
	{
	if(m_Fndrs2Proc[FounderID-1] &  0x01)
		{
		pszFndrReadset = LocateReadset(FounderID);
		m_OutBuffIdx+=sprintf((char *)&m_pszOutBuffer[m_OutBuffIdx],",\"SummedScore-%s:%s\"",pszSkimReadset,pszFndrReadset);
		}
	}

for(FounderID = 1; FounderID <= m_NumFndrs; FounderID++)
	{
	if(m_Fndrs2Proc[FounderID-1] &  0x01)
		{
		pszFndrReadset = LocateReadset(FounderID);
		m_OutBuffIdx+=sprintf((char *)&m_pszOutBuffer[m_OutBuffIdx],",\"RawHap-%s:%s\"",pszSkimReadset,pszFndrReadset);
		}
	}

for(FounderID = 1; FounderID <= m_NumFndrs; FounderID++)
	{
	if(m_Fndrs2Proc[FounderID-1] &  0x01)
		{
		pszFndrReadset = LocateReadset(FounderID);
		m_OutBuffIdx+=sprintf((char *)&m_pszOutBuffer[m_OutBuffIdx],",\"SmthdHap-%s:%s\"",pszSkimReadset,pszFndrReadset);
		}
	}

tsBinCnts *pSkimWinBinCnts;
uint32_t WndFndrCntsIdx;
uint32_t NumWinBins = pSkim->NumWinBins;
pSkimWinBinCnts = &m_pWinBinCnts[pSkim->StartWinBinCntsIdx];
for(WndFndrCntsIdx = 0; WndFndrCntsIdx < NumWinBins; WndFndrCntsIdx++, pSkimWinBinCnts++)
	{
	m_OutBuffIdx+=sprintf((char *)&m_pszOutBuffer[m_OutBuffIdx],"\n\"%s\",%u,%u,%u",
							LocateChrom(pSkimWinBinCnts->ChromID),pSkimWinBinCnts->StartLoci,pSkimWinBinCnts->MultiFounder,pSkimWinBinCnts->NoFounder);

	for(FounderID = 1; FounderID <= m_NumFndrs; FounderID++)
		if(m_Fndrs2Proc[FounderID-1] &  0x01)
			m_OutBuffIdx+=sprintf((char *)&m_pszOutBuffer[m_OutBuffIdx],",%u",pSkimWinBinCnts->Uniques[FounderID-1]);

	for(FounderID = 1; FounderID <= m_NumFndrs; FounderID++)
		if(m_Fndrs2Proc[FounderID-1] &  0x01)
			m_OutBuffIdx+=sprintf((char *)&m_pszOutBuffer[m_OutBuffIdx],",%lld",pSkimWinBinCnts->SummedScores[FounderID-1]);

	for(FounderID = 1; FounderID <= m_NumFndrs; FounderID++)
		if(m_Fndrs2Proc[FounderID-1] &  0x01)
			m_OutBuffIdx+=sprintf((char *)&m_pszOutBuffer[m_OutBuffIdx],",%u",pSkimWinBinCnts->RawFndrHaps[FounderID-1]);
	for(FounderID = 1; FounderID <= m_NumFndrs; FounderID++)
		if(m_Fndrs2Proc[FounderID-1] &  0x01)
			m_OutBuffIdx+=sprintf((char *)&m_pszOutBuffer[m_OutBuffIdx],",%u",pSkimWinBinCnts->SmthdFndrHaps[FounderID-1]);

	if((m_OutBuffIdx + 1000) > m_AllocOutBuff)
		{
		if (!CUtility::RetryWrites(m_hOutFile,m_pszOutBuffer,m_OutBuffIdx))
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
	if (!CUtility::RetryWrites(m_hOutFile,m_pszOutBuffer,m_OutBuffIdx))
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
for(FounderID = 1; FounderID <= m_NumFndrs; FounderID++)
	{
	if(m_Fndrs2Proc[FounderID-1] &  0x01)
		{
		pszFndrReadset = LocateReadset(FounderID);
		sprintf(szOutFile,"%s_Haps.Skim_%s.Fndr_%s.bed",pszRsltsFileBaseName,pszSkimReadset,pszFndrReadset);
		ReportHaplotypesBED(szOutFile,pszFndrReadset,pszSkimReadset,FounderID);
		}
	}

	// generate WIGs enabling raw counts visualisations, one for each founder
for(FounderID = 1; FounderID <= m_NumFndrs; FounderID++)
	{
	if(m_Fndrs2Proc[FounderID-1] &  0x01)
		{
		pszFndrReadset = LocateReadset(FounderID);
		sprintf(szOutFile,"%s_Haps.Skim_%s.Fndr_%s.wig",pszRsltsFileBaseName,pszSkimReadset,pszFndrReadset);
		ReportCountsWIG(szOutFile,pszFndrReadset,pszSkimReadset,FounderID);
		}
	}

// restore state ready for next skim PBA to be loaded and processed
m_NumReadsetNames = NumReadsetNames;
m_NxtszReadsetIdx = NxtszReadsetIdx;
m_NumReadsetNames = NumReadsetNames;
m_NumChromNames = NumChromNames;
m_NxtszChromIdx = NxtszChromIdx;
m_UsedWinBinCnts = 0;
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
m_SkimReadsetID = 0;
return(eBSFSuccess);
}



int
CCallHaplotypes::SmoothBinsHaps(int Ploidy)		// iterate over all window bins, smoothing using adjacent bins to resolve low confidence bins whilst further reducing number of haplotypes in a bin down to this ploidy
{
uint32_t CurChromID;
uint32_t FndrIdx;
uint32_t SkimBinIdx;
uint32_t CurChromSkimBinIdx;
uint32_t SumSkimBinIdx;
uint32_t SmoothStartBinIdx;
uint32_t SmoothEndBinIdx;
uint32_t NumSmoothBins;
uint32_t HapIdx;
uint32_t NumAbsent;
double MeanAllFndrHaps;
uint32_t SumFndrHaps[cMaxFounderReadsets];
tsFndrSmthdMean FndrSmthdMeans[cMaxFounderReadsets];
uint8_t *pRawFndrHap;
uint8_t *pSmthdFndrHap;
tsReadsetMetadata *pSkim;
tsBinCnts *pSkimWinBinCnts;
tsBinCnts *pSumSkimWinBinCnts;

memset(SumFndrHaps,0,sizeof(SumFndrHaps));
pSkim = &m_Readsets[m_SkimReadsetID-1];
pSkimWinBinCnts = &m_pWinBinCnts[pSkim->StartWinBinCntsIdx];
CurChromID = 0;
CurChromSkimBinIdx = 0;
for(SkimBinIdx = 0; SkimBinIdx < pSkim->NumWinBins; SkimBinIdx++,pSkimWinBinCnts++)
	{
	if(CurChromID != pSkimWinBinCnts->ChromID)
		{
		CurChromID = pSkimWinBinCnts->ChromID;
		CurChromSkimBinIdx = SkimBinIdx;
		}
	memset(SumFndrHaps,0,sizeof(SumFndrHaps));
	NumSmoothBins=0;
	MeanAllFndrHaps = 0;
	SmoothStartBinIdx = (SkimBinIdx - CurChromSkimBinIdx) < (m_WinMeanBins)/2 ? CurChromSkimBinIdx : SkimBinIdx - (m_WinMeanBins)/2;
	SmoothEndBinIdx = min(SkimBinIdx + (m_WinMeanBins+1)/2,m_UsedWinBinCnts-1);
	pSumSkimWinBinCnts =  &m_pWinBinCnts[pSkim->StartWinBinCntsIdx + SmoothStartBinIdx];
	for(SumSkimBinIdx = SmoothStartBinIdx; SumSkimBinIdx <= SmoothEndBinIdx; SumSkimBinIdx++,pSumSkimWinBinCnts++)
		{
		if(pSumSkimWinBinCnts->ChromID != CurChromID)
			break;
		NumSmoothBins++;
		pRawFndrHap = pSumSkimWinBinCnts->RawFndrHaps;
		for(FndrIdx = 0; FndrIdx < m_NumFndrs; FndrIdx++,pRawFndrHap++)
			{
			if(!(m_Fndrs2Proc[FndrIdx] & 0x01))
				continue;
			SumFndrHaps[FndrIdx] += (uint32_t)*pRawFndrHap;
			MeanAllFndrHaps += (double)(uint32_t)*pRawFndrHap;
			}
		}

	if(NumSmoothBins)
		{
		memset(FndrSmthdMeans,0,sizeof(tsFndrSmthdMean) * m_NumFndrsProcMap);
		HapIdx = 0;
		NumAbsent = 0;
		for(FndrIdx = 0; FndrIdx < m_NumFndrs; FndrIdx++)
			{
			if(!(m_Fndrs2Proc[FndrIdx] & 0x01))
				continue;
			NumAbsent += pSkimWinBinCnts->RawFndrHaps[FndrIdx] == eHCAbsent ? 0 : 1;
			FndrSmthdMeans[HapIdx].SmthdMean = (double)SumFndrHaps[FndrIdx] / NumSmoothBins;
			FndrSmthdMeans[HapIdx++].FounderID = FndrIdx+1;
			}
		MeanAllFndrHaps /= (NumSmoothBins * m_NumFndrsProcMap);
		if(m_NumFndrsProcMap > 1)
			qsort(FndrSmthdMeans,m_NumFndrsProcMap,sizeof(FndrSmthdMeans[0]),SortFndrSmthdMeans);
		uint32_t PloidyIdx;
		for(PloidyIdx = 0;PloidyIdx < min((uint32_t)Ploidy,m_NumFndrsProcMap);PloidyIdx++)
			{
			pSmthdFndrHap = &pSkimWinBinCnts->SmthdFndrHaps[FndrSmthdMeans[PloidyIdx].FounderID-1];
			if((FndrSmthdMeans[PloidyIdx].SmthdMean * 1.5) >= MeanAllFndrHaps && 
					pSkimWinBinCnts->SummedScores[FndrSmthdMeans[PloidyIdx].FounderID-1]) // needing to identify as to if there was a deletion thus no haplotype even if smoothing!
				*pSmthdFndrHap = eHCHiConf;
			else
				*pSmthdFndrHap = eHCAbsent;
			}
		}
	}
return(0);
}


int
CCallHaplotypes::ChooseHaplotypes(void) // iterate over all the window bins and choose which haplotype to call for that bin interpolating those bins which were initially called as being indeterminate
{
// classify bins into haplotype/absence for each founder within individual bins
ClassifyBinHaps(m_NumFndrsProcMap);

// apply interpolation by smoothing over haplotypes
SmoothBinsHaps(m_Ploidy);
return(0);
}

int
CCallHaplotypes::ReportHaplotypesBED(char *pszOutFile,		// BED file to generate
							char *pszFounder,	// this was the founder haplotype
							char *pszSkim,		// against which this skim readset was aligned
					uint32_t FounderID)			// reporting on haplotypes for this founder
{
uint32_t NumWinBins;
tsBinCnts *pWinBinCnts;

if(m_pszOutBuffer == NULL)
	{
	if((m_pszOutBuffer = new uint8_t[cOutBuffSize]) == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Memory allocation of %u bytes for output buffering failed - %s", cOutBuffSize, strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
	m_AllocOutBuff = cOutBuffSize;
	}
m_OutBuffIdx = 0;

// NOTE: reusing output filehandle and buffering as these are not used concurrently with any other output file generation
if(m_hOutFile != -1)
	{
	close(m_hOutFile);
	m_hOutFile = -1;
	}
#ifdef _WIN32
m_hOutFile = open(pszOutFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE));
#else
	if((m_hOutFile = open64(pszOutFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE)) != -1)
		if(ftruncate(m_hOutFile,0)!=0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"ReportHaplotypesBED: Unable to create/truncate %s - %s",pszOutFile,strerror(errno));
			Reset();
			return(eBSFerrCreateFile);
			}
#endif
if(m_hOutFile < 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ReportHaplotypesBED: Unable to create/truncate %s - %s",pszOutFile,strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}
m_OutBuffIdx = sprintf((char *)m_pszOutBuffer,"track name=\"PBA Skim %s->Fndr %s\" description=\"PBA Skim readset %s aligned against PBA Founder %s\" useScore=0\n", pszSkim,pszFounder, pszSkim,pszFounder);

NumWinBins = m_Readsets[m_SkimReadsetID-1].NumWinBins;
pWinBinCnts = &m_pWinBinCnts[m_Readsets[m_SkimReadsetID-1].StartWinBinCntsIdx];
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
			m_OutBuffIdx+=sprintf((char *)&m_pszOutBuffer[m_OutBuffIdx],"%s\t%u\t%u\n",pszCurChrom,CurStartLoci,CurEndLoci);
		CurEndLoci = 0;
		CurChromID = pWinBinCnts->ChromID;
		pszCurChrom = LocateChrom(CurChromID);
		}

	switch(pWinBinCnts->SmthdFndrHaps[FounderID-1]) {
		case eHCHiConf:			// high confidence that founder has a haplotype
			if(CurEndLoci == 0)
				CurStartLoci = pWinBinCnts->StartLoci;
			CurEndLoci = pWinBinCnts->EndLoci + 1;
			break;

		default:						// have either eHCAbsent or eHCLoConf in founder having a haplotype
			if(CurEndLoci != 0)
				{
				m_OutBuffIdx+=sprintf((char *)&m_pszOutBuffer[m_OutBuffIdx],"%s\t%u\t%u\n",pszCurChrom,CurStartLoci,CurEndLoci);
				CurEndLoci = 0;
				}
			break;
		}

	if((m_OutBuffIdx + 1000) > m_AllocOutBuff)
		{
		if (!CUtility::RetryWrites(m_hOutFile,m_pszOutBuffer,m_OutBuffIdx))
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "ReportHaplotypesBED: Fatal error in RetryWrites()");
			Reset();
			return(eBSFerrFileAccess);
			}
		m_OutBuffIdx = 0;
		}
	}
if(CurEndLoci != 0)
	{
	m_OutBuffIdx+=sprintf((char *)&m_pszOutBuffer[m_OutBuffIdx],"%s\t%u\t%u\n",pszCurChrom,CurStartLoci,CurEndLoci);
	CurEndLoci = 0;
	}

if(m_OutBuffIdx && m_hOutFile != -1)
	{
	if (!CUtility::RetryWrites(m_hOutFile,m_pszOutBuffer,m_OutBuffIdx))
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "ReportHaplotypesBED: Fatal error in RetryWrites()");
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

void
CCallHaplotypes::InitialiseWIGSpan(void) // initialise WIG span vars to values corresponding to no spans having been previously reported
{
m_WIGChromID = 0;
m_WIGRptdChromID = 0;
m_WIGSpanLoci = 0;
m_WIGSpanLen = 0;
m_WIGRptdSpanLen = 0;
m_WIGSpanCnts = 0;
}

int
CCallHaplotypes::CompleteWIGSpan(bool bWrite)				// close off any current WIG span ready to start any subsequent span
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
CCallHaplotypes::AccumWIGBinCnts(uint32_t ChromID,	// accumulate wiggle counts into variableStep spans over this chromosome
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
CCallHaplotypes::AccumWIGCnts(uint32_t ChromID,	// accumulate wiggle counts into variableStep spans over this chromosome
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
CCallHaplotypes::ReportCountsWIG(char *pszOutFile,	// WIG file to generate
							char *pszFounder,	// this was the founder haplotype
							char *pszSkim,		// against which this skim readset was aligned
							uint32_t FounderID) // reporting on exclusive counts for this founder
{
uint32_t NumWinBins;
tsBinCnts *pWinBinCnts;


if(m_pszOutBuffer == NULL)
	{
	if((m_pszOutBuffer = new uint8_t[cOutBuffSize]) == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Memory allocation of %u bytes for output buffering failed - %s", cOutBuffSize, strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
	m_AllocOutBuff = cOutBuffSize;
	}
m_OutBuffIdx = 0;

// NOTE: reusing output filehandle and buffering as these are not used concurrently with any other output file generation
if(m_hOutFile != -1)
	{
	close(m_hOutFile);
	m_hOutFile = -1;
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
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ReportHaplotypesWIG: Unable to create/truncate %s - %s",pszOutFile,strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}
m_OutBuffIdx = sprintf((char *)m_pszOutBuffer,"track name=\"PBA Skim %s->Fndr %s\" description=\"PBA Skim readset %s aligned against PBA Founder %s\" useScore=1\n", pszSkim,pszFounder, pszSkim,pszFounder);

NumWinBins = m_Readsets[m_SkimReadsetID-1].NumWinBins;
pWinBinCnts = &m_pWinBinCnts[m_Readsets[m_SkimReadsetID-1].StartWinBinCntsIdx];

uint32_t WndFndrCntsIdx;
uint32_t CurChromID = 0;
uint32_t CurBinSpan = 0;
uint32_t PrevBinSpan = 0;
InitialiseWIGSpan();
for(WndFndrCntsIdx = 0; WndFndrCntsIdx < NumWinBins; WndFndrCntsIdx++, pWinBinCnts++)
	{
	CurBinSpan = pWinBinCnts->EndLoci-pWinBinCnts->StartLoci+1;
	AccumWIGBinCnts(pWinBinCnts->ChromID,pWinBinCnts->StartLoci + 1,pWinBinCnts->Uniques[FounderID-1],CurBinSpan); // WIG uses 1 based loci instead of the usual UCSC BED 0 based loci
	}
CompleteWIGSpan(true);
if(m_OutBuffIdx && m_hOutFile != -1)
	{
	if (!CUtility::RetryWrites(m_hOutFile,m_pszOutBuffer,m_OutBuffIdx))
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "ReportHaplotypesWIG: Fatal error in RetryWrites()");
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

int32_t					// returned readset identifier (1..n) or < 0 if errors
CCallHaplotypes::LoadPBAFile(char* pszFile,	// load chromosome metadata and PBA data from this file
							uint8_t ReadsetType)	// 0: founder, 1: skim, 2: control
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
memset(pReadsetMetadata,0,sizeof(*pReadsetMetadata));
pReadsetMetadata->ReadsetType = ReadsetType;
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
return(ReadsetID);
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
CCallHaplotypes::LocatePBAfor(int32_t ReadSetID,		// readset identifier 
			 int32_t ChromID)			// chrom identifier
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
	Rslt = ClassifyChromHaplotypes(pWorkQueueEl->ReadsetID,pWorkQueueEl->ChromID,pWorkQueueEl->StartLoci,pWorkQueueEl->EndLoci,pWorkQueueEl->LociIncr,pWorkQueueEl->ChromLen,pWorkQueueEl->AccumBinSize,pWorkQueueEl->NumFndrs,pWorkQueueEl->pPBAs,pWorkQueueEl->pMskPBA);
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
CCallHaplotypes::AlignPBAs(uint32_t ReadsetID,		// this is the readset being aligned to founders - could be a skim or founder
									uint32_t NumFndrs,			// number of founders to be processed against
									int MinKmerSize,			// aligning founders using this minimum Kmer size
									int MaxKmerSize,			// aligning founders using this maximum Kmer size
									uint32_t LociIncr,			// kmer loci starts are incremented by this submultiple of AccumBinSize, used when sampling. Typically 1, 2, 5, 10, 20, 25
									uint32_t AccumBinSize,		// comparing ReadsetID's PBA against founder PBAs then counting using this bin size
									uint32_t NumFndrsProc)		// this many founders in m_Fndrs2Proc
{
uint32_t FounderID;
tsReadsetMetadata *pReadsetMetadata;
tsChromMetadata *pChromMetadata;
tsWorkQueueEl *pWorkQueueEl;
uint32_t CurChromMetadataIdx;
uint8_t *pPBAs[cMaxFounderReadsets+1];							// additional is to allow for the skim readset PBAs
uint8_t *pMskPBA;
char *pszReadset = LocateReadset(ReadsetID);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"AlignPBAs: Starting to align '%s' against %d founders",pszReadset, NumFndrsProc);
m_UsedWinBinCnts = 0;
pReadsetMetadata = &m_Readsets[ReadsetID-1]; 
CurChromMetadataIdx = pReadsetMetadata->StartChromMetadataIdx;
m_MinScoreKmerSize = MinKmerSize;
m_MaxScoreKmerSize = MaxKmerSize;
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
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "AlignPBAs: No PBA for chromosome '%s' in founder '%s', skipping this chromosome", pszChrom, LocateReadset(FounderID));
			break;
			}
		}

	if(m_MaskReadsetID == 0)
		pMskPBA = NULL;
	else
		pMskPBA = LocatePBAfor(m_MaskReadsetID, pChromMetadata->ChromID);

	CurChromMetadataIdx = pChromMetadata->NxtChromMetadataIdx;
	if(FounderID <= NumFndrs || (m_MaskReadsetID != 0  && pMskPBA == NULL))
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
	while(StartLoci < (pChromMetadata->ChromLen - MinKmerSize))
		{
		pWorkQueueEl->ReadsetID = ReadsetID;
		pWorkQueueEl->NumFndrs = NumFndrs;
		pWorkQueueEl->ChromID = pChromMetadata->ChromID;
		pWorkQueueEl->StartLoci = StartLoci;
		pWorkQueueEl->EndLoci = EndLoci;
		pWorkQueueEl->LociIncr = LociIncr;
		pWorkQueueEl->ChromLen = pChromMetadata->ChromLen;
		pWorkQueueEl->AccumBinSize = AccumBinSize;
		pWorkQueueEl->pMskPBA = pMskPBA;
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
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"AlignPBAs: Aligning '%s' against %d founders",pszReadset, NumFndrsProc);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"AlignPBAs: Completed alignment of '%s' against %d founders",pszReadset, NumFndrsProc);
if(m_pWorkQueueEls != NULL)
	{
	delete []m_pWorkQueueEls;
	m_pWorkQueueEls = NULL;	
	m_NumQueueElsProcessed = 0;
	}
return(eBSFSuccess);
}



int
CCallHaplotypes::CountHaplotypes(uint32_t ReadsetID,			// this is the readset being aligned to founders - could be a skim or a founder used as a control
									uint32_t NumFndrs,			// number of founders to be processed against
									int MinKmerSize,			// aligning founders using this minimum Kmer size
									int MaxKmerSize,			// aligning founders using this maximum Kmer size
									uint32_t AccumBinSize)		// comparing ReadsetID's PBA against founder PBAs then counting using this bin size
{
uint32_t FounderIdx;
uint32_t FounderID;
uint32_t NumFndrsProc;
uint32_t PrevNumFndrs;
uint32_t LociIncr;
int64_t WinBinIdx;
tsBinCnts *pWinBinCnts;
tsFndrTotUniqueCnts FndrUniqueCnts[cMaxFounderReadsets];
char *pszReadset = LocateReadset(ReadsetID);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"CountHaplotypes: Aligning '%s' against pool of %d founders, reducing until maximum of %d founders remaining ...",pszReadset,NumFndrs,m_MaxFndrHaps);

memset(m_Fndrs2Proc,0x01,sizeof(m_Fndrs2Proc));
NumFndrsProc = NumFndrs;
while(NumFndrsProc)
	{
	if(NumFndrsProc <= m_MaxFndrHaps)
		LociIncr = 1;
	else
		LociIncr = 5;
	AlignPBAs(ReadsetID,	// this is the readset being aligned to founders - could be a skim or founder
		  NumFndrs,			// number of founders to be processed against
		 MinKmerSize,		// aligning founders using this minimum Kmer size
		 MaxKmerSize,		// aligning founders using this maximum Kmer size
		 LociIncr,			// kmer loci starts are incremented by this submultiple of AccumBinSize, used when sampling. Typically 1, 2, 5, 10, 20, 25
		 AccumBinSize,		// comparing ReadsetID's PBA against founder PBAs then counting using this bin size
		 NumFndrsProc);		// this many founders remaining in m_Fndrs2Proc

// report on number of unique alignments to each founder as a diagnostic aid
	memset(FndrUniqueCnts,0,sizeof(FndrUniqueCnts));
	pWinBinCnts = m_pWinBinCnts;
	FounderIdx = 0;
	for(WinBinIdx = 0; WinBinIdx < m_UsedWinBinCnts; WinBinIdx++,pWinBinCnts++)
		{
		FounderIdx = 0;
		for(FounderID = 1; FounderID <= NumFndrs; FounderID++)
			{
			if(!(m_Fndrs2Proc[FounderID-1] & 0x01))
				continue;
			FndrUniqueCnts[FounderIdx].FounderID = FounderID;
			FndrUniqueCnts[FounderIdx++].TotUniqueCnts += pWinBinCnts->Uniques[FounderID-1];
			}
		}
	if(FounderIdx > 1)
		qsort(FndrUniqueCnts,FounderIdx,sizeof(tsFndrTotUniqueCnts),SortFndrTotUniqueCnts);
	for(FounderIdx = 0; FounderIdx < NumFndrsProc; FounderIdx++)
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"CountHaplotypes: Founder '%s' has %llu unique alignments",LocateReadset(FndrUniqueCnts[FounderIdx].FounderID),FndrUniqueCnts[FounderIdx].TotUniqueCnts);

	if(NumFndrsProc <= m_MaxFndrHaps)	// allowing at most this many founders to be haplotyped on counts
		{
		m_NumFndrsProcMap = NumFndrsProc;
		break;
		}
	// locate members of panel with highest number of uniques
	memset(FndrUniqueCnts,0,sizeof(FndrUniqueCnts));
	pWinBinCnts = m_pWinBinCnts;
	FounderIdx = 0;
	for(WinBinIdx = 0; WinBinIdx < m_UsedWinBinCnts; WinBinIdx++,pWinBinCnts++)
		{
		FounderIdx = 0;
		for(FounderID = 1; FounderID <= NumFndrs; FounderID++)
			{
			if(!(m_Fndrs2Proc[FounderID-1] & 0x01))
				continue;
			FndrUniqueCnts[FounderIdx].FounderID = FounderID;
			FndrUniqueCnts[FounderIdx++].TotUniqueCnts += pWinBinCnts->Uniques[FounderID-1];
			}
		}

	if(FounderIdx > 1)
		qsort(FndrUniqueCnts,FounderIdx,sizeof(tsFndrTotUniqueCnts),SortFndrTotUniqueCnts);

	
	if(FndrUniqueCnts[0].TotUniqueCnts == 0)
		break;

	
	PrevNumFndrs = NumFndrsProc;
	NumFndrsProc = 0;
	memset(m_Fndrs2Proc,0,sizeof(m_Fndrs2Proc));
	for(FounderIdx = 0; FounderIdx < NumFndrs; FounderIdx++)
		{
		if(NumFndrsProc >= m_MaxFndrHaps)
			{
			if(FndrUniqueCnts[FounderIdx].TotUniqueCnts == 0)
				break;
			if(FndrUniqueCnts[FounderIdx].TotUniqueCnts < (FndrUniqueCnts[0].TotUniqueCnts / 50))
				break;
			if((NumFndrsProc + 1) >= PrevNumFndrs)
				break;
			}
		NumFndrsProc++;
		m_Fndrs2Proc[FndrUniqueCnts[FounderIdx].FounderID-1] = 0x01;
		}
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"CountHaplotypes: Completed aligning '%s' against pool of %d founders, reduced down to %d founders",pszReadset,NumFndrs,m_MaxFndrHaps);
m_mtqsort.SetMaxThreads(m_NumThreads);
m_mtqsort.qsort(m_pWinBinCnts,(int64_t)m_UsedWinBinCnts,sizeof(tsBinCnts),SortWinBinCnts);
return(eBSFSuccess);
}


int
CCallHaplotypes::ClassifyChromHaplotypes(uint32_t ReadsetID,			// this readset is used for identifying the readset being aligned to founders - could be a skim or founder
										uint32_t ChromID,				// processing is for this chromosome
										uint32_t StartLoci,				// starting from this loci inclusive
										uint32_t EndLoci,				// and ending at this loci inclusive
										uint32_t LociIncr,				// incrementing loci by this many bases - must be 1, 2, 5, 10, 20, 25 - of which AccumBinSize is a multiple. Used when sampling
										uint32_t ChromLen,				// chromosome length
										uint32_t AccumBinSize,			// accumulating counts into this sized non-overlapping bins
										uint32_t NumFndrs,				// number of founders to be processed against the skim (or founder) PBA at *pPBAs[0]
										uint8_t *pPBAs[],				// pPBAs[0] pts to chromosome skim PBA, followed by ptrs to each of the chromosome founder PBAs
										uint8_t *pMskPBA)				// pts to optional chromosome masking PBA, scoring only for segments contained in mask which are non-zero and where skim and founder PBAs also have an allele.
												// enables a founder to be processed as if a skim, restricting scoring Kmers to be same as if a skim
{
int Rslt;
tsBinCnts WinBinCnts;
uint32_t Loci;
uint32_t NumWinBinCnts = 0;
tsFndrsClassifications InferencedFounders;

memset(&WinBinCnts,0,sizeof(WinBinCnts));
WinBinCnts.ReadsetID = ReadsetID;
WinBinCnts.ChromID = ChromID;
WinBinCnts.StartLoci = StartLoci; 
for(Loci = StartLoci; Loci <=  EndLoci && Loci <= (ChromLen - m_MinScoreKmerSize); Loci+=LociIncr)
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

	if((Rslt = ScorePBAKmer(ReadsetID,&WinBinCnts, Loci, ChromLen, m_MinScoreKmerSize, m_MaxScoreKmerSize,NumFndrs, pPBAs,&InferencedFounders,pMskPBA)) < eBSFSuccess)
		return(Rslt);

	// classify founders from alignment scores
	if(Rslt == 0)	// only 0 if no PBA alignments were accepted for skim against any founder
		WinBinCnts.NoFounder++;
	else
		{
		uint32_t FounderID;
		bool bUnique = false;
		for(FounderID = 1; FounderID <= NumFndrs; FounderID++)
			{
			switch(InferencedFounders.Fndrs[FounderID-1] & 0x03)
				{
				case 0x03:		// uniquely aligned
					WinBinCnts.Uniques[FounderID-1]++;
					bUnique = true;
					break;
				case 0x02:		// shares top billing
				case 0x01:		// was aligned, but not uniquely
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


int				// < 0 if errors, otherwise success
CCallHaplotypes::ScorePBAKmer(uint32_t ReadsetID,	// this readset is used for identifying the readset being aligned to founders - could be a skim or founder
						tsBinCnts *pWinBinCnts,		// bin currently being processed, loci is contained within this bin
						uint32_t Loci,			// processing for Kmer haplotypes starting from this loci
						uint32_t SeqLen,		// sequence or chromosome is this length
						uint32_t MinKmerSize,	// Kmer to score on is of this minimum size
						uint32_t MaxKmerSize,	// Kmer to score on can be extended until this maximum size
						uint32_t NumFndrs,		// number of founders to be processed against the skim (or founder) PBA at *pPBAs[0]
						uint8_t *pPBAs[],		// pts to loci 0 of PBAs in readset order, pPBA[0] pts to the PBA to be scored by aligning against the founder PBAs
						tsFndrsClassifications *pClassifications,		// returned classifications for each processed founder
						uint8_t *pMskPBA)		// pts to optional masking PBA, scoring only for segments contained in mask which are non-zero and where skim and founder PBAs also have an allele.
												// enables a founder to be processed as if a skim, restricting scoring Kmers to be same as if a skim
{
uint32_t FounderID;
uint32_t NumFndrs2Proc;
uint32_t NumFndrAligningAlleles;	// keeping count of number of scoring alleles for this founder so can determine if this founders proportion of scoring alleles reaches the acceptance threshold
tsFndrsClassifications InferencedFounders;
uint8_t *pScoreLoci;
uint8_t *pFndrLoci;
uint32_t AlleleIdx;
uint8_t AlleleMsk;
uint32_t KmerOfs;
uint32_t KmerLoci;

// quick check for skim and all founders having called PBA at specified loci
if( pPBAs[0] == NULL)	// must be a PBA to be scored
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ScorePBAKmer: NULL PBA at pPBA[0]");
	return(-1);
	}
pScoreLoci = pPBAs[0] + Loci;
if(*pScoreLoci == 0)	// will be 0 if no score allele aligned to starting Kmer Loci
	return(0);			// can't score ...

// all founders must have a PBA and alleles at initial Loci
NumFndrs2Proc = 0;
for(FounderID = 1; FounderID <= NumFndrs; FounderID++)
	{
	if(!(m_Fndrs2Proc[FounderID-1] & 0x01))
		continue;
	NumFndrs2Proc++;
	if((pFndrLoci=pPBAs[FounderID]) == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"ScorePBAKmer: NULL PBA founder at pPBA[%u] ...",FounderID);
		return(-1);
		}
	if(*(pFndrLoci+Loci)==0)
		return(0);
	}
if(!NumFndrs2Proc)
	return(0);
// if a mask PBA then check if it has an allele at initial loci
if(pMskPBA != NULL && pMskPBA[Loci] == 0)
	return(0);

// use KmerSize as being a minimum, enabling extension out to a maximum of 10x KmerSize
// if initial KmerSize containing PBAs discovered then allowing a gap which is of maximum 2x KmerSize providing it is followed by another min KmerSize of PBAs
// intent here is to extend lengths to cover skims which are PE with non-overlapping inserts
uint32_t SkimGapLen = 0;
uint32_t SkimPBAsLen = 0;
uint32_t MarkKmerOfs = 0;

for(FounderID = 1; FounderID <= NumFndrs; FounderID++) 
	{
	if(!(m_Fndrs2Proc[FounderID-1] & 0x01))
		continue;
	pScoreLoci = pPBAs[0] + Loci;
	pFndrLoci = pPBAs[FounderID] + Loci;
	KmerOfs = 0;
	MarkKmerOfs = 0;
	SkimPBAsLen = 0;
	SkimGapLen = 0;
	for(KmerLoci = Loci, KmerOfs = 0; KmerLoci <= (SeqLen - MinKmerSize) && KmerOfs <= MaxKmerSize; KmerOfs++,KmerLoci++,pFndrLoci++,pScoreLoci++)
		{
		if(*pFndrLoci == 0)	// will be 0 if no fndr alleles aligned to this Loci within the Kmer
			break;
		if(pMskPBA != NULL && pMskPBA[KmerLoci] == 0)	// masks treated as if a founder
			break;
		if(*pScoreLoci == 0) // skims are allowed to have gaps as long as a gap follows a stretch of PBAs at least KmerSize'd	
			{
			SkimGapLen+=1;
			if(SkimPBAsLen < MinKmerSize || SkimGapLen > MaxKmerSize)
				break;
			}
		else
			{
			if(SkimGapLen > 0)
				{
				SkimGapLen = 0;
				SkimPBAsLen = 0;
				}
			SkimPBAsLen += 1;
			if(SkimPBAsLen >= MinKmerSize)
				MarkKmerOfs = KmerOfs;
			}
		}
	if(MarkKmerOfs < MinKmerSize)			// extension has to be at least original minimal KmerSize
		return(0);
	if(MarkKmerOfs < MaxKmerSize)		// reducing down to minimum extension of any founder or skim
		MaxKmerSize = MarkKmerOfs + 1;
	}
if(MaxKmerSize < MinKmerSize)			// extension has to be at least original minimal KmerSize
	return(0);
MinKmerSize = MaxKmerSize; 

bool bInGap;
int32_t FounderScores[cMaxFounderReadsets];
int32_t *pFounderScore;
memset(FounderScores,0,sizeof(FounderScores));
memset(&InferencedFounders,0,sizeof(InferencedFounders));

pFounderScore = FounderScores;
for(FounderID = 1; FounderID <= NumFndrs; FounderID++,pFounderScore++)
	{
	if(!(m_Fndrs2Proc[FounderID-1] & 0x01))
		continue;
	NumFndrAligningAlleles = 0;		
	pFndrLoci = pPBAs[FounderID] + Loci;
	pScoreLoci = pPBAs[0] + Loci;
	bInGap = false;
	KmerOfs = 0;
	for(KmerLoci = Loci, KmerOfs = 0; KmerLoci <= (SeqLen - MinKmerSize) && KmerOfs < MinKmerSize; KmerOfs++,KmerLoci++,pScoreLoci++,pFndrLoci++)
		{
		if(pMskPBA != NULL && pMskPBA[KmerLoci] == 0)
			break;
		if(*pFndrLoci == 0)	// will be 0 if no alleles aligned to this Loci within the KmerSize
			break;
		if(*pScoreLoci == 0) // if PBA gap for skim then apply a very small gap opening penalty 
			{
			if(!bInGap)
				{
				*pFounderScore -= 1;
				bInGap = true;
				}
			NumFndrAligningAlleles++;
			continue;
			}
		bInGap = false;
		AlleleMsk = 0x03;	
		int32_t FndrBaseAllele;
		int32_t SkimBaseAllele;
		bool bAlleleScored = false;
		bool bAlleleAligned = false;
		for(AlleleIdx = 0; AlleleIdx < 4; AlleleIdx++, AlleleMsk <<= 2)
			{
			SkimBaseAllele =  (*pScoreLoci & AlleleMsk) >> (AlleleIdx * 2);
			if(SkimBaseAllele < m_MinAlignedPBA)
				continue;
			FndrBaseAllele = (*pFndrLoci & AlleleMsk) >> (AlleleIdx * 2);
			if(FndrBaseAllele < m_MinAlignedPBA)
				continue;

			bAlleleAligned = true;
			if(SkimBaseAllele >= m_MinAlignedScoringPBA && FndrBaseAllele >= m_MinAlignedScoringPBA)
				{
				bAlleleScored = true;
				if(SkimBaseAllele == 3 && FndrBaseAllele == 3)
					{
					*pFounderScore += m_AlignedPBAScore; // really reward dirac matches!	
					break;
					}
				}
			}

		if(bAlleleAligned)
			{
			NumFndrAligningAlleles++;
			if(bAlleleScored)
				*pFounderScore += m_AlignedPBAScore;
			}
		else
			*pFounderScore += m_PenaltyPBAScore;
		}

	// did this founder reach minimum proportion of required allele matches for score acceptance?
	if((int32_t)(((NumFndrAligningAlleles * 100) + 99) / MinKmerSize) < m_AcceptKmerMatchPerc)
		*pFounderScore = 0;
	}

// classify from alignment scores
// if highest scorer then classify as 3
// if equal highest then classify as 2
// if not highest or equal highest, and accepted as aligned, then classify as 1
// if not accepted as aligned then classify as 0
uint8_t FndrHighest;
uint8_t FndrNxtHighest;
uint8_t FndrAccepted;
int32_t NumAccepted;
int32_t FndrHighestScore = 0;
int32_t FndrNxtHighestScore = 0;

// locate highest and second highest founder scores
pFounderScore = FounderScores;
for(FounderID = 1; FounderID <= NumFndrs; FounderID++,pFounderScore++)
	{
	if(!(m_Fndrs2Proc[FounderID-1] & 0x01))
		continue;
	if(*pFounderScore <= 0)			
		continue;					
	pWinBinCnts->SummedScores[FounderID-1] += *pFounderScore;

	if(*pFounderScore > FndrHighestScore)
		{
		FndrNxtHighestScore = FndrHighestScore;		// next highest becomes the previous highest
		FndrHighestScore = *pFounderScore;
		}
	else
		if(*pFounderScore <= FndrHighestScore && *pFounderScore > FndrNxtHighestScore)
			FndrNxtHighestScore = *pFounderScore;
	}

// requiring at least a m_MinKmerScoreDiff % differential between highest and next highest to reduce the impact of sequencing/alignment/noise artifacts, if less than m_MinKmerScoreDiff then treat as if equally scoring
if(FndrHighestScore != FndrNxtHighestScore && ((FndrHighestScore * 100) < (FndrNxtHighestScore * (100 + m_MinKmerScoreDiff))))
	FndrHighestScore = FndrNxtHighestScore;

FndrHighest = 0x03;
FndrNxtHighest = 0x02;
FndrAccepted = 0x01;
NumAccepted = 0;
pFounderScore = FounderScores;
for(FounderID = 1; FounderID <= NumFndrs;FounderID++,pFounderScore++)
	{
	InferencedFounders.Fndrs[FounderID-1]= 0;
	if(!(m_Fndrs2Proc[FounderID-1] & 0x01))
		continue;
	if(*pFounderScore <= 0)
		continue;
	NumAccepted++;
	if(FndrHighestScore != 0 && *pFounderScore >= FndrHighestScore && FndrHighestScore > FndrNxtHighestScore)
		InferencedFounders.Fndrs[FounderID-1]= FndrHighest;
	else
		if(*pFounderScore >= FndrNxtHighestScore)
			InferencedFounders.Fndrs[FounderID-1] = FndrNxtHighest;
		else
			InferencedFounders.Fndrs[FounderID-1] = FndrAccepted;
	}
*pClassifications = InferencedFounders;
return(NumAccepted);
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

int						
CCallHaplotypes::ClassifyBinHaps(int MaxHaps)  // classify each window bin as containing a combination of this number of maximum haplotypes in that bin
{
uint32_t FndrIdx;
uint32_t TotUniques;
uint32_t TotNonUniques;
double UniquesMean;
double NonUniquesMean;
uint32_t *pFndrUniques;
uint32_t *pFndrNonUniques;
uint8_t *pRawFndrHap;
uint32_t SkimBinIdx;
tsReadsetMetadata *pSkim;
tsBinCnts *pSkimWinBinCnts;

pSkim = &m_Readsets[m_SkimReadsetID-1];
pSkimWinBinCnts = &m_pWinBinCnts[pSkim->StartWinBinCntsIdx];
for(SkimBinIdx = 0; SkimBinIdx < pSkim->NumWinBins; SkimBinIdx++,pSkimWinBinCnts++)
	{
	memset(pSkimWinBinCnts->RawFndrHaps,(uint8_t)eHCAbsent,sizeof(pSkimWinBinCnts->RawFndrHaps));
	memset(pSkimWinBinCnts->SmthdFndrHaps,(uint8_t)eHCAbsent,sizeof(pSkimWinBinCnts->SmthdFndrHaps));
	// first determine the mean of the uniques over all founders
	TotUniques = 0;
	TotNonUniques = 0;

	pFndrUniques = pSkimWinBinCnts->Uniques;
	pFndrNonUniques = pSkimWinBinCnts->NonUniques;
	for(FndrIdx = 0; FndrIdx < m_NumFndrs; FndrIdx++,pFndrUniques++,pFndrNonUniques++)
		{
		if(!(m_Fndrs2Proc[FndrIdx] & 0x01))
			continue;
		TotUniques += *pFndrUniques;
		TotNonUniques += *pFndrNonUniques;
		}

	UniquesMean = (double)TotUniques/m_NumFndrsProcMap;
	NonUniquesMean = (double)TotNonUniques/m_NumFndrsProcMap;
	pFndrUniques = pSkimWinBinCnts->Uniques;
	pFndrNonUniques = pSkimWinBinCnts->NonUniques;
	pRawFndrHap = pSkimWinBinCnts->RawFndrHaps;
	for(FndrIdx = 0; FndrIdx < m_NumFndrs; FndrIdx++,pFndrUniques++,pFndrNonUniques++,pRawFndrHap++)
		{
		if(!(m_Fndrs2Proc[FndrIdx] & 0x01))
			continue;

		if(*pFndrUniques < 5 && *pFndrNonUniques < 100)
			continue;								// treating as noise counts, so leave ClassFndr as eHCHiConfNone

		if(*pFndrUniques < (UniquesMean/10) && *pFndrNonUniques < (NonUniquesMean/5))
			continue;								// treating as noise counts, so leave ClassFndr as eHCHiConfNone

		if(*pFndrUniques >= 2 && (*pFndrUniques * 5) >= UniquesMean)
			*pRawFndrHap = eHCHiConf;			// this founder has same or higher numbers of unique alignments relative to other founders 
		else
			if((*pFndrUniques * 10) >= UniquesMean)
				*pRawFndrHap = eHCLoConf;		// low confidence that haplotype is present
		}
	}
return(0);			
}



int32_t		// returned chrom identifier, < 1 if unable to accept this chromosome name
CCallHaplotypes::AddChrom(char* pszChrom) // associate unique identifier with this chromosome name
{
int32_t ChromNameIdx;
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


int32_t		// returned chrom identifier, < 1 if unable to locate this chromosome name
CCallHaplotypes::LocateChrom(char* pszChrom) // return unique identifier associated with this chromosome name
{
int32_t ChromNameIdx;
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
CCallHaplotypes::LocateChrom(int32_t ChromID)
{
if(ChromID < 1 || ChromID > m_NumChromNames)
	return(NULL);
return(&m_szChromNames[m_szChromIdx[ChromID-1]]);
}


// NOTE: Readsets are checked for uniqueness as readsets must be unique
int32_t		// returned readset identifier, < 1 if unable to accept this readset name
CCallHaplotypes::AddReadset(char* pszReadset) // associate unique identifier with this readset name
{
int32_t ReadsetNameIdx;
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


int32_t		// returned Readset identifier, < 1 if unable to locate this Readset name
CCallHaplotypes::LocateReadset(char* pszReadset) // return unique identifier associated with this Readset name
{
int32_t ReadsetNameIdx;
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
CCallHaplotypes::LocateReadset(int32_t ReadsetID)
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

// sorting founder unique counts descending, founder ascending
int
CCallHaplotypes::SortFndrTotUniqueCnts(const void* arg1, const void* arg2)
{
tsFndrTotUniqueCnts* pEl1 = (tsFndrTotUniqueCnts*)arg1;
tsFndrTotUniqueCnts* pEl2 = (tsFndrTotUniqueCnts*)arg2;
if(pEl1->TotUniqueCnts > pEl2->TotUniqueCnts)
	return(-1);
if(pEl1->TotUniqueCnts < pEl2->TotUniqueCnts)
	return(1);
if(pEl1->FounderID < pEl2->FounderID)
	return(-1);
if(pEl1->FounderID > pEl2->FounderID)
	return(1);
return(0);
}

// sorting founder smoothed means descending, founder ascending
int
CCallHaplotypes::SortFndrSmthdMeans(const void* arg1, const void* arg2)
{
tsFndrSmthdMean* pEl1 = (tsFndrSmthdMean*)arg1;
tsFndrSmthdMean* pEl2 = (tsFndrSmthdMean*)arg2;
if(pEl1->SmthdMean > pEl2->SmthdMean)
	return(-1);
if(pEl1->SmthdMean < pEl2->SmthdMean)
	return(1);
if(pEl1->FounderID < pEl2->FounderID)
	return(-1);
if(pEl1->FounderID > pEl2->FounderID)
	return(1);
return(0);
}
