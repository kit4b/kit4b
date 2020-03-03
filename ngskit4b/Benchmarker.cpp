// alignbench.cpp : This file contains the 'main' function. Program execution begins and ends there.
// Aligner benchmarker
// Benchmarks aligners
// a) For a given aligner processes the claimed alignments for a given real readset against a reference target
//    and derives an error profile for each claimed alignment in terms of at which base in the alignment there
//    was a mismatch or InDel relative to the reference target. The error profile is recorded for each individual aligned read,
//    and if PE alignments for the pair of reads including the insert or fragment size.
// b) The benchmark application then simulates a set of reads originating at known loci (the ground truth) with the simulated reads
//    having the same error pofiles applied as were present in the claimed alignments.
// c) These simulated reads are then aligned by a aligner undergoing benchmarking
// d) The subsequent alignments are then processed against the ground truth and scored.
//    Scoring is on each individual base, essentally a correctly aligned base is assigned positive points, an incorrectly aligned
//    base is assigned negative points, and an uncalled base is assigned neutral points.
//    Positive, neutral and negative base points are configurable.
//    The score for a given aligner is the sum of points assigned over all bases for the set of simulated reads.
// 
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
#include "Benchmarker.h"

int
BenchmarkProcess(eBMProcMode PMode,			// processing mode
				bool bPrimaryOnly,			// if true then score only primary alignments otherwise score including secondary
				bool bPEReads,				// if true then PE pairs processing otherwise SE reads
				int UnalignedBasePts,		// points to apply for each unaligned base
				int AlignedBasePts,			// points to apply for each base aligned to it's ground truth loci
				int SilentTrimAlignBasePts, // points to apply for each base aligned in a silently trimmed read within ground truth loci range
				int MisalignedBasePts,		// points to apply for each base aligned but not to ground truth loci
				int MaxNumReads,			// maximum number of alignment CIGARs to process or number of simulated reads 
				char *pszObsCIGARsFile,		// observed CIGARs are in this file
				char *pszRefGenomeFile,		// reads are against this target genome file
				char *pszAlignmentsFile,	// input file containing aligned reads (SAM or BAM)
				char *pszSEReads,			// simulated reads are output to this file (SE or PE1 if PE)  (input when scoring as contains ground truth)
				char *pszPE2Reads);			// simulated PE2 reads are output to this file if simulating PE reads (input when scoring as contains ground truth)



#ifdef _WIN32
int BenchmarkAligners(int argc, char* argv[])
{
	// determine my process name
_splitpath(argv[0], NULL, NULL, gszProcName, NULL);
#else
int
BenchmarkAligners(int argc, char** argv)
{
	// determine my process name
CUtility::splitpath((char*)argv[0], NULL, gszProcName);
#endif

int iFileLogLevel;			// level of file diagnostics
int iScreenLogLevel;		// level of file diagnostics
char szLogFile[_MAX_PATH];	// write diagnostics to this file
int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure

eBMProcMode PMode;			// processing mode
bool bPrimaryOnly;			// score only primary alignments otherwise score including secondary
bool bPEReads;				// if true then PE simulations

int NumberOfProcessors;		// number of installed CPUs

int UnalignedBasePts;		// points to apply for each unaligned base
int SilentTrimAlignBasePts; // points to apply for each base aligned in a silently trimmed read within ground truth loci range
int AlignedBasePts;			// points to apply for each base aligned to it's ground truth loci
int MisalignedBasePts;		// points to apply for each base aligned but not to ground truth loci

int MaxNumReads;			// maximum number of alignment CIGARs to process or number of simulated reads 

char szCIGARsFile[_MAX_PATH];	// observed CIGARs are in this file
char szRefGenomeFile[_MAX_PATH];// reads are against this target genome file
char szAlignmentsFile[_MAX_PATH];	// input file containing aligned reads (SAM or BAM)
char szSEReads[_MAX_PATH];		// simulated reads are output to this file (SE or PE1 if PE)
char szPE2Reads[_MAX_PATH];		// simulated PE2 reads are output to this file if simulating PE reads

char szSQLiteDatabase[_MAX_PATH];	// results summaries to this SQLite file
char szExperimentName[cMaxDatasetSpeciesChrom + 1];	// experiment name
char szExperimentDescr[1000];		// describes experiment

struct arg_lit* help = arg_lit0("h", "help", "print this help and exit");
struct arg_lit* version = arg_lit0("v", "version,ver", "print version information and exit");
struct arg_int* FileLogLevel = arg_int0("f", "FileLogLevel", "<int>", "Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file* LogFile = arg_file0("F", "log", "<file>", "diagnostics log file");

struct arg_int* pmode = arg_int0("m", "mode", "<int>", "processing mode: 0 - generate observed CIGARs from alignments, 1:- simulate reads using observed CIGARs, 2: score base alignments using expected CIGARs");
struct arg_file* inalignments = arg_file0("i", "alignments", "<file>", "input alignments (SAM or BAM) when generating observed CIGARs or alignments if scoring using expected CIGARs");

struct arg_lit* scoreprimaryonly = arg_lit0("P", "primary", "set if only primary alignments are to be scored");
struct arg_lit* pereads = arg_lit0("p", "pe",				"set if PE pairs processing otherwise SE reads only");

struct arg_file* refgenome = arg_file0("I", "refgenome", "<file>",	"alignment reference genome suffix array ('ngskit4b index' generated) file");
struct arg_file* reads = arg_file0("o", "reads", "<file>",			"simulated reads to this file SE or PE1 if PE reads, or from if scoring");
struct arg_file* pe2reads = arg_file0("O", "pe2reads", "<file>",	"simulated reads to this file if PE reads, or from if scoring");

struct arg_file* cigars = arg_file0("c", "CIGARs", "<file>",		"observed CIGARs R/W this file");

struct arg_int* alignedbasepts = arg_int0("j", "alignedbasepts", "<int>",		"correctly aligned base points (defaults to 2");
struct arg_int* silenttrimaligdbasepts = arg_int0("J", "silenttrimalignbasepts", "<int>", "silently trimmed alignment within truth bounds base points (defaults to 1)");
struct arg_int* unalignedbasepts = arg_int0("k", "unalignedbasepts", "<int>",	"unaligned base points (defaults to 0)");
struct arg_int* misalignedbasepts = arg_int0("l", "misalignedbasepts", "<int>",	"incorrectly aligned base points (defaults to -2");
struct arg_int* maxnumreads = arg_int0("r", "maxnumreads", "<int>", "process for this number of reads or read pairs if PE (defaults to 10000000");

struct arg_file* summrslts = arg_file0("q", "sumrslts", "<file>", "Output results summary to this SQLite3 database file");
struct arg_str* experimentname = arg_str0("w", "experimentname", "<str>", "experiment name SQLite3 database file");
struct arg_str* experimentdescr = arg_str0("W", "experimentdescr", "<str>", "experiment description SQLite3 database file");

struct arg_end* end = arg_end(200);

	void* argtable[] = { help,version,FileLogLevel,LogFile,
						pmode, inalignments,pereads,refgenome,reads,pe2reads,scoreprimaryonly,cigars,alignedbasepts,silenttrimaligdbasepts,unalignedbasepts,misalignedbasepts,maxnumreads,
						summrslts,experimentname,experimentdescr,end };

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
		gExperimentID = 0;
		gProcessID = 0;
		gProcessingID = 0;
		szSQLiteDatabase[0] = '\0';
		szExperimentName[0] = '\0';
		szExperimentDescr[0] = '\0';


		if (experimentname->count)
		{
			strncpy(szExperimentName, experimentname->sval[0], sizeof(szExperimentName));
			szExperimentName[sizeof(szExperimentName) - 1] = '\0';
			CUtility::TrimQuotedWhitespcExtd(szExperimentName);
			CUtility::ReduceWhitespace(szExperimentName);
		}

		if (strlen(szExperimentName) < 1)
			strcpy(szExperimentName, "N/A");

		if (experimentdescr->count)
		{
			strncpy(szExperimentDescr, experimentdescr->sval[0], sizeof(szExperimentDescr) - 1);
			szExperimentDescr[sizeof(szExperimentDescr) - 1] = '\0';
			CUtility::TrimQuotedWhitespcExtd(szExperimentDescr);
			CUtility::ReduceWhitespace(szExperimentDescr);
		}
		if (strlen(szExperimentDescr) < 1)
			strcpy(szExperimentDescr, "N/A");

		gExperimentID = 0;
		gProcessID = 0;
		gProcessingID = 0;
		szSQLiteDatabase[0] = '\0';

		if (summrslts->count)
		{
			strncpy(szSQLiteDatabase, summrslts->filename[0], sizeof(szSQLiteDatabase) - 1);
			szSQLiteDatabase[sizeof(szSQLiteDatabase) - 1] = '\0';
			CUtility::TrimQuotedWhitespcExtd(szSQLiteDatabase);
			if (strlen(szSQLiteDatabase) < 1)
			{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: After removal of whitespace, no SQLite database specified with '-q<filespec>' option");
				return(1);
			}

			gExperimentID = gSQLiteSummaries.StartExperiment(szSQLiteDatabase, false, true, szExperimentName, szExperimentName, szExperimentDescr);
			if (gExperimentID < 1)
				return(1);
			gProcessID = gSQLiteSummaries.AddProcess((char*)gpszSubProcess->pszName, (char*)gpszSubProcess->pszName, (char*)gpszSubProcess->pszFullDescr);
			if (gProcessID < 1)
				return(1);
			gProcessingID = gSQLiteSummaries.StartProcessing(gExperimentID, gProcessID, (char*)kit4bversion);
			if (gProcessingID < 1)
				return(1);
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "Initialised SQLite database '%s' for results summary collection", szSQLiteDatabase);
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "SQLite database experiment identifier for '%s' is %d", szExperimentName, gExperimentID);
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "SQLite database process identifier for '%s' is %d", (char*)gpszSubProcess->pszName, gProcessID);
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "SQLite database processing instance identifier is %d", gProcessingID);
		}
		else
			szSQLiteDatabase[0] = '\0';

		// ensure all filenames are initialised in case not user specified
		szCIGARsFile[0] = '\0';
		szRefGenomeFile[0] = '\0';
		szAlignmentsFile[0] = '\0';
		szSEReads[0] = '\0';
		szPE2Reads[0] = '\0';
		bPEReads = false;
		UnalignedBasePts = 0;
		AlignedBasePts = 0;
		MisalignedBasePts = 0;
		MaxNumReads = 0;
		SilentTrimAlignBasePts = 0;

		PMode = pmode->count ? (eBMProcMode)pmode->ival[0] : eBMGenCIGARs;
		if (PMode < eBMGenCIGARs || PMode > eBMScore)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Processing mode '-m%d' specified outside of range %d..%d\n", PMode, eBMGenCIGARs, (int)eBMScore);
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
		if (PMode == eBMGenCIGARs || PMode == eBMSimReads)
			{
			if (cigars->count)
				{
				strcpy(szCIGARsFile, cigars->filename[0]);
				CUtility::TrimQuotedWhitespcExtd(szCIGARsFile);
				}
			if (szCIGARsFile[0] == '\0')
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "No observed CIGARs file specified");
				exit(1);
				}
		}
		else
			szCIGARsFile[0] = '\0';

		if(PMode == eBMGenCIGARs || PMode == eBMSimReads)
			{
			if(refgenome->count)
				{
				strcpy(szRefGenomeFile, refgenome->filename[0]);
				CUtility::TrimQuotedWhitespcExtd(szRefGenomeFile);
				}
			if (szRefGenomeFile[0] == '\0')
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "No reference genome file specified");
				exit(1);
				}
			}
		else
			szRefGenomeFile[0] = '\0';

		if (PMode == eBMGenCIGARs || PMode == eBMScore)
			{
			if(inalignments->count)
				{
				strcpy(szAlignmentsFile, inalignments->filename[0]);
				CUtility::TrimQuotedWhitespcExtd(szAlignmentsFile);
				}
			if (szAlignmentsFile[0] == '\0')
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "No input alignment (SAM/BAM) file specified");
				exit(1);
				}
			}
		else
			szAlignmentsFile[0] = '\0';

		bPrimaryOnly = scoreprimaryonly->count ? true : false;
		bPEReads = pereads->count ? true : false;
		MaxNumReads = maxnumreads->count ? maxnumreads->ival[0] : cBMDfltReads;
		if(MaxNumReads > cBMMaxReads)		// silently clamp to be within a reasonable range
			MaxNumReads = cBMMaxReads;
		else
			if (MaxNumReads < cBMMinReads)
				MaxNumReads = cBMMinReads;

		if (PMode == eBMScore)
			{
			AlignedBasePts = alignedbasepts->count ? alignedbasepts->ival[0] : cBMDfltAlignedBasePts;
			if(AlignedBasePts > cBMMaxPts)
				AlignedBasePts = cBMMaxPts;
			else
				if(AlignedBasePts < 1)				// must always be > 0
					AlignedBasePts = 1;
			SilentTrimAlignBasePts = silenttrimaligdbasepts->count ? silenttrimaligdbasepts->ival[0] : cBMDfltSilentTrimAlignBasePts;
			if (SilentTrimAlignBasePts > AlignedBasePts)
				SilentTrimAlignBasePts = AlignedBasePts;
			else
				if (SilentTrimAlignBasePts < 0)				// must always be >= 0
					SilentTrimAlignBasePts = 0;

			MisalignedBasePts = misalignedbasepts->count ? misalignedbasepts->ival[0] : cBMDfltMisalignedBasePts;
			if (MisalignedBasePts < cBMMinPts)
				MisalignedBasePts = cBMMinPts;
			else
				if (MisalignedBasePts >= 0)
					MisalignedBasePts = -1;
			UnalignedBasePts = unalignedbasepts->count ? unalignedbasepts->ival[0] : cBMDfltUnalignedBasePts;
			if (UnalignedBasePts >= AlignedBasePts)
				MisalignedBasePts = AlignedBasePts - 1;
			else
				if (UnalignedBasePts <= MisalignedBasePts)
					MisalignedBasePts = MisalignedBasePts + 1;
			}

		if (PMode == eBMSimReads || PMode == eBMScore)
			{
			if (reads->count)
				{
				strcpy(szSEReads, reads->filename[0]);
				CUtility::TrimQuotedWhitespcExtd(szSEReads);
				}
			if (szSEReads[0] == '\0')
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Simulation or scoring requested but no %s reads file specified", bPEReads? "PE1" : "SE");
				exit(1);
				}

			if(PMode == eBMScore)
				{
				if (pe2reads->count)
					{
					strcpy(szPE2Reads, pe2reads->filename[0]);
					CUtility::TrimQuotedWhitespcExtd(szPE2Reads);
					}
				}
			if (bPEReads)
				{
				if (pe2reads->count)
					{
					strcpy(szPE2Reads, pe2reads->filename[0]);
					CUtility::TrimQuotedWhitespcExtd(szPE2Reads);
					}
				if (szPE2Reads[0] == '\0')
					{
					gDiagnostics.DiagOut(eDLFatal, gszProcName, "Simulation or scoring requested but no PE2 reads file specified");
					exit(1);
					}
				}
			else
				szPE2Reads[0] = '\0';
			}
		else
			{
			szSEReads[0] = '\0';
			szPE2Reads[0] = '\0';
			}

		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Processing parameters:");
		const char* pszDescr;
		switch (PMode) {
			case eBMGenCIGARs:		// generate observed CIGARs from alignments
				pszDescr = "Generate observed CIGARs from alignments";
				break;

			case eBMSimReads:	   // simulate reads using observed CIGARs
				pszDescr = "Simulate reads using observed CIGARs";
				break;

			case eBMScore:			// score base alignments using expected ground truth CIGARs
				pszDescr = "Score base alignments from simulated reads";
				break;
			}


		gDiagnostics.DiagOutMsgOnly(eDLInfo, "Aligner benchmark processing : '%s'", pszDescr);
		if (PMode == eBMGenCIGARs || PMode == eBMSimReads)
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "CIGAR file : '%s'", szCIGARsFile);
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "Aligner benchmark processing for : '%s'", bPEReads ? "PE" : "SE");
		if (PMode == eBMGenCIGARs || PMode == eBMScore)
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "Alignment file (SAM/BAM) : '%s'", szAlignmentsFile);
		if (PMode == eBMGenCIGARs || PMode == eBMSimReads)
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "Indexed reference genome file : '%s'", szRefGenomeFile);
		
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "Maximum number of %s to be processed : %d", bPEReads ? "PE read pairs" : "SE reads", MaxNumReads);
		if(PMode == eBMScore)
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "Scoring is for: %s", bPrimaryOnly ? "Primary alignments only" : "All alignments");
		if(szSEReads[0] != '\0')
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "%s file : '%s'", bPEReads ? "PE1" : "SE", szSEReads);
		if (szPE2Reads[0] != '\0')
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "PE2 file : '%s'", szPE2Reads);


		if(PMode == eBMScore)
			{
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "Aligned base score points: %d", AlignedBasePts);
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "Aligned silently trimmed read base score points: %d", SilentTrimAlignBasePts);
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "Misaligned base score points: %d", MisalignedBasePts);
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "Unaligned base score points: %d", UnalignedBasePts);
			}

		gDiagnostics.DiagOutMsgOnly(eDLInfo, "Experiment name : '%s'", szExperimentName);
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "Experiment description : '%s'", szExperimentDescr);

		if (szExperimentName[0] != '\0')
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "This processing reference: %s", szExperimentName);

		if (gExperimentID > 0)
			{
			int ParamID;
			ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID, ePTText, (int)strlen(szLogFile), "log", szLogFile);
			ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID, ePTInt32, (int)sizeof(PMode), "mode", &PMode);
			ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID, ePTBool, (int)sizeof(bPEReads), "pereads", &bPEReads);
			ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID, ePTBool, (int)sizeof(bPrimaryOnly), "pereads", &bPrimaryOnly);
			ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID, ePTInt32, (int)sizeof(NumberOfProcessors), "cpus", &NumberOfProcessors);
			ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID, ePTText, (int)strlen(szSQLiteDatabase), "sumrslts", szSQLiteDatabase);
			ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID, ePTText, (int)strlen(szExperimentName), "experimentname", szExperimentName);
			ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID, ePTText, (int)strlen(szExperimentDescr), "experimentdescr", szExperimentDescr);
		}

#ifdef _WIN32
		SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
		gStopWatch.Start();
		Rslt = 0;
		Rslt = BenchmarkProcess(PMode,					// processing mode
								bPrimaryOnly,			// if true then score only primary alignments otherwise score including secondary
							    bPEReads,				// if true then PE pair processing otherwise SE reads
								UnalignedBasePts,		// points to apply for each unaligned base
								AlignedBasePts,			// points to apply for each base aligned to it's ground truth loci
								SilentTrimAlignBasePts, // points to apply for each base aligned in a silently trimmed read within ground truth loci range
								MisalignedBasePts,		// points to apply for each base aligned but not to ground truth loci
								MaxNumReads,			// maximum number of alignment CIGARs to process or number of simulated reads 
								szCIGARsFile,			// observed CIGARs are in this file
								szRefGenomeFile,		// reads are against this target genome file
								szAlignmentsFile,		// input file containing aligned reads (SAM or BAM)
								szSEReads,				// simulated reads are output to this file (SE or PE1 if PE)
								szPE2Reads);				// simulated PE2 reads are output to this file if simulating PE reads
		Rslt = Rslt >= 0 ? 0 : 1;
		if (gExperimentID > 0)
		{
			if (gProcessingID)
				gSQLiteSummaries.EndProcessing(gExperimentID, gProcessingID, Rslt);
			gSQLiteSummaries.EndExperiment(gExperimentID);
		}
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
	return 0;
}

TRandomCombined<CRandomMother, CRandomMersenne> RG((int)time(0));

int
BenchmarkProcess(eBMProcMode PMode,	// processing mode
	bool bPrimaryOnly,			// if true then score only primary alignments otherwise score including secondary
	bool bPEReads,				// true if PE pair reads
	int UnalignedBasePts,		// points to apply for each unaligned base
	int AlignedBasePts,			// points to apply for each base aligned to it's ground truth loci
	int SilentTrimAlignBasePts, // points to apply for each base aligned in a silently trimmed read within ground truth loci range
	int MisalignedBasePts,		// points to apply for each base aligned but not to ground truth loci
	int MaxNumReads,			// maximum number of alignment CIGARs to process or number of simulated reads 
	char* pszObsCIGARsFile,		// observed CIGARs are in this file
	char* pszRefGenomeFile,		// alignments are against this target genome file
	char* pszAlignmentsFile,	// input file containing aligned reads (SAM or BAM)
	char* pszSEReads,			// simulated reads are output to this file (SE or PE1 if PE)
	char* pszPE2Reads)			// simulated PE2 reads are output to this file if simulating PE reads
{
int Rslt;
CBenchmark *pBenchmark;
pBenchmark = new CBenchmark;
Rslt = 0;

switch (PMode) {
	case eBMGenCIGARs:
		Rslt = pBenchmark->GenObsCIGARs(bPEReads, MaxNumReads, pszObsCIGARsFile, pszRefGenomeFile, pszAlignmentsFile);
		break;

	case eBMSimReads:
		Rslt = pBenchmark->SimReads(bPEReads, MaxNumReads, pszObsCIGARsFile, pszRefGenomeFile, pszSEReads, pszPE2Reads);
		break;

	case eBMScore:
		Rslt = pBenchmark->Score(bPrimaryOnly,bPEReads, UnalignedBasePts, AlignedBasePts, SilentTrimAlignBasePts, MisalignedBasePts, pszSEReads, pszPE2Reads, pszAlignmentsFile);
		break;
	}


if(pBenchmark != NULL)
	delete pBenchmark;

return(Rslt);
}


CBenchmark::CBenchmark()						// instance constructor
{
m_pGenome = NULL;
m_pAlignments = NULL;
m_pObsErrProfiles = NULL;
m_pGroundTruths = NULL;
m_ppGroundTruthIdx = NULL;
m_pObsCIGARProfFile = NULL;
m_pObsCIGARBuff = NULL;
m_pChromSeqs = NULL;
m_pszSESimReadBuff = NULL;
m_pszPE2SimReadBuff = NULL;
m_pSESimReads = NULL;
m_pPE2SimReads = NULL;
m_hObsSIGARs = -1;
m_hSEReads = -1;
m_hPE2Reads = -1;
Reset();
}

CBenchmark::~CBenchmark()						// instance destructor
{
if(m_pObsCIGARProfFile != NULL)
	delete m_pObsCIGARProfFile;
if(m_pGenome != NULL)
	delete m_pGenome;
if (m_pAlignments != NULL)
	delete m_pAlignments;
if (m_pObsErrProfiles != NULL)
	{
#ifdef _WIN32
	free(m_pObsErrProfiles);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (m_pObsErrProfiles != MAP_FAILED)
		munmap(m_pObsErrProfiles, m_AllocdObsErrProfMem);
#endif
	}

if(m_pGroundTruths != NULL)
	{
#ifdef _WIN32
	free(m_pGroundTruths);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (m_pGroundTruths != MAP_FAILED)
		munmap(m_pGroundTruths, m_AllocdGroundTruthsMem);
#endif
	}

if(m_ppGroundTruthIdx != NULL)
	delete m_ppGroundTruthIdx;

if(m_pObsCIGARBuff != NULL)
	delete m_pObsCIGARBuff;

if(m_pChromSeqs != NULL)
	delete m_pChromSeqs;

if (m_pszSESimReadBuff != NULL)
	delete m_pszSESimReadBuff;

if (m_pszPE2SimReadBuff != NULL)
	delete m_pszPE2SimReadBuff;

if(m_pSESimReads != NULL)
	delete m_pSESimReads;

if(m_pPE2SimReads != NULL)
	delete m_pPE2SimReads;

if (m_hObsSIGARs != -1)
	close(m_hObsSIGARs);

if (m_hSEReads != -1)
	close(m_hSEReads);

if (m_hPE2Reads != -1)
	close(m_hPE2Reads);
}

void
CBenchmark::Reset(void)
{
if (m_pObsCIGARProfFile != NULL)
	{
	delete m_pObsCIGARProfFile;
	m_pObsCIGARProfFile = NULL;
	}

if (m_pGenome != NULL)
	{
	delete m_pGenome;
	m_pGenome = NULL;
	}

if (m_pAlignments != NULL)
	{
	delete m_pAlignments;
	m_pAlignments = NULL;
	}

if (m_pObsErrProfiles != NULL)
	{
#ifdef _WIN32
	free(m_pObsErrProfiles);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (m_pObsErrProfiles != MAP_FAILED)
		munmap(m_pObsErrProfiles, m_AllocdObsErrProfMem);
#endif
	m_pObsErrProfiles = NULL;
	}

if (m_pGroundTruths != NULL)
	{
#ifdef _WIN32
	free(m_pGroundTruths);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (m_pGroundTruths != MAP_FAILED)
		munmap(m_pGroundTruths, m_AllocdGroundTruthsMem);
#endif
	m_pGroundTruths = NULL;
	}

if (m_ppGroundTruthIdx != NULL)
	{
	delete m_ppGroundTruthIdx;
	m_ppGroundTruthIdx = NULL;
	}

if (m_hObsSIGARs != -1)
	{
#ifdef _WIN32
	_commit(m_hObsSIGARs);
#else
	fsync(m_hObsSIGARs);
#endif
	close(m_hObsSIGARs);
	m_hObsSIGARs = -1;
	}

if (m_hSEReads != -1)
	{
#ifdef _WIN32
	_commit(m_hSEReads);
#else
	fsync(m_hSEReads);
#endif
	close(m_hSEReads);
	m_hSEReads = -1;
	}

if (m_hPE2Reads != -1)
	{
#ifdef _WIN32
	_commit(m_hPE2Reads);
#else
	fsync(m_hPE2Reads);
#endif
	close(m_hPE2Reads);
	m_hPE2Reads = -1;
	}

if (m_pObsCIGARBuff != NULL)
	{
	delete m_pObsCIGARBuff;
	m_pObsCIGARBuff = NULL;
	}

if(m_pChromSeqs != NULL)
	{
	delete m_pChromSeqs;
	m_pChromSeqs = NULL;
	}

if (m_pszSESimReadBuff != NULL)
	{
	delete m_pszSESimReadBuff;
	m_pszSESimReadBuff = NULL;
	}

if(m_pszPE2SimReadBuff != NULL)
	{
	delete m_pszPE2SimReadBuff;
	m_pszPE2SimReadBuff = NULL;
	}

if (m_pSESimReads != NULL)
	{
	delete m_pSESimReads;
	m_pSESimReads = NULL;
	}

if (m_pPE2SimReads != NULL)
	{
	delete m_pPE2SimReads;
	m_pPE2SimReads = NULL;
	}

m_bPrimaryOnly = false;
m_bPEReads = false;				// if true then PE pairs processing otherwise SE reads
m_UnalignedBasePts = 0;			// points to apply for each unaligned base
m_AlignedBasePts = 0;			// points to apply for each base aligned to it's ground truth loci
m_SilentTrimAlignBasePts = 0;	// points to apply for each base aligned in a silently trimmed read within ground truth loci range
m_MisalignedBasePts = 0;		// points to apply for each base aligned but not to ground truth loci
m_ObsCIGARBuffLen = 0;
m_AllocdObsCIGARBuffSize = 0;
m_AllocdObsErrProfMem = 0;
m_AllocdObsCIGARBuffSize = 0;

m_TotNumPotentialAlignBases = 0;
m_TotGroundTruthBases = 0;	// total number of potential ground truth bases which could be matched in an alignment 
m_NumGroundTruths = 0;		// loaded this many ground truths loaded from simulated reads
m_UsedGroundTruthsMem = 0;	// loaded ground truths are using this sized memory
m_AllocdGroundTruthsMem = 0;// allocated this size memory to hold observed error profiles
m_UnscoredReads = 0;
m_ScoredReads = 0;

m_NumObsErrProfs = 0;		// loaded this many observed error profiles
m_UsedObsErrProfMem = 0;	// loaded observed error profiles using this sized memory
m_AllocdObsErrProfMem = 0;	// allocated this size memory to hold observed error profiles

m_GenomeLen = 0;
m_GenomeScaledLen = 0;
m_NumChromSeqs = 0;

m_MaxNumReads = 0;				// maximum number of alignment CIGARs to process or number of simulated reads or read pairs 

m_szObsCIGARsFile[0] = '\0';		// observed CIGARs are in this file
m_szGroundTruthFile[0] = '\0';	// simulated reads ground truth (CIGARs and loci for each simulated read)) are in this file
m_szRefGenomeFile[0] = '\0';	// reads are against this target genome file
m_szAlignmentsFile[0] = '\0';	// input file containing aligned reads (SAM or BAM)
m_szSEReads[0] = '\0';			// simulated reads are output to this file (SE or PE1 if PE)
m_szPE2Reads[0] = '\0';			// simulated PE2 reads are output to this file if simulating PE reads

m_AllocOutBuffSize = 0;			// allocated simulated read buffers of this size
m_OutSEBuffIdx = 0;				// currently using this number of chars in pszSESimReadBuff
m_OutPE2BuffIdx = 0;			// currently using this number of chars in pszSESimReadBuff
m_szTargSpecies[0] = '\0';
}

// try loading the indexed genome
int 
CBenchmark::LoadGenome(char* pszRefGenomeFile)		// read alignments were against this targeted genome file
{
int Rslt;
size_t TotLen;
int ChromID;
tsBMChromSeq* pChromSeq;
double LenSCF;			// length scaling factor
int CurScaledStart;
tsSfxHeaderV3 SfxHeader;

if (pszRefGenomeFile == NULL || *pszRefGenomeFile == '\0')
	{
	Reset();
	return(eBSFerrParams);
	}

if (m_pGenome != NULL)
{
	delete m_pGenome;
	m_pGenome = NULL;
}
if ((m_pGenome = new CSfxArray()) == NULL)
{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to instantiate CSfxArray");
	Reset();
	return(eBSFerrObj);
}
if ((Rslt = m_pGenome->Open(pszRefGenomeFile, false, false, false)) != eBSFSuccess)
{
	while (m_pGenome->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal, gszProcName, m_pGenome->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to open input genome suffix array file '%s'", pszRefGenomeFile);
	Reset();
	return(Rslt);
}

// report to user some sfx array metadata as conformation the targeted assembly is correct
strcpy(m_szTargSpecies, m_pGenome->GetDatasetName());

m_pGenome->GetSfxHeader(&SfxHeader);
gDiagnostics.DiagOut(eDLInfo, gszProcName, "Genome Assembly Name: '%s' Descr: '%s' Title: '%s' Version: %d",
	m_szTargSpecies, SfxHeader.szDescription, SfxHeader.szTitle, SfxHeader.Version);

gDiagnostics.DiagOut(eDLInfo, gszProcName, "Loading genome assembly suffix array...");
if ((Rslt = m_pGenome->SetTargBlock(1)) < 0)
{
	while (m_pGenome->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal, gszProcName, m_pGenome->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Fatal: unable to load genome assembly suffix array");
	return(Rslt);
}

TotLen = m_pGenome->GetTotSeqsLen();
// assembly or multifasta sequences now loaded
LenSCF = (double)INT_MAX / (double)TotLen;
m_GenomeScaledLen = TotLen >= (INT64)INT_MAX ? (long)(TotLen * LenSCF) : (long)TotLen;
m_NumChromSeqs = m_pGenome->GetNumEntries();
m_pChromSeqs = new tsBMChromSeq [m_NumChromSeqs];

CurScaledStart = 0;
pChromSeq = m_pChromSeqs;
for (ChromID = 1; ChromID <= m_NumChromSeqs; ChromID++, pChromSeq++)
	{
	pChromSeq->ChromID = ChromID;
	pChromSeq->Len = m_pGenome->GetSeqLen(ChromID);
	pChromSeq->ScaledLen = TotLen >= (INT64)INT_MAX ? (int)(pChromSeq->Len * LenSCF) : pChromSeq->Len;
	pChromSeq->ScaledStart = CurScaledStart;
	CurScaledStart += pChromSeq->ScaledLen;
	}
gDiagnostics.DiagOut(eDLInfo, gszProcName, "Genome assembly suffix array loaded");
return(eBSFSuccess);
}

int
CBenchmark::OpenAlignments(char* pszAlignmentsFile)	// input file containing aligned reads (SAM or BAM)
{
	teBSFrsltCodes Rslt;

	// open SAM/BAM for reading
	if (pszAlignmentsFile == NULL || *pszAlignmentsFile == '\0')
	{
		Reset();
		return(eBSFerrParams);
	}

	if (m_pAlignments != NULL)
	{
		delete m_pAlignments;
		m_pAlignments = NULL;
	}

	if ((m_pAlignments = new CSAMfile) == NULL)
	{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to instantiate class CSAMfile");
		Reset();
		return(eBSFerrInternal);
	}

	if ((Rslt = (teBSFrsltCodes)m_pAlignments->Open(pszAlignmentsFile)) != eBSFSuccess)
	{
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Unable to open SAM/BAM format file %s", pszAlignmentsFile);
		delete m_pAlignments;
		m_pAlignments = NULL;
		Reset();
		return((teBSFrsltCodes)Rslt);
	}
	return(eBSFSuccess);
}



int
CBenchmark::GenObsCIGARs(bool bPEReads,	// true if PE pair processing otherwise SE reads
	int MaxNumReads,		// maximum number of aligned reads or read pair alignments to accept for putative error profile processing
	char* pszObsCIGARsFile,	// write observed error profile CIGARs to this file
	char* pszRefGenomeFile,	// alignments were against this target genome file
	char* pszAlignmentsFile)// input file containing alignments of simulated reads (SAM or BAM)
{
int Rslt;
int NumPutativeAlignments;
int NumAcceptedAlignments;
char* pszLine;
char* pTxt;
int LineLen;
int ReadOfs;
int AlignLoci;
etSeqBase* pReadSeq;
etSeqBase* pTargSeq;
etSeqBase ReadBase;
etSeqBase TargBase;
int CurOpIdx;
int NumMatches;
int NumMismatches;
int CigarIdx;
int NumOps;
int AlignSeqLen;
int MaxProfileSeqLen;
int ProfileSeqLen;
etCIGAROpType TypeOp;
int NumTypeOps;
bool bSloughAlignment;
char szErrProfile[cSRMaxCIGARLen +1 ];
char szObsCIGAR[cSRMaxCIGARLen + 1];
int ObsCIGARLen;
tsBMObsCIGAR* pObsCIGAR;
tsBAMalign *pSAMalign;

Reset();
m_bPEReads = bPEReads;
m_MaxNumReads = MaxNumReads;
strcpy(m_szAlignmentsFile, pszAlignmentsFile);
strcpy(m_szObsCIGARsFile, pszObsCIGARsFile);
strcpy(m_szRefGenomeFile, pszRefGenomeFile);

m_UsedObsErrProfMem = 0;
m_NumObsErrProfs = 0;

m_AllocdObsErrProfMem = (size_t)cSRAllocdObsErrProfMemChunk;
#ifdef _WIN32
m_pObsErrProfiles = (UINT8*)malloc(m_AllocdObsErrProfMem);
if (m_pObsErrProfiles == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenObsCIGARs: Memory allocation of %lld bytes failed", (INT64)m_AllocdObsErrProfMem);
	m_AllocdObsErrProfMem = 0;
	Reset();
	return(eBSFerrMem);
	}
#else
// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
m_pObsErrProfiles = (UINT8*)mmap(NULL, (size_t)m_AllocdObsErrProfMem, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
if (m_pObsErrProfiles == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenObsCIGARs: Memory allocation of %lld bytes through mmap()  failed", (INT64)m_AllocdObsErrProfMem, strerror(errno));
	m_pObsErrProfiles = NULL;
	m_AllocdObsErrProfMem = 0;
	Reset();
	return(eBSFerrMem);
	}
#endif

gDiagnostics.DiagOut(eDLInfo, gszProcName, "Creating emperical error profiles file", pszObsCIGARsFile);
#ifdef _WIN32
if ((m_hObsSIGARs = open(pszObsCIGARsFile, (_O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC), (_S_IREAD | _S_IWRITE))) == -1)
#else
if ((m_hObsSIGARs = open(pszObsCIGARsFile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE)) == -1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to create output file - '%s' - %s",pszObsCIGARsFile, strerror(errno));
	Reset();
	return(eBSFerrOpnFile);
	}
gDiagnostics.DiagOut(eDLInfo, gszProcName, "Opening SAM/BAM alignments file %s from which emperical error profiles will be derived", pszAlignmentsFile);
if ((Rslt = OpenAlignments(m_szAlignmentsFile)) < 0)
	{
	Reset();
	return(Rslt);
	}

m_pObsCIGARBuff = new char [cBMCIGARAllocBuffSize];
m_AllocdObsCIGARBuffSize = cBMCIGARAllocBuffSize;
if(bPEReads)
	m_ObsCIGARBuffLen = sprintf(m_pObsCIGARBuff,"\"ID\",\"SeqLen\",\"PE1 Strand\",\"PE1 CIGAR\",\"PE1 Error Profile\",\"PE Insert Size\",\"PE2 Strand\",\"PE2 CIGAR\",\"PE2 Error Profile\"\n");
else
	m_ObsCIGARBuffLen = sprintf(m_pObsCIGARBuff, "\"ID\",\"SeqLen\",\"Strand\",\"CIGAR\",\"Error Profile\"\n");
CUtility::SafeWrite(m_hObsSIGARs, m_pObsCIGARBuff, m_ObsCIGARBuffLen);
m_ObsCIGARBuffLen = 0;

gDiagnostics.DiagOut(eDLInfo, gszProcName, "Opening alignment targeted assembly/transcriptome indexed file %s for processing", m_szRefGenomeFile);
if((Rslt = LoadGenome(m_szRefGenomeFile)) < 0)
	{
	Reset();
	return(Rslt);
	}
gDiagnostics.DiagOut(eDLInfo, gszProcName, "%s loaded, now processing alignments from %s", m_szRefGenomeFile, pszAlignmentsFile);


etSeqBase ReadSequence[1000];
etSeqBase TargetSequence[1000];

if ((pszLine = new char[cMaxReadLen]) == NULL)				// buffer input lines
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to allocate line buffering");
	Reset();
	return(-1);
	}

if ((pSAMalign = new tsBAMalign) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to allocate memory for alignment tsBAMalign structure");
	delete pszLine;
	Reset();
	return(-1);
	}

// now iterate over all alignments and extract the CIGARs for each alignment
int NumShortReads = 0;
int NumUnableToPair = 0;
int NumCIGARErrs = 0;
int NumCIGARLenErrs = 0;
int NumNonCvtdSoftClips = 0;
int NumUnknownRefs = 0;
int NumFragSizeBounds = 0;
int NumReadLenBounds = 0;
int NumSecAlignments = 0;
int NumNoMateReads = 0;
int NumCIGARSUnknown = 0;
int NumUnmapped = 0;
int NumMissingFeatures = 0;
int NumPEInconsistentStrands = 0;

NumPutativeAlignments = 0;
NumAcceptedAlignments = 0;
MaxProfileSeqLen = cBMMinReadLen;
time_t Then = time(NULL);
time_t Now;
while (Rslt >= eBSFSuccess && (LineLen = m_pAlignments->GetNxtSAMline(pszLine)) > 0)
{
	pszLine[cMaxReadLen - 1] = '\0';
	pTxt = CUtility::TrimWhitespc(pszLine);
	if (*pTxt == '\0')			// simply slough lines which are just whitespace
		continue;
	if (*pTxt == '@')				// only interested in lines with putative alignments
		continue;

	NumPutativeAlignments += 1;
	if (!(NumPutativeAlignments % 100000) || NumPutativeAlignments == 1)
		{
		Now = time(NULL);
		if ((Now - Then) >= 60)
			{
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "Processing %d SAM/BAM CIGAR alignments, putatively accepted %d for use as read simulation error profiles", NumPutativeAlignments, NumAcceptedAlignments);
			Then += 60;
			}
		}

	// primary interest is in the reference chromname, startloci, length
	if ((Rslt = (teBSFrsltCodes)m_pAlignments->ParseSAM2BAMalign(pTxt, pSAMalign, NULL,true)) < eBSFSuccess)
		{
		if (Rslt == eBSFerrFeature)	// not too worried if aligned to feature is missing as some SAMs are missing header features
			{
			NumMissingFeatures++;
			Rslt = eBSFSuccess;
			continue;
			}
		break;
		}

	// check if read has been mapped, if not then slough ...
	if (pSAMalign->refID == -1 || (pSAMalign->flag_nc >> 16) & 0x04)
		{
		NumUnmapped++;
		continue;
		}
	if(pSAMalign->cigar[0] == '*')	// is Cigar is known
		{
		NumCIGARSUnknown++;
		continue;
		}
	if (bPEReads && ((pSAMalign->flag_nc >> 16) & 0x03) != 0x03)// if PE processing then read has to have a mate read
		{
		NumNoMateReads++;
		continue;
		}
	if ((pSAMalign->flag_nc >> 16) & 0x100)						// not interested in any secondary alignments
		{
		NumSecAlignments++;
		continue;
		}


	if(pSAMalign->l_seq > cBMMaxReadLen || pSAMalign->l_seq < cBMMinReadLen)
		{
		NumReadLenBounds++;
		continue;
		}

	// note: tlen can actually be negative (antisense fragment) and excludes soft clipped bases from start/end of fragment
	// further note that if SE alignments, tlen can be 0
	if (abs(pSAMalign->tlen) > cBMMaxFragSize)
		{
		NumFragSizeBounds;
		continue;
		}

	// ensure that the reference is known, if not then silently discard
	if ((pSAMalign->refID = m_pGenome->GetIdent(pSAMalign->szRefSeqName)) < 1)
		{
		NumUnknownRefs++;
		continue;
		}
//Op BAM Description                                          Consumes query    Consumes reference 
//M  0   alignment match(can be a sequence match or mismatch) yes               yes 
//I  1   insertion to the reference                           yes               no 
//D  2   deletion from the reference                          no                yes 
//N  3   skipped region from the reference                    no                yes 
//S  4   soft clipping(clipped sequences present in SEQ)      yes               no 
//H  5   hard clipping(clipped sequences NOT present in SEQ)  no                no 
//P  6   padding(silent deletion from padded reference)       no                no 
//=  7   sequence match                                       yes               yes
//X  8   sequence mismatch                                    yes               yes

	// Not accepting hard trimmed alignments as trimmed off bases are unknown, if soft trimmed then bases are known.
	// Ensuring that CIGAR has no internal soft or hard trimming and does not start/end with eCOPInsert, eCOPDelete,eCOPSkipRegion, or eCOPPadding

	if(pSAMalign->szMateRefSeqName[0] == '=' || !stricmp(pSAMalign->szRefSeqName, pSAMalign->szMateRefSeqName))
		pSAMalign->next_refID = pSAMalign->refID;

	if((NumTypeOps=CoalesceCIGARs(pSAMalign)) <= 0)
		{
		NumCIGARErrs++;
		continue;
		}

	// check if alignment only starts with match or soft trimming, and if soft trimming that next TypeOp is a match
	TypeOp = (etCIGAROpType)(pSAMalign->cigar[0] & 0x0f);
	if(!(TypeOp == eCOPMatch || TypeOp == eCOPSoftClip))
		{
		NumCIGARErrs++;
		continue;
		}
	if(TypeOp == eCOPSoftClip) //	'S'  soft clipping - consumes query sequence only
		{
		if(NumTypeOps == 1)
			{
			NumCIGARErrs++;
			continue;
			}
		TypeOp = (etCIGAROpType)(pSAMalign->cigar[1] & 0x0f);
		if(TypeOp != eCOPMatch)
			{
			NumCIGARErrs++;
			continue;
			}
		
		}

	// check if alignment only ends with match or soft trimming, and if soft trimming 3' that previous TypeOp is a match
	if(NumTypeOps > 1)
		{
		TypeOp = (etCIGAROpType)(pSAMalign->cigar[NumTypeOps-1] & 0x0f);
		if (!(TypeOp == eCOPMatch || TypeOp == eCOPSoftClip))
			{
			NumCIGARErrs++;
			continue;
			}
		if (NumTypeOps == 2 && TypeOp == eCOPSoftClip)
			{
			TypeOp = (etCIGAROpType)(pSAMalign->cigar[0] & 0x0f);
			if (TypeOp != eCOPMatch)
				{
				NumCIGARErrs++;
				continue;
				}
			}
		}

	// ensure that CIGAR contains no internal soft or hard trimming
	if (NumTypeOps >= 3)
		{
		int TypeOpsIdx;
		for (TypeOpsIdx = 1; TypeOpsIdx < (NumTypeOps - 1); TypeOpsIdx++)
			{
			TypeOp = (etCIGAROpType)(pSAMalign->cigar[TypeOpsIdx] & 0x0f);
			if (TypeOp == eCOPHardClip || TypeOp == eCOPSoftClip)
				break;
			}
		if (TypeOpsIdx != NumTypeOps - 1)
			{
			NumCIGARErrs++;
			continue;
			}
		}

	// only accepted alignments for which full sequence is known - sloughed out any hard trimmed alignments	
	AlignSeqLen = m_pAlignments->BAMalignSeq(pSAMalign,	// ptr to tsBAMalign alignment containing packed (2 per byte) sequence to be returned as etSeqBases 
							cBMMaxReadLen,				// maximum length sequence to be returned
							ReadSequence);				// to hold returned sequence
	if(AlignSeqLen != pSAMalign->l_seq) 	// pays to double check :-))
		{
		NumCIGARLenErrs++;
		continue;
		}

	if(AlignSeqLen < MaxProfileSeqLen)
		{
		NumShortReads++;
		continue;
		}

	// read CIGAR is packed, decode into an ASCII string
	ObsCIGARLen = DecodeCIGAR(pSAMalign->flag_nc & 0x0ff, pSAMalign->cigar,sizeof(szObsCIGAR),szObsCIGAR);
	// take read and using the alignment CIGAR as a guide then regenerate a CIGAR which states which bases are exact matches and which mismatch the target sequence
	ReadOfs = 0;
	CigarIdx = 0;
	CurOpIdx = 1;
	NumMatches = 0;
	NumMismatches = 0;
	ProfileSeqLen = 0;
	AlignLoci = pSAMalign->pos;
	szErrProfile[0] = '\0';
	bSloughAlignment = false;
	while ((CurOpIdx = GetClaimOps(pSAMalign, CurOpIdx, &NumOps, &TypeOp)) > 0)
		{
		if(TypeOp == eCOPSoftClip)			// if soft clipping then note that the start ref loci needs to be adjusted back to read start!!!
			{
			if(CurOpIdx == 2)
				{
				AlignLoci -= NumOps;
				if(AlignLoci < 0)			// can't handle soft clipping where ref loci would be negative!
					{
					NumNonCvtdSoftClips++;
					bSloughAlignment = true;
					break;
					}
				}
			TypeOp = eCOPMatch;
			}
		NumOps = AdjustAlignment(TypeOp, NumOps, &ReadOfs, &AlignLoci);
		if(NumOps <= 0)
			{
			NumCIGARLenErrs++;
			bSloughAlignment = true;
			break;
			}

		if(TypeOp == eCOPInsert || TypeOp == eCOPDelete || TypeOp == eCOPSkipRegion || TypeOp == eCOPPadding)
			{
			if (NumMatches > 0)
				{
				ProfileSeqLen += NumMatches;
				CigarIdx += sprintf(&szErrProfile[CigarIdx], "%d=", NumMatches);
				NumMatches = 0;
				}
			else
				if (NumMismatches > 0)
					{
					ProfileSeqLen += NumMismatches;
					CigarIdx += sprintf(&szErrProfile[CigarIdx], "%dX", NumMismatches);
					NumMismatches = 0;
					}
			switch (TypeOp) {
				case eCOPInsert:		//	'I'  insertion relative to target - consumes query sequence only
					ProfileSeqLen += NumOps;
					CigarIdx += sprintf(&szErrProfile[CigarIdx], "%dI", NumOps);
					continue;

				case eCOPDelete:	// 'D'  deletion relative to target - consumes reference sequence only
					CigarIdx += sprintf(&szErrProfile[CigarIdx], "%dD", NumOps);
					continue;

				case eCOPSkipRegion:	//	'N'  skipped region relative to target - intron?  - consumes reference sequence only
					CigarIdx += sprintf(&szErrProfile[CigarIdx], "%dN", NumOps);
					continue;

				case eCOPPadding:		//	'P' padding (silent deletion from padded reference) consumes neither query or reference sequence
					CigarIdx += sprintf(&szErrProfile[CigarIdx], "%dP", NumOps);
					continue;
				}
			}

// can now check read bases against target sequence bases
		pReadSeq = &ReadSequence[ReadOfs];
		m_pGenome->GetSeq(pSAMalign->refID, AlignLoci, TargetSequence, NumOps);
		ReadOfs += NumOps;
		AlignLoci += NumOps;
		pTargSeq = TargetSequence;
		do {
			ReadBase = *pReadSeq++;
			TargBase = *pTargSeq++;
			if((ReadBase & 0x0f) == (TargBase & 0x0f))
				{ // matches
				if (NumMismatches > 0)
					{
					ProfileSeqLen += NumMismatches;
					CigarIdx += sprintf(&szErrProfile[CigarIdx], "%dX", NumMismatches);
					NumMismatches = 0;
					}
				NumMatches++;
				}
			else
				{ // mismatches
				if (NumMatches > 0)
					{
					ProfileSeqLen += NumMatches;
					CigarIdx += sprintf(&szErrProfile[CigarIdx], "%d=", NumMatches);
					NumMatches = 0;
					}
				NumMismatches++;
				}
			}
		while(--NumOps);
		}
	if(bSloughAlignment)
		continue;
	if (NumMatches > 0)
		{
		ProfileSeqLen += NumMatches;
		CigarIdx += sprintf(&szErrProfile[CigarIdx], "%d=", NumMatches);
		NumMatches = 0;
		}
	else
		if (NumMismatches > 0)
			{
			ProfileSeqLen += NumMismatches;
			CigarIdx += sprintf(&szErrProfile[CigarIdx], "%dX", NumMismatches);
			NumMismatches = 0;
			}

	if(ProfileSeqLen > MaxProfileSeqLen)
		MaxProfileSeqLen = ProfileSeqLen;
	else
		if(ProfileSeqLen < MaxProfileSeqLen) // some aligners silently hard trim so need to slough these, current seq len will be shorter than longest observed
			{
			NumShortReads++;
			continue;
			}
	
	if (((sizeof(tsBMObsCIGAR) + cSRAllocdObsErrProfMemChunk/10) >= (m_AllocdObsErrProfMem - m_UsedObsErrProfMem)))
		{
		size_t memreq = m_AllocdObsErrProfMem + (size_t)cSRAllocdObsErrProfMemChunk;
		UINT8* pTmp;
#ifdef _WIN32
		pTmp = (UINT8*)realloc(m_pObsErrProfiles, memreq);
#else
		pTmp = (UINT8*)mremap(m_pObsErrProfiles, m_AllocdObsErrProfMem, memreq, MREMAP_MAYMOVE);
		if (pTmp == MAP_FAILED)
			pTmp = NULL;
#endif
		if (pTmp == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenObsCIGARs: Memory re-allocation to %lld bytes - %s", (INT64)(memreq), strerror(errno));
			delete pszLine;
			delete pSAMalign;
			Reset();
			return(eBSFerrMem);
			}
		m_pObsErrProfiles = pTmp;
		m_AllocdObsErrProfMem = memreq;
		}

	pObsCIGAR = (tsBMObsCIGAR*)&m_pObsErrProfiles[m_UsedObsErrProfMem];
	pObsCIGAR->ReadLen = ProfileSeqLen;
	pObsCIGAR->InsertSize = abs(pSAMalign->tlen);
	pObsCIGAR->Flags = (pSAMalign->flag_nc >> 16) & 0x0000ffff;
	strncpy((char *)&pObsCIGAR->NameCIGARErrProfile, pSAMalign->read_name, pSAMalign->NumReadNameBytes);
	if (pSAMalign->NumReadNameBytes > 3 && 
			pSAMalign->read_name[pSAMalign->NumReadNameBytes - 3] == '/' && (pSAMalign->read_name[pSAMalign->NumReadNameBytes - 2] == '1' || pSAMalign->read_name[pSAMalign->NumReadNameBytes - 2] == '2'))
		{
		pObsCIGAR->NameLen = pSAMalign->NumReadNameBytes - 3;
		pObsCIGAR->NameCIGARErrProfile[pObsCIGAR->NameLen-1] = '\0';
		}
	else
		pObsCIGAR->NameLen = pSAMalign->NumReadNameBytes;
	pObsCIGAR->ErrProfileLen = CigarIdx + 1;
	pObsCIGAR->ObsCIGARlen = ObsCIGARLen + 1;
	strcpy((char*)&pObsCIGAR->NameCIGARErrProfile[pObsCIGAR->NameLen], szObsCIGAR);
	strcpy((char*)&pObsCIGAR->NameCIGARErrProfile[pObsCIGAR->NameLen + pObsCIGAR->ObsCIGARlen], szErrProfile);
	pObsCIGAR->Size = (uint16_t)sizeof(tsBMObsCIGAR) + pObsCIGAR->NameLen + pObsCIGAR->ObsCIGARlen + pObsCIGAR->ErrProfileLen - 1;
	m_UsedObsErrProfMem += pObsCIGAR->Size;
	NumAcceptedAlignments += 1;
	if(NumAcceptedAlignments >= (bPEReads ? MaxNumReads * 2 : MaxNumReads))
		break;
	}
if (NumAcceptedAlignments < (bPEReads ? cBMMinReads * 2 : cBMMinReads))
	{
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Processing of %d SAM/BAM CIGAR alignments completed, only able to putatively accept %d which is fewer than minimum of %d required", NumPutativeAlignments, NumAcceptedAlignments, cBMMinReads);
	delete pszLine;
	delete pSAMalign;
	Reset();
	return(-1);
	}

gDiagnostics.DiagOut(eDLInfo, gszProcName, "Processing of %d SAM/BAM CIGAR alignments completed, putatively accepted %d for use as read simulation error profiles", NumPutativeAlignments, NumAcceptedAlignments);

if(bPEReads)
	{
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Initialising for identification of PE pair association by read names");
	pObsCIGAR = (tsBMObsCIGAR*)m_pObsErrProfiles;
	tsBMObsCIGAR**pSortAssoc = new tsBMObsCIGAR *[NumAcceptedAlignments];
	for(int Idx =0; Idx < NumAcceptedAlignments; Idx++)
		{
		pSortAssoc[Idx] = pObsCIGAR;
		pObsCIGAR = (tsBMObsCIGAR*)((uint8_t *)pObsCIGAR + pObsCIGAR->Size);
		}
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Sorting alignments by read name for PE pair association");
	qsort(pSortAssoc, NumAcceptedAlignments, sizeof(tsBMObsCIGAR*), SortReadPairs);
	// now iterate over read names and only accept those PE pairs where both reads are present

	tsBMObsCIGAR*pNxtObsCIGAR;
	int NumPairs = 0;
	for (int Idx = 0; Idx < NumAcceptedAlignments-1; Idx++)
		{
		pObsCIGAR = pSortAssoc[Idx];
		pNxtObsCIGAR = pSortAssoc[Idx+1];
		if(pObsCIGAR->ReadLen < MaxProfileSeqLen ||	// some may have slipped through before it was known what read lengths were in original sequenced reads prior to any alignment trimming
			pNxtObsCIGAR->ReadLen < MaxProfileSeqLen)
			{
			NumShortReads++;
			continue;
			}
		if(stricmp((char *)pObsCIGAR->NameCIGARErrProfile, (char*)pNxtObsCIGAR->NameCIGARErrProfile))
			{
			NumUnableToPair++;
			continue;
			}
		// have a pair, check on strands, can't be same if PE
		if((pObsCIGAR->Flags & 0x010) == (pNxtObsCIGAR->Flags & 0x010))
			{
			NumPEInconsistentStrands++;
			continue;
			}

		NumPairs += 1;
		m_ObsCIGARBuffLen += sprintf(&m_pObsCIGARBuff[m_ObsCIGARBuffLen],"%d,%d,%c,%s,%s,%d,%c,%s,%s\n", 
			NumPairs, pObsCIGAR->ReadLen,
			(pObsCIGAR->Flags & 0x10) ? '-' : '+',
			&pObsCIGAR->NameCIGARErrProfile[pObsCIGAR->NameLen], 
			&pObsCIGAR->NameCIGARErrProfile[pObsCIGAR->NameLen + pObsCIGAR->ObsCIGARlen],
			pObsCIGAR->InsertSize, (pNxtObsCIGAR->Flags & 0x10) ? '-' : '+', 
			&pNxtObsCIGAR->NameCIGARErrProfile[pNxtObsCIGAR->NameLen], 
			&pNxtObsCIGAR->NameCIGARErrProfile[pNxtObsCIGAR->NameLen + pNxtObsCIGAR->ObsCIGARlen]);
		if ((m_ObsCIGARBuffLen + 1000) >= m_AllocdObsCIGARBuffSize)
			{
			CUtility::SafeWrite(m_hObsSIGARs, m_pObsCIGARBuff, m_ObsCIGARBuffLen);
			m_ObsCIGARBuffLen = 0;
			}
		Idx += 1;
		}
	if (m_ObsCIGARBuffLen)
		{
		CUtility::SafeWrite(m_hObsSIGARs, m_pObsCIGARBuff, m_ObsCIGARBuffLen);
		m_ObsCIGARBuffLen = 0;
		}
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Reported %d PE read pair error profiles to '%s' with read lengths of %dbp", NumPairs, pszObsCIGARsFile, MaxProfileSeqLen);
	delete []pSortAssoc;
	}
else
	{ // SE read error profiles
	pObsCIGAR = (tsBMObsCIGAR*)m_pObsErrProfiles;
	int NumSEs = 0;
	for (int Idx = 0; Idx < NumAcceptedAlignments; Idx++)
		{
		if (pObsCIGAR->ReadLen == MaxProfileSeqLen)	// some may have slipped through before it was known what read lengths were in original sequenced reads prior to any alignment trimming
			{
			NumSEs += 1;
			m_ObsCIGARBuffLen += sprintf(&m_pObsCIGARBuff[m_ObsCIGARBuffLen], "%d,%d,%c,%s,%s\n", NumSEs, pObsCIGAR->ReadLen,
				(pObsCIGAR->Flags & 0x10) ? '-' : '+', &pObsCIGAR->NameCIGARErrProfile[pObsCIGAR->NameLen], &pObsCIGAR->NameCIGARErrProfile[pObsCIGAR->NameLen + pObsCIGAR->ObsCIGARlen]);
			if((m_ObsCIGARBuffLen + 1000) >= m_AllocdObsCIGARBuffSize)
				{
				CUtility::SafeWrite(m_hObsSIGARs, m_pObsCIGARBuff, m_ObsCIGARBuffLen);
				m_ObsCIGARBuffLen = 0;
				}
			}
		else
			NumShortReads++;
		pObsCIGAR = (tsBMObsCIGAR*)((uint8_t*)pObsCIGAR + pObsCIGAR->Size);
		}
	if (m_ObsCIGARBuffLen)
		{
		CUtility::SafeWrite(m_hObsSIGARs, m_pObsCIGARBuff, m_ObsCIGARBuffLen);
		m_ObsCIGARBuffLen = 0;
		}
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Reported %d SE read error profiles to '%s' with read lengths of %dbp", NumSEs, pszObsCIGARsFile, MaxProfileSeqLen);
	}

gDiagnostics.DiagOut(eDLInfo, gszProcName, "Unable to accept reads because -");
gDiagnostics.DiagOutMsgOnly(eDLInfo, "missing features: %u, unmapped: %u, CIGAR missing: %u, CIGARs: %u, CIGAR len : % u",
	NumMissingFeatures, NumUnmapped, NumCIGARSUnknown, NumCIGARErrs, NumCIGARLenErrs);
gDiagnostics.DiagOutMsgOnly(eDLInfo, "no mate: %u, unable to pair: %d, PE strands: %u, not primary: %u",
	NumNoMateReads, NumUnableToPair, NumPEInconsistentStrands, NumSecAlignments);
gDiagnostics.DiagOutMsgOnly(eDLInfo, "short reads: %u, len: %u, frag size: %u",
	NumShortReads, NumReadLenBounds, NumFragSizeBounds);
gDiagnostics.DiagOutMsgOnly(eDLInfo, "unknown ref: %u,  ,soft clip neg loci: %u",
	NumUnknownRefs, NumNonCvtdSoftClips);

delete []pszLine;
delete pSAMalign;
Reset();
return(0);
}


int
CBenchmark::SimReads(bool bPEReads,		// true if PE pair processing otherwise SE reads
		int MaxNumReads,			// maximum number of simulated reads or read pairs
		char* pszCIGARsFile,		// read in observed CIGARs from this file
		char* pszRefGenomeFile,		// reads are against this target genome file
		char* pszSEReads,			// simulated reads are output to this file (SE or PE1 if PE)
		char* pszPE2Reads)			// simulated PE2 reads are output to this file if simulating PE reads
{
int Rslt;
uint8_t SESeq[cBMMaxReadLen];
uint8_t PE2Seq[cBMMaxReadLen];
int PE2StartLoci;
int PE2ChromID;
int SEStartLoci;
int SEChromID;
tsBMPackedCIGARs* pCurCIGAR;
int SeqIdx;
int NumSimReads;
uint32_t CurCIGAR;
int CIGARIdx;
etCIGAROpType CurOpType;
uint8_t* pPermuteBase;
int CurOpCnt;
int ReadSeqIdx;
int CurLoci;

Reset();
m_bPEReads = bPEReads;
m_MaxNumReads = MaxNumReads;
strcpy(m_szObsCIGARsFile, pszCIGARsFile);
strcpy(m_szRefGenomeFile, pszRefGenomeFile);
strcpy(m_szSEReads, pszSEReads);
if(bPEReads)
	strcpy(m_szPE2Reads, pszPE2Reads);

#ifdef _WIN32
m_hSEReads = open(pszSEReads, O_CREATETRUNC);
#else
if ((m_hSEReads = open(pszSEReads, O_RDWR | O_CREAT, S_IREAD | S_IWRITE)) != -1)
if (ftruncate(m_hSEReads, 0) != 0)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "SimReads: Unable to truncate %s - %s", pszSEReads, strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}
#endif

if (m_hSEReads < 0)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "SimReads: unable to create/truncate output file '%s'", pszSEReads);
	Reset();
	return(eBSFerrCreateFile);
	}

if((m_pszSESimReadBuff = new char [cBMCIGARAllocBuffSize]) == NULL)
	{
	Reset();
	return(eBSFerrMem);
	}
m_AllocOutBuffSize = cBMCIGARAllocBuffSize;


if (bPEReads)
{
#ifdef _WIN32
	m_hPE2Reads = open(pszPE2Reads, O_CREATETRUNC);
#else
	if ((m_hPE2Reads = open(pszPE2Reads, O_RDWR | O_CREAT, S_IREAD | S_IWRITE)) != -1)
		if (ftruncate(m_hPE2Reads, 0) != 0)
		{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "SimReads: Unable to truncate %s - %s", pszPE2Reads, strerror(errno));
			Reset();
			return(eBSFerrCreateFile);
		}
#endif

	if (m_hPE2Reads < 0)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "SimReads: unable to create/truncate output file '%s'", pszPE2Reads);
		Reset();
		return(eBSFerrCreateFile);
		}
	if ((m_pszPE2SimReadBuff = new char[cBMCIGARAllocBuffSize]) == NULL)
		{
		Reset();
		return(eBSFerrMem);
		}
	}

if ((Rslt = LoadGenome(pszRefGenomeFile)) != eBSFSuccess)
	{
	Reset();
	return(Rslt);
	}
if((Rslt = InitObsErrProfile(pszCIGARsFile)) != eBSFSuccess)	// read from this observed alignment error profiles file (as generated from actual SAM aligments)))
	{
	Reset();
	return(Rslt);
	}

gDiagnostics.DiagOut(eDLInfo, gszProcName, "SimReads: Starting to simulate %s ...", m_bPEReads ? "PE read pairs" : "SE reads");
// iterate over all the error profiles untill MaxNumReads sequences or pair of sequences have been generated
NumSimReads = 0;
pCurCIGAR = (tsBMPackedCIGARs*)m_pObsErrProfiles;
for(SeqIdx = 0; SeqIdx < MaxNumReads; SeqIdx++)
	{
	int Retries = 1000000;
	while(--Retries)
		{
		// randomly select a target chrom and starting loci
		// iterate over the observed CIGAR and build up a sequence from starting loci with induced errors as per CIGAR
		// output as the PE1 or SE read
		if(!m_bPEReads)
			pCurCIGAR->PEInsertSize = RefSeqConsumedLen(pCurCIGAR->PE1NumCIGAROps, pCurCIGAR->CIGAROpsErrProfOps);
		if((SEStartLoci = ChooseRandStartLoci(pCurCIGAR->PEInsertSize,&SEChromID)) < 0)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "SimReads: unable to locate a random chrom.loci for fragment size %u", pCurCIGAR->PEInsertSize);
			break;
			}
		PE2ChromID = SEChromID;
		
		// iterate over CIGARs
		CurLoci = SEStartLoci;
		CIGARIdx = pCurCIGAR->PE1NumCIGAROps;
		for(ReadSeqIdx = 0; ReadSeqIdx < pCurCIGAR->ReadLen; )
			{
			CurCIGAR = pCurCIGAR->CIGAROpsErrProfOps[CIGARIdx++];
			if (CIGARIdx > pCurCIGAR->PE1NumCIGAROps + pCurCIGAR->PE1NumErrProfOps)
				break;
			CurOpType = (etCIGAROpType)(CurCIGAR & 0x0000f);
			CurOpCnt = (CurCIGAR >> 4) & 0x0fffffff;
			switch(CurOpType) {
				case eCOPMatch:			//	'M' aligned but could be either matching or mismatching, here it is treated as though exactly matching
				case eCOPALignMatch:	//  '=' aligned as exactly matching, consume both query and reference sequence
					m_pGenome->GetSeq(SEChromID,CurLoci, &SESeq[ReadSeqIdx], CurOpCnt);
					CurLoci += CurOpCnt;
					ReadSeqIdx += CurOpCnt;
					continue;
		
				case eCOPALignMismatch:	// 'X' aligned but as a mismatch, consumes both query and reference sequence
					m_pGenome->GetSeq(SEChromID, CurLoci, &SESeq[ReadSeqIdx], CurOpCnt); // get sequence but then permute
					pPermuteBase = &SESeq[ReadSeqIdx];
					for(int PermuteIdx = 0; PermuteIdx < CurOpCnt; PermuteIdx++, pPermuteBase++)
						*pPermuteBase = (*pPermuteBase + 2) & 0x03;
					CurLoci += CurOpCnt;
					ReadSeqIdx += CurOpCnt;
					continue;

				case eCOPSoftClip:				//	'S'  soft clipping - consumes query sequence only
				case eCOPInsert:				//	'I'  insertion relative to target - consumes query sequence only
					m_pGenome->GetSeq(SEChromID, CurLoci, &SESeq[ReadSeqIdx], CurOpCnt); // get sequence but then permute so will not match when aligning
					pPermuteBase = &SESeq[ReadSeqIdx];
					for (int PermuteIdx = 0; PermuteIdx < CurOpCnt; PermuteIdx++, pPermuteBase++)
						*pPermuteBase = (*pPermuteBase + 2) & 0x03;
					ReadSeqIdx += CurOpCnt;		// no changes to target loci
					continue;

				case eCOPSkipRegion:    	    //	'N'  skipped region relative to target - intron?  - consumes reference sequence only
				case eCOPDelete:				// 'D'  deletion relative to target - consumes reference sequence only
					CurLoci += CurOpCnt;		// no change to ReadSeqIdx
					continue;

				case eCOPHardClip:				// 'H'  hard clipping
					m_pGenome->GetSeq(SEChromID, CurLoci, &SESeq[ReadSeqIdx], CurOpCnt); // get sequence but then permute so shouldn't match when aligning
					pPermuteBase = &SESeq[ReadSeqIdx];
					for (int PermuteIdx = 0; PermuteIdx < CurOpCnt; PermuteIdx++, pPermuteBase++)
						*pPermuteBase = (*pPermuteBase + 2) & 0x03;
					ReadSeqIdx += CurOpCnt;		// no changes to target loci
					continue;
				}
			}
		// need to ensure that reads have at most only a couple of indeterminates otherwise alignments likely to be very difficult
		int NumNs = 0;
		for (int BaseIdx = 0; BaseIdx < pCurCIGAR->ReadLen; BaseIdx++)
			if (SESeq[BaseIdx] >= eBaseN)
				NumNs++;
		if (NumNs > pCurCIGAR->ReadLen / 5)
			continue;				// retry for a read without excessive number of indeterminates

		if (m_bPEReads)
			{
			CIGARIdx = pCurCIGAR->PE1NumCIGAROps + pCurCIGAR->PE1NumErrProfOps + pCurCIGAR->PE2NumCIGAROps;
			PE2StartLoci = pCurCIGAR->PEInsertSize + SEStartLoci - RefSeqConsumedLen(pCurCIGAR->PE2NumErrProfOps, &pCurCIGAR->CIGAROpsErrProfOps[CIGARIdx]);
			CurLoci = PE2StartLoci;
			for (ReadSeqIdx = 0; ReadSeqIdx < pCurCIGAR->ReadLen; )
				{
				CurCIGAR = pCurCIGAR->CIGAROpsErrProfOps[CIGARIdx++];
				if (CIGARIdx > pCurCIGAR->PE1NumCIGAROps + pCurCIGAR->PE1NumErrProfOps + pCurCIGAR->PE2NumCIGAROps + pCurCIGAR->PE2NumErrProfOps)
					break;
				CurOpType = (etCIGAROpType)(CurCIGAR & 0x0000f);
				CurOpCnt = (CurCIGAR >> 4) & 0x0fffffff;
				switch (CurOpType) {
					case eCOPMatch:			//	'M' aligned but could be either matching or mismatching, here it is treated as though exactly matching
					case eCOPALignMatch:	//  '=' aligned as exactly matching, consume both query and reference sequence
						m_pGenome->GetSeq(SEChromID, CurLoci, &PE2Seq[ReadSeqIdx], CurOpCnt);
						CurLoci += CurOpCnt;
						ReadSeqIdx += CurOpCnt;
						continue;

					case eCOPALignMismatch:	// 'X' aligned but as a mismatch, consumes both query and reference sequence
						m_pGenome->GetSeq(SEChromID, CurLoci, &PE2Seq[ReadSeqIdx], CurOpCnt); // get sequence but then permute
						pPermuteBase = &PE2Seq[ReadSeqIdx];
						for (int PermuteIdx = 0; PermuteIdx < CurOpCnt; PermuteIdx++, pPermuteBase++)
							*pPermuteBase = (*pPermuteBase + 2) & 0x03;
						CurLoci += CurOpCnt;
						ReadSeqIdx += CurOpCnt;
						continue;

					case eCOPSoftClip:				//	'S'  soft clipping - consumes query sequence only
					case eCOPInsert:				//	'I'  insertion relative to target - consumes query sequence only
						m_pGenome->GetSeq(SEChromID, CurLoci, &PE2Seq[ReadSeqIdx], CurOpCnt); // get sequence but then permute so will not match when aligning
						pPermuteBase = &PE2Seq[ReadSeqIdx];
						for (int PermuteIdx = 0; PermuteIdx < CurOpCnt; PermuteIdx++, pPermuteBase++)
							*pPermuteBase = (*pPermuteBase + 2) & 0x03;
						ReadSeqIdx += CurOpCnt;		// no changes to target loci
						continue;

					case eCOPSkipRegion:    	    //	'N'  skipped region relative to target - intron?  - consumes reference sequence only
					case eCOPDelete:				// 'D'  deletion relative to target - consumes reference sequence only
						CurLoci += CurOpCnt;		// no change to ReadSeqIdx
						continue;

					case eCOPHardClip:				// 'H'  hard clipping
						m_pGenome->GetSeq(SEChromID, CurLoci, &PE2Seq[ReadSeqIdx], CurOpCnt); // get sequence but then permute so shouldn't match when aligning
						pPermuteBase = &PE2Seq[ReadSeqIdx];
						for (int PermuteIdx = 0; PermuteIdx < CurOpCnt; PermuteIdx++, pPermuteBase++)
							*pPermuteBase = (*pPermuteBase + 2) & 0x03;
						ReadSeqIdx += CurOpCnt;		// no changes to target loci
						continue;
					}
				}
			}
		// need to ensure that reads have at most only a couple of indeterminates otherwise alignments likely to be very problematic
		NumNs = 0;
		for(int BaseIdx = 0; BaseIdx < pCurCIGAR->ReadLen; BaseIdx++)
			if(PE2Seq[BaseIdx] >= eBaseN)
				NumNs++;
		if (NumNs > pCurCIGAR->ReadLen / 5)
			continue;				// retry for a read without excessive number of indeterminates
		break;
		}

	if(!Retries)		// if exhusted retries looking for reads without too many indeterminates
		{
		gDiagnostics.DiagOut(eDLWarn, gszProcName, "SimReads: Unable to locate ground truth read with less than 5% of indeterminate bases");
		Reset();
		return(eBSFerrRead);
		}

	if(m_bPEReads)
		{
		// just a check; end loci - start loci should match the fragment length, if not then there is a major problem!
		int InsertSize = CurLoci - SEStartLoci;
		if(InsertSize != pCurCIGAR->PEInsertSize)
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "SimReads: Expected insert size of %d, sim was %d", pCurCIGAR->PEInsertSize, InsertSize);
		}

	NumSimReads += 1;
	// if read strand was '-' then reverse complement before writing out
	if(pCurCIGAR->FlgPE1Strand)
		CSeqTrans::ReverseComplement(pCurCIGAR->ReadLen,SESeq);
	WriteFastaFile(m_bPEReads,false, NumSimReads, pCurCIGAR, SEChromID,SEStartLoci, SESeq);

	if (m_bPEReads)
		{
		// double check on strands - if PE then can't have both reads on same strand!
		if(pCurCIGAR->FlgPE2Strand == pCurCIGAR->FlgPE1Strand)
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "SimReads: both reads are on same strand!");
		if(pCurCIGAR->FlgPE2Strand)
			CSeqTrans::ReverseComplement(pCurCIGAR->ReadLen, PE2Seq);
		WriteFastaFile(m_bPEReads,true, NumSimReads, pCurCIGAR, PE2ChromID, PE2StartLoci, PE2Seq);
		}
	
	pCurCIGAR = (tsBMPackedCIGARs*)((uint8_t *)pCurCIGAR + pCurCIGAR->Size);
	if((uint8_t*)pCurCIGAR >= m_pObsErrProfiles + m_UsedObsErrProfMem)	// recycle error profiles by starting back at first
		pCurCIGAR = (tsBMPackedCIGARs*)m_pObsErrProfiles;
	}

	if (m_OutSEBuffIdx)
		{
		if (!CUtility::SafeWrite(m_hSEReads, m_pszSESimReadBuff, m_OutSEBuffIdx))
			return(eBSFerrFileAccess);
		m_OutSEBuffIdx = 0;
#ifdef _WIN32
		_commit(m_hSEReads);
#else
		fsync(m_hSEReads);
#endif
		close(m_hSEReads);
		m_hSEReads = -1;
		}

	if (m_OutPE2BuffIdx)
		{
		if (!CUtility::SafeWrite(m_hPE2Reads, m_pszPE2SimReadBuff, m_OutPE2BuffIdx))
			return(eBSFerrFileAccess);
		m_OutPE2BuffIdx = 0;
#ifdef _WIN32
		_commit(m_hPE2Reads);
#else
		fsync(m_hPE2Reads);
#endif
		close(m_hPE2Reads);
		m_hPE2Reads = -1;
		}
gDiagnostics.DiagOut(eDLInfo, gszProcName, "SimReads: Generated %d %s", NumSimReads, m_bPEReads ? "PE read pairs":"SE reads");
Reset();
return(Rslt);
}

int					// returns num of chars in decoded packed CIGAR string excluding terminating '\0'
CBenchmark::DecodeCIGAR(int NumPackedOps,			// number of packed CIGAR ops in
	uint32_t* pPackedOPs,		// ptr to packed CIGAR ops
	int MaxCIGARChrs,			// pszCIGAR allocated to hold at most this many chars including terminating '\0'
	char* pszCIGAR)				// decode into this buffer
{
int Idx;
int strCIGAROfs;
etCIGAROpType OpType;
char ChrOp;
int NumOps;
int NumChrs;
int Cnt;

if(NumPackedOps < 1 || pPackedOPs == NULL || MaxCIGARChrs <= 3 || pszCIGAR == NULL)
	return(0);
strCIGAROfs = 0;
for(Idx = 0; Idx < NumPackedOps; Idx++, pPackedOPs++)
	{
	OpType = (etCIGAROpType)(*pPackedOPs & 0x0f);
	NumOps = (*pPackedOPs >> 4) & 0x0fffffff;
	Cnt = NumOps;
	NumChrs = 2;
	while((Cnt /= 10) > 0)
		NumChrs+=1;
	if(NumChrs >= MaxCIGARChrs)
		return(0);

	switch (OpType) {
		case eCOPMatch:			// 'M' aligned but could be either matching or mismatching, consumes query and reference sequence
			ChrOp = 'M';
			break;

		case eCOPALignMatch:	// '=' aligned as exactly matching, consumes both query and reference sequence
			ChrOp = '=';
			break;

		case eCOPALignMismatch:	// 'X' aligned but as a mismatch, consumes both query and reference sequence
			ChrOp = 'X';
			break;
		
		case eCOPInsert:			//'I'  insertion relative to target - consumes query sequence only
			ChrOp = 'I';
			break;

		case eCOPSkipRegion:    	    // 'N'  skipped region relative to target - intron?  - consumes reference sequence only
			ChrOp = 'N';
			break;

		case eCOPDelete:				// 'D'  deletion relative to target - consumes reference sequence only
			ChrOp = 'D';
			break;

		case eCOPSoftClip:				//	'S'  soft clipping - consumes query sequence only
			ChrOp = 'S';
			break;

		case eCOPHardClip:				//	'H'  hard clipping - consumes neither query or reference sequence
			ChrOp = 'H';
			break;

		case eCOPPadding:				//	'P' padding (silent deletion from padded reference) consumes neither query or reference sequence
			ChrOp = 'P';
			break;

		default:		// unrecognised CIGAR operation
			return(-1);
		}
	strCIGAROfs += sprintf(&pszCIGAR[strCIGAROfs],"%d%c", NumOps, ChrOp);
	MaxCIGARChrs -= NumChrs;
	}
return((int)strlen(pszCIGAR));
}

int										// eBSFSuccess if sequence written to file otherwise the error code
CBenchmark::WriteFastaFile(bool bPEReads,	// true if simulations are for PE reads
	bool bPE2,					// true if writing out a PE2 simulated read, false if SE or PE1
	int SeqID,					// primary sequence identifier
	tsBMPackedCIGARs* pCurCIGAR,	// CIGAR error profile 
	int ChromID,				// simulated read is on this chromosome 
	int StartLoci,				// starting at this loci (0 based)
	etSeqBase* pSeq)			// sequence
{
int SeqIdx;
int LineLen;
char Base;
char SEPE;
int *pOutBuffIdx;
int hSimReads;
char *pOutBuff;
char* pChr;
char Strand;
etSeqBase* pBase;
char szReadNameSfx[5];
char szChromName[cMaxDatasetSpeciesChrom+1];
char szCIGAR[100];

if (pCurCIGAR == NULL || (!bPEReads && bPE2) ||
	(!bPE2 && (m_hSEReads == -1 || m_pszSESimReadBuff == NULL)) || (bPE2 && (m_hPE2Reads == -1 || m_pszSESimReadBuff == NULL)))
	return(eBSFerrParams);

if (pCurCIGAR->ReadLen < cBMMinReadLen)
	return(eBSFerrParams);

if (bPE2)
	{
	pOutBuffIdx = &m_OutPE2BuffIdx;
	pOutBuff = m_pszPE2SimReadBuff;
	hSimReads = m_hPE2Reads;
	}
else
	{
	pOutBuffIdx = &m_OutSEBuffIdx;
	pOutBuff = m_pszSESimReadBuff;
	hSimReads = m_hSEReads;
	}
pChr = &pOutBuff[*pOutBuffIdx];

if (*pOutBuffIdx > (m_AllocOutBuffSize - (pCurCIGAR->ReadLen * 2)))
	{
	if (!CUtility::SafeWrite(hSimReads, pOutBuff, *pOutBuffIdx))
		return(eBSFerrFileAccess);
	*pOutBuffIdx = 0;
	pChr = pOutBuff;
	}

	// build descriptor
m_pGenome->GetIdentName(ChromID, cMaxDatasetSpeciesChrom,szChromName);
szChromName[cMaxDatasetSpeciesChrom] = '\0';
if(bPE2)
	{
	Strand = pCurCIGAR->FlgPE2Strand ? '-' : '+';
	strcpy(szReadNameSfx,"PE.2");
	SEPE = '2';
	DecodeCIGAR(pCurCIGAR->PE2NumCIGAROps,&pCurCIGAR->CIGAROpsErrProfOps[pCurCIGAR->PE1NumCIGAROps + pCurCIGAR->PE1NumErrProfOps], sizeof(szCIGAR),szCIGAR);
	}
else
	{
	Strand = pCurCIGAR->FlgPE1Strand ? '-' : '+';
	if(bPEReads)
		{
 		strcpy(szReadNameSfx, "PE.1");
		SEPE = '1';
		}
	else
		{
		strcpy(szReadNameSfx, "SE");
		SEPE = '0';
		}
	DecodeCIGAR(pCurCIGAR->PE1NumCIGAROps, pCurCIGAR->CIGAROpsErrProfOps, sizeof(szCIGAR), szCIGAR);
	}
int Len = sprintf(pChr, ">SR%d%s %c %d %s %d %c %s %d\n", SeqID, szReadNameSfx, SEPE, pCurCIGAR->ReadLen, szChromName, StartLoci+1, Strand,szCIGAR, pCurCIGAR->ID);
*pOutBuffIdx += Len;
pChr += Len;
LineLen = 0;

pBase = pSeq;
for (SeqIdx = 0; SeqIdx < pCurCIGAR->ReadLen; SeqIdx++, pBase += 1)
	{
	switch (*pBase) {
		case eBaseA:
			Base = 'A';
			break;
		case eBaseC:
			Base = 'C';
			break;
		case eBaseG:
			Base = 'G';
			break;
		case eBaseT:
			Base = 'T';
			break;
		default:
			Base = 'N';
			break;
		}
	*pChr++ = Base;
	*pOutBuffIdx += 1;
	LineLen += 1;
	if (LineLen > 70)
		{
		*pChr++ = '\n';
		*pOutBuffIdx += 1;
		LineLen = 0;
		}
	}
if(LineLen)
	{
	*pChr++ = '\n';
	*pOutBuffIdx += 1;
	}
return(eBSFSuccess);
}


int
CBenchmark::LoadGroundTruths(char* pszSESimReads,	// load ground truths for SE or PE1 if PE simulated reads
				char* pszPE2SimReads)	// load ground truths for PE2 if PE simulated reads
{
int Rslt;
int NumDescrs;
int NumEls;
char szDescriptor[301];
int DescrLen;
char szReadName[80];
int SubID;
int ReadLen;
char szRefChrom[80];
uint32_t StartLoci;
char Strand;
char szCIGAR[251];
uint32_t ReadID;
tsBMGroundTruth *pGroundTruth;

if(pszSESimReads == NULL || pszSESimReads[0] == '\0' || (m_bPEReads && (pszPE2SimReads == NULL || pszPE2SimReads[0] == '\0')))
	return(eBSFerrParams);

m_NumGroundTruths = 0;
m_UsedGroundTruthsMem = 0;

m_AllocdGroundTruthsMem = (size_t)cSRAllocdObsErrProfMemChunk;
#ifdef _WIN32
m_pGroundTruths = (UINT8*)malloc(m_AllocdGroundTruthsMem);
if (m_pGroundTruths == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadGroundTruths: Memory allocation of %lld bytes failed", (INT64)m_AllocdGroundTruthsMem);
	m_AllocdGroundTruthsMem = 0;
	Reset();
	return(eBSFerrMem);
	}
#else
// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
m_pGroundTruths = (UINT8*)mmap(NULL, (size_t)m_AllocdGroundTruthsMem, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
if (m_pGroundTruths == MAP_FAILED)
{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadGroundTruths: Memory allocation of %lld bytes through mmap()  failed", (INT64)m_AllocdGroundTruthsMem, strerror(errno));
	m_pGroundTruths = NULL;
	m_AllocdGroundTruthsMem = 0;
	Reset();
	return(eBSFerrMem);
}
#endif

	// parse fasta file containing ground truths
if ((m_pSESimReads = new CFasta()) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to instantiate CFasta");
	Reset();
	return(eBSFerrObj);
	}
if ((Rslt = m_pSESimReads->Open(pszSESimReads)) != eBSFSuccess)
	{
	while (m_pSESimReads->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal, gszProcName, m_pSESimReads->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to open input fasta file containing simulated reads '%s'", pszSESimReads);
	Reset();
	return(Rslt);
	}

// iterate over all reads and parse their descriptor lines as these contain the ground truth for each individual read
// at same time accumulate best case base match score
m_MaxBaseAlignmentScore = 0;
int64_t Score;
NumDescrs = 0;
while ((Rslt = m_pSESimReads->ReadSequence()) > eBSFSuccess)
	{
	if (Rslt == eBSFFastaDescr)		// just read a descriptor line, parse out the feature identifier
		{
		NumDescrs += 1;
		DescrLen = m_pSESimReads->ReadDescriptor(szDescriptor, sizeof(szDescriptor));
		if ((NumEls = sscanf(szDescriptor, "%80s %d %d %80s %u %c %200s %u", szReadName, &SubID, &ReadLen, szRefChrom, &StartLoci, &Strand, szCIGAR, &ReadID)) != 8)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to parse fasta descriptor %d for ground truth of originating loci from '%s'", NumDescrs, pszSESimReads);
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Descriptor: %s", szDescriptor);
			Reset();
			return(eBSFerrFastqDescr);
			}
		m_TotGroundTruthBases += ReadMatchScore(ReadLen, szCIGAR, m_AlignedBasePts, &Score);
		m_MaxBaseAlignmentScore += Score;
		
		if ((cSRAllocdObsErrProfMemChunk / 10) >= (m_AllocdGroundTruthsMem - m_UsedGroundTruthsMem))
		{
			size_t memreq = m_AllocdGroundTruthsMem + (size_t)cSRAllocdObsErrProfMemChunk;
			UINT8* pTmp;
#ifdef _WIN32
			pTmp = (UINT8*)realloc(m_pGroundTruths, memreq);
#else
			pTmp = (UINT8*)mremap(m_pGroundTruths, m_AllocdGroundTruthsMem, memreq, MREMAP_MAYMOVE);
			if (pTmp == MAP_FAILED)
				pTmp = NULL;
#endif
			if (pTmp == NULL)
			{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadGroundTruths: Memory re-allocation to %lld bytes - %s", (INT64)(memreq), strerror(errno));
				Reset();
				return(eBSFerrMem);
			}
			m_pGroundTruths = pTmp;
			m_AllocdGroundTruthsMem = memreq;
		}
		pGroundTruth = (tsBMGroundTruth *)&m_pGroundTruths[m_UsedGroundTruthsMem];
		memset(pGroundTruth,0,sizeof(tsBMGroundTruth));
		pGroundTruth->ID = ++m_NumGroundTruths;
		pGroundTruth->ReadID = ReadID;
		pGroundTruth->NameLen = (int8_t)strlen(szReadName);
		pGroundTruth->ChromNameLen = (int8_t)strlen(szRefChrom);
		pGroundTruth->CIGARLen = (int8_t)strlen(szCIGAR);
		pGroundTruth->FlgPE2 = (SubID == 2) ? 1 : 0;
		pGroundTruth->FlgStrand = Strand == '-' ? 1 : 0;
		pGroundTruth->ReadLen = ReadLen;
		pGroundTruth->StartLoci = StartLoci - 1;		// required as currently 1 based to suit SAM requirements but later will be comparing with a 0 based loci
		strcpy((char *)pGroundTruth->NameChromCIGAR, szReadName);
		strcpy((char*)&pGroundTruth->NameChromCIGAR[pGroundTruth->NameLen + 1], szRefChrom);
		strcpy((char*)&pGroundTruth->NameChromCIGAR[pGroundTruth->NameLen + pGroundTruth->ChromNameLen + 2], szCIGAR);
		pGroundTruth->Size = sizeof(tsBMGroundTruth) + pGroundTruth->NameLen + pGroundTruth->ChromNameLen + pGroundTruth->CIGARLen + 3;
		m_UsedGroundTruthsMem += pGroundTruth->Size;
		}
	}
m_pSESimReads->Close();
delete m_pSESimReads;
m_pSESimReads = NULL;

if (m_bPEReads && pszPE2SimReads != NULL)
	{
	if ((m_pPE2SimReads = new CFasta()) == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to instantiate CFasta");
		Reset();
		return(eBSFerrObj);
		}
	if ((Rslt = m_pPE2SimReads->Open(pszPE2SimReads)) != eBSFSuccess)
		{
		while (m_pPE2SimReads->NumErrMsgs())
			gDiagnostics.DiagOut(eDLFatal, gszProcName, m_pPE2SimReads->GetErrMsg());
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to open input fasta file containing PE2 simulated reads '%s'", pszPE2SimReads);
		Reset();
		return(Rslt);
		}
	NumDescrs = 0;
	while ((Rslt = m_pPE2SimReads->ReadSequence()) > eBSFSuccess)
		{
		if (Rslt == eBSFFastaDescr)		// just read a descriptor line, parse out the feature identifier
			{
			NumDescrs += 1;
			DescrLen = m_pPE2SimReads->ReadDescriptor(szDescriptor, sizeof(szDescriptor));
			if ((NumEls = sscanf(szDescriptor, "%80s %d %d %80s %u %c %80s %u", szReadName, &SubID, &ReadLen, szRefChrom, &StartLoci, &Strand, szCIGAR, &ReadID)) != 8) 
			    {
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to parse fasta descriptor %d for ground truth of originating loci from '%s'", NumDescrs, pszPE2SimReads);
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Descriptor: %s", szDescriptor);
				Reset();
				return(eBSFerrFastqDescr);
				}
			m_TotGroundTruthBases += ReadMatchScore(ReadLen, szCIGAR, m_AlignedBasePts, &Score);
			m_MaxBaseAlignmentScore += Score;

			if ((cSRAllocdObsErrProfMemChunk / 10) >= (m_AllocdGroundTruthsMem - m_UsedGroundTruthsMem))
			{
				size_t memreq = m_AllocdGroundTruthsMem + (size_t)cSRAllocdObsErrProfMemChunk;
				UINT8* pTmp;
#ifdef _WIN32
				pTmp = (UINT8*)realloc(m_pGroundTruths, memreq);
#else
				pTmp = (UINT8*)mremap(m_pGroundTruths, m_AllocdGroundTruthsMem, memreq, MREMAP_MAYMOVE);
				if (pTmp == MAP_FAILED)
					pTmp = NULL;
#endif
				if (pTmp == NULL)
				{
					gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadGroundTruths: Memory re-allocation to %lld bytes - %s", (INT64)(memreq), strerror(errno));
					Reset();
					return(eBSFerrMem);
				}
				m_pGroundTruths = pTmp;
				m_AllocdGroundTruthsMem = memreq;
			}
			pGroundTruth = (tsBMGroundTruth*)&m_pGroundTruths[m_UsedGroundTruthsMem];
			memset(pGroundTruth, 0, sizeof(tsBMGroundTruth));
			pGroundTruth->ID = ++m_NumGroundTruths;
			pGroundTruth->ReadID = ReadID;
			pGroundTruth->NameLen = (int8_t)strlen(szReadName);
			pGroundTruth->ChromNameLen = (int8_t)strlen(szRefChrom);
			pGroundTruth->CIGARLen = (int8_t)strlen(szCIGAR);
			pGroundTruth->FlgPE2 = (SubID == 2) ? 1 : 0;
			pGroundTruth->FlgStrand = Strand == '-' ? 1 : 0;
			pGroundTruth->ReadLen = ReadLen;
			pGroundTruth->StartLoci = StartLoci - 1;
			strcpy((char*)pGroundTruth->NameChromCIGAR, szReadName);
			strcpy((char*)&pGroundTruth->NameChromCIGAR[pGroundTruth->NameLen + 1], szRefChrom);
			strcpy((char*)&pGroundTruth->NameChromCIGAR[pGroundTruth->NameLen + pGroundTruth->ChromNameLen + 2], szCIGAR);
			pGroundTruth->Size = sizeof(tsBMGroundTruth) + pGroundTruth->NameLen + pGroundTruth->ChromNameLen + pGroundTruth->CIGARLen + 3;
			m_UsedGroundTruthsMem += pGroundTruth->Size;
			}
		}
	m_pPE2SimReads->Close();
	delete m_pPE2SimReads;
	m_pPE2SimReads = NULL;
	}

	// have loaded the ground truths, create an index so can do quick binary searches against read names
m_ppGroundTruthIdx = new tsBMGroundTruth * [m_NumGroundTruths];
pGroundTruth = (tsBMGroundTruth*)m_pGroundTruths;
for (uint32_t Idx = 0; Idx < m_NumGroundTruths; Idx++)
	{
	m_ppGroundTruthIdx[Idx] = pGroundTruth;
	pGroundTruth = (tsBMGroundTruth*)((uint8_t*)pGroundTruth + pGroundTruth->Size);
	}
gDiagnostics.DiagOut(eDLInfo, gszProcName, "Sorting ground truths by read name");
qsort(m_ppGroundTruthIdx, m_NumGroundTruths, sizeof(tsBMGroundTruth *), SortGroundTruths);

return(1);
}

tsBMGroundTruth*		// returned ground truth which matches
CBenchmark::LocateGroundTruth(tsBAMalign* pSAMalign)	// this alignment by name
{
int TargPsn;
int IdxLo;
int IdxHi;
int Rslt;
tsBMGroundTruth *pGroundTruth;
IdxLo = 0;
IdxHi = m_NumGroundTruths - 1;
do {
	TargPsn = (IdxLo + IdxHi) / 2;
	pGroundTruth = m_ppGroundTruthIdx[TargPsn];
	if((Rslt = stricmp((char *)pGroundTruth->NameChromCIGAR, pSAMalign->read_name))==0)
		return(pGroundTruth);
	
	if (Rslt > 0)
		IdxHi = TargPsn - 1;
	else
		IdxLo = TargPsn + 1;
	} while (IdxHi >= IdxLo);
return(NULL);
}


//Op BAM Description                                          Consumes query    Consumes reference 
//M  0   alignment match(can be a sequence match or mismatch) yes               yes 
//I  1   insertion to the reference                           yes               no 
//D  2   deletion from the reference                          no                yes 
//N  3   skipped region from the reference                    no                yes 
//S  4   soft clipping(clipped sequences present in SEQ)      yes               no 
//H  5   hard clipping(clipped sequences NOT present in SEQ)  no                no 
//P  6   padding(silent deletion from padded reference)       no                no 
//=  7   sequence match                                       yes               yes
//X  8   sequence mismatch                                    yes               yes
int
CBenchmark::Score(bool bPrimaryOnly,// score only primary alignments otherwise score including secondary
		bool bPEReads,				// if true only PE pair processing otherwise score as if SE reads
		int UnalignedBasePts,		// points to apply for each unaligned base
		int AlignedBasePts,			// points to apply for each base aligned to it's ground truth loci
		int SilentTrimAlignBasePts, // points to apply for each base aligned in a silently trimmed read within ground truth loci range
		int MisalignedBasePts,		// points to apply for each base aligned but not to ground truth loci
		char* pszSEReads,			// input simulated reads which contain ground truths from this file for SE or PE1 if PE
		char* pszPE2Reads,			// input simulated reads which contain ground truths from this file for PE2 if PE
		char* pszAlignmentsFile)	// input file containing alignments of simulated reads (SAM or BAM)
{
int Rslt;
Reset();
m_bPrimaryOnly = bPrimaryOnly;
m_bPEReads = bPEReads;
m_UnalignedBasePts = UnalignedBasePts;
m_AlignedBasePts = AlignedBasePts;
m_SilentTrimAlignBasePts = SilentTrimAlignBasePts;
m_MisalignedBasePts = MisalignedBasePts;
m_MaxNumReads = 0;
strcpy(m_szAlignmentsFile, pszAlignmentsFile);
strcpy(m_szSEReads, pszSEReads);
if(bPEReads)
	strcpy(m_szPE2Reads, pszPE2Reads);
else
	m_szPE2Reads[0] = '\0';
m_MaxBaseAlignmentScore = 0;
m_NumBasesLociCorrect = 0;
m_NumBasesLociIncorrect = 0;
m_NumBasesLociUnclaimed = 0;
memset(m_ReadOverlapHistogram,0,sizeof(m_ReadOverlapHistogram));

if(bPEReads)
	{
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Loading PE ground truth reads from PE1: '%s'", pszSEReads);
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Loading PE ground truth reads from PE2: '%s'", pszPE2Reads);
	}
else
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Loading SE ground truth reads from '%s'", pszSEReads);
if((Rslt = LoadGroundTruths(pszSEReads, pszPE2Reads)) < 1)
	return(Rslt);
gDiagnostics.DiagOut(eDLInfo, gszProcName, "Total of %u ground truth reads loaded having %lld ground truth bases", m_NumGroundTruths, m_TotGroundTruthBases);

gDiagnostics.DiagOut(eDLInfo, gszProcName, "Opening alignments in: '%s'", pszAlignmentsFile);
if ((Rslt = OpenAlignments(pszAlignmentsFile)) != eBSFSuccess)
	return(Rslt);

	// iterate over alignments, looking up alignment read names in the ground truth so each alignment can be scored
char* pszLine;
char* pTxt;
int LineLen;
uint32_t NumMissingFeatures;
uint32_t NumUnmapped;
uint32_t NumNotPrimary;
uint32_t NumCIGARSUnknown;
uint32_t NumNoGroundTruths;
uint32_t NumErrChroms;
uint32_t NumErrStrands;
uint32_t NumLociErrs;
uint32_t NumCorrectChroms;
uint32_t NumErrPE2;
tsBAMalign* pSAMalign;
tsBMGroundTruth *pGroundTruth;
uint32_t NumScoredAlignments;
uint32_t NumPutativeScoredAlignments;
int64_t AlignedPotentialScore;

if ((pszLine = new char[cMaxReadLen]) == NULL)				// buffer input lines
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to allocate line buffering");
	Reset();
	return(-1);
	}

if ((pSAMalign = new tsBAMalign) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to allocate memory for alignment tsBAMalign structure");
	delete pszLine;
	Reset();
	return(-1);
	}

NumMissingFeatures = 0;
NumUnmapped = 0;
NumNotPrimary = 0;
NumCIGARSUnknown = 0;
NumNoGroundTruths = 0;
NumScoredAlignments = 0;
NumPutativeScoredAlignments = 0;
AlignedPotentialScore = 0;
NumErrChroms = 0;
NumErrStrands = 0;
NumLociErrs = 0;
NumCorrectChroms = 0;
NumErrPE2 = 0;

m_TotBaseAlignmentScore = 0;
m_ScoredReads = 0;
m_UnscoredReads = 0;

gDiagnostics.DiagOut(eDLInfo, gszProcName, "Scoring alignments in: '%s'", pszAlignmentsFile);
time_t Then = time(NULL);
time_t Now;
while (Rslt >= eBSFSuccess && (LineLen = m_pAlignments->GetNxtSAMline(pszLine)) > 0)
{
	pszLine[cMaxReadLen - 1] = '\0';
	pTxt = CUtility::TrimWhitespc(pszLine);
	if (*pTxt == '\0')				// simply slough lines which are just whitespace
		continue;
	if (*pTxt == '@')				// only interested in lines with putative alignments
		continue;

	NumPutativeScoredAlignments += 1;
	if (!(NumPutativeScoredAlignments % 100000) || NumPutativeScoredAlignments == 1)
		{
		Now = time(NULL);
		if ((Now - Then) >= 60)
			{
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "Processing %d SAM/BAM CIGAR alignments, scored %d", NumPutativeScoredAlignments, NumScoredAlignments);
			Then += 60;
			}
		}

	// primary interest is in the read name, reference chromname, startloci, length
	if ((Rslt = (teBSFrsltCodes)m_pAlignments->ParseSAM2BAMalign(pTxt, pSAMalign, NULL, true)) < eBSFSuccess)
		{
		if (Rslt == eBSFerrFeature)	// not too worried if the aligned to feature is missing as some SAMs are missing header features
			{						// our interest is as to if the ref chrom is as expected, this is checked later
			NumMissingFeatures++;
			Rslt = eBSFSuccess;
			}
		else
			{
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "Processed %d SAM/BAM CIGAR alignments, missing aligned to reference", NumPutativeScoredAlignments);
			break;
			}
		}

	// check if read has been mapped, if not then slough ...
	if (pSAMalign->refID == -1 || (pSAMalign->flag_nc >> 16) & 0x04)
		{
		m_UnscoredReads++;
		NumUnmapped++;
		continue;
		}
	if (pSAMalign->cigar[0] == '*')	// is Cigar is unknown
		{
		m_UnscoredReads++;
		NumCIGARSUnknown++;
		continue;
		}

	if (m_bPrimaryOnly && (pSAMalign->flag_nc >> 16) & 0x0100)	// only interested in primary reads
		{
		m_UnscoredReads++;
		NumNotPrimary++;
		continue;
		}

	if ((pGroundTruth = LocateGroundTruth(pSAMalign)) == NULL) //fails if read name not a simulated read name
		{
		if(++NumNoGroundTruths == 1)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Alignment read name ('%s') does not match any ground truth read names - alignments nust be from ground truth simulated reads", pSAMalign->read_name);
			delete pszLine;
			delete pSAMalign;
			Reset();
			return(eBSFerrFastaDescr);
			}
		continue;
		}

	// have a ground truth
	m_ScoredReads++;
	pGroundTruth->FlgAligned = true;

	// find out how many bases in alignment could be potentially scored
	int64_t PotentialAlgnBases;
	int64_t PotentialBaseScore;
	
	PotentialAlgnBases = ReadMatchScore(pGroundTruth->ReadLen, (char*)&pGroundTruth->NameChromCIGAR[pGroundTruth->NameLen + pGroundTruth->ChromNameLen + 2], m_AlignedBasePts, &PotentialBaseScore);
	if (PotentialAlgnBases <= 0)
		{
		NumScoredAlignments++;
		m_ReadOverlapHistogram[0] += 1;
		continue;
		}
	AlignedPotentialScore += PotentialBaseScore;
	m_TotNumPotentialAlignBases += PotentialAlgnBases;
		// aligned to same chrom as ground truth?
	if(stricmp((char *)&pGroundTruth->NameChromCIGAR[pGroundTruth->NameLen+1], pSAMalign->szRefSeqName))
		{
		int64_t BaseScore;
		ReadMatchScore(pGroundTruth->ReadLen, (char *)&pGroundTruth->NameChromCIGAR[pGroundTruth->NameLen+ pGroundTruth->ChromNameLen+2], m_MisalignedBasePts, &BaseScore);
		m_TotBaseAlignmentScore += BaseScore;
		m_NumBasesLociIncorrect += PotentialAlgnBases;
		NumScoredAlignments++;
		m_ReadOverlapHistogram[0] += 1;
		pGroundTruth->FlgRefChromErr = true;
		NumErrChroms++;
		continue;
		}
	NumCorrectChroms++;
	if(pGroundTruth->FlgStrand != (((pSAMalign->flag_nc >> 16) & 0x010)== 0x010))
		{
		int64_t BaseScore;
		ReadMatchScore(pGroundTruth->ReadLen, (char*)&pGroundTruth->NameChromCIGAR[pGroundTruth->NameLen + pGroundTruth->ChromNameLen + 2], m_MisalignedBasePts, &BaseScore);
		m_TotBaseAlignmentScore += BaseScore;
		m_NumBasesLociIncorrect += PotentialAlgnBases;
		m_ReadOverlapHistogram[0] += 1;
		NumScoredAlignments++;
		pGroundTruth->FlgStrandErr = true;
		NumErrStrands++;
		continue;
		}
	if(bPEReads && (pGroundTruth->FlgPE2 != (((pSAMalign->flag_nc >> 16) & 0x080) == 0x080)))
		{
		int64_t BaseScore;
		ReadMatchScore(pGroundTruth->ReadLen, (char*)&pGroundTruth->NameChromCIGAR[pGroundTruth->NameLen + pGroundTruth->ChromNameLen + 2], m_MisalignedBasePts, &BaseScore);
		m_TotBaseAlignmentScore += BaseScore;
		m_NumBasesLociIncorrect += PotentialAlgnBases;
		m_ReadOverlapHistogram[0] += 1;
		NumScoredAlignments++;
		pGroundTruth->FlgPE2Err = true;
		NumErrPE2++;
		continue;
		}

	int64_t BaseScore;
	int BasesMatched;
	BasesMatched = ScoreAlignment(pSAMalign, pGroundTruth,&BaseScore);
	m_TotBaseAlignmentScore += BaseScore;
	m_ReadOverlapHistogram[(((int64_t)BasesMatched*100) + 50)/ PotentialAlgnBases]+=1;
	NumScoredAlignments++;
	}
gDiagnostics.DiagOut(eDLInfo, gszProcName, "Completed scoring alignments in: '%s'", pszAlignmentsFile);
gDiagnostics.DiagOut(eDLInfo, gszProcName, "There were a total of %u alignments scored containing %lld potential ground truth matching bases", NumScoredAlignments, m_TotNumPotentialAlignBases);
gDiagnostics.DiagOut(eDLInfo, gszProcName, "%u alignments classified as misaligned due to: chrom: %u, strand: %u, PE mismatch: %u ", NumErrChroms+NumErrStrands+ NumErrPE2,NumErrChroms,NumErrStrands, NumErrPE2);
gDiagnostics.DiagOut(eDLInfo, gszProcName,"Bases loci correct: %d, misaligned: %d, unaligned: %d", m_NumBasesLociCorrect, m_NumBasesLociIncorrect, m_NumBasesLociUnclaimed);
double BaseScoreRate = (double)m_TotBaseAlignmentScore/m_MaxBaseAlignmentScore;
double AlignBaseScoreRate = (double)m_TotBaseAlignmentScore/AlignedPotentialScore;
gDiagnostics.DiagOut(eDLInfo, gszProcName,"Alignment base score rate relative to background ground truth simulation bases is %1.6f", BaseScoreRate);
gDiagnostics.DiagOut(eDLInfo, gszProcName, "Alignment base score rate relative to actually aligned read ground truth simulation bases is %1.6f", AlignBaseScoreRate);

gDiagnostics.DiagOut(eDLInfo, gszProcName, "Aligned reads with ground truth bases percentage overlap histogram:");
for(int Idx = 0; Idx < 101; Idx++)
	gDiagnostics.DiagOutMsgOnly(eDLInfo, "%u%%: %u", Idx, m_ReadOverlapHistogram[Idx]);
delete []pszLine;
delete pSAMalign;
Reset();
return(Rslt);
}



// CIGAR operations recognised are:
// eCOPMatch,			'M' aligned but could be either matching or mismatching, consumes query and reference sequence, advances alignment loci
// eCOPInsert,			'I'  insertion relative to target - consumes query sequence only, does not advance alignment loci
// eCOPDelete,			'D'  deletion relative to target - consumes reference sequence only, advances alignment loci
// eCOPSkipRegion,		'N'  skipped region relative to target - intron?  - consumes reference sequence only, advances alignment loci
// eCOPSoftClip,		'S'  soft clipping - consumes query sequence only
// eCOPHardClip,		'H'  hard clipping - consumes neither query or reference sequence
// eCOPPadding,			'P' padding (silent deletion from padded reference) consumes neither query or reference sequence
// eCOPALignMatch,		'=' aligned as exactly matching, consumes query and reference sequence, advances alignment loci
// eCOPALignMismatch,	'X' aligned but as a mismatch, consumes query and reference sequence, advances alignment loci
// eCOPUnrecognised  unrecognised CIGAR operation
//

int									// returns index of next ClaimOpIdx, 0 if no more claim Ops, -1 if errors 
CBenchmark::GetClaimOps(tsBAMalign* pAlignment,	// claim CIGARs for this alignment
	int ClaimOpIdx,					// get ClaimNumOps, ClaimOpType for this CIGAR operation, starts from 1
	int* pNumOps,					// returned number of operations required
	etCIGAROpType* pTypeOp)			// type of operation
{
	uint32_t Oper;
	int MaxClaimOpIdx;
	if (pNumOps != NULL)
		*pNumOps = 0;
	if (pTypeOp != NULL)
		*pTypeOp = eCOPUnrecognised;
	if (pNumOps == NULL || pTypeOp == NULL || pAlignment == NULL)
		return(-1);
	if(ClaimOpIdx < 1)
		ClaimOpIdx = 1;
	MaxClaimOpIdx = pAlignment->flag_nc & 0x0ff;
	if (ClaimOpIdx > MaxClaimOpIdx)
		return(0);

	Oper = pAlignment->cigar[ClaimOpIdx-1];
	*pTypeOp = (etCIGAROpType)(Oper & 0x0f);
	*pNumOps = (Oper >> 4) & 0x0fffffff;
	return(++ClaimOpIdx);
}

int										    // number of chars consumed in CIGAR, 0 if none consumed, -1 if errors
CBenchmark::ParseExpCIGAR(char* pszCIGAR,	// parse starting from this CIGAR char
	int* pNumOps,						    // returned number of operations required
	etCIGAROpType* pTypeOp)					// type of operation
{
	int NumChrs;
	int NumOps;
	char Chr;

	if (pNumOps != NULL)
		*pNumOps = 0;
	if (pTypeOp != NULL)
		*pTypeOp = eCOPUnrecognised;

	if (pszCIGAR == NULL || pNumOps == NULL || pTypeOp == NULL)
		return(-1);

	if (*pszCIGAR == '\0')
		return(0);

	NumChrs = 0;
	NumOps = 0;
	while ((Chr = *pszCIGAR++) != '\0')
	{
		NumChrs += 1;
		if (Chr >= '0' && Chr <= '9')
		{
			if (NumOps > 100000000)	// hard limit, surely can't have alignments spanning into the Gbp
				return(-1);
			NumOps *= 10;
			NumOps += Chr - '0';
			continue;
		}

		switch (Chr) {
		case '=':			// '=' aligned as exactly matching, consumes both query and reference sequence
			*pTypeOp = eCOPALignMatch;
			break;

		case 'X': case 'x': // 'X' aligned but as a mismatch, consumes both query and reference sequence
			*pTypeOp = eCOPALignMismatch;
			break;

		case 'M': case 'm': // 'M' aligned but could be either matching or mismatching, consumes query and reference sequence
			*pTypeOp = eCOPMatch;
			break;

		case 'I': case 'i':  // 'I' insertion relative to target - consumes query sequence only
			*pTypeOp = eCOPInsert;
			break;

		case 'D': case 'd':	 // 'D'  deletion relative to target - consumes reference sequence only
			*pTypeOp = eCOPDelete;
			break;

		case 'N': case 'n':	 //	'N'  skipped region relative to target - intron?  - consumes reference sequence only
			*pTypeOp = eCOPSkipRegion;
			break;

		case 'S': case 's':  //	'S'  soft clipping - consumes query sequence only
			*pTypeOp = eCOPSoftClip;
			break;

		case 'H': case 'h':	 //	'H'  hard clipping - consumes neither query or reference sequence
			*pTypeOp = eCOPHardClip;
			break;

		case 'P': case 'p':	 //	'P' padding (silent deletion from padded reference) consumes neither query or reference sequence
			*pTypeOp = eCOPPadding;
			break;

		default:
			*pTypeOp = eCOPUnrecognised;
			return(-1);
		}

		*pNumOps = NumOps;
		return(NumChrs);
	}
	return(-1);
}

//Op BAM Description                                          Consumes query    Consumes reference 
//M  0   alignment match(can be a sequence match or mismatch) yes               yes 
//I  1   insertion to the reference                           yes               no 
//D  2   deletion from the reference                          no                yes 
//N  3   skipped region from the reference                    no                yes 
//S  4   soft clipping(clipped sequences present in SEQ)      yes               no 
//H  5   hard clipping(clipped sequences NOT present in SEQ)  no                no 
//P  6   padding(silent deletion from padded reference)       no                no 
//=  7   sequence match                                       yes               yes
//X  8   sequence mismatch                                    yes               yes

// Note that read+target are not adjusted if OpType is for simple match/mismatches!
int
CBenchmark::AdjustAlignment(etCIGAROpType OpType,	// alignment operator type
						int OpCnt,					// number of times operator to be applied
						int* pReadOfs,				// current read offset
						int* pAlignLoci)			// read offset corresponds to this alignment loci
{
switch (OpType) {
	case eCOPMatch:			// 'M' aligned but could be either matching or mismatching, consumes query and reference sequence
	case eCOPALignMatch:	// '=' aligned as exactly matching, consumes both query and reference sequence
	case eCOPALignMismatch:	// 'X' aligned but as a mismatch, consumes both query and reference sequence
		break;				// no change to either target Loci and ReadOfs

	case eCOPInsert:			//'I'  insertion relative to target - consumes query sequence only
		*pReadOfs += OpCnt;		// no changes to target loci
		break;

	case eCOPSkipRegion:    	    // 'N'  skipped region relative to target - intron?  - consumes reference sequence only
	case eCOPDelete:				// 'D'  deletion relative to target - consumes reference sequence only
		*pAlignLoci += OpCnt;		// no change to ExpectReadOfs
		break;

	case eCOPSoftClip:				//	'S'  soft clipping - consumes query sequence only
		*pReadOfs += OpCnt;
		break;

	case eCOPHardClip:				//	'H'  hard clipping - consumes neither query or reference sequence
	case eCOPPadding:				//	'P' padding (silent deletion from padded reference) consumes neither query or reference sequence
		break;

	default:		// unrecognised CIGAR operation
		return(-1);
	}
return(OpCnt);
}

int		// returned number of potential base matches in CIGAR
CBenchmark::ReadMatchScore(int ReadLen,			// read sequence length
	char* pszCIGAR,		// '\0' terminated CIGAR with alignment profile
	int MatchPts,		// pts for each base in alignment profile which is a potential match - could be Pts for match, mismatch, or unaligned
	int64_t* pScore)	// returned number of matches * MatchPts
{
int NumOpTypes;
int NumMatches;
uint32_t OpTypes[cMaxBAMCigarOps];
etCIGAROpType OpType;
if(pScore != NULL)
	*pScore = 0;

// total the number of matches in the CIGAR profile
if ((NumOpTypes = ParseObsCIGAR(pszCIGAR, ReadLen, cMaxBAMCigarOps, OpTypes)) < 1)
		return(0);
NumMatches = 0;
for(int Idx = 0; Idx < NumOpTypes; Idx++)
	{
	OpType = (etCIGAROpType)(OpTypes[Idx] & 0x0f);
	if(OpType == eCOPMatch || OpType == eCOPALignMatch || OpType == eCOPALignMismatch)
		NumMatches += (OpTypes[Idx] >> 4 )& 0x0fffffff;
	}
if(pScore != NULL)
	*pScore = (int64_t)NumMatches * MatchPts;
return(NumMatches);
}

//Op BAM Description                                          Consumes query    Consumes reference 
//M  0   alignment match(can be a sequence match or mismatch) yes               yes 
//I  1   insertion to the reference                           yes               no 
//D  2   deletion from the reference                          no                yes 
//N  3   skipped region from the reference                    no                yes 
//S  4   soft clipping(clipped sequences present in SEQ)      yes               no 
//H  5   hard clipping(clipped sequences NOT present in SEQ)  no                no 
//P  6   padding(silent deletion from padded reference)       no                no 
//=  7   sequence match                                       yes               yes
//X  8   sequence mismatch                                    yes               yes
int				// returns number of correctly aligned base matches
CBenchmark::ScoreAlignment(tsBAMalign* pAlignment,			// claimed alignment
							tsBMGroundTruth *pGroundTruth,	// ground truth for this alignment
							int64_t* pBaseScore)			// returned base alignment score
{
int NumGroundTruthOpTypes;
int NumClaimOpTypes;
int GroundTruthNumOps;
int ClaimNumOps;
uint32_t GroundTruthOps[cMaxBAMCigarOps];

if(pAlignment == NULL || pGroundTruth == NULL || pBaseScore == NULL)
	return(0);
*pBaseScore = 0;

// parse ground truth CIGAR into packed GroundTruthOps
if ((NumGroundTruthOpTypes = ParseObsCIGAR((char*)&pGroundTruth->NameChromCIGAR[pGroundTruth->NameLen + pGroundTruth->ChromNameLen + 2], pGroundTruth->ReadLen, cMaxBAMCigarOps, GroundTruthOps)) < 1)
	return(0);

// if it a silently trimmed claimed sequence then can't really infer which end (could be both) was trimmed
if (pGroundTruth->ReadLen != pAlignment->l_seq)		// check if silently trimmed
	{
	int HardTrimLen = 0;
	int SilentTrimBaseMatches;
	if(eCOPHardClip == (etCIGAROpType)(pAlignment->cigar[0] & 0x0f))
		HardTrimLen = (pAlignment->cigar[0] >> 4) & 0x0fffffff;
	if (eCOPHardClip == (etCIGAROpType)(pAlignment->cigar[(pAlignment->flag_nc & 0x0ff) - 1] & 0x0f))
		HardTrimLen += (pAlignment->cigar[(pAlignment->flag_nc & 0x0ff)-1] >> 4) & 0x0fffffff;
	if (pGroundTruth->ReadLen != pAlignment->l_seq + HardTrimLen) // if not equal then alignment read has been silently trimmed or incorrectly reported
		{
		// does claimed loci range fall within the ground truth loci range, if so then score with m_SilentTrimAlignBasePts otherwise score with m_MisalignedBasePts
		int SilentTrimScore;
		int GroundTruthEnd;
		int64_t BaseScore;
		GroundTruthEnd = pGroundTruth->StartLoci + RefSeqConsumedLen(NumGroundTruthOpTypes,GroundTruthOps) - 1;
		SilentTrimScore = m_SilentTrimAlignBasePts;
		if((pAlignment->pos + 15) >= (int)pGroundTruth->StartLoci && ((int)pAlignment->end - 15) <= GroundTruthEnd)	// allowing a 15bp latitude!
			SilentTrimBaseMatches = ReadMatchScore(pGroundTruth->ReadLen, (char*)&pGroundTruth->NameChromCIGAR[pGroundTruth->NameLen + pGroundTruth->ChromNameLen + 2], m_SilentTrimAlignBasePts, &BaseScore);
		m_NumBasesLociCorrect += SilentTrimBaseMatches;	// have to assume these did match - Ugh, Ugh, and Ugh - because what is the alternative, can't just ignore them
		*pBaseScore = BaseScore;
		return(SilentTrimBaseMatches);
		}
	}


NumClaimOpTypes = pAlignment->flag_nc & 0x0ff;

int ClaimedReadOfs = 0;
int ExpectReadOfs = 0;;
int NumLociCorrect = 0;
int NumLociIncorrect = 0;
int NumLociUnclaimed = 0;
int ClaimOpIdx;
int GroundTruthOpIdx;
int ClaimedLoci;
int GroundTruthLoci;
etCIGAROpType GroundTruthOpType;
etCIGAROpType ClaimOpType;
ClaimOpIdx = 0;
ClaimedLoci = pAlignment->pos;
GroundTruthOpIdx = 0;
GroundTruthLoci = pGroundTruth->StartLoci;
GroundTruthNumOps = 0;
ClaimNumOps = 0;
int ClaimedSeqIdx = 0;
int GroundTruthSeqIdx = 0;
while(GroundTruthSeqIdx < pGroundTruth->ReadLen)
	{
	if (!GroundTruthNumOps)
		{
		if(GroundTruthOpIdx == NumGroundTruthOpTypes)
			break;
		GroundTruthOpType = (etCIGAROpType)(GroundTruthOps[GroundTruthOpIdx] & 0x0f);
		if(GroundTruthOpType == eCOPALignMatch || GroundTruthOpType == eCOPALignMismatch || GroundTruthOpType == eCOPMatch) // simpifying exacts and mismatches
			GroundTruthOpType = eCOPMatch;												  // as if just matching
		else
			if (GroundTruthOpType == eCOPSkipRegion)								    // treat region skipping same
				GroundTruthOpType = eCOPDelete;											// as deletion - consumes target only
		GroundTruthNumOps = (GroundTruthOps[GroundTruthOpIdx++] >> 4) & 0x0fffffff;
		}
	
	if (!ClaimNumOps)
		{
		if (ClaimOpIdx == NumClaimOpTypes)
			break;
		ClaimOpType = (etCIGAROpType)(pAlignment->cigar[ClaimOpIdx] & 0x0f);
		if (ClaimOpType == eCOPALignMatch || ClaimOpType == eCOPALignMismatch || ClaimOpType == eCOPMatch) // simpifying exacts and mismatches
			ClaimOpType = eCOPMatch;										   // as if just matching
		else
			if (ClaimOpType == eCOPSkipRegion)								   // treat region skipping same
				ClaimOpType = eCOPDelete;									   // as deletion - consumes target only
			else
				if (ClaimOpType == eCOPHardClip)								// treat hard clipping same
					ClaimOpType = eCOPSoftClip;									// as soft clipping
		ClaimNumOps = (pAlignment->cigar[ClaimOpIdx++] >> 4) & 0x0fffffff;
		}

	switch(ClaimOpType) {
		case eCOPSoftClip:						//	'S'  claim sequence has been soft clipped
			switch (GroundTruthOpType) {
				case eCOPMatch:					// 'M' aligned but could be either matching or mismatching
					NumLociUnclaimed++;			// not claimed even though ground truth had a match		
					GroundTruthSeqIdx++;		// ground truth is consuming both base and loci
					GroundTruthLoci++;
					break;

				case eCOPSoftClip:				// 'S'  soft clipping
					GroundTruthSeqIdx++;		// ground truth is consuming  both base and loci
					GroundTruthLoci++;			// remember that ground truth loci is from start of read
					break;

				case eCOPDelete:				// ground truth says that there is a deletion in the target
					GroundTruthLoci++;			// increment truth loci but ground truth is not consuming a read base
					break;

				case eCOPPadding:				// silent padding
					break;						// ground truth is not consuming this base and no change in the ground truth loci

				case eCOPInsert:				// 'I' insertion
					GroundTruthSeqIdx++;		// ground truth is consuming this base but no change in the ground truth loci
					break;
				}
			ClaimedSeqIdx++;					// claim is only consuming sequence and not loci
			break;

		case eCOPDelete:						//	'D'  alignment skips over a target subsequence to a new loci in target
			switch (GroundTruthOpType) {
				case eCOPMatch:					// 'M' aligned but could be either matching or mismatching, consumes query and reference sequence
					NumLociUnclaimed++;			// not claimed even though ground truth had a match		
					GroundTruthSeqIdx++;		// ground truth is consuming a base
					GroundTruthLoci++;
					break;

				case eCOPSoftClip:				// 'S'  soft clipping
					GroundTruthSeqIdx++;		// ground truth is consuming this base
					GroundTruthLoci++;			// and loci
					break;

				case eCOPDelete:				// ground truth says that there is a deletion in the read sequence
					GroundTruthLoci += 1;		// increment truth loci but ground truth is not consuming a read base
					break;

				case eCOPPadding:				// silent padding
					break;						// ground truth is not consuming this base and no change in the ground truth loci

				case eCOPInsert:				// 'I' insertion
					GroundTruthSeqIdx++;		// ground truth is consuming this base but no change in the ground truth loci
					break;
				}
			ClaimedLoci++;
			break;

		case eCOPPadding:				// silent padding
			switch (GroundTruthOpType) {
				case eCOPMatch:					// 'M' aligned but could be either matching or mismatching, consumes query and reference sequence
					NumLociUnclaimed++;			// not claimed even though ground truth had a match		
					GroundTruthSeqIdx++;		// ground truth is consuming a base
					GroundTruthLoci++;
					break;

				case eCOPSoftClip:				// 'S'  soft clipping
					GroundTruthSeqIdx++;		// ground truth is consuming this base
					GroundTruthLoci++;			// and loci
					break;

				case eCOPDelete:				// ground truth says that there is a deletion in the read sequence
					GroundTruthLoci++;			// increment truth loci but ground truth is not consuming a read base
					break;

				case eCOPPadding:				// silent padding
					break;						// ground truth is not consuming this base and no change in the ground truth loci

				case eCOPInsert:				// 'I' insertion
					GroundTruthSeqIdx++;		// ground truth is consuming this base but no change in the ground truth loci
					break;
				}
			break;

		case eCOPInsert:				// 'I' insertion
			switch (GroundTruthOpType) {
				case eCOPMatch:					// 'M' aligned but could be either matching or mismatching, consumes query and reference sequence
					NumLociUnclaimed++;			// not claimed even though ground truth had a match		
					GroundTruthSeqIdx++;		// ground truth is consuming a base
					GroundTruthLoci++;
					break;

				case eCOPSoftClip:				// 'S'  soft clipping
					GroundTruthSeqIdx++;		// ground truth is consuming this base and loci
					GroundTruthLoci++;
					break;

				case eCOPDelete:				// ground truth says that there is a deletion in the read sequence
					GroundTruthLoci += 1;		// increment truth loci but ground truth is not consuming a read base
					break;

				case eCOPPadding:				// silent padding
					break;						// ground truth is not consuming this base and no change in the ground truth loci

				case eCOPInsert:				// 'I' insertion
					GroundTruthSeqIdx++;		// ground truth is consuming this base but no change in the ground truth loci
					break;
				}
			ClaimedSeqIdx++;
			break;

		case eCOPMatch:			// 'M' aligned but could be either matching or mismatching, consumes query and reference sequence
			switch (GroundTruthOpType) {
				case eCOPMatch:					// 'M' aligned but could be either matching or mismatching, consumes query and reference sequence
					if(ClaimedLoci == GroundTruthLoci)
						NumLociCorrect++;
					else
						NumLociIncorrect++;
					GroundTruthSeqIdx++;		// ground truth is consuming a base
					GroundTruthLoci++;
					break;

				case eCOPSoftClip:				// 'S'  soft clipping
					NumLociIncorrect++;
					GroundTruthSeqIdx++;		// ground truth is consuming this base
					GroundTruthLoci++;
					break;

				case eCOPDelete:				// ground truth says that there is a deletion in the read sequence
					NumLociIncorrect++;
					GroundTruthLoci++;		// increment truth loci but ground truth is not consuming a read base
					break;

				case eCOPPadding:				// silent padding
					NumLociIncorrect++;
					break;						// ground truth is not consuming this base and no change in the ground truth loci

				case eCOPInsert:				// 'I' insertion
					NumLociIncorrect++;
					GroundTruthSeqIdx++;		// ground truth is consuming this base but no change in the ground truth loci
					break;
				}
			ClaimedLoci++;
			ClaimedSeqIdx++;
			break;
		}	

	GroundTruthNumOps -=1;
	ClaimNumOps-=1;
	}
m_NumBasesLociCorrect += NumLociCorrect;
m_NumBasesLociIncorrect += NumLociIncorrect;
m_NumBasesLociUnclaimed += NumLociUnclaimed;
*pBaseScore = ((int64_t)NumLociCorrect * m_AlignedBasePts) + ((int64_t)NumLociUnclaimed * m_UnalignedBasePts) + ((int64_t)NumLociIncorrect * m_MisalignedBasePts);
return(NumLociCorrect);
}


//Op BAM Description                                          Consumes query    Consumes reference 
//M  0   alignment match(can be a sequence match or mismatch) yes               yes 
//I  1   insertion to the reference                           yes               no 
//D  2   deletion from the reference                          no                yes 
//N  3   skipped region from the reference                    no                yes 
//S  4   soft clipping(clipped sequences present in SEQ)      yes               no 
//H  5   hard clipping(clipped sequences NOT present in SEQ)  no                no 
//P  6   padding(silent deletion from padded reference)       no                no 
//=  7   sequence match                                       yes               yes
//X  8   sequence mismatch                                    yes               yes
int		// for given CIGARs returns total number of bases that these CIGARs consume on aligned to reference sequence
CBenchmark::RefSeqConsumedLen(int NumPackedOps,	  // there are this many packed CIGAR ops in pPackedOps
			uint32_t* pPackedOps) // packed CIGAR ops into this array
{
etCIGAROpType CIGARop;
int RefConsumLen;
int Idx;
if(NumPackedOps < 1 || pPackedOps == NULL || *pPackedOps == 0)
	return(-1);

RefConsumLen = 0;
for(Idx = 0; Idx < NumPackedOps; Idx++, pPackedOps++)
	{
	CIGARop = (etCIGAROpType )(*pPackedOps & 0x0f);
	switch (CIGARop) {
		case eCOPMatch:				//	'M' aligned but could be either matching or mismatching, consumes query and reference sequence
		case eCOPALignMatch:		//  '=' aligned as exactly matching , consumes both query and reference sequence
		case eCOPALignMismatch:		//	'X' aligned but as a mismatch, consumes query and reference sequence
		case eCOPDelete:			//	'D'  deletion relative to target - consumes reference sequence only
		case eCOPSkipRegion:		//	'N'  skipped region relative to target - intron?  - consumes reference sequence only
			RefConsumLen += (*pPackedOps >> 4) & 0x0fffffff;
			continue;

		case eCOPUnrecognised:		// unrecognised CIGAR operation
			return(-1);

		default:					// reference sequence not consumed
			continue;
		}
	}
return(RefConsumLen);
}

int								 // number of accepted CIGAR operations, 0 if unable to accept 
CBenchmark::ParseObsCIGAR(char* pszCIGAR, // CIGAR operations
	int ReadLen,				// length ('=','X','M','D') must be equal to ReadLen
	int MaxCIGAROps,			// pPackedOps can hold at most this many Ops
	uint32_t* pPackedOps)		// pack CIGAR ops into this array, can be NULL if just validating the CIGAR
{
	char Chr;
	int NumChrs;
	uint32_t NumOps;
	int CIGAROps;
	int CIGARReadLen;
	etCIGAROpType CIGARop;

	if (pszCIGAR == NULL || *pszCIGAR == '\0' || ReadLen < cSRMinReadLen || ReadLen > cSRMaxReadLen)
		return(false);

	CIGARReadLen = 0;
	NumOps = 0;
	NumChrs = 0;
	CIGAROps = 0;
	if (pPackedOps != NULL)
		*pPackedOps = 0;

	while ((Chr = *pszCIGAR++) != '\0')
	{
		NumChrs += 1;
		if (NumChrs > cSRMaxCIGARLen)
			return(0);
		if (Chr >= '0' && Chr <= '9')
		{
			if (NumOps > (uint32_t)500000)
				return(0);
			NumOps *= 10;
			NumOps += Chr - '0';
			continue;
		}
		
		if (CIGAROps == MaxCIGAROps)
			return(0);
		switch (Chr) {
		case '=':			// '=' aligned as exactly matching, consumes both query and reference sequence
			CIGARop = eCOPALignMatch;
			CIGARReadLen += NumOps;
			if (CIGARReadLen > ReadLen)
				return(0);
			break;

		case 'X': case 'x': // 'X' aligned but as a mismatch, consumes both query and reference sequence
			CIGARop = eCOPALignMismatch;
			CIGARReadLen += NumOps;
			if (CIGARReadLen > ReadLen)
				return(0);
			break;

		case 'M': case 'm': //	'M' aligned but could be either matching or mismatching, consumes query and reference sequence
			CIGARop = eCOPMatch;
			CIGARReadLen += NumOps;
			if (CIGARReadLen > ReadLen)
				return(0);
			break;

		case 'I': case 'i':  //'I'  insertion relative to target - consumes query sequence only
			CIGARop = eCOPInsert;
			CIGARReadLen += NumOps;
			if (CIGARReadLen > ReadLen)
				return(0);
			break;

		case 'D': case 'd':	 // 'D'  deletion relative to target - consumes reference sequence only
			CIGARop = eCOPDelete;
			break;

		case 'N': case 'n':	 // skipped region in target (intron?)
			CIGARop = eCOPSkipRegion;
			break;

		case 'P': case 'p':	 // 'P' padding (silent deletion from padded reference) consumes neither query or reference sequence
			CIGARop = eCOPPadding;
			break;

		case 'S': case 's':	 // 'S'  soft clipping - consumes query sequence only
			CIGARop = eCOPSoftClip;
			CIGARReadLen += NumOps;
			if (CIGARReadLen > ReadLen)
				return(0);
			break;

		case 'H': case 'h':	 // 'H'  hard clipping - consumes neither query or reference sequence
			CIGARop = eCOPHardClip;
			break;

		default:			// not accepting any other CIGAR operator
			return(0);
		}
		NumOps <<= 4;
		NumOps |= (uint32_t)CIGARop & 0x0f;
		if(pPackedOps != NULL)
			*pPackedOps++ = NumOps;
		CIGAROps += 1;
		NumOps = 0;
	}
	return(CIGARReadLen == ReadLen ? CIGAROps : 0);
}

int				// returns number of alignment inplace coalesced CIGAR ops
CBenchmark::CoalesceCIGARs(tsBAMalign* pBAM)	// in: alignment to coalesce out: updated alignment with coalesced CIGAR ops
{
int NumTypeOps;
int OpIdx;
int NumCoalesced;
etCIGAROpType CIGARop;
uint32_t NumOps;
etCIGAROpType NxtCIGARop;
uint32_t *pPackedOps;
uint32_t* pNxtPackedOps;
if(pBAM == NULL)
	return(-1);
// a quick validation, number of bytes required for CIGARs should be 4 * number of CIGAR ops
NumTypeOps = pBAM->flag_nc & 0x0ff;
if(NumTypeOps != (pBAM->NumCigarBytes / sizeof(uint32_t)))
	return(-1);
if(NumTypeOps == 1)
	return(1);
NumCoalesced = 1;
pPackedOps = pBAM->cigar;
pNxtPackedOps = pPackedOps + 1;
for(OpIdx = 1; OpIdx < NumTypeOps; OpIdx++, pNxtPackedOps++)
	{
	NumOps = (*pPackedOps >> 4) & 0x0fffffff;
	CIGARop = (etCIGAROpType)(*pPackedOps & 0x0f);
	if (CIGARop == eCOPALignMatch || CIGARop == eCOPALignMismatch) // alignment exact match/mismatches are coalesced into a single possible match
		CIGARop = eCOPMatch;
	*pPackedOps = (NumOps << 4) | (uint32_t)CIGARop;
	NxtCIGARop = (etCIGAROpType)(*pNxtPackedOps & 0x0f);
	if (NxtCIGARop == eCOPALignMatch || NxtCIGARop == eCOPALignMismatch)
		NxtCIGARop = eCOPMatch;
	if(NxtCIGARop == CIGARop)
		{
		NumOps += (*pNxtPackedOps >> 4) & 0x0fffffff;
		*pPackedOps = (NumOps <<= 4) | (uint32_t)CIGARop;
		}
	else
		{
		*(++pPackedOps) = (*pNxtPackedOps & 0xfffffff0) | (uint32_t)NxtCIGARop;
		NumCoalesced += 1;
		}
	}
pBAM->flag_nc = (pBAM->flag_nc & 0xfffff000) | NumCoalesced;
return(NumCoalesced);
}

int
CBenchmark::InitObsErrProfile(char* pszInProfFile)	// read from this observed alignment error profiles file (as generated from actual SAM aligments)
{
	int Rslt;
	int NumProcessed;
	int NumFields;
	int CIGARID;
	char PE1Strand;		// '+' sense or watson, '-' antisense or crick
	int PE1NumCIGAROps;
	int PE1NumErrProfOps;
	int ReadLen;
	int PEInsertSize;	// PE reads insert size
	char PE2Strand;		// '+' sense or watson, '-' antisense or crick
	int PE2NumCIGAROps;
	int PE2NumErrProfOps;
	char* pszTxt;

	uint32_t PE1ObsCIGAROps[cMaxBAMCigarOps];
	uint32_t PE2ObsCIGAROps[cMaxBAMCigarOps];

	uint32_t PE1ErrProfOps[cMaxBAMCigarOps * 3];
	uint32_t PE2ErrProfOps[cMaxBAMCigarOps * 3];
	
	tsBMPackedCIGARs* pCurCIGAR;
	if (pszInProfFile == NULL || pszInProfFile[0] == '\0')
	{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "InitObsErrProfile: Observed alignment error profile file not specified");
		Reset();
		return(eBSFerrParams);
	}

	gDiagnostics.DiagOut(eDLInfo, gszProcName, "InitObsErrProfile: Loading observed alignment error profiles from file: %s", pszInProfFile);

	m_UsedObsErrProfMem = 0;
	m_NumObsErrProfs = 0;

	m_AllocdObsErrProfMem = (size_t)cSRAllocdObsErrProfMemChunk;
#ifdef _WIN32
	m_pObsErrProfiles = (UINT8*)malloc(m_AllocdObsErrProfMem);
	if (m_pObsErrProfiles == NULL)
	{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "InitObsErrProfile: Memory allocation of %lld bytes failed", (INT64)m_AllocdObsErrProfMem);
		m_AllocdObsErrProfMem = 0;
		Reset();
		return(eBSFerrMem);
	}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	m_pObsErrProfiles = (UINT8*)mmap(NULL, (size_t)m_AllocdObsErrProfMem, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	if (m_pObsErrProfiles == MAP_FAILED)
	{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "InitObsErrProfile: Memory allocation of %lld bytes through mmap()  failed", (INT64)m_AllocdObsErrProfMem, strerror(errno));
		m_pObsErrProfiles = NULL;
		m_AllocdObsErrProfMem = 0;
		Reset();
		return(eBSFerrMem);
	}
#endif

	if ((m_pObsCIGARProfFile = new CCSVFile) == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "InitObsErrProfile: Unable to instantiate CCSVfile");
		Reset();
		return(eBSFerrObj);
		}

	if ((Rslt = m_pObsCIGARProfFile->Open(pszInProfFile)) != eBSFSuccess)
		{
		while (m_pObsCIGARProfFile->NumErrMsgs())
			gDiagnostics.DiagOut(eDLFatal, gszProcName, m_pObsCIGARProfFile->GetErrMsg());
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "InitObsErrProfile: Unable to open file: %s", pszInProfFile);
		Reset();
		return(Rslt);
		}

	int MinExpFields;
	NumProcessed = 0;
	MinExpFields = m_bPEReads ? 9 : 5;	// if SE then need at least 5 fields; if PE then need at least 9
	while ((Rslt = m_pObsCIGARProfFile->NextLine()) > 0)	// onto next line containing fields
		{
		NumFields = m_pObsCIGARProfFile->GetCurFields();
		
		if (NumFields < MinExpFields)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "InitObsErrProfile: %s contains % fields, expected at least %d", pszInProfFile, NumFields,MinExpFields);
			Reset();
			return(eBSFerrParams);
			}
		if (!NumProcessed && m_pObsCIGARProfFile->IsLikelyHeaderLine())
			continue;
		NumProcessed += 1;
		PE1Strand = '\0';
		ReadLen = 0;
		PE2Strand = '\0';
		PEInsertSize = 0;
		PE1NumCIGAROps = 0;
		PE1NumErrProfOps = 0;
		PE2NumCIGAROps = 0;
		PE2NumErrProfOps = 0;

		m_pObsCIGARProfFile->GetInt(1, &CIGARID);
		m_pObsCIGARProfFile->GetInt(2, &ReadLen);
		if (ReadLen < cSRMinReadLen || ReadLen > cSRMaxReadLen)
			continue;
		m_pObsCIGARProfFile->GetChar(3, &PE1Strand);	// PE1 alignment strand
		m_pObsCIGARProfFile->GetText(4, &pszTxt);		// PE1 CIGAR
		if ((PE1NumCIGAROps = ParseObsCIGAR(pszTxt, ReadLen, cMaxBAMCigarOps, &PE1ObsCIGAROps[0])) <= 0)
			continue;
		m_pObsCIGARProfFile->GetText(5, &pszTxt);		// PE1 error profile
		if ((PE1NumErrProfOps = ParseObsCIGAR(pszTxt, ReadLen, cMaxBAMCigarOps*3, &PE1ErrProfOps[0])) <= 0)
			continue;

		if (m_bPEReads)
			{
			m_pObsCIGARProfFile->GetInt(6, &PEInsertSize);	// PE1-PE2 insert size
			if (abs(PEInsertSize) <= ReadLen/4)
				continue;
			m_pObsCIGARProfFile->GetChar(7, &PE2Strand);	// PE2 alignment strand
			m_pObsCIGARProfFile->GetText(8, &pszTxt);		// PE2 CIGAR
			if ((PE2NumCIGAROps = ParseObsCIGAR(pszTxt, ReadLen, cMaxBAMCigarOps, &PE2ObsCIGAROps[0])) == 0)
				continue;
			m_pObsCIGARProfFile->GetText(9, &pszTxt);		// PE2 CIGAR
			if ((PE2NumErrProfOps = ParseObsCIGAR(pszTxt, ReadLen, cMaxBAMCigarOps*3, &PE2ErrProfOps[0])) == 0)
				continue;
			}

		

		if (((sizeof(tsBMPackedCIGARs) + 10000) >= (m_AllocdObsErrProfMem - m_UsedObsErrProfMem)))
			{
			size_t memreq = m_AllocdObsErrProfMem + (size_t)cSRAllocdObsErrProfMemChunk;
			UINT8* pTmp;
#ifdef _WIN32
			pTmp = (UINT8*)realloc(m_pObsErrProfiles, memreq);
#else
			pTmp = (UINT8*)mremap(m_pObsErrProfiles, m_AllocdObsErrProfMem, memreq, MREMAP_MAYMOVE);
			if (pTmp == MAP_FAILED)
				pTmp = NULL;
#endif
			if (pTmp == NULL)
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "InitObsErrProfile: Memory re-allocation to %lld bytes - %s", (INT64)(memreq), strerror(errno));
				Reset();
				return(eBSFerrMem);
				}
			m_pObsErrProfiles = pTmp;
			m_AllocdObsErrProfMem = memreq;
			}
		pCurCIGAR = (tsBMPackedCIGARs*)&m_pObsErrProfiles[m_UsedObsErrProfMem];
		pCurCIGAR->ID = CIGARID;
		pCurCIGAR->ReadLen = (uint16_t)ReadLen;
		pCurCIGAR->FlgPE1Strand = PE1Strand == '-' ? 1 : 0;
		pCurCIGAR->FlgPE2Strand = PE2Strand == '-' ? 1 : 0;
		pCurCIGAR->PE1NumCIGAROps = (uint8_t)PE1NumCIGAROps;
		pCurCIGAR->PE1NumErrProfOps = (uint8_t)PE1NumErrProfOps;
		pCurCIGAR->PE2NumCIGAROps = (uint8_t)PE2NumCIGAROps;
		pCurCIGAR->PE2NumErrProfOps = (uint8_t)PE2NumErrProfOps;
		pCurCIGAR->PEInsertSize = PEInsertSize;
		memcpy(pCurCIGAR->CIGAROpsErrProfOps, PE1ObsCIGAROps, pCurCIGAR->PE1NumCIGAROps * sizeof(uint32_t));
		memcpy(&pCurCIGAR->CIGAROpsErrProfOps[PE1NumCIGAROps], PE1ErrProfOps, (int64_t)PE1NumErrProfOps * sizeof(uint32_t));

		if(pCurCIGAR->PE2NumCIGAROps)
			{
			memcpy(&pCurCIGAR->CIGAROpsErrProfOps[PE1NumCIGAROps + PE1NumErrProfOps], PE2ObsCIGAROps, PE2NumCIGAROps * sizeof(uint32_t));
			memcpy(&pCurCIGAR->CIGAROpsErrProfOps[PE1NumCIGAROps + PE1NumErrProfOps + PE2NumCIGAROps], PE2ErrProfOps, PE2NumErrProfOps * sizeof(uint32_t));
			}
		pCurCIGAR->Size = (uint16_t)(sizeof(tsBMPackedCIGARs) + ((PE1NumCIGAROps + PE1NumErrProfOps + PE2NumCIGAROps + PE2NumErrProfOps - 1) * (int)sizeof(uint32_t)));

		m_UsedObsErrProfMem += (size_t)pCurCIGAR->Size;
		m_NumObsErrProfs += 1;
		}
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "InitObsErrProfile: Loaded %d observed error profiles", m_NumObsErrProfs);
	delete m_pObsCIGARProfFile;
	m_pObsCIGARProfFile = NULL;
	return(eBSFSuccess);
}



int			// choosen loci on length proportional randomly selected chromosome, -1 if unable to choose a chromosome meeting critera
CBenchmark::ChooseRandStartLoci(int FragSize,		// ensure can generate this sized fragment starting at chrom/loci
								int* pSelChromID)	// returned chrom identifier
{
tsBMChromSeq * pChromSeq;
int64_t Attempts;
int RandChrom;
int TargPsn;
int IdxLo;
int IdxHi;

if(pSelChromID == NULL)
	return(-1);
if(FragSize <= 0)
	FragSize = cBMMaxFragSize;

	// first choose a chrom, longer chroms have a proportional higher chance of being choosen relative to shorter chromosomes
Attempts = (int64_t)m_NumChromSeqs * 1000000; // surely can choose a chrome long enough to accommodate FragSize within this many random choices!
do {
	RandChrom = (int)RG.IRandom(1, m_GenomeScaledLen);
	IdxLo = 0;
	IdxHi = m_NumChromSeqs - 1;
	do {
		TargPsn = (IdxLo + IdxHi) / 2;
		pChromSeq = &m_pChromSeqs[TargPsn];
		if (pChromSeq->ScaledStart <= RandChrom && (pChromSeq->ScaledStart + pChromSeq->ScaledLen) >= RandChrom)
			break;

		if (pChromSeq->ScaledStart > RandChrom)
			IdxHi = TargPsn - 1;
		else
			IdxLo = TargPsn + 1;
	} while (IdxHi >= IdxLo);
// have a putative chrom, ensure it is sufficently long for FragSize otherwise randomly choose another chrom
} while(Attempts-- && pChromSeq->Len < FragSize);
if(pChromSeq->Len < FragSize)
	return(-1);

// selected a chromosome, now randomly choose a start loci on that chromosome
*pSelChromID = pChromSeq->ChromID;
if(pChromSeq->Len == FragSize)
	return(0);

return((int)RG.IRandom(0, pChromSeq->Len - FragSize));
}

// SortReadPairs
// Used to sort read pairs by name ascending with first in pair coming before 2nd in pair
// Return value meaning
// < 0 The element pointed by p1 goes before the element pointed by p2
//	0  The element pointed by p1 is equivalent to the element pointed by p2
// >0 The element pointed by p1 goes after the element pointed by p2
//
// stricmp: -1 if s1 < s2
//           0 if s1 == s2
//           +1 if s1 > s2
// complication is where strlen(s1) != strlen(s2) as terminating '0' compares less than other chrs
// if s1 is "abc" and s2 is "ab" then stricmp(s1,s2) returns -1 whereas stricmp(s2,s1) returns +1
//
int
CBenchmark::SortReadPairs(const void* arg1, const void* arg2)
{
int Cmp;
tsBMObsCIGAR* pEl1 = *(tsBMObsCIGAR**)arg1;
tsBMObsCIGAR* pEl2 = *(tsBMObsCIGAR**)arg2;
//int PE1Len = strlen((char*)pEl1->NameCIGAR);
//int PE2Len = strlen((char*)pEl2->NameCIGAR);
//if (PE1Len != PE2Len)
//	return(PE1Len - PE2Len);
if((Cmp = stricmp((char *)pEl1->NameCIGARErrProfile, (char*)pEl2->NameCIGARErrProfile))!=0)
	return(Cmp);  // Note:  -1 if pEl1 shorter in length than pEl2 but matches up to strlen(pEl1). +1 if pEl1 longer in length than pEl2 but matches up to strlen(pEl2)
if(!(pEl1->Flags & 0x01 && pEl2->Flags & 0x01))	// if either not a PE then just return as matching read names
	return(0);
if(pEl1->Flags & 0x40 && pEl2->Flags & 0x080)	// pEl1 is first in pair, pEl2 is 2nd in pair?
	return(0);									
return(1);					// currently pEl1 is 2nd in pair, need it to be 1st in pair
}

int
CBenchmark::SortGroundTruths(const void* arg1, const void* arg2)
{
	int Cmp;
	tsBMGroundTruth* pEl1 = *(tsBMGroundTruth**)arg1;
	tsBMGroundTruth* pEl2 = *(tsBMGroundTruth**)arg2;
//	int PE1Len = strlen((char*)pEl1->NameChromCIGAR);
//	int PE2Len = strlen((char*)pEl2->NameChromCIGAR);
//	if(PE1Len != PE2Len)
//		return(PE1Len - PE2Len);
	if ((Cmp = stricmp((char*)pEl1->NameChromCIGAR, (char*)pEl2->NameChromCIGAR)) != 0)
		return(Cmp); // Note:  -1 if pEl1 shorter in length than pEl2 but matches up to strlen(pEl1). +1 if pEl1 longer in length than pEl2 but matches up to strlen(pEl2)
	if (!(pEl1->FlgPE2 || pEl2->FlgPE2))	// if either not a PE then just return as must be SE
		return(0);
	if (pEl2->FlgPE2)	// pEl1 is first in pair, pEl2 is 2nd in pair?
		return(0);
	return(1);					// currently pEl1 is 2nd in pair, need it to be 1st in pair
}



