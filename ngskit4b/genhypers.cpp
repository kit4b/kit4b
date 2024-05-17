// genhypers.cpp : Defines the entry point for the console application.
// Locates all hyper elements according to parameterised filtering critera
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
#include "genhypers.h"


int
Process(bool bTargDeps,				// true if process only if any independent src files newer than target
		int PMode,					// processing mode 0: default, 1: summary stats only
		int NumBins,				// when generating length distributions then use this many bins - 0 defaults to using 1000
		int BinDelta,				// when generating length distributions then each bin holds this length delta - 0 defaults to auto determine from NunBins and longest sequence length
		char *pszInputFile,		// bio multialignment (.algn) file to process
		char *pszOutputFile,	// where to write out stats
		char *pszOutputCoreFile, // where to write out the hypercore loci 
		char *pszBiobedFile,	// biobed file containing regional features - exons, introns etc
		int NumIncludeFiles,	// number of include region files
		char **ppszIncludeFiles,// biobed files containing regions to include - default is to exclude none
		int NumExcludeFiles,	// number of exclude region files
		char **ppszExcludeFiles,// biobed file containing regions to exclude - default is to include all
		char *pszUniqueSpeciesSeqs, // ignore alignment block if these species sequences are not unique
		int WindowSize,			// sampling window size
		int NumCoreSpecies,		// number of core species to be in alignment
		char* pszCoreSpecies,	// space or comma separated list of species which must be present in any MAF block alignment before that block will be further processed
		int MinNonCoreSpecies,	// minimum number of species required in an alignment (excluding core species)
		int MinHyperCoreLen,	// minimum hyper core length required
		int MaxHyperColsMismatches,	// hyper cores can have up to this number of columns with at least one mismatches
		int VMismatches,		// number of mismatches in an alignment col to count as a mismatch against MaxHyperColsMismatches
		int MinIdentity,		// minimum identity required when processing hyperconserved
		int NumDistSegs,		// number of match distribution profile segments
		int RegLen,				// regulatory region length - up/dn stream of 5/3' 
		bool bMultipleFeatBits,	// if true then accept alignments in which multiple feature bits are set
		bool bInDelsAsMismatches, // treat InDels as if mismatches
		bool bSloughRefInDels,	// slough columns in which the reference base is an InDel
		bool bFiltLoConfidence,	// filter out low confidence subsequences
		bool bFilt1stLast,		// treat 1st and last subsequences as being low confidence
		int MinSubSeqLen,		// subsequences of less than this length are treated as being low confidence
		int MinIdent,			// treat subsequences of less than this identity as being low confidence
		int NumIncludeChroms,	// number of chromosomes explicitly defined to be included
		char **ppszIncludeChroms,	// ptr to array of reg expressions defining chroms to include - overides exclude
		int NumExcludeChroms,	// number of chromosomes explicitly defined to be excluded
		char **ppszExcludeChroms);	// ptr to array of reg expressions defining chroms to exclude



#ifdef _WIN32
int genhypers(int argc, char* argv[])
{
	// determine my process name
	_splitpath(argv[0], nullptr, nullptr, gszProcName, nullptr);
#else
int
genhypers(int argc, char** argv)
{
	// determine my process name
	CUtility::splitpath((char*)argv[0], nullptr, gszProcName);
#endif
int iScreenLogLevel;		// level of screen diagnostics
int iFileLogLevel;			// level of file diagnostics
char szLogFile[_MAX_PATH];	// write diagnostics to this file

int Rslt;

int NumberOfProcessors;		// number of installed CPUs


bool bTargDeps;						// true if only process if independent files newer than target

int PMode;
int NumBins;				// when generating length distributions then use this many bins - 0 defaults to using 1000
int BinDelta;				// when generating length distributions then each bin holds this length delta - defaults to 1

bool bMultipleFeatBits;
bool bInDelsAsMismatches;
bool bSloughRefInDels;
bool bFiltLoConfidence;

int LenReq;
int Idx;
char szOutputFile[_MAX_PATH];
char szInputFile[_MAX_PATH];
char szOutCoreFile[_MAX_PATH];
char szInputBiobedFile[_MAX_PATH];
int NumIncludeChroms;
char *pszIncludeChroms[cMaxIncludeChroms];
int NumExcludeChroms;
char *pszExcludeChroms[cMaxExcludeChroms];
int NumIncludeFiles;
char *pszIncludeFiles[cMaxExcludeFiles];
int NumExcludeFiles;
char *pszExcludeFiles[cMaxIncludeFiles];
int WindowSize;
char szSpeciesList[cMaxAlignedSpecies * cMaxDatasetSpeciesChrom];
int NumCoreSpecies;
int MinNonCoreSpecies;
int NumSpeciesList;
int RegLen;
int	MinHyperCoreLen;			// minimum hyper core length required
int MaxHyperColsMismatches;		// hyper cores can have at most this number of columns having mismatches
int	MinIdentity;				// minimum identity (50-100) required for hypercores
int DistSegs;					// match distribution profile segments
char szUniqueSpeciesSeqs[512];	// slough blocks in which these species sequences are not unique

char szSQLiteDatabase[_MAX_PATH];	// results summaries to this SQLite file
char szExperimentName[cMaxDatasetSpeciesChrom + 1];			// experiment name
char szExperimentDescr[1000];		// describes experiment

// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_int *ScreenLogLevel=arg_int0("S", "ScreenLogLevel",	"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",	"diagnostics log file");

struct arg_int* pmode = arg_int0("m", "pmode", "<int>", "processing mode 0:default, 1:summary, 2:outspecies");
struct arg_int  *reglen = arg_int0("L","updnstream","<int>",	"length of 5'up or 3'down  stream regulatory region length (default = 2000) 0..1000000");
struct arg_file *infile = arg_file1("i",nullptr,"<file>",			"input from .algn file");
struct arg_file *outfile = arg_file1("o",nullptr,"<file>",			"output to statistics file as CSV");
struct arg_file *outcorefile = arg_file0("O",nullptr,"<file>",		"output hypercore loci to file as CSV");
struct arg_file *inbedfile = arg_file0("r","regions","<file>",		"characterise regions from biobed file");
struct arg_int  *numcorespecies = arg_int0("n","numcorespecies","<int>", "number of core species in alignment, default is num of species in -s<corespecieslist>");
struct arg_int  *minnoncorespecies = arg_int0("N","minnoncorespecies","<int>", "min number of species required in alignment (excludes core species)");
struct arg_str  *specieslist = arg_str0("s","species","<string>","species list, ordered by processing priority, if not specified then defaults to all multialignment species");
struct arg_lit  *multiplefeatbits  = arg_lit0("K","multiplefeatbits",	"single featbit (default) or multiple featbits allowed");
struct arg_lit  *indelsasmismatches  = arg_lit0("j","indelsasmismatches",	"treat InDels same as mismatches (default is to terminate element if InDel)");
struct arg_lit  *sloughrefindels  = arg_lit0("J","sloughrefindels",	"slough columns in which the ref base is an InDel (default is to terminate element if InDel)");
struct arg_lit  *filtloconfidence = arg_lit0("l","filt",		"filter out low confidence subsequences,( slough subsequence < 15mer and first/last subsequence, subsequence must start/end on identical, identity > 70)");
struct arg_int  *minhypercorelen = arg_int0("k","minhypercorelen","<int>",		"minimum (default = 25) required hypercore length (5..1000) - will be maximally extended");
struct arg_int  * maxhypercolsmismatches = arg_int0("X","maxhypercolsmismatches","<int>",	"total number (default = 0) of columns with mismatches allowed in any hypercores (0..500)");
struct arg_int  *minidentity = arg_int0("y","minidentity","<int>",	"minimum percentage identity (default = 100) required in hypercore (50-100)");
struct arg_lit  *chromper = arg_lit0("c","chromper",			"generate stats for each chromosome -= default is for complete genome");
struct arg_int  *numbins = arg_int0("b","numbins","<int>",	"when generating length distributions then use this many bins (defaults to 1000, range 10..10000)");
struct arg_int  *bindelta = arg_int0("B","bindelta","<int>","when generating length distributions then each bin holds this length delta (default 1, range 1,2,5,10,25,50,100,250,500 or 1000)");

struct arg_str  *uniquespeciesseqs = arg_str0("u","uniquespecieseqs","<string>","Ignore alignment blocks in which these species are not unique (use 'any' for all, or '*' for -s<specieslist>)");

struct arg_int  *distsegs = arg_int0("d","distsegs","<int>",	"number of match distribution segments (default 10) used when processing outspecies mode");

struct arg_file *excludefile = arg_filen("E","exclude","<file>",0,cMaxExcludeFiles,	"exclude all regions in biobed file from processing ");
struct arg_file *includefile = arg_filen("I","include","<file>",0,cMaxExcludeFiles,	"include all regions (unless specific regions excluded) in biobed file");
struct arg_int  *windowsize = arg_int0("w","windowsize","<int>","if non-zero then sets fixed size window (in Knt) used to sample along genome");
struct arg_str  *includechroms = arg_strn("z","chromeinclude","<string>",0,cMaxIncludeChroms,"regular expressions defining species.chromosomes to include (overrides exclude) from processing");
struct arg_str  *excludechroms = arg_strn("Z","chromexclude","<string>",0,cMaxExcludeChroms,"regular expressions defining species.chromosomes to exclude from processing");
struct arg_lit  *targdeps = arg_lit0("D","TargDep",				"Generate target file only if missing or older than any of the independent source files");

struct arg_file* summrslts = arg_file0("q", "sumrslts", "<file>", "Output results summary to this SQLite3 database file");
struct arg_str* experimentname = arg_str0("w", "experimentname", "<str>", "experiment name SQLite3 database file");
struct arg_str* experimentdescr = arg_str0("W", "experimentdescr", "<str>", "experiment description SQLite3 database file");

struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,ScreenLogLevel,LogFile,
					
					infile,outfile,outcorefile,inbedfile,
					pmode,numbins,bindelta,uniquespeciesseqs,
					numcorespecies,minnoncorespecies,specieslist,multiplefeatbits,
					indelsasmismatches,sloughrefindels,filtloconfidence,
					minhypercorelen,maxhypercolsmismatches,minidentity,reglen,
					chromper,distsegs,excludefile,includefile,windowsize,includechroms,excludechroms,targdeps,
					summrslts,experimentname,experimentdescr,
					end};

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
		return(1);
		}

	iScreenLogLevel = iFileLogLevel = FileLogLevel->count ? FileLogLevel->ival[0] : eDLInfo;
	if (iFileLogLevel < eDLNone || iFileLogLevel > eDLDebug)
	{
		printf("\nError: FileLogLevel '-l%d' specified outside of range %d..%d\n", iFileLogLevel, eDLNone, eDLDebug);
		return(1);
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
		return(1);
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
	else
		szExperimentName[0] = '\0';

	gExperimentID = 0;
	gProcessID = 0;
	gProcessingID = 0;
	szSQLiteDatabase[0] = '\0';
	szExperimentDescr[0] = '\0';
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

		if (strlen(szExperimentName) < 1)
		{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: After removal of whitespace, no SQLite experiment name specified with '-w<str>' option");
			return(1);
		}
		if (experimentdescr->count)
		{
			strncpy(szExperimentDescr, experimentdescr->sval[0], sizeof(szExperimentDescr) - 1);
			szExperimentDescr[sizeof(szExperimentDescr) - 1] = '\0';
			CUtility::TrimQuotedWhitespcExtd(szExperimentDescr);
		}
		if (strlen(szExperimentDescr) < 1)
		{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: After removal of whitespace, no SQLite experiment description specified with '-W<str>' option");
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
	{
		szSQLiteDatabase[0] = '\0';
		szExperimentDescr[0] = '\0';
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

	NumBins = 1000;
	BinDelta = 1;

	PMode = pmode->count ? pmode->ival[0] : eProcModeStandard;
	if(PMode < eProcModeStandard || PMode > eProcModeOutspecies)
		{
		printf("\nError: Requested processing mode '-x%d' not supported", PMode);
		exit(1);
		}

	MaxHyperColsMismatches = 0;
	MinIdentity = 100;
	if(PMode == eProcModeSummary)
		{	
		szOutCoreFile[0] = '\0';		// where to write out the hypercore loci 
		MinHyperCoreLen =0;				// minimum hyper core length required
		MaxHyperColsMismatches=0;			// hyper cores can have up to this number of columns with mismatches
		MinIdentity=0;					// minimum identity required when processing hyperconserved
		bInDelsAsMismatches=false;		// treat InDels as if mismatches
		bSloughRefInDels =false;		// slough columns in which the reference base is an InDel
		bFiltLoConfidence=false;		// filter out low confidence subsequences
		}
	else	// standard or outspecies processing mode
		{
		MinHyperCoreLen = minhypercorelen->count ? minhypercorelen->ival[0] : cDfltMinHyperCoreLen;
		if(MinHyperCoreLen < cMinHyperCoreLen)
			MinHyperCoreLen = cMinHyperCoreLen;
		else
			if(MinHyperCoreLen > cMaxHyperCoreLen)
				MinHyperCoreLen = cMaxHyperCoreLen;
			
		MaxHyperColsMismatches = maxhypercolsmismatches->count ? maxhypercolsmismatches->ival[0] : cDfltMaxMismatches;
		if(MaxHyperColsMismatches < 0)
			{
			printf("\nMaximum allowed Hypercore columns with mismatches '-X%d' less than 0, assuming you meant to use '-y%d'",MaxHyperColsMismatches,cDfltMaxMismatches);
			MaxHyperColsMismatches = cDfltMaxMismatches;
			}
		else
			if(MaxHyperColsMismatches > cMaxMismatches)
				{
				printf("\nRequested Hypercore columns with mismatches '-X%d' more than than allowed '-X%d', assuming you meant to use '-X%d'",MaxHyperColsMismatches,cMaxMismatches,cMaxMismatches);
				MaxHyperColsMismatches = cDfltMaxMismatches;
				}

		MinIdentity = minidentity->count ? minidentity->ival[0] : cDfltIdentity;
		if(MinIdentity < cMinIdentity)
			{
			printf("\nMinimum identity specified '-y%d' less than %d%%, assuming you meant to use '-y%d'", MinIdentity,cMinIdentity,cMinIdentity);
			MinIdentity = cMinIdentity;
			}
		else
			if(MinIdentity > 100)
				{
				printf("\nMinimum identity specified '-y%d' more than %d, assuming you meant to use '-y%d'", MinIdentity,cMaxIdentity,cMaxIdentity);
				MinIdentity = cMaxIdentity;
				}
		}
	
	if(outcorefile->count)
		strcpy(szOutCoreFile, outcorefile->filename[0]);
	else
		szOutCoreFile[0] = '\0';

	bFiltLoConfidence = filtloconfidence->count ? true : false;
	bInDelsAsMismatches = indelsasmismatches->count ? true : false;
	bSloughRefInDels = sloughrefindels->count ? true : false;
		
	
	NumBins = numbins->count ? numbins->ival[0] : 1000;
	if(NumBins == 0)
		NumBins = 1000;
	if(NumBins < 10 || NumBins > 10000)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Number of bins '-b%d' is not in range 10..10000");
		exit(1);
		}
	BinDelta = bindelta->count ? bindelta->ival[0] : 1;
	if(NumBins == 0)
		NumBins = 1;
	switch(BinDelta) {
		case 1: case 2: case 5: case 10: case 25: case 50: case 100: case 250: case 500: case 1000:
			break;
		default:
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Bin length delta '-B%d' must be either 1,2,5,10,25,50,100,250,500 or 1000");
			exit(1);
		}

	bMultipleFeatBits = multiplefeatbits->count ? true : false;

	strcpy(szInputFile,infile->filename[0]);
	strcpy(szOutputFile,outfile->filename[0]);

	if(inbedfile->count)
		strcpy(szInputBiobedFile,inbedfile->filename[0]);
	else
		szInputBiobedFile[0] = '\0';

	NumIncludeFiles = includefile->count;
	for(Idx=0;Idx < includefile->count; Idx++)
		{
		LenReq = (int)strlen(includefile->filename[Idx]);
		pszIncludeFiles[Idx] = new char [LenReq+1];
		strcpy(pszIncludeFiles[Idx],includefile->filename[Idx]);
		CUtility::TrimQuotes(pszIncludeFiles[Idx]);
		}
	NumExcludeFiles = excludefile->count;
	for(Idx=0;Idx < excludefile->count; Idx++)
		{
		LenReq = (int)strlen(excludefile->filename[Idx]);
		pszExcludeFiles[Idx] = new char [LenReq+1];
		strcpy(pszExcludeFiles[Idx], excludefile->filename[Idx]);
		CUtility::TrimQuotes(pszExcludeFiles[Idx]);
		}

	NumIncludeChroms = includechroms->count;
	for(Idx=0;Idx < includechroms->count; Idx++)
		{
		LenReq = (int)strlen(includechroms->sval[Idx]);
		pszIncludeChroms[Idx] = new char [LenReq+1];
		strcpy(pszIncludeChroms[Idx], includechroms->sval[Idx]);
		CUtility::TrimQuotes(pszIncludeChroms[Idx]);
		}

	NumExcludeChroms = excludechroms->count;
	for(Idx=0;Idx < excludechroms->count; Idx++)
		{
		LenReq = (int)strlen(excludechroms->sval[Idx]);
		pszExcludeChroms[Idx] = new char [LenReq+1];
		strcpy(pszExcludeChroms[Idx], excludechroms->sval[Idx]);
		CUtility::TrimQuotes(pszExcludeChroms[Idx]);
		}

	WindowSize = windowsize->count ? windowsize->ival[0] : 0;
	if(WindowSize < 0)
		{
		printf("\nAllowed sampling result window size '-w%d' less than 0, assuming you meant to use '-w0'", WindowSize);
		WindowSize = 0;
		}

	if(chromper->count && WindowSize == 0) // a hack to force per chromosome stats generation
		WindowSize = 0x1fffffff;

	RegLen = reglen->count ? reglen->ival[0] : cDfltRegLen;
	if(RegLen < cMinRegLen)
		{
		printf("\nRegulatory region length '-L%d' less than minimum %d, assuming you meant to use '-L%d'",RegLen,cMinRegLen,cMinRegLen);
		RegLen = cMinRegLen;
		}
	else
		{
		if(RegLen > cMaxRegLen)
			{
			printf("\nRegulatory region length '-L%d' more than maximum %d, assuming you meant to use '-L%d'",RegLen,cMaxRegLen,cMaxRegLen);
			RegLen = cMaxRegLen;
			}
		}
	
	if(uniquespeciesseqs->count)
		{
		strcpy(szUniqueSpeciesSeqs,uniquespeciesseqs->sval[0]);
		CUtility::TrimQuotes(szUniqueSpeciesSeqs);
		}
	else
		szUniqueSpeciesSeqs[0] = '\0';

	if(!specieslist->count)		// if none specified then default species list to be all present in multialignment file
		{
		CGenHypers *pHypers = new CGenHypers;
		Rslt = pHypers->LoadSpeciesnameList(szInputFile, cMaxAlignedSpecies, sizeof(szSpeciesList), szSpeciesList);
		delete pHypers;
		if (Rslt < 1)
			{
			printf("\nError: Unable to load species list directly from multiple alignment file '%s'\n", szInputFile);
			exit(1);
			}
		}
	else
		strcpy(szSpeciesList,specieslist->sval[0]);
	CUtility::TrimQuotes(szSpeciesList);

	NumSpeciesList = CGenHypers::ParseNumSpecies(szSpeciesList,nullptr);
	if(PMode == eProcModeOutspecies)
		{
		if(NumSpeciesList < 3)
			{
			printf("\nError: At least three (2 core + outspecies) species must be specified in Outspecies mode\n");
			exit(1);
			}
		DistSegs = distsegs->count ? distsegs->ival[0] : cDfltMatchDistSegments;
		if(DistSegs < cMinMatchDistSegments)
			{
			printf("\nWarning: too few distribution segments specified with '-d%d', assume you meant %d segments", DistSegs,cMinMatchDistSegments);
			DistSegs = MinHyperCoreLen;
			}
		else
			if(DistSegs > min(MinHyperCoreLen,cMaxMatchDistSegments))
				{
				printf("\nWarning: too many distribution segments specified with '-d%d', assume you meant %d segments", DistSegs,min(MinHyperCoreLen,cMaxMatchDistSegments));
				DistSegs = min(MinHyperCoreLen,cMaxMatchDistSegments);
				}
		}
	else
		{
		if(NumSpeciesList < 2)		// if multialignment then must be at least 2!
			{
			printf("\nError: At least two species must be specified\n");
			exit(1);
			}
		DistSegs = 0;
		}

	if(PMode != eProcModeOutspecies)
		{
		NumCoreSpecies = numcorespecies->count ? numcorespecies->ival[0] : NumSpeciesList;
		if(NumCoreSpecies < 2)
			NumCoreSpecies = 2;
		if(NumCoreSpecies > NumSpeciesList)
			NumCoreSpecies = NumSpeciesList;

		MinNonCoreSpecies = minnoncorespecies->count ? minnoncorespecies->ival[0] : 0;
		if(MinNonCoreSpecies < 0)
			MinNonCoreSpecies = 0;
		}
	else
		{
		NumCoreSpecies = NumSpeciesList - 1;
		MinNonCoreSpecies = 1;
		}

	if (targdeps->count > 0)
		bTargDeps = true;
	else
		bTargDeps = false;

	if(bTargDeps)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Generate target file only if missing or older than any of the independent source files");
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"bio multialignment (.algn) file to process: '%s'",szInputFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"where to write out stats: '%s'",szOutputFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"biobed file containing regional features: '%s'",szInputBiobedFile);
	switch(PMode) {
		case eProcModeStandard:
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode: standard");
			break;

		case eProcModeSummary:
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode: summary");
			break;

		case eProcModeOutspecies:
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode: outspecies");
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Match distribution profile segments: %d",DistSegs);
			break;
		}

	if(szUniqueSpeciesSeqs[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Ignore alignment blocks where sequences not unique for: '%s'",szUniqueSpeciesSeqs);

	if(PMode != eProcModeSummary)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"where to write out the hypercore loci: '%s'",szOutCoreFile);
	for(Idx = 0; Idx < NumIncludeFiles; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"biobed file containing regions to include: '%s'",pszIncludeFiles[Idx]); 
	for(Idx = 0; Idx < NumExcludeFiles; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"biobed file containing regions to exclude: '%s'",pszExcludeFiles[Idx]);
	if(WindowSize)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"sampling window size: %d",WindowSize);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Number of bins: %d",NumBins);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Bin length delta: %d",BinDelta);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"number of core species to be in alignment: %d",NumCoreSpecies);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"minimum number of non-core species to be in alignment: %d", MinNonCoreSpecies);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"species list: '%s'",	szSpeciesList);
	if(PMode != eProcModeSummary)
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"minimum hyper core length required: %d",MinHyperCoreLen);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"hyper cores can have upto this number of total mismatches: %d",MaxHyperColsMismatches);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"minimum identity required when processing hyperconserved: %d",MinIdentity);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"treat InDels as if mismatches: %s",bInDelsAsMismatches ? "yes" : "no");
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"slough columns in which the reference base is an InDel: %s",bSloughRefInDels ? "yes" : "no");
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"filter out low confidence subsequences: %s",bFiltLoConfidence ? "yes" : "no");
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"treat 1st and last subsequences as being low confidence: %s",bFiltLoConfidence ? "yes" : "no");
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"subsequences of less than this length are treated as being low confidence: %d",bFiltLoConfidence ? 15 : 0);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"treat subsequences of less than this identity as being low confidence: %d",bFiltLoConfidence? 70 : 0);
		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"regulatory region length: %d",RegLen);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"accept alignments in which multiple feature bits are set: %s",bMultipleFeatBits ? "yes" : "no");
	for(Idx = 0; Idx < NumIncludeChroms; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"reg expressions defining chroms to include: '%s'",pszIncludeChroms[Idx]);
	for(Idx = 0; Idx < NumExcludeChroms; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"reg expressions defining chroms to exclude: '%s'",pszExcludeChroms[Idx]); 

	gStopWatch.Start();
#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	Rslt = Process(bTargDeps,				// true if process only if any independent src files newer than target
					PMode,					// processing mode 0: default, 1: summary stats only, 2: outspecies processing
					NumBins,				// when generating length distributions then use this many bins - 0 defaults to using 1000
					BinDelta,				// when generating length distributions then each bin holds this length delta - 0 defaults to auto determine from NunBins and longest sequence length
					szInputFile,		// bio multialignment (.algn) file to process
					szOutputFile,		// where to write out stats
					szOutCoreFile,		// where to write out the hypercore loci 
					szInputBiobedFile,	// biobed file containing regional features - exons, introns etc
					NumIncludeFiles,	// number of include region files
					pszIncludeFiles,	// biobed files containing regions to include - default is to exclude none
					NumExcludeFiles,	// number of exclude region files
					pszExcludeFiles,	// biobed file containing regions to exclude - default is to include all
					szUniqueSpeciesSeqs, // ignore alignment blocks if these species sequences are not unique
					WindowSize,		// sampling window size
					NumCoreSpecies,	// number of core species to be in alignment
					szSpeciesList,		// space or comma separated list of core species
					MinNonCoreSpecies,	// minimum number of non-core species to be in alignment
					MinHyperCoreLen,		// minimum hyper core length required
					MaxHyperColsMismatches,	// hyper cores can have up to this number of columns with mismatches
					1,					// max number of mismatches in an alignment col to count as a mismatch against MaxHyperColsMismatches
					MinIdentity,		// minimum identity required when processing hyperconserved
					DistSegs,			// number of match distribution profile segments
					RegLen,			// regulatory region length - up/dn stream of 5/3' 
					bMultipleFeatBits,	// if true then accept alignments in which multiple feature bits are set
					bInDelsAsMismatches, // treat InDels as if mismatches
					bSloughRefInDels,		// slough columns in which the reference base is an InDel
					bFiltLoConfidence,		// filter out low confidence subsequences
					bFiltLoConfidence,		// treat 1st and last subsequences as being low confidence
					bFiltLoConfidence ? 15 : 0, // subsequences of less than this length are treated as being low confidence
					bFiltLoConfidence? 70 : 0, // treat subsequences of less than this identity as being low confidence
					NumIncludeChroms,	// number of chromosomes explicitly defined to be included
					pszIncludeChroms,	// ptr to array of reg expressions defining chroms to include - overides exclude
					NumExcludeChroms,	// number of chromosomes explicitly defined to be excluded
					pszExcludeChroms);	// ptr to array of reg expressions defining chroms to be excluded

	gStopWatch.Stop();
	Rslt = Rslt < 0 ? 1 : 0;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit Code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	exit(Rslt);
	}
else
	{
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nend of help\n");
	exit(1);
	}
}

int
Process(bool bTargDeps,				// true if process only if any independent src files newer than target
	int ProcMode,					// processing mode 0: default, 1: summary stats only
	int NumBins,				// when generating length distributions then use this many bins - 0 defaults to using 1000
	int BinDelta,				// when generating length distributions then each bin holds this length delta - 0 defaults to auto determine from NunBins and longest sequence length
	char* pszInputFile,		// bio multialignment (.algn) file to process
	char* pszOutputFile,	// where to write out stats
	char* pszOutputCoreFile, // where to write out the hypercore loci 
	char* pszBiobedFile,	// biobed file containing regional features - exons, introns etc
	int NumIncludeFiles,	// number of include region files
	char** ppszIncludeFiles,// biobed files containing regions to include - default is to exclude none
	int NumExcludeFiles,	// number of exclude region files
	char** ppszExcludeFiles,// biobed file containing regions to exclude - default is to include all
	char* pszUniqueSpeciesSeqs, // ignore alignment block if these species sequences are not unique
	int WindowSize,			// sampling window size
	int NumCoreSpecies,		// number of core species to be in alignment
	char* pszCoreSpecies,	// space or comma separated list of species which must be present in any MAF block alignment before that block will be further processed
	int MinNonCoreSpecies,	// minimum number of species required in an alignment (excluding core species)
	int MinHyperCoreLen,		// minimum hyper core length required
	int MaxHyperColsMismatches,	// hyper cores can have upto this number of columns with mismatches
	int VMismatches,		// number of mismatches in an alignment col to count as a mismatch against MaxHyperColsMismatches
	int MinIdentity,		// minimum identity required when processing hyperconserved
	int NumDistSegs,		// number of match distribution profile segments
	int RegLen,				// regulatory region length - up/dn stream of 5/3' 
	bool bMultipleFeatBits,	// if true then accept alignments in which multiple feature bits are set
	bool bInDelsAsMismatches, // treat InDels as if mismatches
	bool bSloughRefInDels,	// slough columns in which the reference base is an InDel
	bool bFiltLoConfidence,	// filter out low confidence subsequences
	bool bFilt1stLast,		// treat 1st and last subsequences as being low confidence
	int MinSubSeqLen,		// subsequences of less than this length are treated as being low confidence
	int MinIdent,			// treat subsequences of less than this identity as being low confidence
	int NumIncludeChroms,	// number of chromosomes explicitly defined to be included
	char** ppszIncludeChroms,	// ptr to array of reg expressions defining chroms to include - overides exclude
	int NumExcludeChroms,	// number of chromosomes explicitly defined to be excluded
	char** ppszExcludeChroms)	// ptr to array of reg expressions defining chroms to exclude
{
int Rslt;
CGenHypers *pGenHypers;
if((pGenHypers = new CGenHypers) == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to instantiate instance of CRNA_DE");
	return(eBSFerrObj);
	}

Rslt = pGenHypers->Process(bTargDeps,			// true if process only if any independent src files newer than target
						ProcMode,				// processing mode 0: default, 1: summary stats only
						NumBins,				// when generating length distributions then use this many bins - 0 defaults to using 1000
						BinDelta,				// when generating length distributions then each bin holds this length delta - 0 defaults to auto determine from NunBins and longest sequence length
						pszInputFile,			// bio multialignment (.algn) file to process
						pszOutputFile,			// where to write out stats
						pszOutputCoreFile,		// where to write out the hypercore loci 
						pszBiobedFile,			// biobed file containing regional features - exons, introns etc
						NumIncludeFiles,		// number of include region files
						ppszIncludeFiles,		// biobed files containing regions to include - default is to exclude none
						NumExcludeFiles,		// number of exclude region files
						ppszExcludeFiles,		// biobed file containing regions to exclude - default is to include all
						pszUniqueSpeciesSeqs,	// ignore alignment block if these species sequences are not unique
						WindowSize,				// sampling window size
						NumCoreSpecies,			// number of core species to be in alignment
						pszCoreSpecies,			// space or comma separated list of species which must be present in any MAF block alignment before that block will be further processed
						MinNonCoreSpecies,		// minimum number of species required in an alignment (excluding core species)
						MinHyperCoreLen,			// minimum hyper core length required
						MaxHyperColsMismatches,	// hyper cores can have up to this number of columns with mismatches
						VMismatches,			// number of mismatches in an alignment col to count as a mismatch against MaxHyperColsMismatches
						MinIdentity,			// minimum identity required when processing hyperconserved
						NumDistSegs,			// number of match distribution profile segments
						RegLen,					// regulatory region length - up/dn stream of 5/3' 
						bMultipleFeatBits,		// if true then accept alignments in which multiple feature bits are set
						bInDelsAsMismatches,	// treat InDels as if mismatches
						bSloughRefInDels,		// slough columns in which the reference base is an InDel
						bFiltLoConfidence,		// filter out low confidence subsequences
						 bFilt1stLast,			// treat 1st and last subsequences as being low confidence
						MinSubSeqLen,			// subsequences of less than this length are treated as being low confidence
						MinIdent,				// treat subsequences of less than this identity as being low confidence
						NumIncludeChroms,		// number of chromosomes explicitly defined to be included
						ppszIncludeChroms,		// ptr to array of reg expressions defining chroms to include - overides exclude
						NumExcludeChroms,		// number of chromosomes explicitly defined to be excluded
						ppszExcludeChroms);		// ptr to array of reg expressions defining chroms to exclude

if(pGenHypers != nullptr)
	delete pGenHypers;
return(Rslt);
}

CGenHypers::CGenHypers()
{
m_pLenRangeClasses = nullptr;
m_pszOutBuffer = nullptr;
Reset();
}

CGenHypers::~CGenHypers()
{
if (m_pszOutBuffer != nullptr)
	delete[]m_pszOutBuffer;
if (m_pLenRangeClasses != nullptr)
	delete[]m_pLenRangeClasses;
}

void
CGenHypers::Reset(void)
{

if(m_pszOutBuffer != nullptr)
	{
	delete []m_pszOutBuffer;
	m_pszOutBuffer = nullptr;
	}
m_AllocOutBuff = 0;
m_OutBuffIdx = 0;

m_NumLenRangeBins = 1000;				// when generating length distributions then use this many bins - 0 defaults to using 1000
m_BinDelta = 1;				// when generating length distributions then each bin holds this length delta - defaults to 1
if(m_pLenRangeClasses != nullptr)
	{
	delete []m_pLenRangeClasses;
	m_pLenRangeClasses = nullptr; // allocated to hold length ranges
	}

m_NumExcludeEls = 0;		// current number of elements in gExcludeChroms
m_pMRA = nullptr;		// pts to most recently accessed or added
m_pLRA = nullptr;		// pts to least recently accessed
memset(m_ExcludeChroms,0,sizeof(m_ExcludeChroms));
}

int
CGenHypers::InitLengthRangeClass(		int NumBins,				// when generating length distributions then use this many bins - 0 defaults to using 1000
		int BinDelta)				// when generating length distributions then each bin holds this length delta - 0 defaults to auto determine from NunBins and longest sequence length
{
tsLenRangeClass *pCurRange;
int Idx;
if(m_pLenRangeClasses != nullptr)
	delete m_pLenRangeClasses;
m_pLenRangeClasses = new tsLenRangeClass [NumBins+1];
pCurRange = m_pLenRangeClasses;
for(Idx = 0; Idx < NumBins; Idx++,pCurRange++)
	{
	pCurRange->ID = Idx+1;
	pCurRange->Min = Idx * BinDelta;
	pCurRange->Max = (Idx + 1) * BinDelta;
	if(BinDelta == 1)
		{
		if(Idx == NumBins-1)
			sprintf(pCurRange->szDescr,"%d plus",pCurRange->Min);
		else
			sprintf(pCurRange->szDescr,"%d",pCurRange->Min);
		}
	else
		{
		if(Idx == NumBins-1)
			sprintf(pCurRange->szDescr,"%d plus",pCurRange->Min);
		else
			sprintf(pCurRange->szDescr,"%d to %d",pCurRange->Min,pCurRange->Max-1);
		}
	}

return(NumBins);
}

// GetLengthRangeClass
// Returns ptr to length range class for specified Length, or nullptr if can't classify into a range
tsLenRangeClass *
CGenHypers::GetLengthRangeClass(int Length)
{
int Idx;
tsLenRangeClass *pRange = m_pLenRangeClasses;

for(Idx = 0; Idx < m_NumLenRangeBins; Idx++,pRange++)
	if(Length >= pRange->Min && Length < pRange->Max)
		return(pRange);
return(&m_pLenRangeClasses[m_NumLenRangeBins -1]);
}

// GetRangeClass
// Returns ptr to length range class for specified range identifier, or nullptr if can't classify into a range
tsLenRangeClass *
CGenHypers::GetRangeClass(int RangeID)
{
if(RangeID < 1 || RangeID > m_NumLenRangeBins)
	return(nullptr);
return(&m_pLenRangeClasses[RangeID-1]);
}

// ParseNumSpecies
// Initialises pProcParams with parsed species names in space or comma delimited list ptd at by pszSpeciesList
// Returns number of species
int
CGenHypers::ParseNumSpecies(char *pszSpeciesList,tsProcParams *pProcParams)
{
// parse out species list
char Chr;
char *pSpecies;
int NumSpecies = 0;
bool InToken = false;
if(pszSpeciesList == nullptr || *pszSpeciesList == '\0')
	return(0);

while(Chr = *pszSpeciesList++)
	{
	if(Chr == '"' || Chr == '\'') // change any single or double quotes into spaces
		Chr = ' ';
	if(isspace(Chr) || Chr==',')
		{
		if(!InToken)			// slough whitespace or ',' if not inside a token parse
			continue;
		InToken = false;
		pszSpeciesList[-1] = '\0';
		if(pProcParams != nullptr)
			{
			strncpy(pProcParams->szSpecies[NumSpecies],pSpecies,cMaxDatasetSpeciesChrom);
			pProcParams->szSpecies[NumSpecies][cMaxDatasetSpeciesChrom-1] = '\0';
			}
		pszSpeciesList[-1] = Chr;
		NumSpecies++;
		if(NumSpecies >= cMaxAlignedSpecies)
			break;
		continue;
		}
	if(!InToken)			// if not already inside token then start token 
		{
		pSpecies = pszSpeciesList-1;
		InToken = true;
		}
	}
if(InToken)
	{
	pszSpeciesList[-1] = '\0';
	if(pProcParams != nullptr)
		{
		strncpy(pProcParams->szSpecies[NumSpecies],pSpecies,cMaxDatasetSpeciesChrom);
		pProcParams->szSpecies[NumSpecies][cMaxDatasetSpeciesChrom-1] = '\0';
		}
	pszSpeciesList[-1] = Chr;
	NumSpecies++;
	}
return(NumSpecies);
}

// ParseUniqueSpeciesSeqs
// Initialises pProcParams with parsed species names in space or comma delimited list ptd at by pszSpeciesList
// Returns number of species
int
CGenHypers::ParseUniqueSpeciesSeqs(char *pszUniqueSpeciesSeqs,tsProcParams *pProcParams)
{
// parse out species list for which sequences must be unique
char Chr;
char *pSpecies;
int NumSpecies = 0;
bool InToken = false;

pProcParams->bAllUniqueSpeciesSeqs = false;
pProcParams->NumUniqueSpeciesSeqs = 0;
if(pszUniqueSpeciesSeqs == nullptr || *pszUniqueSpeciesSeqs == '\0')
	return(0);

if(!stricmp(pszUniqueSpeciesSeqs,"all"))
	{
	pProcParams->bAllUniqueSpeciesSeqs = true;
	return(cMaxAlignedSpecies);
	}

while(Chr = *pszUniqueSpeciesSeqs++)
	{
	if(Chr == '"' || Chr == '\'') // change any single or double quotes into spaces
		Chr = ' ';
	if(isspace(Chr) || Chr==',')
		{
		if(!InToken)			// slough whitespace or ',' if not inside a token parse
			continue;
		InToken = false;
		pszUniqueSpeciesSeqs[-1] = '\0';
		if(pProcParams != nullptr)
			{
			strncpy(pProcParams->UniqueSpeciesSeqs[NumSpecies],pSpecies,cMaxDatasetSpeciesChrom);
			pProcParams->UniqueSpeciesSeqs[NumSpecies][cMaxDatasetSpeciesChrom-1] = '\0';
			}
		pszUniqueSpeciesSeqs[-1] = Chr;
		NumSpecies++;
		if(NumSpecies >= cMaxAlignedSpecies)
			break;
		continue;
		}
	if(!InToken)			// if not already inside token then start token 
		{
		pSpecies = pszUniqueSpeciesSeqs-1;
		InToken = true;
		}
	}
if(InToken)
	{
	pszUniqueSpeciesSeqs[-1] = '\0';
	if(pProcParams != nullptr)
		{
		strncpy(pProcParams->UniqueSpeciesSeqs[NumSpecies],pSpecies,cMaxDatasetSpeciesChrom);
		pProcParams->UniqueSpeciesSeqs[NumSpecies][cMaxDatasetSpeciesChrom-1] = '\0';
		}
	pszUniqueSpeciesSeqs[-1] = Chr;
	NumSpecies++;
	}
pProcParams->NumUniqueSpeciesSeqs = NumSpecies;
return(NumSpecies);
}


// AddExcludeHistory
// Adds a tsExcludeEl to the cached history
// The newly added element will be the MRA
bool
CGenHypers::AddExcludeHistory(int SpeciesID,int ChromID,bool bExclude)
{
tsExcludeEl *pEl;
if(m_NumExcludeEls < cMaxExcludeHistory)
	pEl = &m_ExcludeChroms[m_NumExcludeEls++];
else
	{
	pEl = m_pLRA;		// reuse the least recently accessed element
	m_pLRA = pEl->pPrev;	// 2nd LRA now becomes the LRA
	m_pLRA->pNext = nullptr;
	}
if(m_pMRA != nullptr)
	m_pMRA->pPrev = pEl;
pEl->pNext = m_pMRA;
pEl->pPrev = nullptr;
m_pMRA = pEl;
if(m_pLRA == nullptr)
	m_pLRA = pEl;
pEl->bExclude = bExclude;
pEl->SpeciesID = SpeciesID;
pEl->ChromID = ChromID;
return(bExclude);
}

// LocateExclude
// Locates - starting from the MRA - a tsExcludeEl which matches on SpeciesID and ChromID
// If matches then this tsExcludeEl is made the MRA
// Returns ptr to matching tsExcludeEl if located or nullptr
tsExcludeEl *
CGenHypers::LocateExclude(int SpeciesID,int ChromID)
{
tsExcludeEl *pEl;
pEl = m_pMRA;
while(pEl != nullptr)
	{
	if(pEl->SpeciesID == SpeciesID && pEl->ChromID == ChromID)
		{
		if(m_NumExcludeEls ==1 || m_pMRA == pEl)	// if only, or already the MRA then no need for any relinking 
			return(pEl);

		if(m_pLRA == pEl)						// if was the LRA then the 2nd LRA becomes the LRA
			{
			m_pLRA = pEl->pPrev;
			m_pLRA->pNext = nullptr;
			}
		else									// not the LRA, and not the MRA
			{
			pEl->pPrev->pNext = pEl->pNext;
			pEl->pNext->pPrev = pEl->pPrev;
			}
		m_pMRA->pPrev = pEl;
		pEl->pNext = m_pMRA;
		pEl->pPrev = nullptr;
		m_pMRA = pEl;
		return(pEl);
		}
	pEl = pEl->pNext;
	}
return(nullptr);
}


// ExcludeThisChrom
// Returns true if SpeciesID.ChromosomeID is to be excluded from processing
// ExcludeThisChrom
bool
CGenHypers::ExcludeThisChrom(CMAlignFile *pAlignments,int SpeciesID,int ChromID,tsProcParams *pProcParams)
{
char *pszSpecies;
char *pszChrom;
char szSpeciesChrom[200];

if(!m_RegExprs.HasRegExprs())
	return(false);
// check if this species and chromosome are already known to be included/excluded
tsExcludeEl *pEl;
if((pEl = LocateExclude(SpeciesID,ChromID))!=nullptr)
	return(pEl->bExclude);
// haven't seen this species or chromosome before - or else they have been discarded from history...
pszSpecies = pAlignments->GetSpeciesName(SpeciesID);
pszChrom = pAlignments->GetChromName(ChromID);
sprintf(szSpeciesChrom,"%s.%s",pszSpecies,pszChrom);

// to be excluded or excluded?
return(AddExcludeHistory(SpeciesID,ChromID,!m_RegExprs.Accept(szSpeciesChrom)));
}


// Directly load a space separated list of all species contained in a multialignment file
int
CGenHypers::LoadSpeciesnameList(char* pszAlgn,			 // source bioseq multialignment file
							int MaxSpecies,			 // can accept at most this many multialigned species
							int SpeciesBuffAllocSize, // pszSpeciesBuff is caller allocated to hold a maximum of this many chars including the string null terminator
							char *pszSpeciesBuff)	// all species are written (space separated) into this caller allocated buffer
{
int NumSpecies;
char *pszSpecies;
int CurBuffOfs;
int SpeciesNameLen;
int Rslt;
CMAlignFile* pAlignments;
int SpeciesID;

if ((pAlignments = new CMAlignFile()) == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to create new instance of CMAlignFile");
	return(eBSFerrObj);
	}

if ((Rslt = pAlignments->Open(pszAlgn)) != eBSFSuccess)
	{
	while (pAlignments->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal, gszProcName, pAlignments->GetErrMsg());
	delete pAlignments;
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadSpeciesnameList: Unable to open Algn file %s\n", pszAlgn);
	return(Rslt);
	}
NumSpecies = pAlignments->GetNumSpecies(0);
if (NumSpecies > MaxSpecies)
	{
	delete pAlignments;
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadSpeciesnameList: Too many species, allowed %d actual %d, in multialignment Algn file %s\n", MaxSpecies, NumSpecies, pszAlgn);
	return(eBSFerrSpecies);
	}

	// iterate over all species which are present in the multialignments, copying into pszSpeciesBuff with space separators
CurBuffOfs = 0;
for (SpeciesID = 1; SpeciesID <= NumSpecies; SpeciesID++)
	{
	pszSpecies = pAlignments->GetSpeciesName(SpeciesID);
	SpeciesNameLen = (int)strlen(pszSpecies);
	if((SpeciesNameLen + CurBuffOfs + 1) > SpeciesBuffAllocSize)
		{
		delete pAlignments;
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadSpeciesnameList: Insufficient buffering allowed for species names in multialignment Algn file %s\n", pszAlgn);
		return(eBSFerrSpecies);
		}
	CurBuffOfs += sprintf(&pszSpeciesBuff[CurBuffOfs], "%s ",pszSpecies);
	}
pszSpeciesBuff[CurBuffOfs] = '\0';		// overwrites space char prefixed onto last species name
delete pAlignments;
return(NumSpecies);
}

int 
CGenHypers::ProcessAlignments(char *pszMAF,			 // source bioseq multialignment file
				  tsProcParams *pProcParams) // processing parameters
{
int RefSpeciesID;
int RefChromID;
int BEDChromID;
int PrevRefChromID;
int PrevDispRefChromID;
char *pszRefChrom;
int RefChromOfs;
char RefStrand;
int SpeciesIDs[cMaxAlignedSpecies];
int Rslt;
int RefAlignLen;
CMAlignFile *pAlignments;
int Idx;
int CurBlockID;
bool bLoaded;

if((pAlignments = new CMAlignFile())==nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create new instance of CMAlignFile");
	return(eBSFerrObj);
	}

if((Rslt=pAlignments->Open(pszMAF))!=eBSFSuccess)
	{
	while(pAlignments->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pAlignments->GetErrMsg());
	delete pAlignments;
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open MAF file %s\n",pszMAF);
	return(Rslt);
	}

// ensure all species are represented in multispecies alignment file plus get their species identifiers
for(Idx = 0; Idx < pProcParams->NumSpeciesList; Idx++)
	{
	if((Rslt = SpeciesIDs[Idx] = pAlignments->LocateSpeciesID(pProcParams->szSpecies[Idx]))<1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Species '%s' not represented in %s",pProcParams->szSpecies[Idx],pszMAF);
		Rslt = pAlignments->GetNumSpecies();
		for(Idx=1;Idx<=Rslt;Idx++)
			gDiagnostics.DiagOut(eDLFatal,gszProcName," Represented species: %s",pAlignments->GetSpeciesName(Idx));
		delete pAlignments;
		return(eBSFerrEntry);
		}
	if(!Idx)	// reference species is always the first species in the species name list
		RefSpeciesID = SpeciesIDs[0];
	}

if(pProcParams->bAllUniqueSpeciesSeqs)
	pAlignments->SetAllConfFilt(true);
else
	{
	pAlignments->SetAllConfFilt(false);
	if(pProcParams->NumUniqueSpeciesSeqs)
		{
		for(Idx = 0; Idx < pProcParams->NumUniqueSpeciesSeqs; Idx++)
			{
			if((Rslt = (int)pAlignments->LocateSpeciesID(pProcParams->UniqueSpeciesSeqs[Idx]))<1)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unique sequence species '%s' not represented in %s",pProcParams->UniqueSpeciesSeqs[Idx],pszMAF);
				Rslt = pAlignments->GetNumSpecies();
				for(Idx=1;Idx<=Rslt;Idx++)
					gDiagnostics.DiagOut(eDLFatal,gszProcName," Represented species: %s",pAlignments->GetSpeciesName(Idx));
				delete pAlignments;
				return(eBSFerrEntry);
				}
			pAlignments->SetConfFilt((tSpeciesID)Rslt,true);
			}
		}
	}

// iterate over reference blocks which are sorted by chrom then offset
CurBlockID = 0;
PrevRefChromID = 0;
PrevDispRefChromID = 0;
BEDChromID = 0;
pProcParams->RefSpeciesIdx = 0;		// reference sequence will always be 1st
pProcParams->NxtOutputOffset = pProcParams->WindowSize;

while(CurBlockID >= 0 && ((CurBlockID =						// returned blockid to next start loading from
	LoadContiguousBlocks(RefSpeciesID,	// reference species identifier
			   CurBlockID,			// which block to initially start loading from
			   &bLoaded,			// returned indicator as to if any loaded blocks meet processing requirements
			   &RefChromID,			// returned reference chromosome identifier 
			   &RefStrand,			// returned reference strand
			   &RefAlignLen,		// returned alignment (incl InDels) length
			   &RefChromOfs,		// returned alignment start offset
			   SpeciesIDs,			// input - species of interest identifier array
			   pAlignments,
			   pProcParams)) > 0 || (CurBlockID == eBSFerrAlignBlk && RefAlignLen > 0)))
	{
	if(RefChromID > 0 && RefChromID != PrevDispRefChromID)
		{
		if(PrevDispRefChromID > 0)
			{
			if(pProcParams->PMode == eProcModeSummary)
				ChkOutputSummaryResults(pProcParams->szRefChrom, -1,pProcParams,false,false);
			else
				ChkOutputResults(pProcParams->szRefChrom, -1,pProcParams,false,false);
			}
		PrevDispRefChromID = RefChromID;
		pszRefChrom = pAlignments->GetChromName(RefChromID);

		strcpy(pProcParams->szRefChrom,pszRefChrom);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing chromosome %s",pszRefChrom);
		pProcParams->NxtOutputOffset = pProcParams->WindowSize;
		if(pProcParams->pBiobed != nullptr)
			pProcParams->BEDChromID = pProcParams->pBiobed->LocateChromIDbyName(pProcParams->szRefChrom);
		else
			pProcParams->BEDChromID = 0;
		pProcParams->RefChromID = RefChromID;
		}

	if(!bLoaded)
		continue;

		// not interested if the alignment length would be too short unless extended stats being generated
	if(pProcParams->PMode != eProcModeSummary && RefAlignLen < pProcParams->MinHyperCoreLen)
		continue;
	if(pProcParams->bFiltLoConfidence && RefAlignLen < pProcParams->MinSubSeqLen)
		continue;

		// here we need to normalise the alignments so that there will be no case of all InDels in
		// any column which can occur if none of the processed species is not the reference species
	if((RefAlignLen = NormaliseInDelColumns(pProcParams,RefAlignLen))< pProcParams->MinHyperCoreLen)
		continue;
	if(pProcParams->bFiltLoConfidence && RefAlignLen < pProcParams->MinSubSeqLen)
		continue;

	if(RefChromID != PrevRefChromID)
		PrevRefChromID = RefChromID;

	if(pProcParams->PMode == eProcModeSummary)
		{
		ChkOutputSummaryResults(pProcParams->szRefChrom, RefChromOfs,pProcParams,false,false);
		if(!ProcAlignBlockSummary(RefChromID,RefChromOfs,RefAlignLen,pProcParams))
			break;
		}
	else
		{
		ChkOutputResults(pProcParams->szRefChrom, RefChromOfs,pProcParams,false,false);
		if(!ProcAlignBlock(RefChromID,RefChromOfs,RefAlignLen,pProcParams))
			break;
		}
	}


delete pAlignments;
return(eBSFSuccess);
}

// LoadContiguousBlocks() loads sequences from 1 or more blocks starting at BlockID which are contiguous in the reference sequence.
// Each individual block loaded may differ in the number of relatively aligned species so if any species is absent then its
// corresponding sequence will be set to be eBaseUndef 
int										// returned blockid to next start loading from
CGenHypers::LoadContiguousBlocks(int RefSpeciesID,	// reference species identifier
			   int  BlockID,			// which block to initially start loading from
			   bool *pbLoaded,			// returned indicator as to if any loaded blocks meet processing requirements
			   int *pRefChromID,		// returned reference chromosome identifier 
			   char *pRefStrand,		// returned reference strand
			   int *pRefAlignLen,		// returned alignment (incl InDels) length
			   int *pRefChromOfs,		// returned alignment start offset
			   int *pSpeciesIDs,		// input - species of interest identifier array
			   CMAlignFile *pAlignments,
			   tsProcParams *pProcParams)
{
int CurBlockID;
int PrvBlockID;
int bFirst;
int CurNumAlignments;
int NumSpecies;
int MaxAlignIdxSpecies;
int Idx;
etSeqBase *pSeq;

int	 CurRefChromID;
char CurRefStrand;
int	 CurRefChromOfs;
int	 CurRefChromEndOfs;
int  CurRefAlignLen;

int	 RefChromID;
char RefStrand;
int	 RefChromOfs;
int	 RefChromEndOfs;
int  RefAlignLen;

int RelSpeciesID;
int RelChromID;

bFirst = true;
CurBlockID = BlockID;
PrvBlockID = CurBlockID;
pProcParams->NumSpeciesInAlignment = 0;
pProcParams->MaxAlignIdxSpecies = 0;
RefAlignLen=   0;
CurRefChromID = 0;
CurRefStrand = '+';
CurRefChromOfs = -1;
CurRefChromEndOfs = -1;
CurRefAlignLen = 0;
MaxAlignIdxSpecies = 0;
while(CurBlockID >= 0 && ((CurBlockID = pAlignments->NxtBlock(CurBlockID)) > 0))
	{
	CurRefChromID  = pAlignments->GetRelChromID(CurBlockID,RefSpeciesID);
	CurRefStrand   = pAlignments->GetStrand(CurBlockID,RefSpeciesID);
	CurRefChromOfs = pAlignments->GetRelChromOfs(CurBlockID,RefSpeciesID);
	CurRefChromEndOfs = pAlignments->GetRelChromEndOfs(CurBlockID,RefSpeciesID);
	CurRefAlignLen = pAlignments->GetAlignLen(CurBlockID,RefSpeciesID);
	
			// terminate with blocks in which there are fewer species than the minimum number we need!
	if((CurNumAlignments = pAlignments->GetNumSpecies(CurBlockID)) < pProcParams->NumCoreSpecies)
		break;

	// check if this aligned block would overflow alloc'd memory for holding alignments when block sequences are
	// concatenated with a previously loaded block in the current call to LoadContiguousBlocks().
	// 
	// a copy of each species aligned sequence will be made into pProcParms->pSeqs[SpeciesID]
	if((CurRefAlignLen + RefAlignLen) > pProcParams->MaxSeqAlignLen)
		{
		CurBlockID = PrvBlockID; // force re-read of block next time LoadContiguousBlocks() is called
		break;
		}

	if(bFirst) // if initial block in current call to LoadContiguousBlocks()
		{
		RefChromID =   CurRefChromID;				// noting which chrom, strand, offsets this initial block covers
		RefStrand  =   CurRefStrand;
		RefChromOfs=   CurRefChromOfs;
		RefChromEndOfs=CurRefChromOfs;
		}
	else
		{	// not the initial block in current call to LoadContiguousBlocks()
		if(CurRefChromID != RefChromID ||			// if concatenating sequences then need to ensure the concatenation would be on same chrom as the initial block sequences 
		   CurRefStrand != RefStrand || 
		   (CurRefChromOfs != (RefChromEndOfs + 1)))	// need to ensure that this blocks sequences start immediately following previous blocks end offset 
			{
			CurBlockID = PrvBlockID; // force re-read of block next time
			break;
			}
		}

		// iterate over species aligned in current block, species are assumed to be in priority order
	MaxAlignIdxSpecies = 0;
	int SpeciesExpected = 0;
	for(NumSpecies = Idx = 0; Idx < pProcParams->NumSpeciesList; Idx++) 
		{
		RelSpeciesID = pSpeciesIDs[Idx];
		RelChromID = pAlignments->GetRelChromID(CurBlockID, RelSpeciesID);
		if(RelChromID < 1)			// assume < 1 is because species not in alignment
			{
			if(RelSpeciesID == 1)
				SpeciesExpected++;
			memset(&pProcParams->pSeqs[Idx][RefAlignLen],eBaseUndef,CurRefAlignLen);
			continue;		
			}

		if(ExcludeThisChrom(pAlignments,RelSpeciesID,RelChromID,pProcParams))	// should this chromosome be accepted for processing or sloughed?
			{
			memset(&pProcParams->pSeqs[Idx][RefAlignLen],eBaseUndef,CurRefAlignLen);
			continue;	
			}

			// get sequence for the species and if present - expected to be! - then concatenate on to any previously loaded sequences
		if((pSeq = pAlignments->GetSeq(CurBlockID,RelSpeciesID))==nullptr) 
			{
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadContiguousBlocks: Unexpected missing sequence in alignments");
			break;
			}
		
		memcpy(&pProcParams->pSeqs[Idx][RefAlignLen],pSeq,CurRefAlignLen);
		NumSpecies++;
		MaxAlignIdxSpecies = Idx + 1;
		}
		
	// check if required minimum number of species were in multialignment
	if(NumSpecies < pProcParams->MinAlignSpecies)	
		break;

	if(NumSpecies > pProcParams->NumSpeciesInAlignment)
		pProcParams->NumSpeciesInAlignment = NumSpecies;
	if(MaxAlignIdxSpecies > pProcParams->MaxAlignIdxSpecies)
		pProcParams->MaxAlignIdxSpecies = MaxAlignIdxSpecies;
	RefAlignLen += CurRefAlignLen;
	RefChromEndOfs=CurRefChromEndOfs;
	PrvBlockID = CurBlockID; 
	bFirst = false;
	}

if(bFirst) 
	{
	*pRefChromID  =   CurRefChromID;
	*pRefStrand   =   CurRefStrand;
	*pRefChromOfs =   CurRefChromOfs;
	*pRefAlignLen =   CurRefAlignLen;
	*pbLoaded = false;
	}
else
	{
	*pRefChromID    =RefChromID;
	*pRefStrand     =RefStrand;
	*pRefChromOfs   =RefChromOfs;
	*pRefAlignLen   =RefAlignLen;
	*pbLoaded = true;
	}
return(CurBlockID);
}



// NormaliseInDelColumns
// Because multialignments may be the result of merged alignments resulting in InDels being generated back into the 
// existing reference sequence then if only a subset of species are being processed there could be all InDels in any column
// of the subset sequences. In addition, where there are optional species then there could be eBaseUndefs in InDel columns
// This function will delete all columns which only contain InDels and/or eBaseUndef
// Returns the subset sequence length after eBaseInDel/eBaseUndef columns have been deleted
int
CGenHypers::NormaliseInDelColumns(tsProcParams *pProcParams,int AlignLen)
{
etSeqBase *pSeq;
etSeqBase *pSrcSeq;
etSeqBase SeqBase;
int NumVInDels;
int SeqIdx;
int VIdx;
int FirstInDel=-1;
int NormAlignLen=0;
for(SeqIdx = 0; SeqIdx < AlignLen; SeqIdx++)
	{
	NumVInDels = 0;					// assume no InDels in column
	for(VIdx = 0; VIdx < pProcParams->MaxAlignIdxSpecies; VIdx++)
		{
		pSeq = pProcParams->pSeqs[VIdx];
		SeqBase = pSeq[SeqIdx] & ~cRptMskFlg;
		if(SeqBase == eBaseInDel || SeqBase == eBaseUndef)
			NumVInDels++;
		else
			break;				
		}

	if(NumVInDels == pProcParams->MaxAlignIdxSpecies)	// all were InDels or undefined?
		{
		// mark col for deletion
		for(VIdx = 0; VIdx < pProcParams->MaxAlignIdxSpecies; VIdx++)
			{
			pSeq = pProcParams->pSeqs[VIdx];
			pSeq[SeqIdx] = (etSeqBase)0x0ff;
			}
		if(FirstInDel == -1)		// note idx of first InDel
			FirstInDel = SeqIdx;
		}
	else
		NormAlignLen++;				// accept this column
	}
if(NormAlignLen == AlignLen)	// if no columns to delete
	return(AlignLen);

// have at least one column which is all InDels and is to be deleted
for(VIdx = 0; VIdx < pProcParams->MaxAlignIdxSpecies; VIdx++)
	{
	pSeq = pProcParams->pSeqs[VIdx];
	pSeq = &pSeq[FirstInDel];
	pSrcSeq = pSeq;
	for(SeqIdx = FirstInDel; SeqIdx < AlignLen; SeqIdx++,pSrcSeq++)
		{
		if(*pSrcSeq != (etSeqBase)0x0ff)
			*pSeq++ = *pSrcSeq;
		}
	}
return(NormAlignLen);
}



// ProcAlignBlock
// Process an alignment block which may contain aligned subsequences meeting processing requirements
// Any subsequence (bounded by InDels if bInDelsAsMismatches, or end of block) of longer than MinLen is passed on
// to ProcessAlignment() which will process that subsequence looking for hyperconserved cores 
bool 
CGenHypers::ProcAlignBlock(int RefChromID,	    // reference chromosome
			   int RefChromOfs,		// offset along reference chromosome 
			   int AlignLen,		// alignment length incl any InDels
			   tsProcParams *pProcParams) // global processing parameters
{
int CurSeqLen;
int CurFiltSeqLen;	
int SubRefOfs;
int SubSeqStartIdx;
int NumIdents;
int SeqIdx;
int NumSubSeqs;
int CurNumSeqs;
int VIdx;
bool bAllIdentical; 
etSeqBase *pSeq;
etSeqBase RefBase;
etSeqBase SeqBase;
int NumVInDels;
bool bRefInDel;


if(pProcParams->bFiltLoConfidence && pProcParams->bFilt1stLast)
	NumSubSeqs = CUtility::GetNumSubseqs(AlignLen,		// alignment length incl InDels
							pProcParams->MinAlignSpecies,
							pProcParams->pSeqs);
else
	NumSubSeqs = 0;
CurNumSeqs=0;
CurSeqLen = 0;
CurFiltSeqLen = 0;
SubRefOfs = RefChromOfs;
SubSeqStartIdx = 0;

for(SeqIdx = 0; SeqIdx < AlignLen; SeqIdx++)
	{
	bAllIdentical = true;			// assume all identical in column
	NumVInDels = 0;					// assume no InDel in column
	bRefInDel = false;				// and specifically that the reference base is not an InDel 
	for(VIdx = 0; VIdx < pProcParams->MinAlignSpecies; VIdx++)
		{
		pSeq = pProcParams->pSeqs[VIdx];
		SeqBase = pSeq[SeqIdx] & ~cRptMskFlg;
		if(VIdx == pProcParams->RefSpeciesIdx)
			RefBase = SeqBase;
		if(SeqBase == eBaseInDel)
			{
			NumVInDels++;
			bAllIdentical = false;
			if(VIdx == pProcParams->RefSpeciesIdx)
				bRefInDel = true;
			}
		else
			if(SeqBase == eBaseN || SeqBase != RefBase)
				bAllIdentical = false;
		}
	// if all bases in column were InDels then alert user as this is a dataset processing error
	// should have been filtered out before this function is called..
	if(NumVInDels == pProcParams->MinAlignSpecies && pProcParams->MinAlignSpecies == pProcParams->NumSpeciesInAlignment)
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"Unexpected, all aligned bases were InDels at reference offset %d",RefChromOfs);

	// if in column any base was eBaseInDel then NumVInDels will be > 0
	// if the reference base was an InDel then bRefInDel will be true
	// if in column all identical and a,c,g or t then bAllIdentical will be true
	// if any base was mismatched or not a,c,g or t then bAllIdentical will be false
	if((pProcParams->PMode != eProcModeOutspecies && NumVInDels) && (!pProcParams->bInDelsAsMismatches &&
		(!pProcParams->bSloughRefInDels || (pProcParams->bSloughRefInDels && !bRefInDel))))
		{
		if(CurFiltSeqLen >= pProcParams->MinHyperCoreLen)	// has to be of at least MinLen to be worth further processing
			{
			CurNumSeqs+=1;
			if(CurFiltSeqLen >= pProcParams->MinSubSeqLen)			
				{
				if(pProcParams->bFiltLoConfidence)					// when filtering simply slough any low confidence subsequences
					{
					if(pProcParams->bFilt1stLast && (CurNumSeqs==1 || CurNumSeqs == NumSubSeqs))	// don't process first and last subsequence
						CurFiltSeqLen = 0;
					if(CurFiltSeqLen &&  pProcParams->MinIdent > 0 && ((NumIdents * 100)/CurFiltSeqLen) < pProcParams->MinIdent)
						CurFiltSeqLen = 0;
					}
				if(CurFiltSeqLen >= pProcParams->MinHyperCoreLen)
					ProcessAlignment(RefChromID,SubRefOfs,SubSeqStartIdx,CurFiltSeqLen,pProcParams);
				}
			}
		CurSeqLen = 0;
		CurFiltSeqLen = 0;
		if(RefBase != eBaseInDel)
			RefChromOfs++;
		NumIdents = 0;
		SubSeqStartIdx = SeqIdx + 1;	// mark where next subsequence could start 
		SubRefOfs = RefChromOfs;		// chromosomal offset at which that next subsequence starts
		continue;
		}

	// bInDelsAsMismatches or no InDels in current aligned column - could still be mismatches - ,
	// continue to accept column as part of subsequence
	if(CurSeqLen == 0)		// if first base of putative subsequence...
		{
		if(pProcParams->bFiltLoConfidence && !bAllIdentical) // when filtering, must start on an identical base
			{
			if(RefBase != eBaseInDel)		// advance chromosomal offset  only if not a InDel
				RefChromOfs++;
			SubRefOfs = RefChromOfs;		// chromosomal offset at which that next subsequence starts
			continue;
			}
		NumIdents = 0;
		SubSeqStartIdx = SeqIdx;			// mark where this subsequence has started 
		SubRefOfs = RefChromOfs;			// chromosomal offset at which this subsequence started
		}
	CurSeqLen++;
	if(bAllIdentical)						// mark last identical...
		NumIdents++;
	CurFiltSeqLen = CurSeqLen;
	if(RefBase != eBaseInDel)				// advance chromosomal offset  only if not a InDel
		RefChromOfs++;
	}

// no more bases in this block
if(CurFiltSeqLen >= pProcParams->MinHyperCoreLen)	// has to be at least MinLen to be worth further processing
	{
	if(CurFiltSeqLen >= pProcParams->MinSubSeqLen)			
		{
		if(pProcParams->bFiltLoConfidence) // when filtering simply slough any low confidence subsequences
			{
			if(pProcParams->bFilt1stLast && (!CurNumSeqs || CurNumSeqs == (NumSubSeqs-1)))	// don't process first and last subsequence
				CurFiltSeqLen = 0;
			if(CurFiltSeqLen &&  pProcParams->MinIdent > 0 && ((NumIdents * 100)/CurFiltSeqLen) < pProcParams->MinIdent)
				CurFiltSeqLen = 0;
			}
		if(CurFiltSeqLen >= pProcParams->MinHyperCoreLen)
			ProcessAlignment(RefChromID,SubRefOfs,SubSeqStartIdx,CurFiltSeqLen,pProcParams);
		}
	}

return(true);
}

// ProcAlignBlockSummary
// Process an alignment block which may contain aligned subsequences meeting processing requirements
// Generates stats on the sequences within this alignment block
// Stats: Total number of mismatches
//		  Total number of exact matches
//		  Total numer of ref InDels
//		  Total number of rel InDels
// Above is characterised by region
bool 
CGenHypers::ProcAlignBlockSummary(int RefChromID,	    // reference chromosome
			   int RefChromOfs,		// offset along reference chromosome 
			   int AlignLen,		// alignment length incl any InDels
			   tsProcParams *pProcParams) // global processing parameters
{
int CurSeqLen;
int CurFiltSeqLen;	
int SubRefOfs;
int SubSeqStartIdx;
int SeqIdx;
int NumSubSeqs;
int CurNumSeqs;
int VIdx;

etSeqBase *pSeq;
etSeqBase RefBase;
etSeqBase RelBase;
etSeqBase SeqBase;
int NumInDels;
int NumMatching;

if(pProcParams->bFiltLoConfidence && pProcParams->bFilt1stLast)
	NumSubSeqs = CUtility::GetNumSubseqs(AlignLen,		// alignment length incl InDels
							pProcParams->MinAlignSpecies,
							pProcParams->pSeqs);
else
	NumSubSeqs = 0;
CurNumSeqs=0;
CurSeqLen = 0;
CurFiltSeqLen = 0;
SubRefOfs = RefChromOfs;
SubSeqStartIdx = 0;

for(SeqIdx = 0; SeqIdx < AlignLen; SeqIdx++)
	{
	NumInDels = 0;
	NumMatching = 1;	// 1st base always matches itself!
	for(VIdx = 0; VIdx < pProcParams->MinAlignSpecies; VIdx++)
		{
		pSeq = pProcParams->pSeqs[VIdx];
		SeqBase = pSeq[SeqIdx] & ~cRptMskFlg;
		if(VIdx == pProcParams->RefSpeciesIdx)
			RefBase = SeqBase;

		if(SeqBase == eBaseInDel)
			NumInDels++;
		if(!VIdx)
			RelBase = SeqBase;
		else
			if(RelBase == SeqBase)
				NumMatching++;
		}

	// if all bases in column were InDels then alert user as this is a dataset processing error
	// should have been filtered out before this function is called..
	if(NumInDels ==  pProcParams->MinAlignSpecies && pProcParams->MinAlignSpecies == pProcParams->NumSpeciesInAlignment)
		{
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"Unexpected, all aligned bases were InDels at reference offset %d",RefChromOfs);
		return(false);
		}
	if(RefBase == eBaseInDel)
		{
		ReportSummary(RefChromID,RefChromOfs,0,pProcParams);
		continue;	// don't advance chromofs
		}
	if(NumInDels)	// ref wasn't InDel so if any InDels then they must be in a relative species
		ReportSummary(RefChromID,RefChromOfs,1,pProcParams);
	else
		{
		if((RefBase <= eBaseT) && NumMatching == pProcParams->MinAlignSpecies)		
			ReportSummary(RefChromID,RefChromOfs,2,pProcParams); // match
		else						
			ReportSummary(RefChromID,RefChromOfs,3,pProcParams); // at least one mismatch
		}
	RefChromOfs++;
	}
return(true);
}

int
CGenHypers::ReportSummary(int RefChromID,int RefChromOfs,int ProcMode,tsProcParams *pProcParams)
{
int FeatureBits;
int RegionIdx;
int BitMsk;
int SpliceSiteOverlaps;
int FeatIdx;
int *pStep = pProcParams->pCntStepCnts;

// ensure that the hypercore is in an included region and not part of an excluded region
if(!IncludeFilter(RefChromID,RefChromOfs,RefChromOfs,pProcParams))
		return(eBSFSuccess);

pStep += ProcMode * pProcParams->Regions;
	
if(pProcParams->pBiobed != nullptr)
	{
	if(pProcParams->BEDChromID > 0)
		FeatureBits = pProcParams->pBiobed->GetFeatureBits(pProcParams->BEDChromID,RefChromOfs,RefChromOfs,cRegionFeatBits,pProcParams->UpDnStreamLen);
	else
		FeatureBits = 0;
	RegionIdx = 0;		// default to intergenic if no feature bits set
	if(FeatureBits)		
		{
		BitMsk = cFeatBitCDS;
		for(FeatIdx = 1; FeatIdx < 7; FeatIdx++,BitMsk <<= 1)
			{
			if(BitMsk & FeatureBits)
				{
				if(RegionIdx)		// if already have feature
					return(eBSFSuccess);	// although was sequence of interest, more than one feature bit so can't contribute to stats
				RegionIdx = FeatIdx;
				if(pProcParams->bMultipleFeatBits)	// if multiple features allowed then don't check for any additional
					break;
				}
			}
		}

	// need to remap RegionIdx when incr counts so regions are displayed in a 5'->3' logical order
	switch(RegionIdx) {
		case 0: case 2: case 4: case 6:		// IG,5'UTR, Introns and 3'DS
			pStep[RegionIdx] += 1;	
			break;
		case 1:								// CDS
			pStep[3] += 1;
			break;
		case 3:								// 3'UTR
			pStep[5] += 1;
			break;
		case 5:								// 5'US
			pStep[1] += 1;
			break;
		}
			
	if(RegionIdx != 0)		// if not intergenic then check for splice sites
		{
		SpliceSiteOverlaps = pProcParams->pBiobed->GetSpliceSiteBits(pProcParams->BEDChromID,RefChromOfs,RefChromOfs,cMinSpliceOverlap);
		if(SpliceSiteOverlaps & cIntronExonSpliceSite)
			pStep[7]++;
		if(SpliceSiteOverlaps & cExonIntronSpliceSite)
			pStep[8]++;
		}
	pProcParams->bStatsAvail = true;
	ChkOutputSummaryResults(pProcParams->szRefChrom,RefChromOfs,pProcParams,false,false);	// output results as may be appropriate
	}
else
	{
	*pStep += 1;
	ChkOutputSummaryResults(pProcParams->szRefChrom,RefChromOfs,pProcParams,false,false);	// output results as may be appropriate
	}
return(eBSFSuccess);
}



// OpenBedfile
// Attempts to open specified bedfile
// Returns ptr to opened bedfile or nullptr
CBEDfile *
CGenHypers::OpenBedfile(char *pToOpen)
{
int Rslt;
CBEDfile *pBed;
if(pToOpen != nullptr && pToOpen[0] != '\0')
	{
	if((pBed = (CBEDfile *)new CBEDfile())==nullptr)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBEDfile '%s'",pToOpen);
		return(nullptr);
		}

	if((Rslt = pBed->Open(pToOpen,eBTAnyBed))!=eBSFSuccess)
		{
		while(pBed->NumErrMsgs())
			gDiagnostics.DiagOut(eDLFatal,gszProcName,pBed->GetErrMsg());
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open biobed file %s",pToOpen);
		delete pBed;
		return(nullptr);
		}
	return(pBed);
	}
return(nullptr);
}

//CloseBedfiles
//Closes and deletes all created and opened Biobed files
bool
CGenHypers::CloseBedfiles(tsProcParams *pProcParams)
{
int Idx;
for(Idx=0; Idx < pProcParams->NumIncludes; Idx++)
	{
	if(pProcParams->pIncludes[Idx] != nullptr)
		delete pProcParams->pIncludes[Idx];
	pProcParams->pIncludes[Idx] = nullptr;
	}
for(Idx=0; Idx < pProcParams->NumExcludes; Idx++)
	{
	if(pProcParams->pExcludes[Idx] != nullptr)
		delete pProcParams->pExcludes[Idx];
	pProcParams->pExcludes[Idx] = nullptr;
	}
if(pProcParams->pBiobed != nullptr)
	delete pProcParams->pBiobed;
pProcParams->pBiobed = nullptr;
return(true);
}

// IncludeFilter
// Returns true if subsequence can be processed or false if subsequence is not to be processed
// To be processed a subsequence -
// A) If no Include BED files specified then subsequence is assumed included unless excluded by following rule B).
//    If at least one Include BED file was specified then if subsequence not in any of the Include files then false is returned
// B) If no Exclude files specified then subsequence is returned as Ok (true) to process. If subsequence not in any of the Exclude files
//    then subsequence is returned as Ok (true) to process, otherwise false is returned
bool
CGenHypers::IncludeFilter(int RefChromID,int SubRefOfs,int SubRefEndOfs,tsProcParams *pProcParams)
{
int Idx;
int BEDChromID;
if(pProcParams->NumIncludes)
	{
	for(Idx = 0; Idx < pProcParams->NumIncludes; Idx++)
		{
		if(pProcParams->pIncludes[Idx] == nullptr) // should'nt ever be nullptr but...
			continue;
		if((BEDChromID = pProcParams->pIncludes[Idx]->LocateChromIDbyName(pProcParams->szRefChrom))<1)
			continue;
		if(pProcParams->pIncludes[Idx]->InAnyFeature(BEDChromID,SubRefOfs,SubRefEndOfs))
			break;
		}
	if(Idx == pProcParams->NumIncludes)
		return(false);
	}

if(pProcParams->NumExcludes)
	{
	for(Idx = 0; Idx < pProcParams->NumExcludes; Idx++)
		{
		if(pProcParams->pExcludes[Idx] == nullptr) // should'nt ever be nullptr but...
			continue;
		if((BEDChromID = pProcParams->pExcludes[Idx]->LocateChromIDbyName(pProcParams->szRefChrom))<1)
			continue;
		if(pProcParams->pExcludes[Idx]->InAnyFeature(BEDChromID,SubRefOfs,SubRefEndOfs))
			return(false);
		}
	}
return(true);
}

// Process
// Create alignment stats from files in specified source directory
int
CGenHypers::Process(bool bTargDeps,				// true if process only if any independent src files newer than target
	int PMode,					// processing mode 0: default, 1: summary stats only
	int NumBins,				// when generating length distributions then use this many bins - 0 defaults to using 1000
	int BinDelta,				// when generating length distributions then each bin holds this length delta - 0 defaults to auto determine from NunBins and longest sequence length
	char* pszInputFile,		// bio multialignment (.algn) file to process
	char* pszOutputFile,	// where to write out stats
	char* pszOutputCoreFile, // where to write out the hypercore loci 
	char* pszBiobedFile,	// biobed file containing regional features - exons, introns etc
	int NumIncludeFiles,	// number of include region files
	char** ppszIncludeFiles,// biobed files containing regions to include - default is to exclude none
	int NumExcludeFiles,	// number of exclude region files
	char** ppszExcludeFiles,// biobed file containing regions to exclude - default is to include all
	char* pszUniqueSpeciesSeqs, // ignore alignment block if these species sequences are not unique
	int WindowSize,			// sampling window size
	int NumCoreSpecies,		// number of core species to be in alignment
	char* pszCoreSpecies,	// space or comma separated list of species which must be present in any MAF block alignment before that block will be further processed
	int MinNonCoreSpecies,	// minimum number of species required in an alignment (excluding core species)
	int MinHyperCoreLen,		// minimum hyper core length required
	int MaxHyperColsMismatches,	// hyper cores can have up to this number of columns with at least one mismatches
	int VMismatches,		// number of mismatches in an alignment col to count as a mismatch against MaxHyperColsMismatches
	int MinIdentity,		// minimum identity required when processing hyperconserved
	int NumDistSegs,		// number of match distribution profile segments
	int RegLen,				// regulatory region length - up/dn stream of 5/3' 
	bool bMultipleFeatBits,	// if true then accept alignments in which multiple feature bits are set
	bool bInDelsAsMismatches, // treat InDels as if mismatches
	bool bSloughRefInDels,	// slough columns in which the reference base is an InDel
	bool bFiltLoConfidence,	// filter out low confidence subsequences
	bool bFilt1stLast,		// treat 1st and last subsequences as being low confidence
	int MinSubSeqLen,		// subsequences of less than this length are treated as being low confidence
	int MinIdent,			// treat subsequences of less than this identity as being low confidence
	int NumIncludeChroms,	// number of chromosomes explicitly defined to be included
	char** ppszIncludeChroms,	// ptr to array of reg expressions defining chroms to include - overides exclude
	int NumExcludeChroms,	// number of chromosomes explicitly defined to be excluded
	char** ppszExcludeChroms)	// ptr to array of reg expressions defining chroms to exclude
{
int Rslt;
int Idx;
int *pCntStepCnts;
int NumCnts;
int Regions;
int UpDnStreamLen;
CBEDfile *pBiobed = nullptr;
tsProcParams ProcParams;
int	MismatchScore;
int	MatchScore;
char *pszCSVSpecies;				// to hold comma separated species list
char szInaccessible[_MAX_PATH];

Reset();

if(bTargDeps && (Rslt = CUtility::Chk2TargDepend(szInaccessible,_MAX_PATH,pszOutputFile,pszOutputCoreFile,pszBiobedFile,nullptr)) <= 0)
	{
	if(Rslt)
		{
		gDiagnostics.DiagOutMsgOnly(eDLFatal,"Unable to access source file '%s'",szInaccessible);
		return(Rslt);
		}
	for(Idx=0;Idx<NumIncludeFiles; Idx++)
		{
		Rslt = CUtility::Chk2TargDepend(szInaccessible,_MAX_PATH,pszOutputFile,pszOutputCoreFile,ppszIncludeFiles[Idx],nullptr);
		if(Rslt < 0)
			{
			gDiagnostics.DiagOutMsgOnly(eDLFatal,"Unable to access source include file '%s'",szInaccessible);
			return(Rslt);
			}
		}

	for(Idx=0;Idx<NumExcludeFiles; Idx++)
		{
		Rslt = CUtility::Chk2TargDepend(szInaccessible,_MAX_PATH,pszOutputFile,pszOutputCoreFile,ppszExcludeFiles[Idx],nullptr);
		if(Rslt < 0)
			{
			gDiagnostics.DiagOutMsgOnly(eDLFatal,"Unable to access source exclude file '%s'",szInaccessible);
			return(Rslt);
			}
		}
	}

#ifndef _WIN32
int RegErr;				// regular expression parsing error
char szRegErr[128];			// to hold RegErr as textual representation ==> regerror();
#endif

memset(&ProcParams,0,sizeof(tsProcParams));

// parse out species which must have unique alignment block sequences
if(!stricmp(pszUniqueSpeciesSeqs,"*"))
	ParseUniqueSpeciesSeqs(pszCoreSpecies,&ProcParams);
else
	ParseUniqueSpeciesSeqs(pszUniqueSpeciesSeqs,&ProcParams);

// parse out species list
ProcParams.NumSpeciesList = ParseNumSpecies(pszCoreSpecies,&ProcParams);
pszCSVSpecies = new char[strlen(pszCoreSpecies) * 2];
pszCSVSpecies[0]='\0';
for(Idx = 0; Idx < ProcParams.NumSpeciesList; Idx++)
	{
	if(Idx > 0)
		strcat(pszCSVSpecies,",");
	strcat(pszCSVSpecies,ProcParams.szSpecies[Idx]);
	}
ProcParams.pszSpeciesList = pszCSVSpecies;

for(Idx=0;Idx<NumIncludeFiles; Idx++)
	{
	if((ProcParams.pIncludes[Idx] = OpenBedfile(ppszIncludeFiles[Idx]))==nullptr)
		{
		CloseBedfiles(&ProcParams);
		delete []pszCSVSpecies;
		return(eBSFerrObj);
		}
	ProcParams.NumIncludes++;
	}

for(Idx=0;Idx<NumExcludeFiles; Idx++)
	{
	if((ProcParams.pExcludes[Idx] = OpenBedfile(ppszExcludeFiles[Idx]))==nullptr)
		{
		CloseBedfiles(&ProcParams);
		delete[]pszCSVSpecies;
		return(eBSFerrObj);
		}
	ProcParams.NumExcludes++;
	}


if(pszBiobedFile != nullptr && pszBiobedFile[0] != '\0')
	{
	if((ProcParams.pBiobed = OpenBedfile(pszBiobedFile))==nullptr)
		{
		CloseBedfiles(&ProcParams);
		delete[]pszCSVSpecies;
		return(eBSFerrObj);
		}
	Regions = cNumCharRegs;	
	UpDnStreamLen = RegLen;
	}
else
	{
	Regions = 1;
	UpDnStreamLen = 0;
	ProcParams.pBiobed = nullptr;
	}

m_NumLenRangeBins = NumBins;
m_BinDelta = BinDelta;

if(PMode == eProcModeSummary)
	NumCnts=Regions * 4;	// allows for RefIndels,RelInDels,Matches and Mismatch counts per region
else
	NumCnts=Regions * m_NumLenRangeBins;

if((pCntStepCnts = new int[NumCnts])==nullptr)
	{
	CloseBedfiles(&ProcParams);
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory (%d bytes) for holding alignment statistics",sizeof(int) * NumCnts);
	delete[]pszCSVSpecies;
	return(eBSFerrMem);
	}
memset(pCntStepCnts,0,NumCnts * sizeof(int));

InitLengthRangeClass(NumBins,BinDelta);

ProcParams.bInDelsAsMismatches = bInDelsAsMismatches;
ProcParams.bSloughRefInDels = bSloughRefInDels;			// slough columns in which the reference base is an InDel

ProcParams.bFiltLoConfidence = bFiltLoConfidence;
ProcParams.bFilt1stLast = bFilt1stLast;
ProcParams.MinIdent = MinIdent;
ProcParams.MinSubSeqLen = MinSubSeqLen;
ProcParams.PMode = PMode;
ProcParams.MinAlignSpecies = NumCoreSpecies + MinNonCoreSpecies;
ProcParams.NumCoreSpecies = NumCoreSpecies;
ProcParams.MinNonCoreSpecies = MinNonCoreSpecies;
ProcParams.NumSpeciesInAlignment = 0;
ProcParams.MaxAlignIdxSpecies = 0;
ProcParams.RefSpeciesIdx = 0;
ProcParams.pCntStepCnts = pCntStepCnts;
ProcParams.NumDistSegs = NumDistSegs;
ProcParams.UpDnStreamLen = UpDnStreamLen;
ProcParams.bMultipleFeatBits = bMultipleFeatBits;
ProcParams.Regions = Regions;
ProcParams.NumCnts = NumCnts;
ProcParams.MinHyperCoreLen = MinHyperCoreLen;

ProcParams.MaxHyperColsMismatches = MaxHyperColsMismatches;
ProcParams.MinIdentity = MinIdentity;		// minimum identity required when processing hyperconserved

if(PMode != eProcModeSummary && MinIdentity < 100)
	{
	MismatchScore = (cRandWalk100Score - 1) / (100 - MinIdentity);	// decrease score by this much on each mismatch
	MatchScore = cRandWalk100Score / MinIdentity;					// increase score by this much each time there is a match
	}
else
	{
	MismatchScore = 0;		// 100% identity is an ultra!
	MatchScore = 0;
	}

ProcParams.MismatchScore = MismatchScore;
ProcParams.MatchScore = MatchScore;
ProcParams.VMismatches = VMismatches;

if((Rslt =m_RegExprs.CompileREs(NumIncludeChroms, ppszIncludeChroms,NumExcludeChroms, ppszExcludeChroms)) < eBSFSuccess)
	{
	delete pCntStepCnts;
	CloseBedfiles(&ProcParams);
	delete[]pszCSVSpecies;
	return(eBSFerrMem);
	}


// determine max aligned sequence length for any single species which can be handled
ProcParams.MaxSeqAlignLen = cMAtotMaxSeqAlignLen/ProcParams.NumSpeciesList;
for(Idx = 0; Idx < ProcParams.NumSpeciesList; Idx++)
	{
	if((ProcParams.pSeqs[Idx] = new uint8_t [ProcParams.MaxSeqAlignLen])==nullptr)
		{
		delete pCntStepCnts;
		CloseBedfiles(&ProcParams);
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory (%d bytes) for holding species sequences",ProcParams.MaxSeqAlignLen);
		delete[]pszCSVSpecies;
		return(eBSFerrMem);
		}
	}


ProcParams.hCoreCSVRsltsFile = -1;
Rslt = GenKimura3P(pszInputFile, pszOutputFile, &ProcParams);
if (Rslt == eBSFSuccess)
	{
	if(pszOutputCoreFile != nullptr && pszOutputCoreFile[0] != '\0')
		{
#ifdef _WIN32
		if((ProcParams.hCoreCSVRsltsFile = open(pszOutputCoreFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
		if((ProcParams.hCoreCSVRsltsFile = open(pszOutputCoreFile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create %s - %s",pszOutputCoreFile,strerror(errno));
			close(ProcParams.hRsltsFile);
			delete []pCntStepCnts;
			CloseBedfiles(&ProcParams);
			delete[]pszCSVSpecies;
			return(eBSFerrCreateFile);
			}
		}
	else
		ProcParams.hCoreCSVRsltsFile = -1;
	}

if(Rslt == eBSFSuccess)
	Rslt = ProcessAlignments(pszInputFile,&ProcParams);
if(Rslt == eBSFSuccess && !WindowSize)
	{
	if(PMode == eProcModeSummary)
		OutputSummaryResults( "Alignment",0,&ProcParams,false);
	else
		OutputResults( "Genome",0,&ProcParams,false);
	}
if(pCntStepCnts != nullptr)
	delete []pCntStepCnts;
if(ProcParams.hRsltsFile != -1)
	{
	close(ProcParams.hRsltsFile);
	ProcParams.hRsltsFile = -1;
	}
if(ProcParams.hCoreCSVRsltsFile != -1)
	{
	close(ProcParams.hCoreCSVRsltsFile);
	ProcParams.hCoreCSVRsltsFile = -1;
	}
CloseBedfiles(&ProcParams);
for(Idx = 0; Idx < ProcParams.NumSpeciesList; Idx++)
	{
	if(ProcParams.pSeqs[Idx] != nullptr)
		{
		delete ProcParams.pSeqs[Idx];
		ProcParams.pSeqs[Idx] = nullptr;
		}
	}
delete[]pszCSVSpecies;
return(Rslt);
}



// process alignment
// Process an alignment subsequence which may contain mismatches or 'N's but will only contain InDels if bInDelsAsMismatches is true 
// The subsequence length will be of at least MinLen long
bool
CGenHypers::ProcessAlignment(int RefChromID,		// reference chromosome
				 int ChromOfs,			// chromosome offset (0..n) at which this subsequence started
				 int SeqIdx,			// index (0..n) into pProcParms->pSeqs[] which alignment subsequence starts	
				 int SubSeqLen,			// subsequence length
				 tsProcParams *pProcParams)
{
int NxtSeqIdx;
etSeqBase *pRefBase;
pRefBase = pProcParams->pSeqs[pProcParams->RefSpeciesIdx];
pRefBase += SeqIdx;
while(SubSeqLen >= pProcParams->MinHyperCoreLen)
	{
		// from current SubRefOfs maximally extend to right allowing for MaxMismatches
	NxtSeqIdx = ProcessSubSeq(RefChromID,ChromOfs,SeqIdx,SubSeqLen,pProcParams);
	if(NxtSeqIdx == -1)		// -1 flags that there is no more processing required for this alignment
		break;
	SubSeqLen -= (NxtSeqIdx - SeqIdx);	// update length of remaining sequence in this alignment
	if(SubSeqLen < pProcParams->MinHyperCoreLen)
		break;

	// determine chrom offset at which next
	while(SeqIdx < NxtSeqIdx)
		{
		if((*pRefBase++ & ~cRptMskFlg) != eBaseInDel)
			ChromOfs++;
		SeqIdx++;
		}
	}
return(true);
}

// Reference: Kimura 3 Parameter distance https://en.wikipedia.org/wiki/Genetic_distance
// For each species. relative to the reference species, empirically determine the transistion/transversion1/transversion2 mutation rates.
// Transition:		A<->G,T<->C 
// Transversion1:	T<->A,C<->G
// Transversion2:	T<->G,A<->C
// For each species iterate over all alignment blocks, containing that species, and accumulate counts of each observed ref:species relative nucleotides
// Thus a 4x4 matrix of counts for each species is generated
typedef struct Tag_sKimura3PCnts {
	int32_t SpeciesID;		// counts are for this species
	int64_t NumMatches;		// ref:rel exacts
	int64_t NumTransitions; // ref:rel Transition	A<->G,T<->C
	int64_t NumTransversion1s; // ref:rel Transversion1:	T<->A,C<->G
	int64_t NumTransversion2s; // ref:rel Transversion2:	T<->G,A<->C
	int64_t NumIndeterminate;		// ref:rel is indeterminate - must likely no alignment or a short InDel
	double Distance;		// Kimura3P distance
} tsKimura3PCnts;

int
CGenHypers::GenKimura3P(char* pszMAF,			 // source bioseq multialignment file
						char* pszKimura3PFile,	// write species parameter counts to this output file
	tsProcParams* pProcParams) // processing parameters
{
int RefSpeciesID;
int RefChromID;
int BEDChromID;
int PrevRefChromID;
int PrevDispRefChromID;
char* pszRefChrom;
int RefChromOfs;
char RefStrand;
int SpeciesIdx2IDs[cMaxAlignedSpecies];  // mapping species indexes to species identifiers
int SpeciesID2Idxs[cMaxAlignedSpecies+1];	// mapping species identifiers to species indexes
char szOutFile[_MAX_PATH];
int Rslt;
int RefAlignLen;
CMAlignFile* pAlignments;
int Idx;
int CurBlockID;
bool bLoaded;
tsKimura3PCnts *pKimura3PCnts;
tsKimura3PCnts* pKimura3PCnt;
// compose CSV file name, the file will contain the Kimura3P global genome parameter values
CUtility::AppendFileNameSuffix(szOutFile, pszKimura3PFile, (char*)".Kimura3P.csv", '.');

if (m_pszOutBuffer == nullptr)
{
	m_pszOutBuffer = new uint8_t[cAllOutBuffSize];
	m_AllocOutBuff = cAllOutBuffSize;
}
m_OutBuffIdx = 0;

#ifdef _WIN32
if ((pProcParams->hRsltsFile = open(szOutFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE)) == -1)
#else
if ((pProcParams->hRsltsFile = open(szOutFile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE)) == -1)
#endif
{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenKimura3P: Unable to open or create %s - %s", szOutFile, strerror(errno));
	return(eBSFerrCreateFile);
	}
m_OutBuffIdx = sprintf((char*)m_pszOutBuffer, "\"Species\",\"NumExacts\",\"NumTransitions\",\"NumTransversion1s\",\"NumTransversion2s\",\"NumIndeterminate\"\n");
if (!CUtility::RetryWrites(pProcParams->hRsltsFile, m_pszOutBuffer, m_OutBuffIdx))
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenKimura3P: Fatal error in '%s' RetryWrites()", szOutFile);
	close(pProcParams->hRsltsFile);
	pProcParams->hRsltsFile = -1;
	Reset();
	return(eBSFerrFileAccess);
	}
m_OutBuffIdx = 0;

if ((pAlignments = new CMAlignFile()) == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenKimura3P: Unable to create new instance of CMAlignFile");
	return(eBSFerrObj);
	}

if ((Rslt = pAlignments->Open(pszMAF)) != eBSFSuccess)
	{
	while (pAlignments->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal, gszProcName, pAlignments->GetErrMsg());
	delete pAlignments;
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenKimura3P: Unable to open MAF file %s\n", pszMAF);
	return(Rslt);
	}

// total number of species in the multialignments?
m_NumMAFSpecies = pAlignments->GetNumSpecies();
RefSpeciesID = pAlignments->GetRefSpeciesID();			// get reference species identifer

pKimura3PCnts = new tsKimura3PCnts[m_NumMAFSpecies];
memset(pKimura3PCnts, 0, sizeof(tsKimura3PCnts) * m_NumMAFSpecies);
pKimura3PCnt = pKimura3PCnts;
	// obtain identifiers for each of the species present in the multialignments
for (Idx = 0; Idx < m_NumMAFSpecies; Idx++, pKimura3PCnt++)
	{
	if ((Rslt = SpeciesIdx2IDs[Idx] = pAlignments->LocateSpeciesID(Idx)) < 1)	// mapping species index to species identifier
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenKimura3P: Species indexed as '%d' not represented in %s", Idx, pszMAF);
		delete pAlignments;
		return(eBSFerrEntry);
		}
	SpeciesID2Idxs[Rslt] = Idx;
	pKimura3PCnt->SpeciesID = Rslt;
	}


pAlignments->SetAllConfFilt(false);

	// iterate over reference blocks which are sorted by chrom then offset
CurBlockID = 0;
PrevRefChromID = 0;
PrevDispRefChromID = 0;
BEDChromID = 0;
pProcParams->RefSpeciesIdx = 0;		// reference sequence will always be 1st
pProcParams->NxtOutputOffset = pProcParams->WindowSize;

while (CurBlockID >= 0 && ((CurBlockID =						// returned blockid to next start loading from
		LoadContiguousBlocks(RefSpeciesID,	// reference species identifier
			CurBlockID,			// which block to initially start loading from
			&bLoaded,			// returned indicator as to if any loaded blocks meet processing requirements
			&RefChromID,		// returned reference chromosome identifier 
			&RefStrand,			// returned reference strand
			&RefAlignLen,		// returned alignment (incl InDels) length
			&RefChromOfs,		// returned alignment start offset
			SpeciesIdx2IDs,		// input - species of interest identifier array
			pAlignments,
			pProcParams)) > 0 || (CurBlockID == eBSFerrAlignBlk && RefAlignLen > 0)))
	{
	if (RefChromID > 0 && RefChromID != PrevDispRefChromID)
		{
		PrevDispRefChromID = RefChromID;
		pszRefChrom = pAlignments->GetChromName(RefChromID);

		strcpy(pProcParams->szRefChrom, pszRefChrom);
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenKimura3P: Processing chromosome %s", pszRefChrom);
		pProcParams->NxtOutputOffset = pProcParams->WindowSize;
		if (pProcParams->pBiobed != nullptr)
			pProcParams->BEDChromID = pProcParams->pBiobed->LocateChromIDbyName(pProcParams->szRefChrom);
		else
			pProcParams->BEDChromID = 0;
		pProcParams->RefChromID = RefChromID;
		}

	if (!bLoaded)
		continue;

	int32_t AlignColIdx;
	int32_t SpeciesSeqIdx;
	etSeqBase RefBase;
	etSeqBase *pSeqBases;
	etSeqBase SeqBase;
	for (AlignColIdx = 0; AlignColIdx < RefAlignLen; AlignColIdx++)
		{
		pKimura3PCnt = pKimura3PCnts;
		RefBase = eBaseUndef;
		for (SpeciesSeqIdx = 0; SpeciesSeqIdx < m_NumMAFSpecies; SpeciesSeqIdx++, pKimura3PCnt++)
			{
			pSeqBases = pProcParams->pSeqs[SpeciesSeqIdx];
			SeqBase = pSeqBases[AlignColIdx] & ~cRptMskFlg;
			if (SpeciesSeqIdx == pProcParams->RefSpeciesIdx)
				{
				if(SeqBase > eBaseT)		// skipping over cols where the reference contains a non-canonical base
					break;
				RefBase = SeqBase;
				pKimura3PCnt->NumMatches++;
				continue;
				}
			if(RefBase == eBaseUndef)
				break;

			if (SeqBase > eBaseT)
				{
				pKimura3PCnt->NumIndeterminate++;
				continue;
				}
			switch (RefBase) {
				case eBaseA:
					switch (SeqBase) {
						case eBaseA:
							pKimura3PCnt->NumMatches++;
							break;

						case eBaseC:
							pKimura3PCnt->NumTransversion2s++;
							break;

						case eBaseG:
							pKimura3PCnt->NumTransitions++;
							break;

						case eBaseT:
							pKimura3PCnt->NumTransversion1s++;
							break;
						}
					break;

				case eBaseC:
					switch (SeqBase) {
						case eBaseA:
							pKimura3PCnt->NumTransversion2s++;
							break;

						case eBaseC:
							pKimura3PCnt->NumMatches++;
							break;

						case eBaseG:
							pKimura3PCnt->NumTransversion1s++;
							break;

						case eBaseT:
							pKimura3PCnt->NumTransitions++;
							break;
						}
					break;

				case eBaseG:
					switch (SeqBase) {
						case eBaseA:
							pKimura3PCnt->NumTransitions++;
							break;

						case eBaseC:
							pKimura3PCnt->NumTransversion1s++;
							break;

						case eBaseG:
							pKimura3PCnt->NumMatches++;
							break;

						case eBaseT:
							pKimura3PCnt->NumTransversion2s++;
							break;
						}
					break;

				case eBaseT:
					switch (SeqBase) {
					case eBaseA:
						pKimura3PCnt->NumTransversion1s++;
						break;

					case eBaseC:
						pKimura3PCnt->NumTransitions++;
						break;

					case eBaseG:
						pKimura3PCnt->NumTransversion2s++;
						break;

					case eBaseT:
						pKimura3PCnt->NumMatches++;
						break;
					}
				break;
				}
			}
		}

	if (RefChromID != PrevRefChromID)
		PrevRefChromID = RefChromID;
	}


pKimura3PCnt = pKimura3PCnts;
for(Idx = 0; Idx < m_NumMAFSpecies; Idx++, pKimura3PCnt++)
	{
	char *pszAlignedSpecies = pAlignments->GetSpeciesName(pKimura3PCnt->SpeciesID);
	m_OutBuffIdx += sprintf((char *)&m_pszOutBuffer[m_OutBuffIdx], "\"%s\",%zd,%zd,%zd,%zd,%zd,%f\n",
								pszAlignedSpecies, pKimura3PCnt->NumMatches, pKimura3PCnt->NumTransitions, pKimura3PCnt->NumTransversion1s, pKimura3PCnt->NumTransversion2s, pKimura3PCnt->NumIndeterminate, pKimura3PCnt->Distance);
	if ((m_OutBuffIdx + 1000) > m_AllocOutBuff)
	{	
		if (!CUtility::RetryWrites(pProcParams->hRsltsFile, m_pszOutBuffer, m_OutBuffIdx))
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenKimura3P: Fatal error in '%s' RetryWrites()", szOutFile);
			close(pProcParams->hRsltsFile);
			pProcParams->hRsltsFile = -1;
			Reset();
			return(eBSFerrFileAccess);
			}
		m_OutBuffIdx = 0;
		}
	}

if (m_OutBuffIdx > 0 && pProcParams->hRsltsFile != -1)
	{
	if (!CUtility::RetryWrites(pProcParams->hRsltsFile, m_pszOutBuffer, m_OutBuffIdx))
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenKimura3P: Fatal error in '%s' RetryWrites()", szOutFile);
		close(pProcParams->hRsltsFile);
		pProcParams->hRsltsFile = -1;
		Reset();
		return(eBSFerrFileAccess);
		}
	m_OutBuffIdx = 0;
	}

if (pProcParams->hRsltsFile != -1)
	{
	// commit output file
#ifdef _WIN32
	_commit(pProcParams->hRsltsFile);
#else
	fsync(pProcParams->hRsltsFile);
#endif
	close(pProcParams->hRsltsFile);
	pProcParams->hRsltsFile = -1;
	}

delete pAlignments;
delete []pKimura3PCnts;
return(eBSFSuccess);
}





 
// Process subsequence which starts immediately following the end of the previous subsequence of start of an alignment block
// and which ends immediately prior to the next InDel - if bInDelsAsMismatches false -  or end of alignment block
// Returns the SubSeqStartIdx at which to start a subsequent ProcessSubSeq processing call
int	
CGenHypers::ProcessSubSeq(int RefChromID,		// reference chromosome identifier
			  int ChromOfs,			// chromosome offset (0..n) at which this subsequence starts
			  int SeqIdx,			// index into pProcParams->pSeq[] at which this subsequence starts
			  int MaxLen,			// maximum length of this subsequence
			  tsProcParams *pProcParams)
{
tsLenRangeClass *pRange;
int VIdx;
int *pStep;
etSeqBase *pVSeq;
etSeqBase VBase;
etSeqBase RefBase;
int RegionIdx;
int FeatureBits;
int BitMsk;
int NxtSeqIdx =-1;
int FeatIdx;
int TotColsMismatches = 0;
int ChromOfsEnd = ChromOfs;
int HypercoreLen = 0;			// current reference species core length (includes any refseq InDels into length)
int VMismatches;
int VIndels;
int SpliceSiteOverlaps;
int IdentityScore = cRandWalk100Score;
int CurUltraCoreLen = 0;		// current core with no mismatches
int MaxUltraCoreLen = 0;		// longest core with no mismatches encountered
int RefHyperCoreLen = 0;		// curent reference species core length (excludes any refseq InDels from length)
int WinIdx = 0;
int NumCoreSpecies;
tsDistSeg OGDistProfile[cMaxMatchDistSegments];

etSeqBase OGBase;
int OGMatches;
int OGMismatches;
int OGInDels;
int OGUnaligned;

OGMatches = 0;
OGMismatches = 0;
OGInDels = 0;
OGUnaligned = 0;

if(pProcParams->PMode == eProcModeOutspecies)
	NumCoreSpecies = pProcParams->MinAlignSpecies;
else
	NumCoreSpecies = pProcParams->NumSpeciesInAlignment;

for(HypercoreLen = 0; HypercoreLen < MaxLen ; HypercoreLen++, WinIdx++)
	{
	// determine in the current alignment column the number of mismatches and deletions (missing alignment)
	VMismatches = 0;
	VIndels = 0;

	for(VIdx = 0; VIdx < NumCoreSpecies; VIdx++)
		{
		pVSeq = pProcParams->pSeqs[VIdx];
		VBase = pVSeq[HypercoreLen+SeqIdx] & ~cRptMskFlg;
		
		if(VIdx == pProcParams->RefSpeciesIdx) // if base for reference sequence
			{
			RefBase = VBase;
			// slough columns in which the ref base is an InDel?
			if(RefBase == eBaseInDel)
				{
				if(pProcParams->bSloughRefInDels)	
					break;
				}
			if(VBase == eBaseInDel)
				VIndels++;
			if(VBase == eBaseN || VBase == eBaseInDel)	  // indeterminate ref seq bases and InDels are 
				VMismatches++;							  // treated as if mismatch
			continue;
			}

		// not reference base so can check for missmatch against relative species base
		if(VBase == eBaseInDel)
				VIndels++;
		if(VBase == eBaseInDel || VBase == eBaseN || // treat InDels and indeterminate bases as if mismatches
			RefBase != VBase)		
			VMismatches++;
		}

		// processed column of bases
		// should slough columns in which the ref base is an InDel?
	if(RefBase == eBaseInDel && pProcParams->bSloughRefInDels)	
		continue;

		// if processing outspecies mode and core column consisted of InDels then
		// incr InDel count if there is an alignment onto the out species
	if(VIndels == NumCoreSpecies && pProcParams->PMode == eProcModeOutspecies)
		{
		if(pProcParams->NumSpeciesInAlignment > NumCoreSpecies)
			OGInDels += 1;
		continue;
		}

	// number of mismatches in column now known - if more than allowed then characterise as a hyperconserved mismatch
	if(VMismatches >= pProcParams->VMismatches)
		{
		CurUltraCoreLen = 0;						// any mismatch terminates current ultracore
		if(pProcParams->MismatchScore)				// if scoring mismatches for hypers, will be 0 for ultras... 
			{
			IdentityScore -= pProcParams->MismatchScore;
			// check if the identity has dropped below minimum required
			if(IdentityScore <= 0)
				break;
			}
		if(NxtSeqIdx == -1)						// if 1st mismatch then note where next search for hypercore should start from
			NxtSeqIdx = HypercoreLen+SeqIdx+1;  // if current hypercore not accepted
				// this many total mismatches allowed?
		if(++TotColsMismatches > pProcParams->MaxHyperColsMismatches)
			break;
		}
	else	// mismatches in current column are less than max number allowed
		{
		CurUltraCoreLen++;							// current ultracore increases as bases identical
		if(CurUltraCoreLen > MaxUltraCoreLen)		// is this the longest ultracore in current subsequence?
			MaxUltraCoreLen = CurUltraCoreLen;
		if(pProcParams->MinHyperCoreLen && CurUltraCoreLen >= pProcParams->MinHyperCoreLen) // when it's an ultra of minimal length
			IdentityScore = cRandWalk100Score;			// then assume backup to 100% identity
		else										// else current ultra not at least minimal length
			{
			IdentityScore += pProcParams->MatchScore;	// increment by random walk match score
			if(IdentityScore > cRandWalk100Score)		// can't do any better than 100%...
				IdentityScore = cRandWalk100Score;
			}
		}

	if(pProcParams->PMode == eProcModeOutspecies && pProcParams->NumSpeciesInAlignment > pProcParams->MinAlignSpecies)
		{
		pVSeq = pProcParams->pSeqs[pProcParams->NumSpeciesInAlignment-1];
		OGBase = pVSeq[HypercoreLen+SeqIdx] & ~cRptMskFlg;
		if(OGBase == eBaseUndef)
			OGUnaligned += 1;
		else
			{
			if(RefBase == OGBase && RefBase != eBaseInDel)
				OGMatches+= 1;
			else
				if(RefBase == eBaseInDel || OGBase == eBaseInDel)
					OGInDels += 1;
				else
					OGMismatches += 1;
			}
		}

	if(RefBase != eBaseInDel)
		RefHyperCoreLen++;
	}

// have a sequence which we may be interested in
if(HypercoreLen==MaxLen)// if core was terminated because subsequence ended then no point in subsequent
	NxtSeqIdx = -1;		// checking for more hypercores in same subsequence as these would be internal to this core

// RefHyperCoreLen holds the reference species hypercore length excluding any InDels on the ref species sequence
if(MaxUltraCoreLen < pProcParams->MinHyperCoreLen ||
	RefHyperCoreLen < pProcParams->MinHyperCoreLen)	// if less than required then return where to start next core from
	return(NxtSeqIdx);

ChromOfsEnd = ChromOfs+RefHyperCoreLen-1;	// determine reference chromosomal offset at which sequence ends

// ensure that the hypercore is in an included region and not part of an excluded region
if(!IncludeFilter(RefChromID,ChromOfs,ChromOfsEnd,pProcParams))
		return(NxtSeqIdx);

// classify range
if((pRange = GetLengthRangeClass(RefHyperCoreLen))==nullptr)
	{
	gDiagnostics.DiagOut(eDLWarn,gszProcName,"Unexpected length range classification error - length = %d",RefHyperCoreLen);
	return(NxtSeqIdx);
	}

pStep = pProcParams->pCntStepCnts;
pStep += (pRange->ID-1) * pProcParams->Regions;
RegionIdx = 0;
FeatureBits = 0;
SpliceSiteOverlaps = 0;
if(pProcParams->pBiobed != nullptr)
	{
	if(pProcParams->BEDChromID > 0)
		FeatureBits = pProcParams->pBiobed->GetFeatureBits(pProcParams->BEDChromID,ChromOfs,ChromOfsEnd,cRegionFeatBits,pProcParams->UpDnStreamLen);
	else
		FeatureBits = 0;
	RegionIdx = 0;		// default to intergenic if no feature bits set
	if(FeatureBits)		
		{
		BitMsk = cFeatBitCDS;
		for(FeatIdx = 1; FeatIdx < 7; FeatIdx++,BitMsk <<= 1)
			{
			if(BitMsk & FeatureBits)
				{
				if(RegionIdx)		// if already have feature
					return(NxtSeqIdx);	// although was sequence of interest, more than one feature bit so can't contribute to stats
				RegionIdx = FeatIdx;
				if(pProcParams->bMultipleFeatBits)	// if multiple features allowed then don't check for any additional
					break;
				}
			}
		}

		// need to remap RegionIdx when incr counts so regions are displayed in a 5'->3' logical order
	switch(RegionIdx) {
		case 0: case 2: case 4: case 6:		// IG,5'UTR, Introns and 3'DS
			break;
		case 1:								// CDS
			RegionIdx = 3;
			break;
		case 3:								// 3'UTR
			RegionIdx = 5;
			break;
		case 5:								// 5'US
			RegionIdx = 1;
			break;
		}
	pStep[RegionIdx] += 1;	

	if(RegionIdx != 0)		// if not intergenic then check for splice sites
		{
		SpliceSiteOverlaps = pProcParams->pBiobed->GetSpliceSiteBits(pProcParams->BEDChromID,ChromOfs,ChromOfsEnd,cMinSpliceOverlap);
		if(SpliceSiteOverlaps & cIntronExonSpliceSite)
			pStep[7]++;
		if(SpliceSiteOverlaps & cExonIntronSpliceSite)
			pStep[8]++;
		}
	pProcParams->bStatsAvail = true;
	ChkOutputResults(pProcParams->szRefChrom,ChromOfs,pProcParams,false,false);	// output results as may be appropriate
	}
else
	{
	*pStep += 1;
	ChkOutputResults(pProcParams->szRefChrom,ChromOfs,pProcParams,false,false);	// output results as may be appropriate
	}

memset(OGDistProfile,0,sizeof(OGDistProfile));
if(pProcParams->PMode == eProcModeOutspecies)
	{
	if(pProcParams->NumSpeciesInAlignment == pProcParams->MinAlignSpecies)
		OGUnaligned = RefHyperCoreLen;
	else
		{
		int VMatches;
		int CoreOfs;
		
		int CoreLoci = 0;
		int CoreLenLeft = RefHyperCoreLen;
		int NumSegsLeft = pProcParams->NumDistSegs;
		int SegLoci = CoreLenLeft / NumSegsLeft--;
		int OGSegIdx = 0;
		for(CoreOfs = 0; CoreOfs < HypercoreLen ; CoreOfs++)
			{
			// determine in the current alignment column the number of matches
			VMatches = 0;
			for(VIdx = 0; VIdx < NumCoreSpecies; VIdx++)
				{
				pVSeq = pProcParams->pSeqs[VIdx];
				VBase = pVSeq[CoreOfs+SeqIdx] & ~cRptMskFlg;
				if(VIdx == pProcParams->RefSpeciesIdx) // if base for reference sequence
					{
					RefBase = VBase;
					if(VBase == eBaseN || VBase == eBaseInDel)	  // indeterminate ref seq bases and InDels are 
						break;									  // treated as if mismatch
					}
				else
					if(VBase != RefBase)
						break;
				VMatches += 1;
				}

			if(VMatches == NumCoreSpecies)				// if all core bases in column match then check on outspecies
				{
				pVSeq = pProcParams->pSeqs[pProcParams->NumSpeciesInAlignment-1];
				OGBase = pVSeq[CoreOfs+SeqIdx] & ~cRptMskFlg;
				
				if(OGBase == RefBase)					// profile outgroup bases
					OGDistProfile[OGSegIdx].Matches += 1;
				else
					{
					if(OGBase == eBaseUndef)
						OGDistProfile[OGSegIdx].Unaligned += 1;
					else
						{
						if(OGBase == eBaseInDel)
							OGDistProfile[OGSegIdx].InDels += 1;
						else
							OGDistProfile[OGSegIdx].Mismatches += 1;
						}
					}
				CoreLoci += 1;
				CoreLenLeft -= 1;
				if(NumSegsLeft && CoreLoci >= SegLoci)
					{
					OGSegIdx += 1;
					SegLoci += CoreLenLeft / NumSegsLeft--;
					}
				}
			}
		}
	}

if(!OutputHypercore(pProcParams->szRefChrom,ChromOfs,ChromOfsEnd,FeatureBits | SpliceSiteOverlaps,OGUnaligned,OGMatches,OGMismatches,OGInDels,OGDistProfile,pProcParams))
	return(-1);
if(NxtSeqIdx != -1)
	NxtSeqIdx = HypercoreLen+SeqIdx+1;
return(NxtSeqIdx);
}

// ChkOutputResults
bool
CGenHypers::ChkOutputResults(const char *pszChrom, int ChromOffset, tsProcParams *pProcParams,bool bGenEmptyRows,bool bGenEmptyWindows)
{
if(!pProcParams->WindowSize)
	return(false);
while(ChromOffset < 0 || ChromOffset > pProcParams->NxtOutputOffset)
	{
	if(bGenEmptyWindows || pProcParams->bStatsAvail)
		{
		OutputResults(pszChrom, pProcParams->NxtOutputOffset, pProcParams,bGenEmptyRows);
		if(pProcParams->bStatsAvail)
			memset(pProcParams->pCntStepCnts,0,pProcParams->NumCnts * pProcParams->Regions * sizeof(int));
		}
	pProcParams->bStatsAvail = false;
	pProcParams->NxtOutputOffset += pProcParams->WindowSize;
	if(ChromOffset < 0)
		break;
	}
return(true);
}



bool	
CGenHypers::OutputResults(const char *pszChrom, int ChromOffset, tsProcParams *pProcParams,bool bGenEmptyRows)
{
tsLenRangeClass *pRange;
static bool bOutputHdrFirst = true;
char szLineBuff[2058];
int Idx;
int Steps;
int Instances;
int Len;
int *pStep;
pStep = pProcParams->pCntStepCnts;

if(bOutputHdrFirst)
	{
	bOutputHdrFirst = false;
	if(pProcParams->Regions != 9)
		Len = sprintf(szLineBuff,"\"LenRange\",\"Mismatches\",\"TotInstances\"");
	else
		Len = sprintf(szLineBuff,"\"LenRange\",\"Mismatches\",\"TotInstances\",\"IG\",\"US\",\"5'UTR\",\"CDS\",\"INTRON\",\"3'UTR\",\"DS\",\"5'ExSplice\",\"3'ExSplice\"");
	CUtility::RetryWrites(pProcParams->hRsltsFile,szLineBuff,Len);
	}
pStep = pProcParams->pCntStepCnts;
for(Idx = 0; Idx < m_NumLenRangeBins; Idx++, pStep += pProcParams->Regions)
	{
	if(pProcParams->Regions > 1)
		{
		for(Instances = Steps = 0; Steps < (pProcParams->Regions - 2); Steps++)
			Instances += pStep[Steps];
		}
	else
		Instances = pStep[0];
	pRange = GetRangeClass(Idx+1);
	Len = sprintf(szLineBuff,"\n\"%s\",%d,%d",
						pRange->szDescr,pProcParams->MaxHyperColsMismatches,Instances);
	if(pProcParams->Regions > 1)
		for(Steps = 0; Steps < pProcParams->Regions; Steps++)
			Len += sprintf(&szLineBuff[Len],",%d",pStep[Steps]);
	if(Idx == m_NumLenRangeBins-1)
		Len += sprintf(&szLineBuff[Len],"\n");
	CUtility::RetryWrites(pProcParams->hRsltsFile,szLineBuff,Len);
	}
return(true);
}


bool	
CGenHypers::OutputHypercore(const char *pszChrom, int ChromStartOffset, int ChromEndOffset,
				int FeatureBits,		// feature bits over lapped
				int OGUnaligned,		// number of unaligned bases in outspecies
				int OGMatches,			// number of matching bases in outspecies
				int OGMismatches,		// number of mismatched bases in outspecies
				int OGInDels,			// number of InDels in outspecies
				tsDistSeg SegCnts[],	// array of segment profile counts	
				tsProcParams *pProcParams)
{
static int CoreID = 0;
int Rslt;
char szLineBuff[4096];
int Len;
if(pProcParams->hCoreCSVRsltsFile != -1)
	{
	Len = 0;
	Len += sprintf(&szLineBuff[Len],"%d,\"hypercore\",\"%s\",\"%s\",%d,%d,%d,\"%s\",%d",
			++CoreID,
			pProcParams->szSpecies[pProcParams->RefSpeciesIdx],
			pszChrom,ChromStartOffset,ChromEndOffset,ChromEndOffset-ChromStartOffset+1,
			pProcParams->pszSpeciesList,FeatureBits & (cAnyFeatBits | cOverlaysSpliceSites));
	if(pProcParams->PMode == eProcModeOutspecies)
		{
		Len += sprintf(&szLineBuff[Len],",%d,%d,%d,%d",OGUnaligned,OGMatches,OGMismatches,OGInDels);
		for(int Idx = 0; Idx < pProcParams->NumDistSegs; Idx++)
			Len += sprintf(&szLineBuff[Len],",%d,%d,%d,%d",SegCnts[Idx].Matches,SegCnts[Idx].Mismatches,SegCnts[Idx].InDels,SegCnts[Idx].Unaligned);
		}
	Len += sprintf(&szLineBuff[Len],"\n");
	if((Rslt=write(pProcParams->hCoreCSVRsltsFile,szLineBuff,Len))!=Len)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Write to loci file failed - %s",strerror(errno));
		return(false);
		}
	if(CoreID == 1 || !(CoreID % 500))
		_commit(pProcParams->hCoreCSVRsltsFile);
	}
return(true);
}

// ChkOutputSummaryResults
bool
CGenHypers::ChkOutputSummaryResults(const char *pszChrom, int ChromOffset, tsProcParams *pProcParams,bool bGenEmptyRows,bool bGenEmptyWindows)
{
if(!pProcParams->WindowSize)
	return(false);
while(ChromOffset < 0 || ChromOffset > pProcParams->NxtOutputOffset)
	{
	if(bGenEmptyWindows || pProcParams->bStatsAvail)
		{
		OutputSummaryResults(pszChrom, pProcParams->NxtOutputOffset, pProcParams,bGenEmptyRows);
		if(pProcParams->bStatsAvail)
			memset(pProcParams->pCntStepCnts,0,pProcParams->NumCnts * pProcParams->Regions * sizeof(int));
		}
	pProcParams->bStatsAvail = false;
	pProcParams->NxtOutputOffset += pProcParams->WindowSize;
	if(ChromOffset < 0)
		break;
	}
return(true);
}



bool	
CGenHypers::OutputSummaryResults(const char *pszChrom, int ChromOffset, tsProcParams *pProcParams,bool bGenEmptyRows)
{
static bool bOutputHdrFirst = true;
static char *pszClass[] = {(char *)"RefInDel",(char *)"RelIndel",(char *)"Match",(char *)"Missmatch"};
char szLineBuff[2058];
int Idx;
int Steps;
int Len;
int *pStep;
int64_t Total;
int64_t SumRefBases;
int64_t SumRelBases;

if(bOutputHdrFirst)
	{
	bOutputHdrFirst = false;
	Len = sprintf(szLineBuff,"\"Class\",\"TotInstances\"");
	if(pProcParams->Regions == 9)
		Len += sprintf(&szLineBuff[Len],",\"IG\",\"US\",\"5'UTR\",\"CDS\",\"INTRON\",\"3'UTR\",\"DS\",\"5'ExSplice\",\"3'ExSplice\"");
	CUtility::RetryWrites(pProcParams->hRsltsFile,szLineBuff,Len);
	}

pStep = pProcParams->pCntStepCnts;
for(Idx = 0; Idx < 4; Idx++,pStep += pProcParams->Regions)
	{
	for(Total = Steps = 0; Steps < pProcParams->Regions; Steps++)
		Total += pStep[Steps];
	switch(Idx) {
		case 0:
			SumRelBases = Total;
			break;
		case 1:
			SumRefBases = Total;
			break;
		default:
			SumRelBases += Total;
			SumRefBases += Total;
			break;
		}
#ifdef _WIN32
	Len = sprintf(szLineBuff,"\n\"%s\",%zd",pszClass[Idx],Total);
#else
	Len = sprintf(szLineBuff,"\n\"%s\",%ld",pszClass[Idx],Total);
#endif
	if(pProcParams->Regions > 1)
		{
		for(Steps = 0; Steps < pProcParams->Regions; Steps++)
			Len += sprintf(&szLineBuff[Len],",%d",pStep[Steps]);
		}
	CUtility::RetryWrites(pProcParams->hRsltsFile,szLineBuff,Len);
	}

pStep = pProcParams->pCntStepCnts;
for(Idx = 0; Idx < 2; Idx++)
	{
	switch(Idx) {	// reference sequence
		case 0:
#ifdef _WIN32
			Len=sprintf(szLineBuff,"\n\"RefBases\",%zd",SumRefBases);
#else
			Len=sprintf(szLineBuff,"\n\"RefBases\",%ld",SumRefBases);
#endif
			break;
		case 1:		// relative sequence
#ifdef _WIN32
			Len=sprintf(szLineBuff,"\n\"RelBases\",%zd",SumRelBases);
#else
			Len=sprintf(szLineBuff,"\n\"RelBases\",%ld",SumRelBases);
#endif
			break;
		}

	if(pProcParams->Regions > 1)
		{
		for(Steps = 0; Steps < pProcParams->Regions; Steps++)
			{
			switch(Idx) {
				case 0:	// reference sequence
					SumRefBases = pStep[pProcParams->Regions + Steps] + 
						pStep[(pProcParams->Regions * 2) + Steps] + 
						pStep[(pProcParams->Regions * 3) + Steps];
					break;
				case 1:	// relative sequence
					SumRefBases = pStep[Steps] + 
						pStep[(pProcParams->Regions * 2) + Steps] + 
						pStep[(pProcParams->Regions * 3) + Steps];
					break;
				}
#ifdef _WIN32
			Len += sprintf(&szLineBuff[Len],",%zd",SumRefBases);
#else
			Len += sprintf(&szLineBuff[Len],",%ld",SumRefBases);
#endif
			}
		}
	CUtility::RetryWrites(pProcParams->hRsltsFile,szLineBuff,Len);
	}
return(true);
}



