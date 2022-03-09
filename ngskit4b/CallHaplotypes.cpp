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




int CallHaplotypes(eModeCSH PMode,	// processing mode 0: report imputation haplotype matrix, 1: report both raw and imputation haplotype matrices, 2: additionally generate GWAS allowing visual comparisons, 3: allelic haplotype grouping,4: coverage haplotype grouping, 5: post-process to WIG (default 0)
	        int32_t LimitPrimaryFiles,                // limit number of loaded primary or founder PBA files to this many. 0: no limits, > 0 sets upper limit
			int32_t ExprID,			// assign this experiment identifier for this PBA analysis
			int32_t GrpHapBinSize,           // when grouping into haplotypes (processing mode 3) then use this non-overlapping bin size for accumulation of differentials between all founder PBAs
	        int32_t NumHapGrpPhases,         // number of phases over which to converge group consensus when haplotype clustering into groups - moves the balance between optimal group consensus and processing resource
			int32_t MinCentClustDist,        // haplotype groupings - in processing mode 3/4 only - minimum group centroid clustering distance
			int32_t MaxCentClustDist,        // haplotype groupings - in processing mode 3/4 only - maximum group centroid clustering distance
			int32_t MaxClustGrps,       // haplotype groupings - in processing mode 3 only - targeted maximum number of groups
			int32_t FndrTrim5,			// trim this many aligned PBAs from 5' end of founder aligned segments - reduces false alleles due to sequencing errors
			int32_t FndrTrim3,			// trim this many aligned PBAs from 3' end of founder aligned segments - reduces false alleles due to sequencing errors
			int32_t ProgTrim5,			// trim this many aligned PBAs from 5' end of progeny aligned segments - reduces false alleles due to sequencing errors
			int32_t ProgTrim3,			// trim this many aligned PBAs from 3' end of progeny aligned segments - reduces false alleles due to sequencing errors
			int32_t WWRLProxWindow,		// proximal window size for Wald-Wolfowitz runs test
			int32_t OutliersProxWindow,	// proximal window size for outliers reduction
			char *pszMaskBPAFile,	// optional masking input BPA file, only process BPAs which are intersect of these BPAs and progeny plus founder BPAs
			int NumFounderInputFiles,	// number of input founder file specs
			char *pszFounderInputFiles[],	// names of input founder PBA files (wildcards allowed)
			int NumProgenyInputFiles,	// number of input progeny file specs
			char* pszProgenyInputFiles[],		// names of input progeny PBA files (wildcards allowed)
			char* pszOutFile,		// loci haplotype calls output file (CSV format)
			int	NumIncludeChroms,			// number of chromosome regular expressions to include
			char **ppszIncludeChroms,		// array of include chromosome regular expressions
			int	NumExcludeChroms,			// number of chromosome expressions to exclude
			char **ppszExcludeChroms,		// array of exclude chromosome regular expressions
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

eModeCSH PMode;				// processing mode 0: report imputation haplotype matrix, 1: report both raw and imputation haplotype matrices, 2: additionally generate GWAS allowing visual comparisons,3: allelic haplotype grouping,4: coverage haplotype grouping, 5: post-process to WIG (default 0) 
int32_t LimitPrimaryPBAs;                // limit number of loaded primary or founder PBA files to this many. 0: no limits, > 0 sets upper limit
int32_t ExprID;			    // assign this experiment identifier to this PBA analysis
int32_t FndrTrim5;		    // trim this many aligned PBAs from 5' end of founder aligned segments - reduces false alleles due to sequencing errors
int32_t FndrTrim3;		    // trim this many aligned PBAs from 3' end of founder aligned segments - reduces false alleles due to sequencing errors
int32_t ProgTrim5;		    // trim this many aligned PBAs from 5' end of progeny aligned segments - reduces false alleles due to sequencing errors
int32_t ProgTrim3;		    // trim this many aligned PBAs from 3' end of progeny aligned segments - reduces false alleles due to sequencing errors
int32_t WWRLProxWindow;		// proximal window size for Wald-Wolfowitz runs test
int32_t OutliersProxWindow;	// proximal window size for outliers reduction
int32_t GrpHapBinSize;           // when grouping into haplotypes (processing mode 3) then use this non-overlapping bin size for accumulation of differentials between all founder PBAs
int32_t NumHapGrpPhases;         // number of phases over which to converge group consensus when haplotype clustering into groups - is a balance between optimal group consensus and processing resource
int32_t MinCentClustDist;        // haplotype groupings - in processing mode 3 only - minimum group centroid clustering distance (default 10, min 1, clamped to bin size)
int32_t MaxCentClustDist;        // haplotype groupings - in processing mode 3 only - maximum group centroid clustering distance (default mincentclustdist, clamped to bin size)
int32_t MaxClustGrps;       // haplotype groupings - in processing mode 3 only - targeted maximum number of groups (default 10)

int NumIncludeChroms;
char *pszIncludeChroms[cMaxIncludeChroms];
int NumExcludeChroms;
char *pszExcludeChroms[cMaxExcludeChroms];

int NumFounderInputFiles;		// number of input founder BPA files
char *pszFounderInputFiles[cMaxFounderFileSpecs];		// names of input founder BPA files (wildcards allowed)
int NumProgenyInputFiles;		// number of input progeny BPA files
char *pszProgenyInputFiles[cMaxProgenyFileSpecs];		// names of input progeny BPA files (wildcards allowed)

char szOutFile[_MAX_PATH];		// Windowed haplotype calls output file base name, progeny readset identifier is appended to this base name
char szMaskBPAFile[_MAX_PATH];	// optional masking input BPA file, only process BPAs which are intersect of these BPAs and progeny plus founder BPAs

struct arg_lit *help = arg_lit0 ("h", "help", "print this help and exit");
struct arg_lit *version = arg_lit0 ("v", "version,ver", "print version information and exit");
struct arg_int *FileLogLevel = arg_int0 ("f", "FileLogLevel", "<int>", "Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0 ("F", "log", "<file>", "diagnostics log file");

struct arg_int *pmode = arg_int0 ("m", "mode", "<int>", "processing mode 0: report imputation haplotype matrix, 1: report both raw and imputation haplotype matrices, 2: additionally generate GWAS allowing visual comparisons, 3: allelic haplotype grouping,4: coverage haplotype grouping, 5: post-process to WIG (default 0)");
struct arg_int *limitprimarypbas = arg_int0 ("l", "debug", "<int>", " limit number of loaded primary or founder PBA files to this many. 0: no limits, > 0 sets upper limit (default 0)");

struct arg_int* exprid = arg_int0("e","exprid","<int>","assign this experiment identifier for haplotypes called (default 1");
struct arg_int* fndrtrim5 = arg_int0("y", "fndrtrim5", "<int>", "trim this many aligned PBAs from 5' end of founder aligned segments (default 0)");
struct arg_int* fndrtrim3 = arg_int0("Y", "fndrtrim3", "<int>", "trim this many aligned PBAs from 3' end of founder aligned segments (default 0)");
struct arg_int* progtrim5 = arg_int0("x", "progtrim5", "<int>", "trim this many aligned PBAs from 5' end of progeny aligned segments (default 0)");
struct arg_int* progtrim3 = arg_int0("X", "progtrim3", "<int>", "trim this many aligned PBAs from 3' end of progeny aligned segments (default 0)");

struct arg_int* wwrlproxwindow = arg_int0("w", "wwrlproxwindow", "<int>", "proximal window size for Wald-Wolfowitz runs test (default 1000000)");
struct arg_int* outliersproxwindow = arg_int0("W", "outliersproxwindow", "<int>", "proximal window size for outliers reduction (default 1000000)");

struct arg_int* grphapbinsize = arg_int0("g", "grphapbinsize", "<int>", "haplotype groupings - in processing mode 3/4 only - bin size (default 10000, min 10, will be clamped to actual size)");
struct arg_int* maxclustgrps = arg_int0("G", "maxclustgrps", "<int>", "haplotype groupings - in processing mode 3/4 only - targeted maximum number of groups (default 10)");
struct arg_int* gpphases = arg_int0("p", "gpphases", "<int>", "haplotype groupings - in processing mode 3/4 only - max number of cluster processing phases (default 10)");

struct arg_int* mincentclustdist = arg_int0("c", "mincentclustdist", "<int>", "haplotype groupings - in processing mode 3/4 only - minimum group centroid clustering distance (default 10, min 1, clamped to bin size)");
struct arg_int* maxcentclustdist = arg_int0("C", "maxcentclustdist", "<int>", "haplotype groupings - in processing mode 3/4 only - maximum group centroid clustering distance (default 10000, clamped to bin size)");

struct arg_str  *ExcludeChroms = arg_strn("Z","chromexclude",	"<string>",0,cMaxExcludeChroms,"high priority - regular expressions defining chromosomes to exclude");
struct arg_str  *IncludeChroms = arg_strn("z","chromeinclude",	"<string>",0,cMaxIncludeChroms,"low priority - regular expressions defining chromosomes to include");

	
struct arg_file *founderfiles = arg_filen("I", "founderfiles", "<file>", 0,cMaxFounderFileSpecs,"founder input BPA file(s), wildcards allowed, limit of 200 founder filespecs supported");
struct arg_file *progenyfiles = arg_filen("i", "inprogenyfile", "<file>",0, cMaxProgenyFileSpecs, "progeny input BPA file(s), wildcards allowed, limit of 500 progeny filespecs supported");
struct arg_file *maskbpafile = arg_file0("j", "inmaskbpa", "<file>", "optional masking input BPA file, or haplotype grouping specification file(s) or previously generated haplotype grouping file for post-processing");
struct arg_file *outfile = arg_file1("o", "out", "<file>", "loci haplotype calls output prefix (outputs CSV and WIG format)");
struct arg_int *threads = arg_int0("T","threads","<int>","number of processing threads 0..64 (defaults to 0 which limits threads to maximum of 64 CPU cores)");
struct arg_end *end = arg_end (200);

void *argtable[] = { help,version,FileLogLevel,LogFile,
					pmode,limitprimarypbas,grphapbinsize,gpphases,mincentclustdist,maxcentclustdist,maxclustgrps,exprid,fndrtrim5,fndrtrim3,progtrim5,progtrim3,wwrlproxwindow,outliersproxwindow,
					maskbpafile,IncludeChroms,ExcludeChroms,progenyfiles,founderfiles, outfile,threads,end };

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

		PMode = pmode->count ? (eModeCSH)pmode->ival[0] : eMCSHDefault;
		if(PMode < eMCSHDefault || PMode >= eMCSHPlaceHolder)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Unsupported processing mode '-m%d'\n",PMode);
			exit(1);
			}
		if(PMode <= eMCSHCoverageHapsGrps)
			{
			LimitPrimaryPBAs = limitprimarypbas->count ? limitprimarypbas->ival[0] : 0;
			if(LimitPrimaryPBAs < 0)
				LimitPrimaryPBAs = 0;
			}
		else
			LimitPrimaryPBAs = 0;

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

		NumHapGrpPhases = 0;
		if(PMode == eMCSHAllelicHapsGrps || PMode == eMCSHCoverageHapsGrps) 
			{
			NumHapGrpPhases = gpphases->count ? gpphases->ival[0] : cDfltNumHapGrpPhases;
			if(NumHapGrpPhases < 2) // silently clamping ...
				NumHapGrpPhases = 2;
			else
			if(NumHapGrpPhases > 100)
				NumHapGrpPhases = 100;
			}

		if(IncludeChroms->count)
			{
			for(NumIncludeChroms = Idx = 0; NumIncludeChroms < cMaxIncludeChroms && Idx < IncludeChroms->count; Idx++)
				{
				pszIncludeChroms[Idx] = NULL;
				if(pszIncludeChroms[NumIncludeChroms] == NULL)
					pszIncludeChroms[NumIncludeChroms] = new char[_MAX_PATH];
				strncpy(pszIncludeChroms[NumIncludeChroms], IncludeChroms->sval[Idx], _MAX_PATH);
				pszIncludeChroms[NumIncludeChroms][_MAX_PATH - 1] = '\0';
				CUtility::TrimQuotedWhitespcExtd(pszIncludeChroms[NumIncludeChroms]);
				if(pszIncludeChroms[NumIncludeChroms][0] != '\0')
					NumIncludeChroms++;
				}
			}
		else
			NumIncludeChroms = 0;
		

		if(ExcludeChroms->count)
			{
			for(NumExcludeChroms = Idx = 0; NumExcludeChroms < cMaxExcludeChroms && Idx < ExcludeChroms->count; Idx++)
				{
				pszExcludeChroms[Idx] = NULL;
				if(pszExcludeChroms[NumExcludeChroms] == NULL)
					pszExcludeChroms[NumExcludeChroms] = new char[_MAX_PATH];
				strncpy(pszExcludeChroms[NumExcludeChroms], ExcludeChroms->sval[Idx], _MAX_PATH);
				pszExcludeChroms[NumExcludeChroms][_MAX_PATH - 1] = '\0';
				CUtility::TrimQuotedWhitespcExtd(pszExcludeChroms[NumExcludeChroms]);
				if(pszExcludeChroms[NumExcludeChroms][0] != '\0')
					NumExcludeChroms++;
				}
			}
		else
			NumExcludeChroms = 0;

		if((PMode == eMCSHAllelicHapsGrps || PMode == eMCSHCoverageHapsGrps) && maskbpafile->count == 0)
			{
			GrpHapBinSize = grphapbinsize->count ? grphapbinsize->ival[0] : cDfltGrpHapBinSize;
			if(GrpHapBinSize < 10)        // silently clamping to a more reasonable size!
				GrpHapBinSize = 10;       // upper limit will be actual chromosome size

            MinCentClustDist = mincentclustdist->count ? mincentclustdist->ival[0] : cDfltMinCentClustDist;
			if(MinCentClustDist < 1) 
				MinCentClustDist = 1;
			else
				if(MinCentClustDist > GrpHapBinSize) // silently clamp grouping size to window size-1!
					MinCentClustDist = GrpHapBinSize-1;
			MaxCentClustDist = maxcentclustdist->count ? maxcentclustdist->ival[0] : cDfltMaxCentClustDist;
			if(MaxCentClustDist < MinCentClustDist) 
				MaxCentClustDist = MinCentClustDist;
			else
				if(MaxCentClustDist >= GrpHapBinSize) // silently clamp grouping size to window size-1!
					MaxCentClustDist = GrpHapBinSize-1;
			MaxClustGrps = maxclustgrps->count ? maxclustgrps->ival[0] : cDfltMaxClustGrps;
			if(MaxClustGrps < 1) 
				MaxClustGrps = 1;
			else
				if(MaxClustGrps > cMaxClustGrps) // silently clamp
					MaxClustGrps = cMaxClustGrps;
			}
		else
			{
			GrpHapBinSize = 0;
			MinCentClustDist = 0;
			MaxCentClustDist = 0;
			MaxClustGrps = 0;
			}

		if(PMode < eMCSHAllelicHapsGrps)
			{
			FndrTrim5 = fndrtrim5->count ? fndrtrim5->ival[0] : cDfltFndrTrim5;
			if(FndrTrim5 < 0)
				FndrTrim5 = 0;
			else
				if(FndrTrim5 > 1000)
					FndrTrim5 = 1000;

			FndrTrim3 = fndrtrim3->count ? fndrtrim3->ival[0] : cDfltFndrTrim3;
			if(FndrTrim3 < 0)
				FndrTrim3 = 0;
			else
				if(FndrTrim3 > 1000)
					FndrTrim3 = 1000;

			ProgTrim5 = progtrim5->count ? progtrim5->ival[0] : cDfltProgTrim5;
			if(ProgTrim5 < 0)
				ProgTrim5 = 0;
			else
				if(ProgTrim5 > 1000)
					ProgTrim5 = 1000;

			ProgTrim3 = progtrim3->count ? progtrim3->ival[0] : cDfltProgTrim3;
			if(ProgTrim3 < 0)
				ProgTrim3 = 0;
			else
				if(ProgTrim3 > 1000)
					ProgTrim3 = 1000;

			if(!(PMode == eMCSHAllelicHapsGrps || PMode == eMCSHCoverageHapsGrps))
				{
				WWRLProxWindow = wwrlproxwindow->count ? wwrlproxwindow->ival[0] : cDfltWWRLProxWindow;
				if(WWRLProxWindow < 10000)
					WWRLProxWindow = 0;
				else
					if(WWRLProxWindow > 10000000)	// silently clamping to 10Mbp
						WWRLProxWindow = 10000000;

				OutliersProxWindow = outliersproxwindow->count ? outliersproxwindow->ival[0] : cDlftOutliersProxWindow;
				if(OutliersProxWindow < 10000)
					OutliersProxWindow = 0;
				else
					if(OutliersProxWindow > 10000000)	// silently clamping to 10Mbp
						OutliersProxWindow = 10000000;
				}
			else
				{
				WWRLProxWindow = 0;
				OutliersProxWindow = 0;
				}
			}
		else
			{
			FndrTrim5 = 0;
			FndrTrim3 = 0;
			ProgTrim5 = 0;
			ProgTrim3 = 0;
			WWRLProxWindow = 0;
			OutliersProxWindow = 0;
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

	int MaxAllowedThreads = min(cMaxPBAWorkerThreads,NumberOfProcessors);	// limit to be at most cMaxWorkerThreads
	if((NumThreads = threads->count ? threads->ival[0] : MaxAllowedThreads)==0)
		NumThreads = MaxAllowedThreads;
	if(NumThreads < 0 || NumThreads > MaxAllowedThreads)
		{
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"Warning: Number of threads '-T%d' specified was outside of range %d..%d",NumThreads,1,MaxAllowedThreads);
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"Warning: Defaulting number of threads to %d",MaxAllowedThreads);
		NumThreads = MaxAllowedThreads;
		}

	NumProgenyInputFiles = 0;

	if(PMode < eMCSHAllelicHapsGrps)
		{
		if(progenyfiles->count)
			{
			for(Idx = 0; NumProgenyInputFiles < cMaxProgenyReadsets && Idx < progenyfiles->count; Idx++)
				{
				pszProgenyInputFiles[Idx] = NULL;
				if(pszProgenyInputFiles[NumProgenyInputFiles] == NULL)
					pszProgenyInputFiles[NumProgenyInputFiles] = new char[_MAX_PATH];
				strncpy(pszProgenyInputFiles[NumProgenyInputFiles], progenyfiles->filename[Idx], _MAX_PATH);
				pszProgenyInputFiles[NumProgenyInputFiles][_MAX_PATH - 1] = '\0';
				CUtility::TrimQuotedWhitespcExtd(pszProgenyInputFiles[NumProgenyInputFiles]);
				if(pszProgenyInputFiles[NumProgenyInputFiles][0] != '\0')
					NumProgenyInputFiles++;
				}
			}


		if(!NumProgenyInputFiles)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: After removal of whitespace, there are no input progeny file(s) specified with '-i<filespec>' option)\n");
			exit(1);
			}
		}

	NumFounderInputFiles = 0;
	if(PMode < eMCSHGrpDist2WIG)
		{
		if(founderfiles->count)
			{
			for(NumFounderInputFiles = Idx = 0; NumFounderInputFiles < cMaxFounderFileSpecs && Idx < founderfiles->count; Idx++)
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
		}

	if(maskbpafile->count)
		{
		strcpy (szMaskBPAFile, maskbpafile->filename[0]);
		CUtility::TrimQuotedWhitespcExtd (szMaskBPAFile);
		if(szMaskBPAFile[0] == '\0')
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: After removal of whitespace, no input %s file specified with '-z<filespec>' option)\n",(PMode == eMCSHAllelicHapsGrps || PMode == eMCSHCoverageHapsGrps) ? "haplotype group clustering specification" : "masking BPA");
			exit(1);
			}
		}
	else
		szMaskBPAFile[0] = '\0';

	if(PMode >= eMCSHGrpDist2WIG)
		{
		if(szMaskBPAFile[0] == '\0')
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: No previously generated haplotype grouping file specified with '-z<filespec>' option, required for post-processing centroid distances into WIG output format)\n");
			exit(1);
			}
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
		case eMCSHDefault:
			pszDescr = "Report imputated haplotype matrix";
			break;
		case eMCSHRaw:
			pszDescr = "Report both raw and imputated haplotype matrices";
			break;
		case eMCSHGWAS:
			pszDescr = "Report both raw and imputation haplotype matrices plus generate pre/post imputation GWAS formated files";
			break;
		case eMCSHAllelicHapsGrps:
			pszDescr = "Report allelic haplotype groupings using centroid clustering";
			break;
		case eMCSHCoverageHapsGrps:
			pszDescr = "Report coverage haplotype groupings using centroid clustering";
			break;
		case eMCSHGrpDist2WIG:
			pszDescr = "Post-processing haplotype grouping centroid distances into WIG format";
			break;
		}

	gDiagnostics.DiagOutMsgOnly (eDLInfo, "Calling haplotypes : '%s'", pszDescr);
	
	for(Idx = 0; Idx < NumIncludeChroms; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"reg expressions defining chroms to include: '%s'",pszIncludeChroms[Idx]);
	for(Idx = 0; Idx < NumExcludeChroms; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"reg expressions defining chroms to exclude: '%s'",pszExcludeChroms[Idx]);

	if(LimitPrimaryPBAs > 0)
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "Limiting number of loaded primary or founder PBA files to first: %d",LimitPrimaryPBAs);

	gDiagnostics.DiagOutMsgOnly(eDLInfo, "Experiment identifier: %d", ExprID);
	if(PMode < eMCSHGrpDist2WIG)
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "Trim this many aligned PBAs from 5' end of founder aligned segments: %d", FndrTrim5);
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "Trim this many aligned PBAs from 3' end of founder aligned segments: %d", FndrTrim3);
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "Trim this many aligned PBAs from 5' end of progeny aligned segments: %d", ProgTrim5);
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "Trim this many aligned PBAs from 3' end of progeny aligned segments: %d", ProgTrim3);
		}


	if(PMode < eMCSHAllelicHapsGrps)
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "Proximal window size for Wald-Wolfowitz runs test: %d", WWRLProxWindow);
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "Proximal window size for outliers reduction: %d", OutliersProxWindow);
		if(szMaskBPAFile[0] != '\0')
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "Loading asking BPAs from : '%s'", szMaskBPAFile);
		}
	else
		{
		if(PMode == eMCSHAllelicHapsGrps || PMode == eMCSHCoverageHapsGrps)
			{
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "Number of haplotype grouping centroid clustering phases: %d", NumHapGrpPhases);
			if(szMaskBPAFile[0] != '\0')
				{
				gDiagnostics.DiagOutMsgOnly(eDLInfo, "Loading group centroid clustering specifications from : '%s'", szMaskBPAFile);
				GrpHapBinSize = MinCentClustDist = MaxCentClustDist = MaxClustGrps = 0;
				}
			else
				{
				gDiagnostics.DiagOutMsgOnly(eDLInfo, "Haplotype groupings window size: %d", GrpHapBinSize);
				gDiagnostics.DiagOutMsgOnly(eDLInfo, "Minimum group centroid clustering differential: %d", MinCentClustDist);
				gDiagnostics.DiagOutMsgOnly(eDLInfo, "Maximum group centroid clustering differential: %d", MaxCentClustDist);
				gDiagnostics.DiagOutMsgOnly(eDLInfo, "Targeted maximum number of haplotype groups: %d", MaxClustGrps);
				}
			}
		else
			{
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "Loading previously generated haplotype grouping centroid distances from : '%s'", szMaskBPAFile);
			}
		}

	if(PMode < eMCSHGrpDist2WIG)
		{
		for(Idx = 0; Idx < NumFounderInputFiles; Idx++)
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "Sample or Founder file : '%s'", pszFounderInputFiles[Idx]);

		if(NumProgenyInputFiles)
			for(Idx = 0; Idx < NumProgenyInputFiles; Idx++)
				gDiagnostics.DiagOutMsgOnly(eDLInfo, "Progeny file : '%s'", pszProgenyInputFiles[Idx]);
		}

	gDiagnostics.DiagOutMsgOnly (eDLInfo, "Output file : '%s'", szOutFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"number of threads : %d",NumThreads);

#ifdef _WIN32
	SetPriorityClass (GetCurrentProcess (), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start ();
	Rslt = 0;
	Rslt = CallHaplotypes(PMode,			// processing mode 0: report imputation haplotype matrix, 1: report both raw and imputation haplotype matrices, 2: additionally generate GWAS allowing visual comparisons,3: allelic haplotype grouping,4: coverage haplotype grouping, 5: post-process to WIG (default 0)
		            LimitPrimaryPBAs,                //limit number of loaded primary or founder PBA files to this many. 0: no limits, > 0 sets upper limit			
		            ExprID,					// experiment identifier
					GrpHapBinSize,           // when grouping into haplotypes (processing mode 3) then use this non-overlapping bin size for accumulation of differentials between all founder PBAs
		            NumHapGrpPhases,         // number of phases over which to converge group consensus when haplotype clustering into groups - is a balance between optimal group consensus and processing resource
					MinCentClustDist,        // haplotype groupings - in processing mode 3 only - minimum group centroid clustering distance (default 500, min 1, clamped to bin size)
					MaxCentClustDist,        // haplotype groupings - in processing mode 3 only - maximum group centroid clustering distance (default mincentclustdist, clamped to bin size)
					MaxClustGrps,       // haplotype groupings - in processing mode 3 only - targeted maximum number of groups (default 50)
					FndrTrim5,			// trim this many aligned PBAs from 5' end of founder aligned segments - reduces false alleles due to sequencing errors
					FndrTrim3,			// trim this many aligned PBAs from 3' end of founder aligned segments - reduces false alleles due to sequencing errors
					ProgTrim5,			// trim this many aligned PBAs from 5' end of progeny aligned segments - reduces false alleles due to sequencing errors
					ProgTrim3,			// trim this many aligned PBAs from 3' end of progeny aligned segments - reduces false alleles due to sequencing errors
					WWRLProxWindow,		// proximal window size for Wald-Wolfowitz runs test
					OutliersProxWindow,	// proximal window size for outliers reduction
					szMaskBPAFile,			// optional input masking BPA file or Haplotype group clustering specifications file(s)
					NumFounderInputFiles,	// number of input founder file specs
					pszFounderInputFiles,	// names of input founder PBA files (wildcards allowed)
					NumProgenyInputFiles,		// number of input progeny file specs
					pszProgenyInputFiles,		// names of input progeny PBA files (wildcards allowed)
					szOutFile,				// output to this file
					NumIncludeChroms,			// number of chromosome regular expressions to include
					pszIncludeChroms,		// array of include chromosome regular expressions
					NumExcludeChroms,			// number of chromosome expressions to exclude
					pszExcludeChroms,		// array of exclude chromosome regular expressions
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

int CallHaplotypes(eModeCSH PMode,	// processing mode 0: report imputation haplotype matrix, 1: report both raw and imputation haplotype matrices, 2: additionally generate GWAS allowing visual comparisons,3: allelic haplotype grouping,4: coverage haplotype grouping, 5: post-process to WIG
				   int32_t LimitPrimaryPBAs,                // limit number of loaded primary or founder PBA files to this many. 0: no limits, > 0 sets upper limit
	               int32_t ExprID,			// assign this experiment identifier for this PBA analysis
					int32_t GrpHapBinSize,           // when grouping into haplotypes (processing mode 3) then use this non-overlapping bin size for accumulation of differentials between all founder PBAs
	                int32_t NumHapGrpPhases,         // number of phases over which to converge group consensus when haplotype clustering into groups - is a balance between optimal group consensus and processing resource
					int32_t MinCentClustDist,        // haplotype groupings - in processing mode 3 only - minimum group centroid clustering distance (default 500, min 1, clamped to bin size)
					int32_t MaxCentClustDist,        // haplotype groupings - in processing mode 3 only - maximum group centroid clustering distance (default mincentclustdist, clamped to bin size)
					int32_t MaxClustGrps,       // haplotype groupings - in processing mode 3 only - targeted maximum number of groups (default 50)
	               int32_t FndrTrim5,			// trim this many aligned PBAs from 5' end of founder aligned segments - reduces false alleles due to sequencing errors
				   int32_t FndrTrim3,			// trim this many aligned PBAs from 3' end of founder aligned segments - reduces false alleles due to sequencing errors
				   int32_t ProgTrim5,			// trim this many aligned PBAs from 5' end of progeny aligned segments - reduces false alleles due to sequencing errors
				   int32_t ProgTrim3,			// trim this many aligned PBAs from 3' end of progeny aligned segments - reduces false alleles due to sequencing errors
				   int32_t WWRLProxWindow,		// proximal window size for Wald-Wolfowitz runs test
				   int32_t OutliersProxWindow,	// proximal window size for outliers reduction
				   char* pszMaskBPAFile,	 // optional input masking BPA file or Haplotype group clustering specifications file(s)
				   int NumFounderInputFiles,	// number of input founder file specs
				   char* pszFounderInputFiles[],	// names of input founder PBA files (wildcards allowed)
				   int NumProgenyInputFiles,	// number of input progeny file specs
				   char* pszProgenyInputFiles[],		// names of input progeny PBA files (wildcards allowed)
				   char* pszOutFile,		// Windowed haplotype calls output file (CSV format)
					int	NumIncludeChroms,			// number of chromosome regular expressions to include
					char **ppszIncludeChroms,		// array of include chromosome regular expressions
					int	NumExcludeChroms,			// number of chromosome expressions to exclude
					char **ppszExcludeChroms,		// array of exclude chromosome regular expressions
				   int NumThreads)		// number of worker threads to use
{
int Rslt;
CCallHaplotypes* pCallHaplotypes;

if((pCallHaplotypes = new CCallHaplotypes) == NULL)
	{
	gDiagnostics.DiagOut (eDLFatal, gszProcName, "Unable to instantiate instance of CCallHaplotypes");
	return(eBSFerrObj);
	}
Rslt = pCallHaplotypes->Process(PMode,LimitPrimaryPBAs,ExprID,GrpHapBinSize,NumHapGrpPhases,MinCentClustDist,MaxCentClustDist,MaxClustGrps,FndrTrim5, FndrTrim3, ProgTrim5, ProgTrim3, WWRLProxWindow,OutliersProxWindow,
				pszMaskBPAFile,NumFounderInputFiles,
				pszFounderInputFiles,NumProgenyInputFiles,pszProgenyInputFiles,pszOutFile,NumIncludeChroms,ppszIncludeChroms,NumExcludeChroms,ppszExcludeChroms,NumThreads);
delete pCallHaplotypes;
return(Rslt);
}


CCallHaplotypes::CCallHaplotypes()
{
m_pChromMetadata = NULL;
m_pProgenyFndrAligns = NULL;
m_pWorkQueueEls = NULL;
m_pInBuffer = NULL;
m_pszOutBuffer = NULL;
m_pszWIGBuff = NULL;
m_pAlleleStacks = NULL;
m_pHGBinSpecs = NULL;
m_hOutFile = -1;
m_hWIGOutFile = -1;
m_hInFile = -1;
m_pCSVFile = NULL;
m_pHGCSVFile = NULL;
m_bMutexesCreated = false;
Reset();
}

CCallHaplotypes::~CCallHaplotypes()
{
if(m_hInFile != -1)
	close(m_hInFile);
if(m_hOutFile != -1)
	close(m_hOutFile);
if(m_hWIGOutFile != -1)
    close(m_hWIGOutFile);
if(m_pWorkQueueEls != NULL)
	delete []m_pWorkQueueEls;
if(m_pszOutBuffer != NULL)
	delete []m_pszOutBuffer;
if(m_pszWIGBuff != NULL)
	delete []m_pszWIGBuff;
if(m_pInBuffer != NULL)
	delete []m_pInBuffer;

if(m_pCSVFile != NULL)
	delete m_pCSVFile;
if(m_pHGCSVFile != NULL)
    delete m_pHGCSVFile;

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

if(m_pProgenyFndrAligns != NULL)
	{
#ifdef _WIN32
	free(m_pProgenyFndrAligns);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pProgenyFndrAligns != MAP_FAILED)
		munmap(m_pProgenyFndrAligns, m_AllocdProgenyFndrAlignsMem);
#endif
	}

if(m_pAlleleStacks != NULL)
{
#ifdef _WIN32
	free(m_pAlleleStacks);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pAlleleStacks != MAP_FAILED)
		munmap(m_pAlleleStacks, m_AllocdAlleleStacksMem);
#endif
	}

if(m_pHGBinSpecs != NULL)
{
	tsHGBinSpec* pBin;
	pBin = m_pHGBinSpecs;
	for(int BinIdx = 0; BinIdx < m_UsedHGBinSpecs; BinIdx++, pBin++)
		if(pBin->pHaplotypeGroup != NULL)
			delete[]pBin->pHaplotypeGroup;
#ifdef _WIN32
	free(m_pHGBinSpecs);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pHGBinSpecs != MAP_FAILED)
		munmap(m_pHGBinSpecs, m_AllocdHGBinSpecsMem);
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
if(m_hWIGOutFile != -1)
	{
	close(m_hWIGOutFile);
	m_hWIGOutFile = -1;
	}

if(m_pCSVFile != NULL)
	{
	delete m_pCSVFile;
	m_pCSVFile = NULL;
	}

if(m_pHGCSVFile != NULL)
	{
	delete m_pHGCSVFile;
	m_pHGCSVFile = NULL;
	}

if(m_pWorkQueueEls != NULL)
	{
	delete []m_pWorkQueueEls;
	m_pWorkQueueEls = NULL;
	}
m_AllocWorkQueueEls = 0;
m_TotWorkQueueEls = 0;
m_NumQueueElsProcessed = 0;
m_FastSerialise = 0;

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

if(m_pszWIGBuff != NULL)
	{
	delete []m_pszWIGBuff;
	m_pszWIGBuff = NULL;
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

if(m_pAlleleStacks != NULL)
	{
#ifdef _WIN32
	free(m_pAlleleStacks);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pAlleleStacks != MAP_FAILED)
		munmap(m_pAlleleStacks, m_AllocdAlleleStacksMem);
#endif
	m_pAlleleStacks = NULL;
	}
m_AllocdAlleleStacksMem = 0;
m_UsedAlleleStacks = 0;
m_AllocdAlleleStacks = 0;

if(m_pHGBinSpecs != NULL)
	{
	tsHGBinSpec* pBin;
	pBin = m_pHGBinSpecs;
	for(int BinIdx = 0; BinIdx < m_UsedHGBinSpecs; BinIdx++, pBin++)
		if(pBin->pHaplotypeGroup != NULL)
			delete[]pBin->pHaplotypeGroup;
#ifdef _WIN32
	free(m_pHGBinSpecs);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pHGBinSpecs != MAP_FAILED)
		munmap(m_pHGBinSpecs, m_AllocdHGBinSpecsMem);
#endif
	m_pHGBinSpecs = NULL;
	}
m_UsedHGBinSpecs = 0;
m_AllocdHGBinSpecs = 0;
m_AllocdHGBinSpecsMem = 0;

DeleteMutexes();
if(m_pProgenyFndrAligns != NULL)
	{
#ifdef _WIN32
	free(m_pProgenyFndrAligns);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pProgenyFndrAligns != MAP_FAILED)
		munmap(m_pProgenyFndrAligns, m_AllocdProgenyFndrAlignsMem);
#endif
	m_pProgenyFndrAligns = NULL;
	}
m_UsedProgenyFndrAligns = 0;
m_AllocdProgenyFndrAligns = 0;
m_AllocdProgenyFndrAlignsMem = 0;

m_NumFounders = 0;
m_NumProgenies = 0;

m_ExprID = 0;
m_FndrTrim5 = cDfltFndrTrim5;
m_FndrTrim3 = cDfltFndrTrim3;
m_ProgTrim5 = cDfltProgTrim5;
m_ProgTrim3 = cDfltProgTrim3;
m_GrpHapBinSize = cDfltGrpHapBinSize;
m_NumHapGrpPhases = cDfltNumHapGrpPhases;
m_MinCentClustDist = cDfltMinCentClustDist;
m_MaxCentClustDist = cDfltMaxCentClustDist;
m_MaxClustGrps = cDfltMaxClustGrps;

m_WIGChromID = 0;
m_WIGRptdChromID = 0;
m_WIGSpanLoci = 0;
m_WIGSpanLen = 0;
m_WIGRptdSpanLen = 0;
m_WIGSpanCnts = 0;

m_InNumBuffered = 0;
m_AllocInBuff = 0;

m_OutBuffIdx = 0;	
m_AllocOutBuff = 0;

m_WIGBuffIdx = 0;
m_AllocWIGBuff = 0;

m_LAReadsetNameID = 0;
m_NumReadsetNames = 0;
m_NxtszReadsetIdx = 0;
m_szReadsetNames[0] = '\0';

m_LAChromNameID = 0;
m_NumChromNames = 0;
m_NxtszChromIdx = 0;
m_szChromNames[0] = '\0';

memset(m_Fndrs2Proc,0,sizeof(m_Fndrs2Proc));
m_NumFounders = 0;

m_NumIncludeChroms = 0;
m_NumExcludeChroms = 0;

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
CCallHaplotypes::Process(eModeCSH PMode,	// processing mode 0: report imputation haplotype matrix, 1: report both raw and imputation haplotype matrices, 2: additionally generate GWAS allowing visual comparisons, 3: haplotype grouping, 4: post-process to WIG (default 0)
	            int32_t LimitPrimaryPBAs,                // limit number of loaded primary or founder PBA files to this many. 0: no limits, > 0 sets upper limit
				int32_t ExprID,				// assign this experiment identifier for this PBA analysis
				int32_t GrpHapBinSize,           // when grouping into haplotypes (processing mode 3) then use this non-overlapping bin size for accumulation of differentials between all founder PBAs
	            int32_t NumHapGrpPhases,         // number of phases over which to converge group consensus when haplotype clustering into groups - is a balance between optimal group consensus and processing resource
				int32_t MinCentClustDist,        // haplotype groupings - in processing mode 3 only - minimum group centroid clustering distance (default 500, min 1, clamped to bin size)
				int32_t MaxCentClustDist,        // haplotype groupings - in processing mode 3 only - maximum group centroid clustering distance (default mincentclustdist, clamped to bin size)
				int32_t MaxClustGrps,       // haplotype groupings - in processing mode 3 only - targeted maximum number of groups (default 50)
	            int32_t FndrTrim5,			// trim this many aligned PBAs from 5' end of founder aligned segments - reduces false alleles due to sequencing errors
				int32_t FndrTrim3,			// trim this many aligned PBAs from 3' end of founder aligned segments - reduces false alleles due to sequencing errors
				int32_t ProgTrim5,			// trim this many aligned PBAs from 5' end of progeny aligned segments - reduces false alleles due to sequencing errors
				int32_t ProgTrim3,			// trim this many aligned PBAs from 3' end of progeny aligned segments - reduces false alleles due to sequencing errors
				int32_t WWRLProxWindow,		// proximal window size for Wald-Wolfowitz runs test
				int32_t OutliersProxWindow,	// proximal window size for outliers reduction
				char* pszMaskBPAFile,	    // optional input masking BPA file or Haplotype group clustering specifications file(s)
				int NumFounderInputFiles,	// number of input founder file specs
				char* pszFounderInputFiles[],// names of input founder PBA files (wildcards allowed)
				int NumProgenyInputFiles,	// number of input progeny file specs
				char* pszProgenyInputFiles[],// names of input progeny PBA files (wildcards allowed)
				char* pszOutFile,		    // Windowed haplotype calls output file (CSV format), or WIG if reporting centroid distances
				int	NumIncludeChroms,			// number of chromosome regular expressions to include
				char **ppszIncludeChroms,		// array of include chromosome regular expressions
				int	NumExcludeChroms,			// number of chromosome expressions to exclude
				char **ppszExcludeChroms,		// array of exclude chromosome regular expressions
				int NumThreads)		        // number of worker threads to use
{
int Rslt;
int NumFiles;
int TotNumFiles;
size_t memreq;

Reset();

CreateMutexes();
m_PMode = PMode;
m_LimitPrimaryPBAs = LimitPrimaryPBAs;
m_bAllFndrsLociAligned = true; // replace with actual user parameter!
m_ExprID = ExprID;
m_GrpHapBinSize = GrpHapBinSize;
m_MinCentClustDist = MinCentClustDist;
m_MaxCentClustDist = MaxCentClustDist;
m_MaxClustGrps = MaxClustGrps;
m_NumHapGrpPhases = NumHapGrpPhases;
m_FndrTrim5 = FndrTrim5;
m_FndrTrim3 = FndrTrim3;
m_ProgTrim5 = ProgTrim5;
m_ProgTrim3 = ProgTrim3;
m_WWRLProxWindow = WWRLProxWindow;
m_OutliersProxWindow = OutliersProxWindow;
m_pszRsltsFileBaseName = pszOutFile;
m_NumThreads = NumThreads;

// compile include/exclude chromosome regexpr if user has specified alignments to be filtered by chrom
if(NumIncludeChroms > 0 || NumExcludeChroms > 0)
	if((Rslt = CompileChromRegExprs(NumIncludeChroms,ppszIncludeChroms,NumExcludeChroms,ppszExcludeChroms)) != eBSFSuccess)
		{
		Reset();
		return(Rslt);
		}
m_NumIncludeChroms = NumIncludeChroms;
m_NumExcludeChroms = NumExcludeChroms;
if(NumFounderInputFiles > cMaxFounderFileSpecs)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Can accept at most %d founder wildcarded file specs for processing, %d requested", cMaxFounderFileSpecs, NumFounderInputFiles);
	Reset();
	return(eBSFerrMem);
	}

if(NumProgenyInputFiles > cMaxProgenyFileSpecs)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Can accept at most %d progeny wildcarded file specs for processing, %d requested", cMaxFounderFileSpecs, NumProgenyInputFiles);
	Reset();
	return(eBSFerrMem);
	}

if(pszMaskBPAFile == NULL || pszMaskBPAFile[0] == '\0')
	m_pszMaskBPAFile = NULL;
else
	m_pszMaskBPAFile = pszMaskBPAFile;

if(PMode >= eMCSHGrpDist2WIG)
	{
	Rslt = CSV2WIG(pszMaskBPAFile,pszOutFile);
	Reset();
	return(Rslt);
	}
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
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Memory allocation of %lld bytes through mmap() for chromosome metadata failed - %s", (int64_t)memreq, strerror(errno));
	m_pChromMetadata = NULL;
	Reset();
	return(eBSFerrMem);
	}
#endif
m_AllocdChromMetadataMem = memreq;
m_UsedNumChromMetadata = 0;
m_AllocdChromMetadata = cAllocChromMetadata;

if(!(m_PMode == eMCSHAllelicHapsGrps || m_PMode == eMCSHCoverageHapsGrps))
	{
	memreq = (size_t)cAllocProgenyFndrAligns * sizeof(tsProgenyFndrAligns);
#ifdef _WIN32
	m_pProgenyFndrAligns = (tsProgenyFndrAligns*)malloc(memreq);	// initial and perhaps the only allocation
	if(m_pProgenyFndrAligns == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Initial memory allocation of %lld bytes for chromosome windowed founder counts failed - %s", (int64_t)memreq, strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
#else
	m_pProgenyFndrAligns = (tsProgenyFndrAligns*)mmap(NULL, memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	if(m_pProgenyFndrAligns == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Memory allocation of %lld bytes through mmap() for chromosome windowed founder counts failed - %s", (int64_t)memreq, strerror(errno));
		m_pProgenyFndrAligns = NULL;
		Reset();
		return(eBSFerrMem);
		}
#endif
	m_AllocdProgenyFndrAlignsMem = memreq;
	m_UsedProgenyFndrAligns = 0;
	m_AllocdProgenyFndrAligns = (size_t)cAllocProgenyFndrAligns;


	memreq = (size_t)cAllocAlleleStacks * sizeof(tsAlleleStack);
#ifdef _WIN32
	m_pAlleleStacks = (tsAlleleStack*)malloc(memreq);	// initial and perhaps the only allocation
	if(m_pAlleleStacks == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Initial memory allocation of %lld bytes for allele stacks failed - %s", (int64_t)memreq, strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
#else
	m_pAlleleStacks = (tsAlleleStack*)mmap(NULL, memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	if(m_pAlleleStacks == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Memory allocation of %lld bytes through mmap() for allele stacks failed - %s", (int64_t)memreq, strerror(errno));
		m_pAlleleStacks = NULL;
		Reset();
		return(eBSFerrMem);
		}
#endif
	m_AllocdAlleleStacksMem = memreq;
	m_UsedAlleleStacks = 0;
	m_AllocdAlleleStacks = cAllocAlleleStacks;
	m_AllocdHGBinSpecsMem = 0;
	m_UsedHGBinSpecs = 0;
	m_AllocdHGBinSpecs = 0;
	m_pHGBinSpecs = NULL;
	}
else
	{
	memreq = (size_t)cAllocHGBinSpecs * sizeof(tsHGBinSpec);
#ifdef _WIN32
	m_pHGBinSpecs = (tsHGBinSpec*)malloc(memreq);	// initial and perhaps the only allocation
	if(m_pHGBinSpecs == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Initial memory allocation of %lld bytes for haplotype grouping bin definitions failed - %s", (int64_t)memreq, strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
#else
	m_pHGBinSpecs = (tsHGBinSpec*)mmap(NULL, memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	if(m_pHGBinSpecs == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Memory allocation of %lld bytes through mmap() for haplotype grouping bin definitions failed - %s", (int64_t)memreq, strerror(errno));
		m_pHGBinSpecs = NULL;
		Reset();
		return(eBSFerrMem);
		}
#endif
	m_AllocdHGBinSpecsMem = memreq;
	m_UsedHGBinSpecs = 0;
	m_AllocdHGBinSpecs = cAllocHGBinSpecs;
	m_pProgenyFndrAligns = NULL;
	m_AllocdProgenyFndrAlignsMem = 0;
	m_UsedProgenyFndrAligns = 0;
	m_pAlleleStacks = NULL;
	m_AllocdAlleleStacksMem = 0;
	m_UsedAlleleStacks = 0;
	m_AllocdAlleleStacks = 0;
	}

if((m_pInBuffer = new uint8_t[cInBuffSize]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Memory allocation of %u bytes for input buffering failed - %s", cInBuffSize, strerror(errno));
	Reset();
	return(eBSFerrMem);
	}
m_AllocInBuff = cInBuffSize;
m_InNumBuffered = 0;

gDiagnostics.DiagOut(eDLInfo, gszProcName, "Process: Starting to load pool of founders");
Rslt = eBSFSuccess;		// assume success!
CSimpleGlob glob(SG_GLOB_FULLSORT);
int Idx;
char* pszInFile;
int32_t ReadsetID = 0;
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


	if(LimitPrimaryPBAs > 0 && (TotNumFiles + NumFiles) > LimitPrimaryPBAs)
		{
		NumFiles = LimitPrimaryPBAs - TotNumFiles;
		if(NumFiles == 0)
			break;
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

		if(PMode == eMCSHCoverageHapsGrps)
			ReadsetID = LoadPBACoverage(pszInFile);
		else
			ReadsetID = LoadPBAFile(pszInFile, 0);
		if(ReadsetID <= 0)
			{
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "Process: Errors loading file '%s'",pszInFile);
			Reset();
			return(ReadsetID);
			}
		m_Fndrs2Proc[ReadsetID-1] = 0x01;	// if loaded then assumption is that this founder will be processed
		}
	if(LimitPrimaryPBAs > 0 && TotNumFiles >= LimitPrimaryPBAs)
		break;
	}
m_NumFounders = ReadsetID;
if(m_MaxClustGrps > m_NumFounders)
     m_MaxClustGrps = m_NumFounders;
gDiagnostics.DiagOut(eDLInfo, gszProcName, "Process: Completed loading founders (%d) pool", m_NumFounders);


if(!(m_PMode == eMCSHAllelicHapsGrps || m_PMode == eMCSHCoverageHapsGrps))
	{
	// load the control PBA file if it has been specified
	m_MaskReadsetID = 0;
	if(m_pszMaskBPAFile != NULL)
		{
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Process: Loading control PBA file '%s'", m_pszMaskBPAFile);
		if((ReadsetID = LoadPBAFile(m_pszMaskBPAFile, 2)) <= 0)
			{
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "Process: Errors loading control PBA file");
			Reset();
			return(ReadsetID);
			}
		m_MaskReadsetID = ReadsetID;
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Completed loading control PBA file");
		}

	if((Rslt = GenAlleleStacks(m_NumFounders)) < 1)
		{
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Process: No allele stacks generated");
		Reset();
		return(Rslt);
		}

	if(m_PMode == eMCSHGWAS)
		ReportAnchorsAsGWAS(pszOutFile);
	ReportAnchorsAsCSV(pszOutFile);

	// next is to individually process each progeny PBAs against the founder panel allele stacks
	TotNumFiles = 0;
	m_CurProgenyReadsetID = 0;
	for(Idx = 0; Idx < NumProgenyInputFiles; Idx++)
		{
		glob.Init();
		if(glob.Add(pszProgenyInputFiles[Idx]) < SG_SUCCESS)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to glob progeny '%s", pszProgenyInputFiles[Idx]);
			Reset();
			return(eBSFerrOpnFile);	// treat as though unable to open file
			}
		if((NumFiles = glob.FileCount()) <= 0)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to locate any progeny PBA file matching '%s", pszProgenyInputFiles[Idx]);
			Reset();
			return(eBSFerrOpnFile);	// treat as though unable to open file
			}

		TotNumFiles += NumFiles;
		if(TotNumFiles > cMaxProgenyReadsets)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Can accept at most %d progeny readsets for processing, after wildcard file name expansions there are %d requested", cMaxProgenyReadsets, TotNumFiles);
			Reset();
			return(eBSFerrMem);
			}

		Rslt = eBSFSuccess;
		for(int FileID = 0; Rslt >= eBSFSuccess && FileID < NumFiles; ++FileID)
			{
			pszInFile = glob.File(FileID);
			if((Rslt = ProcessProgenyPBAFile(pszInFile, pszOutFile)) < eBSFSuccess)
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Failed processing progeny PBA file '%s", pszInFile);
				Reset();
				return(Rslt);
				}
			}
		}
	if(m_UsedProgenyFndrAligns == 0)
		{
		gDiagnostics.DiagOut(eDLWarn, gszProcName, "Process: No progeny to  founder allele stack alignments");
		Reset();
		return(eBSFSuccess);
		}

	// sort progeny allele stack overlaps by readset.chrom.loci ascending
	m_mtqsort.SetMaxThreads(m_NumThreads);

	if(m_PMode != eMCSHDefault) // requested to report raw matrix?
		{
		// generate a matrix with chrom.loci (Y) by progeny (X)
		m_mtqsort.qsort(m_pProgenyFndrAligns, (int64_t)m_UsedProgenyFndrAligns, sizeof(tsProgenyFndrAligns), SortProgenyChromLociReadset);
		ReportMatrix(pszOutFile, true);
		}
	m_mtqsort.qsort(m_pProgenyFndrAligns, (int64_t)m_UsedProgenyFndrAligns, sizeof(tsProgenyFndrAligns), SortProgenyFndrAligns);

	for(int32_t ProgIdx = 0; ProgIdx < m_NumProgenies; ProgIdx++)
		{
		if(m_PMode == eMCSHGWAS)
			ReportHaplotypesAsGWAS(pszOutFile, m_ProgenyIDs[ProgIdx], true);
		if(m_PMode != eMCSHDefault)
			ReportHaplotypesByProgeny(pszOutFile, m_ProgenyIDs[ProgIdx], true);

		if(m_WWRLProxWindow > 0)
			{
			ImputeProgenyHeterozygosity(m_WWRLProxWindow, m_ProgenyIDs[ProgIdx]);
			ImputeProgenyHeterozygosity(m_WWRLProxWindow / 5, m_ProgenyIDs[ProgIdx]); // 2nd pass to catch outliers using smaller window
			}

		if(m_OutliersProxWindow > 0)
			ImputeOutliersHaplotypes(m_OutliersProxWindow, m_ProgenyIDs[ProgIdx]);

		if(m_PMode == eMCSHGWAS)
			ReportHaplotypesAsGWAS(pszOutFile, m_ProgenyIDs[ProgIdx]);
		if(m_PMode != eMCSHDefault)
			ReportHaplotypesByProgeny(pszOutFile, m_ProgenyIDs[ProgIdx]);
		}

	// generate a matrix with chrom.loci (Y) by progeny (X)
	m_mtqsort.qsort(m_pProgenyFndrAligns, (int64_t)m_UsedProgenyFndrAligns, sizeof(tsProgenyFndrAligns), SortProgenyChromLociReadset);
	ReportMatrix(pszOutFile);
	}
else
	{
	// optional input masking BPA file or Haplotype group clustering specifications file(s)
	if(m_pszMaskBPAFile != NULL)
		{
		glob.Init();
		if(glob.Add(m_pszMaskBPAFile) < SG_SUCCESS)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to glob haplotype group clustering specification file(s) '%s",m_pszMaskBPAFile);
			Reset();
			return(eBSFerrOpnFile);	// treat as though unable to open file
			}
		if((NumFiles = glob.FileCount()) <= 0)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to locate any haplotype group clustering specification file(s) matching '%s", m_pszMaskBPAFile);
			Reset();
			return(eBSFerrOpnFile);	// treat as though unable to open file
			}
		Rslt = eBSFSuccess;
		for(int FileID = 0; Rslt >= eBSFSuccess && FileID < NumFiles; ++FileID)
			{
			pszInFile = glob.File(FileID);
			if(Rslt = LoadHHGBinSpecs(pszInFile) < 0)
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Failed processing haplotype group clustering specification file '%s", pszInFile);
				Reset();
				return(Rslt);
				}
			}

		if(m_UsedHGBinSpecs == 0)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to parse out any haplotype group clustering specifications from file(s) matching '%s", m_pszMaskBPAFile);
			Reset();
			return(eBSFerrOpnFile);	// treat as though unable to open file
			}

		if(m_UsedHGBinSpecs > 1)
			m_mtqsort.qsort(m_pHGBinSpecs, (int64_t)m_UsedHGBinSpecs, sizeof(tsHGBinSpec), SortHGBinSpecs);

		int32_t CurChromID;
		int32_t HGBinID;
		tsHGBinSpec* pHGBinSpec;
		tsChromMetadata* pChromMetadata;
		CurChromID = 0;
		pHGBinSpec = m_pHGBinSpecs;
		for(HGBinID = 1; HGBinID <= m_UsedHGBinSpecs; HGBinID++, pHGBinSpec++)
			{
			if(pHGBinSpec->ChromID != CurChromID)
				{
				CurChromID = pHGBinSpec->ChromID;
				pChromMetadata = LocateChromMetadataFor(1, pHGBinSpec->ChromID);
				pChromMetadata->HGBinID = HGBinID;
				}
			pHGBinSpec->BinID = HGBinID;
			}
		}

	if((Rslt = GenFounderHaps(m_NumFounders)) < eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Process: Failed generating haplotype groupings");
		Reset();
		return(Rslt);
		}
	}
Reset();
gDiagnostics.DiagOut(eDLInfo, gszProcName, "Completed processing");
return(Rslt);
}


// 
int
CCallHaplotypes::ImputeProgenyHeterozygosity(int32_t MaxDistance, // imputing regional heterozygostic regions having apparent high rates of single haplotype sampling which are within MaxDistance
								int32_t ReadsetID)				// report on this progeny readset only, or if 0 then report on all progeny readsets
{
int32_t CurReadsetID;
int32_t CurChromID;
bool bIsHeterozygote;
int32_t SeqLen;
int32_t NumRuns;
uint8_t ChkHap;
uint8_t PrevHap;
uint32_t NumFa;
uint32_t NumFb;
tsProgenyFndrAligns *pCurPFA;
tsProgenyFndrAligns *pChkPFA;
size_t PFAIdx;
size_t ChkHapStart;
size_t ChkHapEnd;
bool bCheckMe = false;
if(m_UsedProgenyFndrAligns <= 1)
	return(0);

uint32_t NumFndrs[256];

CurReadsetID = 0;
CurChromID = 0;
SeqLen = 0;
NumRuns = 0;
ChkHap = 0;
NumFa = 0;
NumFb = 0;
int NumRandom = 0;
int NumNonrandom = 0;
int NumNotChecked = 0;
pChkPFA = NULL;
pCurPFA = &m_pProgenyFndrAligns[0];
for(PFAIdx = 0; PFAIdx < m_UsedProgenyFndrAligns; PFAIdx++, pCurPFA++)
	{
	if(ReadsetID != 0)
		{
		if(pCurPFA->ReadsetID > ReadsetID)
			break;
		if(pCurPFA->ReadsetID != ReadsetID)
			continue;
		}

	if(pCurPFA->ReadsetID != CurReadsetID || pCurPFA->ChromID != CurChromID)
		{
		CurReadsetID = pCurPFA->ReadsetID;
		CurChromID = pCurPFA->ChromID;
		}

	SeqLen = 0;
	NumRuns = 0;
	ChkHap = 0;
	PrevHap = 0;
	NumFa = 0;
	NumFb = 0;
	uint32_t FaHap = 0;
	uint32_t FbHap = 0;
	memset(NumFndrs, 0, sizeof(NumFndrs));
	ChkHapStart = max(0,(int64_t)PFAIdx-9);
	ChkHapEnd = min((int64_t)m_UsedProgenyFndrAligns-1,(int64_t)PFAIdx+10);
	pChkPFA = &m_pProgenyFndrAligns[ChkHapStart];
	for(; ChkHapStart <= ChkHapEnd; ChkHapStart++,pChkPFA++)
		{
		if(pChkPFA->ReadsetID != CurReadsetID)
			continue;
		if(pChkPFA->NumProgenyFounders == 0)
			continue;
		if(abs(pCurPFA->Loci - pChkPFA->Loci) > MaxDistance)
			continue;

		// must have been at least 1 founder for this progeny..
		// count numbers of each founder
		if(pChkPFA->NumProgenyFounders == 1)
			{
			// only a single founder but which founder?
			ChkHap = 0;
			for(int FndrIdx = 0; ChkHap == 0 && FndrIdx < 4095; FndrIdx++)
				if(BitsVectTest(FndrIdx, pChkPFA->ProgenyFounders))
					ChkHap = FndrIdx+1;

			if(PrevHap == 0 || PrevHap != ChkHap)
				NumRuns++;
			PrevHap = ChkHap;
			NumFndrs[ChkHap-1]++;
			if(NumFndrs[ChkHap - 1] > NumFa)
				{
				NumFa = NumFndrs[ChkHap - 1];
				FaHap = ChkHap;
				}
			else
				if(NumFndrs[ChkHap - 1] > NumFb)
					{
					NumFb = NumFndrs[ChkHap - 1];
					FbHap = ChkHap;
					}
			SeqLen++;
			}
		else // if 2 or more haplotypes present then treat as if 2 sequential haplotyypes - gives a boost to randomness!
			{
			NumRuns+=2;
			SeqLen+=2;
			ChkHap = 0;
			for(int FndrIdx = 0; ChkHap < 2 && FndrIdx < 4095; FndrIdx++)
				if(BitsVectTest(FndrIdx, pChkPFA->ProgenyFounders))
					{
					NumFndrs[FndrIdx]++;
					if(NumFndrs[FndrIdx] > NumFa)
						{
						NumFa = NumFndrs[FndrIdx];
						FaHap = FndrIdx+1;
						}
					else
						if(NumFndrs[FndrIdx] > NumFb)
							{
							NumFb = NumFndrs[FndrIdx];
							FbHap = FndrIdx+1;
							}
					ChkHap++;
					}
			}
		}

	pCurPFA->FaHap = 0;
	pCurPFA->FbHap = 0;
	if(NumRuns >= 3 && SeqLen >= 6) // if possibly sampling a region where there are 2 distinct overlapping haplotypes but sampling results
		// in random single haplotypes being observed then test for random switching between haplotypes
		// and if random switching treat as if heterozygotic region with overlapping haplotypes - Ugh, Ugh and Ugh again. The cons of GBS....
		{
		bIsHeterozygote = m_Stats.IsRandomHaplotypesFaFb(NumFa, NumFb, NumRuns);
		if(bIsHeterozygote)
			{
			pCurPFA->FaHap = FaHap;
			pCurPFA->FbHap = FbHap;
			}
		}
	}

pCurPFA = &m_pProgenyFndrAligns[0];
for(PFAIdx = 0; PFAIdx < m_UsedProgenyFndrAligns; PFAIdx++, pCurPFA++)
	{
	if(ReadsetID != 0)
		{
		if(pCurPFA->ReadsetID > ReadsetID)
			break;
		if(pCurPFA->ReadsetID != ReadsetID)
			continue;
		}

	if(pCurPFA->FaHap != 0 && pCurPFA->FbHap != 0)
		{
		pCurPFA->NumProgenyFounders = 2;
		BitsVectInitialise(0,pCurPFA->ProgenyFounders);
		BitsVectSet(pCurPFA->FaHap-1,pCurPFA->ProgenyFounders);
		BitsVectSet(pCurPFA->FbHap-1,pCurPFA->ProgenyFounders);
		}
	}
return(0);
}


int
CCallHaplotypes::ImputeOutliersHaplotypes(int32_t MaxDistance, // outliers from other called haplotypes which are within MaxDistance
								int32_t ReadsetID)				// report on this progeny readset only, or if 0 then report on all progeny readsets
{
int32_t CurReadsetID;
int32_t CurChromID;
int32_t NxtLociDelta;
int32_t PrevLociDelta;
bool bNxtEqual;
bool bPrevEqual;

tsProgenyFndrAligns *pCurPFA;
tsProgenyFndrAligns *pPrevPFA;
tsProgenyFndrAligns *pNxtPFA;
size_t PFAIdx;

if(m_UsedProgenyFndrAligns <= 1)
	return(0);

CurReadsetID = 0;
CurChromID = 0;

pPrevPFA = NULL;
pNxtPFA = NULL;
PrevLociDelta = 0;
NxtLociDelta = 0;
pCurPFA = &m_pProgenyFndrAligns[0];
for(PFAIdx = 0; PFAIdx < m_UsedProgenyFndrAligns; PFAIdx++, pCurPFA++)
	{
	if(ReadsetID != 0)
		{
		if(pCurPFA->ReadsetID > ReadsetID)
			break;
		if(pCurPFA->ReadsetID != ReadsetID)
			continue;
		}
	if(pCurPFA->ReadsetID != CurReadsetID || pCurPFA->ChromID != CurChromID)
		{
		CurReadsetID = pCurPFA->ReadsetID;
		CurChromID = pCurPFA->ChromID;
		}

	pCurPFA->FaHap = pCurPFA->FbHap = 0;
	if(pCurPFA->NumProgenyFounders == 0 ||	// if previously unable to call then still unable to call!
		pCurPFA->NumProgenyFounders > 1)	// multiploidy calls always accepted with no changes
		continue;
	bNxtEqual = false;
	bPrevEqual = false;
	NxtLociDelta = MaxDistance + 1;
	PrevLociDelta = MaxDistance + 1;
	pNxtPFA = NULL;
	pPrevPFA = NULL;
	if(PFAIdx + 1 < m_UsedProgenyFndrAligns)
		{
		pNxtPFA = &m_pProgenyFndrAligns[PFAIdx + 1];
		if(pNxtPFA->ReadsetID == CurReadsetID &&
		   pNxtPFA->ChromID == CurChromID)
			{
			if((NxtLociDelta = (pNxtPFA->Loci - pCurPFA->Loci)) <= MaxDistance)
				bNxtEqual = BitsVectEqual(pCurPFA->ProgenyFounders, pNxtPFA->ProgenyFounders);
			else
				NxtLociDelta = MaxDistance + 1;
			}
		}

	if(PFAIdx > 0)
		{
		pPrevPFA = &m_pProgenyFndrAligns[PFAIdx - 1];
		if(pPrevPFA->ReadsetID == CurReadsetID &&
		   pPrevPFA->ChromID == CurChromID)
			{
			if((PrevLociDelta = (pCurPFA->Loci - pPrevPFA->Loci)) <= MaxDistance)
				bPrevEqual = BitsVectEqual(pCurPFA->ProgenyFounders, pPrevPFA->ProgenyFounders);
			else
				PrevLociDelta = MaxDistance + 1;
			}
		}

	if(PrevLociDelta == MaxDistance + 1 && NxtLociDelta == MaxDistance + 1)
		continue;
	if(bNxtEqual == true && bPrevEqual == true)
		continue;

	// there is a difference in haplotype calls, change current to be same as nearest of either pPrevPFA or pNxtPFA
	if(PrevLociDelta < NxtLociDelta)
		{
		if(pPrevPFA != NULL && bPrevEqual == false)
			{
			pCurPFA->NumProgenyFounders = pPrevPFA->NumProgenyFounders;
			pCurPFA->ProgenyFounders = pPrevPFA->ProgenyFounders;
			}
		}
	else
		{
		if(pNxtPFA != NULL && bNxtEqual == false)
			{
			pCurPFA->NumProgenyFounders = pNxtPFA->NumProgenyFounders;
			pCurPFA->ProgenyFounders = pNxtPFA->ProgenyFounders;
			}
		}
	}
return(0);
}

int
CCallHaplotypes::ReduceProgenyFounders(int32_t MaxDistance, // where number of founders is more than 2 then attempt to reduce down to a most 2 at any given progeny loci
								int32_t ReadsetID)				// reduce this progeny readset only, or if 0 then all progeny readsets
{
tsProgenyFndrAligns *pCurPFA;
tsProgenyFndrAligns *pNxtPFA;
tsProgenyFndrAligns *pPrvPFA;
size_t PFANxtIdx;
size_t PFAPrvIdx;
size_t PFAIdx;
ts4096Bits CurFounders;
ts4096Bits ExtdFounders;
uint32_t NumIntersect;

int32_t PrevReadsetID;
int32_t PrevChromID;

PrevReadsetID = 0;
PrevChromID = 0;
pCurPFA = &m_pProgenyFndrAligns[0];
for(PFAIdx = 0; PFAIdx < m_UsedProgenyFndrAligns; PFAIdx++, pCurPFA++)
	{
	if(ReadsetID != 0 && pCurPFA->ReadsetID != ReadsetID)
		continue;
	if(pCurPFA->ReadsetID != PrevReadsetID || pCurPFA->ChromID != PrevChromID)
		{
		PrevReadsetID = pCurPFA->ReadsetID;
		PrevChromID = pCurPFA->ChromID;
		}
	if(!pCurPFA->NumProgenyFounders)
		continue;
	if(pCurPFA->NumProgenyFounders <= 2)
		continue;

	// more than 2 putative founders
	// identify which is the unique founder
	CurFounders = pCurPFA->ProgenyFounders;
	tsAlleleStack* pFndrAlleles;
	uint32_t UniqueAllele;
	uint8_t AlleleMsk;

	pFndrAlleles = &m_pAlleleStacks[pCurPFA->AlleleStackID - 1];

	AlleleMsk = 0x03;
	for(UniqueAllele = 0; UniqueAllele < 4; UniqueAllele++, AlleleMsk <<= 2)
		{
		if(pCurPFA->Alleles & AlleleMsk)
			{
			if(pFndrAlleles->NumAlleleFndrs[UniqueAllele] == 1)
				{
				BitsVectClear(CurFounders,pFndrAlleles->Alleles[UniqueAllele]);
				break;
				}
			}
		}
	
	// iterate next/prev until the intersect of current founders has been reduced to <= 2 or > MaxDistance from current loci
	PFAPrvIdx = PFAIdx;
	PFANxtIdx = PFAIdx;
	pNxtPFA = pCurPFA+1;
	pPrvPFA = pCurPFA-1;
	while(PFAPrvIdx != 0 || PFANxtIdx < m_UsedProgenyFndrAligns) 
		{
		if(++PFANxtIdx < m_UsedProgenyFndrAligns)
			{
			if(pNxtPFA->ReadsetID == PrevReadsetID && pNxtPFA->ChromID == PrevChromID && (pNxtPFA->Loci - pCurPFA->Loci) <= MaxDistance)
				{
				// determining the intersect here!
				// extending group of aligned founders through adjacent progeny alignments until one becomes only remaining
				ExtdFounders = pNxtPFA->ProgenyFounders;
				if(NumIntersect = BitsVectIntersect(ExtdFounders, CurFounders) <= 1)
					break;
				pNxtPFA += 1;
				}
			else
				PFANxtIdx = m_UsedProgenyFndrAligns;
			}

		if(PFAPrvIdx-- != 0)
			{
			if(pPrvPFA->ReadsetID == PrevReadsetID && pPrvPFA->ChromID == PrevChromID && (pCurPFA->Loci - pPrvPFA->Loci) <= MaxDistance)
				{
				// determining the intersect here!
				ExtdFounders = pPrvPFA->ProgenyFounders;
				if(NumIntersect = BitsVectIntersect(ExtdFounders, CurFounders) <= 1)
					break;

				pPrvPFA -= 1;
				}
			else
				PFAPrvIdx = 0;
			}
		else
			PFAPrvIdx = 0;
		}
	}
return(0);
}

int 
CCallHaplotypes::ReportHaplotypesByProgeny(char* pszRsltsFileBaseName,		// haplotype results are written to this file base name with '.progeny,m_ExprID.csv' appended
											int32_t ReadsetID,				// report on this progeny readset only, or if 0 then report on all progeny readsets
							  bool bRaw)								// true if reporting on raw haplotypes before any imputing/filtering
{
char szOutFile[_MAX_PATH];
tsProgenyFndrAligns *pCurPFA;
size_t PFAIdx;
int32_t FndrIdx;

int32_t PrevReadsetID;
int32_t PrevChromID;

char *pszProgenyReadset;
char *pszChrom;

if(m_pszOutBuffer == NULL)
	{
	if((m_pszOutBuffer = new uint8_t[cOutBuffSize]) == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "ReportHaplotypesByProgeny: Memory allocation of %u bytes for output buffering failed - %s", cOutBuffSize, strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
	m_AllocOutBuff = cOutBuffSize;
	}
m_OutBuffIdx = 0;

sprintf(szOutFile, "%s.haplotypes.%d.%s.%s.gwas", pszRsltsFileBaseName, m_ExprID, LocateReadset(ReadsetID), bRaw ? "raw" : "imputed");

if(ReadsetID == 0)
	sprintf(szOutFile, "%s.progeny.%d.%s.all.csv", pszRsltsFileBaseName, m_ExprID,bRaw ? "raw" : "imputed");
else
	sprintf(szOutFile, "%s.progeny.%d.%s.%s.csv", pszRsltsFileBaseName, m_ExprID,bRaw ? "raw" : "imputed",LocateReadset(ReadsetID));
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


int32_t 
CCallHaplotypes::ReportAnchorsAsGWAS(char* pszRsltsFileBaseName)		// generate a GWAS format file(s) for viewing marker anchors calls in IGV, useful for looking at associations with coverage etc
{
char szOutFile[_MAX_PATH];
tsAlleleStack* pCurAlleleStack;
char* pszChrom;
int32_t Fndr;
size_t ASIdx;
int32_t AlleleIdx;
uint8_t AlleleMsk;
int32_t Alleles;

if(m_pszOutBuffer == NULL)
	{
	if((m_pszOutBuffer = new uint8_t[cOutBuffSize]) == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "ReportAnchorsAsGWAS: Memory allocation of %u bytes for output buffering failed - %s", cOutBuffSize, strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
	m_AllocOutBuff = cOutBuffSize;
	}
m_OutBuffIdx = 0;
sprintf(szOutFile, "%s.anchors.%d.all.gwas", pszRsltsFileBaseName, m_ExprID);
gDiagnostics.DiagOut(eDLInfo, gszProcName, "Reporting anchors as GWAS to file '%s'", szOutFile);
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
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "ReportAnchorsAsGWAS: Unable to create/truncate %s - %s", szOutFile, strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}
m_OutBuffIdx = sprintf((char*)m_pszOutBuffer, "CHR BP SNP P\n");

int32_t PrevReadsetID = 0;
int32_t PrevChromID = 0;
pCurAlleleStack = &m_pAlleleStacks[0];
for(ASIdx = 0; ASIdx < m_UsedAlleleStacks; ASIdx++, pCurAlleleStack++)
	{
	if(pCurAlleleStack->ChromID != PrevChromID)
		{
		pszChrom = LocateChrom(pCurAlleleStack->ChromID);
		PrevChromID = pCurAlleleStack->ChromID;
		}
	AlleleMsk = 0x03;
	Alleles = 0x00;
	for(AlleleIdx = 3; AlleleIdx >= 0; AlleleIdx--, AlleleMsk <<= 2)
		{
		if(pCurAlleleStack->NumAlleleFndrs[AlleleIdx] == 0)
			continue;
		Alleles |= AlleleMsk;
		}
	if(pCurAlleleStack->NumFndrs == 1)	// if just one founder then which one?
	{
		if(BitsVectTest(0, pCurAlleleStack->ProcFndrs))
			Fndr = 3;
		else
			Fndr = 6;
	}
	else
		Fndr = 9; // both founders
	m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], "%s %d 0x%.4x 0.%d\n", pszChrom, pCurAlleleStack->Loci, Alleles, Fndr);
	if((m_OutBuffIdx + 1000) > m_AllocOutBuff)
		{
		if(!CUtility::RetryWrites(m_hOutFile, m_pszOutBuffer, m_OutBuffIdx))
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "ReportAnchorsAsGWAS: Fatal error in RetryWrites()");
			Reset();
			return(eBSFerrFileAccess);
			}
		m_OutBuffIdx = 0;
		}
	}

if(m_hOutFile != -1)
	{
	if(m_OutBuffIdx)
		{
		if(!CUtility::RetryWrites(m_hOutFile, m_pszOutBuffer, m_OutBuffIdx))
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "ReportAnchorsAsGWAS: Fatal error in RetryWrites()");
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
	}
gDiagnostics.DiagOut(eDLInfo, gszProcName, "Completed reporting anchors as GWAS to file '%s'", szOutFile);
return(eBSFSuccess);
}

int32_t
CCallHaplotypes::ReportAnchorsAsCSV(char* pszRsltsFileBaseName)		// generate a CSV format file(s) containing anchor (potential markers) loci
{
char szOutFile[_MAX_PATH];
char szAlleles[5];
uint32_t BinAlleles;
tsAlleleStack* pCurAlleleStack;
char* pszChrom;
uint32_t FndrIdx;
uint32_t NumFndrs;
size_t ASIdx;
int32_t AlleleIdx;

if(m_pszOutBuffer == NULL)
	{
	if((m_pszOutBuffer = new uint8_t[cOutBuffSize]) == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "ReportAnchorsAsCSV: Memory allocation of %u bytes for output buffering failed - %s", cOutBuffSize, strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
	m_AllocOutBuff = cOutBuffSize;
	}
m_OutBuffIdx = 0;
sprintf(szOutFile, "%s.anchors.%d.all.csv", pszRsltsFileBaseName, m_ExprID);
gDiagnostics.DiagOut(eDLInfo, gszProcName, "Reporting anchors as CVS to file '%s'", szOutFile);
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
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "ReportAnchorsAsCSV: Unable to create/truncate %s - %s", szOutFile, strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}
m_OutBuffIdx = sprintf((char*)m_pszOutBuffer, "\"Chrom\",\"Loci\",\"BinAlleles\",\"Alleles\",\"FndrsAlleleA\",\"FndrsAlleleC\",\"FndrsAlleleG\",\"FndrsAlleleT\"");

int32_t PrevReadsetID = 0;
int32_t PrevChromID = 0;
pCurAlleleStack = &m_pAlleleStacks[0];
for(ASIdx = 0; ASIdx < m_UsedAlleleStacks; ASIdx++, pCurAlleleStack++)
	{
	if(pCurAlleleStack->ChromID != PrevChromID)
		{
		pszChrom = LocateChrom(pCurAlleleStack->ChromID);
		PrevChromID = pCurAlleleStack->ChromID;
		}
	// PBA at a single loci is packed into a single 8bit byte: Allele A in bits 7.6, C in bits 5.4, G in bits 3.2, T in bits 1.0
	BinAlleles = 0;
	for(AlleleIdx = 3; AlleleIdx >= 0; AlleleIdx--)
		{
		if(pCurAlleleStack->NumAlleleFndrs[AlleleIdx] == 0)
			{
			szAlleles[3-AlleleIdx] = '.';
			continue;
			}
		switch(AlleleIdx) {
			case 3:
				szAlleles[3-AlleleIdx] = 'A';
				BinAlleles |= 0x01;
				break;
			case 2:
				szAlleles[3-AlleleIdx] = 'C';
				BinAlleles |= 0x02;
				break;
			case 1:
				szAlleles[3-AlleleIdx] = 'G';
				BinAlleles |= 0x04;
				break;
			case 0:
				szAlleles[3-AlleleIdx] = 'T';
				BinAlleles |= 0x08;
				break;
			}
		}
	szAlleles[4] = '\0';
	m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], "\n%s,%d,%d,\"%s\"", pszChrom, pCurAlleleStack->Loci, BinAlleles, szAlleles);
	for(AlleleIdx = 3; AlleleIdx >= 0; AlleleIdx--)
		{
		if((NumFndrs = pCurAlleleStack->NumAlleleFndrs[AlleleIdx]) == 0)
			{
			m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], ",\"\"");
			continue;
			}
		m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], ",\"");
		for(FndrIdx = 0; NumFndrs && FndrIdx < pCurAlleleStack->NumFndrs; FndrIdx++)
			{
			if(BitsVectTest(FndrIdx, pCurAlleleStack->Alleles[AlleleIdx]))
				{
				m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], "%s", LocateReadset(FndrIdx + 1));
				if(--NumFndrs)
					m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], ":");
				}
			}
		m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], "\"");
		}

	
	if((m_OutBuffIdx + 1000) > m_AllocOutBuff)
		{
		if(!CUtility::RetryWrites(m_hOutFile, m_pszOutBuffer, m_OutBuffIdx))
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "ReportAnchorsAsCSV: Fatal error in RetryWrites()");
			Reset();
			return(eBSFerrFileAccess);
			}
		m_OutBuffIdx = 0;
		}
	}

if(m_hOutFile != -1)
	{
	if(m_OutBuffIdx)
		{
		if(!CUtility::RetryWrites(m_hOutFile, m_pszOutBuffer, m_OutBuffIdx))
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "ReportAnchorsAsCSV: Fatal error in RetryWrites()");
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
	}
gDiagnostics.DiagOut(eDLInfo, gszProcName, "Completed reporting anchors as CVS to file '%s'", szOutFile);
return(eBSFSuccess);
}

int32_t 
CCallHaplotypes::ReportHaplotypesAsGWAS(char* pszRsltsFileBaseName,		// generate a GWAS format file(s) for viewing haplotype calls in IGV, usefull for looking at associations with coverage etc
							  int32_t ReadsetID,						// report on this progeny readset
							  bool bRaw)								// true if reporting on raw haplotypes before any imputing/filtering
{		
char szOutFile[_MAX_PATH];
tsProgenyFndrAligns* pCurPFA;
char *pszFounders;
int32_t FndOfs;
int32_t FndrIdx;
char *pszProgenyReadset;
char *pszChrom;
int32_t Fndr;
size_t PFAIdx;

if(m_pszOutBuffer == NULL)
	{
	if((m_pszOutBuffer = new uint8_t[cOutBuffSize]) == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "ReportHaplotypesAsGWAS: Memory allocation of %u bytes for output buffering failed - %s", cOutBuffSize, strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
	m_AllocOutBuff = cOutBuffSize;
	}
m_OutBuffIdx = 0;
gDiagnostics.DiagOut(eDLInfo, gszProcName, "Reporting haplotypes as GWAS to file '%s'", szOutFile);
sprintf(szOutFile, "%s.haplotypes.%d.%s.%s.gwas", pszRsltsFileBaseName, m_ExprID, LocateReadset(ReadsetID), bRaw ? "raw" : "imputed");
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
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "ReportHaplotypesAsGWAS: Unable to create/truncate %s - %s", szOutFile, strerror(errno));
		Reset();
		return(eBSFerrCreateFile);
	}
m_OutBuffIdx = sprintf((char*)m_pszOutBuffer, "CHR BP SNP P\n");

if((pszFounders = new char[cMaxFounderReadsets*50]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "ReportHaplotypesAsGWAS: Memory allocation of %u bytes failed",cMaxFounderReadsets*50);
	Reset();
	return(eBSFerrMem);
	}
int32_t PrevReadsetID = 0;
int32_t PrevChromID = 0;
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
	FndOfs = 0;
	for(FndrIdx = 0; FndrIdx < m_NumFounders; FndrIdx++)
		{
		if(BitsVectTest(FndrIdx, pCurPFA->ProgenyFounders))
			{
			if(FndOfs != 0)
				pszFounders[FndOfs++] = ':';
			FndOfs+=sprintf(&pszFounders[FndOfs],"%s", LocateReadset(FndrIdx + 1));
			}
		}

	if(pCurPFA->NumProgenyFounders == 1)	// if just one founder then which one?
		{
		if(BitsVectTest(0, pCurPFA->ProgenyFounders))
			Fndr = 3; // visually represents Fa parental haplotype only 
		else
			Fndr = 9; // visually represents Fb parental haplotype only
		}
	else
		Fndr = 1;	// visually represents both parental haplotypes as being present 
	
	m_OutBuffIdx += sprintf((char *)&m_pszOutBuffer[m_OutBuffIdx],"%s %d %s 0.%d\n", pszChrom, pCurPFA->Loci, pszFounders,Fndr);
	if((m_OutBuffIdx + 1000) > m_AllocOutBuff)
		{
		if(!CUtility::RetryWrites(m_hOutFile, m_pszOutBuffer, m_OutBuffIdx))
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "ReportHaplotypesAsGWAS: Fatal error in RetryWrites()");
			delete[]pszFounders;
			Reset();
			return(eBSFerrFileAccess);
			}
		m_OutBuffIdx = 0;
		}
	}
delete[]pszFounders;
if(m_hOutFile != -1)
	{
	if(m_OutBuffIdx)
		{
		if(!CUtility::RetryWrites(m_hOutFile, m_pszOutBuffer, m_OutBuffIdx))
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "ReportHaplotypesAsGWAS: Fatal error in RetryWrites()");
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
	}

gDiagnostics.DiagOut(eDLInfo, gszProcName, "Completed reporting haplotypes as GWAS to file '%s'", szOutFile);
return(eBSFSuccess);
}

int			// returned number of PBAs which are non-conformant
CCallHaplotypes::ValidatePBAs(int32_t Length,
			 uint8_t* pPBAs,
	         bool bSetNoAlleles, // if non-conformant then overwrite *pAlleles to be no alleles present
	         bool bNormalise)    // normalise for coverage, low coverage major alleles (2,0,0,0) are normalised to high coverage (3,0,0,0)
{
int32_t Idx;
int32_t NumErrors = 0;
for(Idx=0; Idx < Length; Idx++,pPBAs++)
	{
	if(!ValidatePBA(pPBAs,bSetNoAlleles,bNormalise))
		NumErrors++;
	}
return(NumErrors);
}

bool 
CCallHaplotypes::ValidatePBA(uint8_t *pAlleles,	// validate that PBA alleles are properly conformant
		         bool bSetNoAlleles, // if non-conformant then overwrite *pAlleles to be no alleles present
	         bool bNormalise)    // normalise for coverage, low coverage major alleles (2,0,0,0) are normalised to high coverage (3,0,0,0)
{
uint32_t AlleleIdx;
uint8_t Allele;
uint8_t AlleleMsk;
int32_t Num2s;
int32_t Num1s;

if(*pAlleles == 0)
	return(true);

Num2s = Num1s = 0;
AlleleMsk = 0x03;
for(AlleleIdx = 0; AlleleIdx < 4; AlleleIdx++, AlleleMsk <<= 2)
	{
	Allele = (*pAlleles & AlleleMsk) >> (AlleleIdx * 2);
	switch(Allele) {
		case 0x03:
			if((*pAlleles & ~AlleleMsk)!=0) // can only be one major ...
				{
				if(bSetNoAlleles)
					*pAlleles = 0;
				return(false);
				}
			else
				return(true);
		case 0x02:
			if(bNormalise && (*pAlleles & ~AlleleMsk)==0) // if low coverage, a major allele is represented as 2,0,0,0 whereas if high coverage then that same allele would have been represented as 3,0,0,0
				{
				*pAlleles = AlleleMsk; // normalisation transforms low coverage major alleles into same representation as high coverage
				return(true);
				}
			Num2s++;
			continue;
		case 0x01:
			Num1s++;
			continue;
		}
	}
if(Num2s > 2 || Num1s > 2 || (Num2s+Num1s > 2))
	{
	if(bSetNoAlleles)
		*pAlleles = 0;
	return(false);
	}
return(true);
}


int
CCallHaplotypes::ReportMatrix(char* pszRsltsFileBaseName,		// matrix results are written to this file base name with 'matrix.csv' appended
			bool bRaw)								// true if reporting on raw haplotypes before any imputing/filtering
{
char szOutFile[_MAX_PATH];
tsProgenyFndrAligns* pCurPFA;
size_t PFAIdx;
int32_t FndrIdx;
int32_t ProgIdx;
int32_t ProgenyHaplotypes[cMaxProgenyReadsets];
uint32_t ProgHaplotype;
int32_t CurChromID;
int32_t CurLoci;
bool bNewLoci;
bool bProgFndrs;

char* pszChrom;

if(m_pszOutBuffer == NULL)
{
	if((m_pszOutBuffer = new uint8_t[cOutBuffSize]) == NULL)
	{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "ReportMatrix: Memory allocation of %u bytes for output buffering failed - %s", cOutBuffSize, strerror(errno));
		Reset();
		return(eBSFerrMem);
	}
	m_AllocOutBuff = cOutBuffSize;
}
gDiagnostics.DiagOut(eDLInfo, gszProcName, "Reporting matrix as CSV to file '%s'", szOutFile);
m_OutBuffIdx = 0;
sprintf(szOutFile, "%s.%d.%s.matrix.csv", pszRsltsFileBaseName,m_ExprID,bRaw ? "raw" : "inputed");
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
pCurPFA = &m_pProgenyFndrAligns[0];
for(PFAIdx = 0; PFAIdx < m_UsedProgenyFndrAligns; PFAIdx++, pCurPFA++)
	{
	if(pCurPFA->ChromID != CurChromID || pCurPFA->Loci != CurLoci)
		{
		if(bNewLoci && bProgFndrs)
			{
			m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], "\n%d,\"%s\",%u", m_ExprID,pszChrom, CurLoci);
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
	m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], "\n%d,\"%s\",%u", m_ExprID,pszChrom, CurLoci);
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
gDiagnostics.DiagOut(eDLInfo, gszProcName, "Completed reporting matrix as CSV to file '%s'", szOutFile);
return(eBSFSuccess);
}


int
CCallHaplotypes::ProcessProgenyPBAFile(char* pszProgenyPBAFile,	// load, process and call haplotypes for this progeny PBA file against previously loaded panel founder PBAs
						char *pszRsltsFileBaseName)			// results are written to this file base name with progeny readset identifier and type appended
{
m_CurProgenyReadsetID = 0;
m_LAChromNameID = 0;
m_LAReadsetNameID = 0;
if((m_CurProgenyReadsetID = LoadPBAFile(pszProgenyPBAFile,1)) <= 0)
	{
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessProgenyPBAFile: Errors loading progeny PBA file '%s'",pszProgenyPBAFile);
	Reset();
	return(-1);
	}
gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessProgenyPBAFile: Completed loading progeny PBA file");
m_ProgenyIDs[m_NumProgenies++] = m_CurProgenyReadsetID;
if(m_pszOutBuffer == NULL)
	{
	if((m_pszOutBuffer = new uint8_t[cOutBuffSize]) == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessProgenyPBAFile: Memory allocation of %u bytes for output buffering failed - %s", cOutBuffSize, strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
	m_AllocOutBuff = cOutBuffSize;
	}
m_OutBuffIdx = 0;

char *pszProgenyReadset = LocateReadset(m_CurProgenyReadsetID);

tsReadsetMetadata *pProgeny = &m_Readsets[m_CurProgenyReadsetID -1];
tsProgenyFndrAligns ProgenyFndrAligns;
int32_t NumFndrUniques;
int32_t NumOverlapping;
uint8_t ProgenyAlleles;
uint8_t Allele;
bool bAlleleAccept;
uint32_t AlleleStackIdx;
tsAlleleStack *pCurAlleleStack;
uint32_t NumAcceptAlleles;
uint32_t NumPotentalFndrs;
uint32_t AlleleIdx;
uint8_t AlleleMsk;

NumFndrUniques = 0;
NumOverlapping = 0;
pCurAlleleStack = m_pAlleleStacks;
for(AlleleStackIdx = 0; AlleleStackIdx < m_UsedAlleleStacks; AlleleStackIdx++, pCurAlleleStack++)
	{
	NumFndrUniques = 0;
	NumPotentalFndrs = 0;
	memset(&ProgenyFndrAligns, 0, sizeof(tsProgenyFndrAligns));
	ProgenyFndrAligns.AlleleStackID = pCurAlleleStack->AlleleStackID;
	ProgenyFndrAligns.ChromID = pCurAlleleStack->ChromID;
	ProgenyFndrAligns.Loci = pCurAlleleStack->Loci;
	ProgenyFndrAligns.ReadsetID = m_CurProgenyReadsetID;
	ProgenyAlleles = LocateReadsetChromLociAlleles(m_CurProgenyReadsetID, pCurAlleleStack->ChromID,pCurAlleleStack->Loci);
	ProgenyFndrAligns.Alleles = ProgenyAlleles;

	if(ProgenyAlleles != 0)		// if coverage
		{
		NumAcceptAlleles = 0;
		bAlleleAccept = false;

		// Note: Allele A in bits 7.6, C in bits 5.4, G in bits 3.2, T in bits 1.0
		AlleleMsk = 0x03;
		for(AlleleIdx = 0; AlleleIdx < 4; AlleleIdx++, AlleleMsk <<= 2) // iterating from allele T (AlleleIdx = 0) through to allele A (AlleleIdx = 3)
			{
			Allele = (ProgenyAlleles & AlleleMsk) >> (AlleleIdx * 2);

			if(Allele == 0)
				continue;

			// accepting progeny with alleles matching:
			// 3 0 0 0	-> single dirac present in at least 0.8 (cScorePBA3MinProp) of alignment coverage
			// 2 0 0 0	-> one allele present in at least 0.70 (cScorePBA2MinLCProp) proportion of alignment low coverage, otherwise there would be more than 1 allele
			// 2 2 0 0	-> two alleles, each present in at least 0.35 (cScorePBA2MinProp) of alignment coverage
			// rational for only accepting certain combinations is to reduce probability of false haplotype calls due to allele noise

			if(Allele == 1) // can't accept very low coverage alleles
				{
				bAlleleAccept = false;
				break;
				}

			if(pCurAlleleStack->NumAlleleFndrs[AlleleIdx]==0) // if no founders having this allele then treat as introgression and treat as if no progeny allele
				{
				bAlleleAccept = false;
				break;
				}

			// totalling number of founder allele which are unique relative to other founders
			if(pCurAlleleStack->NumAlleleFndrs[AlleleIdx] == 1) // later will check to ensure that there was a least 1 unique founder
				NumFndrUniques++;

			NumPotentalFndrs += pCurAlleleStack->NumAlleleFndrs[AlleleIdx]; // when more than 2 (diploid) founders then can't disambiguate so discard

			if(Allele >= 2 && (ProgenyAlleles & ~AlleleMsk) == 0) // if only founder then founder exclusively identified
				{
				NumAcceptAlleles = 1;
				bAlleleAccept = true;
				ProgenyFndrAligns.NumProgenyFounders = BitsVectUnion(ProgenyFndrAligns.ProgenyFounders, pCurAlleleStack->Alleles[AlleleIdx]);
				break;
				}

			if(Allele == 3 && (ProgenyAlleles & ~AlleleMsk)!=0)			// can only accept diracs if exclusive
				{
				bAlleleAccept = false;
				break;
				}

			// accepting any combination of alleles as long as there is no more than 2
			if(++NumAcceptAlleles > 2)					// can only be at most 2 alleles- progeny assumed to be diploid so additional alleles treated as noise
				{
				bAlleleAccept = false;
				break;
				}

			ProgenyFndrAligns.NumProgenyFounders = BitsVectUnion(ProgenyFndrAligns.ProgenyFounders, pCurAlleleStack->Alleles[AlleleIdx]);
			bAlleleAccept = true;
			}

		if(!bAlleleAccept || NumFndrUniques < 1)
			continue;

		if(NumPotentalFndrs > 2)
			continue;

		NumOverlapping++;
		if(AddProgenyFndrAligns(&ProgenyFndrAligns)<=0)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessProgenyPBAFile: Failed allocating progeny to founder allele stack");
			Reset();
			return(eBSFerrMem);	// assuming memory allocation error
			}
		}
	}
	
// remove all PBAs for this progeny as no longer required and memory resources could be limited
tsChromMetadata *pChromMetadata;
uint32_t CurChromMetadataIdx = pProgeny->StartChromMetadataIdx;
for(int32_t ChromIdx = 0; ChromIdx < pProgeny->NumChroms && CurChromMetadataIdx != 0; ChromIdx++)
	{
	pChromMetadata = &m_pChromMetadata[CurChromMetadataIdx - 1];
	if(pChromMetadata->pPBAs != NULL)
#ifdef _WIN32
		free(pChromMetadata->pPBAs);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(pChromMetadata->pPBAs != MAP_FAILED)
		munmap(pChromMetadata->pPBAs, pChromMetadata->ChromLen);
#endif
	pChromMetadata->pPBAs = NULL;
	}
gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessProgenyPBAFile: Located %d progeny loci overlapping %d founder AlleleStacks loci", NumOverlapping, m_UsedAlleleStacks);
return(eBSFSuccess);
}

int
CCallHaplotypes::CompileChromRegExprs(int	NumIncludeChroms,	// number of chromosome regular expressions to include
		char **ppszIncludeChroms,		// array of include chromosome regular expressions
		int	NumExcludeChroms,			// number of chromosome expressions to exclude
		char **ppszExcludeChroms)		// array of exclude chromosome regular expressions
{
int Idx = 0;

#ifndef _WIN32
int RegErr;				// regular expression parsing error
char szRegErr[128];			// to hold RegErr as textual representation ==> regerror();
#endif

#ifdef _WIN32
try {
	for(Idx=0;Idx < NumIncludeChroms;Idx++)
		{
		m_IncludeChromsRE[Idx] = new Regexp();
		m_IncludeChromsRE[Idx]->Parse(ppszIncludeChroms[Idx],false);	// note case insensitive
		}
	}
catch(...)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to process include regexpr chrom '%s'",ppszIncludeChroms[Idx]);
	return(eBSFerrMem);
	}
try {
	for(Idx=0;Idx < NumExcludeChroms;Idx++)
		{
		m_ExcludeChromsRE[Idx] = new Regexp();
		m_ExcludeChromsRE[Idx]->Parse(ppszExcludeChroms[Idx],false);	// note case insensitive
		}
	}
catch(...)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to process exclude regexpr chrom '%s'",ppszExcludeChroms[Idx]);
	return(eBSFerrMem);
	}

#else
for(Idx=0;Idx < NumIncludeChroms;Idx++)
	{

	RegErr=regcomp(&m_IncludeChromsRE[Idx],ppszIncludeChroms[Idx],REG_EXTENDED | REG_ICASE);	// note case insensitive
	if(RegErr)
		{
		regerror(RegErr,&m_IncludeChromsRE[Idx],szRegErr,sizeof(szRegErr));
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to process include chrom '%s' error: %s",ppszIncludeChroms[Idx],szRegErr);
		return(eBSFerrMem);
		}
	}
for(Idx=0;Idx < NumExcludeChroms;Idx++)
	{
	RegErr = regcomp(&m_ExcludeChromsRE[Idx],ppszExcludeChroms[Idx],REG_EXTENDED | REG_ICASE);	// note case insensitive
	if(RegErr)
		{
		regerror(RegErr,&m_ExcludeChromsRE[Idx],szRegErr,sizeof(szRegErr));
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to process exclude chrom '%s' error: %s",ppszExcludeChroms[Idx],szRegErr);
		return(eBSFerrMem);
		}
	}
#endif
return(eBSFSuccess);
}

bool					// true if chrom is accepted, false if chrom not accepted
CCallHaplotypes::AcceptThisChromID(uint32_t ChromID)
{
char *pzChrom;
if(!(m_NumExcludeChroms || m_NumIncludeChroms))
	return(true);
if((pzChrom=LocateChrom(ChromID)) == NULL)
    return(false);
return(AcceptThisChromName(pzChrom));
}

bool					// true if chrom is accepted, false if chrom not accepted
CCallHaplotypes::AcceptThisChromName(char *pszChrom)
{
int IncChromIdx;
int ExclChromIdx;

bool bProcChrom = false;
int MatchesFiltOut = 0;

if(!(m_NumExcludeChroms || m_NumIncludeChroms))
	return(true);

#ifdef _WIN32
RegexpMatch mc;
#else
regmatch_t mc;
int RegErr;					// regular expression parsing error
char szRegErr[128];			// to hold RegErr as textual representation ==> regerror();
#endif

AcquireSerialise();

	// check if to be excluded
bProcChrom = true;
for(ExclChromIdx = 0; ExclChromIdx < m_NumExcludeChroms; ExclChromIdx++)
	{
#ifdef _WIN32
	if(m_ExcludeChromsRE[ExclChromIdx]->Match(pszChrom,&mc))
#else
	if(!regexec(&m_ExcludeChromsRE[ExclChromIdx],pszChrom,1,&mc,0))
#endif
		{
		bProcChrom = false;
		break;
		}
	}

	// to be included?
if(bProcChrom && m_NumIncludeChroms > 0)
	{
	bProcChrom = false;
	for(IncChromIdx = 0; IncChromIdx < m_NumIncludeChroms; IncChromIdx++)
		{
#ifdef _WIN32
		if(m_IncludeChromsRE[IncChromIdx]->Match(pszChrom,&mc))
#else
		if(!regexec(&m_IncludeChromsRE[IncChromIdx],pszChrom,1,&mc,0))
#endif
			{
			bProcChrom = true;
			break;
			}
		}
	}

ReleaseSerialise();

if(!bProcChrom)
	return(false);
return(true);
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
							uint8_t ReadsetType)	// 0: founder, 1: progeny, 2: control
{
int Rslt;
int Version;
char szExperimentID[100];
char szRefAssemblyID[100];
char szReadsetID[100];
uint32_t PrevChromMetadataIdx;
int32_t ChromID;
int32_t ReadsetID;
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

// check file type is PBA
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
if((ReadsetID = AddReadset(szReadsetID, ReadsetType))==0)
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
int32_t ChromLen;

// iterate over all chromosomes
PrevChromMetadataIdx = 0;
while((m_InNumProcessed + 110) <= m_InNumBuffered)
	{
	pBuff = &m_pInBuffer[m_InNumProcessed++];
	ChromNameLen = (int)*pBuff++;
	pszChromName = (char *)pBuff;
	pBuff += 1+ChromNameLen;
	ChromLen = *(int32_t *)pBuff;
	pBuff+=4;
	m_InNumProcessed += ChromNameLen + 1 + sizeof(uint32_t);
	if((m_InNumBuffered - m_InNumProcessed) < (uint32_t)ChromLen && FillInBuffer((uint32_t)ChromLen) == 0)
		break;

	// check if this chromosome is to be retained for further processing
	if(!AcceptThisChromName(pszChromName))
		{
		m_InNumProcessed += ChromLen;
		if((m_InNumBuffered - m_InNumProcessed) < 110 && FillInBuffer(ChromLen) == 0)
			break;
		continue;
		}

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
	pChromMetadata->HGBinID = 0;

//	validate PBA allele composition, earlier releases were retaining very low allele proportions so validating and in-place replacing these incorrect PBAs with non-alignments
	int NumErrs = ValidatePBAs(ChromLen,pBuff,true,true);
//
	memcpy(pChromMetadata->pPBAs, pBuff, ChromLen);
	if(ReadsetType == 0)
		TrimPBAs(m_FndrTrim5, m_FndrTrim3, ChromLen,pChromMetadata->pPBAs);
	else // treating controls as if progeny when trimming
		TrimPBAs(m_ProgTrim5, m_ProgTrim3, ChromLen, pChromMetadata->pPBAs);
	m_InNumProcessed += ChromLen;
	if((m_InNumBuffered - m_InNumProcessed) < 110 && FillInBuffer(ChromLen) == 0)
		break;
	}

close(m_hInFile);
m_hInFile = -1;
return((int32_t)ReadsetID);
}

// trim 5' and 3' aligned segments within the PBAs attempting to reduce sequencing error induced false alleles
int			// returns number of non-trimmed loci in the pPBAs
CCallHaplotypes::TrimPBAs(int32_t Trim5,	// trim 5' this many aligned PBA bases from each aligned segment
		 int32_t Trim3,	// trim 3' this many aligned PBA bases from each aligned segment
		 int32_t PBALen,	// pPBAs contains this many packed base alleles
		uint8_t* pPBAs)		// packed base alleles to be processed
{
int32_t NonTrimmed;
int32_t ToTrim;
int32_t Loci;
uint8_t* pPBA;
if(pPBAs == NULL || PBALen == 0)
	return(0);
if(Trim5 == 0 && Trim3 == 0)
	return(PBALen);
NonTrimmed = PBALen;

// first the 5' trim
if(Trim5)
	{
	ToTrim = Trim5;
	pPBA = pPBAs;
	for(Loci=0;Loci < PBALen;Loci++, pPBA++)
		{
		if(*pPBA == 0)
			{
			ToTrim = Trim5;
			continue;
			}
		if(ToTrim)
			{
			ToTrim-=1;
			NonTrimmed-=1;
			*pPBA = 0;
			}
		}
	}
// then the 3' trim
if(Trim3)
	{
	ToTrim = Trim3;
	pPBA = &pPBAs[PBALen-1];
	for(Loci = 0; Loci < PBALen; Loci++, pPBA--)
		{
		if(*pPBA == 0)
			{
			ToTrim = Trim3;
			continue;
			}
		if(ToTrim)
			{
			NonTrimmed-=1;
			ToTrim -= 1;
			*pPBA = 0;
			}
		}
	}
return(NonTrimmed);
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
CCallHaplotypes::AllocPBAs(int32_t ChromLen)	// allocate memory to hold at least this many packed base alleles
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

if(ReadSetID > m_NumReadsetNames || ReadSetID == 0)
	return(NULL);
pReadsetMetadata = &m_Readsets[ReadSetID-1];
CurChromMetadataIdx = pReadsetMetadata->StartChromMetadataIdx;
for(int32_t ChromIdx = 0; ChromIdx < pReadsetMetadata->NumChroms && CurChromMetadataIdx != 0; ChromIdx++)
	{
	pChromMetadata = &m_pChromMetadata[CurChromMetadataIdx - 1];
	if(pChromMetadata->ChromID == ChromID)
		return(pChromMetadata->pPBAs);
	CurChromMetadataIdx = pChromMetadata->NxtChromMetadataIdx;
	}
return(NULL);
}

tsChromMetadata *								// returned pointer to chromosome metadata
CCallHaplotypes::LocateChromMetadataFor(int32_t ReadSetID,		// readset identifier 
			 int32_t ChromID)			// chrom identifier
{
tsReadsetMetadata *pReadsetMetadata;
tsChromMetadata *pChromMetadata;
uint32_t CurChromMetadataIdx;

if(ReadSetID > m_NumReadsetNames || ReadSetID == 0)
	return(NULL);
pReadsetMetadata = &m_Readsets[ReadSetID-1];
CurChromMetadataIdx = pReadsetMetadata->StartChromMetadataIdx;
for(int32_t ChromIdx = 0; ChromIdx < pReadsetMetadata->NumChroms && CurChromMetadataIdx != 0; ChromIdx++)
	{
	pChromMetadata = &m_pChromMetadata[CurChromMetadataIdx - 1];
	if(pChromMetadata->ChromID == ChromID)
		return(pChromMetadata);
	CurChromMetadataIdx = pChromMetadata->NxtChromMetadataIdx;
	}
return(NULL);
}




size_t									// returned index+1 into  m_pProgenyFndrAligns[] to allocated and initialised ProgenyFndrAligns, 0 if errors
CCallHaplotypes::AddProgenyFndrAligns(tsProgenyFndrAligns *pInitProgenyFndrAligns)	// allocated tsProgenyFndrAligns to be initialised with a copy of pInitProgenyFndrAligns
{
size_t ToAllocdProgenyFndrAligns;
tsProgenyFndrAligns *pProgenyFndrAligns;
size_t memreq;
AcquireSerialise();
if (m_pProgenyFndrAligns == NULL)					// may be NULL first time in
	{
	memreq = cAllocProgenyFndrAligns * sizeof(tsProgenyFndrAligns);
#ifdef _WIN32
	m_pProgenyFndrAligns = (tsProgenyFndrAligns*)malloc((size_t)memreq);
	if (m_pProgenyFndrAligns == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddWinBinCnts: Memory allocation of %lld bytes failed", (int64_t)memreq);
		return(0);
		}
#else
	m_pProgenyFndrAligns = (tsProgenyFndrAligns*)mmap(NULL, (size_t)memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	if (m_pProgenyFndrAligns == MAP_FAILED)
		{
		m_pProgenyFndrAligns = NULL;
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddWinBinCnts: Memory allocation of %lld bytes through mmap()  failed", (int64_t)memreq, strerror(errno));
		return(0);
		}
#endif
	m_AllocdProgenyFndrAlignsMem = memreq;
	m_AllocdProgenyFndrAligns = (size_t)cAllocProgenyFndrAligns;
	m_UsedProgenyFndrAligns = 0;
	}
else
		// needing to allocate more memory?
	if ((m_UsedProgenyFndrAligns) >= m_AllocdProgenyFndrAligns)
		{
		ToAllocdProgenyFndrAligns = m_UsedProgenyFndrAligns + (size_t)cAllocProgenyFndrAligns;
		size_t memreq = ToAllocdProgenyFndrAligns * sizeof(tsProgenyFndrAligns);
#ifdef _WIN32
		pProgenyFndrAligns = (tsProgenyFndrAligns*)realloc(m_pProgenyFndrAligns, memreq);
		if (pProgenyFndrAligns == NULL)
			{
#else
		pProgenyFndrAligns = (tsProgenyFndrAligns*)mremap(m_pProgenyFndrAligns, m_AllocdProgenyFndrAlignsMem, memreq, MREMAP_MAYMOVE);
		if (pProgenyFndrAligns == MAP_FAILED)
			{
#endif
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddProgenyFndrAligns: Memory reallocation to %lld bytes failed - %s", memreq, strerror(errno));
			return(0);
			}
		m_pProgenyFndrAligns = pProgenyFndrAligns;
		m_AllocdProgenyFndrAlignsMem = memreq;
		m_AllocdProgenyFndrAligns = ToAllocdProgenyFndrAligns;
		}
pProgenyFndrAligns = &m_pProgenyFndrAligns[m_UsedProgenyFndrAligns++];
if(pInitProgenyFndrAligns != NULL)
	*pProgenyFndrAligns = *pInitProgenyFndrAligns;
else
	memset(pProgenyFndrAligns,0,sizeof(tsProgenyFndrAligns));
ReleaseSerialise();
return(m_UsedProgenyFndrAligns);
}

void
CCallHaplotypes::AcquireFastSerialise(void)
{
int SpinCnt = 500;
int BackoffMS = 5;

#ifdef _WIN32
while(InterlockedCompareExchange(&m_FastSerialise,1,0)!=0)
	{
	if(SpinCnt -= 1)
		continue;
	CUtility::SleepMillisecs(BackoffMS);
	SpinCnt = 100;
	if(BackoffMS < 500)
		BackoffMS += 2;
	}
#else
while(__sync_val_compare_and_swap(&m_FastSerialise,0,1)!=0)
	{
	if(SpinCnt -= 1)
		continue;
	CUtility::SleepMillisecs(BackoffMS);
	SpinCnt = 100;
	if(BackoffMS < 500)
		BackoffMS += 2;
	}
#endif
}

void
CCallHaplotypes::ReleaseFastSerialise(void)
{
#ifdef _WIN32
InterlockedCompareExchange(&m_FastSerialise,0,1);
#else
__sync_val_compare_and_swap(&m_FastSerialise,1,0);
#endif
}

int32_t									// returned index+1 into m_pHGBinSpecs[] to allocated and initialised allele stack, 0 if errors
CCallHaplotypes::AddHGBinSpec(int32_t ChromID,	    // haplotype group clustering is on this chromosome
	int32_t StartLoci,      // bin starts at this loci
	int32_t NumLoci,	    // covering this many loci
	int32_t MinCentroidDistance, // centroid distance can range from this minimum and up to MaxCentroidDistance when attempting to limit number of haplotype groupings to MaxNumHaplotypeGroups
	int32_t MaxCentroidDistance, // maximum centroid distance which can be used when attempting to meet MaxNumHaplotypeGroups constraint
	int32_t MaxNumHaplotypeGroups) // attempt to constrain number of haplotype groups to be this maximum

{
tsHGBinSpec Tmp;
Tmp.ChromID = ChromID;
Tmp.MaxCentroidDistance = MaxCentroidDistance;
Tmp.NumLoci = NumLoci;
Tmp.StartLoci = StartLoci;
Tmp.MinCentroidDistance = MinCentroidDistance;
Tmp.MaxCentroidDistance = MaxCentroidDistance;
Tmp.MaxNumHaplotypeGroups = MaxNumHaplotypeGroups;
Tmp.Status = 0;
Tmp.pHaplotypeGroup = NULL;
return(AddHGBinSpec(&Tmp));
}

int32_t									// returned index+1 into m_pHGBinSpecs[] to allocated and initialised allele stack, 0 if errors
CCallHaplotypes::AddHGBinSpec(tsHGBinSpec* pInitHGBinSpec)	// allocated tsHGBinSpec to be initialised with a copy of pInitHGBinSpec
{
uint32_t ToAllocdHGBinSpecs;
tsHGBinSpec* pHGBinSpec;
size_t memreq;
AcquireFastSerialise();
if(m_pHGBinSpecs == NULL)					// may be NULL first time in
	{
	memreq = (size_t)cAllocHGBinSpecs * sizeof(tsHGBinSpec);
#ifdef _WIN32
	m_pHGBinSpecs = (tsHGBinSpec*)malloc((size_t)memreq);
	if(m_pHGBinSpecs == NULL)
		{
		ReleaseFastSerialise();
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddHGBinSpec: Memory allocation of %lld bytes failed", (int64_t)memreq);
		return(0);
		}
#else
	m_pHGBinSpecs = (tsHGBinSpec*)mmap(NULL, (size_t)memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	if(m_pHGBinSpecs == MAP_FAILED)
		{
		m_pHGBinSpecs = NULL;
		ReleaseFastSerialise();
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddHGBinSpec: Memory allocation of %lld bytes through mmap()  failed", (int64_t)memreq, strerror(errno));
		return(0);
		}
#endif
	m_AllocdHGBinSpecsMem = memreq;
	m_AllocdHGBinSpecs = cAllocHGBinSpecs;
	m_UsedHGBinSpecs = 0;
	}
else
		// needing to allocate more memory?
	if((m_UsedHGBinSpecs) >= m_AllocdHGBinSpecs)
		{
		ToAllocdHGBinSpecs = m_UsedHGBinSpecs + cReallocHGBinSpecs;
		size_t memreq = (size_t)ToAllocdHGBinSpecs * sizeof(tsHGBinSpec);
#ifdef _WIN32
		pHGBinSpec = (tsHGBinSpec*)realloc(m_pHGBinSpecs, memreq);
		if(pHGBinSpec == NULL)
			{
#else
		pHGBinSpec = (tsHGBinSpec*)mremap(m_pHGBinSpecs, m_AllocdHGBinSpecsMem, memreq, MREMAP_MAYMOVE);
		if(pHGBinSpec == MAP_FAILED)
			{
#endif
			ReleaseFastSerialise();
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddHGBinSpec: Memory reallocation to %lld bytes failed - %s", memreq, strerror(errno));
			return(0);
			}
		m_pHGBinSpecs = pHGBinSpec;
		m_AllocdHGBinSpecsMem = memreq;
		m_AllocdHGBinSpecs = ToAllocdHGBinSpecs;
		}
pHGBinSpec = &m_pHGBinSpecs[m_UsedHGBinSpecs++];
if(pInitHGBinSpec != NULL)
	*pHGBinSpec = *pInitHGBinSpec;
else
	memset(pHGBinSpec, 0, sizeof(tsHGBinSpec));
pHGBinSpec->BinID = m_UsedHGBinSpecs;
ReleaseFastSerialise();
return(m_UsedHGBinSpecs);
}

tsHGBinSpec*    // returned haplotype grouping bin specification, returns NULL if all bins have been iterated
CCallHaplotypes::IterateHGBinSpecs(int32_t BinID)    // previously iterated bin identifier, to access 1st bin then pass 0 as the bin identifier 
{
tsHGBinSpec* pHGBinSpec;

if(BinID < 0) // 1st bin is being requested
    BinID = 0;
AcquireFastSerialise();
while(BinID < m_UsedHGBinSpecs)
	{
	pHGBinSpec = &m_pHGBinSpecs[BinID++];
	if(!(pHGBinSpec->Status & 0x01))
		{
		pHGBinSpec->Status = 0x01;
		m_ProccessingBinID = pHGBinSpec->BinID;
		ReleaseFastSerialise();
		return(pHGBinSpec);
		}
	}
ReleaseFastSerialise();
return(NULL);
}


	// load and parse haplotype grouping bin specification file
    // This file specifies grouping bins for specific regions on chromosomes starting at requested loci of specified length. The number of bins in that region is specified and required maximal clustering distance
    // Thus regions of interest can easily be covered without having to specify explicitly each bin in that region
int								 // eBSFSuccess or error
CCallHaplotypes::LoadHHGBinSpecs(char* pszInFile)  // load grouping specifications from this file

{
int Rslt;
int32_t CurLineNumber;
int32_t NumFields;
int32_t RegionStartLoci;
int32_t NumBins;
int32_t BinSize;
int32_t RegionLen;
int32_t MinCentroidDistance; // centroid distance can range from this minimum and up to MaxCentroidDistance when attempting to limit number of haplotype groupings to MaxNumHaplotypeGroups
int32_t MaxCentroidDistance; // maximum centroid distance which can be used when attempting to meet MaxNumHaplotypeGroups constraint
int32_t MaxNumHaplotypeGroups;// attempt to constrain number of haplotype groups to be this maximum
char *pszChrom;
int32_t ChromID;
int32_t CurChromID;
tsChromMetadata* pChromMetadata;
int32_t NumUnrecognisedChroms;
int32_t NumRegionErrs;
int32_t NumRegions;
int32_t NumHGBinSpecs;

if(m_pCSVFile != NULL) // should be, but better to be sure!
	{
	delete m_pCSVFile;
	m_pCSVFile = NULL;
	}
if((m_pCSVFile = new CCSVFile) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to instantiate CCSVfile");
	Reset();
	return(eBSFerrObj);
	}
m_pCSVFile->SetMaxFields(10); // only actually processing 1st 7 fields: chrom,region start loci, bin size, number of bins, MinCentroidDistance, MaxCentroidDistance,MaxNumHaplotypeGroups

if((Rslt = m_pCSVFile->Open(pszInFile)) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to open file: '%s'", pszInFile);
	Reset();
	return(Rslt);
	}

NumHGBinSpecs = 0;
CurLineNumber = 0;
CurChromID = 0;
pChromMetadata = NULL;
NumRegions = 0;
NumUnrecognisedChroms = 0;
NumRegionErrs = 0;

while((Rslt = m_pCSVFile->NextLine()) > 0)		// onto next line containing fields
	{
	CurLineNumber++;
	if((NumFields = m_pCSVFile->GetCurFields()) < 7)	// must contain at least 7 fields
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Haplotype group bin specification file '%s' expected to contain a minimum of 7 fields, it contains %d at line %d", pszInFile, NumFields, CurLineNumber);
		Reset();
		return(eBSFerrParse);
		}

	// 1st row may be a header row, if so then skip this row
	if(CurLineNumber == 1 && m_pCSVFile->IsLikelyHeaderLine())
		continue;
	m_pCSVFile->GetText(1, &pszChrom);



	if(pszChrom == NULL || pszChrom[0] == '\0')
		{
		if(++NumUnrecognisedChroms <= 100)
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "Unrecognised chromosome/sequence name at line %d in file '%s', region sloughed", CurLineNumber, pszInFile);
		continue;
		}

	if(!AcceptThisChromName(pszChrom)) // may not be accepting this chromosome for processing
		continue;

	// is chromosome known?
	if((ChromID = LocateChrom(pszChrom)) <= 0)
		{
		if(++NumUnrecognisedChroms <= 100)
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "Unrecognised chromosome/sequence name '%s' at line %d in file '%s', region sloughed", pszChrom, CurLineNumber, pszInFile);
		continue;
		}

	// chromosome known so can determine it's size
	if(CurChromID != ChromID)
		{
		pChromMetadata = LocateChromMetadataFor(1, ChromID);
		CurChromID = ChromID;
		}
	if(pChromMetadata->ChromLen < 10)
		{
		if(++NumUnrecognisedChroms <= 100)
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "Chromosome '%s' with length '%d' is smaller than minimum sized bin of 10 at line %d in file '%s', region sloughed", pszChrom, pChromMetadata->ChromLen, CurLineNumber, pszInFile);
		continue;
		}

	m_pCSVFile->GetInt(2, &RegionStartLoci);
	if(RegionStartLoci < 0 || (RegionStartLoci + 10 > pChromMetadata->ChromLen))
		{
		if(++NumRegionErrs <= 100)
			{
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "Chromosome/sequence '%s' has size: %d, referenced in file '%s' at line: %d", pszChrom, pChromMetadata->ChromLen, pszInFile, CurLineNumber);
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "Region start loci '%d' must be in range 0..'%d' (minimum bin size is 10), region sloughed", RegionStartLoci, pChromMetadata->ChromLen - 10);
			}
		continue;
		}

	m_pCSVFile->GetInt(3, &BinSize);
	if(BinSize < 10 || (BinSize + RegionStartLoci > pChromMetadata->ChromLen))    // lower limit of 10 for bin sizes
		{
		if(++NumRegionErrs <= 100)
			{
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "Chromosome/sequence '%s' has size: %d, referenced in file '%s' at line: %d", pszChrom, pChromMetadata->ChromLen, pszInFile, CurLineNumber);
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "Bin size '%d' starting from %d is less than 10 or would extend past end of chromosome, region sloughed", BinSize, RegionStartLoci);
			}
		continue;
		}

	m_pCSVFile->GetInt(4, &NumBins);
	if(NumBins < 1)
		{
		if(++NumRegionErrs <= 100)
			{
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "Chromosome/sequence '%s' has size: %d, referenced in file '%s' at line: %d", pszChrom, pChromMetadata->ChromLen, pszInFile, CurLineNumber);
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "NumBins '%d' must be at least 1, region sloughed", NumBins);
			}
		continue;
		}
	RegionLen = NumBins * BinSize;
	if(RegionLen + RegionStartLoci > pChromMetadata->ChromLen)
		{
		if(++NumRegionErrs <= 100)
			{
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "Chromosome/sequence '%s' has size: %d, referenced in file '%s' at line: %d", pszChrom, pChromMetadata->ChromLen, pszInFile, CurLineNumber);
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "Region size %d (%d bins * %d binsize) + region start loci %d  would extend past end of chromosome, sloughed'", RegionLen, NumBins, BinSize,RegionStartLoci);
			}
		continue;
		}

	m_pCSVFile->GetInt(5, &MinCentroidDistance);
	if(MinCentroidDistance < 1)
		{
		if(++NumRegionErrs <= 100)
			{
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "Chromosome/sequence '%s' has size: %d, referenced in file '%s' at line: %d", pszChrom, pChromMetadata->ChromLen, pszInFile, CurLineNumber);
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "MinCentroidDistance '%d' must be at least 1, region sloughed", MinCentroidDistance);
			}
		continue;
		}

	if(MinCentroidDistance >= BinSize)
		{
		if(++NumRegionErrs <= 100)
			{
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "Chromosome/sequence '%s' has size: %d, referenced in file '%s' at line: %d", pszChrom, pChromMetadata->ChromLen, pszInFile, CurLineNumber);
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "Group clustering MinCentroidDistance '%d' must be less than the binsize '%d', region sloughed", MinCentroidDistance, BinSize);
			}
		continue;
		}


	m_pCSVFile->GetInt(6, &MaxCentroidDistance);
	if(MaxCentroidDistance < MinCentroidDistance)
		{
		if(++NumRegionErrs <= 100)
			{
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "Chromosome/sequence '%s' has size: %d, referenced in file '%s' at line: %d", pszChrom, pChromMetadata->ChromLen, pszInFile, CurLineNumber);
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "MaxCentroidDistance '%d' must be at least MinCentroidDistance '%d', region sloughed", MaxCentroidDistance,MinCentroidDistance);
			}
		continue;
		}

	if(MaxCentroidDistance >= BinSize)
		{
		if(++NumRegionErrs <= 100)
			{
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "Chromosome/sequence '%s' has size: %d, referenced in file '%s' at line: %d", pszChrom, pChromMetadata->ChromLen, pszInFile, CurLineNumber);
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "Group clustering MaxCentroidDistance '%d' must be less than the binsize '%d', region sloughed", MaxCentroidDistance, BinSize);
			}
		continue;
		}

	m_pCSVFile->GetInt(7, &MaxNumHaplotypeGroups);
	if(MaxNumHaplotypeGroups < 1)
		{
		if(++NumRegionErrs <= 100)
			{
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "Chromosome/sequence '%s' has size: %d, referenced in file '%s' at line: %d", pszChrom, pChromMetadata->ChromLen, pszInFile, CurLineNumber);
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "MaxNumHaplotypeGroups '%d' must be at least 1, region sloughed", MaxNumHaplotypeGroups);
			}
		continue;
		}

	if(MaxNumHaplotypeGroups > m_NumFounders) // silently clamp
		MaxNumHaplotypeGroups = m_NumFounders;

	while(NumBins--)
		{
		if(AddHGBinSpec(ChromID, RegionStartLoci, BinSize, MinCentroidDistance,MaxCentroidDistance,MaxNumHaplotypeGroups) < 1)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddHGBinSpec(%d,%d,%d,%d,%d,%d) returned unrecoverable memory allocation error,parsing file '%s' at line %d", ChromID, RegionStartLoci, BinSize, MinCentroidDistance,MaxCentroidDistance,MaxNumHaplotypeGroups, pszInFile, CurLineNumber);
			Reset();
			return(eBSFerrParse);
			}
		RegionStartLoci += BinSize;
		NumHGBinSpecs++;
		}
	NumRegions++;
	}
delete m_pCSVFile;
m_pCSVFile = NULL;
gDiagnostics.DiagOut(eDLInfo, gszProcName, "Parsed and accepted '%d' haplotype group region specifications resulting in '%d' bins for processing, sloughed '%d' regions, in file '%s', total accumulated bins '%d'",NumRegions, NumHGBinSpecs,NumUnrecognisedChroms+NumRegionErrs, pszInFile, m_UsedHGBinSpecs);
return(m_UsedHGBinSpecs);
}


uint32_t									// returned index+1 into m_pAlleleStacks[] to allocated and initialised allele stack, 0 if errors
CCallHaplotypes::AddAlleleStack(tsAlleleStack* pInitAlleleStack)	// allocated tsAlleleStack to be initialised with a copy of pInitAlleleStack
{
uint32_t ToAllocdAlleleStacks;
tsAlleleStack* pAlleleStack;
size_t memreq;
AcquireFastSerialise();
if(m_pAlleleStacks == NULL)					// may be NULL first time in
	{
	memreq = (size_t)cAllocAlleleStacks * sizeof(tsAlleleStack);
#ifdef _WIN32
	m_pAlleleStacks = (tsAlleleStack*)malloc((size_t)memreq);
	if(m_pAlleleStacks == NULL)
		{
		ReleaseFastSerialise();
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddAlleleStack: Memory allocation of %lld bytes failed", (int64_t)memreq);
		return(0);
		}
#else
	m_pAlleleStacks = (tsAlleleStack*)mmap(NULL, (size_t)memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	if(m_pAlleleStacks == MAP_FAILED)
		{
		m_pAlleleStacks = NULL;
		ReleaseFastSerialise();
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddAlleleStack: Memory allocation of %lld bytes through mmap()  failed", (int64_t)memreq, strerror(errno));
		return(0);
		}
#endif
	m_AllocdAlleleStacksMem = memreq;
	m_AllocdAlleleStacks = cAllocAlleleStacks;
	m_UsedAlleleStacks = 0;
	}
else
		// needing to allocate more memory?
	if((m_UsedAlleleStacks) >= m_AllocdAlleleStacks)
		{
		ToAllocdAlleleStacks = m_UsedAlleleStacks + cReallocAlleleStacks;
		size_t memreq = (size_t)ToAllocdAlleleStacks * sizeof(tsAlleleStack);
#ifdef _WIN32
		pAlleleStack = (tsAlleleStack*)realloc(m_pAlleleStacks, memreq);
		if(pAlleleStack == NULL)
			{
#else
		pAlleleStack = (tsAlleleStack*)mremap(m_pAlleleStacks, m_AllocdAlleleStacksMem, memreq, MREMAP_MAYMOVE);
		if(pAlleleStack == MAP_FAILED)
			{
#endif
			ReleaseFastSerialise();
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddAlleleStack: Memory reallocation to %lld bytes failed - %s", memreq, strerror(errno));
			return(0);
			}
		m_pAlleleStacks = pAlleleStack;
		m_AllocdAlleleStacksMem = memreq;
		m_AllocdAlleleStacks = ToAllocdAlleleStacks;
		}
pAlleleStack = &m_pAlleleStacks[m_UsedAlleleStacks++];
if(pInitAlleleStack != NULL)
	*pAlleleStack = *pInitAlleleStack;
else
	memset(pAlleleStack, 0, sizeof(tsAlleleStack));
ReleaseFastSerialise();
return(m_UsedAlleleStacks);
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
									int32_t NumChroms)			// processing this number of chromosomes
{
uint32_t MaxWait;
uint32_t StartQuerySeqIdx;
uint32_t ThreadIdx;
uint32_t StartedInstances;
tsWorkerInstance *pThreadPar;

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
	Sleep(5000);
#else
	sleep(5);
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
int CurBinID;
int CurChromID;
uint32_t NumQueueElsProcessed;
uint32_t MarkNumQueueElsProcessed;
tsWorkQueueEl *pWorkQueueEl;
uint32_t ReqTerminate;
// this thread has started, one more worker thread
AcquireSerialise();
m_NumWorkerInsts++;
ReleaseSerialise();
CurChromID = -1;
NumQueueElsProcessed = 0;
pWorkQueueEl = NULL;
CurBinID = 0;
Rslt = eBSFSuccess;
while(Rslt >= eBSFSuccess) {
	// check if requested to terminate
	AcquireSerialise();
	ReqTerminate = m_ReqTerminate;
	ReleaseSerialise();
	if(ReqTerminate)
		break;

	if(!(m_PMode == eMCSHAllelicHapsGrps || m_PMode == eMCSHCoverageHapsGrps))
		{
		AcquireSerialise();
    	NumQueueElsProcessed = 	m_NumQueueElsProcessed;
	    if(NumQueueElsProcessed == m_TotWorkQueueEls)
		    {
			ReleaseSerialise();
			break;
			}
		m_NumQueueElsProcessed++;
		ReleaseSerialise();
		pWorkQueueEl = &m_pWorkQueueEls[NumQueueElsProcessed];
		Rslt = GenChromAlleleStacks(pWorkQueueEl->ChromID, pWorkQueueEl->ChromLen, pWorkQueueEl->StartLoci, pWorkQueueEl->MaxNumLoci, pWorkQueueEl->NumFndrs, pWorkQueueEl->pFounderPBAs, pWorkQueueEl->pMskPBA);
		}
	else
		{
		tsHGBinSpec* pCurBinSpec;
		if((pCurBinSpec = IterateHGBinSpecs(CurBinID)) == NULL)
			break;
		CurBinID = pCurBinSpec->BinID;

		// locate work queue element corresponding to be processed bin ChromID
		if(pCurBinSpec->ChromID != CurChromID)
			{
			MarkNumQueueElsProcessed = NumQueueElsProcessed; // so can restore if bin chromosome not present in all samples
			pWorkQueueEl = &m_pWorkQueueEls[NumQueueElsProcessed];
			while(NumQueueElsProcessed < m_TotWorkQueueEls && pWorkQueueEl->ChromID != pCurBinSpec->ChromID)
				pWorkQueueEl = &m_pWorkQueueEls[++NumQueueElsProcessed];
			if(NumQueueElsProcessed == m_TotWorkQueueEls)
				{
				pCurBinSpec->Status |= 0x07; // mark as being completed but able to locate chromosome - user specification was for a chromosome not present in all samples
				NumQueueElsProcessed = MarkNumQueueElsProcessed; // attempt to iterate forward onto a bin which does reference a chromsome present in all samples
				continue;
				}
			CurChromID = pCurBinSpec->ChromID;
			}
		Rslt = GenHaplotypeGroups(pCurBinSpec,pWorkQueueEl->NumFndrs, pWorkQueueEl->pFounderPBAs);
		pCurBinSpec->Status = 0x03; // mark processing on this bin as completed
		}
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


int				// < 0 if errors, otherwise success
CCallHaplotypes::GenChromAlleleStacks(int32_t ChromID,	// processing is for this chromosome
										  int32_t ChromSize,		// chromosome is this size
										  int32_t Loci,			// processing for allele stacks starting from this loci
										  int32_t MaxNumLoci,		// processing for this maximum number of loci
										  int32_t NumFndrs,		// number of founders to be processed 1st PBA at *pPBAs[0]
										  uint8_t* pFounderPBAs[],	// pPBAs[] pts to chromosome PBAs, for each of the chromosome founder PBAs
										  uint8_t* pMskPBA)			// pts to optional masking PBA, scoring only for segments contained in mask which are non-zero and where progeny and founder PBAs also have an allele.
																	// enables a founder to be processed as if a progeny, restricting scoring Kmers to be same as if a progeny

{
int32_t FounderIdx;
tsAlleleStack AlleleStack;
uint32_t NumAlleleStacks;
bool bFndrAllele;
uint32_t PrevAlleleIdx;
int32_t NumFndrs2Proc;
uint8_t* pFndrLoci;
uint32_t AlleleIdx;
uint8_t AlleleMsk;
int32_t AlleleLoci;
int32_t FndrIdx;
int32_t DivAlleles;
int32_t LociFounders;
uint8_t FndrBaseAllele;
ts4096Bits Fndrs2Proc;

if(NumFndrs < 1 || NumFndrs > 4096)
	return(-1);

// check that at least 1 founder actually has a PBAs
memset(&Fndrs2Proc, 0, sizeof(Fndrs2Proc));
NumFndrs2Proc = 0;
for(FounderIdx = 0; FounderIdx < NumFndrs; FounderIdx++)
	{
	if(!(m_Fndrs2Proc[FounderIdx] & 0x01))		// not interested as founder not marked for scoring?
		continue;
	if((pFndrLoci = pFounderPBAs[FounderIdx]) == NULL)	// can't score if no PBA for this founder
		continue;
	BitsVectSet(NumFndrs2Proc, Fndrs2Proc);
	NumFndrs2Proc++;							// can use this founder for generating an allelestack
	}
if(NumFndrs2Proc < 1)							// must have at least 1 founder
	return(-2);

int32_t EndLoci = min(Loci + MaxNumLoci, ChromSize) - 1;

memset(&AlleleStack, 0, sizeof(AlleleStack));
AlleleStack.ChromID = ChromID;
AlleleStack.NumFndrs = NumFndrs;
AlleleStack.NumProcFndrs = NumFndrs2Proc;
memcpy(&AlleleStack.ProcFndrs, &Fndrs2Proc, sizeof(ts4096Bits));

bFndrAllele = false;
NumAlleleStacks = 0;
for(AlleleLoci = Loci; AlleleLoci <= EndLoci; AlleleLoci++)
	{
	if(pMskPBA != NULL && pMskPBA[AlleleLoci] == 0)	// skipping over any loci in control which has no PBA
		continue;

// build stack of alleles containing founder alleles at AlleleLoci
// to be accepted as a founder allele stack then each founder must have a dirac allele unique to that founder only
// progeny can then use founder allele stacks to inform as to the lineage  
	AlleleStack.Loci = AlleleLoci;
	memset(AlleleStack.NumAlleleFndrs, 0, sizeof(AlleleStack.NumAlleleFndrs));
	memset(AlleleStack.Alleles, 0, sizeof(AlleleStack.Alleles));
	LociFounders = 0;
	PrevAlleleIdx = 0;
	for(FndrIdx = 0; FndrIdx < NumFndrs; FndrIdx++)
		{
		if(!(m_Fndrs2Proc[FndrIdx] & 0x01))	// skip this founder?
			continue;
		bFndrAllele = false;
		pFndrLoci = pFounderPBAs[FndrIdx] + AlleleLoci;
		if(*pFndrLoci == 0)					
			{
			if(m_bAllFndrsLociAligned)           // must all all founders must have alignment alleles at the AlleleLoci
				{
				bFndrAllele = false;
				break;							// founder has no allele - can't have had an alignment to reference at AlleleLoci
				}
			continue;                          // next founder may have an alignment
			}

		AlleleMsk = 0x03;					// differentiation is at the allele level
		for(AlleleIdx = 0; AlleleIdx < 4; AlleleIdx++, AlleleMsk <<= 2)
			{
			FndrBaseAllele = (*pFndrLoci & AlleleMsk) >> (AlleleIdx * 2);

			if(FndrBaseAllele <= 2)		// no allele, with founders then only interest is in diracs
				continue;				// try next allele

			if(bFndrAllele)				// founders must have a single unique allele
				{
				bFndrAllele = false;
				break;
				}
			bFndrAllele = true;			// accepting this founder as having at least one allele, if more alleles then founder loci treated as having no alleles
			AlleleStack.NumAlleleFndrs[AlleleIdx]++;
			BitsVectSet(FndrIdx, AlleleStack.Alleles[AlleleIdx]);
			}
		if(!bFndrAllele)				// if no alleles of required level for this founder then don't bother with other founders
			break;
		}

	if(!bFndrAllele)				// if no alleles then process next loci
		continue;

		// is there diversity in alleles between all founders - need diversity in order to differentiate founder alleles in progeny
	for(DivAlleles = AlleleIdx = 0; AlleleIdx < 4; AlleleIdx++)
		{
		if(AlleleStack.NumAlleleFndrs[AlleleIdx] == 1)
			DivAlleles++;
		}

	if(DivAlleles < 1) // > NumFndrs2Proc)						// must be at least one founder unique relative to others to be informative...
		continue;

	AddAlleleStack(&AlleleStack);
	NumAlleleStacks++;
	}
return(NumAlleleStacks);
}


int				// < 0 if errors, >=0 length of generated consensus
CCallHaplotypes::GenConsensusPBA(int32_t Loci,			// processing for consensus starting from this loci
			int32_t NumLoci,		// processing for generation of this maximum sized consensus PBAs
			int32_t NumFndrs,		// number of founders to be processed 1st PBA at *pPBAs[0]
			uint8_t* pFounderPBAs[], // pPBAs[] pts to chromosome PBAs, for each of the chromosome founder PBAs
			uint8_t* pConsensusPBAs)     // ptr to preallocated sequence of length MaxNumLoci which is to be updated with consensus PBAs

{
uint32_t AlleleFreq[256];
uint8_t* pFndrLoci;
int32_t AlleleLoci;
int32_t FndrIdx;
int32_t EndLoci;
int32_t ConsensusLen;
uint8_t MaxFeqAllele;

EndLoci = Loci + NumLoci;

ConsensusLen = 0;
for(AlleleLoci = Loci; AlleleLoci < EndLoci; AlleleLoci++,ConsensusLen++)
	{
	memset(AlleleFreq,0,sizeof(AlleleFreq));
	MaxFeqAllele = 0;
	for(FndrIdx = 0; FndrIdx < NumFndrs; FndrIdx++)
		{
		pFndrLoci = pFounderPBAs[FndrIdx] + AlleleLoci;
		AlleleFreq[*pFndrLoci]++;
		if(AlleleFreq[*pFndrLoci] > AlleleFreq[MaxFeqAllele])
			MaxFeqAllele = *pFndrLoci;
		}
	*pConsensusPBAs++ = MaxFeqAllele;
	}
return(ConsensusLen);
}

uint8_t				// concensus PBA for founders in a group of which specified founder is a member at the requested loci
CCallHaplotypes::GenConsensusGroupPBA(int Founder, // consensus for all founders in same group as this founder
	tsHaplotypeGroup *pHaplotypes, // pts to haplotype grouping expected to contain requested chrom and loci
	        int32_t ChromID,        // requiring consensus for this chrom at Loci
	        int32_t Loci,			// requiring consensus at this loci
			int32_t NumFndrs,		// number of founders to be processed 1st PBA at *pPBAs[0]
			uint8_t* pFounderPBAs[]) // pPBAs[] pts to chromosome PBAs, for each of the chromosome founder PBAs
{
int32_t HapGrp;
ts4096Bits* pHapGrp;
uint32_t AlleleFreq[256];
uint8_t* pFndrLoci;
int32_t FndrIdx;
uint8_t MaxFeqAllele;

if(pHaplotypes == NULL)
    return(0);

if(ChromID != pHaplotypes->ChromID || (Loci < pHaplotypes->StartLoci || (Loci >= pHaplotypes->StartLoci + pHaplotypes->NumLoci)))
  return(0);

if(pHaplotypes->MRAGrpChromID != ChromID || pHaplotypes->MRAGrpLoci != Loci) // different loci from that cached? If so then need to clear cache so consensus entries are recalculated
	{
	pHaplotypes->MRAGrpChromID = 0;
	pHaplotypes->MRAGrpLoci = -1;
	memset(pHaplotypes->MRAGrpConsensus,0x0ff,sizeof(pHaplotypes->MRAGrpConsensus));
	}

// groupings in which loci is located now known
// identify in which group the founder is a member
pHapGrp = pHaplotypes->HaplotypeGroup;
for(HapGrp = 0; HapGrp < pHaplotypes->NumHaplotypeGroups; HapGrp++,pHapGrp++)
	{
	if(BitsVectTest(Founder, *pHapGrp))
		break;
	}
if(HapGrp == pHaplotypes->NumHaplotypeGroups) // should have been a member of a group, but check!
	return(0);

if(pHaplotypes->MRAGrpChromID == ChromID && pHaplotypes->MRAGrpLoci == Loci && pHaplotypes->MRAGrpConsensus[HapGrp] != 0x0ff) // with any luck the consensus for group has already been generated and cached
    return(pHaplotypes->MRAGrpConsensus[HapGrp]); // returning the cached consensus

// group has been identified, generate and cache the consensus, consensus will be highest frequency allele within that groups members
pHaplotypes->MRAGrpChromID = ChromID;
pHaplotypes->MRAGrpLoci = Loci;
memset(AlleleFreq,0,sizeof(AlleleFreq));
MaxFeqAllele = 0;
for(FndrIdx = 0; FndrIdx < NumFndrs; FndrIdx++)
	{
	if(!BitsVectTest(FndrIdx, *pHapGrp))
		continue;
	pFndrLoci = pFounderPBAs[FndrIdx] + Loci;
	AlleleFreq[*pFndrLoci]++;
	if(AlleleFreq[*pFndrLoci] > AlleleFreq[MaxFeqAllele])
		MaxFeqAllele = *pFndrLoci;
	}
pHaplotypes->MRAGrpConsensus[HapGrp] = MaxFeqAllele; // cache group concensus
return(MaxFeqAllele);
}


int32_t 
CCallHaplotypes::ReportHaplotypeGroups(int32_t NumHaplotypes)         // number of haplotypes which have been grouped
{
int32_t BinIdx;
tsHGBinSpec* pHGBinSpec;
tsHaplotypeGroup* pHaplotypeGroup;
int32_t CurChromID;
int32_t HaplotypeIdx;
int32_t GrpIdx;
char* pszChrom;
char szOutFile[_MAX_PATH];

int32_t NumGrpMembers[cMaxFounderReadsets];

AcquireSerialise();
if(m_pszOutBuffer == NULL)
	{
	if((m_pszOutBuffer = new uint8_t[cOutBuffSize]) == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "ReportHaplotypeGroups: Memory allocation of %u bytes for output buffering failed - %s", cOutBuffSize, strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
	m_AllocOutBuff = cOutBuffSize;
	}
m_OutBuffIdx = 0;

sprintf(szOutFile, "%s.haplotypegrouping.%d.%d.all.csv", m_pszRsltsFileBaseName, m_ExprID,NumHaplotypes);
gDiagnostics.DiagOut(eDLInfo, gszProcName, "Reporting haplotype grouping to file '%s'", szOutFile);
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
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "ReportHaplotypeGroups: Unable to create/truncate %s - %s", szOutFile, strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}
m_OutBuffIdx = sprintf((char*)m_pszOutBuffer, "\"Chrom\",\"Loci\",\"Len\",\"MinDistance\",\"MaxDistance\",\"MaxHaplotypeGroups\",\"ActualDistance\",\"ActualHaplotypeGroups\"");
for(HaplotypeIdx = 0; HaplotypeIdx < NumHaplotypes; HaplotypeIdx++)
	m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], ",\"%s\"", LocateReadset(HaplotypeIdx+1));
for(HaplotypeIdx = 0; HaplotypeIdx < NumHaplotypes; HaplotypeIdx++)
	m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], ",\"GrpMembers:%d\"", HaplotypeIdx+1);
m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], "\n");

CurChromID = 0;
pHGBinSpec = m_pHGBinSpecs;
for(BinIdx = 0; BinIdx < m_UsedHGBinSpecs; BinIdx++,pHGBinSpec++)
	{
	if((pHGBinSpec->Status & 0x07) != 0x03) // skip any with error status or incomplete processing
		continue;

	if(pHGBinSpec->ChromID != CurChromID)
		{
		CurChromID = pHGBinSpec->ChromID;
		pszChrom = LocateChrom(CurChromID);
		}

	pHaplotypeGroup = pHGBinSpec->pHaplotypeGroup;
	m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], "\"%s\",%d,%d,%d,%d,%d,%d,%d", pszChrom, pHGBinSpec->StartLoci, pHGBinSpec->NumLoci, pHGBinSpec->MinCentroidDistance,pHGBinSpec->MaxCentroidDistance,pHGBinSpec->MaxNumHaplotypeGroups, pHaplotypeGroup->CentroidDistance,pHaplotypeGroup->NumHaplotypeGroups);
	memset(NumGrpMembers, 0, sizeof(NumGrpMembers));
	for(HaplotypeIdx = 0; HaplotypeIdx < NumHaplotypes; HaplotypeIdx++)
		{
		for(GrpIdx = 0; GrpIdx < pHaplotypeGroup->NumHaplotypeGroups; GrpIdx++)
			{
			if(BitsVectTest(HaplotypeIdx, pHaplotypeGroup->HaplotypeGroup[GrpIdx]))
				{
				NumGrpMembers[GrpIdx] += 1;
				m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], ",%d", GrpIdx + 1);
				break;
				}
			}
		}
	for(GrpIdx = 0; GrpIdx < NumHaplotypes; GrpIdx++)
		m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], ",%d", NumGrpMembers[GrpIdx]);

	m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], "\n");
	if((m_OutBuffIdx + 1000) > m_AllocOutBuff)
		{
		if(!CUtility::RetryWrites(m_hOutFile, m_pszOutBuffer, m_OutBuffIdx))
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "ReportHaplotypeGroups: Fatal error in RetryWrites()");
			Reset();
			return(eBSFerrFileAccess);
			}
		m_OutBuffIdx = 0;
		}
	}

if(m_hOutFile != -1)
	{
	if(m_OutBuffIdx)
		{
		if(!CUtility::RetryWrites(m_hOutFile, m_pszOutBuffer, m_OutBuffIdx))
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "ReportAnchorsAsCSV: Fatal error in RetryWrites()");
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
	}
gDiagnostics.DiagOut(eDLInfo, gszProcName, "Completed reporting haplotype grouping to file '%s'", szOutFile);
ReleaseSerialise();
return(eBSFSuccess);
}


tsHaplotypeGroup * // returns allocated ptr to haplotype groups - allocated with 'new', free with 'delete'
CCallHaplotypes::GroupHaplotypes(int32_t ChromID, // grouping is on this chromosome
			          int32_t StartLoci, // grouping starts at this loci
	                     int32_t Length,     // grouping is over this many loci
			             int32_t MinCentroidDistance, // centroid distance can range from this minimum and up to MaxCentroidDistance when attempting to limit number of haplotype groupings to MaxNumHaplotypeGroups
			             int32_t MaxCentroidDistance, // maximum centroid distance which can be used when attempting to meet MaxNumHaplotypeGroups constraint
			             int32_t MaxNumHaplotypeGroups, // attempt to constrain number of haplotype groups to be this maximum
			             uint32_t* pFndrDiffs, // matrix counts of founder differentials
			             int32_t NumFndrs) // number of founders
{
ts4096Bits *pClusterMembers;
ts4096Bits *pMinimisedClusterMembers;
ts4096Bits AssocToAGroup;
int32_t nCol;

uint32_t* pFndrDiff;

int32_t nRow;

int32_t NumHaplotypeGroups;
ts4096Bits CurClusterFounders;
int32_t CurClusterSize;
int32_t TargCentroidDistance;
int32_t TargMaxCentroidDistance;
int32_t TargMinCentroidDistance;
int32_t MinimisedCentroidDistance;
int32_t CurClusterDiffs;
ts4096Bits SelClusterFounders;
int32_t SelClusterSize;
int32_t SelClusterDiffs;
int32_t NumFndrsChkd;

if((pClusterMembers = new ts4096Bits[NumFndrs]) == NULL)
    return(NULL);
if((pMinimisedClusterMembers = new ts4096Bits[NumFndrs]) == NULL)
	{
	delete []pClusterMembers;
	return(NULL);
	}

TargMaxCentroidDistance = MaxCentroidDistance;
TargMinCentroidDistance = MinCentroidDistance;
MinimisedCentroidDistance = -1;
do {
	TargCentroidDistance = (TargMaxCentroidDistance + TargMinCentroidDistance)/2;

	// for each founder determine which other founders are within a specified maximal distance and assign a founder grouping to this group
	// objective is to maximise number of founders in a grouping
	BitsVectInitialise(0, AssocToAGroup);
	BitsVectInitialise(0, SelClusterFounders);
	NumHaplotypeGroups = 0;
	do {
		pFndrDiff = pFndrDiffs;
		SelClusterSize = 0;
		SelClusterDiffs = 0;
		NumFndrsChkd = 0;
		for(nRow = 0; nRow < NumFndrs; nRow++)
			{
			CurClusterSize = 0;
			CurClusterDiffs = 0;
			BitsVectInitialise(0, CurClusterFounders);
			for(nCol = 0; nCol < NumFndrs; nCol++, pFndrDiff++)
				{
				if(BitsVectTest(nCol, AssocToAGroup))
					continue;
				NumFndrsChkd += 1;
				if(*pFndrDiff > (uint32_t)TargCentroidDistance) // if outside of group threshold then can't be member of that group
					continue;
				BitsVectSet(nCol, CurClusterFounders);
				CurClusterSize += 1;
				CurClusterDiffs += *pFndrDiff;
				}

			if(CurClusterSize > SelClusterSize ||
				(CurClusterSize == SelClusterSize && CurClusterDiffs < SelClusterDiffs))
				{
				SelClusterSize = CurClusterSize;
				SelClusterDiffs = CurClusterDiffs;
				SelClusterFounders = CurClusterFounders;
				}
			}
		if(NumFndrsChkd)
			{
			pClusterMembers[NumHaplotypeGroups++] = SelClusterFounders;
			BitsVectUnion(AssocToAGroup, SelClusterFounders);
			}
		}
	while(NumFndrsChkd);

	if(MinCentroidDistance == MaxCentroidDistance)
		break;
	
	if(NumHaplotypeGroups > MaxNumHaplotypeGroups) // if more groups than targeted then need to increase the centroid distance attempting to reduce groups
		TargMinCentroidDistance = TargCentroidDistance + 1;
	else // else if fewer or equal groups than targeted then need to minimise the centroid distance until distance boundary at which number of groups is more than targeted
		{
		if(NumHaplotypeGroups == MaxNumHaplotypeGroups)
			{
			MinimisedCentroidDistance = TargCentroidDistance;
			memcpy(pMinimisedClusterMembers, pClusterMembers, sizeof(ts4096Bits) * MaxNumHaplotypeGroups);
			}
		TargMaxCentroidDistance = TargCentroidDistance - 1;
		}
	}
while(TargMaxCentroidDistance >= TargMinCentroidDistance);

if(MinimisedCentroidDistance > -1 && MinimisedCentroidDistance != TargCentroidDistance)
	{
	TargCentroidDistance = MinimisedCentroidDistance;
	NumHaplotypeGroups = MaxNumHaplotypeGroups;
	memcpy(pClusterMembers,pMinimisedClusterMembers, sizeof(ts4096Bits) * MaxNumHaplotypeGroups);
	}

tsHaplotypeGroup* pGroupings;
size_t MemReq = sizeof(tsHaplotypeGroup) + (sizeof(ts4096Bits) * (NumHaplotypeGroups-1));
pGroupings = (tsHaplotypeGroup *) new uint8_t[MemReq];
memset(pGroupings, 0, MemReq);
pGroupings->Size = (int32_t)MemReq;
pGroupings->ChromID = ChromID;
pGroupings->StartLoci = StartLoci;
pGroupings->NumLoci = Length;
pGroupings->NumFndrs = NumFndrs;
pGroupings->CentroidDistance = TargCentroidDistance;
pGroupings->MinCentroidDistance = MinCentroidDistance;
pGroupings->MaxCentroidDistance = MaxCentroidDistance;
pGroupings->MaxNumHaplotypeGroups = MaxNumHaplotypeGroups;
pGroupings->MRAGrpChromID = -1;
pGroupings->MRAGrpLoci = -1;
pGroupings->NumHaplotypeGroups = NumHaplotypeGroups;
memcpy(pGroupings->HaplotypeGroup, pClusterMembers, sizeof(ts4096Bits) * NumHaplotypeGroups);
delete []pClusterMembers;
delete []pMinimisedClusterMembers;
return(pGroupings);
}



int				// < 0 if errors, otherwise success
CCallHaplotypes::GenHaplotypeGroups(tsHGBinSpec *pHGBinSpec,   // clustering using these specs
										  int32_t NumFndrs,		// number of founders to be processed 1st PBA at *pPBAs[0]
										  uint8_t* pFounderPBAs[])	// pPBAs[] pts to chromosome PBAs, for each of the chromosome founder PBAs
{
int32_t	GroupingPhase;  // grouping can be multiphased, in inital (1) phase there is no imputation for unaligned PBA values, in next (2) phase then unaligned PBA values are replaced with group concensus PBA values

tsHaplotypeGroup *pP1CurHapGroups;

tsHaplotypeGroup* pCurHaplotypeGroups; // pts most recently generated haplotype groups

int32_t EndLoci;
uint8_t* pFndrLoci;
uint8_t *pChkFndrLoci;

int32_t AlleleLoci;

int32_t FndrIdx;
int32_t PrevFndrIdx;
int32_t ChkFndrIdx;

uint32_t *pFndrDiff;
uint32_t* pFndrDiffs;

uint8_t FndrLociPBA;
uint8_t ChkFndrLociPBA;
uint8_t *pConsensusPBA;
uint8_t *pConsensusLociPBA;

pConsensusPBA = new uint8_t[pHGBinSpec->NumLoci];
GenConsensusPBA(pHGBinSpec->StartLoci, pHGBinSpec->NumLoci, NumFndrs, pFounderPBAs, pConsensusPBA);
pFndrDiffs = new uint32_t[NumFndrs * NumFndrs]; // organised as [row][col]

pP1CurHapGroups = NULL;
GroupingPhase = 0; // starting with phase 0 using all founders consensus for missing PBA alignment loci and then onto subsequent refining phases which instead use group membership to derive the group consensus for missing PBA alignment loci
do {
	EndLoci = pHGBinSpec->StartLoci + pHGBinSpec->NumLoci - 1;

		// generating all vs.all matrix
	memset(pFndrDiffs, 0, NumFndrs * NumFndrs * sizeof(uint32_t));
	for(AlleleLoci = pHGBinSpec->StartLoci; AlleleLoci <= EndLoci; AlleleLoci++)
		{
		PrevFndrIdx = -1;
		pFndrDiff = pFndrDiffs;
		pConsensusLociPBA = pConsensusPBA + AlleleLoci - pHGBinSpec->StartLoci;
		for(FndrIdx = 0; FndrIdx < NumFndrs; FndrIdx++)  // iterating founder by row
			{
			pFndrLoci = pFounderPBAs[FndrIdx] + AlleleLoci;

			for(ChkFndrIdx = 0; ChkFndrIdx < NumFndrs; ChkFndrIdx++, pFndrDiff++)  // iterating founder by col
				{
				pChkFndrLoci = pFounderPBAs[ChkFndrIdx] + AlleleLoci;
				switch(GroupingPhase) { // using a switch() instead of a 'if' only because there may be a need for phase specific processing in the future????
					case 0: // inital phase is using consensus over all founders if no alignment PBA
						if((FndrLociPBA = *pFndrLoci) == 0)
							FndrLociPBA = *pConsensusLociPBA;
						if((ChkFndrLociPBA = *pChkFndrLoci) == 0)
							ChkFndrLociPBA = *pConsensusLociPBA;
						break;

					default:  // subsequent phases are refining the group consensus if no alignment PBA
						if(FndrIdx != PrevFndrIdx)
							{
							PrevFndrIdx = FndrIdx;
							if((FndrLociPBA = *pFndrLoci) == 0)
								FndrLociPBA = GenConsensusGroupPBA(FndrIdx, pP1CurHapGroups, pHGBinSpec->ChromID, AlleleLoci, NumFndrs, pFounderPBAs);
							}
						if((ChkFndrLociPBA = *pChkFndrLoci) == 0)
							ChkFndrLociPBA = GenConsensusGroupPBA(ChkFndrIdx, pP1CurHapGroups, pHGBinSpec->ChromID, AlleleLoci, NumFndrs, pFounderPBAs);
						break;
						}

					if(FndrLociPBA != ChkFndrLociPBA)
						*pFndrDiff += 1;
					}
				}
			}

		// have a difference matrix, now goup haplotypes
	pCurHaplotypeGroups = GroupHaplotypes(pHGBinSpec->ChromID, pHGBinSpec->StartLoci, pHGBinSpec->NumLoci,pHGBinSpec->MinCentroidDistance, pHGBinSpec->MaxCentroidDistance,pHGBinSpec->MaxNumHaplotypeGroups, pFndrDiffs, NumFndrs);
	if(pP1CurHapGroups != NULL)
		{
		pP1CurHapGroups->MRAGrpChromID = -1;
		pP1CurHapGroups->MRAGrpLoci = -1;
		memset(pP1CurHapGroups->MRAGrpConsensus, 0, sizeof(pP1CurHapGroups->MRAGrpConsensus));
		if(pP1CurHapGroups->Size == pCurHaplotypeGroups->Size && !memcmp(pCurHaplotypeGroups, pP1CurHapGroups, pCurHaplotypeGroups->Size))
			{
			delete[]pP1CurHapGroups;
			break;
			}
		delete[]pP1CurHapGroups;
		}
	pP1CurHapGroups = pCurHaplotypeGroups;
	}
while(++GroupingPhase < m_NumHapGrpPhases);

pHGBinSpec->pHaplotypeGroup = pCurHaplotypeGroups;

if(pFndrDiffs != NULL)
    delete[]pFndrDiffs;
if(pConsensusPBA != NULL)
    delete []pConsensusPBA;
return(eBSFSuccess);
}


int				// returns number of allele stacks generated
CCallHaplotypes::AlignAlleleStacks(int32_t NumFndrs,			// number of founders to be processed against
									int32_t MaxNumLoci)		// work queue items specify this number of loci for processing by threads
{
int32_t FounderID;
tsReadsetMetadata *pReadsetMetadata;
tsChromMetadata *pChromMetadata;
tsWorkQueueEl *pWorkQueueEl;
uint32_t CurChromMetadataIdx;
uint8_t *pPBAs[cMaxFounderReadsets+1];							// additional is to allow for the progeny readset PBAs
uint8_t *pMskPBA;
int32_t ChromIdx;
int32_t NumUnits;			// NumUnits is the total number of work units to be distributed over available threads for processing chromosomes
int32_t UnitSize;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"AlignAlleleStacks: Starting to generate allele stacks over %d founders ",NumFndrs);
pReadsetMetadata = &m_Readsets[0];		// founders were loaded first so 1st founder will be here
CurChromMetadataIdx = pReadsetMetadata->StartChromMetadataIdx;
if(m_pWorkQueueEls != NULL)				// ensuring only one work queue exists
	delete []m_pWorkQueueEls;
m_pWorkQueueEls = NULL;
m_AllocWorkQueueEls = 0;
pMskPBA = NULL;

// determine maximal number of work queue elements required
CurChromMetadataIdx = pReadsetMetadata->StartChromMetadataIdx;
for(NumUnits = ChromIdx = 0; ChromIdx < pReadsetMetadata->NumChroms && CurChromMetadataIdx != 0; ChromIdx++)
	{
	pChromMetadata = &m_pChromMetadata[CurChromMetadataIdx - 1];
	if(pChromMetadata->ChromLen >= MaxNumLoci)
		NumUnits += (pChromMetadata->ChromLen + MaxNumLoci - 1) / MaxNumLoci;
	NumUnits += 1; // safety - accounting for potential partial AccumBinSizes
	CurChromMetadataIdx = pChromMetadata->NxtChromMetadataIdx;
	}
if(NumUnits == 0)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "AlignAlleleStacks: Number of thread work units is 0");
	return(0);
	}
// sizing work queue to contain a maximum of NumUnits elements
m_AllocWorkQueueEls = NumUnits;
m_pWorkQueueEls = new tsWorkQueueEl[m_AllocWorkQueueEls];
memset(m_pWorkQueueEls, 0, sizeof(tsWorkQueueEl) * m_AllocWorkQueueEls);
m_TotWorkQueueEls = 0;
m_NumQueueElsProcessed = 0;
m_FastSerialise = 0;
pWorkQueueEl = m_pWorkQueueEls;

// iterate along PBAs which are common to all founders, plus mask readset if present, and then initialise elements of work to be undertaken by each thread
CurChromMetadataIdx = pReadsetMetadata->StartChromMetadataIdx;
for(ChromIdx = 0; ChromIdx < pReadsetMetadata->NumChroms && CurChromMetadataIdx != 0; ChromIdx++)
	{
	pChromMetadata = &m_pChromMetadata[CurChromMetadataIdx - 1];

	// all founders, and mask readset if present, must have PBAs for same chromosome
	for(FounderID = 1; FounderID <= NumFndrs; FounderID++)
		{
		if((pPBAs[FounderID-1] = LocatePBAfor(FounderID, pChromMetadata->ChromID)) == NULL)
			{
			char *pszChrom = LocateChrom(pChromMetadata->ChromID);
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "AlignAlleleStacks: No PBA for chromosome '%s' in at least one founder '%s', skipping this chromosome", pszChrom, LocateReadset(FounderID));
			break;
			}
		}

	CurChromMetadataIdx = pChromMetadata->NxtChromMetadataIdx;		// ready for next iteration over the chromosomes
	if(FounderID <= NumFndrs)	// true only if PBAs not common to progeny and all founders, including mask readset
		continue;

	if(m_MaskReadsetID == 0)				// if masking then chromosome PBA for mask also required
		pMskPBA = NULL;
	else
		if((pMskPBA = LocatePBAfor(m_MaskReadsetID, pChromMetadata->ChromID)) == NULL)
			continue;

	// what is the maximal sized unit that a thread would be expected to process in this chromosome
	if(pChromMetadata->ChromLen < MaxNumLoci)	// too much overhead if too many threads processing shorter chromosomes
		UnitSize = pChromMetadata->ChromLen; // a single thread to process this small chromosome
	else
		UnitSize = MaxNumLoci;

	int32_t StartLoci = 0;
	while(StartLoci < pChromMetadata->ChromLen)
		{
		if(m_TotWorkQueueEls == m_AllocWorkQueueEls)
			break;
		pWorkQueueEl->NumFndrs = NumFndrs;
		pWorkQueueEl->ChromID = pChromMetadata->ChromID;
		pWorkQueueEl->StartLoci = StartLoci;
		pWorkQueueEl->ChromLen = pChromMetadata->ChromLen;
		pWorkQueueEl->MaxNumLoci = min(MaxNumLoci, pChromMetadata->ChromLen - StartLoci);
		pWorkQueueEl->pMskPBA = pMskPBA;
		memcpy(pWorkQueueEl->pFounderPBAs,pPBAs,((size_t)NumFndrs * sizeof(uint8_t *)));
		m_TotWorkQueueEls++;
		pWorkQueueEl++;
		StartLoci += UnitSize;
		}
	}

// startup threads
StartWorkerThreads(m_NumThreads,		// there are this many threads in pool
					pReadsetMetadata->NumChroms);		// processing this number of chromosomes
while(!WaitAlignments(60))
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"AlignAlleleStacks: Generating allele stacks over %d founders",NumFndrs);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"AlignAlleleStacks: Completed generating %u allele stacks over %d founders", m_UsedAlleleStacks,NumFndrs);
if(m_pWorkQueueEls != NULL)
	{
	delete []m_pWorkQueueEls;
	m_pWorkQueueEls = NULL;	
	}
m_NumQueueElsProcessed = 0;
m_FastSerialise = 0;
return(m_UsedAlleleStacks);
}



int				
CCallHaplotypes::AlignFounderHaps(int32_t NumFndrs)			// number of founders to be processed against
{
int Rslt;
int32_t FounderID;
tsReadsetMetadata *pReadsetMetadata;
tsChromMetadata *pChromMetadata;
tsWorkQueueEl *pWorkQueueEl;
uint32_t CurChromMetadataIdx;
uint8_t *pPBAs[cMaxFounderReadsets+1];							// additional is to allow for the progeny readset PBAs
uint8_t *pMskPBA;
int32_t ChromIdx;
int32_t NumUnits;			// NumUnits is the total number of work units to be distributed over available threads for processing chromosomes
bool bInitHGBinSpecs;


pReadsetMetadata = &m_Readsets[0];		// founders were loaded first so 1st founder will be here
CurChromMetadataIdx = pReadsetMetadata->StartChromMetadataIdx;
if(m_pWorkQueueEls != NULL)				// ensuring only one work queue exists
	delete []m_pWorkQueueEls;
m_pWorkQueueEls = NULL;
m_AllocWorkQueueEls = 0;
pMskPBA = NULL;

Rslt = eBSFSuccess; // assume no errors!!!

// if cluster group bins not already initialised then use 
bInitHGBinSpecs = m_UsedHGBinSpecs == 0 ? true : false;

// determine maximal number of work queue elements required
CurChromMetadataIdx = pReadsetMetadata->StartChromMetadataIdx;
for(NumUnits = ChromIdx = 0; ChromIdx < pReadsetMetadata->NumChroms && CurChromMetadataIdx != 0; ChromIdx++)
	{
	pChromMetadata = &m_pChromMetadata[CurChromMetadataIdx - 1];
	NumUnits += 1; // safety - accounting for potential partial AccumBinSizes
	CurChromMetadataIdx = pChromMetadata->NxtChromMetadataIdx;
	}
if(NumUnits == 0)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "AlignFounderHaps: Number of thread work units is 0");
	return(-1);
	}
// sizing work queue to contain a maximum of NumUnits elements
m_AllocWorkQueueEls = NumUnits;
m_pWorkQueueEls = new tsWorkQueueEl[m_AllocWorkQueueEls];
memset(m_pWorkQueueEls, 0, sizeof(tsWorkQueueEl) * m_AllocWorkQueueEls);
m_TotWorkQueueEls = 0;
m_NumQueueElsProcessed = 0;
m_FastSerialise = 0;
pWorkQueueEl = m_pWorkQueueEls;

// iterate along PBAs which are common to all founders, plus mask readset if present, and then initialise elements of work to be undertaken by each thread
CurChromMetadataIdx = pReadsetMetadata->StartChromMetadataIdx;
for(ChromIdx = 0; ChromIdx < pReadsetMetadata->NumChroms && CurChromMetadataIdx != 0; ChromIdx++)
	{
	pChromMetadata = &m_pChromMetadata[CurChromMetadataIdx - 1];

	// all founders, and mask readset if present, must have PBAs for same chromosome
	for(FounderID = 1; FounderID <= NumFndrs; FounderID++)
		{
		if((pPBAs[FounderID-1] = LocatePBAfor(FounderID, pChromMetadata->ChromID)) == NULL)
			{
			char *pszChrom = LocateChrom(pChromMetadata->ChromID);
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "AlignFounderHaps: No PBA for chromosome '%s' in at least one founder '%s', skipping this chromosome", pszChrom, LocateReadset(FounderID));
			break;
			}
		}

	CurChromMetadataIdx = pChromMetadata->NxtChromMetadataIdx;		// ready for next iteration over the chromosomes
	if(FounderID <= NumFndrs)	// true only if PBAs not common to progeny and all founders, including mask readset
		continue;

	int32_t StartLoci;
	if(bInitHGBinSpecs) // use user specified bins or create default fixed sizes
		{
		int32_t BinSize;

		// create default fixed sizes
		for(StartLoci = 0; StartLoci < pChromMetadata->ChromLen; StartLoci += m_GrpHapBinSize)
			{
			BinSize = min(pChromMetadata->ChromLen - StartLoci, m_GrpHapBinSize);
			if(BinSize < 10) // any bin smaller is of no value!
				continue;
			if((Rslt = AddHGBinSpec(pChromMetadata->ChromID, StartLoci, BinSize, min(m_MinCentClustDist,BinSize-1),min(m_MaxCentClustDist,BinSize-1),m_MaxClustGrps)) == 0)
				{
				char* pszChrom = LocateChrom(pChromMetadata->ChromID);
				gDiagnostics.DiagOut(eDLInfo, gszProcName, "AlignFounderHaps: Error adding haplotyping group bin specifications on chromosome '%s'", pszChrom);
				if(m_pWorkQueueEls != NULL)
					{
					delete[]m_pWorkQueueEls;
					m_pWorkQueueEls = NULL;
					}
				m_NumQueueElsProcessed = 0;
				return(-1);
				}

			if(StartLoci == 0) // capture bin identifier for 1st bin on current chromosome or sequence
				pChromMetadata->HGBinID = Rslt;
			}
        }

	if(pChromMetadata->HGBinID == 0)   // skipping this chrom if no haplotype group specifications
		continue;

	StartLoci = 0;
	pWorkQueueEl->NumFndrs = NumFndrs;
	pWorkQueueEl->ChromID = pChromMetadata->ChromID;
	pWorkQueueEl->StartLoci = StartLoci;
	pWorkQueueEl->ChromLen = pChromMetadata->ChromLen;
	pWorkQueueEl->MaxNumLoci = pChromMetadata->ChromLen;
	pWorkQueueEl->pMskPBA = pMskPBA;
	memcpy(pWorkQueueEl->pFounderPBAs,pPBAs,((size_t)NumFndrs * sizeof(uint8_t *)));
	m_TotWorkQueueEls++;
	pWorkQueueEl++;
	}

if(m_TotWorkQueueEls > 0)
	{
	m_ProccessingBinID = 0;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"AlignFounderHaps: Starting to generate haplotype groupings for %d founders ",NumFndrs);
// startup threads
	StartWorkerThreads(m_NumThreads,		// there are this many threads in pool
		pReadsetMetadata->NumChroms);		// processing this number of chromosomes
	int WaitSecs = 60;
	while(!WaitAlignments(WaitSecs))
		{
		WaitSecs = min(WaitSecs += 60, 300);
		AcquireFastSerialise();
		int32_t LatestBinID = m_ProccessingBinID;
		ReleaseFastSerialise();
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "AlignFounderHaps: Generating haplotype groupings ... %0.3f%%",(LatestBinID * 100.0)/m_UsedHGBinSpecs);
		}

	gDiagnostics.DiagOut(eDLInfo, gszProcName, "AlignFounderHaps: Completed generating haplotype groupings");
	}
else
    gDiagnostics.DiagOut(eDLInfo, gszProcName, "AlignFounderHaps: Nothing to do, no haplotype groupings on %d founders required", NumFndrs);
if(m_pWorkQueueEls != NULL)
	{
	delete []m_pWorkQueueEls;
	m_pWorkQueueEls = NULL;	
	}
m_NumQueueElsProcessed = 0;
m_FastSerialise = 0;

// report haplotype groupings .....
Rslt = ReportHaplotypeGroups(NumFndrs);
return(Rslt);
}


int
CCallHaplotypes::GenFounderHaps(int32_t NumFndrs)			// number of founders to be processed against
{
int Rslt;
if((Rslt = AlignFounderHaps( NumFndrs)))		// number of founders to be processed against
	return(Rslt);
return(eBSFSuccess);
}

int
CCallHaplotypes::GenAlleleStacks(int32_t NumFndrs,			// number of founders to be processed against
								 int32_t MaxNumLoci)		// work queue items specify this number of loci for processing by threads
{
int Rslt;
uint32_t ASIdx;
tsAlleleStack *pAlleleStack;
if((Rslt = AlignAlleleStacks( NumFndrs,		// number of founders to be processed against
		  MaxNumLoci)) < 1)	// comparing ReadsetID's PBA against founder PBAs then work queue items specify this number of loci for processing by threads
	return(Rslt);
m_mtqsort.SetMaxThreads(m_NumThreads);
if(m_UsedAlleleStacks > 1)
    m_mtqsort.qsort(m_pAlleleStacks, (int64_t)m_UsedAlleleStacks, sizeof(tsAlleleStack), SortAlleleStacks);
pAlleleStack = m_pAlleleStacks;
for(ASIdx = 1; ASIdx <= m_UsedAlleleStacks; ASIdx++, pAlleleStack++)
	pAlleleStack->AlleleStackID = ASIdx;
return(m_UsedAlleleStacks);
}

inline void
CCallHaplotypes::BitsVectSet(uint16_t Bit,		// bit to set, range 0..4095
			ts4096Bits & BitsVect)
{
BitsVect.Bits[Bit / 64] |= ((uint64_t)0x01 << (Bit % 64));
}

inline void
CCallHaplotypes::BitsVectReset(uint16_t Bit,		// bit to reset, range 0..4095
		   ts4096Bits & BitsVect)
{
BitsVect.Bits[Bit / 64] &= ~((uint64_t)0x01 << (Bit % 64));
}


inline bool
CCallHaplotypes::BitsVectTest(uint16_t Bit,		// bit to test, range 0..4095
			ts4096Bits & BitsVect)
{
return(BitsVect.Bits[Bit / 64] & ((uint64_t)0x01 << (Bit % 64)) ? true : false);
}

inline bool
CCallHaplotypes::BitsVectEqual(ts4096Bits & BitsVectA,	// compare for equality
								ts4096Bits & BitsVectB)
{
if(!memcmp(&BitsVectA,&BitsVectB,sizeof(ts4096Bits)))
	return(true);
return(false);
}

inline void
CCallHaplotypes::BitsVectInitialise(bool Set,			// if true then initialse all bits as set, otherwise initialise all bits as reset
				ts4096Bits & BitsVect)
{
uint8_t InitWith;
if(Set)
	InitWith = 0xff;
else
	InitWith = 0;
memset(&BitsVect,InitWith,sizeof(ts4096Bits));
}

inline uint32_t
CCallHaplotypes::BitsVectCount(ts4096Bits& BitsVect)		// count number of set bits
{
uint32_t Count = 0;
uint64_t *pWord = BitsVect.Bits;
uint64_t Bit;
for(uint32_t WordIdx = 0; WordIdx < 64; WordIdx++, pWord++)
	{
	if(*pWord == 0)
		continue;
	for(Bit = 0x01; Bit != 0; Bit <<= 1)
		if(*pWord & Bit)
			Count++;
	}
return(Count);
}

inline uint32_t
CCallHaplotypes::BitsVectUnion(ts4096Bits &BitsVectA, ts4096Bits& BitsVectB)		// union (effective BitsVectA |= BitsVectB) bits in BitsVectA with BitsVectB with BitsVectA updated, returns number of bits set in BitsVectA 
{
uint32_t Count = 0;
uint64_t* pWordA = BitsVectA.Bits;
uint64_t* pWordB = BitsVectB.Bits;
uint64_t Bit;
for(uint32_t WordIdx = 0; WordIdx < 64; WordIdx++, pWordA++, pWordB++)
	{
	*pWordA |= *pWordB;
	if(*pWordA == 0)
		continue;
	for(Bit = 0x01; Bit != 0; Bit <<= 1)
		if(*pWordA & Bit)
			Count++;
	}
return(Count);
}


inline uint32_t 
CCallHaplotypes::BitsVectIntersect(ts4096Bits& BitsVectA, ts4096Bits& BitsVectB)	// intersect (effective BitsVectA &= BitsVectB) of bits in BitsVectA with BitsVectB with BitsVectA updated, returns number of set bits in BitsVectA
{
uint32_t Count = 0;
uint64_t* pWordA = BitsVectA.Bits;
uint64_t* pWordB = BitsVectB.Bits;
uint64_t Bit;
for(uint32_t WordIdx = 0; WordIdx < 64; WordIdx++, pWordA++, pWordB++)
	{
	*pWordA &= *pWordB;
	if(*pWordA == 0)
		continue;
	for(Bit = 0x01; Bit != 0; Bit <<= 1)
		if(*pWordA & Bit)
			Count++;
	}
return(Count);
}

uint32_t 
CCallHaplotypes::BitsVectClear(ts4096Bits& BitsVectA, ts4096Bits& BitsVectB)	    // clear bits in BitsVectA which are set in BitsVectB with BitsVectA updated, returns number of set bits in BitsVectA  
{
uint32_t Count = 0;
uint64_t* pWordA = BitsVectA.Bits;
uint64_t* pWordB = BitsVectB.Bits;
uint64_t Bit;
for(uint32_t WordIdx = 0; WordIdx < 64; WordIdx++, pWordA++, pWordB++)
	{
	*pWordA &= ~(* pWordB);
	if(*pWordA == 0)
		continue;
	for(Bit = 0x01; Bit != 0; Bit <<= 1)
		if(*pWordA & Bit)
			Count++;
	}
return(Count);
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

uint8_t 
CCallHaplotypes::LocateReadsetChromLociAlleles(int32_t ReadsetID,	// return alleles for this readset 
									  int32_t ChromID,		// on this chromosome
									  int32_t Loci)		// at this loci
{
uint8_t *pPBA;
if((pPBA=LocatePBAfor(ReadsetID,ChromID))==NULL)
	return(0);
return(pPBA[Loci]);
}

// NOTE: Readsets are checked for uniqueness as readsets must be unique within a given readetset type
int32_t		// returned readset identifier, < 1 if unable to accept this readset name
CCallHaplotypes::AddReadset(char* pszReadset, // associate unique identifier with this readset name
							uint8_t ReadsetType)	// 0: founder, 1: progeny, 2: control
{
int32_t ReadsetNameIdx;
int ReadsetNameLen;
char Type;
char *pszLAname;
Type = '0' + (char)ReadsetType;

// with any luck the sequence name will be same as the last accessed
if((pszLAname = LocateReadset(m_LAReadsetNameID)) != NULL)
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


int32_t		// returned Readset identifier, < 1 if unable to locate this Readset name
CCallHaplotypes::LocateReadset(char* pszReadset, // return unique identifier associated with this Readset name
							   uint8_t ReadsetType)	// 0: founder, 1: progeny, 2: control
{
int32_t ReadsetNameIdx;
char Type;
char *pszLAReadset;

Type = '0'+ (char)ReadsetType;
// with any luck the Readset name will be same as the last accessed
if((pszLAReadset = LocateReadset(m_LAReadsetNameID)) != NULL)
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

char* 
CCallHaplotypes::LocateReadset(int32_t ReadsetID)
{
int Idx;
Idx = ReadsetID & 0x0fffffff;			// mask out any potential ReadsetType
if(Idx < 1 || Idx > m_NumReadsetNames)
	return(NULL);
return(&(m_szReadsetNames[m_szReadsetIdx[Idx -1]+1])); // skipping lead char which is the ReadsetType
}


int
CCallHaplotypes::LoadPBACoverage(char* pszInPBA)   // file containing PBA
{
int32_t Rslt;
int32_t ReadsetID;
int32_t ChromID;
tsReadsetMetadata* pReadsetMetadata;
tsChromMetadata* pChromMetadata;
uint32_t CurChromMetadataIdx;
FILE* pInStream;
char szInWIG[1000];
char* pszInWIG;
char szLineBuff[1000];
char szChrom[100];
char* pszInBuff;
int32_t CurChromID;
int32_t PrevChromID;
int32_t StartLoci;
int32_t Span;
uint8_t *pAllele;
int32_t StartLociOfs;
int32_t Coverage;
uint64_t TotAllelesLen;
int32_t NumElsParsed;
int32_t WIGtype;    // 0: unknown, 1: fixed steps, 2 variable steps

int32_t LineNumb;

// compose WIG file name from PBA file name
// 
CUtility::AppendFileNameSuffix(szInWIG, pszInPBA, (char*)".covsegs.wig",'.');
pszInWIG = szInWIG;

// ensure WIG file is can be opened for reading
if((pInStream = fopen(pszInWIG,"r"))==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to file %s for reading error: %s",pszInWIG,strerror(errno));
	Reset();
	return(eBSFerrOpnFile);
	}

// load the PBA as per usual load, but then set actual PBAs to be 0
if((ReadsetID = LoadPBAFile(pszInPBA,0)) <= 0)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Errors loading PBA file '%s'",pszInPBA);
	fclose(pInStream);
	Reset();
	return(ReadsetID);
	}
m_Fndrs2Proc[ReadsetID-1] = 0x01;	// if loaded then assumption is that this founder will be processed
pReadsetMetadata = &m_Readsets[ReadsetID-1];
// iterate all PBA chroms and set PBA alleles to be all 0
pReadsetMetadata = &m_Readsets[ReadsetID-1];
CurChromMetadataIdx = pReadsetMetadata->StartChromMetadataIdx;
TotAllelesLen = 0; // used to later on normalise coverage by total length with coverage - in this release not actually used!!
for(ChromID = 1; ChromID <= pReadsetMetadata->NumChroms && CurChromMetadataIdx != 0; ChromID++)
	{
	pChromMetadata = &m_pChromMetadata[CurChromMetadataIdx - 1];
	pAllele = pChromMetadata->pPBAs;
	for(StartLoci=0;StartLoci < pChromMetadata->ChromLen;StartLoci++,pAllele++)
		if(*pAllele != 0)
			{
			*pAllele = 0;
			TotAllelesLen++;
			}
	CurChromMetadataIdx = pChromMetadata->NxtChromMetadataIdx;
	}

// load the WIG and where there is coverage then set the corresponding loci in the PBA to be max(254,log2(coverage) * 10)
// Note that the max limit of 254 is so that when later processing grouping and caching highest frequency coverage (alleles) then 254 can be distinguished from a marker of 255 used if not already cached 
LineNumb = 0;
WIGtype = 0;
pChromMetadata = NULL;
CurChromID = PrevChromID = 0;
Rslt = eBSFSuccess;
while(fgets(szLineBuff, sizeof(szLineBuff) - 1, pInStream) != NULL)
	{
	LineNumb++;
	pszInBuff = CUtility::ReduceWhitespace(szLineBuff);
	if(pszInBuff == NULL || *pszInBuff == '\0' || *pszInBuff == '\n')
		continue;

	if(!strnicmp(pszInBuff, "variableStep",12))
		{
		NumElsParsed = sscanf(pszInBuff, "variableStep chrom=%s span=%d", szChrom, &Span);
		if(NumElsParsed != 2)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "SimWig2PBA: Errors parsing line %d - \"%s\" - in file '%s'",LineNumb,pszInBuff,pszInPBA);
			Rslt = -1;
			break;
			}

		// check for valid chrom, a missing chrom relative to PBA could be that the user has specified a limit on number of chromosomes to process so can't treat as a fatal error
		if((CurChromID = LocateChrom(szChrom)) < 1)
			{
			CurChromID = 0;
			continue;
			}

		if(CurChromID != PrevChromID)
			{
			pChromMetadata = LocateChromMetadataFor(ReadsetID, CurChromID);
			PrevChromID = CurChromID;
			}

		if(Span < 1 || Span > pChromMetadata->ChromLen)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "SimWig2PBA: Span must be between 1..%d \"%s\" at line %d - \"%s\" - in file '%s'",LineNumb,pChromMetadata->ChromLen,pszInBuff,pszInPBA);
			Rslt = -1;
			break;
			}

		WIGtype = 2;
		continue;
		}
	else
		if(!strnicmp(pszInBuff, "fixedStep",10))
			{
			NumElsParsed = sscanf(pszInBuff, "fixedStep chrom=%s start=%d step=%d", szChrom, &StartLoci,&Span);
			if(NumElsParsed != 3)
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "SimWig2PBA: Errors parsing line %d - \"%s\" - in file '%s'",LineNumb,pszInBuff,pszInPBA);
				Rslt = -1;
				break;
				}

			// check for valid chrom, a missing chrom relative to PBA could be that the user has specified a limit on number of chromosomes to process so can't treat as a fatal error
			CurChromID = LocateChrom(szChrom);
			if((CurChromID = LocateChrom(szChrom)) < 1)
				{
				CurChromID = 0;
				continue;
				}

			if(CurChromID != PrevChromID)
				{
				pChromMetadata = LocateChromMetadataFor(ReadsetID, CurChromID);
				PrevChromID = CurChromID;
				}

			if(StartLoci < 1 || (StartLoci + Span - 1) >  pChromMetadata->ChromLen)
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "SimWig2PBA: StartLoci must be between 1..%d, Span %d not extending past end of chromosome \"%s\" at line %d - \"%s\" - in file '%s'",LineNumb,pChromMetadata->ChromLen,Span,pszInBuff,pszInPBA);
				Rslt = -1;
				break;
				}
			WIGtype = 1;
			continue;
			}
	if(CurChromID == 0)
		continue;

	switch(WIGtype) {
			case 0: // unrecognised format, assume some kind of comment or directive, slough
				continue;

			case 1: // fixed step
				NumElsParsed = sscanf(pszInBuff, "%d", &Coverage);
				break;

			case 2: // variable step
				NumElsParsed = sscanf(pszInBuff, "%d %d", &StartLoci, &Coverage);
				break;
		   }
	if(StartLoci < 1 || StartLoci > pChromMetadata->ChromLen) // ensure not writing past end of chromosome
		continue;
	pAllele = &pChromMetadata->pPBAs[StartLoci - 1];
	for(StartLociOfs = 0; StartLociOfs < Span; StartLociOfs++, StartLoci++, pAllele++)
		{
		if(StartLoci > pChromMetadata->ChromLen) // ensure not writing past end of chromosome
			break;
		*pAllele = min((int)log2(Coverage) * 10, 0x0fe);
		}
	}
fclose(pInStream);
return(Rslt >= 0 ? ReadsetID : Rslt);
}

int
CCallHaplotypes::CSV2WIG(char* pszInFile, // file containing haplotype groupings
	char* pszOutFile) // where to write as WIG formated file
{
int Rslt;
int32_t CurLineNumber;
int32_t NumFields;
char* pszChrom;
int32_t CurChromID;
int32_t PrevChromID;
int32_t BinStartLoci;
int32_t BinSize;
int32_t CentroidDistance;
int32_t NumChromBins;
int32_t TotNumChromBins;

if(m_hWIGOutFile != -1)
	{
	close(m_hWIGOutFile);
	m_hWIGOutFile = -1;
	}

if(m_pHGCSVFile != NULL) // should be, but better to be sure!
	{
	delete m_pHGCSVFile;
	m_pHGCSVFile = NULL;
	}

if((m_pHGCSVFile = new CCSVFile) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to instantiate CCSVfile");
	Reset();
	return(eBSFerrObj);
	}
m_pHGCSVFile->SetMaxFields(10); // only actually processing 1st 8 fields: chrom,loci, bin size, min distance, max distance, max haplotypes, actual distance, actual haplotypes

if((Rslt = m_pHGCSVFile->Open(pszInFile)) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to open file: '%s'", pszInFile);
	Reset();
	return(Rslt);
	}

if(m_pszWIGBuff == NULL)
	{
	if((m_pszWIGBuff = new uint8_t[cOutBuffSize]) == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "CSV2WIG: Memory allocation of %u bytes for WIG output buffering failed - %s", cOutBuffSize, strerror(errno));
		Reset();
		return(eBSFerrMem);
		};
	m_AllocWIGBuff = cOutBuffSize;
	}

#ifdef _WIN32
m_hWIGOutFile = open(pszOutFile, (O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC), (_S_IREAD | _S_IWRITE));
#else
if((m_hWIGOutFile = open64(pszOutFile, O_WRONLY | O_CREAT, S_IREAD | S_IWRITE)) != -1)
if(ftruncate(m_hWIGOutFile, 0) != 0)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "CSV2WIG: Unable to create/truncate %s - %s", pszOutFile, strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}
#endif
if(m_hWIGOutFile < 0)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "CSV2WIG: Unable to create/truncate %s - %s", pszOutFile, strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}

CurLineNumber = 0;
BinStartLoci = 0;
BinSize = 0;
CentroidDistance = 0;
PrevChromID = -1;
NumChromBins = 0;
TotNumChromBins = 0;
InitialiseWIGSpan();

while((Rslt = m_pHGCSVFile->NextLine()) > 0)		// onto next line containing fields
	{
	CurLineNumber++;
	if((NumFields = m_pHGCSVFile->GetCurFields()) < 8)	// must contain at least 8 fields
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Haplotype group bin specification file '%s' expected to contain a minimum of 8 fields, it contains %d at line %d", pszInFile, NumFields, CurLineNumber);
		Reset();
		return(eBSFerrParse);
		}

	// 1st row may be a header row, if so then skip this row
	if(CurLineNumber == 1 && m_pHGCSVFile->IsLikelyHeaderLine())
		continue;
	m_pHGCSVFile->GetText(1, &pszChrom);
	if(pszChrom == NULL || pszChrom[0] == '\0')
		{
		gDiagnostics.DiagOut(eDLWarn, gszProcName, "Unrecognised chromosome/sequence name at line %d in file '%s'", CurLineNumber, pszInFile);
		continue;
		}
	CurChromID = AddChrom(pszChrom);
	if(CurChromID != PrevChromID)
		{
		if(NumChromBins > 1)
			CompleteWIGSpan(true);
		NumChromBins = 0;
		PrevChromID = CurChromID;
		}
	m_pHGCSVFile->GetInt(2, &BinStartLoci);
	m_pHGCSVFile->GetInt(3, &BinSize);
	m_pHGCSVFile->GetInt(7, &CentroidDistance);
	NumChromBins++;
	AccumWIGCnts(CurChromID, BinStartLoci+1, BinSize, CentroidDistance); // WIG loci start from 1 not 0
	}
CompleteWIGSpan(true);
delete m_pHGCSVFile;
m_pHGCSVFile = NULL;

if(m_hWIGOutFile != -1)
	{
// commit output file
#ifdef _WIN32
	_commit(m_hWIGOutFile);
#else
	fsync(m_hWIGOutFile);
#endif
	close(m_hWIGOutFile);
	m_hWIGOutFile = -1;
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
char* pszChrom;

	// if existing span then write that span out
if(m_WIGChromID != 0 && m_WIGSpanLen > 0 && m_WIGSpanLoci > 0 && m_WIGSpanCnts > 0)
	{
	// has chrom and/or span changed since previously writing out a span?
	if(m_WIGChromID != m_WIGRptdChromID || m_WIGSpanLen != m_WIGRptdSpanLen)
		{
		pszChrom=LocateChrom(m_WIGChromID);
		m_WIGBuffIdx += sprintf((char *) & m_pszWIGBuff[m_WIGBuffIdx], "variableStep chrom=%s span=%d\n", pszChrom, m_WIGSpanLen);
		m_WIGRptdChromID = m_WIGChromID;
		m_WIGRptdSpanLen = m_WIGSpanLen;
		}
	m_WIGBuffIdx += sprintf((char *) & m_pszWIGBuff[m_WIGBuffIdx], "%d %d\n", m_WIGSpanLoci, (uint32_t)m_WIGSpanCnts);
	}
if((bWrite && m_WIGBuffIdx) || (m_WIGBuffIdx + 500) >  m_AllocWIGBuff)
	{
	if(!CUtility::RetryWrites(m_hWIGOutFile, m_pszWIGBuff, m_WIGBuffIdx))
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"RetryWrites() error");
		return(eBSFerrWrite);
		}
	m_WIGBuffIdx=0;
	}
m_WIGSpanLoci = 0;
m_WIGSpanLen = 0;
m_WIGSpanCnts = 0;
return(eBSFSuccess);
}



int
CCallHaplotypes::AccumWIGCnts(uint32_t ChromID,	// accumulate wiggle counts into variableStep spans over this chromosome
			uint32_t Loci,		// this loci - starts from 1 not 0!
	        uint32_t BinLen,     // span bin is this length
			uint32_t Cnts,		 // has this many counts attributed to the span
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
		m_WIGSpanLen = BinLen;
		m_WIGSpanCnts = Cnts;
		}
	return(eBSFSuccess);
	}

if(m_WIGSpanLen == 0 || m_WIGSpanCnts == 0) // starting a new span?
	{
	m_WIGSpanLoci = Loci;
	m_WIGSpanLen = BinLen;
	m_WIGSpanCnts = (uint64_t)Cnts;
	return(eBSFSuccess);
	}

if(Cnts != m_WIGSpanCnts)
	{
// if less than 5% difference between existing and new span counts then combine adjacent bins into a single span  
	Meanx100 = 100 * (uint32_t)(m_WIGSpanCnts / (uint64_t)m_WIGSpanLen);
	if((Cnts <= 100 && (Cnts * 100) != Meanx100) || (Meanx100 < (Cnts * 95) || Meanx100 >= (Cnts * 105)))
		{
		// write current span out
		if((Rslt = CompleteWIGSpan()) < 0)
			return(Rslt);
		m_WIGSpanLoci = Loci;
		m_WIGSpanLen = BinLen;
		m_WIGSpanCnts = (uint64_t)Cnts;
		return(eBSFSuccess);
		}
	}
m_WIGSpanLen = Loci - m_WIGSpanLoci + BinLen;
return(eBSFSuccess);
}



// sorting by ChromID.StartLoci ascending
int
CCallHaplotypes::SortAlleleStacks(const void* arg1, const void* arg2)
{
tsAlleleStack* pEl1 = (tsAlleleStack*)arg1;
tsAlleleStack* pEl2 = (tsAlleleStack*)arg2;
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

// sorting by ReadsetID.ChromID.StartLoci ascending
int
CCallHaplotypes::SortProgenyFndrAligns(const void* arg1, const void* arg2)
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
return(0);
}

// sorting by ChromID.StartLoci.Readset ascending
int
CCallHaplotypes::SortProgenyChromLociReadset(const void* arg1, const void* arg2)
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
if(pEl1->ReadsetID < pEl2->ReadsetID)
	return(-1);
if(pEl1->ReadsetID > pEl2->ReadsetID)
	return(1);
return(0);
}

// sorting by ChromID.StartLoci.NumLoci.Distance ascending
int
CCallHaplotypes::SortHGBinSpecs(const void* arg1, const void* arg2)
{
tsHGBinSpec* pEl1 = (tsHGBinSpec*)arg1;
tsHGBinSpec* pEl2 = (tsHGBinSpec*)arg2;
if(pEl1->ChromID < pEl2->ChromID)
	return(-1);
if(pEl1->ChromID > pEl2->ChromID)
	return(1);
if(pEl1->StartLoci < pEl2->StartLoci)
	return(-1);
if(pEl1->StartLoci > pEl2->StartLoci)
	return(1);
if(pEl1->NumLoci < pEl2->NumLoci)
	return(-1);
if(pEl1->NumLoci > pEl2->NumLoci)
	return(1);
if(pEl1->MaxCentroidDistance < pEl2->MaxCentroidDistance)
	return(-1);
if(pEl1->MaxCentroidDistance > pEl2->MaxCentroidDistance)
	return(1);
return(0);
}


