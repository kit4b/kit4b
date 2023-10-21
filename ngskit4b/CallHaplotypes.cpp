/*
'kit4b' is being released under the Opensource Software License Agreement (GPLv3)
'kit4b' is Copyright (c) 2019
Please contact Dr Stuart Stephen < stuartjs@g3web.com > if you have any questions regarding 'kit4b'.
*/


#include "stdafx.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#if _WIN32
//#define _CRTDBG_MAP_ALLOC
//#include <crtdbg.h>
#include <process.h>
#include "../libkit4b/commhdrs.h"
#else
#include <sys/mman.h>
#include <pthread.h>
#include "../libkit4b/commhdrs.h"
#endif

#include "./ngskit4b.h"
#include "./CallHaplotypes.h"

#include "../libkit4b/bgzf.h"

int CallHaplotypes(eModeCSH PMode,			// processing mode 0: report imputation haplotype matrix, 1: report both raw and imputation haplotype matrices, 2: additionally generate GWAS allowing visual comparisons, 3: allelic haplotype grouping,4: coverage haplotype grouping, 5: post-process haplotype groupings for DGTs,  6: post-process to WIG,  7 progeny PBAs vs. founder PBAs allelic association scores, 8 founder PBAs vs. founder PBAs allelic association scores, 9 grouping by scores, 10 differential group KMers
			int32_t LimitPrimaryPBAs,		// limit number of loaded primary or founder PBA files to this many. 0: no limits, > 0 sets upper limit
			int32_t ExprID,					// assign this experiment identifier for this PBA analysis
			int32_t SeedRowID,				// generated CSVs will contain monotonically incremented row identifiers seeded with this identifier  
			int32_t AffineGapLen,			// haplotype grouping maximal gap length at which gaps treated same as mismatches (affine scoring); if < 0 then no affine scoring with complete gap scored, if 0 then gaps are not scored, if > 0 then gaps upto and including AffineGapLen are scored
			int32_t GrpHapBinSize,			// when grouping into haplotypes (processing mode 3) then use this non-overlapping bin size for accumulation of differentials between all founder PBAs
			int32_t NumHapGrpPhases,		// number of phases over which to converge group consensus when haplotype clustering into groups - moves the balance between optimal group consensus and processing resource
			int32_t MinCentClustDist,		// haplotype groupings - in processing mode 3/4 only - minimum group centroid clustering distance
			int32_t MaxCentClustDist,		// haplotype groupings - in processing mode 3/4 only - maximum group centroid clustering distance
			int32_t MaxClustGrps,			// haplotype groupings - in processing mode 3/4 only - targeted maximum number of groups
			uint32_t SparseRepPropGrpMbrs,	// only apply sparse representative imputation if number of haplotype group members at least this 
			double SparseRepProp,			// if highest frequency (consensus) allele is 0x00 (no coverage) then if next highest frequency allele is at least this proportion of all members in haplotype group then treat next highest allele as
											// being as being the consensus. Objective to obtain an imputed allele in regions of sparse haplotype coverage
			int32_t MaxReportGrpDGTs,		// when calling group DGTs then, if non-zero - report this many highest scoring DGTs
			int32_t MinDGTGrpMembers,		// groups with fewer than this number of members - in processing mode 5 only -are treated as noise alleles these groups are not used when determining DGT group specific major alleles
			double MinDGTGrpPropTotSamples,	// haplotype groups - in processing mode 5 only -containing less than this proportion of all samples are treated as if containing noise and alleles in these groups are not used when determining DGT group specific major alleles 
			double MinDGTFmeasure,			// only accepting DGT loci with at least this F-measure score
			int32_t KMerSize,				// use this KMer sized sequences when identifying group segregation hammings
			int32_t MinKMerHammings,		// must be at least this hamming separating any two groups
			int32_t KMerNoneCoverage,				// if > 0 then allow upto this number of bp within a KMer to have no coverage, maximum is 50% of KMerSize
			int32_t FndrTrim5,				// trim this many aligned PBAs from 5' end of founder aligned segments - reduces false alleles due to sequencing errors
			int32_t FndrTrim3,				// trim this many aligned PBAs from 3' end of founder aligned segments - reduces false alleles due to sequencing errors
			int32_t ProgTrim5,				// trim this many aligned PBAs from 5' end of progeny aligned segments - reduces false alleles due to sequencing errors
			int32_t ProgTrim3,				// trim this many aligned PBAs from 3' end of progeny aligned segments - reduces false alleles due to sequencing errors
			int32_t WWRLProxWindow,			// proximal window size for Wald-Wolfowitz runs test
			int32_t OutliersProxWindow,		// proximal window size for outliers reduction
			char *pszMaskBPAFile,			// optional masking input BPA file, only process BPAs which are intersect of these BPAs and progeny plus founder BPAs
			char *pszChromFile,				// BED file containing reference assembly chromosome names and sizes
			char *pszAlleleScoreFile,	// CSV file containing allelic association scores previously generated by processing in eMCSHSrcVsRefs or eMCSHRefsVsRefs mode
			int NumFounderInputFiles,		// number of input founder file specs
			char *pszFounderInputFiles[],	// names of input founder PBA files (wildcards allowed)
			int NumProgenyInputFiles,		// number of input progeny file specs
			char* pszProgenyInputFiles[],	// names of input progeny PBA files (wildcards allowed)
			char* pszOutFile,				// loci haplotype calls output file (CSV format)
			int	NumIncludeChroms,			// number of chromosome regular expressions to include
			char **ppszIncludeChroms,		// array of include chromosome regular expressions
			int	NumExcludeChroms,			// number of chromosome expressions to exclude
			char **ppszExcludeChroms,		// array of exclude chromosome regular expressions
			int NumThreads);				// number of worker threads to use

#ifdef _WIN32
int callhaplotypes(int argc, char *argv[])
{
// determine my process name
_splitpath (argv[0], nullptr, nullptr, gszProcName, nullptr);
#else
int
callhaplotypes(int argc, char **argv)
{
// determine my process name
CUtility::splitpath ((char *)argv[0], nullptr, gszProcName);
#endif
int iFileLogLevel;			// level of file diagnostics
int iScreenLogLevel;		// level of file diagnostics
char szLogFile[_MAX_PATH];	// write diagnostics to this file
int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure
int NumberOfProcessors;		// number of installed CPUs
int NumThreads;				// number of threads (0 defaults to number of CPUs or a maximum of cMaxPBAWorkerThreads)

int Idx;

eModeCSH PMode;				// processing mode 
							// 0: report imputation haplotype matrix, 1: report both raw and imputation haplotype matrices, 2: additionally generate GWAS allowing visual comparisons, 
							// 3: allelic haplotype grouping,4: coverage haplotype grouping, 5: post-process haplotype groupings for DGTs,  6: post-process to WIG,  
							// 7 progeny PBAs vs. founder PBAs allelic association scores, 8 founder PBAs vs. founder PBAs allelic association scores, 9 grouping by scores 
							// 10 post-processing haplotype groupings for differential group KMers
							
int32_t LimitPrimaryPBAs;   // limit number of loaded primary or founder PBA files to this many. 0: no limits, > 0 sets upper limit
int32_t ExprID;			    // assign this experiment identifier to this PBA analysis
int32_t SeedRowID;          // generated CSVs will contain monotonically incremented row identifiers seeded with this identifier  
int32_t FndrTrim5;		    // trim this many aligned PBAs from 5' end of founder aligned segments - reduces false alleles due to sequencing errors
int32_t FndrTrim3;		    // trim this many aligned PBAs from 3' end of founder aligned segments - reduces false alleles due to sequencing errors
int32_t ProgTrim5;		    // trim this many aligned PBAs from 5' end of progeny aligned segments - reduces false alleles due to sequencing errors
int32_t ProgTrim3;		    // trim this many aligned PBAs from 3' end of progeny aligned segments - reduces false alleles due to sequencing errors
int32_t WWRLProxWindow;		// proximal window size for Wald-Wolfowitz runs test
int32_t OutliersProxWindow;	// proximal window size for outliers reduction
int32_t AffineGapLen;			// haplotype grouping maximal gap length at which gaps treated same as mismatches (affine scoring); if < 0 then no affine scoring with complete gap scored, if 0 then gaps are not scored, if > 0 then gaps upto and including AffineGapLen are scored
int32_t GrpHapBinSize;           // when grouping into haplotypes (processing mode 3) then use this non-overlapping bin size for accumulation of differentials between all founder PBAs
int32_t NumHapGrpPhases;         // number of phases over which to converge group consensus when haplotype clustering into groups - is a balance between optimal group consensus and processing resource
int32_t MaxReportGrpDGTs;     // when calling group DGTs then, if non-zero - report this many highest scoring DGTs
int32_t MinDGTGrpMembers;        // groups with fewer than this number of members - in processing mode 5 and 10 only -are treated as noise alleles these groups are not used when determining DGT group specific major alleles
double MinDGTGrpPropTotSamples;  // haplotype groups - in processing mode 5 only -containing less than this proportion of all samples are treated as if containing noise and alleles in these groups are not used when determining DGT group specific major alleles 
double MinDGTFmeasure;           // only accepting DGT loci with at least this F-measure score
int32_t MinCentClustDist;        // haplotype groupings - in processing mode 3 only - minimum group centroid clustering distance (default 10, min 1, clamped to bin size)
int32_t MaxCentClustDist;        // haplotype groupings - in processing mode 3 only - maximum group centroid clustering distance (default binsize, clamped to bin size)
int32_t MaxClustGrps;            // haplotype groupings - in processing mode 3 only - targeted maximum number of groups (default 5)
int32_t SparseRepPropGrpMbrs;		// only apply sparse representative imputation if number of haplotype group members at least this 
double SparseRepProp;				// if highest frequency (consensus) allele is 0x00 (no coverage) then if next highest frequency allele is at least this proportion of all members in haplotype group then treat next highest allele as
									// being as being the consensus. Objective to obtain an imputed allele in regions of sparse haplotype coverage
int32_t KMerSize;					// use this KMer sized sequences when identifying group segregation hammings
int32_t MinKMerHammings;			// must be at least this hamming separating any two groups
int32_t KMerNoneCoverage;				// if > 0 then allow upto this number of bp within a KMer to have no coverage, maximum is 50% of KMerSize
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
char szChromFile[_MAX_PATH];	// BED file containing chromosome names and sizes
char szAlleleScoreFile[_MAX_PATH];	// CSV file containing allelic association scores previously generated by processing in eMCSHSrcVsRefs or eMCSHRefsVsRefs mode

struct arg_lit *help = arg_lit0 ("h", "help", "print this help and exit");
struct arg_lit *version = arg_lit0 ("v", "version,ver", "print version information and exit");
struct arg_int *FileLogLevel = arg_int0 ("f", "FileLogLevel", "<int>", "Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0 ("F", "log", "<file>", "diagnostics log file");

struct arg_int *pmode = arg_int0 ("m", "mode", "<int>", "processing mode 0: report imputation haplotype matrix, 1: report both raw and imputation haplotype matrices, 2: additionally generate GWAS allowing visual comparisons, 3: allelic haplotype grouping,4: coverage haplotype grouping, 5: post-process haplotype groupings for DGTs, 6: post-process to WIG, 7 rel PBAs vs. ref PBAs, 8 ref PBAs vs. ref PBAs  allelic association scores,  9 grouping by allelic association scores, 10 post-processing haplotype groupings for differential group KMers (default 0)");
struct arg_int *limitprimarypbas = arg_int0 ("l", "limit", "<int>", " limit number of loaded primary or founder PBA files to this many. 0: no limits, > 0 sets upper limit (default 0)");

struct arg_int* exprid = arg_int0("e","exprid","<int>","assign this experiment identifier for haplotypes called (default 1");
struct arg_int* rowid = arg_int0("E","rowid","<int>","use this monotonically incremented row identifier seed in output CSVs (default 1");

struct arg_int* fndrtrim5 = arg_int0("y", "fndrtrim5", "<int>", "trim this many aligned PBAs from 5' end of founder aligned segments (default 5)");
struct arg_int* fndrtrim3 = arg_int0("Y", "fndrtrim3", "<int>", "trim this many aligned PBAs from 3' end of founder aligned segments (default 5)");
struct arg_int* progtrim5 = arg_int0("x", "progtrim5", "<int>", "trim this many aligned PBAs from 5' end of progeny aligned segments (default 5)");
struct arg_int* progtrim3 = arg_int0("X", "progtrim3", "<int>", "trim this many aligned PBAs from 3' end of progeny aligned segments (default 5)");

struct arg_int* wwrlproxwindow = arg_int0("w", "wwrlproxwindow", "<int>", "proximal window size for Wald-Wolfowitz runs test (default 1000000)");
struct arg_int* outliersproxwindow = arg_int0("W", "outliersproxwindow", "<int>", "proximal window size for outliers reduction (default 1000000)");

struct arg_int* maxreportgrpDGTs = arg_int0("N", "maxreportgrpDGTs", "<int>", "Report at most, this many highest scoring haplotype grouping DGT loci, 0 for no limits, (default 10000000)");
struct arg_int* minDGTgrpmembers = arg_int0("n", "grpDGTmbrs", "<int>", "Minimum DGT haplotype group members, mode 5 and 10 only, (default 10)");
struct arg_dbl* minDGTgrpprptotsamples = arg_dbl0("q", "grpDGTsamples", "<dbl>", "Minimum DGT haplotype group members as proportion of total samples, mode 5 only, (default 0.10)");
struct arg_dbl* minDGTgrpfmeasure = arg_dbl0("Q", "grpDGTfmeasure", "<dbl>", "Only accepting DGT loci with at least this F-measure score (default 0.90)");

struct arg_int* minsparsegrpmbrs = arg_int0("s", "minsparsegrpmbrs", "<int>", "Minimum number of group members required for sparse imputation (default 0 for no sparse imputation)");
struct arg_dbl* minsparseprop = arg_dbl0("S", "minsparseprop", "<dbl>", "Minimum proportion of group members with called alleles before applying sparse imputation (default 0.25)");
struct arg_int* affinegaplen = arg_int0("a", "affinegaplen", "<int>", "haplotype groupings - in processing mode 3/4 only -affine gap length for mismatch scoring (default 3)");
struct arg_int* grphapbinsize = arg_int0("g", "grphapbinsize", "<int>", "haplotype groupings - processing modes 3/4 and 7/8 only - bin size");
struct arg_int* maxclustgrps = arg_int0("G", "maxclustgrps", "<int>", "haplotype groupings - in processing mode 3/4 only - targeted maximum number of groups (default 5)");
struct arg_int* gpphases = arg_int0("p", "gpphases", "<int>", "haplotype groupings - in processing mode 3/4 only - max number of cluster processing phases (default 10)");

struct arg_int* mincentclustdist = arg_int0("d", "mincentclustdist", "<int>", "haplotype groupings - in processing mode 3/4 only - minimum group centroid clustering distance (default 5, min 1, clamped to maxcentclustdist-1)");
struct arg_int* maxcentclustdist = arg_int0("D", "maxcentclustdist", "<int>", "haplotype groupings - in processing mode 3/4 only - maximum group centroid clustering distance (default binsize-1, clamped to bin size-1)");

struct arg_int* kmersize = arg_int0("k", "kmersize", "<int>", "KMers - bp size used to check for group hammings in mode 10 (default 25)");
struct arg_int* minkmerhamming = arg_int0("K", "minkmerhamming", "<int>", "KMers - hammings between any two groups must be at least this in mode 10 (default 2)");
struct arg_int* kmernonecoverage = arg_int0("U", "KMerNoneCoverage", "<int>", "KMers - allows KMers to include this number of sites without coverage in mode 10 (default is all KMer sites must have coverage)");

struct arg_str  *ExcludeChroms = arg_strn("Z","chromexclude",	"<string>",0,cMaxExcludeChroms,"high priority - regular expressions defining chromosomes to exclude");
struct arg_str  *IncludeChroms = arg_strn("z","chromeinclude",	"<string>",0,cMaxIncludeChroms,"low priority - regular expressions defining chromosomes to include");

struct arg_file* chromfile = arg_file1("c", "chromfile", "<file>", "input BED file containing chromosome names and sizes");

struct arg_file* allelescorefile = arg_filen("A", "allelescorefile", "<file>", 0, 1, "input allelic association scores file as generated by 7 rel PBAs vs. ref PBAs or 8 ref PBAs vs. ref PBAs");

struct arg_file *founderfiles = arg_filen("I", "founderfiles", "<file>", 0,cMaxFounderFileSpecs,"founder input BPA file(s), wildcards allowed, limit of 500 founder filespecs supported");
struct arg_file *progenyfiles = arg_filen("i", "inprogenyfile", "<file>",0, cMaxProgenyFileSpecs, "progeny input BPA file(s), wildcards allowed, limit of 500 progeny filespecs supported");
struct arg_file *maskbpafile = arg_file0("j", "inmaskbpa", "<file>", "optional masking input BPA file, or haplotype grouping specification file(s) or previously generated haplotype grouping file for post-processing");
struct arg_file *outfile = arg_file1("o", "out", "<file>", "loci haplotype calls output prefix (outputs CSV and WIG format)");
struct arg_int *threads = arg_int0("T","threads","<int>","number of processing threads 0..64 (defaults to 0 which limits threads to maximum of 64 CPU cores)");
struct arg_end *end = arg_end (200);

void *argtable[] = { help,version,FileLogLevel,LogFile,
					pmode,exprid,rowid,limitprimarypbas,kmersize,minkmerhamming,kmernonecoverage,maxreportgrpDGTs,minDGTgrpmembers,minDGTgrpfmeasure,minDGTgrpprptotsamples,minsparsegrpmbrs,minsparseprop,
					grphapbinsize,gpphases,mincentclustdist,maxcentclustdist,affinegaplen,maxclustgrps,fndrtrim5,fndrtrim3,progtrim5,progtrim3,wwrlproxwindow,outliersproxwindow,
					maskbpafile,chromfile,IncludeChroms,ExcludeChroms,progenyfiles,founderfiles, allelescorefile, outfile,threads,end };

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
		MaxReportGrpDGTs = 0;
		MinDGTGrpMembers = 0; 
		MinDGTGrpPropTotSamples = 0.0;
		MinDGTFmeasure = 0.0;
		MinKMerHammings = cDfltKMerHammings;
		KMerSize = cDfltKMerSize;
		KMerNoneCoverage = 0;

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


		if(PMode == eMCSHKMerGrps)
			{
			KMerSize = kmersize->count ? kmersize->ival[0] : cDfltKMerSize;
			if(KMerSize < 10)			// silent clamping to between 10 and 1000
				KMerSize = 10;
			else
				if (KMerSize > 1000)
					KMerSize = 1000;

			MinKMerHammings = minkmerhamming->count ? minkmerhamming->ival[0] : cDfltKMerHammings;
			if (MinKMerHammings < 1)			// silent clamping to between 1 and 5
				MinKMerHammings = 1;
			else
				if (MinKMerHammings > 5)
					MinKMerHammings = 5;

			
			KMerNoneCoverage = kmernonecoverage->count ? kmernonecoverage->ival[0] : 0;
			if(KMerNoneCoverage < 0)
				KMerNoneCoverage = 0;
			else
				if(KMerNoneCoverage > KMerSize/2)
					KMerNoneCoverage = KMerSize/2;

			MinDGTGrpMembers = minDGTgrpmembers->count ? minDGTgrpmembers->ival[0] : cDfltMinDGTGrpMembers;
			if (MinDGTGrpMembers < 1)
				MinDGTGrpMembers = 1;
			else
				if (MinDGTGrpMembers > 10000)
					MinDGTGrpMembers = 10000;
			}

		if(PMode <= eMCSHRaw)		// need to ensure that at most 4 founders will be processed for progeny parentage 
			{
			if(LimitPrimaryPBAs == 0 || LimitPrimaryPBAs > 4)
				LimitPrimaryPBAs = 4;
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

		SeedRowID = rowid->count ? rowid->ival[0] : 1;
		if(SeedRowID < 1)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Seed row identifier '-E%d' must be in range 1..1000000000\n",SeedRowID);
			exit(1);
			}
		else
			if(SeedRowID > 1000000000)
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Seed row identifier '-E%d' must be in range 1..1000000000\n",SeedRowID);
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

		SparseRepPropGrpMbrs = 0;
		SparseRepProp = 0.0;
		if (PMode == eMCSHAllelicHapsGrps || PMode == eMCSHDGTHapsGrps)
			{
			SparseRepPropGrpMbrs = minsparsegrpmbrs->count ? minsparsegrpmbrs->ival[0] : cDfltSparseRepPropGrpMbrs;
			if (SparseRepPropGrpMbrs < 0)	// silently clamping to reasonable values
				SparseRepPropGrpMbrs = 0;
			else
				if (SparseRepPropGrpMbrs > 100000)
					SparseRepPropGrpMbrs = 100000;
			if (SparseRepPropGrpMbrs > 0)
				SparseRepProp = minsparseprop->count ? minsparseprop->dval[0] : cDfltSparseRepProp;
			if (SparseRepProp < 0.05)    // silently clamping to reasonable values
				SparseRepProp = 0.05;
			else
				if(SparseRepProp > 90.0)
					SparseRepProp = 90.0;
			}


		if(PMode == eMCSHDGTHapsGrps)
			{
			MaxReportGrpDGTs = maxreportgrpDGTs->count ? maxreportgrpDGTs->ival[0] : cDfltMaxReportGrpDGTs;
			if(MaxReportGrpDGTs < 1)
				MaxReportGrpDGTs = 0;
			else
				if(MaxReportGrpDGTs > 2000000000)
					MaxReportGrpDGTs = 2000000000;

			MinDGTGrpMembers = minDGTgrpmembers->count ? minDGTgrpmembers->ival[0] : cDfltMinDGTGrpMembers;
			if(MinDGTGrpMembers < 1)
				MinDGTGrpMembers = 1;
			else
				if(MinDGTGrpMembers > 10000)
					MinDGTGrpMembers = 10000;
			MinDGTGrpPropTotSamples = minDGTgrpprptotsamples->count ? minDGTgrpprptotsamples->dval[0] : cDfltMinDGTGrpPropTotSamples;
			MinDGTGrpPropTotSamples = ((int)(MinDGTGrpPropTotSamples * 1000)) / 1000.0;
			if(MinDGTGrpPropTotSamples < 0.001)
				MinDGTGrpPropTotSamples = 0.001;
			else
				if(MinDGTGrpPropTotSamples > 0.250)
					MinDGTGrpPropTotSamples = 0.250;
			MinDGTFmeasure  = minDGTgrpfmeasure->count ? minDGTgrpfmeasure->dval[0] : cDfltMinDGTFmeasure;
			if(MinDGTFmeasure < 0.5)
				MinDGTFmeasure = 0.5;
			else
				if(MinDGTFmeasure > 0.99)
					MinDGTFmeasure = 0.99;
			}

		if(IncludeChroms->count)
			{
			for(NumIncludeChroms = Idx = 0; NumIncludeChroms < cMaxIncludeChroms && Idx < IncludeChroms->count; Idx++)
				{
				pszIncludeChroms[Idx] = nullptr;
				if(pszIncludeChroms[NumIncludeChroms] == nullptr)
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
				pszExcludeChroms[Idx] = nullptr;
				if(pszExcludeChroms[NumExcludeChroms] == nullptr)
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
			if(GrpHapBinSize < 1000)        // silently clamping to a more reasonable size!
				GrpHapBinSize = 1000;       // upper limit will be actual chromosome size

            MinCentClustDist = mincentclustdist->count ? mincentclustdist->ival[0] : cDfltMinCentClustDist;
			if(MinCentClustDist < 1) 
				MinCentClustDist = 1;
			else
				if(MinCentClustDist >= GrpHapBinSize) // silently clamp grouping size to window size-1!
					MinCentClustDist = GrpHapBinSize-1;
			MaxCentClustDist = maxcentclustdist->count ? maxcentclustdist->ival[0] : GrpHapBinSize-1;
			if(MaxCentClustDist < MinCentClustDist) 
				MaxCentClustDist = MinCentClustDist;
			else
				if(MaxCentClustDist >= GrpHapBinSize) // silently clamp grouping size to window size-1!
					MaxCentClustDist = GrpHapBinSize-1;

			MaxClustGrps = maxclustgrps->count ? maxclustgrps->ival[0] : cDfltMaxClustGrps;
			if(MaxClustGrps < 2) 
				MaxClustGrps = 2;
			else
				if(MaxClustGrps > cMaxClustGrps) // silently clamp
					MaxClustGrps = cMaxClustGrps;

			AffineGapLen = affinegaplen->count ? affinegaplen->ival[0] : cDfltAffineGapLen;
			if (AffineGapLen < 0)		// 0 allowed, means no gap scoring
				AffineGapLen = -1;		// -1 means score over complete gap length
			else
				if (AffineGapLen > cMaxAffineGapLen) // silently clamp
					AffineGapLen = cMaxAffineGapLen;
			}
		else
			{
			if(PMode == eMCSHSrcVsRefs || PMode == eMCSHRefsVsRefs)
				{
				GrpHapBinSize = grphapbinsize->count ? grphapbinsize->ival[0] : cDfltProgFndrBinScoring;
				if (GrpHapBinSize < 10000)        // min bin size supported!
					GrpHapBinSize = 0;			  // if 0 then will be actual chromosome size
				}
			else
				GrpHapBinSize = 0;
			MinCentClustDist = 0;
			MaxCentClustDist = 0;
			MaxClustGrps = 0;
			AffineGapLen = 0;
			}



		FndrTrim5 = 0;
		FndrTrim3 = 0;
		ProgTrim5 = 0;
		ProgTrim3 = 0;

		if(PMode < eMCSHAllelicHapsGrps || PMode == eMCSHSrcVsRefs || PMode == eMCSHKMerGrps || PMode == eMCSHRefsVsRefs)
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

			if(PMode != eMCSHKMerGrps)
				{
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
				}
			}

		if(PMode < eMCSHAllelicHapsGrps)
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


	if (chromfile->count)
		{
		strncpy(szChromFile, chromfile->filename[0], _MAX_PATH);
		szChromFile[_MAX_PATH-1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(szChromFile);
		}
	else
		szChromFile[0] = 0;
	if(szChromFile[0] == 0)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "No BED file containing chromosome names and sizes specified");
		exit(1);
		}

	szAlleleScoreFile[0] = 0;
	if(PMode == eMCSHGroupScores)
		{
		if (allelescorefile->count)
			{
			strncpy(szAlleleScoreFile, allelescorefile->filename[0], _MAX_PATH);
			szAlleleScoreFile[_MAX_PATH - 1] = '\0';
			CUtility::TrimQuotedWhitespcExtd(szAlleleScoreFile);
			}
		if(szAlleleScoreFile[0] == 0)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "No CSV file containing allelic association scores");
			exit(1);
			}
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

	if(PMode < eMCSHAllelicHapsGrps || PMode == eMCSHSrcVsRefs)
		{
		if(progenyfiles->count)
			{
			for(Idx = 0; NumProgenyInputFiles < cMaxProgenyReadsets && Idx < progenyfiles->count; Idx++)
				{
				pszProgenyInputFiles[Idx] = nullptr;
				if(pszProgenyInputFiles[NumProgenyInputFiles] == nullptr)
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
	if(PMode == eMCSHKMerGrps || !(PMode == eMCSHGrpDist2WIG || PMode == eMCSHGroupScores))
		{
		if(founderfiles->count)
			{
			for(NumFounderInputFiles = Idx = 0; NumFounderInputFiles < cMaxFounderFileSpecs && Idx < founderfiles->count; Idx++)
				{
				pszFounderInputFiles[Idx] = nullptr;
				if(pszFounderInputFiles[NumFounderInputFiles] == nullptr)
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
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: After removal of whitespace, no input %s file specified with '-j<filespec>' option)\n",(PMode == eMCSHAllelicHapsGrps || PMode == eMCSHCoverageHapsGrps || PMode == eMCSHDGTHapsGrps || PMode == eMCSHKMerGrps) ? "haplotype group clustering" : "masking BPA");
			exit(1);
			}
		}
	else
		szMaskBPAFile[0] = '\0';

	if(PMode == eMCSHDGTHapsGrps || PMode == eMCSHKMerGrps || PMode == eMCSHGrpDist2WIG)
		{
		if(szMaskBPAFile[0] == '\0')
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: No previously generated haplotype grouping file specified with '-z<filespec>' option)\n");
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
		case eMCSHDGTHapsGrps:
			pszDescr = "post-processing haplotype groupings for differential DGT loci";
			break;
		case eMCSHGrpDist2WIG:
			pszDescr = "Post-processing haplotype grouping centroid distances into WIG format";
			break;

		case eMCSHSrcVsRefs:
			pszDescr = "Source PBAs vs. Reference PBAs allelic association scores";
			break;

		case eMCSHRefsVsRefs:
			pszDescr = "Source PBAs vs. Reference PBAs allelic association scores";
			break;

		case eMCSHGroupScores:
			pszDescr = "Grouping previously generated allelic association scores";

		case eMCSHKMerGrps:
			pszDescr = "Post-processing haplotype groupings for differential group KMers";
			break;
		}

	gDiagnostics.DiagOutMsgOnly (eDLInfo, "Calling haplotypes : '%s'", pszDescr);

	if(PMode != eMCSHKMerGrps)
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "Experiment identifier: %d", ExprID);
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "Monotonically incremented row identifier seed: %d", SeedRowID);
		}

	for(Idx = 0; Idx < NumIncludeChroms; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"reg expressions defining chroms to include: '%s'",pszIncludeChroms[Idx]);
	for(Idx = 0; Idx < NumExcludeChroms; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"reg expressions defining chroms to exclude: '%s'",pszExcludeChroms[Idx]);

	if(LimitPrimaryPBAs > 0)
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "Limiting number of loaded primary or founder PBA files to first: %d",LimitPrimaryPBAs);

	if(PMode < eMCSHAllelicHapsGrps || PMode == eMCSHSrcVsRefs || PMode == eMCSHRefsVsRefs || PMode == eMCSHKMerGrps)
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "Trim this many aligned PBAs from 5' end of founder aligned segments: %d", FndrTrim5);
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "Trim this many aligned PBAs from 3' end of founder aligned segments: %d", FndrTrim3);
		if (PMode != eMCSHKMerGrps)
			{
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "Trim this many aligned PBAs from 5' end of progeny aligned segments: %d", ProgTrim5);
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "Trim this many aligned PBAs from 3' end of progeny aligned segments: %d", ProgTrim3);
			}
		}

	if(PMode == eMCSHGroupScores)
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "Loading association scores from: '%s'", szAlleleScoreFile);
		}
	else
		{
		if (PMode == eMCSHSrcVsRefs || PMode == eMCSHRefsVsRefs)
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "Scoring over bins of this window size: %d", GrpHapBinSize);
		else
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
						if(AffineGapLen < 0)
							gDiagnostics.DiagOutMsgOnly(eDLInfo, "Gaps are full length scored");
						else
							if(AffineGapLen == 0)
								gDiagnostics.DiagOutMsgOnly(eDLInfo, "Gaps are not scored");
							else
								gDiagnostics.DiagOutMsgOnly(eDLInfo, "Gaps are scored for 1st affine bp: %d", AffineGapLen);
						}
					}
				else
					{
					if(PMode == eMCSHDGTHapsGrps)
						{
						if(MaxReportGrpDGTs)
							gDiagnostics.DiagOutMsgOnly(eDLInfo, "Maximum number of highest scoring DGT group loci to report: %d", MaxReportGrpDGTs);
						gDiagnostics.DiagOutMsgOnly(eDLInfo, "Minimum DGT haplotype group members: %d", MinDGTGrpMembers);
						gDiagnostics.DiagOutMsgOnly(eDLInfo, "Minimum DGT haplotype group members as proportion of total samples: %.3f", MinDGTGrpPropTotSamples);
						gDiagnostics.DiagOutMsgOnly(eDLInfo, "Minimum DGT haplotype group allele Log2 F-measure score: %.3f", MinDGTFmeasure);
						}
					else 
						if(PMode == eMCSHKMerGrps)
							{
							gDiagnostics.DiagOutMsgOnly(eDLInfo, "Minimum haplotype group members: %d", MinDGTGrpMembers);
							gDiagnostics.DiagOutMsgOnly(eDLInfo, "KMer size: %dbp", KMerSize);
							gDiagnostics.DiagOutMsgOnly(eDLInfo, "Minimum hammings separating KMer groups: %d", MinKMerHammings);
							gDiagnostics.DiagOutMsgOnly(eDLInfo, "KMers can contain at most this many bp with no coverage: %d", KMerNoneCoverage);
							}
					gDiagnostics.DiagOutMsgOnly(eDLInfo, "Loading previously generated haplotype groupings from : '%s'", szMaskBPAFile);
					}
			}
		}


	if (SparseRepPropGrpMbrs > 0 && (PMode == eMCSHAllelicHapsGrps || PMode == eMCSHDGTHapsGrps || PMode == eMCSHKMerGrps))
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "Process for sparse imputation if number of group members more than: %d", SparseRepPropGrpMbrs);
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "Use alleles for sparse imputation if at least this proportion of all group members: %.3f", SparseRepProp);
		}

	if(PMode < eMCSHGrpDist2WIG || PMode == eMCSHKMerGrps)
		{
		for(Idx = 0; Idx < NumFounderInputFiles; Idx++)
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "Sample or Founder file : '%s'", pszFounderInputFiles[Idx]);

		if(NumProgenyInputFiles && PMode != eMCSHKMerGrps)
			for(Idx = 0; Idx < NumProgenyInputFiles; Idx++)
				gDiagnostics.DiagOutMsgOnly(eDLInfo, "Progeny file : '%s'", pszProgenyInputFiles[Idx]);
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo, "BED containing chromosome names and sizes : '%s'", szChromFile);
	gDiagnostics.DiagOutMsgOnly (eDLInfo, "Output file : '%s'", szOutFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"number of threads : %d",NumThreads);

#ifdef _WIN32
	SetPriorityClass (GetCurrentProcess (), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start ();
	Rslt = CallHaplotypes(PMode,			// processing mode 0: report imputation haplotype matrix, 1: report both raw and imputation haplotype matrices, 2: additionally generate GWAS allowing visual comparisons, 3: allelic haplotype grouping,4: coverage haplotype grouping, 5: post-process haplotype groupings for DGTs,  6: post-process to WIG,  7 progeny PBAs vs. founder PBAs allelic association scores, 8 founder PBAs vs. founder PBAs allelic association scores, 9 grouping by scores, 10 differential group KMers
					ExprID,			         // assign this experiment identifier for this PBA analysis
					SeedRowID,                  // generated CSVs will contain monotonically incremented row identifiers seeded with this identifier  
					LimitPrimaryPBAs,        // limit number of loaded primary or founder PBA files to this many. 0: no limits, > 0 sets upper limit			
					AffineGapLen,			// haplotype grouping maximal gap length at which gaps treated same as mismatches (affine scoring); if < 0 then no affine scoring with complete gap scored, if 0 then gaps are not scored, if > 0 then gaps upto and including m_AffineGapLen are scored
					GrpHapBinSize,           // when grouping into haplotypes (processing mode 3) then use this non-overlapping bin size for accumulation of differentials between all founder PBAs
					NumHapGrpPhases,         // number of phases over which to converge group consensus when haplotype clustering into groups - is a balance between optimal group consensus and processing resource
					MinCentClustDist,        // haplotype groupings - in processing mode 3 only - minimum group centroid clustering distance (default 500, min 1, clamped to bin size)
					MaxCentClustDist,        // haplotype groupings - in processing mode 3 only - maximum group centroid clustering distance (default mincentclustdist, clamped to bin size)
					MaxClustGrps,       // haplotype groupings - in processing mode 3 only - targeted maximum number of groups (default 5)
					SparseRepPropGrpMbrs,	// only apply sparse representative imputation if number of haplotype group members at least this 
					SparseRepProp,			// if highest frequency (consensus) allele is 0x00 (no coverage) then if next highest frequency allele is at least this proportion of all members in haplotype group then treat next highest allele as
					MaxReportGrpDGTs,     // when calling group DGTs then, if non-zero - report this many highest scoring DGTs
					MinDGTGrpMembers,        // groups with fewer than this number of members - in processing mode 5 only -are treated as noise alleles these groups are not used when determining DGT group specific major alleles
					MinDGTGrpPropTotSamples,  // haplotype groups - in processing mode 5 only -containing less than this proportion of all samples are treated as if containing noise and alleles in these groups are not used when determining DGT group specific major alleles 
					MinDGTFmeasure,           // only accepting DGT loci with at least this F-measure score
					KMerSize,					// use this KMer sized sequences when identifying group segregation hammings
					MinKMerHammings,			// must be at least this hamming separating any two groups
					KMerNoneCoverage,				// if > 0 then allow upto this number of bp within a KMer to have no coverage, maximum is 50% of KMerSize
					FndrTrim5,			// trim this many aligned PBAs from 5' end of founder aligned segments - reduces false alleles due to sequencing errors
					FndrTrim3,			// trim this many aligned PBAs from 3' end of founder aligned segments - reduces false alleles due to sequencing errors
					ProgTrim5,			// trim this many aligned PBAs from 5' end of progeny aligned segments - reduces false alleles due to sequencing errors
					ProgTrim3,			// trim this many aligned PBAs from 3' end of progeny aligned segments - reduces false alleles due to sequencing errors
					WWRLProxWindow,		// proximal window size for Wald-Wolfowitz runs test
					OutliersProxWindow,	// proximal window size for outliers reduction
					szMaskBPAFile,			// optional input masking BPA file or Haplotype group clustering specifications file(s)
					szChromFile,			// BED file containing reference assembly chromosome names and sizes
					szAlleleScoreFile,	// CSV file containing allelic association scores previously generated by processing in eMCSHSrcVsRefs or eMCSHRefsVsRefs mode
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

int CallHaplotypes(eModeCSH PMode,	// processing mode 0: report imputation haplotype matrix, 1: report both raw and imputation haplotype matrices, 2: additionally generate GWAS allowing visual comparisons, 3: allelic haplotype grouping,4: coverage haplotype grouping, 5: post-process haplotype groupings for DGTs,  6: post-process to WIG,  7 progeny PBAs vs. founder PBAs allelic association scores, 8 founder PBAs vs. founder PBAs allelic association scores, 9 grouping by scores, 10 differential group KMers
			int32_t ExprID,			         // assign this experiment identifier for this PBA analysis
			int32_t SeedRowID,                  // generated CSVs will contain monotonically incremented row identifiers seeded with this identifier  
			int32_t LimitPrimaryPBAs,                // limit number of loaded primary or founder PBA files to this many. 0: no limits, > 0 sets upper limit
			int32_t AffineGapLen,			// haplotype grouping maximal gap length at which gaps treated same as mismatches (affine scoring); if < 0 then no affine scoring with complete gap scored, if 0 then gaps are not scored, if > 0 then gaps upto and including m_AffineGapLen are scored
			int32_t GrpHapBinSize,           // when grouping into haplotypes (processing mode 3) then use this non-overlapping bin size for accumulation of differentials between all founder PBAs
			int32_t NumHapGrpPhases,         // number of phases over which to converge group consensus when haplotype clustering into groups - moves the balance between optimal group consensus and processing resource
			int32_t MinCentClustDist,        // haplotype groupings - in processing mode 3/4 only - minimum group centroid clustering distance
			int32_t MaxCentClustDist,        // haplotype groupings - in processing mode 3/4 only - maximum group centroid clustering distance
			int32_t MaxClustGrps,            // haplotype groupings - in processing mode 3 only - targeted maximum number of groups
			uint32_t SparseRepPropGrpMbrs,	// only apply sparse representative imputation if number of haplotype group members at least this 
			double SparseRepProp,			// if highest frequency (consensus) allele is 0x00 (no coverage) then if next highest frequency allele is at least this proportion of all members in haplotype group then treat next highest allele as
			int32_t MaxReportGrpDGTs,     // when calling group DGTs then, if non-zero - report this many highest scoring DGTs
			int32_t MinDGTGrpMembers,        // groups with fewer than this number of members - in processing mode 5 only -are treated as noise alleles these groups are not used when determining DGT group specific major alleles
			double MinDGTGrpPropTotSamples,  // haplotype groups - in processing mode 5 only -containing less than this proportion of all samples are treated as if containing noise and alleles in these groups are not used when determining DGT group specific major alleles 
			double MinDGTFmeasure,           // only accepting DGT loci with at least this F-measure score
			int32_t KMerSize,					// use this KMer sized sequences when identifying group segregation hammings
			int32_t MinKMerHammings,			// must be at least this hamming separating any two groups
			int32_t KMerNoneCoverage,				// if > 0 then allow upto this number of bp within a KMer to have no coverage, maximum is 50% of KMerSize
			int32_t FndrTrim5,			// trim this many aligned PBAs from 5' end of founder aligned segments - reduces false alleles due to sequencing errors
			int32_t FndrTrim3,			// trim this many aligned PBAs from 3' end of founder aligned segments - reduces false alleles due to sequencing errors
			int32_t ProgTrim5,			// trim this many aligned PBAs from 5' end of progeny aligned segments - reduces false alleles due to sequencing errors
			int32_t ProgTrim3,			// trim this many aligned PBAs from 3' end of progeny aligned segments - reduces false alleles due to sequencing errors
			int32_t WWRLProxWindow,		// proximal window size for Wald-Wolfowitz runs test
			int32_t OutliersProxWindow,	// proximal window size for outliers reduction
			char *pszMaskBPAFile,	// optional masking input BPA file, only process BPAs which are intersect of these BPAs and progeny plus founder BPAs
			char* pszChromFile,			// BED file containing reference assembly chromosome names and sizes
			char* pszAlleleScoreFile,	// CSV file containing allelic association scores previously generated by processing in eMCSHSrcVsRefs or eMCSHRefsVsRefs mode
			int NumFounderInputFiles,	// number of input founder file specs
			char *pszFounderInputFiles[],	// names of input founder PBA files (wildcards allowed)
			int NumProgenyInputFiles,	// number of input progeny file specs
			char* pszProgenyInputFiles[],		// names of input progeny PBA files (wildcards allowed)
			char* pszOutFile,		// loci haplotype calls output file (CSV format)
			int	NumIncludeChroms,			// number of chromosome regular expressions to include
			char **ppszIncludeChroms,		// array of include chromosome regular expressions
			int	NumExcludeChroms,			// number of chromosome expressions to exclude
			char **ppszExcludeChroms,		// array of exclude chromosome regular expressions
			int NumThreads)		// number of worker threads to use
{
int Rslt;
CCallHaplotypes* pCallHaplotypes;

if((pCallHaplotypes = new CCallHaplotypes) == nullptr)
	{
	gDiagnostics.DiagOut (eDLFatal, gszProcName, "Unable to instantiate instance of CCallHaplotypes");
	return(eBSFerrObj);
	}
Rslt = pCallHaplotypes->Process(PMode,ExprID,SeedRowID,LimitPrimaryPBAs, AffineGapLen,GrpHapBinSize,NumHapGrpPhases,MinCentClustDist,MaxCentClustDist,MaxClustGrps, SparseRepPropGrpMbrs, SparseRepProp,
				MaxReportGrpDGTs,MinDGTGrpMembers,MinDGTGrpPropTotSamples,MinDGTFmeasure, KMerSize, MinKMerHammings,KMerNoneCoverage,
				FndrTrim5, FndrTrim3, ProgTrim5, ProgTrim3, WWRLProxWindow,OutliersProxWindow,
				pszMaskBPAFile,pszChromFile,pszAlleleScoreFile,NumFounderInputFiles,
				pszFounderInputFiles,NumProgenyInputFiles,pszProgenyInputFiles,pszOutFile,NumIncludeChroms,ppszIncludeChroms,NumExcludeChroms,ppszExcludeChroms,NumThreads);
delete pCallHaplotypes;
return(Rslt);
}


CCallHaplotypes::CCallHaplotypes()
{
m_pChromMetadata = nullptr;
m_pProgenyFndrAligns = nullptr;
m_pWorkQueueEls = nullptr;
m_pInBuffer = nullptr;
m_pszOutBuffer = nullptr;
m_pszWIGBuff = nullptr;
m_pAlleleStacks = nullptr;
m_pHGBinSpecs = nullptr;
m_pASBins = nullptr;
m_pASRefGrps = nullptr;
m_pGrpKMerSeqs = nullptr;
m_pDGTLoci = nullptr;
m_pKMerLoci = nullptr;
m_hOutFile = -1;
m_hWIGOutFile = -1;
m_hInFile = -1;
m_pCSVFile = nullptr;
m_pHGCSVFile = nullptr;
m_pAlleleScoresCSVFile = nullptr;

m_pBedFile = nullptr;
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
if(m_pWorkQueueEls != nullptr)
	delete []m_pWorkQueueEls;
if(m_pszOutBuffer != nullptr)
	delete []m_pszOutBuffer;
if(m_pszWIGBuff != nullptr)
	delete []m_pszWIGBuff;
if(m_pInBuffer != nullptr)
	delete []m_pInBuffer;

if(m_pCSVFile != nullptr)
	delete m_pCSVFile;
if(m_pHGCSVFile != nullptr)
    delete m_pHGCSVFile;
if (m_pAlleleScoresCSVFile != nullptr)
	delete m_pAlleleScoresCSVFile;
if (m_pBedFile != nullptr)
	delete m_pBedFile;
if(m_pChromMetadata != nullptr)
	{
	tsCHChromMetadata *pChromMetadata = m_pChromMetadata;
	for(int32_t Idx = 0; Idx < m_UsedNumChromMetadata; Idx++, pChromMetadata++)
		{
		if(pChromMetadata->pPBAs != nullptr)
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

if(m_pProgenyFndrAligns != nullptr)
	{
#ifdef _WIN32
	free(m_pProgenyFndrAligns);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pProgenyFndrAligns != MAP_FAILED)
		munmap(m_pProgenyFndrAligns, m_AllocdProgenyFndrAlignsMem);
#endif
	}

if(m_pAlleleStacks != nullptr)
{
#ifdef _WIN32
	free(m_pAlleleStacks);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pAlleleStacks != MAP_FAILED)
		munmap(m_pAlleleStacks, m_AllocdAlleleStacksMem);
#endif
	}

if (m_pGrpKMerSeqs != nullptr)
{
#ifdef _WIN32
	free(m_pGrpKMerSeqs);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (m_pGrpKMerSeqs != MAP_FAILED)
		munmap(m_pGrpKMerSeqs, m_AllocdKMerSeqMem);
#endif
}


if(m_pDGTLoci != nullptr)
	{
#ifdef _WIN32
	free(m_pDGTLoci);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pDGTLoci != MAP_FAILED)
		munmap(m_pDGTLoci, m_AllocdDGTLociMem);
#endif
	}

if (m_pKMerLoci != nullptr)
	{
#ifdef _WIN32
	free(m_pKMerLoci);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (m_pKMerLoci != MAP_FAILED)
		munmap(m_pKMerLoci, m_AllocdKMerLociMem);
#endif
	} 

if (m_pASBins != nullptr)
{
#ifdef _WIN32
	free(m_pASBins);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (m_pASBins != MAP_FAILED)
		munmap(m_pASBins, m_AllocdASBinsMem);
#endif
}

if (m_pASRefGrps != nullptr)
{
#ifdef _WIN32
	free(m_pASRefGrps);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (m_pASRefGrps != MAP_FAILED)
		munmap(m_pASRefGrps, m_AllocdASRefGrpMem);
#endif
}


if(m_pHGBinSpecs != nullptr)
{
	tsHGBinSpec* pBin;
	pBin = m_pHGBinSpecs;
	for(int32_t BinIdx = 0; BinIdx < m_UsedHGBinSpecs; BinIdx++, pBin++)
		if(pBin->pHaplotypeGroup != nullptr)
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

if(m_pCSVFile != nullptr)
	{
	delete m_pCSVFile;
	m_pCSVFile = nullptr;
	}

if(m_pHGCSVFile != nullptr)
	{
	delete m_pHGCSVFile;
	m_pHGCSVFile = nullptr;
	}

if (m_pAlleleScoresCSVFile != nullptr)
	{
	delete m_pAlleleScoresCSVFile;
	m_pAlleleScoresCSVFile = nullptr;
	}

if(m_pWorkQueueEls != nullptr)
	{
	delete []m_pWorkQueueEls;
	m_pWorkQueueEls = nullptr;
	}
m_AllocWorkQueueEls = 0;
m_TotWorkQueueEls = 0;
m_NumQueueElsProcessed = 0;
m_FastSerialise = 0;

if(m_pInBuffer != nullptr)
	{
	delete []m_pInBuffer;
	m_pInBuffer = nullptr;
	}

if(m_pszOutBuffer != nullptr)
	{
	delete []m_pszOutBuffer;
	m_pszOutBuffer = nullptr;
	}

if(m_pszWIGBuff != nullptr)
	{
	delete []m_pszWIGBuff;
	m_pszWIGBuff = nullptr;
	}

if (m_pBedFile != nullptr)
	{
	delete m_pBedFile;
	m_pBedFile = nullptr;
	}

if(m_pChromMetadata != nullptr)
	{
	tsCHChromMetadata *pChromMetadata = m_pChromMetadata;
	for(int32_t Idx = 0; Idx < m_UsedNumChromMetadata; Idx++, pChromMetadata++)
		{
		if(pChromMetadata->pPBAs != nullptr)
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
	m_pChromMetadata = nullptr;
	}
m_UsedNumChromMetadata = 0;
m_AllocdChromMetadata = 0;
m_AllocdChromMetadataMem = 0;

if(m_pAlleleStacks != nullptr)
	{
#ifdef _WIN32
	free(m_pAlleleStacks);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pAlleleStacks != MAP_FAILED)
		munmap(m_pAlleleStacks, m_AllocdAlleleStacksMem);
#endif
	m_pAlleleStacks = nullptr;
	}
m_AllocdAlleleStacksMem = 0;
m_UsedAlleleStacks = 0;
m_AllocdAlleleStacks = 0;

if(m_pHGBinSpecs != nullptr)
	{
	tsHGBinSpec* pBin;
	pBin = m_pHGBinSpecs;
	for(int32_t BinIdx = 0; BinIdx < m_UsedHGBinSpecs; BinIdx++, pBin++)
		if(pBin->pHaplotypeGroup != nullptr)
			delete[]pBin->pHaplotypeGroup;
#ifdef _WIN32
	free(m_pHGBinSpecs);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pHGBinSpecs != MAP_FAILED)
		munmap(m_pHGBinSpecs, m_AllocdHGBinSpecsMem);
#endif
	m_pHGBinSpecs = nullptr;
	}
m_UsedHGBinSpecs = 0;
m_AllocdHGBinSpecs = 0;
m_AllocdHGBinSpecsMem = 0;

if (m_pGrpKMerSeqs != nullptr)
	{
#ifdef _WIN32
	free(m_pGrpKMerSeqs);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (m_pGrpKMerSeqs != MAP_FAILED)
		munmap(m_pGrpKMerSeqs, m_AllocdKMerSeqMem);
#endif
	m_pGrpKMerSeqs = nullptr;
	}
m_UsedKMerSeqMem = 0;
m_AllocdKMerSeqMem = 0;

if(m_pDGTLoci != nullptr)
	{
#ifdef _WIN32
	free(m_pDGTLoci);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pDGTLoci != MAP_FAILED)
		munmap(m_pDGTLoci, m_AllocdDGTLociMem);
#endif
	m_pDGTLoci = nullptr;
	}
m_UsedDGTLoci = 0;
m_AllocdDGTLoci = 0;
m_AllocdDGTLociMem = 0;

if (m_pKMerLoci != nullptr)
	{
#ifdef _WIN32
	free(m_pKMerLoci);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (m_pKMerLoci != MAP_FAILED)
		munmap(m_pKMerLoci, m_AllocdKMerLociMem);
#endif
	m_pKMerLoci = nullptr;
	}
m_UsedKMerLoci = 0;
m_AllocdKMerLoci = 0;
m_AllocdKMerLociMem = 0;

if (m_pASBins != nullptr)
{
#ifdef _WIN32
	free(m_pASBins);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (m_pASBins != MAP_FAILED)
		munmap(m_pASBins, m_AllocdASBinsMem);
#endif
	m_pASBins = nullptr;
}
m_UsedASBins = 0;
m_AllocdASBins = 0;
m_AllocdASBinsMem = 0;

if (m_pASRefGrps != nullptr)
	{
#ifdef _WIN32
	free(m_pASRefGrps);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (m_pASRefGrps != MAP_FAILED)
		munmap(m_pASRefGrps, m_AllocdASRefGrpMem);
#endif
	m_pASRefGrps = nullptr;
	}
m_UsedASRefGrps = 0;
m_AllocASRefGrps = 0;
m_AllocdASRefGrpMem = 0;

DeleteMutexes();
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

m_NumFounders = 0;
m_NumProgenies = 0;

m_ExprID = 0;
m_InFileOfs = 0;
m_AffineGapLen = cDfltAffineGapLen;
m_FndrTrim5 = cDfltFndrTrim5;
m_FndrTrim3 = cDfltFndrTrim3;
m_ProgTrim5 = cDfltProgTrim5;
m_ProgTrim3 = cDfltProgTrim3;
m_GrpHapBinSize = cDfltGrpHapBinSize;
m_NumHapGrpPhases = cDfltNumHapGrpPhases;
m_MinCentClustDist = cDfltMinCentClustDist;
m_MaxCentClustDist = cDfltMaxCentClustDist;
m_MaxClustGrps = cDfltMaxClustGrps;
m_MaxReportGrpDGTs = cDfltMaxReportGrpDGTs;
m_SparseRepProp = cDfltSparseRepProp;
m_SparseRepPropGrpMbrs = cDfltSparseRepPropGrpMbrs;
m_MinDGTGrpMembers = cDfltMinDGTGrpMembers;
m_MinDGTGrpPropTotSamples = cDfltMinDGTGrpPropTotSamples;
m_MinDGTFmeasure = cDfltMinDGTFmeasure;
m_FbetameasureGrps = cDlftFbetaGrps;
m_WIGChromID = 0;
m_WIGRptdChromID = 0;
m_WIGSpanLoci = 0;
m_WIGSpanLen = 0;
m_WIGRptdSpanLen = 0;
m_WIGSpanCnts = 0;

m_InNumProcessed = 0;
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

memset(m_ChromSizes, 0, sizeof(m_ChromSizes));
m_TotChromSizes = 0;
m_NumChromSizes = 0;
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
if((m_hSerialiseAccess = CreateMutex(nullptr,false,nullptr))==nullptr)
	{
#else
if(pthread_mutex_init (&m_hSerialiseAccess,nullptr)!=0)
	{
#endif
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to create mutex");
	return(eBSFerrInternal);
	}
m_FastSerialise = 0;
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

// loading BED which specifies chrom names and sizes
int		// returning number of chromosomes parsed from BED file and accepted after filtering for wildcards
CCallHaplotypes::LoadChromSizes(char* pszBEDFile) // BED file containing chromosome names and sizes
{
int Rslt;
uint32_t ChromID;
int CurFeatureID;
char szFeatName[cMaxDatasetSpeciesChrom];
char szChromName[cMaxDatasetSpeciesChrom];
int32_t StartLoci;
int32_t EndLoci;
int NumChroms;
int64_t TotChromSizes;


if ((m_pBedFile = new CBEDfile) == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to instantiate CBEDfile");
	return(eBSFerrObj);
	}

if ((Rslt = m_pBedFile->Open(pszBEDFile, eBTAnyBed)) != eBSFSuccess)
	{
	while (m_pBedFile->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal, gszProcName, m_pBedFile->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to open '%s' for processing", pszBEDFile);
	return(eBSFerrOpnFile);
	}

m_TotChromSizes = 0;
m_NumChromSizes = 0;
NumChroms = 0;
TotChromSizes = 0;
CurFeatureID = 0;
while (Rslt == eBSFSuccess && (CurFeatureID = m_pBedFile->GetNextFeatureID(CurFeatureID)) > 0)
	{
	m_pBedFile->GetFeature(CurFeatureID,	// feature instance identifier
		szFeatName,				// where to return feature name
		szChromName,			// where to return chromosome name
		&StartLoci,				// where to return feature start on chromosome (0..n) 
		&EndLoci);				// where to return feature end on chromosome

	// is this chromosome to be be filtered out?
	if (!AcceptThisChromName(szChromName,false))
		{
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "LoadChromSizes: filtered out chromosome '%s'", szChromName);
		continue;
		}
 
	if (StartLoci != 0)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadChromSizes: unable to accept chromosome '%s' because start loci not 0", szChromName);
		delete m_pBedFile;
		m_pBedFile = nullptr;
		return(eBSFerrChrom);
		}

	if (EndLoci < 0 || (size_t)EndLoci >= cAllocPackedBaseAlleles)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadChromSizes: unable to accept chromosome '%s' because end loci not in acceptable range, must be in range 0..%zd", szChromName, cAllocPackedBaseAlleles);
		delete m_pBedFile;
		m_pBedFile = nullptr;
		return(eBSFerrChrom);
		}

	// expecting a single instance of each chromosome 
	if (LocateChrom(szChromName) != 0)	// must be a unique chromosome name
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadChromSizes: multiple instances of chromosome '%s', chromosome names must be unique", szChromName);
		delete m_pBedFile;
		m_pBedFile = nullptr;
		return(eBSFerrChrom);
		}
	// ensure not about to overflow m_ChromSizes[] capacity
	if(NumChroms == cMaxChromNames)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadChromSizes: chromosome '%s' - can only accept a maximum of %d from '%s'", szChromName, NumChroms, pszBEDFile);
		delete m_pBedFile;
		m_pBedFile = nullptr;
		return(eBSFerrChrom);
		}
	ChromID = AddChrom(szChromName);		// new chromosome
	m_ChromSizes[ChromID-1] = EndLoci + 1;	// record it's size (note that EndLoci is actually inclusive, so need to add 1 to obtain length!)
	TotChromSizes += (int64_t)EndLoci + 1;
	NumChroms += 1;
	}
m_NumChromSizes = NumChroms;
m_TotChromSizes = TotChromSizes;
gDiagnostics.DiagOut(eDLInfo, gszProcName, "LoadChromSizes: accepted %d chromosome sizes from '%s'", NumChroms, pszBEDFile);
delete m_pBedFile;
m_pBedFile = nullptr;
return(NumChroms);
}

int
CCallHaplotypes::Process(eModeCSH PMode,	// processing mode 0: report imputation haplotype matrix, 1: report both raw and imputation haplotype matrices, 2: additionally generate GWAS allowing visual comparisons, 3: allelic haplotype grouping,4: coverage haplotype grouping, 5: post-process haplotype groupings for DGTs,  6: post-process to WIG,  7 progeny PBAs vs. founder PBAs allelic association scores, 8 founder PBAs vs. founder PBAs allelic association scores, 9 grouping by scores, 10 differential group KMers
			int32_t ExprID,			         // assign this experiment identifier for this PBA analysis
			int32_t SeedRowID,               // generated CSVs will contain monotonically unique row identifiers seeded with this row identifier  
		    int32_t LimitPrimaryPBAs,        // limit number of loaded primary or founder PBA files to this many. 0: no limits, > 0 sets upper limit
			int32_t AffineGapLen,			 // haplotype grouping maximal gap length at which gaps treated same as mismatches (affine scoring); if < 0 then no affine scoring with complete gap scored, if 0 then gaps are not scored, if > 0 then gaps upto and including AffineGapLen are scored
			int32_t GrpHapBinSize,           // when grouping into haplotypes (processing mode 3) then use this non-overlapping bin size for accumulation of differentials between all founder PBAs
	        int32_t NumHapGrpPhases,         // number of phases over which to converge group consensus when haplotype clustering into groups - moves the balance between optimal group consensus and processing resource
			int32_t MinCentClustDist,        // haplotype groupings - in processing mode 3/4 only - minimum group centroid clustering distance
			int32_t MaxCentClustDist,        // haplotype groupings - in processing mode 3/4 only - maximum group centroid clustering distance
			int32_t MaxClustGrps,            // haplotype groupings - in processing mode 3 only - targeted maximum number of groups
			uint32_t SparseRepPropGrpMbrs,	// only apply sparse representative imputation if number of haplotype group members at least this 
			double SparseRepProp,			// if highest frequency (consensus) allele is 0x00 (no coverage) then if next highest frequency allele is at least this proportion of all members in haplotype group then treat next highest allele as
											// being as being the consensus. Objective to obtain an imputed allele in regions of sparse haplotype coverage
			int32_t MaxReportGrpDGTs,		// when calling group DGTs then, if non-zero - report this many highest scoring DGTs
			int32_t MinDGTGrpMembers,        // groups with fewer than this number of members - in processing mode 5 only -are treated as noise alleles these groups are not used when determining DGT group specific major alleles
			double MinDGTGrpPropTotSamples,  // haplotype groups - in processing mode 5 only -containing less than this proportion of all samples are treated as if containing noise and alleles in these groups are not used when determining DGT group specific major alleles 
			double MinDGTFmeasure,           // only accepting DGT loci with at least this F-measure score
			int32_t KMerSize,					// use this KMer sized sequences when identifying group segregation hammings
			int32_t MinKMerHammings,			// must be at least this hamming separating any two groups
			int32_t KMerNoneCoverage,				// if > 0 then allow upto this number of bp within a KMer to have no coverage, maximum is 50% of KMerSize
			int32_t FndrTrim5,				// trim this many aligned PBAs from 5' end of founder aligned segments - reduces false alleles due to sequencing errors
			int32_t FndrTrim3,				// trim this many aligned PBAs from 3' end of founder aligned segments - reduces false alleles due to sequencing errors
			int32_t ProgTrim5,				// trim this many aligned PBAs from 5' end of progeny aligned segments - reduces false alleles due to sequencing errors
			int32_t ProgTrim3,				// trim this many aligned PBAs from 3' end of progeny aligned segments - reduces false alleles due to sequencing errors
			int32_t WWRLProxWindow,			// proximal window size for Wald-Wolfowitz runs test
			int32_t OutliersProxWindow,		// proximal window size for outliers reduction
			char *pszMaskBPAFile,			// optional masking input BPA file, only process BPAs which are intersect of these BPAs and progeny plus founder BPAs
			char *pszChromFile,				// BED file containing reference assembly chromosome names and sizes
			char *pszAlleleScoreFile,		// CSV file containing allelic association scores previously generated by processing in eMCSHSrcVsRefs or eMCSHRefsVsRefs mode
			int NumFounderInputFiles,		// number of input founder file specs
			char *pszFounderInputFiles[],	// names of input founder PBA files (wildcards allowed)
			int NumProgenyInputFiles,		// number of input progeny file specs
			char* pszProgenyInputFiles[],	// names of input progeny PBA files (wildcards allowed)
			char* pszOutFile,				// loci haplotype calls output file (CSV format)
			int	NumIncludeChroms,			// number of chromosome regular expressions to include
			char **ppszIncludeChroms,		// array of include chromosome regular expressions
			int	NumExcludeChroms,			// number of chromosome expressions to exclude
			char **ppszExcludeChroms,		// array of exclude chromosome regular expressions
			int NumThreads)		// number of worker threads to use
{
int Rslt;
int NumFiles;
int TotNumFiles;
size_t memreq;

Reset();

CreateMutexes();
m_PMode = PMode;
m_ExprID = ExprID;
m_LimitPrimaryPBAs = LimitPrimaryPBAs;
m_bAllFndrsLociAligned = true; // replace with actual user parameter!
m_SeedRowID = SeedRowID;
m_CurRowID = SeedRowID;
m_MaxReportGrpDGTs = MaxReportGrpDGTs;
m_MinDGTGrpMembers = MinDGTGrpMembers;
m_MinDGTGrpPropTotSamples = MinDGTGrpPropTotSamples;
m_MinDGTFmeasure = MinDGTFmeasure;
m_GrpHapBinSize = GrpHapBinSize;
m_MinCentClustDist = MinCentClustDist;
m_MaxCentClustDist = MaxCentClustDist;
m_AffineGapLen = AffineGapLen;
m_MaxClustGrps = MaxClustGrps;
m_SparseRepProp = SparseRepProp;
m_SparseRepPropGrpMbrs = SparseRepPropGrpMbrs;
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
if(Rslt = (m_RegExprs.CompileREs(NumIncludeChroms, ppszIncludeChroms,NumExcludeChroms, ppszExcludeChroms)) < eBSFSuccess)
	{
	Reset();
	return(Rslt);
	}

if ((Rslt = LoadChromSizes(pszChromFile)) < 1) // BED file containing chromosome names and sizes - NOTE chromosomes will be filtered by include/exclude wildcards
	{
	Reset();
	return(Rslt);
	}

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

if(pszMaskBPAFile == nullptr || pszMaskBPAFile[0] == '\0')
	m_pszMaskBPAFile = nullptr;
else
	m_pszMaskBPAFile = pszMaskBPAFile;

if (pszAlleleScoreFile == nullptr || pszAlleleScoreFile[0] == '\0')
	m_pszAlleleScoreFile = nullptr;
else
	m_pszAlleleScoreFile = pszAlleleScoreFile;


if (PMode == eMCSHGroupScores)
{
	Rslt = GroupAlleleScores(pszAlleleScoreFile, pszMaskBPAFile,pszOutFile);
	Reset();
	return(Rslt);
}

if(PMode == eMCSHGrpDist2WIG)
	{
	Rslt = CSV2WIG(pszMaskBPAFile,pszOutFile);
	Reset();
	return(Rslt);
	}


// initial allocation, will be realloc'd as required if more memory required
memreq = (size_t)cAllocChromMetadata * sizeof(tsCHChromMetadata);
#ifdef _WIN32
m_pChromMetadata = (tsCHChromMetadata *)malloc(memreq);	// initial and perhaps the only allocation
if(m_pChromMetadata == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Initial memory allocation of %zd bytes for chromosome metadata failed - %s", (int64_t)memreq, strerror(errno));
	Reset();
	return(eBSFerrMem);
	}
#else
m_pChromMetadata = (tsCHChromMetadata *)mmap(nullptr, memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
if(m_pChromMetadata == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Memory allocation of %zd bytes through mmap() for chromosome metadata failed - %s", (int64_t)memreq, strerror(errno));
	m_pChromMetadata = nullptr;
	Reset();
	return(eBSFerrMem);
	}
#endif
m_AllocdChromMetadataMem = memreq;
m_UsedNumChromMetadata = 0;
m_AllocdChromMetadata = cAllocChromMetadata;

if(PMode == eMCSHDGTHapsGrps)
	{
	Rslt = ProcessGrpLociDGTs(MinDGTGrpMembers,MinDGTGrpPropTotSamples,MinDGTFmeasure,pszMaskBPAFile,NumFounderInputFiles,pszFounderInputFiles,pszOutFile);
	Reset();
	return(Rslt);
	}

if (PMode == eMCSHKMerGrps)
{
	Rslt = GenKMerGrpHammings(KMerSize, MinKMerHammings, KMerNoneCoverage, pszMaskBPAFile, NumFounderInputFiles, pszFounderInputFiles, pszOutFile);
	Reset();
	return(Rslt);
}

if(PMode < eMCSHAllelicHapsGrps)
	{
	memreq = (size_t)cAllocProgenyFndrAligns * sizeof(tsProgenyFndrAligns);
#ifdef _WIN32
	m_pProgenyFndrAligns = (tsProgenyFndrAligns*)malloc(memreq);	// initial and perhaps the only allocation
	if(m_pProgenyFndrAligns == nullptr)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Initial memory allocation of %zd bytes for chromosome windowed founder counts failed - %s",(int64_t)memreq, strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
#else
	m_pProgenyFndrAligns = (tsProgenyFndrAligns*)mmap(nullptr, memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	if(m_pProgenyFndrAligns == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Memory allocation of %zd bytes through mmap() for chromosome windowed founder counts failed - %s", (int64_t)memreq, strerror(errno));
		m_pProgenyFndrAligns = nullptr;
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
	if(m_pAlleleStacks == nullptr)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Initial memory allocation of %zd bytes for allele stacks failed - %s", (int64_t)memreq, strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
#else
	m_pAlleleStacks = (tsAlleleStack*)mmap(nullptr, memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	if(m_pAlleleStacks == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Memory allocation of %zd bytes through mmap() for allele stacks failed - %s", (int64_t)memreq, strerror(errno));
		m_pAlleleStacks = nullptr;
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
	m_pHGBinSpecs = nullptr;
	}
else
	{
	memreq = (size_t)cAllocHGBinSpecs * sizeof(tsHGBinSpec);
#ifdef _WIN32
	m_pHGBinSpecs = (tsHGBinSpec*)malloc(memreq);	// initial and perhaps the only allocation
	if(m_pHGBinSpecs == nullptr)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Initial memory allocation of %zd bytes for haplotype grouping bin definitions failed - %s", (int64_t)memreq, strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
#else
	m_pHGBinSpecs = (tsHGBinSpec*)mmap(nullptr, memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	if(m_pHGBinSpecs == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Memory allocation of %zd bytes through mmap() for haplotype grouping bin definitions failed - %s", (int64_t)memreq, strerror(errno));
		m_pHGBinSpecs = nullptr;
		Reset();
		return(eBSFerrMem);
		}
#endif
	m_AllocdHGBinSpecsMem = memreq;
	m_UsedHGBinSpecs = 0;
	m_AllocdHGBinSpecs = cAllocHGBinSpecs;
	m_pProgenyFndrAligns = nullptr;
	m_AllocdProgenyFndrAlignsMem = 0;
	m_UsedProgenyFndrAligns = 0;
	m_pAlleleStacks = nullptr;
	m_AllocdAlleleStacksMem = 0;
	m_UsedAlleleStacks = 0;
	m_AllocdAlleleStacks = 0;
	}

if((m_pInBuffer = new uint8_t[cInBuffSize]) == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Memory allocation of %u bytes for input buffering failed - %s", cInBuffSize, strerror(errno));
	Reset();
	return(eBSFerrMem);
	}
m_AllocInBuff = cInBuffSize;
m_InNumBuffered = 0;

gDiagnostics.DiagOut(eDLInfo, gszProcName, "Process: Loading metadata for founder pool ...");
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
			ReadsetID = LoadPBAFile(pszInFile, 0, true); // loading readset + chrom metadata only, later the PBAs for each individual chrom will be loaded as needed
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

m_NumProgenies = 0;
m_CurProgenyReadsetID = 0;
if (m_PMode == eMCSHRefsVsRefs)
	{
	if ((Rslt = GenPBAsHomozygosityScores(m_GrpHapBinSize, m_NumFounders, 0, pszOutFile)) < eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Failed processing progeny vs founder PBA files");
		Reset();
		return(Rslt);
		}
	Reset();
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Completed processing");
	return(Rslt);
	}

if(m_PMode == eMCSHSrcVsRefs)
	{
	// individually process each progeny against the founder panel

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

		if(m_NumProgenies + NumFiles > cMaxProgenyReadsets)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Can accept at most %d progeny PBA readsets for processing, after wildcard file name expansions there are %d requested", cMaxProgenyReadsets, TotNumFiles);
			Reset();
			return(eBSFerrMem);
			}

		Rslt = eBSFSuccess;
		for(int FileID = 0; Rslt >= eBSFSuccess && FileID < NumFiles; ++FileID)
			{
			pszInFile = glob.File(FileID);
			ReadsetID = LoadPBAFile(pszInFile, 1, true); // loading readset + chrom metadata only, later the PBAs for each individual chrom will be loaded as needed
			if(ReadsetID <= 0)
				{
				gDiagnostics.DiagOut(eDLInfo, gszProcName, "Process: Errors loading file '%s'",pszInFile);
				Reset();
				return(ReadsetID);
				}
			m_NumProgenies++;
			}
		}

	if ((Rslt = GenPBAsHomozygosityScores(m_GrpHapBinSize, m_NumFounders, m_NumProgenies, pszOutFile)) < eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Failed processing progeny vs founder PBA files");
		Reset();
		return(Rslt);
		}
	Reset();
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Completed processing");
	return(Rslt);
	}

if(!(m_PMode == eMCSHAllelicHapsGrps || m_PMode == eMCSHCoverageHapsGrps || m_PMode == eMCSHDGTHapsGrps))
	{
	// load the control PBA file if it has been specified
	m_MaskReadsetID = 0;
	if(m_pszMaskBPAFile != nullptr)
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
			{
			m_CurRowID = m_SeedRowID;
			ReportHaplotypesByProgeny(pszOutFile, m_ProgenyIDs[ProgIdx]);
			}
		}

	// generate a matrix with chrom.loci (Y) by progeny (X)
	m_mtqsort.qsort(m_pProgenyFndrAligns, (int64_t)m_UsedProgenyFndrAligns, sizeof(tsProgenyFndrAligns), SortProgenyChromLociReadset);
	ReportMatrix(pszOutFile);
	}
else // else (eMCSHAllelicHapsGrps or eMCSHCoverageHapsGrps or eMCSHDGTHapsGrps or eMCSHHyperConserved)
	{
	// optional Haplotype group clustering specifications file(s) or previously generated haplotype groupings file
	if(m_pszMaskBPAFile != nullptr)
		{
		glob.Init();
		if(glob.Add(m_pszMaskBPAFile) < SG_SUCCESS)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to glob haplotype group clustering specification file(s) '%s", m_pszMaskBPAFile);
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

			if(m_UsedHGBinSpecs == 0)
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to parse out any haplotype group clustering specifications from file(s) matching '%s", m_pszMaskBPAFile);
				Reset();
				return(eBSFerrOpnFile);	// treat as though unable to open file
				}

			if(m_UsedHGBinSpecs > 1)
				m_mtqsort.qsort(m_pHGBinSpecs, (int64_t)m_UsedHGBinSpecs, sizeof(tsHGBinSpec), SortHGBinSpecs); // sort ascending by ChromID.StartLoci.NumLoci.Distance ascending

			int32_t CurChromID;
			int32_t HGBinID;
			tsHGBinSpec* pHGBinSpec;
			tsCHChromMetadata* pChromMetadata;
			CurChromID = 0;
			pHGBinSpec = m_pHGBinSpecs;
			for(HGBinID = 1; HGBinID <= m_UsedHGBinSpecs; HGBinID++, pHGBinSpec++)
				{
				pHGBinSpec->BinID = HGBinID;
				if(pHGBinSpec->ChromID != CurChromID)
					{
					CurChromID = pHGBinSpec->ChromID;
					pChromMetadata = LocateChromMetadataFor(1, pHGBinSpec->ChromID);
					pChromMetadata->HGBinID = HGBinID; // first bin on this chromosome 
					}
				}
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
uint32_t NumRuns;
uint8_t ChkHap;
uint8_t PrevHap;
uint32_t NumFa;
uint32_t NumFb;
tsProgenyFndrAligns *pCurPFA;
tsProgenyFndrAligns *pChkPFA;
size_t PFAIdx;
size_t ChkHapStart;
size_t ChkHapEnd;

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
pChkPFA = nullptr;
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
	ChkHapStart = max((size_t)0,(size_t)PFAIdx-9);
	ChkHapEnd = min((int64_t)m_UsedProgenyFndrAligns-1,(int64_t)PFAIdx+10);
	pChkPFA = &m_pProgenyFndrAligns[ChkHapStart];
	for(; ChkHapStart <= ChkHapEnd; ChkHapStart++,pChkPFA++)
		{
		if(pChkPFA->ReadsetID != CurReadsetID)
			continue;
		if(pChkPFA->NumProgenyFounders == 0)
			continue;
		if(abs((int32_t)pCurPFA->Loci - (int32_t)pChkPFA->Loci) > MaxDistance)
			continue;

		// must have been at least 1 founder for this progeny..
		// count numbers of each founder
		if(pChkPFA->NumProgenyFounders == 1)
			{
			// only a single founder but which founder?
			ChkHap = 0;
			for(int FndrIdx = 0; ChkHap == 0 && FndrIdx < cMaxBitVectBits-1; FndrIdx++)
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
			for(int FndrIdx = 0; ChkHap < 2 && FndrIdx < cMaxBitVectBits-1; FndrIdx++)
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

pPrevPFA = nullptr;
pNxtPFA = nullptr;
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
	pNxtPFA = nullptr;
	pPrevPFA = nullptr;
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
		if(pPrevPFA != nullptr && bPrevEqual == false)
			{
			pCurPFA->NumProgenyFounders = pPrevPFA->NumProgenyFounders;
			pCurPFA->ProgenyFounders = pPrevPFA->ProgenyFounders;
			}
		}
	else
		{
		if(pNxtPFA != nullptr && bNxtEqual == false)
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
tsBitsVect CurFounders;
tsBitsVect ExtdFounders;
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

pszProgenyReadset = LocateReadset(ReadsetID);
if(ReadsetID == 0)
	sprintf(szOutFile, "%s.progeny.%u.%s.all.csv", pszRsltsFileBaseName, m_ExprID,bRaw ? "raw" : "imputed");
else
	sprintf(szOutFile, "%s.progeny.%u.%s.%s.csv", pszRsltsFileBaseName, m_ExprID, pszProgenyReadset, bRaw ? "raw" : "imputed");
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
m_OutBuffIdx = sprintf((char*)m_pszOutBuffer, "\"ExprID\",\"RowID\",\"Progeny\",\"Chrom\",\"Loci\"");
for(FndrIdx = 0; FndrIdx < m_NumFounders; FndrIdx++)
	{
	if(m_Fndrs2Proc[FndrIdx] & 0x01)
		m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], ",\"Fndr:%s\"", LocateReadset(FndrIdx+1));
	}
m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], "\n");

m_CurRowID = m_SeedRowID;
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

	m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx],"%u,%u,\"%s\",\"%s\",%u",m_ExprID,m_CurRowID++, pszProgenyReadset, pszChrom, pCurPFA->Loci);
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

if(m_pszOutBuffer == nullptr)
	{
	if((m_pszOutBuffer = new uint8_t[cOutBuffSize]) == nullptr)
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
int32_t FndrIdx;
int32_t NumFndrs;
size_t ASIdx;
int32_t AlleleIdx;


if(m_pszOutBuffer == nullptr)
	{
	if((m_pszOutBuffer = new uint8_t[cOutBuffSize]) == nullptr)
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
m_OutBuffIdx = sprintf((char*)m_pszOutBuffer, "\"ExprID\",\"RowID\",\"Chrom\",\"Loci\",\"BinAlleles\",\"Alleles\",\"FndrsAlleleA\",\"FndrsAlleleC\",\"FndrsAlleleG\",\"FndrsAlleleT\"");

int32_t PrevReadsetID = 0;
int32_t PrevChromID = 0;
pCurAlleleStack = &m_pAlleleStacks[0];
m_CurRowID = m_SeedRowID;
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
	m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], "\n%u,%u,\"%s\",%u,%u,\"%s\"",m_ExprID,m_CurRowID++, pszChrom, pCurAlleleStack->Loci, BinAlleles, szAlleles);
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
uint32_t Fndr;
size_t PFAIdx;

if(m_pszOutBuffer == nullptr)
	{
	if((m_pszOutBuffer = new uint8_t[cOutBuffSize]) == nullptr)
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

if((pszFounders = new char[cMaxFounderReadsets*50]) == nullptr)
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
uint8_t LoMinorAllele;
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
			if ((LoMinorAllele = (*pAlleles & ~AlleleMsk)) != 0) // can only be one major, but could be a very low cover minor ...
				{
				if (LoMinorAllele == 0x01 || LoMinorAllele == 0x04 || LoMinorAllele == 0x10 || LoMinorAllele == 0x40) // when low cover minor then remove that minor
					{
					*pAlleles = AlleleMsk;
					return(true);
					}
				if (bSetNoAlleles)
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
int32_t ProgHaplotype;
int32_t CurChromID;
int32_t CurLoci;
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

sprintf(szOutFile, "%s.%d.%s.matrix.csv", pszRsltsFileBaseName,m_ExprID,bRaw ? "raw" : "inputed");
gDiagnostics.DiagOut(eDLInfo, gszProcName, "Reporting matrix as CSV to file '%s'", szOutFile);

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
m_OutBuffIdx = sprintf((char*)m_pszOutBuffer, "\"ExprID\",\"RowID\",\"Chrom\",\"Loci\"");

for(ProgIdx = 0; ProgIdx < m_NumProgenies; ProgIdx++)
	m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], ",\"Progeny:%s\"", LocateReadset(m_ProgenyIDs[ProgIdx]));

bNewLoci = false;
bProgFndrs = false;
memset(ProgenyHaplotypes,-1, sizeof(ProgenyHaplotypes));
CurLoci = 0;
CurChromID = 0;
m_CurRowID = m_SeedRowID;
pCurPFA = &m_pProgenyFndrAligns[0];
for(PFAIdx = 0; PFAIdx < m_UsedProgenyFndrAligns; PFAIdx++, pCurPFA++)
	{
	if(pCurPFA->ChromID != CurChromID || pCurPFA->Loci != CurLoci)
		{
		if(bNewLoci && bProgFndrs)
			{
			m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], "\n%u,%u,\"%s\",%u", m_ExprID,m_CurRowID++,pszChrom, CurLoci);
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




#ifdef _WIN32
unsigned __stdcall WorkerAlignSelfPBAsInstance(void* pThreadPars)
#else
void* WorkerAlignSelfPBAsInstance(void* pThreadPars)
#endif
{
	int Rslt;
	tsCHWorkerSelfScoreInstance* pPars = (tsCHWorkerSelfScoreInstance*)pThreadPars;			// makes it easier not having to deal with casts!
	CCallHaplotypes* pWorkerInstance = (CCallHaplotypes*)pPars->pThis;

	Rslt = pWorkerInstance->AlignSelfPBAsThread(pPars);
	pPars->Rslt = Rslt;
#ifdef _WIN32
	_endthreadex(0);
	return(eBSFSuccess);
#else
	pthread_exit(&pPars->Rslt);
#endif
}

int 
CCallHaplotypes::AlignSelfPBAs(int32_t NumRefPBAs,int32_t NumSrcPBAs, int32_t ChromID, int32_t ChromLen, int32_t BinsThisChrom, int32_t MaxBinSize, tsCHChromScores *pChromScores,uint8_t *pFndrPBAs[],int MaxThreads)
{
int Rslt = -1;
tsCHWorkerSelfScoreInstance CHWorkerSelfScoreInstances[cMaxPBAWorkerThreads];
int NumThreads;
tsCHWorkerSelfScoreInstance *pThreadPar;
int32_t ThreadIdx;

#ifndef _WIN32
// increase the default thread stack to at least cWorkThreadStackSize
size_t defaultStackSize;
pthread_attr_t threadattr;
pthread_attr_init(&threadattr);
pthread_attr_getstacksize(&threadattr, &defaultStackSize);
if (defaultStackSize < (size_t)cWorkThreadStackSize)
	{
	if ((Rslt = pthread_attr_setstacksize(&threadattr, (size_t)cWorkThreadStackSize)) != 0)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AlignSelfPBAs: pthread_attr_setstacksize(%d) failed, default was %zd", cWorkThreadStackSize, defaultStackSize);
		return(eBSFerrInternal);
		}
	}
#endif

NumThreads = min(NumSrcPBAs > 0 ? NumSrcPBAs : NumRefPBAs, MaxThreads);
m_CurSelfID = 0;
m_NumWorkerInsts = 0;
m_CompletedWorkerInsts = 0;
m_ExpNumWorkerInsts = NumThreads;
pThreadPar = CHWorkerSelfScoreInstances;
for (ThreadIdx = 1; ThreadIdx <= NumThreads; ThreadIdx++, pThreadPar++)
	{
	memset(pThreadPar, 0, sizeof(tsCHWorkerSelfScoreInstance));
#ifdef _WIN32
	pThreadPar->threadHandle = nullptr;
#else
	pThreadPar->threadID = 0;
#endif
	pThreadPar->NumRefPBAs = NumRefPBAs;
	pThreadPar->NumSrcPBAs = NumSrcPBAs;
	pThreadPar->ChromID = ChromID;
	pThreadPar->ChromLen = ChromLen;
	pThreadPar->BinsThisChrom = BinsThisChrom;
	pThreadPar->MaxBinSize = MaxBinSize;
	pThreadPar->pChromScores = pChromScores;
	pThreadPar->pFndrPBAs = pFndrPBAs;
	pThreadPar->ThreadIdx = ThreadIdx;
	pThreadPar->pThis = this;
#ifdef _WIN32
	pThreadPar->threadHandle = (HANDLE)_beginthreadex(nullptr, cWorkThreadStackSize, WorkerAlignSelfPBAsInstance, pThreadPar, 0, &pThreadPar->threadID);
#else
	pThreadPar->threadRslt = pthread_create(&pThreadPar->threadID, &threadattr, WorkerAlignSelfPBAsInstance, pThreadPar);
#endif
	}

// allow threads time to all startup
// check if all threads did actually startup; if they did so then m_NumWorkerInsts will have been incremented to NumInstances
int MaxWait = cMaxWaitThreadsStartup;		// allowing at most cMaxWaitThreadsStartup secs for all threads to startup, polling every second
int StartedInstances;
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
} while (StartedInstances != NumThreads && MaxWait > 0);

if(StartedInstances == NumThreads)				// all started up?
	while(WaitWorkerThreadStatus(1000) != 0);	// returns 0 if all threads started and completed processing, 1 if all threads started but some yet to complete within WaitSecs, 2 if not all threads have started within WaitSecs

pThreadPar = CHWorkerSelfScoreInstances;
for (ThreadIdx = 1; ThreadIdx <= NumThreads; ThreadIdx++, pThreadPar++)
	{
#ifdef _WIN32
	if(pThreadPar->threadHandle == nullptr)
		continue;
	while (WAIT_TIMEOUT == WaitForSingleObject(pThreadPar->threadHandle, 1000))
		{
		}
	CloseHandle(pThreadPar->threadHandle);
	pThreadPar->threadHandle = nullptr;
#else
	struct timespec ts;
	int JoinRlt;
	clock_gettime(CLOCK_REALTIME, &ts);
	ts.tv_sec += 1;
	while ((JoinRlt = pthread_timedjoin_np(pThreadPar->threadID, nullptr, &ts)) != 0)
		{
		ts.tv_sec += 1;
		}

#endif
	if (pThreadPar->Rslt != eBSFSuccess)
		Rslt = pThreadPar->Rslt;
	}
return(Rslt);
}


// thread for scoring an instance of a source PBAs (termed as Self) against all founder PBAs
int
CCallHaplotypes::AlignSelfPBAsThread(tsCHWorkerSelfScoreInstance* pPar)
{
int Rslt;
int BinID;
int EndBinLoci;
int32_t SrcID;
int32_t RefID;
int32_t Loci;
uint8_t SrcPBA;
uint8_t RefPBA;
uint8_t* pSrcPBA;
uint8_t* pSrcLoci;
uint8_t* pRefLoci;
uint32_t ReqTerminate;
int32_t SrcPBAsIdx;
int32_t MaxNumSrcs;
tsCHChromScores* pChromScore;

// one more thread instance has started
AcquireSerialise();
m_NumWorkerInsts++;
ReleaseSerialise();

SrcPBAsIdx = pPar->NumSrcPBAs > 0 ? pPar->NumRefPBAs : 0;		// PBAs layout is such that references start at index 0, followed by sources at index NumRefPBAs, if All vs All
MaxNumSrcs = pPar->NumSrcPBAs > 0 ? pPar->NumSrcPBAs : pPar->NumRefPBAs; // sources are same as references so indexes will be identical
if(pPar->MaxBinSize > pPar->ChromLen)
	pPar->MaxBinSize = pPar->ChromLen;

// iterating over each unprocessed (by either this or another thread) source
for (SrcID = 1; SrcID <= MaxNumSrcs; SrcID++, SrcPBAsIdx++)
	{
	// check if requested to early terminate
	AcquireSerialise();
	ReqTerminate = m_ReqTerminate;
	ReleaseSerialise();
	if (ReqTerminate)
		break;

	// only process a source if no other thread already started to process that source
	AcquireFastSerialise();
	if (SrcID <= m_CurSelfID)	// ensuring source is unprocessed 
		{
		ReleaseFastSerialise();
		continue;				// too late, another thread is processing or completed this source!
		}
	m_CurSelfID = SrcID;		// let other threads know this SrcID is now being processed
	ReleaseFastSerialise();

	// iterate over each of the the references and score the current source relative to these references
	for (RefID = 1; RefID <= pPar->NumRefPBAs; RefID++)
		{
		// starting from first bin at loci 0 on the source
		pChromScore = &pPar->pChromScores[(pPar->NumRefPBAs * pPar->BinsThisChrom * (SrcID - 1)) + ((RefID - 1) * pPar->BinsThisChrom)];

		if(pPar->pFndrPBAs[SrcPBAsIdx] == nullptr || pPar->pFndrPBAs[RefID - 1] == nullptr) // explicitly handle these cases, some sample chromosomes will have no PBAs. 
			{
			for(BinID = 1, Loci = 0; BinID <= pPar->BinsThisChrom && Loci < pPar->ChromLen; BinID++, Loci += pPar->MaxBinSize,pChromScore++)
				{
				memset(pChromScore, 0, sizeof(tsCHChromScores));
				pChromScore->ChromID = pPar->ChromID;
				pChromScore->SrcID = SrcID;
				pChromScore->RefID = RefID;
				pChromScore->BinID = BinID;
				pChromScore->BinLoci = Loci;
				pChromScore->BinSize = (pPar->ChromLen - Loci) > pPar->MaxBinSize ? pPar->MaxBinSize : pPar->ChromLen - Loci;
				pChromScore->ExactScore = 0.0;
				pChromScore->PartialScore = 0.0;
				pChromScore->BinSize = (pPar->ChromLen - Loci) > pPar->MaxBinSize ? pPar->MaxBinSize : pPar->ChromLen - Loci;
				}
			continue;
			}

		pSrcPBA = pPar->pFndrPBAs[SrcPBAsIdx]; // already handled cases of nullptr; unlikely but possible if sample has no PBAs on current chromosome
		pSrcLoci = pSrcPBA;
		pRefLoci = pPar->pFndrPBAs[RefID - 1];
		EndBinLoci = pPar->MaxBinSize - 1;
		BinID = 1;
		for (Loci = 0; Loci < pPar->ChromLen && BinID <= pPar->BinsThisChrom; Loci++)
			{
			if (Loci == 0 || Loci > EndBinLoci)			// initial bin or starting a new bin
				{
				if(Loci > 0)		// Loci > 0 if starting a new bin
					{
					// score just processed bin before starting on new bin
					if(pChromScore->AlignLen > 0)
						{
						// fixed down weighting of 0.5 applied to partial matchings or a match where one allele matched but there was a non-ref allele also present
						pChromScore->PartialScore = (double)(pChromScore->NumExactMatches + (pChromScore->NumPartialMatches + pChromScore->NumNonRefAlleles) / 2) / pChromScore->AlignLen;
						// only scoring exact matches
						pChromScore->ExactScore = (double)pChromScore->NumExactMatches / pChromScore->AlignLen;
						}
					else
						{
						pChromScore->ExactScore = 0.0;
						pChromScore->PartialScore = 0.0;
						}

					pChromScore++;
					EndBinLoci += pPar->MaxBinSize;
					BinID++;
					}
				memset(pChromScore,0,sizeof(tsCHChromScores));
				pChromScore->ChromID = pPar->ChromID;
				pChromScore->SrcID = SrcID;
				pChromScore->RefID = RefID;
				pChromScore->BinID = BinID;
				pChromScore->BinLoci = Loci;
				pChromScore->BinSize = (pPar->ChromLen - Loci) > pPar->MaxBinSize ? pPar->MaxBinSize : pPar->ChromLen - Loci;
				}

			SrcPBA = *pSrcLoci++;
			RefPBA = *pRefLoci++;
			if (RefPBA == 0 || SrcPBA == 0) // both source and reference must have coverage at Loci, if not then iterate to next loci
				continue;

			// enabling partial matching on the basis that fragments are being sequenced. In a diploid both haplotypes at a given loci may not be present after sequencing, especially where there is low coverage - WGS skim reads - or  GBS reads
			// when evaluating same material GBS against WGS then about 40% of the WGS bialleles are present as monoalleles in the GBS 
			if (SrcPBA == RefPBA)	// alleles exactly matching?
				{
				pChromScore->NumExactMatches++;
				if (RefPBA == 0xf0 || RefPBA == 0xcc || RefPBA == 0xc3 || RefPBA == 0x3c || RefPBA == 0x33 || RefPBA == 0x0f)
					pChromScore->NumBiallelicExactMatches++;
				}
			else
				if (SrcPBA & RefPBA) // at least one of the reference alleles present?
					{
					if (~RefPBA & SrcPBA)
						pChromScore->NumNonRefAlleles++;
					else
						pChromScore->NumPartialMatches++;
					}
			pChromScore->AlignLen++;
			}
		// last bin still requires scoring
		if (pChromScore->AlignLen > 0)
			{
			// fixed down weighting of 0.5 applied to partial matchings or a match where one allele matched but there was a non-ref allele also present
			pChromScore->PartialScore = (double)(pChromScore->NumExactMatches + (pChromScore->NumPartialMatches + pChromScore->NumNonRefAlleles) / 2) / pChromScore->AlignLen;

				// only scoring exact matches
			pChromScore->ExactScore = (double)pChromScore->NumExactMatches / pChromScore->AlignLen;
			}
		else
			{
			pChromScore->ExactScore = 0.0;
			pChromScore->PartialScore = 0.0;
			}
		}
	}
if (SrcID <= MaxNumSrcs)
	Rslt = -1;
else
	Rslt = 0;
AcquireSerialise();
m_CompletedWorkerInsts++;
ReleaseSerialise();
return(Rslt);
}



int
CCallHaplotypes::GenPBAsHomozygosityScores(int MaxBinSize,					// binning using these maximal sized bins accross each chromosome, 0 if 1 bin per chromosome
											int32_t NumRefPBAs,				// number of reference PBAs to be processed
											int32_t NumSrcPBAs,				// number of source PBAs to be processed, 0 if none and processing will be reference vs. reference
											char* pszRsltsFileBaseName)		// results are written to this file base name
{
uint8_t** pFndrPBAs;
tsCHChromScores* pChromScores;
tsCHChromScores* pChromScore;

int32_t PBAsID;
int32_t SrcID;

int32_t InitialSrcID;

tsCHReadsetMetadata* pFndrMetadata;
tsCHChromMetadata* pChromMetadata;
uint32_t CurChromID;
int32_t CurChromMetadataIdx;
int32_t ChromIdx;
char* pszChrom;
char szOutScore[_MAX_PATH];

FILE* pOutScoreFile;


// allocate for tsASBin's if not already allocated
if (m_pASBins == nullptr)
	{
	size_t memreq = (size_t)cAllocASBins * sizeof(tsASBin);
#ifdef _WIN32
	m_pASBins = (tsASBin*)malloc(memreq);	// initial and perhaps the only allocation
	if (m_pASBins == nullptr)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Initial memory allocation of %zd bytes for allele score bins failed - %s", (int64_t)memreq, strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
#else
	m_pASBins = (tsASBin*)mmap(nullptr, memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	if (m_pASBins == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Memory allocation of %zd bytes through mmap() for allele score bins failed - %s", (int64_t)memreq, strerror(errno));
		m_pASBins = nullptr;
		Reset();
		return(eBSFerrMem);
		}
#endif
	m_AllocdASBinsMem = memreq;
	m_AllocdASBins = cAllocASBins;
	}
m_UsedASBins = 0;

sprintf(szOutScore, "%s%s", pszRsltsFileBaseName, ".score.csv");

if ((pOutScoreFile = fopen(szOutScore, "w")) == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenPBAsHomozygosityScores: Unable to create/truncate file %s for writing error: %s", szOutScore, strerror(errno));
	Reset();
	return(eBSFerrOpnFile);
	}

int BinsThisChrom = 0;

int32_t ScoreMatrixSize = (NumRefPBAs * (NumSrcPBAs > 0 ? NumSrcPBAs : NumRefPBAs)); // NumSrcPBAs == 0 if all Refs scored against all other Refs
int32_t ScoreChromBinsMatrixSize = 0;

InitialSrcID = NumSrcPBAs > 0 ? NumRefPBAs + 1 : 1;

pFndrPBAs = new uint8_t * [cMaxFounderReadsets + cMaxProgenyReadsets];

pChromScores = nullptr;

char* pszLineBuff = new char[cOutBuffSize];
int32_t LineOfs;

// header row contains the names of those column reference PBAs against which the source PBAs were aligned and scored for homozygosity
LineOfs = sprintf(pszLineBuff, "\"SourcePBA\",\"ReferencePBA\",\"Chrom\",\"Bin\",\"BinLoci\",\"BinSize\",\"AlignLen\",\"NumExactMatches\",\"NumBiallelicExactMatches\",\"NumPartialMatches\",\"NumNonRefAlleles\",\"ExactScore\",\"PartialScore\"\n");
fwrite(pszLineBuff, 1, LineOfs, pOutScoreFile);
fflush(pOutScoreFile);
LineOfs = 0;

gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenPBAsHomozygosityScores: Starting to process %d PBAs against %d reference PBAs", NumSrcPBAs > 0 ? NumSrcPBAs : NumRefPBAs, NumRefPBAs);
pFndrMetadata = &m_Readsets[0];		// references were loaded first so 1st reference will be here
CurChromMetadataIdx = pFndrMetadata->StartChromMetadataIdx;
for (ChromIdx = 0; ChromIdx < pFndrMetadata->NumChroms && CurChromMetadataIdx != 0; ChromIdx++)
	{
	pChromMetadata = &m_pChromMetadata[CurChromMetadataIdx - 1];
	CurChromMetadataIdx = pChromMetadata->NxtChromMetadataIdx;
	CurChromID = pChromMetadata->ChromID;
	pszChrom = LocateChrom(CurChromID);

	if(MaxBinSize == 0)
		BinsThisChrom = 1;
	else
		BinsThisChrom = 1 + ((m_ChromSizes[CurChromID-1] - 1) / MaxBinSize);


	ScoreChromBinsMatrixSize = ScoreMatrixSize * BinsThisChrom; //  required allocation for this chromosome
	if(pChromScores != nullptr)
		delete []pChromScores;
	pChromScores = new tsCHChromScores[ScoreChromBinsMatrixSize];
	memset(pChromScores, 0, sizeof(tsCHChromScores) * ScoreMatrixSize * BinsThisChrom);

	// ensure any previously loaded PBAs for current chromosome are deleted
	if (ChromIdx > 0)
		{
		for (PBAsID = 1; PBAsID <= (NumRefPBAs + NumSrcPBAs); PBAsID++)
			DeleteSampleChromPBAs(PBAsID, pChromMetadata->ChromID);
		}



	// load all PBAs (both references and source PBAs) for current chrom
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenPBAsHomozygosityScores: Loading chromosome PBAs for '%s' from %d PBAs", pszChrom, NumRefPBAs + NumSrcPBAs);
	int WorkerThreadStatus;
	WorkerThreadStatus = StartWorkerLoadChromPBAThreads(m_NumThreads, 1, (NumRefPBAs + NumSrcPBAs), CurChromID, true);
	if (WorkerThreadStatus < 0)
		{
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenPBAsHomozygosityScores: Errors loading chromosome PBAs for '%s' from %d PBAs ....", pszChrom, NumRefPBAs + NumSrcPBAs);
		return(WorkerThreadStatus);
		}

	while (WorkerThreadStatus > 0)
		{
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenPBAsHomozygosityScores: Continuing to load chromosome PBAs for '%s' from %d PBAs ....", pszChrom, NumRefPBAs + NumSrcPBAs);
		WorkerThreadStatus = WaitWorkerThreadStatus(60);
		}
	TerminateLoadChromPBAsThreads(60);
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenPBAsHomozygosityScores: Completed loading of chromosome PBAs for '%s' from %d PBAs", pszChrom, NumRefPBAs + NumSrcPBAs);

	// checking that all PBAs do have current chrom sequences
	for (PBAsID = 1; PBAsID <= (NumRefPBAs + NumSrcPBAs); PBAsID++)
		{
		tsCHChromMetadata* pChromMetadata;

		if ((pChromMetadata = LocateChromMetadataFor(PBAsID, CurChromID)) == nullptr)
			{
			pFndrPBAs[PBAsID - 1] = nullptr;				// nullptr but subsequent processing will check and handle appropriately
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenPBAsHomozygosityScores: No metadata for chromosome '%s' in this PBA '%s', treating as if chromosome loaded but had no alignments ...", pszChrom, LocateReadset(PBAsID));
			}
		else
			{
			pFndrPBAs[PBAsID - 1] = pChromMetadata->pPBAs; // could be assigning nullptr if no chromosome reads but subsequent processing will check and handle appropriately
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenPBAsHomozygosityScores: Metadata for chromosome '%s' exists in this PBA '%s', but had no alignments ...", pszChrom, LocateReadset(PBAsID));
			}
		}

	// all founder PBAs for current chromosome now loaded and normalised for coverage
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenPBAsHomozygosityScores: Completed loading chromosome PBAs for '%s' from %d PBAs, starting alignment scoring ...", pszChrom, NumRefPBAs + NumSrcPBAs);

	int32_t ActualMaxBinSize = MaxBinSize == 0 ? pChromMetadata->ChromLen : MaxBinSize;
	AlignSelfPBAs(NumRefPBAs, NumSrcPBAs, CurChromID, pChromMetadata->ChromLen, BinsThisChrom, ActualMaxBinSize, pChromScores, pFndrPBAs, m_NumThreads);
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenPBAsHomozygosityScores: Completed scoring chromosome PBAs for '%s' from %d PBAs", pszChrom, NumRefPBAs + NumSrcPBAs);

	pChromScore = pChromScores;
	tsCHChromScores* pChromScoreA;
	char *pszRefReadset;
	char *pszSrcReadset;

	int32_t BinIdx;
	int32_t SrcIdx;
	int32_t RefIdx;
	int32_t NumProcSrcPBAs = NumSrcPBAs > 0 ? NumSrcPBAs : NumRefPBAs;

	pChromScore = pChromScores;
	SrcID = NumSrcPBAs > 0 ? NumRefPBAs + 1: 1;
	for (SrcIdx = 0; SrcIdx < NumProcSrcPBAs; SrcIdx++, SrcID++)
		{
		pChromScoreA = &pChromScores[SrcIdx * NumRefPBAs * BinsThisChrom]; // first score for given SrcIdx
		pszSrcReadset = LocateReadset(SrcID);
		for (RefIdx = 0; RefIdx < NumRefPBAs; RefIdx++)
			{
			pszRefReadset = LocateReadset(RefIdx+1);
			int64_t SumRefLengths = 0;
			for (BinIdx = 0; BinIdx < BinsThisChrom; BinIdx++, pChromScore++)
				{
				LineOfs += sprintf(&pszLineBuff[LineOfs], "\"%s\",\"%s\",\"%s\",%d,%d,%d,%d,%d,%d,%d,%d,%.7f,%.7f\n", pszSrcReadset, pszRefReadset,pszChrom, pChromScore->BinID,pChromScore->BinLoci, pChromScore->BinSize,
													pChromScore->AlignLen, pChromScore->NumExactMatches, pChromScore->NumBiallelicExactMatches,pChromScore->NumPartialMatches,pChromScore->NumNonRefAlleles, pChromScore->ExactScore, pChromScore->PartialScore);
				if(LineOfs + 1000 >= cOutBuffSize)
					{
					fwrite(pszLineBuff, 1, LineOfs, pOutScoreFile);
					fflush(pOutScoreFile);
					LineOfs = 0;
					}
				}
			}

		}
	if (LineOfs)
		{
		fwrite(pszLineBuff, 1, LineOfs, pOutScoreFile);
		fflush(pOutScoreFile);
		LineOfs = 0;
		}
	if (pChromScores != nullptr)
		{
		delete[]pChromScores;
		pChromScores = nullptr;
		}

	for (PBAsID = 1; PBAsID <= (NumRefPBAs + NumSrcPBAs); PBAsID++)
		DeleteSampleChromPBAs(PBAsID, pChromMetadata->ChromID);
	}
if (LineOfs)
	fwrite(pszLineBuff, 1, LineOfs, pOutScoreFile);
fflush(pOutScoreFile);
fclose(pOutScoreFile);

if(pFndrPBAs != nullptr)
	delete[]pFndrPBAs;
if(pszLineBuff != nullptr)
	delete[]pszLineBuff;
if (pChromScores != nullptr)
	delete[]pChromScores;

// finished with the score generation
gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenPBAsHomozygosityScores: Completed scoring all source PBAs against all reference PBAs");
return(eBSFSuccess);
}


int
CCallHaplotypes::AlignSelfPBAsKMer(int32_t NumRefPBAs, int32_t KMerSize, int32_t ChromID, int32_t ChromLen, int32_t BinsThisChrom, int32_t MaxBinSize, tsCHChromScores* pChromScores, uint8_t* pFndrPBAs[], int MaxThreads)
{
	int Rslt = -1;
	tsCHWorkerSelfScoreInstance CHWorkerSelfScoreInstances[cMaxPBAWorkerThreads];
	int NumThreads;
	tsCHWorkerSelfScoreInstance* pThreadPar;
	int32_t ThreadIdx;

#ifndef _WIN32
	// increase the default thread stack to at least cWorkThreadStackSize
	size_t defaultStackSize;
	pthread_attr_t threadattr;
	pthread_attr_init(&threadattr);
	pthread_attr_getstacksize(&threadattr, &defaultStackSize);
	if (defaultStackSize < (size_t)cWorkThreadStackSize)
	{
		if ((Rslt = pthread_attr_setstacksize(&threadattr, (size_t)cWorkThreadStackSize)) != 0)
		{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "AlignSelfPBAs: pthread_attr_setstacksize(%d) failed, default was %zd", cWorkThreadStackSize, defaultStackSize);
			return(eBSFerrInternal);
		}
	}
#endif

	NumThreads = min(NumRefPBAs, MaxThreads);
	m_CurSelfID = 0;
	m_NumWorkerInsts = 0;
	m_CompletedWorkerInsts = 0;
	m_ExpNumWorkerInsts = NumThreads;
	pThreadPar = CHWorkerSelfScoreInstances;
	for (ThreadIdx = 1; ThreadIdx <= NumThreads; ThreadIdx++, pThreadPar++)
	{
		memset(pThreadPar, 0, sizeof(tsCHWorkerSelfScoreInstance));
#ifdef _WIN32
		pThreadPar->threadHandle = nullptr;
#else
		pThreadPar->threadID = 0;
#endif
		pThreadPar->NumRefPBAs = NumRefPBAs;
		pThreadPar->NumSrcPBAs = 0;
		pThreadPar->ChromID = ChromID;
		pThreadPar->ChromLen = ChromLen;
		pThreadPar->BinsThisChrom = BinsThisChrom;
		pThreadPar->MaxBinSize = MaxBinSize;
		pThreadPar->pChromScores = pChromScores;
		pThreadPar->pFndrPBAs = pFndrPBAs;
		pThreadPar->ThreadIdx = ThreadIdx;
		pThreadPar->pThis = this;
#ifdef _WIN32
		pThreadPar->threadHandle = (HANDLE)_beginthreadex(nullptr, cWorkThreadStackSize, WorkerAlignSelfPBAsInstance, pThreadPar, 0, &pThreadPar->threadID);
#else
		pThreadPar->threadRslt = pthread_create(&pThreadPar->threadID, &threadattr, WorkerAlignSelfPBAsInstance, pThreadPar);
#endif
	}

	// allow threads time to all startup
	// check if all threads did actually startup; if they did so then m_NumWorkerInsts will have been incremented to NumInstances
	int MaxWait = cMaxWaitThreadsStartup;		// allowing at most cMaxWaitThreadsStartup secs for all threads to startup, polling every second
	int StartedInstances;
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
	} while (StartedInstances != NumThreads && MaxWait > 0);

	if (StartedInstances == NumThreads)				// all started up?
		while (WaitWorkerThreadStatus(1000) != 0);	// returns 0 if all threads started and completed processing, 1 if all threads started but some yet to complete within WaitSecs, 2 if not all threads have started within WaitSecs

	pThreadPar = CHWorkerSelfScoreInstances;
	for (ThreadIdx = 1; ThreadIdx <= NumThreads; ThreadIdx++, pThreadPar++)
	{
#ifdef _WIN32
		if (pThreadPar->threadHandle == nullptr)
			continue;
		while (WAIT_TIMEOUT == WaitForSingleObject(pThreadPar->threadHandle, 1000))
		{
		}
		CloseHandle(pThreadPar->threadHandle);
		pThreadPar->threadHandle = nullptr;
#else
		struct timespec ts;
		int JoinRlt;
		clock_gettime(CLOCK_REALTIME, &ts);
		ts.tv_sec += 1;
		while ((JoinRlt = pthread_timedjoin_np(pThreadPar->threadID, nullptr, &ts)) != 0)
		{
			ts.tv_sec += 1;
		}

#endif
		if (pThreadPar->Rslt != eBSFSuccess)
			Rslt = pThreadPar->Rslt;
	}
	return(Rslt);
}


// thread for scoring an instance of a source PBAs (termed as Self) against all founder PBAs using KMer size subsequences
int
CCallHaplotypes::AlignSelfPBAsKMerThread(tsCHWorkerSelfScoreInstance* pPar)
{
	int Rslt;
	int BinID;
	int EndBinLoci;
	int32_t SrcID;
	int32_t RefID;
	int32_t Loci;
	uint8_t SrcPBA;
	uint8_t RefPBA;
	uint8_t* pSrcPBA;
	uint8_t* pSrcLoci;
	uint8_t* pRefLoci;
	uint32_t ReqTerminate;
	int32_t SrcPBAsIdx;
	int32_t MaxNumSrcs;
	tsCHChromScores* pChromScore;

	// one more thread instance has started
	AcquireSerialise();
	m_NumWorkerInsts++;
	ReleaseSerialise();

	SrcPBAsIdx = pPar->NumSrcPBAs > 0 ? pPar->NumRefPBAs : 0;		// PBAs layout is such that references start at index 0, followed by sources at index NumRefPBAs, if All vs All
	MaxNumSrcs = pPar->NumSrcPBAs > 0 ? pPar->NumSrcPBAs : pPar->NumRefPBAs; // sources are same as references so indexes will be identical
	if (pPar->MaxBinSize > pPar->ChromLen)
		pPar->MaxBinSize = pPar->ChromLen;

	// iterating over each unprocessed (by either this or another thread) source
	for (SrcID = 1; SrcID <= MaxNumSrcs; SrcID++, SrcPBAsIdx++)
	{
		// check if requested to early terminate
		AcquireSerialise();
		ReqTerminate = m_ReqTerminate;
		ReleaseSerialise();
		if (ReqTerminate)
			break;

		// only process a source if no other thread already started to process that source
		AcquireFastSerialise();
		if (SrcID <= m_CurSelfID)	// ensuring source is unprocessed 
		{
			ReleaseFastSerialise();
			continue;				// too late, another thread is processing or completed this source!
		}
		m_CurSelfID = SrcID;		// let other threads know this SrcID is now being processed
		ReleaseFastSerialise();

		// iterate over each of the the references and score the current source relative to these references
		for (RefID = 1; RefID <= pPar->NumRefPBAs; RefID++)
		{
			// starting from first bin at loci 0 on the source
			pChromScore = &pPar->pChromScores[(pPar->NumRefPBAs * pPar->BinsThisChrom * (SrcID - 1)) + ((RefID - 1) * pPar->BinsThisChrom)];

			if (pPar->pFndrPBAs[SrcPBAsIdx] == nullptr || pPar->pFndrPBAs[RefID - 1] == nullptr) // explicitly handle these cases, some sample chromosomes will have no PBAs. 
			{
				for (BinID = 1, Loci = 0; BinID <= pPar->BinsThisChrom && Loci < pPar->ChromLen; BinID++, Loci += pPar->MaxBinSize, pChromScore++)
				{
					memset(pChromScore, 0, sizeof(tsCHChromScores));
					pChromScore->ChromID = pPar->ChromID;
					pChromScore->SrcID = SrcID;
					pChromScore->RefID = RefID;
					pChromScore->BinID = BinID;
					pChromScore->BinLoci = Loci;
					pChromScore->BinSize = (pPar->ChromLen - Loci) > pPar->MaxBinSize ? pPar->MaxBinSize : pPar->ChromLen - Loci;
					pChromScore->ExactScore = 0.0;
					pChromScore->PartialScore = 0.0;
					pChromScore->BinSize = (pPar->ChromLen - Loci) > pPar->MaxBinSize ? pPar->MaxBinSize : pPar->ChromLen - Loci;
				}
				continue;
			}

			pSrcPBA = pPar->pFndrPBAs[SrcPBAsIdx]; // already handled cases of nullptr; unlikely but possible if sample has no PBAs on current chromosome
			pSrcLoci = pSrcPBA;
			pRefLoci = pPar->pFndrPBAs[RefID - 1];
			EndBinLoci = pPar->MaxBinSize - 1;
			BinID = 1;
			for (Loci = 0; Loci < pPar->ChromLen && BinID <= pPar->BinsThisChrom; Loci++)
			{
				if (Loci == 0 || Loci > EndBinLoci)			// initial bin or starting a new bin
				{
					if (Loci > 0)		// Loci > 0 if starting a new bin
					{
						// score just processed bin before starting on new bin
						if (pChromScore->AlignLen > 0)
						{
							// fixed down weighting of 0.5 applied to partial matchings or a match where one allele matched but there was a non-ref allele also present
							pChromScore->PartialScore = (double)(pChromScore->NumExactMatches + (pChromScore->NumPartialMatches + pChromScore->NumNonRefAlleles) / 2) / pChromScore->AlignLen;
							// only scoring exact matches
							pChromScore->ExactScore = (double)pChromScore->NumExactMatches / pChromScore->AlignLen;
						}
						else
						{
							pChromScore->ExactScore = 0.0;
							pChromScore->PartialScore = 0.0;
						}

						pChromScore++;
						EndBinLoci += pPar->MaxBinSize;
						BinID++;
					}
					memset(pChromScore, 0, sizeof(tsCHChromScores));
					pChromScore->ChromID = pPar->ChromID;
					pChromScore->SrcID = SrcID;
					pChromScore->RefID = RefID;
					pChromScore->BinID = BinID;
					pChromScore->BinLoci = Loci;
					pChromScore->BinSize = (pPar->ChromLen - Loci) > pPar->MaxBinSize ? pPar->MaxBinSize : pPar->ChromLen - Loci;
				}

				SrcPBA = *pSrcLoci++;
				RefPBA = *pRefLoci++;
				if (RefPBA == 0 || SrcPBA == 0) // both source and reference must have coverage at Loci, if not then iterate to next loci
					continue;

				// enabling partial matching on the basis that fragments are being sequenced. In a diploid both haplotypes at a given loci may not be present after sequencing, especially where there is low coverage - WGS skim reads - or  GBS reads
				// when evaluating same material GBS against WGS then about 40% of the WGS bialleles are present as monoalleles in the GBS 
				if (SrcPBA == RefPBA)	// alleles exactly matching?
				{
					pChromScore->NumExactMatches++;
					if (RefPBA == 0xf0 || RefPBA == 0xcc || RefPBA == 0xc3 || RefPBA == 0x3c || RefPBA == 0x33 || RefPBA == 0x0f)
						pChromScore->NumBiallelicExactMatches++;
				}
				else
					if (SrcPBA & RefPBA) // at least one of the reference alleles present?
					{
						if (~RefPBA & SrcPBA)
							pChromScore->NumNonRefAlleles++;
						else
							pChromScore->NumPartialMatches++;
					}
				pChromScore->AlignLen++;
			}
			// last bin still requires scoring
			if (pChromScore->AlignLen > 0)
			{
				// fixed down weighting of 0.5 applied to partial matchings or a match where one allele matched but there was a non-ref allele also present
				pChromScore->PartialScore = (double)(pChromScore->NumExactMatches + (pChromScore->NumPartialMatches + pChromScore->NumNonRefAlleles) / 2) / pChromScore->AlignLen;

				// only scoring exact matches
				pChromScore->ExactScore = (double)pChromScore->NumExactMatches / pChromScore->AlignLen;
			}
			else
			{
				pChromScore->ExactScore = 0.0;
				pChromScore->PartialScore = 0.0;
			}
		}
	}
	if (SrcID <= MaxNumSrcs)
		Rslt = -1;
	else
		Rslt = 0;
	AcquireSerialise();
	m_CompletedWorkerInsts++;
	ReleaseSerialise();
	return(Rslt);
}




int
CCallHaplotypes::GenKMerGrpHammings(int32_t KMerSize,		// use this KMer sized sequences when identifying group segregation hammings
	int32_t MinKMerHammings,			// must be at least this hamming between any two group members
	int32_t KMerNoneCoverage,				// if > 0 then allow upto this number of bp within a KMer to have no coverage, maximum is 50% of KMerSize
	char* pszHapGrpFile,             // input, previously generated by 'callhaplotypes', haplotype group file (CSV format)
	int32_t NumRefPBAs,				// number of reference PBAs to be processed
	char* pszFounderInputFiles[],	// names of input founder PBA files (wildcards allowed)
	char* pszRsltsFileBaseName)		// results are written to this file base name
{
int Rslt;
gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenKMerGrpHammings: Initialising for processing of sample PBAs");
Rslt = eBSFSuccess;		// assume success!
CSimpleGlob glob(SG_GLOB_FULLSORT);
int32_t Idx;
char* pszInFile;
int32_t ReadsetID = 0;
int32_t NumFiles;
int32_t TotNumFiles = 0;

if (m_pInBuffer == nullptr)
	{
	if ((m_pInBuffer = new uint8_t[cInBuffSize]) == nullptr)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenKMerGrpHammings: Memory allocation of %u bytes for input buffering failed - %s", cInBuffSize, strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
	m_AllocInBuff = cInBuffSize;
	}
m_InNumBuffered = 0;

if (m_hOutFile != -1)
	{
	close(m_hOutFile);
	m_hOutFile = -1;
	}

// load all samples with out actually allocating memory for each individual chromosome PBAs, but the file offset at which the chromosome PBA starts will be recorded  
for (Idx = 0; Idx < NumRefPBAs; Idx++)
	{
	glob.Init();
	if (glob.Add(pszFounderInputFiles[Idx]) < SG_SUCCESS)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenKMerGrpHammings: Unable to glob '%s' sample PBA files", pszFounderInputFiles[Idx]);
		Reset();
		return(eBSFerrOpnFile);	// treat as though unable to open file
		}
	if ((NumFiles = glob.FileCount()) <= 0)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenKMerGrpHammings: Unable to locate any sample PBA file matching '%s", pszFounderInputFiles[Idx]);
		Reset();
		return(eBSFerrFileName);
		}

	if (m_LimitPrimaryPBAs > 0 && (TotNumFiles + NumFiles) > m_LimitPrimaryPBAs)
		{
		NumFiles = m_LimitPrimaryPBAs - TotNumFiles;
		if (NumFiles == 0)
			break;
		}
	TotNumFiles += NumFiles;

	if (TotNumFiles > cMaxFounderReadsets)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenKMerGrpHammings: Can accept at most %d sample PBA files for processing, after wildcard file name expansions there are %d requested", cMaxFounderReadsets, TotNumFiles);
		Reset();
		return(eBSFerrMem);
		}

	Rslt = eBSFSuccess;
	for (int32_t FileID = 0; Rslt >= eBSFSuccess && FileID < NumFiles; ++FileID)
		{
		pszInFile = glob.File(FileID);
		ReadsetID = LoadPBAFile(pszInFile, 0, true); // loading without allocation for chromosome PBAs
		if (ReadsetID <= 0)
			{
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenKMerGrpHammings: Errors loading file '%s'", pszInFile);
			Reset();
			return(ReadsetID);
			}
		m_Fndrs2Proc[ReadsetID - 1] = 0x01;	// if loaded then assumption is that this founder will be processed
		}
	if (m_LimitPrimaryPBAs > 0 && m_NumReadsetNames >= m_LimitPrimaryPBAs)
		break;
	}
m_NumFounders = TotNumFiles;

// now load the haplotype groupings which were previously generated by 'callhaplotypes' with PMode eMCSHAllelicHapsGrps
glob.Init();
if (glob.Add(pszHapGrpFile) < SG_SUCCESS)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenKMerGrpHammings: Unable to glob haplotype group clustering file(s) '%s", pszHapGrpFile);
	Reset();
	return(eBSFerrOpnFile);	// treat as though unable to open file
	}
if ((NumFiles = glob.FileCount()) <= 0)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenKMerGrpHammings: Unable to locate any haplotype group clustering file(s) matching '%s", pszHapGrpFile);
	Reset();
	return(eBSFerrOpnFile);	// treat as though unable to open file
	}
Rslt = eBSFSuccess;
for (int32_t FileID = 0; Rslt >= eBSFSuccess && FileID < NumFiles; ++FileID)
	{
	pszInFile = glob.File(FileID);
	if (Rslt = LoadHaplotypeGroupings(pszInFile) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenKMerGrpHammings: Failed processing previously generated haplotype group clusters file '%s", pszInFile);
		Reset();
		return(Rslt);
		}
	if (m_UsedHGBinSpecs == 0)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenKMerGrpHammings: Unable to parse out any previously generated haplotype group clusters file from file(s) matching '%s", pszHapGrpFile);
		Reset();
		return(eBSFerrOpnFile);	// treat as though unable to open file
		}

	if (m_UsedHGBinSpecs > 1)
		m_mtqsort.qsort(m_pHGBinSpecs, (int64_t)m_UsedHGBinSpecs, sizeof(tsHGBinSpec), SortHGBinSpecs); // sort ascending by ChromID.StartLoci.NumLoci.Distance ascending
	}

	// now have sample PBA readsets loaded - without any allocations for actual chromosome PBAs - and haplotype grouping bins 
	// iterate over each sample input file, each time processing a single chromosome
char szOutFile[_MAX_PATH];
int32_t SampleID;
int32_t DelSampleID;
int32_t ChromID;
uint8_t** ppPBAs = nullptr;

ppPBAs = new uint8_t * [cMaxFounderReadsets + 1];
memset(ppPBAs, 0, sizeof(uint8_t*) * (cMaxFounderReadsets + 1));

for (ChromID = 1; ChromID <= m_NumChromNames; ChromID++)
	{
	uint32_t ChromSize;
	char* pszChrom = LocateChrom(ChromID);

	gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenKMerGrpHammings: Beginning to load PBAs for %s chromosome for KMer processing .... ", pszChrom);
	int WorkerThreadStatus = StartWorkerLoadChromPBAThreads(m_NumThreads, 1, m_NumFounders, ChromID, true);
	if (WorkerThreadStatus < 0)
		{
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenKMerGrpHammings: Errors loading PBAs for %s chromosome for KMer processing .... ", pszChrom);
		delete[]ppPBAs;
		return(WorkerThreadStatus);
		}
	while (WorkerThreadStatus > 0)
		{
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenKMerGrpHammings: Continuing to load chromosome PBAs for %s .... ", pszChrom);
		WorkerThreadStatus = WaitWorkerThreadStatus(60);
		}
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenKMerGrpHammings: Completed loading chromosome PBAs for %s for KMer processing", pszChrom);
	TerminateLoadChromPBAsThreads(60);

	for (SampleID = 1; SampleID <= m_NumFounders; SampleID++)
		{
		tsCHChromMetadata* pChromMetadata;
		// returned pointer to chromosome metadata
		if ((pChromMetadata = LocateChromMetadataFor(SampleID, ChromID)) == nullptr)
			{
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenKMerGrpHammings: No metadata for chromosome '%s' in at least one founder '%s', skipping this chromosome", pszChrom, LocateReadset(SampleID));
			break;;
			}
		ChromSize = pChromMetadata->ChromLen;
		if ((ppPBAs[SampleID - 1] = pChromMetadata->pPBAs) == nullptr)
			{
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenKMerGrpHammings: No PBA for chromosome '%s' in at least one founder '%s', skipping this chromosome", pszChrom, LocateReadset(SampleID));
			break;
			}
		}

	if (SampleID <= m_NumFounders)
		{
		for (DelSampleID = 1; DelSampleID < SampleID; DelSampleID++)
			DeleteSampleChromPBAs(DelSampleID, ChromID);
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenKMerGrpHammings: Unloaded PBAs for chromosome %s completed, ready for next iteration of chromosome KMer processing", pszChrom);
		continue; // slough this chromosome and try next chromosome
		}

	gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenKMerGrpHammings: Loading all PBAs for chromosome %s completed, now calling KMers", pszChrom);
	if ((Rslt = GenBinKMers(ChromID,          // requiring bin group segregating KMers over this chrom
							ChromSize,		  // chromosome is this size
							KMerSize,			// use this KMer sized sequences when identifying group segregation hammings
							MinKMerHammings,	// must be at least this hammings between any two group members
							KMerNoneCoverage,		// true if all sites in KMers must have coverage
							m_NumFounders,		// number of founders to be processed 1st PBA at *pPBAs[0]
							ppPBAs)) < 0) // ppPBAs[] pts to chromosome PBAs, for each of the chromosome founder PBAs
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenKMerGrpHammings: Fatal error whilst haplotype grouping KMer reporting on chromosome %s", pszChrom);
			Reset();
			return(Rslt);
			}

	gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenKMerGrpHammings: Unloading PBAs for chromosome %s ....", pszChrom);
	for (DelSampleID = 1; DelSampleID <= SampleID; DelSampleID++)
		DeleteSampleChromPBAs(DelSampleID, ChromID);
	}

// reload the PBAs and determine the number of instances of each KMer over the complete genome
for (ChromID = 1; ChromID <= m_NumChromNames; ChromID++)
	{
	int32_t ChromSize;
	char* pszChrom = LocateChrom(ChromID);

	gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenKMerGrpHammings: Beginning to load PBAs for %s chromosome for KMer instance copy number searching .... ", pszChrom);
	int WorkerThreadStatus = StartWorkerLoadChromPBAThreads(m_NumThreads, 1, m_NumFounders, ChromID, true);
	if (WorkerThreadStatus < 0)
	{
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenKMerGrpHammings: Errors loading PBAs for %s chromosome for KMer instance copy number searching .... ", pszChrom);
		delete[]ppPBAs;
		return(WorkerThreadStatus);
	}
	while (WorkerThreadStatus > 0)
	{
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenKMerGrpHammings: Continuing to load chromosome PBAs for %s .... ", pszChrom);
		WorkerThreadStatus = WaitWorkerThreadStatus(60);
	}
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenKMerGrpHammings: Completed loading chromosome PBAs for %s for  KMer instance searching", pszChrom);
	TerminateLoadChromPBAsThreads(60);

	for (SampleID = 1; SampleID <= m_NumFounders; SampleID++)
		{
		tsCHChromMetadata* pChromMetadata;
		// returned pointer to chromosome metadata
		if ((pChromMetadata = LocateChromMetadataFor(SampleID, ChromID)) == nullptr)
			{
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenKMerGrpHammings: No metadata for chromosome '%s' in at least one founder '%s', skipping this chromosome", pszChrom, LocateReadset(SampleID));
			break;
			}
		ChromSize = pChromMetadata->ChromLen;
		if ((ppPBAs[SampleID - 1] = pChromMetadata->pPBAs) == nullptr)
			{
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenKMerGrpHammings: No PBA for chromosome '%s' in at least one founder '%s', skipping this chromosome", pszChrom, LocateReadset(SampleID));
			break;
			}
		}

	if (SampleID <= m_NumFounders)
		{
		for (DelSampleID = 1; DelSampleID < SampleID; DelSampleID++)
			DeleteSampleChromPBAs(DelSampleID, ChromID);
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenKMerGrpHammings: Unloaded PBAs for chromosome %s completed, ready for next iteration of chromosome KMer processing", pszChrom);
		continue; // slough this chromosome and try next chromosome
		}
	
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenKMerGrpHammings: Loading all PBAs for chromosome %s completed, now adding suffix array sequence entries", pszChrom);
	m_PBASfxArray.Reset();
	m_PBASfxArray.SetMaxQSortThreads(m_NumThreads);
	for (SampleID = 1; SampleID <= m_NumFounders; SampleID++)
		{
		char *pszSample = LocateReadset(SampleID);
		if ((Rslt = m_PBASfxArray.AddEntry(ChromID,SampleID,pszSample, ppPBAs[SampleID - 1], ChromSize)) != eBSFSuccess)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenKMerGrpHammings: Fatal error whilst adding entry for '%s' on chromosome %s", pszSample,pszChrom);
			Reset();
			return(Rslt);
			}
		DeleteSampleChromPBAs(SampleID, ChromID); // PBA was copied into suffix array so can now delete
		}

	gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenKMerGrpHammings: Generating suffix array index over chromosome %s ....", pszChrom);
	m_PBASfxArray.SetMaxBaseCmpLen(KMerSize+10);
	m_PBASfxArray.SuffixSort();

	gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenKMerGrpHammings: Completed generating suffix array index over chromosome %s ....", pszChrom);
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenKMerGrpHammings: Checking for any off target KMer sequences from %d potential target KMers ....", m_UsedKMerLoci);
	int32_t NumMatchedChromLoci;
	int32_t KMerMultiInstances;
	int32_t CurKMerChromID;
	char *pszKMerChrom;
	uint32_t KMerIdx;
	tsKMerLoci* pKMerLoci;
	pKMerLoci = m_pKMerLoci;
	int32_t CurGrp;
	uint32_t CurGrpMsk;
	uint32_t GrpsKMerMsk;
	uint32_t NumGrpMembers;
	uint32_t NumGrpMismatches;
	uint8_t* pSeq;
	uint32_t *pNumOffTargets;
	uint32_t NumOffTargets;
	CurKMerChromID = 0;
	for (KMerIdx = 0; KMerIdx < m_UsedKMerLoci; KMerIdx++, pKMerLoci++)
		{
		NumMatchedChromLoci = 0;
		KMerMultiInstances = 0;
		if (CurKMerChromID != pKMerLoci->ChromID)
			{
			CurKMerChromID = pKMerLoci->ChromID;
			pszKMerChrom = LocateChrom(CurKMerChromID);
			}
		CurGrp = 0;
		CurGrpMsk = 0x0001;
		GrpsKMerMsk = pKMerLoci->GrpsKMerMsk;
		pSeq = &m_pGrpKMerSeqs[pKMerLoci->GrpKMerSeqOfs];

		while (GrpsKMerMsk != 0)
			{
			CurGrp += 1;
			if (CurGrpMsk & GrpsKMerMsk)
				{
				NumGrpMembers = *(uint32_t*)pSeq;
				pSeq += sizeof(uint32_t);
				NumGrpMismatches = *(uint32_t*)pSeq;
				pSeq += sizeof(uint32_t);
				pNumOffTargets = (uint32_t*)pSeq;
				pSeq += sizeof(int32_t);
				}
			else
				break;
			
			uint32_t HitChromID;
			uint32_t HitSeqID;
			uint32_t HitLoci;
			int64_t NxtHitIdx;
			int64_t CurHitIdx = 0;
			NumOffTargets = 0;
			while((NxtHitIdx = m_PBASfxArray.IterateExacts(pSeq, pKMerLoci->KMerSize, CurHitIdx,&HitChromID, &HitSeqID,&HitLoci)) != 0)
				{
				if(HitChromID == pKMerLoci->ChromID && HitLoci == pKMerLoci->CurLoci)
					NumMatchedChromLoci++;
				else
					NumOffTargets++;
				CurHitIdx = NxtHitIdx;
				}
			pKMerLoci->TotNumOffTargets += NumOffTargets;
			*pNumOffTargets += NumOffTargets;
			pSeq += pKMerLoci->KMerSize;
			GrpsKMerMsk &= ~CurGrpMsk;
			CurGrpMsk <<= 1;
			}
		}

	m_PBASfxArray.Reset();

	gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenKMerGrpHammings: Unloading PBAs for chromosome %s ....", pszChrom);
	for (DelSampleID = 1; DelSampleID <= SampleID; DelSampleID++)
		DeleteSampleChromPBAs(DelSampleID, ChromID);
	}
if (ppPBAs != nullptr)
	delete[]ppPBAs;

if (m_hOutFile == -1)
	{
	if (m_pszOutBuffer == nullptr)
		{
		if ((m_pszOutBuffer = new uint8_t[cOutBuffSize]) == nullptr)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenKMerGrpHammings: Memory allocation of %u bytes for output buffering failed - %s", cOutBuffSize, strerror(errno));
			Reset();
			return(eBSFerrMem);
			}
		m_AllocOutBuff = cOutBuffSize;
		}
	m_OutBuffIdx = 0;
	sprintf(szOutFile, "%s.Kmers.K%dbp_H%d_U%d.bed", pszRsltsFileBaseName, KMerSize, MinKMerHammings,KMerNoneCoverage);
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenKMerGrpHammings: Reporting KMer grouping loci to file '%s'", szOutFile);
#ifdef _WIN32
	m_hOutFile = open(szOutFile, (O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC), (_S_IREAD | _S_IWRITE));
#else
	if ((m_hOutFile = open64(szOutFile, O_WRONLY | O_CREAT, S_IREAD | S_IWRITE)) != -1)
		if (ftruncate(m_hOutFile, 0) != 0)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenKMerGrpHammings: Unable to create/truncate %s - %s", szOutFile, strerror(errno));
			Reset();
			return(eBSFerrCreateFile);
			}
#endif
	if (m_hOutFile < 0)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenKMerGrpHammings: Unable to create/truncate %s - %s", szOutFile, strerror(errno));
		Reset();
		return(eBSFerrCreateFile);
		}
	}

if (m_hOutFile != -1 && m_UsedKMerLoci)
	{
	m_mtqsort.SetMaxThreads(m_NumThreads);
	m_mtqsort.qsort(m_pKMerLoci, (int64_t)m_UsedKMerLoci, sizeof(tsKMerLoci), SortKMerLoci);
	char* pszChrom;
	uint32_t CurChromID;
	uint32_t KMerIdx;
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenKMerGrpHammings: Writing to file '%s' %u KMers", szOutFile, m_UsedKMerLoci);
	m_CurRowID = 1;
	tsKMerLoci* pKMerLoci;
	pKMerLoci = m_pKMerLoci;
	CurChromID = 0;
	for (KMerIdx = 0; KMerIdx < m_UsedKMerLoci; KMerIdx++, pKMerLoci++)
		{
		if (CurChromID != pKMerLoci->ChromID)
			{
			CurChromID = pKMerLoci->ChromID;
			pszChrom = LocateChrom(CurChromID);
			}
		m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], "%s\t%d\t%u\tGrp%d:%d:%d\n", pszChrom, pKMerLoci->CurLoci, pKMerLoci->CurLoci+ pKMerLoci->KMerSize, pKMerLoci->NumHammingGrps,pKMerLoci->MinReffRelHammingDist, pKMerLoci->MaxReffRelHammingDist);
		if (m_OutBuffIdx + 1000 > m_AllocOutBuff)
			{
			if (!CUtility::RetryWrites(m_hOutFile, m_pszOutBuffer, m_OutBuffIdx))
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenKMerGrpHammings: Fatal error in RetryWrites()");
				Reset();
				return(eBSFerrFileAccess);
				}
			m_OutBuffIdx = 0;
			}
		}
	}

if (m_hOutFile != -1)
	{
	if (m_OutBuffIdx)
		{
		if (!CUtility::RetryWrites(m_hOutFile, m_pszOutBuffer, m_OutBuffIdx))
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenKMerGrpHammings: Fatal error in RetryWrites()");
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


if (m_hOutFile == -1)
{
	if (m_pszOutBuffer == nullptr)
	{
		if ((m_pszOutBuffer = new uint8_t[cOutBuffSize]) == nullptr)
		{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenKMerGrpHammings: Memory allocation of %u bytes for output buffering failed - %s", cOutBuffSize, strerror(errno));
			Reset();
			return(eBSFerrMem);
		}
		m_AllocOutBuff = cOutBuffSize;
	}
	m_OutBuffIdx = 0;
	sprintf(szOutFile, "%s.UniqueKmers.K%dbp_H%d_U%d.bed", pszRsltsFileBaseName, KMerSize, MinKMerHammings, KMerNoneCoverage);
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenKMerGrpHammings: Reporting KMer grouping loci to file '%s'", szOutFile);
#ifdef _WIN32
	m_hOutFile = open(szOutFile, (O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC), (_S_IREAD | _S_IWRITE));
#else
	if ((m_hOutFile = open64(szOutFile, O_WRONLY | O_CREAT, S_IREAD | S_IWRITE)) != -1)
		if (ftruncate(m_hOutFile, 0) != 0)
		{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenKMerGrpHammings: Unable to create/truncate %s - %s", szOutFile, strerror(errno));
			Reset();
			return(eBSFerrCreateFile);
		}
#endif
	if (m_hOutFile < 0)
	{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenKMerGrpHammings: Unable to create/truncate %s - %s", szOutFile, strerror(errno));
		Reset();
		return(eBSFerrCreateFile);
	}
}

if (m_hOutFile != -1 && m_UsedKMerLoci)
{
	m_mtqsort.SetMaxThreads(m_NumThreads);
	m_mtqsort.qsort(m_pKMerLoci, (int64_t)m_UsedKMerLoci, sizeof(tsKMerLoci), SortKMerLoci);
	char* pszChrom;
	uint32_t CurChromID;
	uint32_t KMerIdx;
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenKMerGrpHammings: Writing to file '%s' %u KMers", szOutFile, m_UsedKMerLoci);
	m_CurRowID = 1;
	tsKMerLoci* pKMerLoci;
	pKMerLoci = m_pKMerLoci;
	CurChromID = 0;
	for (KMerIdx = 0; KMerIdx < m_UsedKMerLoci; KMerIdx++, pKMerLoci++)
		{
		if (CurChromID != pKMerLoci->ChromID)
			{
			CurChromID = pKMerLoci->ChromID;
			pszChrom = LocateChrom(CurChromID);
			}
		if(pKMerLoci->TotNumOffTargets > 0)
			continue;
		m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], "%s\t%d\t%u\tGrp%d:%d:%d\n", pszChrom, pKMerLoci->CurLoci, pKMerLoci->CurLoci + pKMerLoci->KMerSize, pKMerLoci->NumHammingGrps, pKMerLoci->MinReffRelHammingDist, pKMerLoci->MaxReffRelHammingDist);
		if (m_OutBuffIdx + 1000 > m_AllocOutBuff)
		{
			if (!CUtility::RetryWrites(m_hOutFile, m_pszOutBuffer, m_OutBuffIdx))
			{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenKMerGrpHammings: Fatal error in RetryWrites()");
				Reset();
				return(eBSFerrFileAccess);
			}
			m_OutBuffIdx = 0;
		}
	}
}

if (m_hOutFile != -1)
{
	if (m_OutBuffIdx)
	{
		if (!CUtility::RetryWrites(m_hOutFile, m_pszOutBuffer, m_OutBuffIdx))
		{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenKMerGrpHammings: Fatal error in RetryWrites()");
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

if (m_hOutFile == -1)
	{
	if (m_pszOutBuffer == nullptr)
		{
		if ((m_pszOutBuffer = new uint8_t[cOutBuffSize]) == nullptr)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenKMerGrpHammings: Memory allocation of %u bytes for output buffering failed - %s", cOutBuffSize, strerror(errno));
			Reset();
			return(eBSFerrMem);
			}
		m_AllocOutBuff = cOutBuffSize;
		}
	m_OutBuffIdx = 0;
	sprintf(szOutFile, "%s.Kmers.K%dbp_H%d_U%d.csv", pszRsltsFileBaseName, KMerSize, MinKMerHammings, KMerNoneCoverage);
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenKMerGrpHammings: Reporting KMer grouping loci to file '%s'", szOutFile);
#ifdef _WIN32
	m_hOutFile = open(szOutFile, (O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC), (_S_IREAD | _S_IWRITE));
#else
	if ((m_hOutFile = open64(szOutFile, O_WRONLY | O_CREAT, S_IREAD | S_IWRITE)) != -1)
		if (ftruncate(m_hOutFile, 0) != 0)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenKMerGrpHammings: Unable to create/truncate %s - %s", szOutFile, strerror(errno));
			Reset();
			return(eBSFerrCreateFile);
			}
#endif
	if (m_hOutFile < 0)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenKMerGrpHammings: Unable to create/truncate %s - %s", szOutFile, strerror(errno));
		Reset();
		return(eBSFerrCreateFile);
		}
	m_OutBuffIdx=sprintf((char *)m_pszOutBuffer,"\"RowID\",\"Chrom\",\"StartLoci\",\"EndLoci\",\"KMerSize\",\"NumGrps\",\"MinHammingDist\",\"MaxHammingDist\",\"CurGrp\",\"NumGrpMembers\",\"GrpErrRate\",\"NumNCs\",\"NumBiallelics\",\"NumIndeterminates\",\"TotNumOffTargets\",\"GrpNumOffTargets\",\"DiplotypeKMerSeq\"\n");
	}

if (m_hOutFile != -1 && m_UsedKMerLoci)
	{
	char szSeq[500];
	char* pszChrom;
	uint32_t CurChromID;
	uint32_t KMerIdx;
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenKMerGrpHammings: Writing to file '%s' %u KMers", szOutFile, m_UsedKMerLoci);
	m_CurRowID = 1;
	tsKMerLoci* pKMerLoci;
	pKMerLoci = m_pKMerLoci;
	CurChromID = 0;
	for (KMerIdx = 0; KMerIdx < m_UsedKMerLoci; KMerIdx++, pKMerLoci++)
		{
		if (CurChromID != pKMerLoci->ChromID)
			{
			CurChromID = pKMerLoci->ChromID;
			pszChrom = LocateChrom(CurChromID);
			}
		
		int32_t CurGrp = 0;
		uint32_t CurGrpMsk = 0x0001;
		uint32_t GrpsKMerMsk = pKMerLoci->GrpsKMerMsk;
		uint8_t *pSeq = &m_pGrpKMerSeqs[pKMerLoci->GrpKMerSeqOfs];
		
		while(GrpsKMerMsk != 0)
			{
			CurGrp += 1;
			if(CurGrpMsk & GrpsKMerMsk)
				{
				int32_t SeqIdx;
				uint8_t Allele;
				int32_t NumNC = 0;
				int32_t NumDiplotype = 0;
				int32_t NumIndeterminate = 0;
				uint32_t NumGrpMembers = *(uint32_t *)pSeq;
				pSeq += sizeof(uint32_t);
				uint32_t NumGrpErrs = *(uint32_t*)pSeq;
				pSeq += sizeof(uint32_t);
				uint32_t NumOffTarget = *(uint32_t *)pSeq;
				pSeq += sizeof(uint32_t);
				uint8_t *pCurAllele = pSeq;
				char* pszBase = szSeq;
				for(SeqIdx = 0; SeqIdx < pKMerLoci->KMerSize; SeqIdx++, pCurAllele++)
					{
					Allele = *pCurAllele;
					switch(Allele) {
						case 0xc0:
							*pszBase++ = 'a';
							*pszBase++ = 'A';
							break;
						case 0xf0:
							*pszBase++ = 'a';
							*pszBase++ = 'C';
							NumDiplotype++;
							break;
						case 0xcc:
							*pszBase++ = 'a';
							*pszBase++ = 'G';
							NumDiplotype++;
							break;
						case 0xc3:
							*pszBase++ = 'a';
							*pszBase++ = 'T';
							NumDiplotype++;
							break;
						
						case 0x30:
							*pszBase++ = 'c';
							*pszBase++ = 'C';
							break;
						case 0x3c:
							*pszBase++ = 'c';
							*pszBase++ = 'G';
							NumDiplotype++;
							break;
						case 0x33:
							*pszBase++ = 'c';
							*pszBase++ = 'T';
							NumDiplotype++;
							break;


						case 0x0c:
							*pszBase++ = 'g';
							*pszBase++ = 'G';
							break;
						case 0x0f:
							*pszBase++ = 'g';
							*pszBase++ = 'T';
							NumDiplotype++;
							break;

						case 0x03:
							*pszBase++ = 't';
							*pszBase++ = 'T';
							break;

						case 0x00:
							*pszBase++ = '-';
							*pszBase++ = '-';
							NumNC++;
							break;


						default:
							*pszBase++ = 'n';
							*pszBase++ = 'N';
							NumIndeterminate++;
							break;
						}
					}
				*pszBase = '\0';
				pSeq += pKMerLoci->KMerSize;
				m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], "%d,\"%s\",%d,%u,%d,%d,%d,%d,%d,%d,%.6f,%d,%d,%d,%u,%u,\"%s\"\n", 
														m_CurRowID++,pszChrom, pKMerLoci->CurLoci, pKMerLoci->CurLoci + pKMerLoci->KMerSize - 1, pKMerLoci->KMerSize,pKMerLoci->NumHammingGrps, pKMerLoci->MinReffRelHammingDist, pKMerLoci->MaxReffRelHammingDist, CurGrp, NumGrpMembers, 
														NumGrpErrs/(double)(NumGrpMembers * pKMerLoci->KMerSize),
														NumNC, NumDiplotype, NumIndeterminate, pKMerLoci->TotNumOffTargets, NumOffTarget,szSeq);
	
				if (m_OutBuffIdx + 10000 > m_AllocOutBuff)
					{
					if (!CUtility::RetryWrites(m_hOutFile, m_pszOutBuffer, m_OutBuffIdx))
						{
						gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenKMerGrpHammings: Fatal error in RetryWrites()");
						Reset();
						return(eBSFerrFileAccess);
						}
					m_OutBuffIdx = 0;
					}
				}
			GrpsKMerMsk &= ~CurGrpMsk;
			CurGrpMsk <<= 1;
			}
		}
	}


if (m_hOutFile != -1)
{
	if (m_OutBuffIdx)
	{
		if (!CUtility::RetryWrites(m_hOutFile, m_pszOutBuffer, m_OutBuffIdx))
		{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenKMerGrpHammings: Fatal error in RetryWrites()");
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
Reset();
gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenKMerGrpHammings: Haplotype grouping KMer reporting completed");
return(0);
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
if(m_pszOutBuffer == nullptr)
	{
	if((m_pszOutBuffer = new uint8_t[cOutBuffSize]) == nullptr)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessProgenyPBAFile: Memory allocation of %u bytes for output buffering failed - %s", cOutBuffSize, strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
	m_AllocOutBuff = cOutBuffSize;
	}
m_OutBuffIdx = 0;

char *pszProgenyReadset = LocateReadset(m_CurProgenyReadsetID);

tsCHReadsetMetadata *pProgeny = &m_Readsets[m_CurProgenyReadsetID -1];
tsProgenyFndrAligns ProgenyFndrAligns;
int32_t NumFndrUniques;
int32_t NumOverlapping;
uint8_t ProgenyAlleles;
uint8_t Allele;
bool bAlleleAccept;
int32_t AlleleStackIdx;
tsAlleleStack *pCurAlleleStack;
int32_t NumAcceptAlleles;
int32_t NumPotentalFndrs;
int32_t AlleleIdx;
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
tsCHChromMetadata *pChromMetadata;
int32_t CurChromMetadataIdx = pProgeny->StartChromMetadataIdx;
for(int32_t ChromIdx = 0; ChromIdx < pProgeny->NumChroms && CurChromMetadataIdx != 0; ChromIdx++)
	{
	pChromMetadata = &m_pChromMetadata[CurChromMetadataIdx - 1];
	if(pChromMetadata->pPBAs != nullptr)
#ifdef _WIN32
		free(pChromMetadata->pPBAs);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(pChromMetadata->pPBAs != MAP_FAILED)
		munmap(pChromMetadata->pPBAs, pChromMetadata->ChromLen);
#endif
	pChromMetadata->pPBAs = nullptr;
	}
gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessProgenyPBAFile: Located %d progeny loci overlapping %d founder AlleleStacks loci", NumOverlapping, m_UsedAlleleStacks);
return(eBSFSuccess);
}


bool					// true if chrom is accepted, false if chrom not accepted
CCallHaplotypes::AcceptThisChromID(uint32_t ChromID)
{
char *pzChrom;
if((pzChrom=LocateChrom(ChromID)) == nullptr)
    return(false);
return(AcceptThisChromName(pzChrom));
}

bool					// true if chrom is accepted, false if chrom not accepted
CCallHaplotypes::AcceptThisChromName(char* pszChrom,   // chromosome name
									 bool bKnown)	// if true then chromosome must have been previously processed and accepted by LoadChromSizes() processing
{
bool bMatch;
if (bKnown)
	return(LocateChrom(pszChrom) < 1 ? false : true);
AcquireSerialise();
bMatch = !m_RegExprs.MatchExcludeRegExpr(pszChrom);
if(bMatch)
    bMatch = m_RegExprs.MatchIncludeRegExpr(pszChrom);
ReleaseSerialise();
return(bMatch);
}




uint32_t				// returns number of unprocessed bytes in buffer
CCallHaplotypes::FillInBuffer(uint32_t MinRequired,uint32_t MaxRequired) // try and fill input buffer with at least MinRequired or refill (if < MinRequired) up to MaxRequired (if 0 then max buffer allocation)
{
int NumRead;
uint32_t UnProcessed;
if(MaxRequired == 0)
    MaxRequired = m_AllocInBuff;
else
    if(MaxRequired < MinRequired)
		MaxRequired = MinRequired;
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
else
	{
	m_InNumBuffered = 0;
	m_InNumProcessed = 0;
	}

// attempt to fill input buffer to it's MaxRequired
do{
	if((NumRead = read(m_hInFile, &m_pInBuffer[m_InNumProcessed], MaxRequired - m_InNumBuffered)) < 0)		
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error %s attempting to read from file", strerror(errno));
		return(0);
		}
	m_InNumBuffered += (uint32_t)NumRead;
	}
while(NumRead > 0 && m_InNumBuffered < MinRequired);
m_InFileOfs = _lseeki64(m_hInFile,0,SEEK_CUR);             // record file offset at which next file read will start from
return(m_InNumBuffered);
}

int32_t					// returned readset identifier (1..n) or < 0 if errors
CCallHaplotypes::LoadPBAFile(char* pszFile,	// load chromosome metadata and PBA data from this file, filters out chroms - AcceptThisChromName()
							uint8_t ReadsetType,	// 0: founder, 1: progeny, 2: control
							 bool bChromMetaOnly)  // load chrom metadata (chrom name,length, file offset at which chrom PBAs start) but don't actually load the chromosome PBAs
{
int Rslt;
int Version;
char szExperimentID[100];
char szRefAssemblyID[100];
char szReadsetID[100];
uint32_t PrevChromMetadataIdx;
uint32_t ChromID;
uint32_t ReadsetID;
tsCHReadsetMetadata *pReadsetMetadata;
tsCHChromMetadata *pChromMetadata;
tsCHChromMetadata *pPrevChromMetadata;
int64_t FileOfsPBA;
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

// attempt to load the readset metadata 
m_InNumBuffered = 0;
m_InNumProcessed = 0;
if((FillInBuffer(500,min(m_AllocInBuff,500u)) == 0) || m_InNumBuffered < 9) // 500 will cover the maximally sized header
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to read at least a partial header from input file '%s'", pszFile);
	return(eBSFerrOpnFile);
	}


// check file type is PBA
if(strncmp((char *)m_pInBuffer,"Type:PbA",8))
	{
	gDiagnostics.DiagOut (eDLFatal, gszProcName, "Input file '%s' exists, unable to parse file type header tag as being a packed base allele file",pszFile);
	return(eBSFerrOpnFile);
	}
// file type in header expected to be followed by a new line, not carriage return then a new line!
if(m_pInBuffer[8] != '\n' && m_pInBuffer[9] != 'V')
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Input file '%s' exists, file type header tag is for a PBA type but this tag is incorrectly terminated. Has file been transformed into ascii?", pszFile);
	return(eBSFerrOpnFile);
	}
m_InNumProcessed = 9;

if(m_InNumBuffered < 500)
	{
	gDiagnostics.DiagOut (eDLFatal, gszProcName, "Input file '%s' exists, unable to read complete header and initial chromosome metadata, expecting file size >= 500, is it a packed base allele file",pszFile);
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

m_InNumProcessed += scanlen + 1;		// header tags were '\n' terminated except for final which was '\0' terminated

m_InFileOfs = m_InNumProcessed;         // 1st chromosome starts immediately following the header metadata
if(_lseeki64(m_hInFile, m_InFileOfs, SEEK_SET) != m_InFileOfs)
	{
	gDiagnostics.DiagOut (eDLFatal, gszProcName, "Input file '%s' exists as a packed base allele file but unable to seek past header metadata",pszFile);
	return(eBSFerrOpnFile);
	}

if((ReadsetID = AddReadset(szReadsetID, ReadsetType))==0)
	{
	gDiagnostics.DiagOut (eDLFatal, gszProcName, "Input file '%s' duplicates the ReadsetID '%s' of a previously loaded readset",pszFile,szReadsetID);
	return(eBSFerrOpnFile);
	}

pReadsetMetadata = &m_Readsets[ReadsetID-1];
memset(pReadsetMetadata,0,sizeof(*pReadsetMetadata));
pReadsetMetadata->ReadsetType = ReadsetType;
pReadsetMetadata->NumChroms = 0;
strcpy(pReadsetMetadata->szFileName, pszFile);
strcpy(pReadsetMetadata->szExperimentID,szExperimentID);
strcpy(pReadsetMetadata->szRefAssemblyID,szRefAssemblyID);
pReadsetMetadata->ReadsetID = ReadsetID;
pReadsetMetadata->StartChromID = 0;
pReadsetMetadata->StartChromMetadataIdx = 0;
pReadsetMetadata->NxtFileChromOfs = m_InFileOfs;

int ChromNameLen;
char *pszChromName;
int32_t ChromLen;

// iterate over all chromosomes
m_InNumProcessed = 0;
m_InNumBuffered = 0;
PrevChromMetadataIdx = 0;
while(FillInBuffer((uint32_t)110,110)==110) // reading chromosome metadata
	{
	pBuff = &m_pInBuffer[m_InNumProcessed++];
	ChromNameLen = (int)*pBuff++;
	pszChromName = (char *)pBuff;
	pBuff += (int64_t)1+ChromNameLen;
	ChromLen = *(int32_t *)pBuff;
	FileOfsPBA = pReadsetMetadata->NxtFileChromOfs + ChromNameLen + 6;
	pReadsetMetadata->NxtFileChromOfs += (int64_t)ChromNameLen + 6 + ChromLen;
	// check if this chromosome is to be retained for further processing
	if(!AcceptThisChromName(pszChromName))
		{
		// not accepting this chromosome
		m_InFileOfs = pReadsetMetadata->NxtFileChromOfs;         // skip over curent chrom's PBAs to start of next chromosome
		if((pReadsetMetadata->NxtFileChromOfs = _lseeki64(m_hInFile, m_InFileOfs, SEEK_SET)) != m_InFileOfs)
			break;
		m_InNumProcessed = 0;
		m_InNumBuffered = 0;
		continue;
		}

		// before accepting chrom then ensure that it's PBA length matches the BED chromosome sizes
	ChromID = AddChrom(pszChromName);
	if (ChromLen != m_ChromSizes[ChromID - 1])
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Input file '%s' with readset '%s' has chromosome '%s' size mismatch - expected size was %d, actual size %d", pszFile, szReadsetID, pszChromName, m_ChromSizes[ChromID - 1], ChromLen);
		return(eBSFerrChrom);
		}

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
	pChromMetadata->HGBinID = 0;
	pChromMetadata->FileOfsPBA = FileOfsPBA;
	pChromMetadata->pPBAs = nullptr;
	if(bChromMetaOnly) // loading chrom metadata only without actually loading the chromosome PBAs?
		{
		m_InFileOfs = pReadsetMetadata->NxtFileChromOfs;         // skip over curent chrom's PBAs to start of next chromosome
		if((pReadsetMetadata->NxtFileChromOfs = _lseeki64(m_hInFile, m_InFileOfs, SEEK_SET)) != m_InFileOfs)
			break;
		m_InNumProcessed = 0;
		m_InNumBuffered = 0;
		continue;
		}
	// loading PBAs
	pChromMetadata->pPBAs = AllocPBAs(ChromLen);
	m_InNumProcessed += ChromNameLen + 5;
	FillInBuffer((uint32_t)ChromLen,ChromLen);
	pBuff = &m_pInBuffer[m_InNumProcessed];


	memcpy(pChromMetadata->pPBAs, pBuff, ChromLen);
	m_InNumProcessed += ChromLen;
	//	validate PBA allele composition, earlier releases were retaining very low allele proportions so validating and in-place replacing these as being non-alignments
	int NumErrs = ValidatePBAs(ChromLen,pChromMetadata->pPBAs,true,true);
	//
	if(ReadsetType == 0)
		TrimPBAs(m_FndrTrim5, m_FndrTrim3, ChromLen,pChromMetadata->pPBAs);
	else // treating controls as if progeny when trimming
		TrimPBAs(m_ProgTrim5, m_ProgTrim3, ChromLen, pChromMetadata->pPBAs);
	}

close(m_hInFile);
m_hInFile = -1;
return((int32_t)ReadsetID);
}

bool  // false: no previous memory allocated for containing the chromosome PBAs, true: allocation was present and has been deleted
CCallHaplotypes::DeleteSampleChromPBAs(int32_t SampleID,   // Sample identifier
									 int32_t ChromID)    // chrom identifier
{
tsCHChromMetadata* pChromMetadata;

// returned pointer to chromosome metadata
if((pChromMetadata = LocateChromMetadataFor(SampleID, ChromID)) == nullptr)
	return(false);
if(pChromMetadata->pPBAs == nullptr) 
	return(false);
#ifdef _WIN32
free(pChromMetadata->pPBAs);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
if(pChromMetadata->pPBAs != MAP_FAILED)
	munmap(pChromMetadata->pPBAs, pChromMetadata->ChromLen);
#endif
pChromMetadata->pPBAs = nullptr;
return(true);
}

uint8_t *  // returned PBA or nullptr if unable to load PBAs for requested sample.chrom 
CCallHaplotypes::LoadSampleChromPBAs(int32_t SampleID,   // Sample identifier
				   int32_t ChromID)    // chrom identifier specifying which PBAs is to be loaded from SampleID file
{
int hInFile;
char* pszChrom;
int64_t ChromSeekOfs;
tsCHChromMetadata* pChromMetadata;
tsCHReadsetMetadata *pReadsetMetadata;
int32_t NumRead;
int32_t NumLoaded;

// returned pointer to chromosome metadata
if((pChromMetadata = LocateChromMetadataFor(SampleID, ChromID)) == nullptr)
	return(nullptr);

pszChrom = LocateChrom(ChromID);

if(pChromMetadata->pPBAs != nullptr) // PBAs may already have been loaded for this chromosome, delete as may be stale
	{
#ifdef _WIN32
	free(pChromMetadata->pPBAs);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(pChromMetadata->pPBAs != MAP_FAILED)
		munmap(pChromMetadata->pPBAs, pChromMetadata->ChromLen);
#endif
	pChromMetadata->pPBAs = nullptr;
	}

if(m_PMode == eMCSHCoverageHapsGrps)
	return(LoadPBAChromCoverage(SampleID, ChromID));

// need to actually load from file
pReadsetMetadata = &m_Readsets[SampleID-1];
#ifdef _WIN32
hInFile = open(pReadsetMetadata->szFileName, O_READSEQ);		// file access is normally sequential..
#else
hInFile = open64(pReadsetMetadata->szFileName, O_READSEQ);		// file access is normally sequential..
#endif
if(hInFile == -1)							// check if file open succeeded
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadSampleChromPBA: Unable to open input file '%s' : %s", pReadsetMetadata->szFileName, strerror(errno));
	return(nullptr);
	}
if((ChromSeekOfs = _lseeki64(hInFile, pChromMetadata->FileOfsPBA, SEEK_SET)) != pChromMetadata->FileOfsPBA)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Input file '%s' opened, _lseek() to file offset %zd failed for chrom '%s', returned %zd, error: '%s'",
										pReadsetMetadata->szFileName, pChromMetadata->FileOfsPBA, ChromSeekOfs, pszChrom, strerror(errno));
	close(hInFile);
	hInFile = -1;
	return(nullptr);
	}
pChromMetadata->pPBAs = AllocPBAs(pChromMetadata->ChromLen);
// attempt to load Chromosomes PBAs from file
NumLoaded = 0;
do{
	if((NumRead = read(hInFile, &pChromMetadata->pPBAs[NumLoaded], pChromMetadata->ChromLen - NumLoaded)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error %s attempting to read chromosome '%s' from file '%s'", strerror(errno),pszChrom,pReadsetMetadata->szFileName);
		close(hInFile);
		hInFile = -1;
		if(pChromMetadata->pPBAs == nullptr) 
			return(nullptr);
#ifdef _WIN32
		free(pChromMetadata->pPBAs);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
		if(pChromMetadata->pPBAs != MAP_FAILED)
			munmap(pChromMetadata->pPBAs, pChromMetadata->ChromLen);
#endif
		pChromMetadata->pPBAs = nullptr;
		return(nullptr);
		}
	NumLoaded += NumRead;
	}
while(NumRead > 0 && NumLoaded < pChromMetadata->ChromLen);
close(hInFile);
hInFile = -1;
//	validate PBA allele composition, earlier releases of 'kalign' were retaining very low allele proportions so validating and in-place replacing these as being non-alignments
int NumErrs = ValidatePBAs(pChromMetadata->ChromLen,pChromMetadata->pPBAs,true,true);
//
if(pReadsetMetadata->ReadsetType == 0)
	TrimPBAs(m_FndrTrim5, m_FndrTrim3, pChromMetadata->ChromLen,pChromMetadata->pPBAs);
else // treating controls as if progeny when trimming
	TrimPBAs(m_ProgTrim5, m_ProgTrim3, pChromMetadata->ChromLen, pChromMetadata->pPBAs);
return(pChromMetadata->pPBAs);
}

// trim 5' and 3' aligned segments within the PBAs attempting to reduce sequencing error induced false alleles
int			// returns number of non-trimmed loci in the pPBAs
CCallHaplotypes::TrimPBAs(uint32_t Trim5,	// trim 5' this many aligned PBA bases from each aligned segment
		 uint32_t Trim3,	// trim 3' this many aligned PBA bases from each aligned segment
		 uint32_t PBALen,	// pPBAs contains this many packed base alleles
		uint8_t* pPBAs)		// packed base alleles to be processed
{
uint32_t NonTrimmed;
uint32_t ToTrim;
uint32_t Loci;
uint8_t* pPBA;
if(pPBAs == nullptr || PBALen == 0)
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
tsCHChromMetadata *pChromMetadata;
size_t memreq;
if (m_pChromMetadata == nullptr)					// may be nullptr first time in
	{
	memreq = cAllocChromMetadata * sizeof(tsCHChromMetadata);
#ifdef _WIN32
	m_pChromMetadata = (tsCHChromMetadata *)malloc((size_t)memreq);
	if (m_pChromMetadata == nullptr)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AllocChromMetadata: Memory allocation of %zd bytes failed",(int64_t)memreq);
		return(eBSFerrMem);
		}
#else
	m_pChromMetadata = (tsCHChromMetadata *)mmap(nullptr, (size_t)memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	if (m_pChromMetadata == MAP_FAILED)
		{
		m_pChromMetadata = nullptr;
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AllocChromMetadata: Memory allocation of %zd bytes through mmap()  failed", (int64_t)memreq, strerror(errno));
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
		size_t memreq = ToAllocdChromMetadata * sizeof(tsCHChromMetadata);
#ifdef _WIN32
		pChromMetadata = (tsCHChromMetadata*)realloc(m_pChromMetadata, memreq);
		if (pChromMetadata == nullptr)
			{
#else
			pChromMetadata = (tsCHChromMetadata*)mremap(m_pChromMetadata, m_AllocdChromMetadataMem, memreq, MREMAP_MAYMOVE);
			if (pChromMetadata == MAP_FAILED)
			{
#endif
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "AllocChromMetadata: Memory reallocation to %zd bytes failed - %s", (int64_t)memreq, strerror(errno));
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
if (pPBAs == nullptr)
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "AllocPBAs Memory allocation of %zd bytes failed",(int64_t)memreq);
#else
pPBAs = (uint8_t*)mmap(nullptr, (size_t)memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
if (pPBAs == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "AllocPBAs: Memory allocation of %zd bytes through mmap()  failed",(int64_t)memreq, strerror(errno));
	pPBAs = nullptr;
	}
#endif
return(pPBAs);
}


uint8_t *								// returned pointer to start of PBA
CCallHaplotypes::LocatePBAfor(int32_t ReadSetID,		// readset identifier 
			 int32_t ChromID)			// chrom identifier
{
tsCHReadsetMetadata *pReadsetMetadata;
tsCHChromMetadata *pChromMetadata;
int32_t CurChromMetadataIdx;

if(ReadSetID > m_NumReadsetNames || ReadSetID == 0)
	return(nullptr);
pReadsetMetadata = &m_Readsets[ReadSetID-1];
CurChromMetadataIdx = pReadsetMetadata->StartChromMetadataIdx;
for(int32_t ChromIdx = 0; ChromIdx < pReadsetMetadata->NumChroms && CurChromMetadataIdx != 0; ChromIdx++)
	{
	pChromMetadata = &m_pChromMetadata[CurChromMetadataIdx - 1];
	if(pChromMetadata->ChromID == ChromID)
		return(pChromMetadata->pPBAs);
	CurChromMetadataIdx = pChromMetadata->NxtChromMetadataIdx;
	}
return(nullptr);
}

tsCHChromMetadata *								// returned pointer to chromosome metadata
CCallHaplotypes::LocateChromMetadataFor(int32_t ReadSetID,		// readset identifier 
			 int32_t ChromID)			// chrom identifier
{
tsCHReadsetMetadata *pReadsetMetadata;
tsCHChromMetadata *pChromMetadata;
int32_t CurChromMetadataIdx;

if(ReadSetID > m_NumReadsetNames || ReadSetID == 0)
	return(nullptr);
pReadsetMetadata = &m_Readsets[ReadSetID-1];
CurChromMetadataIdx = pReadsetMetadata->StartChromMetadataIdx;
for(int32_t ChromIdx = 0; ChromIdx < pReadsetMetadata->NumChroms && CurChromMetadataIdx != 0; ChromIdx++)
	{
	pChromMetadata = &m_pChromMetadata[CurChromMetadataIdx - 1];
	if(pChromMetadata->ChromID == ChromID)
		return(pChromMetadata);
	CurChromMetadataIdx = pChromMetadata->NxtChromMetadataIdx;
	}
return(nullptr);
}




size_t									// returned index+1 into  m_pProgenyFndrAligns[] to allocated and initialised ProgenyFndrAligns, 0 if errors
CCallHaplotypes::AddProgenyFndrAligns(tsProgenyFndrAligns *pInitProgenyFndrAligns)	// allocated tsProgenyFndrAligns to be initialised with a copy of pInitProgenyFndrAligns
{
size_t ToAllocdProgenyFndrAligns;
tsProgenyFndrAligns *pProgenyFndrAligns;
size_t memreq;
AcquireSerialise();
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
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddWinBinCnts: Memory allocation of %zd bytes through mmap()  failed",(int64_t)memreq, strerror(errno));
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
		if (pProgenyFndrAligns == nullptr)
			{
#else
		pProgenyFndrAligns = (tsProgenyFndrAligns*)mremap(m_pProgenyFndrAligns, m_AllocdProgenyFndrAlignsMem, memreq, MREMAP_MAYMOVE);
		if (pProgenyFndrAligns == MAP_FAILED)
			{
#endif
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddProgenyFndrAligns: Memory reallocation to %zd bytes failed - %s", (int64_t)memreq, strerror(errno));
			return(0);
			}
		m_pProgenyFndrAligns = pProgenyFndrAligns;
		m_AllocdProgenyFndrAlignsMem = memreq;
		m_AllocdProgenyFndrAligns = ToAllocdProgenyFndrAligns;
		}
pProgenyFndrAligns = &m_pProgenyFndrAligns[m_UsedProgenyFndrAligns++];
if(pInitProgenyFndrAligns != nullptr)
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


int32_t									// returned index+1 into m_pDGTLoci[] to allocated and initialised DGTLoci, 0 if errors
CCallHaplotypes::AddDGTLoci(int32_t SrcExprID,  // sourced CSV haplotype group experiment identifier
					int32_t SrcRowID,   // sourced CSV haplotype group row identifier
					int32_t ChromID,    // DGT loci is on this chromosome
					int32_t DGTLoci,    // at this loci
					int32_t NumHapGrps, // DGT has this number of actual haplotype groups
					int32_t AlleleGrpIDs[4],   // group identifier for each allele accepted as having the highest F-measure, 0 if no group having accepted F-measure
					double AlleleFMeasures[4],  // F-measure for each allele - 0.0 if no group having an accepted F-measure for that allele
					tsGrpCnts GrpCnts[6])    // counts for 1st 5 groups plus a pseudo-group containing sum of counts for groups after the 1st 5

{
double Xchg;
tsDGTLoci Tmp;
Tmp.SrcExprID = SrcExprID;
Tmp.SrcRowID = SrcRowID;
Tmp.ChromID = ChromID;
Tmp.DGTLoci = DGTLoci;
Tmp.NumHapGrps = NumHapGrps;
memcpy(Tmp.GrpCnts, GrpCnts, 6 * sizeof(tsGrpCnts));
memcpy(Tmp.AlleleGrpIDs,AlleleGrpIDs,4*sizeof(uint32_t));
memcpy(Tmp.AlleleFMeasures,AlleleFMeasures,4*sizeof(double));
// sorting Fmeasures descending 
memcpy(Tmp.AlleleFMeasuresDesc,AlleleFMeasures,4*sizeof(double));
if(Tmp.AlleleFMeasuresDesc[0] < Tmp.AlleleFMeasuresDesc[1])
	{
	Xchg = Tmp.AlleleFMeasuresDesc[0];
	Tmp.AlleleFMeasuresDesc[0] = Tmp.AlleleFMeasuresDesc[1];
	Tmp.AlleleFMeasuresDesc[1] = Xchg;
	}
if(Tmp.AlleleFMeasuresDesc[1] < Tmp.AlleleFMeasuresDesc[2])
	{
	Xchg = Tmp.AlleleFMeasuresDesc[1];
	Tmp.AlleleFMeasuresDesc[1] = Tmp.AlleleFMeasuresDesc[2];
	Tmp.AlleleFMeasuresDesc[2] = Xchg;
	}
if(Tmp.AlleleFMeasuresDesc[2] < Tmp.AlleleFMeasuresDesc[3])
	{
	Xchg = Tmp.AlleleFMeasuresDesc[2];
	Tmp.AlleleFMeasuresDesc[2] = Tmp.AlleleFMeasuresDesc[3];
	Tmp.AlleleFMeasuresDesc[3] = Xchg;
	}
if(Tmp.AlleleFMeasuresDesc[0] < Tmp.AlleleFMeasuresDesc[1])
	{
	Xchg = Tmp.AlleleFMeasuresDesc[0];
	Tmp.AlleleFMeasuresDesc[0] = Tmp.AlleleFMeasuresDesc[1];
	Tmp.AlleleFMeasuresDesc[1] = Xchg;
	}
if(Tmp.AlleleFMeasuresDesc[1] < Tmp.AlleleFMeasuresDesc[2])
	{
	Xchg = Tmp.AlleleFMeasuresDesc[1];
	Tmp.AlleleFMeasuresDesc[1] = Tmp.AlleleFMeasuresDesc[2];
	Tmp.AlleleFMeasuresDesc[2] = Xchg;
	}
if(Tmp.AlleleFMeasuresDesc[0] < Tmp.AlleleFMeasuresDesc[1])
	{
	Xchg = Tmp.AlleleFMeasuresDesc[0];
	Tmp.AlleleFMeasuresDesc[0] = Tmp.AlleleFMeasuresDesc[1];
	Tmp.AlleleFMeasuresDesc[1] = Xchg;
	}
return(AddDGTLoci(&Tmp));
}

int32_t									// returned index+1 into m_pHGBinSpecs[] to allocated and initialised allele stack, 0 if errors
CCallHaplotypes::AddDGTLoci(tsDGTLoci* pInitDGTLoci)	// allocated tsDGTLoci to be initialised with a copy of pInitDGTLoci
{
uint32_t ToAllocdDGTLoci;
tsDGTLoci* pDGTLoci;
size_t memreq;
AcquireFastSerialise();
if(m_pDGTLoci == nullptr)					// may be nullptr first time in
	{
	memreq = (size_t)cAllocHGBinSpecs * sizeof(tsDGTLoci);
#ifdef _WIN32
	m_pDGTLoci = (tsDGTLoci*)malloc((size_t)memreq);
	if(m_pDGTLoci == nullptr)
		{
		ReleaseFastSerialise();
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddDGTLoci: Memory allocation of %zd bytes failed", (int64_t)memreq);
		return(0);
		}
#else
	m_pDGTLoci = (tsDGTLoci*)mmap(nullptr, (size_t)memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	if(m_pDGTLoci == MAP_FAILED)
		{
		m_pDGTLoci = nullptr;
		ReleaseFastSerialise();
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddDGTLoci: Memory allocation of %zd bytes through mmap()  failed", (int64_t)memreq, strerror(errno));
		return(0);
		}
#endif
	m_AllocdDGTLociMem = memreq;
	m_AllocdDGTLoci = cAllocHGBinSpecs;
	m_UsedDGTLoci = 0;
	}
else
		// needing to allocate more memory?
	if((m_UsedDGTLoci) >= m_AllocdDGTLoci)
		{
		ToAllocdDGTLoci = m_UsedDGTLoci + cReallocHGBinSpecs;
		size_t memreq = (size_t)ToAllocdDGTLoci * sizeof(tsDGTLoci);
#ifdef _WIN32
		pDGTLoci = (tsDGTLoci*)realloc(m_pDGTLoci, memreq);
		if(pDGTLoci == nullptr)
			{
#else
		pDGTLoci = (tsDGTLoci*)mremap(m_pDGTLoci, m_AllocdDGTLociMem, memreq, MREMAP_MAYMOVE);
		if(pDGTLoci == MAP_FAILED)
			{
#endif
			ReleaseFastSerialise();
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddDGTLoci: Memory reallocation to %zd bytes failed - %s", (int64_t)memreq, strerror(errno));
			return(0);
			}
		m_pDGTLoci = pDGTLoci;
		m_AllocdDGTLociMem = memreq;
		m_AllocdDGTLoci = ToAllocdDGTLoci;
		}
pDGTLoci = &m_pDGTLoci[m_UsedDGTLoci++];
if(pInitDGTLoci != nullptr)
	*pDGTLoci = *pInitDGTLoci;
else
	memset(pDGTLoci, 0, sizeof(tsDGTLoci));
ReleaseFastSerialise();
return(m_UsedDGTLoci);
}

int64_t									// returned index+1 into m_pGrpKMerSeqs to allocated and copied KMerSeq of size KMerSize, returns 0 if errors
CCallHaplotypes::AddKMerSeq(int32_t NumGrpMembers, // group for which this KMer sequence is representative has this many members
							int32_t NumMismatches,	// over all group members, there were this many mismatches with the consensus KMer sequence
							int32_t KMerSize,	//  consensus KMer sequence to copy is this size
							uint8_t *pKMerSeq)	// copy consensus from here
{
size_t KMerCpyOfs;
uint8_t* pSeq;
size_t memreq;
AcquireFastSerialise();
if (m_pGrpKMerSeqs == nullptr)					// may be nullptr first time in
	{
	memreq = cInitialAllocKMerSeqs;
#ifdef _WIN32
	m_pGrpKMerSeqs = (uint8_t *)malloc((size_t)memreq);
	if (m_pGrpKMerSeqs == nullptr)
		{
		ReleaseFastSerialise();
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddKMerSeq: Memory allocation of %zd bytes failed", (int64_t)memreq);
		return(0);
		}
#else
	m_pGrpKMerSeqs = (uint8_t *)mmap(nullptr, (size_t)memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	if (m_pGrpKMerSeqs == MAP_FAILED)
		{
		m_pGrpKMerSeqs = nullptr;
		ReleaseFastSerialise();
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddKMerSeq: Memory allocation of %zd bytes through mmap()  failed", (int64_t)memreq, strerror(errno));
		return(0);
		}
#endif
	m_AllocdKMerSeqMem = memreq;
	m_UsedKMerSeqMem = 0;
	}
else
		// needing to allocate more memory?
	if ((m_UsedKMerSeqMem + KMerSize + (3*sizeof(int32_t))) >= m_AllocdKMerSeqMem)
		{
		memreq = m_AllocdKMerSeqMem + cReallocKMerSeqs;
#ifdef _WIN32
		pSeq = (uint8_t *)realloc(m_pGrpKMerSeqs, memreq);
		if (pSeq == nullptr)
			{
#else
			pSeq = (uint8_t *)mremap(m_pGrpKMerSeqs, m_AllocdKMerSeqMem, memreq, MREMAP_MAYMOVE);
			if (pSeq == MAP_FAILED)
			{
#endif
				ReleaseFastSerialise();
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddKMerSeq: Memory reallocation to %zd bytes failed - %s", (int64_t)memreq, strerror(errno));
				return(0);
			}
		m_pGrpKMerSeqs = pSeq;
		m_AllocdKMerSeqMem = memreq;
		}
KMerCpyOfs = m_UsedKMerSeqMem;
m_UsedKMerSeqMem += KMerSize + (3*sizeof(int32_t));
pSeq = &m_pGrpKMerSeqs[KMerCpyOfs];
*(int32_t *)pSeq = NumGrpMembers;
pSeq+= sizeof(int32_t);
*(int32_t*)pSeq = NumMismatches;
pSeq += sizeof(int32_t);
*(int32_t*)pSeq = 0;
pSeq += sizeof(int32_t);
memcpy(pSeq, pKMerSeq, KMerSize);
ReleaseFastSerialise();
return(KMerCpyOfs+1);
}


int
CCallHaplotypes::AddKMerLoci(int32_t ChromID,	// KMer is on this chromosome
	int32_t CurLoci,			// starting at this loci
	int32_t KMerSize,			// and for this number of base pairs
	int32_t NumHammingGrps,		// number of actual groups for which hamming differentials were generated
	uint32_t GrpsKMerMsk,			// bit mask of those groups for which are represented in NumHammingGrps (group1 -> bit0, group2 -> bit1, groupN -> bitN-1)
	int64_t GrpKMerSeqOfs,		// offset in m_pGrpKMerSeqs at which 1st group in GrpsKMer sequence starts, 2nd group at GrpKMerSeqOfs+KMerSize ...
	int32_t MinReffRelHammingDist, // the actual minimal hamming between any two groups
	int32_t MaxReffRelHammingDist) // the actual maximal hamming between any two groups
{
tsKMerLoci Tmp;
Tmp.ChromID = ChromID;
Tmp.CurLoci = CurLoci;
Tmp.KMerSize = KMerSize;
Tmp.NumHammingGrps = NumHammingGrps;
Tmp.GrpsKMerMsk = GrpsKMerMsk;
Tmp.GrpKMerSeqOfs = GrpKMerSeqOfs;
Tmp.MinReffRelHammingDist = MinReffRelHammingDist;
Tmp.MaxReffRelHammingDist = MaxReffRelHammingDist;
Tmp.TotNumOffTargets = 0;
return(AddKMerLoci(&Tmp));
}

int32_t									// returned index+1 into m_pKMerLoci[] to allocated and initialised KMers, 0 if errors
CCallHaplotypes::AddKMerLoci(tsKMerLoci* pInitKMerLoci)	// allocated tsKMerLoci to be initialised with a copy of pInitKMerLoci
{
uint32_t ToAllocdKMerLoci;
tsKMerLoci* pKMerLoci;
size_t memreq;
AcquireFastSerialise();
if (m_pKMerLoci == nullptr)					// may be nullptr first time in
	{
	memreq = (size_t)cAllocHGBinSpecs * sizeof(tsKMerLoci);
#ifdef _WIN32
	m_pKMerLoci = (tsKMerLoci*)malloc((size_t)memreq);
	if (m_pKMerLoci == nullptr)
		{
		ReleaseFastSerialise();
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddKMerLoci: Memory allocation of %zd bytes failed", (int64_t)memreq);
		return(0);
		}
#else
	m_pKMerLoci = (tsKMerLoci*)mmap(nullptr, (size_t)memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	if (m_pKMerLoci == MAP_FAILED)
		{
		m_pKMerLoci = nullptr;
		ReleaseFastSerialise();
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddKMerLoci: Memory allocation of %zd bytes through mmap()  failed", (int64_t)memreq, strerror(errno));
		return(0);
		}
#endif
	m_AllocdKMerLociMem = memreq;
	m_AllocdKMerLoci = cAllocHGBinSpecs;
	m_UsedKMerLoci = 0;
	}
else
		// needing to allocate more memory?
	if ((m_UsedKMerLoci) >= m_AllocdKMerLoci)
		{
		ToAllocdKMerLoci = m_UsedKMerLoci + cReallocHGBinSpecs;
		size_t memreq = (size_t)ToAllocdKMerLoci * sizeof(tsKMerLoci);
#ifdef _WIN32
		pKMerLoci = (tsKMerLoci*)realloc(m_pKMerLoci, memreq);
		if (pKMerLoci == nullptr)
			{
#else
			pKMerLoci = (tsKMerLoci*)mremap(m_pKMerLoci, m_AllocdKMerLociMem, memreq, MREMAP_MAYMOVE);
			if (pKMerLoci == MAP_FAILED)
				{
#endif
				ReleaseFastSerialise();
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddKMerLoci: Memory reallocation to %zd bytes failed - %s", (int64_t)memreq, strerror(errno));
				return(0);
				}
			m_pKMerLoci = pKMerLoci;
			m_AllocdKMerLociMem = memreq;
			m_AllocdKMerLoci = ToAllocdKMerLoci;
			}
	pKMerLoci = &m_pKMerLoci[m_UsedKMerLoci++];
if (pInitKMerLoci != nullptr)
	*pKMerLoci = *pInitKMerLoci;
else
	memset(pKMerLoci, 0, sizeof(tsKMerLoci));
ReleaseFastSerialise();
return(m_UsedKMerLoci);
}

int32_t									// returned index+1 into m_pHGBinSpecs[] to allocated and initialised allele stack, 0 if errors
CCallHaplotypes::AddHGBinSpec(int32_t ChromID,	    // haplotype group clustering is on this chromosome
	int32_t ChromSize,			// chromosome size
	int32_t StartLoci,      // bin starts at this loci
	int32_t NumLoci,	    // covering this many loci
	int32_t MinCentroidDistance, // centroid distance can range from this minimum and up to MaxCentroidDistance when attempting to limit number of haplotype groupings to MaxNumHaplotypeGroups
	int32_t MaxCentroidDistance, // maximum centroid distance which can be used when attempting to meet MaxNumHaplotypeGroups constraint
	int32_t MaxNumHaplotypeGroups) // attempt to constrain number of haplotype groups to be this maximum

{
tsHGBinSpec Tmp;
if(StartLoci < 0 || NumLoci < 1 || (StartLoci+ NumLoci) > ChromSize) // sanity check!
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddDGTLoci: Bin loci range (StartLoci: %d NumLoci: %d) extends past chromosome size: %d for ChromID: %d", StartLoci,NumLoci,ChromSize, ChromID);
	return(0);
	}
Tmp.ChromID = ChromID;
Tmp.ChromSize = ChromSize;
Tmp.MaxCentroidDistance = MaxCentroidDistance;
Tmp.NumLoci = NumLoci;
Tmp.StartLoci = StartLoci;
Tmp.MinCentroidDistance = MinCentroidDistance;
Tmp.MaxCentroidDistance = MaxCentroidDistance;
Tmp.MaxNumHaplotypeGroups = MaxNumHaplotypeGroups;
Tmp.ActualHaplotypeGroups = 0;
Tmp.BinID = 0;
Tmp.ProcState = 0; // bin processing state  - 0x00 if unprocessed, 0x01 if currently being processing, 0x03 if processing completed but referenced chromosome not located, 0x07 if processing successfully completed
Tmp.pHaplotypeGroup = nullptr;
return(AddHGBinSpec(&Tmp));
}

int32_t									// returned index+1 into m_pHGBinSpecs[] to allocated and initialised allele stack, 0 if errors
CCallHaplotypes::AddHGBinSpec(tsHGBinSpec* pInitHGBinSpec)	// allocated tsHGBinSpec to be initialised with a copy of pInitHGBinSpec
{
uint32_t ToAllocdHGBinSpecs;
tsHGBinSpec* pHGBinSpec;
size_t memreq;
AcquireFastSerialise();
if(m_pHGBinSpecs == nullptr)					// may be nullptr first time in
	{
	memreq = (size_t)cAllocHGBinSpecs * sizeof(tsHGBinSpec);
#ifdef _WIN32
	m_pHGBinSpecs = (tsHGBinSpec*)malloc((size_t)memreq);
	if(m_pHGBinSpecs == nullptr)
		{
		ReleaseFastSerialise();
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddHGBinSpec: Memory allocation of %zd bytes failed", (int64_t)memreq);
		return(0);
		}
#else
	m_pHGBinSpecs = (tsHGBinSpec*)mmap(nullptr, (size_t)memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	if(m_pHGBinSpecs == MAP_FAILED)
		{
		m_pHGBinSpecs = nullptr;
		ReleaseFastSerialise();
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddHGBinSpec: Memory allocation of %zd bytes through mmap()  failed", (int64_t)memreq, strerror(errno));
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
		if(pHGBinSpec == nullptr)
			{
#else
		pHGBinSpec = (tsHGBinSpec*)mremap(m_pHGBinSpecs, m_AllocdHGBinSpecsMem, memreq, MREMAP_MAYMOVE);
		if(pHGBinSpec == MAP_FAILED)
			{
#endif
			ReleaseFastSerialise();
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddHGBinSpec: Memory reallocation to %zd bytes failed - %s", (int64_t)memreq, strerror(errno));
			return(0);
			}
		m_pHGBinSpecs = pHGBinSpec;
		m_AllocdHGBinSpecsMem = memreq;
		m_AllocdHGBinSpecs = ToAllocdHGBinSpecs;
		}
pHGBinSpec = &m_pHGBinSpecs[m_UsedHGBinSpecs++];
if(pInitHGBinSpec != nullptr)
	*pHGBinSpec = *pInitHGBinSpec;
else
	memset(pHGBinSpec, 0, sizeof(tsHGBinSpec));
pHGBinSpec->BinID = m_UsedHGBinSpecs;
ReleaseFastSerialise();
return(m_UsedHGBinSpecs);
}

tsHGBinSpec*    // returned haplotype grouping bin specification, returns nullptr if all bins have been iterated
CCallHaplotypes::IterateHGBinSpecs(int32_t PrevBinID,    // previously iterated bin identifier, to start from 1st bin then pass 0 as the bin identifier 
								   int32_t ChromID) // only processing bins on this chrom, 0 if processing all bins 
{
tsHGBinSpec* pHGBinSpec;

AcquireFastSerialise();
while(PrevBinID < m_UsedHGBinSpecs)
	{
	pHGBinSpec = &m_pHGBinSpecs[PrevBinID++]; // note: bin specs are unordered so need to iterate over all bins so as to ensure all bins on requested chrom are processed 
	if(ChromID != 0 && pHGBinSpec->ChromID != ChromID) // needing to match on ChromID ..
		{
		if(ChromID > pHGBinSpec->ChromID) // BinIDs ordered by ChromID ascending so checking if iterated past bins for targeted chrom
			break;
		continue;
		}
	// can accept bin spec by ChromID but bin may already have been allocated for processing ...
	if(!(pHGBinSpec->ProcState & 0x07)) // bin processing state  - 0x00 if unprocessed, 0x01 if currently being processing, 0x03 if processing completed but referenced chromosome not located, 0x07 if processing successfully completed
		{
		pHGBinSpec->ProcState |= 0x01; // bin now allocated for processing by caller
		m_ProccessingHGBinID = pHGBinSpec->BinID; // marks last bin accepted for processing, will be used to initialise processing for next chromosome
		ReleaseFastSerialise();
		return(pHGBinSpec);
		}
	}
ReleaseFastSerialise();
return(nullptr);
}


	// load and parse haplotype grouping bin specification file
    // This file specifies grouping bins for specific regions on chromosomes starting at requested loci of specified length. The number of bins in that region is specified and required maximal clustering distance
    // Thus regions of interest can easily be covered without having to specify explicitly each bin in that region
int								 // eBSFSuccess or error
CCallHaplotypes::LoadHHGBinSpecs(char* pszInFile)  // load grouping specifications from this file

{
int Rslt;
uint32_t CurLineNumber;
int32_t NumFields;
int32_t RegionStartLoci;
uint32_t NumBins;
int32_t BinSize;
int32_t RegionLen;
int32_t MinCentroidDistance; // centroid distance can range from this minimum and up to MaxCentroidDistance when attempting to limit number of haplotype groupings to MaxNumHaplotypeGroups
int32_t MaxCentroidDistance; // maximum centroid distance which can be used when attempting to meet MaxNumHaplotypeGroups constraint
int32_t MaxNumHaplotypeGroups;// attempt to constrain number of haplotype groups to be this maximum
char *pszChrom;
int32_t ChromID;
int32_t CurChromID;
tsCHChromMetadata* pChromMetadata;
uint32_t NumUnrecognisedChroms;
uint32_t NumRegionErrs;
uint32_t NumRegions;
uint32_t NumHGBinSpecs;

if(m_pCSVFile != nullptr) // shouldn't have been instantiated, but better to be sure!
	{
	delete m_pCSVFile;
	m_pCSVFile = nullptr;
	}
if((m_pCSVFile = new CCSVFile) == nullptr)
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
pChromMetadata = nullptr;
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



	if(pszChrom == nullptr || pszChrom[0] == '\0')
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

	m_pCSVFile->GetUint(4, &NumBins);
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
	if((RegionLen + RegionStartLoci) > pChromMetadata->ChromLen)
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
		if(AddHGBinSpec(ChromID, pChromMetadata->ChromLen,RegionStartLoci, BinSize, MinCentroidDistance,MaxCentroidDistance,MaxNumHaplotypeGroups) < 1)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddHGBinSpec(%d,%d,%d,%d,%d,%d) returned unrecoverable memory allocation error,parsing file '%s' at line %d", ChromID, pChromMetadata->ChromLen, RegionStartLoci, BinSize, MinCentroidDistance,MaxCentroidDistance,MaxNumHaplotypeGroups, pszInFile, CurLineNumber);
			Reset();
			return(eBSFerrParse);
			}
		RegionStartLoci += BinSize;
		NumHGBinSpecs++;
		}
	NumRegions++;
	}
delete m_pCSVFile;
m_pCSVFile = nullptr;
gDiagnostics.DiagOut(eDLInfo, gszProcName, "Parsed and accepted '%d' haplotype group region specifications resulting in '%d' bins for processing, sloughed '%d' regions, in file '%s', total accumulated bins '%d'",NumRegions, NumHGBinSpecs,NumUnrecognisedChroms+NumRegionErrs, pszInFile, m_UsedHGBinSpecs);
return(m_UsedHGBinSpecs);
}


// load and parse actual haplotype groupings previously generated by 'ngskit4b callhaplotypes'
// some complications with the format of the haplotype groupings file from earlier releases of 'ngskit4b callhaplotypes' are -
// a) field containing the experiment identifier may be missing
// b) field containing the row identifier may be missing
// A check is made for these fields and if missing then defaults are used - experiment identifier: 1
//                                                                        - row identifier: starts at 1
int								 // eBSFSuccess or error
CCallHaplotypes::LoadHaplotypeGroupings(char* pszInFile)  // load previously generated haplotype groupings from this CSV file

{
int Rslt;
uint32_t RowExprID;
uint32_t RowRowID;

uint32_t CurLineNumber;
int32_t NumFields;
int32_t ExpNumFields;
int32_t StartLoci;
int32_t Len;
int32_t MinDistance; // centroid distance can range from this minimum and up to MaxCentroidDistance when attempting to limit number of haplotype groupings to MaxNumHaplotypeGroups
int32_t MaxDistance; // maximum centroid distance which can be used when attempting to meet MaxNumHaplotypeGroups constraint
uint32_t MaxHaplotypeGroups;// attempt to constrain number of haplotype groups to be this maximum
char *pszChrom;
uint32_t ChromID;
uint32_t CurChromID;
tsCHChromMetadata* pChromMetadata;
uint32_t NumUnrecognisedChroms;
uint32_t NumBinSpecErrs;
uint32_t NumAcceptedBins;
uint32_t NumHGBinSpecs;
uint32_t NumSamples;
uint32_t SampleIDIdx;
uint32_t SampleID;
char *pszSampleID;

uint32_t ActualDistance;
uint32_t ActualHaplotypeGroups;

uint32_t RelBinFieldsStart;
bool bExprIDInCSV;
bool bRowIDInCSV;

int32_t *pFieldSampleIDMappings;

if(m_pCSVFile != nullptr) 
	{
	delete m_pCSVFile;
	m_pCSVFile = nullptr;
	}
if((m_pCSVFile = new CCSVFile) == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to instantiate CCSVfile");
	Reset();
	return(eBSFerrObj);
	}
m_pCSVFile->SetMaxFields((cMaxFounderReadsets * 2) + 11); // could be many, many, CSV fields!

if((Rslt = m_pCSVFile->Open(pszInFile)) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to open file: '%s'", pszInFile);
	Reset();
	return(Rslt);
	}

NumSamples = 0;
NumHGBinSpecs = 0;
CurLineNumber = 0;
CurChromID = 0;
pChromMetadata = nullptr;
NumAcceptedBins = 0;
NumUnrecognisedChroms = 0;
NumBinSpecErrs = 0;
ExpNumFields = 0;
pFieldSampleIDMappings = nullptr;
RelBinFieldsStart = 0;
bExprIDInCSV = false;
bRowIDInCSV = false;
RowRowID = 0;
RowExprID = 1;
while((Rslt = m_pCSVFile->NextLine()) > 0)		// onto next line containing fields
	{
	CurLineNumber++;
	if((NumFields = m_pCSVFile->GetCurFields()) < 11)	// must contain at least 9 fields - would be the case if there was a single sample which was grouped into a single group with missing ExprID and RowID
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Haplotype grouping file '%s' expected to contain a minimum of 11 fields, it contains %d at line %d", pszInFile, NumFields, CurLineNumber);
		if(pFieldSampleIDMappings != nullptr)
			delete[]pFieldSampleIDMappings;
		Reset();
		return(eBSFerrParse);
		}
	if(ExpNumFields > 0 && NumFields != ExpNumFields)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Haplotype grouping file '%s' expected to contain %d fields, it contains %d at line %d", pszInFile,ExpNumFields, NumFields, CurLineNumber);
		if(pFieldSampleIDMappings != nullptr)
			delete[]pFieldSampleIDMappings;
		Reset();
		return(eBSFerrParse);
		}

	// 1st row must be a header row, need to parse this out to determine if ExprID, RowID are missing and to obtain the sample identifiers
	// sample identifiers are those starting after field 8 and terminating immediately prior to header field name prefixed as 'GrpMembers:"
	if(CurLineNumber == 1 && m_pCSVFile->IsLikelyHeaderLine())
		{
		NumSamples = 0;
		pFieldSampleIDMappings = new int32_t[cMaxFounderReadsets];
		RelBinFieldsStart = 0;
		m_pCSVFile->GetText(1, &pszSampleID);
		if(!strnicmp(pszSampleID, "ExprID", 6))
			{
			RelBinFieldsStart = 1; // CSV has field for ExprID
			bExprIDInCSV = true;
			}
		m_pCSVFile->GetText(RelBinFieldsStart+1, &pszSampleID);
		if(!strnicmp(pszSampleID, "RowID", 5))
			{
			RelBinFieldsStart++; // CSV has field for ExprID
			bRowIDInCSV = true;
			}

		for(SampleIDIdx = 0; SampleIDIdx < NumFields - (9+RelBinFieldsStart); SampleIDIdx++)
			{
			m_pCSVFile->GetText(RelBinFieldsStart+SampleIDIdx+9, &pszSampleID);
			if(!strnicmp(pszSampleID, "GrpMembers:",11))
				break;
			// all samples must have been loaded as PBAs
			if((SampleID = LocateReadset(pszSampleID, 0)) == 0)
				{
				delete []pFieldSampleIDMappings;
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unrecognised sample identifier '%s' referenced in header line of haplotype grouping file '%s'",pszSampleID, pszInFile);
				Reset();
				return(eBSFerrParse);
				}
			pFieldSampleIDMappings[SampleIDIdx] = SampleID;
			NumSamples++;
			}
		ExpNumFields = NumFields;  // must be consistency, remaining rows should have same number of fields
		continue;
		}
	// past the header, parse out individual bin rows
	RowRowID++; // if RowID actually in the CSV then that value will overwrite this...
	if(bExprIDInCSV)
		{
		m_pCSVFile->GetUint(1, &RowExprID);
		if(bRowIDInCSV)
			m_pCSVFile->GetUint(2, &RowRowID);
		}
	else
		if(bRowIDInCSV)
			m_pCSVFile->GetUint(1, &RowRowID);

	m_pCSVFile->GetText(RelBinFieldsStart+1, &pszChrom);

	if(pszChrom == nullptr || pszChrom[0] == '\0')
		{
		if(++NumUnrecognisedChroms <= 100)
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "Unrecognised chromosome/sequence name at line %d in file '%s', bin sloughed", CurLineNumber, pszInFile);
		continue;
		}

	if(!AcceptThisChromName(pszChrom)) // may not be accepting this chromosome for processing
		continue;

	// is chromosome known?
	if((ChromID = LocateChrom(pszChrom)) <= 0)
		{
		if(++NumUnrecognisedChroms <= 100)
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "Unrecognised chromosome/sequence name '%s' at line %d in file '%s', bin sloughed", pszChrom, CurLineNumber, pszInFile);
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
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "Chromosome '%s' with length '%d' is smaller than minimum sized bin of 10 at line %d in file '%s', bin sloughed", pszChrom, pChromMetadata->ChromLen, CurLineNumber, pszInFile);
		continue;
		}

	m_pCSVFile->GetInt(RelBinFieldsStart+2, &StartLoci);
	if(StartLoci < 0 || (StartLoci + 10 > pChromMetadata->ChromLen))
		{
		if(++NumBinSpecErrs <= 100)
			{
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "Chromosome/sequence '%s' has size: %d, referenced in file '%s' at line: %d", pszChrom, pChromMetadata->ChromLen, pszInFile, CurLineNumber);
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "Start loci '%d' must be in range 0..'%d' (minimum bin size is 10), bin sloughed", StartLoci, pChromMetadata->ChromLen - 10);
			}
		continue;
		}

	m_pCSVFile->GetInt(RelBinFieldsStart+3, &Len);
	if(Len < 10 || (Len + StartLoci > pChromMetadata->ChromLen))    // lower limit of 10 for bin len
		{
		if(++NumBinSpecErrs <= 100)
			{
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "Chromosome/sequence '%s' has size: %d, referenced in file '%s' at line: %d", pszChrom, pChromMetadata->ChromLen, pszInFile, CurLineNumber);
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "Bin Len '%d' starting from %d is less than 10 or would extend past end of chromosome, bin sloughed", Len, StartLoci);
			}
		continue;
		}


	m_pCSVFile->GetInt(RelBinFieldsStart+4, &MinDistance);
	if(MinDistance < 1)
		{
		if(++NumBinSpecErrs <= 100)
			{
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "Chromosome/sequence '%s' has size: %d, referenced in file '%s' at line: %d", pszChrom, pChromMetadata->ChromLen, pszInFile, CurLineNumber);
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "MinDistance '%d' must be at least 1 and <= Len '%d', bin sloughed", MinDistance,Len);
			}
		continue;
		}
	m_pCSVFile->GetInt(RelBinFieldsStart+5, &MaxDistance);
	if(MaxDistance  < MinDistance || MaxDistance > Len)
		{
		if(++NumBinSpecErrs <= 100)
			{
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "Chromosome/sequence '%s' has size: %d, referenced in file '%s' at line: %d", pszChrom, pChromMetadata->ChromLen, pszInFile, CurLineNumber);
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "MaxDistance '%d' must be >= MinDistance '%d' and <= Len '%d', bin sloughed",MaxDistance,MinDistance, Len);
			}
		continue;
		}

	m_pCSVFile->GetUint(RelBinFieldsStart+6, &MaxHaplotypeGroups);


	uint32_t BinSpecID;
	BinSpecID = AddHGBinSpec(ChromID, pChromMetadata->ChromLen,StartLoci,Len,MinDistance,MaxDistance,MaxHaplotypeGroups);
	tsHGBinSpec *pHGBinSpec = &m_pHGBinSpecs[BinSpecID - 1];

	// now can parse out the HapGrp for each sample, then initialise an instance of tsHaplotypeGroup
	m_pCSVFile->GetUint(RelBinFieldsStart+7, &ActualDistance);
	m_pCSVFile->GetUint(RelBinFieldsStart+8, &ActualHaplotypeGroups);

	tsHaplotypeGroup *pHaplotypeGroup;
	size_t MemReq = sizeof(tsHaplotypeGroup) + (sizeof(tsBitsVect) * (ActualHaplotypeGroups-1)); // tsHaplotypeGroup already includes a HaplotypeGroup tsBitsVect
	pHaplotypeGroup = (tsHaplotypeGroup *) new uint8_t[MemReq];
	memset(pHaplotypeGroup, 0, MemReq);
	pHaplotypeGroup->Size = (uint32_t)MemReq;

	uint32_t HapGrp;

	for(SampleIDIdx = 0; SampleIDIdx < NumSamples; SampleIDIdx++)
		{
		m_pCSVFile->GetUint(RelBinFieldsStart+SampleIDIdx+9, &HapGrp);
		BitsVectSet(pFieldSampleIDMappings[SampleIDIdx]-1, pHaplotypeGroup->HaplotypeGroup[HapGrp - 1]);
		}

	pHaplotypeGroup->SrcExprID = RowExprID;
	pHaplotypeGroup->SrcRowID = RowRowID;
	pHaplotypeGroup->CentroidDistance = ActualDistance;
	pHaplotypeGroup->ChromID = ChromID;
	pHaplotypeGroup->ChromSize = pChromMetadata->ChromLen;
	pHaplotypeGroup->MaxCentroidDistance = MaxDistance;
	pHaplotypeGroup->MinCentroidDistance = MinDistance;
	pHaplotypeGroup->MRAGrpChromID = 0;
	pHaplotypeGroup->MRAGrpLoci = -1;
	memset(pHaplotypeGroup->MRAGrpConsensus,0x0ff,sizeof(pHaplotypeGroup->MRAGrpConsensus));
	pHaplotypeGroup->NumFndrs = NumSamples;
	pHaplotypeGroup->MaxHaplotypeGroups = MaxHaplotypeGroups;
	pHaplotypeGroup->ActualHaplotypeGroups = ActualHaplotypeGroups;
	pHaplotypeGroup->NumLoci = Len;
	pHaplotypeGroup->StartLoci = StartLoci;
	pHGBinSpec->ChromID = ChromID;
	pHGBinSpec->ChromSize = pChromMetadata->ChromLen;
	pHGBinSpec->StartLoci = StartLoci;
	pHGBinSpec->NumLoci = Len;
	pHGBinSpec->MinCentroidDistance = MinDistance;
	pHGBinSpec->MaxCentroidDistance = MaxDistance;
	pHGBinSpec->MaxNumHaplotypeGroups = MaxHaplotypeGroups;
	pHGBinSpec->ActualHaplotypeGroups = ActualHaplotypeGroups;
	pHGBinSpec->pHaplotypeGroup = pHaplotypeGroup;
	pHGBinSpec->ProcState = 0x07; // bin processing state  - 0x00 if unprocessed, 0x01 if currently being processing, 0x03 if processing completed but referenced chromosome not located, 0x07 if processing successfully completed
	NumAcceptedBins++;
	}
delete m_pCSVFile;
m_pCSVFile = nullptr;
if(pFieldSampleIDMappings != nullptr)
	delete[]pFieldSampleIDMappings;
gDiagnostics.DiagOut(eDLInfo, gszProcName, "Parsed and accepted '%d' haplotype grouping bins for processing, sloughed '%d' bins, in file '%s', total accumulated bins '%d'", NumAcceptedBins,NumUnrecognisedChroms+NumBinSpecErrs, pszInFile, m_UsedHGBinSpecs);
return(m_UsedHGBinSpecs);
}



// When generating allele scoring based haplotype groupings then the reference groupings are used to determine the source groupings
// Effectively the references are used as proxies for the sources
//
int								 // eBSFSuccess or error
CCallHaplotypes::LoadRefGroupings(char* pszInFile)  // load subset of previously generated haplotype groupings from this CSV file

{
int Rslt;
uint32_t RowExprID;
uint32_t RowRowID;

uint32_t CurLineNumber;
int32_t NumFields;
int32_t ExpNumFields;
int32_t StartLoci;
int32_t Len;
int32_t MinDistance; // centroid distance can range from this minimum and up to MaxCentroidDistance when attempting to limit number of haplotype groupings to MaxNumHaplotypeGroups
int32_t MaxDistance; // maximum centroid distance which can be used when attempting to meet MaxNumHaplotypeGroups constraint
uint32_t MaxHaplotypeGroups;// attempt to constrain number of haplotype groups to be this maximum
char* pszChrom;
uint32_t ChromID;
uint32_t CurChromID;
int32_t CurChromSize;
uint32_t NumUnrecognisedChroms;
uint32_t NumBinSpecErrs;
uint32_t NumAcceptedBins;
uint32_t NumHGBinSpecs;
uint32_t NumSamples;
uint32_t NumUnscoredRefs;
uint32_t SampleIDIdx;
uint32_t SampleID;
char* pszSampleID;

uint32_t ActualDistance;
uint32_t ActualHaplotypeGroups;

uint32_t RelBinFieldsStart;
bool bExprIDInCSV;
bool bRowIDInCSV;

int32_t* pFieldSampleIDMappings;

if (m_pCSVFile != nullptr)
{
	delete m_pCSVFile;
	m_pCSVFile = nullptr;
}
if ((m_pCSVFile = new CCSVFile) == nullptr)
{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadRefGroupings: Unable to instantiate CCSVfile");
	Reset();
	return(eBSFerrObj);
}
m_pCSVFile->SetMaxFields((cMaxFounderReadsets * 2) + 11); // could be many, many, CSV fields!

if ((Rslt = m_pCSVFile->Open(pszInFile)) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadRefGroupings: Unable to open file: '%s'", pszInFile);
	Reset();
	return(Rslt);
	}

NumSamples = 0;
NumUnscoredRefs = 0;
NumHGBinSpecs = 0;
CurLineNumber = 0;
CurChromID = 0;
CurChromSize = 0;
NumAcceptedBins = 0;
NumUnrecognisedChroms = 0;
NumBinSpecErrs = 0;
ExpNumFields = 0;
pFieldSampleIDMappings = nullptr;
RelBinFieldsStart = 0;
bExprIDInCSV = false;
bRowIDInCSV = false;
RowRowID = 0;
RowExprID = 1;
while ((Rslt = m_pCSVFile->NextLine()) > 0)		// onto next line containing fields
	{
	CurLineNumber++;
	if ((NumFields = m_pCSVFile->GetCurFields()) < 11)	// must contain at least 9 fields - would be the case if there was a single sample which was grouped into a single group with missing ExprID and RowID
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadRefGroupings: Haplotype grouping file '%s' expected to contain a minimum of 11 fields, it contains %d at line %d", pszInFile, NumFields, CurLineNumber);
		if (pFieldSampleIDMappings != nullptr)
			delete[]pFieldSampleIDMappings;
		Reset();
		return(eBSFerrParse);
		}
	if (ExpNumFields > 0 && NumFields != ExpNumFields)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadRefGroupings: Haplotype grouping file '%s' expected to contain %d fields, it contains %d at line %d", pszInFile, ExpNumFields, NumFields, CurLineNumber);
		if (pFieldSampleIDMappings != nullptr)
			delete[]pFieldSampleIDMappings;
		Reset();
		return(eBSFerrParse);
		}

	// 1st row must be a header row, need to parse this out to determine if ExprID, RowID are missing and to obtain the sample identifiers
	// sample identifiers are those starting after field 8 and terminating immediately prior to header field name prefixed as 'GrpMembers:"
	if (CurLineNumber == 1 && m_pCSVFile->IsLikelyHeaderLine())
		{
		pFieldSampleIDMappings = new int32_t[cMaxFounderReadsets];
		RelBinFieldsStart = 0;
		m_pCSVFile->GetText(1, &pszSampleID);
		if (!strnicmp(pszSampleID, "ExprID", 6))
			{
			RelBinFieldsStart = 1; // CSV has field for ExprID
			bExprIDInCSV = true;
			}
		m_pCSVFile->GetText(RelBinFieldsStart + 1, &pszSampleID);
		if (!strnicmp(pszSampleID, "RowID", 5))
			{
			RelBinFieldsStart++; // CSV has field for ExprID
			bRowIDInCSV = true;
			}
		for (SampleIDIdx = 0; SampleIDIdx < NumFields - (9 + RelBinFieldsStart); SampleIDIdx++)
			{
			m_pCSVFile->GetText(RelBinFieldsStart + SampleIDIdx + 9, &pszSampleID);
			if (!strnicmp(pszSampleID, "GrpMembers:", 11))
				break;
			if((SampleID = LocateReadset(pszSampleID, 0))==0)		// SampleID will be 0 if this sample was not a reference sample when scoring the source alleles vs reference alleles
				NumUnscoredRefs++;										// keeping tabs on number of samples which are not represented as scored references
			pFieldSampleIDMappings[SampleIDIdx] = SampleID;
			NumSamples++;
			}
		ExpNumFields = NumFields;  // must be consistency, remaining rows should have same number of fields
		continue;
	}
	// past the header, parse out individual bin rows
	RowRowID++; // if RowID actually in the CSV then that value will overwrite this...
	if (bExprIDInCSV)
		{
		m_pCSVFile->GetUint(1, &RowExprID);
		if (bRowIDInCSV)
			m_pCSVFile->GetUint(2, &RowRowID);
		}
	else
		if (bRowIDInCSV)
			m_pCSVFile->GetUint(1, &RowRowID);

	m_pCSVFile->GetText(RelBinFieldsStart + 1, &pszChrom);

	if (pszChrom == nullptr || pszChrom[0] == '\0')
		{
		if (++NumUnrecognisedChroms <= 100)
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "LoadRefGroupings: Unrecognised chromosome/sequence name at line %d in file '%s', bin sloughed", CurLineNumber, pszInFile);
		continue;
		}

	if (!AcceptThisChromName(pszChrom)) // may not be accepting this chromosome for processing
		continue;

	// is chromosome known?
	if ((ChromID = LocateChrom(pszChrom)) <= 0)
		{
		if (++NumUnrecognisedChroms <= 100)
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "LoadRefGroupings: Unrecognised chromosome/sequence name '%s' at line %d in file '%s', bin sloughed", pszChrom, CurLineNumber, pszInFile);
		continue;
		}

	// chromosome known so can determine it's size
	if (CurChromID != ChromID)
		{
		CurChromSize = m_ChromSizes[ChromID - 1];
		CurChromID = ChromID;
		}
	if (CurChromSize < 10)
		{
		if (++NumUnrecognisedChroms <= 100)
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "LoadRefGroupings: Chromosome '%s' with length '%d' is smaller than minimum sized bin of 10 at line %d in file '%s', bin sloughed", pszChrom, CurChromSize, CurLineNumber, pszInFile);
		continue;
		}

	m_pCSVFile->GetInt(RelBinFieldsStart + 2, &StartLoci);
	if (StartLoci < 0 || (StartLoci + 10 > CurChromSize))
		{
		if (++NumBinSpecErrs <= 100)
			{
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "LoadRefGroupings: Chromosome/sequence '%s' has size: %d, referenced in file '%s' at line: %d", pszChrom, CurChromSize, pszInFile, CurLineNumber);
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "LoadRefGroupings: Start loci '%d' must be in range 0..'%d' (minimum bin size is 10), bin sloughed", StartLoci, CurChromSize - 10);
			}
		continue;
		}

	m_pCSVFile->GetInt(RelBinFieldsStart + 3, &Len);
	if (Len < 10 || (Len + StartLoci > CurChromSize))    // lower limit of 10 for bin len
		{
		if (++NumBinSpecErrs <= 100)
			{
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "LoadRefGroupings: Chromosome/sequence '%s' has size: %d, referenced in file '%s' at line: %d", pszChrom, CurChromSize, pszInFile, CurLineNumber);
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "LoadRefGroupings: Bin Len '%d' starting from %d is less than 10 or would extend past end of chromosome, bin sloughed", Len, StartLoci);
			}
		continue;
		}


	m_pCSVFile->GetInt(RelBinFieldsStart + 4, &MinDistance);
	if (MinDistance < 1)
		{
		if (++NumBinSpecErrs <= 100)
			{
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "LoadRefGroupings: Chromosome/sequence '%s' has size: %d, referenced in file '%s' at line: %d", pszChrom, CurChromSize, pszInFile, CurLineNumber);
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "LoadRefGroupings: MinDistance '%d' must be at least 1 and <= Len '%d', bin sloughed", MinDistance, Len);
			}
		continue;
		}
	m_pCSVFile->GetInt(RelBinFieldsStart + 5, &MaxDistance);
	if (MaxDistance  < MinDistance || MaxDistance > Len)
		{
		if (++NumBinSpecErrs <= 100)
			{
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "LoadRefGroupings: Chromosome/sequence '%s' has size: %d, referenced in file '%s' at line: %d", pszChrom, CurChromSize, pszInFile, CurLineNumber);
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "LoadRefGroupings: MaxDistance '%d' must be >= MinDistance '%d' and <= Len '%d', bin sloughed", MaxDistance, MinDistance, Len);
			}
		continue;
		}

	m_pCSVFile->GetUint(RelBinFieldsStart + 6, &MaxHaplotypeGroups);

	// now can parse out the HapGrp for each sample, then initialise an instance of tsHaplotypeGroup
	m_pCSVFile->GetUint(RelBinFieldsStart + 7, &ActualDistance);
	m_pCSVFile->GetUint(RelBinFieldsStart + 8, &ActualHaplotypeGroups);

	tsASRefGrp* pASRefGrp;
	pASRefGrp = AllocASRefGrps(NumSamples);
	uint32_t HapGrp;
	for (SampleIDIdx = 0; SampleIDIdx < NumSamples; SampleIDIdx++)
		{
		if(pFieldSampleIDMappings[SampleIDIdx] == 0)
			continue;
		m_pCSVFile->GetUint(RelBinFieldsStart + SampleIDIdx + 9, &HapGrp);
		pASRefGrp->RefIDGrps[SampleIDIdx].RefID = pFieldSampleIDMappings[SampleIDIdx];
		pASRefGrp->RefIDGrps[SampleIDIdx].Grp = HapGrp;
		}

	pASRefGrp->BinID = RowRowID;
	pASRefGrp->ChromID = ChromID;
	pASRefGrp->ChromSize = CurChromSize;
	pASRefGrp->MaxCentroidDistance = MaxDistance;
	pASRefGrp->MinCentroidDistance = MinDistance;
	pASRefGrp->ActualCentroidDistance = ActualDistance;
	pASRefGrp->NumRefGrps = NumSamples;
	pASRefGrp->MaxHaplotypeGroups = MaxHaplotypeGroups;
	pASRefGrp->ActualHaplotypeGroups = ActualHaplotypeGroups;      // actual number of haplotype groups = ActualHaplotypeGroups;
	pASRefGrp->NumLoci = Len;
	pASRefGrp->StartLoci = StartLoci;
	NumAcceptedBins++;
	}
delete m_pCSVFile;
m_pCSVFile = nullptr;
if (pFieldSampleIDMappings != nullptr)
	delete[]pFieldSampleIDMappings;
m_AcceptedASRefGrps = NumAcceptedBins;

if (m_AcceptedASRefGrps > 1)
	m_mtqsort.qsort(m_pASRefGrps, m_AcceptedASRefGrps, m_pASRefGrps->Size, SortASRefGrp);

// checking for overlapping bins and full chromosome bin coverage
tsASRefGrp* pASRefGrpBin = m_pASRefGrps;
int32_t CurBinLoci;
int32_t ExpNxtBinLoci;
int32_t CurBinSize;
int32_t BinIdx;
CurChromID = 0;


if(NumAcceptedBins == 0)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadRefGroupings: Unable to accept any haplotype groupings in file '%s'", pszInFile);
	return(eBSFerrParse);
	}


// basic validation of grouping bins
CurChromID = 0;
for (BinIdx = 0; BinIdx < m_AcceptedASRefGrps; BinIdx++)
	{
	pASRefGrpBin = (tsASRefGrp*)((uint8_t*)m_pASRefGrps + (BinIdx * pASRefGrpBin->Size));
	if (pASRefGrpBin->ChromID != CurChromID) // first bin this chrom
	{
		if (CurChromID != 0)
			{
			if (ExpNxtBinLoci != CurChromSize)
				break;
			}
		CurChromID = pASRefGrpBin->ChromID;
		CurBinSize = pASRefGrpBin->NumLoci;
		ExpNxtBinLoci = pASRefGrpBin->NumLoci;
		CurChromSize = pASRefGrpBin->ChromSize;
		CurBinLoci = 0;
		continue;
		}
	if (pASRefGrpBin->StartLoci != CurBinLoci)	// new bin loci on same chrom
		{
		if (CurChromSize != pASRefGrpBin->ChromSize)				// same chromosome so expecting to be same size except if last bin on the chromosome
			break;
		if (ExpNxtBinLoci != pASRefGrpBin->StartLoci)				// loci must be same as expected
			break;
		CurBinLoci = pASRefGrpBin->StartLoci;
		CurBinSize = pASRefGrpBin->NumLoci;
		ExpNxtBinLoci += CurBinSize;
		continue;
		}
	break;
	}

if (BinIdx < m_AcceptedASRefGrps || ExpNxtBinLoci != CurChromSize)
{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadRefGroupings: Inconsistencies in file '%s'", pszInFile);
	return(eBSFerrParse);
}


gDiagnostics.DiagOut(eDLInfo, gszProcName, "LoadRefGroupings: Parsed and accepted '%d' reference haplotype grouping bins for processing from file '%s'", NumAcceptedBins, pszInFile, m_AcceptedASRefGrps);
return(NumAcceptedBins);
}


tsASRefGrp *
CCallHaplotypes::LocateASRefGrp(int32_t CurChromID,		// grouping is on this chromosome
	int32_t CurBinLoci)			// in bin starting at this loci;
{
int32_t Idx;
static tsASRefGrp* pMRAASRefGrp = nullptr;	// most recently accessed ASRefGrp
if (m_pASRefGrps == nullptr || m_AcceptedASRefGrps == 0)
	return(nullptr);
if (pMRAASRefGrp == nullptr || pMRAASRefGrp->ChromID != CurChromID || pMRAASRefGrp->StartLoci != CurBinLoci)
	{
	pMRAASRefGrp = m_pASRefGrps;
	for (Idx = 0; Idx < m_AcceptedASRefGrps; Idx++)
		{
		if (pMRAASRefGrp->ChromID == CurChromID && pMRAASRefGrp->StartLoci == CurBinLoci)
			break;
		pMRAASRefGrp = (tsASRefGrp*)((uint8_t*)pMRAASRefGrp + pMRAASRefGrp->Size);
		}
	if (Idx == m_AcceptedASRefGrps)
		{
		pMRAASRefGrp = nullptr;
		return(nullptr);
		}
	}
return(pMRAASRefGrp);
}

tsASRefIDGrp *
CCallHaplotypes::LocateASRefIDGrp(int32_t CurChromID,		// grouping is on this chromosome
				int32_t CurBinLoci,			// in bin starting at this loci
				int32_t RefID)				// grouping for this reference
{
int32_t Idx;
static tsASRefGrp *pMRAASRefGrp = nullptr;	// most recently accessed ASRefGrp
tsASRefIDGrp *pASRefIDGrp = nullptr;
if(m_pASRefGrps == nullptr || m_AcceptedASRefGrps == 0)
	return(nullptr);

if(pMRAASRefGrp == nullptr || pMRAASRefGrp->ChromID != CurChromID || pMRAASRefGrp->StartLoci != CurBinLoci)
	{
	pMRAASRefGrp = m_pASRefGrps;
	for(Idx = 0; Idx < m_AcceptedASRefGrps; Idx++)
		{
		if(pMRAASRefGrp->ChromID == CurChromID && pMRAASRefGrp->StartLoci == CurBinLoci)
			break;
		pMRAASRefGrp = (tsASRefGrp *)((uint8_t *)pMRAASRefGrp + pMRAASRefGrp->Size);
		}
	if(Idx == m_AcceptedASRefGrps)
		{
		pMRAASRefGrp = nullptr;
		return(nullptr);
		}
	}

pASRefIDGrp = pMRAASRefGrp->RefIDGrps;
for (Idx = 0; Idx < pMRAASRefGrp->NumRefGrps; Idx++, pASRefIDGrp++)
	{
	if(pASRefIDGrp->RefID == RefID)
		return(pASRefIDGrp);
	}
return(nullptr);
}

tsASRefGrp *									// returned allocation
CCallHaplotypes::AllocASRefGrps(int NumRefGrps)	// allocated tsASRefGrp allocated to hold this many reference sample groups
{
tsASRefGrp* pASRefGrp;
size_t memreq;
size_t ASRefGrpSize;
ASRefGrpSize = sizeof(tsASRefGrp) + (sizeof(tsASRefIDGrp) * (NumRefGrps - 1));
AcquireFastSerialise();
if (m_pASRefGrps == nullptr)					// may be nullptr first time in
	{
	memreq = (size_t)cAllocASRefGrps * ASRefGrpSize;
#ifdef _WIN32
	m_pASRefGrps = (tsASRefGrp *)malloc(memreq);
	if (m_pASRefGrps == nullptr)
		{
		ReleaseFastSerialise();
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AllocASRefGrp: Memory allocation of %zd bytes failed", (int64_t)memreq);
		return(nullptr);
		}
#else
	m_pASRefGrps = (tsASRefGrp*)mmap(nullptr, memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	if (m_pASRefGrps == MAP_FAILED)
		{
		m_pASRefGrps = nullptr;
		ReleaseFastSerialise();
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AllocASRefGrp: Memory allocation of %zd bytes through mmap()  failed", (int64_t)memreq, strerror(errno));
		return(nullptr);
		}
#endif
	m_AllocdASRefGrpMem = memreq;
	m_AllocASRefGrps = cAllocASRefGrps;
	m_UsedASRefGrps = 0;
	}
else
	// needing to allocate more memory?
	if (((m_UsedASRefGrps + 1) * ASRefGrpSize) >= m_AllocdASRefGrpMem)
		{
		m_AllocASRefGrps = m_UsedASRefGrps + cReallocASRefGrps;
		size_t memreq = (size_t)m_AllocASRefGrps * ASRefGrpSize;
#ifdef _WIN32
		pASRefGrp = (tsASRefGrp*)realloc(m_pASRefGrps, memreq);
		if (pASRefGrp == nullptr)
			{
#else
		pASRefGrp = (tsASRefGrp*)mremap(m_pASRefGrps, m_AllocdASRefGrpMem, memreq, MREMAP_MAYMOVE);
		if (pASRefGrp == MAP_FAILED)
			{
#endif
			ReleaseFastSerialise();
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "AllocASRefGrp: Memory reallocation to %zd bytes failed - %s", (int64_t)memreq, strerror(errno));
			return(0);
			}
		m_pASRefGrps = pASRefGrp;
		m_AllocdASRefGrpMem = memreq;
		m_AllocASRefGrps = m_AllocASRefGrps;
		}
pASRefGrp = (tsASRefGrp *)((uint8_t *)m_pASRefGrps + (m_UsedASRefGrps++ * ASRefGrpSize));
memset(pASRefGrp, 0, ASRefGrpSize);
pASRefGrp->Size = (int32_t)ASRefGrpSize;
ReleaseFastSerialise();
return(pASRefGrp);
}




uint32_t									// returned index+1 into m_pAlleleStacks[] to allocated and initialised allele stack, 0 if errors
CCallHaplotypes::AddAlleleStack(tsAlleleStack* pInitAlleleStack)	// allocated tsAlleleStack to be initialised with a copy of pInitAlleleStack
{
uint32_t ToAllocdAlleleStacks;
tsAlleleStack* pAlleleStack;
size_t memreq;
AcquireFastSerialise();
if(m_pAlleleStacks == nullptr)					// may be nullptr first time in
	{
	memreq = (size_t)cAllocAlleleStacks * sizeof(tsAlleleStack);
#ifdef _WIN32
	m_pAlleleStacks = (tsAlleleStack*)malloc((size_t)memreq);
	if(m_pAlleleStacks == nullptr)
		{
		ReleaseFastSerialise();
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddAlleleStack: Memory allocation of %zd bytes failed",(int64_t)memreq);
		return(0);
		}
#else
	m_pAlleleStacks = (tsAlleleStack*)mmap(nullptr, (size_t)memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	if(m_pAlleleStacks == MAP_FAILED)
		{
		m_pAlleleStacks = nullptr;
		ReleaseFastSerialise();
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddAlleleStack: Memory allocation of %zd bytes through mmap()  failed", (int64_t)memreq, strerror(errno));
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
		if(pAlleleStack == nullptr)
			{
#else
		pAlleleStack = (tsAlleleStack*)mremap(m_pAlleleStacks, m_AllocdAlleleStacksMem, memreq, MREMAP_MAYMOVE);
		if(pAlleleStack == MAP_FAILED)
			{
#endif
			ReleaseFastSerialise();
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddAlleleStack: Memory reallocation to %zd bytes failed - %s", (int64_t)memreq, strerror(errno));
			return(0);
			}
		m_pAlleleStacks = pAlleleStack;
		m_AllocdAlleleStacksMem = memreq;
		m_AllocdAlleleStacks = ToAllocdAlleleStacks;
		}
pAlleleStack = &m_pAlleleStacks[m_UsedAlleleStacks++];
if(pInitAlleleStack != nullptr)
	*pAlleleStack = *pInitAlleleStack;
else
	memset(pAlleleStack, 0, sizeof(tsAlleleStack));
ReleaseFastSerialise();
return(m_UsedAlleleStacks);
}


#ifdef _WIN32
unsigned __stdcall WorkerLoadChromPBAInstance(void * pThreadPars)
#else
void * WorkerLoadChromPBAInstance(void * pThreadPars)
#endif
{
int Rslt;
tsCHWorkerLoadChromPBAsInstance *pPars = (tsCHWorkerLoadChromPBAsInstance *)pThreadPars;			// makes it easier not having to deal with casts!
CCallHaplotypes *pWorkerInstance = (CCallHaplotypes *)pPars->pThis;

Rslt = pWorkerInstance->ProcWorkerLoadChromPBAThread(pPars);
pPars->Rslt = Rslt;
#ifdef _WIN32
_endthreadex(0);
return(eBSFSuccess);
#else
pthread_exit(&pPars->Rslt);
#endif
}


#ifdef _WIN32
unsigned __stdcall WorkerPBAInstance(void * pThreadPars)
#else
void *WorkerPBAInstance(void * pThreadPars)
#endif
{
int Rslt;
tsCHWorkerInstance *pPars = (tsCHWorkerInstance *)pThreadPars;			// makes it easier not having to deal with casts!
CCallHaplotypes *pWorkerInstance = (CCallHaplotypes *)pPars->pThis;

Rslt = pWorkerInstance->ProcWorkerThread(pPars);
pPars->Rslt = Rslt;
#ifdef _WIN32
_endthreadex(0);
return(eBSFSuccess);
#else
pthread_exit(&pPars->Rslt);
#endif
}

// initialise and start pool of worker threads to concurrently load chromosome PBAs
int				// returns 0 if all threads have started processing, 1 if all threads started and some have completed processing, 2 if not all have started
CCallHaplotypes::StartWorkerLoadChromPBAThreads(int32_t NumThreads,		// there are this many threads in pool
										int32_t StartSampleID,				// processing to start from this sample identifer
										int32_t EndSampleID,				// ending with this sample identifier inclusive
										int32_t ChromID,					// loading PBAs for this chromosome
										bool bNormAlleles)					// normalise alleles such that individual alleles can be compared without regard to the proportional coverage (0x22 -> 0x33 as an example)
{
int Rslt = eBSFSuccess;
int32_t FromSampleID;
int32_t NumSamples;
int32_t SamplesThisThread;
int32_t ThreadIdx;
tsCHWorkerLoadChromPBAsInstance *pThreadPar;

#ifndef _WIN32
// increase the default thread stack to at least cWorkThreadStackSize
size_t defaultStackSize;
pthread_attr_t threadattr;
pthread_attr_init(&threadattr);
pthread_attr_getstacksize(&threadattr, &defaultStackSize);
if (defaultStackSize < (size_t)cWorkThreadStackSize)
	{
	if ((Rslt = pthread_attr_setstacksize(&threadattr, (size_t)cWorkThreadStackSize)) != 0)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "StartWorkerLoadChromPBAThreads: pthread_attr_setstacksize(%d) failed, default was %zd", cWorkThreadStackSize, defaultStackSize);
		return(eBSFerrInternal);
		}
	}
#endif


NumSamples = 1 + EndSampleID - StartSampleID;
NumThreads = min(NumSamples, NumThreads);
FromSampleID = StartSampleID;
m_NumWorkerInsts = 0;
m_ExpNumWorkerInsts = NumThreads;
m_CompletedWorkerInsts = 0;
m_ReqTerminate = 0;
pThreadPar = m_WorkerLoadChromPBAInstances;
for (ThreadIdx = 1; ThreadIdx <= NumThreads; ThreadIdx++, pThreadPar++)
	{
	memset(pThreadPar,0,sizeof(tsCHWorkerLoadChromPBAsInstance));
	pThreadPar->bNormAlleles = bNormAlleles;
	pThreadPar->ChromID = ChromID;
	SamplesThisThread = NumSamples / (1 + NumThreads - ThreadIdx);
	NumSamples -= SamplesThisThread;
	pThreadPar->StartSampleID = FromSampleID;
	pThreadPar->EndSampleID = SamplesThisThread - 1 + FromSampleID;
	FromSampleID = pThreadPar->EndSampleID + 1;

	pThreadPar->ThreadIdx = ThreadIdx;
	pThreadPar->pThis = this;
#ifdef _WIN32
	pThreadPar->threadHandle = (HANDLE)_beginthreadex(nullptr, cWorkThreadStackSize, WorkerLoadChromPBAInstance, pThreadPar, 0, &pThreadPar->threadID);
#else
	pThreadPar->threadRslt = pthread_create(&pThreadPar->threadID, &threadattr, WorkerLoadChromPBAInstance, pThreadPar);
#endif
	}

// allow threads time to all startup
Rslt = WaitWorkerThreadStatus(cMaxWaitThreadsStartup);	// waiting untill all threads have started (some may have completed!)
if (Rslt == 2)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "StartWorkerThreads: All %d expected threads did not register as started within allowed %d seconds", NumThreads, cMaxWaitThreadsStartup);
	Rslt = TerminateWorkerThreads();
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "StartWorkerThreads: TerminateWorkThreads(60) returned %d", Rslt);
	return(eBSFerrInternal);
	}
return(Rslt);
}


// initialise and start pool of worker threads
int
CCallHaplotypes::StartWorkerThreads(int32_t NumThreads,		// there are this many threads in pool
									int32_t NumChroms)			// processing this number of chromosomes
{
int Rslt = eBSFSuccess;
int32_t StartQuerySeqIdx;
int32_t ThreadIdx;
tsCHWorkerInstance *pThreadPar;

#ifndef _WIN32
// increase the default thread stack to at least cWorkThreadStackSize
size_t defaultStackSize;
pthread_attr_t threadattr;
pthread_attr_init(&threadattr);
pthread_attr_getstacksize(&threadattr, &defaultStackSize);
if (defaultStackSize < (size_t)cWorkThreadStackSize)
	{
	if((Rslt=pthread_attr_setstacksize(&threadattr, (size_t)cWorkThreadStackSize))!=0)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "StartWorkerThreads: pthread_attr_setstacksize(%d) failed, default was %zd", cWorkThreadStackSize, defaultStackSize);
		return(eBSFerrInternal);
		}
	}
#endif

StartQuerySeqIdx = 0;
m_NumWorkerInsts = 0;
m_CompletedWorkerInsts = 0;
m_ExpNumWorkerInsts = NumThreads;
m_ReqTerminate = 0;
pThreadPar = m_WorkerInstances;
for (ThreadIdx = 1; ThreadIdx <= NumThreads; ThreadIdx++, pThreadPar++)
	{
	memset(pThreadPar,0,sizeof(tsCHWorkerInstance));
	pThreadPar->ThreadIdx = ThreadIdx;
	pThreadPar->pThis = this;
#ifdef _WIN32
	pThreadPar->threadHandle = (HANDLE)_beginthreadex(nullptr, cWorkThreadStackSize, WorkerPBAInstance, pThreadPar, 0, &pThreadPar->threadID);
#else
	pThreadPar->threadRslt = pthread_create(&pThreadPar->threadID, &threadattr, WorkerPBAInstance, pThreadPar);
#endif
	}

// allow threads time to all startup
Rslt = WaitWorkerThreadStatus(cMaxWaitThreadsStartup);	// waiting untill all threads have started (some may have completed!)
if(Rslt == 2)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "StartWorkerThreads: All %d expected threads did not register as started within allowed %d seconds", NumThreads, cMaxWaitThreadsStartup);
	Rslt=TerminateWorkerThreads();
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "StartWorkerThreads: TerminateWorkThreads(60) returned %d", Rslt);
	return(eBSFerrInternal);
	}
 return(Rslt);
}

// stop all threads in worker pool
int	// returns number of threads forced to terminate
CCallHaplotypes::TerminateWorkerThreads(int WaitSecs)				// allow at most this many seconds before force terminating worker pool threads
{
int NumForceTerminated;
uint32_t Idx;
uint32_t NumWorkerInsts;
uint32_t ExpNumWorkerInsts;
uint32_t CompletedWorkerInsts;
tsCHWorkerInstance *pThreadPar;
time_t Then;
time_t Now;

AcquireSerialise();
NumWorkerInsts = m_NumWorkerInsts;
ExpNumWorkerInsts = m_ExpNumWorkerInsts;
CompletedWorkerInsts = m_CompletedWorkerInsts;
ReleaseSerialise();
if(ExpNumWorkerInsts == 0)
	{
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "TerminateWorkerThreads: No worker threads to terminate");
	pThreadPar = m_WorkerInstances;
	AcquireSerialise();
	for (Idx = 0; Idx < ExpNumWorkerInsts; Idx++, pThreadPar += 1)
		memset(pThreadPar, 0, sizeof(tsCHWorkerInstance));
	m_NumWorkerInsts = 0;
	m_ExpNumWorkerInsts = 0;
	m_CompletedWorkerInsts = 0;
	m_ReqTerminate = 0;
	ReleaseSerialise();
	return(0);
	}

// request all worker threads to self terminate
AcquireSerialise();
m_ReqTerminate = 1;
ReleaseSerialise();
if (CompletedWorkerInsts < ExpNumWorkerInsts)
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "TerminateWorkerThreads: Requesting %u worker threads to terminate", ExpNumWorkerInsts - CompletedWorkerInsts);
Then = time(nullptr) + WaitSecs;
NumForceTerminated = 0;
pThreadPar = m_WorkerInstances;
for(Idx = 0; Idx < ExpNumWorkerInsts; Idx++, pThreadPar += 1)
	{
	Now = time(nullptr);
	if(Now >= Then)
		Now = 1;
	else
		Now = Then - Now;

#ifdef WIN32
	if(pThreadPar->threadHandle != nullptr)
		{
		if(WAIT_TIMEOUT == WaitForSingleObject(pThreadPar->threadHandle, (uint32_t)Now * 1000))
			{
			NumForceTerminated += 1;
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "TerminateWorkerThreads: Force terminating thread %u, WaitForSingleObject() returned WAIT_TIMEOUT", pThreadPar->ThreadIdx);
			TerminateThread(pThreadPar->threadHandle,0);
			}
		pThreadPar->threadHandle = nullptr;
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
			pthread_join(pThreadPar->threadID, nullptr);
			}
		pThreadPar->threadID = 0;
		}
#endif
	}

pThreadPar = m_WorkerInstances;
AcquireSerialise();
for(Idx = 0; Idx < ExpNumWorkerInsts; Idx++, pThreadPar += 1)
	memset(pThreadPar,0,sizeof(tsCHWorkerInstance));
m_NumWorkerInsts = 0;
m_ExpNumWorkerInsts = 0;
m_CompletedWorkerInsts = 0;
m_ReqTerminate = 0;
ReleaseSerialise();
return(NumForceTerminated);
}

// stop all threads in worker load chromosome PBAs thread pool
int
CCallHaplotypes::TerminateLoadChromPBAsThreads(int WaitSecs)				// allow at most this many seconds before force terminating worker pool threads
{
int NumForceTerminated;
uint32_t Idx;
uint32_t NumWorkerInsts;
uint32_t ExpNumWorkerInsts;
uint32_t CompletedWorkerInsts;
tsCHWorkerLoadChromPBAsInstance* pThreadPar;
time_t Then;
time_t Now;

AcquireSerialise();
NumWorkerInsts = m_NumWorkerInsts;
ExpNumWorkerInsts = m_ExpNumWorkerInsts;
CompletedWorkerInsts = m_CompletedWorkerInsts;
ReleaseSerialise();
if (ExpNumWorkerInsts == 0)
	{
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "TerminateLoadChromPBAsThreads: No worker threads to terminate");
	pThreadPar = m_WorkerLoadChromPBAInstances;
	AcquireSerialise();
	for (Idx = 0; Idx < ExpNumWorkerInsts; Idx++, pThreadPar += 1)
		memset(pThreadPar, 0, sizeof(tsCHWorkerLoadChromPBAsInstance));
	m_NumWorkerInsts = 0;
	m_ExpNumWorkerInsts = 0;
	m_CompletedWorkerInsts = 0;
	m_ReqTerminate = 0;
	ReleaseSerialise();
	return(0);
	}

// request all worker threads to self terminate
AcquireSerialise();
m_ReqTerminate = 1;
ReleaseSerialise();
if(NumWorkerInsts != CompletedWorkerInsts)
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "TerminateLoadChromPBAsThreads: Requesting %u worker threads to terminate", NumWorkerInsts);
Then = time(nullptr) + WaitSecs;
NumForceTerminated = 0;
pThreadPar = m_WorkerLoadChromPBAInstances;
for (Idx = 0; Idx < ExpNumWorkerInsts; Idx++, pThreadPar += 1)
	{
	Now = time(nullptr);
	if (Now >= Then)
		Now = 1;
	else
		Now = Then - Now;

#ifdef WIN32
	if (pThreadPar->threadHandle != nullptr)
		{
		if (WAIT_TIMEOUT == WaitForSingleObject(pThreadPar->threadHandle, (uint32_t)Now * 1000))
			{
			NumForceTerminated += 1;
			TerminateThread(pThreadPar->threadHandle, 0);
			}
		pThreadPar->threadHandle = nullptr;
		}
#else
	if (pThreadPar->threadID != 0)
		{
		struct timespec ts;
		int JoinRlt;
		void* pExitRslt;
		clock_gettime(CLOCK_REALTIME, &ts);
		ts.tv_sec += Now;
		if ((JoinRlt = pthread_timedjoin_np(pThreadPar->threadID, &pExitRslt, &ts)) != 0)
			{
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "TerminateLoadChromPBAsThreads: Force terminating thread %u, pthread_timedjoin_np() returned %d", pThreadPar->ThreadIdx, JoinRlt);
			NumForceTerminated += 1;
			pthread_cancel(pThreadPar->threadID);
			pthread_join(pThreadPar->threadID, nullptr);
			}
		pThreadPar->threadID = 0;
		}
#endif
	}

pThreadPar = m_WorkerLoadChromPBAInstances;
AcquireSerialise();
for (Idx = 0; Idx < ExpNumWorkerInsts; Idx++, pThreadPar += 1)
	memset(pThreadPar, 0, sizeof(tsCHWorkerLoadChromPBAsInstance));
m_NumWorkerInsts = 0;
m_ExpNumWorkerInsts = 0;
m_CompletedWorkerInsts = 0;
m_ReqTerminate = 0;
ReleaseSerialise();
return(NumForceTerminated);
}

// Using threads to concurrently load chromosome PBAs
// ChromID is specified and range of sample identifiers
int
CCallHaplotypes::ProcWorkerLoadChromPBAThread(tsCHWorkerLoadChromPBAsInstance* pThreadPar)	// worker thread parameters
{
int Rslt;
uint8_t PBA;
uint8_t* pPBAs;
int32_t Ofs;
int32_t CurSampleID;
uint32_t ReqTerminate;
// one more thread instance has started
AcquireSerialise();
m_NumWorkerInsts++;
ReleaseSerialise();

for(CurSampleID = pThreadPar->StartSampleID; CurSampleID <= pThreadPar->EndSampleID; CurSampleID++)
	{
	// check if requested to terminate
	AcquireSerialise();
	ReqTerminate = m_ReqTerminate;
	ReleaseSerialise();
	if(ReqTerminate)
		break;

	// it is possible for a chromosome to be missing in a PBA if no reads were aligned to that chromosome
	// if so then continue loading from remaining samples
	if((pPBAs = LoadSampleChromPBAs(CurSampleID, pThreadPar->ChromID)) == nullptr)
		continue;

	if(pThreadPar->bNormAlleles)
		{
		tsCHChromMetadata *pChromMetadata = LocateChromMetadataFor(CurSampleID, pThreadPar->ChromID);
		for (Ofs = 0; Ofs < pChromMetadata->ChromLen; Ofs++, pPBAs++)
			{
			if ((PBA = *pPBAs) == 0)
				continue;

			if (PBA & 0xc0)
				PBA |= 0xc0;
			if (PBA & 0x30)
				PBA |= 0x30;
			if (PBA & 0x0c)
				PBA |= 0x0c;
			if (PBA & 0x03)
				PBA |= 0x03;
			*pPBAs = PBA;
		}
		}
	}
if(CurSampleID <= pThreadPar->EndSampleID)
	Rslt = -1;
else
	Rslt = eBSFSuccess;
AcquireSerialise();
m_CompletedWorkerInsts++;
ReleaseSerialise();
return(Rslt);
}

int
CCallHaplotypes::ProcWorkerThread(tsCHWorkerInstance *pThreadPar)	// worker thread parameters
{
int Rslt;
int CurBinID;
uint32_t NumQueueElsProcessed;
tsCHWorkQueueEl *pWorkQueueEl;
uint32_t ReqTerminate;
// this thread has started, one more worker thread
AcquireSerialise();
m_NumWorkerInsts++;
ReleaseSerialise();
NumQueueElsProcessed = 0;
pWorkQueueEl = nullptr;
CurBinID = m_ProccessingHGBinID;
Rslt = eBSFSuccess;
while(Rslt >= eBSFSuccess) {
	// check if requested to terminate
	AcquireSerialise();
	ReqTerminate = m_ReqTerminate;
	ReleaseSerialise();
	if(ReqTerminate)
		{
		Rslt = -1;
		break;
		}
	if(m_PMode < eMCSHAllelicHapsGrps) 
		{
		AcquireSerialise();
    	NumQueueElsProcessed = 	m_NumQueueElsProcessed; // multiple work queue elements to be processed, each thread is processing one of these work elements
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
	else // a single work queue element per chromosome - and just one chromosome to be processed, so just one work queue element! threads iterate over haplotype grouping bins for that chromosome
		{
		tsHGBinSpec* pCurBinSpec;
		pWorkQueueEl = &m_pWorkQueueEls[0]; // remember, just a single work queue element ...
		if((pCurBinSpec = IterateHGBinSpecs(CurBinID,pWorkQueueEl->ChromID)) == nullptr)
			break;
		CurBinID = pCurBinSpec->BinID;
		Rslt = GenHaplotypeGroups(pCurBinSpec,pWorkQueueEl->NumFndrs, pWorkQueueEl->pFounderPBAs, pWorkQueueEl->ChromLen);
		AcquireFastSerialise();
		pCurBinSpec->ProcState |= 0x07; // bin processing state  - 0x00 if unprocessed, 0x01 if currently being processing, 0x03 if processing completed but referenced chromosome not located, 0x07 if processing successfully completed
		ReleaseFastSerialise();
		}
	}

AcquireSerialise();
m_CompletedWorkerInsts++;
ReleaseSerialise();
return(Rslt);
}

int				// returns 0 if all threads started and completed processing, 1 if all threads started but some yet to complete within WaitSecs, 2 if not all threads have started within WaitSecs
CCallHaplotypes::WaitWorkerThreadStatus(int WaitSecs)	// // wait at most this many seconds for all threads to start and complete processing before reporting on worker thread status
{
int64_t WaitMS;
if(WaitSecs < 1)
	WaitSecs = 1;
WaitMS = (int64_t)WaitSecs * 1000;

uint32_t CompletedInstances;
uint32_t NumWorkerInsts;
uint32_t ExpNumWorkerInsts;
do {
	CUtility::SleepMillisecs(500);		// poll every 0.5 secs, restricting poll frequency
	AcquireSerialise();
	CompletedInstances = m_CompletedWorkerInsts;
	NumWorkerInsts = m_NumWorkerInsts;
	ExpNumWorkerInsts = m_ExpNumWorkerInsts;
	ReleaseSerialise();
	WaitMS -= 500;
	}
while(ExpNumWorkerInsts > 0 &&		// surely must be some workers expected!
	 CompletedInstances < ExpNumWorkerInsts	 && // have all workers completed 
	 WaitMS > 0);					// within wait period (WaitSecs)

return((NumWorkerInsts == ExpNumWorkerInsts && CompletedInstances == ExpNumWorkerInsts) ? 0 : (NumWorkerInsts == ExpNumWorkerInsts ? 1 : 2));
}


int				// < 0 if errors, otherwise success
CCallHaplotypes::GenChromAlleleStacks(int32_t ChromID,	// processing is for this chromosome
										  int32_t ChromSize,		// chromosome is this size
										  int32_t Loci,			// processing for allele stacks starting from this loci
										  int32_t MaxNumLoci,		// processing for this maximum number of loci
										  int32_t NumFndrs,		// number of founders to be processed 1st PBA at *pPBAs[0]
										  uint8_t* pFounderPBAs[],	// pPBAs[] pts to chromosome PBAs, for each of the chromosome founder PBAs
										  uint8_t* pMskPBA)			// pts to optional masking PBA, scoring only for segments contained in mask which are non-zero and where progeny and founder PBAs also have an allele.
																	// enables a founder to be processed as if a progeny, restricting scoring KMers to be same as if a progeny

{
int32_t FounderIdx;
tsAlleleStack AlleleStack;
uint32_t NumAlleleStacks;
bool bFndrAllele;
int32_t PrevAlleleIdx;
int32_t NumFndrs2Proc;
uint8_t* pFndrLoci;
int32_t AlleleIdx;
uint8_t AlleleMsk;
int32_t AlleleLoci;
int32_t FndrIdx;
uint32_t DivAlleles;
uint32_t LociFounders;
uint8_t FndrBaseAllele;
tsBitsVect Fndrs2Proc;

if(NumFndrs < 1 || NumFndrs > cMaxFounderReadsets)
	return(-1);

// check that at least 1 founder actually has a PBAs
memset(&Fndrs2Proc, 0, sizeof(Fndrs2Proc));
NumFndrs2Proc = 0;
for(FounderIdx = 0; FounderIdx < NumFndrs; FounderIdx++)
	{
	if(!(m_Fndrs2Proc[FounderIdx] & 0x01))		// not interested as founder not marked for scoring?
		continue;
	if((pFndrLoci = pFounderPBAs[FounderIdx]) == nullptr)	// can't score if no PBA for this founder
		continue;
	BitsVectSet(NumFndrs2Proc, Fndrs2Proc);
	NumFndrs2Proc++;							// can use this founder for generating an allelestack
	}
if(NumFndrs2Proc < 1)							// must have at least 1 founder
	return(-2);

int32_t EndLoci = min(Loci + MaxNumLoci, ChromSize);

memset(&AlleleStack, 0, sizeof(AlleleStack));
AlleleStack.ChromID = ChromID;
AlleleStack.NumFndrs = NumFndrs;
AlleleStack.NumProcFndrs = NumFndrs2Proc;
memcpy(&AlleleStack.ProcFndrs, &Fndrs2Proc, sizeof(tsBitsVect));

bFndrAllele = false;
NumAlleleStacks = 0;
for(AlleleLoci = Loci; AlleleLoci < EndLoci; AlleleLoci++)
	{
	if(pMskPBA != nullptr && pMskPBA[AlleleLoci] == 0)	// skipping over any loci in control which has no PBA
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
		if(pFounderPBAs[FndrIdx] == nullptr)
			{
			continue;
			}

		pFndrLoci = pFounderPBAs[FndrIdx] + AlleleLoci;
		if(*pFndrLoci == 0)					
			{
			if(m_bAllFndrsLociAligned)           // all founders must have alignment alleles at the AlleleLoci
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

	if(DivAlleles < 1) 				// must be at least one founder unique relative to others to be informative...
		continue;

	AddAlleleStack(&AlleleStack);
	NumAlleleStacks++;
	}
return(NumAlleleStacks);
}


uint32_t				// the total number of non-consensus alleles over all founders within KMer NumLoci
CCallHaplotypes::GenConsensusPBA(int32_t Loci,			// processing for consensus starting from this loci
			int32_t NumLoci,		// processing for generation of this maximum sized consensus PBAs
			int32_t NumFndrs,		// number of founders to be processed 1st PBA at *pPBAs[0]
			uint8_t* pFounderPBAs[], // pPBAs[] pts to chromosome PBAs, for each of the chromosome founder PBAs
			uint8_t* pConsensusPBAs)     // ptr to preallocated sequence of length NumLoci which is to be updated with consensus PBAs

{
uint32_t AlleleFreq[256];
uint8_t *pConsensusPLoci;
uint8_t* pFndrLoci;
int32_t AlleleLoci;
int32_t FndrIdx;
int32_t EndLoci;
uint32_t NumNonConsensus;
uint8_t MaxFreqAllele;

pConsensusPLoci = pConsensusPBAs;
EndLoci = Loci + NumLoci;
for(AlleleLoci = Loci; AlleleLoci < EndLoci; AlleleLoci++, pConsensusPLoci++)
	{
	memset(AlleleFreq,0,sizeof(AlleleFreq));
	MaxFreqAllele = 0;
	for(FndrIdx = 0; FndrIdx < NumFndrs; FndrIdx++)
		{
		pFndrLoci = pFounderPBAs[FndrIdx] + AlleleLoci;
		AlleleFreq[*pFndrLoci]++;
		if(AlleleFreq[*pFndrLoci] > AlleleFreq[MaxFreqAllele])
			MaxFreqAllele = *pFndrLoci;
		}
	*pConsensusPLoci = MaxFreqAllele;
	}
NumNonConsensus = 0;
for (FndrIdx = 0; FndrIdx < NumFndrs; FndrIdx++)
	{
	pConsensusPLoci = pConsensusPBAs;
	pFndrLoci = pFounderPBAs[FndrIdx] + Loci;
	for(AlleleLoci = 0; AlleleLoci < NumLoci; AlleleLoci++, pFndrLoci++, pConsensusPLoci++)
		if (*pFndrLoci != *pConsensusPLoci)
			NumNonConsensus++;
	}
return(NumNonConsensus);
}

uint8_t				// consensus PBA for founders in a group of which specified founder is a member at the requested loci
CCallHaplotypes::GenFounderConsensusPBA(int32_t Founder, // consensus for all founders in same group as this founder
								tsHaplotypeGroup *pHaplotypes, // pts to haplotype grouping expected to contain requested chrom and loci
								int32_t ChromID,        // requiring consensus for this chrom at Loci
								int32_t ChromSize,		 // chrom size
								int32_t Loci,			// requiring consensus at this loci
								int32_t NumFndrs,		// number of founders to be processed 1st PBA at *pPBAs[0]
								uint8_t* pFounderPBAs[]) // pPBAs[] pts to chromosome PBAs, for each of the chromosome founder PBAs
{
int32_t HapGrp;
tsBitsVect* pHapGrp;
uint32_t AlleleFreq[256];
uint8_t* pFndrLoci;
int32_t FndrIdx;
uint8_t MaxFeqAllele;
uint32_t Max2ndFreqAllele;
int32_t NumGrpMembers;

if (Loci < 0 || ChromSize < 1 || Loci > ChromSize - 1) // sanity check!
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenFounderConsensusPBA: Loci : %d outside of chromosome size: %d for ChromID: %d", Loci, ChromSize, ChromID);
	return(0);
	}

if(pHaplotypes == nullptr)
    return(0);

if(ChromID != pHaplotypes->ChromID || (Loci < pHaplotypes->StartLoci || (Loci >= pHaplotypes->StartLoci + pHaplotypes->NumLoci)))
  return(0);

if(pHaplotypes->MRAGrpChromID != ChromID || pHaplotypes->MRAGrpLoci != Loci) // different loci from that cached? If so then need to clear cache so consensus entries are recalculated
	{
	pHaplotypes->MRAGrpChromID = 0;
	pHaplotypes->MRAGrpLoci = -1;
	memset(pHaplotypes->MRAGrpConsensus,0x0ff,sizeof(pHaplotypes->MRAGrpConsensus));
	}

// identify in which group the founder is a member
pHapGrp = pHaplotypes->HaplotypeGroup;
for(HapGrp = 0; HapGrp < pHaplotypes->ActualHaplotypeGroups; HapGrp++,pHapGrp++)
	{
	if(BitsVectTest(Founder, *pHapGrp))
		break;
	}
if(HapGrp == pHaplotypes->ActualHaplotypeGroups) // should have been a member of a group, but check!
	return(0);

if(pHaplotypes->MRAGrpChromID == ChromID && pHaplotypes->MRAGrpLoci == Loci && pHaplotypes->MRAGrpConsensus[HapGrp] != 0x0ff) // with any luck the consensus for group has already been generated and cached
    return(pHaplotypes->MRAGrpConsensus[HapGrp]); // returning the cached consensus

// group has been identified, generate and cache the consensus, consensus will be highest frequency allele within that groups members
pHaplotypes->MRAGrpChromID = ChromID;
pHaplotypes->MRAGrpLoci = Loci;
memset(AlleleFreq,0,sizeof(AlleleFreq));
MaxFeqAllele = 0;
Max2ndFreqAllele = 0;
NumGrpMembers = 0;
for(FndrIdx = 0; FndrIdx < NumFndrs; FndrIdx++)
	{
	if(!BitsVectTest(FndrIdx, *pHapGrp))
		continue;
	NumGrpMembers += 1;
	pFndrLoci = pFounderPBAs[FndrIdx] + Loci;
	AlleleFreq[*pFndrLoci]++;
	if(AlleleFreq[*pFndrLoci] > AlleleFreq[MaxFeqAllele])
		MaxFeqAllele = *pFndrLoci;
	else
		if (AlleleFreq[*pFndrLoci] > AlleleFreq[Max2ndFreqAllele])
			Max2ndFreqAllele = *pFndrLoci;
	}
// A major issue with consensus alleles arises when alignments are sparse - RNA-seq around exon splice sites plus where there are multiple isoforms, GBS, Skims, InDels are some common examples
// These sparse alignments are evidenced in the haplotype groups as some members having no alleles whilst others do..
// If the proportion of none-allele members to allele members is less than a threshold then it is assumed that the true consensus haplotype allele is the next most frequent allele
if (m_SparseRepProp > 0.0 && MaxFeqAllele == 0x00 && NumGrpMembers >= m_SparseRepPropGrpMbrs && (((double)AlleleFreq[Max2ndFreqAllele] / NumGrpMembers) >= m_SparseRepProp))
	MaxFeqAllele = Max2ndFreqAllele;

pHaplotypes->MRAGrpConsensus[HapGrp] = MaxFeqAllele; // cache group consensus
return(MaxFeqAllele);
}




uint8_t				// concensus PBA for founders in a group of which specified founder is a member at the requested loci
CCallHaplotypes::GenGroupConsensusPBA(int32_t GrpIdx, // consensus for all founders in this group (0 based)
							tsHaplotypeGroup *pHaplotypes, // pts to haplotype grouping expected to contain requested chrom and loci
							int32_t ChromID,        // requiring consensus for this chrom at Loci
							int32_t Loci,			// requiring consensus at this loci
							int32_t NumFndrs,		// number of founders to be processed 1st PBA at *pPBAs[0]
							uint8_t* pFounderPBAs[]) // pPBAs[] pts to chromosome PBAs, for each of the chromosome founder PBAs
{
tsBitsVect* pHapGrp;
int32_t NumGrpMembers;
uint32_t AlleleFreq[256];
uint8_t* pFndrLoci;
int32_t FndrIdx;
uint32_t MaxFreq;
uint8_t MaxFeqAllele;

if(pHaplotypes == nullptr)
    return(0);

if(ChromID != pHaplotypes->ChromID || (Loci < pHaplotypes->StartLoci || (Loci >= pHaplotypes->StartLoci + pHaplotypes->NumLoci)))
  return(0);

if(GrpIdx < 0 || GrpIdx >= pHaplotypes->ActualHaplotypeGroups)
	return(0);

pHapGrp = &pHaplotypes->HaplotypeGroup[GrpIdx];

// group has been identified, generate and cache the consensus, consensus will be highest frequency allele within that groups members
memset(AlleleFreq,0,sizeof(AlleleFreq));
MaxFeqAllele = 0;
MaxFreq = 0;
NumGrpMembers = 0;
for(FndrIdx = 0; FndrIdx < NumFndrs; FndrIdx++)
	{
	if(!BitsVectTest(FndrIdx, *pHapGrp))
		continue;
	NumGrpMembers++;
	pFndrLoci = pFounderPBAs[FndrIdx] + Loci;
	AlleleFreq[*pFndrLoci]++;
	if(AlleleFreq[*pFndrLoci] >= AlleleFreq[MaxFeqAllele])
		{
		MaxFeqAllele = *pFndrLoci;
		MaxFreq = AlleleFreq[*pFndrLoci];
		}
	}

if(NumGrpMembers < m_MinDGTGrpMembers)         // groups with fewer than this number of members are treated as noise alleles these groups are not used when determining DGT group specific major alleles
	return(0xfe);

if(((double)NumGrpMembers/m_NumFounders) < m_MinDGTGrpPropTotSamples) // if group contains a low proportion of total samples then treat these as simply being noise
	return(0xfe); // alleles in these groups are excluded when determining if there are major alleles specific to groups with more significant numbers of members

if(((double)MaxFreq / NumGrpMembers) < m_MinDGTFmeasure) // group has a sufficiently high proportion of samples but do a significant majority of group members share the same allele 
	return(0xff); // flags group as having a mixture of alleles at ChromID.Loci so skip calling group specific DGT allele

return(MaxFeqAllele); // max frequency allele could actually be no alleles - a deletion relative to reference - but will still count when determining if there are two or more groups of which at least one has a group specific allele
}



int32_t				// error or success (>=0) code 
CCallHaplotypes::ProcessGrpLociDGTs(int32_t MinDGTGrpMembers,        // groups with fewer than this number of members are treated as noise alleles and these groups are not used when determining DGT group specific major alleles
									double MinDGTGrpPropTotSamples,  // groups containing less than this proportion of all samples are treated as if containing noise and alleles in these groups are not used when determining DGT group specific major alleles 
									double MinDGTFmeasure,           // only accepting DGT loci with at least this F-measure score
									char* pszHapGrpFile,             // input, previously generated by 'callhaplotypes', haplotype group file (CSV format) 
									int32_t NumSamplePBAsInputFiles,	     // number of input sample PBAs file specs
									char* pszSamplePBAsInputFiles[],	 // names of input sample PBAs files (wildcards allowed)
									char* pszOutFile)		         // write haplotype group loci DGTs to this output file (CSV format)
{
int Rslt;
gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessGrpLociDGTs: Initialising for processing of sample PBAs");
Rslt = eBSFSuccess;		// assume success!
CSimpleGlob glob(SG_GLOB_FULLSORT);
int32_t Idx;
char* pszInFile;
int32_t ReadsetID = 0;
int32_t NumFiles;
int32_t TotNumFiles = 0;

if(m_pInBuffer == nullptr)
	{
	if((m_pInBuffer = new uint8_t[cInBuffSize]) == nullptr)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessGrpLociDGTs: Memory allocation of %u bytes for input buffering failed - %s", cInBuffSize, strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
	m_AllocInBuff = cInBuffSize;
	}
m_InNumBuffered = 0;
// load all samples with out actually allocating memory for each individual chromosome PBAs, but the file offset at which the chromosome PBA starts will be recorded  
for(Idx = 0; Idx < NumSamplePBAsInputFiles; Idx++)
	{
	glob.Init();
	if(glob.Add(pszSamplePBAsInputFiles[Idx]) < SG_SUCCESS)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessGrpLociDGTs: Unable to glob '%s' sample PBA files", pszSamplePBAsInputFiles[Idx]);
		Reset();
		return(eBSFerrOpnFile);	// treat as though unable to open file
		}
	if((NumFiles = glob.FileCount()) <= 0)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessGrpLociDGTs: Unable to locate any sample PBA file matching '%s", pszSamplePBAsInputFiles[Idx]);
		Reset();
		return(eBSFerrFileName);
		}

	if(m_LimitPrimaryPBAs > 0 && (TotNumFiles + NumFiles) > m_LimitPrimaryPBAs)
		{
		NumFiles = m_LimitPrimaryPBAs - TotNumFiles;
		if(NumFiles == 0)
			break;
		}
	TotNumFiles += NumFiles;

	if(TotNumFiles > cMaxFounderReadsets)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessGrpLociDGTs: Can accept at most %d sample PBA files for processing, after wildcard file name expansions there are %d requested", cMaxFounderReadsets, TotNumFiles);
		Reset();
		return(eBSFerrMem);
		}

	Rslt = eBSFSuccess;
	for(int32_t FileID = 0; Rslt >= eBSFSuccess && FileID < NumFiles; ++FileID)
		{
		pszInFile = glob.File(FileID);
		ReadsetID = LoadPBAFile(pszInFile, 0, true); // loading without allocation for chromosome PBAs
		if(ReadsetID <= 0)
			{
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessGrpLociDGTs: Errors loading file '%s'", pszInFile);
			Reset();
			return(ReadsetID);
			}
		m_Fndrs2Proc[ReadsetID - 1] = 0x01;	// if loaded then assumption is that this founder will be processed
		}
	if(m_LimitPrimaryPBAs > 0 && m_NumReadsetNames >= m_LimitPrimaryPBAs)
		break;
	}
m_NumFounders = TotNumFiles;

// now load the haplotype groupings which were previously generated by 'callhaplotypes' with PMode eMCSHAllelicHapsGrps
glob.Init();
if(glob.Add(pszHapGrpFile) < SG_SUCCESS)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessGrpLociDGTs: Unable to glob haplotype group clustering file(s) '%s", pszHapGrpFile);
	Reset();
	return(eBSFerrOpnFile);	// treat as though unable to open file
	}
if((NumFiles = glob.FileCount()) <= 0)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessGrpLociDGTs: Unable to locate any haplotype group clustering file(s) matching '%s", pszHapGrpFile);
	Reset();
	return(eBSFerrOpnFile);	// treat as though unable to open file
	}
Rslt = eBSFSuccess;
for(int32_t FileID = 0; Rslt >= eBSFSuccess && FileID < NumFiles; ++FileID)
	{
	pszInFile = glob.File(FileID);
	if(Rslt = LoadHaplotypeGroupings(pszInFile) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessGrpLociDGTs: Failed processing previously generated haplotype group clusters file '%s", pszInFile);
		Reset();
		return(Rslt);
		}
	if(m_UsedHGBinSpecs == 0)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessGrpLociDGTs: Unable to parse out any previously generated haplotype group clusters file from file(s) matching '%s", pszHapGrpFile);
		Reset();
		return(eBSFerrOpnFile);	// treat as though unable to open file
		}

	if(m_UsedHGBinSpecs > 1)
		m_mtqsort.qsort(m_pHGBinSpecs, (int64_t)m_UsedHGBinSpecs, sizeof(tsHGBinSpec), SortHGBinSpecs); // sort ascending by ChromID.StartLoci.NumLoci.Distance ascending
	}

// now have sample PBA readsets loaded - without any allocations for actual chromosome PBAs - and haplotype grouping bins 
// iterate over each sample input file, each time processing a single chromosome
char szOutFile[_MAX_PATH];
int32_t SampleID;
int32_t DelSampleID;
int32_t ChromID;
uint8_t **ppPBAs = nullptr;

if(m_pszOutBuffer == nullptr)
	{
	if((m_pszOutBuffer = new uint8_t[cOutBuffSize]) == nullptr)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessGrpLociDGTs: Memory allocation of %u bytes for output buffering failed - %s", cOutBuffSize, strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
	m_AllocOutBuff = cOutBuffSize;
	}
m_OutBuffIdx = 0;

if(m_hOutFile == -1)
	{
	sprintf(szOutFile, "%s.HapGrpDGTs.%d.%d.csv", pszOutFile, m_ExprID,m_NumFounders);
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessGrpLociDGTs: Reporting haplotype grouping DGTs to file '%s'", szOutFile);
#ifdef _WIN32
	m_hOutFile = open(szOutFile, (O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC), (_S_IREAD | _S_IWRITE));
#else
	if((m_hOutFile = open64(szOutFile, O_WRONLY | O_CREAT, S_IREAD | S_IWRITE)) != -1)
		if(ftruncate(m_hOutFile, 0) != 0)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessGrpLociDGTs: Unable to create/truncate %s - %s", szOutFile, strerror(errno));
			Reset();
			return(eBSFerrCreateFile);
			}
#endif
	if(m_hOutFile < 0)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessGrpLociDGTs: Unable to create/truncate %s - %s", szOutFile, strerror(errno));
		Reset();
		return(eBSFerrCreateFile);
		}
	m_OutBuffIdx = sprintf((char *)m_pszOutBuffer, "\"ExprID\",\"RowID\",\"SrcExprID\",\"SrcRowID\",\"Chrom\",\"Loci\",\"GrpAllele:A\",\"GrpAllele:C\",\"GrpAllele:G\",\"GrpAllele:T\",\"GrpLog2FMeasure:A\",\"GrpLog2FMeasure:C\",\"GrpLog2FMeasure:G\",\"GrpLog2FMeasure:T\"");
	m_OutBuffIdx += sprintf((char *)&m_pszOutBuffer[m_OutBuffIdx], ",\"ActualHaplotypeGroups\",\"Grp1Mbrs\",\"Grp1Alleles:A\",\"Grp1Alleles:C\",\"Grp1Alleles:G\",\"Grp1Alleles:T\",\"Grp1Alleles:N\"");
	m_OutBuffIdx += sprintf((char *)&m_pszOutBuffer[m_OutBuffIdx], ",\"Grp2Mbrs\",\"Grp2Alleles:A\",\"Grp2Alleles:C\",\"Grp2Alleles:G\",\"Grp2Alleles:T\",\"Grp2Alleles:N\"");
	m_OutBuffIdx += sprintf((char *)&m_pszOutBuffer[m_OutBuffIdx], ",\"Grp3Mbrs\",\"Grp3Alleles:A\",\"Grp3Alleles:C\",\"Grp3Alleles:G\",\"Grp3Alleles:T\",\"Grp3Alleles:N\"");
	m_OutBuffIdx += sprintf((char *)&m_pszOutBuffer[m_OutBuffIdx], ",\"Grp4Mbrs\",\"Grp4Alleles:A\",\"Grp4Alleles:C\",\"Grp4Alleles:G\",\"Grp4Alleles:T\",\"Grp4Alleles:N\"");
	m_OutBuffIdx += sprintf((char *)&m_pszOutBuffer[m_OutBuffIdx], ",\"Grp5Mbrs\",\"Grp5Alleles:A\",\"Grp5Alleles:C\",\"Grp5Alleles:G\",\"Grp5Alleles:T\",\"Grp5Alleles:N\"");
	m_OutBuffIdx += sprintf((char *)&m_pszOutBuffer[m_OutBuffIdx], ",\"GrpPseudoMbrs\",\"GrpPseudoAlleles:A\",\"GrpPseudoAlleles:C\",\"GrpPseudoAlleles:G\",\"GrpPseudoAlleles:T\",\"GrpPseudoAlleles:N\"");
	m_CurRowID = m_SeedRowID;
	}
ppPBAs = new uint8_t * [cMaxFounderReadsets + 1];
memset(ppPBAs, 0, sizeof(uint8_t*)* (cMaxFounderReadsets + 1));

for(ChromID = 1; ChromID <= m_NumChromNames; ChromID++)
	{
	uint32_t ChromSize;
	char* pszChrom = LocateChrom(ChromID);

	gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessGrpLociDGTs: Begining to load PBAs for %s chromosome for DGT processing .... ", pszChrom);
	int WorkerThreadStatus = StartWorkerLoadChromPBAThreads(m_NumThreads, 1, m_NumFounders, ChromID);
	if (WorkerThreadStatus < 0)
		{
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessGrpLociDGTs: Errors loading PBAs for %s chromosome for DGT processing .... ", pszChrom);
		delete []ppPBAs;
		return(WorkerThreadStatus);
		}
	while (WorkerThreadStatus > 0)
		{
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessGrpLociDGTs: Continuing to load chromosome PBAs for %s .... ", pszChrom);
		WorkerThreadStatus = WaitWorkerThreadStatus(60);
		}
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessGrpLociDGTs: Completed loading chromosome PBAs for %s for DGT processing", pszChrom);
	TerminateLoadChromPBAsThreads(60);

	for(SampleID = 1; SampleID <= m_NumFounders; SampleID++)
		{
		tsCHChromMetadata *pChromMetadata;
		// returned pointer to chromosome metadata
		if((pChromMetadata = LocateChromMetadataFor(SampleID, ChromID)) == nullptr)
			{
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessGrpLociDGTs: No metadata for chromosome '%s' in at least one founder '%s', skipping this chromosome", pszChrom, LocateReadset(SampleID));
			break;;
			}
		ChromSize = pChromMetadata->ChromLen;
		if((ppPBAs[SampleID-1] = pChromMetadata->pPBAs) == nullptr)
			{
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessGrpLociDGTs: No PBA for chromosome '%s' in at least one founder '%s', skipping this chromosome", pszChrom, LocateReadset(SampleID));
			break;
			}
		}
	
	if(SampleID <= m_NumFounders)
		{
		for(DelSampleID = 1; DelSampleID < SampleID; DelSampleID++)
			DeleteSampleChromPBAs(DelSampleID, ChromID);
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessGrpLociDGTs: Unloaded PBAs for chromosome %s completed, ready for next iteration of chromosome DGT processing", pszChrom);
		continue; // slough this chromosome and try next chromosome
		}

	gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessGrpLociDGTs: Loading all PBAs for chromosome %s completed, now calling DGTs", pszChrom);
	if((Rslt = GenBinDGTs(ChromID,          // requiring DGTs across this chrom
		ChromSize,			  // chromosome is this size
		m_NumFounders,		  // number of founders to be processed 1st PBA at *pPBAs[0]
		ppPBAs)) < 0) // ppPBAs[] pts to chromosome PBAs, for each of the chromosome founder PBAs
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessGrpLociDGTs: Fatal error whilst haplotype grouping DGT reporting on chromosome %s",pszChrom);
		Reset();
		return(Rslt);
		}
	if(m_OutBuffIdx)
		{
		if(!CUtility::RetryWrites(m_hOutFile, m_pszOutBuffer, m_OutBuffIdx))
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessGrpLociDGTs: Fatal error in RetryWrites()");
			Reset();
			return(eBSFerrFileAccess);
			}
		m_OutBuffIdx = 0;
		}
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessGrpLociDGTs: Unloading PBAs for chromosome %s ....", pszChrom);
	for(DelSampleID = 1; DelSampleID <= SampleID; DelSampleID++)
		DeleteSampleChromPBAs(DelSampleID, ChromID);
	}
if(ppPBAs != nullptr)
	delete []ppPBAs;

char *pszChrom;
uint32_t CurChromID;
uint32_t DGTIdx;
uint32_t AlleleIdx;
uint32_t GrpIdx;
uint32_t ReportTopN = m_MaxReportGrpDGTs;
uint32_t TotGrpUniqueAlleles = m_UsedDGTLoci;

if(ReportTopN && m_UsedDGTLoci > ReportTopN)
	{
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessGrpLociDGTs: Sorting the DGTs prior to reporting the top %u by F-measure scores", ReportTopN);
	m_mtqsort.SetMaxThreads(m_NumThreads);
	m_mtqsort.qsort(m_pDGTLoci, (int64_t)m_UsedDGTLoci, sizeof(tsDGTLoci), SortDGTLociFMeasureDesc);
	m_mtqsort.qsort(m_pDGTLoci, (int64_t)ReportTopN, sizeof(tsDGTLoci), SortDGTLoci);
	}
if(ReportTopN == 0)
	ReportTopN = m_UsedDGTLoci;
ReportTopN = min(ReportTopN, m_UsedDGTLoci);
if(m_UsedDGTLoci)
	{
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessGrpLociDGTs: Writing to file '%s' %u DGTs", szOutFile,ReportTopN);
	m_CurRowID = 1;
	tsDGTLoci* pDGTLoci;
	tsGrpCnts* pDGTGrp;
	pDGTLoci = m_pDGTLoci;
	CurChromID = 0;
	for(DGTIdx =0; DGTIdx < ReportTopN; DGTIdx++, pDGTLoci++)
		{
		if(CurChromID != pDGTLoci->ChromID)
			{
			CurChromID = pDGTLoci->ChromID;
			pszChrom = LocateChrom(CurChromID);
			}
		m_OutBuffIdx+=sprintf((char *) &m_pszOutBuffer[m_OutBuffIdx], "\n%d,%d,%u,%u,\"%s\",%d",m_ExprID,m_CurRowID++,pDGTLoci->SrcExprID,pDGTLoci->SrcRowID, pszChrom, pDGTLoci->DGTLoci);
		for(AlleleIdx = 0; AlleleIdx < 4; AlleleIdx++)
			m_OutBuffIdx += sprintf((char *) &m_pszOutBuffer[m_OutBuffIdx],",%u",pDGTLoci->AlleleGrpIDs[AlleleIdx]);
		for(AlleleIdx = 0; AlleleIdx < 4; AlleleIdx++)
			m_OutBuffIdx += sprintf((char *) &m_pszOutBuffer[m_OutBuffIdx], ",%1.6f", pDGTLoci->AlleleFMeasures[AlleleIdx]);
		m_OutBuffIdx += sprintf((char *) &m_pszOutBuffer[m_OutBuffIdx], ",%u", pDGTLoci->NumHapGrps);
		pDGTGrp = pDGTLoci->GrpCnts;
		for(GrpIdx = 0; GrpIdx < 6; GrpIdx++,pDGTGrp++)
			{
			m_OutBuffIdx += sprintf((char *) &m_pszOutBuffer[m_OutBuffIdx], ",%u", pDGTGrp->NumMembers);
			for(AlleleIdx = 0; AlleleIdx < 5; AlleleIdx++)
				m_OutBuffIdx += sprintf((char *) &m_pszOutBuffer[m_OutBuffIdx], ",%u", pDGTGrp->NumAllele[AlleleIdx]);
			}

		if(m_OutBuffIdx+1000 > m_AllocOutBuff)
			{
			if(!CUtility::RetryWrites(m_hOutFile, m_pszOutBuffer, m_OutBuffIdx))
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessGrpLociDGTs: Fatal error in RetryWrites()");
				Reset();
				return(eBSFerrFileAccess);
				}
			m_OutBuffIdx = 0;
			}
		}	
	}

gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessGrpLociDGTs: Reported the top (by F-measure) %u DGT loci from %u originally identified", ReportTopN,m_UsedDGTLoci);

if(m_hOutFile != -1)
	{
	if(m_OutBuffIdx)
		{
		if(!CUtility::RetryWrites(m_hOutFile, m_pszOutBuffer, m_OutBuffIdx))
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessGrpLociDGTs: Fatal error in RetryWrites()");
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
Reset();
gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessGrpLociDGTs: Haplotype grouping DGT reporting completed");
return(0);
}

int32_t				// error or success (>=0) code 
CCallHaplotypes::GenBinDGTs(int32_t ChromID,          // requiring DGTs across this chrom
							int32_t ChromSize,			// chromosome is this size
							int32_t NumFndrs,		  // number of founders to be processed 1st PBA at *pPBAs[0]
							uint8_t* pFounderPBAs[]) // pPBAs[] pts to chromosome PBAs, for each of the chromosome founder PBAs
{
int32_t BinIdx;
int32_t CurLoci;
int32_t EndLoci;
int32_t GrpIdx;
int32_t AlleleIdx;
int32_t AcceptedGrpAlleles[4];
double AcceptedGrpAllelesFMeasure[4];
uint32_t UinqueGrpAlleles;
uint32_t TotGrpUniqueAlleles;
tsHGBinSpec* pHGBinSpec;
tsHaplotypeGroup* pHaplotypeGroup;
int32_t NumNoiseGrps;                  // number groups chracterised as noise groups
int32_t NumGrpMembers[6];              // number of members in first 5 groups plus a pseudo group to hold sum of group counts in any groups after the 5th
int32_t MaxGrpMembers;                 // maximum number of members in first 5 groups
int32_t GrpMembers;                    // number of group members in current group
int32_t MaxMembersGrpIdx;              // index of a group having the maximum number of members - could be multiple groups, this is the index of the first group 
int32_t NumDGTGrps;                    // DGTs are called over this number of groups which is limited to the 1st 5 groupsuint32_t NumNoiseGrps;                  // number of groups likely contain background noise alleles - having less than m_MinDGTGrpMembers or m_MinDGTGrpPropTotSamples
uint8_t SampleAllele;                   // raw sample PBA allele combination
uint8_t *pSampleLoci;
uint32_t GrpsAlleleCnts[5 * 6];
uint32_t AllGrpsAlleleCnts[5];
uint32_t *pGrpAlleleCnts;
uint8_t GrpRepAllele;
double GrpFbeta2 = m_FbetameasureGrps * m_FbetameasureGrps;
double GrpRecall;
double GrpPrecision;
double GrpFMeasure;

char* pszChrom;


if(m_UsedHGBinSpecs == 0)
	{
	gDiagnostics.DiagOut(eDLWarn, gszProcName, "GenBinDGTs: No haplotype group bins to report on!"); // not a failure, but still warn user
	return(eBSFSuccess);
	}
if((pszChrom = LocateChrom(ChromID)) == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenBinDGTs: Unable to locate chromosome for ChromID: '%d'",ChromID);
	return(eBSFerrChrom);
	}


gDiagnostics.DiagOut(eDLInfo, gszProcName, "Processing haplotype grouping DGTs for chromosome '%s'", pszChrom);

TotGrpUniqueAlleles = 0;
pHGBinSpec = m_pHGBinSpecs;
for(BinIdx = 0; BinIdx < m_UsedHGBinSpecs; BinIdx++, pHGBinSpec++)
	{
	if((pHGBinSpec->ProcState & 0x07) != 0x07) // bin processing state  - 0x00 if unprocessed, 0x01 if currently being processing, 0x03 if processing completed but referenced chromosome not located, 0x07 if processing successfully completed
		continue;
	if(pHGBinSpec->ChromID != ChromID)
		continue;
	pHaplotypeGroup = pHGBinSpec->pHaplotypeGroup;
	// only interested in bins having at least 2 groups
	if(pHaplotypeGroup->ActualHaplotypeGroups < 2)
		continue;

	// determine number of samples in each group
	MaxGrpMembers = 0;
	MaxMembersGrpIdx = 0;
	memset(NumGrpMembers, 0, sizeof(NumGrpMembers));
	for(GrpIdx = 0; GrpIdx < pHaplotypeGroup->ActualHaplotypeGroups; GrpIdx++)
		{
		GrpMembers = BitsVectCount(pHaplotypeGroup->HaplotypeGroup[GrpIdx]);
		if(GrpIdx > 5)
			NumGrpMembers[5] += GrpMembers;
		else
			NumGrpMembers[GrpIdx] = GrpMembers;
		if(GrpIdx < 5 && GrpMembers > MaxGrpMembers) // identify which group has the maximum number of members - could be more than 1 group but this does not matter as goup counts are scaled to the maximum members
			{
			MaxGrpMembers = GrpMembers;
			MaxMembersGrpIdx = GrpIdx;
			}
		}

	// determine number of groups to be characterised as being noise, and if less than 2 groups with meaningful allele counts remaining then skip current bin
	if(MaxGrpMembers < m_MinDGTGrpMembers)
		continue;
	NumDGTGrps = min(5,pHaplotypeGroup->ActualHaplotypeGroups); // could be more than 5 groups but only calling DGTs on the first 5 groups
	NumNoiseGrps = 0;
	for(GrpIdx = 0; GrpIdx < NumDGTGrps; GrpIdx++)
		if(NumGrpMembers[GrpIdx] < m_MinDGTGrpMembers || ((double)NumGrpMembers[GrpIdx]/pHaplotypeGroup->NumFndrs) < m_MinDGTGrpPropTotSamples)
			NumNoiseGrps++;
	if(NumDGTGrps - NumNoiseGrps < 2) // has to be at least two groups remaining after noise groups have been counted
		continue;

	// 2 or more non-noise group sample counts now known
	// iterate over each bin loci and accumulate allele counts for each group
	EndLoci = pHGBinSpec->StartLoci + pHGBinSpec->NumLoci;
	for(CurLoci = pHGBinSpec->StartLoci; CurLoci < EndLoci; CurLoci++)
		{
		memset(AllGrpsAlleleCnts,0,sizeof(AllGrpsAlleleCnts));
		memset(GrpsAlleleCnts,0,sizeof(GrpsAlleleCnts));
		pGrpAlleleCnts = GrpsAlleleCnts;
		pHaplotypeGroup = pHGBinSpec->pHaplotypeGroup;
		GrpRepAllele = 0;
		for(GrpIdx = 0; GrpIdx < pHaplotypeGroup->ActualHaplotypeGroups; GrpIdx++)
			{
			pGrpAlleleCnts = &GrpsAlleleCnts[min(GrpIdx,5)*5];
			for(uint32_t SampleIdx = 0; SampleIdx < (uint32_t)pHaplotypeGroup->NumFndrs; SampleIdx++)
				{
				if(!BitsVectTest(SampleIdx, pHaplotypeGroup->HaplotypeGroup[GrpIdx])) // if sample not member of group then onto next sample
					continue;
				// if using sparse coverage mode then need to generate a representative allele for members of each group at CurLoci
				if(m_SparseRepPropGrpMbrs > 0)
					GrpRepAllele = GenFounderConsensusPBA(SampleIdx, pHaplotypeGroup, ChromID, ChromSize, CurLoci, NumFndrs, pFounderPBAs);
				pSampleLoci = pFounderPBAs[SampleIdx] + CurLoci; // alleles for current sample at CurLoci
				SampleAllele = *pSampleLoci;
				// sum number of alleles - A,C,G,T - any coverage contributes a count for that allele, note that a sample may have no coverage or multiple alleles at a given loci
				if (SampleAllele == 0)     // no coverage
					{
					if (GrpRepAllele == 0)
						{
						pGrpAlleleCnts[4]++;
						AllGrpsAlleleCnts[4]++;
						continue;
						}
					else
						SampleAllele = GrpRepAllele;	// use the representative allele
					}
				// some samples may have multiple alleles - hence testing for all alleles
				if(SampleAllele & 0xc0)   // A
					{
					pGrpAlleleCnts[0]++;
					AllGrpsAlleleCnts[0]++;
					}
				if(SampleAllele & 0x30)   // C
					{
					pGrpAlleleCnts[1]++;
					AllGrpsAlleleCnts[1]++;
					}
				if(SampleAllele & 0x0c)   // G
					{
					pGrpAlleleCnts[2]++;
					AllGrpsAlleleCnts[2]++;
					}
				if(SampleAllele & 0x03)   // T
					{
					pGrpAlleleCnts[3]++;
					AllGrpsAlleleCnts[3]++;
					}
				}
			}
		
		// for each loci and group calculate the recall and precision
		// recall is using the original raw counts for each group as: ThisGrpsAlleleCnts/ThisGrpsNumberOfSamples
		// precision is using log2 scaling to partially normalise for differing group sizes: LOG2(1+(MaxGrpSamples/ThisGrpsSamples))*ThisGprsAlleleCnt
		//
		double GrpScales[5];       // scaling for each group
		double GrpScaleMembers[5]; // scaled members group members
		double GrpScaleAlleles[5 * 5]; // scaled group alleles
		double SumGrpScaleAlleles[5];
		double* pGrpScaleAlleles;
		double* pSumGrpScaleAlleles;
		double* pGrpScales;
		double* pGrpScaleMembers;
		int32_t* pNumGrpMembers;

		memset(GrpScales, 0, 5*sizeof(double));
		memset(GrpScaleMembers, 0, 5*sizeof(double));
		memset(GrpScaleAlleles, 0, 25 * sizeof(double));
		// scaling counts for each group - using log2 scaling to reduce differences between group sizes
		pGrpScales = GrpScales;
		pGrpScaleMembers = GrpScaleMembers;
		pNumGrpMembers = NumGrpMembers;
		pGrpAlleleCnts = &GrpsAlleleCnts[25]; // pseudo groups allele counts
		pSumGrpScaleAlleles = SumGrpScaleAlleles;
		for(AlleleIdx = 0; AlleleIdx < 5; AlleleIdx++, pGrpAlleleCnts++, pSumGrpScaleAlleles++)
			*pSumGrpScaleAlleles = *pGrpAlleleCnts;

		pGrpAlleleCnts = GrpsAlleleCnts;
		pGrpScaleAlleles = GrpScaleAlleles;
		pSumGrpScaleAlleles = SumGrpScaleAlleles;
		for(GrpIdx = 0; GrpIdx < NumDGTGrps; GrpIdx++,pGrpScales++,pNumGrpMembers++,pGrpScaleMembers++)
			{
			if(*pNumGrpMembers == 0)
				continue;
			if(*pNumGrpMembers < m_MinDGTGrpMembers || ((double)*pNumGrpMembers/pHaplotypeGroup->NumFndrs) < m_MinDGTGrpPropTotSamples)
				{
				*pGrpScales = 1.0f;
				*pGrpScaleMembers = *pNumGrpMembers;
				}
			else
				{
				*pGrpScales = log2(1.0 + ((double)MaxGrpMembers / *pNumGrpMembers));
				*pGrpScaleMembers = *pGrpScales * *pNumGrpMembers;
				}
			pSumGrpScaleAlleles = SumGrpScaleAlleles;
			for(AlleleIdx = 0; AlleleIdx < 5; AlleleIdx++,pGrpScaleAlleles++,pGrpAlleleCnts++,pSumGrpScaleAlleles++)
				{
				*pGrpScaleAlleles = *pGrpScales * *pGrpAlleleCnts;
				*pSumGrpScaleAlleles += *pGrpScaleAlleles;
				}
			}

		UinqueGrpAlleles=0;
		memset(AcceptedGrpAlleles, 0, sizeof(AcceptedGrpAlleles));
		memset(AcceptedGrpAllelesFMeasure, 0, sizeof(AcceptedGrpAllelesFMeasure));
		pGrpAlleleCnts = GrpsAlleleCnts;
		pNumGrpMembers = NumGrpMembers;
		pGrpScaleAlleles = GrpScaleAlleles;
		for(GrpIdx = 0; GrpIdx < NumDGTGrps; GrpIdx++,pNumGrpMembers++,pGrpAlleleCnts++,pGrpScaleAlleles++)
			{
			if(*pNumGrpMembers < m_MinDGTGrpMembers || ((double)*pNumGrpMembers/pHaplotypeGroup->NumFndrs) < m_MinDGTGrpPropTotSamples)
				{
				pGrpAlleleCnts += 4;
				pGrpScaleAlleles += 4;
				continue;
				}
			pSumGrpScaleAlleles = SumGrpScaleAlleles;
			for(AlleleIdx = 0; AlleleIdx < 4; AlleleIdx++,pGrpAlleleCnts++,pGrpScaleAlleles++,pSumGrpScaleAlleles++)
				{
				if(AllGrpsAlleleCnts[AlleleIdx] == 0 || *pGrpAlleleCnts == 0 || *pSumGrpScaleAlleles < 1.0)
					continue;

				GrpRecall = (double)*pGrpAlleleCnts / *pNumGrpMembers; // proportion of samples in group with this allele 
				GrpPrecision = *pGrpScaleAlleles / *pSumGrpScaleAlleles; // proportion of samples which have allele in group relative to all samples having that allele
				GrpFMeasure = (1 + GrpFbeta2) * GrpPrecision * GrpRecall / ((GrpFbeta2 * GrpPrecision) + GrpRecall);
				if(GrpFMeasure >= m_MinDGTFmeasure)
					{
					UinqueGrpAlleles += 1;
					AcceptedGrpAlleles[AlleleIdx] = GrpIdx + 1;
					AcceptedGrpAllelesFMeasure[AlleleIdx] = GrpFMeasure;
					}
				}
			}
		if(UinqueGrpAlleles == 0)
			continue;

		tsGrpCnts GrpCnts[6];
		tsGrpCnts *pGrpCnts;

		memset(GrpCnts, 0, sizeof(GrpCnts));
		pGrpAlleleCnts = GrpsAlleleCnts;
		pGrpCnts = GrpCnts;
		for(GrpIdx = 0; GrpIdx < NumDGTGrps; GrpIdx++, pGrpCnts++)
			{
			pGrpCnts->NumMembers = NumGrpMembers[GrpIdx];
			for(AlleleIdx = 0; AlleleIdx < 5; AlleleIdx++,pGrpAlleleCnts++)
				pGrpCnts->NumAllele[AlleleIdx] = *pGrpAlleleCnts;
			}
		if(pHaplotypeGroup->ActualHaplotypeGroups > 5)
			{
			pGrpAlleleCnts = &GrpsAlleleCnts[25]; // pseudo group alleles
			pGrpCnts->NumMembers = NumGrpMembers[5];
			for(AlleleIdx = 0; AlleleIdx < 5; AlleleIdx++,pGrpAlleleCnts++)
				pGrpCnts->NumAllele[AlleleIdx] = *pGrpAlleleCnts;
			}

		if(UinqueGrpAlleles > 0)
			AddDGTLoci(pHaplotypeGroup->SrcExprID, pHaplotypeGroup->SrcRowID, ChromID, CurLoci,pHaplotypeGroup->ActualHaplotypeGroups, AcceptedGrpAlleles, AcceptedGrpAllelesFMeasure, GrpCnts);
		TotGrpUniqueAlleles++;
		}
	}
gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenBinDGTs: There were %d unique loci on '%s' at which DGT alleles were called",TotGrpUniqueAlleles,pszChrom);
if(m_hOutFile != -1)
	{
	if(m_OutBuffIdx)
		{
		if(!CUtility::RetryWrites(m_hOutFile, m_pszOutBuffer, m_OutBuffIdx))
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenBinDGTs: Fatal error in RetryWrites()");
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
	}
return(TotGrpUniqueAlleles);
}

int32_t 
CCallHaplotypes::ReportHaplotypeGroups(int32_t NumHaplotypes)         // number of haplotypes which have been grouped
{
int32_t BinIdx;
int32_t MaxNumHapGrps;
tsHGBinSpec* pHGBinSpec;
tsHaplotypeGroup* pHaplotypeGroup;
int32_t CurChromID;
int32_t HaplotypeIdx;
int32_t GrpIdx;
char* pszChrom;
char szOutFile[_MAX_PATH];

uint32_t *pNumGrpMembers = nullptr;

if(m_UsedHGBinSpecs == 0)
	{
	gDiagnostics.DiagOut(eDLWarn, gszProcName, "ReportHaplotypeGroups:No haplotype group bins to report on!"); // not a failure, but still warn user
	return(eBSFSuccess);
	}

// determine the maximum number of actual haplotype groups over all bins
MaxNumHapGrps = 0;
pHGBinSpec = m_pHGBinSpecs;
for(BinIdx = 0; BinIdx < m_UsedHGBinSpecs; BinIdx++, pHGBinSpec++)
	{
	if((pHGBinSpec->ProcState & 0x07) != 0x07) // skip any with error status or incomplete processing
		continue;
	if(pHGBinSpec->pHaplotypeGroup->ActualHaplotypeGroups > MaxNumHapGrps)
		MaxNumHapGrps = pHGBinSpec->pHaplotypeGroup->ActualHaplotypeGroups;
	if (MaxNumHapGrps > pHGBinSpec->pHaplotypeGroup->MaxHaplotypeGroups)
		gDiagnostics.DiagOut(eDLWarn, gszProcName,"ReportHaplotypeGroupsMaxNumHapGrps %d more than expected maximum %d", MaxNumHapGrps, pHGBinSpec->pHaplotypeGroup->MaxHaplotypeGroups);
	}

if(m_pszOutBuffer == nullptr)
	{
	if((m_pszOutBuffer = new uint8_t[cOutBuffSize]) == nullptr)
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
m_OutBuffIdx = sprintf((char*)m_pszOutBuffer, "\"ExprID\",\"RowID\",\"Chrom\",\"Loci\",\"Len\",\"MinDistance\",\"MaxDistance\",\"MaxHaplotypeGroups\",\"ActualDistance\",\"ActualHaplotypeGroups\"");
for(HaplotypeIdx = 0; HaplotypeIdx < NumHaplotypes; HaplotypeIdx++)
	m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], ",\"%s\"", LocateReadset(HaplotypeIdx+1));
for(HaplotypeIdx = 0; HaplotypeIdx < MaxNumHapGrps; HaplotypeIdx++)
	m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], ",\"GrpMembers:%d\"", HaplotypeIdx+1);
m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], "\n");

pNumGrpMembers = new uint32_t[MaxNumHapGrps];
m_CurRowID = m_SeedRowID;
CurChromID = 0;
pHGBinSpec = m_pHGBinSpecs;
for(BinIdx = 0; BinIdx < m_UsedHGBinSpecs; BinIdx++,pHGBinSpec++)
	{
	if((pHGBinSpec->ProcState & 0x07) != 0x07) // bin processing state  - 0x00 if unprocessed, 0x01 if currently being processing, 0x03 if processing completed but referenced chromosome not located, 0x07 if processing successfully completed
		continue;

	if(pHGBinSpec->ChromID != CurChromID)
		{
		CurChromID = pHGBinSpec->ChromID;
		pszChrom = LocateChrom(CurChromID);
		}

	pHaplotypeGroup = pHGBinSpec->pHaplotypeGroup;
	m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], "%u,%u,\"%s\",%u,%u,%u,%u,%u,%u,%u", m_ExprID,m_CurRowID++,pszChrom, pHGBinSpec->StartLoci, pHGBinSpec->NumLoci, pHGBinSpec->MinCentroidDistance,pHGBinSpec->MaxCentroidDistance,pHGBinSpec->MaxNumHaplotypeGroups, pHaplotypeGroup->CentroidDistance,pHaplotypeGroup->ActualHaplotypeGroups);
	memset(pNumGrpMembers, 0, sizeof(int32_t) * MaxNumHapGrps);
	for(HaplotypeIdx = 0; HaplotypeIdx < NumHaplotypes; HaplotypeIdx++)
		{
		for(GrpIdx = 0; GrpIdx < pHaplotypeGroup->ActualHaplotypeGroups; GrpIdx++)
			{
			if(BitsVectTest(HaplotypeIdx, pHaplotypeGroup->HaplotypeGroup[GrpIdx]))
				{
				pNumGrpMembers[GrpIdx] += 1;
				m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], ",%d", GrpIdx + 1);
				break;
				}
			}
		}
	for(GrpIdx = 0; GrpIdx < MaxNumHapGrps; GrpIdx++)
		m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], ",%d", pNumGrpMembers[GrpIdx]);

	m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], "\n");
	if((m_OutBuffIdx + 1000) > m_AllocOutBuff)
		{
		if(!CUtility::RetryWrites(m_hOutFile, m_pszOutBuffer, m_OutBuffIdx))
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "ReportHaplotypeGroups: Fatal error in RetryWrites()");
			delete []pNumGrpMembers;
			Reset();
			return(eBSFerrFileAccess);
			}
		m_OutBuffIdx = 0;
		}
	}
delete []pNumGrpMembers;
if(m_hOutFile != -1)
	{
	if(m_OutBuffIdx)
		{
		if(!CUtility::RetryWrites(m_hOutFile, m_pszOutBuffer, m_OutBuffIdx))
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "ReportHaplotypeGroups: Fatal error in RetryWrites()");
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
gDiagnostics.DiagOut(eDLInfo, gszProcName, "Completed reporting haplotype groupings to file '%s'", szOutFile);
return(eBSFSuccess);
}

int32_t				// < 0 error, otherwise the total number of non-consensus PBA over all pPBAs: if 0 then all PBAs of ConsensusLen were identical to each other
CCallHaplotypes::GenGrpConsensus(int32_t ConsensusLen,		// PBAs consensus length
			int32_t NumPBAs,					// number of PBAs to be processed, 1st PBA at *pPBAs[0]
			uint8_t* pPBAs[],					// pPBAs[] pts to PBAs over which consensus is be generated
			uint8_t* pConsensusPBAs)			// preallocated (min ConsensusLen) by caller to hold consensus for pPBAs[]
{
int32_t CurOfs;
int32_t PBAsIdx;
uint8_t *pPBA;
uint32_t AlleleFreq[256];
uint8_t MaxFreqAllele;
int32_t NumNonConsensus;

if(ConsensusLen < 1 || NumPBAs < 1 || pConsensusPBAs == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenGrpConsensus: internal parameter error");
	return(eBSFerrParams);
	}

NumNonConsensus = 0;
for (CurOfs = 0; CurOfs < ConsensusLen; CurOfs++, pConsensusPBAs++)
	{
	memset(AlleleFreq, 0, sizeof(AlleleFreq));
	MaxFreqAllele = 0;
	for (PBAsIdx = 0; PBAsIdx < NumPBAs; PBAsIdx++)
		{
		if ((pPBA = pPBAs[PBAsIdx]) == nullptr)
			continue;
		pPBA += CurOfs;
		AlleleFreq[*pPBA]++;
		if (AlleleFreq[*pPBA] > AlleleFreq[MaxFreqAllele])
			MaxFreqAllele = *pPBA;
		}
	*pConsensusPBAs = MaxFreqAllele;
	for (PBAsIdx = 0; PBAsIdx < NumPBAs; PBAsIdx++)
		{
		if((pPBA = pPBAs[PBAsIdx])==nullptr)
			continue;
		pPBA += CurOfs;
		if (*pPBA != MaxFreqAllele)
			NumNonConsensus++;
		}
	}
return(NumNonConsensus);
}

int32_t				// error or success (>=0) code 
CCallHaplotypes::GenBinKMers(int32_t ChromID,          // requiring bin group segregating KMers over this chrom
	int32_t ChromSize,			// chromosome is this size
	int32_t KMerSize,			// use this KMer sized sequences when identifying group segregation homozygosity scores
	int32_t MinKMerHammings,	// must be at least this hammings between any two group KMer sequences of KMerSize base pair
	int32_t KMerNoneCoverage,				// if > 0 then allow upto this number of bp within a KMer to have no coverage, maximum is 50% of KMerSize
	int32_t NumFndrs,			// number of founders to be processed 1st PBA at *pPBAs[0]
	uint8_t* pFounderPBAs[])	// pPBAs[] pts to chromosome PBAs, for each of the chromosome founder PBAs
{
int32_t BinIdx;
int32_t CurLoci;
int32_t EndLoci;
int32_t GrpIdx;
uint32_t TotGrpUniqueAlleles;
uint32_t TotKMerSizeLociGrpsHammingsAccepted;
uint32_t TotGrpKMerSizeLoci;
tsHGBinSpec* pHGBinSpec;
tsHaplotypeGroup* pHaplotypeGroup;
int32_t NumNoiseGrps;                  // number groups chracterised as noise groups
int32_t NumGrpMembers[cMaxClustGrps];              // number of members in first 5 groups plus a pseudo group to hold sum of group counts in any groups after the 5th
int32_t MaxGrpMembers;                 // maximum number of members in first 5 groups
int32_t GrpMembers;                    // number of group members in current group
int32_t MaxMembersGrpIdx;              // index of a group having the maximum number of members - could be multiple groups, this is the index of the first group 
int32_t NumDGTGrps;                    // DGTs are called over this number of groups which is limited to the 1st 5 groupsuint32_t NumNoiseGrps;                  // number of groups likely contain background noise alleles - having less than m_MinDGTGrpMembers or m_MinDGTGrpPropTotSamples

uint32_t SampleIdx;
uint8_t* ppConsensusPBAs[cMaxClustGrps];
char* pszChrom;

memset(ppConsensusPBAs,0,sizeof(ppConsensusPBAs));
if (m_UsedHGBinSpecs == 0)
	{
	gDiagnostics.DiagOut(eDLWarn, gszProcName, "GenBinKMers: No haplotype group bins to report on!"); // not a failure, but still warn user
	return(eBSFSuccess);
	}
if ((pszChrom = LocateChrom(ChromID)) == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenBinKMers: Unable to locate chromosome for ChromID: '%d'", ChromID);
	return(eBSFerrChrom);
	}


gDiagnostics.DiagOut(eDLInfo, gszProcName, "Processing haplotype grouping KMers for chromosome '%s'", pszChrom);

TotGrpUniqueAlleles = 0;
TotKMerSizeLociGrpsHammingsAccepted = 0;
TotGrpKMerSizeLoci = 0;
pHGBinSpec = m_pHGBinSpecs;
for (BinIdx = 0; BinIdx < m_UsedHGBinSpecs; BinIdx++, pHGBinSpec++)
	{
	if ((pHGBinSpec->ProcState & 0x07) != 0x07) // bin processing state  - 0x00 if unprocessed, 0x01 if currently being processing, 0x03 if processing completed but referenced chromosome not located, 0x07 if processing successfully completed
		continue;
	if (pHGBinSpec->ChromID != ChromID)
		continue;
	pHaplotypeGroup = pHGBinSpec->pHaplotypeGroup;
		// only interested in bins having at least 2 groups
	if (pHaplotypeGroup->ActualHaplotypeGroups < 2)
		continue;

		// determine number of samples in each group
	MaxGrpMembers = 0;
	MaxMembersGrpIdx = 0;
	memset(NumGrpMembers, 0, sizeof(NumGrpMembers));
	for (GrpIdx = 0; GrpIdx < pHaplotypeGroup->ActualHaplotypeGroups; GrpIdx++)
		{
		GrpMembers = BitsVectCount(pHaplotypeGroup->HaplotypeGroup[GrpIdx]);
		if (GrpIdx > 5)
			NumGrpMembers[5] += GrpMembers;
		else
			NumGrpMembers[GrpIdx] = GrpMembers;
		if (GrpIdx < 5 && GrpMembers > MaxGrpMembers) // identify which group has the maximum number of members - could be more than 1 group but this does not matter as group counts are scaled to the maximum members
			{
			MaxGrpMembers = GrpMembers;
			MaxMembersGrpIdx = GrpIdx;
			}
		}

		// determine number of groups to be characterised as being noise, and if less than 2 groups with meaningful allele counts remaining then skip current bin
	if (MaxGrpMembers < m_MinDGTGrpMembers)
		continue;

	NumDGTGrps = min(5, pHaplotypeGroup->ActualHaplotypeGroups); // could be more than 5 groups but only calling KMers on the first 5 groups
	NumNoiseGrps = 0;
	for (GrpIdx = 0; GrpIdx < NumDGTGrps; GrpIdx++)
		if (NumGrpMembers[GrpIdx] < m_MinDGTGrpMembers)
			NumNoiseGrps++;
	if (NumDGTGrps - NumNoiseGrps < 2) // has to be at least two groups remaining after noise groups have been counted
		continue;



	// generate consensus for members of each group over full length of bin
	int32_t NumNonConsensus[cMaxClustGrps];
	memset(NumNonConsensus,0,sizeof(NumNonConsensus));
	for (GrpIdx = 0; GrpIdx < pHaplotypeGroup->ActualHaplotypeGroups; GrpIdx++)
		{
		if (ppConsensusPBAs[GrpIdx] != nullptr)
			{
			delete[]ppConsensusPBAs[GrpIdx];
			ppConsensusPBAs[GrpIdx] = nullptr;
			}
		if (NumGrpMembers[GrpIdx] < m_MinDGTGrpMembers)
			continue;
		ppConsensusPBAs[GrpIdx] = new uint8_t [pHGBinSpec->NumLoci];
		uint8_t **ppPBAs = new uint8_t *[NumGrpMembers[GrpIdx]];
		memset(ppPBAs,0,sizeof(uint8_t*)* NumGrpMembers[GrpIdx]);
		int32_t ConsensusIdx = 0;
		for (SampleIdx = 0; SampleIdx < (uint32_t)pHaplotypeGroup->NumFndrs; SampleIdx++)
			{
			if (!BitsVectTest(SampleIdx, pHaplotypeGroup->HaplotypeGroup[GrpIdx])) // if sample not member of current group then onto next sample
				continue;
			
			ppPBAs[ConsensusIdx++] = pFounderPBAs[SampleIdx] + pHGBinSpec->StartLoci;
			}
		ppConsensusPBAs[GrpIdx] = new uint8_t [pHGBinSpec->NumLoci];
		NumNonConsensus[GrpIdx] = GenGrpConsensus(pHGBinSpec->NumLoci, NumGrpMembers[GrpIdx], ppPBAs, ppConsensusPBAs[GrpIdx]);
		delete[]ppPBAs;
		}

	// have the consensus for each group over the full length bin, iterate over the bin looking for KMers which could segregate between bins
	int32_t KMerLoci;
	int32_t KMerEndLoci;										// KMerEndLoci is inclusive
	int32_t KMerSizeLociGrpsIdentical = 0;
	int32_t KMerSizeLociGrpsHammingsAccepted = 0;
	int32_t LociDiff = 1;
	uint8_t* pGrpRepPBAs[cMaxClustGrps];

	EndLoci = pHGBinSpec->StartLoci + pHGBinSpec->NumLoci - 1;	// EndLoci is inclusive
	for (CurLoci = pHGBinSpec->StartLoci; CurLoci <= EndLoci; CurLoci += LociDiff)
		{
		if((EndLoci - CurLoci + 1) < KMerSize)
			break;

		KMerEndLoci = CurLoci + KMerSize - 1;

		for (GrpIdx = 0; GrpIdx < pHaplotypeGroup->ActualHaplotypeGroups; GrpIdx++)
			{
			if (NumGrpMembers[GrpIdx] < m_MinDGTGrpMembers)
				continue;
			pGrpRepPBAs[GrpIdx] = ppConsensusPBAs[GrpIdx] + (CurLoci - pHGBinSpec->StartLoci);
			}
		KMerLoci = KMerEndLoci;

#ifdef USETHISNONCONSENSUS
		// 2 or more non-noise group sample counts now known
		// members of each individual group must have consensus within the Kmer sequence before processing for hammings between groups
	int NumGrpRepNoCoverage[cMaxClustGrps];
	uint8_t GrpSampleAllele;
	LociDiff = 1;
	KMerSizeLociGrpsIdentical = 0;
	KMerSizeLociGrpsHammingsAccepted = 0;
	EndLoci = pHGBinSpec->StartLoci + pHGBinSpec->NumLoci - 1;	// EndLoci is inclusive
	for (CurLoci = pHGBinSpec->StartLoci; CurLoci <= EndLoci; CurLoci+=LociDiff)
		{
		if ((EndLoci - CurLoci + 1) < KMerSize)
			break;
		KMerEndLoci = CurLoci + KMerSize - 1;

		if(KMerNoneCoverage > 0)
			memset(NumGrpRepNoCoverage,0,sizeof(NumGrpRepNoCoverage));
		for(KMerLoci = CurLoci; KMerLoci <= KMerEndLoci; KMerLoci++)
			{
			for (GrpIdx = 0; GrpIdx < pHaplotypeGroup->ActualHaplotypeGroups; GrpIdx++)		// is there agreement intra group - consensus?
				{
				
				if(NumGrpMembers[GrpIdx] < m_MinDGTGrpMembers)
					continue;
				GrpSampleAllele = 0xff;
				for (SampleIdx = 0; SampleIdx < (uint32_t)pHaplotypeGroup->NumFndrs; SampleIdx++)
					{
					if (!BitsVectTest(SampleIdx, pHaplotypeGroup->HaplotypeGroup[GrpIdx])) // if sample not member of current group then onto next sample
						continue;
					pSampleLoci = pFounderPBAs[SampleIdx] + KMerLoci; // alleles for current sample at CurLoci
					SampleAllele = *pSampleLoci;
					if(KMerNoneCoverage == 0 && SampleAllele == 0)
						break;
					if(GrpSampleAllele == 0xff)		// true if 1st sample in group, 1st sample to become representative
						{
						if(KMerNoneCoverage > 0 && SampleAllele == 0)
							{
							NumGrpRepNoCoverage[GrpIdx] += 1;
							if (NumGrpRepNoCoverage[GrpIdx] > KMerNoneCoverage)
								break;
							}
						GrpSampleAllele = SampleAllele;
						if(KMerLoci == CurLoci)
							pGrpRepPBAs[GrpIdx] = pSampleLoci;
						continue;
						}
					if(SampleAllele != GrpSampleAllele)	// does this sample in same group differ from group allele?
						break;
					}
				if(SampleIdx != pHaplotypeGroup->NumFndrs)
					break;
				}

				// were all members of each group identical within that group at KMerLoci?
			if (GrpIdx != pHaplotypeGroup->ActualHaplotypeGroups)
				{
				LociDiff = 1 + KMerLoci - CurLoci;
				break;
				}
				// all individual group members identical at the current KMerLoci, onto next KMerLoci and check if all individual group members still identical with the KMer sequences to be reported
			}

		// if individual group members have consensus then generate the intergroup hammings between the group consensus sequences
		if(GrpIdx == pHaplotypeGroup->ActualHaplotypeGroups && KMerLoci == KMerEndLoci)
#endif
			// if individual group members have consensus then generate the intergroup hammings between the group consensus sequences
		if (GrpIdx == pHaplotypeGroup->ActualHaplotypeGroups)
			{
			uint8_t *pGrpRefKMerSeq;
			uint8_t* pGrpRelKMerSeq;
			uint8_t RefPBA;
			uint8_t RelPBA;
			int32_t RefGrpIdx;
			int32_t RelGrpIdx;
			int32_t RefRelSeqIdx;
			int32_t RefRelHammingDist;
			int32_t MinRefRelHammingDist;
			int32_t MaxRefRelHammingDist;
			int32_t NumLociNoCoverage;
			LociDiff = KMerSize;
			KMerSizeLociGrpsIdentical++;
				// over the KMerSize sequence, each group must be a minimum hamming away from any other group
			RelGrpIdx = 0;
			MinRefRelHammingDist = 0;
			MaxRefRelHammingDist = 0;
			for (RefGrpIdx = 0; RefGrpIdx < pHaplotypeGroup->ActualHaplotypeGroups-1; RefGrpIdx++)
				{
				if (NumGrpMembers[RefGrpIdx] < m_MinDGTGrpMembers)
					continue;
				for (RelGrpIdx = RefGrpIdx+1; RelGrpIdx < pHaplotypeGroup->ActualHaplotypeGroups; RelGrpIdx++)
					{
					if (NumGrpMembers[RelGrpIdx] < m_MinDGTGrpMembers)
						continue;
					pGrpRefKMerSeq = pGrpRepPBAs[RefGrpIdx];
					pGrpRelKMerSeq = pGrpRepPBAs[RelGrpIdx];

						// what is the hamming distance between the ref and rel group?
					RefRelHammingDist = 0;
					NumLociNoCoverage = 0;
					for(RefRelSeqIdx = 0; RefRelSeqIdx < KMerSize; RefRelSeqIdx++, pGrpRefKMerSeq++, pGrpRelKMerSeq++)
						{
						RefPBA = *pGrpRefKMerSeq;
						RelPBA = *pGrpRelKMerSeq;
						if(RefPBA == 0x00 || RelPBA == 0x00)	// no coverage?
							{
							if(KMerNoneCoverage < ++NumLociNoCoverage)			// allowing loci within KMer to have no coverage?
								{
								RefRelHammingDist = 0;
								LociDiff = RefRelSeqIdx + 1;
								break;
								}
							}
						if(RefPBA != RelPBA)
							RefRelHammingDist++;
						}
					
					if(RefRelHammingDist < MinKMerHammings)
						break;

					if (MinRefRelHammingDist == 0 || RefRelHammingDist < MinRefRelHammingDist)
						MinRefRelHammingDist = RefRelHammingDist;
					if (RefRelHammingDist > MaxRefRelHammingDist)
						MaxRefRelHammingDist = RefRelHammingDist;
					}
				if (RefRelHammingDist < MinKMerHammings)
					break;
				}

#ifdef _CHECKVALIDATEKMERSEQS
			if (RefGrpIdx == pHaplotypeGroup->ActualHaplotypeGroups - 1 && RelGrpIdx == pHaplotypeGroup->ActualHaplotypeGroups)
				{
				int32_t ChkIdx;
				int32_t NumChkGrps;
				uint8_t *GrpPSeqs[5];
				NumChkGrps = 0;
				for (ChkIdx = 0; ChkIdx < pHaplotypeGroup->ActualHaplotypeGroups; ChkIdx++)
					{
					if (NumGrpMembers[ChkIdx] < m_MinDGTGrpMembers)
						continue;
					GrpPSeqs[NumChkGrps++] = pGrpRepPBAs[ChkIdx];
					}
				ValidateArrayKMerSeqs(KMerSize, NumChkGrps, GrpPSeqs);
				}
#endif

			if(RefGrpIdx == pHaplotypeGroup->ActualHaplotypeGroups - 1 && RelGrpIdx == pHaplotypeGroup->ActualHaplotypeGroups && RefRelHammingDist >= MinKMerHammings)
				{
				KMerSizeLociGrpsHammingsAccepted++;
				int32_t KMerIdx;
				uint8_t *pSampleLoci;
				uint8_t *pGrpRefKMerLoci;
				size_t CurKMerSeqOfs;
				size_t KMerSeqOfs;
				int NumKMersToValidate;
				int32_t NumMismatches;
				size_t StartKMerSegOfs;
				uint32_t GrpsKMerMsk = 0x0000;
				uint32_t CurGrpKMerMsk = 0x0001;

				NumKMersToValidate = 0;
				for (RefGrpIdx = 0; RefGrpIdx < pHaplotypeGroup->ActualHaplotypeGroups; RefGrpIdx++, CurGrpKMerMsk<<=1)
					{
					if (NumGrpMembers[RefGrpIdx] < m_MinDGTGrpMembers)
						continue;
					pGrpRefKMerSeq = pGrpRepPBAs[RefGrpIdx];

					NumMismatches = 0;
					for (SampleIdx = 0; SampleIdx < (uint32_t)pHaplotypeGroup->NumFndrs; SampleIdx++)
						{
						if (!BitsVectTest(SampleIdx, pHaplotypeGroup->HaplotypeGroup[RefGrpIdx])) // if sample not member of current group then onto next sample
							continue;
						pSampleLoci = pFounderPBAs[SampleIdx] + CurLoci;
						pGrpRefKMerLoci = pGrpRefKMerSeq;
						for(KMerIdx=0; KMerIdx < KMerSize; KMerIdx++, pSampleLoci++, pGrpRefKMerLoci++)
							if(*pGrpRefKMerLoci != *pSampleLoci)
								NumMismatches++;
						}
					CurKMerSeqOfs = AddKMerSeq(NumGrpMembers[RefGrpIdx],NumMismatches,KMerSize, pGrpRefKMerSeq);
					NumKMersToValidate++;
					if(NumKMersToValidate == 1)
						StartKMerSegOfs = CurKMerSeqOfs;
					if(GrpsKMerMsk == 0)
						KMerSeqOfs = CurKMerSeqOfs;
					GrpsKMerMsk |= CurGrpKMerMsk;
					}
#ifdef _CHECKVALIDATEKMERSEQS
				ValidateKMerSeqs(KMerSize, NumDGTGrps - NumNoiseGrps, &m_pGrpKMerSeqs[StartKMerSegOfs-1],MinRefRelHammingDist, MaxRefRelHammingDist);
#endif
				AddKMerLoci(ChromID,CurLoci,KMerSize, NumDGTGrps - NumNoiseGrps, GrpsKMerMsk, KMerSeqOfs-1, MinRefRelHammingDist,MaxRefRelHammingDist);
				}
			}
		}
	TotGrpKMerSizeLoci += KMerSizeLociGrpsIdentical;
	TotKMerSizeLociGrpsHammingsAccepted += KMerSizeLociGrpsHammingsAccepted;
	}

for (GrpIdx = 0; GrpIdx < pHaplotypeGroup->ActualHaplotypeGroups; GrpIdx++)
	{
	if (ppConsensusPBAs[GrpIdx] != nullptr)
		{
		delete[]ppConsensusPBAs[GrpIdx];
		ppConsensusPBAs[GrpIdx] = nullptr;
		}
	}

gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenBinKMers: There are %d potential grouping KMers of size %d with %d KMers having minimum %d hamming differentials between any two groups on '%s'", TotGrpKMerSizeLoci, KMerSize, TotKMerSizeLociGrpsHammingsAccepted, MinKMerHammings, pszChrom);
return(TotKMerSizeLociGrpsHammingsAccepted);
}

bool
CCallHaplotypes::ValidateArrayKMerSeqs(int32_t KMerSize,		// KMers are this length
	int32_t NumGrpSeqs,						// number of KMer sequences to validate, must be at least 1 difference between any 2 sequences
	uint8_t* pSeqs[],						// sequences to be validated
	int32_t MinHammings,				// claimed minimum hammings between any two KMers; if 0 then don't validate the actual hammings
	int32_t MaxHammings)				// claimed maximum hammings between any two KMers; if 0 then don't validate the actual hammings
{
int32_t RefGrp;
int32_t RelGrp;
uint8_t *pRefPBAs;
uint8_t* pRelPBAs;
int32_t LociOfs;
int32_t NumDiffs;
int32_t MinDiffs;
int32_t MaxDiffs;
MinDiffs = 0;
MaxDiffs = 0;
for (RefGrp = 0; RefGrp < NumGrpSeqs - 1; RefGrp++)
	{
	for (RelGrp = RefGrp + 1; RelGrp < NumGrpSeqs; RelGrp++)
		{
		pRefPBAs = pSeqs[RefGrp];
		pRelPBAs = pSeqs[RelGrp];
		NumDiffs = 0;
		for (LociOfs = 0; LociOfs < KMerSize; LociOfs++, pRefPBAs++, pRelPBAs++)
			{
			if (*pRefPBAs != *pRelPBAs)
				NumDiffs++;
			}
		if(NumDiffs == 0)	// expecting some differences between all sequences!
			return(false);
		if(MinHammings >= 1)
			{
			if (MinDiffs == 0 || NumDiffs < MinDiffs)
				MinDiffs = NumDiffs;
			if (NumDiffs > MaxDiffs)
				MaxDiffs = NumDiffs;
			}
		}
	}
if (MinHammings >= 1)
	{
	if (MinDiffs != MinHammings)
		return(false);
	if (MaxDiffs != MaxHammings)
		return(false);
	}
return(true);
}

bool
CCallHaplotypes::ValidateKMerSeqs(int32_t KMerSize,	// KMers are this length
				int32_t NumGrpSeqs,	// number of KMer sequences to validate 
				uint8_t *pGrpSeqs,	// KMer group sequences, each preceded by number of group members as a int32_t
				int32_t MinHammings, // claimed minimum hammings between any two KMers
				int32_t MaxHammings)	// claimed maximum hammings between any two KMers
{
int32_t RefGrp;
int32_t RelGrp;
uint8_t *pRefPBAs;
uint8_t *pRelPBAs;
uint8_t* pRefPBA;
uint8_t* pRelPBA;
int32_t LociOfs;
int32_t NumDiffs;
int32_t MinDiffs;
int32_t MaxDiffs;
bool bChkGrps;
MinDiffs = 0;
MaxDiffs = 0;
pRefPBAs = pGrpSeqs;
if (NumGrpSeqs > 2)
	bChkGrps = true;
else
	bChkGrps = false;

for(RefGrp = 0; RefGrp < NumGrpSeqs-1; RefGrp++,pRefPBAs += KMerSize)
	{
	pRefPBAs += (2*sizeof(uint32_t));
	pRelPBAs = pRefPBAs+KMerSize;
	for(RelGrp = RefGrp+1; RelGrp < NumGrpSeqs; RelGrp++,pRelPBAs += KMerSize)
		{
		pRelPBAs += (2*sizeof(uint32_t));
		pRefPBA = pRefPBAs;
		pRelPBA = pRelPBAs;
		NumDiffs = 0;
		for(LociOfs = 0; LociOfs < KMerSize; LociOfs++, pRefPBA++, pRelPBA++)
			{
			if(*pRefPBA != *pRelPBA)
				NumDiffs++;
			}
		if(NumDiffs == 0)	// expecting some differences!!
			return(false);
		if(MinDiffs == 0 || NumDiffs < MinDiffs)
			MinDiffs = NumDiffs;
		if (NumDiffs > MaxDiffs)
			MaxDiffs = NumDiffs;
		}
	}
if(MinDiffs != MinHammings)
	return(false);
if (MaxDiffs != MaxHammings)
	return(false);
return(true);
}



tsHaplotypeGroup * // returns allocated ptr to haplotype groups - allocated with 'new', free with 'delete'
CCallHaplotypes::GroupHaplotypes(int32_t ChromID, // grouping is on this chromosome
								int32_t ChromSize, // chromosome is this size
								int32_t StartLoci, // grouping starts at this loci
								int32_t Length,     // grouping is over this many loci
								int32_t MinCentroidDistance, // centroid distance can range from this minimum and up to MaxCentroidDistance when attempting to limit number of haplotype groupings to MaxNumHaplotypeGroups
								int32_t MaxCentroidDistance, // maximum centroid distance which can be used when attempting to meet MaxNumHaplotypeGroups constraint
								int32_t MaxNumHaplotypeGroups, // attempt to constrain number of haplotype groups to be this maximum
								uint32_t* pFndrDiffs, // matrix counts of founder differentials
								int32_t NumFndrs) // number of founders
{
tsBitsVect *pClusterMembers;
tsBitsVect *pMinimisedClusterMembers;
tsBitsVect AssocToAGroup;
int32_t nCol;

uint32_t* pFndrDiff;

int32_t nRow;

int32_t NumHaplotypeGroups;
tsBitsVect CurClusterFounders;
int32_t CurClusterSize;
int32_t TargCentroidDistance;
int32_t TargMaxCentroidDistance;
int32_t TargMinCentroidDistance;
int32_t MinimisedCentroidDistance;
int32_t MaximisedNumHaplotypeGroups;
int32_t CurClusterDiffs;
tsBitsVect SelClusterFounders;
tsBitsVect *pSelClusterMembers = nullptr;
tsHaplotypeGroup* pGroupings = nullptr;
size_t MemReq = 0;
int32_t SelClusterSize;
int32_t SelClusterDiffs;
int32_t NumFndrsChkd;

if (StartLoci < 0 || ChromSize < 1 || Length < 1 || (StartLoci + Length) > ChromSize) // sanity check!
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "GroupHaplotypes: StartLoci : %d, Length: %d, outside of chromosome size: %d for ChromID: %d", StartLoci, Length, ChromSize, ChromID);
	return(nullptr);
	}

if((pClusterMembers = new tsBitsVect[NumFndrs]) == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "GroupHaplotypes: Unable to allocate %d ClusterMembers",NumFndrs);
    return(nullptr);
	}
if((pMinimisedClusterMembers = new tsBitsVect[NumFndrs]) == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "GroupHaplotypes: Unable to allocate %d minimised ClusterMembers", NumFndrs);
	delete []pClusterMembers;
	return(nullptr);
	}

TargMaxCentroidDistance = MaxCentroidDistance;
TargMinCentroidDistance = MinCentroidDistance;
MinimisedCentroidDistance = 0xffffffff; // initialise with a sentinel value, later will be checked if has been overwritten  
MaximisedNumHaplotypeGroups = 0;
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
			pSelClusterMembers = &pClusterMembers[NumHaplotypeGroups++];
			*pSelClusterMembers = SelClusterFounders;
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
		if(NumHaplotypeGroups <= MaxNumHaplotypeGroups && NumHaplotypeGroups >= MaximisedNumHaplotypeGroups)
			{
			MinimisedCentroidDistance = TargCentroidDistance; // it may not be possible to converge onto a centroid distance which breaks into a single increment in the number of groups - let's assume 5 groups is the target
			MaximisedNumHaplotypeGroups = NumHaplotypeGroups; // currently at distance 100 there are 4 groups, but when distance is increased to 101 then there are 7 groups 
			memcpy(pMinimisedClusterMembers, pClusterMembers, sizeof(tsBitsVect) * MaxNumHaplotypeGroups); // so need to remember the minimal distance of 4 at which there were 4 groups ... 
			}
		TargMaxCentroidDistance = TargCentroidDistance - 1;
		}
	}
while(TargMaxCentroidDistance >= TargMinCentroidDistance);

if(MinimisedCentroidDistance != 0xffffffff && MinimisedCentroidDistance != TargCentroidDistance)
	{
	TargCentroidDistance = MinimisedCentroidDistance;
	NumHaplotypeGroups = MaximisedNumHaplotypeGroups;
	memcpy(pClusterMembers,pMinimisedClusterMembers, sizeof(tsBitsVect) * NumHaplotypeGroups);
	}

if (NumHaplotypeGroups > MaxNumHaplotypeGroups) // user may need to increase centroid clustering distance!
	{
	gDiagnostics.DiagOut(eDLWarn, gszProcName,"GroupHaplotypes: Number of haplotype groups %d outside max tageted %d groups  allowed at StartLoci : %d, Length: %d, outside of chromosome size: %d for ChromID: %d", NumHaplotypeGroups, MaxNumHaplotypeGroups, StartLoci, Length, ChromSize, ChromID);
	gDiagnostics.DiagOut(eDLWarn, gszProcName,"GroupHaplotypes: Centroid clustering distance, currently %d, needs to be increased so as to reduce number of groups", MaxCentroidDistance);
	}


MemReq = sizeof(tsHaplotypeGroup) + (sizeof(tsBitsVect) * ((size_t)NumHaplotypeGroups-1));
pGroupings = (tsHaplotypeGroup *) new uint8_t[MemReq];
memset(pGroupings, 0, MemReq);
pGroupings->Size = (int32_t)MemReq;
pGroupings->ChromID = ChromID;
pGroupings->ChromSize = ChromSize;
pGroupings->StartLoci = StartLoci;
pGroupings->NumLoci = Length;
pGroupings->NumFndrs = NumFndrs;
pGroupings->CentroidDistance = TargCentroidDistance;
pGroupings->MinCentroidDistance = MinCentroidDistance;
pGroupings->MaxCentroidDistance = MaxCentroidDistance;
pGroupings->MaxHaplotypeGroups = MaxNumHaplotypeGroups;
pGroupings->MRAGrpChromID = -1;
pGroupings->MRAGrpLoci = -1;
pGroupings->ActualHaplotypeGroups = NumHaplotypeGroups;

memcpy(pGroupings->HaplotypeGroup, pClusterMembers, sizeof(tsBitsVect) * NumHaplotypeGroups);
delete []pClusterMembers;
delete []pMinimisedClusterMembers;
return(pGroupings);
}

int				// < 0 if errors, otherwise success
CCallHaplotypes::GenHaplotypeGroups(tsHGBinSpec *pHGBinSpec,		// clustering using these specs
										  int32_t NumFndrs,			// number of founders to be processed 1st PBA at *pPBAs[0]
										  uint8_t* pFounderPBAs[],	// pPBAs[] pts to chromosome PBAs, for each of the chromosome founder PBAs
										  int32_t ChromSize)		// PBAs for current chromosome being processed are this size
{
int Rslt;
uint32_t	GroupingPhase;  // grouping can be multiphased, in inital (1) phase there is no imputation for unaligned PBA values, in next (2) phase then unaligned PBA values are replaced with group concensus PBA values

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

Rslt = eBSFSuccess;
if (pHGBinSpec->StartLoci < 0 || pHGBinSpec->ChromSize < 1 || pHGBinSpec->ChromSize != ChromSize || pHGBinSpec->NumLoci < 1 || (pHGBinSpec->StartLoci + pHGBinSpec->NumLoci) > ChromSize) // sanity check!
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenHaplotypeGroups: pHGBinSpec->StartLoci : %d, pHGBinSpec->NumLoci: %d, outside of chromosome size: %d for ChromID: %d", pHGBinSpec->StartLoci, pHGBinSpec->NumLoci, ChromSize, pHGBinSpec->ChromID);
	return(eBSFerrParams);
	}

uint8_t* pConsensusPBA;
uint8_t* pConsensusLociPBA;
pConsensusPBA = new uint8_t[pHGBinSpec->NumLoci];
GenConsensusPBA(pHGBinSpec->StartLoci, pHGBinSpec->NumLoci, NumFndrs, pFounderPBAs, pConsensusPBA);

pFndrDiffs = new uint32_t[(int64_t)NumFndrs * NumFndrs]; // organised as [row][col]


int ErrCnt = 0;

uint32_t *pFndrGapDiff;
uint32_t *pFndrGapDiffs = new uint32_t[(int64_t)NumFndrs * NumFndrs]; // organised as [row][col]

pCurHaplotypeGroups = nullptr;
pP1CurHapGroups = nullptr;
GroupingPhase = 0; // starting with phase 0 using all founders consensus for missing PBA alignment loci and then onto subsequent refining phases which instead use group membership to derive the group consensus for missing PBA alignment loci
do {
	EndLoci = pHGBinSpec->StartLoci + pHGBinSpec->NumLoci;

		// generating all vs.all matrix
	memset(pFndrDiffs, 0, ((size_t)NumFndrs * NumFndrs) * sizeof(uint32_t));
	memset(pFndrGapDiffs, 0, ((size_t)NumFndrs * NumFndrs) * sizeof(uint32_t));

	for(AlleleLoci = pHGBinSpec->StartLoci; AlleleLoci < EndLoci; AlleleLoci++)
		{
		if(AlleleLoci >= ChromSize)				// shouldn't be ever, ever!
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenHaplotypeGroups: AlleleLoci %d >= ChromSize %d on chrom %d", AlleleLoci, ChromSize, pHGBinSpec->ChromID);

			Rslt = eBSFerrInternal;
			break;
			}
		PrevFndrIdx = -1;
		pFndrDiff = pFndrDiffs;
		pFndrGapDiff = pFndrGapDiffs;
		pConsensusLociPBA = pConsensusPBA + AlleleLoci - pHGBinSpec->StartLoci;

		for(FndrIdx = 0; FndrIdx < NumFndrs; FndrIdx++)  // iterating founder by row
			{
			pFndrLoci = pFounderPBAs[FndrIdx] + AlleleLoci;

			for(ChkFndrIdx = 0; ChkFndrIdx < NumFndrs; ChkFndrIdx++, pFndrDiff++ , pFndrGapDiff++ )  // iterating founder by col
				{
				pChkFndrLoci = pFounderPBAs[ChkFndrIdx] + AlleleLoci;
				switch(GroupingPhase) { // using a switch() instead of a 'if' only because there may be a need for phase specific processing in the future????
					case 0: // inital phase is using consensus over all founders if no alignment PBA
						if((FndrLociPBA = *pFndrLoci) == 0) // if non-aligned (deletion or no coverage?) then use the consensus allele
							FndrLociPBA = *pConsensusLociPBA;

						if((ChkFndrLociPBA = *pChkFndrLoci) == 0) // if non-aligned (deletion or no coverage?) then use the consensus allele
							ChkFndrLociPBA = *pConsensusLociPBA;
						break;

					default:  // subsequent phases are refining the group consensus if no alignment PBA
						if(FndrIdx != PrevFndrIdx)
							{
							PrevFndrIdx = FndrIdx;
							FndrLociPBA = *pFndrLoci;
							if(m_AffineGapLen != 0 && FndrLociPBA  == 0) // if non-aligned (deletion or no coverage?) then generate and use the consensus allele
								FndrLociPBA = GenFounderConsensusPBA(FndrIdx, pP1CurHapGroups, pHGBinSpec->ChromID, pHGBinSpec->ChromSize,AlleleLoci, NumFndrs, pFounderPBAs);
							}
						ChkFndrLociPBA = *pChkFndrLoci;
						if(m_AffineGapLen != 0 && ChkFndrLociPBA == 0)  // if non-aligned (deletion or no coverage?) then generate and use the consensus allele
							ChkFndrLociPBA = GenFounderConsensusPBA(ChkFndrIdx, pP1CurHapGroups, pHGBinSpec->ChromID, pHGBinSpec->ChromSize,AlleleLoci, NumFndrs, pFounderPBAs);
						break;
						}
				if(m_AffineGapLen >= 0 &&					// m_AffineGapLen == -1 if full length base scoring
					(FndrLociPBA == 0 || ChkFndrLociPBA == 0)) 
					{
					if(m_AffineGapLen == 0 || FndrLociPBA == ChkFndrLociPBA)
						{
						*pFndrGapDiff = 0;
						continue;
						}

					*pFndrGapDiff += 1;
					if (*pFndrGapDiff > (uint32_t)m_AffineGapLen)	// m_AffineGapLen must have been >= 0 and still within a relative gap for execution to have reached here 
						continue;							// scoring is only over first m_AffineGapLen of a relative gap
					}
				else 
					*pFndrGapDiff = 0;
				uint8_t DiffCov;
				if (m_PMode == eMCSHCoverageHapsGrps)
					{
					if(FndrLociPBA > ChkFndrLociPBA)
						DiffCov = FndrLociPBA - ChkFndrLociPBA;
					else
						DiffCov = ChkFndrLociPBA - FndrLociPBA;
					if(DiffCov > 10)
						*pFndrDiff += 1;
					}
				else
					if(FndrLociPBA != ChkFndrLociPBA)
						*pFndrDiff += 1;
				}
			}
		}
	if(Rslt != eBSFSuccess)
		break;

		// have a difference matrix, now group haplotypes
	pCurHaplotypeGroups = GroupHaplotypes(pHGBinSpec->ChromID, pHGBinSpec->ChromSize, pHGBinSpec->StartLoci, pHGBinSpec->NumLoci,pHGBinSpec->MinCentroidDistance, pHGBinSpec->MaxCentroidDistance,pHGBinSpec->MaxNumHaplotypeGroups, pFndrDiffs, NumFndrs);
	if(pCurHaplotypeGroups == nullptr) // null if errors whilst grouping
		{
		if (pP1CurHapGroups != nullptr)
			delete [] pP1CurHapGroups;
		pP1CurHapGroups = nullptr;
		Rslt = eBSFerrInternal;
		break;
		}

	if(pP1CurHapGroups == pCurHaplotypeGroups)		// major problem if this occures!
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenHaplotypeGroups: pP1CurHapGroups and pCurHaplotypeGroups pt to same instance!!!!");
		ErrCnt++;	// not useful, but could help if segfaulting?
		Rslt = eBSFerrInternal;
		break;
		}

	if(pP1CurHapGroups != nullptr)
		{
		pP1CurHapGroups->MRAGrpChromID = -1;
		pP1CurHapGroups->MRAGrpLoci = -1;
		memset(pP1CurHapGroups->MRAGrpConsensus, 0, sizeof(pP1CurHapGroups->MRAGrpConsensus));
		if(pP1CurHapGroups->Size == pCurHaplotypeGroups->Size && !memcmp(pCurHaplotypeGroups, pP1CurHapGroups, pCurHaplotypeGroups->Size)) // number of groups and set membership stable?
			{
			delete[]pP1CurHapGroups;
			pP1CurHapGroups = nullptr;
			break;
			}
		delete[]pP1CurHapGroups;
		pP1CurHapGroups = nullptr;
		}
	pP1CurHapGroups = pCurHaplotypeGroups;
	}
while(++GroupingPhase < m_NumHapGrpPhases);

pHGBinSpec->pHaplotypeGroup = pCurHaplotypeGroups;

if(pFndrDiffs != nullptr)
    delete[]pFndrDiffs;
if (pFndrGapDiffs != nullptr)
	delete[]pFndrGapDiffs;
if(pConsensusPBA != nullptr)
    delete []pConsensusPBA;
return(Rslt);
}

int				// returns number of allele stacks generated
CCallHaplotypes::AlignAlleleStacks(int32_t NumFndrs,			// number of founders to be processed against
									int32_t MaxNumLoci)			// work queue items specify this number of loci for processing by threads
{
int32_t FounderID;
tsCHReadsetMetadata *pReadsetMetadata;
tsCHChromMetadata *pChromMetadata;
uint32_t ChromID;
uint32_t CurChromID;
char* pszChrom;
tsCHWorkQueueEl *pWorkQueueEl;
int32_t CurChromMetadataIdx;
uint8_t **ppPBAs;							// additional is to allow for the progeny readset PBAs
uint8_t *pMskPBA;
int32_t ChromIdx;
int32_t MaximalChromLen;
uint32_t NumUnits;			// NumUnits is the total number of work units to be distributed over available threads for processing an individual chromosome
uint32_t UnitSize;
int WorkerThreadStatus;
gDiagnostics.DiagOut(eDLInfo,gszProcName,"AlignAlleleStacks: Starting to generate allele stacks over %d founders ",NumFndrs);
if(m_pWorkQueueEls != nullptr)				// ensuring only one work queue exists
	delete []m_pWorkQueueEls;
m_pWorkQueueEls = nullptr;
m_AllocWorkQueueEls = 0;
pMskPBA = nullptr;

// determine number of work queue elements required to process the maximal length chromosome
NumUnits = 0;
MaximalChromLen = 0;
for(ChromID = 1; ChromID <= (uint32_t)m_NumChromNames; ChromID++)
	{
	int32_t CurChromSize = m_ChromSizes[ChromID - 1];
	if(CurChromSize > MaximalChromLen && CurChromSize >= MaxNumLoci)
		{
		MaximalChromLen = CurChromSize;
		NumUnits = 1 + (MaximalChromLen + MaxNumLoci - 1) / MaxNumLoci;
		}
	}
if(NumUnits == 0)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "AlignAlleleStacks: Number of thread work units is 0");
	return(0);
	}

// sizing work queue to contain NumUnits elements
m_AllocWorkQueueEls = NumUnits;
m_pWorkQueueEls = new tsCHWorkQueueEl[m_AllocWorkQueueEls];
memset(m_pWorkQueueEls, 0, sizeof(tsCHWorkQueueEl) * m_AllocWorkQueueEls);
m_TotWorkQueueEls = 0;
m_NumQueueElsProcessed = 0;
m_FastSerialise = 0;
pWorkQueueEl = m_pWorkQueueEls;

ppPBAs = new uint8_t * [cMaxFounderReadsets + 1];
memset(ppPBAs,0,sizeof(uint8_t *) * (cMaxFounderReadsets+1));
// iterate along PBAs which are common to all founders, plus mask readset if present, and then initialise elements of work to be undertaken by each thread
CurChromID = 0;
pReadsetMetadata = &m_Readsets[0];
CurChromMetadataIdx = pReadsetMetadata->StartChromMetadataIdx;
for(ChromIdx = 0; ChromIdx < pReadsetMetadata->NumChroms && CurChromMetadataIdx != 0; ChromIdx++)
	{
	// ensure any previously loaded PBAs are deleted
	if(CurChromID > 0)
		{
		for(FounderID = 1; FounderID <= NumFndrs; FounderID++)
			DeleteSampleChromPBAs(FounderID, pChromMetadata->ChromID);
		CurChromID = 0;
		}
	pChromMetadata = &m_pChromMetadata[CurChromMetadataIdx - 1];
	CurChromMetadataIdx = pChromMetadata->NxtChromMetadataIdx;		// ready for next iteration over the chromosomes
	CurChromID = pChromMetadata->ChromID;
	pszChrom = LocateChrom(CurChromID);

	gDiagnostics.DiagOut(eDLInfo, gszProcName, "AlignFounderHaps: Loading chromosome PBAs for '%s' over %d founders", pszChrom, m_NumFounders);
	WorkerThreadStatus = StartWorkerLoadChromPBAThreads(m_NumThreads, 1, m_NumFounders, CurChromID);
	if (WorkerThreadStatus < 0)
		{
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "AlignFounderHaps: Errors loading PBAs for %s chromosome", pszChrom);
		delete[]ppPBAs;
		return(WorkerThreadStatus);
		}
	while (WorkerThreadStatus > 0)
		{
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "AlignFounderHaps: Continuing to load chromosome PBAs for %s .... ", pszChrom);
		WorkerThreadStatus = WaitWorkerThreadStatus(60);
		}
	TerminateLoadChromPBAsThreads(60);
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "AlignFounderHaps: Completed loading chromosome PBAs for %s", pszChrom);

	for(FounderID = 1; FounderID <= m_NumFounders; FounderID++)
		{
		tsCHChromMetadata *pChromMetadata;
		// returned pointer to chromosome metadata
		if((pChromMetadata = LocateChromMetadataFor(FounderID, CurChromID)) == nullptr)
			{
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessGrpLociDGTs: No metadata for chromosome '%s' in at least one founder '%s', skipping this chromosome", pszChrom, LocateReadset(FounderID));
			break;
			}

		if((ppPBAs[FounderID-1] = pChromMetadata->pPBAs) == nullptr)
			{
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessGrpLociDGTs: No PBA for chromosome '%s' in at least one founder '%s', skipping this chromosome", pszChrom, LocateReadset(FounderID));
			break;
			}
		}
	
	if(FounderID <= m_NumFounders)
		{
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessGrpLociDGTs: Ready for next iteration of chromosome DGT processing", pszChrom);
		continue; // slough this chromosome and try next chromosome
		}


	if(m_MaskReadsetID == 0)				// if masking then chromosome PBA for mask also required
		pMskPBA = nullptr;
	else
		if((pMskPBA = LocatePBAfor(m_MaskReadsetID, CurChromID)) == nullptr)
			continue;

	// what is the maximal sized unit that a thread would be expected to process in this chromosome
	if(pChromMetadata->ChromLen < MaxNumLoci)	// too much overhead if too many threads processing shorter chromosomes
		UnitSize = pChromMetadata->ChromLen; // a single thread to process this small chromosome
	else
		UnitSize = MaxNumLoci;

	int32_t StartLoci = 0;
	memset(m_pWorkQueueEls, 0, sizeof(tsCHWorkQueueEl) * m_AllocWorkQueueEls);
	pWorkQueueEl = m_pWorkQueueEls;
	m_NumQueueElsProcessed = 0;
	m_TotWorkQueueEls = 0;
	while(StartLoci < pChromMetadata->ChromLen)
		{
		if(m_TotWorkQueueEls == m_AllocWorkQueueEls)
			break;
		pWorkQueueEl->NumFndrs = NumFndrs;
		pWorkQueueEl->ChromID = CurChromID;
		pWorkQueueEl->StartLoci = StartLoci;
		pWorkQueueEl->ChromLen = pChromMetadata->ChromLen;
		pWorkQueueEl->MaxNumLoci = min(MaxNumLoci, pChromMetadata->ChromLen - StartLoci);
		pWorkQueueEl->pMskPBA = pMskPBA;
		memcpy(pWorkQueueEl->pFounderPBAs,ppPBAs,((size_t)NumFndrs * sizeof(uint8_t *)));
		m_TotWorkQueueEls++;
		pWorkQueueEl++;
		StartLoci += UnitSize;
		}

	WorkerThreadStatus = StartWorkerThreads(m_NumThreads, 1);
	if (WorkerThreadStatus < 0)
		{
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "AlignAlleleStacks: Errors generating allele stacks");
		if(ppPBAs != nullptr)
			delete[]ppPBAs;
		return(WorkerThreadStatus);
		}
	while (WorkerThreadStatus > 0)
		{
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "AlignAlleleStacks: Generating allele stacks for chromosome %s .... ", pszChrom);
		WorkerThreadStatus = WaitWorkerThreadStatus(60);
		}
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "AlignAlleleStacks: Completed generating allele stacks for chromosome %s", pszChrom);
	WorkerThreadStatus=TerminateWorkerThreads(60);
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "AlignAlleleStacks: TerminateWorkThreads(60) returned %d", WorkerThreadStatus);
	}

if(CurChromID > 0)
	for(FounderID = 1; FounderID <= NumFndrs; FounderID++)
		DeleteSampleChromPBAs(FounderID, CurChromID);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"AlignAlleleStacks: Completed generating a total of %u allele stacks over %d founders", m_UsedAlleleStacks,NumFndrs);
if(m_pWorkQueueEls != nullptr)
	{
	delete []m_pWorkQueueEls;
	m_pWorkQueueEls = nullptr;	
	}
if(ppPBAs != nullptr)
	delete []ppPBAs;
m_NumQueueElsProcessed = 0;
m_FastSerialise = 0;
return(m_UsedAlleleStacks);
}



int				
CCallHaplotypes::AlignFounderHaps(int32_t NumFndrs)			// number of founders to be processed against
{
int Rslt;
int32_t FounderID;
tsCHReadsetMetadata *pReadsetMetadata;
tsCHChromMetadata *pChromMetadata;
int32_t CurChromID;
char* pszChrom;
int32_t StartLoci;
int32_t BinSize;

tsCHWorkQueueEl *pWorkQueueEl;
int32_t CurChromMetadataIdx;
uint8_t **ppPBAs;
uint8_t *pMskPBA;
int32_t ChromIdx;
uint32_t NumUnits;			// NumUnits is the total number of work units to be distributed over available threads for processing chromosomes
bool bInitHGBinSpecs;


pReadsetMetadata = &m_Readsets[0];		// founders were loaded first so 1st founder will be here
CurChromMetadataIdx = pReadsetMetadata->StartChromMetadataIdx;
if(m_pWorkQueueEls != nullptr)				// ensuring only one work queue exists
    delete[]m_pWorkQueueEls;
m_pWorkQueueEls = nullptr;
m_AllocWorkQueueEls = 0;
pMskPBA = nullptr;

Rslt = eBSFSuccess; // assume no errors!!!

// if cluster group bins not already specified by user then will need to later initialise 
bInitHGBinSpecs = m_UsedHGBinSpecs == 0 ? true : false;

// determine maximal number of work queue elements required, as will be processing 1 chromosome at a time then only 1 work element is required!
CurChromMetadataIdx = pReadsetMetadata->StartChromMetadataIdx;

// sizing work queue to contain a maximum of NumUnits elements
// one work unit corresponds to haplotype processing for a single chromosome at a time, and chromosomes are iterated
NumUnits = 1;
m_AllocWorkQueueEls = NumUnits + 1;  // 1 additional as a safety measure - shouldn't be any accesses beyond the one and only work element but ....
m_pWorkQueueEls = new tsCHWorkQueueEl[m_AllocWorkQueueEls];
memset(m_pWorkQueueEls, 0, sizeof(tsCHWorkQueueEl) * m_AllocWorkQueueEls);
m_TotWorkQueueEls = 0;
m_NumQueueElsProcessed = 0;
m_FastSerialise = 0;
pWorkQueueEl = m_pWorkQueueEls;

if(bInitHGBinSpecs)
	{
	for(ChromIdx = 0; ChromIdx < pReadsetMetadata->NumChroms && CurChromMetadataIdx != 0; ChromIdx++)
		{
		pChromMetadata = &m_pChromMetadata[CurChromMetadataIdx - 1];
		CurChromMetadataIdx = pChromMetadata->NxtChromMetadataIdx;		// ready for next iteration over the chromosomes
		CurChromID = pChromMetadata->ChromID;
		pszChrom = LocateChrom(CurChromID);

		// create default fixed sizes
		for(StartLoci = 0; StartLoci < pChromMetadata->ChromLen; StartLoci += m_GrpHapBinSize)
			{
			BinSize = min(pChromMetadata->ChromLen - StartLoci, m_GrpHapBinSize);
			if(BinSize < 10) // any bin smaller is of no value!
				continue;
			if((Rslt = AddHGBinSpec(CurChromID, pChromMetadata->ChromLen,StartLoci, BinSize, min(m_MinCentClustDist, BinSize - 1), min(m_MaxCentClustDist, BinSize - 1), m_MaxClustGrps)) == 0)
				{
				gDiagnostics.DiagOut(eDLInfo, gszProcName, "AlignFounderHaps: Error adding haplotyping group bin specifications for chromosome '%s'", pszChrom);
				if(m_pWorkQueueEls != nullptr)
					{
					delete[]m_pWorkQueueEls;
					m_pWorkQueueEls = nullptr;
					}
				m_NumQueueElsProcessed = 0;
				return(-1);
				}
			if(StartLoci == 0) // capture bin identifier for 1st bin on current chromosome or sequence
				pChromMetadata->HGBinID = Rslt;
			}
		}
	}


// iterate along PBAs which are common to all founders and then initialise elements of work to be undertaken by each thread

ppPBAs = new uint8_t * [cMaxFounderReadsets + 1];
memset(ppPBAs, 0, sizeof(uint8_t*) * (cMaxFounderReadsets + 1));

m_ProccessingHGBinID = 0;
CurChromID = 0;
CurChromMetadataIdx = pReadsetMetadata->StartChromMetadataIdx;
for(ChromIdx = 0; ChromIdx < pReadsetMetadata->NumChroms && CurChromMetadataIdx != 0; ChromIdx++)
	{
		// ensure any previously loaded PBAs are deleted
	if(CurChromID > 0)
		{
		for(FounderID = 1; FounderID <= NumFndrs; FounderID++)
			DeleteSampleChromPBAs(FounderID, CurChromID);
		CurChromID = 0;
		}

	pChromMetadata = &m_pChromMetadata[CurChromMetadataIdx - 1];
	CurChromMetadataIdx = pChromMetadata->NxtChromMetadataIdx;		// ready for next iteration over the chromosomes
	if(pChromMetadata->HGBinID == 0) // if no haplotype grouping bin specs for this chrom then skip to next chrom
		continue;

	CurChromID = pChromMetadata->ChromID;
	pszChrom = LocateChrom(CurChromID);

	// load this chromosome PBAs from all founders, all founders must have PBAs for this chromosome otherwise chromosome will be skipped

	gDiagnostics.DiagOut(eDLInfo, gszProcName, "AlignFounderHaps: Loading %s on '%s' for %d founders", m_PMode == eMCSHCoverageHapsGrps ? "WIGs" : "PBAs", pszChrom, NumFndrs);
	int WorkerThreadStatus = StartWorkerLoadChromPBAThreads(m_NumThreads, 1, NumFndrs, CurChromID);
	if (WorkerThreadStatus < 0)
		{
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "AlignFounderHaps: Errors loading %s for %s chromosome", m_PMode == eMCSHCoverageHapsGrps ? "WIGs" : "PBAs", pszChrom);
		delete[]ppPBAs;
		return(WorkerThreadStatus);
		}
	while (WorkerThreadStatus > 0)
		{
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "AlignFounderHaps: Continuing to load chromosome PBAs for %s .... ", pszChrom);
		WorkerThreadStatus = WaitWorkerThreadStatus(60);
		}
	TerminateLoadChromPBAsThreads(60);
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "AlignFounderHaps: Completed loading chromosome PBAs for %s", pszChrom);

	for(FounderID = 1; FounderID <= NumFndrs; FounderID++)
		{
		// returned pointer to chromosome metadata
		if((pChromMetadata = LocateChromMetadataFor(FounderID, CurChromID)) == nullptr)
			{
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "AlignFounderHaps: No metadata for chromosome '%s' in at least one founder '%s', skipping this chromosome", pszChrom, LocateReadset(FounderID));
			break;;
			}

		if((ppPBAs[FounderID-1] = pChromMetadata->pPBAs) == nullptr)
			{
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "AlignFounderHaps: No %s for chromosome '%s' in at least one founder '%s', skipping this chromosome", m_PMode == eMCSHCoverageHapsGrps ? "WIGs" : "PBAs", pszChrom, LocateReadset(FounderID));
			break;
			}
		}
	
	if(FounderID <= NumFndrs)
		continue; // slough this chromosome and try next chromosome

	// at least one bin on the current chromosome, initialise the one and only work queue element for this chromosome
	pWorkQueueEl = m_pWorkQueueEls;
	pWorkQueueEl->NumFndrs = NumFndrs;
	pWorkQueueEl->ChromID = CurChromID;
	pWorkQueueEl->StartLoci = 0;
	pWorkQueueEl->ChromLen = pChromMetadata->ChromLen;
	pWorkQueueEl->MaxNumLoci = pChromMetadata->ChromLen;
	pWorkQueueEl->pMskPBA = pMskPBA;
	memcpy(pWorkQueueEl->pFounderPBAs,ppPBAs,((size_t)NumFndrs * sizeof(uint8_t *)));
	m_TotWorkQueueEls = 1;

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"AlignFounderHaps: Starting to generate haplotype groupings on '%s' for %d founders ",pszChrom,NumFndrs);
	WorkerThreadStatus = StartWorkerThreads(m_NumThreads, 1);
	if (WorkerThreadStatus < 0)
		{
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "AlignFounderHaps: Errors Generating haplotype groupings on chrom '%s'", pszChrom);
		if (ppPBAs != nullptr)
			delete[]ppPBAs;
		if (m_pWorkQueueEls != nullptr)
			{	
			delete[]m_pWorkQueueEls;
			m_pWorkQueueEls = nullptr;
			}
		return(WorkerThreadStatus);
		}
	int WaitSecs = 60; // initially wait 1 minute before reporting progress, but increment by 1 minute up to a max of 5 minutes to avoid swamping the logs with progress messages
	while (WorkerThreadStatus > 0)
		{
		WorkerThreadStatus = WaitWorkerThreadStatus(60);
		WaitSecs = min(WaitSecs += 60, 300);
		AcquireFastSerialise();
		int32_t LatestBinID = m_ProccessingHGBinID;
		ReleaseFastSerialise();
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "AlignFounderHaps: Generating haplotype groupings on chromosome '%s' ... %0.3f%% of total haplotype grouping specification bins processed", pszChrom, (LatestBinID * 100.0) / m_UsedHGBinSpecs);
		}
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "AlignFounderHaps: Completed generating all haplotype groupings on chromosome %s", pszChrom);
	WorkerThreadStatus = TerminateWorkerThreads(60);
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "AlignFounderHaps: TerminateWorkThreads(60) returned %d", WorkerThreadStatus);
	}

gDiagnostics.DiagOut(eDLInfo, gszProcName, "AlignFounderHaps: Completed generating haplotype groupings over all chromosomes");

if (ppPBAs != nullptr)
	delete[]ppPBAs;

if(m_pWorkQueueEls != nullptr)
	{
	delete []m_pWorkQueueEls;
	m_pWorkQueueEls = nullptr;	
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
int32_t ASIdx;
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
CCallHaplotypes::BitsVectSet(uint16_t Bit,		// bit to set, range 0..cMaxBitVectBits-1
			tsBitsVect & BitsVect)
{
BitsVect.Bits[Bit / 64] |= ((uint64_t)0x01 << (Bit % 64));
}

inline void
CCallHaplotypes::BitsVectReset(uint16_t Bit,		// bit to reset, range 0..cMaxBitVectBits-1
		   tsBitsVect & BitsVect)
{
BitsVect.Bits[Bit / 64] &= ~((uint64_t)0x01 << (Bit % 64));
}


inline bool
CCallHaplotypes::BitsVectTest(uint16_t Bit,		// bit to test, range 0..cMaxBitVectBits-1
			tsBitsVect & BitsVect)
{
return(BitsVect.Bits[Bit / 64] & ((uint64_t)0x01 << (Bit % 64)) ? true : false);
}

inline bool
CCallHaplotypes::BitsVectEqual(tsBitsVect & BitsVectA,	// compare for equality
								tsBitsVect & BitsVectB)
{
if(!memcmp(&BitsVectA,&BitsVectB,sizeof(tsBitsVect)))
	return(true);
return(false);
}

inline void
CCallHaplotypes::BitsVectInitialise(bool Set,			// if true then initialise all bits as set, otherwise initialise all bits as reset
				tsBitsVect & BitsVect)
{
uint8_t InitWith;
if(Set)
	InitWith = 0xff;
else
	InitWith = 0;
memset(&BitsVect,InitWith,sizeof(tsBitsVect));
}

inline uint32_t
CCallHaplotypes::BitsVectCount(tsBitsVect& BitsVect)		// count number of set bits
{
uint32_t Count = 0;
uint64_t *pWord = BitsVect.Bits;
uint64_t Bit;
for(uint32_t WordIdx = 0; WordIdx < cBitVectWords; WordIdx++, pWord++)
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
CCallHaplotypes::BitsVectUnion(tsBitsVect &BitsVectA, tsBitsVect& BitsVectB)		// union (effective BitsVectA |= BitsVectB) bits in BitsVectA with BitsVectB with BitsVectA updated, returns number of bits set in BitsVectA 
{
uint32_t Count = 0;
uint64_t* pWordA = BitsVectA.Bits;
uint64_t* pWordB = BitsVectB.Bits;
uint64_t Bit;
for(uint32_t WordIdx = 0; WordIdx < cBitVectWords; WordIdx++, pWordA++, pWordB++)
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
CCallHaplotypes::BitsVectIntersect(tsBitsVect& BitsVectA, tsBitsVect& BitsVectB)	// intersect (effective BitsVectA &= BitsVectB) of bits in BitsVectA with BitsVectB with BitsVectA updated, returns number of set bits in BitsVectA
{
uint32_t Count = 0;
uint64_t* pWordA = BitsVectA.Bits;
uint64_t* pWordB = BitsVectB.Bits;
uint64_t Bit;
for(uint32_t WordIdx = 0; WordIdx < cBitVectWords; WordIdx++, pWordA++, pWordB++)
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
CCallHaplotypes::BitsVectClear(tsBitsVect& BitsVectA, tsBitsVect& BitsVectB)	    // clear bits in BitsVectA which are set in BitsVectB with BitsVectA updated, returns number of set bits in BitsVectA  
{
uint32_t Count = 0;
uint64_t* pWordA = BitsVectA.Bits;
uint64_t* pWordB = BitsVectB.Bits;
uint64_t Bit;
for(uint32_t WordIdx = 0; WordIdx < cBitVectWords; WordIdx++, pWordA++, pWordB++)
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

int32_t		// returned chrom identifier, 0 if unable to accept this chromosome name
CCallHaplotypes::AddChrom(char* pszChrom) // associate unique identifier with this chromosome name
{
int32_t ChromNameIdx;
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


int32_t		// returned chrom identifier, 0 if unable to locate this chromosome name
CCallHaplotypes::LocateChrom(char* pszChrom) // return unique identifier associated with this chromosome name
{
int32_t ChromNameIdx;
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
CCallHaplotypes::LocateChrom(int32_t ChromID)
{
if(ChromID < 1 || ChromID > m_NumChromNames)
	return(nullptr);
return(&m_szChromNames[m_szChromIdx[ChromID-1]]);
}

uint8_t 
CCallHaplotypes::LocateReadsetChromLociAlleles(int32_t ReadsetID,	// return alleles for this readset 
									  int32_t ChromID,		// on this chromosome
									  int32_t Loci)		// at this loci
{
uint8_t *pPBA;
if((pPBA=LocatePBAfor(ReadsetID,ChromID))==nullptr)
	return(0);
	
return(pPBA[Loci]);
}

// NOTE: Readsets are checked for uniqueness as readsets must be unique within a given readetset type
int32_t		// returned readset identifier, 0 if unable to accept this readset name
CCallHaplotypes::AddReadset(char* pszReadset, // associate unique identifier with this readset name
							uint8_t ReadsetType)	// 0: founder, 1: progeny, 2: control
{
int32_t ReadsetNameIdx;
int ReadsetNameLen;
char Type;
char *pszLAname;
Type = '0' + (char)ReadsetType;
if(m_NumReadsetNames == 0)
	{
	memset(m_NumReadsetTypes,0,sizeof(m_NumReadsetTypes));
	memset(m_FndrIDMappings,0,sizeof(m_FndrIDMappings));
	}

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
m_NumReadsetTypes[ReadsetType]++;
if(ReadsetType == 0)			// founder types are special, need to maintain mappings of readset identifiers
	m_FndrIDMappings[m_NumReadsetTypes[ReadsetType]-1]= m_LAReadsetNameID;
return(m_LAReadsetNameID);
}

int32_t	// returned readset identifier index for founder type readsets, -1 if unable to locate readset identifier
CCallHaplotypes::LocateReadsetIDIdx(int32_t ReadsetID,		// requiring idx for this readset identifier
					uint8_t ReadsetType)	// 0: founder, 1: progeny, 2: control
{
int32_t ReadsetIDIdx;
if(ReadsetType != 0 || m_NumReadsetNames == 0 || m_NumReadsetTypes[ReadsetType] == 0)
	return(-1);
// currently a linear search but expecting only a few hundred readset identifers to be founder types so shouldn't be nuch of an overhead
for(ReadsetIDIdx = 0; ReadsetIDIdx < m_NumReadsetTypes[ReadsetType]; ReadsetIDIdx++)
	if(ReadsetID == m_FndrIDMappings[ReadsetIDIdx])
		return(ReadsetIDIdx);
return(-1);
}

int32_t		// returned Readset identifier, 0 if unable to locate this Readset name
CCallHaplotypes::LocateReadset(char* pszReadset, // return unique identifier associated with this Readset name
							   uint8_t ReadsetType)	// 0: founder, 1: progeny, 2: control
{
int32_t ReadsetNameIdx;
char Type;
char *pszLAReadset;

if(m_NumReadsetNames > 0 && m_NumReadsetTypes[ReadsetType] == 0)
	return(0);

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

char* 
CCallHaplotypes::LocateReadset(int32_t ReadsetID)
{
int32_t Idx;
Idx = ReadsetID & 0x00ffffff;			// mask out any potential ReadsetType
if(Idx < 1 || Idx > m_NumReadsetNames)
	return(nullptr);
return(&(m_szReadsetNames[m_szReadsetIdx[Idx -1]+1])); // skipping lead char which is the ReadsetType
}


int
CCallHaplotypes::LoadPBACoverage(char* pszInWIGFile)   // file containing PBA, file name extension will be replaced with 'covsegs.wig' which will be expected to be name of file containing the WIG coverage
{
int32_t ReadsetID;
FILE* pInStream;
char szInWIG[1000];
char* pszInWIG;

// compose WIG file name from PBA file name and check if file can be opened
CUtility::AppendFileNameSuffix(szInWIG, pszInWIGFile, (char*)".covsegs.wig", '.');
pszInWIG = szInWIG;
if((pInStream = fopen(pszInWIG, "r")) == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadPBACoverage: Unable to open WIG file %s for reading, error: %s", pszInWIG, strerror(errno));
	Reset();
	return(eBSFerrOpnFile);
	}
fclose(pInStream);

// load PBA file, but skip allocation and loading the actual PBAs
if((ReadsetID = LoadPBAFile(pszInWIGFile, 0, true)) <= 0)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadPBACoverage: Errors loading WIG file '%s'", pszInWIGFile);
	fclose(pInStream);
	Reset();
	return(ReadsetID);
	}
m_Fndrs2Proc[ReadsetID - 1] = 0x01;	// if loaded then assumption is that this founder will be processed
return(ReadsetID);
}


// Note: an assumption is that within the WIG input file ordering is by chromosome ascending
uint8_t * // loaded PBAs for requested chromosome or nullptr if unable to load
CCallHaplotypes::LoadPBAChromCoverage(int32_t ReadsetID, // loading chromosome coverage for this readset and mapping coverage as though PBAs
				int32_t ChromID)    // coverage is for this chromosome
{
int Rslt;

char szInWIG[_MAX_PATH];
char* pszInWIG;
char szLineBuff[200];
tsCHReadsetMetadata* pReadsetMetadata;
tsCHChromMetadata* pChromMetadata;
char* pszReadset;
FILE* pInStream;
uint32_t LineNumb;
uint32_t WIGtype;    // 0: unknown, 1: fixed steps, 2 variable steps
int32_t CurChromID;
uint32_t Coverage;
int32_t StartLoci;
int32_t Span;
uint8_t* pAllele;
int32_t StartLociOfs;
char *pszChrom;
char szChrom[100];
char* pszInBuff;
uint32_t NumElsParsed;

// compose WIG file name from PBA file name
pReadsetMetadata = &m_Readsets[ReadsetID - 1];
CUtility::AppendFileNameSuffix(szInWIG, pReadsetMetadata->szFileName, (char*)".covsegs.wig", '.');
pszInWIG = szInWIG;

if((pszReadset = LocateReadset(ReadsetID)) == nullptr)
     return(nullptr);
if(ChromID == 1)
     pReadsetMetadata->NxtFileChromOfs = 0;

if((pszChrom = LocateChrom(ChromID)) == nullptr)
    return(nullptr);

pChromMetadata = LocateChromMetadataFor(ReadsetID, ChromID);
if(pChromMetadata->pPBAs == nullptr)
	{
	if((pChromMetadata->pPBAs = AllocPBAs(pChromMetadata->ChromLen)) == nullptr)
		return(nullptr);
	}
memset(pChromMetadata->pPBAs, 0, pChromMetadata->ChromLen);

// WIG file can now be opened and coverage for requested chromosome parsed and mapped as though PBAs
if((pInStream = fopen(pszInWIG,"r"))==nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadPBACoverage: Unable to open WIG file %s for reading, error: %s",pszInWIG,strerror(errno));
	return(nullptr);
	}
#if _WIN32
_fseeki64(pInStream, pReadsetMetadata->NxtFileChromOfs, SEEK_SET);
#else
fseeko64(pInStream, pReadsetMetadata->NxtFileChromOfs, SEEK_SET);
#endif
// load the WIG and where there is coverage then set the corresponding loci in the PBA to be max(254,log2(coverage) * 10)
// Note that the max limit of 254 is so that when later processing grouping and caching highest frequency coverage (alleles) then 254 can be distinguished from a marker of 255 used if not already cached
// Note: an assumption is that within the WIG input file ordering is by chromosome ascending
LineNumb = 0;
WIGtype = 0;
CurChromID = 0;
Rslt = eBSFSuccess;

// allocate temporary buffer to hold actual counts without any transformations, this buffer is used for smothing prior to transformation
uint32_t * pRawCnt;
uint32_t *pRawCnts = nullptr;
if((pRawCnts = new uint32_t [pChromMetadata->ChromLen]) == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadPBACoverage: Error allocating buffer hold actual coverage counts for chromosome of length %d prior to smoothing", pChromMetadata->ChromLen);
	return(nullptr);
	}
memset(pRawCnts,0,sizeof(uint32_t) * pChromMetadata->ChromLen);

while(fgets(szLineBuff, sizeof(szLineBuff) - 1, pInStream) != nullptr)
	{
	LineNumb++;
	pszInBuff = CUtility::ReduceWhitespace(szLineBuff);
	if(pszInBuff == nullptr || *pszInBuff == '\0' || *pszInBuff == '\n')
		continue;

	if(pszInBuff[0] == 'v' && !strnicmp(pszInBuff, "variableStep", 12))
		{
		NumElsParsed = sscanf(pszInBuff, "variableStep chrom=%s span=%d", szChrom, &Span);
		if(NumElsParsed != 2)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadPBACoverage: Errors parsing WIG line %d - \"%s\" - in WIG file '%s'",LineNumb,pszInBuff,pszInWIG);
			Rslt = -1;
			break;
			}

		// check for valid chrom, a missing chrom relative to PBA could be that the user has specified a limit on number of chromosomes to process so can't treat as a fatal error
		if((CurChromID = LocateChrom(szChrom)) < 1)
			{
			CurChromID = 0;
			continue;
			}

		if(CurChromID > ChromID) // onto coverage for next chromosome?
			break;

		if(Span < 1 || Span > pChromMetadata->ChromLen)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadPBACoverage: Span must be between 1..%d \"%s\" at line %d - \"%s\" - for chromosome '%s' in WIG file '%s'",pChromMetadata->ChromLen,LineNumb,pszInBuff,pszChrom,pszInWIG);
			Rslt = -1;
			break;
			}
#if _WIN32
		pReadsetMetadata->NxtFileChromOfs = _ftelli64(pInStream);
#else
		pReadsetMetadata->NxtFileChromOfs = ftello64(pInStream);
#endif

		WIGtype = 2;
		continue;
		}
	else
		if(pszInBuff[0] == 'f' && !strnicmp(pszInBuff, "fixedStep",10))
			{
			NumElsParsed = sscanf(pszInBuff, "fixedStep chrom=%s start=%d step=%d", szChrom, &StartLoci,&Span);
			if(NumElsParsed != 3)
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadPBACoverage: Errors parsing line %d - \"%s\" - in WIG file '%s'",LineNumb,pszInBuff,pszInWIG);
				Rslt = -1;
				break;
				}

			// check for valid chrom, a missing chrom relative to PBA could be that the user has specified a limit on number of chromosomes to process so can't treat as a fatal error
			if((CurChromID = LocateChrom(szChrom)) < 1)
				{
				CurChromID = 0;
				continue;
				}

			if(CurChromID > ChromID) // onto coverage for next chromosome?
				break;

			if(StartLoci < 1 || (StartLoci + Span - 1) >  pChromMetadata->ChromLen)
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadPBACoverage: StartLoci must be between 1..%d, Span %d extending past end of chromosome \"%s\" at line %d - \"%s\" - in WIG file '%s'",pChromMetadata->ChromLen,Span,pszChrom,LineNumb, pszInBuff,pszInWIG);
				Rslt = -1;
				break;
				}
#if _WIN32
		pReadsetMetadata->NxtFileChromOfs = _ftelli64(pInStream);
#else
		pReadsetMetadata->NxtFileChromOfs = ftello64(pInStream);
#endif
			WIGtype = 1;
			continue;
			}
	if(CurChromID == 0) // 0 if yet to parse out coverage for targeted chromosome
		{
#if _WIN32
		pReadsetMetadata->NxtFileChromOfs = _ftelli64(pInStream);
#else
		pReadsetMetadata->NxtFileChromOfs = ftello64(pInStream);
#endif
		continue;
		}

	switch(WIGtype) {
			case 0: // unrecognised format, assume some kind of comment or directive, slough
				CurChromID = 0;
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
	pRawCnt = &pRawCnts[StartLoci - 1];
	for(StartLociOfs = 0; StartLociOfs < Span; StartLociOfs++, StartLoci++, pRawCnt++)
		{
		if(StartLoci > pChromMetadata->ChromLen) // ensure not writing past end of chromosome
			break;
		*pRawCnt = Coverage;
		}
	}
// apply smoothing using a window size which is 10% of bin size
uint32_t CurWinSize;
uint32_t MaxWinSize;
uint32_t WinMeanCoverage;
uint32_t WinSumCoverage;
uint32_t *pWinBeginRawCnt;
uint32_t * pWinEndRawCnt;
MaxWinSize = max(10,m_GrpHapBinSize / 10);
CurWinSize = 0;
pWinBeginRawCnt = pRawCnts;
pWinEndRawCnt =pRawCnts;
WinSumCoverage = 0;
pAllele = pChromMetadata->pPBAs;
for (StartLociOfs = 0; StartLociOfs < pChromMetadata->ChromLen; StartLociOfs++,pAllele++, pWinBeginRawCnt++, pWinEndRawCnt++)
	{
	if(CurWinSize < MaxWinSize)
		{
		CurWinSize += 1;
		pWinBeginRawCnt = pRawCnts;
		}
	else
		WinSumCoverage -= *pWinBeginRawCnt;
	WinSumCoverage += *pWinEndRawCnt;
	WinMeanCoverage = WinSumCoverage / CurWinSize;
	*pAllele = min((int)log2(WinMeanCoverage) * 10, 0x0fe);
	}
fclose(pInStream);
if(pRawCnts != nullptr)
	delete []pRawCnts;

if(Rslt < 0)
	{
	DeleteSampleChromPBAs(ReadsetID, CurChromID);
	pChromMetadata->pPBAs = nullptr;
	}
return(pChromMetadata->pPBAs);
}

int
CCallHaplotypes::CSV2WIG(char* pszInFile, // file containing haplotype groupings
	char* pszOutFile) // where to write as WIG formated file
{
int Rslt;
int32_t CurLineNumber;
int32_t NumFields;
char* pszChrom;
uint32_t CurChromID;
uint32_t PrevChromID;
int32_t BinStartLoci;
int32_t BinSize;
int32_t CentroidDistance;
uint32_t NumChromBins;
uint32_t TotNumChromBins;

if(m_hWIGOutFile != -1)
	{
	close(m_hWIGOutFile);
	m_hWIGOutFile = -1;
	}

if(m_pHGCSVFile != nullptr) // should be, but better to be sure!
	{
	delete m_pHGCSVFile;
	m_pHGCSVFile = nullptr;
	}

if((m_pHGCSVFile = new CCSVFile) == nullptr)
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

if(m_pszWIGBuff == nullptr)
	{
	if((m_pszWIGBuff = new uint8_t[cOutBuffSize]) == nullptr)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "CSV2WIG: Memory allocation of %u bytes for WIG output buffering failed - %s", cOutBuffSize, strerror(errno));
		Reset();
		return(eBSFerrMem);
		};
	m_AllocWIGBuff = cOutBuffSize;
	m_WIGBuffIdx = 0;
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
	if(pszChrom == nullptr || pszChrom[0] == '\0')
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
	if (BinStartLoci < 0)
		{
		gDiagnostics.DiagOut(eDLWarn, gszProcName, "Negative loci %d at line %d in file '%s'", BinStartLoci, CurLineNumber, pszInFile);
		continue;
		}
	m_pHGCSVFile->GetInt(3, &BinSize);
	if (BinSize < 1)
		{
		gDiagnostics.DiagOut(eDLWarn, gszProcName, "Bin size %d less than 1 at line %d in file '%s'", BinSize, CurLineNumber, pszInFile);
		continue;
		}
	m_pHGCSVFile->GetInt(7, &CentroidDistance);
	if (CentroidDistance < 1)
		{
		gDiagnostics.DiagOut(eDLWarn, gszProcName, "Centroid distancee %d less than 1 at line %d in file '%s'", CentroidDistance, CurLineNumber, pszInFile);
		continue;
		}

	NumChromBins++;
	AccumWIGCnts(CurChromID, BinStartLoci+1, BinSize, CentroidDistance); // WIG loci start from 1 not 0
	}
CompleteWIGSpan(true);
delete m_pHGCSVFile;
m_pHGCSVFile = nullptr;

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


int64_t									// returned index+1 into m_pASBins[] to allocated and initialised allele scores, 0 if errors
CCallHaplotypes::AddASBin(int32_t SrcID,			// identifies the source PBA which was aligned against
	int32_t RefID,				// this reference PBA scoring allelic differences
	int32_t ChromID,			// scoring was on this chromosome
	int32_t ChromSize,			// chromosome is this size
	int32_t BinLoci,			// scored bin starts at this loci
	int32_t BinSize,			// bin is this size
	double ScoreExact,			// scoring for exact matches
	double ScorePartial)		// scoring both eact and partial matches
{
tsASBin Tmp;
if (BinLoci < 0 || BinSize < 1 || (BinLoci + BinSize) > ChromSize) // sanity check!
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddASLoci: Bin loci range (BinLoci: %d BinSize: %d) extends past chromosome size: %d for ChromID: %d", BinLoci, BinSize, ChromSize, ChromID);
	return(0);
	}
Tmp.SrcID = SrcID;
Tmp.RefID = RefID;
Tmp.ChromID = ChromID;
Tmp.ChromSize = ChromSize;
Tmp.BinLoci = BinLoci;
Tmp.BinSize = BinSize;
Tmp.ScoreExact = ScoreExact;
Tmp.ScorePartial = ScorePartial;
Tmp.ProcState = 0; // bin processing state  - 0x00 if unprocessed, 0x01 if currently being processing, 0x03 if processing completed but referenced chromosome not located, 0x07 if processing successfully completed
return(AddASBin(&Tmp));
}

int64_t									// returned index+1 into m_pASBins[] to allocated and initialised allele scores, 0 if errors
CCallHaplotypes::AddASBin(tsASBin* pInitASBin)	// allocated tsASBin to be initialised with a copy of pInitASBin
{
int64_t ToAllocdASBins;
tsASBin* pASBin;
size_t memreq;
AcquireFastSerialise();
if (m_pASBins == nullptr)					// may be nullptr first time in
	{
	memreq = (size_t)cAllocASBins * sizeof(tsASBin);
#ifdef _WIN32
	m_pASBins = (tsASBin*)malloc((size_t)memreq);
	if (m_pASBins == nullptr)
		{
		ReleaseFastSerialise();
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddASBin: Memory allocation of %zd bytes failed", (int64_t)memreq);
		return(0);
		}
#else
	m_pASBins = (tsASBin*)mmap(nullptr, (size_t)memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	if (m_pASBins == MAP_FAILED)
		{
		m_pASBins = nullptr;
		ReleaseFastSerialise();
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddASBin: Memory allocation of %zd bytes through mmap()  failed", (int64_t)memreq, strerror(errno));
		return(0);
		}
#endif
	m_AllocdASBinsMem = memreq;
	m_AllocdASBins = cAllocASBins;
	m_UsedASBins = 0;
	}
else
		// needing to allocate more memory?
if ((m_UsedASBins) >= m_AllocdASBins)
	{
	ToAllocdASBins = m_UsedASBins + (int64_t)cReallocASBins;
	size_t memreq = (size_t)ToAllocdASBins * sizeof(tsASBin);
#ifdef _WIN32
	pASBin = (tsASBin*)realloc(m_pASBins, memreq);
	if (pASBin == nullptr)
		{
#else
		pASBin = (tsASBin*)mremap(m_pASBins, m_AllocdASBinsMem, memreq, MREMAP_MAYMOVE);
		if (pASBin == MAP_FAILED)
			{
#endif
			ReleaseFastSerialise();
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddASBin: Memory reallocation to %zd bytes failed - %s", (int64_t)memreq, strerror(errno));
			return(0);
			}
	m_pASBins = pASBin;
	m_AllocdASBinsMem = memreq;
	m_AllocdASBins = ToAllocdASBins;
	}
pASBin = &m_pASBins[m_UsedASBins++];
if (pInitASBin != nullptr)
	*pASBin = *pInitASBin;
else
	memset(pASBin, 0, sizeof(tsASBin));
pASBin->BinID = m_UsedASBins;
ReleaseFastSerialise();
return(m_UsedASBins);
}

tsASBin*    // returned allele scoring bin specification, returns nullptr if all bins have been iterated
CCallHaplotypes::IterateASBins(int32_t PrevBinID,    // previously iterated bin identifier, to start from 1st bin then pass 0 as the bin identifier 
	int32_t ChromID) // only processing bins on this chrom, 0 if processing all bins 
{
tsASBin* pASBin;

AcquireFastSerialise();
while (PrevBinID < m_UsedASBins)
	{
	pASBin = &m_pASBins[PrevBinID++]; // note: bin specs are unordered so need to iterate over all bins so as to ensure all bins on requested chrom are processed 
	if (ChromID != 0 && pASBin->ChromID != ChromID) // needing to match on ChromID ..
		{
		if (ChromID > pASBin->ChromID) // BinIDs ordered by ChromID ascending so checking if iterated past bins for targeted chrom
			break;
		continue;
		}
	// can accept bin spec by ChromID but bin may already have been allocated for processing ...
	if (!(pASBin->ProcState & 0x07)) // bin processing state  - 0x00 if unprocessed, 0x01 if currently being processing, 0x03 if processing completed but referenced chromosome not located, 0x07 if processing successfully completed
		{
		pASBin->ProcState |= 0x01; // bin now allocated for processing by caller
		m_ProccessingASBinID = pASBin->BinID; // marks last bin accepted for processing, will be used to initialise processing for next chromosome
		ReleaseFastSerialise();
		return(pASBin);
		}
	}
ReleaseFastSerialise();
return(nullptr);
}


int
CCallHaplotypes::GroupAlleleScores(char* pszInFile, // file containing allele association scores
	char *pszWGSHapGrpFile,	// WGS haplotype groupings file
	char* pszOutFile) // where to write grouping CSV file
{
int32_t Rslt;
int64_t CurLineNumber;
int64_t ASBinIdx;
int32_t NumFields;

char *pszSrcPBA;
int32_t SrcID;
char *pszRefPBA;
int32_t RefID;

char* pszChrom;
uint32_t CurChromID;
uint32_t PrevChromID;
int32_t BinStartLoci;
int32_t BinSize;
int32_t CentroidDistance;
uint32_t TotNumChromBins;

double ExactScore;
double PartialScore;

if (m_hOutFile != -1)
	{
	close(m_hOutFile);
	m_hOutFile = -1;
	}

if (m_pAlleleScoresCSVFile != nullptr) // should be, but better to be sure!
	{
	delete m_pAlleleScoresCSVFile;
	m_pAlleleScoresCSVFile = nullptr;
	}

if ((m_pAlleleScoresCSVFile = new CCSVFile) == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "GroupAlleleScores: Unable to instantiate CCSVfile");
	Reset();
	return(eBSFerrObj);
	}
m_pAlleleScoresCSVFile->SetMaxFields(14); // only processing 1st 13 fields, allowing 1 extra so can detect if more thasn 13

if ((Rslt = m_pAlleleScoresCSVFile->Open(pszInFile)) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "GroupAlleleScores: Unable to open file: '%s'", pszInFile);
	Reset();
	return(Rslt);
	}

// allocate for tsASBin's if not already allocated
if(m_pASBins == nullptr)
	{
	size_t memreq = (size_t)cAllocASBins * sizeof(tsASBin);
#ifdef _WIN32
	m_pASBins = (tsASBin*)malloc(memreq);	// initial and perhaps the only allocation
	if (m_pASBins == nullptr)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "GroupAlleleScores: Initial memory allocation of %zd bytes for allele score bins failed - %s", (int64_t)memreq, strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
#else
	m_pASBins = (tsASBin*)mmap(nullptr, memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	if (m_pASBins == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "GroupAlleleScores: Memory allocation of %zd bytes through mmap() for allele score bins failed - %s", (int64_t)memreq, strerror(errno));
		m_pASBins = nullptr;
		Reset();
		return(eBSFerrMem);
		}
#endif
	m_AllocdASBinsMem = memreq;
	m_AllocdASBins = cAllocASBins;
	}
m_UsedASBins = 0;
ASBinIdx = 0;
CurLineNumber = 0;
BinStartLoci = 0;
BinSize = 0;
CentroidDistance = 0;
SrcID = 0;
RefID = 0;
PrevChromID = -1;
TotNumChromBins = 0;
// allele association scores file structure
// "SourcePBA","ReferencePBA","Chrom","Bin","BinLoci","BinSize","AlignLen","NumExactMatches","NumBiallelicExactMatches","NumPartialMatches","NumNonRefAlleles","ExactScore","PartialScore"
while ((Rslt = m_pAlleleScoresCSVFile->NextLine()) > 0)		// onto next line containing fields
	{
	CurLineNumber++;
	if ((NumFields = m_pAlleleScoresCSVFile->GetCurFields()) < 13)	// must contain at least 13 fields
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "GroupAlleleScores: Allele association score file '%s' expected to contain a minimum of 13 fields, it contains %d at line %zd", pszInFile, NumFields, CurLineNumber);
		Reset();
		return(eBSFerrParse);
		}

	// 1st row may be a header row, if so then skip this row
	if (CurLineNumber == 1 && m_pAlleleScoresCSVFile->IsLikelyHeaderLine())
		continue;

	Rslt = m_pAlleleScoresCSVFile->GetText(3, &pszChrom);
	if (Rslt < 0 || pszChrom == nullptr || pszChrom[0] == '\0')
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "GroupAlleleScores: Unable to parse chromosome name at line %zd in file '%s'", CurLineNumber, pszInFile);
		Reset();
		return(eBSFerrParse);
		}

	if (!AcceptThisChromName(pszChrom))
		continue;

	Rslt = m_pAlleleScoresCSVFile->GetText(1, &pszSrcPBA);
	if (Rslt < 0 || pszSrcPBA == nullptr || pszSrcPBA[0] == '\0')
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "GroupAlleleScores: Unable to parse source identifier at line %zd in file '%s'", CurLineNumber, pszInFile);
		Reset();
		return(eBSFerrParse);
		}


	if((SrcID = LocateReadset(pszSrcPBA, 1)) == 0)
		SrcID = AddReadset(pszSrcPBA,1);

	Rslt = m_pAlleleScoresCSVFile->GetText(2, &pszRefPBA);
	if (Rslt < 0 || pszRefPBA == nullptr || pszRefPBA[0] == '\0')
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "GroupAlleleScores: Unable to parse reference identifier at line %zd in file '%s'", CurLineNumber, pszInFile);
		Reset();
		return(eBSFerrParse);
		}

	if ((RefID = LocateReadset(pszRefPBA, 0)) == 0)
		RefID = AddReadset(pszRefPBA, 0);


	CurChromID = AddChrom(pszChrom);
	if (CurChromID != PrevChromID)
		{
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "GroupAlleleScores: Loading bin scores for chromosome '%s' at line %zd from file '%s'", pszChrom, CurLineNumber, pszInFile);
		PrevChromID = CurChromID;
		}

	Rslt = m_pAlleleScoresCSVFile->GetInt(5, &BinStartLoci);
	if (Rslt < 0 || BinStartLoci < 0)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "GroupAlleleScores: Error whilst parsing bin loci %d at line %zd in file '%s'", BinStartLoci, CurLineNumber, pszInFile);
		Reset();
		return(eBSFerrParse);
		}
	if(BinStartLoci >= m_ChromSizes[CurChromID -1])
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "GroupAlleleScores: Bin start loci %d not within chromosome '%s' size %d at line %zd in file '%s'", BinStartLoci, m_ChromSizes[CurChromID - 1], pszChrom, CurLineNumber, pszInFile);
		Reset();
		return(eBSFerrParse);
		}

	Rslt = m_pAlleleScoresCSVFile->GetInt(6, &BinSize);
	if (Rslt < 0 || BinSize < 1)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "GroupAlleleScores: Error whilst parsing bin size %d at line %zd in file '%s'", BinSize, CurLineNumber, pszInFile);
		Reset();
		return(eBSFerrParse);
		}

	if (BinStartLoci + BinSize  > m_ChromSizes[CurChromID - 1])
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "GroupAlleleScores: Bin end loci %d not within chromosome '%s' size %d at line %zd in file '%s'", BinStartLoci + BinSize - 1, m_ChromSizes[CurChromID - 1], pszChrom, CurLineNumber, pszInFile);
		Reset();
		return(eBSFerrParse);
		}

	if ((Rslt = m_pAlleleScoresCSVFile->GetDouble(12, &ExactScore)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "GroupAlleleScores: Error while parsing ExactScore at line %zd in file '%s'", ExactScore, CurLineNumber, pszInFile);
		Reset();
		return(eBSFerrParse);
		}

	if (ExactScore < 0.0 || ExactScore > 1.0)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "GroupAlleleScores: ExactScore %f not in range 0..1 at line %zd in file '%s'", ExactScore, CurLineNumber, pszInFile);
		Reset();
		return(eBSFerrParse);
		}

	if ((Rslt = m_pAlleleScoresCSVFile->GetDouble(13, &PartialScore)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "GroupAlleleScores: Error while parsing PartialScore at line %zd in file '%s'", ExactScore, CurLineNumber, pszInFile);
		Reset();
		return(eBSFerrParse);
		}
	if (PartialScore < 0.0 || PartialScore > 1.0)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "GroupAlleleScores: PartialScore %f not in range 0.0 .. 1.0 at line %zd in file '%s'", PartialScore, CurLineNumber, pszInFile);
		Reset();
		return(eBSFerrParse);
		}
	if (PartialScore < ExactScore)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "GroupAlleleScores: PartialScore %f must be at least same as ExactScore %f at line %zd in file '%s'", PartialScore, ExactScore, CurLineNumber, pszInFile);
		Reset();
		return(eBSFerrParse);
		}
	

	if((ASBinIdx = AddASBin(SrcID,RefID,CurChromID, m_ChromSizes[CurChromID - 1], BinStartLoci, BinSize, ExactScore,PartialScore)) < 1)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "GroupAlleleScores: Internal errors processing successfully parsed allele scored bin whilst at line %zd in file '%s'", PartialScore, ExactScore, CurLineNumber, pszInFile);
		Reset();
		return((int32_t)ASBinIdx);
		}
	}
delete m_pAlleleScoresCSVFile;
m_pAlleleScoresCSVFile = nullptr;


	// sorting by ChromID.BinLoci.SrcID ascending, then ScorePartial descending, then RefID ascending
if (m_UsedASBins > 1)
	{
	m_mtqsort.SetMaxThreads(m_NumThreads);
	m_mtqsort.qsort(m_pASBins, m_UsedASBins, sizeof(tsASBin), SortASBin);
	}
else
	if(m_UsedASBins < 1)
		{
		gDiagnostics.DiagOut(eDLWarn, gszProcName, "GroupAlleleScores: Nothing to do, no parsed allele score bins accepted from file '%s'", pszInFile);
		return(eBSFerrParse);
		}

// basic validation of source allele score bins
// Should be same number of bins at each loci as there are source scored samples
// checking for overlapping bins and full chromosome bin coverage
tsASBin* pASBin = m_pASBins;
int32_t CurBinLoci;
int32_t ExpNxtBinLoci;
int32_t CurBinSize;
int32_t CurChromSize;
int32_t NumBinsAtLoci;
int64_t BinIdx;
int32_t ExpBinsAtLoci = m_NumReadsetTypes[0] * m_NumReadsetTypes[1];
int32_t AnError = 0;
CurChromID = 0;
for (BinIdx = 0; BinIdx < m_UsedASBins; BinIdx++, pASBin++)
	{
	if(pASBin->ChromID != CurChromID) // first bin this chrom
		{
		if(CurChromID != 0)
			{
			if(NumBinsAtLoci != ExpBinsAtLoci)
				break;
			if(ExpNxtBinLoci != CurChromSize)
				break;
			}
		CurChromID = pASBin->ChromID;
		CurBinSize = pASBin->BinSize;
		ExpNxtBinLoci = pASBin->BinSize;
		CurChromSize = pASBin->ChromSize;
		CurBinLoci = 0;
		NumBinsAtLoci = 1;
		continue;
		}
	if(pASBin->BinLoci != CurBinLoci)	// new bin loci on same chrom
		{
		if(CurChromSize != pASBin->ChromSize)				// same chromosome so expecting to be same size except if last bin on the chromosome
			break;
		if(NumBinsAtLoci != ExpBinsAtLoci)			// expecting same numbers of bins at each loci
			break;
		if(ExpNxtBinLoci != pASBin->BinLoci)				// loci must be same as expected
			break;
		CurBinLoci = pASBin->BinLoci;
		CurBinSize = pASBin->BinSize;
		ExpNxtBinLoci += CurBinSize;
		NumBinsAtLoci = 1;
		continue;
		}
	NumBinsAtLoci++;
	if (CurChromSize != pASBin->ChromSize)				// same chromosome so expecting to be same size
		break;
	if (NumBinsAtLoci > ExpBinsAtLoci)			// expecting same numbers of bins at each loci
		break;
	if(CurBinSize != pASBin->BinSize)
		break;
	}

if (BinIdx < m_UsedASBins || NumBinsAtLoci != ExpBinsAtLoci || ExpNxtBinLoci != CurChromSize)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "GroupAlleleScores: Inconsistencies in file '%s'", pszWGSHapGrpFile);
	return(eBSFerrParse);
	}

// 
// load the WGS haplotype groupings file
if((Rslt = LoadRefGroupings(pszWGSHapGrpFile)) < eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "GroupAlleleScores: Errors parsing WGS reference groupings from file '%s'", pszWGSHapGrpFile);
	return(Rslt);
	}

// checking that WGS grouping bins and allele score bins have same start/sizes
ExpBinsAtLoci = m_NumReadsetTypes[0] * m_NumReadsetTypes[1];
pASBin = m_pASBins;
tsASRefGrp *pASRefGrpBin = m_pASRefGrps;
for (BinIdx = 0; BinIdx < m_UsedASBins; BinIdx += ExpBinsAtLoci, pASBin+= ExpBinsAtLoci, pASRefGrpBin = (tsASRefGrp *)((uint8_t *)pASRefGrpBin + pASRefGrpBin->Size))
	{
	if(pASBin->ChromID != pASRefGrpBin->ChromID || pASBin->BinLoci != pASRefGrpBin->StartLoci || pASBin->BinSize != pASRefGrpBin->NumLoci)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "GroupAlleleScores: Mismatch between WGS grouping bins and Allele scoring bins");
		return(eBSFerrParse);
		}
	}

Rslt = ProcessAlleleScoreBins(pszOutFile);

return(Rslt);
}


int32_t
CCallHaplotypes::ProcessAlleleScoreBins(char* pszOutFile) // where to write allele score bins grouping CSV file
{
int32_t Rslt;

if(m_UsedASBins < 1)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessAlleleScoredBins: Nothing to do, no allele scored bins available for processing");
	Reset();
	return(eBSFerrInternal);
	}

if (m_hOutFile != -1)
	{
	close(m_hOutFile);
	m_hOutFile = -1;
	}
if (m_pszOutBuffer == nullptr)
	{
	if ((m_pszOutBuffer = new uint8_t[cOutBuffSize]) == nullptr)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessAlleleScoredBins: Memory allocation of %u bytes for grouping output buffering failed - %s", cOutBuffSize, strerror(errno));
		Reset();
		return(eBSFerrMem);
		};
	m_AllocOutBuff = cOutBuffSize;
	}
m_OutBuffIdx = 0;

#ifdef _WIN32
m_hOutFile = open(pszOutFile, (O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC), (_S_IREAD | _S_IWRITE));
#else
if ((m_hOutFile = open64(pszOutFile, O_WRONLY | O_CREAT, S_IREAD | S_IWRITE)) != -1)
	if (ftruncate(m_hOutFile, 0) != 0)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessAlleleScoredBins: Unable to create/truncate %s - %s", pszOutFile, strerror(errno));
		Reset();
		return(eBSFerrCreateFile);
		}
#endif
if (m_hOutFile < 0)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessAlleleScoredBins: Unable to create/truncate %s - %s", pszOutFile, strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}
Rslt = eBSFSuccess; // a big assumption!!!

// iterate over bins and retain the highest scoring for each ChromID.BinLoci.SrcID
tsASBin *pASBin = m_pASBins;
tsASBin* pASBinMH = nullptr;
int32_t* pRefIDGrpMbrs;
int32_t RefIDidx;
int32_t GrpIdx;
int RefIdx;
int64_t BinIdx;
int CurChromID = 0;
int CurBinLoci = 0;
int CurBinSize = 0;
int32_t CurSrcID = 0;
double CurScorePartial = 0.0;
double CurScoreExact = 0.0;
for(BinIdx = 0; BinIdx < m_UsedASBins; BinIdx++, pASBin++)
	{
	if (pASBin->ChromID == CurChromID && pASBin->BinLoci == CurBinLoci && pASBin->SrcID == CurSrcID)	// if same Chrom.BinLoci.SrcID
		{
		if (pASBin->ScorePartial < CurScorePartial || (pASBin->ScorePartial == CurScorePartial && pASBin->ScoreExact <= CurScoreExact)) // using ScoreExact as the tiebreaker
			{
			pASBin->ProcState &= ~0x40;
			continue;
			}
		// have a higher scoring bin than previously highest scorer
		pASBinMH->ProcState &= ~0x40;
		pASBin->ProcState |= 0x40;
		CurScorePartial = pASBin->ScorePartial;
		CurScoreExact = pASBin->ScoreExact;
		pASBinMH = pASBin;
		continue;
		}
	// only reaching here if starting a new Chrom.BinLoci.CurSrcID
	// until otherwise proven this will be the highest scoring
	CurChromID = pASBin->ChromID;
	CurBinLoci = pASBin->BinLoci;
	CurSrcID = pASBin->SrcID;
	CurScorePartial = pASBin->ScorePartial;
	CurScoreExact = pASBin->ScoreExact;
	pASBin->ProcState |= 0x40;					// flag that this is the highest scoring bin
	pASBinMH = pASBin;
	}


// output header
m_OutBuffIdx = sprintf((char*)&m_pszOutBuffer[0], "\"Chrom\",\"BinLoci\",\"BinSize\"");
for (RefIdx = 0; RefIdx < m_NumReadsetTypes[0]; RefIdx++)
	m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], ",\"%s\"", LocateReadset(m_FndrIDMappings[RefIdx]));
if (!CUtility::RetryWrites(m_hOutFile, m_pszOutBuffer, m_OutBuffIdx))
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessAlleleScoreBins: Fatal error in RetryWrites()");
	Reset();
	return(eBSFerrFileAccess);
	}
m_OutBuffIdx = 0;


// reporting the number of reference instances for each highest scoring src bins at ChromID.BinLoci
RefIDidx = 0;
GrpIdx = 0;
RefIdx = 0;

pRefIDGrpMbrs = new int32_t[m_NumReadsetTypes[0]];
memset(pRefIDGrpMbrs,0,sizeof(int32_t) * m_NumReadsetTypes[0]);
pASBin = m_pASBins;
CurBinLoci = 0;
CurChromID = 0;
CurBinSize = 0;
for (BinIdx = 0; BinIdx < m_UsedASBins; BinIdx++, pASBin++)
	{
	if (CurChromID != pASBin->ChromID || CurBinLoci != pASBin->BinLoci)	// starting a new chromosome or bin loci changed
		{
		if (CurChromID > 0)
			{
			m_OutBuffIdx += sprintf((char *) & m_pszOutBuffer[m_OutBuffIdx], "\n\"%s\",%d,%d", LocateChrom(CurChromID), CurBinLoci,CurBinSize);
			for(RefIdx = 0; RefIdx < m_NumReadsetTypes[0]; RefIdx++)
				{
				m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx],",%d", pRefIDGrpMbrs[RefIdx]);

				if(m_OutBuffIdx + 1000 > m_AllocOutBuff)
					{
					if (!CUtility::RetryWrites(m_hOutFile, m_pszOutBuffer, m_OutBuffIdx))
						{
						gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessAlleleScoreBins: Fatal error in RetryWrites()");
						if (pRefIDGrpMbrs != nullptr)
							delete[]pRefIDGrpMbrs;
						Reset();
						return(eBSFerrFileAccess);
						}
					m_OutBuffIdx = 0;
					}
				}
			}
		memset(pRefIDGrpMbrs, 0, sizeof(int32_t) * m_NumReadsetTypes[0]);
		CurChromID = pASBin->ChromID;
		CurBinLoci = pASBin->BinLoci;
		CurBinSize = pASBin->BinSize;
		}
	if (pASBin->ProcState & 0x40)
		{
		pRefIDGrpMbrs[LocateReadsetIDIdx(pASBin->RefID,0)] += 1;
		}

	}
if (CurChromID > 0)
	{
	m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], "\n\"%s\",%d,%d", LocateChrom(CurChromID), CurBinLoci,CurBinSize);
	for (RefIdx = 0; RefIdx < m_NumReadsetTypes[0]; RefIdx++)
		{
		m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], ",%d", pRefIDGrpMbrs[RefIdx]);
		if (m_OutBuffIdx + 1000 > m_AllocOutBuff)
			{
			if (!CUtility::RetryWrites(m_hOutFile, m_pszOutBuffer, m_OutBuffIdx))
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessAlleleScoreBins: Fatal error in RetryWrites()");
				if (pRefIDGrpMbrs != nullptr)
					delete[]pRefIDGrpMbrs;
				Reset();
				return(eBSFerrFileAccess);
				}
			m_OutBuffIdx = 0;
			}
		}
	}
if (m_OutBuffIdx > 0 && m_hOutFile != -1)
	{
	if (!CUtility::RetryWrites(m_hOutFile, m_pszOutBuffer, m_OutBuffIdx))
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessAlleleScoreBins: Fatal error in RetryWrites()");
		if (pRefIDGrpMbrs != nullptr)
			delete[]pRefIDGrpMbrs;
		Reset();
		return(eBSFerrFileAccess);
		}
	m_OutBuffIdx = 0;
	}
if(pRefIDGrpMbrs != nullptr)
	delete []pRefIDGrpMbrs;
if (m_hOutFile != -1)
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

// report the actual highest scored reference for each ChromID.BinLoci.SrcID. Highest scored bin will have ProcState bit 5 (0x40) set
char szSrcRefs[_MAX_PATH];
CUtility::AppendFileNameSuffix(szSrcRefs, pszOutFile,(char *)".srcrefs.csv",'.');

#ifdef _WIN32
m_hOutFile = open(szSrcRefs, (O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC), (_S_IREAD | _S_IWRITE));
#else
if ((m_hOutFile = open64(szSrcRefs, O_WRONLY | O_CREAT, S_IREAD | S_IWRITE)) != -1)
if (ftruncate(m_hOutFile, 0) != 0)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessAlleleScoredBins: Unable to create/truncate %s - %s", szSrcRefs, strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}
#endif
if (m_hOutFile < 0)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessAlleleScoredBins: Unable to create/truncate %s - %s", szSrcRefs, strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}
Rslt = eBSFSuccess; // a big assumption!!!


// output header
pASBin = m_pASBins;
m_OutBuffIdx = sprintf((char*)&m_pszOutBuffer[0], "\"Chrom\",\"BinLoci\",\"BinSize\"");
for (BinIdx = 0; BinIdx < m_NumReadsetTypes[1]; BinIdx++, pASBin+= m_NumReadsetTypes[0])
	m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], ",\"%s\"", LocateReadset(pASBin->SrcID));
if (!CUtility::RetryWrites(m_hOutFile, m_pszOutBuffer, m_OutBuffIdx))
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessAlleleScoreBins: Fatal error in RetryWrites()");
	Reset();
	return(eBSFerrFileAccess);
	}
m_OutBuffIdx = 0;
CurBinLoci = 0;
CurChromID = 0;
CurBinSize = 0;
pASBin = m_pASBins;
for (BinIdx = 0; BinIdx < m_UsedASBins; BinIdx++, pASBin++)
	{
	if (!(pASBin->ProcState & 0x40))   // only interested in highest scoring bins
		continue;
	if (CurChromID != pASBin->ChromID)
		{
		CurChromID = pASBin->ChromID;
		CurBinLoci = pASBin->BinLoci;
		CurBinSize = pASBin->BinSize;
		m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], "\n\"%s\",%d,%d", LocateChrom(CurChromID), CurBinLoci, CurBinSize);
		}
	else
		if (CurBinLoci != pASBin->BinLoci)	// same chromosome but bin loci changed
			{
			CurBinLoci = pASBin->BinLoci;
			CurBinSize = pASBin->BinSize;
			m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], "\n\"%s\",%d,%d", LocateChrom(CurChromID), CurBinLoci, CurBinSize);
			}
	m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], ",\"%s\"",LocateReadset(pASBin->RefID));
	if (m_OutBuffIdx + 1000 > m_AllocOutBuff)
		{
		if (!CUtility::RetryWrites(m_hOutFile, m_pszOutBuffer, m_OutBuffIdx))
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessAlleleScoreBins: Fatal error in RetryWrites()");
			Reset();
			return(eBSFerrFileAccess);
			}
		}
	}

if (m_OutBuffIdx > 0 && m_hOutFile != -1)
	{
	if (!CUtility::RetryWrites(m_hOutFile, m_pszOutBuffer, m_OutBuffIdx))
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessAlleleScoreBins: Fatal error in RetryWrites()");
		if (pRefIDGrpMbrs != nullptr)
			delete[]pRefIDGrpMbrs;
		Reset();
		return(eBSFerrFileAccess);
		}
	m_OutBuffIdx = 0;
	}
if (m_hOutFile != -1)
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

// generate group association mappings - sources are associated to same groups as the references that sources scored to 
CUtility::AppendFileNameSuffix(szSrcRefs, pszOutFile, (char*)".srcrefs.grps.csv", '.');

#ifdef _WIN32
m_hOutFile = open(szSrcRefs, (O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC), (_S_IREAD | _S_IWRITE));
#else
if ((m_hOutFile = open64(szSrcRefs, O_WRONLY | O_CREAT, S_IREAD | S_IWRITE)) != -1)
if (ftruncate(m_hOutFile, 0) != 0)
{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessAlleleScoredBins: Unable to create/truncate %s - %s", szSrcRefs, strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
}
#endif
if (m_hOutFile < 0)
{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessAlleleScoredBins: Unable to create/truncate %s - %s", szSrcRefs, strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
}
Rslt = eBSFSuccess; // a big assumption!!!


// output header
pASBin = m_pASBins;
m_OutBuffIdx = sprintf((char*)&m_pszOutBuffer[0], "\"ExprID\",\"RowID\",\"Chrom\",\"Loci\",\"Len\",\"MinDistance\",\"MaxDistance\",\"MaxHaplotypes\",\"ActualDistance\",\"ActualHaplotypes\"");
for (BinIdx = 0; BinIdx < m_NumReadsetTypes[1]; BinIdx++, pASBin+= m_NumReadsetTypes[0])
	m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], ",\"%s\"", LocateReadset(pASBin->SrcID));
m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], ",\"GrpMembers:1\",\"GrpMembers:2\",\"GrpMembers:3\",\"GrpMembers:4\",\"GrpMembers:5\"");

if (!CUtility::RetryWrites(m_hOutFile, m_pszOutBuffer, m_OutBuffIdx))
{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessAlleleScoreBins: Fatal error in RetryWrites()");
	Reset();
	return(eBSFerrFileAccess);
}
m_OutBuffIdx = 0;
CurBinLoci = 0;
CurChromID = 0;
CurBinSize = 0;
pASBin = m_pASBins;
int32_t CurRefGrp;
int32_t GrpInstDist[100];

tsASRefGrp* pMRAASRefGrp = nullptr;
tsASRefIDGrp* pASRefIDGrp = nullptr;
for (BinIdx = 0; BinIdx < m_UsedASBins; BinIdx++, pASBin++)
	{
	if (!(pASBin->ProcState & 0x40))   // only interested in highest scoring bins
		continue;
	if (CurChromID != pASBin->ChromID || CurBinLoci != pASBin->BinLoci || CurBinSize != pASBin->BinSize)
		{
		if(CurChromID != 0)
			for(GrpIdx = 0; GrpIdx < 5; GrpIdx++)
				m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx],",%d",GrpInstDist[GrpIdx]);
		memset(GrpInstDist,0,sizeof(GrpInstDist));
		CurChromID = pASBin->ChromID;
		CurBinLoci = pASBin->BinLoci;
		CurBinSize = pASBin->BinSize;
		pMRAASRefGrp = LocateASRefGrp(CurChromID, CurBinLoci);
		m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], "\n1,%d,\"%s\",%d,%d", pMRAASRefGrp->BinID,LocateChrom(CurChromID), CurBinLoci, CurBinSize);
		m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], ",%d,%d,%d,%d,%d", pMRAASRefGrp->MinCentroidDistance, pMRAASRefGrp->MaxCentroidDistance, pMRAASRefGrp->MaxHaplotypeGroups, pMRAASRefGrp->ActualCentroidDistance, pMRAASRefGrp->ActualHaplotypeGroups);
		}

	if((pASRefIDGrp = LocateASRefIDGrp(CurChromID, CurBinLoci, pASBin->RefID))==nullptr)
		CurRefGrp = 1;
	else
		CurRefGrp = pASRefIDGrp->Grp;
	GrpInstDist[CurRefGrp-1]++;

	m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], ",%d", CurRefGrp);
	if (m_OutBuffIdx + 1000 > m_AllocOutBuff)
		{
		if (!CUtility::RetryWrites(m_hOutFile, m_pszOutBuffer, m_OutBuffIdx))
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessAlleleScoreBins: Fatal error in RetryWrites()");
			Reset();
			return(eBSFerrFileAccess);
			}
		}
	}
if (CurChromID != 0)
	for (GrpIdx = 0; GrpIdx < 5; GrpIdx++)
		m_OutBuffIdx += sprintf((char*)&m_pszOutBuffer[m_OutBuffIdx], ",%d", GrpInstDist[GrpIdx]);
if (m_OutBuffIdx > 0 && m_hOutFile != -1)
	{
	if (!CUtility::RetryWrites(m_hOutFile, m_pszOutBuffer, m_OutBuffIdx))
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessAlleleScoreBins: Fatal error in RetryWrites()");
		if (pRefIDGrpMbrs != nullptr)
			delete[]pRefIDGrpMbrs;
		Reset();
		return(eBSFerrFileAccess);
		}
	m_OutBuffIdx = 0;
	}
if (m_hOutFile != -1)
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
return(Rslt);
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
CCallHaplotypes::AccumWIGCnts(int32_t ChromID,	// accumulate wiggle counts into variableStep spans over this chromosome
			int32_t Loci,		// this loci - starts from 1 not 0!
	        int32_t BinLen,     // span bin is this length
			uint32_t Cnts,		 // has this many counts attributed to the span
			int32_t MaxSpanLen) // allow WIG spans to be this maximal length
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

// sorting by ChromID.Loci ascending
int
CCallHaplotypes::SortASRefGrp(const void* arg1, const void* arg2)
{
	tsASRefGrp* pEl1 = (tsASRefGrp*)arg1;
	tsASRefGrp* pEl2 = (tsASRefGrp*)arg2;
	if (pEl1->ChromID < pEl2->ChromID)					// sort ascending
		return(-1);
	if (pEl1->ChromID > pEl2->ChromID)
		return(1);
	if (pEl1->StartLoci < pEl2->StartLoci)
		return(-1);
	if (pEl1->StartLoci > pEl2->StartLoci)
		return(1);
	return(0);
}


// sorting by ChromID.BinLoci.SrcID ascending, then ScorePartial descending, ScoreExact descending, then RefID ascending
int
CCallHaplotypes::SortASBin(const void* arg1, const void* arg2)
{
	tsASBin* pEl1 = (tsASBin*)arg1;
	tsASBin* pEl2 = (tsASBin*)arg2;
	if (pEl1->ChromID < pEl2->ChromID)					// sort ascending
		return(-1);
	if (pEl1->ChromID > pEl2->ChromID)
		return(1);
	if (pEl1->BinLoci < pEl2->BinLoci)
		return(-1);
	if (pEl1->BinLoci > pEl2->BinLoci)
		return(1);
	if (pEl1->SrcID < pEl2->SrcID)
		return(-1);
	if (pEl1->SrcID > pEl2->SrcID)
		return(1);

	if (pEl1->ScorePartial < pEl2->ScorePartial)		// sort descending
		return(1);
	if (pEl1->ScorePartial > pEl2->ScorePartial)
		return(-1);
	if (pEl1->ScoreExact < pEl2->ScoreExact)			// sort descending
		return(1);
	if (pEl1->ScoreExact > pEl2->ScoreExact)
		return(-1);

	if (pEl1->RefID < pEl2->RefID)						// back to sort ascending
		return(-1);
	if (pEl1->RefID > pEl2->RefID)
		return(1);

	return(0);
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

// sorting by DGTLoci by FMeasure descending
// Note that at any loci there could be - but unlikely - upto 4 F-measures, one for each allele
int
CCallHaplotypes::SortDGTLociFMeasureDesc(const void* arg1, const void* arg2)
{
tsDGTLoci* pEl1 = (tsDGTLoci*)arg1;
tsDGTLoci* pEl2 = (tsDGTLoci*)arg2;
if(pEl1->AlleleFMeasuresDesc[0] > pEl2->AlleleFMeasuresDesc[0])
	return(-1);
if(pEl1->AlleleFMeasuresDesc[0] < pEl2->AlleleFMeasuresDesc[0])
	return(1);
if(pEl1->AlleleFMeasuresDesc[1] < 0.1 && pEl2->AlleleFMeasuresDesc[1] < 0.1) // assumes any Fmeasure < 0.1 represents non-acceptance
	return(0);
if(pEl1->AlleleFMeasuresDesc[1] > pEl2->AlleleFMeasuresDesc[1])
	return(-1);
if(pEl1->AlleleFMeasuresDesc[1] < pEl2->AlleleFMeasuresDesc[1])
	return(1);
if(pEl1->AlleleFMeasuresDesc[2] < 0.1 && pEl2->AlleleFMeasuresDesc[2] < 0.1)
	return(0);
if(pEl1->AlleleFMeasuresDesc[2] > pEl2->AlleleFMeasuresDesc[2])
	return(-1);
if(pEl1->AlleleFMeasuresDesc[2] < pEl2->AlleleFMeasuresDesc[2])
	return(1);
if(pEl1->AlleleFMeasuresDesc[3] < 0.1 && pEl2->AlleleFMeasuresDesc[3] < 0.1)
	return(0);
if(pEl1->AlleleFMeasuresDesc[3] > pEl2->AlleleFMeasuresDesc[3])
	return(-1);
if(pEl1->AlleleFMeasuresDesc[3] < pEl2->AlleleFMeasuresDesc[3])
	return(1);
return(0);
}

// sorting by ChromID.DGTLoci ascending
int
CCallHaplotypes::SortDGTLoci(const void* arg1, const void* arg2)
{
tsDGTLoci* pEl1 = (tsDGTLoci*)arg1;
tsDGTLoci* pEl2 = (tsDGTLoci*)arg2;
if(pEl1->ChromID < pEl2->ChromID)
	return(-1);
if(pEl1->ChromID > pEl2->ChromID)
	return(1);
if(pEl1->DGTLoci < pEl2->DGTLoci)
	return(-1);
if(pEl1->DGTLoci > pEl2->DGTLoci)
	return(1);
return(0);
}

// sorting by ChromID.KMerLoci ascending
int
CCallHaplotypes::SortKMerLoci(const void* arg1, const void* arg2)
{
	tsKMerLoci* pEl1 = (tsKMerLoci*)arg1;
	tsKMerLoci* pEl2 = (tsKMerLoci*)arg2;
	if (pEl1->ChromID < pEl2->ChromID)
		return(-1);
	if (pEl1->ChromID > pEl2->ChromID)
		return(1);
	if (pEl1->CurLoci < pEl2->CurLoci)
		return(-1);
	if (pEl1->CurLoci > pEl2->CurLoci)
		return(1);
	return(0);
}

