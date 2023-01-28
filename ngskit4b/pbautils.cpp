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
#include "./pbautils.h"

int Process(ePBAuMode PMode,	// processing mode: ePBAu2Fasta PBA to Fasta, ePBAu2PBA Fasta to PBA, ePBAu2PBAConcordance concordance over PBA samples, ePBAu2WIGConcordance concordance over WIG samples, ePBAu2AVCF allelic variant VCF, ePBAu2GVCF allelic genotype VCF,ePBAu2DVCF deletion genotype VCF,ePBAu2BED generate BED containing all transcribed regions
	int32_t LimitPBAs,			// limit number of loaded PBA files to this many. 1 .. cMaxPBAReadsets
	int32_t PBAsTrim5,			// trim this many aligned PBAs from 5' end of aligned segments - reduces false alleles due to sequencing errors
	int32_t PBAsTrim3,			// trim this many aligned PBAs from 3' end of aligned segments - reduces false alleles due to sequencing errors
	double GTPropNAThres,		// when genotyping VCF then proportion of non-aligned calls over all samples must be < this threshold
	double GTPropHetThres,		// when genotyping VCF then proportion of heterozygous calls must be >= this threshold
	char* pszRefAssemb,			// reference assembly identifier
	char* pszRefAssembFile,		// PBA file containing reference assembly sequences
	char* pszChromFile,			// BED file containing reference assembly chromosome names and sizes
	char* pszSeqID,				// sequence identifier
	char* pszExprID,			// experiment identifier
	int NumInputFiles,			// number of input founder file specs
	char *pszInputFiles[],		// names of input founder PBA files (wildcards allowed)
	char *szOutFile,			// output to this file
	char *pszROIFile,			// optional BED file containing regions of interest
	char *pszGTSampleFiltFile,	// optional file containing those readseq/sample names for whom only their genotypes are to be reported
	int NumIncludeChroms,		// number of chromosome regular expressions to include
	char *pszIncludeChroms[],	// array of include chromosome regular expressions
	int NumExcludeChroms,		// number of chromosome expressions to exclude
	char *pszExcludeChroms[],	// array of exclude chromosome regular expressions
	int NumThreads);			// number of worker threads to use

#ifdef _WIN32
int pbautils(int argc, char* argv[])
{
	// determine my process name
_splitpath(argv[0], nullptr, nullptr, gszProcName, nullptr);
#else
int
pbautils(int argc, char** argv)
{
	// determine my process name
CUtility::splitpath((char*)argv[0], nullptr, gszProcName);
#endif
int iFileLogLevel;			// level of file diagnostics
int iScreenLogLevel;		// level of file diagnostics
char szLogFile[_MAX_PATH];	// write diagnostics to this file
char szChromFile[_MAX_PATH];	// BED file containing chromosome names and sizes
char szRefAssembFile[_MAX_PATH];	// PBA file containing reference assembly sequences

int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure
int NumberOfProcessors;		// number of installed CPUs
int NumThreads;				// number of threads (0 defaults to number of CPUs or a maximum of cMaxPBAWorkerThreads)

int Idx;

ePBAuMode PMode;			// processing mode: ePBAu2Fasta PBA to Fasta, ePBAu2PBA Fasta to PBA, ePBAu2PBAConcordance concordance over PBA samples, ePBAu2WIGConcordance concordance over WIG samples, ePBAu2AVCF allelic variant VCF, ePBAu2GVCF allelic genotype VCF, ePBAu2DVCF deletion genotype VC,FePBAu2BED generate BED containing all transcribed regions
int32_t LimitPBAs;			// limit number of loaded PBA files to this many. 0: no limits, > 0 sets upper limit
int32_t PBAsTrim5;			// trim this many aligned PBAs from 5' end of aligned segments - reduces false alleles due to sequencing errors
int32_t PBAsTrim3;			// trim this many aligned PBAs from 3' end of aligned segments - reduces false alleles due to sequencing errors
double GTPropNAThres;		// when genotyping VCF then proportion of non-aligned calls over all samples must be < this threshold
double GTPropHetThres;		// when genotyping VCF then proportion of heterozygous calls must be >= this threshold
char szExprID[cMaxRefAssembName + 1];
char szSeqID[cMaxRefAssembName+1];
char szRefAssemb[cMaxRefAssembName+1]; 
int NumIncludeChroms;
char* pszIncludeChroms[cMaxIncludeChroms];
int NumExcludeChroms;
char* pszExcludeChroms[cMaxExcludeChroms];

int NumInputFiles;		// number of input files
char* pszInputFiles[cMaxWildCardFileSpecs];		// names of input files (wildcards allowed)
char szOutFile[_MAX_PATH];		// output converted to this file
char szROIFile[_MAX_PATH];		// optional BED file containing regions of interest
char szGTSampleFiltFile[_MAX_PATH];	// optional file containing those readseq/sample names for whom only their genotypes are to be reported
struct arg_lit* help = arg_lit0("h", "help", "print this help and exit");
struct arg_lit* version = arg_lit0("v", "version,ver", "print version information and exit");
struct arg_int* FileLogLevel = arg_int0("f", "FileLogLevel", "<int>", "Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file* LogFile = arg_file0("F", "log", "<file>", "diagnostics log file");

struct arg_int* pmode = arg_int0("m", "mode", "<int>", "processing mode: 0 PBA to Fasta, 1 Fasta to PBA, 2 concordance over PBA samples, 3 concordance over WIG samples, 4 allelic variant VCF, 5 allelic genotype VCF, 6 deletion genotype VCF, 7 segment BED");
struct arg_int* limitpbas = arg_int0("l", "limitpbas", "<int>", " limit number of loaded PBA files to this many. 0: no limits, > 0 sets upper limit (default 0)");
struct arg_int* pbastrim5 = arg_int0("x", "trim5", "<int>", "trim this many aligned PBAs from 5' end of aligned segments (default 0, range 0..100)");
struct arg_int* pbastrim3 = arg_int0("X", "trim3", "<int>", "trim this many aligned PBAs from 3' end of aligned segments (default 0, range 0..100)");

struct arg_dbl* mingthet = arg_dbl0("y", "mingthet", "<dbl>", "(mode 5,6) minimum proportion required genotyping aligned sample heterozygosity or (mode 7) min proportion of samples for segment initiation");
struct arg_dbl* maxgtpropna = arg_dbl0("Y", "maxgtpropna", "<dbl>", "(mode 5,6) maximum proportion allowed genotyping non-aligned samples  or (mode 7) min proportion of samples for segment extension");


struct arg_str* exprid = arg_str0("e", "exprid", "<str>", "assign this experiment identifier, mode 1 only (default is 'E000000')");
struct arg_str* seqid = arg_str0("s", "seqid", "<str>", "sequence identifier, mode 1 only (default is 'S000000')");
struct arg_str* refassemb = arg_str0("r", "refassemb", "<str>", "reference assembly, mode 1 only (default is 'Wm82v4')");
struct arg_str* ExcludeChroms = arg_strn("Z", "chromexclude", "<string>", 0, cMaxExcludeChroms, "high priority - regular expressions defining chromosomes to exclude");
struct arg_str* IncludeChroms = arg_strn("z", "chromeinclude", "<string>", 0, cMaxIncludeChroms, "low priority - regular expressions defining chromosomes to include");
struct arg_file* infiles = arg_filen("i", "infiles", "<file>", 0, cMaxWildCardFileSpecs, "input file(s), wildcards allowed, limit of 200 filespecs supported");
struct arg_file* chromfile = arg_file1("c", "chromfile", "<file>", "input BED file containing chromosome names and sizes");

struct arg_file* refassembfile = arg_file0("R", "refassembfile", "<file>", "reference PBA file, required when generating VCF from PBA files");
struct arg_file* roifile = arg_file0("C", "roifile", "<file>", "BED file containing regions of interest - optional - when generating VCFs from PBA files");
struct arg_file* gtfiltsamplefile = arg_file0("g", "gtfiltsampleile", "<file>", "text file containing names (spaces,quotes,or line separated) of readsets/samples - optional - to actually report in GT VCF file");
struct arg_file* outfile = arg_file1("o", "out", "<file>", "output to this file");
struct arg_int* threads = arg_int0("T", "threads", "<int>", "number of processing threads 0..8 (defaults to 0 which limits threads to maximum of 8 CPU cores)");

struct arg_end* end = arg_end(200);
void* argtable[] = {help,version,FileLogLevel,LogFile,
					pmode,limitpbas, pbastrim5,pbastrim3,mingthet,maxgtpropna, refassemb,exprid,seqid,chromfile,refassembfile,
					ExcludeChroms,IncludeChroms,infiles,outfile,roifile,gtfiltsamplefile,threads,end};
char** pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc, (char**)argv, &pAllArgs);
if(argerrors >= 0)
argerrors = arg_parse(argerrors, pAllArgs, argtable);

/* special case: '--help' takes precedence over error reporting */
if(help->count > 0)
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
if(version->count > 0)
	{
	printf("\n%s %s Version %s\n", gszProcName, gpszSubProcess->pszName, kit4bversion);
	return(1);
	}

if(!argerrors)
	{
	if(FileLogLevel->count && !LogFile->count)
		{
		printf("\nError: FileLogLevel '-f%d' specified but no logfile '-F<logfile>\n'", FileLogLevel->ival[0]);
		exit(1);
		}

	iScreenLogLevel = iFileLogLevel = FileLogLevel->count ? FileLogLevel->ival[0] : eDLInfo;
	if(iFileLogLevel < eDLNone || iFileLogLevel > eDLDebug)
		{
		printf("\nError: FileLogLevel '-l%d' specified outside of range %d..%d\n", iFileLogLevel, eDLNone, eDLDebug);
		exit(1);
		}

	if(LogFile->count)
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

	PMode = pmode->count ? (ePBAuMode)pmode->ival[0] : ePBAu2Fasta;
	if(PMode < ePBAu2Fasta || PMode >= ePBAuPlaceholder)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Unsupported processing mode '-m%d'\n", PMode);
		exit(1);
		}
	if (PMode == ePBAu2PBAConcordance || PMode == ePBAu2BED || PMode == ePBAu2WIGConcordance)
		{
		LimitPBAs = limitpbas->count ? limitpbas->ival[0] : cMaxPBAReadsets;
		if (LimitPBAs < 1 || LimitPBAs > cMaxPBAReadsets)
			LimitPBAs = cMaxPBAReadsets;
		}
	else
		LimitPBAs = 1;

	PBAsTrim5 = pbastrim5->count ? pbastrim5->ival[0] : 0;
	if (PBAsTrim5 < 0)			// silently force to be in range 0.100
		PBAsTrim5 = 0;
	else
		if (PBAsTrim5 > 100)
			PBAsTrim5 = 100;

	PBAsTrim3 = pbastrim3->count ? pbastrim3->ival[0] : 0;
	if (PBAsTrim3 < 0)			// silently force to be in range 0.100
		PBAsTrim3 = 0;
	else
		if (PBAsTrim3 > 100)
			PBAsTrim3 = 100;

	GTPropNAThres = 0.0;
	GTPropHetThres = 0.0;
	if (PMode == ePBAu2BED)
		{
		GTPropNAThres = maxgtpropna->count ? maxgtpropna->dval[0] : cDfltGTPropThres; // repurposing for minimum proportion aligned samples segment initiation 
		if(GTPropNAThres < 0.001)	// silently force to be in accepted range
			GTPropNAThres = 0.001;
		else
			if(GTPropNAThres > 0.50)
				GTPropNAThres = 0.50;
			
		GTPropHetThres = mingthet->count ? mingthet->dval[0] : GTPropNAThres/2.0; // repurposing for minimum proportion aligned samples segment continuation
		if(GTPropHetThres < 0.001)	// silently force to be in accepted range
			GTPropHetThres = 0.001;
		else
			if(GTPropHetThres > GTPropNAThres)
				GTPropHetThres = GTPropNAThres;
		}
	else
		if(PMode == ePBAu2GVCF || PMode == ePBAu2DVCF)
			{
			GTPropNAThres = maxgtpropna->count ? maxgtpropna->dval[0] : cDfltGTPropThres;
			if(GTPropNAThres < 0.00)	// silently force to be in accepted range
				GTPropNAThres = 0.00;
			else
				if(GTPropNAThres > 0.25)
					GTPropNAThres = 0.25;
			
			GTPropHetThres = mingthet->count ? mingthet->dval[0] : cDfltGTPropThres;
			if(GTPropHetThres < 0.00)	// silently force to be in accepted range
				GTPropHetThres = 0.00;
			else
				if(GTPropHetThres > 0.25)
					GTPropHetThres = 0.25;
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
	int MaxAllowedThreads = min((int)cMaxPBAutilityThreads, NumberOfProcessors);	// limit to be at most cMaxWorkerThreads
	if ((NumThreads = threads->count ? threads->ival[0] : MaxAllowedThreads) == 0)
		NumThreads = MaxAllowedThreads;
	if (NumThreads < 0 || NumThreads > MaxAllowedThreads)
		{
		gDiagnostics.DiagOut(eDLWarn, gszProcName, "Warning: Number of threads '-T%d' specified was outside of range %d..%d", NumThreads, 1, MaxAllowedThreads);
		gDiagnostics.DiagOut(eDLWarn, gszProcName, "Warning: Defaulting number of threads to %d", MaxAllowedThreads);
		NumThreads = MaxAllowedThreads;
		}

	szRefAssembFile[0] = '\0';
	if (PMode == ePBAu2PBA || PMode == ePBAu2AVCF || PMode == ePBAu2GVCF || PMode == ePBAu2DVCF)
		{
		if (refassemb->count)
			{
			strncpy(szRefAssemb, refassemb->sval[0], cMaxRefAssembName);
			szRefAssemb[cMaxRefAssembName] = '\0';
			CUtility::TrimQuotedWhitespcExtd(szRefAssemb);
			CUtility::CleanText(szRefAssemb);
			}
		else
			szRefAssemb[0] = '\0';
		if (szRefAssemb[0] == '\0')
			strcpy(szRefAssemb, cDfltRefAssemb);
		szRefAssemb[cMaxRefAssembName] = '\0';

		if (exprid->count)
			{
			strncpy(szExprID, exprid->sval[0], cMaxRefAssembName);
			szExprID[cMaxRefAssembName] = '\0';
			CUtility::TrimQuotedWhitespcExtd(szExprID);
			CUtility::CleanText(szExprID);
			}
		else
			szExprID[0] = '\0';
		if (szExprID[0] == '\0')
			strcpy(szExprID, cDfltExprID);
		
		if (seqid->count)
			{
			strncpy(szSeqID, seqid->sval[0], cMaxRefAssembName);
			szSeqID[cMaxRefAssembName] = '\0';
			CUtility::TrimQuotedWhitespcExtd(szSeqID);
			CUtility::CleanText(szSeqID);
			}
		else
			szSeqID[0] = '\0';
		if (szSeqID[0] == '\0')
			strcpy(szSeqID, cDfltSeqID);
		szSeqID[cMaxRefAssembName] = '\0';
		if(PMode == ePBAu2AVCF || PMode == ePBAu2GVCF || PMode == ePBAu2DVCF)
			{
			if (refassembfile->count)
				{
				strncpy(szRefAssembFile, refassembfile->filename[0], _MAX_PATH);
				szRefAssembFile[_MAX_PATH - 1] = '\0';
				CUtility::TrimQuotedWhitespcExtd(szRefAssembFile);
				}
			else
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "No reference PBA assembly file specified, required for VCF generation");
				exit(1);
				}
			}
		else
			szRefAssembFile[0] = '\0';
		}
	else
		{
		szRefAssembFile[0] = '\0';
		szRefAssemb[0] = '\0';
		szExprID[0] = '\0';
		szSeqID[0] = '\0';
		}

	if (chromfile->count)
		{
		strncpy(szChromFile, chromfile->filename[0], _MAX_PATH);
		szChromFile[_MAX_PATH-1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(szChromFile);
		}
	else
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "No BED file containing chromosome names and sizes specified");
		exit(1);
		}

	NumIncludeChroms = 0;
	if (IncludeChroms->count)
		{
		for (NumIncludeChroms = Idx = 0; NumIncludeChroms < cMaxIncludeChroms && Idx < IncludeChroms->count; Idx++)
			{
			pszIncludeChroms[Idx] = nullptr;
			if (pszIncludeChroms[NumIncludeChroms] == nullptr)
				pszIncludeChroms[NumIncludeChroms] = new char[_MAX_PATH];
			strncpy(pszIncludeChroms[NumIncludeChroms], IncludeChroms->sval[Idx], _MAX_PATH);
			pszIncludeChroms[NumIncludeChroms][_MAX_PATH - 1] = '\0';
			CUtility::TrimQuotedWhitespcExtd(pszIncludeChroms[NumIncludeChroms]);
			if (pszIncludeChroms[NumIncludeChroms][0] != '\0')
				NumIncludeChroms++;
			}
		}
	

	NumExcludeChroms = 0;
	if (ExcludeChroms->count)
		{
		for (NumExcludeChroms = Idx = 0; NumExcludeChroms < cMaxExcludeChroms && Idx < ExcludeChroms->count; Idx++)
			{
			pszExcludeChroms[Idx] = nullptr;
			if (pszExcludeChroms[NumExcludeChroms] == nullptr)
				pszExcludeChroms[NumExcludeChroms] = new char[_MAX_PATH];
			strncpy(pszExcludeChroms[NumExcludeChroms], ExcludeChroms->sval[Idx], _MAX_PATH);
			pszExcludeChroms[NumExcludeChroms][_MAX_PATH - 1] = '\0';
			CUtility::TrimQuotedWhitespcExtd(pszExcludeChroms[NumExcludeChroms]);
			if (pszExcludeChroms[NumExcludeChroms][0] != '\0')
				NumExcludeChroms++;
			}
		}

	NumInputFiles = 0;
	if (infiles->count)
		{
		for (Idx = 0; NumInputFiles < cMaxWildCardFileSpecs && Idx < infiles->count; Idx++)
			{
			pszInputFiles[Idx] = nullptr;
			if (pszInputFiles[NumInputFiles] == nullptr)
				pszInputFiles[NumInputFiles] = new char[_MAX_PATH];
			strncpy(pszInputFiles[NumInputFiles], infiles->filename[Idx], _MAX_PATH);
			pszInputFiles[NumInputFiles][_MAX_PATH - 1] = '\0';
			CUtility::TrimQuotedWhitespcExtd(pszInputFiles[NumInputFiles]);
			if (pszInputFiles[NumInputFiles][0] != '\0')
				NumInputFiles++;
			}
		}

	if (!NumInputFiles)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: After removal of whitespace, there are no input file(s) specified with '-i<filespec>' option)\n");
		exit(1);
		}

	strcpy(szOutFile, outfile->filename[0]);
	CUtility::TrimQuotedWhitespcExtd(szOutFile);
	if(szOutFile[0] == '\0')
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "No output file specified");
		exit(1);
		}


	szROIFile[0] = '\0';
	szGTSampleFiltFile[0] = '\0';
	if((PMode == ePBAu2AVCF || PMode == ePBAu2GVCF || PMode == ePBAu2DVCF))
		{
		if(roifile->count)
			{
			strncpy(szROIFile, roifile->filename[0], _MAX_PATH);
			szROIFile[_MAX_PATH-1] = '\0';
			CUtility::TrimQuotedWhitespcExtd(szROIFile);
			}

		if((PMode == ePBAu2GVCF  || PMode == ePBAu2DVCF) && gtfiltsamplefile)
			{
			strncpy(szGTSampleFiltFile, gtfiltsamplefile->filename[0], _MAX_PATH);
			szGTSampleFiltFile[_MAX_PATH-1] = '\0';
			CUtility::TrimQuotedWhitespcExtd(szGTSampleFiltFile);
			}
		}

	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Processing parameters:");
	const char* pszDescr;
	switch(PMode) {
		case ePBAu2Fasta:
			pszDescr = "generating Fasta file from PBA file";
			break;
		case ePBAu2PBA:
			pszDescr = "generating PBA file from Fasta file";
			break;
		case ePBAu2PBAConcordance:
			pszDescr = "generating concordance CSV file from PBAs samples";
			break;
		case ePBAu2WIGConcordance:
			pszDescr = "generating concordance CSV file from WIG samples";
			break;
		case ePBAu2AVCF:
			pszDescr = "generate allelic VCF file";
			break;
		case ePBAu2GVCF:
			pszDescr = "generate genotype VCF file";
			break;
		case ePBAu2DVCF:
			pszDescr = "generate a deletion VCF";
			break;
		case ePBAu2BED:
			pszDescr = "generate BED containing all transcribed regions";
			break;
	}

	gDiagnostics.DiagOutMsgOnly(eDLInfo, "PBA utilities : '%s'", pszDescr);

	if (PMode == ePBAu2PBAConcordance || PMode == ePBAu2BED || PMode == ePBAu2WIGConcordance)
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "Limit number of processed PBAs or WIGs to : %d", LimitPBAs);

	gDiagnostics.DiagOutMsgOnly(eDLInfo, "Trim 5' segments by : %dbp", PBAsTrim5);
	gDiagnostics.DiagOutMsgOnly(eDLInfo, "Trim 3' segments by : %dbp", PBAsTrim3);

	if(PMode == ePBAu2GVCF || PMode == ePBAu2DVCF)
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "Maximum proportion unaligned samples : %.3f", GTPropNAThres);
		if(PMode == ePBAu2GVCF)
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "Minimum proportion heterozygosity between aligned samples : %.3f", GTPropHetThres);
		else
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "Minimum proportion unaligned samples : %.3f", GTPropHetThres);
		}

	if(PMode == ePBAu2BED)
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "Minimum proportion aligned samples segment initiation : %.3f", GTPropNAThres);
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "Minimum proportion aligned samples segment continuation : %.3f", GTPropHetThres);
		}

	if(PMode == ePBAu2AVCF || PMode == ePBAu2GVCF || PMode == ePBAu2DVCF)
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "PBA reference assembly : '%s'", szRefAssembFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo, "BED containing chromosome names and sizes : '%s'", szChromFile);
	if (PMode == ePBAu2PBA || PMode == ePBAu2WIGConcordance)
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "Reference assembly : '%s'", szRefAssemb);
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "Experiment identifier : '%s'", szExprID);
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "Sequence identifier : '%s'", szSeqID);
		}

	if(szROIFile[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "Regions of interest file : '%s'", szROIFile);

	if(szGTSampleFiltFile[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "Only report on GT samples in this file : '%s'", szGTSampleFiltFile);

	for (Idx = 0; Idx < NumIncludeChroms; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "reg expressions defining chroms to include: '%s'", pszIncludeChroms[Idx]);
	for (Idx = 0; Idx < NumExcludeChroms; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "reg expressions defining chroms to exclude: '%s'", pszExcludeChroms[Idx]);

	for (Idx = 0; Idx < NumInputFiles; Idx++)
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "Input file spec (%d) : '%s'", Idx+1, pszInputFiles[Idx]);

	gDiagnostics.DiagOutMsgOnly(eDLInfo, "Output to file : '%s'", szOutFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo, "number of threads : %d", NumThreads);
#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = 0;
	Rslt = Process(PMode,	// processing mode: ePBAu2Fasta PBA to Fasta, ePBAu2PBA Fasta to PBA, ePBAu2PBAConcordance concordance over PBA samples, ePBAu2WIGConcordance concordance over WIG samples, ePBAu2AVCF allelic variant VCF, ePBAu2GVCF allelic genotype VCF, ePBAu2DVCF deletion genotype VCF, ePBAu2BED generate BED containing all transcribed regions
					LimitPBAs,			// limit number of loaded PBA files to this many. 1 .. cMaxPBAReadsets
					PBAsTrim5,			// trim this many aligned PBAs from 5' end of aligned segments - reduces false alleles due to sequencing errors
					PBAsTrim3,			// trim this many aligned PBAs from 3' end of aligned segments - reduces false alleles due to sequencing errors
					GTPropNAThres,		// when genotyping VCF then proportion of non-aligned calls over all samples must be < this threshold
					GTPropHetThres,		// when genotyping VCF then proportion of heterozygous calls must be >= this threshold
					szRefAssemb,		// reference assembly
					szRefAssembFile,	// PBA file containing reference assembly sequences
					szChromFile,		// BED file containing chromosome names and sizes
					szSeqID,			// sequence identifier
					szExprID,			// experiment identifier
					NumInputFiles,		// number of input founder file specs
					pszInputFiles,		// names of input founder PBA files (wildcards allowed)
					szOutFile,			// output to this file
					szROIFile,			// optional BED file containing regions of interest
					szGTSampleFiltFile,	// optional file containing those readseq/sample names for whom only their genotypes are to be reported
					NumIncludeChroms,	// number of chromosome regular expressions to include
					pszIncludeChroms,	// array of include chromosome regular expressions
					NumExcludeChroms,	// number of chromosome expressions to exclude
					pszExcludeChroms,	// array of exclude chromosome regular expressions
					NumThreads);		// number of worker threads to use

	Rslt = Rslt >= 0 ? 0 : 1;
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
}

int Process(ePBAuMode PMode,	// processing mode: ePBAu2Fasta PBA to Fasta, ePBAu2PBA Fasta to PBA, ePBAu2PBAConcordance concordance over PBA samples, ePBAu2WIGConcordance concordance over WIG samples, ePBAu2AVCF allelic variant VCF, ePBAu2GVCF allelic genotype VCF, ePBAu2DVCF deletion genotype VCF,ePBAu2BED generate BED containing all transcribed regions
	int32_t LimitPBAs,			// limit number of loaded PBA files to this many. 1 .. cMaxPBAReadsets
	int32_t PBAsTrim5,			// trim this many aligned PBAs from 5' end of aligned segments - reduces false alleles due to sequencing errors
	int32_t PBAsTrim3,			// trim this many aligned PBAs from 3' end of aligned segments - reduces false alleles due to sequencing errors
	double GTPropNAThres,		// when genotyping VCF then proportion of non-aligned calls over all samples must be < this threshold
	double GTPropHetThres,		// when genotyping VCF then proportion of heterozygous calls must be >= this threshold
	char* pszRefAssemb,			// reference assembly identifier
	char* pszRefAssembFile,		// PBA file containing reference assembly sequences
	char* pszChromFile,			// BED file containing chromosome names and sizes
	char* pszSeqID,				// sequence identifier
	char* pszExprID,			// experiment identifier
	int32_t NumInputFiles,			// number of input founder file specs
	char* pszInputFiles[],		// names of input founder PBA files (wildcards allowed)
	char* pszOutFile,			// output to this file
	char *pszROIFile,			// optional BED file containing regions of interest
	char *pszGTSampleFiltFile,	// optional file containing those readseq/sample names for whom only their genotypes are to be reported
	int NumIncludeChroms,		// number of chromosome regular expressions to include
	char* pszIncludeChroms[],	// array of include chromosome regular expressions
	int NumExcludeChroms,		// number of chromosome expressions to exclude 
	char* pszExcludeChroms[],	// array of exclude chromosome regular expressions
	int NumThreads)				// number of worker threads to use
{
int Rslt;
CPBAutils *pPBAutils;
if((pPBAutils = new CPBAutils) == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to instantiate instance of CPBAutils");
	return(eBSFerrObj);
	}
Rslt = pPBAutils->Process(PMode, LimitPBAs, PBAsTrim5,PBAsTrim3,GTPropNAThres,GTPropHetThres,pszRefAssemb,pszRefAssembFile, pszChromFile, pszSeqID, pszExprID,
					NumInputFiles, pszInputFiles, pszOutFile, pszROIFile,pszGTSampleFiltFile,
					NumIncludeChroms, pszIncludeChroms, NumExcludeChroms, pszExcludeChroms, NumThreads);
delete pPBAutils;
return(Rslt);
}

CPBAutils::CPBAutils()	// constructor
{
m_hInFile = -1;
m_hOutFile = -1;
m_hWIGOutFile = -1;

m_pChromMetadata = nullptr;
m_pInBuffer = nullptr;
m_pOutBuffer = nullptr;
m_pszWIGBuff = nullptr;
m_pBedFile = nullptr;
m_pROIFile = nullptr;
m_pReportGTSample = nullptr;
m_bMutexesCreated = false;
}

CPBAutils::~CPBAutils()	// destructor
{
if (m_hInFile != -1)
	close(m_hInFile);
if (m_hOutFile != -1)
	close(m_hInFile);
if (m_hWIGOutFile != -1)
	close(m_hWIGOutFile);

if (m_pInBuffer != nullptr)
	delete[]m_pInBuffer;
if (m_pOutBuffer != nullptr)
	delete[]m_pOutBuffer;
if (m_pszWIGBuff != nullptr)
	delete []m_pszWIGBuff;
if (m_pBedFile != nullptr)
	delete m_pBedFile;
if(m_pROIFile != nullptr)
	delete m_pROIFile;

if(m_pReportGTSample != nullptr)
	delete []m_pReportGTSample;

if(m_pChromMetadata != nullptr)
	{
	tsPUChromMetadata* pChromMetadata = m_pChromMetadata;
	for (uint32_t Idx = 0; Idx < m_UsedNumChromMetadata; Idx++, pChromMetadata++)
		{
		if (pChromMetadata->pPBAs != nullptr)
#ifdef _WIN32
			free(pChromMetadata->pPBAs);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
			if (pChromMetadata->pPBAs != MAP_FAILED)
				munmap(pChromMetadata->pPBAs, pChromMetadata->ChromLen);
#endif
		}

#ifdef _WIN32
	free(m_pChromMetadata);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (m_pChromMetadata != MAP_FAILED)
		munmap(m_pChromMetadata, m_AllocdChromMetadataMem);
#endif
	}
if(m_bMutexesCreated)
	DeleteMutexes();
}

void
CPBAutils::Reset(void)
{
if (m_hInFile != -1)
	{
	close(m_hInFile);
	m_hInFile = -1;
	}
if (m_hOutFile != -1)
	{
	close(m_hOutFile);
	m_hOutFile = -1;
	}

if (m_hWIGOutFile != -1)
	{
	close(m_hWIGOutFile);
	m_hWIGOutFile = -1;
	}

if (m_pInBuffer != nullptr)
	{
	delete[]m_pInBuffer;
	m_pInBuffer = nullptr;
	}

if (m_pOutBuffer != nullptr)
	{
	delete[]m_pOutBuffer;
	m_pOutBuffer = nullptr;
	}

if (m_pszWIGBuff != nullptr)
	{
	delete[]m_pszWIGBuff;
	m_pszWIGBuff = nullptr;
	}

if (m_pBedFile != nullptr)
	{
	delete m_pBedFile;
	m_pBedFile = nullptr;
	}

if(m_pROIFile != nullptr)
	{
	delete m_pROIFile;
	m_pROIFile = nullptr;
	}

if(m_pReportGTSample != nullptr)
	{
	delete []m_pReportGTSample;
	m_pReportGTSample = nullptr;
	}

if (m_pChromMetadata != nullptr)
	{
	tsPUChromMetadata* pChromMetadata = m_pChromMetadata;
	for (uint32_t Idx = 0; Idx < m_UsedNumChromMetadata; Idx++, pChromMetadata++)
		{
		if (pChromMetadata->pPBAs != nullptr)
#ifdef _WIN32
			free(pChromMetadata->pPBAs);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
		if (pChromMetadata->pPBAs != MAP_FAILED)
				munmap(pChromMetadata->pPBAs, pChromMetadata->ChromLen);
#endif
		}

#ifdef _WIN32
	free(m_pChromMetadata);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (m_pChromMetadata != MAP_FAILED)
		munmap(m_pChromMetadata, m_AllocdChromMetadataMem);
#endif
	m_pChromMetadata = nullptr;
	}

m_FastSerialise = 0;

m_UsedNumChromMetadata = 0;
m_AllocdChromMetadata = 0;
m_AllocdChromMetadataMem = 0;

m_InNumBuffered = 0;
m_AllocInBuff = 0;

m_OutBuffIdx = 0;
m_AllocOutBuff = 0;

m_NumReadsetIDs = 0;

m_LAReadsetNameID = 0;
m_NumReadsetNames = 0;
m_NxtszReadsetIdx = 0;
m_szReadsetNames[0] = '\0';

m_LAChromNameID = 0;
m_NumChromNames = 0;
m_NxtszChromIdx = 0;
m_szChromNames[0] = '\0';

m_NumIncludeChroms = 0;
m_NumExcludeChroms = 0;

m_szRefAssemb[0] = '\0';
m_szRefAssembFile[0] = '\0';

m_NumChromSizes = 0;

m_szOutFile[0] = '\0';

m_NumThreads = cMaxPBAutilityThreads;
m_PBAsTrim5 = 0;
m_PBAsTrim3 = 0;

m_GTPropNAThres = 0.0;
m_GTPropHetThres = 0.0;

m_WIGChromID = 0;
m_WIGRptdChromID = 0;
m_WIGSpanLoci = 0;
m_WIGSpanLen = 0;
m_WIGRptdSpanLen = 0;
m_WIGSpanCnts = 0;
m_WIGBuffIdx = 0;
m_AllocWIGBuff = 0;

m_MRA_ROIChromID = 0;
m_MRA_ChromID = 0;

if (m_bMutexesCreated)
	DeleteMutexes();
m_bMutexesCreated = false;
}

void
CPBAutils::InitialiseWIGSpan(void) // initialise WIG span vars to values corresponding to no spans having been previously reported
{
	m_WIGChromID = 0;
	m_WIGRptdChromID = 0;
	m_WIGSpanLoci = 0;
	m_WIGSpanLen = 0;
	m_WIGRptdSpanLen = 0;
	m_WIGSpanCnts = 0;
}

int
CPBAutils::CompleteWIGSpan(bool bWrite)				// close off any current WIG span ready to start any subsequent span
{
	char* pszChrom;

	// if existing span then write that span out
	if (m_WIGChromID != 0 && m_WIGSpanLen > 0 && m_WIGSpanLoci > 0 && m_WIGSpanCnts > 0)
	{
		// has chrom and/or span changed since previously writing out a span?
		if (m_WIGChromID != m_WIGRptdChromID || m_WIGSpanLen != m_WIGRptdSpanLen)
		{
			pszChrom = LocateChrom(m_WIGChromID);
			m_WIGBuffIdx += sprintf((char*)&m_pszWIGBuff[m_WIGBuffIdx], "variableStep chrom=%s span=%d\n", pszChrom, m_WIGSpanLen);
			m_WIGRptdChromID = m_WIGChromID;
			m_WIGRptdSpanLen = m_WIGSpanLen;
		}
		m_WIGBuffIdx += sprintf((char*)&m_pszWIGBuff[m_WIGBuffIdx], "%d %d\n", m_WIGSpanLoci, (uint32_t)m_WIGSpanCnts);
	}
	if ((bWrite && m_WIGBuffIdx) || (m_WIGBuffIdx + 500) > m_AllocWIGBuff)
	{
		if (!CUtility::RetryWrites(m_hWIGOutFile, m_pszWIGBuff, m_WIGBuffIdx))
		{
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "RetryWrites() error");
			return(eBSFerrWrite);
		}
		m_WIGBuffIdx = 0;
	}
	m_WIGSpanLoci = 0;
	m_WIGSpanLen = 0;
	m_WIGSpanCnts = 0;
	return(eBSFSuccess);
}



int
CPBAutils::AccumWIGCnts(uint32_t ChromID,	// accumulate wiggle counts into variableStep spans over this chromosome
	uint32_t Loci,		// this loci - starts from 1 not 0!
	uint32_t BinLen,     // span bin is this length
	uint32_t Cnts,		 // has this many counts attributed to the span
	uint32_t MaxSpanLen) // allow WIG spans to be this maximal length
{
	int Rslt;
	uint32_t Meanx100;
	if (ChromID != m_WIGChromID || m_WIGSpanLen >= MaxSpanLen || Cnts == 0)		// onto a different chromosome, or current span is at maximal length?
	{
		if (m_WIGChromID != 0)
		{
			if ((Rslt = CompleteWIGSpan()) < 0)
				return(Rslt);
		}
		if (Cnts > 0)
		{
			m_WIGChromID = ChromID;
			m_WIGSpanLoci = Loci;
			m_WIGSpanLen = BinLen;
			m_WIGSpanCnts = Cnts;
		}
		return(eBSFSuccess);
	}

	if (m_WIGSpanLen == 0 || m_WIGSpanCnts == 0) // starting a new span?
	{
		m_WIGSpanLoci = Loci;
		m_WIGSpanLen = BinLen;
		m_WIGSpanCnts = (uint64_t)Cnts;
		return(eBSFSuccess);
	}

	if (Cnts != m_WIGSpanCnts)
	{
		// if less than 5% difference between existing and new span counts then combine adjacent bins into a single span  
		Meanx100 = 100 * (uint32_t)(m_WIGSpanCnts / (uint64_t)m_WIGSpanLen);
		if ((Cnts <= 100 && (Cnts * 100) != Meanx100) || (Meanx100 < (Cnts * 95) || Meanx100 >= (Cnts * 105)))
		{
			// write current span out
			if ((Rslt = CompleteWIGSpan()) < 0)
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

int 
CPBAutils::Process(ePBAuMode PMode,	// processing mode: ePBAu2Fasta PBA to Fasta, ePBAu2PBA Fasta to PBA, ePBAu2PBAConcordance concordance over PBA samples, ePBAu2WIGConcordance concordance over WIG samples, ePBAu2AVCF allelic variant VCF, ePBAu2GVCF allelic genotype VCF, ePBAu2DVCF deletion genotype VCF
	int32_t LimitPBAs,			// limit number of loaded PBA files to this many. 1...cMaxPBAReadsets
	int32_t PBAsTrim5,			// trim this many aligned PBAs from 5' end of aligned segments - reduces false alleles due to sequencing errors
	int32_t PBAsTrim3,			// trim this many aligned PBAs from 3' end of aligned segments - reduces false alleles due to sequencing errors
	double GTPropNAThres,		// when genotyping VCF then proportion of non-aligned calls over all samples must be < this threshold
	double GTPropHetThres,		// when genotyping VCF then proportion of heterozygous calls must be >= this threshold
	char* pszRefAssemb,			// reference assembly identifier
	char* pszRefAssembFile,		// PBA file containing reference assembly sequences
	char* pszChromFile,			// BED file containing chromosome names and sizes
	char* pszSeqID,				// sequence identifier
	char* pszExprID,			// experiment identifier
	int NumInputFiles,			// number of input founder file specs
	char* pszInputFiles[],		// names of input founder PBA files (wildcards allowed)
	char* pszOutFile,			// output to this file
	char *pszROIFile,			// optional BED file containing regions of interest
	char *pszGTSampleFiltFile,	// optional file containing those readseq/sample names for whom only their genotypes are to be reported
	int NumIncludeChroms,		// number of chromosome regular expressions to include
	char* pszIncludeChroms[],	// array of include chromosome regular expressions
	int NumExcludeChroms,		// number of chromosome expressions to exclude
	char* pszExcludeChroms[],	// array of exclude chromosome regular expressions
	int NumThreads)				// number of worker threads to use
{
int Rslt;
Reset();
CreateMutexes();
m_PBAsTrim5 = PBAsTrim5;
m_PBAsTrim3 = PBAsTrim3;
m_GTPropNAThres = GTPropNAThres;
m_GTPropHetThres = GTPropHetThres;
strcpy(m_szRefAssemb, pszRefAssemb);
strcpy(m_szOutFile, pszOutFile);
if(pszRefAssembFile != nullptr && pszRefAssembFile[0] != '\0')
	strcpy(m_szRefAssembFile, pszRefAssembFile);

// compile include/exclude chromosome regexpr if user has specified alignments to be filtered by chrom
if(Rslt = (m_RegExprs.CompileREs(NumIncludeChroms, pszIncludeChroms,NumExcludeChroms, pszExcludeChroms)) < eBSFSuccess)
	{
	Reset();
	return(Rslt);
	}

// if specified then load regions of interest from file
m_pROIFile = nullptr;
if(pszROIFile != nullptr && pszROIFile[0] != '\0')
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading regions of interest file '%s'", pszROIFile);
	if((m_pROIFile = new CBEDfile()) == nullptr)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBEDfile");
		Reset();
		return(eBSFerrObj);
		}
	if((Rslt=m_pROIFile->Open(pszROIFile))!=eBSFSuccess)
		{
		while(m_pROIFile->NumErrMsgs())
			gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pROIFile->GetErrMsg());
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open regions of interest file '%s'",pszROIFile);
		Reset();
		return(Rslt);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Regions of interest file '%s' loaded", pszROIFile);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing chromosome sizes from file '%s'", pszChromFile);
if ((Rslt = LoadChromSizes(pszChromFile)) < 1) // BED file containing chromosome names and sizes - NOTE chromosomes will be filtered by include/exclude wildcards
	{
	Reset();
	return(Rslt);
	}

m_NumThreads = NumThreads;

	// initial allocation, will be realloc'd if more memory required
size_t memreq = (size_t)cAllocChromMetadata * sizeof(tsPUChromMetadata);
#ifdef _WIN32
m_pChromMetadata = (tsPUChromMetadata*)malloc(memreq);	// initial and perhaps the only allocation
if (m_pChromMetadata == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Initial memory allocation of %zd bytes for chromosome metadata failed - %s", (int64_t)memreq, strerror(errno));
	Reset();
	return(eBSFerrMem);
	}
#else
m_pChromMetadata = (tsPUChromMetadata*)mmap(nullptr, memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
if (m_pChromMetadata == MAP_FAILED)
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

if ((m_pInBuffer = new uint8_t[cInBuffSize]) == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Memory allocation of %u bytes for input buffering failed - %s", cInBuffSize, strerror(errno));
	Reset();
	return(eBSFerrMem);
	}
m_AllocInBuff = cInBuffSize;
m_InNumBuffered = 0;
int32_t Idx;
char* pszInFile;
int32_t ReadsetID = 0;
int32_t RefReadsetID = 0;
int32_t NumFiles = 0;
int32_t TotNumFiles = 0;
CSimpleGlob glob(SG_GLOB_FULLSORT);
if (PMode <= ePBAu2PBA)
	{
	glob.Init();
	if (glob.Add(pszInputFiles[0]) < SG_SUCCESS)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to glob '%s' input file spec", pszInputFiles[0]);
		Reset();
		return(eBSFerrOpnFile);	// treat as though unable to open file
		}
	if ((NumFiles = glob.FileCount()) < 1)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to glob '%s' input file spec", pszInputFiles[0]);
		Reset();
		return(eBSFerrOpnFile);	// treat as though unable to open file
		}

	pszInFile = glob.File(0);
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Process: Starting to load input %s file", PMode == ePBAu2Fasta ? "PBA" : "Fasta");

	if (PMode == ePBAu2Fasta)
		{
		ReadsetID = LoadPBAFile(pszInFile, 0, true); // loading readset + chrom metadata only, later the PBAs for each individual chrom will be loaded on demand
		if (ReadsetID <= 0)
			{
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "Process: Errors loading file '%s'", pszInFile);
			Reset();
			return(ReadsetID);
			}
		Rslt = PBA2Fasta(1);
		}
	else
		Rslt = Fasta2PBA(pszRefAssemb, pszSeqID, pszExprID, pszInFile, m_szOutFile);
	Reset();
	return(Rslt);
	}

if (PMode == ePBAu2AVCF)
	{
	glob.Init();
	if (glob.Add(pszRefAssembFile) < SG_SUCCESS)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to glob '%s' reference PBA file spec", pszRefAssembFile);
		Reset();
		return(eBSFerrOpnFile);	// treat as though unable to open file
		}
	if ((NumFiles = glob.FileCount()) <= 0)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to locate any reference PBA file matching '%s'", pszRefAssembFile);
		Reset();
		return(eBSFerrOpnFile);	// treat as though unable to open file
		}
	if (NumFiles > 1)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Only one (%d located) reference PBA file matching '%s' allowed", NumFiles, pszRefAssembFile);
		Reset();
		return(eBSFerrOpnFile);	// treat as though unable to open file
		}
	pszInFile = glob.File(0);
	RefReadsetID = LoadPBAFile(pszInFile, 0, true); // loading readset + chrom metadata only, later the PBAs for each individual chrom will be loaded on demand
	if (RefReadsetID <= 0)
		{
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Process: Errors loading reference PBA file '%s'", pszInFile);
		Reset();
		return(RefReadsetID);
		}

	glob.Init();
	if (glob.Add(pszInputFiles[0]) < SG_SUCCESS)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to glob '%s' file spec", pszInputFiles[0]);
		Reset();
		return(eBSFerrOpnFile);	// treat as though unable to open file
		}
	if ((NumFiles = glob.FileCount()) <= 0)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to locate any file matching '%s'", pszInputFiles[0]);
		Reset();
		return(eBSFerrOpnFile);	// treat as though unable to open file
		}
	if (NumFiles > 1)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Only one (%d located) file matching '%s' allowed", NumFiles, pszInputFiles[0]);
		Reset();
		return(eBSFerrOpnFile);	// treat as though unable to open file
		}
	pszInFile = glob.File(0);
	ReadsetID = LoadPBAFile(pszInFile, 0, true); // loading readset + chrom metadata only, later the PBAs for each individual chrom will be loaded on demand
	if (ReadsetID <= 0)
		{
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Process: Errors loading file '%s'", pszInFile);
		Reset();
		return(ReadsetID);
		}

	Rslt = GenAllelicVCF(RefReadsetID,		// reference PBA assembly
						   ReadsetID);		// PBA file with allelic differences
	Reset();
	return(Rslt);
	}

if (PMode == ePBAu2GVCF || PMode == ePBAu2DVCF)
	{
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Process: Starting to load input PBA files");
	Rslt = eBSFSuccess;		// assume success!
	TotNumFiles = 0;
	glob.Init();
	if (glob.Add(pszRefAssembFile) < SG_SUCCESS)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to glob '%s' reference PBA file spec", pszRefAssembFile);
		Reset();
		return(eBSFerrOpnFile);	// treat as though unable to open file
		}
	if ((NumFiles = glob.FileCount()) <= 0)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to locate any reference PBA file matching '%s'", pszRefAssembFile);
		Reset();
		return(eBSFerrOpnFile);	// treat as though unable to open file
		}
	if (NumFiles > 1)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Only one (%d located) reference PBA file matching '%s' allowed", NumFiles, pszRefAssembFile);
		Reset();
		return(eBSFerrOpnFile);	// treat as though unable to open file
		}
	pszInFile = glob.File(0);
	RefReadsetID = LoadPBAFile(pszInFile, 0, true); // loading readset + chrom metadata only, later the PBAs for each individual chrom will be loaded on demand
	if (RefReadsetID <= 0)
		{
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Process: Errors loading reference PBA file '%s'", pszInFile);
		Reset();
		return(RefReadsetID);
		}

	for (Idx = 0; Idx < NumInputFiles; Idx++)
		{
		glob.Init();
		if (glob.Add(pszInputFiles[Idx]) < SG_SUCCESS)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to glob '%s' input file spec", pszInputFiles[Idx]);
			Reset();
			return(eBSFerrOpnFile);	// treat as though unable to open file
			}
		if ((NumFiles = glob.FileCount()) <= 0)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to locate any input file matching '%s'", pszInputFiles[Idx]);
			Reset();
			return(eBSFerrFileName);
			}

		TotNumFiles += NumFiles;

		Rslt = eBSFSuccess;
		for (int FileID = 0; Rslt >= eBSFSuccess && FileID < NumFiles; ++FileID)
			{
			pszInFile = glob.File(FileID);
			ReadsetID = LoadPBAFile(pszInFile, 0, true); // loading readset + chrom metadata only, later the PBAs for each individual chrom will be loaded on demand
			if (ReadsetID <= 0)
				{
				gDiagnostics.DiagOut(eDLInfo, gszProcName, "Process: Errors loading file '%s'", pszInFile);
				Reset();
				return(ReadsetID);
				}
			if (ReadsetID == LimitPBAs)
				break;
			}
		if (ReadsetID == LimitPBAs)
			{
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "Process: Loaded maximum supported (%d) number of PBA files", LimitPBAs);
			break;
			}
		}
	m_NumReadsetIDs = ReadsetID;

	if(pszGTSampleFiltFile != nullptr && (Rslt = InitialiseGTSampleFiltering(pszGTSampleFiltFile)) < 0)
		{
		Reset();
		return(Rslt);
		}

	Rslt = GenGenotypeVCF(PMode,			// will be either ePBAu2GVCF allelic genotype VCF or ePBAu2DVCF deletion genotype VCF
						RefReadsetID,		// reference PBA assembly
						   m_NumReadsetIDs-1,			// number of PBA samples to report on genotypes relative to reference
						   GTPropNAThres,		// when genotyping VCF then proportion of non-aligned calls over all samples must be < this threshold
							GTPropHetThres);		// when genotyping VCF then proportion of heterozygous calls must be >= this threshold
	Reset();
	return(Rslt);
	}

	// concordance analysis
gDiagnostics.DiagOut(eDLInfo, gszProcName, "Process: Starting to load input PBA files");
Rslt = eBSFSuccess;		// assume success!
TotNumFiles = 0;
for (Idx = 0; Idx < NumInputFiles; Idx++)
	{
	glob.Init();
	if (glob.Add(pszInputFiles[Idx]) < SG_SUCCESS)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to glob '%s' input file spec", pszInputFiles[Idx]);
		Reset();
		return(eBSFerrOpnFile);	// treat as though unable to open file
		}
	if ((NumFiles = glob.FileCount()) <= 0)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to locate any input file matching '%s'", pszInputFiles[Idx]);
		Reset();
		return(eBSFerrFileName);
		}

	TotNumFiles += NumFiles;

	Rslt = eBSFSuccess;
	for (int FileID = 0; Rslt >= eBSFSuccess && FileID < NumFiles; ++FileID)
		{
		pszInFile = glob.File(FileID);
		if (PMode == ePBAu2PBA)
			ReadsetID = LoadPBACoverage(pszInFile);		 // will load readset + chrom metadata only from PBA after checking that coverage WIG file actually exists
		else
			ReadsetID = LoadPBAFile(pszInFile, 0, true); // loading readset + chrom metadata only, later the PBAs for each individual chrom will be loaded on demand
		if (ReadsetID <= 0)
			{
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "Process: Errors loading file '%s'", pszInFile);
			Reset();
			return(ReadsetID);
			}
		if (ReadsetID == LimitPBAs)
			break;
		}
	if (ReadsetID == LimitPBAs)
		{
		gDiagnostics.DiagOut(eDLWarn, gszProcName, "Process: Loaded maximum supported (%d) number of PBA files", LimitPBAs);
		break;
		}
	}
m_NumReadsetIDs = ReadsetID;
if(PMode == ePBAu2BED)
	Rslt = GeneratePBASegsBed(m_GTPropNAThres,m_GTPropHetThres,cMinSegLen);
else
	Rslt = GeneratePBAConcordance();
return(Rslt);
}

// check if assumed ascii string only contains ascii chars
bool
CPBAutils::CheckIfAscii(int Len,	// length to check
			 uint8_t* pAscii)	// assumed ascii - should only contain printable characters
{
if(Len > 0)
	{
	while(Len--)
		if(!isascii((char) * pAscii++))
			return(false);
	}
return(true);
}

int32_t
CPBAutils::InitialiseGTSampleFiltering(char *pszGTSampleFiltFile) // pszFiltFile contains sample names which are to be reported on in the generated GT VCF file, filter file contains 1 name per line
{
int Num2Accept;
int32_t NameLen;
int32_t SampleID;
char *pTxt;
char Chr;
int32_t LineNum;
int32_t LineLen;
FILE *pStream;
char szSampleName[200];
m_pReportGTSample = nullptr;
if(pszGTSampleFiltFile == nullptr || pszGTSampleFiltFile[0] == '\0')
	return(0);

m_pReportGTSample = new uint8_t[m_NumReadsetIDs];
memset(m_pReportGTSample,0,m_NumReadsetIDs); // default - assume no samples will be accepted for reporting

if((pStream = fopen(pszGTSampleFiltFile,"r"))==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "InitialiseGTSampleFiltering: Unable to access '%s' error: %s", pszGTSampleFiltFile,strerror(errno));
	return(eBSFerrOpnFile);
	}
LineNum = 0;
while(fgets((char *)m_pInBuffer,m_AllocInBuff-1,pStream)!= NULL)
	{
	LineNum += 1;
	LineLen = (int32_t)strlen((char *)m_pInBuffer);
	if(!CheckIfAscii(LineLen,m_pInBuffer))
		{
		fclose(pStream);
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "InitialiseGTSampleFiltering: Unable to parse '%s' at line", pszGTSampleFiltFile,LineNum);
		return(eBSFerrParse);
		}

	pTxt = (char *)m_pInBuffer;
	pTxt = CUtility::TrimQuotedWhitespcExtd(pTxt);
	NameLen = 0;
	while((Chr = *pTxt++) != '\0')
		{
		if(!isspace(Chr) && Chr != '\'' && Chr != '"' && Chr != ',' && Chr != ';')
			{
			if(NameLen < sizeof(szSampleName) - 1)
				szSampleName[NameLen++] = Chr;
			continue;
			}
		if(NameLen == 0)
			continue;
		szSampleName[NameLen] = '\0';
		NameLen = 0;
		if((SampleID = LocateReadset(szSampleName,0)) <= 1)
			continue;
		m_pReportGTSample[SampleID-1] |= 0x01;
		}
	if(NameLen > 0)
		{
		szSampleName[NameLen] = '\0';
		if((SampleID = LocateReadset(szSampleName,0)) <= 1)
			continue;
		m_pReportGTSample[SampleID-1] |= 0x01;
		}
	}
fclose(pStream);

// count number of accepted samples
Num2Accept = 0;
for(SampleID = 1; SampleID < m_NumReadsetIDs; SampleID++)
	if(m_pReportGTSample[SampleID] & 0x01)
		Num2Accept++;
if(Num2Accept == 0)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "InitialiseGTSampleFiltering: No readset sample names in '%s' specified which match any loaded readsets - nothing to do!", pszGTSampleFiltFile);
	return(eBSFerrOpnFile);
	}
else
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "InitialiseGTSampleFiltering: Reporting on %d genotyped samples out of %d to be processed ", Num2Accept,m_NumReadsetIDs-1);
return(Num2Accept);
}

// loading BED which specifies chrom names and sizes
int		// returning number of chromosomes parsed from BED file and accepted after filtering for wildcards
CPBAutils::LoadChromSizes(char* pszBEDFile) // BED file containing chromosome names and sizes
{
int Rslt;
uint32_t ChromID;
int CurFeatureID;
char szFeatName[cMaxDatasetSpeciesChrom];
char szChromName[cMaxDatasetSpeciesChrom];
int32_t StartLoci;
int32_t EndLoci;
int NumChroms;

if ((m_pBedFile = new CBEDfile) == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to instantiate CBEDfile");
	return(eBSFerrObj);
	}

if ((Rslt = m_pBedFile->Open(pszBEDFile, eBTAnyBed)) != eBSFSuccess)
	{
	while (m_pBedFile->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal, gszProcName, m_pBedFile->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadChromSizes: Unable to process '%s' as BED", pszBEDFile);
	return(eBSFerrOpnFile);
	}

NumChroms = 0;
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
		return(eBSFerrChrom);
		}

	if (EndLoci < 0 || (size_t)EndLoci >= cAllocPackedBaseAlleles)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadChromSizes: unable to accept chromosome '%s' because end loci not in acceptable range, must be in range 0..%zd", szChromName, cAllocPackedBaseAlleles);
		return(eBSFerrChrom);
		}

	// expecting a single instance of each chromosome 
	if (LocateChrom(szChromName) != 0)	// must be a unique chromosome name
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadChromSizes: multiple instances of chromosome '%s', chromosome names must be unique", szChromName);
		return(eBSFerrChrom);
		}

	ChromID = AddChrom(szChromName);		// new chromosome
	m_ChromSizes[ChromID-1] = EndLoci+1;	// record it's size
	NumChroms += 1;
	}
m_NumChromSizes = NumChroms;
gDiagnostics.DiagOut(eDLInfo, gszProcName, "LoadChromSizes: accepted %d chromosome sizes from '%s'", NumChroms, pszBEDFile);
delete m_pBedFile;
m_pBedFile = nullptr;
return(NumChroms);
}


int
CPBAutils::GeneratePBAConcordance(void)	// generate allelic concordance across all sample PBAs
{
if (m_pOutBuffer == nullptr)
	{
	if ((m_pOutBuffer = new uint8_t[cOutBuffSize]) == nullptr) // an overkill in terms of buffer size!
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Memory allocation of %u bytes for output buffering failed - %s", cOutBuffSize, strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
	m_AllocOutBuff = cOutBuffSize;
	}
m_OutBuffIdx = 0;

if (m_pszWIGBuff == nullptr)
{
	if ((m_pszWIGBuff = new uint8_t[cOutBuffSize]) == nullptr) // an overkill in terms of buffer size!
	{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Memory allocation of %u bytes for WIG output buffering failed - %s", cOutBuffSize, strerror(errno));
		Reset();
		return(eBSFerrMem);
	}
	m_AllocWIGBuff = cOutBuffSize;
}
m_WIGBuffIdx = 0;

CUtility::AppendFileNameSuffix(m_szWIGFile, m_szOutFile, (char *)".wig", '.');

#ifdef _WIN32
m_hOutFile = open(m_szOutFile, (O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC), (_S_IREAD | _S_IWRITE));
#else
if ((m_hOutFile = open64(m_szOutFile, O_WRONLY | O_CREAT, S_IREAD | S_IWRITE)) != -1)
if (ftruncate(m_hOutFile, 0) != 0)
{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to create/truncate %s - %s", m_szOutFile, strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
}
#endif
if (m_hOutFile < 0)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to create/truncate %s - %s", m_szOutFile, strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}

#ifdef _WIN32
m_hWIGOutFile = open(m_szWIGFile, (O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC), (_S_IREAD | _S_IWRITE));
#else
if ((m_hWIGOutFile = open64(m_szWIGFile, O_WRONLY | O_CREAT, S_IREAD | S_IWRITE)) != -1)
if (ftruncate(m_hWIGOutFile, 0) != 0)
{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to create/truncate %s - %s", m_szWIGFile, strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
}
#endif
if (m_hWIGOutFile < 0)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to create/truncate %s - %s", m_szWIGFile, strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}

m_OutBuffIdx = sprintf((char *)m_pOutBuffer, "\"NumSamples\",\"Chrom\",\"Length\",\"Loci Full Coverage\",\"Loci Min 50%% Coverage\",\"Loci No Coverage\",\"Loci Full Concordance\",\"Loci 90%% Concordance\",\"PolyAlleles\",\"MonoAlleles\",\"Loci Full Mixture\"\n");
if (!CUtility::RetryWrites(m_hOutFile, m_pOutBuffer, m_OutBuffIdx))
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Fatal error in RetryWrites()");
	Reset();
	return(eBSFerrFileAccess);
	}
m_OutBuffIdx = 0;

gDiagnostics.DiagOut(eDLInfo, gszProcName, "Process: Completed loading %d input PBA files\n\tStarting concordance analysis for each chromosome", m_NumReadsetIDs);
uint32_t ChromID;
uint32_t CurSampleID;
uint32_t NumNoCoverage;
uint32_t NumFullCoverage;
uint32_t Num50Coverage;
uint32_t NumConcordantLoci;
uint32_t NumNearConcordantLoci;
uint32_t CurConcordantAlleles;
uint32_t CurMonoAlleles;
uint32_t CurPolyAlleles;
uint32_t NumMonoAlleles;
uint32_t NumPolyAlleles;
uint8_t CurSampleAlleles;
uint8_t RefSampleAllele;
uint8_t **pSamplePBAs;
pSamplePBAs = new uint8_t * [m_NumReadsetIDs];
memset(pSamplePBAs, 0, sizeof(uint8_t*)* m_NumReadsetIDs);
InitialiseWIGSpan();
for (ChromID = 1; (int)ChromID < m_NumChromNames; ChromID++)
	{
	if(ChromID != 1)
		CompleteWIGSpan(true);
	LoadChromPBAs(ChromID);			// load PBAs for this chromosome
	uint32_t ChromLen;
	uint32_t NumNoAlleles;
	tsPUChromMetadata* pChromMetadata;
	pChromMetadata = LocateChromMetadataFor(1, ChromID);			// major assumption is that each sample has same length chromosomes as 1st sample
	char* pszChrom = LocateChrom(pChromMetadata->ChromID);
	ChromLen = pChromMetadata->ChromLen;
	NumNoCoverage = 0;
	NumFullCoverage = 0;
	Num50Coverage = 0;
	NumConcordantLoci = 0;
	NumNearConcordantLoci = 0;
	NumMonoAlleles = 0;
	NumPolyAlleles = 0;
	uint8_t* pPBAs;

	// process for intersection of samples with read coverage at each loci
	uint32_t CurLoci = 0;
	for (CurLoci = 0; CurLoci < ChromLen; CurLoci++)
		{
		NumNoAlleles = 0;
		CurConcordantAlleles = 0;
		CurMonoAlleles = 0;
		CurPolyAlleles = 0;
		for (CurSampleID = 1; (int)CurSampleID <= m_NumReadsetIDs; CurSampleID++)
			{
			// cache the PBAs for each sample
			if (CurLoci == 0)
				{
				if ((pSamplePBAs[CurSampleID-1] = LocatePBAfor(CurSampleID, pChromMetadata->ChromID)) == nullptr)
					{
					gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to locate PBAs for '%s' in sample %d",pszChrom, CurSampleID);
					delete[]pSamplePBAs;
					return(eBSFerrChrom);
					}
				}
			pPBAs = pSamplePBAs[CurSampleID-1];
			CurSampleAlleles = pPBAs[CurLoci];
			if (CurSampleID == 1)	// treating as if the reference, concordance will be the proportion of other samples with exactly the same allele
				RefSampleAllele = CurSampleAlleles;
			if (CurSampleAlleles == RefSampleAllele)
				CurConcordantAlleles++;
			if (CurSampleAlleles == 0)
				NumNoAlleles++;
			else
				switch (CurSampleAlleles) {
					case 0xc0: case 0x30: case 0x0c: case 0x03:
					case 0x80: case 0x40: case 0x20: case 0x10: case 0x08: case 0x04: case 0x02: case 0x01:
						CurMonoAlleles++;
						break;
					default:
						CurPolyAlleles++;
						break;
					}

			}

		if (NumNoAlleles == 0)  // will be 0 if all samples at CurLoci had coverage
			{
			NumFullCoverage++;
			if (CurConcordantAlleles == m_NumReadsetIDs) // full concordance if full coverage and all alleles exactly the same
				{
				NumConcordantLoci += 1;
				if (CurPolyAlleles)
					NumPolyAlleles++;
				else
					NumMonoAlleles++;
				}
			else
				if((((int)CurConcordantAlleles * 100) / 90) >= m_NumReadsetIDs) // near concordant if at least 90% exactly the same as the assumed reference
					NumNearConcordantLoci += 1;
			}
		else
			if (NumNoAlleles == m_NumReadsetIDs)
				NumNoCoverage++;
			else
				if (((int)NumNoAlleles * 2) > m_NumReadsetIDs)
					Num50Coverage++;
		AccumWIGCnts(ChromID, CurLoci + 1, 1, m_NumReadsetIDs - NumNoAlleles); // WIG loci start from 1 not 0
		}

	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Process: Samples:%d Chrom:%s Length: %d, loci coverage - full: %d +50%%: %d none: %d, loci concordance - full: %d +90%%: %d, alleles - poly: %d mono: %d", m_NumReadsetIDs, pszChrom, ChromLen, NumFullCoverage, Num50Coverage, NumNoCoverage, NumConcordantLoci, NumNearConcordantLoci, NumPolyAlleles, NumMonoAlleles, NumFullCoverage - NumConcordantLoci);

	m_OutBuffIdx = sprintf((char*)m_pOutBuffer, "%d,\"%s\",%d,%d,%d,%d,%d,%d,%d,%d,%d\n", m_NumReadsetIDs, pszChrom, ChromLen, NumFullCoverage, Num50Coverage, NumNoCoverage, NumConcordantLoci,NumNearConcordantLoci, NumPolyAlleles, NumMonoAlleles, NumFullCoverage- NumConcordantLoci);
	if (!CUtility::RetryWrites(m_hOutFile, m_pOutBuffer, m_OutBuffIdx))
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Fatal error in RetryWrites()");
		Reset();
		return(eBSFerrFileAccess);
		}
	m_OutBuffIdx = 0;
	DeleteAllChromPBAs();
	}
CompleteWIGSpan(true);
if (pSamplePBAs != nullptr)
	delete[]pSamplePBAs;

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

if (m_hWIGOutFile != -1)
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

int
CPBAutils::GenerateWIGConcordance(void)	// generate coverage depth concordance across all sample WIGs
{
	if (m_pOutBuffer == nullptr)
	{
		if ((m_pOutBuffer = new uint8_t[cOutBuffSize]) == nullptr) // an overkill in terms of buffer size!
		{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Memory allocation of %u bytes for output buffering failed - %s", cOutBuffSize, strerror(errno));
			Reset();
			return(eBSFerrMem);
		}
		m_AllocOutBuff = cOutBuffSize;
	}
	m_OutBuffIdx = 0;

	if (m_pszWIGBuff == nullptr)
	{
		if ((m_pszWIGBuff = new uint8_t[cOutBuffSize]) == nullptr) // an overkill in terms of buffer size!
		{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Memory allocation of %u bytes for WIG output buffering failed - %s", cOutBuffSize, strerror(errno));
			Reset();
			return(eBSFerrMem);
		}
		m_AllocWIGBuff = cOutBuffSize;
	}
	m_WIGBuffIdx = 0;

	CUtility::AppendFileNameSuffix(m_szWIGFile, m_szOutFile, (char *)".wig", '.');

#ifdef _WIN32
	m_hOutFile = open(m_szOutFile, (O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC), (_S_IREAD | _S_IWRITE));
#else
	if ((m_hOutFile = open64(m_szOutFile, O_WRONLY | O_CREAT, S_IREAD | S_IWRITE)) != -1)
		if (ftruncate(m_hOutFile, 0) != 0)
		{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to create/truncate %s - %s", m_szOutFile, strerror(errno));
			Reset();
			return(eBSFerrCreateFile);
		}
#endif
	if (m_hOutFile < 0)
	{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to create/truncate %s - %s", m_szOutFile, strerror(errno));
		Reset();
		return(eBSFerrCreateFile);
	}

#ifdef _WIN32
	m_hWIGOutFile = open(m_szWIGFile, (O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC), (_S_IREAD | _S_IWRITE));
#else
	if ((m_hWIGOutFile = open64(m_szWIGFile, O_WRONLY | O_CREAT, S_IREAD | S_IWRITE)) != -1)
		if (ftruncate(m_hWIGOutFile, 0) != 0)
		{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to create/truncate %s - %s", m_szWIGFile, strerror(errno));
			Reset();
			return(eBSFerrCreateFile);
		}
#endif
	if (m_hWIGOutFile < 0)
	{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to create/truncate %s - %s", m_szWIGFile, strerror(errno));
		Reset();
		return(eBSFerrCreateFile);
	}

	m_OutBuffIdx = sprintf((char*)m_pOutBuffer, "\"NumSamples\",\"Chrom\",\"Length\",\"Loci Full Coverage\",\"Loci Min 50%% Coverage\",\"Loci No Coverage\",\"Loci Full Concordance\",\"Loci 90%% Concordance\",\"PolyAlleles\",\"MonoAlleles\",\"Loci Full Mixture\"\n");
	if (!CUtility::RetryWrites(m_hOutFile, m_pOutBuffer, m_OutBuffIdx))
	{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Fatal error in RetryWrites()");
		Reset();
		return(eBSFerrFileAccess);
	}
	m_OutBuffIdx = 0;

	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Process: Completed loading %d input PBA files\n\tStarting concordance analysis for each chromosome", m_NumReadsetIDs);
	uint32_t ChromID;
	uint32_t CurSampleID;
	uint32_t NumNoCoverage;
	uint32_t NumFullCoverage;
	uint32_t Num50Coverage;
	uint32_t NumConcordantLoci;
	uint32_t NumNearConcordantLoci;
	uint32_t CurConcordantAlleles;
	uint32_t CurMonoAlleles;
	uint32_t CurPolyAlleles;
	uint32_t NumMonoAlleles;
	uint32_t NumPolyAlleles;
	uint8_t CurSampleAlleles;
	uint8_t RefSampleAllele;
	uint8_t** pSamplePBAs;
	pSamplePBAs = new uint8_t * [m_NumReadsetIDs];
	memset(pSamplePBAs, 0, sizeof(uint8_t*) * m_NumReadsetIDs);
	InitialiseWIGSpan();
	for (ChromID = 1; (int)ChromID < m_NumChromNames; ChromID++)
	{
		if (ChromID != 1)
			CompleteWIGSpan(true);
		LoadChromPBAs(ChromID);			// load PBAs for this chromosome
		uint32_t ChromLen;
		uint32_t NumNoAlleles;
		tsPUChromMetadata* pChromMetadata;
		pChromMetadata = LocateChromMetadataFor(1, ChromID);			// major assumption is that each sample has same length chromosomes as 1st sample
		char* pszChrom = LocateChrom(pChromMetadata->ChromID);
		ChromLen = pChromMetadata->ChromLen;
		NumNoCoverage = 0;
		NumFullCoverage = 0;
		Num50Coverage = 0;
		NumConcordantLoci = 0;
		NumNearConcordantLoci = 0;
		NumMonoAlleles = 0;
		NumPolyAlleles = 0;
		uint8_t* pPBAs;

		// process for intersection of samples with read coverage at each loci
		uint32_t CurLoci = 0;
		for (CurLoci = 0; CurLoci < ChromLen; CurLoci++)
		{
			NumNoAlleles = 0;
			CurConcordantAlleles = 0;
			CurMonoAlleles = 0;
			CurPolyAlleles = 0;
			for (CurSampleID = 1; (int)CurSampleID <= m_NumReadsetIDs; CurSampleID++)
				{
				if (CurLoci == 0)
					{
					if ((pSamplePBAs[CurSampleID - 1] = LoadPBAChromCoverage(CurSampleID, pChromMetadata->ChromID)) == nullptr)
						{
						gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to locate PBAs for '%s' in sample %d", pszChrom, CurSampleID);
						delete[]pSamplePBAs;
						return(eBSFerrChrom);
						}
					}
				pPBAs = pSamplePBAs[CurSampleID - 1];
				CurSampleAlleles = pPBAs[CurLoci];
				if (CurSampleID == 1)	// treating as if the reference, concordance will be the proportion of other samples with exactly the same allele
					RefSampleAllele = CurSampleAlleles;

// here is where the coverage (scaled 0..254) needs to be segmented into groups relative to the mean coverage over all samples
				if (CurSampleAlleles == RefSampleAllele)
					CurConcordantAlleles++;
				if (CurSampleAlleles == 0)
					NumNoAlleles++;
				else
					switch (CurSampleAlleles) {
					case 0xc0: case 0x30: case 0x0c: case 0x03:
					case 0x80: case 0x40: case 0x20: case 0x10: case 0x08: case 0x04: case 0x02: case 0x01:
						CurMonoAlleles++;
						break;
					default:
						CurPolyAlleles++;
						break;
					}

				}

			if (NumNoAlleles == 0)  // will be 0 if all samples at CurLoci had coverage
			{
				NumFullCoverage++;
				if (CurConcordantAlleles == m_NumReadsetIDs) // full concordance if full coverage and all alleles exactly the same
				{
					NumConcordantLoci += 1;
					if (CurPolyAlleles)
						NumPolyAlleles++;
					else
						NumMonoAlleles++;
				}
				else
					if ((((int)CurConcordantAlleles * 100) / 90) >= m_NumReadsetIDs) // near concordant if at least 90% exactly the same as the assumed reference
						NumNearConcordantLoci += 1;
			}
			else
				if (NumNoAlleles == m_NumReadsetIDs)
					NumNoCoverage++;
				else
					if (((int)NumNoAlleles * 2) > m_NumReadsetIDs)
						Num50Coverage++;
			AccumWIGCnts(ChromID, CurLoci + 1, 1, m_NumReadsetIDs - NumNoAlleles); // WIG loci start from 1 not 0
		}

		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Process: Samples:%d Chrom:%s Length: %d, loci coverage - full: %d +50%%: %d none: %d, loci concordance - full: %d +90%%: %d, alleles - poly: %d mono: %d", m_NumReadsetIDs, pszChrom, ChromLen, NumFullCoverage, Num50Coverage, NumNoCoverage, NumConcordantLoci, NumNearConcordantLoci, NumPolyAlleles, NumMonoAlleles, NumFullCoverage - NumConcordantLoci);

		m_OutBuffIdx = sprintf((char*)m_pOutBuffer, "%d,\"%s\",%d,%d,%d,%d,%d,%d,%d,%d,%d\n", m_NumReadsetIDs, pszChrom, ChromLen, NumFullCoverage, Num50Coverage, NumNoCoverage, NumConcordantLoci, NumNearConcordantLoci, NumPolyAlleles, NumMonoAlleles, NumFullCoverage - NumConcordantLoci);
		if (!CUtility::RetryWrites(m_hOutFile, m_pOutBuffer, m_OutBuffIdx))
		{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Fatal error in RetryWrites()");
			Reset();
			return(eBSFerrFileAccess);
		}
	m_OutBuffIdx = 0;
	DeleteAllChromPBAs();
	}
CompleteWIGSpan(true);
if (pSamplePBAs != nullptr)
	delete[]pSamplePBAs;

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

if (m_hWIGOutFile != -1)
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

// generate coverage (usually, but not restricted to, RNA-seq alignments) segments and output to BED
// intent is that this may allow discovery of novel transcripts and if currently annotated features require extension or are not actually transcribed
int
CPBAutils::GeneratePBASegsBed(double MinCovThres,	// segment coverage initiation if previously not meeting coverage threshold and now at least this proportion of all PBA samples have coverage
							double ContCovThres,	// coverage continuation if coverage stays above this proportional threshold and,
							int32_t MinLen)			// coverage continues for at least this many bp
{
char szBEDFile[_MAX_PATH+1];

if (m_pOutBuffer == nullptr)
	{
	if ((m_pOutBuffer = new uint8_t[cOutBuffSize]) == nullptr) // an overkill in terms of buffer size!
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Memory allocation of %u bytes for output buffering failed - %s", cOutBuffSize, strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
	m_AllocOutBuff = cOutBuffSize;
	}
m_OutBuffIdx = 0;

CUtility::AppendFileNameSuffix(szBEDFile, m_szOutFile, (char *)".bed", '.');
gDiagnostics.DiagOut(eDLInfo, gszProcName, "Generating coverage segments, segment initiated if %.3f of all samples have coverage, coverage stays above %.3f of all samples and extends at least %dbp", MinCovThres, ContCovThres,MinLen);
gDiagnostics.DiagOut(eDLInfo, gszProcName,"Coverage segments will be written to '%s' in BED format",szBEDFile);
#ifdef _WIN32
m_hOutFile = open(szBEDFile, (O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC), (_S_IREAD | _S_IWRITE));
#else
if ((m_hOutFile = open64(szBEDFile, O_WRONLY | O_CREAT, S_IREAD | S_IWRITE)) != -1)
if (ftruncate(m_hOutFile, 0) != 0)
{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to create/truncate '%s' - %s", szBEDFile, strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
}
#endif
if (m_hOutFile < 0)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to create/truncate '%s' - %s", szBEDFile, strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}

if (!CUtility::RetryWrites(m_hOutFile, m_pOutBuffer, m_OutBuffIdx))
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Fatal error in RetryWrites()");
	Reset();
	return(eBSFerrFileAccess);
	}
m_OutBuffIdx = 0;

gDiagnostics.DiagOut(eDLInfo, gszProcName, "Process: Completed loading %d input PBA files\n\tStarting segment discovery for each chromosome", m_NumReadsetIDs);
uint32_t ChromID;
uint32_t ChromLen;
uint32_t CurLoci;
uint32_t CurSampleID;
uint32_t NumbLociSamplesWithCoverage;
double PropLociSamplesWithCoverage;
int32_t CurSegStartLoci;
int32_t CurSegLen;
uint8_t **pSamplePBAs;
uint8_t* pPBAs;
tsPUChromMetadata* pChromMetadata;
pSamplePBAs = new uint8_t * [m_NumReadsetIDs];
uint64_t AllChromsSegsLen;
uint64_t NumAllChromSegs;
uint32_t CurChromSegsLen;
uint32_t CurNumChromSegs;
AllChromsSegsLen = 0;
NumAllChromSegs = 0;
for (ChromID = 1; (int)ChromID < m_NumChromNames; ChromID++)
	{
	LoadChromPBAs(ChromID);			// load PBAs for this chromosome
	pChromMetadata = LocateChromMetadataFor(1, ChromID);			// major assumption is that each sample has same length chromosomes as 1st sample
	char* pszChrom = LocateChrom(pChromMetadata->ChromID);
	gDiagnostics.DiagOut(eDLInfo, gszProcName,"Processing chrom '%s'",pszChrom);

	ChromLen = pChromMetadata->ChromLen;
	memset(pSamplePBAs, 0, sizeof(uint8_t*)* m_NumReadsetIDs);
	for (CurSampleID = 1; (int)CurSampleID <= m_NumReadsetIDs; CurSampleID++)
		{
		if ((pSamplePBAs[CurSampleID-1] = LocatePBAfor(CurSampleID, pChromMetadata->ChromID)) == nullptr)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to locate PBAs for '%s' in sample %d",pszChrom, CurSampleID);
			delete[]pSamplePBAs;
			return(eBSFerrChrom);
			}
		}
	gDiagnostics.DiagOut(eDLInfo, gszProcName,"Processing chrom '%s'. %d sample PBAs loaded",pszChrom,m_NumReadsetIDs);
	NumbLociSamplesWithCoverage = 0;
	CurSegLen = 0;
	CurSegStartLoci = 0;
	CurChromSegsLen = 0;
	CurNumChromSegs = 0;
	// process for intersection of samples with read coverage at each loci
	for (CurLoci = 0; CurLoci < ChromLen; CurLoci++)
		{
		NumbLociSamplesWithCoverage = 0;
		for (CurSampleID = 1; (int)CurSampleID <= m_NumReadsetIDs; CurSampleID++)
			{
			pPBAs = pSamplePBAs[CurSampleID-1];
			if(pPBAs[CurLoci] != 0)
				NumbLociSamplesWithCoverage++;;
			}

		PropLociSamplesWithCoverage = (double)(NumbLociSamplesWithCoverage+1)/m_NumReadsetIDs;
		if(CurSegLen == 0 && PropLociSamplesWithCoverage < MinCovThres)
			continue;
		if (CurLoci == ChromLen -1 ||  PropLociSamplesWithCoverage < ContCovThres)	// dropped below the min?
			{
			if(CurSegLen < MinLen)
				{
				CurSegLen = 0;
				continue;
				}
			// accepting coverage, starts at DeemedCoverageStartLoci and is of length CurSegLen
			CurNumChromSegs++;
			m_OutBuffIdx += sprintf((char *)&m_pOutBuffer[m_OutBuffIdx], "%.80s\t%d\t%d\t%.80s.%d\t0\t+\t0\t%d\t0\t1\t%d,\t0\n", pszChrom, CurSegStartLoci, CurSegLen, pszChrom, CurNumChromSegs, CurSegLen, CurSegLen);
			if((m_OutBuffIdx + 10000) > m_AllocOutBuff)
				{
				if (!CUtility::RetryWrites(m_hOutFile, m_pOutBuffer, m_OutBuffIdx))
					{
					gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Fatal error in RetryWrites()");
					Reset();
					return(eBSFerrFileAccess);
					}
				m_OutBuffIdx = 0;
				}

			CurChromSegsLen+=CurSegLen;
			CurSegLen = 0;
			continue;
			}
		if(CurSegLen ==0)	// if just starting segment coverage then record loci at which it starts
			CurSegStartLoci = CurLoci;
		CurSegLen++;
		}

	NumAllChromSegs += CurNumChromSegs;
	AllChromsSegsLen += CurChromSegsLen;
	gDiagnostics.DiagOut(eDLInfo, gszProcName,"Completed chrom '%s', discovered %u segments totaling %u in length",pszChrom,CurNumChromSegs,CurChromSegsLen);
	DeleteAllChromPBAs();
	}
gDiagnostics.DiagOut(eDLInfo, gszProcName,"Completed processing, over all chroms discovered %zd segments totaling %zd in length",NumAllChromSegs,AllChromsSegsLen);
if (pSamplePBAs != nullptr)
	delete[]pSamplePBAs;
if(m_OutBuffIdx)
	{
	if (!CUtility::RetryWrites(m_hOutFile, m_pOutBuffer, m_OutBuffIdx))
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Fatal error in RetryWrites()");
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

return(eBSFSuccess);
}

bool 
CPBAutils::AcceptROIChromID(int32_t ChromID)	// testing if this chromosome can be mapped to any ROI chromosome
{
char *pszMRAChrom;

if(m_pROIFile == nullptr)			// if no ROIs then all chroms will be accepted
	return(true);

if(ChromID != m_MRA_ChromID)
	{
	pszMRAChrom = LocateChrom(ChromID);
	if((m_MRA_ROIChromID = m_pROIFile->LocateChromIDbyName(pszMRAChrom)) < 1)
		return(false);
	m_MRA_ChromID = ChromID;
	}
return(true);
}

bool 
CPBAutils::AcceptROILoci(int32_t ChromID,	// chromosome to map onto ROI chromosome
					uint32_t Loci)	// testing if this loci is within any ROI on the mapped ChromID
{
char *pszMRAChrom;

if(m_pROIFile == nullptr)			// if no ROIs then all chrom.loci will be accepted
	return(true);

if(ChromID != m_MRA_ChromID)
	{
	pszMRAChrom = LocateChrom(ChromID);
	if((m_MRA_ROIChromID = m_pROIFile->LocateChromIDbyName(pszMRAChrom)) < 1)
		return(false);
	m_MRA_ChromID = ChromID;
	}
return(m_pROIFile->InAnyFeature(m_MRA_ROIChromID,Loci,Loci));
}

int
CPBAutils::LoadPBACoverage(char* pszInWIGFile)   // file containing PBA, file name extension will be replaced with 'coverage.wig' which will be expected to be name of file containing the WIG coverage
{
int32_t ReadsetID;
FILE* pInStream;
char szInWIG[1000];
char* pszInWIG;

// compose WIG file name from PBA file name and check if file can be opened
CUtility::AppendFileNameSuffix(szInWIG, pszInWIGFile, (char*)".covsegs.wig", '.');
pszInWIG = szInWIG;
if ((pInStream = fopen(pszInWIG, "r")) == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadPBACoverage: Unable to open WIG file %s for reading, error: %s", pszInWIG, strerror(errno));
	Reset();
	return(eBSFerrOpnFile);
	}
fclose(pInStream);

// load PBA file, but skip allocation and loading the actual PBAs
if ((ReadsetID = LoadPBAFile(pszInWIGFile, 0, true)) <= 0)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadPBACoverage: Errors loading WIG file '%s'", pszInWIGFile);
	fclose(pInStream);
	Reset();
	return(ReadsetID);
	}
return(ReadsetID);
}


// Note: an assumption is that within the WIG input file ordering is by chromosome ascending
uint8_t* // loaded PBAs for requested chromosome or nullptr if unable to load
CPBAutils::LoadPBAChromCoverage(uint32_t ReadsetID, // loading chromosome coverage for this readset and mapping coverage as though PBAs
	uint32_t ChromID)    // coverage is for this chromosome
{
int Rslt;

char szInWIG[_MAX_PATH];
char* pszInWIG;
char szLineBuff[1000];
tsPUReadsetMetadata* pReadsetMetadata;
tsPUChromMetadata* pChromMetadata;
char* pszReadset;
FILE* pInStream;
uint32_t LineNumb;
uint32_t WIGtype;    // 0: unknown, 1: fixed steps, 2 variable steps
uint32_t CurChromID;
uint32_t Coverage;
uint32_t StartLoci;
uint32_t Span;
uint8_t* pAllele;
uint32_t StartLociOfs;
char* pszChrom;
char szChrom[100];
char* pszInBuff;
uint32_t NumElsParsed;

// compose WIG file name from PBA file name
pReadsetMetadata = &m_Readsets[ReadsetID - 1];
CUtility::AppendFileNameSuffix(szInWIG, pReadsetMetadata->szFileName, (char*)".covsegs.wig", '.');
pszInWIG = szInWIG;

if ((pszReadset = LocateReadset(ReadsetID)) == nullptr)
	return(nullptr);
if (ChromID == 1)
	pReadsetMetadata->NxtFileChromOfs = 0;

if ((pszChrom = LocateChrom(ChromID)) == nullptr)
	return(nullptr);

pChromMetadata = LocateChromMetadataFor(ReadsetID, ChromID);
if (pChromMetadata->pPBAs == nullptr)
	{
	if ((pChromMetadata->pPBAs = AllocPBAs(pChromMetadata->ChromLen)) == nullptr)
		return(nullptr);
	}
memset(pChromMetadata->pPBAs, 0, pChromMetadata->ChromLen);

// WIG file can now be opened and coverage for requested chromosome parsed and mapped as though PBAs
if ((pInStream = fopen(pszInWIG, "r")) == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadPBACoverage: Unable to open WIG file %s for reading, error: %s", pszInWIG, strerror(errno));
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

while (fgets(szLineBuff, sizeof(szLineBuff) - 1, pInStream) != nullptr)
	{
	LineNumb++;
	pszInBuff = CUtility::ReduceWhitespace(szLineBuff);
	if (pszInBuff == nullptr || *pszInBuff == '\0' || *pszInBuff == '\n')
		continue;

	if (pszInBuff[0] == 'v' && !strnicmp(pszInBuff, "variableStep", 12))
		{
		NumElsParsed = sscanf(pszInBuff, "variableStep chrom=%s span=%d", szChrom, &Span);
		if (NumElsParsed != 2)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadPBACoverage: Errors parsing WIG line %d - \"%s\" - in WIG file '%s'", LineNumb, pszInBuff, pszInWIG);
			Rslt = -1;
			break;
			}

		// check for valid chrom, a missing chrom relative to PBA could be that the user has specified a limit on number of chromosomes to process so can't treat as a fatal error
		if ((CurChromID = LocateChrom(szChrom)) < 1)
			{
			CurChromID = 0;
			continue;
			}

		if (CurChromID > ChromID) // onto coverage for next chromosome?
			break;

		if (Span < 1 || Span > pChromMetadata->ChromLen)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadPBACoverage: Span must be between 1..%d \"%s\" at line %d - \"%s\" - for chromosome '%s' in WIG file '%s'", pChromMetadata->ChromLen, LineNumb, pszInBuff, pszChrom, pszInWIG);
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
		if (pszInBuff[0] == 'f' && !strnicmp(pszInBuff, "fixedStep", 10))
			{
			NumElsParsed = sscanf(pszInBuff, "fixedStep chrom=%s start=%d step=%d", szChrom, &StartLoci, &Span);
			if (NumElsParsed != 3)
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadPBACoverage: Errors parsing line %d - \"%s\" - in WIG file '%s'", LineNumb, pszInBuff, pszInWIG);
				Rslt = -1;
				break;
				}

			// check for valid chrom, a missing chrom relative to PBA could be that the user has specified a limit on number of chromosomes to process so can't treat as a fatal error
			if ((CurChromID = LocateChrom(szChrom)) < 1)
				{
				CurChromID = 0;
				continue;
				}

			if (CurChromID > ChromID) // onto coverage for next chromosome?
				break;

			if (StartLoci < 1 || (StartLoci + Span - 1) >  pChromMetadata->ChromLen)
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadPBACoverage: StartLoci must be between 1..%d, Span %d extending past end of chromosome \"%s\" at line %d - \"%s\" - in WIG file '%s'", pChromMetadata->ChromLen, Span, pszChrom, LineNumb, pszInBuff, pszInWIG);
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
	if (CurChromID == 0) // 0 if yet to parse out coverage for targeted chromosome
		{
#if _WIN32
		pReadsetMetadata->NxtFileChromOfs = _ftelli64(pInStream);
#else
		pReadsetMetadata->NxtFileChromOfs = ftello64(pInStream);
#endif
		continue;
		}

	switch (WIGtype) {
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
	if (StartLoci < 1 || StartLoci > pChromMetadata->ChromLen) // ensure not writing past end of chromosome
		continue;
	pAllele = &pChromMetadata->pPBAs[StartLoci - 1];
	for (StartLociOfs = 0; StartLociOfs < Span; StartLociOfs++, StartLoci++, pAllele++)
		{
		if (StartLoci > pChromMetadata->ChromLen) // ensure not writing past end of chromosome
			break;
		*pAllele = min((int)log2(Coverage) * 10, 0x0fe);
		}
}
fclose(pInStream);
if (Rslt < 0)
	{
	DeleteSampleChromPBAs(ReadsetID, CurChromID);
	pChromMetadata->pPBAs = nullptr;
	}
return(pChromMetadata->pPBAs);
}

int
CPBAutils::LoadChromPBAs(uint32_t ChromID,			// load PBAs for this chromosome
			  uint32_t StartSampleID,				// chromosome PBAs for this sample through 
			  uint32_t EndSampleID)					// to this sample inclusive
{
int Rslt;
uint8_t* pPBAs;
uint32_t CurSampleID;
if (EndSampleID == 0)
	EndSampleID = m_NumReadsetIDs;

for (CurSampleID = StartSampleID; CurSampleID <= EndSampleID; CurSampleID++)
	{
	if ((pPBAs = LoadSampleChromPBAs(CurSampleID, ChromID)) == nullptr)
		break;
}
if (CurSampleID <= EndSampleID)
	Rslt = -1;
else
	Rslt = eBSFSuccess;
return(Rslt);
}

uint8_t
CPBAutils::LocateReadsetChromLociAlleles(uint32_t ReadsetID,	// return alleles for this readset 
	uint32_t ChromID,		// on this chromosome
	uint32_t Loci)		// at this loci
{
	uint8_t* pPBA;
	if ((pPBA = LocatePBAfor(ReadsetID, ChromID)) == nullptr)
		return(0);

	return(pPBA[Loci]);
}


// PBA at a single loci is packed into a single 8bit byte: Allele A in bits 7.6, C in bits 5.4, G in bits 3.2, T in bits 1.0
// following allele scoring thresholds are extracted from KAligner.h, you must ensure that these thresholds are updated each time KAligner is updated
//const double cScorePBA3MinProp = 0.75;		// score PBA as 3 if allele proportion of all counts is >= this threshold and coverage is >= 5
//const double cScorePBA2MinProp = 0.35;		// score PBA as 2 if allele proportion of all counts is >= this threshold and coverage is >= 5
//const double cScorePBA1MinProp = 0.20;		// score PBA as 1 if allele proportion of all counts is >= this threshold and coverage is >= 5
// when coverage is less than 5 (usually the case with progeny reads) and thus confidence in alleles is reduced then scores are reduced
//const double cScorePBA2MinLCProp = 0.70;		// score PBA as 2 if allele proportion of all counts is >= this threshold and coverage is < 5
//const double cScorePBA1MinLCProp = 0.30;		// score PBA as 1 if allele proportion of all counts is >= this threshold and coverage is < 5

uint8_t
CPBAutils::ConsensusHaploid(uint8_t PBA) // returns consensus haploid allele as PBA from a possibly polyallelic PBA
{
if(PBA == 0)	// no alleles
	return(0);

// process for unambiguous A,C,G or T - values 3
if((PBA & 0x0c0) == 0x0c0)
	return(0x0c0);
if((PBA & 0x030) == 0x030)
	return(0x030);
if((PBA & 0x0c)==0x0c)
	return(0x0c);
if((PBA & 0x03)==0x03)
	return(0x03);

// try for allele with value 2
if((PBA & 0x0c0) == 0x080)
	return(0x0c0);
if((PBA & 0x030) == 0x020)
	return(0x030);
if((PBA & 0x0c)==0x08)
	return(0x0c);
if((PBA & 0x03)==0x02)
	return(0x03);

// bottom of barrel, alleles are low coverage
// try for allele with value 1
if(PBA & 0x0c0)
	return(0x0c0);
if(PBA & 0x030)
	return(0x030);
if(PBA & 0x0c)
	return(0x0c);
return(0x03); // PBA contained at least 1 allele so to have reached here then must have been T with value 1
}

uint8_t
CPBAutils::ConsensusDiploid(uint8_t PBA) // returns consensus diploid alleles as PBA from a possibly polyallelic PBA
{
uint8_t Mono;
uint8_t Di;
if(PBA == 0)
	return(0);
Mono = ConsensusHaploid(PBA);
if((PBA & Mono) == Mono) // if top value was 3 then there was at least 5 reads and at least 75% of these were this allele so accept as only a single allele present
	return(Mono);

if((PBA & ~Mono) == 0)	// any other alleles present?
	return(Mono);		// just a single allele

Di = ConsensusHaploid(PBA & ~Mono);
return(Mono | Di);
}


etSeqBase
CPBAutils::Alleles2Base(uint8_t PBA) // identify and return consensus base from a possibly diallelic PBA
{
etSeqBase AltBases[2];
m_NumPBAbases++;

switch (PBA) {
	case 0xc0:
	case 0x80:
	case 0x40:
	case 0x90:
	case 0x84:
	case 0x81:
		m_NumMonoAllelic++;
		return(eBaseA);

	case 0x30:
	case 0x20:
	case 0x10:
	case 0x50:
	case 0x24:
	case 0x21:
		m_NumMonoAllelic++;
		return(eBaseC);

	case 0x0c:
	case 0x08:
	case 0x04:
	case 0x18:
	case 0x48:
	case 0x09:
		m_NumMonoAllelic++;
		return(eBaseG);

	case 0x03:
	case 0x02:
	case 0x01:
	case 0x42:
	case 0x12:
	case 0x0A:
		m_NumMonoAllelic++;
		return(eBaseT);

	case 0xA0:	// eBaseA or eBaseC
	case 0x60:
		AltBases[0] = eBaseA;
		AltBases[1] = eBaseC;
		break;

	case 0x88:  // eBaseA or eBaseG
	case 0x44:
		AltBases[0] = eBaseA;
		AltBases[1] = eBaseG;
		break;

	case 0x82:  // eBaseA or eBaseT
	case 0x41:
		AltBases[0] = eBaseA;
		AltBases[1] = eBaseT;
		break;

	case 0x28:  // eBaseC or eBaseG
	case 0x14:
		AltBases[0] = eBaseC;
		AltBases[1] = eBaseG;
		break;

	case 0x22:  // eBaseC or eBaseT
	case 0x11:
		AltBases[0] = eBaseC;
		AltBases[1] = eBaseT;
		break;

	case 0x06:	// eBaseG or eBaseT
	case 0x05:
		AltBases[0] = eBaseG;
		AltBases[1] = eBaseT;
		break;

	default: // can't discriminate with any confidence
		m_NumNonAllelic++;
		return(eBaseN);
	}
m_NumDiAllelic++;
return(AltBases[m_NumDiAllelic & 0x01]); // not using random selection as selection needs to be deterministic - reproducible between runs
}

char
CPBAutils::PBAFastaBase(uint8_t PBA)  	// returns Fasta base char - 'a','c','g','t' or 'n' - as being the consensus allele in PBA
{
char FastaChar;
etSeqBase Base;

Base = Alleles2Base(PBA);
switch (Base) {
	case eBaseA:
		FastaChar = 'a';
		break;
	case eBaseC:
		FastaChar = 'c';
		break;
	case eBaseG:
		FastaChar = 'g';
		break;
	case eBaseT:
		FastaChar = 't';
		break;
	default:
		FastaChar = 'N';
	}
return(FastaChar);
}


int32_t 
CPBAutils::GenAllelicVCF(int32_t RefPBAID,	// reference PBAs identifier
						int32_t AllelicPBAID)	// PBAs with allelic differences
{

int32_t Rslt = 0;
etSeqBase RefBase;
uint8_t* pRefPBAs;
uint8_t* pAllelicPBAs;
int32_t ChromID;
uint32_t ChromIdx;
uint32_t ChromSize;
uint32_t Loci;
uint8_t RefAlleles;
uint8_t AAleles;
uint32_t NumChromDiffs;
uint32_t TotChromsDiffs;
uint32_t NumChromNA;
uint32_t NumChromDiffNA;
uint32_t NumChromMatches;
char *pszChromName;

tsPUReadsetMetadata* pReadsetMetadata;
tsPUChromMetadata* pChromMetadata;

#ifdef _WIN32
m_hOutFile = open(m_szOutFile, (O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC), (_S_IREAD | _S_IWRITE));
#else
if ((m_hOutFile = open64(m_szOutFile, O_WRONLY | O_CREAT, S_IREAD | S_IWRITE)) != -1)
if (ftruncate(m_hOutFile, 0) != 0)
{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenAllelicVCF: Unable to create/truncate %s - %s", m_szOutFile, strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
}
#endif
if (m_hOutFile < 0)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenAllelicVCF: Unable to create/truncate %s - %s", m_szOutFile, strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}

if (m_pOutBuffer == nullptr)
	{
	if ((m_pOutBuffer = new uint8_t[cOutBuffSize]) == nullptr) // an overkill in terms of buffer size!
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenAllelicVCF: Memory allocation of %u bytes for output buffering failed - %s", cOutBuffSize, strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
	m_AllocOutBuff = cOutBuffSize;
	}
m_OutBuffIdx = 0;

pReadsetMetadata = &m_Readsets[RefPBAID - 1];

m_OutBuffIdx = sprintf((char *)m_pOutBuffer,"##fileformat=VCFv4.1\n##source=pbautils%s\n##reference=%s\n##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n",
				  kit4bversion, (char *)pReadsetMetadata->szRefAssemblyID);
m_OutBuffIdx += sprintf((char *)&m_pOutBuffer[m_OutBuffIdx],"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");

if (!CUtility::RetryWrites(m_hOutFile, m_pOutBuffer, m_OutBuffIdx))
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenAllelicVCF: Fatal error in RetryWrites()");
	Reset();
	return(eBSFerrFileAccess);
	}
m_OutBuffIdx = 0;

// iterate over chromosomes
TotChromsDiffs = 0;
for(ChromIdx = 1; ChromIdx <= pReadsetMetadata->NumChroms; ChromIdx++)
	{
	pChromMetadata = LocateChromMetadataFor(RefPBAID, ChromIdx);
	ChromID = pChromMetadata->ChromID;
	if(!AcceptROIChromID(ChromID))
		continue;
	ChromSize = pChromMetadata->ChromLen;

	if((pRefPBAs = LoadSampleChromPBAs(RefPBAID, ChromID)) == nullptr)
		break;

	if((pAllelicPBAs = LoadSampleChromPBAs(AllelicPBAID, ChromID)) == nullptr)
		break;
	pszChromName = LocateChrom(ChromID);
	NumChromDiffs = 0;
	NumChromNA = 0;
	NumChromDiffNA = 0;
	NumChromMatches = 0;
	for(Loci = 0; Loci < ChromSize; Loci++,pRefPBAs++,pAllelicPBAs++)
		{
		if(!AcceptROILoci(ChromID,Loci))
			continue;
		RefAlleles = ConsensusHaploid(*pRefPBAs); // processing is dependent on the reference PBA being a consensus PBA, force to be consensus!
		AAleles = *pAllelicPBAs;
		if(RefAlleles == 0 || AAleles == 0) // need both PBAs to have alleles
			{
			NumChromNA++;
			if(RefAlleles || AAleles) // check if one was aligned in which case there was a difference in coverage
				NumChromDiffNA++;
			continue;
			}

		if((RefAlleles & 0x0aa) == (AAleles & 0x0aa)) // matching as dirac or major alleles treated as if exactly matching
			{
			NumChromMatches++;
			continue;
			}

		// AllelicPBAID has at least one allele which differs from the reference
		NumChromDiffs++;
		TotChromsDiffs++;
		// generate VCF...
		// determine which alleles differ
		int BaseIdx;
		int AltsIdx;
		int SumFreqs;
		char szALTs[100];
		char szALTsFreq[100];
		uint32_t AlleleMsk;
		int32_t Shl;
		SumFreqs = 0;
		AltsIdx = 0;
		AlleleMsk = 0x0c0;
		Shl = 6;
		for(BaseIdx = 0; AAleles != 0 && BaseIdx < 4; BaseIdx++, AlleleMsk >>= 2, Shl-=2)
			{
			if(!(AAleles & AlleleMsk))			// if no allele then try next allele
				continue;

			if (AltsIdx > 0)
				{
				szALTs[AltsIdx] = ',';
				szALTsFreq[AltsIdx++] = ',';
				}
			szALTs[AltsIdx] = CSeqTrans::MapBase2Ascii(BaseIdx);
			switch ((AAleles & AlleleMsk) >> Shl) {
				case 1:
					szALTsFreq[AltsIdx] = '1';
					SumFreqs += 1;
					break;
				case 2:
					szALTsFreq[AltsIdx] = '4';
					SumFreqs += 2;
					break;
				case 3:
					szALTsFreq[AltsIdx] = '5';
					SumFreqs += 3;
					break;
				}
			AltsIdx++;
			szALTs[AltsIdx] = '\0';
			szALTsFreq[AltsIdx] = '\0';
			AAleles &= ~AlleleMsk;
			}
		
		AlleleMsk = 0xc0;
		for(RefBase = 0; RefBase < 4; RefBase++,AlleleMsk >>= 2)
			{
			if(AlleleMsk & RefAlleles)
				break;
			}

		// for reasons I do not understand - others have complained! - IGV does not accept GVFs with multiple alleles if one matches the reference base.
		// but is happy if the reference base is set to be indeterminate 'N'. Go figure ....
		m_OutBuffIdx += sprintf((char *)&m_pOutBuffer[m_OutBuffIdx], "%s\t%u\tSNP%d\t%c\t%s\t%d\tPASS\tAF=%s;DP=%d\n",
					pszChromName, Loci + 1, TotChromsDiffs, 'N', // CSeqTrans::MapBase2Ascii(RefBase),
					szALTs, 100, szALTsFreq, 
					SumFreqs);		// have no reliable coverage depth so simply default to sum of frequencies

		if(m_OutBuffIdx+1000 > m_AllocOutBuff)
			{
			if (!CUtility::RetryWrites(m_hOutFile, m_pOutBuffer, m_OutBuffIdx))
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenAllelicVCF: Fatal error in RetryWrites()");
				Reset();
				return(eBSFerrFileAccess);
				}
			m_OutBuffIdx = 0;
			}
		}
	if(m_OutBuffIdx)
		{
		if (!CUtility::RetryWrites(m_hOutFile, m_pOutBuffer, m_OutBuffIdx))
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenAllelicVCF: Fatal error in RetryWrites()");
			Reset();
			return(eBSFerrFileAccess);
			}
		m_OutBuffIdx = 0;
		}
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenAllelicVCF: '%s' - %dbp : homo alleles: %d, hetro alleles: %d, non-aligned: %d, differential non-aligned: %d", pszChromName, ChromSize, NumChromMatches, NumChromDiffs, NumChromNA,NumChromDiffNA);
	DeleteSampleChromPBAs(RefPBAID, ChromID);
	DeleteSampleChromPBAs(AllelicPBAID, ChromID);
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


bool
CPBAutils::ReportGTSample(int32_t SampleID) // returns true if this samples genotype is to be reported
{
if(m_pReportGTSample == NULL || (m_pReportGTSample[SampleID] & 0x01))
	return(true);
return(false);
}

int32_t 
CPBAutils::GenGenotypeVCF(ePBAuMode PMode, // will be either ePBAu2GVCF allelic genotype VCF or ePBAu2DVCF deletion genotype VCF
						int32_t RefSampleID,	// reference identifier
						int32_t NumSamples,		// number of samples to genotype relative to reference
						double GTPropNAThres,		// when genotyping VCF then proportion of non-aligned calls over all samples must be < this threshold
						double GTPropHetThres)		// when genotyping VCF then proportion of heterozygous calls must be >= this threshold	
{
int32_t Rslt = 0;
etSeqBase RefBase;
uint8_t* pRefPBAs;
uint8_t *pAllelicPBAs;
uint8_t **ppChromPBAs;
char *pszSampleHaps;
int32_t SampleHapsIdx;
int32_t ChromID;
uint32_t ChromIdx;
uint32_t ChromSize;
uint32_t Loci;
uint8_t RefAlleles;
uint8_t SampleAlleles;
uint8_t AltAlleles;
int SampleID;
char *pszChromName;
char szALTs[2048];
uint32_t TotGTIDs;
uint32_t ChromGTIDs;
bool bRptThisLoci = false;
uint32_t LastReportedLoci = 0;
uint32_t ReportUpToLoci = 0;
uint32_t AlleleCnts[256]; // counts of all possible allele combinations

tsPUReadsetMetadata* pReadsetMetadata;
tsPUChromMetadata* pChromMetadata;

#ifdef _WIN32
m_hOutFile = open(m_szOutFile, (O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC), (_S_IREAD | _S_IWRITE));
#else
if ((m_hOutFile = open64(m_szOutFile, O_WRONLY | O_CREAT, S_IREAD | S_IWRITE)) != -1)
if (ftruncate(m_hOutFile, 0) != 0)
{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenGenotypeVCF: Unable to create/truncate %s - %s", m_szOutFile, strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
}
#endif
if (m_hOutFile < 0)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenGenotypeVCF: Unable to create/truncate %s - %s", m_szOutFile, strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}

if (m_pOutBuffer == nullptr)
	{
	if ((m_pOutBuffer = new uint8_t[cOutBuffSize]) == nullptr) // an overkill in terms of buffer size!
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenGenotypeVCF: Memory allocation of %u bytes for output buffering failed - %s", cOutBuffSize, strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
	m_AllocOutBuff = cOutBuffSize;
	}
m_OutBuffIdx = 0;

pReadsetMetadata = &m_Readsets[RefSampleID - 1];

m_OutBuffIdx = sprintf((char *)m_pOutBuffer,"##fileformat=VCFv4.1\n##source=pbautils%s\n##reference=%s\n##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Unphased genotypes\">\n",
				  kit4bversion, (char *)pReadsetMetadata->szRefAssemblyID);
m_OutBuffIdx += sprintf((char *)&m_pOutBuffer[m_OutBuffIdx],"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
for(SampleID = RefSampleID+1; SampleID <= NumSamples+1; SampleID++)
	{
	if(!ReportGTSample(SampleID))
		continue;
	char *pszSampleName = LocateReadset(SampleID);
	m_OutBuffIdx += sprintf((char *)&m_pOutBuffer[m_OutBuffIdx],"\t%s",pszSampleName);
	}
m_OutBuffIdx += sprintf((char *)&m_pOutBuffer[m_OutBuffIdx],"\n");

if (!CUtility::RetryWrites(m_hOutFile, m_pOutBuffer, m_OutBuffIdx))
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenGenotypeVCF: Fatal error in RetryWrites()");
	Reset();
	return(eBSFerrFileAccess);
	}
m_OutBuffIdx = 0;

ppChromPBAs = new uint8_t * [NumSamples];
pszSampleHaps = new char [NumSamples * 10];

// iterate over chromosomes
ChromGTIDs = 0;
TotGTIDs = 0;
for(ChromIdx = 1; ChromIdx <= pReadsetMetadata->NumChroms; ChromIdx++)
	{
	ChromGTIDs = 0;
	pChromMetadata = LocateChromMetadataFor(RefSampleID, ChromIdx);
	ChromID = pChromMetadata->ChromID;
	if(!AcceptROIChromID(ChromID))
		continue;
	pszChromName = LocateChrom(ChromID);
	ChromSize = pChromMetadata->ChromLen;
	
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenGenotypeVCF: Processing chrom: '%s'",pszChromName);
	if((pRefPBAs = LoadSampleChromPBAs(RefSampleID, ChromID)) == nullptr)
		break;

	uint8_t *pPBA = pRefPBAs;
		// normalise reference alleles to be consensus monoallelic - VCF assumes that alignments were to a consensus reference
	for(Loci = 0; Loci < ChromSize; Loci++,pPBA++)
		*pPBA = ConsensusHaploid(*pPBA);


	// load chrom PBAs for each sample
	for(SampleID = RefSampleID+1; SampleID <= NumSamples+1; SampleID++)
		{
		if((ppChromPBAs[SampleID-2] = LoadSampleChromPBAs(SampleID, ChromID)) == nullptr)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Fatal error loading chrom PBAs");
			if(ppChromPBAs != nullptr)
				delete []ppChromPBAs;
			if(pszSampleHaps != nullptr)
				delete []pszSampleHaps;
			Reset();
			return(eBSFerrChrom);
			}

			// normalise these to be ConsensusDiploid
		pPBA = ppChromPBAs[SampleID-2];
		for(Loci = 0; Loci < ChromSize; Loci++,pPBA++)
			*pPBA = ConsensusDiploid(*pPBA);
		}


	// both reference and all samples have been normalised over the current chromosomes full length
	for(Loci = 0; Loci < ChromSize; Loci++,pRefPBAs++)
		{
		if(!AcceptROILoci(ChromID,Loci))
			{
			// flush out any genotypes yet to be reported
			if(m_OutBuffIdx)
				{
				if (!CUtility::RetryWrites(m_hOutFile, m_pOutBuffer, m_OutBuffIdx))
					{
					gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenGenotypeVCF: Fatal error in RetryWrites()");
					DeleteAllChromPBAs();
					if(ppChromPBAs != nullptr)
						delete []ppChromPBAs;
					if(pszSampleHaps != nullptr)
						delete []pszSampleHaps;
					Reset();
					return(eBSFerrFileAccess);
					}
				m_OutBuffIdx = 0;
				}
			continue;
			}

		if((RefAlleles = *pRefPBAs)==0) // reference must have coverage before comparing with samples
			continue;

		// if all samples same then no differential alleles at this loci
		uint8_t PrevSampleAlleles = 0;
		uint32_t NumDiffAlleles = 0;
		uint32_t NumMissing = 0;

		memset(AlleleCnts,0,sizeof(AlleleCnts));
		AltAlleles = 0;
		for(SampleID = RefSampleID+1; SampleID <= NumSamples+1; SampleID++)
			{
			pAllelicPBAs = ppChromPBAs[SampleID-2];
			SampleAlleles = pAllelicPBAs[Loci];
			AlleleCnts[SampleAlleles]++;
			if(SampleAlleles != 0)
				AltAlleles |= SampleAlleles;
			}



		bRptThisLoci = false;
		bool bCoverageInRange = false;
		double CovProportion = 1.0 - (double)AlleleCnts[0] / NumSamples;
		if(PMode == ePBAu2DVCF)
			{
			if(CovProportion >= GTPropNAThres)			// deletion processing needs to have min coverage of at least GTPropNAThres
				bCoverageInRange = true;
			}
		else			// else ePBAu2GVCF
			{
			if(CovProportion >= (1.0 - GTPropNAThres)) //  if only processing for alleles then coverage must be at least (1.0 - GTPropNAThres)
				bCoverageInRange = true;
			}

		if(bCoverageInRange) // bCoverageInRange will be true if meeting requirements for proportions of aligned samples at current loci
			{
			// need highest frequency and sum of other frequencies to determine proportion of heterozygosity between samples and then compare against GTPropHetThres
			// if PMode == ePBAu2DVCF then missing coverage is treated as if allelic variation!
			uint32_t HighAlleleCnt;
			uint32_t LowAlleleCnts;
			uint32_t CntsIdx;

			HighAlleleCnt = 0;
			LowAlleleCnts = 0;
			for(CntsIdx = PMode == ePBAu2GVCF ? 1 : 0; CntsIdx <= 255; CntsIdx++) // ePBAu2GVCF excludes non-aligned samples (assumed deletions)
				{
				if(AlleleCnts[CntsIdx] > HighAlleleCnt) // new highest allelic count?
					{
					LowAlleleCnts += HighAlleleCnt;		// previous highest count is back summed into the lowest 
					HighAlleleCnt = AlleleCnts[CntsIdx];
					}
				else
					LowAlleleCnts += AlleleCnts[CntsIdx];
				}

			// now know highest frequency allelic count and sum of all others so can determine if sufficient heterozygosity
			if(PMode == ePBAu2GVCF) // ePBAu2GVCF excludes deletions, only using alleles when determining heterozygosity
				{
				if ((double)LowAlleleCnts / (NumSamples - AlleleCnts[0]) >= GTPropHetThres)	// if very little heterozygosity in samples excluding non-aligned then not forced to report
					bRptThisLoci = true;
				}
			else  // ePBAu2DVCF heterozygosity includes deletions as well as allelic
				{
				if ((double)LowAlleleCnts / NumSamples >= GTPropHetThres)	// if very little heterozygosity in samples including non-aligned then not forced to report
					bRptThisLoci = true;;
				}
			}

		if(!bRptThisLoci && Loci > ReportUpToLoci)
			{
			ReportUpToLoci = 0;
			continue;
			}

		uint32_t PrevToReport;
		if(bRptThisLoci)  // if some samples have alleles which differ from the reference and/or multiple samples have no alignments at the current loci and will be reported
			{			  // report 3bp loci up/down from this loci which has been identified as needed to be reported
			if(ReportUpToLoci == 0)
				{
				PrevToReport = min((uint32_t)4,Loci - LastReportedLoci); // ensure that a max of 3 previous unreported loci are also being reported
				if(PrevToReport > 1)
					{
					Loci -= PrevToReport;
					pRefPBAs -= PrevToReport;
					ReportUpToLoci = Loci+3;
					continue;
					}
				}
			ReportUpToLoci = Loci+3;  // ensure that will always report next 10 even if no heterozygosity or unaligned samples in these 10 loci
			}

		LastReportedLoci = Loci;
		AltAlleles &= ~RefAlleles;				// remove consensus ref allele from alternative alleles
		SampleHapsIdx = 0;
		pszSampleHaps[0] = '\0';
		for(SampleID = RefSampleID+1; SampleID <= NumSamples+1; SampleID++)
			{
			if(!ReportGTSample(SampleID))
				continue;
			pAllelicPBAs = ppChromPBAs[SampleID-2];
			SampleAlleles = pAllelicPBAs[Loci];

			// reference had coverage, if sample has no coverage then sample has a deletion - or just had no sequenced read covering that loci!!
			if(SampleAlleles == 0) 
				{
				SampleHapsIdx += sprintf(&pszSampleHaps[SampleHapsIdx],"\t./.");
				continue;
				}

			if(RefAlleles == SampleAlleles) // does sample exactly match reference?
				{
				SampleHapsIdx += sprintf(&pszSampleHaps[SampleHapsIdx],"\t0/0");
				continue;
				}

			// sample has at least one allele which doesn't match the reference
			uint32_t Msk;
			uint32_t AltIdx1;
			int32_t AltIdx2;
			AltIdx1 = 0;
			Msk = 0x0c0;
			if(RefAlleles & SampleAlleles)	// does at least one of the sample allele(s) match the reference?
				{
				do
					{
					if(Msk & RefAlleles)
						continue;
					if(!(Msk & AltAlleles))
						continue;
					AltIdx1++;
					if(Msk & SampleAlleles)
						{
						SampleHapsIdx += sprintf(&pszSampleHaps[SampleHapsIdx],"\t0/%d",AltIdx1);
						break;
						}
					}
				while(Msk >>= 2);
				continue;
				}

			// sample alleles do not include the reference allele
			AltIdx1 = -1;
			AltIdx2   = 0;
			Msk = 0x0c0;
			do
				{
				if(!(Msk & AltAlleles))
					continue;
				AltIdx2++;
				if(Msk & SampleAlleles)
					{
					if(AltIdx1 == -1)
						AltIdx1 = AltIdx2;
					SampleAlleles &= ~Msk;
					}
				}
			while((SampleAlleles != 0) && (Msk >>= 2));
			SampleHapsIdx += sprintf(&pszSampleHaps[SampleHapsIdx],"\t%d/%d",AltIdx1,AltIdx2);
			continue;
			}

		uint8_t AltBase;
		uint8_t AltsMsk = 0xc0;
		uint32_t ALTsIdx = 0;
		if(AltAlleles == 0)
			ALTsIdx += sprintf(&szALTs[ALTsIdx],"N");
		else
			{
			for(AltBase = 0; AltBase < 4; AltBase++,AltsMsk >>= 2)
				{
				if(AltsMsk & AltAlleles)
					{
					if(ALTsIdx != 0)
						ALTsIdx += sprintf(&szALTs[ALTsIdx],",");
					ALTsIdx += sprintf(&szALTs[ALTsIdx],"%c",CSeqTrans::MapBase2Ascii(AltBase));
					}
				}
			}


		AltsMsk = 0xc0;
		for(RefBase = 0; RefBase < 4; RefBase++,AltsMsk >>= 2)
			{
			if(AltsMsk & RefAlleles)
				break;
			}

		TotGTIDs += 1;
		ChromGTIDs += 1;
		m_OutBuffIdx += sprintf((char *)&m_pOutBuffer[m_OutBuffIdx], "%s\t%u\tGT%u\t%c\t%s\t.\tPASS\t.\tGT%s\n",
						pszChromName, Loci + 1, TotGTIDs, CSeqTrans::MapBase2Ascii(RefBase),
						szALTs,pszSampleHaps);

		// ensure no buffer overflows
		// every sample requires 4 chars per loci plus an overhead (allowing 500 per loci for fixed fields)
		if((m_OutBuffIdx + 500 + (NumSamples * 4)) >= m_AllocOutBuff)
			{
			if (!CUtility::RetryWrites(m_hOutFile, m_pOutBuffer, m_OutBuffIdx))
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenGenotypeVCF: Fatal error in RetryWrites()");
				DeleteAllChromPBAs();
				if(ppChromPBAs != nullptr)
					delete []ppChromPBAs;
				if(pszSampleHaps != nullptr)
					delete []pszSampleHaps;
				Reset();
				return(eBSFerrFileAccess);
				}
			m_OutBuffIdx = 0;
			}
		}
	if(m_OutBuffIdx)
		{
		if (!CUtility::RetryWrites(m_hOutFile, m_pOutBuffer, m_OutBuffIdx))
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenGenotypeVCF: Fatal error in RetryWrites()");
			DeleteAllChromPBAs();
			if(ppChromPBAs != nullptr)
				delete []ppChromPBAs;
			if(pszSampleHaps != nullptr)
				delete []pszSampleHaps;
			Reset();
			return(eBSFerrFileAccess);
			}
		m_OutBuffIdx = 0;
		}
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenGenotypeVCF: '%s' (%dbp) : reported GTs - chrom: %u, total: %u",pszChromName,ChromSize,ChromGTIDs,TotGTIDs);
	DeleteAllChromPBAs();
	}
gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenGenotypeVCF: Over all chromosomes there are a total of %u reported GTs",TotGTIDs);
if(ppChromPBAs != nullptr)
	delete []ppChromPBAs;
if(pszSampleHaps != nullptr)
	delete []pszSampleHaps;
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

int
CPBAutils::PBA2Fasta(uint32_t ReadSetID)
{
uint8_t Alleles;
uint32_t BaseIdx;

int LineLen = 0;
char FastaBase;
uint8_t* pPBAs;
int32_t ChromID;
uint32_t ChromIdx;

tsPUReadsetMetadata* pReadsetMetadata;
tsPUChromMetadata* pChromMetadata;

m_NumPBAbases = 0;
m_NumDiAllelic = 0;
m_NumMonoAllelic = 0;
m_NumNonAllelic = 0;


if (m_pOutBuffer == nullptr)
	{
	if ((m_pOutBuffer = new uint8_t[cOutBuffSize]) == nullptr)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "PBA2Fasta: Memory allocation of %u bytes for output buffering failed - %s", cOutBuffSize, strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
	m_AllocOutBuff = cOutBuffSize;
	}
m_OutBuffIdx = 0;

pReadsetMetadata = &m_Readsets[ReadSetID - 1];
m_OutBuffIdx = 0;

#ifdef _WIN32
m_hOutFile = open(m_szOutFile, (O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC), (_S_IREAD | _S_IWRITE));
#else
if ((m_hOutFile = open64(m_szOutFile, O_WRONLY | O_CREAT, S_IREAD | S_IWRITE)) != -1)
	if (ftruncate(m_hOutFile, 0) != 0)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "PBA2Fasta: Unable to create/truncate %s - %s", m_szOutFile, strerror(errno));
		Reset();
		return(eBSFerrCreateFile);
		}
#endif
if (m_hOutFile < 0)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "PBA2Fasta: Unable to create/truncate %s - %s", m_szOutFile, strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}

// iterate over chromosomes
for(ChromIdx = 1; ChromIdx <= pReadsetMetadata->NumChroms; ChromIdx++)
	{
	pChromMetadata = LocateChromMetadataFor(ReadSetID, ChromIdx);
	ChromID = pChromMetadata->ChromID;
	if((pPBAs = LoadSampleChromPBAs(ReadSetID, ChromID)) == nullptr)
		break;
	if (ChromIdx > 1 && LineLen > 0)
		m_pOutBuffer[m_OutBuffIdx++] = '\n';
	m_OutBuffIdx += sprintf((char *)&m_pOutBuffer[m_OutBuffIdx], ">%s\n", LocateChrom(ChromID));
	LineLen = 0;
	for (BaseIdx = 0; BaseIdx < pChromMetadata->ChromLen; BaseIdx++, pPBAs++)
		{
		Alleles = *pPBAs;
		FastaBase = PBAFastaBase(Alleles);
		m_pOutBuffer[m_OutBuffIdx++] = FastaBase;
		if (++LineLen == 70)
			{
			m_pOutBuffer[m_OutBuffIdx++] = '\n';
			LineLen = 0;
			}
		if ((m_OutBuffIdx + 1000) > m_AllocOutBuff)
			{
			if (!CUtility::RetryWrites(m_hOutFile, m_pOutBuffer, m_OutBuffIdx))
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "PBA2Fasta: Fatal error in RetryWrites()");
				Reset();
				return(eBSFerrFileAccess);
				}
			m_OutBuffIdx = 0;
			}
		}
	DeleteSampleChromPBAs(ReadSetID, ChromID);
	}

if (m_OutBuffIdx)
	{
	if (!CUtility::RetryWrites(m_hOutFile, m_pOutBuffer, m_OutBuffIdx))
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "PBA2Fasta: Fatal error in RetryWrites()");
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
gDiagnostics.DiagOut(eDLInfo, gszProcName, "PBA2Fasta: Total of %zd bases written to Fasta file - %s", m_NumPBAbases, m_szOutFile);
gDiagnostics.DiagOut(eDLInfo, gszProcName, "PBA2Fasta: Total of %zd monoallelic bases", m_NumMonoAllelic);
gDiagnostics.DiagOut(eDLInfo, gszProcName, "PBA2Fasta: Total of %zd diallelic bases", m_NumDiAllelic);
gDiagnostics.DiagOut(eDLInfo, gszProcName, "PBA2Fasta: Total of %zd nonallelic bases", m_NumNonAllelic);
return(0);
}

int
CPBAutils::CreateMutexes(void)
{
if (m_bMutexesCreated)
	return(eBSFSuccess);

#ifdef _WIN32
if ((m_hSerialiseAccess = CreateMutex(nullptr, false, nullptr)) == nullptr)
	{
#else
if (pthread_mutex_init(&m_hSerialiseAccess, nullptr) != 0)
	{
#endif
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Fatal: unable to create mutex");
	return(eBSFerrInternal);

	}
m_FastSerialise = 0;
m_bMutexesCreated = true;
return(eBSFSuccess);
}

void
CPBAutils::DeleteMutexes(void)
{
	if (!m_bMutexesCreated)
		return;
#ifdef _WIN32
	CloseHandle(m_hSerialiseAccess);
#else
	pthread_mutex_destroy(&m_hSerialiseAccess);
#endif
	m_bMutexesCreated = false;
}

void
CPBAutils::AcquireSerialise(void)
{
	uint32_t WaitRslt;
#ifdef _WIN32
	WaitRslt = WaitForSingleObject(m_hSerialiseAccess, INFINITE);
	if (WaitRslt != WAIT_OBJECT_0)
	{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Fatal:WaitForSingleObject() returned error %u", WaitRslt);
		Reset();
		exit(1);
	}
#else
	pthread_mutex_lock(&m_hSerialiseAccess);
#endif
}

void
CPBAutils::ReleaseSerialise(void)
{
#ifdef _WIN32
	ReleaseMutex(m_hSerialiseAccess);
#else
	pthread_mutex_unlock(&m_hSerialiseAccess);
#endif
}

void
CPBAutils::AcquireFastSerialise(void)
{
	int SpinCnt = 500;
	int BackoffMS = 5;

#ifdef _WIN32
	while (InterlockedCompareExchange(&m_FastSerialise, 1, 0) != 0)
	{
		if (SpinCnt -= 1)
			continue;
		CUtility::SleepMillisecs(BackoffMS);
		SpinCnt = 100;
		if (BackoffMS < 500)
			BackoffMS += 2;
	}
#else
	while (__sync_val_compare_and_swap(&m_FastSerialise, 0, 1) != 0)
	{
		if (SpinCnt -= 1)
			continue;
		CUtility::SleepMillisecs(BackoffMS);
		SpinCnt = 100;
		if (BackoffMS < 500)
			BackoffMS += 2;
	}
#endif
}

void
CPBAutils::ReleaseFastSerialise(void)
{
#ifdef _WIN32
	InterlockedCompareExchange(&m_FastSerialise, 0, 1);
#else
	__sync_val_compare_and_swap(&m_FastSerialise, 1, 0);
#endif
}

// NOTE: Readsets are checked for uniqueness as readsets must be unique within a given readetset type
uint32_t		// returned readset identifier, 0 if unable to accept this readset name
CPBAutils::AddReadset(char* pszReadset, // associate unique identifier with this readset name
	uint8_t ReadsetType)	// 0: founder, 1: progeny, 2: control
{
	uint32_t ReadsetNameIdx;
	int ReadsetNameLen;
	char Type;
	char* pszLAname;
	Type = '0' + (char)ReadsetType;

	// with any luck the sequence name will be same as the last accessed
	if ((pszLAname = LocateReadset(m_LAReadsetNameID)) != nullptr)
	{
		if (pszLAname[-1] == Type &&	// prefixed with ReadsetType
			!stricmp(pszReadset, pszLAname))
			return(0); // non-unique within ReadsetType!
	}

	// iterate over all known readsets in case this readset to add is a duplicate
	for (ReadsetNameIdx = 0; ReadsetNameIdx < m_NumReadsetNames; ReadsetNameIdx++)
	{
		pszLAname = &m_szReadsetNames[m_szReadsetIdx[ReadsetNameIdx]];
		if (*pszLAname++ != Type)
			continue;
		if (!stricmp(pszReadset, pszLAname))
		{
			m_LAReadsetNameID = ReadsetNameIdx + 1;
			return(0); // non-unique within ReadsetType!
		}
	}

	// not a duplicate
	ReadsetNameLen = (int)strlen(pszReadset);
	if ((m_NxtszReadsetIdx + ReadsetNameLen + 2) > (int)sizeof(m_szReadsetNames))
		return(0);		// no more space to hold reeadset name, treating as if a non-unique!
	if (m_NumReadsetNames == cMaxPBAReadsets)
		return(0);		// unable to hold any more readsets, treating as if a non-unique!

	m_szReadsetIdx[m_NumReadsetNames++] = m_NxtszReadsetIdx;
	pszLAname = &m_szReadsetNames[m_NxtszReadsetIdx];
	*pszLAname++ = Type;
	strcpy(pszLAname, pszReadset);
	m_NxtszReadsetIdx += ReadsetNameLen + 2;
	m_LAReadsetNameID = m_NumReadsetNames;
	return(m_LAReadsetNameID);
}

uint32_t		// returned Readset identifier, 0 if unable to locate this Readset name
CPBAutils::LocateReadset(char* pszReadset, // return unique identifier associated with this Readset name
	uint8_t ReadsetType)	// 0: founder, 1: progeny, 2: control
{
	uint32_t ReadsetNameIdx;
	char Type;
	char* pszLAReadset;

	Type = '0' + (char)ReadsetType;
	// with any luck the Readset name will be same as the last accessed
	if ((pszLAReadset = LocateReadset(m_LAReadsetNameID)) != nullptr)
	{
		if (pszLAReadset[-1] == Type &&	// prefixed with ReadsetType
			!stricmp(pszReadset, pszLAReadset))
			return(m_LAReadsetNameID);
	}

	// iterate over all known Readsets
	for (ReadsetNameIdx = 0; ReadsetNameIdx < m_NumReadsetNames; ReadsetNameIdx++)
	{
		pszLAReadset = &m_szReadsetNames[m_szReadsetIdx[ReadsetNameIdx]];
		if (*pszLAReadset++ != Type)		// readset must be that of ReadsetType
			continue;
		if (!stricmp(pszReadset, pszLAReadset))
		{
			m_LAReadsetNameID = ReadsetNameIdx + 1;
			return(m_LAReadsetNameID);
		}
	}
	return(0);
}

char*
CPBAutils::LocateReadset(uint32_t ReadsetID)
{
	uint32_t Idx;
	Idx = ReadsetID & 0x0fffffff;			// mask out any potential ReadsetType
	if (Idx < 1 || Idx > m_NumReadsetNames)
		return(nullptr);
	return(&(m_szReadsetNames[m_szReadsetIdx[Idx - 1] + 1])); // skipping lead char which is the ReadsetType
}



uint32_t		// returned chrom identifier, 0 if unable to accept this chromosome name
CPBAutils::AddChrom(char* pszChrom) // associate unique identifier with this chromosome name
{
int32_t ChromNameIdx;
int ChromNameLen;
char* pszLAname;

// with any luck the sequence name will be same as the last accessed
if ((pszLAname = LocateChrom(m_LAChromNameID)) != nullptr)
	if (!stricmp(pszChrom, pszLAname))
		return(m_LAChromNameID);

// iterate over all known chroms in case this chrom to add is a duplicate
for (ChromNameIdx = 0; ChromNameIdx < m_NumChromNames; ChromNameIdx++)
	if (!stricmp(pszChrom, &m_szChromNames[m_szChromIdx[ChromNameIdx]]))
	{
	m_LAChromNameID = ChromNameIdx + 1;
	return(m_LAChromNameID);
	}

// chrom is not a duplicate
ChromNameLen = (int)strlen(pszChrom);
if ((m_NxtszChromIdx + ChromNameLen + 1) > (int)sizeof(m_szChromNames))
	return(0);
if (m_NumChromNames == cMaxChromNames)
	return(0);

m_szChromIdx[m_NumChromNames++] = m_NxtszChromIdx;
strcpy(&m_szChromNames[m_NxtszChromIdx], pszChrom);
m_NxtszChromIdx += ChromNameLen + 1;
m_LAChromNameID = m_NumChromNames;
return(m_LAChromNameID);
}


uint32_t		// returned chrom identifier, 0 if unable to locate this chromosome name
CPBAutils::LocateChrom(char* pszChrom) // return unique identifier associated with this chromosome name
{
	int32_t ChromNameIdx;
	char* pszLAname;

	// with any luck the sequence name will be same as the last accessed
	if ((pszLAname = LocateChrom(m_LAChromNameID)) != nullptr)
		if (!stricmp(pszChrom, pszLAname))
			return(m_LAChromNameID);

	// iterate over all known chroms
	for (ChromNameIdx = 0; ChromNameIdx < m_NumChromNames; ChromNameIdx++)
		if (!stricmp(pszChrom, &m_szChromNames[m_szChromIdx[ChromNameIdx]]))
		{
			m_LAChromNameID = ChromNameIdx + 1;
			return(m_LAChromNameID);
		}
	return(0);
}

char*
CPBAutils::LocateChrom(uint32_t ChromID)
{
	if (ChromID < 1 || (int32_t)ChromID > m_NumChromNames)
		return(nullptr);
	return(&m_szChromNames[m_szChromIdx[ChromID - 1]]);
}


bool					// true if chrom is accepted, false if chrom not accepted
CPBAutils::AcceptThisChromID(uint32_t ChromID)
{
char* pzChrom;
if ((pzChrom = LocateChrom(ChromID)) == nullptr)
	return(false);
return(AcceptThisChromName(pzChrom));
}


bool					// true if chrom is accepted, false if chrom not accepted
CPBAutils::AcceptThisChromName(char* pszChrom,   // chromosome name
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


int
CPBAutils::AllocChromMetadata(void)
{
	uint32_t ToAllocdChromMetadata;
	tsPUChromMetadata* pChromMetadata;
	size_t memreq;
	if (m_pChromMetadata == nullptr)					// may be nullptr first time in
	{
		memreq = cAllocChromMetadata * sizeof(tsPUChromMetadata);
#ifdef _WIN32
		m_pChromMetadata = (tsPUChromMetadata*)malloc((size_t)memreq);
		if (m_pChromMetadata == nullptr)
		{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "AllocChromMetadata: Memory allocation of %zd bytes failed", (int64_t)memreq);
			return(eBSFerrMem);
		}
#else
		m_pChromMetadata = (tsPUChromMetadata*)mmap(nullptr, (size_t)memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
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
			size_t memreq = ToAllocdChromMetadata * sizeof(tsPUChromMetadata);
#ifdef _WIN32
			pChromMetadata = (tsPUChromMetadata*)realloc(m_pChromMetadata, memreq);
			if (pChromMetadata == nullptr)
			{
#else
			pChromMetadata = (tsPUChromMetadata*)mremap(m_pChromMetadata, m_AllocdChromMetadataMem, memreq, MREMAP_MAYMOVE);
			if (pChromMetadata == MAP_FAILED)
			{
#endif
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "AllocChromMetadata: Memory reallocation to %zd bytes failed - %s", (int64_t)memreq, strerror(errno));
				return(eBSFerrMem);
			}
			m_pChromMetadata = pChromMetadata;
			m_AllocdChromMetadataMem = memreq;
			m_AllocdChromMetadata = ToAllocdChromMetadata;
			}
	return(eBSFSuccess);
		}

uint8_t*
CPBAutils::AllocPBAs(uint32_t ChromLen)	// allocate memory to hold at least this many packed base alleles
{
	uint8_t* pPBAs;
	size_t memreq;
	memreq = (size_t)ChromLen;	// no safety margin!
#ifdef _WIN32
	pPBAs = (uint8_t*)malloc((size_t)memreq);
	if (pPBAs == nullptr)
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AllocPBAs Memory allocation of %zd bytes failed", (int64_t)memreq);
#else
	pPBAs = (uint8_t*)mmap(nullptr, (size_t)memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	if (pPBAs == MAP_FAILED)
	{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AllocPBAs: Memory allocation of %zd bytes through mmap()  failed", (int64_t)memreq, strerror(errno));
		pPBAs = nullptr;
	}
#endif
	return(pPBAs);
}


uint8_t*								// returned pointer to start of PBA
CPBAutils::LocatePBAfor(uint32_t ReadSetID,		// readset identifier 
	uint32_t ChromID)			// chrom identifier
{
	tsPUReadsetMetadata* pReadsetMetadata;
	tsPUChromMetadata* pChromMetadata;
	uint32_t CurChromMetadataIdx;

	if (ReadSetID > m_NumReadsetNames || ReadSetID == 0)
		return(nullptr);
	pReadsetMetadata = &m_Readsets[ReadSetID - 1];
	CurChromMetadataIdx = pReadsetMetadata->StartChromMetadataIdx;
	for (uint32_t ChromIdx = 0; ChromIdx < pReadsetMetadata->NumChroms && CurChromMetadataIdx != 0; ChromIdx++)
	{
		pChromMetadata = &m_pChromMetadata[CurChromMetadataIdx - 1];
		if (pChromMetadata->ChromID == ChromID)
			return(pChromMetadata->pPBAs);
		CurChromMetadataIdx = pChromMetadata->NxtChromMetadataIdx;
	}
	return(nullptr);
}

tsPUChromMetadata*								// returned pointer to chromosome metadata
CPBAutils::LocateChromMetadataFor(uint32_t ReadSetID,		// readset identifier 
	uint32_t ChromID)			// chrom identifier
{
	tsPUReadsetMetadata* pReadsetMetadata;
	tsPUChromMetadata* pChromMetadata;
	uint32_t CurChromMetadataIdx;

	if (ReadSetID > m_NumReadsetNames || ReadSetID == 0)
		return(nullptr);
	pReadsetMetadata = &m_Readsets[ReadSetID - 1];
	CurChromMetadataIdx = pReadsetMetadata->StartChromMetadataIdx;
	for (uint32_t ChromIdx = 0; ChromIdx < pReadsetMetadata->NumChroms && CurChromMetadataIdx != 0; ChromIdx++)
	{
		pChromMetadata = &m_pChromMetadata[CurChromMetadataIdx - 1];
		if (pChromMetadata->ChromID == ChromID)
			return(pChromMetadata);
		CurChromMetadataIdx = pChromMetadata->NxtChromMetadataIdx;
	}
	return(nullptr);
}



uint32_t				// returns number of unprocessed bytes in buffer
CPBAutils::FillInBuffer(uint32_t MinRequired, uint32_t MaxRequired) // try and fill input buffer with at least MinRequired or refill (if < MinRequired) up to MaxRequired (if 0 then max buffer allocation)
{
	int NumRead;
	uint32_t UnProcessed;
	if (MaxRequired == 0)
		MaxRequired = m_AllocInBuff;
	else
		if (MaxRequired < MinRequired)
			MaxRequired = MinRequired;
	UnProcessed = m_InNumBuffered - m_InNumProcessed;

	// if already filled to MinRequired then no further action required
	if (MinRequired <= UnProcessed)
		return(UnProcessed);

	// copy down any bytes which are yet to be processed
	if (UnProcessed) {
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
	do {
		if ((NumRead = read(m_hInFile, &m_pInBuffer[m_InNumProcessed], MaxRequired - m_InNumBuffered)) < 0)
		{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error %s attempting to read from file", strerror(errno));
			return(0);
		}
		m_InNumBuffered += (uint32_t)NumRead;
	} while (NumRead > 0 && m_InNumBuffered < MinRequired);
	m_InFileOfs = _lseeki64(m_hInFile, 0, SEEK_CUR);             // record file offset at which next file read will start from
	return(m_InNumBuffered);
}


// Fasta2PBA
// Convert Fasta into PBA
int
CPBAutils::Fasta2PBA(char *pszExperimentID, // experiment identifier
	char *pszReferenceID,	  // reference identifier
	char *pszReadsetID,		  // readset identifier
	char* pszInFile,		  // Fasta format input file
	char *pszOutFile)		  // PBA format output file
{
	CFasta Fasta;
	uint8_t* pPackedBaseAlleles;
	size_t AvailBuffSize;
	char szChromName[cBSFSourceSize];
	char szDescription[cBSFDescriptionSize];
	uint32_t SeqLen;
	int Descrlen;
	bool bChromSeq;
	int Rslt;
	int ChromID;
	uint32_t PBAIdx;
	uint32_t ChromPBAlen;
	uint32_t ChromLenIdx;
	uint32_t ChromIndeterminates;
	uint32_t NumPBAChroms;
	size_t AssemblyChromIndeterminates;
	size_t AssembSeqLen;

	if ((Rslt = Fasta.Open(pszInFile, true)) != eBSFSuccess)
		{
		if (Rslt != eBSFerrNotFasta)
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Fasta2PBA: Unable to open '%s' [%s] %s", pszInFile, Fasta.ErrText((teBSFrsltCodes)Rslt), Fasta.GetErrMsg());
		return(Rslt);
		}

	m_AllocInBuff = (size_t)cMaxAXAllocBuffChunk * 16;	// allocating to hold each fasta input sequence, will be reallocated as required to hold complete sequences
	// note malloc is used as can then simply realloc to expand as may later be required
	if ((m_pInBuffer = (uint8_t*)malloc(m_AllocInBuff)) == nullptr)
	{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Fasta2PBA:- Unable to allocate memory (%zd bytes) for sequence buffer", m_AllocInBuff);
		Fasta.Close();
		return(eBSFerrMem);
	}
	m_InNumBuffered = 0;

	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Fasta2PBA: Loading Fasta sequences from %s..", pszInFile);


#ifdef _WIN32
	m_hOutFile = open(pszOutFile, (O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC), (_S_IREAD | _S_IWRITE));
#else
	if ((m_hOutFile = open64(pszOutFile, O_WRONLY | O_CREAT, S_IREAD | S_IWRITE)) != -1)
		if (ftruncate(m_hOutFile, 0) != 0)
		{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to create/truncate %s - %s", pszOutFile, strerror(errno));
			Reset();
			return(eBSFerrCreateFile);
		}
#endif
	if (m_hOutFile < 0)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "ReportHaplotypesByProgeny: Unable to create/truncate %s - %s", pszOutFile, strerror(errno));
		Reset();
		return(eBSFerrCreateFile);
		}

	// Write out the PBA header
	// This header contains a series of '\n' separated tagname:values
	// Following the header are a variable number of chromosomes and associated sequences
m_InNumBuffered = sprintf((char*)m_pInBuffer, "Type:%s\nVersion:1\nExperimentID:%s\nReferenceID:%s\nReadsetID:%s", "PbA", pszExperimentID, pszReferenceID, pszReadsetID);
m_InNumBuffered += 1;
if (!CUtility::RetryWrites(m_hOutFile, m_pInBuffer, m_InNumBuffered))
	{
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "RetryWrites() error");
	return(eBSFerrWrite);
	}
m_InNumBuffered = 0;
bChromSeq = false;
ChromPBAlen = 0;
ChromIndeterminates = 0;
AssemblyChromIndeterminates = 0;
AssembSeqLen = 0;
NumPBAChroms = 0;
while ((Rslt = SeqLen = Fasta.ReadSequence(&m_pInBuffer[m_InNumBuffered], (int)(m_AllocInBuff - m_InNumBuffered), true, false)) > eBSFSuccess)
	{
	if (SeqLen == eBSFFastaDescr)		// just read a descriptor line
		{
		if (bChromSeq)				// not the first descriptor so write out previous with it's associated sequence
			{
			// check if chromosome is to retained
			if(AcceptThisChromName(szChromName))
				{
				// before accepting chrom then ensure that it's sequence length matches the BED chromosome sizes
				ChromID = LocateChrom(szChromName);
				if (ChromPBAlen != m_ChromSizes[ChromID - 1])
					{
					gDiagnostics.DiagOut(eDLFatal, gszProcName, "Input file '%s' has chromosome '%s' size mismatch - expected size was %u, actual size %u", pszInFile, szChromName, m_ChromSizes[ChromID - 1], ChromPBAlen);
					return(eBSFerrChrom);
					}

				*(uint32_t*)&m_pInBuffer[ChromLenIdx] = ChromPBAlen;
				if (!CUtility::RetryWrites(m_hOutFile, m_pInBuffer, m_InNumBuffered))
					{
					gDiagnostics.DiagOut(eDLInfo, gszProcName, "RetryWrites() error");
					return(eBSFerrWrite);
					}
				gDiagnostics.DiagOut(eDLInfo, gszProcName, "Chrom '%s' is length %d and has %d indeterminate bases", szChromName, ChromPBAlen, ChromIndeterminates);
				AssembSeqLen += (size_t)ChromPBAlen;
				AssemblyChromIndeterminates += (size_t)ChromIndeterminates;
				NumPBAChroms++;
				}
			}
		ChromIndeterminates = 0;
		m_InNumBuffered = 0;
		Descrlen = Fasta.ReadDescriptor(szDescription, cBSFDescriptionSize);
		sscanf(szDescription, " %s[ ,]", szChromName);
		pPackedBaseAlleles = &m_pInBuffer[m_InNumBuffered];
		*pPackedBaseAlleles = (uint8_t)strlen(szChromName);
		strcpy((char*)&pPackedBaseAlleles[1], szChromName);
		pPackedBaseAlleles += 2 + *pPackedBaseAlleles;
		*(uint32_t*)pPackedBaseAlleles = 0;
		ChromLenIdx = (uint32_t)(pPackedBaseAlleles - m_pInBuffer);
		pPackedBaseAlleles += 4;
		m_InNumBuffered = (uint32_t)(pPackedBaseAlleles - m_pInBuffer);
		ChromPBAlen = 0;
		continue;		
		}

			// Fasta to PBA alleles
	pPackedBaseAlleles = &m_pInBuffer[m_InNumBuffered];
	for (PBAIdx = 0; PBAIdx < SeqLen; PBAIdx++,pPackedBaseAlleles++)
		{
		switch (*pPackedBaseAlleles & ~cRptMskFlg) {
			case eBaseA:
				*pPackedBaseAlleles = 0xc0;
				break;
			case eBaseC:
				*pPackedBaseAlleles = 0x30;
				break;
			case eBaseG:
				*pPackedBaseAlleles = 0x0c;
				break;
			case eBaseT:
				*pPackedBaseAlleles = 0x03;
				break;
			default:
				ChromIndeterminates++;
				*pPackedBaseAlleles = 0x00;
				break;
			}
		}

		// increase buffering as may be required to handle ReadSequence() on next iteration still extending current sequence rather than a new fasta descriptor
	bChromSeq = true;
	m_InNumBuffered += SeqLen;
	ChromPBAlen += SeqLen;
	AvailBuffSize = m_AllocInBuff - m_InNumBuffered;
	if (AvailBuffSize < (size_t)(cMaxAXAllocBuffChunk / 8))
		{
		size_t NewSize = (size_t)cMaxAXAllocBuffChunk + m_AllocInBuff;
		uint8_t* pTmp;
		if ((pTmp = (uint8_t*)realloc(m_pInBuffer, NewSize)) == nullptr)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessFastaFile:- Unable to reallocate memory (%u bytes) for sequence buffer", (uint32_t)NewSize);
			return(eBSFerrMem);
			}
		m_pInBuffer = pTmp;
		m_AllocInBuff = (uint32_t)NewSize;
		}
	}

if (Rslt == eBSFSuccess && bChromSeq && m_InNumBuffered)	
	{
	// check if chromosome is to retained
	if(AcceptThisChromName(szChromName))
		{
		// before accepting chrom then ensure that it's sequence length matches the BED chromosome sizes
		ChromID = LocateChrom(szChromName);
		if (ChromPBAlen != m_ChromSizes[ChromID - 1])
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Input file '%s' has chromosome '%s' size mismatch - expected size was %u, actual size %u", pszInFile, szChromName, m_ChromSizes[ChromID - 1], ChromPBAlen);
			return(eBSFerrChrom);
			}

		*(uint32_t*)&m_pInBuffer[ChromLenIdx] = ChromPBAlen;
		if (!CUtility::RetryWrites(m_hOutFile, m_pInBuffer, m_InNumBuffered))
			{
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "RetryWrites() error");
			return(eBSFerrWrite);
			}
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Chrom '%s' is length %d and has %d indeterminate bases", szChromName, ChromPBAlen, ChromIndeterminates);
		AssembSeqLen += (size_t)ChromPBAlen;
		AssemblyChromIndeterminates += (size_t)ChromIndeterminates;
		NumPBAChroms++;
		}
	}
if(Rslt == eBSFSuccess)
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Generated PBA contains %d sequences totaling %dbp with %d indeterminate bases",NumPBAChroms, AssembSeqLen, AssemblyChromIndeterminates);
else
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Fasta2PBA: Processing errors..");
	return(Rslt);
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
gDiagnostics.DiagOut(eDLInfo, gszProcName, "Fasta2PBA: Completed writing to PBA file %s..", pszOutFile);
return(Rslt);
}



int32_t					// returned readset identifier (1..n) or < 0 if errors
CPBAutils::LoadPBAFile(char* pszFile,	// load chromosome metadata and PBA data from this file, filters out chroms - AcceptThisChromName()
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
	tsPUReadsetMetadata* pReadsetMetadata;
	tsPUChromMetadata* pChromMetadata;
	tsPUChromMetadata* pPrevChromMetadata;
	int64_t FileOfsPBA;
	uint8_t* pBuff;
	int scanlen;
	int NumTags;

	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Loading '%s' file", pszFile);
#ifdef _WIN32
	m_hInFile = open(pszFile, O_READSEQ);		// file access is normally sequential..
#else
	m_hInFile = open64(pszFile, O_READSEQ);		// file access is normally sequential..
#endif
	if (m_hInFile == -1)							// check if file open succeeded
	{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to open input file '%s' : %s", pszFile, strerror(errno));
		return(eBSFerrOpnFile);
	}

	// attempt to load the readset metadata 
	m_InNumBuffered = 0;
	m_InNumProcessed = 0;
	if ((FillInBuffer(500u, min(m_AllocInBuff, 500u)) == 0) || m_InNumBuffered < 9) // 500 will cover the maximally sized header
	{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to read at least a partial header from input file '%s'", pszFile);
		return(eBSFerrOpnFile);
	}


	// check file type is PBA
	if (strncmp((char*)m_pInBuffer, "Type:PbA\n", 9))
	{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Input file '%s' exists, unable to parse file type header tag, is it a packed base allele file", pszFile);
		return(eBSFerrOpnFile);
	}
	m_InNumProcessed = 9;

	if (m_InNumBuffered < 500)
	{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Input file '%s' exists, unable to read complete header and initial chromosome metadata, expecting file size >= 500, is it a packed base allele file", pszFile);
		return(eBSFerrOpnFile);
	}

	// parse out tagnames and associated values
	pBuff = &m_pInBuffer[m_InNumProcessed];
	NumTags = sscanf((char*)pBuff, "Version:%d\nExperimentID:%99[^\n]\nReferenceID:%99[^\n]\nReadsetID:%99[^\n]%n", &Version, szExperimentID, szRefAssemblyID, szReadsetID, &scanlen);
	if (NumTags != 4 || Version != 1 || szExperimentID[0] == '\0' || szRefAssemblyID[0] == '\0' || szReadsetID[0] == '\0')
	{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Input file '%s' exists as a packed base allele file but inconsistencies in header tag values", pszFile);
		return(eBSFerrOpnFile);
	}
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Metadata ExperimentID:%s, ReferenceID:%s, ReadsetID:%s", szExperimentID, szRefAssemblyID, szReadsetID);

	m_InNumProcessed += scanlen + 1;		// header tags were '\n' terminated except for final which was '\0' terminated

	m_InFileOfs = m_InNumProcessed;         // 1st chromosome starts immediately following the header metadata
	if (_lseeki64(m_hInFile, m_InFileOfs, SEEK_SET) != m_InFileOfs)
	{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Input file '%s' exists as a packed base allele file but unable to seek past header metadata", pszFile);
		return(eBSFerrOpnFile);
	}

	if ((ReadsetID = AddReadset(szReadsetID, ReadsetType)) == 0)
	{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Input file '%s' duplicates the ReadsetID '%s' of a previously loaded readset", pszFile, szReadsetID);
		return(eBSFerrOpnFile);
	}

	pReadsetMetadata = &m_Readsets[ReadsetID - 1];
	memset(pReadsetMetadata, 0, sizeof(*pReadsetMetadata));
	pReadsetMetadata->ReadsetType = ReadsetType;
	pReadsetMetadata->NumChroms = 0;
	strcpy(pReadsetMetadata->szFileName, pszFile);
	strcpy(pReadsetMetadata->szExperimentID, szExperimentID);
	strcpy(pReadsetMetadata->szRefAssemblyID, szRefAssemblyID);
	pReadsetMetadata->ReadsetID = ReadsetID;
	pReadsetMetadata->StartChromID = 0;
	pReadsetMetadata->StartChromMetadataIdx = 0;
	pReadsetMetadata->NxtFileChromOfs = m_InFileOfs;

	int ChromNameLen;
	char* pszChromName;
	uint32_t ChromLen;

	// iterate over all chromosomes
	m_InNumProcessed = 0;
	m_InNumBuffered = 0;
	PrevChromMetadataIdx = 0;
	while (FillInBuffer((uint32_t)110, 110) == 110) // reading chromosome metadata
	{
		pBuff = &m_pInBuffer[m_InNumProcessed++];
		ChromNameLen = (int)*pBuff++;
		pszChromName = (char*)pBuff;
		pBuff += 1 + ChromNameLen;
		ChromLen = *(int32_t*)pBuff;
		FileOfsPBA = pReadsetMetadata->NxtFileChromOfs + ChromNameLen + 6;
		pReadsetMetadata->NxtFileChromOfs += ChromNameLen + 6 + ChromLen;
		
		// check if this chromosome is to be retained for further processing
		if (!AcceptThisChromName(pszChromName))
		{
			// not accepting this chromosome
			m_InFileOfs = pReadsetMetadata->NxtFileChromOfs;         // skip over curent chrom's PBAs to start of next chromosome
			if ((pReadsetMetadata->NxtFileChromOfs = _lseeki64(m_hInFile, m_InFileOfs, SEEK_SET)) != m_InFileOfs)
				break;
			m_InNumProcessed = 0;
			m_InNumBuffered = 0;
			continue;
		}

		// before accepting chrom then ensure that it's PBA length matches the BED chromosome sizes
		ChromID = AddChrom(pszChromName);
		if (ChromLen != m_ChromSizes[ChromID - 1])
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Input file '%s' has chromosome '%s' size mismatch - expected size was %u, actual size %u", pszFile, pszChromName, m_ChromSizes[ChromID - 1], ChromLen);
			return(eBSFerrChrom);
			}

		if ((Rslt = AllocChromMetadata()) != eBSFSuccess)
			return(Rslt);

		pChromMetadata = &m_pChromMetadata[m_UsedNumChromMetadata++];
		if (PrevChromMetadataIdx != 0)
			{
			pPrevChromMetadata = &m_pChromMetadata[PrevChromMetadataIdx - 1];
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
		if (bChromMetaOnly) // loading chrom metadata only without actually loading the chromosome PBAs?
			{
			m_InFileOfs = pReadsetMetadata->NxtFileChromOfs;         // skip over curent chrom's PBAs to start of next chromosome
			if ((pReadsetMetadata->NxtFileChromOfs = _lseeki64(m_hInFile, m_InFileOfs, SEEK_SET)) != m_InFileOfs)
				break;
			m_InNumProcessed = 0;
			m_InNumBuffered = 0;
			continue;
			}
		// loading PBAs
		pChromMetadata->pPBAs = AllocPBAs(ChromLen);
		m_InNumProcessed += ChromNameLen + 5;
		FillInBuffer((uint32_t)ChromLen, ChromLen);
		pBuff = &m_pInBuffer[m_InNumProcessed];

		//	validate PBA allele composition, earlier releases were retaining very low allele proportions so validating and in-place replacing these as being non-alignments
		int NumErrs = ValidatePBAs(ChromLen, pBuff, true, true);
		//
		memcpy(pChromMetadata->pPBAs, pBuff, ChromLen);
		if(m_PBAsTrim5 > 0 || m_PBAsTrim3 > 0)		// trimming aligned segment boundaries as these may have higher polyallele rates
			TrimPBAs(m_PBAsTrim5, m_PBAsTrim3, ChromLen, pChromMetadata->pPBAs);
		m_InNumProcessed += ChromLen;
		}

	close(m_hInFile);
	m_hInFile = -1;
	return((int32_t)ReadsetID);
}

int			// returned number of PBAs which are non-conformant
CPBAutils::ValidatePBAs(uint32_t Length,
	uint8_t* pPBAs,
	bool bSetNoAlleles, // if non-conformant then overwrite *pAlleles to be no alleles present
	bool bNormalise)    // normalise for coverage, low coverage major alleles (2,0,0,0) are normalised to high coverage (3,0,0,0)
{
	uint32_t Idx;
	int32_t NumErrors = 0;
	for (Idx = 0; Idx < Length; Idx++, pPBAs++)
	{
		if (!ValidatePBA(pPBAs, bSetNoAlleles, bNormalise))
			NumErrors++;
	}
	return(NumErrors);
}

// trim 5' and 3' aligned segments within the PBAs attempting to reduce sequencing error induced false alleles
int			// returns number of non-trimmed loci in the pPBAs
CPBAutils::TrimPBAs(uint32_t Trim5,	// trim 5' this many aligned PBA bases from each aligned segment
	uint32_t Trim3,	// trim 3' this many aligned PBA bases from each aligned segment
	uint32_t PBALen,	// pPBAs contains this many packed base alleles
	uint8_t* pPBAs)		// packed base alleles to be processed
{
	uint32_t NonTrimmed;
	uint32_t ToTrim;
	uint32_t Loci;
	uint8_t* pPBA;
	if (pPBAs == nullptr || PBALen == 0)
		return(0);
	if (Trim5 == 0 && Trim3 == 0)
		return(PBALen);
	NonTrimmed = PBALen;

	// first the 5' trim
	if (Trim5)
	{
		ToTrim = Trim5;
		pPBA = pPBAs;
		for (Loci = 0; Loci < PBALen; Loci++, pPBA++)
		{
			if (*pPBA == 0)
			{
				ToTrim = Trim5;
				continue;
			}
			if (ToTrim)
			{
				ToTrim -= 1;
				NonTrimmed -= 1;
				*pPBA = 0;
			}
		}
	}
	// then the 3' trim
	if (Trim3)
	{
		ToTrim = Trim3;
		pPBA = &pPBAs[PBALen - 1];
		for (Loci = 0; Loci < PBALen; Loci++, pPBA--)
		{
			if (*pPBA == 0)
			{
				ToTrim = Trim3;
				continue;
			}
			if (ToTrim)
			{
				NonTrimmed -= 1;
				ToTrim -= 1;
				*pPBA = 0;
			}
		}
	}
	return(NonTrimmed);
}

void
CPBAutils::DeleteAllChromPBAs(void) // delete all currently loaded PBAs - all sample chroms
{
if (m_pChromMetadata != nullptr)
	{
	tsPUChromMetadata* pChromMetadata = m_pChromMetadata;
	for (uint32_t Idx = 0; Idx < m_UsedNumChromMetadata; Idx++, pChromMetadata++)
		{
		if (pChromMetadata->pPBAs != nullptr)
#ifdef _WIN32
			free(pChromMetadata->pPBAs);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
		if (pChromMetadata->pPBAs != MAP_FAILED)
				munmap(pChromMetadata->pPBAs, pChromMetadata->ChromLen);
#endif
		pChromMetadata->pPBAs = nullptr; 
		}
	}
}

bool  // false: no previous memory allocated for containing the chromosome PBAs, true: allocation was present and has been deleted
CPBAutils::DeleteSampleChromPBAs(uint32_t SampleID,   // Sample identifier
	uint32_t ChromID)    // chrom identifier
{
	tsPUChromMetadata* pChromMetadata;

	// returned pointer to chromosome metadata
	if ((pChromMetadata = LocateChromMetadataFor(SampleID, ChromID)) == nullptr)
		return(false);
	if (pChromMetadata->pPBAs == nullptr)
		return(false);
#ifdef _WIN32
	free(pChromMetadata->pPBAs);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (pChromMetadata->pPBAs != MAP_FAILED)
		munmap(pChromMetadata->pPBAs, pChromMetadata->ChromLen);
#endif
	pChromMetadata->pPBAs = nullptr;
	return(true);
}

uint8_t*  // returned PBA or nullptr if unable to load PBAs for requested sample.chrom 
CPBAutils::LoadSampleChromPBAs(uint32_t SampleID,   // Sample identifier
	uint32_t ChromID)    // chrom identifier specifying which PBAs is to be loaded from SampleID file
{
	int hInFile;
	char* pszChrom;
	int64_t ChromSeekOfs;
	tsPUChromMetadata* pChromMetadata;
	tsPUReadsetMetadata* pReadsetMetadata;
	uint32_t NumRead;
	uint32_t NumLoaded;

	// returned pointer to chromosome metadata
	if ((pChromMetadata = LocateChromMetadataFor(SampleID, ChromID)) == nullptr)
		return(nullptr);

	pszChrom = LocateChrom(ChromID);

	if (pChromMetadata->pPBAs != nullptr) // PBAs may already have been loaded for this chromosome, delete as may be stale
	{
#ifdef _WIN32
		free(pChromMetadata->pPBAs);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
		if (pChromMetadata->pPBAs != MAP_FAILED)
			munmap(pChromMetadata->pPBAs, pChromMetadata->ChromLen);
#endif
		pChromMetadata->pPBAs = nullptr;
	}

	// need to actually load from file
	pReadsetMetadata = &m_Readsets[SampleID - 1];
#ifdef _WIN32
	hInFile = open(pReadsetMetadata->szFileName, O_READSEQ);		// file access is normally sequential..
#else
	hInFile = open64(pReadsetMetadata->szFileName, O_READSEQ);		// file access is normally sequential..
#endif
	if (hInFile == -1)							// check if file open succeeded
	{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadSampleChromPBA: Unable to open input file '%s' : %s", pReadsetMetadata->szFileName, strerror(errno));
		return(nullptr);
	}
	if ((ChromSeekOfs = _lseeki64(hInFile, pChromMetadata->FileOfsPBA, SEEK_SET)) != pChromMetadata->FileOfsPBA)
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
	do {
		if ((NumRead = read(hInFile, &pChromMetadata->pPBAs[NumLoaded], pChromMetadata->ChromLen - NumLoaded)) < 0)
		{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error %s attempting to read chromosome '%s' from file '%s'", strerror(errno), pszChrom, pReadsetMetadata->szFileName);
			close(hInFile);
			hInFile = -1;
			if (pChromMetadata->pPBAs == nullptr)
				return(nullptr);
#ifdef _WIN32
			free(pChromMetadata->pPBAs);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
			if (pChromMetadata->pPBAs != MAP_FAILED)
				munmap(pChromMetadata->pPBAs, pChromMetadata->ChromLen);
#endif
			pChromMetadata->pPBAs = nullptr;
			return(nullptr);
		}
		NumLoaded += NumRead;
	} while (NumRead > 0 && NumLoaded < pChromMetadata->ChromLen);
	close(hInFile);
	hInFile = -1;
	//	validate PBA allele composition, earlier releases were retaining very low allele proportions so validating and in-place replacing these as being non-alignments
	int NumErrs = ValidatePBAs(pChromMetadata->ChromLen, pChromMetadata->pPBAs, true, true);
	//
	if (m_PBAsTrim5 > 0 || m_PBAsTrim3 > 0)		// trimming aligned segment boundaries as these may have higher polyallele rates
		TrimPBAs(m_PBAsTrim5, m_PBAsTrim3, pChromMetadata->ChromLen, pChromMetadata->pPBAs);
	return(pChromMetadata->pPBAs);
}

bool
CPBAutils::ValidatePBA(uint8_t* pAlleles,	// validate that PBA alleles are properly conformant
	bool bSetNoAlleles, // if non-conformant then overwrite *pAlleles to be no alleles present
	bool bNormalise)    // normalise for coverage, low coverage major alleles (2,0,0,0) are normalised to high coverage (3,0,0,0)
{
	uint32_t AlleleIdx;
	uint8_t Allele;
	uint8_t AlleleMsk;
	int32_t Num2s;
	int32_t Num1s;

	if (*pAlleles == 0)
		return(true);

	Num2s = Num1s = 0;
	AlleleMsk = 0x03;
	for (AlleleIdx = 0; AlleleIdx < 4; AlleleIdx++, AlleleMsk <<= 2)
	{
		Allele = (*pAlleles & AlleleMsk) >> (AlleleIdx * 2);
		switch (Allele) {
		case 0x03:
			if ((*pAlleles & ~AlleleMsk) != 0) // can only be one major ...
			{
				if (bSetNoAlleles)
					*pAlleles = 0;
				return(false);
			}
			else
				return(true);
		case 0x02:
			if (bNormalise && (*pAlleles & ~AlleleMsk) == 0) // if low coverage, a major allele is represented as 2,0,0,0 whereas if high coverage then that same allele would have been represented as 3,0,0,0
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
	if (Num2s > 2 || Num1s > 2 || (Num2s + Num1s > 2))
	{
		if (bSetNoAlleles)
			*pAlleles = 0;
		return(false);
	}
	return(true);
}


