/*
This toolkit is a source base clone of 'BioKanga' release 4.4.2 (https://github.com/csiro-crop-informatics/biokanga) and contains
significant source code changes enabling new functionality and resulting process parameterisation changes. These changes have resulted in
incompatibility with 'BioKanga'.

Because of the potential for confusion by users unaware of functionality and process parameterisation changes then the modified source base
and resultant compiled executables have been renamed to 'kit4b' - K-mer Informed Toolkit for Bioinformatics.
The renaming will force users of the 'BioKanga' toolkit to examine scripting which is dependent on existing 'BioKanga'
parameterisations so as to make appropriate changes if wishing to utilise 'kit4b' parameterisations and functionality.

'kit4b' is being released under the Opensource Software License Agreement (GPLv3)
'kit4b' is Copyright (c) 2019, 2020
Please contact Dr Stuart Stephen < stuartjs@g3web.com > if you have any questions regarding 'kit4b'.

Original 'BioKanga' copyright notice has been retained and immediately follows this notice..
*/
/*
 * CSIRO Open Source Software License Agreement (GPLv3)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENSE for the complete license information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */


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

const char* cpszArtef5Seq = "ACACTCTTTCCCTACACGACGCTGTTCCATCT";	// default artifact seq for 5' simulated read ends (Illumina Single End Adapter 1)
const char* cpszArtef3Seq = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"; // default artefact seq for 3' simulated read ends (Illumina Single End Sequencing Primer)


int
Process(etSRPMode PMode,		// processing mode
		etSRSEMode SEMode,	// induced sequencer error rate mode
		bool bPEgen,		// true if paired ends are to be generated
		int PEmin,			// PE minimum insert
		int PEmax,			// PE maximum insert
		double PropRandReads, // generate completely random reads at this rate
		int DistCluster,	// cluster generated reads into windows of this median length, 0 if no clustering
		double SeqErrRate,	// dynamic sequencing errors are to be induced into generated sequences at this rate
		bool bSeqErrProfile,// true if to generate composite sequencer errors with uniform profile (default is Illumina 3' skewed)
		int SNPrate,		// simulate SNPs at this rate per million bases
		int InDelSize,		// simulated InDel size range
		double InDelRate,	// simulated InDel rate per read
		etSRFMode FMode,		// output format
		int NumThreads,		// number of worker threads to use
		char Strand,		// generate for this strand '+' or '-' or for both '*'
		int NumReads,		// number of reads required (will be doubled if paired end reads)
		int ReadLen,		// read lengths
		double Artef5Rate,			// rate (0..1) at which to insert 5' artefact sequences
		int NumArtef5Seqs,			// number of user specified 5' artefact sequences
		char *pszArtef5Seqs[], // 5' artefact sequences
		double Artef3Rate,			// rate (0..1) at which to insert 3' artefact sequences
		int NumArtef3Seqs,			// number of user specified 3' artefact sequences
		char *pszArtef3Seqs[], // 5' artefact sequences
		int CutMin,			// min cut length
		int CutMax,			// max cut length
		bool bDedupe,		// true if unique read sequences only to be generated
		int Region,			// Process regions 0:ALL,1:CDS,2:5'UTR,3:3'UTR,4:Introns,5:5'US,6:3'DS,7:Intergenic (default = ALL)
		int UpDnRegLen,		// if processing regions then up/down stream regulatory length
		char *pszFeatFile,	// optionally generate transcriptome reads from features or genes in this BED file
		char *pszInFile,	// input from this bioseq assembly
		char *pszProfFile,	// input from this profile site preferences file
		char *pszOutPEFile, // output partner paired end simulated reads to this file
		char *pszOutFile,	// output simulated reads to this file
		char *pszOutSNPs);   // output simulated SNP loci to this file


#ifdef _WIN32
int SimReads(int argc, char* argv[])
{
// determine my process name
_splitpath(argv[0],nullptr,nullptr,gszProcName,nullptr);
#else
int SimReads(int argc, char** argv)
{
// determine my process name
CUtility::splitpath((char *)argv[0],nullptr,gszProcName);
#endif
int iScreenLogLevel;		// level of screen diagnostics
int iFileLogLevel;			// level of file diagnostics
char szLogFile[_MAX_PATH];	// write diagnostics to this file
int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure

int Idx;
int LenReq;

etSRPMode PMode;					// processing mode
etSRFMode FMode;					// multifasta
char Strand;				// generate for this strand '+' or '-' or for both '*'
int NumReads;				// number of reads required
int ReadLen;				// read lengths
int CutMin;					// min length
int CutMax;					// max cut length
bool bDedupe;				// true if unique reads only to be generated

int NumberOfProcessors;		// number of installed CPUs
int NumThreads;				// number of threads (0 defaults to number of CPUs)
int SNPrate;				// generate SNPs at this rate per million bases

etSRBEDRegion Region;					// Process regions 0:ALL,1:CDS,2:5'UTR,3:3'UTR,4:Introns,5:5'US,6:3'DS,7:Intergenic (default = ALL)
int UpDnRegLen;				// if processing regions then up/down stream regulatory length

double PropRandReads;		// proportion of reads to be generated with completely random sequences

etSRSEMode SEMode;			// induced sequencer error rate mode
int DistCluster;			// cluster generated reads into windows of this median length, 0 if no clustering
double SeqErrRate;			// dynamic sequencer error rate
int InDelSize;				// simulated InDel size range
double InDelRate;			// simulated InDel rate per read
bool bSeqErrProfile;		// true if to generate composite sequencer errors with uniform profile (default is Illumina 3' skewed)
bool bPEgen;				// true if paired ends are to be generated
int PEmin;					// PE minimum insert
int PEmax;					// PE maximum insert
char szInFile[_MAX_PATH];	// input from this bioseq assembly
char szProfFile[_MAX_PATH];// input from this profile site preferences file
char szOutFile[_MAX_PATH];	// output simulated reads to this file
char szOutPEFile[_MAX_PATH];// output partner paired end simulated reads to this file
char szSNPFile[_MAX_PATH];	// output simulated SNPs to this file
char szFeatFile[_MAX_PATH]; // optional BED feature or gene file if generating transcriptome reads

double Artef5Rate;			// rate (0..1) at which to insert 5' artefact sequences
int NumArtef5Seqs;			// number of user specified 5' artefact sequences
char *pszArtef5Seqs[cMaxArtefSeqs]; // 5' artefact sequences
double Artef3Rate;			// rate (0..1) at which to insert 3' artefact sequences
int NumArtef3Seqs;			// number of user specified 3' artefact sequences
char *pszArtef3Seqs[cMaxArtefSeqs]; // 5' artefact sequences


char szSQLiteDatabase[_MAX_PATH];	// results summaries to this SQLite file
char szExperimentName[cMaxDatasetSpeciesChrom+1];			// experiment name
char szExperimentDescr[1000];		// describes experiment

// command line args
struct arg_lit  *help    = arg_lit0("h","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *pmode = arg_int0("m","mode","<int>",		    "processing mode: 0 - random start and end, 1 - Profiled start with random end sites, 2 - random start with profiled end sites,  3 - profiled start with profiled end sites (0 - default)");
struct arg_int *fmode = arg_int0("M","format","<int>",		    "output format: 0 - multifasta wrapped, 1 - multifasta non-wrapped");
struct arg_int *numreads = arg_int0("n","nreads","<int>",	    "number of reads required, minimum 5, maximum 500000000 (default = 10000000)");


struct arg_int *region    = arg_int0("G","genomicregion","<int>","Process regions 0:ALL,1:CDS,2:5'UTR,3:3'UTR,4:Introns,5:5'US,6:3'DS,7:Intergenic (default = ALL)");
struct arg_int *updnreglen    = arg_int0("U","updnreglen","<int>","Up/Dn stream regulatory region length, used if processing regions (default is 2000)");


struct arg_int *generrmode = arg_int0("g","generrmode","<int>", "simulate sequencer error modes: 0 - no errors, 1 - induce fixed num errors per read, 2 - static profile, 3 - dynamic according to '-z<rate>' (0 - default)");
struct arg_dbl *seqerrrate = arg_dbl0("z","seqerrs","<dbl>",	"simulate composite sequencer errors with induced error mean rate in range 0-0.20%");
struct arg_lit *seqerrprof = arg_lit0("Z","seqerrprofile",      "generate composite sequencer errors with uniform profile (default is Illumina 3' skewed)");

struct arg_int *indelsize = arg_int0("x","indelsize","<int>",   "simulate micro-InDel size range: 1..9 of this length maximum (default is 3)");
struct arg_dbl *indelrate = arg_dbl0("X","indelrate","<dbl>",	"simulate micro-InDels with mean rate per read in range 0 - 100% of all reads (default is 0)");

struct arg_str *strand=arg_str0("s", "strand","<str>",          "generate for this strand '+' or '-' only (default is both)");
struct arg_int *readlen = arg_int0("l","length","<int>",	    "read lengths (default = 100, max is 100000)");
struct arg_int *cutmin = arg_int0("c","cutmin","<int>",		    "min cut length (minimum = read length)");
struct arg_int *cutmax = arg_int0("C","cutmax","<int>",		    "max cut length (maximum = 100000)");
struct arg_file *infile = arg_file1("i","in","<file>",			"input from this raw multifasta or bioseq assembly");
struct arg_file *inmnase = arg_file0("I","inprofile","<file>",	"input from this profile site preferences file");
struct arg_int *distcluster = arg_int0("D","distcluster","<int>","distribute generated reads as clusters into windows of this median length (0..300), 0 if no clustering");
struct arg_file *featfile = arg_file0("t","featfile","<file>",	"use features or genes in this BED file to generate transcriptome reads from target genome");
struct arg_file *outfile = arg_file1("o","out","<file>",		"output simulated (or N/1 if paired) reads to this file");
struct arg_file *outpefile = arg_file0("O","outpe","<file>",	"output simulated (N/2) paired end reads to this file");
struct arg_file *outsnpfile = arg_file0("u","outsnp","<file>",	"output simulated SNP loci to this BED file, if no SNP rate specified then defaults to 1000 per Mbp");
struct arg_int *threads = arg_int0("T","threads","<int>",		"number of processing threads 0..128 (defaults to 0 which sets threads to number of CPU cores)");
struct arg_lit  *dedupe = arg_lit0("d","dedupe",                "generate unique read sequences only");
struct arg_int *snprate = arg_int0("N","snprate","<int>",       "generate SNPs at the specified rate per Mb (default = 0, range 1..20000)");

struct arg_lit  *pegen = arg_lit0("p","pegen",				    "generate paired end reads");
struct arg_int *pemin = arg_int0("j","pemin","<int>",           "generate paired end reads with minimum fragment lengths (default is 2x read length)");
struct arg_int *pemax = arg_int0("J","pemax","<int>",           "generate paired end reads with maximum  fragment lengths  (default is min fragment length plus read length)");

struct arg_dbl *proprandreads = arg_dbl0("R","randreads","<dbl>","proportion (0 to 0.9000) of reads to be generated with random sequences not likely to be in target genome (default = 0.0)");

struct arg_dbl *artif5rate = arg_dbl0("a","artif5rate","<dbl>",	"randomly induce sequencer adaptor/linker artefacts at 5' end of sequences rate (default 0.0, max 0.9)");
struct arg_str *artif5str = arg_strn("A","artif5str","<string>",0,cMaxArtefSeqs,"5' artefacts sequence(s) (default is 'ACACTCTTTCCCTACACGACGCTGTTCCATCT'");

struct arg_dbl *artif3rate = arg_dbl0("b","artif3rate","<dbl>",	"randomly induce sequencer adaptor/linker artefacts at 3' end of sequences rate (default 0.0, max 0.9)");
struct arg_str *artif3str = arg_strn("A","artif3str","<string>",0,cMaxArtefSeqs,"3' artefacts sequences (default is 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT'");


struct arg_file *summrslts = arg_file0("q","sumrslts","<file>",		"Output results summary to this SQLite3 database file");
struct arg_str *experimentname = arg_str0("w","experimentname","<str>",		"experiment name SQLite3 database file");
struct arg_str *experimentdescr = arg_str0("W","experimentdescr","<str>",	"experiment description SQLite3 database file");

struct arg_end *end = arg_end(200);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					pmode,region,updnreglen,pegen,pemin,pemax,fmode,numreads,proprandreads,snprate,generrmode,distcluster,seqerrrate,seqerrprof,
					artif5rate,artif5str,artif3rate,artif3str,
					indelsize,indelrate,strand,readlen,cutmin,cutmax,dedupe,featfile,
					infile,inmnase,outpefile,outfile,outsnpfile,summrslts,
					experimentname,experimentdescr,
					threads,
					end};

char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

/* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
		{
		printf("\n%s %s %s, Version %s\nOptions ---\n", gszProcName,gpszSubProcess->pszName,gpszSubProcess->pszFullDescr,kit4bversion);
		arg_print_syntax(stdout,argtable,"\n");
		arg_print_glossary(stdout,argtable,"  %-25s %s\n");
		printf("\nNote: Parameters can be entered into a parameter file, one parameter per line.");
		printf("\n      To invoke this parameter file then precede its name with '@'");
		printf("\n      e.g. %s %s @myparams.txt\n",gszProcName,gpszSubProcess->pszName);
		printf("\nPlease report any issues regarding usage of %s at https://github.com/kit4b/issues\n\n",gszProcName);
		return(1);
		}

	/* special case: '--version' takes precedence error reporting */
if (version->count > 0)
		{
		printf("\n%s %s Version %s\n",gszProcName,gpszSubProcess->pszName,kit4bversion);
		return(1);
		}

if (!argerrors)
	{
	if(FileLogLevel->count && !LogFile->count)
		{
		printf("\nError: FileLogLevel '-f%d' specified but no logfile '-F<logfile>'",FileLogLevel->ival[0]);
		exit(1);
		}

	iScreenLogLevel = iFileLogLevel = FileLogLevel->count ? FileLogLevel->ival[0] : eDLInfo;
	if(iFileLogLevel < eDLNone || iFileLogLevel > eDLDebug)
		{
		printf("\nError: FileLogLevel '-l%d' specified outside of range %d..%d",iFileLogLevel,eDLNone,eDLDebug);
		exit(1);
		}
	if(LogFile->count)
		{
		strncpy(szLogFile,LogFile->filename[0],_MAX_PATH);
		szLogFile[_MAX_PATH-1] = '\0';
		}
	else
		{
		iFileLogLevel = eDLNone;
		szLogFile[0] = '\0';
		}

	// now that log parameters have been parsed then initialise diagnostics log system
	if(!gDiagnostics.Open(szLogFile,(etDiagLevel)iScreenLogLevel,(etDiagLevel)iFileLogLevel,true))
		{
		printf("\nError: Unable to start diagnostics subsystem\n");
		if(szLogFile[0] != '\0')
			printf(" Most likely cause is that logfile '%s' can't be opened/created\n",szLogFile);
		return(1);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Subprocess %s Version %s starting",gpszSubProcess->pszName,kit4bversion);
	gExperimentID = 0;
	gProcessID = 0;
	gProcessingID = 0;
	szSQLiteDatabase[0] = '\0';
	szExperimentName[0] = '\0';
	szExperimentDescr[0] = '\0';

	if(experimentname->count)
		{
		strncpy(szExperimentName,experimentname->sval[0],sizeof(szExperimentName));
		szExperimentName[sizeof(szExperimentName)-1] = '\0';
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
	if(summrslts->count)
		{
		strncpy(szSQLiteDatabase,summrslts->filename[0],sizeof(szSQLiteDatabase)-1);
		szSQLiteDatabase[sizeof(szSQLiteDatabase)-1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(szSQLiteDatabase);
		if(strlen(szSQLiteDatabase) < 1)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no SQLite database specified with '-q<filespec>' option");
			return(1);
			}

		if(strlen(szExperimentName) < 1)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no SQLite experiment name specified with '-w<str>' option");
			return(1);
			}
		if(experimentdescr->count)
			{
			strncpy(szExperimentDescr,experimentdescr->sval[0],sizeof(szExperimentDescr)-1);
			szExperimentDescr[sizeof(szExperimentDescr)-1] = '\0';
			CUtility::TrimQuotedWhitespcExtd(szExperimentDescr);
			}
		if(strlen(szExperimentDescr) < 1)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no SQLite experiment description specified with '-W<str>' option");
			return(1);
			}

		gExperimentID = gSQLiteSummaries.StartExperiment(szSQLiteDatabase,false,true,szExperimentName,szExperimentName,szExperimentDescr);
		if(gExperimentID < 1)
			return(1);
		gProcessID = gSQLiteSummaries.AddProcess((char *)gpszSubProcess->pszName,(char *)gpszSubProcess->pszName,(char *)gpszSubProcess->pszFullDescr);
		if(gProcessID < 1)
			return(1);
		gProcessingID = gSQLiteSummaries.StartProcessing(gExperimentID,gProcessID,(char *)kit4bversion);
		if(gProcessingID < 1)
			return(1);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Initialised SQLite database '%s' for results summary collection",szSQLiteDatabase);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"SQLite database experiment identifier for '%s' is %d",szExperimentName,gExperimentID);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"SQLite database process identifier for '%s' is %d",(char *)gpszSubProcess->pszName,gProcessID);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"SQLite database processing instance identifier is %d",gProcessingID);
		}
	else
		{
		szSQLiteDatabase[0] = '\0';
		szExperimentDescr[0] = '\0';
		}

	szFeatFile[0] = '\0';
	szSNPFile[0] = '\0';

	PMode = (etSRPMode)(pmode->count ? pmode->ival[0] : eSRPMStandard);
	if(PMode < eSRPMStandard || PMode >= eSRPMplaceholder)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Processing mode '-m%d' specified outside of range %d..%d",PMode,(int)eSRPMStandard,(int)eSRPMplaceholder-1);
		exit(1);
		}

	Region = eSRMEGRAny;
	UpDnRegLen = 0;

	InDelRate = 0.0f;
	InDelSize = 0;

	bPEgen = pegen->count ? true : false;
	if(bPEgen && PMode != eSRPMStandard)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Processing mode '-m%d' not supported when generating paired ends",PMode);
		exit(1);
		}

	if(bPEgen && outpefile->count == 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Paired end processing '-p' requested but no partner paired end output file specifed with '-O<file>'");
		exit(1);
		}

	FMode = (etSRFMode)(fmode->count ? fmode->ival[0] : eSRFMFasta);
	if(FMode < eSRFMFasta || FMode >= eSRFMplaceholder)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Output format '-m%d' specified outside of range %d..%d",PMode,(int)eSRFMFasta,(int)eSRFMplaceholder-1);
		exit(1);
		}

	if(bPEgen && !(FMode == eSRFMNWFasta || FMode == eSRFMFasta))
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Output format '-M%d' not supported (only fasta output supported) when generating paired ends",FMode);
		exit(1);
		}

	PropRandReads = proprandreads->count ? proprandreads->dval[0] : 0.0;
	if(PropRandReads < 0.0 || PropRandReads >= 0.9000)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Proportion of random reads '-R%f' must in range 0 to 0.9000",PropRandReads);
		exit(1);
		}

	SNPrate = snprate->count ? snprate->ival[0] : 0;
	if(SNPrate < 0 || SNPrate > 20000)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: simulated SNP rate '-N%d' specified outside of range 0..20000",SNPrate);
		exit(1);
		}

	if(outsnpfile->count)
		{
		strncpy(szSNPFile,outsnpfile->filename[0],_MAX_PATH);
		szSNPFile[_MAX_PATH-1] = '\0';
		if(SNPrate == 0)
			{
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Info: SNP file output requested but no SNP rate specified, defaulting SNP rate to 1000 per Mbp");
			SNPrate = 1000;
			}
		}
	else
		szSNPFile[0] = '\0';

	DistCluster = distcluster->count ? distcluster->ival[0] : 0;
	if(DistCluster < 0 || DistCluster > cMaxDistClusterLen)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: distributed cluster size '-D%d' must be in range 0 to %d",DistCluster,cMaxDistClusterLen);
		exit(1);
		}

	InDelRate = indelrate->count ? indelrate->dval[0] : 0.0f;
	if(InDelRate < 0.0 || InDelRate > 1.0)
		{
		printf("\nError: simulated InDel size range '-X%f' must be in range 0.0 to 1.0",InDelRate);
		exit(1);
		}
	if(InDelRate > 0.0)
		{
		InDelSize = indelsize->count ? indelsize->ival[0] :3;
		if(InDelSize < 1 || InDelSize > 9)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: simulated InDel size range '-x%d' specified outside of range 1..9",InDelSize);
			exit(1);
			}
		}
	else
		InDelSize = 0;

	SEMode = (etSRSEMode)(generrmode->count ? generrmode->ival[0] : eSRSEPnone);
	if(SEMode < eSRSEPnone || SEMode >= eSRSEPplaceholder)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: simulated error rate mode '-m%d' specified outside of range %d..%d",SEMode,(int)eSRSEPnone,(int)eSRSEPplaceholder-1);
		exit(1);
		}
	if(SEMode == eSRSEPdyn || SEMode == eSRSEPfixerrs)
		{
		SeqErrRate = SEMode == eSRSEPdyn ? 0.01 : 5.0;
		double MaxSeqErrRate = SEMode == eSRSEPdyn ? 0.20 : 30;
		SeqErrRate = seqerrrate->count ? seqerrrate->dval[0] : SeqErrRate;
		if(SeqErrRate < 0.0 || SeqErrRate > MaxSeqErrRate)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: sequencer error rate '-z%f' must be in range 0.0 to %f",SeqErrRate,MaxSeqErrRate);
			exit(1);
			}
		}
	else
		SeqErrRate = -1;

	if((FMode < eSRFMNWFasta) && SEMode != eSRSEPnone)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: sequencer error rate mode '-g%d' only allowed if generating multifasta output with '-M2/3/4'",SEMode);
		exit(1);
		}

	bSeqErrProfile = seqerrprof->count ? true : false;
	if(bSeqErrProfile && SEMode == eSRSEPnone)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Uniform profile '-Z' requested but no error rate mode specified with '-g<mode>'");
		exit(1);
		}

	Strand = strand->count ? *(char *)strand->sval[0] : '*';
	if(!(Strand == '+' || Strand == '-' || Strand == '*'))
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Strand specified '-s%c' must be one of '+', '-' or '*'",Strand);
		exit(1);
		}

	NumReads = numreads->count ? numreads->ival[0] : cSRDfltNumReads;
	if(NumReads < cSRMinNumReads || NumReads > cSRMaxNumReads)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Number of reads '-n%d' specified outside of range %d..%d",NumReads,cSRMinNumReads,cSRMaxNumReads);
		exit(1);
		}

	ReadLen = readlen->count ? readlen->ival[0] : cSRDfltReadLen;
	if(ReadLen < cSRMinReadLen || ReadLen > cSRMaxReadLen)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Read length '-a%d' specified outside of range %d..%d",ReadLen,cSRMinReadLen,cSRMaxReadLen);
		exit(1);
		}

		CutMin = cutmin->count ? cutmin->ival[0] : ReadLen;
		if(CutMin < ReadLen || CutMin > cSRMaxCutLen)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Minimum cut length '-c%d' must be in range %d..%d",CutMin,ReadLen,cSRMaxCutLen);
			exit(1);
			}

		CutMax = cutmax->count ? cutmax->ival[0] : CutMin;
		if(CutMax < CutMin || CutMax > cSRMaxCutLen)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Maximum cut length '-C%d' must be in range %d..%d",CutMax,CutMin,cSRMaxCutLen);
			exit(1);
			}


		if(bPEgen)
			{
			PEmin = pemin->count ? pemin->ival[0] : ReadLen * 2;
			if(PEmin < cMinPEFragLen || PEmin > cMaxPEFragLen)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: paired end minimum fragment length '-j%d' specified outside of range %d..%d",PEmin,cMinPEFragLen,cMaxPEFragLen);
				exit(1);
				}
			PEmax = min(PEmin + ReadLen,cMaxPEFragLen);
			PEmax = pemax->count ? pemax->ival[0] : PEmax;
			if(PEmax < PEmin || PEmax > cMaxPEFragLen)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: paired end maximum fragment length '-J%d' specified outside of range %d..%d",PEmax,PEmin,cMaxPEFragLen);
				exit(1);
				}

			if(PEmin <= CutMin)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: paired end min fragment length (%d) must be more than min read length (%d)",PEmin,CutMin);
				exit(1);
				}

			}
		else
			{
			PEmin = 0;
			PEmax = 0;
			}

	bDedupe = dedupe->count ? true : false;

	if(bPEgen && bDedupe)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Deduping not supported when generating paired ends");
		exit(1);
		}

	if(PMode != eSRPMStandard && inmnase->count)
		{
		strncpy(szProfFile,inmnase->filename[0],_MAX_PATH);
		szProfFile[_MAX_PATH-1] = '\0';
		}
	else
		szProfFile[0] = '\0';

	strcpy(szInFile,infile->filename[0]);

	strcpy(szOutFile,outfile->filename[0]);

	if(featfile->count)
		{
		strncpy(szFeatFile,featfile->filename[0],_MAX_PATH);
		szFeatFile[_MAX_PATH-1] = '\0';
		FMode = eSRFMFasta;
		}
	else
		szFeatFile[0] = '\0';

	if(Region != eSRMEGRAny && szFeatFile[0] == '\0')
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Filtering on region specified with '-G%d' but no BED feature file specified with '-t<file>'",Region);
		exit(1);
		}

	if(bPEgen)
		strcpy(szOutPEFile,outpefile->filename[0]);
	else
		szOutPEFile[0] = '\0';


	int MaxArtefLen = min(ReadLen - 5,cMaxArtefSeqLen);
	if(artif5rate->count ||  artif3rate->count)
		{
		if(!(FMode == eSRFMNWFasta || FMode == eSRFMFasta ) || bPEgen || Region != eSRMEGRAny || bDedupe)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Sorry artefact sequence processing not supported with currently selected options - only single end fasta is supported!");
			exit(1);
			}
		}


	if(FMode == eSRFMNWFasta || FMode == eSRFMFasta)
		{
		Artef5Rate = artif5rate->count ? artif5rate->dval[0] : 0.0;
		if(Artef5Rate < 0.0 || Artef5Rate > 1.0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: 5' artefact rate '-a%1.2f' must be in the range 0.0 to 1.0",Artef5Rate);
			exit(1);
			}
		}
	else
		Artef5Rate = 0.0;

	if(Artef5Rate > 0.0)
		{
		NumArtef5Seqs = artif5str->count;
		if(NumArtef5Seqs > 0)
			{
			for(Idx=0;Idx < artif5str->count; Idx++)
				{
				LenReq = (int)strlen(artif5str->sval[Idx]);
				pszArtef5Seqs[Idx] = new char [LenReq+1];
				strcpy(pszArtef5Seqs[Idx],artif5str->sval[Idx]);
				if((LenReq = CSimReads::TrimSeqChrs(pszArtef5Seqs[Idx])) < 1 || LenReq > MaxArtefLen)
					{
					if(LenReq == -1)
						gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: 5' artefact sequence '-A\"%s\"' contains non-sequence char",artif5str->sval[Idx]);
					else
						{
						if(LenReq == 0)
							gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: 5' artefact sequence '-A\"%s\"' after filtering is empty",artif5str->sval[Idx]);
						else
							gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: 5' artefact sequence '-A\"%s\"' is too long (must be <= %d bp)",pszArtef5Seqs[Idx],MaxArtefLen);
						}

					exit(1);
					}
				}
			}
		else
			{
			pszArtef5Seqs[0] = new char [cMaxArtefSeqLen+1];
			strncpy(pszArtef5Seqs[0], cpszArtef5Seq,MaxArtefLen);
			pszArtef5Seqs[0][MaxArtefLen] = '\0';
			NumArtef5Seqs = 1;
			}
		}
	else
		{
		NumArtef5Seqs = 0;
		pszArtef5Seqs[0] = nullptr;
		}

	if(FMode == eSRFMNWFasta || FMode == eSRFMFasta)
		{
		Artef3Rate = artif3rate->count ? artif3rate->dval[0] : 0.0;
		if(Artef3Rate < 0.0 || Artef3Rate > 1.0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: 5' artefact rate '-b%1.2f' must be in the range 0.0 to 1.0",Artef3Rate);
			exit(1);
			}
		}
	else
		Artef3Rate = 0.0;



	if(Artef3Rate > 0.0)
		{
		NumArtef3Seqs = artif3str->count;
		if(NumArtef3Seqs > 0)
			{
			for(Idx=0;Idx < artif3str->count; Idx++)
				{
				LenReq = (int)strlen(artif3str->sval[Idx]);
				pszArtef3Seqs[Idx] = new char [LenReq+1];
				strcpy(pszArtef3Seqs[Idx],artif3str->sval[Idx]);
				if((LenReq = CSimReads::TrimSeqChrs(pszArtef3Seqs[Idx])) < 1 || LenReq > MaxArtefLen)
					{
					if(LenReq == -1)
						gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: 3' artefact sequence '-A\"%s\"' contains non-sequence char",artif3str->sval[Idx]);
					else
						{
						if(LenReq == 0)
							gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: 3' artefact sequence '-A\"%s\"' after filtering is empty",artif3str->sval[Idx]);
						else
							gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: 3' artefact sequence '-A\"%s\"' is too long (must be <= %d bp)",pszArtef3Seqs[Idx],MaxArtefLen);
						}

					exit(1);
					}
				}
			}
		else
			{
			pszArtef3Seqs[0] = new char [cMaxArtefSeqLen+1];
			strncpy(pszArtef3Seqs[0], cpszArtef3Seq,MaxArtefLen);
			pszArtef3Seqs[0][MaxArtefLen] = '\0';
			NumArtef3Seqs = 1;
			}
		}
	else
		{
		NumArtef3Seqs = 0;
		pszArtef3Seqs[0] = nullptr;
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
	int MaxAllowedThreads = min(cMaxWorkerThreads,NumberOfProcessors);	// limit to be at most cMaxWorkerThreads
	if((NumThreads = threads->count ? threads->ival[0] : MaxAllowedThreads)==0)
		NumThreads = MaxAllowedThreads;
	if(NumThreads < 0 || NumThreads > MaxAllowedThreads)
		{
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"Warning: Number of threads '-T%d' specified was outside of range %d..%d",NumThreads,1,MaxAllowedThreads);
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"Warning: Defaulting number of threads to %d",MaxAllowedThreads);
		NumThreads = MaxAllowedThreads;
		}

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing parameters:");

	const char *pszDescr;
	switch(PMode) {
		case eSRPMStandard:
			pszDescr = "random start and random end sites";
			break;
		case eSRPMProfRand:
			pszDescr = "profiled start with random end sites";
			break;
		case eSRPMRandProf:
			pszDescr = "profiled start with random end sites";
			break;
		case eSRPMProfProf:
			pszDescr = "profiled start with profiled end sites";
			break;
		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"processing mode is : '%s'",pszDescr);

	switch(FMode) {
		case eSRFMFasta:
			pszDescr = "Multifasta, wrapped read sequences";
			break;
		case eSRFMNWFasta:
			pszDescr = "Multifasta, non-wrapped read sequences";
			break;
		}


	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Generate reads as : '%s'", bPEgen ? "paired ends" : "single ended");
	if(bPEgen)
		gDiagnostics.DiagOutMsgOnly(eDLInfo," Paired end fragment length: min %d, max %d", PEmin, PEmax);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output format is : '%s'",pszDescr);

	switch(SEMode) {
		case eSRSEPnone:
			pszDescr = "No induced simulated errors";
			break;
		case eSRSEPfixerrs:
			pszDescr = "simulate fixed number of errors in each read at rate specified by '-z<rate>'";
			break;
		case eSRSEPstatic:
			pszDescr = "use internal static profile";
			break;
		case eSRSEPdyn:
			pszDescr = "dynamic according to '-z<rate>'";
			break;
		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Simulated sequencer error mode : '%s'",pszDescr);

	if(SNPrate > 0)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Simulate SNPs at this rate per million bases : %d",SNPrate);
	if(szSNPFile[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Write simulated SNP loci to this BED file : '%s'",szSNPFile);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Generate random reads in this Proportion: %0.3f",PropRandReads);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Cluster reads into these sized windows : %d", DistCluster);

	if(SEMode == eSRSEPdyn)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Induce composite sequence errors at mean rate of : %f",SeqErrRate);
	else
		if(SEMode == eSRSEPfixerrs)
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Induce fixed number of sequence errors : %d",(int)SeqErrRate);
	if(SEMode != eSRSEPnone)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Induce sequence errors profile : '%s'",bSeqErrProfile ? "uniform" : "3' skewed");

	if(InDelRate == 0.0f)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Do not induce InDels");
	else
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Induce InDels at mean rate per read of : %1.3f",InDelRate);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"InDels will have a size range of : 1 to %d",InDelSize);
		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Simulate for this strand : '%c'",Strand);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Number of %s reads required: %d",bPEgen ?  "paired end" : "single ended", NumReads);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Read lengths: %d",ReadLen);

	if(Artef5Rate > 0.0 && NumArtef5Seqs > 0)
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Induce 5' artefacts at this rate: %1.2f",Artef5Rate);
		for(Idx = 0; Idx < NumArtef5Seqs; Idx++)
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"   5' artefact (%d) : '%s'",Idx+1,pszArtef5Seqs[Idx]);
		}
	else
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Induce 5' artefacts at this rate: %1.2f",Artef5Rate);

	if(Artef3Rate > 0.0 && NumArtef3Seqs > 0)
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Induce 3' artefacts at this rate: %1.2f",Artef3Rate);
		for(Idx = 0; Idx < NumArtef3Seqs; Idx++)
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"   3' artefact (%d) : '%s'",Idx+1,pszArtef3Seqs[Idx]);
		}
	else
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Induce 3' artefacts at this rate: %1.2f",Artef3Rate);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"min cut length: %d",CutMin);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"max cut length: %d",CutMax);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"unique read sequences only: %s",bDedupe ? "yes" : "no");


	gDiagnostics.DiagOutMsgOnly(eDLInfo,"input bioseq assembly file: '%s'",szInFile);
	if(PMode != eSRPMStandard)
		{
		if(szProfFile[0] == '\0')
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"static profiled site preferences derived from site octamer");
		else
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"input profiled site preferences file: '%s'",szProfFile);
		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"output simulated reads file: '%s'",szOutFile);
	if(bPEgen)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"output simulated paired read partners to file: %s", szOutPEFile);

	if(szFeatFile[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"input BED feature or gene file: '%s'",szFeatFile);

	if(Region != eSRMEGRAny)
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Process for region: %s",CSimReads::Region2Txt(Region));
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Up/Dn stream regulatory length: %d",UpDnRegLen);
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"number of threads : %d",NumThreads);

	if(gExperimentID > 0)
		{
		int ParamID;
		int Idy;
		int bTrueInt = 1;
		int bFalseInt = 0;

		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTText,(int)strlen(szLogFile),"log",szLogFile);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,sizeof(PMode),"mode",&PMode);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,sizeof(FMode),"format",&FMode);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTText,1,"strand",&Strand);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,sizeof(NumReads),"nreads",&NumReads);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,sizeof(ReadLen),"length",&ReadLen);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,sizeof(CutMin),"cutmin",&CutMin);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,sizeof(CutMax),"cutmax",&CutMax);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,sizeof(bTrueInt),"dedupe",bDedupe  == true ? &bTrueInt : &bFalseInt);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,sizeof(NumThreads),"threads",&NumThreads);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,sizeof(SNPrate),"snprate",&SNPrate);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,sizeof(Region),"genomicregion",&Region);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,sizeof(UpDnRegLen),"updnreglen",&UpDnRegLen);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTDouble,sizeof(PropRandReads),"randreads",&PropRandReads);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,sizeof(DistCluster),"distcluster",&DistCluster);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTDouble,sizeof(SeqErrRate),"mode",&SeqErrRate);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,sizeof(InDelSize),"indelsize",&InDelSize);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,sizeof(bTrueInt),"seqerrprofile",bSeqErrProfile == true ? &bTrueInt : &bFalseInt);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,sizeof(PMode),"indelrate",&InDelRate);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,sizeof(bTrueInt),"pegen",bPEgen == true ? &bTrueInt : &bFalseInt);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,sizeof(PEmin),"pemin",&PEmin);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,sizeof(PEmax),"pemax",&PEmax);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTDouble,sizeof(Artef5Rate),"artif5rate",&Artef5Rate);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,sizeof(NumArtef5Seqs),"NumArtef5Seqs",&NumArtef5Seqs);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTDouble,sizeof(Artef3Rate),"artif3rate",&Artef3Rate);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,sizeof(NumArtef3Seqs),"NumArtef3Seqs",&NumArtef3Seqs);

		for(Idy = 0; Idy < NumArtef5Seqs; Idy++)
			ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTText,(int)strlen(pszArtef5Seqs[Idy]),"artif5str",pszArtef5Seqs[Idy]);
		for(Idy = 0; Idy < NumArtef3Seqs; Idy++)
			ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTText,(int)strlen(pszArtef3Seqs[Idy]),"artif3str",pszArtef3Seqs[Idy]);

		if(szInFile[0] != '\0')
			ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTText,(int)strlen(szInFile),"in",szInFile);
		if(szProfFile[0] != '\0')
			ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTText,(int)strlen(szProfFile),"inprofile",szProfFile);
		if(szOutFile[0] != '\0')
			ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTText,(int)strlen(szOutFile),"out",szOutFile);
		if(szOutPEFile[0] != '\0')
			ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTText,(int)strlen(szOutPEFile),"outpe",szOutPEFile);
		if(szSNPFile[0] != '\0')
			ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTText,(int)strlen(szSNPFile),"outsnp",szSNPFile);
		if(szFeatFile[0] != '\0')
			ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTText,(int)strlen(szFeatFile),"featfile",szFeatFile);

		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,sizeof(NumThreads),"threads",&NumThreads);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,sizeof(NumberOfProcessors),"cpus",&NumberOfProcessors);

		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTText,(int)strlen(szSQLiteDatabase),"sumrslts",szSQLiteDatabase);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTText,(int)strlen(szExperimentName),"experimentname",szExperimentName);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTText,(int)strlen(szExperimentDescr),"experimentdescr",szExperimentDescr);
		}

#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = Process((etSRPMode)PMode,SEMode,bPEgen,PEmin,PEmax,PropRandReads,
			DistCluster,SeqErrRate,bSeqErrProfile,SNPrate,
						InDelSize,InDelRate,FMode,NumThreads,Strand,NumReads,
						ReadLen,Artef5Rate,NumArtef5Seqs,pszArtef5Seqs,Artef3Rate,NumArtef3Seqs,pszArtef3Seqs,
						CutMin,CutMax,bDedupe,Region,UpDnRegLen,szFeatFile,szInFile,szProfFile,szOutPEFile,szOutFile,szSNPFile);
	Rslt = Rslt >=0 ? 0 : 1;
	if(gExperimentID > 0)
		{
		if(gProcessingID)
			gSQLiteSummaries.EndProcessing(gExperimentID, gProcessingID,Rslt);
		gSQLiteSummaries.EndExperiment(gExperimentID);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	return(Rslt);
	}
else
	{
	printf("\n%s %s %s, Version %s\n", gszProcName,gpszSubProcess->pszName,gpszSubProcess->pszFullDescr,kit4bversion);
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nUse '-h' to view option and parameter usage\n");
	return(1);
	}
return 0;
}

int
Process(etSRPMode PMode,		// processing mode
		etSRSEMode SEMode,	// induced sequencer error rate mode
		bool bPEgen,		// true if paired ends are to be generated
		int PEmin,			// PE minimum insert
		int PEmax,			// PE maximum insert
		double PropRandReads, // generate completely random reads at this rate
		int DistCluster,	// cluster generated reads into windows of this median length, 0 if no clustering
		double SeqErrRate,	// dynamic sequencing errors are to be induced into generated sequences at this rate
		bool bSeqErrProfile,// true if to generate composite sequencer errors with uniform profile (default is Illumina 3' skewed)
		int SNPrate,		// simulate SNPs at this rate per million bases
		int InDelSize,		// simulated InDel size range
		double InDelRate,	// simulated InDel rate per read
		etSRFMode FMode,		// output format
		int NumThreads,		// number of worker threads to use
		char Strand,		// generate for this strand '+' or '-' or for both '*'
		int NumReads,		// number of reads required (will be doubled if paired end reads)
		int ReadLen,		// read lengths
		double Artef5Rate,			// rate (0..1) at which to insert 5' artefact sequences
		int NumArtef5Seqs,			// number of user specified 5' artefact sequences
		char *pszArtef5Seqs[], // 5' artefact sequences
		double Artef3Rate,			// rate (0..1) at which to insert 3' artefact sequences
		int NumArtef3Seqs,			// number of user specified 3' artefact sequences
		char *pszArtef3Seqs[], // 5' artefact sequences
		int CutMin,			// min cut length
		int CutMax,			// max cut length
		bool bDedupe,		// true if unique read sequences only to be generated
		int Region,			// Process regions 0:ALL,1:CDS,2:5'UTR,3:3'UTR,4:Introns,5:5'US,6:3'DS,7:Intergenic (default = ALL)
		int UpDnRegLen,		// if processing regions then up/down stream regulatory length
		char *pszFeatFile,	// optionally generate transcriptome reads from features or genes in this BED file
		char *pszInFile,	// input from this bioseq assembly
		char *pszProfFile,	// input from this profile site preferences file
		char *pszOutPEFile, // output partner paired end simulated reads to this file
		char *pszOutFile,	// output simulated reads to this file
		char *pszOutSNPs)   // output simulated SNP loci to this file
{
int Rslt;
CSimReads *pSimReads;

if((pSimReads = new CSimReads) == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate instance of CSimReads");
	return(-1);
	}

Rslt = pSimReads->GenSimReads(PMode, SEMode, bPEgen, PEmin, PEmax, PropRandReads, DistCluster,
		SeqErrRate,	bSeqErrProfile,	SNPrate, InDelSize,	InDelRate, FMode,
		NumThreads,	Strand,	NumReads, ReadLen,	Artef5Rate,	NumArtef5Seqs, pszArtef5Seqs,Artef3Rate,NumArtef3Seqs,	
		pszArtef3Seqs,	CutMin,	CutMax,	bDedupe,Region,	UpDnRegLen,	pszFeatFile,pszInFile,pszProfFile,pszOutPEFile,	pszOutFile,	pszOutSNPs);
delete pSimReads;
return(Rslt);
}


