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

// kingsax.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#if _WIN32
#include "../libkit4b/commhdrs.h"
#else

#include "../libkit4b/commhdrs.h"
#endif

#include "ngskit4b.h"

const size_t cGbpSeqSize = 1000000000;	// define our definition of a Gbp
const int cMaxSimGenomeGbp = 1000;      // when simulating genomes then allowing for genomes of upto 1 Tbp which would require ~ 6TB of memory!
const int cNumSimChroms = 10000;		// simulated genomes will be split over this many equal sized chromosomes
const int cMinBlkEls = 1;				// user can specify suffix blocks down to this minimum size (GB)
const int cDlftBlkEls = 0;				// default depends on available memory and will default to be at most 75% of physical memory
const int cMaxBlkEls = 40;				// user can specify suffix blocks up to this maximum size (GB)

const int cMaxInFileSpecs = 100;			// user can specify upto this many input files

const int cMaxDupEntries = 10;				// report 1st 10 duplicate entry names

const int cMaxAXAllocBuffChunk = 0x00ffffff;	// buffer for fasta sequences is realloc'd in this sized chunks

int
CreateBioseqSuffixFile(int Mode,
						int MaxKMers,			// downstream this suffix array will only be used for processing of these sized maximal K-mers, only advantage is that shorter k-mers are quicker to sort than longer k-mers
						int MinSeqLen,			// only accept for indexing sequences which are at least this length
						int SimGenomeSize,		// if 1..120 then simulating indexing of a genome of this size in Gbp.
						int MaxThreads,			// max threads
						bool bSOLiD,				// true if to process for colorspace (SOLiD)
						int NumInputFiles,			// number of input file specs
						char *pszInputFiles[],		// names of input files (wildcards allowed)
						char *pszDestSfxFile,	// output suffix array to this file
						char *pszRefSpecies,	// species assembly which is being indexed
						char *pszDescr,			// describes assembly
						char *pszTitle);        // title to use


CSfxArray *m_pSfxFile;				// suffix array file being created


#ifdef _WIN32
int kingsax(int argc, char* argv[])
{
// determine my process name
_splitpath(argv[0],nullptr,nullptr,gszProcName,nullptr);
#else
int
kingsax(int argc, char** argv)
{
// determine my process name
CUtility::splitpath((char *)argv[0],nullptr,gszProcName);
#endif
int iScreenLogLevel;		// level of screen diagnostics
int iFileLogLevel;			// level of file diagnostics
char szLogFile[_MAX_PATH];	// write diagnostics to this file

int Rslt;

char szOutputFileSpec[_MAX_PATH];
int MinSeqLen;								// only accept for indexing sequences which are at least this length
int SimGenomeSize;							// if 1..1000 then simulating indexing of a genome of this size in Gbp.
int NumInputFileSpecs;						// number of input file specs
char *pszInputFileSpecs[cMaxInFileSpecs];	// input files
char szDescription[cMBSFFileDescrLen];
char szTitle[cMBSFShortFileDescrLen];
char szRefSpecies[cMaxDatasetSpeciesChrom];
int iMode;									// processing mode
bool bSOLiD;								// colorspace (SOLiD) generation
int NumberOfProcessors;						// number of installed CPUs
int NumThreads;								// number of threads (0 defaults to number of CPUs)
int MaxKMers;								// downstream this suffix array will only be used for processing of these sized maximal K-mers, only advantage is that shorter k-mers are quicker to sort than longer k-mers

char szSQLiteDatabase[_MAX_PATH];	// results summaries to this SQLite file
char szExperimentName[cMaxDatasetSpeciesChrom+1];			// experiment name
char szExperimentDescr[1000];		// describes experiment

// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *Mode=arg_int0("m", "mode",	"<int>",			"Processing mode, 0=standard, 1=bisulphite index, 2=simulated genome (default 0)");
struct arg_lit  *solid = arg_lit0("C","colorspace",             "Generate for colorspace (SOLiD)");
struct arg_int *simgenomesize=arg_int0("s", "simgenomesize",	"<int>","Simulated genome size in Gbp (default 5, range 1..1000");

struct arg_int *minseqlen=arg_int0("l", "minseqlen",			"<int>","Do not accept for indexing sequences less than this length (default 50, range 1..1000000)");
struct arg_int *maxkmers = arg_int0("k", "maxkmers", "<int>",	"when sorting compare up to this length KMer for equality (default 100000, range 10..1000000)");

struct arg_file *infiles = arg_filen("i",nullptr,"<file>",0,cMaxInFileSpecs,	"input from wildcarded kingss or fasta files");
struct arg_file *OutFile = arg_file0("o",nullptr,"<file>",			"output suffix array file");
struct arg_str *Descr = arg_str0("d","descr","<string>",		"full description");
struct arg_str *Title = arg_str0("t","title","<string>",		"short title");
struct arg_str *RefSpecies = arg_str1("r","ref","<string>",		"reference species");
struct arg_int *threads = arg_int0("T","threads","<int>",		"number of processing threads 0..128 (defaults to 0 which sets threads to number of CPU cores)");
struct arg_file *summrslts = arg_file0("q","sumrslts","<file>",		"Output results summary to this SQLite3 database file");
struct arg_str *experimentname = arg_str0("w","experimentname","<str>",		"experiment name SQLite3 database file");
struct arg_str *experimentdescr = arg_str0("W","experimentdescr","<str>",	"experiment description SQLite3 database file");
struct arg_end *end = arg_end(200);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					summrslts,experimentname,experimentdescr,
					Mode,maxkmers,minseqlen,simgenomesize,solid,infiles,OutFile,RefSpecies,Descr,Title,
					threads,end};

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
		exit(1);
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

	MaxKMers = 0;
	MinSeqLen = 0;
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

	iMode = Mode->count ? Mode->ival[0] : 0;
	if(iMode < 0 || iMode > 2)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: processing mode '-m%d' must be specified in range %d..%d",iMode,0,1);
		exit(1);
		}

	if(iMode == 2)
		{
		SimGenomeSize = simgenomesize->count ? simgenomesize->ival[0] : 5;
		if(SimGenomeSize < 1 || iMode > cMaxSimGenomeGbp)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: simulated genome size in Gbp '-s%d' must be specified in range 1..%d",SimGenomeSize,cMaxSimGenomeGbp);
			exit(1);
			}
		}
	else
		SimGenomeSize = 0;

	MaxKMers = maxkmers->count ? maxkmers->ival[0] : 0;
	if(MaxKMers == 0)
		MaxKMers = 100000;
	if(MaxKMers < 10 || MaxKMers > cMaxReadLen)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Maximum KMer sort length, -k%d, must be in the range 10..%dbp", MaxKMers, cMaxReadLen);
		exit(1);
		}

	bSOLiD = solid->count ? true : false;

	int Idx;

	if(iMode != 2)
		{
		MinSeqLen = minseqlen->count ? minseqlen->ival[0] : 50;
		if(MinSeqLen < 1)
			MinSeqLen = 1;
		else
			if(MinSeqLen > 1000000)
				MinSeqLen = 1000000;

		for(NumInputFileSpecs=Idx=0;NumInputFileSpecs < cMaxInFileSpecs && Idx < infiles->count; Idx++)
			{
			pszInputFileSpecs[Idx] = nullptr;
			if(pszInputFileSpecs[NumInputFileSpecs] == nullptr)
				pszInputFileSpecs[NumInputFileSpecs] = new char [_MAX_PATH];
			strncpy(pszInputFileSpecs[NumInputFileSpecs],infiles->filename[Idx],_MAX_PATH);
			pszInputFileSpecs[NumInputFileSpecs][_MAX_PATH-1] = '\0';
			CUtility::TrimQuotedWhitespcExtd(pszInputFileSpecs[NumInputFileSpecs]);
			if(pszInputFileSpecs[NumInputFileSpecs][0] != '\0')
				NumInputFileSpecs++;
			}

		if(!NumInputFileSpecs)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, no input file(s) specified with '-i<filespec>' option)\n");
			exit(1);
			}
		}
	else
		{
		NumInputFileSpecs = 0;
		pszInputFileSpecs[0] = nullptr;
		}

	if(OutFile->count)
		{
		strcpy(szOutputFileSpec,OutFile->filename[0]);
		CUtility::TrimQuotedWhitespcExtd(szOutputFileSpec);
		}
	else
		szOutputFileSpec[0] = '\0';

	if(iMode != 2 && szOutputFileSpec[0] == '\0')
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: no output file(s) specified with '-o<filespec>' option)\n");
		exit(1);
		}

	strncpy(szRefSpecies,RefSpecies->sval[0],cMaxDatasetSpeciesChrom);
	szRefSpecies[cMaxDatasetSpeciesChrom-1] = '\0';

	if(!Title->count)
		strcpy(szTitle,szRefSpecies);
	else
		{
		strncpy(szTitle,Title->sval[0],cMBSFShortFileDescrLen);
		szTitle[cMBSFShortFileDescrLen-1] = '\0';
		}

	if(!Descr->count)
		strcpy(szDescription,szRefSpecies);
	else
		{
		strncpy(szDescription,Descr->sval[0],cMBSFFileDescrLen);
		szDescription[cMBSFFileDescrLen-1] = '\0';
		}

// show user current resource limits
#ifndef _WIN32
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Resources: %s",CUtility::ReportResourceLimits());
#endif

	int AvailGBMemory;
#ifdef _WIN32
	SYSTEM_INFO SystemInfo;
	GetSystemInfo(&SystemInfo);
	NumberOfProcessors = SystemInfo.dwNumberOfProcessors;

	MEMORYSTATUSEX status;
	status.dwLength = sizeof(status);
	GlobalMemoryStatusEx(&status);
	AvailGBMemory = (int)(status.ullTotalPhys / 0x040000000);	// set to be gigabytes
#else
	NumberOfProcessors = sysconf(_SC_NPROCESSORS_CONF);

	long pages = sysconf(_SC_PHYS_PAGES);
	long page_size = sysconf(_SC_PAGE_SIZE);
	AvailGBMemory = (int)(((int64_t)pages * (int64_t)page_size)/0x040000000);
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
	gDiagnostics.DiagOutMsgOnly(eDLInfo, "Sorting suffixes of max KMer length : %d", MaxKMers);
	if(iMode == 2)
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Simulating indexing of a genome of size: %dGbp",SimGenomeSize);
		if(AvailGBMemory < (SimGenomeSize * 7))		// allowing for 6 bytes per base plus overheads
			gDiagnostics.DiagOutMsgOnly(eDLWarn,"May be memory allocation problems when loading or indexing this simulated genome, available memory: %dGB",AvailGBMemory);
		}

	if(bSOLiD)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Process for colorspace (SOLiD)");

	if(iMode != 2)
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Accepting for indexing sequences of length at least: %dbp",MinSeqLen);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Process input files as: '%s'",iMode == 0 ? "standard" : "bisulphite");
		for(Idx=0; Idx < NumInputFileSpecs; Idx++)
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input source file spec: '%s'",pszInputFileSpecs[Idx]);
		}

	if(szOutputFileSpec[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output to suffix array file: '%s'",szOutputFileSpec);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Reference species: '%s'",szRefSpecies);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Title text: '%s'",szTitle);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Descriptive text: '%s'",szDescription);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Number of threads : %d",NumThreads);

	if(szExperimentName[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"This processing reference: %s",szExperimentName);

#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = CreateBioseqSuffixFile(iMode, MaxKMers,MinSeqLen,SimGenomeSize,NumThreads,bSOLiD,NumInputFileSpecs,pszInputFileSpecs,szOutputFileSpec,szRefSpecies,szDescription,szTitle);
	Rslt = Rslt >=0 ? 0 : 1;
	if(gExperimentID > 0)
		{
		if(gProcessingID)
			gSQLiteSummaries.EndProcessing(gExperimentID,gProcessingID,Rslt);
		gSQLiteSummaries.EndExperiment(gExperimentID);
		}
	gStopWatch.Stop();
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	exit(Rslt);
	}
else
	{
	printf("\n%s %s %s, Version %s\n", gszProcName,gpszSubProcess->pszName,gpszSubProcess->pszFullDescr,kit4bversion);
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nUse '-h' to view option and parameter usage\n");
	exit(1);
	}
}


void InitpSfx(void)
{
m_pSfxFile = nullptr;
}

void Reset(void)
{
if(m_pSfxFile != nullptr)
	{
	delete m_pSfxFile;
	m_pSfxFile = nullptr;
	}
}



// ProcessFastaFile
// Parse input fasta format file into a biosequence suffix array file
int
ProcessFastaFile(int MinSeqLen,		// only accept for indexing sequences of at least this length
				char *pszFile)
{
CFasta Fasta;
uint8_t *pSeqBuff;
uint8_t *pMskBase;
uint32_t MskIdx;
size_t BuffOfs;
size_t AllocdBuffSize;
size_t AvailBuffSize;
char szName[cBSFSourceSize];
char szDescription[cBSFDescriptionSize];
uint32_t SeqLen;
int Descrlen;
bool bFirstEntry;
bool bEntryCreated;
int Rslt;
int SeqID;
uint32_t NumSeqsUnderlength;

if((Rslt=Fasta.Open(pszFile,true))!=eBSFSuccess)
	{
	if(Rslt != eBSFerrNotFasta)
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile: Unable to open '%s' [%s] %s",pszFile,Fasta.ErrText((teBSFrsltCodes)Rslt),Fasta.GetErrMsg());
	return(Rslt);
	}

AllocdBuffSize = (size_t)cMaxAXAllocBuffChunk * 16;
// note malloc is used as can then simply realloc to expand as may later be required
if((pSeqBuff = (uint8_t *)malloc(AllocdBuffSize)) == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile:- Unable to allocate memory (%u bytes) for sequence buffer",(uint32_t)AllocdBuffSize);
	Fasta.Close();
	return(eBSFerrMem);
	}
AvailBuffSize = AllocdBuffSize;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"ProcessFastaFile:- Adding %s..",pszFile);

bFirstEntry = true;
bEntryCreated = false;
SeqID = 0;
BuffOfs = 0;
NumSeqsUnderlength = 0;
while((Rslt = SeqLen = Fasta.ReadSequence(&pSeqBuff[BuffOfs],(int)min(AvailBuffSize,(size_t)cMaxAXAllocBuffChunk),true,false)) > eBSFSuccess)
	{
	if(SeqLen == eBSFFastaDescr)		// just read a descriptor line
		{
		SeqID++;
		if(bEntryCreated)				// add any previous entry
			{
			if(BuffOfs < (size_t)MinSeqLen)
				NumSeqsUnderlength += 1;
			else
				if((Rslt=m_pSfxFile->AddEntry(szName,pSeqBuff,(uint32_t)BuffOfs)) < eBSFSuccess)
					{
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile - error %d %s",Rslt,m_pSfxFile->GetErrMsg());
					break;
					}
			Rslt = eBSFSuccess;
			}
		Descrlen = Fasta.ReadDescriptor(szDescription,cBSFDescriptionSize);
		// An assumption - will one day bite real hard - is that the
		// fasta descriptor line starts with some form of unique identifier.
		// Use this identifier as the entry name.
		if(sscanf(szDescription," %s[ ,]",szName)!=1)
			sprintf(szName,"%s.%d",pszFile,++SeqID);

		bFirstEntry = false;
		bEntryCreated = true;
		BuffOfs = 0;
		continue;
		}
	else
		if(bFirstEntry)	// if there was no descriptor then dummy up one...
			{
			SeqID++;
			sprintf(szName,"%s.%d",pszFile,SeqID);
			strcpy(szDescription,"No Description provided");
			bFirstEntry = false;
			bEntryCreated = true;
			}

	// remove any repeat masking flags so that sorts can actually sort
	// if run of more than 25 Ns and at least 5 Ns to end of buffer then randomly mutate
	// every 13th N
	//	e.g <25Ns>r<12Ns>r<12Ns> where r is a pseudorandom base
	pMskBase = &pSeqBuff[BuffOfs];
	int SeqNs = 0;
	for(MskIdx = 0; MskIdx < SeqLen; MskIdx++,pMskBase++)
		{
		*pMskBase &= ~cRptMskFlg;
		if(*pMskBase == eBaseN && (MskIdx+5) < SeqLen)
			{
			if(++SeqNs > 25 &&
				pMskBase[1] == eBaseN &&
				pMskBase[2] == eBaseN &&
				pMskBase[3] == eBaseN &&
				pMskBase[4] == eBaseN)
				{
				if(!(SeqNs % 13))	// mutate every 13th
					*pMskBase = rand() % 4;
				}
			}
		else
			SeqNs = 0;
		}

	BuffOfs += SeqLen;
	AvailBuffSize -= SeqLen;
	if(AvailBuffSize < (size_t)(cMaxAXAllocBuffChunk / 8))
		{
		size_t NewSize = (size_t)cMaxAXAllocBuffChunk + AllocdBuffSize;
		uint8_t *pTmp;
		if((pTmp = (uint8_t *)realloc(pSeqBuff,NewSize))==nullptr)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile:- Unable to reallocate memory (%u bytes) for sequence buffer",(uint32_t)NewSize);
			return(eBSFerrMem);
			}
		pSeqBuff = pTmp;
		AllocdBuffSize = NewSize;
		AvailBuffSize = AllocdBuffSize - BuffOfs;
		}
	}

if(Rslt >= eBSFSuccess && bEntryCreated && BuffOfs > 0)			// close entry
	{
	if(BuffOfs < (size_t)MinSeqLen)
		{
		NumSeqsUnderlength += 1;
		Rslt = eBSFSuccess;
		}
	else
		{
		if((Rslt=m_pSfxFile->AddEntry(szName,pSeqBuff,(uint32_t)BuffOfs)) < eBSFSuccess)
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile - error %d %s",Rslt,m_pSfxFile->GetErrMsg());
		else
			Rslt = eBSFSuccess;
		}
	}
if(pSeqBuff != nullptr)
	free(pSeqBuff);
if(NumSeqsUnderlength > 0)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"ProcessFastaFile - %u sequences not accepted for indexing as length under %dbp ",NumSeqsUnderlength,MinSeqLen);
return(Rslt);
}


// simulates the loading of fasta sequences totalling SimGenomeSize Gbp
// the sequences will be randomly (near random!) generated into cNumSimChroms chromosomes of equal length
int
ProcessSimGenome(int SimGenomeSize)		// simulate genome of this total sequence size split into cNumSimChroms chromosomes of random length
{
int Rslt;
char szName[80];
int CurChromID;
uint8_t *pSeqBuff;
etSeqBase *pBase;
size_t BuffOfs;
size_t BuffRandRelOfs;
size_t AllocdBuffSize;
size_t AvailBuffSize;
size_t SimChromSize;
size_t SimGenomeSizebp;
size_t RemainingSimGenomeSizebp;
int NumSimChroms;

TRandomCombined<CRandomMother,CRandomMersenne> RGseeds((int)SimGenomeSize); // looking for reproducible pseudorandom numbers so seed with targeted genome size

SimGenomeSizebp = (size_t)SimGenomeSize * cGbpSeqSize;
NumSimChroms = cNumSimChroms;

AllocdBuffSize = SimGenomeSizebp / (size_t)(NumSimChroms - 1); // a little larger than actually required
// note malloc is used as can then simply realloc to expand as may later be required
if((pSeqBuff = (uint8_t *)malloc(AllocdBuffSize)) == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessSimGenome:- Unable to allocate memory (%u bytes) for simulated sequence buffer",(uint32_t)AllocdBuffSize);
	return(eBSFerrMem);
	}
AvailBuffSize = AllocdBuffSize;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"ProcessSimGenome:- Simulating and adding %d chromosome sequences ..",NumSimChroms);

Rslt = eBSFSuccess;
RemainingSimGenomeSizebp = SimGenomeSizebp;
for(CurChromID = 1; CurChromID <= NumSimChroms; CurChromID++)
	{
	SimChromSize = RemainingSimGenomeSizebp / (NumSimChroms - (CurChromID - 1));
	RemainingSimGenomeSizebp -= SimChromSize;
	// the first chromosome contains all randomly generated bases
	// subsequent chromosomes contain random variations on the first chromosome so as to ensure there are repeated regions
	if(CurChromID == 1)
		{
		pBase = (etSeqBase *)pSeqBuff;
		for(BuffOfs = 0; BuffOfs < SimChromSize; BuffOfs++,pBase++)
			{
			*pBase = (etSeqBase) RGseeds.IRandom(0,3);
			}
		}
	else
		{
		BuffRandRelOfs = (size_t)RGseeds.IRandom(1,1000);
		pBase = (etSeqBase *)pSeqBuff;
		for(BuffOfs = 0; BuffOfs < SimChromSize; BuffOfs+=BuffRandRelOfs,pBase+=BuffRandRelOfs)
			{
			*pBase = (etSeqBase) RGseeds.IRandom(0,3);
			BuffRandRelOfs = (size_t)RGseeds.IRandom(10,1000);
			}
		}
	BuffOfs = SimChromSize;

	sprintf(szName,"SimChrom%d",CurChromID);
	if((Rslt=m_pSfxFile->AddEntry(szName,pSeqBuff,(uint32_t)BuffOfs)) < eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessSimGenome - adding %s of size %zd (cumulative total bases %zd) error %d %s",szName,BuffOfs,(int64_t)SimGenomeSizebp - RemainingSimGenomeSizebp,Rslt,m_pSfxFile->GetErrMsg());
		break;
		}
	if(!(CurChromID % 1000))
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Added chromosome: %s ..",szName);
	}
if(Rslt >= eBSFSuccess)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Completed adding %d simulated chromosome sequences totaling %dGbp",NumSimChroms,SimGenomeSize);

if(pSeqBuff != nullptr)
	free(pSeqBuff);
return(Rslt);
}

size_t
FileSeqLen(char* pszFastaFile)			// input file containing fasta sequences
{
	int Rslt;
	int SeqLen;
	size_t TotSeqLen;
	bool bInSeq;
	int NumAccepted;
	int NumProcessed;

	CFasta* pFasta = nullptr;
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Starting to process source fasta file '%s'", pszFastaFile);

	if ((pFasta = new CFasta()) == nullptr)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to create CFasta object");
		return((size_t)eBSFerrObj);
		}

	if ((Rslt = pFasta->Open(pszFastaFile, true)) < eBSFSuccess)
		{
		while (pFasta->NumErrMsgs())
			gDiagnostics.DiagOut(eDLFatal, gszProcName, pFasta->GetErrMsg());
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to open fasta file '%s'", pszFastaFile);
		delete pFasta;
		return((size_t)Rslt);
		}

	TotSeqLen = 0;
	NumAccepted = 0;
	NumProcessed = 0;
	char szFeature[128];
	char szBEDFeature[256];
	int DescrLen;
	bInSeq = false;
	while ((Rslt = SeqLen = pFasta->ReadSequence(nullptr, 0)) > eBSFSuccess)
		{
		if (SeqLen == eBSFFastaDescr)		// just read a descriptor line, parse out the feature identifier
			{
			DescrLen = pFasta->ReadDescriptor(szBEDFeature, sizeof(szBEDFeature));
			sscanf(szBEDFeature, " %s[ ,]", szFeature);
			szFeature[35] = '\0';
			bInSeq = false;
			continue;
			}
		bInSeq = true;
		NumProcessed += 1;
		TotSeqLen += SeqLen;
		continue;
		}

	delete pFasta;
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Sequences processed %d ", NumProcessed);
	return(TotSeqLen);
}

int
CreateBioseqSuffixFile(int Mode,
						int MaxKMers,			// downstream this suffix array will only be used for processing of these sized maximal K-mers, only advantage is that shorter k-mers are quicker to sort than longer k-mers
						int MinSeqLen,			// only accept for indexing sequences which are at least this length
						int SimGenomeSize,		// if 1..1000 then simulating indexing of a genome of this size in Gbp.
					   int MaxThreads,			// max threads
					   bool bSOLiD,				// true if to process for colorspace (SOLiD)
						int NumInputFiles,			// number of input file specs
						char *pszInputFiles[],		// names of input files (wildcards allowed)
						char *pszDestSfxFile,	// output suffix array to this file
						char *pszRefSpecies,	// species assembly which is being indexed
						char *pszDescr,			// describes assembly
						char *pszTitle)        // title to use
{
int Rslt;
int Idx;
int NumGlobs;
int64_t SumFileSizes;
uint32_t DupEntries[cMaxDupEntries];
int NumDupEntries;
char szDupEntry[100];

InitpSfx();

CSimpleGlob glob(SG_GLOB_FULLSORT);

if(Mode == 2)
	SumFileSizes = (int64_t)SimGenomeSize * cGbpSeqSize;
else
	{
	// determine crude estimate of total genome size
	for(Idx = 0; Idx < NumInputFiles; Idx++)
		if (glob.Add(pszInputFiles[Idx]) < SG_SUCCESS)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to glob '%s",pszInputFiles[Idx]);
			Reset();
			return(eBSFerrOpnFile);
			}


	SumFileSizes = 0;
	for (NumGlobs = 0; NumGlobs < glob.FileCount(); NumGlobs += 1)
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Will process in this order: %d '%s", NumGlobs+1,glob.File(NumGlobs));

		SumFileSizes += FileSeqLen(glob.File(NumGlobs));
		}

	if(NumGlobs == 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to locate any input files");
		Reset();
		return(eBSFerrOpnFile);
		}
	}

m_pSfxFile = new CSfxArray;
if(m_pSfxFile == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"CreateBioseqSuffixFile - error Unable to create instance of CSfxArray");
	Reset();
	return(eBSFerrMem);
	}

m_pSfxFile->SetMaxQSortThreads(MaxThreads);

if(Mode == 2 && (pszDestSfxFile == nullptr || pszDestSfxFile[0]=='\0'))
	Rslt=m_pSfxFile->Open(false,bSOLiD);
else
	Rslt=m_pSfxFile->Open(pszDestSfxFile,true,Mode==1 ? true : false,bSOLiD);

if(Rslt !=eBSFSuccess)
	{
	while(m_pSfxFile->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pSfxFile->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"CreateBioseqSuffixFile, unable to create '%s' - error %d %s",pszDestSfxFile,Rslt,m_pSfxFile->GetErrMsg());
	Reset();
	return(Rslt);
	}

if((Rslt=m_pSfxFile->SetDescription(pszDescr)) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"CreateBioseqSuffixFile: Unable to set description '%s' into %s",pszDescr,pszDestSfxFile);
	Reset();
	return(Rslt);
	}
if((Rslt=m_pSfxFile->SetTitle(pszTitle)) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"CreateBioseqSuffixFile: Unable to set title '%s' into %s",pszTitle,pszDestSfxFile);
	Reset();
	return(Rslt);
	}

if((Rslt = m_pSfxFile->SetDatasetName(pszRefSpecies)) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"CreateBioseqSuffixFile: Unable to set dataset %s",pszRefSpecies);
	Reset();
	return(Rslt);
	}
m_pSfxFile->SetInitalSfxAllocEls(SumFileSizes);	// just a hint which is used for initial allocations by suffix processing
m_pSfxFile->SetMaxBaseCmpLen(MaxKMers);			// suffix array will only be used for processing of these maximally sized K-mers
Rslt = eBSFSuccess;

if(Mode == 2)			// Mode == 2 if simulating the genome sequence and indexing same
	{
	Rslt = ProcessSimGenome(SimGenomeSize);
	if(Rslt < eBSFSuccess)
		{
		m_pSfxFile->Close(false);
		Reset();
		return(Rslt);
		}
	}
else
	{
	for (int n = 0; Rslt >= eBSFSuccess &&  n < glob.FileCount(); ++n)
		{
		// try opening as a fasta file
		Rslt = ProcessFastaFile(MinSeqLen,glob.File(n));
		if(Rslt < eBSFSuccess)
			{
			m_pSfxFile->Close(false);
			Reset();
			return(Rslt);
			}

		// check for duplicate entry names
		if((NumDupEntries = m_pSfxFile->ChkDupEntries(cMaxDupEntries,&DupEntries[0])) > 0)
			{
			while(NumDupEntries--)
				{
				m_pSfxFile->GetIdentName(DupEntries[NumDupEntries],sizeof(szDupEntry)-1,szDupEntry); // get sequence name for specified entry identifier
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"CreateBioseqSuffixFile, duplicate sequence entry name '%s' in file '%s'",szDupEntry,glob.File(n));
				}
			m_pSfxFile->Close(false);
			Reset();
			return(-1);
			}
		}
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"CreateBioseqSuffixFile: sorting suffix array...");
if((Rslt = m_pSfxFile->Finalise()) < eBSFSuccess)
	Rslt = m_pSfxFile->Close();
gDiagnostics.DiagOut(eDLInfo,gszProcName,"CreateBioseqSuffixFile: completed...");
Reset();
return(Rslt);
}

