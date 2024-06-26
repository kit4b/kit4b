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
// Extracts fasta sequences from a multifasta file
// Sequences to be extracted are identified by their descriptors

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

CUtility RegExpr;                      // used for regular expression processing    

const int cMaxExtractDescrs = 20;					// at most this number of extraction descriptors can be processed

int
fastaextractproc(int PMode,				// processing mode - 0 by matching descriptor, 1 subsample by selecting every Nth sequence, 2 extracting sequences for each feature defined in this BED
		int NthSel,						// if > 0 then select every Nth descriptor.sequence to be extracted
		int XSense,						// extraction sense - 0 - original sense, 1 - reverse only, 2 - complement only 3 - reverse complement
		int NumDescrs,					// number of descriptors
		char *pszDescrs[],				// extract sequences with these descriptors
		char* pszBED,					// in PMode 2, extracting sequences for each feature defined in this BED
		char *pszInFasta,				// extract from this multifasta file
		char *pszOutFasta);				// output extracted sequences to this multifasta file



#ifdef _WIN32
int fastaextract(int argc, char* argv[])
{
// determine my process name
_splitpath(argv[0],nullptr,nullptr,gszProcName,nullptr);
#else
int 
fastaextract(int argc, char** argv)
{
// determine my process name
CUtility::splitpath((char *)argv[0],nullptr,gszProcName);
#endif

int iScreenLogLevel;		// level of screen diagnostics
int iFileLogLevel;			// level of file diagnostics
char szLogFile[_MAX_PATH];	// write diagnostics to this file
int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure

int PMode;					// processing mode
int NthSel;					// select every Nth sequence
int XSense;					// extraction sense - 0 - original sense, 1 - reverse only, 2 - complement only 3 - reverse complement

char szInBEDFile[_MAX_PATH];	// in PMode 2, extracting sequences for each feature defined in this BED
char szInFile[_MAX_PATH];		// read from this file
char szOutFile[_MAX_PATH];		// write extracted sequences to this file

int NumExtractDescrs;
char *pszExtractDescrs[cMaxExtractDescrs+1];

char szSQLiteDatabase[_MAX_PATH];	// results summaries to this SQLite file
char szExperimentName[cMaxDatasetSpeciesChrom+1];			// experiment name
char szExperimentDescr[1000];		// describes experiment

// command line args
struct arg_lit  *help    = arg_lit0("h","help",                 "Print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"Print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"Diagnostics log file");

struct arg_int *mode = arg_int0("m","mode","<int>",				"Processing mode: 0 extract by descriptor, 1 extract every Nth sequence, 2 extract subsequences defined in BED file starting at chrom.loci");
struct arg_int *nthsel = arg_int0("n","nth","<int>",			"Extract every Nth descriptor 1 .. Nth");
struct arg_str *extractdescrs = arg_strn("e","extract","<str>",0,cMaxExtractDescrs,		"Extract sequences with these descriptors (wildcards allowed)");
struct arg_file* inputbedfile = arg_file0("b", "bed", "<file>",	"Extract sequences for each feature defined in this Input BED file");
struct arg_file *inputfile = arg_file1("i","in","<file>",		"Input file containing sequences to extract");
struct arg_file *outfile = arg_file1("o","out","<file>",		"Output extracted sequences to this file");

struct arg_int *xsense = arg_int0("x","xsense","<int>",           "Extracted sequences as: 0 - original sense, 1 - reverse only, 2 - complement only 3 - reverse complement");

struct arg_file *summrslts = arg_file0("q","sumrslts","<file>",					"Output results summary to this SQLite3 database file");
struct arg_str *experimentname = arg_str0("w","experimentname","<str>",			"experiment name SQLite3 database file");
struct arg_str *experimentdescr = arg_str0("W","experimentdescr","<str>",		"experiment description SQLite3 database file");

struct arg_end *end = arg_end(200);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					summrslts,experimentname,experimentdescr,
	                mode,nthsel,extractdescrs,  xsense, inputbedfile, inputfile, outfile,
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
		printf("\nError: FileLogLevel '-f%d' specified but no logfile '-F<logfile>\n'",FileLogLevel->ival[0]);
		return(1);
		}

	iScreenLogLevel = iFileLogLevel = FileLogLevel->count ? FileLogLevel->ival[0] : eDLInfo;
	if(iFileLogLevel < eDLNone || iFileLogLevel > eDLDebug)
		{
		printf("\nError: FileLogLevel '-l%d' specified outside of range %d..%d\n",iFileLogLevel,eDLNone,eDLDebug);
		return(1);
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
		gProcessID = gSQLiteSummaries.AddProcess((char *)gszProcName,(char *)gszProcName,(char *)szExperimentDescr);
		if(gProcessID < 1)
			return(1);
		gProcessingID = gSQLiteSummaries.StartProcessing(gExperimentID,gProcessID,(char *)kit4bversion);
		if(gProcessingID < 1)
			return(1);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Initialised SQLite database '%s' for extraction results summary",szSQLiteDatabase);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"SQLite database experiment identifier for '%s' is %d",szExperimentName,gExperimentID);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"SQLite database process identifier for '%s' is %d",(char *)gszProcName,gProcessID);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"SQLite database processing instance identifier is %d",gProcessingID);
		}
	else
		{
		szSQLiteDatabase[0] = '\0';
		szExperimentDescr[0] = '\0';
		}

	PMode = mode->count ? mode->ival[0] : 0;
	if(PMode < 0 || PMode > 2)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Processing mode '-m%d' must be in range 0..2",PMode);
		return(1);
		}

	
	if(PMode == 1)
		{
		NthSel = nthsel->count ? nthsel->ival[0] : 0;
		if(NthSel <= 0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Sample every Nth sequence requested but '-n<th>' not specified as %d",NthSel);
			return(1);
			}
		}
	else
		NthSel = 0;

	if(PMode == 0 && extractdescrs->count)
		{
		int Idx;
		int LenExtractDescr;
		NumExtractDescrs = 0;
		for(Idx=0;Idx < extractdescrs->count; Idx++)
			{
			pszExtractDescrs[Idx] = nullptr;
			LenExtractDescr = (int)strlen(extractdescrs->sval[Idx]);
			if(pszExtractDescrs[NumExtractDescrs] == nullptr)
				pszExtractDescrs[Idx] = new char [LenExtractDescr+1];
			strcpy(pszExtractDescrs[NumExtractDescrs],extractdescrs->sval[Idx]);
			CUtility::TrimQuotes(pszExtractDescrs[NumExtractDescrs]);
			if(pszExtractDescrs[NumExtractDescrs][0] != '\0')
				NumExtractDescrs++;
			}

		if(!NumExtractDescrs)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: After removal of whitespace, descriptors specified with '-e<descriptor>' option)\n");
			exit(1);
			}
		}
	else
		{
		NumExtractDescrs = 0;
		pszExtractDescrs[0] = nullptr;
		}

	szInBEDFile[0] = '\0';
	if (PMode == 2)
		{
		if(inputbedfile->count > 1)
			{
			strcpy(szInBEDFile, inputbedfile->filename[0]);
			CUtility::TrimQuotedWhitespcExtd(szInBEDFile);
			}
		if (szInBEDFile[0] == '\0')
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: No input BED file specified containing feature specifications using '-b<file>' option)\n");
			exit(1);
			}
		}

	XSense = xsense->count ? xsense->ival[0] : 0;
	if(XSense < 0 || XSense > 3)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Extraction sense '-x%d' must be in range 0..3",XSense);
		return(1);
		}
	
	strcpy(szInFile,inputfile->filename[0]);
	CUtility::TrimQuotedWhitespcExtd(szInFile);

	strcpy(szOutFile,outfile->filename[0]);
	CUtility::TrimQuotedWhitespcExtd(szOutFile);

// show user current resource limits
#ifndef _WIN32
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Resources: %s",CUtility::ReportResourceLimits());
#endif

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing parameters:");
	const char *pszDescr;
	switch(PMode) {
		case 0:
			pszDescr = "Extract fasta sequence(s) matching descriptors";
			break;
		case 1:
			pszDescr = "Sample every Nth sequence";
			break;
		case 2:
			pszDescr = "Extract fasta subsequence(s) as specified in BED formated file";
			break;
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode is : '%s'",pszDescr);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Experiment name : '%s'",szExperimentName);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Experiment description : '%s'",szExperimentDescr);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"SQLite database file: '%s'",szSQLiteDatabase);

	switch(PMode) {
		case 0:
			for(int Idx = 0; Idx < NumExtractDescrs; Idx++)
				gDiagnostics.DiagOutMsgOnly(eDLInfo,"Descriptor (%d) : '%s'",Idx+1,pszExtractDescrs[Idx]);
			break;
		case 1:
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Sample every Nth sequence: %d",NthSel);
			break;
		case 2:
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "BED formated file soecifying features: '%s'", szInBEDFile);
			break;
		}

	switch(XSense) {
			case 0:
				pszDescr = "original";
				break; 
			case 1:
				pszDescr = "reverse";
				break; 
			case 2:
				pszDescr = "complement";
				break; 
			case 3:
				pszDescr = "reverse complement";
				break;
			}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Sequences extracted will be in : %s sense",pszDescr);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Extract sequences from file : '%s'",szInFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Write extracted features to file: '%s'",szOutFile);

	if(gExperimentID > 0)
		{
		int ParamID;
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTText,(int)strlen(szLogFile),"log",szLogFile);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,sizeof(PMode),"mode",&PMode);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,sizeof(XSense),"xsense",&XSense);

	    ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,sizeof(NumExtractDescrs),"NumExtractDescrs",&NumExtractDescrs);
		switch(PMode) {
			case 0:
				for(int Idx=0; Idx < NumExtractDescrs; Idx++)
					ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTText,(int)strlen(pszExtractDescrs[Idx]),"extractdescrs",pszExtractDescrs[Idx]);
				break;
			case 1:
				ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTInt32,sizeof(NthSel),"nthsel",&NthSel);
				break;
			case 2:
				ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID, ePTText, (int)strlen(szInFile), "bed", szInBEDFile);
				break;
			}

		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTText,(int)strlen(szOutFile),"out",szOutFile);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTText,(int)strlen(szInFile),"in",szInFile);

		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTText,(int)strlen(szSQLiteDatabase),"sumrslts",szSQLiteDatabase);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTText,(int)strlen(szExperimentName),"experimentname",szExperimentName);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID,ePTText,(int)strlen(szExperimentDescr),"experimentdescr",szExperimentDescr);
		}


#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();

	Rslt = fastaextractproc(PMode,NthSel,XSense,NumExtractDescrs,pszExtractDescrs, szInBEDFile,szInFile,szOutFile);
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


return 0;
}

const size_t cMaxAllocBuffChunk = 0x0ffffff;		// allocate for input sequences buffering in these sized chunks

CBEDfile *m_pBedFile; // Extract sequences for each feature defined in this Input BED file
int m_hFEInFile;	// file handle for opened multifasta input file
int m_hFEOutFile;	// file handle for opened multifasta output file

size_t m_SeqBuffLen;		// number of bases currently buffered in m_pSeqBuff
size_t m_AllocdSeqBuffMem;  // size of memory currently allocated to m_pSeqBuff
uint8_t *m_pSeqBuff;			// buffers sequences as read from file

int m_NumExtractDescrs;		// number of extract descriptors
regex* m_pRegexDescrs[100]; // to hold precompiled regular expressions

void
FEReset(void)
{
if(m_pBedFile != nullptr)
	{
	delete m_pBedFile;
	m_pBedFile = nullptr;
	}

if(m_hFEInFile != -1)
	{
	close(m_hFEInFile);
	m_hFEInFile = -1;
	}
if(m_hFEOutFile != -1)
	{
#ifdef _WIN32
	_commit(m_hFEOutFile);
#else
	fsync(m_hFEOutFile);
#endif
	close(m_hFEOutFile);
	m_hFEOutFile = -1;
	}

if(m_pSeqBuff != nullptr)
	{
#ifdef _WIN32
	free(m_pSeqBuff);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pSeqBuff != MAP_FAILED)
		munmap(m_pSeqBuff,m_AllocdSeqBuffMem);
#endif
	m_pSeqBuff = nullptr;
	}

m_AllocdSeqBuffMem = 0;
m_SeqBuffLen = 0;
for (int Idx = 0; Idx < m_NumExtractDescrs; Idx++)
	delete m_pRegexDescrs[Idx];
m_NumExtractDescrs = 0;
}

void
FEInit(void)
{
m_hFEInFile = -1;
m_hFEOutFile = -1;
m_pBedFile = nullptr;
m_pSeqBuff = nullptr;
m_pRegexDescrs[0] = nullptr;
m_NumExtractDescrs = 0;
FEReset();
}

etSeqBase *
AllocSeqBuff(size_t SeqLen)				// allocate for at least this sequence length
{
size_t memreq;
etSeqBase *pTmp;

if(m_pSeqBuff != nullptr && m_AllocdSeqBuffMem >= SeqLen)
	return(m_pSeqBuff);

if(m_pSeqBuff == nullptr)
	{
	memreq = max(SeqLen,(size_t)cMaxAllocBuffChunk);
#ifdef _WIN32
	m_pSeqBuff = (etSeqBase *) malloc(SeqLen);	// initial and perhaps the only allocation
	if(m_pSeqBuff == nullptr)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AllocSeqBuff: Memory allocation of %zd bytes - %s",(int64_t)SeqLen,strerror(errno));
		return(nullptr);
		}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	m_pSeqBuff = (etSeqBase *)mmap(nullptr,memreq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pSeqBuff == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AllocSeqBuff: Memory allocation of %zd bytes through mmap()  failed - %s",(int64_t)memreq,strerror(errno));
		m_pSeqBuff = nullptr;
		return(nullptr);
		}
#endif
	}
else
	{
	memreq = SeqLen + cMaxAllocBuffChunk;
#ifdef _WIN32
	pTmp = (etSeqBase *) realloc(m_pSeqBuff,memreq);
#else
	pTmp = (etSeqBase *)mremap(m_pSeqBuff,m_AllocdSeqBuffMem,memreq,MREMAP_MAYMOVE);
	if(pTmp == MAP_FAILED)
		pTmp = nullptr;
#endif
	if(pTmp == nullptr)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AllocSeqBuff: Memory re-allocation to %zd bytes - %s",(int64_t)memreq,strerror(errno));
		return(nullptr);
		}
	m_pSeqBuff = pTmp;
	}
m_AllocdSeqBuffMem = memreq;
return(m_pSeqBuff);
}

int
WriteSeqFile(int XSense,						// extraction sense - 0 - original sense, 1 - reverse only, 2 - complement only 3 - reverse complement
			char *pszDescr,int SeqLen,etSeqBase *pSeq)
{
char Chr;
int LineLen;
int BuffIdx;
int SeqIdx;
uint8_t *pszOutBuff;

pszOutBuff = new uint8_t [0x07ffff];

switch(XSense) {
	case 0:		// original sense
		break;
	case 1:     // reverse only
		CSeqTrans::ReverseSeq(SeqLen,pSeq);
		break;
	case 2:		// complement only
		CSeqTrans::ComplementStrand(SeqLen,pSeq);
		break;
	case 3:		// reverse and complement
		CSeqTrans::ReverseComplement(SeqLen,pSeq);
		break;
	}

LineLen = 0;
		// write out descriptor to file
BuffIdx = sprintf((char *)pszOutBuff,">%s\n",pszDescr);
for(SeqIdx = 0; SeqIdx < SeqLen; SeqIdx++,pSeq++,LineLen++)
	{
	if(LineLen > 75)
		{
		pszOutBuff[BuffIdx++] = '\n';
		LineLen = 0;
		}
	switch(*pSeq & 0x07) {
		case eBaseA:
			Chr = 'A';
			break;
		case eBaseC:
			Chr = 'C';
			break;
		case eBaseG:
			Chr = 'G';
			break;
		case eBaseT:
			Chr = 'T';
			break;
		default:
			Chr = 'N';
			break;
		}
	pszOutBuff[BuffIdx++] = Chr;
	if(BuffIdx > (0x07ffff - 0x3ff))
		{
		CUtility::RetryWrites(m_hFEOutFile,pszOutBuff,BuffIdx);
		BuffIdx = 0;
		}
	}
pszOutBuff[BuffIdx++] = '\n';
CUtility::RetryWrites(m_hFEOutFile,pszOutBuff,BuffIdx);
delete []pszOutBuff;
return(SeqLen);
}

int
fastaextractproc(int PMode,				// processing mode - 0 by matching descriptor, 1 subsample by selecting every Nth sequence, 2 extracting sequences for each feature defined in this BED
		int NthSel,						// if > 0 then select every Nth descriptor.sequence to be extracted
		int XSense,						// extraction sense - 0 - original sense, 1 - reverse only, 2 - complement only 3 - reverse complement
		int NumDescrs,					// number of descriptor reqular expressions (in pszDescrs[]) to match on
		char *pszDescrs[],				// extract sequences with these descriptors
		char *pszBED,					// in PMode 2, extracting sequences for each feature defined in this BED
		char *pszInFasta,				// extract from this multifasta file
        char *pszOutFasta)				// output extracted sequences to this multifasta file
{
int Rslt;
int NumMatches;
int SeqID;
int Descrlen;
int SeqLen;
size_t AvailBuffSize;
bool bExtract;
char szInDescription[cBSFDescriptionSize];
char *pszDescriptor;
int CurChromID;
char *pszChromName;
char szChromName[cMaxDatasetSpeciesChrom];

CFasta Fasta;

FEInit();

if(PMode == 0 && NumDescrs >= 1)
	{
	if ((Rslt = RegExpr.CompileIncludeREs(NumDescrs, pszDescrs)) < eBSFSuccess)
		return(Rslt);
	}
else if (PMode == 2)
	{
	if ((m_pBedFile = new CBEDfile) == nullptr)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to instantiate CBEDfile");
		return(eBSFerrObj);
		}

	if ((Rslt = m_pBedFile->Open(pszBED, eBTAnyBed)) != eBSFSuccess)
		{
		while (m_pBedFile->NumErrMsgs())
			gDiagnostics.DiagOut(eDLFatal, gszProcName, m_pBedFile->GetErrMsg());
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to open '%s' for processing", pszBED);
		return(eBSFerrOpnFile);
		}
	}

#ifdef _WIN32
if((m_hFEOutFile = open(pszOutFasta, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((m_hFEOutFile = open(pszOutFasta,O_RDWR | O_CREAT |O_TRUNC, S_IREAD | S_IWRITE))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create or truncate %s - %s",pszOutFasta,strerror(errno));
	FEReset();
	return(eBSFerrCreateFile);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Output results file created/truncated: '%s'",pszOutFasta);

if((Rslt=Fasta.Open(pszInFasta,true))!=eBSFSuccess)
	{
	if(Rslt != eBSFerrNotFasta)
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile: Unable to open '%s' [%s] %s",pszInFasta,Fasta.ErrText((teBSFrsltCodes)Rslt),Fasta.GetErrMsg());
	return(Rslt);
	}

if(m_pSeqBuff == nullptr)				// if not already allocated then allocate to hold cMaxAllocBuffChunk bases 
	{
	SeqLen = cMaxAllocBuffChunk;
	if(AllocSeqBuff(SeqLen) == nullptr)
		{
		Rslt = eBSFerrMem;
		Fasta.Close();
		return(Rslt);
		}
	}

NumMatches = 0;
SeqID = 0;
m_SeqBuffLen = 0;
AvailBuffSize = m_AllocdSeqBuffMem;
bExtract = false;
while((Rslt = SeqLen = Fasta.ReadSequence(&m_pSeqBuff[m_SeqBuffLen],(int)min(AvailBuffSize,(size_t)cMaxAllocBuffChunk),true,false)) > eBSFSuccess)
	{
	if(SeqLen == eBSFFastaDescr)		// just read a descriptor line
		{
		if(bExtract && m_SeqBuffLen)
			{
			if(PMode != 2)
				WriteSeqFile(XSense,szInDescription,(int)m_SeqBuffLen,m_pSeqBuff);
			else
				{
				char szFeatChrom[cMaxDatasetSpeciesChrom + 1];
				char szFeatName[cMaxGeneNameLen + 1];
				int FeatStart;
				int FeatEnd;
				int NthFeature = 0;
				int CurFeatID = 0;
				while ((CurFeatID = m_pBedFile->LocateFeatureIDinRangeOnChrom(CurChromID, 0, (int)m_SeqBuffLen, ++NthFeature)) > 0)
					{
					Rslt = m_pBedFile->GetFeature(CurFeatID, szFeatName, szFeatChrom, &FeatStart, &FeatEnd);
					WriteSeqFile(XSense, szFeatName, FeatEnd-FeatStart+1, &m_pSeqBuff[FeatStart]);
					NumMatches++;
					}
				}
			}
		else
			if(bExtract)
				gDiagnostics.DiagOut(eDLWarn,gszProcName,"Matching descriptor '%s' but no sequence in file - '%s'",szInDescription,pszOutFasta);
		bExtract = false;
		m_SeqBuffLen = 0;
		AvailBuffSize = m_AllocdSeqBuffMem;
		bExtract = false;
		SeqID++;

		Descrlen = Fasta.ReadDescriptor(szInDescription,cBSFDescriptionSize-10);

		switch(PMode) {
			case 0:		// processing for sequences with matching descriptors
				if (NumDescrs > 0 && !RegExpr.MatchIncludeRegExpr(szInDescription))
					continue;
				break;

			case 1:		// sampling every Nth sequence, includes 1st sequence
				if (NthSel > 1 && (SeqID % NthSel))
					continue;
				break;

			case 2:		// processing for chrom or sequence identifer as specified in a BED file
				// chrom name expected to be prefix in descriptor and terminated by first whitespace
				char Chr;
				pszChromName = szChromName;
				pszDescriptor = szInDescription;
				for (int Idx = 0; Idx < sizeof(szChromName)-1; Idx++)
					{
					Chr = *pszDescriptor++;
					if(isspace(Chr))
						break;
					*pszChromName++ = Chr;
					}
				*pszChromName = '\0';
				if((CurChromID = m_pBedFile->LocateChromIDbyName(szChromName)) < 1)
					continue;
				break;
			}

		if(PMode != 2)
			NumMatches += 1;
		bExtract = true;				// flag that the sequence is to be extracted
		m_SeqBuffLen = 0;
		AvailBuffSize = m_AllocdSeqBuffMem;
		continue;
		}

	if(!bExtract)						// slough this sequence?
		continue;

	// accumulate complete sequence prior to writing to file, this enables subsequence extraction as needed for the BED processing
	m_SeqBuffLen += SeqLen;
	AvailBuffSize -= SeqLen;
	if(AvailBuffSize < (size_t)(cMaxAllocBuffChunk / 8))
		{
		if(AllocSeqBuff(m_AllocdSeqBuffMem + SeqLen) == nullptr)
			{
			Rslt = eBSFerrMem;
			Fasta.Close();
			return(Rslt);
			}
		AvailBuffSize = m_AllocdSeqBuffMem - m_SeqBuffLen;
		}
	}
if(Rslt < eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Encountered parse errors");
	while(Fasta.NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,Fasta.GetErrMsg());
	}
if(bExtract && m_SeqBuffLen)					// was previous sequence to be extracted?
	{
	if (PMode != 2)
		WriteSeqFile(XSense, szInDescription, (int)m_SeqBuffLen, m_pSeqBuff);
	else
		{
		char szFeatChrom[cMaxDatasetSpeciesChrom + 1];
		char szFeatName[cMaxGeneNameLen + 1];
		int FeatStart;
		int FeatEnd;
		int NthFeature = 0;
		int CurFeatID = 0;
		while ((CurFeatID = m_pBedFile->LocateFeatureIDinRangeOnChrom(CurChromID, 0, (int)m_SeqBuffLen, ++NthFeature)) > 0)
			{
			Rslt = m_pBedFile->GetFeature(CurFeatID, szFeatName, szFeatChrom, &FeatStart, &FeatEnd);
			WriteSeqFile(XSense, szFeatName, FeatEnd - FeatStart + 1, &m_pSeqBuff[FeatStart]);
			NumMatches++;
			}
		}
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Extracted %d sequences",NumMatches);
FEReset();
return(Rslt);
}