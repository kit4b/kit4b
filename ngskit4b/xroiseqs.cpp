// xroiseqs.cpp : Defines the entry point for the console application.
// Extract ROI (as specified in BED file) sequences from a suffix assembly file and output these into a multifasta file
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
#include "./ngskit4b.h"

const int cAllocSeqLen = 1000000;	 // allocate for sequences in this sized chunks

typedef enum eSamplingProcMode {
	eSMPdefault = 0				// currently default processing only
} etSamplingProcMode;

int ParseRegions(char *pszRegionList);
char *Regions2Txt(int Regions);
int TrimQuotes(char *pszTxt);

int Process(int PMode,						// currently just the single default which is the complete sequence for each ROI
			char *pszROIFile,		// input BED file containing regions to be sampled
			char *pszSamplesFile,	// output ROI sequences to this multifasta file
			char *pszSfxFile);	// source ROI sequences from this suffix array file


#ifdef _WIN32
int xroiseqs(int argc, char* argv[])
{
// determine my process name
_splitpath(argv[0],NULL,NULL,gszProcName,NULL);
#else
int 
xroiseqs(int argc, char **argv)
{
// determine my process name
CUtility::splitpath((char *)argv[0],NULL,gszProcName);
#endif
int iScreenLogLevel;		// level of screen diagnostics
int iFileLogLevel;			// level of file diagnostics
char szLogFile[_MAX_PATH];	// write diagnostics to this file

int Rslt;

int PMode;						// currently just the single default which is the complete sequence for each ROI
char szROIFile[_MAX_PATH];		// input BED file containing regions to be sampled
char szSamplesFile[_MAX_PATH];	// output ROI sequences to this multifasta file
char szSfxFile[_MAX_PATH];	// source ROI sequences from this fasta assembly file

// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_int *ScreenLogLevel=arg_int0("S", "ScreenLogLevel",	"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int  *procmode = arg_int0("m","procmode","<int>",	"processing mode 0:default");
struct arg_file* roifile = arg_file1("i", "roifile", "<file>",			"input BED file containing ROI specifying regions containing fasta sequences to be extracted");
struct arg_file* samplesfile = arg_file1("o", "samplesfile","<file>",	"output ROI sequences to this multifasta file");
struct arg_file *sfxfile = arg_file1("I","sfxfile","<file>",			"source ROI sequences from this suffix array assembly file");

struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,ScreenLogLevel,LogFile,
					procmode,roifile,samplesfile,sfxfile,
					end};

char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

    /* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
        {
		printf("Usage: %s ", gszProcName);
        arg_print_syntax(stdout,argtable,"\n");
        arg_print_glossary(stdout,argtable,"  %-25s %s\n");
		exit(1);
        }

    /* special case: '--version' takes precedence error reporting */
if (version->count > 0)
        {
		printf("\n%s: Version: %s\n",gszProcName, kit4bversion);
		exit(1);
        }

if (!argerrors)
	{
	iScreenLogLevel = ScreenLogLevel->count ? ScreenLogLevel->ival[0] : eDLInfo;
	if(iScreenLogLevel < eDLNone || iScreenLogLevel > eDLDebug)
		{
		printf("\nError: ScreenLogLevel '-S%d' specified outside of range %d..%d",iScreenLogLevel,eDLNone,eDLDebug);
		exit(1);
		}
	if(FileLogLevel->count && !LogFile->count)
		{
		printf("\nError: FileLogLevel '-f%d' specified but no logfile '-F<logfile>'",FileLogLevel->ival[0]);
		exit(1);
		}

	iFileLogLevel = FileLogLevel->count ? FileLogLevel->ival[0] : eDLInfo;
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

	PMode = procmode->count ? procmode->ival[0] : eSMPdefault;
	if(PMode < eSMPdefault || PMode > eSMPdefault)
		PMode = eSMPdefault;

	strncpy(szROIFile,roifile->filename[0],sizeof(szROIFile)-1);
	CUtility::TrimQuotedWhitespcExtd (szROIFile);
	if(szROIFile[0] == '\0')
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: After removal of whitespace, no input ROI file specified with '-i<filespec>' parameter)\n");
		exit(1);
		}

	strncpy(szSamplesFile,samplesfile->filename[0],sizeof(szSamplesFile)-1);
	CUtility::TrimQuotedWhitespcExtd (szSamplesFile);
	if(szSamplesFile[0] == '\0')
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: After removal of whitespace, no output ROI sequences file specified with '-o<filespec>' parameter)\n");
		exit(1);
		}

	strncpy(szSfxFile,sfxfile->filename[0],sizeof(szSfxFile)-1);
	CUtility::TrimQuotedWhitespcExtd (szSfxFile);
	if(szSfxFile[0] == '\0')
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: After removal of whitespace, no input suffix assembly file specified with '-I<filespec>' option)\n");
		exit(1);
		}


		// now that command parameters have been parsed then initialise diagnostics log system
	if(!gDiagnostics.Open(szLogFile,(etDiagLevel)iScreenLogLevel,(etDiagLevel)iFileLogLevel,true))
		{
		printf("\nError: Unable to start diagnostics subsystem.");
		if(szLogFile[0] != '\0')
			printf(" Most likely cause is that logfile '%s' can't be opened/created",szLogFile);
		exit(1);
		}

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %s Processing parameters:", kit4bversion);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode: Default");
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input suffix assembly file: '%s'",szSfxFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output ROI fasta sequences to this file: '%s'",szSamplesFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input ROI BED file: '%s'",szROIFile);
#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = Process(PMode,						// currently just the single default which is the complete sequence for each ROI
					szROIFile,		// input BED file containing regions to be sampled
					szSamplesFile,	// output ROI sequences to this multifasta file
					szSfxFile);	// source ROI sequences from this suffix assembly file);
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


// ParseRegions
// Parses space or comma delimited list of regions in which
// 1 == Intergenic, 2 == US, 3 == 5'UTR, 4 == CDS, 5 == Intron, 6 == 3'UTR, 7 == DS
//
// Returns bitmap of regions or -1 if parse errors
// If no region specified then assumes all regions are selected
int
ParseRegions(char *pszRegionList)
{
// parse out region list
char Chr;
int Region = 0;
if(pszRegionList == NULL || *pszRegionList == '\0')
	return(cFeatBitIG | cFeatBitUpstream | cFeatBit5UTR | cFeatBitCDS | cFeatBitIntrons | cFeatBit3UTR | cFeatBitDnstream);

while(Chr = *pszRegionList++) {
	if(isspace(Chr) || Chr == ',')		// accept spaces and commas as separators
		continue;
	switch(Chr) {
		case '1':						// intergenic to be filtered
			Region |= cFeatBitIG;
			break;
		case '2':						// 5'US to be filtered
			Region |= cFeatBitUpstream;
			break;
		case '3':						// 5'UTR to be filtered
			Region |= cFeatBit5UTR;
			break;
		case '4':
			Region |= cFeatBitCDS;		// CDS to be filtered
			break;
		case '5':
			Region |=  cFeatBitIntrons;	// any intronic to be filtered
			break;
		case '6':
			Region |=  cFeatBit3UTR;	// any 3'UTR to be filtered
			break;
		case '7':
			Region |=  cFeatBitDnstream;	// any 3'DS to be filtered 	
			break;
		default:
			return(-1);
		}
	}
return(Region);
}

// Regions2Txt
// Returns textual representation of regions
char *
Regions2Txt(int Regions)
{
static char szRegions[200];
if(!Regions)
	return((char *)"All");
if(Regions & cFeatBitIG || Regions == 0)
	strcpy(szRegions,"Intergenic");
else
	szRegions[0] = '\0';
if(Regions & cFeatBitUpstream)
	{
	if(szRegions[0] != '\0')
		strcat(szRegions,", ");
	strcat(szRegions,"5'US");
	}
if(Regions & cFeatBit5UTR)
	{
	if(szRegions[0] != '\0')
		strcat(szRegions,", ");
	strcat(szRegions,"5'UTR");
	}
if(Regions & cFeatBitCDS)
	{
	if(szRegions[0] != '\0')
		strcat(szRegions,", ");
	strcat(szRegions,"CDS");
	}
if(Regions & cFeatBitIntrons)
	{
	if(szRegions[0] != '\0')
		strcat(szRegions,", ");
	strcat(szRegions,"Introns");
	}
if(Regions & cFeatBit3UTR)
	{
	if(szRegions[0] != '\0')
		strcat(szRegions,", ");
	strcat(szRegions,"3'UTR");
	}
if(Regions & cFeatBitDnstream)
	{
	if(szRegions[0] != '\0')
		strcat(szRegions,", ");
	strcat(szRegions,"3'DS");
	}
return(szRegions);
}

// TrimQuotes
// Removes any leading and trailing quotes from pszTxt
int
TrimQuotes(char *pszTxt)
{
char *pDst = pszTxt;
char Chr;
int Len = 0;
while((Chr = *pszTxt++))	
	{
	if((!Len || *pszTxt == '\0') && (Chr == '"' || Chr == '\''))
		continue;
	*pDst++ = Chr;
	Len++;
	}
*pDst = Chr;
return(Len);
}




int
Process(int PMode,				// currently just the single default which is the complete sequence for each ROI
		char *pszROIFile,		// input BED file containing regions to be sampled
		char *pszSamplesFile,	// output ROI sequences to this multifasta file
		char *pszSfxFile)	// source ROI sequences from this suffix array file
{
int Rslt;
CBEDfile *pBED;
CSfxArray *pSfxArray;

char szCurChrom[cMaxDatasetSpeciesChrom];
int AllocSeqLen;
etSeqBase *pSeqBuffer;
int hRsltsFile;

pSeqBuffer = NULL;
AllocSeqLen = 0;


if((pBED = new CBEDfile())==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"\nUnable to instantiate CBEDfile");
	return(eBSFerrObj);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading BED/GFF3 file %s",pszROIFile);
if((Rslt = pBED->Open(pszROIFile,eBTAnyBed))!=eBSFSuccess)
	{
	while(pBED->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pBED->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open BED/GFF3 file %s",pszROIFile);
	delete pBED;
	return(Rslt);
	}

if((pSfxArray = new CSfxArray())==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"\nUnable to instantiate CSfxArray");
	return(eBSFerrObj);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading suffix assembly file %s",pszSfxFile);
if((Rslt = pSfxArray->Open(pszSfxFile))!=eBSFSuccess)
	{
	while(pSfxArray->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pSfxArray->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open suffix assembly file '%s'",pSfxArray);
	delete pBED;
	delete pSfxArray;
	return(Rslt);
	}

#ifdef _WIN32
if((hRsltsFile = open(pszSamplesFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((hRsltsFile = open(pszSamplesFile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create %s - %s",pszSamplesFile,strerror(errno));
		delete pBED;
		delete pSfxArray;
		return(eBSFerrCreateFile);
		}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Output ROI fasta file created: '%s'",pszSamplesFile);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Iterating over ROIs sequences ...");
// all input files now opened
if((pSeqBuffer = new uint8_t[cAllocSeqLen]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory for buffering input fasta sequences");
	delete pBED;
	delete pSfxArray;
	close(hRsltsFile);
	return(eBSFerrMem);
	}
AllocSeqLen = cAllocSeqLen;

// iterate over all the BED features

char szFeatureName[cMaxDatasetSpeciesChrom];
char szChrom[cMaxDatasetSpeciesChrom];
int SfxEntryID;
int FeatLen;
int FeatStart;
int FeatEnd;
int CurFeatureID;
int NxtFeatureID;
int NumFeatures;
char *pszSeqBuffer;
int SeqLen;
int AllocszSeqLen;
int SeqBuffOfs;
int NumWarnings;

if((pszSeqBuffer = new char[cAllocSeqLen]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory for buffering output ROI sequences");
	delete pBED;
	delete pSfxArray;
	close(hRsltsFile);
	delete []pSeqBuffer;
	return(eBSFerrMem);
	}
AllocszSeqLen = cAllocSeqLen;
SeqBuffOfs = 0;
NumWarnings = 0;
Rslt = eBSFSuccess;
szCurChrom[0] = '\0';
CurFeatureID = 0;
SfxEntryID = 0;
NumFeatures = 0;
while ((NxtFeatureID = pBED->GetNextFeatureID(CurFeatureID)) > 0)
	{
	CurFeatureID = NxtFeatureID;
	if((Rslt = pBED->GetFeature(CurFeatureID,szFeatureName,szChrom,&FeatStart,&FeatEnd)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to locate expected feature ROI in input BED");
		break;
		}
	FeatLen = FeatEnd - FeatStart + 1;
	if(pSeqBuffer == NULL || FeatLen > AllocSeqLen)
		{
		if(pSeqBuffer != NULL)
			delete []pSeqBuffer;
		AllocSeqLen = 0;
		if((pSeqBuffer = new uint8_t [FeatLen + cAllocSeqLen])==NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory for buffering");
			Rslt = eBSFerrMem;
			break;
			}
		AllocSeqLen = FeatLen + cAllocSeqLen;
		}
	if(SfxEntryID < 1 || szCurChrom[0] == '\0' || stricmp(szCurChrom, szChrom))
		{
		if((SfxEntryID = pSfxArray->GetIdent(szChrom)) < 1)
			{
			SfxEntryID = 0;
			if(stricmp(szCurChrom, szChrom) && NumWarnings++ < 100)
				gDiagnostics.DiagOut(eDLWarn,gszProcName,"Unable to locate ROI chromosome '%s' in suffix assembly file '%s'",szChrom,pszSfxFile);
			strcpy(szCurChrom,szChrom);
			continue;
			}
		}
	
	if((SeqLen=pSfxArray->GetSeq(SfxEntryID,FeatStart-1,pSeqBuffer,FeatLen))!=FeatLen)
		{
		szCurChrom[0] = '\0';
		if(NumWarnings++ < 100)
			gDiagnostics.DiagOut(eDLWarn,gszProcName,"Incomplete ROI '%s' sequence from chromosome '%s' loaded from suffix assembly file '%s'",szFeatureName,szChrom,pszSfxFile);
		continue;
		}

	if(NumFeatures++ == 0)
		SeqBuffOfs += sprintf(&pszSeqBuffer[SeqBuffOfs], ">%s %s:%d:%d", szFeatureName, szChrom, FeatStart, SeqLen);
	else
		SeqBuffOfs += sprintf(&pszSeqBuffer[SeqBuffOfs], "\n>%s %s:%d:%d", szFeatureName, szChrom, FeatStart, SeqLen);

	int NumCols;
	int SeqOfs = 0;
	int NxtFastaCol = 0;
	while(SeqLen)
		{
		NumCols = SeqLen > 70 ? 70 : SeqLen;
		if(NumCols > 70)
			NumCols = 70;
		SeqBuffOfs += sprintf(&pszSeqBuffer[SeqBuffOfs],"\n");
		CSeqTrans::MapSeq2Ascii(&pSeqBuffer[SeqOfs],NumCols,&pszSeqBuffer[SeqBuffOfs]);
		SeqBuffOfs += NumCols;
		SeqLen -= NumCols;
		SeqOfs += NumCols;
		if((SeqBuffOfs + 1000) > AllocszSeqLen)
			{
			CUtility::RetryWrites(hRsltsFile,pszSeqBuffer,SeqBuffOfs);
			SeqBuffOfs = 0;
			}
		}
	}

if (Rslt == eBSFSuccess && SeqBuffOfs)
	CUtility::RetryWrites(hRsltsFile,pszSeqBuffer,SeqBuffOfs);

if(hRsltsFile != -1)
	{
#ifdef _WIN32
	_commit(hRsltsFile);
#else
	fsync(hRsltsFile);
#endif
	close(hRsltsFile);
	}
if(pSeqBuffer != NULL)
	delete []pSeqBuffer;
if(pszSeqBuffer != NULL)
	delete []pszSeqBuffer;
if(pBED != NULL)
	delete pBED;
if(pSfxArray != NULL)
	delete pSfxArray;
return(Rslt);
}

