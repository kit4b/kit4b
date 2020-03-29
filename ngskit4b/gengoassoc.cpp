/*
This toolkit is a source base clone of 'BioKanga' release 4.4.2 (https://github.com/csiro-crop-informatics/biokanga) and contains
significant source code changes enabling new functionality and resulting process parameterisation changes. These changes have resulted in
incompatibilty with 'BioKanga'.

Because of the potentential for confusion by users unaware of functionality and process parameterisation changes then the modified source base
and resultant compiled executables have been renamed to 'kit4b' - K-mer Informed Toolkit for Bioinformatics.
The renaming will force users of the 'BioKanga' toolkit to examine scripting which is dependent on existing 'BioKanga'
parameterisations so as to make appropriate changes if wishing to utilise 'kit4b' parameterisations and functionality.

'kit4b' is being released under the Opensource Software License Agreement (GPLv3)
'kit4b' is Copyright (c) 2019, 2020
Please contact Dr Stuart Stephen < stuartjs@g3web.com > if you have any questions regarding 'kit4b'.

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

#include "./ngskit4b.h"
int
Process(etGOAssocParseType Type,		// expected association file type
		char *pszGOAssocFile,			// source GO association file to parse
		char *pszBIOGOAssocFile,		// biogoassoc file to generate
		char *pszNameMapFile,			// optional gene name mapping file
		char *pszFiltFile); 			// optional gene name filter file


#ifdef _WIN32
int gengoassoc(int argc, char* argv[])
{
	// determine my process name
	_splitpath(argv[0], NULL, NULL, gszProcName, NULL);
#else
int
gengoassoc(int argc, char** argv)
{
	// determine my process name
	CUtility::splitpath((char*)argv[0], NULL, gszProcName);
#endif
int iScreenLogLevel;		// level of screen diagnostics
int iFileLogLevel;			// level of file diagnostics
char szLogFile[_MAX_PATH];	// write diagnostics to this file

int Rslt;
etGOAssocParseType iType;			// expected association file type
char szGOAssocFile[_MAX_PATH];		// source GO association file
char szBIOGOAssocFile[_MAX_PATH];	// generated biogoassoc file
char szNameMapFile[_MAX_PATH];		// optional gene name mapping file
char szFiltFile[_MAX_PATH];			// optional gene name filter file


// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_int *ScreenLogLevel=arg_int0("S", "ScreenLogLevel",	"<int>","Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *Type=arg_int0("t", "type",	"<int>",			"association file type 0:UCSC (default), 1:GO standard, 2:GO TAIR, 3:GO Flybase, 4:PGSB GFF3");
struct arg_file *GOAssocFile = arg_file1("i",NULL,"<file>",	    "association file, gene to GO terms ");
struct arg_file *BIOGOAssocFile = arg_file1("o",NULL,"<file>",  "output to biogoassoc file");
struct arg_file *NameMapFile = arg_file0("I",NULL,"<file>",     "optional, gene mapping file");
struct arg_file *FiltFile = arg_file0("r",NULL,"<file>",		"optional, gene filtering file");

struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,ScreenLogLevel,LogFile,
					Type,GOAssocFile,BIOGOAssocFile,NameMapFile,FiltFile,
					end};

char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

/* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
        {
		printf("\n%s %s %s, Version %s\nOptions ---\n", gszProcName, gpszSubProcess->pszName, gpszSubProcess->pszFullDescr, kit4bversion);
        arg_print_syntax(stdout,argtable,"\n");
        arg_print_glossary(stdout,argtable,"  %-25s %s\n");
		printf("\nNote: Parameters can be entered into a parameter file, one parameter per line.");
		printf("\n      To invoke this parameter file then precede its name with '@'");
		printf("\n      e.g. %s @myparams.txt\n",gszProcName);
		printf("\nPlease report any issues regarding usage of %s to stuart.stephen@csiro.au\n\n",gszProcName);
		exit(1);
        }
    /* special case: '--version' takes precedence error reporting */
if (version->count > 0)
        {
		printf("\n%s Version %s\n",gszProcName,kit4bversion);
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

	iType = (etGOAssocParseType)(Type->count ? Type->ival[0] :eGOAPMUCSU);
	if(iType < eGOAPMUCSU || iType > eGOAPGFF)
		{
		printf("\nError: Unsupported GO association file type '-x%d' specified",(int)iType);
		exit(1);
		}
	strcpy(szGOAssocFile,GOAssocFile->filename[0]);
	strcpy(szBIOGOAssocFile,BIOGOAssocFile->filename[0]);

	if(NameMapFile->count)
		strcpy(szNameMapFile,NameMapFile->filename[0]);
	else
		szNameMapFile[0] = '\0';

	if(FiltFile->count)
		strcpy(szFiltFile,FiltFile->filename[0]);
	else
		szFiltFile[0] = '\0';

			// now that command parameters have been parsed then initialise diagnostics log system
	if(!gDiagnostics.Open(szLogFile,(etDiagLevel)iScreenLogLevel,(etDiagLevel)iFileLogLevel,true))
		{
		printf("\nError: Unable to start diagnostics subsystem.");
		if(szLogFile[0] != '\0')
			printf(" Most likely cause is that logfile '%s' can't be opened/created",szLogFile);
		exit(1);
		}

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %s",kit4bversion);

	char *pTxt;
	switch(iType) {
		case eGOAPMUCSU:
			pTxt = (char *)"parse file as UCSC";
			break;
		case eGOAPGO:
			pTxt = (char *)"parse file as GO annotations";
			break;
		case eGOAPGOTAIR:
			pTxt = (char *)"parse file as TAIR GO annotations (reduce isoforms e.g. strip off any suffixed '.[1-9]')";
			break;
		case  eGOAPGOFB:
			pTxt = (char *)"parse file as GO Flybase annotations";
			break;
		case eGOAPGFF:
			pTxt = (char*)"parse file as PGSB GFF3 with note GO associations";
			break;
		}

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"expected file type: '%s'",pTxt);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"input gene to GO association file: '%s'",szGOAssocFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"output biogoassoc file: '%s'",szBIOGOAssocFile);
	if(szNameMapFile[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"gene name mapping file: '%s'",szNameMapFile);
	if(szFiltFile[0] != '\0')
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"gene name filter file: '%s'",szFiltFile);

	gStopWatch.Start();
#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	Rslt = Process(iType,szGOAssocFile,szBIOGOAssocFile,szNameMapFile,szFiltFile);
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
Process(etGOAssocParseType Type,	// expected GO association type
		char *pszGOAssocFile,		// source GO association file to parse
		char *pszBIOGOAssocFile,	// biogoassoc file to generate
		char *pszNameMapFile,		// optional gene name mapping file
		char *pszFiltFile)			// optional gene name filter file
{
int Rslt;

CGOAssocs *pGOAssoc = NULL;

if(pszGOAssocFile == NULL || *pszGOAssocFile == '\0' ||
	pszBIOGOAssocFile == NULL || *pszBIOGOAssocFile == '\0')
	return(eBSFerrParams);

pGOAssoc = new CGOAssocs();
if((Rslt = pGOAssoc->Open(pszBIOGOAssocFile,true))!=eBSFSuccess)
	{
	while(pGOAssoc->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pGOAssoc->GetErrMsg());
	delete pGOAssoc;
	return(Rslt);
	}

if((Rslt = pGOAssoc->Parse(Type,pszGOAssocFile,pszNameMapFile,pszFiltFile)) < eBSFSuccess)
	{
	while(pGOAssoc->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pGOAssoc->GetErrMsg());
	delete pGOAssoc;
	return(Rslt);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"NumGeneMaps: %d NumGeneGOitems: %d GOGenesCnt: %d",
						pGOAssoc->GetNumGeneMaps(),
						pGOAssoc->GetNumGeneGOitems(),
						pGOAssoc->GetNumGenes());

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Number of genes associated to BIO: %d CEL: %d MOL: %d",
						pGOAssoc->GetNumGOGenes(eONTBiological),
						pGOAssoc->GetNumGOGenes(eONTCellular),
						pGOAssoc->GetNumGOGenes(eONTMolecular));

Rslt = pGOAssoc->Close(true);
delete pGOAssoc;
return(Rslt);
}


