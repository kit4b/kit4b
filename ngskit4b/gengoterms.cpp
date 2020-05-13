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

const int cMaxObjectIDLen = 50;			// max expected length of any object identifier
const int cMaxObjectSymLen = 50;		// max expected length of any object symbol
const int cMaxQualLen = 25;				// qualifier max expected length
const int cMaxAspectLen = 5;			// max expected length of any aspect (should be 1 but who knows!)
const int cMaxGOIDLen = 50;				// max expected length of any GO term ID


int
Process(char *pszOBOFile,				// GO ontology file to parse
		char *pszGOTermFile);			// biogoterm file to generate

#ifdef _WIN32
int gengoterms(int argc, char* argv[])
{
	// determine my process name
	_splitpath(argv[0], NULL, NULL, gszProcName, NULL);
#else
int
gengoterms(int argc, char** argv)
{
	// determine my process name
	CUtility::splitpath((char*)argv[0], NULL, gszProcName);
#endif
int iScreenLogLevel;		// level of screen diagnostics
int iFileLogLevel;			// level of file diagnostics
char szLogFile[_MAX_PATH];	// write diagnostics to this file

int Rslt;
char szOBOFile[_MAX_PATH];		// source GO ontology file
char szGOTermFile[_MAX_PATH];	// generated biogoterm file

struct arg_lit  *help    = arg_lit0("h","help",                 "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_file *OBOFile = arg_file1("i",NULL,"<file>",			"source OBO ontology file");
struct arg_file *GOTermFile = arg_file1("o",NULL,"<file>",		"output to biogoterm file");

struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					OBOFile,GOTermFile,
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
	if(FileLogLevel->count && !LogFile->count)
		{
		printf("\nError: FileLogLevel '-f%d' specified but no logfile '-F<logfile>\n'",FileLogLevel->ival[0]);
		exit(1);
		}

	iScreenLogLevel = iFileLogLevel = FileLogLevel->count ? FileLogLevel->ival[0] : eDLInfo;
	if(iFileLogLevel < eDLNone || iFileLogLevel > eDLDebug)
		{
		printf("\nError: FileLogLevel '-l%d' specified outside of range %d..%d\n",iFileLogLevel,eDLNone,eDLDebug);
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

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %s",kit4bversion);

	strcpy(szOBOFile,OBOFile->filename[0]);
	strcpy(szGOTermFile,GOTermFile->filename[0]);

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing parameters:");

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"input OBO ontology file: '%s'",szOBOFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"output biogoterm file: '%s'",szGOTermFile);
	gStopWatch.Start();
#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	Rslt = Process(szOBOFile,szGOTermFile);
	Rslt = Rslt >=0 ? 0 : 1;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	exit(Rslt);
	}
else
	{
	printf("\n%s Generate GO and PO terms, Version %s\n",gszProcName,kit4bversion);
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nUse '-h' to view option and parameter usage\n");
	exit(1);
	}
return 0;
}

int
Process(char *pszOBOFile,		// GO ontology file to parse
		char *pszGOTermFile)	// biogoterm file to generate
{
int Rslt;

FILE *pOBOStream = NULL;
CGOTerms *pGOTerms = NULL;

if(pszOBOFile == NULL || *pszOBOFile == '\0' ||
	pszGOTermFile == NULL || *pszGOTermFile == '\0')
	return(eBSFerrParams);

if((pOBOStream = fopen(pszOBOFile,"r"))==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open GO ontology file %s error: %s",pszOBOFile,strerror(errno));
	return(eBSFerrOpnFile);
	}

pGOTerms = new CGOTerms();
Rslt = pGOTerms->Open(pszGOTermFile,true);
if(Rslt == eBSFSuccess)
	Rslt = pGOTerms->Parse(pOBOStream);
if(Rslt == eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Total GO Terms - Cellular: %d, Biological: %d, Molecular: %d, Plant Anatomical: %d, Plant Developmental: %d",
		pGOTerms->NumGOTerms(eONTCellular),pGOTerms->NumGOTerms(eONTBiological),pGOTerms->NumGOTerms(eONTMolecular),pGOTerms->NumGOTerms(eONTPlantAnatomical),pGOTerms->NumGOTerms(eONTPlantDev));
	Rslt = pGOTerms->Close(true);
	}
else
	{
	while(pGOTerms->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pGOTerms->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to process '%s'",pszGOTermFile);
	}
fclose(pOBOStream);
delete pGOTerms;
return(Rslt);
}
