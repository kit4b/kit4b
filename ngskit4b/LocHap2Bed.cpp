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

#include "./ngskit4b.h"
#include "LocHap2Bed.h"

int Process (eModeLocHap Mode,				// processing mode
			 int MinCoverage,				// must be at least this coverage at SNP site
			 double MinAlleleProp,			// putative allele must be at least this proportion of total site read coverage
			 double PValueThres,			// only accept SNP alleles which have a PValue <= this threshold
			 char* pszTrackName,			// track name
			 char* pszAssemblyName,			// UCSC assembly name - for SARS-CoV-2 it is "wuhCor1"
			 char* pszExperimentDescr,		// describes experiment
			 int NumInWhitelist,			// number of white listed file names
			 char** ppszWhitelisted,		// names of white listed files
			 int NumInBlackList,			// number of black listed files
			 char** ppsBlacklisted,			// names of black listed files
			 int NumSNPFiles,				// number of input SNP files
			 char** ppszSNPFiles,			// input SNP files
			 char* pszOutFile);				// output SNPs to this UCSC Personal Genome SNP format file

#ifdef _WIN32
int LocHap2Bed (int argc, char *argv[])
{
	// determine my process name
	_splitpath (argv[0], NULL, NULL, gszProcName, NULL);
#else
int
LocHap2Bed (int argc, char **argv)
{
	// determine my process name
	CUtility::splitpath ((char *)argv[0], NULL, gszProcName);
#endif
	int Len;
	int iFileLogLevel;			// level of file diagnostics
	int iScreenLogLevel;		// level of file diagnostics
	char szLogFile[_MAX_PATH];	// write diagnostics to this file
	int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure
	int Idx;
	int NumberOfProcessors;		// number of installed CPUs

	eModeLocHap PMode;			// processing mode
	int MinCoverage;				// must be at least this coverage at SNP site
	double MinAlleleProp;			// putative allele must be at least this proportion of total site read coverage
	double PValueThres;				// only accept SNP alleles which have a PValue <= this threshold

	int NumSNPFiles;			// number of input SNP files
	char* ppszSNPFiles[cLHMaxSetMembers + 1];  // input SNP files
	char szOutBedFile[_MAX_PATH];		// output in UCSC BED format to this file

	int NumInWhitelist;				// number of white listed file names
	char * ppszWhitelisted[cLHMaxSetMembers];	// names of white listed files
	int NumInBlackList;					// number of black listed files
	char * ppsBlacklisted[cLHMaxSetMembers];	// names of black listed files

	char szTrackName[cMaxDatasetSpeciesChrom + 1];	// UCSC track name
	char szAssemblyName[cMaxDatasetSpeciesChrom + 1];	// UCSC assembly name
	char szExperimentDescr[cMaxDatasetSpeciesChrom + 1];	// describes experiment


	struct arg_lit *help = arg_lit0 ("h", "help", "print this help and exit");
	struct arg_lit *version = arg_lit0 ("v", "version,ver", "print version information and exit");
	struct arg_int *FileLogLevel = arg_int0 ("f", "FileLogLevel", "<int>", "Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
	struct arg_file *LogFile = arg_file0 ("F", "log", "<file>", "diagnostics log file");

	struct arg_int *pmode = arg_int0 ("m", "mode", "<int>", "processing mode: 0 default and currently only processing mode externally exposed");
	struct arg_file* snpfiles = arg_filen("i", "insnps", "<file>", 1, cLHMaxSetMembers, "Load SNPs from disnp or trisnp file(s) (can contain regular expressions)");

	struct arg_str *whitelisted = arg_strn ("w", "whitelist", "<whitelist>", 0, cLHMaxSetMembers, "names of files which are whitelisted (can contain regular expressions)");
	struct arg_str *blacklisted = arg_strn ("b", "blacklist", "<blacklist>", 0, cLHMaxSetMembers, "names of files which are blacklisted (can contain regular expressions)");

	struct arg_file *outbed = arg_file1 ("o", "outbed", "<file>", "output in UCSC BED format to this file");
	struct arg_str *experimentdescr = arg_str1 ("e", "experiment", "<str>", "UCSC track experiment description");
	struct arg_str *assemblyname = arg_str1 ("a", "assembly", "<str>", "UCSC assembly name");
	struct arg_str *trackname = arg_str1 ("t", "track", "<str>", "UCSC track name");
	struct arg_dbl *pvaluethres = arg_dbl0 ("p", "pvalue", "<dbl>", "SNP maximum allele PValue threshold (default 0.05, range 0.01..0.25)");
	struct arg_dbl *minalleleprop = arg_dbl0 ("P", "minalleleprop", "<dbl>", "SNP minimum allele proportion of loci site coverage (default 0.1, range 0.01..0.95)");
	struct arg_int *mincoverage = arg_int0 ("c", "mincoverage", "<int>", "SNP minimum loci site coverage (default 20, range 5..10000)");

	struct arg_end *end = arg_end (200);

	void *argtable[] = { help,version,FileLogLevel,LogFile,
						pmode,pvaluethres,minalleleprop,mincoverage,snpfiles,whitelisted,blacklisted,outbed,experimentdescr,assemblyname,trackname,end };

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
		if (!gDiagnostics.Open (szLogFile, (etDiagLevel)iScreenLogLevel, (etDiagLevel)iFileLogLevel, true))
		{
			printf ("\nError: Unable to start diagnostics subsystem\n");
			if (szLogFile[0] != '\0')
				printf (" Most likely cause is that logfile '%s' can't be opened/created\n", szLogFile);
			exit (1);
		}

		gDiagnostics.DiagOut (eDLInfo, gszProcName, "Subprocess %s Version %s starting", gpszSubProcess->pszName, kit4bversion);
		gExperimentID = 0;
		gProcessID = 0;
		gProcessingID = 0;
		szExperimentDescr[0] = '\0';

		PMode = pmode->count ? (eModeLocHap)pmode->ival[0] : eMLHdefault;
		if (PMode < eMLHdefault || PMode >= eMLHPlaceHolder)
		{
			gDiagnostics.DiagOut (eDLFatal, gszProcName, "Error: Processing mode '-m%d' specified outside of range %d..%d\n", PMode, eMLHdefault, eMLHPlaceHolder-1);
			exit (1);
		}

		strncpy (szExperimentDescr, experimentdescr->sval[0], sizeof (szExperimentDescr) - 1);
		szExperimentDescr[sizeof (szExperimentDescr) - 1] = '\0';
		CUtility::TrimQuotedWhitespcExtd (szExperimentDescr);
		CUtility::ReduceWhitespace (szExperimentDescr);
		if ((Len = (int)strlen (szExperimentDescr)) < 3 || Len > 80)
		{
			gDiagnostics.DiagOut (eDLFatal, gszProcName, "Error: Expected UCSC experiment description length to be in range of 3..80\n");;
			exit (1);
		}

		strncpy (szAssemblyName, assemblyname->sval[0], sizeof (szAssemblyName) - 1);
		szAssemblyName[sizeof (szAssemblyName) - 1] = '\0';
		CUtility::TrimQuotedWhitespcExtd (szAssemblyName);
		CUtility::ReduceWhitespace (szAssemblyName);
		if ((Len = (int)strlen (szAssemblyName)) < 3 || Len > 50)
		{
			gDiagnostics.DiagOut (eDLFatal, gszProcName, "Error: Expected UCSC assembly name length to be in range of 3..50\n");
			exit (1);
		}

		strncpy (szTrackName, trackname->sval[0], sizeof (szTrackName) - 1);
		szTrackName[sizeof (szTrackName) - 1] = '\0';
		CUtility::TrimQuotedWhitespcExtd (szTrackName);
		CUtility::ReduceWhitespace (szTrackName);
		if ((Len = (int)strlen (szTrackName)) < 3 || Len > 50)
		{
			gDiagnostics.DiagOut (eDLFatal, gszProcName, "Error: Expected UCSC track name length to be in range of 3..50\n");
			exit (1);
		}

		PValueThres = pvaluethres->count ? pvaluethres->dval[0] : 0.05;
		if (PValueThres <= 0.0)
			PValueThres = 0.01;
		else
			if (PValueThres > 0.25)
				PValueThres = 0.25;

		MinAlleleProp = minalleleprop->count ? minalleleprop->dval[0] : 0.10;
		if (MinAlleleProp <= 0.0)
			MinAlleleProp = 0.01;
		else
			if (MinAlleleProp > 0.95)
				MinAlleleProp = 0.95;

		MinCoverage = mincoverage->count ? mincoverage->ival[0] : 20;
		if (MinCoverage <= 5)
			MinCoverage = 5;
		else
			if (MinCoverage > 10000)
				MinCoverage = 10000;

		// show user current resource limits
#ifndef _WIN32
		gDiagnostics.DiagOut (eDLInfo, gszProcName, "Resources: %s", CUtility::ReportResourceLimits ());
#endif

#ifdef _WIN32
		SYSTEM_INFO SystemInfo;
		GetSystemInfo (&SystemInfo);
		NumberOfProcessors = SystemInfo.dwNumberOfProcessors;
#else
		NumberOfProcessors = sysconf (_SC_NPROCESSORS_CONF);
#endif

		if(!snpfiles->count)
		{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: No input SNP file(s) specified with '-i<filespec>' option)");
			exit(1);
		}

		memset(ppszSNPFiles,0,sizeof(ppszSNPFiles));
		char szTmpFile[_MAX_PATH+1];
		CSimpleGlob glob(SG_GLOB_FULLSORT);
		for(NumSNPFiles = Idx = 0; NumSNPFiles < cLHMaxSetMembers && Idx < snpfiles->count; Idx++)
			{
			strncpy(szTmpFile, snpfiles->filename[Idx],sizeof(szTmpFile));
			szTmpFile[_MAX_PATH] = '\0';
			CUtility::TrimQuotedWhitespcExtd(szTmpFile);
			if(szTmpFile[0] == '\0')
				continue;
			glob.Init();
			if(glob.Add(szTmpFile) >= SG_SUCCESS)
				{
				for(int n = 0; Rslt >= eBSFSuccess && n < glob.FileCount(); ++n)
					{
					if(NumSNPFiles >= cLHMaxSetMembers)
						{
						gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Too many input SNP files specified with '-i<filespec>' option");
						exit(1);
						}
					if(ppszSNPFiles[NumSNPFiles] == NULL)
						ppszSNPFiles[NumSNPFiles] = new char[_MAX_PATH];
					strncpy(ppszSNPFiles[NumSNPFiles], glob.File(n), _MAX_PATH);
					ppszSNPFiles[NumSNPFiles++][_MAX_PATH - 1] = '\0';
					}
				}
			}

		if(!NumSNPFiles)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: After wildcard expansion andd/or removal of whitespace, no input SNP file(s) specified with '-i<filespec>' option");
			exit(1);
			}


		strcpy (szOutBedFile, outbed->filename[0]);
		CUtility::TrimQuotedWhitespcExtd (szOutBedFile);
		if (szOutBedFile[0] == '\0')
		{
			gDiagnostics.DiagOut (eDLFatal, gszProcName, "No UCSC BED output file specified");
			exit (1);
		}

		memset (ppszWhitelisted, 0, sizeof (ppszWhitelisted));
		memset(ppsBlacklisted, 0, sizeof(ppsBlacklisted));
		NumInWhitelist = 0;
		NumInBlackList = 0;
		if (whitelisted->count)
		{
			for (int Idx = 0; NumInWhitelist < cLHMaxSetMembers && Idx < whitelisted->count; Idx++)
			{
				ppszWhitelisted[Idx] = NULL;
				if (ppszWhitelisted[NumInWhitelist] == NULL)
					ppszWhitelisted[NumInWhitelist] = new char[_MAX_PATH];
				strncpy (ppszWhitelisted[NumInWhitelist], whitelisted->sval[Idx], _MAX_PATH);
				ppszWhitelisted[NumInWhitelist][_MAX_PATH - 1] = '\0';
				CUtility::TrimQuotedWhitespcExtd (ppszWhitelisted[NumInWhitelist]);
				if (ppszWhitelisted[NumInWhitelist][0] != '\0')
					NumInWhitelist++;
			}
		}

		if (blacklisted->count)
		{
			for (int Idx = 0; NumInBlackList < cLHMaxSetMembers && Idx < blacklisted->count; Idx++)
			{
				ppsBlacklisted[Idx] = NULL;
				if (ppsBlacklisted[NumInBlackList] == NULL)
					ppsBlacklisted[NumInBlackList] = new char[_MAX_PATH];
				strncpy (ppsBlacklisted[NumInBlackList], blacklisted->sval[Idx], _MAX_PATH);
				ppsBlacklisted[NumInBlackList][_MAX_PATH - 1] = '\0';
				CUtility::TrimQuotedWhitespcExtd (ppsBlacklisted[NumInBlackList]);
				if (ppsBlacklisted[NumInBlackList][0] != '\0')
					NumInBlackList++;
			}
		}


		gDiagnostics.DiagOut (eDLInfo, gszProcName, "Processing parameters:");
		const char *pszDescr;
		switch (PMode) {
			case eMLHdefault:
				pszDescr = "input SNPs are those generated by kalign";
				break;
		}



		gDiagnostics.DiagOutMsgOnly (eDLInfo, "disnp or trisnp localised haplotype to UCSC BED format conversion : '%s'", pszDescr);

		gDiagnostics.DiagOutMsgOnly (eDLInfo, "SNP minimum loci site coverage : '%d'", MinCoverage);
		gDiagnostics.DiagOutMsgOnly (eDLInfo, "SNP minimum allele proportion of loci site coverage : '%f'", MinAlleleProp);

		gDiagnostics.DiagOutMsgOnly (eDLInfo, "SNP allele PValue threshold : '%f'", PValueThres);
		gDiagnostics.DiagOutMsgOnly (eDLInfo, "UCSC assembly name : '%s'", szAssemblyName);
		gDiagnostics.DiagOutMsgOnly (eDLInfo, "UCSC track name : '%s'", szTrackName);
		gDiagnostics.DiagOutMsgOnly (eDLInfo, "UCSC track experiment description : '%s'", szExperimentDescr);

		for(Idx = 0; Idx < NumSNPFiles; Idx++)
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "Input SNP file (%d) : '%s'", Idx + 1, ppszSNPFiles[Idx]);
		// check file extension, if '.csv' then generate output formated for csv instead of the default BED format
		int Len;
		if ((Len = (int)strlen (szOutBedFile)) >= 5 && !stricmp (&szOutBedFile[Len - 4], ".csv"))
			gDiagnostics.DiagOutMsgOnly (eDLInfo, gszProcName, "Output to '%s' is in CSV format - NOT CURRENTLY SUPPORTED", szOutBedFile);
		else
			gDiagnostics.DiagOutMsgOnly (eDLInfo, gszProcName, "Output to '%s' is in UCSC BED format", szOutBedFile);

		gDiagnostics.DiagOutMsgOnly (eDLInfo, "Set operation : '%s'", pszDescr);

		if (NumInWhitelist)
			gDiagnostics.DiagOutMsgOnly (eDLInfo, "Number of whitelisted file names : '%d'", NumInWhitelist);

		if (NumInBlackList)
			gDiagnostics.DiagOutMsgOnly (eDLInfo, "Number of blacklisted file names : '%d'", NumInBlackList);


		if (szExperimentDescr[0] != '\0')
			gDiagnostics.DiagOutMsgOnly (eDLInfo, "Experiment description: %s", szExperimentDescr);

#ifdef _WIN32
		SetPriorityClass (GetCurrentProcess (), BELOW_NORMAL_PRIORITY_CLASS);
#endif
		gStopWatch.Start ();
		Rslt = 0;
		Rslt = Process (PMode,					// processing mode
						MinCoverage,				// must be at least this coverage at SNP site
						MinAlleleProp,				// putative allele must be at least this proportion of total site read coverage
						PValueThres,				// only accept SNP alleles which have a PValue <= this threshold
						szTrackName,				// track name
						szAssemblyName,				// UCSC assembly name - for SARS-CoV-2 it is "wuhCor1"
						szExperimentDescr,			// describes experiment
						NumInWhitelist,				// number of white listed file names
						ppszWhitelisted,			// names of white listed files
						NumInBlackList,				// number of black listed files
						ppsBlacklisted,				// names of black listed files
						NumSNPFiles,				// number of input SNP files
						ppszSNPFiles,				// input SNP files
						szOutBedFile);				// output in UCSC BED format to this file
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
	return 0;
}


int Process (eModeLocHap Mode,				// processing mode
			 int MinCoverage,				// must be at least this coverage at SNP site
			 double MinAlleleProp,			// putative allele must be at least this proportion of total site read coverage
			 double PValueThres,			// only accept SNP alleles which have a PValue <= this threshold
			 char* pszTrackName,			// track name
			 char* pszAssemblyName,		// UCSC assembly name - for SARS-CoV-2 it is "wuhCor1"
			 char* pszExperimentDescr,		// describes experiment
			 int NumInWhitelist,			// number of white listed file names
			 char** ppszWhitelisted,		// names of white listed files
			 int NumInBlackList,			// number of black listed files
			 char** ppsBlacklisted,		// names of black listed files
			 int NumSNPFiles,				// number of input SNP files
			 char** ppszSNPFiles,			// input SNP files
			 char* pszOutFile)				// output SNPs to this UCSC Personal Genome SNP format file
{
int Rslt;
CLocHap2Bed *pCLocHap2Bed;
if ((pCLocHap2Bed = new CLocHap2Bed) == NULL)
	{
	gDiagnostics.DiagOut (eDLFatal, gszProcName, "Unable to instantiate instance of CLocHap2Bed");
	return(eBSFerrInternal);
	}

Rslt = pCLocHap2Bed->Process (Mode, MinCoverage, MinAlleleProp, PValueThres, pszTrackName, pszAssemblyName, pszExperimentDescr, 
								  NumInWhitelist, ppszWhitelisted, NumInBlackList, ppsBlacklisted, NumSNPFiles, ppszSNPFiles, pszOutFile);
if (pCLocHap2Bed != NULL)
	delete pCLocHap2Bed;
return(Rslt);
}

CLocHap2Bed::CLocHap2Bed (void)
{
	m_pCSV = NULL;
	m_pszLineBuff = NULL;
	m_hOutFile = -1;
	m_ppszSNPFiles = NULL;
	Reset ();
}


CLocHap2Bed::~CLocHap2Bed (void)
{
if (m_pCSV != NULL)
	delete m_pCSV;
if (m_hOutFile != -1)
	close (m_hOutFile);
if (m_pszLineBuff != NULL)
	delete[]m_pszLineBuff;
}

void
CLocHap2Bed::Reset (void)
{
if (m_pCSV != NULL)
	{
	delete m_pCSV;
	m_pCSV = NULL;
	}
if (m_hOutFile != -1)
	{
	if (m_LineBuffOffs && m_pszLineBuff != NULL)
		CUtility::SafeWrite (m_hOutFile, m_pszLineBuff, m_LineBuffOffs);
#ifdef _WIN32
	_commit (m_hOutFile);
#else
	fsync (m_hOutFile);
#endif
	close (m_hOutFile);
	m_hOutFile = -1;
	}

if (m_pszLineBuff != NULL)
	{
	delete[]m_pszLineBuff;
	m_pszLineBuff = NULL;
	}
m_AllocdLineBuff = 0;
m_LineBuffOffs = 0;
m_NumSNPFiles = 0;
m_szOutBEDFile[0] = '\0';
}

int 
CLocHap2Bed::Process2BED(void)						// processing to BED format
{
int Rslt;
int NumFields;
int NumElsParsed;
bool bBiSNPs;
int ExpNumSNPs;
int FileIdx;
int FldIdx;
int SNPIdx;
int BaseIdx;
char BaseChr;
int* pCnt;
int AlleleCombinationIdx;
int NumAlleleCombinations;
m_pCSV->SetMaxFields(cCSVLHMaxFields);
char *pszSNPFile;
char *pTxt;
tsLocalisedHaplotype LocalisedHaplotype;
tsLHCnts *pLHCnts;
CStats Stats;
double PValue;
m_LocalSeqErrRate = cDfltLocalSeqErrRate;

// generate BED header
m_LineBuffOffs = sprintf(m_pszLineBuff,"track name=\"%s\" description=\"%s\" useScore=1\n", m_szTrackName, m_szDescription);
CUtility::SafeWrite(m_hOutFile, m_pszLineBuff, m_LineBuffOffs);
m_LineBuffOffs = 0;

// iterate over all called SNP files and process each one
bBiSNPs = false;
ExpNumSNPs = 0;
for(FileIdx = 0; FileIdx < m_NumSNPFiles; FileIdx++)
	{
	pszSNPFile = m_ppszSNPFiles[FileIdx];
	if(!AcceptThisFile(pszSNPFile))
		{
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Not processing file: %s", pszSNPFile);
		continue;
		}
	if((Rslt = m_pCSV->Open(pszSNPFile)) != eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to open file: %s", pszSNPFile);
		Reset();
		return(Rslt);
		}
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Processing file: %s", pszSNPFile);
	NumElsParsed = 0;
	while((Rslt = m_pCSV->NextLine()) > 0)	// onto next line
		{
		NumFields = m_pCSV->GetCurFields();
		if(NumFields != 37 && NumFields != 92)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Expected 37 fields if disnp, or 92 fields if trisnp. Processed line with %d fields from '%s'", NumFields, pszSNPFile);
			Reset();
			return(eBSFerrFileType);
			}

		if(!NumElsParsed && m_pCSV->IsLikelyHeaderLine()) // slough header line
			continue;

		if(ExpNumSNPs == 0)
			{
			if(NumFields == 37)
				{
				ExpNumSNPs = 2;
				bBiSNPs = true;
				NumAlleleCombinations = 16;
				}
			else
				{
				ExpNumSNPs = 3;
				bBiSNPs = false;
				NumAlleleCombinations = 64;
				}
			}
		else
			{
			if(NumFields == 37 && !bBiSNPs)
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Previous files processed as %s, this file '%s' is a %s. Can't mix file types!", 
										bBiSNPs ? "bisnps" : "trisnps", bBiSNPs ? "trisnps" : "bisnps",pszSNPFile);
				Reset();
				return(eBSFerrFileType);
				}
			}


		NumElsParsed += 1;
		memset(&LocalisedHaplotype,0,sizeof(LocalisedHaplotype));
		m_pCSV->GetText(4, &pTxt);
		// a special case: UCSC references the reference SARS-CoV-2 as being NC_O45512v2,not NC_O45512.2
		// so change NC_O45512.2 to be NC_O45512v2
		if(!stricmp(pTxt, "NC_045512.2"))
			strcpy(LocalisedHaplotype.szChrom, "NC_045512v2");
		else
			strncpy(LocalisedHaplotype.szChrom,pTxt, cMaxGeneNameLen);
		LocalisedHaplotype.szChrom[cMaxGeneNameLen] = '\0';
		FldIdx = 5;
		pLHCnts = &LocalisedHaplotype.LHCnts[0];
		for(SNPIdx = 0; SNPIdx < ExpNumSNPs; SNPIdx++, pLHCnts++)
			{
			m_pCSV->GetInt(FldIdx++, &pLHCnts->Loci);
			m_pCSV->GetChar(FldIdx++, &BaseChr);
			switch(BaseChr) {
				case 'a': case 'A':
					pLHCnts->Base = eBaseA;
					break;
				case 'c': case 'C':
					pLHCnts->Base = eBaseC;
					break;
				case 'g': case 'G':
					pLHCnts->Base = eBaseG;
					break;
				case 't': case 'T':
					pLHCnts->Base = eBaseT;
					break;
				default:
					pLHCnts->Base = eBaseN;
					break;
				}
			for(BaseIdx = 0; BaseIdx <= 4; BaseIdx++)
				{
				m_pCSV->GetInt(FldIdx++, &pLHCnts->Cnts[BaseIdx]);
				pLHCnts->CovBases += pLHCnts->Cnts[BaseIdx];
				}
			}
		m_pCSV->GetInt(FldIdx++, &LocalisedHaplotype.DepthCnt);
		if(LocalisedHaplotype.DepthCnt < m_MinCoverage)		// have to have at least this many reads which were covering the haplotype
			continue;
		m_pCSV->GetInt(FldIdx++, &LocalisedHaplotype.AntisenseCnt);
		m_pCSV->GetInt(FldIdx++, &LocalisedHaplotype.NumHaplotypes);
		pCnt = LocalisedHaplotype.AlleleCombCnts;
		for(AlleleCombinationIdx = 0; AlleleCombinationIdx < NumAlleleCombinations; AlleleCombinationIdx++,pCnt++)
			{
			m_pCSV->GetInt(FldIdx++, pCnt);
			LocalisedHaplotype.SumAlleleCombCnts += *pCnt;
			}
		pLHCnts = &LocalisedHaplotype.LHCnts[0];
		for(SNPIdx = 0; SNPIdx < ExpNumSNPs; SNPIdx++, pLHCnts++)
			{
			for(BaseIdx = 0; BaseIdx < 4; BaseIdx++)
				{
				if(pLHCnts->Cnts[BaseIdx] == 0)
					continue;

				if(pLHCnts->Cnts[BaseIdx] / (double)pLHCnts->CovBases < m_MinAlleleProp)
					{
					pLHCnts->Cnts[BaseIdx] = 0;
					continue;
					}
				PValue = 1.0 - Stats.Binomial(pLHCnts->CovBases, pLHCnts->Cnts[BaseIdx], m_LocalSeqErrRate);
				if(PValue <= m_PValueThres)
					pLHCnts->bIsSNP = true;
				else
					pLHCnts->Cnts[BaseIdx] = 0;
				}
			}


		pCnt = LocalisedHaplotype.AlleleCombCnts;
		for(AlleleCombinationIdx = 0; AlleleCombinationIdx < NumAlleleCombinations; AlleleCombinationIdx++, pCnt++)
			{
			if(*pCnt == 0)
				continue;

			if(*pCnt / (double)LocalisedHaplotype.SumAlleleCombCnts < m_MinAlleleProp)
				{
				*pCnt = 0;
				continue;
				}
			PValue = 1.0 - Stats.Binomial(LocalisedHaplotype.SumAlleleCombCnts, *pCnt, m_LocalSeqErrRate);
			if(PValue < m_PValueThres)
				LocalisedHaplotype.NumAllelicTypes++;
			else
				*pCnt = 0;
			}

		if(LocalisedHaplotype.NumAllelicTypes < 1)
			continue;
		int StartLoci = LocalisedHaplotype.LHCnts[0].Loci;
		int EndLoci = LocalisedHaplotype.LHCnts[ExpNumSNPs-1].Loci;
		int LocHapLen = 1+EndLoci - StartLoci;
		int Score;
		Score = min(1000,(LocalisedHaplotype.NumAllelicTypes * 100) + 50);

		int RGB[3];
		switch(LocalisedHaplotype.NumAllelicTypes) {
			case 1:
				RGB[0] = 0;
				RGB[1] = 0;
				RGB[2] = 127;
				break;
			case 2:
				RGB[0] = 0;
				RGB[1] = 0;
				RGB[2] = 255;
				break;
			case 3:
				RGB[0] = 0;
				RGB[1] = 64;
				RGB[2] = 255;
				break;
			case 4:
				RGB[0] = 0;
				RGB[1] = 127;
				RGB[2] = 255;
				break;
			case 5:
				RGB[0] = 0;
				RGB[1] = 255;
				RGB[2] = 255;
				break;
			case 6:
				RGB[0] = 64;
				RGB[1] = 0;
				RGB[2] = 255;
				break;
			case 7:
				RGB[0] = 64;
				RGB[1] = 0;
				RGB[2] = 0;
				break;
			case 8:
				RGB[0] = 127;
				RGB[1] = 0;
				RGB[2] = 0;
				break;
			default:
				RGB[0] = 255;
				RGB[1] = 0;
				RGB[2] = 0;
				break;
			}
		m_LineBuffOffs += sprintf(&m_pszLineBuff[m_LineBuffOffs],"%s\t%d\t%d\t%dht\t%d\t.\t%d\t%d\t%d,%d,%d\t1\t%d\t%d\n",
								 LocalisedHaplotype.szChrom, StartLoci, EndLoci+1, 
								 LocalisedHaplotype.NumAllelicTypes, Score, StartLoci, EndLoci+1,
								 RGB[0], RGB[1], RGB[2], LocHapLen,0);
		if((m_LineBuffOffs + 1000) > m_AllocdLineBuff)
			{
			CUtility::SafeWrite(m_hOutFile, m_pszLineBuff, m_LineBuffOffs);
			m_LineBuffOffs = 0;
			}
		}
	m_pCSV->Close();
	}
if(m_LineBuffOffs)
	{
	CUtility::SafeWrite(m_hOutFile, m_pszLineBuff, m_LineBuffOffs);
	m_LineBuffOffs = 0;
	}

#ifdef _WIN32
_commit(m_hOutFile);
#else
fsync(m_hOutFile);
#endif
close(m_hOutFile);
m_hOutFile = -1;
return(0);
}

int 
CLocHap2Bed::Process2CSV (void)						// processing to CSV format
{
return(0);
}


int	
CLocHap2Bed::CompileChromRegExprs(int NumREIncludeFiles,	// number of file name regular expressions to include in processing
								  char** ppszREIncludeFiles,	// array of file name regular expressions to include in processing
								  int	NumREExcludeFiles,	// number of file name expressions to exclude from processing
								  char** ppszREExcludeFiles)	// array of file name regular expressions to exclude from processing
{
int Idx;

#ifndef _WIN32
int RegErr;				// regular expression parsing error
char szRegErr[128];		// to hold RegErr as textual representation ==> regerror();
#endif
m_NumREIncludeFiles= NumREIncludeFiles;
m_NumREExcludeFiles= NumREExcludeFiles;

#ifdef _WIN32
if(NumREIncludeFiles)
{
try {
	for(Idx = 0; Idx < NumREIncludeFiles; Idx++)
		{
		m_REIncludeFiles[Idx] = new Regexp();
		m_REIncludeFiles[Idx]->Parse(ppszREIncludeFiles[Idx], false);	// note case insensitive
		}
	}
catch(...)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to process include regexpr file name '%s'", ppszREIncludeFiles[Idx]);
	return(eBSFerrMem);
	}
}

if(NumREExcludeFiles)
{
try {
	for(Idx = 0; Idx < NumREExcludeFiles; Idx++)
		{
		m_REExcludeFiles[Idx] = new Regexp();
		m_REExcludeFiles[Idx]->Parse(ppszREExcludeFiles[Idx], false);	// note case insensitive
		}
	}
catch(...)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to process exclude regexpr file '%s' name", ppszREExcludeFiles[Idx]);
	return(eBSFerrMem);
	}
}
#else
if(NumREIncludeFiles)
{
for(Idx = 0; Idx < NumREIncludeFiles; Idx++)
	{
	RegErr = regcomp(&m_REIncludeFiles[Idx], ppszREIncludeFiles[Idx], REG_EXTENDED | REG_ICASE);	// note case insensitive
	if(RegErr)
		{
		regerror(RegErr, &m_REIncludeFiles[Idx], szRegErr, sizeof(szRegErr));
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to process include file regexpr '%s' error: %s", ppszREIncludeFiles[Idx], szRegErr);
		return(eBSFerrMem);
		}
	}
}
if(NumREExcludeFiles)
{
for(Idx = 0; Idx < NumREExcludeFiles; Idx++)
	{
	RegErr = regcomp(&m_REExcludeFiles[Idx], ppszREExcludeFiles[Idx], REG_EXTENDED | REG_ICASE);	// note case insensitive
	if(RegErr)
		{
		regerror(RegErr, &m_REExcludeFiles[Idx], szRegErr, sizeof(szRegErr));
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to process exclude file regexpr '%s' error: %s", ppszREExcludeFiles[Idx], szRegErr);
		return(eBSFerrMem);
		}
	}
}
#endif
return(eBSFSuccess);
}



bool					// true if file is accepted for processing, false if file not accepted
CLocHap2Bed::AcceptThisFile(char *pszFileName)
{
int IncFileIdx;
int ExclFileIdx;
bool bProcFile = false;
int MatchesFiltOut = 0;

if(!(m_NumREExcludeFiles || m_NumREIncludeFiles))
	return(true);

#ifdef _WIN32
RegexpMatch mc;
#else
regmatch_t mc;
int RegErr;					// regular expression parsing error
char szRegErr[128];			// to hold RegErr as textual representation ==> regerror();
#endif

		// check if to be excluded
bProcFile = true;
for(ExclFileIdx = 0; ExclFileIdx < m_NumREExcludeFiles; ExclFileIdx++)
	{
#ifdef _WIN32
	if(m_REExcludeFiles[ExclFileIdx]->Match(pszFileName, &mc))
#else
	if(!regexec(&m_REExcludeFiles[ExclFileIdx], pszFileName, 1, &mc, 0))
#endif
		{
		bProcFile = false;
		break;
		}
	}

	// to be included?
if(bProcFile && m_NumREIncludeFiles > 0)
	{
	bProcFile = false;
	for(IncFileIdx = 0; IncFileIdx < m_NumREIncludeFiles; IncFileIdx++)
		{
#ifdef _WIN32
		if(m_REIncludeFiles[IncFileIdx]->Match(pszFileName, &mc))
#else
		if(!regexec(&m_REIncludeFiles[IncFileIdx], pszFileName, 1, &mc, 0))
#endif
			{
			bProcFile = true;
			break;
			}
		}
	}

return(bProcFile);
}


int 
CLocHap2Bed::Process(eModeLocHap Mode,				// processing mode
								  int MinCoverage,				// must be at least this coverage at SNP site
								  double MinAlleleProp,			// putative allele must be at least this proportion of total site read coverage
								  double PValueThres,			// only accept SNP alleles which have a PValue <= this threshold
								  char* pszTrackName,			// track name
								  char* pszAssemblyName,		// UCSC assembly name - for SARS-CoV-2 it is "wuhCor1"
								  char* pszExperimentDescr,		// describes experiment
								  int NumInWhitelist,			// number of white listed file names
								  char** ppszWhitelisted,		// names of white listed files
								  int NumInBlackList,			// number of black listed files
								  char** ppsBlacklisted,		// names of black listed files
								  int NumSNPFiles,				// number of input SNP files
								  char** ppszSNPFiles,			// input SNP files
								  char* pszOutFile)				// output SNPs to this UCSC Personal Genome SNP format file
{
int Rslt;
int Len;
Reset();
m_NumSNPFiles = NumSNPFiles;
m_ppszSNPFiles = ppszSNPFiles;
strcpy (m_szOutBEDFile, pszOutFile);
strcpy (m_SpecAssemblyName, pszAssemblyName);
strcpy (m_szDescription, pszExperimentDescr);
strcpy (m_szTrackName, pszTrackName);
m_PValueThres = PValueThres;
m_MinCoverage = MinCoverage;
m_MinAlleleProp = MinAlleleProp;

if((Rslt=CompileChromRegExprs(NumInWhitelist, ppszWhitelisted, NumInBlackList, ppsBlacklisted)) != eBSFSuccess)
	return(Rslt);

if ((m_pszLineBuff = new char[cLHAllocLineBuffSize]) == NULL)
	{
	gDiagnostics.DiagOut (eDLFatal, gszProcName, "Process: Unable to allocate memory for output line buffering -- %s", strerror (errno));
	return(eBSFerrMem);
	}
m_AllocdLineBuff = cLHAllocLineBuffSize;

// check file extension, if '.csv' then generate output formated for CSV instead of the default BED format
if ((Len = (int)strlen (pszOutFile)) >= 5 && !stricmp (&pszOutFile[Len - 4], ".csv"))
	m_ReportFormat = eRHTCsv;
else
	m_ReportFormat = eRHTBed;

#ifdef _WIN32
if ((m_hOutFile = open (pszOutFile, (_O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC), (_S_IREAD | _S_IWRITE))) == -1)
#else
if ((m_hOutFile = open (pszOutFile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE)) == -1)
#endif
	{
	gDiagnostics.DiagOut (eDLFatal, gszProcName, "Process: Unable to create/truncate output file - '%s' - %s", pszOutFile, strerror (errno));
	Reset ();
	return(eBSFerrOpnFile);
	}

if ((m_pCSV = new CCSVFile) == NULL)
	{
	gDiagnostics.DiagOut (eDLFatal, gszProcName, "Unable to instantiate CCSVfile");
	Reset ();
	return(eBSFerrObj);
	}

switch (Mode) {
	case eRHTBed:
		Rslt = Process2BED();
		break;
	case eRHTCsv:
		Rslt = Process2CSV();
		break;
}

Reset ();
return(Rslt);
}
