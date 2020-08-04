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
#include "sarscov2ml.h"

int Process (eModeSC2 Mode,					// processing mode
	 char *pszMatrixFile,				// input matrix file
	 char *pszIsolateClassFile,				// input isolate classification file
	 char *pszOutFile);						// output feature classifications file

#ifdef _WIN32
int sarscov2ml(int argc, char *argv[])
{
	// determine my process name
	_splitpath (argv[0], NULL, NULL, gszProcName, NULL);
#else
int
sarscov2ml(int argc, char **argv)
{
	// determine my process name
	CUtility::splitpath ((char *)argv[0], NULL, gszProcName);
#endif
	int iFileLogLevel;			// level of file diagnostics
	int iScreenLogLevel;		// level of file diagnostics
	char szLogFile[_MAX_PATH];	// write diagnostics to this file
	int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure
	int NumberOfProcessors;		// number of installed CPUs

	 eModeSC2 PMode;					// processing mode
	 char szMatrixFile[_MAX_PATH]; // input matrix file
	 char szIsolateClassFile[_MAX_PATH]; // input isolate classification file
	 char szOutFile[_MAX_PATH];			 // output feature associations file

	struct arg_lit *help = arg_lit0 ("h", "help", "print this help and exit");
	struct arg_lit *version = arg_lit0 ("v", "version,ver", "print version information and exit");
	struct arg_int *FileLogLevel = arg_int0 ("f", "FileLogLevel", "<int>", "Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
	struct arg_file *LogFile = arg_file0 ("F", "log", "<file>", "diagnostics log file");

	struct arg_int *pmode = arg_int0 ("m", "mode", "<int>", "processing mode: 0 default and currently only processing mode externally exposed");
	struct arg_file *matrixfile = arg_file1("i", "in", "<file>", "Load matrix from this file");
	struct arg_file *isolateclassfile = arg_file0("I","isolateclass", "<file>", "Load isolate classification from this file");
	struct arg_file *outfile = arg_file1 ("o", "out", "<file>", "output file");
	struct arg_end *end = arg_end (200);

	void *argtable[] = { help,version,FileLogLevel,LogFile,
						pmode,matrixfile,isolateclassfile,outfile,end };

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

		PMode = pmode->count ? (eModeSC2)pmode->ival[0] : eMSC2default;
		if (PMode < eMSC2default || PMode >= eMSC2PlaceHolder)
		{
			gDiagnostics.DiagOut (eDLFatal, gszProcName, "Error: Processing mode '-m%d' specified outside of range %d..%d\n", PMode, eMSC2default, eMSC2PlaceHolder-1);
			exit (1);
		}


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

		strcpy (szMatrixFile, matrixfile->filename[0]);
		CUtility::TrimQuotedWhitespcExtd (szMatrixFile);
		if (szMatrixFile[0] == '\0')
			{
			gDiagnostics.DiagOut (eDLFatal, gszProcName, "No matrix file specified");
			exit (1);
			}

		if(isolateclassfile->count)
			{
			strcpy (szIsolateClassFile, isolateclassfile->filename[0]);
			CUtility::TrimQuotedWhitespcExtd (szIsolateClassFile);
			if (szIsolateClassFile[0] == '\0')
				{
				gDiagnostics.DiagOut (eDLFatal, gszProcName, "No isolate classification file specified");
				exit (1);
				}
			}
		else
			szIsolateClassFile[0] = '\0';

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
			case eMSC2default:
				pszDescr = "Training";
				break;
		}



		gDiagnostics.DiagOutMsgOnly (eDLInfo, "Processing matrix : '%s'", pszDescr);
		
		gDiagnostics.DiagOutMsgOnly (eDLInfo, "Input matrix file : '%s'", szMatrixFile);

		if(szIsolateClassFile[0] != '\0')
			gDiagnostics.DiagOutMsgOnly (eDLInfo, "Input isolate classification file : '%s'", szIsolateClassFile);

		gDiagnostics.DiagOutMsgOnly (eDLInfo, "Output file : '%s'", szOutFile);


#ifdef _WIN32
		SetPriorityClass (GetCurrentProcess (), BELOW_NORMAL_PRIORITY_CLASS);
#endif
		gStopWatch.Start ();
		Rslt = 0;
		Rslt = Process (PMode,					// processing mode
						szMatrixFile,		// input matrix file
						szIsolateClassFile,		// input association file
						szOutFile);				// output feature associations file
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

int Process(eModeSC2 Mode,					// processing mode
			char* pszMatrixFile,				// input matrix file
			char* pszIsolateClassFile,				// input isolate classification file
			char* pszOutFile)						// output feature classifications file
{
int Rslt;
CSarsCov2ML *pCSarsCov2ML;
if ((pCSarsCov2ML = new CSarsCov2ML) == NULL)
	{
	gDiagnostics.DiagOut (eDLFatal, gszProcName, "Unable to instantiate instance of CSarsCov2ML");
	return(eBSFerrInternal);
	}
Rslt = pCSarsCov2ML->Process(Mode,pszMatrixFile,pszIsolateClassFile,pszOutFile);

if (pCSarsCov2ML != NULL)
	delete pCSarsCov2ML;
return(Rslt);
}



CSarsCov2ML::CSarsCov2ML(void)
{
m_pMatrix = NULL;
Reset();
}

CSarsCov2ML::~CSarsCov2ML(void)
{
if(m_pMatrix != NULL)
	delete []m_pMatrix;
}

void
CSarsCov2ML::Reset(void) // reset state back to that immediately following class instantiation
{
if(m_pMatrix != NULL)
	{
	delete []m_pMatrix;
	m_pMatrix = NULL;
	}
m_NumCols = 0;								// matrix has this number of columns
m_NumRows = 0;								// matrix has this number of rows

m_NumFeatureNames = 0;
m_szFeatureNames[0] = '\0';
m_szFeatureIdx[0] = 0;
m_NxtszFeatureIdx = 0;

m_NumReadsetNames = 0;
m_szReadsetNames[0] = '\0';
m_szReadsetIdx[0] = 0;
m_NxtszReadsetIdx = 0;
memset(m_Classifications,0,sizeof(m_Classifications));
memset(m_ReadsetClassification,0,sizeof(m_ReadsetClassification));
m_TotClassified = 0;
m_NumClassificationNames = 0;
m_NxtszClassificationIdx = 0;
m_szClassificationNames[0] = '\0';
m_nCombination = 0;
m_rCombination = 0;
}


int 
CSarsCov2ML::LoadMatrixDimensions(char* pszMatrixFile, // matrix file to dimension
					uint32_t *pCols,	 // number of columns - includes extra column holding row name identifiers
					uint32_t *pRows)	 // number of rows - includes extra row holding column name identifiers
{
int Rslt;
int NumFields;
int ExpFields;
int LineNum;
CCSVFile *pCSV;
if(pCols != NULL)
	*pCols = 0;
if(pRows != NULL)
	*pRows = 0;
if((pCSV = new CCSVFile) == NULL)
	{
	gDiagnostics.DiagOut (eDLFatal, gszProcName, "Unable to instantiate instance of CCSVFile");
	Reset();
	return(eBSFerrInternal);
	}
pCSV->SetMaxFields(cMaxMatrixCols);
if((Rslt = pCSV->Open(pszMatrixFile)) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to open file: %s", pszMatrixFile);
	Reset();
	return(Rslt);
	}
ExpFields = 0;
NumFields = 0;
LineNum = 0;
while((Rslt = pCSV->NextLine()) > 0)	// onto next line
	{
	LineNum++;
	NumFields = pCSV->GetCurFields();
	if(ExpFields > 0 && NumFields != ExpFields)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Inconsistency in number of fields (expected %d, parsed %d) near line %d in file: '%s'", ExpFields,NumFields,LineNum, pszMatrixFile);
		delete pCSV;
		Reset();
		return(eBSFerrFieldCnt);
		}

	if(LineNum == 1)
		{
		if(pCSV->IsLikelyHeaderLine())		// parse header line as this contains the features (loci)
			{
			ExpFields = NumFields;
			continue;
			}
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Expected header row at line 1 in file: '%s'",pszMatrixFile);
		delete pCSV;
		Reset();
		return(eBSFerrParse);
		}
	continue;
	}
if(pCols != NULL)
	*pCols = NumFields;
if(pRows != NULL)
	*pRows = LineNum;
delete pCSV;
return(eBSFSuccess);
}

int
CSarsCov2ML::LoadMatrix(char* pszMatrixFile)	// matrix file to load
{
int Rslt;
int NumFields;
int LineNum;
int FieldIdx;
char *pTxt;
int FeatID;
int ReadsetID;
int Value;
uint32_t *pValue;
CCSVFile *pCSV;

if(m_pMatrix != NULL)
	{
	delete []m_pMatrix;
	m_pMatrix = NULL;
	}
m_NumCols = 0;
m_NumRows = 0;

if((Rslt=LoadMatrixDimensions(pszMatrixFile,&m_NumCols,&m_NumRows))!=eBSFSuccess)
	return(Rslt);

if((m_pMatrix = new uint32_t[(size_t)m_NumCols * m_NumRows]) == NULL)
	{
	Reset();
	return(eBSFerrMem);
	}
memset(m_pMatrix,0,(size_t)m_NumCols * m_NumRows * sizeof(*m_pMatrix));

if((pCSV = new CCSVFile) == NULL)
	{
	gDiagnostics.DiagOut (eDLFatal, gszProcName, "Unable to instantiate instance of CCSVFile");
	return(eBSFerrInternal);
	}

pCSV->SetMaxFields(m_NumCols+1);
if((Rslt = pCSV->Open(pszMatrixFile)) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to open file: %s", pszMatrixFile);
	Reset();
	return(Rslt);
	}
gDiagnostics.DiagOut(eDLInfo, gszProcName, "Processing file: %s", pszMatrixFile);
LineNum = 0;
pValue = &m_pMatrix[1];
while((Rslt = pCSV->NextLine()) > 0)	// onto next line
	{
	LineNum++;
	NumFields = pCSV->GetCurFields();
	if(NumFields != m_NumCols)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Inconsistency in number of fields (expected %d, parsed %d) near line %d in file: '%s'", m_NumCols+1,NumFields,LineNum, pszMatrixFile);
		delete pCSV;
		Reset();
		return(eBSFerrFieldCnt);
		}
	if(LineNum == 1)
		{
		if(pCSV->IsLikelyHeaderLine())		// parse header line as this contains the features (loci)
			{
			for(FieldIdx = 2; FieldIdx <= NumFields; FieldIdx++)
				{
				pCSV->GetText(FieldIdx, &pTxt);
				if((FeatID = AddFeature(pTxt)) < 1)
					{
					gDiagnostics.DiagOut(eDLFatal, gszProcName, "Failed adding feature '%s' from file: '%s'",pTxt,pszMatrixFile);
					delete pCSV;
					Reset();
					return(FeatID);
					}
				*pValue++ = FeatID;
				}
			continue;
			}
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Expected header row at line 1 in file: '%s'",pszMatrixFile);
		delete pCSV;
		Reset();
		return(eBSFerrParse);
		}

	// parse out the readset name
	pCSV->GetText(1, &pTxt);
	if((ReadsetID = AddReadset(pTxt)) < 1)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Failed adding readset '%s' from file: '%s'",pTxt,pszMatrixFile);
		delete pCSV;
		Reset();
		return(ReadsetID);
		}
	*pValue++ = ReadsetID;

	// now for the features values
	for(FieldIdx = 2; FieldIdx <= NumFields; FieldIdx++)
		{
		pCSV->GetInt(FieldIdx, &Value);
		*pValue++ = Value;
		}

	}

delete pCSV;
return(eBSFSuccess);
}


int 
CSarsCov2ML::LoadClassifications(char* pszClassFile)	// classifications file to load
{
int Rslt;
int LineNum;
int ReadsetID;
uint32_t ClassificationID;
char *pTxt;
CCSVFile *pCSV;
if((pCSV = new CCSVFile) == NULL)
	{
	gDiagnostics.DiagOut (eDLFatal, gszProcName, "Unable to instantiate instance of CCSVFile");
	return(eBSFerrInternal);
	}
if((Rslt = pCSV->Open(pszClassFile)) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to open file: %s", pszClassFile);
	Reset();
	return(Rslt);
	}
gDiagnostics.DiagOut(eDLInfo, gszProcName, "Processing file: %s", pszClassFile);
LineNum = 0;
while((Rslt = pCSV->NextLine()) > 0)	// onto next line
	{
	LineNum++;
	if(LineNum == 1 && pCSV->IsLikelyHeaderLine())		// skip header line
		continue;
	pCSV->GetText(1, &pTxt);			// expected to contain the readset prefix
	if((ReadsetID = LocateReadset(pTxt)) < 1)
		{
		gDiagnostics.DiagOut(eDLWarn, gszProcName, "Unable to locate matching readset for '%s' from file: '%s'",pTxt,pszClassFile);
		continue;
		}

	pCSV->GetText(2, &pTxt);			// expected to contain the associated training classification
	if(pTxt == NULL || pTxt[0] == '\0' || !stricmp("missing",pTxt))	// unclassified?
		continue;

	if((ClassificationID = AddClassification(pTxt)) < 1)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to add classification for '%s' from file: '%s'",pTxt,pszClassFile);
		delete pCSV;
		Reset();
		return(eBSFerrEntryCreate);
		}
	m_TotClassified += 1;
	m_Classifications[ClassificationID-1].NumReadsets += 1;
	m_ReadsetClassification[ReadsetID-1].ReadsetID = ReadsetID;
	m_ReadsetClassification[ReadsetID-1].ClassificationID = ClassificationID;
	}
for(ClassificationID = 0; ClassificationID < m_NumClassificationNames; ClassificationID++)
	m_Classifications[ClassificationID].Proportion = m_Classifications[ClassificationID].NumReadsets / (double)m_TotClassified;
delete pCSV;


return(eBSFSuccess);
}

uint32_t		// returned feature (loci) identifier, 0 if unable to accept this feature name
CSarsCov2ML::AddFeature(char* pszFeature) // associate unique identifier with this feature name
{
uint32_t FeatureNameIdx;
uint32_t FeatureNameLen;

// iterate over all known Features in case this Feature to add is a duplicate
for(FeatureNameIdx = 0; FeatureNameIdx < m_NumFeatureNames; FeatureNameIdx++)
	if(!stricmp(pszFeature, &m_szFeatureNames[m_szFeatureIdx[FeatureNameIdx]]))
		return(FeatureNameIdx + 1);

// Feature is not a duplicate
FeatureNameLen = (int)strlen(pszFeature);
if((m_NxtszFeatureIdx +FeatureNameLen + 1) > (int)sizeof(m_szFeatureNames))
	return(eBSFerrMaxEntries);
if(m_NumFeatureNames == cMaxMatrixCols)
	return(0);

m_szFeatureIdx[m_NumFeatureNames] = m_NxtszFeatureIdx;
strcpy(&m_szFeatureNames[m_NxtszFeatureIdx], pszFeature);
m_NxtszFeatureIdx += FeatureNameLen + 1;
return(++m_NumFeatureNames);
}

char* 
CSarsCov2ML::LocateFeature(uint32_t FeatureID)
{
if(FeatureID < 1 || FeatureID > m_NumFeatureNames)
	return(NULL);
return(&m_szFeatureNames[m_szFeatureIdx[FeatureID-1]]);
}

uint32_t 
CSarsCov2ML::LocateFeature(char *pszFeaturePrefix)		// match on this prefix, may be full length
{
int Len;
uint32_t FeatureNameIdx;
if(pszFeaturePrefix == NULL || pszFeaturePrefix[0] == '\0')
	return(0);
Len = (int)strlen(pszFeaturePrefix);

for(FeatureNameIdx = 0; FeatureNameIdx < m_NumFeatureNames; FeatureNameIdx++)
	if(!strnicmp(pszFeaturePrefix, &m_szFeatureNames[m_szFeatureIdx[FeatureNameIdx]],Len))
		return(FeatureNameIdx + 1);
return(0);
}

uint32_t		// returned readset identifier, 0 if unable to accept this readset name
CSarsCov2ML::AddReadset(char* pszReadset) // associate unique identifier with this readset name
{
uint32_t ReadsetIdx;
int ReadsetNameLen;

// iterate over all known readset in case this readset to add is a duplicate
for(ReadsetIdx = 0; ReadsetIdx < m_NumReadsetNames; ReadsetIdx++)
	if(!stricmp(pszReadset, &m_szReadsetNames[m_szReadsetIdx[ReadsetIdx]]))
		return(ReadsetIdx + 1);

// readset is not a duplicate
ReadsetNameLen = (int)strlen(pszReadset);
if(((size_t)m_NxtszReadsetIdx + ReadsetNameLen + 1) > sizeof(m_szReadsetNames))
	return(0);
if(m_NumReadsetNames == cMaxMatrixRows)
	return(0);

m_szReadsetIdx[m_NumReadsetNames] = m_NxtszReadsetIdx;
strcpy(&m_szReadsetNames[m_NxtszReadsetIdx], pszReadset);
m_NxtszReadsetIdx += ReadsetNameLen + 1;
return(++m_NumReadsetNames);
}

char* // returned ptr to readset name 
CSarsCov2ML::LocateReadset(uint32_t ReadsetID) // readset name identifier
{
if(ReadsetID < 1 || ReadsetID > m_NumReadsetNames)
	return(NULL);
return(&m_szReadsetNames[m_szReadsetIdx[ReadsetID - 1]]);
}

uint32_t 
CSarsCov2ML::LocateReadset(char *pszReadsetPrefix)		// match on this prefix, may be full length
{
int Len;
uint32_t ReadsetIdx;
if(pszReadsetPrefix == NULL || pszReadsetPrefix[0] == '\0')
	return(0);
Len = (int)strlen(pszReadsetPrefix);

for(ReadsetIdx = 0; ReadsetIdx < m_NumReadsetNames; ReadsetIdx++)
	if(!strnicmp(pszReadsetPrefix, &m_szReadsetNames[m_szReadsetIdx[ReadsetIdx]],Len))
		return(ReadsetIdx + 1);
return(0);
}

uint32_t		// returned classification identifier, 0 if unable to accept this classification name
CSarsCov2ML::AddClassification(char* pszClassification) // associate unique identifier with this classification name
{
uint32_t ClassificationIdx;
int ClassificationNameLen;

// iterate over all known Classifications in case this Classification to add is a duplicate
for(ClassificationIdx = 0; ClassificationIdx < m_NumClassificationNames; ClassificationIdx++)
	if(!stricmp(pszClassification, &m_szClassificationNames[m_Classifications[ClassificationIdx].szClassificationIdx]))
		return(ClassificationIdx + 1);

// readset is not a duplicate
ClassificationNameLen = (int)strlen(pszClassification);
if(((size_t)m_NxtszClassificationIdx + ClassificationNameLen + 1) > sizeof(m_szClassificationNames))
	return(eBSFerrMaxEntries);
if(m_NumClassificationNames == cMaxClassifications)
	return(0);

m_Classifications[m_NumClassificationNames].szClassificationIdx = m_NxtszClassificationIdx;
strcpy(&m_szClassificationNames[m_NxtszClassificationIdx], pszClassification);
m_NxtszClassificationIdx += ClassificationNameLen + 1;
return(++m_NumClassificationNames);
}

char* // returned ptr to Classification name 
CSarsCov2ML::LocateClassification(uint32_t ClassificationID) // Classification name identifier
{
if(ClassificationID < 1 || ClassificationID > m_NumClassificationNames)
	return(NULL);
return(&m_szClassificationNames[m_Classifications[ClassificationID - 1].szClassificationIdx]);
}



bool
CSarsCov2ML::Init_nCr(int n,			// initialise for n total elements
		 int r)			// from which will be drawing combinations of r elements
{
int Idx;
if(r <= 1 || r > cMaxR_nCr || n < r || n > cMaxN_nCr)
	{
	m_nCombination = 0;
	m_rCombination = 0;
	return(false);
	}

m_nCombination = n;
m_rCombination = r;
for(Idx = r-1; Idx >= 0; Idx--)
	m_nCrCombinations[Idx]=Idx;
return(true);
}


bool	// false if all combinations iterated, true if next combination generated
CSarsCov2ML::Iter_nCr(void)				// iterator for drawing next combination of r elements from n total
{
int Idx;
if(m_nCombination == 0 || m_rCombination == 0)
	return(false);

if(m_nCrCombinations[0] >= m_nCombination - m_rCombination)
	return(false);	// exhausted all combinations

for(Idx = m_rCombination-1; Idx >= 0; Idx--)
	{
	if(m_nCrCombinations[Idx] < m_nCombination - m_rCombination + Idx)
		{
		m_nCrCombinations[Idx] += 1;
		break;
		}
	}
while(Idx != m_rCombination-1)
	{
	m_nCrCombinations[Idx+1] = m_nCrCombinations[Idx] + 1;
	Idx++;
	}
return(true);
}




int
CSarsCov2ML::RunKernel(uint32_t NumLinkedFeatures,	// require this many features to be linked
				  uint32_t MinPropRows,				// require at least this many rows to show same linkage
				  uint32_t FeatType)				// linkage is between these minimum feature types
{

if(NumLinkedFeatures < 1 || NumLinkedFeatures > cMaxR_nCr)
	return(0);

uint32_t ClassIdx;
uint32_t RowIdx;
uint32_t ColIdx;
uint32_t *pCell;
uint32_t ClassifiedAs;
uint32_t NumRowsClassified;
uint32_t TotColCnts;
CStats Stats;
// iterate over feature loci and discover which loci are proportionally over/under represented relative to classification
ClassIdx = 0;
memset(m_TopLinkages,0,sizeof(m_TopLinkages));

if(0)
	{
	for(ColIdx = 1; ColIdx < m_NumCols; ColIdx++)
		{
		TotColCnts = 0;
		NumRowsClassified = 0;
		memset(m_ColClassifiedThresCnts,0,sizeof(m_ColClassifiedThresCnts));
		memset(m_ColClassifiedBelowCnts,0,sizeof(m_ColClassifiedBelowCnts));
		for(RowIdx = 1; RowIdx < m_NumRows; RowIdx++)
			{
			pCell = &m_pMatrix[RowIdx*m_NumCols];
			if(m_ReadsetClassification[*pCell-1].ReadsetID == 0)	// 0 if unclassified
				continue;
			NumRowsClassified += 1;
			ClassifiedAs = m_ReadsetClassification[*pCell-1].ClassificationID;
			pCell += ColIdx;
			if(*pCell >= FeatType)
				{
				m_ColClassifiedThresCnts[ClassifiedAs-1] += *pCell;
				TotColCnts += 1;
				}
			else
				m_ColClassifiedBelowCnts[ClassifiedAs-1] += 1;
			}
		if(NumRowsClassified < 50)		// currently just an arbitrary threshold - 
			continue;
		if(TotColCnts < NumRowsClassified/100)	// currently just an arbitrary threshold - 
			continue;
		double Prob = Stats.FishersExactTest(m_ColClassifiedThresCnts[0],m_ColClassifiedBelowCnts[0],m_ColClassifiedThresCnts[1],m_ColClassifiedBelowCnts[1]);
		if(Prob <= 0.01)
			{
			m_TopLinkages[ClassIdx].Prob = Prob;
			m_TopLinkages[ClassIdx].ColIdx = ColIdx;
			ClassIdx+=1;
			if(ClassIdx > cMaxN_nCr)
				break;
			}
		}
	}
else   // else looking for linkages which are not classified
	{
	for(ColIdx = 1; ColIdx < m_NumCols; ColIdx++)
		{
		TotColCnts = 0;
		NumRowsClassified = 0;
		memset(m_ColClassifiedThresCnts,0,sizeof(m_ColClassifiedThresCnts));
		memset(m_ColClassifiedBelowCnts,0,sizeof(m_ColClassifiedBelowCnts));
		for(RowIdx = 1; RowIdx < m_NumRows; RowIdx++)		// counting rows with that feature
			{
			pCell = &m_pMatrix[RowIdx*m_NumCols];
			NumRowsClassified += 1;
			pCell += ColIdx;
			if(*pCell >= FeatType)
				{
				m_ColClassifiedThresCnts[0] += *pCell;
				TotColCnts += 1;
				}
			else
				m_ColClassifiedBelowCnts[0] += 1;
			}
		if(NumRowsClassified < 50)		// currently just an arbitrary threshold - 
			continue;
		if(TotColCnts < 50)	// currently just an arbitrary threshold 
			continue;
		double Prob = 1.0 - Stats.Binomial(NumRowsClassified,m_ColClassifiedThresCnts[0],0.01);
		if(Prob <= 0.05)
			{
			m_TopLinkages[ClassIdx].Prob = Prob;
			m_TopLinkages[ClassIdx].ColIdx = ColIdx;
			ClassIdx+=1;
			if(ClassIdx > cMaxN_nCr)
				break;
			}
		}
	}


if(ClassIdx < NumLinkedFeatures)
	return(0);

if(ClassIdx > cMaxN_nCr)	// can't handle more!
	ClassIdx = cMaxN_nCr;

// just playing :-)
uint32_t NumFeatsLinked;
uint32_t NumRowsLinked;
uint32_t *pLinkedCell;


Init_nCr(ClassIdx,NumLinkedFeatures);
do
	{
	// work to do here!
	NumRowsLinked = 0;
	for(RowIdx = 1; RowIdx < m_NumRows; RowIdx++)
		{
		pCell = &m_pMatrix[RowIdx*m_NumCols];
		if(m_ReadsetClassification[*pCell-1].ReadsetID == 0)	// 0 if unclassified
			continue;

		NumFeatsLinked = 0;
		for(uint32_t IdxR = 0; IdxR < m_rCombination; IdxR++)
			{
			pLinkedCell = pCell + m_TopLinkages[m_nCrCombinations[IdxR]].ColIdx;
			if(*pLinkedCell >= FeatType)
				NumFeatsLinked++;
			}
		if(NumFeatsLinked < NumLinkedFeatures)
			continue;
		NumRowsLinked++;
		}

	if(NumRowsLinked >= MinPropRows)		
		{
		int BuffIdx;
		char szBuff[10000];
		BuffIdx = sprintf(szBuff,"Have %d rows (thres %u) with %u features linked by this minimum feature type %d -->",NumRowsLinked,MinPropRows,NumLinkedFeatures, FeatType);
		NumFeatsLinked = 0;
		for(uint32_t IdxR = 0; IdxR < m_rCombination; IdxR++)
			{
			BuffIdx += sprintf(&szBuff[BuffIdx]," %s", LocateFeature(m_TopLinkages[m_nCrCombinations[IdxR]].ColIdx));
			NumFeatsLinked++;
			}
		gDiagnostics.DiagOut(eDLInfo, gszProcName,szBuff);
		}
	}
while(Iter_nCr());

uint32_t *pLinkedCells[10];
uint32_t CurLinkagesIdx;
uint32_t NxtLinkagesIdx;
uint32_t NumLinked;
uint32_t CheckMe;
CheckMe = 0;
for(CurLinkagesIdx = 0; CurLinkagesIdx < ClassIdx-1; CurLinkagesIdx++)
	{
	for(NxtLinkagesIdx = CurLinkagesIdx+1; NxtLinkagesIdx < ClassIdx; NxtLinkagesIdx++)
		{
		NumLinked = 0;
		for(RowIdx = 1; RowIdx < m_NumRows; RowIdx++)
			{
			pCell = &m_pMatrix[RowIdx*m_NumCols];
			if(m_ReadsetClassification[*pCell-1].ReadsetID == 0)	// 0 if unclassified
				continue;
			pLinkedCells[0] = pCell + m_TopLinkages[CurLinkagesIdx].ColIdx;
			if(*pLinkedCells[0] < FeatType)
				continue;
			pLinkedCells[1] = pCell + m_TopLinkages[NxtLinkagesIdx].ColIdx;
			if(*pLinkedCells[1] < FeatType)
				continue;
			NumLinked+=1;
			}
		if(NumLinked > 10)
			{
			char *pszFeat1 = LocateFeature(m_TopLinkages[CurLinkagesIdx].ColIdx);
			char *pszFeat2 = LocateFeature(m_TopLinkages[NxtLinkagesIdx].ColIdx);
//			gDiagnostics.DiagOut(eDLInfo, gszProcName, "%s and %s linked %d times with which %d",pszFeat1,pszFeat2,NumLinked,FeatType);
			CheckMe++;
			}
		}
	}
return(ClassIdx);
}

int 
CSarsCov2ML::Process(eModeSC2 Mode,					// processing mode
				   char* pszMatrixFile,			// input matrix file
				   char* pszIsolateClassFile,			// input isolate classification file
				   char* pszOutFile)					// output feature classifications file
{
int Rslt;
// load the input matrix containing row isolates and column loci
if((Rslt=LoadMatrix(pszMatrixFile))!=eBSFSuccess)
	return(Rslt);

// if present then parse in the classification file
if(pszIsolateClassFile != NULL && pszIsolateClassFile[0] != '\0')
	if((Rslt=LoadClassifications(pszIsolateClassFile))!=eBSFSuccess)
		return(Rslt);
int K2 = RunKernel(2,30,3);
int K3 = RunKernel(3,20,3);
int K4 = RunKernel(4,20,3);
int K5 = RunKernel(5,10,3);
int K6 = RunKernel(6,5,3);
int K7 = RunKernel(7,3,3);
int K8 = RunKernel(8,3,3);
int K9 = RunKernel(9,2,3);
Reset();
return(0);
}

