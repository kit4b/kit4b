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
	uint32_t NumLinkedFeatures,				// require this many features to be linked
	uint32_t MinLinkedRows,					// require at least this many rows to show same linkage
	uint32_t FeatClassValue,				// linkage is between these minimum feature class values
	 char *pszMatrixFile,					// input matrix file
	 char *pszIsolateClassFile,				// input isolate classification file
	 char *pszOutFile);						// output feature classifications file

#ifdef _WIN32
int sarscov2ml(int argc, char *argv[])
{
	// determine my process name
	_splitpath (argv[0], nullptr, nullptr, gszProcName, nullptr);
#else
int
sarscov2ml(int argc, char **argv)
{
	// determine my process name
	CUtility::splitpath ((char *)argv[0], nullptr, gszProcName);
#endif
	int iFileLogLevel;			// level of file diagnostics
	int iScreenLogLevel;		// level of file diagnostics
	char szLogFile[_MAX_PATH];	// write diagnostics to this file
	int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure
	int NumberOfProcessors;		// number of installed CPUs

	 eModeSC2 PMode;					// processing mode
	 uint32_t NumLinkedFeatures;		// require this many features to be linked
	uint32_t MinLinkedRows;				// require at least this many rows to show same linkage
	uint32_t FeatClassValue;			// linkage is between these minimum feature class values
	 char szMatrixFile[_MAX_PATH];		// input matrix file
	 char szIsolateClassFile[_MAX_PATH]; // input isolate classification file
	 char szOutFile[_MAX_PATH];			 // output feature associations file

	struct arg_lit *help = arg_lit0 ("h", "help", "print this help and exit");
	struct arg_lit *version = arg_lit0 ("v", "version,ver", "print version information and exit");
	struct arg_int *FileLogLevel = arg_int0 ("f", "FileLogLevel", "<int>", "Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
	struct arg_file *LogFile = arg_file0 ("F", "log", "<file>", "diagnostics log file");

	struct arg_int *pmode = arg_int0 ("m", "mode", "<int>", "processing mode: 0 discover feature linkages");

	struct arg_int *numlinkedfeatures = arg_int0 ("l", "NumLinkedFeatures", "<int>", "require this many features to be linked (default 5, range 1..50)");
	struct arg_int *minlinkedrows = arg_int0 ("r", "MinLinkedRows", "<int>", "require at least this many rows to show same linkage (default 50, range 2..100000)");
	struct arg_int *featclassvalue = arg_int0 ("c", "FeatClassValue", "<int>", "linkage is between these minimum feature class values (default 3, range 0..1000)");

	struct arg_file *matrixfile = arg_file1("i", "in", "<file>", "Load matrix from this file");
	struct arg_file *isolateclassfile = arg_file0("I","isolateclass", "<file>", "Load isolate feature classifications from this file");
	struct arg_file *outfile = arg_file1 ("o", "out", "<file>", "output file");
	struct arg_end *end = arg_end (200);

	void *argtable[] = { help,version,FileLogLevel,LogFile,
						pmode,numlinkedfeatures,minlinkedrows,featclassvalue,matrixfile,isolateclassfile,outfile,end };

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

		NumLinkedFeatures = numlinkedfeatures->count ? numlinkedfeatures->ival[0] : 5;
		if (NumLinkedFeatures < 1 || NumLinkedFeatures > cMaxR_nCr)
			{
			gDiagnostics.DiagOut (eDLFatal, gszProcName, "Error: Number of linked features '-l%d' specified outside of range %d..%d\n", NumLinkedFeatures, 1, cMaxR_nCr);
			exit (1);
			}

		MinLinkedRows = minlinkedrows->count ? minlinkedrows->ival[0] : 50;
		if (MinLinkedRows < 2 || MinLinkedRows > 100000)
		{
			gDiagnostics.DiagOut (eDLFatal, gszProcName, "Error: Number of rows with same property or linkage '-r%d' specified outside of range %d..%d\n", MinLinkedRows, 2, 100000);
			exit (1);
		}

		FeatClassValue = featclassvalue->count ? featclassvalue->ival[0] : 3;
		if (FeatClassValue < 0 || FeatClassValue > 1000)
		{
			gDiagnostics.DiagOut (eDLFatal, gszProcName, "Error: Minimum feature class value '-c%d' specified outside of range %d..%d\n", FeatClassValue, 0, 1000);
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
				pszDescr = "locate linkages between features";
				break;

			}

		gDiagnostics.DiagOutMsgOnly (eDLInfo, "Processing is to : '%s'", pszDescr);
		gDiagnostics.DiagOutMsgOnly (eDLInfo, "Linkage between this number of features : %d", NumLinkedFeatures);
		gDiagnostics.DiagOutMsgOnly (eDLInfo, "Minimum number samples with same linkage : %d", MinLinkedRows);
		gDiagnostics.DiagOutMsgOnly (eDLInfo, "Minimum feature class value : %d", FeatClassValue);

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
						NumLinkedFeatures,		// require this many features to be linked
						MinLinkedRows,				// require at least this many rows to show same linkage
						FeatClassValue,			// linkage is between these minimum feature types
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

int Process(eModeSC2 Mode,						// processing mode
			uint32_t NumLinkedFeatures,			// require this many features to be linked
			uint32_t MinPropRows,				// require at least this many rows to show same linkage
			uint32_t FeatClassValue,			// linkage is between these minimum feature class values
			char* pszMatrixFile,				// input matrix file
			char* pszIsolateClassFile,			// input isolate classification file
			char* pszOutFile)					// output feature classifications file
{
int Rslt;
CSarsCov2ML *pCSarsCov2ML;
if ((pCSarsCov2ML = new CSarsCov2ML) == nullptr)
	{
	gDiagnostics.DiagOut (eDLFatal, gszProcName, "Unable to instantiate instance of CSarsCov2ML");
	return(eBSFerrInternal);
	}
Rslt = pCSarsCov2ML->Process(Mode,NumLinkedFeatures,MinPropRows,FeatClassValue,pszMatrixFile,pszIsolateClassFile,pszOutFile);

if (pCSarsCov2ML != nullptr)
	delete pCSarsCov2ML;
return(Rslt);
}



CSarsCov2ML::CSarsCov2ML(void)
{
m_pMatrix = nullptr;
m_hOutFile = -1;
m_pOutBuffer = nullptr;
Reset();
}

CSarsCov2ML::~CSarsCov2ML(void)
{
if(m_pMatrix != nullptr)
	delete []m_pMatrix;
if(m_pOutBuffer != nullptr)
	delete []m_pOutBuffer;
if(m_hOutFile != -1)
	close(m_hOutFile);
}

void
CSarsCov2ML::Reset(void) // reset state back to that immediately following class instantiation
{
if(m_hOutFile != -1)
	{
	if(m_pOutBuffer != nullptr && m_OutBuffIdx)
		{
		CUtility::RetryWrites(m_hOutFile,m_pOutBuffer,m_OutBuffIdx);
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
if(m_pOutBuffer != nullptr)
	{
	delete []m_pOutBuffer;
	m_pOutBuffer = nullptr;
	}
m_OutBuffIdx = 0;
m_AllocOutBuff = 0;

if(m_pMatrix != nullptr)
	{
	delete []m_pMatrix;
	m_pMatrix = nullptr;
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
m_NthCombination = 0;
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
if(pCols != nullptr)
	*pCols = 0;
if(pRows != nullptr)
	*pRows = 0;
if((pCSV = new CCSVFile) == nullptr)
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
if(pCols != nullptr)
	*pCols = NumFields;
if(pRows != nullptr)
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

if(m_pMatrix != nullptr)
	{
	delete []m_pMatrix;
	m_pMatrix = nullptr;
	}
m_NumCols = 0;
m_NumRows = 0;

if((Rslt=LoadMatrixDimensions(pszMatrixFile,&m_NumCols,&m_NumRows))!=eBSFSuccess)
	return(Rslt);

if((m_pMatrix = new uint32_t[(size_t)m_NumCols * m_NumRows]) == nullptr)
	{
	Reset();
	return(eBSFerrMem);
	}
memset(m_pMatrix,0,(size_t)m_NumCols * m_NumRows * sizeof(*m_pMatrix));

if((pCSV = new CCSVFile) == nullptr)
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
uint32_t NumMissing;
uint32_t ClassificationID;
char *pTxt;
CCSVFile *pCSV;
if((pCSV = new CCSVFile) == nullptr)
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
NumMissing = 0;
while((Rslt = pCSV->NextLine()) > 0)	// onto next line
	{
	LineNum++;
	if(LineNum == 1 && pCSV->IsLikelyHeaderLine())		// skip header line
		continue;
	pCSV->GetText(1, &pTxt);			// expected to contain the readset prefix
	if((ReadsetID = LocateReadset(pTxt)) < 1)
		{
		gDiagnostics.DiagOut(eDLWarn, gszProcName, "Unable to locate matching entry in matrix for '%s' from file: '%s'",pTxt,pszClassFile);
		NumMissing++;
		continue;
		}

	pCSV->GetText(2, &pTxt);			// expected to contain the associated training classification
	if(pTxt == nullptr || pTxt[0] == '\0' || !stricmp("missing",pTxt))	// unclassified?
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
if(NumMissing)
	gDiagnostics.DiagOut(eDLWarn, gszProcName, "Unable to locate a total of %u readsets",NumMissing);
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
	return(nullptr);
return(&m_szFeatureNames[m_szFeatureIdx[FeatureID-1]]);
}

uint32_t 
CSarsCov2ML::LocateFeature(char *pszFeaturePrefix)		// match on this prefix, may be full length
{
int Len;
uint32_t FeatureNameIdx;
if(pszFeaturePrefix == nullptr || pszFeaturePrefix[0] == '\0')
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
	return(nullptr);
return(&m_szReadsetNames[m_szReadsetIdx[ReadsetID - 1]]);
}

uint32_t 
CSarsCov2ML::LocateReadset(char *pszReadsetPrefix)		// match on this prefix, may be full length
{
int Len;
uint32_t ReadsetIdx;
if(pszReadsetPrefix == nullptr || pszReadsetPrefix[0] == '\0')
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
	return(nullptr);
return(&m_szClassificationNames[m_Classifications[ClassificationID - 1].szClassificationIdx]);
}



bool
CSarsCov2ML::Init_nCr(int n,			// initialise for n total elements
		 int r)			// from which will be drawing combinations of r elements
{
int Idx;
if(r < 1 || r > cMaxR_nCr || n < r || n > cMaxN_nCr)
	{
	m_nCombination = 0;
	m_rCombination = 0;
	return(false);
	}

m_nCombination = n;
m_rCombination = r;
for(Idx = r-1; Idx >= 0; Idx--)
	m_nCrCombinations[Idx]=Idx;
m_NthCombination = 0;	// incremented with each call to Iter_nCr()
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

if(m_NthCombination++)	// if not first then iterate next combination
	{
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
	}
return(true);
}




int
CSarsCov2ML::RunKernel(eModeSC2 Mode,					// processing mode
					uint32_t NumLinkedFeatures,			// require this many features to be linked
					uint32_t MinLinkedRows,				// require at least this many rows to show same linkage
					uint32_t FeatClassValue)			// linkage is between these minimum feature class values
{

if(NumLinkedFeatures < 1 || NumLinkedFeatures > cMaxR_nCr)
	return(0);

uint32_t ClassIdx;
uint32_t RowIdx;
uint32_t ColIdx;
uint32_t *pCell;
uint32_t NumRowsClassified;
uint32_t LinkedRows;
CStats Stats;
// iterate over feature loci and discover which loci are proportionally over/under represented relative to classification
ClassIdx = 0;
memset(m_TopLinkages,0,sizeof(m_TopLinkages));

for(ColIdx = 1; ColIdx < m_NumCols; ColIdx++)
	{
	LinkedRows = 0;
	NumRowsClassified = 0;
	memset(m_ColClassifiedThresCnts,0,sizeof(m_ColClassifiedThresCnts));
	memset(m_ColClassifiedBelowCnts,0,sizeof(m_ColClassifiedBelowCnts));
	for(RowIdx = 1; RowIdx < m_NumRows; RowIdx++)		// counting rows with that feature
		{
		pCell = &m_pMatrix[RowIdx*m_NumCols];
		NumRowsClassified += 1;
		pCell += ColIdx;
		if(*pCell >= FeatClassValue)
			{
			m_ColClassifiedThresCnts[0] += 1;
			LinkedRows += 1;
			}
		else
			m_ColClassifiedBelowCnts[0] += 1;
		}
	if(LinkedRows < MinLinkedRows)
		continue;
	if(m_ColClassifiedThresCnts[0] == 0)
		continue;

	if(ClassIdx >= cMaxN_nCr)		// need to prune?
		{
		// attempt to find a linkage which is lower than this columns LinkedRows and if located then replace
		tsScoredCol *pLinkage = m_TopLinkages;
		uint32_t Idx;
		for(Idx = 0; Idx < ClassIdx; Idx++,pLinkage++)
			{
			if(pLinkage->LinkedRows < LinkedRows)
				{
				pLinkage->LinkedRows = LinkedRows;
				pLinkage->ColIdx = ColIdx;
				break;
				}
			}
		if(Idx == (int)cMaxN_nCr)
			continue;
		}
	m_TopLinkages[ClassIdx].LinkedRows = LinkedRows;
	m_TopLinkages[ClassIdx].ColIdx = ColIdx;
	ClassIdx+=1;
	}

if(ClassIdx < NumLinkedFeatures)
	return(0);

if(ClassIdx > cMaxN_nCr)	// can't handle more!
	ClassIdx = cMaxN_nCr;

qsort(m_TopLinkages, ClassIdx, sizeof(tsScoredCol), SortTopLinkages);

Init_nCr(ClassIdx,NumLinkedFeatures);
InitSarsCov2MLThreads(NumLinkedFeatures,MinLinkedRows,FeatClassValue,20);


return(ClassIdx);
}

int
CSarsCov2ML::ProcThreadML(tsThreadML *pPars)
{
// local per thread
int NumReportedLinkages;
uint32_t RowIdx;
uint32_t NumFeatsLinked;
uint32_t NumRowsLinked;
uint32_t *pLinkedCell;
uint32_t *pCell;
char *pszLoci;
uint32_t nCrCombinations[cMaxR_nCr];			// indexes for each potential nCr combination instance
uint32_t ClassifiedAsCnts[cMaxClassifications];

NumReportedLinkages = 0;
AcquireSerialise();
while(Iter_nCr())
	{
	memcpy(nCrCombinations,m_nCrCombinations,sizeof(nCrCombinations));
	ReleaseSerialise();
	if(m_NumClassificationNames)
		memset(ClassifiedAsCnts,0,sizeof(ClassifiedAsCnts));
	NumRowsLinked = 0;
	for(RowIdx = 1; RowIdx < m_NumRows; RowIdx++)
		{
		pCell = &m_pMatrix[RowIdx*m_NumCols];
		if(m_ReadsetClassification[*pCell-1].ReadsetID == 0)	// 0 if unlinked
			continue;

		NumFeatsLinked = 0;
		for(uint32_t IdxR = 0; IdxR < m_rCombination; IdxR++)
			{
			pLinkedCell = pCell + m_TopLinkages[nCrCombinations[IdxR]].ColIdx;
			if(*pLinkedCell >= pPars->FeatClassValue)
				NumFeatsLinked++;
			}
		if(NumFeatsLinked < pPars->NumLinkedFeatures)
			continue;
		if(m_NumClassificationNames)
			ClassifiedAsCnts[m_ReadsetClassification[*pCell-1].ClassificationID-1]++;
		NumRowsLinked++;
		}

	AcquireSerialise();
	if(NumRowsLinked >= pPars->MinLinkedRows)
		{
		m_OutBuffIdx += sprintf(&m_pOutBuffer[m_OutBuffIdx],"%u,%u,%u,%u,%d",m_NumRows,NumRowsLinked,pPars->MinLinkedRows,pPars->NumLinkedFeatures, pPars->FeatClassValue);
		NumFeatsLinked = 0;
		for(uint32_t IdxR = 0; IdxR < m_rCombination; IdxR++)
			{
			pszLoci = LocateFeature(m_TopLinkages[nCrCombinations[IdxR]].ColIdx);
			if(!strncmp(pszLoci,"Loci:",5))	// if loci using 'Loci:' as a chrom/seq name then strip this off as only used if processing a single unnamed sequence as in SARS-Cov-2
				pszLoci+=5;
			m_OutBuffIdx += sprintf(&m_pOutBuffer[m_OutBuffIdx],",\"%s\"", pszLoci);
			NumFeatsLinked++;
			}
		if(m_NumClassificationNames)
			{
			for(uint32_t IdxR = 0; IdxR < m_NumClassificationNames; IdxR++)
				{
				m_OutBuffIdx += sprintf(&m_pOutBuffer[m_OutBuffIdx],",%u", ClassifiedAsCnts[IdxR]);
				NumFeatsLinked++;
				}
			}
		m_OutBuffIdx += sprintf(&m_pOutBuffer[m_OutBuffIdx],"\n");
		if(m_OutBuffIdx + 1000 > m_AllocOutBuff)
			{
			CUtility::RetryWrites(m_hOutFile,m_pOutBuffer,m_OutBuffIdx);
			m_OutBuffIdx = 0;
			}
		NumReportedLinkages+=1;
		}
	}
ReleaseSerialise();
return(NumReportedLinkages);
}

#ifdef _WIN32
unsigned __stdcall ProcessMLThread(void * pThreadPars)
#else
void *ProcessMLThread(void * pThreadPars)
#endif
{
int Rslt;
tsThreadML *pPars = (tsThreadML *)pThreadPars;			// makes it easier not having to deal with casts!
CSarsCov2ML *pSarsCov2ML = (CSarsCov2ML *)pPars->pThis;
Rslt = pSarsCov2ML->ProcThreadML(pPars);
pPars->Rslt = Rslt;
#ifdef _WIN32
_endthreadex(0);
return(eBSFSuccess);
#else
pthread_exit(&pPars->Rslt);
#endif
}

int
CSarsCov2ML::InitSarsCov2MLThreads(	uint32_t NumLinkedFeatures,			// require this many features to be linked
				uint32_t MinLinkedRows,				// require at least this many rows to show same linkage
				uint32_t FeatClassValue,			// linkage is between these minimum feature class values
				int NumThreads)			// use this many threads
{
tsThreadML *pThreads;
tsThreadML *pThread;

if(CreateMutexes()!=eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Failed to create thread synchronisation mutexes");
	Reset();
	return(cBSFSyncObjErr);
	}

if ((pThreads = new tsThreadML[NumThreads]) == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "InitSarsCov2MLThreads: Memory allocation for thread context failed");
	return(eBSFerrMem);
	}
memset(pThreads,0,sizeof(tsThreadML) * NumThreads);
int ThreadIdx;
pThread = pThreads;
for (ThreadIdx = 1; ThreadIdx <= NumThreads; ThreadIdx++, pThread++)
	{
	pThread->ThreadIdx = ThreadIdx;
	pThread->pThis = this;

	pThread->NumLinkedFeatures=NumLinkedFeatures;
	pThread->MinLinkedRows=MinLinkedRows;
	pThread->FeatClassValue=FeatClassValue;

#ifdef _WIN32
	pThread->threadHandle = (HANDLE)_beginthreadex(nullptr, 0x0fffff, ProcessMLThread, pThread, 0, &pThread->threadID);
#else
	pThread->threadRslt = pthread_create(&pThread->threadID, nullptr, ProcessMLThread, pThread);
#endif
	}

// allow threads a few seconds to startup
#ifdef _WIN32
Sleep(5000);
#else
sleep(5);
#endif

pThread = pThreads;
for (ThreadIdx = 0; ThreadIdx < NumThreads; ThreadIdx++, pThread++)
	{
#ifdef _WIN32
	while (WAIT_TIMEOUT == WaitForSingleObject(pThread->threadHandle, 60000))
		{
		AcquireSerialise();

		ReleaseSerialise();
		};
	CloseHandle(pThread->threadHandle);
#else
	struct timespec ts;
	int JoinRlt;
	clock_gettime(CLOCK_REALTIME, &ts);
	ts.tv_sec += 60;
	while ((JoinRlt = pthread_timedjoin_np(pThread->threadID, nullptr, &ts)) != 0)
		{
		AcquireSerialise();

		ReleaseSerialise();
		ts.tv_sec += 60;
		}
#endif
	}

if(pThreads != nullptr)
	delete []pThreads;
DeleteMutexes();
return(0);
}


int 
CSarsCov2ML::Process(eModeSC2 Mode,						// processing mode
					uint32_t NumLinkedFeatures,			// require this many features to be linked
					uint32_t MinLinkedRows,				// require at least this many rows to show same linkage
					uint32_t FeatClassValue,			// linkage is between these minimum feature class values
					char* pszMatrixFile,				// input matrix file
					char* pszIsolateClassFile,			// input isolate classification file
					char* pszOutFile)					// output feature classifications file
{
int Rslt;
m_PMode = Mode;
strcpy(m_szMatrixFile,pszMatrixFile);
strcpy(m_szIsolateClassFile,pszIsolateClassFile);
strcpy(m_szOutFile,pszOutFile);
m_MinRowsClassified = MinLinkedRows;



// load the input matrix containing row isolates and column loci
if((Rslt=LoadMatrix(pszMatrixFile))!=eBSFSuccess)
	return(Rslt);

// if present then parse in the classification file
if(pszIsolateClassFile != nullptr && pszIsolateClassFile[0] != '\0')
	{
	if((Rslt=LoadClassifications(pszIsolateClassFile))!=eBSFSuccess)
		return(Rslt);
	}
else   // no classification file so default to a single unclassified class
	{
	uint32_t ClassificationID = AddClassification((char *)"Defaulted");
	m_Classifications[0].NumReadsets = m_NumReadsetNames;
	for(uint32_t ReadSetID = 0; ReadSetID < m_NumReadsetNames; ReadSetID++)
		{
		m_ReadsetClassification[ReadSetID].ReadsetID = ReadSetID+1;
		m_ReadsetClassification[ReadSetID].ClassificationID = ClassificationID;
		}
	m_Classifications[0].Proportion = 1.0;
	}


#ifdef _WIN32
m_hOutFile = open(pszOutFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE));
#else
if((m_hOutFile = open64(pszOutFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE)) != -1)
	if(ftruncate(m_hOutFile,0)!=0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create/truncate %s - %s",pszOutFile,strerror(errno));
		Reset();
		return(eBSFerrCreateFile);
		}
#endif
if(m_hOutFile < 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create/truncate %s - %s",pszOutFile,strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}

m_pOutBuffer = new char [cMaxAllocOutBuff];
m_AllocOutBuff = cMaxAllocOutBuff;

m_OutBuffIdx = sprintf(m_pOutBuffer,"\"TotSamples\",\"LinkedSamples\",\"MinSamples\",\"FeatLinks\",\"MinFeatValue\"");
for(uint32_t LinkIdx = 0; LinkIdx < NumLinkedFeatures; LinkIdx++)
	m_OutBuffIdx += sprintf(&m_pOutBuffer[m_OutBuffIdx],",\"L%u\"",LinkIdx);
for(uint32_t Idx = 0; Idx < m_NumClassificationNames; Idx++)
	m_OutBuffIdx += sprintf(&m_pOutBuffer[m_OutBuffIdx],",\"%s (%u)\"", LocateClassification(Idx+1),m_Classifications[Idx].NumReadsets);
m_OutBuffIdx += sprintf(&m_pOutBuffer[m_OutBuffIdx],"\n");
CUtility::RetryWrites(m_hOutFile,m_pOutBuffer,m_OutBuffIdx);
m_OutBuffIdx = 0;

int K = RunKernel(Mode,NumLinkedFeatures,MinLinkedRows,FeatClassValue);

if(m_OutBuffIdx)
	{
	CUtility::RetryWrites(m_hOutFile,m_pOutBuffer,m_OutBuffIdx);
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
Reset();
return(0);
}


int
CSarsCov2ML::CreateMutexes(void)
{
if(m_bMutexesCreated)
	return(eBSFSuccess);

#ifdef _WIN32
if((m_hMtxIterReads = CreateMutex(nullptr,false,nullptr))==nullptr)
	{
#else
if(pthread_mutex_init (&m_hMtxIterReads,nullptr)!=0)
	{
#endif
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to create mutex");
	return(eBSFerrInternal);
	}

m_bMutexesCreated = true;
return(eBSFSuccess);
}

void
CSarsCov2ML::DeleteMutexes(void)
{
if(!m_bMutexesCreated)
	return;
#ifdef _WIN32
CloseHandle(m_hMtxIterReads);
#else
pthread_mutex_destroy(&m_hMtxIterReads);
#endif
m_bMutexesCreated = false;
}


void
CSarsCov2ML::AcquireSerialise(void)
{
#ifdef _WIN32
WaitForSingleObject(m_hMtxIterReads,INFINITE);
#else
pthread_mutex_lock(&m_hMtxIterReads);
#endif
}

void
CSarsCov2ML::ReleaseSerialise(void)
{
#ifdef _WIN32
ReleaseMutex(m_hMtxIterReads);
#else
pthread_mutex_unlock(&m_hMtxIterReads);
#endif
}


int
CSarsCov2ML::SortTopLinkages(const void* arg1, const void* arg2)
{
tsScoredCol* pEl1 = (tsScoredCol*)arg1;
tsScoredCol* pEl2 = (tsScoredCol*)arg2;

if(pEl1->ColIdx < pEl2->ColIdx)
	return(-1);
if(pEl1->ColIdx > pEl2->ColIdx)
	return(1);
return(0);
}


