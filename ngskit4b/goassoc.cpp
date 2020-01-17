// goassoc.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
/*
This toolkit is a source base clone of 'BioKanga' release 4.4.2 (https://github.com/csiro-crop-informatics/biokanga) and contains
significant source code changes enabling new functionality and resulting process parameterisation changes. These changes have resulted in
incompatibilty with 'BioKanga'.

Because of the potentential for confusion by users unaware of functionality and process parameterisation changes then the modified source base
and resultant compiled executables have been renamed to 'kit4b' - K-mer Informed Toolkit for Bioinformatics.
The renaming will force users of the 'BioKanga' toolkit to examine scripting which is dependent on existing 'BioKanga'
parameterisations so as to make appropriate changes if wishing to utilise 'kit4b' parameterisations and functionality.

'kit4b' is being released under the Opensource Software License Agreement (GPLv3)
'kit4b' is Copyright (c) 2019
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

#ifndef _WIN32
#include <config.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const int cMaxGOterms = 100000;		// assume may need to handle this many GOterms when generating GraphViz dot file
const int cElGenes2Alloc = 50000;	// allocate for unordered gene lists in this many increments
const int cAllocUniqueGenes = 2000;  // allocate for ordered unique gene lists in this many increments

const double cMinDotThres  = 0.000001; // min threshold at which terms with P values above this are filtered out from GraphViz dot file
const double cDfltDotThres = 0.01;	// default threshold at which terms with P values above this are filtered out from GraphViz dot file
const double cMaxDotThres  = 0.5;	// max threshold at which terms with P values above this are filtered out from GraphViz dot file
 
// processing modes
typedef enum TAG_ePMode {
	ePMdefault,					// default is for hits CSV file
	ePMDESeqPVal,				// process DESeq (use PVal) generated file
	ePMDESeqPValAdj,			// process DESeq (use PValAdj) generated file
	ePMplaceholder				// used to set the enumeration range
	} etPMode;

// which P-Value to use when generating GraphViz dot output
typedef enum eDotPValue {
	eDPVunweighted = 0,				// unweighted and no multiple term corection
	eDPVweighted,					// weighted but with no multiple term correction
	eDPVunweightedMTC,				// unweighted and multiple term corection
	eDPVweightedNTC					// weighted and multiple term correction
	} etDotPValue;

typedef enum eMTCMode {
	ePMNoCorrection = 0,		// no multiple corrections
	ePMBonferroni,				// Bonferroni correction
	ePMHolm,					// Step-down 
	ePMBenjaminiHochberg		// Benjamini and Hochberg FDR
	} etMTCMode;

typedef enum eTestMeth {
	eTMFisher = 0,				// default is to use fisher exact
	eTMChiSquare				// optionally could use chi-square
	} etTestMeth;


typedef struct TAG_sElementGene {
	int Weighting;						// weighting to associate with this element
	int Length;							//length of this element						
	char szGeneName[cMaxFeatNameLen+1]; // element or gene name
} tsElementGene;

typedef struct TAG_sTermStats {
	int TotNumGenes;				// total number of genes, some may not be associated with any GO term
	int SumTotGeneWeightings;		// sum of all gene weightings, some may not be associated with any GO term		
	int NumGenesNotInBED;			// number of genes not located in BED 
	int SumGenesNotInBEDWeightings;// sum of all genes not in BED weightings	
	int NumGenesNotOnStrand;		// number of genes not on requested strand
	int SumGenesNotOnStrandWeightings;  // sum of all genes not on requested strand weightings		
	int NumGenesNotInGO;			// number of genes with no GO:Term association
	int SumGenesNotInGOWeightings;	// sum of all genes not in GO:Term association weightings		
	int NumGenesInGO;				// number of genes which were associated with GO:Terms
	int SumNumGenesInGOWeightings;	// sum of all gene weightings which were associated with GO:Terms		
} tsTermStats;


char *MultipleTestingMode2Txt(int Mode);
char *Ontologies2Txt(etOntologies Ontology);
char *TestMeth2Txt(etTestMeth TestMeth);

int
Process(etPMode PMode,			// processing mode
		double DESeqThres,		// cutoff PVal when processing DESeq generated file
		etMTCMode MTC,			// multiple test correction mode
		bool bCanonicalise,     // canonicalise sample name isoforms by removing any numerical suffix in range '.[0-99] (TAIR genes)
		etTestMeth TestMeth,	// hypothesis test method
		etOntologies Ontology,  // which root ontology to process
		bool bMultGeneHits,		// allow multiple gene hits by sample genes 
		int MinLen,				// process sample elements of >= this minimim length
		int MaxLen,				// process sample elements of <= this maximum length
		int MinWeight,			// process sample elements of >= this minimim weight
		int MaxWeight,			// process sample elements of <= this maximum weight
		bool bProp,				// propagate counts from GO:Term into parent terms
		bool bBkgndLen,			// background counts proportional to gene lengths
		char Strand,			// background counts are for which strand '*', '+' or '-'
		char *pszBED,			// gene BED file
		char *pszGoAssoc,		// gene to ID association file
		char *pszGOTerms,		// GO ontology file
		char *pszPopulationGenes, // file containing background population genes
		char *pszHitsFile,		// file containing sample hits
		char *pszResultsFile,	// results file to generate
		char *pszGraphVizDot,	// GraphViz dot to generate
		char *pszRootGOTerm,	// optionally treat this term as the root term when generating GraphViz
		double dblDotThres,		// filter out terms from GraphViz dot file with higher P-values
		etDotPValue DotPValue);	// which P-Value to use when generating GraphViz dot file

int
SetBkgndCnts(etMTCMode MTC,		// multiple test correction mode
		bool bCanonicalise,     // canonicalise sample name isoforms by removing any numerical suffix in range '.[0-99] (TAIR ge
 		etOntologies Ontology, // which root ontology to process
		bool bProp,				// propagate counts from GO:Term into parent terms
		bool bBkgndLen,			// background counts proportional to gene lengths
		char Strand,			// background counts are for which strand '*', '+' or '-'
		int PopulationCnt,		// number of genes in background population 
		tsElementGene *pPopulationGenes, // genes to treat as being the background population
		CBEDfile *pBED,			// gene BED file
		CGOAssocs *pAssocs,		// gene to ID association file
		CGOTerms *pGOTerms,		// GO ontology file
		tsTermStats *pStats);	// returned processing stats


int
SetSampleCnts(etMTCMode MTC,	// multiple test correction mode
	    bool bCanonicalise,     // canonicalise sample name isoforms by removing any numerical suffix in range '.[0-99] (TAIR ge
  		etOntologies Ontology, // which root ontology to process
		bool bProp,				// propagate counts from GO:Term into parent terms
		char Strand,			// background counts are for which strand '*', '+' or '-'
		int SampleCnt,			// number of genes in sample 
		tsElementGene *pSampleGenes, // genes + counts to treat as being the sample
		CBEDfile *pBED,			// gene BED file
		CGOAssocs *pAssocs,		// gene to ID association file
		CGOTerms *pGOTerms,		// GO ontology file
		tsTermStats *pStats);	// returned processing stats

int GenPValuesChiSquare(etOntologies Ontology,CGOTerms *pGOTerms,tsTermStats *pPopulationStats,tsTermStats *pSampleStats);
int GenPValuesFishersExactTest(etOntologies Ontology,CGOTerms *pGOTerms,tsTermStats *pPopulationStats,tsTermStats *pSampleStats);

int MTCNone(etOntologies Ontology,CGOTerms *pGOTerms);
int MTCBonferroni(etOntologies Ontology,CGOTerms *pGOTerms);
int MTCHolm(etOntologies Ontology,CGOTerms *pGOTerms);
int MTCBenjaminiHochberg(etOntologies Ontology,CGOTerms *pGOTerms);

int
ParseCSVHits(char *pszHitsFile,			// file to parse
			CBEDfile* pBED,				// hits must be to genes in this file
		  int MinLen,					// filter out elements less than this length
		  int MaxLen,					// filter out elements more than this length
		  int MinWeight,				// filter out elements less than this number of hits or weighting
		  int MaxWeight,				// filter out elements more than this number of hits or weighting
		  int *pNumEls,					// where to return number of elements in returned list				
		  tsElementGene **ppRetElements);	// returned list of elements

int
ParseDESeq(etPMode PMode,				// procesisng mode
		   double Thres,				// only process for PVals <= this threshold
		   char *pszDESeqFile,			// file to parse
		  int *pNumEls,					// where to return number of elements in returned list				
		  tsElementGene **ppRetElements);	// returned list of elements

char *DotPValueTxt(etDotPValue DotPValue);
int ReduceGeneList(bool bKeepDups,int NumGenes,tsElementGene *pGeneList,int *pNumUniqueGenes,tsElementGene **ppUniqueGeneList);

int
OutputResults(etOntologies Ontology,
			  tsTermStats *pPopulationStats,
			  tsTermStats *pSampleStats,
			  char *pszResultsFile,
			  char *pszGraphViz,
			  CGOTerms *pGOTerms,
			  char *pszRootGOTerm,
			  double dblDotThres,
			  etDotPValue DotPValue);

static int SortGenes(const void *arg1, const void *arg2);

#ifdef _WIN32
int goassoc(int argc, char* argv[])
{
	// determine my process name
	_splitpath(argv[0], NULL, NULL, gszProcName, NULL);
#else
int
goassoc(int argc, char** argv)
{
	// determine my process name
	CUtility::splitpath((char*)argv[0], NULL, gszProcName);
#endif
int iScreenLogLevel;		// level of screen diagnostics
int iFileLogLevel;			// level of file diagnostics
char szLogFile[_MAX_PATH];	// write diagnostics to this file

int Rslt;

etPMode PMode;					// processing mode
etOntologies iOntology;			// ontologies to process: 1 Cellular, 2 Biological, 3 Molecular
double DESeqThres;				// only process DESeq PVal or PValAdj <= this threshold
int iMTC;						// Multiple Test Correction mode
int iTestMeth;					// test method
bool bMultGeneHits;				// allow multiple sample hits onto population gene
int iMinLen;					// process sample elements of >= this minimim length
int iMaxLen;					// process sample elements of <= this maximum length

int iMinWeight;					// process sample elements of >= this minimim weight
int iMaxWeight;					// process sample elements of <= this maximum weight

bool bProp;						// propagate counts from GO:Term into parent terms
bool bBkgndLen;					// background counts proportional to gene lengths
int  iBkgndStrand;				// background strand 0==both, 1=='+', 2=='-'
char cStrand;					// background strand

char szBED[_MAX_PATH];			// gene BED file
char szGoAssoc[_MAX_PATH];		// gene to ID association file
char szGOTerms[_MAX_PATH];		// GO ontology file
char szElGenes[_MAX_PATH];		// background population genes file
char szHitsFile[_MAX_PATH];		// file containing hits
char szResultsFile[_MAX_PATH];	// results file to generate
char szGraphVizDot[_MAX_PATH];	// GraphViz dot to generate
char szRootGOTerm[100];			// root term when generating GraphViz dot output
double dblDotThres;				// filter out any probs > than this when generating GraphViz dot file
int iDotPValue;					// which P-Value to use when generating GraphViz dot file
bool bCanonicalise;				// canonicalise sample name isoforms by removing any numerical suffix in range '.[0-99] (TAIR genes)

char szSQLiteDatabase[_MAX_PATH];	// results summaries to this SQLite file
char szExperimentName[cMaxDatasetSpeciesChrom + 1];			// experiment name
char szExperimentDescr[1000];		// describes experiment

// command line args
struct arg_lit  *help    = arg_lit0("hH","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *pmode = arg_int0("m","mode","<int>",		    "processing mode: 0 - hits, 1 - DESeq PVal, 2 - DESeq PValAdj (default = 0)");
struct arg_dbl *deseqthres = arg_dbl0("d","deseqthres","<dbl>",	"only process DESeq PVal or PValAdj <= this threshold (default = 0.05");

struct arg_lit  *canonicalise = arg_lit0("c","canonicalise",	"canonicalise sample name isoforms by removing any numerical suffix in range '.[0-99] (TAIR genes)");

struct arg_int  *Ontology = arg_int0("r","ontology","<int>",	"Ontologies: (default)1:Cellular 2:Biological 3:Molecular");
struct arg_int  *TestMeth = arg_int0("t","test","<int>",		"Test method (default) 0:Chi-square 1:Fisher exact 2:Hypergeometric");

struct arg_int  *mtc = arg_int0("x","MTC","<int>",				"MTC (default)0:None 1:Bonferroni 2:Holm 3:Benjamini and Hochberg FDR");
struct arg_lit  *MultGeneHits = arg_lit0("X","multhits",		"allow multiple sample hits onto population gene");

struct arg_int  *MinLen = arg_int0("n","minlen","<int>",		"process sample elements of >= this minimim length");
struct arg_int  *MaxLen = arg_int0("N","maxlen","<int>",		"process sample elements of <= this maximum length");

struct arg_int  *MinWeight = arg_int0("w","minweight","<int>",	"process sample elements of >= this minimim weighting");
struct arg_int  *MaxWeight = arg_int0("W","maxweight","<int>",	"process sample elements of <= this maximum weighting");

struct arg_lit  *Prop = arg_lit0("p","propagate",				"propagate counts from GO:Terms into their parent terms");
struct arg_lit  *BkgndLen = arg_lit0("l","bkgndlen",			"background counts are proportional to gene lengths");
struct arg_int  *BkgndStrand = arg_int0("s","strand","<int>",	"background strand default 0==both,1=='+',2=='-'");

struct arg_file *BED = arg_file1("b","bedgenefile","<file>",		"input bioBED gene file");
struct arg_file *GoAssoc = arg_file1("i","goassocfile","<file>",			"input GO associations file");
struct arg_file *GOTerms = arg_file1("I","gotermsfile","<file>",			"input GO terms file");
struct arg_file *ElGenes = arg_file0("P","popgenesfile","<file>",			"input population genes .csv file");
struct arg_file *HitsFile = arg_file1("g","samplehitsfile","<file>",		"input sample file - gene+weighting+len (hits) .csv file");
struct arg_file *ResultsFile = arg_file1("o","rsltsfile","<file>",		"output results file");
struct arg_file *GraphVizDot = arg_file0("O","graphviz","<file>",		"output into GraphViz dot file");
struct arg_str  *RootGOTerm = arg_str0("T","dotrootterm","<string>",		"treat this term as the root when generating GraphViz dot file");
struct arg_dbl  *DotThres = arg_dbl0("j","dotprobthres","<dbl>",	"dot - only for probs of less or equal this threshold");
struct arg_int  *DotPValue = arg_int0("J","dotpvalue","<int>",		"dot - 0: unweighted 1: weighted 2: unweighted + MTC 3: weighted + MTC");

struct arg_file* summrslts = arg_file0("q", "sumrslts", "<file>", "Output results summary to this SQLite3 database file");
struct arg_str* experimentname = arg_str0("w", "experimentname", "<str>", "experiment name SQLite3 database file");
struct arg_str* experimentdescr = arg_str0("W", "experimentdescr", "<str>", "experiment description SQLite3 database file");

struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					summrslts,experimentname,experimentdescr,
					pmode,deseqthres,canonicalise,Ontology,mtc,TestMeth,MultGeneHits,MinLen,MaxLen,MinWeight,MaxWeight,Prop,BkgndLen,BkgndStrand,BED,GoAssoc,
					GOTerms,ElGenes,HitsFile,ResultsFile,GraphVizDot,RootGOTerm,DotThres,DotPValue,
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

	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Subprocess %s Version %s starting", gpszSubProcess->pszName, kit4bversion);
	gExperimentID = 0;
	gProcessID = 0;
	gProcessingID = 0;
	szSQLiteDatabase[0] = '\0';
	szExperimentName[0] = '\0';
	szExperimentDescr[0] = '\0';

	if (experimentname->count)
	{
		strncpy(szExperimentName, experimentname->sval[0], sizeof(szExperimentName));
		szExperimentName[sizeof(szExperimentName) - 1] = '\0';
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
	if (summrslts->count)
	{
		strncpy(szSQLiteDatabase, summrslts->filename[0], sizeof(szSQLiteDatabase) - 1);
		szSQLiteDatabase[sizeof(szSQLiteDatabase) - 1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(szSQLiteDatabase);
		if (strlen(szSQLiteDatabase) < 1)
		{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: After removal of whitespace, no SQLite database specified with '-q<filespec>' option");
			return(1);
		}

		if (strlen(szExperimentName) < 1)
		{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: After removal of whitespace, no SQLite experiment name specified with '-w<str>' option");
			return(1);
		}
		if (experimentdescr->count)
		{
			strncpy(szExperimentDescr, experimentdescr->sval[0], sizeof(szExperimentDescr) - 1);
			szExperimentDescr[sizeof(szExperimentDescr) - 1] = '\0';
			CUtility::TrimQuotedWhitespcExtd(szExperimentDescr);
		}
		if (strlen(szExperimentDescr) < 1)
		{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: After removal of whitespace, no SQLite experiment description specified with '-W<str>' option");
			return(1);
		}

		gExperimentID = gSQLiteSummaries.StartExperiment(szSQLiteDatabase, false, true, szExperimentName, szExperimentName, szExperimentDescr);
		if (gExperimentID < 1)
			return(1);
		gProcessID = gSQLiteSummaries.AddProcess((char*)gpszSubProcess->pszName, (char*)gpszSubProcess->pszName, (char*)gpszSubProcess->pszFullDescr);
		if (gProcessID < 1)
			return(1);
		gProcessingID = gSQLiteSummaries.StartProcessing(gExperimentID, gProcessID, (char*)kit4bversion);
		if (gProcessingID < 1)
			return(1);
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Initialised SQLite database '%s' for results summary collection", szSQLiteDatabase);
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "SQLite database experiment identifier for '%s' is %d", szExperimentName, gExperimentID);
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "SQLite database process identifier for '%s' is %d", (char*)gpszSubProcess->pszName, gProcessID);
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "SQLite database processing instance identifier is %d", gProcessingID);
	}
	else
	{
		szSQLiteDatabase[0] = '\0';
		szExperimentDescr[0] = '\0';
	}


	PMode = (etPMode)(pmode->count ? pmode->ival[0] : ePMdefault);
	if(PMode < ePMdefault || PMode >= ePMplaceholder)
		{
		printf("\nError: Processing mode '-m%d' specified outside of range %d..%d",PMode,0,(int)ePMplaceholder-1);
		exit(1);
		}

	if(PMode != ePMdefault)
		{
		DESeqThres = deseqthres->count ? deseqthres->dval[0] : 0.05f;
		if(DESeqThres < 0.0001f || DESeqThres > 0.25f)
			{
			printf("\nError: DESeq threshold '-d%f' outside of range 0.0001 to 0.25",DESeqThres);
			exit(1);
			}
		}
	else
		DESeqThres = 0.0f;

	bCanonicalise = canonicalise->count ? true : false;

	iOntology = (etOntologies)(Ontology->count ? Ontology->ival[0] : eONTCellular);
	if(iOntology < eONTCellular || iOntology > eONTMolecular)
		{
		printf("\nError: Requested ontology '-O%d' not supported",iOntology);
		exit(1);
		}

	iMTC = mtc->count ? mtc->ival[0] : ePMNoCorrection;
	if(iMTC < ePMNoCorrection || iMTC > ePMBenjaminiHochberg)
		{
		printf("\nError: Requested Multiple Testing Correction Mode '-x%d' not supported",iMTC);
		exit(1);
		}

	iTestMeth = TestMeth->count ? TestMeth->ival[0] : eTMFisher;
	if(iTestMeth < eTMFisher || iTestMeth > eTMChiSquare)
		{
		printf("\nError: Requested hypothesis test method '-t%d' not supported",iTestMeth);
		exit(1);
		}

	bMultGeneHits = MultGeneHits->count ? true : false;

	iMinLen = MinLen->count ? MinLen->ival[0] : 1;
	if(iMinLen < 1)
		{
		printf("\nError: Requested minimum sample element length '-n%d' not supported",iMinLen);
		exit(1);
		}
	iMaxLen = MaxLen->count ? MaxLen->ival[0] : 1000000000;
	if(iMaxLen < 1)
		{
		printf("\nError: Requested max sample element length '-N%d' not supported",iMaxLen);
		exit(1);
		}
	if(iMaxLen < iMinLen)
		{
		printf("\nError: Requested max sample element length '-N%d' less than minimum '-n%d' requested",iMaxLen,iMinLen);
		exit(1);
		}

	iMinWeight = MinWeight->count ? MinWeight->ival[0] : 1;
	if(iMinWeight < 1)
		{
		printf("\nError: Requested minimum sample weight '-w%d' not supported",iMinWeight);
		exit(1);
		}
	iMaxWeight = MaxWeight->count ? MaxWeight->ival[0] : 1000000000;
	if(iMaxWeight < 1)
		{
		printf("\nError: Requested max sample weight '-W%d' not supported",iMaxWeight);
		exit(1);
		}
	if(iMaxWeight < iMinWeight)
		{
		printf("\nError: Requested max sample weight '-W%d' less than minimum '-w%d' requested",iMaxWeight,iMinWeight);
		exit(1);
		}

	iBkgndStrand = BkgndStrand->count ? BkgndStrand->ival[0] : 0;
	if(iBkgndStrand < 0 || iBkgndStrand > 2)
		{
		printf("\nError: Requested background strand '-s%d' not supported",iBkgndStrand);
		exit(1);
		}
	switch(iBkgndStrand) {
		case 0:
			cStrand = '*';
			break;
		case 1:
			cStrand = '+';
			break;
		case 2:
			cStrand = '-';
			break;
		}

	bProp = Prop->count ? true : false;
	bBkgndLen = BkgndLen->count ? true : false;

	dblDotThres = DotThres->count ? DotThres->dval[0] : cDfltDotThres;
	if(dblDotThres < cMinDotThres)
		{
		printf("\nError: Requested dot threshold must be at least %f",cMinDotThres);
		exit(1);
		}
	else
		if(dblDotThres >= cMaxDotThres)
				{
				printf("\nError: Requested dot threshold must be at less than %f",cMaxDotThres);
				exit(1);
				}

	iDotPValue = DotPValue->count ? DotPValue->ival[0] : eDPVunweightedMTC;
	if(iDotPValue < eDPVunweighted || iDotPValue > eDPVweightedNTC)
		{
		printf("\nError: Requested dot P-Value not supported");
		exit(1);
		}

	if(RootGOTerm->count)
		strcpy(szRootGOTerm,RootGOTerm->sval[0]);
	else
		szRootGOTerm[0] = '\0';

	if(ElGenes->count)
		strcpy(szElGenes,ElGenes->filename[0]);
	else
		szElGenes[0] = '\0';
	strcpy(szBED,BED->filename[0]);
	strcpy(szGoAssoc,GoAssoc->filename[0]);
	strcpy(szGOTerms,GOTerms->filename[0]);
	strcpy(szHitsFile,HitsFile->filename[0]);
	strcpy(szResultsFile,ResultsFile->filename[0]);
	if(GraphVizDot->count)
		strcpy(szGraphVizDot,GraphVizDot->filename[0]);
	else
		szGraphVizDot[0] = '\0';

			// now that command parameters have been parsed then initialise diagnostics log system
	if(!gDiagnostics.Open(szLogFile,(etDiagLevel)iScreenLogLevel,(etDiagLevel)iFileLogLevel,true))
		{
		printf("\nError: Unable to start diagnostics subsystem.");
		if(szLogFile[0] != '\0')
			printf(" Most likely cause is that logfile '%s' can't be opened/created",szLogFile);
		exit(1);
		}

// show user current resource limits
#ifndef _WIN32
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Resources: %s",CUtility::ReportResourceLimits());
#endif

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %s Processing parameters:",kit4bversion);

	const char *pszDescr;
	switch(PMode) {
		case ePMdefault:
			pszDescr = "process CSV hits file";
			break;
		case ePMDESeqPVal:
			pszDescr = "process DESeq (use PVal) generated file";
			break;
		case ePMDESeqPValAdj:
			pszDescr = "process DESeq (use PValAdj) generated file";
			break;
		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Processing mode is : '%s'",pszDescr);

	if(PMode != ePMdefault)
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"DESeq threshold: %f",DESeqThres);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Ontology: %s",Ontologies2Txt(iOntology));
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Canonicalise gene names: %s",bCanonicalise ? "Yes" : "No");
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Hypothesis test method: %s",TestMeth2Txt((etTestMeth)iTestMeth));
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Multiple Test Correction Mode: %s",MultipleTestingMode2Txt(iMTC));
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Allow multiple sample hits per population gene: %s",bMultGeneHits ? "Yes":"No");
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Propagate counts from child GO:Terms into parent GO:Terms: %s",bProp ? "Yes" : "No");

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum length sample hits: %d",iMinLen);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Maximum length sample hits: %d",iMaxLen);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Minimum sample weight: %d",iMinWeight);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Maximum sample weight: %d",iMaxWeight);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"background counts proportional to gene lengths: %s",bBkgndLen ? "Yes" : "No");
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"background counts on: %c",cStrand);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"gene BED population file: '%s'",szBED);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"gene to ID association file: '%s'",szGoAssoc);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"GO ontology file: '%s'",szGOTerms);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"file containing population genes: '%s'",szElGenes[0] == '\0' ? "NONE" : szElGenes);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"file containing hits: '%s'",szHitsFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"results file to generate: '%s'",szResultsFile);
	if(szGraphVizDot[0] != '\0')
		{
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"GraphViz dot file to generate: '%s'",szGraphVizDot);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Filter out terms with P-values more than: '%f'",dblDotThres);
		gDiagnostics.DiagOutMsgOnly(eDLInfo,"Use this P-Value in generated dot file: '%s'",DotPValueTxt((etDotPValue)iDotPValue));
		if(szRootGOTerm[0] != '\0')
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Use this term as the root: '%s'",szRootGOTerm);
		}

	if (gExperimentID > 0)
	{
		int ParamID;
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID, ePTText, (int)strlen(szLogFile), "log", szLogFile);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID, ePTInt32, sizeof(PMode), "pmode", &PMode);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID, ePTDouble, sizeof(DESeqThres), "deseqthres", &DESeqThres);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID, ePTBool, sizeof(bCanonicalise), "canonicalise", &bCanonicalise);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID, ePTInt32, sizeof(iOntology), "ontology", &iOntology);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID, ePTInt32, sizeof(iTestMeth), "test", &iTestMeth);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID, ePTInt32, sizeof(iMTC), "mtc", &iMTC);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID, ePTBool, sizeof(bMultGeneHits), "multhits", &bMultGeneHits);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID, ePTBool, sizeof(bProp), "propagate", &bProp);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID, ePTInt32, sizeof(iMinLen), "minlen", &iMinLen);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID, ePTInt32, sizeof(iMaxLen), "maxlen", &iMaxLen);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID, ePTInt32, sizeof(iMinWeight), "minweight", &iMinWeight);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID, ePTInt32, sizeof(iMaxWeight), "maxweight", &iMaxWeight);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID, ePTBool, sizeof(bBkgndLen), "bkgndlen", &bBkgndLen);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID, ePTInt32, sizeof(cStrand), "strand", &cStrand);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID, ePTText, (int)strlen(szBED), "bedgenefile", szBED);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID, ePTText, (int)strlen(szGoAssoc), "goassocfile", szGoAssoc);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID, ePTText, (int)strlen(szGOTerms), "gotermsfile", szGOTerms);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID, ePTText, (int)strlen(szElGenes), "popgenesfile", szElGenes);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID, ePTText, (int)strlen(szHitsFile), "samplehitsfile", szHitsFile);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID, ePTText, (int)strlen(szResultsFile), "szResultsFile", szResultsFile);

		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID, ePTText, (int)strlen(szGraphVizDot), "graphviz", szGraphVizDot);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID, ePTText, (int)strlen(szRootGOTerm), "dotrootterm", szRootGOTerm);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID, ePTDouble, sizeof(dblDotThres), "dotprobthres", &dblDotThres);
		ParamID = gSQLiteSummaries.AddParameter(gExperimentID, gProcessingID, ePTInt32, sizeof(iDotPValue), "dotpvalue", &iDotPValue);
	}
	gStopWatch.Start();
#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	Rslt = Process(PMode,DESeqThres,(etMTCMode)iMTC,bCanonicalise,(etTestMeth)iTestMeth,iOntology,bMultGeneHits,iMinLen,iMaxLen,iMinWeight,iMaxWeight,bProp,bBkgndLen,cStrand,szBED,szGoAssoc,szGOTerms,szElGenes,szHitsFile,szResultsFile,szGraphVizDot,szRootGOTerm,dblDotThres,(etDotPValue)iDotPValue);
	gStopWatch.Stop();
	if (gExperimentID > 0)
		{
		if (gProcessingID)
			gSQLiteSummaries.EndProcessing(gExperimentID, gProcessingID, Rslt);
		gSQLiteSummaries.EndExperiment(gExperimentID);
		}
	Rslt = Rslt < 0 ? 1 : 0;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit Code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	exit(Rslt);
	}
else
	{
	printf("\n%s Generate GO Analysis, Version %s\n",gszProcName,kit4bversion);
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nUse '-h' to view option and parameter usage\n");
	exit(1);
	}
}

char *
TestMeth2Txt(etTestMeth TestMeth)
{
switch(TestMeth) {
	case eTMFisher:			// default is to use fisher exact
		return((char *)"Fisher exact");
	case eTMChiSquare:		// or optionally chi-square
		return((char *)"Chi-square");
	default:
		break;
	}
return((char *)"Unsupported");
}


// which P-Value to use when generating GraphViz dot output
char *
DotPValueTxt(etDotPValue DotPValue)
{
switch(DotPValue) {
	case eDPVunweighted:				// unweighted and no multiple term corection
		return((char *)"Unweighted and no multiple term corection");
	case eDPVweighted:					// weighted but with no multiple term correction
		return((char *)"Weighted and no multiple term correction");
	case eDPVunweightedMTC:				// unweighted and multiple term corection
		return((char *)"Unweighted and but with multiple term corection");
	case eDPVweightedNTC:					// weighted and multiple term correction
		return((char *)"Weighted and and with multiple term correction");
	default:
		break;
	}
return((char *)"Unsupported");
}

char *
MultipleTestingMode2Txt(int Mode)
{
char *pszMode;
switch(Mode) {
	case ePMNoCorrection: 
		pszMode = (char *)"no multiple corrections"; 
		break;
	case ePMBonferroni: 
		pszMode = (char *)"Bonferroni correction"; 
		break;
	case ePMHolm: 
		pszMode = (char *)"Holm Step-down";  
		break;
	case ePMBenjaminiHochberg: 
		pszMode=(char *)"Benjamini and Hochberg FDR"; 
		break;
	default: 
		pszMode=(char *)"Unknown"; 
		break;
	}
return(pszMode);
}

char *
Ontologies2Txt(etOntologies Ontology)
{
char *pszOntology;
switch(Ontology) {
	case eONTCellular:
		pszOntology = (char *)"Cellular component";
		break;
	case eONTBiological:
		pszOntology = (char *)"Biological process";
		break;
	case eONTMolecular:
		pszOntology =(char *) "Molecular function";
		break;
	default:
		pszOntology =(char *) "Unrecognised";
		break;
	}
return(pszOntology);
}

int
Process(etPMode PMode,			// processing mode
		double DESeqThres,		// cutoff PVal when processing DESeq generated file
		etMTCMode MTC,			// multiple test correction mode
		bool bCanonicalise,     // canonicalise sample name isoforms by removing any numerical suffix in range '.[0-99] (TAIR genes)
		etTestMeth TestMeth,	// hypothesis test method
		etOntologies Ontology, // which root ontology to process
		bool bMultGeneHits,		// allow multiple gene hits by sample genes 
		int MinLen,				// process sample elements of >= this minimim length
		int MaxLen,				// process sample elements of <= this maximum length
		int MinWeight,			// process sample elements of >= this minimim weight
		int MaxWeight,			// process sample elements of <= this maximum weight
		bool bProp,				// propagate counts from GO:Term into parent terms
		bool bBkgndLen,			// background counts proportional to gene lengths
		char Strand,			// background counts are for which strand '*', '+' or '-'
		char *pszBED,			// gene BED file
		char *pszGoAssoc,		// gene to ID association file
		char *pszGOTerms,		// GO ontology file
		char *pszPopulationGenes, // file containing background population genes
		char *pszHitsFile,		// file containing sample genes
		char *pszResultsFile,	// results file to generate
		char *pszGraphVizDot,	// GraphViz dot to generate
		char *pszRootGOTerm,	// treat this term as the root term
		double dblDotThres,		// filter out P-values > DotThres from generated .dot file
		etDotPValue DotPValue)	// which P-Value to use when generating GraphViz dot file
{
int Rslt;
CGOAssocs *pAssocs = NULL;
CGOTerms *pGOTerms = NULL;
CBEDfile *pBED = NULL;
tsElementGene *pPopulationGenes = NULL;
int AllocdPopulationGenes = 0;
int PopulationGeneCnt = 0;
tsElementGene *pSampleGenes = NULL;
int AllocdSampleGenes = 0;
int SampleGeneCnt = 0;

int NumFields = 0;
int GeneWeighting = 0;
int GeneLength = 0;
char *pszGene = NULL;


if(pszGoAssoc == NULL || *pszGoAssoc == '\0' ||
   pszGOTerms == NULL || *pszGOTerms == '\0' ||
   pszBED == NULL || *pszBED == '\0' ||
   pszHitsFile == NULL || *pszHitsFile == '\0' ||
	pszResultsFile == NULL || *pszResultsFile == '\0')
	return(eBSFerrParams);


pBED = new CBEDfile();
if((Rslt = pBED->Open(pszBED))!=eBSFSuccess)
	{
	while(pBED->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pBED->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open BED gene file '%s'",pszBED);
	delete pBED;
	return(Rslt);
	}

pAssocs = new CGOAssocs();
if((Rslt = pAssocs->Open(pszGoAssoc))!=eBSFSuccess)
	{
	while(pAssocs->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pAssocs->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open gene to GO:ID association file '%s'",pszGoAssoc);
	delete pAssocs;
	delete pBED;
	return(Rslt);
	}

pGOTerms = new CGOTerms();
if((Rslt = pGOTerms->Open(pszGOTerms))!=eBSFSuccess)
	{
	while(pGOTerms->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pGOTerms->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open GO:Term ontology file '%s'",pszGOTerms);
	delete pAssocs;
	delete pBED;
	delete pGOTerms;
	return(Rslt);
	}

pPopulationGenes = NULL;
if(pszPopulationGenes != NULL && pszPopulationGenes[0] != '\0')
	if((Rslt = ParseCSVHits(pszPopulationGenes, pBED, MinLen, MaxLen, MinWeight, MaxWeight,&PopulationGeneCnt, &pPopulationGenes)) < 0)
		{
		if(pPopulationGenes != NULL)
			{
			delete pPopulationGenes;
			pPopulationGenes = NULL;
			}
		delete pAssocs;
		delete pBED;
		delete pGOTerms;
		return(Rslt);
		}

switch(PMode) {
	case ePMdefault:		// default is for hits CSV file
		Rslt = ParseCSVHits(pszHitsFile,pBED,MinLen,MaxLen,MinWeight,MaxWeight,&SampleGeneCnt,&pSampleGenes);
		break;
	default:				// process DESeq generated file
		Rslt = ParseDESeq(PMode,DESeqThres,pszHitsFile,&SampleGeneCnt,&pSampleGenes);
		break;
	}

if(!SampleGeneCnt || Rslt < eBSFSuccess)
	{
	if(!SampleGeneCnt && Rslt >= eBSFSuccess)
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Nothing to do - no samples meeting filtering criteria");

	if(pSampleGenes != NULL)
		delete pSampleGenes;
	if (pPopulationGenes != NULL)
		delete pPopulationGenes;
	delete pAssocs;
	delete pBED;
	delete pGOTerms;
	return(Rslt);
	}

tsTermStats PopulationStats;
tsTermStats SampleStats;
int NumUniqueElGenes = 0;
tsElementGene *pUniquePopGeneList = NULL;
int NumUniqueSampleGenes = 0;
tsElementGene *pUniqueSampleGeneList = NULL;

if(PopulationGeneCnt && pPopulationGenes != NULL)
	Rslt = ReduceGeneList(false,PopulationGeneCnt,pPopulationGenes,&NumUniqueElGenes,&pUniquePopGeneList);
if(Rslt == eBSFSuccess)
   Rslt = ReduceGeneList(bMultGeneHits,SampleGeneCnt,pSampleGenes,&NumUniqueSampleGenes,&pUniqueSampleGeneList);

if(Rslt == eBSFSuccess)
	Rslt = SetBkgndCnts(MTC,bCanonicalise,Ontology,bProp,bBkgndLen,Strand,NumUniqueElGenes,pUniquePopGeneList,pBED,pAssocs,pGOTerms,&PopulationStats);
if(Rslt >= eBSFSuccess)
	Rslt = SetSampleCnts(MTC,bCanonicalise,Ontology,bProp,Strand,NumUniqueSampleGenes,pUniqueSampleGeneList,pBED,pAssocs,pGOTerms,&SampleStats);

pBED->Close();
delete pBED;
pAssocs->Close();
delete pAssocs;

// can now generate the unadjusted (not corrected for multiple testing) P values
if(Rslt >= eBSFSuccess)
	{
	switch(TestMeth) {
		case eTMFisher:			// default is to use fisher exact 
			Rslt = GenPValuesFishersExactTest(Ontology,pGOTerms,&PopulationStats,&SampleStats);
			break;
		case eTMChiSquare:		// optionly to use chi-square
			Rslt = GenPValuesChiSquare(Ontology,pGOTerms,&PopulationStats,&SampleStats);
			break;
		}
	}

if(Rslt >= eBSFSuccess)
	switch(MTC) {
		case ePMNoCorrection:		// no multiple test corrections
			MTCNone(Ontology,pGOTerms);
			break;
		case ePMBonferroni:			// Bonferroni correction
			MTCBonferroni(Ontology,pGOTerms);
			break;

		case ePMHolm:				// Step-down
			MTCHolm(Ontology,pGOTerms);
			break;

		case ePMBenjaminiHochberg:	// Benjamini and Hochberg FDR
			MTCBenjaminiHochberg(Ontology,pGOTerms);
			break;
		};

// now time to generate results...
if(Rslt >= eBSFSuccess)
	OutputResults(Ontology,&PopulationStats,&SampleStats,
					pszResultsFile,pszGraphVizDot,pGOTerms,pszRootGOTerm,dblDotThres,DotPValue);
	
pGOTerms->Close();
delete pGOTerms;
if(pUniquePopGeneList != NULL)
	delete pUniquePopGeneList;
if(pUniqueSampleGeneList!=NULL)
	delete pUniqueSampleGeneList;
return(Rslt);
}

// ParseDESeq
// hits file expected to contain csv delimited output as generated by DESeq
// The gene names are in column 2, PVal in column 8, and PValAdj in column 9
// First row contains column titles
int
ParseDESeq(etPMode PMode,				// procesisng mode
		   double Thres,				// only process for PVals <= this threshold
		   char *pszDESeqFile,			// file to parse
		  int *pNumEls,					// where to return number of elements in returned list				
		  tsElementGene **ppRetElements)// returned list of elements
{
int Rslt;
int NumFields;
int LineNum;
char *pszPVal;
char *pszGene;
double PVal;
char *pszTerm;
CCSVFile *pCSVDESeq;
int SampleGeneCnt;
int AllocdSampleGenes;
tsElementGene *pSampleGenes;
tsElementGene *pTmp;
int AboveThresCnt;
int NAcnt;
int InfCnt;

*pNumEls = 0;
*ppRetElements = NULL;
pSampleGenes = NULL;
AllocdSampleGenes = 0;
SampleGeneCnt= 0;
AboveThresCnt = 0;
NAcnt = 0;
InfCnt = 0;

pCSVDESeq = new CCSVFile();
if((Rslt = pCSVDESeq->Open(pszDESeqFile))!=eBSFSuccess)
	{
	while(pCSVDESeq->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pCSVDESeq->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open sample genes hit file '%s'",pszDESeqFile);
	delete pCSVDESeq;
	return(Rslt);
	}

while((Rslt=pCSVDESeq->NextLine()) > eBSFSuccess)
	{		
	NumFields = Rslt;
	LineNum = pCSVDESeq->GetLineNumber();
	if(NumFields < 9)			// unlikely to be DESeq file if less than 10 fields
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Too few fields at line %d, not a DESeq generated file: %s",LineNum,NumFields,pszDESeqFile);
		if(pSampleGenes != NULL)
			{
			delete pSampleGenes;
			pSampleGenes = NULL;
			}
		SampleGeneCnt = 0;
		Rslt = eBSFerrFileAccess;
		break;
		}
	if(LineNum == 1)			// 1st line expected to contain column titles, slough 
		continue;

	Rslt = pCSVDESeq->GetText(2,&pszGene);
	if(Rslt < eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to parse gene name at line %d from %s",LineNum,pszDESeqFile);
		if(pSampleGenes != NULL)
			{
			delete pSampleGenes;
			pSampleGenes = NULL;
			}
		SampleGeneCnt = 0;
		break;
		}
	if(!stricmp("NA",pszGene))	// if DESeq outputs gene names (and all other fields in this row) as "NA" if for any reason it can't determine the differential expression
		{
		NAcnt += 1;
		continue;				// simply slough this row
		}

	Rslt = pCSVDESeq->GetText(PMode == ePMDESeqPVal ? 8 : 9,&pszPVal);
	if(Rslt < eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to parse %s value at line %d from %s",PMode == ePMDESeqPVal ? "PVal" : "PValAdj",LineNum,pszDESeqFile);
		if(pSampleGenes != NULL)
			{
			delete pSampleGenes;
			pSampleGenes = NULL;
			}
		SampleGeneCnt = 0;
		break;
		}

	if(!stricmp("NA",pszPVal) || 
		!stricmp("-Inf",pszPVal) ||	// need to check for under/overflows from DESeq
		!stricmp("+Inf",pszPVal))
		{
		InfCnt += 1;
		continue;				// simply slough this row
		}

	PVal = strtod(pszPVal,&pszTerm);		// strtod is safer than atof in that under/overflows etc can be detected
	if(*pszTerm != '\0')
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to parse %s value at line %d from %s as a double",PMode == ePMDESeqPVal ? "PVal" : "PValAdj",LineNum,pszDESeqFile);
		if(pSampleGenes != NULL)
			{
			delete pSampleGenes;
			pSampleGenes = NULL;
			}
		SampleGeneCnt = 0;
		Rslt = eBSFerrParse;
		break;
		}

	if(PVal > Thres)				// only interested if PVal or PValAdj <= Thres
		{
		AboveThresCnt += 1;
		continue;
		}

	if(pSampleGenes == NULL || SampleGeneCnt >= AllocdSampleGenes)
		{
		AllocdSampleGenes += cElGenes2Alloc;
		if((pTmp = new tsElementGene[AllocdSampleGenes])==NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Memory allocation for sample genes failed");
			if(pSampleGenes != NULL)
				{
				delete pSampleGenes;
				pSampleGenes = NULL;
				}
			SampleGeneCnt = 0;
			Rslt = eBSFerrMem;
			break;
			}
		if(SampleGeneCnt)
			{
			memcpy(pTmp,(char *)pSampleGenes,sizeof(tsElementGene) * SampleGeneCnt);
			delete pSampleGenes;
			}
		pSampleGenes = pTmp;
		}
	strncpy(pSampleGenes[SampleGeneCnt].szGeneName,pszGene,cMaxFeatNameLen);
	pSampleGenes[SampleGeneCnt].szGeneName[cMaxFeatNameLen] = '\0';
	pSampleGenes[SampleGeneCnt].Length= 1;
	pSampleGenes[SampleGeneCnt++].Weighting = 1;
	}
if(Rslt >= eBSFSuccess)
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Accepted %d entries - another %d were above threshold, %d were NA and %d were 'Inf'", SampleGeneCnt,AboveThresCnt,NAcnt,InfCnt);

delete pCSVDESeq;
*ppRetElements = pSampleGenes;
*pNumEls = SampleGeneCnt;
return(Rslt);
}

// ParseCSVHits
// CSV hits file contains rows: gene name, weighting
// weighting is optional
int
ParseCSVHits(char *pszHitsFile,			// file to parse
			CBEDfile* pBED,				// hits must be to genes in this file
		  int MinLen,					// filter out elements less than this length
		  int MaxLen,					// filter out elements more than this length
		  int MinWeight,				// filter out elements less than this number of hits or weight
		  int MaxWeight,				// filter out elements more than this number of hits or weight
		  int *pNumEls,					// where to return number of elements in returned list				
		  tsElementGene **ppRetElements)	// returned list of elements
{
int Rslt;
int LineNum;
int GeneLength;
int SampleGeneCnt;
int AllocdSampleGenes;
tsElementGene *pSampleGenes;
tsElementGene *pTmp;
CCSVFile *pCSVSample;
char *pszGene;
int NumFields;
int GeneWeighting;
int CurFeatureID;

*pNumEls = 0;
*ppRetElements = NULL;
pSampleGenes = NULL;
AllocdSampleGenes = 0;
SampleGeneCnt= 0;

pCSVSample = new CCSVFile();
if((Rslt = pCSVSample->Open(pszHitsFile))!=eBSFSuccess)
	{
	while(pCSVSample->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pCSVSample->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open sample genes hit file '%s'",pszHitsFile);
	delete pCSVSample;
	return(Rslt);
	}
LineNum = 0;
while((Rslt=pCSVSample->NextLine()) > eBSFSuccess)
	{
	if (++LineNum == 1)			// 1st line expected to contain column titles, slough 
		continue;
	NumFields = Rslt;
	if(NumFields == 71)			// 71 fields if csv was generated by 'kit4b rnade'
		Rslt = pCSVSample->GetText(2, &pszGene);	// 2nd field is expected to be the gene name
	else
		Rslt = pCSVSample->GetText(1, &pszGene);	// first field is expected to be the gene name
	if (Rslt < eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process", "Unable to read gene name line %d from %s", pCSVSample->GetLineNumber(), pszHitsFile);
		if (pSampleGenes != NULL)
			{
			delete pSampleGenes;
			pSampleGenes = NULL;
			}
		SampleGeneCnt = 0;
		break;
		}
	// ensure gene is known ..
	CurFeatureID = pBED->LocateFeatureIDbyName(pszGene);
	if (CurFeatureID < 1)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to locate sample gene '%s' at line %d in %s", pszGene, pCSVSample->GetLineNumber(), pszHitsFile);
		continue;
		}
	GeneLength = pBED->GetFeatLen(CurFeatureID);
	if (GeneLength < MinLen || GeneLength > MaxLen)
		continue;


	if (NumFields >= 2) // 2nd field is optional, if present then expected to be the sample rating weighting or number of hits
		{
		if (NumFields == 71)
			pCSVSample->GetInt(4, &GeneWeighting);
		else
			Rslt = pCSVSample->GetInt(2, &GeneWeighting);
		if (Rslt < eBSFSuccess)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process", "Unable to read gene weighting at line %d from %s", pCSVSample->GetLineNumber(), pszHitsFile);
			if (pSampleGenes != NULL)
				{
				delete pSampleGenes;
				pSampleGenes = NULL;
				}
			SampleGeneCnt = 0;
			break;
			}
		if (GeneWeighting == 0)
			continue;
		}
	else
		GeneWeighting = 1;

	if (GeneWeighting < MinWeight || GeneWeighting > MaxWeight)
		continue;

	if(pSampleGenes == NULL || SampleGeneCnt >= AllocdSampleGenes)
		{
		AllocdSampleGenes += cElGenes2Alloc;
		if((pTmp = new tsElementGene[AllocdSampleGenes])==NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Memory allocation for %d sample genes failed",AllocdSampleGenes);
			if(pSampleGenes != NULL)
				{
				delete pSampleGenes;
				pSampleGenes = NULL;
				}
			SampleGeneCnt = 0;
			Rslt = eBSFerrMem;
			break;
			}

		if(SampleGeneCnt)
			{
			memcpy(pTmp,(char *)pSampleGenes,sizeof(tsElementGene) * SampleGeneCnt);
			delete pSampleGenes;
			}
		pSampleGenes = pTmp;
		}
	strncpy(pSampleGenes[SampleGeneCnt].szGeneName,pszGene,cMaxFeatNameLen);
	pSampleGenes[SampleGeneCnt].szGeneName[cMaxFeatNameLen] = '\0';
	pSampleGenes[SampleGeneCnt].Length= GeneLength;
	pSampleGenes[SampleGeneCnt++].Weighting = GeneWeighting;
	}

delete pCSVSample;
*ppRetElements = pSampleGenes;
*pNumEls = SampleGeneCnt;
return(eBSFSuccess);
}


// MTCNone
// No multiple test corrections applied
int
MTCNone(etOntologies Ontology,CGOTerms *pGOTerms)
{
tsGOTermCnts *pCnts;
int TotNumTermCnts = pGOTerms->GetNumTermCnts(Ontology);
for(int Cnt = 1; Cnt <= TotNumTermCnts; Cnt++)
	{
	pCnts = pGOTerms->GetTermCnts(Cnt);
	if(pCnts->RootOntology != Ontology ||
		!pCnts->NumSampleGenes)						// not all terms will have sample gene hits	
		continue;
	pCnts->HGMTCProbK = pCnts->HGProbK;
	pCnts->HGMTCWeightedProbK = pCnts->HGWeightedProbK;
	}
return(eBSFSuccess);
}

// MTCBonferroni
// Bonferroni - multiply Pr(K) for each term by number of genes in population which have associated terms
int
MTCBonferroni(etOntologies Ontology,CGOTerms *pGOTerms)
{
int Cnt;
double Fact;
tsGOTermCnts *pCnts;
int TotNumTermCnts;

TotNumTermCnts = pGOTerms->GetNumTermCnts(Ontology);
Fact = (double)pGOTerms->GetNumSampledTermCnts(Ontology,1);

for(Cnt = 1; Cnt <= TotNumTermCnts; Cnt++)
	{
	pCnts = pGOTerms->GetTermCnts(Cnt);
	if(pCnts->RootOntology != Ontology ||
		!pCnts->NumSampleGenes)						// not all terms will have sample gene hits	
		continue;
	pCnts->HGMTCProbK = pCnts->HGProbK * Fact;
	if(pCnts->HGMTCProbK > 1.0)			// ensure doesn't exceed 1.0
		pCnts->HGMTCProbK = 1.0;
	pCnts->HGMTCWeightedProbK = pCnts->HGWeightedProbK * Fact;
	if(pCnts->HGMTCWeightedProbK > 1.0)
		pCnts->HGMTCWeightedProbK = 1.0;
	}
return(eBSFSuccess);
}

// MTCHolm
// Holm - sort ascending the Pr(K) for each term (rank inversely n..1)
// Then multiply Pr(K) * min(Rank,
int
MTCHolm(etOntologies Ontology,CGOTerms *pGOTerms)
{
tsGOTermCnts *pCnts;
int TotNumTermCnts;
unsigned int CurRank;
unsigned int MaxRank;
double Fact;
int Cnt;
pGOTerms->SortCntsHGProbK();
TotNumTermCnts = pGOTerms->GetNumTermCnts(Ontology);
MaxRank = pGOTerms->GetNumSampledTermCnts(Ontology,1); // number of terms which have at least 1 sample
for(CurRank=MaxRank,Cnt = 1; Cnt <= TotNumTermCnts; Cnt++)
	{
	pCnts = pGOTerms->GetTermCnts(Cnt);
	if(pCnts->RootOntology != Ontology ||
		!pCnts->NumSampleGenes)						// not all terms will have sample gene hits	
		continue;
	Fact = (double)CurRank;
	pCnts->HGMTCProbK = pCnts->HGProbK * Fact;
	if(pCnts->HGMTCProbK > 1.0)
		pCnts->HGMTCProbK = 1.0;
	CurRank--;
	}

pGOTerms->SortCntsHGWeightedProbK();
for(CurRank=MaxRank,Cnt = 1; Cnt <= TotNumTermCnts; Cnt++)
	{
	pCnts = pGOTerms->GetTermCnts(Cnt);
	if(pCnts->RootOntology != Ontology ||
		!pCnts->NumSampleGenes)						// not all terms will have sample gene hits	
		continue;
	Fact = (double)CurRank;
	pCnts->HGMTCWeightedProbK = pCnts->HGWeightedProbK * Fact;
	if(pCnts->HGMTCWeightedProbK > 1.0)
		pCnts->HGMTCWeightedProbK = 1.0;
	CurRank--;
	}
return(eBSFSuccess);
}


// BenjaminiHochberg
// Benjamini Hochberg False Discovery Rate
// sort ascending the Pr(K) for each term (rank 1..n)
// Multiply each term Pr(K) by NumSampleGenes / rank
int
MTCBenjaminiHochberg(etOntologies Ontology,CGOTerms *pGOTerms)
{
int Cnt;
tsGOTermCnts *pCnts;
int TotNumTermCnts;
double Fact;
unsigned int MaxRank;
unsigned int CurRank;
pGOTerms->SortCntsHGProbK();
TotNumTermCnts = pGOTerms->GetNumTermCnts(Ontology);
MaxRank = pGOTerms->GetNumSampledTermCnts(Ontology,1);
for(CurRank=1,Cnt = 1; Cnt <= TotNumTermCnts; Cnt++)
	{
	pCnts = pGOTerms->GetTermCnts(Cnt);
	if(pCnts->RootOntology != Ontology ||
		!pCnts->NumSampleGenes)						// not all terms will have sample gene hits	
		continue;

	Fact = (double)MaxRank/CurRank;
	pCnts->HGMTCProbK = pCnts->HGProbK * Fact;
	if(pCnts->HGMTCProbK > 1.0)
		pCnts->HGMTCProbK = 1.0;
	CurRank++;
	}

pGOTerms->SortCntsHGWeightedProbK();
for(CurRank=1,Cnt = 1; Cnt <= TotNumTermCnts; Cnt++)
	{
	pCnts = pGOTerms->GetTermCnts(Cnt);
	if(pCnts->RootOntology != Ontology ||
		!pCnts->NumSampleGenes)						// not all terms will have sample gene hits	
		continue;
	Fact = (double)MaxRank/CurRank;
	pCnts->HGMTCWeightedProbK = pCnts->HGWeightedProbK * Fact;
	if(pCnts->HGMTCWeightedProbK > 1.0)
		pCnts->HGMTCWeightedProbK = 1.0;
	CurRank++;
	}
return(eBSFSuccess);
}



int
GenPValuesChiSquare(etOntologies Ontology,CGOTerms *pGOTerms,tsTermStats *pPopulationStats,tsTermStats *pSampleStats)
{
CStats *pStats = new CStats();
double QVal;
int Cells[4];
tsGOTermCnts *pCnts;
int NumTermCnts = pGOTerms->GetNumTermCnts(Ontology);

unsigned int TElements;     // number of genes in sample set
unsigned int KElements;		// number of genes in sample set which were annotated to term
unsigned int N1Elements;    // number of genes in population annotated to term
unsigned int N2Elements;    // number of genes in population not annotated to term (n1+n2 == population)
unsigned int TWElements;     // sum of all weights genes in sample set
unsigned int KWElements;	 // sum of all weights genes in sample set which were annotated to term
unsigned int N1WElements;    // sum of all weights genes in population annotated to term
unsigned int N2WElements;    // sum of all weights genes in population not annotated to term (n1+n2 == population)

for(int Cnt = 1; Cnt <= NumTermCnts; Cnt++)
	{
	pCnts = pGOTerms->GetTermCnts(Cnt);
	if(!pCnts->NumSampleGenes)
		continue;
	TElements = pSampleStats->NumGenesInGO;
	KElements = pCnts->NumSampleGenes;
	N1Elements = pCnts->NumBkgdGenes;
	N2Elements = pPopulationStats->NumGenesInGO - N1Elements;

	Cells[0] = TElements - KElements;
	Cells[1] = KElements;
	Cells[2] = N2Elements;
	Cells[3] = N1Elements;
	if((QVal = pStats->CalcChiSqr(2,2,Cells)) < 0.0)	// less than 0.0 if any expected count is less than 5
		pCnts->HGProbK = 2.0;							// down stream processing can recognise this value as being > 1.0 and thus an error
	else
		pCnts->HGProbK = pStats->ChiSqr2PVal(1,QVal);
	
	TWElements = pSampleStats->SumNumGenesInGOWeightings;
	KWElements = pCnts->SampleCnt;
	N1WElements = pCnts->BkgdCnt;
	N2WElements = pPopulationStats->SumNumGenesInGOWeightings - N1WElements;

	Cells[0] = TWElements - KWElements;
	Cells[1] = KWElements;
	Cells[2] = N2WElements;
	Cells[3] = N1WElements;
	if((QVal = pStats->CalcChiSqr(2,2,Cells)) < 0.0)	// less than 0.0 if any expected count is less than 5
		pCnts->HGWeightedProbK = 2.0;							// down stream processing can recognise this value as being > 1.0 and thus an error
	else
		pCnts->HGWeightedProbK = pStats->ChiSqr2PVal(1,QVal);
	}
return(eBSFSuccess);
}

int
GenPValuesFishersExactTest(etOntologies Ontology,CGOTerms *pGOTerms,tsTermStats *pPopulationStats,tsTermStats *pSampleStats)
{
CStats *pStats = new CStats();

tsGOTermCnts *pCnts;
int NumTermCnts = pGOTerms->GetNumTermCnts(Ontology);

unsigned int TElements;     // number of genes in sample set
unsigned int KElements;		// number of genes in sample set which were annotated to term
unsigned int N1Elements;    // number of genes in population annotated to term
unsigned int N2Elements;    // number of genes in population not annotated to term (n1+n2 == population)
unsigned int TWElements;     // sum of all weights genes in sample set
unsigned int KWElements;	 // sum of all weights genes in sample set which were annotated to term
unsigned int N1WElements;    // sum of all weights genes in population annotated to term
unsigned int N2WElements;    // sum of all weights genes in population not annotated to term (n1+n2 == population)

for(int Cnt = 1; Cnt <= NumTermCnts; Cnt++)
	{
	pCnts = pGOTerms->GetTermCnts(Cnt);
	if(!pCnts->NumSampleGenes)
		continue;
	TElements = pSampleStats->NumGenesInGO;
	KElements = pCnts->NumSampleGenes;
	N1Elements = pCnts->NumBkgdGenes;
	N2Elements = pPopulationStats->NumGenesInGO - N1Elements;
	pCnts->HGProbK = pStats->FishersExactTest(TElements - KElements, KElements, N2Elements, N1Elements);

	TWElements = pSampleStats->SumNumGenesInGOWeightings;
	KWElements = pCnts->SampleCnt;
	N1WElements = pCnts->BkgdCnt;
	N2WElements = pPopulationStats->SumNumGenesInGOWeightings - N1WElements;
	pCnts->HGWeightedProbK = pStats->FishersExactTest(TWElements - KWElements, KWElements, N2WElements, N1WElements);
	}
delete pStats;
return(eBSFSuccess);
}


int
SetBkgndCnts(etMTCMode MTC,		// multiple test correction mode
		bool bCanonicalise,     // canonicalise sample name isoforms by removing any numerical suffix in range '.[0-99] (TAIR ge
 		etOntologies Ontology, // which root ontology to process
		bool bProp,				// propagate counts from GO:Term into parent terms
		bool bBkgndLen,			// background counts proportional to gene lengths
		char Strand,			// background counts are for which strand '*', '+' or '-'
		int PopulationCnt,		// number of genes in background population 
		tsElementGene *pPopulationGenes, // genes to treat as being the background population
		CBEDfile *pBED,			// gene BED file
		CGOAssocs *pAssocs,		// gene to ID association file
		CGOTerms *pGOTerms,		// GO ontology file
		tsTermStats *pStats)	// returned processing stats
{
int Rslt;
int CurFeatureID;
char szGene[cMaxFeatNameLen];
char szChrom[cMaxDatasetSpeciesChrom];
int GeneStart;
int GeneEnd;
int GeneLen;
char GeneStrand;
int NumGOIDs;
int Cnt;
char *pszGOID;
int NumFeatures;
int TotGOIDsCounted;
char *pszGOIDs[1000];
char szGOIDs[10000];
char *pszConcatGOIDs;
int CntGOIDs;

int CurGeneIdx;
int CurGeneWeighting;
bool bMore;
bool bPopIsAll;
tsGOTerm *pTerm;
int NameInst;
int NumRedundantNames;
int NumGenes;
int NumGenesInGO;

memset((char *)pStats,0,sizeof(tsTermStats));

TotGOIDsCounted = 0;
CurFeatureID = 0;
CurGeneIdx = 0;
NumRedundantNames = 0;
NumGenes = 0;
NumGenesInGO = 0;
Rslt = eBSFSuccess;
bPopIsAll = (PopulationCnt == 0 || pPopulationGenes == NULL) ? true : false;
if(bPopIsAll)
	NumFeatures = pBED->GetNumFeatures();
	
bMore = true;
while(bMore)
	{
	if(bPopIsAll)		// use all genes as the background population?
		{
		if((CurFeatureID = pBED->GetNextFeatureID(CurFeatureID)) <= 0)
			break;
		NumGenes += 1;
		NameInst = pBED->GetNameInst(CurFeatureID);	// could be multiple instances so just accept the first as being the cannonical
		if(NameInst > 1)
			{
			NumRedundantNames += 1;
			continue;
			}
		pStats->TotNumGenes++;
		CurGeneWeighting = 1;
		}
	else				// use specified genes as the background population
		{
		if(CurGeneIdx >= PopulationCnt)
			break;
		CurGeneWeighting = pPopulationGenes[CurGeneIdx].Weighting;
		pStats->TotNumGenes++;
		NumGenes += 1;
		CurFeatureID = pBED->LocateFeatureIDbyName(pPopulationGenes[CurGeneIdx++].szGeneName);
		if(CurFeatureID <= 0) // not an error if can't locate gene
			{
			pStats->NumGenesNotInBED++;
			pStats->SumGenesNotInBEDWeightings += CurGeneWeighting;
			continue;
			}
		}

	Rslt = pBED->GetFeature(CurFeatureID,szGene,szChrom,&GeneStart,&GeneEnd,NULL,&GeneStrand);

	GeneLen = GeneEnd - GeneStart;
	if(bBkgndLen)
		CurGeneWeighting *= GeneLen;
	pStats->SumTotGeneWeightings += CurGeneWeighting;

	if(Strand != '*' && Strand != GeneStrand)
		{
		pStats->NumGenesNotOnStrand++;
		pStats->SumGenesNotOnStrandWeightings += CurGeneWeighting;
		continue;
		}
	
	if((NumGOIDs = pAssocs->GetNumGOIDs(szGene)) <= 0)
		{
		if(bCanonicalise && pAssocs->TrimNameIso(szGene))
			NumGOIDs = pAssocs->GetNumGOIDs(szGene);
		}
	CntGOIDs = 0;
	if(NumGOIDs > 0)
		{
		pszConcatGOIDs = szGOIDs;
		*pszConcatGOIDs = '\0';
		for(Cnt = 0; Cnt < NumGOIDs; Cnt++)
			{
			pszGOID = pAssocs->GetGOID(szGene,Cnt+1);
			if(pszGOID != NULL)
				{
				pTerm = pGOTerms->LocateGOID(pszGOID);
				if(pTerm == NULL || pTerm->RootOntology != Ontology)
					continue;
				strcpy(pszConcatGOIDs,pszGOID);
				pszGOIDs[CntGOIDs++] = pszConcatGOIDs;
				pszConcatGOIDs += strlen(pszGOID) + 1;
				}
			}
		}
	if(CntGOIDs)
		{
		if((Rslt=pGOTerms->AddCount(Ontology,bProp,false,CurGeneWeighting,CntGOIDs,pszGOIDs)) < eBSFSuccess)
			break;
		if(Rslt > 0)
			{
			pStats->NumGenesInGO++;
			pStats->SumNumGenesInGOWeightings += CurGeneWeighting;
			TotGOIDsCounted += Rslt;
			NumGenesInGO += 1;
			}
		else
			CntGOIDs = 0;
		}
	if(!CntGOIDs)
		{
		pStats->NumGenesNotInGO++;
		pStats->SumGenesNotInGOWeightings += CurGeneWeighting;
		}
	}
if(bPopIsAll)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processed %d background genes of which %d were redundant and %d had GO associations",NumGenes,NumRedundantNames,NumGenesInGO);

return(Rslt >= eBSFSuccess ? TotGOIDsCounted : Rslt);
}



int
SetSampleCnts(etMTCMode MTC,		// multiple test correction mode
		bool bCanonicalise,     // canonicalise sample name isoforms by removing any numerical suffix in range '.[0-99] (TAIR ge
  		etOntologies Ontology, // which root ontology to process
		bool bProp,				// propagate counts from GO:Term into parent terms
		char Strand,			// background counts are for which strand '*', '+' or '-'
		int SampleCnt,			// number of genes in sample 
		tsElementGene *pSampleGenes, // genes + counts to treat as being the sample
		CBEDfile *pBED,			// gene BED file
		CGOAssocs *pAssocs,		// gene to ID association file
		CGOTerms *pGOTerms,		// GO ontology file
		tsTermStats *pStats)	// returned processing stats
{
int Rslt;
int CurFeatureID;
char szGene[cMaxFeatNameLen];
char szChrom[cMaxDatasetSpeciesChrom];
int GeneStart;
int GeneEnd;
int GeneLen;
char GeneStrand;
int NumGOIDs;
int Cnt;
char *pszGOID;
int TotGOIDsCounted;

char *pszGOIDs[1000];
char szGOIDs[10000];
char *pszConcatGOIDs;
int CntGOIDs;

int CurGeneIdx;
int CurGeneWeighting;
bool bMore;
tsGOTerm *pTerm;

memset((char *)pStats,0,sizeof(tsTermStats));
CurFeatureID = 0;
CurGeneIdx = 0;
TotGOIDsCounted = 0;
bMore = true;
Rslt = eBSFSuccess;
while(bMore)
	{
	if(CurGeneIdx >= SampleCnt)
		break;
	CurGeneWeighting = pSampleGenes[CurGeneIdx].Weighting;
	pStats->TotNumGenes++;
	pStats->SumTotGeneWeightings += CurGeneWeighting;
	CurFeatureID = pBED->LocateFeatureIDbyName(pSampleGenes[CurGeneIdx++].szGeneName);
	if(CurFeatureID <= 0) // not a fatal error if can't locate gene
		{
		pStats->NumGenesNotInBED++;
		pStats->SumGenesNotInBEDWeightings += CurGeneWeighting;
		continue;
		}
	Rslt = pBED->GetFeature(CurFeatureID,szGene,szChrom,&GeneStart,&GeneEnd,NULL,&GeneStrand);
	if(Strand != '*' && (Strand != GeneStrand))
		{
		pStats->NumGenesNotOnStrand++;
		pStats->SumGenesNotOnStrandWeightings += CurGeneWeighting;
		continue;
		}
	GeneLen = GeneEnd - GeneStart;
	if((NumGOIDs = pAssocs->GetNumGOIDs(szGene)) <= 0)
		{
		if(bCanonicalise && pAssocs->TrimNameIso(szGene))
			NumGOIDs = pAssocs->GetNumGOIDs(szGene);
		}

	CntGOIDs = 0;
	if(NumGOIDs > 0)
		{
		pszConcatGOIDs = szGOIDs;
		*pszConcatGOIDs = '\0';
		CntGOIDs = 0;
		for(Cnt = 0; Cnt < NumGOIDs; Cnt++)
			{
			pszGOID = pAssocs->GetGOID(szGene,Cnt+1);
			if(pszGOID != NULL)
				{
				pTerm = pGOTerms->LocateGOID(pszGOID);
				if(pTerm == NULL || pTerm->RootOntology != Ontology)
					continue;

				strcpy(pszConcatGOIDs,pszGOID);
				pszGOIDs[CntGOIDs++] = pszConcatGOIDs;
				pszConcatGOIDs += strlen(pszGOID) + 1;
				}
			}
		}
	if(CntGOIDs)
		{
		if((Rslt=pGOTerms->AddCount(Ontology,bProp,true,CurGeneWeighting,CntGOIDs,pszGOIDs)) < eBSFSuccess)
			break;
		if(Rslt > 0)
			{
			pStats->NumGenesInGO++;
			pStats->SumNumGenesInGOWeightings += CurGeneWeighting;
			TotGOIDsCounted += Rslt;
			}
		else
			CntGOIDs = 0;
		}
	if(!CntGOIDs)
		{
		pStats->NumGenesNotInGO++;
		pStats->SumGenesNotInGOWeightings += CurGeneWeighting;
		}
	}
return(Rslt >= eBSFSuccess ? TotGOIDsCounted : Rslt);
}

// ReduceGeneList
// Sorts and then (if bKeepDups id false) reduces genes in the list by collapsing genes with same name into one instance
// When collapsing uses the maximum of all collapsed genes
// Assumes that any filtering for length and weight will have already been performed
int
ReduceGeneList(bool bKeepDups,int NumGenes,tsElementGene *pGeneList,int *pNumUniqueGenes,tsElementGene **ppUniqueGeneList)
{
int Idx;
tsElementGene *pGene;
tsElementGene *pPrevGene;
tsElementGene *pUniqueGene;
tsElementGene *pTmp;


int AllocdUniqueGenes;
int UniqueGenesCnt;
tsElementGene *pUniqueGeneList;

*ppUniqueGeneList = NULL;
*pNumUniqueGenes = 0;

if(!NumGenes || pGeneList == NULL)
	return(eBSFerrParams);

// firstly sort genes by gene name ascending
if(NumGenes > 1)
	qsort(pGeneList,NumGenes,sizeof(tsElementGene),SortGenes);

pGene = pGeneList;
pPrevGene = NULL;
pUniqueGene = NULL;
UniqueGenesCnt = 0;
AllocdUniqueGenes = 0;
pUniqueGeneList = NULL;
for(Idx = 0; Idx < NumGenes; Idx++,pGene++)
	{
	if(pPrevGene == NULL || bKeepDups || stricmp(pGene->szGeneName,pPrevGene->szGeneName))
		{
		pPrevGene = pGene;
		if(pUniqueGeneList == NULL || UniqueGenesCnt >= AllocdUniqueGenes)
			{
			AllocdUniqueGenes += cAllocUniqueGenes;
			if((pTmp = new tsElementGene[AllocdUniqueGenes])==NULL)
				{
				if(pUniqueGeneList != NULL)
					delete pUniqueGeneList;
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory for unique genes list");
				return(eBSFerrMem);
				}
			if(UniqueGenesCnt)
				{
				memcpy(pTmp,(char *)pUniqueGeneList,sizeof(tsElementGene) * UniqueGenesCnt);
				delete pUniqueGeneList;
				}
			pUniqueGeneList = pTmp;
			}
		pUniqueGene = &pUniqueGeneList[UniqueGenesCnt++];
		strcpy(pUniqueGene->szGeneName,pGene->szGeneName);
		pUniqueGene->Weighting = pGene->Weighting;
		}
	else
		{	
		if(pUniqueGene->Weighting < pGene->Weighting)
			pUniqueGene->Weighting = pGene->Weighting;
		}
	}
*ppUniqueGeneList = pUniqueGeneList;
*pNumUniqueGenes = UniqueGenesCnt;
return(eBSFSuccess);
}


// grey1..100
char *
MapProbNodeAttribs(double Prb,double dblDotThres)
{
static char szAttrib[200];
int Level;
if(Prb > dblDotThres)
	return((char *)"");

Level = (int)floor(50.0 + (50*Prb/dblDotThres));
if(Level > 100)
	Level = 100;
sprintf(szAttrib,"color=grey%1d style=filled",Level);
return(szAttrib);
}

char *
GetRootAttrib(tsGOTerm *pParent)
{
switch(pParent->RootOntology) {
	case eONTCellular:	// Cellular component
		return((char *)"CELL");
	case eONTBiological:// Biological process
		return((char *)"BIO");
	case eONTMolecular:	// Molecular function
		return((char *)"MOL");
	default:
		break;
	}
return((char *)" ");
}

double
GetDotPValue(tsGOTermCnts *pCnts,etDotPValue DotPValue)
{
switch(DotPValue) {
	case eDPVunweighted:				// unweighted and no multiple term corection
		return(pCnts->HGProbK);
	case eDPVweighted:					// weighted but with no multiple term correction
		return(pCnts->HGWeightedProbK);
	case eDPVunweightedMTC:				// unweighted and multiple term corection
		return(pCnts->HGMTCProbK);
	case eDPVweightedNTC:				// weighted and multiple term correction
		return(pCnts->HGMTCWeightedProbK);
	default:
		break;
	}
return(1.0);
}

int
GenGOdot(int ParentTermID,FILE *pDot,CGOTerms *pGOTerms, int *pTermStates,double dblDotThres,etDotPValue DotPValue)
{
int Rslt;
tsGOTerm *pParent;
tsGOTerm *pChild;
tsGOTermCnts *pCnts;
double ProbK;
int NthChild;
int ChildTermID;
int *pParentTermState;
int *pChildTermState;
char *pszParentGOTermID;
char *pszParentGOName;
char *pszChildGOTermID;
char *pszChildGOName;

pParent = pGOTerms->GetGOTermByID(ParentTermID);
pszParentGOTermID = pGOTerms->GetGOIDforTerm(pParent);
pszParentGOName = pGOTerms->GetGONameforTerm(pParent);
pParentTermState = &pTermStates[ParentTermID];

pCnts = pGOTerms->GetExistingTermCnts(pParent);
if(pCnts != NULL && pCnts->NumSampleGenes)		
	{
	if(*pParentTermState == 0)	// if first time for this term
		{
		*pParentTermState = 0x01;
		ProbK = GetDotPValue(pCnts,DotPValue);
		if(ProbK <= dblDotThres)
			{
			fprintf(pDot,"\t N%d [%s label = \"%s\\n%s\\n%s\\nP(K)=%g\"];\n",ParentTermID,MapProbNodeAttribs(ProbK,dblDotThres),GetRootAttrib(pParent),pszParentGOTermID,pszParentGOName,ProbK);
			*pParentTermState = 0x03;
			}
		}
	}
else
	*pParentTermState = 0x01;

NthChild = 1;
while((pChild = pGOTerms->LocateNthChild(pParent,NthChild++))!=NULL)
	{
	ChildTermID = pChild->GOTermID;
	pChildTermState = &pTermStates[ChildTermID];
	pCnts = pGOTerms->GetExistingTermCnts(pChild);
	if(pCnts != NULL && pCnts->NumSampleGenes)
		{
		ProbK = GetDotPValue(pCnts,DotPValue);
		if(ProbK <= dblDotThres)
			{
			if(*pParentTermState & 0x02)
				{
				pszChildGOTermID = pGOTerms->GetGOIDforTerm(pChild);
				pszChildGOName = pGOTerms->GetGONameforTerm(pChild);
				fprintf(pDot,"\t N%d -> N%d;\n",ParentTermID,ChildTermID);
				}
			}
		}

	if(!(*pChildTermState & 0x04))
		{		
		if(Rslt = GenGOdot(ChildTermID,pDot,pGOTerms,pTermStates,dblDotThres,DotPValue) < 0)
			return(Rslt);
		}
	*pChildTermState |= 0x04;
	}
return(0);
}

//
// OutputResults
// Writes results out to file
int
OutputResults(etOntologies Ontology,
  			  tsTermStats *pPopulationStats,
			  tsTermStats *pSampleStats,
			  char *pszResultsFile,
			  char *pszGraphViz,
			  CGOTerms *pGOTerms,
			  char *pszRootGOTerm,
			  double dblDotThres,
			  etDotPValue DotPValue)	// which P-Value to use when generating GraphViz dot file
{
FILE *pRslts;
FILE *pDot;

int NumTermCnts;
tsGOTerm *pTerm;
char *pszGOTermID;
char *pszGOName;
tsGOTermCnts *pCnts;
int Cnt;
double Expected;
double Enrichment;
double WExpected;
double WEnrichment;

if((pRslts = fopen(pszResultsFile,"w+"))==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open/create results file %s error: %s",pszResultsFile,strerror(errno));
	return(eBSFerrOpnFile);
	}
fprintf(pRslts,"\"GO:Term\",\"Name\",\"P(K)\",\"Pw(K)\",\"Pmtc(K)\",\"Pmtcw(K)\",\"PopGenes\",\"WPopGenes\",\"SampleGenes\",\"WSampleGenes\",\"ExpectedGenes\",\"WExpectedGenes\",\"Enrichment\",\"WEnrichment\",\"TotPopulation\",\"WTotPopulation\",\"TotSample\",\"WTotSample\",\"Ontology Class\"\n");
pGOTerms->SortCntsHGMTCProbK();
NumTermCnts = pGOTerms->GetNumTermCnts(Ontology);
for(Cnt = 1; Cnt <= NumTermCnts; Cnt++)
	{	
	pCnts = pGOTerms->GetTermCnts(Cnt);
	if(!pCnts->NumSampleGenes)		
		continue;
	pTerm = pGOTerms->GetGOTermByID(pCnts->TermID);
	pszGOTermID = pGOTerms->GetGOIDforTerm(pTerm);
	pszGOName = pGOTerms->GetGONameforTerm(pTerm);

	Expected = ((double)(pSampleStats->NumGenesInGO * pCnts->NumBkgdGenes)) / pPopulationStats->NumGenesInGO;
	Enrichment = pCnts->NumSampleGenes / (Expected > 0.0 ? Expected : 0.001);

	WExpected = ((double)(pSampleStats->SumNumGenesInGOWeightings * pCnts->BkgdCnt)) / pPopulationStats->SumNumGenesInGOWeightings;
	WEnrichment = pCnts->SampleCnt / (WExpected > 0.0 ? WExpected : 0.001);


	fprintf(pRslts,"\"%s\",\"%s\",%g,%g,%g,%g,%d,%d,%d,%d,%1.4f,%1.4f,%1.4f,%1.4f,%d,%d,%d,%d,\"%s\"\n",
		pszGOTermID,pszGOName,
		pCnts->HGProbK,pCnts->HGWeightedProbK,pCnts->HGMTCProbK,pCnts->HGMTCWeightedProbK,
		pCnts->NumBkgdGenes,pCnts->BkgdCnt,pCnts->NumSampleGenes,pCnts->SampleCnt,
		Expected,WExpected,Enrichment,WEnrichment,
		pPopulationStats->NumGenesInGO,pPopulationStats->SumNumGenesInGOWeightings,pSampleStats->NumGenesInGO,pSampleStats->SumNumGenesInGOWeightings,
		pGOTerms->RootOntology2Txt((etOntologies)pTerm->RootOntology));
	}
fclose(pRslts);


int RootID;
int *pTermState;

// now for the GraphViz .dot file
if(pszGraphViz == NULL || pszGraphViz[0] == '\0')
	return(eBSFSuccess);

if((pDot = fopen(pszGraphViz,"w+"))==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open/create results file %s error: %s",pszGraphViz,strerror(errno));
	return(eBSFerrOpnFile);
	}
fprintf(pDot,"digraph GOgraph {");
if((pTermState = new int [cMaxGOterms]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory for GraphViz term states");
	fclose(pDot);
	return(eBSFerrMem);
	}
memset(pTermState,0,sizeof(int) * cMaxGOterms);

if(pszRootGOTerm != NULL && pszRootGOTerm[0] != '\0')
	{
	if((pTerm = pGOTerms->LocateGOID(pszRootGOTerm))!= NULL)
		{
		if(pTerm->RootOntology == Ontology)
			GenGOdot(pTerm->GOTermID,pDot,pGOTerms,pTermState,dblDotThres,DotPValue);
		else
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Requested root term %s not in selected ontology class",pszRootGOTerm);
			fclose(pDot);
			delete pTermState;
			return(eBSFerrGOID);
			}
		}
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Requested root term %s not located",pszRootGOTerm);
	fclose(pDot);
	delete pTermState;
	return(eBSFerrGOID);
	}
else
	{
	RootID = pGOTerms->GetRootGOTermID(Ontology);
	GenGOdot(RootID,pDot,pGOTerms,pTermState,dblDotThres,DotPValue);
	}

fprintf(pDot,"}\n");
fclose(pDot);
delete pTermState;
return(eBSFSuccess);
}


static 
int SortGenes(const void *arg1, const void *arg2)
{
tsElementGene *pGene1 = (tsElementGene *)arg1;
tsElementGene *pGene2 = (tsElementGene *)arg2;
char cGene1 = tolower(pGene1->szGeneName[0]);
char cGene2 = tolower(pGene2->szGeneName[0]);
if(cGene1 < cGene2)
	return(-1);
if(cGene1 > cGene2)
	return(1);
return(stricmp(pGene1->szGeneName,pGene2->szGeneName));
}
