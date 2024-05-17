/*
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
#include "./CGenMAFAlgn.h"



// CreateMAlignFile
int		// Create a multialignment file from files (either axt or mfa format) in specified source directory
	ProcCreateMAlignFile(char* pszChromLens, char* pszSrcDirPath, char* pszDestAlignFile,
	char* pszDescr, char* pszTitle,
	char* pszRefSpecies, char* pszRelSpecies, bool bIsAXT, bool bSwapRefRel);

int
	ProcGenbioDataPointsFile(char* pszMAF, char* pszDataPointsFile, char* pszDescr, char* pszTitle);

#ifdef _WIN32
int genmafalgn(int argc, char* argv[])
{
	// determine my process name
	_splitpath(argv[0], nullptr, nullptr, gszProcName, nullptr);
#else
int
genmafalgn(int argc, char** argv)
{
	// determine my process name
	CUtility::splitpath((char*)argv[0], nullptr, gszProcName);
#endif
	int iScreenLogLevel;		// level of screen diagnostics
	int iFileLogLevel;			// level of file diagnostics
	char szLogFile[_MAX_PATH];	// write diagnostics to this file

	bool bIsAXT;
	bool bSwapRefRel;
	int Rslt;
	int iMode;

	char szOutputFileSpec[_MAX_PATH];
	char szInputFileSpec[_MAX_PATH];
	char szChromLens[_MAX_PATH];
	char szRefSpecies[80];
	char szRelSpecies[80];
	char szDescription[cMBSFFileDescrLen];
	char szTitle[cMBSFShortFileDescrLen];


	// command line args
	struct arg_lit* help = arg_lit0("hH", "help", "print this help and exit");
	struct arg_lit* version = arg_lit0("v", "version,ver", "print version information and exit");
	struct arg_int* FileLogLevel = arg_int0("f", "FileLogLevel", "<int>", "Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
	struct arg_int* ScreenLogLevel = arg_int0("S", "ScreenLogLevel", "<int>", "Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
	struct arg_file* LogFile = arg_file0("F", "FileLogLevel", "<file>", "diagnostics log file");

	struct arg_lit* IsAXT = arg_lit0("x", "axt", "MAF (default) or AXT two species source files");
	struct arg_file* ChromLens = arg_file0("c", "chromlens", "<file>", "input file containing species.chrom lengths, required for .axt processing");
	struct arg_int* Mode = arg_int0("m", "mode", "<int>", "processing mode: 0 (default) MAF or AXT -> multialign, 1 multialign -> datapoints");
	struct arg_file* InFile = arg_file1("i", nullptr, "<file>", "input from .maf or axt files (wildcards accepted)");
	struct arg_file* OutFile = arg_file1("o", nullptr, "<file>", "output as multialign biosequence file");
	struct arg_str* Descr = arg_str1("d", "descr", "<string>", "full description");
	struct arg_str* Title = arg_str1("t", "title", "<string>", "short title");
	struct arg_str* RefSpecies = arg_str1("r", "ref", "<string>", "reference species ");
	struct arg_str* RelSpecies = arg_str0("R", "rel", "<string>", "relative species (axt only) ");
	struct arg_lit* SwapRefRel = arg_lit0("X", "exchange", "exchange ref and rel species - only applies to AXT alignments");


	struct arg_end* end = arg_end(20);

	void* argtable[] = { help,version,FileLogLevel,ScreenLogLevel,LogFile,SwapRefRel,IsAXT,ChromLens,Mode,InFile,OutFile,Descr,Title,RefSpecies,RelSpecies,end };

	char** pAllArgs;
	int argerrors;
	argerrors = CUtility::arg_parsefromfile(argc, (char**)argv, &pAllArgs);
	if (argerrors >= 0)
		argerrors = arg_parse(argerrors, pAllArgs, argtable);

	/* special case: '--help' takes precedence over error reporting */
	if (help->count > 0)
	{
		printf("\n%s ", gszProcName);
		arg_print_syntax(stdout, argtable, "\n");
		arg_print_glossary(stdout, argtable, "  %-25s %s\n");
		printf("\nNote: Parameters can be entered into a parameter file, one parameter per line.");
		printf("\n      To invoke this parameter file then precede its name with '@'");
		printf("\n      e.g. %s @myparams.txt\n\n", gszProcName);
		exit(1);
	}

	/* special case: '--version' takes precedence error reporting */
	if (version->count > 0)
	{
		printf("\n%s Version %s", gszProcName, kit4bversion);
		exit(1);
	}


	if (!argerrors)
	{
		iScreenLogLevel = ScreenLogLevel->count ? ScreenLogLevel->ival[0] : eDLInfo;
		if (iScreenLogLevel < eDLNone || iScreenLogLevel > eDLDebug)
		{
			printf("\nError: ScreenLogLevel '-S%d' specified outside of range %d..%d", iScreenLogLevel, eDLNone, eDLDebug);
			exit(1);
		}
		if (FileLogLevel->count && !LogFile->count)
		{
			printf("\nError: FileLogLevel '-f%d' specified but no logfile '-F<logfile>'", FileLogLevel->ival[0]);
			exit(1);
		}

		iFileLogLevel = FileLogLevel->count ? FileLogLevel->ival[0] : eDLInfo;
		if (iFileLogLevel < eDLNone || iFileLogLevel > eDLDebug)
		{
			printf("\nError: FileLogLevel '-l%d' specified outside of range %d..%d", iFileLogLevel, eDLNone, eDLDebug);
			exit(1);
		}
		if (LogFile->count)
		{
			strncpy(szLogFile, LogFile->filename[0], _MAX_PATH);
			szLogFile[_MAX_PATH - 1] = '\0';
		}
		else
		{
			iFileLogLevel = eDLNone;
			szLogFile[0] = '\0';
		}

		bIsAXT = IsAXT->count ? true : false;
		bSwapRefRel = SwapRefRel->count ? true : false;

		strcpy(szInputFileSpec, InFile->filename[0]);
		strcpy(szOutputFileSpec, OutFile->filename[0]);
		if (bIsAXT)
		{
			if (ChromLens->count)
			{
				strncpy(szChromLens, ChromLens->filename[0], _MAX_PATH);
				szChromLens[_MAX_PATH - 1] = '\0';
			}
			else
			{
				printf("\nError: Processing .axt files but no chromosome length file '-c<chromlen_file>' specified");
				exit(1);
			}
		}
		else
			szChromLens[0] = '\0';

		iMode = Mode->count ? Mode->ival[0] : 0;

		if (Title->count < 1)
		{
			printf("No short title specified with '-t' parameter\n");
			exit(1);
		}
		if (Descr->count < 1)
		{
			printf("No descriptive text specified with '-d' parameter\n");
			exit(1);
		}

		strncpy(szTitle, Title->sval[0], sizeof(szTitle));
		szTitle[cMBSFShortFileDescrLen - 1] = '\0';
		strncpy(szDescription, Descr->sval[0], sizeof(szDescription));
		szDescription[cMBSFFileDescrLen - 1] = '\0';


		if (RefSpecies->count)
			strcpy(szRefSpecies, RefSpecies->sval[0]);
		else
			szRefSpecies[0] = '\0';

		if (RelSpecies->count)
			strcpy(szRelSpecies, RelSpecies->sval[0]);
		else
			szRelSpecies[0] = '\0';

		if (iMode == 0)
		{
			if (szRefSpecies[0] == '\0')
			{
				printf("\nError: In Mode 0 (creating biosequence multiple alignment file) reference species must be specified as '-r RefSpecies'");
				exit(1);
			}

			if (bIsAXT && szRelSpecies[0] == '\0')
			{
				printf("\nError: In Mode 0 (creating biosequence multiple alignment file) and processing .AXT files, relative species must be specified as '-R RelSpecies'");
				exit(1);
			}
		}

		// now that command parameters have been parsed then initialise diagnostics log system
		if (!gDiagnostics.Open(szLogFile, (etDiagLevel)iScreenLogLevel, (etDiagLevel)iFileLogLevel, true))
		{
			printf("\nError: Unable to start diagnostics subsystem.");
			if (szLogFile[0] != '\0')
				printf(" Most likely cause is that logfile '%s' can't be opened/created", szLogFile);
			exit(1);
		}

		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Version: %s Processing parameters:", kit4bversion);


		gStopWatch.Start();
#ifdef _WIN32
		SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
		switch (iMode) {
		case 0:	// creating bioseq multialignment file
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "Version: %s Processing parameters:", kit4bversion);
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "Mode: 0 (creating bioseq .ALGN multiple alignment file from .AXT or MAF files)");
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "Output file: '%s'", szOutputFileSpec);
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "Source files: '%s'", szInputFileSpec);
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "Source files are expected to be %s alignments", bIsAXT ? "AXT" : "MAF");
			if (bIsAXT)
				gDiagnostics.DiagOutMsgOnly(eDLInfo, "Chromosome lengths file: '%s'", szChromLens);
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "Reference Species: %s", szRefSpecies);
			if (bIsAXT)
				gDiagnostics.DiagOutMsgOnly(eDLInfo, "Relative Species: %s", szRelSpecies);
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "Exchange Ref/Rel alignments: %s", bSwapRefRel ? "yes" : "no");
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "Title text: %s", szTitle);
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "Descriptive text: %s", szDescription);
			Rslt = ProcCreateMAlignFile(szChromLens, szInputFileSpec, szOutputFileSpec, szDescription, szTitle, szRefSpecies, szRelSpecies, bIsAXT, bSwapRefRel);
			break;

		case 1:	// creating bioseq data points file from bioseq multialignment file
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "Version: %s Processing parameters:", kit4bversion);
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "Mode: 1 (creating bioseq data points .DPS file from bioseq .ALGN file)");
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "Output file: '%s'", szOutputFileSpec);
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "Source file: '%s'", szInputFileSpec);
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "Title text: %s", szTitle);
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "Descriptive text: %s", szDescription);
			Rslt = ProcGenbioDataPointsFile(szInputFileSpec, szOutputFileSpec, szDescription, szTitle);
			break;

		default:
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "requested processing mode %d not supported\n", iMode);
			exit(1);
		}

		gStopWatch.Stop();
		Rslt = Rslt >= 0 ? 0 : 1;
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Exit code: %d Total processing time: %s", Rslt, gStopWatch.Read());
		exit(Rslt);
	}
	else
	{
		arg_print_errors(stdout, end, gszProcName);
		arg_print_syntax(stdout, argtable, "\nend of help\n");
		exit(1);
	}
	exit(Rslt);
}

int		// Create a multialignment file from files (either axt or mfa format) in specified source directory
ProcCreateMAlignFile(char* pszChromLens, char* pszSrcDirPath, char* pszDestAlignFile,
	char* pszDescr, char* pszTitle,
	char* pszRefSpecies, char* pszRelSpecies, bool bIsAXT, bool bSwapRefRel)
{
int Rslt;
CGenMAFAlgn* pCGenMAFAlgn;
if ((pCGenMAFAlgn = new CGenMAFAlgn) == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to instantiate instance of CRNA_DE");
	return(eBSFerrObj);
	}
Rslt = pCGenMAFAlgn->CreateMAlignFile(pszChromLens, pszSrcDirPath, pszDestAlignFile,pszDescr,pszTitle,pszRefSpecies,pszRelSpecies,bIsAXT,bSwapRefRel);
if (pCGenMAFAlgn != nullptr)
	delete pCGenMAFAlgn;
return(Rslt);
}


int
ProcGenbioDataPointsFile(char* pszMAF, char* pszDataPointsFile, char* pszDescr, char* pszTitle)
{
int Rslt;
CGenMAFAlgn* pCGenMAFAlgn;
if ((pCGenMAFAlgn = new CGenMAFAlgn) == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to instantiate instance of CRNA_DE");
	return(eBSFerrObj);
	}
Rslt = pCGenMAFAlgn->GenbioDataPointsFile(pszMAF, pszDataPointsFile, pszDescr, pszTitle);
if (pCGenMAFAlgn != nullptr)
	delete pCGenMAFAlgn;
return(Rslt);
}

CGenMAFAlgn::CGenMAFAlgn()
{
m_pAlignValidate = nullptr;
}

CGenMAFAlgn::~CGenMAFAlgn()
{
if (m_pAlignValidate != nullptr)
	delete m_pAlignValidate;
}

void 
CGenMAFAlgn::Reset(void)
{
if (m_pAlignValidate != nullptr)
	{
	delete m_pAlignValidate;
	m_pAlignValidate = nullptr;
	}

}

bool
CGenMAFAlgn::IsBase(char Base)
{
	if (Base == 'a' || Base == 'A' ||
		Base == 'c' || Base == 'C' ||
		Base == 'g' || Base == 'G' ||
		Base == 't' || Base == 'T')
		return(true);
	return(false);
}

char*
CGenMAFAlgn::StripWS(char* pChr)
{
	while (*pChr != '\0' && isspace(*pChr))
		pChr++;
	return(pChr);
}

int
CGenMAFAlgn::GenbioDataPointsFile(char* pszMAF, char* pszDataPointsFile, char* pszDescr, char* pszTitle)
{
int Rslt;
CMAlignFile* pAlignFile;
Reset();
if((m_pAlignValidate = new CAlignValidate)==nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to instantiate instance of CAlignValidate");
	return(eBSFerrObj);
	}

if ((pAlignFile = new CMAlignFile) == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to instantiate instance of CMAlignFile");
	return(eBSFerrObj);
	}
Rslt = pAlignFile->MultiAlignAll2DataPts(pszMAF, pszDataPointsFile, pszDescr, pszTitle);
if (Rslt < 0)
{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenbioDataPointsFile: Errors creating '%s' from '%s'\n", pszDataPointsFile, pszMAF);
	while (pAlignFile->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error-> %s", pAlignFile->GetErrMsg());
}

delete pAlignFile;
return(Rslt);
}

// ProcessThisFile
// Parse input multialignment file. Format is expected to be either MAF or AXT
// as specified by (tsProcParams *)pParams.
int
CGenMAFAlgn::ProcessThisFile(char* pszFile,
	void* pParams)			// will be tsProcParams
{
	int Rslt;
	int hFile;
	unsigned long Elapsed;
	unsigned long CurElapsed;
	CStopWatch CurTime;
	gzFile pgzFile;

	char* pNxtChr;
	char Chr;
	int iBuffLen;							// number of buffered chrs in szBuffer
	int LineLen;							// current number of chars in pszLineBuffer;
	int64_t CurLineNumb;					// current line being processed 
	tsProcParams* pProcParams = (tsProcParams*)pParams;
	bool bIsAXT = pProcParams->bIsAXT;
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Processing File: %s", pszFile);
	pgzFile = nullptr;
	hFile = -1;

	int NameLen = (int)strlen(pszFile);
	if (NameLen >= 4 && !stricmp(".gz", &pszFile[NameLen - 3]))
	{
		if ((pgzFile = gzopen(pszFile, "r")) == nullptr)
		{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to open %s as a gzip'd file - %s", pszFile, strerror(errno));
			return(eBSFerrOpnFile);
		}
		if (gzbuffer(pgzFile, cMARawBuffSize * 2) != 0)
		{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to set gzbuffer size to %d", cgzAllocInBuffer);
			gzclose(pgzFile);
			return(eBSFerrMem);
		}
	}
	else
	{
#ifdef _WIN32
		hFile = open(pszFile, _O_RDONLY | _O_SEQUENTIAL);
#else
		hFile = open(pszFile, O_RDWR);
#endif
		if (hFile == -1)
		{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessThisFile: Unable to open input file for processing - '%s' - %s",
				pszFile, strerror(errno));
			return(eBSFerrOpnFile);
		}
	}

	Rslt = pProcParams->pAlignFile->StartFile(pszFile);
	if (Rslt != eBSFSuccess)
		return(Rslt);
	
	CurTime.Start();
	Elapsed = CurTime.ReadUSecs();
	CurLineNumb = 0;
	LineLen = 0;
	while (Rslt >= eBSFSuccess && ((iBuffLen = pgzFile != nullptr ? gzread(pgzFile, pProcParams->pRawBuffer, cMARawBuffSize) : read(hFile, pProcParams->pRawBuffer, cMARawBuffSize)) > 0))
	{
		pNxtChr = pProcParams->pRawBuffer;
		int BuffIdx = 0;
		while (BuffIdx++ < iBuffLen)
		{
			Chr = *pNxtChr++;
			if (Chr == '\r')		// slough CR, use NL as EOL 
				continue;
			if (Chr == '\n')		// EOL
			{
			CurLineNumb++;
			if (!LineLen)
				continue;
			if((CurLineNumb % 100000) == 0)
				{
				CurElapsed = CurTime.ReadUSecs();
				if ((CurElapsed - Elapsed) > 60)
					{
					Elapsed = CurElapsed;
					gDiagnostics.DiagOut(eDLInfo, gszProcName, "Parsing line %zd in file '%s'", CurLineNumb, pszFile);
					}
				}
			pProcParams->pszLineBuffer[LineLen] = 0;
			if (bIsAXT)
				Rslt = ProcessAXTline(LineLen, pProcParams);
			else
				Rslt = ProcessMAFline(LineLen, pProcParams);
			if (Rslt < 0)
			{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessThisFile: Errors processing %s line", bIsAXT ? "AXT" : "MAF");
				return(Rslt);
			}
			LineLen = 0;
			continue;
			}
			if (!LineLen)		// strip any leading space
			{
				if (isspace(Chr))
					continue;
			}
			pProcParams->pszLineBuffer[LineLen++] = Chr;
		}
	}
	if (LineLen)
	{
		pProcParams->pszLineBuffer[LineLen] = 0;
		if (bIsAXT)
			Rslt = ProcessAXTline(LineLen, pProcParams);
		else
			Rslt = ProcessMAFline(LineLen, pProcParams);
		if (Rslt < 0)
		{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessThisFile: Errors processing %s line", bIsAXT ? "AXT" : "MAF");
			return(Rslt);
		}
	}
	pProcParams->pAlignFile->EndAlignBlock();
	pProcParams->pAlignFile->EndFile();

	pgzFile != nullptr ? gzclose(pgzFile) : close(hFile);
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessThisFile: Total alignment blocks accepted: %d", pProcParams->pAlignFile->GetNumBlocks());
	return(Rslt);
}

// ProcessMAFline
// Process line from assumed MAF formated file
//
//Each sequence line begins with "s" and contains the following required fields: 
//src -- Name of one of the source sequences included in the alignment. For sequences that are resident in a browser assembly, the form database.chromosome allows automatic creation of links to other assemblies. Non-browser sequences are typically referenced by the species name alone. Species names must not contain spaces: concatenate multi-word names or replace the space with an underscore. 
//start -- Start of the aligning region in the source sequence, using zero-based position coordinates. If the strand value is "-", this field defines the start relative to the reverse-complemented source sequence. 
//size -- Size of the aligning region in the source sequence. This is equal to the number of non-dash characters in the text field (see below). 
//strand -- "+" or "-". If the value is "-", the sequence aligns to the reverse-complemented source. 
//srcSize -- Size of the entire source sequence, not just the portions involved in the alignment. 
//text -- Nucleotides (or amino acids) in the alignment are represented by upper-case letters; repeats are shown in lower case. Insertions are indicated by "-". 
//
// Notes:
// MAF start coordinates are zero based (0..n) but when displayed in the UCSC browser they are one based (1..n+1)
int
CGenMAFAlgn::ProcessMAFline(int LineLen, tsProcParams * pProcParams)
{
	static int32_t Score;						// alignment line score or pass value
	static int32_t RefChromID;					// reference chromosome
	static bool bRefChromNxt;				// true if next chromosome is the first in block - assume this is the reference
	char szSpecies[cMaxDatasetSpeciesChrom];
	char szChrom[cMaxDatasetSpeciesChrom];
	uint32_t Start;
	uint32_t Len;
	uint32_t AlignLen;
	int32_t ChromLen;
	char cStrand;
	int32_t Psn;
	int32_t Cnt;
	int Rslt = eBSFSuccess;
	char* pszLine = pProcParams->pszLineBuffer;
	switch (*pszLine) {
	case '#':					// parameter line or header
		if (pszLine[1] == '#')	// '##' is a header line, should have "##maf varname=value ..." format
		{
			pszLine = StripWS(pszLine + 6);
			if (*pszLine != '\0')
				Rslt = pProcParams->pAlignFile->AddHeader(pszLine);
		}
		else						// '#' is a parameter line containing parameters used to run the alignment program...
		{
			pszLine = StripWS(pszLine + 1);
			if (*pszLine != '\0')
				Rslt = pProcParams->pAlignFile->AddParameters(pszLine);
		}
		return(Rslt);

	case 'a':		// new alignment block starting
		pszLine = StripWS(pszLine + 1);
		if (*pszLine != '\0')
		{
			if (!sscanf(pszLine, "score=%d", &Score))	// note: only interested in integral part of score (scores are floating points)
				Score = 0;
		}
		else
			Score = 0;

		Rslt = pProcParams->pAlignFile->StartAlignBlock(Score);		// start new block, any existing will be autoclosed
		if (Rslt < eBSFSuccess)
		{
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "ProcessMAFline: Unable to start alignment block - %d - %s", Rslt, CErrorCodes::ErrText((teBSFrsltCodes)Rslt));
			break;
		}
		bRefChromNxt = true;
		break;

	case 's':		// species sequence
		pszLine = StripWS(pszLine + 1);

		Cnt = sscanf(pszLine, "%s %u %u %c %u %n",
			szSpecies, &Start, &Len, &cStrand, &ChromLen, &Psn);
		if (Cnt != 5)
		{
			pszLine[25] = '\0';
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "ProcessMAFline: Unable to parse species sequence params: %s", pszLine);
			return(eBSFerrParse);
		}

		char* pChr = szSpecies;
		char Chr;
		while (Chr = *pChr++)
		{
			switch (Chr) {
			case '.':
			case '_':
				pChr[-1] = '\0';
				strncpy(szChrom, pChr, sizeof(szChrom));
				szChrom[cMaxDatasetSpeciesChrom - 1] = '\0';
				break;
			default:
				continue;
			}
			break;
		}
		szSpecies[cMaxDatasetSpeciesChrom - 1] = '\0';
		pszLine = StripWS(pszLine + Psn);
		AlignLen = (uint32_t)strlen(pszLine);
		if ((Rslt = pProcParams->pAlignFile->AddAlignSeq(szSpecies, szChrom, (int)ChromLen, (int)Start, (int)AlignLen, cStrand, pszLine)) != eBSFSuccess)
		{
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "ProcessMAFline: Unable to add sequence for %s.%s ofs: %u len: %u\n", szSpecies, szChrom, Start, AlignLen);
			while (pProcParams->pAlignFile->NumErrMsgs())
				gDiagnostics.DiagOut(eDLWarn, gszProcName, "ProcessMAFline: Error-> %s", pProcParams->pAlignFile->GetErrMsg());
		}
		break;
	}
	return(Rslt);
}


void
CGenMAFAlgn::RevSeq(int SeqLenInclInDels, char* pszSeq)
{
	int Cnt;
	char tmp;
	char* pRight = pszSeq + SeqLenInclInDels - 1;
	for (Cnt = 0; Cnt < SeqLenInclInDels / 2; Cnt++)
	{
		tmp = *pszSeq;
		*pszSeq++ = *pRight;
		*pRight-- = tmp;
	}
}

// ReverseComplement
bool // Inplace reverse complement pszSeq
CGenMAFAlgn::ReverseComplement(int SeqLenInclInDels, char* pszSeq)
{
	int Cnt;
	RevSeq(SeqLenInclInDels, pszSeq);
	for (Cnt = 0; Cnt < SeqLenInclInDels; Cnt++, pszSeq++)
	{
		switch (*pszSeq) {
		case 'a':
			*pszSeq = 't';
			break;
		case 'A':
			*pszSeq = 'T';
			break;
		case 'c':
			*pszSeq = 'g';
			break;
		case 'C':
			*pszSeq = 'G';
			break;
		case 'g':
			*pszSeq = 'c';
			break;
		case 'G':
			*pszSeq = 'C';
			break;
		case 't':
			*pszSeq = 'a';
			break;
		case 'T':
			*pszSeq = 'A';
			break;
		case 'u':
			*pszSeq = 'a';
			break;
		case 'U':
			*pszSeq = 'A';
			break;
		case '\0':
			return(false);
		default:
			break;
		}
	}
	return(true);
}

//AXT file Structure
//Each alignment block in an axt file contains three lines: a summary line and 2 sequence lines. Blocks are separated from one another by blank lines. 

//1. Summary line 

// 0 chr19 3001012 3001075 chr11 70568380 70568443 - 3500
// The summary line contains chromosomal position and size information about the alignment.
// It consists of 9 required fields: 
// Alignment number -- The alignment numbering starts with 0 [0..numalignments-1]
// Chromosome (primary organism) 
// Alignment start (primary organism) -- The first base is numbered 1 (in MAF co-ordinates are 0 based!). 
// Alignment end (primary organism) -- The end base is included. 
// Chromosome (aligning organism) 
// Alignment start (aligning organism) 
// Alignment end (aligning organism) 
// Strand (aligning organism)
// Blastz score -- Different blastz scoring matrices are used for different organisms
//
// Note:  -- If the strand value is "-", the values of the aligning organism's start and end fields
//           are relative to the reverse-complemented coordinates of its chromosome. 
int
CGenMAFAlgn::ProcessAXTline(int LineLen, tsProcParams * pProcParams)
{
	int Rslt;
	int SeqLenInclInDels;
	int RefChromID;
	int RelChromID;
	int RefChromLen;
	int RelChromLen;

	static bool b1stAlign = true;
	static int32_t AlignNum;
	static int32_t Start1;
	static int32_t End1;
	static int32_t Start2;
	static int32_t End2;
	static int32_t Len1;
	static int32_t Len2;
	static char cStrand;
	static char szChrom1[cMaxDatasetSpeciesChrom];
	static char szChrom2[cMaxDatasetSpeciesChrom];
	static int32_t Score;
	char* pszLine = pProcParams->pszLineBuffer;
	static int BlockID = 0;

	// alignment number as a numeric starts a new block
	if (isdigit(*pszLine))
	{
		// close any currently opened block
		b1stAlign = true;
		sscanf(pszLine, "%*d %s %d %d %s %d %d %c %d", // note: only interested in integral part of score (scores are floating points)
			szChrom1, &Start1, &End1, szChrom2, &Start2, &End2, &cStrand, &Score);
		Len1 = End1 - Start1 + 1;
		Len2 = End2 - Start2 + 1;

		//----- hack to get around the factor that there could be zillions of scafolds in some assemblies which will exceed the limit on number of chromosomes (currently 20000)
#ifdef _WIN32
		if ((szChrom2[0] == 's' || szChrom2[0] == 'S') && !strnicmp("SCAFFOLD", szChrom2, 8)) // if chromosome name starts with "SCAFFOLD" then make all these names the same
#else
		if ((szChrom2[0] == 's' || szChrom2[0] == 'S') && !strncasecmp("SCAFFOLD", szChrom2, 8)) // if chromosome name starts with "SCAFFOLD" then make all these names the same

#endif
			szChrom2[8] = '\0';				// by simply truncating
		// ----- end of hack
		if ((Rslt = pProcParams->pAlignFile->StartAlignBlock(Score)) < eBSFSuccess)	// start a new block
			return(Rslt);

		BlockID = Rslt;
		return(eBSFSuccess);
	}

	// sequence base starts alignment
	Rslt = eBSFSuccess;
	if (IsBase(*pszLine))
	{
		SeqLenInclInDels = (int)strlen(pszLine);
		if (b1stAlign)
		{
			strcpy(pProcParams->pszAXTBuffer, pszLine);
			b1stAlign = false;
			return(eBSFSuccess);
		}
		// not 1st alignment so pszLine pts to 2nd alignment sequence
		RefChromID = m_pAlignValidate->GetChromID(pProcParams->pszRefSpecies, szChrom1);
		RefChromLen = m_pAlignValidate->GetChromLen(RefChromID);
		if (RefChromLen < 0)
			RefChromLen = 0;
		RelChromID = m_pAlignValidate->GetChromID(pProcParams->pszRelSpecies, szChrom2);
		RelChromLen = m_pAlignValidate->GetChromLen(RelChromID);
		if (RelChromLen < 0)
			RelChromLen = 0;

		if (pProcParams->bSwapRefRel)	// need to exchange ref and rel?
		{
			if (cStrand == '-')
			{
				ReverseComplement(SeqLenInclInDels, pszLine);
				ReverseComplement(SeqLenInclInDels, pProcParams->pszAXTBuffer);
			}

			Rslt = pProcParams->pAlignFile->AddAlignSeq(pProcParams->pszRelSpecies, szChrom2, RelChromLen, Start2 - 1, SeqLenInclInDels, '+', pszLine);
			Rslt = pProcParams->pAlignFile->AddAlignSeq(pProcParams->pszRefSpecies, szChrom1, RefChromLen, Start1 - 1, SeqLenInclInDels, cStrand, pProcParams->pszAXTBuffer); // note that AXT have 1..n, we use 0..n locus co-ordinates

		}
		else	// not exchanging....
		{

			Rslt = pProcParams->pAlignFile->AddAlignSeq(pProcParams->pszRefSpecies, szChrom1, RefChromLen, Start1 - 1, SeqLenInclInDels, '+', pProcParams->pszAXTBuffer); // note that AXT have 1..n, we use 0..n locus co-ordinates
			Rslt = pProcParams->pAlignFile->AddAlignSeq(pProcParams->pszRelSpecies, szChrom2, RelChromLen, Start2 - 1, SeqLenInclInDels, cStrand, pszLine);

		}

		if (Rslt < 0)
		{
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "Unable to add AXT Alignment!");
			while (pProcParams->pAlignFile->NumErrMsgs())
				gDiagnostics.DiagOut(eDLWarn, gszProcName, "ProcessAXTline: Error-> %s", pProcParams->pAlignFile->GetErrMsg());
		}
	}
	return(Rslt);
}

// CreateMAlignFile
int		// Create a multialignment file from files (either axt or mfa format) in specified source directory
CGenMAFAlgn::CreateMAlignFile(char* pszChromLens, char* pszSrcDirPath, char* pszDestAlignFile,
	char* pszDescr, char* pszTitle,
	char* pszRefSpecies, char* pszRelSpecies, bool bIsAXT, bool bSwapRefRel)
{
int Rslt;
tsProcParams ProcParams;

memset(&ProcParams, 0, sizeof(tsProcParams));
if ((m_pAlignValidate = new CAlignValidate) == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to instantiate instance of CAlignValidate");
	return(eBSFerrObj);
	}

if ((ProcParams.pAlignFile = new CMAlignFile) == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to instantiate instance of CMAlignFile");
	return(eBSFerrObj);
	}

ProcParams.pRawBuffer = new char[cMARawBuffSize];
if (ProcParams.pRawBuffer == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to allocate memory (%d bytes) for raw buffer", cMARawBuffSize);
	return(eBSFerrMem);
	}

if ((ProcParams.pszLineBuffer = new char[cMALineSize]) == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to allocate memory (%d bytes) for line buffer", cMALineSize);
	delete ProcParams.pRawBuffer;
	return(eBSFerrMem);
	}

if (bIsAXT)
	{
	if ((ProcParams.pszAXTBuffer = new char[cMALineSize]) == nullptr)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to allocate memory (%d bytes) for AXT buffers", cMALineSize);
		delete ProcParams.pAlignFile;
		delete ProcParams.pRawBuffer;
		delete ProcParams.pszLineBuffer;
		return(eBSFerrMem);
		}

	if ((Rslt = m_pAlignValidate->ProcessChroms(pszChromLens)) != eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to process file containing chromosome lengths from %s", pszChromLens);
		delete ProcParams.pAlignFile;
		delete ProcParams.pRawBuffer;
		delete ProcParams.pszLineBuffer;
		if (ProcParams.pszAXTBuffer == nullptr)
			delete ProcParams.pszAXTBuffer;
		return(Rslt);
		}
	}

ProcParams.bSwapRefRel = bSwapRefRel;
ProcParams.bIsAXT = bIsAXT;
ProcParams.pszRefSpecies = pszRefSpecies;
ProcParams.pszRelSpecies = pszRelSpecies;

if ((Rslt = ProcParams.pAlignFile->Open(pszDestAlignFile, eMAPOCreate, pszRefSpecies)) != eBSFSuccess)
	{
	while (ProcParams.pAlignFile->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal, gszProcName, ProcParams.pAlignFile->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Failed to create output file '%s'", pszDestAlignFile);
	}
if (Rslt == eBSFSuccess)
	{
	ProcParams.pAlignFile->SetDescription(pszDescr);
	ProcParams.pAlignFile->SetTitle(pszTitle);
	CSimpleGlob glob(SG_GLOB_FULLSORT);
	if (glob.Add(pszSrcDirPath) >= SG_SUCCESS)
		{
		for (int n = 0; Rslt >= eBSFSuccess && n < glob.FileCount(); ++n)
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "Will process in this order: %d '%s", n + 1, glob.File(n));

		for (int n = 0; Rslt >= eBSFSuccess && n < glob.FileCount(); ++n)
			{
			Rslt = ProcessThisFile(glob.File(n), &ProcParams);
			if (Rslt < eBSFSuccess)
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Errors processing '%s", glob.File(n));
				break;
				}
			}
		}
	else
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to glob '%s", pszSrcDirPath);
		Rslt = eBSFerrOpnFile;	// treat as though unable to open file
		}
	if (Rslt >= eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Completed processing, now commiting file header");
		if ((Rslt = ProcParams.pAlignFile->Close(true)) != eBSFSuccess)
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Errors whilst commiting file header");
		}
	else
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Completed processing but with errors");
}
delete ProcParams.pAlignFile;
m_pAlignValidate->Reset();
delete ProcParams.pRawBuffer;
delete ProcParams.pszLineBuffer;
if (ProcParams.pszAXTBuffer == nullptr)
	delete ProcParams.pszAXTBuffer;
return(Rslt);
}

