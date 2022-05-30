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
#include "seghaplotypes.h"
#include "pangenome.h"


int Process(eModePG PMode,			// processing mode
			char* pszPrefix,		// descriptor prefix
			int BinSizeKbp,			// Wiggle score is number of alignments over these sized bins
			char* pszInFile,		// input fasta or SAM file
			char* pszOutFile);		// output to this file

#ifdef _WIN32
int pangenome(int argc, char *argv[])
{
	// determine my process name
	_splitpath (argv[0], NULL, NULL, gszProcName, NULL);
#else
int
pangenome(int argc, char **argv)
{
	// determine my process name
	CUtility::splitpath ((char *)argv[0], NULL, gszProcName);
#endif
	int iFileLogLevel;			// level of file diagnostics
	int iScreenLogLevel;		// level of file diagnostics
	char szLogFile[_MAX_PATH];	// write diagnostics to this file
	int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure
	int NumberOfProcessors;		// number of installed CPUs

	 eModePG PMode;					// processing mode
	 int BinSizeKbp;	// Wiggle score is number of alignments over these sized bins
	 char szInFile[_MAX_PATH];		// input fasta or SAM file
	 char szPrefix[_MAX_PATH];		// use this prefix
	 char szOutFile[_MAX_PATH];		 // output file

	struct arg_lit *help = arg_lit0 ("h", "help", "print this help and exit");
	struct arg_lit *version = arg_lit0 ("v", "version,ver", "print version information and exit");
	struct arg_int *FileLogLevel = arg_int0 ("f", "FileLogLevel", "<int>", "Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
	struct arg_file *LogFile = arg_file0 ("F", "log", "<file>", "diagnostics log file");

	struct arg_int *pmode = arg_int0 ("m", "mode", "<int>", "processing mode: 0 prefix fasta descriptor, 1 filtering SAM for target prefixes, 2 generate Wiggle all alignments in bin, 3 generate Wiggle unique loci in bin");
	struct arg_int *binsizekbp = arg_int0 ("b", "binsizekbp", "<int>", "if generating Wiggle then maximum Kbp bin size (default 10, range 1..1000");

	struct arg_str *prefix = arg_str1("p","prefix","<str>", "prefix to apply (alphanumeric only, length limited to a max of 10 chars)");
	struct arg_file *infile = arg_file1("i", "in", "<file>", "input file");
	struct arg_file *outfile = arg_file1 ("o", "out", "<file>", "output file");
	struct arg_end *end = arg_end (200);

	void *argtable[] = { help,version,FileLogLevel,LogFile,
						pmode,prefix,binsizekbp,infile,outfile,end };

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

		PMode = pmode->count ? (eModePG)pmode->ival[0] : eMPGDefault;
		if (PMode < eMPGDefault || PMode >= eMPGPlaceHolder)
			{
			gDiagnostics.DiagOut (eDLFatal, gszProcName, "Error: Processing mode '-m%d' specified outside of range %d..%d\n", PMode, eMPGDefault, eMPGPlaceHolder-1);
			exit (1);
			}

		if(PMode >= eMPGWiggleUniqueLoci)
			{
			BinSizeKbp = binsizekbp->count ? binsizekbp->ival[0] : cDfltWiggleBinSize;
			if(BinSizeKbp < cMinWiggleBinSize || BinSizeKbp > cMaxWiggleBinSize)
				{
				gDiagnostics.DiagOut (eDLFatal, gszProcName, "Error: Windowed maximum bin size specified as '-b%d'Kbp , outside of range %d..%d\n", BinSizeKbp, cMinWiggleBinSize,cMaxWiggleBinSize);
				exit (1);
				}
			}
		else
			BinSizeKbp = 0;


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

		szPrefix[0] = '\0';

		if(prefix->count)
			{
			strcpy (szPrefix, prefix->sval[0]);
			CUtility::TrimQuotedWhitespcExtd (szPrefix);
			}
		if (szPrefix[0] == '\0')
			{
			gDiagnostics.DiagOut (eDLFatal, gszProcName, "No prefix specified");
			exit (1);
			}
		if(strlen(szPrefix) > cMaxLenPrefix)
			{
			gDiagnostics.DiagOut (eDLFatal, gszProcName, "Prefix \"%s\" length must <= %d chars",szPrefix,cMaxLenPrefix);
			exit (1);
			}
		// prefix must be alpha-numeric only
		char Chr;
		char *pChr;
		pChr = szPrefix;
		while(Chr = *pChr++)
			{
			if(!isalnum(Chr))
				{
				gDiagnostics.DiagOut (eDLFatal, gszProcName, "Prefix \"%s\" may only contain alpha-numeric chars",szPrefix);
				exit (1);
				}
			}

		strcpy (szInFile, infile->filename[0]);
		CUtility::TrimQuotedWhitespcExtd (szInFile);
		if (szInFile[0] == '\0')
			{
			gDiagnostics.DiagOut (eDLFatal, gszProcName, "No input file specified");
			exit (1);
			}

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
			case eMPGDefault:
				pszDescr = "Prefixing fasta descriptors with specified prefix";
				break;
			case eMPGFilterPrefix:
				pszDescr = "Filter SAM alignments to targets with specified prefix";
				break;

			case eMPGWiggleUniqueLoci:
				pszDescr = "Generate UCSC Wiggle format file using unique loci counts";
				break;

			case eMPGWiggleAll:
				pszDescr = "Generate UCSC Wiggle format file using all alignment counts";
				break;
		}



		gDiagnostics.DiagOutMsgOnly (eDLInfo, "Pangenome processing : '%s'", pszDescr);
		if(BinSizeKbp > 0)
			gDiagnostics.DiagOutMsgOnly (eDLInfo, "Counts accumulated into maximal sized bins of : %dKbp'", BinSizeKbp);	
		if(szPrefix[0] != '\0')
			gDiagnostics.DiagOutMsgOnly (eDLInfo, "Prefix : '%s'", szPrefix);
		gDiagnostics.DiagOutMsgOnly (eDLInfo, "Input file : '%s'", szInFile);
		gDiagnostics.DiagOutMsgOnly (eDLInfo, "Output file : '%s'", szOutFile);


#ifdef _WIN32
		SetPriorityClass (GetCurrentProcess (), BELOW_NORMAL_PRIORITY_CLASS);
#endif
		gStopWatch.Start ();
		Rslt = 0;
		Rslt = Process (PMode,			// processing mode
						szPrefix,		// descriptor prefix
						BinSizeKbp,		// bin size in Kbp
						szInFile,		// input fasta or SAM file
						szOutFile);		// output to this file
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


CPangenome::CPangenome()
{
m_pInBuffer = NULL;
m_pOutBuffer = NULL;
m_pSAMfile=NULL;
m_pBAMalignment =NULL;
m_pSAMloci = NULL;
m_hInFile = -1;			// input file handle
m_hOutFile = -1;			// output file handle
Reset();
}

CPangenome::~CPangenome()
{
if(m_hInFile != -1)
	close(m_hInFile);
if(m_hOutFile != -1)
	close(m_hOutFile);
if(m_pInBuffer != NULL)
	delete []m_pInBuffer;
if(m_pOutBuffer != NULL)
	delete []m_pOutBuffer;
if(m_pSAMfile!=NULL)
	delete m_pSAMfile;

if(m_pBAMalignment !=NULL)
	delete m_pBAMalignment;

if (m_pSAMloci != NULL)
	{
#ifdef _WIN32
	free(m_pSAMloci);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (m_pSAMloci != MAP_FAILED)
		munmap(m_pSAMloci, m_AllocdSAMlociMem);
#endif
	}
}

void
CPangenome::Reset(void)	// resets class instance state back to that immediately following instantiation
{
if(m_hInFile != -1)
	{
	close(m_hInFile);
	m_hInFile = -1;
	}
if(m_hOutFile != -1)
	{
	close(m_hOutFile);
	m_hOutFile = -1;
	}

if(m_pInBuffer != NULL)
	{
	delete []m_pInBuffer;
	m_pInBuffer = NULL;
	}

if(m_pOutBuffer != NULL)
	{
	delete []m_pOutBuffer;
	m_pOutBuffer = NULL;
	}

if (m_pSAMloci != NULL)
	{
#ifdef _WIN32
	free(m_pSAMloci);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (m_pSAMloci != MAP_FAILED)
		munmap(m_pSAMloci, m_AllocdSAMlociMem);
#endif
	m_pSAMloci = NULL;
	}

m_InBuffIdx = 0;
m_AllocInBuff = 0;

m_OutBuffIdx = 0;	
m_AllocOutBuff = 0;

m_CurNumSAMloci = 0;
m_AllocdSAMloci = 0;
m_AllocdSAMlociMem = 0;

m_LASeqNameID = 0;
m_NumSeqNames = 0;
m_NxtszSeqNameIdx = 0;
m_szSeqNames[0] = '\0';
m_szSeqNameIdx[0] = 0;
}

int
CPangenome::PrefixFasta(char* pszPrefix,		// descriptor prefix
			char* pszInFile,		// input fasta or SAM file
			char* pszOutFile)		// output to this file
{
bool bNL;

int NumRead;
int PrefixLen;
uint32_t Idx;
uint8_t InChr;
uint8_t *pInChr;
uint8_t *pOutChr;

Reset();

#ifdef _WIN32
m_hInFile = open(pszInFile, O_READSEQ );		// file access is normally sequential..
#else
m_hInFile = open64(pszInFile, O_READSEQ );		// file access is normally sequential..
#endif
if(m_hInFile == -1)							// check if file open succeeded
	{
	gDiagnostics.DiagOut (eDLFatal, gszProcName, "Unable to open input file '%s' : %s",pszInFile,strerror(errno));
	Reset();
	return(eBSFerrOpnFile);
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

if((m_pInBuffer = new uint8_t[cAllocPGBuffInSize]) == NULL)
	{
	gDiagnostics.DiagOut (eDLFatal, gszProcName, "Unable to allocate memory for input file buffering");
	Reset();
	return(eBSFerrMem);
	}
m_AllocInBuff = cAllocPGBuffInSize;

if((m_pOutBuffer = new uint8_t[cAllocPGBuffOutSize]) == NULL)
	{
	gDiagnostics.DiagOut (eDLFatal, gszProcName, "Unable to allocate memory for output file buffering");
	delete []m_pInBuffer;
	Reset();
	return(eBSFerrMem);
	}
m_AllocOutBuff = cAllocPGBuffOutSize;

PrefixLen = (int)strlen(pszPrefix);
m_InBuffIdx = 0;
m_OutBuffIdx = 0;
pOutChr = m_pOutBuffer;
bNL = true; // 1st line in fasta may be a descriptor with no preceding NL so pretend there was a preceding NL  
while((NumRead = (int)read(m_hInFile, m_pInBuffer, (int)m_AllocInBuff)) > 0)
	{
	m_InBuffIdx = NumRead;
	pInChr = m_pInBuffer;
	for(Idx = 0; Idx < m_InBuffIdx; Idx++)
		{
		if((InChr = *pInChr++) == '\0')
			break;
		switch(InChr) {
			case '\0': case '\t': case '\n': case '\r':	// whitespace or line endings are accepted
				break;
			default:
				if(InChr >= 0x20 && InChr <= 0x7f)
					break;
			gDiagnostics.DiagOut (eDLFatal, gszProcName, "Expected Fasta '%s' to only contain ascii chars",pszInFile);
			Reset();
			return(eBSFerrFastqChr);
			}
		*pOutChr++ = InChr;
		m_OutBuffIdx++;
		if(InChr == '\n')
			bNL = true;
		else
			{
			if(bNL == true && InChr == '>')
				{
				strcpy((char *)pOutChr,pszPrefix);
				pOutChr += PrefixLen;
				*pOutChr++ = cTagSHTerm1;
				*pOutChr++ = cTagSHTerm2;
				m_OutBuffIdx += PrefixLen+1;
				}
			bNL = false;
			}
		
		if(m_OutBuffIdx >= (m_AllocOutBuff - PrefixLen - 2))
			{
			CUtility::RetryWrites(m_hOutFile, m_pOutBuffer, m_OutBuffIdx);	
			m_OutBuffIdx = 0;
			pOutChr = m_pOutBuffer;
			}
		}
	if(InChr == '\0')
		{
		*pOutChr++ = '\0';
		m_OutBuffIdx++;
		break;
		}
	}

if(m_OutBuffIdx)
	CUtility::RetryWrites(m_hOutFile, m_pOutBuffer, m_OutBuffIdx);	

// commit output file
#ifdef _WIN32
_commit(m_hOutFile);
#else
fsync(m_hOutFile);
#endif
Reset();
return(eBSFSuccess);
}


int
CPangenome::FilterSAM(char* pszPrefix,			// target prefix used for filtering from
			char* pszInFile,		// input SAM file with matching alignments
			char* pszOutFile)		// output to this SAM file
{
bool bNL;
int NumRead;
uint32_t NumPrefixed;
int PrefixLen;
uint32_t Idx;
uint8_t InChr;
uint8_t *pInChr;
uint8_t *pOutChr;

char szAssembRef[120];
char szSeqRef[120];
uint32_t SeqRefLen;
char szPrefix[cMaxSHLenPrefix + 3];		// allowing for 2 terminating chars plus '\0'

Reset();

strcpy(szPrefix,pszPrefix);
PrefixLen = (int)strlen(szPrefix);
szPrefix[PrefixLen++] = cTagSHTerm1;
szPrefix[PrefixLen++] = cTagSHTerm2;
szPrefix[PrefixLen] = '\0';

#ifdef _WIN32
m_hInFile = open(pszInFile, O_READSEQ );		// file access is normally sequential..
#else
m_hInFile = open64(pszInFile, O_READSEQ );		// file access is normally sequential..
#endif
if(m_hInFile == -1)							// check if file open succeeded
	{
	gDiagnostics.DiagOut (eDLFatal, gszProcName, "Unable to open input file '%s' : %s",pszInFile,strerror(errno));
	Reset();
	return(eBSFerrOpnFile);
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

if((m_pInBuffer = new uint8_t[cAllocPGBuffInSize]) == NULL)
	{
	gDiagnostics.DiagOut (eDLFatal, gszProcName, "Unable to allocate memory for input file buffering");
	Reset();
	return(eBSFerrMem);
	}
m_AllocInBuff = cAllocPGBuffInSize;

if((m_pOutBuffer = new uint8_t[cAllocPGBuffOutSize]) == NULL)
	{
	gDiagnostics.DiagOut (eDLFatal, gszProcName, "Unable to allocate memory for output file buffering");
	Reset();
	return(eBSFerrMem);
	}
m_AllocOutBuff = cAllocPGBuffOutSize;

m_InBuffIdx = 0;
m_OutBuffIdx = 0;
pOutChr = m_pOutBuffer;
bNL = false;
uint32_t RecordStartIdx;
uint8_t *pChr;
uint32_t SAMreclen = 0;
RecordStartIdx = 0;
NumPrefixed = 0;
while((NumRead = (int)read(m_hInFile, m_pInBuffer, (int)m_AllocInBuff)) > 0)
	{
	m_InBuffIdx = NumRead;
	pInChr = m_pInBuffer;
	for(Idx = 0; Idx < m_InBuffIdx; Idx++)
		{
		InChr = *pInChr++;
		switch(InChr) {
			case '\0': case '\t': case '\n': case '\r':	// whitespace or line endings are accepted
				break;
			default:
				if(InChr >= 0x20 && InChr <= 0x7f)
					break;
			gDiagnostics.DiagOut (eDLFatal, gszProcName, "Expected SAM '%s' to only contain ASCII chars",pszInFile);
			Reset();
			return(eBSFerrFastqChr);
			}

		*pOutChr++ = InChr;
		m_OutBuffIdx++;
		SAMreclen++;
		if(SAMreclen > 20000)
			{
			gDiagnostics.DiagOut (eDLFatal, gszProcName, "Expected SAM '%s' to contain individual records of less than 20K chars in length",pszInFile);
			Reset();
			return(eBSFerrFastqChr);
			}
		if(InChr == '\0' || InChr == '\n')	// have a complete SAM record?
			{
			pChr = &m_pOutBuffer[RecordStartIdx];
			if(*pChr == '@')		// is SAM record a header?
				{
				if(pChr[1] == 'S' && pChr[2] == 'Q')	// only retaining sequence names with specified prefix
					{
					if(sscanf((char*)&pChr[3], "\tAS:%s\tSN:%s\tLN:%u", szAssembRef, szSeqRef, &SeqRefLen) == 3)
						{
						if(strnicmp(szSeqRef, szPrefix, PrefixLen))
							{
							SAMreclen = 0;
							m_OutBuffIdx = RecordStartIdx;
							pOutChr = &m_pOutBuffer[RecordStartIdx];
							continue;
							}
					// sequence name has specified prefix, strip prefix off and accept this header record
						sprintf((char*)&pChr[3], "\tAS:%s\tSN:%s\t:LN:%u\n", szAssembRef, &szSeqRef[PrefixLen], SeqRefLen);
						m_OutBuffIdx -= (PrefixLen - 1);
						}
					}
				SAMreclen = 0;
				RecordStartIdx = m_OutBuffIdx;
				pOutChr = &m_pOutBuffer[RecordStartIdx];
				}
			else    // not a header record so assume it's an alignment
				{
				if(SAMreclen <= 20)			// could be a near empty record in which case slough this record
					{
					SAMreclen = 0;
					m_OutBuffIdx = RecordStartIdx;
					pOutChr = &m_pOutBuffer[RecordStartIdx];
					continue;
					}

			// skip over 1st 2 tabs until target sequence name
				pChr = &m_pOutBuffer[RecordStartIdx];
				int NumTabs = 0;
				do{
					if(*pChr++ == '\t')
						NumTabs++;
				}
				while(--SAMreclen && NumTabs < 2);

				if(NumTabs != 2 || strnicmp((char*)pChr, szPrefix, PrefixLen)){
					SAMreclen = 0;
					m_OutBuffIdx = RecordStartIdx;
					pOutChr = &m_pOutBuffer[RecordStartIdx];
					continue;
				}

			// have an alignment target prefixed by szPrefix
				NumPrefixed++;
				// strip prefix off
				do{
					*pChr++ = (InChr = pChr[PrefixLen]);
				}
				while(!(InChr == '\0' || InChr == '\n'));
				m_OutBuffIdx -= PrefixLen;
				RecordStartIdx = m_OutBuffIdx;
				pOutChr = &m_pOutBuffer[RecordStartIdx];
				SAMreclen = 0;
				}
			}

			if(m_OutBuffIdx == RecordStartIdx && m_OutBuffIdx >= (m_AllocOutBuff - 100000))	// assumes max length SAM record of 100K chars
			{
				CUtility::RetryWrites(m_hOutFile, m_pOutBuffer, m_OutBuffIdx);
				m_OutBuffIdx = 0;
				RecordStartIdx = 0;
				pOutChr = m_pOutBuffer;
				SAMreclen = 0;
			}
		}
		if(InChr == '\0'){
			*pOutChr++ = '\0';
			RecordStartIdx++;
			break;
		}
	}

	if(RecordStartIdx)
		CUtility::RetryWrites(m_hOutFile, m_pOutBuffer, RecordStartIdx);

	// commit output file
#ifdef _WIN32
	_commit(m_hOutFile);
#else
	fsync(m_hOutFile);
#endif
	Reset();
	szPrefix[PrefixLen - 2] = '\0';
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Filtered %u alignments with target prefixed by '%s' into '%s'", NumPrefixed, szPrefix, pszOutFile);
	return(eBSFSuccess);
}


int Process(eModePG PMode,			// processing mode
			char* pszPrefix,		// descriptor prefix
			int BinSizeKbp,			// bin size in Kbp
			char* pszInFile,		// input fasta or SAM file
			char* pszOutFile)		// output to this file
{
int Rslt;
CPangenome* pPangenome;

if((pPangenome = new CPangenome) == NULL)
	{
	gDiagnostics.DiagOut (eDLFatal, gszProcName, "Unable to instantiate instance of CPangenome");
	return(eBSFerrObj);
	}
Rslt = pPangenome->Process(PMode,pszPrefix,BinSizeKbp,pszInFile,pszOutFile);
delete pPangenome;
return(Rslt);
}

int 
CPangenome::Process(eModePG PMode,			// processing mode
			char* pszPrefix,		// descriptor prefix
			int BinSizeKbp,			// bin size in Kbp
			char* pszInFile,		// input fasta or SAM file
			char* pszOutFile)		// output to this file
{
int NumRead;
uint8_t Chr;
int Rslt = eBSFerrParams;
uint8_t *m_pInBuffer;

if(PMode < eMPGWiggleUniqueLoci)		// wiggles input can be BAM (binary) so can't just check for ascii format input
	{
	if((m_pInBuffer = new uint8_t[cAllocAsciiChkSize]) == NULL)
		{
		gDiagnostics.DiagOut (eDLFatal, gszProcName, "Unable to allocate memory for input file ASCII only check");
		return(eBSFerrMem);
		}
	m_AllocInBuff = cAllocAsciiChkSize;

	// input file currently expected to be ascii only, binary format (BAM, gz'd etc., currently can't be processed)
	// read in up to 4M bytes and if any not ascii then report to user
	#ifdef _WIN32
	m_hInFile = open(pszInFile, O_READSEQ );		// file access is normally sequential..
	#else
	m_hInFile = open64(pszInFile, O_READSEQ );		// file access is normally sequential..
	#endif
	if(m_hInFile == -1)							// check if file open succeeded
		{
		gDiagnostics.DiagOut (eDLFatal, gszProcName, "Unable to open input file '%s' : %s",pszInFile,strerror(errno));
		Reset();
		return(eBSFerrOpnFile);
		}

	if((NumRead = (int)read(m_hInFile, m_pInBuffer, (int)cAllocAsciiChkSize)) <= 1000)	// expecting input file to be sized at least 1000
		{
		gDiagnostics.DiagOut (eDLFatal, gszProcName, "Expected '%s to contain at least 1000 chars, read() returned %d",pszInFile,NumRead);
		Reset();
		return(eBSFerrOpnFile);
		}
	close(m_hInFile);
	m_hInFile = -1;
	while(--NumRead)
		{
		switch(Chr = m_pInBuffer[NumRead]) {
			case '\0':		// can accept as long as not within 1st 1000 chars, some processes place a NULL char at end of ascii file
				if(NumRead < 1000)
					break;
				continue;
			case ' ': case '\t': case '\n': case '\r':	// whitespace or line endings are accepted
				continue;
			default:
				if(Chr >= 0x20 && Chr <= 0x7f)
					continue;
				break;
			}

		gDiagnostics.DiagOut (eDLFatal, gszProcName, "Expected '%s' to contain only ascii chars (BAM or gz'd input files not currently supported)",pszInFile);
		Reset();
		return(eBSFerrOpnFile);
		}
	delete []m_pInBuffer;
	m_pInBuffer = NULL;
	m_AllocInBuff = 0;
	}

switch(PMode){
	case eMPGDefault:				// prefixing fasta descriptors
		Rslt = PrefixFasta(pszPrefix,		// descriptor prefix
				pszInFile,		// input fasta or SAM file
				pszOutFile);		// output to this file
		break;

	case eMPGFilterPrefix:
		Rslt = FilterSAM(pszPrefix,		// target prefix used for filtering from
						pszInFile,		// input SAM file with matching alignments
						pszOutFile);	// output to this SAM file
		break;

	case eMPGWiggleAll:			// generate UCSC Wiggle format file using all alignments in bins
	case eMPGWiggleUniqueLoci:			// generate UCSC Wiggle format file using only unique alignments in bins
		Rslt = GenBinnedWiggle(PMode,			// processing mode
			pszPrefix,		// target prefix used for filtering from
			BinSizeKbp,					// Wiggle score is number of alignments over this sized bin
			pszInFile,					// alignments are in this SAM/BAM file 
			pszOutFile);				// write out Wiggle to this file
		break;
	}

return(Rslt);
}


// generate UCSC Wiggle file from SAM/BAM alignment input file
// wiggle file contains smoothed alignment density along each alignment targeted sequence
int	
CPangenome::GenBinnedWiggle(eModePG PMode,			// processing mode
			char *pszName,	// wiggle track name
			uint32_t BinSizeKbp,	// Wiggle score is number of alignments over this sized bins
			char* pszInFile,		// alignments are in this SAM/BAM file 
			char* pszOutFile)		// write out Wiggle to this file
{
teBSFrsltCodes Rslt;

int LineLen;
char *pszTxt;
int WindowSize = BinSizeKbp * 1000;
int TargID;
tsSAMloci *pSAMloci;

// open SAM for reading
if(pszInFile == NULL || *pszInFile == '\0')
	return(eBSFerrParams);

if((m_pSAMfile = new CSAMfile) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenBinnedWiggle: Unable to instantiate class CSAMfile");
	return(eBSFerrInternal);
	}

if((m_pBAMalignment = new tsBAMalign) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenBinnedWiggle: Unable to instantiate tsBAMalign");
	Reset();
	return(eBSFerrInternal);
	}

if((m_pInBuffer = new uint8_t[cMaxReadLen * 3]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenBinnedWiggle: Unable to allocate buffers");
	Reset();
	return(eBSFerrMem);
	}
m_AllocInBuff = cMaxReadLen * 3;

if((m_pOutBuffer = new uint8_t[cAllocPGBuffOutSize]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenBinnedWiggle: Unable to allocate buffers");
	Reset();
	return(eBSFerrMem);
	}
m_AllocOutBuff = cAllocPGBuffOutSize;

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

m_OutBuffIdx = sprintf((char *)m_pOutBuffer,"name=\"Coverage %s\" type=wiggle_0\n",pszName);

m_AllocdSAMlociMem = (size_t)cAllocNumSAMloci * sizeof(tsSAMloci);
#ifdef _WIN32
m_pSAMloci = (tsSAMloci*)malloc(m_AllocdSAMlociMem);
if (m_pSAMloci == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenBinnedWiggle: Memory allocation of %zd bytes failed", (int64_t)m_AllocdSAMlociMem);
	m_AllocdSAMlociMem = 0;
	Reset();
	return(eBSFerrMem);
	}
#else
// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
m_pSAMloci = (tsSAMloci*)mmap(NULL, m_AllocdSAMlociMem, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
if (m_pSAMloci == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenBinnedWiggle: Memory allocation of %zd bytes through mmap()  failed", (int64_t)m_AllocdSAMlociMem, strerror(errno));
	m_pSAMloci = NULL;
	m_AllocdSAMlociMem = 0;
	Reset();
	return(eBSFerrMem);
	}
#endif
m_AllocdSAMloci = cAllocNumSAMloci;

if((Rslt = (teBSFrsltCodes)m_pSAMfile->Open(pszInFile)) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenBinnedWiggle: Unable to open SAM/BAM format file %s",pszInFile);
	Reset();
	return((teBSFrsltCodes)Rslt);
	}

uint32_t NumMissingFeatures = 0;
size_t NumAlignmentsProc = 0;
size_t NumAcceptedAlignments = 0;
uint32_t NumUnmapped = 0;

time_t Then = time(NULL);
time_t Now;
while(Rslt >= eBSFSuccess && (LineLen = m_pSAMfile->GetNxtSAMline((char *)m_pInBuffer)) > 0)
	{
	m_pInBuffer[LineLen] = '\0';
	pszTxt = CUtility::TrimWhitespc((char *)m_pInBuffer);
	if (*pszTxt == '\0')			// simply slough lines which are just whitespace
		continue;
	if (*pszTxt == '@')				// only interested in lines with putative alignments
		continue;

	NumAlignmentsProc += 1;
	if (!(NumAlignmentsProc % 100000) || NumAlignmentsProc == 1)
		{
		Now = time(NULL);
		if ((Now - Then) >= 60)
			{
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "Accepted %zd SAM/BAM alignments", NumAcceptedAlignments);
			Then += 60;
			}
		}

	// primary interest is in the reference chrom name, startloci, length
	if ((Rslt = (teBSFrsltCodes)m_pSAMfile->ParseSAM2BAMalign(pszTxt, m_pBAMalignment, NULL,true)) < eBSFSuccess)
		{
		if (Rslt == eBSFerrFeature)	// not too worried if aligned to feature is missing as some SAMs are missing header features
			{
			NumMissingFeatures++;
			Rslt = eBSFSuccess;
			continue;
			}
		break;
		}

	// check if read has been mapped, if not then slough ...
	if (m_pBAMalignment->refID == -1 || (m_pBAMalignment->flag_nc >> 16) & 0x04)
		{
		NumUnmapped++;
		continue;
		}

	// can now access the alignment loci info
	// hard or soft clipping is currently of no interest as using a large sliding bin and counting loci in that bin
	if(m_CurNumSAMloci > m_AllocdSAMloci) // needing to realloc?
		{
		size_t memreq = m_AllocdSAMlociMem + ((size_t)cAllocNumSAMloci * sizeof(tsSAMloci));
		uint8_t* pTmp;
#ifdef _WIN32
		pTmp = (uint8_t*)realloc(m_pSAMloci, memreq);
#else
		pTmp = (uint8_t*)mremap(m_pSAMloci, m_AllocdSAMlociMem, memreq, MREMAP_MAYMOVE);
		if (pTmp == MAP_FAILED)
			pTmp = NULL;
#endif
		if (pTmp == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "GenBinnedWiggle: Memory re-allocation to %zd bytes - %s", (int64_t)(memreq), strerror(errno));
			Reset();
			return(eBSFerrMem);
			}
		m_pSAMloci = (tsSAMloci *)pTmp;
		m_AllocdSAMlociMem = memreq;
		m_AllocdSAMloci += (size_t)cAllocNumSAMloci;
		}
	pSAMloci = &m_pSAMloci[m_CurNumSAMloci++];
	TargID = AddTargSeqName(m_pBAMalignment->szRefSeqName);
	pSAMloci->TargID = TargID;
	pSAMloci->TargLoci = m_pBAMalignment->pos;
	pSAMloci->Cnt = 1;
	NumAcceptedAlignments++;
	}

if(NumAcceptedAlignments == 0)		// ugh, no alignments!
	{
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenBinnedWiggle: No alignments accepted for processing!");
	Reset();
	return(eBSFSuccess);			// not an error!
	}
else
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenBinnedWiggle: Accepted %zd total alignments for binning",NumAcceptedAlignments);


// sort SAMloci by TargID.TargLoci ascending
if(m_CurNumSAMloci > 1)
	qsort(m_pSAMloci, m_CurNumSAMloci, sizeof(tsSAMloci), SortSAMTargLoci);

uint32_t WindowCnts;
uint32_t NumUniqueLoci;
size_t LociIdx;
tsSAMloci *pUniqueSamLoci;
tsSAMloci *pWinStartSAMloci;

// could be multiple alignments to same loci so reduce such that each loci is unique with a count of alignments to that loci unless unique loci in which case only 1 count attributed
pSAMloci = m_pSAMloci;
pUniqueSamLoci = pSAMloci++;
pUniqueSamLoci->Cnt = 1;
NumUniqueLoci = 1;
for(LociIdx = 1; LociIdx < m_CurNumSAMloci; LociIdx++, pSAMloci++)
	{
	if(pSAMloci->TargID == pUniqueSamLoci->TargID && pSAMloci->TargLoci == pUniqueSamLoci->TargLoci)
		{
		if(PMode == eMPGWiggleAll) // if all to be counted then counting multiple alignments to this loci
			pUniqueSamLoci->Cnt++; 
		}
	else // onto another loci
		{
		pUniqueSamLoci++;
		*pUniqueSamLoci = *pSAMloci;
		pUniqueSamLoci->Cnt = 1;
		NumUniqueLoci++;
		}
	}
m_CurNumSAMloci = NumUniqueLoci;
gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenBinnedWiggle: Alignments were to %zd unique loci",m_CurNumSAMloci);

// iterate over all loci and accumulate counts within each bin sized to be WindowSize bp
uint32_t BinStart;
uint32_t BinEnd;
WindowCnts = 0;
BinStart = (m_pSAMloci->TargLoci / (uint32_t)WindowSize) * (uint32_t)WindowSize;
BinEnd = BinStart + (uint32_t)WindowSize - 1; 
pSAMloci = m_pSAMloci;
pWinStartSAMloci = pSAMloci;
for(LociIdx = 0; LociIdx < m_CurNumSAMloci; LociIdx++,pSAMloci++)
	{
	if(pSAMloci->TargID == pWinStartSAMloci->TargID && pSAMloci->TargLoci <= BinEnd)
		WindowCnts += pSAMloci->Cnt;
	else // outside of an existing bin into which counts were being accumulated
		{
		if(pSAMloci->TargID != pWinStartSAMloci->TargID)	// special case if new loci on a different chromosome to prev as last bin needs to be truncated to last known loci on previous chromosome
			BinEnd = pUniqueSamLoci->TargLoci;
		GenBinnedCoverage(pWinStartSAMloci->TargID,BinStart,BinEnd,WindowCnts);
		WindowCnts = 0;

		BinStart = (pSAMloci->TargLoci / (uint32_t)WindowSize) * (uint32_t)WindowSize;
		BinEnd = BinStart + (uint32_t)WindowSize - 1; 

		pWinStartSAMloci = pSAMloci;
		}

	pUniqueSamLoci = pSAMloci;		// recording last loci known to be within WindowSize 
	}

if(WindowCnts)
	GenBinnedCoverage(pWinStartSAMloci->TargID,BinStart,pUniqueSamLoci->TargLoci,WindowCnts);

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
Reset();
return(eBSFSuccess);
}

int		// bin score assigned
CPangenome::GenBinnedCoverage(int TargID,	// coverage is on this targeted chrom/seq
				uint32_t BinStart,		// coverage starts at this loci inclusive
				uint32_t BinEnd,		// coverage ends at this loci inclusive
				uint32_t BinCnts)		// total number of coverage counts
{
uint32_t BinLen;
BinLen = 1 + BinEnd - BinStart; // inclusive!
if((m_OutBuffIdx + 1000) > m_AllocOutBuff)
	{
	CUtility::RetryWrites(m_hOutFile,m_pOutBuffer,m_OutBuffIdx);
	m_OutBuffIdx = 0;
	}
m_OutBuffIdx += sprintf((char *)&m_pOutBuffer[m_OutBuffIdx],"variableStep chrom=%s span=%d\n%d %d\n",LocateTargSeqName(TargID),BinLen,BinStart + 1,BinCnts); // Wiggle uses 1-start coordinate system

return((int)BinCnts);
}

int		// returned sequence name identifier, < 1 if unable to accept this chromosome name
CPangenome::AddTargSeqName(char* pszSeqName) // associate unique identifier with this sequence name
{
int SeqNameIdx;
int SeqNameLen;
char *pszLAname;

// with any luck the sequence name will be same as the last accessed
if((pszLAname = LocateTargSeqName(m_LASeqNameID)) != NULL)
	if(!stricmp(pszSeqName,pszLAname))
		return(m_LASeqNameID);


// iterate over all known sequence names in case this name to add is a duplicate
for(SeqNameIdx = 0; SeqNameIdx < m_NumSeqNames; SeqNameIdx++)
	if(!stricmp(pszSeqName, &m_szSeqNames[m_szSeqNameIdx[SeqNameIdx]]))
		{
		m_LASeqNameID = SeqNameIdx + 1;
		return(m_LASeqNameID);
		}

// sequence name is not a duplicate
SeqNameLen = (int)strlen(pszSeqName);
if((m_NxtszSeqNameIdx + SeqNameLen + 1) > (int)sizeof(m_szSeqNames))
	return(eBSFerrMaxEntries);
if(m_NumSeqNames == cMaxSeqNames)
	return(eBSFerrMaxEntries);

m_szSeqNameIdx[m_NumSeqNames++] = m_NxtszSeqNameIdx;
strcpy(&m_szSeqNames[m_NxtszSeqNameIdx], pszSeqName);
m_NxtszSeqNameIdx += SeqNameLen + 1;
m_LASeqNameID = m_NumSeqNames;
return(m_NumSeqNames);
}

char*							// returned sequence name
CPangenome::LocateTargSeqName(int SeqNameID)	// identifier returned by call to AddTargSeqName
{
if(SeqNameID < 1 || SeqNameID > m_NumSeqNames)
	return(NULL);
return(&m_szSeqNames[m_szSeqNameIdx[SeqNameID-1]]);
}

// SortSAMTargLoci
// Sort m_pSAMloci by ascending Targ.Loci
int
CPangenome::SortSAMTargLoci(const void* arg1, const void* arg2)
{
tsSAMloci* pEl1 = (tsSAMloci*)arg1;
tsSAMloci* pEl2 = (tsSAMloci*)arg2;

if(pEl1->TargID < pEl2->TargID)
	return(-1);
if(pEl1->TargID > pEl2->TargID)
	return(1);

if(pEl1->TargLoci < pEl2->TargLoci)
	return(-1);
if(pEl1->TargLoci > pEl2->TargLoci)
	return(1);
return(0);
}

