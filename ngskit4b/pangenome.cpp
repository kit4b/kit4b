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
#include "pangenome.h"


int Process(eModePG PMode,			// processing mode
			char* pszPrefix,		// descriptor prefix
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
	 char szInFile[_MAX_PATH];		// input fasta or SAM file
	 char szPrefix[_MAX_PATH];		// use this prefix
	 char szOutFile[_MAX_PATH];		 // output file

	struct arg_lit *help = arg_lit0 ("h", "help", "print this help and exit");
	struct arg_lit *version = arg_lit0 ("v", "version,ver", "print version information and exit");
	struct arg_int *FileLogLevel = arg_int0 ("f", "FileLogLevel", "<int>", "Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
	struct arg_file *LogFile = arg_file0 ("F", "log", "<file>", "diagnostics log file");

	struct arg_int *pmode = arg_int0 ("m", "mode", "<int>", "processing mode: 0 prefix fasta descriptor, 1 filtering SAM for target prefixes");
	struct arg_str *prefix = arg_str1("p","prefix","<str>", "prefix to apply");
	struct arg_file *infile = arg_file1("i", "in", "<file>", "input file");
	struct arg_file *outfile = arg_file1 ("o", "out", "<file>", "output file");
	struct arg_end *end = arg_end (200);

	void *argtable[] = { help,version,FileLogLevel,LogFile,
						pmode,prefix,infile,outfile,end };

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

		strcpy (szPrefix, prefix->sval[0]);
		CUtility::TrimQuotedWhitespcExtd (szPrefix);
		if (szPrefix[0] == '\0')
			{
			gDiagnostics.DiagOut (eDLFatal, gszProcName, "No prefix specified");
			exit (1);
			}
		if(strlen(szPrefix) > 6)
			{
			gDiagnostics.DiagOut (eDLFatal, gszProcName, "Prefix \"%s\" must <= 6 in length",szPrefix);
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
			pChr[-1] = tolower(Chr);
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
		}



		gDiagnostics.DiagOutMsgOnly (eDLInfo, "Pangenome processing : '%s'", pszDescr);
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


int
PrefixFasta(char* pszPrefix,		// descriptor prefix
			char* pszInFile,		// input fasta or SAM file
			char* pszOutFile)		// output to this file
{
bool bNL;
int m_hInFile;			// input file handle
int m_hOutFile;			// output file handle
int NumRead;
int PrefixLen;
uint32_t Idx;
uint8_t InChr;
uint8_t *pInChr;
uint8_t *pOutChr;
uint32_t m_InBuffIdx;	// currently buffering this many input bytes
size_t m_AllocInBuff;	// m_pInBuffer allocated to hold this many input bytes
uint8_t *m_pInBuffer;	// allocated for buffering input

uint32_t m_OutBuffIdx;	// currently buffering this many output bytes
size_t m_AllocOutBuff;	// m_pOutBuffer allocated to hold this many input bytes
uint8_t *m_pOutBuffer;	// allocated for buffering output

#ifdef _WIN32
m_hInFile = open(pszInFile, O_READSEQ );		// file access is normally sequential..
#else
m_hInFile = open64(pszInFile, O_READSEQ );		// file access is normally sequential..
#endif
if(m_hInFile == -1)							// check if file open succeeded
	{
	gDiagnostics.DiagOut (eDLFatal, gszProcName, "Unable to open input file '%s' : %s",pszInFile,strerror(errno));
	return(eBSFerrOpnFile);
	}

#ifdef _WIN32
m_hOutFile = open(pszOutFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE));
#else
if((m_hOutFile = open64(pszOutFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE)) != -1)
	if(ftruncate(m_hOutFile,0)!=0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create/truncate %s - %s",pszOutFile,strerror(errno));
		close(m_hInFile);
		return(eBSFerrCreateFile);
		}
#endif
if(m_hOutFile < 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create/truncate %s - %s",pszOutFile,strerror(errno));
	close(m_hInFile);
	return(eBSFerrCreateFile);
	}

if((m_pInBuffer = new uint8_t[cAllocPGBuffInSize]) == NULL)
	{
	gDiagnostics.DiagOut (eDLFatal, gszProcName, "Unable to allocate memory for input file buffering");
	close(m_hInFile);
	close(m_hOutFile);
	return(eBSFerrMem);
	}
m_AllocInBuff = cAllocPGBuffInSize;

if((m_pOutBuffer = new uint8_t[cAllocPGBuffOutSize]) == NULL)
	{
	gDiagnostics.DiagOut (eDLFatal, gszProcName, "Unable to allocate memory for output file buffering");
	delete []m_pInBuffer;
	close(m_hInFile);
	close(m_hOutFile);
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
				*pOutChr++ = '|';
				m_OutBuffIdx += PrefixLen+1;
				}
			bNL = false;
			}
		
		if(m_OutBuffIdx >= (m_AllocOutBuff - PrefixLen - 2))
			{
			CUtility::SafeWrite(m_hOutFile, m_pOutBuffer, m_OutBuffIdx);	
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
	CUtility::SafeWrite(m_hOutFile, m_pOutBuffer, m_OutBuffIdx);	

// close input file
close(m_hInFile);
// commit and close output file
#ifdef _WIN32
_commit(m_hOutFile);
#else
fsync(m_hOutFile);
#endif
close(m_hOutFile);
delete []m_pInBuffer;
delete []m_pOutBuffer;
return(eBSFSuccess);
}


int
FilterSAM(char* pszPrefix,			// target prefix used for filtering from
			char* pszInFile,		// input SAM file with matching alignments
			char* pszOutFile)		// output to this SAM file
{
bool bNL;
int m_hInFile;			// input file handle
int m_hOutFile;			// output file handle
int NumRead;
int PrefixLen;
uint32_t Idx;
uint8_t InChr;
uint8_t *pInChr;
uint8_t *pOutChr;
uint32_t m_InBuffIdx;	// currently buffering this many input bytes
size_t m_AllocInBuff;	// m_pInBuffer allocated to hold this many input bytes
uint8_t *m_pInBuffer;	// allocated for buffering input

uint32_t m_OutBuffIdx;	// currently buffering this many output bytes
size_t m_AllocOutBuff;	// m_pOutBuffer allocated to hold this many input bytes
uint8_t *m_pOutBuffer;	// allocated for buffering output

char szAssembRef[120];
char szSeqRef[120];
uint32_t SeqRefLen;
char szPrefix[10];
strcpy(szPrefix,pszPrefix);
PrefixLen = (int)strlen(szPrefix);
szPrefix[PrefixLen++] = '|';
szPrefix[PrefixLen] = '\0';

#ifdef _WIN32
m_hInFile = open(pszInFile, O_READSEQ );		// file access is normally sequential..
#else
m_hInFile = open64(pszInFile, O_READSEQ );		// file access is normally sequential..
#endif
if(m_hInFile == -1)							// check if file open succeeded
	{
	gDiagnostics.DiagOut (eDLFatal, gszProcName, "Unable to open input file '%s' : %s",pszInFile,strerror(errno));
	return(eBSFerrOpnFile);
	}

#ifdef _WIN32
m_hOutFile = open(pszOutFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE));
#else
if((m_hOutFile = open64(pszOutFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE)) != -1)
	if(ftruncate(m_hOutFile,0)!=0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create/truncate %s - %s",pszOutFile,strerror(errno));
		close(m_hInFile);
		return(eBSFerrCreateFile);
		}
#endif
if(m_hOutFile < 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create/truncate %s - %s",pszOutFile,strerror(errno));
	close(m_hInFile);
	return(eBSFerrCreateFile);
	}

if((m_pInBuffer = new uint8_t[cAllocPGBuffInSize]) == NULL)
	{
	gDiagnostics.DiagOut (eDLFatal, gszProcName, "Unable to allocate memory for input file buffering");
	close(m_hInFile);
	close(m_hOutFile);
	return(eBSFerrMem);
	}
m_AllocInBuff = cAllocPGBuffInSize;

if((m_pOutBuffer = new uint8_t[cAllocPGBuffOutSize]) == NULL)
	{
	gDiagnostics.DiagOut (eDLFatal, gszProcName, "Unable to allocate memory for output file buffering");
	delete []m_pInBuffer;
	close(m_hInFile);
	close(m_hOutFile);
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
while((NumRead = (int)read(m_hInFile, m_pInBuffer, (int)m_AllocInBuff)) > 0)
	{
	m_InBuffIdx = NumRead;
	pInChr = m_pInBuffer;
	for(Idx = 0; Idx < m_InBuffIdx; Idx++)
		{
		InChr = *pInChr++;
		*pOutChr++ = InChr;
		m_OutBuffIdx++;
		SAMreclen++;
		if(InChr == '\0' || InChr == '\n')	// have a complete SAM record?
			{
			pChr = &m_pOutBuffer[RecordStartIdx];
			if(*pChr == '@')		// is SAM record a header?
				{
				if(pChr[1] == 'S' && pChr[2] == 'Q')	// only retaining sequence names with specified prefix
					{
					if(sscanf((char *)&pChr[3],"\tAS:%s\tSN:%s\tLN:%u",szAssembRef,szSeqRef,&SeqRefLen)==3)
						{
						if(strnicmp(szSeqRef, szPrefix, PrefixLen))
							{
							SAMreclen = 0;
							m_OutBuffIdx = RecordStartIdx;
							pOutChr = &m_pOutBuffer[RecordStartIdx];
							continue;
							}
						// sequence name has specified prefix, strip prefix off and accept this header record
						sprintf((char *)&pChr[3],"\tAS:%s\tSN:%s\t:LN:%u\n",szAssembRef,&szSeqRef[PrefixLen],SeqRefLen);
						m_OutBuffIdx -= (PrefixLen-1);
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
				do {
					if(*pChr++ == '\t')
						NumTabs++;
					}
				while(--SAMreclen && NumTabs < 2);

				if(NumTabs != 2 || strnicmp((char *)pChr, szPrefix, PrefixLen))
					{
					SAMreclen = 0;
					m_OutBuffIdx = RecordStartIdx;
					pOutChr = &m_pOutBuffer[RecordStartIdx];
					continue;
					}

				// strip prefix off
				do {
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
			CUtility::SafeWrite(m_hOutFile, m_pOutBuffer, m_OutBuffIdx);	
			m_OutBuffIdx = 0;
			RecordStartIdx = 0;
			pOutChr = m_pOutBuffer;
			SAMreclen = 0;
			}
		}
	if(InChr == '\0')
		{
		*pOutChr++ = '\0';
		RecordStartIdx++;
		break;
		}
	}

if(RecordStartIdx)
	CUtility::SafeWrite(m_hOutFile, m_pOutBuffer, RecordStartIdx);	

// close input file
close(m_hInFile);
// commit and close output file
#ifdef _WIN32
_commit(m_hOutFile);
#else
fsync(m_hOutFile);
#endif
close(m_hOutFile);
delete []m_pInBuffer;
delete []m_pOutBuffer;
return(eBSFSuccess);
}


int Process(eModePG PMode,			// processing mode
			char* pszPrefix,		// descriptor prefix
			char* pszInFile,		// input fasta or SAM file
			char* pszOutFile)		// output to this file
{
int Rslt = eBSFerrParams;

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
	}


return(Rslt);
}