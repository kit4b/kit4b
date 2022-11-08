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
#ifdef _WIN32
#include "./commhdrs.h"
#else
#include "./commhdrs.h"
#endif

CFasta::CFasta(void)
{
m_hFile = -1;
m_gzFile = NULL;
m_pszFastqSeq = NULL;
m_pszFastqSeqQ = NULL;
memset(m_FastaBlocks, 0, sizeof(m_FastaBlocks));
m_pCurFastaBlock = NULL;
Cleanup();
}

CFasta::~CFasta(void)
{
if (m_hFile >= 0)
	close(m_hFile);

if (m_gzFile != NULL)
	gzclose(m_gzFile);

if (m_pszFastqSeq != NULL)
	delete[]m_pszFastqSeq;

if (m_pszFastqSeqQ != NULL)
	delete[]m_pszFastqSeqQ;

for (int32_t Idx = 0; Idx < cNumFastaBlocks; Idx++)
	{
	if (m_FastaBlocks[Idx].pBlock != NULL)
		delete[]m_FastaBlocks[Idx].pBlock;
	}
}

void
CFasta::Cleanup(void)
{
if(m_hFile >= 0)
	{
	if(!m_bRead)
		{
#ifdef _WIN32
		_commit(m_hFile);
#else
		fsync(m_hFile);
#endif
		}
	close(m_hFile);
	m_hFile = -1;
	}
if(m_gzFile != NULL)
	{
	gzclose(m_gzFile);
	m_gzFile = NULL;
	}
if (m_pszFastqSeq != NULL)
	{
	delete []m_pszFastqSeq;
	m_pszFastqSeq = NULL;
	}
if (m_pszFastqSeqQ != NULL)
	{
	delete []m_pszFastqSeqQ;
	m_pszFastqSeqQ = NULL;
	}
for (int32_t Idx = 0; Idx < cNumFastaBlocks; Idx++)
	{
	if (m_FastaBlocks[Idx].pBlock != NULL)
		delete []m_FastaBlocks[Idx].pBlock;
	}
memset(m_FastaBlocks, 0, sizeof(m_FastaBlocks));
m_pCurFastaBlock = NULL;
m_FileDescrOfs = 0;
m_FileReadDescrOfs = 0;
m_szDescriptor[0] = '\0';
m_bDescrAvail = false;
m_bForceRetSeq = false;
m_DescriptorLen = 0;
m_szFile[0] = '\0';
m_CurLineLen = 0;
m_bRead = true;
m_bIsGZ = false;
m_bIsFastQ = false;
m_FastqSeqIdx = 0; 
m_FastqSeqLen = 0;
m_CurFastQParseLine = 0;
m_NumMissingSequences = 0;
}

uint64_t
CFasta::InitialFileSize(void)				// file size when initially opened for reading
{
return(m_StatFileSize);
}

int32_t
CFasta::Open(char *pszFile,						// fasta or fastq file path+name to open
			 bool Read,							// TRUE if opening for read, FALSE for write
			 uint32_t BufferSize)			// use this size buffer for staging
{
int32_t Rslt;
if(pszFile == NULL || *pszFile == '\0')
	return(eBSFerrParams);

// force buffer size  to be within acceptable range
if(BufferSize < cMinStageBuffSize)
	BufferSize = cMinStageBuffSize;
else
	if(BufferSize > cMaxStageBuffSize)
		BufferSize = cMaxStageBuffSize;

Cleanup();	

if (m_pszFastqSeq == NULL)
	{
	if ((m_pszFastqSeq = new char[cMaxFastQSeqLen+10]) == NULL)
		return(eBSFerrMem);
	}
if (m_pszFastqSeqQ == NULL)
	{
	if((m_pszFastqSeqQ = new char[cMaxFastQSeqLen+10]) == NULL)
		return(eBSFerrMem);
	}

#ifdef _WIN32
struct _stat64 st;
if(!_stat64(pszFile,&st))
#else
struct stat64 st;
if(!stat64(pszFile,&st))
#endif
	m_StatFileSize = (int64_t)st.st_size;
else
	m_StatFileSize = 0;

if(Read)		// read
	{
	strcpy(m_szFile,pszFile);

	// if file has extension of ".gz' then assume that this file has been compressed and needs processing with gzopen/gzread/gzclose
	int32_t NameLen = (int32_t)strlen(pszFile);
	if(NameLen >= 4 && !stricmp(".gz",&pszFile[NameLen-3]))
		{
		if((m_gzFile = gzopen(pszFile,"r"))==NULL)
			{
			AddErrMsg("CFasta::Open","Unable to open %s as a gzip'd file - %s",pszFile,strerror(errno));
			Rslt = eBSFerrOpnFile;
			Cleanup();
			return(Rslt);
			}

		if(gzbuffer(m_gzFile,cgzAllocInBuffer)!=0)
			{
			AddErrMsg("CFasta::Open","Unable to set gzbuffer size to %d",cgzAllocInBuffer);
			Rslt = eBSFerrMem;
			Cleanup();
			return(Rslt);
			}

		if(!gzdirect(m_gzFile))			// false if file has been compressed
			m_bIsGZ = true;
		else
			m_bIsGZ = false;
		}
	else
		{
		m_bIsGZ = false;
#ifdef _WIN32
		m_hFile = open(pszFile, O_READSEQ );		// file access is normally sequential..
#else
		m_hFile = open64(pszFile, O_READSEQ );		// file access is normally sequential..
#endif
		if(m_hFile == -1)							// check if file open succeeded
			{
			AddErrMsg("CFasta::Open","Unable to open %s - %s",pszFile,strerror(errno));
			Rslt =eBSFerrOpnFile;
			Cleanup();
			return(eBSFerrOpnFile);
			}
		}

	// Check if file looks like it contains sequence data, either as a fasta or fastq file
	if(Read && (Rslt = CheckIsFasta()) != eBSFSuccess)
		{
		AddErrMsg("CFasta::Open","File %s exists but not a fasta or fastq file",pszFile);
		Cleanup();
		return(Rslt);
		}

	}
else			// write
	{
#ifdef _WIN32
	m_hFile = open(pszFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE));
#else
	if((m_hFile = open64(pszFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE)) != -1)
		if(ftruncate(m_hFile,0)!=0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",pszFile,strerror(errno));
			return(eBSFerrCreateFile);
			}
#endif
	if(m_hFile < 0)
		{
		AddErrMsg("CFasta::Open","Unable to open %s - %s",pszFile,strerror(errno));
		Rslt = eBSFerrCreateFile ;
		Cleanup();
		return(Rslt);
		}
	m_StatFileSize = 0;
	}

if((uint64_t)BufferSize > (m_StatFileSize+100))			// a little additional never hurts!
	BufferSize = (uint32_t)(m_StatFileSize + 100);

memset(m_FastaBlocks, 0, sizeof(m_FastaBlocks));

m_pCurFastaBlock = &m_FastaBlocks[0];
m_pCurFastaBlock->pBlock = new uint8_t[BufferSize];
if (m_pCurFastaBlock->pBlock == NULL)
	{
	AddErrMsg("CFasta::Open", "Memory allocation of %d bytes for %s- %s", BufferSize, pszFile, strerror(errno));
	Cleanup();			// closes opened file..
	return(eBSFerrMem);
	}
m_pCurFastaBlock->AllocSize = BufferSize;
if (m_bRead && cNumFastaBlocks > 1)
	{
	for (int32_t Idx = 1; Idx < cNumFastaBlocks; Idx++)
		{
		m_FastaBlocks[Idx].pBlock = new uint8_t[BufferSize];
		if (m_FastaBlocks[Idx].pBlock == NULL)
			{
			AddErrMsg("CFasta::Open", "Memory allocation of %d bytes for %s- %s", BufferSize, pszFile, strerror(errno));
			Cleanup();			// closes opened file..
			return(eBSFerrMem);
			}
		m_FastaBlocks[Idx].AllocSize = BufferSize;
		}
	}

if((Rslt=Reset())!=eBSFSuccess)
	return(Rslt);
m_bRead = Read;
return(eBSFSuccess);
}


// FastaEstimateSizes
// Generally suitable only for NGS read dataset size estimates
const int32_t cEstChrsBuff = 0x03fffff;					// base estimates on the first 4M chars of file

uint32_t								// returns estimated number of sequences in fata/fastq file
CFasta::FastaEstSizes(char *pszFile,				// fasta or fastq file path+name to estimate sizes
			  int64_t *pFileSize,						// file is this size on disk
			  int32_t *pEstMaxDescrLen,				// with estimated maximum descriptor length
			  int32_t *pEstMeanDescrLen,				// estimated mean descriptor length
			  int32_t *pEstMaxSeqLen,					// and estimated maximum sequence length
			  int32_t *pEstMeanSeqLen,				// estimated mean sequence length
			  int32_t *pEstScoreSchema)				// guestimated scoring schema - 0: no scoring, 1: Solexa, 2: Illumina 1.3+, 3: Illumina 1.5+, 4: Illumina 1.8+ or could be Sanger 
{
int64_t FileSize;
uint32_t NumSeqs;
char MaxScoreChr;		// max score char in any of the reads, used for reasonable guestimate of scoring schema utilised
char MinScoreChr;		// min score char in any of the reads, used for reasonable guestimate of scoring schema utilised

gzFile gzdFile;
int32_t hFile;
bool bIsGZ;
long FileOfs;

int32_t BuffSize;
int32_t NumInBuff;
uint8_t *pBuff;
char Chr;
char *pChr;
int32_t NumChrsParsed;
int32_t NumChrsParsedThisLine;

bool bInDescrLine;
bool bInSeqLine;
int32_t NumLines2Slough;

int32_t CurseqLen;
int32_t MaxSeqLen;
int32_t TotSeqLen;
int32_t MaxDescrLen;
int32_t TotDescrLen;
uint32_t EstNumSeqs;

MaxScoreChr = 0;
MinScoreChr = 0;

if(pEstMaxDescrLen != NULL)
	*pEstMaxDescrLen = 0;
if(pEstMaxSeqLen != NULL)
	*pEstMaxSeqLen = 0;
if(pEstMeanDescrLen != NULL)
	*pEstMeanDescrLen = 0;
if(pEstMeanSeqLen != NULL)
	*pEstMeanSeqLen = 0;
if(pFileSize != NULL)
	*pFileSize = 0;

BuffSize = cEstChrsBuff;

#ifdef _WIN32
struct _stat64 st;
if(!_stat64(pszFile,&st))
#else
struct stat64 st;
if(!stat64(pszFile,&st))
#endif
	FileSize = (int64_t)st.st_size;
else
	FileSize = 0;

if(pFileSize != NULL)
	*pFileSize = FileSize;

if(FileSize == 0)		// 0 if file not readable or if 0 length
	{
	gDiagnostics.DiagOut(eDLWarn,gszProcName,"FastaEstSizes: Unable to estimate sizes for file '%s', does file exist and not 0 length, or is it readable",pszFile);
	return(0);
	}

if(FileSize < 10)		// arbitary minimum Fasta file size, allows for short descriptor + 1 base sequence...
	{
	gDiagnostics.DiagOut(eDLWarn,gszProcName,"FastaEstSizes: Unable to estimate sizes for file '%s', file exists but is only %d long",pszFile,FileSize);
	return(0);
	}

if((pBuff = new uint8_t [BuffSize+10])==NULL)			// allows for trailing '\0's
	{
	gDiagnostics.DiagOut(eDLWarn,gszProcName,"FastaEstSizes: Unable to allocate %d memory for file '%s'",BuffSize,pszFile);
	return(0);
	}


int32_t NameLen = (int32_t)strlen(pszFile);
if(NameLen >= 4 && !stricmp(".gz",&pszFile[NameLen-3]))
	{
	if((gzdFile = gzopen(pszFile,"r"))==NULL)
		{
		delete []pBuff;
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"FastaEstSizes: Unable to open gzip file '%s'",pszFile);
		return(0);
		}

	if(!gzdirect(gzdFile))			// false if file has been compressed
		bIsGZ = true;
	else
		bIsGZ = false;
	gzseek(gzdFile,0,SEEK_SET);
	NumInBuff = gzread(gzdFile,pBuff,BuffSize);
	FileOfs = gzoffset(gzdFile);
	gzclose(gzdFile);
	}
else
	{
	if((hFile = open(pszFile,O_READSEQ))== -1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"FastaEstSizes: Unable to open %s - %s",pszFile,strerror(errno));
		delete []pBuff;
		return(0);
		}
	NumInBuff = read(hFile,pBuff,BuffSize);
	close(hFile);
	bIsGZ = false;
	}

if(NumInBuff < 10)			// what kind of fasta file would be this small :-)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"FastaEstSizes: Unable to estimate sizes for file '%s', file exists but is only %d long",pszFile,NumInBuff);
	delete []pBuff;
	return(0);
	}
pBuff[NumInBuff] = 0;

pChr = (char *)pBuff;
NumChrsParsed = 0;
NumChrsParsedThisLine = 0;
MaxSeqLen = 0;
MaxDescrLen = 0;
NumSeqs = 0;
TotDescrLen = 0;
TotSeqLen = 0;
CurseqLen = 0;
bInDescrLine = false;
bInSeqLine = false;					
NumLines2Slough = 0;
while((Chr = *pChr++) != '\0')
	{
	NumChrsParsed += 1;

	if(Chr == '\r' || Chr == '\n')				// end of current line?
		{
		while(*pChr == '\r' || *pChr == '\n')	// slough repeated CR/LF's
			{
			pChr++;
			NumChrsParsed += 1;
			}
		if(NumLines2Slough)
			{
			NumLines2Slough -= 1;
			NumChrsParsedThisLine = 0;
			continue;
			}

		if(bInDescrLine)
			{
			TotDescrLen += NumChrsParsedThisLine;
			if(NumChrsParsedThisLine > MaxDescrLen)
				MaxDescrLen = NumChrsParsedThisLine;
			CurseqLen = 0;
			bInSeqLine = true;
			bInDescrLine = false;
			}
		else
			{
			if(bInSeqLine)
				{
				TotSeqLen += NumChrsParsedThisLine;
				CurseqLen += NumChrsParsedThisLine;
				if(CurseqLen > MaxSeqLen)
					MaxSeqLen = CurseqLen;
				Chr = *pChr;
				if(Chr == '>' || Chr == '@' || Chr == '+' || Chr == '\0')
					{
					bInSeqLine = false;
					NumSeqs += 1;
					}
				}
			}
		NumChrsParsedThisLine = 0;
		continue;
		}

	if (NumLines2Slough)
		{
		if (NumLines2Slough == 1)	// if sloughing the quality scores then determine min/max score chars so a guestimate can be made of the scoring schema utilised
			{
			if (Chr > MaxScoreChr)
				MaxScoreChr = Chr;
			if (MinScoreChr == 0 || MinScoreChr > Chr)
				MinScoreChr = Chr;
			}
		continue;
		}

	NumChrsParsedThisLine += 1;
	if(bInDescrLine || bInSeqLine)
		continue;

	if(NumChrsParsedThisLine==1)
		{
		if(Chr == '>' || Chr == '@' || Chr == '+')			// starting descriptor line for fasta or fastq; or separator descriptor line between sequence and phred scores in fastq
			{
			if (Chr != '+')
				bInDescrLine = true;
			else
				NumLines2Slough = 2;
			continue;
			}
		}
	}

delete []pBuff;

if(!bIsGZ)
	{
	if((int64_t)NumChrsParsed == FileSize)
		EstNumSeqs = NumSeqs;
	else
		EstNumSeqs = (uint32_t)(((NumSeqs + 1) * FileSize)/(uint64_t)NumChrsParsed);	// better to slightly over estimate than underestimate..
	}
else
	{
	double ComprRatio = (double)NumChrsParsed/FileOfs;
	int64_t UncomprFileSize = (int64_t)(FileSize * ComprRatio);
	if((int64_t)NumChrsParsed >= UncomprFileSize)
		EstNumSeqs = NumSeqs;
	else
		EstNumSeqs = (uint32_t)(((NumSeqs + 1) * UncomprFileSize)/(uint64_t)NumChrsParsed);	// better to slightly over estimate than underestimate..
	}


if(pEstMaxDescrLen != NULL)
	*pEstMaxDescrLen = MaxDescrLen;
if(pEstMaxSeqLen != NULL)
	*pEstMaxSeqLen = MaxSeqLen;

if(pEstMeanDescrLen != NULL)
	*pEstMeanDescrLen = TotDescrLen/NumSeqs;
if(pEstMeanSeqLen != NULL)
	*pEstMeanSeqLen = TotSeqLen/NumSeqs;

if (pEstScoreSchema != NULL)
	{
	if (MaxScoreChr >= 'K')		// solexa, illumina 1.3+ and illuminia 1.5+ all have scores >= 'K'
		{
		if (MinScoreChr <= '?')		// solexa starts from ':' on up to 'h'
			*pEstScoreSchema = 1;
		else
			if (MinScoreChr <= 'A')		// illumina 1.3+ starts from '@' on up to 'h'
				*pEstScoreSchema = 2;
			else
				*pEstScoreSchema = 3;	// illumina 1.5+ starts from 'B' on up to 'h'
		}
	else
		if (MaxScoreChr == 0)		// no associated scoring if no score lines processed
			*pEstScoreSchema = 0;
		else
			*pEstScoreSchema = 4;	// illumina 1.8+ starts from '#' on up to 'J'
	}
return(EstNumSeqs);
}

// CheckIsFasta
// Crude check on presumed fasta or fastq file contents
// Reads 1st 1Mbp of contents and checks if contents are consistent with fasta or fastq file formats
// Also handles Fastq with SOLiD 2 base representation
const int32_t cChkFastaSize = 0x100000;

int32_t
CFasta::CheckIsFasta(void)
{
char *pszBuff;
int32_t BuffCnt;
char *pChr;
char Base;
int32_t BuffIdx;
int32_t CurLineNum;		// line currently being processed
int32_t ChrPsn;			// chr psn 1..n in line being processed
bool bIsFastQSOLiD;	// some fastq files (from NCBA SRA SRP000191) have SOLiD 2 base representation sequences
int32_t NumSeqChrs;		// cur number of sequence line nucleotides processed
int32_t NumAmbiguousBases; // number of ambiguous bases in current fasta sequence
int32_t NumQualChrs;	// cur number of quality line scores processed
int32_t FileState;		// 0 unknown, 1 in fasta, 2 in fasta descr, 3 in fasta seq, 4 in fastq, 5 in fastq seq, 6 waiting for '+', 7 in '+', 8 quality values
bool bSkipEOL;		// true if remainder of current line not to be parsed
bool bIscsfasta;	// true if file has been determined to be SOLiD csfasta format
bool bIsfasta;		// true if file has been determined to be basespace

if(m_gzFile == NULL && m_hFile == -1)
	return(eBSFerrFileClosed);

if((pszBuff = new char [cChkFastaSize]) == NULL)
	return(eBSFerrMem);

if(m_gzFile != NULL)
	{
	gzseek(m_gzFile,0,SEEK_SET);
	BuffCnt = gzread(m_gzFile,pszBuff,cChkFastaSize-1);
	gzseek(m_gzFile,0,SEEK_SET);
	}
else
	{
	if(_lseeki64(m_hFile,0,SEEK_SET)!=0)		
		return(eBSFerrFileAccess);
	BuffCnt = read(m_hFile,pszBuff,cChkFastaSize-1);
	if(_lseeki64(m_hFile,0,SEEK_SET)!=0)	
		{
		delete pszBuff;	
		return(eBSFerrFileAccess);
		}
	}

if(BuffCnt < 4)		// not interested in fasta or fastq sequences which are too short!
	{
	AddErrMsg("CFasta::CheckIsFasta","Errors whilst reading file - '%s' file too short, less than 4 chrs", m_szFile);
	delete pszBuff;	
	return(eBSFerrNotFasta);
	}

// quick check in case user has specified a bioseq file
if((pszBuff[0] == 'b' && pszBuff[1] == 'i' && pszBuff[2] == 'o') ||
   (pszBuff[0] == 's' && pszBuff[1] == 'f' && pszBuff[2] == 'x'))
	{
	AddErrMsg("CFasta::CheckIsFasta","Not a fasta file, could be a bioseq/biosfx file - '%s'", m_szFile);
	delete pszBuff;	
	return(eBSFerrNotFasta);
	}

FileState = 0;
pChr = pszBuff;
CurLineNum = 1;
bIsFastQSOLiD = false;
NumSeqChrs = 0;
NumQualChrs = 0;
NumAmbiguousBases = 0;
ChrPsn = 0;
bSkipEOL = false;
bIscsfasta = false;
bIsfasta = false; 
for(BuffIdx = 0; BuffIdx < BuffCnt; BuffIdx++,pChr++)
	{
	if(bSkipEOL == true)
		{
		if(*pChr != '\n')
			continue;
		bSkipEOL = false;
		}

	if(*pChr == '\r')			// simply slough CR's - NL's are the line terminators of interest
		{
		ChrPsn = 0;
		continue;
		}	

	if(*pChr == '\n')
		ChrPsn = 0;
	else
		ChrPsn += 1;
	if(!isspace(*pChr) && (*pChr < 0x20 || (uint8_t)*pChr > 0x7f))
		{
		AddErrMsg("CFasta::CheckIsFasta","Errors whilst reading file near line %d#%d - '%s' - illegal chr 0x%x", CurLineNum,ChrPsn, m_szFile,*pChr);
		delete pszBuff;	
		return(eBSFerrNotFasta);
		}
	
	switch(FileState) {
		case 0:		// file type not yet determined
			if(isspace(*pChr))			// be very tolerant of whitespace on reads...
				{
				if(*pChr == '\n')
					CurLineNum += 1;
				continue;
				}
			switch(Base = tolower(*pChr)) {
				case '>':		// indicates in fasta, but could be csfasta if no header present
					FileState = 2; // and reading descriptor line
					continue;
				case '@':		// indicates in fastq
					bIsfasta = true;
					FileState = 4; // and reading sequence identifier line
					continue;
				case 'a': case 'c': case 'g': case 't': case 'u': case 'n': case '-':
					// if an ambiguity sequence char then providing not too many, accept as most likely fasta
					case 'k':		// G || T 
					case 'm':		// A || C
					case 'r':		// A || G
					case 'y':		// C || T
					case 's':		// C || G
					case 'w':		// A || T
					case 'b':		// C || G || T
					case 'v':		// A || C || G
					case 'h':       // A || C || T
					case 'd':		// A || G || T
					case 'x':		// A || C || G || T
					case '.':		// not (A || C || G || T)
					bIsfasta = true;
					FileState = 3;	// assume fasta and reading sequence
					NumAmbiguousBases = 0;
					continue;
				case '#':		// SOLiD files csfasta files may contain a header with header lines starting with '#"
					bSkipEOL = true; // simply slough this line
					bIscsfasta = true;
					continue;

				default:		
					break;
				}
			AddErrMsg("CFasta::CheckIsFasta","Errors whilst reading file near line %d#%d - '%s' - unexpected chr '%c'",CurLineNum,ChrPsn, m_szFile,*pChr);
			delete pszBuff;	
			return(eBSFerrNotFasta);

		case 1:	// presumed fasta and expecting either descriptor or sequence
			if(isspace(*pChr))
				{
				if(*pChr == '\n')
					CurLineNum += 1;
				continue;
				}
			switch(tolower(*pChr)) {
				case '>':
					FileState = 2;	// descriptor
					continue;
				case 'a': case 'c': case 'g': case 't':
				case 'u': case 'n': case '-':
				case 'k': case 'm':	case 'r': case 'y':	case 's': case 'w':	case 'b': case 'v':	case 'h': case 'd':	case 'x': case '.':
					FileState = 3;	// sequence
					NumAmbiguousBases = 0;
					continue;
				default:			// any other char is an error
					break;
				}
			AddErrMsg("CFasta::CheckIsFasta","Errors whilst reading file near line %d#%d - '%s' - illegal chr '%c'",CurLineNum,ChrPsn, m_szFile,*pChr);
			delete pszBuff;	
			return(eBSFerrNotFasta);

		case 2:	// presumed fasta file and reading the descriptor line
			if(*pChr != '\n')		// anything is accepted within descriptors
				continue;
			CurLineNum += 1;
			FileState = 1;	// back to either descriptor or sequence
			continue;

		case 3: // presumed fasta or csfasta file and reading the sequence line
			if(*pChr == '\n')
				{
				CurLineNum += 1;
				NumSeqChrs = 0;
				FileState = 1;	// back to either descriptor or sequence
				continue;
				}
			if(isspace(*pChr))	// pays to be tolerant...
				continue;
			if(bIscsfasta)
				{
				switch(*pChr) {
					case '0': case '1': case '2': case '3': case '4': case '.':  // '4' and '.' represents an indeterminate call
						NumSeqChrs += 1;
						continue;
					default:
						break;
					}
				}
			else
				switch(tolower(*pChr)) {
				case 'a': case 'c': case 'g': case 't':	case 'u': case 'n':
						if(NumAmbiguousBases)
							NumAmbiguousBases--;
						bIsfasta = true;
						continue;

					case '-': case 'k': case 'm':	case 'r': case 'y':	case 's': case 'w':	case 'b': case 'v':	case 'h': case 'd':	case 'x': case '.':
						NumSeqChrs += 1;
						if(NumSeqChrs > 50 && NumAmbiguousBases++ > (NumSeqChrs / 2))	// this limit is purely arbitrary
							{
							AddErrMsg("CFasta::CheckIsFasta","Too many ambiguous bases in assumed fasta sequence near line %d#%d - '%s' - non-sequence chr '%c'",CurLineNum,ChrPsn, m_szFile,*pChr);
							break;
							}
						bIsfasta = true;
						continue;

					case '0': case '1': case '2': case '3': case '4':   // could be colorspace with no header???
						if(bIsfasta)	
							break;
						NumSeqChrs += 1;
						bIscsfasta = true;
						continue;

					default:
						break;
					}
			AddErrMsg("CFasta::CheckIsFasta","Errors whilst reading file near line %d#%d - '%s' - non-sequence chr '%c'",CurLineNum,ChrPsn, m_szFile,*pChr);
			delete pszBuff;	
			return(eBSFerrNotFasta);
		
		case 4:	// presumed fastq file and seq identifier line
			if(*pChr != '\n')		// anything is accepted within sequence identifiers
				continue;
			CurLineNum += 1;
			FileState = 5;	// next is the sequence
			NumSeqChrs = 0;
			bIsFastQSOLiD = false;
			continue;

		case 5: // presumed fastq file and reading the sequence line
			if(*pChr == '\n')
				{
				FileState = 6;	// onto the '+' line
				CurLineNum += 1;
				continue;
				}
			if(isspace(*pChr))	// tolerant of spaces
				continue;
			
			switch(Base = tolower(*pChr)) {
				case 'a': case 'c': case 'g': case 't': case 'u': case 'n':
				case 'k': case 'm':	case 'r': case 'y':	case 's': case 'w':	case 'b': case 'v':	case 'h': case 'd':	case 'x': case '.':
					if(bIsFastQSOLiD && Base != '.')
						{
						AddErrMsg("CFasta::CheckIsFasta","Errors whilst reading file near line %d#%d - '%s' - expected SOLiD base representation", CurLineNum,ChrPsn,m_szFile);
						delete pszBuff;	
						return(eBSFerrNotFasta);
						}
					NumSeqChrs+=1;
					continue;
					
				case '0': case '1': case '2': case '3': case '4':		// SOLiD representation...
					if(NumSeqChrs == 1)
						bIsFastQSOLiD = true;
					if(bIsFastQSOLiD)
						{
						NumSeqChrs += 1;
						continue;
						}
					AddErrMsg("CFasta::CheckIsFasta","Errors whilst reading file near line %d#%d - '%s' - unexpected SOLiD base representation", CurLineNum,ChrPsn,m_szFile);
					delete pszBuff;	
					return(eBSFerrNotFasta);
						
				default:
					break;
				}
			AddErrMsg("CFasta::CheckIsFasta","Errors whilst reading file near line %d#%d - '%s' - non-sequence chr '%d'", CurLineNum,ChrPsn,m_szFile,*pChr);
			delete []pszBuff;	
			return(eBSFerrNotFasta);

		case 6:	// presumed fastq file and reading the duplicate seq ident line
			if(isspace(*pChr))		// need to see a '+'
				{
				if(*pChr == '\n')
					CurLineNum += 1;
				continue;
				}
			if(*pChr == '+')
				{
				FileState = 7;
				continue;
				}
			AddErrMsg("CFasta::CheckIsFasta","Errors whilst reading file near line %d#%d - '%s' - expected '+' not chr '%c'", CurLineNum,ChrPsn,m_szFile,*pChr);
			delete []pszBuff;	
			return(eBSFerrNotFasta);

		case 7:					// processing fastq duplicate seq identifier line
			if(*pChr != '\n')		// anything is accepted within seq identifier line
				continue;
			CurLineNum += 1;
			FileState = 8;	// onto quality
			NumQualChrs = 0;
			continue;

		case 8: // presumed fastq file and reading the quality line
			if(*pChr != '\n')
				{
				if(isspace(*pChr))		// be very tolerant of whitespace
					continue;
				NumQualChrs += 1;
				continue;
				}
			if(NumQualChrs != NumSeqChrs)
				{
				AddErrMsg("CFasta::CheckIsFasta","Errors whilst reading file near line %d - '%s' - sequence length (%d) not same as quality length (%d)", CurLineNum,m_szFile,NumSeqChrs,NumQualChrs);
				delete []pszBuff;	
				return(eBSFerrNotFasta);
				}
			CurLineNum += 1;
			NumSeqChrs = 0;
			NumQualChrs = 0;
			FileState = 9;			// expecting new fastq block
			continue;
		case 9:						// assumed fastq file and expecting a new block to start
			if(*pChr == '\n')		// be tolerant of newlines
				{
				CurLineNum += 1;
				continue;
				}
			if(*pChr == '@')
				{
				FileState = 4;
				continue;
				}
			AddErrMsg("CFasta::CheckIsFasta","Errors whilst reading file near line %d#%d - '%s' - expected '@' not chr '%c'", CurLineNum,ChrPsn,m_szFile,*pChr);
			delete []pszBuff;	
			return(eBSFerrNotFasta);
		}
	}
delete []pszBuff;	
m_bIsFastQ = FileState >= 4 ? true : false;
m_bIscsfasta = bIscsfasta;
return(eBSFSuccess);	
}

// IsFastQ
// Returns true if currently opened file is in fastq format
bool
CFasta::IsFastq(void)
{
return(m_bIsFastQ);
}

// IsSOLiD
// Returns true if currently opened file is in SOLiD format
bool
CFasta::IsSOLiD(void)
{
return(m_bIscsfasta);
}

// Reset
// Resets context to that immediately following an Open() with file seek to FileOfs
// NOTE: gzip library can't handle file offsets which are greater than 2^31 - 1
int32_t
CFasta::Reset(int64_t FileOfs)
{
if(m_hFile == -1 && m_gzFile == NULL)
	return(eBSFerrClosed);		
int64_t SeekPsn;
if(m_hFile != -1)
	SeekPsn = _lseeki64(m_hFile,FileOfs,SEEK_SET);
else
	SeekPsn = gzseek(m_gzFile,(long)FileOfs,SEEK_SET);
if(SeekPsn != FileOfs)
	{
	AddErrMsg("CFasta::Reset","Seek failed to offset %d on %s - %s",FileOfs,m_szFile,strerror(errno));
	Cleanup();
	return(eBSFerrFileAccess);
	}

m_pCurFastaBlock = &m_FastaBlocks[0];
m_pCurFastaBlock->FileOfs = 0;
m_pCurFastaBlock->BuffCnt = 0;
m_pCurFastaBlock->BuffIdx = 0;

m_FileReadDescrOfs = 0;
m_FileDescrOfs = 0;
m_CurLineLen = 0;
m_szDescriptor[0] = '\0';
m_bDescrAvail = false;
m_bForceRetSeq = false;
m_DescriptorLen = 0;
m_CurFastQParseLine = 0;
return(eBSFSuccess);
}


// ReadSequence
// Returns upto Max2Read bases from fasta or fastq input file
// Sloughs all whitespace and non-alpha characters within sequences
// If a csfasta file then sloughs lines starting with '#'
//
// Will return -
//    eBSFSuccess		fasta file has been completely processed
//	  eBSFFastaDescr    a descriptor line has been just processed, call GetDescriptor()
//	  n					number of sequence symbols returned
//    < 0				error
// 
//  
uint8_t CFasta::m_SOLiDmap[5][5] = {
	{'a','c','g','t','n'},	 
	{'c','a','t','g','n'},		
	{'g','t','a','c','n'},
	{'t','g','c','a','n'},
	{'n','n','n','n','n'}};		


int32_t						// returns actual number read (eBSFSuccess == EOF,eBSFFastaDescr == End of current sequence, descriptor line now available)
CFasta::ReadSequence(void *pRetSeq,		// where to return sequence, can be NULL if only interested in the sequence length
					 int32_t Max2Ret,		// max to return, ignored if pRetSeq is NULL
					 bool bSeqBase,		// if false then return ascii, if true then return as etSeqBase
					 bool RptMskUpperCase)	// default is false, UCSC softmasked use lowercase when repeat masked
{
bool bInDescriptor;		// true whilst processing descriptor or fastq sequence identifier characters
bool bMoreToDo;
char Chr;
int64_t FileOfs;			// will contain file offset corresponding to last buffer read from file
int32_t SeqLen = 0;
char *pAscii = (char *)pRetSeq;
int32_t Rslt;
bool bSloughEOL;	// if true then skip to end of current line
int32_t PrevSOLiDbase;

if (m_gzFile == NULL && m_hFile == -1 || m_pCurFastaBlock == NULL || m_pCurFastaBlock->pBlock == NULL)
	return(eBSFerrClosed);
if(!m_bRead)
	return(eBSFerrRead);
if(pRetSeq != NULL && !Max2Ret)
	return(eBSFerrParams);		

if(m_bIsFastQ)
	{
	if(m_FastqSeqIdx != m_FastqSeqLen)
		{
		if(pRetSeq == NULL)
			{
			if(Max2Ret == 0)	// user only interested in the actual sequence length?
				SeqLen = m_FastqSeqLen - m_FastqSeqIdx;
			else
				SeqLen = min(Max2Ret,m_FastqSeqLen - m_FastqSeqIdx); 
			m_FastqSeqIdx += SeqLen;
			return(SeqLen);
			}

		SeqLen = min(Max2Ret,m_FastqSeqLen - m_FastqSeqIdx); 
		if(pRetSeq != NULL && bSeqBase)
			Ascii2Sense((char *)&m_pszFastqSeq[m_FastqSeqIdx],SeqLen,(etSeqBase *)pRetSeq,RptMskUpperCase);
		else
			{
			strncpy((char *)pRetSeq,&m_pszFastqSeq[m_FastqSeqIdx],SeqLen);
			((char *)pRetSeq)[SeqLen] = '\0';
			}
		m_FastqSeqIdx += SeqLen;
		return(SeqLen);
		}
	m_FastqSeqIdx = 0;
	if((Rslt = ParseFastQblockQ()) <= eBSFSuccess)
		return(Rslt);
	m_bDescrAvail = true;
	return(Rslt);
	}

m_bDescrAvail = false;
bInDescriptor = false;
bMoreToDo = true;
bSloughEOL = false;
PrevSOLiDbase = 0;

while(bMoreToDo) {
	if (m_pCurFastaBlock->BuffIdx >= m_pCurFastaBlock->BuffCnt)	// time to refill m_pBuffer with another (up to) m_BuffSize chars?
		{
		if(m_gzFile != NULL)
			{
			FileOfs = gztell(m_gzFile);
			m_pCurFastaBlock->BuffCnt = gzread(m_gzFile, m_pCurFastaBlock->pBlock, m_pCurFastaBlock->AllocSize);
			}
		else
			{
			FileOfs = _lseeki64(m_hFile,0,SEEK_CUR);
			m_pCurFastaBlock->BuffCnt = read(m_hFile, m_pCurFastaBlock->pBlock, m_pCurFastaBlock->AllocSize);
			}

		if (m_pCurFastaBlock->BuffCnt <= 0)
			break;
		m_pCurFastaBlock->BuffIdx = 0;
		m_pCurFastaBlock->FileOfs = FileOfs;
		}

	while (m_pCurFastaBlock->BuffIdx < m_pCurFastaBlock->BuffCnt)
		{
		Chr = m_pCurFastaBlock->pBlock[m_pCurFastaBlock->BuffIdx++];
		// ensure reading an ascii text file - only allow whitespace and chrs >= 0x20 and <= 0x7f
		// note that if within a descriptor line then chars > 0x7f are tolerated but will be substituted with '?' 
		if(bInDescriptor && (uint8_t)Chr > 0x7f)
			Chr = '?';
		
		if((Chr < 0x20 || (uint8_t)Chr > 0x7f) && !(Chr == '\n' || Chr == '\r' || Chr == '\t' || Chr == ' '))
			{
			AddErrMsg("CFasta::ReadSequence","Errors whilst reading file - '%s' - illegal chr 0x%x", m_szFile,Chr);
			return(eBSFerrFileAccess);
			}
		if(bSloughEOL)							// slough to end of line?
			{
			if(Chr == '\n')			
				bSloughEOL = false;
			continue;
			}

		if(bInDescriptor)						// parse out Descriptor lines
			{
			// slough any descriptor prefix consisting of spaces or tab chars
			if(m_DescriptorLen == 0 && (Chr == ' ' || Chr == '\t'))
				continue;

			if((Chr == '\n' || Chr == '\r'))	// end of Descriptor?
				{
				m_szDescriptor[m_DescriptorLen] = '\0';
				m_bDescrAvail = true;
				return(eBSFFastaDescr);		 // let caller know current sequence has ended and a new descriptor line is available
				}
			if(m_DescriptorLen < cMaxFastaDescrLen)
				m_szDescriptor[m_DescriptorLen++] = Chr;
			continue;
			}
		
		switch(Chr) {
			case '#':					// if not in a descriptor line then can accept if csfasta
				if(m_bIscsfasta == true)
					bSloughEOL = true;	
				continue;


			case '>':					// start of a new fasta Descriptor
				if(SeqLen == 0 && m_bForceRetSeq)       // ensuring that there will always be a simulated sequence of length 1 returned if no actual sequence between two descriptors
					{
					m_pCurFastaBlock->BuffIdx--;		// effective pushback of descriptor..
					m_bForceRetSeq = false;
					if(pRetSeq != NULL)
						*(uint8_t *)pRetSeq = bSeqBase ? eBaseN : 'N';
					return(1);
					}
				if(SeqLen)				// if sequence avail then return sequence, parse descriptor on next call..
					{
					m_bForceRetSeq = false;
					m_pCurFastaBlock->BuffIdx--;		// effective pushback of descriptor..
					if(pRetSeq != NULL && bSeqBase)
						return(Ascii2Sense((char *)pRetSeq,SeqLen,(etSeqBase *)pRetSeq,RptMskUpperCase));
					return(SeqLen);
					}
				bInDescriptor = true;	// parse out descriptor
				m_DescriptorLen = 0;
				m_szDescriptor[0] = '\0';
				m_FileDescrOfs = m_pCurFastaBlock->FileOfs + m_pCurFastaBlock->BuffIdx - 1; // note file offset at which descriptor starts
				continue;


			default:
				if(m_bIscsfasta)		// csfasta processing
					{
					if(pRetSeq != NULL) {
						switch(Chr) {
							case 'a': case 'A':
								PrevSOLiDbase = eBaseA;
								continue;
							case 'c': case 'C':
								PrevSOLiDbase = eBaseC;
								continue;
							case 'g': case 'G':
								PrevSOLiDbase = eBaseG;
								continue;
							case 't': case 'T':
								PrevSOLiDbase = eBaseT;
								continue;
							case 'n': case 'N':			// should never see a SOLiD sequence starting on an indeterminate base
								PrevSOLiDbase = eBaseT; // use a default - will be correct 25% of the time :-)
								continue;

							case '0': case '1': case '2': case '3':
								Chr = m_SOLiDmap[PrevSOLiDbase][Chr-'0'];
								switch(Chr) {
									case 'a': case 'A':
										PrevSOLiDbase = eBaseA;
										break;
									case 'c': case 'C':
										PrevSOLiDbase = eBaseC;
										break;
									case 'g': case 'G':
										PrevSOLiDbase = eBaseG;
										break;
									case 't': case 'T':
										PrevSOLiDbase = eBaseT;
										break;
									case 'n': case 'N':  // although a known indeterminate, take a 25% chance that the indeterminate
										PrevSOLiDbase = eBaseA; // was an eBaseA 
										break;			 
									}
								break;

							case '4':			// although a known indeterminate, take a 25% chance that the indeterminate
							case '.':			// was an eBaseA, otherwise the rest of sequence would be
								Chr = 'n';		// treated as all bases were indeterminate
								PrevSOLiDbase = eBaseA;
								break;
							}
						}
					else
						Chr = 'n';
					}
				if(Chr != '-' && !isalpha(Chr))		// slough any non-alpha except for '-' outside of Descriptors
					continue;

				if(pRetSeq != NULL)
					{
					*pAscii++ = Chr;
					if(++SeqLen < Max2Ret)
						*pAscii='\0';
					}
				else
					SeqLen++;

				if(pRetSeq == NULL || Max2Ret == 0 || SeqLen < Max2Ret)
					continue;

				m_bForceRetSeq = false;
				if(pRetSeq != NULL && bSeqBase)
					return(Ascii2Sense((char *)pRetSeq,SeqLen,(etSeqBase *)pRetSeq,RptMskUpperCase));
				return(SeqLen);
			}
		}
	}


if (m_pCurFastaBlock->BuffCnt < 0)
	{
	AddErrMsg("CFasta::ReadSequence","Errors whilst reading file - '%s' - %s", m_szFile, strerror(errno));
	return(eBSFerrFileAccess);
	}
if(SeqLen)
	{
	m_bForceRetSeq = false;
	if(pRetSeq != NULL && bSeqBase)
		return(Ascii2Sense((char *)pRetSeq,SeqLen,(etSeqBase *)pRetSeq,RptMskUpperCase));
	return(SeqLen);
	}
if(bInDescriptor)
	{
	m_bForceRetSeq = false;
	m_bDescrAvail = true;
	m_szDescriptor[m_DescriptorLen] = '\0';
	return(eBSFFastaDescr);
	}
return(eBSFSuccess);
}


// LocateDescriptor
// locates a descriptor whose text matches that specified
// If bFromStart==false then the search starts from the previous descriptor (or start of file if first call)
// If bFromStart==true then the search starts from the start of file
// Returns eBSFSuccess if descriptor is located
int32_t													// eBSFSuccess if matched and descriptor is available
CFasta::LocateDescriptor(char *pszPrefixToMatch,	// prefix to match or NULL if match any
						 bool bFromStart)			// true: from start, false: from last descriptor or sequence
{
char Buffer[16000];
int32_t Cnt;
int32_t CmpLen;
int32_t FilePsn;

if(m_gzFile == NULL && m_hFile == -1 )
	return(eBSFerrClosed);
if(!m_bRead)
	return(eBSFerrRead);

if(bFromStart)
	{
	m_pCurFastaBlock->BuffCnt = 0;
	m_pCurFastaBlock->BuffIdx = 0;
	m_bDescrAvail = false;
	m_FastqSeqLen = 0;
	m_FastqSeqIdx = 0;
	m_FastqSeqQLen = 0;

	if(m_gzFile != NULL)
		FilePsn = gzseek(m_gzFile,0,SEEK_SET);
	else
		FilePsn = (int32_t)_lseeki64(m_hFile,0,SEEK_SET);
	if(FilePsn != 0)
		{
		AddErrMsg("CFasta::LocateDescriptor","Seek failed to offset 0 on %s - %s",m_szFile,strerror(errno));
		return(eBSFerrFileAccess);
		}
	}

if(pszPrefixToMatch != NULL && pszPrefixToMatch[0] != '\0')
	CmpLen = (int32_t)strlen(pszPrefixToMatch);
else
	CmpLen = 0;
Cnt = 0;
while(m_bDescrAvail || (Cnt = ReadSequence(Buffer,sizeof(Buffer)))>=0)
	{
	m_bDescrAvail = false;
	if(Cnt > 0 && Cnt != eBSFFastaDescr)			// slough sequences, looking for descriptor lines...
		continue;	
	if(!CmpLen || !strnicmp(m_szDescriptor,pszPrefixToMatch,CmpLen))
		{
		m_bDescrAvail = true;
		return(eBSFSuccess);
		}
	}
return(eBSFerrFastaDescr);
}


// ParseFastQblockQ
// Parses a fastq block (seq identifier + sequence + quality scores)
//			fastq can also contain SOLiD 2 base representation
// Returns eBSFSuccess if EOF
//         eBSFFastaDescr if fastq block was parsed
//         standard error codes if any errors
int32_t						
CFasta::ParseFastQblockQ(void)	
{
char Chr;
int64_t FileOfs;			// will contain file offset corresponding to last buffer read from file
int32_t ParseState;
int32_t SeqLen = 0;
bool bIsFastQSOLiD;	// some fastq files (from NCBA SRA SRP000191) have SOLiD sequences
char PrvBase;		// used if decoding SOLiD sequences
if (m_gzFile == NULL && m_hFile == -1 || m_pCurFastaBlock == NULL || m_pCurFastaBlock->pBlock == NULL)
	return(eBSFerrClosed);
if(!m_bRead)
	return(eBSFerrRead);

m_FastqSeqLen = 0;
m_FastqSeqIdx = 0;
m_FastqSeqQLen = 0;
bIsFastQSOLiD = false;
ParseState = 0;
while(ParseState < 6) {
	if (m_pCurFastaBlock->BuffIdx >= m_pCurFastaBlock->BuffCnt)
		{
		if(m_gzFile != NULL)
			{
			FileOfs = gztell(m_gzFile);
			m_pCurFastaBlock->BuffCnt = gzread(m_gzFile, m_pCurFastaBlock->pBlock, m_pCurFastaBlock->AllocSize);
			}
		else
			{
			FileOfs = _lseeki64(m_hFile,0,SEEK_CUR);
			m_pCurFastaBlock->BuffCnt = read(m_hFile, m_pCurFastaBlock->pBlock, m_pCurFastaBlock->AllocSize);
			}

		if (m_pCurFastaBlock->BuffCnt <= 0)
			break;
		m_pCurFastaBlock->BuffIdx = 0;
		m_pCurFastaBlock->FileOfs = FileOfs;
		}
	while (ParseState < 6 && m_pCurFastaBlock->BuffIdx < m_pCurFastaBlock->BuffCnt)
		{
		Chr = m_pCurFastaBlock->pBlock[m_pCurFastaBlock->BuffIdx++];
		if(Chr == '\n')
			m_CurFastQParseLine += 1;

		// ensure reading an ascii text file - only allow whitespace and chrs >= 0x20 and <= 0x7f
		if(!isspace(Chr) &&  (Chr < 0x20 || (uint8_t)Chr > 0x7f))
			{
			AddErrMsg("CFasta::ParseFastQblockQ","Errors whilst reading file - '%s' - near line %d, illegal chr 0x%x", m_szFile,m_CurFastQParseLine+1,Chr);
			return(eBSFerrFastqChr);
			}

		switch(ParseState) {
			case 0:		// expecting to read new block
				if(isspace(Chr)) // slough whitespace between blocks
					continue;
				if(Chr == '@')	// marks start of sequence identifier (descriptor) line
					{
					ParseState = 1;
					m_DescriptorLen = 0;
					m_FastqSeqLen = 0;
					m_FastqSeqQLen = 0;
					m_FileDescrOfs = m_pCurFastaBlock->FileOfs + m_pCurFastaBlock->BuffIdx - 1; // note file offset at which descriptor starts
					continue;
					}
				AddErrMsg("CFasta::ParseFastQblockQ","Errors whilst reading fastq file - '%s' - near line %d, expected '@' sequence identifer line, instead read 0x%x", m_szFile,m_CurFastQParseLine+1,Chr);
				return(eBSFerrFastqSeqID);

			case 1:		// parsing sequence identifier
				if(Chr == '\n' || Chr == '\r')	// end of Descriptor?
					{
					m_szDescriptor[m_DescriptorLen] = '\0';
					ParseState = 2;					// next expecting the sequence
					bIsFastQSOLiD = false;
					}
				else
					if(m_DescriptorLen < cMaxFastaDescrLen)
						m_szDescriptor[m_DescriptorLen++] = Chr;
				continue;

			case 2:		// parsing sequence - bit of a hack here as need to be able to handle the SOLiD sequences as well as the usual 'a','c','g','t' and 'n's
				if(Chr == '\n' || Chr == '\r')	// end of sequence?
					{
					if(m_FastqSeqLen == 0)		
						continue;
					m_pszFastqSeq[m_FastqSeqLen] = '\0';
					ParseState = 3;					// next expecting duplicate sequence identifier
					continue;
					}

				if(m_FastqSeqLen >= cMaxFastQSeqLen)	// hacky... simply silently truncate overlength sequences
					continue;

				switch(Chr) {
					case '+':			// a very special case! very rarely the sequence may have been completely filtered out so need to check for this and if no sequence then pretend there was a single indeterminate base
										// if we were to be strict and simply treat as an error then this could upset downstream processing whereby PE reads were being processed and a mate read was being expected
										// a warning is written to log
						if(m_FastqSeqLen == 0)
							{
							m_pszFastqSeq[m_FastqSeqLen++] = 'N';
							m_pszFastqSeq[m_FastqSeqLen] = '\0';
							m_pCurFastaBlock->BuffIdx--;			// push back so '+' will be reparsed as if true duplicate sequence identifier
							if(m_NumMissingSequences++ == 0)
								gDiagnostics.DiagOut(eDLWarn,gszProcName,"ParseFastQblockQ: At least one instance of missing sequence in file '%s' - first is near line %d, simulating sequence of length 1", m_szFile,m_CurFastQParseLine+1);
							ParseState = 3;		
							continue;
							}
						AddErrMsg("CFasta::ParseFastQblockQ","sequence char (0x%2x %c) error whilst reading fastq file - '%s' - near line %d at chr position %d", (int32_t)Chr,Chr, m_szFile,m_CurFastQParseLine+1,m_FastqSeqLen);
						return(eBSFerrFastqSeq);

					case 'a': case 'A': case 'c': case 'C': case 'g': case 'G': case 't': case 'T': case 'n': case 'N':
						if(bIsFastQSOLiD)		// if already determined that SOLiD is being parsed then not expecting acgtn's
							{
							AddErrMsg("CFasta::ParseFastQblockQ","Unexpected SOLiD sequence type error whilst reading fastq file - '%s' - near line %d", m_szFile,m_CurFastQParseLine+1);
							return(eBSFerrFastqSeq);
							}
						break;

					case '0': case '1': case '2': case '3': case '4': case '.':  // '4' and '.' represents an indeterminate call
						if(m_FastqSeqLen > 1 && !bIsFastQSOLiD)
							{
							AddErrMsg("CFasta::ParseFastQblockQ","Unexpected SOLiD sequence type error whilst reading fastq file - '%s' - near line %d", m_szFile,m_CurFastQParseLine+1);
							return(eBSFerrFastqSeq);
							}

						if(m_FastqSeqLen == 1 && !bIsFastQSOLiD)
							{
							bIsFastQSOLiD = true;
							m_FastqSeqLen = 0;
							PrvBase = tolower(m_pszFastqSeq[0]);
							}
						else
							{
							if((PrvBase = m_pszFastqSeq[m_FastqSeqLen-1])=='n')
								PrvBase = 'a'; // making previous base as not being indeterminate gives a 25% chance remainder of sequence will be correct
							}

						if(Chr == '.' || Chr == '4')	
							Chr = 'n';
						else
							switch(PrvBase) {
								case 'a':
									switch(Chr) {
										case '0': Chr = 'a'; break;
										case '1': Chr = 'c'; break;
										case '2': Chr = 'g'; break;
										case '3': Chr = 't'; break;
										}
									break;
								case 'c':
									switch(Chr) {
										case '1': Chr = 'a'; break;
										case '0': Chr = 'c'; break;
										case '3': Chr = 'g'; break;
										case '2': Chr = 't'; break;
										}
									break;
								case 'g':
									switch(Chr) {
										case '2': Chr = 'a'; break;
										case '3': Chr = 'c'; break;
										case '0': Chr = 'g'; break;
										case '1': Chr = 't'; break;
										}
									break;
								case 't':
									switch(Chr) {
										case '3': Chr = 'a'; break;
										case '2': Chr = 'c'; break;
										case '1': Chr = 'g'; break;
										case '0': Chr = 't'; break;
										}
									break;
								}
						break;
	
					default:
						AddErrMsg("CFasta::ParseFastQblockQ","sequence char (0x%2x %c) error whilst reading fastq file - '%s' - near line %d at chr position %d", (int32_t)Chr,Chr, m_szFile,m_CurFastQParseLine+1,m_FastqSeqLen);
						return(eBSFerrFastqSeq);
					}
					
				m_pszFastqSeq[m_FastqSeqLen++] = Chr;
				continue;

			case 3:		// parsing '+' duplicate sequence identifier
				if(Chr == '\n' || Chr == '\r')
					continue;
				if(Chr == '+')
					{
					ParseState = 4;
					continue;
					}
				AddErrMsg("CFasta::ParseFastQblockQ","Errors whilst reading fastq file - '%s' - near line %d, expected '+' sequence identifer line, instead read 0x%x", m_szFile,m_CurFastQParseLine+1,Chr);
				return(eBSFerrFastqDescr);


			case 4:		// parsing duplicate sequence identifier
				if(!(Chr == '\n' || Chr == '\r'))	// slough duplicate, end of identifier?
					continue;
				ParseState = 5;			// next should be the quality scores
				continue;

			case 5:		// parsing quality scores
				if(Chr == '\n' || Chr == '\r')	// end of quality scores?
					{
					if(m_FastqSeqQLen == 0)
						{
						if(m_FastqSeqLen > 1)
							continue;
						m_pszFastqSeqQ[m_FastqSeqQLen++] = 'a';
						}

					m_pszFastqSeqQ[m_FastqSeqQLen] = '\0';
					ParseState = 6;					// block processed
					}
				else
					{
					if(m_FastqSeqQLen >= cMaxFastQSeqLen)	// silently slough overlength quality scores (same as overlength sequences)
						continue;
					if(bIsFastQSOLiD && !m_FastqSeqQLen)    // slough 1st quality score on SOLiD reads as that applies to the 1st bases which is a ligator sequence base which is also sloughed
						{
						bIsFastQSOLiD = false;
						continue;
						}
					m_pszFastqSeqQ[m_FastqSeqQLen++] = Chr;
					}
				continue;
			}
		}
	}
if (m_pCurFastaBlock->BuffCnt < 0)
	{
	AddErrMsg("CFasta::ParseFastQblockQ","Errors whilst reading fastq file - '%s' near line %d, - %s", m_szFile, m_CurFastQParseLine+1, strerror(errno));
	return(eBSFerrFileAccess);
	}

if(ParseState == 0)
	return(eBSFSuccess);

if(ParseState == 5)
	{
	if(m_FastqSeqQLen > 0)	// need to allow for last block having quality scores not being terminated by a CR/LF as last block could be EOF terminated
		{
		m_pszFastqSeqQ[m_FastqSeqQLen] = '\0';
		ParseState = 6;					// treat as if block was processed
		}
	else
		if(m_FastqSeqLen == 1 && m_FastqSeqQLen == 0)	// could have been a filtered out sequence which is still present as a fastq sequence block
			{
			m_pszFastqSeqQ[m_FastqSeqQLen++] = 'a';
			m_pszFastqSeqQ[m_FastqSeqQLen] = '\0';
			ParseState = 6;					// treat as if block was processed
			}
	}

// check all elements were present and not empty
if(ParseState != 6 || m_DescriptorLen == 0 || m_FastqSeqLen == 0 && m_FastqSeqQLen == 0)
	{
	AddErrMsg("CFasta::ParseFastQblockQ","Errors whilst reading fastq file block, empty elements - '%s' near line %d", m_szFile,m_CurFastQParseLine+1);
	AddErrMsg("CFasta::ParseFastQblockQ","Parse state: %d, Descriptor length: %d, Seq length: %d, Qual length: %d", ParseState,m_DescriptorLen,m_FastqSeqLen,m_FastqSeqQLen);
	return(eBSFerrFileAccess);
	}
if(m_FastqSeqLen != m_FastqSeqQLen)
	{
	AddErrMsg("CFasta::ParseFastQblockQ","Expected number of quality scores (%d) to be same as number of bases (%d) - '%s' near line %d", m_FastqSeqQLen,m_FastqSeqLen,m_szFile,m_CurFastQParseLine+1);
	return(eBSFerrFileAccess);
	}
return(eBSFFastaDescr);
}



// ReadSubsequence
// Reads a subsubsequence from fasta starting at SeqOfs within the next sequence to be located in the fasta file
int32_t										// number actually read
CFasta::ReadSubsequence(etSeqBase *pSeq,	// where to return subsequence (NO TERMINATING eBaseEOS)
				 uint32_t MaxLen,			// reads upto MaxLen sequence bases from file
				 uint32_t SeqOfs,			// relative offset in sequence at which to read
				 bool bFromStart)			// true: from start of fasta file, false: from start of currently buffered sequence
{
char Buffer[0x2fff];
uint32_t CurSeqOfs = 0;
uint32_t SeqLen = 0;
int32_t Cnt;
int32_t CpyLen;
int32_t CpyFromIdx;

if (m_gzFile == NULL && m_hFile == -1 || m_pCurFastaBlock == NULL || m_pCurFastaBlock->pBlock == NULL)
	return(eBSFerrClosed);

if(!m_bRead)
	return(eBSFerrRead);

if(bFromStart) // if from start of source Fasta file
	{
	int32_t FileOfs;
	m_pCurFastaBlock->BuffCnt = 0;
	m_pCurFastaBlock->BuffIdx = 0;
	m_bDescrAvail = false;
	if(m_gzFile != NULL)
		FileOfs = gzseek(m_gzFile,0,SEEK_SET);
	else
		FileOfs = (int32_t)_lseeki64(m_hFile,0,SEEK_SET);

	if(FileOfs!=0)
		{
		AddErrMsg("CFasta::ReadSubsequence","Seek failed to offset 0 on %s - %s",m_szFile,strerror(errno));
		return(eBSFerrFileAccess);
		}
	}

m_bDescrAvail = false;
CurSeqOfs = 0;
while(SeqLen < MaxLen && (Cnt = ReadSequence(Buffer,sizeof(Buffer),true)) > 0)
	{
	if(Cnt == eBSFFastaDescr)		// slough descriptors, looking for a sequence...
		{
		if(SeqLen)		// any descriptor terminates the currently accumulating sequence
			{
			m_bDescrAvail = true;
			break;
			}
		continue; // no currently accumulating sequence so bypass descriptor
		}
	if(Cnt + CurSeqOfs < (int32_t)SeqOfs)
		{
		CurSeqOfs += Cnt;
		continue;
		}
	if(CurSeqOfs < SeqOfs)
		CpyFromIdx = SeqOfs - CurSeqOfs;
	else
		CpyFromIdx = 0;
	CpyLen = MaxLen - SeqLen;
	if(CpyLen > Cnt)
		CpyLen = Cnt;
	memmove(pSeq,&Buffer[CpyFromIdx],CpyLen);
	pSeq += CpyLen;
	SeqLen += CpyLen;
	}
return(SeqLen);
}

char 
CFasta::Base2Chr(etSeqBase Base)		 // returns ASCII representation of Base
{
Base &= 0x0f;
switch(Base) {
	case eBaseA:
		return('a');
	case eBaseC:
		return('c');
	case eBaseG:
		return('g');
	case eBaseT:
		return('t');
	}
return('n');
}

etSeqBase 
CFasta::Chr2Base(char Base)			 // returns etSeqBase representation of ASCII base
{
switch(Base) {
	case 'a': case 'A':
		return(eBaseA);
	case 'c': case 'C':
		return(eBaseC);
	case 'g': case 'G':
		return(eBaseG);
	case 't': case 'T': case 'u': case 'U':
		return(eBaseT);
	}
return(eBaseN);
}

// Ascii2Sense
// Translates ascii into etSeqBase's
// Caller can override assumption that lowercase represents softmasked repeats
int32_t
CFasta::Ascii2Sense(char *pAscii,		// expected to be '\0' terminated, or SeqLen long
					int32_t MaxSeqLen,		// maximal sized sequence that pSeq can hold				
					etSeqBase *pSeq,	// where to return translated bases
					bool RptMskUpperCase) // true if bases are softmasked as uppercase instead of default lowercase
{
char Base;
int32_t SeqLen = 0;
while(MaxSeqLen-- && (Base = *pAscii++)!='\0')
	{
	SeqLen++;
	switch(Base) {
		case 'a':
			*pSeq++ = eBaseA | (RptMskUpperCase ? 0 : cRptMskFlg);
			continue;
		case 'A':
			*pSeq++ = eBaseA  | (RptMskUpperCase ? cRptMskFlg : 0);
			continue;
		case 'c':
			*pSeq++ = eBaseC  | (RptMskUpperCase ? 0 : cRptMskFlg);
			continue;
		case 'C':
			*pSeq++ =  eBaseC | (RptMskUpperCase ? cRptMskFlg : 0);
			continue;
		case 'g':
			*pSeq++ =  eBaseG  | (RptMskUpperCase ? 0 : cRptMskFlg);
			continue;
		case 'G':
			*pSeq++ =  eBaseG | (RptMskUpperCase ? cRptMskFlg : 0);
			continue;
		case 't': case 'u': 
			*pSeq++ =  eBaseT  | (RptMskUpperCase ? 0 : cRptMskFlg);
			continue;
		case 'T': case 'U':
			*pSeq++ =  eBaseT | (RptMskUpperCase ? cRptMskFlg : 0);
			continue;

		case '-':
			*pSeq++ =  eBaseInDel;
			continue;

		default: // 'N'
			*pSeq++ = eBaseN;
			continue;
		}
	}

return(SeqLen);
}

// ReadDescriptor
// Copies last parsed descriptor (truncates at MaxLen) into user supplied buffer.
// Note: ReadSequence() is expected to return a sequence subsequent to a ReadDescriptor() call.
//       Some fasta/fastq filtering processes can delete the complete sequence leaving descriptors immediately followed by another descriptor!
//       Calls to ReadSequence() with user expecting a sequence to be returned would then fail, so a flag is set forcing ReadSequence() to always return at least 1bp immediately after a call to ReadDescriptor!!!
int32_t							// returns strlen of available descriptor or 0 if none
CFasta::ReadDescriptor(char *pszDescriptor,int32_t MaxLen)
{
if(m_gzFile == NULL && m_hFile == -1)
	return(eBSFerrClosed);

if(!m_bRead)
	return(eBSFerrRead);

if(!m_DescriptorLen || !m_bDescrAvail)
	return(eBSFerrParams);

strncpy(pszDescriptor,m_szDescriptor,MaxLen);
if(m_DescriptorLen >= (uint32_t)MaxLen)
	{
	pszDescriptor[MaxLen-1] = '\0';
	return(MaxLen-1);
	}
m_FileReadDescrOfs = m_FileDescrOfs;
m_bForceRetSeq = true; // ensures a dummy sequence of length 1 will always be returned if a sequence missing between two descriptors - Ugh, Ugh and Ugh again
return(m_DescriptorLen);
}

// ReadQualityValues
// Copies last parsed sequence quality scores (truncates at MaxLen) into user supplied buffer.
int32_t							// returns strlen of available quality scores or 0 if none
CFasta::ReadQValues(char *pszValues,int32_t MaxLen)
{
if(m_gzFile == NULL && m_hFile == -1)
	return(eBSFerrClosed);

if(!m_bRead)
	return(eBSFerrRead);

if(!m_bDescrAvail)
	return(eBSFerrParams);

if(m_FastqSeqQLen)
	{
	strncpy(pszValues,m_pszFastqSeqQ,MaxLen);
	if(m_FastqSeqQLen >= MaxLen)
		{
		pszValues[MaxLen-1] = '\0';
		return(MaxLen-1);
		}
	}
return(m_FastqSeqQLen);
}


// returns file offset at which descriptor returned by ReadDescriptor() was parsed from
int64_t 
CFasta::GetDescrFileOfs(void)				
{
return(m_FileReadDescrOfs);
}


// Write
// Write ascii (assumed to be a IUB symbol) IUB char to file
// Note that whitespace is sloughed w/o any error being returned
// Note that any other nonalpha chr is treated as an error
int32_t 
CFasta::Write(char Symbol)
{
int32_t NumWritten;

if(m_hFile == -1)
	return(eBSFerrClosed);

if(m_bRead)
	return(eBSFerrWrite);

if(isspace(Symbol))
	return(eBSFSuccess);

if(!isalpha(Symbol))
	return(eBSFerrFastaSymb);

if (m_pCurFastaBlock->BuffIdx >= m_pCurFastaBlock->AllocSize - 2)		// 2 is safety margin to account for any nl's etc
	{
	NumWritten = write(m_hFile, m_pCurFastaBlock->pBlock, m_pCurFastaBlock->BuffIdx);
	if (NumWritten != m_pCurFastaBlock->BuffIdx)
		{
		AddErrMsg("CFasta::Write", "Errors whilst writing %d bytes - '%s' - %s", m_pCurFastaBlock->BuffIdx, m_szFile, strerror(errno));
		return(eBSFerrFileAccess);
		}
	m_pCurFastaBlock->BuffIdx = 0;
	}
if(m_CurLineLen >= cMaxGenFastaLineLen)
	{
	m_pCurFastaBlock->pBlock[m_pCurFastaBlock->BuffIdx++] = '\n';
	m_CurLineLen = 0;
	}
m_pCurFastaBlock->pBlock[m_pCurFastaBlock->BuffIdx++] = Symbol;
m_CurLineLen++;
return(eBSFSuccess);
}

// Write Cnt ascii (assumed IUB symbols w/o any checking) IUB chars to file
int32_t
CFasta::Write(char *pSymbols,uint32_t Cnt)
{
int32_t Rslt = eBSFSuccess;
if(m_bRead)
	return(eBSFerrWrite);

if(m_hFile == -1)
	return(eBSFerrClosed);

while(Cnt-- && Rslt == eBSFSuccess)
	Rslt = Write(*pSymbols++);
return(Rslt);
}

// Write fasta descriptor line to file
// 
int32_t
CFasta::WriteDescriptor(char *pszDescriptor)
{
size_t Len;
int32_t NumWritten;
if(m_bRead)
	return(eBSFerrWrite);

if(m_hFile == -1)
	return(eBSFerrClosed);

if(pszDescriptor == NULL || *pszDescriptor == '\0')
	return(eBSFerrParams);
// truncate overlength descriptor lines
Len = strlen(pszDescriptor);
if(Len > cMaxFastaDescrLen)
	Len = cMaxFastaDescrLen;

if (m_pCurFastaBlock->BuffIdx >= m_pCurFastaBlock->AllocSize - (int32_t)Len - 16)  // 16 is a safety margin to allow for nl's etc
	{
	NumWritten = write(m_hFile, m_pCurFastaBlock->pBlock, m_pCurFastaBlock->BuffIdx);
	if (NumWritten != m_pCurFastaBlock->BuffIdx)
		{
		AddErrMsg("CFasta::WriteDescriptor", "Errors whilst writing %d bytes - '%s' - %s", m_pCurFastaBlock->BuffIdx, m_szFile, strerror(errno));
		return(eBSFerrFileAccess);
		}
	m_pCurFastaBlock->BuffIdx = 0;
	}

if(m_CurLineLen)
	m_pCurFastaBlock->pBlock[m_pCurFastaBlock->BuffIdx++] = '\n';

if(*pszDescriptor != '>')
	m_pCurFastaBlock->pBlock[m_pCurFastaBlock->BuffIdx++] = '>';
strncpy((char *)&m_pCurFastaBlock->pBlock[m_pCurFastaBlock->BuffIdx], pszDescriptor, Len);
m_pCurFastaBlock->BuffIdx += (int32_t)Len;
m_pCurFastaBlock->pBlock[m_pCurFastaBlock->BuffIdx++] = '\n';
m_CurLineLen=0;
return(eBSFSuccess);
}

int32_t
CFasta::Close(void)
{
if(!m_bRead)
	{
	if(m_hFile == -1)		// if already closed then treat as success
		return(eBSFSuccess);

	if (m_pCurFastaBlock->BuffIdx)
		if (m_pCurFastaBlock->BuffIdx != write(m_hFile, m_pCurFastaBlock->pBlock, m_pCurFastaBlock->BuffIdx))
			{
			AddErrMsg("CFasta::Close", "Errors whilst writing %d bytes - '%s' - %s", m_pCurFastaBlock->BuffIdx, m_szFile, strerror(errno));
			return(eBSFerrFileAccess);
			}
	}

Cleanup();
return(eBSFSuccess);
}
