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
#include <sys/resource.h>
#endif

uint16_t			// generated 16bit hash over the lowercased chromosome name; hashes as 0 if pszName == null or is empty
CUtility::GenHash16(char *pszName)	// name will be lowercased whilst hashing
{
int Hash;
char Chr;
Hash = 19937;			// circular prime as hash seed

if(pszName == NULL || pszName[0] == '\0')
	return(0);

while(Chr = *pszName++)
	{
	Hash = (Hash ^ (int)tolower(Chr)) * 3119;	// a circular prime
	Hash ^= (Hash >> 13);
	Hash &= 0x0ffff;
	}
if(Hash == 0)			// 0 reserved as an error indicator so force hash to be 19937
	Hash = 19937;
return((uint16_t)(Hash & 0x0ffff));
}

int						// // generated 24bit hash over the lowercased chromosome name; hashes as 0 if pszName == null or is empty
CUtility::GenHash24(char *pszName) // name will be lowercased whilst hashing
{
int Hash;
char Chr;
Hash = 37199;			// circular prime as hash seed

if(pszName == NULL || pszName[0] == '\0')
	return(0);

while(Chr = *pszName++)
	{
	Hash = (Hash ^ (int)tolower(Chr)) * 19937;	// a circular prime
	Hash ^= (Hash >> 21);
	Hash &= 0x0ffffff;
	}
if(Hash == 0)			// 0 reserved as an error indicator so force hash to be 1
	Hash += 1;
return(Hash);
}

// ClassifyFileType
// Attempt to classify the alignment file as one of CSV, BED or SAM/BAM from it's initial 8k char contents
// Currently processes CSV, BED and SAM format file types
// Assumes must be SAM if initial lines have at least one prefixed by a '@' followed by a 2 letter record type code 
//
etClassifyFileType
CUtility::ClassifyFileType(char *pszFileName)
{
int hFile;
gzFile gz;
BGZF* pInBGZF;
int BuffLen;
int BuffIdx;
uint8_t Buffer[cFileClassifyBuffLen];
uint8_t *pBuff;
bool bStartNL;
bool bSkipEOL;
uint8_t Chr;
int NumLines;
int FieldCnt;
int TabCnt;
int CommaCnt;
int FldLen;
bool bInQuotes;
int LikelyCSV;
int LikelyBED;
int LikelySAM;
int LikelyNonCSVSAMBED;
bool bSeenSAMHdrs;
int FileNameLen;
bool bGZd;

FileNameLen = (int)strlen(pszFileName);
bGZd = false;
if(FileNameLen >= 4)
	{
	if(!stricmp(&pszFileName[FileNameLen-3],".gz"))
		bGZd = true;
	else
		{
		if(FileNameLen >= 5 && !stricmp(&pszFileName[FileNameLen-4],".bam"))
			{
			hFile = open(pszFileName,O_READSEQ);
			if(hFile == -1)
				return(eCFTopenerr);
			// BAM will using BGZF compression ..
			if((pInBGZF = bgzf_dopen(hFile, "r"))==NULL)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"ClassifyFileType: unable to initialise for BGZF processing on file '%s'",pszFileName);
				close(hFile);
				return(eCFTopenerr);
				}
			hFile = -1;

			// try reading the header, bgzf_read will confirm it does start with "BAM\1" ....
			if((BuffLen = (int)bgzf_read(pInBGZF,Buffer,100)) < 100)		// will be < 100 if errors ...
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"ClassifyFileType: Not a BAM format file '%s'",pszFileName);
				bgzf_close(pInBGZF);
				return(eCFTopenerr);
				}
			bgzf_close(pInBGZF);
			return(eCFTSAM);
			}
		}
	}

// now can try to actually open file and read in first cFileClassifyBuffLen chars
if(bGZd)
	{
	gz = gzopen(pszFileName,"rb");
	if(gz == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Open: unable to open for reading gzip'd file '%s'",pszFileName);
		return(eCFTopenerr);
		}
	BuffLen = gzread(gz,Buffer,sizeof(Buffer)-1);
	gzclose(gz);
	}
else
	{
	hFile = open(pszFileName,O_READSEQ);
	if(hFile == -1)
		return(eCFTopenerr);
	// read the 1st cFileTypeBuffLen into buffer
	BuffLen = read(hFile,Buffer,sizeof(Buffer)-1);
	close(hFile);
	}

if(BuffLen < cMinFileClassifyLen)	// an arbitary lower limit!
	return(eCFTlenerr);

Buffer[BuffLen] = '\0';
pBuff = Buffer;
NumLines = 0;
LikelyCSV = 0;
LikelyBED = 0;
LikelySAM = 0;
LikelyNonCSVSAMBED = 0;
BuffIdx = 0;
bStartNL = true;
bSeenSAMHdrs = false;
while(Chr = *pBuff++)
	{
	BuffIdx += 1;
	if(bStartNL)
		{
		FieldCnt = 0;
		TabCnt = 0;
		CommaCnt = 0;
		FldLen = 0;
		bInQuotes = false;
		bStartNL = false;
		bSkipEOL = false;
		NumLines += 1;
		}
	if(Chr == '\n' || Chr == '\r')			// if at end of line
		{
		bStartNL = true;
		bSkipEOL = false;
		if(FieldCnt < 3)					// BED can have down to 3 fields, CSV alignment and SAM should have more
			continue;

		if(!bInQuotes)
			{
			if(CommaCnt >= 3 && CommaCnt > TabCnt)		// if at least as many commas as tabs as assumed field separators then most likely a CSV file
				LikelyCSV += 10;
			else							// if more tabs than commas then could be either BED or SAM
				{
				if(bSeenSAMHdrs)
					{
					LikelyBED += 5;
					LikelySAM += 10;			// SAM would be distinguished by it's header lines starting with '@"
					}
				else
					{
					LikelyBED += 20;
					LikelySAM += 5;
					}
				}
			}
		continue;
		}
	
	if(bSkipEOL)
		continue;

	if(!FieldCnt && !FldLen && (Chr == ' ' || Chr == '\t'))		// simply slough all leading whitespaces before intial field starts
		continue;

	// nested quotes are potentially a problem; currently quotes are simply sloughed
	if(Chr == '\'' || Chr == '"')
		{
		bInQuotes = !bInQuotes;
		continue;
		}


	if(!FieldCnt && !bInQuotes && Chr == '@' || Chr == '>')
		{
		if(Chr == '@')				// if SAM then header line(s) should be present and can be expected to start with "@HD", "@SQ", "@RG", "@PG", "@CO" 
			{		
			if(BuffIdx <  (BuffLen - 3))
				{
				if(((*pBuff == 'H' && pBuff[1] == 'D') ||
					(*pBuff == 'S' && pBuff[1] == 'Q') ||
					(*pBuff == 'R' && pBuff[1] == 'G') ||
					(*pBuff == 'P' && pBuff[1] == 'G') ||
					(*pBuff == 'C' && pBuff[1] == 'O')) &&
					(pBuff[2] == ' ' || pBuff[2] == '\t' ))
					{
					bSeenSAMHdrs = true;
					LikelyNonCSVSAMBED = 0;
					LikelySAM += 10000;
					bSkipEOL = true;
					continue;
					}
				else
					{
					if(!bSeenSAMHdrs)					// if no SAM headers parsed then could easily be a fastq...
						{
						LikelyNonCSVSAMBED += 50;
						bSkipEOL = true;
						continue;
						}
					}
				}
			}
		if(Chr == '>')									// if at start of line then could easily be fasta...
			LikelyNonCSVSAMBED += 50;
		}

	switch(Chr) {
		case ' ':			// simply slough spaces
			continue;

		case ',':			// if comma then likely is a csv, but could still be BED if in optional fields 9 (itemRgb) onwards 
			if(TabCnt < 8 && FieldCnt >= TabCnt)
				{
				FieldCnt += 1;	
				CommaCnt += 1;
				FldLen = 0;
				}
			break;

		case '\t':			// tabs are in BED and SAM as field separators, but could also be present in CSV as spacers
			if(CommaCnt < 3 && FieldCnt >= CommaCnt)
				{
				FieldCnt += 1;	
				TabCnt += 1;
				FldLen = 0;
				}
			break;

		default:			// any other char is assumed to be part of an actual field value
			FldLen += 1;
			break;
		}
	}

if(LikelyNonCSVSAMBED >= 250 || (LikelyCSV < 10 && LikelyBED < 10 && LikelySAM < 500))
	return(eCFTunknown);

if(LikelyCSV >= LikelyBED && LikelyCSV >= LikelySAM)
	return(eCFTCSV);
if(LikelyBED >= LikelySAM)
	return(eCFTBED);
return(eCFTSAM);	
}

bool						// true if requested number of bytes was written to pOutFile
CUtility::SafeWrite(void *pOutFile,		// either gzFile or file handle
					void *pBuff,		// pts to buffer containing data to write
					size_t Max2Write,	// write this many bytes from pBuff
					bool bgzOut)		// true if pOutFile pts to a gzFile, false if pOutFile is file handle
{
if(pOutFile == NULL || pBuff == NULL || Max2Write == 0)
	return(false);
if(bgzOut)
	return(SafeWrite_gz(*(gzFile *)pOutFile,pBuff,Max2Write));
return(RetryWrites(*(int *)pOutFile,pBuff,Max2Write));
}

bool	// true if requested number of bytes have been compressed to pgzFile
CUtility::SafeWrite_gz(gzFile pgzFile,void *pBuff,size_t Max2Write)
{
int Written;
int BuffLen;
uint8_t *pByte = (uint8_t *)pBuff;
Written = 0;
while(Written >= 0 && Max2Write > 0)
	{
	BuffLen = (int)min(Max2Write,(size_t)(cgzAllocOutBuffer * 9) / 10);			// limit writes to max of 90% cgzAllocOutBuffer bytes at a time
	Written = gzwrite(pgzFile,pBuff,BuffLen);
	if(Written <= 0)
		{
		printf("SafeWrite_gz: gzwrite %d bytes failed",BuffLen);
		break;
		}
	
	pByte += Written;
	BuffLen -= Written;
	Max2Write -= Written;
	Written = 0;
	}
return(Max2Write == 0 ? true : false);
}


bool				// true if requested number of bytes was written to hOutFile 
CUtility::RetryWrites(int hOutFile,void *pBuff,size_t Max2Write,
					  int MaxSecs)		// retry writes for at most this number of secs until complete Max2Write bytes written 
{
int Secs2Write;
int Written;
int BuffLen;
uint8_t *pByte = (uint8_t *)pBuff;
Written = 0;
Secs2Write = 0;
while(Written >= 0 && Max2Write > 0)
	{
	BuffLen = (int)min(Max2Write,(size_t)0x03fffffff);			// limit writes to max of 1GB at a time
	while(BuffLen) {
		Written = write(hOutFile, pByte, BuffLen);
		if(Written == -1 && !(errno == EINTR || errno == EAGAIN || errno == EWOULDBLOCK))
			{
			printf("RetryWrites: failed - %s",strerror(errno));
			return(false);
			}
		if(Written > 0)		
			{
			pByte += Written;
			BuffLen -= Written;
			Max2Write -= (size_t)Written;
			}

		if(Written <= 0 && Max2Write)
			{
			Secs2Write += 1;
			if(Secs2Write > MaxSecs)
				{
				printf("RetryWrites: failed - timed out");
				return(false);
				}
#ifdef _WIN32
			Sleep(1000);	// if write() was blocked or interrupted then no point in hammering the system, back off for 1 sec
#else
			sleep(1);
#endif
			}
		else
			Secs2Write = 0;
		Written = 0;
		}
	}
return(Max2Write == 0 ? true : false);
}

char *			// returns current and maximum resource limits ( getrlimit ), returns NULL if WIN32
CUtility::ReportResourceLimits(void)
{
#ifdef WIN32		// currently not bothering with windows resources
	return(NULL);
#else
int Ofs;
int Idx;
int PhysMemGB;
int NumberOfProcessors;
char szCurLimit[50];
char szHardLimit[50];
static char szResourceLimits[500];
struct rlimit RLimit;
Ofs = 0;
NumberOfProcessors = sysconf(_SC_NPROCESSORS_ONLN);
PhysMemGB = (int)(((((uint64_t)sysconf(_SC_PHYS_PAGES) * (uint64_t)sysconf(_SC_PAGESIZE)) + 0x01fffffff)) / ((uint64_t)0x03fffffff));

Ofs = sprintf(szResourceLimits,"cores: %d physmem: %d (GB) ", NumberOfProcessors, PhysMemGB);

for(Idx = 0; Idx < 5; Idx++)
	{
	switch(Idx) {
		case 0:        // address space
			getrlimit(RLIMIT_AS,&RLimit);
			Ofs += sprintf(&szResourceLimits[Ofs],"virtual mem:");
			break;
		case 1:        // max locked memory
			getrlimit(RLIMIT_MEMLOCK,&RLimit);
			Ofs += sprintf(&szResourceLimits[Ofs],"locked mem:");
			break;
		case 2:        // data seg size
			getrlimit(RLIMIT_DATA,&RLimit);
			Ofs += sprintf(&szResourceLimits[Ofs],"data seg size:");
			break;
		case 3:        // threads
			getrlimit(RLIMIT_NPROC,&RLimit);
			Ofs += sprintf(&szResourceLimits[Ofs],"threads:");
			break;
		case 4:        // stack size
			getrlimit(RLIMIT_STACK,&RLimit);
			Ofs += sprintf(&szResourceLimits[Ofs],"stack size:");
			break;			
		}
	if(RLimit.rlim_cur == RLIM_INFINITY)
		strcpy(szCurLimit,"Unlimited");
	else
		sprintf(szCurLimit,"%lu",RLimit.rlim_cur);
	if(RLimit.rlim_max == RLIM_INFINITY)
		strcpy(szHardLimit,"Unlimited");
	else
		sprintf(szHardLimit,"%lu",RLimit.rlim_max);
	Ofs += sprintf(&szResourceLimits[Ofs]," %s (cur) %s (max) ", szCurLimit, szHardLimit);
	}
return(szResourceLimits);
#endif
}

// GetNumSubseqs
// Returns number of subsequences - no filtering - common w/o InDels
// between the sequences in pProcParams->pSeq
// Note: see GetNumFilteredSubseqs() for filtered processed equivalent
int
CUtility::GetNumSubseqs(int AlignLen,		// alignment length incl InDels
					   int NumSeqs2Align,	// number of sequences
						etSeqBase *pSeqs[1]) // pts to array of ptrs to the sequences to process
{
int SeqIdx;
int Cnt;
int VIdx;
int NumSeqs;
bool bInDel = false;
int CurSeqLen = 0;
etSeqBase *pSeq;
etSeqBase SeqBase;
NumSeqs = 0;
for(SeqIdx = Cnt = 0; Cnt < AlignLen; Cnt++, SeqIdx++)
	{
	for(VIdx = 0; VIdx < NumSeqs2Align; VIdx++)
		{
		pSeq = pSeqs[VIdx];
		SeqBase = pSeq[SeqIdx] & ~cRptMskFlg;
		if(SeqBase == eBaseInDel || SeqBase == eBaseN)
			bInDel = true;
		}

	if(bInDel)
		{
		if(CurSeqLen > 0)		// if previously was processing a subsequence
			NumSeqs++;
		CurSeqLen = 0;
		bInDel = false;
		continue;
		}

	// no InDels in any of the required species, must be part of an aligned subsequence
	CurSeqLen++;
	}

if(CurSeqLen > 0)
	NumSeqs++;
return(NumSeqs);
}


// GetFilteredSubSeq
// Returns first subsequence meeting filtering requirements
// Filtering:
// Left and right flanking bases of subsequence are removed untill base which is identical over all species in alignment
// If subsequence length is less than ReqMinLen then subsequence is sloughed
// If overall identity within subsequence is less than ReqMinIdent then subsequence is sloughed
int								// returned subsequence length
CUtility::GetFilteredSubSeq(int AlignLen, // remaining alignment length incl InDels
		  int StartIdx,			// initial starting index into pRef[]/pRel[]
		  int *pFirstIdx,		// returned index into pRef[]/pRel[] of first base in returned subsequence
		  int *pLastIdx,		// returned index into pRef[]/pRel[] of last base in returned subsequence
		  etSeqBase *pRef,		// reference sequence
		  etSeqBase *pRel,		// relative sequence
		  int ReqMinLen,		// subsequence must be of at least this length
		  int ReqMinIdent,		// and have at least this overall identity (0..100%)
		  bool bStartEndIdent)  // and start/end on base that is identical
{
etSeqBase RefBase;
etSeqBase RelBase;
int Idx;
int SeqLen;
int SeqIdent;
int NumIdent;
int FirstMatchIdx = -1;
int LastMatchIdx = -1;
int SubSeqCnt = 0;
if(StartIdx < 0 || 
   AlignLen <= 0 || ReqMinLen < 1 || AlignLen < ReqMinLen || 
   pRef == NULL || pRel == NULL ||
   pFirstIdx == NULL || pLastIdx == NULL ||
   ReqMinIdent < 0 || ReqMinIdent > 100)
	return(0);

pRef = &pRef[StartIdx];
pRel = &pRel[StartIdx];
for(Idx = 0; Idx < AlignLen; Idx++)
	{
	RefBase = *pRef++ & ~cRptMskFlg;
	RelBase = *pRel++ & ~cRptMskFlg;

	if(RefBase > eBaseT || RelBase > eBaseT)	// treat any base not a,c,g,t as an InDel
		{
		if(FirstMatchIdx >= 0)		// if at least one identical base then it might be a subsequence...
			{
			SeqLen = LastMatchIdx - FirstMatchIdx + 1;
			if(SeqLen >= ReqMinLen)	// if subsequence of minimum required length
				{
				if(ReqMinIdent)			// does it need to have a minimum identity?
					SeqIdent = (NumIdent * 100) / SeqLen;
				if(!ReqMinIdent || SeqIdent >= ReqMinIdent)
					{
					// meets filtering requirements!!!!
					*pFirstIdx = StartIdx + FirstMatchIdx;
					*pLastIdx = StartIdx + LastMatchIdx;
					return(SeqLen);
					}
				}
			LastMatchIdx = FirstMatchIdx = -1;
			}
		continue;
		}

	// not an InDel, must be part of an aligned subsequence
	if(RefBase == RelBase)			// if identical...
		{
		if(FirstMatchIdx == -1)		// if first then note where identity occured
			{
			FirstMatchIdx = Idx;
			NumIdent = 0;
			}
		NumIdent++;
		LastMatchIdx = Idx;			// note where last identity occured
		continue;
		}

	if(!bStartEndIdent)
		{
		if(FirstMatchIdx == -1)		// if first base in subsequence then note where occured
			{
			FirstMatchIdx = Idx;
			NumIdent = 0;
			}
		LastMatchIdx = Idx;			// note where last identity occured
		}
	}

if(FirstMatchIdx >= 0)		// if at least one identical base then it might be a subsequence...
	{
	SeqLen = LastMatchIdx - FirstMatchIdx + 1;
	if(SeqLen >= ReqMinLen)	// if subsequence of minimum required length
		{
		if(ReqMinIdent)			// does it need to have a minimum identity?
			SeqIdent = (NumIdent * 100) / SeqLen;
		if(!ReqMinIdent || SeqIdent >= ReqMinIdent)
			{
			// meets filtering requirements!!!!
			*pFirstIdx = StartIdx + FirstMatchIdx;
			*pLastIdx = StartIdx + LastMatchIdx;
			return(SeqLen);
			}
		}
	}
*pFirstIdx = -1;
*pLastIdx = -1;
return(0);
}



// GetNumFilteredSubseqs
// Returns number of subsequences in alignment block which would meet filtering requirements
// as processed by GetFilteredSubSeq()
int
CUtility::GetNumFilteredSubseqs(int AlignLen,	// alignment length incl InDels
		   etSeqBase *pRef,		// reference sequence
		   etSeqBase *pRel,		// relative sequence
		  int ReqMinLen,		// subsequence must be of at least this length
		  int ReqMinIdent,		// and have at least this overall identity (0..100%)
		  bool bStartEndIdent)  // and start/end on base that is identical
{
etSeqBase RefBase;
etSeqBase RelBase;
int Idx;
int SeqLen;
int SeqIdent;
int NumIdent;
int FirstMatchIdx = -1;
int LastMatchIdx = -1;
int SubSeqCnt = 0;

if(AlignLen <= 0 || ReqMinLen < 1 || AlignLen < ReqMinLen || 
   pRef == NULL || pRel == NULL || 
   ReqMinIdent < 0 || ReqMinIdent > 100)
	return(0);

for(Idx = 0; Idx < AlignLen; Idx++)
	{
	RefBase = *pRef++ & ~cRptMskFlg;
	RelBase = *pRel++ & ~cRptMskFlg;

	if(RefBase > eBaseT || RelBase > eBaseT)	// treat any base not a,c,g,t as an InDel
		{
		if(FirstMatchIdx >= 0)		// if at least one identical base then it might be a subsequence...
			{
			SeqLen = LastMatchIdx - FirstMatchIdx + 1;
			if(SeqLen >= ReqMinLen)	// if subsequence of minimum required length
				{
				if(ReqMinIdent)			// does it need to have a minimum identity?
					SeqIdent = (NumIdent * 100) / SeqLen;
				if(!ReqMinIdent || SeqIdent >= ReqMinIdent)
					SubSeqCnt += 1;
				}
			LastMatchIdx = FirstMatchIdx = -1;
			}
		continue;
		}

	// not an InDel, must be part of an aligned subsequence
	if(RefBase == RelBase)			// if identical...
		{
		if(FirstMatchIdx == -1)		// if first then note where identity occured
			{
			FirstMatchIdx = Idx;
			NumIdent = 0;
			}
		NumIdent++;
		LastMatchIdx = Idx;			// note where last identity occured
		continue;
		}

	if(!bStartEndIdent)
		{
		if(FirstMatchIdx == -1)		// if first base in subsequence then note where occured
			{
			FirstMatchIdx = Idx;
			NumIdent = 0;
			}
		LastMatchIdx = Idx;			// note where last identity occured
		}
	}

if(FirstMatchIdx >= 0)		// if at least one identical base then it might be a subsequence...
	{
	SeqLen = LastMatchIdx - FirstMatchIdx + 1;
	if(SeqLen >= ReqMinLen)	// if subsequence of minimum required length
		{
		if(ReqMinIdent)			// does it need to have a minimum identity?
			SeqIdent = (NumIdent * 100) / SeqLen;
		if(!ReqMinIdent || SeqIdent >= ReqMinIdent)
			SubSeqCnt += 1;
		}
	}
return(SubSeqCnt);
}


// TrimQuotes
// Removes (inplace) any leading and trailing quotes from pszTxt
int
CUtility::TrimQuotes(char *pszTxt)
{
char *pDst = pszTxt;
char Chr;
int Len = 0;
while((Chr = *pszTxt++))	
	{
	if(Chr == '"' || Chr == '\'')
		continue;
	*pDst++ = Chr;
	Len++;
	}
*pDst = Chr;
return(Len);
}

// TrimWhitespc
// Returns ptr to 1st non-whitespace char in buffer ptd at by pszText
// Will also have trimmed down returned text by replacing trailing whitespace with '\0' terminator 
// Will return NULL if pszText is empty after trimming
char *
CUtility::TrimWhitespc(char *pszText)
{
char Chr;
char *pChr;
// trim any leading/trailing whitespace
while(Chr = *pszText++)
	{
	if(!isspace(Chr))
		break;
	}
pszText-=1;
if(Chr == '\0')		// all whitespace?
	return(NULL);

pChr = pszText;
while(Chr = *pChr++);	// locate terminator '\0'
pChr -= 2;				// backup to last chr
while(Chr = *pChr--)	// backup as long as non whitespace
	if(!isspace(Chr))
		break;
pChr[2] = '\0';
return(pszText);
}

// TrimWhitespcExtd
// Removes any leading/trailing whitespace and
// copies trimmed text back into original text buffer (pszText)
char *
CUtility::TrimWhitespcExtd(char *pszText)
{
char *pszTrimed;
char *pszMov;
if(pszText == NULL || pszText[0] == '\0')
	return(pszText);
if((pszTrimed = TrimWhitespc(pszText)) == NULL)
	{
	pszText[0] = '\0';
	return(pszText);
	}
pszMov = pszText; 
while (*pszMov++ = *pszTrimed++);
return(pszText);
}

// TrimQuotedWhitespcExtd
// Removes leading/trailing quotes and then any leading/trailing whitespace and
// copies trimmed text back into original text buffer (pszText)
char *
CUtility::TrimQuotedWhitespcExtd(char *pszText)
{
if(pszText == NULL || pszText[0] == '\0')
	return(pszText);
if(!TrimQuotes(pszText))
	return(pszText);
return(TrimWhitespcExtd(pszText));
}

// ReduceWhitespace
// Removes leading/trailing whitespace and replaces any internal whitespace with a single space
char *
CUtility::ReduceWhitespace(char *pszText)
{
char *pDst;
char *pSrc;
char Chr;
if(pszText == NULL || pszText[0] == '\0')
	return(pszText);
TrimWhitespcExtd(pszText);
if(pszText[0] == '\0')
	return(pszText);
pSrc = pszText;
pDst = pszText;
while(Chr = *pSrc++)
	{
	if(!isspace(Chr))
		{
		*pDst++ = Chr;
		continue;
		}
	if(pDst[-1] != ' ')
		*pDst++ = ' ';
	}
*pDst = '\0';
return(pszText);
}


// cleantext
// inplace cleans potentially troublesome chars from text by replacing these chars with a single space char
char *
CUtility::CleanText(char* pText)
{
char Chr;
char *pTxt;
if((pTxt = pText) == NULL)
	return(NULL);

while((Chr = *pTxt++) != '\0')
	{
	if(Chr < 0x20 || Chr == '%' || Chr == '\'' || Chr == '"')
		pTxt[-1] = ' ';
	}
return(pText);
}


const int cMaxNumOptions = 60100;		// allow at most this many options in option file, increased from original 1024 because snpmarkers needs to handle up to 20000 (max of 3 params per isolate) different isolates 
const int cMaxOptionLen  = 1024;		// any single option can be of this length
const int cMaxOptionLineLen = 8196;		// allow for option lines to be upto this length, one line can contain multiple options
const int cMaxTotOptionsLen = (cMaxNumOptions * 101);	// allow for options to total at most this many chars - assumes options average 100 chars incl '\0' separators

// iterate over options checking if some options are contained in an external option file]
// name of external option file is preceded by '@' i.e. '@filename' specifies to initially read parameters from this file
int												// returned cnt of args ptd to by retargv
CUtility::arg_parsefromfile(int argc,			// cnt of args ptd to by *argv[]
							char *argv[],		// pts to array of args
							char **retargv[])	// returned array of args
{
char *pszOptn;
char *pszOptnFile;
int NumOptns;
char szOptnLineBuff[cMaxOptionLineLen];
char szOption[cMaxOptionLen];				// to hold a single option as parsed from szOptnLineBuff
static char szOptnsBuff[cMaxTotOptionsLen];
static char *pszOptns[cMaxNumOptions];		// handles at most cMaxNumOptions argv options
char Chr;
char *pChr;
char *pOpt;
bool bInQuotes;
bool bInParam;
bool bEOL;

int NxtOptnsBuffIdx;
FILE *pOptStream;

if(retargv == NULL || argc > cMaxNumOptions)
	return(-1);
*retargv = NULL;

NumOptns = 0;
NxtOptnsBuffIdx = 0;
for(int OptnCnt = 0; OptnCnt < argc; OptnCnt++)
	{
	if(NumOptns == cMaxNumOptions)
		break;
	if(*argv[OptnCnt] == '@')
		{
		pszOptnFile = TrimWhitespc(argv[OptnCnt]+1);
		if(pszOptnFile == NULL)
			{
			printf("No options file specified following '@' switch");
			return(-1);
			}

		// open options file
		if((pOptStream = fopen(pszOptnFile,"r"))==NULL)
			{
			printf("Unable to open options file '%s'\nError: %s",pszOptnFile,strerror(errno));
			return(-1);
			}

		// parse each line into one or more options using whitespace as the option separator except where option has been quoted
		while(fgets(szOptnLineBuff,sizeof(szOptnLineBuff),pOptStream)!= NULL)
			{
			// trim any leading/trailing whitespace
			pszOptn = TrimWhitespc(szOptnLineBuff);
			if(pszOptn == NULL || pszOptn[0] == '\0' ||		// slough empty lines or lines with 1st non-whitespace being '#' or ';' or "//" as comment lines
				pszOptn[0] == '#' || pszOptn[0] == ';' || (pszOptn[0] == '/' && pszOptn[1] == '/'))
				continue;

			// parse out each option delimited by whitespace
			pChr = pszOptn;
			pOpt = szOption;
			bInQuotes = false;
			bInParam = false;
			bEOL = false;
			while(!bEOL) {
				Chr = *pChr++;
				if(Chr == 0x16)						// special fudge, on some systems '-' is represented as 0x16
					Chr = '-';
				*pOpt++ = Chr;
				switch(Chr) {
					case 0x00:						// EOL
						bEOL = true;
						break;

					case '"':						// quotes
					case '\'':
						bInQuotes = !bInQuotes;		// toggle inquote indicator
						bInParam = true;			// treat as still in parameter until whitespace or EOL
						continue;

					case ' ':
					case '\t':
						if(bInQuotes)				// allow whitespace inside quotes
							continue;
						if(!bInParam)				// if between parameters then slough whitespace
							{
							pOpt -= 1;
							continue;
							}
						pOpt[-1] = '\0';			// terminate current option
						break;

					default:						// any other character is a parameter chr
						bInParam = true;
						continue;
					}
			
				pszOptns[NumOptns] = &szOptnsBuff[NxtOptnsBuffIdx];
				strcpy(pszOptns[NumOptns++],szOption);
				NxtOptnsBuffIdx += (int)strlen(szOption) + 1;
				pOpt = szOption;
				bInQuotes = false;
				bInParam = false;
				if(NumOptns == cMaxNumOptions)
					break;
				}
			if(NumOptns == cMaxNumOptions)
				break;
			}

		fclose(pOptStream);
		}
	else
		{
		pszOptns[NumOptns] = &szOptnsBuff[NxtOptnsBuffIdx];
		strcpy(pszOptns[NumOptns++],argv[OptnCnt]);
		NxtOptnsBuffIdx += (int)strlen(argv[OptnCnt]) + 1;
		if(NumOptns == cMaxNumOptions)
			break;
		}
	}
*retargv = pszOptns;
return(NumOptns);
}

// ChkTargDepend
// Check list of source files against a target file
// If target file does not exist, or any source file is newer than target then returns 1
// If target file is newer or same age as all the source files then returns 0
// If any source file is specified but does not exist then returns -1 and copies file name into pszInaccessible
// Use NULL or empty source file name to terminate list of source files
int
CUtility::ChkTargDepend(char *pszInaccessible,	// where to copy name of the first detected inaccessible source file (if any)
						int MaxLen,				// max length to copy
						char *pszTarg,
						 ...)					// var number of source files, use NULL or empty (pszSource[0] == '\0') to terminate 
{
int StatRslt;
char *pszSrc;
struct stat TargStat;
struct stat SrcStat;
va_list argptr;
va_start(argptr, pszTarg);

if(pszInaccessible != NULL && MaxLen > 0)
	*pszInaccessible = '\0';

StatRslt = stat(pszTarg,&TargStat);

if(StatRslt < 0)	// < 0 if unable to access target file
	return(1);

while(((pszSrc = va_arg( argptr, char *)) != NULL) && pszSrc[0] != '\0')
   {
   if(*pszSrc == '\0')						  // unspecified independent src files are simply sloughed
	   continue;
   if((StatRslt = stat(pszSrc,&SrcStat)) < 0) // expect independent src file to at least be accessible
		{
		if(pszInaccessible != NULL && MaxLen > 1)
			{
			strncpy(pszInaccessible,pszSrc,MaxLen);
			pszInaccessible[MaxLen-1] = '\0';
			}
		return(-1);							  // treat inaccessible independents as errors
		}  
	if(TargStat.st_mtime < SrcStat.st_mtime)   // if target older than dependent source then return 1 
		return(1);
   }
va_end( argptr );              /* Reset variable arguments.      */
return(0);	// target exists and is no older than any independent source file
}

// Chk2TargDepend
// Check list of source files against two target files
// If target file does not exist, or any source file is newer than target then returns 1
// If target file is newer or same age as all the source files then returns 0
// If any source file is specified but does not exist then returns -1 and copies file name into pszInaccessible
// Use NULL or empty source file name to terminate list of source files
int
CUtility::Chk2TargDepend(char *pszInaccessible,	// where to copy name of the first detected inaccessible source file (if any)
						int MaxLen,				// max length to copy
						char *pszTarg1,char *pszTarg2,
						 ...)					// var number of source files, use NULL or empty (pszSource[0] == '\0') to terminate
{
int StatRslt;
char *pszSrc;
struct stat TargStat1;
struct stat TargStat2;
struct stat SrcStat;
va_list argptr;
va_start(argptr, pszTarg2);

if(pszInaccessible != NULL && MaxLen > 0)
	*pszInaccessible = '\0';

StatRslt = stat(pszTarg1,&TargStat1);
if(StatRslt < 0)	// < 0 if unable to access target 1 file
	return(1);

if(pszTarg2 != NULL && *pszTarg2)
	{
	StatRslt = stat(pszTarg2,&TargStat2);
	if(StatRslt < 0)	// < 0 if unable to access target 2 file
		return(1);
	}

while(((pszSrc = va_arg( argptr, char *)) != NULL) && pszSrc[0] != '\0')
   {
   if(*pszSrc == '\0')						  // unspecified independent src files are simply sloughed
	   continue;
   if((StatRslt = stat(pszSrc,&SrcStat)) < 0) // expect independent src file to at least be accessible
		{
		if(pszInaccessible != NULL && MaxLen > 1)
			{
			strncpy(pszInaccessible,pszSrc,MaxLen);
			pszInaccessible[MaxLen-1] = '\0';
			}
		return(-1);							  // treat inaccessible independents as errors
		}
   if(TargStat1.st_mtime < SrcStat.st_mtime)   // if target 1 older than dependent source then return 1 
		return(1);
	if(pszTarg2 != NULL && *pszTarg2 && TargStat2.st_mtime < SrcStat.st_mtime)   // if target 2 older than dependent source then return 1 
		return(1);
   }
va_end( argptr );              /* Reset variable arguments.      */
return(0);	// target exists and is no older than any independent source file
}
	
// splits the path at the last '\' (windows) or '/' (linux), returning the prefix as the directory and the suffix as the filename
void 
CUtility::splitpath(char *pszFullPath, char *pszDir,
				 char *pFname)
{
int PathLen = (int)strlen(pszFullPath);
char *p, *end;
if(pszDir)
		*pszDir = 0;
if(pFname)
	*pFname = 0;

	 /* look for end of directory part */
end = NULL;
for (p = pszFullPath; *p; p++) 
	if (*p == '/' || *p == '\\') 
		end = p + 1;

if(end)  /* got a directory */
	{
	if (pszDir)
		{
		 memmove( pszDir, pszFullPath, (end - pszFullPath) * sizeof(char));
		pszDir[end - pszFullPath] = 0;
		}
	pszFullPath = end;
	}
if(pFname != NULL && *pszFullPath)
	{
	do {
		*pFname++ = *pszFullPath;
		}
	while(*pszFullPath++);
	}
}

// utility function for appending a filename suffix onto a base file name starting from the last occurrence of a marker char in the base file name
int				// returned length of suffixed file name
CUtility::AppendFileNameSuffix(char *pszRetSuffixedName,	// buffer allocated by caller to hold suffixed name
					 char *pszBaseFileName,		// base file name prefix
					 char *pszSuffix,			// suffix to append onto base file name at last separator marker char in base file name
					 char MarkerChar)			// separator marker char in base name at which to append suffix from 
{
if(pszRetSuffixedName == NULL)
	return(0);
pszRetSuffixedName[0] = '\0';
if(pszBaseFileName == NULL || pszBaseFileName[0] == '\0')
	{
	if(pszSuffix != NULL && pszSuffix[0] != '\0')
		{
		strcpy(pszRetSuffixedName, pszSuffix);
		return((int)strlen(pszRetSuffixedName));
		}
	return(0);
	}

int Idx;

strcpy(pszRetSuffixedName, pszBaseFileName);
int BaseLen = (int)strlen(pszRetSuffixedName);
char *pChr = &pszRetSuffixedName[BaseLen-1];
for(Idx = BaseLen; Idx > 0; Idx-=1, pChr-=1)
	{
	if (*pChr == MarkerChar)
		{
		*pChr = '\0';
		break;
		}
	if(*pChr == '\\' || *pChr == '/')	// don't suffix onto any path components
		break;
	}
if(pszSuffix != NULL && pszSuffix[0] != '\0')
	strcat(pszRetSuffixedName, pszSuffix);
return((int)strlen(pszRetSuffixedName));
}

void 
CUtility::SleepMillisecs(uint32_t milliseconds) // cross-platform sleep function
{
#ifdef WIN32
	Sleep(milliseconds);
#elif _POSIX_C_SOURCE >= 199309L
	struct timespec ts;
	ts.tv_sec = milliseconds / 1000;
	ts.tv_nsec = (milliseconds % 1000) * 1000000;
	nanosleep(&ts, NULL);
#else
	usleep(milliseconds * 1000);
#endif
}

