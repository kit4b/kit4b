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
// This version can handle DNA in both colour and base space
//
//
#include "stdafx.h"
#ifdef _WIN32
#include <process.h>
#include "./commhdrs.h"
#else
#include <sys/mman.h>
#include <pthread.h>
#include "./commhdrs.h"
#endif

#include "./SfxArray.h"

static int gMaxBaseCmpLen = cMaxReadLen;   // used to limit the number of bases being compared for length in QSortSeqCmp32() and QSortSeqCmp40()
static int QSortSeqCmp32(const void *p1,const void *p2);
static int QSortSeqCmp40(const void *p1,const void *p2);
static int QSortEntryNames(const void *p1,const void *p2);

static uint8_t *gpSfxArray = NULL;
static etSeqBase *gpSeq = NULL;

// Suffix array elements can be sized as either 4 or 5 bytes dependent on the total length of concatenated sequences
// If total length is less than 4G then can use 4 byte elements, if longer then will use 5 byte elements
inline 
int64_t SfxOfsToLoci(int SfxElSize,			// size in bytes of suffix element - expected to be either 4 or 5
				void *pSfx,					// pts to 1st element of suffix array
				int64_t Ofs)					// offset to suffix element
{
uint64_t Loci;
uint8_t *pSfxEls = (uint8_t *)pSfx;
Ofs *= SfxElSize;
Loci = (uint64_t)*(uint32_t *)&pSfxEls[Ofs];
if(SfxElSize == 5)
	Loci |= ((uint64_t)pSfxEls[Ofs+4]) << 32;
return(Loci);
}


// constructor
CSfxArray::CSfxArray(void)
{
m_pEntriesBlock = NULL;
m_pSfxBlock = NULL;
m_pBisulfateBases = NULL;
m_pOccKMerClas = NULL;
m_pCoreKMers = NULL;
m_pKMerCntDist = NULL;
m_hFile = -1;
m_bThreadActive = false;
m_AllocSfxBlockMem = 0;
m_AllocEntriesBlockMem = 0;
m_AllocSfxBlockMem = 0;
m_AllocBisulfiteMem = 0;
m_EstSfxEls = 0;
m_CoreKMerLen = 0;
m_NumCoreKMers = 0;
m_AllocdCoreKMersMem = 0;
m_MaxIter = cDfltMaxIter;
m_bInMemSfx = false;
m_MaxQSortThreads = cDfltSortThreads;
m_MTqsort.SetMaxThreads(m_MaxQSortThreads);
m_MaxSfxBlockEls = cMaxAllowConcatSeqLen;
m_CASSeqFlags = 0;
m_MaxKMerCnts = 0;
m_pKMerCntDist = NULL;
m_CASSeqFlags = 0;
m_AllocEntriesBlockMem = 0;
m_AllocSfxBlockMem = 0;
m_AllocBisulfiteMem = 0;
m_bTermThread = false;
m_ReqBlockRslt = eBSFSuccess;
m_ReqBlockID = 0;
m_szFile[0] = '\0';
m_bHdrDirty = false;
m_bCreate = false;
m_bBisulfite = false;
m_bColorspace = false;
m_MaxIter = cDfltMaxIter;
m_bV3File = false;
m_bInMemSfx = false;
m_MaxMMExploreInDel = cMaxMMExploreInDel;
m_MaxInDelLen = cMaxMicroInDelLen;
m_MinInDelSeqLen = cMinInDelSeqLen;
m_MaxSfxBlockEls = cMaxAllowConcatSeqLen;
m_MaxKMerOccs = 0;
m_AllocOccKMerClasMem = 0;
m_OccKMerLen = 0;
m_CoreKMerLen = 0;
m_NumCoreKMers = 0;
m_AllocdCoreKMersMem = 0;


gMaxBaseCmpLen = cMaxReadLen;

#ifdef _WIN32
m_threadID = 0;
m_JobMutex = NULL;
m_JobReqEvent = NULL;
m_JobAckEvent = NULL;
#endif
}

//destructor
CSfxArray::~CSfxArray(void)
{
if(m_bThreadActive)	// if processing reads against sfx blocks then need to terminate the readahead background thread
	{
	m_bTermThread = true;
	m_ReqBlockRslt = (teBSFrsltCodes)1;
#ifdef _WIN32
	SetEvent(m_JobReqEvent);
	WaitForSingleObject(m_threadHandle,cSigTermWaitSecs * 1000);
	CloseHandle(m_threadHandle);
	CloseHandle(m_JobMutex);
	CloseHandle(m_JobReqEvent);
	CloseHandle(m_JobAckEvent);
#else
	pthread_cond_signal(&m_JobReqEvent);
	pthread_join(m_threadID,NULL);
	pthread_mutex_destroy(&m_JobMutex);
	pthread_cond_destroy(&m_JobReqEvent);
	pthread_cond_destroy(&m_JobAckEvent);
#endif
	}

if(m_pSfxBlock != NULL)
	{
#ifdef _WIN32
	free(m_pSfxBlock);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pSfxBlock != MAP_FAILED)
		munmap(m_pSfxBlock,m_AllocSfxBlockMem);
#endif
	}

if(m_pEntriesBlock != NULL)
	{
#ifdef _WIN32
	free(m_pEntriesBlock);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pEntriesBlock != MAP_FAILED)
		munmap(m_pEntriesBlock,m_AllocEntriesBlockMem);
#endif
	}

if(m_pBisulfateBases != NULL)
	{
#ifdef _WIN32
	free(m_pBisulfateBases);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pBisulfateBases != MAP_FAILED)
		munmap(m_pBisulfateBases,m_AllocBisulfiteMem);
#endif
	}

if(m_pOccKMerClas != NULL)
	{
#ifdef _WIN32
	free(m_pOccKMerClas);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pOccKMerClas != MAP_FAILED)
		munmap(m_pOccKMerClas,m_AllocOccKMerClasMem);
#endif
	}

if(m_pCoreKMers != NULL)
	{
#ifdef _WIN32
	free(m_pCoreKMers);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (m_pCoreKMers != MAP_FAILED)
		munmap(m_pCoreKMers, m_AllocdCoreKMersMem);
#endif
	}

if(m_hFile != -1)
	close(m_hFile);

}

int
CSfxArray::Reset(bool bFlush)			// reset state back to that immediately following instantiation
{
int Rslt;

if(m_bThreadActive)	// if still processing reads against sfx blocks then need to terminate the readahead background thread
	{
	m_bTermThread = true;
	m_ReqBlockRslt = (teBSFrsltCodes)1;
#ifdef _WIN32
	SetEvent(m_JobReqEvent);
	WaitForSingleObject(m_threadHandle,cSigTermWaitSecs * 1000);
	CloseHandle(m_threadHandle);
	m_threadHandle = NULL;
	CloseHandle(m_JobMutex);
	m_JobMutex = NULL;
	CloseHandle(m_JobReqEvent);
	m_JobReqEvent = NULL;
	CloseHandle(m_JobAckEvent);
	m_JobAckEvent = NULL;
#else
	pthread_cond_signal(&m_JobReqEvent);
	pthread_join(m_threadID,NULL);
	pthread_mutex_destroy(&m_JobMutex);
	pthread_cond_destroy(&m_JobReqEvent);
	pthread_cond_destroy(&m_JobAckEvent);
#endif
	m_bThreadActive = false;
	}

if(m_hFile != -1)
	{
	if(bFlush)
		Rslt = Flush2Disk();
	close(m_hFile);
	m_hFile = -1;
	}

memset(&m_SfxHeader,0,sizeof(m_SfxHeader));

if(m_pSfxBlock != NULL)
	{
#ifdef _WIN32
	free(m_pSfxBlock);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pSfxBlock != MAP_FAILED)
		munmap(m_pSfxBlock,m_AllocSfxBlockMem);
#endif
	m_pSfxBlock = NULL;
	}

if(m_pEntriesBlock != NULL)
	{
#ifdef _WIN32
	free(m_pEntriesBlock);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pEntriesBlock != MAP_FAILED)
		munmap(m_pEntriesBlock,m_AllocEntriesBlockMem);
#endif
	m_pEntriesBlock = NULL;
	}

if(m_pBisulfateBases != NULL)
	{
#ifdef _WIN32
	free(m_pBisulfateBases);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pBisulfateBases != MAP_FAILED)
		munmap(m_pBisulfateBases,m_AllocBisulfiteMem);
#endif
	m_pBisulfateBases = NULL;
	}

if(m_pOccKMerClas != NULL)
	{
#ifdef _WIN32
	free(m_pOccKMerClas);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pOccKMerClas != MAP_FAILED)
		munmap(m_pOccKMerClas,m_AllocOccKMerClasMem);
#endif
	m_pOccKMerClas = NULL;
	}

if (m_pCoreKMers != NULL)
	{
#ifdef _WIN32
	free(m_pCoreKMers);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (m_pCoreKMers != MAP_FAILED)
		munmap(m_pCoreKMers, m_AllocdCoreKMersMem);
#endif
	m_pCoreKMers = NULL;
	}

m_MaxKMerCnts = 0;
m_pKMerCntDist = NULL;
m_CASSeqFlags = 0;
m_AllocEntriesBlockMem = 0;
m_AllocSfxBlockMem = 0;
m_AllocBisulfiteMem = 0;
m_bTermThread = false;
m_ReqBlockRslt = eBSFSuccess;
m_ReqBlockID = 0;
m_szFile[0] = '\0';
m_bHdrDirty = false;
m_bCreate = false;
m_bBisulfite = false;
m_bColorspace = false;
m_MaxIter = cDfltMaxIter;
m_bV3File = false;
m_bInMemSfx = false;
m_MaxMMExploreInDel = cMaxMMExploreInDel;
m_MaxInDelLen = cMaxMicroInDelLen;
m_MinInDelSeqLen = cMinInDelSeqLen;
m_MaxSfxBlockEls = cMaxAllowConcatSeqLen;
m_MaxKMerOccs = 0;
m_AllocOccKMerClasMem = 0;
m_OccKMerLen = 0;
m_CoreKMerLen = 0;
m_NumCoreKMers = 0;
m_AllocdCoreKMersMem = 0;
return(eBSFSuccess);
}


void
CSfxArray::SetMaxQSortThreads(int MaxThreads)			// sets maximum number of threads to use in multithreaded qsorts
{
m_MaxQSortThreads = MaxThreads;
m_MTqsort.SetMaxThreads(MaxThreads);
}

int						// returns the previously utilised MaxBaseCmpLen
CSfxArray::SetMaxBaseCmpLen(int MaxBaseCmpLen)		// sets maximum number of bases which need to be compared for equality in multithreaded qsorts, will be clamped to be in range 10..cMaxReadLen
{
int PrevMaxBaseCmpLen = gMaxBaseCmpLen;
if(MaxBaseCmpLen > cMaxReadLen)
	MaxBaseCmpLen = cMaxReadLen;
else
	if(MaxBaseCmpLen < 10)
		MaxBaseCmpLen = 10;
gMaxBaseCmpLen = MaxBaseCmpLen;
return(PrevMaxBaseCmpLen);
}

void
CSfxArray::SetInitalSfxAllocEls(int64_t EstSfxEls)		// user estimate of the initial number of suffix elements to allocate
{
m_EstSfxEls = EstSfxEls;
}


// InitHdr
// Initialise header with default values
void
CSfxArray::InitHdr(void)
{
memset(&m_SfxHeader,0,sizeof(tsSfxHeaderV3));
m_SfxHeader.Magic[0] = 's';
m_SfxHeader.Magic[1] = 'f';
m_SfxHeader.Magic[2] = 'x';
m_SfxHeader.Magic[3] = '5';
m_SfxHeader.Version = cSFXVersion;	        // file structure version
m_SfxHeader.FileLen = sizeof(tsSfxHeaderV3);	// current file length (nxt write psn)
m_SfxHeader.szDatasetName[0] = '\0';
m_SfxHeader.szDescription[0] = '\0';
m_SfxHeader.szTitle[0] = '\0';
m_SfxHeader.Attributes = 0;
m_bHdrDirty = true;
}

// Hdr2Disk
// Writes file header out to disk after any appropriate endian processing
teBSFrsltCodes
CSfxArray::Hdr2Disk(void)
{
tsSfxHeaderV3 *pHdr;
int WrtLen;
if(!m_bHdrDirty)
	return(eBSFSuccess);
WrtLen = sizeof(tsSfxHeaderV3);
if(m_bBisulfite)
	m_SfxHeader.Attributes |= 0x01;
else
	m_SfxHeader.Attributes &= ~0x01;

if(m_bColorspace)
	m_SfxHeader.Attributes |= 0x02;
else
	m_SfxHeader.Attributes &= ~0x02;

pHdr = &m_SfxHeader;
if(_lseeki64(m_hFile,0,SEEK_SET) ||
		write(m_hFile,pHdr,WrtLen)!=WrtLen)
	{
	AddErrMsg("CSfxArray::Flush2Disk","Unable to write file header to disk on file %s - error %s",m_szFile,strerror(errno));
	Reset(false);
	return(eBSFerrFileAccess);
	}
m_bHdrDirty = false;
return(eBSFSuccess);
}

bool
CSfxArray::IsSOLiD(void)
{
return(m_bColorspace);
}

// TransformToColorspace
// transformation of basespace and colorspace (SOLiD)
//
uint8_t CSfxArray::CvrtSOLiDmap[5][7] = {
	{0,1,2,3,4,(uint8_t)eBaseEOS,(uint8_t)eBaseEOG},	// PrevBase==etBaseA
	{1,0,3,2,4,(uint8_t)eBaseEOS,(uint8_t)eBaseEOG},	// PrevBase==etBaseC
	{2,3,0,1,4,(uint8_t)eBaseEOS,(uint8_t)eBaseEOG},	// PrevBase==etBaseG
	{3,2,1,0,4,(uint8_t)eBaseEOS,(uint8_t)eBaseEOG},	// PrevBase==etBaseT
	{0,1,2,3,4,(uint8_t)eBaseEOS,(uint8_t)eBaseEOG},	// PrevBase==etBaseN
	};

int
CSfxArray::TransformToColorspace(uint8_t *pSrcBases,	// basespace sequence to transform
					int64_t SeqLen,				// sequence length
					uint8_t *pDstColorspace,		// where to write transformed
					bool bHiNibble,				// if false then write into bits 0..3, otherwise bits 4..7, of destination
					etSeqBase Seed)				// use this base as the seed base immediately preceding the first base
{
etSeqBase PrvBase;
etSeqBase CurBase;
uint8_t ColorSpace;
if(pSrcBases == NULL || SeqLen == 0 || pDstColorspace == NULL || Seed > eBaseN)
	return(eBSFerrParams);
PrvBase = Seed;
while(SeqLen--)
	{
	CurBase = *pSrcBases++ & 0x0f;						// conversion is on the low nibble
	*pDstColorspace &= bHiNibble ? 0x0f : 0xf0;
	if(!(CurBase == eBaseEOS || CurBase == eBaseEOG) && PrvBase <= 4)
		{
		if(CurBase > eBaseN)
			CurBase = eBaseN;
		else
			CurBase &= ~cRptMskFlg;				// remove any repeat masking
		ColorSpace = CvrtSOLiDmap[PrvBase][CurBase];
		*pDstColorspace++ |= bHiNibble ? ColorSpace << 4 : ColorSpace;
		}
	else
		*pDstColorspace++ |= bHiNibble ? CurBase << 4 : CurBase;
	PrvBase = CurBase;
	}
return(eBSFSuccess);
}

int
CSfxArray::TransformToBasespace(uint8_t *pSrcColorspace,	// colorspace sequence to transform
					int64_t SeqLen,				// sequence length
					uint8_t *pDstBases,			// where to write transformed
					bool bHiNibble,				// if false then write into bits 0..3, otherwise bits 4..7, of destination
					etSeqBase Seed)				// use this base as the seed base immediately preceding the first base
{
etSeqBase PrvBase;
etSeqBase CurBase;
uint8_t Base;
if(pSrcColorspace == NULL || SeqLen == 0 || pDstBases == NULL || Seed > eBaseN)
	return(eBSFerrParams);
PrvBase = Seed;
while(SeqLen--)
	{
	CurBase = *pSrcColorspace++ & 0x0f;						// conversion is on the low nibble
	*pDstBases &= bHiNibble ? 0x0f : 0xf0;
	if(CurBase <= eBaseN && PrvBase <= eBaseN)
		Base = CvrtSOLiDmap[PrvBase][CurBase];
	else
		Base = CurBase;
	*pDstBases++ |= bHiNibble ? Base << 4 : Base;
	PrvBase = Base;
	}
return(eBSFSuccess);
}

// obtain a copy of the header for external diagnostics
tsSfxHeaderV3 *
CSfxArray::GetSfxHeader(tsSfxHeaderV3 *pCopyTo)
{
if(pCopyTo != NULL)
	{
	memmove(pCopyTo,&m_SfxHeader,sizeof(tsSfxHeaderV3));
	return(pCopyTo);
	}
return(&m_SfxHeader);
}

teBSFrsltCodes
CSfxArray::SfxBlock2Disk(void)
{
teBSFrsltCodes Rslt;
uint8_t *pBisBases;
uint8_t *pBases;
uint64_t Idx;

uint64_t WrtLen;
// only need to write if something to write!
if(m_pSfxBlock == NULL || m_pSfxBlock->BlockID == 0 || m_pSfxBlock->ConcatSeqLen == 0 || m_pSfxBlock->NumEntries == 0)
	return(eBSFSuccess);

// if bisulfite processing then need to map all occurences of eBaseT to be eBaseC, and all occurances of eBaseA to be eBaseG
// do the sort, and then restore back to original values. Means that the memory requirements are greatly increased....
//CSAIS SAIS;
if(m_bBisulfite)
	{
	pBisBases = m_pBisulfateBases;
	pBases = m_pSfxBlock->SeqSuffix;
	for(Idx = 0; Idx < m_pSfxBlock->ConcatSeqLen; Idx++,pBases++,pBisBases++)
		{
		if(*pBases == eBaseT)
			*pBisBases = eBaseC;
		else
			{
			if(*pBases == eBaseA)
				*pBisBases = eBaseG;
			else
				*pBisBases = *pBases;
			}
		}
	if(m_bColorspace)
		{
		TransformToColorspace(m_pBisulfateBases,m_pSfxBlock->ConcatSeqLen,m_pBisulfateBases);
		TransformToColorspace(m_pSfxBlock->SeqSuffix,m_pSfxBlock->ConcatSeqLen,m_pSfxBlock->SeqSuffix);
		}

	QSortSeq((int64_t)m_pSfxBlock->ConcatSeqLen,m_pBisulfateBases,m_pSfxBlock->SfxElSize,(void *)&m_pSfxBlock->SeqSuffix[m_pSfxBlock->ConcatSeqLen]);
	}
else
	{
	if(m_bColorspace)
		TransformToColorspace(m_pSfxBlock->SeqSuffix,m_pSfxBlock->ConcatSeqLen,m_pSfxBlock->SeqSuffix);
	QSortSeq(m_pSfxBlock->ConcatSeqLen,m_pSfxBlock->SeqSuffix,m_pSfxBlock->SfxElSize,(void *)&m_pSfxBlock->SeqSuffix[m_pSfxBlock->ConcatSeqLen]);
	}

if (m_bColorspace)	// set hi nibbles of sequence to be original sequence
	TransformToBasespace(m_pSfxBlock->SeqSuffix, m_pSfxBlock->ConcatSeqLen, m_pSfxBlock->SeqSuffix, true);

if (!m_bInMemSfx)
	{
	// set block size and file offset for suffix block into header
	m_SfxHeader.NumSfxBlocks = 1;
	m_SfxHeader.SfxBlockSize = sizeof(tsSfxBlock) + m_pSfxBlock->ConcatSeqLen - 1 + ((size_t)m_pSfxBlock->ConcatSeqLen * m_pSfxBlock->SfxElSize);
	m_SfxHeader.SfxBlockOfs = m_SfxHeader.FileLen;

	// now write...
	WrtLen = sizeof(tsSfxBlock) + m_pSfxBlock->ConcatSeqLen - 1;
	if((Rslt=ChunkedWrite(m_SfxHeader.FileLen,(uint8_t *)m_pSfxBlock,WrtLen))!=eBSFSuccess)
		{
		AddErrMsg("CSfxArray::SfxBlock2Disk","Unable to write suffix block sequence to disk");
		Reset(false);
		return(Rslt);
		}

	m_SfxHeader.FileLen += WrtLen;

	WrtLen = (int64_t)m_pSfxBlock->SfxElSize * m_pSfxBlock->ConcatSeqLen;
	if((Rslt=ChunkedWrite(m_SfxHeader.FileLen,(uint8_t *)&m_pSfxBlock->SeqSuffix[m_pSfxBlock->ConcatSeqLen],WrtLen))!=eBSFSuccess)
		{
		AddErrMsg("CSfxArray::SfxBlock2Disk","Unable to write suffix block array to disk");
		Reset(false);
		return(Rslt);
		}
	m_SfxHeader.FileLen += WrtLen;

	m_pSfxBlock->BlockID = 0;
	m_pSfxBlock->NumEntries = 0;
	m_pSfxBlock->ConcatSeqLen = 0;
	}
return(eBSFSuccess);
}

teBSFrsltCodes
CSfxArray::Entries2Disk(void)
{
teBSFrsltCodes Rslt;

// any entries to write?
if(m_pEntriesBlock == NULL || !m_pEntriesBlock->NumEntries)
	{
	if(m_pEntriesBlock != NULL)
		{
#ifdef _WIN32
		free(m_pEntriesBlock);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
		if(m_pEntriesBlock != MAP_FAILED)
			munmap(m_pEntriesBlock,m_AllocEntriesBlockMem);
#endif
		m_pEntriesBlock = NULL;
		}
	m_SfxHeader.EntriesSize = 0;
	m_SfxHeader.EntriesOfs = 0;
	m_bHdrDirty = true;
	return(eBSFSuccess);
	}

// now write...
m_pEntriesBlock->MaxEntries = m_pEntriesBlock->NumEntries;
m_SfxHeader.EntriesSize = sizeof(tsSfxEntriesBlock) + (sizeof(tsSfxEntry) * (m_pEntriesBlock->NumEntries - 1));

if(!m_bInMemSfx)
	{
	m_SfxHeader.EntriesOfs = m_SfxHeader.FileLen;
	if((Rslt=ChunkedWrite(m_SfxHeader.EntriesOfs,(uint8_t *)m_pEntriesBlock,(int64_t)m_SfxHeader.EntriesSize))!=eBSFSuccess)
		{
		AddErrMsg("CSfxArray::AddEntry","Unable to write entries block to disk");
		Reset(false);
		return(Rslt);
		}
	m_SfxHeader.FileLen += m_SfxHeader.EntriesSize;
	}
m_bHdrDirty = true;
return(eBSFSuccess);
}


// Disk2Hdr
// Loads in header from disk, checks if header was created using older release and if so then updates to latest
teBSFrsltCodes
CSfxArray::Disk2Hdr(char *pszFile)
{
uint8_t HdrVer[8];
uint32_t Version;
tsSfxHeaderVv SfxHeaderVv;

if(_lseeki64(m_hFile,0,SEEK_SET)!=0)			// header is at start of file
	{
	AddErrMsg("CSfxArray::Disk2Hdr","Seek failed to offset 0 on %s - %s",pszFile,strerror(errno));
	Reset(false);			// closes opened file..
	return(eBSFerrFileAccess);
	}

// start by reading in the 1st 8 bytes, this covers the version
if(8 != read(m_hFile,&HdrVer,8))
	{
	AddErrMsg("CSfxArray::Disk2Hdr","Read of file header failed on %s - %s",pszFile,strerror(errno));
	Reset(false);			// closes opened file..
	return(eBSFerrFileAccess);
	}

// header in memory, validate it as being either a V3 or V4 file header
if(tolower(HdrVer[0]) != 's' ||
	tolower(HdrVer[1]) != 'f' ||
	tolower(HdrVer[2]) != 'x' ||
	(tolower(HdrVer[3]) < '3' || tolower(HdrVer[3]) > '5'))
	{
	AddErrMsg("CSfxArray::Disk2Hdr","%s opened but invalid magic signature - not a 'kit4b index' generated suffix array file",pszFile);
	Reset(false);			// closes opened file..
	return(eBSFerrNotBioseq);
	}

// can we handle this file structure version?
Version = *(uint32_t *)&HdrVer[4];
if(Version < cSFXVersionBack || Version > cSFXVersion)
	{
	AddErrMsg("CSfxArray::Disk2Hdr","%s opened as 'kit4b index' generated suffix array file but structure version %d is incompatiable with this release %d of kit4b",
				pszFile, Version, cSFXVersion);
	Reset(false);			// closes opened file..
	return(eBSFerrFileVer);
	}

m_bV3File = Version == cSFXVersionBack ? true : false;

if(_lseeki64(m_hFile,0,SEEK_SET)!=0)			// header is at start of file
	{
	AddErrMsg("CSfxArray::Disk2Hdr","Seek failed to offset 0 on %s - %s",pszFile,strerror(errno));
	Reset(false);			// closes opened file..
	return(eBSFerrFileAccess);
	}

if(m_bV3File)
	{
	if(sizeof(tsSfxHeaderVv) != read(m_hFile,&SfxHeaderVv,sizeof(tsSfxHeaderVv)))
		{
		AddErrMsg("CSfxArray::Disk2Hdr","Read of V3 file header failed on %s - %s",pszFile,strerror(errno));
		Reset(false);			// closes opened file..
		return(eBSFerrFileAccess);
		}
	memcpy(&m_SfxHeader,&SfxHeaderVv,sizeof(tsSfxHeaderVv));
	m_SfxHeader.Magic[3] = '5';
	m_SfxHeader.Version = cSFXVersion;
	memcpy(&m_SfxHeader.szDescription,&SfxHeaderVv.szDescription,sizeof(SfxHeaderVv.szDescription));
	memcpy(&m_SfxHeader.szTitle,&SfxHeaderVv.szTitle,sizeof(SfxHeaderVv.szTitle));
	}
else
	{
	if(sizeof(tsSfxHeaderV3) != read(m_hFile,&m_SfxHeader,sizeof(tsSfxHeaderV3)))
		{
		AddErrMsg("CSfxArray::Disk2Hdr","Read of V4 file header failed on %s - %s",pszFile,strerror(errno));
		Reset(false);			// closes opened file..
		return(eBSFerrFileAccess);
		}
	}


m_bBisulfite = m_SfxHeader.Attributes & 0x01 ? true : false;
m_bColorspace = m_SfxHeader.Attributes & 0x02 ? true : false;
m_bHdrDirty = false;
return(eBSFSuccess);
}

// Disk2Entries
// Loads any entries into memory from file
teBSFrsltCodes
CSfxArray::Disk2Entries(void)
{
teBSFrsltCodes Rslt;

uint32_t Idx;
tsSfxEntryV3 *pEntryV3;
tsSfxEntry *pEntry;
tsSfxEntriesBlockV3 *pSfxEntriesBlockV3;

if(m_pEntriesBlock != NULL)
	{
#ifdef _WIN32
	free(m_pEntriesBlock);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pEntriesBlock != MAP_FAILED)
		munmap(m_pEntriesBlock,m_AllocEntriesBlockMem);
#endif
	m_pEntriesBlock = NULL;
	}

// any entries to load?
if(m_SfxHeader.EntriesOfs == 0 || m_SfxHeader.EntriesSize == 0)
	return(eBSFSuccess);

if(m_bV3File)
	{
	if((pSfxEntriesBlockV3 = (tsSfxEntriesBlockV3 *)new uint8_t [m_SfxHeader.EntriesSize])==NULL)
		{
		AddErrMsg("CSfxArray::Disk2Entries","unable to allocate %d bytes for holding V3 entries block",m_SfxHeader.EntriesSize);
		Reset(false);			// closes opened file..
		return(eBSFerrMem);
		}

	if((Rslt=ChunkedRead(m_SfxHeader.EntriesOfs,(uint8_t *)pSfxEntriesBlockV3,m_SfxHeader.EntriesSize))!=eBSFSuccess)
		{
		AddErrMsg("CSfxArray::Disk2Entries","unable to load V3 entries block of length %d from offset %d",m_SfxHeader.EntriesSize,m_SfxHeader.EntriesOfs);
		Reset(false);			// closes opened file..
		return(eBSFerrFileAccess);
		}

	m_SfxHeader.EntriesSize = sizeof(tsSfxEntriesBlock) + (pSfxEntriesBlockV3->NumEntries * sizeof(tsSfxEntry));

#ifdef _WIN32
	m_pEntriesBlock = (tsSfxEntriesBlock *) malloc(m_SfxHeader.EntriesSize);
	if(m_pEntriesBlock==NULL)
		{
		AddErrMsg("CSfxArray::Disk2Entries","unable to allocate %u bytes for holding entries block",m_SfxHeader.EntriesSize);
		delete pSfxEntriesBlockV3;
		Reset(false);			// closes opened file..
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and seems to have issues if more than 2GB allocation
	m_pEntriesBlock = (tsSfxEntriesBlock *)mmap(NULL,m_SfxHeader.EntriesSize, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pEntriesBlock == MAP_FAILED)
		{
		AddErrMsg("CSfxArray::Disk2Entries","unable to allocate %u bytes for holding entries block",m_SfxHeader.EntriesSize);
		m_pEntriesBlock = NULL;
		delete pSfxEntriesBlockV3;
		Reset(false);			// closes opened file..
		return(eBSFerrMem);
		}
#endif
	m_pEntriesBlock->NumEntries = pSfxEntriesBlockV3->NumEntries;
	m_pEntriesBlock->MaxEntries = pSfxEntriesBlockV3->NumEntries;
	pEntryV3 = pSfxEntriesBlockV3->Entries;
	pEntry = m_pEntriesBlock->Entries;
	for(Idx = 0; Idx < (uint32_t)m_pEntriesBlock->MaxEntries; Idx++, pEntryV3++,pEntry++)
		{
		pEntry->EntryID = pEntryV3->EntryID;						// identifies each entry (1..n), unique within this suffix file
		pEntry->fBlockID = pEntryV3->fBlockID;					// which SfxBlock contains sequence for this entry
		memmove(pEntry->szSeqName,pEntryV3->szSeqName,sizeof(pEntry->szSeqName));	// entry name - typically a chromosome name
		pEntry->NameHash = pEntryV3->NameHash;					// hash over szSeqName
		pEntry->SeqLen = pEntryV3->SeqLen;						// sequence length - excludes sequence terminator eSeqEOS
		pEntry->StartOfs = (int64_t)(uint64_t)pEntryV3->StartOfs;						// offset into concatenated sequences of where this sequence starts
		pEntry->EndOfs = (int64_t)(uint64_t)pEntryV3->EndOfs;						// offset into concatenated sequences of where this sequence ends
		}
	delete []pSfxEntriesBlockV3;
	}
else
	{ 

#ifdef _WIN32
	m_pEntriesBlock = (tsSfxEntriesBlock *) malloc(m_SfxHeader.EntriesSize);
	if(m_pEntriesBlock==NULL)
		{
		AddErrMsg("CSfxArray::Disk2Entries","unable to allocate %u bytes for holding entries block",m_SfxHeader.EntriesSize);
		Reset(false);			// closes opened file..
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and seems to have issues if more than 2GB allocation
	m_pEntriesBlock = (tsSfxEntriesBlock *)mmap(NULL,m_SfxHeader.EntriesSize, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pEntriesBlock == MAP_FAILED)
		{
		AddErrMsg("CSfxArray::Disk2Entries","unable to allocate %u bytes for holding entries block",m_SfxHeader.EntriesSize);
		m_pEntriesBlock = NULL;
		Reset(false);			// closes opened file..
		return(eBSFerrMem);
		}
#endif

	if((Rslt=ChunkedRead(m_SfxHeader.EntriesOfs,(uint8_t *)m_pEntriesBlock,m_SfxHeader.EntriesSize))!=eBSFSuccess)
		{
		AddErrMsg("CSfxArray::Disk2Entries","unable to load entries block of length %d from offset %d",m_SfxHeader.EntriesSize,m_SfxHeader.EntriesOfs);
		Reset(false);			// closes opened file..
		return(eBSFerrFileAccess);
		}
	}

return(eBSFSuccess);
}


// ChunkedWrite
// Seeks to specified 64bit file offset and writes to disk as chunks of no more than INT_MAX/16
teBSFrsltCodes
CSfxArray::ChunkedWrite(int64_t WrtOfs,uint8_t *pData,int64_t WrtLen)
{
int BlockLen;
if(_lseeki64(m_hFile,WrtOfs,SEEK_SET) != WrtOfs)
	{
	AddErrMsg("CSfxArray::ChunkedWrite","Unable to seek to %zd on file %s - error %s",WrtOfs,m_szFile,strerror(errno));
	Reset(false);
	return(eBSFerrFileAccess);
	}

while(WrtLen)
	{
	BlockLen = WrtLen > (int64_t)(INT_MAX/16) ? (INT_MAX/16) : (int)WrtLen;
	WrtLen -= BlockLen;
	if(write(m_hFile,pData,BlockLen)!=BlockLen)
		{
		AddErrMsg("CSfxArray::ChunkedWrite","Unable to write to disk on file %s - error %s",m_szFile,strerror(errno));
		Reset(false);
		return(eBSFerrFileAccess);
		}
	pData += BlockLen;
	}
return(eBSFSuccess);
}

// ChunkedRead
// Seeks to specified 64bit file offset and reads from disk as chunks of no more than INT_MAX/32
// This allows checks for thread termination requests at reasonable time intervals
// returns (teBSFrsltCodes)1 if m_bTermThread was set by another thread
teBSFrsltCodes
CSfxArray::ChunkedRead(int64_t RdOfs,uint8_t *pData,int64_t RdLen)
{
uint32_t BlockLen;
uint32_t ActualRdLen;
if(_lseeki64(m_hFile,RdOfs,SEEK_SET) != RdOfs)
	{
	AddErrMsg("CSfxArray::ChunkedRead","Unable to seek to %zd on file %s - error %s",RdOfs,m_szFile,strerror(errno));
	return(eBSFerrFileAccess);
	}
if(m_bTermThread)
	return((teBSFrsltCodes)1);

while(RdLen)
	{
	BlockLen = RdLen > (int64_t)(INT_MAX/2) ? (INT_MAX/2) : (uint32_t)RdLen;
	RdLen -= BlockLen;
	if((ActualRdLen = read(m_hFile,pData,BlockLen))!=BlockLen)
		{
		AddErrMsg("CSfxArray::ChunkedRead","Unable to read from disk file %s - error %s",m_szFile,strerror(errno));
		return(eBSFerrFileAccess);
		}
	pData += BlockLen;
	if(m_bTermThread)
		return((teBSFrsltCodes)1);
	}
return(eBSFSuccess);
}

// Flush2Disk
// Writes to disk any outstanding suffix block, then the entries block, and finally the file header
teBSFrsltCodes
CSfxArray::Flush2Disk(void)			// flush and commit to disk
{
teBSFrsltCodes Rslt;
int64_t ReallocSize;

if (!m_bInMemSfx && m_hFile == -1)
	return(eBSFSuccess);

if(!m_bCreate)
	return(eBSFerrParams);

// now that actual sequence length is known then need to allocate for suffix array elements
if(m_pSfxBlock->ConcatSeqLen)
	{
	if(m_pSfxBlock->ConcatSeqLen < cThres8ByteSfxEls)			// if less than cThres8ByteSfxEls then can use 4 bytes per suffix element otherwise need to go to 5 bytes
		m_pSfxBlock->SfxElSize = 4;
	else
		m_pSfxBlock->SfxElSize = 5;

	ReallocSize = (int64_t)sizeof(tsSfxBlock) + m_pSfxBlock->ConcatSeqLen + (m_pSfxBlock->ConcatSeqLen * m_pSfxBlock->SfxElSize);

		// m_pSfxBlock almost certainly needs to be extended
	if(ReallocSize > (int64_t)m_AllocSfxBlockMem) 
		{
		tsSfxBlock *pRealloc;
#ifdef _WIN32
		pRealloc = (tsSfxBlock *)realloc(m_pSfxBlock,(size_t)ReallocSize);
#else
		pRealloc = (tsSfxBlock *)mremap(m_pSfxBlock,m_AllocSfxBlockMem,(size_t)ReallocSize,MREMAP_MAYMOVE);
		if(pRealloc == MAP_FAILED)
			pRealloc = NULL;
#endif
		if(pRealloc == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: SfxBlock memory re-allocation to %zd bytes - %s",(int64_t)ReallocSize,strerror(errno));
			return(eBSFerrMem);
			}

		m_pSfxBlock = pRealloc;
		m_AllocSfxBlockMem = ReallocSize;
		}
	}

if((Rslt=SfxBlock2Disk())!=eBSFSuccess) // sorting is immediately prior to writing to disk
	return(Rslt);
if((Rslt=Entries2Disk())!=eBSFSuccess)
	return(Rslt);

if (!m_bInMemSfx)
	{
	if((Rslt=Hdr2Disk())!=eBSFSuccess)
		return(Rslt);
	}
return(eBSFSuccess);
}

// Memory resident suffix array
// User creates the suffix array in memory and can immediately access for alignments etc without needing to write to file
int
CSfxArray::Open(bool bBisulfite,
				  bool bColorspace)
{
Reset(false);						// reset context in case any previously accessed file still opened
InitHdr();
m_bCreate = true;
m_bBisulfite = bBisulfite;
m_bColorspace = bColorspace;
m_bInMemSfx = true;
return(eBSFSuccess);
}

// Open()
// Opens specified suffix file pszFile
// Option to create or truncate pszFile
// Option to process for bisulfite indexing
// Option to process for colorspace
int
CSfxArray::Open(char *pszFile,
					  bool bCreate,
					  bool bBisulfite,
					  bool bColorspace)
{
teBSFrsltCodes Rslt;
if(pszFile == NULL || *pszFile == '\0') // validate parameters
	return(eBSFerrParams);

Reset(false);				// reset context in case file still opened

#ifdef _WIN32
if(!bCreate)
	m_hFile = open(pszFile, O_READSEQ ); // file access is normally sequential..
else
	m_hFile = open(pszFile, ( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE) );
#else
if(!bCreate)
	m_hFile = open64(pszFile, O_READSEQ ); // file access is normally sequential..
else
	if((m_hFile = open64(pszFile,O_WRONLY | O_CREAT, S_IREAD | S_IWRITE))!=1)
	  if(ftruncate(m_hFile,0)!=0)
			{
			AddErrMsg("CSfxArray::Open","Unable to truncate %s - %s",pszFile,strerror(errno));
			Rslt = bCreate ? eBSFerrCreateFile : eBSFerrOpnFile;
			return(Rslt);
			}

#endif

if(m_hFile == -1)					// check if file open succeeded
	{
	AddErrMsg("CSfxArray::Open","Unable to open %s - %s",pszFile,strerror(errno));
	Rslt = bCreate ? eBSFerrCreateFile : eBSFerrOpnFile;
	Reset(false);
	return(Rslt);
	}

strncpy(m_szFile,pszFile,_MAX_PATH);
m_szFile[_MAX_PATH-1] = '\0';


if(bCreate)
	{
	InitHdr();
	if((Rslt = Hdr2Disk()) != eBSFSuccess)
		{
		Reset(false);			// closes opened file..
		return(Rslt);
		}
	m_bCreate = true;
	m_bBisulfite = bBisulfite;
	m_bColorspace = bColorspace;
	}
else // else opening existing file
	{
	// load header
	if((Rslt=Disk2Hdr(pszFile))!=eBSFSuccess)
		{
		Reset(false);			// closes opened file..
		return(eBSFerrFileAccess);
		}
	// check if file was appropriate for bisulfite
	if(m_bBisulfite != bBisulfite)
		{
		AddErrMsg("CSfxArray::Open","BiosfxV3 file '%s' %s bisulfite index",pszFile,
			m_bBisulfite ? "was generated with" : "does not contain");

		Reset(false);			// closes opened file..
		return(eBSFerrFileType);
		}
	// check if file was appropriate for colorspace
	if(m_bColorspace != bColorspace)
		{
		AddErrMsg("CSfxArray::Open","BiosfxV3 file '%s' %s colorspace index",pszFile,
			m_bColorspace ? "was generated with" : "does not contain");

		Reset(false);			// closes opened file..
		return(eBSFerrFileType);
		}

	// load any directory entries
	if((Rslt=Disk2Entries()) < eBSFSuccess)
		{
		Reset(false);			// closes opened file..
		return(eBSFerrFileAccess);
		}

	// allocate suffix block memory
#ifdef _WIN32
	m_pSfxBlock = (tsSfxBlock *) malloc((size_t)m_SfxHeader.SfxBlockSize);
	if(m_pSfxBlock == NULL)
		{
		AddErrMsg("CSfxArray::Open","Fatal: unable to allocate %zd bytes contiguous memory for index",(int64_t)m_SfxHeader.SfxBlockSize);
		Reset(false);
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and seems to have issues if more than 2GB allocation
	m_pSfxBlock = (tsSfxBlock *)mmap(NULL,m_SfxHeader.SfxBlockSize, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pSfxBlock == MAP_FAILED)
		{
		AddErrMsg("CSfxArray::Open","Fatal: unable to allocate block memory");
		m_pSfxBlock = NULL;
		Reset(false);
		return(eBSFerrMem);
		}
#endif
	m_AllocSfxBlockMem = (size_t)m_SfxHeader.SfxBlockSize;
	m_pSfxBlock->BlockID = 0;
	m_pSfxBlock->NumEntries = 0;
	m_pSfxBlock->ConcatSeqLen = 0;

	// initialise thread for loading suffix blocks as a background task
	// this initialisation must be treated as being a transacted task - either completes or is rolled back
#ifdef _WIN32
	if((m_JobMutex = CreateMutex(NULL,false,NULL))==NULL)
#else
	if(pthread_mutex_init (&m_JobMutex,NULL)!=0)
#endif
		{
		AddErrMsg("CSfxArray::Open","Fatal: unable to create m_JobMutex");
		Reset(false);			// closes opened file..
		return(eBSFerrInternal);
		}

#ifdef _WIN32
	if((m_JobReqEvent = CreateEvent(NULL,false,false,NULL))==NULL)
#else
	if(pthread_cond_init(&m_JobReqEvent,NULL)!=0)
#endif
		{
		AddErrMsg("CSfxArray::Open","Fatal: unable to create m_JobReqEvent");
#ifdef _WIN32
		CloseHandle(m_JobMutex);
		m_JobMutex = NULL;
#else
		pthread_mutex_destroy(&m_JobMutex);
#endif
		Reset(false);			// closes opened file..
		return(eBSFerrInternal);
		}

#ifdef _WIN32
	if((m_JobAckEvent = CreateEvent(NULL,false,false,NULL))==NULL)
#else
	if(pthread_cond_init(&m_JobAckEvent,NULL)!=0)
#endif
		{
		AddErrMsg("CSfxArray::Open","Fatal: unable to create m_JobAckEvent");
#ifdef _WIN32
		CloseHandle(m_JobMutex);
		m_JobMutex = NULL;
		CloseHandle(m_JobReqEvent);
		m_JobReqEvent = NULL;
#else
		pthread_mutex_destroy(&m_JobMutex);
		pthread_cond_destroy(&m_JobReqEvent);
#endif
		Reset(false);			// closes opened file..
		return(eBSFerrInternal);
		}

m_CASSeqFlags = 0;

	// now to create/start the actual background readahead thread
#ifdef _WIN32
	if((m_threadHandle = (HANDLE)_beginthreadex(NULL,0x0fffff,ThreadedPrereadBlocks,this,0,&m_threadID))==NULL)
		{
		AddErrMsg("CSfxArray::Open","Fatal: unable to create background readahead thread");
		CloseHandle(m_JobMutex);
		m_JobMutex = NULL;
		CloseHandle(m_JobReqEvent);
		m_JobReqEvent = NULL;
		CloseHandle(m_JobAckEvent);
		m_JobAckEvent = NULL;
		Reset(false);			// closes opened file..
		return(eBSFerrInternal);
		}
	m_bThreadActive = true;
	ResumeThread(m_threadHandle);
#else
	m_bThreadActive = true;			// assume thread will be created, set false if creation fails
	if((m_threadRslt =	pthread_create (&m_threadID , NULL , ThreadedPrereadBlocks , this ))!=0)
		{
		m_bThreadActive = false;
		AddErrMsg("CSfxArray::Open","Fatal: unable to create background readahead thread");
		Reset(false);			// closes opened file..
		return(eBSFerrInternal);
		}

#endif
	// allow a few seconds for thread to have started and be initialised ready to load blocks
#ifdef _WIN32
	Sleep(2000);
#else
	sleep(2);
#endif
	// thread started, request that initial block be preloaded
	if(m_SfxHeader.NumSfxBlocks > 0)
		{
		m_ReqBlockID = 1;
		m_ReqBlockRslt = (teBSFrsltCodes)1;
#ifdef _WIN32
		SetEvent(m_JobReqEvent);
#else
		pthread_cond_signal(&m_JobReqEvent);
#endif
		}
	}

return(eBSFSuccess);
}

// Fasta2InMemSfxArray
// Parse input fasta format file into an in-memory suffix array
// Suffix array can then be accessed as if read from a pre-generated suffix array file
int
CSfxArray::Fasta2InMemSfxArray(char *pszFile,	// generate in-memory suffix array from this fasta file
	uint32_t MinSeqLen,				// skipping sequences which are less than this length
	uint32_t MaxSeqLen,				// skipping sequences which are longer than this length
	bool bBisulfite,				// true if Bisulfite
	bool bColorspace)				// true if colorspace
{
CFasta Fasta;
uint8_t *pSeqBuff;
uint8_t *pMskBase;
uint32_t MskIdx;
size_t BuffOfs;
size_t AllocdBuffSize;
size_t AvailBuffSize;
char szName[cBSFSourceSize];
char szDescription[cBSFDescriptionSize];
uint32_t SeqLen;
int Descrlen;
bool bFirstEntry;
bool bEntryCreated;
int Rslt;
int SeqID;
uint32_t NumSeqsUnderlength;
uint32_t NumSeqsOverlength;

gDiagnostics.DiagOut(eDLInfo, gszProcName, "Fasta2InMemSfxArray: generating in-memory suffix array from fasta format file '%s'", pszFile);

// if not specified then default MinSeqLen to be cMinAllowSeqLen and MaxSeqLen to be cMaxAllowSeqLen
if(MinSeqLen == 0)
	MinSeqLen = cMinAllowInMemSeqLen;
if(MaxSeqLen == 0)
	MaxSeqLen = cMaxAllowInMemSeqLen;

if(pszFile == NULL || pszFile[0] == '\0')
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Fasta2InMemSfxArray: fasta filename not specified");
	return(eBSFerrParams);
	}

if(MinSeqLen < cMinAllowInMemSeqLen || MaxSeqLen > cMaxAllowInMemSeqLen ||
	MinSeqLen > MaxSeqLen)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Fasta2InMemSfxArray: MinSeqLen (%u) and/or MaxSeqLen (%u) parameter error", MinSeqLen, MaxSeqLen);
	return(eBSFerrParams);
	}

// expecting input file to be in fasta format
if ((Rslt = Fasta.Open(pszFile, true)) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Fasta2InMemSfxArray: Unable to open '%s' [%s] %s", pszFile, Fasta.ErrText((teBSFrsltCodes)Rslt), Fasta.GetErrMsg());
	return(Rslt);
	}

// ensure initiated for in memory suffix processing
Open(bBisulfite, bColorspace);

AllocdBuffSize = (size_t)cMaxMemSfxSeqAlloc;
// note malloc is used as can then simply realloc to expand as may later be required
if ((pSeqBuff = (uint8_t *)malloc(AllocdBuffSize)) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Fasta2InMemSfxArray: Unable to allocate memory (%u bytes) for sequence buffer", (uint32_t)AllocdBuffSize);
	Fasta.Close();
	return(eBSFerrMem);
	}
AvailBuffSize = AllocdBuffSize;

bFirstEntry = true;
bEntryCreated = false;
SeqID = 0;
BuffOfs = 0;
NumSeqsUnderlength = 0;
NumSeqsOverlength = 0;
while ((Rslt = SeqLen = Fasta.ReadSequence(&pSeqBuff[BuffOfs], (int)min(AvailBuffSize, (size_t)cMaxMemSfxSeqAlloc), true, false)) > eBSFSuccess)
	{
	if (SeqLen == eBSFFastaDescr)		// just read a descriptor line
		{
		SeqID++;
		if (bEntryCreated)				// add any previous entry
			{
			if (BuffOfs < (size_t)MinSeqLen || BuffOfs >(size_t)MaxSeqLen)
				{
				if(BuffOfs < (size_t)MinSeqLen)
					NumSeqsUnderlength += 1;
				else
					NumSeqsOverlength += 1;
				}
			else
				if ((Rslt = AddEntry(szName, pSeqBuff, (uint32_t)BuffOfs)) < eBSFSuccess)
					{
					gDiagnostics.DiagOut(eDLFatal, gszProcName, "Fasta2InMemSfxArray: AddEntry() error %d %s", Rslt, GetErrMsg());
					if (pSeqBuff != NULL)
						free(pSeqBuff);
					Reset(false);
					return(Rslt);
					}
			Rslt = eBSFSuccess;
			}
		Descrlen = Fasta.ReadDescriptor(szDescription, cBSFDescriptionSize);
			// An assumption - will one day bite real hard - is that the
			// fasta descriptor line starts with some form of unique identifier.
			// Use this identifier as the entry name.
		if (sscanf(szDescription, " %s[ ,]", szName) != 1)
			sprintf(szName, "%s.%d", pszFile, ++SeqID);

		bFirstEntry = false;
		bEntryCreated = true;
		BuffOfs = 0;
		continue;
		}
	else
		if (bFirstEntry)	// if there was no descriptor then dummy up one...
			{
			SeqID++;
			sprintf(szName, "%s.%d", pszFile, SeqID);
			strcpy(szDescription, "No Description provided");
			bFirstEntry = false;
			bEntryCreated = true;
			}

		// remove any repeat masking flags so that sorts can actually sort
		// if run of more than 25 Ns and at least 5 Ns to end of buffer then randomly mutate
		// every 13th N
		//	e.g <25Ns>r<12Ns>r<12Ns> where r is a pseudorandom base
	pMskBase = &pSeqBuff[BuffOfs];
	int SeqNs = 0;
	for (MskIdx = 0; MskIdx < SeqLen; MskIdx++, pMskBase++)
		{
		*pMskBase &= ~cRptMskFlg;
		if (*pMskBase == eBaseN && (MskIdx + 5) < SeqLen)
			{
			if (++SeqNs > 25 &&
					pMskBase[1] == eBaseN &&
					pMskBase[2] == eBaseN &&
					pMskBase[3] == eBaseN &&
					pMskBase[4] == eBaseN)
				{
				if (!(SeqNs % 13))	// mutate every 13th
					*pMskBase = rand() % 4;
				}
			}
		else
			SeqNs = 0;
		}

	BuffOfs += SeqLen;
	AvailBuffSize -= SeqLen;
	if (AvailBuffSize < (size_t)(cMaxMemSfxSeqAlloc / 8))
		{
		size_t NewSize = (size_t)cMaxMemSfxSeqAlloc + AllocdBuffSize;
		uint8_t *pTmp;
		if ((pTmp = (uint8_t *)realloc(pSeqBuff, NewSize)) == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Fasta2InMemSfxArray: Unable to reallocate memory (%u bytes) for sequence buffer", (uint32_t)NewSize);
			if (pSeqBuff != NULL)
				free(pSeqBuff);
			Reset(false);
			return(eBSFerrMem);
			}
		pSeqBuff = pTmp;
		AllocdBuffSize = NewSize;
		AvailBuffSize = AllocdBuffSize - BuffOfs;
		}
	}

if (Rslt >= eBSFSuccess && bEntryCreated && BuffOfs > 0)			// close entry
	{
	if (BuffOfs < (size_t)MinSeqLen || BuffOfs >(size_t)MaxSeqLen)
		{
		if (BuffOfs < (size_t)MinSeqLen)
			NumSeqsUnderlength += 1;
		else
			NumSeqsOverlength += 1;
		Rslt = eBSFSuccess;
		}
	else
		{
		if ((Rslt = AddEntry(szName, pSeqBuff, (uint32_t)BuffOfs)) < eBSFSuccess)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Fasta2InMemSfxArray: AddEntry error %d %s", Rslt, GetErrMsg());
			if (pSeqBuff != NULL)
				free(pSeqBuff);
			Reset(false);
			return(Rslt);
			}
		else
			Rslt = eBSFSuccess;
		}
	}
if (pSeqBuff != NULL)
	free(pSeqBuff);

if(NumSeqsUnderlength)
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Fasta2InMemSfxArray: NumSeqsUnderlength %u", NumSeqsUnderlength);
if (NumSeqsOverlength)
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Fasta2InMemSfxArray: NumSeqsOverlength %u", NumSeqsOverlength);

if (m_pEntriesBlock->NumEntries == 0)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Fasta2InMemSfxArray: No sequences accepted for indexing");
	Close(false);
	return(eBSFerrEntry);
	}

int NumDupEntries;
uint32_t DupEntries[10];
char szDupEntry[100];

// check for any duplicate entry names
if ((NumDupEntries = ChkDupEntries(10, &DupEntries[0])) > 0)
	{
	while (NumDupEntries--)
		{
		GetIdentName(DupEntries[NumDupEntries], sizeof(szDupEntry) - 1, szDupEntry); // get sequence name for specified entry identifier
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Fasta2InMemSfxArray: duplicate sequence entry name '%s' in file '%s'", szDupEntry, pszFile);
		}
	Close(false);
	return(eBSFerrFastqSeqID);
	}

gDiagnostics.DiagOut(eDLInfo, gszProcName, "Fasta2InMemSfxArray:  sorting suffix array...");
Finalise();
gDiagnostics.DiagOut(eDLInfo, gszProcName, "Fasta2InMemSfxArray:  completed sorting suffix array");
return(eBSFSuccess);
}




int
CSfxArray::Close(bool bFlush)
{
return(Reset(bFlush));
}


teBSFrsltCodes
CSfxArray::SetDatasetName(char *pszDataset)
{
if (!m_bInMemSfx && m_hFile == -1)
	return(eBSFerrClosed);
if(!m_bCreate)
	return(eBSFerrRead);
if(pszDataset == NULL || pszDataset[0] == '\0')
	return(eBSFerrParams);
strncpy((char *)m_SfxHeader.szDatasetName,pszDataset,cMaxDatasetSpeciesChrom);
m_SfxHeader.szDatasetName[cMaxDatasetSpeciesChrom-1] = '\0';
return(eBSFSuccess);
}

char *
CSfxArray::GetDatasetName(void)
{
if (!m_bInMemSfx && m_hFile == -1)
	return(NULL);
if(m_bCreate)
	return(NULL);
return((char *)m_SfxHeader.szDatasetName);
}

teBSFrsltCodes
CSfxArray::SetDescription(char *pszDescription)
{
if(!m_bInMemSfx && m_hFile == -1)
	return(eBSFerrClosed);
if(!m_bCreate)
	return(eBSFerrRead);
if(pszDescription == NULL || pszDescription[0] == '\0')
	return(eBSFerrParams);
strncpy((char *)m_SfxHeader.szDescription,pszDescription,cMBSFFileDescrLen);
m_SfxHeader.szDescription[cMBSFFileDescrLen-1] = '\0';
return(eBSFSuccess);
}

teBSFrsltCodes
CSfxArray::GetDescription(int MaxLen,char *pszDescription)
{
if (!m_bInMemSfx && m_hFile == -1)
	return(eBSFerrClosed);
if(m_bCreate)
	return(eBSFerrWrite);
strncpy(pszDescription,(char *)m_SfxHeader.szDescription,MaxLen);
pszDescription[MaxLen-1] = '\0';
return(eBSFSuccess);
}


teBSFrsltCodes
CSfxArray::SetTitle(char *pszTitle)
{
if (!m_bInMemSfx && m_hFile == -1)
	return(eBSFerrClosed);
if(!m_bCreate)
	return(eBSFerrRead);
if(pszTitle == NULL || pszTitle[0] == '\0')
	return(eBSFerrParams);
strncpy((char *)m_SfxHeader.szTitle,pszTitle,cMBSFShortFileDescrLen);
m_SfxHeader.szTitle[cMBSFShortFileDescrLen-1] = '\0';
return(eBSFSuccess);
}


teBSFrsltCodes
CSfxArray::GetTitle(int MaxLen,char *pszTitle)
{
if (!m_bInMemSfx && m_hFile == -1)
	return(eBSFerrClosed);
if(m_bCreate)
	return(eBSFerrWrite);
strncpy(pszTitle,(char *)m_SfxHeader.szTitle,MaxLen);
pszTitle[MaxLen-1] = '\0';
return(eBSFSuccess);
}

// set maximum iterations limit when processing non-unique matching subsequences
int								// previous value of m_MaxIter
CSfxArray::SetMaxIter(int MaxIter)
{
int Prev = m_MaxIter;
m_MaxIter = MaxIter > 0 ? MaxIter : 0;
m_MaxKMerOccs = m_MaxIter;
return(Prev);
}

int								// current value of m_MaxIter
CSfxArray::GetMaxIter(void)
{
return(m_MaxIter);
}

// AddEntry
// Adds new entry and it's associated sequence
teBSFrsltCodes
CSfxArray::AddEntry(char *pszSeqIdent,	// sequence identifier, typically chromosome name
				etSeqBase   *pSeq,			// sequence to add to suffix array
				uint32_t SeqLen,				// sequence length excluding any eBaseEOS
				uint16_t Flags)				// any user specified flags
{
uint32_t Idx;
etSeqBase  *pTmpSeq;
tsSfxEntry *pTmpEntry;


if(SeqLen < 1 || SeqLen > cMaxAllowSeqLen)
	{
	AddErrMsg("CSfxArray::AddEntry","SeqLen %d not in range 1..%u",SeqLen,cMaxAllowSeqLen);
	Reset(false);
	return(eBSFerrParams);
	}

if(m_pEntriesBlock == NULL)					// will be NULL until at least one entry has been added
	{
	if(m_pBisulfateBases != NULL)
		{
#ifdef _WIN32
		free(m_pBisulfateBases);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
		if(m_pBisulfateBases != MAP_FAILED)
			munmap(m_pBisulfateBases,m_AllocBisulfiteMem);
#endif
		m_pBisulfateBases = NULL;
		}
	m_AllocBisulfiteMem = 0;

	if(m_pSfxBlock != NULL)
		{
#ifdef _WIN32
		free(m_pSfxBlock);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
		if(m_pSfxBlock != MAP_FAILED)
			munmap(m_pSfxBlock,m_AllocSfxBlockMem);
#endif
		m_pSfxBlock = NULL;
		}
	m_AllocSfxBlockMem = 0;

	// initially allocate for default number of entries
	// will realloc as required and, later when writing entries to file, will trim back to actual
	m_SfxHeader.EntriesSize = sizeof(tsSfxEntriesBlock) + (sizeof(tsSfxEntry) * (cAllocSfxEntries-1));
	m_SfxHeader.EntriesOfs = 0;
#ifdef _WIN32
	m_pEntriesBlock = (tsSfxEntriesBlock *) malloc(m_SfxHeader.EntriesSize);
	if(m_pEntriesBlock==NULL)
		{
		AddErrMsg("CSfxArray::AddEntry","unable to allocate %u bytes for holding entries block",m_SfxHeader.EntriesSize);
		Reset(false);			// closes opened file..
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and seems to have issues if more than 2GB allocation
	m_pEntriesBlock = (tsSfxEntriesBlock *)mmap(NULL,m_SfxHeader.EntriesSize, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pEntriesBlock == MAP_FAILED)
		{
		AddErrMsg("CSfxArray::AddEntry","unable to allocate %u bytes for holding entries block",m_SfxHeader.EntriesSize);
		m_pEntriesBlock = NULL;
		Reset(false);			// closes opened file..
		return(eBSFerrMem);
		}
#endif
	memset(m_pEntriesBlock,0,sizeof(tsSfxEntriesBlock));
	m_pEntriesBlock->MaxEntries = cAllocSfxEntries;
	m_AllocEntriesBlockMem = m_SfxHeader.EntriesSize;
	m_SfxHeader.NumSfxBlocks = 0;

	if(m_EstSfxEls < (int64_t)cReallocBlockEls)
		m_EstSfxEls = (int64_t)cReallocBlockEls;

	m_AllocSfxBlockMem = (uint64_t)sizeof(tsSfxBlock) + (size_t)max(m_EstSfxEls,(uint64_t)SeqLen + 10);
#ifdef _WIN32
	m_pSfxBlock = (tsSfxBlock *) malloc((size_t)m_AllocSfxBlockMem);
	if(m_pSfxBlock == NULL)
		{
		AddErrMsg("CSfxArray::AddEntry","Fatal: unable to allocate suffix block memory %zu",m_AllocSfxBlockMem);
		m_AllocSfxBlockMem = 0;
		Reset();
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and seems to have issues if more than 2GB allocation
	m_pSfxBlock = (tsSfxBlock *)mmap(NULL,m_AllocSfxBlockMem, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pSfxBlock == MAP_FAILED)
		{
		AddErrMsg("CSfxArray::AddEntry","Fatal: unable to allocate suffix block memory %zu",m_AllocSfxBlockMem);
		m_pSfxBlock = NULL;
		m_AllocSfxBlockMem = 0;
		Reset(false);			// closes opened file..
		return(eBSFerrMem);
		}
#endif
	memset(m_pSfxBlock,0,sizeof(tsSfxBlock));
	m_pSfxBlock->BlockID = 1;

	if(m_bBisulfite)
		{
		m_AllocBisulfiteMem = sizeof(tsSfxBlock) + max(cReallocBlockEls,(uint64_t)SeqLen + 10);
#ifdef _WIN32
		m_pBisulfateBases = (uint8_t *) malloc((size_t)m_AllocBisulfiteMem);
		if(m_pBisulfateBases == NULL)
			{
			AddErrMsg("CSfxArray::AddEntry","Fatal: unable to allocate Bisulfite memory %zu",m_AllocBisulfiteMem);
			m_AllocBisulfiteMem = 0;
			Reset();
			return(eBSFerrMem);
			}
#else
		// gnu malloc is still in the 32bit world and seems to have issues if more than 2GB allocation
		m_pBisulfateBases = (uint8_t *)mmap(NULL,(size_t)m_AllocBisulfiteMem, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
		if(m_pBisulfateBases == MAP_FAILED)
			{
			AddErrMsg("CSfxArray::AddEntry","Fatal: unable to allocate Bislfite memory %zu",m_AllocBisulfiteMem);
			m_pBisulfateBases = NULL;
			m_AllocBisulfiteMem = 0;
			Reset(false);			// closes opened file..
			return(eBSFerrMem);
			}
#endif
		}
	}
else	// else already at least one entry
	{
	// here is where to check if this new entry would cause the current block to exceed the max allowed block size
	if((m_pSfxBlock->ConcatSeqLen + SeqLen) > m_MaxSfxBlockEls)
		{
		AddErrMsg("CSfxArray::AddEntry","Fatal: Total concatenated sequence length (%zu) is more than maximum (%zu) supported",m_pSfxBlock->ConcatSeqLen + SeqLen,m_MaxSfxBlockEls);
		Reset(false);			// closes opened file..
		return(eBSFerrMem);
		}

	// ensure not about to exceed max number of allowed entries
	if(m_pEntriesBlock->NumEntries >= cMaxSfxEntries)
		{
		AddErrMsg("CSfxArray::AddEntry","Reached max number (%d) of supported entries",cMaxSfxEntries);
		Reset(false);			// closes opened file..
		return(eBSFerrMaxEntries);
		}

	// check if the entries block needs to be extended
	if(m_pEntriesBlock->NumEntries == m_pEntriesBlock->MaxEntries)
		{
		tsSfxEntriesBlock *pRealloc;
		size_t ReallocSize = sizeof(tsSfxEntriesBlock) + (sizeof(tsSfxEntry) * (m_pEntriesBlock->MaxEntries + cAllocSfxEntries-1));
#ifdef _WIN32
		pRealloc = (tsSfxEntriesBlock *)realloc(m_pEntriesBlock,ReallocSize);
#else
		pRealloc = (tsSfxEntriesBlock *)mremap(m_pEntriesBlock,m_AllocEntriesBlockMem,ReallocSize,MREMAP_MAYMOVE);
		if(pRealloc == MAP_FAILED)
			pRealloc = NULL;
#endif
		if(pRealloc == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Memory re-allocation to %d bytes - %s",ReallocSize,strerror(errno));
			return(eBSFerrMem);
			}

		m_pEntriesBlock = pRealloc;
		m_pEntriesBlock->MaxEntries += cAllocSfxEntries;
		m_SfxHeader.EntriesSize = (uint32_t)ReallocSize;
		m_AllocEntriesBlockMem = ReallocSize;
		}

	// check if the m_pSfxBlock needs to be extended
	if((m_pSfxBlock->ConcatSeqLen + SeqLen + 10) > m_AllocSfxBlockMem) // 10 is to allow for appended eBaseEOS's and slight safety margin
		{
		tsSfxBlock *pRealloc;
		int64_t ReallocSize = m_AllocSfxBlockMem + max(cReallocBlockEls,((uint64_t)SeqLen) + 10);	

#ifdef _WIN32
		pRealloc = (tsSfxBlock *)realloc(m_pSfxBlock,(size_t)ReallocSize);
#else
		pRealloc = (tsSfxBlock *)mremap(m_pSfxBlock,(size_t)m_AllocSfxBlockMem,(size_t)ReallocSize,MREMAP_MAYMOVE);
		if(pRealloc == MAP_FAILED)
			pRealloc = NULL;
#endif
		if(pRealloc == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: SfxBlock memory re-allocation to %zd bytes - %s",(int64_t)ReallocSize,strerror(errno));
			return(eBSFerrMem);
			}

		m_pSfxBlock = pRealloc;
		m_AllocSfxBlockMem = ReallocSize;
		}

	if(m_bBisulfite)
		{
			// check if the m_pSfxBlock needs to be extended
		if((m_pSfxBlock->ConcatSeqLen + SeqLen + 10) > m_AllocBisulfiteMem) // 10 is to allow a slight safety margin
			{
			uint8_t *pRealloc;
			int64_t ReallocSize = m_AllocBisulfiteMem + max(cReallocBlockEls,(uint64_t)SeqLen + 10);
#ifdef _WIN32
			pRealloc = (uint8_t *)realloc(m_pBisulfateBases,(size_t)ReallocSize);
#else
			pRealloc = (uint8_t *)mremap(m_pBisulfateBases,(size_t)m_AllocBisulfiteMem,(size_t)ReallocSize,MREMAP_MAYMOVE);
			if(pRealloc == MAP_FAILED)
				pRealloc = NULL;
#endif
			if(pRealloc == NULL)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Bisulfite memory re-allocation to %zd bytes - %s",ReallocSize,strerror(errno));
				return(eBSFerrMem);
				}

			m_pBisulfateBases = pRealloc;
			m_AllocBisulfiteMem = ReallocSize;
			}
		}

	}

pTmpEntry = &m_pEntriesBlock->Entries[m_pEntriesBlock->NumEntries++];
pTmpEntry->EntryID = m_pEntriesBlock->NumEntries;
strncpy((char *)pTmpEntry->szSeqName,pszSeqIdent,sizeof(pTmpEntry->szSeqName));
pTmpEntry->szSeqName[sizeof(pTmpEntry->szSeqName)-1] = '\0';
pTmpEntry->NameHash = CUtility::GenHash16(pszSeqIdent);

pTmpEntry->fBlockID = (m_pSfxBlock->BlockID & 0x0ff) | Flags << 8;
pTmpEntry->SeqLen = SeqLen;
pTmpEntry->StartOfs = m_pSfxBlock->ConcatSeqLen;
pTmpEntry->EndOfs = m_pSfxBlock->ConcatSeqLen + SeqLen - 1;
pTmpSeq = &m_pSfxBlock->SeqSuffix[m_pSfxBlock->ConcatSeqLen];
for(Idx = 0; Idx < SeqLen; Idx++,pTmpSeq++,pSeq++)
	*pTmpSeq = *pSeq & ~cRptMskFlg;
*pTmpSeq = eBaseEOS;
m_pSfxBlock->NumEntries+=1;
m_pSfxBlock->ConcatSeqLen += SeqLen + 1;

return(eBSFSuccess);
}


// finalise and commit to disk after all entries added with AddEntry()
int
CSfxArray::Finalise(void)
{
int Rslt;
if ((Rslt = Flush2Disk()) == eBSFSuccess)
	{
	if (!m_bInMemSfx)
		Reset(false);
	}
return(Rslt);
}


// ThreadedPrereadBlocks
// Background thread waits for job requests which could be either to load a specified m_ReqBlockID from disk or to self terminate
// Job requests are signaled through m_JobReqEvent, and job completion notified back to the requestor through m_JobAckEvent
#ifdef _WIN32
unsigned __stdcall CSfxArray::ThreadedPrereadBlocks(void * pThreadPars)
#else
void *CSfxArray::ThreadedPrereadBlocks(void * pThreadPars)
#endif
{
CSfxArray *pSfxArray = (CSfxArray *)pThreadPars;
teBSFrsltCodes Rslt;
bool bTermThread;
int64_t BlockSize;
int ReqBlockID;
bTermThread = false;		// set true if main thread requests this thread to terminate

do {
#ifdef _WIN32
	switch(WaitForSingleObject(pSfxArray->m_JobReqEvent,cSigWaitSecs*1000)) {
		case WAIT_TIMEOUT:			// periodically wakeup incase m_JobReqEvent signal was missed
		case WAIT_OBJECT_0:			// m_JobReqEvent was signalled
			break;
		default:					// any other indicates that object has been abandoned
			bTermThread = true;		// treat as a terminate request
			continue;
		}
	// need to serialise access to shared job resources
	switch(WaitForSingleObject(pSfxArray->m_JobMutex,cSigWaitSecs*1000)) {
		case WAIT_TIMEOUT:			// timed out, loop back
			continue;
		case WAIT_OBJECT_0:			// gained access
			break;
		default:					// any other indicates that object has been abandoned
			bTermThread = true;		// treat as a terminate request
			continue;
		}

	// terminate requests have priority
	if(bTermThread = pSfxArray->m_bTermThread)
		{
		SetEvent(pSfxArray->m_JobAckEvent);
		ReleaseMutex(pSfxArray->m_JobMutex);
		break;
		}
#else
	struct timespec abstime;
	pthread_mutex_lock(&pSfxArray->m_JobMutex);
	while(!pSfxArray->m_bTermThread && pSfxArray->m_ReqBlockID == 0)
		{
		clock_gettime(CLOCK_REALTIME,&abstime);
		abstime.tv_sec += cSigWaitSecs;
		pthread_cond_timedwait(&pSfxArray->m_JobReqEvent,&pSfxArray->m_JobMutex,&abstime);
		}
	// terminate requests have priority
	if(bTermThread = pSfxArray->m_bTermThread)
		{
		pthread_cond_signal(&pSfxArray->m_JobAckEvent);
		pthread_mutex_unlock(&pSfxArray->m_JobMutex);
		break;
		}
#endif
	if((ReqBlockID = pSfxArray->m_ReqBlockID) > 0)	// fgnd has requested block to be loaded?
		{
		if(pSfxArray->m_pSfxBlock->BlockID == ReqBlockID)	// is already loaded?
			{
			ReqBlockID = 0;
			pSfxArray->m_ReqBlockID = 0;
			pSfxArray->m_ReqBlockRslt = eBSFSuccess;
#ifdef _WIN32
			SetEvent(pSfxArray->m_JobAckEvent);
			ReleaseMutex(pSfxArray->m_JobMutex);
#else
			pthread_cond_signal(&pSfxArray->m_JobAckEvent);
			pthread_mutex_unlock(&pSfxArray->m_JobMutex);
#endif
			}

		// if request was not preloaded then need to load from disk
		if(ReqBlockID) // != 0 if need to load
			{
			int64_t BlkOfs;
			BlkOfs = 0;
			BlockSize = pSfxArray->m_SfxHeader.SfxBlockSize - BlkOfs;
			if((Rslt=pSfxArray->ChunkedRead(pSfxArray->m_SfxHeader.SfxBlockOfs,((uint8_t *)pSfxArray->m_pSfxBlock) + BlkOfs,BlockSize)) < eBSFSuccess)
				{
				pSfxArray->m_ReqBlockID = 0;
				pSfxArray->m_ReqBlockRslt = Rslt;	// error reading from disk, let main thread know
#ifdef _WIN32
				SetEvent(pSfxArray->m_JobAckEvent);
				ReleaseMutex(pSfxArray->m_JobMutex);
#else
				pthread_cond_signal(&pSfxArray->m_JobAckEvent);
				pthread_mutex_unlock(&pSfxArray->m_JobMutex);
#endif
				continue;
				}
			if(Rslt == 1)		// chunked read may have exited early because m_bTermThread was set true
				{
				pSfxArray->m_ReqBlockID = 0;
				bTermThread = true;
#ifdef _WIN32
				SetEvent(pSfxArray->m_JobAckEvent);
				ReleaseMutex(pSfxArray->m_JobMutex);
#else
				pthread_cond_signal(&pSfxArray->m_JobAckEvent);
				pthread_mutex_unlock(&pSfxArray->m_JobMutex);
#endif
				break;
				}

			ReqBlockID = 0;		// requested block now loaded
			pSfxArray->m_ReqBlockID = 0;
			pSfxArray->m_ReqBlockRslt = eBSFSuccess;
#ifdef _WIN32
			SetEvent(pSfxArray->m_JobAckEvent);
			ReleaseMutex(pSfxArray->m_JobMutex);
#else
			pthread_cond_signal(&pSfxArray->m_JobAckEvent);
			pthread_mutex_unlock(&pSfxArray->m_JobMutex);
#endif
			}
		}
	else
		{
#ifdef _WIN32
		ReleaseMutex(pSfxArray->m_JobMutex);
#else
		pthread_mutex_unlock(&pSfxArray->m_JobMutex);
#endif
		}
	}
while(!bTermThread);
// have been requested to terminate this thread
pSfxArray->m_ReqBlockRslt = eBSFSuccess;
#ifdef _WIN32
_endthreadex(0);
return(eBSFSuccess); // unreached, but keeps compilers happy!
#else
pthread_exit(&pSfxArray->m_ReqBlockRslt);
#endif
}


// Disk2SfxBlock
teBSFrsltCodes
CSfxArray::Disk2SfxBlock(int BlockID)
{
teBSFrsltCodes Rslt;


if(m_pSfxBlock == NULL || !m_bThreadActive)
	return(eBSFerrInternal);

if(BlockID < 1 || m_SfxHeader.NumSfxBlocks == 0 || (uint32_t)BlockID > m_SfxHeader.NumSfxBlocks)
	return(eBSFerrParams);


do {
#ifdef _WIN32
	// gain lock, could take some time if background thread loading large SfxBlock from disk
	switch(WaitForSingleObject(m_JobMutex,cSigMaxWaitSecs * 1000)) {
		case WAIT_OBJECT_0:		// gained the mutex
			break;
		case WAIT_TIMEOUT:		// timed out, treat as a file access error
			return(eBSFerrFileAccess);
		default:				// error..
			return(cBSFSyncObjErr);
		}
#else
	pthread_mutex_lock(&m_JobMutex);
#endif
	// already loaded?
	if(m_pSfxBlock->BlockID == BlockID)
		break;

	// block not loaded, or is in readahead
	// request preload thread to make available the required block
	m_ReqBlockID = BlockID;
#ifdef _WIN32
	SetEvent(m_JobReqEvent);
	ReleaseMutex(m_JobMutex);
	WaitForSingleObject(m_JobAckEvent,cSigWaitSecs * 1000);
#else
	struct timespec abstime;
	clock_gettime(CLOCK_REALTIME,&abstime);
	abstime.tv_sec += cSigWaitSecs;
	pthread_cond_signal(&m_JobReqEvent);
	pthread_cond_timedwait(&m_JobAckEvent,&m_JobMutex,&abstime);
	pthread_mutex_unlock(&m_JobMutex);
#endif
	}
while(1);
Rslt = m_ReqBlockRslt;
#ifdef _WIN32
ReleaseMutex(m_JobMutex);
#else
pthread_mutex_unlock(&m_JobMutex);
#endif
return(Rslt);
}

int
CSfxArray::Next(int PrevBlockID)
{
if(PrevBlockID < 0 || (uint32_t)PrevBlockID >= m_SfxHeader.NumSfxBlocks)
	return(0);
return(PrevBlockID + 1);
}

// SetTargBlock
// specifies which block is to be loaded for subsequent sequence matches
int
CSfxArray::SetTargBlock(int BlockID)
{
return(Disk2SfxBlock(BlockID));
}

int										// if non-zero then returned number of identifiers 
CSfxArray::ChkDupEntries(int MaxIdents,		// maximum number of identifers to return in pIdents (caller allocates to hold returned identifiers)
					  uint32_t *pIdents)		// checks if there are duplicate entry names and reports identifier
{
int NumDups;
bool bInDup;
tsSfxEntry **ppEntries;
tsSfxEntry *pCurEntry;
tsSfxEntry *pPrevEntry;

uint32_t Idx;

if(m_pEntriesBlock == NULL || MaxIdents > 0 && pIdents == NULL) 
	return(eBSFerrEntry);
if(m_pEntriesBlock->NumEntries < 2)
	return(0);

if((ppEntries = new tsSfxEntry * [m_pEntriesBlock->NumEntries]) == NULL)
	{
	AddErrMsg("CSfxArray::ChkDupEntries","unable to allocate memory for sorting %d entry pointers",m_pEntriesBlock->NumEntries);
	return(eBSFerrMem);
	}

for(Idx = 0; Idx < m_pEntriesBlock->NumEntries; Idx++)
	ppEntries[Idx] = &m_pEntriesBlock->Entries[Idx];

// now sort on entries names
m_MTqsort.qsort(ppEntries,m_pEntriesBlock->NumEntries,sizeof(tsSfxEntry *),QSortEntryNames);

// iterate looking for duplicates
NumDups = 0;
bInDup = false;
pPrevEntry = ppEntries[0];
for(Idx = 1; Idx < m_pEntriesBlock->NumEntries && NumDups < MaxIdents; Idx++)
	{
	pCurEntry = ppEntries[Idx];
	if(stricmp((char *)pCurEntry->szSeqName,(char *)pPrevEntry->szSeqName))
		{
		pPrevEntry = pCurEntry;
		bInDup = false;
		continue;
		}
	if(!bInDup && MaxIdents && NumDups < MaxIdents)
		{
		*pIdents++ = pPrevEntry->EntryID;
		NumDups += 1;
		bInDup = true;
		}
	else
		if(!MaxIdents)
			{
			NumDups = 1;
			break;
			}
	}
delete []ppEntries;
return(NumDups);
}


uint16_t											// returned flags prior to this call
CSfxArray::SetResetIdentFlags(uint32_t EntryID,	// atomic set/reset flags for this entry
					uint16_t SetFlags,			// set these flags then
					uint16_t ResetFlags)			// reset these flags
{
uint32_t fBlockID;
uint16_t Flags;
if(SetFlags == 0 && ResetFlags == 0)
	return(GetIdentFlags(EntryID));

if(m_pEntriesBlock == NULL || EntryID < 1 || EntryID > m_pEntriesBlock->NumEntries)
	return(eBSFerrEntry);
SerialiseBaseFlags();
fBlockID = m_pEntriesBlock->Entries[EntryID-1].fBlockID;
Flags = (fBlockID >> 8) & 0x0ffff;
Flags |= SetFlags;
Flags &= ~ResetFlags;
m_pEntriesBlock->Entries[EntryID-1].fBlockID = ((uint32_t)Flags << 8) | fBlockID & 0x0ff;
ReleaseBaseFlags();
return((fBlockID >> 8) & 0x0ffff);
}

uint16_t					// returned flags
CSfxArray::GetIdentFlags(uint32_t EntryID)
{
uint32_t Flags;
if(m_pEntriesBlock == NULL || EntryID < 1 || EntryID > m_pEntriesBlock->NumEntries)
	return(eBSFerrEntry);
SerialiseBaseFlags();
Flags = m_pEntriesBlock->Entries[EntryID-1].fBlockID;
ReleaseBaseFlags();
return((Flags >> 8) & 0x0ffff);
}


void
CSfxArray::InitAllIdentFlags(uint16_t Flags)		// initialise all entries to have these flags
{
uint32_t EntryID;
uint32_t fBlockID;
uint32_t SetFlags;
tsSfxEntry *pSfxEntry;
if(m_pEntriesBlock == NULL || m_pEntriesBlock->NumEntries < 1)
	return;
SetFlags = (uint32_t)Flags << 8;
pSfxEntry = m_pEntriesBlock->Entries;
SerialiseBaseFlags();
for(EntryID = 0; EntryID < m_pEntriesBlock->NumEntries; EntryID++,pSfxEntry++)
	{
	fBlockID = pSfxEntry->fBlockID & 0x0ff;
	fBlockID |= SetFlags;
	pSfxEntry->fBlockID = fBlockID;
	}
ReleaseBaseFlags();
return;
}


int
CSfxArray::GetIdentName(uint32_t EntryID,int MaxLen,char *pszSeqIdent)
{
if(m_pEntriesBlock == NULL || EntryID < 1 || (uint32_t)EntryID > m_pEntriesBlock->NumEntries)
	return(eBSFerrEntry);
strncpy(pszSeqIdent,(const char *)m_pEntriesBlock->Entries[EntryID-1].szSeqName,MaxLen);
pszSeqIdent[MaxLen-1] = '\0';
return(eBSFSuccess);
}

int		// returns identifer for specified sequence name
CSfxArray::GetIdent(char *pszSeqIdent)
{
uint32_t Idx;
if(m_pEntriesBlock == NULL || m_pEntriesBlock->NumEntries < 1)
	return(eBSFerrEntry);
for(Idx = 0; Idx < m_pEntriesBlock->NumEntries; Idx++)
	{
	if(!stricmp(pszSeqIdent,(const char *)m_pEntriesBlock->Entries[Idx].szSeqName))
		return((int)Idx+1);
	}
return(eBSFerrEntry);
}

int
CSfxArray::GetNumEntries(void)
{
if(m_pEntriesBlock == NULL)
	return(eBSFerrEntry);
return(m_pEntriesBlock->NumEntries);
}

uint32_t
CSfxArray::GetSeqLen(uint32_t EntryID)
{
tsSfxEntry *pEntry;
if(m_pEntriesBlock == NULL || EntryID < 1 || (uint32_t)EntryID > m_pEntriesBlock->NumEntries)
	return(0);
pEntry = &m_pEntriesBlock->Entries[EntryID-1];
return(pEntry->SeqLen);
}

// returns total length of all sequences
uint64_t
CSfxArray::GetTotSeqsLen(void)
{
uint64_t Idx;
uint64_t TotLen;
tsSfxEntry *pEntry;
if(m_pEntriesBlock == NULL || !m_pEntriesBlock->NumEntries)
	return((size_t)eBSFerrEntry);
TotLen = 0;
pEntry = &m_pEntriesBlock->Entries[0];
for(Idx = 0; Idx < m_pEntriesBlock->NumEntries; Idx++,pEntry++)
	TotLen += pEntry->SeqLen;
return(TotLen);
}

// returns minimum length of any sequence in currently loaded entries block
uint32_t 
CSfxArray::GetMinSeqLen(void)							
{
uint64_t Idx;
uint32_t MinLen;
tsSfxEntry *pEntry;
if(m_pEntriesBlock == NULL || !m_pEntriesBlock->NumEntries)
	return(eBSFerrEntry);
MinLen = 0;
pEntry = &m_pEntriesBlock->Entries[0];
for(Idx = 0; Idx < m_pEntriesBlock->NumEntries; Idx++,pEntry++)
	if(MinLen == 0 || pEntry->SeqLen < MinLen)
		MinLen = pEntry->SeqLen;
return(MinLen);
}

// returns maximum length of any sequence in currently loaded entries block
uint32_t 
CSfxArray::GetMaxSeqLen(void)							
{
uint64_t Idx;
uint32_t MaxLen;
tsSfxEntry *pEntry;
if(m_pEntriesBlock == NULL || !m_pEntriesBlock->NumEntries)
	return(eBSFerrEntry);
MaxLen = 0;
pEntry = &m_pEntriesBlock->Entries[0];
for(Idx = 0; Idx < m_pEntriesBlock->NumEntries; Idx++,pEntry++)
	if(pEntry->SeqLen > MaxLen)
		MaxLen = pEntry->SeqLen;
return(MaxLen);
}

void
CSfxArray::SerialiseBaseFlags(void)
{
int SpinCnt = 500;
int BackoffMS = 5;

#ifdef _WIN32
while(InterlockedCompareExchange(&m_CASSeqFlags,1,0)!=0)
	{
	if(SpinCnt -= 1)
		continue;
	CUtility::SleepMillisecs(BackoffMS);
	SpinCnt = 100;
	if(BackoffMS < 500)
		BackoffMS += 2;
	}
#else
while(__sync_val_compare_and_swap(&m_CASSeqFlags,0,1)!=0)
	{
	if(SpinCnt -= 1)
		continue;
	CUtility::SleepMillisecs(BackoffMS);
	SpinCnt = 100;
	if(BackoffMS < 500)
		BackoffMS += 2;
	}
#endif
}

void
CSfxArray::ReleaseBaseFlags(void)
{
#ifdef _WIN32
InterlockedCompareExchange(&m_CASSeqFlags,0,1);
#else
__sync_val_compare_and_swap(&m_CASSeqFlags,1,0);
#endif
}

int 
CSfxArray::GetBaseFlags(uint32_t EntryID,		// identifies sequence containing loci flags to be returned - flags returned are in bits 0..3
				uint32_t Loci)				// offset within sequence of base flags to return
{
int Rslt;
int Flags;
tsSfxEntry *pEntry;
etSeqBase *pSeq;
if(m_pEntriesBlock == NULL || EntryID < 1 || (uint32_t)EntryID > m_pEntriesBlock->NumEntries || m_bColorspace)
	return(eBSFerrEntry);
pEntry = &m_pEntriesBlock->Entries[EntryID-1];

if(Loci >= pEntry->SeqLen)	// requested start offset must be less or equal than the end of sequence
	return(eBSFerrInternal);

// ensure SfxBlock containing entry sequence is loaded into memory
if(m_pSfxBlock == NULL || m_pSfxBlock->BlockID != (pEntry->fBlockID & 0x0ff))
	if((Rslt=Disk2SfxBlock(pEntry->fBlockID  & 0x0ff))!=eBSFSuccess)
		return(Rslt);

pSeq = &m_pSfxBlock->SeqSuffix[pEntry->StartOfs];
SerialiseBaseFlags();
Flags = pSeq[Loci];
ReleaseBaseFlags();
return((Flags >> 4) & 0x0f);
}

int									// NOTE: returns original flags (EntryID flags in bits 3..0) immediately prior to the set set/reset
CSfxArray::SetBaseFlags(uint32_t EntryID,	// identifies sequence containing loci flags to be set
				uint32_t Loci,				// offset within sequence of flags to set
	            int SetFlgs,			// set flags in bits 0..3
				int ResetFlgs)			// reset flags in bits 0..3
{
int Rslt;
int Flags;
int PriorFlags;
tsSfxEntry *pEntry;
etSeqBase *pSeq;
if(m_pEntriesBlock == NULL || EntryID < 1 || (uint32_t)EntryID > m_pEntriesBlock->NumEntries || m_bColorspace)
	return(eBSFerrEntry);
pEntry = &m_pEntriesBlock->Entries[EntryID-1];

if(Loci >= pEntry->SeqLen)	// requested start offset must be less or equal than the end of sequence
	return(eBSFerrInternal);

// ensure SfxBlock containing entry sequence is loaded into memory
if(m_pSfxBlock == NULL || m_pSfxBlock->BlockID != (pEntry->fBlockID & 0x0ff))
	if((Rslt=Disk2SfxBlock((pEntry->fBlockID & 0x0ff)))!=eBSFSuccess)
		return(Rslt);

pSeq = &m_pSfxBlock->SeqSuffix[pEntry->StartOfs];

SetFlgs <<= 4;
SetFlgs &= 0x0f0;
ResetFlgs <<= 4;
ResetFlgs &= 0x0f0;

SerialiseBaseFlags();
PriorFlags = Flags = pSeq[Loci];
Flags &= ~ResetFlgs;
Flags |= SetFlgs;
pSeq[Loci] = Flags;
ReleaseBaseFlags();

return((PriorFlags >> 4) & 0x0f);
}

int									// NOTE: returns original flags (EntryID1 flags in bits 3..0 and EntryID2 in bits 7..4) immediately prior to the set set/reset
CSfxArray::SetBaseFlags(int EntryID1,	// identifies 1st sequence containing loci flags to be set
				uint32_t Loci1,			// offset within sequence of flags to set
				int EntryID2,			// identifies 2nd sequence containing loci flags to be set
				uint32_t Loci2,			// offset within sequence of flags to set
	            int SetFlgs,			// set flags in bits 0..3
				int ResetFlgs)			// reset flags in bits 0..3
{
int Rslt;
int Flags;
int PriorFlags1;
int PriorFlags2;
tsSfxEntry *pEntry1;
tsSfxEntry *pEntry2;
etSeqBase *pSeq1;
etSeqBase *pSeq2;

if(m_pEntriesBlock == NULL || EntryID1 < 1 || (uint32_t)EntryID1 > m_pEntriesBlock->NumEntries || EntryID2 < 1 || (uint32_t)EntryID2 > m_pEntriesBlock->NumEntries|| m_bColorspace)
	return(eBSFerrEntry);


pEntry1 = &m_pEntriesBlock->Entries[EntryID1-1];
if(Loci1 >= pEntry1->SeqLen)	// requested start offset must be less or equal than the end of sequence
	return(eBSFerrInternal);

pEntry2 = &m_pEntriesBlock->Entries[EntryID2-1];
if(Loci2 >= pEntry2->SeqLen)	// requested start offset must be less or equal than the end of sequence
	return(eBSFerrInternal);


// ensure SfxBlock containing entry sequence is loaded into memory
// currently only a single SfxBlock is supported so only need to check one of the entries
if(m_pSfxBlock == NULL || m_pSfxBlock->BlockID != (pEntry1->fBlockID & 0x0ff))
	if((Rslt=Disk2SfxBlock(pEntry1->fBlockID & 0x0ff))!=eBSFSuccess)
		return(Rslt);

SetFlgs <<= 4;
SetFlgs &= 0x0f0;
ResetFlgs <<= 4;
ResetFlgs &= 0x0f0;
pSeq1 = &m_pSfxBlock->SeqSuffix[pEntry1->StartOfs];
pSeq2 = &m_pSfxBlock->SeqSuffix[pEntry2->StartOfs];

SerialiseBaseFlags();
PriorFlags1 = Flags = pSeq1[Loci1];
Flags &= ~ResetFlgs;
Flags |= SetFlgs;
pSeq1[Loci1] = Flags;

PriorFlags2 = Flags = pSeq2[Loci2];
Flags &= ~ResetFlgs;
Flags |= SetFlgs;
pSeq2[Loci2] = Flags;
ReleaseBaseFlags();

PriorFlags1 = (PriorFlags1 >> 4) & 0x0f;
PriorFlags1 |= (PriorFlags2 & 0x0f0);
return(PriorFlags1);
}

// GetBase
// Returns base for EntryID at Loci
// Base is returned in bits 0..3
int
CSfxArray::GetBase(int EntryID,			// identifies sequence containing base to be returned
				uint32_t Loci)					// offset within sequence of base to return
{
int Rslt;
tsSfxEntry *pEntry;
etSeqBase *pSeq;
if(m_pEntriesBlock == NULL || EntryID < 1 || (uint32_t)EntryID > m_pEntriesBlock->NumEntries)
	return(eBSFerrEntry);
pEntry = &m_pEntriesBlock->Entries[EntryID-1];

if(Loci >= pEntry->SeqLen)	// requested start offset must be less or equal than the end of sequence
	return(eBSFerrInternal);

// ensure SfxBlock containing entry sequence is loaded into memory
if(m_pSfxBlock == NULL || m_pSfxBlock->BlockID != (pEntry->fBlockID & 0x0ff))
	if((Rslt=Disk2SfxBlock(pEntry->fBlockID & 0x0ff))!=eBSFSuccess)
		return(Rslt);

pSeq = &m_pSfxBlock->SeqSuffix[pEntry->StartOfs];
if(m_bColorspace)
	return((pSeq[Loci] >> 4) & 0x0f);
return(pSeq[Loci] & 0x0f);
}


// GetSeq
// Returns a copy of the requested subsequence into pRetSeq from the
// sequence at EntryID starting at offset Loci and of maximal length Len
uint32_t					// returned sequence length (may be shorter than that requested) or 0 if errors
CSfxArray::GetSeq(int EntryID,uint32_t Loci,etSeqBase *pRetSeq,uint32_t Len)
{
int Rslt;
uint32_t RetLen;
tsSfxEntry *pEntry;
etSeqBase *pSeq;
if(m_pEntriesBlock == NULL || EntryID < 1 || (uint32_t)EntryID > m_pEntriesBlock->NumEntries)
	return(0);
pEntry = &m_pEntriesBlock->Entries[EntryID-1];

if(Loci >= pEntry->SeqLen)	// any sequence to return?
	return(0);

// ensure SfxBlock containing entry sequence is loaded into memory
if(m_pSfxBlock == NULL || m_pSfxBlock->BlockID != (pEntry->fBlockID & 0x0ff))
	if((Rslt=Disk2SfxBlock(pEntry->fBlockID & 0x0ff))!=eBSFSuccess)
		return(0);

pSeq = &m_pSfxBlock->SeqSuffix[pEntry->StartOfs];
pSeq = &pSeq[Loci];
Len = RetLen = min(Len,pEntry->SeqLen - Loci);
if(m_bColorspace)
	while(Len--)
		*pRetSeq++ = (*pSeq++ >> 4) & 0x0f;
else
	while(Len--)
		*pRetSeq++ = *pSeq++ & 0x0f;
return(RetLen);
}

// GetPtrSeq
// Returns a ptr to the sequence at EntryID starting at offset Loci
etSeqBase *
CSfxArray::GetPtrSeq(int EntryID,uint32_t Loci)
{
int Rslt;
tsSfxEntry *pEntry;
etSeqBase *pSeq;
if(m_pEntriesBlock == NULL || EntryID < 1 || (uint32_t)EntryID > m_pEntriesBlock->NumEntries)
	return(NULL);
pEntry = &m_pEntriesBlock->Entries[EntryID-1];
if(Loci >= pEntry->SeqLen)	// requested start offset must be less or equal than the end of sequence
	return(NULL);

// ensure SfxBlock containing entry sequence is loaded into memory
if(m_pSfxBlock == NULL || m_pSfxBlock->BlockID != (pEntry->fBlockID & 0x0ff))
	if((Rslt=Disk2SfxBlock(pEntry->fBlockID & 0x0ff))!=eBSFSuccess)
		return(NULL);
pSeq = &m_pSfxBlock->SeqSuffix[pEntry->StartOfs];
return(&pSeq[Loci]);
}


// GetColorspaceSeq
// Returns a copy of the requested subsequence in color space into pRetSeq from the
// sequence at EntryID starting at offset Loci and of length Len
uint32_t			// returned sequence length (may be shorter than that requested) or 0 if errors
CSfxArray::GetColorspaceSeq(int EntryID,uint32_t Loci,etSeqBase *pRetSeq,uint32_t Len)
{
uint32_t Rslt;
uint32_t RetLen;
tsSfxEntry *pEntry;
etSeqBase *pSeq;
if(m_pEntriesBlock == NULL || EntryID < 1 || (uint32_t)EntryID > m_pEntriesBlock->NumEntries || !m_bColorspace)
	return(0);
pEntry = &m_pEntriesBlock->Entries[EntryID-1];
if(Loci >= pEntry->SeqLen)	// requested start offset must be less or equal than the end of sequence
	return(0);

// ensure SfxBlock containing entry sequence is loaded into memory
if(m_pSfxBlock == NULL || m_pSfxBlock->BlockID != (pEntry->fBlockID & 0xff))
	if((Rslt=Disk2SfxBlock(pEntry->fBlockID & 0x0ff))!=eBSFSuccess)
		return(0);
pSeq = &m_pSfxBlock->SeqSuffix[pEntry->StartOfs];
RetLen = Len = min(Len,pEntry->SeqLen - Loci);
pSeq = &pSeq[Loci];
while(Len--)
	*pRetSeq++ = *pSeq++ & 0x0f;
return(RetLen);
}

// ExactMatchLen
// Returns the exact match length, up to MaxMatchLen, 3' to pProbe and pTarg 
// Extends  probe against target taking into account the eBaseEOS terminator
int													// length of exact match
CSfxArray::ExactMatchLen(etSeqBase *pProbe,		// determine exactly matching length between probe
							etSeqBase *pTarg,		// and target sequences
							int MaxMatchLen)		// for up to this length
{
int MatchLen;
uint8_t El1;
uint8_t El2;
if(pProbe == NULL || pTarg == NULL || MaxMatchLen == 0)
	return(0);
for(MatchLen=0; MatchLen < MaxMatchLen; MatchLen++)
	{
	El2 = *pTarg++ & 0x0f;
	if(El2 > eBaseT)
		break;
	El1 = *pProbe++ & 0x0f;
	if(El1 > eBaseT)
		break;
	if(El1 != El2)
		break;
	}
return(MatchLen);
}


// CmpProbeTarg
// Compares probe against target taking into account the eBaseEOS terminator
int
CSfxArray::CmpProbeTarg(etSeqBase *pEl1,etSeqBase *pEl2,int Len)
{
int Psn;
uint8_t El1;
uint8_t El2;
for(Psn=0; Psn < Len; Psn++)
	{
	El2 = *pEl2++ & 0x0f;
	if(El2 == eBaseEOS)
		return(-1);
	El1 = *pEl1++ & 0x0f;
	if(El1 > El2)
		return(1);
	if(El1 < El2)
		return(-1);
	}
return(0);
}


// BSCmpProbeTarg
// Compares probe against target taking into account the eBaseEOS terminator
// For bisulfite processing treats either pEl1(C or T) as matching pEl2(C or T), and pEl1(A or G) as matching pEl2(A or G)
int
CSfxArray::BSCmpProbeTarg(etSeqBase *pEl1,etSeqBase *pEl2,int Len)
{
int Psn;
uint8_t El1;
uint8_t El2;
for(Psn=0; Psn < Len; Psn++)
	{
	El1 = *pEl1++ & 0x0f;
	El2 = *pEl2++ & 0x0f;
	switch(El2) {
		case eBaseA:
		case eBaseG:
			if(El1 == eBaseA || El1 == eBaseG)
				continue;
			if(El1 == eBaseC || El1 == eBaseT)
				return(-1);
			return(1);

		case eBaseC:
		case eBaseT:
			if(El1 == eBaseC || El1 == eBaseT)
				continue;
			return(1);

		case eBaseN:
			if(El1 == eBaseN)
				continue;
			return(-1);

		case eBaseEOS:
			return(-1);

		default:
			return(1);
		}
	}
return(0);
}


// GetBisBase
// Uses cnts of probe vs target mismatches at target eBaseC and eBaseG to determine
// most likely bisulfite base mapping for that probe read
etSeqBase
CSfxArray::GetBisBase(int TargMatchLen,etSeqBase *pTargBase,etSeqBase *pProbeBase)
{
int Idx;
int CntA = 0;
int CntT = 0;

for(Idx = 0; Idx < TargMatchLen; Idx++,pTargBase++,pProbeBase++)
	{
	if(*pTargBase == eBaseEOS)		// mustn't match across entry sequences
		break;
	if(*pProbeBase == *pTargBase)
		continue;

	// mismatch...
	if(*pTargBase == eBaseC || *pTargBase == eBaseG)
		{
		if(*pProbeBase == eBaseA && *pTargBase == eBaseG)
			CntA++;
		else
			if(*pProbeBase == eBaseT && *pTargBase == eBaseC)
				CntT++;
		}
	}
if(CntA >= CntT)
	return(eBaseA);
return(eBaseT);
}


// MapChunkHit2Entry
// Maps the chunk hit loci to the relevant sequence entry
// If many entries expected then this mapping would be a good candidate for optimisation!
tsSfxEntry *
CSfxArray::MapChunkHit2Entry(uint64_t ChunkOfs)
{

uint32_t CurBlockID;

CurBlockID = m_pSfxBlock->BlockID;

tsSfxEntry *pProbe;
int64_t Lo,Mid,Hi;	// search limits

Lo = 0;
Hi = m_pEntriesBlock->NumEntries-1;
while(Hi >= Lo) {
	Mid = (Hi + Lo)/2;
	pProbe = &m_pEntriesBlock->Entries[Mid];
	if((pProbe->fBlockID & 0x0ff) > CurBlockID)
		{
		Hi = Mid - 1;
		continue;
		}
	else
		if((pProbe->fBlockID & 0x0ff) < CurBlockID)
			{
			Lo = Mid + 1;
			continue;
			}

	// block matches
	if(pProbe->StartOfs <= ChunkOfs && pProbe->EndOfs >= ChunkOfs)
		return(pProbe);

	if(pProbe->StartOfs > ChunkOfs)
		{
		Hi = Mid - 1;
		continue;
		}
	else
		if(pProbe->EndOfs < ChunkOfs)
			{
			Lo = Mid + 1;
			continue;
			}
	}

return(NULL);
}

const int cMaxBisProbeLen = 8196;						// max allowed bisulfate probe sequence length

// Bisulfite entry point
int64_t						// returned match index or 0 if no matches or errors
CSfxArray::BSIterateExacts(etSeqBase *pProbeSeq,		// probe
 						 uint32_t ProbeLen,			// probe length
						 int64_t PrevHitIdx,		// 0 if starting new sequence, otherwise set to return value of previous successful iteration return
 						 uint32_t *pTargEntryID,	// if match then where to return suffix entry (chromosome) matched on
						 uint32_t *pHitLoci)		// if match then where to return loci
{
int64_t HitIdx;
int Idx;
tsSfxEntry *pEntry;
etSeqBase BisBase;
etSeqBase *pSeq;
etSeqBase *pProbe;
etSeqBase *pBisSeq;
etSeqBase BisProbe[cMaxBisProbeLen];	// used to hold copy of probe with eBaseA's mapped to eBaseG's and eBaseT's mapped to eBaseC's

if(!m_bBisulfite)
	return(-1);
if(ProbeLen > cMaxBisProbeLen)
	return(-1);

pBisSeq = BisProbe;
pSeq = pProbeSeq;
BisBase = eBaseN;

while((HitIdx = IterateExacts(pSeq,		// probe
 					 ProbeLen,			// probe length
					 PrevHitIdx,		// 0 if starting new sequence, otherwise set to return value of previous successful iteration return
 					 pTargEntryID,		// if match then where to return suffix entry (chromosome) matched on
					 pHitLoci))>0)		// if match then where to return loci
	{
	// check original probe against real sequence here
	pEntry = &m_pEntriesBlock->Entries[*pTargEntryID-1];
	pSeq = &m_pSfxBlock->SeqSuffix[pEntry->StartOfs + *pHitLoci];
	pProbe = pProbeSeq;
	BisBase = GetBisBase(ProbeLen,pSeq,pProbe);
	for(Idx=0; Idx < (int)ProbeLen; Idx++,pProbe++,pSeq++)
		{
		if(*pProbe == *pSeq)
			continue;
		if(*pSeq == eBaseC || *pSeq == eBaseG)
			{
			switch(BisBase) {
				case eBaseA:		// only allow A
					if(*pProbe == eBaseA && *pSeq == eBaseG)
						continue;
					break;
				case eBaseT:		// only allow T
					if(*pProbe == eBaseT && *pSeq == eBaseC)
						continue;
					break;
				}
			}
		break;
		}
	if(Idx == ProbeLen)
		return(HitIdx);
	PrevHitIdx = HitIdx;
	}
return(HitIdx);
}

// locate all sequences whereby -
//
// a) Check if starting sequence K-mer index is not the lowest index for that identical K-mer sequence
//    If not first then iterate until K-mer no longer matches, this will be the starting K-mer sequence
// b) Iterate over K-mers accumulating counts of identical K-mers attributing to the originating cultivar
//    When the K-mer sequence changes then pass counts back via the callback function
// Iterate all sequences starting at the specified suffix element index and finishing at StartSfxIdx + SfxIdxRange - 1
//

// GenKMerCultThreadRange
// For each thread instance and given starting suffix index determine the ending suffix index that the thread should process
// so as to spread the processing reasonably uniformally between the threads  
int													// < 0 if errors, 0 if no processing required by thread, 1 if *pSfxIdxRange initilised to number of elements to process 
CSfxArray::GenKMerCultThreadRange(int KMerLen,	// processing K-Mers of this length
					  int ThreadInst,				// starting and range required for this thread instance
					  int NumThreads,				// total number of threads which will be used for generating cultivar K-Mers
					  int64_t StartSfxIdx,			// given starting suffix index (0..N) for the thread instance
					  int64_t *pEndSfxIdx)			// returned end suffix (inclusive) index to be processed by this thread
{
int Ofs;
int64_t TargPsn1;
int64_t TargPsn2;
etSeqBase *pEl1;
etSeqBase *pEl2;
int64_t NumSfxEls;
int64_t PutativeRange;
int64_t PutativeEndSfxEl;
void *pSfxArray;
etSeqBase *pTarg;


if(NumThreads < 1 || NumThreads > 64 ||
	ThreadInst < 1 || ThreadInst > NumThreads ||
	StartSfxIdx < 0 || pEndSfxIdx == NULL)
	return(-1);
*pEndSfxIdx = 0;
if(m_pSfxBlock == NULL)
	return(-1);
NumSfxEls = m_pSfxBlock->ConcatSeqLen;
if(StartSfxIdx >= NumSfxEls)
	return(-1);
pTarg = (etSeqBase *)&m_pSfxBlock->SeqSuffix[0];
pSfxArray = (void *)&m_pSfxBlock->SeqSuffix[m_pSfxBlock->ConcatSeqLen];

if((NumSfxEls - StartSfxIdx) < 50000)		// if less than 50K putative K_Mers then let the current thread handle them all
	{
	*pEndSfxIdx = NumSfxEls-1;
	return(1);
	}
PutativeRange = (NumSfxEls - StartSfxIdx) / (1 + NumThreads - ThreadInst);
PutativeEndSfxEl = StartSfxIdx + PutativeRange;
if((PutativeEndSfxEl + 25000) >= NumSfxEls)	// any other thread shouldn't have just 25K to process, let this thread handle remaining
	{
	*pEndSfxIdx = NumSfxEls-1;
	return(1);
	}

TargPsn1 = SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,PutativeEndSfxEl);
PutativeEndSfxEl += 1;
for(;PutativeEndSfxEl < NumSfxEls; PutativeEndSfxEl++)
	{
	TargPsn2 = SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,PutativeEndSfxEl);
	if((TargPsn2 + KMerLen) >= NumSfxEls)
		{
		*pEndSfxIdx = PutativeEndSfxEl;
		return(1);
		}
	pEl1 = &pTarg[TargPsn1];
	pEl2 = &pTarg[TargPsn2];
	for(Ofs=0; Ofs < KMerLen; Ofs++,pEl1++, pEl2++)
		{
		if((*pEl1 & 0x0f) > eBaseT || (*pEl1 & 0x0f) != (*pEl2 & 0x0f))	// must only contain ACGT
			{
			*pEndSfxIdx = PutativeEndSfxEl;
			return(1);
			}
		}
	}

return(0);	// nothing for this thread to process
}


int64_t						// < 0 if errors, otherwise the total number of identified K-Mers over both strands (could be more than reported if MinCultivars > 1
CSfxArray::GenKMerCultsCnts(bool bSenseOnly,			// true if sense strand only processing, default is to process both sense and antisense
							int64_t StartSfxIdx,			// starting suffix index
							int64_t EndSfxIdx,			// finish processing at this suffix index, inclusive - if 0 then process all remaining
							int PrefixKMerLen,			// report on K-Mers having this prefix sequence length
						    int SuffixKMerLen,			// and allow for the K-mers containing suffixes of this length (can be 0) 
							int MinCultivars,			// only report if K-Mers present in at least this many different cultivars (0 if must be present in all cultivars)
							int MaxHomozygotic,			// only report prefixes if K-Mer suffixes are homozygotic between a maximum of this many cultivars, if 0 then no check for homozygotic suffixes
							void *pThis,				// callers class instance
						    int (* pCallback)(void *pThis,tsKMerCultsCnts *pCultsCnts)) // callback on K-Mers
{
int Rslt;
tsKMerCultsCnts KMerCultsCnts;					// accepted prefix sequence counts
int64_t NumSfxEntries;
int64_t NumSfxEls;
int64_t NumKMersLocated;
int TargSeqLen;
int CultivarIdx;
int64_t TargPsn1;
int64_t TargPsn2;
int64_t ConcatSeqLen;
tsSfxEntry *pEntry1;
tsSfxEntry *pEntry2;
etSeqBase *pEl1;
etSeqBase *pEl2;
int Ofs;
int64_t HomoSfxIdx;
int64_t TargHomoPsn;
etSeqBase *pElPrefix;
tsSfxEntry *pHomoEntry;
int CultivarsHomozygotic[cMaxCultivars];
int NumCultivarsHomozygotic;
etSeqBase AntisenseKMer[cTotCultivarKMerLen+1];
uint32_t HitEntryID;
uint32_t HitLoci;
int64_t PrevHitIdx;
int64_t CurHitIdx;

etSeqBase *pTarg;			// target sequence
void *pSfxArray;			// target sequence suffix array

// check prefix and suffix lengths within acceptable ranges
TargSeqLen = PrefixKMerLen + SuffixKMerLen;
if(PrefixKMerLen < cMinCultivarPreSufLen || SuffixKMerLen < 0 ||
   PrefixKMerLen > cMaxCultivarPreSufLen || SuffixKMerLen > cMaxCultivarPreSufLen ||
   MinCultivars < 0 ||
   TargSeqLen > cTotCultivarKMerLen ||
   StartSfxIdx < 0 || EndSfxIdx < 0)
	return(-1);


// ensure suffix loaded for iteration
if(m_pSfxBlock == NULL)
	return(-1);

if(m_pSfxBlock->NumEntries > cMaxCultivars)
	return(-1);

if((uint32_t)MinCultivars > m_pSfxBlock->NumEntries)
	return(-1);

if(MinCultivars == 0)
	MinCultivars = m_pSfxBlock->NumEntries;

if(StartSfxIdx > (int64_t)m_pSfxBlock->ConcatSeqLen)
	return(-1);

if(EndSfxIdx == 0 || EndSfxIdx >= (int64_t)m_pSfxBlock->ConcatSeqLen)
	EndSfxIdx = m_pSfxBlock->ConcatSeqLen - 1;

pTarg = (etSeqBase *)&m_pSfxBlock->SeqSuffix[0];
pSfxArray = (void *)&m_pSfxBlock->SeqSuffix[m_pSfxBlock->ConcatSeqLen];
ConcatSeqLen = m_pSfxBlock->ConcatSeqLen;
NumSfxEntries = m_pSfxBlock->NumEntries;
NumSfxEls = m_pSfxBlock->ConcatSeqLen;

Rslt = 0;
NumKMersLocated = 0;

do {
		// find initial K-Mer sequence which is the first instance of that sequence
	KMerCultsCnts.NumCultivars = 0;
	for(;StartSfxIdx <= EndSfxIdx;StartSfxIdx++)
		{
		TargPsn1 = SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,StartSfxIdx);
		// potential starting sequence, check it's length, only containing ACGT is at least as long as the prefix+suffix required
		if((TargPsn1 + TargSeqLen) >= ConcatSeqLen)	
			continue;				// keep looking..
		pEl1 = &pTarg[TargPsn1];
		for(Ofs=0; Ofs < TargSeqLen; Ofs++,pEl1++)
			if((*pEl1 & 0x0f) > eBaseT)	// must only contain ACGT
				break;
		if(Ofs == TargSeqLen)
			break;						// meets criteria required length with out any bases other than ACGT
		}
	if(StartSfxIdx > EndSfxIdx)			// if unable to locate initial K-mer sequence
		break;							// will break out of enclosing so..while and return number of unique K-Mers located

	// if also suffix processing then check that this that the suffixes for this prefix are homozygotic over at most MaxHomozygotic cultivars
	if(SuffixKMerLen && MaxHomozygotic)
		{
		memset(CultivarsHomozygotic,0,sizeof(CultivarsHomozygotic));
		NumCultivarsHomozygotic = 0;

		pHomoEntry = MapChunkHit2Entry(TargPsn1);
		CultivarsHomozygotic[pHomoEntry->EntryID - 1] = 1;
		NumCultivarsHomozygotic = 1;
		pElPrefix = &pTarg[TargPsn1];

		// antisense so need to process against the prefix reverse complemented
		if(!bSenseOnly)
			{
			memmove(AntisenseKMer,pElPrefix,TargSeqLen);
			CSeqTrans::ReverseComplement(TargSeqLen,AntisenseKMer);
			AntisenseKMer[TargSeqLen]= eBaseEOS;
			PrevHitIdx = 0;
			while((CurHitIdx = IterateExacts(AntisenseKMer,TargSeqLen,PrevHitIdx,&HitEntryID,&HitLoci)) > 0)
				{
				PrevHitIdx = CurHitIdx;
				if(CultivarsHomozygotic[HitEntryID-1] == 0)
					{
					NumCultivarsHomozygotic += 1;
					CultivarsHomozygotic[HitEntryID-1] = 1;
					}
				}
			}

		// check how many other cultivars have this exact matching prefix + suffix
		for(HomoSfxIdx = StartSfxIdx+1;HomoSfxIdx <= EndSfxIdx;HomoSfxIdx++)
			{
			TargHomoPsn = SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,HomoSfxIdx);
			if((TargHomoPsn + TargSeqLen) >= ConcatSeqLen)	
				break;
			pEl1 = pElPrefix;
			pEl2 = &pTarg[TargHomoPsn];
			for(Ofs=0; Ofs < TargSeqLen; Ofs++,pEl1++,pEl2++)
				if((*pEl1 & 0x0f) != (*pEl2 & 0x0f))	
					break;								
			if(Ofs < PrefixKMerLen || (Ofs < TargSeqLen && ((*pEl2 & 0x0f) > eBaseT)))
				break;												// onto a new prefix

			// still on same prefix, what about the suffix?
			if(Ofs < TargSeqLen)									
				{
				// check if remainder of suffix contains any no-cannonical bases
				for(; Ofs < TargSeqLen; Ofs++,pEl2++)
					if((*pEl2 & 0x0f) > eBaseT)	
						break;
				if(Ofs < TargSeqLen)
					continue;
				// suffix has changed but still same prefix
				pElPrefix = &pTarg[TargHomoPsn];					// keep checking with same prefix but the new suffix
				if(!bSenseOnly)										// check antisense?
					{
					memmove(AntisenseKMer,pElPrefix,TargSeqLen);
					CSeqTrans::ReverseComplement(TargSeqLen,AntisenseKMer);
					AntisenseKMer[TargSeqLen]= eBaseEOS;
					PrevHitIdx = 0;
					while((CurHitIdx = IterateExacts(AntisenseKMer,TargSeqLen,PrevHitIdx,&HitEntryID,&HitLoci)) > 0)
						{
						PrevHitIdx = CurHitIdx;
						if(CultivarsHomozygotic[HitEntryID-1] == 0)
							{
							NumCultivarsHomozygotic += 1;
							CultivarsHomozygotic[HitEntryID-1] = 1;
							}
						}
					}
				continue;
				}

			// prefix + suffix are homozygotic 
			pHomoEntry = MapChunkHit2Entry(TargHomoPsn);			// prefix + suffix still identical
			if(CultivarsHomozygotic[pHomoEntry->EntryID-1] == 0)
				{
				NumCultivarsHomozygotic += 1;
				CultivarsHomozygotic[pHomoEntry->EntryID-1] = 1;
				}
			}
		if(NumCultivarsHomozygotic > MaxHomozygotic)
			{
			StartSfxIdx = HomoSfxIdx;
			continue;
			}
		}
 

	// have an initial starting K-Mer, initialise cultivar counts and entry identifiers ...
	memset(&KMerCultsCnts,0,sizeof(KMerCultsCnts));	
	KMerCultsCnts.TotCultivars = m_pSfxBlock->NumEntries;
	for(CultivarIdx = 0; CultivarIdx < (int)KMerCultsCnts.TotCultivars; CultivarIdx++)
		KMerCultsCnts.CultCnts[CultivarIdx].EntryID = CultivarIdx + 1;
	// add the starting K-Mer sequence
	pEntry1 = MapChunkHit2Entry(TargPsn1);
	pEl1 = &pTarg[TargPsn1];

	memmove(KMerCultsCnts.KMerSeq,pEl1,TargSeqLen);
	KMerCultsCnts.KMerSeq[TargSeqLen] = eBaseEOS;
	KMerCultsCnts.NumCultivars = 1;
	KMerCultsCnts.SenseCnts = 1;
	KMerCultsCnts.CultCnts[pEntry1->EntryID-1].SenseCnts = 1;

	if(StartSfxIdx++ == EndSfxIdx)
		break;

	for(;StartSfxIdx <= EndSfxIdx;StartSfxIdx++)
		{
		TargPsn2 = SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,StartSfxIdx);
		if((TargPsn2 + TargSeqLen) >= ConcatSeqLen)			// putative would be too short?
			break;
		// not too short but may contain bases other than ACGT so check for these whilst comparing against the initial K_Mer sequence
		pEl1 = KMerCultsCnts.KMerSeq;						
		pEl2 = &pTarg[TargPsn2];
		for(Ofs=0; Ofs < PrefixKMerLen; Ofs++,pEl1++,pEl2++)	// bases must be identical over prefix between initial and putative
			if((*pEl1 & 0x0f) != (*pEl2 & 0X0f))				// if not identical then will need to start a new prefix
				break;
		if(Ofs < PrefixKMerLen)
			break;												// start a new prefix after callbacking with current prefix counts

		// prefix still identical, if suffix  processing requested then ensure suffix only contains cannonical ACGT bases
		if(SuffixKMerLen)
			{
			for(Ofs = 0; Ofs < SuffixKMerLen; Ofs++,pEl2++)
				if((*pEl2 & 0x0f) > eBaseT)
					break;
			if(Ofs < SuffixKMerLen)
				continue;										// unable to extend suffix out to required length so keep looking with same prefix
			}

		// putative K-Mer has an identical prefix and, if also suffix processing, the suffix only contains ACGTs
		// accumulate counts
		pEntry2 = MapChunkHit2Entry(TargPsn2);
		if(KMerCultsCnts.CultCnts[pEntry2->EntryID-1].SenseCnts == 0 && KMerCultsCnts.CultCnts[pEntry2->EntryID-1].AntisenseCnts == 0)
			KMerCultsCnts.NumCultivars += 1;
		KMerCultsCnts.CultCnts[pEntry2->EntryID-1].SenseCnts += 1;
		KMerCultsCnts.SenseCnts += 1;
		}

	// if at least one cultivar with counts then may be sufficent cultivars to report via the callback
	if(KMerCultsCnts.NumCultivars)
		{
		NumKMersLocated += 1;
		// if checking antisense also then revcpl the prefix and look for exact matches on this
		if(!bSenseOnly)			// true if sense strand only processing, default is to process both sense and antisense
			AntisenseKMerCultsCnts(PrefixKMerLen,SuffixKMerLen,&KMerCultsCnts);

		// if number of cultivars in which K-Mers with identical sequences discovered is more than minimum required then report back to caller 
		if(KMerCultsCnts.NumCultivars >= (uint32_t)MinCultivars)
			{
			Rslt = (*pCallback)(pThis,&KMerCultsCnts);
			if(Rslt < 0)
				return(Rslt);
			}
		KMerCultsCnts.NumCultivars = 0;
		}
	}
while(StartSfxIdx <= EndSfxIdx);
if(Rslt < 0)
	return(Rslt);
if(KMerCultsCnts.NumCultivars)
	{
	NumKMersLocated += 1;
	// if checking antisense also then revcpl the prefix and look for exact matches on this
	if(!bSenseOnly)			// true if sense strand only processing, default is to process both sense and antisense
		AntisenseKMerCultsCnts(PrefixKMerLen,SuffixKMerLen,&KMerCultsCnts);

	// if number of cultivars in which K-Mers with identical sequences discovered is more than minimum required then report back to caller 
	if(KMerCultsCnts.NumCultivars >= (uint32_t)MinCultivars)
		{
		Rslt = (*pCallback)(pThis,&KMerCultsCnts);
		if(Rslt < 0)
			return(Rslt);
		}
	}
return(NumKMersLocated);
}

int64_t																// returned number of antisense K-Mers identified
CSfxArray::AntisenseKMerCultsCnts(int PrefixKMerLen,				// report on K-Mers having this prefix sequence length
						    int SuffixKMerLen,					// and allow for the K-mers containing suffixes of this length (can be 0) 
						   tsKMerCultsCnts *pKMerCultsCnts)		// update these cultivar counts
{
int Ofs;
etSeqBase *pSeq;
uint32_t HitEntryID;
uint32_t HitLoci;
int64_t TargPsn;
etSeqBase *pTarg;			// target sequence
void *pSfxArray;			// target sequence suffix array
int64_t PrevHitIdx;
int64_t CurHitIdx;
int64_t NumAntisenseKMers;

etSeqBase AntisensePrefix[cTotCultivarKMerLen+1];

pTarg = (etSeqBase *)&m_pSfxBlock->SeqSuffix[0];
pSfxArray = (void *)&m_pSfxBlock->SeqSuffix[m_pSfxBlock->ConcatSeqLen];

// antisense so need to process against the prefix reverse complemented
memmove(AntisensePrefix,pKMerCultsCnts->KMerSeq,PrefixKMerLen);
CSeqTrans::ReverseComplement(PrefixKMerLen,AntisensePrefix);
AntisensePrefix[PrefixKMerLen]= eBaseEOS;

// locate 1st instance of the antisense K-Mer if one exists
PrevHitIdx = 0;
CurHitIdx;
NumAntisenseKMers = 0;
while((CurHitIdx = IterateExacts(AntisensePrefix,PrefixKMerLen,PrevHitIdx,&HitEntryID,&HitLoci)) > 0)
	{
	PrevHitIdx = CurHitIdx;
	// may have a hit on the prefix but also need to ensure that the suffix will be sufficiently long
	if(SuffixKMerLen > 0)
		{
		TargPsn = SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,CurHitIdx-1);
		pSeq = &pTarg[TargPsn + PrefixKMerLen];
		for(Ofs = 0; Ofs < SuffixKMerLen; Ofs++,pSeq++)
			if((*pSeq & 0x0f) > eBaseT)
				break;
		if(Ofs < SuffixKMerLen)
			continue;
		}
	// accepted this antisense hit, record cnt
	if(pKMerCultsCnts->CultCnts[HitEntryID-1].SenseCnts == 0 && pKMerCultsCnts->CultCnts[HitEntryID-1].AntisenseCnts == 0)
		pKMerCultsCnts->NumCultivars += 1;
	pKMerCultsCnts->CultCnts[HitEntryID-1].AntisenseCnts += 1;
	pKMerCultsCnts->AntisenseCnts += 1;
	NumAntisenseKMers += 1;
	}
return(NumAntisenseKMers);
}

int64_t						// returned number of cultivar sequences
CSfxArray::LocateMultiCultivarMarkers(int PrefixKMerLen,	// prefix K-mer length
						   int SuffixKMerLen,			// suffix K-mer length
						   int NumCultivars,			// number of cultivars in pCultivars
						   uint32_t *pEntryIDs,			// sequences of interest will be on these entries
						   void *pThis,					// callers class instance
						   int (* pCallback)(void *pThis,uint32_t EntryID,etSeqBase *pSeq)) // callback on cultivar marker sequences
{
uint32_t CultivarFlags[cMaxCultivars];
int64_t NumSfxEntries;
int64_t NumSfxEls;
int64_t SfxElsIdx1;
int64_t SfxElsIdx2;
int TargSeqLen;
int64_t MarkSfxElsIdx;
int CultivarIdx;
int NumCultivarsPrefixed;
int64_t TargPsn1;
int64_t TargPsn2;
int64_t NxtTargPsn;
int64_t ConcatSeqLen;
tsSfxEntry *pEntry1;
tsSfxEntry *pEntry2;
tsSfxEntry *pProvisionalEntry;
etSeqBase *pProvisionalSeq;
etSeqBase *pEl1;
etSeqBase *pEl2;
int Ofs;

etSeqBase *pTarg;			// target sequence
void *pSfxArray;			// target sequence suffix array

if(NumCultivars < 1 || NumCultivars > cMaxCultivars || pEntryIDs == NULL ||
   PrefixKMerLen < cMinCultivarPreSufLen || SuffixKMerLen < cMinCultivarPreSufLen ||
   PrefixKMerLen > cMaxCultivarPreSufLen || SuffixKMerLen > cMaxCultivarPreSufLen)
	return(0);

// ensure suffix loaded for iteration
if(m_pSfxBlock == NULL)
	return(0);

// check that all entry identifiers are known ...
for(CultivarIdx = 0; CultivarIdx < NumCultivars; CultivarIdx++)
	if(pEntryIDs[CultivarIdx] < 1 || pEntryIDs[CultivarIdx] > m_pSfxBlock->NumEntries)
		return(0);

pTarg = (etSeqBase *)&m_pSfxBlock->SeqSuffix[0];
pSfxArray = (void *)&m_pSfxBlock->SeqSuffix[m_pSfxBlock->ConcatSeqLen];
ConcatSeqLen = m_pSfxBlock->ConcatSeqLen;
NumSfxEntries = m_pSfxBlock->NumEntries;
NumSfxEls = m_pSfxBlock->ConcatSeqLen;

 // start processing
// first identify the prefixes which are common to the cultivars
TargSeqLen = PrefixKMerLen + SuffixKMerLen;
MarkSfxElsIdx = 0;
NumCultivarsPrefixed = 0;
for(SfxElsIdx1 = 0; SfxElsIdx1 < (NumSfxEls-1); SfxElsIdx1++)
	{
	TargPsn1 = SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,SfxElsIdx1);
	// slough if remaining concatenated sequences too short, or current sequence is early terminated with a conctenation 'N'
	if((TargPsn1 + TargSeqLen) >= ConcatSeqLen)	
		continue;
	pEl1 = &pTarg[TargPsn1];
	for(Ofs=0; Ofs < TargSeqLen; Ofs++,pEl1++)
		if((*pEl1 & 0x0f) > eBaseT)
			break;
	if(Ofs < TargSeqLen)
		continue;

	// current sequence is long enough to explore
	pEl1 = &pTarg[TargPsn1];
	pEntry1 = NULL;
	if(NumCultivarsPrefixed)				// only clear cultivar flags if it is known at least one has flags set
		{
		memset(CultivarFlags,0,sizeof(CultivarFlags));
		NumCultivarsPrefixed = 0;
		}
	pProvisionalEntry = NULL;
	pProvisionalSeq = NULL;
	NxtTargPsn = TargPsn1;
	for(SfxElsIdx2 = SfxElsIdx1 + 1; SfxElsIdx2 < NumSfxEls; SfxElsIdx2++)
		{
		TargPsn1 = NxtTargPsn;
		TargPsn2 = SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,SfxElsIdx2);
		// slough if remaining concatenated sequences too short, or current sequence is early terminated with a conctenation 'N'
		if((TargPsn2 + TargSeqLen) >= ConcatSeqLen)				
			continue;
		pEl2 = &pTarg[TargPsn2];
		for(Ofs=0; Ofs < TargSeqLen; Ofs++,pEl2++)
			if((*pEl2 & 0x0f) > eBaseT)
				break;
		if(Ofs < TargSeqLen)
			continue;

		// check if prefix still matching
		pEl1 = &pTarg[TargPsn1];
		pEl2 = &pTarg[TargPsn2];
		for(Ofs=0; Ofs < PrefixKMerLen; Ofs++,pEl1++,pEl2++)
			if((*pEl1 & 0x0f) != (*pEl2 & 0x0f))
				break;
		if(Ofs != PrefixKMerLen)								// break if no longer sharing the same prefix
			break;

		// sharing same prefix
		NxtTargPsn = TargPsn2;

		// sharing same prefix, check now if also sharing same suffix 
		for(Ofs=0; Ofs < SuffixKMerLen; Ofs++,pEl1++,pEl2++)
			if((*pEl1 & 0x0f) != (*pEl2 & 0x0f))
				break;
		if(Ofs == SuffixKMerLen)								// exactly matching?
			{
			// if exclusively on same cultivar then will accept, otherwise slough
			if(pEntry1 == NULL)
				pEntry1 = MapChunkHit2Entry(TargPsn1);
			pEntry2 = MapChunkHit2Entry(TargPsn2);
			if(pEntry1->EntryID != pEntry2->EntryID)			// will slough as not cultivar specific
				{
				pProvisionalEntry = NULL;

				continue;
				}
			// will need to provisionally accept...
			if(pProvisionalEntry != NULL)
				{
				pProvisionalEntry = pEntry1;
				pProvisionalSeq = &pTarg[TargPsn1];
				}
			}

		// different suffix
		// if provisionally accepted exists then that sequence can be reported
		if(pProvisionalEntry != NULL)
			{
			if(pCallback(pThis,pProvisionalEntry->EntryID,pProvisionalSeq))
				break;
			pProvisionalEntry = NULL;
			pProvisionalSeq = NULL;
			}
		if(NumCultivarsPrefixed < NumCultivars)
			{
			if(pEntry1 == NULL)
				{
				pEntry1 = MapChunkHit2Entry(SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,SfxElsIdx1));
				NumCultivarsPrefixed += 1;
				CultivarFlags[pEntry1->EntryID-1] |= 0x01;	
				}
			pEntry2 = MapChunkHit2Entry(SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,SfxElsIdx2));
			if(!(CultivarFlags[pEntry2->EntryID - 1] & 0x01))	// if previously unmarked then can mark this cultivar as having a prefix match
				{
				NumCultivarsPrefixed += 1;
				CultivarFlags[pEntry2->EntryID-1] |= 0x01;		
				}
			}
		}

	// if provisionally accepted exists then that sequence can be reported
	if(pProvisionalEntry != NULL)
		{
		if(pCallback(pThis,pProvisionalEntry->EntryID,pProvisionalSeq))
			break;
		pProvisionalEntry = NULL;
		pProvisionalSeq = NULL;
		}

	SfxElsIdx1 = SfxElsIdx2;
	}

return(0);		// no more hits
}


int64_t 
CSfxArray::IterateExacts(etSeqBase *pProbeSeq,	// probe
 						 uint32_t ProbeLen,			// probe length
						 int64_t PrevHitIdx,			// 0 if starting new sequence, otherwise set to return value of previous successful iteration return
 						 uint32_t *pTargEntryID,		// if match then where to return suffix entry (chromosome) matched on
						 uint32_t *pHitLoci,			// if match then where to return loci
						 uint32_t MaxExtdLen,			// attempt to extend exact match of ProbeLen out to MaxExtdLen and report extended match length in *pHitExtdLen
						int *pHitExtdLen)			// was able to extend core out to this length
{
int Cmp;
int64_t TargPsn;
tsSfxEntry *pEntry;
etSeqBase *pEl1;
etSeqBase *pEl2;

etSeqBase *pTarg;			// target sequence
void *pSfxArray;			// target sequence suffix array
int64_t SfxLen;				// number of suffixs in pSfxArray

if(MaxExtdLen == 0 || pHitExtdLen == NULL)
	return(IterateExacts(pProbeSeq,ProbeLen,PrevHitIdx,pTargEntryID,pHitLoci));

*pHitLoci = 0;
*pTargEntryID = 0;
*pHitExtdLen = 0;

// ensure suffix loaded for iteration and prev hit was not the last!
if(m_pSfxBlock == NULL || (uint64_t)PrevHitIdx >= m_pSfxBlock->ConcatSeqLen)
	return(0);

pTarg = (etSeqBase *)&m_pSfxBlock->SeqSuffix[0];
pSfxArray = (void *)&m_pSfxBlock->SeqSuffix[m_pSfxBlock->ConcatSeqLen];
SfxLen = m_pSfxBlock->ConcatSeqLen;

if(!PrevHitIdx)
	{
	// locate first exact match
	if((TargPsn = LocateFirstExact(pProbeSeq,ProbeLen,pTarg,m_pSfxBlock->SfxElSize,pSfxArray,0,0,SfxLen-1)) == 0)
		return(0);	// no match
	TargPsn -= 1;

	pEl2 = &pTarg[SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,TargPsn)];
	pEl1 = pProbeSeq;
	*pHitExtdLen = ExactMatchLen(pEl1,pEl2,MaxExtdLen);
	pEntry = MapChunkHit2Entry(SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,TargPsn));
	*pTargEntryID = pEntry->EntryID;
	*pHitLoci = (uint32_t)(SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,TargPsn) - pEntry->StartOfs);
	return(TargPsn + 1);
	}

  // check if probe matches next suffix
pEl2 = &pTarg[SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,PrevHitIdx)];
pEl1 = pProbeSeq;

if(m_bBisulfite)
	Cmp = BSCmpProbeTarg(pEl1,pEl2,ProbeLen);
else
	Cmp = CmpProbeTarg(pEl1,pEl2,ProbeLen);
if(!Cmp)
	{
	*pHitExtdLen = ExactMatchLen(pEl1,pEl2,MaxExtdLen);
	pEntry = MapChunkHit2Entry(SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,PrevHitIdx));
	*pTargEntryID = pEntry->EntryID;
	*pHitLoci = (uint32_t)(SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,PrevHitIdx) - pEntry->StartOfs);
	return(PrevHitIdx+1);
	}

return(0);		// no more hits


}


int64_t						// returned hit index (1..n) or 0 if no hits
CSfxArray::IterateExacts(etSeqBase *pProbeSeq,// probe
 						 uint32_t ProbeLen,		// probe length
						 int64_t PrevHitIdx,		// 0 if starting new sequence, otherwise set to return value of previous successful iteration return
 						 uint32_t *pTargEntryID,	// if match then where to return suffix entry (chromosome) matched on
						 uint32_t *pHitLoci,		// if match then where to return loci
						 uint32_t *pTargSeqLen,	// optionally update with the matched target sequence length
						 etSeqBase **ppTargSeq) // optionally update with ptr to start of the target sequence, exact match will have started at &ppTargSeq[*pHitLoci]
{
int Cmp;
int64_t TargPsn;
tsSfxEntry *pEntry;
etSeqBase *pEl1;
etSeqBase *pEl2;

etSeqBase *pTarg;			// target sequence
void *pSfxArray;			// target sequence suffix array
int64_t SfxLen;				// number of suffixs in pSfxArray


*pHitLoci = 0;
*pTargEntryID = 0;

if(pTargSeqLen != NULL)
	*pTargSeqLen = 0;
if(ppTargSeq != NULL)
	*ppTargSeq = NULL;

// ensure suffix loaded for iteration and prev hit was not the last!
if(m_pSfxBlock == NULL || (uint64_t)PrevHitIdx >= m_pSfxBlock->ConcatSeqLen)
	return(0);

pTarg = (etSeqBase *)&m_pSfxBlock->SeqSuffix[0];
pSfxArray = (void *)&m_pSfxBlock->SeqSuffix[m_pSfxBlock->ConcatSeqLen];
SfxLen = m_pSfxBlock->ConcatSeqLen;

if(!PrevHitIdx)
	{
	// locate first exact match
	if((TargPsn = LocateFirstExact(pProbeSeq,ProbeLen,pTarg,m_pSfxBlock->SfxElSize,pSfxArray,0,0,SfxLen-1)) == 0)
		return(0);	// no match
	TargPsn -= 1;
	pEntry = MapChunkHit2Entry(SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,TargPsn));
	*pTargEntryID = pEntry->EntryID;
	*pHitLoci = (uint32_t)(SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,TargPsn) - pEntry->StartOfs);
	
	if(pTargSeqLen != NULL)
		*pTargSeqLen = pEntry->SeqLen;
	if(ppTargSeq != NULL)
		*ppTargSeq = &pTarg[pEntry->StartOfs];

	return(TargPsn + 1);
	}

  // check if probe matches next suffix
pEl2 = &pTarg[SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,PrevHitIdx)];
pEl1 = pProbeSeq;

if(m_bBisulfite)
	Cmp = BSCmpProbeTarg(pEl1,pEl2,ProbeLen);
else
	Cmp = CmpProbeTarg(pEl1,pEl2,ProbeLen);
if(!Cmp)
	{
	pEntry = MapChunkHit2Entry(SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,PrevHitIdx));
	*pTargEntryID = pEntry->EntryID;
	*pHitLoci = (uint32_t)(SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,PrevHitIdx) - pEntry->StartOfs);

	if(pTargSeqLen != NULL)
		*pTargSeqLen = pEntry->SeqLen;
	if(ppTargSeq != NULL)
		*ppTargSeq = &pTarg[pEntry->StartOfs];

	return(PrevHitIdx+1);
	}

return(0);		// no more hits
}

int64_t						// returned hit index (1..n) or 0 if no hits
CSfxArray::IterateExactsRange(etSeqBase *pProbeSeq,	// probe sequence
				uint32_t ProbeLen,		// probe length
				int64_t PrevHitIdx,		// 0 if starting new sequence, otherwise set to return value of previous successful iteration return
				uint32_t TargEntryID,		// accepted hits must be onto this target suffix entry (chromosome)
				uint32_t StartLoci,		// accepted hit loci must be in the range StartLoci ... EndLoci inclusive
				uint32_t EndLoci,			// accepted hit loci must be in the range StartLoci ... EndLoci inclusive
				uint32_t *pHitLoci,		// if match then where to return loci
				uint32_t *pTargSeqLen,	// optionally update with the matched target sequence length
				etSeqBase **ppTargSeq)	// optionally update with ptr to start of the target sequence, exact match will have started at &ppTargSeq[*pHitLoci]
{
	int Cmp;

	int64_t TargPsn;
	uint32_t HitLoci;
	tsSfxEntry *pEntry;
	etSeqBase *pEl1;
	etSeqBase *pEl2;

	etSeqBase *pTarg;			// target sequence
	void *pSfxArray;			// target sequence suffix array
	int64_t SfxLen;				// number of suffixs in pSfxArray
	*pHitLoci = 0;

	if (pTargSeqLen != NULL)
		*pTargSeqLen = 0;
	if (ppTargSeq != NULL)
		*ppTargSeq = NULL;

	if(TargEntryID == 0)
		return(0);

	// ensure suffix loaded for iteration and prev hit was not the last!
	if (m_pSfxBlock == NULL || (uint64_t)PrevHitIdx >= m_pSfxBlock->ConcatSeqLen)
		return(0);

	pTarg = (etSeqBase *)&m_pSfxBlock->SeqSuffix[0];
	pSfxArray = (void *)&m_pSfxBlock->SeqSuffix[m_pSfxBlock->ConcatSeqLen];
	SfxLen = m_pSfxBlock->ConcatSeqLen;

	if (!PrevHitIdx)
		{
		// locate first exact match
		if ((TargPsn = LocateFirstExact(pProbeSeq, ProbeLen, pTarg, m_pSfxBlock->SfxElSize, pSfxArray, 0, 0, SfxLen - 1)) == 0)
			return(0);	// no match
		TargPsn -= 1;
		pEntry = MapChunkHit2Entry(SfxOfsToLoci(m_pSfxBlock->SfxElSize, pSfxArray, TargPsn));
		if(pEntry->EntryID == TargEntryID)
			{
			HitLoci = (uint32_t)(SfxOfsToLoci(m_pSfxBlock->SfxElSize, pSfxArray, TargPsn) - pEntry->StartOfs);
			if(HitLoci >= StartLoci && HitLoci <= EndLoci)
				{
				*pHitLoci = HitLoci;
				if (pTargSeqLen != NULL)
					*pTargSeqLen = pEntry->SeqLen;
				if (ppTargSeq != NULL)
					*ppTargSeq = &pTarg[pEntry->StartOfs];
				return(TargPsn + 1);
				}
			}
		PrevHitIdx = TargPsn + 1;
		}

do {
	// check if probe matches next suffix
	pEl2 = &pTarg[SfxOfsToLoci(m_pSfxBlock->SfxElSize, pSfxArray, PrevHitIdx)];
	pEl1 = pProbeSeq;

	if (m_bBisulfite)
		Cmp = BSCmpProbeTarg(pEl1, pEl2, ProbeLen);
	else
		Cmp = CmpProbeTarg(pEl1, pEl2, ProbeLen);
	if (!Cmp)
		{
		pEntry = MapChunkHit2Entry(SfxOfsToLoci(m_pSfxBlock->SfxElSize, pSfxArray, PrevHitIdx));
		if (pEntry->EntryID == TargEntryID)
			{
			HitLoci = (uint32_t)(SfxOfsToLoci(m_pSfxBlock->SfxElSize, pSfxArray, PrevHitIdx) - pEntry->StartOfs);
			if (HitLoci >= StartLoci && HitLoci <= EndLoci)
				{
				*pHitLoci = HitLoci;
				if (pTargSeqLen != NULL)
					*pTargSeqLen = pEntry->SeqLen;
				if (ppTargSeq != NULL)
					*ppTargSeq = &pTarg[pEntry->StartOfs];
				return(PrevHitIdx + 1);
				}
			}
		PrevHitIdx += 1;
		}
	}
while(!Cmp);
return(0);		// no more hits
}


int											// number of target sequences which were prequalified
CSfxArray::PreQualTargs(uint32_t ProbeEntryID,	// probe sequence entry identifier used to determine if a self hit to be sloughed
			int ProbeSeqLen,				// probe sequence length
			etSeqBase *pProbeSeq,			// prequalify with cores from this probe sequence
			int QualCoreLen,				// prequalifying core lengths
			 int DeltaCoreOfs,				// offset core windows of QualCoreLen along the probe sequence when checking for overlaps
			int TargLenDiffBp,				// 0 if disabled, accepted prequalified targets must have length within this differential (bp) of probe sequence length (used in transcriptome processing)
			int MaxPreQuals,				// max number of target sequences to pre-qualify
			tsQualTarg *pQualTargs)			// into this prequalified list
{
etSeqBase *pHomo;
int HomoIdx;
int HomoBaseCnts[4];
int MaxAcceptHomoCnt;
int NumQualSeqs;
int QualSeqIdx;
int ProbeOfs;
int LastProbeOfs;
etSeqBase *pCoreSeq;
etSeqBase *pTargSeq;
uint32_t TargSeqLen;
int64_t PrevHitIdx;
int64_t NxtHitIdx;
uint32_t TargEntryID;
uint32_t HitLoci;
tsQualTarg *pQualTarg;
tsQualTarg *pFiltCntQualHits;

int TargScoreLen;
int QuickScore;

memset(pQualTargs,0,MaxPreQuals * sizeof(tsQualTarg));
NumQualSeqs = 0;
PrevHitIdx = 0;
LastProbeOfs = ProbeSeqLen - (QualCoreLen + 100);
pCoreSeq = pProbeSeq;
MaxAcceptHomoCnt = (QualCoreLen * 90) / 100; // if any core contains more than 90% of the same base then treat as being a near homopolymer core and slough
for(ProbeOfs = 0; ProbeOfs < LastProbeOfs; ProbeOfs+=DeltaCoreOfs, pCoreSeq+=DeltaCoreOfs)
	{
	// with PacBio reads most homopolymer runs are actual artefact inserts so don't bother processing these homopolymer runs for cores
	if(MaxAcceptHomoCnt > 0)
		{
		HomoBaseCnts[0] = HomoBaseCnts[1] = HomoBaseCnts[2] = HomoBaseCnts[3] = 0;
		for(pHomo = pCoreSeq, HomoIdx = 0; HomoIdx < QualCoreLen; HomoIdx+=1, pHomo += 1)
			HomoBaseCnts[*pHomo & 0x03] += 1;
		if(HomoBaseCnts[0] > MaxAcceptHomoCnt || HomoBaseCnts[1] > MaxAcceptHomoCnt || HomoBaseCnts[2] > MaxAcceptHomoCnt || HomoBaseCnts[3] > MaxAcceptHomoCnt)
			continue;
		}

	PrevHitIdx = 0;
	while((NxtHitIdx = IterateExacts(pCoreSeq,QualCoreLen,PrevHitIdx,&TargEntryID,&HitLoci,&TargSeqLen,&pTargSeq))!=0)
		{
		PrevHitIdx = NxtHitIdx;
		if(ProbeEntryID == TargEntryID)
			continue;
		if((uint32_t)TargLenDiffBp > 0 && ((TargSeqLen > (ProbeSeqLen + (uint32_t)TargLenDiffBp)) || (TargSeqLen < (ProbeSeqLen - (uint32_t)TargLenDiffBp))))
			continue;

	   TargScoreLen = min((uint32_t)150, TargSeqLen - HitLoci); // attempting to quick SW score over up to 150bp with a minimum of 100bp, this range including the QualCoreLen exactly matching core
	   if(TargScoreLen < 100)
			continue;

		// see if exactly matching core can be extended, if extension is out to at least 60% of targeted scoring length then will accept as not needing the cost of a SW
		int QAnchorIdx;
		int MatchBaseLen;
		int OverlapIdentity;
		etSeqBase *pQAnchor;
		etSeqBase *pTAnchor;
		pTargSeq = &pTargSeq[HitLoci];
		QAnchorIdx = QualCoreLen;
		pQAnchor = &pCoreSeq[QAnchorIdx];
		pTAnchor = &pTargSeq[QAnchorIdx];
		for(MatchBaseLen = QualCoreLen; QAnchorIdx < TargScoreLen; QAnchorIdx++,pQAnchor++,pTAnchor++,MatchBaseLen++)
			{
			if((*pQAnchor & 0x07) != (*pTAnchor & 0x07))
				break;
			}
		OverlapIdentity = (MatchBaseLen * 100) / TargScoreLen;  
		if(OverlapIdentity < 60) // if able to get at least 60% of the TargScoreLen exactly matching then not worth the overhead of quickscore processing - just accept!
			{
			QuickScore = QuickScoreOverlap(TargScoreLen - MatchBaseLen,pQAnchor,pTAnchor);
			OverlapIdentity = ((MatchBaseLen + ((TargScoreLen - MatchBaseLen) * QuickScore)/100)*100)/TargScoreLen;
			}
		else
			QuickScore = 0;

	   if(OverlapIdentity < 60) 
			continue;

		pQualTarg = pQualTargs;
		if(NumQualSeqs > 0)
			{
			for(QualSeqIdx = 0; QualSeqIdx < NumQualSeqs; QualSeqIdx++, pQualTarg++)
				{
				if(pQualTarg->TargEntryID == TargEntryID)
					{
					if(pQualTarg->Hits < 0x0ff)
						pQualTarg->Hits += 1;
					break;
					}
				}
			if(QualSeqIdx < NumQualSeqs || QualSeqIdx == MaxPreQuals)
				continue;
			}
		pQualTarg->Hits = 1;
		pQualTarg->TargEntryID = TargEntryID;
		NumQualSeqs += 1;
		}
	}
if(NumQualSeqs < 50)
	return(NumQualSeqs);

int NumFiltCntQualSeqs;
pQualTarg = pQualTargs;
pFiltCntQualHits = pQualTarg;

for(NumFiltCntQualSeqs = QualSeqIdx = 0; QualSeqIdx < NumQualSeqs; QualSeqIdx++, pQualTarg++)
	{
	if(pQualTarg->Hits <= 2)
		continue;
	if(NumFiltCntQualSeqs != QualSeqIdx)
		*pFiltCntQualHits = *pQualTarg;
	pFiltCntQualHits += 1;
	NumFiltCntQualSeqs += 1;
	}
return(NumFiltCntQualSeqs);
}


const int cMinQuickScoreSeedCoreLen = 8;					// expecting QuickScoreOverlap SeedCoreLen to be at least this many bp and
const int cMaxQuickScoreSeedCoreLen = 100;					// no longer than this length
const int cMinQuickScoreOverlapLen = 16;					// QuickScoreOverlap sequence lengths must be at least this many bp and
const int cMaxQuickScoreOverlapLen = 1000;                   // no longer than this maximal bp length

int															// estimated identity (0..100) of overlap between pProbe and pTarg derived by dividing final score by the MatchScore and then normalising by length
CSfxArray::QuickScoreOverlap(int SeedCoreLen,				// initial seed core exact matching length between probe and target was this many bp 
									int SeqLen,				// quick score over this many remaining bases in pProbe and pTarg following the exactly matching SeedCoreLen, must be in the range cMaxQuickScoreOverlapLen..cMaxQuickScoreOverlapLen  
									etSeqBase *pProbe,		// probe sequence scoring onto
									etSeqBase *pTarg,		// this target sequence
								    int MatchScore,			// exact match score  ((in any Pacbio sequence then expecting ~85% of all base calls to be correct), note that between any 2 sequences then the relative error rate doubles!
									int InDelPenalty,		// insertions and deletions much more likely than substitutions (in any Pacbio sequence then expecting ~14% of all error events to be insertions)
							        int SubstitutePenalty)  // base call subs are relatively rare (in any Pacbio sequence then expecting ~1% of all error events to be substitutions)
{
bool bNoScores;
etSeqBase *pP1;
etSeqBase *pT1;
etSeqBase PBase;
etSeqBase TBase;
int HiScore;
int BandScores[cMaxQuickScoreOverlapLen + 1];
int *pCurScore;
int DiagScore; 
int PutativeMatchScore;
int PutativeInsertScore;
int PutativeDeleteScore;

int IdxP1;
int IdxT1;
int NxtHiIdxT1;
int HiIdxT1;
int LoIdxT1;

if(SeedCoreLen < cMinQuickScoreSeedCoreLen ||  SeedCoreLen > cMaxQuickScoreSeedCoreLen ||
	SeqLen < cMinQuickScoreOverlapLen ||  SeqLen > cMaxQuickScoreOverlapLen ||
	pProbe == NULL || pTarg == NULL)
	return(0);

DiagScore = SeedCoreLen * MatchScore;
pCurScore = BandScores;
for(IdxP1 = 0; IdxP1 < min(SeedCoreLen,SeqLen); IdxP1++,pCurScore++)
	{
	*pCurScore = DiagScore;
	DiagScore -= InDelPenalty;
	}
if(IdxP1 < SeqLen+1)
	memset(&BandScores[IdxP1],0,(SeqLen+1 - IdxP1) * sizeof(BandScores[0]));

LoIdxT1 = 0;
pP1 = pProbe;
NxtHiIdxT1 = min(SeedCoreLen+1, SeqLen);
for(IdxP1 = 0; IdxP1 < SeqLen; IdxP1 += 1, pP1 += 1)
	{
	pCurScore = &BandScores[LoIdxT1+1];
	DiagScore = BandScores[LoIdxT1];
	if(BandScores[LoIdxT1] > InDelPenalty)
		BandScores[LoIdxT1] -= InDelPenalty;
	else
		BandScores[LoIdxT1] = 0;
	PBase = *pP1 & 0x07;
	bNoScores = true;
	HiIdxT1 = NxtHiIdxT1;
	pT1 = &pTarg[LoIdxT1+1];
	for(IdxT1 = LoIdxT1; IdxT1 < HiIdxT1; IdxT1 += 1, pT1 += 1, pCurScore++)
		{
		if(*pCurScore == 0 && DiagScore == 0 && pCurScore[-1] == 0)
			{
			if(LoIdxT1 == IdxT1)
				LoIdxT1 += 1;
			continue;
			}
		NxtHiIdxT1 = min(IdxT1+2,SeqLen);
		bNoScores = false;
		TBase = *pT1 & 0x07;
		if(PBase == TBase)			// exactly matches
			PutativeMatchScore = DiagScore + MatchScore;
		else                       // either mismatch or InDel
			PutativeMatchScore = DiagScore - SubstitutePenalty; // assume it was a mismatch
		PutativeInsertScore = *pCurScore - InDelPenalty;			// score as if it was an insertion
		PutativeDeleteScore = pCurScore[-1] - InDelPenalty;		// score as if it was a deletion
		DiagScore = *pCurScore;		
		if(PutativeMatchScore > 0 && PutativeMatchScore >= PutativeInsertScore && PutativeMatchScore >= PutativeDeleteScore)
			*pCurScore = PutativeMatchScore;
		else
			if(PutativeInsertScore > 0 && PutativeInsertScore >= PutativeMatchScore && PutativeInsertScore >= PutativeDeleteScore)
				*pCurScore = PutativeInsertScore;
			else
				if(PutativeDeleteScore > 0 && PutativeDeleteScore >= PutativeMatchScore && PutativeDeleteScore >= PutativeInsertScore)
					*pCurScore = PutativeDeleteScore;
				else
					*pCurScore = 0;
		}
	if(bNoScores)
		return(0);
	}

HiScore = 0;
pCurScore = &BandScores[1];
for(IdxP1 = 0; IdxP1 < SeqLen; IdxP1 += 1, pCurScore++)
	if(*pCurScore > HiScore)
		HiScore = *pCurScore;
if(HiScore > 0)
	HiScore = (int)(((double)(HiScore + MatchScore-1) / (double)MatchScore) * (100.0 / (double)(SeqLen + SeedCoreLen)));
return(HiScore);
}

// QuickScoreOverlap expected to return a count of less than 20 for random sequence overlaps, and >= 40 for real overlapping sequences even with the PacBio error profiles  
int                                         //  estimated identity (0..100) of overlap between pProbe and pTarg, much quicker but far less accurate, than the SW banded function as scoring only based on the number of 4-mers matching within a 16bp window along the two overlapping sequences
CSfxArray::QuickScoreOverlap(int SeqLen,	// both probe and target are of this minimum length (must be at least 16bp), will be clamped to 200bp max
				  etSeqBase *pProbe,		// scoring overlap of probe sequence onto
				  etSeqBase *pTarg)		    // this target sequence
{
	int Matches;
	int Ofs;
	int Idx;
	int MaxBases2Cmp;
	uint32_t Targ16Bases;
	uint32_t Targ4Bases;
	uint32_t Probe4Bases;

	if (SeqLen < 16 || pProbe == NULL || pTarg == NULL)
		return(0);
	if(SeqLen > 200)
		SeqLen = 200;

	Matches = 0;
	Targ16Bases = 0;
	Probe4Bases = 0;

	for (Idx = 0; Idx < 4; Idx++)  // probe sliding window is 4bp packed into 8bits
		{
		Probe4Bases <<= 2;
		Probe4Bases |= *pProbe++ & 0x03;
		}

	for (Idx = 0; Idx < 10; Idx++)  // target window initially starts out as 10bp but will be increased out to a max of 16bp packed into 32bits
		{
		Targ16Bases <<= 2;
		Targ16Bases |= *pTarg++ & 0x03;
		}


	for (Ofs = 0; Ofs < SeqLen; Ofs+=4) // iterating, sliding windows by 4bp, over full sequence length
		{
		Targ4Bases = Targ16Bases;
		MaxBases2Cmp = min(Ofs + 6, 12);
		for (Idx = 0; Idx < MaxBases2Cmp; Idx++)
			{
			if (Probe4Bases == (Targ4Bases & 0x0ff))
				{
				Matches += 1;
				break;
				}
			Targ4Bases >>= 2;
			}

		for (Idx = Ofs; Idx < (Ofs + 4) && Idx < SeqLen; Idx++)  // slide target and probe windows to right by at most 4bp
			{
			Targ16Bases <<= 2;
			Targ16Bases |= *pTarg++ & 0x03;
			Probe4Bases <<= 2;
			Probe4Bases |= *pProbe++ & 0x03;
			}
		Probe4Bases &= 0x0ff;           // probe window is 4bp only!
		}
Matches *= 4;       // 4 bases per match
if(Matches > SeqLen)
	Matches = SeqLen;
return((Matches * 100) / SeqLen);
}


int64_t						// returned hit idex (1..n) or 0 if no hits
CSfxArray::IteratePacBio(etSeqBase *pProbeSeq,				// probe sequence
 									uint32_t ProbeLen,			// probe sequence length
 									 uint32_t SeedCoreLen,		// using this seed core length
									 uint32_t SloughEntryID,		// if > 0 then hits to this entry (normally would be the probe sequence entry identifier) are to be sloughed
									 uint32_t MinTargLen,			// accepted hit target sequences must be at least this length
									 uint32_t MaxTargLen,			// if > 0 then accepted hit target sequences must be no longer than this length
									 int64_t PrevHitIdx,			// 0 if starting new sequence, otherwise set to return value of previous successful iteration return
 									 uint32_t *pTargEntryID,		// if match then where to return suffix entry (chromosome) matched on
									 uint32_t *pHitLoci,			// if match then where to return loci
									  int NumPreQuals,			// number of pre-qualified sequences in pQualTargs
									tsQualTarg *pQualTargs,		// holds prequalified target sequence identifiers
									 int PacBioMinKmersExtn)	// accepting as putative overlap if extension matches at least this many cPacBiokExtnKMerLen (currently 4bp)
{
int Cmp;
bool bFirst;
tsSfxEntry *pEntry;

uint32_t HitLoci;
etSeqBase *pTarg;			// target sequence
void *pSfxArray;			// target sequence suffix array
int64_t SfxLen;				// number of suffixs in pSfxArray
int64_t TargLoci;
int TargKMerLen;

int PreQualIdx;
tsQualTarg *pQualTarg;

etSeqBase *pQAnchor;
etSeqBase *pTAnchor;
int MatchBaseLen;
int OverlapIdentity;
int QAnchorIdx;
int QuickScore;

uint32_t AcceptExactExtdCoreLen = ((cMaxPacBioSeedCoreLen+1)/2); // if core can be simply extended out to AcceptExactExtdCoreLen with exact matching then immediately accept; will be 50bp as cMaxPacBioSeedCoreLen is currently 100bp

etSeqBase *pEl2;


*pHitLoci = 0;
*pTargEntryID = 0;

*pTargEntryID = 0;
*pHitLoci = 0;

// ensure suffix loaded for iteration and prev hit was not the last!
// also require that the probe length must be at least 100 + SeedCoreLen so matching subseqs can be explored
if(SeedCoreLen < cMinPacBioSeedCoreLen || SeedCoreLen > cMaxPacBioSeedCoreLen || ProbeLen < SeedCoreLen || m_pSfxBlock == NULL || (uint64_t)PrevHitIdx >= m_pSfxBlock->ConcatSeqLen)
	return(0);

if(SeedCoreLen >= AcceptExactExtdCoreLen)	// if looking with seed cores of at least 50bp (cMaxPacBioSeedCoreLen currently is 100bp)  then will not bother to seed extend. Alignments should be reasonably high confidence.
	{
	int64_t NxtHitIdx;
	while((NxtHitIdx = IterateExacts(pProbeSeq, SeedCoreLen, PrevHitIdx, pTargEntryID,	pHitLoci)) > 0)
		{
		if(SloughEntryID != 0 && SloughEntryID == *pTargEntryID)
			{
			PrevHitIdx = NxtHitIdx;
			continue;
			}

		pEntry = &m_pEntriesBlock->Entries[*pTargEntryID - 1];
		if(MinTargLen >  pEntry->SeqLen || (MaxTargLen > 0 && MaxTargLen < pEntry->SeqLen))
			{
			PrevHitIdx = NxtHitIdx;
			continue;
			}

		if(NumPreQuals != 0 && pQualTargs != NULL)
			{
			pQualTarg = pQualTargs;
			for(PreQualIdx = 0; PreQualIdx < NumPreQuals; PreQualIdx+=1, pQualTarg+=1)
				{
				if(pQualTarg->TargEntryID == *pTargEntryID)
					break;
				}
			if(PreQualIdx == NumPreQuals)
				{
				PrevHitIdx = NxtHitIdx;
				continue;
				}
			}

		return(NxtHitIdx);
		}	
	*pTargEntryID = 0;
	*pHitLoci = 0;
	return(0);
	}

// will be maximally seed extending potentially out to cPacBioSeedCoreExtn (currently 100bp) so ensure sufficient sequence remaing in prob e sequence to seed extend
if(ProbeLen < cPacBioSeedCoreExtn)
	return(0);

pTarg = (etSeqBase *)&m_pSfxBlock->SeqSuffix[0];
pSfxArray = (void *)&m_pSfxBlock->SeqSuffix[m_pSfxBlock->ConcatSeqLen];
SfxLen = m_pSfxBlock->ConcatSeqLen;
bFirst = false;

if(!PrevHitIdx)	// if locate first exact match using SeedCoreLen
	{
	// check if this core sequence over-occurs or does not exist
	if(OverOccKMerClas(SeedCoreLen,pProbeSeq) != 1)
		return(0);

	if((PrevHitIdx = LocateFirstExact(pProbeSeq,SeedCoreLen,pTarg,m_pSfxBlock->SfxElSize,pSfxArray,0,0,SfxLen-1)) == 0)
		return(0);	// no match
	TargLoci = SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,PrevHitIdx-1);

	if((pEntry = MapChunkHit2Entry(TargLoci))==NULL)
		return(0);

	if(!(pEntry->EntryID == SloughEntryID || MinTargLen > pEntry->SeqLen  || (MaxTargLen > 0 && MaxTargLen < pEntry->SeqLen)))
		{
		if(NumPreQuals != 0 && pQualTargs != NULL)
			{
			bFirst = false;
			pQualTarg = pQualTargs;
			for(PreQualIdx = 0; PreQualIdx < NumPreQuals; PreQualIdx+=1, pQualTarg+=1)
				{
				if(pEntry->EntryID == pQualTarg->TargEntryID)
					{
					bFirst = true;
					break;
					}
				}
			}
		else
			bFirst = true;
		if(bFirst)
			{
			HitLoci = (uint32_t)(TargLoci - pEntry->StartOfs);
			pEl2 = &pTarg[TargLoci];
			}
		}
	else
		bFirst = false;
	}

TargKMerLen = cPacBiokExtnKMerLen;  // when extending then using number of matching cPacBiokExtnKMerLen (currently 4bp) within extension whereby each matching 4-mer is at most 8 bp displaced due to assumed InDels
while(1)
	{
	if(!bFirst)
		{
		TargLoci = SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,PrevHitIdx);
		if((pEntry = MapChunkHit2Entry(TargLoci))==NULL)
			return(0);

		pEl2 = &pTarg[TargLoci];
		
		if((Cmp = CmpProbeTarg(pProbeSeq,pEl2,SeedCoreLen)) != 0)	// only 0 if still matching on the seed core
			return(0);

		if(pEntry->EntryID == SloughEntryID || MinTargLen > pEntry->SeqLen || (MaxTargLen > 0 && MaxTargLen < pEntry->SeqLen))
			{
			PrevHitIdx+=1;
			continue;
			}

		if(NumPreQuals != 0 && pQualTargs != NULL)
			{
			pQualTarg = pQualTargs;
			for(PreQualIdx = 0; PreQualIdx < NumPreQuals; PreQualIdx+=1, pQualTarg+=1)
				{
				if(pQualTarg->TargEntryID == pEntry->EntryID)
					break;
				}
			if(PreQualIdx == NumPreQuals)
				{
				PrevHitIdx+=1;
				continue;
				}
			}

		HitLoci = (uint32_t)(TargLoci - pEntry->StartOfs);
		}

	if(pEntry->SeqLen < (HitLoci + cPacBioSeedCoreExtn + 10) || MinTargLen > pEntry->SeqLen || (MaxTargLen > 0 && MaxTargLen < pEntry->SeqLen))
		{
		if(!bFirst)
			PrevHitIdx += 1;
		else
			bFirst = false;
		continue;
		}

	// see if the target sequence, using the initial matching core can be simply extended whilst still matching
	QAnchorIdx = SeedCoreLen;
	pQAnchor = &pProbeSeq[QAnchorIdx];
	pTAnchor = &pEl2[QAnchorIdx];
	for(MatchBaseLen = (int)SeedCoreLen; QAnchorIdx < (int)AcceptExactExtdCoreLen; QAnchorIdx++,pQAnchor++,pTAnchor++,MatchBaseLen++)
		{
		if((*pQAnchor & 0x07) != (*pTAnchor & 0x07))
			break;
		}
	if((uint32_t)MatchBaseLen < AcceptExactExtdCoreLen)			// no need to try extending if exactly matching extended core is at least cMaxPacBioSeedExtn (currently 50bp)
		{
		QuickScore = QuickScoreOverlap(cPacBioSeedCoreExtn - MatchBaseLen,pQAnchor,pTAnchor);  // QuickScore is a count of the number of bases covered by exactly matching, non-overlapping, 4-mers between probe and target in the extension 
	//		QuickScore = QuickScoreOverlap(min(10,MatchBaseLen),cPacBioSeedCoreExtn - MatchBaseLen,pQAnchor,pTAnchor);
		OverlapIdentity = ((MatchBaseLen + ((cPacBioSeedCoreExtn - MatchBaseLen) * QuickScore)/100)*100)/cPacBioSeedCoreExtn; // OverlapIdentity is accounting for the exactly matching prefix extended core

		if(OverlapIdentity < (int)AcceptExactExtdCoreLen)  
			{
			if(!bFirst)
				PrevHitIdx += 1;
			else
				bFirst = false;
			continue;
			}
		}

	*pTargEntryID = pEntry->EntryID;
	*pHitLoci = HitLoci;
	return(bFirst ? PrevHitIdx : PrevHitIdx+1);
	}
return(0);		// no more hits
}


int										// number (upto Lim) of non-canonical bases in pSeq 
CSfxArray::NumNonCanonicals(int Lim, // process for at most this many non-canonical bases
					uint32_t SeqLen,		// pSeq is of this max length
					etSeqBase *pSeq)	// pSeq to process for count of non-canonical
{
etSeqBase Base;
int Cnt;
if(Lim < 0 || SeqLen == 0 || pSeq == NULL)
	return(0);
Cnt = 0;
while(SeqLen--)
	{
	Base = *pSeq++ & 0x0f;
	if(Base >= eBaseN)
		{
		if(Base > eBaseN)
			break;
		Cnt++;
		if(Cnt >= Lim)
			break;
		}
	}
return(Cnt);
}

//
// LocateSfxHammings
// Returns minimum Hammings discovered from source squence to targeted assembly
//
int						// < 0 if errors, otherwise minimum hamming of probe to target
CSfxArray::LocateSfxHammings(int RHamm,				// restricted hammings limit
                         bool bSAHammings,				// if true then hammings on both sense and antisense required
						 int KMerLen,					// hammings for k-mers of this length
						 int SampleNth,					// sample every Nth k-mer
						 int SrcSeqLen,					// generate hammings over sequence of this length
						 uint32_t MinCoreDepth,		// initially explore cores to at least this depth and only continue exploring if depth would be no more than MaxCoreDepth
						 uint32_t MaxCoreDepth,		// explore cores to at most this depth
						 uint32_t EntryID,  				// KMer sequences are to be from this suffix entry
						 uint32_t EntryLoci,				// and starting at this loci
				 		 int IntraInterBoth,	    // 0: hammings over both intra (same sequence as probe K-mer drawn from) and inter (different sequences to that from which probe K-mer drawn), 1: Intra only, 2: Inter only
						 uint8_t *pHammings,				// returned hammings
						 int NumAllocdIdentNodes,		// memory has been allocated by caller for holding upto this many tsIdentNodes
						 tsIdentNode *pAllocsIdentNodes) // memory allocated by caller for holding tsIdentNodes
{
etSeqBase *pSrcSeq;			// hamming k-mers from this sequence
int Rslt;
int Ofs;
uint32_t EntryLen;
int Idx;
etSeqBase *pN;
etSeqBase *pHSrcSeq;
etSeqBase NonNSeq[501];
int NumNs;
int NIdxs[9];
int NIters;
int CurNIter;
int LowestHamming;

// quick sanity checks
if(EntryID == 0 || RHamm < 1 || KMerLen < 10 || KMerLen > 500 || SrcSeqLen < KMerLen || pHammings == NULL)
	return(eBSFerrInternal);
EntryLen = GetSeqLen(EntryID);
if(EntryLen == 0 || ((EntryLen - EntryLoci) < (uint32_t)SrcSeqLen))
	return(eBSFerrInternal);
pSrcSeq = GetPtrSeq(EntryID,EntryLoci);
if(pSrcSeq == NULL)
	return(eBSFerrInternal);

// iterate over each K-mer instance and save the Hamming for that instance
if(SampleNth <= 0)
	SampleNth = 1;
Rslt = 0;
for(Ofs = 0; Ofs <= (SrcSeqLen - KMerLen); Ofs += SampleNth,EntryLoci+=SampleNth,pSrcSeq += SampleNth)
	{
	// check on number of indeterminate bases in current K-mer, if between 1 and 4 then will iteratively substitute these with cannonical bases 
	// if more than 4 then will treat the current K-mer as if it is Hamming 0 away from any other K-mer in the target 
	NumNs = 0;
	NIters = 0;
	pN = pSrcSeq;
	for(Idx = 0; NumNs <= 4 && Idx < KMerLen; Idx++,pN++)
		if((*pN & 0x0f) >= eBaseN)
			NIdxs[NumNs++] = Idx;
	if(Idx == KMerLen && NumNs >= 1 && NumNs <= 4)
		{
		memcpy(NonNSeq,pSrcSeq,KMerLen);
		pHSrcSeq = 	NonNSeq;
		switch(NumNs) {
			case 1: NIters = 4; break;
			case 2: NIters = 4*4; break;
			case 3: NIters = 4*4*4; break;
			case 4: NIters = 4*4*4*4; break;
			default: NIters = 0; pHSrcSeq = pSrcSeq; break;
			}
		}
	else
		pHSrcSeq = pSrcSeq;

	if(NumNs <= 4)	// only process if not too many non-canonicals
		{
		CurNIter = 0;
		LowestHamming = 0x7f;
		do {
			if(NumNs)
				{
				for(Idx = 0; Idx < NumNs; Idx++)
					NonNSeq[NIdxs[Idx]] = (CurNIter >> (2 * Idx)) & 0x03;
				CurNIter += 1;
				}
			Rslt = LocateHamming(0,1,2*MinCoreDepth,10*MaxCoreDepth,bSAHammings,pHSrcSeq,KMerLen,EntryID,EntryLoci,IntraInterBoth,NumAllocdIdentNodes,pAllocsIdentNodes);
			if(Rslt > 1 && RHamm > 1)
				{
				Rslt = LocateHamming(2,2,2*MinCoreDepth,8*MaxCoreDepth,bSAHammings,pHSrcSeq,KMerLen,EntryID,EntryLoci,IntraInterBoth,NumAllocdIdentNodes,pAllocsIdentNodes);
				if(Rslt > 2 && RHamm > 2)
					{
					Rslt = LocateHamming(3,3,2*MinCoreDepth,4*MaxCoreDepth,bSAHammings,pHSrcSeq,KMerLen,EntryID,EntryLoci,IntraInterBoth,NumAllocdIdentNodes,pAllocsIdentNodes);
					if(Rslt > 3 && RHamm > 3)
						{
						Rslt = LocateHamming(4,4,MinCoreDepth,2*MaxCoreDepth,bSAHammings,pHSrcSeq,KMerLen,EntryID,EntryLoci,IntraInterBoth,NumAllocdIdentNodes,pAllocsIdentNodes);
						if(Rslt > 4 && RHamm > 4)
							Rslt = LocateHamming(5,RHamm,MinCoreDepth,MaxCoreDepth,bSAHammings,pHSrcSeq,KMerLen,EntryID,EntryLoci,IntraInterBoth,NumAllocdIdentNodes,pAllocsIdentNodes);
						}
					}
				}
			if(Rslt < 0)
				return(Rslt);
			if((Rslt & 0x7f) < LowestHamming)
				LowestHamming = Rslt & 0x7f;
			if(LowestHamming == 0)
				break;
			}
		while(NumNs > 0 && CurNIter < NIters);
		if(LowestHamming > 20)
			LowestHamming = 20;
		*pHammings = LowestHamming;
		}
	else
		{
		*pHammings = 0;
		Rslt = 0;
		}
	pHammings += SampleNth;
	}
return(Rslt < 0 ? Rslt : 0);
}

//
// LocateHammings
// Returns minimum Hammings discovered from source squence to targeted assembly
//
int						// < 0 if errors, otherwise minimum hamming of probe to target
CSfxArray::LocateHammings(int RHamm,					// restricted hammings limit
                         bool bSAHammings,				// if true then hammings on both sense and antisense required
						 int KMerLen,					// hammings for k-mers of this length
						 int SampleNth,					// sample every Nth k-mer
						 int SrcSeqLen,					// sequence containing k-mers is of this length
						 uint32_t MinCoreDepth,		// initially explore cores to at least this depth and only continue exploring if depth would be no more than MaxCoreDepth
						 uint32_t MaxCoreDepth,		// explore cores to at most this depth
						 etSeqBase *pSrcSeq,			// hamming k-mers from this sequence
						 uint32_t SeqEntry,  				// sequence was from this suffix entry (0 if probe not from same assembly)
						 uint32_t SeqLoci,				// and starting at this loci
				 		int IntraInterBoth,	    // 0: hammings over both intra (same sequence as probe K-mer drawn from) and inter (different sequences to that from which probe K-mer drawn), 1: Intra only, 2: Inter only
						 uint8_t *pHammings,				// returned hammings
						 int NumAllocdIdentNodes,		// memory has been allocated by caller for holding upto this many tsIdentNodes
						 tsIdentNode *pAllocsIdentNodes) // memory allocated by caller for holding tsIdentNodes
{
int Rslt;
int Ofs;
int Idx;
etSeqBase *pN;
etSeqBase *pHSrcSeq;
etSeqBase NonNSeq[501];
int NumNs;
int NIdxs[9];
int NIters;
int CurNIter;
int LowestHamming;

// quick sanity checks
if(RHamm < 1 || KMerLen < 10 || KMerLen > 500 || SrcSeqLen < KMerLen || pSrcSeq == NULL || pHammings == NULL)
	return(eBSFerrInternal);

// iterate over each K-mer instance and save the Hamming for that instance
if(SampleNth <= 0)
	SampleNth = 1;
Rslt = 0;
for(Ofs = 0; Ofs <= (SrcSeqLen - KMerLen); Ofs+=SampleNth,pSrcSeq+=SampleNth)
	{
	// check on number of indeterminate bases in current K-mer, if between 1 and 4 then will iteratively substitute these with cannonical bases 
	// if more than 4 then will treat the current K-mer as if it is Hamming 0 away from any other K-mer in the target 
	NumNs = 0;
	NIters = 0;
	pN = pSrcSeq;
	for(Idx = 0; NumNs <= 4 && Idx < KMerLen; Idx++,pN++)
		if((*pN & 0x0f) >= eBaseN)
			NIdxs[NumNs++] = Idx;
	if(Idx == KMerLen && NumNs >= 1 && NumNs <= 4)
		{
		memcpy(NonNSeq,pSrcSeq,KMerLen);
		pHSrcSeq = 	NonNSeq;
		switch(NumNs) {
			case 1: NIters = 4; break;
			case 2: NIters = 4*4; break;
			case 3: NIters = 4*4*4; break;
			case 4: NIters = 4*4*4*4; break;
			default: NIters = 0; pHSrcSeq = pSrcSeq; break;
			}
		}
	else
		pHSrcSeq = pSrcSeq;
	if(NumNs <= 4)	// only process if not too many non-canonicals
		{
		CurNIter = 0;
		LowestHamming = 0x7f;
		do {
			if(NumNs)
				{
				for(Idx = 0; Idx < NumNs; Idx++)
					NonNSeq[NIdxs[Idx]] = (CurNIter >> (2 * Idx)) & 0x03;
				CurNIter += 1;
				}
			Rslt = LocateHamming(0,1,2*MinCoreDepth,10*MaxCoreDepth,bSAHammings,pHSrcSeq,KMerLen,SeqEntry,SeqLoci+Ofs,IntraInterBoth,NumAllocdIdentNodes,pAllocsIdentNodes);
			if(Rslt > 1 && RHamm > 1)
				{
				Rslt = LocateHamming(2,2,2*MinCoreDepth,8*MaxCoreDepth,bSAHammings,pHSrcSeq,KMerLen,SeqEntry,SeqLoci+Ofs,IntraInterBoth,NumAllocdIdentNodes,pAllocsIdentNodes);
				if(Rslt > 2 && RHamm > 2)
					{
					Rslt = LocateHamming(3,3,2*MinCoreDepth,4*MaxCoreDepth,bSAHammings,pHSrcSeq,KMerLen,SeqEntry,SeqLoci+Ofs,IntraInterBoth,NumAllocdIdentNodes,pAllocsIdentNodes);
					if(Rslt > 3 && RHamm > 3)
						{
						Rslt = LocateHamming(4,4,MinCoreDepth,2*MaxCoreDepth,bSAHammings,pHSrcSeq,KMerLen,SeqEntry,SeqLoci+Ofs,IntraInterBoth,NumAllocdIdentNodes,pAllocsIdentNodes);
						if(Rslt > 4 && RHamm > 4)
							Rslt = LocateHamming(5,RHamm,MinCoreDepth,MaxCoreDepth,bSAHammings,pHSrcSeq,KMerLen,SeqEntry,SeqLoci+Ofs,IntraInterBoth,NumAllocdIdentNodes,pAllocsIdentNodes);
						}
					}
				}
			if(Rslt < 0)
				return(Rslt);
			if(Rslt < LowestHamming)
				LowestHamming = Rslt & 0x7f;
			if(LowestHamming == 0)
				break;
			}
		while(NumNs > 0 && CurNIter < NIters);

		*pHammings = LowestHamming;
		}
	else
		{
		*pHammings = 0;
		Rslt = 0;
		}
	pHammings += SampleNth;
	}
return(Rslt < 0 ? Rslt : 0);
}


//
// LocateHamming
// Returns minimum Hamming discovered from source probe squence to targeted assembly
//
int						// < 0 if errors, otherwise minimum hamming of probe to target
CSfxArray::LocateHamming(int RHammMin,				// process for Hammings of at least this minimum
						 int RHammMax,					// process for Hammings upto this limit
						 uint32_t MinCoreDepth,			// initially explore cores to at least this depth and only continue exploring if depth would be no more than MaxCoreDepth
						 uint32_t MaxCoreDepth,			// explore cores to at most this depth
                         bool bSAHammings,				// if true then hammings on both sense and antisense required
						 etSeqBase *pProbeSeq,			// probe sequence
						 uint32_t ProbeLen,				// probe length
 						 uint32_t ProbeEntry,  			// probe was from this suffix entry (0 if probe not from same assembly)
						 uint32_t ProbeLoci,				// and starting at this loci
				 		 int IntraInterBoth,	    // 0: hammings over both intra (same sequence as probe K-mer drawn from) and inter (different sequences to that from which probe K-mer drawn), 1: Intra only, 2: Inter only
						 int NumAllocdIdentNodes,		// memory has been allocated by caller for holding upto this many tsIdentNodes
						 tsIdentNode *pAllocsIdentNodes) // memory allocated by caller for holding tsIdentNodes
{
uint8_t RevCplSeq1Kbp[0x03ff];

uint8_t *pRevCplSeq;

bool bDoAntisense;
int CurProbeSegOfs;
int CurProbeSegLen;
int CoreLen;
int SegIdx;
uint32_t IterCnt;
int CurHamming;
int MaxSegIdx;

etSeqBase *pProbeBase;
etSeqBase *pTargBase;

int64_t TargIdx;

int CurSubCnt;
int PatIdx;
int64_t TargSeqLeftIdx;
uint32_t ProbeSeqLeftIdx;
int TargMatchLen;
int64_t PutativeTargLoci;

uint8_t ProbeBase;
uint8_t TargBase;

uint32_t HitLoci;
tsSfxEntry *pEntry;
int Cmp;

int64_t LastTargIdx;
uint32_t NumCopies;

etSeqBase *pTarg;			// target sequence
void *pSfxArray;			// target sequence suffix array
int64_t SfxLen;				// number of suffixs in pSfxArray


// ensure suffix loaded for iteration!
if(m_pSfxBlock == NULL || m_pSfxBlock->ConcatSeqLen == 0)
	return(-1);

pTarg = (etSeqBase *)&m_pSfxBlock->SeqSuffix[0];
SfxLen = m_pSfxBlock->ConcatSeqLen;
pSfxArray = (void *)&m_pSfxBlock->SeqSuffix[SfxLen];
if(bSAHammings)
	{
	if(ProbeLen >= sizeof(RevCplSeq1Kbp))
		{
		if((pRevCplSeq = (uint8_t *)malloc(ProbeLen+1))==NULL)
			return(eBSFerrMem);
		}
	else
		pRevCplSeq = RevCplSeq1Kbp;
	}
else
	pRevCplSeq = NULL;
if(RHammMin == 0)
	{
	// first try for a hamming of 0! Could strike it lucky!
	int64_t PrevHitIdx = 0;
	uint32_t TargEntryID;
	while((PrevHitIdx = IterateExacts(pProbeSeq,ProbeLen,PrevHitIdx,&TargEntryID,&HitLoci)) > 0)
		{
		if(ProbeEntry == TargEntryID && ProbeLoci == HitLoci)  // not interested in self hits!
			continue;

		if(ProbeEntry != 0 && IntraInterBoth != 0)
			{
			if((IntraInterBoth == 1 && ProbeEntry != TargEntryID) ||
				(IntraInterBoth == 2 && ProbeEntry == TargEntryID))
				continue;
			}

		if(pRevCplSeq != NULL && (ProbeLen >= sizeof(RevCplSeq1Kbp)))
			free(pRevCplSeq);
		return(0);
		}

	if(bSAHammings)
		{
		memmove(pRevCplSeq,pProbeSeq,ProbeLen);
		CSeqTrans::ReverseComplement(ProbeLen,pRevCplSeq);
		PrevHitIdx = 0;
		while((PrevHitIdx = IterateExacts(pRevCplSeq,ProbeLen,PrevHitIdx,&TargEntryID,&HitLoci)) > 0)
			{
			if(ProbeEntry != 0 && IntraInterBoth != 0)
				{
				if((IntraInterBoth == 1 && ProbeEntry != TargEntryID) ||
					(IntraInterBoth == 2 && ProbeEntry == TargEntryID))
					continue;
				}
	
			if(pRevCplSeq != NULL && (ProbeLen >= sizeof(RevCplSeq1Kbp)))
				free(pRevCplSeq);
			return(0);
			}
		}
	}

if(RHammMax == 0)		// if caller was only interested in hammings of 0 then any other hammings must be of at least 1
	{
	if(pRevCplSeq != NULL && (ProbeLen >= sizeof(RevCplSeq1Kbp)))
		free(pRevCplSeq);
	return(1);
	}

// expected Hamming is at least 1
CoreLen = ProbeLen/(RHammMax + 1);
CurHamming = ProbeLen/CoreLen;
bDoAntisense = false;
do
	{
	CurProbeSegOfs = 0;
	CurProbeSegLen = CoreLen;
	MaxSegIdx = ProbeLen/CoreLen;

	for(SegIdx = 0; SegIdx < MaxSegIdx; SegIdx++, CurProbeSegOfs += CoreLen)
		{
		TargIdx = LocateFirstExact(&pProbeSeq[CurProbeSegOfs],CoreLen,pTarg,m_pSfxBlock->SfxElSize,pSfxArray,0,0,SfxLen-1);

		if(TargIdx == 0) // if no segment match then onto next segment
			continue;
		TargIdx -= 1;
		IterCnt = 0;
		NumCopies = 0;
		while(CurHamming > RHammMin && IterCnt < MaxCoreDepth)
			{
			if(IterCnt++) {				// if iterating subsequent (to LocateFirstExact) targets
				// ensure not about to iterate past end of suffix array!
				if((TargIdx + 1) >= (int64_t)m_pSfxBlock->ConcatSeqLen || (SfxOfsToLoci((int32_t)m_pSfxBlock->SfxElSize,pSfxArray,TargIdx+1) + CoreLen) > (int64_t)m_pSfxBlock->ConcatSeqLen)
					break;

				if(IterCnt == MinCoreDepth  && !NumCopies)
					{
					// check how many more exact copies there are of the current probe subsequence, if too many then don't bother exploring these
					LastTargIdx = LocateLastExact(&pProbeSeq[CurProbeSegOfs],CoreLen,pTarg,m_pSfxBlock->SfxElSize,pSfxArray,0,TargIdx-1,SfxLen-1);
					NumCopies = LastTargIdx > 0 ? (uint32_t)(1 + LastTargIdx - TargIdx) : 0;
					if(MaxCoreDepth && NumCopies > MaxCoreDepth)		// only checking at the MinCoreDepth iteration allows a little slack
						break;										// try next core segment
					}

				pTargBase = &pTarg[SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,TargIdx+1)]; // then check that target core is still matching
				pProbeBase = &pProbeSeq[CurProbeSegOfs];
					{
					int Ofs;
					uint8_t El1;
					uint8_t El2;
					etSeqBase *pEl1= pProbeBase;
					etSeqBase *pEl2 = pTargBase;
					Cmp = 0;
					for(Ofs=0; Ofs < CoreLen; Ofs++)
						{
						El2 = *pEl2++ & 0x0f;
						if(El2 == eBaseEOS || El2 == eBaseN)
							{
							Cmp = -1;
							break;
							}
						El1 = *pEl1++ & 0x0f;
						if(El1 > El2)
							{
							Cmp = 1;
							break;
							}
						if(El1 < El2 || El1 == eBaseN)
							{
							Cmp = -1;
							break;
							}
						}

					if(Cmp != 0)				// will be non-zero if target no longer matches
						break;					// try next core segment
					TargIdx += 1;
					}
				}


			if(SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,TargIdx) < (uint32_t)CurProbeSegOfs)
				continue;

			TargMatchLen = ProbeLen;
			TargSeqLeftIdx = SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,TargIdx) - CurProbeSegOfs;
			ProbeSeqLeftIdx = 0;
			PutativeTargLoci = TargSeqLeftIdx;

			// ensure comparisons are within start/end range of target sequence/assembly
			if(!TargMatchLen || ((uint64_t)TargSeqLeftIdx + (uint32_t)TargMatchLen) > m_pSfxBlock->ConcatSeqLen)
				continue;

			// ensure this target segment not already processed as an extension of a previous segment
			if(SegIdx > 0)
				{
				int DupIdx;
				int CurDupSegLen = 0;
				int CurDupSegOfs = 0;
				pTargBase = &pTarg[TargSeqLeftIdx];
				pProbeBase = pProbeSeq;
				for(DupIdx = 0; DupIdx < SegIdx; DupIdx++, CurDupSegOfs += CoreLen, pProbeBase += CoreLen, pTargBase+=CoreLen)
					{
					Cmp = CmpProbeTarg(pProbeBase,pTargBase,CoreLen);
					if(Cmp==0)
						break;
					}
				if(DupIdx < SegIdx)
					continue;
				}

			// now do the matching
			pTargBase = &pTarg[TargSeqLeftIdx]; // then check that target is still matching
			pProbeBase = pProbeSeq;
			CurSubCnt = 0;
			for(PatIdx = 0; PatIdx < TargMatchLen; PatIdx++,pTargBase++,pProbeBase++)
				{
				TargBase = *pTargBase & 0x0f;
				ProbeBase = *pProbeBase & 0x0f;
				if(TargBase == eBaseEOS)		// mustn't match across entry sequences
					break;
				if(ProbeBase == eBaseN || ProbeBase == TargBase) // treating probe indeterminate bases as if they would match thus reducing the Hamming
					continue;

				// execution here only if mismatch
				if(++CurSubCnt > CurHamming)
					break;
				}
			if(PatIdx != TargMatchLen)
				continue;

			// check for self matches

			if(CurSubCnt == 0 && !bDoAntisense && ProbeEntry > 0)
				{
				pEntry = MapChunkHit2Entry(TargSeqLeftIdx);

				if(ProbeEntry == pEntry->EntryID)
					{
					HitLoci = (uint32_t)(TargSeqLeftIdx - pEntry->StartOfs);
					if(ProbeLoci == HitLoci)
						continue;
					}
				if(ProbeEntry != 0 && IntraInterBoth != 0)
					{
					if((IntraInterBoth == 1 && ProbeEntry != pEntry->EntryID) ||
						(IntraInterBoth == 2 && ProbeEntry == pEntry->EntryID))
						continue;
					}
				}

			if(CurSubCnt < CurHamming)
				CurHamming = CurSubCnt;
			}

		if(CurHamming <= RHammMin)
			break;
		}

	if((CurHamming > RHammMin) && bSAHammings)
		{
		memmove(pRevCplSeq,pProbeSeq,ProbeLen);
		pProbeSeq = pRevCplSeq;
		CSeqTrans::ReverseComplement(ProbeLen,pProbeSeq);
		bSAHammings = false;
		ProbeEntry = 0;
		bDoAntisense = true;
		}
	else
		bDoAntisense = false;
	}
while((CurHamming > RHammMin) && bDoAntisense);
if(pRevCplSeq != NULL && (ProbeLen >= sizeof(RevCplSeq1Kbp)))
	free(pRevCplSeq);
return(CurHamming);
}

//
// LocateApproxUniques
// Locates probe to target matches with mismatches
// If m_bBisulfite is set true then matches using bisulfite rules
//
int						// < 0 if errors, 0 if no matches, 1 if a unique match, 2 if multiple matches
CSfxArray::LocateApproxUniques(int AccumReadHits,		// how many reads have already been matched, must be 0 or 1
						 etSeqBase *pProbeSeq,			// probe
 						 etSeqBase *pPatternSeq,		// contains pattern to match with, etBaseN represents wildcard bases to match against any in the target
														// will be updated on return with etBaseN's changed to actual subsitution bases  - NULL if all bases can be wildcarded
						 uint32_t ProbeLen,				// probe, and also pattern, length
						 int ExpMismatches,				// expected wildcard matches
 						 uint32_t *pTargEntryID,			// if unique match then where to return suffix entry (chromosome) matched on
						 uint32_t *pHitLoci,				// if unique match then where to return loci
						 int CurMaxIter)				// max allowed iterations per subsegmented sequence when matching that subsegment
{
int CurProbeSegOfs;
int CurProbeSegLen;
int SegIdx;
int IterCnt;

etSeqBase *pProbeBase;
etSeqBase *pPatternBase;
etSeqBase *pTargBase;

int64_t TargIdx;
int64_t LastTargIdx;
uint32_t NumCopies;

int CurSubCnt;
int PatIdx;
int64_t TargSeqLeftIdx;
int ProbeSeqLeftIdx;
int TargMatchLen;
int64_t PutativeTargLoci;
int CurInstances;

etSeqBase BisBase;
uint8_t ProbeBase;
uint8_t TargBase;

int Cmp;

etSeqBase *pTarg;			// target sequence
void *pSfxArray;			// target sequence suffix array
int64_t SfxLen;				// number of suffixs in pSfxArray
tsSfxEntry *pEntry;
*pHitLoci = 0;
*pTargEntryID = 0;

// ensure suffix loaded for iteration!
if(m_pSfxBlock == NULL || m_pSfxBlock->ConcatSeqLen == 0)
	return(0);

pTarg = (etSeqBase *)&m_pSfxBlock->SeqSuffix[0];
pSfxArray = (void *)&m_pSfxBlock->SeqSuffix[m_pSfxBlock->ConcatSeqLen];
SfxLen = m_pSfxBlock->ConcatSeqLen;

CurProbeSegOfs = 0;
CurInstances = AccumReadHits;

for(SegIdx = 0; SegIdx <= ExpMismatches; SegIdx++,CurProbeSegOfs += CurProbeSegLen)
	{
	CurProbeSegLen = (ProbeLen - CurProbeSegOfs) / (1 + (ExpMismatches - SegIdx));

	TargIdx = LocateFirstExact(&pProbeSeq[CurProbeSegOfs],CurProbeSegLen,pTarg,m_pSfxBlock->SfxElSize,pSfxArray,0,0,SfxLen-1);
	if(TargIdx == 0) // if no segment match then onto next segment
		continue;
	TargIdx -= 1;
	IterCnt = 0;
	NumCopies = 0;
	while((!CurMaxIter || IterCnt < CurMaxIter) && CurInstances < 2)
		{
		if(IterCnt++) {				// if iterating subsequent (to LocateFirstExact) targets
			// ensure not about to iterate past end of suffix array!
			if((TargIdx + 1) >= (int64_t)m_pSfxBlock->ConcatSeqLen || (SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,TargIdx+1) + CurProbeSegLen) >  (int64_t)m_pSfxBlock->ConcatSeqLen)
				break;

			if(IterCnt == 100 && !NumCopies)
				{
				// check how many more exact copies there are of the current probe subsequence, if too many then don't bother exploring these
				LastTargIdx = LocateLastExact(&pProbeSeq[CurProbeSegOfs],CurProbeSegLen,pTarg,m_pSfxBlock->SfxElSize,pSfxArray,0,TargIdx-1,SfxLen-1);
				NumCopies = LastTargIdx > 0 ? (uint32_t)(1 + LastTargIdx - TargIdx) : 0;
				if(CurMaxIter && NumCopies > (uint32_t)CurMaxIter)		// only checking at the 100th iteration allows a little slack
					break;										// try next core segment
				}

			pTargBase = &pTarg[SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,TargIdx+1)]; // then check that target is still matching
			pProbeBase = &pProbeSeq[CurProbeSegOfs];
			if(m_bBisulfite)
				Cmp = BSCmpProbeTarg(pProbeBase,pTargBase,CurProbeSegLen);
			else
				Cmp = CmpProbeTarg(pProbeBase,pTargBase,CurProbeSegLen);
			if(Cmp != 0)				// != 0 if target no longer matches
				break;
			TargIdx += 1;
			}

		TargMatchLen = ProbeLen;
		TargSeqLeftIdx = SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,TargIdx) - CurProbeSegOfs;
		ProbeSeqLeftIdx = 0;
		PutativeTargLoci = TargSeqLeftIdx;

		// ensure comparisons are within start/end range of target sequence/assembly
		if(!TargMatchLen || TargSeqLeftIdx < 0 || (TargSeqLeftIdx + TargMatchLen) > (int) m_pSfxBlock->ConcatSeqLen)
			continue;

		// ensure this target segment not already processed as an extension of a previous segment
		if(SegIdx > 0)
			{
			int DupIdx;
			int CurDupSegLen;
			int CurDupSegOfs = 0;
			pTargBase = &pTarg[TargSeqLeftIdx];
			pProbeBase = pProbeSeq;
			for(DupIdx = 0; DupIdx < SegIdx; DupIdx++,CurDupSegOfs += CurDupSegLen,pProbeBase += CurDupSegLen,pTargBase+=CurDupSegLen)
				{
				CurDupSegLen = (ProbeLen - CurDupSegOfs) / (1 + (ExpMismatches - DupIdx));
				if(m_bBisulfite)
					Cmp = BSCmpProbeTarg(pProbeBase,pTargBase,CurDupSegLen);
				else
					Cmp = CmpProbeTarg(pProbeBase,pTargBase,CurDupSegLen);
				if(Cmp==0)
					break;
				}
			if(DupIdx < SegIdx)
				continue;
			}

		// now do the matching
		pTargBase = &pTarg[TargSeqLeftIdx]; // then check that target is still matching
		pProbeBase = pProbeSeq;
		pPatternBase = pPatternSeq;
		CurSubCnt = 0;
		if(m_bBisulfite)
			BisBase = GetBisBase(TargMatchLen,pTargBase,pProbeBase);
		for(PatIdx = 0; PatIdx < TargMatchLen; PatIdx++,pTargBase++,pProbeBase++)
			{
			if(pPatternBase != NULL)
				pPatternBase += 1;
			TargBase = *pTargBase & 0x0f;
			ProbeBase = *pProbeBase & 0x0f;
			if(TargBase == eBaseEOS)		// mustn't match across entry sequences
				break;
			if(ProbeBase == TargBase)
				continue;

			if(m_bBisulfite)
				{
				if(TargBase == eBaseC || TargBase == eBaseG)
					{
					switch(BisBase) {
						case eBaseA:		// only allow A
							if(ProbeBase == eBaseA && TargBase == eBaseG)
								continue;
							break; // mismatch
						case eBaseT:		// only allow T
							if(ProbeBase == eBaseT && TargBase == eBaseC)
								continue;
							break; // mismatch
						}
					}
				}
			// execution here only if mismatch
			if(pPatternBase != NULL && pPatternBase[-1] != eBaseN)
				break;
			if(++CurSubCnt > ExpMismatches)
				break;
			}
		if(PatIdx != TargMatchLen)
			continue;

		if(ExpMismatches != CurSubCnt)
			continue;

		if(++CurInstances == 1)		// if 1st instance then record target entry loci
			{
			pEntry = MapChunkHit2Entry(PutativeTargLoci);
			if(pEntry == NULL)	// should never happen!
				continue;
			*pTargEntryID = pEntry->EntryID;
			*pHitLoci = (uint32_t)(PutativeTargLoci - pEntry->StartOfs);
			}
		}
	if(CurInstances > 1)
		break;
	}
return(CurInstances);
}


int						// < 0 if errors, too many matches, or sequence locally non-unique to ProbeEntry , 0 if no matches, > 0 total number of matches  (see genzygosity.cpp)
CSfxArray::LocateAllNearMatches(int MaxTotMM,			// max number of mismatches allowed
						 int CoreLen,					// core window length
						 int CoreDelta,					// core window offset increment (1..n)
						 int MaxNumCoreSlides,			// max number of times to slide core on each strand
						 int MaxMatches,				// if more than 0 and encountering more matches than this limit then return match count -1
 					 	 int NumEntries,				// number of entries in pEntryMatch
						 uint32_t *pEntryMatches,			// return number of matches to each entry in this array
						 uint32_t ProbeEntry,				// probe was from this entry
						 uint32_t ProbeOfs,				// starting at this offset
						 int ProbeLen,					// probe length
						 etSeqBase *pProbeSeq,			// probe
						 int CurMaxIter,				// max allowed iterations per subsegmented sequence when matching that subsegment
						 int NumAllocdIdentNodes,		// memory has been allocated by caller for holding upto this many tsIdentNodes
						 tsIdentNode *pAllocsIdentNodes) // memory allocated by caller for holding tsIdentNodes
{
int CurCoreSegOfs;				// current core segment relative start
int IterCnt;					// count iterator for current segment target matches
int CurMMCnt;					// total number of mismatches for current target sequence being processed
int CoreMMCnt;					// total number of mismatches for current target core being processed
int TotalHitInstances;			// total number of matches
etSeqBase *pProbeBase;
etSeqBase *pTargBase;
char CurStrand;
int64_t TargIdx;

tsIdentNode *pHashArray[cHashEntries+1];		// hash array holding ptrs to identifier nodes
tsIdentNode *pIdentNodes = pAllocsIdentNodes;	// identifier nodes
tsIdentNode *pCurIdentNode;
tsIdentNode *pNewIdentNode;
int CurNumIdentNodes;
int Hash;
uint32_t HitLoci;
uint32_t TargSeqID;
uint32_t NumTargSeqProc;

bool bFirstIter;				// set false after the first subsequence core returned by LocateFirstExact has been processed

uint32_t PatIdx;
int64_t TargSeqLeftIdx;
uint32_t ProbeSeqLeftIdx;
int TargMatchLen;

int64_t LastTargIdx;
uint32_t NumCopies;

etSeqBase BisBase;
uint8_t ProbeBase;
uint8_t TargBase;

int Cmp;
int CurNumCoreSlides;
int CurCoreDelta;
tsHitLoci *pCurHit;

etSeqBase *pTarg;			// target sequence
void *pSfxArray;			// target sequence suffix array
int64_t SfxLen;				// number of suffixs in pSfxArray
tsSfxEntry *pEntry;
BisBase=eBaseN;

// ensure suffix array block loaded for iteration!
if(m_pSfxBlock == NULL || m_pSfxBlock->ConcatSeqLen == 0)
	return(eBSFerrInternal);

pTarg = (etSeqBase *)&m_pSfxBlock->SeqSuffix[0];
pSfxArray = (void *)&m_pSfxBlock->SeqSuffix[m_pSfxBlock->ConcatSeqLen];
SfxLen = m_pSfxBlock->ConcatSeqLen;

pCurHit = NULL;
CurStrand = '+';			// initially treat probe as not reverse complemented, e.g. '+' matches

// counts of matches to each suffix array entry are being recorded so initialise these counts to be 0
if(pEntryMatches != NULL && NumEntries >= 1)
	memset(pEntryMatches,0,sizeof(uint32_t) * NumEntries);

NumTargSeqProc = 0;
CurNumIdentNodes = 0;
CurNumCoreSlides = 0;
TotalHitInstances = 0;
do
	{
	CurCoreDelta = CoreDelta;
	CurNumCoreSlides = 0;
	memset(pHashArray,0,sizeof(pHashArray));
	CurNumIdentNodes = 0;
	for(CurCoreSegOfs = 0;
		CurNumCoreSlides < MaxNumCoreSlides &&
		CurCoreSegOfs <= (ProbeLen - CoreLen) &&
		CurCoreDelta > CoreLen/3 &&
		CurNumIdentNodes < NumAllocdIdentNodes;	// can only allow up to cMaxNumIdentNodes to be saved
	    CurNumCoreSlides += 1,
	    CurCoreSegOfs += CurCoreDelta)
		{
		if(CurNumCoreSlides >= MaxNumCoreSlides)
			break;
		if((CurCoreSegOfs + CoreLen + CurCoreDelta) > ProbeLen)
			CurCoreDelta = ProbeLen - (CurCoreSegOfs + CoreLen);

		TargIdx = LocateFirstExact(&pProbeSeq[CurCoreSegOfs],CoreLen,pTarg,m_pSfxBlock->SfxElSize,pSfxArray,0,0,SfxLen-1);
		if(TargIdx == 0)        // 0 if no core segment matches
			continue;			// try for match on next core segment after shifting core to right

		TargIdx -= 1;
		IterCnt = 0;
		NumCopies = 0;
		bFirstIter = true;		// set false after the first subsequence core returned by LocateFirstExact has been processed
		while(!CurMaxIter || IterCnt < CurMaxIter)
			{
			if(CurNumIdentNodes >= NumAllocdIdentNodes)
				break;

			if(!bFirstIter) {

				// ensure not about to iterate past end of suffix array!
				if((TargIdx + 1) >= (int64_t)m_pSfxBlock->ConcatSeqLen || (SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,TargIdx+1) + CoreLen) >  (int64_t)m_pSfxBlock->ConcatSeqLen)
					break;

				if(IterCnt == 100 && !NumCopies)
					{
					// check how many more exact copies there are of the current probe subsequence, if too many then don't bother exploring these
					LastTargIdx = LocateLastExact(&pProbeSeq[CurCoreSegOfs],CoreLen,pTarg,m_pSfxBlock->SfxElSize,pSfxArray,0,TargIdx-1,SfxLen-1);
					NumCopies = LastTargIdx > 0 ? (uint32_t)(1 + LastTargIdx - TargIdx) : 0;
					if(CurMaxIter && NumCopies > (uint32_t)CurMaxIter)		// only checking at the 100th iteration allows a little slack
						break;										// try next core segment
					}


				// check that this new putative core is still matching
				pTargBase = &pTarg[SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,TargIdx+1)];
				pProbeBase = &pProbeSeq[CurCoreSegOfs];
				if(m_bBisulfite)
					Cmp = BSCmpProbeTarg(pProbeBase,pTargBase,CoreLen);
				else
					{
					int Ofs;
					uint8_t El1;
					uint8_t El2;
					etSeqBase *pEl1= pProbeBase;
					etSeqBase *pEl2 = pTargBase;
					Cmp = 0;
					for(Ofs=0; Ofs < CoreLen; Ofs++)
						{
						El2 = *pEl2++ & 0x0f;
						if(El2 == eBaseEOS)
							{
							Cmp = -1;
							break;
							}
						El1 = *pEl1++ & 0x0f;
						if(El1 > El2)
							{
							Cmp = 1;
							break;
							}
						if(El1 < El2)
							{
							Cmp = -1;
							break;
							}
						}
					}

				if(Cmp != 0)				// will be non-zero if target no longer matches
					break;					// try next core segment
				TargIdx += 1;
				}

			bFirstIter = false;
			if(SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,TargIdx) < (uint32_t)CurCoreSegOfs)
				continue;

			TargMatchLen = ProbeLen;
			TargSeqLeftIdx = SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,TargIdx) - CurCoreSegOfs;
			ProbeSeqLeftIdx = 0;

			// ensure comparisons are still within start/end range of target sequence/assembly
			if(!TargMatchLen || (TargSeqLeftIdx + TargMatchLen) > (int64_t)m_pSfxBlock->ConcatSeqLen)
				continue;

			// check if target already processed
			TargSeqID = (uint32_t)(1 + SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,TargIdx) - (uint32_t)CurCoreSegOfs);
			Hash = (TargSeqID & cHashEntries);
			if((pCurIdentNode = pHashArray[Hash])==NULL)
				{
				pCurIdentNode = pHashArray[Hash] = &pIdentNodes[CurNumIdentNodes++];
				pCurIdentNode->TargSeqID = TargSeqID;
				pCurIdentNode->pNxt = NULL;
				}
			else
				{
				while(pCurIdentNode != NULL && pCurIdentNode->TargSeqID != TargSeqID)
					pCurIdentNode = pCurIdentNode->pNxt;
				if(pCurIdentNode != NULL)							// will be non-null if have already processed this sequence
					continue;
				pNewIdentNode = &pIdentNodes[CurNumIdentNodes++];
				pNewIdentNode->TargSeqID = TargSeqID;
				pNewIdentNode->pNxt = pHashArray[Hash];
				pHashArray[Hash] = pNewIdentNode;
				}
			NumTargSeqProc += 1;
			IterCnt += 1;

			// now do the matching allowing for missmatches
				{
				pTargBase = &pTarg[TargSeqLeftIdx];
				pProbeBase = pProbeSeq;
				CurMMCnt = 0;
				CoreMMCnt = 0;
				bool bPairMM = false;
				for(PatIdx = 0; PatIdx < (uint32_t)TargMatchLen; PatIdx++,pTargBase++,pProbeBase++)
					{
					TargBase = *pTargBase & 0x0f;
					ProbeBase = *pProbeBase & 0x0f;
					if(TargBase == eBaseEOS)		// mustn't match across entry sequences
						break;

					if(ProbeBase == TargBase)
						continue;

					// execution here only if mismatch
					if(++CurMMCnt > MaxTotMM)
						break;
					}
				if(PatIdx != TargMatchLen)
					continue;

				// check for self matches
				pEntry = MapChunkHit2Entry(TargSeqLeftIdx);

				if(ProbeEntry == pEntry->EntryID)
					{
					HitLoci = (uint32_t)(TargSeqLeftIdx - pEntry->StartOfs);
					if(ProbeOfs == HitLoci)
						continue;
					// not a self match, so subsequence not unique to this sfx entry
					if(CurStrand == '-')									// restore probe sequence if had started processing '-' strand
						CSeqTrans::ReverseComplement(ProbeLen,pProbeSeq);
					return(-1);
					}

				// match onto another sfx entry, check if too many match instances
				if(MaxMatches > 0 && TotalHitInstances == MaxMatches)
					{
					if(CurStrand == '-')									// restore probe sequence if had started processing '-' strand
						CSeqTrans::ReverseComplement(ProbeLen,pProbeSeq);
					return(-1);
					}

				// have accepted this match so record which entry the match was on
				if(pEntry !=NULL && pEntry->EntryID <= (uint32_t)NumEntries)
					pEntryMatches[pEntry->EntryID-1] += 1;
				TotalHitInstances += 1;
				}
			}
		}
	if(CurStrand == '+')
		{
		CSeqTrans::ReverseComplement(ProbeLen,pProbeSeq);
		CurStrand = '-';
		NumTargSeqProc = 0;
		}
	else
		break;
	}
while(1);

if(CurStrand == '-')									// restore probe sequence if had started processing '-' strand
	CSeqTrans::ReverseComplement(ProbeLen,pProbeSeq);

return(TotalHitInstances);				// no errors, return total number of hit instances to all entries
}

int										// < 0 if errors, 0 if no matches to any other chroms allowing up to MaxTotMM, 1 if matches to any other chrom
CSfxArray::MatchesOtherChroms( int ProbeChromID,	// probe is from this chrom
				    int MaxTotMM,		// allow for at most this number of missmatches
				    int ProbeLen,		// probe sequence is this length
					etSeqBase *pProbe)	// check for matches from this probe sequence
{
int64_t PrevHitIdx;
int64_t NxtHitIdx;
uint32_t TargHitID;
uint32_t TargHitLoci;
int CoreLen;
int CurCoreSegOfs;
int NumMM;
int Idx;
etSeqBase *pStartCoreBase;
etSeqBase *pHitCoreSeq;
etSeqBase *pTarg;			// target sequence
void *pSfxArray;			// target sequence suffix array
int64_t SfxLen;				// number of suffixs in pSfxArray
tsSfxEntry *pEntry;

pTarg = (etSeqBase *)&m_pSfxBlock->SeqSuffix[0];
pSfxArray = (void *)&m_pSfxBlock->SeqSuffix[m_pSfxBlock->ConcatSeqLen];
SfxLen = m_pSfxBlock->ConcatSeqLen;
CoreLen = ProbeLen / (1+MaxTotMM);

if(CoreLen < 8)				// have to have a minimum core otherwise may as well do a linear search!
	CoreLen = 8;				

for(CurCoreSegOfs = 0; CurCoreSegOfs < (ProbeLen - CoreLen); CurCoreSegOfs += CoreLen)
	{
	if((CurCoreSegOfs + CoreLen) > ProbeLen)
		CurCoreSegOfs = ProbeLen - CoreLen;
	pStartCoreBase = &pProbe[CurCoreSegOfs];
	PrevHitIdx = 0;
	NxtHitIdx = 0;
	while((NxtHitIdx = IterateExacts(pStartCoreBase,CoreLen,PrevHitIdx,&TargHitID,&TargHitLoci)) > 0)
		{
		PrevHitIdx = NxtHitIdx;
		if(TargHitID == ProbeChromID)			// only interested in hits onto other chroms
			continue;
		if(TargHitLoci < (uint32_t)CurCoreSegOfs)					// must be able to extend core left
			continue;

		pEntry = &m_pEntriesBlock->Entries[TargHitID-1];
		if((TargHitLoci + (uint32_t)ProbeLen - (uint32_t)CurCoreSegOfs) > pEntry->SeqLen)					// must be able to extend core right
			continue;

		// check flanks and total number of missmatches, if <= MaxTotMM then return as having matched at least one other chrom
		pHitCoreSeq = &m_pSfxBlock->SeqSuffix[pEntry->StartOfs + TargHitLoci];
		pTarg = pHitCoreSeq - CurCoreSegOfs;
		pHitCoreSeq = pProbe;
		NumMM = 0;
		for(Idx = 0; NumMM > MaxTotMM && Idx < ProbeLen; Idx++, pTarg++, pHitCoreSeq++)
			{
			if((*pTarg & 0x0f) == (*pHitCoreSeq & 0x0f))
				continue;
			NumMM += 1;
			}
		if(NumMM <= MaxTotMM)
			return(1);
		}
	}

return(0);
}


int										// < 0 if errors, 0 if no matches to any other chroms not in pProbeChromIDs[] allowing up to MaxTotMM, 1 if matches to any other chrom
CSfxArray::MatchesOtherChroms(int NumProbeIDs, // number of chroms in pProbeChromIDs
					int *pProbeChromIDs,	// probe is from one of these chroms
				    int MaxTotMM,		// allow for at most this number of missmatches
				    int ProbeLen,		// probe sequence is this length
					etSeqBase *pProbe)	// check for matches from this probe sequence
{
int ProbeIdx;
int *pProbeID;

int64_t PrevHitIdx;
int64_t NxtHitIdx;
uint32_t TargHitID;
uint32_t TargHitLoci;
int CoreLen;
int CurCoreSegOfs;
int NumMM;
int Idx;
etSeqBase *pStartCoreBase;
etSeqBase *pHitCoreSeq;
etSeqBase *pTarg;			// target sequence
void *pSfxArray;			// target sequence suffix array
int64_t SfxLen;				// number of suffixs in pSfxArray
tsSfxEntry *pEntry;

pTarg = (etSeqBase *)&m_pSfxBlock->SeqSuffix[0];
pSfxArray = (void *)&m_pSfxBlock->SeqSuffix[m_pSfxBlock->ConcatSeqLen];
SfxLen = m_pSfxBlock->ConcatSeqLen;
CoreLen = ProbeLen / (1+MaxTotMM);

if(CoreLen < 8)				// have to have a minimum core otherwise may as well do a linear search!
	CoreLen = 8;				

for(CurCoreSegOfs = 0; CurCoreSegOfs < (ProbeLen - CoreLen); CurCoreSegOfs += CoreLen)
	{
	if((CurCoreSegOfs + CoreLen) > ProbeLen)
		CurCoreSegOfs = ProbeLen - CoreLen;
	pStartCoreBase = &pProbe[CurCoreSegOfs];
	PrevHitIdx = 0;
	NxtHitIdx = 0;
	while((NxtHitIdx = IterateExacts(pStartCoreBase,CoreLen,PrevHitIdx,&TargHitID,&TargHitLoci)) > 0)
		{
		PrevHitIdx = NxtHitIdx;
		pProbeID = pProbeChromIDs;
		for(ProbeIdx = 0; ProbeIdx < NumProbeIDs; ProbeIdx++,pProbeID++)
			if(*pProbeID == TargHitID)
				break;
		if(ProbeIdx < NumProbeIDs)
			continue;

		if(TargHitLoci < (uint32_t)CurCoreSegOfs)					// must be able to extend core left
			continue;

		pEntry = &m_pEntriesBlock->Entries[TargHitID-1];
		if((TargHitLoci + (uint32_t)ProbeLen - (uint32_t)CurCoreSegOfs) > pEntry->SeqLen)					// must be able to extend core right
			continue;

		// check flanks and total number of missmatches, if <= MaxTotMM then return as having matched at least one other chrom
		pHitCoreSeq = &m_pSfxBlock->SeqSuffix[pEntry->StartOfs + TargHitLoci];
		pTarg = pHitCoreSeq - CurCoreSegOfs;
		pHitCoreSeq = pProbe;
		NumMM = 0;
		for(Idx = 0; NumMM > MaxTotMM && Idx < ProbeLen; Idx++, pTarg++, pHitCoreSeq++)
			{
			if((*pTarg & 0x0f) == (*pHitCoreSeq & 0x0f))
				continue;
			NumMM += 1;
			}
		if(NumMM <= MaxTotMM)
			return(1);
		}
	}

return(0);
}

//
// LocateAnyAlignment
// Utilises a sliding core
// Locates cores containing at most CoreMismatches
int						// < 0 if errors, 0 if no alignmentse, 1 if at least one alignment
CSfxArray::IfAnyAlignments(uint32_t ReadID,					// identifies this read
								int MaxTotMM,			        // max number of mismatches allowed
								int CoreLen,					// core window length
								int CoreDelta,					// core window offset increment (1..n)
								int MaxNumCoreSlides,			// max number of times to slide core on each strand
								eALStrand Align2Strand,			// align to this strand
								etSeqBase *pProbeSeq,			// probe sequence
								int ProbeLen,					// probe length
								int CurMaxIter,				// max allowed iterations per subsegmented sequence when matching that subsegment
								int NumAllocdIdentNodes,		// memory has been allocated by caller for holding up to this many tsIdentNodes
								tsIdentNode *pAllocsIdentNodes) // memory allocated by caller for holding tsIdentNodes
{
int CurCoreSegOfs;				// current core segment relative start
int IterCnt;					// count iterator for current segment target matches
int CurMMCnt;					// total number of mismatches for current target sequence being processed
int CoreMMCnt;					// total number of mismatches for current target core being processed
int LowMMCnt;					// least number of substitutions thus far required for match against target
int NxtLowMMCnt;				// next best least number of substitutions thus far required for match against target
int LowHitInstances;			// number of hit instances thus far for hits with LowMMCnt mismatches
etSeqBase *pProbeBase;
etSeqBase *pTargBase;
char CurStrand;
int64_t TargIdx;

tsIdentNode *pHashArray[cHashEntries + 1];		// hash array holding ptrs to identifier nodes
tsIdentNode *pIdentNodes = pAllocsIdentNodes;	// identifier nodes
tsIdentNode *pCurIdentNode;
tsIdentNode *pNewIdentNode;
int CurNumIdentNodes;
int Hash;

uint32_t TargSeqID;
uint32_t NumTargSeqProc;

bool bFirstIter;				// set false after the first subsequence core returned by LocateFirstExact has been processed

int64_t LastTargIdx;
uint32_t NumCopies;

uint32_t PatIdx;
int64_t TargSeqLeftIdx;
int64_t ProbeSeqLeftIdx;
int TargMatchLen;

etSeqBase BisBase;
uint8_t ProbeBase;
uint8_t TargBase;

int Cmp;
int CurNumCoreSlides;
int CurCoreDelta;
tsHitLoci *pCurHit;

etSeqBase *pTarg;			// target sequence
void *pSfxArray;			// target sequence suffix array
int64_t SfxLen;				// number of suffixs in pSfxArray
BisBase = eBaseN;

// ensure suffix array block loaded for iteration!
if (m_pSfxBlock == NULL || m_pSfxBlock->ConcatSeqLen == 0)
	return(eBSFerrInternal);

pTarg = (etSeqBase *)&m_pSfxBlock->SeqSuffix[0];
pSfxArray = (void *)&m_pSfxBlock->SeqSuffix[m_pSfxBlock->ConcatSeqLen];
SfxLen = m_pSfxBlock->ConcatSeqLen;


LowHitInstances = 0;
LowMMCnt =  MaxTotMM + 1;
NxtLowMMCnt = LowMMCnt;

pCurHit = NULL;


NumTargSeqProc = 0;
CurNumIdentNodes = 0;
CurNumCoreSlides = 0;

if (Align2Strand == eALSCrick)
	{
	if (m_bColorspace)
		CSeqTrans::ReverseSeq(ProbeLen, pProbeSeq);
	else
		CSeqTrans::ReverseComplement(ProbeLen, pProbeSeq);
	CurStrand = '-';
	}
else
	CurStrand = '+';

do
	{
	CurCoreDelta = CoreDelta;
	CurNumCoreSlides = 0;
	memset(pHashArray, 0, sizeof(pHashArray));
	CurNumIdentNodes = 0;
	for (CurCoreSegOfs = 0;
			CurNumCoreSlides < MaxNumCoreSlides &&
			CurCoreSegOfs <= (ProbeLen - CoreLen) &&
			CurCoreDelta > CoreLen / 3 &&
			CurNumIdentNodes < NumAllocdIdentNodes;	// can only allow upto cMaxNumIdentNodes to be saved
	CurNumCoreSlides += 1,
	CurCoreSegOfs += CurCoreDelta)
		{
		if (CurNumCoreSlides >= MaxNumCoreSlides)
			break;
		if ((CurCoreSegOfs + CoreLen + CurCoreDelta) > ProbeLen)
			CurCoreDelta = ProbeLen - (CurCoreSegOfs + CoreLen);

		TargIdx = LocateFirstExact(&pProbeSeq[CurCoreSegOfs], CoreLen, pTarg, m_pSfxBlock->SfxElSize, pSfxArray, 0, 0, SfxLen - 1);
		if (TargIdx == 0)        // 0 if no core segment matches
			continue;			// try for match on next core segment after shifting core to right

		TargIdx -= 1;
		IterCnt = 0;
		NumCopies = 0;
		bFirstIter = true;		// set false after the first subsequence core returned by LocateFirstExact has been processed
		while (!CurMaxIter || IterCnt < CurMaxIter)
			{
			if (CurNumIdentNodes >= NumAllocdIdentNodes)
				break;

			if (!bFirstIter)
				{

					// ensure not about to iterate past end of suffix array!
				if ((TargIdx + 1) >= (int64_t)m_pSfxBlock->ConcatSeqLen || (SfxOfsToLoci(m_pSfxBlock->SfxElSize, pSfxArray, TargIdx + 1) + CoreLen) >(int64_t)m_pSfxBlock->ConcatSeqLen)
					break;

				if (IterCnt == 100 && !NumCopies)
					{
						// check how many more exact copies there are of the current probe subsequence, if too many then don't bother exploring these
					LastTargIdx = LocateLastExact(&pProbeSeq[CurCoreSegOfs], CoreLen, pTarg, m_pSfxBlock->SfxElSize, pSfxArray, 0, TargIdx - 1, SfxLen - 1);
					NumCopies = LastTargIdx > 0 ? (uint32_t)(1 + LastTargIdx - TargIdx) : 0;
					if (CurMaxIter && NumCopies > (uint32_t)CurMaxIter)		// only checking at the 100th iteration allows a little slack
						break;										// try next core segment
					}

					// check that this new putative core is still matching
				pTargBase = &pTarg[SfxOfsToLoci(m_pSfxBlock->SfxElSize, pSfxArray, TargIdx + 1)];
				pProbeBase = &pProbeSeq[CurCoreSegOfs];
				if (m_bBisulfite)
					Cmp = BSCmpProbeTarg(pProbeBase, pTargBase, CoreLen);
				else
				{
					int Ofs;
					uint8_t El1;
					uint8_t El2;
					etSeqBase *pEl1 = pProbeBase;
					etSeqBase *pEl2 = pTargBase;
					Cmp = 0;
					for (Ofs = 0; Ofs < CoreLen; Ofs++)
					{
						El2 = *pEl2++ & 0x0f;
						if (El2 == eBaseEOS)
						{
							Cmp = -1;
							break;
						}
						El1 = *pEl1++ & 0x0f;
						if (El1 > El2)
						{
							Cmp = 1;
							break;
						}
						if (El1 < El2)
						{
							Cmp = -1;
							break;
						}
					}
				}

				if (Cmp != 0)				// will be non-zero if target no longer matches
					break;					// try next core segment
				TargIdx += 1;
				}

			bFirstIter = false;
			if (SfxOfsToLoci(m_pSfxBlock->SfxElSize, pSfxArray, TargIdx) < (uint32_t)CurCoreSegOfs)
				continue;

			TargMatchLen = ProbeLen;
			TargSeqLeftIdx = SfxOfsToLoci(m_pSfxBlock->SfxElSize, pSfxArray, TargIdx) - CurCoreSegOfs;
			ProbeSeqLeftIdx = 0;

			// ensure comparisons are still within start/end range of target sequence/assembly
			if (!TargMatchLen || ((uint64_t)TargSeqLeftIdx + (uint32_t)TargMatchLen) > m_pSfxBlock->ConcatSeqLen)
				continue;

				// check if target already processed
			TargSeqID = (uint32_t)(1 + SfxOfsToLoci(m_pSfxBlock->SfxElSize, pSfxArray, TargIdx) - (uint32_t)CurCoreSegOfs);
			Hash = (TargSeqID & cHashEntries);
			if ((pCurIdentNode = pHashArray[Hash]) == NULL)
				{
				pCurIdentNode = pHashArray[Hash] = &pIdentNodes[CurNumIdentNodes++];
				pCurIdentNode->TargSeqID = TargSeqID;
				pCurIdentNode->pNxt = NULL;
				}
			else
				{
				while (pCurIdentNode != NULL && pCurIdentNode->TargSeqID != TargSeqID)
					pCurIdentNode = pCurIdentNode->pNxt;
				if (pCurIdentNode != NULL)							// will be non-null if have already processed this sequence
					continue;
				pNewIdentNode = &pIdentNodes[CurNumIdentNodes++];
				pNewIdentNode->TargSeqID = TargSeqID;
				pNewIdentNode->pNxt = pHashArray[Hash];
				pHashArray[Hash] = pNewIdentNode;
				}
			NumTargSeqProc += 1;
			IterCnt += 1;

				// now do the matching allowing for missmatches
				// if colorspace then need to allow for adjacent apparent mismatches actually being the result of a single substitution relative to the target
				// if colorspace and a standalone, without an adjacent, mismatch then this is supposedly a sequencing error, still treat as though a mismatch
				{
					pTargBase = &pTarg[TargSeqLeftIdx];
					pProbeBase = pProbeSeq;
					CurMMCnt = 0;
					CoreMMCnt = 0;
					bool bPairMM = false;
					if (m_bBisulfite)
						BisBase = GetBisBase(TargMatchLen, pTargBase, pProbeBase);
					for (PatIdx = 0; PatIdx < (uint32_t)TargMatchLen; PatIdx++, pTargBase++, pProbeBase++)
					{
						TargBase = *pTargBase & 0x0f;
						ProbeBase = *pProbeBase & 0x0f;
						if (TargBase == eBaseEOS)		// mustn't match across entry sequences
							break;

						if (m_bColorspace)		// in colorspace, unpaired mismatches are sequencer errors but still treat as if a substitution
						{
							if (!bPairMM && TargBase != ProbeBase)
							{
								if (PatIdx < ((uint32_t)TargMatchLen - 1))
								{
									if ((pTargBase[1] & 0x0f) == (pProbeBase[1] & 0x0f))
									{
										if (++CurMMCnt > MaxTotMM)
											break;
										if (CurMMCnt >= NxtLowMMCnt)
											break;
										continue;
									}
								}
								// accept as being a mismatch
								bPairMM = true;
							}
							else
							{
								bPairMM = false;
								continue;
							}
						}
						else				// in basespace
						{
							if (ProbeBase == TargBase)
								continue;

							if (m_bBisulfite)
							{
								if (TargBase == eBaseC || TargBase == eBaseG)
								{
									switch (BisBase)
									{
										case eBaseA:		// only allow A
											if (ProbeBase == eBaseA && TargBase == eBaseG)
												continue;
											break; // mismatch
										case eBaseT:		// only allow T
											if (ProbeBase == eBaseT && TargBase == eBaseC)
												continue;
											break; // mismatch
									}
								}
							}
						}

						// execution here only if mismatch
						if (++CurMMCnt > MaxTotMM)
							break;
						if (CurMMCnt >= NxtLowMMCnt)
							break;
					}
				if (PatIdx != TargMatchLen)
					continue;

				// processing continues here only if number of mismatches over whole probe accepted
				return(eHRhits);
				}
			}
		}
	if (CurStrand == '+' && Align2Strand == eALSboth)	// if just processed watson '+' strand then will need to process crick or '-' strand if processing both strands
		{
		if (m_bColorspace)
			CSeqTrans::ReverseSeq(ProbeLen, pProbeSeq);
		else
			CSeqTrans::ReverseComplement(ProbeLen, pProbeSeq);
		CurStrand = '-';
		Align2Strand = eALSCrick;
		NumTargSeqProc = 0;
		}
	else									// either processed both or just the crick strand so no further processing required
		Align2Strand = eALSnone;			// two passes max - first on the '+' strand with optional 2nd on '-' strand
	}
while(Align2Strand != eALSnone);

if (CurStrand == '-')									// restore probe sequence if had started processing '-' strand
	{
	if (m_bColorspace)
		CSeqTrans::ReverseSeq(ProbeLen, pProbeSeq);
	else
		CSeqTrans::ReverseComplement(ProbeLen, pProbeSeq);
	}
return(eHRnone);				// no putative hits
}


int											// if < eBSFSuccess then errors, 0 if unable to adaptively trim, >= MinTrimLen if able to adaptively trim
CSfxArray::AdaptiveTrim(uint32_t SeqLen,	// untrimmed probe and target sequence are both this bp long
			etSeqBase *pProbeSeq,			// end trimming this probe sequence
			etSeqBase *pTargSeq,			// when aligned against this target sequence
			uint32_t ChromID,				// target sequence is part of this chromosome sequence
			uint32_t ChromSeqLen,			// chrom sequence is this length
			uint32_t MinTrimLen,			// accepted trimmed probe sequence must be at least this minimum length
			uint32_t MaxMM,					// and contain at most this many mismatches proportional to a 100bp trimmed sequence match
			uint32_t MinFlankMatches,		// 5' and 3' flanks of trimmed probe must start/end with at least this many exactly matching bases
			uint32_t *pTrimSeqLen,			// trimmed probe is this length
			uint32_t *pTrimStart,			// accepted trimmed probe has this many bp trimmed from start (5') of probe
			uint32_t *pTrimEnd,				// accepted trimmed probe has this many bp trimmed from end (3') of probe
			uint32_t *pTrimMMs)				// after trimming there are this many mismatches in the trimmed probe
{
uint32_t Ofs;
uint32_t NumRegions;
uint32_t NumRegionsMinExact;
uint32_t CurRegionIdx;
uint32_t CurStartRegionIdx;
uint32_t FirstStartRegionIdx;
uint32_t LastStartRegionIdx;
uint32_t FirstEndRegionIdx;
uint32_t LastEndRegionIdx;
tsATRegion ATRegions[cMaxATSeqLen];
tsATRegion *pCurATRegion;
tsATRegion *pCurStartATRegion;

etSeqBase *pProbeBase;
etSeqBase *pTargBase;
bool bCurMM;

if(pTrimSeqLen != NULL)
	*pTrimSeqLen = 0;
if(pTrimStart != NULL)
	*pTrimStart = 0;
if(pTrimEnd != NULL)
	*pTrimEnd	= 0;
if(pTrimMMs != NULL)
	*pTrimMMs = 0;

// some simple parameter validation
if(SeqLen < cMinATSeqLen || SeqLen > cMaxATSeqLen || 
	pProbeSeq == NULL || pTargSeq == NULL ||
	MinTrimLen < cMinATTrimmedLen || MinTrimLen > SeqLen ||
	MaxMM > (((cMaxATMM * SeqLen) + 99) / 100) || MinFlankMatches > cMaxATMaxFlankMatches)
	return(eBSFerrParams);
if(MinFlankMatches == 0)					// ensuring adaptively trimmed sequence always start/end on an exact match
	MinFlankMatches = 1;

pProbeBase = pProbeSeq;
pTargBase = pTargSeq;

if(MinTrimLen == SeqLen) // if MinTrimLen same as SeqLen then simply need to total up number of mismatches
	{
	uint32_t MaxAllowedMMs = ((SeqLen * MaxMM)+99)/100;
	uint32_t SeqMMs = 0;
	for(Ofs=0; Ofs < SeqLen; Ofs+=1, pProbeBase+=1, pTargBase+=1)
		{
		bCurMM = (*pProbeBase & 0x0f) == (*pTargBase & 0x0f) ? false : true;
		if(bCurMM)
			{
			if(++SeqMMs > MaxAllowedMMs)
				break;
			if(Ofs < MinFlankMatches || (SeqLen - Ofs) < MinFlankMatches)
				{
				SeqMMs = MaxAllowedMMs+1;
				break;
				}
			}
		}
	if(SeqMMs <= MaxAllowedMMs)
		{
		if(pTrimSeqLen != NULL)
			*pTrimSeqLen = SeqLen;
		if(pTrimMMs != NULL)
			*pTrimMMs = SeqMMs;
		return(SeqLen);
		}
	return(0);
	}

// determine number and sizes of regions containing either mismatches or exact matches
pProbeBase = pProbeSeq;
pTargBase = pTargSeq;
NumRegions = 0;
pCurATRegion = NULL;
bCurMM = false;
NumRegionsMinExact = 0;
uint32_t NumExactInRegions = 0;
for(Ofs=0; Ofs < SeqLen; Ofs+=1, pProbeBase+=1, pTargBase+=1)
	{
	bCurMM = (*pProbeBase & 0x0f) == (*pTargBase & 0x0f) ? false : true;
	if(!bCurMM)
		NumExactInRegions += 1;
	if(pCurATRegion == NULL || bCurMM != pCurATRegion->flgMM ? true : false)
		{
		pCurATRegion = &ATRegions[NumRegions++];
		pCurATRegion->flgMM = bCurMM ? 1 : 0;
		pCurATRegion->RegionLen = 1;
		pCurATRegion->RegionOfs = Ofs;
		}
	else
		{
		pCurATRegion->RegionLen += 1;
		if(pCurATRegion->RegionLen == cMinATExactLen && !bCurMM)
			NumRegionsMinExact += 1;
		}
	}

if(!NumRegionsMinExact)
	return(0);

// for each region determine if that region would be suitable for starting and/or ending a trimmed sequence
pCurATRegion = ATRegions;
CurRegionIdx = 0;
FirstStartRegionIdx = 0;
LastStartRegionIdx = 0;
FirstEndRegionIdx = 0;
LastEndRegionIdx = 0;
for(CurRegionIdx = 0; CurRegionIdx < NumRegions; CurRegionIdx++,pCurATRegion++)
	{
	if(pCurATRegion->flgMM == 0 && pCurATRegion->RegionLen >= MinFlankMatches) // putative trimmed sequences must start/end with at least MinFlankMatches exact matches
		{
		if(pCurATRegion->RegionOfs <= SeqLen - MinTrimLen)
			{
			pCurATRegion->flgTrim5 = 1;
			LastStartRegionIdx = CurRegionIdx+1;
			if(FirstStartRegionIdx == 0)
				FirstStartRegionIdx = LastStartRegionIdx;
			}
		else
			pCurATRegion->flgTrim5 = 0;
		if(pCurATRegion->RegionOfs + pCurATRegion->RegionLen >= (uint16_t)MinTrimLen)
			{
			pCurATRegion->flgTrim3 = 1;
			LastEndRegionIdx = CurRegionIdx+1;
			if(FirstEndRegionIdx == 0)
				FirstEndRegionIdx = LastEndRegionIdx;
			}
		else
			pCurATRegion->flgTrim3 = 0;
		}
	else
		{
		pCurATRegion->flgTrim5 = 0;
		pCurATRegion->flgTrim3 = 0;
		}
	}

if(FirstStartRegionIdx == 0 || FirstEndRegionIdx == 0)
	return(0);

// iterate over each putative starting region and maximally extend that region into ending regions providing that number of mismatches would be accepted
uint32_t CurTrimmedSeqLen;
uint32_t CurTrimmedSeqMM;
uint32_t CurTrimmedRegionIdx;

uint32_t BestTrimmedSeqLen;
uint32_t BestTrimmedSeqMM;
uint32_t BestTrimStart;
uint32_t BestTrimEnd;

BestTrimmedSeqLen = 0;
BestTrimmedSeqMM = 0;
BestTrimStart = 0;
BestTrimEnd = 0;
pCurStartATRegion = &ATRegions[FirstStartRegionIdx-1];
for(CurStartRegionIdx = FirstStartRegionIdx-1; CurStartRegionIdx < LastStartRegionIdx; CurStartRegionIdx++,pCurStartATRegion++)
	{
	if(pCurStartATRegion->flgTrim5 == 0)   // find region to start trimming from
		continue;
	CurTrimmedSeqLen = 0;
	CurTrimmedSeqMM = 0;
	pCurATRegion = pCurStartATRegion;
	CurTrimmedRegionIdx = CurStartRegionIdx;
	while(CurTrimmedRegionIdx++ < LastEndRegionIdx)
		{
		CurTrimmedSeqLen +=  pCurATRegion->RegionLen;
		if(pCurATRegion->flgMM)
			{
			if(MaxMM == 0)
				break;
			CurTrimmedSeqMM += pCurATRegion->RegionLen;
			if((MaxMM+1.0) / 100.0 <= (double)CurTrimmedSeqMM/(SeqLen - pCurStartATRegion->RegionOfs))
				break;
			}
		else
			{
			if(BestTrimmedSeqLen == 0)
				{
				BestTrimStart = pCurStartATRegion->RegionOfs;
				BestTrimEnd = SeqLen - (BestTrimStart + CurTrimmedSeqLen);
				BestTrimmedSeqLen = CurTrimmedSeqLen;
				BestTrimmedSeqMM = 0;
				pCurATRegion += 1;
				continue;
				}
			}

		if(CurTrimmedSeqLen < MinTrimLen || !pCurATRegion->flgTrim3)
			{
			pCurATRegion += 1;
			continue;
			}

		pCurATRegion += 1;

		// check if the currently trimmed sequence length would be acceptable for mismatches
		if((MaxMM+1.0) / 100.0 <= (double)CurTrimmedSeqMM/CurTrimmedSeqLen)
			continue;

		// prepared to accept the trimmed sequence but only accept if longer than any previously accepted trimmed sequence
		if(BestTrimmedSeqLen < CurTrimmedSeqLen || (BestTrimmedSeqLen == CurTrimmedSeqLen && (BestTrimmedSeqMM == 0 || CurTrimmedSeqMM < BestTrimmedSeqMM)))
			{
			BestTrimStart = pCurStartATRegion->RegionOfs;
			BestTrimEnd = SeqLen - (BestTrimStart + CurTrimmedSeqLen);
			BestTrimmedSeqLen = CurTrimmedSeqLen;
			BestTrimmedSeqMM  = CurTrimmedSeqMM;
			}
		}
	}

if(BestTrimmedSeqLen >= MinTrimLen)
	{
	if(pTrimSeqLen != NULL)
		*pTrimSeqLen = BestTrimmedSeqLen;
	if(pTrimStart != NULL)
		*pTrimStart = BestTrimStart;
	if(pTrimEnd != NULL)
		*pTrimEnd	= BestTrimEnd;
	if(pTrimMMs != NULL)
		*pTrimMMs = BestTrimmedSeqMM;
	return(BestTrimmedSeqLen);
	}
return(0);
}

//
// LocateCoreMultiples
// Utilises a sliding core
// Locates cores containing at most CoreMismatches
//
// Locates probe to target matches with mismatches
// If m_bBisulfite is set true then matches using bisulfite rules
//
int						// < 0 if errors, 0 if no matches or change, 1 if mumber of matches accepted, 2 MMDelta criteria not met, 3 too many match instances
CSfxArray::LocateCoreMultiples(uint32_t ExtdProcFlags,	// flags indicating if lower levels need to do any form of extended processing with this specific read...
						 uint32_t ReadID,					// identifies this read
						 int MinChimericLen,			// minimum chimeric length as a percentage (0 to disable, otherwise 50..99) of probe sequence
						 int MaxTotMM,			        // max number of mismatches allowed
						 int CoreLen,					// core window length
						 int CoreDelta,					// core window offset increment (1..n)
						 int MaxNumCoreSlides,			// max number of times to slide core on each strand
						 int MMDelta,					// minimum (1..n) mismatch difference between the best and next best core alignment
 					 	 eALStrand Align2Strand,		// align to this strand
  						 int *pLowHitInstances,			// In/Out number of match instances for lowest number of mismatches thus far for this read
						 int *pLowMMCnt,				// In/Out lowest number of mismatches thus far for this read
						 int *pNxtLowMMCnt,				// In/Out next to lowest number of mismatches thus far for this read
						 etSeqBase *pProbeSeq,			// probe
						 int ProbeLen,					// probe length
						 int MaxHits,					// (IN) process for at most this number of hits
						 tsHitLoci *pHits,				// where to return hits (at most MaxHits)
						 int CurMaxIter,				// max allowed iterations per subsegmented sequence when matching that subsegment
						 int NumAllocdIdentNodes,		// memory has been allocated by caller for holding up to this many tsIdentNodes
						 tsIdentNode *pAllocsIdentNodes) // memory allocated by caller for holding tsIdentNodes
{
int CurCoreSegOfs;				// current core segment relative start
int IterCnt;					// count iterator for current segment target matches
int CurMMCnt;					// total number of mismatches for current target sequence being processed
int CoreMMCnt;					// total number of mismatches for current target core being processed
int LowMMCnt;					// least number of substitutions thus far required for match against target
int NxtLowMMCnt;				// next best least number of substitutions thus far required for match against target
int LowHitInstances;			// number of hit instances thus far for hits with LowMMCnt mismatches
etSeqBase *pProbeBase;
etSeqBase *pTargBase;
char CurStrand;
int64_t TargIdx;

uint32_t MinProbeChimericLen;
char BestChimericStrand;
int BestChimericLen;
uint16_t BestTrim5Flank;
uint16_t BestTrim3Flank;
int BestMaxChimericMMs;

tsIdentNode *pHashArray[cHashEntries+1];		// hash array holding ptrs to identifier nodes
tsIdentNode *pIdentNodes = pAllocsIdentNodes;	// identifier nodes
tsIdentNode *pCurIdentNode;
tsIdentNode *pNewIdentNode;
int CurNumIdentNodes;
int Hash;

uint32_t TargSeqID;
uint32_t NumTargSeqProc;

bool bFirstIter;				// set false after the first subsequence core returned by LocateFirstExact has been processed

uint32_t PatIdx;
int64_t TargSeqLeftIdx;
int64_t ProbeSeqLeftIdx;
int TargMatchLen;

etSeqBase BisBase;
uint8_t ProbeBase;
uint8_t TargBase;

int Cmp;
int CurNumCoreSlides;
int CurCoreDelta;
tsHitLoci *pCurHit;

etSeqBase *pTarg;			// target sequence
void *pSfxArray;			// target sequence suffix array
int64_t SfxLen;				// number of suffixs in pSfxArray
tsSfxEntry *pEntry;
BisBase=eBaseN;

// force MinChimericLen to be either 0, or in the range 15..99. if MinChimericLen is < 15 or > 99 then treat as if a normal full length match required 
if(MinChimericLen >= 15 && MinChimericLen <= 99)
	{
	MinProbeChimericLen = max(CoreLen,(MinChimericLen * ProbeLen) / 100);	// clamp to be at least the core length!
	}
else
	MinProbeChimericLen = 0;

// ensure suffix array block loaded for iteration!
if(m_pSfxBlock == NULL || m_pSfxBlock->ConcatSeqLen == 0)
	return(eBSFerrInternal);

// if already matching exactly and > MaxHits then can't improve, treat as same or new LowMMCnt but with too many multiple instances
if(*pLowHitInstances > MaxHits && *pLowMMCnt == 0)
   return(eHRHitInsts);

// if already matching exactly and less than required mismatch deltas then can't improve
if(*pLowHitInstances >= 1 && *pLowMMCnt == 0 && (*pNxtLowMMCnt - *pLowMMCnt) < MMDelta)
	return(eHRMMDelta);

pTarg = (etSeqBase *)&m_pSfxBlock->SeqSuffix[0];
pSfxArray = (void *)&m_pSfxBlock->SeqSuffix[m_pSfxBlock->ConcatSeqLen];
SfxLen = m_pSfxBlock->ConcatSeqLen;


if(*pLowHitInstances <= 0 || *pLowMMCnt < 0 || *pNxtLowMMCnt < 0)	// if never seen any previous matches then ensure substitution counts are initialised
	{
	LowHitInstances = *pLowHitInstances = 0;
	LowMMCnt = *pLowMMCnt = MaxTotMM + MMDelta + 1;
	NxtLowMMCnt = *pNxtLowMMCnt = LowMMCnt;
	}
else
	{
	LowHitInstances = *pLowHitInstances;
	LowMMCnt = *pLowMMCnt;
	NxtLowMMCnt = *pNxtLowMMCnt;
	}

if(LowHitInstances < MaxHits)
	pCurHit = &pHits[LowHitInstances];
else
	pCurHit = NULL;


NumTargSeqProc = 0;
CurNumIdentNodes = 0;
CurNumCoreSlides = 0;

if(Align2Strand == eALSCrick)
	{
	if(m_bColorspace)
		CSeqTrans::ReverseSeq(ProbeLen,pProbeSeq);
	else
		CSeqTrans::ReverseComplement(ProbeLen,pProbeSeq);
	CurStrand = '-';
	}
else
	CurStrand = '+';

BestChimericStrand = '?';
BestChimericLen = 0;
BestTrim5Flank = 0;
BestTrim3Flank = 0;
BestMaxChimericMMs = 0;

do
	{
	CurCoreDelta = CoreDelta;
	CurNumCoreSlides = 0;
	memset(pHashArray,0,sizeof(pHashArray));
	CurNumIdentNodes = 0;
	for(CurCoreSegOfs = 0;
		CurNumCoreSlides < MaxNumCoreSlides &&
		CurCoreSegOfs <= (ProbeLen - CoreLen) &&
		CurCoreDelta > CoreLen/3 &&
		CurNumIdentNodes < NumAllocdIdentNodes;	// can only allow upto cMaxNumIdentNodes to be saved
	    CurNumCoreSlides += 1,
	    CurCoreSegOfs += CurCoreDelta)
		{
		if(CurNumCoreSlides >= MaxNumCoreSlides)
			break;
		if((CurCoreSegOfs + CoreLen + CurCoreDelta) > ProbeLen)
			CurCoreDelta = ProbeLen - (CurCoreSegOfs + CoreLen);

		if (m_pOccKMerClas != NULL && CoreLen == m_OccKMerLen && (OverOccKMerClas(CoreLen, &pProbeSeq[CurCoreSegOfs]) != 1))	// 1 if number of K-mers of CoreLen is within range 1..max
			continue;

		TargIdx = LocateFirstExact(&pProbeSeq[CurCoreSegOfs],CoreLen,pTarg,m_pSfxBlock->SfxElSize,pSfxArray,0,0,SfxLen-1);
		if(TargIdx == 0)        // 0 if no core segment matches
			continue;			// try for match on next core segment after shifting core to right

		TargIdx -= 1;
		IterCnt = 0;
		bFirstIter = true;		// set false after the first subsequence core returned by LocateFirstExact has been processed
		while(!CurMaxIter || IterCnt < CurMaxIter)
			{
			if(CurNumIdentNodes >= NumAllocdIdentNodes)
				break;

			if(!bFirstIter) {

				// ensure not about to iterate past end of suffix array!
				if((TargIdx + 1) >= (int64_t)m_pSfxBlock->ConcatSeqLen || (SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,TargIdx+1) + CoreLen) >  (int64_t)m_pSfxBlock->ConcatSeqLen)
					break;

				// check that this new putative core is still matching
				pTargBase = &pTarg[SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,TargIdx+1)];
				pProbeBase = &pProbeSeq[CurCoreSegOfs];
				if(m_bBisulfite)
					Cmp = BSCmpProbeTarg(pProbeBase,pTargBase,CoreLen);
				else
					{
					int Ofs;
					uint8_t El1;
					uint8_t El2;
					etSeqBase *pEl1= pProbeBase;
					etSeqBase *pEl2 = pTargBase;
					Cmp = 0;
					for(Ofs=0; Ofs < CoreLen; Ofs++)
						{
						El2 = *pEl2++ & 0x0f;
						if(El2 == eBaseEOS)
							{
							Cmp = -1;
							break;
							}
						El1 = *pEl1++ & 0x0f;
						if(El1 > El2)
							{
							Cmp = 1;
							break;
							}
						if(El1 < El2)
							{
							Cmp = -1;
							break;
							}
						}
					}

				if(Cmp != 0)				// will be non-zero if target no longer matches
					break;					// try next core segment
				TargIdx += 1;
				}

			bFirstIter = false;
			if(SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,TargIdx) < (uint32_t)CurCoreSegOfs)
				continue;

			TargMatchLen = ProbeLen;
			TargSeqLeftIdx = SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,TargIdx) - CurCoreSegOfs;
			ProbeSeqLeftIdx = 0;

			pEntry = MapChunkHit2Entry(TargSeqLeftIdx);

			// ensure comparisons are still within start/end range of target sequence/assembly
			if(pEntry == NULL || !TargMatchLen || ((uint64_t)TargSeqLeftIdx + (uint32_t)TargMatchLen - 1) > pEntry->EndOfs)
				continue;

			// check if target already processed
			TargSeqID = (uint32_t)(1 + SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,TargIdx) - (uint32_t)CurCoreSegOfs);
			Hash = (TargSeqID & cHashEntries);
			if((pCurIdentNode = pHashArray[Hash])==NULL)
				{
				pCurIdentNode = pHashArray[Hash] = &pIdentNodes[CurNumIdentNodes++];
				pCurIdentNode->TargSeqID = TargSeqID;
				pCurIdentNode->pNxt = NULL;
				}
			else
				{
				while(pCurIdentNode != NULL && pCurIdentNode->TargSeqID != TargSeqID)
					pCurIdentNode = pCurIdentNode->pNxt;
				if(pCurIdentNode != NULL)							// will be non-null if have already processed this sequence
					continue;
				pNewIdentNode = &pIdentNodes[CurNumIdentNodes++];
				pNewIdentNode->TargSeqID = TargSeqID;
				pNewIdentNode->pNxt = pHashArray[Hash];
				pHashArray[Hash] = pNewIdentNode;
				}
			NumTargSeqProc += 1;
			IterCnt += 1;

			// now do the matching allowing for mismatches
			// if colorspace then need to allow for adjacent apparent mismatches actually being the result of a single substitution relative to the target
			// if colorspace and a standalone, without an adjacent, mismatch then this is supposedly a sequencing error, still treat as though a mismatch


			if(MinProbeChimericLen > 0) // this chimeric processing is a naive implementation and would be a prime target for optimisation if time becomes available
				{
				pTargBase = &pTarg[TargSeqLeftIdx];
				pProbeBase = pProbeSeq;
				CurMMCnt = 0;
				CoreMMCnt = 0;
				bool bPairMM = false;
				bool bInMatches;

				etSeqBase *pProbeSubSeq;
				etSeqBase *pTargSubSeq;

				int CurChimericLen;
				int MaxChimericMMs;
				int MaxChimericLen;
				uint16_t Trim5Flank;
				uint16_t Trim3Flank;


				pTargSubSeq = pTargBase;
				pProbeSubSeq = pProbeBase;
				MaxChimericLen = 0;
				MaxChimericMMs = 0;
				CurChimericLen = 0;
				bInMatches = false;

				int AdaptiveTrimRslt;
				uint32_t TrimmedSeqLen;
				uint32_t TrimmedLeftFlankLen;
				uint32_t TrimmedRightFlankLen;
				uint32_t TrimmedMM;
				uint32_t TargSeqLen;
				TargSeqLen = GetSeqLen(TargSeqID);
				AdaptiveTrimRslt = AdaptiveTrim(TargMatchLen,pProbeSeq,&pTarg[TargSeqLeftIdx], TargSeqID,TargSeqLen,MinProbeChimericLen,MaxTotMM,3,&TrimmedSeqLen,&TrimmedLeftFlankLen,&TrimmedRightFlankLen,&TrimmedMM);
				MaxChimericLen = TrimmedSeqLen; 				
				MaxChimericMMs = TrimmedMM;
				Trim5Flank = TrimmedLeftFlankLen;
				Trim3Flank = TrimmedRightFlankLen;

				if(MaxChimericLen < (int)MinProbeChimericLen)
					continue;

				if(MaxChimericLen > BestChimericLen ||	// if at least as long as any previous match then this is a new unique putative hit
					(MaxChimericLen == BestChimericLen && MaxChimericMMs < BestMaxChimericMMs))
					{
					if(BestChimericLen > 0 && MaxChimericLen > BestChimericLen)
						LowMMCnt = MaxChimericMMs + MMDelta + 1;

					BestChimericLen = MaxChimericLen;
					BestTrim5Flank = Trim5Flank;
					BestTrim3Flank = Trim3Flank;
					BestMaxChimericMMs = MaxChimericMMs;
					BestChimericStrand = CurStrand;
					pCurHit = pHits;
					LowHitInstances = 1;

					NxtLowMMCnt = LowMMCnt;
					LowMMCnt = MaxChimericMMs;

					pEntry = MapChunkHit2Entry(TargSeqLeftIdx + Trim5Flank);
					pCurHit->FlgChimeric = 1;
					pCurHit->FlgInDel = 0;
					pCurHit->FlgInsert = 0;
					pCurHit->FlgSplice = 0;
					pCurHit->FlgNonOrphan = 0;
					memset(pCurHit->Seg,0,sizeof(pCurHit->Seg));
					pCurHit->Seg[0].Strand = CurStrand;
					if(CurStrand == '+')
						{
						pCurHit->Seg[0].TrimLeft = Trim5Flank;
						pCurHit->Seg[0].TrimRight = Trim3Flank;
						}
					else
						{
						pCurHit->Seg[0].TrimLeft = Trim3Flank;
						pCurHit->Seg[0].TrimRight = Trim5Flank;
						}
					pCurHit->Seg[0].ChromID = pEntry->EntryID;
					pCurHit->Seg[0].MatchLoci = (uint32_t)(TargSeqLeftIdx - pEntry->StartOfs);
					pCurHit->Seg[0].MatchLen = ProbeLen;
					pCurHit->Seg[0].Mismatches = MaxChimericMMs;
					pCurHit->Seg[0].TrimMismatches = MaxChimericMMs;
					pCurHit->BisBase = BisBase;
					}
				else							//	MaxChimericMMs >= LowMMCnt
					if(MaxChimericLen == BestChimericLen && MaxChimericMMs == BestMaxChimericMMs)	// multiple instances with this number of mismatches?
						{
						LowHitInstances += 1;
						if(pCurHit != NULL && LowHitInstances <= MaxHits)
							{
							pCurHit += 1;
							pEntry = MapChunkHit2Entry(TargSeqLeftIdx + Trim5Flank);
							pCurHit->FlgChimeric = 1;
							pCurHit->FlgInDel = 0;
							pCurHit->FlgInsert = 0;
							pCurHit->FlgSplice = 0;
							pCurHit->FlgNonOrphan = 0;
							memset(pCurHit->Seg,0,sizeof(pCurHit->Seg));
							pCurHit->Seg[0].Strand = CurStrand;
							if(CurStrand == '+')
								{
								pCurHit->Seg[0].TrimLeft = Trim5Flank;
								pCurHit->Seg[0].TrimRight = Trim3Flank;
								}
							else
								{
								pCurHit->Seg[0].TrimLeft = Trim3Flank;
								pCurHit->Seg[0].TrimRight = Trim5Flank;
								}
							pCurHit->Seg[0].ChromID = pEntry->EntryID;
							pCurHit->Seg[0].MatchLoci = (uint32_t)(TargSeqLeftIdx - pEntry->StartOfs);
							pCurHit->Seg[0].MatchLen = ProbeLen;
							pCurHit->Seg[0].Mismatches = MaxChimericMMs;
							pCurHit->Seg[0].TrimMismatches = MaxChimericMMs;
							pCurHit->BisBase = BisBase;
							}
						}
					else						//	CurMMCnt > LowMMCnt
						{
						if(MaxChimericLen == BestChimericLen && MaxChimericMMs < NxtLowMMCnt) // is this an instance with the next fewest mismatches?
							NxtLowMMCnt = MaxChimericMMs;
						}

				if(MaxChimericLen == ProbeLen && LowHitInstances > MaxHits && LowMMCnt == 0)
					break;
				}
			else   // standard, non-chimeric, full read match required processing
				{
				pTargBase = &pTarg[TargSeqLeftIdx];
				pProbeBase = pProbeSeq;
				CurMMCnt = 0;
				CoreMMCnt = 0;
				bool bPairMM = false;
				if(m_bBisulfite)
					BisBase = GetBisBase(TargMatchLen,pTargBase,pProbeBase);

				for(PatIdx = 0; PatIdx < (uint32_t)TargMatchLen; PatIdx++,pTargBase++,pProbeBase++)
					{
					TargBase = *pTargBase & 0x0f;
					ProbeBase = *pProbeBase & 0x0f;
					if(TargBase == eBaseEOS)		// mustn't match across entry sequences
						break;

					if(m_bColorspace)		// in colorspace, unpaired mismatches are sequencer errors but still treat as if a substitution
						{
						if(!bPairMM && TargBase != ProbeBase)
							{
							if(PatIdx < ((uint32_t)TargMatchLen-1))
								{
								if((pTargBase[1] & 0x0f) == (pProbeBase[1] & 0x0f))
									{
									if(++CurMMCnt > MaxTotMM)
										break;
									if(CurMMCnt >= NxtLowMMCnt)
										break;
									continue;
									}
								}
							// accept as being a mismatch
							bPairMM = true;
							}
						else
							{
							bPairMM = false;
							continue;
							}
						}
					else				// in basespace
						{
						if(ProbeBase == TargBase)
							continue;

						if(m_bBisulfite)
							{
							if(TargBase == eBaseC || TargBase == eBaseG)
								{
								switch(BisBase) {
									case eBaseA:		// only allow A
										if(ProbeBase == eBaseA && TargBase == eBaseG)
											continue;
										break; // mismatch
									case eBaseT:		// only allow T
										if(ProbeBase == eBaseT && TargBase == eBaseC)
											continue;
										break; // mismatch
									}
								}
							}
						}

					// execution here only if mismatch
					if(++CurMMCnt > MaxTotMM)
						break;
					if(CurMMCnt >= NxtLowMMCnt)
						break;
					}
				if(PatIdx != TargMatchLen)
					continue;

				// processing continues here only if number of mismatches over whole probe accepted
				if(CurMMCnt < LowMMCnt)	// if fewer mismatches than any previous match then this is a new unique putative hit
					{
					pCurHit = pHits;
					LowHitInstances = 1;
					NxtLowMMCnt = LowMMCnt;
					LowMMCnt = CurMMCnt;
					pEntry = MapChunkHit2Entry(TargSeqLeftIdx);
					pCurHit->FlgChimeric = 0;
					pCurHit->FlgInDel = 0;
					pCurHit->FlgInsert = 0;
					pCurHit->FlgSplice = 0;
					pCurHit->FlgNonOrphan = 0;
					memset(pCurHit->Seg,0,sizeof(pCurHit->Seg));
					pCurHit->Seg[0].Strand = CurStrand;
					pCurHit->Seg[0].ChromID = pEntry->EntryID;
					pCurHit->Seg[0].MatchLoci = (uint32_t)(TargSeqLeftIdx - pEntry->StartOfs);
					pCurHit->Seg[0].MatchLen = ProbeLen;
					pCurHit->Seg[0].Mismatches = CurMMCnt;
					pCurHit->Seg[0].TrimMismatches = CurMMCnt;
					pCurHit->BisBase = BisBase;
					}
				else							//	CurMMCnt >= LowMMCnt
					if(CurMMCnt == LowMMCnt)	// multiple instances with this number of mismatches?
						{
						LowHitInstances += 1;
						if(pCurHit != NULL && LowHitInstances <= MaxHits)
							{
							pCurHit += 1;
							pEntry = MapChunkHit2Entry(TargSeqLeftIdx);
							pCurHit->FlgChimeric = 0;
							pCurHit->FlgInDel = 0;
							pCurHit->FlgInsert = 0;
							pCurHit->FlgSplice = 0;
							pCurHit->FlgNonOrphan = 0;
							memset(pCurHit->Seg,0,sizeof(pCurHit->Seg));
							pCurHit->Seg[0].Strand = CurStrand;
							pCurHit->Seg[0].ChromID = pEntry->EntryID;
							pCurHit->Seg[0].MatchLoci = (uint32_t)(TargSeqLeftIdx - pEntry->StartOfs);
							pCurHit->Seg[0].MatchLen = ProbeLen;
							pCurHit->Seg[0].Mismatches = CurMMCnt;
							pCurHit->Seg[0].TrimMismatches = CurMMCnt;
							pCurHit->BisBase = BisBase;
							}
						}
					else						//	CurMMCnt > LowMMCnt
						{
						if(CurMMCnt < NxtLowMMCnt) // is this an instance with the next fewest mismatches?
							NxtLowMMCnt = CurMMCnt;
						}
				if(LowHitInstances > MaxHits && LowMMCnt == 0)
					break;
				}
			}
		if(LowHitInstances > MaxHits && LowMMCnt == 0)
			{
			Align2Strand = eALSnone;					// can't improve so early exit
			break;
			}
		}
	if(CurStrand == '+' && Align2Strand == eALSboth)	// if just processed watson '+' strand then will need to process crick or '-' strand if processing both strands
		{
		if(m_bColorspace)
			CSeqTrans::ReverseSeq(ProbeLen,pProbeSeq);
		else
			CSeqTrans::ReverseComplement(ProbeLen,pProbeSeq);
		CurStrand = '-';
		Align2Strand = eALSCrick;
		NumTargSeqProc = 0;
		}
	else									// either processed both or just the crick strand so no further processing required
		Align2Strand = eALSnone;			// two passes max - first on the '+' strand with optional 2nd on '-' strand
	}
while(!(LowHitInstances > MaxHits && LowMMCnt == 0) && Align2Strand != eALSnone);

if(CurStrand == '-')									// restore probe sequence if had started processing '-' strand
		if(m_bColorspace)
			CSeqTrans::ReverseSeq(ProbeLen,pProbeSeq);
		else
			CSeqTrans::ReverseComplement(ProbeLen,pProbeSeq);

// no changes to LowMMCnt or LowHitInstances?
if(*pLowMMCnt == LowMMCnt && *pLowHitInstances == LowHitInstances)
	{
	if(*pNxtLowMMCnt > NxtLowMMCnt)		// NxtLowMMCnt and LowHitInstances same but perhaps NxtLowMMCnt was reduced
		{
		*pNxtLowMMCnt    = NxtLowMMCnt;
		if((NxtLowMMCnt - *pLowMMCnt) < MMDelta)
			return(eHRMMDelta);			        // report hits but MMDelta criteria not met
		return(eHRRMMDelta);					// report that the only change was that the hamming delta was reduced but still accepted
		}
	return(eHRnone);					// return no change
	}

*pLowMMCnt       = LowMMCnt;
*pLowHitInstances = LowHitInstances;
*pNxtLowMMCnt    = NxtLowMMCnt;

if(*pLowHitInstances >= 1 && (*pNxtLowMMCnt - *pLowMMCnt) < MMDelta)
	return(eHRMMDelta);			// report hits but MMDelta criteria not met


if(*pLowHitInstances > MaxHits) // MMDelta requirement met but too many multiple hits?
   return(eHRHitInsts);

return(eHRhits);				// MMDelta requirement met and within the multiple hits limit
}



teBSFrsltCodes
CSfxArray::InitOverOccKMers(int KMerLen,	// will be processing for over occurring KMers of this length (max cMaxKmerLen)
				int MaxKMerOccs)		// which if there are more than MaxKMerOccs instances will be classified as an over-occurrence
{
m_MaxIter = m_MaxKMerOccs = MaxKMerOccs;
if(KMerLen < 1 || KMerLen > cMaxKmerLen || MaxKMerOccs < 1)
	return(eBSFerrParams);

SerialiseBaseFlags();
if(m_pOccKMerClas != NULL && KMerLen == m_OccKMerLen)
	{
	memset(m_pOccKMerClas,0,m_AllocOccKMerClasMem);
	ReleaseBaseFlags();
	return(eBSFSuccess);
	}

if(m_pOccKMerClas != NULL)
	{
	delete m_pOccKMerClas;
	m_pOccKMerClas = NULL;
	}

m_OccKMerLen = KMerLen;
m_AllocOccKMerClasMem = (size_t)0x01 << ((KMerLen-1) * 2);

#ifdef _WIN32
m_pOccKMerClas = (uint8_t *) malloc(m_AllocOccKMerClasMem);
if(m_pOccKMerClas==NULL)
	{
	AddErrMsg("CSfxArray::InitOverOccKMers","unable to allocate %zu bytes for over occurring KMers",m_AllocOccKMerClasMem);
	Reset(false);			// closes opened file..
	ReleaseBaseFlags();
	return(eBSFerrMem);
	}
#else
// gnu malloc is still in the 32bit world and seems to have issues if more than 2GB allocation
m_pOccKMerClas = (uint8_t *)mmap(NULL,m_AllocOccKMerClasMem, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
if(m_pOccKMerClas == MAP_FAILED)
	{
	AddErrMsg("CSfxArray::InitOverOccKMers","unable to allocate %zu bytes for over occurring KMers",m_AllocOccKMerClasMem);
	m_pOccKMerClas = NULL;
	Reset(false);			// closes opened file..
	ReleaseBaseFlags();
	return(eBSFerrMem);
	}
#endif
memset(m_pOccKMerClas,0,m_AllocOccKMerClasMem);
ReleaseBaseFlags();
return(eBSFSuccess);
}

int												// 0: no instances, 1: number occurrences <= m_MaxKMerOccs, 2: number occurrences > m_MaxKMerOccs
CSfxArray::OverOccKMerClas(int KMerLen,		// KMer length, checked and must be equal to m_OccKMerLen
				etSeqBase *pSeq)				// return over occurrence classification for this sequence
{
etSeqBase *pBase;
size_t PackedSeqIdx;
int ClasShf;
int OccKMerClas;
uint32_t NumCopies;
int Idx;

if(KMerLen > cMaxKmerLen || m_pOccKMerClas == NULL || KMerLen != m_OccKMerLen)		// NOTE - now classifying both over or none occurrences of all K-mers even if KMerLen > cMaxKmerLen
	{ 
	NumCopies = NumExactKMers(m_MaxKMerOccs, KMerLen, pSeq);
	if (NumCopies == 0)
		OccKMerClas = 0;
	else
		if (NumCopies <= m_MaxKMerOccs)
			OccKMerClas = 1;
		else
			OccKMerClas = 2;
	return(OccKMerClas);
	}

pBase = pSeq;
PackedSeqIdx = 0;
ClasShf = (*pBase++ & 0x03) * 2;
for(Idx = 1; Idx < KMerLen; Idx++,pBase++)
	{
	PackedSeqIdx <<= 2;
	PackedSeqIdx |= *pBase & 0x03;
	}
if((OccKMerClas = (m_pOccKMerClas[PackedSeqIdx] >> ClasShf) & 0x03) > 0)
	return(OccKMerClas - 1);

NumCopies = NumExactKMers(m_MaxKMerOccs,KMerLen,pSeq);
if(NumCopies == 0)
	OccKMerClas = 1;
else
	if(NumCopies <= m_MaxKMerOccs)	 
		OccKMerClas = 2;
	else
		OccKMerClas = 3;	

// NOTE: no serialisation of access to m_pOccKMerClas required as access is at the byte level with overwriting 2 bits only
// If another thread accesses the same byte then there may be a conflict but this may very occasionally result in two threads
// independently regenerating the same classification if that Kmer not previously classified
uint8_t Tmp = m_pOccKMerClas[PackedSeqIdx] & ~(0x03 << ClasShf);
Tmp |= (OccKMerClas << ClasShf); 
m_pOccKMerClas[PackedSeqIdx] = Tmp;
return(OccKMerClas - 1);
}


int 
CSfxArray::InitKMerCountDists(uint32_t MaxKMerCnts,		// // max KMer target count + 1 to allow for K-mers with 0 counts, range will be 0..MakKMerCnts-1 inclusive
	uint32_t *pKMerCntDist)						// record count distributions into this array sized MaxKMerCnts
{
m_MaxKMerCnts = MaxKMerCnts;
m_pKMerCntDist = pKMerCntDist;
memset(pKMerCntDist,0,sizeof(uint32_t) * MaxKMerCnts);
return(eBSFSuccess);
}


//
// LocateQuerySeqs
// Utilises a sliding core
int						// < 0 if errors, 0 if no matches or change, 1 if mumber of matches accepted, 2 MMDelta criteria not met, 3 too many match instances
CSfxArray::LocateQuerySeqs(uint32_t QuerySeqID,			// identifies this query sequence
	etSeqBase* pProbeSeq,			// probe
	uint32_t ProbeLen,				// probe length
	uint32_t CoreLen,				// core window length
	uint32_t CoreDelta,				// core window offset increment (1..n)
	eALStrand Align2Strand,		// align to this strand
	uint32_t MaxHits,				// (IN) process for at most this number of hits
	tsQueryAlignNodes* pHits,		// where to return hits (at most MaxHits)
	uint32_t CurMaxIter,				// max allowed iterations per subsegmented sequence when matching that subsegment
	int BaseMatchPts,				// award this many points for matching bases when extending 5' and 3' flanks
	int BaseMismatchPts)			// penalise this many points for mismatching bases when extending 5' and 3' flanks

{
	uint32_t CurCoreSegOfs;				// current core segment relative start
	uint32_t IterCnt;					// count iterator for current segment target matches
	uint32_t NumMatches;					// number of hit instances thus far
	char CurStrand;
	int64_t TargIdx;
	uint32_t NumCopies;
	int64_t TargSeqLeftIdx;

	uint32_t TargSeqID;
	uint32_t TargSeqLoci;
	uint32_t CurNumCoreSlides;
	uint32_t CurCoreDelta;
	tsQueryAlignNodes* pCurHit;
	tsQueryAlignNodes* pPrevHit;
	uint32_t SeqIdx;
	uint32_t Flank3Len;
	uint8_t TargBase;
	uint32_t Flank5Len;
	uint8_t ProbeBase;
	etSeqBase* pProbeBase;
	etSeqBase* pTargBase;
	uint32_t RptTargIdx;

	uint32_t TargHash[256];		// matched target identifiers and match sense are hashed: bit 0 == sense, bits 1..7 low order 0..6 of target identifier
	uint32_t HashIdx;

	etSeqBase* pTarg;			// target sequence
	void* pSfxArray;			// target sequence suffix array
	int64_t SfxLen;				// number of suffixs in pSfxArray
	tsSfxEntry* pEntry;

	// ensure suffix array block loaded for iteration!
	if (m_pSfxBlock == NULL || m_pSfxBlock->ConcatSeqLen == 0)
		return(eBSFerrInternal);

	pTarg = (etSeqBase*)&m_pSfxBlock->SeqSuffix[0];
	pSfxArray = (void*)&m_pSfxBlock->SeqSuffix[m_pSfxBlock->ConcatSeqLen];
	SfxLen = m_pSfxBlock->ConcatSeqLen;

	memset(TargHash, 0, sizeof(TargHash));

	NumMatches = 0;

	pCurHit = NULL;

	CurNumCoreSlides = 0;

	if (Align2Strand == eALSCrick)
	{
		CSeqTrans::ReverseComplement(ProbeLen, pProbeSeq);
		CurStrand = '-';
	}
	else
		CurStrand = '+';

	do
	{
		CurCoreDelta = CoreDelta;
		CurNumCoreSlides = 0;
		for (CurCoreSegOfs = 0;
			CurCoreSegOfs < (uint32_t)(ProbeLen - CoreLen);
			CurNumCoreSlides += 1,
			CurCoreSegOfs += CurCoreDelta)
		{
			if (((int)CurCoreSegOfs + CoreLen + CurCoreDelta) > ProbeLen)
				CurCoreDelta = ProbeLen - (CurCoreSegOfs + CoreLen);

			if (m_pOccKMerClas != NULL && CoreLen == m_OccKMerLen && (OverOccKMerClas(CoreLen, &pProbeSeq[CurCoreSegOfs]) != 1))	// 1 if number of K-mers of CoreLen is within range 1..max
				continue;													// outside of range so try next core

			// find first exactly matching and how many are exactly matching
			NumCopies = NumExactKMers(CurMaxIter, CoreLen, &pProbeSeq[CurCoreSegOfs], &TargIdx);

			// counting number of occurrences for each KMer core?
			if (m_pKMerCntDist != NULL && m_MaxKMerCnts > 0)
			{
				if (NumCopies >= m_MaxKMerCnts)
					NumCopies = m_MaxKMerCnts - 1;
#ifdef _WIN32
				InterlockedIncrement((volatile unsigned int*)&m_pKMerCntDist[NumCopies]);
#else
				__sync_fetch_and_add(&m_pKMerCntDist[NumCopies], 1);
#endif
			}
			if (NumCopies == 0 || NumCopies >= CurMaxIter || TargIdx == 0)	// if too many exact matches for this core then hopefully there will be a manageable number of matches on the next core
				continue;
			uint64_t LastTargIdx = TargIdx + NumCopies - 1;

			TargIdx -= 1;
			IterCnt = 0;
			while (IterCnt++ < NumCopies)
			{
				TargSeqLeftIdx = SfxOfsToLoci(m_pSfxBlock->SfxElSize, pSfxArray, TargIdx++);
				pEntry = MapChunkHit2Entry(TargSeqLeftIdx);
				TargSeqID = pEntry->EntryID;
				TargSeqLoci = (uint32_t)(TargSeqLeftIdx - pEntry->StartOfs);
				HashIdx = (TargSeqID & 0x07f) << 1 | (CurStrand == '+' ? 0 : 1);

				// have a putative core, check if this core overlaps with an already discovered core extension on to the same targeted sequence
				if ((SeqIdx = TargHash[HashIdx]) != 0)			// only bother checking for match contained within a previous if it is known there was at least one previous core match for this sequence 
				{
					do {
						pCurHit = &pHits[SeqIdx - 1];
						if (pCurHit->TargSeqID == TargSeqID)
						{
							// current core match is aligning to same target sequence and strand as a previously accepted aligned core match
							// if current core is contained within a prev alignment then continue checking
							if (CurCoreSegOfs < pCurHit->QueryStartOfs || ((CurCoreSegOfs + CoreLen) > (pCurHit->QueryStartOfs + pCurHit->AlignLen)))
								continue;
							if (TargSeqLoci < pCurHit->TargSeqLoci || ((TargSeqLoci + CoreLen) > (pCurHit->TargSeqLoci + pCurHit->AlignLen)))
								continue;
							if ((CurCoreSegOfs - pCurHit->QueryStartOfs) == (TargSeqLoci - pCurHit->TargSeqLoci))
								break;
						}
					} while ((SeqIdx = pCurHit->NxtHashMatch) != 0);
					if (SeqIdx != 0)			// current core part of a previous core extension, try for another targeted sequence
						continue;
				}

				// now maximally extend out on the flanks, scoring exacts/mismatches and retaining 5'/3' extensions with highest scores
				// an affine style score threshold is used when determining how many bp to extend a flank plus the mismatch rate must be no more than 15%
				// first try 5' extending
				int Num5MMs;
				int Num3MMs;
				int Extn5Score;
				int Extn3Score;
				int ExtnFloor;

				Num5MMs = 0;
				Num3MMs = 0;
				Extn5Score = CoreLen * BaseMatchPts;
				ExtnFloor = Extn5Score - (4 * BaseMatchPts); // on 5' extensions, when floor is reached then at least 4 exact matches required to recover back to a potential new high scoring extension
				int HiScoreScore = Extn5Score;
				int HiScoreExtnLen = 0;
				int HiScoreMM = 0;
				Flank5Len = 0;
				pTargBase = &pTarg[TargSeqLeftIdx - 1];
				pProbeBase = &pProbeSeq[CurCoreSegOfs - 1];
				for (SeqIdx = min(TargSeqLoci, CurCoreSegOfs); SeqIdx > 0; SeqIdx--, pTargBase--, pProbeBase--)
				{
					TargBase = *pTargBase & 0x07;
					ProbeBase = *pProbeBase & 0x07;
					if (TargBase > eBaseN || ProbeBase > eBaseN)		// mustn't match inter targ or probe sequences
						break;
					if (TargBase == eBaseN || ProbeBase == eBaseN || ProbeBase != TargBase)
					{
						if (Extn5Score >= ExtnFloor)
						{
							if (Extn5Score > ExtnFloor)
								Extn5Score -= BaseMismatchPts;		// heavily penalise mismatches until floor, aka like affine gap scoring
							if (Extn5Score < ExtnFloor)
								Extn5Score = ExtnFloor;
						}
						Num5MMs += 1;
					}
					else
						Extn5Score += BaseMatchPts;			// exact matches are rewarded
					Flank5Len += 1;
					if (Extn5Score >= HiScoreScore && ((Num5MMs * 100) / (Flank5Len + CoreLen)) <= 15) // 15% max mismatch rate allowed on extension
					{
						HiScoreScore = Extn5Score;
						HiScoreExtnLen = Flank5Len;
						HiScoreMM = Num5MMs;
						ExtnFloor = HiScoreScore - (4 * BaseMatchPts);
					}
				}
				Flank5Len = HiScoreExtnLen;
				Num5MMs = HiScoreMM;

				// now try 3' flank extending
				Flank3Len = 0;
				Num3MMs = 0;
				Extn3Score = CoreLen * BaseMatchPts;
				ExtnFloor = Extn3Score - (4 * BaseMatchPts); // on 3' extensions, when floor is reached then at least 4 exact matches required to recover back to a potential new high scoring extension
				HiScoreScore = Extn3Score;
				HiScoreExtnLen = 0;
				HiScoreMM = 0;

				pTargBase = &pTarg[TargSeqLeftIdx + CoreLen];
				pProbeBase = &pProbeSeq[CurCoreSegOfs + CoreLen];
				for (SeqIdx = CurCoreSegOfs + (uint32_t)CoreLen, RptTargIdx = TargSeqLoci + (uint32_t)CoreLen; SeqIdx < ProbeLen && RptTargIdx < pEntry->SeqLen; SeqIdx++, RptTargIdx++, pTargBase++, pProbeBase++)
				{
					TargBase = *pTargBase & 0x07;
					ProbeBase = *pProbeBase & 0x07;
					if (TargBase > eBaseN || ProbeBase > eBaseN)		// mustn't match inter targ or probe sequences
						break;

					if (TargBase == eBaseN || ProbeBase == eBaseN || ProbeBase != TargBase)
					{
						if (Extn3Score > ExtnFloor)
							Extn3Score -= BaseMismatchPts;		// heavily penalise mismatches until floor, aka like affine gap scoring
						if (Extn3Score < ExtnFloor)
							Extn3Score = ExtnFloor;
						Num3MMs += 1;
					}
					else
						Extn3Score += BaseMatchPts;			// exact matches are rewarded
					Flank3Len += 1;
					if (Extn3Score >= HiScoreScore && ((Num3MMs * 100) / (Flank3Len + CoreLen)) <= 15) // 15% max mismatch rate allowed on extension
					{
						HiScoreScore = Extn3Score;
						HiScoreExtnLen = Flank3Len;
						HiScoreMM = Num3MMs;
						ExtnFloor = HiScoreScore - (4 * BaseMatchPts);
					}
				}
				Flank3Len = HiScoreExtnLen;
				Num3MMs = HiScoreMM;

				// have a maximally extended core
				pCurHit = &pHits[NumMatches];
				pCurHit->AlignNodeID = NumMatches;
				pCurHit->QueryID = QuerySeqID;
				pCurHit->QueryStartOfs = CurCoreSegOfs - Flank5Len;
				pCurHit->TargSeqLoci = (uint32_t)(TargSeqLoci - Flank5Len);
				pCurHit->AlignLen = Flank5Len + CoreLen + Flank3Len;
				pCurHit->NumMismatches = Num5MMs + Num3MMs;
				pCurHit->FlgStrand = CurStrand == '+' ? 0 : 1;
				pCurHit->TargSeqID = TargSeqID;

				pCurHit->Flg2Rpt = 0;
				pCurHit->FlgFirst2tRpt = 0;
				pCurHit->FlgScored = 0;
				pCurHit->FlgRedundant = 0;
				pCurHit->HiScorePathNextIdx = 0;
				pCurHit->HiScore = 0;

				tsQueryAlignNodes* pPrevHit;
				if ((SeqIdx = TargHash[HashIdx]) != 0)
				{
					do {
						pPrevHit = &pHits[SeqIdx - 1];
						if (!pPrevHit->FlgRedundant && (pPrevHit->FlgStrand == pCurHit->FlgStrand && pPrevHit->TargSeqID == pCurHit->TargSeqID))
						{
							// does latest (pCurHit) query aligned subsequence intersect with the PrevHit?
							// if so then choose longest overlap, mark other as being redundant so it can later be discarded
							// note that 1bp non-overlap is being treated as if an overlap
							if (!((pCurHit->QueryStartOfs + pCurHit->AlignLen + CoreLen - 1) < pPrevHit->QueryStartOfs ||
								pCurHit->QueryStartOfs > (pPrevHit->QueryStartOfs + pPrevHit->AlignLen + CoreLen - 1)))
							{
								// check for obvious containment of a much smaller core within a larger core, the much smaller core is discarded
								if ((pCurHit->AlignLen + CoreLen / 2) <= pPrevHit->AlignLen)
									pCurHit->FlgRedundant = true;
								else
									if (pCurHit->AlignLen >= (pPrevHit->AlignLen + CoreLen / 2))
										pPrevHit->FlgRedundant = true;

								// check for intersect on target
								if (!((pCurHit->TargSeqLoci + pCurHit->AlignLen + CoreLen - 1) < pPrevHit->TargSeqLoci ||
									pCurHit->TargSeqLoci > (pPrevHit->TargSeqLoci + pPrevHit->AlignLen + CoreLen - 1)))
								{
									// have intersect on both query and target - need to choose!
									// currently choosing the longest alignment length and
									// as a tiebreaker then the lowest number of mismatches followed by discarding pCurHit
									if (pCurHit->AlignLen < pPrevHit->AlignLen)
										pCurHit->FlgRedundant = true;
									else
										if (pCurHit->AlignLen > pPrevHit->AlignLen)
											pPrevHit->FlgRedundant = true;
										else
											if (pCurHit->NumMismatches >= pPrevHit->NumMismatches)
												pCurHit->FlgRedundant = true;
											else
												pPrevHit->FlgRedundant = true;

								}
							}
						}
					} while ((SeqIdx = pPrevHit->NxtHashMatch) != 0);
				}
				if (pCurHit->FlgRedundant)
					continue;
				NumMatches += 1;
				pCurHit->NxtHashMatch = TargHash[HashIdx];
				TargHash[HashIdx] = NumMatches;
				if (NumMatches >= MaxHits)
					break;
			}
			if (NumMatches >= MaxHits)
			{
				Align2Strand = eALSnone;					// can't improve so early exit
				break;
			}
		}
		if (CurStrand == '+' && Align2Strand == eALSboth)	// if just processed watson '+' strand then will need to process crick or '-' strand if processing both strands
		{
			CSeqTrans::ReverseComplement(ProbeLen, pProbeSeq);
			CurStrand = '-';
			Align2Strand = eALSCrick;
		}
		else									// either processed both or just the crick strand so no further processing required
			Align2Strand = eALSnone;			// two passes max - first on the '+' strand with optional 2nd on '-' strand
	} while (!(NumMatches >= MaxHits) && Align2Strand != eALSnone);

	if (CurStrand == '-')									// restore probe sequence if had started processing '-' strand
		CSeqTrans::ReverseComplement(ProbeLen, pProbeSeq);

	// iterate and remove all alignment nodes marked as being redundant
	uint32_t NonRedundant = 0;
	if (NumMatches > 0)
	{
#if _ALIGNNODESDIAG_
		pCurHit = pPrevHit = pHits;
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "Query (%.8d) has %.5d matches", pCurHit->QueryID, NumMatches);
		for (SeqIdx = 0; SeqIdx < NumMatches; SeqIdx++, pCurHit++)
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "P %c,%.3d,%.3d,%.5d", pCurHit->FlgStrand ? '-' : '+', pCurHit->QueryStartOfs, pCurHit->AlignLen, pCurHit->TargSeqLoci);
#endif
		if (NumMatches == 1)
			return(pHits->FlgRedundant ? 0 : 1);
		pCurHit = pPrevHit = pHits;
		for (SeqIdx = 0; SeqIdx < NumMatches; SeqIdx++, pCurHit++)
		{
			if (pCurHit->FlgRedundant)
				continue;
			if (pPrevHit != pCurHit)
				*pPrevHit = *pCurHit;
			pPrevHit++;
			NonRedundant++;
		}
#if _ALIGNNODESDIAG_
		pCurHit = pPrevHit = pHits;
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "Query (%.8d) has %.5d non-redundant matches", pCurHit->QueryID, NonRedundant);
		for (SeqIdx = 0; SeqIdx < NonRedundant; SeqIdx++, pCurHit++)
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "F %c,%.3d,%.3d,%.5d", pCurHit->FlgStrand ? '-' : '+', pCurHit->QueryStartOfs, pCurHit->AlignLen, pCurHit->TargSeqLoci);
#endif
	}
	return(NonRedundant);				// number of alignment nodes non-redundant, MMDelta requirement met and within the multiple hits limit
}


// LocateBestMatches
// Locate, at most MaxHits, alignments having no more than MaxTotMM mismatches, additional matches are sloughed
// 
int						// < 0 if errors, 0 if no matches, 1..MaxHits, or MaxHits+1 if additional matches have been sloughed
CSfxArray::LocateBestMatches(uint32_t ReadID,			// identifies this read
						 int MaxTotMM,			        // return matches having at most this number of mismatches
						 int CoreLen,					// core window length
						 int CoreDelta,					// core window offset increment (1..n)
						 int MaxNumCoreSlides,			// max number of times to slide core on each strand
						 eALStrand Align2Strand,		// align to this strand
 						 etSeqBase *pProbeSeq,			// probe
						 int ProbeLen,					// probe length
						 int MaxHits,					// return at most this number of best hits
						 int *pHitInstances,			// returned number of match instances in pHits
						 tsHitLoci *pHits,				// where to return best hits (at most MaxHits)
						 int CurMaxIter,				// max allowed iterations per subsegmented sequence when matching that subsegment
						 int NumAllocdIdentNodes,		// memory has been allocated by caller for holding up to this many tsIdentNodes
						 tsIdentNode *pAllocsIdentNodes) // memory allocated by caller for holding tsIdentNodes
{
int CurCoreSegOfs;				// current core segment relative start
int IterCnt;					// count iterator for current segment target matches
int CurMMCnt;					// total number of mismatches for current target sequence being processed
int CoreMMCnt;					// total number of mismatches for current target core being processed
int LowHitInstances;			// number of hit instances thus far for hits with LowMMCnt mismatches
etSeqBase *pProbeBase;
etSeqBase *pTargBase;
char CurStrand;
int64_t TargIdx;

tsIdentNode *pHashArray[cHashEntries+1];		// hash array holding ptrs to identifier nodes
tsIdentNode *pIdentNodes = pAllocsIdentNodes;	// identifier nodes
tsIdentNode *pCurIdentNode;
tsIdentNode *pNewIdentNode;
int CurNumIdentNodes;
int Hash;

uint32_t TargSeqID;
uint32_t NumTargSeqProc;

bool bFirstIter;				// set false after the first subsequence core returned by LocateFirstExact has been processed

int64_t LastTargIdx;
uint32_t NumCopies;

uint32_t PatIdx;
int64_t TargSeqLeftIdx;
int64_t ProbeSeqLeftIdx;
int TargMatchLen;

etSeqBase BisBase;
uint8_t ProbeBase;
uint8_t TargBase;

int Cmp;
int CurNumCoreSlides;
int CurCoreDelta;
tsHitLoci *pCurHit;

bool bMatchesSloughed;

etSeqBase *pTarg;			// target sequence
void *pSfxArray;			// target sequence suffix array
int64_t SfxLen;				// number of suffixs in pSfxArray
tsSfxEntry *pEntry;
BisBase=eBaseN;

// ensure suffix array block loaded for iteration!
if(m_pSfxBlock == NULL || m_pSfxBlock->ConcatSeqLen == 0)
	return(eBSFerrInternal);

pTarg = (etSeqBase *)&m_pSfxBlock->SeqSuffix[0];
pSfxArray = (void *)&m_pSfxBlock->SeqSuffix[m_pSfxBlock->ConcatSeqLen];
SfxLen = m_pSfxBlock->ConcatSeqLen;


LowHitInstances = 0;
NumTargSeqProc = 0;
CurNumIdentNodes = 0;
CurNumCoreSlides = 0;

if(pHitInstances != NULL)
	*pHitInstances = 0;
bMatchesSloughed = false;

if(Align2Strand == eALSCrick)
	{
	if(m_bColorspace)
		CSeqTrans::ReverseSeq(ProbeLen,pProbeSeq);
	else
		CSeqTrans::ReverseComplement(ProbeLen,pProbeSeq);
	CurStrand = '-';
	}
else
	CurStrand = '+';

do
	{
	CurCoreDelta = CoreDelta;
	CurNumCoreSlides = 0;
	memset(pHashArray,0,sizeof(pHashArray));
	CurNumIdentNodes = 0;
	for(CurCoreSegOfs = 0;
		CurNumCoreSlides < MaxNumCoreSlides &&
		CurCoreSegOfs <= (ProbeLen - CoreLen) &&
		CurCoreDelta > CoreLen/3 &&
		CurNumIdentNodes < NumAllocdIdentNodes;	// can only allow upto cMaxNumIdentNodes to be saved
	    CurNumCoreSlides += 1,
	    CurCoreSegOfs += CurCoreDelta)
		{
		if(CurNumCoreSlides >= MaxNumCoreSlides)
			break;
		if((CurCoreSegOfs + CoreLen + CurCoreDelta) > ProbeLen)
			CurCoreDelta = ProbeLen - (CurCoreSegOfs + CoreLen);

		TargIdx = LocateFirstExact(&pProbeSeq[CurCoreSegOfs],CoreLen,pTarg,m_pSfxBlock->SfxElSize,pSfxArray,0,0,SfxLen-1);
		if(TargIdx == 0)        // 0 if no core segment matches
			continue;			// try for match on next core segment after shifting core to right

		TargIdx -= 1;
		IterCnt = 0;
		NumCopies = 0;
		bFirstIter = true;		// set false after the first subsequence core returned by LocateFirstExact has been processed
		while(!CurMaxIter || IterCnt < CurMaxIter)
			{
			if(CurNumIdentNodes >= NumAllocdIdentNodes)
				break;

			if(!bFirstIter) {

				// ensure not about to iterate past end of suffix array!
				if((TargIdx + 1) >= (int64_t)m_pSfxBlock->ConcatSeqLen || (SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,TargIdx+1) + CoreLen) >  (int64_t)m_pSfxBlock->ConcatSeqLen)
					break;

				if(IterCnt == 100 && !NumCopies)
					{
					// check how many more exact copies there are of the current probe subsequence, if too many then don't bother exploring these
					LastTargIdx = LocateLastExact(&pProbeSeq[CurCoreSegOfs],CoreLen,pTarg,m_pSfxBlock->SfxElSize,pSfxArray,0,TargIdx-1,SfxLen-1);
					NumCopies = LastTargIdx > 0 ? (uint32_t)(1 + LastTargIdx - TargIdx) : 0;
					if(CurMaxIter && NumCopies > (uint32_t)CurMaxIter)		// only checking at the 100th iteration allows a little slack
						break;										// try next core segment
					}

				// check that this new putative core is still matching
				pTargBase = &pTarg[SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,TargIdx+1)];
				pProbeBase = &pProbeSeq[CurCoreSegOfs];
				if(m_bBisulfite)
					Cmp = BSCmpProbeTarg(pProbeBase,pTargBase,CoreLen);
				else
					{
					int Ofs;
					uint8_t El1;
					uint8_t El2;
					etSeqBase *pEl1= pProbeBase;
					etSeqBase *pEl2 = pTargBase;
					Cmp = 0;
					for(Ofs=0; Ofs < CoreLen; Ofs++)
						{
						El2 = *pEl2++ & 0x0f;
						if(El2 == eBaseEOS)
							{
							Cmp = -1;
							break;
							}
						El1 = *pEl1++ & 0x0f;
						if(El1 > El2)
							{
							Cmp = 1;
							break;
							}
						if(El1 < El2)
							{
							Cmp = -1;
							break;
							}
						}
					}

				if(Cmp != 0)				// will be non-zero if target no longer matches
					break;					// try next core segment
				TargIdx += 1;
				}

			bFirstIter = false;
			if(SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,TargIdx) < (uint32_t)CurCoreSegOfs)
				continue;

			TargMatchLen = ProbeLen;
			TargSeqLeftIdx = SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,TargIdx) - CurCoreSegOfs;
			ProbeSeqLeftIdx = 0;

			// ensure comparisons are still within start/end range of target sequence/assembly
			if(!TargMatchLen || ((uint64_t)TargSeqLeftIdx + (uint32_t)TargMatchLen) > m_pSfxBlock->ConcatSeqLen)
				continue;

			// check if target already processed
			TargSeqID = (uint32_t)(1 + SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,TargIdx) - (uint32_t)CurCoreSegOfs);
			Hash = (TargSeqID & cHashEntries);
			if((pCurIdentNode = pHashArray[Hash])==NULL)
				{
				pCurIdentNode = pHashArray[Hash] = &pIdentNodes[CurNumIdentNodes++];
				pCurIdentNode->TargSeqID = TargSeqID;
				pCurIdentNode->pNxt = NULL;
				}
			else
				{
				while(pCurIdentNode != NULL && pCurIdentNode->TargSeqID != TargSeqID)
					pCurIdentNode = pCurIdentNode->pNxt;
				if(pCurIdentNode != NULL)							// will be non-null if have already processed this sequence
					continue;
				pNewIdentNode = &pIdentNodes[CurNumIdentNodes++];
				pNewIdentNode->TargSeqID = TargSeqID;
				pNewIdentNode->pNxt = pHashArray[Hash];
				pHashArray[Hash] = pNewIdentNode;
				}
			NumTargSeqProc += 1;
			IterCnt += 1;

			// now do the matching allowing for mismatches
			// if colorspace then need to allow for adjacent apparent mismatches actually being the result of a single substitution relative to the target
			// if colorspace and a standalone, without an adjacent, mismatch then this is supposedly a sequencing error, still treat as though a mismatch
				{
				pTargBase = &pTarg[TargSeqLeftIdx];
				pProbeBase = pProbeSeq;
				CurMMCnt = 0;
				CoreMMCnt = 0;
				bool bPairMM = false;
				if(m_bBisulfite)
					BisBase = GetBisBase(TargMatchLen,pTargBase,pProbeBase);
				for(PatIdx = 0; PatIdx < (uint32_t)TargMatchLen; PatIdx++,pTargBase++,pProbeBase++)
					{
					TargBase = *pTargBase & 0x0f;
					ProbeBase = *pProbeBase & 0x0f;
					if(TargBase == eBaseEOS)		// mustn't match across entry sequences
						break;

					if(m_bColorspace)		// in colorspace, unpaired mismatches are sequencer errors but still treat as if a substitution
						{
						if(!bPairMM && TargBase != ProbeBase)
							{
							if(PatIdx < ((uint32_t)TargMatchLen-1))
								{
								if((pTargBase[1] & 0x0f) == (pProbeBase[1] & 0x0f))
									{
									if(++CurMMCnt > MaxTotMM)
										break;
									continue;
									}
								}
							// accept as being a mismatch
							bPairMM = true;
							}
						else
							{
							bPairMM = false;
							continue;
							}
						}
					else				// in basespace
						{
						if(ProbeBase == TargBase)
							continue;

						if(m_bBisulfite)
							{
							if(TargBase == eBaseC || TargBase == eBaseG)
								{
								switch(BisBase) {
									case eBaseA:		// only allow A
										if(ProbeBase == eBaseA && TargBase == eBaseG)
											continue;
										break; // mismatch
									case eBaseT:		// only allow T
										if(ProbeBase == eBaseT && TargBase == eBaseC)
											continue;
										break; // mismatch
									}
								}
							}
						}

					// execution here only if mismatch
					if(++CurMMCnt > MaxTotMM)
						break;
					}
				if(PatIdx != TargMatchLen)
					continue;

				// need to ensure only the hits with lowest number of mismatches are being retained
				pCurHit = NULL;
				if(LowHitInstances)
					{
					int BestIdx;
					tsHitLoci *pNxtHit;
					if(LowHitInstances == MaxHits)	// if already at max allowed to be returned then at least one is going to be sloughed
						bMatchesSloughed = true;
					pNxtHit = &pHits[0]; 
					for(BestIdx = 0; BestIdx < LowHitInstances; BestIdx++,pNxtHit++)
						{
						if(pNxtHit->Seg[0].Mismatches > CurMMCnt)		// if this hit has more missmatches than current target sequence then insert
							{
							pCurHit = pNxtHit;
							if((BestIdx + 1) < MaxHits)
								memmove(&pNxtHit[1],pNxtHit,sizeof(tsHitLoci) * (LowHitInstances - BestIdx));
							break;
							}
						}
					if(BestIdx == LowHitInstances && LowHitInstances < MaxHits)	// not inserting, may be able to append
						pCurHit = pNxtHit;
					}
				else
					pCurHit = &pHits[0];			// must have been first hit

				if(pCurHit != NULL)
					{
					pEntry = MapChunkHit2Entry(TargSeqLeftIdx);
					pCurHit->FlgChimeric = 0;
					pCurHit->FlgInDel = 0;
					pCurHit->FlgInsert = 0;
					pCurHit->FlgSplice = 0;
					pCurHit->FlgNonOrphan = 0;
					memset(pCurHit->Seg,0,sizeof(pCurHit->Seg));
					pCurHit->Seg[0].Strand = CurStrand;
					pCurHit->Seg[0].ChromID = pEntry->EntryID;
					pCurHit->Seg[0].MatchLoci = (uint32_t)(TargSeqLeftIdx - pEntry->StartOfs);
					pCurHit->Seg[0].MatchLen = ProbeLen;
					pCurHit->Seg[0].Mismatches = CurMMCnt;
					pCurHit->Seg[0].TrimMismatches = CurMMCnt;
					pCurHit->BisBase = BisBase;
					if(LowHitInstances < MaxHits)
						LowHitInstances += 1;
					else
						{
						pCurHit = &pHits[LowHitInstances-1]; 
						MaxTotMM = pCurHit->Seg[0].Mismatches;
						}
					}
				}
			}
		if(LowHitInstances == MaxHits && MaxTotMM == 0 && !bMatchesSloughed)
			{
			Align2Strand = eALSnone;					// can't improve so early exit
			break;
			}
		}
	if(CurStrand == '+' && Align2Strand == eALSboth)	// if just processed watson '+' strand then will need to process crick or '-' strand if processing both strands
		{
		if(m_bColorspace)
			CSeqTrans::ReverseSeq(ProbeLen,pProbeSeq);
		else
			CSeqTrans::ReverseComplement(ProbeLen,pProbeSeq);
		CurStrand = '-';
		Align2Strand = eALSCrick;
		NumTargSeqProc = 0;
		}
	else									// either processed both or just the crick strand so no further processing required
		Align2Strand = eALSnone;			// two passes max - first on the '+' strand with optional 2nd on '-' strand
	}
while(!(LowHitInstances == MaxHits && MaxTotMM == 0 && !bMatchesSloughed) && Align2Strand != eALSnone);

if(CurStrand == '-')									// restore probe sequence if had started processing '-' strand
	if(m_bColorspace)
		CSeqTrans::ReverseSeq(ProbeLen,pProbeSeq);
	else
		CSeqTrans::ReverseComplement(ProbeLen,pProbeSeq);
if(pHitInstances != NULL)
	*pHitInstances = LowHitInstances;
if(LowHitInstances == 0)
	return(0);			// no hits
return(bMatchesSloughed ? LowHitInstances+1 : LowHitInstances);		// accepted hits and may have sloughed some
}

int						// < 0 if errors, eHRnone if no matches or too many, eHRhits if mumber of matches accepted
CSfxArray::LocateSpliceJuncts(uint32_t ExtdProcFlags,			// flags indicating if lower levels need to do any form of extended processing with this specific read...
						 uint32_t ReadID,					// identifies this read
						 int MaxSpliceJunctLen,		// junction has to be no longer than this length
						 int MaxTotMM,					// max number of mismatches allowed
						 int CoreLen,					// core window length
						 eALStrand Align2Strand,		// align to this strand
  						 int *pLowHitInstances,			// In/Out number of match instances for lowest number of mismatches thus far for this read
						 int *pLowMMCnt,				// In/Out lowest number of mismatches thus far for this read
						 int *pNxtLowMMCnt,				// In/Out next to lowest number of mismatches thus far for this read
						 etSeqBase *pProbeSeq,			// probe
						 int ProbeLen,					// probe length
						 int MaxHits,					// (IN) process for at most this number of hits (currently will always only allow a single hit)
						 tsHitLoci *pHits,				// where to return any hit
 						 int *pScore,					// where to return score
						 int CurMaxIter,				// max allowed iterations per subsegmented sequence when matching that subsegment
						 int NumAllocdIdentNodes,		// memory has been allocated by caller for holding upto this many tsIdentNodes
						 tsIdentNode *pAllocsIdentNodes) // memory allocated by caller for holding tsIdentNodes
{
int Phase;						// current processing phase: 0 - look for InDel anchored 5' to 3', 1 - look for splice junction anchored 3' to 5'
int CurCoreSegOfs;				// current core segment relative start
int IterCnt;					// count iterator for current segment target matches
etSeqBase *pProbeBase;
etSeqBase *pTargBase;
char CurStrand;
int64_t TargIdx;

uint32_t NumTargSeqProc;

bool bFirstIter;				// set false after the first subsequence core returned by LocateFirstExact has been processed


int64_t TargSeqLeftIdx;
int64_t ProbeSeqLeftIdx;
int TargMatchLen;

int SpliceJunctLenLimit;
uint32_t NumCopies;

etSeqBase BisBase;

tsHitLoci RightHitLoci;
tsHitLoci LeftHitLoci;

int RsltRightSplice;
int RsltLeftSplice;
int BestScoreInstances;

int Cmp;
tsHitLoci *pCurHit;
tsSfxEntry *pEntry;
etSeqBase *pTarg;			// target sequence
void *pSfxArray;			// target sequence suffix array
int64_t SfxLen;				// number of suffixs in pSfxArray
BisBase=eBaseN;

// ensure suffix array block loaded for iteration!
if(m_pSfxBlock == NULL || m_pSfxBlock->ConcatSeqLen == 0)
	return(eBSFerrInternal);

pTarg = (etSeqBase *)&m_pSfxBlock->SeqSuffix[0];
pSfxArray = (void *)&m_pSfxBlock->SeqSuffix[m_pSfxBlock->ConcatSeqLen];
SfxLen = m_pSfxBlock->ConcatSeqLen;

if(MaxTotMM > cMaxJunctAlignMM)		// silently clamp
	MaxTotMM = cMaxJunctAlignMM;

pCurHit = NULL;

NumTargSeqProc = 0;
*pLowHitInstances = 0;
*pLowMMCnt = 0;
*pNxtLowMMCnt = 0;
memset(pHits,0,sizeof(tsHitLoci));
BestScoreInstances = 0;

if(Align2Strand == eALSCrick)
	{
	if(m_bColorspace)
		CSeqTrans::ReverseSeq(ProbeLen,pProbeSeq);
	else
		CSeqTrans::ReverseComplement(ProbeLen,pProbeSeq);
	CurStrand = '-';
	}
else
	CurStrand = '+';

do
	{
	Phase = 0;
	for(Phase = 0; Phase < 2; Phase+=1)
		{
		switch(Phase) {
			case 0:
				CurCoreSegOfs = 0;
				break;
			case 1:
				CurCoreSegOfs = ProbeLen - CoreLen;
				break;
			}

		if (m_pOccKMerClas != NULL && CoreLen == m_OccKMerLen && (OverOccKMerClas(CoreLen, &pProbeSeq[CurCoreSegOfs]) != 1))	// 1 if number of K-mers of CoreLen is within range 1..max
			continue;

		TargIdx = LocateFirstExact(&pProbeSeq[CurCoreSegOfs],CoreLen,pTarg,m_pSfxBlock->SfxElSize,pSfxArray,0,0,SfxLen-1);
		if(TargIdx == 0)        // 0 if no core segment matches
			continue;			// try for match on next core segment after shifting core to right
		TargIdx -= 1;
		IterCnt = 0;
		NumCopies = 0;
		bFirstIter = true;		// set false after the first subsequence core returned by LocateFirstExact has been processed
		while(!CurMaxIter || IterCnt < CurMaxIter)
			{
			if(!bFirstIter) {

				// ensure not about to iterate past end of suffix array!
				if((TargIdx + 1) >= (int64_t)m_pSfxBlock->ConcatSeqLen || (SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,TargIdx+1) + (Phase == 0 ? ProbeLen : CoreLen)) >=  (int64_t)m_pSfxBlock->ConcatSeqLen)
					break;

				// check that this new putative core is still matching
				pTargBase = &pTarg[SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,TargIdx+1)];
				pProbeBase = &pProbeSeq[CurCoreSegOfs];
				if(m_bBisulfite)
					Cmp = BSCmpProbeTarg(pProbeBase,pTargBase,CoreLen);
				else
					{
					int Ofs;
					uint8_t El1;
					uint8_t El2;
					etSeqBase *pEl1= pProbeBase;
					etSeqBase *pEl2 = pTargBase;
					Cmp = 0;
					for(Ofs=0; Ofs < CoreLen; Ofs++)
						{
						El2 = *pEl2++ & 0x0f;
						if(El2 == eBaseEOS)
							{
							Cmp = -1;
							break;
							}
						El1 = *pEl1++ & 0x0f;
						if(El1 > El2)
							{
							Cmp = 1;
							break;
							}
						if(El1 < El2)
							{
							Cmp = -1;
							break;
							}
						}
					}

				if(Cmp != 0)				// will be non-zero if target no longer matches
					break;					// try next core segment
				TargIdx += 1;
				}

			bFirstIter = false;

			if(SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,TargIdx) < (uint32_t)CurCoreSegOfs)
				continue;

			TargMatchLen = ProbeLen;
			TargSeqLeftIdx = SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,TargIdx) - CurCoreSegOfs;
			ProbeSeqLeftIdx = 0;

			// ensure comparisons are still within start/end range of target sequence/assembly
			if((TargSeqLeftIdx + ProbeLen) >= (int64_t)m_pSfxBlock->ConcatSeqLen)
				continue;

			pEntry = MapChunkHit2Entry(SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,TargIdx));
			if(pEntry == NULL)
				continue;
			if(TargSeqLeftIdx < (int64_t)pEntry->StartOfs || ((TargSeqLeftIdx +  ProbeLen) > (int64_t)pEntry->EndOfs))
				continue;

			NumTargSeqProc += 1;
			IterCnt += 1;

			pTargBase = &pTarg[TargSeqLeftIdx];
			pProbeBase = pProbeSeq;

				// explore for a splice junction to right if matched on core at 5' end of probe (phase 0)
			if(Phase == 0)
				{
				SpliceJunctLenLimit = (int)(m_pSfxBlock->ConcatSeqLen - TargSeqLeftIdx);
				if(SpliceJunctLenLimit > (cMinJunctAlignSep + cMinJunctSegLen))
					{
					SpliceJunctLenLimit -= (cMinJunctAlignSep + cMinJunctSegLen);
					if(SpliceJunctLenLimit > MaxSpliceJunctLen)
						SpliceJunctLenLimit = MaxSpliceJunctLen;

					RsltRightSplice = ExploreSpliceRight(CurStrand,			// aligning on this strand
						SpliceJunctLenLimit,				// any junction not allowed to be longer than this
						MaxTotMM,							// can be at most this many mismatches in total
						CoreLen,							// core length used
						ProbeLen,							// length of probe excluding any eBaseEOS
						pProbeBase,							// pts to probe sequence
						TargSeqLeftIdx,						// pTarg corresponds to this suffix index
 						m_pSfxBlock->ConcatSeqLen,				// total length of target
						pTargBase,							// pts to target sequence
						&RightHitLoci);						// where to return alignment loci

					if(RsltRightSplice > 0 &&  RightHitLoci.Score >= pHits->Score)
						{
						if(RightHitLoci.Score == pHits->Score)
							{
							if(pHits->Seg[0].MatchLoci == RightHitLoci.Seg[0].MatchLoci)
								continue;
							if(++BestScoreInstances > MaxHits)
								continue;
							}
						else
							BestScoreInstances = 0;
						pHits[BestScoreInstances++] = RightHitLoci;
						}
					}
				}

			// explore for a splice junction to the left if matched on core at 3' end of probe (phase 1)
			if(Phase == 1 && (TargSeqLeftIdx >= (uint32_t)(CurCoreSegOfs + cMinJunctSegLen)))
				{
				SpliceJunctLenLimit = min((int32_t)TargSeqLeftIdx, (int32_t)MaxSpliceJunctLen);
				if(SpliceJunctLenLimit >= (cMinJunctAlignSep + cMinJunctSegLen))
					{
					SpliceJunctLenLimit -= cMinJunctSegLen;

					RsltLeftSplice = ExploreSpliceLeft(CurStrand,							// aligning on this strand
						SpliceJunctLenLimit,		// any junction has to be no longer than this length
						MaxTotMM,				// can be at most this many mismatches in total
						CoreLen,					// core length used
						ProbeLen,				// length of probe excluding any eBaseEOS
						pProbeBase,				// pts to probe sequence
						TargSeqLeftIdx,	        // pTargBase corresponds to this suffix index
						m_pSfxBlock->ConcatSeqLen,		// total length of target
						pTargBase,				// pts to target sequence
						&LeftHitLoci);			// where to return alignment loci

					if(RsltLeftSplice > 0 &&  LeftHitLoci.Score >= pHits->Score)
						{
						if(LeftHitLoci.Score == pHits->Score)
							{
							if(pHits->Seg[0].MatchLoci == LeftHitLoci.Seg[0].MatchLoci)
								continue;
							if(++BestScoreInstances > MaxHits)
								continue;
							}
						else
							BestScoreInstances = 0;;
						pHits[BestScoreInstances++] = LeftHitLoci;
						}
					}
				}
			}
		if(BestScoreInstances >= 1 && pHits->Score >= cMaxScore) // unlikely to do any better than this so early exit...
			{
			Align2Strand = eALSnone;
			break;
			}
		}


	if(CurStrand == '+' && Align2Strand == eALSboth)	// if just processed watson '+' strand then will need to process crick or '-' strand if processing both strands
		{
		if(m_bColorspace)
			CSeqTrans::ReverseSeq(ProbeLen,pProbeSeq);
		else
			CSeqTrans::ReverseComplement(ProbeLen,pProbeSeq);
		CurStrand = '-';
		Align2Strand = eALSCrick;
		NumTargSeqProc = 0;
		}
	else									// either processed both or just the crick strand so no further processing required
		Align2Strand = eALSnone;			// two passes max - first on the '+' strand with optional 2nd on '-' strand
	}
while(!(BestScoreInstances >= 1 && pHits->Score >= cMaxScore) && Align2Strand != eALSnone);

if(CurStrand == '-')									// restore probe sequence if had started processing '-' strand
		if(m_bColorspace)
			CSeqTrans::ReverseSeq(ProbeLen,pProbeSeq);
		else
			CSeqTrans::ReverseComplement(ProbeLen,pProbeSeq);


if(BestScoreInstances == 0)
	return(eHRnone);			// return no matches

int Idx;

if(pHits->Score > cMaxScore)
	pHits->Score = cMaxScore;

pCurHit = pHits;
for(Idx = 0; Idx < min(MaxHits,BestScoreInstances); Idx++,pCurHit++)
	{
	if((pEntry = MapChunkHit2Entry(pCurHit->Seg[0].MatchLoci))==NULL)
		return(eHRnone);

	pCurHit->Seg[0].ChromID = pEntry->EntryID;
	pCurHit->Seg[0].MatchLoci -= pEntry->StartOfs;
	if(pCurHit->Seg[1].MatchLoci > 0)
		{
		if((pEntry = MapChunkHit2Entry(pCurHit->Seg[1].MatchLoci))==NULL)
			return(eHRnone);
		pCurHit->Seg[1].ChromID = pEntry->EntryID;
		pCurHit->Seg[1].MatchLoci -= pEntry->StartOfs;
		}
	}
*pScore = pHits->Score;
*pLowHitInstances = BestScoreInstances;
*pLowMMCnt = pHits->Seg[0].Mismatches + pHits->Seg[1].Mismatches;
*pNxtLowMMCnt = *pLowMMCnt+2;

return(BestScoreInstances <= MaxHits ? eHRhits : eHRnone);
}

int						// < 0 if errors, eHRnone if none or too many matches or change, eHRhits if mumber of matches accepted
CSfxArray::LocateInDels(uint32_t ExtdProcFlags,			// flags indicating if lower levels need to do any form of extended processing with this specific read...
						 uint32_t ReadID,					// identifies this read
						 int microInDelLen,			// microInDel length maximum
						 int MaxTotMM,					// max number of mismatches allowed
						 int CoreLen,					// core window length
 					 	 eALStrand Align2Strand,		// align to this strand
  						 int *pLowHitInstances,			// In/Out number of match instances for lowest number of mismatches thus far for this read
						 int *pLowMMCnt,				// In/Out lowest number of mismatches thus far for this read
						 int *pNxtLowMMCnt,				// In/Out next to lowest number of mismatches thus far for this read
						 etSeqBase *pProbeSeq,			// probe
						 int ProbeLen,					// probe length
						 int MaxHits,					// (IN) process for at most this number of hits (currently this will always be 1)
						 tsHitLoci *pHits,				// where to return any hit
						 int *pScore,					// where to return InDel score
						 int CurMaxIter,				// max allowed iterations per subsegmented sequence when matching that subsegment
						 int NumAllocdIdentNodes,		// memory has been allocated by caller for holding upto this many tsIdentNodes
						 tsIdentNode *pAllocsIdentNodes) // memory allocated by caller for holding tsIdentNodes
{
int Phase;						// current processing phase: 0 - look for InDel anchored 5' to 3', 1 - look for InDel anchored 3' to 5'
int CurCoreSegOfs;				// current core segment relative start
int IterCnt;					// count iterator for current segment target matches
etSeqBase *pProbeBase;
etSeqBase *pTargBase;
char CurStrand;
int64_t TargIdx;
uint32_t NumCopies;

uint32_t NumTargSeqProc;

bool bFirstIter;				// set false after the first subsequence core returned by LocateFirstExact has been processed

int64_t TargSeqLeftIdx;
uint32_t ProbeSeqLeftIdx;
int TargMatchLen;

etSeqBase BisBase;

tsHitLoci RightHitLoci;
tsHitLoci LeftHitLoci;
int RsltRightInDel;
int RsltLeftInDel;
int BestScoreInstances;

int Cmp;
tsHitLoci *pCurHit;
tsSfxEntry *pEntry;
etSeqBase *pTarg;			// target sequence
void *pSfxArray;			// target sequence suffix array
int64_t SfxLen;				// number of suffixs in pSfxArray
BisBase=eBaseN;

// ensure suffix array block loaded for iteration!
if(m_pSfxBlock == NULL || m_pSfxBlock->ConcatSeqLen == 0)
	return(eBSFerrInternal);

pTarg = (etSeqBase *)&m_pSfxBlock->SeqSuffix[0];
pSfxArray = (void *)&m_pSfxBlock->SeqSuffix[m_pSfxBlock->ConcatSeqLen];
SfxLen = m_pSfxBlock->ConcatSeqLen;


pCurHit = NULL;

NumTargSeqProc = 0;
*pLowHitInstances = 0;
*pLowMMCnt = 0;
*pNxtLowMMCnt = 0;
memset(pHits,0,sizeof(tsHitLoci));
BestScoreInstances = 0;

if(Align2Strand == eALSCrick)
	{
	if(m_bColorspace)
		CSeqTrans::ReverseSeq(ProbeLen,pProbeSeq);
	else
		CSeqTrans::ReverseComplement(ProbeLen,pProbeSeq);
	CurStrand = '-';
	}
else
	CurStrand = '+';

do
	{
	Phase = 0;
	for(Phase = 0; Phase < 2; Phase+=1)
		{
		switch(Phase) {
			case 0:
				CurCoreSegOfs = 0;
				break;
			case 1:
				CurCoreSegOfs = ProbeLen - CoreLen;
				break;
			}

		if (m_pOccKMerClas != NULL && CoreLen == m_OccKMerLen && (OverOccKMerClas(CoreLen, &pProbeSeq[CurCoreSegOfs]) != 1))	// 1 if number of K-mers of CoreLen is within range 1..max
			continue;

		TargIdx = LocateFirstExact(&pProbeSeq[CurCoreSegOfs],CoreLen,pTarg,m_pSfxBlock->SfxElSize,pSfxArray,0,0,SfxLen-1);
		if(TargIdx == 0)        // 0 if no core segment matches
			continue;			// try for match on next core segment after shifting core to right
		TargIdx -= 1;
		IterCnt = 0;
		NumCopies = 0;
		bFirstIter = true;		// set false after the first subsequence core returned by LocateFirstExact has been processed
		while(!CurMaxIter || IterCnt < CurMaxIter)
			{
			if(!bFirstIter) {

				// ensure not about to iterate past end of suffix array!
				if((TargIdx + 1) >= (int64_t)m_pSfxBlock->ConcatSeqLen || (SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,TargIdx+1) + CoreLen) >  (int64_t)m_pSfxBlock->ConcatSeqLen)
					break;

				// check that this new putative core is still matching
				pTargBase = &pTarg[SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,TargIdx+1)];
				pProbeBase = &pProbeSeq[CurCoreSegOfs];
				if(m_bBisulfite)
					Cmp = BSCmpProbeTarg(pProbeBase,pTargBase,CoreLen);
				else
					{
					int Ofs;
					uint8_t El1;
					uint8_t El2;
					etSeqBase *pEl1= pProbeBase;
					etSeqBase *pEl2 = pTargBase;
					Cmp = 0;
					for(Ofs=0; Ofs < CoreLen; Ofs++)
						{
						El2 = *pEl2++ & 0x0f;
						if(El2 == eBaseEOS)
							{
							Cmp = -1;
							break;
							}
						El1 = *pEl1++ & 0x0f;
						if(El1 > El2)
							{
							Cmp = 1;
							break;
							}
						if(El1 < El2)
							{
							Cmp = -1;
							break;
							}
						}
					}

				if(Cmp != 0)				// will be non-zero if target no longer matches
					break;					// try next core segment
				TargIdx += 1;
				}

			bFirstIter = false;
			if(SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,TargIdx) < (uint32_t)CurCoreSegOfs)
				continue;

			TargMatchLen = ProbeLen;
			TargSeqLeftIdx = SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,TargIdx) - CurCoreSegOfs;
			ProbeSeqLeftIdx = 0;
			pEntry = MapChunkHit2Entry(SfxOfsToLoci(m_pSfxBlock->SfxElSize,pSfxArray,TargIdx));
			if(pEntry == NULL)
				continue;
			if(TargSeqLeftIdx < (int64_t)pEntry->StartOfs || ((TargSeqLeftIdx +  ProbeLen - 1) > (int64_t)pEntry->EndOfs))
				continue;

			// ensure comparisons are still within start/end range of target sequence/assembly
			if((TargSeqLeftIdx + ProbeLen) > (int64_t)m_pSfxBlock->ConcatSeqLen)
				continue;

			NumTargSeqProc += 1;
			IterCnt += 1;

			pTargBase = &pTarg[TargSeqLeftIdx];
			pProbeBase = pProbeSeq;

			// explore to the right if phase 0
			if(Phase == 0)
				{
				RsltRightInDel = ExploreInDelMatchRight(ExtdProcFlags,			// flags indicating if lower levels need to do any form of extended processing with this specific read...
					ReadID,								// identifies this read
					CurStrand,							// aligning on this strand
					microInDelLen,						// microInDel can be at most this length
					MaxTotMM,							// can be at most this many mismatches in total including any InDel
					ProbeLen,							// length of probe excluding any eBaseEOS
					pProbeBase,							// pts to probe sequence
					pEntry,								// target at TargSeqLeftIdx is in this sfx entry
					TargSeqLeftIdx,						// pTarg corresponds to this suffix index
 					pTargBase,							// pts to target sequence
					&RightHitLoci);						// where to return alignment loci

				if(RsltRightInDel > 0 &&  RightHitLoci.Score >= pHits->Score)
					{

					if(RightHitLoci.Score == pHits->Score)
						{
						if(pHits->Seg[0].MatchLoci == RightHitLoci.Seg[0].MatchLoci)
							continue;
						if(++BestScoreInstances > MaxHits)
							continue;
						}
					else
						BestScoreInstances = 0;
					pHits[BestScoreInstances++] = RightHitLoci;
					}
				}

			// explore to the left if phase 1
			if(Phase == 1)
				{
				RsltLeftInDel = ExploreInDelMatchLeft(ExtdProcFlags,			// flags indicating if lower levels need to do any form of extended processing with this specific read...
						 ReadID,					// identifies this read
					CurStrand,							// aligning on this strand
					microInDelLen,			// microInDel can be at most this length
					MaxTotMM,				// can be at most this many mismatches in total
					ProbeLen,				// length of probe excluding any eBaseEOS
					pProbeBase,				// pts to probe sequence
					pEntry,					// target at TargSeqLeftIdx is in this sfx entry
					TargSeqLeftIdx,			// pTargBase corresponds to this suffix index
					pTargBase,				// pts to target sequence
					&LeftHitLoci);			// where to return alignment loci


				if(RsltLeftInDel > 0 &&  LeftHitLoci.Score >= pHits->Score)
					{
					if(LeftHitLoci.Score == pHits->Score)
						{
						if(pHits->Seg[0].MatchLoci == LeftHitLoci.Seg[0].MatchLoci)
							continue;
						if(++BestScoreInstances > MaxHits)
							continue;

						}
					else
						BestScoreInstances = 0;
					pHits[BestScoreInstances++] = LeftHitLoci;
					}
				}
			}
		if(BestScoreInstances >= 1 && pHits->Score >= cMaxScore)							// unlikely to do any better than this so early exit...
			{
			Align2Strand = eALSnone;
			break;
			}
		}


	if(CurStrand == '+' && Align2Strand == eALSboth)	// if just processed watson '+' strand then will need to process crick or '-' strand if processing both strands
		{
		if(m_bColorspace)
			CSeqTrans::ReverseSeq(ProbeLen,pProbeSeq);
		else
			CSeqTrans::ReverseComplement(ProbeLen,pProbeSeq);
		CurStrand = '-';
		Align2Strand = eALSCrick;
		NumTargSeqProc = 0;
		}
	else									// either processed both or just the crick strand so no further processing required
		Align2Strand = eALSnone;			// two passes max - first on the '+' strand with optional 2nd on '-' strand
	}
while(!(BestScoreInstances >= 1 && pHits->Score >= cMaxScore) && Align2Strand != eALSnone);

if(CurStrand == '-')									// restore probe sequence if had started processing '-' strand
		if(m_bColorspace)
			CSeqTrans::ReverseSeq(ProbeLen,pProbeSeq);
		else
			CSeqTrans::ReverseComplement(ProbeLen,pProbeSeq);

if(BestScoreInstances == 0)
	return(eHRnone);			// return no matches

int Idx;
if(pHits->Score > cMaxScore)
	pHits->Score = cMaxScore;

pCurHit = pHits;
tsSfxEntry *pEntry2;

for(Idx = 0; Idx < min(MaxHits,BestScoreInstances); Idx++,pCurHit++)
	{
	if(pCurHit->Seg[0].Strand != '-' &&  pCurHit->Seg[0].Strand != '+')
		{
		pCurHit->Seg[0].Strand = '?';
		}
	if((pEntry = MapChunkHit2Entry(pCurHit->Seg[0].MatchLoci))==NULL)
		return(eHRnone);
	if((pEntry2 = MapChunkHit2Entry(pCurHit->Seg[1].MatchLoci))==NULL)
		return(eHRnone);
	if(pEntry->EntryID != pEntry2->EntryID)			// InDels must only be accepted if intra-chromosome; Not sure as to how but occasionally ExploreInDelMatchLeft/Right
		return(eHRnone);							// are identifying intra-chrom InDels in the Wheat GSS assembly, not observed in any other assembly
													// This should be investigated so this check is a short term work around...

	pCurHit->Seg[0].ChromID = pEntry->EntryID;
	pCurHit->Seg[0].MatchLoci -= pEntry->StartOfs;

	if(pCurHit->Seg[1].MatchLoci > 0)
		{
		pCurHit->Seg[1].ChromID = pEntry2->EntryID;
		pCurHit->Seg[1].MatchLoci -= pEntry2->StartOfs;
		}
	}
*pScore = pHits->Score;
*pLowHitInstances = min(MaxHits,BestScoreInstances);
*pLowMMCnt = pHits->Seg[0].Mismatches + pHits->Seg[1].Mismatches;
*pNxtLowMMCnt = *pLowMMCnt+2;

return(BestScoreInstances <= MaxHits ? eHRhits : eHRnone);
}

//
// AlignReads
// Multiphase align reads
int														// < 0 if errors, 0 if no matches or change, 1 if mumber of matches accepted, 2 MMDelta criteria not met, 3 too many match instances
CSfxArray::AlignReads(uint32_t ExtdProcFlags,			// flags indicating if lower levels need to do any form of extended processing with this specific read...
						 uint32_t ReadID,					// identifies this read
						 int MinChimericLen,			// minimum chimeric length as a percentage (0 to disable, otherwise 15..99) of probe sequence
						 int MaxTotMM,					// max number of mismatches allowed
						 int CoreLen,					// core window length
						 int CoreDelta,					// core window offset increment (1..n)
						 int MaxNumCoreSlides,			// limit number of core slides to this many per strand over the read length
						 int MinCoreLen,				// do not reduce core below this length
						 int MMDelta,					// minimum (1..n) mismatch difference between the best and next best core alignment
 					 	 eALStrand Align2Strand,		// align to this strand
						 int microInDelLen,				// microInDel length maximum
						 int MaxSpliceJunctLen,			// junction has to be no longer than this length
  						 int *pLowHitInstances,			// In/Out number of match instances for lowest number of mismatches thus far for this read
						 int *pLowMMCnt,				// In/Out lowest number of mismatches thus far for this read
						 int *pNxtLowMMCnt,				// In/Out next to lowest number of mismatches thus far for this read
						 etSeqBase *pProbeSeq,			// probe
						 int ProbeLen,					// probe length
						 int MaxHits,					// (IN) process for at most this number of hits
						 tsHitLoci *pHits,				// where to return hits (at most MaxHits)
 						 int NumAllocdIdentNodes,		// memory has been allocated by caller for holding upto this many tsIdentNodes
						 tsIdentNode *pAllocsIdentNodes) // memory allocated by caller for holding tsIdentNodes

{
int Rslt;
int AllowMM;
int CL;
int BestInDelScore;
int SpliceCore;

if(MaxTotMM > 0)
	{
	for(AllowMM = 0; AllowMM <= MaxTotMM; AllowMM++)
		{
		CL = ProbeLen / (AllowMM+MMDelta);
		if(CL <= CoreLen)
			break;
		Rslt = LocateCoreMultiples(ExtdProcFlags,ReadID,0,AllowMM,CL,CL,MaxNumCoreSlides,MMDelta,Align2Strand,pLowHitInstances,
								pLowMMCnt,pNxtLowMMCnt,pProbeSeq,ProbeLen,MaxHits,pHits,m_MaxIter,NumAllocdIdentNodes,pAllocsIdentNodes);
		if(Rslt != 0)
			return(Rslt);
		if(AllowMM > MaxTotMM)
			break;
		}
	}
else
	AllowMM = 0;

if(AllowMM <= MaxTotMM)
	{
	Rslt = LocateCoreMultiples(ExtdProcFlags,ReadID,0,MaxTotMM,CoreLen,CoreDelta,MaxNumCoreSlides,MMDelta,Align2Strand,pLowHitInstances,
								pLowMMCnt,pNxtLowMMCnt,pProbeSeq,ProbeLen,MaxHits,pHits,m_MaxIter,NumAllocdIdentNodes,pAllocsIdentNodes);
	if(Rslt != 0)
		return(Rslt);
	}

// Note: only interested in unique splice junctions or InDels
if(Rslt == 0 && microInDelLen > 0)
	{
	// NOTE:
	//    This is a rather compute intensive function, in order to get some reasonable throughput-
	//	  a) the CoreLen is increased to be a min of either (ProbeLen-1)/2 or (CoreLen * 2) so as to reduce the number of internal processing iterations
	//    b) the MaxTotMM is clamped to be at most cMaxJunctAlignMM
// changed - smaller cores are likely to have more hits
	SpliceCore = CoreLen;
	Rslt = LocateInDels(ExtdProcFlags,ReadID,microInDelLen,MaxTotMM > cMaxMicroInDelMM ? cMaxMicroInDelMM : MaxTotMM,SpliceCore,Align2Strand,pLowHitInstances,
								pLowMMCnt,pNxtLowMMCnt,pProbeSeq,ProbeLen,1,pHits,&BestInDelScore,m_MaxIter,NumAllocdIdentNodes,pAllocsIdentNodes);
	if(Rslt != 0)
		return(Rslt);
	}

if(Rslt == 0 && MaxSpliceJunctLen > 0)
	// NOTE:
	//    This is a rather compute intensive function, in order to get some reasonable throughput-
	//	  a) the CoreLen is increased to be a min of either (ProbeLen-1)/2 or (CoreLen * 2) so as to reduce the number of internal processing iterations
	//    b) the MaxTotMM is clamped to be at most cMaxJunctAlignMM
	{
// changed - smaller cores are likely to have more hits, was: SpliceCore = min((CoreLen * 2), (ProbeLen - 1) / 2);
	SpliceCore = CoreLen;
	Rslt=LocateSpliceJuncts(ExtdProcFlags,ReadID,MaxSpliceJunctLen,MaxTotMM > cMaxJunctAlignMM ? cMaxJunctAlignMM : MaxTotMM,
						SpliceCore,Align2Strand,pLowHitInstances,
								pLowMMCnt,pNxtLowMMCnt,pProbeSeq,ProbeLen,1,pHits,&BestInDelScore,m_MaxIter,NumAllocdIdentNodes,pAllocsIdentNodes);
	if(Rslt != 0)
		return(Rslt);
	}

if(MinChimericLen > 0)
	{
	CoreLen = max(MinCoreLen,ProbeLen/(MaxTotMM+4));
	CoreDelta = max(ProbeLen/(MaxNumCoreSlides-1),CoreLen);
	Rslt = LocateCoreMultiples(ExtdProcFlags,ReadID,MinChimericLen,MaxTotMM,CoreLen,CoreDelta,MaxNumCoreSlides,MMDelta,Align2Strand,pLowHitInstances,
								pLowMMCnt,pNxtLowMMCnt,pProbeSeq,ProbeLen,MaxHits,pHits,m_MaxIter,NumAllocdIdentNodes,pAllocsIdentNodes);
	return(Rslt);
	}

return(0);
}



int64_t			// index+1 in pSfxArray of first exactly matching probe or 0 if no match
CSfxArray::LocateFirstExact(etSeqBase *pProbe,  // pts to probe sequence
				  int ProbeLen,					// probe length to exactly match over
				  etSeqBase *pTarg,				// target sequence
				  int SfxElSize,				// size in bytes of suffix element - expected to be either 4 or 5
				  void *pSfxArray,				// target sequence suffix array
				  int64_t TargStart,				// position in pTarg (0..n) corresponding to start of suffix array
				  int64_t SfxLo,					// low index in pSfxArray
				  int64_t SfxHi,					// high index in pSfxArray
				  bool bIgnoreKMerCores)	     // true to bypass checking for CoreKMers starts
{
etSeqBase *pEl1;
etSeqBase *pEl2;
uint8_t El1;
uint8_t El2;

int CmpRslt;
int Ofs;
int64_t Mark;
int64_t TargPsn;

if(!bIgnoreKMerCores && ProbeLen == m_CoreKMerLen && m_pCoreKMers != NULL)
	{
	TargPsn = 0;
	NumExactKMers(0,ProbeLen, pProbe,&TargPsn);
	return(TargPsn);
	}

do {
	pEl1 = pProbe;
	TargPsn = ((int64_t)SfxLo + SfxHi) / 2L;
	pEl2 = &pTarg[SfxOfsToLoci(SfxElSize,pSfxArray,TargPsn + TargStart)];
	if(m_bBisulfite)
		CmpRslt = BSCmpProbeTarg(pEl1,pEl2,ProbeLen);
	else
		{
		CmpRslt = 0;
		for(Ofs=0; Ofs < ProbeLen; Ofs++)
			{
			El2 = *pEl2++ & 0x0f;
			if(El2 == eBaseEOS)
				{
				CmpRslt = -1;
				break;
				}
			El1 = *pEl1++ & 0x0f;
			if(El1 > El2)
				{
				CmpRslt = 1;
				break;
				}
			if(El1 < El2)
				{
				CmpRslt = -1;
				break;
				}
			}
		}

	if(!CmpRslt)	// if a match then may not be the lowest indexed match
		{
		if(TargPsn == 0 || SfxLo == TargPsn) // check if already lowest
			return(TargPsn + 1);
		// iterate until lowest located
		while(1) {
			if(CmpRslt == 0)
				{
				Mark = TargPsn;
				if(Mark == 0)
					return(Mark+1);
				SfxHi = TargPsn - 1;
				}
			TargPsn = ((int64_t)SfxLo + SfxHi) / 2L;
			pEl2 = &pTarg[SfxOfsToLoci(SfxElSize,pSfxArray,TargPsn + TargStart)];
			if(m_bBisulfite)
				CmpRslt = BSCmpProbeTarg(pEl1,pEl2,ProbeLen);
			else
				{
				pEl1 = pProbe;
				CmpRslt = 0;
				for(Ofs=0; Ofs < ProbeLen; Ofs++)
					{
					El2 = *pEl2++ & 0x0f;
					if(El2 == eBaseEOS)
						{
						CmpRslt = -1;
						break;
						}
					El1 = *pEl1++ & 0x0f;
					if(El1 > El2)
						{
						CmpRslt = 1;
						break;
						}
					if(El1 < El2)
						{
						CmpRslt = -1;
						break;
						}
					}
				}
			if(CmpRslt == 0)				// 0 if still matching
				continue;
			SfxLo = TargPsn + 1;
			if(SfxLo == Mark)
				return(Mark+1);
			}
		}

	if(CmpRslt < 0)
		{
		if(TargPsn == 0)
			break;
		SfxHi = TargPsn - 1;
		}
	else
		SfxLo = TargPsn+1;
	}
while(SfxHi >= SfxLo);

return(0);	// unable to locate any instance of pProbe
}

// iterate over all KMers of KMerLen in suffix array and update each KMer record
int
CSfxArray::InitialiseCoreKMers(int KMerLen)			// initialise for KMer cores of this length
{
etSeqBase *pCurKMer;
void *pSfxArray;			// target sequence suffix array
int64_t SfxLen;				// number of suffixs in pSfxArray

// ensure suffix array block loaded for iteration!
if (m_pSfxBlock == NULL || m_pSfxBlock->ConcatSeqLen == 0)
	return(eBSFerrInternal);

pCurKMer = (etSeqBase *)&m_pSfxBlock->SeqSuffix[0];
pSfxArray = (void *)&m_pSfxBlock->SeqSuffix[m_pSfxBlock->ConcatSeqLen];
SfxLen = m_pSfxBlock->ConcatSeqLen;

// limited to KMer cores in the range of 1 through to cMaxKmerLen inclusive
if(KMerLen < 1 || KMerLen > cMaxKmerLen)
	return(eBSFerrParams);
SerialiseBaseFlags();
if (!(KMerLen == m_CoreKMerLen && m_pCoreKMers != NULL))
	{
	m_CoreKMerLen = 0;
	if (m_pCoreKMers != NULL)
		{
#ifdef _WIN32
		free(m_pCoreKMers);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
		if (m_pCoreKMers != MAP_FAILED)
			munmap(m_pCoreKMers, m_AllocdCoreKMersMem);
#endif
		m_pCoreKMers = NULL;
		}
	m_NumCoreKMers = (size_t)0x04 << ((KMerLen - 1) * 2);
	m_AllocdCoreKMersMem = m_NumCoreKMers * sizeof(tsCoreKMer);
#ifdef _WIN32
	m_pCoreKMers = (tsCoreKMer *)malloc(m_AllocdCoreKMersMem);
	if (m_pCoreKMers == NULL)
		{
		AddErrMsg("CSfxArray::InitOverOccKMers", "unable to allocate %zu bytes for core KMers", m_AllocdCoreKMersMem);
		ReleaseBaseFlags();
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and seems to have issues if more than 2GB allocation
	m_pCoreKMers = (tsCoreKMer *)mmap(NULL, m_AllocdCoreKMersMem, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	if (m_pCoreKMers == MAP_FAILED)
		{
		AddErrMsg("CSfxArray::InitOverOccKMers", "unable to allocate %zu bytes for core KMers", m_AllocdCoreKMersMem);
		m_pCoreKMers = NULL;
		ReleaseBaseFlags();
		return(eBSFerrMem);
		}
#endif
	}
memset(m_pCoreKMers, 0, m_AllocdCoreKMersMem);
m_CoreKMerLen = KMerLen;
ReleaseBaseFlags();

return(eBSFSuccess);
}

uint32_t			// number of exactly matching KMers up to specified MaxToMatch inclusive, or 0 if no match, MaxToMatch+1 if more than MaxToMatch Kmers
CSfxArray::NumExactKMers(uint32_t MaxToMatch,	// if 0 then no limit on matches, otherwise match up to this number of instances and if more then return MaxToMatch+1 
					int KMerLen,				// KMer length to exactly match over
					etSeqBase *pKMerSeq,		// pts to KMer sequence
					int64_t *pFirstTargIdx)		// optionally pts to returned first target index
{
uint32_t NumCopies;
int Idx;
int64_t KMerCoreIdx;
int SfxElSize;
etSeqBase *pTarg;
void *pSfxArray;
int64_t SfxLen;
int64_t FirstTargIdx;
int64_t LastTargIdx;
int64_t LimitTargIdx;

if(pFirstTargIdx != NULL)
	*pFirstTargIdx = 0;

// pKMerSeq may already have been processed for exact matches so can quickly lookup the number of exact matches
if (KMerLen == m_CoreKMerLen && m_pCoreKMers != NULL)
	{
	// generate index from bases in pKMerSeq
	etSeqBase Base;
	pTarg = pKMerSeq;
	KMerCoreIdx = 0;
	for(Idx = 0; Idx < KMerLen; Idx++)
		{
		Base = *pTarg++ & 0x0f;
		if(Base > eBaseT)
			{
			KMerCoreIdx = -1;
			break;
			}
		KMerCoreIdx |= (((int64_t)Base & 0x03) << (Idx * 2));
		}
	if (KMerCoreIdx != -1)
		{
		SerialiseBaseFlags();
		if (m_pCoreKMers[KMerCoreIdx].FirstCore != 0)		// will be -1 if was previously processed for first and there was no core KMer in suffix array
			{
			NumCopies = m_pCoreKMers[KMerCoreIdx].NumCores;
			if (pFirstTargIdx != NULL)
				*pFirstTargIdx = m_pCoreKMers[KMerCoreIdx].FirstCore;
			ReleaseBaseFlags();
			if (MaxToMatch > 0 && NumCopies > MaxToMatch)
				return(MaxToMatch + 1);
			return(NumCopies);
			}
		ReleaseBaseFlags();
		}
	}
else
	KMerCoreIdx = -1;

SfxElSize = m_pSfxBlock->SfxElSize;
pTarg = (etSeqBase *)&m_pSfxBlock->SeqSuffix[0];
pSfxArray = (void *)&m_pSfxBlock->SeqSuffix[m_pSfxBlock->ConcatSeqLen];
SfxLen = m_pSfxBlock->ConcatSeqLen;
FirstTargIdx = LocateFirstExact(pKMerSeq,KMerLen,pTarg,m_pSfxBlock->SfxElSize,pSfxArray,0,0,SfxLen-1,true);
if(FirstTargIdx < 1)		// not even one match?
	{
	if(KMerCoreIdx != -1)
		{
		SerialiseBaseFlags();
		m_pCoreKMers[KMerCoreIdx].FirstCore = -1;
		m_pCoreKMers[KMerCoreIdx].NumCores = 0;
		ReleaseBaseFlags();
		}
	return(0); // not a single match!
	}

if (pFirstTargIdx != NULL)
	*pFirstTargIdx = FirstTargIdx; // NOTE: LocateFirstExact() returns FirstTargIdx + 1 if able to locate the first exactly matching, otherwise 0 was returned!!!!!

// have at least one match
if(KMerCoreIdx != -1) 
	{
	// exactly how many matches are there?
	// NOTE: clamps number of matches to be no more than 0x0fffffff0
	LastTargIdx = LocateLastExact(pKMerSeq, KMerLen, pTarg, m_pSfxBlock->SfxElSize, pSfxArray, 0, FirstTargIdx -1, SfxLen - 1,true); // LocateLastExact() returns LastTargIdx + 1 if able to locate the last exactly matching, otherwise 0 was returned 
	SerialiseBaseFlags();
	m_pCoreKMers[KMerCoreIdx].FirstCore = FirstTargIdx;
	NumCopies = m_pCoreKMers[KMerCoreIdx].NumCores = LastTargIdx > 0 ? (uint32_t)min((uint64_t)(1 + LastTargIdx - FirstTargIdx),(uint64_t)0x0fffffff0) : 0;
	m_pCoreKMers[KMerCoreIdx].NumCores = NumCopies;
	ReleaseBaseFlags();
	if (MaxToMatch > 0 && NumCopies > MaxToMatch)
		return(MaxToMatch + 1);
	return(NumCopies);
	}

if(MaxToMatch > 0)
	LimitTargIdx = min(SfxLen-1,(int64_t)((uint64_t)FirstTargIdx + MaxToMatch + 1));
else
	LimitTargIdx = SfxLen-1;
LastTargIdx = LocateLastExact(pKMerSeq,KMerLen,pTarg,m_pSfxBlock->SfxElSize,pSfxArray,0,FirstTargIdx - 1,LimitTargIdx,true);
NumCopies = LastTargIdx > 0 ? (uint32_t)(1 + LastTargIdx - FirstTargIdx) : 0;

if(MaxToMatch > 0 && NumCopies > MaxToMatch)
	NumCopies = MaxToMatch+1;
return(NumCopies);
}

int64_t			// index+1 in pSfxArray of last exactly matching probe or 0 if no match
CSfxArray::LocateLastExact(etSeqBase *pProbe, // pts to probe sequence
				  int ProbeLen,					// probe length to exactly match over
				  etSeqBase *pTarg,				// target sequence
				  int SfxElSize,				// size in bytes of suffix element - expected to be either 4 or 5
				  void *pSfxArray,				// target sequence suffix array
				  int64_t TargStart,				// position in pTarg (0..n) corresponding to start of suffix array
				  int64_t SfxLo,					// low index in pSfxArray
				  int64_t SfxHi,					// high index in pSfxArray
				  bool bIgnoreKMerCores)	     // true to bypass checking for CoreKMers starts
{
etSeqBase *pEl1;
etSeqBase *pEl2;
uint8_t El1;
uint8_t El2;

int CmpRslt;
int Ofs;
int64_t Mark;
int64_t TargPsn;
int64_t SfxHiMax = SfxHi;
int64_t SfxLoMax = SfxLo;

if (!bIgnoreKMerCores && ProbeLen == m_CoreKMerLen && m_pCoreKMers != NULL)
	{
	Mark = NumExactKMers(0, ProbeLen, pProbe, &TargPsn);
	if(Mark > 0)
		return(TargPsn + Mark - 1);
	else
		return(0);
	}

do {
	pEl1 = pProbe;
	TargPsn = ((int64_t)SfxLo + SfxHi) / 2L;
	pEl2 = &pTarg[SfxOfsToLoci(SfxElSize,pSfxArray,TargPsn + TargStart)];
	if(m_bBisulfite)
		CmpRslt = BSCmpProbeTarg(pEl1,pEl2,ProbeLen);
	else
		{
		CmpRslt = 0;
		for(Ofs=0; Ofs < ProbeLen; Ofs++)
			{
			El2 = *pEl2++ & 0x0f;
			if(El2 == eBaseEOS)
				{
				CmpRslt = -1;
				break;
				}
			El1 = *pEl1++ & 0x0f;
			if(El1 > El2)
				{
				CmpRslt = 1;
				break;
				}
			if(El1 < El2)
				{
				CmpRslt = -1;
				break;
				}
			}
		}

	if(!CmpRslt)	// if a match then may not be the highest indexed match
		{
		if(TargPsn == SfxHiMax || SfxHi == TargPsn) // check if already highest
			return(TargPsn + 1);
		// iterate until highest located
		while(1) {
			if(CmpRslt == 0)
				{
				Mark = TargPsn;
				if(Mark == SfxHi)
					return(Mark+1);
				SfxLo = TargPsn + 1;
				}
			TargPsn = ((int64_t)SfxLo + SfxHi) / 2L;
			pEl2 = &pTarg[SfxOfsToLoci(SfxElSize,pSfxArray,TargPsn + TargStart)];
			if(m_bBisulfite)
				CmpRslt = BSCmpProbeTarg(pEl1,pEl2,ProbeLen);
			else
				{
				pEl1 = pProbe;
				CmpRslt = 0;
				for(Ofs=0; Ofs < ProbeLen; Ofs++)
					{
					El2 = *pEl2++ & 0x0f;
					if(El2 == eBaseEOS)
						{
						CmpRslt = -1;
						break;
						}
					El1 = *pEl1++ & 0x0f;
					if(El1 > El2)
						{
						CmpRslt = 1;
						break;
						}
					if(El1 < El2)
						{
						CmpRslt = -1;
						break;
						}
					}
				}
			if(CmpRslt == 0)				// 0 if still matching
				continue;
			SfxHi = TargPsn - 1;
			if(SfxHi == Mark)
				return(Mark+1);
			}
		}

	if(CmpRslt < 0)
		{
		if(TargPsn == 0)
			break;
		SfxHi = TargPsn - 1;
		}
	else
		SfxLo = TargPsn+1;
	}
while(SfxHi >= SfxLo);

return(0);	// unable to locate any instance of pProbe
}


int CSfxArray::CompareExtdMatches( const void *arg1, const void *arg2 )
{
tsSfxMatch *pEl1 = (tsSfxMatch *)arg1;
tsSfxMatch *pEl2 = (tsSfxMatch *)arg2;
if(pEl1->ProbePsn == pEl2->ProbePsn )
	{
	if(pEl1->TargPsn < pEl2->TargPsn)
		return(-1);
	if(pEl1->TargPsn > pEl2->TargPsn)
		return(1);
	if(pEl1->MatchLen < pEl2->MatchLen)
		return(-1);
	if(pEl1->MatchLen > pEl2->MatchLen)
		return(1);
	return(0);
	}

return(pEl1->ProbePsn < pEl2->ProbePsn ? -1 : 1);
}


// NEW CODE FOLLOWS


CAlignOpCode::CAlignOpCode(void)
{
m_MaxAlignOpCodes = cMaxAlignOpcodesLen;
m_CurNumAlignOpCodes = 0;
}

CAlignOpCode::~CAlignOpCode()
{
}

void
CAlignOpCode::Reset(void)
{
m_CurNumAlignOpCodes = 0;
}

void
CAlignOpCode::Clone(CAlignOpCode &CloneFrom)		// clone from
{
m_CurNumAlignOpCodes = CloneFrom.m_CurNumAlignOpCodes;				// current number of opcodes + associated values in m_pAlignOpCodes
if(m_CurNumAlignOpCodes)
	memmove(m_AlignOpCodes,CloneFrom.m_AlignOpCodes,m_CurNumAlignOpCodes);	// to hold alignment opcodes + associated values for lowest scored alignment
}

bool
CAlignOpCode::Init(int Len,uint8_t *pOpCodes)
{
if(Len > cMaxAlignOpcodesLen)
	return(false);
memmove(m_AlignOpCodes,pOpCodes,Len);
m_CurNumAlignOpCodes = Len;
return(true);
}

bool
CAlignOpCode::Append(int Len,uint8_t *pOpCodes)
{
if((m_CurNumAlignOpCodes + Len) > cMaxAlignOpcodesLen)
	return(false);
memmove(&m_AlignOpCodes[m_CurNumAlignOpCodes],pOpCodes,Len);
m_CurNumAlignOpCodes += Len;
return(true);
}

bool
CAlignOpCode::Append(CAlignOpCode &AppendFrom)
{
if((m_CurNumAlignOpCodes + AppendFrom.m_CurNumAlignOpCodes) > cMaxAlignOpcodesLen)
	return(false);
memmove(&m_AlignOpCodes[m_CurNumAlignOpCodes],AppendFrom.m_AlignOpCodes,AppendFrom.m_CurNumAlignOpCodes);
m_CurNumAlignOpCodes += AppendFrom.m_CurNumAlignOpCodes;
return(true);
}

// truncate at the Nth mismatch
bool
CAlignOpCode::TruncNthMM(int NthMM)
{
uint8_t *pAlignOpCodes;
int Ofs;
if(m_CurNumAlignOpCodes == 0 || !NthMM)
	return(false);
pAlignOpCodes = m_AlignOpCodes;
for(Ofs = 0;Ofs < m_CurNumAlignOpCodes; Ofs++,pAlignOpCodes++)
	{
	if(*pAlignOpCodes != eAOPMismatch)
		continue;
	if(--NthMM)
		{
		m_CurNumAlignOpCodes = Ofs;
		return(true);
		}
	}
return(false);
}

// AppendNthMM
// Locate the Nth instance of a mismatch and then replace following opcodes
bool
CAlignOpCode::AppendNthMM(int NthMM,CAlignOpCode &AppendFrom)
{
uint8_t *pAlignOpCodes;
int Ofs;
if(m_CurNumAlignOpCodes == 0 || !NthMM)
	return(false);
pAlignOpCodes = m_AlignOpCodes;
for(Ofs = 0;Ofs < m_CurNumAlignOpCodes; Ofs++,pAlignOpCodes++)
	{
	if(*pAlignOpCodes != eAOPMismatch)
		continue;
	if(--NthMM)
		break;
	}
if(Ofs == m_CurNumAlignOpCodes)
	return(false);
if((Ofs + AppendFrom.m_CurNumAlignOpCodes) > m_MaxAlignOpCodes)
	return(false);
m_CurNumAlignOpCodes = Ofs + AppendFrom.m_CurNumAlignOpCodes;
if(AppendFrom.m_CurNumAlignOpCodes)
	memmove(pAlignOpCodes,AppendFrom.m_AlignOpCodes,AppendFrom.m_CurNumAlignOpCodes);
return(false);
}

eAlignOPcodes										// opcode popped, or eAOPempty if opcodes empty
CAlignOpCode::PopAlignOpCode(uint32_t *pCounts)		// opcode's associated count or value
{
uint8_t *pAlignOpCodes;
uint8_t Val;
uint32_t Counts;
*pCounts = 0;
if(!m_CurNumAlignOpCodes)
	return(eAOPempty);
pAlignOpCodes = &m_AlignOpCodes[m_CurNumAlignOpCodes-1];
Counts = 0;
while(m_CurNumAlignOpCodes--)
	{
	Val = *pAlignOpCodes--;
	if(Val & 0x080)					// all opcodes have bit 7 set, values have bit 7 reset
		break;
	Counts <<= 7;
	Counts |= Val;
	}
*pCounts = Counts;
return((eAlignOPcodes)Val);
}


bool										// false if no room to push requested opcode, otherwise true
CAlignOpCode::PushAlignOpCode(eAlignOPcodes OpCode,	// opcode to push
				uint32_t Counts)			// opcode associated count or value
{
uint8_t *pAlignOpCodes;
uint32_t Len;
int RemainDepth;
RemainDepth = m_MaxAlignOpCodes - m_CurNumAlignOpCodes;
if(!RemainDepth)
	return(false);
pAlignOpCodes = &m_AlignOpCodes[m_CurNumAlignOpCodes];
switch(OpCode) {
	case eAOPMismatch:				// probe base does not match target sequence base, probe base is in next byte bits 0..3 and targ sequence base is in bits 4..7
		if(RemainDepth < 2)			// need at least 2
			return(false);
		*pAlignOpCodes++ = (uint8_t)OpCode;
		*pAlignOpCodes++ = (uint8_t)(Counts & 0x07f);
		m_CurNumAlignOpCodes += 2;
		break;
	case eAOPMatch:					// probe bases exactly matches target sequence bases for length specified in following byte(s)
	case eAOPInsert:				// probe has inserted bases relative to target, length of insertion is specified in following byte(s)
	case eAOPDelete:				// probe has deleted bases relative to target, length of insertion is specified in following byte(s)
	case eAOPSkip:					// skip this number of bases in target, number of bases to skip is in specified in following byte(s)
		Len = Counts;				// check that there is sufficent space to hold opcode + value
		while(Len > 0)
			{
			RemainDepth -= 1;
			Len >>= 7;
			}
		if(RemainDepth < 1)			// need at least another 1 for opcode
			return(false);
		*pAlignOpCodes++ = (uint8_t)OpCode;
		m_CurNumAlignOpCodes += 1;
		while(Counts > 0)
			{
			*pAlignOpCodes++ = (uint8_t)(Counts & 0x07f);
			Counts >>= 7;
			m_CurNumAlignOpCodes += 1;
			}
		break;
		}
return(true);
}

bool										// success or otherwise
CAlignOpCode::PushAlignMismatch(etSeqBase ProbeBase,	// probe base
				  etSeqBase RefBase)    // target or reference base
{
return(PushAlignOpCode(eAOPMismatch,((RefBase & 0x07) << 3) | (ProbeBase & 0x07)));
}

// AlignPartnerRead
// Have been able to unquely align one read out of a pair, now need to align the other read
// if not bPairStrand
//		if PE1 was to sense strand then expect PE2 on the antisense strand downstream towards the 3' end of targeted chrom:		b3primeExtend=true,bAntisense=true
//		if PE1 was to antisense strand then expect PE2 on the sense strand upstream towards the 5' end of targeted chrom:		b3primeExtend=false,bAntisense=false
//		if PE2 was to sense strand then expect PE1 on the antisense strand downstream towards the 3' end of targeted chrom:		b3primeExtend=true,bAntisense=true
//		if PE2 was to antisense strand then expect PE1 on the sense strand upstream towards the 5' end of targeted chrom:		b3primeExtend=false,bAntisense=false

// if bPairStrand
//		if PE1 was to sense strand then expect PE2 on the sense strand downstream towards the 3' end of targeted chrom:			b3primeExtend=true,bAntisense=false
//		if PE1 was to antisense strand then expect PE2 on the antisense strand upstream towards the 5' end of targeted chrom:	b3primeExtend=false,bAntisense=true
//		if PE2 was to sense strand then expect PE1 on the sense strand upstream towards the 5' end of targeted chrom:			b3primeExtend=false,bAntisense=false
//		if PE2 was to antisense strand then expect PE2 on the antisense strand downstream towards the 3' end of targeted chrom:	b3primeExtend=true,bAntisense=true

int									// -1 if errors, 0 if no match, 1 if mismatches only, 2 InDel and insertion into probe, 3 InDel and deletion from probe
CSfxArray::AlignPairedRead(bool b3primeExtend,	// if false extend towards 5' of targeted chrom, if true extend towards 3' of targeted chrom
			bool bAntisense,		// if false extend with read sense, if true extend with read RevCpl'd (antisense)
			uint32_t ChromID,		// accepted aligned read was on this chromosome
			uint32_t StartLoci,		// accepted aligned read started at this loci
			uint32_t EndLoci,		// and ending at this loci
			int MinInsertSize,	  // expecting partner to align at least this distance away from accepted aligned read (inclusive of read lengths)
			int MaxInsertSize,	  // but no more than this distance away (inclusive of read length)
			int MaxAllowedMM,	  // any accepted alignment can have at most this many mismatches
			int MinHamming,		  // and must be at least this Hamming away from the next best putative alignment
			int ReadLen,		  // length of read excluding any eBaseEOS
			int MinChimericLen,		// minimum chimeric length as a percentage (0 to disable, otherwise 15..99) of probe sequence
			int CoreLen,			// core window length, 0 to disable
			int CoreDelta,			// core window offset increment (1..n)
			int MaxNumCoreSlides,	// max number of times to slide core
		    etSeqBase *pRead,	  // pts to 5' start of read sequence
		    tsHitLoci *pAlign)	  // where to return any paired read alignment loci
{
int LowestNumMM;
int LowestMMcnt;
int AlignLoci;
uint32_t PutChromLen;
etSeqBase *pPutChromSeq;
uint32_t PutHitLoci;
uint32_t PutCoreOfs;
int64_t PutHitIdx;
int64_t NxtPutHitIdx;
uint32_t StartPutTargLoci;
uint32_t EndPutTargLoci;

int32_t MinPutLen;

uint32_t ChromSeqLen;
etSeqBase ReadSeq[cMaxReadLen+1];

if(pAlign != nullptr)
	memset(pAlign,0,sizeof(tsHitLoci));

// double check that the Start/End loci are within the chromosome length
if((ChromSeqLen = GetSeqLen(ChromID)) == 0)
	return(-1);
if(StartLoci >= EndLoci ||
	EndLoci >= ChromSeqLen)
	return(-1);

// calc actual minimum alignment length that would be accepted after any 5/3' chimeric trimming
if(CoreLen > 0 && MinChimericLen >= 15 && MinChimericLen <= 99) // if MinChimericLen > 0 then has to be in range 15% .. 99% otherwise will treat as if no chimeric trimming required
	{
	MinPutLen = ((ReadLen * MinChimericLen) + 50) / 100;
	if(CoreLen > MinPutLen)
		CoreLen = MinPutLen;
	}
else
	MinPutLen = ReadLen;
if (MinPutLen == ReadLen)
	{
	MinChimericLen = 0;
	CoreLen = 0;
	}

// ensure insert sizes are reasonable, min can't be more than max
if(MinInsertSize > MaxInsertSize)
	return(0);
// force minimum has to be a least the untrimmed read length
 if(MinInsertSize < ReadLen)
	{
	MaxInsertSize += ReadLen - MinInsertSize;
	MinInsertSize = ReadLen;
	}

// check if any putative extension would still be within targeted ChromSeqLen
if(b3primeExtend)
	{
	if((uint32_t)(StartLoci + MinInsertSize) >= ChromSeqLen)	// MaxInsertSize could extend over chromosome boundary so will be explicitly later checking for this...
		return(0);
	StartPutTargLoci = StartLoci + MinInsertSize - ReadLen;
	EndPutTargLoci = min(ChromSeqLen-ReadLen,(uint32_t)(StartLoci + MaxInsertSize - ReadLen));
	}
else // else must be 5' extension
	{
	if(EndLoci < (uint32_t)MinInsertSize)
		return(0);
	if(EndLoci <= (uint32_t)MaxInsertSize)
		StartPutTargLoci = 0;
	else
		StartPutTargLoci = EndLoci - MaxInsertSize;
	EndPutTargLoci = EndLoci - MinInsertSize;
	}

if((pPutChromSeq = GetPtrSeq(ChromID))==NULL)
	return(0);

memmove(ReadSeq,pRead,ReadLen);
ReadSeq[ReadLen] = eBaseEOS;
pRead = ReadSeq;

// if aligning antisense then revcpl the read sequence
if(bAntisense)
	{
	if(m_bColorspace)
		CSeqTrans::ReverseSeq(ReadLen,ReadSeq);
	else
		CSeqTrans::ReverseComplement(ReadLen,ReadSeq);
	}

LowestNumMM = MaxAllowedMM+2;
LowestMMcnt = 0;

int Rslt;
uint32_t TrimSeqLen;
uint32_t Trim5Flank;
uint32_t Trim3Flank;
uint32_t MaxChimericMMs;

uint32_t PrevBestMaxChimericMMs = MaxAllowedMM + 1;
if((EndPutTargLoci - StartPutTargLoci) >= 1000)		// if a relatively large insert size (likely with WGS RNA alignments) then use suffix array to locate exactly match cores which can serve as loci for adaptive trimming the read against
	{
	for(PutCoreOfs = 0; ((int)PutCoreOfs + CoreLen) <= ReadLen; PutCoreOfs += (uint32_t)CoreDelta)
		{
		PutHitIdx = 0;
		NxtPutHitIdx = 0;
		while ((NxtPutHitIdx = IterateExactsRange(&ReadSeq[PutCoreOfs],CoreLen, PutHitIdx,ChromID, StartPutTargLoci, EndPutTargLoci,&PutHitLoci,&PutChromLen,&pPutChromSeq)) > 0)
			{
			if(PutCoreOfs > PutHitLoci || (PutHitLoci + ReadLen - PutCoreOfs) >= ChromSeqLen)
				continue;

			Rslt = AdaptiveTrim(ReadLen,ReadSeq,&pPutChromSeq[PutHitLoci - PutCoreOfs], ChromID,ChromSeqLen,MinPutLen,MaxAllowedMM,3,&TrimSeqLen,&Trim5Flank,&Trim3Flank,&MaxChimericMMs);
			if(Rslt > (int)MinPutLen || (Rslt == (int)MinPutLen && MaxChimericMMs < PrevBestMaxChimericMMs))
				{
				 PrevBestMaxChimericMMs = MaxChimericMMs;
				 MinPutLen = Rslt;
				 AlignLoci = PutHitLoci - PutCoreOfs;
				 pAlign->FlgChimeric = MinPutLen == ReadLen ? 0 : 1;
				 pAlign->FlgInDel = 0;
				 pAlign->FlgInsert = 0;
				 pAlign->FlgSplice = 0;
				 pAlign->FlgNonOrphan = 0;
				 pAlign->Seg[0].Strand = bAntisense ? '-' : '+';
				 if (!bAntisense)
					{
					 pAlign->Seg[0].TrimLeft = Trim5Flank;
					 pAlign->Seg[0].TrimRight = Trim3Flank;
					}
				 else
					{
					 pAlign->Seg[0].TrimLeft = Trim3Flank;
					 pAlign->Seg[0].TrimRight = Trim5Flank;
					}
				 pAlign->Seg[0].ChromID = ChromID;
				 pAlign->Seg[0].ReadOfs = 0;
				 pAlign->Seg[0].MatchLoci = (uint32_t)(PutHitLoci - PutCoreOfs);
				 pAlign->Seg[0].MatchLen = ReadLen;
				 pAlign->Seg[0].Mismatches = MaxChimericMMs;
				 pAlign->Seg[0].TrimMismatches = MaxChimericMMs;
				}
			PutHitIdx = NxtPutHitIdx;
			}
		}
	}
else  // else if relatively small maximum insert size ( < 1000bp ) as is likely with WGS or GBS DNA then quicker to do a linear scan for target loci against which the read can be matched
	{
	for(PutHitLoci = StartPutTargLoci; PutHitLoci <= EndPutTargLoci; PutHitLoci++)
		{
		Rslt = AdaptiveTrim(ReadLen, ReadSeq, &pPutChromSeq[PutHitLoci], ChromID,ChromSeqLen,MinPutLen, MaxAllowedMM, 3, &TrimSeqLen, &Trim5Flank, &Trim3Flank, &MaxChimericMMs);
		if (Rslt > (int)MinPutLen || (Rslt == (int)MinPutLen && MaxChimericMMs < PrevBestMaxChimericMMs))
			{
			PrevBestMaxChimericMMs = MaxChimericMMs;
			MinPutLen = Rslt;
			AlignLoci = PutHitLoci;
			pAlign->FlgChimeric = MinPutLen == ReadLen ? 0 : 1;
			pAlign->FlgInDel = 0;
			pAlign->FlgInsert = 0;
			pAlign->FlgSplice = 0;
			pAlign->FlgNonOrphan = 0;
			pAlign->Seg[0].Strand = bAntisense ? '-' : '+';
			if (!bAntisense)
			{
				pAlign->Seg[0].TrimLeft = Trim5Flank;
				pAlign->Seg[0].TrimRight = Trim3Flank;
			}
			else
			{
				pAlign->Seg[0].TrimLeft = Trim3Flank;
				pAlign->Seg[0].TrimRight = Trim5Flank;
			}
			pAlign->Seg[0].ChromID = ChromID;
			pAlign->Seg[0].ReadOfs = 0;
			pAlign->Seg[0].MatchLoci = PutHitLoci;
			pAlign->Seg[0].MatchLen = ReadLen;
			pAlign->Seg[0].Mismatches = MaxChimericMMs;
			pAlign->Seg[0].TrimMismatches = MaxChimericMMs;
			if(MinPutLen == ReadLen && MaxChimericMMs == 0)
				break;
			}
		}
	}
return(PrevBestMaxChimericMMs <= (uint32_t)MaxAllowedMM ? 1 : 0);
}


int									// -1 if errors, 0 if no match, 1 if mismatches only, 2 InDel and insertion into probe, 3 InDel and deletion from probe
CSfxArray::ExploreSpliceRight(char CurStrand,					// aligning on this strand
			int MaxSpliceJunctLen,		// junction has to be no longer than this length
			int MaxTotMM,			// can be at most this many mismatches in total
			int CoreLen,			// core length used
 		   int ProbeLen,			// length of probe excluding any eBaseEOS
		   etSeqBase *pProbe,		// pts to 5' start of probe sequence
		   int64_t TargOfs,			// pTarg corresponds to this suffix offset
		   int64_t TargLen,			// max length of target to be explored
		   etSeqBase *pTarg,		// pts to target sequence
		   tsHitLoci *pHit)			// where to return hit loci
{
uint32_t CurScore;					// current score for putative splice being explored

tsHitLoci CurSplice;				// current splice junction loci

uint32_t Idx;
etSeqBase *pP,*pT;
etSeqBase *pCurP;

etSeqBase ProbeBase,TargBase;
int CurNumMismatches;

uint32_t MaxSegLen;

uint32_t TmpTargLen;
int CurMismatches;
int CurSpliceSeqLen;
etSeqBase Donor[2];
etSeqBase Accept[2];

int64_t TargSeqLen;

int ProbeMMOfss[cMaxPutInDelOfss];			// offsets of all mismatches in probe which may need to be explored as putative InDel starts
int MMIdx;									// current index into ProbeMMOfss
int TotMM;									// total mismatches to be explored as putative InDels

memset(pHit,0,sizeof(tsHitLoci));
if((TargOfs + ProbeLen + cMinJunctAlignSep) > TargLen)
	return(0);

memset(&CurSplice,0,sizeof(tsHitLoci));

if(MaxTotMM > cMaxJunctAlignMM)		// silently clamp
	MaxTotMM = cMaxJunctAlignMM;

TargSeqLen = TargLen;

pP = pProbe + CoreLen;
pT = pTarg + CoreLen;

CurNumMismatches = 0;

for(Idx = CoreLen; Idx < (uint32_t)ProbeLen && CurNumMismatches <= max(MaxTotMM,cMaxJunctAlignMM * 5); Idx++,pT+=1,pP+=1)
	{
	ProbeBase = *pP & 0x07;
	TargBase = *pT & 0x07;
	if(TargBase > eBaseN || ProbeBase > eBaseN)	// treat as no match any EOS
		return(0);

	if(ProbeBase == TargBase && ProbeBase <= eBaseT)
		continue;

		// have a mismatch
	ProbeMMOfss[CurNumMismatches] = Idx;		// record probe offset at which this mismatch occurred as this is where putative splice junctions may be subsequently explored
	CurNumMismatches += 1;
	}

if(CurNumMismatches < (cMaxJunctAlignMM * 4) ||		    // only bother exploring putative splice junctions if too many mismatches
	cMinJunctSegLen	> (ProbeLen - ProbeMMOfss[0]))	// and the length of sequence which could be in a spliced segment is worth the effort
	{
	if(CurNumMismatches > MaxTotMM)
		return(0);

	memset(pHit->Seg,0,sizeof(pHit->Seg));
	pHit->Seg[0].MatchLen = ProbeLen;
	pHit->Seg[0].MatchLoci = TargOfs;
	pHit->Seg[0].Mismatches = CurNumMismatches;
	pHit->Seg[0].TrimMismatches = CurNumMismatches;
	pHit->Seg[0].ReadOfs = 0;
	pHit->Seg[0].Strand = CurStrand;
	pHit->Score = cBaseScore + (ProbeLen * cScoreMatch) - (CurNumMismatches * cScoreMismatch);
	return(1);
	}

TotMM = min(CurNumMismatches,MaxTotMM);  // can accept at most this many mismatches in right segment

etSeqBase *pTStart;			// pts to start of putative right target segment
etSeqBase *pTEnd;			// pts to end of putative right target segment
etSeqBase *pTDonor;			// pts to putative splice donor site (GT) to right of left target segment
int MaxHashDiff;
int ProbeSegHash;
int MinTargSegHash;
int MaxTargSegHash;
int TargSegHash;

// explore putative junctions and if so then choose the highest scoring junction as the one to report
// look right on target and ensure that start of any potential right splice segment will be on same chrom...
pT = &pTarg[ProbeMMOfss[TotMM-1]];
for(Idx = 0; Idx < (uint32_t)(cMinJunctAlignSep+cMinJunctSegLen); Idx++,pT+=1)
	{
	if((*pT & 0x07) > eBaseN)
		return(0);
	}
int	InHash = 0;
int OutHash = 0;

// iterate each putative splice site (mismatch) in probe
for(MMIdx = 0; MMIdx <= TotMM && cMinJunctSegLen	< (ProbeLen - ProbeMMOfss[MMIdx]); MMIdx++)
	{
	if(CurSplice.Score >= cMaxScore)	// unlikely to do any better than this so early exit from loop
		break;
	MaxSegLen = ProbeLen - ProbeMMOfss[MMIdx];
	pCurP = &pProbe[ProbeMMOfss[MMIdx]];			// for this putative splice junction, pCurP pts at 1st base of right probe segment
	pTDonor = &pTarg[ProbeMMOfss[MMIdx]];			// if splice junction then could be a donor (GT) site


	pTStart = pTDonor + cMinJunctAlignSep;			// if splice junction then could be start of right segment
	MaxHashDiff = 4 * (MaxTotMM - MMIdx);

	// calculate hash for probe right segment
	pP = pCurP;
	ProbeSegHash = 100000;
	for(Idx = 0; Idx < MaxSegLen; Idx++)
		ProbeSegHash += *pP++ & 0x07;
	MinTargSegHash = ProbeSegHash - MaxHashDiff; // will only bother exploring full segment if target hash between MinProbeSegHash and MaxProbeSegHash
	MaxTargSegHash = ProbeSegHash + MaxHashDiff;
	// calculate initial hash for target segment - 1
	TargSegHash = 100000;
	pTEnd = pTStart;
	for(Idx = 0; Idx < MaxSegLen-1; Idx++)
		{
		if((*pTEnd & 0x07) > eBaseN)
			break;
		TargSegHash += *pTEnd++ & 0x07;
		}
	if(Idx < (MaxSegLen-1))
		break;

	// now iterate over target, but only bother fully exploring if target segment hash in range of min/max hash
	for(CurSpliceSeqLen = cMinJunctAlignSep;CurSpliceSeqLen < MaxSpliceJunctLen - (int)MaxSegLen; CurSpliceSeqLen++,pTStart++,pTEnd++)
		{
		if((TargBase = *pTEnd & 0x07) > eBaseN)
			break;
		TargSegHash += TargBase;
		if(TargSegHash < MinTargSegHash || TargSegHash > MaxTargSegHash)
			{
			TargSegHash -= (*pTStart & 0x07);
			InHash += 1;
			continue;
			}
		OutHash += 1;
		// hash in range, now worth exploring this putative right splice segment
		TargSegHash -= (*pTStart & 0x07);
		pT = pTStart;
		pP = pCurP;
		TmpTargLen = (uint32_t)(TargSeqLen - (TargOfs + ProbeMMOfss[MMIdx] + CurSpliceSeqLen + 1));
		if(TmpTargLen < MaxSegLen)
			break;

		CurMismatches = 0;
		for(Idx = 0; Idx < MaxSegLen && (MMIdx + CurMismatches) < MaxTotMM; Idx++,pT+=1,pP+=1)
			{
			ProbeBase = *pP & 0x07;
			TargBase = *pT & 0x07;
			if(ProbeBase > eBaseN || TargBase > eBaseN)
				break;

			if(ProbeBase == TargBase && ProbeBase <= eBaseT)
				continue;

			CurMismatches += 1;
			}
		if(Idx != MaxSegLen)
			{
			if(ProbeBase > eBaseN || TargBase > eBaseN)
				break;
			continue;
			}
		CurScore = cBaseScore + (ProbeLen * cScoreMatch) - (((MMIdx + CurMismatches) * cScoreMismatch) + ((CurSpliceSeqLen/1000) *  cSpliceLen));
		// in colorspace need to get basespace for comparison against the cannonical donor/acceptor sites
		if(m_bColorspace)
			{
			Donor[0] = (*pTDonor >> 4) & 0x07;
			Donor[1] = (pTDonor[1] >> 4) & 0x07;
			Accept[0] = (pTStart[-1] >> 4) & 0x07;
			Accept[1] = (pTStart[-2] >> 4) & 0x07;
			}
		else
			{
			Donor[0] = *pTDonor & 0x07;
			Donor[1] = pTDonor[1] & 0x07;
			Accept[0] = pTStart[-1] & 0x07;
			Accept[1] = pTStart[-2] & 0x07;
			}

		// following splice scoring is a little convoluted because we really are not too sure as to what strand a read really originated from
		// especially if Illumina so we first try the strand canonical splicesites and if no match then try the antisense splicesites but give
		// higher score if on the expected strand
		if(CurStrand == '+')
			{
			if(((Donor[0] == eBaseG) && (Donor[1] == eBaseT)) && (Accept[0] == eBaseG && Accept[1] == eBaseA))
				CurScore += cSpliceDonorAccept;
			else
				if(((Donor[0] == eBaseC) && (Donor[1] == eBaseT)) && (Accept[0] == eBaseC && Accept[1] == eBaseA))
					CurScore += cSpliceDonorAccept/2;
			}
		else
			{
			if(((Donor[0] == eBaseC) && (Donor[1] == eBaseT)) && (Accept[0] == eBaseC && Accept[1] == eBaseA))
				CurScore += cSpliceDonorAccept;
			else
				if(((Donor[0] == eBaseG) && (Donor[1] == eBaseT)) && (Accept[0] == eBaseG && Accept[1] == eBaseA))
					CurScore += cSpliceDonorAccept/2;
			}

		if(CurScore > CurSplice.Score)					// exploring right
			{
			memset(CurSplice.Seg,0,sizeof(CurSplice.Seg));
			CurSplice.Seg[0].MatchLen = ProbeMMOfss[MMIdx];
			CurSplice.Seg[0].MatchLoci = TargOfs;		// target suffix, not chrom loci
			CurSplice.Seg[0].Mismatches = MMIdx;
			CurSplice.Seg[0].TrimMismatches = MMIdx;
			CurSplice.Seg[0].Strand = CurStrand;

			CurSplice.Seg[1].MatchLen = ProbeLen - ProbeMMOfss[MMIdx];
			CurSplice.Seg[1].MatchLoci =	TargOfs + ProbeMMOfss[MMIdx] + CurSpliceSeqLen;		// target suffix, not chrom loci
			CurSplice.Seg[1].Mismatches = CurMismatches;
			CurSplice.Seg[1].TrimMismatches = CurMismatches;
			CurSplice.Seg[1].ReadOfs = ProbeMMOfss[MMIdx];
			CurSplice.Seg[1].Strand = CurStrand;
			CurSplice.Score = CurScore;
			CurSplice.FlgSplice = 1;
			}
		}
	}

if(CurSplice.Score == 0)	//  will be 0 if no splice junction accepted
	return(0);
*pHit = CurSplice;
return(3);
}

// ExploreSpliceLeft
// Normally this method is only invoked if unable to align using LocateCoreMultiples, e.g. exact or near exact matches, which
// would not have aligned reads in which there were micro InDels (small InDels in range of 1..microInDelLen, currently defined max is 20)
// micro InDels are only explored if there are more than cMaxMMExploreInDel (currently defined as 7) mismatches between the read and any other subsequence
// of same length in the targeted genome. The flanking length also needs to be at least cMinInDelSeqLen (currently defined as 7).
// Insertions (or deletions from the targeted sequence) into the read sequence are first scored, then deletions (insertions into target) from
// the read are scored with the lowest scored InDel chosen. Currently scoring is 2 for a mismatch, 10 for starting an InDel and 1 for InDel extension per base extended.
//
int									// -1 if errors, 0 if no match, 1 if mismatches only, 2 InDel and insertion into probe, 3 InDel and deletion from probe
CSfxArray::ExploreSpliceLeft(char CurStrand,					// aligning on this strand
			int MaxSpliceJunctLen,		// junction has to be no longer than this length
			int MaxTotMM,			// can be at most this many mismatches in total
			int CoreLen,			// core length used
 		   int ProbeLen,			// length of probe excluding any eBaseEOS
		   etSeqBase *pProbe,		// pts to 5' start of probe sequence
		   int64_t TargOfs,			// pTarg corresponds to this suffix offset
		   int64_t TargLen,			// max length of target to be explored
		   etSeqBase *pTarg,		// pts to target sequence, 5' 1st base if bRight, 3' last base if !bRight
		   tsHitLoci *pHit)			// where to return hit loci
{
uint32_t CurScore;					// current score for putative InDel being explored
tsHitLoci CurSplice;				// current splice loci

uint32_t Idx;
etSeqBase *pP,*pT;
etSeqBase *pCurP;

etSeqBase ProbeBase,TargBase;
int CurNumMismatches;


uint32_t TmpTargLen;
int CurMismatches;
int CurSpliceSeqLen;
etSeqBase Donor[2];
etSeqBase Accept[2];

int64_t TargSeqLen;
uint32_t MaxSegLen;

int ProbeMMOfss[cMaxPutInDelOfss];			// offsets of all mismatches in probe which may need to be explored as putative InDel starts
int MMIdx;									// current index into ProbeMMOfss
int TotMM;									// total mismatches to be explored as putative InDels

memset(pHit,0,sizeof(tsHitLoci));

if(TargOfs < (uint32_t)(cMinJunctAlignSep + cMinJunctSegLen))
	return(0);

memset(&CurSplice,0,sizeof(tsHitLoci));

if(MaxTotMM > cMaxJunctAlignMM)		// silently clamp
	MaxTotMM = cMaxJunctAlignMM;

						// when exploring left then need to set pProbe to pt at 3' end of probe
pProbe += ProbeLen - 1;
pTarg += ProbeLen - 1;
TargSeqLen = TargLen;

pP = pProbe - CoreLen;
pT = pTarg - CoreLen;

// count mismatches and record where these are, if too many mismatches then will explore for an InDel starting at each mismatch
CurNumMismatches = 0;
for(Idx = CoreLen; Idx < (uint32_t)ProbeLen && CurNumMismatches <= max(MaxTotMM,cMaxJunctAlignMM * 5); Idx++,pT-=1,pP-=1)
	{
	ProbeBase = *pP & 0x07;
	TargBase = *pT & 0x07;

	if(TargBase > eBaseN || ProbeBase > eBaseN)	// treat as no match any EOS
		return(0);

	if(ProbeBase == TargBase && ProbeBase <= eBaseT)
		continue;

		// have a mismatch
 	ProbeMMOfss[CurNumMismatches] = Idx;		// record probe offset at which this mismatch occured as this is where putative InDels may be subsequently explored
	CurNumMismatches += 1;
	}

if(CurNumMismatches < (cMaxJunctAlignMM * 4) ||		    // only bother exploring putative splice junctions if too many mismatches
	cMinJunctSegLen	> (ProbeLen - ProbeMMOfss[0]))	// and the length of sequence which could be in a spliced segment is worth the effort
	{
	if(CurNumMismatches > MaxTotMM)
		return(0);
	memset(pHit->Seg,0,sizeof(pHit->Seg));
	pHit->Seg[0].MatchLen = ProbeLen;
	pHit->Seg[0].MatchLoci = TargOfs;
	pHit->Seg[0].Mismatches = CurNumMismatches;
	pHit->Seg[0].TrimMismatches = CurNumMismatches;
	pHit->Seg[0].ReadOfs = 0;
	pHit->Seg[0].Strand = CurStrand;
	pHit->Score = cBaseScore + (ProbeLen * cScoreMatch) - (CurNumMismatches * cScoreMismatch);
	return(1);
	}
TotMM = min(CurNumMismatches,MaxTotMM); // can accept at most this many mismatches

etSeqBase *pTStart;			// pts to start of putative left target segment
etSeqBase *pTEnd;			// pts to end of putative left target segment
etSeqBase *pTDonor;			// pts to putative splice donor site (GT)
int MaxHashDiff;
int ProbeSegHash;
int MinTargSegHash;
int MaxTargSegHash;
int TargSegHash;

// explore putative junctions and if so then choose the highest scoring junction as the one to report

// look left on target and ensure that start of any potential splice segment will be on same chrom...
pT = &pTarg[-1 * ProbeMMOfss[TotMM]];
for(Idx = 0; Idx < (uint32_t)(cMinJunctAlignSep+cMinJunctSegLen); Idx++,pT-=1)
	{
	if((*pT & 0x07) > eBaseN)
		return(0);
	}

// iterate each putative splice site (mismatch) in probe
for(MMIdx = 0; MMIdx <= TotMM && cMinJunctSegLen	< (ProbeLen - ProbeMMOfss[MMIdx]); MMIdx++)
	{
	if(CurSplice.Score >= cMaxScore)	// unlikely to do any better than this so early exit from loop
		break;
	MaxSegLen = ProbeLen - ProbeMMOfss[MMIdx];
	pCurP = &pProbe[-1 * ProbeMMOfss[MMIdx]];			// for this putative splice junction, pCurP pts at 1st base of right probe segment
	pTDonor = &pTarg[-1 * ProbeMMOfss[MMIdx]];			// if splice junction then could be a donor (GT) site

	pTStart = pTDonor - cMinJunctAlignSep;		// if splice junction then could be start of right segment

	MaxHashDiff = 4 * (MaxTotMM - MMIdx);

	// calculate hash for probe left segment
	pP = pCurP;
	ProbeSegHash = 100000;
	for(Idx = 0; Idx < MaxSegLen; Idx++)
		ProbeSegHash += *pP-- & 0x07;
	MinTargSegHash = ProbeSegHash - MaxHashDiff; // will only bother exploring full segment if target hash between MinProbeSegHash and MaxProbeSegHash
	MaxTargSegHash = ProbeSegHash + MaxHashDiff;

	// calculate initial hash for target segment - 1
	TargSegHash = 100000;
	pTEnd = pTStart;
	for(Idx = 0; Idx < MaxSegLen-1; Idx++)
		{
		if((*pTEnd & 0x07) > eBaseN)
			break;
		TargSegHash += *pTEnd-- & 0x07;
		}
	if(Idx < (MaxSegLen-1))
		break;

	// now iterate over target, but only bother fully exploring if target segment hash in range of min/max hash
	for(CurSpliceSeqLen = cMinJunctAlignSep;CurSpliceSeqLen < MaxSpliceJunctLen - (int)MaxSegLen; CurSpliceSeqLen++,pTStart--,pTEnd--)
		{
		if((TargBase = *pTEnd & 0x07) > eBaseN)
			break;
		TargSegHash += TargBase;
		if(TargSegHash < MinTargSegHash || TargSegHash > MaxTargSegHash)
			{
			TargSegHash -= (*pTStart & 0x07);
			continue;
			}
		// hash in range, now worth exploring this putative right splice segment
		TargSegHash -= (*pTStart & 0x07);
		pT = pTStart;
		pP = pCurP;


		TmpTargLen = (uint32_t)(TargOfs - (CurSpliceSeqLen));
		if(TmpTargLen < 1)
			break;

		CurMismatches = 0;
		for(Idx = 0; Idx < MaxSegLen && (MMIdx + CurMismatches) < MaxTotMM; Idx++,pT-=1,pP-=1)
			{
			ProbeBase = *pP & 0x07;
			TargBase = *pT & 0x07;
			if(ProbeBase > eBaseN || TargBase > eBaseN)
				break;

			if(ProbeBase == TargBase && ProbeBase <= eBaseT)
				continue;

			CurMismatches += 1;
			}
		if(Idx != MaxSegLen)
			{
			if(ProbeBase > eBaseN || TargBase > eBaseN)
				break;
			continue;
			}
		CurScore = cBaseScore + (ProbeLen * cScoreMatch) - (((MMIdx + CurMismatches) * cScoreMismatch) + ((CurSpliceSeqLen/1000) *  cSpliceLen));

		if(m_bColorspace)
			{
			Accept[0] = (*pTDonor >> 4) & 0x07;
			Accept[1] = (pTDonor[-1] >> 4) & 0x07;
			Donor[0] = (pTStart[1] >> 4) & 0x07;
			Donor[1] = (pTStart[2] >> 4) & 0x07;
			}
		else
			{
			Accept[0] = *pTDonor & 0x07;
			Accept[1] = pTDonor[-1] & 0x07;
			Donor[0] = pTStart[1] & 0x07;
			Donor[1] = pTStart[2] & 0x07;
			}
		// following splice scoring is a little convoluted because we really are not too sure as to what strand a read really originated from
		// especially if Illumina so we first try the strand cannonical splicesites and if no match then try the antisense splicesites but give
		// higher score if on the expected strand
		if(CurStrand == '+')
			{
			if(((Donor[0] == eBaseG) && (Donor[1] == eBaseT)) && (Accept[0] == eBaseG && Accept[1] == eBaseA))
				CurScore += cSpliceDonorAccept;
			else
				if(((Donor[0] == eBaseC) && (Donor[1] == eBaseT)) && (Accept[0] == eBaseC && Accept[1] == eBaseA))
					CurScore += cSpliceDonorAccept/2;
			}
		else
			{
			if(((Donor[0] == eBaseC) && (Donor[1] == eBaseT)) && (Accept[0] == eBaseC && Accept[1] == eBaseA))
				CurScore += cSpliceDonorAccept;
			else
				if(((Donor[0] == eBaseG) && (Donor[1] == eBaseT)) && (Accept[0] == eBaseG && Accept[1] == eBaseA))
					CurScore += cSpliceDonorAccept/2;
			}

		if(CurScore > CurSplice.Score)					// exploring left
			{
			memset(CurSplice.Seg,0,sizeof(CurSplice.Seg));
			CurSplice.Seg[0].MatchLen = ProbeLen - ProbeMMOfss[MMIdx];
			CurSplice.Seg[0].MatchLoci =	TargOfs - (CurSpliceSeqLen);		// target suffix, not chrom loci
			CurSplice.Seg[0].Mismatches = CurMismatches;
			CurSplice.Seg[0].TrimMismatches = CurMismatches;
			CurSplice.Seg[0].Strand = CurStrand;

			CurSplice.Seg[1].MatchLen = ProbeMMOfss[MMIdx];
			CurSplice.Seg[1].MatchLoci = CurSplice.Seg[0].MatchLoci	+ CurSplice.Seg[0].MatchLen + CurSpliceSeqLen;		// target suffix, not chrom loci
			CurSplice.Seg[1].Mismatches = MMIdx;
			CurSplice.Seg[1].TrimMismatches = MMIdx;
			CurSplice.Seg[1].ReadOfs = CurSplice.Seg[0].MatchLen;
			CurSplice.Seg[1].Strand = CurStrand;
			CurSplice.Score = CurScore;
			CurSplice.FlgSplice = 1;
			}
		}
	}

if(CurSplice.Score == 0)	// will be less than 0xffff if splice junction located
	return(0);

// have accepted this splice junction
*pHit = CurSplice;
return(3);
}


// ExploreInDelMatchRight
// Normally this method is only invoked if unable to align using LocateCoreMultiples, e.g. exact or near exact matches, which
// would not have aligned reads in which there were micro InDels (small InDels in range of 1..microInDelLen, currently defined max is 20)
// micro InDels are only explored if there are more than cMaxMMExploreInDel (currently defined as 7) mismatches between the read and any other subsequence
// of same length in the targeted genome. The flanking length also needs to be at least cMinInDelSeqLen (currently defined as 7).
// Insertions (or deletions from the targeted sequence) into the read sequence are first scored, then deletions (insertions into target) from
// the read are scored with the lowest scored InDel chosen. Currently scoring is 2 for a mismatch, 10 for starting an InDel and 1 for InDel extension per base extended.
//
int									// -1 if errors, 0 if no match, 1 if mismatches only, 2 InDel and insertion into probe, 3 InDel and deletion from probe
CSfxArray::ExploreInDelMatchRight(uint32_t ExtdProcFlags, // flags indicating if lower levels need to do any form of extended processing with this specific read...
						 uint32_t ReadID,					  // identifies this read
						 char CurStrand,					// aligning on this strand
			int microInDelLen,		   		// microInDel can be upto (inclusive) this length
			int MaxTotMM,			// can be at most this many mismatches in total
 		   int ProbeLen,			// length of probe excluding any eBaseEOS
		   etSeqBase *pProbe,		// pts to 5' start of probe sequence
		   tsSfxEntry *pEntry,		// target at TargSeqLeftIdx is in this sfx entry
		   int64_t TargOfs,			// pTarg corresponds to this suffix offset
		   etSeqBase *pTarg,		// pts to target sequence
		   tsHitLoci *pHit)			// where to return hit loci
{
int CurScore;						// current score for putative InDel being explored

tsHitLoci CurInsert;				// current insert loci
tsHitLoci CurDelete;				// current delete loci

uint32_t Idx;
etSeqBase *pP,*pT;
etSeqBase ProbeBase,TargBase;
int CurNumMismatches;

uint32_t TmpProbeLen;
uint32_t TmpTargLen;
int CurInDelMismatches;
int CurInDelLen;

uint32_t TargSeqLen;

int ProbeMMOfss[cMaxPutInDelOfss];			// offsets of all mismatches in probe which may need to be explored as putative InDel starts
int MMIdx;									// current index into ProbeMMOfss
int TotMM;									// total mismatches to be explored as putative InDels

memset(pHit,0,sizeof(tsHitLoci));
if(TargOfs < (int64_t)pEntry->StartOfs || ((TargOfs + ProbeLen - 1) > (int64_t)pEntry->EndOfs))
	return(0);

TargSeqLen = (uint32_t)(pEntry->SeqLen -  (TargOfs - pEntry->StartOfs));

memset(&CurInsert,0,sizeof(tsHitLoci));
memset(&CurDelete,0,sizeof(tsHitLoci));

if(MaxTotMM > cMaxPutInDelOfss)		// silently clamp so as to ensure ProbeMMOfss never overflowed
	MaxTotMM = cMaxPutInDelOfss;

pP = pProbe;
pT = pTarg;

// count mismatches and record where these are, if too many mismatches then will explore for an InDel starting at each mismatch
CurNumMismatches = 0;
for(Idx = 0; Idx < (uint32_t)ProbeLen && CurNumMismatches <= max(MaxTotMM,cMaxMMExploreInDel); Idx++,pT+=1,pP+=1)
	{
	ProbeBase = *pP & 0x07;
	TargBase = *pT & 0x07;

	if(TargBase > eBaseN || ProbeBase > eBaseN)	// treat as no match any EOS
		return(0);

	if(ProbeBase == TargBase && ProbeBase <= eBaseT)
		continue;


		// have a mismatch
	ProbeMMOfss[CurNumMismatches] = Idx;		// record probe offset at which this mismatch occurred as this is where putative InDels may be subsequently explored
	CurNumMismatches += 1;
	}

if(CurNumMismatches < cMaxMMExploreInDel ||		    // only bother exploring putative InDels if too many mismatches
	cMinInDelSeqLen	> (ProbeLen - ProbeMMOfss[0]))	// and the length of sequence which could be InDel'd is worth the effort
	{
	if(CurNumMismatches > MaxTotMM)
		return(0);
	memset(pHit->Seg,0,sizeof(pHit->Seg));
	pHit->Seg[0].MatchLen = ProbeLen;
	pHit->Seg[0].MatchLoci = TargOfs;
	pHit->Seg[0].Mismatches = CurNumMismatches;
	pHit->Seg[0].TrimMismatches = CurNumMismatches;
	pHit->Seg[0].Strand = CurStrand;
	pHit->Score = cBaseScore + (ProbeLen * cScoreMatch) - (CurNumMismatches * cScoreMismatch);
	return(1);
	}



// too many mismatches so now explore for insertions into probe starting from each mismatch
TotMM = min(MaxTotMM,CurNumMismatches);
for(MMIdx = 0; MMIdx <= TotMM && cMinInDelSeqLen	< (ProbeLen - ProbeMMOfss[MMIdx]); MMIdx++)
	{
	// microInDels are to be limited in length
	for(CurInDelLen = 1; CurInDelLen <= microInDelLen; CurInDelLen++)
		{
		CurScore = cBaseScore + (ProbeLen * cScoreMatch) - (((CurInDelLen - 1) * cScoreInDelExtn) + cScoreInDelOpn);

		pP = &pProbe[ProbeMMOfss[MMIdx]+CurInDelLen];
		pT = &pTarg[ProbeMMOfss[MMIdx]];

		TmpProbeLen = ProbeLen - (ProbeMMOfss[MMIdx] + CurInDelLen);
		if(TmpProbeLen < cMinInDelSeqLen)
			break;
		TmpTargLen = TargSeqLen - ProbeMMOfss[MMIdx];
		if(TmpTargLen < TmpProbeLen)
			break;
		CurInDelMismatches = 0;
		for(Idx = 0; Idx < TmpProbeLen &&
					 (MMIdx + CurInDelMismatches) <= MaxTotMM;
					  Idx++,pT+=1,pP+=1)
			{
			ProbeBase = *pP & 0x07;
			TargBase = *pT & 0x07;
			if(ProbeBase > eBaseN || TargBase > eBaseN)
				break;

			if(ProbeBase == TargBase && ProbeBase <= eBaseT)
				continue;

			CurInDelMismatches += 1;
			CurScore -= cScoreMismatch;
			if(TmpProbeLen < (uint32_t)(cMinInDelSeqLen * CurInDelMismatches))
				break;
			}
	    if(Idx != TmpProbeLen)
			continue;

		if(CurScore > CurInsert.Score)					// explore right, insertion into probe...
			{
			memset(CurInsert.Seg,0,sizeof(CurInsert.Seg));
			CurInsert.Seg[0].MatchLen = ProbeMMOfss[MMIdx];
			CurInsert.Seg[0].MatchLoci =	TargOfs;		// target suffix, not chrom loci
			CurInsert.Seg[0].Mismatches = MMIdx;
			CurInsert.Seg[0].TrimMismatches = MMIdx;
			CurInsert.Seg[0].Strand = CurStrand;

			CurInsert.Seg[1].MatchLen = TmpProbeLen;
			CurInsert.Seg[1].MatchLoci =	CurInsert.Seg[0].MatchLoci + CurInsert.Seg[0].MatchLen;		// target suffix, not chrom loci
			CurInsert.Seg[1].Mismatches = CurInDelMismatches;
			CurInsert.Seg[1].TrimMismatches = CurInDelMismatches;
			CurInsert.Seg[1].ReadOfs = CurInsert.Seg[0].MatchLen + CurInDelLen;
			CurInsert.Seg[1].Strand = CurStrand;
			CurInsert.Score = CurScore;
			CurInsert.FlgChimeric = 0;
			CurInsert.FlgInDel = 1;
			CurInsert.FlgInsert = 1;
			CurInsert.FlgSplice = 0;
			CurInsert.FlgNonOrphan = 0;
			}
		}
	}

// now know if there was a putative insertion, check if there is a putative deletion and if so then choose the highest scoring InDel as the one to report
for(MMIdx = 0; MMIdx <= TotMM && cMinInDelSeqLen	< (ProbeLen - ProbeMMOfss[MMIdx]); MMIdx++)
	{
	for(CurInDelLen = 1; CurInDelLen <= microInDelLen; CurInDelLen++)
		{
		CurScore = cBaseScore + (ProbeLen * cScoreMatch) - (((CurInDelLen - 1) * cScoreInDelExtn) + cScoreInDelOpn);

		pP = &pProbe[ProbeMMOfss[MMIdx]];
		pT = &pTarg[ProbeMMOfss[MMIdx]+CurInDelLen];
		TmpProbeLen = ProbeLen - ProbeMMOfss[MMIdx];
		if(TmpProbeLen < cMinInDelSeqLen)
			break;
		TmpTargLen = TargSeqLen - (ProbeMMOfss[MMIdx]+CurInDelLen);
		if(TmpTargLen < TmpProbeLen)
			break;
		CurInDelMismatches = 0;
		for(Idx = 0; Idx < TmpProbeLen && (MMIdx + CurInDelMismatches) <= MaxTotMM; Idx++,pT+=1,pP+=1)
			{
			ProbeBase = *pP & 0x07;
			TargBase = *pT & 0x07;
			if(ProbeBase > eBaseN || TargBase > eBaseN)
				break;

			if(ProbeBase == TargBase && ProbeBase <= eBaseT)
				continue;

			CurInDelMismatches += 1;
			CurScore -= cScoreMismatch;
			if(TmpProbeLen < (uint32_t)(cMinInDelSeqLen * CurInDelMismatches))
				break;
			}
	    if(Idx != TmpProbeLen)
			continue;

		if(CurScore > CurDelete.Score) // explore right, deletion from probe...
			{
			memset(CurDelete.Seg,0,sizeof(CurDelete.Seg));
			CurDelete.Seg[0].MatchLen = ProbeMMOfss[MMIdx];
			CurDelete.Seg[0].MatchLoci =	TargOfs;		// target suffix, not chrom loci
			CurDelete.Seg[0].Mismatches = MMIdx;
			CurDelete.Seg[0].TrimMismatches = MMIdx;
			CurDelete.Seg[0].Strand = CurStrand;

			CurDelete.Seg[1].MatchLen = TmpProbeLen;
			CurDelete.Seg[1].MatchLoci =	TargOfs + CurDelete.Seg[0].MatchLen + CurInDelLen;		// target suffix, not chrom loci
			CurDelete.Seg[1].Mismatches = CurInDelMismatches;
			CurDelete.Seg[1].TrimMismatches = CurInDelMismatches;
			CurDelete.Seg[1].ReadOfs = CurDelete.Seg[0].MatchLen;
			CurDelete.Seg[1].Strand = CurStrand;
			CurDelete.Score = CurScore;
			CurDelete.FlgChimeric = 0;
			CurDelete.FlgInDel = 1;
			CurDelete.FlgSplice = 0;
			CurDelete.FlgInsert = 0;
			CurDelete.FlgNonOrphan = 0;
			}
		}
	}

if(CurDelete.Score == 0x00 && CurInsert.Score == 0x00)	// at least one will be > 0 if microInDel located
	return(0);

// have accepted this microInDel
if(CurDelete.Score > CurInsert.Score)
	{
	*pHit = CurDelete;
	return(3);
	}
*pHit = CurInsert;
return(2);
}

// ExploreInDelMatchLeft
// Normally this method is only invoked if unable to align using LocateCoreMultiples, e.g. exact or near exact matches, which
// would not have aligned reads in which there were micro InDels (small InDels in range of 1..microInDelLen, currently defined max is 20)
// micro InDels are only explored if there are more than cMaxMMExploreInDel (currently defined as 7) mismatches between the read and any other subsequence
// of same length in the targeted genome. The flanking length also needs to be at least cMinInDelSeqLen (currently defined as 7).
// Insertions (or deletions from the targeted sequence) into the read sequence are first scored, then deletions (insertions into target) from
// the read are scored with the lowest scored InDel chosen. Currently scoring is 2 for a mismatch, 10 for starting an InDel and 1 for InDel extension per base extended.
//
int									// -1 if errors, 0 if no match, 1 if mismatches only, 2 InDel and insertion into probe, 3 InDel and deletion from probe
CSfxArray::ExploreInDelMatchLeft(uint32_t ExtdProcFlags,			// flags indicating if lower levels need to do any form of extended processing with this specific read...
						 uint32_t ReadID,		// identifies this read
						 char CurStrand,	// aligning on this strand
			int microInDelLen,		   		// microInDel can be upto (inclusive) this length
			int MaxTotMM,			// can be at most this many mismatches in total
 		   int ProbeLen,			// length of probe excluding any eBaseEOS
		   etSeqBase *pProbe,		// pts to 5' start of probe sequence
		   tsSfxEntry *pEntry,		// target at TargSeqLeftIdx is in this sfx entry
		   int64_t TargOfs,			// pTarg corresponds to this suffix offset
		   etSeqBase *pTarg,		// pts to target sequence
		   tsHitLoci *pHit)			// where to return hit loci
{
int CurScore;						// current score for putative InDel being explored

tsHitLoci CurInsert;				// current insert loci
tsHitLoci CurDelete;				// current delete loci

int Idx;
etSeqBase *pP,*pT;
etSeqBase ProbeBase,TargBase;
int CurNumMismatches;

uint32_t TmpProbeLen;
int CurInDelMismatches;
int CurInDelLen;

uint32_t TargSeqLen;

int ProbeMMOfss[cMaxPutInDelOfss];			// offsets of all mismatches in probe which may need to be explored as putative InDel starts
int MMIdx;									// current index into ProbeMMOfss
int TotMM;									// total mismatches to be explored as putative InDels

memset(pHit,0,sizeof(tsHitLoci));
if(TargOfs < (int64_t)pEntry->StartOfs || ((TargOfs + ProbeLen - 1) > (int64_t)pEntry->EndOfs))
	return(0);


TargSeqLen = (uint32_t)((TargOfs - pEntry->StartOfs) + ProbeLen);

memset(&CurInsert,0,sizeof(tsHitLoci));
memset(&CurDelete,0,sizeof(tsHitLoci));

if(MaxTotMM > cMaxPutInDelOfss)		// silently clamp so as to ensure ProbeMMOfss never overflowed
	MaxTotMM = cMaxPutInDelOfss;

					// if exploring left then need to set pProbe to pt at 3' end of probe
pP = pProbe + ProbeLen - 1;
pT = pTarg  + ProbeLen - 1;

// count mismatches and record where these are, if too many mismatches then will explore for an InDel starting at each mismatch
CurNumMismatches = 0;
for(Idx = (ProbeLen-1); Idx >= 0 && CurNumMismatches <= max(MaxTotMM,cMaxMMExploreInDel); Idx--,pT-=1,pP-=1)
	{
	ProbeBase = *pP & 0x07;
	TargBase = *pT & 0x07;
	if(TargBase > eBaseN || ProbeBase > eBaseN)	// treat as no match any EOS
		return(0);

	if(ProbeBase == TargBase && ProbeBase <= eBaseT)
		continue;

		// have a mismatch
	ProbeMMOfss[CurNumMismatches] = Idx;		// record probe offset at which this mismatch occurred as this is where putative InDels may be subsequently explored
	CurNumMismatches += 1;
	}



if(CurNumMismatches < cMaxMMExploreInDel ||		    // only bother exploring putative InDels if too many mismatches
	cMinInDelSeqLen	> ProbeMMOfss[0])	// and the length of sequence which could be InDel'd is worth the effort
	{
	if(CurNumMismatches > MaxTotMM)
		return(0);
	memset(pHit->Seg,0,sizeof(pHit->Seg));
	pHit->Seg[0].MatchLen = ProbeLen;
	pHit->Seg[0].MatchLoci = TargOfs;
	pHit->Seg[0].Mismatches = CurNumMismatches;
	pHit->Seg[0].TrimMismatches = CurNumMismatches;
	pHit->Seg[0].Strand = CurStrand;
	pHit->Score = cBaseScore + (ProbeLen * cScoreMatch) - (CurNumMismatches * cScoreMismatch);
	return(1);
	}

// too many mismatches so now explore for InDels as insertions into probe starting from each mismatch
TotMM = min(MaxTotMM,CurNumMismatches);
CurScore = 0;
for(MMIdx = 0; MMIdx <= TotMM && cMinInDelSeqLen	< ProbeMMOfss[MMIdx]; MMIdx++)
	{
	// microInDels are to be limited in length
	for(CurInDelLen = 1; CurInDelLen <= microInDelLen; CurInDelLen++)
		{
		CurScore = cBaseScore + (ProbeLen * cScoreMatch) - (((CurInDelLen - 1) * cScoreInDelExtn) + cScoreInDelOpn);
		CurScore -= (MMIdx * cScoreMismatch);
		if(CurScore < CurInsert.Score)
			break;

		pP = &pProbe[ProbeMMOfss[MMIdx] - CurInDelLen];
		pT = &pTarg[ProbeMMOfss[MMIdx]];

		TmpProbeLen = ProbeMMOfss[MMIdx] - (CurInDelLen-1);

		if(TmpProbeLen < cMinInDelSeqLen)
			break;

		CurInDelMismatches = 0;
		for(Idx = 0; Idx < (int)TmpProbeLen &&
				 (MMIdx + CurInDelMismatches) <= MaxTotMM;
					  Idx++,pT-=1,pP-=1)
			{
			ProbeBase = *pP & 0x07;
			TargBase = *pT & 0x07;


			if(ProbeBase > eBaseN || TargBase > eBaseN)
				break;

			if(ProbeBase == TargBase && ProbeBase <= eBaseT)
				continue;

			CurInDelMismatches += 1;
			CurScore -= cScoreMismatch;

			if(TmpProbeLen < (uint32_t)(cMinInDelSeqLen * CurInDelMismatches))
				break;
			}

	    if(Idx != TmpProbeLen)
			continue;


		if(CurScore > CurInsert.Score) // explore left, insertion into probe...
			{
			memset(CurInsert.Seg,0,sizeof(CurInsert.Seg));
			CurInsert.Seg[0].MatchLen = TmpProbeLen;
			CurInsert.Seg[0].MatchLoci =	(uint32_t)(TargOfs + CurInDelLen);		// target suffix, not chrom loci
			CurInsert.Seg[0].Mismatches = CurInDelMismatches;
			CurInsert.Seg[0].TrimMismatches = CurInDelMismatches;
			CurInsert.Seg[0].Strand = CurStrand;

			CurInsert.Seg[1].MatchLen = ProbeLen - (TmpProbeLen + CurInDelLen);
			CurInsert.Seg[1].MatchLoci =	CurInsert.Seg[0].MatchLoci  + CurInsert.Seg[0].MatchLen;		// target suffix, not chrom loci
			CurInsert.Seg[1].Mismatches = MMIdx;
			CurInsert.Seg[1].ReadOfs = CurInsert.Seg[0].MatchLen + CurInDelLen;
			CurInsert.Seg[1].Strand = CurStrand;
			CurInsert.Score = CurScore;
			CurInsert.FlgChimeric = 0;
			CurInsert.FlgInDel = 1;
			CurInsert.FlgInsert = 1;
			CurInsert.FlgSplice = 0;
			CurInsert.FlgNonOrphan = 0;
			}
		}
	}

// now know if there was a putative insertion, check if there is a putative deletion and if so then choose the lowest scoring InDel as the one to report
CurScore = 0;
for(MMIdx = 0; MMIdx <= TotMM && cMinInDelSeqLen	< ProbeMMOfss[MMIdx]; MMIdx++)
	{
	for(CurInDelLen = 1; CurInDelLen <= microInDelLen; CurInDelLen++)
		{
		CurScore = cBaseScore + (ProbeLen * cScoreMatch) - (((CurInDelLen - 1) * cScoreInDelExtn) + cScoreInDelOpn);
		CurScore -= (MMIdx * cScoreMismatch);
		if(CurScore < CurDelete.Score)
			break;

		pP = &pProbe[ProbeMMOfss[MMIdx]];
		pT = &pTarg[ProbeMMOfss[MMIdx]-CurInDelLen];

		TmpProbeLen = ProbeMMOfss[MMIdx] + 1;

		if(TmpProbeLen < cMinInDelSeqLen)
			break;
		if((uint32_t)CurInDelLen > TargOfs)
			break;

		CurInDelMismatches = 0;
		for(Idx = 0; Idx < (int)TmpProbeLen && (MMIdx + CurInDelMismatches) <= MaxTotMM; Idx++,pT-=1,pP-=1)
			{
			ProbeBase = *pP & 0x07;
			TargBase = *pT & 0x07;
			if(ProbeBase > eBaseN || TargBase > eBaseN)
				break;

			if(ProbeBase == TargBase && ProbeBase <= eBaseT)
				continue;

			CurInDelMismatches += 1;
			CurScore -= cScoreMismatch;
			if(TmpProbeLen < (uint32_t)(cMinInDelSeqLen * CurInDelMismatches))
				break;
			}
	    if(Idx != TmpProbeLen)
			continue;


		if(CurScore > CurDelete.Score) // explore left, deletion from probe...
			{
			memset(CurDelete.Seg,0,sizeof(CurDelete.Seg));
			CurDelete.Seg[0].MatchLen = TmpProbeLen;
			CurDelete.Seg[0].MatchLoci =	(uint32_t)(TargOfs - CurInDelLen);		// target suffix, not chrom loci
			CurDelete.Seg[0].Mismatches = CurInDelMismatches;
			CurDelete.Seg[0].Strand = CurStrand;

			CurDelete.Seg[1].MatchLen = ProbeLen - TmpProbeLen;
			CurDelete.Seg[1].MatchLoci = CurDelete.Seg[0].MatchLoci + TmpProbeLen + CurInDelLen; // target suffix, not chrom loci
			CurDelete.Seg[1].Mismatches = MMIdx;
			CurDelete.Seg[1].ReadOfs = TmpProbeLen;
			CurDelete.Seg[1].Strand = CurStrand;
			CurDelete.Score = CurScore;
			CurDelete.FlgChimeric = 0;
			CurDelete.FlgInDel = 1;
			CurDelete.FlgInsert = 0;
			CurDelete.FlgSplice = 0;
			CurDelete.FlgNonOrphan = 0;
			}
		}
	}

if(CurDelete.Score == 0x00 && CurInsert.Score == 0x00)	// at least one will be > 0 if microInDel located
	return(0);

// have accepted this microInDel
if(CurDelete.Score > CurInsert.Score)
	{
	*pHit = CurDelete;
	return(3);
	}
*pHit = CurInsert;
return(2);
}


int
CSfxArray::QSortSeq(int64_t SeqLen,		// total concatenated sequence length
						etSeqBase *pSeq,	// pts to start of concatenated sequences
						int SfxElSize,		// suffix element size (will be either 4 or 5)
						void *pArray)		// allocated to hold suffix elements
{
int64_t Idx;
gpSeq = pSeq;

switch(SfxElSize) {
	case 4:	// 4 bytes per element
			{
			uint32_t *pIdx = (uint32_t *)pArray;
			gpSfxArray = (uint8_t *)pIdx;
			for(Idx = 0; Idx < SeqLen; Idx++, pIdx++)
				*pIdx = (uint32_t)Idx;
			m_MTqsort.qsort(gpSfxArray,SeqLen,4,QSortSeqCmp32);
			}
		break;
	case 5:		// 5 bytes per element
			{
			uint8_t *pIdx = (uint8_t *)pArray;
			gpSfxArray = pIdx;
			gpSeq = pSeq;
			for(Idx = 0; Idx < SeqLen; Idx++, pIdx++)
				{
				*(uint32_t *)pIdx = Idx & 0x0ffffffff;
				pIdx += 4;
				*pIdx = (Idx >> 32) & 0x00ff;
				}
			m_MTqsort.qsort(gpSfxArray,SeqLen,5,QSortSeqCmp40);
			}
		break;
	default:			// any other element size is unsupported
		return(-1);
	}
return(0);
}

// QSortSeqCmp32
// qsorts suffix elements whereby each element occupies 32bits, 4 bytes, and is an offset into gpSeq[]
static int QSortSeqCmp32(const void *p1,const void *p2)
{
etSeqBase *pSeq1;
etSeqBase *pSeq2;
etSeqBase b1;
etSeqBase b2;
int MaxCmpLen;

pSeq1 = &gpSeq[*(uint32_t *)p1];
pSeq2 = &gpSeq[*(uint32_t *)p2];

// compare seqs for at most gMaxBaseCmpLen bases
MaxCmpLen = gMaxBaseCmpLen + 10;		// allow an additional 10bp extension as a small safety margin in case user miscalculated comparison K-mer sizes!
while(MaxCmpLen--)
	{
	if((b1 = (*pSeq1++ & 0x0f)) != (b2 = (*pSeq2++ & 0x0f)))
		return(b1 < b2 ? -1 : 1);
	if (b1 == eBaseEOS || b2 == eBaseEOS)
		break;
	}
return(0);
}

// QSortSeqCmp40
// qsorts suffix elements whereby each element occupies 40bits, 5 bytes, and is an offset into gpSeq[]
static int QSortSeqCmp40(const void *p1,const void *p2)
{
etSeqBase *pSeq1;
etSeqBase *pSeq2;
uint8_t *pP;
int64_t Ofs1;
int64_t Ofs2;
etSeqBase b1;
etSeqBase b2;
int MaxCmpLen;
pP = (uint8_t *)p1;
Ofs1 = (int64_t)*(uint32_t *)pP;
Ofs1 |= ((int64_t)pP[4] << 32);
pP = (uint8_t *)p2;
Ofs2 = (int64_t)*(uint32_t *)pP;
Ofs2 |= ((int64_t)pP[4] << 32);

pSeq1 = &gpSeq[Ofs1];
pSeq2 = &gpSeq[Ofs2];

// compare seqs for at most gMaxBaseCmpLen bases
MaxCmpLen = gMaxBaseCmpLen + 10;		// allow an additional 10bp extension as a small safety margin in case user miscalculated comparison K-mer sizes!
while(MaxCmpLen--)
	{
	if((b1 = (*pSeq1++ & 0x0f)) != (b2 = (*pSeq2++ & 0x0f)))
		return(b1 < b2 ? -1 : 1);
	if(b1 == eBaseEOS || b2 == eBaseEOS)
		break;
	}
return(0);
}


static int QSortEntryNames(const void *p1,const void *p2)
{
tsSfxEntry *pE1 = *(tsSfxEntry **)p1;
tsSfxEntry *pE2 = *(tsSfxEntry **)p2;

return(stricmp((char *)pE1->szSeqName,(char *)pE2->szSeqName));
}
