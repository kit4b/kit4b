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
#include <process.h>
#include "./commhdrs.h"
#else
#include <sys/mman.h>
#include <pthread.h>
#include "./commhdrs.h"
#endif

#include "./CPBASfxArray.h"

static int gMaxBaseCmpLen = cMaxReadLen;   // used to limit the number of bases being compared for length in QSortSeqCmp32() and QSortSeqCmp40()
static int CPBAQSortSeqCmp32(const void* p1, const void* p2);
static int CPBAQSortSeqCmp40(const void* p1, const void* p2);
static int CPBAQSortEntryNames(const void* p1, const void* p2);

static uint8_t* gpPBASfxArray = NULL;
static uint8_t* gpPBASeq = NULL;

CPBASfxArray::CPBASfxArray()
{
m_pSfxSeqs = nullptr;
m_pEntries = nullptr;
m_pSfxIdx = nullptr;
m_AllocSfxIdxMem = 0;
m_AllocEntriesMem = 0;
m_AllocSfxSeqsMem = 0;
m_CurUsedSfxSeqsMem = 0;
m_CurNumEntries = 0;
m_AllocEntries = 0;
m_MaxQSortThreads = 16;
gMaxBaseCmpLen = cMaxReadLen;
m_SfxElSize = 4;
SetMaxQSortThreads(cDfltSortThreads);
}

CPBASfxArray::~CPBASfxArray()
{
if (m_pSfxSeqs != NULL)
	{
#ifdef _WIN32
	free(m_pSfxSeqs);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (m_pSfxSeqs != MAP_FAILED)
		munmap(m_pSfxSeqs, m_AllocSfxSeqsMem);
#endif
	}

if (m_pEntries != NULL)
	{
#ifdef _WIN32
	free(m_pEntries);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (m_pEntries != MAP_FAILED)
		munmap(m_pEntries, m_AllocEntriesMem);
#endif
	}

if (m_pSfxIdx != NULL)
	{
#ifdef _WIN32
	free(m_pSfxIdx);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (m_pSfxIdx != MAP_FAILED)
		munmap(m_pSfxIdx, m_AllocSfxIdxMem);
#endif
	}
}

void 
CPBASfxArray::Reset(void)
{
if (m_pSfxSeqs != NULL)
	{
#ifdef _WIN32
	free(m_pSfxSeqs);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (m_pSfxSeqs != MAP_FAILED)
		munmap(m_pSfxSeqs, m_AllocSfxSeqsMem);
#endif
	m_pSfxSeqs = nullptr;
	}
m_AllocSfxSeqsMem = 0;

if (m_pEntries != NULL)
	{
#ifdef _WIN32
	free(m_pEntries);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (m_pEntries != MAP_FAILED)
		munmap(m_pEntries, m_AllocEntriesMem);
#endif
	m_pEntries = nullptr;
	}
m_AllocEntriesMem = 0;

if (m_pSfxIdx != NULL)
	{
#ifdef _WIN32
	free(m_pSfxIdx);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (m_pSfxIdx != MAP_FAILED)
		munmap(m_pSfxIdx, m_AllocSfxIdxMem);
#endif
	m_pSfxIdx = nullptr;
	}
m_AllocSfxIdxMem = 0;

m_CurUsedSfxSeqsMem = 0;
m_CurNumEntries = 0;
m_AllocEntries = 0;
gMaxBaseCmpLen = cMaxReadLen;

m_SfxElSize = 4;
SetMaxQSortThreads(cDfltSortThreads);
}

// Suffix array elements can be sized as either 4 or 5 bytes dependent on the total length of concatenated sequences
// If total length is less than 4G then can use 4 byte elements, if longer then will use 5 byte elements
inline
int64_t SfxOfsToLoci(int SfxElSize,	// size in bytes of suffix element - expected to be either 4 or 5
	void* pSfx,						// pts to 1st element of suffix array
	int64_t Ofs)					// offset to suffix element
{
uint64_t Loci;
uint8_t* pSfxEls = (uint8_t*)pSfx;
Ofs *= SfxElSize;
Loci = (uint64_t) * (uint32_t*)&pSfxEls[Ofs];
if (SfxElSize == 5)
	Loci |= ((uint64_t)pSfxEls[Ofs + 4]) << 32;
return(Loci);
}

// AddEntry
// Adds new entry and it's associated sequence
teBSFrsltCodes
CPBASfxArray::AddEntry(uint32_t UserID,		// user supplied identifier
	uint32_t SeqID,							// sequence id, typically related to the textual pszSeqIdent
	char* pszSeqIdent,						// sequence identifier, typically chromosome name
	uint8_t* pSeq,							// PBA sequence to add to suffix array
	int32_t SeqLen)							// PBA sequence length
{
uint8_t* pTmpSeq;
tsPBASfxEntry* pTmpEntry;
size_t ReallocSize;

if (SeqLen < cMinAllowInMemPBASeqLen || SeqLen > cMaxAllowInMemPBASeqLen)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddEntry: SeqLen %u not in range 1..%u", SeqLen, cMaxAllowInMemPBASeqLen);
	Reset();
	return(eBSFerrParams);
	}

if (m_pEntries == NULL)					// will be NULL until at least one entry has been added
	{
	if (m_pSfxSeqs != NULL)				// if no entries then shouldn't be any sequences, but better to be sure!
		{
#ifdef _WIN32
		free(m_pSfxSeqs);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
		if (m_pSfxSeqs != MAP_FAILED)
			munmap(m_pSfxSeqs, m_AllocSfxSeqsMem);
#endif
		m_pSfxSeqs = NULL;
		}
	m_AllocSfxSeqsMem = 0;
	m_CurUsedSfxSeqsMem = 0;

		// initially allocate for default number of entries
	m_CurNumEntries = 0;
	m_AllocEntries = 0;
	m_AllocEntriesMem = 0;
	ReallocSize = sizeof(tsPBASfxEntry) * cAllocPBASfxEntries;
#ifdef _WIN32
	m_pEntries = (tsPBASfxEntry *)malloc(ReallocSize);
	if (m_pEntries == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddEntry: unable to allocate %zu bytes for holding entries", ReallocSize);
		Reset();
		return(eBSFerrMem);
		}
#else
		// gnu malloc is still in the 32bit world and seems to have issues if more than 2GB allocation
	m_pEntries = (tsPBASfxEntry*)mmap(NULL, ReallocSize, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	if (m_pEntries == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddEntry: unable to allocate %zd bytes for holding entries", ReallocSize);
		m_pEntries = NULL;
		Reset();
		return(eBSFerrMem);
		}
#endif
	memset(m_pEntries, 0, ReallocSize);
	m_AllocEntriesMem = (int64_t)ReallocSize;
	m_AllocEntries = cAllocPBASfxEntries;
	}

if(m_pSfxSeqs == nullptr)
	{
	ReallocSize = (size_t)cAllocSfxSeqs;

#ifdef _WIN32
	m_pSfxSeqs = (uint8_t *)malloc((size_t)ReallocSize);
	if (m_pSfxSeqs == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddEntry: unable to allocate suffix memory %zd", ReallocSize);
		Reset();
		return(eBSFerrMem);
		}
#else
		// gnu malloc is still in the 32bit world and seems to have issues if more than 2GB allocation
	m_pSfxSeqs = (uint8_t *)mmap(NULL, ReallocSize, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	if (m_pSfxSeqs == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddEntry: unable to allocate suffix block memory %zd", ReallocSize);
		m_pSfxSeqs = NULL;
		Reset();
		return(eBSFerrMem);
		}
#endif
	m_AllocSfxSeqsMem = ReallocSize;
	m_CurUsedSfxSeqsMem = 0;
	memset(m_pSfxSeqs, 0, ReallocSize);
	}


// ensure not about to exceed max number of allowed entries
if ((m_CurNumEntries+1) > cMaxPBASfxEntries)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddEntry: Reached max number (%d) of supported entries", cMaxPBASfxEntries);
	Reset();
	return(eBSFerrMaxEntries);
	}

	// check if adding new entry sequence would cause the currently used m_CurUsedSfxSeqsMem to overflow cMaxAllowConcatPBASeqLen
if ((m_CurUsedSfxSeqsMem + SeqLen + 10) > cMaxAllowConcatPBASeqLen)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddEntry: Total concatenated sequence length (%zd) is more than maximum (%zd) supported", m_CurUsedSfxSeqsMem + SeqLen + 10, cMaxAllowConcatPBASeqLen);
	Reset();
	return(eBSFerrMem);
	}

	// check if entries needs to be extended
if ((m_CurNumEntries+1) >= m_AllocEntries)
	{
	int32_t ReallocEntriesTo;
	tsPBASfxEntry* pRealloc;
	ReallocEntriesTo = min(cMaxPBASfxEntries, m_CurNumEntries + cAllocPBASfxEntries + 1);
	ReallocSize = (sizeof(tsPBASfxEntry) * ReallocEntriesTo);
#ifdef _WIN32
	pRealloc = (tsPBASfxEntry*)realloc(m_pEntries, ReallocSize);
#else
	pRealloc = (tsPBASfxEntry*)mremap(m_pEntries, (size_t)m_AllocEntriesMem, ReallocSize, MREMAP_MAYMOVE);
	if (pRealloc == MAP_FAILED)
		pRealloc = NULL;
#endif
	if (pRealloc == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddEntry: Memory reallocation from %zd bytes to %zd  - %s", m_AllocEntriesMem, (int64_t)ReallocSize, strerror(errno));
		Reset();
		return(eBSFerrMem);
		}

	m_pEntries = pRealloc;
	m_AllocEntries = ReallocEntriesTo;
	m_AllocEntriesMem = (int64_t)ReallocSize;
	}

		// check if m_pSfxSeqs needs to be extended
if ((m_CurUsedSfxSeqsMem + SeqLen + 10) >= m_AllocSfxSeqsMem) // 10 is to allow for appended sequence separators (0xff)
	{
	uint8_t* pRealloc;
	ReallocSize = (size_t)min(cMaxAllowConcatPBASeqLen, cReallocSfxSeqs + m_AllocSfxSeqsMem + SeqLen + 10);

#ifdef _WIN32
	pRealloc = (uint8_t*)realloc(m_pSfxSeqs, ReallocSize);
#else
	pRealloc = (uint8_t*)mremap(m_pSfxSeqs, (size_t)m_AllocSfxSeqsMem, ReallocSize, MREMAP_MAYMOVE);
	if (pRealloc == MAP_FAILED)
		pRealloc = NULL;
#endif
	if (pRealloc == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddEntry: SfxSeq memory reallocation from %zd to %zd bytes - %s", m_AllocSfxSeqsMem, (int64_t)ReallocSize, strerror(errno));
		Reset();
		return(eBSFerrMem);
		}

	m_pSfxSeqs = pRealloc;
	m_AllocSfxSeqsMem = ReallocSize;
	}

pTmpEntry = &m_pEntries[m_CurNumEntries++];
pTmpEntry->EntryID = m_CurNumEntries;
pTmpEntry->UserID = UserID;
pTmpEntry->SeqID = SeqID;
strncpy((char*)pTmpEntry->szSeqName, pszSeqIdent, sizeof(pTmpEntry->szSeqName));
pTmpEntry->szSeqName[sizeof(pTmpEntry->szSeqName) - 1] = '\0';
pTmpEntry->NameHash = CUtility::GenHash16(pszSeqIdent);
pTmpEntry->SeqLen = SeqLen;
pTmpEntry->StartOfs = m_CurUsedSfxSeqsMem;
pTmpEntry->EndOfs = m_CurUsedSfxSeqsMem + SeqLen - 1;
pTmpSeq = &m_pSfxSeqs[m_CurUsedSfxSeqsMem];
memcpy(pTmpSeq, pSeq, SeqLen);
pTmpSeq += SeqLen;
*pTmpSeq = 0xff;
m_CurUsedSfxSeqsMem += SeqLen + 1;
if(m_CurUsedSfxSeqsMem < cThres8BytePBASfxEls)			// if less than cThres8BytePBASfxEls then can use 4 bytes per suffix element otherwise need to go to 5 bytes
	m_SfxElSize = 4;
else
	m_SfxElSize = 5;
return(eBSFSuccess);
}

int
CPBASfxArray::SuffixSort(void)
{
int64_t ReallocSize;
if(m_CurUsedSfxSeqsMem == 0)
	return(0);

ReallocSize = (m_CurUsedSfxSeqsMem + 100) * m_SfxElSize;

if (m_pSfxIdx == nullptr)
	{
#ifdef _WIN32
	m_pSfxIdx = (void *)malloc((size_t)ReallocSize);
	if (m_pSfxIdx == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "SuffixSort: SfxIdx memory allocation to %zd bytes - %s", ReallocSize, strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and seems to have issues if more than 2GB allocation
	m_pSfxIdx = (void *)mmap(NULL, (size_t)ReallocSize, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	if (m_pSfxIdx == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "SuffixSort: SfxIdx memory allocation to %zd bytes - %s", ReallocSize, strerror(errno));
		m_pSfxIdx = NULL;
		Reset();
		return(eBSFerrMem);
		}
#endif
	m_AllocSfxIdxMem = ReallocSize;
	}
else
	if (ReallocSize > m_AllocSfxIdxMem)
		{
		void* pRealloc;
#ifdef _WIN32
		pRealloc = (void*)realloc(m_pSfxIdx,(size_t)ReallocSize);
#else
		pRealloc = (void*)mremap(m_pSfxIdx, (size_t)m_AllocSfxIdxMem, (size_t)ReallocSize, MREMAP_MAYMOVE);
		if (pRealloc == MAP_FAILED)
			pRealloc = NULL;
#endif
		if (pRealloc == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "SuffixSort: SfxIdx memory re-allocation from %zd to %zd bytes - %s", m_AllocSfxIdxMem, ReallocSize, strerror(errno));
			return(eBSFerrMem);
			}

		m_pSfxIdx = pRealloc;
		m_AllocSfxIdxMem = ReallocSize;
		}
QSortSeq(m_CurUsedSfxSeqsMem, m_pSfxSeqs, m_SfxElSize, m_pSfxIdx);
return(eBSFSuccess);
}



void
CPBASfxArray::SetMaxQSortThreads(int MaxThreads)			// sets maximum number of threads to use in multithreaded qsorts
{
m_MaxQSortThreads = MaxThreads;
m_MTqsort.SetMaxThreads(MaxThreads);
}

int						// returns the previously utilised MaxBaseCmpLen
CPBASfxArray::SetMaxBaseCmpLen(int MaxBaseCmpLen)		// sets maximum number of bases which need to be compared for equality in multithreaded qsorts, will be clamped to be in range 10..cMaxReadLen
{
	int PrevMaxBaseCmpLen = gMaxBaseCmpLen;
	if (MaxBaseCmpLen > cMaxReadLen)
		MaxBaseCmpLen = cMaxReadLen;
	else
		if (MaxBaseCmpLen < 10)
			MaxBaseCmpLen = 10;
	gMaxBaseCmpLen = MaxBaseCmpLen;
	return(PrevMaxBaseCmpLen);
}

// Compares probe against target taking into account any 0xff terminator
int
CPBASfxArray::CmpProbeTarg(uint8_t* pEl1, uint8_t* pEl2, int Len)
{
int Ofs;
uint8_t El1;
uint8_t El2;
for (Ofs = 0; Ofs < Len; Ofs++)
	{
	El2 = *pEl2++;
	if (El2 == 0xff)
		return(-1);
	El1 = *pEl1++;
	if (El1 > El2)
		return(1);
	if (El1 < El2)
		return(-1);
	}
return(0);
}

// MapHit2Entry
// Maps the hit offset to the relevant sequence entry
// If many entries expected then this mapping would be a good candidate for optimisation!
tsPBASfxEntry*
CPBASfxArray::MapHit2Entry(int64_t Ofs)
{
tsPBASfxEntry* pProbe;
int64_t Lo, Mid, Hi;	// search limits

Lo = 0;
Hi = m_CurNumEntries - 1;
while (Hi >= Lo) {
	Mid = (Hi + Lo) / 2;
	pProbe = &m_pEntries[Mid];

	if (pProbe->StartOfs <= Ofs && pProbe->EndOfs >= Ofs)
		return(pProbe);

	if (pProbe->StartOfs > Ofs)
		{
		Hi = Mid - 1;
		continue;
		}
	else
		if (pProbe->EndOfs < Ofs)
			{
			Lo = Mid + 1;
			continue;
			}
	}
return(NULL);
}

int64_t			// returned index (1..n) into m_pSfxIdx[] of first suffix sequence exactly matching pProbe, 0 if none matching
CPBASfxArray::LocateFirstExact(uint8_t *pProbe,  // pts to probe sequence
	int ProbeLen,					// probe length to exactly match over
	int SfxElSize,					// size in bytes of suffix element - expected to be either 4 or 5
	int64_t TargStart,				// position in pTarg (0..n) corresponding to start of suffix array
	int64_t SfxLo,					// low index in pSfxArray
	int64_t SfxHi)					// high index in pSfxArray
{
uint8_t* pEl1;
uint8_t* pEl2;
uint8_t El1;
uint8_t El2;

int CmpRslt;
int Ofs;
int64_t Mark;
int64_t TargPsn;

do {
	pEl1 = pProbe;
	TargPsn = ((int64_t)SfxLo + SfxHi) / 2L;
	pEl2 = &m_pSfxSeqs[SfxOfsToLoci(SfxElSize, m_pSfxIdx, TargPsn + TargStart)];
	CmpRslt = 0;
	for (Ofs = 0; Ofs < ProbeLen; Ofs++)
		{
		El2 = *pEl2++;
		if (El2 == 0xff)
			{
			CmpRslt = -1;
			break;
			}
		El1 = *pEl1++;
		if (El1 > El2)
			{
			CmpRslt = 1;
			break;
			}
		if (El1 < El2)
			{
			CmpRslt = -1;
			break;
			}
		}

	if (!CmpRslt)	// if a match then may not be the lowest indexed match
		{
		if (TargPsn == 0 || SfxLo == TargPsn) // check if already lowest
			return(TargPsn + 1);
		// iterate until lowest located
		while (1) {
			if (CmpRslt == 0)
				{
				Mark = TargPsn;
				if (Mark == 0)
					return(Mark + 1);
				SfxHi = TargPsn - 1;
				}
			TargPsn = ((int64_t)SfxLo + SfxHi) / 2L;
			pEl2 = &m_pSfxSeqs[SfxOfsToLoci(SfxElSize, m_pSfxIdx, TargPsn + TargStart)];

			pEl1 = pProbe;
			CmpRslt = 0;
			for (Ofs = 0; Ofs < ProbeLen; Ofs++)
				{
				El2 = *pEl2++;
				if (El2 == 0xff)
					{
					CmpRslt = -1;
					break;
					}
				El1 = *pEl1++;
				if (El1 > El2)
					{
					CmpRslt = 1;
					break;
					}
				if (El1 < El2)
					{
					CmpRslt = -1;
					break;
					}
				}
				
			if (CmpRslt == 0)				// 0 if still matching
				continue;
			SfxLo = TargPsn + 1;
			if (SfxLo == Mark)
				return(Mark + 1);
		}
	}

	if (CmpRslt < 0)
		{
		if (TargPsn == 0)
			break;
		SfxHi = TargPsn - 1;
		}
	else
		SfxLo = TargPsn + 1;
	} 
while (SfxHi >= SfxLo);
return(0);	// unable to locate any instance of pProbe
}



int64_t		// returned index (1..n) into m_pSfxIdx[] of first suffix sequence exactly matching pProbe, 0 if none matching
CPBASfxArray::LocateExact(uint8_t* pProbe,  // pts to probe sequence
			int ProbeLen)					// probe length to exactly match over
{
return(LocateFirstExact(pProbe,ProbeLen,m_SfxElSize,0,0, m_CurUsedSfxSeqsMem - 1));
}


bool		// returns true if the probe sequence exactly matches the sequence in claimed entry starting at claimed loci, otherwise false
CPBASfxArray::IsValidatedClaimedExact(uint8_t *pProbeSeq, 
						int32_t ProbeLen,		// probe length
						tsPBASfxEntry* pEntry,	// claimed match was on this entry sequence 
						uint32_t HitLoci)		// starting at this loci in entry sequence
{
uint8_t* pClaimSeq;
int32_t Idx;

pClaimSeq = &m_pSfxSeqs[pEntry->StartOfs + HitLoci];

for (Idx = 0; Idx < ProbeLen - 1; Idx++, pProbeSeq++, pClaimSeq++)
	{
	if (*pProbeSeq == 0xff || *pClaimSeq == 0xff)
		return(false);
	if (*pProbeSeq != *pClaimSeq)
		return(false);
	}
return(true);
}

int64_t						// returned index (1..n) into m_pSfxIdx[] of first suffix sequence exactly matching pProbe, 0 if none matching
CPBASfxArray::IterateExacts(uint8_t* pProbeSeq,// probe
				int32_t ProbeLen,		// probe length
				int64_t PrevHitIdx,		// 0 if starting new sequence, otherwise set to return value of previous successful iteration return
				uint32_t* pUserID,		// if match then where to return UserID in matched entry
				uint32_t* pSeqID,		// if match then where to return SeqID in matched entry
				uint32_t* pHitLoci)		// if match then where to return loci
{
int Cmp;
int64_t TargPsn;
tsPBASfxEntry* pEntry;
uint8_t* pEl1;
uint8_t* pEl2;

uint8_t* pTarg;			// target sequence
void* pSfxArray;		// target sequence suffix array
int64_t SfxLen;			// number of suffixs in pSfxArray

if(pUserID != nullptr)
	*pUserID = 0;
if (pSeqID != nullptr)
	*pSeqID = 0;
if (pHitLoci != nullptr)
	*pHitLoci = 0;

	// ensure prev hit was not the last!
if (PrevHitIdx >= m_CurUsedSfxSeqsMem)
	return(0);

pTarg = m_pSfxSeqs;
pSfxArray = m_pSfxIdx;
SfxLen = m_CurUsedSfxSeqsMem;

if (!PrevHitIdx)
	{
		// locate first exact match
	if ((TargPsn = LocateFirstExact(pProbeSeq, ProbeLen, m_SfxElSize, 0, 0, m_CurUsedSfxSeqsMem - 1)) == 0)
		return(0);	// no match
	TargPsn -= 1;
	pEntry = MapHit2Entry(SfxOfsToLoci(m_SfxElSize, pSfxArray, TargPsn));
	if(pUserID != nullptr)
		*pUserID = pEntry->UserID;
	if (pSeqID != nullptr)
		*pSeqID = pEntry->SeqID;
	if (pHitLoci != nullptr)
		*pHitLoci = (uint32_t)(SfxOfsToLoci(m_SfxElSize, pSfxArray, TargPsn) - pEntry->StartOfs);
#ifdef _CHECKVALIDATEKMERSEQS
	bool bIsValidated = IsValidatedClaimedExact(pProbeSeq, ProbeLen, pEntry, (uint32_t)(SfxOfsToLoci(m_SfxElSize, pSfxArray, TargPsn) - pEntry->StartOfs));
#endif
	return(TargPsn + 1);
	}

	// check if probe matches next suffix
pEl2 = &pTarg[SfxOfsToLoci(m_SfxElSize, pSfxArray, PrevHitIdx)];
pEl1 = pProbeSeq;

Cmp = CmpProbeTarg(pEl1, pEl2, ProbeLen);
if (!Cmp)
	{
	pEntry = MapHit2Entry(SfxOfsToLoci(m_SfxElSize, pSfxArray, PrevHitIdx));
	if (pUserID != nullptr)
		*pUserID = pEntry->UserID;
	if (pSeqID != nullptr)
		*pSeqID = pEntry->SeqID;
	if(pHitLoci != nullptr)
		*pHitLoci = (uint32_t)(SfxOfsToLoci(m_SfxElSize, pSfxArray, PrevHitIdx) - pEntry->StartOfs);
#ifdef _CHECKVALIDATEKMERSEQS
	bool bIsValidated = IsValidatedClaimedExact(pProbeSeq, ProbeLen, pEntry, (uint32_t)(SfxOfsToLoci(m_SfxElSize, pSfxArray, PrevHitIdx) - pEntry->StartOfs));
#endif
	return(PrevHitIdx + 1);
	}

return(0);		// no more hits
}

int
CPBASfxArray::QSortSeq(int64_t SeqLen,		// total concatenated sequence length
	uint8_t *pSeq,	// pts to start of concatenated sequences
	int SfxElSize,		// suffix element size (will be either 4 or 5)
	void* pArray)		// allocated to hold suffix elements
{
int64_t Idx;
gpPBASeq = pSeq;

switch (SfxElSize) {
	case 4:	// 4 bytes per element
		{
		uint32_t* pIdx = (uint32_t*)pArray;
		gpPBASfxArray = (uint8_t*)pIdx;
		for (Idx = 0; Idx < SeqLen; Idx++, pIdx++)
			*pIdx = (uint32_t)Idx;
		m_MTqsort.qsort(gpPBASfxArray, SeqLen, 4, CPBAQSortSeqCmp32);
		}
		break;
	case 5:		// 5 bytes per element
		{
		uint8_t* pIdx = (uint8_t*)pArray;
		gpPBASfxArray = pIdx;
		gpPBASeq = pSeq;
		for (Idx = 0; Idx < SeqLen; Idx++, pIdx++)
			{
			*(uint32_t*)pIdx = Idx & 0x0ffffffff;
			pIdx += 4;
			*pIdx = (Idx >> 32) & 0x00ff;
			}
		m_MTqsort.qsort(gpPBASfxArray, SeqLen, 5, CPBAQSortSeqCmp40);
		}
		break;
	default:			// any other element size is unsupported
		return(-1);
	}
return(0);
}

// QSortSeqCmp32
// qsorts suffix elements whereby each element occupies 32bits, 4 bytes, and is an offset into gpSeq[]
static int CPBAQSortSeqCmp32(const void* p1, const void* p2)
{
	uint8_t* pSeq1;
	uint8_t* pSeq2;
	uint8_t b1;
	uint8_t b2;
	int MaxCmpLen;

	pSeq1 = &gpPBASeq[*(uint32_t*)p1];
	pSeq2 = &gpPBASeq[*(uint32_t*)p2];

	// compare seqs for at most gMaxBaseCmpLen bases
	MaxCmpLen = gMaxBaseCmpLen + 1;		// allow an additional 1bp extension as a small safety margin in case user miscalculated comparison K-mer sizes!
	while (MaxCmpLen--)
	{
		if ((b1 = *pSeq1++) != (b2 = *pSeq2++))
			return(b1 < b2 ? -1 : 1);
		if(b1 == 0xff || b2 == 0xff)
			break;
	}
	return(0);
}

// QSortSeqCmp40
// qsorts suffix elements whereby each element occupies 40bits, 5 bytes, and is an offset into gpSeq[]
static int CPBAQSortSeqCmp40(const void* p1, const void* p2)
{
	uint8_t* pSeq1;
	uint8_t* pSeq2;
	uint8_t* pP;
	int64_t Ofs1;
	int64_t Ofs2;
	uint8_t b1;
	uint8_t b2;
	int MaxCmpLen;
	pP = (uint8_t*)p1;
	Ofs1 = (int64_t) * (uint32_t*)pP;
	Ofs1 |= ((int64_t)pP[4] << 32);
	pP = (uint8_t*)p2;
	Ofs2 = (int64_t) * (uint32_t*)pP;
	Ofs2 |= ((int64_t)pP[4] << 32);

	pSeq1 = &gpPBASeq[Ofs1];
	pSeq2 = &gpPBASeq[Ofs2];

	// compare seqs for at most gMaxBaseCmpLen bases
	MaxCmpLen = gMaxBaseCmpLen + 1;		// allow an additional 1bp extension as a small safety margin in case user miscalculated comparison K-mer sizes!
	while (MaxCmpLen--)
		{
		if ((b1 = *pSeq1++) != (b2 = *pSeq2++))
			return(b1 < b2 ? -1 : 1);
		if (b1 == 0xff || b2 == 0xff)
			break;
		}
	return(0);
}


static int CPBAQSortEntryNames(const void* p1, const void* p2)
{
	tsPBASfxEntry* pE1 = *(tsPBASfxEntry**)p1;
	tsPBASfxEntry* pE2 = *(tsPBASfxEntry**)p2;

	return(stricmp((char*)pE1->szSeqName, (char*)pE2->szSeqName));
}
