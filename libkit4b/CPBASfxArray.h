#pragma once
#include "./commdefs.h"


const uint32_t cMinAllowInMemPBASeqLen = 0x04;		// minimum allowed length of any individual sequence allowed when constructing in-memory suffix array
const int32_t cMaxAllowInMemPBASeqLen = 0x7fffff00;		// maximum allowed length of any individual sequence allowed when constructing in-memory suffix array, slightly less than 2Gbp
const int64_t cMaxAllowConcatPBASeqLen = 1000000000000;		// max supported concatenation length of all sequences (must fit within 40bits)

const int64_t cAllocSfxSeqs = ((int64_t)cMaxAllowInMemPBASeqLen * 10);	// initial allocation of memory for m_pSfxSeqs in this sized chunk
const int64_t cReallocSfxSeqs = ((int64_t)cMaxAllowInMemPBASeqLen * 5);					// realloc memory for m_pSfxSeqs in this sized chunk

const int64_t cThres8BytePBASfxEls = 4000000000;			// if concatenated sequence length >= this threshold then use 5bytes per suffix element instead of 4 when creating suffix index


const int32_t cMaxPBASfxEntries = 1000000;				// supporting a maximum of this many entries
const int32_t cAllocPBASfxEntries = 10000;			// allocating/reallocating for this many entries

#pragma pack(1)
// each entry for sequences is described by the following fixed size structure
typedef struct TAG_sPBASfxEntry {
	int32_t EntryID;					// identifies each entry (1..n), unique within this suffix file
	int32_t UserID;						// user supplied identifier when entry was created with call to AddEntry()
	int32_t SeqID;						// user supplied sequence id, typically related to the textual pszSeqIdent
	uint8_t szSeqName[cMaxDatasetSpeciesChrom];	// entry name - typically a chromosome name
	uint16_t NameHash;					// hash over szSeqName
	int32_t SeqLen;						// sequence length - excludes any sequence terminator 0xff
	int64_t StartOfs;					// offset into concatenated sequences of where this sequence starts
	int64_t EndOfs;						// offset into concatenated sequences of where this sequence ends
} tsPBASfxEntry;
#pragma pack()


class CPBASfxArray : public CErrorCodes, protected CEndian
{

	int32_t m_SfxElSize;						// number of bytes per suffix array element = will be either 4 or 5

	int64_t m_CurUsedSfxSeqsMem;					// currently using this much of allocation memory size in m_pSfxSeqs (is also the number of elements in suffix array)
	int64_t m_AllocSfxSeqsMem;					// allocation memory size	
	uint8_t *m_pSfxSeqs;						// to hold concatenated PBA sequences

	int32_t m_CurNumEntries;					// current number of entries in m_pEntries[]
	int32_t m_AllocEntries;					// allocation is for this many entries
	int64_t m_AllocEntriesMem;					// allocation memory size
	tsPBASfxEntry *m_pEntries;					// allocated to hold entry descriptors

	int64_t m_AllocSfxIdxMem;					// allocation memory size
	void *m_pSfxIdx;							// allocated to hold sorted suffix offsets into m_pSfxSeqs[]

	int m_MaxQSortThreads;						// max number of threads to use when sorting
	CMTqsort m_MTqsort;							// multithreaded qsort

	tsPBASfxEntry* MapHit2Entry(int64_t Ofs);	// Maps the hit offset in m_pSfxSeqs[] to the relevant sequence entry

	bool		// returns true if the probe sequence exactly matches the sequence in claimed entry starting at claimed loci, otherwise false
		IsValidatedClaimedExact(uint8_t* pProbeSeq,
			int32_t ProbeLen,		// probe length
			tsPBASfxEntry* pEntry,	// claimed match was on this entry sequence 
			uint32_t HitLoci);		// starting at this loci in entry sequence

	int64_t										// index+1 in pSfxArray of first exactly matching probe or 0 if no match
		LocateFirstExact(uint8_t* pProbe,  // pts to probe sequence
			int ProbeLen,					// probe length to exactly match over
			int SfxElSize,					// size in bytes of suffix element - expected to be either 4 or 5
			int64_t TargStart,				// position in pTarg (0..n) corresponding to start of suffix array
			int64_t SfxLo,					// low index in pSfxArray
			int64_t SfxHi);					// high index in pSfxArray

	int64_t					// returned index (1..n) into m_pSfxIdx[] of first suffix sequence exactly matching pProbe, 0 if none matching
		LocateExact(uint8_t* pProbe,		// pts to probe sequence
			int ProbeLen);					// probe length to exactly match over

	int CmpProbeTarg(uint8_t* pEl1, uint8_t* pEl2, int Len);

public:
	CPBASfxArray();								// constructor
	~CPBASfxArray();							// destructor

	void Reset(void);							// reset instance state to that immediately following construction

	void SetMaxQSortThreads(int MaxThreads);			// sets maximum number of threads to use in multithreaded qsorts
	
	int											// returns the previously utilised MaxBaseCmpLen
		SetMaxBaseCmpLen(int MaxBaseCmpLen);	// sets maximum number of bases which need to be compared for equality in multithreaded qsorts, will be clamped to be in range 10..(5*cMaxReadLen)

	int	SuffixSort(void);						// construct suffix sorted index over m_pSfxSeqs[]
	int	QSortSeq(int64_t SeqLen,				// total concatenated sequence length
		uint8_t *pSeq,							// pts to start of concatenated sequences
		int SfxElSize,							// suffix element size (will be either 4 or 5)
		void* pArray);							// allocated to hold suffix elements

	teBSFrsltCodes
		AddEntry(uint32_t UserID,					// user supplied identifier
			uint32_t SeqID,							// sequence id, typically related to the textual pszSeqIdent 
			char* pszSeqIdent,						// sequence identifier, typically chromosome name
			uint8_t* pSeq,							// PBA sequence to add to suffix array
			int32_t SeqLen);						// PBA sequence length


	int64_t						// returned index (1..n) into m_pSfxIdx[] of first suffix sequence exactly matching pProbe, 0 if none matching
		IterateExacts(uint8_t* pProbeSeq,// probe
			int32_t ProbeLen,		// probe length
			int64_t PrevHitIdx,		// 0 if starting new sequence, otherwise set to return value of previous successful iteration return
			uint32_t* pUserID=nullptr,		// if match then where to return UserID in matched entry
			uint32_t* pSeqID = nullptr,		// if match then where to return SeqID in matched entry
			uint32_t* pHitLoci = nullptr);		// if match then where to return loci
		

};

