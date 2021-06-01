#pragma once

const uint32_t cMaxStoredSeqs = 0x7ffffff0;			// can store up to this many sequences
const uint32_t cSeqIDMsk = 0x7fffffff;				// bit 31 reserved as a user defined flag bit so sequence identifiers will have bit31 masked off by CSeqStore functions
const uint32_t cMinSeqStoreLen = 0x01;				// any individual sequence must be >= this min length
const uint32_t cMaxSeqStoreLen = 0xfffffff0;			// any individual sequence must be <= this max length
const uint32_t cAllocNumSeqHdrs = 50000;				// allocate for this incremental number of SeqHdrs
const uint32_t cAllocDescrSeqMem = (cAllocNumSeqHdrs * 1000);	// allocate for this incremental descr + seq bytes


typedef uint32_t tSeqID;		// sequence identifiers are 32bit


#pragma pack(1)
typedef struct TAG_sSeqHdr {
	tSeqID SeqID;			// uniquely identifies sequence
	uint32_t Flags;			// any flags associated with this sequence
	uint8_t DescrLen;			// sequence descriptor is this length
	uint32_t SeqLen;			// sequence is this length
	uint64_t DescrSeqOfs;		// offset into descriptor followed by sequence	
	} tsSeqHdr;
#pragma pack()

class CSeqStore
{
	uint32_t m_NumStoredSeqs;			// this many sequences currently stored
	uint32_t m_MinSeqLen;				// minimum length of any currently stored sequence
	uint32_t m_MaxSeqLen;				// maximum length of any currently stored sequence
	size_t m_TotStoredSeqsLen;		// total length of all currently stored sequences
	size_t m_UsedSeqHdrMem;			// tsSeqHdrs occupying this much memory 
	size_t m_AllocdSeqHdrSize;		// currently this much memory allocated to hold tsSeqHdrs
	tsSeqHdr *m_pSeqHdrs;			// allocated to hold tsSeqHdrs
	size_t m_UsedDescrSeqsMem;		// concatenated descriptors plus sequences currently require this much memory 
	size_t m_AllocdDescrSeqsSize;	// currently this much memory has been allocated to hold concatenated descriptors plus sequences
	uint8_t *m_pDescrSeqs;			// allocated to hold concatenated descriptors plus sequences

	uint32_t m_NumSeqDescrIdxd;		// currently this many sequence descriptors have been indexed
	size_t m_AllocdSeqDescrIdxSize;	// currently this much memory allocated to hold m_NumSeqDescrIdxd
	uint32_t *m_pSeqDescrIdx;         // allocated to hold index used to quickly retrieve sequence identifiers for given descriptor

	CMTqsort m_MTqsort;					// multi-threaded qsort
	static int SortSeqDescr(const void *arg1, const void *arg2);

	tSeqID								// sequence identifier of sequence matching on same descriptor
		LocateSeqID(char *pszDescr);	// must match on this descriptor
	
public:
	CSeqStore();
	~CSeqStore();
	int Reset(void);

	tSeqID			// identifier by which this sequence can later be retrieved (0 if unable to add sequence)
		AddSeq(uint32_t Flags,			// any flags associated with this sequence
			   char *pszDescr,			// sequence descriptor
			   uint32_t SeqLen,			// sequence is this length
			   etSeqBase *pSeq);		// sequence to add

	int GenSeqDescrIdx(void);			// generate index over all sequence descriptors, index used when retrieving sequence identifier for given descriptor by GetSeqID()


	tSeqID								// sequence identifier of sequence matching on same descriptor
		GetSeqID(char *pszDescr);		// must match on this descriptor

	uint32_t								// returned number of currently stored sequences
			GetNumSeqs(void);			// get number of currently stored sequences

	size_t								// returned total length of all currently stored sequences
			GetTotalLen(void);			// get total length of all currently stored sequences

	uint32_t								// returned max length of any currently stored sequence
			GetMaxSeqLen(void);			// get max length of any currently stored sequence

	uint32_t								// returned min length of any currently stored sequence
			GetMinSeqLen(void);			// get min length of any currently stored sequence

	uint32_t								// returned flags
			GetFlags(tSeqID SeqID);		// get flags associated with sequence identified by SeqID

	uint32_t								// returned sequence length
			GetLen(tSeqID SeqID);		// get sequence length identified by SeqID

	uint32_t								// returned flags
			GetDescr(tSeqID SeqID,		// get descriptor associated with sequence identified by SeqID
					int MaxDescrLen,	// return at most this many descriptor chars copied into
					char *pszDescr);	// this returned sequence descriptor

	uint32_t								// previous flags
			SetFlags(tSeqID SeqID,		// set flags associated with sequence identified by SeqID
					uint32_t Flags);		// flags to be associated with this sequence

	uint32_t								// returned sequence copied into pRetSeq is this length; 0 if errors
			GetSeq(tSeqID SeqID,		// get sequence identified by SeqID
				uint32_t StartOfs,		// starting from this sequence offset
				uint32_t MaxSeqSize,		// return at most this many sequence bases  
				etSeqBase *pRetSeq);	// copy sequence into this caller allocated buffer 
};

