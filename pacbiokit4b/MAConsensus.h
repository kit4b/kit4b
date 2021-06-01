#pragma once

// following group are limits are reasonably arbitrary - could be increased with little risk but memory requirements could be an issue
const uint32_t cMaxNumRefSeqs = 0x0ffffff;		// can process up to this many reference sequences (16M)
const uint64_t cMaxTotRefSeqLens = 0x07ffffffff;	// which total in length to no more than this many bases (28Gbp)
const uint32_t cMinRefSeqLen = 100;				// only accepting individual reference sequences of at least this length (100bp)
const uint32_t cMaxRefSeqLen = 0x00fffffff;		// only accepting individual reference sequences no longer than this length (256Mbp)

#pragma pack(1)

// multiple alignment columns
typedef struct TAG_sMAlignConCol {
	uint32_t RefID;			// column is in this reference sequence
	uint64_t ColIdx;			// index + 1 in m_pMACols[] at which this alignment column starts
    uint64_t NxtColIdx;		// index + 1 in m_pMACols[] at which next alignment column starts, 0 if this is the last alignment column in RefID reference sequence
    uint64_t PrevColIdx;		// index + 1 in m_pMACols[] at which previous alignment column starts, 0 if this is the first alignment column in RefID reference sequence
	uint32_t Ofs;				// alignment col is for this reference sequence relative offset (1..n)
	uint16_t  Extn;			// if > 0 then relative offset to RefLoci of an inserted column 
	uint32_t  Depth;			// currently with this total coverage 
	uint8_t   ConsBase;       // initialised with the reference base, later updated when all alignments known with the consensus base
	uint32_t  BaseCnts[eBaseInDel+1];	// counts for number of bases, indexed by eBaseA..eBaseInDel  	
} tsMAlignConCol;


typedef struct TAG_sMARefSeq {
	uint32_t RefID;					// uniquely identifies this reference sequence
	uint32_t SeqLen;					// reference sequence is this length
	uint64_t StartColIdx;				// sequence starts, inclusive, from this this alignment column in m_pMACols[]
	uint64_t EndColIdx;				// sequence ends, inclusive, at this this alignment column in m_pMACols[]
	} tsMARefSeq;


#pragma pack()

class CMAConsensus
{

	bool m_bStartedMultiAlignments;     // set true after 1st call to StartMultiAlignments(), reset after call to GenMultialignConcensus() 

	uint32_t m_NumRefSeqs;			// consensus is over this many reference sequences
	uint32_t m_AllocdRefSeqs;			// m_pRefSeqs currently allocated to hold this many tsRefSeqs
	tsMARefSeq *m_pRefSeqs;			// allocated to hold reference sequence detail

	uint64_t m_MATotRefSeqLen;		// can accept up to total reference sequence length 
	uint64_t m_MACurTotRefSeqLen;		// actual current total reference sequence length
	uint64_t m_MACurCols;				// actual current number of multialignment columns used
	uint64_t m_AllocMACols;			// this number of columns have been allocated
	size_t m_AllocMAColsSize;		// allocated memory size for multialignment cols
	tsMAlignConCol *m_pMACols;				// pts to allocated memory for multialignment cols 

	tsMAlignConCol *						// inserted column or NULL if errors inserting
		InsertCol(uint64_t PrevCol);	// allocate and insert new column after PrevCol

public:
	CMAConsensus();
	~CMAConsensus();

	void Reset(void);							// reset state back to that immediately following instantiation

	int Init(uint32_t NumRefSeqs,					// max number of reference sequences which will be added
			 uint64_t TotRefSeqLen);				// max total bases of all reference sequences which will be added

	uint32_t								// returned reference sequence indentifier to be used when adding alignments to this reference sequence
		AddRefSeq(uint32_t SeqLen,		// reference sequence to which other sequences are aligned is this length
					etSeqBase *pRefSeq); // reference sequence 

	int												// eBSFSuccess or error
		AddMultiAlignment(uint32_t RefSeqID,			// alignment is against this sequence
					  uint32_t RefStartOfs,			// alignment starts at this reference sequence offset (1..n)
					  uint32_t RefEndOfs,				// alignment ends at this reference sequence offset inclusive
					  uint32_t ProbeStartOfs,			// alignment starts at this probe sequence offset (1..n)
					  uint32_t ProbeEndOfs,			// alignment ends at this probe sequence offset inclusive
					  etSeqBase *pProbeSeq,			// alignment probe sequence
					  uint32_t NumMAAlignOps,			// number of alignment operators
					   tMAOp *pMAAlignOps);			// alignment operators

	int	GenMultialignConcensus(void);          // generate multiple alignment consensus over all multialignment columns at m_pMACols

	int											// number of bases in consensus 
		GetConsensus(uint32_t RefSeqID);			// for this reference sequence identifier

	int											// number of bases returned in pRetBases 
		GetConsensus(uint32_t RefSeqID,			// consensus sequence for this reference identifier
					    uint32_t RefStartOfs,		// starting offset 1..N
						uint32_t RefLen,			// return at most this many bases
						etSeqBase *pRetBases);  // where to return bases

};

