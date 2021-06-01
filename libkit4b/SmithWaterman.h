#pragma once
#include "./commdefs.h"

const uint32_t cSWMinProbeOrTargLen = 5;			// require probe and target lengths to be at least this number of bp
const uint32_t cSWMinCells = (cSWMinProbeOrTargLen * cSWMinProbeOrTargLen);	// not worth the effort if < this number of cells! minimum query and target length is 5bp
const uint32_t cSWMaxProbeOrTargLen = 1000000;  // require either probe or target length to be no longer than this limit
const uint64_t cSWMaxCells = ((uint64_t)cSWMaxProbeOrTargLen * cSWMaxProbeOrTargLen/10);	// 10^10 is getting big, big! lets hope this limit is never reached

// default smith-waterman scores
const int cSWDfltMatchScore = 1;		// score for matching bases
const int cSWDfltMismatchPenalty = -1;	// mismatch penalty
const int cSWDfltGapOpenPenalty = -3;	// gap opening penalty
const int cSWDfltGapExtnPenalty = -1;	// gap extension penalty
const int cSWDfltDlyGapExtn = 2;        // apply gap extension penalty for gaps of at least 2bp (must be in range 1..63)
const int cSWDfltProgPenaliseGapExtn = 0; // if non-zero then progressively increment gap extension penalty for gaps of at least this length, 0 to disable, used for PacBio style error profiles


typedef uint32_t tSWTrcBckCell;		// traceback cells are 32 bits
const  uint32_t cSWScoreMsk   =     0x001fffff; // score is clamped to 2^21-1  max
const uint32_t  cSWInDelLenMsk   =     0x07e00000; // these 6 bits used to hold the current gap opened length (clamped to be no more than 63)
const int cSWInDelLenShf = 21;						// shift factor to move InDel length into LSBs or from LSBs into cSWInDelLenMsk
const  uint32_t cSWGapOpnFlg  =     0x08000000; // set in cell if gap has been opened
const  uint32_t cSWTrcBckMatchFlg = 0x10000000; // set if probe and target base were exactly matching
const  uint32_t cSWTrcBckMsk  =     0xE0000000; // 3 bits to hold traceback direction (0=none,001=diag,010=down,100=left)
const  uint32_t cSWTrcBckDiagFlg =  0x20000000; // traceback diagonally, bases align, could be exact match (cSWTrcBckMatchFlg set) or mismatched bases
const  uint32_t cSWTrcBckDownFlg =  0x40000000; // traceback down, base insert into target (or deletion from probe)
const  uint32_t cSWTrcBckLeftFlg =  0x80000000; // traceback left, base insert into probe (or deletion from target)

#pragma pack(1)
typedef struct TAG_sSWColBand {
	uint32_t TrcBckCellPsn;				// cells in this band have been allocated from cells starting at m_pTrcBckCells[TrcBckCellPsn-1], 0 if none allocated
	uint32_t StartTargBasePsn;			// cells in this probe column start at this target base position 1..m_TargLen
	uint32_t EndTargBasePsn;			    // cells in this probe column end at this target  1..m_TargLen 
} tsSWColBand;
#pragma pack()

class CSmithWaterman
{
	bool m_bAligned;					// set true following successful Needleman-Wunsch scoring alignment
	bool m_bBanded;						// set true if banded processing
	tSWTrcBckCell *m_pTrcBckCells;		// array used for holding traceback cells
	size_t m_TrcBckCellsAllocdSize;		// allocation size
	size_t m_TrcBckCellsAllocd;			// actual cells in alloc'd m_pTrcBckCells (in tSWTrcBckCell's)
	size_t m_TrcBckCellsUsed;			// current cells used

	uint32_t m_ColBandsUsed;				// number of column active bands used - should be same as m_ProbeLen
	uint32_t m_ColBandsAllocd;			// number of column active bands allocated
	tsSWColBand *m_pColBands;			// each column, or probe base, in the matrix will contain a band of active target base cells

	etSeqBase *m_pProbe;				// alloc'd probe sequence memory
	uint32_t m_ProbeAllocd;				// actual m_pProbe alloc'd size (in teSWSeqBases's)
	uint32_t m_ProbeLen;					// current probe length

	etSeqBase *m_pTarg;					// alloc'd target sequence memory
	uint32_t m_TargAllocd;					// actual m_pTarg alloc'd size (in teSWSeqBases's)
	uint32_t m_TargLen;						// current targ length

	int m_MatchScore;					// score for matching bases (0..100)
	int m_MismatchPenalty;				// mismatch penalty (-100..0)
	int m_GapOpenPenalty;				// gap opening penalty (-100..0)
	int m_GapExtnPenalty;				// gap extension penalty (-100..0)
	int m_DlyGapExtn;					// delayed gap penalties, only apply gap extension penalty if gap at least this length
    int m_ProgPenaliseGapExtn;			// if non-zero then progressively increment gap extension penalty for gaps of at least this length, 0 to disable, used for PacBio style error profiles

	uint32_t m_SWBandInitial;				// expecting initial overlap for max scoring alignment path to have non-overlap flank of at most this many bases
	double m_SWPathLenDiff;				// restrain path length differential between query and probe for max score be at most this proportion of query

	int m_PeakScore;					// peak score in any cell
	uint32_t m_NumBasesAligned;				// this many bases (exact and subs) were aligned between probe and target
	uint32_t m_NumBasesExact;				// of the m_NumBasesAligned, this many were exact matches, remainder were subs
	uint32_t m_NumProbeInserts;				// number of inserted bases into the probe relative to the target
	uint32_t m_NumTargInserts;				// number of inserted bases into the target relative to the probe

	uint32_t m_ProbeAlignStartOfs;			// offset in the probe at which the alignment starts
	uint32_t m_TargAlignStartOfs;			    // offset in the target at which the alignment starts

	uint32_t m_PeakProbeIdx;					// highest scoring cell is at this matrix column or probe ofs
	uint32_t m_PeakTargIdx;					// highest scoring cell is at this matrix row or target ofs 

	inline tSWTrcBckCell *						// returns ptr to newly allocated cell  or NULL if errors
			AllocBandCell(uint32_t ProbeBasePsn,			// current probe base position 1..m_ProbeLen
							  uint32_t TargBasePsn);		// current target base position 1..m_TargLen


	inline tSWTrcBckCell *						// returns ptr to cell 
			DerefBandCell(uint32_t ProbeBasePsn,	// current probe base position 1..m_ProbeLen
							  uint32_t TargBasePsn);		// current target base position 1..m_TargLen

	inline tSWTrcBckCell *					// returns ptr to cell  or NULL if errors
			DerefBandCellLeft(uint32_t ProbeBasePsn,	// current probe base position 1..m_ProbeLen
							  uint32_t TargBasePsn);		// current target base position 1..m_TargLen


	inline tSWTrcBckCell *					// returns ptr to cell  or NULL if errors
			DerefBandCellDiag(uint32_t ProbeBasePsn,	// current probe base position 1..m_ProbeLen
							  uint32_t TargBasePsn);		// current target base position 1..m_TargLen

	inline tSWTrcBckCell *									// returns ptr to cell  or NULL if errors
			DerefBandCellDown(uint32_t ProbeBasePsn,	// current probe base position 1..m_ProbeLen
							  uint32_t TargBasePsn);		// current target base position 1..m_TargLen

	


public:
	CSmithWaterman(void);
	~CSmithWaterman(void);

	void Reset(void);

	bool SetScores(int MatchScore= cSWDfltMatchScore,			// score for match
				int MismatchPenalty  = cSWDfltMismatchPenalty,	// penalty for mismatch
				int GapOpenPenalty  = cSWDfltGapOpenPenalty,	// penalty for opening a gap
				int GapExtnPenalty  = cSWDfltGapOpenPenalty,	// penalty if extending already opened gap
				int DlyGapExtn = cSWDfltDlyGapExtn,				// delayed gap penalties, only apply gap extension penalty if gap at least this length
				int ProgPenaliseGapExtn = cSWDfltProgPenaliseGapExtn);	// if non-zero then progressively increment gap extension penalty for gaps of at least this length, 0 to disable, used for PacBio style error profiles


	bool SetProbe( uint32_t Len,etSeqBase *pSeq);					// set probe sequence to use in subsequent alignments
	bool SetTarg( uint32_t Len,etSeqBase *pSeq);					// set target sequence to use in subsequent alignments

	int // smith-waterman style local alignment, returns highest score
		Align(bool bBanded = false,					// true (currently experimental) to use banded or constrained SW 
				 uint32_t MaxStartNonOverlap = 20,	// if banded then initial non-overlapping path expected to be at most this many bp, 0 for no limits
				 double MaxPathLenDiff = 0.10);		// if banded then path length differential between query and probe for max score be at most this proportion of query length

	int											// returned total alignment length between probe and target including InDels
		GetAlignStats(uint32_t *pNumAlignedBases=NULL,// returned number of bases aligning between probe and target
				 uint32_t *pNumExactBases=NULL,          // of the aligning bases there were this many exact matches, remainder were substitutions
				 uint32_t *pNumProbeInsertBases=NULL,    // this many bases were inserted into the probe relative to the target
				 uint32_t *pNumTargInsertBases=NULL,		// this many bases were inserted into the target relative to the probe
				 uint32_t *pProbeStartOfs=NULL,			// alignment starts at this probe offset (1 based)
				 uint32_t *pTargStartOfs=NULL);			// alignment starts at this target offset (1 based)
	
	int GetProbeStartOfs(void);	// get offset (1..n) in probe at which alignment starts
	int GetTargStartOfs(void);  // get offset (1..n) in target at which alignment starts
	int GetNumAlignedBases(void);	    // get number of bases which align (exactly plus subs), excluding InDels 

	int GetProbeAlign( uint32_t Len, etSeqBase *pBuff); // get probe alignment
	int GetTargAlign( uint32_t Len, etSeqBase *pBuff);	 // get target alignment

	int						// < 0 errors, 0 - no anchors, 1 - anchors returned (note that retuned anchor offsets are 1 based inclusive)
		GetAnchors(uint32_t MinAnchorLen,	// anchors must be of at least this length
							   uint32_t *pProbeOfs5,	// returned probe 5' anchor starts at this probe sequence offset 
							   uint32_t *pTargOfs5,	// returned target 5' anchor starts at this target sequence offset 
							   uint32_t *pProbeOfs3,	// returned probe 3' anchor ends at this probe sequence offset
							   uint32_t *pTargOfs3);	// returned target 3' anchor ends at this target sequence offset

	int DumpScores(char *pszFile,		// dump Smith-Waterman matrix to this csv file
					char Down = '<',	// use this char to represent cell down link representing base inserted into target relative to probe
					char Left = '^',	// use this char to represent cell left link representing base inserted into probe relative to target
				    char Diag = '\\');	// use this char to represent cell diagonal  representing matching base either exact or mismatch

	};
