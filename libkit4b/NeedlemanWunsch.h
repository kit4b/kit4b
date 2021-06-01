#pragma once
#include "./commdefs.h"

const uint32_t cNWMinProbeOrTargLen = 5;			// require probe and target lengths to be at least this number of bp
const uint32_t cNWMinCells = (cNWMinProbeOrTargLen * cNWMinProbeOrTargLen);	// not worth the effort if > this number of cells! minimum query and target length is 5bp
const uint32_t cNWMaxCells = 2000000000;		// 2 billion cells is getting big, big! restriction is that must be < 31bits, probelen * targlen limit
const uint32_t cNWMaxProbeOrTargLen = (cNWMaxCells/cNWMinProbeOrTargLen);  // require either probe or target length to be no longer than this limit

// default Needleman-Wunsch scores
const int cNWDfltMatchScore = 1;			// score for matching bases
const int cNWDfltMismatchScore = -1;		// mismatch penalty
const int cNWDfltGapOpenScore = -3;			// gap opening penalty
const int cNWDfltGapExtScore = -1;			// gap extension penalty

typedef unsigned int tNWTrcBckCell;		  // traceback cells are 32 bits
const uint32_t cNWScoreMsk     =   0x001fffff; // score is clamped to +/- 2^20-1 max with sign in cNWScoreNegFlg
const uint32_t cNWScoreNegFlg  =   0x00100000; // treat score as negative if this bit is set
const uint32_t cNWSpareBits=       0x07e00000; // these 6 bits currently not used
const uint32_t cNWGapOpnFlg  =     0x08000000; // set in cell if gap has been opened
const uint32_t cNWTrcBckMatchFlg = 0x10000000; // set if probe and target base were exactly matching
const uint32_t cNWTrcBckMsk  =     0xE0000000; // 3 bits to hold traceback direction (0=none,001=diag,010=up,100=left)
const uint32_t cNWTrcBckDiagFlg =  0x20000000; // traceback diagonally, bases align, could be exact match (cSWTrcBckMatchFlg set) or mismatched bases
const uint32_t cNWTrcBckDownFlg =  0x40000000; // traceback down, bases inserted into target (or deletion from probe)
const uint32_t cNWTrcBckLeftFlg =  0x80000000; // traceback left, bases inserted into probe (or deletion from target)


#pragma pack(1)
typedef struct TAG_sNWColBand {
	uint32_t TrcBckCellPsn;				// cells in this band have been allocated from cells starting at m_pTrcBckCells[TrcBckCellPsn-1], 0 if none allocated
	uint32_t StartTargBasePsn;			// cells in this probe column start at this target base position 1..m_TargLen
	uint32_t EndTargBasePsn;			    // cells in this probe column end at this target  1..m_TargLen 
} tsNWColBand;
#pragma pack()


class CNeedlemanWunsch
{
	bool m_bAligned;					// set true following successful Needleman-Wunsch scoring alignment
	bool m_bBanded;						// set true if banded processing

	tNWTrcBckCell *m_pTrcBckCells;		// array used for holding traceback cells
	size_t m_TrcBckCellsAllocdSize;		// allocation size
	uint32_t m_TrcBckCellsAllocd;			// actual cells in alloc'd m_pTrcBckCells (in tNWTrcBckCell's)
	uint32_t m_TrcBckCellsUsed;			// current cells used

	uint32_t m_ColBandsUsed;				// number of column active bands used - should be same as m_ProbeLen
	uint32_t m_ColBandsAllocd;			// number of column active bands allocated
	tsNWColBand *m_pColBands;			// each column, or probe base, in the matrix will contain a band of active target base cells

	etSeqBase *m_pProbe;				// alloc'd probe sequence memory
	uint32_t m_ProbeAllocd;				// actual m_pProbe alloc'd size (in teSWSeqBases's)
	uint32_t m_ProbeLen;					// current probe length

	etSeqBase *m_pTarg;					// alloc'd target sequence memory
	uint32_t m_TargAllocd;				// actual m_pTarg alloc'd size (in teSWSeqBases's)
	uint32_t m_TargLen;					// current targ length

	int m_MatchScore;					// score for matching bases
	int m_MismatchScore;				// mismatch penalty
	int m_GapOpenScore;					// gap opening penalty
	int m_GapExtScore;					// gap extension penalty

	uint32_t m_NWBandInitial;				// expecting initial overlap for max scoring alignment path to have non-overlap flank of at most this many bases
	double m_NWPathLenDiff;				// restrain path length differential between query and probe for max score be at most this proportion of query

	int m_PeakScore;					// peak score in any cell
	uint32_t m_NumBasesAligned;			// this many bases (exact and subs) were aligned between probe and target
	uint32_t m_NumBasesExact;				// of the m_NumBasesAligned, this many were exact matches, remainder were subs
	uint32_t m_NumProbeInserts;			// number of inserted bases into the probe relative to the target
	uint32_t m_NumTargInserts;			// number of inserted bases into the target relative to the probe

	uint32_t m_ProbeAlignStartOfs;		// offset in the probe at which the alignment starts
	uint32_t m_TargAlignStartOfs;			// offset in the target at which the alignment starts

	uint32_t m_PeakProbeIdx;				// highest scoring cell is at this matrix column or probe ofs
	uint32_t m_PeakTargIdx;				// highest scoring cell 

	tNWTrcBckCell *						// returns ptr to newly allocated cell  or NULL if errors
			AllocBandCell(uint32_t ProbeBasePsn,			// current probe base position 1..m_ProbeLen
							  uint32_t TargBasePsn);		// current target base position 1..m_TargLen


	tNWTrcBckCell *						// returns ptr to cell 
			DerefBandCell(uint32_t ProbeBasePsn,	// current probe base position 1..m_ProbeLen
							  uint32_t TargBasePsn);		// current target base position 1..m_TargLen

	tNWTrcBckCell *					// returns ptr to cell  or NULL if errors
			DerefBandCellLeft(uint32_t ProbeBasePsn,	// current probe base position 1..m_ProbeLen
							  uint32_t TargBasePsn);		// current target base position 1..m_TargLen


	tNWTrcBckCell *					// returns ptr to cell  or NULL if errors
			DerefBandCellDiag(uint32_t ProbeBasePsn,	// current probe base position 1..m_ProbeLen
							  uint32_t TargBasePsn);		// current target base position 1..m_TargLen

	tNWTrcBckCell *									// returns ptr to cell  or NULL if errors
			DerefBandCellDown(uint32_t ProbeBasePsn,	// current probe base position 1..m_ProbeLen
							  uint32_t TargBasePsn);		// current target base position 1..m_TargLen


public:
	CNeedlemanWunsch(void);
	~CNeedlemanWunsch(void);

	void Reset(void);

    bool SetScores(int MatchScore = cNWDfltMatchScore,		// set scoring
					int MismatchScore = cNWDfltMismatchScore,
					int GapOpenScore = cNWDfltGapOpenScore,
					int GapExtScore = cNWDfltGapExtScore);

	bool SetProbe(uint32_t Len,etSeqBase *pSeq);
	bool SetTarg(uint32_t Len,etSeqBase *pSeq);

	int Align(void);			// Needleman-Wunsch style global alignment, returns highest score

	int											// returned total alignment length between probe and target including InDels
		GetAlignStats(uint32_t *pNumAlignedBases=NULL,// returned number of bases aligning between probe and target
				 uint32_t *pNumExactBases=NULL,          // of the aligning bases there were this many exact matches, remainder were substitutions
				 uint32_t *pNumProbeInsertBases=NULL,    // this many bases were inserted into the probe relative to the target
				 uint32_t *pNumTargInsertBases=NULL);	// this many bases were inserted into the target relative to the probe
	
	int GetNumAlignedBases(void);	    // get number of bases which align (exactly plus subs), excluding InDels 
	int GetProbeAlign(uint32_t Len, etSeqBase *pBuff); // get probe alignment
	int GetTargAlign(uint32_t Len, etSeqBase *pBuff);	 // get target alignment

	int DumpScores(char *pszFile,		// dump Needleman-Wunsch matrix to this csv file
					char Down = '<',	// use this char to represent cell down link representing base inserted into target relative to probe
					char Left = '^',	// use this char to represent cell left link representing base inserted into probe relative to target
				    char Diag = '\\');	// use this char to represent cell diagonal  representing matching base either exact or mismatch

    };


