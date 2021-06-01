#pragma once
/*
 * CSIRO Open Source Software License Agreement (GPLv3)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENSE for the complete license information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */

#include "../libkit4b/commdefs.h"
#include "SSW.h"

const int cAllocTargetSeqSize = 100000;				// allocation buffer size to hold target sequences of at least this length
const int cAllocTProbetSeqSize = 100000;			// allocation buffer size to hold probe sequences of at least this length
const int cDfltMinTargetSeqLen    = 1000;			// minimum target sequence length accepted  
const int cDfltMinTargetOvlpLen   = 500;			// requiring target sequences to end overlap by at least this many bp or else be fully contained or containing  
const int cDfltMaxArtefactDev = 75;					// default percentage deviation from the mean over a 1Kbp alignment window allowed when classing overlaps as being artefactual
const int cDfltSSWDInitiatePathOfs = 250;			// default is to require SW paths to have started within this many bp on either the probe or target - effectively anchoring the SW
const int cDfltSSWDOvlpFloat = 100;					// allowed float in bp on SW overlaps

const int cMaxProbePBSSWs = 100;						// explore with SW at most this many probe alignments against target sequences
const int cPBSSWSummaryTargCoreHitCnts = 100;		// summary core hit counts on at most this many targets

#pragma pack(1)

typedef enum TAG_ePBSSWOverlapClass {
	ePBSSWOLCOverlapping = 0,	// probe classified as overlapping target, either 5' or 3'
	ePBSSWOLCcontains,			// probe completely contains the target
	ePBSSWOLCcontained,			// probe is completely contained within the target
	ePBSSWOLCartefact			// probe contains a subsequence of target, classifying as an artefact overlap and not further processed
} ePBSSWOverlapClass;

// seed core hits 
typedef struct TAG_sPBSSWCoreHit {
	uint32_t TargSeqID;				// hit was onto this target  node
	uint32_t ProbeOfs;                // hit was from this probe offset
	uint32_t TargOfs;					// onto this target offset
	uint32_t HitLen;					// hit was of this length
	uint32_t WinHits;					// number of core hits onto target relative to this core which are within a window of probelen
	uint8_t flgRevCpl:1;				// 1 if core sequence was revcpl'd before matching
	uint8_t flgMulti:1;				// 1 if core sequence was target multiloci and this instance to not be further processed
	uint8_t flgClustered:1;			// 1 if core hit identified as part of a cluster of hits
	} tsPBSSWACoreHit;

typedef struct TAG_sPBSSWCoreHitCnts {
	uint32_t TargSeqID;				// node identifier for hit target sequence	
	uint32_t	STargStartOfs;			// lowest target offset for any sense hit from probe
	uint32_t	STargEndOfs;			// highest target offset for any sense hit from probe
	uint32_t	ATargStartOfs;			// lowest target offset for any antisense hit from probe
	uint32_t	ATargEndOfs;			// highest target offset for any antisense hit from probe
	uint32_t	SProbeStartOfs;			// lowest probe offset for any sense hit onto target
	uint32_t	SProbeEndOfs;			// highest probe offset for any sense hit onto target
	uint32_t	AProbeStartOfs;			// lowest probe offset for any antisense hit onto target
	uint32_t	AProbeEndOfs;			// highest probe offset for any antisense hit onto target
	uint32_t NumSHits;				// number of hits onto target sequence from sense probe
	uint32_t NumAHits;				// number of hits onto target sequence from antisense probe
} sPBSSWCoreHitCnts;

typedef struct TAG_sPBSSWInstance {
	uint32_t InstanceID;					// each instance is uniquely identified
	uint8_t  FlgActive:1;                 // 1 if this instance is currently actively performing an alignment, new alignments can only be initiated if FlgActive == 0
	uint32_t NumCoreHits;					// currently this many core hits in m_pCoreHits
	uint32_t AllocdCoreHits;				// m_pCoreHits currently allocated to hold at most this many core hits
	size_t AllocdCoreHitsSize;			// m_pCoreHits current allocation size
	tsPBSSWACoreHit *pCoreHits;			// allocated to hold all core hits	
	uint32_t NumTargCoreHitCnts;			// current number of summary target core hit counts in TargCoreHitCnts
	sPBSSWCoreHitCnts TargCoreHitCnts[cPBSSWSummaryTargCoreHitCnts]; // top targets by core hit counts

	uint32_t ProbeSeqLen;					// current probe sequence length
	uint32_t AllocdProbeSeqSize;			// current allocation size for buffered probe sequence in pProbeSeq 	
	etSeqBase *pProbeSeq;				// allocated to hold the current probe sequence

	uint32_t TargSeqLen;					// current target sequence length
	uint32_t AllocdTargSeqSize;			// current allocation size for buffered target sequence in pTargSeq 	
	etSeqBase *pTargSeq;				// allocated to hold the current target sequence

	CMTqsort *pmtqsort;					// muti-threaded qsort
	CSSW *pSW;							// Smith-waterman class instance used by this alignment instance
} tsPBSSWInstance;

#pragma pack()

class CSWAlign
{

	int m_SWMatchScore;						// SW score for matching bases (0..100)
	int m_SWMismatchPenalty;				// SW mismatch penalty (-100..0)
	int m_SWGapOpenPenalty;					// SW gap opening penalty (-100..0)
	int m_SWGapExtnPenalty;					// SW gap extension penalty (-100..0)
	int m_SWProgExtnPenaltyLen;				// only apply gap extension penalty if gap at least this length (1..63) - use if aligning PacBio
	int m_SWDlyGapExtn;						// delayed gap penalties, only apply gap extension penalty if gap at least this length
	int m_SWProgPenaliseGapExtn;			// if non-zero then progressively increment gap extension penalty for gaps of at least this length, 0 to disable, used for PacBio style error profiles
	int m_SWAnchorLen;						// identified first and last anchors in alignment to be of at least this length

	int m_CPMatchScore;						// path class, ClassifyPath(), score for matching bases (0..100)
	int m_CPMismatchPenalty;				// path class mismatch penalty (-100..0)
	int m_CPGapOpenPenalty;					// path class gap opening penalty (-100..0)
	int m_CPGapExtnPenalty;					// path class gap extension penalty (-100..0)
	int m_MaxInitiatePathOfs;				// if non-zero then only allow new paths to start if within that offset (0 to disable) on either probe or target - effectively an anchored SW

	uint32_t m_MinOverlapLen;					// the putative overlap would be of at least this length
	uint32_t m_MaxSeedCoreDepth;				// only further extend a seed core if there are no more than this number of matching cores in all targeted sequences
	uint32_t m_DeltaCoreOfs;					// offset core windows of coreSeqLen along the probe sequence when checking for overlaps 
	uint32_t m_SeedCoreLen;					// putative overlaps are explored if there are cores of at least this length in any putative overlap
	uint32_t m_MinNumCores;					// and if the putative overlap contains at least this many cores
	uint32_t m_MaxAcceptHitsPerSeedCore;		// limit accepted hits per seed core to no more this many
	uint32_t m_MinNumSeedCores;				// require at least this many seed cores between overlapping scaffold sequences
	uint32_t m_OverlapFloat;					// allow up to this much float on overlaps to account for the PacBio error profile
	uint32_t m_MinPBSeqLen;					// individual target PacBio sequences must be of at least this length
	uint32_t m_MinPBSeqOverlap;				// any overlap of a PacBio onto a target PacBio must be of at least this many bp to be considered for contributing towards error correction (defaults to 5Kbp) 
	uint32_t m_MaxArtefactDev;				// classify overlaps as artefactual if sliding window of 1Kbp over any overlap deviates by more than this percentage from the overlap mean



	uint32_t m_DfltMaxProbeSeqLen;	// initially allocate for this length probe sequence to be aligned, will be realloc'd as may be required

	uint32_t m_NumSWInstances;			// number of SW instances currently used
	uint32_t m_AllocdSWInstances;			// allocated to hold at most this many instances
	tsPBSSWInstance *m_pSWInstances;	// allocated to hold SW instances

	uint32_t m_MaxTargSeqLen;				// longest target sequence length loaded
	CSfxArray *m_pSfxArray;			// targeted sequences are loaded/indexed into this suffix array

	int					// returns 0 if core overlapped (uses a non-exhaustive search) a previously added core, index 1..N of just added core hit or -1 if errors
	AddCoreHit(tsPBSSWInstance *pInstance,		// using this instance
				bool bRevCpl,					// true if core sequence was revcpl'd before matching
			   uint32_t ProbeOfs,                 // hit started at this probe offset
			   uint32_t TargSeqID,               // probe core matched onto this target sequence
			   uint32_t TargOfs,                  // probe core matched starting at this target loci
			   uint32_t HitLen);					// hit was of this length

	int									// returns index 1..N of core hits remaining or -1 if errors
	RemoveAddedCoreHits(tsPBSSWInstance *pInstance,		// using this instance
				int NumToRemove);                   // removing the last NumToRemove AddCoreHit() added

	int
	IdentifyCoreHits(tsPBSSWInstance *pInstance,// using this instance
					 bool bRevCpl,			// true if probe sequence to be reverse complemented
					 uint32_t MinTargLen = 1,		// minimum accepted target length
					 uint32_t MaxTargLen = 0);		// maximum accepted target length

	int LoadTargFastaFile(int MinSeqLen,		// only accept for indexing sequences of at least this length
				char *pszFile,					// file containing sequences
				int Flags = 0);					// flags are user defined

	bool m_bLocksCreated;					// set true when locks created
	int CreateLocks(void);
	void AcquireLock(bool bExclusive);
	void ReleaseLock(bool bExclusive);
	void DeleteLocks(void);

#ifdef _WIN32
	SRWLOCK m_hRwLock;
#else
	pthread_rwlock_t m_hRwLock;
#endif

static int SortCoreHitsByTargProbeOfs(const void *arg1, const void *arg2);
static int SortCoreHitsByProbeTargOfs(const void *arg1, const void *arg2);
static int SortCoreHitsDescending(const void *arg1, const void *arg2);

public:
	CSWAlign();
	~CSWAlign();

	void Reset(void);
	int Initialise(uint32_t MaxSWAInstances,	// initialise for this many instances
			uint32_t MinOverlapLen,			// the putative overlap would be of at least this length
			uint32_t MaxSeedCoreDepth,		// only further extend a seed core if there are no more than this number of matching cores in all targeted sequences
			uint32_t DeltaCoreOfs,			// offset core windows of coreSeqLen along the probe sequence when checking for overlaps 
			uint32_t CoreSeqLen,				// putative overlaps are explored if there are cores of at least this length in any putative overlap
			uint32_t MinNumCores,				// and if the putative overlap contains at least this many cores
			uint32_t MaxAcceptHitsPerSeedCore, // limit accepted hits per seed core to no more this many
			uint32_t DfltMaxProbeSeqLen);		// initially allocate for this length probe sequence to be aligned, will be realloc'd as may be required


	int									// returns number of target sequences loaded and indexed, or error result if < 0
		LoadTargetSeqs(int MinSeqLen,	// only accepting target sequences of at least this length
			char *pszTargSeqsFile,		// load target sequences from this file
			int NumThreads = 4);		// use at most this number of threads when indexing target sequences

	uint32_t			// returned alignment instance identifier or 0 if errors
			InitInstance(void);        // initialise instance - uses scores etc., as previously set with SetScores() and/or SetCPScores() and/or SetMaxInitiatePathOfs()

	bool SetScores(int MatchScore = cSSWDfltMatchScore,			// score for match
				   int MismatchPenalty = cSSWDfltMismatchPenalty,	// penalty for mismatch
				   int GapOpenPenalty = cSSWDfltGapOpenPenalty,	// penalty for opening a gap
				   int GapExtnPenalty = cSSWDfltGapExtnPenalty,	// penalty if extending already opened gap
				   int DlyGapExtn = cSSWDfltDlyGapExtn,			// delayed gap penalties, only apply gap extension penalty if gap at least this length
				   int ProgPenaliseGapExtn = cSSWDfltProgPenaliseGapExtn,	// if non-zero then progressively increment gap extension penalty for gaps of at least this length, 0 to disable, used for PacBio style error profiles
				   int AnchorLen = cSSWDfltAnchorLen);				// identified first and last anchors in alignment to be of at least this length

	bool SetCPScores(int MatchScore = cSSWDfltMatchScore,		// ClassifyPath() score for match
					 int MismatchPenalty = cSSWDfltMismatchPenalty,	// ClassifyPath() penalty for mismatch
					 int GapOpenPenalty = cSSWDfltGapOpenPenalty,	// ClassifyPath() penalty for opening a gap
					 int GapExtnPenalty = cSSWDfltGapExtnPenalty);	// ClassifyPath() penalty if extending already opened gap

	bool SetMaxInitiatePathOfs(int MaxInitiatePathOfs = cDfltSSWDInitiatePathOfs,	// require SW paths to have started within this many bp (0 to disable) on either the probe or target - effectively anchoring the SW 
								int OverlapFloat = cDfltSSWDOvlpFloat);		// with this overlap float 

	int
			AlignProbeSeq(uint32_t SWAInstance,			// using this alignment instance
							uint32_t ProbeSeqLen,			// sequence to align is this length
							etSeqBase *pProbeSeq,      // probe sequence to align
							uint32_t MinTargLen = 1,          // aligned to targets must be at least this long
							uint32_t MaxTargLen = 0,			// and if > 0 then target no longer than this many bp
							bool bSenseOnly = false,   // true if to align probe sense only, false to align both sense and antisense
	    					tsSSWCell *pRetMatched = NULL);    // optional (if not NULL) returned match detail
		
};

