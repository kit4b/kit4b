#pragma once
/*
 * CSIRO Open Source Software License Agreement (GPLv3)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENSE for the complete license information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */

#include "./pacbiocommon.h"

const int cMaxInFileSpecs = 100;			// user can specify upto this many input files
const int cMinSeqLen = 500;					// minimum sequence length accepted for processing
const int cMaxDupEntries = 10;				// report 1st 10 duplicate entry names
const int cMaxAllocBuffChunk = 0x00ffffff;	// buffer for fasta sequences is realloc'd in this sized chunks

const int cFlgLCSeq = 0x01;					// sequence is a low confidence read (PacBio) containing many sequence errors
const int cFlgHCSeq = 0x02;					// sequence a high confidence sequence with few sequence errors
const int cflgAlgnd = 0x04;					// sequence is either being aligned or aligned

const int cDfltMinHCSeqLen = 1000;			// default min high confidence sequence length
const int cMinHCSeqLen = 100;				// min high confidence sequence length
const int cDfltMinHCSeqOverlap = 1000;		// default min high confidence sequence overlap length onto PacBio sequence

const int cMinMaxArtefactDev = 1;			// user can specify down to this minimum or 0 to disable
const int cDfltMaxArtefactDev = 20;			// default percentage deviation from the mean allowed when classing overlaps as being artefactual
const int cDfltScaffMaxArtefactDev = 10;		// but when scaffolding with error corrected reads there should be less percentage deviation from the mean
const int cMaxMaxArtefactDev = 25;			// user can specify up to this maximum 


const int cMaxPacBioErrCorLen = 250000;		// allowing for error corrected read sequences of up to this length
const int cMaxPacBioMAFLen = (cMaxPacBioErrCorLen * 100);	// allowing for multialignment format buffering of up to this length

typedef enum TAG_ePBPMode {								// processing mode
	ePBPMErrCorrect,									// error correct
	ePBPMConsensus,										// generate consensus from previously generated multiple alignments
	ePBPMOverlapDetail									// generate overlap detail from previously generated consensus sequences
	} etPBPMode;


#pragma pack(1)

typedef struct TAG_sPBECCKMerOfs {					// no KMer in antisense strand if both MinOfs and MaxOfs are 0
	uint32_t MinOfs;										// at least one antisense KMer is located between this minimum and
	uint32_t MaxOfs;                                      // and this maximum offset (sub 1) in the antisense sequence
} tsPBECCKMerOfs;

// seed core hits 
typedef struct TAG_sPBECCCoreHit {
	uint32_t ProbeNodeID;				// core hit was from this probe  node 
	uint32_t TargNodeID;				// hit was onto this target  node
	uint32_t ProbeOfs;                // hit was from this probe offset
	uint32_t TargOfs;					// onto this target offset
	uint32_t HitLen;					// hit was of this length
	uint32_t WinHits;					// number of core hits onto target relative to this core which are within a window of probelen
	uint8_t flgRevCpl:1;				// 1 if core sequence was revcpl'd before matching
	uint8_t flgMulti:1;				// 1 if core sequence was target multiloci and this instance to not be further processed
	uint8_t flgClustered:1;			// 1 if core hit identified as part of a cluster of hits
	} tsPBECCCoreHit;

typedef struct TAG_sPBECCScaffNode {
	uint32_t NodeID;					// uniquely identifies this node
	uint32_t VertexID;				// assembly graph vertex identifier
	uint32_t EntryID;					// suffix array entry identifier for indexed sequence
	uint32_t RefSeqID;				// consensus class identifier for this sequence
	uint32_t SeqLen;					// length in bp of this scaffolding node sequence
	uint8_t flgCurProc:1;			// sequence is currently being processed
	uint8_t flgContained:1;			// sequence is fully contained within at least one other sequence
	uint8_t flgContains:1;			// sequence fully contains at least one other sequence
	uint8_t flgUnderlength:1;        // sequence is under length
	uint8_t flgHCseq:1;				// loaded as a high confidence (non-PacBio) sequence
} tsPBECCScaffNode;

typedef struct TAG_sPBECCCoreHitCnts {
	uint32_t TargNodeID;				// node identifier for hit target sequence	
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
} sPBECCCoreHitCnts;

typedef struct TAG_sThreadPBECContigs {
	int ThreadIdx;					// uniquely identifies this thread
	void *pThis;					// will be initialised to pt to class instance
#ifdef _WIN32
	HANDLE threadHandle;			// handle as returned by _beginthreadex()
	uint32_t threadID;				// identifier as set by _beginthreadex()
#else
	int threadRslt;					// result as returned by pthread_create ()
	pthread_t threadID;				// identifier as set by pthread_create ()
#endif
	int Rslt;						// processing result

	CMTqsort *pmtqsort;				// muti-threaded qsort

	CSSW *pSW;						// Smith-Waterman class

	bool bSelfHits;					// if true then processing for self hits with objective of detecting retained hairpins

	bool bRevCpl;					// true if probe sequence is to be revcpl when looking for overlaps
	uint32_t MaxSeedCoreDepth;		// only further extend a seed core if there are no more than this number of matching cores in all targeted sequences
	uint32_t DeltaCoreOfs;			// offset core windows of coreSeqLen along the probe sequence when checking for overlaps 
	uint32_t CoreSeqLen;				// putative overlaps are explored if there are cores of at least this length in any putative overlap
	uint32_t MinNumCores;				// and if the putative overlap contains at least this many cores
	uint32_t MaxAcceptHitsPerSeedCore; // limit accepted hits per seed core to no more this many
	uint32_t MinPBSeqLen;				// only process PacBio sequences which are at least this length

	uint32_t NumTargCoreHitCnts;		// current number of summary target core hit counts in TargCoreHitCnts
	sPBECCCoreHitCnts TargCoreHitCnts[cSummaryTargCoreHitCnts+1]; // top targets by core hit counts
	uint32_t NumCoreHits;				// currently this many core hits in m_pCoreHits
	uint32_t AllocdCoreHits;				// m_pCoreHits currently allocated to hold at most this many core hits
	size_t AllocdCoreHitsSize;		// m_pCoreHits current allocation size
	tsPBECCCoreHit *pCoreHits;			// allocated to hold all core hits	

	uint32_t AllocdProbeSeqSize;		// current allocation size for buffered probe sequence in pProbeSeq 	
	etSeqBase *pProbeSeq;			// allocated to hold the current probe sequence

	uint32_t AllocdTargSeqSize;		// current allocation size for buffered target sequence in pTargSeq 	
	etSeqBase *pTargSeq;			// allocated to hold the current target sequence

	uint32_t AlignErrMem;				// number of times alignments failed because of memory allocation errors
	uint32_t AlignExcessLen;			// number of times alignments failed because length of probe * target was excessive
} tsThreadPBECContigs;

#ifdef USETHISCODE
	uint32_t ErrCorBuffIdx;			// index into m_szErrCorLineBuff at which to next copy a corrected sequence
	uint32_t AllocdErrCorLineBuff;	// allocation size for m_pszErrCorLineBuff
	char *pszErrCorLineBuff;		// allocated buffering for error corrected and scored sequences
	uint32_t MultiAlignBuffIdx;		// index into m_pszMultiAlignLineBuff at which to write next multialignment
	uint32_t AllocdMultiAlignLineBuff;// allocation size for m_pszMultiAlignLineBuff
	char *pszMultiAlignLineBuff;	// allocated buffering for multialignments
#endif

#pragma pack()


class CPBECContigs
{
	etPBPMode m_PMode;						// processing mode

	uint32_t m_DeltaCoreOfs;					// offset by this many bp the core windows of coreSeqLen along the probe sequence when checking for overlaps
	uint32_t m_MaxSeedCoreDepth;				// only further process a seed core if there are no more than this number of matching cores in all targeted sequences

	uint32_t m_MinSeedCoreLen;				// use seed cores of this length when identifying putative overlapping scaffold sequences
	uint32_t m_MinNumSeedCores;				// require at least this many seed cores between overlapping scaffold sequences

	int m_SWMatchScore;						// SW score for matching bases (0..100)
	int m_SWMismatchPenalty;				// SW mismatch penalty (-100..0)
	int m_SWGapOpenPenalty;					// SW gap opening penalty (-100..0)
	int m_SWGapExtnPenalty;					// SW gap extension penalty (-100..0)
	int m_SWProgExtnPenaltyLen;				// only apply gap extension penalty if gap at least this length (1..63) - use if aligning PacBio
	uint32_t m_MaxArtefactDev;				// classify overlaps as artefactual if sliding window of 1Kbp over any overlap deviates by more than this percentage from the overlap mean
	
	uint32_t m_NumOverlapProcessed;			// number of PacBio reads processed for overlapping other PacBio reads
	uint32_t m_ProvOverlapping;               // number of high confidence sequences overlapping onto at least one contig to be error corrected
	uint32_t m_ProvContained;					// number of high confidence sequences provisionally totally contained within at least one contig to be error corrected
	uint32_t m_ProvArtefact;					// number of high confidence sequences provisionally only partially, likely an alignment artefact, contained within at least one contig to be error corrected
	uint32_t m_ProvSWchecked;					// number of times SW used to identify overlaps

	uint32_t m_OverlapFloat;					// allow up to this much float on overlaps to account for the PacBio error profile
	uint32_t m_MinContigLen;					// individual contigs must be of at least this length
	uint32_t m_MinHCSeqLen;					// only accepting hiconfidence sequences of at least this length (defaults to 1Kbp)

	char m_szContigFile[_MAX_PATH];			// input contigs file
	char m_szHiConfFile[_MAX_PATH];			// input hiconfidence sequences file		

	char m_szErrCorFile[_MAX_PATH];			// name of file into which write error corrected and scored sequences
	int m_hErrCorFile;						// file handle for writing error corrected and scored sequences

	int ErrCorBuffIdx;						// index into m_szErrCorLineBuff at which to next copy a corrected sequence
	int AllocdErrCorLineBuff;				// allocation size for m_pszErrCorLineBuff
	char *pszErrCorLineBuff;				// allocated buffering for error corrected and scored sequences
	int MultiAlignBuffIdx;					// index into m_pszMultiAlignLineBuff at which to write next multialignment
	int AllocdMultiAlignLineBuff;			// allocation size for m_pszMultiAlignLineBuff
	char *pszMultiAlignLineBuff;			// allocated buffering for multialignments


	int m_NumThreads;							// maximum number of worker threads to use

	uint32_t m_MaxTargSeqLen;						// max length of any targeted contig sequence to be error corrected
	uint32_t m_NumPBScaffNodes;					// m_pPBScaffNodes currently holds many contig nodes to be error corrected
	uint32_t m_AllocdPBScaffNodes;				// m_pPBScaffNodes allocated to hold this many contig nodes
	tsPBECCScaffNode *m_pPBScaffNodes;				// allocated to hold contig nodes
	uint32_t *m_pMapEntryID2NodeIDs;				// used to map from suffix array entry identifiers to the corresponding contig node identifier
	CSfxArray *m_pSfxArray;					// contigs to error correct are indexed in this suffix array

	uint32_t m_NumHiConfSeqs;						// m_pSeqStore currently holds many hi confidence sequences to be used for error correcting
    uint32_t m_MaxHiConfSeqLen;					// max length of any high confidence sequence
	CSeqStore  *m_pSeqStore;					// high confidence sequences used for errror correcting are loaded into this sequence store

	CMAConsensus *m_pMAConsensus;               // class to call consensus bases in contigs

	void Init(void);							// initialise state to that immediately following construction
	void Reset(bool bSync);						// reset state, if bSync true then fsync before closing output file handles

	int64_t EstSumSeqLens(int NumTargFiles,char **pszTargFiles);		// guestimate and return total sequence length by simply summing the lengths of each file - likely to grossly over estimate

	int LoadSeqs(int MinSeqLen,int NumTargFiles,
				char **pszTargFiles,			// parse, and index sequences in this file into memory resident suffix array; file expected to contain either fasta or fastq sequences
				bool bSeqStore = false,			// false is to load into suffix array, true to load into sequence store
				int Flags = cFlgLCSeq);			// which by default are low confidence PacBio read sequences

	int ProcessFastaFile(int MinSeqLen,			// only accept for indexing sequences of at least this length
				char *pszFile,					// file containing sequences
				bool bSeqStore = false,			// default is to load into suffix array, set true to load into sequence store
				int Flags = cFlgLCSeq);			// which by default are low confidence PacBio read sequences

	int InitiateECContigs(int NumECThreads);	// initiate contig error correction using this many threads

	int IdentifyCoreHits(uint32_t HiConfSeqID,	// identify all overlaps of this probe sequence HiConfSeqID onto target sequences
				uint32_t MinTargLen,				// accepted target hit sequences must be at least this length
				uint32_t MaxTargLen,				// and if > 0 then accepted targets no longer than this length
				tsThreadPBECContigs *pPars);		// thread specific

	int					// returns 0 if core overlapped (uses a non-exhaustive search) a previously added core, index 1..N of just added core hit or -1 if errors
		AddCoreHit(uint32_t HiConfSeqID,			// core hit was from this probe  
			   bool bRevCpl,					// true if core sequence was revcpl'd before matching
			   uint32_t ProbeOfs,                 // hit started at this probe offset
			   uint32_t TargNodeID,               // probe core matched onto this target scaffold node
			   uint32_t TargOfs,                  // probe core matched starting at this target loci
			   uint32_t HitLen,					// hit was of this length
               tsThreadPBECContigs *pPars);		// thread specific

	uint32_t										// returned tsPBScaffNode node identifier
		MapEntryID2NodeID(uint32_t EntryID);		// suffix array entry identifier

	CMTqsort m_mtqsort;				// muti-threaded qsort

static int SortLenDescending(const void *arg1, const void *arg2); // Sort scaffolding nodes by length descending
static int SortCoreHitsByTargProbeOfs(const void *arg1, const void *arg2);
static int SortCoreHitsByProbeTargOfs(const void *arg1, const void *arg2);
static int SortCoreHitsDescending(const void *arg1, const void *arg2);

	bool m_bMutexesCreated;			// will be set true if synchronisation mutexes have been created
	int CreateMutexes(void);
	void DeleteMutexes(void);

#ifdef _WIN32
	alignas(4)	volatile unsigned int m_CASSerialise; // used with synchronous compare and swap (CAS) for serialising access - replaces AcquireSerialise() as much more efficient
	alignas(4)	volatile unsigned int m_CASLock; // used with synchronous compare and swap (CAS) for serialising access -  - replaces AcquireLock(True) as much more efficient
#else
	__attribute__((aligned(4))) volatile unsigned int m_CASSerialise; // used with synchronous compare and swap (CAS) for serialising access - replaces AcquireSerialise() as much more efficient
	__attribute__((aligned(4))) volatile unsigned int m_CASLock; // used with synchronous compare and swap (CAS) for serialising access -  - replaces AcquireLock(True) as much more efficient
#endif
	void AcquireCASSerialise(void);
	void ReleaseCASSerialise(void);
	void AcquireCASLock(void);
	void ReleaseCASLock(void);


public:
	CPBECContigs();
	~CPBECContigs();

	int ThreadPBECContigs(tsThreadPBECContigs *pThreadPar);

	int
	Process(etPBPMode PMode,	// processing mode
		int DeltaCoreOfs,			// offset by this many bp the core windows of coreSeqLen along the probe sequence when checking for overlaps
		int MaxSeedCoreDepth,		// only further process a seed core if there are no more than this number of matching cores in all targeted sequences
		int MinSeedCoreLen,			// use seed cores of this length when identifying putative overlapping scaffold sequences
		int MinNumSeedCores,        // require at least this many seed cores between overlapping scaffold sequences
		int SWMatchScore,			// score for matching bases (0..50)
		int SWMismatchPenalty,		// mismatch penalty (-50..0)
		int SWGapOpenPenalty,		// gap opening penalty (-50..0)
		int SWGapExtnPenalty,		// gap extension penalty (-50..0)
		int SWProgExtnPenaltyLen,	// progressive gap scoring then only apply gap extension score if gap at least this length (0..63) - use if aligning PacBio
		int MaxArtefactDev,			// classify overlaps as artefactual if sliding window of 1Kbp over any overlap deviates by more than this percentage from the overlap mean
		int MinContigLen,			// only accepting contigs of at least this length (defaults to 10Kbp)
		int MinHCSeqLen,			// only accepting hiconfidence reads of at least this length (defaults to 1Kbp)
		char *pszContigFile,			// input multifasta contig file
		char *pszHiConfFile,			// input hiconfidence file
	    char *pszErrCorFile,		// name of file into which write error corrected contigs
		int NumThreads);			// maximum number of worker threads to use

};




