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

#pragma once

const int cMinBlitzCoreLen = 8;			// minimum allowed core or seed length
const int cDfltBlitzCoreLen = 20;		// default core or seed length
const int cMaxBlitzCoreLen = 100;		// max allowed core or seed length

const int cDfltBlitzMinQueryLenAlignedPct = 75;  // to be accepted a query sequence must align over at least this percentage (1..100) of it's length onto target

const int cMinBlitzPathScore = 25;		// user specified min allowed path score before that path will be reported
const int cDfltBlitzPathScore = 75;		// default minimum path score before that path will be reported
const int cMaxBlitzPathScore = 100000;		// user specified max minimum allowed path score before that path will be reported

const int cMaxBlitzOverlapFloat = 8;		// allowing for overlap float of at most this many bases, needed because the cores with extensions are independent and may be overextended

const int cDfltBlitzMaxPathsToReport = 1;	// by default report at most this many scoring paths for any query sequence


const int cDfltBlitzExactMatchScore = 1;		// default score for an exact match
const int cMinBlitzExactMatchScore = 1;			// user specified minimum score for an exact match
const int cMaxBlitzExactMatchScore = 10;		// user specified maximum score for an exact match

const int cDfltBlitzMismatchScore = 2;		// default cost for mismatches when scoring
const int cMinBlitzMismatchScore = 1;		// user specified minimum cost for mismatches when scoring
const int cMaxBlitzMismatchScore = 10;		// user specified maximum cost for mismatches when scoring

// using affine gap scoring but limiting the gap extension cost to just the first 100bp
const int cDfltBlitzGapOpenScore = 5;	// default cost for opening path gap when scoring path
const int cMinBlitzGapOpenScore = 1;		// user specified minimum cost for opening path gap when scoring path
const int cMaxBlitzGapOpenScore = 50;	// user specified maximum cost for opening path gap when scoring path

const int cGapBlitzExtendCost = 1;		// cost for extending gap per 10bp extension when scoring path
const int cGapBlitzExtendCostLimit = 10;	// clamp gap extension cost to be no more than this
const int cGapBlitzMaxLength = 500000;   // treat any gaps longer than this length as being not on same path - allows for RNA-seq with skipped exons

const int cMinBlitzCoreDelta = 1;		// minimum allowed core shift delta in bp
const int cMaxBlitzCoreDelta = 50;		// max allowed core shift delta in bp

const int cDfltBlitzSensCoreIters  = 1500;	// default sensitivity core explore depth 
const int cMoreBlitzSensCoreIters  = 2000;	// more sensitivity core explore depth
const int cUltraBlitzSensCoreIters = 3000;	// ultra sensitivity core explore depth
const int cMinBlitzSensCoreIters   = 750;	// min sensitivity core explore depth

const int cMinBlitzOccKMerDepth = 100;       // user can override the core exploration search depth from this minimum 
const int cMaxBlitzOccKMerDepth = 20000;		// up to this maximum 

const int cMaxBlitzQuerySeqIdentLen = 80;		// allow for fasta sequence identifiers of upto this length 
const int cMaxBlitzDescrLen = 128;				// allow for fasta descriptors (incl identifiers) of up to this length
const int cSAMBlitztruncSeqLen = 16000;			// SAM format is really designed for short read sequences, down stream apps may have problems handling alignments with long query sequences so slough sequences longer than this length if SAM output
const int cAllocBlitzQuerySeqLen = 0x0200000;	// initially allocate to hold a query sequence of up to this length (2Mbp), will be realloc'd if needed for longer query sequences
const int cMaxBlitzQuerySeqLen = (cAllocBlitzQuerySeqLen * 8);    // can handle query sequences of up to this maximal length (16Mbp), longer sequences will be truncated to this length and the user warned
const int cMaxBlitzReadAheadQuerySeqs = 50000;	// read ahead and enqueue up to at most this many query sequences

const int cNumBlitzAllocdAlignNodes = 200000;  // allow each query sequence to have up to this many aligned subsequences

const int cAlignBlitzRprtBufferSize = 500000; // buffer for buffering alignment results ready to write to file

#pragma pack(1)

typedef enum TAG_eBLZPMode {
	eBLZPMdefault,						// default processing mode
	eBLZPMplaceholder					// used as a placeholder and sets the range of these enumerations
}etBLZPMode;


typedef enum TAG_eBLZSensitivity {
	eBLZSdefault = 0,					// default processing sensitivity
	eBLZSMoreSens,						// more sensitive - slower
	eBLZSUltraSens,					// ultra sensitive - much slower
	eBLZSLessSens,						// less sensitive - quicker
	eBLZSplaceholder					// used as a placeholder and sets the range of these enumerations
}etBLZSensitivity;


typedef enum TAG_eBLZRsltsFomat {
	eBLZRsltsPSL = 0,	// default results format is PSL
	eBLZRsltsPSLX,		// results format is PSLX
	eBLZRsltsMAF,		// default results format is MAF
	eBLZRsltsBED,		// results as BED
	eBLZRsltsSQLite,	// results as SQLite database
	eBLZRsltsSAM,		// results in SAM format
	eBLZRsltsplaceholder   // used as a placeholder and flags the range of these enumerations
}etBLZRsltsFomat;


typedef struct TAG_sQuerySeq {
    int SeqID;						// monotonically increasing unique sequence identifier
	char szQueryIdent[cMaxBlitzQuerySeqIdentLen+1];	// fasta identifier
	int QuerySeqLen;				// query sequence length
	uint8_t *pQuerySeq;				// allocated to hold sequence 
} tsQuerySeq;

typedef struct TAG_sLoadQuerySeqsThreadPars {
	int ThreadIdx;					// uniquely identifies this thread
	void *pThis;					// will be initialised to pt to CBlitz instance

#ifdef _WIN32
	HANDLE threadHandle;			// handle as returned by _beginthreadex()
	unsigned int threadID;			// identifier as set by _beginthreadex()
#else
	int threadRslt;					// result as returned by pthread_create ()
	pthread_t threadID;				// identifier as set by pthread_create ()
#endif
	bool bIsSAMOutput;				// will be true if aligning for SAM output format with PE processing capability
	bool bIsSAMPE;					// will be true if aligning SAM PE
	int *pRslt;						// write intermediate result codes to this location
	int Rslt;						// returned result code
} tsLoadQuerySeqsThreadPars;

typedef struct TAG_sThreadQuerySeqsPars {
	int ThreadIdx;					// uniquely identifies this thread
	void *pThis;					// will be initialised to pt to CBlitz instance

#ifdef _WIN32
	HANDLE threadHandle;			// handle as returned by _beginthreadex()
	unsigned int threadID;			// identifier as set by _beginthreadex()
#else
	int threadRslt;					// result as returned by pthread_create ()
	pthread_t threadID;				// identifier as set by pthread_create ()
#endif
	bool bIsSAMOutput;				// will be true if aligning for SAM output format with PE processing capability
	bool bIsSAMPE;					// will be true if aligning SAM PE
	uint32_t NumAllocdAlignNodes;				// number of allocated alignment nodes
	tsQueryAlignNodes *pAllocdAlignNodes;	// allocated to hold aligned subsequences
	tsQueryAlignNodes **ppFirst2Rpts;		// allocated to hold ptrs to alignment nodes which are marked as being FlgFirst2tRpt
	uint32_t NumAllocdAlignNodesPE2;				// number of allocated alignment nodes
	tsQueryAlignNodes *pAllocdAlignNodesPE2;	// allocated to hold aligned PE2 subsequences
	tsQueryAlignNodes **ppFirst2RptsPE2;		// allocated to hold ptrs to PE2 alignment nodes which are marked as being FlgFirst2tRpt
	int *pRslt;						// write intermediate result codes to this location
	int Rslt;						// returned result code
} tsThreadQuerySeqsPars;

#pragma pack()


class CBlitz
{
	etBLZPMode m_ProcMode;			// processing mode
	int m_SampleNthRawRead;			// sample every Nth raw read (or read pair) for processing (1..10000)
	etBLZSensitivity m_Sensitivity;  // alignment sensitivity 0 - standard, 1 - high, 2 - very high, 3 - low sensitivity
	eALStrand m_AlignStrand;		// align on to watson, crick or both strands of target
	int m_CoreLen;					// use this core length as the exactly matching seed length to be 5' and 3' extended whilst no more than m_MaxSubRate
	int m_CoreDelta;				// offset cores by this many bp
	int m_MaxInsertLen;				// SAM output, accept observed insert sizes of at most this (default = 100000)
	int m_QueryLenAlignedPct;    	// only report alignment paths if the percentage of total aligned bases to the query sequence length is at least this percentage (1..100)
	int m_MaxIter;					// max allowed iterations (depth) per subsegmented sequence (core) when matching that subsegment
	int m_MaxOccKMerDepth;			// maximum depth to explore over-occurring core K-mers
	int m_MinPathScore;				// to be reported alignment paths on any target sequence must score at least this value
	int m_MaxPathsToReport;			// report at most this many alignment paths for any query
	int m_AlignPathID;				// alignment path identifier

	int m_MismatchScore;			// decrease score by this for each mismatch bp
	int m_ExactMatchScore;			// increase score by this for each exactly matching bp
	int m_GapOpenScore;				// decrease score by this for each gap open

	int m_MaxQuerySeqLen;			// query sequences over this length will be sloughed, required because SAM format is assuming relatively short NGS reads as query sequences

	etBLZRsltsFomat m_RsltsFormat;	// output results format
	char *m_pszInputFile;			// name of input file containing query sequences, or PE1 if PE processing when SAM output
	char *m_pszInputFilePE2;		// name of input file containing PE2 query sequences if PE processing when SAM output

	char *m_pszSfxFile;				// target as suffix array
	char *m_pszOutFile;				// where to write alignments
	uint32_t m_ReportedPaths;			// total number of aligned paths reported
	uint32_t m_QueriesPaths;			// this many query sequences had at least one reported path
	uint32_t m_NumQueriesProc;		// total number of query sequences processed

	CSQLitePSL *m_pSQLitePSL;		// used when outputting PSL rows directly into SQLite database
	int m_ExprID;					// experiment identifier allocated when initialising database with experiment details

	int m_hOutNonAlignedFilePE1;	// results output file handle for non-aligned PE1 reads
	int m_hOutNonAlignedFilePE2;	// results output file handle for non-aligned PE2 reads

	int m_hOutFile;					// results output file handle
	int m_szLineBuffIdx;			// offset into m_pszLineBuff at which to next write
	char *m_pszLineBuff;			// allocated to hold output line buffering
	CSfxArray *m_pSfxArray;			// suffix array holds genome of interest
	char m_szTargSpecies[cMaxDatasetSpeciesChrom+1]; // suffix array was generated over this targeted species

	int m_TotSeqIDs;				// total number of query sequences which have been parsed and enqueued
	int m_NumQuerySeqs;				// number of query sequences currently in m_pQuerySeqs
	int m_NxtQuerySeqIdx;			// index into m_pQuerySeqs[] at which to dequeue the next query sequence
	int m_AllocdQuerySeqs;			// number of query sequences allocated
	tsQuerySeq *m_pQuerySeqs;		// allocated to hold array of query sequences

	uint32_t *m_pKmerOccsDist;		// allocated to hold Kmer count distributions (up to cMaxOccKMerDepth counts)
	
	uint8_t m_TermBackgoundThreads; // if non-zero then all background threads are to immediately terminate processing


	int m_NumThreads;				// number of worker threads to use

	bool m_bAllQuerySeqsLoaded;			// set true when all query sequences have been parsed and loaded
	teBSFrsltCodes m_LoadQuerySeqsRslt;	// set with exit code from background query sequences load thread, read after checking if m_bAllQuerySeqsLoaded has been set
	int m_ThreadLoadQuerySeqsRslt;		// returned by query sequence loading thread
#ifdef _WIN32
	unsigned int m_ThreadLoadQuerySeqsID;	// query sequences loading thread identifier
#else
	pthread_t m_ThreadLoadQuerySeqsID;
#endif


	CMTqsort m_mtqsort;				// muti-threaded qsort

	void Init(void);			// initialise state to that immediately following construction
	void Reset(bool bSync);		// reset state, if bSync true then fsync before closing output file handles

	int InitQuerySeqThreads(int NumThreads,			// use this many threads
							int AlignNodes);			// each thread is allocatd this many subsequence alignment nodes


	int InitLoadQuerySeqs(void);		// query sequences are loaded asynchronously to the alignments

	bool m_bMutexesCreated;			// will be set true if synchronisation mutexes have been created
	int CreateMutexes(void);
	void DeleteMutexes(void);
	void AcquireSerialise(void);
    void ReleaseSerialise(void);
    void AcquireSerialiseMH(void);
	void ReleaseSerialiseMH(void);
	void AcquireLock(bool bExclusive);
	void ReleaseLock(bool bExclusive);


#ifdef _WIN32
	HANDLE m_hMtxIterReads;
	HANDLE m_hMtxMHReads;
	SRWLOCK m_hRwLock;
	HANDLE m_hThreadLoadQuerySeqs;
#else
	pthread_mutex_t m_hMtxIterReads;
	pthread_mutex_t m_hMtxMHReads;
	pthread_rwlock_t m_hRwLock;
#endif

	static int SortQueryAlignNodes(const void *arg1, const void *arg2); // Sort alignment nodes by TargSeqID.QueryID.FlgStrand.QueryStartOfs.TargStartOfs ascending
	static int SortHighScoreDescend(const void *arg1, const void *arg2); // Sort alignment nodes which are the first in path by score descending

public:
	CBlitz();
	~CBlitz();
	int
	Process(char *pszProcessName,		// name of process requesting blitz alignments
			etBLZPMode PMode,			// processing mode
			int SampleNthRawRead,		// sample every Nth raw read (or read pair) for processing (1..10000)
			char *pszExprName,				// experiment name
			char *pszExprDescr,				// experiment description
			char *pszParams,				// string containing blitz parameters
			bool KMerDist,				// true if K_mer counts distributions to be reported
			etBLZSensitivity Sensitivity,	// sensitivity 0 - standard, 1 - high, 2 - very high, 3 - low sensitivity
			eALStrand AlignStrand,			// align on to watson, crick or both strands of target
			int MismatchScore,				// decrease score by this for each mismatch bp
			int ExactMatchScore,			// increase score by this for each exactly matching bp
			int GapOpenScore,				// decrease score by this for each gap open
			int  CoreLen,					// use this core length as the exactly matching seed length to be 5' and 3' extended whilst no more than m_MaxSubRate
			int  CoreDelta,					// offset cores by this many bp
			int MaxInsertLen,				// SAM output, accept observed insert sizes of at most this (default = 50000)
			int MaxOccKMerDepth,			// maximum depth to explore over-occurring core K-mers
			int  MinPathScore,				// only report alignment paths on any target sequence if the path score is >= this minimum score
			int QueryLenAlignedPct,			// only report alignment paths if the percentage of total aligned bases to the query sequence length is at least this percentage (1..100)
			int  MaxPathsToReport,			// report at most this many alignment paths for any query
			etBLZRsltsFomat RsltsFormat,	// output results format
			char *pszInputFile,				// name of input file containing query sequences (PE1 if PE processing when output format is SAM)
			char *pszInputFilePE2,			// name of input file containing PE2 query sequences (only applies if output format is SAM)
			char *pszSfxFile,				// target as suffix array
			char *pszOutFile,				// where to write alignments
			int NumThreads);				// number of worker threads to use

	int ProcLoadQuerySeqsFile(tsLoadQuerySeqsThreadPars *pPars);
	int ProcLoadSAMQuerySeqsFile(tsLoadQuerySeqsThreadPars *pPars);
	int ProcAlignQuerySeqs(tsThreadQuerySeqsPars *pPars);
	int ProcAlignSAMQuerySeqsSE(tsThreadQuerySeqsPars *pPars);	// align as SE only reporting alignments in SAM file format
	int ProcAlignSAMQuerySeqsPE(tsThreadQuerySeqsPars *pPars);	// align as PE only reporting alignments in SAM file format
	int	ReportNonAligned(char *pszDescPE1, int LenSeqPE1, uint8_t *pSeqPE1, char *pszDescPE2 = NULL, int LenSeqPE2 = 0, uint8_t *pSeqPE2 = NULL);

	teBSFrsltCodes // When SAM aligning then will be aligning sequencing reads which are length truncated to be no longer than cSAMtruncSeqLen bases
		LoadSAMRawReads(bool bIsPairReads,	// true if paired end processing - PE1 reads in pszPE1File and PE2 reads in pszPE2File
						char *pszPE1File,	// process PE1 reads from this file
						char *pszPE2File);	// optionally process PE2 reads from this file

	int										// returned enqueued query identifier
		EnqueueQuerySeq(char *pszQueryIdent,    // query identifier
			int QuerySeqLen,				// query sequence length
			uint8_t *pQuerySeq,				// query sequence
			char *pszQueryIdentPE2 = NULL,    // PE2 query identifier parsed from fasta descriptor
			int QuerySeqLenPE2 = 0,				// PE2 query sequence length
			uint8_t *pQuerySeqPE2 = NULL);       // PE2 query sequence


	int										// number of sequences returned, 0 if none to be dequeued, < 0 if errors
		DequeueQuerySeq(int MaxLenQueryIdent,			// maximum length query identifier
			int *pSeqID,					// returned sequence identifier
			char *pszQueryIdent,			// where to return query identifier
			int *pQuerySeqLen,				// where to return query sequence length
			uint8_t **pDequeuedQuerySeq,		// where to return ptr to dequeued query sequence, caller is responsible for deleting memory allocated to hold the returned sequence (delete *pDequeuedQuerySeq)
			int *pSeqIDPE2 = NULL,			// if SAM PE: returned PE2 sequence identifier
			char *pszQueryIdentPE2 = NULL,	// if SAM PE: where to return PE2 query identifier
			int *pQuerySeqLenPE2 = NULL,	// if SAM PE: where to return PE2 query sequence length
			uint8_t **pDequeuedQuerySeqPE2 = NULL);	// if SAM PE: where to return ptr to dequeued PE2 query sequence, caller is responsible for deleting memory allocated to hold the returned sequence (delete *pDequeuedQuerySeqPE2)

	int
		BlocksAlignStats(uint32_t *pMatches,	// returned number of bases that match that aren't repeats
					uint32_t *pmisMatches,	// returned number of bases that don't match
					uint32_t *prepMatches,	// returned number of bases that match but are part of repeats
					uint32_t *pnCount,		// returned number of 'N' bases
				char  Strand,				// query sequence strand, '+' or '-')
				uint8_t *pQuerySeq,			// the query sequence
				uint32_t qSize,				// Query sequence size
				uint32_t TargSeqID,			// CSfxArray sequence identifier
				uint32_t tSize,				// Target sequence size 
				uint32_t TargPathStartOfs,	// at this starting offset
				uint32_t TargPathEndOfs,		// ending at this offset (inclusive)
				uint32_t NumPathNodes,		// number of alignment nodes in alignment path
				int SortedPathIdx,
				tsQueryAlignNodes *pAlignNodes,		// alignment nodes
				tsQueryAlignNodes **ppFirst2Rpts); // allocated to hold ptrs to alignment nodes which are marked as being FlgFirst2tRpt

	// reporting KMer counts distributions as CSV
	int	ReportKMerDist(uint32_t KMerLen,		// KMer count distributions are for this length K-mer 
			uint32_t MaxKMerCnts,				// max number of target occurrences which ranges from 0 to MaxKMerCnts inclusive
			uint32_t *pKmerOccsDist);			// array of length MaxKMerCnts + 1 which holds number of query K-mers having target occurrences

	uint32_t							// returned number of high scoring path head nodes
		IdentifyHighScorePaths(int MinPathScore,			// only report paths having at least this minimum score
			int  MaxPathsToReport,		// report at most this many alignment paths for any query
			uint32_t QueryLen,			// query length
			uint32_t NumNodes,			// number of alignment nodes
			tsQueryAlignNodes *pAlignNodes,		// alignment nodes
			tsQueryAlignNodes **ppFirst2Rpts);	// allocated to hold ptrs to alignment nodes which are marked as being FlgFirst2tRpt

	uint32_t	// returns number of nodes in characterised path
		CharacterisePath(uint32_t QueryLen,			// query length
				uint8_t* pQuerySeq,					// the query sequence as etSeqBase's
				uint32_t SortedPathIdx,				// index of paths head node
				tsQueryAlignNodes *pAlignNodes,		// alignment nodes
				tsQueryAlignNodes **ppFirst2Rpts,	// allocated to hold ptrs to alignment nodes which are marked as being FlgFirst2tRpt
				uint32_t *pQueryPathEndOfs = NULL,			// query path ends at this offset
				uint32_t *pTargPathEndOfs = NULL,			// target path ends at this offset
				uint32_t *pqNumInsert = NULL,			// qNumInsert, Number of inserts in query
				uint32_t *pqBaseInsert = NULL,		// qBaseInsert, Number of bases inserted in query
				uint32_t *ptNumInsert = NULL,			// tNumInsert, Number of inserts in target
				uint32_t *tBaseInsert = NULL);		// tBaseInsert, Number of bases inserted in target

	// Frequently there will be gaps between nodes and these gaps need to be closed through interpolation
	uint32_t			// number of nodes after consolidation
		ConsolidateNodes(bool bSense,		// if true then query sequence is aligning antisense to target
				uint32_t QueryLen,					// query length
				uint8_t* pQuerySeq,					// the query sequence as etSeqBase's
				uint32_t SortedPathIdx,				// alignment path starts from this node - 0 based 
				tsQueryAlignNodes* pAlignNodes,		// alignment nodes
				tsQueryAlignNodes** ppFirst2Rpts);	// allocated to hold ptrs to alignment nodes which are marked as being FlgFirst2tRpt


	int			// number of reported paths 
		Report( uint32_t MinPathScore,			// only report paths having at least this minimum score
				uint32_t  MaxPathsToReport,		// report at most this many alignment paths for any query
				char *pszQuerySeqIdent,		// query sequence name 
				uint32_t QueryLen,			// query length
				uint8_t *pQuerySeq,			// the query sequence
				uint32_t NumNodes,			// number of alignment nodes
				tsQueryAlignNodes *pAlignNodes, // alignment nodes
				tsQueryAlignNodes **ppFirst2Rpts);	// allocated to hold ptrs to alignment nodes which are marked as being FlgFirst2tRpt




	int	// reporting alignment as SQLite PSL format 
		ReportAsSQLitePSL(uint32_t Matches,				// Number of bases that match that aren't repeats
						uint32_t misMatches,			// Number of bases that don't match
						uint32_t repMatches,			// Number of bases that match but are part of repeats
						uint32_t nCount,				// Number of 'N' bases
						uint32_t	qNumInsert,			// Number of inserts in query
						uint32_t qBaseInsert,			// Number of bases inserted in query
						uint32_t tNumInsert,			// Number of inserts in target
						uint32_t tBaseInsert,			// Number of bases inserted in target
						char  Strand,				// query sequence strand, '+' or '-'
						char *pszQuerySeqIdent,     // this query sequence
						uint32_t qSize,				// Query sequence size
						uint32_t qStart,				// Alignment start position in query
						uint32_t qEnd,				// Alignment end position in query
						char *pszTargName,			// aligning to this target
						uint32_t tSize,				// Target sequence size 
						uint32_t TargPathStartOfs,	// at this starting offset
						uint32_t TargPathEndOfs,		// ending at this offset (inclusive)
						uint32_t NumPathNodes,		// number of alignment nodes in alignment path
						int SortedPathIdx,
						tsQueryAlignNodes *pAlignNodes,		// alignment nodes
						tsQueryAlignNodes **ppFirst2Rpts); // allocated to hold ptrs to alignment nodes which are marked as being FlgFirst2tRpt


	int		// reporting alignment as PSL format
			ReportAsPSL(uint32_t Matches,			// Number of bases that match that aren't repeats
					uint32_t misMatches,			// Number of bases that don't match
					uint32_t repMatches,			// Number of bases that match but are part of repeats
					uint32_t nCount,				// Number of 'N' bases
					uint32_t	qNumInsert,			// Number of inserts in query
					uint32_t qBaseInsert,			// Number of bases inserted in query
					uint32_t tNumInsert,			// Number of inserts in target
					uint32_t tBaseInsert,			// Number of bases inserted in target
					char  Strand,				// query sequence strand, '+' or '-'
					char *pszQuerySeqIdent,     // this query sequence
					uint32_t qSize,				// Query sequence size
					uint32_t qStart,				// Alignment start position in query
					uint32_t qEnd,				// Alignment end position in query
					char *pszTargName,			// aligning to this target
					uint32_t tSize,				// Target sequence size 
					uint32_t TargPathStartOfs,	// at this starting offset
					uint32_t TargPathEndOfs,		// ending at this offset (inclusive)
					uint32_t NumPathNodes,		// number of alignment nodes in alignment path
					int SortedPathIdx,
					tsQueryAlignNodes *pAlignNodes,		// alignment nodes
					tsQueryAlignNodes **ppFirst2Rpts);  // allocated to hold ptrs to alignment nodes which are marked as being FlgFirst2tRpt

	int		// reporting alignment as PSLX format
			ReportAsPSLX(uint32_t Matches,			// Number of bases that match that aren't repeats
					uint32_t misMatches,			// Number of bases that don't match
					uint32_t repMatches,			// Number of bases that match but are part of repeats
					uint32_t nCount,				// Number of 'N' bases
					uint32_t	qNumInsert,			// Number of inserts in query
					uint32_t qBaseInsert,			// Number of bases inserted in query
					uint32_t tNumInsert,			// Number of inserts in target
					uint32_t tBaseInsert,			// Number of bases inserted in target
					char  Strand,				// query sequence strand, '+' or '-'
					char *pszQuerySeqIdent,     // this query sequence
					uint32_t qSize,				// Query sequence size
					uint8_t *pQuerySeq,			// the query sequence
					uint32_t qStart,				// Alignment start position in query
					uint32_t qEnd,				// Alignment end position in query
					uint32_t TargSeqID,			// CSfxArray sequence identifier
					char *pszTargName,			// aligning to this target
					uint32_t tSize,				// Target sequence size 
					uint32_t TargPathStartOfs,	// at this starting offset
					uint32_t TargPathEndOfs,		// ending at this offset (inclusive)
					uint32_t NumPathNodes,		// number of alignment nodes in alignment path
					int SortedPathIdx,
					tsQueryAlignNodes *pAlignNodes,		// alignment nodes
					tsQueryAlignNodes **ppFirst2Rpts);  // allocated to hold ptrs to alignment nodes which are marked as being FlgFirst2tRpt

	int		// reporting alignment as SAM format
		ReportAsSAM(uint32_t Flags,		// use as reported SAM flags
			char  Strand,				// query sequence strand, '+' or '-'
			uint32_t AlignScore,			// alignment has this score
			char *pszQuerySeqIdent,     // this query sequence
			uint32_t qSize,				// Query sequence size
			uint8_t *pQuerySeq,			// the query sequence
			uint32_t qStart,				// Alignment start position in query
			uint32_t qEnd,				// Alignment end position in query
			char *pszTargName,			// aligning to this target
			uint32_t tSize,				// Target sequence size 
			uint32_t TargPathStartOfs,	// at this starting offset
			uint32_t TargPathEndOfs,		// ending at this offset (inclusive)
			char RNEXT,					// Reference sequence name of the primary alignment of the NEXT read in the template, '*' if unknown, '=' if PE and both ends align to same reference
			uint32_t PNEXT,				// 1-based Position of the primary alignment of the NEXT read in the template. Set as 0 when the information is unavailable
			int TLEN,					// signed template length, If all segments are mapped to the same reference, the unsigned observed template length equals the number of bases from the leftmost mapped base to the rightmost mapped base
			uint32_t NumPathNodes,		// number of alignment nodes in alignment path
			int SortedPathIdx,
			tsQueryAlignNodes *pAlignNodes,		// alignment nodes
			tsQueryAlignNodes **ppFirst2Rpts);  // allocated to hold ptrs to alignment nodes which are marked as being FlgFirst2tRpt

	int		// reporting alignment as MAF format
			ReportAsMAF(int PathScore,			// score for this path
					uint32_t Matches,				// Number of bases that match that aren't repeats
					uint32_t misMatches,			// Number of bases that don't match
					uint32_t repMatches,			// Number of bases that match but are part of repeats
					uint32_t nCount,				// Number of 'N' bases
					uint32_t	qNumInsert,			// Number of inserts in query
					uint32_t qBaseInsert,			// Number of bases inserted in query
					uint32_t tNumInsert,			// Number of inserts in target
					uint32_t tBaseInsert,			// Number of bases inserted in target
					char  Strand,				// query sequence strand, '+' or '-'
					char *pszQuerySeqIdent,     // this query sequence
					uint32_t qSize,				// Query sequence size
					uint8_t *pQuerySeq,			// the query sequence
					uint32_t qStart,				// Alignment start position in query
					uint32_t qEnd,				// Alignment end position in query
					uint32_t TargSeqID,			// CSfxArray sequence identifier
					char *pszTargName,			// aligning to this target
					uint32_t tSize,				// Target sequence size 
					uint32_t TargPathStartOfs,	// at this starting offset
					uint32_t TargPathEndOfs,		// ending at this offset (inclusive)
					uint32_t NumPathNodes,		// number of alignment nodes in alignment path
					int SortedPathIdx,
					tsQueryAlignNodes *pAlignNodes,		// alignment nodes
					tsQueryAlignNodes **ppFirst2Rpts);  // allocated to hold ptrs to alignment nodes which are marked as being FlgFirst2tRpt
	
	int			// reporting alignment as BED format
			ReportAsBED(char *pszQuerySeqIdent,     // this query sequence
					char  Strand,				// query sequence strand, '+' or '-'
					uint32_t AlignScore,			// alignment has this score
					char *pszTargName,			// aligning to this target
					uint32_t NumPathNodes,		// number of alignment nodes in alignment path
					uint32_t TargPathStartOfs,	// at this starting offset
					uint32_t TargPathEndOfs,		// ending at this offset (inclusive)
					int SortedPathIdx,
					tsQueryAlignNodes *pAlignNodes,		// alignment nodes
					tsQueryAlignNodes **ppFirst2Rpts); // allocated to hold ptrs to alignment nodes which are marked as being FlgFirst2tRpt

	uint32_t	IdentifyHighScorePaths(uint32_t QueryLen,		// query length
			uint32_t TargSeqLen,					// targeted sequence length
			bool bStrand,						// reporting for paths on this strand - false if sense, true if antisense
			uint32_t NumNodes,					// reporting on this number of nodes starting from StartNodeIdx
			uint32_t StartNodeIdx,				// report for nodes starting at this node index (1..NumNodes) which is expected to be the first alignment node of a new target sequence
			tsQueryAlignNodes *pAlignNodes,		// alignment nodes
			uint32_t MinPathScore,				// only report those series having at least this score
			uint32_t  MaxPathsToReport);		// report at most this many alignment paths for any query

	uint32_t									// returned best score for paths starting at pAlignNodes[ExploreNodeIdx]
		HighScoreSW(uint32_t QueryLen,			// query length
			uint32_t TargSeqLen,				// targeted sequence length
 			bool bStrand,						// scoring for series on this strand - false if sense, true if antisense
			uint32_t ExploreNodeIdx,			// node to be explored for maximally scored path
			uint32_t NumNodes,					// total number of alignment nodes 
			tsQueryAlignNodes *pAlignNodes);	// alignment nodes
			
};

