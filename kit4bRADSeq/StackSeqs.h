#pragma once

const int cMinP1StackDepth = 5;		// user can specify P1 stack as must contain at least this many reads
const int cDfltP1StackDepth = 10;	// by default P1 stack must contain at least this many reads
const int cMaxP1StackDepth = 1000;	// user can specify at most P1 stack to contain at least this many reads

const int cDfltP1StackEnd = -1;		// by default stack end float is set as being 5% of mean read length
const int cMinP1StackEnd = 0;		// no P1 stack end float allowed
const int cMaxP1StackEnd = 20;		// max P1 stack end float allowed

const int cMinP1StackSubRate = 0;	// min user specified P1 stack substitution rate 
const int cDfltP1StackSubRate = 1;	// default P1 stack substitution rate
const int cMaxP1StackSubRate = 10;	// max user specified P1 stack substitution rate

const int cDfltP2MinOvrl = -1;		// by default P2 minimum read overlap will be set to be 50% of mean P2 read length
const int cMinP2MinOvrl = 10;		// user can specify P2 minimum read overlap down to this many bp
const int cMaxP2MinOvrl = 100;		// user can specify at most this P2 minimum read overlap 

const int cMinP2MaxOvrlSubRate = 0;		// min user specified P2 substitution rate 
const int cDfltP2MaxOvrlSubRate = 5;	// default P2 substitution rate
const int cMaxP2MaxOvrlSubRate = 10;	// max user specified P2 substitution rate


const int cAllocRawSeqLen = 50000;	// allow for P1/P2 read sequences loaded from file of up to this length

const uint32_t cMaxSfxBlkEls = 4000000000;	// construct suffix block arrays with 4byte array elements if no more than this, otherwise use 5 byte array elements

const int64_t cMaxConcatSeqLen = (int64_t)cMaxSfxBlkEls * 50; // limit concatenated read sequences length to total no more than this length   
const size_t cMinConcatSeqLen = (size_t)100000;	// always allocate at least this for holding concatenated read sequences

const int cReallocNumReads = 2000000;  // incrementally allocate/realloc for this number of reads
const int cReallocConcatSeqs = (cReallocNumReads * 100); // incrementally allocate for concatenated sequences in this sized memory increments

const int cAbsMinCoreLen = 12;	// core length is 1/3 of minimum overlap but is constrained to be in the range cAbsMinCoreLen to cAbsMaxCoreLen
const int cAbsMaxCoreLen = 50;	// always limit cores to be <= maximum core length
const int cChkIterDepth = 100;  // if more than this many matches of core then check if there are too many matching cores

const uint8_t cCSeqBOS = 0xf0;     // marks beginning of concatenated sequences, 1st sequence immediately follows
const uint8_t cCSeqSep = 0x80;	 // separators between concatenated sequences will always have bit 7 set
const uint8_t cCSeqEOS = 0xff;	 // marks end of concatenated sequences

#pragma pack(1)
typedef enum TAG_ePMode {
	ePMDefault = 0,				// default sensitivity processing mode
	ePMMoreSens,				// more sensitive - slower
	ePMUltraSens,				// ultra sensitive - much slower
	ePMLessSens,				// less sensitive - quicker
	ePMplaceholder				// used to set the enumeration range
} etPMode;


typedef struct TAG_sxRdsSfxDirEl {
	uint32_t BlockID;				// identifies (1..n) this suffix block within this file
    uint32_t SeqID;				// SeqID for initial read sfx'd in this block
	uint64_t NumSuffixEls;	    // number of suffix elements
	uint32_t ElSize;				// each element in this block is this many bytes long - currently will be either 4 or 5
	uint64_t SfxOfs;				// file offset at which suffix array starts
	uint64_t SizeOfSfx;			// total memory (and disk space) required for holding suffix array
	uint32_t NumReads;			// number of reads sfx'd in this block
	} tsxRdsSfxDirEl;

typedef struct TAG_sRdsSrcFile {
	uint32_t  NumReads;					// number of reads loaded from this source file
	uint32_t	SrcFileID;					// uniquely identifies source file (1..N)
	uint8_t   SrcFileName[_MAX_PATH];		// reads source file name
} tsRdsSrcFile;

typedef struct TAG_sRdsSfxHdr {
	uint32_t NumSrcFiles;			// number of source files from which reads were loaded
	uint64_t SumReadLens;		    // sum total of all read lengths
	uint64_t ConcatSeqLen;		// length of all concatenated sequences excluding initial cCSeqBOS and final cCSeqEOS
	tsRdsSrcFile RdsSrcFiles[cMaxSrcFiles]; // reads were loaded from these source files
	} tsRdsSfxHdr;


// concatenated sequences start with tsASeqSep followed by read sequences separated by tsASeqSep with last tsASeqSep set to 0xFFFFFFFFFF
// currently size of this sequence separator is expected to be 6 bytes
typedef struct TAG_sASeqSep {
	uint8_t  SeqSep;		// cCSeqBOS if starting 1st sequence, cCSeqSep if starting an intermediate sequence, cCSeqEOS if previous sequence was the final sequence
	uint32_t SeqIDlo;		// 40bit unique sequence identifier with bit 7 set on each byte so these can be distinguished from sequence base bytes
	uint8_t  SeqIDhi;		// unique sequence identifier requires 5 bytes
	} tsASeqSep;

// suffixed sequence (initially reads but could be the merge of multiple reads as contigs are assembled)
typedef struct TAG_tsSfxdSeq {
	uint32_t SeqID;		// unique sequence identifier
	uint64_t ConcatSeqOfs; // offset in concatenated sequences (m_pConcatSeqs) at which this sequence starts
	uint32_t ReadLen;		// length of this sequence
	uint8_t Flags;		// holds various combinations of processing state flags for this sequence CAUTION: if changed larger than 1byte then need to synchronise access
} tsSfxdSeq;

typedef struct TAG_sArtefactMap {			// used to record best mapping for probe sequence on to a target whilst attempting to determine consensus basses
	uint32_t   TargSeqID;						// mapping was on to this target
	tsSfxdSeq *pTargSfxdSeq;				// this is the target
	etSeqBase *pTargLeftStart;				// mapping starts from this leftmost base target 
	int ProbeLeftOfs;						// and this leftmost 5' probe base
	int TargLeftOfs;						// and this leftmost 5' target base
	int CurOverlayLen;						// overlay is of this length
	int CurMMCnt;							// and had this many mismatches
} tsArtefactMap;


typedef struct TAG_sSeqStarts {
	uint32_t PairID;					// uniquely identifies this pair
	uint16_t P1SeqLen;				// P1 sequence length
	uint16_t P2SeqLen;				// P2 sequence length
	size_t P1SeqOfs;				// P1 sequence starts at this offset + 1 in m_pP1Seqs2Assemb
	size_t P2SeqOfs;				// P2 sequence starts at this offset + 1 in m_pP2Seqs2Assemb
} tsSeqStarts;

// tsSfxdSeq Flags
const uint8_t cFlagNA = 0x02;			// this read is not to be processed for assignment to any contig
const uint8_t cFlagMerged = 0x04;	    // this sequence has been merged with another sequence
const uint8_t cFlagDup = 0x08;		// sequence identified as being a duplicate of another


#pragma pack()

class CStackSeqs
{

	bool m_bIsPairedEndProc;	// if true then processing is for paired end P1 plus P2, if false then single ended P1 processing
	int m_NumThreads;			// max number of processing threads to use
	CMTqsort m_MTqsort;			// multithreaded sorting

	int m_hOutCtgsFile;			// contigs output file
	int m_hOutVCFfile;			// VCF output file

	uint32_t m_AllocdPathIDs;		// m_pPathIDs currently allocated to hold this many path identifiers
	uint32_t m_NumPathIDs;		// number of path identifiers currently in m_pPathIDs
	uint32_t *m_pPathIDs;			// allocated to hold path identifiers

	int m_NumRawFiles;			// number of raw reads files processed
	int m_NumRdsFiles;			// number of preprocessed (kangar) reads files processed

	int m_hInFile;				// input file handle
	char m_szInFile[_MAX_PATH];	// input file

	char m_szCtgDescr[80];		// contig descriptor prefix


	uint32_t m_TotSeqsParsed;		// total number of sequences parsed
    uint32_t m_TotSeqsUnderLen;	// total number of sequences filtered out because underlength
	uint32_t m_TotSeqsExcessNs;	// total number of sequences filtered out because too many Ns

	uint32_t m_TotP1Seqs2Assemb;	// original number of sequences to be assembled after length filtering
	uint32_t m_NumP1Seqs2Assemb;	// number of sequences in m_pP1Seqs2Assemb
	size_t m_P1Seqs2AssembLen;	// current length of sequences in m_pP1Seqs2Assemb
	uint32_t m_TotP2Seqs2Assemb;	// original number of sequences to be assembled after length filtering
	uint32_t m_NumP2Seqs2Assemb;	// number of sequences in m_pP2Seqs2Assemb
	size_t m_P2Seqs2AssembLen;	// current length of sequences in m_pP2Seqs2Assemb

	size_t m_AllocMemP1Seqs2Assemb;	// memory currently allocated for P1 m_pP1Seqs2Assemb 
	size_t m_AllocMemP2Seqs2Assemb;	// memory currently allocated for P2 m_pP2Seqs2Assemb 
	uint8_t *m_pP1Seqs2Assemb;		// holds P1 sequences (always in basespace) to assemble into contigs, concatenated with uint32_t lengths prepended to each sequence
	uint8_t *m_pP2Seqs2Assemb;		// holds P2 sequences (always in basespace) concatenated with uint32_t lengths prepended to each sequence



	uint32_t m_AllocSeqStarts;			// number of allocated read sequence starts
	size_t	m_AllocMemSeqStarts;		// memory allocated for holding the sequence starts
	uint32_t m_NumSeqStarts;				// number of sequence starts
	tsSeqStarts *m_pSeqStarts;			// pts to array of sequence starts


	int m_MeanReadLen;			// mean length of all reads

	int m_SfxElSize;			// suffix element size - either 4, <= cMaxSfxBlkEls, or 5 if > cMaxSfxBlkEls 
	int64_t m_NumSuffixEls;		// number of elements in suffix array
	int64_t m_AllocMemSfx;	    // allocated memory size for suffix array
	uint32_t *m_pSuffixArray;     // to hold suffix array for concatenated read/contig sequences

	size_t m_CurMaxMemWorkSetBytes;     // currently set max working set in bytes, need to convert to pages when setting working set
	uint32_t m_WinPageSize;				// windows memory page size in bytes (0 if process not on Windows)
	size_t m_BaseWinMinMem;				// windows base min working set memory in bytes when process initially started
	size_t m_BaseWinMaxMem;				// windows base max working set memory in bytes when process initially started

	teBSFrsltCodes AddSeq(int P1SeqLen,		// P1 sequence length
						uint8_t *pP1Seq,		// ptr to P1 sequence
						int P2SeqLen,		// P2 sequence length
						uint8_t *pP2Seq);		// ptr to P1 sequence

	tsRdsSfxHdr m_P1RdsSfxHdr;	// concatenated P1 reads header
	tsRdsSfxHdr m_P2RdsSfxHdr;	// concatenated P2 reads header

	teBSFrsltCodes								// returns number of concatenated sequences or if < 0 then error code
		GenSfxdSeqs(void);			// generates concatenated sequences into m_pConcatSeqs with m_pSfxdSeqs holding offsets into m_pConcatSeqs ready for sorting

	teBSFrsltCodes GenP1RdsSfx(void);		// generate suffix array over concatenated P1 read sequences as either 4 or 5byte sfx els

	// ChunkedRead
	// Seeks to specified 64bit file offset and reads from disk as chunks of no more than INT_MAX/32
	teBSFrsltCodes ChunkedRead(int64_t RdOfs,uint8_t *pData,int64_t RdLen);
	teBSFrsltCodes ChunkedRead(int hFile,char *pszFile,int64_t RdOfs,uint8_t *pData,int64_t RdLen);

	static inline uint64_t		// unpacks the 5 bytes ptd at by p5Bytes and returns as uint64_t 
		Unpack5(uint8_t *p5Bytes)
		{
		// ensures that only 5 bytes are actually accessed, can't just access as a uint64_t and mask retain bottom 40 bits...
		return((uint64_t)((uint64_t)p5Bytes[4] << 32 | *(uint32_t *)p5Bytes));
		}

	static inline uint8_t *Pack5(uint64_t Value, uint8_t *p5Bytes)
		{
		*(uint32_t *)p5Bytes = (uint32_t)Value;
		p5Bytes[4] = (uint8_t)(Value >> 32);
		return(p5Bytes + 5);
		}

 	static int SfxSortFunc(const void *arg1, const void *arg2);
	static int Sfx5SortFunc(const void *arg1, const void *arg2);
	static int ReadsSortFunc(const void *arg1, const void *arg2);

	CStopWatch m_StopWatch;

	// serialisations and locks
	bool m_bMutexesCreated;						// set true if mutexes and rwlocks created/initialised
	void AcquireSerialise(void);
	void ReleaseSerialise(void);
	void AcquireLock(bool bExclusive);				// defaults as read only lock
	void ReleaseLock(bool bExclusive);
#ifdef _WIN32
	HANDLE m_hMtxIterReads;
	HANDLE m_hMtxMHReads;
	SRWLOCK m_hRwLock;
	static unsigned __stdcall ThreadedReadsAssemb(void * pThreadPars);
#else
	pthread_mutex_t m_hMtxIterReads;
	pthread_mutex_t m_hMtxMHReads;
	pthread_rwlock_t m_hRwLock;
	static void *ThreadedReadsAssemb(void * pThreadPars);
#endif

	int CreateMutexes(void);
	void DeleteMutexes(void);

	bool SetMaxMemWorkSetSize(size_t Bytes);

	teBSFrsltCodes AllocReadsMemory(size_t P1ReqAllocSize,size_t P2ReqAllocSize); // allocate memory as may be required

public:
	CStackSeqs(void);
	~CStackSeqs(void);

	void Reset(bool bSync = true);
	int Init(void);

	// Transforms 31bit identifier into a 40bit identifier such that the MSB of each byte is set to 1
	// Objective is to allow the identifier to be easily distinguished from ordinary sequence data which has the MSB reset
	static inline uint64_t IDtoXForm(uint32_t ID)
		{
		uint64_t XFormID;
		if(ID > 0x07fffffff)
			return( 0x0ffffffff);
		XFormID = 0x8080808080 | ((ID & 0x07f) | (ID & 0x3f80) << 1 | (ID & 0x01fc000) << 2 | (ID & 0x0fe00000) << 3 | (ID & 0x07f0000000) << 4);  
		return(XFormID);
		}

	// XFormToID
	// Transforms a transformed 40bit identifier (generated by IDtoXForm()) back into it's original 31bit identifier
	static inline uint32_t XFormToID(uint64_t XFormID)
		{
		uint32_t ID;
		ID = (uint32_t)((XFormID & 0x07f) | (XFormID & 0x07f00) >> 1 | (XFormID & 0x07f0000) >> 2 | (XFormID & 0x07f000000) >> 3 | (XFormID & 0x07f00000000) >> 4);  
		return(ID);
		}


	void SetCtgDescr(char *pszCtgDescr);			// set contig descriptor prefix
	void SetPairedEndProc(bool bIsPairedEndProc);	// set processing is for paired end

	void SetNumThreads(int maxThreads);

	teBSFrsltCodes
		LoadRawReads(int MaxNs,				// filter out input sequences having higher than this number of indeterminate bases per 100bp (default is 1, range 0..10)
					int Trim5,				// trim this number of 5' bases from input sequences (default is 0, range 0..20)
					int Trim3,				// trim this number of 3' bases from input sequences (default is 0, range 0..20)
					int MinSeqLen,		    // filter out input sequences (after any trimming) which are less than this length (default is 50bp, range 30..10000)
					int FileID,				// uniquely identifies source file P1 and P2
					char *pszP1File,		// process from this P1 file
					char *pszP2File);		// optionally process from this P2 file


	teBSFrsltCodes 
		Process(etPMode PMode,						// processing sensitivity mode
				int P1StackEnd,						// P1 stack maximum end float (default is 1, range 0..10)");
				int P1StackDepth,					// P1 stack minimum depth (default is 10, range 5..1000)");
				int P1StackSubRate,					// P1 stack maximum substitution rate (default is 1%, range 0..10%)");
				int P2MinOvrl,						// P2 minimum read overlap (default is 30% of read length, range 10..100)");
				int P2MaxOvrlSubRate,				// P2 maximum read overlap substitution rate (default is 5%%, range 0..10%%)");
				char *pszCtgDescr,					// generated contig descriptor prefix 
				char *pszOutCtgsFile,				// assembled contigs written to this file
				char *pszOutVCFfile,				// output Variant Call Format (VCF 4.1) to this file
				int NumThreads);					// max number of worker threads to use

};

