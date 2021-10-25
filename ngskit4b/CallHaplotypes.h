#pragma once

const int32_t cMaxPBAWorkerThreads = 64;			// limiting number of worker threads to be a max of this many

const int32_t cMaxFounderFileSpecs = 100;			// can accept at most this max number of input wildcarded founder file specs
const int32_t cMaxProgenyFileSpecs = 500;			// can specify a max of this number of input wildcarded progeny file specs

const int32_t cMaxFounderReadsets = 1000;			// after input wildcard expansion then can accept at most this max number of founder readsets
const int32_t cMaxProgenyReadsets = 5000;			// after input wildcard expansion then can specify at most this max number of progeny readsets

const int32_t cMaxChromNames = 100000;				// max number of unique chromsome names allowed
const int32_t cMaxChromMetadata = 1000000;			// allowing for a maximum of this many chromosomes/contigs over all founders
const int32_t cAllocChromMetadata = 100000;			// allocate chrom metadata in this sized increments
const size_t cAllocPackedBaseAlleles = 0x3fffffff;	// allocate packed base alleles to this maximal size, 1 allocation per chromosome per readset
const uint32_t cInBuffSize = 0x7fffffff;			// buffer size for file reads
const uint32_t cOutBuffSize = 0x6fffffff;			// buffer size for file writes
const uint32_t cDfltFndrTrim5 = 10;					// trim this many aligned PBAs from 5' end of founder aligned segments - reduces false alleles due to sequencing errors
const uint32_t cDfltFndrTrim3 = 10;					// trim this many aligned PBAs from 3' end of founder aligned segments - reduces false alleles due to sequencing errors
const uint32_t cDfltProgTrim5 = 5;					// trim this many aligned PBAs from 5' end of progeny aligned segments - reduces false alleles due to sequencing errors
const uint32_t cDfltProgTrim3 = 5;						// trim this many aligned PBAs from 3' end of progeny aligned segments - reduces false alleles due to sequencing errors

// following allele scoring thresholds are extracted from KAligner.h, you must ensure that these thresholds are updated each time KAligner is updated
const double cScorePBA3MinProp = 0.75;		// score PBA as 3 if allele proportion of all counts is >= this threshold and coverage is >= 5
const double cScorePBA2MinProp = 0.35;		// score PBA as 2 if allele proportion of all counts is >= this threshold and coverage is >= 5
const double cScorePBA1MinProp = 0.20;		// score PBA as 1 if allele proportion of all counts is >= this threshold and coverage is >= 5
// when coverage is less than 5 (usually the case with progeny reads) and thus confidence in alleles is reduced then scores are reduced
const double cScorePBA2MinLCProp = 0.70;		// score PBA as 2 if allele proportion of all counts is >= this threshold and coverage is < 5
const double cScorePBA1MinLCProp = 0.30;		// score PBA as 1 if allele proportion of all counts is >= this threshold and coverage is >= 5

const uint32_t cAllocProgenyFndrAligns = 1000000; // initial allocation for progeny to founder allele stacks alignments

const uint32_t cAllocAlleleStacks = 1000000;		// initial allocation for allele stacks - number of stacks
const uint32_t cReallocAlleleStacks = (cAllocAlleleStacks/2);	// realloc allele stacks with this sized increments - number of stacks


typedef enum TAG_eModeCSH {
	eMCSHDefault = 0,		// default is to process founder allele stacks
	eMCSHPlaceHolder		// used to mark end of processing modes
	}eModeCSH;

typedef enum TAG_eHapClass	// confidence in a haplotype for given founder being present or absent within a given bin
{
	eHCAbsent = 0,		// high confidence that founder has no haplotype
	eHCLoConf,			// low confidence that founder has a haplotype
	eHCHiConf			// high confidence that founder has a haplotype
} eHapClass;

#pragma pack(1)

typedef struct TAG_s256Bits {
	uint8_t NumBits;		// total number of bits over all 64bit words in Bits256
	uint64_t Bits[4];	// 256bits encoded into 4 x 64bit words
} ts256Bits;


typedef struct TAG_sAlleleStack {	// stacked informative - all dirac alleles and not all alleles the same in stack -, e.g alleles at the given chrom and loci
	uint32_t ChromID;		// stack is on this chromosome
	uint32_t Loci;			// and at this loci
	uint8_t NumFndrs;		// number of founders
	uint8_t NumProcFndrs;	// number of founders actually processed - some founders may be marked as not for processing
	ts256Bits ProcFndrs;	// bitmap of founders actually processed
	uint8_t NumAlleleFndrs[4];	// number of founders mapped in corresponding Alleles[] 
	ts256Bits Alleles[4];	// bitmaps (founderA in bit0) for each founder allele contained in this stack - max of 256 founders per allele 
} tsAlleleStack;

typedef struct TAG_sProgenyFndrAligns {
	uint32_t ReadsetID;				// identifies the progeny readset
	uint32_t ChromID;				// stack is on this chromosome
	uint32_t Loci;					// and at this loci
	uint8_t Alleles;				// progeny has these alleles present at the allelestack chrom.loci
	uint8_t NumProgenyFounders;		// number of potential progeny founders in ProgenyParents bitmap
	ts256Bits ProgenyFounders;		// bitmap of potential progeny founders - progeny has a minor/major allele shared with corresponding founder
} tsProgenyFndrAligns;

typedef struct TAG_sChromMetadata
{
uint32_t ReadsetID;			// chromosome is from this readset
uint32_t ChromMetadataIdx;	// index of this metadata
uint32_t NxtChromMetadataIdx; // chromosome metadata in a given readset are forward linked
uint32_t ChromID;			// chromosome identifier
uint32_t ChromLen;			// has this many loci bases
uint8_t *pPBAs;				// pts to memory allocation holding packed base alleles for this chromosome
} tsChromMetadata;

typedef struct TAG_sReadsetMetadata
{
uint32_t ReadsetID;				// identifies from which readset these packed base alleles were generated
char szExperimentID[100];		// sequencing experiment
char szRefAssemblyID[100];		// alignments were against this target assembly
uint8_t ReadsetType;			// 0: founder, 1: progeny, 2: control
uint32_t NumChroms;				// readset has this number of chromosomes
uint32_t StartChromMetadataIdx;	// index of starting chrom metadata for this readset
uint32_t StartChromID;			// starting chrom for this readset
} tsReadsetMetadata;

typedef struct TAG_sWorkQueueEl 
{
uint32_t ChromID;				// processing is for this chromosome
uint32_t StartLoci;				// starting from this loci inclusive
uint32_t ChromLen;				// chromosome length
uint32_t MaxNumLoci;			// max number of loci to be processed - starting from StartLoci - to be processed in this queue element
uint32_t NumFndrs;				// number of founders to be processed against the progeny (or founder) PBA at *pPBAs[0]
uint8_t *pFounderPBAs[cMaxFounderReadsets]; // ptrs to each of the chromosome founder PBAs indexed in FounderID ascending order
uint8_t *pMskPBA;				// pts to optional chromosome masking PBA, scoring only for segments contained in mask which are non-zero and where progeny and founder PBAs also have an allele.
								// enables a founder to be processed as if a progeny, restricting scoring Kmers to be same as if a progeny
} tsWorkQueueEl;

typedef struct TAG_sWorkerInstance {
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
} tsWorkerInstance;

#pragma pack()

class CCallHaplotypes
{
	CStats m_Stats;					// used for determining statistical significance of individual founder unique counts relative to other founder counts

	uint32_t m_FndrTrim5;			// trim this many aligned PBAs from 5' end of founder aligned segments - reduces false alleles due to sequencing errors
	uint32_t m_FndrTrim3;			// trim this many aligned PBAs from 3' end of founder aligned segments - reduces false alleles due to sequencing errors
	uint32_t m_ProgTrim5;			// trim this many aligned PBAs from 5' end of progeny aligned segments - reduces false alleles due to sequencing errors
	uint32_t m_ProgTrim3;			// trim this many aligned PBAs from 3' end of progeny aligned segments - reduces false alleles due to sequencing errors

	uint64_t m_FndrsProcMap;			// bitmap (founder1 in bit position 0) corresponding to founders which are to be processed, if bit is not set then the corresponding founder is not to be processed
	uint32_t m_NumFndrsProcMap;			// number of founders in m_FndrsProcMap
	uint32_t m_NumFounders;				// number of founders
	uint32_t m_NumProgenies;			// number of progenies
	uint32_t m_MaskReadsetID;			// masking readset identifier, 0 if no masking readset
	int32_t m_CurProgenyReadsetID;		// current progeny readset identifier
	int32_t m_ProgenyIDs[cMaxProgenyReadsets];	// array of progeny readset identifiers in order of loading
	int32_t m_LAReadsetNameID;			// name identifier last returned by AddReadsetName()
	int32_t m_NumReadsetNames;			// number of readsets names currently in m_szReadsets
	uint32_t m_NxtszReadsetIdx;			// current concatenated (names separated by '\0') of all readset names in m_szReadsets, note that the readset type (founder 0, progeny 1 or control 2) is prepended to each name
	char m_szReadsetNames[(cMaxFounderReadsets + cMaxProgenyReadsets + 1) * (cMaxDatasetSpeciesChrom + 1)];	// used to hold concatenated readset names, each separated by '\0', allowing for 1 progeny readset
	uint32_t m_szReadsetIdx[cMaxFounderReadsets + cMaxProgenyReadsets + 1];	// array of indexes into m_szReadsetNames giving the starts of each founder or progeny  name. , allowing for a control readset name
	tsReadsetMetadata m_Readsets[cMaxFounderReadsets + cMaxProgenyReadsets + 1];	// array of all readset metadata

	int32_t m_LAChromNameID;						// last accessed chromosome identifier from call to AddChrom()
	int32_t m_NumChromNames;						// number of chromosome names currently in m_szChromNames
	int32_t m_NxtszChromIdx;						// current concatenated (names separated by '\0') of all chromosome names in m_szChromNames
	char m_szChromNames[cMaxChromNames * cMaxDatasetSpeciesChrom];	// used to hold concatenated chromosome names, each separated by '\0'
	uint32_t m_szChromIdx[cMaxChromNames];		// array of indexes into m_szChromNames giving the starts of each chromosome name

	uint32_t m_UsedNumChromMetadata;	// current number of chrom metadata used 
	uint32_t m_AllocdChromMetadata;		// current allocation is for this many 
	size_t m_AllocdChromMetadataMem;	// current mem allocation size for m_pChromMetadata
	tsChromMetadata *m_pChromMetadata;	// allocated to hold metadata for all founder chromosomes

	uint32_t m_UsedProgenyFndrAligns;			// number of actually used progeny to founder alignments
	uint32_t m_AllocdProgenyFndrAligns;			// number of allocated founder alignments
	size_t m_AllocdProgenyFndrAlignsMem;		// current mem allocation size for m_pProgenyFndrAligns
	tsProgenyFndrAligns*m_pProgenyFndrAligns;	// allocated to hold progeny allele stack alignments

	uint32_t m_UsedAlleleStacks;		// number of actually used allele stacks
	uint32_t m_AllocdAlleleStacks;		// number of allocated allele stacks
	size_t m_AllocdAlleleStacksMem;		// current mem allocation size for m_pAlleleStacks
	tsAlleleStack *m_pAlleleStacks;		// allocated to hold allele stacks

	uint32_t m_AllocWorkQueueEls;		// work queue allocated to hold at most this number of elements
	uint32_t m_TotWorkQueueEls;			// total number of work queue elements to be processed
#ifdef WIN32
	alignas(4) volatile uint32_t m_NumQueueElsProcessed;	// number thus far processed or in progress
#else
	__attribute__((aligned(4))) volatile uint32_t m_NumQueueElsProcessed;	// number thus far processed or in progress
#endif
	tsWorkQueueEl *m_pWorkQueueEls;	// thread work queue elements - currently really a list that is iterated over but in future may be adapted to become a queue

	uint8_t m_Fndrs2Proc[cMaxFounderReadsets];	// array of founders which are to be processed, indexed by FounderID-1. If LSB is set then that founder is marked for processing

	uint32_t m_InNumProcessed;	// this number of input bytes have been processed
	uint32_t m_InNumBuffered;	// currently buffering this many input bytes
	uint32_t m_AllocInBuff;		// m_pInBuffer allocated to hold this many input bytes
	uint8_t *m_pInBuffer;		// allocated for buffering input

	uint32_t m_OutBuffIdx;		// currently buffering this many output bytes
	uint32_t m_AllocOutBuff;	// m_pszOutBuffer allocated to hold this many output bytes
	uint8_t *m_pszOutBuffer;	// allocated for buffering output

	int m_hInFile;				// input file handle
	int m_hOutFile;				// file handle for writing results

	char *m_pszMaskBPAFile;		// optional input masking BPA file, only process BPAs which are intersect of these BPAs and progeny plus founder BPAs

#ifdef WIN32
	alignas(4)	volatile unsigned int m_ReqTerminate; // used with synchronous compare and swap (CAS) for serialising access to thread early terminate request
	alignas(4)	volatile unsigned int m_NumThreads; // used with synchronous compare and swap (CAS) for serialising access to number of worker threads to use
	alignas(4)	volatile unsigned int m_NumWorkerInsts; // used with synchronous compare and swap (CAS) for serialising access to actual number of worker instances
	alignas(4)	volatile unsigned int m_CompletedWorkerInsts; // used with synchronous compare and swap (CAS) for serialising access to  number of completed worker instances
#else
	__attribute__((aligned(4))) volatile unsigned int m_ReqTerminate; // used with synchronous compare and swap (CAS) for serialising access to thread early terminate request
	__attribute__((aligned(4))) volatile unsigned int m_NumThreads; // used with synchronous compare and swap (CAS) for serialising access to number of worker threads to use
	__attribute__((aligned(4))) volatile unsigned int m_NumWorkerInsts; // used with synchronous compare and swap (CAS) for serialising access to actual number of worker instances
	__attribute__((aligned(4))) volatile unsigned int m_CompletedWorkerInsts; // used with synchronous compare and swap (CAS) for serialising access to  number of completed worker instances
#endif

	tsWorkerInstance m_WorkerInstances[cMaxPBAWorkerThreads];	// to hold all worker instance thread parameters

	bool m_bMutexesCreated;		// will be set true if synchronisation mutexes have been created
	int CreateMutexes(void);
	void DeleteMutexes(void);

#ifdef _WIN32
	HANDLE m_hSerialiseAccess;
#else
	pthread_mutex_t m_hSerialiseAccess;
#endif
	void AcquireSerialise(void);
	void ReleaseSerialise(void);

	// initialise and start pool of worker threads
	int	StartWorkerThreads(uint32_t NumThreads,		// there are this many threads in pool
							uint32_t NumChroms);		// processing this number of chromosomes

	bool WaitAlignments(int WaitSecs);	// allow at most this many seconds for pool of worker threads to complete PBA scoring

	int	TerminateWorkerThreads(int WaitSecs = 120);	// allow at most this many seconds before force terminating threads


	int32_t					// returned readset identifier (1..n) or < 1 if errors
		LoadPBAFile(char* pszFile,		// load chromosome metadata and PBA data from this file
				uint8_t ReadsetType);	// 0: founder, 1: progeny, 2: control

	// trim 5' and 3' aligned segments within the PBAs attempting to reduce sequencing error induced false alleles
	int			// returns number of non-trimmed loci in the pPBAs
		TrimPBAs(uint32_t Trim5,	// trim 5' this many aligned PBA bases from each aligned segment
								  uint32_t Trim3,	// trim 3' this many aligned PBA bases from each aligned segment
								  uint32_t PBALen,	// pPBAs contains this many packed base alleles
								  uint8_t* pPBAs);	// packed base alleles to be processed

	uint8_t LocateReadsetChromLociAlleles(int32_t ReadsetID,	// return alleles for this readset 
										int32_t ChromID,		// on this chromosome
										uint32_t Loci);			// at this loci

	int32_t		// returned chrom identifier, < 1 if unable to accept this chromosome name
		AddChrom(char* pszChrom); // associate unique identifier with this chromosome name

	char*	// returns chromosome name associated with 
		LocateChrom(int32_t ChromID); // chromosome identifier

	int32_t		// returned chrom identifier, < 1 if unable to locate this chromosome name
		LocateChrom(char* pszChrom); // return unique identifier associated with this chromosome name

	int32_t		// returned readset name identifier, < 1 if unable to accept this readset name
		AddReadset(char* pszReadset,		// associate unique identifier with this readset name, readset names must be unique within the readset type
				   uint8_t ReadsetType);	// 0: founder, 1: progeny, 2: control

	char*							// returned readset name
		LocateReadset(int32_t ReadsetID);	// identifier returned by call to AddReadsetName

	int32_t		// returned readset name identifier, 0 if unable to locate existing readset
		LocateReadset(char* pszReadset, // associate unique identifier with this readset name
					  uint8_t ReadsetType);	// 0: founder, 1: progeny, 2: control

	uint8_t *								// returned pointer to start of PBA
		LocatePBAfor(int32_t ReadSetID,	// readset identifier 
			 int32_t ChromID);				// chrom identifier

	uint32_t								// number buffered
		FillInBuffer(uint32_t MinRequired); // try and fill input buffer with at least this many bytes reading from currently opened input file handle

	uint8_t *AllocPBAs(uint32_t ChromLen);	// allocate to hold at least this many packed base alleles

	int AllocChromMetadata(void);			// allocate for additional chromosome metadata

	uint32_t									// returned index+1 into  m_pProgenyFndrAligns[] to allocated and initialised ProgenyFndrAligns, 0 if errors
		AddProgenyFndrAligns(tsProgenyFndrAligns* pInitProgenyFndrAligns);	// allocated tsProgenyFndrAligns to be initialised with a copy of pInitProgenyFndrAligns


	uint32_t									// returned index+1 into m_pAlleleStacks[] to allocated and initialised allele stack, 0 if errors
		AddAlleleStack(tsAlleleStack* pInitAlleleStack);	// allocated tsAlleleStack to be initialised with a copy of pInitAlleleStack

	int ProcessProgenyPBAFile(char* pszProgenyPBAFile,	// load, process and call haplotypes for this progeny PBA file against previously loaded panel founder PBAs
							char *pszRsltsFileBaseName);		// results are written to this file base name with progeny readset identifier and type appended

	int				// < 0 if errors, otherwise success
		GenChromAlleleStacks(uint32_t ChromID,	// processing is for this chromosome
												  uint32_t ChromSize,		// chromosome is this size
												  uint32_t Loci,			// processing for allele stacks starting from this loci
												  uint32_t MaxNumLoci,		// processing for this maximum number of loci
												  uint32_t NumFndrs,		// number of founders to be processed 1st PBA at *pPBAs[0]
												  uint8_t* pFounderPBAs[],	// pPBAs[] pts to chromosome PBAs, for each of the chromosome founder PBAs
												  uint8_t* pMskPBA);			// pts to optional masking PBA, scoring only for segments contained in mask which are non-zero and where progeny and founder PBAs also have an allele.
																			// enables a founder to be processed as if a progeny, restricting scoring Kmers to be same as if a progeny

	int
		AlignAlleleStacks(uint32_t NumFndrs,			// number of founders to be processed against
					uint32_t MaxNumLoci = 100000);		// comparing ReadsetID's PBA against founder PBAs then work queue items specify this number of loci for processing by threads

	int GenAlleleStacks(uint32_t NumFndrs,		// number of founders to be processed against
						uint32_t MaxNumLoci = 100000);		// work queue items specify this number of loci for processing by threads

	double ChiSqr(int NumRows, int32_t* pExp,int32_t *pObs); // currently not used!

	int	ReportMatrix(char* pszRsltsFileBaseName);		// matrix results are written to this file base name with '.SummaryMatrix.csv' appended

	int ReportHaplotypesByProgeny(char* pszRsltsFileBaseName,		// haplotype results are written to this file base name with 'ProgenyID.csv' appended
								  uint32_t ReadsetID = 0);				// report on this progeny readset only, or if 0 then report on all progeny readsets
							

	void	Bits256Shl1(ts256Bits& Bits256);
	void	Bits256Set(uint8_t Bit,		// bit to set, range 0..255
				   ts256Bits& Bits256);
	void	Bits256Reset(uint8_t Bit,		// bit to reset, range 0..255
					 ts256Bits& Bits256);
	bool Bits256Test(uint8_t Bit,		// bit to test, range 0..255
					ts256Bits& Bits256);
	void Bits256Initialise(bool Set,			// if true then initialse all bits as set, otherwise initialise all bits as reset
										   ts256Bits& Bits256);
	uint32_t Bits256Count(ts256Bits& Bits256);		// count number of set bits

	uint32_t Bits256Combine(ts256Bits& Bits256A, ts256Bits& Bits256B);		// combine (effective |= ) bits in Bits256A with Bits256B with Bits256A updated 

	CMTqsort m_mtqsort;				// multi-threaded qsort
	static int SortAlleleStacks(const void* arg1, const void* arg2);
	static int SortProgenyFndrAligns(const void* arg1, const void* arg2);
	static int SortProgenyChromLociReadset(const void* arg1, const void* arg2);

public:
	CCallHaplotypes();
	~CCallHaplotypes();
	void Reset(void);	// resets class instance state back to that immediately following instantiation
	int Process(eModeCSH PMode,	// processing mode 0: call observed haplotype counts, mode 1: GBS processing
				char* pszTrackName,		// track name
				char* pszTrackDescr,	// track descriptor
				uint32_t FndrTrim5,			// trim this many aligned PBAs from 5' end of founder aligned segments - reduces false alleles due to sequencing errors
				uint32_t FndrTrim3,			// trim this many aligned PBAs from 3' end of founder aligned segments - reduces false alleles due to sequencing errors
				uint32_t ProgTrim5,			// trim this many aligned PBAs from 5' end of progeny aligned segments - reduces false alleles due to sequencing errors
				uint32_t ProgTrim3,			// trim this many aligned PBAs from 3' end of progeny aligned segments - reduces false alleles due to sequencing errors
				char* pszMaskBPAFile,	// optional masking input BPA file, only process BPAs which are intersect of these BPAs and progeny plus founder BPAs
				int NumFounderInputFiles,	// number of input founder file specs
				char* pszFounderInputFiles[],	// names of input founder PBA files (wildcards allowed)
				int NumProgenyInputFiles,	// number of input progeny file specs
				char* pszProgenyInputFiles[],		// names of input progeny PBA files (wildcards allowed)
				char* pszOutFile,		// Windowed haplotype calls output file (CSV format)
				int NumThreads);		// number of worker threads to use

	int ProcWorkerThread(tsWorkerInstance *pThreadPar);	// worker thread parameters
};


