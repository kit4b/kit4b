#pragma once

const int32_t cDfltAccumBinSize = 100;					// counts accumulated into this sized non-overlapping bins - bin size is in Kbp, must be multiple of 10Kbp
const int32_t cDfltNumMeanBins = 10;				// means derived over this number of sliding bins
const int32_t cMaxPBAWorkerThreads = 32;				// limiting number of worker threads to be a max of this many
const int32_t cDfltKmerSize = 100;						// default Kmer PBA alignment length for Kmer scoring
const int32_t cDfltAcceptKmerMatchPerc = 97;			// default accept as Kmer match if this percentage of base alleles match
const int32_t cMaxFounderReadsets = 16;					// can process at most this max number of founder readsets
const int32_t cMaxSkimReadsets = 100;					// can specify a max of this number of input skim file specs (each file spec can contain wildcards)
const int32_t cMaxChromNames = 100000;					// max number of unique chromsome names allowed
const int32_t cMaxChromMetadata = 1000000;				// allowing for a maximum of this many chromosomes/contigs over all founders
const int32_t cAllocChromMetadata = 100000;				// allocate chrom metadata in this sized increments
const size_t cAllocPackedBaseAlleles = 0x3fffffff;		// allocate packed base alleles to this maximal size, 1 allocation per chromosome per readset
const size_t cAllocWinBinCnts = 10000000;				// allocate windowed founder counts in this sized increments
const int32_t cDfltRawCntsBinSize = 10000;				// raw counts are for this sized bins
const uint32_t cInBuffSize = 0x7fffffff;				// buffer size for file reads
const uint32_t cOutBuffSize = 0x7fffffff;				// buffer size for file writes

// following allele scoring thresholds are extracted from KAligner.h, you must ensure that these thresholds are updated each time KAligner is updated
const double cScorePBA3MinProp = 0.75;		// score PBA as 3 if allele proportion of all counts is >= this threshold and coverage is >= 5
const double cScorePBA2MinProp = 0.35;		// score PBA as 2 if allele proportion of all counts is >= this threshold and coverage is >= 5
const double cScorePBA1MinProp = 0.20;		// score PBA as 1 if allele proportion of all counts is >= this threshold and coverage is >= 5
// when coverage is less than 5 (usually the case with skim reads) and thus confidence in alleles is reduced then scores are reduced
const double cScorePBA2MinLCProp = 0.70;		// score PBA as 2 if allele proportion of all counts is >= this threshold and coverage is < 5
const double cScorePBA1MinLCProp = 0.30;		// score PBA as 1 if allele proportion of all counts is >= this threshold and coverage is >= 5

typedef enum TAG_eModeSH {
	eMSHDefault = 0,		// default is to generate segmentations using bin counts of unique loci only
	eMSHSegAll,			// generate segmentations using bin counts of all alignment loci incl non-unique
	eMSHPlaceHolder			// used to mark end of processing modes
	}eModeSH;

#pragma pack(1)
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
bool bIsSkim;					// true if this was a skim readset, false if a founder
uint32_t NumChroms;				// readset has this number of chromosomes
uint32_t StartChromMetadataIdx;	// index of starting chrom metadata for this readset
uint32_t StartChromID;			// starting chrom for this readset
uint32_t LRAWinBinCntIdx;		// index into m_pWinBinCnts[] of LRA windowed bin counts
uint32_t StartWinBinCntsIdx;	// windowed bin counts for this founder start at m_pWinBinCnts[StartWinBinCntsIdx]
uint32_t NumWinBins;			// this number of window bins are used by this readset
} tsReadsetMetadata;

typedef struct TAG_sBinCnts
{
uint32_t ReadsetID;				// counts are for this readset
uint32_t ChromID;				// called on this chromosome
uint8_t RawHaplotypeClass;		// haplotype class raw called: 0: indeterminate or missing counts, 1: Fa, 2: Fb, 3: Fa+Fb
uint8_t HaplotypeClass;			// haplotype class after the RawHaplotypeClass has been confidence processed
uint8_t HaplotypeConf;			// level of confidence in RawHaplotypeClass, 0 low, 1 medium, 2 high - very arbitrary but derived from differences in ChiSqrs 
uint32_t StartLoci;				// window starting at this loci
uint32_t EndLoci;				// window ending at this loci inclusive
uint32_t MultiFounder;			// number of loci in window with multiple founders having same highest scores
uint32_t NoFounder;				// number of loci in window having no skim PBA
uint32_t NonUniques[cMaxFounderReadsets];	// Kmers shared - non-unique -with these other founders summed in window over StartLoci through to EndLoci
uint32_t Uniques[cMaxFounderReadsets];	// Kmers uniquely aligning to this founder summed in window over StartLoci through to EndLoci
} tsBinCnts;

typedef struct TAG_sWorkQueueEl 
{
uint32_t ReadsetID;				// readset being aligned to founders - could be a skim or founder
uint32_t ChromID;				// processing is for this chromosome
uint32_t StartLoci;				// starting from this loci inclusive
uint32_t EndLoci;				// and ending at this loci inclusive
uint32_t ChromLen;				// chromosome length
uint32_t AccumBinSize;			// counts are accumulated into this sized (bp) non-overlapping bins
uint32_t NumFndrs;				// number of founders to be processed against the skim (or founder) PBA at *pPBAs[0]
uint8_t *pPBAs[cMaxFounderReadsets]; // pPBAs[0] pts to chromosome skim PBA, followed by ptrs to each of the chromosome founder PBAs
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
	uint32_t m_NumFndrs;				// number of founders
	uint32_t m_LAReadsetNameID;			// name identifier last returned by AddReadsetName()
	uint32_t m_NumReadsetNames;			// number of readsets names currently in m_szReadsets
	uint32_t m_NxtszReadsetIdx;		// current concatenated (names separated by '\0') of all readset names in m_szReadsets
	char m_szReadsetNames[(cMaxFounderReadsets + 1) * cMaxDatasetSpeciesChrom];	// used to hold concatenated readset names, each separated by '\0', allowing for 1 skim readset
	uint32_t m_szReadsetIdx[cMaxFounderReadsets + 1];	// array of indexes into m_szFounders giving the starts of each founder name. , allowing for 1 skim readset
	tsReadsetMetadata m_Readsets[cMaxFounderReadsets + 1];	// array of all readset metadata. , allowing for 1 skim readset

	int m_NChroms;									// speeds up debugging - limit processing: 0 All chromosomes, 1..n 1st to nTh chromosome
	uint32_t m_LAChromNameID;						// last accessed chromosome identifier from call to AddChrom()
	uint32_t m_NumChromNames;						// number of chromosome names currently in m_szChromNames
	uint32_t m_NxtszChromIdx;						// current concatenated (names separated by '\0') of all chromosome names in m_szChromNames
	char m_szChromNames[cMaxChromNames * cMaxDatasetSpeciesChrom];	// used to hold concatenated chromosome names, each separated by '\0'
	uint32_t m_szChromIdx[cMaxChromNames];		// array of indexes into m_szChromNames giving the starts of each chromosome name

	uint32_t m_UsedNumChromMetadata;	// current number of chrom metadata used 
	uint32_t m_AllocdChromMetadata;		// current allocation is for this many 
	size_t m_AllocdChromMetadataMem;	// current mem allocation size for m_pChromMetadata
	tsChromMetadata *m_pChromMetadata;	// allocated to hold metadata for all founder chromosomes

	uint32_t m_UsedWinBinCnts;			// number of actually used bin counts
	uint32_t m_AllocdWinBinCnts;		// number of allocated windowed bin counts
	size_t m_AllocdWinBinCntsMem;		// current mem allocation size for m_pWinBinCnts
	tsBinCnts *m_pWinBinCnts;			// allocated to hold bin counts

	uint32_t m_TotWorkQueueEls;		// total number of workqueue elements to be processed
#ifdef WIN32
	alignas(4) volatile uint32_t m_NumQueueElsProcessed;	// number thus far processed or in progress
#else
	__attribute__((aligned(4))) volatile uint32_t m_NumQueueElsProcessed;	// number thus far processed or in progress
#endif
	tsWorkQueueEl *m_pWorkQueueEls;	// thread work queue elements - currently really a list that is iterated over but in future may be adapted to become a queue

	uint32_t m_WinMeanBins;			// derive haplotype calls from mean of counts in a sliding window containing this many bins
	bool m_bTrimmedMean;		// mean is a trimmed mean with highest and lowest bin counts in sliding window removed
	bool m_bIsMonoploid;		// false: diploid, true: monoploid
	uint32_t m_AccumBinSize;	// counts are accumulated over this sized (bp) bins
	int32_t m_ScoreKmerSize;	// Kmer PBA alignment length for Kmer scoring
	int32_t m_AcceptKmerMatchPerc; // only accept as Kmer match if at least this percentage of base alleles match

	uint32_t m_InNumProcessed;	// this number of input bytes have been processed
	uint32_t m_InNumBuffered;	// currently buffering this many input bytes
	uint32_t m_AllocInBuff;		// m_pInBuffer allocated to hold this many input bytes
	uint8_t *m_pInBuffer;		// allocated for buffering input

	uint32_t m_OutBuffIdx;		// currently buffering this many output bytes
	uint32_t m_AllocOutBuff;	// m_pOutBuffer allocated to hold this many output bytes
	uint8_t *m_pOutBuffer;		// allocated for buffering output

	int m_hInFile;				// input file handle
	int m_hOutFile;				// output file handle

	int m_hROIInFile;			// file handle for read in ROI
	int m_hROIOutFile;			// file handle for writing out haplotypes in ROI
	char *m_pszROIInFile;		// defining regions of input file
	char *m_pszROIOutFile;		// Regions of interest haplotype calls output file (CSV format)

	// some globals currently used for assessing the proportions of bins likely to contain both founders
	uint32_t m_PropExpAboveObs;	// number of bins where proportion of expected Fa+Fb for Fa and Fb combined was above the observed skim proportion of Fa+Fb 
	uint32_t m_PropExpEqualObs;	// number of bins where proportion of expected Fa+Fb for Fa and Fb combined was equal to the observed skim proportion of Fa+Fb
	uint32_t m_PropExpBelowObs;	// number of bins where proportion of expected Fa+Fb for Fa and Fb combined was below the observed skim proportion of Fa+Fb

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


	int
		LoadPBAFile(char* pszFile,bool bIsSkim = false);	// load chromosome metadata and PBA data from this file


	uint32_t		// returned chrom identifier, 0 if unable to accept this chromosome name
		AddChrom(char* pszChrom); // associate unique identifier with this chromosome name

	char*	// returns chromosome name associated with 
		LocateChrom(uint32_t ChromID); // chromosome identifier

	uint32_t		// returned chrom identifier, 0 if unable to locate this chromosome name
		LocateChrom(char* pszChrom); // return unique identifier associated with this chromosome name

	uint32_t		// returned readset name identifier, 0 if unable to accept this readset name
		AddReadset(char* pszReadset);		// associate unique identifier with this readset name, readset names must be unique including the skim readset!

	char*							// returned readset name
		LocateReadset(uint32_t ReadsetID);// identifier returned by call to AddReadsetName

	uint32_t		// returned readset name identifier, 0 if unable to locate existing readset
		LocateReadset(char* pszReadset); // associate unique identifier with this readset name

	uint8_t *								// returned pointer to start of PBA
		LocatePBAfor(uint32_t ReadSetID,	// readset identifier 
			 uint32_t ChromID);				// chrom identifier

	uint32_t								// number buffered
		FillInBuffer(uint32_t MinRequired); // try and fill input buffer with at least this many bytes reading from currently opened input file handle

	uint8_t *AllocPBAs(uint32_t ChromLen);	// allocate to hold at least this many packed base alleles

	int AllocChromMetadata(void);			// allocate for additional chromosome metadata

	tsBinCnts *							// returned ptr to allocated and initialised bin counts
		AddWinBinCnts(tsBinCnts *pInitBinCnts);	// allocated tsWinBinCnts to be initialised with a copy of pInitWinBinCnts

	int ProcessSkimPBAFile(char* pszSkimPBAFile,	// load, process and call haplotypes for this skim PBA file against previously loaded panel founder PBAs
							char *pszRsltsFileBaseName);		// results are written to this file base name with skim readset identifier and type appended

	int CountHaplotypes(uint32_t ReadsetID,		// this is the readset being aligned to founders - could be a skim or founder
							uint32_t NumFndrs,		// number of founders to be processed against
							int KmerSize,			// aligning founders using this Kmer size
							uint32_t AccumBinSize = 100000);	// comparing Readset's PBA against founder PBAs then counting most probable haplotypes using non-overlapping bins of this size

	int
		ClassifyChromHaplotypes(uint32_t ReadsetID,				// this readset is used for identifying the readset being aligned to founders - could be a skim or founder
								uint32_t ChromID,				// processing is for this chromosome
								uint32_t StartLoci,				// starting from this loci inclusive
								uint32_t EndLoci,				// and ending at this loci inclusive
								uint32_t ChromLen,				// chromosome length
								uint32_t AccumBinSize,			// accumulating counts into this sized (bp) non-overlapping bins
								uint32_t NumFndrs,				// number of founders to be processed against the skim (or founder) PBA at *pPBAs[0]
								uint8_t *pPBAs[]);				// pPBAs[0] pts to chromosome skim PBA, followed by ptrs to each of the chromosome founder PBAs

	tsBinCnts *			// returned WinBinCnts if located in pReadset
			LocateWinBinCnts(tsBinCnts *pWinBinCnts,	// chrom and start loci in this bin are to be located
								tsReadsetMetadata *pReadset);	// in this readset

	int64_t				// < 0 if errors, returned packed Kmer scores, (2 bits per founder, founder 1 in bits 0..1, founder 31 in bits 60..61) scored as being as being present at this skim loci
		ScoreSkimLoci(uint32_t ReadsetID,				// this readset is used for identifying the readset being aligned to founders - could be a skim or founder
					uint32_t Loci,				// processing for Kmer haplotypes starting from this loci
						uint32_t SeqLen,		// sequence or chromosome is this length
						uint32_t KmerSize,		// Kmer to score on is of this size
						uint32_t NumFndrs,				// number of founders to be processed against the skim (or founder) PBA at *pPBAs[0]
						uint8_t *pPBA[]);		// pts to loci 0 of PBAs in readset order

	int32_t				// returned Kmer score resulting from PBA alignment of pFndrA against pFndrB
		ScoreFounderLoci(uint32_t Loci,			// processing for Kmer haplotypes starting from this loci
						uint32_t SeqLen,		// sequence or chromosome is this length
						uint32_t KmerSize,		// Kmer to score on is of this size
						uint8_t *pAFndrPBA,		// PBA for FndrA
						uint8_t *pBFndrPBA);	// PBA for FndrB

	int ChooseHaplotypes(void); // iterate over all the window bins and choose which haplotype to call for that bin interpolating those bins which were initially called as being indeterminate

	// each returned packed haplotype occupies 4 bits, allowing for up to a total of 16 haplotype states
	//			0 - unable to call as either there are no founder and/or skim counts
	//			1 - Fa as likely the only haplotype
	//			2 - Fb as likely the only haplotype
	//			3 - Fa + Fb as likely both haplotypes present
	uint64_t						
		ChooseHaplotype(uint64_t HMMBinHaplotypes,		// 16 packed previously called haplotypes for use in as HMM transition probabilities, bits 0..3 contain previously called, through to bits 60..63 as oldest called (32 state)
						tsBinCnts* pSkimWinBinCnts);	// bin counts for skim
				 
	double ChiSqr(int NumRows, int32_t* pExp,int32_t *pObs);

	int ReportHaplotypesBED(char *pszOutFile,	// report haplotypes for specified founder to BED pszOutFile
							char *pszFounder,	// this was the founder haplotype
							char *pszSkim,		// against which this skim readset was aligned
							uint32_t HaplotypeMsk=0x01);	// reporting on haplotypes matching this mask, 0x01 for Fa, 0x02 for Fb, 0x04 for Fc, 0x010 for Fd

	CMTqsort m_mtqsort;				// multi-threaded qsort
	static int SortWinBinCnts(const void* arg1, const void* arg2);

public:
	CCallHaplotypes();
	~CCallHaplotypes();
	void Reset(void);	// resets class instance state back to that immediately following instantiation
	int Process(eModeSH PMode,		// processing mode
			int NChroms,			// debugging - limit processing: 0 All chromosomes, 1..n 1st to nTh chromosome
			char *pszTrackName,		// track name
			char *pszTrackDescr,	// track descriptor
			int KmerSize,			// aligning founders using this Kmer size
			int AcceptKmerMatchPerc, // only accept as Kmer match if this percentage of base alleles match
			uint32_t AccumBinSize,	 // accumulating counts into this sized (bp) non-overlapping bins
			int WinMeanBins,		// derive haplotype calls from mean of counts in a sliding window containing this many bins
			bool bTrimmedMean,		// mean is a trimmed mean with highest and lowest bin counts in sliding window removed
			bool bIsMonoploid,		// false: diploid, true: monoploid
			int NumFounderInputFiles,	// number of input founder file specs
			char *pszFounderInputFiles[],	// names of input founder PBA files (wildcards allowed)
			int NumSkimInputFiles,	// number of input skim file specs
			char* pszSkimInputFiles[],		// names of input skim PBA files (wildcards allowed)
			char *pszROIInFile,		// defining regions of input file
			char *pszROIOutFile,	// Regions of interest haplotype calls output file (CSV format)
			char* pszOutFile,		// output to this file
			int NumThreads);		// number of worker threads to use

	int ProcWorkerThread(tsWorkerInstance *pThreadPar);	// worker thread parameters
};


