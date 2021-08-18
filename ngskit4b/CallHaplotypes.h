#pragma once

const int32_t cDfltMinAlignedScoringPBA = 2;		// to be accepted as aligned and scored then PBAs in both founder and skim must be at least this PBA level
const int32_t cDfltMinAlignedPBA = 1;				// to be accepted as aligned and not penalised, PBAs in both founder and skim must be at least this PBA level
const int32_t cDfltAlignedPBAScore = 10;			// if PBA accepted for scoring (>= m_MinAlignedScoringPBA) then add this score to the summed Kmer alignment score
const int32_t cDfltPenaltyPBAScore = (-3 * cDfltAlignedPBAScore);	// if PBA not accepted as aligned ( < m_MinAlignedPBA) then apply this penalty to the summed Kmer alignment score
const int32_t cDfltMinKmerScoreDiff = 2;			// must be at least this Kmer percentage score differential between highest scoring founders before accepting a founder as being a uniquely highest scoring

const int32_t cDfltAccumBinSize = 100;					// counts accumulated into this sized non-overlapping bins - bin size is in Kbp
const int32_t cDfltNumMeanBins = 20;					// means derived over this number of sliding bins
const int32_t cMaxPBAWorkerThreads = 64;				// limiting number of worker threads to be a max of this many
const int32_t cMinKmerSize = 20;						// allow Kmers down to this size
const int32_t cDfltMinKmerSize = 50;					// default minimum Kmer PBA alignment length for Kmer scoring
const int32_t cMaxKmerSize = 2000;						// maximum Kmer PBA alignment length for Kmer scoring, default is for 10x specified minimum kmer size 

const int32_t cDfltAcceptKmerMatchPerc = 97;			// default accept as Kmer match if this percentage of base alleles match
const int32_t cMaxFounderReadsets = 500;				// can accept at most this max number of founder readsets in the pool, this limit is rather arbitrary!
const int32_t cMaxSkimReadsets = 1000;					// can specify a max of this number of input skim file specs (each file spec can contain wildcards)
const int32_t cMaxChromNames = 100000;					// max number of unique chromsome names allowed
const int32_t cMaxChromMetadata = 1000000;				// allowing for a maximum of this many chromosomes/contigs over all founders
const int32_t cAllocChromMetadata = 100000;				// allocate chrom metadata in this sized increments
const size_t cAllocPackedBaseAlleles = 0x3fffffff;		// allocate packed base alleles to this maximal size, 1 allocation per chromosome per readset
const size_t cAllocWinBinCnts = 1000000;				// allocate windowed founder counts in this sized increments
const int32_t cDfltRawCntsBinSize = 10000;				// raw counts are for this sized bins
const uint32_t cInBuffSize = 0x7fffffff;				// buffer size for file reads
const uint32_t cOutBuffSize = 0x6fffffff;				// buffer size for file writes


// following allele scoring thresholds are extracted from KAligner.h, you must ensure that these thresholds are updated each time KAligner is updated
const double cScorePBA3MinProp = 0.75;		// score PBA as 3 if allele proportion of all counts is >= this threshold and coverage is >= 5
const double cScorePBA2MinProp = 0.35;		// score PBA as 2 if allele proportion of all counts is >= this threshold and coverage is >= 5
const double cScorePBA1MinProp = 0.20;		// score PBA as 1 if allele proportion of all counts is >= this threshold and coverage is >= 5
// when coverage is less than 5 (usually the case with skim reads) and thus confidence in alleles is reduced then scores are reduced
const double cScorePBA2MinLCProp = 0.70;		// score PBA as 2 if allele proportion of all counts is >= this threshold and coverage is < 5
const double cScorePBA1MinLCProp = 0.30;		// score PBA as 1 if allele proportion of all counts is >= this threshold and coverage is >= 5

typedef enum TAG_eModeSH {
	eMSHDefault = 0,		// default is to generate segmentations using bin counts of unique loci only
	eMSHAlleles,			// report allele sites present in a genome PBA
	eMSHPlaceHolder			// used to mark end of processing modes
	}eModeSH;

typedef enum TAG_eHapClass	// confidence in a haplotype for given founder being present or absent within a given bin
{
	eHCAbsent = 0,		// high confidence that founder has no haplotype
	eHCLoConf,			// low confidence that founder has a haplotype
	eHCHiConf			// high confidence that founder has a haplotype
} eHapClass;

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
uint8_t ReadsetType;			// 0: founder, 1: skim, 2: control
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
uint32_t StartLoci;				// window starting at this loci
uint32_t EndLoci;				// window ending at this loci inclusive
uint32_t MultiFounder;			// number of loci in window with multiple founders having same highest scores
uint32_t NoFounder;				// number of loci in window having no skim PBA
uint8_t SmthdFndrHaps[cMaxFounderReadsets];	// smoothed -classification of individual founder haplotype presence in this bin - eHCHiConfNone,eHCHiConfPresent 
uint8_t RawFndrHaps[cMaxFounderReadsets];	// raw - sans smoothing -classification of individual founder haplotype presence in this bin - eHCAbsent, eHCLoConf, eHCHiConf
uint32_t NonUniques[cMaxFounderReadsets];	// Kmers shared - non-unique -with these other founders summed in window over StartLoci through to EndLoci
uint32_t Uniques[cMaxFounderReadsets];		// Kmers uniquely aligning to each founder summed in window over StartLoci through to EndLoci
size_t SummedScores[cMaxFounderReadsets];  // summed scores for each founder in window
} tsBinCnts;

typedef struct TAG_sFndrTotUniqueCnts
{
	uint32_t FounderID;		// this founder
	size_t TotUniqueCnts; // has accumulated this total number of unique counts
} tsFndrTotUniqueCnts;

typedef struct TAG_sFndrSmthdMean
{
	uint32_t FounderID;		// this founder
	double SmthdMean;		// has this raw haplotypes smoothed mean
} tsFndrSmthdMean;

typedef struct TAG_sFndrsClassifications	// sized to hold all founder classifications
{
uint8_t Fndrs[cMaxFounderReadsets];			// holds each founders classification 
} tsFndrsClassifications;

typedef struct TAG_sWorkQueueEl 
{
uint32_t ReadsetID;				// readset being aligned to founders - could be a skim or founder
uint32_t ChromID;				// processing is for this chromosome
uint32_t StartLoci;				// starting from this loci inclusive
uint32_t EndLoci;				// and ending at this loci inclusive
uint32_t LociIncr;				// incrementing loci by this many bases - must be 1, 2, 5, 10, 20 - of which AccumBinSize is a multiple. Used when sampling
uint32_t ChromLen;				// chromosome length
uint32_t AccumBinSize;			// counts are accumulated into this sized (bp) non-overlapping bins
uint32_t NumFndrs;				// number of founders to be processed against the skim (or founder) PBA at *pPBAs[0]
uint8_t *pPBAs[cMaxFounderReadsets]; // pPBAs[0] pts to chromosome skim PBA, followed by ptrs to each of the chromosome founder PBAs
uint8_t *pMskPBA;				// pts to optional chromosome masking PBA, scoring only for segments contained in mask which are non-zero and where skim and founder PBAs also have an allele.
								// enables a founder to be processed as if a skim, restricting scoring Kmers to be same as if a skim

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
	int32_t m_MinAlignedScoringPBA;		// to be accepted as aligned and scored then PBAs in both founder and skim must be at least this PBA level
	int32_t m_MinAlignedPBA;		// to be accepted as aligned ( < m_MinScoringPBA) with no penalty applied then PBAs in both founder and skim must be at least this PBA level
	int32_t m_AlignedPBAScore;		// if PBA accepted for scoring (>= m_MinAlignedScoringPBA) then add this score to the summed Kmer alignment score
	int32_t m_PenaltyPBAScore;		// if PBA not accepted as aligned ( < m_MinAlignedPBA) then apply this penalty to the summed Kmer alignment score
	int32_t m_AcceptKmerMatchPerc;	// only accept Kmer alignment if at least this percentage of base alleles matched at >= m_MinAlignedPBA
	int32_t m_MinKmerScoreDiff;		// must be at least this Kmer percentage score differential between highest scoring founders before accepting a founder as being a uniquely highest scoring

	uint64_t m_FndrsProcMap;			// bitmap (founder1 in bit position 0) corresponding to founders which are to be processed, if bit is not set then the corresponding founder is not to be processed
	uint32_t m_NumFndrsProcMap;			// number of founders in m_FndrsProcMap
	uint32_t m_NumFndrs;				// number of founders
	uint32_t m_Ploidy;					// targeted species is this ploidy - monoploid, diploid, triploid ... 
	uint32_t m_MaskReadsetID;			// masking readset identifier, 0 if no masking readset
	int32_t m_SkimReadsetID;			// skim readset identifier

	int32_t m_LAReadsetNameID;			// name identifier last returned by AddReadsetName()
	int32_t m_NumReadsetNames;			// number of readsets names currently in m_szReadsets
	uint32_t m_NxtszReadsetIdx;		// current concatenated (names separated by '\0') of all readset names in m_szReadsets
	char m_szReadsetNames[(cMaxFounderReadsets + 1) * cMaxDatasetSpeciesChrom];	// used to hold concatenated readset names, each separated by '\0', allowing for 1 skim readset
	uint32_t m_szReadsetIdx[cMaxFounderReadsets + 1];	// array of indexes into m_szFounders giving the starts of each founder name. , allowing for 1 skim readset
	tsReadsetMetadata m_Readsets[cMaxFounderReadsets + 1];	// array of all readset metadata. , allowing for 1 skim readset

	int m_NChroms;									// speeds up debugging - limit processing: 0 All chromosomes, 1..n 1st to nTh chromosome
	int32_t m_LAChromNameID;						// last accessed chromosome identifier from call to AddChrom()
	int32_t m_NumChromNames;						// number of chromosome names currently in m_szChromNames
	int32_t m_NxtszChromIdx;						// current concatenated (names separated by '\0') of all chromosome names in m_szChromNames
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

	uint32_t m_WinMeanBins;		// derive haplotype calls from mean of counts in a sliding window containing this many bins
	uint32_t m_MaxFndrHaps;		// process for reduction down to this maximum number of called haplotypes
	uint8_t m_Fndrs2Proc[cMaxFounderReadsets];	// array of founders which are to be processed, indexed by FounderID-1. If LSB is set then that founder is marked for processing
	uint32_t m_AccumBinSize;	// counts are accumulated over this sized (bp) bins
	int32_t m_MinScoreKmerSize;	// minimum Kmer PBA alignment length for Kmer scoring
	int32_t m_MaxScoreKmerSize;	// m_MinScoreKmerSize can be extended out to this maximal size

	uint32_t m_InNumProcessed;	// this number of input bytes have been processed
	uint32_t m_InNumBuffered;	// currently buffering this many input bytes
	uint32_t m_AllocInBuff;		// m_pInBuffer allocated to hold this many input bytes
	uint8_t *m_pInBuffer;		// allocated for buffering input

	uint32_t m_OutBuffIdx;		// currently buffering this many output bytes
	uint32_t m_AllocOutBuff;	// m_pszOutBuffer allocated to hold this many output bytes
	uint8_t *m_pszOutBuffer;	// allocated for buffering output

	uint32_t m_WIGChromID;				// current WIG span is on this chromosome
	uint32_t m_WIGRptdChromID;			// previously reported WIG chrom
	uint32_t m_WIGSpanLoci;				// current WIG span starts at this loci
	uint32_t m_WIGSpanLen;				// current span is this length
	uint32_t m_WIGRptdSpanLen;			// previously reported WIG span length
	uint64_t m_WIGSpanCnts;				// current span has accumulated this many counts

	int m_hInFile;				// input file handle
	int m_hOutFile;				// file handle for writing results

	char *m_pszMaskBPAFile;		// optional input masking BPA file, only process BPAs which are intersect of these BPAs and skim plus founder BPAs
	

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


	int32_t					// returned readset identifier (1..n) or < 1 if errors
		LoadPBAFile(char* pszFile,		// load chromosome metadata and PBA data from this file
				uint8_t ReadsetType);	// 0: founder, 1: skim, 2: control

	int PBAReportAlleles(char* pszPBAFile,		// load and report allele details present in this PBA file
						char *pszRsltsFileBaseName);	// results are written to this file base name with skim readset identifier and type appended

	int32_t		// returned chrom identifier, < 1 if unable to accept this chromosome name
		AddChrom(char* pszChrom); // associate unique identifier with this chromosome name

	char*	// returns chromosome name associated with 
		LocateChrom(int32_t ChromID); // chromosome identifier

	int32_t		// returned chrom identifier, < 1 if unable to locate this chromosome name
		LocateChrom(char* pszChrom); // return unique identifier associated with this chromosome name

	int32_t		// returned readset name identifier, < 1 if unable to accept this readset name
		AddReadset(char* pszReadset);		// associate unique identifier with this readset name, readset names must be unique including the skim readset!

	char*							// returned readset name
		LocateReadset(int32_t ReadsetID);// identifier returned by call to AddReadsetName

	int32_t		// returned readset name identifier, 0 if unable to locate existing readset
		LocateReadset(char* pszReadset); // associate unique identifier with this readset name

	uint8_t *								// returned pointer to start of PBA
		LocatePBAfor(int32_t ReadSetID,	// readset identifier 
			 int32_t ChromID);				// chrom identifier

	uint32_t								// number buffered
		FillInBuffer(uint32_t MinRequired); // try and fill input buffer with at least this many bytes reading from currently opened input file handle

	uint8_t *AllocPBAs(uint32_t ChromLen);	// allocate to hold at least this many packed base alleles

	int AllocChromMetadata(void);			// allocate for additional chromosome metadata

	tsBinCnts *							// returned ptr to allocated and initialised bin counts
		AddWinBinCnts(tsBinCnts *pInitBinCnts);	// allocated tsWinBinCnts to be initialised with a copy of pInitWinBinCnts

	int ProcessSkimPBAFile(char* pszSkimPBAFile,	// load, process and call haplotypes for this skim PBA file against previously loaded panel founder PBAs
							char *pszRsltsFileBaseName);		// results are written to this file base name with skim readset identifier and type appended

	int
		AlignPBAs(uint32_t ReadsetID,		// this is the readset being aligned to founders - could be a skim or founder
							uint32_t NumFndrs,			// number of founders to be processed against
							int MinKmerSize,			// aligning founders using this minimum Kmer size
							int MaxKmerSize,			// aligning founders using this maximum Kmer size
							uint32_t LociIncr,			// kmer loci starts are incremented by this submultiple of AccumBinSize, used when sampling. Typically 1, 2, 5, 10, 20, 25
							uint32_t AccumBinSize,		// comparing ReadsetID's PBA against founder PBAs then counting using this bin size
							uint32_t NumFndrsProc);		// this many founders in m_Fndrs2Proc

	int CountHaplotypes(uint32_t ReadsetID,		// this is the readset being aligned to founders - could be a skim or founder
							uint32_t NumFndrs,		// number of founders to be processed against
							int MinKmerSize,			// aligning founders using this minimum Kmer size
							int MaxKmerSize,			// aligning founders using this maximum Kmer size
							uint32_t AccumBinSize = 100000);	// comparing Readset's PBA against founder PBAs then counting most probable haplotypes using non-overlapping bins of this size

	int
		ClassifyChromHaplotypes(uint32_t ReadsetID,				// this readset is used for identifying the readset being aligned to founders - could be a skim or founder
								uint32_t ChromID,				// processing is for this chromosome
								uint32_t StartLoci,				// starting from this loci inclusive
								uint32_t EndLoci,				// and ending at this loci inclusive
								uint32_t LociIncr,				// incrementing loci by this many bases - must be 1, 2, 5, 10, 20 - of which AccumBinSize is a multiple. Used when sampling
								uint32_t ChromLen,				// chromosome length
								uint32_t AccumBinSize,			// accumulating counts into this sized (bp) non-overlapping bins
								uint32_t NumFndrs,				// number of founders to be processed against the skim (or founder) PBA at *pPBAs[0]
								uint8_t *pPBAs[],				// pPBAs[0] pts to chromosome skim PBA, followed by ptrs to each of the chromosome founder PBAs
								uint8_t *pMskPBA = NULL);		// pts to optional (NULL if no mask) chromosome masking PBA, scoring only for segments contained in mask which are non-zero and where skim and founder PBAs also have an allele.
																// enables a founder to be processed as if a skim, restricting scoring Kmers to be same as if a skim


	tsBinCnts *			// returned WinBinCnts if located in pReadset
			LocateWinBinCnts(tsBinCnts *pWinBinCnts,	// chrom and start loci in this bin are to be located
								tsReadsetMetadata *pReadset);	// in this readset

	int				// < 0 if errors, otherwise success
		ScorePBAKmer(uint32_t ReadsetID,		// this readset is used for identifying the readset being aligned to founders - could be a skim or founder
					tsBinCnts *pWinBinCnts,		// bin currently being processed, loci is contained within this bin

					uint32_t Loci,				// processing for Kmer haplotypes starting from this loci
					uint32_t SeqLen,		// sequence or chromosome is this length
					uint32_t MinKmerSize,		// Kmer to score on is of this minimum size
					uint32_t MaxKmerSize,		// Kmer to score on can be extended until this maximum size
					uint32_t NumFndrs,		// number of founders to be processed against the skim (or founder) PBA at *pPBAs[0]
					uint8_t *pPBA[],		// pts to loci 0 of PBAs in readset order
					tsFndrsClassifications *pClassifications,		// returned classifications for each processed founder
					uint8_t *pMskPBA = NULL);		// pts to optional masking PBA, scoring only for segments contained in mask which are non-zero and where skim and founder PBAs also have an allele.
												// enables a founder to be processed as if a skim, restricting scoring Kmers to be same as if a skim


	int ChooseHaplotypes(void); // iterate over all the window bins and choose which haplotype to call for that bin interpolating those bins which were initially called as being indeterminate

	int	ClassifyBinHaps(int MaxHaps = 2);	// classify each window bin as containing a combination of this number of maximum haplotypes in that bin

	int SmoothBinsHaps(int Ploidy = 2);		// iterate over all window bins, smoothing using adjacent bins to resolve low confidence bins whilst further reducing number of haplotypes in a bin down to this ploidy 


	double ChiSqr(int NumRows, int32_t* pExp,int32_t *pObs); // currently not used!

	int ReportHaplotypesBED(char *pszOutFile,	// report haplotypes for specified founder to BED pszOutFile
							char *pszFounder,	// this was the founder haplotype
							char *pszSkim,		// against which this skim readset was aligned
							uint32_t FounderID);// reporting on haplotypes for this founder

	void InitialiseWIGSpan(void);				// initialise WIG span vars to values corresponding to no spans having been previously reported
	int CompleteWIGSpan(bool bWrite = false);	// close off any current WIG span ready to start any subsequent span, if bWrite is true then write to disk 
	int AccumWIGCnts(uint32_t ChromID,	// accumulate wiggle counts into variableStep spans over this chromosome
			uint32_t Loci,		// this loci - WIG uses 1 based loci instead of the usual UCSC BED 0 based loci
			uint32_t Cnts,		// has this many counts attributed
			uint32_t MaxSpanLen = 1000000); // allow WIG spans to be this maximal length
	int AccumWIGBinCnts(uint32_t ChromID,	// accumulate wiggle counts into variableStep spans over this chromosome
			uint32_t Loci,		// bin counts starting from this loci - WIG uses 1 based loci instead of the usual UCSC BED 0 based loci
			uint32_t Cnts,		 // bin has this many counts attributed
			uint32_t BinLen);	// bin is this length

	int
		ReportCountsWIG(char *pszOutFile,	// WIG file
							char *pszFounder,	// this was the founder haplotype
							char *pszSkim,		// against which this skim readset was aligned
							uint32_t FounderID); // reporting on exclusive counts for this founder

	CMTqsort m_mtqsort;				// multi-threaded qsort
	static int SortWinBinCnts(const void* arg1, const void* arg2);
	static int SortFndrTotUniqueCnts(const void* arg1, const void* arg2);
	static int SortFndrSmthdMeans(const void* arg1, const void* arg2);

public:
	CCallHaplotypes();
	~CCallHaplotypes();
	void Reset(void);	// resets class instance state back to that immediately following instantiation
	int Process(eModeSH PMode,		// processing mode
			int NChroms,			// debugging - limit processing: 0 All chromosomes, 1..n 1st to nTh chromosome
			char *pszTrackName,		// track name
			char *pszTrackDescr,	// track descriptor
			 int MinKmerSize,				// aligning founders using this minimum Kmer size
			 int MaxKmerSize,				// aligning founders using this maximum Kmer size
			int32_t MinAlignedScoringPBA,		// to be accepted as aligned and scored then PBAs in both founder and skim must be at least this PBA level
			int32_t MinAlignedPBA,		// to be accepted as aligned ( < m_MinScoringPBA) with no penalty applied then PBAs in both founder and skim must be at least this PBA level
			int32_t AlignedPBAScore,		// if PBA accepted for scoring (>= m_MinAlignedScoringPBA) then add this score to the summed Kmer alignment score
			int32_t PenaltyPBAScore,		// if PBA not accepted as aligned ( < m_MinAlignedPBA) then apply this penalty to the summed Kmer alignment score
			int32_t AcceptKmerMatchPerc,	// only accept Kmer alignment if at least this percentage of base alleles matched at >= m_MinAlignedPBA
			int32_t MinKmerScoreDiff,		// must be at least this Kmer percentage score differential between highest scoring founders before accepting a founder as being a uniquely highest scoring
			uint32_t AccumBinSize,	 // accumulating counts into this sized (bp) non-overlapping bins
			int WinMeanBins,		// derive haplotype calls from mean of counts in a sliding window containing this many bins
			int MaxFndrHaps,		// process for reduction down to this maximum number of called haplotypes
			int Ploidy,				// smoothed haplotypes for organism having this ploidy
			char *pszMaskBPAFile,	// optional input masking BPA file, only process BPAs which are intersect of these BPAs and skim plus founder BPAs
			int NumFounderInputFiles,	// number of input founder file specs
			char *pszFounderInputFiles[],	// names of input founder PBA files (wildcards allowed)
			int NumSkimInputFiles,	// number of input skim file specs
			char* pszSkimInputFiles[],		// names of input skim PBA files (wildcards allowed)
			char* pszOutFile,		// output to this file
			int NumThreads);		// number of worker threads to use

	int ProcWorkerThread(tsWorkerInstance *pThreadPar);	// worker thread parameters
};


