#pragma once
// This class is targeting the analysis of concordance and/or discrepancies between generated DGTs and QTLs


const int32_t cMaxPBAWorkerThreads = 128;			// limiting number of worker threads to be a max of this many

const int32_t cMaxPBAFileSpecs = 500;				// can accept at most this max number of input wildcarded PBA file specs

const int32_t cMaxPBAFiles = 8000;					// after input wildcard expansion then can accept at most this max number of PBA files

const int cMaxIncludeChroms = 20;					// max number of include chromosomes regular expressions
const int cMaxExcludeChroms = 20;					// max number of exclude chromosomes regular expressions

const int32_t cMaxChromNames = 100000;				// max number of unique chromsome names allowed
const int32_t cMaxChromMetadata = 1000000;			// allowing for a maximum of this many chromosomes/contigs over all founders
const int32_t cAllocChromMetadata = 100000;			// allocate chrom metadata in this sized increments
const size_t cAllocPackedBaseAlleles = 0x3fffffff;	// allocate packed base alleles to this maximal size, 1 allocation per chromosome per readset
const uint32_t cInBuffSize = 0x7fffffff;			// buffer size for file reads
const uint32_t cOutBuffSize = 0x6fffffff;			// buffer size for file writes
const int32_t cDfltFndrTrim5 = 5;					// default is for no trimming - can trim PBAs from 5' end of founder aligned segments - reduces false alleles due to sequencing errors
const int32_t cDfltFndrTrim3 = 5;					// default is for no trimming - can trim PBAs from 3' end of founder aligned segments - reduces false alleles due to sequencing errors
const int32_t cDfltProgTrim5 = 5;					// default is for no trimming - can trim PBAs 5' end of progeny aligned segments - reduces false alleles due to sequencing errors
const int32_t cDfltProgTrim3 = 5;					// default is for no trimming - can trim PBAs 3' end of progeny aligned segments - reduces false alleles due to sequencing errors

const double cDfltMinCoverage = 0.8;				// if coverage < this threshold then class as being low coverage
const double cDfltHomozPropThres = 0.95;			// if proportion of samples in Grp1Prop is >= this proportion then then characterise as homozygous

const int32_t cMaxDGTQTLAllelesHashes = 32771;		// max number of unique tsDGTQTLAlleles hashes, chosen to be highest prime < 2^15
const int cAllocDGTQTLAlleles = 1000000;			// alloc size

const int cMaxWaitThreadsStartup = 120;				// allowing all threads to initialise and report starting up within this many seconds = m_NumWorkerInsts == m_ExpNumWorkerInsts
const int cWorkThreadStackSize = 0x01ffff;			// threads created with stacks this size - no recursive functions so an overkill, but who knows, memory is cheap!

typedef enum TAG_eModeDGTA {
	eDGTADefault = 0,						// processing mode 0: reporting QTL loci instances only
	eDGTADGT,								// processing mode 1: reporting DGT or QTL loci instances
	eDGTAPlaceHolder								// used to mark end of processing modes
}eModeDGTA;


#pragma pack(1)
typedef struct TAG_sCHChromMetadata
{
	int32_t ReadsetID;			// chromosome is from this readset
	int32_t ChromMetadataIdx;	// index of this metadata
	int32_t NxtChromMetadataIdx; // chromosome metadata in a given readset are forward linked
	int32_t ChromID;			// chromosome identifier
	int32_t ChromLen;			// has this many loci bases
	int32_t HGBinID;            // initial haplotype grouping bin identifier for this chromosome
	int64_t FileOfsPBA;         // PBAs for this chromosome starts at this file offset
	uint8_t* pPBAs;				// pts to memory allocation holding packed base alleles for this chromosome
} tsCHChromMetadata;

typedef struct TAG_sCHReadsetMetadata
{
	int32_t ReadsetID;				// identifies from which readset these packed base alleles were generated
	char szExperimentID[100];		// sequencing experiment
	char szRefAssemblyID[100];		// alignments were against this target assembly
	char szFileName[_MAX_PATH];     // readset loaded from this path+file
	uint8_t ReadsetType;			// 0: primary or founder, 1: progeny, 2: control
	int32_t NumChroms;				// readset has this number of chromosomes
	int32_t StartChromMetadataIdx;	// index of starting chrom metadata for this readset
	int32_t StartChromID;			// starting chrom for this readset
	int64_t NxtFileChromOfs;        // file offset at which the next file chromosome metadata can be read from 
} tsCHReadsetMetadata;

typedef struct TAG_sDGTQTLAlleles {
	uint32_t Idx;			// index of this instance in m_pDGTQTLAlleles[]
	uint32_t Nxt;			// m_pDGTQTLAlleles[Nxt-1] of next instance having same hash as this instance, 0 if no more instances with same hash
	int32_t ChromID;		// on this chromosome
	int32_t Loci;			// at this loci
	uint32_t flgDGT:1;		// set if instance contains a DGT 
	uint32_t flgQTL : 1;	// set if instance contains a QTL
	uint32_t flgSamplesLoCov : 1;		// set if samples have low proportional coverage at this loci
	uint32_t  flgQTLRefMismatch:1;		// set if the QTL ref allele does not match the reference assembly allele
	uint32_t flgSamplesRefMismatch:1;	// set if the samples highest frequency diplotype does not exactly match the reference assembly allele
	uint32_t flgSamplesHomoz:1;			// set if samples homozygosity >= default of 0.9
	uint32_t flgSamplesHeteroz : 1;		// set if samples highest frequency allele is < 4x next highest frequency allele
	uint32_t flgSamplesPolyAllelic :1;	// set if flgSamplesHomoz true and stacked samples were polyallelic
	uint8_t flgSamplesMonoAllelic:1;	// set  if flgSamplesHomoz true and stacked samples were monoallelic
	uint8_t RefAlleles;		// reference assembly having these alleles present = note 1 allele if consensus but could be 2 if diplotype
	uint8_t DGTAlleles[4];	// DGTs having these alleles present, indexed by alleles A..T with value being that group identified by the respective alelle
	uint8_t QTLAlleles[2];	// QTLs having these alleles present, indexed by Ref=0 and Alt=1
} tsDGTQTLAlleles;

typedef struct TAG_sCHWorkerLoadChromPBAsInstance {
	int ThreadIdx;					// uniquely identifies this thread
	void* pThis;					// will be initialised to pt to class instance
#ifdef _WIN32
	HANDLE threadHandle;			// handle as returned by _beginthreadex()
	uint32_t threadID;				// identifier as set by _beginthreadex()
#else
	int threadRslt;					// result as returned by pthread_create ()
	pthread_t threadID;				// identifier as set by pthread_create ()
#endif
	int32_t StartSampleID;         // processing to start from this sample identifer
	int32_t EndSampleID;           // ending with this sample identifier inclusive
	int32_t ChromID;               // loading PBAs for this chromosome
	bool bNormAlleles;				// true to normalise alleles such that individual alleles can be compared without regard to the proportional coverage (0x22 -> 0x33 as an example)
	int Rslt;						// processing result
} tsCHWorkerLoadChromPBAsInstance;

typedef struct TAG_sCHWorkerInstance {
	int ThreadIdx;					// uniquely identifies this thread
	void* pThis;					// will be initialised to pt to class instance
#ifdef _WIN32
	HANDLE threadHandle;			// handle as returned by _beginthreadex()
	uint32_t threadID;				// identifier as set by _beginthreadex()
#else
	int threadRslt;					// result as returned by pthread_create ()
	pthread_t threadID;				// identifier as set by pthread_create ()
#endif
	int Rslt;						// processing result
} tsCHWorkerInstance;

#pragma pack()


class CDGTvQTLs
{
	eModeDGTA m_PMode;							// processing mode
	char* m_pszRsltsFileBaseName;				// output results file basename

	uint32_t m_FndrTrim5;						// trim this many aligned PBAs from 5' end of founder aligned segments - reduces false alleles due to sequencing errors
	uint32_t m_FndrTrim3;						// trim this many aligned PBAs from 3' end of founder aligned segments - reduces false alleles due to sequencing errors
	uint32_t m_ProgTrim5;						// trim this many aligned PBAs from 5' end of progeny aligned segments - reduces false alleles due to sequencing errors
	uint32_t m_ProgTrim3;						// trim this many aligned PBAs from 3' end of progeny aligned segments - reduces false alleles due to sequencing errors

	int32_t m_LAChromNameID;					// last accessed chromosome identifier from call to AddChrom()
	int32_t m_NumChromNames;					// number of chromosome names currently in m_szChromNames
	int32_t m_NxtszChromIdx;					// current concatenated (names separated by '\0') of all chromosome names in m_szChromNames
	char m_szChromNames[cMaxChromNames * cMaxDatasetSpeciesChrom];	// used to hold concatenated chromosome names, each separated by '\0'
	int32_t m_szChromIdx[cMaxChromNames];		// array of indexes into m_szChromNames giving the starts of each chromosome name

	int32_t m_ChromSizes[cMaxChromNames];		// array of chromosome sizes indexed by ChromNameID-1
	int32_t m_NumChromSizes;					// number of chrom sizes accepted from chrom name+sizes BED file - should be same as m_NumChromNames!!!
	int64_t m_TotChromSizes;					// sum of all chrom sizes

	int32_t m_UsedNumChromMetadata;				// current number of chrom metadata used 
	int32_t m_AllocdChromMetadata;				// current allocation is for this many 
	size_t m_AllocdChromMetadataMem;			// current mem allocation size for m_pChromMetadata
	tsCHChromMetadata* m_pChromMetadata;		// allocated to hold metadata for all chromosomes

	int32_t m_LAReadsetNameID;					// name identifier last returned by AddReadsetName()
	int32_t m_NumReadsetNames;					// number of readsets names currently in m_szReadsets
	int32_t m_NxtszReadsetIdx;					// current concatenated (names separated by '\0') of all readset names in m_szReadsets, note that the readset type (founder 0, progeny 1 or control 2) is prepended to each name
	char m_szReadsetNames[(cMaxPBAFiles + 1) * (cMaxDatasetSpeciesChrom + 1)];	// used to hold concatenated readset names, each separated by '\0', allowing for 1 reference readset
	int32_t m_NumReadsetTypes[256];				// total numbers of each readset types - currently only 0: founder, 1: progeny, 2: control utilized but in future may incorporate more types
	int32_t m_FndrIDMappings[cMaxPBAFiles + 1];	// holds founder identifiers in order of addition via AddReadset
	int32_t m_szReadsetIdx[cMaxPBAFiles + 1];	// array of indexes into m_szReadsetNames giving the starts of each name. , allowing for a reference readset name
	tsCHReadsetMetadata m_Readsets[cMaxPBAFiles + 1];	// array of all readset metadata


	CBEDfile* m_pBedFile;						// BED file containing reference assembly chromosome names and sizes
	CUtility m_RegExprs;						// regular expression processing

	int32_t m_NumFounders;						// number of founders
	uint32_t m_Fndrs2Proc[cMaxPBAFiles+1];		// array of samples which are to be processed, indexed by SampleID-1. If LSB is set then that founder is marked for processing
	uint32_t m_InNumProcessed;					// this number of input bytes have been processed
	uint32_t m_InNumBuffered;					// currently buffering this many input bytes
	uint32_t m_AllocInBuff;						// m_pInBuffer allocated to hold this many input bytes
	uint8_t* m_pInBuffer;						// allocated for buffering input

	uint32_t m_OutBuffIdx;						// currently buffering this many output bytes
	uint32_t m_AllocOutBuff;					// m_pszOutBuffer allocated to hold this many output bytes
	uint8_t* m_pszOutBuffer;					// allocated for buffering output

	int m_hInFile;								// input file handle
	int64_t m_InFileOfs;						// input file next read will start from this file offset
	int m_hOutFile;								// file handle for writing results

	CCSVFile* m_pInCSVFile;						// holds input DGT or QTL CSV file

	uint32_t m_NumDGTsLoaded;					// number of DGTs parsed and loaded
	uint32_t m_NumQTLsLoaded;					// number of QTLs parsed and loaded

	uint32_t m_UsedDGTQTLAlleles;				// number of instances, combines DGT and QTL where same chrom.loci intersects
	uint32_t m_AllocDGTQTLAlleles;				// instances allocated
	size_t m_AllocDGTQTLAllelesMem;				// memory allocated
	tsDGTQTLAlleles *m_pDGTQTLAlleles;			// will be allocated to hold instances of DGTQTLAlleles, these instances will be linked on same hash from m_DGTQTLAllelesHashes[hash] through tsDGTQTLAlleles.NxtOfs 
	uint32_t m_DGTQTLAllelesHashes[cMaxDGTQTLAllelesHashes];	// instances of tsDGTQTLAlleles are hashed on tsDGTQTLAlleles.ChromID andtsDGTQTLAlleles.Loci, instances with same hashes are linked by tsDGTQTLAlleles.NxtOfs

	uint32_t m_NumQTLInstances;					// number of QTLs characterised
	uint32_t m_NumDGTInstances;					// number of DGTs characterised
	uint32_t m_NumDGTQTLbothInstances;			// number of both a DGT and QTL present, intersect, at same loci
	uint32_t m_NumQTLRefMismatched;				// number of times the QTL claimed reference base mismatches the reference assembly base
	uint32_t m_NumSamplesLoCov;					// number of DGT/QTLs having insufficient coverage to further characterise
	// m_NumSamplesHeteroz is dependent on samples having sufficient coverage
	uint32_t m_NumSamplesHeteroz;				// number of times samples had highest frequency allele less than 4x the next highest
	// m_NumSamplesHomoz is dependent on samples not being heterozygous
	uint32_t m_NumSamplesHomoz;					// number of times samples loci have homozygosity >= 0.9
	// the following m_NumSamplesRefMismatched, m_NumSampleBiAllelic are dependent on samples having sufficient coverage and being homozygous
	uint32_t m_NumSamplesRefMismatched;			// number of times samples with homozygosity had alleles matching the reference assembly allele
	uint32_t m_NumSamplesMonoAllelic;			// number of times samples which were homozygous were monoallelic
	uint32_t m_NumSamplesPolyAllelic;			// number of times samples which were homozygous were polyallelic

	double m_MinCoverage;						// if coverage < this threshold then class as being low coverage
	double m_HomozPropThres;					// group 1 proportion of samples is >= than this proportion then characterise as homozygous

	tsDGTQTLAlleles*							// returned instance or nullptr if unable to locate
		LocateDGTQTLAlleles(int32_t ChromID,	// instance on this chrom
			int32_t Loci);						// and at this loci

	tsDGTQTLAlleles*							// returned instance or nullptr if unable to locate
		LocateDGTQTLAlleles(int32_t ChromID);	// instance on this chrom and at the lowest loci (5')


	int32_t									// eBSFSuccess or error
		AddDGTQTLAlleles(tsDGTQTLAlleles* pInitAlleles,bool bIsQTL = false); // add instance of tsDGTQTLAlleles

	int	LoadChromSizes(char* pszBEDFile);						// BED file containing chromosome names and sizes

	int32_t									// returned chrom identifier, 0 if unable to accept this chromosome name
		AddChrom(char* pszChrom); // associate unique identifier with this chromosome name

	char*	// returns chromosome name associated with 
		LocateChrom(int32_t ChromID);		// chromosome identifier

	int32_t								// returned chrom identifier, 0 if unable to locate this chromosome name
		LocateChrom(char* pszChrom);		// return unique identifier associated with this chromosome name

	tsCHChromMetadata*						// returned pointer to chromosome metadata
		LocateChromMetadataFor(int32_t ReadSetID,		// readset identifier 
			int32_t ChromID);				// chrom identifier

	bool										// true if chrom is accepted for processing, false if chrom not accepted
		AcceptThisChromID(uint32_t ChromID);

	bool										// true if chrom is accepted for processing, false if chrom not accepted
		AcceptThisChromName(char* pszChrom,		// chromosome name
			bool bKnown = true);	// if true then chromosome must have been previously processed and accepted by LoadChromSizes() processing - size is known!

	int32_t		// returned readset name identifier, 0 if unable to accept this readset name
		AddReadset(char* pszReadset,		// associate unique identifier with this readset name, readset names must be unique within the readset type
			uint8_t ReadsetType);	// 0: founder, 1: progeny, 2: control

	char*									// returned readset name
		LocateReadset(int32_t ReadsetID);	// identifier returned by call to AddReadsetName

	int32_t		// returned readset name identifier, 0 if unable to locate existing readset
		LocateReadset(char* pszReadset,		// associate unique identifier with this readset name
			uint8_t ReadsetType);	// 0: founder, 1: progeny, 2: control

	int32_t	// returned readset identifier index for founder type readsets, -1 if unable to locate readset identifier
		LocateReadsetIDIdx(int32_t ReadsetID,		// requiring idx for this readset identifier
			uint8_t ReadsetType);			// 0: founder, 1: progeny, 2: control

	uint8_t*								// returned pointer to start of PBA
		LocatePBAfor(int32_t ReadSetID,	// readset identifier 
			int32_t ChromID);				// chrom identifier

	bool m_bMutexesCreated;						// will be set true if synchronisation mutexes have been created
	int CreateMutexes(void);
	void DeleteMutexes(void);

#ifdef WIN32
	HANDLE m_hSerialiseAccess;
	alignas(4) volatile uint32_t m_FastSerialise;	// interlocked access to founder stack allocator -  - AcquireFastSerialise()
	alignas(4)	volatile unsigned int m_ReqTerminate;	// used with synchronous compare and swap (CAS) for serialising access to thread early terminate request
	alignas(4)	volatile unsigned int m_NumThreads;		// used with synchronous compare and swap (CAS) for serialising access to number of worker threads to use
	alignas(4)	volatile unsigned int m_NumWorkerInsts; // used with synchronous compare and swap (CAS) for serialising access to actual number of worker instances
	alignas(4)	volatile unsigned int m_ExpNumWorkerInsts; // used with synchronous compare and swap (CAS) for serialising access to expected number of worker instances
	alignas(4)	volatile unsigned int m_CompletedWorkerInsts; // used with synchronous compare and swap (CAS) for serialising access to  number of completed worker instances
#else
	pthread_mutex_t m_hSerialiseAccess;
	__attribute__((aligned(4))) volatile uint32_t m_FastSerialise;	// fast serialised access to founder stack allocator - AcquireFastSerialise()
	__attribute__((aligned(4))) volatile unsigned int m_ReqTerminate;	// used with synchronous compare and swap (CAS) for serialising access to thread early terminate request
	__attribute__((aligned(4))) volatile unsigned int m_NumThreads;		// used with synchronous compare and swap (CAS) for serialising access to number of worker threads to use
	__attribute__((aligned(4))) volatile unsigned int m_NumWorkerInsts; // used with synchronous compare and swap (CAS) for serialising access to actual number of worker instances
	__attribute__((aligned(4))) volatile unsigned int m_ExpNumWorkerInsts;	// used with synchronous compare and swap (CAS) for serialising access to expected number of worker instances
	__attribute__((aligned(4))) volatile unsigned int m_CompletedWorkerInsts; // used with synchronous compare and swap (CAS) for serialising access to  number of completed worker instances
#endif

	void AcquireSerialise(void);
	void ReleaseSerialise(void);
	void AcquireFastSerialise(void);
	void ReleaseFastSerialise(void);

	tsCHWorkerInstance m_WorkerInstances[cMaxPBAWorkerThreads];	// to hold all worker instance thread parameters
	tsCHWorkerLoadChromPBAsInstance m_WorkerLoadChromPBAInstances[cMaxPBAWorkerThreads];	// to hold all worker instance thread parameters for loading chromosome PBAs

	int // returns 0 if all threads have started processing, 1 if all threads started and some have completed processing, 2 if not all have started
		StartWorkerLoadChromPBAThreads(int32_t NumThreads,		// there are this many threads in pool
			int32_t StartSampleID,	// processing to start from this sample identifer
			int32_t EndSampleID,	// ending with this sample identifier inclusive
			int32_t ChromID,		// loading PBAs for this chromosome
			bool bNormAlleles = false);	// normalise alleles such that individual alleles can be compared for presence without regard to differences in proportional coverage (0x22 -> 0x33 as an example)

	int				// returns 0 if all threads started and completed processing, 1 if all threads started but some yet to complete within WaitSecs, 2 if not all threads have started within WaitSecs
		WaitWorkerThreadStatus(int WaitSecs);	// // wait at most this many seconds for all threads to start and complete processing before reporting on worker thread status

	// stop all threads in worker pool
	int	// returns number of threads forced to terminate
		TerminateWorkerThreads(int WaitSecs = 120);	// allow at most this many seconds before force terminating threads

	// stop all threads in worker load chromosome PBAs thread pool
	int	// returns number of threads forced to terminate
		TerminateLoadChromPBAsThreads(int WaitSecs = 120); // allow at most this many seconds before force terminating worker pool threads

	int32_t				// error or success (>=0) code 
		ProcessDGTs(char* pszHapGrpDGTFile,             // input, previously generated by 'callhaplotypes', haplotype group DGT file (CSV format) 
			int32_t NumSamplePBAsInputFiles,	     // number of input sample PBAs file specs
			char* pszSamplePBAsInputFiles[],		// names of input sample PBAs files (wildcards allowed)
			char* pszOutFile);						// write results to this output file (CSV format)


	int32_t										// returned readset identifier (1..n) or < 1 if errors
		LoadPBAFile(char* pszFile,				// load chromosome metadata and PBA data from this file
			uint8_t ReadsetType,		// 0: founder, 1: progeny, 2: control
			bool bChromMetaOnly = false);// if true then load chrom metadata (chrom name,length, file offset at which chrom PBAs start) but don't actually load the chromosome PBAs

	uint8_t*   // loaded PBAs or nullptr if errors loading PBAs for SampleID.ChromID
		LoadSampleChromPBAs(int32_t SampleID,	// Sample identifier
			int32_t ChromID);		// chrom identifier specifying which PBAs is to be loaded from SampleID file

	bool  // false: no previous memory allocated for containing the chromosome PBAs, true: allocation was present and has been deleted
		DeleteSampleChromPBAs(int32_t SampleID,// Sample identifier
			int32_t ChromID);	// chrom identifier

	bool ValidatePBA(uint8_t* pAlleles,				// validate that PBA alleles are properly conformant
		bool bSetNoAlleles = true,					// if non-conformant then overwrite *pAlleles to be no alleles present
		bool bNormalise = true);					// normalise for coverage, low coverage major alleles (2,0,0,0) are normalised to high coverage (3,0,0,0)

	int			// returned number of PBAs which are non-conformant
		ValidatePBAs(int32_t Length, uint8_t* pPBAs,// validate sequence of PBAs
			bool bSetNoAlleles = true,				// if non-conformant then overwrite *pAlleles to be no alleles present
			bool bNormalise = true);					// normalise for coverage, low coverage major alleles (2,0,0,0) are normalised to high coverage (3,0,0,0)

	uint32_t								// number buffered
		FillInBuffer(uint32_t MinRequired, uint32_t MaxRequired = 0); // try and fill input buffer with at least MinRequired or refill (if < MinRequired) up to MaxRequired (if 0 then max buffer allocation)


	// trim 5' and 3' aligned segments within the PBAs attempting to reduce sequencing error induced false alleles
	int			// returns number of non-trimmed loci in the pPBAs
		TrimPBAs(uint32_t Trim5,			// trim 5' this many aligned PBA bases from each aligned segment
			uint32_t Trim3,				// trim 3' this many aligned PBA bases from each aligned segment
			uint32_t PBALen,			// pPBAs contains this many packed base alleles
			uint8_t* pPBAs);			// packed base alleles to be processed

	uint8_t
		LocateReadsetChromLociAlleles(int32_t ReadsetID,	// return alleles for this readset 
			int32_t ChromID,		// on this chromosome
			int32_t Loci);			// at this loci

	uint8_t* AllocPBAs(int32_t ChromLen);	// allocate to hold at least this many packed base alleles
	int AllocChromMetadata(void);			// allocate for additional chromosome metadata

	int LoadDGTs(char *pszDGTsFile);		// CSV DGTs file containing chrom.loci and alleles for each group tag
	
	int LoadQTLs(char *pszQTLsFile);		// CSV QTLs file containing chrom.loci and alleles for each QTL

	int
		ProcessDGTQTLsPBAs(eModeDGTA PMode);	// actual processing of DGTs/QTLs against reference and sample PBAs

	int
		AnalyseInstance(eModeDGTA PMode,				// processing mode
			tsDGTQTLAlleles* pDGTQTLAlleles,	// processing this instance
			int32_t NumPBAs,					// against this number of PBAs 
			uint8_t** ppPBAs);					// ptrs to loaded PBAs

	const char*  // returns character representation of a diplotype alleles
		Diplotype2Txt(uint8_t Alleles);

	CMTqsort m_mtqsort;				// multi-threaded qsorts
	static int SortDGTQTLAlleles(const void* arg1, const void* arg2);

public:
	CDGTvQTLs();
	~CDGTvQTLs();
	void Reset(void);					// resets class instance state back to that immediately following instantiation
	
	int Process(eModeDGTA PMode,		// processing mode
		double MinCoverage,				// if no coverage >= this threshold then class as being low coverage
		double HomozPropThres,			// if proportion of samples in Grp1Prop is >= this proportion then characterise as homozygous		char* pszAssembRefFile,			// contains original reference assembly, diplotype PBA, against which samples were aligned
		char* pszAssembRefFile,			// contains original reference assembly, PBA, against which samples were aligned
		char* pszChromFile,				// BED file, contains chromosome names and sizes
		char* pszDGTsFile,				// file containing DGT loci and allele sample groupings
		char* pszQTLsFile,				// file containing QTL SNP loci
		int NumPBAFiles,				// number of input PBA file specs, these are the source PBA files used during the grouping generation
		char** ppszPBAFiles,			// names of input PBA files (wildcards allowed)
		int	NumIncludeChroms,			// number of chromosome regular expressions to include
		char** ppszIncludeChroms,		// array of include chromosome regular expressions
		int	NumExcludeChroms,			// number of chromosome expressions to exclude
		char** ppszExcludeChroms,		// array of exclude chromosome regular expressions
		char* pszRsltsFileBaseName,		// write results to this file base name - will be suffixed with result types
		int NumThreads);				// number of worker threads to use

	int ProcWorkerThread(tsCHWorkerInstance* pThreadPar);	// worker thread parameters
	int ProcWorkerLoadChromPBAThread(tsCHWorkerLoadChromPBAsInstance* pThreadPar);	// worker thread parameters
};

