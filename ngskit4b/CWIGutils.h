#pragma once

const int32_t cMaxWildCardFileSpecs = 200;			// can accept at most this max number of input wildcarded file specs
const int32_t cMaxWIGReadsets = 4000;				// currently allowing for at most this number of WIG readsets
const int32_t cMaxIncludeChroms = 20;				// max number of include chromosomes regular expressions
const int32_t cMaxExcludeChroms = 20;				// max number of exclude chromosomes regular expressions

const int32_t cMaxChromNames = 100000;				// max number of unique chromsome names allowed
const int32_t cMaxChromMetadata = 1000000;			// allowing for a maximum of this many chromosomes/contigs over all founders
const int32_t cAllocChromMetadata = 100000;			// allocate chrom metadata in this sized increments

const uint32_t cMaxWIGutilityThreads = 4;			// relatively few threads likely to be required
const size_t cAllocPackedBaseAlleles = 0x3fffffff;	// allocate packed base alleles to this maximal size, 1 allocation per chromosome per readset

typedef enum TAG_eWIGuMode {	// WIG processing modes
	eWIGu2DEseq = 0,  // default is for processing WIGs into DEseq counts ready for DE analysis
	eWIGu2TODEseq,	  // process WIG into DE transcript turnover pseudo counts ready for DE analysis
	eWIGuPlaceholder  // used as a placeholder to mark number of processing modes
	} eWIGuMode;

#pragma pack(1)

typedef struct TAG_sWUChromMetadata
{
	uint32_t ReadsetID;			// chromosome is from this readset
	uint32_t ChromMetadataIdx;	// index of this metadata
	uint32_t NxtChromMetadataIdx; // chromosome metadata in a given readset are forward linked
	uint32_t ChromID;			// chromosome identifier
	uint32_t ChromLen;			// has this many loci bases
	int64_t FileOfs;			// coverage counts for this chromosome start at this file offset
	uint32_t* pCnts;			// pts to memory allocation holding coverage counts for this chromosome
} tsWUChromMetadata;

typedef struct TAG_sWUReadsetMetadata
{
	uint32_t ReadsetID;				// identifies from which readset these packed base alleles were generated
	char szFileName[_MAX_PATH];     // WIG counts loaded from this path+file
	uint32_t NumChroms;				// readset has this number of chromosomes
	uint32_t StartChromMetadataIdx;	// index of starting chrom metadata for this readset
	uint32_t StartChromID;			// starting chrom for this readset
} tsWUReadsetMetadata;

typedef struct TAG_sChromROICnts {
	uint32_t NumLociCoverage;		// number of exon+intron loci in region with actual coverage
	uint32_t NumExonLociCoverage;	// number of exon only loci in region with actual coverage
	uint64_t TotCnts;				// total counts in region (includes both exon and intronic counts)
	uint64_t Tot50Cnts;				// total counts over 1st 50% of region (includes both exon and intronic counts)
	uint64_t TotExonCnts;			// total counts within exons only
	uint64_t Tot50ExonCnts;			// total counts over 1st 50% of region within exons only
	double ExonTurnOver;			// exon turnover ratio
	double TurnOver;				// exon+intron turn over ratio
	} tsChromROICnts;

typedef struct TAG_sChromROI {
	int32_t ReadsetID;				// ROI applies to this readset
	uint32_t *pCnts;				// loaded counts for chromosome currently being processed
	tsChromROICnts ROICnts;			// current ROI
	} tsChromROI;

#pragma pack()

class CWIGutils {
	eWIGuMode m_PMode;				// processing mode: eWIGu2DEseq WIG to DEseq, eWIGu2TODEseq WIG to turnover DEseq

	int32_t m_NumReadsetIDs;			// total number of input readsets loaded for processing
	uint32_t m_LAReadsetNameID;			// name identifier last returned by AddReadsetName()
	uint32_t m_NumReadsetNames;			// number of readsets names currently in m_szReadsets
	uint32_t m_NxtszReadsetIdx;			// current concatenated (names separated by '\0') of all readset names in m_szReadsets, note that the readset type (founder 0, progeny 1 or control 2) is prepended to each name
	char m_szReadsetNames[cMaxWIGReadsets * (cMaxDatasetSpeciesChrom + 1)];	// used to hold concatenated readset names, each separated by '\0'
	uint32_t m_szReadsetIdx[cMaxWIGReadsets];	// array of indexes into m_szReadsetNames giving the starts of each name
	tsWUReadsetMetadata m_Readsets[cMaxWIGReadsets];	// array of all readset metadata

	int32_t m_LAChromNameID;						// last accessed chromosome identifier from call to AddChrom()
	int32_t m_NumChromNames;						// number of chromosome names currently in m_szChromNames
	int32_t m_NxtszChromIdx;						// current concatenated (names separated by '\0') of all chromosome names in m_szChromNames
	char m_szChromNames[cMaxChromNames * cMaxDatasetSpeciesChrom];	// used to hold concatenated chromosome names, each separated by '\0'
	uint32_t m_szChromIdx[cMaxChromNames];			// array of indexes into m_szChromNames giving the starts of each chromosome name
	uint32_t m_ChromSizes[cMaxChromNames];			// array of chromosome sizes indexed by ChromNameID-1
	uint32_t m_NumChromSizes;						// number of chrom sizes accepted from chrom name+sizes BED file - should be same as m_NumChromNames!!!
	uint32_t m_UsedNumChromMetadata;	// current number of chrom metadata used 
	uint32_t m_AllocdChromMetadata;		// current allocation is for this many 
	size_t m_AllocdChromMetadataMem;	// current mem allocation size for m_pChromMetadata
	tsWUChromMetadata* m_pChromMetadata;	// allocated to hold metadata for all founder chromosomes

	CBEDfile* m_pBedFile;	// BED file containing reference assembly chromosome names and sizes

	char m_szOutFile[_MAX_PATH];	// write converted format into this output file 
	CBEDfile *m_pROIFile;			// BED file containing regions of interest
	CUtility m_RegExprs;            // regular expression processing

	int32_t m_MRA_ROIChromID;	// most recently accessed region of interest chrom identifier (ROI BED chrom identifier) 
	int32_t m_MRA_ChromID;		// most recently accessed chrom identifier (internal chrom identifier)

// loading BED which specifies chrom names and sizes
	int		// returning number of chromosomes parsed from BED file and accepted after filtering for wildcards
		LoadChromSizes(char* pszBEDFile); // BED file containing chromosome names and sizes

			// Note: an assumption is that within the WIG input file ordering is by chromosome ascending
	uint32_t* // loaded WIG counts for requested chromosome or nullptr if unable to load
		LoadChromCoverage(uint32_t ReadsetID, // loading chromosome coverage for this readset
							uint32_t ChromID);    // coverage is for this chromosome

	int32_t				// returned readset identifier (1..n) or < 0 if errors
			InitialiseMetadata(char *pszReadset,	  // readset name
								char* pszInWIGFile);   // initialise chromosome metadata from file containing WIG coverage

	char* LocateReadset(uint32_t ReadsetID);

	uint32_t		// returned Readset identifier, 0 if unable to locate this Readset name
		LocateReadset(char* pszReadset, // return unique identifier associated with this Readset name
					uint8_t ReadsetType);	// 0: founder, 1: progeny, 2: control

	tsWUChromMetadata*								// returned pointer to chromosome metadata
		LocateChromMetadataFor(uint32_t ReadSetID,		// readset identifier 
			uint32_t ChromID);			// chrom identifier

	// NOTE: Readsets are checked for uniqueness as readsets must be unique within a given readetset type
	uint32_t		// returned readset identifier, 0 if unable to accept this readset name
		AddReadset(char* pszReadset, // associate unique identifier with this readset name
			uint8_t ReadsetType);	// 0: founder, 1: progeny, 2: control

	uint32_t		// returned chrom identifier, 0 if unable to accept this chromosome name
		AddChrom(char* pszChrom); // associate unique identifier with this chromosome name

	uint32_t		// returned chrom identifier, 0 if unable to locate this chromosome name
		LocateChrom(char* pszChrom); // return unique identifier associated with this chromosome name

	char* LocateChrom(uint32_t ChromID); // returns chrom name

	bool					// true if chrom is accepted, false if chrom not accepted
		AcceptThisChromName(char* pszChrom,   // chromosome name
									 bool bKnown);	// if true then chromosome must have been previously processed and accepted by LoadChromSizes() processing


	int AllocChromMetadata(void);			// allocate for additional chromosome metadata

	uint32_t*	AllocCnts(uint32_t ChromLen);	// allocate memory to hold at least this many coverage counts

	void	DeleteAllChromCnts(void); // delete all currently loaded PBAs - all sample chroms



	bool  // false: no previous memory allocated for containing the chromosome PBAs, true: allocation was present and has been deleted
		DeleteSampleChromCnts(uint32_t SampleID,   // Sample identifier
					uint32_t ChromID);    // chrom identifier

	bool m_bMutexesCreated;		// will be set true if synchronisation mutexes have been created
	int CreateMutexes(void);
	void DeleteMutexes(void);

#ifdef _WIN32
	HANDLE m_hSerialiseAccess;
	alignas(4) volatile uint32_t m_FastSerialise;	// interlocked access to founder stack allocator -  - AcquireFastSerialise()
	alignas(4)	volatile unsigned int m_NumThreads; // used with synchronous compare and swap (CAS) for serialising access to number of worker threads to use
#else
	pthread_mutex_t m_hSerialiseAccess;
	__attribute__((aligned(4))) volatile uint32_t m_FastSerialise;	// fast serialised access to founder stack allocator - AcquireFastSerialise()
	__attribute__((aligned(4))) volatile unsigned int m_NumThreads; // used with synchronous compare and swap (CAS) for serialising access to number of worker threads to use
#endif
	void AcquireSerialise(void);
	void ReleaseSerialise(void);
	void AcquireFastSerialise(void);
	void ReleaseFastSerialise(void);


public:
	CWIGutils();	// constructor
	~CWIGutils();	// destructor

	void Reset(void);	// reset state back to that immediately following instantiation

	int Process(eWIGuMode PMode,	// processing mode: eWIGu2DEseq WIG to DEseq, eWIGu2TODEseq WIG to turnover DEseq
		int32_t NumInputFiles,		// number of input WIG file specs
		char* pszInputFiles[],		// names of input WIG files (wildcards allowed)
		char* pszOutFile,			// output to this file
		char* pszChromFile,			// BED file containing reference assembly chromosome names and sizes
		char* pszROIFile,			// BED file containing regions/genes of interest
		int NumIncludeChroms,		// number of chromosome regular expressions to include
		char *pszIncludeChroms[],	// array of include chromosome regular expressions
		int NumExcludeChroms,		// number of chromosome expressions to exclude
		char *pszExcludeChroms[],	// array of exclude chromosome regular expressions
		int32_t NumThreads);		// number of worker threads to use
    };

