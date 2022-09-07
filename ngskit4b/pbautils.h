#pragma once

const int32_t cMaxRefAssembName = 60;				// user specified reference assembly names can be at most this length excluding trailing '\0'
const int32_t cMaxWildCardFileSpecs = 200;			// can accept at most this max number of input wildcarded file specs
const int32_t cMaxPBAReadsets = 4000;				// currently allowing for at most this number of readsets
const int32_t cMaxIncludeChroms = 20;				// max number of include chromosomes regular expressions
const int32_t cMaxExcludeChroms = 20;				// max number of exclude chromosomes regular expressions
const int32_t cMaxChromNames = 100000;				// max number of unique chromsome names allowed
const int32_t cMaxChromMetadata = 1000000;			// allowing for a maximum of this many chromosomes/contigs over all founders
const int32_t cAllocChromMetadata = 100000;			// allocate chrom metadata in this sized increments
const size_t cAllocPackedBaseAlleles = 0x3fffffff;	// allocate packed base alleles to this maximal size, 1 allocation per chromosome per readset
const uint32_t cInBuffSize = 0x7fffffff;			// buffer size for file reads
const uint32_t cOutBuffSize = 0x6fffffff;			// buffer size for file writes
const int cMaxAXAllocBuffChunk = 0x00ffffff;		// buffer for fasta sequences is realloc'd in this sized chunks 
const uint32_t cMaxPBAutilityThreads = 16;			// relatively few threads likely to be required

const char cDfltRefAssemb[] = { "Wm82.v2" };			// default reference assembly if none specified
const char cDfltSeqID[] = { "S000000" };				// default sequence identifier if none specified
const char cDfltExprID[] = { "E000000" };				// default experiment identifier if none specified


// currently only the default processing mode
typedef enum TAG_ePBAuMode {
	ePBAu2Fasta = 0,	// generate Fasta assembly from PBA
	ePBAu2PBA,			// generate PBA from Fasta
	ePBAu2PBAConcordance,	// report on degree of allelic concordance - proportions of loci where all samples had alignments and all samples with exactly same alleles
	ePBAu2WIGConcordance,	// report on degree of coverage depth concordance - proportions of coverage depth (WIG) where all samples had alignments and all samples with near same coverage
	ePBAu2VCF,// output a VCF containing allelic differences between a PBA and the original reference assembly from which the PBA was generated
	ePBAuPlaceholder		// used as a placeholder to mark number of processing modes
} ePBAuMode;

#pragma pack(1)

typedef struct TAG_sChromMetadata
{
	uint32_t ReadsetID;			// chromosome is from this readset
	uint32_t ChromMetadataIdx;	// index of this metadata
	uint32_t NxtChromMetadataIdx; // chromosome metadata in a given readset are forward linked
	uint32_t ChromID;			// chromosome identifier
	uint32_t ChromLen;			// has this many loci bases
	uint32_t HGBinID;            // initial haplotype grouping bin identifier for this chromosome
	int64_t FileOfsPBA;         // PBAs for this chromosome starts at this file offset
	uint8_t* pPBAs;				// pts to memory allocation holding packed base alleles for this chromosome
} tsChromMetadata;

typedef struct TAG_sReadsetMetadata
{
	uint32_t ReadsetID;				// identifies from which readset these packed base alleles were generated
	char szExperimentID[100];		// sequencing experiment
	char szRefAssemblyID[100];		// alignments were against this target assembly
	char szFileName[_MAX_PATH];     // readset loaded from this path+file
	uint8_t ReadsetType;			// 0: primary or founder, 1: progeny, 2: control
	uint32_t NumChroms;				// readset has this number of chromosomes
	uint32_t StartChromMetadataIdx;	// index of starting chrom metadata for this readset
	uint32_t StartChromID;			// starting chrom for this readset
	int64_t NxtFileChromOfs;        // file offset at which the next file chromosome metadata can be read from 
} tsReadsetMetadata;


#pragma pack()

class CPBAutils
{
	char m_szRefAssembFile[_MAX_PATH];	// PBA file containing reference assembly sequences
	char m_szRefAssemb[cMaxRefAssembName+1];    // assembly name to use when converting Fasta to PBA
	char m_szOutFile[_MAX_PATH];				// write converted format into this output file 

	int m_hInFile;				// input file handle
	int64_t m_InFileOfs;        // input file next read will start from this file offset
	int m_hOutFile;				// file handle for writing results

	int32_t m_NumReadsetIDs;			// total number of input readsets loaded for processing
	uint32_t m_LAReadsetNameID;			// name identifier last returned by AddReadsetName()
	uint32_t m_NumReadsetNames;			// number of readsets names currently in m_szReadsets
	uint32_t m_NxtszReadsetIdx;			// current concatenated (names separated by '\0') of all readset names in m_szReadsets, note that the readset type (founder 0, progeny 1 or control 2) is prepended to each name
	char m_szReadsetNames[cMaxPBAReadsets * (cMaxDatasetSpeciesChrom + 1)];	// used to hold concatenated readset names, each separated by '\0'
	uint32_t m_szReadsetIdx[cMaxPBAReadsets];	// array of indexes into m_szReadsetNames giving the starts of each name
	tsReadsetMetadata m_Readsets[cMaxPBAReadsets];	// array of all readset metadata

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
	tsChromMetadata* m_pChromMetadata;	// allocated to hold metadata for all founder chromosomes

	uint32_t m_InNumProcessed;	// this number of input bytes have been processed
	uint32_t m_InNumBuffered;	// currently buffering this many input bytes
	uint32_t m_AllocInBuff;		// m_pInBuffer allocated to hold this many input bytes
	uint8_t* m_pInBuffer;		// allocated for buffering input

	uint32_t m_OutBuffIdx;		// currently buffering this many output bytes
	uint32_t m_AllocOutBuff;	// m_pszOutBuffer allocated to hold this many output bytes
	uint8_t* m_pOutBuffer;	// allocated for buffering output

	int64_t m_NumPBAbases;		// number of PBA called Fasta consensus bases
	int64_t m_NumDiAllelic;		// number of PBA called Fasta diallelic consensus bases
	int64_t m_NumMonoAllelic;	// number of PBA called Fasta monoallelic consensus bases
	int64_t m_NumNonAllelic;	// number of PBA called Fasta indeterminate consensus bases

	uint32_t m_NumIncludeChroms;				// number of RE include chroms
	uint32_t m_NumExcludeChroms;				// number of RE exclude chroms

	uint32_t m_PBAsTrim5;			// trim this many aligned PBAs from 5' end of aligned segments - reduces false alleles due to sequencing errors
	uint32_t m_PBAsTrim3;			// trim this many aligned PBAs from 3' end of aligned segments - reduces false alleles due to sequencing errors

	uint32_t m_WIGChromID;				// current WIG span is on this chromosome
	uint32_t m_WIGRptdChromID;			// previously reported WIG chrom
	uint32_t m_WIGSpanLoci;				// current WIG span starts at this loci
	uint32_t m_WIGSpanLen;				// current span is this length
	uint32_t m_WIGRptdSpanLen;			// previously reported WIG span length
	uint64_t m_WIGSpanCnts;				// current span has accumulated this many counts

	uint32_t m_WIGBuffIdx;		// currently buffering this many output bytes
	uint32_t m_AllocWIGBuff;	// m_pszOutBuffer allocated to hold this many output bytes
	uint8_t* m_pszWIGBuff;	    // allocated for buffering output
	int m_hWIGOutFile;			// file handle for writing results
	char m_szWIGFile[_MAX_PATH];// write WIG format coverage into this output file 

	CBEDfile* m_pBedFile;	// BED file containing reference assembly chromosome names and sizes
	CUtility m_RegExprs;            // regular expression processing

	bool					// true if chrom is accepted, false if chrom not accepted
		AcceptThisChromID(uint32_t ChromID);

	bool					// true if chrom is accepted, false if chrom not accepted
		AcceptThisChromName(char* pszChrom,   // chromosome name
							bool bKnown = true);	// if true then chromosome must have been previously processed and accepted by LoadChromSizes() processing - size is known!

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


	void Reset(void);		// resets class instance state back to that immediately following instantiation

	uint32_t		// returned chrom identifier, 0 if unable to accept this chromosome name
		AddChrom(char* pszChrom); // associate unique identifier with this chromosome name

	char*	// returns chromosome name associated with 
		LocateChrom(uint32_t ChromID); // chromosome identifier

	uint32_t		// returned chrom identifier, 0 if unable to locate this chromosome name
		LocateChrom(char* pszChrom); // return unique identifier associated with this chromosome name

	uint32_t		// returned readset name identifier, 0 if unable to accept this readset name
		AddReadset(char* pszReadset,		// associate unique identifier with this readset name, readset names must be unique within the readset type
			uint8_t ReadsetType);	// 0: founder, 1: progeny, 2: control

	char*							// returned readset name
		LocateReadset(uint32_t ReadsetID);	// identifier returned by call to AddReadsetName

	uint32_t		// returned readset name identifier, 0 if unable to locate existing readset
		LocateReadset(char* pszReadset, // associate unique identifier with this readset name
			uint8_t ReadsetType);	// 0: founder, 1: progeny, 2: control

	etSeqBase
		Alleles2Base(uint8_t PBA);           // identify and return consensus base from a possibly diallelic PBA

	char
		PBAFastaBase(uint8_t PBA); 	// returns Fasta base char - 'a','c','g','t' or 'n' - as being the consensus allele in PBA

	int
		PBA2Fasta(uint32_t ReadSetID);

	uint8_t*								// returned pointer to start of PBA
		LocatePBAfor(uint32_t ReadSetID,	// readset identifier 
			uint32_t ChromID);				// chrom identifier

	uint8_t LocateReadsetChromLociAlleles(uint32_t ReadsetID,	// return alleles for this readset 
		uint32_t ChromID,		// on this chromosome
		uint32_t Loci);			// at this loci

	uint8_t*   // loaded PBAs or NULL if errors loading PBAs for SampleID.ChromID
		LoadSampleChromPBAs(uint32_t SampleID,   // Sample identifier
			uint32_t ChromID);    // chrom identifier specifying which PBAs is to be loaded from SampleID file

	int
		LoadChromPBAs(uint32_t ChromID,	// load PBAs for this chromosome
			uint32_t StartSampleID = 1,	// chromosome PBAs for this sample through 
			uint32_t EndSampleID = 0);	// to this sample inclusive, 0 if from StartSampleID through to last sample

	bool  // false: no previous memory allocated for containing the chromosome PBAs, true: allocation was present and has been deleted
		DeleteSampleChromPBAs(uint32_t SampleID,   // Sample identifier
			uint32_t ChromID);    // chrom identifier

	tsChromMetadata*								// returned pointer to chromosome metadata
		LocateChromMetadataFor(uint32_t ReadSetID,		// readset identifier 
			uint32_t ChromID);			// chrom identifier

	uint8_t* AllocPBAs(uint32_t ChromLen);	// allocate to hold at least this many packed base alleles

	int AllocChromMetadata(void);			// allocate for additional chromosome metadata

	uint32_t				// returns number of unprocessed bytes in buffer
		FillInBuffer(uint32_t MinRequired, uint32_t MaxRequired); // try and fill input buffer with at least MinRequired or refill (if < MinRequired) up to MaxRequired (if 0 then max buffer allocation)

	int32_t					// returned readset identifier (1..n) or < 0 if errors
		LoadPBAFile(char* pszFile,	// load chromosome metadata and PBA data from this file, filters out chroms - AcceptThisChromName()
			uint8_t ReadsetType,	// 0: founder, 1: progeny, 2: control
			bool bChromMetaOnly);  // load chrom metadata (chrom name,length, file offset at which chrom PBAs start) but don't actually load the chromosome PBAs

	int
		Fasta2PBA(char* pszExperimentID, // experiment identifier
			char* pszReferenceID,	  // reference identifier
			char* pszReadsetID,		  // readset identifier
			char* pszInFile,		  // Fasta format input file
			char* pszOutFile);		  // PBA format output file

	bool ValidatePBA(uint8_t* pAlleles,		// validate that PBA alleles are properly conformant
		bool bSetNoAlleles = true, // if non-conformant then overwrite *pAlleles to be no alleles present
		bool bNormalise = true);    // normalise for coverage, low coverage major alleles (2,0,0,0) are normalised to high coverage (3,0,0,0)


	// trim 5' and 3' aligned segments within the PBAs attempting to reduce sequencing error induced false alleles
	int			// returns number of non-trimmed loci in the pPBAs
		TrimPBAs(uint32_t Trim5,	// trim 5' this many aligned PBA bases from each aligned segment
			uint32_t Trim3,			// trim 3' this many aligned PBA bases from each aligned segment
			uint32_t PBALen,		// pPBAs contains this many packed base alleles
			uint8_t* pPBAs);		// packed base alleles to be processed

	int			// returned number of PBAs which are non-conformant
		ValidatePBAs(uint32_t Length, uint8_t* pPBAs, // validate sequence of PBAs
			bool bSetNoAlleles = true, // if non-conformant then overwrite *pAlleles to be no alleles present
			bool bNormalise = true);    // normalise for coverage, low coverage major alleles (2,0,0,0) are normalised to high coverage (3,0,0,0)

	void InitialiseWIGSpan(void);				// initialise WIG span vars to values corresponding to no spans having been previously reported
	int CompleteWIGSpan(bool bWrite = false);	// close off any current WIG span ready to start any subsequent span, if bWrite is true then write to disk 
	int AccumWIGCnts(uint32_t ChromID,	// accumulate wiggle counts into variableStep spans over this chromosome
		uint32_t Loci,		// span is starting at this loci
		uint32_t BinLen,       // span bin is this length
		uint32_t Cnts,		// has this many counts attributed to the bin
		uint32_t MaxSpanLen = 10000000); // allow reported variableStep WIG spans to be continuous up this maximal length and then start another span

	int  // simulating a in-memory PBA with PBA loci alleles replaced by coverage obtained from a WIG coverage file generated at the same time as the PBA file - enables haplotype coverage grouping
		LoadPBACoverage(char* pszInPBA);   // file containing PBA, file name extension will be replaced with 'coverage.wig' which will be expected to be name of file containing the WIG format coverage

	uint8_t* LoadPBAChromCoverage(uint32_t ReadsetID, // loading chromosome coverage for this readset and mapping coverage as though PBAs
		uint32_t ChromID);    // coverage is on this chromosome

	int32_t GeneratePBAConcordance(void);	// generate allelic concordance across all sample PBAs

	int32_t GenerateWIGConcordance(void);	// generate coverage depth concordance across all sample WIGs

	int32_t	LoadChromSizes(char* pszBEDFile); // BED file containing chromosome names and sizes

	int32_t GenAllelicDiffs(int32_t RefPBAID,	// reference PBAs identifier
						int32_t AllelicPBAID);	// PBAs with allelic differences

public:
	CPBAutils();	// constructor
	~CPBAutils();	// destructor

	int Process(ePBAuMode PMode,	// ePBAu2Fasta: PBA to Fasta, ePBAu2PBA: Fasta to PBA, ePBAu2PBAConcordance: concordance over all PBA samples, ePBAu2WIGConcordance: concordance over all WIG samples
		int32_t LimitPBAs,			// limit number of loaded PBA files to this many. 1...cMaxPBAReadsets
		int32_t PBAsTrim5,			// trim this many aligned PBAs from 5' end of aligned segments - reduces false alleles due to sequencing errors
		int32_t PBAsTrim3,			// trim this many aligned PBAs from 3' end of aligned segments - reduces false alleles due to sequencing errors
		char* pszRefAssemb,			// reference assembly identifier
		char* pszRefAssembFile,		// PBA file containing reference assembly sequences
		char* pszChromFile,			// BED file containing chromosome names and sizes
		char* pszSeqID,				// sequence identifier
		char* pszExprID,			// experiment identifier
		int NumInputFiles,			// number of input founder file specs
		char* pszInputFiles[],		// names of input founder PBA files (wildcards allowed)
		char* pszOutFile,			// output to this file
		int NumIncludeChroms,		// number of chromosome regular expressions to include
		char* pszIncludeChroms[],	// array of include chromosome regular expressions
		int NumExcludeChroms,		// number of chromosome expressions to exclude
		char* pszExcludeChroms[],	// array of exclude chromosome regular expressions
		int NumThreads);			// number of worker threads to use
};

