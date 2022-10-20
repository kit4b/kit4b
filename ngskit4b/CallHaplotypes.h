#pragma once

const int32_t cMaxPBAWorkerThreads = 64;			// limiting number of worker threads to be a max of this many

const int32_t cMaxFounderFileSpecs = 200;			// can accept at most this max number of input wildcarded founder file specs
const int32_t cMaxProgenyFileSpecs = 500;			// can specify a max of this number of input wildcarded progeny file specs

const int32_t cMaxFounderReadsets = 4000;			// after input wildcard expansion then can accept at most this max number of founder readsets
const int32_t cMaxProgenyReadsets = 4000;			// after input wildcard expansion then can specify at most this max number of progeny readsets

const int cMaxIncludeChroms = 20;		// max number of include chromosomes regular expressions
const int cMaxExcludeChroms = 20;		// max number of exclude chromosomes regular expressions

const int32_t cMaxChromNames = 100000;				// max number of unique chromsome names allowed
const int32_t cMaxChromMetadata = 1000000;			// allowing for a maximum of this many chromosomes/contigs over all founders
const int32_t cAllocChromMetadata = 100000;			// allocate chrom metadata in this sized increments
const size_t cAllocPackedBaseAlleles = 0x3fffffff;	// allocate packed base alleles to this maximal size, 1 allocation per chromosome per readset
const uint32_t cInBuffSize = 0x7fffffff;			// buffer size for file reads
const uint32_t cOutBuffSize = 0x6fffffff;			// buffer size for file writes
const int32_t cDfltFndrTrim5 = 0;					// default is for no trimming - can trim PBAs from 5' end of founder aligned segments - reduces false alleles due to sequencing errors
const int32_t cDfltFndrTrim3 = 0;					// default is for no trimming - can trim PBAs from 3' end of founder aligned segments - reduces false alleles due to sequencing errors
const int32_t cDfltProgTrim5 = 0;					// default is for no trimming - can trim PBAs 5' end of progeny aligned segments - reduces false alleles due to sequencing errors
const int32_t cDfltProgTrim3 = 0;					// default is for no trimming - can trim PBAs 3' end of progeny aligned segments - reduces false alleles due to sequencing errors

const int32_t cDfltWWRLProxWindow = 1000000;			// default for proximal window size for Wald-Wolfowitz runs test
const int32_t cDlftOutliersProxWindow = 1000000;		// default for proximal window size for outliers reduction

// PBA at a single loci is packed into a single 8bit byte: Allele A in bits 7.6, C in bits 5.4, G in bits 3.2, T in bits 1.0
// following allele scoring thresholds are extracted from KAligner.h, you must ensure that these thresholds are updated each time KAligner is updated
const double cScorePBA3MinProp = 0.75;		// score PBA as 3 if allele proportion of all counts is >= this threshold and coverage is >= 5
const double cScorePBA2MinProp = 0.35;		// score PBA as 2 if allele proportion of all counts is >= this threshold and coverage is >= 5
const double cScorePBA1MinProp = 0.20;		// score PBA as 1 if allele proportion of all counts is >= this threshold and coverage is >= 5
// when coverage is less than 5 (usually the case with progeny reads) and thus confidence in alleles is reduced then scores are reduced
const double cScorePBA2MinLCProp = 0.70;		// score PBA as 2 if allele proportion of all counts is >= this threshold and coverage is < 5
const double cScorePBA1MinLCProp = 0.30;		// score PBA as 1 if allele proportion of all counts is >= this threshold and coverage is < 5

const uint32_t cAllocProgenyFndrAligns = 10000000; // initial allocation for progeny to founder allele stacks alignments

const uint32_t cAllocAlleleStacks = 1000000;		// initial allocation for allele stacks - number of stacks
const uint32_t cReallocAlleleStacks = (cAllocAlleleStacks/2);	// realloc allele stacks with this sized increments - number of stacks

const uint32_t cAllocHGBinSpecs = 500000;  // initial allocation is to hold this many haplotype group clustering specifications
const uint32_t cReallocHGBinSpecs = (cAllocHGBinSpecs/2);	// realloc haplotype group clustering in this sized increments - number of bins


const int32_t cDfltGrpHapBinSize = 10000;  // default haplotype clustering over these sized KMer non-overlapping windows
const int32_t cDfltMinCentClustDist = 10;  // default haplotype clustering minimum group centroid clustering differential
const int32_t cDfltMaxCentClustDist = 10000;  // default haplotype clustering maximum group centroid clustering differential
const int32_t cDfltMaxClustGrps = 5;    // default number of haplotype groups
const int32_t cMaxClustGrps = 20;        // can specify up to this maximum of haplotype groups, pointless if more!
const int32_t cDfltNumHapGrpPhases = 10; // default is a max of 10 phases in which to converge group consensus when haplotype clustering into groups - a balance between optimal group consensus and processing resources.  

const int32_t cDfltMaxReportGrpQGLs = 10000000; // default is to report this many highest scoring F-measure group QGLs 
const int32_t cDfltMinQGLGrpMembers = 10;      // haplotype groups with fewer than this number of members are considered as if containing noise and alleles in these groups are not used when determining QGL group specific major alleles    
const double cDfltMinQGLGrpPropTotSamples = 0.10; // haplotype groups containing less than this proportion of total samples are treated as if containing noise and alleles in these groups are not used when determining QGL group specific major alleles 
const double cDfltMinQGLFmeasure = 0.90;   // only accepting QGL loci with at least this F-measure score
const double cDlftFbetaGrps = 1.0;         // default Fbeta-measure

const int32_t cDfltSparseRepPropGrpMbrs = 0;  // only apply sparse representative imputation if number of haplotype group members at least this, default is for no sparse imputations
const double cDfltSparseRepProp = 0.25;	   // if highest frequency (consensus) allele is 0x00 (no coverage) then if next highest frequency allele is at least this proportion of all members in haplotype group then treat next highest allele as
										   // being as being the consensus. Objective to obtain an imputed allele in regions of sparse haplotype coverage

typedef enum TAG_eModeCSH {
	eMCSHDefault = 0,		// processing mode 0: report progeny imputation haplotype matrix
	eMCSHRaw,               // report both raw and imputation haplotype matrices
	eMCSHGWAS,				// additionally generate GWAS allowing visual comparisons
	eMCSHAllelicHapsGrps,   // generating allelic haplotype groupings
	eMCSHCoverageHapsGrps,  // generating coverage haplotype groupings
	eMCSHQGLHapsGrps,       // post-processing haplotype groupings for differential QGL loci
	eMCSHGrpDist2WIG,       // post-processing haplotype grouping centroid distances into WIG format
	eMCSHPlaceHolder		// used to mark end of processing modes
	}eModeCSH;

typedef enum TAG_eHapClass	// confidence in a haplotype for given founder being present or absent within a given bin
{
	eHCAbsent = 0,		// high confidence that founder has no haplotype
	eHCLoConf,			// low confidence that founder has a haplotype
	eHCHiConf			// high confidence that founder has a haplotype
} eHapClass;

#pragma pack(1)

typedef struct TAG_sGrpCnts
{
	uint32_t NumMembers; // number of members in this group
	uint32_t NumAllele[5]; // number of members having this allele, indexed by 0=A,1=C,2=G,3=T,4=N
} tsGrpCnts;

typedef struct TAG_sQGLLoci { // loci which have been identified as group QGLs
	uint32_t SrcExprID;  // sourced CSV haplotype group experiment identifier
	uint32_t SrcRowID;   // sourced CSV haplotype group row identifier
	uint32_t ChromID;    // QGL loci is on this chromosome
	uint32_t QGLLoci;    // at this loci
	uint32_t NumHapGrps; // QGL has this number of actual haplotype groups
	uint32_t AlleleGrpIDs[4];   // group identifier for each allele accepted as having the highest F-measure, 0 if no group having accepted F-measure
	double AlleleFMeasures[4];  // F-measure for each allele - 0.0 if no group having an accepted F-measure for that allele
	double AlleleFMeasuresDesc[4];  // F-measure sorted descending without regard to allele ordering - used when filtering for N loci having highest F-measures 
	tsGrpCnts GrpCnts[6];    // counts for 1st 5 groups plus a pseudo-group containing sum of counts for groups after the 1st 5 
	} tsQGLLoci;

typedef struct TAG_s4096Bits {
	uint64_t Bits[64];		// 4096 bits encoded into 64 x 64bit words
} ts4096Bits;

typedef struct TAG_sHaplotypeGroup {
	uint32_t Size;           // this instance is this size (bytes)
	uint32_t SrcExprID;      // haplotype group initialised from CSV row having this ExprID
	uint32_t SrcRowID;      // haplotype group initialised from CSV row having this RowID
	uint32_t ChromID;	    // group is for this chromosome
	uint32_t StartLoci;		// group starting from this loci
	uint32_t NumLoci;	    // covering this many loci
	uint32_t MinCentroidDistance; // centroid distance can range from this minimum and up to MaxCentroidDistance when attempting to limit number of haplotype groupings to MaxNumHaplotypeGroups
	uint32_t MaxCentroidDistance; // maximum centroid distance which can be used when attempting to meet MaxNumHaplotypeGroups constraint
	uint32_t MaxHaplotypeGroups; // attempt to constrain number of haplotype groups to be this maximum
	uint32_t CentroidDistance; // actual centroid distance used when limiting number of haplotype groupings to MaxNumHaplotypeGroups
	uint32_t MRAGrpChromID;	     // most recently updated group consensus PBAs was on this chrom - GenConsensusGroupBPA()
	uint32_t MRAGrpLoci;      // most recently updated group consensus PBAs was at this loci - GenConsensusGroupBPA()
	uint8_t MRAGrpConsensus[256];    // most recently returned haplotype group consensus allele - GenConsensusGroupBPA()
	uint32_t NumFndrs;       // groups contain a total of this many founders
	uint32_t NumHaplotypeGroups;      // contains this many haplotype groups
	ts4096Bits HaplotypeGroup[1]; // could be (likely!) multiple groups
	} tsHaplotypeGroup;

typedef struct TAG_sHGBinSpec { // haplotype group clustering is for this bin specification
	uint32_t ProcState;        // bin processing state  - 0x00 if unprocessed, 0x01 if currently being processing, 0x03 if processing completed but referenced chromosome not located, 0x07 if processing successfully completed       
	uint32_t ChromID;	    // haplotype group clustering is on this chromosome
	uint32_t BinID;          // uniquely identifies bin - globally unique 
	uint32_t StartLoci;      // bin starts at this loci
	uint32_t NumLoci;	    // covering this many loci
	uint32_t MinCentroidDistance; // centroid distance can range from this minimum and up to MaxCentroidDistance when attempting to limit number of haplotype groupings to MaxNumHaplotypeGroups
	uint32_t MaxCentroidDistance; // maximum centroid distance which can be used when attempting to meet MaxNumHaplotypeGroups constraint
	uint32_t MaxNumHaplotypeGroups; // attempt to constrain number of haplotype groups to be this maximum
	tsHaplotypeGroup* pHaplotypeGroup; // haplotype groupings for this bin
	} tsHGBinSpec;

typedef struct TAG_sAlleleStack {	// stacked informative - all dirac alleles and not all alleles the same in stack -, e.g founder alleles differ at the given chrom and loci
	uint32_t AlleleStackID; // uniquely identifies this allele stack instance
	uint32_t ChromID;		// stack is on this chromosome
	uint32_t Loci;			// and at this loci
	uint32_t NumFndrs;		// number of founders
	uint32_t NumProcFndrs;	// number of founders actually processed - some founders may be marked as not for processing
	ts4096Bits ProcFndrs;	// bitmap of founders actually processed
	uint32_t NumAlleleFndrs[4];	// number of founders mapped in corresponding Alleles[], Alleles[0] allele T through to Alleles[3] allele A 
	ts4096Bits Alleles[4];	// bitmaps (founderA in bit0) for each founder allele contained in this stack - max of 4000 founders per allele 
} tsAlleleStack;

typedef struct TAG_FndrHapGroup {   // founder haplotype group - founders which have been grouped having same sequence of haplotype alleles
	uint32_t FndHapID;  // uniquely identifies this instance
	uint8_t GrpChked : 1; // set when group members have been checked as sharing same sequence {Len} of allele 
	uint32_t ChromID;	// haplotype is on this chrom
	uint32_t Loci;		// starting from this loci
	uint32_t Len;        // is this length
	uint32_t NumFndrs;   // shared by this number of founders
	ts4096Bits Fndrs;    // bitmap of sharing founders
	} tsFndrHapGroup;

typedef struct TAG_sProgenyFndrAligns {
	uint32_t AlleleStackID; // progeny alignments derived from this  allele stack instance
	uint32_t ReadsetID;				// identifies the progeny readset
	uint32_t ChromID;				// stack is on this chromosome
	uint32_t Loci;					// and at this loci
	uint32_t FaHap;					// 0 : none assigned 
	uint32_t FbHap;					// 0: none assigned 
	uint8_t Alleles;				// progeny has these alleles present at the allelestack chrom.loci
	uint32_t NumProgenyFounders;	// number of potential progeny founders in ProgenyParents bitmap
	ts4096Bits ProgenyFounders;		// bitmap of potential progeny founders - progeny has a minor/major allele shared with corresponding founder
} tsProgenyFndrAligns;

typedef struct TAG_sCHChromMetadata
{
uint32_t ReadsetID;			// chromosome is from this readset
uint32_t ChromMetadataIdx;	// index of this metadata
uint32_t NxtChromMetadataIdx; // chromosome metadata in a given readset are forward linked
uint32_t ChromID;			// chromosome identifier
uint32_t ChromLen;			// has this many loci bases
uint32_t HGBinID;            // initial haplotype grouping bin identifier for this chromosome
int64_t FileOfsPBA;         // PBAs for this chromosome starts at this file offset
uint8_t *pPBAs;				// pts to memory allocation holding packed base alleles for this chromosome
} tsCHChromMetadata;

typedef struct TAG_sCHReadsetMetadata
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
} tsCHReadsetMetadata;

typedef struct TAG_sCHWorkQueueEl 
{
uint32_t ChromID;				// processing is for this chromosome
uint32_t StartLoci;				// starting from this loci inclusive
uint32_t ChromLen;				// chromosome length
uint32_t MaxNumLoci;				// max number of loci to be processed - starting from StartLoci - to be processed in this queue element
uint32_t NumFndrs;				// number of founders to be processed against the progeny (or founder) PBA at *pPBAs[0]
uint8_t *pFounderPBAs[cMaxFounderReadsets]; // ptrs to each of the chromosome founder PBAs indexed in FounderID ascending order
uint8_t *pMskPBA;				// pts to optional chromosome masking PBA, scoring only for segments contained in mask which are non-zero and where progeny and founder PBAs also have an allele.
								// enables a founder to be processed as if a progeny, restricting scoring Kmers to be same as if a progeny
} tsCHWorkQueueEl;

typedef struct TAG_sCHWorkerInstance {
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
} tsCHWorkerInstance;

typedef struct TAG_sCHWorkerLoadChromPBAsInstance {
	int ThreadIdx;					// uniquely identifies this thread
	void *pThis;					// will be initialised to pt to class instance
#ifdef _WIN32
	HANDLE threadHandle;			// handle as returned by _beginthreadex()
	uint32_t threadID;				// identifier as set by _beginthreadex()
#else
	int threadRslt;					// result as returned by pthread_create ()
	pthread_t threadID;				// identifier as set by pthread_create ()
#endif
	uint32_t StartSampleID;         // processing to start from this sample identifer
	uint32_t EndSampleID;           // ending with this sample identifier inclusive
	uint32_t ChromID;               // loading PBAs for this chromosome
	int Rslt;						// processing result
} tsCHWorkerLoadChromPBAsInstance;


#pragma pack()

class CCallHaplotypes
{
	eModeCSH m_PMode;				// processing mode 0: report imputation haplotype matrix, 1: report both raw and imputation haplotype matrices, 2: additionally generate GWAS allowing visual comparisons, 3: allelic haplotype grouping,4: coverage haplotype grouping, 5: post-process haplotype groupings for QGLs,  6: post-process to WIG
	uint32_t m_LimitPrimaryPBAs;    // limit number of loaded primary or founder PBA files to this many. 0: no limits, > 0 sets upper limit

	uint32_t m_MaxReportGrpQGLs;    // when calling group QGLs then, if non-zero - report this many highest scoring QGLs
	uint32_t m_ExprID;				// assign this experiment identifier for this PBA analysis
	uint32_t m_SeedRowID;           // generated CSVs will contain monotonically unique row identifiers seeded with this row identifier  
	uint32_t m_CurRowID;            // current row identifier, increment after each row assignment  

	uint32_t m_FndrTrim5;			// trim this many aligned PBAs from 5' end of founder aligned segments - reduces false alleles due to sequencing errors
	uint32_t m_FndrTrim3;			// trim this many aligned PBAs from 3' end of founder aligned segments - reduces false alleles due to sequencing errors
	uint32_t m_ProgTrim5;			// trim this many aligned PBAs from 5' end of progeny aligned segments - reduces false alleles due to sequencing errors
	uint32_t m_ProgTrim3;			// trim this many aligned PBAs from 3' end of progeny aligned segments - reduces false alleles due to sequencing errors
	bool m_bAllFndrsLociAligned;    // true if all founders at a given loci must have a PBA alignment for founders to be processed at that alignment - if false then only one founder need have a PBA alignment
	int32_t m_WWRLProxWindow;		// proximal window size for Wald-Wolfowitz runs test
	int32_t m_OutliersProxWindow;	// proximal window size for outliers reduction

	uint32_t m_NumHapGrpPhases;          // number of grouping phases - using more phases converges the group consensus alleles but at the cost of processing times, optimal balance seems to be at 3 phases but this will be subject to change  
	uint32_t m_GrpHapBinSize;            // when grouping into haplotypes (processing mode 3) then use this non-overlapping bin size for accumulation of differentials between all founder PBAs
	uint32_t m_MinCentClustDist;         // haplotype groupings - in processing mode 3 only - minimum group centroid clustering distance
	uint32_t m_MaxCentClustDist;         // haplotype groupings - in processing mode 3 only - maximum group centroid clustering distance
	uint32_t m_MaxClustGrps;             // haplotype groupings - in processing mode 3 only - targeted maximum number of groups

	uint32_t m_MinQGLGrpMembers;         // groups with fewer than this number of members are treated as noise alleles these groups are not used when determining QGL group specific major alleles
	double m_MinQGLGrpPropTotSamples;  // haplotype groups containing less than this proportion of all samples are treated as if containing noise and alleles in these groups are not used when determining QGL group specific major alleles 
	double m_MinQGLFmeasure;            // only accepting QGL loci with at least this F-measure score
	double m_FbetameasureGrps;         // default Fbeta-measure

	uint32_t m_SparseRepPropGrpMbrs;		// only apply sparse representative imputation if number of haplotype group members at least this 
	double m_SparseRepProp;				// if highest frequency (consensus) allele is 0x00 (no coverage) then if next highest frequency allele is at least this proportion of all members in haplotype group then treat next highest allele as
										// being as being the consensus. Objective to obtain an imputed allele in regions of sparse haplotype coverage

	uint32_t m_WIGChromID;				// current WIG span is on this chromosome
	uint32_t m_WIGRptdChromID;			// previously reported WIG chrom
	uint32_t m_WIGSpanLoci;				// current WIG span starts at this loci
	uint32_t m_WIGSpanLen;				// current span is this length
	uint32_t m_WIGRptdSpanLen;			// previously reported WIG span length
	uint64_t m_WIGSpanCnts;				// current span has accumulated this many counts

	uint32_t m_NumFndrsProcMap;			// number of founders in m_FndrsProcMap
	uint32_t m_NumFounders;				// number of founders
	uint32_t m_NumProgenies;				// number of progenies
	uint32_t m_MaskReadsetID;			// masking readset identifier, 0 if no masking readset
	uint32_t m_CurProgenyReadsetID;		// current progeny readset identifier
	uint32_t m_ProgenyIDs[cMaxProgenyReadsets];	// array of progeny readset identifiers in order of loading
	uint32_t m_LAReadsetNameID;			// name identifier last returned by AddReadsetName()
	uint32_t m_NumReadsetNames;			// number of readsets names currently in m_szReadsets
	uint32_t m_NxtszReadsetIdx;			// current concatenated (names separated by '\0') of all readset names in m_szReadsets, note that the readset type (founder 0, progeny 1 or control 2) is prepended to each name
	char m_szReadsetNames[(cMaxFounderReadsets + cMaxProgenyReadsets + 1) * (cMaxDatasetSpeciesChrom + 1)];	// used to hold concatenated readset names, each separated by '\0', allowing for 1 progeny readset
	uint32_t m_szReadsetIdx[cMaxFounderReadsets + cMaxProgenyReadsets + 1];	// array of indexes into m_szReadsetNames giving the starts of each founder or progeny  name. , allowing for a control readset name
	tsCHReadsetMetadata m_Readsets[cMaxFounderReadsets + cMaxProgenyReadsets + 1];	// array of all readset metadata

	uint32_t m_LAChromNameID;						// last accessed chromosome identifier from call to AddChrom()
	uint32_t m_NumChromNames;						// number of chromosome names currently in m_szChromNames
	uint32_t m_NxtszChromIdx;						// current concatenated (names separated by '\0') of all chromosome names in m_szChromNames
	char m_szChromNames[cMaxChromNames * cMaxDatasetSpeciesChrom];	// used to hold concatenated chromosome names, each separated by '\0'
	uint32_t m_szChromIdx[cMaxChromNames];		// array of indexes into m_szChromNames giving the starts of each chromosome name

	uint32_t m_ChromSizes[cMaxChromNames];	// array of chromosome sizes indexed by ChromNameID-1
	uint32_t m_NumChromSizes;				// number of chrom sizes accepted from chrom name+sizes BED file - should be same as m_NumChromNames!!!

	uint32_t m_UsedNumChromMetadata;	// current number of chrom metadata used 
	uint32_t m_AllocdChromMetadata;		// current allocation is for this many 
	size_t m_AllocdChromMetadataMem;	// current mem allocation size for m_pChromMetadata
	tsCHChromMetadata *m_pChromMetadata;	// allocated to hold metadata for all founder chromosomes

	size_t m_UsedProgenyFndrAligns;			// number of actually used progeny to founder alignments
	size_t m_AllocdProgenyFndrAligns;			// number of allocated founder alignments
	size_t m_AllocdProgenyFndrAlignsMem;		// current mem allocation size for m_pProgenyFndrAligns
	tsProgenyFndrAligns*m_pProgenyFndrAligns;	// allocated to hold progeny allele stack alignments

	uint32_t m_UsedAlleleStacks;		// number of actually used allele stacks
	uint32_t m_AllocdAlleleStacks;		// number of allocated allele stacks
	size_t m_AllocdAlleleStacksMem;		// current mem allocation size for m_pAlleleStacks
	tsAlleleStack *m_pAlleleStacks;		// allocated to hold founder allele stacks
	uint32_t m_ProccessingBinID;         // most recent bin undergoing processing - processing may be incomplete
	uint32_t m_UsedHGBinSpecs;		    // number of actually used haplotype grouping bins
	uint32_t m_AllocdHGBinSpecs;		// number of allocated haplotype grouping bins
	size_t m_AllocdHGBinSpecsMem;		// current mem allocation size for m_pHGBinSpecs
	tsHGBinSpec* m_pHGBinSpecs;         // allocated to hold haplotype grouping bin specifications


	uint32_t m_UsedQGLLoci;		        // number of actually used QGLLoci
	uint32_t m_AllocdQGLLoci;		    // number of allocated QGLLoci
	size_t m_AllocdQGLLociMem;		        // current mem allocation size for m_pQGLLoci
	tsQGLLoci * m_pQGLLoci;             // allocated to hold haplotype grouping QGL loci


	uint32_t m_AllocWorkQueueEls;		// work queue allocated to hold at most this number of elements
	uint32_t m_TotWorkQueueEls;			// total number of work queue elements to be processed
#ifdef WIN32
	alignas(4) volatile uint32_t m_NumQueueElsProcessed;	// number of work queue elements processed
	alignas(4) volatile uint32_t m_FastSerialise;	// interlocked access to founder stack allocator -  - AcquireFastSerialise()
#else
	__attribute__((aligned(4))) volatile uint32_t m_NumQueueElsProcessed;	// number of work queue elements processed
	__attribute__((aligned(4))) volatile uint32_t m_FastSerialise;	// fast serialised access to founder stack allocator - AcquireFastSerialise()
#endif
	tsCHWorkQueueEl *m_pWorkQueueEls;	// thread work queue elements - currently really a list that is iterated over but in future may be adapted to become a queue

	uint32_t m_Fndrs2Proc[cMaxFounderReadsets];	// array of founders which are to be processed, indexed by FounderID-1. If LSB is set then that founder is marked for processing

	uint32_t m_InNumProcessed;	// this number of input bytes have been processed
	uint32_t m_InNumBuffered;	// currently buffering this many input bytes
	uint32_t m_AllocInBuff;		// m_pInBuffer allocated to hold this many input bytes
	uint8_t *m_pInBuffer;		// allocated for buffering input

	uint32_t m_OutBuffIdx;		// currently buffering this many output bytes
	uint32_t m_AllocOutBuff;	// m_pszOutBuffer allocated to hold this many output bytes
	uint8_t *m_pszOutBuffer;	// allocated for buffering output

	uint32_t m_WIGBuffIdx;		// currently buffering this many output bytes
	uint32_t m_AllocWIGBuff;	// m_pszOutBuffer allocated to hold this many output bytes
	uint8_t *m_pszWIGBuff;	    // allocated for buffering output

	int m_hInFile;				// input file handle
	int64_t m_InFileOfs;        // input file next read will start from this file offset
	int m_hOutFile;				// file handle for writing results
	int m_hWIGOutFile;			// file handle for writing results
	CCSVFile *m_pCSVFile;       // holds optional user specified haplotype group specifications in processing mode 3
	CCSVFile* m_pHGCSVFile;        // holds input previously generated haplotype grouping file for  post-processing

	char* m_pszRsltsFileBaseName;	// output results file basename
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

	tsCHWorkerInstance m_WorkerInstances[cMaxPBAWorkerThreads];	// to hold all worker instance thread parameters
	tsCHWorkerLoadChromPBAsInstance m_WorkerLoadChromPBAInstances[cMaxPBAWorkerThreads];	// to hold all worker instance thread parameters for loading chromosome PBAs

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
#else
	pthread_mutex_t m_hSerialiseAccess;
#endif
	void AcquireSerialise(void);
	void ReleaseSerialise(void);
	void AcquireFastSerialise(void);
	void ReleaseFastSerialise(void);

	// initialise and start pool of worker threads
	int	StartWorkerThreads(uint32_t NumThreads,		// there are this many threads in pool
							int32_t NumChroms);		// processing this number of chromosomes

	int StartWorkerLoadChromPBAThreads(uint32_t NumThreads,		// there are this many threads in pool
									   uint32_t StartSampleID,        // processing to start from this sample identifer
									   uint32_t EndSampleID,           // ending with this sample identifier inclusive
									   uint32_t ChromID);       // loading PBAs for this chromosome

	bool WaitAlignments(int WaitSecs);	// allow at most this many seconds for pool of worker threads to complete PBA scoring

	int	TerminateWorkerThreads(int WaitSecs = 120);	// allow at most this many seconds before force terminating threads

	tsHGBinSpec*    // returned haplotype grouping bin specification, returns NULL if all bins on chromosome have been iterated
		IterateHGBinSpecs(uint32_t PrevBinID,    // previously iterated bin identifier, to start from 1st bin then pass 0 as the bin identifier 
						  uint32_t ChromID=0); // only processing bins on this chrom, 0 if processing all bins

	// load and parse haplotype grouping bin specification file
	int								 // eBSFSuccess or error
		LoadHHGBinSpecs(char* pszInFile);  // load grouping specifications from this file

		// load and parse actual haplotype groupings previously generated by 'ngskit4b callhaplotypes'
	int								 // eBSFSuccess or error
		LoadHaplotypeGroupings(char* pszInFile);  // load previously generated haplotype groupings from this CSV file

	int32_t					// returned readset identifier (1..n) or < 1 if errors
		LoadPBAFile(char* pszFile,		// load chromosome metadata and PBA data from this file
					uint8_t ReadsetType,	// 0: founder, 1: progeny, 2: control
					bool bChromMetaOnly = false);  // if true then load chrom metadata (chrom name,length, file offset at which chrom PBAs start) but don't actually load the chromosome PBAs

	uint8_t *   // loaded PBAs or NULL if errors loading PBAs for SampleID.ChromID
		LoadSampleChromPBAs(uint32_t SampleID,   // Sample identifier
						uint32_t ChromID);    // chrom identifier specifying which PBAs is to be loaded from SampleID file

	bool  // false: no previous memory allocated for containing the chromosome PBAs, true: allocation was present and has been deleted
		DeleteSampleChromPBAs(uint32_t SampleID,   // Sample identifier
											   uint32_t ChromID);    // chrom identifier

	// trim 5' and 3' aligned segments within the PBAs attempting to reduce sequencing error induced false alleles
	int			// returns number of non-trimmed loci in the pPBAs
		TrimPBAs(uint32_t Trim5,	// trim 5' this many aligned PBA bases from each aligned segment
								  uint32_t Trim3,	// trim 3' this many aligned PBA bases from each aligned segment
								  uint32_t PBALen,	// pPBAs contains this many packed base alleles
								  uint8_t* pPBAs);	// packed base alleles to be processed

	uint8_t LocateReadsetChromLociAlleles(uint32_t ReadsetID,	// return alleles for this readset 
										uint32_t ChromID,		// on this chromosome
										uint32_t Loci);			// at this loci

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

	uint8_t *								// returned pointer to start of PBA
		LocatePBAfor(uint32_t ReadSetID,	// readset identifier 
			 uint32_t ChromID);				// chrom identifier

	tsCHChromMetadata *								// returned pointer to chromosome metadata
		LocateChromMetadataFor(uint32_t ReadSetID,		// readset identifier 
			 uint32_t ChromID);			// chrom identifier

	uint32_t								// number buffered
		FillInBuffer(uint32_t MinRequired,uint32_t MaxRequired =0); // try and fill input buffer with at least MinRequired or refill (if < MinRequired) up to MaxRequired (if 0 then max buffer allocation)


	uint8_t *AllocPBAs(uint32_t ChromLen);	// allocate to hold at least this many packed base alleles

	int AllocChromMetadata(void);			// allocate for additional chromosome metadata

	size_t									// returned index+1 into  m_pProgenyFndrAligns[] to allocated and initialised ProgenyFndrAligns, 0 if errors
		AddProgenyFndrAligns(tsProgenyFndrAligns* pInitProgenyFndrAligns);	// allocated tsProgenyFndrAligns to be initialised with a copy of pInitProgenyFndrAligns


	uint32_t									// returned index+1 into m_pAlleleStacks[] to allocated and initialised allele stack, 0 if errors
		AddAlleleStack(tsAlleleStack* pInitAlleleStack);	// allocated tsAlleleStack to be initialised with a copy of pInitAlleleStack

	int32_t									// returned index+1 into m_pHGBinSpecs[] to allocated and initialised allele stack, 0 if errors
	AddQGLLoci(tsQGLLoci* pInitQGLLoci);	// allocated tsQGLLoci to be initialised with a copy of pInitQGLLoci

	int32_t									// returned index+1 into m_pHGBinSpecs[] to allocated and initialised allele stack, 0 if errors
		AddQGLLoci(	uint32_t SrcExprID,  // sourced CSV haplotype group experiment identifier
					uint32_t SrcRowID,   // sourced CSV haplotype group row identifier
					uint32_t ChromID,    // QGL loci is on this chromosome
					uint32_t QGLLoci,    // at this loci
				    uint32_t NumHapGrps, // QGL has this number of actual haplotype groups
					uint32_t GrpID[4],   // group identifier for each allele accepted as having the highest F-measure, 0 if no group having accepted F-measure
					double FMeasure[4],  // F-measure for each allele - 0.0 if no group having an accepted F-measure for that allele
					tsGrpCnts GrpCnts[6]);    // counts for 1st 5 groups plus a pseudo-group containing sum of counts for groups after the 1st 5

	int32_t									// returned index+1 into m_pHGBinSpecs[] to allocated and initialised allele stack, 0 if errors
		AddHGBinSpec(tsHGBinSpec* pInitHGBinSpec);	// allocated tsHGBinSpec to be initialised with a copy of pInitHGBinSpec

	int32_t									// returned index+1 into m_pHGBinSpecs[] to allocated and initialised allele stack, 0 if errors
		AddHGBinSpec(uint32_t ChromID,	    // haplotype group clustering is on this chromosome
			uint32_t StartLoci,      // bin starts at this loci
			uint32_t NumLoci,	    // covering this many loci
			uint32_t MinCentroidDistance, // centroid distance can range from this minimum and up to MaxCentroidDistance when attempting to limit number of haplotype groupings to MaxNumHaplotypeGroups
			uint32_t MaxCentroidDistance, // maximum centroid distance which can be used when attempting to meet MaxNumHaplotypeGroups constraint
			uint32_t MaxNumHaplotypeGroups); // attempt to constrain number of haplotype groups to be this maximum



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

	
	int				// < 0 if errors, >=0 length of generated consensus
		GenConsensusPBA(uint32_t Loci,			// processing for allele stacks starting from this loci
			uint32_t NumLoci,		// processing for generation of this sized consensus PBAs
			uint32_t NumFndrs,		// number of founders to be processed 1st PBA at *pPBAs[0]
			uint8_t* pFounderPBAs[], // pPBAs[] pts to chromosome PBAs, for each of the chromosome founder PBAs
			uint8_t* pConsensusPBAs);     // ptr to preallocated sequence of length MaxNumLoci which is to be updated with consensus PBAs

	uint8_t				// concensus PBA for founders in a group of which specified founder is a member at the requested loci
		GenFounderConsensusPBA(uint32_t Founder, // consensus for all founders in same group as this founder
			tsHaplotypeGroup* pHaplotypes, // pts to haplotype grouping expected to contain the chrom and loci
			uint32_t ChromID,        // requiring consensus for this chrom at Loci
			uint32_t Loci,			// requiring consensus at this loci
			uint32_t NumFndrs,		// number of founders to be processed 1st PBA at *pPBAs[0]
			uint8_t* pFounderPBAs[]); // pPBAs[] pts to chromosome PBAs, for each of the chromosome founder PBAs

	uint8_t				// concensus PBA for founders in a group of which specified founder is a member at the requested loci
		GenGroupConsensusPBA(uint32_t GrpIdx, // consensus for all founders in this group
			tsHaplotypeGroup* pHaplotypes, // pts to haplotype grouping expected to contain the chrom and loci
			uint32_t ChromID,        // requiring consensus for this chrom at Loci
			uint32_t Loci,			// requiring consensus at this loci
			uint32_t NumFndrs,		// number of founders to be processed 1st PBA at *pPBAs[0]
			uint8_t* pFounderPBAs[]); // pPBAs[] pts to chromosome PBAs, for each of the chromosome founder PBAs


	int32_t				// error or success (>=0) code 
		ProcessGrpLociQGLs(uint32_t MinQGLGrpMembers,        // groups with fewer than this number of members are treated as noise alleles and these groups are not used when determining QGL group specific major alleles
						   double MinQGLGrpPropTotSamples,  // groups containing less than this proportion of all samples are treated as if containing noise and alleles in these groups are not used when determining QGL group specific major alleles 
						   double MinQGLFmeasure,           // only accepting QGL loci with at least this F-measure score
						   char* pszHapGrpFile,             // input, previously generated by 'callhaplotypes', haplotype group file (CSV format) 
						   uint32_t NumSamplePBAsInputFiles,	     // number of input sample PBAs file specs
						   char* pszSamplePBAsInputFiles[],	 // names of input sample PBAs files (wildcards allowed)
						   char* pszOutFile);		         // write haplotype group loci QGLs to this output file (CSV format)


	int32_t				// error or success (>=0) code 
		GenBinQGLs(uint32_t ChromID,          // requiring QGLs across this chrom
			uint32_t NumFndrs,		  // number of founders to be processed 1st PBA at *pPBAs[0]
			uint8_t* pFounderPBAs[]); // pPBAs[] pts to chromosome PBAs, for each of the chromosome founder PBAs


	tsHaplotypeGroup * // returns allocated ptr to haplotype groups - allocated with 'new', caller must free with 'delete'
		GroupHaplotypes(uint32_t ChromID, // grouping is on this chromosome
			          uint32_t StartLoci, // grouping starts at this loci
	                     uint32_t Length,     // grouping is over this many loci
			             int32_t MinCentroidDistance, // centroid distance can range from this minimum and up to MaxCentroidDistance when attempting to limit number of haplotype groupings to MaxNumGroups
			             int32_t MaxCentroidDistance, // maximum centroid distance which can be used when attempting to meet MaxNumGroups constraint
			             uint32_t MaxNumGroups, // attempt to constrain number of haplotype groups to be this maximum
			             uint32_t* pFndrDiffs, // matrix counts of founder differentials
			             uint32_t NumFndrs); // number of founders

	int				// < 0 if errors, otherwise success
		GenHaplotypeGroups(tsHGBinSpec* pHGBinSpec,   // clustering using these specs
			uint32_t NumFndrs,		// number of founders to be processed 1st PBA at *pPBAs[0]
			uint8_t* pFounderPBAs[]);	// pPBAs[] pts to chromosome PBAs, for each of the chromosome founder PBAs

	int
		AlignAlleleStacks(uint32_t NumFndrs,			// number of founders to be processed against
					uint32_t MaxNumLoci = 100000);		// comparing ReadsetID's PBA against founder PBAs then work queue items specify this number of loci for processing by threads

	int
		AlignFounderHaps(uint32_t NumFndrs);			// number of founders to be processed against

	int GenAlleleStacks(uint32_t NumFndrs,		// number of founders to be processed against
						uint32_t MaxNumLoci = 100000);		// work queue items specify this number of loci for processing by threads

	int GenFounderHaps(uint32_t NumFndrs);		// number of founders to be processed against)

	int	ReportMatrix(char* pszRsltsFileBaseName,		// matrix results are written to this file base name with '.SummaryMatrix.csv' appended
		bool bRaw = false);								// true if reporting on raw haplotypes before any imputing/filtering

	int ImputeOutliersHaplotypes(int32_t MaxDistance,	// imputing outliers from other called haplotypes which are within MaxDistance
							uint32_t ReadsetID);			// report on this progeny readset only, or if 0 then report on all progeny readsets

	int
		ReduceProgenyFounders(uint32_t MaxDistance, // where number of founders is more than 2 then attempt to reduce down to a most 2 at any given progeny loci
								uint32_t ReadsetID);				// reduce this progeny readset only, or if 0 then all progeny readsets

	int
		ImputeProgenyHeterozygosity(int32_t MaxDistance, // imputing regional heterozygostic regions having apparent high rates of single haplotype sampling which are within MaxDistance
								uint32_t ReadsetID);				// report on this progeny readset only, or if 0 then report on all progeny readsets

	int ReportHaplotypesByProgeny(char* pszRsltsFileBaseName,		// haplotype results are written to this file base name with 'ProgenyID.csv' appended
								  uint32_t ReadsetID,				// report on this progeny readset only, or if 0 then report on all progeny readsets
	                              bool bRaw = false);								// true if reporting on raw haplotypes before any imputing/filtering

	bool ValidatePBA(uint8_t* pAlleles,		// validate that PBA alleles are properly conformant
		bool bSetNoAlleles = true, // if non-conformant then overwrite *pAlleles to be no alleles present
		bool bNormalise = true);    // normalise for coverage, low coverage major alleles (2,0,0,0) are normalised to high coverage (3,0,0,0)

	int			// returned number of PBAs which are non-conformant
		ValidatePBAs(uint32_t Length, uint8_t* pPBAs, // validate sequence of PBAs
			bool bSetNoAlleles = true, // if non-conformant then overwrite *pAlleles to be no alleles present
			bool bNormalise=true);    // normalise for coverage, low coverage major alleles (2,0,0,0) are normalised to high coverage (3,0,0,0)


	void	BitsVectSet(uint16_t Bit,		// bit to set, range 0..4095
				   ts4096Bits& BitsVect);
	void	BitsVectReset(uint16_t Bit,		// bit to reset, range 0..4095
					 ts4096Bits& BitsVect);
	bool BitsVectEqual(ts4096Bits & BitsVectA,	// compare for equality
								ts4096Bits & BitsVectB);
	bool BitsVectTest(uint16_t Bit,		// bit to test, range 0..4095
					ts4096Bits& BitsVect);
	void BitsVectInitialise(bool Set,			// if true then initialse all bits as set, otherwise initialise all bits as reset
										   ts4096Bits& BitsVect);
	uint32_t BitsVectCount(ts4096Bits& BitsVect);		// count number of set bits

	uint32_t BitsVectUnion(ts4096Bits& BitsVectA, ts4096Bits& BitsVectB);		// union (effective BitsVectA | BitsVectB) of bits in BitsVectA with BitsVectB with BitsVectA updated, returns number of set bits in BitsVectA  
	uint32_t BitsVectIntersect(ts4096Bits& BitsVectA, ts4096Bits& BitsVectB);	// intersect (effective BitsVectA & BitsVectB) of bits in BitsVectA with BitsVectB with BitsVectA updated, returns number of set bits in BitsVectA  
	uint32_t BitsVectClear(ts4096Bits& BitsVectA, ts4096Bits& BitsVectB);	    // clear bits in BitsVectA which are set in BitsVectB with BitsVectA updated, returns number of set bits in BitsVectA  

	CMTqsort m_mtqsort;				// multi-threaded qsort
	static int SortAlleleStacks(const void* arg1, const void* arg2);
	static int SortProgenyFndrAligns(const void* arg1, const void* arg2);
	static int SortProgenyChromLociReadset(const void* arg1, const void* arg2);
	static int SortHGBinSpecs(const void* arg1, const void* arg2);
	static int SortQGLLociFMeasureDesc(const void* arg1, const void* arg2);
	static int SortQGLLoci(const void* arg1, const void* arg2);

	int32_t ReportHaplotypesAsGWAS(char* pszRsltsFileBaseName,		// generate a GWAS format file(s) for viewing haplotype calls in IGV, useful for looking at associations with coverage etc
								   uint32_t ReadsetID,				// report on this progeny readset
								   bool bRaw = false);				// true if reporting on raw haplotypes before any imputing/filtering

	int32_t ReportAnchorsAsGWAS(char* pszRsltsFileBaseName);		// generate a GWAS format file(s) for viewing marker anchors calls in IGV, useful for looking at associations with coverage etc
	int32_t ReportAnchorsAsCSV(char* pszRsltsFileBaseName);		// generate a CSV format file(s) containing anchor (potential markers) loci

	int32_t ReportHaplotypeGroups(uint32_t NumHaplotypes);         // number of haplotypes which have been grouped 

	int	LoadChromSizes(char* pszBEDFile); // BED file containing chromosome names and sizes

	int  CSV2WIG(char* pszInFile, // input file containing haplotype groupings
			char* pszOutFile); // where to write as WIG formated file

	int  // simulating a in-memory PBA with PBA loci alleles replaced by coverage obtained from a WIG coverage file generated at the same time as the PBA file - enables haplotype coverage grouping
		LoadPBACoverage(char* pszInPBA);   // file containing PBA, file name extension will be replaced with 'coverage.wig' which will be expected to be name of file containing the WIG format coverage

	uint8_t *LoadPBAChromCoverage(uint32_t ReadsetID, // loading chromosome coverage for this readset and mapping coverage as though PBAs
						uint32_t ChromID);    // coverage is on this chromosome

	void InitialiseWIGSpan(void);				// initialise WIG span vars to values corresponding to no spans having been previously reported
	int CompleteWIGSpan(bool bWrite = false);	// close off any current WIG span ready to start any subsequent span, if bWrite is true then write to disk 
	int AccumWIGCnts(uint32_t ChromID,	// accumulate wiggle counts into variableStep spans over this chromosome
			uint32_t Loci,		// span is starting at this loci
		    uint32_t BinLen,       // span bin is this length
			uint32_t Cnts,		// has this many counts attributed to the bin
			uint32_t MaxSpanLen = 10000000); // allow reported variableStep WIG spans to be continuous up this maximal length and then start another span

	CStats m_Stats;					// used for determining statistical significance of individual founder unique counts relative to other founder counts

public:
	CCallHaplotypes();
	~CCallHaplotypes();
	void Reset(void);	// resets class instance state back to that immediately following instantiation
	int Process(eModeCSH PMode,	// processing mode 0: report imputation haplotype matrix, 1: report both raw and imputation haplotype matrices, 2: additionally generate GWAS allowing visual comparisons, 3: allelic haplotype grouping,4: coverage haplotype grouping, 5: post-process haplotype groupings for QGLs,  6: post-process to WIG
			int32_t ExprID,			         // assign this experiment identifier for this PBA analysis
			int32_t SeedRowID,                  // generated CSVs will contain monotonically unique row identifiers seeded with this row identifier  
			int32_t LimitPrimaryPBAs,        // limit number of loaded primary or founder PBA files to this many. 0: no limits, > 0 sets upper limit
			int32_t GrpHapBinSize,           // when grouping into haplotypes (processing mode 3) then use this non-overlapping bin size for accumulation of differentials between all founder PBAs
	        int32_t NumHapGrpPhases,         // number of phases over which to converge group consensus when haplotype clustering into groups - moves the balance between optimal group consensus and processing resource
			int32_t MinCentClustDist,        // haplotype groupings - in processing mode 3/4 only - minimum group centroid clustering distance
			int32_t MaxCentClustDist,        // haplotype groupings - in processing mode 3/4 only - maximum group centroid clustering distance
			int32_t MaxClustGrps,            // haplotype groupings - in processing mode 3 only - targeted maximum number of groups
			uint32_t SparseRepPropGrpMbrs,	// only apply sparse representative imputation if number of haplotype group members at least this 
			double SparseRepProp,			// if highest frequency (consensus) allele is 0x00 (no coverage) then if next highest frequency allele is at least this proportion of all members in haplotype group then treat next highest allele as
										// being as being the consensus. Objective to obtain an imputed allele in regions of sparse haplotype coverage
			int32_t MaxReportGrpQGLs,     // when calling group QGLs then, if non-zero - report this many highest scoring QGLs
			int32_t MinQGLGrpMembers,        // groups with fewer than this number of members - in processing mode 5 only -are treated as noise alleles these groups are not used when determining QGL group specific major alleles
			double MinQGLGrpPropTotSamples,  // haplotype groups - in processing mode 5 only -containing less than this proportion of all samples are treated as if containing noise and alleles in these groups are not used when determining QGL group specific major alleles 
			double MinQGLFmeasure,           // only accepting QGL loci with at least this F-measure score
			int32_t FndrTrim5,			// trim this many aligned PBAs from 5' end of founder aligned segments - reduces false alleles due to sequencing errors
			int32_t FndrTrim3,			// trim this many aligned PBAs from 3' end of founder aligned segments - reduces false alleles due to sequencing errors
			int32_t ProgTrim5,			// trim this many aligned PBAs from 5' end of progeny aligned segments - reduces false alleles due to sequencing errors
			int32_t ProgTrim3,			// trim this many aligned PBAs from 3' end of progeny aligned segments - reduces false alleles due to sequencing errors
			int32_t WWRLProxWindow,		// proximal window size for Wald-Wolfowitz runs test
			int32_t OutliersProxWindow,	// proximal window size for outliers reduction
			char *pszMaskBPAFile,	// optional masking input BPA file, only process BPAs which are intersect of these BPAs and progeny plus founder BPAs
			char* pszChromFile,			// BED file containing reference assembly chromosome names and sizes
			int NumFounderInputFiles,	// number of input founder file specs
			char *pszFounderInputFiles[],	// names of input founder PBA files (wildcards allowed)
			int NumProgenyInputFiles,	// number of input progeny file specs
			char* pszProgenyInputFiles[],		// names of input progeny PBA files (wildcards allowed)
			char* pszOutFile,		// loci haplotype calls output file (CSV format)
			int	NumIncludeChroms,			// number of chromosome regular expressions to include
			char **ppszIncludeChroms,		// array of include chromosome regular expressions
			int	NumExcludeChroms,			// number of chromosome expressions to exclude
			char **ppszExcludeChroms,		// array of exclude chromosome regular expressions
			int NumThreads);		// number of worker threads to use

	int ProcWorkerThread(tsCHWorkerInstance* pThreadPar);	// worker thread parameters
	int ProcWorkerLoadChromPBAThread(tsCHWorkerLoadChromPBAsInstance* pThreadPar);	// worker thread parameters

};


