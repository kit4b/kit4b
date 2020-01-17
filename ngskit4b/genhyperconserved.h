#pragma once

const int cMAtotMaxSeqAlignLen = 0x0fffffff; // total (over all aligned species) max seq length that can be buffered in concatenated seqs

const int cMinCoreLen = 1;				// allow core lengths to be specified down to cMinCoreLen
const int cDfltMinCoreLen = 50;			// if core lengths not specified then default to cDfltMinCoreLen
const int cMaxMinCoreLen = 10000;		// minimum core lengths can be specified up to this length

const int cMinIdentity = 50;			// accept down to 50% identity
const int cMaxIdentity = 100;			// can't do any better than 100% identity!
const int cDfltIdentity = 90;			// use 90% identity if non specified and hypercores are to be processed
const int cRandWalk100Score = 10000;	// random walk score for 100% identity

const int cMaxMismatches = 500;			// total number of mismatches allowed in any hypercore before terminating that hypercore
const int cDfltMaxMismatches = 100;		// default if non specified

const int cMinMismatchHistLen = 10;		// minimal sized window over which mismatch history is maintained
const int cMaxMismatchHistLen = 500;	// maximal sized window over which mismatch history is maintained

const int cMaxIncludeFiles = 10;		// maximum number of include region filter files
const int cMaxExcludeFiles = 10;		// maximum number of exclude region filter files

const int cMaxIncludeChroms = 20;		// max number of include chromosomes regular expressions
const int cMaxExcludeChroms = 20;		// max number of exclude chromosomes regular expressions

const int cMinMatchDistSegments = 4;	// min match distribution profile segments (applies if eProcModeOutspecies)
const int cDfltMatchDistSegments = 10;	// default match distribution profile segments (applies if eProcModeOutspecies)
const int cMaxMatchDistSegments = 100;	// max match distribution profile segments (applies if eProcModeOutspecies)

typedef enum eProcMode {
	eProcModeStandard = 0,				// default processing
	eProcModeSummary,					// summary processing
	eProcModeOutspecies					// same as default but with outspecies species additional detail
} etProcMode;

typedef struct TAG_sLenRangeClass {
	int ID;					// uniquely identifies this range
	int Min;				// minimum length in this range
	int Max;				// maximum length in this range
	char szDescr[10];			// descriptive text
} tsLenRangeClass;


typedef struct TAG_sExcludeSpeciesChrom {
	tSpeciesID SpeciesID;		// which species
	tChromID ChromID;			// chromosome is not to be processed
} tsExcludeSpeciesChrom;


typedef struct TAG_sDistSeg {
	int Matches;			// number of exact matches in this segment
	int Mismatches;			// number of mismatches
	int InDels;				// total number of InDel bases
	int Unaligned;			// number of unaligned bases
} tsDistSeg;

typedef struct TAG_sProcParams
{
	int ProcMode;					// processing mode 0: default, 1: summary stats, 2: outspecies processing
	bool bStatsAvail;
	char* pszSpeciesList;			// comma separated species list starting with reference species
	int NumSpeciesList;				// number of species in species list
	int RefSpeciesIdx;				// current index into pSeqs for the reference species
	char szSpecies[cMaxAlignedSpecies][cMaxDatasetSpeciesChrom];	// species names of interest - other species are sloughed
									// only alignments with species included will be processed
									// first species is the reference species
	int NumCoreSpecies;				// number of core species (in species list priority order) required in an alignment
	int MinAlignSpecies;			// minimum number of species required in an alignment - will be >= NumCoreSpecies
	int MaxNumStatsSpecies;			// max number of species to accumulate stats for
	int NumSpeciesInAlignment;		// actual number of sequences in current alignment
	int MaxAlignIdxSpecies;			// pSeqs[MaxAlignIdxSpecies-1] is last actual alignment sequence
	char szRefChrom[cMaxDatasetSpeciesChrom]; // current reference chromosome
	int RefChromID;					// current reference chromosome identifier
	bool bFiltLoConfidence;			// true if low confidence subseqences to be filtered out
	bool bFilt1stLast;				// true if treat 1st and last subsequences as being low confidence
	int MinIdent;					// treat subsequences of less than this identity as being low confidence
	int MinSubSeqLen;				// subsequences of less than this length are treated as being low confidence
	bool bInDelsAsMismatches;		// treat InDels as if mismatches (one mismatch == one InDel column)
	bool bSloughRefInDels;			// slough columns in which the reference base is an InDel
	int	WindowSize;					// sampling window size
	int NxtOutputOffset;			// when to next output results - i.e end of current window
	int NumDistSegs;				// number of match distribution profile segments
	int SeqOfs;						// offset into sequences
	int SeqLen;						// sequence length
	etSeqBase* pSeqs[cMaxAlignedSpecies];  // ptrs to each sequence
	int MaxSeqAlignLen;				// max length alignment which can be processed (how much mem was alloc'd to pSeq[n])			
	int MinHyperLen;				// minimum hyper core length required
	int MinUltraLen;				// minimum ultra core length required
	int MinCoreLen;					// maximum of either MinHyperLen or MinUltraLen
	int MaxHyperMismatches;			// hyper cores can have upto this number of total mismatches
	int VMismatches;				// number of mismatches in an alignment col to count as a missmatch against MaxMismatches
	int* pCntStepCnts;				// array of stats counters
	int NumCnts;					// number of steps in pCntStepCnts
	int Regions;					// number of regions per step
	CBEDfile* pBiobed;				// if not NULL then opened biobed file for regional characteristics
	int BEDChromID;					// BED chromosome identifier corresponding to RefChromID
	int NumIncludes;				// number of biobed files containing regions to include
	int NumExcludes;				// number of biobed files containing regions to exclude
	CBEDfile* pIncludes[cMaxIncludeFiles];	// if opened biobed files for regions to include - all other regions are to be excluded
	CBEDfile* pExcludes[cMaxExcludeFiles];	// if opened biobed files for regions to exclude 
	int UpDnStreamLen;				// up/dn stream regional length when characterising
	bool bMultipleFeatBits;			// if false then stats only generated if a single feature bit is set - e.g if both exons and introns overlapped then no stat generated
	int hRsltsFile;					// write stats results into this CSV file
	int hCoreCSVRsltsFile;			// write hypercore loci into this CSV file
	int NumIncludeChroms;			// number of chromosomes explicitly defined to be included
	char** ppszIncludeChroms;		// ptr to array of reg expressions defining chroms to include - overides exclude
	int NumExcludeChroms;			// number of chromosomes explicitly defined to be excluded
	char** ppszExcludeChroms;		// ptr to array of reg expressions defining chroms to include
#ifdef _WIN32
	Regexp* IncludeChromsRE[cMaxIncludeChroms];	// compiled regular expressions
	Regexp* ExcludeChromsRE[cMaxExcludeChroms];
#else
	regex_t IncludeChromsRE[cMaxIncludeChroms];	// compiled regular expressions
	regex_t ExcludeChromsRE[cMaxExcludeChroms];
#endif

	bool bAllUniqueSpeciesSeqs;		// true if all species sequences must not overlap any other sequence
	int NumUniqueSpeciesSeqs;		// number of species in UniqueSpeciesSeqs
	char UniqueSpeciesSeqs[cMaxAlignedSpecies][cMaxDatasetSpeciesChrom];	// species for which sequences in alignment blocks must not overlap with sequence in any other block


	int MinIdentity;				// minimum required identity (50-100%) from which MismatchScore and MatchScore were derived  
	int MismatchScore;				// score to reduce random walk score by for mismatches - processing terminates if score drops below 0
	int MatchScore;					// score to increase random walk score by for matches - random walk score limited to cRandWalk100Score;
} tsProcParams;

int ParseNumSpecies(char* pszSpeciesList, tsProcParams* pProcParams);
int ParseUniqueSpeciesSeqs(char* pszUniqueSpeciesSeqs, tsProcParams* pProcParams);

bool ProcAlignBlock(int RefChromID, int RefChromOfs, int AlignLen, tsProcParams* pProcParams);
bool ProcAlignBlockSummary(int RefChromID, int RefChromOfs, int AlignLen, tsProcParams* pProcParams);
bool ProcessAlignment(int RefChromID, int ChromOfs, int SeqIdx, int SubSeqLen, tsProcParams* pProcParams);
int ProcessSubSeq(int RefChromID, int ChromOfs, int SeqIdx, int MaxLen, tsProcParams* pProcParams);
bool ChkOutputResults(const char* pszChrom, int ChromOffset, tsProcParams* pProcParams, bool bGenEmptyRows, bool bGenEmptyWindows);
bool OutputResults(const char* pszChrom, int ChromOffset, tsProcParams* pProcParams, bool bGenEmptyRows);
bool ChkOutputSummaryResults(const char* pszChrom, int ChromOffset, tsProcParams* pProcParams, bool bGenEmptyRows, bool bGenEmptyWindows);
bool OutputSummaryResults(const char* pszChrom, int ChromOffset, tsProcParams* pProcParams, bool bGenEmptyRows);
bool OutputHypercore(const char* pszChrom, int ChromStartOffset, int ChromEndOffset, int FeatureBits,
	int OGUnaligned, int OGMatches, int OGMismatches, int OGInDels, tsDistSeg SegCnts[], tsProcParams* pProcParams);
char* ChkSpeciesChromWellFormed(char* pszSpeciesChroms);
int ReportSummary(int RefChromID, int RefChromOfs, int ProcMode, tsProcParams* pProcParams);
bool IncludeFilter(int RefChromID, int SubRefOfs, int SubRefEndOfs, tsProcParams* pProcParams);

int										// returned blockid to next start loading from
LoadContiguousBlocks(int RefSpeciesID,	// reference species identifier
	int  BlockID,			// which block to initially start loading from
	bool* pbLoaded,			// returned indicator as to if any loaded blocks meet processing requirements
	int* pRefChromID,		// returned reference chromosome identifier 
	char* pRefStrand,		// returned reference strand
	int* pRefAlignLen,		// returned alignment (incl InDels) length
	int* pRefChromOfs,		// returned alignment start offset
	int* pSpeciesIDs,		// input - species of interest identifier array
	CMAlignFile* pAlignments,
	tsProcParams* pProcParams);

