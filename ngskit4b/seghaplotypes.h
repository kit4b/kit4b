#pragma once

const int cMaxSHLenPrefix = 10;							// founder or Tag name prefixes must be no longer than this number of chars
const char cTagSHTerm1 = '|';							// tags are terminated by 2 chrs, this is the first. Preceding tag name must be alpha-numeric only and a max of cMaxSHLenPrefix in length
const char cTagSHTerm2 = '#';							// and this is the second

const int cMaxSHFounders = 5;							// allowing for at most this many founders (quad-parental + 1)
const int cAllocSHBuffInSize = 0x3fffffff;				// allocating buffering for input file
const int cAllocSHBuffOutSize = cAllocSHBuffInSize;		// allocating buffering for output files
const int cAllocSHNumSAMloci = 0x0ffffff;				// allocate/realloc for this many SAMloci
const int cMaxSHSeqNames = 0x0ffffff;					// can accept at most this many unique sequence names

const int cMinSHBinSize = 1;							// counting number of alignment read loci within sliding window (smoothing) of this Kbp size 
const int cDfltSHBinSize = 10;							// counting number of alignment read loci within sliding window (smoothing) of this Kbp size 
const int cMaxSHBinSize = 1000;							// counting number of alignment read loci within sliding window (smoothing) of Kbp this size

const int cBEDNoScore = 100;							// if not scoring bins then use this value as the score

typedef enum TAG_eModeSH {
	eMSHDefault = 0,		// default is to generate segmentations using bin counts of unique loci only
	eMSHSegAll,			// generate segmentations using bin counts of all alignment loci incl non-unique
	eMSHPlaceHolder			// used to mark end of processing modes
	}eModeSH;

#pragma pack(1)

typedef struct TAG_tsSHSAMloci
	{
	int TargID;			// target sequence name identifier (founder tag has been removed)
	int FounderID;		// founder tag identifier
	uint32_t TargLoci;	// where alignment starts on target - no adjustments for trimming, loci is that recorded in SAM record
	uint32_t Cnt;		// number of alignments at this loci
	} tsSHSAMloci;

typedef struct TAG_tsSHBin
	{
	uint32_t TargID;		// target sequence name identifier
	uint32_t TargLoci;		// bin starts at this loci
	uint32_t BinLen;		// bin is this length
	uint8_t fCalled:1;		// set if haplotypes for this bin have been called
	uint8_t fInfer:1;		// haplotype call was inferenced from an adjacent bin
	uint32_t CalledHaplotype[cMaxSHFounders]; // called haplotype score
	uint32_t RawCnts[cMaxSHFounders];		// founder raw counts, indexed by FounderID - 1 
	uint32_t SmoothedCnts[cMaxSHFounders];	// founder smoothed counts, indexed by FounderID - 1
	} tsSHBin;

typedef struct TAG_tsTargSeq
	{
	uint32_t fAligned:1;	// flags that there has been at least one accepted alignment to this target sequence
	uint32_t TargSeqID;		// target sequence name identifier
	uint32_t TargSeqLen;	// target sequence length
	uint32_t TargSeqNameOfs; // offset into m_szSeqNames at which target sequence name starts
	uint32_t NumBins;		 // this many bins allocated to hold counts for alignments to this target sequence
	uint32_t BinsOfs;		 // allocated bins start at this offset in m_pBins[]
	} tsTargSeq;

#pragma pack()

class CSegHaplotypes
{
	uint32_t m_BinSizeKbp;			// each bin is at most this many Kbp in size
	uint32_t m_AllocdBins;			// number of allocated bins
	tsSHBin *m_pBins;				// allocated to hold bins for maximal sized targeted sequence

	uint32_t m_NumAlignedTargSeqs;	// number of target sequences which were aligned to at least once
	uint32_t m_LASeqNameID;			// sequence name identifier last returned by AddTargSeqName()
	uint32_t m_NumSeqNames;			// number of sequence names currently in m_szSeqNames
	uint32_t m_NxtszSeqNameIdx;		// current concatenated (names separated by '\0') of all chromosome names in m_szSeqNames
	char m_szSeqNames[cMaxSHSeqNames * (cMaxDatasetSpeciesChrom/2)];	// used to hold concatenated chromosome names, each separated by '\0'
	tsTargSeq m_TargSeqs[cMaxSHSeqNames];	// one entry for each target sequence

	uint32_t m_LAFounderID;			// name identifier last returned by AddFounderName()
	uint32_t m_NumFounders;			// number of founder names currently in m_szFounders
	uint32_t m_NxtszFounderIdx;		// current concatenated (names separated by '\0') of all founder names in m_szFounders
	char m_szFounders[cMaxSHFounders * (cMaxSHLenPrefix + 1)];	// used to hold concatenated founder names, each separated by '\0'
	uint32_t m_szFounderIdx[cMaxSHFounders];	// array of indexes into m_szFounders giving the starts of each founder name

	uint32_t m_InBuffIdx;	// currently buffering this many input bytes
	size_t m_AllocInBuff;	// m_pInBuffer allocated to hold this many input bytes
	uint8_t *m_pInBuffer;	// allocated for buffering input

	uint32_t m_OutBuffIdx;		// currently buffering this many output bytes
	uint32_t m_AllocOutBuff;	// m_pOutBuffer allocated to hold this many output bytes
	uint8_t *m_pOutBuffer;		// allocated for buffering output

	int m_hOutFile;			// output file handle

	CSAMfile *m_pSAMfile;			// processing alignments from SAM/BAM file 
	tsBAMalign *m_pBAMalignment;	// SAM/BAM alignment record

	size_t m_CurNumSAMloci;			// number of loci currently accepted
	size_t m_AllocdSAMloci;			// m_pSAMlociMem can hold at most this many  
	size_t m_AllocdSAMlociMem;		// allocated memory for holding accepted alignment loci
	tsSHSAMloci *m_pSAMloci;		// pts to memory dynamically allocated as needed to hold all accepted alignment loci

	int ApplySmoothing(void);	// Cnts in immediately adjacent bin counts are weighted at 0.5 and these weighted counts added to the current bin cnts as smoothed bin counts


	int								// total number of bins which have been assigned a haplotype
		IdentifySegments(int MinCoverage,		// requiring at least this min number of counts to make call for any bin
						 bool bInterpolate = false, // if true - if bin counts less than minimum coverage then interpolate from adjacent bins which were called
						 bool bDontScore = false);	// if true then mark bin haplotype with score of cBEDNoScore without actually scoring according to bin counts

	int genBED(char *pszTrackName,		// track name
				char *pszTrackDescr);	// track descriptor


	int	
		GenBinnedSegments(eModeSH PMode,			// processing mode
			char *pszTrackName,		// track name
			char *pszTrackDescr,	// track descriptor
			bool bDontScore,			// don't score segment bins
			uint32_t BinSizeKbp,	// Wiggle score is number of alignments over this sized bins
			char* pszInFile,		// alignments are in this SAM/BAM file 
			char* pszOutFile);		// write out BED to this file

	int		// bin score assigned 
		GenBinnedCoverage(int TargID,	// coverage is on this targeted chrom/seq
				uint32_t BinStart,		// coverage starts at this loci inclusive
				uint32_t BinEnd,		// coverage ends at this loci inclusive
				uint32_t BinCnts);		// total number of coverage counts

	uint32_t							// returned sequence identifier, 0 if unable to accept this sequence name
		AddTargSeqName(char* pszSeqName, // associate unique identifier with this sequence name
						int SeqLen);	 // sequence is this length

	char*							// returned sequence name
		LocateTargSeqName(uint32_t SeqID);	// identifier returned by call to AddTargSeqName

	tsTargSeq *							// returned sequence detail
		LocateTargSeq(uint32_t SeqNameID);	// identifier returned by call to AddTargSeqName

	tsTargSeq *							// returned sequence detail
		LocateTargSeq(char *pszSeqName);	// sequence name

	uint32_t		// returned founder name identifier, 0 if unable to accept this founder name
		AddFounder(char* pszFounder); // associate unique identifier with this founder name

	char*							// returned founder name
		LocateFounder(uint32_t FounderID);// identifier returned by call to AddFounderName

	uint32_t		// returned founder name identifier, 0 if unable to accept this founder name
		LocateFounder(char* pszFounder); // associate unique identifier with this founder name

	
	int			// returned len of parsed founder name, 0 if unable to parse out a founder (add 2 to len to obtain start of non-prefixed original string 
		ParseFounder(char* pszIn);			// input null terminated string assumed to contain founder name within first cMaxSHLenPrefix chars


	// SortSAMTargLoci
// Sort m_pSAMloci by ascending Targ.Loci
	static int SortSAMTargLoci(const void* arg1, const void* arg2);


// Sort m_pSAMloci by ascending Founder.Targ.Loci
	static int SortSAMFounderTargLoci(const void* arg1, const void* arg2);

public:
	CSegHaplotypes();
	~CSegHaplotypes();
	void Reset(void);	// resets class instance state back to that immediately following instantiation

	int Process(eModeSH PMode,			// processing mode
			char *pszTrackName,			// track name
			char *pszTrackDescr,		// track descriptor
			bool bDontScore,				// don't score segment bins
			int BinSizeKbp,				// using this sized bins
			char* pszInFile,			// input SAM file
			char* pszOutFile);			// output to this file
};


