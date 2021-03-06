#pragma once

const int cMaxLenPrefix = 10;							// prefixes must be no longer than this number of chars
const int cAllocPGBuffInSize = 0x3fffffff;				// allocating buffering for input file
const int cAllocPGBuffOutSize = cAllocPGBuffInSize;		// allocating buffering for output files
const int cAllocAsciiChkSize = 0x03fffff;				// check if first 4MB of input file is ascii - only ascii files are parsed
const int cAllocNumSAMloci = 0x0ffffff;					// allocate/realloc for this many SAMloci
const int cMaxSeqNames = 0x0ffffff;						// can accept at most this many unique sequence names

const int cMinWiggleBinSize = 1;					// Wiggle scores are generated by counting number of alignment read loci within sliding window (smoothing) of this Kbp size 
const int cDfltWiggleBinSize = 10;				// Wiggle scores are generated by counting number of alignment read loci within sliding window (smoothing) of this Kbp size 
const int cMaxWiggleBinSize = 1000;					// Wiggle scores are generated by counting number of alignment read loci within sliding window (smoothing) of Kbp this size

typedef enum TAG_eModePG {
	eMPGDefault = 0,		// default is to prefix fasta descriptors with originating genome identifier
	eMPGFilterPrefix,		// filter SAM alignments removing alignments to targets which are not prefixed
	eMPGWiggleUniqueLoci,	// generate UCSC Wiggle file using bin counts of unique loci only
	eMPGWiggleAll,			// generate UCSC Wiggle file using bin counts of all alignment loci incl non-unique
	eMPGPlaceHolder			// used to mark end of processing modes
	}eModePG;

#pragma pack(1)

typedef struct TAG_tsSAMloci
	{
	int TargID;		// target name identifier
	uint32_t TargLoci;	// where alignment starts on target - no adjustments for trimming, loci is that recorded in SAM record
	uint32_t Cnt;		// number of alignments at this loci
	} tsSAMloci;

#pragma pack()

class CPangenome
{

	int m_LASeqNameID;			// sequence name identifier last returned by AddTargSeqName()
	int m_NumSeqNames;			// number of sequence names currently in m_szSeqNames
	int m_NxtszSeqNameIdx;		// current concatenated (names separated by '\0') of all chromosome names in m_szSeqNames
	char m_szSeqNames[cMaxSeqNames * (cMaxDatasetSpeciesChrom/2)];	// used to hold concatenated chromosome names, each separated by '\0'
	uint32_t m_szSeqNameIdx[cMaxSeqNames];		// array of indexes into m_szSeqNames giving the starts of each sequence name

	uint32_t m_InBuffIdx;	// currently buffering this many input bytes
	size_t m_AllocInBuff;	// m_pInBuffer allocated to hold this many input bytes
	uint8_t *m_pInBuffer;	// allocated for buffering input

	uint32_t m_OutBuffIdx;		// currently buffering this many output bytes
	uint32_t m_AllocOutBuff;	// m_pOutBuffer allocated to hold this many output bytes
	uint8_t *m_pOutBuffer;		// allocated for buffering output

	int m_hInFile;			// input file handle
	int m_hOutFile;			// output file handle

	CSAMfile *m_pSAMfile;			// processing alignments from SAM/BAM file 
	tsBAMalign *m_pBAMalignment;	// SAM/BAM alignment record

	size_t m_CurNumSAMloci;			// number of loci currently accepted
	size_t m_AllocdSAMloci;			// m_pSAMlociMem can hold at most this many  
	size_t m_AllocdSAMlociMem;		// allocated memory for holding accepted alignment loci
	tsSAMloci *m_pSAMloci;		// pts to memory dynamically allocated as needed to hold all accepted alignment loci

	int
	PrefixFasta(char* pszPrefix,	// descriptor prefix
			char* pszInFile,		// input fasta or SAM file
			char* pszOutFile);		// output to this file

	int
	FilterSAM(char* pszPrefix,			// target prefix used for filtering from
			char* pszInFile,		// input SAM file with matching alignments
			char* pszOutFile);		// output to this SAM file

	int	
		GenBinnedWiggle(eModePG PMode,			// processing mode
			char *pszName,	// wiggle track name
			uint32_t BinSizeKbp,	// Wiggle score is number of alignments over this sized bins
			char* pszInFile,		// alignments are in this SAM/BAM file 
			char* pszOutFile);		// write out Wiggle to this file

	int		// bin score assigned 
		GenBinnedCoverage(int TargID,	// coverage is on this targeted chrom/seq
				uint32_t BinStart,		// coverage starts at this loci inclusive
				uint32_t BinEnd,		// coverage ends at this loci inclusive
				uint32_t BinCnts);		// total number of coverage counts

	int		// returned sequence identifier, < 1 if unable to accept this sequence name
		AddTargSeqName(char* pszSeqName);	// associate unique identifier with this sequence name

	char*							// returned sequence name
		LocateTargSeqName(int SeqID);	// identifier returned by call to AddTargSeqName

	// SortSAMTargLoci
// Sort m_pSAMloci by ascending Targ.Loci
	static int SortSAMTargLoci(const void* arg1, const void* arg2);

public:
	CPangenome();
	~CPangenome();
	void Reset(void);	// resets class instance state back to that immediately following instantiation

	int Process(eModePG PMode,			// processing mode
			char* pszPrefix,		// descriptor prefix
			int BinSizeKbp,	// Wiggle score is number of alignments over this sized bins
			char* pszInFile,		// input fasta or SAM file
			char* pszOutFile);		// output to this file
};
