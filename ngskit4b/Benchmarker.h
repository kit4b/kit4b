#pragma once

// base weightings will be silently clamped to be in the range of cBMMinPts..cBMMaxPts
// defaults were chosen to heavily favor correctly aligned bases, and that it is preferable to have unaligned bases rather than misaligned bases
const int cBMMinBaseWtg = 0;					// base weighting minimum
const int cBMMaxBaseWtg = 100;					// base weighting max
const int cBMDfltAlignedBaseWtg = 50;			// default points to apply for each base aligned to it's ground truth loci
const int cBMDfltSilentTrimAlignBaseWtg = 25;	// default weighting for if read was silently trimmed and falls within ground truth loci range then credit this weighting - Ugh, Ugh and Ugh again 
const int cBMDfltUnalignedBaseWtg = 20;			// default weighting to apply for each unaligned base
const int cBMDfltMisalignedBaseWtg = 0;			// default weighting to apply for each base aligned but not to ground truth loci

// silently clamping number of reads, SE or PE pairs, to be within following ranges
const int cBMMaxReads = 100000000;		// max number of reads or read pairs to process for MAGIC simulations
const int cBMDfltLimitReads = 5000000;	// restrict number of raw reads which are to be aligned
const int cBMDfltReads = 2000000;		// default number of reads or read pairs to process for MAGIC simulations
const int cBMMinReads = 100;			// min number of reads or read pairs to process for MAGIC simulations
const int cBMMinReadLen = 50;			// only processing for read sequences in range cBMMinReadLen
const int cBMMaxReadLen = 1000;			// upto cBMMaxReadLen inclusive

const int cBMMaxFragSize = 200000;				// can handle PE fragment sizes upto this inclusive size

const int cBMCIGARAllocBuffSize = 10000000; // allocation size for buffering CIGAR CSV file read/writes and simulated reads

typedef enum _TAG_eBMProcModes {
	eBMLimitReads = 0,	// limit number of raw reads which are to be aligned
	eBMGenCIGARs,		// generate observed CIGARs from alignments
	eBMSimReads,		// simulate reads using observed CIGARs
	eBMScore			// score base alignments using ground truth CIGARs and loci
} eBMProcMode;

#pragma pack(1)

typedef struct TAG_sBMObsCIGAR {
	uint16_t Size;			// total size of this structure instance
	uint8_t NameLen;		// read name length
	uint8_t ObsCIGARlen;	// observed CIGAR length
	uint16_t ErrProfileLen;		// error profile for read
	uint16_t Flags;			// same as flags in SAM format specification
	int32_t ReadLen;		// read length
	int32_t InsertSize;		// if PE pair then insert size
	uint8_t NameCIGARErrProfile[1];	// read name concatenated with observed CIGAR from read and then the read error profile
} tsBMObsCIGAR;

typedef struct TAG_sBMChromSeq {
	int ChromID;			// chromosome identifier in loaded suffix array genome
	int Len;				// actual length
	int ScaledLen;			// chrom length scaled such that the sum of all chrom scaled lengths is less than INT_MAX
	int ScaledStart;		// scaled start of this chrom
} tsBMChromSeq;

typedef struct TAG_sBMPackedCIGARs {
	uint32_t ID;						// observed CIGAR identifier
	uint16_t Size;						// actual size of this tsBMPackedCIGARs instance in bytes 
	uint16_t ReadLen;					// read length
	uint8_t FlgPE1Strand : 1;			// PE1 strand - 0 if watson, 1 if crick
	uint8_t FlgPE2Strand : 1;			// PE2 strand - 0 if watson, 1 if crick
	uint16_t PE1NumCIGAROps;			// number of PE1 CIGAR operators
	uint16_t PE1NumErrProfOps;			// number of PE1 error profile operators
	uint16_t PE2NumCIGAROps;			// number of PE2 CIGAR operators
	uint16_t PE2NumErrProfOps;			// number of PE1 error profile operators
	uint32_t PEInsertSize;				// if PE then insert size
	uint32_t CIGAROpsErrProfOps[1];		// place holder for read CIGAR + error profile operations; PE1 first followed by PE2
} tsBMPackedCIGARs;

typedef struct TAG_sBMGroundTruth {
	uint32_t ID;						// uniquely identifies this ground truth instance
	uint32_t ReadID;					// ground truth identifier expected to be unique for each read SE or read pair
	uint16_t Size;						// actual size of this tsBMGroundTruth instance in bytes
	int8_t NameLen;						// read name length, excludes terminating '\0'
	int8_t ChromNameLen;				// chromosome name length, excludes terminating '\0'
	int8_t CIGARLen;					// unpacked CIGAR string length, excludes terminating '\0'
	uint32_t StartLoci;					// simulated read started at this loci on ChromName
	uint8_t FlgStrand : 1;				// strand - 0 if watson, 1 if crick
	uint8_t FlgPE2 : 1;					// 0 if SE or PE1 of a PE, 1 if PE2
	uint8_t FlgAligned : 1;				// set if alignment read with same name and matching SE or PE was discovered - picks up multialignments to same ground truth
	uint8_t FlgStrandErr : 1;			// aligned to incorrect strand
	uint8_t FlgRefChromErr : 1;			// aligned to incorrect ref chrom
	uint8_t FlgPE2Err : 1;				// aligned pair end is incorrect - alignment PE2 does not match ground truth PE2
	uint32_t ReadLen;					// ground truth is for a simulated read of this length
	uint32_t PotentialBasesAligning;	// potentially this ground truth read has this many bases aligning to the target
	uint8_t NameChromCIGAR[1];			// read name concatenated chromosome name concatenated with CIGAR
} tsBMGroundTruth;


#pragma pack()

class CBenchmark {
	bool m_bPrimaryOnly;			// if true then only score primary read alignments otherwise score all including secondary
	bool m_bPEReads;				// if true then PE pair only processing otherwise treating all reads as if SE reads
	int m_UnalignedBasePts;			// points to apply for each unaligned base
	int m_AlignedBasePts;			// points to apply for each base aligned to it's ground truth loci
	int m_SilentTrimAlignBasePts;	// points to apply for each base aligned in a silently trimmed read within ground truth loci range
	int m_MisalignedBasePts;		// points to apply for each base aligned but not to ground truth loci

	int64_t m_TotNumPotentialAlignBases;	// number of match bases in actual alignments which potentially could have been aligned to ground truth
	int64_t m_NumBasesLociCorrect;			// total number of bases aligned correctly to ground truth loci
	int64_t m_NumBasesLociIncorrect;		// total number of bases aligned incorrectly to ground truth loci
	int64_t m_NumBasesLociUnclaimed;		// total number of ground truth bases which were not aligned
	int64_t m_NumSilentTrimBaseMatches;		// total number of bases accepted as aligned even though in reads which have been silently trimmed - neither soft or hard trimmed

	uint32_t m_ReadOverlapHistogram[101];	// histogram of read alignments by percentile proportion of read bases in reads overlapping with ground truth read bases

	int m_MaxNumReads;				// maximum number of alignment CIGARs to process or number of simulated reads or read pairs 

	char m_szObsCIGARsFile[_MAX_PATH];	// observed CIGARs are in this file
	char m_szGroundTruthFile[_MAX_PATH];	// simulated reads ground truth (CIGARs and loci for each simulated read)) are in this file
	char m_szRefGenomeFile[_MAX_PATH];// reads are against this target genome file
	char m_szAlignmentsFile[_MAX_PATH];	// input file containing aligned reads (SAM or BAM)
	char m_szSEReads[_MAX_PATH];		// simulated reads are output to this file (SE or PE1 if PE)
	char m_szPE2Reads[_MAX_PATH];		// simulated PE2 reads are output to this file if simulating PE reads
	char m_szResultsFile[_MAX_PATH];	// benchmarking m2 results appended to this CSV file
	char m_szExperimentDescr[cMaxDatasetSpeciesChrom + 1];	// experiment descriptor by which benchmarking results can be identified in szResultsFile

	INT64 m_GenomeLen;				// total genome length including concatenators
	INT32 m_GenomeScaledLen;		// sum of all chrom scaled lengths, will always be less than INT_MAX
	int m_NumChromSeqs;				// number of chromosomes loaded
	tsBMChromSeq *m_pChromSeqs;		// pts to chromseqs array

	int m_hSEReads;					// simulated reads are output to this file handle (SE or PE1 if PE)
	int m_hPE2Reads;				// simulated reads are output to this file handle (PE2)

	int m_AllocOutBuffSize;			// allocated simulated read buffers of this size
	int m_OutSEBuffIdx;				// currently using this number of chars in pszSESimReadBuff
	char *m_pszSESimReadBuff;			// allocated to buffer SE or PE1 simulated reads prior to disk write
	int m_OutPE2BuffIdx;			// currently using this number of chars in pszSESimReadBuff
	char* m_pszPE2SimReadBuff;		// allocated to buffer PE2 simulated reads prior to disk write

	int m_hObsSIGARs;				// observed CIGARs CSV file handle
	int m_ObsCIGARBuffLen;			// current ObsCIGARBuff buffering
	int m_AllocdObsCIGARBuffSize;	// m_pObsCIGARBuff allocation size	
	char * m_pObsCIGARBuff;			// allocated for buffering CIGARs

	char m_szTargSpecies[cMaxDatasetSpeciesChrom+1];	// reference genome species
	CSfxArray *m_pGenome;				// reference genome
	CSAMfile *m_pAlignments;			// alignments against the reference genome
	CCSVFile *m_pObsCIGARProfFile;		// observed CIGAR alignment profiles

	CFasta *m_pSESimReads;				// when scoring used to process simulated SE or PE1 reads for ground truth
	CFasta* m_pPE2SimReads;				// when scoring used to process simulated SE or PE1 reads for ground truth

	uint32_t m_NumObsErrProfs;		// loaded this many observed error profiles
	size_t m_UsedObsErrProfMem;		// loaded observed error profiles are using this sized memory
	size_t m_AllocdObsErrProfMem;	// allocated this size memory to hold observed error profiles
	uint8_t* m_pObsErrProfiles;		// pts to memory allocated

	uint32_t m_UnscoredReads;		// number of aligned simulation reads which were unable to be scored - unaligned, missing CIGAR etc
	uint32_t m_ScoredReads;			// number of aligned simulation which were scored

	uint32_t m_NumGroundTruthReads;		// loaded this many ground truths loaded from simulated reads
	INT64 m_TotGroundTruthReadBases;	// total number of sequence bases in loaded ground truth reads
	size_t m_UsedGroundTruthsMem;	// loaded ground truths are using this sized memory
	size_t m_AllocdGroundTruthsMem;	// allocated this size memory to hold observed error profiles
	uint8_t *m_pGroundTruths;		// allocated to hold ground truths loaded from simulated reads
	tsBMGroundTruth **m_ppGroundTruthIdx;	// allocated to hold ground truth index sorted by read names ascending


	void Reset(void);					// reset instance state back to that immediately following instanciation
	int	LoadGenome(char* pszRefGenomeFile);		// read alignments were against this targeted genome file
	int OpenAlignments(char* pszAlignmentsFile);// input file containing aligned reads (SAM or BAM)

	int		// returned number of potential base matches in CIGAR
		PotentialMatchBases(int ReadLen,			// read sequence length
						char *pszCIGAR,				// '\0' terminated CIGAR with alignment profile
						bool bTreatSoftHardAsMatches=false);  // if true then treat soft and hard clipped as if matches (used when scoring aligned reads)

	int    // returns number of claimed base matches which were ground truth bases
		ActualMatchBases(tsBAMalign* pAlignment,			// claimed alignment
			tsBMGroundTruth* pGroundTruth);	// ground truth for this alignment
			

	int
		AdjustAlignment(etCIGAROpType OpType,	// alignment operator type
			int OpCnt,							// number of times operator to be applied
			int* pReadOfs,						// current read offset
			int* pAlignLoci);					// read offset corresponds to this alignment loci

	int										    // number of chars consumed in CIGAR, 0 if none consumed, -1 if errors
		ParseExpCIGAR(char* pszCIGAR,			// parse starting from this CIGAR char
			int* pNumOps,						// returned number of operations required
			etCIGAROpType* pTypeOp);			// type of operation

	int					// returns num of chars in decoded packed CIGAR string excluding terminating '\0'
		DecodeCIGAR(int NumPackedOps,			// number of packed CIGAR ops in
					uint32_t *pPackedOPs,		// ptr to packed CIGAR ops
					int MaxCIGARChrs,			// pszCIGAR allocated to hold at most this many chars including terminating '\0'
					char *pszCIGAR,			// decode into this buffer
					bool bTreatSoftHardAsMatches=false);  // if true then treat soft and hard clipped as if matches (used when generating simulated reads)

	int									// returns index of next ClaimOpIdx, 0 if no more claim Ops, -1 if errors 
		GetClaimOps(tsBAMalign* pAlignment,	// claim CIGAR in this alignment
			int ClaimOpIdx,			// get ClaimNumOps, ClaimOpType for this CIGAR operation, starts from 1
			int* pNumOps,			// returned number of operations required
			etCIGAROpType* pTypeOp);		// type of operation

	int							  // number of accepted CIGAR operations, 0 if unable to accept 
		ParseObsCIGAR(char* pszCIGAR, // CIGAR operations
			int ReadLen, // length ('=','X','M','D') must be equal to ReadLen
			int MaxCIGAROps,	// pPackedOps can hold at mos this many Ops
			uint32_t* pPackedOps);	// pack CIGAR ops into this array, can be NULL if just validating the CIGAR

	int		// for given CIGARs returns total number of bases that these CIGARs consume on aligned to reference sequence
		RefSeqConsumedLen(int NumPackedOps,	  // there are this many packed CIGAR ops in pPackedOps
			uint32_t* pPackedOps, // packed CIGAR ops into this array
			bool bTreatSoftHardAsMatches=false);	// if true then treat soft and hard clipped as if matches (used when generating simulated reads)

	int				// returns number of alignment inplace coalesced CIGAR ops
		CoalesceCIGARs(tsBAMalign *pBAM);	// in: alignment to coalesce out: updated alignment with coalesced CIGAR ops

	int	InitObsErrProfile(char* pszInProfFile);	// read from this observed alignment error profiles file (as generated from actual SAM aligments)

	int	LoadGroundTruths(char* pszSESimReads,		// load ground truths for SE or PE1 if PE simulated reads
			char* pszPE2SimReads,					// load ground truths for PE2 if PE simulated reads
			bool bTreatSoftHardAsMatches=false);	// if true then treat soft and hard clipped as if matches

	tsBMGroundTruth*		// returned ground truth which matches
		LocateGroundTruth(tsBAMalign* pSAMalign,	// this alignment by name
							bool b2ndInPair);		// if true (assumes PE ground truths) then locate PE2

	int			// choosen loci on length proportional randomly selected chromosome, -1 if unable to choose a chromosome meeting critera
		ChooseRandStartLoci(int FragSize,		// ensure can generate this sized fragment starting at chrom/loci
			int* pSelChromID);	// returned chrom identifier

	int										// eBSFSuccess if sequence written to file otherwise the error code
		WriteFastaFile(bool bPEReads,	// true if simulations are for PE reads
			bool bPE2,			// true if writing out a PE2 simulated reads false if SE or PE1
			int SeqID,					// primary sequence identifier
			tsBMPackedCIGARs* pCurCIGAR,	// CIGAR error profile 
			int ChromID,				// simulated read is on this chromosome 
			int StartLoci,				// starting at this loci (0 based)
			etSeqBase* pSeq);			// sequence

	static int SortReadPairs(const void* arg1, const void* arg2);
	static int SortGroundTruths(const void* arg1, const void* arg2);

public:
	CBenchmark();						// instance constructor
	~CBenchmark();						// instance destructor

	int
		GenLimitReads(bool bPEReads,	// true if PE pair processing otherwise SE reads
			int MaxNumReads,			// restrict reads to be at most this many SE reads or PE read pairs
			char* pszInSEReads,			// SE or PE1 reads are input from this file
			char* pszInPE2Reads,		// PE2 reads are input from this file
			char* pszOutSEReads,		// SE or PE1 reads are output to this file
			char* pszOutPE2Reads);		// PE2 reads are output to this file

	int
		GenObsCIGARs(bool bPEReads,		// true if PE pair processing otherwise SE reads
				int MaxNumReads,		// maximum number of aligned reads or read pair alignments to process
				char* pszObsCIGARsFile,	// write observed CIGARs to this file
				char* pszRefGenomeFile,	// alignments were against this target genome file
				char* pszAlignmentsFile);	// input file containing alignments of simulated reads (SAM or BAM)

	int
		SimReads(bool bPEReads,		// true if PE pair processing otherwise SE reads
			int MaxNumReads,			// maximum number of simulated reads or read pairs
			char* pszCIGARsFile,		// read in observed CIGARs from this file
			char* pszGroundTruthFile,	// simulated reads ground truth (MAGIC and loci for each simulated read)) are writen to this file
			char* pszSEReads,			// simulated reads are output to this file (SE or PE1 if PE)
			char* pszPE2Reads);			// simulated PE2 reads are output to this file if simulating PE reads


	int
		Score(bool bScorePrimaryOnly,		// score only primary alignments otherwise score including secondary
			bool bScoreMatedPE,			// if true then both mates of a PE must have been aligned for alignment to be scored
			bool bPEReads,				// true if PE pair processing
			int UnalignedBasePts,		// points to apply for each unaligned base
			int AlignedBasePts,			// points to apply for each base aligned to it's ground truth loci
			int SilentTrimAlignBasePts, // points to apply for each base aligned in a silently trimmed read within ground truth loci range
			int MisalignedBasePts,		// points to apply for each base aligned but not to ground truth loci
			char* pszResultsFile,		// benchmarking m2 results appended to this CSV file
			char* pszExperimentDescr,	// experiment descriptor by which benchmarking results can be identified in szResultsFile
			char* pszControlAligner,	// control aligner generating error profile from which simulated reads were generated 
			char* pszScoredAligner,		// aligner aligning simulated reads and which was scored
			char* pszSEReads,			// input simulated reads which contain ground truths from this file for SE or PE1 if PE
			char* pszPE2Reads,			// input simulated reads which contain ground truths from this file for PE2 if PE
			char* pszAlignmentsFile);   // input file containing alignments of simulated reads (SAM or BAM)

};