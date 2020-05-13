/*
This toolkit is a source base clone of 'BioKanga' release 4.4.2 (https://github.com/csiro-crop-informatics/biokanga) and contains
significant source code changes enabling new functionality and resulting process parameterisation changes. These changes have resulted in
incompatibility with 'BioKanga'.

Because of the potential for confusion by users unaware of functionality and process parameterisation changes then the modified source base
and resultant compiled executables have been renamed to 'kit4b' - K-mer Informed Toolkit for Bioinformatics.
The renaming will force users of the 'BioKanga' toolkit to examine scripting which is dependent on existing 'BioKanga'
parameterisations so as to make appropriate changes if wishing to utilise 'kit4b' parameterisations and functionality.

'kit4b' is being released under the Opensource Software License Agreement (GPLv3)
'kit4b' is Copyright (c) 2019, 2020
Please contact Dr Stuart Stephen < stuartjs@g3web.com > if you have any questions regarding 'kit4b'.
 */

#pragma once


const int cMaxArtefSeqs = 20;			// allow at most this number of artefact sequences
const int cMaxArtefSeqLen = 40;			// artefact sequences can be at most this length


const int cSRDfltNumReads = 10000000;		// default number of reads
const int cSRMinNumReads =  5;			// minimum number of reads
const int cSRMaxNumReads =  500000000;	// maximum number of reads
const int cSRMinReadLen  =  20;			// minimum read length
const int cSRDfltReadLen= 100;			// default read length
const int cSRMaxReadLen = 1000;			// max simulated read length
const int cSRMaxCIGARLen = 200;			// limit any CIGAR string to be at most this length

const int cSRMinChromLen = 1000;		// reads simulated from chromosomes of at least this length

const int cUpdnstream = 2000;			// when regional filtering this is the length of the up/down stream region

const int cSRMinCutLen   =  cSRMinReadLen;	// minimum cut length
const int cSRMaxCutLen   =  cSRMaxReadLen;	// maximum cut length

const int cMinPEFragLen     = 75;		// paired end minimum fragment length
const int cMaxPEFragLen     = 100000;	// paired end maximum fragment length

const int cMaxDistClusterLen = 300;		// max clustered read bin size

const int cMaxProfileBatchSize = 500000;  // max number of end profiled reads to simulate per thread before checkpointing to disk
const int cMaxBatchSize = 5000000;		// max number of defaulted Hamming and non-profiled reads to simulate per thread before checkpointing to disk

#define NUCONLYMSK (~cRptMskFlg & 0x0f)	// hiorder bits used as attributes - bit 5 used to flag subsequence already selected
#define SSSELECTED 0x010				// used as an attribute to flag subsequence starting this loci already selected

const int cSRMaxAllocBuffChunk = 0x07fffff;			// buffer fasta
const int cSRAllocNumChroms = 4096;					// allocate for chrom seqs in this many increments

const int cSRAllocdObsErrProfMemChunk = 0x07ffffff;		// allocate for observed error profiles in this sized memory chunks 

// processing modes
typedef enum TAG_eSRPMode {
	eSRPMStandard,					// default - standard random start and end
	eSRPMProfRand,					// profiled start with random end sites
	eSRPMRandProf,					// random start with profiled end sites
	eSRPMProfProf,					// profiled start with profiled end sites
	eSRPMplaceholder				// used to set the enumeration range
	} etSRPMode;

typedef enum TAG_eSRBEDRegion {
	eSRMEGRAny = 0,		// process any region
	eSRMEGRCDS,			// part of feature overlaps CDS
	eSRMEGR5UTR,			// part of feature overlaps 5'UTR
	eSRMEGR3UTR,			// part of feature overlaps 3'UTR
	eSRMEGRIntrons,		// part of feature overlaps Intron
	eSRMEGRUpstream,		// part of feature overlaps 5'upstream regulatory
	eSRMEGRDnstream,		// part of feature overlaps 3'downstream regulatory
	eSRMEGRIntergenic		// part of feature overlaps intergenic
} etSRBEDRegion;


// output format
typedef enum TAG_eSRFMode {
	eSRFMFasta,				// multifasta, sequences wrapped if longer than 79bp
	eSRFMNWFasta,			// multifasta, sequences non-wrapped even if longer than 79bp
	eSRFMplaceholder		// used to set the enumeration range
	} etSRFMode;

// simulated error rate mode selection
typedef enum TAG_eSRSEMode {
	eSRSEPnone,				// no simulated errors
	eSRSEPfixerrs,			// simulate fixed number of errors in each read
	eSRSEPstatic,				// use internal static profile
	eSRSEPdyn,				// dynamic according to '-z<rate>'
	eSRSEPplaceholder			// used to set the enumeration range
	} etSRSEMode;

#pragma pack(1)

typedef struct TAG_sPackedCIGARS {
	uint32_t ID;						// observed CIGAR identifier
	uint16_t Size;						// actual size of this tsPackedCIGARS instance in bytes 
	uint16_t ReadLen;					// read length
	uint8_t FlgPE1Strand:1;				// PE1 strand - 0 if watson, 1 if crick
	uint8_t FlgPE2Strand : 1;			// PE2 strand - 0 if watson, 1 if crick
	uint8_t PE1NumCIGAROps;				// number of PE1 CIGAR operators
	uint8_t PE2NumCIGAROps;				// number of PE2 CIGAR operators
	uint32_t PEInsertSize;				// if PE then insert size
	uint32_t CIGAROps[1];	// place holder for read CIGAR operations; PE1 first followed by PE2
} tsPackedCIGARS;

typedef struct TAG_sSRChromSeq {
	int ChromSeqID;			// m_pChromSeqs[ChromSeqID-1] of this tsChromSeq instance
	int ChromID;
	char Strand;			// '*' if chrom seq is for genome assembly, '+' or '-' if seq is for a feature or gene
	int RelDensity;			// relative (0 to 1000), to other reads, compared with rand(0,1000) and if less then no read generated in current iteration
	int ScaledStart;		// scaled start of this chrom
	int ScaledLen;			// chrom length scaled such that the sum of all chrom scaled lengths is less than INT_MAX
	char szChromName[cMaxDatasetSpeciesChrom];
	int Len;				// actual length
	UINT32 SeqOfs;			// offset at which this chroms sequence starts
} tsSRChromSeq;


#pragma pack(4)
typedef struct TAG_sSRSimRead {
	int ChromSeqID;			// simulated read is on this m_pChromSeqs[]
	int Status;				// 0 if yet to be reported on, 1 if reported, 2 if not reported because dup read seq, 3 if not reported because from same loci
	int FlgPE2:1;			// if paired reads generation then set if this is the 3' end
	struct TAG_sSRSimRead *pPartner; // if paired reads generation then points to partner of this read
	int ChromID;			// chromosome identifier
	int Strand;				// on this strand (0 if '+', 1 if '-')
	int StartLoci;			// starts at this loci
	int EndLoci;			// ends at this loci
	int Len;				// and was originally of this length
	int InDelLen;			// and has had a deletion (<0) or an insertion (>0) of this length
	int Lenx;				// resulting after the InDel in the simulated read being of this length
	etSeqBase *pSeq;		// pts to start of this simulated reads sequence
	} tsSRSimRead;

// Induced errors will be heavily 3' biased simulating Illumina read substitution profiles
typedef struct TAG_sSRInduceErrProf {
	double Proportion;			// proportion of reads (0 - 100.0) with
	int NumSubs;				// this number of sequencer errors
} tsSRInduceErrProf;

typedef struct TAG_sReadRprt {
	etSeqBase FwdSeq[cMaxReadLen + 1];
	etSeqBase RevSeq[cMaxReadLen + 1];
	char szColorspace[cMaxReadLen + 2];
	char szLineBuff[(cMaxReadLen + 2) * 2];
	int LineLen;
	int NumCols;
} sReadRprt;

typedef struct TAG_sSRWorkerPars {
	int ThreadIdx;					// index of this thread (1..m_NumThreads)
	void *pThis;					// will be initialised to pt to CSimReads instance
#ifdef _WIN32
	HANDLE threadHandle;	// handle as returned by _beginthreadex()
	unsigned int threadID;	// identifier as set by _beginthreadex()
#else
	int threadRslt;			// result as returned by pthread_create ()
	pthread_t threadID;		// identifier as set by pthread_create ()
#endif
    int Rslt;				// processing result code
	int RandSeed;			// random generator seed for this worker
	etSRPMode PMode;			// processing mode
	int Region;				//  Process regions 0:ALL,1:CDS,2:5'UTR,3:3'UTR,4:Introns,5:5'US,6:3'DS,7:Intergenic (default = ALL)
	int UpDnRegLen;			// if processing regions then up/down stream regulatory length
	int NumReqReads;		// number of reads to be generated by this worker thread
	int NumGenReads;		// number of actually generated reads
	bool bDedupe;			// if true then slough any reads with a Hamming of 0
	int ReadLen;			// read length
	int CutMin;				// min cut length
	int CutMax;				// max cut length
	char Strand;			// generate for this strand '+' or '-' or for both '*'
	bool bPEgen;			// true if paired ends are to be generated
	bool bMaxIters;			// set true by thread if exhusted attempts to find chrom from which reads can be simulated
	int PEmin;				// PE minimum fragment
	int PEmax;				// PE maximum fragment
	tsSRSimRead *pReads;		// to hold all reads generated by this worker thread - will have been preallocated to hold NumReqReads
} tsSRWorkerPars;
#pragma pack()

class CSimReads
{
	etSRPMode m_PMode;				// processing mode
	etSRFMode m_FMode;				// output format
	int m_hOutFile;					// output results file handle
	char *m_pszOutFile;				// output file name

	int m_hOutPEFile;					// output results paired end file handle
	char *m_pszOutPEFile;				// output paired end file name

	CCSVFile *m_pProfCSV;			// used to hold profile site preferences whilst loading into m_pProfSel
	double *m_pProfSel;				// allocated array of profile site selection preferences (0.0..1.0) indexed by sequence octamers

	int m_PropRandReads;			// what proportion ( * 1 million) are to be generated as completely random reads
	int m_DistCluster;				// cluster generated reads into windows of this median length, 0 if no clustering
	int m_InDelSize;				// simulated InDel size range
	double m_InDelRate;				// simulated InDel rate per read


	double m_Artef5Rate;			// rate (0..1) at which to insert 5' artefact sequences
	int m_NumArtef5Seqs;			// number of user specified 5' artefact sequences
	char **m_ppszArtef5Seqs;		// 5' artefact sequences
	double m_Artef3Rate;			// rate (0..1) at which to insert 3' artefact sequences
	int m_NumArtef3Seqs;			// number of user specified 3' artefact sequences
	char **m_ppszArtef3Seqs;		// 5' artefact sequences

	int m_Artef5SeqLens[cMaxArtefSeqs];				// sequence lengths for each 5' artefact
	etSeqBase m_Artef5Seqs[cMaxArtefSeqs][cMaxArtefSeqLen+1];	// each 5' artefact sequence
	int m_Artef3SeqLens[cMaxArtefSeqs];				// sequence lengths for each 3' artefact
	etSeqBase m_Artef3Seqs[cMaxArtefSeqs][cMaxArtefSeqLen+1];	// each 3' artefact sequence

	CBioSeqFile *m_pBioSeqFile;		// genome assembly
	CBEDfile *m_pBEDFile;			// optional if simulating transcript reads

	char m_szSpecies[cMaxDatasetSpeciesChrom+1];		// species (title) as retrieved from bioseqfile
	int m_NumChromSeqs;				// number of chromosomes loaded
	int m_AllocdChromSeqs;			// number allocated
	size_t m_AllocdChromSize;		// allocd mem size for m_pChromSeqs
	tsSRChromSeq *m_pChromSeqs;		// pts to chromseqs array

	UINT8 *m_pGenomeSeq;			// allocated to hold concatenated (separated by eBaseEOSs) chromosome sequences
	size_t m_AllocdGenomeMem;		// allocd mem size for m_pGenomeSeq
	INT64 m_GenomeLen;				// total genome length including concatenators
	INT32 m_GenomeScaledLen;		// sum of all chrom scaled lengths, will always be less than INT_MAX

	tsSRSimRead *m_pSimReads;			// allocated to hold simulated reads - NOTE: mmap/malloc used because of GNU malloc limitations
	INT64 m_AllocdMemReads;			// actual memory allocation size used when allocating m_pSimReads
	int m_NumReadsAllocd;			// this many reads have been allocated for in m_pSimReads
	int m_TotReqReads;				// number of reads required to be simulated - will be 2x user requested number if simulating paired end reads
	int m_CurNumGenReads;			// current number of generated reads - updated every N reads generated by worker threads

	int m_MaxFastaLineLen;			// wrap sequences in multifasta output files if line is longer than this many bases
	int SimInDels(tsSRSimRead *pSimRead,int *pReadLen,etSeqBase *pRead);
	int	SimArtefacts(bool b3ArtefSeq,		// if false then 5' artefact else if true then 3' artefact
			int ReadLen,etSeqBase *pRead);

	int // number of substitutions inplace induced into this read
		SimSeqErrors(int SeqLen,etSeqBase *pRead);

	int // number of substitutions inplace induced into this read
		SimSeqRand(int SeqLen,etSeqBase *pRead);


#ifdef _WIN32
	HANDLE m_hMtxDedupe;
#else
	pthread_mutex_t m_hMtxDedupe;
#endif

	static int SortSimReads(const void *arg1, const void *arg2);

public:
	CSimReads();
	~CSimReads();

	void Reset(bool bSync);

	static int SortSimLoci(const void *arg1, const void *arg2);
	static int	TrimSeqChrs(char *pszTxt);	// trims quote marks, space/tabs and validates sequence as only containing acgtu
	static char *Region2Txt(etSRBEDRegion Region);

	int
		GenSimReads(etSRPMode PMode,		// processing mode
				etSRSEMode SEMode,	// induced sequencer error rate mode
				bool bPEgen,		// true if paired ends are to be generated
				int PEmin,			// PE minimum fragment length
				int PEmax,			// PE maximum fragment length
				double PropRandReads, // generate completely random reads at this rate
				int DistCluster,	// cluster generated reads into windows of this median length, 0 if no clustering
				double SeqErrRate,	// dynamic sequencing errors are to be induced into generated sequences at this rate
				bool bSeqErrProfile,// true if to generate composite sequencer errors with uniform profile (default is Illumina 3' skewed)
				int SNPrate,		// simulate SNPs at this rate per million bases
				int InDelSize,		// simulated InDel size range
				double InDelRate,	// simulated InDel rate per read
				etSRFMode FMode,		// output format
				int NumThreads,		// number of worker threads to use
				char Strand,		// generate for this strand '+' or '-' or for both '*'
				int NumReads,		// number of reads required (will be 2x this number if generating paired ends)
				int ReadLen,		// read lengths
				double Artef5Rate,			// rate (0..1) at which to insert 5' artefact sequences
				int NumArtef5Seqs,			// number of user specified 5' artefact sequences
				char *pszArtef5Seqs[], // 5' artefact sequences
				double Artef3Rate,			// rate (0..1) at which to insert 3' artefact sequences
				int NumArtef3Seqs,			// number of user specified 3' artefact sequences
				char *pszArtef3Seqs[], // 5' artefact sequences
				int CutMin,			// min cut length
				int CutMax,			// max cut length
				bool bDedupe,		// true if unique read sequences only to be generated
				int Region,			// Process regions 0:ALL,1:CDS,2:5'UTR,3:3'UTR,4:Introns,5:5'US,6:3'DS,7:Intergenic (default = ALL)
				int UpDnRegLen,		// if processing regions then up/down stream regulatory length
				char *pszFeatFile,	// optionally generate transcriptome reads from features or genes in this BED file
				char *pszInFile,	// input from this raw multifasta or bioseq assembly
				char *pszProfFile,	// input from this profile site preferences file
				char *pszOutPEFile, // output partner paired end simulated reads to this file
				char *pszOutFile,	// output simulated reads to this file
				char *pszOutSNPs);   // output simulated SNP loci to this file

	int LoadFasta(size_t *pTotLen,int MinChromLen,char *pszFastaFile);
	int LoadBioseq(size_t *pTotLen,int MinChromLen,char *pszBioSeqFile);
	int LoadGenome(int MinChromLen,			// warn if loaded chromosomes are less than this length; may be too short to sample reads from
			char *pszBioSeqFile);

	int
	LoadTranscriptome(char *pszBioSeqFile,			// genome assembly
				char *pszFeatFile);				// BED file containing features or genes

	int
	SimulateSNPs(char *pszOutSNPs,   // output simulated SNP loci to this file
			int SNPrate);		 // required SNPs per million bases

	
	int	CmpSeqs(int Len, etSeqBase* pSeq1, etSeqBase* pSeq2); // Returns the number of differences between two sequences of length Len
	bool IsDupSeq(int Len, etSeqBase* pSeq1, etSeqBase* pSeq2); // sequences the same?
	bool IsDupSeqStrands(int Len, etSeqBase* pReadSeq1, etSeqBase* pReadSeq2); // both the forward and reverse complements are compared for equivilence
	void ShowDupSeqs(int ReadLen, tsSRSimRead* pRead1, tsSRSimRead* pRead2); // display duplicated sequences

	int
		ReportReads(bool bPEgen,			// true if paired end simulated reads being simulated
		    int Region,			// Process regions 0:ALL,1:CDS,2:5'UTR,3:3'UTR,4:Introns,5:5'US,6:3'DS,7:Intergenic (default = ALL)
			int NomReadLen,			// reads nominal length
			int NumPrevReported,	// number of reads previously reported
			int MaxReads,			// maximum number of reads remaining to be reported
			bool bDedupe,			// true if reads are to be deduped
			int AvailReads);			// number of reads to be reported on from m_pSimReads[]

	int
		LocateRevCplSeq(int Len,etSeqBase *pProbe,int NumSortedReads);

	int			// index of exactly matching probe or -1 if no match
		LocateFirstExact(etSeqBase *pProbe,				// pts to probe sequence
				 int ProbeLen,					// probe length to exactly match over
				  int IdxLo,					// low index in m_pSimReads
				  int IdxHi);					// high index in m_pSimReads

// generate '+' strand index from K-mer of length SeqLen starting at pSeq
	int GenPSeqIdx(int SeqLen,etSeqBase *pSeq);	
// generate '-' strand index from K-mer of length SeqLen starting at pSeq
int GenMSeqIdx(int SeqLen,etSeqBase *pSeq);

int
InitProfSitePrefs(char *pszInProfFile);	// read from this profile site selectivity file (for MNase, generated by MNaseSitePred process), if NULL then static profiling

	int ThreadSimReads(void * pThreadPars);
};

