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

Original 'BioKanga' copyright notice has been retained and immediately follows this notice..
*/
/*
 * CSIRO Open Source Software License Agreement (GPLv3)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENSE for the complete license information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */

#pragma once

const uint32_t cMinRSSeqLen = 16;						// raw sequences must be at least this length after any end timming to be acceptable for processing
const uint32_t cMaxRSSeqLen = cMaxReadLen;			// can accept read sequences upto this length

const int cDfltPhredScore = 30;						// if no quality scores available in input reads then use this as the default quality score

const int cMinBasePhredScore = 10;					// if processing for minimum mean Phred scores then this is the minimum which can be specified, and any individual base must have at least this Phred

const uint32_t cMaxContamSeqs = 250;					// max of this many contaminant sequences allowed
const uint32_t cReAllocContamSeqs = 64;				// initially allocate then realloc as needed for this many contaminate sequences

const int cMaxKMerLen = 12;							// can process maximal sized K-mers of this length

const int cDfltContamSubRate = 1;					// default allowed contamimamt substitution rate (bases per 25bp of contaminant overlap, 1st 15bp of overlap no subs allowed)
const int cDfltMinContamLen = 5;					// default is to accept contaminant overlaps if overlaps of at least this many bases 

const int cMinDupSeeds  = 100000;					// user specified minimum number of reads to use as as duplicate seeds and off target alignments
const int cDfltDupSeeds = 5000000;					// default number of reads to use as as duplicate seeds and off target alignments
const int cMaxDupSeeds = 25000000;				    // user specified maximum number of reads to use as as duplicate seeds and off target alignments
const int cRSDMaxDupInstances = 2500;				// report instance counts up this max number of instances

const int cRSDMaxWorkerThreads = 64;				// can handle at most 64 threads
const int cRSDMaxInFileSpecs = 125;					// allow up to this many input files each for PE1 (or single ended) and PE2, 250 total max

const uint32_t cHashMask =       0x0fffff;			// hashing is over 10 bases
const uint32_t cMaxHashArrayEntries = (cHashMask + 1);        // alloc hash array to hold this many entries, must be at least 1 + maximal sized sized hash
const uint32_t cReallocPackedWrds = (((cMaxRSSeqLen+3)/4) * 1000);	// if needing to realloc memory for holding packed sampled reads then realloc this many additional words 

// processing mode enumerations
typedef enum TAG_eRSDMode
	{
	eRSDindependent = 0,							// process as independent readsets
	eRSDpooled										// process all readsets pooled 
	} etRSDMode;

#pragma pack(1)


typedef struct TAG_sSampledSeq {
	uint32_t NxtSeq;			// index+1 into m_pSampledSeqs[] at which next sampled sequence with same hash starts or 0 if no same hashed sequence 
	uint32_t NumInstances;		// number of instances of this sequence, or if PE then both sequences concatenated
	uint32_t NumRevCplInstances;	// number of instances where the match was with a reverse complemented SE or PE probe only
	uint16_t PE1ReadLen;		// PE1 sequence length
	uint16_t PE2ReadLen;		// if PE processing then PE2 sequence length
	uint16_t NumPackedWrds;	// number of UINT32s in PackedSeqs
	uint16_t Flags:16;		// to hold sundry flags
	uint32_t PackedSeqs[1];	// to contain the packed sequences, 16 bases per uint32_t
	} tsSampledSeq;

typedef struct TAG_sSeqCharacteristics {
	int64_t NumReads;				// number of reads processed
	int MeanReadLen;			// mean read length
	int MinReadLen;				// minimum read length
	int MaxReadLen;				// maximum read length
	int64_t NumSampled;			// number of reads sampled for duplicates and off target alignments
	int64_t NumCntd;				// number of reads which were counted in duplicate counts
	int64_t NumNoneCntd;			// number of reads which were processed for duplicate counts but not contributing towards these counts
	int64_t NotProcUL;			// number of reads not accepted for processing because after any end trimming they were < cMinRawSeqLen
	int64_t NotProcNs;			// number of reads not accepted for duplicate processing because they contained indeterminate bases
	int64_t NotProcQS;			// number of reads not accepted for duplicate processing because of low quality scoring
	int MinScore;				// minimum phred or quality score
	int MeanQScore;				// mean quality score
	int MaxScore;				// maximum phred or quality score
	} tsSeqCharacteristics;

typedef struct TAG_sInReadsFile {
	int FileID;					// uniquely identifies this raw reads file
	int ProcByThreadIdx;		// index of thread which has started to process this reads file; 0 if no thread yet to start processing
	char szFileName[_MAX_PATH];	 // reads were processed from this file
	int QSSchema;				// quality scoring schema - guestimated scoring schema - 0: no scoring, 1: Solexa, 2: Illumina 1.3+, 3: Illumina 1.5+, 4: Illumina 1.8+ or could be Sanger 
	int EstNumReads;			// estimated number of reads in this file
	int EstMeanReadLen;			// estimated mean read length
	tsSeqCharacteristics SeqCharacteristics;		// sequence characteristics for current reads file
} tsInReadsFile;


typedef struct TAG_sThreadNGSQCPars {
	int ThreadIdx;					// index of this thread (1..m_NumThreads)
	void *pThis;					// will be initialised to pt to CReadStats instance

#ifdef _WIN32
	HANDLE threadHandle;			// handle as returned by _beginthreadex()
	unsigned int threadID;			// identifier as set by _beginthreadex()
#else
	int threadRslt;					// result as returned by pthread_create ()
	pthread_t threadID;				// identifier as set by pthread_create ()
#endif
	int Rslt;						// returned result code
	int ProcessingID;				// processing instance identifier, used if processing eRSDindependent to identify output file instances
	int PE1RawReadLen;				// current length of read being buffered in PE1RawReadsBuff
	uint8_t PE1RawReadsBuff[cMaxRSSeqLen + 16];		// holds PE1 raw read sequence (16 is for a small safety margin!)
	int PE2RawReadLen;				// current length of read being buffered in PE2RawReadsBuff
	uint8_t PE2RawReadsBuff[cMaxRSSeqLen + 16];		// holds PE2 raw read sequence
	int PE1QScoresLen;				// current length of quality scores buffered in PE1QScoresBuff
	uint8_t PE1QScoresBuff[cMaxRSSeqLen + 16];		// holds PE1 raw read quality scores
	int PE2QScoresLen;				// current length of quality scores buffered in PE2QScoresBuff
	uint8_t PE2QScoresBuff[cMaxRSSeqLen + 16];		// holds PE2 raw read quality scores

	int NumInputFilesProcessed;		// number of input read files processed by this thread
	int64_t TotNumSEReads;			// total number of SE reads processed by this thread
	int64_t TotNumPEReads;			// total number of PE reads processed by this thread
	tsSeqCharacteristics SeqCharacteristics; // sequence characteristics for all reads processed by this thread
	uint32_t KMerCntOfs[cMaxRSSeqLen * cMaxKMerLen];  // each thread buffers offsets into m_pKMerCnts[] untill all K-mers in a read have been identified then updates m_pKMerCnts as an atomic block 
} tsThreadNGSQCPars;

typedef struct TAG_sThreadIndependentNGSQCPars {
	int ThreadIdx;					// index of this thread (1..m_NumThreads)
#ifdef _WIN32
	HANDLE threadHandle;			// handle as returned by _beginthreadex()
	unsigned int threadID;			// identifier as set by _beginthreadex()
#else
	int threadRslt;					// result as returned by pthread_create ()
	pthread_t threadID;				// identifier as set by pthread_create ()
#endif
	bool bProcCompleted;			// set true when this readset processing is complete
	int Rslt;						// returned result code
	etRSDMode PMode;				// processing mode; eRSDindependent or eRSDpooled
	int ProcessingID;				// processing instance identifier, used if processing eRSDindependent to identify output file instances
	bool bStrand;					// true if read strand specific distributions
	int Trim5;						// trim this number of bases from 5' end of reads when loading the reads
	int Trim3;						// trim this number of bases from 3' end of reads when loading the reads
	int MaxKMerLen;					// processing is for upto this KMer length inclusive
	int KMerCCC;					// concordance correlation coefficient measure KMer length
	int MaxContamSubRate;			// max allowed contamimant substitution rate (bases per 25bp of contaminant overlap, 1st 15bp of overlap no subs allowed)
	int MinContamLen;				// accept contaminant overlaps if overlap at least this many bases 
	int ReqMaxDupSeeds;				// requested to sample for this many duplicate seeds and off target alignments
	int  MinPhredScore;				// only accept reads for duplicate and KMer processing if their minimum Phred score is at least this threshold
	int NumThreads;					// number of worker threads to use
	bool bAffinity;					// thread to core affinity
	int NumPE1InputFiles;			// number of PE1 input files
	char **ppszInPE1files;			// input PE1 5' read files
	int NumPE2InputFiles;			// number of PE2 input files
	char **ppszInPE2files;		    // input PE2 3' read files
	char *pszContaminantFile;		// contaminants fasta file
	char *pszOutDistFile;			// where to write distributions CSV file
	char *pszOutHTMLFile;			//  where to write distributions HTML5 file
} tsThreadIndependentNGSQCPars;


#pragma pack()

class CReadStats
{
	etRSDMode m_PMode;			// processing mode; // processing mode; eRSDindependent or eRSDpooled

	int m_MaxKMerLen;			// processing is for upto this KMer length inclusive
	int m_KMerCCC;				// concordance correlation coefficient measure KMer length

	uint32_t *m_pBaseNs;			// to hold indeterminates Ns counts at each offset 5' to 3' along read length 
	uint32_t *m_pScores;			// to hold Phred scores at each offset 5' to 3' along read length
	int m_EstMaxSeqLen;			// estimated maximum read lengths 
	int m_AllocdMaxReadLen;		// current K-mer count allocations are for these maximal length reads
	int m_KMerCntsEls;			// total number of count elements per base
	size_t m_AllocdKMerCntsMem;	// allocated memory for K-mer counts
	uint32_t *m_pKMerCnts;		// to contain all K-mer counts along lengths of reads

	int m_MinReadLen;			// minimum length read processed
	int m_MaxReadLen;			// maximum length read processed
	uint32_t m_ReadLenDist[cMaxRSSeqLen+2]; // read length distributions, includes count m_ReadLenDist[cMaxRSSeqLen+1] of reads which were truncated to cMaxRSSeqLen
	uint64_t m_ProbNoReadErrDist[100];	// probabilities of read being error free distributions 

	int m_hReadLenDistRptFile;		// file handle for read length distributions report file
	char m_szReadLenDistRptFile[_MAX_PATH];	//  read length distributions report file


	int m_hKMerDistRptFile;		// file handle for Kmer distributions report file
	char m_szKMerDistRptFile[_MAX_PATH];	//  Kmer distributions report file

	int m_hPearsonDistRptFile;	// file handle for Pearson distribution report file
	char m_szPearsonDistRptFile[_MAX_PATH];	//  Pearson distributions report file

	int m_hQScoreDistRptFile;	// file handle for Phred quality score distribution report file
	char m_szQScoreDistRptFile[_MAX_PATH];	// Phred quality score distribution report file

	int m_hErrFreeReadDistRptFile;	// file handle for proportional error free read distribution report file
	char m_szErrFreeReadDistRptFile[_MAX_PATH];	// proportional error free read distribution report file


	int m_hDuplicatesDistRptFile;	// file handle for duplicate reads distribution report file
	char m_szDuplicatesDistRptFile[_MAX_PATH];	// duplicate reads distribution report file

	char m_szSVGFile[_MAX_PATH];  // SVG file for graph output

	bool m_bStrand;				// true if read strand specific distributions
	int m_Trim5;				// trim this number of bases from 5' end of reads when loading the reads
	int m_Trim3;				// trim this number of bases from 3' end of reads when loading the reads
	int m_ReqMaxDupSeeds;		// requested to sample for this many duplicate seeds and off target alignments
	int	m_ActMaxDupSeeds;		// actually sampled for this many duplicate seeds and off target alignments
	int m_MinMeanPhredScore;	// only accept reads for duplicate and KMer processing if their mean Phred score is at least this threshold
	int m_NumThreads;			// number of worker threads to use
	bool m_bAffinity;			// thread to core affinity
	int m_NumPE1InputFiles;		// number of PE1 input files
	char **m_ppszInPE1files;	// input PE1 5' read files
	int m_NumPE2InputFiles;		// number of PE2 input files
	char **m_ppszInPE2files;	// input PE2 3' read files

	int m_EstPE1MeanReadLen;	// PE1 estimated mean read length
	int m_EstPE2MeanReadLen;	// PE1 estimated mean read length

	char *m_pszContaminantFile;	// contaminants fasta file
	int m_hContamRptFile;		// file handle for contaminants ovelap report file
	char m_szContamRptFile[_MAX_PATH];	// contaminants ovelap report file
	int m_MaxContamSubRate;		// max allowed contamimant substitution rate (bases per 25bp of contaminant overlap, 1st 15bp of overlap no subs allowed)
	int m_MinContamLen;			// accept contaminant overlaps if overlap at least this many bases 
	CContaminants *m_pContaminates; // contaminants class for processing of contaminate sequences

	uint32_t m_NumChkdPE1ContamHits;	// number of PE1 sequences checked for contaminant hits
	uint32_t m_NumPE1ContamHits;	// number of PE1 sequences with contaminate hits
	uint32_t m_NumChkdPE2ContamHits;	// number of PE2 sequences checked for contaminant hits
	uint32_t m_NumPE2ContamHits;	// number of PE2 sequences with contaminate hits

	uint32_t m_NumSeqHashes;		// number of hashes in m_pSeqHashes
	uint32_t *m_pSeqHashes;		// array, indexed by sampled sequence hashes, holding offsets into packed sequence samples

	uint32_t m_AllocdSampledSeqWrds;  // this many sampled sequence words have been allocated
	uint32_t m_UsedSampledSeqWrds;	// this many sampled sequence words have been used 
	size_t m_AllocdSampledSeqMem;	// allocation for this many bytes
	uint32_t *m_pSampledSeqs;		// allocated to hold sampled sequences


	char *m_pszOutDistFile;		// where to write distributions CSV file
	char *m_pszOutHTMLFile;		//  where to write distributions HTML5 file

	bool m_bPEProc;				// false if SE, true if PE processing
	int m_NumInFiles;			// total number of input files in m_InReadsFiles
	int m_NumPE1InFiles;		// this number of PE1 input files in m_InReadsFiles
	int m_NumPE2InFiles;		// this number of PE2 input files in m_InReadsFiles
	tsInReadsFile m_InReadsFiles[cRSDMaxInFileSpecs * 2];  // to hold all input read file details ( 2* allows for PE reads)

	bool m_bTerminate;			// set true if processing threads are to terminate
	
	bool m_bMutexesCreated;			// set true if mutexes and rwlocks created/initialised
#ifdef _WIN32
	CRITICAL_SECTION m_hSCritSect;
	CRITICAL_SECTION m_hSCritSectScores;
	CRITICAL_SECTION m_hSCritSectKMers;
	HANDLE m_hMtxMHReads;
	SRWLOCK m_hRwLock;
#else
	pthread_spinlock_t m_hSpinLock;
	pthread_spinlock_t m_hSpinLockScores;
	pthread_spinlock_t m_hSpinLockKMers;
	pthread_mutex_t m_hMtxMHReads;
	pthread_rwlock_t m_hRwLock;
#endif

	void Init(void);			// initialisation
	void Reset(void);			// reset state back to that immediately following initialisation

	int CreateMutexes(void);			// create and initialise serialisation mutexes and locks
	void DeleteMutexes(void);			// delete serialisation mutexes and locks

	void AcquireSerialise(void);		// serialise access when determining seed sequences
	void ReleaseSerialise(void);
	void AcquireSerialiseScores(void);	// serialise access when updating score counts
	void ReleaseSerialiseScores(void);
	void AcquireSerialiseKMers(void);	// serialise access when updating KMer counts
	void ReleaseSerialiseKMers(void);
	void AcquireLock(bool bExclusive = false);		// defaults as read only lock
	void ReleaseLock(bool bExclusive = false);


	int
		LoadContamiantsFile(char *pszFile);	// loads contamiants fasta file into memory resident m_pContamSfx

	int				// 0 if contaminate duplicate by either name or sequence
		AddContaminant(char *pszName, etSeqBase *pSeq, int SeqLen);

	int				// < 0 if errors, otherwise the number of contaminent sequences
			FinaliseContaminants(void);

	int
		AddContamHash(uint32_t Hash, int SeqOf, int ContamSeq);

	int				// <0 if errors, 0 if no matches, >0 at least one contaminant sequence overlapping
		LocateContaminentMatch(int SeqLen,				// targ sequence is of this length
								etSeqBase *pSeq,		// target sequence		
								bool bPE2 = false);		// false if processing SE/PE1 read, true if PE2

	int				// returns 0 if not accepted as read instance, 1 if this is the first instance of an accepted sampled read, 2..N if multiple instances exist
		AddReadInst(tsThreadNGSQCPars *pThread, // thread specific processing state and context
				int PE1ReadLen,			// number of bases in PE1 read
				uint8_t *pPE1RawRead,		// PE1 read sequence
				int PE2ReadLen,			// number of bases in PE2 read
				uint8_t *pPE2RawRead);	// PE2 read sequence

	teBSFrsltCodes
		AnalyseReads(tsThreadNGSQCPars *pThread, // thread specific processing state and context
				int PE1ReadLen,			// number of bases in PE1 read
				uint8_t *pPE1RawRead,		// PE1 read sequence
				int PE1QScoreLen,		// number of quality scores in PE1 read (0 if no associated quality scores)
				uint8_t *pPE1QScores,		// PE1 read quality scores
				tsInReadsFile *pPE1File,// file containing PE1 or SE reads
				int PE2ReadLen = 0,		// number of bases in PE2 read
				uint8_t *pPE2RawRead = nullptr,	// PE2 read sequence
				int PE2QScoreLen = 0,		// number of quality scores in PE2 read (0 if no associated quality scores)
				uint8_t *pPE2QScores = nullptr,	// PE2 read quality scores
				tsInReadsFile *pPE2File = nullptr);	// file containing PE2

	int				// minimum score at any offset within the sequence
		AccumQScores(tsThreadNGSQCPars *pThread, // thread specific processing state and context
				int QSSchema,					// quality scoring schema - guestimated scoring schema - 0: no scoring, 1: Solexa, 2: Illumina 1.3+, 3: Illumina 1.5+, 4: Illumina 1.8+ or Sanger 
				int ReadLen,					// number of bases in read
				uint8_t *pRawRead,				// read sequence
				uint8_t *pQScores = nullptr);			// read quality scores (nullptr if no associated quality scores)

	// accumulate K-mer counts at each base offset along the read length
	int				
		AccumKMers(tsThreadNGSQCPars *pThread, // thread specific processing state and context
				int ReadLen,					// number of bases in read
				uint8_t *pSeq);					// read sequence

	teBSFrsltCodes
		ProcessReads(tsThreadNGSQCPars *pThread, // thread specific processing state and context
				tsInReadsFile *pPE1File, // file containing PE1 or SE reads
				tsInReadsFile *pPE2File);			 // file containing PE2 reads if PE processing

	double									// returned Pearson
		Pearsons(int KMerLen,				// KMer length for which counts were accumulated 
					uint32_t *pCtrlKMerCnts,  // KMer counts for control
					uint32_t *pExprKMerCnts);	// KMer counts for experiment	

public:
	CReadStats();
	~CReadStats();

	int ProcNGSQC(tsThreadNGSQCPars *pThread);		// read quality analytics thread processing

	int
		ProcessReadsetDist(etRSDMode PMode,			// processing mode; eRSDindependent or eRSDpooled
					int ProcessingID,				// processing instance identifier, used if processing eRSDindependent to identify output file instances
					bool bStrand,					// true if read strand specific distributions
					int Trim5,						// trim this number of bases from 5' end of reads when loading the reads
					int Trim3,						// trim this number of bases from 3' end of reads when loading the reads
					int MaxKMerLen,					// processing is for upto this KMer length inclusive
					int KMerCCC,					// concordance correlation coefficient measure KMer length
					int MaxContamSubRate,			// max allowed contamimant substitution rate (bases per 25bp of contaminant overlap, 1st 15bp of overlap no subs allowed)
					int MinContamLen,				// accept contaminant overlaps if overlap at least this many bases 
					int ReqMaxDupSeeds,				// requested to sample for this many duplicate seeds and off target alignments
					int  MinPhredScore,				// only accept reads for duplicate and KMer processing if their minimum Phred score is at least this threshold
					int NumThreads,					// number of worker threads to use
					bool bAffinity,					// thread to core affinity
					int NumPE1InputFiles,			// number of PE1 input files
					char *pszInPE1files[],			// input PE1 5' read files
					int NumPE2InputFiles,			// number of PE2 input files
					char *pszInPE2files[],		    // input PE2 3' read files
					char *pszContaminantFile,		// contaminants fasta file
					char *pszOutDistFile,			// where to write distributions CSV file
					char *pszOutHTMLFile);			//  where to write distributions HTML5 file



	};

