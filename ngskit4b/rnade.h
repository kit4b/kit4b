#pragma once

const double cNormCntsScale = 0.0;			// default normalise experiment to control counts is autoscaling

const double cClampFoldChange = 25.0;		// clamp fold changes to be at most cClampFoldChange

const int cMaxInFileSpecs = 100;				// allow at most a total of this many wildcarded control or experiment input read alignment loci files

const int cMinNumBins= 5;					// min number of bins
const int cDfltNumBins = 10;				// default number of bins
const int cMaxNumBins = 200;				// max number of bins

const int cMinStartLociThres = 1;			// user specified min number of unique start loci
const int cDfltStartLociThres = 5;			// default min number of unique start loci
const int cMaxStartLociThres = 200;			// user specified max number of unique start loci

const int cMinFeatCntThres = 1;				// user specified min feature read counts
const int cDfltFeatCntThres = 10;			// default min feature read counts
const int cMaxFeatCntThres = 200;			// user specified max min feature read counts

const int cMaxCoalesceWinLen = 20;			// maximum counts coalescing window length
const int cDfltCoalesceWinLen = 1;			// default counts coalescing window length

const int cMaxFeats2ProcAlloc = 200;		// allocate at most this many FeatureIDs to a thread for processing at any time

const int cWrtBinBuffSize = 0x03fffff;		// use a 4 MB write buffer size
const int cWrtStatBuffSize = 0x0fffff;		// use a 1 MB write buffer size

const int cAlignReadsLociInitalAlloc  = 30000000;	// initial allocation to hold this many read alignment loci
const int cAlignReadsLociRealloc	  = 15000000;	// realloc read alignments allocation in this sized increments
const int cDataBuffAlloc = 0x0fffffff;		// allocation size to hold gene features

// there are some assumptions here in terms of the max number of aligned read loci within the max length transcripts
// one day these assumptions may be shown to be invalid but let's see....
const int cMaxAssumTransLen = 2000000;					// assume no transcribed region will be longer than this
const int cMaxAssumTransLoci = cMaxAssumTransLen/10;	// assume very long transcribed regions will be low abundance reads and number of unique read aligned loci will be at most 10% of cMaxAssumTransLen

const int cMaxConfidenceIterations = 10000;	// max number of iterations when calculating confidence intervals and PValues

const int cMaxExclZones = 1000;			// max allowed number of exclusion zones within which reads are to be excluded

const int cPossion1SeqLen = 10000;
const int cPossion2SeqLen = 20000;
const int cPossion3SeqLen = 40000;
const int cPossion4SeqLen = 80000;
const int cPossion5SeqLen = 100000;
const int cPossion6SeqLen = 200000;
const int cPossion7SeqLen = 400000;
const int cPossion8SeqLen = 800000;
const int cPossion9SeqLen = 1000000;
const int cPossion10SeqLen = 2000000;


// following thresholds are used for characterisation of differential transcription state
// characterisation is into high, moderate, low and none with none including both negatively correlated and corelections less than cLoPeason
const double cHiPearsonThres = 0.8;					// >= this Pearson if control and experiment highly correlated in their reads distributions
const double cModPearsonThes = 0.5;					// >= this Pearson if control and experiment moderately correlated
const double cLoPearsonThres = 0.3;					// >= this Pearson if control and experiment have low correlated
const double cNoPearsonThres = cLoPearsonThres;		// < this Pearson if control and experiment have either negative or no correlated

const double cNoFoldChange = 1.25;				// if less than this fold change then characterise as none
const double cLoFoldChange = 1.50;				// if less than this fold change then characterise as low
const double cModFoldChange = 1.75;				// if less than this fold change then characterise as moderate
const double cHiFoldChange = cModFoldChange;	// if equal or more than this fold change then characterise as hi

// There are multiple processing phases
typedef enum TAG_eProcPhase {
	ePPInit = 0,								// initial initialisation
	ePPLoadFeatures,							// loading features of interest from BED file
	ePPLoadFeatClass,							// loading gene feature classifications
	ePPLoadExclZones,							// loading loci of reads which are to be excluded from any processing
	ePPLoadReads,								// loading read alignments from disk
	ePPReducePCRartifacts,						// reducing library PCR excessive read count artifacts
	ePPCoalesceReads,							// coalescing reads
	ePPNormLibCnts,								// normalise library counts to be very nearly the same
	ePPAllocDEmem,								// allocating memory for usage by threads and to hold feature DE fold changes etc
	ePPDDd,										// attempting to determine transcripts DE status for both counts and Pearson with PValue on transcript counts
	ePPReport,									// report to file transcript status, fold changes and Pearson
	ePPCompleted								// processing completed
	} etProcPhase;

typedef enum TAG_ePearsonScore {
	ePSIndeterminate = 0,						// insufficent counts, or other filtering criteria, by which Pearson can be calculated
	ePSNone,									// control and experiment have either negative or correlated below that of cLoPearsonThres
	ePSLow,										// lowly correlated Pearsons (cLoPearsonThres)
	ePSMod,										// moderately correlated Pearsons (cModPearsonThres)
	ePSHi										// highly correlated Pearsons (cHiPearsonThres)
} etPearsonScore;

typedef enum TAG_eCntsScore {
	eDEIndeterminate = 0,						// insufficent counts, or other filtering criteria, by which foldchange can be calculated
	eDEHi = 1,									// high changes in DE counts (cHiFoldChange)
	eDESMod,									// moderate changes in DE counts (cModFoldChange)
	eDSElow,									// low change in DE lowly correlated Pearsons (cLoFoldChange)
	eDESNone,									// no or very little change in DE counts (cNoFoldChange)
} etCntsScore;

// processing modes
typedef enum TAG_eDEPMode {
	eDEPMdefault = 0,				// Standard sensitivity (2500 iterations)
	eDEPMMoreSens,				// More sensitive (slower 5000 iterations)
	eDEPMUltraSens,				// Ultra sensitive (very slow 10000 iterations)
	eDEPMLessSens,				// Less sensitive (quicker 1000 iterations)
	eDEPMplaceholder				// used to set the enumeration range
	} etDEPMode;


typedef enum eDEBEDRegion {
	eDEMEGRAny = 0,				// process any region
	eDEMEGRExons,					// only process exons
	eDEMEGRIntrons,				// only process introns
	eDEMEGRCDS,					// only process CDSs
	eDEMEGUTR,					// only process UTRs
	eDEMEG5UTR,					// only process 5'UTRs
	eDEMEG3UTR					// only process 3'UTRs
} etDEBEDRegion;

// strand processing modes
typedef enum TAG_eStrandProc {
		eStrandDflt,			// default is to ignore the strand
		eStrandWatson,			// process for Watson or sense
		eStrandCrick,			// process for Crick or antisense
		eStrandPlaceholder
} etStrandProc;

#pragma pack(1)

typedef struct TAG_sExclZone {
	int RegionID;			// identifies this loci exclusion
	int ChromID;			// loci on this chrom, as generated by IDtoChrom()
	char Strand;		    // loci on this strand
	int StartLoci;			// loci starts inclusive
	int EndLoci;			// loci ends inclusive
} tsExclZone;

typedef struct TAG_sRefIDChrom {
	uint32_t ChromID;			// uniquely identifies chromosome
	uint32_t Hash;			// hash on chromosome name - GenHash24()
	uint32_t Next;			// used to link chroms with same Hash 
	char szChromName[cMaxDatasetSpeciesChrom+1];	// chromosome name
} tsRefIDChrom;

typedef struct TAG_sAlignReadLoci {
	uint8_t ExprFlag:1;			// 0 if aligned read from control, 1 if from experiment
	uint8_t Sense:1;				// 0 if '-' or antisense crick, 1 if '+' or sense watson
	uint8_t FileID:5;				// from which reads alignment file was this alignment loaded
	uint32_t NormCnts;			// coalesced counts with control vs experiment library size normalisation
	uint32_t ArtCnts;				// used as temp count storage during artefact reduction processing
	uint32_t AlignHitIdx;			// current read hit index + 1 for this read
	uint32_t ChromID;				// identifies hit chromosome - ChromToID(pszChrom)
	uint32_t Loci;				// 5' start loci on chromosome of this alignment
	uint32_t ReadLen;				// read length
} tsAlignReadLoci;

typedef struct TAG_sAlignLociInstStarts {
	uint32_t Bin;					// instance counts are for this bin
	uint32_t RelLoci;				// instance counts are at this loci
	uint32_t NumCtrlStarts;		// number of control start loci instances at this loci instance
	uint32_t NumExprStarts;		// number of experiment start loci instances at this loci instance
} tsAlignLociInstStarts ;

typedef struct TAG_sAlignBin {
	uint32_t Bin;					// poisson counts are from this bin (1..n)
	uint32_t BinRelStartLoci;		// bin starts at this relative loci
	uint32_t BinRelEndLoci;		// bin starts at this relative loci
	uint32_t NumCtrlInstStarts;	// number of control start loci instances in this bin
	uint32_t NumExprInstStarts;	// number of experiment start loci instances in this bin
	uint32_t ControlCnts;			// original total control counts in this bin over control start loci instances
	uint32_t ExperimentCnts;		// original total experiment counts in this bin over experiment start loci instances
	uint32_t ControlCoverage;     // coverage apportioned to this bin (control reads likely to overlap multiple bins)
	uint32_t ExperimentCoverage;  // coverage apportioned to this bin (experiment reads likely to overlap multiple bins)
	uint32_t ControlPoissonCnts;	 // after poisson control counts
	uint32_t ExperimentPoissonCnts;// after poisson experiment counts
} tsAlignBin;

typedef struct TAG_sFeatDE {
		char szFeatName[128];	// feature name
		int FeatLen;			// feature transcribed length
		int NumExons;			// number of exons
		int UserClass;			// user classification for this feature
		int DEscore;			// transformed DE score (range 0..9)
		int CntsScore;			// score for DE cnts (range 0..4)
		int PearsonScore;		// score for Pearson (range 0..4)
		int CtrlCnts;			// control read counts
		int ExprCnts;			// expression read counts
		int SumCtrlExprCnts;	// sum control + expression counts
		double PValueMedian;	// median P-value (0..1.0)
		double PValueLow95;		// low 95 P-value (0..1.0)
		double PValueHi95;		// high 95 P-value (0..1.0)
		double ObsFoldChange;	// observed fold change (0..n) if less than 1.0 then negative fold change
		double FoldMedian;		// fold change median (0..50.0) if less than 1.0 then negative fold change
		double FoldLow95;		// fold change low 95 percentile (0..50.0) if less than 1.0 then negative fold change
		double FoldHi95;		// fold change high 95 percentile (0..50.0) if less than 1.0 then negative fold change
		double PearsonObs;		// Pearson observed	(-1.0 to 1.0)
		double PearsonMedian;	// Pearson median (-1.0 to 1.0)
		double PearsonLow95;	// Pearson low 95 percentile (-1.0 to 1.0)
		double PearsonHi95;		// Pearson high 95 percentile (-1.0 to 1.0)
		int TotCtrlStartLoci;	// total number of unique control start loci in this transcript
		int TotExprStartLoci;	// total number of unique experiment start loci in this transcript
		int BinsShared;			// number of bins with both control and experiment counts
		int BinsExclCtrl;		// number of bins holding counts exclusive to control
		int BinsExclExpr;		// number of bins holding counts exclusive to experiment
		uint32_t BinsCtrlDepth[cMaxNumBins];	// to hold read coverage depth in each control bin
		uint32_t BinsExprDepth[cMaxNumBins];	// to hold read coverage depth in each experiment bin
} tsFeatDE;


// each thread has it's own instance of the following
typedef struct TAG_ThreadInstData {
	uint32_t ThreadInst;				// uniquely identifies this thread instance
	void *pThis;					// will be initialised to pt to class instance

#ifdef _WIN32
	HANDLE threadHandle;			// handle as returned by _beginthreadex()
	unsigned int threadID;			// identifier as set by _beginthreadex()
#else
	int threadRslt;					// result as returned by pthread_create ()
	pthread_t threadID;				// identifier as set by pthread_create ()
#endif
	CStats *pStats;						// used for ChiSquare processing
	teBSFrsltCodes Rslt;				// thread processing completed result - eBSFSuccess if no errors
	int FeatureID;						// currently processing this feature

	uint32_t CurFeatLen;					// current feature length over which counts are being processed
	int CurRegionLen;					// region length for curent gene or feature being processed

	double *pPValues;					// to hold PValues for transcript
	double *pFeatFoldChanges;			// to hold feature fold changes (Ctrl + 0.001)/(Expr + 0.001)) whilst determining feature confidence interval
	double *pPearsons;					 // prealocated to hold Pearsons whilst determining confidence interval
	uint32_t NumBinsWithLoci;				// pAlignBins currently contains this number of bin instances with at least one aligned loci instance
	tsAlignBin *pAlignBins;				// preallocated to hold alignment bin cnts reserved for this thread

	tsAlignBin *pPoissonAlignBins;		// to hold alignment bins with poisson applied
	uint32_t NumBinInstStarts;				// m_pBinInstsCnts currently contains this many experiment and control start instances
	tsAlignLociInstStarts *pBinLociInstStarts;	// preallocated pts to list of control and experiment start instance counts for all bins
	int MaxFeats2Proc;					// max number of FeatureIDs which can be allocated for processing by this thread by GetFeats2Proc() into Feats2Proc[]
	int NumFeats2Proc;					// process this number of features in Feats2Proc[]
	int NumFeatsProcessed;				// number features processed currently from Feats2Proc
	int Feats2Proc[cMaxFeats2ProcAlloc];	// these are the feature identifiers for the features to be processed
	CSimpleRNG *pSimpleRNG;				// random generator exclusively for use by this thread
} tsThreadInstData;

#pragma pack()

class CRNA_DE {

	bool m_bDEMutexesCreated;				// true if mutexes for serialisation have been created
	int m_NumDEThreads;							// number of processing threads
	int m_MaxConfidenceIterations;				// max number of iterations over start loci when inducing random poisson noise
	tsThreadInstData *m_pThreadsInstData;		// all allocated thread instance data

	etDEPMode m_DEPMode;					// processing mode
	etProcPhase m_ProcessingPhase;		// current processing phase
	bool m_bFiltNonaligned;				// true if only features having at least one read aligned are to be be reported

	int m_NumExclZones;					// total number of read exclusion zones loaded
	int m_NumExclReads;					// total number of reads which were excluded because they overlaid an exclusion zone
	tsExclZone *m_pExclZones;			// pts to exclusion zones to which overlaying reads are to be excluded from processing

	uint32_t m_LimitAligned;				// for test/evaluation can limit number of reads parsed to be no more than this number (0 for no limit)
	int m_CoWinLen;						// counts coalescing window length

	tsAlignReadLoci *m_pCtrlAlignReadLoci; // memory allocated to hold control read alignment loci, reads are written contiguously into this memory
	uint32_t m_AllocdCtrlAlignReadsLoci;			  // how many instances of control tsAlignReadLoci have been allocated
	size_t m_AllocdCtrlAlignReadsMem;			 // size of allocated memory
	uint32_t m_CurNumCtrlAlignReadsLoci;			  // m_pAlignReadLoci currently contains a total of this many control read alignment loci
	uint64_t m_CurSumCtrlReadsLen;				  // current summed control reads length

	tsAlignReadLoci *m_pExprAlignReadLoci;	// memory allocated to hold experiment read alignment loci, reads are written contiguously into this memory
	uint32_t m_AllocdExprAlignReadsLoci;			// how instances of experiment tsAlignReadLoci have been allocated
	size_t m_AllocdExprAlignReadsMem;			 // size of allocated memory
	uint32_t m_CurNumExprAlignReadsLoci;			// m_pAlignReadLoci currently contains a total of this many experiment read alignment loci
	uint64_t m_CurSumExprReadsLen;				// current summed control reads length

	uint32_t m_NumLoadedCtrlReads;			// total number of control reads actually loaded prior to any  coalescing and library size normalisation
	uint32_t m_NumLoadedExprReads;			// total number of expression reads loaded loaded prior to any coalescing and library size normalisation
	uint32_t m_MeanLenCtrlReads;				// mean length of loaded control reads
	uint32_t m_MeanLenExprReads;				// mean length of loaded experiment reads


	uint32_t m_NumNormCtrlReads;				// total number of control reads after library size normalisation
	uint32_t m_NumNormExprReads;				// total number of expression reads after library size normalisation

	uint32_t m_NumPoissonNormCtrlReads;		// total number of control reads after poisson noise
	uint32_t m_NumPoissonNormExprReads;		// total number of expression after after poisson noise

	uint32_t m_NumFeaturesLoaded;				// total number of features loaded for processing
	uint32_t m_NumFeaturesProcessed;			// current number of features processed

	tsAlignBin *m_pAlignBins;		// memory allocated to hold alignment bin cnts
	uint32_t m_AllocdAlignBins;			// how instances of tsAlignBin have been allocated
	tsAlignLociInstStarts *m_pBinInstsStarts;	// pts to list of control and experiment start instance counts for all bins
	uint32_t m_AllocBinInstStarts = 0;			// m_pBinInstsCnts is currently allocated to hold at most this number of start loci instance counts

	int m_FeatsPerThread;					// number of features to be processed as a block by each thread

	double *m_pPearsons;					// to hold Pearsons whilst determining confidence interval
	size_t m_AllocNumPearsons;				// number allocated


	tsThreadInstData *m_pThreadInst;		// pts to current thread instance

	char m_DEAlignStrand;				// process for reads on this strand only
	char m_FeatStrand;					// process for genes or features on this strand only
	etDEBEDRegion m_Region;				// process for this genomic region only
	int m_NumBins;						// Bin regions into this many non-overlapping bins

	double m_LibSizeNormExpToCtrl;		// factor by which experiment counts can be normalised to that of control counts (accounts for library size ratio)

	CBEDfile *m_pBEDFeatFile;

	size_t m_AllocBinWrtBuff = 0;		// memory allocated for output bin write buffers
	size_t m_AllocStatsWrtBuff = 0;		// memory allocated for output stats write buffers
	size_t m_WrtBinBuffOfs = 0;			// offset at which to next write
	size_t m_WrtStatsBuffOfs = 0;		// offset at which to next write
	uint8_t *m_pWrtBinBuff;				// used to buffer output writes to bin counts file
	uint8_t *m_pWrtStatsBuff;				// used to buffer output writes to stats file
	int m_hOutStatsFile;				// output stats results file
	bool m_bWrtStatHdr;					// true if stats file requires header
	int m_hOutBinFile;					// output bin counts file
	bool m_bWrtBinHdr;					// true if bin counts file requires header row

	double *m_pFeatFoldChanges;			// to hold feature fold changes (Ctrl + 0.001)/(Expr + 0.001)) whilst determining feature confidence interval
	size_t m_AllocNumFeatFoldChanges;	// number allocated

	double *m_pPValues;					// to hold feature PValues
	size_t m_AllocNumPValues;			// number allocated

	tsAlignBin *m_pPoissonAlignBins;	 // to hold bin poisson counts whilst determining pearson confidence interval
	size_t m_AllocNumPoissonCnts;		 // number allocated

	CSimpleRNG m_SimpleRNG;				// used to generate Poisson distributed random cnts

	int m_NumFeatsDEd;					// number of features processed into tsFeatDE's
	int m_AllocdFeatsDEd;				// number alloc'd
	tsFeatDE *m_pFeatDEs;				// allocated to hold features processed for DE

	int m_NumNoneDEd;					// number of features which are candidates for determining as being none DE'd in phase ePPAllocDEmem

	int m_MinFeatCntThres;				// minimum feature count threshold, control or experiment, required (1 to 200, defaults to 20)
	int m_MinStartLociThres;			// minimum feature unique start loci threshold, control or experiment, required (1 to 200, defaults to 10)


	tsRefIDChrom *m_pChroms = NULL;
	int m_NumChromsAllocd = 0;
	int m_CurNumChroms = 0;
	int m_MRAChromID = 0;
	uint32_t *m_pChromHashes;			// allocated to provide hash index for chroms

	int m_LastFeatureAllocProc = 0;					// set to the last FeatureID allocated to a thread for processing, 0 if starting feature allocation, -1 if all features have been allocated

	int m_Possion1[cPossion1SeqLen];
	int m_Possion2[cPossion2SeqLen];
	int m_Possion3[cPossion3SeqLen];
	int m_Possion4[cPossion4SeqLen];
	int m_Possion5[cPossion5SeqLen];
	int m_Possion6[cPossion6SeqLen];
	int m_Possion7[cPossion7SeqLen];
	int m_Possion8[cPossion8SeqLen];
	int m_Possion9[cPossion9SeqLen];
	int m_Possion10[cPossion10SeqLen];




	int		DECreateMutexes(void);
	void	DEDeleteMutexes(void);
	inline void	DEAcquireSerialise(void);
	inline void ReleaseSerialise(void);
	
	double ClampFoldChange(double Scale);

	int	LoadGeneFeatClasses(char *pszFeatClassFile); // // Loads gene or feature classification from a CSV or tab delimited file
	int	LoadExclZones(char *pszExclZonesFile); // Loads read start/end loci for which any covering reads are to be excluded
	teBSFrsltCodes 	LoadGeneFeatures(char Strand,			// features on this strand
					 char *pszInFeatFile);	// from this BED file

	int	LoadAlignedReads(bool bExpr,	// false if loading control, true if loading experiment
			char Strand,		// process for this strand '+' or '-' or for both '*'
			int FType,				// input element file format: 0 - auto, 1 - CSV, 2 - BED, 3 - SAM/SAM
			int FileID,			// uniquely identifies this input file
			char *pszInAlignFile);	// load aligned reads from this file)

	int
	LoadAlignedReadFiles(char Strand,		// process for this strand '+' or '-' or for either '*'
				int FType,					// input element file format: 0 - auto, 1 - CSV, 2 - BED, 3 - SAM/BAM
				int NumInputControlSpecs,	// number of input file specs
				char **pszInControlFiles,	// input control aligned reads files
				int NumInputExperimentSpecs,// number of input file specs
				char **pszInExperimentFiles);	// input experiment aligned reads files

	teBSFrsltCodes				// Load aligned reads from SAM/BAM alignment files
		LoadAlignedReadsSAM(bool bIsExperiment,		// false if control file, true if experiment
		int FileID,						// uniquely identifies this file
		char* pszInFile,
		char FiltStrand);

	teBSFrsltCodes 	LoadAlignedReadsBED(bool bIsExperiment,		// false if control file, true if experiment
				int FileID,						// uniquely identifies this file
				char *pszInFile,
				char FiltStrand);

	teBSFrsltCodes 	LoadAlignedReadsCSV(bool bIsExperiment,		// false if control file, true if experiment
				int FileID,						// uniquely identifies this file
				char *pszInFile,
				char FiltStrand);



	teBSFrsltCodes ProcessReads4DE(void);

	teBSFrsltCodes						// if > eBSFSuccess then read was silently sloughed because it was in an exclusion loci
		AddReadHit(int FileID,				// parsed from this file
			bool bIsExperiment,		// true if this is an experimental read
			char *pszChrom,			// hit to this chrom
			char Strand,			// on this strand
			int StartLoci,			// starts at this loci
			int ReadLen);			// and is of this length

	int	LocateStartAlignment(char Strand,uint32_t ChromID,uint32_t StartLoci,uint32_t EndLoci, uint32_t NumAlignReadsLoci,tsAlignReadLoci *pAlignReadLoci);


	int CoalesceReadAlignments(int WinLen,				// coalescing window length 1..20
						   bool bSamesense,			// if true then only coalesce reads with same sense
						   bool bExperiment);		// if true then coalesce experiment read loci otherwise control read loci

	teBSFrsltCodes ReducePCRartifacts(int FlankLen,				// 5' and 3' flank length over which the background mean rate is to be measured
				int ArtifactCntsThres);				// if counts at any loci are >= this threshold then process for PCR artifact reduction

	uint32_t										// reduced total counts by this
		ReducePCR(int FlankLen,						// 5' and 3' flank length over which the background mean rate is to be measured
				int ArtifactCntsThres,				// if counts at any loci are >= this threshold then process for PCR artifact reduction
				uint32_t NumAlignReadsLoci,tsAlignReadLoci *pAlignReadLoci);

	teBSFrsltCodes NormaliseLibraryCounts(void);	// normalise library sizes

	void InitPoissonSeqs(tsThreadInstData *pThreadInst); // GenPoissonSeqs
	int RandPoisson(tsThreadInstData *pThreadInst,int Lambda);
	teBSFrsltCodes GenPoissonLibraryCounts(void);   // generate poisson noise over all library counts

	double									// returned Pearson sample correlation coefficient
		Pearsons(tsAlignBin *pAlignBins);		// bins containing alignment counts

		double								    // returned median PValue
		PearsonsPValue(tsThreadInstData *pThreadInst,	// thread instance
			double Pearson,					// observed pearson
			int MaxPerms,					// maximum number of permutions
			tsAlignBin *pAlignBins,			// bins containing alignment counts
			double *pPValueLow95,			// returned low 95 percentile
			double *pPValueHi95,			// returned upper 95 percentile
			double *pLow95,					// returned low 95 percentile
			double *pHi95,					// returned upper 95 percentile
			double *pMedian,				// returned median
			double *pFeatLow95,				// returned low 95 percentile for feature
			double *pFeatHi95,				// returned upper 95 percentile for feature
			double *pFeatMedian);			// returned median for feature

	double										// returned Pearson
		PoissonPearsons(tsAlignBin *pAlignBins);		// bins containing alignment counts

	double								// returns *pUpper - *pLower
		ConfInterval95(int N,				// number of bins containing at least one count
			   double Pearson,		// Pearsons r
			   double *pUpper,		// returned upper for a confidence interval of 95
			   double *pLower);		// returned lower for a confidence interval of 95

	double								// returns *pUpper - *pLower
		ConfInterval99(int N,		// number of bins containing at least one count
			   double Pearson,		// Pearsons r
			   double *pUpper,		// returned upper for a confidence interval of 95
			   double *pLower);		// returned lower for a confidence interval of 95
	double	poz (double	z);						/* returns cumulative probability from -oo to z */
	double z2r(double z); // Convert Fisher z' to Pearson r
	double r2z(double r); // Convert Pearson r to Fisher z
	
	int
		AddDEPearsons(tsThreadInstData *pThreadInst,
				char *pszFeatName,        	// counts DE and Pearson is for this feature
				int NumExons,				// feature has this number of exons
				int UserClass);				// user classification for this feature



	teBSFrsltCodes
		GenBinAlignStarts(tsThreadInstData *pThreadInst,	// thread instance
					  uint32_t RegionOfs,			// StartLoci is at this region offset
					  char *pszChrom,			// alignments are on this chrom
					  uint32_t StartLoci,			// must start on or after this loci
					  uint32_t EndLoci);			// must start on or before this loci


	teBSFrsltCodes
		AddAlignBinCnts(tsThreadInstData *pThreadInst,	// thread instance data
				int RelLoci,			// relative offset within  pThreadInst->CurRegionLen from which these cnts were derived
				int MeanControlReadLen,	// control reads were of this mean length
				bool bSense,				// 1 if aligned sense, 0 if antisense
				int ControlCnts,		// attribute this many control counts (alignments) to relevant bin(s)
				int MeanExperimentReadLen,	// experiment reads were of this mean length
				int ExperimentCnts);		// attribute this many experiment counts (alignments) to relevant bin(s)

	int ReportDEandPearsons(void);
	int ReportDEandPearsonBinCounts(void);

	char *IDtoChrom(uint32_t ChromID);			// returns ptr to chrom for ChromID
	uint32_t					  // unique identifier or 0 if chromosome not known
		ChromToID(char *pszChrom, // get unique chromosome identifier for this chrom
				bool bAdd);	  // if true and chromosome not already known then add and return it's identifier



	int  // Returns 0 if Ctrl within Delta of Expr, -1 if Ctrl more than Delta below Expr, 1 if Ctrl more than Delta above Expr
		CmpLoose(double Delta,double Ctrl, double Expr);

	bool												// true if at least one FeatureID to be processed, false if all FeatureIDs have been already allocated for processing
		GetFeats2Proc(tsThreadInstData *pThreadInst);	// thread instance data
	int	ProcessFeature(tsThreadInstData *pThreadInst);	// instance data for this individual thread - could be multiple threads calling this method
	teBSFrsltCodes	GenBinStarts(tsThreadInstData *pThreadInst,	// thread instance - could be multiple threads calling this method
				  uint32_t RegionOfs,			// StartLoci is at this rlative offset within the current region
				  uint32_t ChromID,			// alignments are on this chrom
				  uint32_t StartLoci,			// must start on or after this loci
				  uint32_t EndLoci,			// must start on or before this loci
				  uint32_t NumAlignReadsLoci,
				  tsAlignReadLoci *pAlignReadsLoci);

	tsAlignLociInstStarts *
		UpdateBinLociInstStarts(tsThreadInstData *pThreadInst,			// thread instance
					uint32_t Bin,								// starts are in this bin
					int RelLoci,							// starts at this relative (to start of binned transcript) loci
					uint32_t CtrlStarts,						// this number of control reads start at RelLoci
					uint32_t ExprStarts);						// this number of experiment reads start at RelLoci

	teBSFrsltCodes AddAlignBinCnts(tsThreadInstData *pThreadInst,	// thread instance data
			uint32_t RelLoci,			// relative offset within  pThreadInst->CurRegionLen from which these cnts were derived
			uint32_t MeanControlReadLen,	// control reads were of this mean length
			bool bSense,				// 1 if aligned sense, 0 if antisense
			uint32_t ControlCnts,		// attribute this many control counts (alignments) to relevant bin(s)
			uint32_t MeanExperimentReadLen,	// experiment reads were of this mean length
			uint32_t ExperimentCnts);		// attribute this many experiment counts (alignments) to relevant bin(s)

	static int SortAlignments(const void *arg1, const void *arg2);
	static int SortDoubles(const void *arg1, const void *arg2);
	static int SortDEScore(const void *arg1, const void *arg2);
	static int SortFoldMedian(const void *arg1, const void *arg2);

	CMTqsort m_mtqsort;					// multithreaded quicksort

#ifdef _WIN32
	CRITICAL_SECTION m_hSCritSect;	// used to serialise
	unsigned __stdcall ThreadedDEproc(void * pThreadPars);
#else
	pthread_spinlock_t m_hSpinLock;
	void *ThreadedDEproc(void * pThreadPars);
#endif

public:
	CRNA_DE();	// constructor
	~CRNA_DE();	// destructor

	void Reset(void);

	static char *Region2Txt(etDEBEDRegion Region);
	static char ReportStrand(etStrandProc StrandProc);
	int	ProcThreadRNADE(tsThreadInstData* pThreadPar);	// worker thread parameters
	teBSFrsltCodes Process(etDEPMode PMode,					// processing mode
						int NumThreads,						// number of threads (0 defaults to number of CPUs)
						int  CoWinLen,						// counts coalescing window length
						int ArtifactCntsThres,				// if counts at any loci are >= this threshold then process for PCR artifact reduction
						uint32_t LimitAligned,				// for test/evaluation can limit number of reads parsed to be no more than this number (0 for no limit)
						bool bFiltNonaligned,				// true if only features having at least one read aligned are to be be reported
						char AlignStrand,					// process for reads aligning to this strand only
						char FeatStrand,					// process for genes or features on this strand only
						etDEBEDRegion Region,					// process for this genomic region only
						int	NumBins,						// number of non-overlapping count bins
						int MinFeatCntThres,				// minimum feature count threshold, control or experiment, required (1 to 200, defaults to 20)
						int MinStartLociThres,				// minimum feature unique start loci threshold, control or experiment, required (1 to 200, defaults to 10)
						double NormCntsScale,				// counts normalisation scale factor
						int FType,							// input element file format: 0 - auto, 1 - CSV, 2 - BED, 3 - SAM/BAM
						int NumInputControlSpecs,			// number of input file specs
						char **pszInControlFiles,			// input control aligned reads files
						int NumInputExperimentSpecs,		// number of input file specs
						char **pszInExperimentFiles,		// input experiment aligned reads files
						char *pszInFeatFile,				// input gene or feature BED file
						char *pszFeatClassFile,				// classify genes or features from this file
						char *pszExclZonesFile,				// exclude reads overlaying zone loci specified from this file
						char *pszOutfile,					// output into this file
						char *pszBinCountsFile);			// output bin counts to this file
	};
