#pragma once

const int cNumAllocEls = 10000;			// alloc for tsChromEls in this many increments
const int cMaxNumBedFiles = 20;			// can handle at most this many source BED files
const int cMaxNumChromFeats = 200;		// can handle at most this many chromosomes

const int cMaxIncludeChroms = 20;		// max number of include chromosomes regular expressions
const int cMaxExcludeChroms = 20;		// max number of exclude chromosomes regular expressions

// BED regions of interest
typedef enum eBEDRegion {
	eMEGRAny = 0,				// process any region
	eMEGRIntergenic,			// only process intergenic
	eMEGRExons,					// only process exons
	eMEGRIntrons,				// only process introns
	eMEGRCDS,					// only process CDSs
	eMEGUTR,					// only process UTRs
	eMEG5UTR,					// only process 5'UTRs
	eMEG3UTR					// only process 3'UTRs
} etBEDRegion;

#pragma pack(1)
typedef struct TAG_sBED {
	uint8_t BedID;				// uniquely identifies this BED file
	char szBedFile[_MAX_PATH];	// BED file name
	uint32_t NumEls;				// number of elements loaded from this file
} tsBED;

typedef struct TAG_sFeatureEl {
	uint8_t BedID;				// BED file from which this feature was loaded
	uint32_t FeatureID;			// feature start/ends share same feature identifier
	uint32_t ChromID;				// chrom on which this feature is located
	uint32_t Loci;				// feature starts or ends at this loci (0..N)
	uint8_t FlgStart:1;			// will be 1 if feature start
	uint8_t FlgEnd:1;				// will be 1 if feature end
	char Strand;				// feature is on this strand
} tsFeatureEl;

typedef struct TAG_sChromFeatures {
	uint32_t ChromID;					// uniquely identifies this chromosome
	char szChrom[cMaxDatasetSpeciesChrom];
	int StartLoci;
	int EndLoci;
	int NumFeatures;						// current number of tsFeatureEls
	int  NumFeaturesAlloc;					// number of features allocated
	size_t NumFeaturesAllocMem;
	tsFeatureEl *pFeatures;					// features on this chromosome			
} tsChromFeatures;
#pragma pack()

class CChromFeatures
{
	char szOutFile[_MAX_PATH];			// name of BED file to contain generated merged features
	int m_hOutFile;						// file handle for generated BED file to contain merged features
	char m_LineBuff[8196];				// used for buffering output BED format merged features
	int m_LineBuffIdx;					// current index into m_LineBuff

	CBEDfile *m_pBedFile;				// BED file from which elements are currently being loaded
	tsBED *m_pCurBedFile;				// current bed
	int m_NumBedFiles;					// number of BED files from which elements were loaded
	tsBED m_BedFiles[cMaxNumBedFiles];  // elements were loaded from these BED files

	int m_NumChroms;					// number of chromosomes currently loaded
	tsChromFeatures *m_pCurLocChrom;	// last located chrom
	tsChromFeatures m_Chroms[cMaxNumChromFeats]; // hold chroms

	tsChromFeatures *LocateChrom(char *pszChrom); 

	CUtility m_RegExprs;            // regular expression processing

	char m_szFiltChrom[_MAX_PATH];	// used to cache last chrom processed	
	bool m_bFiltChrom;				// and it's filtered status

	void DumpFeat(char *pszMsg,tsFeatureEl *pFeat);

	int ReportMergedFeat(uint32_t FeatID,				// uniquely identifies this merged feature
						 char *pszChrom,			// feature is on this chrom
						 uint32_t MergeStartLoci,		// starting at this loci
						 uint32_t MergeEndLoci,		// ending at this loci
						 char Strand);				// and is on this strand

	static int SortFeatures(const void *arg1, const void *arg2);

public:
	CChromFeatures(void);
	~CChromFeatures(void);
	void Reset(void);
	bool ExcludeThisChrom(char *pszChrom);			 // returns true if chrom to be excluded from processing
	int SetChromFilters(int NumIncludeChroms,		 // number of chromosome regular expressions to include
						char **ppszIncludeChroms,	 // array of include chromosome regular expressions
						int NumExcludeChroms,		 // number of chromosome expressions to exclude
						char **ppszExcludeChroms);	 // array of exclude chromosome regular expressions);

	int AddFeature(char *pszBedFile,			// feature is from this BED file
						   char *pszChrom,		// and is on this chromosome
						   int StartLoci,		// starts at this loci
						   int EndLoci,			// and finishes at this loci
						   char Strand);		// and is on this strand

	int MergeFeatures(char *pszOutFile,			// file to which merged features are to be written
							  bool bStrandSpecific,		// if true then merging is strand specific
							  int MaxSep,		// merge features which are separated by up to this separation, if 0 then don't merge
							  int MinLen);		// only report merged features if they span at least this distance

	teBSFrsltCodes LoadBED(etBEDRegion Region,	// functional region of interest
				char ROIStrand,					// region on this strand
				 char *pszInFile);				// UCSC BED containing regions
	
};
