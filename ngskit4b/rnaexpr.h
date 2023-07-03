#pragma once


const int32_t cMaxSampleRows = 100000;		// expecting at most this number of RNA sample rows
const int32_t cMaxRNAESamples = 8000;		// limit number of RNA samples (RNASampleNames) to be at most this many
const int32_t cMaxWGSSamples =  2000;		// limit number of WGS samples (WGSSampleNames) to be at most this many
const int32_t cMaxMaterials =  2000;		// limit number of materials (MaterialNames) to be at most this many
const int32_t cMaxRNAEFeatures = 1000000;	// limit number of features (FeatureNames) to be at most this many
const int32_t cMaxChromNames = 100000;		// max number of unique chromsome names allowed
const int32_t cMaxSampleRefNameLen = 100;	// sample reference names can be at most this length
const int32_t cMaxFeatRefNameLen = 200;		// needs to be quite long as may be a concatenation of chrom + loci + species/cultivar
const int32_t cMaxMaterialNameLen = 100;	// material names can be at most this length
const int32_t cMaxRNAEWorkerThreads = 10;	// limit number of threads to at most this many

#define FEATVALUETYPE int32_t		// currently processing is for feature values representing integer counts, future processing may include other types

typedef enum TAG_eModeRNAE {
	eMRNAEDefault = 0,		// processing mode 0: check RNA replicate samples for expression and allelic inconsistencies
	eMRNAEHomozScores,		// check biological sample replicates for allelic homozygosity score inconsistencies
	eMRNAEPlaceHolder		// used to mark end of processing modes
	}eModeRNAE;

#pragma pack(1)

typedef struct TAG_sSampleFeats {
	size_t Nxt;            // offset to next SampleFeat in linked list 
	int32_t SampleRefOfs;  // offset in m_pSampleFeatures of this sample reference name
	int32_t Features[1];   // feature values for this sample
	} tsSampleFeats;

typedef struct TAG_sSampleFeatValue {
	int32_t SampleID;      // identifies the sample
	int32_t FeatID;        // identifies the feature
	int32_t Value;         // SampleID.FeatID has this value
	} tsSampleFeatValue;

typedef struct TAG_sRNAGBSWGSMaterial {
	uint32_t RNASampleID;	// RNA sample name identifier
	uint32_t WGSSampleID;	// WGS sample name identifier
	uint32_t MaterialID;	// Material name identifier
	uint32_t LocationID;	// RNA location
	uint32_t EntryBookID;	// RNA entry name
	uint32_t SampleNum;			// RNA replicate number - currently just 1 or 2
	uint32_t PlotNum;			// RNA plot number
	uint32_t RangeNum;			// RNA range number
	uint32_t RowNum;				// RNA row number
	} tsRNAGBSWGSMaterial;

#pragma pack()

class CRNAExpr
{
	uint32_t m_LARNAMetaNameID;			// name identifier last returned by AddRNAMetaName()
	uint32_t m_NumRNAMetaNames;			// number of RNAMetaNames names currently in m_szRNAMetaNames
	uint32_t m_NxtszRNAMetaNameIdx;			// current concatenated (names separated by '\0') size of all Meta names in m_szRNAMetaNames - offset at which to add new name
	char m_szRNAMetaNames[cMaxRNAESamples * (2*cMaxSampleRefNameLen + 1)];	// used to hold concatenated RNAMetaName (Location,EntryBook names, each separated by '\0'
	uint32_t m_szRNAMetaNamesIdx[2*cMaxRNAESamples];	// array of indexes into m_szRNAMetaNames giving the start offsets of each name - location or entry boo;

	uint32_t m_LARNASampleNameID;			// name identifier last returned by AddRNASampleName()
	uint32_t m_NumRNASampleNames;			// number of RNASampleNames names currently in m_szRNASampleNames
	uint32_t m_NxtszRNASampleNameIdx;			// current concatenated (names separated by '\0') size of all SampleName names in m_szRNASampleNames - offset at which to add new name
	char m_szRNASampleNames[cMaxRNAESamples * (cMaxSampleRefNameLen + 1)];	// used to hold concatenated RNASampleName names, each separated by '\0'
	uint32_t m_szRNASampleNamesIdx[cMaxRNAESamples];	// array of indexes into m_szRNASampleNames giving the start offsets of each name

	uint32_t m_LAWGSSampleNameID;			// name identifier last returned by AddWGSSampleName()
	uint32_t m_NumWGSSampleNames;			// number of WGSSampleNames names currently in m_szWGSSampleNames
	uint32_t m_NxtszWGSSampleNameIdx;			// current concatenated (names separated by '\0') size of all SampleName names in m_szWGSSampleNames - offset at which to add new name
	char m_szWGSSampleNames[cMaxWGSSamples * (cMaxSampleRefNameLen + 1)];	// used to hold concatenated WGSSampleName names, each separated by '\0'
	uint32_t m_szWGSSampleNamesIdx[cMaxWGSSamples];	// array of indexes into m_szWGSSampleNames giving the start offsets of each name

	uint32_t m_LAMaterialNameID;			// name identifier last returned by AddMaterialName()
	uint32_t m_NumMaterialNames;			// number of MaterialNames names currently in m_szMaterialNames
	uint32_t m_NxtszMaterialNameIdx;			// current concatenated (names separated by '\0') size of all MaterialName names in m_szMaterialNames - offset at which to add new name
	char m_szMaterialNames[cMaxMaterials * (cMaxMaterialNameLen + 1)];	// used to hold concatenated MaterialName names, each separated by '\0'
	uint32_t m_szMaterialNamesIdx[cMaxMaterials];	// array of indexes into m_szMaterialNames giving the start offsets of each name


	int32_t m_LAChromNameID;						// last accessed chromosome identifier from call to AddChrom()
	int32_t m_NumChromNames;						// number of chromosome names currently in m_szChromNames
	int32_t m_NxtszChromNameIdx;					// current concatenated (names separated by '\0') size of all chromosome names in m_szChromNames - offset at which to add new name
	char m_szChromNames[cMaxChromNames * (cMaxDatasetSpeciesChrom + 1)];	// used to hold concatenated chromosome names, each separated by '\0'
	uint32_t m_szChromNamesIdx[cMaxChromNames];			// array of indexes into m_szChromNames giving the starts of each chromosome name

	int32_t m_LAFeatureNameID;						// last accessed feature name identifier from call to AddFeatureName()
	int32_t m_NumFeatureNames;						// number of feature names currently in m_szFeatureNames
	int32_t m_NxtszFeatureNameIdx;					// current concatenated (names separated by '\0') size of all feature names in m_szFeatureNames - offset at which to add new name
	char m_szFeatureNames[cMaxRNAEFeatures * (cMaxFeatRefNameLen + 1)];	// used to hold concatenated chromosome names, each separated by '\0'
	uint32_t m_szFeatureNamesIdx[cMaxRNAEFeatures];			// array of indexes into m_szChromNames giving the starts of each chromosome name

	CCSVFile *m_pRNAGBSWGSFile;						// load RNA/GBS/WIG sample names and RNA associated metadata using this CSV parser
	CCSVFile *m_pInCntsFile;						// load expression level counts using this CSV parser

	CCSVFile *m_pInWGSvWGSFile;						// load WGS vs. WGS scores using this CSV parser
	CCSVFile *m_pInRNAvWGSFile;						// load RNA vs. WGS scores from this file
	CCSVFile* m_pInRNAvRNAFile;						// load RNA vs. WGS scores from this file

	size_t m_UsedRNACntsMem;						// actually using this much memory in m_pRNACntsMem for holding currently loaded RNA feature counts
	size_t m_AllocdRNACntsMem;						// allocated total size for m_pRNACntsMem
	uint8_t *m_pRNACntsMem;							// allocated to hold RNA feature counts (each count is a int32_t)

	size_t m_UsedRNAhomozScoresMem;					// actually using this much memory in m_pRNAhomozScoresMem for holding currently loaded RNA scores
	size_t m_AllocdRNAhomozScoresMem;				// allocated total size for m_pRNAhomozScoresMem
	uint8_t *m_pRNAhomozScoresMem;					// allocated to hold RNA scores (each score is a double)

	size_t m_UsedRNAvRNAhomozScoresMem;					// actually using this much memory in m_pRNAvRNAhomozScoresMem for holding currently loaded RNAvRNA scores
	size_t m_AllocdRNAvRNAhomozScoresMem;				// allocated total size for m_pRNAvRNAhomozScoresMem
	uint8_t* m_pRNAvRNAhomozScoresMem;					// allocated to hold RNAvRNA scores

	size_t m_UsedWGShomozScoresMem;					// actually using this much memory in m_pWGShomozScoresMem for holding currently loaded WGS scores
	size_t m_AllocdWGShomozScoresMem;				// allocated total size for m_pWGShomozScoresMem
	uint8_t *m_pWGShomozScoresMem;					// allocated to hold WGS scores

	// as generated by LoadRNACntsFile()
	int32_t m_NumbHdrRNAvCntsMappings;					// number of header RNA sample identifiers mappings in m_HdrRNAvCntsMappings
	int32_t m_HdrRNAvCntsMappings[cMaxSampleRows];		// mappings between order of header RNA sample names in RNA file and RNA sample names in materials file

	// as generated by LoadRNAvsWGSScoresFile()
	int32_t m_NumbHdrRNAvWGSMappings;					// number of header RNA sample identifiers mappings in m_HdrRNAvWGSMappings
	int32_t m_HdrRNAvWGSMappings[cMaxSampleRows];		// mappings between order of header RNA sample names in RNA file and RNA sample names in materials file
	int32_t m_NumbRowRNAvWGSMappings;				// number of row RNA sample identifiers mappings in m_RowRNAvWGSMappings
	int32_t m_RowRNAvWGSMappings[cMaxSampleRows];		// mappings between order of row RNA sample names in RNA file and RNA sample names in materials file

	// as generated by LoadWGSvsWGSScoresFile()
	int32_t m_NumbHdrWGSvWGSMappings;					// number of header WGS sample identifier mappings in m_HdrWGSvWGSMappings
	int32_t m_HdrWGSvWGSMappings[cMaxSampleRows];		// mappings between order of header WGS sample names in score files and WGS sample names in materials file
	int32_t m_NumbRowWGSvWGSMappings;					// number of header WGS sample identifier mappings in m_RowWGSvWGSMappings
	int32_t m_RowWGSvWGSMappings[cMaxSampleRows];		// mappings between order of row WGS sample names in score files and WGS sample names in materials file

	// as generated by LoadRNAvsRNAScoresFile()
	int32_t m_NumbHdrRNAvRNAMappings;				// number of header RNAvRNA sample identifier mappings in m_HdrRNAvRNAMappings
	int32_t m_HdrRNAvRNAMappings[cMaxSampleRows];	// mappings between order of header  RNAvRNA sample names in score files and RNA sample names in materials file
	int32_t m_NumbRowRNAvRNAMappings;				// number of row RNA sample identifiers mappings in m_RowRNAvRNAMappings
	int32_t m_RowRNAvRNAMappings[cMaxSampleRows];	// mappings between order of row RNA sample names in RNA vs. RNA file and RNA sample names in materials file




	tsRNAGBSWGSMaterial m_RNAGBSWGSMaterials[cMaxRNAESamples];	// associated RNA/GBS/WGS identifiers and materials, indexed by RNASampleID-1

	void Reset(void);								// reset instance state back to that immediately following instantiation

	uint32_t		// returned chrom identifier, 0 if unable to accept this chromosome name
		AddChromName(char* pszChromName); // associate unique identifier with this chromosome name

	char*	// returns chromosome name associated with 
		LocateChromName(uint32_t ChromID); // chromosome identifier

	uint32_t		// returned chrom identifier, 0 if unable to locate this chromosome name
		LocateChromID(char* pszChrom); // return unique identifier associated with this chromosome name

	uint32_t		// returned RNA metadata name identifier, 0 if unable to accept this metadata name
		AddRNAMetaName(char* pszRNAMetaName); // associate unique identifier with this chromosome name

	char*	// returns RNA metadata name associated with 
		LocateRNAMetaName(uint32_t RNAMetaNameID); // metadata identifier

	uint32_t		// returned RNA metadata name identifier, 0 if unable to locate this metadata name
		LocateRNAMetaNameID(char* pszRNAMetaName); // return unique identifier associated with this metadata name


	uint32_t		// returned material identifier, 0 if unable to accept this material name
		AddMaterialName(char* pszMaterialName); // associate unique identifier with this material name

	char*	// returns material name
		LocateMaterialName(uint32_t MaterialNameID); // material identifier

	uint32_t		// returned material identifier, 0 if unable to locate this material name
		LocateMaterialNameID(char* pszMaterialName); // return unique identifier associated with this material name

	char*	// returns material name corresponding to WGS identifier
		LocateWGSMaterialName(uint32_t WGSNameID); // WGS identifier

	uint32_t		// returned WGS sample name identifier, 0 if unable to accept this WGS name
		AddWGSSampleName(char* pszWGSSampleName); // associate unique identifier with this WGS name

	char*	// returns WGS sample name
		LocateWGSSampleName(uint32_t WGSSampleID); // WGSSample identifier

	uint32_t		// returned WGS sample identifier, 0 if unable to locate this WGS sample name
		LocateWGSSampleNameID(char* pszWGSSampleName); // return unique identifier associated with this WGS name


	uint32_t		// returned SampleName identifier, 0 if unable to accept this SampleName name
		AddRNASampleName(char* pszSampleName);	// associate unique identifier with this SampleName, SampleNames must be unique

	char*							// returned SampleName name
		LocateRNASampleName(uint32_t SampleNameID);	// identifier returned by call to AddSampleName

	uint32_t		// returned SampleName identifier, 0 if unable to locate existing SampleName
		LocateRNASampleNameID(char* pszSampleName); // associate unique identifier with this SampleName


	uint32_t		// returned FeatureName name identifier, 0 if unable to accept this FeatureName
		AddFeatureName(char* pszFeatureName);	// associate unique identifier with this FeatureName, FeatureName names must be unique

	char*							// returned SampleName
		LocateFeatureName(uint32_t SampleNameID);	// identifier returned by call to AddSampleName

	uint32_t		// returned FeatureName identifier, 0 if unable to locate existing FeatureName
		LocateFeatureNameID(char* pszFeatureName); // associate unique identifier with this FeatureName

	int LoadRNACntsFile(char* pszInCntsFile,			// parse and load counts from this file
						bool bNormaliseCnts = true);	// normalise individual RNA feature counts, aka library size, to maximal total counts of any replicate

	int LoadWGSvsWGSScoresFile(char* pszInWGSvWGSFile);		// parse and load WGS vs. WGS homozygosity scores

	int LoadRNAvsWGSScoresFile(char* pszInRNAvWGSFile);		// parse and load RNA vs. WGS homozygosity scores

	int LoadRNAvsRNAScoresFile(char* pszInRNAvWGSFile);		// parse and load RNA vs. WGS homozygosity scores
	
	int LoadMaterialRNAGBSWGSFile(char* pszMaterialRNAGBSWGSFile); // Load sample material mappings

	int	AllocateRNAcnts(int EstNumFeatureCnts);	// initially allocate for this estimated total number of RNA feature counts, memory will be reallocd to hold more if subsequently required

	int ReallocRNAcnts(int NumFeatCntsPerRow);	// realloc as may be required to have room for at least one more row of RNA feature counts

	int AllocateRNAhomozScores(int EstRNAhomozScores);	// initially allocate for this estimated total number of RNA vs. WGS homozygosity scores, memory will be reallocd to hold more if subsequently required

	int ReallocRNAhomozScores(int32_t NumScoresPerRow); // realloc as may be required to have room for at least one more row of RNA homozygosity scores 

	int AllocateWGShomozScores(int EstWGShomozScores);	// initially allocate for this estimated total number of WGS vs. WGS homozygosity scores, memory will be reallocd to hold more if subsequently required

	int ReallocWGShomozScores(int32_t NumScoresPerRow); // realloc as may be required to have room for at least one more row of WGS homozygosity scores

	int AllocateRNAvRNAhomozScores(int EstWGShomozScores);	// initially allocate for this estimated total number of RNAvRNA homozygosity scores, memory will be reallocd to hold more if subsequently required

	int ReallocRNAvRNAhomozScores(int32_t NumScoresPerRow); // realloc as may be required to have room for at least one more row of RNAvRNA homozygosity scores

	double	// returned correlation coefficient, used for finding maximal replicate pairing pearson 'r's
		PearsonCorrelationRNAvCnts(int32_t StartRow,	// starting from this row (1..N)
			int32_t EndRow,								// ending with this inclusive row (N)
			int32_t xSampleID,							// x sample (column)
			int32_t ySampleID);							// y sample (column)

	double	// returned correlation coefficient
		PearsonCorrelationRNAvRNA(int32_t StartRow,		// starting from this row (1..N)
			int32_t EndRow,								// ending with this inclusive row (N)
			int32_t xSampleID,							// x sample (column)
			int32_t ySampleID);							// y sample (column)

	double	// returned correlation coefficient
		PearsonCorrelationRNAvWGS(int32_t StartRow,		// starting from this row (1..N)
			int32_t EndRow,								// ending with this inclusive row (N)
			int32_t xSampleID,							// x sample (column)
			int32_t ySampleID);							// y sample (column)

	double	// returned correlation coefficient
		PearsonCorrelationWGSvWGS(int32_t StartRow,		// starting from this row (1..N)
			int32_t EndRow,								// ending with this inclusive row (N)
			int32_t xSampleID,							// x sample (column)
			int32_t ySampleID);							// y sample (column)

	int	// checking for expression level counts inconsistencies - if replicate 1 and replicate 2 from same biological sample then expecting expression level counts to track over the features
		GenExprCntsPearsons(eModeRNAE PMode,			// processing mode
				  char* pszMaterialRNAGBSWGSFile,		// load RNA/GBS/WIG sample names and RNA associated metadata
				  char* pszInCntsFile,					// load coverage counts from this file
				  char* pszOutRslts);					// write results to this file, will be suffixed appropriately

	int	// checking for RNA replicate homozoygosity inconsistencies - if replicate 1 and replicate 2 from same biological sample then expecting homozygosity to track
		GenRNAvRNAPearsons(eModeRNAE PMode,			// processing mode
			char* pszMaterialRNAGBSWGSFile,		// load RNA/GBS/WIG sample names and RNA associated metadata
			char* pszInHomozygosityFile,		// load homozygosity from this file
			char* pszOutRsltsFile);				// write results to this file, will be suffixed appropriately


public:
	CRNAExpr();	// constructor
	~CRNAExpr();	// destructor

	int
		Process(eModeRNAE PMode,			// processing mode
				char* pszMaterialRNAGBSWGSFile, // load RNA/GBS/WIG sample names and RNA associated metadata
				char* pszInCntsFile,		// load coverage counts from this file
				char* pszInWGSvWGSFile,		// load WGS vs. WGS homozygosity scores from this file
				char* pszInRNAvWGSFile,		// load RNA vs. WGS homozygosity scores from this file
				char* pszInRNAvRNAFile,		// load RNA vs. RNA homozygosity scores from this file
				char* pszOutRsltsFile);			// write results to this file, will be suffixed appropriately
 };

