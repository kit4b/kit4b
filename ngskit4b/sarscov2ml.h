#pragma once

const int cMaxMatrixCols = 30000;		// allowing for matrix containing at most this many loci columns - SARS-CoV-2 has just under 30Kbp
const int cMaxMatrixRows = 50000;		// allowing for matrix to contain this many readsets or isolates
const int cMaxClassifications = 100;	// allowing for this many classifications
const int cMaxR_nCr = 50;				// allowing for at most this many r elements as a combination to be drawn from n total elements
const int cMaxN_nCr = 1000;				// allowing for at most this many n total elements from which r elements can be drawn from as a combination

const uint32_t DfltMinRowsClassified = 50;	// default minimum number of rows (samples) containing feature

const int cMaxAllocOutBuff = 0x0ffffff;	// output buffering for this many chars

typedef enum TAG_eModeSC2 {
	eMSC2default = 0,		// default processing is to locate linkages between features
	eMSC2PlaceHolder		// used to mark end of processing modes
} eModeSC2;

#pragma pack(1)
typedef struct TAG_sScoredCol
	{
	double Prob;				// PValue of feature significance
	uint32_t ColIdx;			// matrix (feature) column index 
	} tsScoredCol;

typedef struct TAG_sTrainClassify {
	uint32_t ReadsetID;			// readset classification
	uint32_t ClassificationID;	// associated classification
} tsTrainClassify;

typedef struct TAG_sClassifications {
	uint32_t ClassificationID;	// classification
	uint32_t szClassificationIdx; // classification name starts at m_szClassificationNames[szClassificationIdx]
	uint32_t NumReadsets;		// number of readsets with this classification
	double Proportion;		// readsets with this classification are this proportion of all classified readsets
} tsClassifications;
#pragma pack()

class CSarsCov2ML
{

	eModeSC2 m_PMode;							// processing mode
	char m_szMatrixFile[_MAX_PATH];				// input matrix file
	char m_szIsolateClassFile[_MAX_PATH];		// input isolate classification file
	char m_szOutFile[_MAX_PATH];				// output linked features
	int m_hOutFile;								// opened output file handle

	int m_OutBuffIdx;							// currently output buffering this many chars
	int m_AllocOutBuff;							// output buffer can hold a max of this many chars
	char *m_pOutBuffer;							// allocated for output buffering


	uint32_t m_NumFeatureNames;						// number of chromosome names currently in m_szChromNames
	uint32_t m_NxtszFeatureIdx;						// current concatenated (names separated by '\0') of all Feature names in m_szFeatureNames 
	char m_szFeatureNames[cMaxMatrixCols * (cMaxDatasetSpeciesChrom/2)];	// used to hold concatenated Feature names, each separated by '\0'
	uint32_t m_szFeatureIdx[cMaxMatrixCols];		// array of indexes into m_szFeatureNames giving the starts of each Feature name

	uint32_t m_NumReadsetNames;						// number of readset names currently in m_szReadsetNames
	uint32_t m_NxtszReadsetIdx;						// current concatenated (names separated by '\0') of all readset names in m_szReadsetNames 
	char m_szReadsetNames[cMaxMatrixRows * (cMaxDatasetSpeciesChrom / 2)];	// used to hold concatenated readset names, each separated by '\0'
	uint32_t m_szReadsetIdx[cMaxMatrixRows];	// array of indexes into m_szReadsetNames giving the starts of each readset name



	uint32_t m_NumClassificationNames;						// number of Classification names currently in m_szClassificationNames
	uint32_t m_NxtszClassificationIdx;						// current concatenated (names separated by '\0') of all Classification names in m_szClassificationNames 
	char m_szClassificationNames[cMaxClassifications * cMaxDatasetSpeciesChrom];	// used to hold concatenated Classification names, each separated by '\0'
	tsClassifications m_Classifications[cMaxClassifications];	// classifications

	uint32_t m_TotClassified;						  // total number of readsets classified
	tsTrainClassify m_ReadsetClassification[cMaxMatrixRows]; // individual readsets and their classifications

	uint32_t m_NumCols;								// matrix has this number of columns including row name identifiers in 1st column
	uint32_t m_NumRows;								// matrix has this number of rows including column name identifiers in 1st row
	uint32_t *m_pMatrix;							// allocated to hold matrix [m_NumRows][m_NumCols]

	uint32_t m_nCombination;						// from n total elements
	uint32_t m_rCombination;						// will be drawing combinations of r elements
	uint32_t m_nCrCombinations[cMaxR_nCr];			// indexes for each potential nCr combination instance

	tsScoredCol m_TopLinkages[cMaxN_nCr];			// columns containing features which potentially have linkages
	uint32_t m_ColClassifiedThresCnts[cMaxClassifications];	// classification thresholds
	uint32_t m_ColClassifiedBelowCnts[cMaxClassifications];

	uint32_t m_MinRowsClassified;					// at least this number of rows (samples) must have been characterised

	bool Init_nCr(int n,							// initialise for n total elements
					int r);							// from which will be drawing combinations of r elements

	bool											// false if all combinations iterated, true if next combination generated
			Iter_nCr(void);							// iterator for drawing next combination of r elements from n total


	int LoadMatrixDimensions(char* pszMatrixFile,	// matrix file to dimension
					uint32_t *pCols,				// number of columns - includes extra column holding row name identifiers
					uint32_t *pRows);				// number of rows - includes extra row holding column name identifiers

	int	LoadMatrix(char* pszMatrixFile); // matrix file to load

	int LoadClassifications(char* pszClassFile);	// classifications file to load

	uint32_t		// returned feature (loci) identifier, 0 if unable to accept this feature name
			AddFeature(char* pszFeature); // associate unique identifier with this feature name
	char* // returned ptr to feature name
			LocateFeature(uint32_t FeatureID);	// feature identifier

	uint32_t LocateFeature(char *pszFeaturePrefix);	// match on this prefix, may be full length

	uint32_t LocateReadset(char *pszReadsetPrefix);	// match on this prefix, may be full length

	uint32_t		// returned readset identifier, 0 if unable to accept this readset name
			AddReadset(char* pszReadset); // associate unique identifier with this readset name
	char* // returned ptr to readset name 
			LocateReadset(uint32_t ReadsetID); // readset name identifier

	uint32_t		// returned classification identifier, 0 if unable to accept this classification name
			AddClassification(char* pszClassification); // associate unique identifier with this classification name
	char* // returned ptr to classification
			LocateClassification(uint32_t ClassificationID);	// Classification identifier

public:
	CSarsCov2ML(void);
	~CSarsCov2ML(void);

	void Reset(void);							// reset state back to that immediately following class instantiation


	int Process(eModeSC2 Mode,						// processing mode
					uint32_t NumLinkedFeatures,			// require this many features to be linked
					uint32_t MinPropRows,				// require at least this many rows to show same linkage
					uint32_t FeatClassValue,			// linkage is between these minimum feature class values
					char* pszMatrixFile,			// input matrix file
					char* pszIsolateClassFile,		// input isolate classification file
					char* pszOutFile);				// output feature classifications file

	int	RunKernel(eModeSC2 Mode,						// processing mode
					uint32_t NumLinkedFeatures,			// require this many features to be linked
					uint32_t MinPropRows,				// require at least this many rows to show same linkage
					uint32_t FeatClassValue);			// linkage is between these minimum feature class values

};

