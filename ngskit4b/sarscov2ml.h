#pragma once

const int cMaxMatrixCols = 30000;		// allowing for matrix containing at most this many loci columns - SARS-CoV-2 has just under 30Kbp
const int cMaxMatrixRows = 50000;		// allowing for matrix to contain this many readsets or isolates
const int cMaxClassifications = 100;	// allowing for this many classifications

typedef enum TAG_eModeSC2 {
	eMSC2default = 0,		// no modes, currently just a default!
	eMSC2PlaceHolder			// used to mark end of processing modes
} eModeSC2;

#pragma pack(1)
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


	int Process(eModeSC2 Mode,					// processing mode
			char* pszMatrixFile,				// input matrix file
			char* pszIsolateClassFile,			// input isolate classification file
			char* pszOutFile);					// output feature classifications file

	int	RunKernel(uint32_t NumLinkedFeatures,	// require this many features to be linked
				  uint32_t MinPropRows,			// require at least this many rows to show same linkage
				  uint32_t FeatType);			// linkage is between these minimum feature types

};

