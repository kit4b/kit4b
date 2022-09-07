#pragma once

const int32_t cMaxGMLDWorkerThreads = 30;       // allowing at most this number of worker threads
const int32_t cMaxSampleRefNameLen = 100;        // sample reference names can be at most this length
const int32_t cMaxFeatRefNameLen = 200;       // needs to be quite long as may be a concatenation of chrom + loci + species/cultivar

typedef enum TAG_eModeGMLD {
	eGMLDDefault = 0,  // transposition, columns to be features, rows to be samples
	eGMLDPlaceHolder   // acts as a placeholder, sets max-1 range of all enumerations
	} teModeGMLD;


typedef enum TAG_eTypeGMLD {
	eTGMLDDefault = 0,  // default input file format type
	eTGMLDPlaceHolder   // acts as a placeholder, sets max-1 range of all enumerations
	} teTypeGMLD;

typedef enum TAG_eReductGMLD {
	eRGMLDDefault = 0,  // default feature reduction
	eRGMLDPlaceHolder   // acts as a placeholder, sets max-1 range of all enumerations
	} teReductGMLD;

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
#pragma pack()

class CGenMLdatasets {
	teModeGMLD m_PMode;  // processing mode
	teTypeGMLD m_FType;          // input sample feature file format type
	teReductGMLD m_RMode;          // feature reduction mode


	int32_t m_NumSampleRefs; // number of sample reference names contained in m_pszSampleRefs
	size_t m_UsedSampleRefsMem; // actually used this allocation memory 
	size_t m_AllocdSampleRefsMem; // current allocation size for sample references
	uint8_t* m_pszSampleRefsMem;  // allocated to hold concatenated sample reference names, each separated by '\0'

	int32_t m_NumFeatRefs; // number of feature reference names contained in m_pszFeatureRefs
	size_t m_UsedFeatRefsMem; // actually used allocation size
	size_t m_AllocdFeatRefsMem; // current allocation size for feature references
	uint8_t* m_pszFeatRefsMem;  // allocated to hold concatenated feature reference names, each separated by '\0'


	int32_t m_SampleFeatsSize;    // each sample feature requires this much memory
	size_t m_UsedSampleFeatsMem; // actually using this much memory in m_pSampleFeats for holding currently loaded sample features
	size_t m_AllocdSampleFeatsMem; // allocated total size for m_pSampleFeats
	uint8_t *m_pSampleFeatsMem; // allocated to hold sample features

	CCSVFile *m_pSFFile; // input CSV file containing sample features

	int32_t ReduceSampleFeatures(teReductGMLD RMode);       // feature reduction mode
	int32_t LoadSampleFeatures(teTypeGMLD FType,            // input sample feature file format type
						   char* pszInSampleFeats);     // load sample features from this file
	int32_t AssociateSampleLabels(char* pszInSampleLabels); // associate samples with these labels

	int32_t AllocateMemSamples(int NumSamples,  // initially allocate for this number of samples, memory will be reallocd to hold more if subsequently required
					   int NumFeatures); // initially allocate for this number of features
	int32_t ReallocSamples(void);      // realloc as may be required to have room for at least one more sample reference name and associated sample features

	int32_t AddSampleRefName(char *pszSampleRef);           // add sample reference name to m_pszSampleRefsMem and returns the ID 1..N
	int32_t AddFeatRefName(char* pszFeatRef);               // add feature reference name to m_pszFeatsRefsMem and returns the ID 1..N
	int32_t AddSampleFeatValue(int32_t SampleID,  // sample identifier
				   int32_t FeatID,    // feature identifier
				   int32_t Value);    // SampleID.FeatureID has this value
public:
	CGenMLdatasets();       // constructor
	~CGenMLdatasets();      // destructor
	void Reset(void);       // reset back to instantiation

	int Process(teModeGMLD PMode,                  // processing mode
				teTypeGMLD FType,                  // input sample feature file format type
				teReductGMLD RMode,                  // feature reduction mode
				char* pszInSampleFeats,     // input sample feature file
				char* pszInSampleLabels,    // input sample labels (classes) file
				char* pszOutMLdataset,     // output ML dataset - row per sample, column per feature, last column containing sample labels
				int32_t	NumThreads);		// maximum number of worker threads to use

	};

