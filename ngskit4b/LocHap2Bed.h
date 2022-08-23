#pragma once

const int cLHMaxSetMembers = 1000;	// a max of this many members in a set
const int cCSVLHMaxFields = 1000;	// max of this number of fields in a csv snp/di/tri file
const int cLHAllocLineBuffSize = 0x0ffffff;	// allocate line buffer for output which can hold this many chars
const double cDfltLocalSeqErrRate = 0.02;	// default localised sequence error rate
const int cMaxREIncludeFiles = cLHMaxSetMembers;	// upto this many regular expressions for files which can be processed
const int cMaxREExcludeFiles = cLHMaxSetMembers;	// upto this many regular expressions for files which can be processed


typedef enum TAG_eModeLocHap {
	eMLHdefault = 0,		// no modes, currently just a default!
	eMLHPlaceHolder			// used to mark end of processing modes
} eModeLocHap;

typedef enum TAG_eRHTFormat {			// report format
	eRHTBed = 0,						// report in BED format
	eRHTCsv,							// report in CSV format
	eRHTPlaceHolder						// currently a placeholder
} eRHTFormat;

typedef enum TAG_etLHSetOp {
	eLHSONone = 0,					// no Set operation
	eLHSOUnion,						// union of SetA and SetB
	eLHSOIntersect,					// intersect of SetA and SetB
	eLHSOCplUnion,					// complement of union SetA and SetB
	eLHSOCplIntersect,				// complement of intersect SetA and SetB
	eLHSOSubtract,					// SetA subtract SetB
	eLHSOPlaceHolder					// just a place holder to mark end of enumeration
} etLHSetOp;

#pragma pack(1)
typedef struct TAG_sCultivar {
	char szName[cMaxDatasetSpeciesChrom + 1];	// name of this cultivar
	etSeqBase CalledBase;						// base called as being the major allele for this cultivar
	int Score;									// score
	int TotalBaseCnts;							// total number of base counts
	int BaseCnts[5];							// counts for each base
	uint8_t flgSNPInferenced : 1;				// 1 if allele was inferenced 
	uint8_t flgInSetA : 1;						// 1 if this cultivar is a member of SetA
	uint8_t flgInSetB : 1;						// 1 if this cultivar is a member of SetB
} tsCultivar;

typedef struct TAG_sLHCnts {
	bool bIsSNP;		// true if at least one of the alleles at this loci is classified as being a SNP
	int Loci;			// loci in haplotype
	int CovBases;		// covering bases at loci
	etSeqBase Base;		// reference base at loci
	int Cnts[5];		// number of counts for each base
} tsLHCnts;

typedef struct TAG_sLocalisedHaplotype {
	char szChrom[cMaxGeneNameLen+1];	// haplotype is on this chromosome
	int NumHaplotypes;	// number of haplotypes which were called by kalign
	int DepthCnt;		// number of reads containing the localised haplotypes
	int AntisenseCnt;	// number of antisense reads containing localised haplotypes
	tsLHCnts LHCnts[3];	// counts for each haplotype loci in this localised haplotype
	int NumAllelicTypes;	// number of accepted allelic combinations
	int SumAlleleCombCnts;	// sum of all AlleleCombCnts; 
	int AlleleCombCnts[64];	// counts for each combination of haplotype alleles
} tsLocalisedHaplotype;
#pragma pack()

class CLocHap2Bed
{
	eRHTFormat m_ReportFormat;						// processing results are to be reported in this file format
	int m_MinCoverage;							// must be at least this read coverage at SNP site
	double m_MinAlleleProp;						// putative allele must be at least this proportion of total site read coverage
	double m_PValueThres;						// acceptance PValue
	double m_LocalSeqErrRate;					// local sequencing error rate
	int m_NumSNPFiles;							// number of input SNP files
	char** m_ppszSNPFiles;						// input SNP files
	char m_szOutBEDFile[_MAX_PATH];				// output in UCSC BED format to this file
	char m_TargAssemblyName[cMaxDatasetSpeciesChrom + 1]; // alignments were against this targeted assembly
	char m_SpecAssemblyName[cMaxDatasetSpeciesChrom + 1]; // in generated pgSNP file use this as the assembly name
	char m_szDescription[_MAX_PATH];
	char m_szTrackName[100];
	int m_NumCultivars;							// number of cultivars being processed
	tsCultivar m_Cultivars[cMaxCultivars];		// cultivar specific metadata
	CCSVFile *m_pCSV;							// used to load SNP calls
	int m_hOutFile;								// output eRHTFormat'd file handle

	int m_LineBuffOffs;							// offset in m_pszLineBuff at which to next write buffered lines
	int m_AllocdLineBuff;						// 	m_pszLineBuff allocated to a max of this many chars			
	char *m_pszLineBuff;						// line buffer memory to hold buffered output

	CUtility m_RegExprs;            // regular expression processing

	int
		Process2BED(void);						// processing to BED format
	int
		Process2CSV(void);						// processing to CSV format - currently unimplemented

	bool					// true if file is accepted for processing, false if file not accepted
		AcceptThisFile(char* pszFileName);

public:
	CLocHap2Bed(void);
	~CLocHap2Bed (void);

	void Reset (void);							// re-initialise

	int Process (eModeLocHap Mode,				// processing mode
				 int MinCoverage,				// must be at least this coverage at SNP site
				 double MinAlleleProp,			// putative allele must be at least this proportion of total site read coverage
				 double PValueThres,			// only accept SNP alleles which have a PValue <= this threshold
				 char *pszTrackName,			// track name
				 char *pszAssemblyName,			// UCSC assembly name - for SARS-CoV-2 it is "wuhCor1"
				 char *pszExperimentDescr,		// describes experiment
				 int NumInWhitelist,			// number of white listed file names
				 char **ppszWhitelisted,		// names of white listed files
				 int NumInBlackList,			// number of black listed files
				 char **ppszBlacklisted,		// names of black listed files
				 int NumSNPFiles,				// number of input SNP files
				 char** ppszSNPFiles,			// input SNP files
				 char *pszOutFile);				// output SNPs to this UCSC Personal Genome SNP format file
};

