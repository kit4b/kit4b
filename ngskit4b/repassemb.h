#pragma once
// repurpose assembly using low coverage SNPs from a subspecies to
// modify a reference assembly
// Useful when attempting to discover differential haplotype regions
// in F4's relative to unassembled founders
//
const int cMaxSNPChromNames = 5000000;		// can only handle a max of this many chromosome names - memory is statically allocated
const int cAllocNumSNPs = 10000000;			// allocate for SNPs in this sized blocks
const int cAllocNumASegments = 10000000;	// allocate for aligned segments in this sized blocks
const int cAlignNumSSNPfields = 23;			// if generated by 'kit4b kalign' then there will be this many CSV fields
const int cAlignNumSSNPXfields = 23 + (3 * 65);	// if generated by 'kit4b kalign' for all 3 read counts frame shifts (3*4^3) + fixed reference frame shifts (3*3) then there will be this many CSV fields
const int cSSNPMarkerNumFields = 4 + 9;		// if generated by 'kit4b snpmarkers' then there will be a minimum of this number of CSVfields
const int cAllocAssembFastaMem = 0x03fffffff;		// allocating memory buffering for assembly sequences in chunks of this size 

const uint32_t cDfltMinSiteCoverage = 20;		// must be at least this number of covering bases at SNP site to be accepted
const double cDfltMinAlleleProp = 0.5;			// major allele must be at least this proportion of all counts (incl ref base) at SNP site to be accepted
const double cDfltPValueThres = 0.05;			// and PValue must be <= this threshold

#pragma pack(1)
typedef enum TAG_eRPPMode{
	eRPMdefault,		// default processing mode is for SNP major allele replacements
	eRPMNonalignSegs,	// processing for non-aligned segments - these are assumed to be deletions from target assembly
	eRPMplaceholder	// placeholder
} etRPPMode;


typedef struct TAG_sASegment // processing for non-aligned segments - these are assumed to be deletions from target assembly
{
uint32_t ChromID;		// aligned segment is on this chromosome
uint32_t Loci;			// starting at this loci inclusive
uint32_t Len;			// and is this length (number bases aligned)
} tsASegment;

typedef struct TAG_sSNPSSite // default processing mode is for SNP major allele replacements
{
	uint32_t ChromID;		// SNP is on the chromosome
	uint32_t Loci;			// at this loci
	etSeqBase RefBase;		// reference assembly base
	etSeqBase RepBase;		// use this base when repurposing the assembly - if eBaseN then do not repurpose
	uint32_t AlleleBaseCnts[5]; // counts of major allelic base called over all readset alignments
} tsSNPSite;

typedef struct TAG_sChromSites
{
	uint32_t ChromID;		// SNPs are on this chromosome
	uint32_t SNPSiteIdx;	// index into m_pSNPSites for 1st site on this chromosome
	uint32_t NumSNPSites;	// number of SNP sites on this chromosome
} tsChromSites;

#pragma pack()

class CRepAssemb
{
	char m_szTargAssemblyName[cMaxDatasetSpeciesChrom + 1]; // alignments were against this targeted assembly
	int m_NumSNPFiles;								// re-purposing from this number of kalign generated SNP files
	char m_szCurChrom[cMaxDatasetSpeciesChrom+1];	// processing SNP on this chromosome
	uint32_t m_CurChromID;						// current chromosome identifier

	int m_LAChromNameID;						// last accessed chromosome identifier from call to AddChrom()
	int m_NumChromNames;						// number of chromosome names currently in m_szChromNames
	int m_NxtszChromIdx;						// current concatenated (names separated by '\0') of all chromosome names in m_szChromNames
	int m_NumChromSites;						// number of chromosomes with at least one accepted SNP site
	char m_szChromNames[cMaxSNPChromNames * (cMaxDatasetSpeciesChrom/2)];	// used to hold concatenated chromosome names, each separated by '\0'
	uint32_t m_szChromIdx[cMaxSNPChromNames];		// array of indexes into m_szChromNames giving the starts of each chromosome name
	tsChromSites m_ChromSites[cMaxSNPChromNames];	// chromosome sites

	CCSVFile* m_pCSV;							// used to load SNP calls

	uint32_t m_MinCoverage;						// must be at least this number of covering bases at SNP site to be accepted
	double m_MinAlleleProp;						// major allele must be at least this proportion of all alleles (incl ref base) at SNP site to be accepted
	double m_PValueThres;						// and PValue must be <= this threshold

	size_t m_SeqBuffIdx;						// currently buffering this many bases
	size_t m_AllocSeqBuffMem;					// m_pSeqBuffer allocated to this size (realloc'd if required)
	uint8_t *m_pSeqBuffer;						// allocated to buffer input/output fasta sequences 
	int m_hInFile;								// input assembly file handle
	int m_hOutFile;								// output assembly file handle
	int m_hOutRepSNPsFile;						// output consensus sites file handle

	uint32_t m_NumSortedSites;					// number of sorted sites enabling a binary search over this number of sites
	uint32_t m_NumSNPSites;						// number of SNP sites accepted
	uint32_t m_NumAllocdSNPSites;				// number of SNP sites currently allocated
	size_t m_AllocdSNPSitesMem;					// allocd memory size for holding m_NumAllocdSNPSites
	tsSNPSite *m_pSNPSites;						// pts to allocd SNP sites

	uint32_t m_MaxASegmentLen;					// max length aligned segment in m_pASegments
	uint32_t m_NumASegments;					// number of aligned segments accepted
	uint32_t m_NumAllocdASegments;				// number of aligned segments currently allocated
	size_t m_AllocdASegmentsMem;				// allocd memory size for holding m_NumAllocdASegments
	tsASegment *m_pASegments;					// pts to allocated aligned segments
					

	int		// returned chrom identifier, < 1 if unable to accept this chromosome name
		AddChrom(char* pszChrom); // associate unique identifier with this chromosome name

	char*	// returns chromosome name associated with 
		LocateChrom(int ChromID); // chromosome identifier

	int		// returned chrom identifier, < 1 if unable to locate this chromosome name
		LocateChrom(char* pszChrom); // return unique identifier associated with this chromosome name

	int MergeASegs(char* pszChrom,		// segment is on this chromosome
			uint32_t Loci,				// starting at this loci inclusive
			uint32_t Len,				// segment is this length
			uint32_t Score = 1);		// 0 if segment is an unaligned segment

	tsChromSites * // returned ptr
		LocateChromSites(char* pszChrom);	// sites are to be those on this chromosome

	int						// returned total number of SNPs parsed and accepted
		LoadKalignSNPs(char *pszSNPFile);	// load and parse kalign generated SNPs from wild-carded file specification

	int						// eBSFSuccess or error
		LoadKalignASegs(char *pszASegsFiles);	// load and parse kalign generated aligned segments from this - can be wild-carded


	int
		LoadAssembly(char* pszAssembly); // load assembly sequences in this fasta file into memory (m_pSeqBuffer will be allocated to hold raw fasta sequences)

	int									// number of bases replaced
		ReplaceSNPBases(char *pszChrom,	// chromosome name
			 uint32_t SeqLen,			// max of this number of bases in sequence
			 char *pszSeq);				// sequence

	int									// number of bases replaced
		ReplaceASegmentBases(char *pszChrom, // chromosome name
			 uint32_t SeqLen,			// max of this number of bases in fasta sequence length - includes any whitespace
			 char *pszSeq);				// sequence

	int ReportSNPdist(char *pszRepAssembFile);	// base file name used for reporting major allele SNPs used when purposing assembly file

	CStats m_Stats;								// basic statistics

	tsSNPSite *LocateSite(uint32_t ChromID, uint32_t Loci);	// locate existing site using a binary search over sorted sites
	tsASegment *LocateASegment(uint32_t ChromID, uint32_t Loci);	// locate a tsASegment containing the requested ChromID.Loci
	static int SortSNPChromLoci(const void* arg1, const void* arg2);
	static int SortASegmentChromIDLoci(const void* arg1, const void* arg2);

public:
	CRepAssemb(void);			// constructor

	~CRepAssemb(void);			// destructor
	
	void Reset(void);				// resets instance state back to that immediately following instantiation

	int ProccessAssembly(etRPPMode PMode,	// processing mode
				uint32_t MinCoverage,			// must be at least this number of covering bases at SNP site to be accepted
				double MinAlleleProp,			// major allele must be at least this proportion of all alleles (incl ref base) at SNP site to be accepted
				double PValueThres,				// SNP PValue must be <= this threshold
				char *pszProcFile,		// input kalign generated SNPs or alignment segments file name, filename can be wildcarded
				char *pszAssembFile,	// input assembly to repurpose
				char *pszRepAssembFile);	// write out repurposed assembly to this file

};

