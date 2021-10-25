#pragma once

const uint32_t cMaxChromMappings = 10000; // maximum of this many chromosome size and name mappings
const uint32_t cMaxChromNames = (cMaxChromMappings * 2);	// can process at most this many chromosome mappings
const uint32_t cOutBuffSize = 0x0ffffff;	// allocate this sized buffer for outputting results
const int32_t cMaxFounderReadsets = 3;				// max number of founders
const int32_t cMaxProgenyReadsets = 1000;			// max number of F4s expected to be present in GBS input file

const uint32_t cAllocProgenyFndrAligns = 1000000; // initial allocation for progeny to founder allele stacks alignments

typedef enum TAG_eModeGBSMapSNPs
{
	eMGBSMDefault = 0, // default is to map SNP GBS to PBA GBS haplotypes
	eMGBSMConsistency, // Compare 2 PBA GBS matrices for haplotype consistency
	eMGBSMPlaceholder  // used to mark end of enumeration range
} eModeGBSMapSNPs;

#pragma pack(1)
typedef struct TAG_s256Bits {
	uint8_t NumBits;		// total number of bits over all 64bit words in Bits256
	uint64_t Bits[4];	// 256bits encoded into 4 x 64bit words
} ts256Bits;

typedef struct TAG_sProgenyFndrAligns {
	uint32_t ReadsetID;				// identifies the progeny readset
	uint32_t ChromID;				// stack is on this chromosome
	uint32_t Loci;					// and at this loci
	uint32_t Source;				// identifies source - 0: conversion of SNP to PBA haplotypes, 1: matrix A, 2: matrix B when comparing 2 matrices for haplotype consistency
	uint8_t Alleles;				// progeny has these alleles present at the allelestack chrom.loci
	uint8_t NumProgenyFounders;		// number of potential progeny founders in ProgenyParents bitmap
	ts256Bits ProgenyFounders;		// bitmap of potential progeny founders - progeny has a minor/major allele shared with corresponding founder
} tsProgenyFndrAligns;

typedef struct TAG_sChromMapping
{
	uint32_t RefChromID;	// unique identifier
	uint32_t AliasChromID;	// alias for reference chrom
	uint32_t Size;			// chromosome size (bp)
} tsChromMapping;


#pragma pack()

class CGBSmapSNPs {
	eModeGBSMapSNPs m_PMode;	// processing mode
	char *m_pszInNMFile;	// processing chromosome name mapping from this file
	char *m_pszInGBSFile;	// processing GBS SNP calls from this file
	char *m_pszOutFile;		// windowed GBS SNP calls output file

	int32_t m_LAReadsetNameID;			// name identifier last returned by AddReadsetName()
	int32_t m_NumReadsetNames;			// number of readsets names currently in m_szReadsets
	uint32_t m_NxtszReadsetIdx;			// current concatenated (names separated by '\0') of all readset names in m_szReadsets, note that the readset type (founder 0, progeny 1 or control 2) is prepended to each name
	char m_szReadsetNames[(cMaxProgenyReadsets + cMaxFounderReadsets) * (cMaxDatasetSpeciesChrom + 1)];	// used to hold concatenated readset names, each separated by '\0', allowing for 1 progeny readset
	uint32_t m_szReadsetIdx[cMaxProgenyReadsets + cMaxFounderReadsets];	// array of indexes into m_szReadsetNames giving the starts of each founder or progeny  name. , allowing for a control readset name
	int32_t m_NumFounders;					// number of founders
	uint8_t m_Fndrs2Proc[cMaxFounderReadsets];	// array of founders which are to be processed, indexed by FounderID-1. If LSB is set then that founder is marked for processing
	int32_t m_FndrIDs[cMaxFounderReadsets];	// founder readset identifiers for Fa..Fn

	int32_t m_NumProgenies;				// number of F4s
	int32_t	m_ProgenyIDs[cMaxProgenyReadsets];		// F4 readset identifiers

	uint32_t m_WIGChromID;				// current WIG span is on this chromosome
	uint32_t m_WIGRptdChromID;			// previously reported WIG chrom
	uint32_t m_WIGSpanLoci;				// current WIG span starts at this loci
	uint32_t m_WIGSpanLen;				// current span is this length
	uint32_t m_WIGRptdSpanLen;			// previously reported WIG span length
	uint64_t m_WIGSpanCnts;				// current span has accumulated this many counts

	uint32_t m_LAChromNameID;	// last accessed chromosome identifier from call to AddChrom()
	uint32_t m_NumChromNames;	// number of chromosome names currently in m_szChromNames
	uint32_t m_NxtszChromIdx;	// current concatenated (names separated by '\0') of all chromosome names in m_szChromNames
	char m_szChromNames[cMaxChromNames * cMaxDatasetSpeciesChrom];	// used to hold concatenated chromosome names, each separated by '\0'
	uint32_t m_szChromIdx[cMaxChromNames];	// array of indexes into m_szChromNames giving the starts of each chromosome name

	CCSVFile *m_pInNMFile;		// to contain chromosome name mapping and chromosome sizes
	CCSVFile *m_pInGBSFile;		// processing GBS SNP calls from this file
	int m_hOutFile;			// file handle for writing results in CSV format
	
	uint32_t m_NumChromMappings;		// number of chromosome name mappings loaded
	tsChromMapping *m_pChromMappings;	// fixed allocation (cMaxChromMappings) to hold chromosome size and name mappings

	uint32_t m_UsedProgenyFndrAligns;			// number of actually used progeny to founder alignments
	uint32_t m_AllocdProgenyFndrAligns;			// number of allocated founder alignments
	size_t m_AllocdProgenyFndrAlignsMem;		// current mem allocation size for m_pProgenyFndrAligns
	tsProgenyFndrAligns* m_pProgenyFndrAligns;	// allocated to hold progeny allele stack alignments

	uint32_t m_OutBuffIdx;		// currently buffering this many output bytes
	uint32_t m_AllocOutBuff;		// m_pOutBuffer allocated to hold this many output bytes
	uint8_t *m_pszOutBuffer;		// allocated for buffering output

	// inplace modification of sample identifier ensuring that sample identifiers have a prefix of 'S0' where original identifier starts with 'S[1-9]'
	 int32_t		// updated strlen(pszSampleID)
		MakeConsistentSampleName(char* pszSampleID);

	int32_t		// returned readset name identifier, < 1 if unable to accept this readset name
		AddReadset(char* pszReadset,		// associate unique identifier with this readset name, readset names must be unique within the readset type
				   uint8_t ReadsetType);	// 0: founder, 1: progeny, 2: control

	char*							// returned readset name
		LocateReadset(int32_t ReadsetID);	// identifier returned by call to AddReadsetName

	int32_t		// returned readset name identifier, 0 if unable to locate existing readset
		LocateReadset(char* pszReadset, // associate unique identifier with this readset name
					  uint8_t ReadsetType);	// 0: founder, 1: progeny, 2: control

	uint32_t					// returned chrom identifier, 0 if unable to accept this chromosome name
		AddChrom(char* pszChrom); // associate unique identifier with this chromosome name

	char*	// returns chromosome name associated with 
		LocateChrom(uint32_t ChromID); // chromosome identifier

	uint32_t		// returned chrom identifier, 0 if unable to locate this chromosome name
		LocateChrom(char* pszChrom); // return unique identifier associated with this chromosome name

	tsChromMapping *
		LocateAliasChromMapping(char* pszChrom);	// locates and returns chromosome mapping for a alias name which matches the prefix of this chromosome 

	int	LoadNM(char* pszInNMFile);		// processing chromosome name mapping from this file);
	int LoadGBSSNPs(uint32_t MatrixID,			// identifies matrix source
					char* pszInGBSFile);	// processing chromosome name mapping from this file)

	int	LoadHaplotypeMatrix(uint32_t MatrixID,			// identifies matrix source
				char* pszInGBSFile);						// processing chromosome name mapping from this file)

	void Reset(void);			// reset class states back to that immediately following class instantiation 

	uint32_t									// returned index+1 into  m_pProgenyFndrAligns[] to allocated and initialised ProgenyFndrAligns, 0 if errors
		AddProgenyFndrAligns(tsProgenyFndrAligns* pInitProgenyFndrAligns);	// allocated tsProgenyFndrAligns to be initialised with a copy of pInitProgenyFndrAligns

	uint8_t			// returned PBA alleles
		SNPs2Alleles(char* pszSNPs,	// translate char representation of major/minor SNPs into it's packed byte allelic representation
				bool bMajorOnly = false);	// true if only major/major single allele to be returned as PBA

	int	ReportMatrix(char* pszRsltsFileBaseName);		// matrix results are written to this file base name with '.SummaryMatrix.csv' appended

	int ReportConsistency(char* pszRsltsFileBaseName);		// results are written to this file base name with '.consistency.csv' appended);
	
	int ReportHaplotypesByProgeny(char* pszRsltsFileBaseName,		// haplotype results are written to this file base name with 'ProgenyID.csv' appended
								  uint32_t ReadsetID = 0);				// report on this progeny readset only, or if 0 then report on all progeny readsets

	void	Bits256Shl1(ts256Bits& Bits256);
	void	Bits256Set(uint8_t Bit,		// bit to set, range 0..255
					   ts256Bits& Bits256);
	void	Bits256Reset(uint8_t Bit,		// bit to reset, range 0..255
						 ts256Bits& Bits256);
	bool Bits256Test(uint8_t Bit,		// bit to test, range 0..255
					 ts256Bits& Bits256);
	void Bits256Initialise(bool Set,			// if true then initialse all bits as set, otherwise initialise all bits as reset
						   ts256Bits& Bits256);
	uint32_t Bits256Count(ts256Bits& Bits256);		// count number of set bits

	uint32_t Bits256Combine(ts256Bits& Bits256A, ts256Bits& Bits256B);		// combine (effective |= ) bits in Bits256A with Bits256B with Bits256A updated 

	CMTqsort m_mtqsort;				// multi-threaded qsort
	static int SortProgenyFndrAligns(const void* arg1, const void* arg2);
	static int SortProgenyChromLociReadset(const void* arg1, const void* arg2);

public:
	CGBSmapSNPs();
	~CGBSmapSNPs();

	int			// success or otherwise ( >= 0 success, < 0 if processing failed)
	Process(eModeGBSMapSNPs PMode,			// processing mode
							char *pszInNMFile,		// processing chromosome name mapping from this file
							char *pszInGBSFile,	// processing GBS SNP calls from this file
							char *pszOutFile);	// windowed GBS haplotype matrix output file

};

