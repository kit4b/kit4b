#pragma once

const uint32_t cMaxChromMappings = 10000; // maximum of this many chromosome size and name mappings
const uint32_t cMaxChromNames = (cMaxChromMappings * 2);	// can process at most this many chromosome mappings
const uint32_t cMaxAllocRsltsBuff = 0x0ffffff;	// allocate this sized buffer for outputting results
const int32_t cMaxGBSF4s = 2000;			// max number of F4s expected to be present in GBS input file

typedef enum TAG_eModeGBSMapSNPs
{
eMGBSMdefault = 0, // default is to map from SNP call loci into PBA windowed bins

} eModeGBSMapSNPs;

#pragma pack(1)
typedef struct TAG_sChromMapping
{
	uint32_t RefChromID;	// unique identifier
	uint32_t AliasChromID;	// alias for reference chrom
	uint32_t Size;			// chromosome size (bp)
	uint32_t StartBinID;	// starting windowed bin identifier
	uint32_t EndBinID;		// final bin identifier (inclusive)
} tsChromMapping;

typedef struct TAG_sWinBin
{
uint32_t BinID;				// uniquely identifies this bin
uint32_t RefChromID;		// bin is on this chromosome
uint32_t StartLoci;			// bin starts at this loci
uint32_t EndLoci;			// bin ends at this loci inclusive
uint32_t LociCounts;		// has this many SNP loci
uint32_t CntNA;				// count of loci having F4 SNPs which were either missing (NA) or didn't match either Fa or Fb
uint32_t CntExclFa;			// count of loci having F4 SNPs exclusively matching Fa
uint32_t CntExclFb;			// count of loci having F4 SNPs exclusively matching Fb
uint32_t CntFaFb;			// count of loci having F4 SNPs matching both Fa and Fb
} tsWinBin;

#pragma pack()

class CGBSmapSNPs {
	eModeGBSMapSNPs m_PMode;	// processing mode
	int m_BinSize;			// SNP loci counts are accumulated into this sized (bp) non-overlapping bins
	char *m_pszGBSName;		// column name referencing the GBS readset which is to be processed
	char *m_pszInNMFile;	// processing chromosome name mapping from this file
	char *m_pszInGBSFile;	// processing GBS SNP calls from this file
	char *m_pszOutFile;		// windowed GBS SNP calls output file

	uint32_t m_WIGChromID;				// current WIG span is on this chromosome
	uint32_t m_WIGRptdChromID;			// previously reported WIG chrom
	uint32_t m_WIGSpanLoci;				// current WIG span starts at this loci
	uint32_t m_WIGSpanLen;				// current span is this length
	uint32_t m_WIGRptdSpanLen;			// previously reported WIG span length
	uint64_t m_WIGSpanCnts;				// current span has accumulated this many counts

	char m_szFndrA[cMaxDatasetSpeciesChrom];		// founder A name
	char m_szFndrB[cMaxDatasetSpeciesChrom];		// founder B name

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

	uint32_t m_NumWinBins;		// total number of window bins used
	tsWinBin *m_pWinBins;		// single allocation to hold m_NumWinBins windowed bins

	uint32_t m_OutBuffIdx;		// currently buffering this many output bytes
	uint32_t m_AllocOutBuff;		// m_pOutBuffer allocated to hold this many output bytes
	uint8_t *m_pszOutBuffer;		// allocated for buffering output

	uint32_t					// returned chrom identifier, 0 if unable to accept this chromosome name
		AddChrom(char* pszChrom); // associate unique identifier with this chromosome name

	char*	// returns chromosome name associated with 
		LocateChrom(uint32_t ChromID); // chromosome identifier

	uint32_t		// returned chrom identifier, 0 if unable to locate this chromosome name
		LocateChrom(char* pszChrom); // return unique identifier associated with this chromosome name

	tsChromMapping *
		LocateAliasChromMapping(char* pszChrom);	// locates and returns chromosome mapping for a alias name which matches the prefix of this chromosome 

	int	LoadNM(char* pszInNMFile);		// processing chromosome name mapping from this file);
	int LoadGBSSNPs(char* pszGBSName,	// column name referencing the GBS readset which is to be processed
					char* pszInGBSFile);	// processing chromosome name mapping from this file)

	int ReportBinCounts(char* pszOutFile, // reporting bin counts to this file as CSV
						char* pszGBSName);				// counts are for this GBS readset

	int ReportCountsWIG(char* pszOutFile,
							 bool bFb,						// false if founder Fa, true if founder Fb
							 char* pszGBSName,				// column name referencing the GBS readset which is to be processed
							 char *pszFounder);				// founder

	void InitialiseWIGSpan(void);				 // initialise WIG span vars to values corresponding to no spans having been previously reported
	int CompleteWIGSpan(bool bWrite = false);	// close off any current WIG span ready to start any subsequent span, if bWrite is true then write to disk
	int AccumWIGBinCnts(uint32_t ChromID,	// accumulate wiggle counts into variableStep spans over this chromosome
			uint32_t Loci,		// bin counts starting from this loci - WIG uses 1 based loci instead of the usual UCSC BED 0 based loci
			uint32_t Cnts,		// bin has this many counts attributed
			uint32_t BinLen);	// bin is this length
	int AccumWIGCnts(uint32_t ChromID,	// accumulate wiggle counts into variableStep spans over this chromosome
			uint32_t Loci,		// this loci  - WIG uses 1 based loci instead of the usual UCSC BED 0 based loci
			uint32_t Cnts,		 // has this many counts attributed
			uint32_t MaxSpanLen = 1000000); // allow WIG spans to be this maximal length

	void Reset(void);			// reset class states back to that immediately following class instantiation 
public:
	CGBSmapSNPs();
	~CGBSmapSNPs();

	int			// success or otherwise ( >= 0 success, < 0 if processing failed)
	Process(eModeGBSMapSNPs PMode,			// processing mode
							int BinSize,					// SNP loci counts are accumulated into this sized (Kbp) non-overlapping bins
							char *pszGBSName,				// column name referencing the GBS readset which is to be processed
							char *pszInNMFile,		// processing chromosome name mapping from this file
							char *pszInGBSFile,	// processing GBS SNP calls from this file
							char *pszOutFile);		// windowed GBS SNP calls output file

};

