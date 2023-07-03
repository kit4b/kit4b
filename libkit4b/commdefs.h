#pragma once

#ifndef _WIN32
//#define min(a,b) ((a) < (b) ? (a) : (b))
//#define max(a,b) ((a) > (b) ? (a) : (b))
#endif

#ifdef _WIN32
// flags when opening existing file for read/write with mainly sequential access
#define O_READORWRITESEQ  (_O_RDWR | _O_BINARY | _O_SEQUENTIAL)

// flags when opening existing file for read with mainly sequential access
#define O_READSEQ  (_O_RDONLY | _O_BINARY | _O_SEQUENTIAL)

// flags when opening existing file for read/write only with mainly random access
#define O_READORWRITERAND  (_O_RDWR | _O_BINARY)

// flags when creating new file, or truncating existing file, with mainly sequential writes
#define O_CREATETRUNC ( _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE)

#else
// flags when opening existing file for read/write with mainly sequential access
#define O_READORWRITESEQ  (O_RDWR)

// flags when opening existing file for read with mainly sequential access
#define O_READSEQ  (O_RDONLY)

// flags when opening existing file for read/write only with mainly random access
#define O_READORWRITERAND  (O_RDWR)

// flags when creating new file, or truncating existing file, with mainly sequential writes
#define O_CREATETRUNC ( O_RDWR | O_CREAT),(S_IREAD | S_IWRITE)

// string functions
#define stricmp strcasecmp
#define strnicmp strncasecmp

// file functions
#define _lseeki64 lseek64

#define _telli64(fd) lseek64(fd,0,SEEK_CUR)

#define	_commit(fd) fsync(fd);

// some name buffer sizes
#ifndef _MAX_FNAME
#define _MAX_FNAME 256
#endif

#ifndef _MAX_PATH
#define _MAX_PATH 260			// was previously 256 but to maintain compatibility with windows changed to 260
#endif

#endif


const unsigned int cBSFTypeAny = 0;				// any type of biosequence file is accepted
const unsigned int cBSFTypeSeq = 1;				// biosequence file contains sequences (equiv to multifasta)
const unsigned int cBSFTypeMultiAlign = 2;		// biosequence file contains multialignments (equiv to UCSC .mfa)
const unsigned int cBSFTypeSfx = 3;				// biosequence file contains suffix arrays
const unsigned int cBSFTypeRlst = 4;			// biosequence file contains intermediate results
const unsigned int cBSFTypeHash = 5;			// biosequence file contains hashed sequences
const unsigned int cBSFTypeDataPts = 6;			// biosequence file contains visualisation data points
const unsigned int cBSFTypeFeat = 7;			// biosequence file contains feature locii
const unsigned int cBSFTypeGOAssoc = 8;			// biosequence file contains GO associations
const unsigned int cBSFTypeGOTerms = 9;			// biosequence file contains GO terms
const unsigned int cBSFTypeFMIndexes = 10;		// biosequence file contains compressed suffix FMIndexes
const unsigned int cBSFTypeError = 128;			// unknown biosequence file type or file header not read

const unsigned int cMBSFFileDescrLen = 1024;	// max length of file description field in file header including '\0'
const unsigned int cMBSFShortFileDescrLen = 64; // max length of short file description field in file header including '\0'

#pragma pack(1)

// nucleotide base - fits in a 4bit nibble
typedef enum eSeqBase {
	eBaseA=0,			// A 
	eBaseC,				// C
	eBaseG,				// G
	eBaseT,				// T or U
	eBaseN,				// N
	eBaseUndef,			// base is undefined - no alignment
	eBaseInDel,			// '-' reserved for use with alignments 
	eBaseEOS,			// end of sequence marker
	eBaseEOG = 0x0f  	// end of genome marker
	} teSeqBases;
typedef uint8_t etSeqBase; // type used to hold teSeqBases
// following is combined with each eSeqBase
const uint8_t cRptMskFlg = 0x08;			// used to flag base as part of a repeat masked sequence
const uint8_t cMarkMskFlg = 0x10;			// used to flag that this base has been marked resulting from internal processing
const uint8_t cMarkConChngFlg = 0x20;			// used to flag that this base has been marked as a changed consensus base 

const int cMaxAutoSeq2Ascii = 0x7fff;			// maximum number of seq2ascii autotranslated


typedef enum eFuncRegion {
	eFRIntergenic = 0,			// 0: intergenic region
	eFRUpstream,				// 1: 5' upstream regulatory region
	eFR5UTR,					// 2: 5'UTR
	eFRCDS,						// 3: CDS
	eFRIntronic,				// 4: Intronic
	eFR3UTR,					// 5: 3'UTR
	eFRDnstream,				// 6: 3' downstream regulatory region
	eFRAuto						// 7: any functional region
} teFuncRegion;


// ontologies classes
typedef enum eOntologies {
	eONTUndefined = 0,	// undefined - error	
	eONTCellular,		// Cellular component
	eONTBiological,		// Biological process
	eONTMolecular,		// Molecular function
	eONTPlantAnatomical,	// plant anatomy
	eONTPlantDev		// plant structure developmental growth
	} etOntologies;

#pragma pack()

// typedefs to make it easy to change common identifier types
typedef int32_t tGeneID;			// gene identifier
typedef int32_t tChromID;			// chromosome identifier
typedef int32_t tAlignBlockID;		// alignment block identifier
typedef int32_t tSpeciesID;			// species identifier

#ifdef _SHORTMAXREADLEN_
const int cMaxReadLen = 2048;				// can process reads (if more than this length then must be sequences!) of up to this length
#else
const int cMaxReadLen = (0x30000);			// reads (in future expected for PacBio fastq and Illumina moleculo fastq format) can be up to this length so need to be able to handle them
#endif

// releases previous to 2.75.x had lower max gene/dataset/species/chrom name limits than later releases
const int cMaxGeneNameLenV3 = 52;				// max gene or element name length (including  trailing '\0')
const int cMaxDatasetSpeciesChromV3 = 36;		// max length of any dataset, species or chromosome name (including  trailing '\0')
const int cMaxDescrIDLenV3 = 50;				// descriptor needs to uniquely identify reads within 1st cMaxDescrIDLen chars
												// descriptors will be truncated at the 1st whitespace char

// 2.75.0 and later name length limits have been increased
const int cMaxGeneNameLen = 81;					// max gene or element name length (including  trailing '\0')
const int cMaxDatasetSpeciesChrom = 81;			// max length of any dataset, species or chromosome name (including  trailing '\0')
const int cMaxDescrIDLen = 80;					// NGS read descriptors needs to uniquely identify reads within 1st cMaxDescrIDLen chars
												// descriptors will be truncated at the 1st whitespace char

extern class CDiagnostics gDiagnostics;			// all diagnostics should be generated through this class
extern char gszProcName[_MAX_FNAME];			// process name using this library

#pragma pack(1)
typedef struct TAG_sBSFRdsHdr {
	uint8_t Magic[4];			 		// magic chars to identify this file as a biosequence reads file
	int64_t RdsOfs;							// where concatenated reads start on disk
	int64_t TotReadsLen;						// total memory space required to hold descriptors + packed read/quality sequences
	int32_t Version;							// file structure version
	int32_t NumRds;							// number of reads in this file after filtering
	int8_t PMode;								// processing mode used when rds file created
	int8_t QMode;								// original fastq quality scoring method
	int8_t Trim5;								// trimmed this many bases off leading sequence 5' end
	int8_t Trim3;								// trimmed this many bases off trailing sequence 3' end
	int8_t NumFiles;							// number of source files in FileNames[]
	int8_t FileNames[8196];					// source file names concatenated with '\0' separators - if too many to fit then overflow names simply sloughed
	int32_t OrigNumReads;						// number of reads originally processed before any deduping or other filtering
	uint32_t FlagsK:1;						// processing flag - bit representing '-k' flag
	uint32_t FlagsCS:1;						// processing flag - bit representing that reads were processed from SOLiD (colorspace)	
	uint32_t FlagsPR:1;						// processing flags - bit representing that contained reads are paired				
} tsBSFRdsHdr;



// V5 raw reads structure
typedef struct TAG_sRawReadV5 {
	uint32_t ReadID;				// unique read identifier
	uint32_t PairReadID;			// this read's paired raw read identifier (0 if not a paired read) bit32 reset if 5' fwd read, bit32 set if 3' read of pair
	uint16_t NumReads;			// number of source reads merged into this read
	uint8_t FileID;				// identifies file from which this read was parsed
	uint8_t DescrLen;				// descriptor length - starts at &Read[0]
	uint16_t ReadLen;				// length of following packed read and quality scores (starts at &Read[DescrLen])
	uint8_t  Read[1];				// descriptor followed by packed read + phred score (read in bits 0..2, phred in 3..7)
} tsRawReadV5;

// V6 raw reads structure
// Only difference from V5 is that NumReads range has been expanded to 32bits from 16bits in V5
typedef struct TAG_sRawReadV6 {
	uint32_t ReadID;				// unique read identifier
	uint32_t PairReadID;			// this read's paired raw read identifier (0 if not a paired read) bit32 reset if 5' fwd read, bit32 set if 3' read of pair
	uint32_t NumReads;			// number of source reads merged into this read
	uint8_t FileID;				// identifies file from which this read was parsed
	uint8_t DescrLen;				// descriptor length - starts at &Read[0]
	uint16_t ReadLen;				// length of following packed read and quality scores (starts at &Read[DescrLen])
	uint8_t  Read[1];				// descriptor followed by packed read + phred score (read in bits 0..2, phred in 3..7)
} tsRawReadV6;

const int cMaxWorkerThreads = 128;			// limiting max number of threads to this many

#pragma pack()
