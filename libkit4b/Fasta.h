#pragma once
#include "./commdefs.h"

/*
Fastq scoring schema (from http://en.wikipedia.org/wiki/FASTQ_format )

SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS.....................................................
..........................XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX......................
...............................IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII......................
.................................JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ......................
..LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL....................................................
!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
|                         |    |        |                              |                     |
33                        59   64       73                            104                   126
0........................26...31.......40
-5....0........9.............................40
0........9.............................40
3.....9.............................40
0.2......................26...31........41

S - Sanger        Phred+33,  raw reads typically (0, 40)
X - Solexa        Solexa+64, raw reads typically (-5, 40)
I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)
J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40)
with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator (bold)
(Note: See discussion above).
L - Illumina 1.8+ Phred+33,  raw reads typically (0, 41)
*/

const uint32_t cDfltStageBuffSize = 0x03ffffff;  // 64M buffer as default
const uint32_t cMaxStageBuffSize = 0x03fffffff;	 // 1GB buffer as maximum
const uint32_t cMinStageBuffSize = 0x0fffff;	 // 1M buffer as minimum

const int32_t cNumFastaBlocks = 1;						 // when readahead disk buffering implemented then set to be 2 or more

const uint32_t cMaxGenFastaLineLen = 79;		// limit generated Fasta lines to this length
const uint32_t cMaxFastaDescrLen   = 8192;	    // Fasta descriptor lines can be concatenated..
const uint32_t cMaxFastQSeqLen = cMaxReadLen;	// Allow FastQ sequences to be of this max length

const int32_t cgzAllocInBuffer = 0x1ffffff;				// gz processing input buffer size
const int32_t cgzAllocOutBuffer = 0x1ffffff;			// gz processing output buffer size

#pragma pack(1)
typedef struct TAG_sFastaBlock
	{
	int64_t  FileOfs;				// file offset from where fasta block (m_pCurBuffer) last read from file
	int32_t BuffIdx;				// process next char from pBlock[BuffIdx]
	int32_t BuffCnt;				// block currently loaded with this many chars
	int32_t AllocSize;			// block was allocated to buffer at most this many chars in Fasta[]
	uint8_t *pBlock;				// allocated to hold a block of fasta file content
	} tsFastaBlock;
#pragma pack()

class CFasta : public CErrorCodes
{
	int32_t m_hFile;				// opened for write fasta
	gzFile m_gzFile;			// opened for read (could be compressed) fasta or fastq
	char m_szFile[_MAX_PATH];	// to hold fasta file path+name
	uint64_t m_StatFileSize;		// file size as returned by stat() when file initially opened
	bool m_bIsGZ;				// true if processing a gz compressed file
	bool m_bIsFastQ;			// true if processing a fastq file
	bool m_bIscsfasta;			// sequences are in SOLiD csfasta format
	static uint8_t m_SOLiDmap[5][5]; // used for transforming from SOLiD colorspace into basespace
	bool m_bRead;				// TRUE if reading fasta file, FALSE if write to fasta file

	tsFastaBlock *m_pCurFastaBlock;    // buffered fasta block currently being processed
	tsFastaBlock m_FastaBlocks[cNumFastaBlocks];    // allow for at most cNumFastaBlocks buffered fasta file blocks; currently not implemented but in future will allow for readahead of blocks

	bool m_bForceRetSeq;        // true if 1bp dummy sequence is to be returned even if no sequence is actually available
	bool m_bDescrAvail;			// true if NEW descriptor available, reset by ReadDescriptor()
	char m_szDescriptor[cMaxFastaDescrLen+1];	// to hold last descriptor parsed
	char *m_pszFastqSeq;	    // to hold last FastQ sequence line parsed
	char *m_pszFastqSeqQ;	    // to hold last FastQ quality line parsed
	int32_t m_FastqSeqLen;
	int32_t m_FastqSeqIdx;
	int32_t	m_FastqSeqQLen;
	int32_t m_CurFastQParseLine;	// current FastQ line being processed
	int32_t m_NumMissingSequences;	// total number of fastq records which were empty with no sequence present

	int64_t m_FileDescrOfs;		// file offset from where last descriptor was parsed
	int64_t m_FileReadDescrOfs;   // offset of last descriptor returned by ReadDescriptor()
	uint32_t m_DescriptorLen;
	uint32_t m_CurLineLen;

	int32_t CheckIsFasta(void);		// checks if file contents are likely to be fasta or fastq format
	int32_t	ParseFastQblockQ(void); // Parses a fastq block (seq identifier + sequence + quality scores)

public:
	CFasta(void);
	~CFasta(void);
	void Cleanup(void);
	int32_t Reset(int64_t FileOfs = 0l);				// reset context to that following an Open() with option to start processing at FileOfs
												// NOTE: gzip library can't handle file offsets which are greater than 2^31 - 1

	int32_t Open(char *pszFile,bool Read = true,uint32_t BufferSize = cDfltStageBuffSize);
	bool IsFastq(void);							// true if opened file is in fastq format
	bool IsSOLiD(void);							// true if opened file is in SOLiD or colorspace format
	uint64_t InitialFileSize(void);				// file size when initially opened for reading

	uint32_t										// returns estimated number of sequences in fasta/fastq file
		FastaEstSizes(char *pszFile,			// fasta or fastq file path+name to estimate sizes
			  int64_t *pFileSize = NULL,			// file is this size on disk
			  int32_t *pEstMaxDescrLen = NULL,	// with estimated maximum descriptor length
			  int32_t *pEstMeanDescrLen = NULL,	// estimated mean descriptor length
			  int32_t *pEstMaxSeqLen = NULL,		// and estimated maximum sequence length
			  int32_t *pEstMeanSeqLen = NULL,		// estimated mean sequence length
			  int32_t *pEstScoreSchema = NULL);	// guestimated scoring schema - 0: no scoring, 1: Solexa, 2: Illumina 1.3+, 3: Illumina 1.5+, 4: Illumina 1.8+ or could be Sanger 
	
	int32_t Close(void);							 // close opened fasta file
	int32_t	LocateDescriptor(char *pszPrefixToMatch=NULL,// prefix to match or NULL if match any
						 bool bFromStart=false);	 // true: from start, false: from next descriptor

	int32_t											// number actually read
		ReadSubsequence(etSeqBase *pSeq,		// where to return subsequence (NO TERMINATING eBaseEOS)
					 uint32_t MaxLen,		// reads upto MaxLen sequence bases from file
					 uint32_t SeqOfs=0,		// relative offset in sequence at which to read
					 bool bFromStart=false);	// true: from fasta file start, false: from next sequence

	int32_t ReadSequence(void *pRetSeq = NULL,		// where to return sequence, can be NULL if only interested in the sequence length
					 int32_t Max2Ret = 0,			// max to return, ignored if pRetSeq is NULL
					 bool bSeqBase = true,		// if false then return ascii, if true then return as etSeqBase
					 bool RptMskUpperCase = false);	// default is false, UCSC softmasked use lowercase when repeat masked


	int32_t ReadQValues(char *pRetQValues,	// where to return quality values
					 int32_t Max2Ret);			// max to return

	int32_t ReadDescriptor(char *pszDescriptor,int32_t MaxLen); // copies last descriptor processed into pszDescriptor and returns copied length
	int64_t GetDescrFileOfs(void);				// returns file offset at which descriptor returned by ReadDescriptor() was parsed from

	static int32_t Ascii2Sense(char *pAscii,		  // expected to be '\0' terminated, or SeqLen long
					int32_t MaxSeqLen,				  // maximal sized sequence that pSeq can hold
					etSeqBase *pSeq,			  // where to return sequence
					bool RptMskUpperCase=false);  // true if bases are softmasked as uppercase instead of default lowercase

	static char Base2Chr(etSeqBase Base);		 // returns ascii representation of Base
	static etSeqBase Chr2Base(char Base);			 // returns etSeqBase representation of ascii base

	int32_t Write(char Symbol);							// writes sequence char to fasta file
	int32_t Write(char *pSymbols,uint32_t Cnt);		// writes Cnt sequence chars to fasta file
	int32_t WriteDescriptor(char *pszDescriptor);		// terminates any current sequence and writes fasta file descriptor starting new sequence
};

