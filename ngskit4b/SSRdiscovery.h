/*
This toolkit is a source base clone of 'BioKanga' release 4.4.2 (https://github.com/csiro-crop-informatics/biokanga) and contains
significant source code changes enabling new functionality and resulting process parameterisation changes. These changes have resulted in
incompatibilty with 'BioKanga'.

Because of the potentential for confusion by users unaware of functionality and process parameterisation changes then the modified source base
and resultant compiled executables have been renamed to 'kit4b' - K-mer Informed Toolkit for Bioinformatics.
The renaming will force users of the 'BioKanga' toolkit to examine scripting which is dependent on existing 'BioKanga'
parameterisations so as to make appropriate changes if wishing to utilise 'kit4b' parameterisations and functionality.

'kit4b' is being released under the Opensource Software License Agreement (GPLv3)
'kit4b' is Copyright (c) 2019
Please contact Dr Stuart Stephen < stuartjs@g3web.com > if you have any questions regarding 'kit4b'.
 */

#pragma once

const size_t cMaxAllocBuffChunk = 0x0ffffff;		// allocate for input sequences buffering in these sized chunks
const int cMaxAllocRptSSRs = 0x07fffff;				// SSR reporting buffer size

const int cMinRepElLen = 1;				// minimum element K-mer length
const int cDfltMinRepElLen = 2;			// default minimum element K-mer length
const int cDfltMaxRepElLen = 5;			// default maximum element k-mer length
const int cMaxRepElLen = 20;			// maximum element k-mer length

const int cMinTandemRpts = 2;			// min number of tandem repeats
const int cDfltMinTandemRpts = 5;		// default min tandem repeats
const int cDfltMaxTandemRpts = 10;		// defualt max tandem repeats
const int cMaxTandemRpts = 50;			// max number of tandem repeats

const int cMinSSRFlankLen = 25;			// minimum user specified SSR flanking lengths
const int cDfltSSRFlankLen = 100;		// default SSR flanking lengths
const int cMaxSSRFlankLen = 500;		// maximum user specified SSR flanking lengths

// reporting of SSR loci is in one of the following formats
typedef enum eRptSSRsFromat {
	eRFCsv = 0,			// report as CSV
	eRFBed,				// BED file format
	eRFSam				// SAM file format
	} teRptSSRsFromat;

#pragma pack(1)
typedef struct TAG_sKMerDist {
	UINT32 Cnt;							// number of occurances of SSR with this repeating K-mer element
	UINT32 TandemRpts[1];				// number of tandem repeats (extended to user specified max repeats)
	} tsKMerDist;
#pragma pack()

class CSSRDiscovery
{
	CStopWatch m_CurTime;			// used for progress messaging
	teRptSSRsFromat m_RptSSRsFormat;	// format in which to report SSRs
	int m_hOutFile;				// output processing results to this file
	int m_hOutKMerFreqFile;		// output repeating element KMer freqs to this file 

	int m_IdxRptSSRs;			// current index into m_pszRptSSRsBuff for buffered SSRs reporting
	char *m_pszRptSSRsBuff;		// allocated to buffer SSRs reporting

	size_t m_SeqBuffLen;		// number of bases currently buffered in m_pSeqBuff
	size_t m_AllocdSeqBuffMem;  // size of memory currently allocated to m_pSeqBuff
	UINT8 *m_pSeqBuff;			// buffers sequences as read from file

	int m_MinRepElLen;			// identify repeating elements of this minimum length
	int m_MaxRepElLen;			// ranging upto this maximum length
	int m_MinTandemRpts;		// minimum number of tandem repeats
	int m_MaxTandemRpts;		// maximum number of repeats
	int m_SSRFlankLen;			// SSR flanking lengths

	UINT32 m_TotNumAcceptedSSRs;	// total number of SSRs accepted
	UINT32 m_TotNumExcessiveTandemSSRs;    // total number of putative SSRs not accepted because m_MaxTandemRpts

	UINT32 m_TotNumAcceptedKmerSSRs[cMaxRepElLen+1];

	int m_KMerFreqLen;				// if > 0 then counting this length KMer
	int m_NumKMers;					// number of KMers
	int m_SizeOfKMerDist;			// size of each tsKMerDist allowing for max number of tandem repeats requested by user
	size_t m_AllocdKMerFreqMem;		// size of memory currently allocated to m_pKMerFreq
	tsKMerDist *m_pKMerDist;		// allocated to hold KMer sequence instance counts; only used if single KMer element length being processed

	int ProcessBioseqFile(char *pszFile);	// load and process a bioseq file for SSRs
	int ProcessFastaFile(char *pszFile);	// load and process a multifasta file for SSRs

	int Report(int RepElLen,			// identified SSR contains repeat elements of this length
			   int NumTandemEls,		// each repeat element is repeated this many times
			   INT64 SSRStartOfs,		// repeat element starts at this offset within the targeted sequence
				char *pszDescr,			// descriptor for the targeted sequence
				char *pszInFile,		// sequence parsed from this file
			   INT64 TargSeqLen,		// targeted sequence is this length
			   etSeqBase *pTargSeq);	// targeted sequence within which the SSR has been located

	int
		ReportCSV(int RepElLen,			// identified SSR contains repeat elements of this length
			int NumTandemEls,			// each repeat element is repeated this many times
			INT64 SSRStartOfs,			// repeat element starts at this offset within the targeted sequence
			char *pszDescr,				// descriptor for the targeted sequence
			char *pszInFile,			// sequence parsed from this file
			INT64 TargSeqLen,			// targeted sequence is this length
			etSeqBase *pTargSeq);		// targeted sequence within which the SSR has been located

	int
		ReportBED(int RepElLen,			// identified SSR contains repeat elements of this length
			int NumTandemEls,			// each repeat element is repeated this many times
			INT64 SSRStartOfs,			// repeat element starts at this offset within the targeted sequence
			char *pszDescr,				// descriptor for the targeted sequence
			char *pszInFile,			// sequence parsed from this file
			INT64 TargSeqLen,			// targeted sequence is this length
			etSeqBase *pTargSeq);		// targeted sequence within which the SSR has been located

	int
		ReportSAM(int RepElLen,			// identified SSR contains repeat elements of this length
			int NumTandemEls,			// each repeat element is repeated this many times
			INT64 SSRStartOfs,			// repeat element starts at this offset within the targeted sequence
			char *pszDescr,				// descriptor for the targeted sequence
			char *pszInFile,			// sequence parsed from this file
			INT64 TargSeqLen,			// targeted sequence is this length
			etSeqBase *pTargSeq);		// targeted sequence within which the SSR has been located

	int	ReportKMers(char *pszKMerFreqFile);	// report SSR repeating element K-mer frequencies to this file

	int	ReportProgress(bool bForce = false);	// let user know that there is processing activity, normally progress reportde evry 60 sec unless bForce set true

	int
		IdentifySSRs(
			 char *pszDescr,			// descriptor for the targeted sequence
			 char *pszInFile,			// sequence parsed from this file
			 INT64 TargSeqLen,			// sequence length of targeted sequence within which to search for SSRs
			 etSeqBase *pTargSeq);		// identify SSRs in this targeted sequence

	etSeqBase *AllocSeqBuff(size_t SeqLen);	// allocate for at least this sequence length

	int CntKMer(int KMerLen,			// count this KMer
				int Rpts,				// tandem repeat counts
				etSeqBase *pSeq);		// sequence 

public:
	CSSRDiscovery(void);
	~CSSRDiscovery(void);

	void Init(void);						// initialisation
	void Reset(void);					// resets state back to that imediately following initialisation
	int
		Process(int PMode,			// procesisng mode - currently unused..
			teRptSSRsFromat RptSSRsFormat,	// report SSRs in this file format
			int MinRepElLen,	// identify repeating elements of this minimum length
			int MaxRepElLen,		// ranging upto this maximum length
			int MinTandemRpts,		// minimum number of tandem repeats
			int MaxTandemRpts,		// maximum number of repeats
			int SSRFlankLen,		// SSR flanking sequences lengths to report
			int NumInFileSpecs,		// number of input, could be wildcarded, file specs
			char *pszInFile[],		// files to be processed
			char *pszKMerFreqFile,	// optional, output element KMer freq to this file
			char *pszOutFile);		// SSRs to this file

};



