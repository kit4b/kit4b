#pragma once

const int cMARawBuffSize = 0x1fffffff;
const int cMALineSize = 0x0fffffff;


typedef struct TAG_sProcParams
{
	bool bIsAXT;			// true if source files are assumed to be AXT, otherwise (false) they are MAFs
	bool bSwapRefRel;		// exchange reference and relative sequences
	char* pszRefSpecies;	// reference species name
	char* pszRelSpecies;	// relative species name
	char* pRawBuffer;		// pts to raw char buffer from which MAF/AXT lines are processed
	char* pszLineBuffer;	// pts to line buffer
	char* pszAXTBuffer;		// pts to line buffer for use by AXT processing 
	CMAlignFile* pAlignFile; // bioseq multialignment file
} tsProcParams;

class CGenMAFAlgn
{
	CAlignValidate* m_pAlignValidate;

	bool	IsBase(char Base);
	char*	StripWS(char* pChr);
	
	int	ProcessThisFile(char* pszFile,
			void* pParams);			// will be tsProcParams
	int	ProcessMAFline(int LineLen, tsProcParams* pProcParams);
	bool // Inplace reverse complement pszSeq
		ReverseComplement(int SeqLenInclInDels, char* pszSeq);
	void	RevSeq(int SeqLenInclInDels, char* pszSeq);

	int	ProcessAXTline(int LineLen, tsProcParams* pProcParams);

public:
	CGenMAFAlgn();
	~CGenMAFAlgn();
	void Reset(void);

	int	GenbioDataPointsFile(char* pszMAF, char* pszDataPointsFile, char* pszDescr, char* pszTitle);
	int	// Create a multialignment file from files (either axt or mfa format) in specified source directory
		CreateMAlignFile(char* pszChromLens, char* pszSrcDirPath, char* pszDestAlignFile,
			char* pszDescr, char* pszTitle,
			char* pszRefSpecies, char* pszRelSpecies, bool bIsAXT, bool bSwapRefRel);

};

