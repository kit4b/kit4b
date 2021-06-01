#pragma once
#include "./commdefs.h"

const uint32_t cBSFVersion = 10;					// increment each time the file strucure for bioseq or suffix is changed
const uint32_t cBSFVersionBack = 10;				// minimum version that is this version is still compatible with

const uint32_t cBSFMaxDirEntries = 20000000;	// entry directory can contain upto this many entries
const uint32_t cBSFDataBuffSize     = 0x07ffff; // buffer size used when creating entries
const uint32_t cAllocEntriesIncr = 0x0fffff;	// alloc this size when reallocing directory
const uint32_t cBSFSourceSize = 255;			 // max string size for source of each entry 					
const uint32_t cBSFDescriptionSize = 4095;		 // max string size for describing each entry

const uint32_t cMaxFileReadChunk = 2000000000;	 // read in these chunk sized blocks to overcome 2G read limit

typedef int32_t tBSFEntryID;							// entry identifiers

#pragma pack(1)

typedef enum eDataType {
	eUndefDataType = 0,
	eSeqBaseType,				// data is etSeqBase packed to 4bits/base
	eAsciiType,					// data is 7bit ascii
	eBinDataType,				// data is raw 8bit binary
	eDataPts					// data is visualisation data points
} teDataType;
typedef uint8_t etDataType;

typedef struct TAG_sBSFDirEntry {
	int64_t DataPsn;						    // absolute file psn at which biosequence data starts
	uint32_t Size;							// size in bytes of this instance when concatenated
	int32_t EntryID;							// 1..n uniquely identifies this entry instance
	int32_t NameInst;							// name instance (1..n) -- can have multiple entries with duplicate names
	uint32_t DataLen;							// data length when unpacked
	uint16_t Hash;							// hash on szName
	uint8_t DType:4;							// data type - teDataType
	uint8_t MskSense:1;						// 1 if uppercase represents repeat masked, 0 if lowercase represents repeat masked
	char szName[1];							// place holder for entry name + '\0' + optional description + '\0'
} tsBSFDirEntry;
#pragma pack()

#pragma pack(8)
typedef struct TAG_sBSFHeader {
	uint8_t Magic[4];			 		// magic chars to identify this file as a biosequence file
	int64_t FileLen;							// current file length (write psn for newly created entries)
	int64_t SeqOfs;							// where concatenated sequences start on disk
	int64_t SeqSize;							// disk/memory space required to hold concatenated sequences
	int64_t DirEntriesOfs;					// where directory entries start on disk
	int64_t DirEntryIDIdxOfs;					// where index by EntryID into directory entries start on disk
	int64_t DirEntryNameIdxOfs;				// where index by name into directory entries start on disk
	int32_t Type;								// biosequence file type
	int32_t Version;							// file structure version
	int32_t MaxEntries;						// max number of entries supported 
	int32_t NumEntries;						// actual number of directory entries
	int32_t DirEntriesSize;					// actual disk/memory space required to hold directory entries
	char szDatasetName[cMaxDatasetSpeciesChrom]; // dataset name
	char szDescription[cMBSFFileDescrLen];	// describes contents of file
	char szTitle[cMBSFShortFileDescrLen];	// short title by which this file can be distingished from other files in dropdown lists etc
} tsBSFHeader;

#pragma pack()


class CBioSeqFile : protected CEndian, public CConformation
{
	int m_hFile;							   // opened/created file handle
	char m_szFile[_MAX_PATH+1];				   // file name	as opened/created
	bool m_bHdrDirty;						   // TRUE if current header needs to be written to disk 
	bool m_bCreate;								// TRUE if file opened in create mode
	tsBSFHeader m_FileHdr;					   // file header

	int m_AllocDirEntriesSize;					// disk/memory space required to hold concatenated directory entries
	tsBSFDirEntry *m_pDirEntries;				// pts to  array of tsBSFDirectory's
	tsBSFDirEntry **m_ppDirEntryNames;			// sorted (by name) array of ptrs into m_pDirEntries
	tsBSFDirEntry **m_ppDirEntryIDs;			// sorted (by EntryID) array of ptrs into m_pDirEntries

	tsBSFDirEntry *m_pCreatingEntry;		   // entry currently being created, NULL if none
	uint32_t m_DataBuffLen;				   // number of bytes of data currently in m_pDataBuff
	uint32_t m_DataBuffNibs;			   // number of nibbles currently in m_pDataBuff
	uint8_t *m_pDataBuff;				   // data buffer for use whilst entries are created
	uint64_t m_DataWrtPsn;						// file offset at which to next write contents of m_pDataBuff
	teBSFrsltCodes LoadEntries(void);			// load entries directory into memory
	teBSFrsltCodes ReadDisk(int64_t DiskOfs,int Len,void *pTo); // reads disk block into memory
	teBSFrsltCodes Hdr2Disk(void);			// write header to disk
	teBSFrsltCodes Disk2Hdr(char *pszSeqFile,int FileType); // read header from disk

protected:
	void InitHdr(uint32_t Type);			// initialise file header with biosequence file type
	teBSFrsltCodes Flush2Disk(void);					   // flush and commit to disk
	tsBSFDirEntry *LocateEntry(tBSFEntryID EntryID); // locate entry by EntryID
	teBSFrsltCodes SortEntries(void);
	unsigned short GenNameHash(char *pszName);
	static int SortEntryNames(const void *arg1, const void *arg2); // used to sort by entry name->EntryID


public:
	CBioSeqFile(void);
	virtual ~CBioSeqFile(void);
	teBSFrsltCodes Reset(bool bFlush = true);		// reset state back to that immediately following instantiation before file opened/created
	virtual teBSFrsltCodes Open(char *pszSeqFile,	// specifies file to open or create
				uint32_t Type = eSeqBaseType,	// biosequence file type 
				bool bCreate = false);				// create file if it doesn't already exist, truncate if it does

	teBSFrsltCodes SetDescription(char *pszDescription);
	teBSFrsltCodes GetDescription(int MaxLen,char *pszDescription);
	teBSFrsltCodes SetTitle(char *pszTitle);
	teBSFrsltCodes GetTitle(int MaxLen,char *pszTitle);

	int GetType(void);
	teBSFrsltCodes Close(void);
	teBSFrsltCodes SetDatasetName(char *pszDataset);	// sets file dataset name
	char *GetDatasetName(void);						// returns file dataset name
	int LocateEntryIDbyName(char *pszName,			// entry name to locate
							int Ith = 1);			// Ith instance to locate (1..n)
	teBSFrsltCodes Exists(tBSFEntryID EntryID);	// returns eBSFSuccess if specified entry exists
	int NumEntries(void);						// returns number of entries
	tBSFEntryID Next(tBSFEntryID Current);		// returns next entry identifer after current
	tBSFEntryID Prev(tBSFEntryID Current);		// returns previous entry identifer before current

	int GetName(tBSFEntryID EntryID,int MaxLen,char *pszName); // get entry name
	int GetDescription(tBSFEntryID EntryID,int MaxLen,char *pDescr); // get entry description
	teBSFrsltCodes GetNameDescription(tBSFEntryID EntryID,int MaxNameLen,char *pszName,int MaxDescrLen,char *pszDescr);
	
	etDataType GetDataType(tBSFEntryID EntryID);
	uint32_t GetDataLen(tBSFEntryID EntryID);		// get entry data length (returns 0 if no associated data or data unavailable)
	int GetData(tBSFEntryID EntryID,etDataType ReqDataType,uint32_t Ofs,uint8_t *pBuffer,uint32_t MaxLen);

	tBSFEntryID CreateEntry(char *pszName, char *pszDescription,etDataType DataType, bool RptMskUpperCase=false);

	teBSFrsltCodes AddData(uint32_t DataLen, uint8_t *pData);		// add data to the currently created entry
	teBSFrsltCodes SealEntry(void);							// seal currently created entry and commits to file - no more data can be added to that entry						

	int GetShuffledBases(int ChromID,	// reference chromosome identifier
					  int ChromOfs,			// start on reference chromosome (0..n) 
					  int NumPts,			// number of points to return
					  uint8_t *pToRet,	// where to return etSeqBases
					  etSeqBase MissingMarker);	// value to use if missing bases

	void ShuffleDinucs(etSeqBase *pSeq,int SeqLen);
	bool Dump2XML(char *pszXMLfile,  uint32_t MaxDumpSeqLen);

	static int  LocateBase(etSeqBase Probe, int TargLen,etSeqBase *pTarget);
	static int LocateSequence(int ProbeLen, etSeqBase *pProbe,int TargLen,etSeqBase *pTarget);

	static int GetNumMasked(int ProbeLen, etSeqBase *pProbe);

	static int
		PackBases(uint32_t SeqLen,	// unpacked sequence length 
		  etSeqBase *pUnpacked,				// pts to sequence to pack from
		  uint32_t NibbleOfs,			// nibble to start pack into (0..SeqLen-1)
		  uint8_t *pPacked);			// where to pack into
	static int
		UnpackBases(uint32_t SeqLen,	// unpacked sequence length 
		  etSeqBase *pUnpacked,				// pts to sequence to unpack into
		  uint32_t NibbleOfs,			// nibble to start unpack from (0..SeqLen-1)
		  uint8_t *pPacked);				// where to unpack from


	static uint32_t ClassComplexity(etSeqBase *pSeq,
							 uint32_t Len,
							 uint32_t ThresPercent);

	static uint32_t ScoreComplexity(etSeqBase *pSeq,
							 uint32_t Len);

	static uint32_t								// number of aligned bases
		QuickAlignRight(uint32_t AllowedMismatches,	// total allowed mismatches
		   uint32_t MaxMismatchSeqLen,	// max number of mismatches allowed in any run 
		   uint32_t AllowedProbeInDels,	// total allowed indel events on probe
		   uint32_t AllowedTargInDels,	// total allowed indel events on target
		   uint32_t MaxInDelExtension,	// max length of any InDel extension 
		   uint32_t ProbeLen,			// remaining probe length
		   uint8_t *pProbe,
		   uint32_t TargLen,			// remaining target length
		   uint8_t *pTarg);

	static uint32_t						// number of aligned bases
		QuickAlignLeft(uint32_t AllowedMismatches,	// total allowed mismatches
		   uint32_t MaxMismatchSeqLen,	// max number of mismatches allowed in any run 
		   uint32_t AllowedProbeInDels,	// total allowed indel events on probe
		   uint32_t AllowedTargInDels,	// total allowed indel events on target
		   uint32_t MaxInDelExtension,	// max length of any InDel extension 
		   uint32_t ProbeLen,			// remaining probe length
		   uint8_t *pProbe,
		   uint32_t TargLen,			// remaining target length
		   uint8_t *pTarg);

	static char * MakeXMLsafe(char *pszStr);
};

