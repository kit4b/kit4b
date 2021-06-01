#pragma once

const int cMaxOjCompLen = cMaxGeneNameLen;
const int cMaxLinkageEvidence = 101;			// allowing upto this length of linkage evidence
const int cAllocAGPentries = 10000;				// alloc/realloc for AGP entries in this many increments

// The sequencing status of the component. These typically correspond to keywords in the International Sequence Database (GenBank/EMBL/DDBJ) submission
typedef enum TAG_eAGPseqStatus {
	eAGPSSA=0,			// Active Finishing
	eAGPSSD,		// Draft HTG (often phase1 and phase2 are called Draft, whether or not they have the draft keyword).
	eAGPSSF,		// Finished HTG (phase 3)
	eAGPSSG,		// Whole Genome Finishing
	eAGPSSN,		// gap with specified size
	eAGPSSO,		// Other sequence (typically means no HTG keyword)
	eAGPSSP,		// Pre Draft
	eAGPSSU,		// gap of unknown size, typically defaulting to predefined values.
	eAGPSSW,		// WGS contig
} teAGPseqStatus;

// The combination of gap type and linkage (column 8b) indicates whether the gap is captured or uncaptured. In some cases, the gap types are assigned a biological value (e.g. centromere).
typedef enum TAG_eAGPgapType {
	AGPfragment = 0,      // gap between two sequence contigs (also called a ‘sequence gap’).
	AGPclone,			  // a gap between two clones that do not overlap.
	AGPcontig,			  // a gap between clone contigs (also called a "layout gap").
	AGPcentromere,		  // a gap inserted for the centromere.
	AGPshort_arm,		  // a gap inserted at the start of an acrocentric chromosome.
	AGPheterochromatin,	  // a gap inserted for an especially large region of heterochromatic sequence (may also include the centromere).
	AGPtelomere,		  // a gap inserted for the telomere.
	AGPrepeat,			  // an unresolvable repeat.
	AGPscaffold,		  // a gap between two sequence contigs in a scaffold (superscaffold or ultra-scaffold)
	AGPGapUndefined,	  // gap type is unknown - a catchall for types not defined in the AGP 1.1 and 2.0 specification
	} teAGPgapType;

typedef enum TAG_eAGPcompOrient {
	eAGPOplus = 0,	// plus
	eAGPOminus,		// minus
	eAGPOunknown,	// unknown
	eAGPna			// na
} teAGPcompOrient;

typedef struct TAG_sAGPgap {
		uint32_t GapLen;		// gap length 
		teAGPgapType GapType;		// gap type
		bool bLinkage;		// true if there is evidence of linkage between the adjacent lines
		char szEvidence[cMaxLinkageEvidence];	// evidence of linkage
} tsAGPgap;

typedef struct TAG_sAGPcomp {
	uint32_t Start;		// start offset in component
	uint32_t End;			// end offset in component
	teAGPcompOrient Orientation;	// component orientation relative to object
	char szCompID[cMaxOjCompLen];			// component identifier
} tsAGPcomp;


typedef struct TAG_sAGPentry {
	uint32_t EntryID;			// uniquely identifies this entry (1..m_NumAGPentries)
	char szObjIdent[cMaxOjCompLen];	// identifer of object being assembled
	uint32_t Start;			// object start loci (1 based)
	uint32_t End;				// object end loci
	uint32_t PartNum;			// identifies this entry line for object (1..n)
	teAGPseqStatus Type;		// component type - A,D,F,G,N,O,P,U,W
	union {
		tsAGPgap Gap;		// if Type is N or U
		tsAGPcomp Comp;		// if Type is not N or U
		};	
} tsAGPentry;


class CAGPs
{
	uint32_t m_NumAGPentries;		// current number of AGP entries ptd at by m_pAGPentries
	uint32_t m_AllocAGPentries;		// number of AGP entries allocated
	size_t m_AllocdAGPMem;		// size of memory allocated to hold AGP entries
	tsAGPentry *m_pAGPentries;  // used to hold parsed AGP entries

	char *TrimWhitespace(char *pTxt); // trim leading/trailing whitespace

	static int SortAGPentries(const void *pEl1,const void *pEl2);

public:
	CAGPs(void);
	~CAGPs(void);
	void Reset(void);
	int LoadAGPs(char *pszAGPFile);		 // load AGPs from specified file
	tsAGPentry *Entry(uint32_t EntryID);	 // returns specified entry
	tsAGPentry *LocateCompID(char *pszCompID); // returns entry with specified component identifier
	tsAGPentry *Next(tsAGPentry *pPrev); // returns next entry (NULL if none) after pPrev. If pPrev == NULL then returns first

};
