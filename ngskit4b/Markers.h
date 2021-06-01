/*
This toolkit is a source base clone of 'BioKanga' release 4.4.2 (https://github.com/csiro-crop-informatics/biokanga) and contains
significant source code changes enabling new functionality and resulting process parameterisation changes. These changes have resulted in
incompatibility with 'BioKanga'.

Because of the potential for confusion by users unaware of functionality and process parameterisation changes then the modified source base
and resultant compiled executables have been renamed to 'kit4b' - K-mer Informed Toolkit for Bioinformatics.
The renaming will force users of the 'BioKanga' toolkit to examine scripting which is dependent on existing 'BioKanga'
parameterisations so as to make appropriate changes if wishing to utilise 'kit4b' parameterisations and functionality.

'kit4b' is being released under the Opensource Software License Agreement (GPLv3)
'kit4b' is Copyright (c) 2019, 2020
Please contact Dr Stuart Stephen < stuartjs@g3web.com > if you have any questions regarding 'kit4b'.

Original 'BioKanga' copyright notice has been retained and immediately follows this notice..
*/
/*
 * CSIRO Open Source Software License Agreement (GPLv3)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENSE for the complete license information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */

#pragma once
#include "../libkit4b/commdefs.h"

const int cMaxLenName = 100;			// accept species or chrom/contig names of at most this length
const int cMaxMarkerSpecies = 20000;	// allow at most 20000 different species or cultivars
const int cMaxSeqID = 100000000;		// allow at most 10^8 different target species sequences

// following two constants are used when determining species specific SNPs
const int cAltMaxBases = 1;			// must be no more than this number of bases in alternative species
const double cMinBaseThres = 0.50;  // and species base must be at least this proportion of total bases

const double cDfltLSER = 0.02;			// if local sequencing error rate is unknown then use this default
const double cFloorLSER = 0.001;		// use a floor to prevent divide by 0 errors whilst processing

const uint16_t cFlgSNPcnts = 0x01;		// alignment base counts parsed from SNP file
const uint16_t cFlgImputCnts = 0x02;      // alignment base counts imputed from coverage and assumed to be the target reference base
const uint16_t cFlgAlignCnts = 0x04;		// alignment base counts imputed from alignment sequences


#pragma pack(1)
typedef struct TAG_sAlignLoci {
	int64_t AlignID;				// uniquely identifies this alignment instance (1..N)
	uint16_t TargSpeciesID;		// identifies aligned to cultivar or species
	uint32_t TargSeqID;			// identifies aligned to sequence - could be a chrom/contig/transcript
	uint32_t TargLoci;			// loci within SeqID at which SNPs observed
	uint8_t TargRefBase;			// loci is this reference base
	uint16_t ProbeSpeciesID;		// identifies cultivar or species with sequences aligning to TargSpecies
	uint8_t FiltLowTotBases:1;	// 1 if this alignment has fewer TotBases than reporting threshold
	uint16_t NumSpeciesWithCnts;	// number of species at this loci which have TotBases >= reporting threshold
	uint32_t TotBases;			// sum of base counts in ProbeBaseCnts
	uint8_t CultSpecBase;			// cultivar specific base allowing identification of this cultivar
	double LSER;				// Cultivar local sequencing error rate
	uint16_t Flags;				// any loci associated flags
	uint32_t ProbeBaseCnts[5];	// indexed by A,C,G,T,N : number instances probe base aligned to TargRefBase 
} tsAlignLoci;

const int cAllocAlignLoci = 100000000;	// initially allocate to hold this many alignments
const int cReAllocAlignPerc = 120;	    // if required to realloc then proportionally realloc by this percentage

typedef struct TAG_sSNPSSpecies {
	uint16_t SpeciesID;	// uniquely identifies this species instance
	uint8_t szSpecies[cMaxLenName+1];	// species name (allow for trailing '\0' terminator)
	uint8_t IsRefSpecies:1;	// set if this is a reference or target species, reset if a probe species
} tsSNPSSpecies;


typedef struct TAG_sSeqName {
	uint32_t SeqID;			// uniquely identifies this sequence name (1..N)
	uint8_t Len;				// length of this tsSeqName instance 
	uint64_t NxtSeqOfs;		// offset into m_pAllocSeqNames at which next sequence with same name hash starts
	uint8_t szSeqName[1];		// sequence name including terminating '\0'
} tsSeqName;
const size_t cAllocSeqNames = 5000000;	// allocate to incrementally hold this many sequence names
const size_t cAllocMemSeqNames = (sizeof(tsSeqName) + cMaxLenName) * cAllocSeqNames; // allocate in this sized increments memory (m_pAllocSeqNames) for holding sequence names
const size_t cAllocMinDiffSeqNames = (sizeof(tsSeqName) + cMaxLenName) * 100; // reallocate if less than this many bytes remaining 
#pragma pack()

class CMarkers
{
	tsSNPSSpecies *m_pCurSpecies;			// currently processed species
	uint16_t m_NumSpecies;						// current number of species in m_Species (also includes the reference species)
	uint16_t m_RefSpeciesID;						// identifier for species identified as being the reference species
	tsSNPSSpecies m_Species[cMaxMarkerSpecies+1];	// array of currently known species/cultivars - 1 additional to account for reference species (m_Species[0])

	uint32_t m_NumSeqNames;			// currently there are this many sequence names
	uint64_t m_UsedMemSeqNames;		// memory currently used for sequence names
	uint64_t m_AllocMemSeqNames;		// memory allocated for sequence names
	tsSeqName *m_pAllocSeqNames;	// allocated to hold sequence names

	uint64_t  m_AllocMemSeqNameIDsOfs;	// memory allocated for sequence m_pSeqNameIDsOfs
	uint64_t *m_pAllocSeqNameIDsOfs;		// allocated to hold sequence identifiers to sequence offsets in m_pAllocSeqNames

	uint32_t m_UsedNameHashArray;			// currently using this number of entries in the SeqNameHashArray
	uint64_t *m_pSeqNameHashArray;		// allocated to hold offsets into m_pAllocSeqNames for sequence name hashes

	uint8_t m_szCurSeqName[cMaxLenName+1];	// holds last processed sequence name
	uint32_t m_CurSeqNameID;				// and it's corresponding identifer

	int64_t m_NumSSNPLoci;				// number of loci hold SNPs called in SNP files
	int64_t m_UsedAlignLoci;				// currently using this many alignment loci
	int64_t m_AllocAlignLoci;				// allocated to hold this many alignment loci
	size_t m_AllocMemAlignLoci;			// allocation memory size
	tsAlignLoci *m_pAllocAlignLoci;		// allocated to hold alignment loci 

	int m_hOutFile;						// file handle for outputting snpmarkers 
	char *m_pszBuff;					// snpmarkers buffer
	int m_BuffIdx;						// currently this many chars buffered

	int m_hBEDOutFile;					// file handle for outputting VCF 
	char *m_pszBEDBuff;					// VCF buffer
	int m_BEDBuffIdx;					// currently this many chars buffered

	double m_LSER;						// if local sequencing error rate unknown then use this default

	int m_NumThreads;					// max number of threads

	int64_t AddLoci(uint16_t TargSpeciesID,		// reads were aligned to this cultivar or species
				uint32_t TargSeqID,		// alignments to this sequence - could be a chrom/contig/transcript - from pszSpecies
				uint32_t TargLoci,			// loci within target sequence at which SNPs observed
				etSeqBase TargRefBase,	// loci has this reference base
				uint16_t ProbeSpeciesID,	// reads were aligned from this cultivar or species
				uint32_t ProbeCntA,		// number instances probe base A aligned to TargRefBase 
				uint32_t ProbeCntC,		// number instances probe base C aligned to TargRefBase
				uint32_t ProbeCntG,		// number instances probe base G aligned to TargRefBase
				uint32_t ProbeCntT,		// number instances probe base T aligned to TargRefBase
				uint32_t ProbeCntN,		// number instances probe base N aligned to TargRefBase
				double LSER,			// local sequencing error rate
				uint16_t Flags);			// any loci associated flags


	bool m_bMutexesCreated;			// will be set true if synchronisation mutexes have been created
	int CreateMutexes(void);
	void DeleteMutexes(void);
	void AcquireSerialise(void);
	void ReleaseSerialise(void);
	void AcquireSerialiseMH(void);
	void ReleaseSerialiseMH(void);
	void AcquireLock(bool bExclusive);
	void ReleaseLock(bool bExclusive);

#ifdef _WIN32
	HANDLE m_hMtxIterReads;
	HANDLE m_hMtxMHReads;
	SRWLOCK m_hRwLock;
	HANDLE m_hThreadLoadQuerySeqs;
#else
	pthread_mutex_t m_hMtxIterReads;
	pthread_mutex_t m_hMtxMHReads;
	pthread_rwlock_t m_hRwLock;
#endif

	bool m_bSorted;									// set true if alignments sorted
	CMTqsort m_MTqsort;								// using multithreaded sorting
	static int QSortAlignSeqLociSpecies(const void *arg1, const void *arg2); // qsorts alignment loci by TargSeqID,TargLoci,ProbeSpeciesID ascending

public:
	CMarkers(void);
	~CMarkers(void);
	void Reset(void);	// clears all allocated resources

	int Init(int NumThreads = 1); //Initialise resources for specified number of threads; currently not actually multithreaded but may be in future

	int64_t		// qsorts alignment loci by TargSeqID,TargLoci,ProbeSpeciesID ascending
		SortTargSeqLociSpecies(void);

	int 
		LoadSNPFile(int MinBases,			// accept SNPs with at least this number covering bases
					  double MaxPValue,		// accept SNPs with at most this P-value
					  double SNPMmajorPC,	// only accept SNP for processing if major allele >= this proportion of total allele counts
					  char *pszRefSpecies,	// this is the reference species 
					  char *pszProbeSpecies, // this species reads were aligned to the reference species from which SNPs were called 
					  char *pszSNPFile);	// SNP file to parse and load

	int64_t 
		AddImputedAlignments(int MinBases,		// must be at least this number of reads covering the SNP loci
					  char *pszRefSpecies,		// this is the reference species 
					char *pszProbeSpecies,		// this species reads were aligned to the reference species from which SNPs were called 
					char *pszAlignFile,			// file containing alignments (CSV,BED, or SAM)
					int FType = 0,				// input alignment file format: 0 - auto, 1 - CSV, 2 - BED, 3 - SAM)
					bool bSeqs = false,			// if alignment file contains the read sequence then impute bases from the actual sequences	
					int EstNumSeqs = 0,			// estimated number of sequences (0 if no estimate)
					int EstSeqLen = 0);			// estimated mean sequence length (0 if no estimate)		

	int64_t	// Add simulated alignments for when no actual alignments were available, just SNP calls
		AddSimulatedAlignments(int MinBases,			// using this as the simulated number of reads covering the SNP loci
					 char *pszRefSpecies,				// this is the reference species 
					char *pszProbeSpecies);				// this species reads were aligned to the reference species from which SNPs were called

	uint16_t									// returned species identifer (1..cMaxSpecies)
		AddSpecies(char *pszSpecies,bool IsRefSpecies = false);	// cultivar or species

	char *								// returned species name corresponding to SpeciesID
		SpeciesIDtoName(uint16_t SpeciesID); // species identifier

	uint16_t NameToSpeciesID(char *pszSpecies); // returned species identifier for specified name, returns 0 if name not previously added with AddSpecies)

	bool								// true if a reference or target species
		IsRefSpecies(uint16_t SpeciesID);	// species to test

	uint32_t								// returned sequence identifier (1..cMaxSeqID)
		AddTargSeq(char *pszSeqName);	// sequence name - could be chrom, contig, transcript name

	uint32_t NameToSeqID(char *pszSeqName); // returned sequence identifier for specified name, returns 0 if name not previously added with AddTargSeq)

	char *								// returned sequence name
		SeqIDtoName(uint32_t SeqID);		// sequence identifier for which name is to be returned

	int PreAllocEstSNPs(int64_t EstNumSNPS);	// preallocate memory for this many estimated SNP loci
	int PreAllocImputedSNPs(int NumbIsolates);	// when all known SNPs have been loaded then can allocate for additional imputed SNPs using number of isolates

	int64_t AddLoci(char *pszTargSpecies,	// reads were aligned to this cultivar or species
				char *pszTargSeq,		// alignments to this sequence - could be a chrom/contig/transcript - from pszSpecies
				uint32_t TargLoci,		// loci within target sequence at which SNPs observed
				etSeqBase TargRefBase,	// loci has this reference base
				char *pszProbeSpecies,	// reads were aligned from this cultivar or species
				uint32_t ProbeCntA,		// number instances probe base A aligned to TargRefBase 
				uint32_t ProbeCntC,		// number instances probe base C aligned to TargRefBase
				uint32_t ProbeCntG,		// number instances probe base G aligned to TargRefBase
				uint32_t ProbeCntT,		// number instances probe base T aligned to TargRefBase
				uint32_t ProbeCntN,		// number instances probe base N aligned to TargRefBase
				double LSER);			// local sequencing error rate

	int IdentSpeciesSpec(int AltMaxCnt,	// max count allowed for base being processed in any other species, 0 if no limit
						int MinCnt,		// min count required for base being processed in species
						double SNPMmajorPC,		// to be processed major putative SNP base must be at least this percentage of total
						int MinSpeciesWithCnts = 0,			// must be at least this number of species with base counts more than MinSpeciesTotCntThres - 0 if no limit 
						int MinSpeciesTotCntThres = 0);		// individual species must have at least this number of total bases at SNP loci - 0 if no threshold

	int64_t NumAlignLoci(void);					// returns current number of alignment/SNP loci

	int64_t											// number of markers reported
		Report(char *pszRefGenome,			    // reference genome assembly against which other species were aligned
			int NumRelGenomes,					// number of relative genome names
			char *pszRelGenomes[],				// relative genome names
			char *pszReportFile,				// report to this file
			int MinSpeciesWithCnts = 0,			// must be at least this number of species with base counts more than MinSpeciesTotCntThres - 0 if no limit 
			int MinSpeciesTotCntThres = 0,  	// individual species must have at least this number of total bases at SNP loci - 0 if no threshold
			bool bSloughRefOnly = false,		// do not report if no inter-cultivar SNP marker, i.e if cultivars all same with the polymorthic site relative to reference only
			bool bSloughNonHetero = false);		// do not report unless all cultivars are relative heterozygotic - no two cultivars have same base 
};


