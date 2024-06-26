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

#include "stdafx.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#if _WIN32
#include <process.h>
#include "../libkit4b/commhdrs.h"
#else
#include <sys/mman.h>
#include <pthread.h>
#include "../libkit4b/commhdrs.h"
#endif
#include "Markers.h"

CMarkers::CMarkers(void)
{
m_hOutFile = -1; 
m_pszBuff = nullptr;
m_hBEDOutFile = -1;
m_pszBEDBuff = nullptr;
m_pAllocSeqNames = nullptr;
m_pAllocSeqNameIDsOfs = nullptr;
m_pSeqNameHashArray = nullptr;
m_pAllocAlignLoci = nullptr;
m_bMutexesCreated = false;
Reset();
}


CMarkers::~CMarkers(void)
{
Reset();
}


// this process is not currently not actually multithreaded but one day, real soon, it will be multithreaded :-)
int
CMarkers::Init(int NumThreads)	//Initialise resources for specified number of threads
{
int Rslt;
m_NumThreads = NumThreads;
if((Rslt=CreateMutexes())!= eBSFSuccess)
	return(Rslt);
return(eBSFSuccess);
}

// Reset
// Release all allocated resources
void 
CMarkers::Reset(void)	// clears all allocated resources
{
if(m_pAllocSeqNames != nullptr)
	{
#ifdef _WIN32
	free(m_pAllocSeqNames);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pAllocSeqNames != MAP_FAILED)
		munmap(m_pAllocSeqNames,m_AllocMemSeqNames);
#endif
	m_pAllocSeqNames = nullptr;
	}

if(m_pAllocSeqNameIDsOfs != nullptr)
	{
#ifdef _WIN32
	free(m_pAllocSeqNameIDsOfs);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pAllocSeqNameIDsOfs != MAP_FAILED)
		munmap(m_pAllocSeqNameIDsOfs,m_AllocMemSeqNameIDsOfs);
#endif
	m_pAllocSeqNameIDsOfs = nullptr;
	}

if(m_pAllocAlignLoci != nullptr)
	{
#ifdef _WIN32
	free(m_pAllocAlignLoci);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pAllocAlignLoci != MAP_FAILED)
		munmap(m_pAllocAlignLoci,m_AllocMemAlignLoci);
#endif
	m_pAllocAlignLoci = nullptr;
	}

if(m_pSeqNameHashArray != nullptr)
	{
	delete m_pSeqNameHashArray;
	m_pSeqNameHashArray = nullptr;
	}

if(m_hOutFile != -1)
	{
	close(m_hOutFile);
	m_hOutFile = -1;
	}

if(m_pszBuff != nullptr)
	{
	delete m_pszBuff;
	m_pszBuff = nullptr;
	}

if(m_hBEDOutFile != -1)
	{
	close(m_hBEDOutFile);
	m_hBEDOutFile = -1;
	}

if(m_pszBEDBuff != nullptr)
	{
	delete m_pszBEDBuff;
	m_pszBEDBuff = nullptr;
	}

m_BuffIdx = 0;
m_BEDBuffIdx = 0;

DeleteMutexes();
m_bMutexesCreated = false;
m_NumSpecies = 0; 
m_RefSpeciesID = 0;
m_NumSeqNames = 0;		
m_UsedMemSeqNames = 0;		
m_AllocMemSeqNames = 0;	
m_UsedNameHashArray = 0;	
m_AllocMemSeqNameIDsOfs = 0;
m_UsedAlignLoci = 0;	
m_AllocAlignLoci = 0;	
m_AllocMemAlignLoci = 0;	
memset(m_Species,0,sizeof(m_Species));
m_pCurSpecies = nullptr;
m_szCurSeqName[0] = '\0';
m_CurSeqNameID = 0;
m_bSorted = false; 
m_LSER = cDfltLSER;
m_NumThreads = 1;
}

uint16_t											// returned species identifier (1..cMaxSpecies)
CMarkers::AddSpecies(char *pszSpecies,bool IsRefSpecies)	// cultivar or species name
{
int SpeciesID;
tsSNPSSpecies *pSpecies;

// linear search should be ok as normally only expect a few species involved in marker processing
if(m_NumSpecies > 0)
	{
	if(m_pCurSpecies != nullptr && !stricmp((char *)m_pCurSpecies->szSpecies,pszSpecies))
		{
		if(IsRefSpecies && m_RefSpeciesID != m_pCurSpecies->SpeciesID)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddSpecies: Attempted to add multiple reference species (%s)",pszSpecies);
			return(0);
			}
		return(m_pCurSpecies->SpeciesID);
		}
	pSpecies = &m_Species[0];
	for(SpeciesID = 0; SpeciesID < m_NumSpecies; SpeciesID++,pSpecies++)
		{
		if(!stricmp((char *)pSpecies->szSpecies,pszSpecies))
			{
			if(IsRefSpecies && m_RefSpeciesID != pSpecies->SpeciesID)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddSpecies: Attempted to add multiple reference species (%s)",pszSpecies);
				return(0);
				}
			m_pCurSpecies = pSpecies;
			return(pSpecies->SpeciesID);
			}
		}
	}

if(m_NumSpecies == cMaxMarkerSpecies+1)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddSpecies: Attempted to add more (%s) than the max (%d) supported species names",pszSpecies,cMaxMarkerSpecies);
	return(0);
	}
pSpecies = &m_Species[m_NumSpecies++];
pSpecies->SpeciesID = m_NumSpecies;
pSpecies->IsRefSpecies = IsRefSpecies ? 1 : 0;
if(IsRefSpecies)
	m_RefSpeciesID = pSpecies->SpeciesID;
strncpy((char *)pSpecies->szSpecies,pszSpecies,sizeof(pSpecies->szSpecies)-1);
pSpecies->szSpecies[sizeof(pSpecies->szSpecies)-1] = '\0';
m_pCurSpecies = pSpecies;
return(pSpecies->SpeciesID);
}

char *
CMarkers::SpeciesIDtoName(uint16_t SpeciesID)
{
if(SpeciesID < 1 || SpeciesID > m_NumSpecies)
	return(nullptr);
return((char *)&m_Species[SpeciesID-1].szSpecies[0]);
}

bool
CMarkers::IsRefSpecies(uint16_t SpeciesID)
{
if(SpeciesID < 1 || SpeciesID > m_NumSpecies || m_RefSpeciesID < 1)
	return(false);
return(m_RefSpeciesID == SpeciesID);
}

uint16_t					// returned species identifier for specified name, returns 0 if name not previously added with AddSpecies)
CMarkers::NameToSpeciesID(char *pszSpecies)
{
int SpeciesID;
tsSNPSSpecies *pSpecies;

if(m_NumSpecies == 0)
	return(0);

// linear search should be ok as normally only expect a few species involved in marker processing
if(m_pCurSpecies != nullptr && !stricmp((char *)m_pCurSpecies->szSpecies,pszSpecies))
	return(m_pCurSpecies->SpeciesID);

pSpecies = &m_Species[0];
for(SpeciesID = 0; SpeciesID < m_NumSpecies; SpeciesID++,pSpecies++)
	{
	if(!stricmp((char *)pSpecies->szSpecies,pszSpecies))
		{
		m_pCurSpecies = pSpecies;
		return(pSpecies->SpeciesID);
		}
	}
return(0);
}

char *								// returned sequence name
CMarkers::SeqIDtoName(uint32_t SeqID)	// sequence identifier (1..m_NumSeqNames) for which name is to be returned
{
tsSeqName *pSeqName;
uint8_t *pSeqNames;
uint64_t Ofs;
if(SeqID < 1 || SeqID > m_NumSeqNames)
	return(nullptr);
Ofs = m_pAllocSeqNameIDsOfs[SeqID-1];
pSeqNames = (uint8_t *)m_pAllocSeqNames;
pSeqName = (tsSeqName *)&pSeqNames[Ofs];
return((char *)pSeqName->szSeqName);
}

uint32_t 
CMarkers::NameToSeqID(char *pszSeqName) // returned sequence identifier for specified name, returns 0 if name not previously added with AddTargSeq)
{
uint16_t Hash;
uint64_t SeqNameOfs;
uint64_t NxtSeqOfs;
tsSeqName *pSeqName;

if(m_pSeqNameHashArray == nullptr || m_pAllocSeqNames == nullptr || m_NumSeqNames == 0)
	return(0);
// may have struck it lucky and sequence name same as previously added ...
if(m_szCurSeqName[0] != '\0' && !stricmp(pszSeqName,(char *)m_szCurSeqName))
	return(m_CurSeqNameID);
// hash sequence name and use as index into hash array
Hash = CUtility::GenHash16(pszSeqName);
SeqNameOfs = m_pSeqNameHashArray[Hash];	// SeqNameOfs will be 0 if this is a new hash not previously processed

uint8_t *pSeqNames = (uint8_t *)m_pAllocSeqNames;	// used as a convenience instead of requiring many casts when subsequently deriving pSeqName

if(SeqNameOfs == 0)	// no sequence name with this hash previously added?
	return(0);

// seen same hash previously, iterate along names with same hash and check to see if name already known
NxtSeqOfs = SeqNameOfs;
while(SeqNameOfs != 0)
	{
	pSeqName = (tsSeqName *)&pSeqNames[SeqNameOfs - 1];
	if(!stricmp((char *)pSeqName->szSeqName,pszSeqName))
		{
		strcpy((char *)m_szCurSeqName,pszSeqName);
		m_CurSeqNameID = pSeqName->SeqID;
		return(m_CurSeqNameID);
		}
	SeqNameOfs = pSeqName->NxtSeqOfs;
	}
return(0);
}

uint32_t										// returned sequence identifier (1..cMaxSeqID)
CMarkers::AddTargSeq(char *pszSeqName)	// sequence name - could be a chrom, contig, transcript name
{
int SeqNameLen;
uint16_t Hash;
uint64_t SeqNameOfs;
uint64_t NxtSeqOfs;
tsSeqName *pSeqName;

// allocate for 16bit name hashes
if(m_pSeqNameHashArray == nullptr)
	{
	if((m_pSeqNameHashArray = new uint64_t [0x010001])==nullptr)	// allowing for case whereby hash was generated exceeding 16bits!!!!!
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddTargSeq: Memory allocation of %zd bytes - %s",(int64_t)0x010001 * sizeof(uint64_t),strerror(errno));
		return(eBSFerrMem);
		}
	m_UsedNameHashArray = 0;
	memset((size_t *)m_pSeqNameHashArray,0,(size_t)(sizeof(uint64_t) * 0x010001));
	}

if(m_pAllocSeqNames == nullptr)		// initial allocation?
	{
	size_t memreq = cAllocMemSeqNames;

#ifdef _WIN32
	m_pAllocSeqNames = (tsSeqName *) malloc(memreq);	// initial and perhaps the only allocation
	if(m_pAllocSeqNames == nullptr)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddTargSeq: Memory allocation of %zd bytes - %s",(int64_t)memreq,strerror(errno));
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	m_pAllocSeqNames = (tsSeqName *)mmap(nullptr,memreq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pAllocSeqNames == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddTargSeq: Memory allocation of %zd bytes through mmap()  failed - %s",(int64_t)memreq,strerror(errno));
		m_pAllocSeqNames = nullptr;
		return(eBSFerrMem);
		}
#endif
	m_AllocMemSeqNames = (uint64_t)cAllocMemSeqNames;
	m_UsedMemSeqNames = 0;
	m_NumSeqNames = 0;
	}
else
	{
	if(m_AllocMemSeqNames <= (cAllocMinDiffSeqNames + m_UsedMemSeqNames))	// play safe and increase allocation?
		{
		size_t memreq;
		memreq = cAllocMemSeqNames + (size_t)m_AllocMemSeqNames;
#ifdef _WIN32
		pSeqName = (tsSeqName *) realloc(m_pAllocSeqNames,memreq);
#else
		pSeqName = (tsSeqName *)mremap(m_pAllocSeqNames,m_AllocMemSeqNames,memreq,MREMAP_MAYMOVE);
		if(pSeqName == MAP_FAILED)
			pSeqName = nullptr;
#endif
		if(pSeqName == nullptr)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddTargSeq: Sequence names memory re-allocation to %zd from %zd bytes - %s",(int64_t)memreq,m_AllocMemSeqNames,strerror(errno));
			return(eBSFerrMem);
			}
		m_pAllocSeqNames = pSeqName;
		m_AllocMemSeqNames = memreq;
		}
	}


if(m_pAllocSeqNameIDsOfs == nullptr)		// initial allocation?
	{
	size_t memreq = cAllocSeqNames * sizeof(uint64_t);

#ifdef _WIN32
	m_pAllocSeqNameIDsOfs = (uint64_t *) malloc(memreq);	// initial and perhaps the only allocation
	if(m_pAllocSeqNameIDsOfs == nullptr)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddTargSeq: Memory allocation to %zd bytes - %s",(int64_t)memreq,strerror(errno));
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	m_pAllocSeqNameIDsOfs = (uint64_t *)mmap(nullptr,memreq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pAllocSeqNameIDsOfs == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddTargSeq: Memory allocation of %zd bytes through mmap()  failed - %s",(int64_t)memreq,strerror(errno));
		m_pAllocSeqNameIDsOfs = nullptr;
		return(eBSFerrMem);
		}
#endif
	m_AllocMemSeqNameIDsOfs = (uint64_t)memreq;
	}
else
	{
	if(m_AllocMemSeqNameIDsOfs <= ((m_NumSeqNames + 10) * sizeof(uint64_t)))	// play safe and increase allocation?
		{
		uint64_t *pRealloc;
		size_t memreq;
		memreq = (size_t)(m_AllocMemSeqNames + (cAllocSeqNames * sizeof(uint64_t)));
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"AddTargSeq: memory re-allocation to %zd from %zd bytes",(int64_t)memreq,m_AllocMemSeqNameIDsOfs);

#ifdef _WIN32
		pRealloc = (uint64_t *) realloc(m_pAllocSeqNameIDsOfs,memreq);
#else
		pRealloc = (uint64_t *)mremap(m_pAllocSeqNameIDsOfs,m_AllocMemSeqNameIDsOfs,memreq,MREMAP_MAYMOVE);
		if(pRealloc == MAP_FAILED)
			pRealloc = nullptr;
#endif
		if(pRealloc == nullptr)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddTargSeq: Sequence name idex memory re-allocation to %zd from %zd bytes - %s",(int64_t)memreq,m_AllocMemSeqNameIDsOfs,strerror(errno));
			return(eBSFerrMem);
			}
		m_pAllocSeqNameIDsOfs = (uint64_t *)pRealloc;
		m_AllocMemSeqNameIDsOfs = memreq;
		}
	}

// may have struck it lucky and sequence name same as previously added ...
if(m_szCurSeqName[0] != '\0' && !stricmp(pszSeqName,(char *)m_szCurSeqName))
	return(m_CurSeqNameID);

// hash sequence name and use as index into hash array
Hash = CUtility::GenHash16(pszSeqName);
SeqNameOfs = m_pSeqNameHashArray[Hash];	// SeqNameOfs will be 0 if this is a new hash not previously processed

uint8_t *pSeqNames = (uint8_t *)m_pAllocSeqNames;	// used as a convenience instead of requiring many casts when subsequently deriving pSeqName

if(SeqNameOfs == 0)	// no sequence name with this hash previously added?
	{
	m_UsedNameHashArray += 1;
	NxtSeqOfs = 0;
	}
else
	{
	// seen same hash previously, iterate along names with same hash and check to see if name already known
	NxtSeqOfs = SeqNameOfs;
	while(SeqNameOfs != 0)
		{
		pSeqName = (tsSeqName *)&pSeqNames[SeqNameOfs - 1];
		if(!stricmp((char *)pSeqName->szSeqName,pszSeqName))
			{
			strcpy((char *)m_szCurSeqName,pszSeqName);
			m_CurSeqNameID = pSeqName->SeqID;
			return(m_CurSeqNameID);
			}
		SeqNameOfs = pSeqName->NxtSeqOfs;
		}
	}
pSeqName = (tsSeqName *)&pSeqNames[m_UsedMemSeqNames];
SeqNameLen = (int)min(cMaxLenName-1,(int)strlen(pszSeqName));
pSeqName->Len = (uint8_t)(sizeof(tsSeqName) + SeqNameLen);
strncpy((char *)pSeqName->szSeqName,pszSeqName,SeqNameLen);
pSeqName->szSeqName[SeqNameLen] = '\0';
pSeqName->SeqID = m_NumSeqNames + 1;
pSeqName->NxtSeqOfs = NxtSeqOfs;
m_pSeqNameHashArray[Hash] = m_UsedMemSeqNames + 1;
m_pAllocSeqNameIDsOfs[m_NumSeqNames] = m_UsedMemSeqNames;
m_UsedMemSeqNames += pSeqName->Len;
m_NumSeqNames += 1;
return(pSeqName->SeqID);
}

int
CMarkers::PreAllocEstSNPs(int64_t EstNumSNPS)	// estimating will be required to process this many SNP loci
{
size_t memreq; 
int64_t AllocdAlignLoci;

AllocdAlignLoci = ((99 + EstNumSNPS) * (int64_t)125) / (int64_t)100;  // allowing for an extra 25% to reduce probability of a realloc being subsequently required if estimate was a little low
memreq = (size_t)AllocdAlignLoci * sizeof(tsAlignLoci); 
#ifdef _WIN32
m_pAllocAlignLoci = (tsAlignLoci *) malloc(memreq);	// initial and perhaps the only allocation
if(m_pAllocAlignLoci == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"PreAllocSNPs: Memory allocation of %zd bytes - %s",(int64_t)memreq,strerror(errno));
	return(eBSFerrMem);
	}
#else
// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
m_pAllocAlignLoci = (tsAlignLoci *)mmap(nullptr,memreq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
if(m_pAllocAlignLoci == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"PreAllocSNPs: Memory allocation of %zd bytes through mmap()  failed - %s",(int64_t)memreq,strerror(errno));
	m_pAllocAlignLoci = nullptr;
	return(eBSFerrMem);
	}
#endif
m_AllocMemAlignLoci = memreq;
m_UsedAlignLoci = 0;
m_AllocAlignLoci = AllocdAlignLoci;
memset(m_pAllocAlignLoci,0,memreq);
return(eBSFSuccess);
}

int
CMarkers::PreAllocImputedSNPs(int NumbIsolates)	// calculate the number of actually required SNP loci from number of reported SNPs and number of isolates
{
int64_t NumRefSNPs;
int64_t TotalRefImputed;
tsAlignLoci *pCurLoci;
uint32_t TargLoci;
uint32_t TargSeqID;
uint32_t TargSpeciesID;


if(m_pAllocAlignLoci == nullptr || m_UsedAlignLoci < 1) // must be at least 1 known SNP loci before any can be imputed!
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "PreAllocImputedSNPs: need at least one known SNP loci before any can be imputed!");
	return(eBSFerrParams);
	}

// ensure sorting of known SNPs has been completed as will be iterating over these in order to calculate additional imputed SNP loci structure instances required
SortTargSeqLociSpecies();
// determine number of distinct SNP loci
pCurLoci = m_pAllocAlignLoci;
TargLoci = pCurLoci->TargLoci;
TargSeqID = pCurLoci->TargSeqID;
TargSpeciesID = pCurLoci->TargSpeciesID;
NumRefSNPs = 1;
for(int64_t Idx = 0; Idx < m_UsedAlignLoci; Idx++, pCurLoci++)
	if(pCurLoci->TargLoci != TargLoci || TargSeqID != pCurLoci->TargSeqID || TargSpeciesID != pCurLoci->TargSpeciesID)
		{
		TargLoci = pCurLoci->TargLoci;
		TargSeqID = pCurLoci->TargSeqID;
		TargSpeciesID = pCurLoci->TargSpeciesID;
		NumRefSNPs++;
		}
TotalRefImputed = (NumRefSNPs * NumbIsolates);

if ((TotalRefImputed + 100) > m_AllocAlignLoci)	// play safe when increasing allocation
	{
	size_t memreq;
	memreq = (TotalRefImputed+100) * sizeof(tsAlignLoci);
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Memory re-allocation to %zd from %zd bytes for estimated %zd SNPs", (int64_t)memreq, (int64_t)m_AllocMemAlignLoci, (int64_t)TotalRefImputed);

#ifdef _WIN32
	pCurLoci = (tsAlignLoci*)realloc(m_pAllocAlignLoci, memreq);
#else
	pCurLoci = (tsAlignLoci*)mremap(m_pAllocAlignLoci, m_AllocMemAlignLoci, memreq, MREMAP_MAYMOVE);
	if (pCurLoci == MAP_FAILED)
		pCurLoci = nullptr;
#endif
	if (pCurLoci == nullptr)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "PreAllocImpunedSNPs: Memory re-allocation to %zd bytes - %s", (int64_t)memreq, strerror(errno));
		return(eBSFerrMem);
		}
	m_pAllocAlignLoci = pCurLoci;
	m_AllocMemAlignLoci = memreq;
	m_AllocAlignLoci = TotalRefImputed + 100;
	}
return(eBSFSuccess);
}

int64_t 
CMarkers::NumAlignLoci(void)					// returns current number of alignment/SNP loci
{
return(m_UsedAlignLoci);
}

int64_t 
CMarkers::AddLoci(uint16_t TargSpeciesID,	// reads were aligned to this cultivar or species
				uint32_t TargSeqID,		// alignments to this sequence - could be a chrom/contig/transcript - from pszSpecies
				uint32_t TargLoci,		// loci within target sequence at which SNPs observed
				etSeqBase TargRefBase,	// loci has this reference base
				uint16_t ProbeSpeciesID,	// reads were aligned from this cultivar or species
				uint32_t ProbeCntA,		// number instances probe base A aligned to TargRefBase 
				uint32_t ProbeCntC,		// number instances probe base C aligned to TargRefBase
				uint32_t ProbeCntG,		// number instances probe base G aligned to TargRefBase
				uint32_t ProbeCntT,		// number instances probe base T aligned to TargRefBase
				uint32_t ProbeCntN,		// number instances probe base U aligned to TargRefBase
				double LSER,			// local sequencing error rate
				uint16_t Flags)			// any loci associated flags
{
tsAlignLoci *pLoci;

if(m_pAllocAlignLoci == nullptr)		// initial allocation?
	{
	size_t memreq = cAllocAlignLoci * sizeof(tsAlignLoci);

#ifdef _WIN32
	m_pAllocAlignLoci = (tsAlignLoci *) malloc(memreq);	// initial and perhaps the only allocation
	if(m_pAllocAlignLoci == nullptr)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddLoci: Memory allocation of %zd bytes - %s",(int64_t)memreq,strerror(errno));
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	m_pAllocAlignLoci = (tsAlignLoci *)mmap(nullptr,memreq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pAllocAlignLoci == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddLoci: Memory allocation of %zd bytes through mmap()  failed - %s",(int64_t)memreq,strerror(errno));
		m_pAllocAlignLoci = nullptr;
		return(eBSFerrMem);
		}
#endif
	m_AllocMemAlignLoci = memreq;
	m_UsedAlignLoci = 0;
	m_AllocAlignLoci = cAllocAlignLoci;
	}
else
	{
	if(m_AllocAlignLoci <= (m_UsedAlignLoci + 100))	// play safe when increasing allocation
		{
		size_t memreq;
		size_t AllocTo;
		AllocTo = ((size_t)cReAllocAlignPerc * m_AllocAlignLoci)/100;
		memreq = ((AllocTo + 1) * sizeof(tsAlignLoci));
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"AddLoci: memory re-allocation to %zd from %zd bytes",(int64_t)memreq,(int64_t)m_AllocMemAlignLoci);

#ifdef _WIN32
		pLoci = (tsAlignLoci *) realloc(m_pAllocAlignLoci,memreq);
#else
		pLoci = (tsAlignLoci *)mremap(m_pAllocAlignLoci,m_AllocMemAlignLoci,memreq,MREMAP_MAYMOVE);
		if(pLoci == MAP_FAILED)
			pLoci = nullptr;
#endif
		if(pLoci == nullptr)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddLoci: Memory re-allocation to %zd bytes - %s",(int64_t)memreq,strerror(errno));
			return(eBSFerrMem);
			}
		m_pAllocAlignLoci = pLoci;
		m_AllocMemAlignLoci = memreq;
		m_AllocAlignLoci = AllocTo; 
		}
	}

pLoci = &m_pAllocAlignLoci[m_UsedAlignLoci++];
pLoci->AlignID = m_UsedAlignLoci;
pLoci->TargSpeciesID = TargSpeciesID;
pLoci->TargSeqID = TargSeqID;
pLoci->TargLoci = TargLoci;
pLoci->TargRefBase = TargRefBase;
pLoci->Flags = Flags;
pLoci->ProbeSpeciesID = ProbeSpeciesID;
pLoci->ProbeBaseCnts[0] = ProbeCntA;
pLoci->ProbeBaseCnts[1] = ProbeCntC;
pLoci->ProbeBaseCnts[2] = ProbeCntG;
pLoci->ProbeBaseCnts[3] = ProbeCntT;
pLoci->ProbeBaseCnts[4] = ProbeCntN;
pLoci->LSER = LSER;
m_bSorted = false;
return(m_UsedAlignLoci);
}

// AddLoci
// Add loci on target which has identified SNP when aligned to from probe sequences
int64_t  
CMarkers::AddLoci(char *pszTargSpecies,	// reads were aligned to this cultivar or species
				char *pszTargSeq,		// alignments to this sequence - could be a chrom/contig/transcript - from pszSpecies
				uint32_t TargLoci,			// loci within target sequence at which SNPs observed
				etSeqBase TargRefBase,	// loci has this reference base
				char *pszProbeSpecies,	// reads were aligned from this cultivar or species
				uint32_t ProbeCntA,		// number instances probe base A aligned to TargRefBase 
				uint32_t ProbeCntC,		// number instances probe base C aligned to TargRefBase
				uint32_t ProbeCntG,		// number instances probe base G aligned to TargRefBase
				uint32_t ProbeCntT,		// number instances probe base T aligned to TargRefBase
				uint32_t ProbeCntN,		// number instances probe base U aligned to TargRefBase
				double LSER)				// local sequencing error rate
{
uint16_t TargID;
uint16_t ProbeID;
uint32_t TargSeqID;
int64_t Rslt;

TargID = AddSpecies(pszTargSpecies,true);
if(TargID == 0 || TargID != m_RefSpeciesID)
	return(-1);

ProbeID = AddSpecies(pszProbeSpecies,false);
if(ProbeID == 0 || ProbeID == m_RefSpeciesID)
	return(-1);

TargSeqID = AddTargSeq(pszTargSeq);
if(TargSeqID == 0)
	return(-1);

Rslt = AddLoci(TargID,TargSeqID,TargLoci,TargRefBase,ProbeID,ProbeCntA,ProbeCntC,ProbeCntG,ProbeCntT,ProbeCntN, LSER,cFlgSNPcnts);
return(Rslt);
}

int 
CMarkers::LoadSNPFile(int MinBases,			// accept SNPs with at least this number covering bases
					  double MaxPValue,		// accept SNPs with at most this P-value
					  double SNPMmajorPC,	// only accept SNP for processing if major allele >= this proportion of total allele counts
					  char *pszRefSpecies,	// this is the reference species 
					  char *pszProbeSpecies, // this species reads were aligned to the reference species from which SNPs were called 
					  char *pszSNPFile)		// SNP file to parse and load
{
int Rslt;
int64_t Rslt64;
int NumFields;
int NumElsParsed;

char *pszRefSeq;
etSeqBase RefBase;
int StartLoci;
int Mismatches;
char *pszRefBase;
int BaseCnt[5];
double LSER;

int CoveringBases;
int NumFilteredOut;
int FilteredCovBases;
int FilteredPValue;
int FilteredMajor;

CCSVFile *pCSV = new CCSVFile;
if(pCSV == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CCSVfile");
	return(eBSFerrObj);
	}

if((Rslt=pCSV->Open(pszSNPFile))!=eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open SNP file: %s",pszSNPFile);
	delete pCSV;
	return(Rslt);
	}

// fix for reported issue whereby all SNPs were filtered out from this SNP file
// and thus probe or query was not being registered
uint16_t TargID;
uint16_t ProbeID;
uint32_t TargSeqID;

TargID = AddSpecies(pszRefSpecies, true);
if (TargID == 0 || TargID != m_RefSpeciesID)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to register species: %s", pszRefSpecies);
	delete pCSV;
	return(eBSFerrFeature);
	}

ProbeID = AddSpecies(pszProbeSpecies, false);
if (ProbeID == 0 || ProbeID == m_RefSpeciesID)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to register species: %s", pszProbeSpecies);
	delete pCSV;
	return(eBSFerrFeature);
	}

double PValue;

NumElsParsed = 0;
NumFilteredOut = 0;
FilteredCovBases = 0;
FilteredPValue = 0;
FilteredMajor = 0;

while((Rslt=pCSV->NextLine()) > 0)	// onto next line containing fields
	{
	NumFields = pCSV->GetCurFields();
	if(!NumElsParsed && (NumFields >= 23 && pCSV->IsLikelyHeaderLine())) // expecting at least 23 fields so only then check for header line
		continue;

	NumElsParsed += 1;

	pCSV->GetText(4,&pszRefSeq);
	pCSV->GetInt(5,&StartLoci);
	pCSV->GetInt(12,&Mismatches);
	pCSV->GetText(13,&pszRefBase);
	pCSV->GetInt(14,&BaseCnt[0]);
	pCSV->GetInt(15,&BaseCnt[1]);
	pCSV->GetInt(16,&BaseCnt[2]);
	pCSV->GetInt(17,&BaseCnt[3]);
	pCSV->GetInt(18,&BaseCnt[4]);

	TargSeqID = AddTargSeq(pszRefSeq);
	if (TargSeqID == 0)
		continue;

	// apply any filtering
	pCSV->GetInt(11, &CoveringBases);
	if (CoveringBases < MinBases)
		{
		NumFilteredOut += 1;
		FilteredCovBases += 1;
		continue;
		}

	pCSV->GetDouble(10, &PValue);
	if (PValue > MaxPValue)
		{
		NumFilteredOut += 1;
		FilteredPValue += 1;
		continue;
		}

	if((double)Mismatches/CoveringBases < SNPMmajorPC)
		{
		NumFilteredOut += 1;
		FilteredMajor += 1;
		continue;
		}

	pCSV->GetDouble(19, &LSER);
	if(LSER < cFloorLSER)		// use a floor to prevent potential div errors later in processing
		LSER = cFloorLSER;

	switch(pszRefBase[0]) {
		case 'a': case 'A':
			RefBase = eBaseA;
			BaseCnt[0] = CoveringBases - Mismatches;
			break;
		case 'c': case 'C':
			RefBase = eBaseC;
			BaseCnt[1] = CoveringBases - Mismatches;
			break;
		case 'g': case 'G':
			RefBase = eBaseG;
			BaseCnt[2] = CoveringBases - Mismatches;
			break;
		case 't': case 'T':
			RefBase = eBaseT;
			BaseCnt[3] = CoveringBases - Mismatches;
			break;
		default:
			RefBase = eBaseN;
			BaseCnt[4] = CoveringBases - Mismatches;
			break;
		}

	int BaseIdx;
	for(BaseIdx = 0; BaseIdx < 5; BaseIdx++)
		{
		if(BaseIdx == RefBase)
			continue;
		if(((double)BaseCnt[BaseIdx] + 1)/CoveringBases >= SNPMmajorPC)
			break;
		}
	if(BaseIdx > 4)
		{
		NumFilteredOut += 1;
		FilteredMajor += 1;
		continue;
		}

	Rslt64 = AddLoci(pszRefSpecies,		// reads were aligned to this cultivar or species
				pszRefSeq,		// alignments to this sequence - could be a chrom/contig/transcript - from pszSpecies
				StartLoci,		// loci within target sequence at which SNPs observed
				RefBase,		// loci has this reference base
				pszProbeSpecies,// reads were aligned from this cultivar or species
				BaseCnt[0],		// number instances probe base A aligned to TargRefBase 
				BaseCnt[1],		// number instances probe base C aligned to TargRefBase
				BaseCnt[2],		// number instances probe base G aligned to TargRefBase
				BaseCnt[3],		// number instances probe base T aligned to TargRefBase
				BaseCnt[4],		// number instances probe base U aligned to TargRefBase
				LSER);			// local sequencing error rate

	if(Rslt64 < 1)
		{
		if(pCSV != nullptr)
			delete pCSV;
		return(-1);
		}
	}
if(pCSV != nullptr)
	delete pCSV;
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Parsed %d SNPs from file: %s",NumElsParsed,pszSNPFile);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Accepted %d SNPs, filtered out - %d high P-Value, %d low coverage, %d major allele",
					NumElsParsed - NumFilteredOut, FilteredPValue, FilteredCovBases, FilteredMajor);
m_NumSSNPLoci = m_UsedAlignLoci;
return(NumElsParsed - NumFilteredOut);
}

// 
// AddImputedAlignments
// Add alignments for species where no snp was called but other species do have snp called
// The call could be because there were none or insufficient reads covering the loci, or there was coverage but no snp!
int64_t 
CMarkers::AddImputedAlignments(int MinBases,			// must be at least this number of reads covering the SNP loci
					  char *pszRefSpecies,				// this is the reference species 
					char *pszProbeSpecies,				// this species reads were aligned to the reference species from which SNPs were called 
					char *pszAlignFile,					// file containing alignments
					int FType,							// input alignment file format: 0 - auto, 1 - CSV, 2 - BED, 3 - SAM)
					bool bSeqs,							// if alignment file contains the read sequence then impute bases from the actual sequences	
					int EstNumSeqs,						// estimated number of sequences (0 if no estimate)
					int EstSeqLen)						// estimated mean sequence length (0 if no estimate)           			
{
int Rslt;
int64_t Rslt64;
int NumEls;
uint16_t RefSpeciesID;
uint16_t ProbeSpeciesID;
int MinLength = 50;
int MaxLength = 1000;
uint16_t ImputFlags;
CHyperEls *pHypers;

if((ProbeSpeciesID = NameToSpeciesID(pszProbeSpecies)) < 1)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to locate identifier for probe species '%s'",pszProbeSpecies);
	return(eBSFerrInternal);
	}
if((RefSpeciesID = NameToSpeciesID(pszRefSpecies)) < 1)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to locate identifier for reference species species '%s'",pszRefSpecies);
	return(eBSFerrInternal);
	}


if((pHypers = new CHyperEls)==nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CHyperEls");
	return(eBSFerrObj);
	}

if(EstNumSeqs != 0)
	if((Rslt = pHypers->PreAllocMem(EstNumSeqs,EstSeqLen)) != eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to prealloc allocate memory for sequences");
		delete pHypers;
		return(eBSFerrObj);
		}

etClassifyFileType FileType;

if(FType == 0)
	FileType = CUtility::ClassifyFileType(pszAlignFile);
else
	FileType = (etClassifyFileType)(FType - 1);

ImputFlags = cFlgImputCnts;
switch(FileType) {
	case eCFTopenerr:		// unable to open file for reading
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open file: '%s'",pszAlignFile);
		delete pHypers;
		return(eBSFerrOpnFile);

	case eCFTlenerr:		// file length is insufficent to classify type
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Unable to classify file type (insuffient data points): '%s'",pszAlignFile);
		delete pHypers; 
		return(eBSFerrFileAccess);

	case eCFTunknown:		// unable to reliably classify
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Unable to reliably classify file type: '%s'",pszAlignFile);
		delete pHypers; 
		return(eBSFerrFileType);

	case eCFTCSV:			// file has been classified as being CSV
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Parsing CSV file (%d..%d): '%s'",MinLength,MaxLength,pszAlignFile);
		if((Rslt = pHypers->ParseCSVFileElements(pszAlignFile,MinLength,MaxLength,eCSVFdefault)) < 0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Parse errors in CSV file (%d..%d): '%s'",MinLength,MaxLength,pszAlignFile);
			delete pHypers;
			return(Rslt);
			}
		break;

	case eCFTBED:			// file has been classified as being BED
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Parsing BED file (%d..%d): '%s'",MinLength,MaxLength,pszAlignFile);
		if((Rslt = pHypers->ParseBEDFileElements(pszAlignFile,MinLength,MaxLength)) < 0)
			{	
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Parse errors in BED file (%d..%d): '%s'",MinLength,MaxLength,pszAlignFile);
			delete pHypers; 
			return(Rslt);
			}
		break;

	case eCFTSAM:			// file has been classified as being SAM
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Parsing SAM file (%d..%d): '%s'",MinLength,MaxLength,pszAlignFile);
		if((Rslt = pHypers->ParseSAMFileElements(pszAlignFile,MinLength,MaxLength,bSeqs)) < 0)
			{	
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Processing errors in SAM/BAM file (%d..%d): '%s'",MinLength,MaxLength,pszAlignFile);
			delete pHypers; 
			return(Rslt);
			}
		ImputFlags = cFlgAlignCnts;
		break;

	default:
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to classify file type: '%s'",pszAlignFile);
		delete pHypers; 
		return(eBSFerrFileType);
	}


NumEls = pHypers->NumEls();
if(NumEls == 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"No elements with length range %d..%d in file: '%s'",MinLength,MaxLength,pszAlignFile);
	delete pHypers;
	return(Rslt);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loaded and parsed %d elements",NumEls);

// try to find overlaps on to SNPs called on other species where this probe species is not represented!
bool bProbeAligned;
int64_t AlignIdx;
int64_t UsedAlignLoci;
tsAlignLoci *pAlign;
uint32_t CurTargSeqID;
uint32_t CurTargLociLoci;
uint32_t PrevTargSeqID;
int NumOverlapping;
int HyperChromID;
char *pszTargSeq;
etSeqBase TargRefBase;

int64_t TotNonOverlapping;
int64_t TotOverlapping;

uint32_t ProbeCntA;		// number instances probe base A aligned to TargRefBase 
uint32_t ProbeCntC;		// number instances probe base C aligned to TargRefBase
uint32_t ProbeCntG;		// number instances probe base G aligned to TargRefBase
uint32_t ProbeCntT;		// number instances probe base T aligned to TargRefBase
uint32_t ProbeCntN;		// number instances probe base U aligned to TargRefBase

TotOverlapping = 0;
TotNonOverlapping = 0;
AlignIdx = 0;
bProbeAligned = false;
PrevTargSeqID = 0;
CurTargSeqID = 0;
pszTargSeq = nullptr;
UsedAlignLoci = m_NumSSNPLoci;
for(AlignIdx = 0; AlignIdx < UsedAlignLoci; AlignIdx++)
	{
	pAlign = &m_pAllocAlignLoci[AlignIdx];			// m_pAllocAlignLoci could be realloc'd so best to take the address each iteration....
	if(AlignIdx == 0)
		{
		CurTargSeqID = pAlign->TargSeqID;
		CurTargLociLoci = pAlign->TargLoci;
		TargRefBase = pAlign->TargRefBase;
		if(pAlign->ProbeSpeciesID == ProbeSpeciesID)
			bProbeAligned = true;
		else
			bProbeAligned = false;
		continue;
		}


	if(bProbeAligned == false && (CurTargSeqID != pAlign->TargSeqID || CurTargLociLoci != pAlign->TargLoci))
		{
				// no snp called for probe species, check if there were reads aligned to reference
		if(CurTargSeqID != PrevTargSeqID || pszTargSeq == nullptr)
			{
			pszTargSeq = SeqIDtoName(CurTargSeqID);
			HyperChromID = pHypers->GetChromID(pszTargSeq);
			PrevTargSeqID = CurTargSeqID;
			}

		ProbeCntA = 0;
		ProbeCntC = 0;
		ProbeCntG = 0;
		ProbeCntT = 0;
		ProbeCntN = 0;
		if(HyperChromID < 1)		// will be < 1 if no target sequence alignments in alignment file
			NumOverlapping = 0;
		else
			{
			NumOverlapping = pHypers->LocateLociBaseCnts(HyperChromID,CurTargLociLoci,&ProbeCntA,&ProbeCntC,&ProbeCntG,&ProbeCntT,&ProbeCntN);

			if(NumOverlapping == 0 || (ProbeCntA == 0 && ProbeCntC == 0 && ProbeCntG == 0 && ProbeCntT == 0 && ProbeCntN == 0))
				{
				ProbeCntA = 0;
				ProbeCntC = 0;
				ProbeCntG = 0;
				ProbeCntT = 0;
				ProbeCntN = 0;
				switch(TargRefBase) {
					case eBaseA:
						ProbeCntA = NumOverlapping;
						break;
					case eBaseC:
						ProbeCntC = NumOverlapping;
						break;
					case eBaseG:
						ProbeCntG = NumOverlapping;
						break;
					case eBaseT:
						ProbeCntT = NumOverlapping;
						break;
					case eBaseN:
						ProbeCntN = NumOverlapping;
						break;
					}
				}
			}


		if(NumOverlapping == 0)
			TotNonOverlapping += 1;
		else
			TotOverlapping += 1;

		// add number overlapping to a new loci record with refbase count set to NumOverlapping
		Rslt64 = AddLoci(RefSpeciesID,	// reads were aligned to this cultivar or species
				CurTargSeqID,			// alignments to this sequence - could be a chrom/contig/transcript - from pszSpecies
				CurTargLociLoci,		// loci within target sequence at which SNPs observed
				TargRefBase,			// loci has this reference base
				ProbeSpeciesID,			// reads were aligned from this cultivar or species
				ProbeCntA,				// number instances probe base A aligned to TargRefBase 
				ProbeCntC,				// number instances probe base C aligned to TargRefBase
				ProbeCntG,				// number instances probe base G aligned to TargRefBase
				ProbeCntT,				// number instances probe base T aligned to TargRefBase
				ProbeCntN,				// number instances probe base U aligned to TargRefBase
				m_LSER,					// defaulting local sequencing error rate
				ImputFlags);			// user flag to indicate these are imputed counts, not from the SNP file
		if(Rslt64 < 1)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddImputedAlignments: AddLoci returned error %d",(int)Rslt64);
			delete pHypers;
			return(Rslt64);			
			}

		pAlign = &m_pAllocAlignLoci[AlignIdx];
		}

	if(CurTargSeqID != pAlign->TargSeqID || CurTargLociLoci != pAlign->TargLoci)
		{
		CurTargSeqID = pAlign->TargSeqID;
		CurTargLociLoci = pAlign->TargLoci;
		TargRefBase = pAlign->TargRefBase;
		bProbeAligned = false;
		}
	if(pAlign->ProbeSpeciesID == ProbeSpeciesID)
		bProbeAligned = true;
	}

if(!bProbeAligned)
	{
	if(CurTargSeqID != PrevTargSeqID || pszTargSeq == nullptr)
		{
		pszTargSeq = SeqIDtoName(CurTargSeqID);
		HyperChromID = pHypers->GetChromID(pszTargSeq);
		PrevTargSeqID = CurTargSeqID;
		}
	ProbeCntA = 0;
	ProbeCntC = 0;
	ProbeCntG = 0;
	ProbeCntT = 0;
	ProbeCntN = 0;
	if(HyperChromID < 1)
		NumOverlapping = 0;
	else
		{
		NumOverlapping = pHypers->LocateLociBaseCnts(HyperChromID,CurTargLociLoci,&ProbeCntA,&ProbeCntC,&ProbeCntG,&ProbeCntT,&ProbeCntN);

		if(NumOverlapping == 0 || (ProbeCntA == 0 && ProbeCntC == 0 && ProbeCntG == 0 && ProbeCntT == 0 && ProbeCntN == 0))
			{
			ProbeCntA = 0;
			ProbeCntC = 0;
			ProbeCntG = 0;
			ProbeCntT = 0;
			ProbeCntN = 0;
			switch(TargRefBase) {
				case eBaseA:
					ProbeCntA = NumOverlapping;
					break;
				case eBaseC:
					ProbeCntC = NumOverlapping;
					break;
				case eBaseG:
					ProbeCntG = NumOverlapping;
					break;
				case eBaseT:
					ProbeCntT = NumOverlapping;
					break;
				case eBaseN:
					ProbeCntN = NumOverlapping;
					break;
				}
			}
		}
	if(NumOverlapping < MinBases)
		NumOverlapping = 0;
	if(NumOverlapping == 0)
		TotNonOverlapping += 1;
	else
		TotOverlapping += 1;

		// add number overlapping to a new loci record with refbase count set to NumOverlapping
	Rslt64 = AddLoci(RefSpeciesID,	// reads were aligned to this cultivar or species
				CurTargSeqID,			// alignments to this sequence - could be a chrom/contig/transcript - from pszSpecies
				CurTargLociLoci,		// loci within target sequence at which SNPs observed
				TargRefBase,			// loci has this reference base
				ProbeSpeciesID,			// reads were aligned from this cultivar or species
				ProbeCntA,				// number instances probe base A aligned to TargRefBase 
				ProbeCntC,				// number instances probe base C aligned to TargRefBase
				ProbeCntG,				// number instances probe base G aligned to TargRefBase
				ProbeCntT,				// number instances probe base T aligned to TargRefBase
				ProbeCntN,				// number instances probe base U aligned to TargRefBase
				m_LSER,					// use default local sequencing error rate
				ImputFlags);			// user flag to indicate these are imputed counts, not from the SNP file

	if(Rslt64 < 1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddImputedAlignments: AddLoci returned error %d",(int)Rslt64);
		delete pHypers;
		return(Rslt64);
		}	
	pAlign = &m_pAllocAlignLoci[AlignIdx];
	}
delete pHypers;
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Added %d loci with alignments, %d loci with no alignments",TotOverlapping,TotNonOverlapping);
return(TotOverlapping);
}



// 
// AddSimulatedAlignments
// Add simulated alignments for when no actual alignments were available, just SNP calls
int64_t 
CMarkers::AddSimulatedAlignments(int MinBases,			// using this as the simulated number of reads covering the SNP loci
					  char *pszRefSpecies,				// this is the reference species 
					char *pszProbeSpecies)				// this species reads were aligned to the reference species from which SNPs were called 
{
int64_t Rslt64;
uint16_t RefSpeciesID;
uint16_t ProbeSpeciesID;
uint16_t ImputFlags;
int64_t TotSimulated = 0;

if((ProbeSpeciesID = NameToSpeciesID(pszProbeSpecies)) < 1)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to locate identifier for probe species '%s'",pszProbeSpecies);
	return(eBSFerrInternal);
	}
if((RefSpeciesID = NameToSpeciesID(pszRefSpecies)) < 1)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to locate identifier for reference species species '%s'",pszRefSpecies);
	return(eBSFerrInternal);
	}

ImputFlags = cFlgImputCnts;

int64_t AlignIdx;
int64_t UsedAlignLoci;
tsAlignLoci *pAlign;

uint32_t ProbeCntA;		// number instances probe base A aligned to TargRefBase 
uint32_t ProbeCntC;		// number instances probe base C aligned to TargRefBase
uint32_t ProbeCntG;		// number instances probe base G aligned to TargRefBase
uint32_t ProbeCntT;		// number instances probe base T aligned to TargRefBase
uint32_t ProbeCntN;		// number instances probe base U aligned to TargRefBase
uint32_t SimAlignLoci;	
uint32_t SimTargSeqID;
uint8_t SimRefBase;
bool bSNPCalled = false;
AlignIdx = 0;
UsedAlignLoci = m_NumSSNPLoci;
pAlign = &m_pAllocAlignLoci[AlignIdx];
SimAlignLoci = pAlign->TargLoci;
SimTargSeqID = pAlign->TargSeqID;
SimRefBase = pAlign->TargRefBase;
for(AlignIdx = 0; AlignIdx < UsedAlignLoci; AlignIdx++)
	{
	pAlign = &m_pAllocAlignLoci[AlignIdx];			// m_pAllocAlignLoci could be realloc'd so best to take the address each iteration....
	if(pAlign->TargLoci == SimAlignLoci && SimTargSeqID == pAlign->TargSeqID)			// simulate alignment just once for each SNP loci on a target sequence
		{
		if(pAlign->ProbeSpeciesID == ProbeSpeciesID)	// already actually called a SNP for probe species?
			{
			bSNPCalled = true;
			continue;
			}
		if(AlignIdx != UsedAlignLoci-1)
			continue;
		}

	// target loci changed
	if(bSNPCalled)			// already called a SNP at previous alignment loci so don't need to simulate
		{
		SimAlignLoci = pAlign->TargLoci;
		SimTargSeqID = pAlign->TargSeqID;
		SimRefBase = pAlign->TargRefBase;
		AlignIdx--;
		bSNPCalled = false;
		continue;
		}

	// didn't call a SNP at SimAlignLoci so need to simulate
	ProbeCntA = 0;				// number instances probe base A aligned to TargRefBase 
	ProbeCntC = 0;				// number instances probe base C aligned to TargRefBase
	ProbeCntG = 0;				// number instances probe base G aligned to TargRefBase
	ProbeCntT = 0;				// number instances probe base T aligned to TargRefBase
	ProbeCntN = 0;				// number instances probe base U aligned to TargRefBase
	switch(SimRefBase) {
		case eBaseA:
			ProbeCntA = MinBases;
			break;
		case eBaseC:
			ProbeCntC = MinBases;
			break;
		case eBaseG:
			ProbeCntG = MinBases;
			break;
		case eBaseT:
			ProbeCntT = MinBases;
			break;
		case eBaseN:
			ProbeCntN = MinBases;
			break;
			}

			// add number overlapping to a new loci record with refbase count set to NumOverlapping
	Rslt64 = AddLoci(RefSpeciesID,	// reads were aligned to this cultivar or species
			SimTargSeqID,			// alignments to this sequence - could be a chrom/contig/transcript - from pszSpecies
			SimAlignLoci,			// loci within target sequence at which SNPs observed
			SimRefBase,				// loci has this reference base
			ProbeSpeciesID,			// reads were aligned from this cultivar or species
			ProbeCntA,				// number instances probe base A aligned to TargRefBase 
			ProbeCntC,				// number instances probe base C aligned to TargRefBase
			ProbeCntG,				// number instances probe base G aligned to TargRefBase
			ProbeCntT,				// number instances probe base T aligned to TargRefBase
			ProbeCntN,				// number instances probe base U aligned to TargRefBase
			m_LSER,					// defaulting local sequencing error rate
			ImputFlags);			// user flag to indicate these are imputed counts, not from the SNP file
	if(Rslt64 < 1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddSimulatedAlignments: AddLoci returned error %d",(int)Rslt64);
		return(Rslt64);			
		}
	TotSimulated++;
	if(AlignIdx == UsedAlignLoci-1)
		continue;
	pAlign = &m_pAllocAlignLoci[AlignIdx];			// m_pAllocAlignLoci could be realloc'd so best to take the address each iteration....
	SimAlignLoci = pAlign->TargLoci;
	SimTargSeqID = pAlign->TargSeqID;
	SimRefBase = pAlign->TargRefBase;
	AlignIdx--;
	bSNPCalled = false;
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Added %d simulated reference loci",TotSimulated);
return(TotSimulated);
}

int
CMarkers::IdentSpeciesSpec(int AltMaxCnt,				// max count allowed for base being processed in any other species, 0 if no limit
						int MinCnt,						// min count required for base being processed in species
						double SNPMmajorPC,				// to be processed major putative SNP base must be at least this proportion of total
						int MinSpeciesWithCnts,			// must be at least this number of species with base counts more than MinSpeciesTotCntThres - 0 if no limit 
						int MinSpeciesTotCntThres)		// individual species must have at least this number of total bases at SNP loci - 0 if no threshold

{
int64_t AlignIdx;
tsAlignLoci *pAlign;
tsAlignLoci *pAlignSpecies;
tsAlignLoci *pAlignSpeciesA;
int NumSpecies = m_NumSpecies-1; // note: m_NumSpecies includes reference species name so need to subtract 1 to get number of cultivars/isolates
int SpeciesIdx;
int SpeciesIdxA;
int BaseIdx;
bool bAcceptSpec;
etSeqBase BestAcceptBase;
double BestAcceptConf;
double CurConf;
int NumSpeciesWithCnts;

SortTargSeqLociSpecies();

time_t Then = time(nullptr);
time_t Now;
int64_t NumLociProc = 0;
int64_t PutSNPLoci = m_UsedAlignLoci / (int64_t)(m_NumSpecies - 1);
pAlign = &m_pAllocAlignLoci[0];
for(AlignIdx = 0; AlignIdx < m_UsedAlignLoci; AlignIdx += NumSpecies, pAlign += NumSpecies)
	{
	if ((NumLociProc % 10000) == 0)
		{
		Now = time(nullptr);
		if ((Now - Then) >= 60)
			{
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "Filtering: Processed %zd from %zd putative SNP loci", NumLociProc, PutSNPLoci);
			Then += 60;
			}
		}
	NumLociProc++;
	NumSpeciesWithCnts = 0;
	pAlignSpecies = pAlign;
	for(SpeciesIdx = 0; SpeciesIdx < NumSpecies; SpeciesIdx++,pAlignSpecies += 1)
		{
		pAlignSpecies->TotBases = pAlignSpecies->ProbeBaseCnts[0]+pAlignSpecies->ProbeBaseCnts[1]+pAlignSpecies->ProbeBaseCnts[2]+pAlignSpecies->ProbeBaseCnts[3]+pAlignSpecies->ProbeBaseCnts[4];
		if(pAlignSpecies->TotBases == 0 || (uint32_t)MinSpeciesTotCntThres > pAlignSpecies->TotBases)
			{
			pAlignSpecies->CultSpecBase = eBaseN;
			if(pAlignSpecies->LSER < cFloorLSER)
				pAlignSpecies->LSER = m_LSER;
			pAlignSpecies->FiltLowTotBases = 1;
			continue;
			}
		pAlignSpecies->FiltLowTotBases = 0;
		NumSpeciesWithCnts += 1;

		// if proportion of major SNP base above a threshold then check if any of the other species have any bases
		BestAcceptBase = eBaseN;
		BestAcceptConf = 0.0;
		for(BaseIdx = 0; BaseIdx < 4; BaseIdx++)
			{
			CurConf = ((pAlignSpecies->ProbeBaseCnts[BaseIdx] + 1) / (double)pAlignSpecies->TotBases);
			if(pAlignSpecies->ProbeBaseCnts[BaseIdx] >= (uint32_t)MinCnt && (CurConf >= SNPMmajorPC))
				{
				bAcceptSpec = true;
				pAlignSpeciesA = pAlign;
				for(SpeciesIdxA = 0; SpeciesIdxA < NumSpecies; SpeciesIdxA++,pAlignSpeciesA += 1)
					{
					if(SpeciesIdxA == SpeciesIdx)
						continue;

					if(AltMaxCnt > 0 && pAlignSpeciesA->ProbeBaseCnts[BaseIdx] >= (uint32_t)AltMaxCnt)
						{
						bAcceptSpec = false;
						break;
						}
					}

				if(bAcceptSpec)
					{
					if(CurConf > BestAcceptConf)
						{
						BestAcceptBase = BaseIdx;
						BestAcceptConf = (pAlignSpecies->ProbeBaseCnts[BaseIdx] / (double)pAlignSpecies->TotBases);
						}
					}
				}

			pAlignSpecies->CultSpecBase = BestAcceptBase;
			if(pAlignSpecies->LSER < cFloorLSER)
				pAlignSpecies->LSER = m_LSER;
			}
		}
	pAlignSpecies = pAlign;
	for(SpeciesIdx = 0; SpeciesIdx < NumSpecies; SpeciesIdx++,pAlignSpecies += 1)
		pAlignSpecies->NumSpeciesWithCnts = NumSpeciesWithCnts; 
	}
gDiagnostics.DiagOut(eDLInfo, gszProcName, "Filtered: Processed %zd from %zd putative SNP loci", NumLociProc, PutSNPLoci);
return(0);
}

const int cRptBuffSize = 0x03fffff;          // will allocate this sized reporting buffer
int64_t 
CMarkers::Report(char *pszRefGenome,		// reference genome assembly against which other species were aligned
			int NumRelGenomes,				// number of relative genome names
			char *pszRelGenomes[],			// relative genome names
			char *pszReportFile,			// report to this file
			int MinSpeciesWithCnts,			// must be at least this number of species with base counts more than MinSpeciesTotCntThres - 0 if no limit 
			int MinSpeciesTotCntThres,		// individual species must have at least this number of total bases at SNP loci - 0 if no limit
			bool bSloughRefOnly,			// do not report if no inter-cultivar SNP marker, i.e if cultivars all same with the polymorphic site relative to reference only
			bool bSloughNonHetero)			// do not report unless all cultivars are relative heterozygotic - no two cultivars have same base 
{
uint32_t Idx;
char cBase;
tsAlignLoci *pAlign;
tsAlignLoci *pTmpAlign;
uint16_t ChkIdx;
char *pszTargSeq;
uint32_t NumCultivars;
uint32_t PrevTargSeqID;
int64_t NumSloughed;
int64_t NumReported;
char szBEDReportFile[_MAX_PATH+1];

SortTargSeqLociSpecies();

#ifdef _WIN32
m_hOutFile = open(pszReportFile,O_CREATETRUNC );
#else
if((m_hOutFile = open(pszReportFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
	if(ftruncate(m_hOutFile,0)!=0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",pszReportFile,strerror(errno));
			Reset();
			return(eBSFerrCreateFile);
			}
#endif

if(m_hOutFile < 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output file '%s'",pszReportFile);
	Reset();
	return(eBSFerrCreateFile);
	}

if((m_pszBuff = new char [cRptBuffSize]) == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Report: unable to allocate memory for buffering reports");
	Reset();
	return(eBSFerrMem);
	}
NumCultivars = m_NumSpecies - 1;

// only generating a BED if maximum of 4 cultivars and bSloughNonHetero
if(NumCultivars <= 4 && bSloughNonHetero)
	{
	int NameLen = (int)strlen(pszReportFile);
	strcpy(szBEDReportFile,pszReportFile);
	if(NameLen > 4 && !stricmp(&szBEDReportFile[NameLen-4],".csv"))
		szBEDReportFile[NameLen-4] = '\0';
	strcat(szBEDReportFile,".bed");
#ifdef _WIN32
	m_hBEDOutFile = open(szBEDReportFile,O_CREATETRUNC );
#else
	if((m_hBEDOutFile = open(szBEDReportFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
		if(ftruncate(m_hBEDOutFile,0)!=0)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",szBEDReportFile,strerror(errno));
				Reset();
				return(eBSFerrCreateFile);
				}
#endif

	if(m_hBEDOutFile < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output file '%s'",szBEDReportFile);
		Reset();
		return(eBSFerrCreateFile);
		}

	if((m_pszBEDBuff = new char [cRptBuffSize]) == nullptr)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Report: unable to allocate memory for buffering BED");
		Reset();
		return(eBSFerrMem);
		}

	m_BEDBuffIdx = sprintf(m_pszBEDBuff, "track name=\"SNP Dirac loci\" description=\"SNP Dirac Loci\" useScore=0\n");
	}


m_BuffIdx = sprintf(m_pszBuff,"\"%s:TargSeq\",\"Loci\",\"TargBase\",\"NumSpeciesWithCnts\"",pszRefGenome);
for(Idx = 0; Idx < NumCultivars; Idx++)
	{
	m_BuffIdx += sprintf(&m_pszBuff[m_BuffIdx],",\"%s:CntsSrc\",\"%s:Base\",\"%s:LSER\",\"%s:BaseCntTot\",\"%s:BaseCntA\",\"%s:BaseCntC\",\"%s:BaseCntG\",\"%s:BaseCntT\",\"%s:BaseCntN\"",
				pszRelGenomes[Idx],pszRelGenomes[Idx],pszRelGenomes[Idx],pszRelGenomes[Idx],pszRelGenomes[Idx],pszRelGenomes[Idx],pszRelGenomes[Idx],pszRelGenomes[Idx],pszRelGenomes[Idx]);
	if((m_BuffIdx + (cRptBuffSize/4)) > cRptBuffSize)
		{
		CUtility::RetryWrites(m_hOutFile, m_pszBuff, m_BuffIdx);
		m_BuffIdx = 0;
		}
	}
if(m_BuffIdx > 0)
	{
	CUtility::RetryWrites(m_hOutFile,m_pszBuff,m_BuffIdx);
	m_BuffIdx = 0;
	}

pAlign = &m_pAllocAlignLoci[0];
PrevTargSeqID = 0;
pszTargSeq = nullptr;
NumSloughed = 0;
NumReported = 0;
time_t Then = time(nullptr);
time_t Now;
int64_t LociIdx;
int64_t NumLociProc = 0;
int64_t PutSNPLoci = m_UsedAlignLoci / (int64_t)NumCultivars;
for(LociIdx = 0; LociIdx < m_UsedAlignLoci; LociIdx += NumCultivars)
	{
	if((NumLociProc % 10000)==0)
		{
		Now = time(nullptr);
		if ((Now - Then) >= 60)
			{
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "Reporting: Processed %zd putative SNP loci accepting %zd", NumLociProc, NumReported);
			Then += 60;
			}
		}
	NumLociProc++;
	pAlign = &m_pAllocAlignLoci[LociIdx];
	if(pAlign->NumSpeciesWithCnts < MinSpeciesWithCnts)
		{
		NumSloughed++;
		continue;
		}

	if(bSloughRefOnly) // user may have requested that only variants between the cultivars are of interest; if any cultivar has a indeterminate base then slough
		{
		pTmpAlign = pAlign;
		for(ChkIdx = 0; ChkIdx < NumCultivars; ChkIdx++,pTmpAlign++)
			{
			if(pTmpAlign->CultSpecBase == eBaseN)
				break;
			}
		if(ChkIdx != NumCultivars)	// at least 1 cultivar was indeterminate?
			{
			NumSloughed += 1;
			continue;
			}
		}

	// user may have requested that only variants between the cultivars are of interest; if all cultivars have same variant, even if different to reference, then slough
	if(bSloughNonHetero)
		{
		pTmpAlign = pAlign;
		uint8_t CultBase = 0x0ff;
		for(ChkIdx = 0; ChkIdx < NumCultivars; ChkIdx++,pTmpAlign++)
			{
			if(!pTmpAlign->FiltLowTotBases)
				{
				if(CultBase == 0x0ff)
					CultBase = pTmpAlign->CultSpecBase;
				else
					if(CultBase != pTmpAlign->CultSpecBase)
						break;
				}
			}
		if(ChkIdx == NumCultivars)	// all cultivar bases same or below filtering threshold?
			{
			NumSloughed += 1;
			continue;
			}
		}
	char RefBase;
	pTmpAlign = pAlign;
	switch(pTmpAlign->TargRefBase) {
		case eBaseA:
			RefBase = 'A';
			break;
		case eBaseC:
			RefBase = 'C';
			break;
		case eBaseG:
			RefBase = 'G';
			break;
		case eBaseT:
			RefBase = 'T';
			break;
		default:
			RefBase = 'N';
			break;
		}

	if(pszTargSeq == nullptr || pTmpAlign->TargSeqID != PrevTargSeqID)
		{
		if((pszTargSeq = SeqIDtoName(pTmpAlign->TargSeqID))==nullptr)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Report: unable to locate target sequence name for identifier: %u", pTmpAlign->TargSeqID);
			Reset();
			return(eBSFerrInternal);
			}
		PrevTargSeqID = pTmpAlign->TargSeqID;
		}
	m_BuffIdx += sprintf(&m_pszBuff[m_BuffIdx], "\n\"%s\",%d,\"%c\",%d", pszTargSeq, pTmpAlign->TargLoci, RefBase, pTmpAlign->NumSpeciesWithCnts);
		// if also reporting in BED format
	if(m_hBEDOutFile != -1)
		m_BEDBuffIdx += sprintf(&m_pszBEDBuff[m_BEDBuffIdx],"%s\t%u\t%u\t", pszTargSeq, pTmpAlign->TargLoci, pTmpAlign->TargLoci + 1);

	for(ChkIdx = 0; ChkIdx < (uint16_t)NumCultivars; ChkIdx++, pTmpAlign++)
		{
		switch(pTmpAlign->CultSpecBase) {
				case eBaseA:
					cBase = 'A';
					break;
				case eBaseC:
					cBase = 'C';
					break;
				case eBaseG:
					cBase = 'G';
					break;
				case eBaseT:
					cBase = 'T';
					break;
				case eBaseN:
					cBase = 'N';
					break;
				}
		m_BuffIdx += sprintf(&m_pszBuff[m_BuffIdx],",\"%c\",\"%c\",%f,%u,%u,%u,%u,%u,%u",(pTmpAlign->Flags & cFlgSNPcnts) ? 'S' : 'I',	cBase,
						pTmpAlign->LSER, pTmpAlign->TotBases, pTmpAlign->ProbeBaseCnts[0], pTmpAlign->ProbeBaseCnts[1], pTmpAlign->ProbeBaseCnts[2], pTmpAlign->ProbeBaseCnts[3], pTmpAlign->ProbeBaseCnts[4]);

		if(m_hBEDOutFile != -1)
			m_BEDBuffIdx += sprintf(&m_pszBEDBuff[m_BEDBuffIdx],"%c%c",cBase,ChkIdx==NumCultivars-1 ? '\n':',');
		}

	if((m_BuffIdx + (cRptBuffSize / 4)) > cRptBuffSize)
		{
		CUtility::RetryWrites(m_hOutFile,m_pszBuff,m_BuffIdx);
		m_BuffIdx = 0;
		}

	if(m_hBEDOutFile != -1 && ((m_BEDBuffIdx + (cRptBuffSize / 4)) > cRptBuffSize))
			{
			CUtility::RetryWrites(m_hBEDOutFile,m_pszBEDBuff,m_BEDBuffIdx);
			m_BEDBuffIdx = 0;
			}

	NumReported += 1;
	}

if(m_hOutFile != -1)
	{
	if(m_BuffIdx && m_pszBuff != nullptr)
		{
		CUtility::RetryWrites(m_hOutFile,m_pszBuff,m_BuffIdx);
		m_BuffIdx = 0;
		}
#ifdef _WIN32
	_commit(m_hOutFile);
#else
	fsync(m_hOutFile);
#endif
	close(m_hOutFile);
	m_hOutFile = -1;

	if(m_pszBuff != nullptr)
		{
		delete m_pszBuff;
		m_pszBuff = nullptr;
		}
	}

if(m_hBEDOutFile != -1)
	{
	if(m_BEDBuffIdx && m_pszBEDBuff != nullptr)
		{
		CUtility::RetryWrites(m_hBEDOutFile,m_pszBEDBuff,m_BEDBuffIdx);
		m_BEDBuffIdx = 0;
		}
#ifdef _WIN32
	_commit(m_hBEDOutFile);
#else
	fsync(m_hBEDOutFile);
#endif
	close(m_hBEDOutFile);
	m_hBEDOutFile = -1;
	if(m_pszBEDBuff != nullptr)
		{
		delete m_pszBEDBuff;
		m_pszBEDBuff = nullptr;
		}
	}

gDiagnostics.DiagOut(eDLInfo, gszProcName, "Reported: Processed %zd from putative SNP loci accepting %zd", NumLociProc, NumReported);
return(NumReported);
}

int
CMarkers::CreateMutexes(void)
{
	if (m_bMutexesCreated)
		return(eBSFSuccess);

#ifdef _WIN32
	InitializeSRWLock(&m_hRwLock);
#else
	if (pthread_rwlock_init(&m_hRwLock, nullptr) != 0)
	{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Fatal: unable to create rwlock");
		return(eBSFerrInternal);
	}
#endif

#ifdef _WIN32
	if ((m_hMtxIterReads = CreateMutex(nullptr, false, nullptr)) == nullptr)
	{
#else
	if (pthread_mutex_init(&m_hMtxIterReads, nullptr) != 0)
	{
		pthread_rwlock_destroy(&m_hRwLock);
#endif
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Fatal: unable to create mutex");
		return(eBSFerrInternal);
	}

#ifdef _WIN32
	if ((m_hMtxMHReads = CreateMutex(nullptr, false, nullptr)) == nullptr)
	{
#else
	if (pthread_mutex_init(&m_hMtxMHReads, nullptr) != 0)
	{
#endif
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Fatal: unable to create mutex");
#ifdef _WIN32
		CloseHandle(m_hMtxIterReads);
#else
		pthread_rwlock_destroy(&m_hRwLock);
		pthread_mutex_destroy(&m_hMtxIterReads);
#endif
		return(eBSFerrInternal);
	}

	m_bMutexesCreated = true;
	return(eBSFSuccess);
	}

void
CMarkers::DeleteMutexes(void)
{
	if (!m_bMutexesCreated)
		return;
#ifdef _WIN32
	CloseHandle(m_hMtxIterReads);
	CloseHandle(m_hMtxMHReads);
#else
	pthread_mutex_destroy(&m_hMtxIterReads);
	pthread_mutex_destroy(&m_hMtxMHReads);
	pthread_rwlock_destroy(&m_hRwLock);
#endif
	m_bMutexesCreated = false;
}

void
CMarkers::AcquireSerialise(void)
{
#ifdef _WIN32
	WaitForSingleObject(m_hMtxIterReads, INFINITE);
#else
	pthread_mutex_lock(&m_hMtxIterReads);
#endif
}

void
CMarkers::ReleaseSerialise(void)
{
#ifdef _WIN32
	ReleaseMutex(m_hMtxIterReads);
#else
	pthread_mutex_unlock(&m_hMtxIterReads);
#endif
}

void
CMarkers::AcquireSerialiseMH(void)
{
#ifdef _WIN32
	WaitForSingleObject(m_hMtxMHReads, INFINITE);
#else
	pthread_mutex_lock(&m_hMtxMHReads);
#endif
}

void
CMarkers::ReleaseSerialiseMH(void)
{
#ifdef _WIN32
	ReleaseMutex(m_hMtxMHReads);
#else
	pthread_mutex_unlock(&m_hMtxMHReads);
#endif
}

void
CMarkers::AcquireLock(bool bExclusive)
{
#ifdef _WIN32
	if (bExclusive)
		AcquireSRWLockExclusive(&m_hRwLock);
	else
		AcquireSRWLockShared(&m_hRwLock);
#else
	if (bExclusive)
		pthread_rwlock_wrlock(&m_hRwLock);
	else
		pthread_rwlock_rdlock(&m_hRwLock);
#endif
}

void
CMarkers::ReleaseLock(bool bExclusive)
{

#ifdef _WIN32
	if (bExclusive)
		ReleaseSRWLockExclusive(&m_hRwLock);
	else
		ReleaseSRWLockShared(&m_hRwLock);
#else
	pthread_rwlock_unlock(&m_hRwLock);
#endif
}


static tsAlignLoci *gpAllocAlignLoci;		// as allocated to hold alignment loci;

int64_t
CMarkers::SortTargSeqLociSpecies(void)	
{
// check if anything to sort ....
if(m_pAllocAlignLoci == nullptr || m_UsedAlignLoci == 0)
	{
	m_bSorted = false;
	return(0);
	}

if (m_UsedAlignLoci == 1)
	{
	m_bSorted = true;
	return(1);
	}
	
if(!m_bSorted)
	{
	// trying with 16 threads to see what throughput improvement we gain over a single threaded approach
	gpAllocAlignLoci = m_pAllocAlignLoci;
	m_MTqsort.qsort(gpAllocAlignLoci, m_UsedAlignLoci, sizeof(tsAlignLoci), QSortAlignSeqLociSpecies);
	m_bSorted = true;
	}
return(m_UsedAlignLoci);
}

// QSortAlignSeqLociSpecies
// qsorts alignment loci by TargSeqID,TargLoci,ProbeSpeciesID ascending
int CMarkers::QSortAlignSeqLociSpecies(const void *p1,const void *p2)
{
tsAlignLoci *pAlign1;
tsAlignLoci *pAlign2;

pAlign1 = (tsAlignLoci *)p1;
pAlign2 = (tsAlignLoci *)p2;
if(pAlign1->TargSeqID > pAlign2->TargSeqID)
	return(1);
if(pAlign1->TargSeqID < pAlign2->TargSeqID)
	return(-1);
if(pAlign1->TargLoci > pAlign2->TargLoci)
	return(1);
if(pAlign1->TargLoci < pAlign2->TargLoci)
	return(-1);
if(pAlign1->ProbeSpeciesID > pAlign2->ProbeSpeciesID)
	return(1);
if(pAlign1->ProbeSpeciesID < pAlign2->ProbeSpeciesID)
	return(-1);
return(0);
}