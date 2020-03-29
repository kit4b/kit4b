/*
This toolkit is a source base clone of 'BioKanga' release 4.4.2 (https://github.com/csiro-crop-informatics/biokanga) and contains
significant source code changes enabling new functionality and resulting process parameterisation changes. These changes have resulted in
incompatibilty with 'BioKanga'.

Because of the potentential for confusion by users unaware of functionality and process parameterisation changes then the modified source base
and resultant compiled executables have been renamed to 'kit4b' - K-mer Informed Toolkit for Bioinformatics.
The renaming will force users of the 'BioKanga' toolkit to examine scripting which is dependent on existing 'BioKanga'
parameterisations so as to make appropriate changes if wishing to utilise 'kit4b' parameterisations and functionality.

'kit4b' is being released under the Opensource Software License Agreement (GPLv3)
'kit4b' is Copyright (c) 2019, 2020
Please contact Dr Stuart Stephen < stuartjs@g3web.com > if you have any questions regarding 'kit4b'.

Orginal 'BioKanga' copyright notice has been retained and immediately follows this notice..
*/
/*
 * CSIRO Open Source Software License Agreement (GPLv3)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENSE for the complete license information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */

#include "stdafx.h"
#ifdef _WIN32
#include "./commhdrs.h"
#else
#include "./commhdrs.h"
#endif

const int cNumAllocChroms = 1000;		// allocate for chromosomes in this many instances
const int cNumAllocLoci	  = 50000;		// allocate for loci in this many instances


CFeatLoci::CFeatLoci(void)
{
m_NumLoci = 0;
m_AllocLoci = 0;
m_NumChroms = 0;
m_AllocChroms = 0;
m_pLocii = NULL;
m_pChroms = NULL;
}

CFeatLoci::~CFeatLoci(void)
{
if(m_pLocii != NULL)
	delete m_pLocii;
if(m_pChroms != NULL)
	delete m_pChroms;

}

// GenNameHash
// Generates a 16bit hash on specified name
// This hash can then be used to quickly eliminate probe names which can't match a target name by comparing hashes
unsigned short
CFeatLoci::GenNameHash(char *pszName)
{
unsigned short Hash = 0;
unsigned short MSB;
char chr;
while(chr = *pszName++)
	{
	MSB = (Hash & 0x0c000) >> 13;
	Hash <<= 2;
	Hash ^= tolower(chr);
	Hash ^= MSB;
	}
return(Hash);
}


int 
CFeatLoci::GetNumLoci(void)
{
return(m_NumLoci);
}

int
CFeatLoci::LocateChrom(char *pszChrom)
{
int Idx;
tsChrom *pChrom;
unsigned short Hash;
if(m_NumChroms > 0 && (pChrom = m_pChroms)!=NULL)
	{
	Hash = GenNameHash(pszChrom);
	for(Idx = 0; Idx < m_NumChroms; Idx++,pChrom++)
		{
		if(pChrom->Hash == Hash && !stricmp(pszChrom,pChrom->szName))
			return(pChrom->ChromID);
		}
	}
return(0);	// couldn't locate
}

int
CFeatLoci::AddChrom(char *pszChrom)
{
int ChromID;
tsChrom *pChrom;

if((ChromID = LocateChrom(pszChrom)) > 0)
	return(ChromID);

if(m_pChroms == NULL || m_NumChroms >= m_AllocChroms)
	{
	if((pChrom = new tsChrom[m_AllocChroms + cNumAllocChroms])==NULL)
		return(-1);

	if(m_pChroms != NULL)
		{
		if(m_NumChroms)
			memmove(pChrom,m_pChroms,sizeof(tsLoci) * m_NumChroms);
		delete m_pChroms;
		}
	m_AllocChroms += cNumAllocChroms;
	m_pChroms = pChrom;
	}
pChrom = &m_pChroms[m_NumChroms++];
pChrom->ChromID = m_NumChroms;
strcpy(pChrom->szName,pszChrom);
pChrom->Hash = GenNameHash(pszChrom);
return(m_NumChroms);
}

int 
CFeatLoci::AddLoci(int RefID,char *pszChrom, int StartLoci, int EndLoci)
{
int Rslt;
tsLoci *pLoci;
int ChromID;
if((Rslt = ChromID = AddChrom(pszChrom)) < 1)
	return(Rslt);

if(m_pLocii == NULL || m_NumLoci >= m_AllocLoci)
	{
	if((pLoci = new tsLoci[m_AllocLoci + cNumAllocLoci])==NULL)
		return(-1);

	if(m_pLocii != NULL)
		{
		if(m_NumLoci)
			memmove(pLoci,m_pLocii,sizeof(tsLoci) * m_NumLoci);
		delete m_pLocii;
		}
	m_AllocLoci += cNumAllocLoci;
	m_pLocii = pLoci;
	}
pLoci = &m_pLocii[m_NumLoci++];
pLoci->RefID = RefID;
pLoci->ChromID = ChromID;
pLoci->StartLoci = StartLoci;
pLoci->EndLoci = EndLoci;
return(m_NumLoci);
}

// ClusterLoci
// Any loci which is within MaxDistance of another loci is assumed to be part of a larger clustered loci 
int
CFeatLoci::ClusterLoci(int MaxDistance)
{
int NumEls;
int Idx;
int ChromID;
tsLoci *pLoci;
tsLoci *pNxt;
if(m_pLocii == NULL || m_NumLoci < 2)
	return(m_NumLoci);

// firstly sort loci
qsort(m_pLocii,m_NumLoci,sizeof(tsLoci),SortLoci);
if(MaxDistance > 0)	// need to cluster?
	{
	// iterate over all loci, any loci within MaxDistance is combined into a single loci
	pLoci = m_pLocii;
	NumEls = 1;
	for(pNxt = &m_pLocii[1],Idx = 1; Idx < m_NumLoci; Idx++,pNxt++)
		{
		if(pNxt->ChromID != pLoci->ChromID)	// different chromosome?
			{
			pLoci += 1;
			if(pLoci != pNxt)
				*pLoci = *pNxt;
			NumEls += 1;
			continue;
			}

		// same chromosome
		// check if overlaid loci
		if(pLoci->EndLoci >= pNxt->StartLoci)
			{
			if(pNxt->EndLoci > pLoci->EndLoci)
				pLoci->EndLoci = pNxt->EndLoci;
			continue;
			}

		// check if separated by more than MaxDistance
		if((pLoci->EndLoci + MaxDistance) < pNxt->StartLoci)
			{
			pLoci += 1;
			if(pLoci != pNxt)
				*pLoci = *pNxt;
			NumEls += 1;
			continue;
			}

		// separated by less than MaxDistance so combine
		pLoci->EndLoci = pNxt->EndLoci;
		}
	m_NumLoci = NumEls;
	}

// now process for the starting loci on each chromosome
ChromID = -1;
for(pLoci = m_pLocii,Idx = 0; Idx < m_NumLoci; Idx++,pLoci++)
	{
	if(pLoci->ChromID != ChromID)	// different chromosome?
		{
		m_pChroms[pLoci->ChromID-1].LociID = Idx + 1;
		ChromID = pLoci->ChromID;
		}
	}

return(m_NumLoci);
}

// Locate1stLoci
// Locate 1st loci on specified chromosome
int
CFeatLoci::Locate1stLoci(int ChromID)
{
if(m_pLocii == NULL || ChromID < 1 || ChromID > m_NumChroms)
	return(-1);
return(m_pChroms[ChromID-1].LociID);
}

int
CFeatLoci::Locate1stLoci(char *pszChrom)
{
int ChromID;
if((ChromID = LocateChrom(pszChrom)) < 1)
	return(ChromID);
return(Locate1stLoci(ChromID));
}

// returns next loci on same chromosome as that for LociID (1..n)
int 
CFeatLoci::LocateNxtLociChrom(int LociID)
{
tsLoci *pLoci;
int ChromID;
if(m_pLocii == NULL || LociID < 1 || LociID >= m_NumLoci)
	return(-1);
pLoci = &m_pLocii[LociID-1];
ChromID = pLoci->ChromID;
if(ChromID == pLoci[1].ChromID)
	return(LociID + 1);
return(0);
}


int 
CFeatLoci::GetLoci(int LociID,int *pRefID,char **ppszChrom,int *pStartLoci,int *pEndLoci)
{
tsLoci *pLoci;
if(m_pLocii == NULL || LociID < 1 || LociID > m_NumLoci)
	return(-1);
pLoci = &m_pLocii[LociID-1];
if(pRefID != NULL)
	*pRefID = pLoci->RefID;
if(ppszChrom != NULL)
	*ppszChrom = m_pChroms[pLoci->ChromID-1].szName;
if(pStartLoci != NULL)
	*pStartLoci = pLoci->StartLoci;
if(pEndLoci != NULL)
	*pEndLoci = pLoci->EndLoci;
return(LociID);
}

// SortLoci
// Used to sort by ChromID ---> Start ---> End
int 
CFeatLoci::SortLoci( const void *arg1, const void *arg2)
{
tsLoci *pEl1 = (tsLoci *)arg1;
tsLoci *pEl2 = (tsLoci *)arg2;

if(pEl1->ChromID < pEl2->ChromID)
	return(-1);
if(pEl1->ChromID > pEl2->ChromID)
	return(1);
if(pEl1->StartLoci < pEl2->StartLoci)
	return(-1);
if(pEl1->StartLoci > pEl2->StartLoci)
	return(1);
if(pEl1->EndLoci < pEl2->EndLoci)
	return(-1);
if(pEl1->EndLoci > pEl2->EndLoci)
	return(1);
return(0);
}
