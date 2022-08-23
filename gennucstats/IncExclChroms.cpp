#include "stdafx.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#if _WIN32
#include "../libkit4b/commhdrs.h"
#else
#include "../libkit4b/commhdrs.h"
#endif

#include "IncExclChroms.h"


CIncExclChroms::CIncExclChroms()
{
m_AllocdIncExclChroms = 0;
m_IncExclChroms = 0;
m_pIncExclChroms = NULL;
m_pLastmatch = NULL;
}

CIncExclChroms::~CIncExclChroms()
{
if(m_pIncExclChroms != NULL)
	delete m_pIncExclChroms;
}

void CIncExclChroms::Reset(void)
{
if(m_pIncExclChroms != NULL)
	{
	delete m_pIncExclChroms;
	m_pIncExclChroms = NULL;
	}
m_AllocdIncExclChroms = 0;
m_IncExclChroms = 0;
m_pLastmatch = NULL;
}

int
CIncExclChroms::InitChromExclusion(int	NumIncludeChroms,	// number of chromosome regular expressions to include
		char **ppszIncludeChroms,	// array of include chromosome regular expressions
		int	NumExcludeChroms,		// number of chromosome expressions to exclude
		char **ppszExcludeChroms)	// array of exclude chromosome regular expressions
{
int Rslt;
Reset();
Rslt = m_RegExprs.CompileREs(NumIncludeChroms,ppszIncludeChroms,NumExcludeChroms,ppszExcludeChroms);
if(Rslt < eBSFSuccess)
	Reset();
return(Rslt);
}

//IncludeThisChrom
//Returns 0 if to be excluded, > 1 as the chromid when to be included or < 0 if error
//
int
CIncExclChroms::IncludeThisChrom(char *pszChrom)
{
int ChromIdx;
bool bProcChrom;
tsIncExclChrom *pIncExcl;

// if no include/exclude then chromosome is to be included
if(!m_RegExprs.HasRegExprs())
	{
	if((pIncExcl = LocateChrom(pszChrom))!=NULL)
		return(pIncExcl->ChromID);
	return(AddChrom(pszChrom,true));
	}

// optimisation is to check if chromosome previously known to be included/excluded
if((pIncExcl = LocateChrom(pszChrom))!=NULL)
	return((int)pIncExcl->bInc ? pIncExcl->ChromID : 0);

bProcChrom = m_RegExprs.Accept(pszChrom);
ChromIdx = AddChrom(pszChrom,bProcChrom);
return(bProcChrom ? ChromIdx : 0);
}

tsIncExclChrom *
CIncExclChroms::LocateChrom(char *pszChrom)
{
tsIncExclChrom *pTmp;
int Idx;
if(m_pIncExclChroms == NULL || m_IncExclChroms == 0)
	return(NULL);
if(m_pLastmatch != NULL && !stricmp(m_pLastmatch->szChrom,pszChrom))
	return(m_pLastmatch);

pTmp = m_pIncExclChroms;
for(Idx = 0; Idx < m_IncExclChroms; Idx++, pTmp++)
	if(!stricmp(pTmp->szChrom,pszChrom))
		{
		m_pLastmatch = pTmp;
		return(pTmp);
		}
return(NULL);
}

char *
CIncExclChroms::LocateChrom(int ChromID)	// returns chrom corresponding to specified ChromID
{
if(m_pIncExclChroms == NULL || ChromID < 1 || ChromID > m_IncExclChroms)
	return(NULL);
return(m_pIncExclChroms[ChromID-1].szChrom);
}


int					// returns total number of chromosomes added, becomes the chrom identifier! ( < 0 if error)
CIncExclChroms::AddChrom(char *pszChrom,bool bProcChrom) // caches chrom processing state
{
tsIncExclChrom *pTmp;
if(m_pIncExclChroms == NULL || m_IncExclChroms >= m_AllocdIncExclChroms)
	{
	if(m_pIncExclChroms == NULL)
		{
		m_IncExclChroms = 0;
		m_AllocdIncExclChroms = 0;
		}
	if((pTmp = new tsIncExclChrom [m_AllocdIncExclChroms + cAllocIncExclChroms])==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory (%d bytes) for holding IncExcl chromosomes",sizeof(tsIncExclChrom) * (m_AllocdIncExclChroms + cAllocIncExclChroms));
		return(eBSFerrMem);
		}
	if(m_pIncExclChroms != NULL)
		{
		memcpy(pTmp,m_pIncExclChroms,sizeof(tsIncExclChrom) * m_IncExclChroms);
		delete m_pIncExclChroms;
		}
	m_AllocdIncExclChroms += cAllocIncExclChroms;
	m_pIncExclChroms = pTmp;
	}
pTmp = &m_pIncExclChroms[m_IncExclChroms++];
pTmp->ChromID = m_IncExclChroms;
pTmp->bInc = bProcChrom;
strncpy(pTmp->szChrom,pszChrom,cMaxDatasetSpeciesChrom);
pTmp->szChrom[cMaxDatasetSpeciesChrom-1] = '\0';
m_pLastmatch = pTmp;
return(m_IncExclChroms);
}
