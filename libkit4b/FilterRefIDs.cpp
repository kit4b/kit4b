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
#ifdef _WIN32
#include "./commhdrs.h"
#else
#include "./commhdrs.h"
#endif



CFilterRefIDs::CFilterRefIDs(void)
{
m_pFilterRefIDs = NULL;
m_NumFilterRefIDs = 0;
m_AllocdFilterRefIDs = 0;
}

CFilterRefIDs::~CFilterRefIDs(void)
{
if(m_pFilterRefIDs != NULL)
	delete m_pFilterRefIDs;
}

void
CFilterRefIDs::Reset(void)
{
if(m_pFilterRefIDs != NULL)
	{
	delete m_pFilterRefIDs;
	m_pFilterRefIDs= NULL;
	}
m_NumFilterRefIDs = 0;
m_AllocdFilterRefIDs = 0;
}


int
CFilterRefIDs::Open(char *pszFile)
{
int Rslt;
int *pRefIDs;
int RefID;
int NumFields;

if((Rslt=CCSVFile::Open(pszFile))!=eBSFSuccess)
	{
	AddErrMsg("CFilterRefIDs::Open","Unable to open '%s'",pszFile);
	return(Rslt);
	}

// now parse each RefID which is to be filtered out
while((Rslt=NextLine()) > 0)	// onto next line containing fields
	{
	NumFields = GetCurFields();
	if(NumFields < 1)
		{
		AddErrMsg("CFilterRefIDs::Open","Expected at least 1 fields in '%s', GetCurFields() returned '%d'",pszFile,NumFields);
		return(eBSFerrFieldCnt);
		}
	if(!m_NumFilterRefIDs && IsLikelyHeaderLine())
		continue;
	GetInt(1,&RefID);

	if(m_pFilterRefIDs == NULL || m_NumFilterRefIDs == m_AllocdFilterRefIDs)
		{
		if((pRefIDs = new int [m_AllocdFilterRefIDs + cFiltRefIDsAllocChunk])==NULL)
			{
			AddErrMsg("CFilterRefIDs::Open","Unable to alloc memory for filtered RefIDs");
			return(eBSFerrMem);
			}

		if(m_pFilterRefIDs != NULL)
			{
			if(m_NumFilterRefIDs)
				memmove(pRefIDs,m_pFilterRefIDs,sizeof(int) * m_NumFilterRefIDs);
			delete m_pFilterRefIDs;
			}
		else
			{
			m_AllocdFilterRefIDs = 0;
			m_NumFilterRefIDs = 0;
			}
		m_AllocdFilterRefIDs += cFiltRefIDsAllocChunk;
		m_pFilterRefIDs = pRefIDs;
		}
	m_pFilterRefIDs[m_NumFilterRefIDs++] = RefID;
	}

CCSVFile::Close();
if(m_NumFilterRefIDs > 1)
	qsort(m_pFilterRefIDs,m_NumFilterRefIDs,sizeof(int),SortRefIDs);

return(m_NumFilterRefIDs);
}

// returns true if RefID is in RefIDs loaded
bool 
CFilterRefIDs::Locate(int RefID)
{
int *pProbe;
int Lo,Mid,Hi;	// search limits
Lo = 0; Hi = m_NumFilterRefIDs-1;
while(Hi >= Lo) {
	Mid = (Hi + Lo)/2;
	pProbe = &m_pFilterRefIDs[Mid];
	if(RefID == *pProbe)
		return(true);
	if(RefID < *pProbe)	
		{
		Hi = Mid - 1;
		continue;
		}
	Lo = Mid + 1;
	}
return(false);
}


int 
CFilterRefIDs::SortRefIDs( const void *arg1, const void *arg2)
{
return(*(int *)arg1 - *(int *)arg2);
}



