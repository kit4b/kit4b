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

CSimReads::CSimReads()
{
m_pBioSeqFile = NULL;
m_pBEDFile = NULL;
m_pChromSeqs = NULL;
m_pProfSel = NULL;			// allocated array of profile site selection preferences (0.0..1.0) indexed by sequence octamers
m_pProfCSV = NULL;			// used to hold DNase site preferences whilst loading into m_pProfSel
m_hOutFile = -1;
m_hOutPEFile = -1;
m_pSimReads = NULL;			// allocated to hold simulated reads
m_pGenomeSeq = NULL;		// allocated to hold concatenated (separated by eBaseEOSs) chromosome sequences

Reset(false);
}


CSimReads::~CSimReads()
{
Reset(false);
}


// TrimSeqChrs
// Removes any quotes and whitespace (space and tabs only) from pszTxt
// Also checks for legal base chars 'acgt'
// Returns -1 if illegal char, 0 if empty sequence, otherwise the sequence length
int
CSimReads::TrimSeqChrs(char *pszTxt)
{
char *pDst = pszTxt;
char Chr;
int Len = 0;
while((Chr = *pszTxt++))
	{
	if(Chr == '\0' || (Chr == '\t'  || Chr == ' '  || Chr == '"' || Chr == '\''))
		continue;
	switch(Chr) {
		case 'a': case 'A':
			Chr = 'a';
			break;
		case 'c': case 'C':
			Chr = 'c';
			break;
		case 'g': case 'G':
			Chr = 'g';
			break;
		case 't': case 'T': case 'u': case 'U':
			Chr = 't';
			break;
		default:
			return(-1);
		}
	*pDst++ = Chr;
	Len++;
	}
*pDst = Chr;
return(Len);
}

char *
CSimReads::Region2Txt(etSRBEDRegion Region)
{
switch(Region) {
	case eSRMEGRAny:		// process any region
		return((char *)"All");

	case eSRMEGRCDS:	// only process CDS
		return((char *)"CDS");

	case eSRMEGR5UTR:	// only process 5' UTR
		return((char *)"5' UTR");

	case eSRMEGR3UTR:	// only process 3' UTR
		return((char *)"3' UTR");

	case eSRMEGRIntrons:		// only process Introns
		return((char *)"Introns");

	case eSRMEGRUpstream:		// only process 5'upstream regulatory
		return((char *)"5'US");

	case eSRMEGRDnstream:		// only process 3'upstream regulatory
		return((char *)"3'DS");

	case eSRMEGRIntergenic:		// only process intergenic
		return((char *)"Intergenic");

	default:
		break;
	}
return((char *)"Unsupported");
}

TRandomCombined<CRandomMother,CRandomMersenne> RGseeds((int)time(0));


int
CSimReads::SimInDels(tsSRSimRead *pSimRead,int *pReadLen,etSeqBase *pRead)
{
int Idx;
etSeqBase InsertBases[20];
bool bInsert;
double Thres;
int InDelSize;
int InDelPsn;
if(m_InDelRate == 0.0 || m_InDelSize == 0)
	return(0);
Thres = RGseeds.Random();	// pick a number, any number, between 0 and 1 inclusive
if(Thres > m_InDelRate)		// InDel this read?
	return(0);

bInsert = RGseeds.IRandom(0,1) == 0 ? false : true;			// pick insert or delete
if(bInsert)
	{
	for(Idx = 0; Idx < m_InDelSize; Idx++)
		InsertBases[Idx] = (etSeqBase) RGseeds.IRandom(0,3);
	InDelSize = (int)RGseeds.IRandom(1,m_InDelSize);            // pick the insertion size
	InDelPsn = (int)RGseeds.IRandom(0,*pReadLen-InDelSize);     // pick the insert position
	memcpy(&pRead[InDelPsn+InDelSize],&pRead[InDelPsn],*pReadLen-(InDelPsn+InDelSize));
	memcpy(&pRead[InDelPsn],InsertBases,InDelSize);
	pSimRead->InDelLen = InDelSize;
	return(-1 * InDelSize);
	}
else
	{
	InDelSize = (int)RGseeds.IRandom(1,m_InDelSize);            // pick the deletion size
	InDelPsn = (int)RGseeds.IRandom(0,*pReadLen-InDelSize);     // pick the deletion position
	memcpy(&pRead[InDelPsn],&pRead[InDelPsn+InDelSize],*pReadLen-(InDelPsn+InDelSize));
	pSimRead->InDelLen = -1 * InDelSize;
	return(InDelSize);
	}
}

int
CSimReads::SimArtefacts(bool b3ArtefSeq,		// if false then 5' artefact else if true then 3' artefact
			int ReadLen,etSeqBase *pRead)
{
int NumArtSeqs;
double ArtefRate;

int ArtefSeqIdx;
int ArtefSeqLen;
int ArtefactLen;
etSeqBase *pArtefSeq;

double Thres;

if(b3ArtefSeq)
	{
	ArtefRate = m_Artef3Rate;
	NumArtSeqs = m_NumArtef3Seqs;
	}
else
	{
	ArtefRate = m_Artef5Rate;
	NumArtSeqs = m_NumArtef5Seqs;
	}

if(NumArtSeqs == 0 || ArtefRate == 0.0)
	return(ReadLen);

Thres = RGseeds.Random();	// pick a number, any number, between 0 and 1 inclusive
if(ArtefRate < 1.0 && Thres > ArtefRate)	// 5' artefact this read?
	return(0);

if(NumArtSeqs > 1)
	ArtefSeqIdx = RGseeds.IRandom(0,NumArtSeqs);			// pick artefact sequence
else
	ArtefSeqIdx = 0;
if(b3ArtefSeq)
	{
	pArtefSeq = m_Artef3Seqs[ArtefSeqIdx];
	ArtefSeqLen = m_Artef3SeqLens[ArtefSeqIdx];
	}
else
	{
	pArtefSeq = m_Artef5Seqs[ArtefSeqIdx];
	ArtefSeqLen = m_Artef5SeqLens[ArtefSeqIdx];
	}

ArtefactLen = min(ArtefSeqLen,ReadLen - 10);
ArtefactLen = RGseeds.IRandom(1,ArtefSeqLen);

if(b3ArtefSeq)
	memcpy(&pRead[ReadLen - ArtefactLen],pArtefSeq,ArtefactLen);
else
	{
	memcpy(&pRead[ArtefactLen],pRead,ReadLen - ArtefactLen);
	memcpy(pRead,&pArtefSeq[ArtefSeqLen - ArtefactLen],ArtefactLen);
	}

return(ReadLen);
}


// static profile - sequencing error rates as initially hardcoded
// updated (15th August 2013) to better represent the observed error rates with current Illumina 100bp read sets
// assumes mean substitution rate per 100bp of 1%
// =POISSON.DIST(NumSubs,MeanExpectedSubs=1,FALSE)
tsSRInduceErrProf StaticErrProfile[] =
	{
	{ 0.367879, 0},	// proportion of reads with no substitutions
	{ 0.367879, 1},	// proportion of reads  with 1 substitution
	{ 0.183944, 2},	// proportion of reads   with 2 substitutions
	{ 0.061313, 3},	// proportion of reads  with 3 substitutions
	{ 0.015328, 4},	// proportion of reads  with 4 substitutions
	{ 0.003066, 5},	// proportion of reads   with 5 substitutions
	{ 0.000511, 6},	// proportion of reads  with 6 substitutions
	{ 0.000073, 7},	// proportion of reads  with 7 substitutions
	{ -1.0f, 8 }	// remainder of reads  with 8 substitutions
	};
#define cNumProfEls (int)(sizeof(StaticErrProfile)/sizeof(tsSRInduceErrProf))

// dynamic profile - initialised according to user specified induced error rates
tsSRInduceErrProf DynErrProfile[cNumProfEls];


// Illumina cumulative error profile along length of read with moderate increase in subs at 5' start of reads but most subs are down at the 3' end of reads
int IlluminaSpatialDist[20] = { 40,55,64,72,80,88,96,104,112,121,131,142,156,174,197,228,270,325,400,500 };
#define cNumIlluminaSpatialDist (int)(sizeof(IlluminaSpatialDist)/sizeof(int))  

etSRSEMode m_SEMode;				// simulated sequencer error mode
double m_DynProbErr = 0.01;		// default as being a 1% error rate
bool m_bUniformDist = true;		// default as being a uniform distribution of induced errors
int m_PrvProfSeqLen = -1;		// set to read length for which DynErrProfile has been generated
int m_InducedErrDist[cNumProfEls];	// to hold induced error count distribution
int m_InducedErrPsnDist[2000];  // to hold read sequence psn induced error count


// CAUTION: uses the static profile for
// determining composite read distribution
int // number of substitutions inplace induced into this read
CSimReads::SimSeqErrors(int SeqLen,etSeqBase *pRead)
{
int Idx;
int RandSubs;
int SubBase;
int NumSubs2Induce;
etSeqBase *pSeq;
double Thres;

UINT8 *pSubd;

if(m_SEMode == eSRSEPnone)
	return(0);

pSeq = pRead;
for(Idx = 0; Idx < SeqLen; Idx++,pSeq++)
	*pSeq = *pSeq & ~cRptMskFlg;

// determine how many errors are to be induced into the current read
tsSRInduceErrProf *pProfile;
switch(m_SEMode) {
	case eSRSEPdyn:
		if(m_PrvProfSeqLen != SeqLen)
			{
			pProfile = &DynErrProfile[0];
			double CurThres = pow(1.0-m_DynProbErr,(double)SeqLen);
			double AccThres = 0.0;
			for(Idx = 0; Idx < (cNumProfEls-1); Idx++,pProfile++)
				{
				pProfile->NumSubs = Idx;
				pProfile->Proportion = CurThres;
				AccThres += CurThres;
				CurThres = (1-AccThres)/2;
				}
			pProfile->NumSubs = cNumProfEls-1;
			pProfile->Proportion = -1.0;		// flags last
			m_PrvProfSeqLen = SeqLen;
			}
		pProfile = &DynErrProfile[0];
		break;
	case eSRSEPstatic:
		pProfile = StaticErrProfile;
		break;
	case eSRSEPfixerrs:
		if((int)m_DynProbErr <= 0)
			{
			m_InducedErrDist[0] += 1;
			return(0);
			}
		pProfile = NULL;
		break;
	}

if(pProfile != NULL)
	{
	Thres = RGseeds.Random();	// pick a number, any number, between 0 and 1 inclusive
	for(Idx = 0; Idx < cNumProfEls; Idx++, pProfile++)
		{
		if(pProfile->Proportion <= 0.0 || pProfile->Proportion >= Thres)
			break;
		Thres -= pProfile->Proportion;
		}
	if(!pProfile->NumSubs)
		{
		m_InducedErrDist[0] += 1;
		return(0);
		}
	}

if(pProfile != NULL)
	NumSubs2Induce = pProfile->NumSubs;
else
	NumSubs2Induce = (int)m_DynProbErr;

pSubd = new UINT8 [SeqLen+10];
memset(pSubd,0,SeqLen);
RandSubs = 0;

int Psn;
while(RandSubs < NumSubs2Induce)
	{
	if(m_bUniformDist)
		Idx = RGseeds.IRandom(0,SeqLen);
	else
		{
		int MinPsn;
		int MaxPsn;
		int DistIdx;
		Psn = RGseeds.IRandom(0,IlluminaSpatialDist[cNumIlluminaSpatialDist-1]);
if(Psn < 0 || Psn > 500)
	printf("\n Check psn");
		for(DistIdx = 0; DistIdx < cNumIlluminaSpatialDist-1; DistIdx++)
			{
			if(Psn <= IlluminaSpatialDist[DistIdx])
				break;
			}
		MinPsn = (DistIdx * SeqLen) / cNumIlluminaSpatialDist;
		if(DistIdx == cNumIlluminaSpatialDist - 1)
			MaxPsn = SeqLen-1;
		else
			MaxPsn = (MinPsn + (SeqLen / cNumIlluminaSpatialDist)) - 1;

		Psn = RGseeds.IRandom(MinPsn,MaxPsn);
if(Psn < MinPsn || Psn > MaxPsn)
	printf("\n Check psn A");
		Idx=min(SeqLen-1,Psn);
		}
	if(pSubd[Idx] != 0)
		continue;

	pSeq = &pRead[Idx];
	do {
		SubBase = RGseeds.IRandom(0,3);
		}
	while(SubBase == *pSeq);
	*pSeq = SubBase | cRptMskFlg;
	pSubd[Idx] = 1;
	m_InducedErrPsnDist[Idx] += 1;
	RandSubs += 1;
	}
m_InducedErrDist[NumSubs2Induce] += 1;
delete pSubd;
return(NumSubs2Induce);
}

//
// Randomises all bases in this sequence making it highly unlikely that it can be aligned
int // number of substitutions inplace induced into this read
CSimReads::SimSeqRand(int SeqLen,etSeqBase *pRead)
{
int Idx;
int SubBase;
etSeqBase *pSeq;
if(m_SEMode == eSRSEPnone)
	return(0);
pSeq = pRead;
for(Idx = 0; Idx < SeqLen; Idx++,pSeq++)
	{
	*pSeq = *pSeq & ~cRptMskFlg;
	while(1)
		{
		SubBase = RGseeds.IRandom(0,3);
		if(SubBase != *pSeq)
			{
			*pSeq = SubBase;
			break;
			}
		}
	};
return(SeqLen);
}


void
CSimReads::Reset(bool bSync)
{
if(m_hOutFile != -1)
	{
	if(bSync)
#ifdef _WIN32
		_commit(m_hOutFile);
#else
		fsync(m_hOutFile);
#endif
	close(m_hOutFile);
	m_hOutFile = -1;
	}
if(m_hOutPEFile != -1)
	{
	if(bSync)
#ifdef _WIN32
		_commit(m_hOutPEFile);
#else
		fsync(m_hOutPEFile);
#endif
	close(m_hOutPEFile);
	m_hOutPEFile = -1;
	}
if(m_pBioSeqFile != NULL)
	{
	delete m_pBioSeqFile;
	m_pBioSeqFile = NULL;
	}
if(m_pBEDFile != NULL)
	{
	delete m_pBEDFile;
	m_pBEDFile = NULL;
	}

if(m_pProfSel != NULL)
	{
	delete m_pProfSel;
	m_pProfSel = NULL;
	}
if(m_pProfCSV != NULL)
	{
	delete m_pProfCSV;
	m_pProfCSV = NULL;
	}

if(m_pChromSeqs != NULL)
	{
#ifdef _WIN32
	free(m_pChromSeqs);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pChromSeqs != MAP_FAILED)
		munmap(m_pChromSeqs,m_AllocdChromSize);
#endif
	m_pChromSeqs = NULL;
	}

if(m_pSimReads != NULL)
	{
#ifdef _WIN32
	free(m_pSimReads);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pSimReads != MAP_FAILED)
		munmap(m_pSimReads,m_AllocdMemReads);
#endif
	m_pSimReads = NULL;
	}

if(m_pGenomeSeq != NULL)
	{
#ifdef _WIN32
	free(m_pGenomeSeq);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pGenomeSeq != MAP_FAILED)
		munmap(m_pGenomeSeq,m_AllocdGenomeMem);
#endif
	m_pGenomeSeq = NULL;
	}

m_AllocdMemReads = 0;
m_NumReadsAllocd = 0;

m_GenomeLen = 0;
m_AllocdGenomeMem = 0;
m_GenomeScaledLen = 0;
m_szSpecies[0] = '\0';
m_NumChromSeqs = 0;
m_AllocdChromSeqs = 0;
m_AllocdChromSize = 0;



m_DynProbErr = 0.01;	// default as being a 1% error rate
m_bUniformDist = false;		// default as being a uniform distribution of induced errors
m_InDelSize = 3;				// simulated InDel size range
m_InDelRate = 0.0;				// simulated InDel rate per read
m_PrvProfSeqLen = -1;  // set to read length for which DynErrProfile has been generated
m_TotReqReads = 0;
m_CurNumGenReads = 0;
m_MaxFastaLineLen = 79;
memset(m_InducedErrDist,0,sizeof(m_InducedErrDist));	// to hold induced error count distribution
memset(m_InducedErrPsnDist,0,sizeof(m_InducedErrPsnDist));	// to hold read sequence psn induced error count
}

// generate '+' strand index from K-mer of length SeqLen starting at pSeq
int
CSimReads::GenPSeqIdx(int SeqLen,etSeqBase *pSeq)
{
int Idx;
int SeqIdx;
int Base;
for(Idx=SeqIdx=0; Idx < SeqLen; Idx++,pSeq++)
	{
	Base = *pSeq & NUCONLYMSK;
	if(Base > eBaseT)
		return(-1);
	SeqIdx <<= 2;
	SeqIdx |= Base;
	}
return(SeqIdx);
}

// generate '-' strand index from K-mer of length SeqLen starting at pSeq
int
CSimReads::GenMSeqIdx(int SeqLen,etSeqBase *pSeq)
{
int Idx;
int SeqIdx;
int Base;

pSeq += SeqLen-1;
for(Idx=SeqIdx=0; Idx < SeqLen; Idx++,pSeq--)
	{
	Base = *pSeq & NUCONLYMSK;
	if(Base > eBaseT)
		return(-1);
	switch(Base) {
		case 0:
			Base = 3;
			break;
		case 1:
			Base = 2;
			break;
		case 2:
			Base = 1;
			break;
		case 3:
			Base = 0;
			break;
		}
	SeqIdx <<= 2;
	SeqIdx |= Base;
	}
return(SeqIdx);
}

int
CSimReads::InitProfSitePrefs(char *pszInProfFile)	// read from this profile site selectivity file (for MNase, generated by MNaseSitePred process), if NULL then static profiling
{
int Rslt;
int NumProcessed;
int NumFields;
int OctIdx;
char *pszOctamer;
etSeqBase Octamer[9];
double SitePref;

if(pszInProfFile == NULL || pszInProfFile[0] == '\0')
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Static profile site selection preferences based on site octamer");
else
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading profile site selection preferences from file: %s",pszInProfFile);

m_pProfSel = (double *) new double [0x010000];	// to hold 4^8 (octamer) site preferences
if(m_pProfSel == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Memory allocation for profile site preferences of 64k doubles failed");
	Reset(false);
	return(eBSFerrMem);
	}

if(pszInProfFile == NULL || pszInProfFile[0] == '\0')
	{
	for(OctIdx = 0; OctIdx < 0x010000; OctIdx++)
		m_pProfSel[OctIdx] = max(0.0000001f,(double)OctIdx/(double)0x0ffff);
	return(eBSFSuccess);
	}

for(OctIdx = 0; OctIdx < 0x010000; OctIdx++)
	m_pProfSel[OctIdx] = 0.0f;

if((m_pProfCSV = new CCSVFile) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CCSVfile");
	Reset(false);
	return(eBSFerrObj);
	}

if((Rslt=m_pProfCSV->Open(pszInProfFile))!=eBSFSuccess)
	{
	while(m_pProfCSV->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pProfCSV->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open file: %s",pszInProfFile);
	Reset(false);
	return(Rslt);
	}

NumProcessed = 0;
while((Rslt=m_pProfCSV->NextLine()) > 0)	// onto next line containing fields
	{
	NumFields = m_pProfCSV->GetCurFields();
	if(NumFields < 4)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"file: %s contains % fields, expected at least 4",pszInProfFile,NumFields);
		Reset(false);
		return(eBSFerrParams);
		}
	if(!NumProcessed && m_pProfCSV->IsLikelyHeaderLine())
		continue;
	NumProcessed += 1;
	m_pProfCSV->GetText(1,&pszOctamer);
	memcpy(Octamer,CSeqTrans::MapAscii2Sense(pszOctamer),8);
	OctIdx = GenPSeqIdx(8,Octamer);
	m_pProfCSV->GetDouble(4,&SitePref);
	m_pProfSel[OctIdx] = SitePref;
	}

delete m_pProfCSV;
m_pProfCSV = NULL;

return(eBSFSuccess);
}







int
CSimReads::LoadFasta(size_t *pTotLen,int MinChromLen,char *pszFastaFile)
{
CFasta Fasta;
unsigned char *pSeqBuff;
char szName[cBSFSourceSize];
char szDescription[cBSFDescriptionSize];
int SeqLen;
int Descrlen;
bool bFirstEntry;
bool bEntryCreated;
int Rslt;
int SeqID;
tsSRChromSeq *pChromSeq;
UINT8 *pSeq;
int ChromID;
int SeqOfs;
int NumLenWarnings;
size_t TotLen;
tsSRChromSeq *pTmp;
size_t memreq;

m_GenomeLen = 0;
m_AllocdGenomeMem = 0;
TotLen = 0;
m_NumChromSeqs = 0;
m_GenomeLen = 0;
ChromID = 0;

m_AllocdChromSize = sizeof(tsSRChromSeq) * cSRAllocNumChroms;
#ifdef _WIN32
m_pChromSeqs = (tsSRChromSeq *) malloc((size_t)m_AllocdChromSize);
if(m_pChromSeqs == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadFasta: Memory allocation of %lld bytes failed",(INT64)m_AllocdChromSize);
	m_AllocdChromSize = 0;
	Reset(false);
	return(eBSFerrMem);
	}
#else
// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
m_pChromSeqs = (tsSRChromSeq *)mmap(NULL,(size_t)m_AllocdChromSize, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
if(m_pChromSeqs == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadFasta: Memory allocation of %lld bytes through mmap()  failed",(INT64)m_AllocdChromSize,strerror(errno));
	m_pChromSeqs = NULL;
	m_AllocdChromSize = 0;
	Reset(false);
	return(eBSFerrMem);
	}
#endif
m_AllocdChromSeqs = cSRAllocNumChroms;
memset(m_pChromSeqs,0,m_AllocdChromSize);

m_AllocdGenomeMem = cSRMaxAllocBuffChunk;		// an initial allocation to hold the assembly sequence, will be extended as may be required
#ifdef _WIN32
m_pGenomeSeq = (UINT8 *) malloc((size_t)m_AllocdGenomeMem);
if(m_pGenomeSeq == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadFasta: Memory allocation of %lld bytes failed",(INT64)m_AllocdGenomeMem);
	Reset(false);
	return(eBSFerrMem);
	}
#else
// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
m_pGenomeSeq = (UINT8 *)mmap(NULL,(size_t)m_AllocdGenomeMem, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
if(m_pGenomeSeq == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadFasta: Memory allocation of %lld bytes through mmap()  failed",(INT64)m_AllocdGenomeMem,strerror(errno));
	m_pGenomeSeq = NULL;
	Reset(false);
	return(eBSFerrMem);
	}
#endif


if((pSeqBuff = new unsigned char [cSRMaxAllocBuffChunk]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile:- Unable to allocate memory (%d bytes) for sequence buffer",cSRMaxAllocBuffChunk);
	Reset(false);
	return(eBSFerrMem);
	}

if((Rslt=Fasta.Open(pszFastaFile,true))!=eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open '%s' [%s] %s",pszFastaFile,Fasta.ErrText((teBSFrsltCodes)Rslt),Fasta.GetErrMsg());
	Reset(false);
	delete pSeqBuff;
	return(Rslt);
	}
bFirstEntry = true;
bEntryCreated = false;
SeqID = 0;
ChromID = 0;
TotLen = 0;

NumLenWarnings = 0;
while((Rslt = SeqLen = Fasta.ReadSequence(pSeqBuff,cSRMaxAllocBuffChunk-1,true,false)) > eBSFSuccess)
	{
	if(bFirstEntry || SeqLen == eBSFFastaDescr)		// just read a descriptor line
		{
		if(!bFirstEntry)
			{
			pChromSeq->Len = SeqOfs;
			if(pChromSeq->Len < MinChromLen && NumLenWarnings++ < 10)
				gDiagnostics.DiagOut(eDLWarn,gszProcName,"LoadFasta: Chrom '%s' of length %d may be too short from which to reliably sample simulated reads...",pChromSeq->szChromName,pChromSeq->Len);
			m_pGenomeSeq[pChromSeq->SeqOfs+SeqOfs] = eBaseEOS;
			m_GenomeLen += (size_t)SeqOfs + 1;
			}

		if(ChromID == m_AllocdChromSeqs)
			{
			memreq = sizeof(tsSRChromSeq) * ((size_t)m_AllocdChromSeqs + cSRAllocNumChroms);
#ifdef _WIN32
			pTmp = (tsSRChromSeq *) realloc(m_pChromSeqs,memreq);
#else
			pTmp = (tsSRChromSeq *)mremap(m_pChromSeqs,m_AllocdChromSize,memreq,MREMAP_MAYMOVE);
			if(pTmp == MAP_FAILED)
				pTmp = NULL;
#endif
			if(pTmp == NULL)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadFasta: Memory re-allocation to %lld bytes - %s",(INT64)(memreq),strerror(errno));
				delete pSeqBuff;
				return(eBSFerrMem);
				}
			m_pChromSeqs = pTmp;
			m_AllocdChromSize = memreq;
			m_AllocdChromSeqs += cSRAllocNumChroms;
			}

     	pChromSeq = &m_pChromSeqs[ChromID++];
		pChromSeq->ChromSeqID = ChromID;
		pChromSeq->ChromID = ChromID;
		m_NumChromSeqs = ChromID;
		Descrlen = Fasta.ReadDescriptor(szDescription,cBSFDescriptionSize);
		if(SeqLen != eBSFFastaDescr || sscanf(szDescription," %s[ ,]",szName)!=1)
			sprintf(szName,"%s.%d",pszFastaFile,ChromID);
		strncpy(pChromSeq->szChromName,szName,sizeof(pChromSeq->szChromName));
		pChromSeq->szChromName[50] = '\0';	  // limiting chrom names to be at most 50 chars as additional detail will later be appended		
		pChromSeq->Strand = '*';
		pChromSeq->RelDensity = 1000;
		pChromSeq->SeqOfs = (UINT32)m_GenomeLen;
		bEntryCreated = true;
		bFirstEntry = false;
		SeqOfs = 0;
		if(SeqLen == eBSFFastaDescr)
			continue;
		}

	// realloc m_pGenomeSeq as and when required
	if((m_GenomeLen + SeqOfs + SeqLen + 100) >= (INT64)m_AllocdGenomeMem)
		{
		memreq = m_AllocdGenomeMem + SeqLen + cSRMaxAllocBuffChunk;
#ifdef _WIN32
		pSeq = (UINT8 *) realloc(m_pGenomeSeq,memreq);
#else
		pSeq = (UINT8 *)mremap(m_pGenomeSeq,m_AllocdGenomeMem,memreq,MREMAP_MAYMOVE);
		if(pSeq == MAP_FAILED)
			pSeq = NULL;
#endif
		if(pSeq == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadFasta: Memory re-allocation to %lld bytes - %s",(INT64)(memreq),strerror(errno));
			delete pSeqBuff;
			return(eBSFerrMem);
			}
		m_pGenomeSeq = pSeq;
		m_AllocdGenomeMem = memreq;
		}

	pSeq = &m_pGenomeSeq[pChromSeq->SeqOfs+SeqOfs];
	memcpy(pSeq,pSeqBuff,SeqLen);
	SeqOfs += SeqLen;
	TotLen += SeqLen;
	}

if(Rslt < eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile [%s] %s",Fasta.ErrText((teBSFrsltCodes)Rslt),Fasta.GetErrMsg());
	delete pSeqBuff;
	return(false);
	}
if(bEntryCreated)
	{
	pChromSeq->Len = SeqOfs;
	if(pChromSeq->Len < MinChromLen && NumLenWarnings++ < 10)
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"LoadFasta: Chrom '%s' of length %d may be too short from which to reliably sample simulated reads...",pChromSeq->szChromName,pChromSeq->Len);

	m_pGenomeSeq[pChromSeq->SeqOfs+SeqOfs] = eBaseEOS;
	m_GenomeLen += SeqOfs + 1;
	}
delete pSeqBuff;
*pTotLen = TotLen;
return(eBSFSuccess);
}

int
CSimReads::LoadBioseq(size_t *pTotLen,int MinChromLen,char *pszBioSeqFile)
{
int Rslt;
size_t TotLen;
int ChromID;
etSeqBase *pSeq;
tsSRChromSeq *pChromSeq;
int Len;
int NumLenWarnings;

*pTotLen = 0;
if((m_pBioSeqFile = new CBioSeqFile) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBioSeqFile object");
	return(eBSFerrObj);
	}
if((Rslt = m_pBioSeqFile->Open(pszBioSeqFile))!=eBSFSuccess)
	{
	if(Rslt == eBSFerrFileAccess)
		{
		delete m_pBioSeqFile;
		m_pBioSeqFile = NULL;
		return(Rslt);
		}

	while(m_pBioSeqFile->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pBioSeqFile->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open assembly sequence file '%s'",pszBioSeqFile);
	delete m_pBioSeqFile;
	m_pBioSeqFile = NULL;
	Reset(false);
	return(Rslt);
	}

m_pBioSeqFile->GetTitle(sizeof(m_szSpecies),m_szSpecies);

m_GenomeLen = 0;
m_AllocdGenomeMem = 0;
TotLen = 0;
m_NumChromSeqs = 0;
ChromID = 0;
m_AllocdChromSeqs = m_pBioSeqFile->NumEntries();

m_AllocdChromSize = sizeof(tsSRChromSeq) * m_AllocdChromSeqs;
#ifdef _WIN32
m_pChromSeqs = (tsSRChromSeq *) malloc((size_t)m_AllocdChromSize);
if(m_pChromSeqs == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadAssembly: Memory allocation of %lld bytes failed",(INT64)m_AllocdChromSize);
	Reset(false);
	return(eBSFerrMem);
	}
#else
// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
m_pChromSeqs = (tsSRChromSeq *)mmap(NULL,(size_t)m_AllocdChromSize, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
if(m_pChromSeqs == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadAssembly: Memory allocation of %lld bytes through mmap()  failed",(INT64)m_AllocdChromSize,strerror(errno));
	m_pChromSeqs = NULL;
	Reset(false);
	return(eBSFerrMem);
	}
#endif

memset(m_pChromSeqs,0,m_AllocdChromSize);

// determine total sequence length required for allocating to hold genome as one concatenated sequence
NumLenWarnings = 0;
while((ChromID = m_pBioSeqFile->Next(ChromID)) > 0) {
	pChromSeq = &m_pChromSeqs[ChromID-1];
	pChromSeq->ChromID = ChromID;
	m_pBioSeqFile->GetName(ChromID,sizeof(pChromSeq->szChromName),pChromSeq->szChromName);
	pChromSeq->Strand = '*';
	pChromSeq->RelDensity = 1000;
	Len = m_pBioSeqFile->GetDataLen(ChromID);
	if(Len < MinChromLen && NumLenWarnings++ < 10)
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"LoadAssembly: Chrom '%s' of length %d may be too short from which to reliably sample simulated reads...",pChromSeq->szChromName,Len);
	m_GenomeLen += (INT64)Len + 1;		// allow for a concatenation separator
	pChromSeq->Len = Len;
	pChromSeq->SeqOfs = 0;
	pChromSeq->ChromSeqID = ++m_NumChromSeqs;
	}
m_AllocdGenomeMem = (size_t)m_GenomeLen;
#ifdef _WIN32
m_pGenomeSeq = (UINT8 *) malloc((size_t)m_AllocdGenomeMem);
if(m_pGenomeSeq == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadAssembly: Memory allocation of %lld bytes failed",(INT64)m_AllocdGenomeMem);
	Reset(false);
	return(eBSFerrMem);
	}
#else
// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
m_pGenomeSeq = (UINT8 *)mmap(NULL,(size_t)m_AllocdGenomeMem, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
if(m_pGenomeSeq == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadAssembly: Memory allocation of %lld bytes through mmap()  failed",(INT64)m_AllocdGenomeMem,strerror(errno));
	m_pGenomeSeq = NULL;
	Reset(false);
	return(eBSFerrMem);
	}
#endif

ChromID = 0;
TotLen = 0;
pSeq = m_pGenomeSeq;
while((ChromID = m_pBioSeqFile->Next(ChromID)) > 0) {
	pChromSeq = &m_pChromSeqs[ChromID-1];
	pChromSeq->SeqOfs = (UINT32)(UINT64)(pSeq - m_pGenomeSeq);
	TotLen += pChromSeq->Len;
	m_pBioSeqFile->GetData(ChromID,eSeqBaseType,0,pSeq,pChromSeq->Len);
	pSeq += pChromSeq->Len;
	*pSeq++ = eBaseEOS;		// concatenator separator
	}
delete m_pBioSeqFile;
m_pBioSeqFile = NULL;
*pTotLen = TotLen;
return(eBSFSuccess);
}

int
CSimReads::LoadGenome(int MinChromLen,			// warn if loaded chromosomes are less than this length; may be too short to sample reads from
			char *pszBioSeqFile)
{
int Rslt;
size_t TotLen;
int ChromID;
tsSRChromSeq *pChromSeq;
double LenSCF;			// length scaling factor
int CurScaledStart;

// assume a bioseq file, if not a bioseq file then try to load as a raw multifasta
if((Rslt = LoadBioseq(&TotLen,MinChromLen,pszBioSeqFile)) != eBSFSuccess)
	{
	if(Rslt != eBSFerrFileAccess)
		return(Rslt);

	if((Rslt = LoadFasta(&TotLen,MinChromLen,pszBioSeqFile)) != eBSFSuccess)
		return(Rslt);
	}

// assembly or multifasta sequences now loaded
LenSCF = (double)INT_MAX/(double)TotLen;
m_GenomeScaledLen = TotLen >= (INT64)INT_MAX ? (long)(TotLen * LenSCF) : (long)TotLen;
CurScaledStart = 0;
pChromSeq = &m_pChromSeqs[0];
for(ChromID = 0; ChromID < m_NumChromSeqs; ChromID++,pChromSeq++)
	{
	pChromSeq->ScaledLen = TotLen >= (INT64)INT_MAX ? (int)(pChromSeq->Len * LenSCF) : pChromSeq->Len;
	pChromSeq->ScaledStart = CurScaledStart;
	CurScaledStart += pChromSeq->ScaledLen;
	}
return(eBSFSuccess);
}


int
CSimReads::LoadTranscriptome(char *pszBioSeqFile,			// genome assembly
				char *pszFeatFile)				// BED file containing features or genes
{
int Rslt;
bool bTranscribed;
int Len;
size_t TotLen;
int ChromID;
int FeatID;
tsSRChromSeq *pChromSeq;
etSeqBase *pSeq;
double LenSCF;			// length scaling factor
int CurScaledStart;

if((m_pBioSeqFile = new CBioSeqFile) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBioSeqFile object");
	return(eBSFerrObj);
	}
if((Rslt = m_pBioSeqFile->Open(pszBioSeqFile))!=eBSFSuccess)
	{
	while(m_pBioSeqFile->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pBioSeqFile->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open assembly sequence file '%s'",pszBioSeqFile);
	Reset(false);
	return(Rslt);
	}

m_pBioSeqFile->GetTitle(sizeof(m_szSpecies),m_szSpecies);

// open feature file
if((m_pBEDFile = new CBEDfile) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBEDfile object");
	return(eBSFerrObj);
	}
if((Rslt = m_pBEDFile->Open(pszFeatFile))!=eBSFSuccess)
	{
	while(m_pBEDFile->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pBEDFile->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open feature/gene BED file '%s'",pszFeatFile);
	Reset(false);
	return(Rslt);
	}

// check if BED contains features as genes (transcribed lengths) or simply features as being regions in the genome with no transcriptional information
bTranscribed = m_pBEDFile->ContainsGeneDetail();

// alloc for total number of features, sum of transcribed length if genes or sum of feature lengths if regions
m_AllocdChromSeqs = m_pBEDFile->GetNumFeatures();
if((m_pChromSeqs = new tsSRChromSeq [m_AllocdChromSeqs])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadTranscriptome: unable to allocate memory (%d bytes) for %d Feature Seqs",sizeof(tsSRChromSeq) * m_AllocdChromSeqs,m_AllocdChromSeqs);
	Reset(false);
	return(eBSFerrMem);
	}
memset(m_pChromSeqs,0,sizeof(tsSRChromSeq) * m_AllocdChromSeqs);

// determine total sequence length required for allocating to hold all features or genes as one concatenated sequence
m_GenomeLen = 0;
TotLen = 0;
m_NumChromSeqs = 0;

int CurFeatScore;
int MaxFeatScore;
char CurStrand;

MaxFeatScore = 0;
FeatID = 0;
while((FeatID = m_pBEDFile->GetNextFeatureID(FeatID)) > 0) {
	pChromSeq = &m_pChromSeqs[FeatID-1];
	pChromSeq->ChromID = FeatID;
	m_pBEDFile->GetFeature(FeatID,pChromSeq->szChromName,NULL,NULL,NULL,&CurFeatScore,&CurStrand);
	pChromSeq->Strand = CurStrand;
	CurFeatScore = min(CurFeatScore,1000);
	MaxFeatScore = max(MaxFeatScore,CurFeatScore);
	if(bTranscribed)
		Len = m_pBEDFile->GetTranscribedLen(FeatID);			// get transcribed length for gene
	else
		Len = m_pBEDFile->GetFeatLen(FeatID);
	m_GenomeLen += (INT64)Len + 1;		// allow for a concatenation separator
	pChromSeq->Len = Len;
	pChromSeq->SeqOfs = 0;
	pChromSeq->ChromSeqID = ++m_NumChromSeqs;
	}

FeatID = 0;
while((FeatID = m_pBEDFile->GetNextFeatureID(FeatID)) > 0) {
	pChromSeq = &m_pChromSeqs[FeatID-1];
	m_pBEDFile->GetFeature(FeatID,NULL,NULL,NULL,NULL,&CurFeatScore);
	CurFeatScore = min(1000,CurFeatScore);
	if(MaxFeatScore == 0)			// if no scores for any feature then assume all are to be equally sampled
		pChromSeq->RelDensity = 1000;
	else
		pChromSeq->RelDensity = CurFeatScore == 0 ? 0 : 1 + ((1000 * CurFeatScore) / MaxFeatScore);
	}

// now know the feature sequence lengths, alloc to hold all sequences
#ifdef _WIN32
m_pGenomeSeq = (UINT8 *) malloc((size_t)m_GenomeLen);
if(m_pGenomeSeq == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadTranscriptome: Memory allocation of %lld bytes failed",(INT64)m_GenomeLen);
	Reset(false);
	return(eBSFerrMem);
	}
#else
// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
m_pGenomeSeq = (UINT8 *)mmap(NULL,(size_t)m_GenomeLen, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
if(m_pGenomeSeq == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadAssembly: Memory allocation of %lld bytes through mmap()  failed",(INT64)m_GenomeLen,strerror(errno));
	m_pGenomeSeq = NULL;
	Reset(false);
	return(eBSFerrMem);
	}
#endif

int CurRegionLen;
int NumExons;
int StartLoci;
int EndLoci;
int Idx;
int ExonLen;
int NumWarns = 0;
int CurFeatureID = 0;
int PrevFeatureChromID = 0;
TotLen = 0;
pSeq = m_pGenomeSeq;

char szGenomeChromName[128];
int FeatureChromID;
int GenomeChromID;

while((CurFeatureID = m_pBEDFile->GetNextFeatureID(CurFeatureID)) > 0) {
	FeatureChromID = m_pBEDFile->GetFeatureChromID(CurFeatureID);
	if(FeatureChromID != PrevFeatureChromID)
		{
		m_pBEDFile->GetChromosome(FeatureChromID,szGenomeChromName);
		if((GenomeChromID = m_pBioSeqFile->LocateEntryIDbyName(szGenomeChromName)) <= 0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadTranscriptome: Unable to load locate feature '%s' chrom '%s'",pChromSeq->szChromName,szGenomeChromName);
			Reset(false);
			return(eBSFerrChrom);
			}
		PrevFeatureChromID = FeatureChromID;
		}

	pChromSeq = &m_pChromSeqs[CurFeatureID-1];
	pChromSeq->SeqOfs = (UINT32)(UINT64)(pSeq - m_pGenomeSeq);
	TotLen += pChromSeq->Len;
	CurRegionLen = 0;
	if(pChromSeq->Len > 0)
		{
		if(bTranscribed)
			NumExons = m_pBEDFile->GetNumExons(CurFeatureID);
		else
			NumExons = 1;
		for(Idx = 1; Idx <= NumExons; Idx++)
			{
			StartLoci = m_pBEDFile->GetExonStart(CurFeatureID,Idx);
			EndLoci   = m_pBEDFile->GetExonEnd(CurFeatureID,Idx);
			ExonLen = 1 + EndLoci - StartLoci;
			CurRegionLen += ExonLen;
			if((Rslt=m_pBioSeqFile->GetData(GenomeChromID,eSeqBaseType,StartLoci,pSeq,ExonLen))<=0)
				{
				if(NumWarns++ < 10)
					gDiagnostics.DiagOut(eDLWarn,gszProcName,"LoadTranscriptome: Unable to load sequence for feature '%s' from chrom '%s'",pChromSeq->szChromName,szGenomeChromName);
				else
					{
					Reset(false);
					return(eBSFerrMem);
					}
				}
			pSeq += ExonLen;
			}
		}

	if(pChromSeq->Len != CurRegionLen)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadTranscriptome: Expected sum of exon lengths (%d) to equal transcript length (%d) for element/gene '%s'",CurRegionLen,pChromSeq->Len,pChromSeq->szChromName);
		Reset(false);
		return(eBSFerrInternal);
		}
	*pSeq++ = eBaseEOS;		// concatenator separator
	}

delete m_pBioSeqFile;
m_pBioSeqFile = NULL;

LenSCF = (double)INT_MAX/(double)TotLen;
m_GenomeScaledLen = TotLen >= (INT64)INT_MAX ? (long)(TotLen * LenSCF) : (long)TotLen;
CurScaledStart = 0;
pChromSeq = &m_pChromSeqs[0];
for(ChromID = 0; ChromID < m_NumChromSeqs; ChromID++,pChromSeq++)
	{
	pChromSeq->ScaledLen = TotLen >= (INT64)INT_MAX ? (int)(pChromSeq->Len * LenSCF) : pChromSeq->Len;
	pChromSeq->ScaledStart = CurScaledStart;
	CurScaledStart += pChromSeq->ScaledLen;
	}
return(eBSFSuccess);
}


// SimulateSNPs
// Simulate SNPs in the genome at the specified rate per million bases
// Reports their loci
int
CSimReads::SimulateSNPs(char *pszOutSNPs,   // output simulated SNP loci to this file
			int SNPrate)		 // required SNPs per million bases
{
int hFile;
int BuffOfs;
char szLineBuff[8196];
char szChromName[128];
int Ofs;
int SNPiD;
int NumChromSNPs;
etSeqBase *pSeq;
etSeqBase PrevBase;
etSeqBase SNPbase;

int ChromID;
tChromID AssembChromID;
int ChromLoci;

tsSRChromSeq *pChromSeq;

if(SNPrate == 0)
	return(eBSFSuccess);

if(pszOutSNPs != NULL && pszOutSNPs[0] != '\0')
	{
#ifdef _WIN32
	hFile = open(pszOutSNPs,O_CREATETRUNC );
#else
	if((hFile = open(pszOutSNPs,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
		if(ftruncate(hFile,0)!=0)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",pszOutSNPs,strerror(errno));
				Reset(false);
				return(eBSFerrCreateFile);
				}
#endif

	if(hFile < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output file '%s'",pszOutSNPs);
		Reset(false);
		return(eBSFerrCreateFile);
		}
	BuffOfs = sprintf(szLineBuff,"track type=bed name=\"SimSNPs\" description=\"Simulated SNPS\"\n");
	CUtility::SafeWrite(hFile,szLineBuff,BuffOfs);
	}
else
	hFile = -1;

TRandomCombined<CRandomMother,CRandomMersenne> RG((int)time(0));
BuffOfs = 0;
SNPiD = 0;
for(ChromID = 0; ChromID < m_NumChromSeqs; ChromID++,pChromSeq++)
	{
	pChromSeq = &m_pChromSeqs[ChromID];
	pSeq = &m_pGenomeSeq[pChromSeq->SeqOfs];
	NumChromSNPs = 1 + (int)(pChromSeq->Len * ((double)SNPrate/1000000));
	while(NumChromSNPs)
		{
		Ofs = (int)RG.IRandom(0,pChromSeq->Len-1);
		pSeq = &m_pGenomeSeq[pChromSeq->SeqOfs+Ofs];
		if((PrevBase = *pSeq) > eBaseT)
			continue;
		if((SNPbase = (int)RG.IRandom(0,3))==PrevBase)
			continue;
		*pSeq = SNPbase;
		NumChromSNPs -= 1;

		if(hFile != -1)
			{
			// if transcriptome read generation then remap the transcript relative SNP loci back to the assembly chrom loci
			if(m_pBEDFile != NULL)
				{
				m_pBEDFile->MapTransOfs2Loci(pChromSeq->ChromID,Ofs,NULL,&AssembChromID,&ChromLoci);
				m_pBEDFile->GetChromosome(AssembChromID,szChromName);
				}
			else
				{
				ChromLoci = Ofs;
				strcpy(szChromName,pChromSeq->szChromName);
				}

			SNPiD += 1;
			BuffOfs+=sprintf(&szLineBuff[BuffOfs],"%s\t%d\t%d\tSNP_%d\t%d\t+\n",
					szChromName,ChromLoci,ChromLoci+1,SNPiD,0);

			if(BuffOfs + 200 > sizeof(szLineBuff))
				{
				CUtility::SafeWrite(hFile,szLineBuff,BuffOfs);
				BuffOfs = 0;
				}
			}
		}
	}


if(hFile != -1 && BuffOfs > 0)
	{
	CUtility::SafeWrite(hFile,szLineBuff,BuffOfs);
	}
if(hFile != -1)
	close(hFile);
return(eBSFSuccess);
}

// CmpSeqs
// Returns the number of differences between two sequences of length Len
int
CSimReads::CmpSeqs(int Len,etSeqBase *pSeq1,etSeqBase *pSeq2)
{
int Diffs = 0;
etSeqBase Base1;
etSeqBase Base2;
int Idx;
for(Idx=0; Idx < Len; Idx++)
	if((Base1 = (*pSeq1++ & NUCONLYMSK)) != (Base2 = (*pSeq2++ & NUCONLYMSK)))
		Diffs += 1;
return(Diffs);
}

bool
CSimReads::IsDupSeq(int Len,etSeqBase *pSeq1,etSeqBase *pSeq2)
{
return(CmpSeqs(Len,pSeq1,pSeq2) == 0 ? true : false);
}



// both the forward and reverse complements are compared for equivilence
bool
CSimReads::IsDupSeqStrands(int Len,etSeqBase *pReadSeq1,etSeqBase *pReadSeq2)
{
etSeqBase *pSeq1 = pReadSeq1;
etSeqBase *pSeq2 = pReadSeq2;
etSeqBase Base1;
etSeqBase Base2;
int Idx;

// first assume both reads will ultimately align onto the forward (watson) strand
pSeq1 = pReadSeq1;
pSeq2 = pReadSeq2;
for(Idx=0; Idx < Len; Idx++)
	if((Base1 = (*pSeq1++ & NUCONLYMSK)) != (Base2 = (*pSeq2++ & NUCONLYMSK)))
		break;
if(Idx == Len)
	return(true);

// now assume that pReadSeq1 stays watson but that pReadSeq2 will be reverse complemented to align onto the crick strand
pSeq1 = pReadSeq1;
pSeq2 = &pReadSeq2[Len-1];
for(Idx=0; Idx < Len; Idx++)
	{
	Base1 = (*pSeq1++ & NUCONLYMSK);
	Base2 = (*pSeq2-- & NUCONLYMSK);
	switch(Base2) {
		case eBaseA: Base2=eBaseT; break;
		case eBaseC: Base2=eBaseG; break;
		case eBaseG: Base2=eBaseC; break;
		case eBaseT: Base2=eBaseA; break;
		}
	if(Base1 != Base2)
		break;
	}
if(Idx == Len)
	return(true);

// finally assume that pReadSeq2 stays watson but that pReadSeq1 will be reverse complemented to align onto the crick strand
pSeq1 = &pReadSeq1[Len-1];
pSeq2 = pReadSeq2;
for(Idx=0; Idx < Len; Idx++)
	{
	Base1 = (*pSeq1-- & NUCONLYMSK);
	Base2 = (*pSeq2++ & NUCONLYMSK);
	switch(Base1) {
		case eBaseA: Base1=eBaseT; break;
		case eBaseC: Base1=eBaseG; break;
		case eBaseG: Base1=eBaseC; break;
		case eBaseT: Base1=eBaseA; break;
		}
	if(Base1 != Base2)
		return(false);
	}
return(true);
}

//
void
CSimReads::ShowDupSeqs(int ReadLen,tsSRSimRead *pRead1,tsSRSimRead *pRead2)
{
int Idx;
etSeqBase *pSeq;
UINT8 Seq[2000];
char szSeq[2000];
pSeq = pRead1->pSeq;
for(Idx = 0; Idx < ReadLen; Idx++,pSeq++)
	Seq[Idx] = *pSeq & NUCONLYMSK;
CSeqTrans::MapSeq2Ascii(Seq,ReadLen,szSeq);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"chrom: %1.2d, loci: %1.9d seq: '%s'",
			 pRead1->ChromID,pRead1->StartLoci,szSeq);
pSeq = pRead2->pSeq;
for(Idx = 0; Idx < ReadLen; Idx++,pSeq++)
	Seq[Idx] = *pSeq & NUCONLYMSK;
CSeqTrans::MapSeq2Ascii(Seq,ReadLen,szSeq);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"chrom: %1.2d, loci: %1.9d seq: '%s'",
			 pRead2->ChromID,pRead2->StartLoci,szSeq);
}



int
CSimReads::ReportReads(bool bPEgen,	// true if paired end simulated reads being simulated
		    int Region,				// Process regions 0:ALL,1:CDS,2:5'UTR,3:3'UTR,4:Introns,5:5'US,6:3'DS,7:Intergenic (default = ALL)
			int NomReadLen,			// reads nominal length
			int NumPrevReported,	// number of reads previously reported
			int MaxReads,			// maximum number of reads remaining to be reported
			bool bDedupe,			// true if reads are to be deduped
			int AvailReads)			// number of reads to be reported on from m_pSimReads[]
{
tsSRSimRead *pRead;
tsSRSimRead *pPairRead;
tsSRSimRead *pPrevRead;
tsSRChromSeq *pChromSeq;
etSeqBase *pReadSeq;
char *pszIsRand;
int ReadOfs;
int ReadsIdx;
int ReadLenRem;

int NumOfDups;
int NumReported;

char szCIGAR[128];				// to hold a simulated reads SAM CIGAR

sReadRprt *pPEreads;			// to hold each either SE or PE read being reported

szCIGAR[0] = '\0';


// check if any reads to report
if(m_pSimReads == NULL || MaxReads < 1 || AvailReads < 1)
	return(0);

if((pPEreads = new sReadRprt[2])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ReportReads: Memory allocation error");
	return(eBSFerrMem);
	}

if(bDedupe && AvailReads > 1)
	{
	NumOfDups = 0;
	// for detection of the duplicate read sequences, the sequences need to be sorted
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Deduping %d read sequences, sorting...",AvailReads);
	qsort(m_pSimReads,AvailReads,sizeof(tsSRSimRead),SortSimReads);
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Sorting completed, checking for duplicates");
	// firstly find and mark all watson '+' strand duplicates
	pPrevRead = m_pSimReads;
	pRead = &m_pSimReads[1];
	for(ReadsIdx = 1; ReadsIdx < AvailReads; ReadsIdx++, pRead++)
		{
		if(pPrevRead->Len == pRead->Len && IsDupSeq(pRead->Len,pPrevRead->pSeq,pRead->pSeq))
			{
			if(pRead->Status > 1)
				pPrevRead->Status = 2;
			pRead->Status = 2;
			if(pRead->pPartner != NULL)
				pRead->pPartner->Status = 2;
			NumOfDups+=1;
			}
		else
			{
			if(pPrevRead->Len == pRead->Len)
				{
				int Diffs = CmpSeqs(pRead->Len,pPrevRead->pSeq,pRead->pSeq);
				if(Diffs == 0)
					{
					if(pRead->Status > 1)
						pPrevRead->Status = 2;
					pRead->Status = 2;		// currently treat these as if duplicates
					if(pRead->pPartner != NULL)
						pRead->pPartner->Status = 2;
					NumOfDups+=1;
					}
				}
			pPrevRead = pRead;
			}
		}
	// next, for those sequences not already marked as duplicates, see if after revcpl into the crick they match one still on watson
	// iterate and revcpl each non-dup sequence
	// search for this crick sequence
	// if none found this is a truely unique sequence
	// if found then mark all instances as being duplicates
	int RevIdx;
	pRead = m_pSimReads;
	for(ReadsIdx = 0; ReadsIdx < AvailReads; ReadsIdx++, pRead++)
		{
		if((RevIdx=LocateRevCplSeq(pRead->Len,pRead->pSeq,AvailReads)) >= 0)
			{
			ShowDupSeqs(pRead->Len,&m_pSimReads[RevIdx],pRead);
			pRead->Status = 2;
			if(pRead->pPartner != NULL)
				pRead->pPartner->Status = 2;
			NumOfDups+=1;
			}
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Marked %d reads out of %d as duplicate sequences",NumOfDups,AvailReads);
	}

// if transcriptome read generation then remap the transcript relative loci back to the assembly chrom loci
if(Region == eSRMEGRAny && m_pBEDFile != NULL)
	{
	tChromID ChromID;
	int ChromLoci;

	pRead = m_pSimReads;
	NumReported = 0;
	for(ReadsIdx = 0; ReadsIdx < AvailReads; ReadsIdx++, pRead++)
		{
		if(pRead->Status > 0)	// slough if read was determined as being a duplicate
			continue;
		pChromSeq = &m_pChromSeqs[pRead->ChromSeqID-1];
		m_pBEDFile->MapTransOfs2Loci(pRead->ChromID,pRead->StartLoci,NULL,&ChromID,&ChromLoci);
		m_pBEDFile->GetChromosome(ChromID,pChromSeq->szChromName);
		pRead->StartLoci = ChromLoci;
		m_pBEDFile->MapTransOfs2Loci(pRead->ChromID,pRead->EndLoci,NULL,NULL,&ChromLoci);
		pRead->ChromID = (int)ChromID;
		pRead->EndLoci = ChromLoci;
		}
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Writing simulated reads to file: '%s'",m_pszOutFile);


memset(pPEreads,0,2 * sizeof(sReadRprt));
pRead = m_pSimReads;
NumReported = 0;

for(ReadsIdx = 0; ReadsIdx < AvailReads; ReadsIdx++, pRead++)
	{
	if(pRead->Status > 0)	// slough if read already reported on, or if a duplicate
		continue;
	int Idx;
	int Pair1ReadLen;
	int Pair2ReadLen;

	pReadSeq = pRead->pSeq;
	Pair1ReadLen = pRead->Len;
	for(Idx = 0; Idx < (Pair1ReadLen + m_InDelSize); Idx++,pReadSeq++)	// extra read length is in case InDels being simulated and upto 10 bases deleted
		{
		if((*pReadSeq & NUCONLYMSK) > eBaseN)
			break;
		pPEreads[0].FwdSeq[Idx] = *pReadSeq & NUCONLYMSK;
		}
	if(Idx != (Pair1ReadLen + m_InDelSize))
		{
		pRead->Status = 2;
		if(bPEgen == true)
			pRead->pPartner->Status = 2;
		continue;
		}

	if(bPEgen == true)
		{
		pPairRead =  pRead->pPartner;
		Pair2ReadLen = pPairRead->Len;
		pReadSeq = pPairRead->pSeq;
		for(Idx = 0; Idx < (Pair2ReadLen + m_InDelSize); Idx++,pReadSeq++)	// extra read length is in case InDels being simulated and upto 10 bases deleted
			{
			if((*pReadSeq & NUCONLYMSK) > eBaseN)
				break;
			pPEreads[1].FwdSeq[Idx] = *pReadSeq & NUCONLYMSK;
			}
		if(Idx != (Pair2ReadLen + m_InDelSize))
			{
			pPairRead->Status = 2;
			pRead->Status = 2;
			continue;
			}
		}

	pReadSeq = pPEreads[0].FwdSeq;

	NumReported += 1;
	pChromSeq = &m_pChromSeqs[pRead->ChromSeqID-1];
	switch(m_FMode) {
		case eSRFMFasta:
		case eSRFMNWFasta:
			{
			int NumSubs;
			int InDelSize;

			if(m_FMode == eSRFMFasta)
				m_MaxFastaLineLen = 79;
			else
				m_MaxFastaLineLen = cSRMaxReadLen;

			pReadSeq = pPEreads[0].FwdSeq;
			if(pRead->Strand)
				{
				memcpy(pPEreads[0].RevSeq,pReadSeq,Pair1ReadLen+m_InDelSize);
				CSeqTrans::ReverseComplement(Pair1ReadLen+m_InDelSize,pPEreads[0].RevSeq);
				pReadSeq = pPEreads[0].RevSeq;
				}

			SimArtefacts(false,Pair1ReadLen,pReadSeq);
			SimArtefacts(true,Pair1ReadLen,pReadSeq);

			if(!m_PropRandReads || (m_PropRandReads <= RGseeds.IRandom(0,1000000)))
				{
				InDelSize = SimInDels(pRead,&Pair1ReadLen,pReadSeq);
				NumSubs = SimSeqErrors(Pair1ReadLen,pReadSeq);
				pszIsRand = (char *)"lcl";
				}
			else
				{
				InDelSize = 0;
				NumSubs = SimSeqRand(Pair1ReadLen,pReadSeq);
				pszIsRand = (char *)"lcr";
				}

			// CIGAR chars recognised/used are:
			//	=  matching
			//	X  mismatch
			//	M either matching or mismatching
			//	I  insertion relative to target - score after skipping I bases in read
			//	D  deletion relative to target - score after skipping D bases in target
			//	N  skipped region relative to target - score after skipping N bases in target
			//	S  soft clipping - score after skipping S bases in both read and target
			//	H  hard clipping - score starting from first base in sequence
			//
			// TLEN template length
			// >chr1.2003.2102.+.1234567 100:25M65M:0

			pPEreads[0].LineLen+=sprintf(&pPEreads[0].szLineBuff[pPEreads[0].LineLen],">%s|%1.8d|%s|%d|%d|%d|%c|%d|%d\n",pszIsRand,
					NumPrevReported+NumReported,pChromSeq->szChromName,pRead->StartLoci,pRead->EndLoci+InDelSize,Pair1ReadLen,pRead->Strand ? '-' : '+',NumSubs,InDelSize);

			ReadOfs = 0;
			ReadLenRem = Pair1ReadLen;
			while(ReadLenRem)
				{
				pPEreads[0].NumCols = ReadLenRem > m_MaxFastaLineLen ? m_MaxFastaLineLen : ReadLenRem;
				CSeqTrans::MapSeq2UCAscii(&pReadSeq[ReadOfs],pPEreads[0].NumCols,&pPEreads[0].szLineBuff[pPEreads[0].LineLen]);
				pPEreads[0].LineLen += pPEreads[0].NumCols;
				pPEreads[0].LineLen += sprintf(&pPEreads[0].szLineBuff[pPEreads[0].LineLen],"\n");
				ReadLenRem -= pPEreads[0].NumCols;
				ReadOfs += pPEreads[0].NumCols;
				}

			if(bPEgen == true)
				{
				NumReported += 1;
				pReadSeq = pPEreads[1].FwdSeq; 
				if(pPairRead->Strand)
					{
					memcpy(pPEreads[1].RevSeq,pReadSeq,Pair2ReadLen+m_InDelSize);
					CSeqTrans::ReverseComplement(Pair2ReadLen+m_InDelSize,pPEreads[1].RevSeq);
					pReadSeq = pPEreads[1].RevSeq;
					}

				SimArtefacts(false,Pair2ReadLen,pReadSeq);
				SimArtefacts(true,Pair2ReadLen,pReadSeq);

				if(!m_PropRandReads || (m_PropRandReads <= RGseeds.IRandom(0,1000000)))
					{
					InDelSize = SimInDels(pPairRead,&Pair2ReadLen,pReadSeq);
					NumSubs = SimSeqErrors(Pair2ReadLen,pReadSeq);
					pszIsRand = (char *)"lcl";
					}
				else
					{
					InDelSize = 0;
					NumSubs = SimSeqRand(Pair2ReadLen,pReadSeq);
					pszIsRand = (char *)"lcr";
					}

				pPEreads[1].LineLen+=sprintf(&pPEreads[1].szLineBuff[pPEreads[1].LineLen],">%s|%1.8d|%s|%d|%d|%d|%c|%d|%d\n",pszIsRand,
						NumPrevReported+NumReported,pChromSeq->szChromName,pPairRead->StartLoci,pPairRead->EndLoci+InDelSize,Pair2ReadLen,pPairRead->Strand ? '-' : '+',NumSubs,InDelSize);

				ReadOfs = 0;
				ReadLenRem = Pair2ReadLen;
				while(ReadLenRem)
					{
					pPEreads[1].NumCols = ReadLenRem > m_MaxFastaLineLen ? m_MaxFastaLineLen : ReadLenRem;
					CSeqTrans::MapSeq2UCAscii(&pReadSeq[ReadOfs],pPEreads[1].NumCols,&pPEreads[1].szLineBuff[pPEreads[1].LineLen]);
					pPEreads[1].LineLen += pPEreads[1].NumCols;
					pPEreads[1].LineLen += sprintf(&pPEreads[1].szLineBuff[pPEreads[1].LineLen],"\n");
					ReadLenRem -= pPEreads[1].NumCols;
					ReadOfs += pPEreads[1].NumCols;
					}
				}
			break;
			}

		}

	if((pPEreads[0].LineLen + 2) > sizeof(pPEreads[0].szLineBuff) / 2)
		{
		if(write(m_hOutFile,pPEreads[0].szLineBuff,pPEreads[0].LineLen) != pPEreads[0].LineLen)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Errors whilst writing %d bytes - '%s' - %s",pPEreads[0].LineLen, m_pszOutFile, strerror(errno));
			Reset(false);
			return(eBSFerrFileAccess);
			}
		pPEreads[0].LineLen=0;
		}
	if(bPEgen == true)
		{
		if((pPEreads[1].LineLen + 2) > sizeof(pPEreads[1].szLineBuff) / 2)
			{
			if(write(m_hOutPEFile,pPEreads[1].szLineBuff,pPEreads[1].LineLen) != pPEreads[1].LineLen)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Errors whilst writing %d bytes - '%s' - %s",pPEreads[1].LineLen, m_pszOutPEFile, strerror(errno));
				delete pPEreads;
				Reset(false);
				return(eBSFerrFileAccess);
				}
			pPEreads[1].LineLen=0;
			}
		}

	pRead->Status = 1;						// mark this read as having being reported
	if(bPEgen)
		pPairRead->Status = 1;

	if(NumPrevReported + NumReported == MaxReads)
		break;
	}

if(pPEreads[0].LineLen && write(m_hOutFile,pPEreads[0].szLineBuff,pPEreads[0].LineLen) != pPEreads[0].LineLen)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Errors whilst writing %d bytes - '%s' - %s",pPEreads[0].LineLen, m_pszOutFile, strerror(errno));
	delete pPEreads;
	Reset(false);
	return(eBSFerrFileAccess);
	}
#ifdef _WIN32
_commit(m_hOutFile);
#else
fsync(m_hOutFile);
#endif

if(bPEgen == true)
	{
	if(pPEreads[1].LineLen && write(m_hOutPEFile,pPEreads[1].szLineBuff,pPEreads[1].LineLen) != pPEreads[1].LineLen)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Errors whilst writing %d bytes - '%s' - %s",pPEreads[1].LineLen, m_pszOutPEFile, strerror(errno));
		delete pPEreads;
		Reset(false);
		return(eBSFerrFileAccess);
		}
	#ifdef _WIN32
	_commit(m_hOutPEFile);
	#else
	fsync(m_hOutPEFile);
	#endif
	}
if(pPEreads != NULL)
	delete pPEreads;
return(NumReported);
}



#ifdef _WIN32
unsigned __stdcall WorkerThread(void * pThreadPars)
#else
void * WorkerThread(void * pThreadPars)
#endif
{
int Rslt = 0;
tsSRWorkerPars *pPars = (tsSRWorkerPars *)pThreadPars; // makes it easier not having to deal with casts!
CSimReads *pThis = (CSimReads *)pPars->pThis;
Rslt = pThis->ThreadSimReads(pPars);
pPars->Rslt = Rslt;
#ifdef _WIN32
_endthreadex(0);
return(eBSFSuccess);
#else
pthread_exit(NULL);
#endif
}

int
CSimReads::GenSimReads(etSRPMode PMode,		// processing mode
		etSRSEMode SEMode,	// induced sequencer error rate mode
		bool bPEgen,		// true if paired ends are to be generated
		int PEmin,			// PE minimum fragment length
		int PEmax,			// PE maximum fragment length
		double PropRandReads, // generate completely random reads at this rate
		int DistCluster,	// cluster generated reads into windows of this median length, 0 if no clustering
		double SeqErrRate,	// dynamic sequencing errors are to be induced into generated sequences at this rate
		bool bSeqErrProfile,// true if to generate composite sequencer errors with uniform profile (default is Illumina 3' skewed)
		int SNPrate,		// simulate SNPs at this rate per million bases
		int InDelSize,		// simulated InDel size range
		double InDelRate,	// simulated InDel rate per read
		etSRFMode FMode,		// output format
		int NumThreads,		// number of worker threads to use
		char Strand,		// generate for this strand '+' or '-' or for both '*'
		int NumReads,		// number of reads required (will be 2x this number if generating paired ends)
		int ReadLen,		// read lengths
		double Artef5Rate,			// rate (0..1) at which to insert 5' artefact sequences
		int NumArtef5Seqs,			// number of user specified 5' artefact sequences
		char *pszArtef5Seqs[], // 5' artefact sequences
		double Artef3Rate,			// rate (0..1) at which to insert 3' artefact sequences
		int NumArtef3Seqs,			// number of user specified 3' artefact sequences
		char *pszArtef3Seqs[], // 5' artefact sequences
		int CutMin,			// min cut length
		int CutMax,			// max cut length
		bool bDedupe,		// true if unique read sequences only to be generated
		int Region,			// Process regions 0:ALL,1:CDS,2:5'UTR,3:3'UTR,4:Introns,5:5'US,6:3'DS,7:Intergenic (default = ALL)
		int UpDnRegLen,		// if processing regions then up/down stream regulatory length
		char *pszFeatFile,	// optionally generate transcriptome reads from features or genes in this BED file
		char *pszInFile,	// input from this raw multifasta or bioseq assembly
		char *pszProfFile,	// input from this profile site preferences file
		char *pszOutPEFile, // output partner paired end simulated reads to this file
		char *pszOutFile,	// output simulated reads to this file
		char *pszOutSNPs)   // output simulated SNP loci to this file
{
int Rslt;
int ReadsOfs;
int NumReadsReq;
int ReadsCnt;
tsSRWorkerPars WorkerThreads[cMaxWorkerThreads];
tsSRWorkerPars *pCurThread;
int ReadsPerThread;
int ThreadIdx;
bool bFirst;
int MaxReadsPerBatch;		// process at most this number of simulated reads per batch
int ReportedReads;			// number of reads reported on by last call to ReportReads()
int TotReportedReads;		// total number of reads reported on by all batches processed by ReportReads()
int ExhustedChroms;			// number of threads which exhusted attempts to find chroms from which reads can be simulated
int CurNumGenReads;
int PrevNumGenReads;
int MinChromLen;

Reset(false);

RGseeds.RandomInit((int)time(NULL));

m_PMode = PMode;
m_FMode = FMode;
m_SEMode = SEMode;
m_PropRandReads = (int)(PropRandReads * 1000000.0);
m_DynProbErr = SeqErrRate;
m_bUniformDist = bSeqErrProfile;
m_InDelSize = InDelSize;
m_InDelRate = InDelRate;
m_DistCluster = DistCluster;
m_TotReqReads = bPEgen ? NumReads * 2 : NumReads;

m_Artef5Rate = Artef5Rate;				// rate (0..1) at which to insert 5' artefact sequences
m_NumArtef5Seqs = NumArtef5Seqs;		// number of user specified 5' artefact sequences
m_ppszArtef5Seqs = pszArtef5Seqs;		// 5' artefact sequences
m_Artef3Rate = Artef3Rate;				// rate (0..1) at which to insert 3' artefact sequences
m_NumArtef3Seqs = NumArtef3Seqs;		// number of user specified 3' artefact sequences
m_ppszArtef3Seqs = pszArtef3Seqs;		// 5' artefact sequences
m_Artef5SeqLens[0] = 0;
m_Artef3SeqLens[0] = 0;
for(int Idx = 0; Idx < m_NumArtef5Seqs; Idx++)
	{
	m_Artef5SeqLens[Idx] = (int)strlen(pszArtef5Seqs[Idx]);
	CSeqTrans::MapAscii2Sense(pszArtef5Seqs[Idx],m_Artef5SeqLens[Idx],m_Artef5Seqs[Idx]);
	}
for(int Idx = 0; Idx < m_NumArtef3Seqs; Idx++)
	{
	m_Artef3SeqLens[Idx] = (int)strlen(pszArtef3Seqs[Idx]);
	CSeqTrans::MapAscii2Sense(pszArtef3Seqs[Idx],m_Artef3SeqLens[Idx],m_Artef3Seqs[Idx]);
	}

if(PMode != etSRPMode::eSRPMStandard)
	{
	if((Rslt=InitProfSitePrefs(pszProfFile))!=eBSFSuccess)
		{
		Reset(false);
		return(Rslt);
		}
	}

// estimate minimum length chroms from which reads can be reliably sampled
if(bPEgen)
	MinChromLen = PEmax;
else
	MinChromLen = CutMax;
MinChromLen += 20;
if(MinChromLen < cSRMinChromLen)
	MinChromLen = cSRMinChromLen;

if(Region != eSRMEGRAny || pszFeatFile == NULL || pszFeatFile[0] == '\0')
	{
	if((Rslt=LoadGenome(MinChromLen,pszInFile))!=eBSFSuccess)	// need to load complete assembly into memory
		{
		Reset(false);
		return(Rslt);
		}
	}
else
	{
	if((Rslt=LoadTranscriptome(pszInFile,pszFeatFile))!=eBSFSuccess)	// need to load transcriptome into memory
		{
		Reset(false);
		return(Rslt);
		}
	}

// if filtering by region...
if(Region != eSRMEGRAny)
	{
	// open feature file
	if((m_pBEDFile = new CBEDfile) == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBEDfile object");
		return(eBSFerrObj);
		}

	if((Rslt = m_pBEDFile->Open(pszFeatFile))!=eBSFSuccess)
		{
		while(m_pBEDFile->NumErrMsgs())
			gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pBEDFile->GetErrMsg());
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open feature/gene BED file '%s'",pszFeatFile);
		Reset(false);
		return(Rslt);
		}

	 if(!m_pBEDFile->ContainsGeneDetail())
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"BED file '%s' does not contain gene feature regions",pszFeatFile);
		Reset(false);
		return(eBSFerrEntry);
		}
	}

if(SNPrate > 0)				// simulate SNPs?
	{
	if((Rslt = SimulateSNPs(pszOutSNPs,SNPrate))!=eBSFSuccess)
		{
		Reset(false);
		return(Rslt);
		}
	}

#ifdef _WIN32
m_hOutFile = open(pszOutFile,O_CREATETRUNC );
#else
if((m_hOutFile = open(pszOutFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
    if(ftruncate(m_hOutFile,0)!=0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",pszOutFile,strerror(errno));
			Reset(false);
			return(eBSFerrCreateFile);
			}
#endif

if(m_hOutFile < 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output file '%s'",pszOutFile);
	Reset(false);
	return(eBSFerrCreateFile);
	}
m_pszOutFile = pszOutFile;

if(bPEgen)
	{
#ifdef _WIN32
	m_hOutPEFile = open(pszOutPEFile,O_CREATETRUNC );
#else
	if((m_hOutPEFile = open(pszOutPEFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
		if(ftruncate(m_hOutPEFile,0)!=0)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",pszOutPEFile,strerror(errno));
				Reset(false);
				return(eBSFerrCreateFile);
				}
#endif

	if(m_hOutPEFile < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output file '%s'",pszOutPEFile);
		Reset(false);
		return(eBSFerrCreateFile);
		}
	m_pszOutPEFile = pszOutPEFile;
	}
else
	{
	m_pszOutPEFile = NULL;
	m_hOutPEFile = -1;
	}

// Don't bother checkpointing (write to file in batches every N reads generated) unless profiling
// Checkpointing is only used if the processing is expected to take some time to complete
MaxReadsPerBatch = cMaxProfileBatchSize;
MaxReadsPerBatch *= NumThreads;


// Allocate to hold all reads if bDedupe is TRUE even though they will be checkpointed, and written to disk, every cChkNumReads. This is because
// when deduping the reads the deduping needs to be over all reads and not just the reads in the current checkpointed batch
if(bDedupe)			// if will be deduping then generate some extra reads assuming that a few will be dups and thus be removed
	{
	NumReadsReq = (int)(((INT64)m_TotReqReads * 120)/100);
	m_NumReadsAllocd = NumReadsReq + 100;
	}
else
	{
	NumReadsReq = m_TotReqReads;
	m_NumReadsAllocd = MaxReadsPerBatch + 100;
	}

m_AllocdMemReads = (INT64)m_NumReadsAllocd * sizeof(tsSRSimRead);
#ifdef _WIN32
m_pSimReads = (tsSRSimRead *) malloc((size_t)m_AllocdMemReads);
if(m_pSimReads == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Memory allocation of %lld bytes for simulated reads failed",(INT64)m_AllocdMemReads);
	Reset(false);
	return(eBSFerrMem);
	}
#else
// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
m_pSimReads = (tsSRSimRead *)mmap(NULL,(size_t)m_AllocdMemReads, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
if(m_pSimReads == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Memory allocation of %lld bytes for simulated reads through mmap()  failed",(INT64)m_AllocdMemReads,strerror(errno));
	m_pSimReads = NULL;
	Reset(false);
	return(eBSFerrMem);
	}
#endif

memset(m_pSimReads,0,(size_t)m_AllocdMemReads);

#ifdef _WIN32
if((m_hMtxDedupe = CreateMutex(NULL,false,NULL))==NULL)
#else
if(pthread_mutex_init (&m_hMtxDedupe,NULL)!=0)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to create mutex");
	Reset(false);
	return(eBSFerrInternal);
	}

CurNumGenReads = 0;
PrevNumGenReads = 0;
ReadsOfs = 0;
TotReportedReads = 0;
bFirst =true;
TRandomCombined<CRandomMother,CRandomMersenne> RGseeds((int)time(0));
do {
	// initialise all worker thread parameters and start the threads
	if(!bDedupe)
		{
		ReadsOfs = 0;
		ReadsCnt = min((m_TotReqReads - TotReportedReads),MaxReadsPerBatch);
		}
	else
		{
		ReadsCnt = min((int)(((INT64)(m_TotReqReads - TotReportedReads) * 105)/100),MaxReadsPerBatch);
		if(ReadsCnt < 10000 && (ReadsOfs + 10000) < m_NumReadsAllocd)
			ReadsCnt = 10000;
		}
	NumReadsReq = 0;
	ExhustedChroms = 0;
	for(ThreadIdx = 0; ThreadIdx < NumThreads; ThreadIdx++)
		{
		pCurThread = &WorkerThreads[ThreadIdx];
		ReadsPerThread = ReadsCnt/(NumThreads-ThreadIdx);
		ReadsCnt -= ReadsPerThread;
		memset(pCurThread,0,sizeof(tsSRWorkerPars));
		pCurThread->pThis = this;
		pCurThread->RandSeed = RGseeds.IRandom(1,INT_MAX);
		pCurThread->ThreadIdx = ThreadIdx;
		pCurThread->bDedupe = bDedupe;
		pCurThread->ReadLen = ReadLen;
		pCurThread->Region = Region;
		pCurThread->UpDnRegLen = UpDnRegLen;
		pCurThread->CutMax = CutMax;
		pCurThread->CutMin = CutMin;

		pCurThread->bPEgen = bPEgen;
		pCurThread->PEmin = PEmin;
		pCurThread->PEmax = PEmax;
		pCurThread->NumGenReads=0;
		pCurThread->NumReqReads=ReadsPerThread;
		NumReadsReq += ReadsPerThread;
		pCurThread->pReads = &m_pSimReads[ReadsOfs];
		ReadsOfs += ReadsPerThread;
		pCurThread->PMode=PMode;
		pCurThread->Strand=Strand;
		pCurThread->bMaxIters = false;

	#ifdef _WIN32
		pCurThread->threadHandle = (HANDLE)_beginthreadex(NULL,0x0fffff,WorkerThread,pCurThread,0,&pCurThread->threadID);
	#else
		pCurThread->threadRslt =	pthread_create(&pCurThread->threadID , NULL , WorkerThread , pCurThread );
	#endif
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Randomly selecting up to %s%d reads...",bFirst?" ":" another ",NumReadsReq);
	bFirst = false;


	// wait for all threads to terminate - be patient, could be a long, long wait if dynamic Hamming dist determinations
    PrevNumGenReads = 0;
	for(ThreadIdx = 0; ThreadIdx < NumThreads; ThreadIdx++)
		{
		pCurThread = &WorkerThreads[ThreadIdx];

		// report on number of reads generated every 60 secs assuming that at least 5 reads were generated in that time period
#ifdef _WIN32
		while(WAIT_TIMEOUT == WaitForSingleObject( pCurThread->threadHandle, 60000 * 10))
			{
			WaitForSingleObject(m_hMtxDedupe,INFINITE);
			CurNumGenReads = m_CurNumGenReads;
			ReleaseMutex(m_hMtxDedupe);
			if(CurNumGenReads > (PrevNumGenReads+1))
				{
				gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: %u reads generated",CurNumGenReads);
				PrevNumGenReads = CurNumGenReads;
				}
			}
		CloseHandle(pCurThread->threadHandle);
#else
		struct timespec ts;
		int JoinRlt;
		clock_gettime(CLOCK_REALTIME, &ts);
		ts.tv_sec += 60 * 10;
		while((JoinRlt = pthread_timedjoin_np(pCurThread->threadID, NULL, &ts)) != 0)
			{
			pthread_mutex_lock(&m_hMtxDedupe);
			CurNumGenReads = m_CurNumGenReads;
			pthread_mutex_unlock(&m_hMtxDedupe);
			if(CurNumGenReads > (PrevNumGenReads+1))
				{
				gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: %u reads generated",CurNumGenReads);
				PrevNumGenReads = CurNumGenReads;
				}
			ts.tv_sec += 60;
			}

#endif
		if(pCurThread->bMaxIters)
			ExhustedChroms += 1;
		}

	if(ExhustedChroms < NumThreads)
		{
		ReportedReads = ReportReads(bPEgen,		// true if paired end simulated reads being simulated
		        Region,						// Hamming regional processing?
				ReadLen,					// read length
			    TotReportedReads,			// number of reads thus far reported on
				m_TotReqReads,				// maximum number of reads required to be reported on
				bDedupe,					// true if reads are to be deduped
				ReadsOfs);					// number of reads to report on from this batch
		TotReportedReads += ReportedReads;
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Current block processing: %1.7d total: %1.8d",ReportedReads,TotReportedReads);
		}
	else
		{
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"Unable to generate simulated reads as assembly chromosomes not of sufficent length from which to randomly select reads");
		break;
		}
	}
while(TotReportedReads < m_TotReqReads && ReadsOfs < m_NumReadsAllocd);

#ifdef _WIN32
CloseHandle(m_hMtxDedupe);
#else
pthread_mutex_destroy(&m_hMtxDedupe);
#endif
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Total simulated reads generated: %1.8d",TotReportedReads);

if(m_hOutFile != -1)
	{
#ifdef _WIN32
	_commit(m_hOutFile);
#else
	fsync(m_hOutFile);
#endif
	close(m_hOutFile);
	m_hOutFile = -1;
	}
if(m_hOutPEFile != -1)
	{
#ifdef _WIN32
	_commit(m_hOutPEFile);
#else
	fsync(m_hOutPEFile);
#endif
	close(m_hOutPEFile);
	m_hOutPEFile = -1;
	}

if(Region != eSRMEGRAny)
	{
	if(m_pBEDFile != NULL)
		{
		delete m_pBEDFile;
		m_pBEDFile = NULL;
		}
	}

Reset(true);
return(eBSFSuccess);
}



int 
CSimReads::ThreadSimReads(void * pThreadPars)
{
tsSRWorkerPars *pWorkerPars = (tsSRWorkerPars *)pThreadPars;
etSeqBase ReadSeq[cMaxReadLen+1];
etSeqBase *pReadSeq;
etSeqBase Nuc;
tsSRSimRead *pRead;

int RandChrom;

int BEDChromID;
int Features;

int RandCutSite1;
int RandCutLen;
int RandCutSite2;
int RandStrand;				// 0 if '+', 1 if '-'

int PEFragSize;
int PERandCutSite1;
int PERandCutLen;
int PERandCutSite2;
int PERandStrand;

int OctIdx;
int CurNumGenReads;

double RandCutProb;
double ProfScore;
tsSRChromSeq *pChromSeq;
etSeqBase *pSeq;
etSeqBase *pSeqSite;

int TargPsn;
int IdxLo;
int IdxHi;
UINT64 NumChromIters;

CurNumGenReads = 0;
NumChromIters = 0;
pRead = pWorkerPars->pReads;
pRead += pWorkerPars->NumGenReads;
pWorkerPars->bMaxIters = false;
TRandomCombined<CRandomMother,CRandomMersenne> RG(pWorkerPars->RandSeed);
while(pWorkerPars->NumGenReads < pWorkerPars->NumReqReads)
	{
	if(NumChromIters++ > ((UINT64)m_NumChromSeqs * 50))
		{
		pWorkerPars->bMaxIters = true;			// flag this thread is terminating because it has exhusted attempts to find chrom from which read can be simulated
		break;
		}
		// first randomly choose chromosome
	RandChrom = (int)RG.IRandom(1,m_GenomeScaledLen);

	IdxLo = 0;
	IdxHi = m_NumChromSeqs-1;
	do {
		TargPsn = (IdxLo + IdxHi) / 2;
		pChromSeq = &m_pChromSeqs[TargPsn];
		if(pChromSeq->ScaledStart <= RandChrom && (pChromSeq->ScaledStart + pChromSeq->ScaledLen) >= RandChrom)
			break;

		if(pChromSeq->ScaledStart > RandChrom)
			IdxHi = TargPsn - 1;
		else
			IdxLo = TargPsn+1;
		}
	while(IdxHi >= IdxLo);

	if(pWorkerPars->Strand != '*' && pChromSeq->Strand != '*' && pWorkerPars->Strand != pChromSeq->Strand)
		continue;

	// randomly choose cut length?
	if(pWorkerPars->CutMin != pWorkerPars->CutMax)
		{
		RandCutLen = (int)RG.IRandom(pWorkerPars->CutMin,pWorkerPars->CutMax);
		if(pWorkerPars->bPEgen)
			PERandCutLen = (int)RG.IRandom(pWorkerPars->CutMin,pWorkerPars->CutMax);
		else
			PERandCutLen = RandCutLen;
		}
	else
		{
		RandCutLen = pWorkerPars->CutMin;
		PERandCutLen = RandCutLen;
		}

		// if paired end then choose the fragment length
	if(pWorkerPars->bPEgen)
		{
		PEFragSize = (int)RG.IRandom(pWorkerPars->PEmin,pWorkerPars->PEmax);
		if(PEFragSize < (min(RandCutLen,PERandCutLen) + 1))		// allow for user simulating overlapped paired reads as required for AllpathsLG
			continue;
		}
	else
		PEFragSize = 0;

	// skip any extremely short chromosomes - could be a contig?
	if(pWorkerPars->bPEgen)
		{
		if(pChromSeq->Len < (PEFragSize + 20))
			continue;
		}
	else
		if(pChromSeq->Len < (pWorkerPars->CutMax + 20))
			continue;

	if(pChromSeq->RelDensity < 1)				// don't bother to generate reads if relative density too low
		continue;

	if(pChromSeq->RelDensity < 1000)				// if less than 1000 then don't accept this putative read if RelDensity < rand(0,999)
		{
		if(pChromSeq->RelDensity < (int)RG.IRandom(0,999))
			continue;
		}

	// try and replicate, very crude attempt!, both the variance in transcript levels and the clumped distribution of reads within transcripts
	// as observed with RNASeq sequenced read distributions
	// firstly the transcript level is simply a linear function of the first 8bp at the start of the transcript sequence
	if(m_DistCluster > 0 && pChromSeq->RelDensity == 1000)
		{
		int TransLev = 0;
		TransLev = GenPSeqIdx(8,&m_pGenomeSeq[pChromSeq->SeqOfs]);
		// overall transcript level now determined
		if(RG.IRandom(0,0x010000) >= TransLev)
			continue;
		}

	// randomly choose strand?
	if(pWorkerPars->Strand == '*')
		RandStrand = (int)RG.IRandom(0,1);
	else
		RandStrand = pWorkerPars->Strand == '-' ? 1 : 0;
	if(pWorkerPars->bPEgen)
		PERandStrand = RandStrand == 0 ? 1 : 0;

	// randomly choose initial cut sites
	// note that specified range is such that '+' and '-' start/end loci after allowing for cut lengths will always be on the chromosome
	if(pWorkerPars->bPEgen)
		{
		if(RandStrand == 0)
			{
			RandCutSite1 = (int)RG.IRandom(7,pChromSeq->Len - (PEFragSize + 7));
			PERandCutSite1 = RandCutSite1 + PEFragSize - PERandCutLen;
			}
		else
			{
			RandCutSite1 = (int)RG.IRandom(7 + PEFragSize - RandCutLen,pChromSeq->Len - RandCutLen - 7);
			PERandCutSite1 = RandCutSite1 - (PEFragSize - RandCutLen) - 1;
			}
		}
	else
		{
		RandCutSite1 = (int)RG.IRandom(7,pChromSeq->Len - (RandCutLen + 7));
		PERandCutSite1 = RandCutSite1;
		}

	RandCutSite2 = RandCutSite1 + RandCutLen;
	PERandCutSite2 = PERandCutSite1 + PERandCutLen;


	// filter by region here
	if(pWorkerPars->Region != eSRMEGRAny)
		{
		if(m_pBEDFile != NULL)
			{
			if((BEDChromID = m_pBEDFile->LocateChromIDbyName(pChromSeq->szChromName)) < 1)
				continue;
			Features = m_pBEDFile->GetFeatureBits(BEDChromID,RandCutSite1+(RandCutLen/2),RandCutSite1+(RandCutLen/2),cRegionFeatBits,pWorkerPars->UpDnRegLen);
			switch(pWorkerPars->Region) {
				case eSRMEGRCDS:			// part of feature overlaps CDS
					if(!(Features & cFeatBitCDS))
						continue;
					break;
				case eSRMEGR5UTR:			// part of feature overlaps 5'UTR
					if(!(Features & cFeatBit5UTR) || (Features & cFeatBitCDS))
						continue;
					break;
				case eSRMEGR3UTR:			// part of feature overlaps 3'UTR
					if(!(Features & cFeatBit3UTR) || (Features & (cFeatBitCDS | cFeatBit5UTR)))
						continue;
					break;
				case eSRMEGRIntrons:		// part of feature overlaps Intron
					if(!(Features & cFeatBitIntrons) || (Features & (cFeatBitCDS | cFeatBit5UTR | cFeatBit3UTR)))
						continue;
					break;
				case eSRMEGRUpstream:		// part of feature overlaps 5'upstream regulatory
					if(!(Features & cFeatBitUpstream) || (Features & (cFeatBitCDS | cFeatBit5UTR | cFeatBit3UTR | cFeatBitIntrons)))
						continue;
					break;
				case eSRMEGRDnstream:		// part of feature overlaps 3'downstream regulatory
					if(!(Features & cFeatBitDnstream) || (Features & (cFeatBitCDS | cFeatBit5UTR | cFeatBit3UTR | cFeatBitIntrons | cFeatBitUpstream)))
						continue;
					break;
				case eSRMEGRIntergenic:	// part of feature overlaps intergenic
					if(Features & (cFeatBitCDS | cFeatBit5UTR | cFeatBit3UTR | cFeatBitIntrons | cFeatBitUpstream | cFeatBitDnstream))
						continue;
					break;
				}
			}
		}

		// now for the clumping factor
		// clumps are distributed along the length of the transcript into non-overlapping clustering bins of width m_DistCluster
	if(m_DistCluster && pChromSeq->Len > m_DistCluster)
		{
		int ClumpProb;
		int ClumpOfs;
		ClumpOfs = (RandCutSite1 / m_DistCluster) * m_DistCluster; // ClumpOfs is the start of bin, of m_DistCluster len, which would contain RandCutSite
		if(((ClumpOfs/m_DistCluster) % 5) == 2)					   // every 5th bin starting at the second is the bin containing the clustered reads
			{
			pSeq = &m_pGenomeSeq[pChromSeq->SeqOfs  + ClumpOfs];	// bases at the start of clustering bin used to determine the prob of reads in that bin
			ClumpProb = GenPSeqIdx(5,pSeq);
			}
		else
			ClumpProb = 0;											// reads outside of clustering bins are given a no chance...
		if(RG.IRandom(1,1024) > ClumpProb)
				continue;
		}

	if(pWorkerPars->PMode == etSRPMode::eSRPMProfRand || pWorkerPars->PMode == etSRPMode::eSRPMProfProf)
		{
		// randomly choose site1 cut prob
		RandCutProb = RG.Random();
		pSeq =  &m_pGenomeSeq[pChromSeq->SeqOfs];
		pSeqSite = &pSeq[RandCutSite1-4];
		if(!RandStrand)
			OctIdx = GenPSeqIdx(8,pSeqSite);
		else
			OctIdx = GenMSeqIdx(8,pSeqSite);
		if(OctIdx < 0)		// -1 if pSeqSite contained 'n'
			continue;
		ProfScore = m_pProfSel[OctIdx];
		if(ProfScore < RandCutProb)
			continue;
		}

	pSeq = &m_pGenomeSeq[pChromSeq->SeqOfs];

	if(pWorkerPars->PMode == etSRPMode::eSRPMProfProf || pWorkerPars->PMode == etSRPMode::eSRPMRandProf)
		{
		// randomly choose site2 cut prob
		RandCutProb = RG.Random();
		pSeqSite = &pSeq[RandCutSite2-4];
		if(!RandStrand)
			OctIdx = GenPSeqIdx(8,pSeqSite);
		else
			OctIdx = GenMSeqIdx(8,pSeqSite);
		if(OctIdx < 0)		// -1 if seq contained 'n'
			continue;
		ProfScore = m_pProfSel[OctIdx];
		if(ProfScore < RandCutProb)
			continue;
		}

		// ensure that this sequence will be exactly alignable - slough if any contained 'n's
	pSeqSite = &pSeq[RandCutSite1];
	pReadSeq = ReadSeq;
	for(OctIdx = RandCutSite1; OctIdx < RandCutSite2; OctIdx++,pSeqSite++)
		{
		if((Nuc = (*pSeqSite & NUCONLYMSK)) > eBaseT)
			break;
		*pReadSeq++ = Nuc;
		}
	if(OctIdx < RandCutSite2)
		continue;
	*pReadSeq = eBaseEOS;

	// if generating paired ends then check if partner read contains any 'n's - if so then slough
	if(pWorkerPars->bPEgen)
		{
		pSeqSite = &pSeq[PERandCutSite1];
		for(OctIdx = PERandCutSite1; OctIdx < PERandCutSite2; OctIdx++,pSeqSite++)
			{
			if((Nuc = (*pSeqSite & NUCONLYMSK)) > eBaseT)
				break;
			}
		if(OctIdx < PERandCutSite2)
			continue;
		}

	// we have a sequence starting at RandCutSite1 and ending at RandCutSite2-1 which is of length RandCutLen
	// if non-duplicates required then mark subsequence as selected
	if(pWorkerPars->bDedupe)
		{
		pSeqSite = &pSeq[RandCutSite1];
#ifdef _WIN32
		WaitForSingleObject(m_hMtxDedupe,INFINITE);
#else
		pthread_mutex_lock(&m_hMtxDedupe);
#endif
		if(*pSeqSite & SSSELECTED)	// check if this subsequence already selected
			{
#ifdef _WIN32
			ReleaseMutex(m_hMtxDedupe);
#else
			pthread_mutex_unlock(&m_hMtxDedupe);
#endif
			continue;
			}
		*pSeqSite |= SSSELECTED;
#ifdef _WIN32
		ReleaseMutex(m_hMtxDedupe);
#else
		pthread_mutex_unlock(&m_hMtxDedupe);
#endif
		}

	pRead->Status = 0;
	pRead->ChromSeqID = pChromSeq->ChromSeqID;
	pRead->ChromID = pChromSeq->ChromID;
	pRead->Strand = RandStrand;
	pRead->StartLoci = RandCutSite1;
	pRead->EndLoci = RandCutSite2 - 1;
	pRead->Len = RandCutLen;
	pRead->pSeq = &pSeq[RandCutSite1];
	pWorkerPars->NumGenReads += 1;
	pRead->FlgPE2 = 0;
	if(pWorkerPars->bPEgen)
		pRead->pPartner = &pRead[1];
	else
		pRead->pPartner = NULL;
	pRead += 1;

	if(pWorkerPars->bPEgen)
		{
		pRead->pPartner = &pRead[-1];
		pRead->FlgPE2 = 1;
		pRead->Status = 0;
		pRead->ChromSeqID = pChromSeq->ChromSeqID;
		pRead->ChromID = pChromSeq->ChromID;
		pRead->Strand = PERandStrand;
		pRead->StartLoci = PERandCutSite1;
		pRead->EndLoci = PERandCutSite2 - 1;
		pRead->Len = PERandCutLen;
		pRead->pSeq = &pSeq[pRead->StartLoci];
		pWorkerPars->NumGenReads += 1;
		pRead += 1;
		}

	CurNumGenReads += 1;
	NumChromIters = 0;

	// time to let main thread know that some progress is being made?
	if(CurNumGenReads >= 500)
		{
#ifdef _WIN32
		WaitForSingleObject(m_hMtxDedupe,INFINITE);
#else
		pthread_mutex_lock(&m_hMtxDedupe);
#endif

		m_CurNumGenReads += CurNumGenReads;
		CurNumGenReads = 0;
#ifdef _WIN32
		ReleaseMutex(m_hMtxDedupe);
#else
		pthread_mutex_unlock(&m_hMtxDedupe);
#endif
		}
	}
#ifdef _WIN32
_endthreadex(0);
return(eBSFSuccess);
#else
pthread_exit(NULL);
#endif
}





int
CSimReads::LocateRevCplSeq(int Len,etSeqBase *pProbe,int NumSortedReads)
{
int Rslt;
etSeqBase RevSeq[cMaxReadLen+1];
for(int Idx = 0; Idx < Len; Idx++,pProbe++)
	RevSeq[Idx] = *pProbe & NUCONLYMSK;
CSeqTrans::ReverseComplement(Len,RevSeq);
Rslt = LocateFirstExact(RevSeq,Len,0,NumSortedReads-1);
return(Rslt);
}



int			// index of exactly matching probe or -1 if no match
CSimReads::LocateFirstExact(etSeqBase *pProbe,				// pts to probe sequence
				 int ProbeLen,					// probe length to exactly match over
				  int IdxLo,					// low index in m_pSimReads
				  int IdxHi)					// high index in m_pSimReads
{
etSeqBase *pEl1;
etSeqBase *pEl2;

int CmpRslt;

int TargPsn;
int LowPsn;

pEl1 = pProbe;
do {
	TargPsn = (IdxLo + IdxHi) / 2;
	pEl2 = m_pSimReads[TargPsn].pSeq;
	CmpRslt = CmpSeqs(ProbeLen,pEl1,pEl2);

	if(!CmpRslt)	// if have a match but might not be the lowest indexed match
		{
		if(TargPsn == 0 || IdxLo == TargPsn) // check if lowest
			return(TargPsn);
		LowPsn = LocateFirstExact(pProbe,ProbeLen,IdxLo,TargPsn - 1);
		return(LowPsn < 0 ? TargPsn : LowPsn);
		}

	if(CmpRslt < 0)
		IdxHi = TargPsn - 1;
	else
		IdxLo = TargPsn+1;
	}
while(IdxHi >= IdxLo);

return(-1);	// unable to locate any instance of pProbe
}



// SortSimLoci
// Sort simulated reads by chrom,loci,len
int
CSimReads::SortSimLoci(const void *arg1, const void *arg2)
{
tsSRSimRead *pEl1 = (tsSRSimRead *)arg1;
tsSRSimRead *pEl2 = (tsSRSimRead *)arg2;
if(pEl1->ChromID < pEl2->ChromID)
	return(-1);
if(pEl1->ChromID > pEl2->ChromID)
	return(1);
if(pEl1->StartLoci < pEl2->StartLoci)
	return(-1);
if(pEl1->StartLoci > pEl2->StartLoci)
	return(1);
if(pEl1->Len < pEl2->Len)
	return(-1);
if(pEl1->Len > pEl2->Len)
	return(1);
return(0);
}


// SortSimReads
// Sort simulated reads by sequence
int
CSimReads::SortSimReads(const void *arg1, const void *arg2)
{
tsSRSimRead *pEl1 = (tsSRSimRead *)arg1;
tsSRSimRead *pEl2 = (tsSRSimRead *)arg2;
etSeqBase *pSeq1 = pEl1->pSeq;
etSeqBase *pSeq2 = pEl2->pSeq;
etSeqBase Base1;
etSeqBase Base2;
int Idx;
int Len = min(pEl1->Len,pEl2->Len);
for(Idx=0; Idx < Len; Idx++)
	if((Base1 = (*pSeq1++ & NUCONLYMSK)) != (Base2 = (*pSeq2++ & NUCONLYMSK)))
		return(Base1 - Base2);
return(0);
}

