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


CSmithWaterman::CSmithWaterman(void)
{
m_pTrcBckCells = NULL;
m_pColBands = NULL;
m_pProbe = NULL;
m_pTarg = NULL;
Reset();
}

CSmithWaterman::~CSmithWaterman(void)
{
Reset();
}

// Reset
// Resets instance to that immediately following that of instance construction
// Note that all buffers will be deallocated..
void
CSmithWaterman::Reset(void)
{
if(m_pTrcBckCells != NULL)
	{
#ifdef _WIN32
	free(m_pTrcBckCells);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pTrcBckCells != MAP_FAILED)
		munmap(m_pTrcBckCells,m_TrcBckCellsAllocdSize);
#endif
	}
m_TrcBckCellsAllocdSize = 0;
m_TrcBckCellsAllocd = 0;	
m_TrcBckCellsUsed = 0;		

if(m_pColBands != NULL)			// each column in the matrix will contain a band of active cells
	{
	delete m_pColBands;
	m_pColBands = NULL;
	}
m_ColBandsAllocd = 0;
m_ColBandsUsed = 0;

if(m_pProbe != NULL)
	{
	delete m_pProbe;
	m_pProbe = NULL;
	}

m_ProbeAllocd = 0;	// actual m_pProbe alloc'd size (in teSWSeqBases's)
m_ProbeLen = 0;			// current probe length

if(m_pTarg != NULL)				// alloc'd target sequence memory
	{
	delete m_pTarg;
	m_pTarg = NULL;
	}
m_TargAllocd = 0;			// actual m_pTarg alloc'd size (in teSWSeqBases's)
m_TargLen = 0;				// current targ length

m_MatchScore = cSWDfltMatchScore;				
m_MismatchPenalty = cSWDfltMismatchPenalty;		
m_GapOpenPenalty = cSWDfltGapOpenPenalty;			
m_GapExtnPenalty = cSWDfltGapExtnPenalty;			
m_DlyGapExtn = cSWDfltDlyGapExtn;
m_ProgPenaliseGapExtn = cSWDfltProgPenaliseGapExtn;
m_SWBandInitial = 0;
m_SWPathLenDiff = 1.0;

m_PeakScore = 0;
m_PeakProbeIdx = 0;
m_PeakTargIdx = 0;
m_NumBasesAligned = 0;
m_NumBasesExact = 0;
m_NumProbeInserts = 0;
m_NumTargInserts = 0;
m_ProbeAlignStartOfs = 0;
m_TargAlignStartOfs = 0;
m_bBanded = false;
m_bAligned = false;
}

// SetScores
// Set match, mismatch, gap opening and gap extension scores
bool 
CSmithWaterman::SetScores(int MatchScore,	// score for match
				int MismatchPenalty,		// penalty for mismatch
				int GapOpenPenalty,			// penalty for opening a gap
				int GapExtnPenalty,			// penalty if extending already opened gap
				int DlyGapExtn,				// delayed gap penalties, only apply gap extension penalty if gap at least this length
				int ProgPenaliseGapExtn)	// if non-zero then progressively increment gap extension penalty for gaps of at least this length, 0 to disable, used for PacBio style error profiles
{
if(MatchScore <= 0 || MatchScore > 100 || DlyGapExtn < 1 || DlyGapExtn > 63 || ProgPenaliseGapExtn < 0 || ProgPenaliseGapExtn > 63 || MismatchPenalty < -100 || MismatchPenalty > 0 || GapOpenPenalty < -100 || GapOpenPenalty > 0 || GapExtnPenalty < -100 ||  GapExtnPenalty > 0)
	return(false);
m_MatchScore = MatchScore;				
m_MismatchPenalty = MismatchPenalty;		
m_GapOpenPenalty = GapOpenPenalty;			
m_GapExtnPenalty = GapExtnPenalty;			
m_DlyGapExtn = DlyGapExtn;				
if(ProgPenaliseGapExtn > 0 && ProgPenaliseGapExtn < DlyGapExtn)
	ProgPenaliseGapExtn = DlyGapExtn;	
m_ProgPenaliseGapExtn = ProgPenaliseGapExtn;
m_ColBandsUsed = 0;
m_TrcBckCellsUsed = 0;
m_NumBasesAligned = 0;
m_NumBasesExact = 0;
m_NumProbeInserts = 0;
m_NumTargInserts = 0;
m_ProbeAlignStartOfs = 0;
m_TargAlignStartOfs = 0;
m_PeakProbeIdx = 0;
m_PeakTargIdx = 0;
m_bAligned = false;
return(true);
}

// SetProbe
// Set probe sequence to use in subsequent alignments
bool 
CSmithWaterman::SetProbe(uint32_t Len,etSeqBase *pSeq)
{
if(Len < cSWMinProbeOrTargLen || Len > cSWMaxProbeOrTargLen || pSeq == NULL || *pSeq > eBaseN) 	// can't be bothered with very short or very long probes!
	return(false);

if(m_pProbe == NULL || m_ProbeAllocd < Len)
	{
	if(m_pProbe != NULL)
		delete m_pProbe;
	m_ProbeAllocd = Len + 100;
	m_pProbe = new etSeqBase [m_ProbeAllocd];
	if(m_pProbe == NULL)
		return(false);
	}
memmove(m_pProbe,pSeq,Len);
m_ProbeLen = Len;
m_TrcBckCellsUsed = 0;
m_ColBandsUsed = 0;
m_NumBasesAligned = 0;
m_NumBasesExact = 0;
m_NumProbeInserts = 0;
m_NumTargInserts = 0;
m_ProbeAlignStartOfs = 0;
m_TargAlignStartOfs = 0;
m_PeakProbeIdx = 0;
m_PeakTargIdx = 0;
m_bAligned = false;
return(true);
}

// SetTarg
// Set target sequence to use in subsequent alignments
bool 
CSmithWaterman::SetTarg( uint32_t Len,etSeqBase *pSeq)
{
if(Len < cSWMinProbeOrTargLen || Len > cSWMaxProbeOrTargLen || pSeq == NULL || *pSeq > eBaseN)	// can't be bothered with very short or very long targets!
	return(false);
if(m_pTarg == NULL || m_TargAllocd < Len)
	{
	if(m_pTarg != NULL)
		delete m_pTarg;
	m_TargAllocd = Len + 100;
	m_pTarg = new etSeqBase [m_TargAllocd];
	if(m_pTarg == NULL)
		return(false);
	}
memmove(m_pTarg,pSeq,Len);
m_TargLen = Len;
m_TrcBckCellsUsed = 0;
m_ColBandsUsed = 0;
m_NumBasesAligned = 0;
m_NumBasesExact = 0;
m_NumProbeInserts = 0;
m_NumTargInserts = 0;
m_ProbeAlignStartOfs = 0;
m_TargAlignStartOfs = 0;
m_PeakProbeIdx = 0;
m_PeakTargIdx = 0;
m_bAligned = false;
return(true);
}


// Align
// Align probe SetProbe() against target SetTarg() using scores specified via SetScores() 
// using smith-waterman dynamic programing with optimisations
int						// smith-waterman style local alignment, returns highest score
CSmithWaterman::Align(bool bBanded,			// true (currently experimental) to use banded or constrained SW 
				 uint32_t MaxStartNonOverlap,	// if banded then initial non-overlapping path expected to be at most this many bp, 0 for no limits
				 double MaxPathLenDiff)		// if banded then path length differential between query and probe for max score be at most this proportion of query length
{
uint32_t IdxP;							// current index into m_Probe[]
uint32_t IdxT;							// current index into m_Targ[]
uint64_t NumCells;						// m_ProbeLen * m_TargLen - total number of cells
int NewScore;
int DiagScore;							// putative diagonal score
int LeftScore;							// putative left score
int DownScore;							// putative down score
int PrevScore;							// score in back referenced cell
tSWTrcBckCell DiagDir;					// to hold back reference direction as diagonal  

etSeqBase *pProbe;						// current &m_Probe[IdxP]
etSeqBase ProbeBase;					// current m_Probe[IdxP]

etSeqBase *pTarg;						// current &m_Targ[IdxT]
etSeqBase TargBase;						// current m_Targ[IdxT]

int LeftInDelLen;
int DownInDelLen;

tSWTrcBckCell *pCell;
tSWTrcBckCell *pPrevCell;
bool bMatch;
uint32_t StartIdxT;
uint32_t EndIdxT;
uint32_t DeltaIdxT;

m_bAligned = false;
if(m_ProbeLen < cSWMinProbeOrTargLen || m_TargLen < cSWMinProbeOrTargLen || ((uint64_t)m_ProbeLen * (uint64_t)m_TargLen > (int64_t)cSWMaxCells))
	return(eBSFerrMaxEntries);	
m_bBanded = bBanded;
if(bBanded)
	{
	if(MaxStartNonOverlap < 5)
		MaxStartNonOverlap = 5;
	m_SWBandInitial = MaxStartNonOverlap;
	if(MaxPathLenDiff < 0.05)
		MaxPathLenDiff = 0.05;
	else
		if(MaxPathLenDiff > 1.0)
			MaxPathLenDiff = 1.0;
	m_SWPathLenDiff = MaxPathLenDiff;
	}
else
	{
	m_SWBandInitial = 0;
	MaxPathLenDiff = 0.0;
	}

NumCells = ((uint64_t)m_ProbeLen * m_TargLen);
if(m_bBanded)
	NumCells /= 10;

if(m_pTrcBckCells == NULL || m_TrcBckCellsAllocd < NumCells || ((uint64_t)m_TrcBckCellsAllocd > (uint64_t)NumCells * 5))
	{
	NumCells += 100;						// small overallocation as a saftety margin
	if(m_pTrcBckCells != NULL)
		{
#ifdef _WIN32
		free(m_pTrcBckCells);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
		if(m_pTrcBckCells != MAP_FAILED)
			munmap(m_pTrcBckCells,m_TrcBckCellsAllocdSize);
#endif
		}
	m_TrcBckCellsAllocd = NumCells;
	m_TrcBckCellsAllocdSize = sizeof(tSWTrcBckCell) * m_TrcBckCellsAllocd;
#ifdef _WIN32
	m_pTrcBckCells = (tSWTrcBckCell *) malloc(m_TrcBckCellsAllocdSize);
	if(m_pTrcBckCells == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to allocate %I64d bytes contiguous memory for traceback cells",(int64_t)m_TrcBckCellsAllocdSize);
		m_TrcBckCellsAllocdSize = 0;
		m_TrcBckCellsAllocd = 0;
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and seems to have issues if more than 2GB allocation
	m_pTrcBckCells = (tSWTrcBckCell *)mmap(NULL,m_TrcBckCellsAllocdSize, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pTrcBckCells == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to allocate %I64d bytes contiguous memory for traceback cells",(int64_t)m_TrcBckCellsAllocdSize);
		m_TrcBckCellsAllocdSize = 0;
		m_TrcBckCellsAllocd = 0;
		return(eBSFerrMem);
		}
#endif

	}
m_TrcBckCellsUsed = 0;

if(m_bBanded)
	{
	if(m_pColBands == NULL || m_ColBandsAllocd < m_ProbeLen || m_ColBandsAllocd > m_ProbeLen * 2)
		{
		if(m_pColBands != NULL)
			delete m_pColBands;
		m_ColBandsAllocd = m_ProbeLen + 1000;
		m_pColBands = new tsSWColBand [m_ColBandsAllocd];
		if(m_pColBands == NULL)
			return(false);
		}
	}

// cell defaults are score = 0, no backreferences, no gap extension
memset(m_pTrcBckCells,0,NumCells * sizeof(tSWTrcBckCell));

m_PeakScore = 0;
m_NumBasesAligned = 0;
m_NumBasesExact = 0;
m_ProbeAlignStartOfs = 0;
m_TargAlignStartOfs = 0;
m_PeakProbeIdx = 0;
m_PeakTargIdx = 0;
pProbe = m_pProbe;
for(IdxP = 0; IdxP < m_ProbeLen; IdxP++)
	{
	ProbeBase = *pProbe++ & ~cRptMskFlg;
	if(m_bBanded && (m_SWBandInitial > 0 && m_SWPathLenDiff > 0.0))
		{
		DeltaIdxT = (int)(((int64_t)m_TargLen * IdxP) / ((int64_t)m_ProbeLen * (1.0/m_SWPathLenDiff)));
		if(DeltaIdxT < m_SWBandInitial)
			DeltaIdxT = m_SWBandInitial;
		if(DeltaIdxT > m_TargLen)
			DeltaIdxT = m_TargLen;

		StartIdxT = (int)(((int64_t)m_TargLen * IdxP) / m_ProbeLen);
		EndIdxT = StartIdxT;
		if(StartIdxT >= DeltaIdxT)
			StartIdxT -= DeltaIdxT;
		EndIdxT += DeltaIdxT;
		if(EndIdxT > (int)m_TargLen)
			EndIdxT = m_TargLen;
		}
	else
		{
		StartIdxT = 0;
		EndIdxT = m_TargLen;
		}

	if(!m_bBanded)
		pCell = &m_pTrcBckCells[((uint64_t)IdxP * m_TargLen) + (uint64_t)StartIdxT];
	pTarg = &m_pTarg[StartIdxT];
	for(IdxT = StartIdxT; IdxT < EndIdxT; IdxT++,pCell++)
		{
		if(m_bBanded)
			pCell = AllocBandCell(IdxP+1,IdxT+1);

		TargBase = *pTarg++ & ~cRptMskFlg;
		bMatch = ProbeBase == TargBase;
		LeftInDelLen = 0;
		DownInDelLen = 0;

		// calc the 3 scores (DiagScore, LeftScore, DownScore) so can determine which highest scoring path direction to take

		// diagonal is either MatchScore or MismatchScore added to prev diagonal score
		if(IdxT > 0 && IdxP > 0)
			{
			if(!m_bBanded)
				pPrevCell = pCell - m_TargLen - 1; 
			else
				pPrevCell = DerefBandCellDiag(IdxP+1,IdxT+1);
			if(pPrevCell != NULL)
					{
					PrevScore = (int)(*pPrevCell & cSWScoreMsk);
					DiagScore = PrevScore + (bMatch ? m_MatchScore : m_MismatchPenalty);
					DiagDir = (DiagScore > 0) ? cSWTrcBckDiagFlg : 0;
					if(DiagScore > (int)cSWScoreMsk)
						DiagScore = cSWScoreMsk;
					}
				else
					DiagScore = 0;
			}
		else
			{
			DiagScore = bMatch ? m_MatchScore : 0;
			*pCell = DiagScore | (bMatch ? cSWTrcBckMatchFlg : 0);
			continue;
			}

		// leftscore is either GapExtnPenalty (if gap already opened and gap at least m_DlyGapExtn) or GapOpenScore added to prev left score
		if(!m_bBanded)
			pPrevCell = pCell - m_TargLen;
		else
			pPrevCell = DerefBandCellLeft(IdxP+1,IdxT+1);

		int GapExtnPenalty;
		if(pPrevCell!=NULL)
			{
			PrevScore = (int)(*pPrevCell & cSWScoreMsk);
			GapExtnPenalty = 0;
			if(*pPrevCell & cSWGapOpnFlg)  // was there a gap previously opened?
				{
				LeftInDelLen =  1 + ((*pPrevCell & cSWInDelLenMsk) >> cSWInDelLenShf); // extending existing gap
				if(LeftInDelLen > 63)		// clamp as only 6bits avail to hold length
					LeftInDelLen = 63;

				if(m_GapExtnPenalty !=0  && LeftInDelLen >= m_DlyGapExtn)
					{
					if(m_ProgPenaliseGapExtn == 0 || LeftInDelLen < m_ProgPenaliseGapExtn)   
						GapExtnPenalty = m_GapExtnPenalty;	
					else
						GapExtnPenalty = m_GapExtnPenalty * 2; 
					}
				LeftScore = PrevScore + GapExtnPenalty; 
				}
			else
				{
				LeftInDelLen = 1;			// opening a gap
				if(m_DlyGapExtn == 1)
					GapExtnPenalty = m_GapExtnPenalty;
				if(m_ProgPenaliseGapExtn == 1)
					GapExtnPenalty += m_GapExtnPenalty;
				LeftScore = PrevScore + m_GapOpenPenalty + GapExtnPenalty;
				}
			LeftInDelLen <<= cSWInDelLenShf;
			}
		else
			{
			LeftInDelLen = 0;
			LeftScore = 0;
			}

		// down score is either GapExtnPenalty (if gap already opened) or GapOpenScore added to prev down score
		if(!m_bBanded)
			pPrevCell = pCell - 1;
		else
			pPrevCell = DerefBandCellDown(IdxP+1,IdxT+1);
		if(pPrevCell != NULL)
			{
			PrevScore = (int)(*pPrevCell & cSWScoreMsk);
			GapExtnPenalty = 0;
			if(*pPrevCell & cSWGapOpnFlg)  // was there a gap previously opened?
				{
				DownInDelLen =  1 + ((*pPrevCell & cSWInDelLenMsk) >> cSWInDelLenShf); // extending gap
				if(DownInDelLen > 63)			// clamp as only 6bits avail to hold length
					DownInDelLen = 63;

				if(m_GapExtnPenalty != 0 && DownInDelLen >= m_DlyGapExtn)
					{
					if(m_ProgPenaliseGapExtn == 0 || DownInDelLen < m_ProgPenaliseGapExtn)   
						GapExtnPenalty = m_GapExtnPenalty;	
					else
						GapExtnPenalty = m_GapExtnPenalty * 2;  
					}
				DownScore = PrevScore + GapExtnPenalty; 
				}
			else
				{
				DownInDelLen = 1;
				if(m_DlyGapExtn == 1)
					GapExtnPenalty = m_GapExtnPenalty;
				if(m_ProgPenaliseGapExtn == 1)
					GapExtnPenalty += m_GapExtnPenalty;
				DownScore = PrevScore + m_GapOpenPenalty + GapExtnPenalty;
				}
			DownInDelLen <<= cSWInDelLenShf;
			}
		else
			DownScore = 0;
		
		// if no score was > 0 then set cell score + flags to 0
		if(DiagScore <= 0 && DownScore <= 0 && LeftScore <= 0)
			{
			*pCell = 0;	
			continue;		
			}

		// select highest score into cell together with traceback and gap opened flag..
		if(DiagScore >= DownScore && DiagScore >= LeftScore)
			{
			*pCell = DiagScore | DiagDir | (bMatch ? cSWTrcBckMatchFlg : 0);
			NewScore = DiagScore;
			}
		else
			if(DownScore >= LeftScore)
				{			
				*pCell = DownScore | cSWTrcBckDownFlg | (bMatch ? 0 : (cSWGapOpnFlg | (uint32_t)DownInDelLen));
				NewScore = DownScore;
				}
			else
				{			
				*pCell = LeftScore | cSWTrcBckLeftFlg | (bMatch ? 0 : (cSWGapOpnFlg | (uint32_t)LeftInDelLen));
				NewScore = LeftScore;
				}

		if(NewScore > m_PeakScore)
			{
			m_PeakScore = NewScore;
			m_PeakProbeIdx = IdxP;
			m_PeakTargIdx = IdxT;
			}
		}
	}
m_bAligned = true;
return(m_PeakScore);
}



// GetNumAlignedBases
// get number of bases which align (exactly plus subs), excluding InDels
// Also internally updates m_ProbeAlignStart and m_TargAlignStart
int
CSmithWaterman::GetNumAlignedBases(void)	 // get number of bases which align (exactly plus subs), excluding InDels
{
tSWTrcBckCell TrcBckDir;
tSWTrcBckCell *pPeakCell;
uint32_t ProbeIdx;
uint32_t TargIdx;

if(!m_bAligned || m_PeakProbeIdx == 0 || m_PeakTargIdx == 0)
	return(0);

if(m_NumBasesAligned > 0)
	return(m_NumBasesAligned);

m_ProbeAlignStartOfs = m_PeakProbeIdx;
m_TargAlignStartOfs = m_PeakTargIdx;
m_NumBasesAligned = 0;
m_NumBasesExact = 0;
m_NumProbeInserts = 0;
m_NumTargInserts = 0;
if(m_bBanded)
	{
	ProbeIdx = m_PeakProbeIdx+1;
	TargIdx =  m_PeakTargIdx+1;
	pPeakCell = DerefBandCell(ProbeIdx,TargIdx);
	}
else
	pPeakCell = &m_pTrcBckCells[((uint64_t)m_PeakProbeIdx * m_TargLen) +  m_PeakTargIdx];

do {
	switch(TrcBckDir = (*pPeakCell & cSWTrcBckMsk)) {
		case cSWTrcBckDiagFlg:		// back on the diagonal, either a match or mismatch
			m_NumBasesAligned += 1;
			if(*pPeakCell & cSWTrcBckMatchFlg)
				m_NumBasesExact += 1;
			if(m_bBanded)
				{
				ProbeIdx -= 1;
				TargIdx -= 1;
				pPeakCell = DerefBandCell(ProbeIdx,TargIdx);
				}
			else
				pPeakCell -= m_TargLen + 1;
			m_ProbeAlignStartOfs -= 1;
			m_TargAlignStartOfs -= 1;
			break;


		case cSWTrcBckLeftFlg:			// left, insertion into probe or deletion from target
			if(m_bBanded)
				{
				ProbeIdx -= 1;
				pPeakCell = DerefBandCell(ProbeIdx,TargIdx);
				}
			else
				pPeakCell -= m_TargLen;		// treating as insertion into probe
			m_ProbeAlignStartOfs -= 1;
			m_NumProbeInserts += 1;
			break;

		case cSWTrcBckDownFlg:			// down, insertion into target or deletion  
			if(m_bBanded)
				{
				TargIdx -= 1;
				pPeakCell = DerefBandCell(ProbeIdx,TargIdx);
				}
			else
				pPeakCell -= 1;
			m_TargAlignStartOfs -= 1;
			m_NumTargInserts += 1;
			break;

		default:					// cell has no trace back so must be final cell in path
			m_NumBasesAligned += 1;
			if(*pPeakCell & cSWTrcBckMatchFlg)
				m_NumBasesExact += 1;
			break;		
		}
	}
while(TrcBckDir);
if(m_NumBasesAligned > 0)
	{
	m_ProbeAlignStartOfs += 1;
	m_TargAlignStartOfs += 1;
	}
return(m_NumBasesAligned);
}



int									   // offset or -1 if errors
CSmithWaterman::GetProbeStartOfs(void) // get offset (0..n) in probe at which alignment starts
{
if(!GetNumAlignedBases())
	return(-1);
return(m_ProbeAlignStartOfs);
}

int										// psn or -1 if errors
CSmithWaterman::GetTargStartOfs(void)  // get offset (0..n) in target at which alignment starts
{
if(!GetNumAlignedBases())
	return(-1);
return(m_TargAlignStartOfs);
}

int						// < 0 errors, 0 - no anchors, 1 - anchors returned (note that returned anchor offsets are 1 based inclusive)
CSmithWaterman::GetAnchors(uint32_t MinAnchorLen,	// anchors must be of at least this length
						   uint32_t *pProbeOfs5,	// returned probe 5' anchor starts at this probe sequence offset + 1
						   uint32_t *pTargOfs5,	// returned target 5' anchor starts at this target sequence offset + 1
						   uint32_t *pProbeOfs3,	// returned probe 3' anchor ends at this probe sequence offset + 1
						   uint32_t *pTargOfs3)	// returned target 3' anchor ends at this target sequence offset + 1
{
tSWTrcBckCell TrcBckDir;
tSWTrcBckCell *pPeakCell;
uint32_t ProbeIdx;
uint32_t TargIdx;
bool bHaveAnchors;
uint32_t CurAnchorLen;
uint32_t ProbeOfs5;
uint32_t TargOfs5;
uint32_t ProbeOfs3;
uint32_t TargOfs3;

if(GetNumAlignedBases() < (int)MinAnchorLen)
	{
	if(pProbeOfs5 != NULL)
		*pProbeOfs5 = 0;
	if(pTargOfs5 != NULL)
		*pTargOfs5 = 0;
	if(pProbeOfs3 != NULL)
		*pProbeOfs3 = 0;
	if(pTargOfs3 != NULL)
		*pTargOfs3 = 0;
	return(0);
	}

bHaveAnchors = false;
CurAnchorLen = 0;
ProbeOfs5 = 0;
TargOfs5 = 0;
ProbeOfs3 = 0;
TargOfs3 = 0;
ProbeIdx = m_PeakProbeIdx+1;
TargIdx =  m_PeakTargIdx+1;
if(m_bBanded)
	pPeakCell = DerefBandCell(ProbeIdx,TargIdx);
else
	pPeakCell = &m_pTrcBckCells[((uint64_t)m_PeakProbeIdx * m_TargLen) +  m_PeakTargIdx];

CurAnchorLen = 0;
do
	{
	TrcBckDir = *pPeakCell  & cSWTrcBckMsk;
	switch(TrcBckDir) {
		case cSWTrcBckDiagFlg:			// back on the diagonal, could be part of an anchor
			CurAnchorLen += 1;
			if(CurAnchorLen >= MinAnchorLen)
				{
				ProbeOfs5 = ProbeIdx;
				TargOfs5 = TargIdx;

				if(!bHaveAnchors)		// first anchor discovered must be the 3' anchor
					{
					ProbeOfs3 = ProbeIdx + CurAnchorLen - 1;
					TargOfs3 = TargIdx + CurAnchorLen - 1;
					bHaveAnchors = true;
					}
				}
				
			ProbeIdx -= 1;
			TargIdx -= 1;
			if(m_bBanded)
				pPeakCell = DerefBandCell(ProbeIdx,TargIdx);
			else
				pPeakCell -= m_TargLen + 1;
			break;

		case cSWTrcBckLeftFlg:			// left, treating as insertion into probe so not part of an anchor
			CurAnchorLen = 0; 
			ProbeIdx -= 1;
			if(m_bBanded)
				pPeakCell = DerefBandCell(ProbeIdx,TargIdx);
			else
				pPeakCell -= m_TargLen;
			break;

		case cSWTrcBckDownFlg:		     // down, treating as insertion into target so not part of an anchor
			CurAnchorLen = 0; 
			TargIdx -= 1;
			if(m_bBanded)
				pPeakCell = DerefBandCell(ProbeIdx,TargIdx);
			else
				pPeakCell -= 1;
			break;

		default:					// no direction, must be final cell so treat as if part of an anchor
			CurAnchorLen += 1; 
			if(CurAnchorLen >= MinAnchorLen)
				{
				ProbeOfs5 = ProbeIdx;
				TargOfs5 = TargIdx;

				if(!bHaveAnchors)		// first anchor discovered must be the 3' anchor
					{
					ProbeOfs3 = ProbeIdx + CurAnchorLen - 1;
					TargOfs3 = TargIdx + CurAnchorLen - 1;
					bHaveAnchors = true;
					}
				}
			break;	
		}
	}
while(TrcBckDir);
if(bHaveAnchors)
	{
	if(pProbeOfs5 != NULL)
		*pProbeOfs5 = ProbeOfs5;
	if(pTargOfs5 != NULL)
		*pTargOfs5 = TargOfs5;
	if(pProbeOfs3 != NULL)
		*pProbeOfs3 = ProbeOfs3;
	if(pTargOfs3 != NULL)
		*pTargOfs3 = TargOfs3;
	}

return(bHaveAnchors ? 1 : 0);
}


int 
CSmithWaterman::GetProbeAlign(uint32_t Len, etSeqBase *pBuff) // get probe alignment
{
tSWTrcBckCell TrcBckDir;
tSWTrcBckCell *pPeakCell;
etSeqBase *pProbe;
etSeqBase *pProbeStart;
uint32_t ProbeIdx;
uint32_t TargIdx;

if(!GetNumAlignedBases())
	return(-1);
if(Len < (m_NumBasesAligned + m_NumProbeInserts + m_NumTargInserts))
	return(-1);

if(m_bBanded)
	{
	ProbeIdx = m_PeakProbeIdx+1;
	TargIdx =  m_PeakTargIdx+1;
	pPeakCell = DerefBandCell(ProbeIdx,TargIdx);
	}
else
	pPeakCell = &m_pTrcBckCells[((uint64_t)m_PeakProbeIdx * m_TargLen) +  m_PeakTargIdx];
pProbeStart = m_pProbe + m_ProbeAlignStartOfs;
pProbe = pProbeStart + m_NumBasesAligned + m_NumProbeInserts - 1;
pBuff += m_NumBasesAligned + m_NumProbeInserts + m_NumTargInserts - 1;

do
	{
	TrcBckDir = *pPeakCell  & cSWTrcBckMsk;
	switch(TrcBckDir) {
		case cSWTrcBckDiagFlg:			// back on the diagonal, report base
			*pBuff-- = *pProbe--; 
			if(m_bBanded)
				{
				ProbeIdx -= 1;
				TargIdx -= 1;
				pPeakCell = DerefBandCell(ProbeIdx,TargIdx);
				}
			else
				pPeakCell -= m_TargLen + 1;
			break;

		case cSWTrcBckLeftFlg:			// left, treating as insertion into probe so report base
			*pBuff-- = *pProbe--; 
			if(m_bBanded)
				{
				ProbeIdx -= 1;
				pPeakCell = DerefBandCell(ProbeIdx,TargIdx);
				}
			else
				pPeakCell -= m_TargLen;
			break;

		case cSWTrcBckDownFlg:		     // down, treating as insertion into target so report as InDel
			*pBuff-- = eBaseInDel; 
			if(m_bBanded)
				{
				TargIdx -= 1;
				pPeakCell = DerefBandCell(ProbeIdx,TargIdx);
				}
			else
				pPeakCell -= 1;
			break;

		default:					// no direction, must be final cell so report probe base
			*pBuff = *pProbe; 
			break;	
		}
	}
while(TrcBckDir);
return(m_NumBasesAligned + m_NumProbeInserts + m_NumTargInserts);
}

int 
CSmithWaterman::GetTargAlign(uint32_t Len, etSeqBase *pBuff) // get target alignment
{
tSWTrcBckCell TrcBckDir;
tSWTrcBckCell *pPeakCell;
etSeqBase *pTarg;
etSeqBase *pTargStart;
uint32_t ProbeIdx;
uint32_t TargIdx;

if(!GetNumAlignedBases())
	return(-1);
if(Len < (m_NumBasesAligned + m_NumProbeInserts + m_NumTargInserts))
	return(-1);

if(m_bBanded)
	{
	ProbeIdx = m_PeakProbeIdx+1;
	TargIdx =  m_PeakTargIdx+1;
	pPeakCell = DerefBandCell(ProbeIdx,TargIdx);
	}
else
	pPeakCell = &m_pTrcBckCells[((uint64_t)m_PeakProbeIdx * m_TargLen) +  m_PeakTargIdx];
pTargStart = m_pTarg + m_TargAlignStartOfs;
pTarg = pTargStart + m_NumBasesAligned + m_NumTargInserts - 1;
pBuff += m_NumBasesAligned + m_NumProbeInserts + m_NumTargInserts - 1;

do
	{
	TrcBckDir = *pPeakCell  & cSWTrcBckMsk;
	switch(TrcBckDir) {
		case cSWTrcBckDiagFlg:		// back on the diagonal, report base
			*pBuff-- = *pTarg--; 
			if(m_bBanded)
				{
				ProbeIdx -= 1;
				TargIdx -= 1;
				pPeakCell = DerefBandCell(ProbeIdx,TargIdx);
				}
			else
				pPeakCell -= m_TargLen + 1;
			break;

		case cSWTrcBckLeftFlg:			// left, treating as insertion into probe so report InDel
			*pBuff-- = eBaseInDel; 
			if(m_bBanded)
				{
				ProbeIdx -= 1;
				pPeakCell = DerefBandCell(ProbeIdx,TargIdx);
				}
			else			
				pPeakCell -= m_TargLen;
			break;

		case cSWTrcBckDownFlg:		     // down, treating as insertion into target so report as target base
			*pBuff-- = *pTarg--; 
			if(m_bBanded)
				{
				TargIdx -= 1;
				pPeakCell = DerefBandCell(ProbeIdx,TargIdx);
				}
			else
				pPeakCell -= 1;
			break;

		default:				// no direction, must be final cell so report target base
			*pBuff = *pTarg; 
			break;	

		}
	}
while(TrcBckDir);
return(m_NumBasesAligned + m_NumProbeInserts + m_NumTargInserts);
}


int												// returned total alignment length between probe and target including InDels
CSmithWaterman::GetAlignStats(uint32_t *pNumAlignedBases,			// returned number of bases aligning between probe and target
				 uint32_t *pNumExactBases,           // of the aligning bases there were this many exact matches, remainder were substitutions
				 uint32_t *pNumProbeInsertBases,     // this many bases were inserted into the probe relative to the target
				 uint32_t *pNumTargInsertBases,		// this many bases were inserted into the target relative to the probe
				 uint32_t *pProbeStartOfs,			// alignment starts at this probe offset + 1
				 uint32_t *pTargStartOfs)			// alignment starts at this target offset + 1
{
if(!GetNumAlignedBases())
	return(-1);
if(pNumAlignedBases != NULL)
	*pNumAlignedBases = m_NumBasesAligned;
if(pNumExactBases != NULL)
	*pNumExactBases = m_NumBasesExact;
if(pNumProbeInsertBases != NULL)
	*pNumProbeInsertBases = m_NumProbeInserts;
if(pNumTargInsertBases != NULL)
	*pNumTargInsertBases = m_NumTargInserts;
if(pProbeStartOfs != NULL)
	*pProbeStartOfs = m_ProbeAlignStartOfs;
if(pTargStartOfs != NULL)
	*pTargStartOfs = m_TargAlignStartOfs;
return(m_NumBasesAligned + m_NumProbeInserts + m_NumTargInserts);
}

//
// Helper functions for Band cell dereferencing

tSWTrcBckCell *						// returns ptr to newly allocated cell  or NULL if errors
CSmithWaterman::AllocBandCell(uint32_t ProbeBasePsn,	// current probe base position 1..m_ProbeLen
							  uint32_t TargBasePsn)		// current target base position 1..m_TargLen 
{
tsSWColBand *pSWColBand;
tSWTrcBckCell *pTrcBckCell;
#ifdef _DEBUG
if(ProbeBasePsn == 0 || ProbeBasePsn > m_ProbeLen || TargBasePsn == 0 || TargBasePsn > m_TargLen)
	return(NULL);
if((ProbeBasePsn - 1) > m_ColBandsUsed || ProbeBasePsn < m_ColBandsUsed)
	return(NULL);
#endif

if(m_TrcBckCellsUsed == m_TrcBckCellsAllocd)			// need to allocate more?
	{
	size_t ReallocSize;
	uint32_t EstReqCells;
	tSWTrcBckCell *pRealloc;
	EstReqCells = (uint32_t)(((double)(m_ProbeLen * 2) / ProbeBasePsn) * m_TrcBckCellsAllocd); 
	ReallocSize = sizeof(tSWTrcBckCell) * EstReqCells;
#ifdef _WIN32
	pRealloc = (tSWTrcBckCell *)realloc(m_pTrcBckCells,ReallocSize);
#else
	pRealloc = (tSWTrcBckCell *)mremap(m_pTrcBckCells,m_TrcBckCellsAllocdSize,ReallocSize,MREMAP_MAYMOVE);
	if(pRealloc == MAP_FAILED)
		pRealloc = NULL;
#endif
	if(pRealloc == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AllocBandCell: traceback cell memory re-allocation to %I64d bytes - %s",(int64_t)ReallocSize,strerror(errno));
		return(NULL);
		}
	m_TrcBckCellsAllocdSize = ReallocSize;
	m_TrcBckCellsAllocd = EstReqCells;
	m_pTrcBckCells = pRealloc;
	}

pSWColBand = &m_pColBands[ProbeBasePsn - 1];
pTrcBckCell = &m_pTrcBckCells[m_TrcBckCellsUsed++];
*pTrcBckCell = 0;
if(ProbeBasePsn > m_ColBandsUsed)
	{
	m_ColBandsUsed += 1;
	pSWColBand->StartTargBasePsn = TargBasePsn;
	pSWColBand->EndTargBasePsn = TargBasePsn;
	pSWColBand->TrcBckCellPsn = (uint32_t)(uint64_t)m_TrcBckCellsUsed;
	return(pTrcBckCell);
	}
if(pSWColBand->EndTargBasePsn != TargBasePsn - 1)
	{
	m_TrcBckCellsUsed -= 1;
	return(NULL);
	}
pSWColBand->EndTargBasePsn = TargBasePsn;
return(pTrcBckCell);
}

tSWTrcBckCell *					// returns ptr to cell or NULL if errors
CSmithWaterman::DerefBandCell(uint32_t ProbeBasePsn,	// current probe base position 1..m_ProbeLen
							  uint32_t TargBasePsn)		// current target base position 1..m_TargLen
{
uint32_t BandIdx;
tsSWColBand *pSWColBand;
#ifdef _DEBUG
if(m_ColBandsUsed == 0 || ProbeBasePsn == 0 || ProbeBasePsn > m_ColBandsUsed || TargBasePsn == 0 || TargBasePsn > m_TargLen)
	return(NULL);
#endif
pSWColBand = &m_pColBands[ProbeBasePsn - 1];
#ifdef _DEBUG
if(TargBasePsn < pSWColBand->StartTargBasePsn || TargBasePsn > pSWColBand->EndTargBasePsn || pSWColBand->TrcBckCellPsn == 0)
	return(NULL);
#endif
BandIdx = pSWColBand->TrcBckCellPsn - 1 + (TargBasePsn - pSWColBand->StartTargBasePsn);
return(&m_pTrcBckCells[BandIdx]);
}

tSWTrcBckCell *					// returns ptr to cell  or NULL if errors
CSmithWaterman::DerefBandCellLeft(uint32_t ProbeBasePsn,	// current probe base position 1..m_ProbeLen
							  uint32_t TargBasePsn)		// current target base position 1..m_TargLen
{
uint32_t BandIdx;
tsSWColBand *pSWColBand;
#ifdef _DEBUG
if(m_ColBandsUsed == 0 || ProbeBasePsn < 2 || ProbeBasePsn > m_ColBandsUsed || TargBasePsn == 0 || TargBasePsn > m_TargLen)
	return(NULL);
#endif
pSWColBand = &m_pColBands[ProbeBasePsn - 2];
if(TargBasePsn < pSWColBand->StartTargBasePsn || TargBasePsn > pSWColBand->EndTargBasePsn || pSWColBand->TrcBckCellPsn == 0)
	return(NULL);
BandIdx = pSWColBand->TrcBckCellPsn - 1 + (TargBasePsn - pSWColBand->StartTargBasePsn);
return(&m_pTrcBckCells[BandIdx]);
}

tSWTrcBckCell *					// returns ptr to cell  or NULL if errors
CSmithWaterman::DerefBandCellDiag(uint32_t ProbeBasePsn,	// current probe base position 1..m_ProbeLen
							  uint32_t TargBasePsn)		// current target base position 1..m_TargLen
{
uint32_t BandIdx;
tsSWColBand *pSWColBand;
#ifdef _DEBUG
if(m_ColBandsUsed == 0 || ProbeBasePsn < 2 || ProbeBasePsn > m_ColBandsUsed || TargBasePsn < 2 || TargBasePsn > m_TargLen)
	return(NULL);
#endif
TargBasePsn -= 1;
pSWColBand = &m_pColBands[ProbeBasePsn - 2];
if(TargBasePsn < pSWColBand->StartTargBasePsn || TargBasePsn > pSWColBand->EndTargBasePsn || pSWColBand->TrcBckCellPsn == 0)
	return(NULL);
BandIdx = pSWColBand->TrcBckCellPsn - 1 + (TargBasePsn - pSWColBand->StartTargBasePsn);
return(&m_pTrcBckCells[BandIdx]);
}

tSWTrcBckCell *					// returns ptr to cell  or NULL if errors
CSmithWaterman::DerefBandCellDown(uint32_t ProbeBasePsn,	// current probe base position 1..m_ProbeLen
							  uint32_t TargBasePsn)		// current target base position 1..m_TargLen
{
uint32_t BandIdx;
tsSWColBand *pSWColBand;
#ifdef _DEBUG
if(m_ColBandsUsed == 0 || ProbeBasePsn == 0 || ProbeBasePsn > m_ColBandsUsed || TargBasePsn == 0 || TargBasePsn > m_TargLen)
	return(NULL);
#endif
pSWColBand = &m_pColBands[ProbeBasePsn - 1];
if(TargBasePsn < pSWColBand->StartTargBasePsn || TargBasePsn > pSWColBand->EndTargBasePsn || pSWColBand->TrcBckCellPsn == 0)
	return(NULL);
BandIdx = pSWColBand->TrcBckCellPsn - 1 +  (TargBasePsn - pSWColBand->StartTargBasePsn);
return(&m_pTrcBckCells[BandIdx - 1]);
}


int
CSmithWaterman::DumpScores(char *pszFile,		// dump Smith-Waterman matrix to this csv file
					char Down,	// use this char to represent cell down link representing base inserted into target relative to probe
					char Left,	// use this char to represent cell left link representing base inserted into probe relative to target
				    char Diag)	// use this char to represent cell diagonal  representing matching base either exact or mismatch

{
char *pszBuff;
char *pRow;
uint32_t BuffIdx;
uint32_t TargIdx;
uint32_t ProbeIdx;
tSWTrcBckCell *pCell;
tSWTrcBckCell TrcBckDir;
int TrcBckScore;
int hDumpFile;

uint32_t EstDumpLen;

if(m_bBanded)
	{
	gDiagnostics.DiagOut(eDLWarn,gszProcName,"Banded DumpScores not currently supported");
	return(eBSFerrParams);
	}


if(!m_bAligned || m_PeakProbeIdx == 0 || m_PeakTargIdx == 0 || pszFile == NULL || pszFile[0] == '\0')
	return(eBSFerrParams);

EstDumpLen = (uint32_t)min((uint64_t)m_TargLen * m_ProbeLen * 10,(uint64_t)0x03ffff);
if((pszBuff = new char [EstDumpLen])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate dump buffer memory %d bytes",EstDumpLen);
	return(eBSFerrMem);
	}

#ifdef _WIN32
if((hDumpFile = open(pszFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((hDumpFile = open(pszFile,O_RDWR | O_CREAT | O_TRUNC,S_IREAD | S_IWRITE))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create %s - %s",pszFile,strerror(errno));
	return(eBSFerrCreateFile);
	}

// target bases along top as columns, probe bases vertically as rows
pRow = pszBuff;
for(BuffIdx = TargIdx = 0; TargIdx < m_TargLen; TargIdx++)
	{
	BuffIdx += sprintf(&pszBuff[BuffIdx],",\" \",\"%c\"", CSeqTrans::MapBase2Ascii(m_pTarg[TargIdx]));
	if(BuffIdx + 1000 > EstDumpLen)
		{
		CUtility::RetryWrites(hDumpFile,pszBuff,BuffIdx);
		BuffIdx = 0;
		}
	}
pCell = m_pTrcBckCells;
for(ProbeIdx = 0; ProbeIdx < m_ProbeLen; ProbeIdx++)
	{
	BuffIdx += sprintf(&pszBuff[BuffIdx],"\n%c",CSeqTrans::MapBase2Ascii(m_pProbe[ProbeIdx]));
	for(TargIdx = 0; TargIdx < m_TargLen; TargIdx++,pCell++)
		{
		TrcBckDir = *pCell & cSWTrcBckMsk;
		TrcBckScore = *pCell & cSWScoreMsk;
		if(TrcBckDir == cSWTrcBckDiagFlg)
			BuffIdx += sprintf(&pszBuff[BuffIdx],",\"%c\",%d",Diag,TrcBckScore);
		else 
			{
			if(TrcBckDir == cSWTrcBckDownFlg)
				BuffIdx += sprintf(&pszBuff[BuffIdx],",\"%c\",%d",Down,TrcBckScore);
			else
				{
				if(TrcBckDir == cSWTrcBckLeftFlg)
					BuffIdx += sprintf(&pszBuff[BuffIdx],",\"%c\",%d",Left,TrcBckScore);
				else
					BuffIdx += sprintf(&pszBuff[BuffIdx],",\" \",%d",TrcBckScore);
				}
			}
		if(BuffIdx + 1000 > EstDumpLen)
			{
			CUtility::RetryWrites(hDumpFile,pszBuff,BuffIdx);
			BuffIdx = 0;
			}
		}
	}
if(BuffIdx)
	CUtility::RetryWrites(hDumpFile,pszBuff,BuffIdx);
#ifdef _WIN32
_commit(hDumpFile);
#else
fsync(hDumpFile);
#endif
close(hDumpFile);
return(eBSFSuccess);
}



