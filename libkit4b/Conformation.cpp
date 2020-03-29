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

CConformation::CConformation(void)
{
m_pOctStructParams = NULL;
m_pDiStructParams = NULL;
m_NumDiStructParams = 0;
}

CConformation::~CConformation(void)
{
if(m_pOctStructParams != NULL)
	delete m_pOctStructParams;
if(m_pDiStructParams != NULL)
	delete m_pDiStructParams;
}

// StructParamsLoaded
// Returns true if structural parameters loaded
bool
CConformation::StructParamsLoaded(void)
{
return(m_pOctStructParams == NULL ? false : true);
}

// LoadStructOctamersParams
// Load structural parameters from file, two different sets of structural parameters exist - one (the most comprehensive) is for dimers, the other is for octamers
// For octamers the file layout is that of one comma separated set of parameters per unique octamer line with optional heading line
// Octamer,Twist,Roll,Tilt,Rise,Slide,Shift,3-StepTwist,3-StepRoll,3-StepSlide,3-StepShift,Energy,
// MinorGroove,RMSD,Q-Twist,Q+Twist,Q-Roll,Q+Roll,3Q-Twist,3Q+Twist,3Q-Roll,3Q+Roll, ORChID (hydroxyl radical cleavage values)
//
teBSFrsltCodes 
CConformation::LoadStructOctamersParams(char *pszStructParamsFile) // load structural parameters from file
{
FILE *pParamsStream;
int LineNum;
int NumOctParams;
char szLineBuff[512];
int Cnt;
char Octamer[9];
etSeqBase octbases[8];
int OctIdx;
double Twist,Roll,Tilt,Rise,Slide,Shift,TriStepTwist,TriStepRoll,TriStepSlide,TriStepShift;
double Energy,MinorGroove,RMSD,QminusTwist,QplusTwist,QminusRoll,QplusRoll;
double TriQminusTwist,TriQplusTwist,TriQminusRoll,TriQplusRoll,ORChID;

tsOctStructParam *pStruct1;
char chr, *pDst, *pSrc;

if(m_pOctStructParams == NULL)
	{	
	m_pOctStructParams = (tsOctStructParam *) new UINT8[cOctStructParamAllocSize];
	if(m_pOctStructParams == NULL)
		{
		AddErrMsg("CConformation::LoadStructParams","Unable to allocate memory to hold structural parameters");
		return(eBSFerrMem);
		}
	memset(m_pOctStructParams,0,cOctStructParamAllocSize);
	}
m_szOctStructParamFile[0] = '\0';


if((pParamsStream = fopen(pszStructParamsFile,"r"))==NULL)
	{
	AddErrMsg("CConformation::LoadStructParams","Unable to open parameters file %s error: %s",pszStructParamsFile,strerror(errno));
	delete m_pOctStructParams;
	m_pOctStructParams = NULL;
	return(eBSFerrOpnFile);
	}


LineNum = 0;
NumOctParams = 0;
int NxScanPsn;
int ParamLines = 0;
while(fgets(szLineBuff,sizeof(szLineBuff),pParamsStream)!= NULL)
	{
	LineNum++;
	if(strlen(szLineBuff) < 5)	// simply slough lines which are too short to contain anything worth parsing
		continue;

	// strip any whitespace and quotes
	pDst = pSrc = szLineBuff;
	while(chr = *pSrc++)
		if(!isspace(chr) && chr != '\'' && chr != '"')
			*pDst++ = chr;
	*pDst = '\0';
	if(szLineBuff[0] == '\0')
		continue;
	ParamLines += 1;
	// if unable to initially parse the first line then assume it's a header line and slough
	 Cnt = sscanf(szLineBuff,"%8[^,],%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf %n",
			Octamer,&Twist,&Roll,&Tilt,&Rise,&Slide,&Shift,&TriStepTwist,&TriStepRoll,&TriStepSlide,&TriStepShift, &NxScanPsn);
	 if(Cnt < 11 && ParamLines == 1)
		 continue;

	 if(Cnt == 11)
		 Cnt += sscanf(&szLineBuff[NxScanPsn],",%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
				&Energy,&MinorGroove,&RMSD,&QminusTwist,&QplusTwist,&QminusRoll,&QplusRoll,&TriQminusTwist,&TriQplusTwist,&TriQminusRoll,&TriQplusRoll, &ORChID);

	 if(Cnt != 23)
		{
		AddErrMsg("CConformation::LoadStructParams","Error parsing structural parameters file %s at line %d, expected 23 but only parsed %d parameters\n%s\n",pszStructParamsFile,LineNum,Cnt,szLineBuff);
		fclose(pParamsStream);
		delete m_pOctStructParams;
		m_pOctStructParams = NULL;
		return(eBSFerrStructParm); 
		}

	CSeqTrans::MapAscii2Sense(Octamer,8,octbases);
	if((OctIdx = StructParamIdx(octbases)) < 0)
		{
		AddErrMsg("CConformation::LoadStructParams","Octamer in file %s at line %d contains unrecognised bases (only 'a','c','g','t' accepted)\n%s\n",pszStructParamsFile,LineNum,szLineBuff);
		fclose(pParamsStream);
		delete m_pOctStructParams;
		m_pOctStructParams = NULL;
		return((teBSFrsltCodes)OctIdx); 
		}
	pStruct1	 = &m_pOctStructParams[OctIdx];
	if(pStruct1->Param.twist != 0)
		{
		fclose(pParamsStream);
		delete m_pOctStructParams;
		AddErrMsg("CConformation::LoadStructParams","Duplicate octamer structural properties parsed at line %d in file %s\n",LineNum,pszStructParamsFile);
		return(eBSFerrStructParm);
		}

	pStruct1->Param.energy = (int)(Energy * 10000.0);
	pStruct1->Param.minorgroove = (int)(MinorGroove * 10000.0);
	pStruct1->Param.twist =  (int)(Twist * 10000.0);
	pStruct1->Param.roll = (int)(Roll * 10000.0);
	pStruct1->Param.tilt = (int)(Tilt * 10000.0);
	pStruct1->Param.rise = (int)(Rise * 10000.0);
	pStruct1->Param.slide = (int)(Slide * 10000.0);
	pStruct1->Param.shift = (int)(Shift * 10000.0);
	pStruct1->Param.rmsd = (int)(RMSD * 10000.0);
	pStruct1->Param.orchid = (int)(ORChID * 10000.0);
	NumOctParams++;
	}
fclose(pParamsStream);

if(NumOctParams != cNumParamOctamers)
	AddErrMsg("CConformation::LoadStructParams","Warning, not all octamers in '%s' have structural properties\n",pszStructParamsFile);
strncpy(m_szOctStructParamFile,pszStructParamsFile,_MAX_PATH-1);
m_szOctStructParamFile[_MAX_PATH-1] = '\0';

// inference mid-step major groove from twist + rise
int StepIdx;
int SumTwists;
int SumRises;
double MajorGroove;

memset(octbases,0,sizeof(octbases));
pStruct1 = m_pOctStructParams;
for(OctIdx = 0; OctIdx < cNumParamOctamers; OctIdx++, pStruct1++)
	{
	for(StepIdx=0; StepIdx < 8; StepIdx++)
		octbases[StepIdx] = 0x03 & (OctIdx >> (StepIdx*2));
	SumTwists = 0;
	SumRises = 0;
	for(StepIdx=1; StepIdx <= 7; StepIdx++)		// NOTE: conformational steps are interbase so octamers contain 7 steps
		{
		SumTwists += StructValue(eSStwist,StepIdx,8,octbases,0);
		SumRises += StructValue(eSSrise,StepIdx,8,octbases,0);
		}

	MajorGroove = (double)SumRises * 10000.0 * 360.0/(double)SumTwists;
	
	pStruct1->Param.majorgroove = (int)MajorGroove - pStruct1->Param.minorgroove;
	} 
return(eBSFSuccess);
}

// dimer structural parameters, generated from 'DiProDB: a database for dinucleotide properties', http://diprodb.fli-leibniz.de
teBSFrsltCodes 
CConformation::LoadStructDimersParams(char *pszStructParamsFile) // load structural parameters from file
{
FILE *pParamsStream;
int LineNum;
int NumFields;
char szLineBuff[512];

tsDiStructParam DiStructParam;
int DiMerIdx;

tsDiStructParam *pDiStructParam;
int NumDiStructParams;
char chr, *pDst, *pSrc;

if(m_pDiStructParams == NULL)
	{	
	m_pDiStructParams = (tsDiStructParam *) new unsigned char[cDiStructParamAllocSize];
	if(m_pDiStructParams == NULL)
		{
		AddErrMsg("CConformation::LoadStructDimersParams","Unable to allocate memory to hold structural parameters");
		return(eBSFerrMem);
		}
	memset(m_pDiStructParams,0,cDiStructParamAllocSize);
	}
m_szDiStructParamFile[0] = '\0';
m_NumDiStructParams = 0;
NumDiStructParams = 0;

if((pParamsStream = fopen(pszStructParamsFile,"r"))==NULL)
	{
	AddErrMsg("CConformation::LoadStructParams","Unable to open dimer structural parameters file %s error: %s",pszStructParamsFile,strerror(errno));
	delete m_pDiStructParams;
	m_pDiStructParams = NULL;
	return(eBSFerrOpnFile);
	}

pDiStructParam = m_pDiStructParams;
LineNum = 0;
int NxScanPsn;
int ParamLines = 0;
while(fgets(szLineBuff,sizeof(szLineBuff),pParamsStream)!= NULL)
	{
	LineNum++;
	if(strlen(szLineBuff) < 10)	// simply slough lines which are too short to contain anything worth parsing
		continue;

	// strip any whitespace and quotes
	pDst = pSrc = szLineBuff;
	while(chr = *pSrc++)
		if(!isspace(chr) && chr != '\'' && chr != '"')
			*pDst++ = chr;
	*pDst = '\0';
	if(szLineBuff[0] == '\0')
		continue;
	ParamLines += 1;

	// expecting at least 18 fields in CSV rows
	pSrc = szLineBuff; 
	// parse out expected structural parameter identifier followed by parameter name
	NumFields = sscanf(pSrc,"%d,%80[^,],%n",&DiStructParam.ID,DiStructParam.szName,&NxScanPsn);
	if(NumFields != 2)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to parse dimer structural parameter identifier and name from line %d in '%s' file",LineNum,pszStructParamsFile);
		delete m_pDiStructParams;
		m_pDiStructParams = NULL;
		fclose(pParamsStream);
		return(eBSFerrConfParam);
		}
	// check that the structural parameter is unique and monotonically incrementing
	for(DiMerIdx = 1; DiMerIdx <= NumDiStructParams; DiMerIdx++)
		if(DiStructParam.ID != DiMerIdx)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Dimer structural parameter identifier %d parsed from line %d in '%s' file is not monotonically incrementing from 1 to 124",DiStructParam.ID,LineNum,pszStructParamsFile);
			delete m_pDiStructParams;
			m_pDiStructParams = NULL;
			fclose(pParamsStream);
			return(eBSFerrConfParam);
			}	
	// next expecting the parameter values for each dimer - 16 dimers AA..TT
	pSrc += NxScanPsn;
	for(DiMerIdx=0;DiMerIdx < 16; DiMerIdx++)
		{
		NumFields = sscanf(pSrc,"%lf,%n",&DiStructParam.DiValue[DiMerIdx],&NxScanPsn);
		if(NumFields != 1)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to parse dimer structural parameter value for dimer %d from line %d in '%s' file",DiMerIdx, LineNum,pszStructParamsFile);
			delete m_pDiStructParams;
			m_pDiStructParams = NULL;
			fclose(pParamsStream);
			return(eBSFerrConfParam);
			}
		pSrc += NxScanPsn;
		}

	// finally, treat all remaining fields as supplementary
	if(*pSrc != '\0')
		strncpy(DiStructParam.szSupFields,pSrc,sizeof(DiStructParam.szSupFields)-1);
	*pDiStructParam++ = DiStructParam; 
	NumDiStructParams++;
	}
if(NumDiStructParams != 124)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to parse all expected 124 dimer structural parameters from '%s' file",pszStructParamsFile);
	delete m_pDiStructParams;
	m_pDiStructParams = NULL;
	fclose(pParamsStream);
	return(eBSFerrConfParam);
	}
fclose(pParamsStream);
strncpy(m_szDiStructParamFile,pszStructParamsFile,_MAX_PATH-1);
m_szDiStructParamFile[_MAX_PATH-1] = '\0';
m_NumDiStructParams = NumDiStructParams;
return(eBSFSuccess);
}

// reproducible psudeo-randomise conformation characteristics
//
int
CConformation::PsudeoRandomise(void)
{
int Len;
int OctIdx;
tsOctStructParam tmp;
CRandomMersenne Random(1); // random sequence must be reproducible so seed required
for(Len = cNumParamOctamers; Len > 1; Len--)
	{				
    OctIdx       = Random.IRandom(0,Len-1);
    tmp         = m_pOctStructParams[OctIdx];
    m_pOctStructParams[OctIdx]   = m_pOctStructParams[Len-1];
    m_pOctStructParams[Len-1] = tmp;
	}
return(eBSFSuccess);
}


// determine index to use for assumed octamer sequence
// eBSFerrStructStep returned if any base in octamer is indeterminate - 'N'
int
CConformation::StructParamIdx(etSeqBase *pOctamer)		// sequence
{
etSeqBase Base;
int OctIdx = 0;
int Len = 8;
if(pOctamer == NULL || m_pOctStructParams == NULL)
	return((int)eBSFerrParams);

while(Len--)
	{
	OctIdx <<= 2;
	Base = *pOctamer++  & ~cRptMskFlg;
	if(Base > eBaseT)
		return((int)eBSFerrStructStep);								// unrecognised base
	OctIdx |= Base;
	}
return(OctIdx);
}

teBSFrsltCodes
CConformation::GetSequenceConformation(teOctStructStats Param,		// which structural parameter value to return
				 unsigned int iStartOfs, // initial starting offset (0..n) in pSeq
				  unsigned int iNumSteps,		  // number of steps (0 for all) to process starting at pSeq[iStartPsn]|pSeq[iStartPsn+1]
  				  unsigned int SeqLen,			  // total length of sequence
				  etSeqBase *pSeq,				  // sequence to be processed
				  int *pRetConfValue,			  // where to return conformation
				  int UndefBaseValue)			  // value to return for undefined or indeterminate ('N') bases 
{
unsigned int Step;
unsigned int LastStep;
int IdxRetVal;

if(SeqLen < 8 || iStartOfs >= SeqLen - 1 || 
   iStartOfs + iNumSteps >= SeqLen ||
   pSeq == NULL || m_pOctStructParams == NULL)
	return(eBSFerrParams);

if(iNumSteps == 0)
	iNumSteps = SeqLen - iStartOfs - 1;
LastStep = iStartOfs + iNumSteps;

for(IdxRetVal =0,Step = iStartOfs+1; Step <= LastStep; IdxRetVal++,Step++)
	pRetConfValue[IdxRetVal] = StructValue(Param,Step,SeqLen,pSeq,UndefBaseValue);
return(eBSFSuccess);
}

int
CConformation::StructValue(teOctStructStats Param,		// which structural parameter value to return
			unsigned int Step,			// which step in sequence to return structural value for
			unsigned int SeqLen,		// total length of sequence
			etSeqBase *pSeq,			// sequence to be processed
			int UndefBaseValue)			// value to return for undefined or indeterminate ('N') bases 
{
tsOctStructParam *pStruct;
etSeqBase *pOctamer;
etSeqBase octamer[8];
int OctIdx;
unsigned int InterpStep;
int ParamOfs;
bool bContinue;
bool bLeft;
int Iters;
INT64 InterpValue;
bContinue = false;
int	Cnt;
int	OctOfs;
int	SeqOfs;

switch(Param) {
	case eSSenergy:				// minimal energy int(energy * 10000) e.g. -408.2632 ==> -4082632
		ParamOfs = offsetof(tsOctStructParam,Param.energy);
		break;

	case eSSminorgroove:				// groove int(dimensions * 10000) e.g 10.784 ==> 107840
		ParamOfs = offsetof(tsOctStructParam,Param.minorgroove);
		break;

	case eSSmajorgroove:				
		ParamOfs = offsetof(tsOctStructParam,Param.majorgroove);
		break;

	case eSStwist:					// twist int(angle * 10000) e.g 	37.6262 ==> 376262
		ParamOfs = offsetof(tsOctStructParam,Param.twist);
		break;

	case eSSroll:					// roll int(angle * 10000) e.g 	2.4139 ==> 24139
		ParamOfs = offsetof(tsOctStructParam,Param.roll);
		break;

	case eSStilt:					// tilt int(tilt * 10000) e.g 	-0.022 ==> -220
		ParamOfs = offsetof(tsOctStructParam,Param.tilt);
		break;

	case eSSrise:					// rise int(rise * 10000) e.g 	3.1409 ==> 31409
		ParamOfs = offsetof(tsOctStructParam,Param.rise);
		break;

	case eSSslide:					// slide int(slide * 10000) e.g 	-0.0968 ==> -968	
		ParamOfs = offsetof(tsOctStructParam,Param.slide);
		break;

	case eSSshift:					// shift int(shift * 10000) e.g 	0.0645 ==> 645
		ParamOfs = offsetof(tsOctStructParam,Param.shift);
		break;

	case eSSrmsd:					// rmsd int(rmsd * 10000) e.g 	0.3078 ==> 3078
		ParamOfs = offsetof(tsOctStructParam,Param.rmsd);
		break;

	case eSSORChidVal:
		ParamOfs = offsetof(tsOctStructParam,Param.orchid);
		break;

	default:
		return(eBSFerrParams);
	};

// if can use the midstep of an octamer then do so...
if(Step > 3 && Step < SeqLen - 3)
	{
	pOctamer = pSeq + Step - 4;
	if((OctIdx = StructParamIdx(pOctamer)) < 0)
		return(UndefBaseValue);
	else
		{
		pStruct = &m_pOctStructParams[OctIdx];
		return(*(int *)(((unsigned char *)pStruct)+ParamOfs));
		}
	}

	// not a octamer midstep, need to interpolate
	// step    cnt   toofs fromofs
	// 1	   5     3		3
	// 2       6     2		2
	// 3       7     1      1
	//
	// 5       7     0		SeqLen - 8
	// 6       6	 0      SeqLen - 7
	// 7       5	 0		SeqLen - 6
InterpStep = Step;
memset(octamer,eBaseA,8);
if(InterpStep < 4)
	{
	Cnt = 4 + InterpStep;
	OctOfs = 4 - InterpStep;
	SeqOfs = 0;
	}
else
	{
	InterpStep = 8 - (SeqLen - InterpStep);						// adjust Step to 5..7
	Cnt = 12 - InterpStep;
	OctOfs = 0;
	SeqOfs = SeqLen - (Cnt + 1);
	}
memmove(&octamer[OctOfs],&pSeq[SeqOfs],Cnt);

//		Interpolate(Step,octamer,pRetStructParam);
if(InterpStep > 4)								// sequence to left is known
	{
	bLeft = true;
	InterpStep -= 4;						
	}
else										// else sequence to right is known
	{
	bLeft = false;
	InterpStep = 4 - InterpStep;	
	}

switch(InterpStep) {
	case 1: 
		Iters = 4;
		break;
	case 2:
		Iters = 16;
		break;
	default:
		Iters = 64;
		break;
	}

InterpValue=0;
for(Cnt = 0; Cnt < Iters; Cnt++)
	{
	OctIdx = StructParamIdx(octamer);
	if(OctIdx < 0)
		return(UndefBaseValue);
	pStruct = &m_pOctStructParams[OctIdx];
	InterpValue += *(int *)(((unsigned char *)pStruct)+ParamOfs);;

	if(!bLeft)
		{
		octamer[0]++;
		if(octamer[0] > eBaseT)
			{
			octamer[0] = eBaseA;
			octamer[1]++;
			if(octamer[1] > eBaseT)
				{
				octamer[1] = eBaseA;
				octamer[2]++;
				}
			}
		}
	else
		{
		octamer[7]++;
		if(octamer[7] > eBaseT)
			{
			octamer[7] = eBaseA;
			octamer[6]++;
			if(octamer[6] > eBaseT)
				{
				octamer[6] = eBaseA;
				octamer[5]++;
				}
			}
		}
	}
InterpValue /= Iters;
return((int)InterpValue);
}



	