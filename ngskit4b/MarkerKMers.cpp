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

#include "./ngskit4b.h"
#include "MarkerKMers.h"

static tsPutMarker *gpPutativeMarkers = nullptr;		// used when sorting putative marker sequences
static size_t gPutativeMarkerSeqLen = 0;			// length of putative marker sequences
static size_t gPutMarkerSize = 0;					// size of a tsPutMarker including the putative sequence

CMarkerKMers::CMarkerKMers(void)
{
m_pSfxArray = nullptr;
m_pMarkerBuff = nullptr;
m_pPutMarkers = nullptr;
m_pPutMarkersIndex = nullptr;
m_hOutFile = -1;
#ifdef _WIN32
InitializeSRWLock(&m_hRwLock);
#else
pthread_rwlock_init( &m_hRwLock,nullptr);
#endif

#ifdef _WIN32
InitializeCriticalSectionAndSpinCount(&m_hSCritSect,1000);
#else
pthread_spin_init(&m_hSpinLock,PTHREAD_PROCESS_PRIVATE);
#endif
Reset();
}


CMarkerKMers::~CMarkerKMers(void)
{
if(m_hOutFile != -1)
	close(m_hOutFile);
if(m_pSfxArray != nullptr)
	delete m_pSfxArray;
if(m_pMarkerBuff != nullptr)
	delete m_pMarkerBuff;
if(m_pPutMarkers != nullptr)
	{
#ifdef _WIN32
	free(m_pPutMarkers);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pPutMarkers != MAP_FAILED)
		munmap(m_pPutMarkers,m_AllocPutMarkersSize);
#endif
	}
if(m_pPutMarkersIndex != nullptr)
	{
#ifdef _WIN32
	free(m_pPutMarkersIndex);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pPutMarkersIndex != MAP_FAILED)
		munmap(m_pPutMarkersIndex,m_AllocPutMarkersIndexSize);
#endif
	}

#ifndef _WIN32
pthread_rwlock_destroy(&m_hRwLock);
pthread_spin_destroy(&m_hSpinLock);
#endif
}

void
CMarkerKMers::Reset(bool bSync)
{
if(m_pSfxArray != nullptr)
	{
	delete m_pSfxArray;
	m_pSfxArray = nullptr;
	}

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

if(m_pMarkerBuff != nullptr)
	{
	delete m_pMarkerBuff;
	m_pMarkerBuff = nullptr;
	}

if(m_pPutMarkers != nullptr)
	{
#ifdef _WIN32
	free(m_pPutMarkers);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pPutMarkers != MAP_FAILED)
		munmap(m_pPutMarkers,m_AllocPutMarkersSize);
#endif
	m_pPutMarkers = nullptr;
	}

if(m_pPutMarkersIndex != nullptr)
	{
#ifdef _WIN32
	free(m_pPutMarkersIndex);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pPutMarkersIndex != MAP_FAILED)
		munmap(m_pPutMarkersIndex,m_AllocPutMarkersIndexSize);
#endif
	m_pPutMarkersIndex = nullptr;
	}

m_szDataset[0] = '\0';
m_szMarkerFile[0] = '\0';
m_NumSfxEntries = 0;
m_NumPrefixKMers = 0;					
m_TotSenseCnts = 0;			
m_TotAntisenseCnts = 0;	
m_MarkerID = 0;		
m_MarkerBuffOfs = 0;			
m_AllocMarkerBuffSize = 0;	
m_PutMarkerSize = 0;
m_NumPutMarkers = 0;
m_AllocPutMarkersSize = 0;
m_AllocPutMarkersIndexSize = 0;
memset(m_AllCultivars,0,sizeof(m_AllCultivars));
}


// serialise access to critical code
inline void
CMarkerKMers::EnterCritSect(void)
{
int SpinCnt = 5000;
#ifdef _WIN32
while(!TryEnterCriticalSection(&m_hSCritSect))
	{
	if(SpinCnt -= 1)
		continue;
	SwitchToThread();
	SpinCnt = 500;
	}
#else
while(pthread_spin_trylock(&m_hSpinLock)==EBUSY)
	{
	if(SpinCnt -= 1)
		continue;
	sched_yield();
	SpinCnt = 500;
	}
#endif
}

inline void
CMarkerKMers::LeaveCritSect(void)
{
#ifdef _WIN32
LeaveCriticalSection(&m_hSCritSect);
#else
pthread_spin_unlock(&m_hSpinLock);
#endif
}

// AcquireLock
// Aquire lock, either exclusive or read only, on class shared instance vars
void 
CMarkerKMers::AcquireLock(bool bExclusive)
{
#ifdef _WIN32
if(bExclusive)
	AcquireSRWLockExclusive(&m_hRwLock);
else
	AcquireSRWLockShared(&m_hRwLock);
#else
if(bExclusive)
	pthread_rwlock_wrlock(&m_hRwLock);
else
	pthread_rwlock_rdlock(&m_hRwLock);
#endif
}


// ReleaseLock
// Release lock, must be exclusive or read only as used when AcquireLock was used
void
CMarkerKMers::ReleaseLock(bool bExclusive)
{
#ifdef _WIN32
if(bExclusive)
	ReleaseSRWLockExclusive(&m_hRwLock);
else
	ReleaseSRWLockShared(&m_hRwLock);
#else
pthread_rwlock_unlock(&m_hRwLock);
#endif
}

	
int
CMarkerKMers::ReportMarker(tsPutMarker *pMarker)
{
int TruncMarkerLen;
char szMarkerSeq[cMaxKMerLen+1];				// additional for '\0' terminator

if(m_hOutFile == -1 || pMarker == nullptr || m_pMarkerBuff == nullptr)	// better safe than sorry...
	return(0);

TruncMarkerLen = m_PrefixLen > cMaxKMerLen ? cMaxKMerLen : m_PrefixLen;
CSeqTrans::MapSeq2UCAscii(pMarker->MarkerSeq,TruncMarkerLen,szMarkerSeq);
AcquireLock(true);
if((m_MarkerBuffOfs + (2 * cMaxKMerLen)) > m_AllocMarkerBuffSize)
	{
	CUtility::RetryWrites(m_hOutFile,m_pMarkerBuff,m_MarkerBuffOfs);
	m_MarkerBuffOfs = 0;
	}
m_MarkerID += 1;
m_MarkerBuffOfs += sprintf((char *)&m_pMarkerBuff[m_MarkerBuffOfs],">KMerM%d %d|%d|%d|%u|%u\n%s\n",m_MarkerID,TruncMarkerLen,m_KMerLen,pMarker->NumCultivars,pMarker->SenseCnts,pMarker->AntisenseCnts,szMarkerSeq);

ReleaseLock(true);
return(eBSFSuccess);
}


int64_t														// returns number of K-Mers processed
CMarkerKMers::GetKMerProcProgress(int64_t *pTotSenseCnts,	// number of K-Mers accepted on sense strand
					int64_t *pTotAntisenseCnts)				// number of K-Mers accepted on antisense strand
{
int64_t NumPutMarkers;
EnterCritSect();
if(pTotSenseCnts != nullptr)
	*pTotSenseCnts = m_TotSenseCnts;
if(pTotAntisenseCnts != nullptr)
	*pTotAntisenseCnts = m_TotAntisenseCnts;
NumPutMarkers = m_NumPutMarkers;
LeaveCritSect();
return(NumPutMarkers);
}

int
CMarkerKMers::LocKMers(etPMode PMode,			// processing mode - defaults to 0
		  int KMerLen,					// this length K-mers
		  int PrefixLen,				// inter-cultivar shared prefix length
		  int SuffixLen,				// cultivar specific suffix length
		  int MinWithPrefix,			// minimum number of cultivars required to have the shared prefix
		  int MaxHomozygotic,			// only report prefixes if K-Mer suffixes are homozygotic between a maximum of this many cultivars, if 1 then no other cultivars
		  char *pszSfxPseudoGenome,		// contains pregenerated suffix over psuedochromosomes for each cultivar
		  char *pszMarkerFile,			// output potential markers to this file
		  int NumThreads)				// max number of threads allowed
{
int Rslt;
int64_t NumPutativePrefixKMers;
int64_t TotSenseCnts;
int64_t TotAntisenseCnts;

int EntryID;
char szSfxEntryName[100];
tsCultivar *pCultivar;

Reset();

m_PMode = PMode;
m_KMerLen = KMerLen;
m_PrefixLen = PrefixLen;
m_SuffixLen = SuffixLen;
m_MinWithPrefix = MinWithPrefix;
m_MaxHomozygotic = MaxHomozygotic;
m_PutMarkerSize = (int)sizeof(tsPutMarker) + m_PrefixLen - 1;

m_NumThreads = NumThreads;
strncpy(m_szMarkerFile,pszMarkerFile,sizeof(m_szMarkerFile));
m_szMarkerFile[sizeof(m_szMarkerFile)-1] = '\0';

// allocate buffers
if((m_pMarkerBuff = new uint8_t [cMarkerSeqBuffSize])==nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Unable to allocate (%d bytes) for marker buffering",cMarkerSeqBuffSize);
	Reset();
	return(eBSFerrMem);
	}
m_AllocMarkerBuffSize = cMarkerSeqBuffSize;


	// allocate initial putative marker sequence memory
m_AllocPutMarkersSize = (size_t)cAllocNumPutativeSeqs * m_PutMarkerSize;
#ifdef _WIN32
m_pPutMarkers = (tsPutMarker *) malloc(m_AllocPutMarkersSize);
if(m_pPutMarkers == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to allocate %zd bytes contiguous memory for putative marker sequences",(int64_t)m_AllocPutMarkersSize);
	m_AllocPutMarkersSize = 0;
	Reset(false);
	return(eBSFerrMem);
	}
#else
// gnu malloc is still in the 32bit world and seems to have issues if more than 2GB allocation
m_pPutMarkers = (tsPutMarker *)mmap(nullptr,m_AllocPutMarkersSize, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
if(m_pPutMarkers == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to allocate %zd bytes contiguous memory for putative marker sequences",(int64_t)m_AllocPutMarkersSize);
	m_AllocPutMarkersSize = 0;
	m_pPutMarkers = nullptr;
	Reset(false);
	return(eBSFerrMem);
	}
#endif
m_NumPutMarkers = 0;

if((m_pSfxArray = new CSfxArray)==nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Unable to instantiate instance of CSfxArray");
	Reset();
	return(eBSFerrObj);
	}

if((Rslt = m_pSfxArray->Open(pszSfxPseudoGenome))!=eBSFSuccess)
	{
	while(m_pSfxArray->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pSfxArray->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open input bioseq suffix array file '%s'",pszSfxPseudoGenome);
	Reset();
	return(Rslt);
	}

// report to user some sfx array metadata for user conformation the targeted assembly is correct
strcpy(m_szDataset,m_pSfxArray->GetDatasetName());
if(m_pSfxArray->IsSOLiD())
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to process SOLiD colorspace suffix array file '%s'",pszSfxPseudoGenome);
	Reset();
	return(eBSFerrRefDataset);
	}
tsSfxHeaderV3 SfxHeader;
m_pSfxArray->GetSfxHeader(&SfxHeader);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Psuedo-assembly Name: '%s' Descr: '%s' Title: '%s' Version: %d",
					 m_szDataset,SfxHeader.szDescription,SfxHeader.szTitle,SfxHeader.Version);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Assembly block size: %zu",SfxHeader.SfxBlockSize);

m_NumSfxEntries = m_pSfxArray->GetNumEntries();
if(m_NumSfxEntries < 1)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"No pseudo-chroms in '%s'",pszSfxPseudoGenome);
	Reset();
	return(eBSFerrEntry);
	}

if(m_NumSfxEntries < 2 || m_NumSfxEntries > cMaxCultivars)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Number of pseudo-chroms (%d) in '%s' must be between 2 and %d",m_NumSfxEntries,pszSfxPseudoGenome,cMaxCultivars);
	Reset();
	return(eBSFerrEntry);
	}

// ensure assembly and suffix array has been fully loaded into memory 
int CurBlockID = 1;		// currently only single block supported
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading genome assembly suffix array...");
if((Rslt=m_pSfxArray->SetTargBlock(CurBlockID))<0)
	{
	while(m_pSfxArray->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pSfxArray->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to load assembly suffix array");
	return(Rslt);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Assembly suffix array loaded");


// now report whilst populating m_AllCultivars[]
pCultivar = m_AllCultivars;
for(EntryID = 1; EntryID <= m_NumSfxEntries; EntryID++, pCultivar++)
	{
	pCultivar->Status = 0;
	m_pSfxArray->GetIdentName(EntryID,sizeof(szSfxEntryName),szSfxEntryName);
	pCultivar->EntryID = EntryID;
	strncpy(pCultivar->szEntryName,szSfxEntryName,sizeof(pCultivar->szEntryName));
	pCultivar->szEntryName[sizeof(pCultivar->szEntryName)-1] = '\0';
	pCultivar->EntryLen = m_pSfxArray->GetSeqLen(EntryID);
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"   Processing K-mers against pseudo-chromosome '%s'",szSfxEntryName);
	}

if(m_MinWithPrefix == 0 || m_MinWithPrefix > m_NumSfxEntries)
	m_MinWithPrefix = m_NumSfxEntries;

// looks good to go so create/truncate output marker sequence file
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Truncating/Creating output multi-fasta marker file '%s'",m_szMarkerFile);
#ifdef _WIN32
m_hOutFile = open(m_szMarkerFile,O_CREATETRUNC );
#else
if((m_hOutFile = open(m_szMarkerFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
	if(ftruncate(m_hOutFile,0)!=0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",m_szMarkerFile,strerror(errno));
			Reset();
			return(eBSFerrCreateFile);
			}
#endif

if(m_hOutFile < 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output file '%s'",m_szMarkerFile);
	m_hOutFile = -1;
	Reset();
	return(eBSFerrCreateFile);
	}

// initialise and startup K-mer processing worker threads
tsKMerThreadPars WorkerThreads[cMaxWorkerThreads];			// allow for max possible user configured number of threads
int ThreadIdx;
memset(WorkerThreads,0,sizeof(WorkerThreads));
int NumActiveThreads;
int64_t StartSfxIdx;
int64_t EndSfxIdx;

// partition the processing over multiple threads
StartSfxIdx = 0;
for(NumActiveThreads = 0; NumActiveThreads < m_NumThreads; NumActiveThreads++)
	{
	if((Rslt = m_pSfxArray->GenKMerCultThreadRange(PrefixLen,NumActiveThreads+1,m_NumThreads,StartSfxIdx,&EndSfxIdx))<1)
		break;
	WorkerThreads[NumActiveThreads].StartSfxIdx = StartSfxIdx;
	WorkerThreads[NumActiveThreads].EndSfxIdx = EndSfxIdx;
	StartSfxIdx = EndSfxIdx + 1;
	}
if(Rslt < 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to partition processing between threads");
	Reset();
	return(Rslt);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Process: Partitioning processing between %d threads",NumActiveThreads);
for(ThreadIdx = 0; ThreadIdx < NumActiveThreads; ThreadIdx++)
	{
	WorkerThreads[ThreadIdx].ThreadIdx = ThreadIdx + 1;
	WorkerThreads[ThreadIdx].pThis = this;
#ifdef _WIN32
	WorkerThreads[ThreadIdx].threadHandle = (HANDLE)_beginthreadex(nullptr,0x0fffff,KMerThreadStart,&WorkerThreads[ThreadIdx],0,&WorkerThreads[ThreadIdx].threadID);
#else
	WorkerThreads[ThreadIdx].threadRslt =	pthread_create (&WorkerThreads[ThreadIdx].threadID , nullptr , KMerThreadStart , &WorkerThreads[ThreadIdx] );
#endif
	}

// allow threads a few seconds to startup
#ifdef _WIN32
	Sleep(10000);
#else
	sleep(10);
#endif

// let user know that this K-mer processing process is working hard...
NumPutativePrefixKMers = GetKMerProcProgress(&TotSenseCnts,&TotAntisenseCnts);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress - putative prefix K-Mers: %zd",NumPutativePrefixKMers);

// wait for all threads to have completed
for(ThreadIdx = 0; ThreadIdx < NumActiveThreads; ThreadIdx++)
	{
#ifdef _WIN32
	while(WAIT_TIMEOUT == WaitForSingleObject( WorkerThreads[ThreadIdx].threadHandle, 60000 * 10))
		{
		NumPutativePrefixKMers = GetKMerProcProgress(&TotSenseCnts,&TotAntisenseCnts);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress - putative prefix K-Mers: %zd",NumPutativePrefixKMers);
		}
	CloseHandle( WorkerThreads[ThreadIdx].threadHandle);
#else
	struct timespec ts;
	int JoinRlt;
	clock_gettime(CLOCK_REALTIME, &ts);
	ts.tv_sec += 60 * 10;
	while((JoinRlt = pthread_timedjoin_np(WorkerThreads[ThreadIdx].threadID, nullptr, &ts)) != 0)
		{
		NumPutativePrefixKMers = GetKMerProcProgress(&TotSenseCnts,&TotAntisenseCnts);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress - putative prefix K-Mers: %zd",NumPutativePrefixKMers);
		ts.tv_sec += 60;
		}
#endif
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Completed - putative prefix K-Mers: %zd",m_NumPutMarkers);

if(m_pSfxArray != nullptr)						// releases memory!
	{
	delete m_pSfxArray;
	m_pSfxArray = nullptr;
	}

if(gProcessingID > 0)
	{
	gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"Kmers",ePTInt64,sizeof(m_NumPutMarkers),"NumPutMarkers",&m_NumPutMarkers);
	gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"Kmers",ePTInt64,sizeof(m_TotSenseCnts),"TotSenseCnts",&m_TotSenseCnts);
	gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"Kmers",ePTInt64,sizeof(m_TotAntisenseCnts),"TotAntisenseCnts",&m_TotAntisenseCnts);
	}

if(m_NumPutMarkers == 0 || m_NumPutMarkers > 0x07fffffff)
	{
	if(m_NumPutMarkers == 0)
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Nothing to do, no putative prefix K-Mers");
	else
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Too many (%zd) putative prefix K-Mers to process",m_NumPutMarkers);
	if(m_hOutFile != -1)
		{
		if(m_MarkerBuffOfs)
			CUtility::RetryWrites(m_hOutFile,m_pMarkerBuff,m_MarkerBuffOfs);
#ifdef _WIN32
		_commit(m_hOutFile);
#else
		fsync(m_hOutFile);
#endif
		close(m_hOutFile);
		m_hOutFile = -1;
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Closed output multi-fasta marker file '%s'",m_szMarkerFile);
		}
	Reset();
	return(m_NumPutMarkers == 0 ? eBSFSuccess : eBSFerrMaxEntries);
	}

// at least one putative marker sequences identified
// allocate and initialise for index over putative  markers
m_AllocPutMarkersIndexSize = (size_t)m_NumPutMarkers * sizeof(uint32_t);
#ifdef _WIN32
m_pPutMarkersIndex = (uint32_t *) malloc(m_AllocPutMarkersIndexSize);
if(m_pPutMarkersIndex == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to allocate %zd bytes contiguous memory for putative marker sequence index",(int64_t)m_AllocPutMarkersIndexSize);
	m_AllocPutMarkersIndexSize = 0;
	Reset(false);
	return(eBSFerrMem);
	}
#else
// gnu malloc is still in the 32bit world and seems to have issues if more than 2GB allocation
m_pPutMarkersIndex = (uint32_t *)mmap(nullptr,m_AllocPutMarkersIndexSize, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
if(m_pPutMarkersIndex == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to allocate %zd bytes contiguous memory for putative marker sequences",(int64_t)m_AllocPutMarkersIndexSize);
	m_AllocPutMarkersIndexSize = 0;
	m_pPutMarkersIndex = nullptr;
	Reset(false);
	return(eBSFerrMem);
	}
#endif

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Identifying redundant (duplicates, antisense, overlaps) putative markers ...");

uint32_t TotRedundant;
uint32_t NumOverlapping;
uint32_t NumDuplicates;
uint32_t NumAntisense;
uint32_t Idx;
uint32_t *pPutMarkerIdx;
tsPutMarker *pMarker;
tsPutMarker *pAntisense;
pPutMarkerIdx = m_pPutMarkersIndex;
for(Idx = 0; Idx < m_NumPutMarkers; Idx++, pPutMarkerIdx++)
	*pPutMarkerIdx = Idx;

TotRedundant = 0;
NumOverlapping = 0;
NumDuplicates = 0;
NumAntisense = 0;

// if more than 1 then need to sort the marker sequences
if(m_NumPutMarkers > 1)
	{

	tsPutMarker *pNxtMarker;
	CMTqsort mtqsort;
	mtqsort.SetMaxThreads(NumThreads);
	gpPutativeMarkers = m_pPutMarkers;
	gPutativeMarkerSeqLen = m_PrefixLen;
	gPutMarkerSize = m_PutMarkerSize;
	mtqsort.qsort(m_pPutMarkersIndex,m_NumPutMarkers,sizeof(uint32_t),SortPutativeSeqs);
	pPutMarkerIdx = m_pPutMarkersIndex;
	for(Idx = 0; Idx < m_NumPutMarkers; Idx++, pPutMarkerIdx++)
		{
		pMarker = (tsPutMarker *)((uint8_t *)m_pPutMarkers + (*pPutMarkerIdx * (uint64_t)m_PutMarkerSize));
		pMarker->MarkerID = Idx + 1;
		}

	// filter for duplicates, overlaps and antisense marker sequences
	pPutMarkerIdx = m_pPutMarkersIndex;
	for(Idx = 0; Idx < m_NumPutMarkers; Idx++, pPutMarkerIdx++)
		{
		pMarker = (tsPutMarker *)((uint8_t *)m_pPutMarkers + (*pPutMarkerIdx * (uint64_t)m_PutMarkerSize));
		
		// check for duplicate and mark first instance of that duplicate
		if(Idx < (m_NumPutMarkers - 1))
			{
			pNxtMarker = (tsPutMarker *)((uint8_t *)m_pPutMarkers + (pPutMarkerIdx[1] * (uint64_t)m_PutMarkerSize));;
			if(!memcmp(pMarker->MarkerSeq,pNxtMarker->MarkerSeq,m_PrefixLen))
				{
				pMarker->Flags |= cMarkerDupFlg;
				NumDuplicates += 1;
				}
			}

		// check for overlap and mark the overlapping instance
		if(IsMarkerOverlapping(pMarker))
			{
			pMarker->Flags |= cMarkerOvlFlg;
			NumOverlapping += 1;
			}

		// if antisense processing then check for other marker antisense to current and mark that other marker if current not duplicate nor overlapping
		if(m_PMode == ePMSenseAntiKMers && (pAntisense = GetMarkerAntisense(pMarker)) != nullptr)
			{
			if(pAntisense->MarkerID < pMarker->MarkerID || pMarker->Flags & cMarkerAntiFlg)	// if antisense already processed or current sense marked antisense then don't re-mark
				continue;
			NumAntisense += 1;
			if(pMarker->Flags & 0x0ff)						// mark current as antisense if already marked
				pMarker->Flags |= cMarkerAntiFlg;
			else
				pAntisense->Flags |= cMarkerAntiFlg;
			}
		}

	// sort the markers by NumCultivars and sense/antisense counts as that is the order in which markers will be reported
	mtqsort.qsort(m_pPutMarkersIndex,m_NumPutMarkers,sizeof(uint32_t),SortNumCultivarsCnts);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Identified redundant putative markers, duplicates: %u, antisense: %u, overlaps: %u",NumDuplicates,NumAntisense,NumOverlapping);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Writing non-redundant markers to file: '%s'",m_szMarkerFile);

pPutMarkerIdx = m_pPutMarkersIndex;
for(Idx = 0; Idx < m_NumPutMarkers; Idx++, pPutMarkerIdx++)
	{
	pMarker = (tsPutMarker *)((uint8_t *)m_pPutMarkers + (*pPutMarkerIdx * m_PutMarkerSize));
	if(pMarker->Flags & 0x0ff)
		continue;
	ReportMarker(pMarker);
	}

if(m_hOutFile != -1)
	{
	if(m_MarkerBuffOfs)
		CUtility::RetryWrites(m_hOutFile,m_pMarkerBuff,m_MarkerBuffOfs);
#ifdef _WIN32
	_commit(m_hOutFile);
#else
	fsync(m_hOutFile);
#endif
	close(m_hOutFile);
	m_hOutFile = -1;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Closed output multi-fasta marker file '%s'",m_szMarkerFile);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Number of reported markers: %u",m_MarkerID);


Reset();
return(eBSFSuccess);
}


// check for overlaps, by m_PrefixLen-1 bases, of one potential prefix marker onto any other potential prefix marker
bool													// true if overlapping onto at least one other prefix marker
CMarkerKMers::IsMarkerOverlapping(tsPutMarker *pMarker)	// marker to check if overlapping onto any other sequence
{
int Cmp;
tsPutMarker *pPutMarker;
etSeqBase *pEl1;
etSeqBase *pEl2;
uint8_t *pByte = (uint8_t *)m_pPutMarkers;
uint32_t TargPsn;
uint32_t Psn;
int64_t SfxHi = (int64_t)(uint64_t)m_NumPutMarkers-1;
int64_t SfxLo = 0;
pEl1 = &pMarker->MarkerSeq[1];							// checking if pMarker is overlapping with 1 base 5' overhang  
do {
	TargPsn = (uint32_t)(((uint64_t)(SfxLo + SfxHi)) / 2L);
	Psn = m_pPutMarkersIndex[TargPsn];
	pPutMarker = (tsPutMarker *)&pByte[(uint64_t)Psn * (uint32_t)m_PutMarkerSize];
	pEl2 = pPutMarker->MarkerSeq;
	Cmp = memcmp(pEl1,pEl2,m_PrefixLen-1);
	if(!Cmp)
		{
		if(pMarker->MarkerID != pPutMarker->MarkerID)	// accept overlap if not to self		
			return(true);
		// self hit so check if overlapping either bracketing K-Mers
		if(TargPsn > 0)
			{
			Psn = m_pPutMarkersIndex[TargPsn-1];
			pPutMarker = (tsPutMarker *)&pByte[(uint64_t)Psn * (uint32_t)m_PutMarkerSize];
			pEl2 = pPutMarker->MarkerSeq;
			if(!memcmp(pEl1,pEl2,m_PrefixLen-1))
				return(true);
			}
		if(TargPsn < m_NumPutMarkers-1)
			{
			Psn = m_pPutMarkersIndex[TargPsn+1];
			pPutMarker = (tsPutMarker *)&pByte[(uint64_t)Psn * (uint32_t)m_PutMarkerSize];
			pEl2 = pPutMarker->MarkerSeq;
			if(!memcmp(pEl1,pEl2,m_PrefixLen-1))
				return(true);
			}
		return(false);
		}
	if(Cmp < 1)
		SfxHi = (int64_t)(uint64_t)TargPsn - 1;
	else
		SfxLo = (int64_t)(uint64_t)TargPsn + 1;
	}
while(SfxHi >= SfxLo);
return(false);
}

tsPutMarker *											// nullptr if no others antisense, else the other antisense marker
CMarkerKMers::GetMarkerAntisense(tsPutMarker *pMarker)	// marker to check if this is antisense to any other marker sequence
{
int Cmp;
etSeqBase AntisenseSeq[cMaxKMerLen];

tsPutMarker *pPutMarker;
etSeqBase *pEl1;
etSeqBase *pEl2;
uint8_t *pByte = (uint8_t *)m_pPutMarkers;
uint32_t TargPsn;
uint32_t Psn;
int64_t SfxHi = (int64_t)(uint64_t)m_NumPutMarkers-1;
int64_t SfxLo = 0;

memmove(AntisenseSeq,&pMarker->MarkerSeq,m_PrefixLen);
CSeqTrans::ReverseComplement(m_PrefixLen,AntisenseSeq);
pEl1 = AntisenseSeq;
do {
	TargPsn = (uint32_t)(((uint64_t)(SfxLo + SfxHi)) / 2L);
	Psn = m_pPutMarkersIndex[TargPsn];
	pPutMarker = (tsPutMarker *)&pByte[(uint64_t)Psn * (uint32_t)m_PutMarkerSize];
	pEl2 = pPutMarker->MarkerSeq;
	Cmp = memcmp(pEl1,pEl2,m_PrefixLen);
	if(!Cmp)
		{
		if(pMarker->MarkerID != pPutMarker->MarkerID)		
			return(pPutMarker);
		// self hit so check if overlapping either bracketing K-Mers
		if(TargPsn > 0)
			{
			Psn = m_pPutMarkersIndex[TargPsn-1];
			pPutMarker = (tsPutMarker *)&pByte[(uint64_t)Psn * (uint32_t)m_PutMarkerSize];
			pEl2 = pPutMarker->MarkerSeq;
			if(!memcmp(pEl1,pEl2,m_PrefixLen))
				return(pPutMarker);
			}
		if(TargPsn < m_NumPutMarkers-1)
			{
			Psn = m_pPutMarkersIndex[TargPsn+1];
			pPutMarker = (tsPutMarker *)&pByte[(uint64_t)Psn * (uint32_t)m_PutMarkerSize];
			pEl2 = pPutMarker->MarkerSeq;
			if(!memcmp(pEl1,pEl2,m_PrefixLen))
				return(pPutMarker);
			}
		return(nullptr);
		}

	if(Cmp < 1)
		SfxHi = (int64_t)(uint64_t)TargPsn - 1;
	else
		SfxLo = (int64_t)(uint64_t)TargPsn + 1;
	}
while(SfxHi >= SfxLo);
return(nullptr);
}

// Thread startup
#ifdef WIN32
unsigned int __stdcall CMarkerKMers::KMerThreadStart(void *args)
{
#else
void * CMarkerKMers::KMerThreadStart(void *args)
{
#endif
tsKMerThreadPars *pArgs = (tsKMerThreadPars *)args;
pArgs->pThis->LocateSharedPrefixKMers(pArgs);
#ifdef WIN32
ExitThread(1);
#else
return nullptr;
#endif
}

// Call back to notify of a putative marker
int 
CMarkerKMers::MarkersCallback(void *pThis,tsKMerCultsCnts *pCultsCnts)
{
int Idx;
uint8_t *pBase;
tsPutMarker *pPutMarker;
CMarkerKMers *pMarkerKMers;
etSeqBase *pMarkerSeq;
etSeqBase *pSrc;
if(pThis == nullptr || pCultsCnts == nullptr)		// better safe than risking a segfault...
	return(-1);
pMarkerKMers = (CMarkerKMers *)pThis;

pMarkerKMers->EnterCritSect();

// need to realloc memory to hold more marker sequences? 
if(((pMarkerKMers->m_NumPutMarkers + 1) * pMarkerKMers->m_PutMarkerSize) > pMarkerKMers->m_AllocPutMarkersSize)
	{
	size_t ReallocSize;
	tsPutMarker *pRealloc;
	ReallocSize = pMarkerKMers->m_AllocPutMarkersSize + ((size_t)cAllocNumPutativeSeqs * (size_t)pMarkerKMers->m_PutMarkerSize);

#ifdef _WIN32
	pRealloc = (tsPutMarker *)realloc(pMarkerKMers->m_pPutMarkers,(size_t)ReallocSize);
#else
	pRealloc = (tsPutMarker *)mremap(pMarkerKMers->m_pPutMarkers,pMarkerKMers->m_AllocPutMarkersSize,(size_t)ReallocSize,MREMAP_MAYMOVE);
	if(pRealloc == MAP_FAILED)
		pRealloc = nullptr;
#endif
	if(pRealloc == nullptr)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"MarkersCallback: putative marker sequences memory re-allocation to %zd bytes - %s",(int64_t)ReallocSize,strerror(errno));
		pMarkerKMers->LeaveCritSect();
		return(eBSFerrMem);
		}
	pMarkerKMers->m_pPutMarkers = pRealloc;
	pMarkerKMers->m_AllocPutMarkersSize = ReallocSize;
	}

pBase = (uint8_t *)pMarkerKMers->m_pPutMarkers;
pPutMarker = (tsPutMarker *)(pBase + ((size_t)pMarkerKMers->m_PutMarkerSize * pMarkerKMers->m_NumPutMarkers));
pMarkerSeq = pPutMarker->MarkerSeq;
pSrc = pCultsCnts->KMerSeq;
for(Idx = 0; Idx < pMarkerKMers->m_PrefixLen; Idx++,pMarkerSeq++,pSrc++)
	*pMarkerSeq = *pSrc & 0x0f;
pPutMarker->SenseCnts = (uint32_t)(pCultsCnts->SenseCnts  > 0x0ffffffff ? 0x0ffffffff : pCultsCnts->SenseCnts);
pPutMarker->AntisenseCnts = (uint32_t)(pCultsCnts->AntisenseCnts > 0x0ffffffff ? 0x0ffffffff : pCultsCnts->AntisenseCnts);
pPutMarker->NumCultivars = pCultsCnts->NumCultivars;
pPutMarker->Flags = 0;
pPutMarker->MarkerID = ++pMarkerKMers->m_NumPutMarkers;
pMarkerKMers->m_TotSenseCnts += pCultsCnts->SenseCnts;
pMarkerKMers->m_TotAntisenseCnts += pCultsCnts->AntisenseCnts;
pMarkerKMers->LeaveCritSect();
return(0);
}

int
CMarkerKMers::LocateSharedPrefixKMers(tsKMerThreadPars *pPars)
{
m_pSfxArray->GenKMerCultsCnts(m_PMode == ePMNSenseKMers ? true : false,pPars->StartSfxIdx,pPars->EndSfxIdx,m_PrefixLen,m_SuffixLen,m_MinWithPrefix,m_MaxHomozygotic,this,MarkersCallback);
return(eBSFSuccess);
}

// SortPutativeSeqs
// Sort putative marker sequences ascending
int
CMarkerKMers::SortPutativeSeqs(const void *arg1, const void *arg2)
{
size_t Ofs1;
size_t Ofs2;
uint8_t *pByte;
etSeqBase *pEl1;
etSeqBase *pEl2;

tsPutMarker *pPM1;
tsPutMarker *pPM2;
Ofs1 = gPutMarkerSize * *(uint32_t *)arg1;
Ofs2 = gPutMarkerSize * *(uint32_t *)arg2;

pByte = (uint8_t *)gpPutativeMarkers;
pPM1 = (tsPutMarker *)(pByte + Ofs1);
pPM2 = (tsPutMarker *)(pByte + Ofs2);
pEl1 = pPM1->MarkerSeq;
pEl2 = pPM2->MarkerSeq;
return(memcmp(pEl1,pEl2,gPutativeMarkerSeqLen));
}

// sort the markers by NumCultivars and sense/antisense counts descending
int
CMarkerKMers::SortNumCultivarsCnts(const void *arg1, const void *arg2)
{
size_t Ofs1;
size_t Ofs2;
uint8_t *pByte;
uint64_t PM1Cnts;
uint64_t PM2Cnts;

tsPutMarker *pPM1;
tsPutMarker *pPM2;
Ofs1 = gPutMarkerSize * *(uint32_t *)arg1;
Ofs2 = gPutMarkerSize * *(uint32_t *)arg2;

pByte = (uint8_t *)gpPutativeMarkers;
pPM1 = (tsPutMarker *)(pByte + Ofs1);
pPM2 = (tsPutMarker *)(pByte + Ofs2);

if(pPM1->NumCultivars > pPM2->NumCultivars)
	return(-1);
if(pPM1->NumCultivars < pPM2->NumCultivars)
	return(1);
PM1Cnts = pPM1->SenseCnts + pPM1->AntisenseCnts;
PM2Cnts = pPM2->SenseCnts + pPM2->AntisenseCnts;
if(PM1Cnts > PM2Cnts)
	return(-1);
if(PM1Cnts < PM2Cnts)
	return(1);
return(0);
}