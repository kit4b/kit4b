/*
This toolkit is a source base clone of 'BioKanga' release 4.4.2 (https://github.com/csiro-crop-informatics/biokanga) and contains
significant source code changes enabling new functionality and resulting process parameterisation changes. These changes have resulted in
incompatibility with 'BioKanga'.

Because of the potential for confusion by users unaware of functionality and process parameterisation changes then the modified source base
and resultant compiled executables have been renamed to 'kit4b' - K-mer Informed Toolkit for Bioinformatics.
The renaming will force users of the 'BioKanga' toolkit to examine scripting which is dependent on existing 'BioKanga'
parameterisations so as to make appropriate changes if wishing to utilise 'kit4b' parameterisations and functionality.

'kit4b' is being released under the Opensource Software License Agreement (GPLv3)
'kit4b' is Copyright (c) 2019
Please contact Dr Stuart Stephen < stuartjs@g3web.com > if you have any questions regarding 'kit4b'.

Original 'BioKanga' copyright notice has been retained and immediately follows this notice..
*/
/*
 * CSIRO Open Source Software License Agreement (GPLv3)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENSE for the complete license information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */

// KAligner.cpp : contains the CKAligner class implementation
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

#define _DISNPS_ 1    // experimental, exploring the utility of generating DiSNPs and TriSNPs

#include "ngskit4b.h"
#include "KAligner.h"

#include "../libkit4b/bgzf.h"
#include "../libkit4b/SAMfile.h"

sNAR CKAligner::m_NARdesc[] = {{eNARUnaligned,(char *)"NA",(char *)"Not processed for alignment"},
	{eNARAccepted,(char *)"AA",(char *)"Alignment accepted"},
	{eNARNs,(char *)"EN",(char *)"Excessive indeterminate (Ns) bases"},
	{eNARNoHit,(char *)"NL",(char *)"No potential alignment loci"},
	{eNARMMDelta,(char *)"MH",(char *)"Mismatch delta (minimum Hamming) criteria not met"},
	{eNARMultiAlign,(char *)"ML",(char *)"Aligned to multiloci"},
	{eNARTrim,(char *)"ET",(char *)"Excessively end trimmed"},
	{eNARSpliceJctn,(char *)"OJ",(char *)"Aligned as orphaned splice junction"},
	{eNARmicroInDel,(char *)"OM",(char *)"Aligned as orphaned microInDel"},
	{eNARPCRdup,(char *)"DP",(char *)"Duplicate PCR"},
	{eNARNonUnique,(char *)"DS",(char *)"Duplicate read sequence"},
	{eNARChromFilt,(char *)"FC",(char *)"Aligned to filtered target sequence"},
	{eNARRegionFilt,(char *)"PR",(char *)"Aligned to a priority region"},
	{eNARPEInsertMin,(char *)"UI",(char *)"PE under minimum insert size"},
	{eNARPEInsertMax,(char *)"OI",(char *)"PE over maximum insert size"},
	{eNARPENoHit,(char *)"UP",(char *)"PE partner not aligned"},
	{eNARPEStrand,(char *)"IS",(char *)"PE partner aligned to inconsistent strand"},
	{eNARPEChrom,(char *)"IT",(char *)"PE partner aligned to different target sequence"},
	{eNARPEUnalign,(char *)"NP",(char *)"PE alignment not accepted"},
{eNARLociConstrained,(char *)"LC",(char *)"Alignment violated loci base constraints"}};

CKAligner::CKAligner(void)
{
Init();
}


CKAligner::~CKAligner(void)
{
Reset(false);
}


int
CKAligner::Align(etPMode PMode,			// processing mode
		char *pszExperimentName,		// labels this experiment
		uint32_t SampleNthRawRead,		// sample every Nth raw read (or read pair) for processing (1..10000)
		etFQMethod Quality,				// quality scoring for fastq sequence files
		bool bSOLiD,					// if true then processing in colorspace
		bool bBisulfite,				// if true then process for bisulfite methylation patterning
		etPEproc PEproc,				// paired reads alignment processing mode
		int PairMinLen,					// accept paired end alignments with apparent length of at least this (default = 100)
		int PairMaxLen,					// accept paired end alignments with apparent length of at most this (default = 300)
		bool bPairStrand,				// accept paired ends if on same strand
		bool bPEcircularised,			// experimental - true if processing for PE spaning circularised fragments
		bool bPEInsertLenDist,			// experimental - true if stats file to include PE insert length distributions for each transcript
		eALStrand AlignStrand,			// align on to watson, crick or both strands of target
		int MinChimericLen,				// minimum chimeric length as a percentage (0 to disable, otherwise 15..99) of probe sequence length: negative if chimeric diagnostics to be reported
		bool bChimericRpt,				// report chimeric trimming detail for individual reads (default is not to report)
		int microInDelLen,				// microInDel length maximum
		int SpliceJunctLen,				// maximum splice junction length when aligning RNAseq reads
		int MinSNPreads,				// must be at least this number of reads covering any loci before processing for SNPs at this loci
		double QValue,					// QValue controlling FDR (Benjamini�Hochberg) SNP prediction
		double SNPNonRefPcnt,			// only process for SNP if more/equal than this percentage number of reads are non-ref at putative SNP loci (defaults to 25.0) 
		int MarkerLen,					// marker sequences of this length
		double MarkerPolyThres,			// maximum allowed marker sequence base polymorphism independent of centroid SNP (default 0.333, range 0.0 to 0.5)
		int PCRartefactWinLen,			// if >= 0 then window size to use when attempting to reduce the number of  PCR differential amplification artefacts (reads stacking to same loci)
		etMLMode MLMode,				// multimatch loci reads processing mode
		int MaxMLmatches,				// accept at most this number of multimatched alignments for any read
		bool bClampMaxMLmatches,		// accept as if MaxMLmatches even if more than this number of MaxMLmatches multimached alignments
		bool bLocateBestMatches,		// align for best matches, not just those 1H better than any other, upto at most MaxMLmatches
		int MaxNs,					    // allow at most this number of indeterminate eBaseNs in read before deeming as nonalignable
		int MinEditDist,				// any matches must have at least this edit distance to the next best match
		int MaxSubs,					// maximum number of substitutions allowed per 100bp of actual read length
		int Trim5,						// trim this number of bases from 5' end of reads when loading the reads
		int Trim3,						// trim this number of bases from 3' end of reads when loading the reads
		int MinAcceptReadLen,					// only accepting reads for alignment if at least this length after any end trimming
		int MaxAcceptReadLen,					// only accepting reads for alignment if no longer than this length after any end trimming
		int MinFlankExacts,				// trim matched reads on 5' and 3' flanks until at least this number of exactly matching bases in flanks
		int PCRPrimerCorrect,			// initially align with MaxSubs+PCRPrimerCorrect subs allowed but then correct substitutions in 5' 12bp until overall sub rate within MaxSubs
		int MaxRptSAMSeqsThres,			// report all SAM chroms or sequences if number of reference chroms <= this limit (defaults to 10000)
		etFMode FMode,					// output format mode
		teSAMFormat SAMFormat,			// if SAM output format then could be SAM or BAM compressed dependent on the file extension used
		int SitePrefsOfs,				// offset read start sites when processing  octamer preferencing, range -100..100
		int NumThreads,					// number of worker threads to use
		char *pszTrackTitle,			// track title if output format is UCSC BED
		int NumPE1InputFiles,			// number of input PE1 or single ended file specs
		char *pszPE1InputFiles[],		// names of input files (wildcards allowed unless processing paired ends) containing raw reads
		int NumPE2InputFiles,			// number of input PE2 file specs
		char *pszPE2InputFiles[],		// optional raw paired reads are in these files
		char *pszPriorityRegionFile,	// optional high priority BED file contains prioritised region loci
		bool bFiltPriorityRegions,		// true if non-priority alignments to be filtered out 
		char *pszOutFile,				// where to write alignments
		char *pszSNPFile,				// Output SNPs (CSV or VCF format) to this file (default is to output file name with '.snp' appended)
		char *pszMarkerFile,			// Output markers to this file
		char *pszSNPCentroidFile,		// Output SNP centroids (CSV format) to this file (default is for no centroid processing)
		char *pszSfxFile,				// target as suffix array
		char *pszStatsFile,				// aligner induced substitutions stats file
		char *pszMultiAlignFile,		// file to contain reads which are aligned to multiple locations
		char *pszNoneAlignFile,			// file to contain reads which were non-alignable
		char *pszSitePrefsFile,			// file to contain aligned reads octamer preferencing
		char *pszLociConstraintsFile,	// loci base constraints file
		char *pszContamFile,			// file containing contaminant sequences
		int	NumIncludeChroms,			// number of chromosome regular expressions to include
		char **ppszIncludeChroms,		// array of include chromosome regular expressions
		int	NumExcludeChroms,			// number of chromosome expressions to exclude
		char **ppszExcludeChroms)		// array of exclude chromosome regular expressions
{
int Rslt = 0;
int BuffLen = 0;
int BuffOfs = 0;
int SeqIdx = 0;
char szPEInsertDistFile[_MAX_PATH];
char szOutBAIFile[_MAX_PATH];
Init(); 

m_bPackedBaseAlleles = (FMode == eFMPBA) ? true : false;
if(pszExperimentName == nullptr || pszExperimentName[0] == '\0')
	strcpy(m_szExperimentName,"Unspecified");
else
	strcpy(m_szExperimentName,pszExperimentName);		// labels this experiment

if(MinChimericLen > 0)					// too confusing if trimming chimeric and then PCR primer trimming or flank exact trimming. Chimeric trimming should handle both PCR and flank exacts
	PCRPrimerCorrect = 0;

m_SampleNthRawRead = SampleNthRawRead;

m_pszTrackTitle = pszTrackTitle;

m_pszOutFile = pszOutFile;

if(SAMFormat == etSAMFBAM)
	{
	strcpy(szOutBAIFile,pszOutFile);
	strcat(szOutBAIFile,".bai");
	}
else
	szOutBAIFile[0] = '\0';
m_pszOutBAIFile = szOutBAIFile;

m_pszSNPRsltsFile = pszSNPFile;
if(m_FMode != eFMbed && pszSNPFile != nullptr && pszSNPFile[0] != '\0')
	{
	// check if file extension is VCF, if so then output SNPs as VCF instead of the default CSV
	int SNPFileNameLen;
	m_bXCSVFrameShifts = false;
	SNPFileNameLen = (int)strlen(m_pszSNPRsltsFile);

	if(!stricmp(".vcf",&pszSNPFile[SNPFileNameLen-4]))
		m_bSNPsVCF = true;
	else
		if(!stricmp(".csvx", &pszSNPFile[SNPFileNameLen - 5]))	// special case - an extension of 'csvx' flags that generated SNPs csv is to be in extended format containing codon frameshifts
			{
			m_bXCSVFrameShifts = true;
			pszSNPFile[SNPFileNameLen-1] = '\0';
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "Generating extended SNPs CSV containing codon frame shifts into '%s'", pszSNPFile);
			}
	}
m_pszSNPCentroidFile = pszSNPCentroidFile;
m_NumPE2InputFiles = NumPE2InputFiles;

m_ppszPE2InputFiles = pszPE2InputFiles;
m_pszMultiAlignFile = pszMultiAlignFile;
m_pszNoneAlignFile = pszNoneAlignFile;
m_pszStatsFile = pszStatsFile;
m_pszSfxFile = pszSfxFile;

if(bPEInsertLenDist && m_pszStatsFile != nullptr && m_pszStatsFile[0] != '\0')
	{
	strcpy(szPEInsertDistFile,pszStatsFile);
	strcat(szPEInsertDistFile,".peinserts.csv");
	m_pszPEInsertDistFile = szPEInsertDistFile;
	}
else
	{
	bPEInsertLenDist = false;
	m_pszPEInsertDistFile = nullptr;
	}

m_pszSitePrefsFile = pszSitePrefsFile;

m_NumPE1InputFiles = NumPE1InputFiles;
m_ppszPE1InputFiles = pszPE1InputFiles;
m_bFiltPriorityRegions = bFiltPriorityRegions;

m_bPEcircularised = bPEcircularised;		// experimental - true if processing for PE spaning circularised fragments
m_bPEInsertLenDist = bPEInsertLenDist;		// true if stats file to include PE insert length distributions for each transcript

m_PMode = PMode;
m_FMode = FMode;
m_SAMFormat = SAMFormat;
m_PEproc = PEproc;
m_QMethod = Quality;
m_NumThreads = NumThreads;
m_bBisulfite = bBisulfite;
m_MaxMLmatches = MaxMLmatches;
m_bClampMaxMLmatches = bClampMaxMLmatches;
m_bLocateBestMatches = bLocateBestMatches; 
m_MLMode = MLMode;
m_MaxNs = MaxNs;
m_MaxSubs = MaxSubs;
m_MinEditDist = MinEditDist;

m_PairMinLen = PairMinLen;
m_PairMaxLen = PairMaxLen;
m_bPairStrand = bPairStrand;

if(PCRPrimerCorrect > 0)
	m_InitalAlignSubs = min(MaxSubs + PCRPrimerCorrect,cMaxAllowedSubs);
else
	m_InitalAlignSubs = MaxSubs;

m_AlignStrand = AlignStrand;
m_MinChimericLen = abs(MinChimericLen);
m_bReportChimerics = bChimericRpt && MLMode == eMLdefault ? true : false;
m_microInDelLen = microInDelLen;
m_SpliceJunctLen = SpliceJunctLen;
m_MinSNPreads = MinSNPreads;
m_SNPNonRefPcnt = SNPNonRefPcnt/100.0;
m_MinCoverageSegments = min(m_MinSNPreads,cMinCoverageSegments);

m_QValue = QValue;
m_Marker5Len = MarkerLen/2;
m_Marker3Len = MarkerLen - 1 - m_Marker5Len;
m_MarkerPolyThres = MarkerPolyThres;
m_pszMarkerFile = pszMarkerFile;

m_Trim5 = Trim5;						// trim this number of bases from 5' end of reads when loading the reads
m_Trim3 = Trim3;						// trim this number of bases from 3' end of reads when loading the reads

m_MinAcceptReadLen = MinAcceptReadLen;	// only accepting reads for alignment if at least this length after any end trimming
m_MaxAcceptReadLen = MaxAcceptReadLen;	// only accepting reads for alignment if no longer than this length after any end trimming

m_SitePrefsOfs = SitePrefsOfs;			// offset read start sites when processing  octamer preferencing, range -100..100

m_MaxRptSAMSeqsThres = MaxRptSAMSeqsThres;

if(CreateMutexes()!=eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Failed to create thread synchronisation mutexes");
	Reset(false);
	return(cBSFSyncObjErr);
	}

m_mtqsort.SetMaxThreads(NumThreads);

// load contaminants if user has specified a contaminant sequence file
if(pszContamFile != nullptr && pszContamFile[0] != '\0')
	{
	if((m_pContaminants = new CContaminants)==nullptr)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CContaminants");
		Reset(false);
		return(eBSFerrObj);
		}
	if((Rslt = m_pContaminants->LoadContaminantsFile(pszContamFile)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to load contaminate sequences file '%s'",pszContamFile);
		Reset(false);
		return(Rslt);
		}
	if(Rslt == 0)
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"No contaminant sequences loaded from '%s'",pszContamFile);
	}

// compile include/exclude chromosome regular expressions
if(NumIncludeChroms > 0 || NumExcludeChroms > 0)
	if((Rslt = CompileChromRegExprs(NumIncludeChroms,ppszIncludeChroms,NumExcludeChroms,ppszExcludeChroms)) != eBSFSuccess)
		{
		Reset(false);
		return(Rslt);
		}

// if specified then load high priority regions from BED file
m_pPriorityRegionBED = nullptr;
if(pszPriorityRegionFile != nullptr && pszPriorityRegionFile[0] != '\0')
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading high priority regions BED file '%s'", pszPriorityRegionFile);
	if((m_pPriorityRegionBED = new CBEDfile()) == nullptr)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBEDfile");
		Reset(false);
		return(eBSFerrObj);
		}
	if((Rslt=m_pPriorityRegionBED->Open(pszPriorityRegionFile))!=eBSFSuccess)
		{
		while(m_pPriorityRegionBED->NumErrMsgs())
			gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pPriorityRegionBED->GetErrMsg());
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open high priority regions BED file '%s'",pszPriorityRegionFile);
		Reset(false);
		return(Rslt);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"High priority regions BED file '%s' loaded", pszPriorityRegionFile);
	}


// open bioseq file containing suffix array for targeted assembly to align reads against
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading suffix array file '%s'", pszSfxFile);
if((m_pSfxArray = new CSfxArray()) == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CSfxArray");
	Reset(false);
	return(eBSFerrObj);
	}
if((Rslt=m_pSfxArray->Open(pszSfxFile,false,bBisulfite,bSOLiD))!=eBSFSuccess)
	{
	while(m_pSfxArray->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pSfxArray->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open input bioseq suffix array file '%s'",pszSfxFile);
	Reset(false);
	return(Rslt);
	}

// report to user some sfx array metadata as conformation the targeted assembly is correct
strcpy(m_szTargSpecies,m_pSfxArray->GetDatasetName());
m_bIsSOLiD = m_pSfxArray->IsSOLiD();
tsSfxHeaderV3 SfxHeader;
m_pSfxArray->GetSfxHeader(&SfxHeader);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Genome Assembly Name: '%s' Descr: '%s' Title: '%s' Version: %d",
					 m_szTargSpecies,SfxHeader.szDescription,SfxHeader.szTitle,SfxHeader.Version);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Assembly has blocks: %d, max block size: %zu",SfxHeader.NumSfxBlocks,SfxHeader.SfxBlockSize);

// if user has specified that there are constraints on loci bases then need to load the loci constraints file
if(pszLociConstraintsFile != nullptr && pszLociConstraintsFile[0] != '\0')
	{
	if((Rslt=LoadLociConstraints(pszLociConstraintsFile)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to load loci constraints");
		Reset(false);
		return(Rslt);
		}
	}


// restrict the max core iterations according to the requested sensitivity
int MaxIter;
switch(PMode) {
	case ePMdefault:			// default processing mode
		MaxIter = cDfltKASensCoreIters;
		break;
	case ePMMoreSens:			// more sensitive - slower
		MaxIter = cMoreKASensCoreIters;
		break;
	case ePMUltraSens:			// ultra sensitive - much slower
		MaxIter = cUltraKASensCoreIters;
		break;
	case ePMLessSens:			// less sensitive - quicker
	default:
		MaxIter = cMinKASensCoreIters;
	}
m_pSfxArray->SetMaxIter(MaxIter);

// reads are loaded asynchronously to the alignment processing
if((Rslt=InitiateLoadingReads()) < eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Failed to load reads");
	Reset(false);
	return(Rslt);
	}

// if user wants to utilise reads aligning to multiple loci then need to make an initial alloc for holding these..
if(m_MLMode > eMLrand && m_MLMode < eMLall)
	{
	size_t memreq = (size_t)sizeof(tsReadHit) * cAllocMultihits;

#ifdef _WIN32
	m_pMultiHits = (tsReadHit *) malloc(memreq);	// initial and perhaps the only allocation

	if(m_pMultiHits == nullptr)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Memory allocation of %zd bytes - %s",(int64_t)memreq,strerror(errno));
		Reset(false);
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	m_pMultiHits = (tsReadHit *)mmap(nullptr,memreq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pMultiHits == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Memory allocation of %zd bytes through mmap()  failed - %s",(int64_t)memreq,strerror(errno));
		m_pMultiHits = nullptr;
		Reset(false);
		return(eBSFerrMem);
		}
#endif
	m_AllocdMultiHits = cAllocMultihits;
	m_AllocdMultiHitsMem = memreq;
	m_NumMultiHits = 0;
	m_NumUniqueMultiHits = 0;
	m_NumProvMultiAligned = 0;
	}

// if outputting multiloci all then need to allocate memory for these
if(m_MLMode >= eMLall || m_FMode == eFMsamAll)
	{
	size_t memreq = (size_t)(sizeof(tsReadHit) + 150) * 5 * cAllocMultihits;	// read sizes not known yet so assume 100bp reads plus long descriptors plus many multiloci - realloc'd as may be required

#ifdef _WIN32
	m_pMultiAll = (tsReadHit *) malloc(memreq);	// initial and perhaps the only allocation

	if(m_pMultiAll == nullptr)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Memory allocation of %zd bytes for multiAll - %s",(int64_t)memreq,strerror(errno));
		Reset(false);
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	m_pMultiAll = (tsReadHit *)mmap(nullptr,memreq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pMultiAll == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Memory allocation of %zd bytes for multiAll through mmap()  failed - %s",(int64_t)memreq,strerror(errno));
		m_pMultiAll = nullptr;
		Reset(false);
		return(eBSFerrMem);
		}
#endif
	m_AllocMultiAllMem = memreq;
	m_NumMultiAll = 0;
	m_NxtMultiAllOfs = 0;
	}

// reasonably confident that there will be some results to report so create/trunc all required result files
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Creating/truncating result files..");
if((Rslt=CreateOrTruncResultFiles())!=eBSFSuccess)
	{
	// problems.. need to ensure all background threads (at this stage should only be the reads loading or sfx loading thread) are cleanly terminated
	m_TermBackgoundThreads = 1;	// need to immediately self-terminate?
	m_pSfxArray->Reset(false);
#ifdef _WIN32
	if(m_hThreadLoadReads != nullptr)
		{
		while(WAIT_TIMEOUT == WaitForSingleObject(m_hThreadLoadReads, 5000))
			{
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: waiting for reads load thread to terminate");
			}
		CloseHandle(m_hThreadLoadReads);
		m_hThreadLoadReads = nullptr;
		}
#else
	if(m_ThreadLoadReadsID != 0)
		{
		struct timespec ts;
		int JoinRlt;
		clock_gettime(CLOCK_REALTIME, &ts);
		ts.tv_sec += 5;
		while((JoinRlt = pthread_timedjoin_np(m_ThreadLoadReadsID, nullptr, &ts)) != 0)
			{
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: waiting for reads load thread to terminate");
			ts.tv_sec += 60;
			}
		m_ThreadLoadReadsID = 0;
		}
#endif
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: reads load thread terminated");
	Reset(false);
	return(Rslt);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Creating/truncating result files completed");
m_CurReadsSortMode = eRSMReadID;			// reads were loaded and assigned ascending read identifiers so that is their initial sort order

// locate all read matches
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Aligning in %s...",bSOLiD ? "colorspace" : "basespace");
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Aligning for %s cored matches...",bBisulfite ? "bisulfite" : "normal");

Rslt = LocateCoredApprox(MinEditDist,m_InitalAlignSubs);

if(Rslt < eBSFSuccess)
	{
	Reset(false);
	return(Rslt);
	}
	
if(m_bReportChimerics)
	{
	char szChimericsFile[_MAX_PATH];
	strcpy(szChimericsFile,m_pszOutFile);
	strcat(szChimericsFile,".chimericseqs.csv");
	ReportChimerics(szChimericsFile);
	}

// if autodetermining max subs that were allowed from actual reads then let user know what the average read length was
size_t TotReadsLen = 0;
tsReadHit *pReadHit;
int AvReadsLen;
int MinReadsLen = -1;
int MaxReadsLen = 0;
pReadHit = m_pReadHits;
for(uint32_t RIdx = 0; RIdx < m_NumReadsLoaded; RIdx++)
	{
	TotReadsLen += pReadHit->ReadLen;
	if(MinReadsLen > pReadHit->ReadLen || MinReadsLen == -1)
		MinReadsLen = pReadHit->ReadLen;
	if(MaxReadsLen < pReadHit->ReadLen)
		MaxReadsLen = pReadHit->ReadLen;
	pReadHit = (tsReadHit *)((uint8_t *)pReadHit + sizeof(tsReadHit) + pReadHit->ReadLen + pReadHit->DescrLen);
	}
AvReadsLen = (int)(TotReadsLen/m_NumReadsLoaded);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Average length of all reads was: %d (min: %d, max: %d)",AvReadsLen,MinReadsLen,MaxReadsLen);
m_MaxReadsLen = MaxReadsLen;
m_MinReadsLen = MinReadsLen;
m_AvReadsLen = AvReadsLen;
if(gProcessingID > 0)
	{
	gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"ReadLen",ePTInt32,sizeof(AvReadsLen),"MeanLen",&AvReadsLen);
	gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"ReadLen",ePTInt32,sizeof(MinReadsLen),"MinLen",&MinReadsLen);
	gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"ReadLen",ePTInt32,sizeof(MaxReadsLen),"MaxLen",&MaxReadsLen);
	}

if(m_InitalAlignSubs != 0)
	{
	int MeanSubs;
	int MinSubs;
	int MaxSubs;
	MeanSubs = max(1,(AvReadsLen * m_InitalAlignSubs)/100);
	MinSubs = max(1,(MinReadsLen * m_InitalAlignSubs)/100);
	MaxSubs = max(1,(MaxReadsLen * m_InitalAlignSubs)/100);

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Typical allowed aligner induced substitutions was: %d (min: %d, max: %d)", MeanSubs, MinSubs, MaxSubs);
	if(gProcessingID > 0)
		{
		gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"AllowedSubs",ePTInt32,sizeof(MeanSubs),"MeanSubs",&MeanSubs);
		gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"AllowedSubs",ePTInt32,sizeof(MinSubs),"MinSubs",&MinSubs);
		gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"AllowedSubs",ePTInt32,sizeof(MaxSubs),"MaxSubs",&MaxSubs);
		}
	}
else
	if(gProcessingID > 0)
		gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"AllowedSubs",ePTInt32,sizeof(MaxSubs),"MeanSubs",&m_InitalAlignSubs);
 
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Provisionally accepted %d aligned reads (%d uniquely, %d aligning to multiloci) aligning to a total of %d loci", m_TotAcceptedAsAligned,m_TotAcceptedAsUniqueAligned,m_TotAcceptedAsMultiAligned,m_TotLociAligned);

m_OrigNumReadsLoaded = m_NumReadsLoaded;		// make a copy of actual number of reads loaded as if m_MLMode >= eMLall then will be overwritten with number of multialigned loci 
if(m_MLMode >= eMLall)		// a little involved as need to reuse m_pReadHits ptrs so sorting and reporting of multi-SAM hits will be same as if normal processing...
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Treating accepted %d multialigned reads as uniquely aligned %d source reads in subsequent processing",m_TotAcceptedAsMultiAligned,m_TotLociAligned - m_TotAcceptedAsUniqueAligned);
	if(m_pReadHits != nullptr)
		{
#ifdef _WIN32
		free(m_pReadHits);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
		if(m_pReadHits != MAP_FAILED)
			munmap(m_pReadHits,m_AllocdReadHitsMem);
#endif
		m_pReadHits = nullptr;
		}
	m_pReadHits = m_pMultiAll;
	m_pMultiAll = nullptr;
	m_AllocdReadHitsMem = m_AllocMultiAllMem;
	m_UsedReadHitsMem = m_NxtMultiAllOfs;
	m_NumReadsLoaded = m_NumMultiAll;
	m_FinalReadID = m_NumMultiAll;
	m_AllocMultiAllMem = 0;
	m_NxtMultiAllOfs = 0;
	m_NumMultiAll = 0;
	if(m_ppReadHitsIdx != nullptr)
		{
		delete []m_ppReadHitsIdx;
		m_ppReadHitsIdx = nullptr;
		m_AllocdReadHitsIdx = 0;
		}
	m_ppReadHitsIdx = 0;
	}

// if PE processing then try assign partners within insert size constraints
if(PEproc != ePEdefault)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Paired end association and partner alignment processing started..");
	if((Rslt=ProcessPairedEnds(PEproc,MinEditDist,PairMinLen,PairMaxLen,bPairStrand,m_InitalAlignSubs)) < eBSFSuccess)
		{
		Reset(false);
		return(Rslt);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Paired end association and partner alignment processing completed..");
	}

// try and assign multimatch read loci?
// only applies if SE processing and non-random assignment of a single multiloci loci
if(PEproc == ePEdefault && m_MLMode > eMLrand && m_MLMode != eMLall)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Multialignment processing started..");
	if((Rslt = AssignMultiMatches()) < eBSFSuccess)
		{
		Reset(false);
		return(Rslt);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Multialignment processing completed");
	}

IdentifyConstraintViolations(PEproc != ePEdefault);

// if requested then attempt to reduce the number of  PCR differential amplification artefacts (reads stacking to same loci)
if(PEproc == ePEdefault && PCRartefactWinLen >= 0)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing to reduce PCR differential amplification artefacts processing started..");
	if((Rslt=ReducePCRduplicates(PCRartefactWinLen)) < eBSFSuccess)
		{
		Reset(false);
		return(Rslt);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"PCR differential amplification artefacts processing completed");
	}

if(PCRPrimerCorrect > 0) 
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"PCR 5' Primer correction processing started..");
	if((Rslt=PCR5PrimerCorrect(m_MaxSubs)) < eBSFSuccess)
		{
		Reset(false);
		return(Rslt);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"PCR 5' Primer correction processing completed");
	}

if(MinFlankExacts > 0) // autotrim aligned reads flanks?
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Autotrim aligned read flank processing started..");
	if((Rslt=AutoTrimFlanks(MinFlankExacts)) < eBSFSuccess)
		{
		Reset(false);
		return(Rslt);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Autotrim aligned read flank processing completed");
	}

// if splice junctions being processed then check for orphans and remove these
if(PEproc == ePEdefault && SpliceJunctLen > 0)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Removal of orphan splice junction processing started..");
	if((Rslt=RemoveOrphanSpliceJuncts(SpliceJunctLen)) < eBSFSuccess)
		{
		Reset(false);
		return(Rslt);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Removal of orphan splice junction processing completed");
	}

// if processing for microInDels then need to check for orphans and remove these
if(PEproc == ePEdefault && microInDelLen > 0)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Removal of orphan microInDels processing started..");
	if((Rslt=RemoveOrphanMicroInDels(microInDelLen)) < eBSFSuccess)
		{
		Reset(false);
		return(Rslt);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Removal of orphan microInDels processing completed");
	}

BuffLen = 0;
BuffOfs = 0;
SeqIdx;

// now apply any chromosome filtering that user has specified
if(NumExcludeChroms || NumIncludeChroms)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Filtering aligned reads by chromosome started..");
	if((Rslt=FiltByChroms())<eBSFSuccess)
		{
		Reset(false);
		return(Rslt);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Filtering aligned reads by chromosome completed");
	}

// apply priority region filtering if requested
if(m_pPriorityRegionBED != nullptr && m_bFiltPriorityRegions)
	FiltByPriorityRegions();

// user interested in the nonaligned?
if(m_hNoneAlignFile != -1 || m_gzNoneAlignFile != nullptr)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Reporting of non-aligned reads started..");
	if((Rslt=ReportNoneAligned())<eBSFSuccess)
		{
		Reset(false);
		return(Rslt);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Reporting of non-aligned reads completed");
	}

// user interested in the multialigned?
// these only include those reads which would otherwise have been accepted but aligned to multiple loci
if(m_hMultiAlignFile != -1 || m_gzMultiAlignFile != nullptr)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Reporting of multialigned reads started..");
	if((Rslt=ReportMultiAlign()) < eBSFSuccess)
		{
		Reset(false);
		return(Rslt);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Reporting of multialigned reads completed");
	}

// all processing to accept alignments completed, can now report basic stats
if((Rslt=ReportAlignStats()) < eBSFSuccess)
	{
	Reset(false);
	return(Rslt);
	}

Rslt = eBSFSuccess;
if(!m_bPackedBaseAlleles)
	{
	if(m_NARAccepted && m_hSitePrefsFile != -1)
		ProcessSiteProbabilites(m_SitePrefsOfs);

	// now time to write out the read hits
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Reporting of aligned result set started...");
	if(FMode >= eFMsam)
		{
		if(m_hJctOutFile != -1 || m_hIndOutFile != -1)	// even though SAM for read alignments, splice and indels are reported as BED format
			{
			Rslt = WriteReadHits(PEproc == ePEdefault ? false : true);
			if(Rslt < eBSFSuccess)
				{
				Reset(false);
				return(Rslt);
				}
			}
		Rslt = WriteBAMReadHits(FMode,SAMFormat,PEproc == ePEdefault ? false : true,6);	// default to compression level 6
		}
	else
		Rslt = WriteReadHits(PEproc == ePEdefault ? false : true);

	if(Rslt < eBSFSuccess)
		{
		Reset(false);
		return(Rslt);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Reporting of aligned result set completed");

	if(m_bPEInsertLenDist && m_NARAccepted)
		ReportPEInsertLenDist();

	if(Rslt >= eBSFSuccess && m_NARAccepted && m_hStatsFile > 0 && m_MaxAlignLen > 0)
		Rslt = WriteBasicCountStats();

	if(Rslt >= eBSFSuccess && m_NARAccepted && m_hStatsFile > 0 && m_MaxAlignLen > 0)
		Rslt = ReportTargHitCnts();

	if(Rslt >= eBSFSuccess  && m_NARAccepted && m_hSitePrefsFile > 0)
		Rslt = WriteSitePrefs();

	m_TotNumSNPs = 0;
	if(Rslt >= eBSFSuccess && m_NARAccepted && MinSNPreads > 0 && m_hSNPfile != -1)
		{
		bool bMarkers;
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing for SNPs and writing out SNPs to file '%s",m_pszSNPRsltsFile);
		if(m_hMarkerFile != -1)
			{
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing for Markers and writing out marker sequences to file '%s",m_pszMarkerFile);
			bMarkers = true;
			}
		else
			bMarkers = false;
		Rslt = ProcessSNPs();			// track title if output format is to be UCSC BED, will have '_SNPs' appended
		if(Rslt >= eBSFSuccess)
			{
			if(bMarkers)
				gDiagnostics.DiagOut(eDLInfo,gszProcName,"Marker processing completed with %d marker sequences writtten to file '%s",m_MarkerID,m_pszMarkerFile);
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"SNP processing completed with %d putative SNPs discovered",m_TotNumSNPs);
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"There are %zd aligned loci bases which are covered by %zd read bases with mean coverage of %1.2f",m_LociBasesCovered,m_LociBasesCoverage,m_LociBasesCoverage/(double)m_LociBasesCovered);
			}
		}

	if(gProcessingID != 0)
		gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"SNPs",ePTInt32,sizeof(AvReadsLen),"Cnt",&m_TotNumSNPs);
	}
else
	Rslt = ProcessSNPs();
Reset(Rslt >= eBSFSuccess ? true : false);
return(Rslt);
}

void
CKAligner::Init(void)
{
m_hInFile = -1;
m_hOutFile = -1;

m_hBAIFile = -1;

m_hIndOutFile = -1;
m_hJctOutFile = -1;
m_hStatsFile = -1;
m_hInsertLensFile = -1;
m_hSitePrefsFile = -1;
m_hNoneAlignFile = -1;
m_hMultiAlignFile = -1;
m_hSNPfile = -1;
m_hDiSNPfile = -1;
m_hTriSNPfile = -1;
m_hMarkerFile = -1;	
m_hSNPCentsfile = -1;
m_hWIGSpansFile = -1;
m_hPackedBaseAllelesFile = -1;
m_bPackedBaseAlleles = false;
m_gzOutFile = nullptr;
m_gzIndOutFile = nullptr;
m_gzJctOutFile = nullptr;
m_gzNoneAlignFile = nullptr;
m_gzMultiAlignFile = nullptr;
m_gzSNPfile = nullptr;
m_gzSNPCentsfile = nullptr;

m_bgzOutFile = false;
m_bgzNoneAlignFile = false;
m_bgzMultiAlignFile = false;
m_AllocdCovSegBuff = 0;
m_CovSegBuffIdx = 0;
m_pszCovSegBuff = nullptr;
m_pReadHits = nullptr;
m_ppReadHitsIdx = nullptr;
m_pMultiHits = nullptr;
m_pMultiAll = nullptr;
m_pLociPValues = nullptr;
m_pPackedBaseAlleles = nullptr;
m_pSfxArray = nullptr;
m_pPriorityRegionBED = nullptr;
m_pAllocsIdentNodes = nullptr;
m_pAllocsMultiHitLoci = nullptr;
m_pAllocsMultiHitBuff = nullptr;
m_pChromSNPs = nullptr;
m_pszLineBuff = nullptr;
m_pLenDist = nullptr;
m_pSNPCentroids = nullptr;
m_pConstraintLoci = nullptr; 
m_pContaminants = nullptr;
m_WIGChromID = 0;
m_WIGSpanLoci = 0;
m_WIGSpanLen = 0;
m_WIGSpanCnts = 0;
m_NumConstraintLoci = 0;
m_NumConstrainedChroms = 0;
m_ConstrainedChromIDs[0] = 0;
m_bSNPsVCF = false;
m_LociBasesCovered = 0;
m_LociBasesCoverage = 0;
m_PrevSizeOf = 0;
m_NumMultiAll = 0;
m_NxtMultiAllOfs = 0;
m_AllocMultiAllMem = 0;
m_szLineBuffIdx = 0;
m_MaxMLmatches = 0;
m_bClampMaxMLmatches = false;
m_bLocateBestMatches = false;
m_TotAllocdIdentNodes = 0;
m_PerThreadAllocdIdentNodes = 0;
m_AllocdReadHitsMem = 0;
m_UsedReadHitsMem = 0;
m_NumReadsLoaded = 0;
m_UsedReadHitsMem = 0;
m_FinalReadID = 0;
m_NumDescrReads = 0;
m_FinalReadID = 0;
m_AllocdReadHitsIdx = 0;
m_AllocdMultiHits = 0;
m_AllocdMultiHitsMem = 0;
m_NumMultiHits = 0;
m_NumUniqueMultiHits = 0;
m_NumProvMultiAligned = 0;
m_MaxAlignLen = 0;
m_MaxNs = 0;
m_MaxRptSAMSeqsThres = 0;
m_bPEcircularised = false;
m_bPEInsertLenDist = false;	
m_bIsSOLiD = false;
m_bBisulfite = false;
m_AlignStrand = eALSnone;
m_MinChimericLen = 0;
m_microInDelLen = 0;
m_SpliceJunctLen = 0;
m_AllocLociPValuesMem = 0;
m_NumLociPValues = 0;
m_QValue = 0.0;
m_MinSNPreads = 0;
m_SNPNonRefPcnt = 0.0; 
m_MaxDiSNPSep = cDfltMaxDiSNPSep;
m_MarkerID = 0;	
m_Marker5Len = 0;
m_Marker3Len = 0;
m_MarkerPolyThres = 0;
m_NumReadsProc = 0;
m_NxtReadProcOfs = 0;
m_ElimPlusTrimed = 0;
m_ElimMinusTrimed = 0;
m_PEproc = ePEdefault;
m_TotNonAligned = 0;
m_NumSloughedNs = 0;
m_TotAcceptedAsAligned = 0;
m_TotLociAligned = 0;
m_TotAcceptedAsUniqueAligned = 0;
m_TotAcceptedAsMultiAligned = 0;
m_TotNotAcceptedDelta = 0;
m_TotAcceptedHitInsts = 0;
m_SitePrefsOfs = cDfltRelSiteStartOfs;
m_pszOutFile = nullptr;
m_szIndRsltsFile[0] = '\0';
m_szJctRsltsFile[0] = '\0';
m_szCoverageSegmentsFile[0] = '\0';
m_szDiSNPFile[0] = '\0';
m_pszSNPRsltsFile = nullptr;
m_pszSNPCentroidFile = nullptr;
m_pszSitePrefsFile = nullptr;
m_MaxReadsLen = 0;
m_MinReadsLen = 0;
m_AvReadsLen = 0;
m_NARAccepted = 0;
m_TermBackgoundThreads = 0;
m_CurClusterFrom = 0;
m_bFiltPriorityRegions = false;
m_SAMFormat = etSAMFformat;
m_CurReadsSortMode = eRSMunsorted;
m_ThreadCoredApproxRslt = 0;
m_bXCSVFrameShifts = false;
m_AllocPackedBaseAllelesMem = 0;
m_NumPackedBaseAlleles = 0;

#ifdef _WIN32
m_hThreadLoadReads = nullptr;
#else
m_ThreadLoadReadsID = 0;
#endif
memset(m_AlignQSubDist,0,sizeof(m_AlignQSubDist));
memset(m_AlignMSubDist,0,sizeof(m_AlignMSubDist));
memset(m_MultiHitDist,0,sizeof(m_MultiHitDist));
memset(&m_FileHdr,0,sizeof(m_FileHdr));
m_bMutexesCreated = false;
}

void
CKAligner::Reset(bool bSync)			// if bSync true then fsync before closing output file handles
{
m_TermBackgoundThreads = 0x01;	// need to require any background threads to self-terminate
if(m_hInFile != -1)
	{
	close(m_hInFile);
	m_hInFile = -1;
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

if(m_hBAIFile != -1)
	{
	if(bSync)
#ifdef _WIN32
		_commit(m_hBAIFile);
#else
		fsync(m_hBAIFile);
#endif
	close(m_hBAIFile);
	m_hBAIFile = -1;
	}

if(m_hJctOutFile != -1)
	{
	if(bSync)
#ifdef _WIN32
		_commit(m_hJctOutFile);
#else
		fsync(m_hJctOutFile);
#endif
	close(m_hJctOutFile);
	m_hJctOutFile = -1;
	}
if(m_hIndOutFile != -1)
	{
	if(bSync)
#ifdef _WIN32
		_commit(m_hIndOutFile);
#else
		fsync(m_hIndOutFile);
#endif
	close(m_hIndOutFile);
	m_hIndOutFile = -1;
	}

if(m_hStatsFile != -1)
	{
	if(bSync)
#ifdef _WIN32
		_commit(m_hStatsFile);
#else
		fsync(m_hStatsFile);
#endif
	close(m_hStatsFile);
	m_hStatsFile = -1;
	}

if(m_hInsertLensFile != -1)
	{
	if(bSync)
#ifdef _WIN32
		_commit(m_hInsertLensFile);
#else
		fsync(m_hInsertLensFile);
#endif
	close(m_hInsertLensFile);
	m_hInsertLensFile = -1;
	}

if(m_hNoneAlignFile != -1)
	{
	if(bSync)
#ifdef _WIN32
		_commit(m_hNoneAlignFile);
#else
		fsync(m_hNoneAlignFile);
#endif
	close(m_hNoneAlignFile);
	m_hNoneAlignFile = -1;
	}
if(m_hMultiAlignFile != -1)
	{
	if(bSync)
#ifdef _WIN32
		_commit(m_hMultiAlignFile);
#else
		fsync(m_hMultiAlignFile);
#endif
	close(m_hMultiAlignFile);
	m_hMultiAlignFile = -1;
	}
if(m_hSitePrefsFile != -1)
	{
	if(bSync)
#ifdef _WIN32
		_commit(m_hSitePrefsFile);
#else
		fsync(m_hSitePrefsFile);
#endif
	close(m_hSitePrefsFile);
	m_hSitePrefsFile = -1;
	}

if(m_hSNPfile != -1)
	{
	if(bSync)
#ifdef _WIN32
		_commit(m_hSNPfile);
#else
		fsync(m_hSNPfile);
#endif
	close(m_hSNPfile);
	m_hSNPfile = -1;
	}

if(m_hDiSNPfile != -1)
	{
	if(bSync)
#ifdef _WIN32
		_commit(m_hDiSNPfile);
#else
		fsync(m_hDiSNPfile);
#endif
	close(m_hDiSNPfile);
	m_hDiSNPfile = -1;
	}
if(m_hTriSNPfile != -1)
	{
	if(bSync)
#ifdef _WIN32
		_commit(m_hTriSNPfile);
#else
		fsync(m_hTriSNPfile);
#endif
	close(m_hTriSNPfile);
	m_hTriSNPfile = -1;
	}

if(m_hWIGSpansFile != -1)
	{
	if(bSync)
#ifdef _WIN32
		_commit(m_hWIGSpansFile);
#else
		fsync(m_hWIGSpansFile);
#endif
	close(m_hWIGSpansFile);
	m_hWIGSpansFile = -1;
	}

if(m_hPackedBaseAllelesFile != -1)
	{
	if(bSync)
#ifdef _WIN32
		_commit(m_hPackedBaseAllelesFile);
#else
		fsync(m_hPackedBaseAllelesFile);
#endif
	close(m_hPackedBaseAllelesFile);
	m_hPackedBaseAllelesFile = -1;
	}

if(m_hMarkerFile != -1)
	{
	if(bSync)
#ifdef _WIN32
		_commit(m_hMarkerFile);
#else
		fsync(m_hMarkerFile);
#endif
	close(m_hMarkerFile);
	m_hMarkerFile = -1;
	}

if(m_hSNPCentsfile != -1)
	{
	if(bSync)
#ifdef _WIN32
		_commit(m_hSNPCentsfile);
#else
		fsync(m_hSNPCentsfile);
#endif
	close(m_hSNPCentsfile);
	m_hSNPCentsfile = -1;
	}

if(m_pSNPCentroids != nullptr)
	{
	delete []m_pSNPCentroids;
	m_pSNPCentroids = nullptr;
	}


if(m_pConstraintLoci != nullptr)
	{
	delete []m_pConstraintLoci;
	m_pConstraintLoci = nullptr;
	}

if(m_gzOutFile != nullptr)
	{
	gzclose(m_gzOutFile);
	m_gzOutFile = nullptr;
	}

if(m_gzIndOutFile != nullptr)
	{
	gzclose(m_gzIndOutFile);
	m_gzIndOutFile = nullptr;
	}

if(m_gzJctOutFile != nullptr)
	{
	gzclose(m_gzJctOutFile);
	m_gzJctOutFile = nullptr;
	}

if(m_gzNoneAlignFile != nullptr)
	{
	gzclose(m_gzNoneAlignFile);
	m_gzNoneAlignFile = nullptr;
	}

if(m_gzMultiAlignFile != nullptr)
	{
	gzclose(m_gzMultiAlignFile);
	m_gzMultiAlignFile = nullptr;
	}

if(m_gzSNPfile != nullptr)
	{
	gzclose(m_gzSNPfile);
	m_gzSNPfile = nullptr;
	}

if(m_gzSNPCentsfile != nullptr)
	{
	gzclose(m_gzSNPCentsfile);
	m_gzSNPCentsfile = nullptr;
	}

if(m_pszLineBuff != nullptr)
	{
	delete []m_pszLineBuff;
	m_pszLineBuff = nullptr;
	}
if(m_pAllocsIdentNodes != nullptr)
	{
	delete []m_pAllocsIdentNodes;
	m_pAllocsIdentNodes = nullptr;
	}

if(m_pAllocsMultiHitBuff != nullptr)
	{
	delete []m_pAllocsMultiHitBuff;
	m_pAllocsMultiHitBuff = nullptr;
	}

if(m_pAllocsMultiHitLoci != nullptr)
	{
	delete []m_pAllocsMultiHitLoci;
	m_pAllocsMultiHitLoci = nullptr;
	}
if(m_pSfxArray != nullptr)
	{
	delete m_pSfxArray;
	m_pSfxArray = nullptr;
	}
if(m_pPriorityRegionBED != nullptr)
	{
	delete m_pPriorityRegionBED;
	m_pPriorityRegionBED = nullptr;
	}

if(m_pReadHits != nullptr)
	{
#ifdef _WIN32
	free(m_pReadHits);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pReadHits != MAP_FAILED)
		munmap(m_pReadHits,m_AllocdReadHitsMem);
#endif
	m_pReadHits = nullptr;
	}

if(m_pMultiAll != nullptr)
	{
#ifdef _WIN32
	free(m_pMultiAll);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pMultiAll != MAP_FAILED)
		munmap(m_pMultiAll,m_AllocMultiAllMem);
#endif
	m_pMultiAll = nullptr;
	}

if(m_ppReadHitsIdx != nullptr)
	{
	delete []m_ppReadHitsIdx;
	m_ppReadHitsIdx = nullptr;
	}
if(m_pMultiHits != nullptr)
	{
#ifdef _WIN32
	free(m_pMultiHits);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pMultiHits != MAP_FAILED)
		munmap(m_pMultiHits,m_AllocdMultiHitsMem);
#endif
	m_pMultiHits = nullptr;
	}

if(m_pLociPValues != nullptr)
	{
#ifdef _WIN32
	free(m_pLociPValues);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pLociPValues != MAP_FAILED)
		munmap(m_pLociPValues,m_AllocPackedBaseAllelesMem);
#endif
	m_pLociPValues = nullptr;
	}

if(m_pPackedBaseAlleles != nullptr)
	{
#ifdef _WIN32
	free(m_pPackedBaseAlleles);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pPackedBaseAlleles != MAP_FAILED)
		munmap(m_pPackedBaseAlleles, m_AllocPackedBaseAllelesMem);
#endif
	m_pPackedBaseAlleles = nullptr;
	}

if(m_pChromSNPs != nullptr)
	{
	delete []m_pChromSNPs;
	m_pChromSNPs = nullptr;
	}

if(m_pLenDist != nullptr)
	{
	delete []m_pLenDist;
	m_pLenDist = nullptr;
	}

if(m_pContaminants != nullptr)
	{
	delete m_pContaminants;
	m_pContaminants = nullptr;
	}

if(m_pszCovSegBuff != nullptr)
	{
	delete []m_pszCovSegBuff;
	m_pszCovSegBuff = nullptr;
	}

DeleteMutexes();

Init();
m_TermBackgoundThreads = 0x0;	// can startup any background thread processing
}


// Load alignment loci base constraints from CSV file
// Expected file format is -
// <chrom>,<startloci>,<endloci>,<baselist>
// Whereby -
// chrom  must the name of a targeted chrom/sequence
// startloci is the start loci on the chrom
// endloci is the end loci on the chrom
// baselist names one or more bases which a read must align with for that aligned read to be accepted
//          if the baselist is '.' then no substitutions allowed at that loci
// examples -
// Chr1A,100,100,RA       reads aligning over Chr1A.100 will only be accepted if the read base is A or the target base
// Chr1A,100,100,C        reads aligning over Chr1A.100 will only be accepted if the read base is C
// Chr1A,100,100,CT       reads aligning over Chr1A.100 will only be accepted if the read base is C or T
// Chr1A,100,107,R        reads aligning over Chr1A.100 to Chr1A.107 will only be accepted if the read bases matches the target bases
// 

teBSFrsltCodes
CKAligner::LoadLociConstraints(char *pszLociConstraints)	// load loci constraints from file
{
int Rslt;
int NumLines;
int NumFields;
int CSVLineNum;

int Idx;

CCSVFile *pInFile = new CCSVFile;

char *pszTargChrom;
char szPrevTargChrom[cMaxGeneNameLen];

int ChromID;
int StartLoci;
int EndLoci;
uint8_t Constraint;

char *pszBases;
char *pBase;
char Base;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading loci base constraints from CSV file '%s' ...",pszLociConstraints);

if((Rslt=pInFile->Open(pszLociConstraints)) !=eBSFSuccess)
	{
	while(pInFile->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName, pInFile->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open '%s' for processing",pszLociConstraints);
	delete pInFile;
	return(eBSFerrOpnFile);
	}

if(m_pConstraintLoci == nullptr)
	m_pConstraintLoci = new tsConstraintLoci [cMaxConstrainedLoci];


NumLines = 0;
NumFields = 0;
szPrevTargChrom[0] = 0;
ChromID = 0;
while((Rslt= pInFile->NextLine()) > 0)	// onto next line containing fields
	{
	NumLines += 1;
	NumFields = pInFile->GetCurFields();
	if(NumFields < 4)
		{
		CSVLineNum = pInFile->GetLineNumber();
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected at least 4 fields at line %d in '%s', GetCurFields() returned '%d'",CSVLineNum,pszLociConstraints,NumFields);
		pInFile->Close();
		delete pInFile;
		return(eBSFerrFieldCnt);
		}

	if(NumLines == 1)		// slough 1st line if it is a title line, assuming title line contains only text values ... 
		{
		if(pInFile->IsLikelyHeaderLine())
			continue;
		}

	pInFile->GetText(1,&pszTargChrom);
	pInFile->GetInt(2,&StartLoci);
	pInFile->GetInt(3,&EndLoci);
	pInFile->GetText(4,&pszBases);

	// get chrom identifier
	if(ChromID == 0 || szPrevTargChrom[0] == 0 || stricmp(pszTargChrom,szPrevTargChrom))
		{
		strncpy(szPrevTargChrom,pszTargChrom,sizeof(szPrevTargChrom)-1);
		szPrevTargChrom[sizeof(szPrevTargChrom)-1] = '\0';
		if((ChromID = m_pSfxArray->GetIdent(szPrevTargChrom)) <= 0)
			{
			CSVLineNum = pInFile->GetLineNumber();
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to find matching indexed identifier for '%s' at line %d in '%s'",szPrevTargChrom,CSVLineNum,pszLociConstraints);
			pInFile->Close();
			delete pInFile;
			return(eBSFerrFieldCnt);
			}
		
		}

	if(StartLoci < 0 || StartLoci > EndLoci)
		{
		CSVLineNum = pInFile->GetLineNumber();
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Start loci must be >= 0 and <= end loci for '%s' at line %d in '%s'",szPrevTargChrom,CSVLineNum,pszLociConstraints);
		pInFile->Close();
		delete pInFile;
		return(eBSFerrFieldCnt);
		}

	if((uint32_t)EndLoci >= m_pSfxArray->GetSeqLen(ChromID))
		{
		CSVLineNum = pInFile->GetLineNumber();
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"End loci must be > targeted sequence length for '%s' at line %d in '%s'",szPrevTargChrom,CSVLineNum,pszLociConstraints);
		pInFile->Close();
		delete pInFile;
		return(eBSFerrFieldCnt);
		}

	Constraint = 0;
	pBase = pszBases;
	while((Base = *pBase++) != '\0')
		{
		switch(Base) {
			case 'a': case 'A':
				Constraint |= 0x01;
				break;
			case 'c': case 'C':
				Constraint |= 0x02;
				break;
			case 'g': case 'G':
				Constraint |= 0x04;
				break;
			case 't': case 'T':
				Constraint |= 0x08;
				break;
			case 'r': case 'R':
				Constraint |= 0x10;
				break;
			case ' ': case '\t':
				continue;
			default:
				CSVLineNum = pInFile->GetLineNumber();
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Illegal base specifiers for '%s' at line %d in '%s'",szPrevTargChrom,CSVLineNum,pszLociConstraints);
				pInFile->Close();
				delete pInFile;
				return(eBSFerrFieldCnt);
			}
		}
	if(Constraint == 0)
		{
		CSVLineNum = pInFile->GetLineNumber();
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Illegal base specifiers for '%s' at line %d in '%s'",szPrevTargChrom,CSVLineNum,pszLociConstraints);
		pInFile->Close();
		delete pInFile;
		return(eBSFerrFieldCnt);
		}

	// chrom, start, end and bases parsed
	for(Idx = 0; Idx < m_NumConstrainedChroms; Idx++)
		{
		if(m_ConstrainedChromIDs[Idx] == ChromID)
			break;
		}
	if(Idx == m_NumConstrainedChroms)
		{
		if(m_NumConstrainedChroms == cMaxConstrainedChroms)
			{
			CSVLineNum = pInFile->GetLineNumber();
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Number of constrained chroms would be more than max (%d) allowed for '%s' at line %d in '%s'",cMaxConstrainedChroms,szPrevTargChrom,CSVLineNum,pszLociConstraints);
			pInFile->Close();
			delete pInFile;
			return(eBSFerrFieldCnt);
			}
		m_ConstrainedChromIDs[Idx] = ChromID;
		m_NumConstrainedChroms += 1;
		}

	if(m_NumConstraintLoci == cMaxConstrainedLoci)
		{
		CSVLineNum = pInFile->GetLineNumber();
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Number of constrained loci would be more than max (%d) allowed for '%s' at line %d in '%s'",cMaxConstrainedLoci,szPrevTargChrom,CSVLineNum,pszLociConstraints);
		pInFile->Close();
		delete pInFile;
		return(eBSFerrFieldCnt);
		}
	m_pConstraintLoci[m_NumConstraintLoci].ChromID = ChromID;
	m_pConstraintLoci[m_NumConstraintLoci].StartLoci = StartLoci;
	m_pConstraintLoci[m_NumConstraintLoci].EndLoci = EndLoci;
	m_pConstraintLoci[m_NumConstraintLoci++].Constraint = Constraint;
	}

pInFile->Close();
delete pInFile;

// sort the constraints by chromid.start.end ascending
if(m_NumConstraintLoci > 1)
	m_mtqsort.qsort(m_pConstraintLoci,m_NumConstraintLoci,sizeof(tsConstraintLoci),SortConstraintLoci);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Completed loading %d loci base constraints for %d target sequences from CSV file '%s'",m_NumConstraintLoci,m_NumConstrainedChroms, pszLociConstraints);
return((teBSFrsltCodes)m_NumConstraintLoci);
}

teBSFrsltCodes
CKAligner::Disk2Hdr(char *pszRdsFile)			// read from disk and validate header
{
if(_lseeki64(m_hInFile,0,SEEK_SET)!=0)			// read in header..
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Seek failed to offset 0 - %s",pszRdsFile,strerror(errno));
	Reset(false);			// closes opened files..
	return(eBSFerrFileAccess);
	}

if(sizeof(tsBSFRdsHdr) != read(m_hInFile,&m_FileHdr,sizeof(tsBSFRdsHdr)))
	return(eBSFerrNotBioseq);

// header read, validate it as being a reads file header
if(tolower(m_FileHdr.Magic[0]) != 'b' ||
	tolower(m_FileHdr.Magic[1]) != 'i' ||
	tolower(m_FileHdr.Magic[2]) != 'o' ||
	tolower(m_FileHdr.Magic[3]) != 'r')
	return(eBSFerrNotBioseq);

	// can we handle this version?
if(m_FileHdr.Version < cBSFRdsVersionBack || m_FileHdr.Version > cBSFRdsVersion)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"%s opened as a preprocessed reads file - expected between version %d and %d, file version is %d",
					pszRdsFile,cBSFRdsVersionBack,cBSFRdsVersion,m_FileHdr.Version);
	Reset(false);			// closes opened files..
	return(eBSFerrFileVer);
	}

return(eBSFSuccess);
}

// NOTE: will only return SNP bases, e.g. aligned read sequence must not contain microInDel or span splice junctions
inline eSeqBase
CKAligner::AdjAlignSNPBase(tsReadHit *pReadHit,	// aligned read
		   uint32_t ChromID,			// read expected to have aligned to this chromosome
			uint32_t Loci)            // base to be returned is at this alignment loci, base will be complemented if antisense alignment
{
uint32_t AdjStartLoc;
uint32_t AdjEndLoc;
etSeqBase Base;
uint8_t *pBases;
tsSegLoci *pSeg;
tsHitLoci *pHit;

pHit = &pReadHit->HitLoci.Hit;
if(pHit->FlgInDel || pHit->FlgSplice)
	return(eBaseEOS);

pSeg = &pReadHit->HitLoci.Hit.Seg[0];
if(pSeg->ChromID != ChromID)
	return(eBaseEOS);

if(pSeg->Strand == '+')
	{
	AdjStartLoc=((uint32_t)pSeg->MatchLoci + pSeg->TrimLeft);
	AdjEndLoc=((uint32_t)pSeg->MatchLoci + (pSeg->MatchLen - pSeg->TrimRight - 1));
	}
else
	{
	AdjStartLoc=((uint32_t)pSeg->MatchLoci + pSeg->TrimRight);
	AdjEndLoc=((uint32_t)pSeg->MatchLoci + (pSeg->MatchLen - pSeg->TrimLeft - 1));
	}

if(AdjStartLoc > Loci || AdjEndLoc < Loci)
	return(eBaseEOS);

pBases = &pReadHit->Read[pReadHit->DescrLen+1];

if(pSeg->Strand == '+')
	Base = pBases[Loci - pSeg->MatchLoci]  & 0x07;
else
	{
	Base = pBases[pSeg->MatchLoci + pSeg->MatchLen - Loci - 1] & 0x07;
	switch(Base) {
		case eBaseA: Base = eBaseT; break;
		case eBaseC: Base = eBaseG; break;
		case eBaseG: Base = eBaseC; break;
		case eBaseT: Base = eBaseA; break;
		default:
			break;
		}
	}
return((eSeqBase)Base);
}

inline uint32_t
CKAligner::AdjStartLoci(tsSegLoci *pSeg)
{
if(pSeg->Strand == '+')
	return((uint32_t)pSeg->MatchLoci + pSeg->TrimLeft);
else
	return((uint32_t)pSeg->MatchLoci + pSeg->TrimRight);
}

inline uint32_t
CKAligner::AdjEndLoci(tsSegLoci *pSeg)
{
if(pSeg->Strand == '+')
	return((uint32_t)pSeg->MatchLoci + (pSeg->MatchLen - pSeg->TrimRight - 1));
else
	return((uint32_t)pSeg->MatchLoci  + (pSeg->MatchLen - pSeg->TrimLeft - 1));
}

inline uint32_t
CKAligner::AdjHitLen(tsSegLoci *pSeg)
{
return((uint32_t)pSeg->MatchLen - pSeg->TrimLeft - pSeg->TrimRight);
}

inline uint32_t
CKAligner::AdjAlignStartLoci(tsHitLoci *pHit)
{
tsSegLoci *pSeg = &pHit->Seg[0];
if(pSeg->Strand == '+')
	return((uint32_t)pSeg->MatchLoci + pSeg->TrimLeft);
else
	return((uint32_t)pSeg->MatchLoci + pSeg->TrimRight);
}

inline uint32_t
CKAligner::AdjAlignEndLoci(tsHitLoci *pHit)
{
tsSegLoci *pSeg;
if(pHit->FlgInDel || pHit->FlgSplice)
	pSeg = &pHit->Seg[1];
else
	pSeg = &pHit->Seg[0];
if(pSeg->Strand == '+')
	return((uint32_t)pSeg->MatchLoci + (pSeg->MatchLen - pSeg->TrimRight - 1));
else
	return((uint32_t)pSeg->MatchLoci  + (pSeg->MatchLen - pSeg->TrimLeft - 1));
}

inline uint32_t
CKAligner::AdjAlignHitLen(tsHitLoci *pHit)
{
uint32_t HitLen;
tsSegLoci *pSeg = &pHit->Seg[0];
HitLen = (uint32_t)pSeg->MatchLen - pSeg->TrimLeft - pSeg->TrimRight;
if(pHit->FlgInDel || pHit->FlgSplice)
	{
	pSeg = &pHit->Seg[1];
	HitLen += (uint32_t)pSeg->MatchLen - pSeg->TrimLeft - pSeg->TrimRight;
	}
return(HitLen);
}

inline int
CKAligner::AdjAlignMismatches(tsHitLoci *pHit)
{
int Mismatches;
tsSegLoci *pSeg = &pHit->Seg[0];
Mismatches = pSeg->TrimMismatches;
if(pHit->FlgInDel || pHit->FlgSplice)
	{
	pSeg = &pHit->Seg[1];
	Mismatches += pSeg->TrimMismatches;
	}
return(Mismatches);
}

// AutoTrimFlanks
// Intent is that this will be useful for delimiting RNAseq reads covering exon boundaries
// Autotrimmed aligned reads must be at least 50% of their untrimmed length or they will be discarded; exception is that if paired end processing then
// trimming of these reads is limited so that at least 1/3rd of the read is retained as a central core
int
CKAligner::AutoTrimFlanks(int MinFlankExacts) // Autotrim back aligned read flanks until there are at least MinFlankExacts exactly matching bases in the flanks
{
int Rslt;
tsReadHit *pReadHit;
uint32_t SeqIdx;
uint32_t Idx;
int MinTrimmedLen;
uint8_t *pSeq;
uint8_t ReadSeq[cMaxFastQSeqLen + 10000];
uint8_t TargSeq[cMaxFastQSeqLen + 10000];

if((Rslt=SortReadHits(eRSMHitMatch,false)) < eBSFSuccess)
	{
	Reset(false);
	return(Rslt);
	}

m_ElimPlusTrimed = 0;
m_ElimMinusTrimed = 0;
if(MinFlankExacts > 0)	// do  3' and 5' autotrim? Note that currently can't trim multiseg hits
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Starting 5' and 3' flank sequence autotrim processing...");
	etSeqBase *pHitSeq;
	pReadHit = nullptr;
	int ExactLen;
	int LeftOfs;
	int RightOfs;
	uint32_t MatchLen;
	int TrimMismatches;


	while((pReadHit = IterReads(pReadHit))!=nullptr)
		{
		pReadHit->HitLoci.FlagTR = 0;
		if(pReadHit->NAR != eNARAccepted || pReadHit->HitLoci.FlagSegs == 1 || pReadHit->HitLoci.Hit.FlgChimeric == 1)
			continue;
		MatchLen = pReadHit->HitLoci.Hit.Seg[0].MatchLen;
		if(MatchLen != pReadHit->ReadLen) // will only be different if read already determined to be spliced or containing microInDel
			{
			pReadHit->NumHits = 0;
			if(pReadHit->HitLoci.Hit.Seg[0].Strand == '+')		// treating as if trimmed
				m_ElimPlusTrimed += 1;
			else
				m_ElimMinusTrimed += 1;
			continue;
			}

		MinTrimmedLen = (MatchLen+1)/2;							// post trimming the alignment must be at least this length
		if(MinTrimmedLen < 15)
			MinTrimmedLen = 15;

				// copy read into ReadSeq masking any hi-order phred scores (in bits 4..7) and repeat mask (bit 3) out
		pHitSeq = &pReadHit->Read[pReadHit->DescrLen+1];
		pSeq = ReadSeq;
		for(SeqIdx = 0; SeqIdx < MatchLen; SeqIdx++)
			*pSeq++ = *pHitSeq++ & 0x07;

		if(m_bIsSOLiD)
			{
			uint32_t Loci = AdjStartLoci(&pReadHit->HitLoci.Hit.Seg[0]);
			if(pReadHit->HitLoci.Hit.Seg[0].Strand == '+')
				Loci += 1;
			if(m_pSfxArray->GetColorspaceSeq(pReadHit->HitLoci.Hit.Seg[0].ChromID,
									Loci,
									TargSeq,MatchLen) == 0) // get colorspace sequence
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: GetColorspaceSeq() for chrom %d, Loci %d, matchlen %d failed",
										pReadHit->HitLoci.Hit.Seg[0].ChromID,Loci,MatchLen);
				Reset(false);
				return(eBSFerrMem);
				}

			if(pReadHit->HitLoci.Hit.Seg[0].Strand == '-')
				CSeqTrans::ReverseSeq(MatchLen,TargSeq);

			// convert target assembly sequence read into colorspace
			pSeq = ReadSeq;
			uint8_t PrvBase = *pSeq;
			for(SeqIdx = 1; SeqIdx <= MatchLen; SeqIdx++,pSeq++)
				{
				*pSeq = SOLiDmap[PrvBase][pSeq[1]];
				PrvBase = pSeq[1];
				}
			MatchLen -= 1;
			}
		else
			{
			// get basespace sequence
			m_pSfxArray->GetSeq(pReadHit->HitLoci.Hit.Seg[0].ChromID,(uint32_t)(uint64_t)pReadHit->HitLoci.Hit.Seg[0].MatchLoci,TargSeq,(uint32_t)MatchLen);
			if(pReadHit->HitLoci.Hit.Seg[0].Strand == '-')
				CSeqTrans::ReverseComplement(MatchLen,TargSeq);
			}

		// trim from 5' towards 3'
		pSeq = ReadSeq;
		pHitSeq = TargSeq;
		ExactLen = 0;
		TrimMismatches = 0;
		int PEmincore;
		if(m_PEproc == ePEdefault)
			PEmincore = MatchLen;
		else
			PEmincore = MatchLen / 3;
		for(Idx = 0; Idx <= (MatchLen-(uint32_t)MinTrimmedLen) && Idx < (uint32_t)PEmincore; Idx++,pSeq++,pHitSeq++)
			{
			if(*pSeq != *pHitSeq)
				{
				ExactLen = 0;
				TrimMismatches += 1;
				continue;
				}
			ExactLen += 1;
			if(ExactLen == MinFlankExacts)
				break;
			}

		if(m_PEproc == ePEdefault)
			{
			// if can't trim to a 5' flank of at least MinFlankExacts then this read is to be sloughed
			if((Idx + MinTrimmedLen) > MatchLen || ExactLen < MinFlankExacts)
				{
				pReadHit->NumHits = 0;
				pReadHit->NAR = eNARTrim;
				if(pReadHit->HitLoci.Hit.Seg[0].Strand == '+')
					m_ElimPlusTrimed += 1;
				else
					m_ElimMinusTrimed += 1;
				continue;
				}
			}
		LeftOfs = Idx - (MinFlankExacts - 1);

		// trim from 3' back towards 5'
		pHitSeq = &TargSeq[MatchLen-1];
		pSeq = &ReadSeq[MatchLen-1];
		ExactLen = 0;
		if(m_PEproc == ePEdefault)
			PEmincore = 0;
		else
			PEmincore = (MatchLen * 2) / 3;

		for(Idx = MatchLen-1; Idx >= (uint32_t)(LeftOfs+MinTrimmedLen) && Idx > (uint32_t)PEmincore; Idx--,pSeq--,pHitSeq--)
			{
			if(*pSeq != *pHitSeq)
				{
				ExactLen = 0;
				TrimMismatches += 1;
				continue;
				}
			ExactLen += 1;
			if(ExactLen == MinFlankExacts)
				break;
			}

		if(m_PEproc == ePEdefault)
			{
			// if can't trim to a 3' flank of at least MinFlankExacts then this read is to be sloughed
			if(ExactLen != MinFlankExacts || Idx < (uint32_t)(LeftOfs + MinTrimmedLen))
				{
				pReadHit->NumHits = 0;
				pReadHit->NAR = eNARTrim;
				if(pReadHit->HitLoci.Hit.Seg[0].Strand == '+')
					m_ElimPlusTrimed += 1;
				else
					m_ElimMinusTrimed += 1;
				continue;
				}
			}
		RightOfs = Idx + MinFlankExacts;

		// left and right offset in the read at which the exact matching flanks start is now known
		if(m_bIsSOLiD)
			{
			MatchLen += 1;
			RightOfs += 1;
			}

		pReadHit->HitLoci.Hit.Seg[0].TrimLeft = LeftOfs;
		pReadHit->HitLoci.Hit.Seg[0].TrimRight = MatchLen - RightOfs;
		if(LeftOfs || (MatchLen - RightOfs)) // any 5' or 3' trimming was required?
			{
			pReadHit->HitLoci.Hit.Seg[0].TrimMismatches = pReadHit->HitLoci.Hit.Seg[0].Mismatches - TrimMismatches;
			pReadHit->HitLoci.FlagTR = 1;
			}
		else
			pReadHit->HitLoci.FlagTR = 0;
		}

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Finished 5' and 3' flank sequence autotriming, %d plus strand and %d minus strand aligned reads removed",m_ElimPlusTrimed,m_ElimMinusTrimed);

	if((Rslt=SortReadHits(eRSMHitMatch,false,true)) < eBSFSuccess)
		{
		Reset(false);
		return(Rslt);
		}

	if(gProcessingID)
		{
		gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"Filtering",ePTInt32,sizeof(m_ElimPlusTrimed),"5' trimmed",&m_ElimPlusTrimed);
		gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"Filtering",ePTInt32,sizeof(m_ElimMinusTrimed),"3' trimmed",&m_ElimMinusTrimed);
		}
	}
return(eBSFSuccess);
}

// ReportChimerics
// Report reads which have been chimerically aligned
const int cChimericSeqBuffLen = 0x07fffff;  // use this sized buffer when reporting the chimeric sequences

int
CKAligner::ReportChimerics(char *pszChimericSeqFile)			// report chimerically trimmed read sequences to file pszChimericSeqFile
{
int hChimerics;
char Strand;
char szChromName[128];
uint32_t Ofs;
uint32_t PrevChromID;

char *pszLineBuff;
int BuffIdx;
tsReadHit *pCurReadHit;
tsReadHit *pNxtReadHit;
uint32_t TrimLeft;
uint32_t TrimRight;
uint32_t NumLeftRightTrimmed;
uint32_t SeqIdx;
uint32_t NumTrimmed;
uint32_t NumLeftTrimmed;
uint32_t NumRightTrimmed;
uint32_t MatchLen;
uint8_t *pSeq;
uint8_t *pChimericSeq;
uint8_t *pTo;
int CopyLen;
uint32_t NumReads;

if(pszChimericSeqFile == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: No chimeric sequences file specified");
	return(eBSFerrCreateFile);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Starting to report chimeric flank trimmed read sequences to file '%s' ...", pszChimericSeqFile);
PrevChromID = 0;
hChimerics = -1;
pszLineBuff = nullptr;
if(pszChimericSeqFile != nullptr && pszChimericSeqFile[0] != '\0')
	{
#ifdef _WIN32
	hChimerics = open(pszChimericSeqFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE) );
#else
	if((hChimerics = open(pszChimericSeqFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE))!=-1)
		if(ftruncate(hChimerics,0)!=0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate chimeric sequences file %s - %s",pszChimericSeqFile,strerror(errno));
			return(eBSFerrCreateFile);
			}
#endif

	if(hChimerics < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate chimeric sequences file '%s'",pszChimericSeqFile);
		return(eBSFerrCreateFile);
		}	
	if((pszLineBuff = new char [cChimericSeqBuffLen]) == nullptr)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to allocate memory %d for chimeric sequences buffering",cChimericSeqBuffLen);
		return(eBSFerrMem);
		}	
	BuffIdx = sprintf(pszLineBuff,"\"Chrom\",\"Loci\",\"Strand\",\"ReadDescr\",\"5'TrimLen\",\"5'TrimSeq\",\"AlignLen\",\"AlignSeq\",\"3'TrimLen\",\"3'TrimSeq\"\n");
	if(!CUtility::RetryWrites(hChimerics,pszLineBuff,BuffIdx))
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"RetryWrites() error");
		return(eBSFerrWrite);
		}
	BuffIdx = 0;
	}
NumTrimmed = 0;
NumLeftTrimmed = 0;
NumRightTrimmed = 0;
NumLeftRightTrimmed = 0;

pCurReadHit = m_pReadHits;
pTo = (uint8_t *)pCurReadHit;
NumReads = 0;
while(pCurReadHit != nullptr) {
	NumReads += 1;
	CopyLen = sizeof(tsReadHit) + pCurReadHit->DescrLen + pCurReadHit->ReadLen;
	if(pCurReadHit->ReadID != m_FinalReadID)
		pNxtReadHit = (tsReadHit *)((uint8_t *)pCurReadHit + CopyLen);
	else
		pNxtReadHit = nullptr;
	pCurReadHit->HitLoci.FlagTR = 0;
	if(pCurReadHit->NAR != eNARAccepted || pCurReadHit->HitLoci.Hit.FlgChimeric == 0) // only reporting reads which were accepted as being chimeric aligned
		{
		pCurReadHit = pNxtReadHit;
		continue;
		}

	TrimLeft = pCurReadHit->HitLoci.Hit.Seg[0].TrimLeft;
	TrimRight = pCurReadHit->HitLoci.Hit.Seg[0].TrimRight;

	if(TrimLeft == 0 && TrimRight == 0)	// if accepted as a chimeric then should have had at least one flank to be trimmed, treat as full length match ...
		{
		pCurReadHit->HitLoci.Hit.FlgChimeric = 0;
		pCurReadHit = pNxtReadHit;
		continue;
		}
	if(TrimLeft > 0)
		NumLeftTrimmed += 1;
	if(TrimRight > 0)
		NumRightTrimmed += 1;
	if(TrimLeft > 0 && TrimRight > 0)
		NumLeftRightTrimmed += 1;

	MatchLen = pCurReadHit->ReadLen - (TrimLeft + TrimRight);
	pSeq = &pCurReadHit->Read[pCurReadHit->DescrLen+1];
	pChimericSeq = pSeq + TrimLeft;
	
	if(hChimerics != -1)
		{
		Strand = pCurReadHit->HitLoci.Hit.Seg[0].Strand;
		Ofs = AdjAlignStartLoci(&pCurReadHit->HitLoci.Hit);
		if(PrevChromID == 0 || pCurReadHit->HitLoci.Hit.Seg[0].ChromID != PrevChromID)
			{
			m_pSfxArray->GetIdentName(pCurReadHit->HitLoci.Hit.Seg[0].ChromID,sizeof(szChromName),szChromName);
			PrevChromID = pCurReadHit->HitLoci.Hit.Seg[0].ChromID;
			}
		BuffIdx += sprintf(&pszLineBuff[BuffIdx],"\"%s\",%d,\"%c\",",szChromName,Ofs,Strand);
		}
	if(TrimLeft > 0)
		{
		if(hChimerics != -1)
			{
			BuffIdx += sprintf(&pszLineBuff[BuffIdx],"\"%s\",%d,\"",(char *)pCurReadHit->Read,TrimLeft);
			CSeqTrans::MapSeq2Ascii(pSeq,TrimLeft,&pszLineBuff[BuffIdx]);
			BuffIdx += TrimLeft;
			BuffIdx += sprintf(&pszLineBuff[BuffIdx],"\",%d,\"",MatchLen);
			CSeqTrans::MapSeq2Ascii(pChimericSeq,MatchLen,&pszLineBuff[BuffIdx]);
			BuffIdx += MatchLen;
			BuffIdx += sprintf(&pszLineBuff[BuffIdx],"\",%d,\"",TrimRight);
			if(TrimRight)
				CSeqTrans::MapSeq2Ascii(&pChimericSeq[MatchLen],TrimRight,&pszLineBuff[BuffIdx]);
			BuffIdx += TrimRight;
			BuffIdx += sprintf(&pszLineBuff[BuffIdx],"\"\n");
			}
		for(SeqIdx = 0; SeqIdx < MatchLen; SeqIdx++)
			*pSeq++ = *pChimericSeq++;
		}
	else
		if(hChimerics != -1)
			{
			BuffIdx += sprintf(&pszLineBuff[BuffIdx],"\"%s\",0,\"\",%d,\"",(char *)pCurReadHit->Read,MatchLen);
			CSeqTrans::MapSeq2Ascii(pSeq,MatchLen,&pszLineBuff[BuffIdx]);
			BuffIdx += MatchLen;
			BuffIdx += sprintf(&pszLineBuff[BuffIdx],"\",%d,\"",TrimRight);
			if(TrimRight)
				CSeqTrans::MapSeq2Ascii(&pSeq[MatchLen],TrimRight,&pszLineBuff[BuffIdx]);
			BuffIdx += TrimRight;
			BuffIdx += sprintf(&pszLineBuff[BuffIdx],"\"\n");
			}

	if(hChimerics != -1 && BuffIdx > (cChimericSeqBuffLen / 2))
		{
		if(!CUtility::RetryWrites(hChimerics,pszLineBuff,BuffIdx))
			{
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"RetryWrites() error");
			return(eBSFerrWrite);
			}
		BuffIdx = 0;
		}

	pCurReadHit = pNxtReadHit;
	NumTrimmed += 1;
	}

if(hChimerics != -1)
	{
	if(BuffIdx)
		if(!CUtility::RetryWrites(hChimerics,pszLineBuff,BuffIdx))
			{
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"RetryWrites() error");
			return(eBSFerrWrite);
			}

#ifdef _WIN32
	_commit(hChimerics);
#else
	fsync(hChimerics);
#endif
	close(hChimerics);
	hChimerics = -1;
	}
if(pszLineBuff != nullptr)
	delete []pszLineBuff;
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Completed reporting to file '%s'  %u chimeric flanks, %u 5', %u 3', %u both 5' and 3' flanks trimmed",pszChimericSeqFile,NumTrimmed,NumLeftTrimmed,NumRightTrimmed, NumLeftRightTrimmed);
return(eBSFSuccess);
}



// PCR5PrimerCorrect
// Intent is that this will be useful for aligning with tight substitution constraints whereby end PCR hexamer artefacts dominate the base substitution profile
// Substitutions in the first 5' KLen bases are progressively corrected, using the targeted sequence as template
int 
CKAligner::PCR5PrimerCorrect(int MaxAllowedSubRate,	// after corrections overall sub rate for read must be no more than this 
					int KLen) // progressively correct substitutions in first Klen 5' read bases - assumed to be PCR random primer artefacts - until overall read meets MaxAllowedSubRate
{
int Rslt;
tsReadHit *pReadHit;
uint8_t *pSeq;
int MaxMMs;
int Ofs;
int CurMMs;
etSeqBase Base;
int NumCorrectedReads;
int NumCorrectedBases;
int NumUnacceptedReads;
etSeqBase *pHitSeq;
pReadHit = nullptr;
uint32_t MatchLen;
uint8_t TargSeq[cMaxFastQSeqLen+ 10000];


if(m_bIsSOLiD || KLen < 1)
	return(0);

if((Rslt=SortReadHits(eRSMHitMatch,false)) < eBSFSuccess)
	{
	Reset(false);
	return(Rslt);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Starting PCR 5' %dbp primer correction processing targeting substitution rate of %d ...",KLen, MaxAllowedSubRate);

NumCorrectedReads = 0;
NumCorrectedBases = 0;
NumUnacceptedReads = 0;
while((pReadHit = IterReads(pReadHit))!=nullptr)
	{
	if(pReadHit->NAR != eNARAccepted || pReadHit->HitLoci.FlagSegs == 1)
		continue;

	MaxMMs = ((MaxAllowedSubRate * pReadHit->ReadLen)+50) / 100;

	if(pReadHit->LowMMCnt <= MaxMMs)				// if already meeting targeted sub rate then no correction required
		continue;

	MatchLen = pReadHit->HitLoci.Hit.Seg[0].MatchLen;
	if(MatchLen != pReadHit->ReadLen)     // should usually be the same but some earlier processing function may have modified the length so don't attempt correction
		continue;

	// get target sequence
	m_pSfxArray->GetSeq(pReadHit->HitLoci.Hit.Seg[0].ChromID,(uint32_t)(uint64_t)pReadHit->HitLoci.Hit.Seg[0].MatchLoci,TargSeq,(uint32_t)MatchLen);
	if(pReadHit->HitLoci.Hit.Seg[0].Strand == '-')
		CSeqTrans::ReverseComplement(MatchLen,TargSeq);

	// progressively look for mismatches over first KLen of read and correct until meeting  MaxAllowedSubRate
	pSeq =  &pReadHit->Read[pReadHit->DescrLen+1];
	pHitSeq = TargSeq;
	CurMMs = pReadHit->LowMMCnt;
	for(Ofs = 0; Ofs < (int)KLen; Ofs++,pSeq++,pHitSeq++)
		{
		if((*pSeq & 0x07) != *pHitSeq)					
			{
			if(--CurMMs <= MaxMMs)
				break;
			}
		}
	if(CurMMs <= MaxMMs)
		{
		pSeq =  &pReadHit->Read[pReadHit->DescrLen+1];
		pHitSeq = TargSeq;
		CurMMs = pReadHit->LowMMCnt;
		for(Ofs = 0; Ofs < (int)KLen; Ofs++,pSeq++,pHitSeq++)
			{
			if(((Base = *pSeq) & 0x07) != *pHitSeq)					// required an aligner induced sub?
				{
				*pSeq = (Base & 0xf8) | *pHitSeq;
				NumCorrectedBases += 1;
				if(--CurMMs <= MaxMMs)
					break;
				}
			}
		pReadHit->LowMMCnt = CurMMs;
		pReadHit->HitLoci.Hit.Seg[0].Mismatches = CurMMs;
		NumCorrectedReads += 1;
		}
	else
		{
		pReadHit->NumHits = 0;
		pReadHit->NAR = eNARNoHit;
		NumUnacceptedReads += 1;
		}
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Completed PCR 5' primer correction, %d reads with %d bases corrected, %d reads with excessive substitutions rejected",NumCorrectedReads,NumCorrectedBases, NumUnacceptedReads);

if((Rslt=SortReadHits(eRSMHitMatch,false,true)) < eBSFSuccess)
	{
	Reset(false);
	return(Rslt);
	}

if(gProcessingID)
	{
	gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"PCR5PrimerCorrect",ePTInt32,sizeof(NumCorrectedReads),"NumCorrectedReads",&NumCorrectedReads);
	gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"PCR5PrimerCorrect",ePTInt32,sizeof(NumCorrectedBases),"NumCorrectedBases",&NumCorrectedBases);
	gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"PCR5PrimerCorrect",ePTInt32,sizeof(NumUnacceptedReads),"NumUnacceptedReads",&NumUnacceptedReads);
	}

return(eBSFSuccess);
}


//
int
CKAligner::NumUniqueAlignedReads(void)		// return count of accepted uniquely aligned reads
{
tsReadHit *pReadHit = nullptr;
int NumUniques = 0;
while((pReadHit = IterReads(pReadHit))!=nullptr)
	{
	if(pReadHit->NAR == eNARAccepted)
		NumUniques += 1;
	}
return(NumUniques);
}

int
CKAligner::DedupeNonuniqueReads(void)	// dedupe reads
{
int NumSubDups;

tsReadHit *pReadHit1;
tsReadHit *pReadHitMark;
tsReadHit *pReadHit = nullptr;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Deduping aligned reads which after substitutions are no longer unique");
SortReadHits(eRSMHitMatch,false);
NumSubDups = 0;
while((pReadHit = IterSortedReads(pReadHit))!=nullptr)
	{
	if(pReadHit->NAR != eNARAccepted)		// only interested in reads accepted as aligned
		continue;

	pReadHit1 = pReadHit;			// start iterating forward checking for multiple probe hits onto same locus
	pReadHitMark = pReadHit;
	while((pReadHit1 = IterSortedReads(pReadHit1)) != nullptr)
		{
		if(pReadHit1->NAR != eNARAccepted)	// only interested in reads with a single hits
			break;

		if(	pReadHit->HitLoci.Hit.Seg[0].ChromID == pReadHit1->HitLoci.Hit.Seg[0].ChromID && // iterated probe hitting same locus as original probe?
			AdjStartLoci(&pReadHit->HitLoci.Hit.Seg[0]) == AdjStartLoci(&pReadHit1->HitLoci.Hit.Seg[0]) &&
			pReadHit->HitLoci.Hit.Seg[0].Strand == pReadHit1->HitLoci.Hit.Seg[0].Strand)
			{
			// the iterated probe will be of equal or lower level than original probe
			// simply treat iterated probe as though it had multiple matches
			pReadHit1->NumHits = 0;
			pReadHit1->LowHitInstances = 0;
			pReadHit1->NAR = eNARNonUnique;
			pReadHitMark = pReadHit1;
			NumSubDups += 1;
			}
		else
			break;
		}
	pReadHit = pReadHitMark;
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Deduping completed, removed %d non-unique reads",NumSubDups);

if(NumSubDups)
	SortReadHits(eRSMHitMatch,false,true);

if(gProcessingID)
	gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"Filtering",ePTInt32,sizeof(NumSubDups),"Nonunique",&NumSubDups);

return(eBSFSuccess);
}


//
// PCR results in differential read duplications with some sequence templates being many fold duplicated relative to others
// Differential artifact reduction is very crude!
// Each unique loci to which a read is aligned is iterated, then the number of unique read alignment sites up/down stream within WinLen bases is counted and this count is then used
// to limit the number of reads attributed to the loci currently being iterated.
// The limit is set to be a function of the maximum of the up/downstream unique aligned loci within WinLen...
int
CKAligner::ReducePCRduplicates(int WinLen)		// remove potential PCR artefacts
{
int Rslt;
int NumSubDups;
int CurChromID;
int CurLen;
int CurStart;
uint8_t CurStrand;
int UpUniques;
int DnUniques;
int LimitDups;
int PropWin;

tsReadHit *pReadHit1;
tsReadHit *pReadHitMark;
tsReadHit *pReadHit = nullptr;

// sort reads by match loci
if((Rslt=SortReadHits(eRSMHitMatch,false)) < eBSFSuccess)
	{
	Reset(false);
	return(Rslt);
	}

NumSubDups = 0;
while((pReadHit = IterSortedReads(pReadHit))!=nullptr)
	{
	if(pReadHit->NAR != eNARAccepted)		// only interested in reads with unique hits
		continue;

	if(WinLen > 0)
		{
		UpUniques = NumUpUniques(pReadHit,WinLen,true);
		DnUniques = NumDnUniques(pReadHit,WinLen,true);
		LimitDups = max(UpUniques,DnUniques);
		PropWin = (int)(((double)LimitDups/WinLen) * 100.0);
		if(PropWin < 5)
			LimitDups = 1;
		else
			if(PropWin <= 10)
				LimitDups = 2;
			else
				if(PropWin <= 20)
					LimitDups = 3;
				else
					if(PropWin <= 40)
						LimitDups = 4;
					else
						if(PropWin <= 60)
							LimitDups = 5;
						else
							if(PropWin <= 80)
								LimitDups = 10;
							else
								LimitDups = 50;
		}
	else
		LimitDups = 0;

	pReadHit1 = pReadHit;			// start iterating forward checking for multiple probe hits onto same locus
	pReadHitMark = pReadHit;
	CurChromID = pReadHit->HitLoci.Hit.Seg[0].ChromID;
	CurStart = AdjStartLoci(&pReadHit->HitLoci.Hit.Seg[0]);
	CurLen = AdjHitLen(&pReadHit->HitLoci.Hit.Seg[0]);
	CurStrand = pReadHit->HitLoci.Hit.Seg[0].Strand;
	while((pReadHit1 = IterSortedReads(pReadHit1)) != nullptr)
		{
		if(pReadHit1->NAR != eNARAccepted)	// only interested in reads with a unique hits
			continue;

		if(CurChromID == pReadHit1->HitLoci.Hit.Seg[0].ChromID && // iterated probe hitting same locus as original probe?
			CurStart == AdjStartLoci(&pReadHit1->HitLoci.Hit.Seg[0]) &&
			CurStrand == pReadHit1->HitLoci.Hit.Seg[0].Strand)
			{
			if(CurLen != AdjHitLen(&pReadHit1->HitLoci.Hit.Seg[0]))
				continue;
			if(LimitDups > 0)
				{
				LimitDups -= 1;
				continue;
				}

			pReadHit1->NumHits = 0;
			pReadHit1->LowHitInstances = 0;
			pReadHit1->NAR = eNARPCRdup;
			pReadHitMark = pReadHit1;
			NumSubDups += 1;
			}
		else
			break;
		}
	pReadHit = pReadHitMark;
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Removed %d potential PCR artefact reads",NumSubDups);
if(NumSubDups)
	SortReadHits(eRSMHitMatch,false,true);
if(gProcessingID)
	gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"Filtering",ePTInt32,sizeof(NumSubDups),"PCRartefacts",&NumSubDups);
return(eBSFSuccess);
}


int
CKAligner::RemoveOrphanSpliceJuncts(int SpliceJunctLen)	// remove unsupported orphan splice junctions
{
int Idx;
if(SpliceJunctLen > 0)
	{
	SortReadHits(eRSMHitMatch,false);
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Filtering out any orphan (unsupported by >= 2 reads) splice junction reads");
	int NumSpliceJuncs = 0;
	int	NumSpliceAccepted = 0;
	int	NumSpliceNotAccepted = 0;
	tsSegJuncts *pSegJuncts = nullptr;
	tsSegJuncts *pJunct;
	tsSegJuncts *pNxtJunct;
	tsReadHit *pCurHit = nullptr;
	tsReadHit *pReadHit = nullptr;
	while((pReadHit = IterReads(pReadHit))!=nullptr)
		{
		if(pReadHit->NAR != eNARAccepted || pReadHit->HitLoci.Hit.FlgSplice == 0)
			continue;
		pCurHit = pReadHit;
		NumSpliceJuncs += 1;
		}
	if(NumSpliceJuncs > 1)
		{
		// allocate for this number of splice junctions, then sort and dedupe
		pSegJuncts = new tsSegJuncts[NumSpliceJuncs];
		if(pSegJuncts == nullptr)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to allocate memory (%d bytes) for splice junction processing",NumSpliceJuncs * sizeof(tsSegJuncts));
			Reset(false);
			return(eBSFerrMem);
			}
		pJunct = pSegJuncts;
		pReadHit = nullptr;
		while((pReadHit = IterReads(pReadHit))!=nullptr)
			{
			if(pReadHit->NAR != eNARAccepted || pReadHit->HitLoci.Hit.FlgSplice == 0)
				continue;
			pJunct->pRead = pReadHit;
			pJunct->ChromID = pReadHit->HitLoci.Hit.Seg[0].ChromID;
			pJunct->Starts = AdjEndLoci(&pReadHit->HitLoci.Hit.Seg[0]);
			pJunct->Ends = AdjStartLoci(&pReadHit->HitLoci.Hit.Seg[1]);
			pJunct->Cnt = 1;
			pJunct += 1;
			}
		// now sort ascending by chrom, start, end
		m_mtqsort.qsort(pSegJuncts,NumSpliceJuncs,sizeof(tsSegJuncts),SortSegJuncts);
		// iterate and determine multiplicity of junctions
		pJunct = pSegJuncts;
		pNxtJunct = &pSegJuncts[1];
		for(Idx = 0; Idx < (NumSpliceJuncs-1); Idx++,pJunct++,pNxtJunct++)
			{
			if(pJunct->ChromID == pNxtJunct->ChromID &&
			   (pJunct->Starts <= (pNxtJunct->Starts + 3) && pJunct->Starts >= (pNxtJunct->Starts - 3)) &&
			   (pJunct->Ends <= (pNxtJunct->Ends + 3) && pJunct->Ends >= (pNxtJunct->Ends - 3)))
				{
				pJunct->pRead->HitLoci.Hit.FlgNonOrphan = 1;
				pNxtJunct->pRead->HitLoci.Hit.FlgNonOrphan = 1;
				}
			}
		// junctions not supported by at least two reads are to be treated as simply unaligned
		pJunct = pSegJuncts;
		for(Idx = 0; Idx < NumSpliceJuncs; Idx++,pJunct++)
			{
			if(pJunct->pRead->HitLoci.Hit.FlgNonOrphan != 1)
				{
				pJunct->pRead->NAR = eNARSpliceJctn;
				pJunct->pRead->NumHits = 0;	// treat as unaligned
				pJunct->pRead->LowHitInstances = 0;
				NumSpliceNotAccepted += 1;
				}
			else
				NumSpliceAccepted += 1;
			}
		delete []pSegJuncts;
		}
	else
		if(NumSpliceJuncs > 0 && pCurHit != nullptr)
			{
			pCurHit->NAR = eNARSpliceJctn;
			pCurHit->NumHits = 0;	// treat as unaligned
			pCurHit->LowHitInstances = 0;
			NumSpliceNotAccepted += 1;
			}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"From %d reads with putative splice sites %d orphans were removed", NumSpliceJuncs,NumSpliceNotAccepted);
	if(gProcessingID)
		gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"Filtering",ePTInt32,sizeof(NumSpliceNotAccepted),"OrphanSpliceJuncts",&NumSpliceNotAccepted);
	if(NumSpliceNotAccepted)
		SortReadHits(eRSMHitMatch,false,true);
	}

return(eBSFSuccess);
}

int
CKAligner::RemoveOrphanMicroInDels(int microInDelLen) // remove any unsupported orphan microInDels
{
int Idx;
int NumInDelJuncs = 0;
int	NumInDelsAccepted = 0;
int	NumInDelsNotAccepted = 0;

tsSegJuncts *pSegJuncts = nullptr;
tsSegJuncts *pJunct;
tsSegJuncts *pNxtJunct;
tsReadHit *pCurHit = nullptr;
tsReadHit *pReadHit = nullptr;

if(microInDelLen > 0)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Filtering out any orphan microIndel reads");
	SortReadHits(eRSMHitMatch,false);
	while((pReadHit = IterReads(pReadHit))!=nullptr)
		{
		if(pReadHit->NAR != eNARAccepted || pReadHit->HitLoci.Hit.FlgInDel == 0)
			continue;
		pCurHit = pReadHit;
		NumInDelJuncs += 1;
		}
	if(NumInDelJuncs > 1)
		{
		// allocate for this number of InDel junctions, then sort and dedupe
		pSegJuncts = new tsSegJuncts[NumInDelJuncs];
		if(pSegJuncts == nullptr)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to allocate memory (%d bytes) for microInDel processing",NumInDelJuncs * sizeof(tsSegJuncts));
			Reset(false);
			return(eBSFerrMem);
			}
		pJunct = pSegJuncts;
		pReadHit = nullptr;
		while((pReadHit = IterReads(pReadHit))!=nullptr)
			{
			if(pReadHit->NAR != eNARAccepted || pReadHit->HitLoci.Hit.FlgInDel == 0)
				continue;
			pJunct->pRead = pReadHit;
			pJunct->ChromID = pReadHit->HitLoci.Hit.Seg[0].ChromID;
			pJunct->Starts = AdjEndLoci(&pReadHit->HitLoci.Hit.Seg[0]);
			pJunct->Ends = AdjStartLoci(&pReadHit->HitLoci.Hit.Seg[1]);
			pJunct->Cnt = 1;
			pJunct += 1;
			}
		// now sort ascending by chrom, start, end
		m_mtqsort.qsort(pSegJuncts,NumInDelJuncs,sizeof(tsSegJuncts),SortSegJuncts);
		// iterate and determine multiplicity of junctions
		pJunct = pSegJuncts;
		pNxtJunct = &pSegJuncts[1];
		for(Idx = 0; Idx < (NumInDelJuncs-1); Idx++,pJunct++,pNxtJunct++)
			{
			if(pJunct->ChromID == pNxtJunct->ChromID &&
			   (pJunct->Starts <= (pNxtJunct->Starts + 3) && pJunct->Starts >= (pNxtJunct->Starts - 3)) &&
			   (pJunct->Ends <= (pNxtJunct->Ends + 3) && pJunct->Ends >= (pNxtJunct->Ends - 3)))				{
				pJunct->pRead->HitLoci.Hit.FlgNonOrphan = 1;
				pNxtJunct->pRead->HitLoci.Hit.FlgNonOrphan = 1;
				}
			}
		// InDels not supported by at least two reads are to be treated as simply unaligned
		pJunct = pSegJuncts;
		for(Idx = 0; Idx < NumInDelJuncs; Idx++,pJunct++)
			{
			if(pJunct->pRead->HitLoci.Hit.FlgNonOrphan != 1)
				{
				pJunct->pRead->NAR = eNARmicroInDel;
				pJunct->pRead->NumHits = 0;	// treat as unaligned
				pJunct->pRead->LowHitInstances = 0;
				NumInDelsNotAccepted += 1;
				}
			else
				NumInDelsAccepted += 1;
			}
		delete []pSegJuncts;
		}
	else
		if(NumInDelJuncs > 0 && pCurHit != nullptr)
			{
			pCurHit->NAR = eNARmicroInDel;
			pCurHit->NumHits = 0;	// treat as unaligned
			pCurHit->LowHitInstances = 0;
			NumInDelsNotAccepted += 1;
			}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"From %d reads with putative microIndels %d orphans were removed", NumInDelJuncs,NumInDelsNotAccepted);
	if(gProcessingID)
		gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"Filtering",ePTInt32,sizeof(NumInDelsNotAccepted),"OrphanInDels",&NumInDelsNotAccepted);

	if(NumInDelsNotAccepted)
		SortReadHits(eRSMHitMatch,false,true);
	}
return(eBSFSuccess);
}



int								// -1: base not meeting constraints, 0: chromID has no constraints, 1: ChromID constrained but base accepted
CKAligner::AcceptBaseConstraint(uint32_t ChromID,			// base aligned to this chrom/sequence
								  uint32_t Loci,			// aligned to this loci
								  etSeqBase Base)		// base in read
{
int Idx;
etSeqBase TargBase;
int Rslt;
tsConstraintLoci *pConstraintLoci;

if(m_NumConstrainedChroms == 0 || m_NumConstraintLoci == 0 || m_pConstraintLoci == nullptr) // there may be no constraints!
	return(0);

Rslt = 0;		// assuming no constraints on the targeted ChromID

// current implementation not expecting too many constraints, and only a few constrained chroms, so just a simple linear search for matching constraints
pConstraintLoci = m_pConstraintLoci;
for(Idx = 0; Idx < m_NumConstraintLoci; Idx++,pConstraintLoci++)
	{
	if(pConstraintLoci->ChromID < ChromID)
		continue;
	if(pConstraintLoci->ChromID > ChromID)
		break;

	// there is at least one constraint on targeted chrom
	Rslt = 1;
	if(pConstraintLoci->StartLoci <= Loci && pConstraintLoci->EndLoci >= Loci)
		{
		if(pConstraintLoci->Constraint & 0x10)	// accept if reads base same as targets base
			{
			TargBase = m_pSfxArray->GetBase(ChromID,Loci);
			if(TargBase == Base)
				continue;
			}

		if(pConstraintLoci->Constraint & 0x01 && Base == eBaseA)
			continue;
		if(pConstraintLoci->Constraint & 0x02 && Base == eBaseC)
			continue;
		if(pConstraintLoci->Constraint & 0x04 && Base == eBaseG)
			continue;
		if(pConstraintLoci->Constraint & 0x08 && Base == eBaseT)
			continue;
		return(-1);
		}
	}
return(Rslt);
}

bool												  // true if read alignment meets any loci base constraints
CKAligner::AcceptLociConstraints(tsReadHit *pReadHit)   // read alignment to check
{
int Rslt;
int Idx;
int SeqIdx;
uint32_t ChromID;
uint8_t *pHitSeq;
uint8_t *pSeq;
uint8_t ReadSeq[cMaxSeqLen+1];
etSeqBase Base;

uint32_t StartLoci;
uint32_t EndLoci; 
uint32_t CurLoci;

if(pReadHit == nullptr || pReadHit->NAR != eNARAccepted)	  // only interested in reads which have, thus far, been accepted as being aligned
	return(true);

if(m_NumConstrainedChroms == 0 || m_NumConstraintLoci == 0 || m_pConstraintLoci == nullptr) // there may be no constraints!
	return(true);

ChromID = pReadHit->HitLoci.Hit.Seg[0].ChromID;

// check if aligned chrom has any constraint loci
for(Idx = 0; Idx < m_NumConstrainedChroms; Idx++)
	if(ChromID == m_ConstrainedChromIDs[Idx])
		break;
if(Idx == m_NumConstrainedChroms)
	return(true);


// copy read into ReadSeq masking any hi-order phred scores (in bits 4..7) and repeat mask (bit 3) out
pHitSeq = &pReadHit->Read[pReadHit->DescrLen+1];
pSeq = ReadSeq;
for(SeqIdx = 0; SeqIdx < pReadHit->ReadLen; SeqIdx++)
	*pSeq++ = *pHitSeq++ & 0x07;

// if read was aligned antisense to target then reverse complement
if(pReadHit->HitLoci.Hit.Seg[0].Strand == '-')
	CSeqTrans::ReverseComplement(pReadHit->ReadLen,ReadSeq);

StartLoci = AdjStartLoci(&pReadHit->HitLoci.Hit.Seg[0]);
EndLoci = AdjEndLoci(&pReadHit->HitLoci.Hit.Seg[0]);
SeqIdx = pReadHit->HitLoci.Hit.Seg[0].ReadOfs + pReadHit->HitLoci.Hit.Seg[0].TrimLeft;
for(CurLoci = StartLoci; CurLoci <= EndLoci; CurLoci++, SeqIdx++)
	{
	Base = ReadSeq[SeqIdx];
	Rslt = AcceptBaseConstraint(ChromID,CurLoci,Base);
	if(Rslt != 1)
		return(Rslt == -1 ? false : true);
	}

if(pReadHit->HitLoci.FlagSegs)
	{
	StartLoci = AdjStartLoci(&pReadHit->HitLoci.Hit.Seg[1]);
	EndLoci = AdjEndLoci(&pReadHit->HitLoci.Hit.Seg[1]);
	SeqIdx = pReadHit->HitLoci.Hit.Seg[1].ReadOfs + pReadHit->HitLoci.Hit.Seg[1].TrimLeft;
	for(CurLoci = StartLoci; CurLoci <= EndLoci; CurLoci++, SeqIdx++)
		{
		Base = ReadSeq[SeqIdx];
		Rslt = AcceptBaseConstraint(ChromID,CurLoci,Base);
		if(Rslt != 1)
			return(Rslt == -1 ? false : true);
		}
	}
return(true);
}

// Identify and mark with eNARLociConstrained any currently accepted alignments which violate a loci base constraint
int						// number of new alignments identified as violating a loci base constraint
CKAligner::IdentifyConstraintViolations(bool bPEread) // true if processing for PE's, false if processing for SE
{
int NumIdentified;
tsReadHit *pReadHit;
tsReadHit *pPEReadHit;
bool bIsPE2;
if(m_NumConstrainedChroms == 0 || m_NumConstraintLoci == 0 || m_pConstraintLoci == nullptr)
	return(0);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Identifying %s loci base constraint violations ...",bPEread ? "PE" : "SE");
SortReadHits(eRSMHitMatch,false);

NumIdentified = 0;
pReadHit = nullptr;
while((pReadHit = IterSortedReads(pReadHit))!=nullptr)
	{
	if(pReadHit->NAR == eNARAccepted && !AcceptLociConstraints(pReadHit))
		{
		pReadHit->NAR = eNARLociConstrained;
		pReadHit->NumHits = 0;
		pReadHit->LowHitInstances = 0;
		NumIdentified += 1;
		}

	if(bPEread && pReadHit->NAR == eNARLociConstrained)		// if PE processsing then ensure that the partner read also marked 
		{
		// locate partner read for current read
		// if current read is PE1 then PE2 will be next read, if PE2 then PE1 will be previous read
		bIsPE2 = pReadHit->PairReadID & 0x80000000 ? true : false;	// top bit set if PE2
		if(bIsPE2)	
			pPEReadHit = (tsReadHit *)((uint8_t *)pReadHit - pReadHit->PrevSizeOf);   
		else
			pPEReadHit = (tsReadHit *)((uint8_t *)pReadHit + sizeof(tsReadHit) + pReadHit->ReadLen + pReadHit->DescrLen);
		if(pPEReadHit->NAR != eNARLociConstrained)
			{
			pPEReadHit->NAR = eNARLociConstrained;
			pPEReadHit->NumHits = 0;
			pPEReadHit->LowHitInstances = 0;
			NumIdentified += 1;
			}
		}
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Identified %d %s loci base constraint violations",NumIdentified,bPEread ? "PE" : "SE");
if(NumIdentified > 0)
	SortReadHits(eRSMHitMatch,false,true);
return(NumIdentified);
}



bool					// true if chrom is accepted, false if chrom not accepted
CKAligner::AcceptThisChromID(uint32_t ChromID)
{
char szChromName[128];
bool bProcChrom = false;
int MatchesFiltOut = 0;

m_pSfxArray->GetIdentName(ChromID,sizeof(szChromName),szChromName);

AcquireSerialise();
if((bProcChrom = !m_RegExprs.MatchExcludeRegExpr(szChromName)) == true)
    bProcChrom = m_RegExprs.MatchIncludeRegExpr(szChromName);

ReleaseSerialise();

return(bProcChrom);
}


// returns
// > 0: accepted with return value the fragment (insert) size
//   0: not paired end alignment, both ends not uniquely aligned with NumHits == 1
// -1: alignment strands not consistent
// -2: aligned to different chromosomes
// -3: both ends on filtered chrom
// -4: 5' end on filtered chrom
// -5: 3' end on filtered chrom
// -6: fragment < PairMinLen
// -7: fragment > PairMaxLen
int				// > 0: accepted with return value the fragment (insert) size
CKAligner::AcceptProvPE(int PairMinLen,		// only accept paired reads with a combined sequence length of at least this
			 int PairMaxLen,				// only accept paired reads with a combined sequence length of no more than this
			 bool bPairStrand,				// accept paired ends if on same strand
			 tsReadHit *pFwdReadHit,
			 tsReadHit *pRevReadHit)
{
uint32_t FwdChrom;
uint8_t FwdStrand;
uint32_t FwdProbeHitLen;
uint32_t FwdStartLoci;
uint32_t FwdEndLoci;

uint32_t RevChrom;
uint8_t RevStrand;
uint32_t RevProbeHitLen;
uint32_t RevStartLoci;
uint32_t RevEndLoci;

// to be a pair then expecting both ends to have been aligned
if(!(pFwdReadHit->NumHits == 1 && pRevReadHit->NumHits == 1))
	return(0);

FwdChrom = pFwdReadHit->HitLoci.Hit.Seg[0].ChromID;
FwdStrand = pFwdReadHit->HitLoci.Hit.Seg[0].Strand;
FwdProbeHitLen = AdjHitLen(&pFwdReadHit->HitLoci.Hit.Seg[0]);
FwdStartLoci = AdjStartLoci(&pFwdReadHit->HitLoci.Hit.Seg[0]);
if(pFwdReadHit->HitLoci.Hit.FlgInDel || pFwdReadHit->HitLoci.Hit.FlgSplice)
	{
	FwdEndLoci = AdjEndLoci(&pFwdReadHit->HitLoci.Hit.Seg[1]);
	FwdProbeHitLen += AdjHitLen(&pFwdReadHit->HitLoci.Hit.Seg[1]);
	}
else
	FwdEndLoci = AdjEndLoci(&pFwdReadHit->HitLoci.Hit.Seg[0]);
RevChrom = pRevReadHit->HitLoci.Hit.Seg[0].ChromID;
RevStrand = pRevReadHit->HitLoci.Hit.Seg[0].Strand;
RevProbeHitLen = AdjHitLen(&pRevReadHit->HitLoci.Hit.Seg[0]);
RevStartLoci = AdjStartLoci(&pRevReadHit->HitLoci.Hit.Seg[0]);
if(pRevReadHit->HitLoci.Hit.FlgInDel || pRevReadHit->HitLoci.Hit.FlgSplice)
	{
	RevEndLoci = AdjEndLoci(&pRevReadHit->HitLoci.Hit.Seg[1]);
	RevProbeHitLen += AdjHitLen(&pRevReadHit->HitLoci.Hit.Seg[1]);
	}
else
	RevEndLoci = AdjEndLoci(&pRevReadHit->HitLoci.Hit.Seg[0]);

bool bFwdChrom = AcceptThisChromID(FwdChrom);
if(FwdChrom != RevChrom)									
	{
	bool bRevChrom = AcceptThisChromID(RevChrom);
	if(bRevChrom == true && bFwdChrom == true)
		return(-2);
	if(bRevChrom == false && bFwdChrom == false)
		return(-3);
	if(bFwdChrom == false)
		return(-4);
	return(-5);
	}
else
	if(bFwdChrom == false)
		return(-3);

return(PEInsertSize(PairMinLen,PairMaxLen,bPairStrand,FwdStrand,FwdStartLoci,FwdEndLoci,RevStrand,RevStartLoci,RevEndLoci));
}

// calculates PE alignment insert size
// returns
// > 0: accepted with return value the fragment (insert) size
//   0: not paired end alignment, both ends not uniquely aligned with NumHits == 1
// -1: alignment strands not consistent
// -2: aligned to different chromosomes
// -3: both ends on filtered chrom
// -4: 5' end on filtered chrom
// -5: 3' end on filtered chrom
// -6: fragment < PairMinLen
// -7: fragment > PairMaxLen
int										// returned PE insert size, <= 0 if errors
CKAligner::PEInsertSize(int PairMinLen,	// only accept paired reads with a combined sequence length of at least this
			 int PairMaxLen,			// only accept paired reads with a combined sequence length of no more than this
			 bool bPairStrand,			// accept paired ends if on same strand		
			 uint8_t PE1Strand,			// PE1 aligned on to this strand
			 uint32_t PE1StartLoci,		// PE read starts at this loci
			 uint32_t PE1EndLoci,			// PE1 read ends at this loci
			 uint8_t PE2Strand,			// PE2 aligned on to this strand
			 uint32_t PE2StartLoci,		// PE2 read starts at this loci
			 uint32_t PE2EndLoci)			// PE2 read ends at this loci
{
int SeqFragLen;

if((bPairStrand && PE1Strand != PE2Strand) ||				
	(!bPairStrand && PE1Strand == PE2Strand))
		return(-1);

// if processing for circulars then 
//		if fwd is on '+' then distance = FwdEndLoci - RevStartLoci
//		if fwd is on '-' then distance = RevStartLoci - FwdEndLoci
// else
//		if fwd is on '+' then distance = RevEndLoci - FwdStartLoci
//		if fwd is on '-' then distance = FwdStartLoci - RevStartLoci
if(m_bPEcircularised)
	{
	if(PE1Strand == '+')
		SeqFragLen = 1 + (int)PE1EndLoci - (int)PE2StartLoci;
	else
		SeqFragLen = 1 + (int)PE2StartLoci - (int)PE1EndLoci;
	}
else
	SeqFragLen = 1 + max(PE1EndLoci, PE2EndLoci) - min(PE1StartLoci, PE2StartLoci);

if(SeqFragLen < 0)			// treat as inconsistent strand if fragment length is negative 
	return(-1);

if(SeqFragLen < PairMinLen)
	return(-6);

if(SeqFragLen > PairMaxLen)
	return(-7);

// can accept
return(SeqFragLen);
}

#ifdef _WIN32
unsigned __stdcall KProcessPairedEndsThread(void * pThreadPars)
#else
void *KProcessPairedEndsThread(void * pThreadPars)
#endif
{
	int Rslt;
	tsPEThreadPars *pPars = (tsPEThreadPars *)pThreadPars;			// makes it easier not having to deal with casts!
	CKAligner *pKAligner = (CKAligner *)pPars->pThis;
	Rslt = pKAligner->ProcessPairedEnds(pPars);
	pPars->Rslt = Rslt;
#ifdef _WIN32
	_endthreadex(0);
	return(eBSFSuccess);
#else
pthread_exit(&pPars->Rslt);
#endif
}

//ProcessPairedEnds
// Matches read pairs and if one read is aligned but partner is not due to muiltihit loci then
// attempts to find a unique loci for that multihit read within the expected approximate pair distance range
// If user has optionally requested then will accept reads which are orphaned (other PE not aligned) or PEs where PE1 and PE2 align to separate chroms/contigs 
int
CKAligner::ProcessPairedEnds(etPEproc PEproc, // paired reads alignment processing mode
				  int MinEditDist, // accepted alignments must be at least this Hamming away from other putative alignments
				  int PairMinLen,  // only accept paired reads with a combined sequence length of at least this
				  int PairMaxLen,  // only accept paired reads with a combined sequence length of no more than this
				  bool bPairStrand,	// accept paired ends if on same strand
				  int MaxSubs)	   // aligned reads can have at most this many subs per 100bp
{
int Rslt;
uint32_t PairReadIdx;

int UnderLenPairs = 0;
int OverLenPairs = 0;

int PartnerUnpaired = 0;
int PartnerPaired = 0;
int UnalignedPairs = 0;
int UniquePairCnt = 0;
int NumFilteredByChrom = 0;

int AcceptedNumPaired = 0;
int AcceptedNumSE = 0;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Generating paired reads index over %d paired reads", m_NumReadsLoaded/2);
SortReadHits(eRSMPairReadID,false,true);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Starting to associate Paired End reads to be within insert size range ...");

if(m_pLenDist != nullptr)
	{
	delete []m_pLenDist;
	m_pLenDist = nullptr;
	}
if((m_pLenDist = new int[cPairMaxLen+1])==nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate %d bytes memory for paired read sequence length distributions",
								(int)sizeof(int) * (cPairMaxLen+1));
	Reset(false);
	return(eBSFerrMem);
	}
memset(m_pLenDist,0,sizeof(int) * (cPairMaxLen+1));

time_t Started = time(0);
uint32_t PrevPairReadIdx = 0;
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processed putative 0 pairs, accepted 0");

// split number of pairs to be processed amongst number of worker threads


int ThreadIdx;
int NumThreadsUsed;
tsPEThreadPars WorkerThreads[cMaxWorkerThreads];
memset(WorkerThreads, 0, sizeof(WorkerThreads));
uint32_t NumPairsThisThread;
PairReadIdx = 0;
NumThreadsUsed = 0;
for (ThreadIdx = 0; ThreadIdx < m_NumThreads; ThreadIdx++, NumThreadsUsed += 1)
	{
	WorkerThreads[ThreadIdx].ThreadIdx = ThreadIdx + 1;
	WorkerThreads[ThreadIdx].pThis = this;
	WorkerThreads[ThreadIdx].bPairStrand = bPairStrand;
	WorkerThreads[ThreadIdx].MaxSubs = MaxSubs;
	WorkerThreads[ThreadIdx].MinEditDist = MinEditDist;
	WorkerThreads[ThreadIdx].StartPairIdx = PairReadIdx;
	NumPairsThisThread = ((m_NumReadsLoaded / 2)- PairReadIdx) / (m_NumThreads - ThreadIdx);
	WorkerThreads[ThreadIdx].NumPairsToProcess = NumPairsThisThread;
	PairReadIdx += NumPairsThisThread;
	WorkerThreads[ThreadIdx].PairMaxLen = PairMaxLen;
	WorkerThreads[ThreadIdx].PairMinLen = PairMinLen;
	WorkerThreads[ThreadIdx].PEproc = PEproc;
#ifdef _WIN32
	WorkerThreads[ThreadIdx].threadHandle = (HANDLE)_beginthreadex(nullptr, 0x0fffff, KProcessPairedEndsThread, &WorkerThreads[ThreadIdx], 0, &WorkerThreads[ThreadIdx].threadID);
#else
	WorkerThreads[ThreadIdx].threadRslt = pthread_create(&WorkerThreads[ThreadIdx].threadID, nullptr, KProcessPairedEndsThread, &WorkerThreads[ThreadIdx]);
#endif
	if((WorkerThreads[ThreadIdx].StartPairIdx + WorkerThreads[ThreadIdx].NumPairsToProcess) == m_NumReadsLoaded / 2)
		{
		NumThreadsUsed = ThreadIdx + 1;
		break;
		}
	}

// allow threads a few seconds to startup
#ifdef _WIN32
Sleep(5000);
#else
sleep(5);
#endif

uint32_t ReportProgressSecs;
ReportProgressSecs = 60;

// wait for all threads to have completed
for (ThreadIdx = 0; ThreadIdx < NumThreadsUsed; ThreadIdx++)
	{
#ifdef _WIN32
	while (WAIT_TIMEOUT == WaitForSingleObject(WorkerThreads[ThreadIdx].threadHandle, (DWORD)ReportProgressSecs * 1000))
		{
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Progress: Still associating Paired Ends ...");
		}
	CloseHandle(WorkerThreads[ThreadIdx].threadHandle);
#else
	struct timespec ts;
	int JoinRlt;
	clock_gettime(CLOCK_REALTIME, &ts);
	ts.tv_sec += ReportProgressSecs;
	while ((JoinRlt = pthread_timedjoin_np(WorkerThreads[ThreadIdx].threadID, nullptr, &ts)) != 0)
		{
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Progress: Still associating Paired Ends ...");
		ts.tv_sec += ReportProgressSecs;
		}

#endif
	AcceptedNumPaired += WorkerThreads[ThreadIdx].AcceptedNumPaired;
	PartnerPaired += WorkerThreads[ThreadIdx].PartnerPaired;
	PartnerUnpaired += WorkerThreads[ThreadIdx].PartnerUnpaired;

	UnderLenPairs += WorkerThreads[ThreadIdx].UnderLenPairs;
	OverLenPairs += WorkerThreads[ThreadIdx].OverLenPairs;
	NumFilteredByChrom += WorkerThreads[ThreadIdx].NumFilteredByChrom;
	UnalignedPairs += WorkerThreads[ThreadIdx].UnalignedPairs;
	AcceptedNumSE += WorkerThreads[ThreadIdx].AcceptedNumSE;
	PartnerUnpaired += WorkerThreads[ThreadIdx].PartnerUnpaired;
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Completed association of Paired End reads from %u pairs, accepted %u pairs", m_NumReadsLoaded /2,AcceptedNumPaired);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"From %d Paired End pairs there were %d accepted (of which %d pairs were from recovered orphans)",
			m_NumReadsLoaded / 2,AcceptedNumPaired,PartnerPaired);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"%d Paired End pairs unrecoverable as still orphan partnered",PartnerUnpaired - PartnerPaired);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"%d Paired End pairs were under length, %d over length, not accepted as being paired",UnderLenPairs,OverLenPairs);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"%d Paired End aligned pairs were filtered out by chromosome",NumFilteredByChrom);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"%d Paired End pairs have neither end uniquely aligned",UnalignedPairs);
if(PEproc == ePEuniqueSE || PEproc == ePEorphanSE)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"%d Paired End reads were unable to be associated with partner read and accepted as if SE aligned",AcceptedNumSE);

if(gProcessingID > 0)
	{
	uint32_t NumPairsLoaded;
	NumPairsLoaded = m_NumReadsLoaded/2;
	gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"PEProcessing",ePTUint32,sizeof(NumPairsLoaded),"Loaded",&NumPairsLoaded);
	gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"PEProcessing",ePTUint32,sizeof(AcceptedNumPaired),"Accepted",&AcceptedNumPaired);
	gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"PEProcessing",ePTUint32,sizeof(PartnerPaired),"RecoveredOrphans",&PartnerPaired);
	gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"PEProcessing",ePTUint32,sizeof(UnderLenPairs),"UnderInsertSize",&UnderLenPairs);
	gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"PEProcessing",ePTUint32,sizeof(OverLenPairs),"OverInsertSize",&OverLenPairs);
	gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"PEProcessing",ePTUint32,sizeof(NumFilteredByChrom),"FilteredByChrom",&NumFilteredByChrom);
	gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"PEProcessing",ePTUint32,sizeof(UnalignedPairs),"NonUniqueEnds",&UnalignedPairs);
	gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"PEProcessing",ePTUint32,sizeof(AcceptedNumSE),"AcceptedNumSE",&AcceptedNumSE);
	}

if(m_hStatsFile != -1)
	{
	int hGlobalPEInsertDist;
	char szGlobalPEInsertDist[_MAX_PATH];
	char szLineBuff[(cMaxReadLen + 1000)*2];
	int BuffIdx = 0;

	CUtility::AppendFileNameSuffix(szGlobalPEInsertDist, m_pszStatsFile, (char*)".GlobalPEInsertDist.csv", '.');
#ifdef _WIN32
	hGlobalPEInsertDist = open(szGlobalPEInsertDist, (O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC), (_S_IREAD | _S_IWRITE));
#else
	if ((hGlobalPEInsertDist = open(szGlobalPEInsertDist, O_WRONLY | O_CREAT, S_IREAD | S_IWRITE)) != -1)
		if (ftruncate(hGlobalPEInsertDist, 0) != 0)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to truncate global PE insert size distribution file '%s' - %s", szGlobalPEInsertDist, strerror(errno));
			return(eBSFerrCreateFile);
			}
#endif

	if (hGlobalPEInsertDist < 0)
		{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to truncate global PE insert size distribution file '%s'", szGlobalPEInsertDist);
			return(eBSFerrCreateFile);
		}

	for(PairReadIdx = 0; PairReadIdx <= (uint32_t)cPairMaxLen; PairReadIdx++)
		{
		BuffIdx += sprintf(&szLineBuff[BuffIdx],"%d,%d\n",PairReadIdx,m_pLenDist[PairReadIdx]);
		if(BuffIdx + 100 > sizeof(szLineBuff))
			{
			if(!CUtility::RetryWrites(hGlobalPEInsertDist,szLineBuff,BuffIdx))
				{
				gDiagnostics.DiagOut(eDLInfo,gszProcName,"RetryWrites() error");
				return(eBSFerrWrite);
				}
			BuffIdx = 0;
			}
		}
	if(BuffIdx)
		{
		if(!CUtility::RetryWrites(hGlobalPEInsertDist,szLineBuff,BuffIdx))
			{
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"RetryWrites() error");
			return(eBSFerrWrite);
			}
		BuffIdx = 0;
		}
#ifdef _WIN32
	_commit(hGlobalPEInsertDist);
#else
	fsync(hGlobalPEInsertDist);
#endif
	close(hGlobalPEInsertDist);
	}

if((Rslt=SortReadHits(eRSMHitMatch,false,true)) < eBSFSuccess)
	{
	Reset(false);
	return(Rslt);
	}
return(eBSFSuccess);
}




int
CKAligner::ProcessPairedEnds(tsPEThreadPars *pPars)	   // paired reads parameters
{
	int Rslt;
	uint32_t PrevChromID;
	int SeqFragLen;
	uint32_t PairReadIdx;
	tsReadHit *pFwdReadHit;
	tsReadHit *pRevReadHit;
	uint32_t OrphStartLoci;
	uint32_t OrphEndLoci;
	bool bFwdUnaligned;
	bool bRevUnaligned;
	int ProbeLen;
	int MaxTotMM;
	int CoreLen;
	int MaxNumSlides;
	int CoreDelta;

	int NumPaired = 0;
	int NegPairs = 0;
	int UnderLenPairs = 0;
	int OverLenPairs = 0;
	int LongestSeqFragLen = 0;
	int OverPairMaxLen = 0;

	bool bPartnerPaired;
	int PartnerUnpaired = 0;
	int PartnerPaired = 0;
	int UnalignedPairs = 0;
	int UniquePairCnt = 0;
	int NumFilteredByChrom = 0;

	int AcceptedNumPaired = 0;
	int AcceptedNumSE = 0;

	uint32_t CurPairIdx;

	uint8_t ReadSeq[cMaxSeqLen + 1];

	tsHitLoci HitLoci;
	uint8_t *pHitSeq;
	uint8_t *pSeq;
	uint32_t SeqIdx;
	bool b3primeExtend;
	bool bAntisense;
	PrevChromID = -1;

	uint32_t PrevPairReadIdx = 0;
	for (CurPairIdx = pPars->StartPairIdx; CurPairIdx < (pPars->StartPairIdx + pPars->NumPairsToProcess); CurPairIdx += 1)
		{
		PairReadIdx = CurPairIdx * 2;
		pFwdReadHit = m_ppReadHitsIdx[PairReadIdx];
		pRevReadHit = m_ppReadHitsIdx[PairReadIdx + 1];
		pFwdReadHit->FlgPEAligned = 0;
		pRevReadHit->FlgPEAligned = 0;
		bFwdUnaligned = (pFwdReadHit->NAR == eNARNs || pFwdReadHit->NAR == eNARNoHit || pFwdReadHit->NAR == eNARUnaligned);
		bRevUnaligned = (pRevReadHit->NAR == eNARNs || pRevReadHit->NAR == eNARNoHit || pRevReadHit->NAR == eNARUnaligned);
		// at least one of the ends must have been accepted as being aligned in order to process as PE
		if (!(pFwdReadHit->NAR == eNARAccepted || pRevReadHit->NAR == eNARAccepted))
			{
			UnalignedPairs += 1;
			continue;
			}

		// at least one end had a unique hit accepted, but if other end had no potential hits at all then
		// can't recover into an aligned PE with both ends aligned
		if (pPars->PEproc == ePEunique && (bFwdUnaligned || bRevUnaligned))
			{
			pFwdReadHit->NumHits = 0;
			pRevReadHit->NumHits = 0;
			pFwdReadHit->LowHitInstances = 0;
			pRevReadHit->LowHitInstances = 0;
			if (pFwdReadHit->NAR == eNARAccepted)
				pFwdReadHit->NAR = eNARPENoHit;
			if (pRevReadHit->NAR == eNARAccepted)
				pRevReadHit->NAR = eNARPENoHit;
			PartnerUnpaired += 1;
			continue;
			}

		// even if both ends have been accepted as aligned it could be that these alignments can't be accepted as a PE - may be to different chroms, or overlength etc
		if (pFwdReadHit->NAR == eNARAccepted && pRevReadHit->NAR == eNARAccepted)
			{
			// both ends were accepted as being aligned but can these be accepted as PE within allowed insert size range ...
			SeqFragLen = AcceptProvPE(pPars->PairMinLen, pPars->PairMaxLen, pPars->bPairStrand, pFwdReadHit, pRevReadHit);
			if (SeqFragLen > 0)
				{
				// this paired end alignment has been accepted
				pFwdReadHit->FlgPEAligned = 1;
				pRevReadHit->FlgPEAligned = 1;
				LongestSeqFragLen = max(LongestSeqFragLen, SeqFragLen);
				AcquireSerialise();
				m_pLenDist[SeqFragLen] += 1;
				ReleaseSerialise();
				AcceptedNumPaired += 1;
				continue;	// onto next pair
				}

			// although ends aligned, am unable to accept as a PE alignment
			switch (SeqFragLen) {
			case -1:			// alignment strands not consistent
				pFwdReadHit->NAR = eNARPEStrand;
				pRevReadHit->NAR = eNARPEStrand;
				break;

			case -2:			// aligned to different chromosomes
				pFwdReadHit->NAR = eNARPEChrom;
				pRevReadHit->NAR = eNARPEChrom;
				break;

			case -3:			// both ends were to a filtered chrom so can't use either end as an anchor
				NumFilteredByChrom += 1;
				pFwdReadHit->NumHits = 0;
				pRevReadHit->NumHits = 0;
				pFwdReadHit->LowHitInstances = 0;
				pRevReadHit->LowHitInstances = 0;
				pFwdReadHit->NAR = eNARChromFilt;
				pRevReadHit->NAR = eNARChromFilt;
				continue;	// try next pair

			case -4:		// 5' end to filtered chrom so can't be used as an anchor
				pFwdReadHit->NAR = eNARChromFilt;
				pFwdReadHit->LowHitInstances = 0;
				pFwdReadHit->NumHits = 0;
				break;

			case -5:		// 3' end to filtered chrom so can't be used as an anchor
				pRevReadHit->NAR = eNARChromFilt;
				pRevReadHit->NumHits = 0;
				pRevReadHit->LowHitInstances = 0;
				break;

			case -6:		// under min insert size
				pFwdReadHit->NAR = eNARPEInsertMin;
				pRevReadHit->NAR = eNARPEInsertMin;
				break;

			case -7:		// over max insert size
				pFwdReadHit->NAR = eNARPEInsertMax;
				pRevReadHit->NAR = eNARPEInsertMax;
				break;
			}

			// if not allowed to orphan recover or treat as SE alignments then try next pair of reads
			if (pPars->PEproc == ePEunique)
			{
				pFwdReadHit->NumHits = 0;
				pRevReadHit->NumHits = 0;
				pFwdReadHit->LowHitInstances = 0;
				pRevReadHit->LowHitInstances = 0;
				if (pFwdReadHit->NAR == eNARAccepted)
					pFwdReadHit->NAR = eNARPENoHit;
				if (pRevReadHit->NAR == eNARAccepted)
					pRevReadHit->NAR = eNARPENoHit;
				PartnerUnpaired += 1;
				continue;
			}
		}

		// at least one end was uniquely aligning although not accepted as a PE
		PartnerUnpaired += 1;
		if (pPars->PEproc == ePEorphan || pPars->PEproc == ePEorphanSE)		// allowed to try and recover?
		{
			if (pFwdReadHit->NumHits == 1 && !bRevUnaligned)
			{
				if (AcceptThisChromID(pFwdReadHit->HitLoci.Hit.Seg[0].ChromID))
				{
					// first try using the 5' alignment as an anchor, if that doesn't provide a pair then will later try using the 3' as the anchor
					bPartnerPaired = false;
					pHitSeq = &pRevReadHit->Read[pRevReadHit->DescrLen + 1];   // forward read uniquely aligned so will try to align the other read within insert distance
					pSeq = ReadSeq;
					for (SeqIdx = 0; SeqIdx < pRevReadHit->ReadLen; SeqIdx++)
						*pSeq++ = *pHitSeq++ & 0x07;

					// AlignPartnerRead
					// Have been able to unquely align one read out of a pair, now need to align the other read
					// if not bPairStrand
					//		if PE1 was to sense strand then expect PE2 on the antisense strand downstream towards the 3' end of targeted chrom:		b3primeExtend=true,bAntisense=true
					//		if PE1 was to antisense strand then expect PE2 on the sense strand upstream towards the 5' end of targeted chrom:		b3primeExtend=false,bAntisense=false
					// if bPairStrand
					//		if PE1 was to sense strand then expect PE2 on the sense strand downstream towards the 3' end of targeted chrom:			b3primeExtend=true,bAntisense=false
					//		if PE1 was to antisense strand then expect PE2 on the antisense strand upstream towards the 5' end of targeted chrom:	b3primeExtend=false,bAntisense=true
					b3primeExtend = pFwdReadHit->HitLoci.Hit.Seg[0].Strand == '+' ? true : false;
					if (pPars->bPairStrand)
						bAntisense = pFwdReadHit->HitLoci.Hit.Seg[0].Strand == '+' ? false : true;
					else
						bAntisense = pFwdReadHit->HitLoci.Hit.Seg[0].Strand == '+' ? true : false;
					if (m_bPEcircularised)
						b3primeExtend = !b3primeExtend;
					
					bool bCheckMe;
					if(pFwdReadHit->ReadLen == pFwdReadHit->HitLoci.Hit.Seg[0].MatchLen)
						bCheckMe = true;
					else
						bCheckMe = false;
					OrphStartLoci = AdjStartLoci(&pFwdReadHit->HitLoci.Hit.Seg[0]);
					OrphEndLoci = AdjEndLoci(&pFwdReadHit->HitLoci.Hit.Seg[0]);

					// note: MaxSubs is specified by user as being per 100bp of read length, e.g. if user specified '-s5' and a read is 200bp then 
					// 10 mismatches will be allowed for that specific read
					ProbeLen = pRevReadHit->ReadLen;
					MaxTotMM = m_MaxSubs == 0 ? 0 : max(1, (int)(0.5 + ((ProbeLen - 1) * m_MaxSubs) / 100.0));

					if (MaxTotMM > cMaxTotAllowedSubs)		// irrespective of length allow at most this many subs
						MaxTotMM = cMaxTotAllowedSubs;

					// The window core length is set to be read length / (subs+1) for minimum Hamming difference of 1, and
					// to be read length / (subs+2) for minimum Hamming difference of 2
					// The window core length is clamped to be at least m_MinCoreLen
					CoreLen = max(m_MinCoreLen, ProbeLen / (m_MinEditDist == 1 ? MaxTotMM + 1 : MaxTotMM + 2));
					MaxNumSlides = max(1, ((m_MaxNumSlides * ProbeLen) + 99) / 100);
					CoreDelta = max(ProbeLen / m_MaxNumSlides - 1, CoreLen);

					Rslt = m_pSfxArray->AlignPairedRead(b3primeExtend, bAntisense,
						pFwdReadHit->HitLoci.Hit.Seg[0].ChromID,	  // accepted aligned read was on this chromosome
						OrphStartLoci,			// accepted aligned read started at this loci
						OrphEndLoci,			// and ending at this loci
						pPars->PairMinLen,				// expecting partner to align at least this distance away from accepted aligned read
						pPars->PairMaxLen,				// but no more than this distance away
						pPars->MaxSubs,				// any accepted alignment can have at most this many mismatches
						pPars->MinEditDist,			// and must be at least this Hamming away from the next best putative alignment
						ProbeLen,		  // length of read excluding any eBaseEOS
						m_MinChimericLen,	// minimum chimeric length as a percentage (0 to disable, otherwise 50..99) of probe sequence
						CoreLen,			// core window length, 0 to disable
						CoreDelta,			// core window offset increment (1..n)
						MaxNumSlides,		// max number of times to slide core
						ReadSeq,			// pts to 5' start of read sequence
						&HitLoci);			// where to return any paired read alignment loci

					if (Rslt == 1)
						{
						SeqFragLen = PEInsertSize(pPars->PairMinLen, pPars->PairMaxLen, pPars->bPairStrand, pFwdReadHit->HitLoci.Hit.Seg[0].Strand, OrphStartLoci, OrphEndLoci, HitLoci.Seg[0].Strand, AdjStartLoci(&HitLoci.Seg[0]), AdjEndLoci(&HitLoci.Seg[0]));
						if (SeqFragLen <= 0)
							Rslt = 0;
						}

					if (Rslt == 1)
						{
						// with 5' anchor was able to find an alignment within the min/max insert size and it is known chrom accepted
						pRevReadHit->HitLoci.Hit = HitLoci;
						pRevReadHit->NumHits = 1;
						pRevReadHit->LowMMCnt = HitLoci.Seg[0].Mismatches;
						pRevReadHit->LowHitInstances = 1;
						// this paired end alignment has been accepted
						pFwdReadHit->FlgPEAligned = 1;
						pRevReadHit->FlgPEAligned = 1;
						pFwdReadHit->NAR = eNARAccepted;
						pRevReadHit->NAR = eNARAccepted;
						LongestSeqFragLen = max(LongestSeqFragLen, SeqFragLen);
						m_pLenDist[SeqFragLen] += 1;
						AcceptedNumPaired += 1;
						PartnerPaired += 1;
						continue;	// try next pair
						}
					}
				else
					if (pFwdReadHit->NAR == eNARAccepted)
					{
						pFwdReadHit->NumHits = 0;
						pFwdReadHit->LowHitInstances = 0;
						pFwdReadHit->NAR = eNARChromFilt;
					}
			}


			if (pRevReadHit->NumHits == 1 && !bFwdUnaligned)
			{
				if (AcceptThisChromID(pRevReadHit->HitLoci.Hit.Seg[0].ChromID))
				{
					pHitSeq = &pFwdReadHit->Read[pFwdReadHit->DescrLen + 1];    // reverse read uniquely aligned so will try to align the other read within insert distance
					pSeq = ReadSeq;
					for (SeqIdx = 0; SeqIdx < pFwdReadHit->ReadLen; SeqIdx++)
						*pSeq++ = *pHitSeq++ & 0x07;
					// AlignPartnerRead
					// Have been able to unquely align one read out of a pair, now need to align the other read
					// if not bPairStrand
					//		if PE2 was to sense strand then expect PE1 on the antisense strand downstream towards the 3' end of targeted chrom:		b3primeExtend=true,bAntisense=true
					//		if PE2 was to antisense strand then expect PE1 on the sense strand upstream towards the 5' end of targeted chrom:		b3primeExtend=false,bAntisense=false
					// if bPairStrand
					//		if PE2 was to sense strand then expect PE1 on the sense strand upstream towards the 5' end of targeted chrom:			b3primeExtend=false,bAntisense=false
					//		if PE2 was to antisense strand then expect PE2 on the antisense strand downstream towards the 3' end of targeted chrom:	b3primeExtend=true,bAntisense=true
					b3primeExtend = pRevReadHit->HitLoci.Hit.Seg[0].Strand == '+' ? true : false;
					bAntisense = pRevReadHit->HitLoci.Hit.Seg[0].Strand == '+' ? true : false;
					if (pPars->bPairStrand)
					{
						b3primeExtend = !b3primeExtend;
						bAntisense = !bAntisense;
					}
					if (m_bPEcircularised)
						b3primeExtend = !b3primeExtend;

					bool bCheckMe;
					if(pRevReadHit->ReadLen == pRevReadHit->HitLoci.Hit.Seg[0].MatchLen)
						bCheckMe = true;
					else
						bCheckMe = false;
					OrphStartLoci = AdjStartLoci(&pRevReadHit->HitLoci.Hit.Seg[0]);
					OrphEndLoci = AdjEndLoci(&pRevReadHit->HitLoci.Hit.Seg[0]);

					// note: MaxSubs is specified by user as being per 100bp of read length, e.g. if user specified '-s5' and a read is 200bp then 
					// 10 mismatches will be allowed for that specific read
					ProbeLen = pFwdReadHit->ReadLen;
					MaxTotMM = m_MaxSubs == 0 ? 0 : max(1, (int)(0.5 + ((ProbeLen - 1) * m_MaxSubs) / 100.0));

					if (MaxTotMM > cMaxTotAllowedSubs)		// irrespective of length allow at most this many subs
						MaxTotMM = cMaxTotAllowedSubs;

					// The window core length is set to be read length / (subs+1) for minimum Hamming difference of 1, and
					// to be read length / (subs+2) for minimum Hamming difference of 2
					// The window core length is clamped to be at least m_MinCoreLen
					CoreLen = max(m_MinCoreLen, pFwdReadHit->ReadLen / (m_MinEditDist == 1 ? MaxTotMM + 1 : MaxTotMM + 2));
					MaxNumSlides = max(1, ((m_MaxNumSlides * ProbeLen) + 99) / 100);
					CoreDelta = max(ProbeLen / m_MaxNumSlides - 1, CoreLen);


					Rslt = m_pSfxArray->AlignPairedRead(b3primeExtend, bAntisense,
						pRevReadHit->HitLoci.Hit.Seg[0].ChromID,	  // accepted aligned read was on this chromosome
						OrphStartLoci,			// accepted aligned read started at this loci
						OrphEndLoci,		  // and ending at this loci
						pPars->PairMinLen,			// expecting partner to align at least this distance away from accepted aligned read
						pPars->PairMaxLen,		// but no more than this distance away
						pPars->MaxSubs,				// any accepted alignment can have at most this many mismatches
						pPars->MinEditDist,			// and must be at least this Hamming away from the next best putative alignment
						ProbeLen,		  // length of read excluding any eBaseEOS
						m_MinChimericLen,	// minimum chimeric length as a percentage (0 to disable, otherwise 50..99) of probe sequence
						CoreLen,			// core window length, 0 to disable
						CoreDelta,			// core window offset increment (1..n)
						MaxNumSlides,		// max number of times to slide core
						ReadSeq,	  // pts to 5' start of read sequence
						&HitLoci);	  // where to return any paired read alignment loci

					if (Rslt == 1)
						{
						SeqFragLen = PEInsertSize(pPars->PairMinLen, pPars->PairMaxLen, pPars->bPairStrand, HitLoci.Seg[0].Strand, AdjStartLoci(&HitLoci.Seg[0]), AdjEndLoci(&HitLoci.Seg[0]), pRevReadHit->HitLoci.Hit.Seg[0].Strand, OrphStartLoci, OrphEndLoci);
						if (SeqFragLen <= 0)
							Rslt = 0;
						}

					if (Rslt == 1)
					{
						// with 3' anchor was able to find an alignment within the min/max insert size
						pFwdReadHit->HitLoci.Hit = HitLoci;
						pFwdReadHit->LowMMCnt = HitLoci.Seg[0].Mismatches;
						pFwdReadHit->NumHits = 1;
						pFwdReadHit->LowHitInstances = 1;
						// this paired end alignment has been accepted
						pFwdReadHit->FlgPEAligned = 1;
						pRevReadHit->FlgPEAligned = 1;
						pFwdReadHit->NAR = eNARAccepted;
						pRevReadHit->NAR = eNARAccepted;
						LongestSeqFragLen = max(LongestSeqFragLen, SeqFragLen);
						AcquireSerialise();
						m_pLenDist[SeqFragLen] += 1;
						ReleaseSerialise();
						AcceptedNumPaired += 1;
						PartnerPaired += 1;
						continue;	// try next pair
					}
				}
				else
					if (pRevReadHit->NAR == eNARAccepted)
					{
						pFwdReadHit->NumHits = 0;
						pFwdReadHit->LowHitInstances = 0;
						pFwdReadHit->NAR = eNARChromFilt;
					}
			}
		}

		// unable to accept as PE
		if (pFwdReadHit->NAR == eNARChromFilt || pRevReadHit->NAR == eNARChromFilt)
			NumFilteredByChrom += 1;
		if (pFwdReadHit->NAR == eNARPEInsertMin || pRevReadHit->NAR == eNARPEInsertMin)
			UnderLenPairs += 1;
		if (pFwdReadHit->NAR == eNARPEInsertMax || pRevReadHit->NAR == eNARPEInsertMax)
			OverLenPairs += 1;

		if (!(pPars->PEproc == ePEorphanSE || pPars->PEproc == ePEuniqueSE))
		{
			pFwdReadHit->NumHits = 0;
			pFwdReadHit->LowHitInstances = 0;
			pRevReadHit->NumHits = 0;
			pRevReadHit->LowHitInstances = 0;
			if (pFwdReadHit->NAR == eNARAccepted)
				pFwdReadHit->NAR = eNARPENoHit;
			if (pRevReadHit->NAR == eNARAccepted)
				pRevReadHit->NAR = eNARPENoHit;
			continue;
		}

		// allowed to accept as being SE if was able to uniquely align
		bool bChromFilt;
		if (pFwdReadHit->NumHits == 1)
			bChromFilt = AcceptThisChromID(pFwdReadHit->HitLoci.Hit.Seg[0].ChromID);
		else
			bChromFilt = false;
		if (pFwdReadHit->NumHits != 1 || !bChromFilt)
		{
			pFwdReadHit->NumHits = 0;
			pFwdReadHit->LowHitInstances = 0;
			if (pFwdReadHit->NAR == eNARAccepted)
				pFwdReadHit->NAR = bChromFilt ? eNARChromFilt : eNARPEUnalign;
		}
		else
		{
			pFwdReadHit->NAR = eNARAccepted;
			AcceptedNumSE += 1;
		}

		if (pRevReadHit->NumHits == 1)
			bChromFilt = AcceptThisChromID(pRevReadHit->HitLoci.Hit.Seg[0].ChromID);
		else
			bChromFilt = false;

		if (pRevReadHit->NumHits != 1 || !bChromFilt)
		{
			pRevReadHit->NumHits = 0;
			pRevReadHit->LowHitInstances = 0;
			if (pRevReadHit->NAR == eNARAccepted)
				pRevReadHit->NAR = bChromFilt ? eNARChromFilt : eNARPEUnalign;
		}
		else
		{
			pRevReadHit->NAR = eNARAccepted;
			AcceptedNumSE += 1;
		}
	}
pPars->UnalignedPairs = UnalignedPairs;
pPars->AcceptedNumPaired = AcceptedNumPaired;
pPars->AcceptedNumSE = AcceptedNumSE;
pPars->PartnerPaired = PartnerPaired;
pPars->PartnerUnpaired = PartnerUnpaired;
pPars->NumFilteredByChrom = NumFilteredByChrom;
pPars->UnderLenPairs = UnderLenPairs;
pPars->OverLenPairs = OverLenPairs;
pPars->Rslt = 0;
return(0);
}


int
CKAligner::ReportAlignStats(void)		// report basic alignment statistics
{
int Rslt;
uint32_t NumUniques = 0;
uint32_t NumPlusHits = 0;
uint32_t NumNoMatches = 0;
uint32_t NumChimeric = 0;
uint32_t NumMultiMatches = 0;
uint32_t NumHamming = 0;
tsReadHit *pReadHit = nullptr;
bool bSimReads = false;
uint32_t NumReads1EdgeAligned = 0;
uint32_t NumReads2EdgeAligned = 0;
uint32_t NumReadsMisaligned = 0;
uint32_t NumIndels = 0;
uint32_t NumSpliced = 0;
uint32_t PrevTargEntry = 0;
uint32_t NumTrimmed = 0;

uint32_t NARUnaligned = 0;				// read has yet to be aligned
uint32_t NARAccepted = 0;					// read has been accepted as being aligned
uint32_t NARNs = 0;						// not accepted because contains too many indeterminate bases
uint32_t NARNoHit = 0;					// not accepted as aligned because unable to find any potential hits
uint32_t NARMMDelta = 0;					// not accepted as aligned because MMDelta criteria not met
uint32_t NARMultiAlign = 0;				// not accepted as aligned because was multiloci aligned
uint32_t NARTrim = 0;						// not accepted as aligned because aligned read excessively trimmed
uint32_t NARSpliceJctn = 0;				// not accepted as aligned because aligned read orphan splice junction
uint32_t NARmicroInDel = 0;				// not accepted as aligned because aligned read orphan microInDel
uint32_t NARPCRdup = 0;					// not accepted as aligned because aligned read was a PCR duplicate
uint32_t NARNonUnique = 0;				// not accepted as aligned because aligned read was a duplicate sequence
uint32_t NARChromFilt = 0;				// not accepted as aligned because aligned to a filtered chrom
uint32_t NARRegionFilt = 0;				// not accepted as aligned because not aligned to a priority region
uint32_t NARPEInsertMin = 0;				// not accepted as aligned because both PE ends align but less than minimum insert size
uint32_t NARPEInsertMax = 0;				// not accepted as aligned because both PE ends align but more than maximum insert size
uint32_t NARPENoHit = 0;					// not accepted as aligned because PE aligning and although able to align this read was unable to align partner end
uint32_t NARPEStrand = 0;					// not accepted as aligned because PE aligning and although able to align this read other read was aligned to inconsistent strand
uint32_t NARPEChrom = 0;					// not accepted as aligned because PE aligning and although able to align this read other read was aligned to different chromosome
uint32_t NARPEUnalign = 0;				// not accepted as aligned because PE aligning and unable to accept this alignment
uint32_t NARLociConstrained = 0;			// not accepted as aligned because alignment violated loci base constraints
if((Rslt=SortReadHits(eRSMHitMatch,false)) < eBSFSuccess)
	{
	Reset(false);
	return(Rslt);
	}
while((pReadHit = IterReads(pReadHit))!=nullptr)
	{
	switch(pReadHit->NAR) {
		case eNARUnaligned:				// read has yet to be aligned
			NARUnaligned += 1;
			break;

		case eNARAccepted:									// accepted hit, accumulate strand count
			NARAccepted += 1;
			NumUniques += 1;
			if(pReadHit->HitLoci.Hit.Seg[0].Strand == '+')
				NumPlusHits += 1;
			if(pReadHit->HitLoci.Hit.FlgInDel)
				NumIndels += 1;
			if(pReadHit->HitLoci.Hit.FlgChimeric)
				NumChimeric += 1;
			if(pReadHit->HitLoci.Hit.FlgSplice)
				NumSpliced += 1;
			if(pReadHit->HitLoci.FlagTR)
				NumTrimmed += 1;
			break;

		case eNARNs:					// not accepted because contains too many indeterminate bases
//			NARNs += 1;
//			NumNoMatches += 1;
			break;
		case eNARNoHit:					// not accepted as aligned because unable to find any potential hits
			NARNoHit += 1;
			NumNoMatches += 1;
			break;
		case eNARMMDelta:				// not accepted as aligned because MMDelta (Hamming) criteria not met
			NARMMDelta += 1;
			NumHamming += 1;
			break;
		case eNARMultiAlign:			// not accepted as aligned because was multiloci aligned
			NARMultiAlign += 1;
			NumMultiMatches+= 1;
			break;
		case eNARTrim:					// not accepted as aligned because aligned read excessively trimmed
			NARTrim += 1;
			break;
		case eNARSpliceJctn:			// not accepted as aligned because aligned read orphan splice junction
			NARSpliceJctn += 1;
			break;
		case eNARmicroInDel:			// not accepted as aligned because aligned read orphan microInDel
			NARmicroInDel += 1;
			break;
		case eNARPCRdup:				// not accepted as aligned because aligned read was a PCR duplicate
			NARPCRdup += 1;
			break;
		case eNARNonUnique:				// not accepted as aligned because aligned read was a duplicate sequence
			NARNonUnique += 1;
			break;
		case eNARChromFilt:				// not accepted as aligned because aligned to a filtered chrom
			NARChromFilt += 1;
			break;
		case eNARRegionFilt:			// not accepted as aligned because not aligned to a priority region
			NARRegionFilt += 1;
			break;
		case eNARPEInsertMin:			// not accepted as aligned because both PE ends align but less than minimum insert size
			NARPEInsertMin += 1;
			break;
		case eNARPEInsertMax:			// not accepted as aligned because both PE ends align but more than maximum insert size
			NARPEInsertMax += 1;
			break;
		case eNARPENoHit:				// not accepted as aligned because PE aligning and although able to align this read was unable to align partner end
			NARPENoHit += 1;
			break;
		case eNARPEStrand:				// not accepted as aligned because PE aligning and although able to align this read other read was aligned to inconsistent strand
			NARPEStrand += 1;
			break;
		case eNARPEChrom:				// not accepted as aligned because PE aligning and although able to align this read other read was aligned to different chromosome
			NARPEChrom += 1;
			break;
		case eNARPEUnalign:				// not accepted as aligned because PE aligning and unable to accept this alignment
			NARPEUnalign += 1;
			break;
		case eNARLociConstrained:		// not accepted as aligned because alignment violated loci base constraints
			NARLociConstrained += 1;
			break;
		default:							
			break;
		}
	}

if(m_MLMode >= eMLall)
	NumNoMatches = m_TotNonAligned;
NumNoMatches += m_NumSloughedNs;
NARNs = m_NumSloughedNs;
gDiagnostics.DiagOut(eDLInfo,gszProcName,"From %u source reads there are %u accepted alignments, %u on '+' strand, %u on '-' strand", m_OrigNumReadsLoaded,NumUniques,NumPlusHits,NumUniques-NumPlusHits);
if(bSimReads)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"There are %u (%u 2 edge, %u 1 edge) high confidence aligned simulated reads with %u misaligned",NumReads2EdgeAligned + NumReads1EdgeAligned,NumReads2EdgeAligned,NumReads1EdgeAligned,NumReadsMisaligned);
if(NumIndels)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Of the accepted aligned reads, %u contained microInDels", NumIndels);
if(NumSpliced)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Of the accepted aligned reads, %u contained splice junctions", NumSpliced);
if(NumChimeric)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Of the accepted aligned reads, %u were chimeric", NumChimeric);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"A further %u multiloci aligned reads could not accepted as hits because they were unresolvable",NumMultiMatches);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"A further %u aligned reads were not accepted as hits because of insufficient Hamming edit distance",NumHamming);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"A further %u '+' and %u '-' strand aligned reads not accepted because of flank trimming (%d were trimmed) requirements",m_ElimPlusTrimed,m_ElimMinusTrimed,NumTrimmed);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Unable to align %u source reads of which %d were not aligned as they contained excessive number of indeterminate 'N' bases",NumNoMatches,NARNs);


gDiagnostics.DiagOut(eDLInfo,gszProcName,"Read nonalignment reason summary:");
gDiagnostics.DiagOut(eDLInfo,gszProcName,"   %u (%s) %s",NARUnaligned,m_NARdesc[eNARUnaligned].pszNAR, m_NARdesc[eNARUnaligned].pszNARdescr);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"   %u (%s) %s",NARAccepted,m_NARdesc[eNARAccepted].pszNAR, m_NARdesc[eNARAccepted].pszNARdescr);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"   %u (%s) %s",NARNs,m_NARdesc[eNARNs].pszNAR, m_NARdesc[eNARNs].pszNARdescr);

if(m_MLMode >= eMLall)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"   %u (%s) %s",NumNoMatches,m_NARdesc[eNARNoHit].pszNAR, m_NARdesc[eNARNoHit].pszNARdescr);
else
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"   %u (%s) %s",NARNoHit,m_NARdesc[eNARNoHit].pszNAR, m_NARdesc[eNARNoHit].pszNARdescr);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"   %u (%s) %s",NARMMDelta,m_NARdesc[eNARMMDelta].pszNAR, m_NARdesc[eNARMMDelta].pszNARdescr);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"   %u (%s) %s",NARMultiAlign,m_NARdesc[eNARMultiAlign].pszNAR, m_NARdesc[eNARMultiAlign].pszNARdescr);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"   %u (%s) %s",NARTrim,m_NARdesc[eNARTrim].pszNAR, m_NARdesc[eNARTrim].pszNARdescr);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"   %u (%s) %s",NARSpliceJctn,m_NARdesc[eNARSpliceJctn].pszNAR, m_NARdesc[eNARSpliceJctn].pszNARdescr);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"   %u (%s) %s",NARmicroInDel,m_NARdesc[eNARmicroInDel].pszNAR, m_NARdesc[eNARmicroInDel].pszNARdescr);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"   %u (%s) %s",NARPCRdup,m_NARdesc[eNARPCRdup].pszNAR, m_NARdesc[eNARPCRdup].pszNARdescr);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"   %u (%s) %s",NARNonUnique,m_NARdesc[eNARNonUnique].pszNAR, m_NARdesc[eNARNonUnique].pszNARdescr);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"   %u (%s) %s",NARChromFilt,m_NARdesc[eNARChromFilt].pszNAR, m_NARdesc[eNARChromFilt].pszNARdescr);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"   %u (%s) %s",NARRegionFilt,m_NARdesc[eNARRegionFilt].pszNAR, m_NARdesc[eNARRegionFilt].pszNARdescr);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"   %u (%s) %s",NARPEInsertMin,m_NARdesc[eNARPEInsertMin].pszNAR, m_NARdesc[eNARPEInsertMin].pszNARdescr);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"   %u (%s) %s",NARPEInsertMax,m_NARdesc[eNARPEInsertMax].pszNAR, m_NARdesc[eNARPEInsertMax].pszNARdescr);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"   %u (%s) %s",NARPENoHit,m_NARdesc[eNARPENoHit].pszNAR, m_NARdesc[eNARPENoHit].pszNARdescr);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"   %u (%s) %s",NARPEStrand,m_NARdesc[eNARPEStrand].pszNAR, m_NARdesc[eNARPEStrand].pszNARdescr);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"   %u (%s) %s",NARPEChrom,m_NARdesc[eNARPEChrom].pszNAR, m_NARdesc[eNARPEChrom].pszNARdescr);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"   %u (%s) %s",NARPEUnalign,m_NARdesc[eNARPEUnalign].pszNAR, m_NARdesc[eNARPEUnalign].pszNARdescr);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"   %u (%s) %s",NARLociConstrained,m_NARdesc[eNARLociConstrained].pszNAR, m_NARdesc[eNARLociConstrained].pszNARdescr);

if(gProcessingID > 0)
	{
	uint32_t Cricks;
	Cricks = NumUniques-NumPlusHits;
	if(m_MLMode >= eMLall)
		gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"Alignments",ePTUint32,sizeof(m_OrigNumReadsLoaded),"NumLoaded",&m_OrigNumReadsLoaded);
	else
		gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"Alignments",ePTUint32,sizeof(m_NumReadsLoaded),"NumLoaded",&m_NumReadsLoaded);
	gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NumUniques),"AcceptedAligned",&NumUniques);
	gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NumPlusHits),"AcceptedSense",&NumPlusHits);
	gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"Alignments",ePTUint32,sizeof(Cricks),"AcceptedAntisense",&Cricks);
	gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NumChimeric),"AcceptedChimeric",&NumChimeric);
	gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NumIndels),"AcceptedIndels",&NumIndels);
	gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NumSpliced),"AcceptedSpliced",&NumSpliced);
	gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NumMultiMatches),"RjctMultiMatches",&NumMultiMatches);
	gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NumHamming),"RjctNumHamming",&NumHamming);
	gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"Alignments",ePTUint32,sizeof(m_ElimPlusTrimed),"RjctPlusTrim",&m_ElimPlusTrimed);
	gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"Alignments",ePTUint32,sizeof(m_ElimMinusTrimed),"RjctMinusTrim",&m_ElimMinusTrimed);
	gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NumTrimmed),"NumTrimmed",&NumTrimmed);
	gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NumNoMatches),"Unalignable",&NumNoMatches);
	gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"Alignments",ePTUint32,sizeof(m_NumSloughedNs),"RjctExcessNs",&m_NumSloughedNs);


	gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NARUnaligned),"NARUnaligned",&NARUnaligned);
	gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NARAccepted),"NARAccepted",&NARAccepted);
	gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NARNs),"NARNs",&NARNs);

	if(m_MLMode >= eMLall)
		gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NumNoMatches),"NARNoHit",&NumNoMatches);
	else
		gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NARNoHit),"NARNoHit",&NARNoHit);
	
	gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NARMMDelta),"NARMMDelta",&NARMMDelta);
	gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NARMultiAlign),"NARMultiAlign",&NARMultiAlign);
	gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NARTrim),"NARTrim",&NARTrim);
	gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NARSpliceJctn),"NARSpliceJctn",&NARSpliceJctn);
	gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NARmicroInDel),"NARmicroInDel",&NARmicroInDel);
	gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NARPCRdup),"NARPCRdup",&NARPCRdup);
	gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NARNonUnique),"NARNonUnique",&NARNonUnique);
	gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NARChromFilt),"NARChromFilt",&NARChromFilt);
	gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NARRegionFilt),"NARRegionFilt",&NARRegionFilt);
	gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NARPEInsertMin),"NARPEInsertMin",&NARPEInsertMin);
	gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NARPEInsertMax),"NARPEInsertMax",&NARPEInsertMax);
	gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NARPENoHit),"NARPENoHit",&NARPENoHit);
	gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NARPEChrom),"NARPEChrom",&NARPEChrom);
	gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NARPEUnalign),"NARPEUnalign",&NARPEUnalign);
	gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"Alignments",ePTUint32,sizeof(NARLociConstrained),"NARLociConstrained",&NARLociConstrained);
	}

m_NARAccepted = NARAccepted;		// further downstream processing may be interested in number of reads being reported as accepted
return(eBSFSuccess);
}


int
CKAligner::ReportNoneAligned(void)
{
int ReadLen;
int SeqOfs;
int	NxtFastaCol;
int NumCols;
tsReadHit *pReadHit;
uint32_t SeqIdx;
etSeqBase Sequence[cMaxFastQSeqLen+ 10000];	// to hold sequence (sans quality scores) for current read
uint8_t *pSeqVal;
etSeqBase *pSeq;

int NumNoMatches = 0;
int NumMultiMatches = 0;
int NumHamming = 0;

int LineLen = 0;
// user interested in the non-alignable reads?
// these only include those which had no alignment at all
if(m_hNoneAlignFile != -1 || m_gzNoneAlignFile != nullptr)
	{
	SortReadHits(eRSMHitMatch,false);
	pReadHit = nullptr;
	while((pReadHit = IterSortedReads(pReadHit))!=nullptr)
		{
		if(!(pReadHit->NAR == eNARNs || pReadHit->NAR == eNARNoHit))
			continue;
		ReadLen = pReadHit->ReadLen;
		pSeqVal = &pReadHit->Read[pReadHit->DescrLen+1];
		pSeq = Sequence;
		for(SeqIdx = 0; SeqIdx < pReadHit->ReadLen; SeqIdx++,pSeq++,pSeqVal++)
			*pSeq = (*pSeqVal & 0x07);

		LineLen += sprintf(&m_pszLineBuff[LineLen],">lcl|na|%d %s %d|%d|%d\n",pReadHit->ReadID,pReadHit->Read,pReadHit->ReadID,pReadHit->NumReads,pReadHit->ReadLen);
		if((LineLen + (2 * cMaxFastQSeqLen)) > cAllocLineBuffSize)
			{
			if(m_hNoneAlignFile != -1)
				CUtility::RetryWrites(m_hNoneAlignFile,m_pszLineBuff,LineLen);
			else
				CUtility::SafeWrite_gz(m_gzNoneAlignFile,m_pszLineBuff,LineLen);
			LineLen = 0;
			}
		SeqOfs = 0;
		NxtFastaCol = 0;
		while(ReadLen)
			{
			NumCols = ReadLen > 70 ? 70 : ReadLen;
			if((NumCols + NxtFastaCol) > 70)
				NumCols = 70 - NxtFastaCol;
			CSeqTrans::MapSeq2Ascii(&Sequence[SeqOfs],NumCols,&m_pszLineBuff[LineLen]);
			LineLen += NumCols;
			NxtFastaCol += NumCols;
			SeqOfs += NumCols;
			ReadLen -= NumCols;
			if(!ReadLen || NxtFastaCol >= 70)
				{
				LineLen += sprintf(&m_pszLineBuff[LineLen],"\n");
				NxtFastaCol = 0;
				}

			if((LineLen + (2 * cMaxFastQSeqLen)) > cAllocLineBuffSize)
				{
				if(m_hNoneAlignFile != -1)
					CUtility::RetryWrites(m_hNoneAlignFile,m_pszLineBuff,LineLen);
				else
					CUtility::SafeWrite_gz(m_gzNoneAlignFile,m_pszLineBuff,LineLen);
				LineLen = 0;
				}
			}
		}
	if(LineLen)
		{
		if(m_hNoneAlignFile != -1)
			CUtility::RetryWrites(m_hNoneAlignFile,m_pszLineBuff,LineLen);
		else
			CUtility::SafeWrite_gz(m_gzNoneAlignFile,m_pszLineBuff,LineLen);
		}

	if(m_hNoneAlignFile != -1)
		{
#ifdef _WIN32
		_commit(m_hNoneAlignFile);
#else
		fsync(m_hNoneAlignFile);
#endif
		close(m_hNoneAlignFile);
		m_hNoneAlignFile = -1;
		}
	else
		{
		gzclose(m_gzNoneAlignFile);
		m_gzNoneAlignFile = nullptr;
		}
	}
return(eBSFSuccess);
}

int
CKAligner::ReportMultiAlign(void)
{
int ReadLen;
int SeqOfs;
int	NxtFastaCol;
int NumCols;
tsReadHit *pReadHit;
uint32_t SeqIdx;
etSeqBase Sequence[cMaxFastQSeqLen+ 10000];	// to hold sequence (sans quality scores) for current read
uint8_t *pSeqVal;
etSeqBase *pSeq;
int LineLen = 0;

if(m_hMultiAlignFile != -1 || m_gzMultiAlignFile != nullptr)
	{
	SortReadHits(eRSMHitMatch,false);
	pReadHit = nullptr;
	while((pReadHit = IterSortedReads(pReadHit))!=nullptr)
		{
		if(pReadHit->NAR != eNARMultiAlign)
			continue;

		ReadLen = pReadHit->ReadLen;
		pSeqVal = &pReadHit->Read[pReadHit->DescrLen+1];
		pSeq = Sequence;
		for(SeqIdx = 0; SeqIdx < pReadHit->ReadLen; SeqIdx++,pSeq++,pSeqVal++)
			*pSeq = (*pSeqVal & 0x07);

		LineLen += sprintf(&m_pszLineBuff[LineLen],">lcl|ml|%d %s %d|%d|%d\n",pReadHit->ReadID,pReadHit->Read,pReadHit->ReadID,pReadHit->NumReads,pReadHit->ReadLen);
		if((LineLen + (2 * cMaxFastQSeqLen)) > cAllocLineBuffSize)
			{
			if(m_hMultiAlignFile != -1)
				CUtility::RetryWrites(m_hMultiAlignFile,m_pszLineBuff,LineLen);
			else
				CUtility::SafeWrite_gz(m_gzMultiAlignFile,m_pszLineBuff,LineLen);
			LineLen = 0;
			}
		SeqOfs = 0;
		NxtFastaCol = 0;
		while(ReadLen)
			{
			NumCols = ReadLen > 70 ? 70 : ReadLen;
			if((NumCols + NxtFastaCol) > 70)
				NumCols = 70 - NxtFastaCol;
			CSeqTrans::MapSeq2Ascii(&Sequence[SeqOfs],NumCols,&m_pszLineBuff[LineLen]);
			LineLen += NumCols;
			NxtFastaCol += NumCols;
			SeqOfs += NumCols;
			ReadLen -= NumCols;
			if(!ReadLen || NxtFastaCol >= 70)
				{
				LineLen += sprintf(&m_pszLineBuff[LineLen],"\n");
				NxtFastaCol = 0;
				}

			if((LineLen + (2 * cMaxFastQSeqLen)) > cAllocLineBuffSize)
				{
				if(m_hMultiAlignFile != -1)
					CUtility::RetryWrites(m_hMultiAlignFile,m_pszLineBuff,LineLen);
				else
					CUtility::SafeWrite_gz(m_gzMultiAlignFile,m_pszLineBuff,LineLen);
				LineLen = 0;
				}
			}
		}
	if(LineLen)
		{
		if(m_hMultiAlignFile != -1)
			CUtility::RetryWrites(m_hMultiAlignFile,m_pszLineBuff,LineLen);
		else
			CUtility::SafeWrite_gz(m_gzMultiAlignFile,m_pszLineBuff,LineLen);
		}

	if(m_hMultiAlignFile != -1)
		{
#ifdef _WIN32
		_commit(m_hMultiAlignFile);
#else
		fsync(m_hMultiAlignFile);
#endif
		close(m_hMultiAlignFile);
		m_hMultiAlignFile = -1;
		}
	else
		{
		gzclose(m_gzMultiAlignFile);
		m_gzMultiAlignFile = nullptr;
		}
	}
return(eBSFSuccess);
}



int
CKAligner::FiltByChroms(void)
{
tBSFEntryID PrevTargEntry;
tBSFEntryID ExcludeTargEntryID;
char szChromName[128];
bool bProcChrom = false;
int MatchesFiltOut = 0;
tsReadHit *pReadHit;

SortReadHits(eRSMHitMatch,false);

if(m_RegExprs.HasRegExprs())
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Now filtering matches by chromosome");
	PrevTargEntry = 0;
	pReadHit = nullptr;
	while((pReadHit = IterSortedReads(pReadHit))!=nullptr)
		{
		if(pReadHit->NAR == eNARAccepted)
			{
			if(pReadHit->HitLoci.Hit.Seg[0].ChromID != (uint32_t)PrevTargEntry)
				{
				m_pSfxArray->GetIdentName(pReadHit->HitLoci.Hit.Seg[0].ChromID,sizeof(szChromName),szChromName);
				PrevTargEntry = pReadHit->HitLoci.Hit.Seg[0].ChromID;
				ExcludeTargEntryID = (tBSFEntryID)-1;
				}
			else
				{
				if(pReadHit->HitLoci.Hit.Seg[0].ChromID == (uint32_t)ExcludeTargEntryID)
					{
					if(pReadHit->NAR == eNARAccepted)
						MatchesFiltOut += 1;
					pReadHit->NAR = eNARChromFilt;
					pReadHit->NumHits = 0;
					pReadHit->LowHitInstances = 0;
					}
				continue;
				}

			// to be included?
			if((bProcChrom = !m_RegExprs.MatchExcludeRegExpr(szChromName)) == true)
				bProcChrom = m_RegExprs.MatchIncludeRegExpr(szChromName);

			if(!bProcChrom)
				{
				ExcludeTargEntryID = pReadHit->HitLoci.Hit.Seg[0].ChromID;
				if(pReadHit->NAR  == eNARAccepted)
					MatchesFiltOut += 1;
				pReadHit->NAR = eNARChromFilt;
				pReadHit->NumHits = 0;
				pReadHit->LowHitInstances = 0;
				}
			else
				ExcludeTargEntryID = -1;
			}
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Filtering by chromosome completed - removed %d  matches",MatchesFiltOut);

	if(MatchesFiltOut)
		SortReadHits(eRSMHitMatch,false,true);
	if(gProcessingID > 0)
		gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"Filtered",ePTInt32,sizeof(MatchesFiltOut),"Chroms",&MatchesFiltOut);
		
	}
return(eBSFSuccess);
}

// remove hits not aligned into a priority region
int
CKAligner::FiltByPriorityRegions(void) // remove hits not aligned into a priority region
{
char szPriorityChromName[100];
int PriorityChromID;

uint32_t MatchesFiltOut = 0;
uint32_t MatchesFiltIn = 0;

tsReadHit *pReadHit;
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Now filtering matches by prioritorised regions");
SortReadHits(eRSMHitMatch,false);
pReadHit = nullptr;
while((pReadHit = IterSortedReads(pReadHit))!=nullptr)
	{
	if(pReadHit->NAR == eNARAccepted)
		{
		// check if hit loci within region designated as being a priority exact matching region
		if(m_pSfxArray->GetIdentName(pReadHit->HitLoci.Hit.Seg[0].ChromID,sizeof(szPriorityChromName),szPriorityChromName)!=eBSFSuccess)
			{
			MatchesFiltOut += 1;
			pReadHit->NAR = eNARRegionFilt;
			pReadHit->NumHits = 0;
			pReadHit->HitLoci.Hit.Seg[0].ChromID = 0;
			pReadHit->HitLoci.Hit.Seg[1].ChromID = 0;
			continue;
			}
		if((PriorityChromID = m_pPriorityRegionBED->LocateChromIDbyName(szPriorityChromName)) < 1)
			{
			MatchesFiltOut += 1;
			pReadHit->NAR = eNARRegionFilt;
			pReadHit->NumHits = 0;
			pReadHit->HitLoci.Hit.Seg[0].ChromID = 0;
			pReadHit->HitLoci.Hit.Seg[1].ChromID = 0;
			continue;
			}

		if(!m_pPriorityRegionBED->InAnyFeature(PriorityChromID,(int)pReadHit->HitLoci.Hit.Seg[0].MatchLoci,
										(int)(pReadHit->HitLoci.Hit.Seg[0].MatchLoci + pReadHit->HitLoci.Hit.Seg[0].MatchLen-1)))
			{
			MatchesFiltOut += 1;
			pReadHit->NAR = eNARRegionFilt;
			pReadHit->NumHits = 0;
			pReadHit->HitLoci.Hit.Seg[0].ChromID = 0;
			pReadHit->HitLoci.Hit.Seg[1].ChromID = 0;
			}
		else
			MatchesFiltIn += 1;
		}
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Filtering by prioritised regions completed - retained %u, removed %u matches",MatchesFiltIn,MatchesFiltOut);

if(MatchesFiltOut)
	SortReadHits(eRSMHitMatch,false,true);

if(gProcessingID > 0)
	{
	gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"Filtered",ePTInt32,sizeof(MatchesFiltIn),"RegionRetained",&MatchesFiltIn);
	gSQLiteSummaries.AddResult(gExperimentID, gProcessingID,(char *)"Filtered",ePTInt32,sizeof(MatchesFiltOut),"RegionRemoved",&MatchesFiltOut);
	}
return(eBSFSuccess);
}


int
CKAligner::WriteBasicCountStats(void)
{
int Idx;
int SeqOfs;
char szLineBuff[(cMaxReadLen + 1000) * 8];
int BuffIdx;
int PhredBand;

if(m_hStatsFile > 0 && m_MaxAlignLen > 0)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Writing out basic count stats to file");

	// report multiple hit distribution if MLMode not simply the default
	BuffIdx = 0;
	if(m_MLMode > eMLdefault)
		{
		BuffIdx = sprintf(szLineBuff,"\"Multihit distribution\"\n,");
		for(Idx=0; Idx < m_MaxMLmatches; Idx++)
			{
			BuffIdx += sprintf(&szLineBuff[BuffIdx],",%d",Idx+1);
			if((BuffIdx + 512) > sizeof(szLineBuff))
				{
				CUtility::RetryWrites(m_hStatsFile,szLineBuff,BuffIdx);
				BuffIdx = 0;
				}
			}
		BuffIdx += sprintf(&szLineBuff[BuffIdx],"\n,,\"Instances\"");
		for(Idx=0; Idx < m_MaxMLmatches; Idx++)
			{
			BuffIdx += sprintf(&szLineBuff[BuffIdx],",%d",m_MultiHitDist[Idx]);
			if((BuffIdx + 512) > sizeof(szLineBuff))
				{
				CUtility::RetryWrites(m_hStatsFile,szLineBuff,BuffIdx);
				BuffIdx = 0;
				}
			}
		BuffIdx += sprintf(&szLineBuff[BuffIdx],"\n");
		}

	BuffIdx += sprintf(&szLineBuff[BuffIdx],"\"Phred Score Instances\"\n,\"Psn\"");
	for(SeqOfs = 0; SeqOfs < m_MaxAlignLen; SeqOfs++)
		{
		BuffIdx += sprintf(&szLineBuff[BuffIdx],",%d",SeqOfs+1);
		if((BuffIdx + 512) > sizeof(szLineBuff))
			{
			CUtility::RetryWrites(m_hStatsFile,szLineBuff,BuffIdx);
			BuffIdx = 0;
			}
		}

	for(PhredBand = 0; PhredBand <= 3; PhredBand++)
		{
		switch(PhredBand) {
			case 0:
				BuffIdx += sprintf(&szLineBuff[BuffIdx],"\n,\"Phred 0..9\"");
				break;
			case 1:
				BuffIdx += sprintf(&szLineBuff[BuffIdx],"\n,\"Phred 10..19\"");
				break;
			case 2:
				BuffIdx += sprintf(&szLineBuff[BuffIdx],"\n,\"Phred 20..29\"");
				break;
			case 3:
				BuffIdx += sprintf(&szLineBuff[BuffIdx],"\n,\"Phred 30+\"");
				break;
			}

		for(SeqOfs = 0; SeqOfs < m_MaxAlignLen; SeqOfs++)
			{
			BuffIdx += sprintf(&szLineBuff[BuffIdx],",%d",m_AlignQSubDist[PhredBand][SeqOfs].QInsts);
			if((BuffIdx + 512) > sizeof(szLineBuff))
				{
				CUtility::RetryWrites(m_hStatsFile,szLineBuff,BuffIdx);
				BuffIdx = 0;
				}
			}
		}

	BuffIdx += sprintf(&szLineBuff[BuffIdx],"\n\"Aligner Induced Subs\"\n,\"Psn\"");
	for(SeqOfs = 0; SeqOfs < m_MaxAlignLen; SeqOfs++)
		{
		BuffIdx += sprintf(&szLineBuff[BuffIdx],",%d",SeqOfs+1);
		if((BuffIdx + 512) > sizeof(szLineBuff))
			{
			CUtility::RetryWrites(m_hStatsFile,szLineBuff,BuffIdx);
			BuffIdx = 0;
			}
		}

	for(PhredBand = 0; PhredBand <= 3; PhredBand++)
		{
		switch(PhredBand) {
			case 0:
				BuffIdx += sprintf(&szLineBuff[BuffIdx],"\n,\"Phred 0..8\"");
				break;
			case 1:
				BuffIdx += sprintf(&szLineBuff[BuffIdx],"\n,\"Phred 9..19\"");
				break;
			case 2:
				BuffIdx += sprintf(&szLineBuff[BuffIdx],"\n,\"Phred 20..29\"");
				break;
			case 3:
				BuffIdx += sprintf(&szLineBuff[BuffIdx],"\n,\"Phred 30+\"");
				break;
			}

		for(SeqOfs = 0; SeqOfs < m_MaxAlignLen; SeqOfs++)
			{
			BuffIdx += sprintf(&szLineBuff[BuffIdx],",%d",m_AlignQSubDist[PhredBand][SeqOfs].Subs);
			if((BuffIdx + 512) > sizeof(szLineBuff))
				{
				CUtility::RetryWrites(m_hStatsFile,szLineBuff,BuffIdx);
				BuffIdx = 0;
				}
			}
		}

	BuffIdx += sprintf(&szLineBuff[BuffIdx],"\n\"Multiple substitutions\"\n,\"NumSubs\"");
	for(SeqOfs = 0; SeqOfs < m_MaxAlignLen; SeqOfs++)
		{
		BuffIdx += sprintf(&szLineBuff[BuffIdx],",%d",SeqOfs);
		if((BuffIdx + 512) > sizeof(szLineBuff))
			{
			CUtility::RetryWrites(m_hStatsFile,szLineBuff,BuffIdx);
			BuffIdx = 0;
			}
		}
	BuffIdx += sprintf(&szLineBuff[BuffIdx],"\n,\"Instances\"");
	for(SeqOfs = 0; SeqOfs < m_MaxAlignLen; SeqOfs++)
		{
		BuffIdx += sprintf(&szLineBuff[BuffIdx],",%d",m_AlignMSubDist[SeqOfs]);
		if((BuffIdx + 512) > sizeof(szLineBuff))
			{
			CUtility::RetryWrites(m_hStatsFile,szLineBuff,BuffIdx);
			BuffIdx = 0;
			}
		}
	BuffIdx += sprintf(&szLineBuff[BuffIdx],"\n");
	CUtility::RetryWrites(m_hStatsFile,szLineBuff,BuffIdx);
	BuffIdx = 0;
	}
return(eBSFSuccess);
}


bool								// true if file to be generated compressed with gzopen/gzwrite/gzclose
CKAligner::FileReqWriteCompr(char *pszFile) // If last 3 chars of file name is ".gz" then this file is assumed to require compression
{
int Len;
if(pszFile == nullptr || pszFile[0] == '\0')
	return(false);
if((Len = (int)strlen(pszFile)) < 4)
	return(false);
return(stricmp(".gz",&pszFile[Len-3]) == 0 ? true : false);
}

// Create, or if file already exists then truncate, the multitude of results files which user may have requested
// If a primary output (m_pszOutFile, m_pszNoneAlignFile or m_pszMultiAlignFile) file name ends with '.gz' then will create file with gzopen ready for output compression
// If output is to SAM or BAM then CSAMfile will handle the file processing
int
CKAligner::CreateOrTruncResultFiles(void)
{
int LineLen;

if(m_bPackedBaseAlleles)
	{
#ifdef _WIN32
	m_hPackedBaseAllelesFile = open(m_pszOutFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE) );
#else
	if((m_hPackedBaseAllelesFile = open(m_pszOutFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE))!=-1)
		if(ftruncate(m_hPackedBaseAllelesFile,0)!=0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate loci base classification binary output file  %s - %s",m_pszOutFile,strerror(errno));
			return(eBSFerrCreateFile);
			}
#endif
	if(m_hPackedBaseAllelesFile < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate loci base classification binary output file '%s'",m_pszOutFile);
		return(eBSFerrCreateFile);
		}

	CUtility::AppendFileNameSuffix(m_szCoverageSegmentsFile, m_pszOutFile, (char*)".covsegs.wig",'.');

#ifdef _WIN32
	m_hWIGSpansFile = open(m_szCoverageSegmentsFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE) );
#else
	if((m_hWIGSpansFile = open(m_szCoverageSegmentsFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE))!=-1)
		if(ftruncate(m_hWIGSpansFile,0)!=0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate coverage segments file  %s - %s",m_szCoverageSegmentsFile,strerror(errno));
			return(eBSFerrCreateFile);
			}
#endif
	if(m_hWIGSpansFile < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate coverage segments output file '%s'",m_szCoverageSegmentsFile);
		return(eBSFerrCreateFile);
		}
	return(eBSFSuccess);
	}

// determine which primary output files need to be generated as compressed
m_bgzOutFile = FileReqWriteCompr(m_pszOutFile);
m_bgzNoneAlignFile = FileReqWriteCompr(m_pszNoneAlignFile);
m_bgzMultiAlignFile = FileReqWriteCompr(m_pszMultiAlignFile);


// create/truncate all output reporting files
if(m_FMode < eFMsam)
	{
	if(!m_bgzOutFile)
		{
#ifdef _WIN32
		m_hOutFile = open(m_pszOutFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE));
#else
		if((m_hOutFile = open(m_pszOutFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE))!=-1)
			if(ftruncate(m_hOutFile,0)!=0)
					{
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",m_pszOutFile,strerror(errno));
					return(eBSFerrCreateFile);
					}
#endif
		if(m_hOutFile < 0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output file '%s'",m_pszOutFile);
			return(eBSFerrCreateFile);
			}
		}
	else
		{
		m_gzOutFile = gzopen(m_pszOutFile,"wb");
		if(m_gzOutFile == nullptr)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output file '%s'",m_pszOutFile);
			return(eBSFerrCreateFile);
			}
		gzbuffer(m_gzOutFile,cAllocLineBuffSize);		// large buffer to reduce number of writes required
		}
	}
else
	{
	m_hOutFile = -1;
	m_gzOutFile = nullptr;
	}

if((m_pszLineBuff = new char [cAllocLineBuffSize])==nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to allocate memory (%d bytes) for output line buffer",cAllocLineBuffSize);
	return(eBSFerrMem);
	}

if(m_FMode == eFMbed && m_MLMode == eMLall)
	{
	LineLen = sprintf(m_pszLineBuff,"track type=bed name=\"%s\" description=\"%s\"\n",m_pszTrackTitle,m_pszTrackTitle);
	if(!m_bgzOutFile)
		CUtility::RetryWrites(m_hOutFile,m_pszLineBuff,LineLen);
	else
		CUtility::SafeWrite_gz(m_gzOutFile,m_pszLineBuff,LineLen);
	LineLen = 0;
	}
LineLen = 0;
m_szLineBuffIdx = 0;

// If BED output requested then separate BEDs for InDels and splice junctions
if(m_microInDelLen && m_FMode == eFMbed)
	{
	strcpy(m_szIndRsltsFile,m_pszOutFile);
	if(m_bgzOutFile)								// if compressing the primary alignment results file then remove the ".gz' file suffix 
		m_szJctRsltsFile[strlen(m_pszOutFile)-3] = '\0';
	strcat(m_szIndRsltsFile,".ind");
#ifdef _WIN32
	m_hIndOutFile = open(m_szIndRsltsFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE) );
#else
	if((m_hIndOutFile = open(m_szIndRsltsFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE))!=-1)
		if(ftruncate(m_hIndOutFile,0)!=0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate microInDel file %s - %s",m_szIndRsltsFile,strerror(errno));
			return(eBSFerrCreateFile);
			}
#endif
	if(m_hIndOutFile < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate microInDel output file '%s'",m_szIndRsltsFile);
		return(eBSFerrCreateFile);
		}
	}

if(m_SpliceJunctLen && (m_FMode == eFMbed || m_FMode == eFMsam || m_FMode == eFMsamAll))
	{
	strcpy(m_szJctRsltsFile,m_pszOutFile);
	if(m_bgzOutFile)								// if compressing the primary alignment results file then remove the ".gz' file suffix 
		m_szJctRsltsFile[strlen(m_pszOutFile)-3] = '\0';
	strcat(m_szJctRsltsFile,".jct");
#ifdef _WIN32
	m_hJctOutFile = open(m_szJctRsltsFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE) );
#else
	if((m_hJctOutFile = open(m_szJctRsltsFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE))!=-1)
		if(ftruncate(m_hJctOutFile,0)!=0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate  splice junct  %s - %s",m_szJctRsltsFile,strerror(errno));
			return(eBSFerrCreateFile);
			}
#endif
	if(m_hJctOutFile < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate splice junct output file '%s'",m_szJctRsltsFile);
		return(eBSFerrCreateFile);
		}
	}


if(m_pszSitePrefsFile != nullptr && m_pszSitePrefsFile[0] != '\0')
	{
#ifdef _WIN32
	m_hSitePrefsFile = open(m_pszSitePrefsFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE) );
#else
	if((m_hSitePrefsFile = open(m_pszSitePrefsFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE))!=-1)
	   if(ftruncate(m_hSitePrefsFile,0)!=0)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",m_pszSitePrefsFile,strerror(errno));
				return(eBSFerrCreateFile);
				}
#endif

	if(m_hSitePrefsFile < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate site preferencing file '%s'",m_pszSitePrefsFile);
		return(eBSFerrCreateFile);
		}
	}
else
	m_hSitePrefsFile = -1;


if(m_pszSNPRsltsFile != nullptr && m_pszSNPRsltsFile[0] != '\0' && m_MinSNPreads > 0 && m_FMode <= eFMsam)
	{
#ifdef _WIN32
	m_hSNPfile = open(m_pszSNPRsltsFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE) );
#else
	if((m_hSNPfile = open(m_pszSNPRsltsFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE))!=-1)
		if(ftruncate(m_hSNPfile,0)!=0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate SNP file  %s - %s",m_pszSNPRsltsFile,strerror(errno));
			return(eBSFerrCreateFile);
			}
#endif
	if(m_hSNPfile < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate SNP output file '%s'",m_pszSNPRsltsFile);
		return(eBSFerrCreateFile);
		}

	CUtility::AppendFileNameSuffix(m_szCoverageSegmentsFile, m_pszSNPRsltsFile, (char*)".covsegs.wig",'.');

#ifdef _WIN32
	m_hWIGSpansFile = open(m_szCoverageSegmentsFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE) );
#else
	if((m_hWIGSpansFile = open(m_szCoverageSegmentsFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE))!=-1)
		if(ftruncate(m_hWIGSpansFile,0)!=0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate coverage segments file  %s - %s",m_szCoverageSegmentsFile,strerror(errno));
			return(eBSFerrCreateFile);
			}
#endif
	if(m_hWIGSpansFile < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate coverage segments output file '%s'",m_szCoverageSegmentsFile);
		return(eBSFerrCreateFile);
		}


	if(m_pszSNPCentroidFile != nullptr && m_pszSNPCentroidFile[0] != '\0')
		{
#ifdef _WIN32
		m_hSNPCentsfile = open(m_pszSNPCentroidFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE) );
#else
		if((m_hSNPCentsfile = open(m_pszSNPCentroidFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE))!=-1)
			if(ftruncate(m_hSNPCentsfile,0)!=0)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate SNP file  %s - %s",m_pszSNPCentroidFile,strerror(errno));
				return(eBSFerrCreateFile);
				}
#endif
		if(m_hSNPCentsfile < 0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate SNP output file '%s'",m_pszSNPCentroidFile);
			return(eBSFerrCreateFile);
			}
		}
	else
		m_hSNPCentsfile = -1;

#ifdef _DISNPS_
	CUtility::AppendFileNameSuffix(m_szDiSNPFile, m_pszSNPRsltsFile, (char*)".disnp.csv",'.');
	CUtility::AppendFileNameSuffix(m_szTriSNPFile, m_pszSNPRsltsFile, (char*)".trisnp.csv", '.');

	if(m_szDiSNPFile[0] != '\0')
		{
#ifdef _WIN32
		m_hDiSNPfile = open(m_szDiSNPFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE) );
#else
		if((m_hDiSNPfile = open(m_szDiSNPFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE))!=-1)
			if(ftruncate(m_hDiSNPfile,0)!=0)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate DiSNP file  %s - %s",m_szDiSNPFile,strerror(errno));
				return(eBSFerrCreateFile);
				}
#endif
		if(m_hDiSNPfile < 0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate DiSNP output file '%s'",m_szDiSNPFile);
			return(eBSFerrCreateFile);
			}
		}
	else
		m_hDiSNPfile = -1;

	if(m_szTriSNPFile[0] != '\0')
		{
#ifdef _WIN32
		m_hTriSNPfile = open(m_szTriSNPFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE) );
#else
		if((m_hTriSNPfile = open(m_szTriSNPFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE))!=-1)
			if(ftruncate(m_hTriSNPfile,0)!=0)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate TriSNP file  %s - %s",m_szTriSNPFile,strerror(errno));
				return(eBSFerrCreateFile);
				}
#endif
		if(m_hTriSNPfile < 0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate TriSNP output file '%s'",m_szTriSNPFile);
			return(eBSFerrCreateFile);
			}
		}
	else
		m_hTriSNPfile = -1;
#endif

	if(m_pszMarkerFile != nullptr && m_pszMarkerFile[0] != '\0')
		{
#ifdef _WIN32
		m_hMarkerFile = open(m_pszMarkerFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE) );
#else
		if((m_hMarkerFile = open(m_pszMarkerFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE))!=-1)
			if(ftruncate(m_hMarkerFile,0)!=0)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate SNP file  %s - %s",m_pszMarkerFile,strerror(errno));
				return(eBSFerrCreateFile);
				}
#endif
		if(m_hMarkerFile < 0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate SNP output file '%s'",m_pszMarkerFile);
			return(eBSFerrCreateFile);
			}
		}
	else
		m_hMarkerFile = -1;
	}
else
	{
	m_hSNPfile = -1;
	m_hSNPCentsfile = -1;
	m_hMarkerFile = -1;
	m_hDiSNPfile = -1;
	m_hTriSNPfile = -1;
	};

// none-aligned fasta reads file to be also generated?
if(m_pszNoneAlignFile != nullptr && m_pszNoneAlignFile[0] != '\0')
	{
	if(!m_bgzNoneAlignFile)
		{
#ifdef _WIN32
		m_hNoneAlignFile = open(m_pszNoneAlignFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE) );
#else
		if((m_hNoneAlignFile = open(m_pszNoneAlignFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE))!=-1)
			if(ftruncate(m_hNoneAlignFile,0)!=0)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate nonaligned %s - %s",m_pszNoneAlignFile,strerror(errno));
				return(eBSFerrCreateFile);
				}
#endif

		if(m_hNoneAlignFile < 0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output nonaligned file '%s'",m_pszNoneAlignFile);
			return(eBSFerrCreateFile);
			}
		}
	else
		{
		m_gzNoneAlignFile = gzopen(m_pszNoneAlignFile,"wb");
		if(m_gzNoneAlignFile == nullptr)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate nonaligned file '%s'",m_pszNoneAlignFile);
			return(eBSFerrCreateFile);
			}
		gzbuffer(m_gzNoneAlignFile,cAllocLineBuffSize);		// large buffer to reduce number of writes required
		}
	}
else
	{
	m_hNoneAlignFile = -1;
	m_gzNoneAlignFile = nullptr;
	}

// multialigned fasta reads file to be also generated?
if(m_pszMultiAlignFile != nullptr && m_pszMultiAlignFile[0] != '\0')
	{
	if(!m_bgzMultiAlignFile)
		{
#ifdef _WIN32
		m_hMultiAlignFile = open(m_pszMultiAlignFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE) );
#else
		if((m_hMultiAlignFile = open(m_pszMultiAlignFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE))!=-1)
			if(ftruncate(m_hMultiAlignFile,0)!=0)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate multiple alignment file %s - %s",m_pszMultiAlignFile,strerror(errno));
				return(eBSFerrCreateFile);
				}
#endif

		if(m_hMultiAlignFile < 0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate multiple alignment file '%s'",m_pszMultiAlignFile);
			return(eBSFerrCreateFile);
			}
		}
	else
		{
		m_gzMultiAlignFile = gzopen(m_pszMultiAlignFile,"wb");
		if(m_gzMultiAlignFile == nullptr)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate multialigned file '%s'",m_pszMultiAlignFile);
			return(eBSFerrCreateFile);
			}
		gzbuffer(m_gzMultiAlignFile,cAllocLineBuffSize);		// large buffer to reduce number of writes required
		}
	}
else
	{
	m_hMultiAlignFile = -1;
	m_bgzMultiAlignFile = false;
	}

// substitution stats to be also generated? Note that these can't be generated for InDels or splice junctions
if(m_pszStatsFile != nullptr && m_pszStatsFile[0] != '\0')
	{
#ifdef _WIN32
	m_hStatsFile = open(m_pszStatsFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE) );
#else
	if((m_hStatsFile = open(m_pszStatsFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE))!=-1)
		if(ftruncate(m_hStatsFile,0)!=0)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate stats %s - %s",m_pszStatsFile,strerror(errno));
				return(eBSFerrCreateFile);
				}
#endif

	if(m_hStatsFile < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output stats file '%s'",m_pszStatsFile);
		return(eBSFerrCreateFile);
		}
	}
else
	m_hStatsFile = -1;


// PE insert lengths to be also generated? Note that these can't be generated for InDels or splice junctions
if(m_pszPEInsertDistFile != nullptr && m_pszPEInsertDistFile[0] != '\0')
	{
#ifdef _WIN32
	m_hInsertLensFile = open(m_pszPEInsertDistFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE) );
#else
	if((m_hInsertLensFile = open(m_pszPEInsertDistFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE))!=-1)
		if(ftruncate(m_hInsertLensFile,0)!=0)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate stats %s - %s",m_pszPEInsertDistFile,strerror(errno));
				return(eBSFerrCreateFile);
				}
#endif

	if(m_hInsertLensFile < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output stats file '%s'",m_pszPEInsertDistFile);
		return(eBSFerrCreateFile);
		}
	}
else
	m_hInsertLensFile = -1;

return(eBSFSuccess);
}

int
CKAligner::CompileChromRegExprs(int	NumIncludeChroms,	// number of chromosome regular expressions to include
		char **ppszIncludeChroms,		// array of include chromosome regular expressions
		int	NumExcludeChroms,			// number of chromosome expressions to exclude
		char **ppszExcludeChroms)		// array of exclude chromosome regular expressions
{
return(m_RegExprs.CompileREs(NumIncludeChroms, ppszIncludeChroms,NumExcludeChroms, ppszExcludeChroms));
}

#ifdef _WIN32
unsigned __stdcall KLoadReadFilesThread(void * pThreadPars)
#else
void *KLoadReadFilesThread(void * pThreadPars)
#endif
{
int Rslt;
tsLoadReadsThreadPars *pPars = (tsLoadReadsThreadPars *)pThreadPars;			// makes it easier not having to deal with casts!
CKAligner *pKAligner = (CKAligner *)pPars->pThis;
Rslt = pKAligner->ProcLoadReadFiles(pPars);
pPars->Rslt = Rslt;
#ifdef _WIN32
_endthreadex(0);
return(eBSFSuccess);
#else
pthread_exit(&pPars->Rslt);
#endif
}

int
CKAligner::InitiateLoadingReads(void)
{
tsLoadReadsThreadPars ThreadPars;
memset(&ThreadPars, 0, sizeof(ThreadPars));
// initiate loading the reads
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading reads from file...");
m_ThreadLoadReadsRslt = -1;


ThreadPars.pRslt = &m_ThreadLoadReadsRslt;
ThreadPars.pThis = this;
ThreadPars.SampleNthRawRead = m_SampleNthRawRead;
ThreadPars.Rslt = 0;


#ifdef _WIN32
m_hThreadLoadReads = ThreadPars.threadHandle = (HANDLE)_beginthreadex(nullptr,0x0fffff,KLoadReadFilesThread,&ThreadPars,0,&m_ThreadLoadReadsID);
if (m_hThreadLoadReads == nullptr)
	{
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Failed to initiate read loading thread ...");
	return(-1);
	}
#else
int ThreadRslt = ThreadPars.threadRslt = pthread_create (&m_ThreadLoadReadsID , nullptr , KLoadReadFilesThread , &ThreadPars);
if (ThreadRslt != 0)
	{
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Failed to initiate read loading thread ...");
	return(-1);
	}
#endif

// allow a few seconds for threads to actually startup
#ifdef _WIN32
Sleep(2000);
#else
sleep(2);
#endif

// wait a few additional seconds, if major problems with loading reads then these should show up very quickly
#ifdef _WIN32
if(WaitForSingleObject(m_hThreadLoadReads, 1000) != WAIT_TIMEOUT)
	{
	CloseHandle(m_hThreadLoadReads);
	m_hThreadLoadReads = nullptr;
	return(m_ThreadLoadReadsRslt);
	}
#else
struct timespec ts;
int JoinRlt;
clock_gettime(CLOCK_REALTIME, &ts);
ts.tv_sec += 1;
if((JoinRlt = pthread_timedjoin_np(m_ThreadLoadReadsID, nullptr, &ts)) == 0)
	{
	m_ThreadLoadReadsID = 0;
	return(m_ThreadLoadReadsRslt);
	}
#endif
return(eBSFSuccess);
}


#ifdef _WIN32
unsigned __stdcall KAssignMultiMatchesThread(void * pThreadPars)
#else
void *KAssignMultiMatchesThread(void * pThreadPars)
#endif
{
int Rslt;
tsClusterThreadPars *pPars = (tsClusterThreadPars *)pThreadPars; // makes it easier not having to deal with casts!
CKAligner *pThis = (CKAligner *)pPars->pThis;
Rslt = pThis->ProcAssignMultiMatches(pPars);
pPars->Rslt = Rslt;
#ifdef _WIN32
_endthreadex(0);
return(eBSFSuccess); // unreached, but keeps compilers happy!
#else
pthread_exit(&pPars->Rslt);
#endif
}

int
CKAligner::RunClusteringThreads(int NumThreads)
{
int ThreadIdx;
uint32_t ClusterStartIdx;
tsClusterThreadPars ClusterThreads[cMaxWorkerThreads];

ClusterStartIdx = 0;
memset(ClusterThreads,0,sizeof(ClusterThreads));
for(ThreadIdx = 0; ThreadIdx < NumThreads; ThreadIdx++)
	{
	ClusterThreads[ThreadIdx].ThreadIdx = ThreadIdx + 1;
	ClusterThreads[ThreadIdx].pThis = this;
#ifdef _WIN32
	ClusterThreads[ThreadIdx].threadHandle = (HANDLE)_beginthreadex(nullptr,0x0fffff,KAssignMultiMatchesThread,&ClusterThreads[ThreadIdx],0,&ClusterThreads[ThreadIdx].threadID);
#else
	ClusterThreads[ThreadIdx].threadRslt =	pthread_create (&ClusterThreads[ThreadIdx].threadID , nullptr , KAssignMultiMatchesThread , &ClusterThreads[ThreadIdx] );
#endif
	}

for(ThreadIdx = 0; ThreadIdx < NumThreads; ThreadIdx++)
	{
#ifdef _WIN32
	while(WAIT_TIMEOUT == WaitForSingleObject( ClusterThreads[ThreadIdx].threadHandle, 60000))
		{
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Progress: Still clustering ...");
		}
	CloseHandle( ClusterThreads[ThreadIdx].threadHandle);
#else
	struct timespec ts;
	int JoinRlt;
	clock_gettime(CLOCK_REALTIME, &ts);
	ts.tv_sec += 60;
	while((JoinRlt = pthread_timedjoin_np(ClusterThreads[ThreadIdx].threadID, nullptr, &ts)) != 0)
		{
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Progress: Still clustering ...");
		ts.tv_sec += 60;
		}
#endif
	}
return(eBSFSuccess);
}

// GetClusterStartEnd
// Clustering threads call this function to obtain the next from, until inclusive indexes into m_pMultiHits[] to be clustered
// Returns number of
int												// returns 0 if finished clustering or cnt of multihits to be processed by this thread
CKAligner::GetClusterStartEnd(uint32_t *pMatchFrom,			// cluster from this inclusive index
					uint32_t *pMatchUntil)		// until this inclusive index
{
int Rslt;
uint32_t NumLeft2Cluster;
uint32_t Num4Thread;
AcquireSerialise();
NumLeft2Cluster = m_NumMultiHits - m_CurClusterFrom;
if(NumLeft2Cluster > 0)
	{
	if(NumLeft2Cluster < 100)	// if < 100 yet to be processed then give it all to the one thread
		Num4Thread = NumLeft2Cluster;
	else
		{
		Num4Thread = min(2000u,(uint32_t)m_NumThreads + (NumLeft2Cluster / (uint32_t)m_NumThreads));
		Num4Thread = min(Num4Thread,NumLeft2Cluster);
		}

	*pMatchFrom = m_CurClusterFrom;
	*pMatchUntil = m_CurClusterFrom + Num4Thread - 1;
	m_CurClusterFrom = min(1 + *pMatchUntil,m_NumMultiHits);
	Rslt = Num4Thread;
	}
else
	Rslt = 0;
ReleaseSerialise();
return(Rslt);
}


int
CKAligner::ProcAssignMultiMatches(tsClusterThreadPars *pPars)
{
int Rslt;
uint32_t MatchFrom;
uint32_t MatchUntil;
uint32_t HitIdx;
uint32_t Score;
uint32_t ClustHitIdx;
tsReadHit *pCurHit;
tsReadHit *pClustHit;
tsReadHit *pPrevProcCurHit;
int Overlap;
uint32_t ClustEndLoci;

while((Rslt = GetClusterStartEnd(&MatchFrom,&MatchUntil)) > 0)
	{
	pCurHit = &m_pMultiHits[MatchFrom];
	pPrevProcCurHit = nullptr;
	for(HitIdx=MatchFrom;HitIdx <= MatchUntil; HitIdx++,pCurHit++)
		{
		if(!pCurHit->HitLoci.FlagMH)			// only interested in assigning reads which align to multiple hit loci
			continue;							// reads which map to a single hit loci do not require processing

		if(pPrevProcCurHit != nullptr &&
					AdjStartLoci(&pPrevProcCurHit->HitLoci.Hit.Seg[0]) == AdjStartLoci(&pCurHit->HitLoci.Hit.Seg[0]) &&
					AdjHitLen(&pPrevProcCurHit->HitLoci.Hit.Seg[0]) == AdjHitLen(&pCurHit->HitLoci.Hit.Seg[0]) &&
					pPrevProcCurHit->HitLoci.Hit.Seg[0].Strand == pCurHit->HitLoci.Hit.Seg[0].Strand &&
					pPrevProcCurHit->HitLoci.Hit.Seg[0].ChromID == pCurHit->HitLoci.Hit.Seg[0].ChromID)
			{
			pCurHit->HitLoci.Hit.Score = pPrevProcCurHit->HitLoci.Hit.Score;
			continue;
			}

		pPrevProcCurHit = nullptr;
		pCurHit->HitLoci.Hit.Score = 0;
		pClustHit = pCurHit;
		ClustHitIdx = HitIdx;
		while(ClustHitIdx-- > 0)				// checking for clustering upstream of current hit loci
			{
			pClustHit -= 1;

			// can't cluster with reads on a different chrom!
			if(pClustHit->HitLoci.Hit.Seg[0].ChromID != pCurHit->HitLoci.Hit.Seg[0].ChromID)
				break;

			// can't cluster if much too far away
			if((AdjStartLoci(&pCurHit->HitLoci.Hit.Seg[0]) - AdjStartLoci(&pClustHit->HitLoci.Hit.Seg[0])) >= (int)m_MaxReadsLen)
				break;
			// only cluster if >= cClustMultiOverLap
			ClustEndLoci = AdjEndLoci(&pClustHit->HitLoci.Hit.Seg[0]);
			if((int)ClustEndLoci < (AdjStartLoci(&pCurHit->HitLoci.Hit.Seg[0]) + cClustMultiOverLap))
				continue;
			Overlap = min(AdjHitLen(&pCurHit->HitLoci.Hit.Seg[0]),(int)ClustEndLoci - AdjStartLoci(&pCurHit->HitLoci.Hit.Seg[0]));

			if((m_MLMode == eMLuniq && pClustHit->HitLoci.FlagMH) ||
				(pCurHit->HitLoci.Hit.Score & cUniqueClustFlg && (pCurHit->HitLoci.Hit.Score & ~cUniqueClustFlg) >= 0x01fff))
				continue;

			if(pClustHit->HitLoci.Hit.Seg[0].Strand == pCurHit->HitLoci.Hit.Seg[0].Strand && pClustHit->ReadID != pCurHit->ReadID)
				{
				if(!pClustHit->HitLoci.FlagMH)	// clustering to a unique aligned reads has much higher priority than to other multialigned reads
					{
					Score = 1 + (Overlap * cClustUniqueScore)/cClustScaleFact;
					if(pCurHit->HitLoci.Hit.Score & cUniqueClustFlg)
						Score += pCurHit->HitLoci.Hit.Score & ~cUniqueClustFlg;
					if(Score > 0x01fff)			// clamp upstream scores to be no more than 0x01fff so still room for dnstream scores
						Score = 0x01fff;
					pCurHit->HitLoci.Hit.Score = (uint16_t)(Score | cUniqueClustFlg);	// flag that this score is because now clustering to uniquely aligned reads
					if(Score == 0x01fff)
						break;					// this is a good candidate loci
					}
				else							// never seen a uniquely aligned so still clustering to other non-unique aligned reads
					{
					if(!(pCurHit->HitLoci.Hit.Score & cUniqueClustFlg))
						{
						Score = 1 + (Overlap * cClustMultiScore)/cClustScaleFact;
						Score += pCurHit->HitLoci.Hit.Score & ~cUniqueClustFlg;
						if(Score > 0x01fff)			// clamp upstream scores to be no more than 0x01fff so still room for dnstream scores
							Score = 0x01fff;
						pCurHit->HitLoci.Hit.Score = Score;
						}
					}
				}
			}

		// now cluster downstream
		pClustHit = pCurHit;
		ClustHitIdx = HitIdx;
		while(++ClustHitIdx < m_NumMultiHits)				// checking for clustering downstream of current hit loci
			{
			pClustHit += 1;
			// can't cluster with reads on a different chrom!
			if(pClustHit->HitLoci.Hit.Seg[0].ChromID != pCurHit->HitLoci.Hit.Seg[0].ChromID)
				break;

			// can't score if no overlap of at least cClustMultiOverLap
			ClustEndLoci = AdjEndLoci(&pCurHit->HitLoci.Hit.Seg[0]);
			if(AdjStartLoci(&pClustHit->HitLoci.Hit.Seg[0]) > (ClustEndLoci - (uint32_t)cClustMultiOverLap))
				break;

			Overlap = min(AdjHitLen(&pClustHit->HitLoci.Hit.Seg[0]),(int)ClustEndLoci - AdjStartLoci(&pClustHit->HitLoci.Hit.Seg[0]));

			if((m_MLMode == eMLuniq && pClustHit->HitLoci.FlagMH) ||
				(pCurHit->HitLoci.Hit.Score & cUniqueClustFlg && (pCurHit->HitLoci.Hit.Score & ~cUniqueClustFlg) >= 0x03fff))
				continue;

			if(pClustHit->HitLoci.Hit.Seg[0].Strand == pCurHit->HitLoci.Hit.Seg[0].Strand && pClustHit->ReadID != pCurHit->ReadID)
				{
				if(!pClustHit->HitLoci.FlagMH)	// clustering to a unique aligned reads has much higher priority than to other multialigned reads
					{
					Score = 1 + (Overlap * cClustUniqueScore)/cClustScaleFact;
					if(pCurHit->HitLoci.Hit.Score & cUniqueClustFlg)
						Score += pCurHit->HitLoci.Hit.Score & ~cUniqueClustFlg;
					if(Score > 0x3fff)
						Score = 0x3fff;
					pCurHit->HitLoci.Hit.Score = (uint16_t)(Score | cUniqueClustFlg);
					if(Score == 0x3fff)
						break;
					}
				else
					{
					if(!(pCurHit->HitLoci.Hit.Score & cUniqueClustFlg))
						{
						Score = 1 + (Overlap * cClustMultiScore)/cClustScaleFact;
						Score += pCurHit->HitLoci.Hit.Score & ~cUniqueClustFlg;
						if(Score > 0x3fff)
							Score = 0x3fff;
						pCurHit->HitLoci.Hit.Score = (uint16_t)Score;
						}
					}
				}
			}
		pPrevProcCurHit = pCurHit;
		}
	}
pPars->Rslt = 1;
#ifdef _WIN32
_endthreadex(0);
return(eBSFSuccess);
#else
pthread_exit(&pPars->Rslt);
#endif
}

// AssignMultiMatches
// Use clustering (within a sliding window) to determine which matching loci should be assigned to those reads with multiple hits
// This is multipass clustering process as need to cluster multimatches with uniquely matched before clustering with other multimatched.
int
CKAligner::AssignMultiMatches(void) // false to cluster with uniques, true to cluster with multimatches
{
uint32_t HitIdx;
uint32_t ClustHitIdx;
uint32_t Distance;
int NumAssigned;
int NumUnlocated;
tsReadHit *pCurHit;
tsReadHit *pClustHit;

int ClusterUniqueAssigned;
int ClusterAllAssigned;

if(m_MLMode <= eMLrand)		// nothing to do?
	return(0);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Assigning %d reads which aligned to multiple loci to a single loci",m_NumProvMultiAligned);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Sorting...");
SortReadHits(eRSMReadID,false);
m_mtqsort.qsort(m_pMultiHits,m_NumMultiHits,sizeof(tsReadHit),SortMultiHits);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Sorting completed, now clustering...");

RunClusteringThreads(m_NumThreads);

// sort now by ascending ReadID and descending scores
// and assign the read match with the highest score to that read
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Assigning...");
m_mtqsort.qsort(m_pMultiHits,m_NumMultiHits,sizeof(tsReadHit),SortMultiHitReadIDs);
uint32_t CurReadID = 0;
uint32_t BestScore;
uint32_t NxtBestScore;

tsReadHit *pAssign2Read;
NumAssigned = 0;
NumUnlocated = 0;
ClusterUniqueAssigned = 0;
ClusterAllAssigned = 0;
pCurHit = m_pMultiHits;
for(HitIdx=0;HitIdx < m_NumMultiHits; HitIdx++,pCurHit++)
	{
	if(!pCurHit->HitLoci.FlagMH)	 // only interested in reads with multiple hits
		continue;
	if(CurReadID == pCurHit->ReadID) // only interested in the first for each read
		continue;
	CurReadID = pCurHit->ReadID;
	NumAssigned += 1;

	// only interested if best score for read is at least cMHminScore and that score is at least 2x next best score for read
	BestScore = (uint32_t)(pCurHit->HitLoci.Hit.Score & ~cUniqueClustFlg);
	if(BestScore < cMHminScore)
		continue;

	if((pCurHit->HitLoci.Hit.Score & cUniqueClustFlg) == (pCurHit[1].HitLoci.Hit.Score & cUniqueClustFlg))
		{
		NxtBestScore = (uint32_t)(pCurHit[1].HitLoci.Hit.Score & ~cUniqueClustFlg);
		if(BestScore < (NxtBestScore * 2))
			continue;
		}

	// accept this loci as the loci as being the putative read hit loci for this multihit
	pCurHit->HitLoci.FlagMHA = 1;
	if(pCurHit->HitLoci.Hit.Score & cUniqueClustFlg)
		{
		ClusterUniqueAssigned += 1;
		pCurHit->HitLoci.FlagHL = eHLclustunique;
		}
	else
		{
		ClusterAllAssigned += 1;
		pCurHit->HitLoci.FlagHL = eHLclustany;
		}
	}

// need to ensure that as a result of the assignments no assigned multialigned read is now actually an orphan
// orphans are those assigned as eHLclustany if within the window there are no other assigned multireads...
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Checking for orphans (unclustered) from %d putative assignments..",NumAssigned);
m_mtqsort.qsort(m_pMultiHits,m_NumMultiHits,sizeof(tsReadHit),SortMultiHits);
NumUnlocated = 0;
ClusterUniqueAssigned = 0;
ClusterAllAssigned = 0;
int PutativeAssignments = NumAssigned;
bool bAcceptMulti;
NumAssigned = 0;
pCurHit = m_pMultiHits;
for(HitIdx=0;HitIdx < m_NumMultiHits; HitIdx++,pCurHit++)
	{
	if(!pCurHit->HitLoci.FlagMHA)			// only interested in multihit reads which have been assigned to a single loci
		continue;

	bAcceptMulti = false;
	if(pCurHit->HitLoci.FlagHL == eHLclustany) // if was clustered with another multihit then ensure that the other multihit still in play..
		{
		pClustHit = pCurHit;
		ClustHitIdx = HitIdx;
		while(ClustHitIdx-- > 0)				// checking for clustering upstream of current hit loci
			{
			pClustHit -= 1;
			Distance = AdjStartLoci(&pCurHit->HitLoci.Hit.Seg[0]) - AdjStartLoci(&pClustHit->HitLoci.Hit.Seg[0]);
			if(Distance > (uint32_t)(cClustMultiOverLap + (int)AdjHitLen(&pClustHit->HitLoci.Hit.Seg[0]))) // finish if much too far away
				break;
			if(pClustHit->HitLoci.Hit.Seg[0].ChromID != pCurHit->HitLoci.Hit.Seg[0].ChromID)
				break;							// can't cluster with reads on a different chrom!
			if(!pClustHit->HitLoci.FlagMH || pClustHit->HitLoci.FlagMHA == 1)
				{
				bAcceptMulti = true;
				break;
				}
			}

		if(!bAcceptMulti)
			{
			// now cluster downstream
			pClustHit = pCurHit;
			ClustHitIdx = HitIdx;
			while(++ClustHitIdx < m_NumMultiHits)				// checking for clustering downstream of current hit loci
				{
				pClustHit += 1;
				Distance = AdjStartLoci(&pClustHit->HitLoci.Hit.Seg[0]) - AdjStartLoci(&pCurHit->HitLoci.Hit.Seg[0]);
				if(Distance > (uint32_t)(cClustMultiOverLap + (int)AdjHitLen(&pCurHit->HitLoci.Hit.Seg[0])))				// can't score if much too far away
					break;

				if(pClustHit->HitLoci.Hit.Seg[0].ChromID != pCurHit->HitLoci.Hit.Seg[0].ChromID)
					break;
				if(!pClustHit->HitLoci.FlagMH || pClustHit->HitLoci.FlagMHA == 1)
					{
					bAcceptMulti = true;
					break;
					}
				}
			}
		if(!bAcceptMulti)
			pCurHit->HitLoci.FlagMHA = 0;
		}
	else
		bAcceptMulti = true;

	if(bAcceptMulti == true)
		{
		if((pAssign2Read = LocateRead(pCurHit->ReadID))==nullptr)
			{
			if(NumUnlocated++ < 10)
				gDiagnostics.DiagOut(eDLWarn,gszProcName,"A read (ReadID = %d) was not located in call to LocateRead()",pCurHit->ReadID);
			continue;
			}
		pAssign2Read->HitLoci = pCurHit->HitLoci;
		pAssign2Read->NumHits = 1;
		pAssign2Read->NAR = eNARAccepted;
		pAssign2Read->LowHitInstances = 1;
		if(pCurHit->HitLoci.FlagHL == eHLclustunique)
			ClusterUniqueAssigned += 1;
		else
			if(pCurHit->HitLoci.FlagHL == eHLclustany)
				ClusterAllAssigned += 1;
		NumAssigned += 1;
		}
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Clustering completed, removed %d unclustered orphans from %d putative resulting in %d (%d clustered near unique, %d clustered near other multiloci reads) multihit reads accepted as assigned",
					 PutativeAssignments-NumAssigned,PutativeAssignments,NumAssigned,ClusterUniqueAssigned,ClusterAllAssigned);
SortReadHits(eRSMHitMatch,false,true);
return(0);
}

// Median calculation
uint32_t
CKAligner::MedianInsertLen(uint32_t NumInserts,				// number of insert lengths in pInsertLens
				uint32_t *pInsertLens)		// insert lengths
{
uint32_t Idx;
uint32_t Idy;
uint32_t NumUnder;
uint32_t NumOver;
uint32_t *pInsertLen;
uint32_t *pInsertLen1;
uint32_t Median;
uint32_t CurMedian;
uint32_t LoMedian;
uint32_t HiMedian;

if(NumInserts == 0 || pInsertLens == nullptr)
	return(0);
if(NumInserts == 1)
	return(*pInsertLens);
if(NumInserts == 2)
	return((pInsertLens[0] + pInsertLens[1])/2);

LoMedian = 0;
HiMedian = 0;
pInsertLen = pInsertLens;
for(Idx = 0; Idx < NumInserts; Idx++,pInsertLen++)
	{
	if(*pInsertLen < LoMedian || LoMedian == 0)
		LoMedian = *pInsertLen;
	if(*pInsertLen > HiMedian || HiMedian == 0)
		HiMedian = *pInsertLen;
	}
if(HiMedian == LoMedian)
	return(HiMedian);

pInsertLen = pInsertLens;
Median = *pInsertLen++;
for(Idx = 1; Idx < NumInserts; Idx++,pInsertLen++)
	{
	CurMedian = *pInsertLen;
	if(CurMedian < LoMedian || CurMedian > HiMedian)
		continue;
	
	NumUnder = 0;
	NumOver = 0;
	pInsertLen1 = pInsertLens;
	for(Idy = 0; Idy < NumInserts; Idy++,pInsertLen1++)
		{
		if(Idy == Idx || *pInsertLen1 == CurMedian)
			continue;
		if(*pInsertLen1 < CurMedian)
			NumUnder += 1;
		else
			NumOver += 1;
		}
	if(NumUnder == NumOver)
		return(CurMedian);
	if(NumUnder < NumOver)
		LoMedian = CurMedian;
	else
		HiMedian = CurMedian;
	}
return((LoMedian + HiMedian)/2);
}

// Report PE insert length distributions for each transcript or assembly contig/sequence
int
CKAligner::ReportPEInsertLenDist(void)
{
uint32_t PE1Start;
uint32_t PE2Start;
uint32_t PE1Len;
uint32_t PE2Len;
uint32_t TLen;
uint32_t LenDist[53];				// holds counts for insert lengths from 100 to 600 in 10bp increments plus under 100bp and over 600bp
char szDistBuff[0x7fff];
int BuffOfs;
int Idx;
uint32_t NumPEs;
uint64_t SumInsertLens;
uint32_t MedianInsert;
uint32_t CurTargID;
uint32_t NumTargIDs;
char szTargChromName[128];
uint32_t TargChromLen;
tsReadHit *pPE1ReadHit;
tsReadHit *pPE2ReadHit;
uint32_t *pInsertLens;

if(m_hInsertLensFile == -1 || !m_bPEInsertLenDist)		// only process PE insert distributions if actually aligning PEs and user requested the insert distributions
	return(0);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Reporting PE insert lengths for each transcript or assembly sequence, sorting reads");
SortReadHits(eRSMPEHitMatch,false);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Completed sort");

pInsertLens = new uint32_t [1000000 + 10];	// big assumption that target transcript or chromosome will have no more than 1M paired end hits - check is made and only 1st 1M inserts are processed!

BuffOfs = sprintf(szDistBuff,"\"TargSeq\",\"TargLen\",\"TotPEs\",\"MedianInsertLen\",\"MeanInsertLen\",\"Insert <100bp\"");
for(Idx = 0; Idx < 50; Idx++)
	BuffOfs += sprintf(&szDistBuff[BuffOfs],",\"Insert %dbp\"",100 + (10 * Idx));
BuffOfs += sprintf(&szDistBuff[BuffOfs],",\"Insert >600bp\"\n");
CUtility::RetryWrites(m_hInsertLensFile,szDistBuff,BuffOfs);
BuffOfs = 0;

pPE1ReadHit = nullptr;
CurTargID = 0;
NumTargIDs = 0;
NumPEs = 0;
SumInsertLens = 0;
memset(LenDist,0,sizeof(LenDist));
memset(pInsertLens,0,sizeof(uint32_t) * 1000001);
while((pPE1ReadHit = IterSortedReads(pPE1ReadHit))!=nullptr)
	{
	if((pPE2ReadHit = IterSortedReads(pPE1ReadHit)) == nullptr) // expecting PE2 to immediately follow PE1
		break;

	// both ends of pair must be accepted as aligned and both must be aligned as PE
	if(!(pPE1ReadHit->NAR == eNARAccepted && pPE1ReadHit->FlgPEAligned &&  pPE2ReadHit->NAR == eNARAccepted && pPE2ReadHit->FlgPEAligned))
		break;

	// both pPE1ReadHit and pPE2ReadHit have been accepted as being PE aligned
	if(pPE1ReadHit->HitLoci.Hit.Seg[0].ChromID != (uint32_t)CurTargID)  // now processing a different transcript or assembly sequence?
		{
		if(NumTargIDs && NumPEs)
			{
			MedianInsert = MedianInsertLen(min(1000000u,NumPEs),pInsertLens);
			BuffOfs += sprintf(&szDistBuff[BuffOfs],"\"%s\",%u,%u,%u,%u",szTargChromName,TargChromLen,NumPEs,MedianInsert,(uint32_t)(SumInsertLens/NumPEs));
			for(Idx = 0; Idx < 52; Idx++)
				BuffOfs += sprintf(&szDistBuff[BuffOfs],",%d",LenDist[Idx]);
			BuffOfs += sprintf(&szDistBuff[BuffOfs],"\n");
			if((BuffOfs + 4096) > sizeof(szDistBuff))
				{
				CUtility::RetryWrites(m_hInsertLensFile,szDistBuff,BuffOfs);
				BuffOfs = 0;
				}
			}

		CurTargID = pPE1ReadHit->HitLoci.Hit.Seg[0].ChromID;
		m_pSfxArray->GetIdentName(CurTargID,sizeof(szTargChromName),szTargChromName);
		TargChromLen = m_pSfxArray->GetSeqLen(CurTargID);
		NumTargIDs += 1;
		NumPEs = 0;
		SumInsertLens = 0;
		memset(LenDist,0,sizeof(LenDist));
		memset(pInsertLens,0,sizeof(uint32_t) * 1000001);
		}

	if(NumPEs >= 1000000)
		{
		pPE1ReadHit = pPE2ReadHit;
		continue;
		}

	PE1Start = AdjAlignStartLoci(&pPE1ReadHit->HitLoci.Hit) + 1;
	PE2Start = AdjAlignStartLoci(&pPE2ReadHit->HitLoci.Hit) + 1;
	PE1Len = AdjAlignHitLen(&pPE1ReadHit->HitLoci.Hit);
	PE2Len = AdjAlignHitLen(&pPE2ReadHit->HitLoci.Hit);
	if(PE1Start <= PE2Start)
		TLen = (PE2Start - PE1Start) + PE2Len;
	else
		TLen = (PE1Start - PE2Start) + PE1Len;

	SumInsertLens += TLen;
	if(TLen < 100)
		LenDist[0] += 1;
	else
		{
		if(TLen > 600)
			LenDist[51] += 1;
		else
			LenDist[1 + ((TLen - 99)/10)] += 1;
		}
	if(NumPEs < 1000000)
		pInsertLens[NumPEs] = TLen;
	NumPEs += 1;
	pPE1ReadHit = pPE2ReadHit;
	}

if(NumTargIDs)
	{
	if(NumPEs)
		{
		MedianInsert = MedianInsertLen(min((uint32_t)1000000,NumPEs),pInsertLens);
		BuffOfs += sprintf(&szDistBuff[BuffOfs],"\"%s\",%u,%u,%u,%u",szTargChromName,TargChromLen,NumPEs,MedianInsert,(uint32_t)(SumInsertLens/NumPEs));
		for(Idx = 0; Idx < 52; Idx++)
			BuffOfs += sprintf(&szDistBuff[BuffOfs],",%u",LenDist[Idx]);
		BuffOfs += sprintf(&szDistBuff[BuffOfs],"\n");
		}

	if(BuffOfs)
		CUtility::RetryWrites(m_hInsertLensFile,szDistBuff,BuffOfs);
	}

delete []pInsertLens;
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Reporting PE insert lengths for %u transcripts or assembly sequences completed",NumTargIDs);
SortReadHits(eRSMHitMatch,false);
return(NumTargIDs);
}

// Report number of accepted alignments onto each targeted transcript or assembly contig/sequence
int
CKAligner::ReportTargHitCnts(void)
{
char szDistBuff[0x7fff];
char szTrimer[4];
uint32_t TriMerCnts[4*4*4];	// cnts for each possible initial starting tri-mer; if containing an indeterminate then counted separately
uint32_t Indeterminates;	// cnts of initial starting tri-mers which contain at least 1 indeterminate

bool bIsIndeterminate;

int TriIdx;
int SeqIdx;
int Base;
uint8_t* pSeqVal;

uint32_t NumLociUniques;
uint32_t CurLociUnique;
uint32_t LociUnique;

int BuffOfs;
uint32_t NumHits;
uint32_t CurTargID;
uint32_t NumTargIDs;
char szTargChromName[128];
uint32_t TargChromLen;
tsReadHit *pReadHit;

if(m_hStatsFile == -1)
	return(0);

int hTargHitCnts;
char szTargHitCnts[_MAX_PATH];
int BuffIdx = 0;

CUtility::AppendFileNameSuffix(szTargHitCnts, m_pszStatsFile, (char*)".AlignCntsDist.csv", '.');
#ifdef _WIN32
hTargHitCnts = open(szTargHitCnts, (O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC), (_S_IREAD | _S_IWRITE));
#else
if ((hTargHitCnts = open(szTargHitCnts, O_WRONLY | O_CREAT, S_IREAD | S_IWRITE)) != -1)
if (ftruncate(hTargHitCnts, 0) != 0)
{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to truncate global PE insert size distribution file '%s' - %s", szTargHitCnts, strerror(errno));
	return(eBSFerrCreateFile);
}
#endif

if (hTargHitCnts < 0)
{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to truncate global PE insert size distribution file '%s'", szTargHitCnts);
	return(eBSFerrCreateFile);
}


gDiagnostics.DiagOut(eDLInfo,gszProcName,"Reporting accepted read alignment counts on to targeted transcripts or sequences, sorting reads");
SortReadHits(eRSMHitMatch,false);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Completed sort");

BuffOfs = sprintf(szDistBuff,"\"FeatID\",\"TargSeq\",\"TargLen\",\"NumHits\",\"RPKM\",\"NumUniqueLoci\"");
for (TriIdx = 0; TriIdx < 64; TriIdx++)
	{
	int TrimerMsk;
	int TrimerIdx;
	for(TrimerIdx = 0, TrimerMsk = 0x30; TrimerIdx < 3; TrimerIdx += 1, TrimerMsk >>= 2)
		{
		switch((TriIdx & TrimerMsk) >> (2 * (2 - TrimerIdx))) {
			case 0x00:
				szTrimer[TrimerIdx] = 'A';
				break;
			case 0x01:
				szTrimer[TrimerIdx] = 'C';
				break;
			case 0x02:
				szTrimer[TrimerIdx] = 'G';
				break;
			default:
				szTrimer[TrimerIdx] = 'T';
				break;
			}
		}
	szTrimer[3] = '\0';
	BuffOfs += sprintf(&szDistBuff[BuffOfs], ",\"%s\"",szTrimer);
	}
BuffOfs += sprintf(&szDistBuff[BuffOfs],",Indeterminates\n");

if(!CUtility::RetryWrites(hTargHitCnts,szDistBuff,BuffOfs))
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"RetryWrites() error");
	return(eBSFerrWrite);
	}
BuffOfs = 0;

pReadHit = nullptr;
CurTargID = 0;
uint32_t ExpTargID = 0;
TargChromLen = 0;
NumTargIDs = 0;
NumHits = 0;
memset(TriMerCnts,0,sizeof(TriMerCnts));
Indeterminates = 0;
TriIdx = 0;
NumLociUniques = 0;
int64_t NumReads = 0;
int64_t NumReadsAligned = 0;
double RPKM = 0.0;

// how many reads in total? Needed for RPKM calculations
while ((pReadHit = IterSortedReads(pReadHit)) != nullptr)
	NumReads += 1;

pReadHit = nullptr;
while((pReadHit = IterSortedReads(pReadHit))!=nullptr)
	{
	// read must have been accepted as aligned
	if(pReadHit->NAR != eNARAccepted)
		continue;

	// read accepted as aligned
	if(pReadHit->HitLoci.Hit.Seg[0].ChromID != (uint32_t)CurTargID)  // now processing a different transcript or assembly sequence?
		{
		if(CurTargID != 0)
			{
			RPKM = ((double)NumHits * 1000.0f);
			RPKM /= (double)TargChromLen;
			RPKM *= 1000000.0f / (double)NumReads;

			BuffOfs += sprintf(&szDistBuff[BuffOfs], "%u,\"%s\",%u,%u,%f,%u", CurTargID,szTargChromName, TargChromLen, NumHits, RPKM, NumLociUniques);
			for (TriIdx = 0; TriIdx < 64; TriIdx++)
				BuffOfs += sprintf(&szDistBuff[BuffOfs], ",%1.4f", NumHits == 0 ? 0.0 : (double)TriMerCnts[TriIdx] / (double)NumHits);
			BuffOfs += sprintf(&szDistBuff[BuffOfs], ",%u\n", Indeterminates);
			if (((size_t)BuffOfs + 4096) > sizeof(szDistBuff))
				{
				if(!CUtility::RetryWrites(hTargHitCnts, szDistBuff, BuffOfs))
					{
					gDiagnostics.DiagOut(eDLInfo,gszProcName,"RetryWrites() error");
					return(eBSFerrWrite);
					}
				BuffOfs = 0;
				}
			ExpTargID = CurTargID + 1;
			}
		else
			ExpTargID = 1;
		CurTargID = pReadHit->HitLoci.Hit.Seg[0].ChromID;
		while(CurTargID != ExpTargID) // checking for targets w/o hits and reporting these as having no hits
			{
			m_pSfxArray->GetIdentName(ExpTargID, sizeof(szTargChromName), szTargChromName);
			TargChromLen = m_pSfxArray->GetSeqLen(ExpTargID);
			BuffOfs += sprintf(&szDistBuff[BuffOfs], "%u,\"%s\",%u,0,0.0,0", ExpTargID,szTargChromName, TargChromLen);
			for (TriIdx = 0; TriIdx < 64; TriIdx++)
				BuffOfs += sprintf(&szDistBuff[BuffOfs], ",0.0");
			BuffOfs += sprintf(&szDistBuff[BuffOfs], ",0\n");
			if (((size_t)BuffOfs + 4096) > sizeof(szDistBuff))
				{
				if(!CUtility::RetryWrites(hTargHitCnts, szDistBuff, BuffOfs))
					{
					gDiagnostics.DiagOut(eDLInfo,gszProcName,"RetryWrites() error");
					return(eBSFerrWrite);
					}
				BuffOfs = 0;
				}
			ExpTargID += 1;
			}
		// CurTargID will always be same as ExpTargID by now

		m_pSfxArray->GetIdentName(CurTargID,sizeof(szTargChromName),szTargChromName);
		TargChromLen = m_pSfxArray->GetSeqLen(CurTargID);
		NumTargIDs += 1;
		NumHits = 0;
		NumLociUniques = 0;
		memset(TriMerCnts,0,sizeof(TriMerCnts));
		Indeterminates = 0;
		}


	bIsIndeterminate = false;
	pSeqVal = &pReadHit->Read[pReadHit->DescrLen + 1];
	TriIdx = 0;
	for (SeqIdx = 0; SeqIdx < 3; SeqIdx++, pSeqVal++)
		{
		Base = (*pSeqVal & 0x07);		// not interested in any soft masking within a SAM sequence
		if(Base > eBaseT)
			{
			bIsIndeterminate = true;
			break;
			}
		TriIdx <<= 2;
		TriIdx |= (Base&0x03);
		}
	if(bIsIndeterminate)
		Indeterminates += 1;
	else
		TriMerCnts[TriIdx] += 1;
	LociUnique = AdjAlignStartLoci(&pReadHit->HitLoci.Hit);
	if(NumLociUniques == 0 || LociUnique != CurLociUnique)
		{
		NumLociUniques += 1;
		CurLociUnique = LociUnique;
		}
	NumHits += 1;
	}

if(NumHits)
	{
	RPKM = ((double)NumHits * 1000.0f);
	RPKM /= (double)TargChromLen;
	RPKM *= 1000000.0f / (double)NumReads;
	BuffOfs += sprintf(&szDistBuff[BuffOfs],"%u,\"%s\",%u,%u,%f,%u", CurTargID,szTargChromName,TargChromLen,NumHits, RPKM,NumLociUniques);
	for(TriIdx = 0; TriIdx < (4*4*4); TriIdx++)
		BuffOfs += sprintf(&szDistBuff[BuffOfs], ",%1.4f", NumHits == 0 ? 0.0 : (double)TriMerCnts[TriIdx]/(double)NumHits);
	BuffOfs += sprintf(&szDistBuff[BuffOfs], ",%u\n", Indeterminates);
	if(!CUtility::RetryWrites(hTargHitCnts,szDistBuff,BuffOfs))
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"RetryWrites() error");
		return(eBSFerrWrite);
		}
	BuffOfs = 0;
	}
ExpTargID = m_pSfxArray->GetNumEntries();
while (CurTargID != ExpTargID) // checking for targets w/o hits and reporting these as having no hits
	{
	m_pSfxArray->GetIdentName(CurTargID+1, sizeof(szTargChromName), szTargChromName);
	TargChromLen = m_pSfxArray->GetSeqLen(CurTargID+1);
	BuffOfs += sprintf(&szDistBuff[BuffOfs], "%u,\"%s\",%u,0,0.0,0", CurTargID + 1,szTargChromName, TargChromLen);
	for (TriIdx = 0; TriIdx < 64; TriIdx++)
		BuffOfs += sprintf(&szDistBuff[BuffOfs], ",0");
	BuffOfs += sprintf(&szDistBuff[BuffOfs], ",0\n");
	if (((size_t)BuffOfs + 4096) > sizeof(szDistBuff))
		{
		if(!CUtility::RetryWrites(hTargHitCnts, szDistBuff, BuffOfs))
			{
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"RetryWrites() error");
			return(eBSFerrWrite);
			}
		BuffOfs = 0;
		}
	CurTargID += 1;
	}
if (BuffOfs)
	{
	if(!CUtility::RetryWrites(hTargHitCnts, szDistBuff, BuffOfs))
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"RetryWrites() error");
		return(eBSFerrWrite);
		}
	BuffOfs = 0;
	}

#ifdef _WIN32
_commit(hTargHitCnts);
#else
fsync(hTargHitCnts);
#endif
close(hTargHitCnts);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Completed reporting read alignment counts on to %d targeted transcripts or sequences",NumTargIDs);
return(NumTargIDs);
}



// Write results as BAM or SAM format
int
CKAligner::WriteBAMReadHits(etFMode ProcMode,	   // eFMsam or eFMsamAll
							teSAMFormat SAMFormat, // if SAM output format then could be SAM or BAM compressed dependent on the file extension used
							bool bPEProc,		   // true if processing paired ends
							int ComprLev)		   // BAM to be BGZF compressed at the requested level (0..9)
{
int Rslt;
tsBAMalign *pBAMalign;			// to hold each SAM or BAM alignment as it is constructed
eSAMFileType FileType;
CSAMfile *pSAMfile;

char szChromName[128];
int NumAlignedToSeqs;
int NumSeqsInHdr;
uint32_t NumReportedBAMreads;
uint32_t PrevNumReportedBAMreads;
tsReadHit *pReadHit;
tBSFEntryID PrevTargEntry;
int NumChroms;
bool bRptAllChroms;
uint32_t ChromSeqLen;

int ChromID;
uint32_t CurChromID;
uint16_t EntryFlags;

if((pBAMalign = new tsBAMalign) == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "WriteBAMReadHits: Unable to allocate memory for instance of tsBAMalign");
	return(eBSFerrInternal);
	}


if((pSAMfile = new CSAMfile) == nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"WriteBAMReadHits: Unable to instantiate class CSAMfile");
	return(eBSFerrInternal);
	}

switch(SAMFormat) {
	case etSAMFformat:			// output SAM
		if(m_bgzOutFile)
			FileType = eSFTSAMgz;
		else
			FileType = eSFTSAM;
		break;
	case etSAMFBAM:				// output as BAM compressed with bgzf
		FileType = eSFTBAM_BAI;
		break;
	}

if((Rslt = pSAMfile->Create(FileType,m_pszOutFile,ComprLev,(char *)kit4bversion)) < eBSFSuccess)
	{
	delete pBAMalign;
	delete pSAMfile;
	return(Rslt);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Sorting alignments by ascending chrom.loci");
SortReadHits(eRSMHitMatch,false);

// mark entries (chroms/contigs/sequences) for which there is at least one alignment so only these marked entries
// will be written to the SAM or BAM header
pReadHit = nullptr;
PrevTargEntry = 0;
m_PrevSAMTargEntry = 0;
CurChromID = 0;
NumChroms = m_pSfxArray->GetNumEntries();
bRptAllChroms = m_MaxRptSAMSeqsThres >= NumChroms ? true : false;

// identify chroms to be reported
// only reporting those which have accepted alignments unless reporting all chromosomes even if not all have alignments
while((pReadHit = IterSortedReads(pReadHit))!=nullptr)
	{
	if(pReadHit->NAR == eNARAccepted)
		{
		if(CurChromID == 0 || CurChromID != pReadHit->HitLoci.Hit.Seg[0].ChromID)
			{
			CurChromID = pReadHit->HitLoci.Hit.Seg[0].ChromID;
			if(CurChromID > 0)
				m_pSfxArray->SetResetIdentFlags(CurChromID,0x01,0x00);
			}
		}
	}

NumAlignedToSeqs = 0;
NumSeqsInHdr = 0;
for(ChromID = 1; ChromID <= NumChroms; ChromID++)
	{
	EntryFlags = m_pSfxArray->GetIdentFlags(ChromID);
	if(EntryFlags & 0x01 || bRptAllChroms)
		{
		m_pSfxArray->GetIdentName(ChromID,sizeof(szChromName),szChromName);
		ChromSeqLen=m_pSfxArray->GetSeqLen(ChromID);
		if((Rslt = pSAMfile->AddRefSeq(m_szTargSpecies,szChromName,ChromSeqLen)) < 1)
			{
			delete pBAMalign;
			delete pSAMfile;
			return(Rslt);
			}
		NumSeqsInHdr += 1;
		if(EntryFlags & 0x01)
			NumAlignedToSeqs += 1;
		}
	}
pSAMfile->StartAlignments();

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Header written with references to %d sequences of which %d have at least 1 alignments",NumSeqsInHdr, NumAlignedToSeqs);
NumAlignedToSeqs = 0;
NumSeqsInHdr = 0;

// now write out each alignment in BAM format
pReadHit = nullptr;

PrevTargEntry = 0;
m_PrevSAMTargEntry = 0;
PrevNumReportedBAMreads = 0;
NumReportedBAMreads = 0;
int RefID;
int BAMRefID;
tsReadHit *pNxtRead; 
bool bLastAligned;
int  ReadIs;				// 0 if SE, 1 if PE1 of a PE, 2 if PE2 of a PE


gDiagnostics.DiagOut(eDLInfo,gszProcName,"Reported %s %u read alignments",SAMFormat == etSAMFformat ? "SAM" : "BAM",NumReportedBAMreads);
RefID = 0;
time_t Started = time(0);
while((pReadHit = IterSortedReads(pReadHit))!=nullptr)
	{
	if(pReadHit->NAR == eNARAccepted || ProcMode == eFMsamAll)
		{
		if(!bPEProc)
			ReadIs = 0;
		else
			ReadIs = pReadHit->PairReadID & 0x080000000 ? 0x02 : 0x01;

		if(pReadHit->NAR == eNARAccepted)
			{
			if(pReadHit->HitLoci.Hit.Seg[0].ChromID != (uint32_t)m_PrevSAMTargEntry)
				{
				m_pSfxArray->GetIdentName(pReadHit->HitLoci.Hit.Seg[0].ChromID,sizeof(m_szSAMTargChromName),m_szSAMTargChromName);
				m_PrevSAMTargEntry = pReadHit->HitLoci.Hit.Seg[0].ChromID;
				}
			BAMRefID = 0;
			}
		else   // else also reporting reads not accepted as being aligned
			{
			BAMRefID = -1;
			m_szSAMTargChromName[0] = '*';
			m_szSAMTargChromName[1] = '\0';
			}

		if ((Rslt = ReportBAMread(pReadHit, BAMRefID, ReadIs,pBAMalign)) < eBSFSuccess)
			{
			delete pBAMalign;
			return(Rslt);
			}
		strcpy(pBAMalign->szRefSeqName,m_szSAMTargChromName);

		// look ahead to check if current read is the last accepted aligned read
		bLastAligned = false;
		if(pReadHit->NAR == eNARAccepted)
			{
			pNxtRead = IterSortedReads(pReadHit);
			if(pNxtRead == nullptr || pNxtRead->NAR != eNARAccepted)
				bLastAligned = true;
			}

		if((Rslt = pSAMfile->AddAlignment(pBAMalign,bLastAligned)) < eBSFSuccess)
			return(Rslt);

		NumReportedBAMreads += 1;

		if(NumReportedBAMreads > (PrevNumReportedBAMreads + 50000))
			{
			time_t Now = time(0);
			unsigned long ElapsedSecs = (unsigned long) (Now - Started);
			if(ElapsedSecs >= 60)
				{
				gDiagnostics.DiagOut(eDLInfo,gszProcName,"Reported %s %u read alignments",SAMFormat == etSAMFformat ? "SAM" : "BAM",NumReportedBAMreads);
				Started = Now;
				}
			PrevNumReportedBAMreads = NumReportedBAMreads;
			}
		// user may be interested in the distribution of the aligner induced substitutions, after any auto-trimming of flanks,
		// along the length of the reads and how this distribution relates to the quality scores
		if(m_hStatsFile != -1)
			WriteSubDist(pReadHit);
		}
	}

delete pBAMalign;
pSAMfile->Close();
delete pSAMfile;
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Completed reporting %s %u read alignments",SAMFormat == etSAMFformat ? "SAM" : "BAM",NumReportedBAMreads);
return(0);
}

// BAM index
// uint8_t magic[4];    // "BAI\1"
// uint32_t n_rf;       // number of reference sequences following
//    uint32_t n_bin;   // number of distinct bins for current reference sequence
//        uint32_t bin; // distinct bin
//        uint32_t chunks; // number of chunks following
//            uint64_t chumk_beg;		// virtual file offset at which chunk starts
//            uint64_t chumk_end;		// virtual file offset at which chunk ends
//    uint32_t n_intv;  // number of 16kb intervals for linear index
//        uint64_t ioffset;   // virtual file offset of first alignment in interval


// following BAM bin functions are copied from the specification at http://samtools.sourceforge.net/SAMv1.pdf
/* calculate bin given an alignment covering [beg,end) (zero-based, half-close-half-open) */
int CKAligner::BAMreg2bin(int beg, int end)
{
--end;
if (beg>>14 == end>>14) return ((1<<15)-1)/7 + (beg>>14);
if (beg>>17 == end>>17) return ((1<<12)-1)/7 + (beg>>17);
if (beg>>20 == end>>20) return ((1<<9)-1)/7 + (beg>>20);
if (beg>>23 == end>>23) return ((1<<6)-1)/7 + (beg>>23);
if (beg>>26 == end>>26) return ((1<<3)-1)/7 + (beg>>26);
return 0;
}

/* calculate the list of bins that may overlap with region [beg,end) (zero-based) */
#define MAX_BIN (((1<<18)-1)/7)
int CKAligner::BAMreg2bins(int beg, int end, uint16_t *plist)
{
int i = 0, k;
--end;
plist[i++] = 0;
for (k = 1 + (beg>>26); k <= 1 + (end>>26); ++k) plist[i++] = k;
for (k = 9 + (beg>>23); k <= 9 + (end>>23); ++k) plist[i++] = k;
for (k = 73 + (beg>>20); k <= 73 + (end>>20); ++k) plist[i++] = k;
for (k = 585 + (beg>>17); k <= 585 + (end>>17); ++k) plist[i++] = k;
for (k = 4681 + (beg>>14); k <= 4681 + (end>>14); ++k) plist[i++] = k;
return i;
}

int
CKAligner::ReportBAMread(tsReadHit *pReadHit,	// read to report
			int RefID,							// read aligns to this BAM refID (-1) if unaligned
			int  ReadIs,						// 0 if SE, 1 if PE1 of a PE, 2 if PE2 of a PE
			tsBAMalign *pBAMalign)				// BAM alignment to return
{
tsReadHit *pPEReadHit;
char szSEQName[128];
char *pszSEQName;
char *pszPEQName;
char *pszQName;
char *pszRNext;

int SeqIdx;
etSeqBase Sequence[cMaxFastQSeqLen+ 10000];	// to hold sequence (sans quality scores) for current read

int SumScores;
uint8_t *pSeqVal;
etSeqBase *pSeq;
char *pQScore;
char ExchScore;

int Flags;
int MAPQ;
int GapLen;
int QNameLen;

int CigarIdx;
int Seg0RightTrimLen;
int Seg0LeftTrimLen;
int Seg0Hitlen;
int Seg1Hitlen;
int Seg1RightTrimLen;
int Seg1LeftTrimLen;

int SEStart;
int PEStart;
int SELen;
int PELen;
int PNext;
int TLen;

int LimitGapErrs;
if(pReadHit == nullptr)
	return(eBSFerrInternal);

pBAMalign->block_size = 0; 	
pBAMalign->refID = 0; 
pBAMalign->pos = 0; 
pBAMalign->end = 0; 
pBAMalign->bin_mq_nl = 0; 
pBAMalign->flag_nc = 0; 
pBAMalign->l_seq = 0; 
pBAMalign->next_refID = 0; 
pBAMalign->next_pos = 0; 
pBAMalign->tlen = 0; 
pBAMalign->NumReadNameBytes = 0; 
pBAMalign->NumCigarBytes = 0; 
pBAMalign->NumSeqBytes = 0; 
pBAMalign->NumAux = 0; 
pBAMalign->szRefSeqName[0] = 0; 
pBAMalign->szMateRefSeqName[0] = 0; 
pBAMalign->read_name[0] = 0;
pBAMalign->cigar[0] = 0;
pBAMalign->seq[0] = 0;
pBAMalign->qual[0] = 0;
pBAMalign->auxData[0].array_type = 0; 
pBAMalign->auxData[0].NumVals = 0; 
pBAMalign->auxData[0].tag[0] = 0;
pBAMalign->auxData[0].value[0] = 0;
pBAMalign->auxData[0].val_type = 0;

memcpy(szSEQName,(char *)pReadHit->Read,pReadHit->DescrLen+1);
pszSEQName = szSEQName;
pszQName = pszSEQName;

pPEReadHit = nullptr;
Flags = 0;
SEStart = 0;
SELen = 0;
PNext = 0;
TLen = 0;
pszRNext = (char *)"*";
pszPEQName = nullptr;

switch(ReadIs) {
	case 0:			// reads have been processed as SE reads
		if(pReadHit->NAR == eNARAccepted)
			{	
			Flags = pReadHit->HitLoci.Hit.Seg[0].Strand == '+' ? 0 : cSAMFlgAS;	// accepted as aligned, flag which strand alignment was to
			SEStart = AdjAlignStartLoci(&pReadHit->HitLoci.Hit);
			SELen = AdjAlignHitLen(&pReadHit->HitLoci.Hit);
			}
		else
			Flags = cSAMFlgUnmapped;		// read was unaligned
		break;

	case 1:			// PE1 of a PE pair
		Flags = cSAMFlgReadPaired | cSAMFlgReadPairMap | cSAMFlgPE1;
		if(pReadHit->NAR == eNARAccepted)
			{
			Flags |= pReadHit->HitLoci.Hit.Seg[0].Strand == '+' ? 0 : cSAMFlgAS;	// accepted as aligned, flag which strand alignment was to
			SEStart = AdjAlignStartLoci(&pReadHit->HitLoci.Hit);
			SELen = AdjAlignHitLen(&pReadHit->HitLoci.Hit);
			}
		else
			Flags |= cSAMFlgUnmapped;

		pPEReadHit = (tsReadHit *)((uint8_t *)pReadHit + sizeof(tsReadHit) + pReadHit->ReadLen + pReadHit->DescrLen);
		if(pReadHit->FlgPEAligned && pPEReadHit->FlgPEAligned && pPEReadHit->NAR == eNARAccepted)
			{
			Flags |= pPEReadHit->HitLoci.Hit.Seg[0].Strand == '+' ? 0 : cSAMFlgMateAS;	// accepted as aligned, flag which strand alignment was to

			if(pReadHit->NAR == eNARAccepted)
				{
				pszRNext = (char *)"=";
				PEStart = AdjAlignStartLoci(&pPEReadHit->HitLoci.Hit);
				PELen = AdjAlignHitLen(&pPEReadHit->HitLoci.Hit);

				if(SEStart <= PEStart)
					TLen = (PEStart - SEStart) + PELen;
				else
					TLen = (SEStart - PEStart) + SELen;
				}
			}
		else
			Flags |= cSAMFlgMateUnmapped;
		break;

	case 2:			// this read is PE2 of a PE
		Flags = cSAMFlgReadPaired | cSAMFlgReadPairMap | cSAMFlgPE2;
		if(pReadHit->NAR == eNARAccepted)
			{
			Flags |= pReadHit->HitLoci.Hit.Seg[0].Strand  == '+' ? 0 : cSAMFlgAS;	// accepted as aligned, flag which strand alignment was to
			SEStart = AdjAlignStartLoci(&pReadHit->HitLoci.Hit);
			SELen = AdjAlignHitLen(&pReadHit->HitLoci.Hit);
			}
		else
			Flags |= cSAMFlgUnmapped;

		pPEReadHit = (tsReadHit *)((uint8_t *)pReadHit - pReadHit->PrevSizeOf);
		if(pReadHit->FlgPEAligned && pPEReadHit->FlgPEAligned && pPEReadHit->NAR == eNARAccepted)
			{
			Flags |= pPEReadHit->HitLoci.Hit.Seg[0].Strand  == '+' ? 0 : cSAMFlgMateAS;	// accepted as aligned, flag which strand alignment was to
			if(pReadHit->NAR == eNARAccepted)
				{
				pszRNext = (char *)"=";
				PEStart = AdjAlignStartLoci(&pPEReadHit->HitLoci.Hit);
				PELen = AdjAlignHitLen(&pPEReadHit->HitLoci.Hit);
				if(SEStart <= PEStart)
					TLen = (PEStart - SEStart) + PELen;
				else
					TLen = (SEStart - PEStart) + SELen;		
				}
			}
		else
			Flags |= cSAMFlgMateUnmapped;
		break;
	}

LimitGapErrs = 10;

pSeqVal = &pReadHit->Read[pReadHit->DescrLen+1];
pSeq = Sequence;
pQScore = (char *)pBAMalign->qual;
SumScores = 0;
for(SeqIdx = 0; SeqIdx < pReadHit->ReadLen; SeqIdx++,pSeq++,pQScore++,pSeqVal++)
	{
	*pSeq = (*pSeqVal & 0x07);		// not interested in any soft masking within a SAM sequence
	SumScores += (*pSeqVal >> 4) & 0x0f;	// any quality scores would be in bits 4..7
	*pQScore = (char)(33 + ((((*pSeqVal >> 4) & 0x0f) * 40))/15);
	}

if(SumScores == 0)		// if there were no associated quality scores
	memset(pBAMalign->qual,0x0ff,pReadHit->ReadLen);
else
	{
	*pQScore = '\0';
	if(pReadHit->NAR == eNARAccepted && pReadHit->HitLoci.Hit.Seg[0].Strand != '+')	// sequence is reversed so quality scores also should be
		{
		pQScore -= 1;
		for(SeqIdx = 0; SeqIdx < (pReadHit->ReadLen/2); SeqIdx++,pQScore--)
			{
			ExchScore = (char)pBAMalign->qual[SeqIdx];
			pBAMalign->qual[SeqIdx] = *pQScore;
			*pQScore = ExchScore;
			}
		}
	}

MAPQ = 254;		// SAM recommendation is that no mapping quality should be the maximum of 255 ????

CigarIdx = 0;
if(pReadHit->NAR == eNARAccepted)
	{
	Seg0RightTrimLen = pReadHit->HitLoci.Hit.Seg[0].TrimRight;
	Seg0LeftTrimLen = pReadHit->HitLoci.Hit.Seg[0].TrimLeft;
	Seg0Hitlen = AdjHitLen(&pReadHit->HitLoci.Hit.Seg[0]);
	Seg1Hitlen = 0;
	Seg1RightTrimLen = 0;
	Seg1LeftTrimLen = 0;

	if(pReadHit->HitLoci.Hit.Seg[0].Strand == '+')
		{
		if(Seg0LeftTrimLen != 0)
			pBAMalign->cigar[CigarIdx++] = Seg0LeftTrimLen << 4 | 4;		// 'S'
		pBAMalign->cigar[CigarIdx++] = Seg0Hitlen << 4 | 0;					// 'M'
		if(Seg0RightTrimLen > 0)
			pBAMalign->cigar[CigarIdx++] = Seg0RightTrimLen << 4 | 4;		// 'S'
		}
	else
		{
		if(Seg0RightTrimLen != 0)
			pBAMalign->cigar[CigarIdx++] = Seg0RightTrimLen << 4 | 4;		// 'S'
		pBAMalign->cigar[CigarIdx++] = Seg0Hitlen << 4 | 0;					// 'M'	
		if(Seg0LeftTrimLen > 0)
			pBAMalign->cigar[CigarIdx++] = Seg0LeftTrimLen << 4 | 4;		// 'S'
		}

	if(pReadHit->HitLoci.FlagSegs != 0)
		{
		Seg1Hitlen = AdjHitLen(&pReadHit->HitLoci.Hit.Seg[1]);
		Seg1LeftTrimLen = pReadHit->HitLoci.Hit.Seg[1].TrimLeft;
		Seg1RightTrimLen = pReadHit->HitLoci.Hit.Seg[1].TrimRight;
		if(pReadHit->HitLoci.Hit.FlgSplice == 1)  // splice
			{
			MAPQ -= 20;
			GapLen = (int)pReadHit->HitLoci.Hit.Seg[1].MatchLoci - (int)(pReadHit->HitLoci.Hit.Seg[0].MatchLoci + pReadHit->HitLoci.Hit.Seg[0].MatchLen);
			pBAMalign->cigar[CigarIdx++] = GapLen << 4 | 3;		// 'N'
			}
		else   // else if an InDel
			{
			MAPQ -= 10;
			if(pReadHit->HitLoci.Hit.FlgInsert)
				{
				GapLen = pReadHit->ReadLen -
							((pReadHit->HitLoci.Hit.Seg[0].MatchLen - pReadHit->HitLoci.Hit.Seg[0].TrimRight) +
								(pReadHit->HitLoci.Hit.Seg[1].MatchLen  - pReadHit->HitLoci.Hit.Seg[1].TrimLeft));
				if(LimitGapErrs-- > 0 && GapLen <= 0)
					gDiagnostics.DiagOut(eDLInfo,gszProcName,"Check - an apparent insertion GapLen of %d",GapLen);
				pBAMalign->cigar[CigarIdx++] = GapLen << 4 | 1;		// 'I'
				}
			else
				{
				GapLen = (int)pReadHit->HitLoci.Hit.Seg[1].MatchLoci -
					(int)(pReadHit->HitLoci.Hit.Seg[0].MatchLoci + pReadHit->HitLoci.Hit.Seg[0].MatchLen);
				if(LimitGapErrs-- > 0 && GapLen <= 0)
					gDiagnostics.DiagOut(eDLInfo,gszProcName,"Check - an apparent deletion GapLen of %d, Seg[1].MatchLoci:%d, Seg[0].MatchLoci: %d,Seg[0].MatchLen: %d",GapLen,
										pReadHit->HitLoci.Hit.Seg[1].MatchLoci,pReadHit->HitLoci.Hit.Seg[0].MatchLoci,pReadHit->HitLoci.Hit.Seg[0].MatchLen);
				pBAMalign->cigar[CigarIdx++] = abs(GapLen) << 4 | 2;		// 'D'
				}
			}

		if(pReadHit->HitLoci.Hit.Seg[0].Strand == '+')
			{
			if(Seg1LeftTrimLen != 0)
				pBAMalign->cigar[CigarIdx++] = Seg1LeftTrimLen << 4 | 4;		// 'S'
			pBAMalign->cigar[CigarIdx++] = Seg1Hitlen << 4 | 0;		// 'M'
			if(Seg1RightTrimLen > 0)
				pBAMalign->cigar[CigarIdx++] = Seg1RightTrimLen << 4 | 4;		// 'S'
			}
		else
			{
			if(Seg1RightTrimLen != 0)
				pBAMalign->cigar[CigarIdx++] = Seg1RightTrimLen << 4 | 4;		// 'S'
			pBAMalign->cigar[CigarIdx++] = Seg1Hitlen << 4 | 0;		// 'M'
			if(Seg1LeftTrimLen > 0)
				pBAMalign->cigar[CigarIdx++] = Seg1LeftTrimLen << 4 | 4;		// 'S'
			}
		}

	if(ReadIs > 0 && pReadHit->FlgPEAligned && pPEReadHit->FlgPEAligned && pPEReadHit->NAR == eNARAccepted)
		PNext = AdjStartLoci(&pPEReadHit->HitLoci.Hit.Seg[0]);
	else
		PNext = -1;
	MAPQ = max(1,(int)(MAPQ * ((double)(Seg0Hitlen + Seg1Hitlen) / pReadHit->ReadLen)));
	if(MAPQ > 254)
		MAPQ = 254;
	QNameLen = 1 + (int)strlen(pszQName);
	pBAMalign->NumReadNameBytes = QNameLen;
	strcpy(pBAMalign->read_name,pszQName);
	pBAMalign->NumCigarBytes = CigarIdx * sizeof(uint32_t);
	pBAMalign->flag_nc = Flags << 16 | CigarIdx;
	pBAMalign->refID = RefID;
	pBAMalign->next_refID = PNext == -1 ? -1 : RefID;
	pBAMalign->pos = SEStart;
	pBAMalign->end = SEStart+SELen-1;
	pBAMalign->bin_mq_nl = BAMreg2bin(pBAMalign->pos,SEStart+SELen) << 16 | MAPQ << 8 | QNameLen;
	pBAMalign->next_pos = PNext;
	if(TLen >= 0)
		pBAMalign->tlen = TLen;
	else
		pBAMalign->tlen = -1 * abs(TLen);
	pBAMalign->l_seq = pReadHit->ReadLen;
	pBAMalign->NumSeqBytes = (pReadHit->ReadLen + 1)/2;
	pBAMalign->NumAux = 0;
	}
else   // treating as being unaligned
	{
	QNameLen = 1 + (int)strlen(pszQName);
	pBAMalign->NumReadNameBytes = QNameLen;
	strcpy(pBAMalign->read_name,pszQName);
	pBAMalign->NumCigarBytes = 4;
	pBAMalign->cigar[0] = pReadHit->ReadLen << 4 | 0;
	pBAMalign->flag_nc = (Flags << 16) | 0x01;
	pBAMalign->refID = -1;
	pBAMalign->next_refID = -1;
	pBAMalign->pos = -1;
	pBAMalign->end = 0;
	pBAMalign->bin_mq_nl = (128 << 8) | QNameLen;
	pBAMalign->next_pos = -1;
	pBAMalign->tlen = 0;
	pBAMalign->l_seq = pReadHit->ReadLen;
	pBAMalign->NumSeqBytes = (pReadHit->ReadLen + 1)/2;
	pBAMalign->NumAux = 1;
	pBAMalign->auxData[0].tag[0] = 'Y';
	pBAMalign->auxData[0].tag[1] = 'U';
	pBAMalign->auxData[0].val_type = 'Z';
	pBAMalign->auxData[0].NumVals = 1;
	pBAMalign->auxData[0].array_type = 'A';		// not actually required because value type is a string
	strcpy((char *)pBAMalign->auxData[0].value,m_NARdesc[pReadHit->NAR].pszNAR);
	}

if(pReadHit->NAR == eNARAccepted && pReadHit->HitLoci.Hit.Seg[0].Strand != '+')   // 1.1.6 seems that downstream applications are expecting the read sequences to be reverse complemented if sequence was mapped to Crick strand????
	CSeqTrans::ReverseComplement(pReadHit->ReadLen,Sequence);

uint8_t Byte;
int Ofs;

etSeqBase Base;

Byte = 0;
Ofs = 0;
do
	{
	Base = Sequence[Ofs];
	switch(Base) {       // `=ACMGRSVTWYHKDBN' --> [0, 15];
		case eBaseA:
			Byte |= 1;
			break;
		case eBaseC:
			Byte |= 2;
			break;
		case eBaseG:
			Byte |= 4;
			break;
		case eBaseT:
			Byte |= 8;
			break;
		default:
			Byte |= 15;
		}
	if(!(Ofs & 0x01))
		Byte <<= 4;

	if((Ofs & 0x01) || Ofs == pReadHit->ReadLen-1)
		{
		pBAMalign->seq[Ofs/2] = Byte;
		Byte = 0;
		}
	}
while(++Ofs < pReadHit->ReadLen);

return(eBSFSuccess);
}




// ReplaceTabs
char *
CKAligner::ReplaceTabs(char *pszTabTxt) // Inplace replacement of any tabs with a single space char
{
char Chr;
char *pTxt = pszTabTxt;
while(Chr = *pTxt++)
	if(Chr == '\t')
		pTxt[-1]=' ';
return(pszTabTxt);
}

int
CKAligner::AppendStr(char *pszBuff,	// write to this buffer
		  char LeadSep,		// if != '\0' then prefix with this separator (usually ',' or '\t')
		  char Quote,		// if != '\0' then quote string with this char (usually single or double quote char)
		  char *pStr,		// '\0' terminated string
		  char TrailSep)	// if != '\0' then suffix with this separator (usually ',' or '\t' or '\n')
{
int Len = 0;
if(LeadSep != '\0')
	{
	*pszBuff++ = LeadSep;
	Len += 1;
	}
if(Quote != '\0')
	{
	*pszBuff++ = Quote;
	Len += 1;
	}
while(*pStr)
	{
	*pszBuff++ = *pStr++;
	Len += 1;
	}
if(Quote != '\0')
	{
	*pszBuff++ = Quote;
	Len += 1;
	}
if(TrailSep != '\0')
	{
	*pszBuff++ = TrailSep;
	Len += 1;
	}
*pszBuff = '\0';
return(Len);
}

int
CKAligner::AppendChrs(char *pszBuff,	// write to this buffer
		  char LeadSep,		// if != '\0' then prefix with this separator (usually ',' or '\t')
		  char Quote,		// if != '\0' then quote chars with this char (usually single or double quote char)
		  int NumChrs,		// number of chars to append
		  char *Chrs,		// pts to chars to append
		  char TrailSep)	// if != '\0' then suffix with this separator (usually ',' or '\t' or '\n')
{
int Len = 0;
if(LeadSep != '\0')
	{
	*pszBuff++ = LeadSep;
	Len += 1;
	}
if(Quote != '\0')
	{
	*pszBuff++ = Quote;
	Len += 1;
	}
if(NumChrs)
	{
	while(NumChrs--)
		{
		*pszBuff++ = *Chrs++;
		Len += 1;
		}
	}
if(Quote != '\0')
	{
	*pszBuff++ = Quote;
	Len += 1;
	}
if(TrailSep != '\0')
	{
	*pszBuff++ = TrailSep;
	Len += 1;
	}
*pszBuff = '\0';
return(Len);
}


// very fast version of uitoa
int							// length written
CKAligner::AppendUInt(char *pszBuff,	// write to this buffer
		  char LeadSep,		// if > '\0' then prefix with this separator (usually ',' or '\t')
		  uint32_t Value,
		  char TrailSep)	// if > '\0' then suffix with this separator (usually ',' or '\t' or '\n')
{
int Len = 0;
char *pChr;
char *pMark;
char Tmp;
if(LeadSep != '\0')
	{
	*pszBuff++ = LeadSep;
	Len += 1;
	}

if(Value)
	{
	pChr = pszBuff;
	while(Value)
		{
		*pChr++ = '0' + (char)(Value%10);
		Value/=10;
		Len += 1;
		}
	pMark = pChr;
	*pChr-- = '\0';
	while(pszBuff < pChr)
		{
		Tmp = *pChr;
		*pChr-- = *pszBuff;
		*pszBuff++ = Tmp;
		}
	pszBuff = pMark;
	}
else
	{
	Len += 1;
	*pszBuff++ = '0';
	}
if(TrailSep)
	{
	Len += 1;
	*pszBuff++ = TrailSep;
	}
*pszBuff = '\0';
return(Len);
}

// user may be interested in the distribution of the aligner induced substitutions, after any auto-trimming of flanks,
// along the length of the reads and how this distribution relates to the quality scores
int
CKAligner::WriteSubDist(tsReadHit *pReadHit)
{
int QScoreIdx;
int NumMSubs;
int SeqIdx;
uint8_t *pSeqVal;
etSeqBase *pAssembSeq;
tsSegLoci *pSeg;

etSeqBase AssembSeq[cMaxFastQSeqLen+ 10000];	// to hold targeted genome assembly sequence

if(m_hStatsFile == -1 || pReadHit == nullptr || pReadHit->NAR != eNARAccepted || pReadHit->HitLoci.FlagSegs !=0) // if read was segmented because of InDel or splice junction then can't handle, simply slough
	return(eBSFSuccess);

pSeg = &pReadHit->HitLoci.Hit.Seg[0];
if(pSeg->ChromID == 0)
	return(eBSFSuccess);
if(pSeg->Strand == '\0')	// default strand to be sense if not specified
	pSeg->Strand = '+';

m_pSfxArray->GetSeq(pSeg->ChromID,AdjStartLoci(pSeg),AssembSeq,AdjHitLen(pSeg));	// get sequence for entry starting at offset and of length len
if(pSeg->Strand == '-')
	CSeqTrans::ReverseComplement(AdjHitLen(pSeg),AssembSeq);
		
if(m_MaxAlignLen < pReadHit->ReadLen)
	m_MaxAlignLen = pReadHit->ReadLen;

pSeqVal = &pReadHit->Read[pReadHit->DescrLen+1];
pSeqVal += pSeg->ReadOfs + pSeg->TrimLeft;
pAssembSeq = AssembSeq;
NumMSubs = 0;
for(SeqIdx = pSeg->ReadOfs + pSeg->TrimLeft; SeqIdx < (pReadHit->ReadLen - pSeg->TrimRight); SeqIdx++, pSeqVal++, pAssembSeq++)
	{
	QScoreIdx = (*pSeqVal >> 4) & 0x0f;
	// note that Phred scores were scaled to fit within 4bits (((Phred + 2)* 15)/40) by genreads
	if(QScoreIdx <= 3)		//Phred 0..8?
		QScoreIdx = 0;
	else
		if(QScoreIdx <= 7)	// Phred 9..19?
			QScoreIdx = 1;
		else
			if(QScoreIdx <= 11) // Phred 20..29?
				QScoreIdx = 2;
			else
				QScoreIdx = 3;	// Phred 30+
	m_AlignQSubDist[QScoreIdx][SeqIdx].QInsts += 1;
	if((*pSeqVal & 0x07) != (*pAssembSeq & 0x07))
		{
		m_AlignQSubDist[QScoreIdx][SeqIdx].Subs += 1;
		NumMSubs += 1;
		}
	}
if(m_MaxMSubDist < NumMSubs)
	m_MaxMSubDist = NumMSubs;
m_AlignMSubDist[NumMSubs] += 1;

return(eBSFSuccess);
}


int
CKAligner::WriteReadHits(bool bPEProc)		   // true if processing paired ends	
{
const char *pszBsMap;
int LineLen;
char szChromName[128];
int SeqIdx;
etSeqBase ReadSeq[cMaxFastQSeqLen+ 10000];	// to hold sequence (sans quality scores) for current read
uint8_t *pSeqVal;
etSeqBase *pReadSeq;

tsReadHit *pReadHit;
tBSFEntryID PrevTargEntry;

if(m_FMode > eFMbed && (m_hJctOutFile == -1 && m_hIndOutFile == -1))
	return(0);

m_MaxAlignLen = 0;

if(m_FMode == eFMbed)
	{
	LineLen = sprintf(m_pszLineBuff,"track type=bed name=\"%s\" description=\"%s\"\n",m_pszTrackTitle,m_pszTrackTitle);
	if(!m_bgzOutFile)
		{
		if(!CUtility::RetryWrites(m_hOutFile,m_pszLineBuff,LineLen))
			{
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"RetryWrites() error");
			return(eBSFerrWrite);
			}
		}
	else
		CUtility::SafeWrite_gz(m_gzOutFile,m_pszLineBuff,LineLen);
	}

if(m_hJctOutFile != -1)
	{
	LineLen = sprintf(m_pszLineBuff,"track type=bed name=\"JCT_%s\" description=\"%s\"\n",m_pszTrackTitle,m_pszTrackTitle);
	if(!CUtility::RetryWrites(m_hJctOutFile,m_pszLineBuff,LineLen))
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"RetryWrites() error");
		return(eBSFerrWrite);
		}
	}

if(m_hIndOutFile != -1)
	{
	LineLen = sprintf(m_pszLineBuff,"track type=bed name=\"IND_%s\" description=\"%s\"\n",m_pszTrackTitle,m_pszTrackTitle);
	if(!CUtility::RetryWrites(m_hIndOutFile,m_pszLineBuff,LineLen))
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"RetryWrites() error");
		return(eBSFerrWrite);
		}
	}
LineLen = 0;

pReadHit = nullptr;
LineLen = 0;
PrevTargEntry = 0;
const char *pszAlignType;

bool bPrevInDelSeg = false;
bool bPrevJunctSeg = false;
bool bPrevAlignSeg = false;
bool bSkipBEDformat = false;;

while((pReadHit = IterSortedReads(pReadHit))!=nullptr)
	{
	bSkipBEDformat = false;
	if(pReadHit->NAR == eNARAccepted)
		{
		if(pReadHit->HitLoci.Hit.FlgInDel)
			pszAlignType = "ari";
		else
			if(pReadHit->HitLoci.Hit.FlgSplice)
				pszAlignType = "arj";
			else
				{
				pszAlignType = "ar";
				if(m_FMode > eFMbed)
					bSkipBEDformat = true;
				}

		if(pReadHit->HitLoci.Hit.FlgInDel)
			pszAlignType = "ari";
		else
			if(pReadHit->HitLoci.Hit.FlgSplice)
				pszAlignType = "arj";
			else
				{
				pszAlignType = "ar";
				if(m_FMode > eFMbed)
					bSkipBEDformat = true;
				}

		tsSegLoci *pSeg;
		int SegIdx;

		pSeg = &pReadHit->HitLoci.Hit.Seg[0];
		if(pSeg->ChromID != (uint32_t)PrevTargEntry)
			{
			m_pSfxArray->GetIdentName(pSeg->ChromID,sizeof(szChromName),szChromName);
			PrevTargEntry = pSeg->ChromID;
			}

		if(pSeg->Strand == '\0')	// default strand to be sense if not specified
			pSeg->Strand = '+';

		int Score = (int)min(1000.0,(999 * m_OctSitePrefs[pSeg->Strand == '+' ? 0 : 1][pReadHit->SiteIdx].RelScale));

		if(m_FMode >= eFMbed)
			{
			if(pReadHit->HitLoci.FlagSegs==0)
				{
				if(bPrevInDelSeg || bPrevJunctSeg)
					{
					if(LineLen > 0)
						{
						if(bPrevInDelSeg && m_hIndOutFile != -1)
							{
							if(!CUtility::RetryWrites(m_hIndOutFile,m_pszLineBuff,LineLen))
								{
								gDiagnostics.DiagOut(eDLInfo,gszProcName,"RetryWrites() error");
								return(eBSFerrWrite);
								}
							}
						else
							if(bPrevJunctSeg && m_hJctOutFile != -1)
								{
								if(!CUtility::RetryWrites(m_hJctOutFile,m_pszLineBuff,LineLen))
									{
									gDiagnostics.DiagOut(eDLInfo,gszProcName,"RetryWrites() error");
									return(eBSFerrWrite);
									}
								}
						LineLen = 0;
						}
					bPrevInDelSeg = false;
					bPrevJunctSeg = false;
					}
				pSeg = &pReadHit->HitLoci.Hit.Seg[0];

				if(!bSkipBEDformat)
					LineLen+=sprintf(&m_pszLineBuff[LineLen],"%s\t%d\t%d\t%s\t%d\t%c\n",
								szChromName,AdjStartLoci(pSeg),AdjEndLoci(pSeg) + 1,pszAlignType,Score,pSeg->Strand);
				bPrevAlignSeg = true;
				}
			else // segmented alignment
				{
				if(bPrevAlignSeg)
					{
					if(LineLen > 0)
						{
						if(m_FMode == eFMbed)
							{
							if(!m_bgzOutFile)
								CUtility::RetryWrites(m_hOutFile,m_pszLineBuff,LineLen);
							else
								CUtility::SafeWrite_gz(m_gzOutFile,m_pszLineBuff,LineLen);
							}
						LineLen = 0;
						}
					bPrevAlignSeg = false;
					}
				if(pReadHit->HitLoci.Hit.FlgInDel)
					{
					if(bPrevJunctSeg && m_hJctOutFile != -1)
						{
						if(LineLen > 0)
							if(!CUtility::RetryWrites(m_hJctOutFile,m_pszLineBuff,LineLen))
								{
								gDiagnostics.DiagOut(eDLInfo,gszProcName,"RetryWrites() error");
								return(eBSFerrWrite);
								}
						LineLen = 0;
						}
					bPrevJunctSeg = false;
					bPrevInDelSeg = true;
					}
				else
					{
					if(bPrevInDelSeg && m_hIndOutFile != -1)
						{
						if(LineLen > 0)
							{
							if(!CUtility::RetryWrites(m_hIndOutFile,m_pszLineBuff,LineLen))
								{
								gDiagnostics.DiagOut(eDLInfo,gszProcName,"RetryWrites() error");
								return(eBSFerrWrite);
								}
							LineLen = 0;
							}
						}
					bPrevInDelSeg = false;
					bPrevJunctSeg = true;
					}

				int AjAlignStartLoci;
				int AjAlignEndLoci;

				if(!bSkipBEDformat)
					{
					AjAlignStartLoci = AdjAlignStartLoci(&pReadHit->HitLoci.Hit);
					AjAlignEndLoci = AdjAlignEndLoci(&pReadHit->HitLoci.Hit);
					LineLen+=sprintf(&m_pszLineBuff[LineLen],"%s\t%d\t%d\t%s\t%d\t%c\t%d\t%d\t0\t2\t%d,%d\t0,%d\n",
									szChromName,AjAlignStartLoci,AjAlignEndLoci+1,pszAlignType,Score,pSeg->Strand,AjAlignStartLoci,AjAlignEndLoci+1,
									AdjHitLen(pSeg),AdjHitLen(&pSeg[1]),AdjStartLoci(&pSeg[1])-AdjStartLoci(pSeg));
					}
				}

			if((cAllocLineBuffSize - LineLen) < 1000)
				{
				if(bPrevAlignSeg)
					{
					if(m_FMode == eFMbed)
						{
						if(!m_bgzOutFile)
							CUtility::RetryWrites(m_hOutFile,m_pszLineBuff,LineLen);
						else
							CUtility::SafeWrite_gz(m_gzOutFile,m_pszLineBuff,LineLen);
						}
					LineLen = 0;
					bPrevAlignSeg = false;
					}
				else
					{
					if(bPrevInDelSeg && m_hIndOutFile != -1)
						CUtility::RetryWrites(m_hIndOutFile,m_pszLineBuff,LineLen);
					else
						if(bPrevJunctSeg && m_hJctOutFile != -1)
							CUtility::RetryWrites(m_hJctOutFile,m_pszLineBuff,LineLen);
					LineLen = 0;
					bPrevInDelSeg = false;
					bPrevJunctSeg = false;
					}
				}
			continue;
			}

		if(!m_bIsSOLiD)
			{
			pSeqVal = &pReadHit->Read[pReadHit->DescrLen+1];
			pReadSeq = ReadSeq;
			for(SeqIdx = 0; SeqIdx < pReadHit->ReadLen; SeqIdx++,pReadSeq++,pSeqVal++)
				*pReadSeq = (*pSeqVal & 0x07);
			}
		bPrevAlignSeg = true;
		for(SegIdx = 0; SegIdx < 2; SegIdx++)
			{
			pSeg = &pReadHit->HitLoci.Hit.Seg[SegIdx];
			if(pSeg->ChromID == 0)
				continue;
			if(pSeg->Strand == '\0')	// default strand to be sense if not specified
				pSeg->Strand = '+';
			if(pSeg->ChromID != (uint32_t)PrevTargEntry)
				{
				m_pSfxArray->GetIdentName(pSeg->ChromID,sizeof(szChromName),szChromName);
				PrevTargEntry = pSeg->ChromID;
				}

			LineLen += AppendUInt(&m_pszLineBuff[LineLen],0,pReadHit->ReadID,',');
			LineLen += AppendStr(&m_pszLineBuff[LineLen],0,'"',(char *)pszAlignType,',');
			LineLen += AppendStr(&m_pszLineBuff[LineLen],0,'"',m_szTargSpecies,',');
			LineLen += AppendStr(&m_pszLineBuff[LineLen],0,'"',szChromName,',');

			LineLen += AppendUInt(&m_pszLineBuff[LineLen],0,AdjStartLoci(pSeg),',');
			LineLen += AppendUInt(&m_pszLineBuff[LineLen],0,AdjEndLoci(pSeg),',');
			LineLen += AppendUInt(&m_pszLineBuff[LineLen],0,AdjHitLen(pSeg),',');
			LineLen += AppendChrs(&m_pszLineBuff[LineLen],0,'"',1,(char *)&pSeg->Strand,',');
			LineLen += AppendUInt(&m_pszLineBuff[LineLen],0,Score,',');
			LineLen += AppendUInt(&m_pszLineBuff[LineLen],0,0,'\0');

			if(m_bBisulfite) {
				switch(pReadHit->HitLoci.Hit.BisBase) {
					case eBaseT:
						pszBsMap = "TC:C";
						break;
					case eBaseA:
						pszBsMap = "AG:T";
						break;
					default:
						pszBsMap = "?:?";
					}
				}
			else
				pszBsMap = "N/A";

			LineLen += AppendUInt(&m_pszLineBuff[LineLen],',',pReadHit->NumReads,',');
			LineLen += AppendUInt(&m_pszLineBuff[LineLen],0,pSeg->TrimMismatches,',');
			LineLen += AppendStr(&m_pszLineBuff[LineLen],0,'"',(char *)pszBsMap,',');
			LineLen += AppendStr(&m_pszLineBuff[LineLen],0,'"',(char *)pReadHit->Read,'\0');

			LineLen += AppendStr(&m_pszLineBuff[LineLen],',','"',CSeqTrans::MapSeq2Ascii(&ReadSeq[pSeg->ReadOfs+pSeg->TrimLeft],AdjHitLen(pSeg)),0);

			LineLen += AppendStr(&m_pszLineBuff[LineLen],0,0,(char *)"\n",0);
			if(LineLen + ((cMaxFastQSeqLen * 2) + 1024) > cAllocLineBuffSize)
				{
				if(!m_bgzOutFile)
					CUtility::RetryWrites(m_hOutFile,m_pszLineBuff,LineLen);
				else
					CUtility::SafeWrite_gz(m_gzOutFile,m_pszLineBuff,LineLen);
				LineLen = 0;
				}
			}

		// user may be interested in the distribution of the aligner induced substitutions, after any auto-trimming of flanks,
		// along the length of the reads and how this distribution relates to the quality scores
		if(m_hStatsFile != -1 && m_FMode <= eFMbed)
			WriteSubDist(pReadHit);
		}
	}
if(LineLen)
	{
	if(bPrevAlignSeg && m_FMode <= eFMbed)
		{
		if(!m_bgzOutFile)
			CUtility::RetryWrites(m_hOutFile,m_pszLineBuff,LineLen);
		else
			CUtility::SafeWrite_gz(m_gzOutFile,m_pszLineBuff,LineLen);

		LineLen = 0;
		bPrevAlignSeg = false;
		}
	else
		{
		if(bPrevInDelSeg && m_hIndOutFile != -1)
			CUtility::RetryWrites(m_hIndOutFile,m_pszLineBuff,LineLen);
		else
			if(bPrevJunctSeg && m_hJctOutFile != -1)
				CUtility::RetryWrites(m_hJctOutFile,m_pszLineBuff,LineLen);
		LineLen = 0;
		bPrevInDelSeg = false;
		bPrevJunctSeg = false;
		}
	}
return(eBSFSuccess);
}


int
CKAligner::AddMultiHit(tsReadHit *pReadHit)
{
int NumMultiHits;
int CopyLen;
size_t HitLen;
tsReadHit *pMultiHit;

HitLen = sizeof(tsReadHit) + pReadHit->ReadLen + pReadHit->DescrLen; 
#ifdef _WIN32
DWORD WaitRslt = WaitForSingleObject(m_hMtxMultiMatches,INFINITE);
#else
pthread_mutex_lock(&m_hMtxMultiMatches);
#endif
if(m_NxtMultiAllOfs + (2 * HitLen) >= m_AllocMultiAllMem)
	{
	size_t memreq = m_AllocMultiAllMem + (cAllocMultihits * HitLen);
#ifdef _WIN32
	pMultiHit = (tsReadHit *) realloc(m_pMultiAll,memreq);
#else
	pMultiHit = (tsReadHit *)mremap(m_pMultiAll,m_AllocMultiAllMem,memreq,MREMAP_MAYMOVE);
	if(pMultiHit == MAP_FAILED)
		pMultiHit = nullptr;
#endif
	if(pMultiHit == nullptr)
		{
#ifdef _WIN32
ReleaseMutex(m_hMtxMultiMatches);
#else
pthread_mutex_unlock(&m_hMtxMultiMatches);
#endif
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddMultiHit: Memory re-allocation to %zd bytes - %s",memreq,strerror(errno));
		return(eBSFerrMem);
		}
	m_pMultiAll = pMultiHit;
	m_AllocMultiAllMem = memreq;
	}
pMultiHit = (tsReadHit *)((uint8_t *)m_pMultiAll + m_NxtMultiAllOfs);

CopyLen = sizeof(tsReadHit) + pReadHit->DescrLen + pReadHit->ReadLen;
memcpy(pMultiHit,pReadHit,CopyLen);

m_NxtMultiAllOfs += CopyLen;
m_NumMultiAll += 1;
pMultiHit->ReadID = m_NumMultiAll;
NumMultiHits = m_NumMultiAll;
#ifdef _WIN32
ReleaseMutex(m_hMtxMultiMatches);
#else
pthread_mutex_unlock(&m_hMtxMultiMatches);
#endif
return((int)NumMultiHits);
}

int				// normally NumHits, but will be actual number of hits if unable to accept any of the loci hit because of chromosome filtering
CKAligner::WriteHitLoci(tsThreadMatchPars *pThreadPars,tsReadHit *pReadHit,int NumHits,tsHitLoci *pHits)
{
int Rslt;
tsHitLoci *pHit;
tsSegLoci *pSeg;

int ReadHitBuffIdx;						// index into output szReadHits

tBSFEntryID PrevTargEntry;
uint8_t ReadHit[sizeof(tsReadHit) + cMaxKADescrLen + cMaxFastQSeqLen + 10000];
tsReadHit *pMultiHit;

m_MaxAlignLen = 0;
PrevTargEntry = 0;
ReadHitBuffIdx = 0;

// check if hits are to chromosomes which are to be retained
if(NumHits > 0 && m_RegExprs.HasRegExprs())
	{
	tsHitLoci *pAcceptHit;
	int AcceptHits;
	AcceptHits = 0;
	pHit = pHits;
	pAcceptHit = pHit;
	for(int HitIdx = 0; HitIdx < NumHits; HitIdx++,pHit++)
		{
		pSeg = &pHit->Seg[0];
		if(AcceptThisChromID(pSeg->ChromID))
			{
			AcceptHits += 1;
			if(pAcceptHit != pHit)
				*pAcceptHit++ = *pHit;
			}
		}
	NumHits = AcceptHits;
	if(m_FMode != eFMsamAll && NumHits == 0)
		return(0);
	}

PrevTargEntry = 0;
pHit = pHits;

pMultiHit = (tsReadHit *)&ReadHit;
memmove(&ReadHit,pReadHit,sizeof(tsReadHit) + pReadHit->ReadLen + pReadHit->DescrLen);
if(NumHits == 0)
	{
	pMultiHit->NumHits = 0;
	NumHits = 1;
	}
else
	pMultiHit->NumHits = 1;


for(int HitIdx = 0; HitIdx < NumHits; HitIdx++,pHit++)
	{
	if(pMultiHit->NumHits > 0)
		{
		memmove(&pMultiHit->HitLoci.Hit,pHit,sizeof(tsHitLoci));
		pMultiHit->HitLoci.FlagSegs = (pHit->FlgInDel == 1 || pHit->FlgSplice == 1) ? 1 : 0;
		}
	else
		pMultiHit->HitLoci.FlagSegs = 0;

	if((Rslt = AddMultiHit(pMultiHit)) < eBSFSuccess)
		return(Rslt);
	}

return(NumHits);
}

void
CKAligner::InitialiseWIGSpan(void) // initialise WIG span vars to values corresponding to no spans having been previously reported
{
m_WIGChromID = 0;
m_WIGRptdChromID = 0;
m_WIGSpanLoci = 0;
m_WIGSpanLen = 0;
m_WIGRptdSpanLen = 0;
m_WIGSpanCnts = 0;
}

int
CKAligner::CompleteWIGSpan(bool bWrite)				// close off any current WIG span ready to start any subsequent span
{
char szChromName[cMaxDatasetSpeciesChrom + 1];

	// if existing span then write that span out
if(m_WIGChromID != 0 && m_WIGSpanLen > 0 && m_WIGSpanLoci > 0 && m_WIGSpanCnts > 0)
	{
	// has chrom and/or span changed since previously writing out a span?
	if(m_WIGChromID != m_WIGRptdChromID || m_WIGSpanLen != m_WIGRptdSpanLen)
		{
		m_pSfxArray->GetIdentName(m_WIGChromID, sizeof(szChromName), szChromName);
		m_CovSegBuffIdx += sprintf(&m_pszCovSegBuff[m_CovSegBuffIdx],"variableStep chrom=%s span=%d\n",szChromName,m_WIGSpanLen);
		m_WIGRptdChromID = m_WIGChromID;
		m_WIGRptdSpanLen = m_WIGSpanLen;
		}
	m_CovSegBuffIdx += sprintf(&m_pszCovSegBuff[m_CovSegBuffIdx],"%d %d\n",m_WIGSpanLoci,(uint32_t)((m_WIGSpanCnts + m_WIGSpanLen-1)/m_WIGSpanLen));
	}
if((bWrite && m_CovSegBuffIdx) || (m_CovSegBuffIdx + 500) >  m_AllocdCovSegBuff)
	{
	if(!CUtility::RetryWrites(m_hWIGSpansFile, m_pszCovSegBuff, m_CovSegBuffIdx))
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"RetryWrites() error");
		return(eBSFerrWrite);
		}
	m_CovSegBuffIdx=0;
	}
m_WIGSpanLoci = 0;
m_WIGSpanLen = 0;
m_WIGSpanCnts = 0;
return(eBSFSuccess);
}



int
CKAligner::AccumWIGCnts(uint32_t ChromID,	// accumulate wiggle counts into variableStep spans over this chromosome
			uint32_t Loci,		// this loci - starts from 1 not 0!
			uint32_t Cnts,		 // has this many counts attributed
			uint32_t MaxSpanLen) // allow WIG spans to be this maximal length
{
	int Rslt;
	uint32_t Meanx100;
	if(ChromID != m_WIGChromID || m_WIGSpanLen >= MaxSpanLen || Cnts == 0)		// onto a different chromosome, or current span is at maximal length?
	{
		if(m_WIGChromID != 0)
		{
			if((Rslt = CompleteWIGSpan()) < 0)
				return(Rslt);
		}
		if(Cnts > 0)
		{
			m_WIGChromID = ChromID;
			m_WIGSpanLoci = Loci;
			m_WIGSpanLen = 1;
			m_WIGSpanCnts = Cnts;
		}
		return(eBSFSuccess);
	}

	if(m_WIGSpanLen == 0 || m_WIGSpanCnts == 0)
	{
		m_WIGSpanLoci = Loci;
		m_WIGSpanLen = 1;
		m_WIGSpanCnts = (uint64_t)Cnts;
		return(eBSFSuccess);
	}

	Meanx100 = 100 * (uint32_t)(m_WIGSpanCnts / (uint64_t)m_WIGSpanLen);
	if((Cnts <= 5 && (Cnts * 100) != Meanx100) || (Meanx100 < (Cnts * 75) || Meanx100 >= (Cnts * 125)))
	{
		// write current span out
		if((Rslt = CompleteWIGSpan()) < 0)
			return(Rslt);
		m_WIGSpanLoci = Loci;
		m_WIGSpanLen = 1;
		m_WIGSpanCnts = Cnts;
		return(eBSFSuccess);
	}
	m_WIGSpanCnts += Cnts;
	m_WIGSpanLen = Loci - m_WIGSpanLoci + 1;
	return(eBSFSuccess);
}


// OutputSNPs
// Currently can't process for SNPs in InDels or splice junctions
// FDR: Benjamini�Hochberg
// QValue == acceptable FDR e.g. 0.05% or 0.01%
// PValue == 1.0 - Stats.Binomial(TotBasesInColumn,NumBasesInColMismatching+1,GlobalSeqErrRate);
// PValueIdx == sorted index of PValue and locus pairs, 1..k
// Generate PValues for all alignment columns meeting minimum constraints into an array of structures containing column loci and associated PValues
// Sort array of structures ascending on PValues
// Iterate array 1 to k and accept as SNPs those elements with PValues < (PValueIdx/k) * QValue
int
CKAligner::OutputSNPs(void)
{
int Rslt;
double PValue;
double GlobalSeqErrRate;
double LocalSeqErrRate;
tsSNPcnts* pSNP;
uint32_t Loci;

uint32_t PrevScaledLociCoverage;
uint32_t	CoverageSegLen;
uint32_t	StartCoverageSegLoci;

int Idx;
int NumSNPs;
int TotBases;
char szChromName[cMaxDatasetSpeciesChrom + 1];
int LineLen;
double Proportion;
double AdjPValue;
int RelRank;
tsLociPValues* pLociPValues;
size_t memreq;
tsSNPcnts* pSNPWinL;
tsSNPcnts* pSNPWinR;
uint32_t LocalBkgndRateWindow;
uint32_t LocalBkgndRateWinFlank;
uint32_t LocalTotMismatches;
uint32_t LocalTotMatches;
uint32_t LocTMM;
uint32_t LocTM;
CStats Stats;
uint8_t *pPackedBaseAlleles;
uint8_t PackedBaseAlleles;

tsMonoSNP sMonoSNP;
tsDiSNP sDiSNP;
tsTriSNP sTriSNP;

tsSNPcnts PrevSNPs[3]; // to hold counts for prev 2 plus currrent SNP - used when reporting Di/TriSNPs

int CurDiSNPLoci;
int PrevDiSNPLoci;
int DiSNPBuffIdx;
char szDiSNPs[4000];

int CurTriSNPLoci;
int PrevTriSNPLoci;
int FirstTriSNPLoci;
int TriSNPBuffIdx;
char szTriSNPs[4000];
int TotNumDiSNPs;
int TotNumTriSNPs;

uint8_t SNPFlanks[9];
uint8_t* pSNPFlank;
int SNPFlankIdx;
int SNPCentroidIdx;
uint8_t Base;
tsSNPCentroid* pCentroid;
int NumCovReads;

m_pSfxArray->GetIdentName(m_pChromSNPs->ChromID, sizeof(szChromName), szChromName);
if(m_pszCovSegBuff == nullptr)
	{
	m_pszCovSegBuff = new char [cAllocCovSegSize];
	m_AllocdCovSegBuff = cAllocCovSegSize;
	}
m_CovSegBuffIdx = 0;

if(!m_bPackedBaseAlleles)			
	{
	if (m_pLociPValues == nullptr)					// will be nullptr first time in
		{
		memreq = cAllocLociPValues * sizeof(tsLociPValues);
#ifdef _WIN32
		m_pLociPValues = (tsLociPValues*)malloc((size_t)memreq);
		if (m_pLociPValues == nullptr)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "OutputSNPs: Memory allocation of %zd bytes failed", (int64_t)memreq);
			return(eBSFerrMem);
			}
#else
		// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
		m_pLociPValues = (tsLociPValues*)mmap(nullptr, (size_t)memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
		if (m_pLociPValues == MAP_FAILED)
			{
			m_pLociPValues = nullptr;
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "OutputSNPs: Memory allocation of %zd bytes through mmap()  failed", (int64_t)memreq, strerror(errno));
			return(eBSFerrMem);
			}
#endif
		m_AllocLociPValuesMem = memreq;
		m_NumLociPValues = 0;
		}
	}
else
	{ // this code block contains the PBA processing
// NOTE: initial allocation size is likely to be sufficent for any single chromosome
// but the reallocation to larger size processing is retained in case there is, in future, a need to realloc! 
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenPBAChrom: Processing chrom '%s'", szChromName);
	if (m_pPackedBaseAlleles == nullptr)					// will be nullptr first time in
		{
		memreq = max((size_t)0x07ffffff0,((size_t)m_pChromSNPs->ChromLen + 0x0fff) * 2);  // allocating a little extra for the header
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenPBAChrom: Initial memory allocation is for %zd bytes", memreq);
#ifdef _WIN32
		m_pPackedBaseAlleles = (uint8_t*)malloc(memreq);
		if (m_pPackedBaseAlleles == nullptr)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "OutputSNPs: Memory allocation of %zd bytes failed", (int64_t)memreq);
			return(eBSFerrMem);
			}
#else
		m_pPackedBaseAlleles = (uint8_t*)mmap(nullptr, memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
		if (m_pPackedBaseAlleles == MAP_FAILED)
			{
			m_pPackedBaseAlleles = nullptr;
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "OutputSNPs: Initial memory allocation of %zd bytes through mmap()  failed", (int64_t)memreq, strerror(errno));
			return(eBSFerrMem);
			}
#endif
		m_AllocPackedBaseAllelesMem = memreq;

		// as this the first chromosome being called then generate header
		// This header contains a series of '\n' separated tagname:values 
		// Following the header are a variable number of chromosomes
		m_NumPackedBaseAlleles = sprintf((char *)m_pPackedBaseAlleles,"Type:%s\nVersion:1\nExperimentID:%s\nReferenceID:%s\nReadsetID:%s","PbA",m_szExperimentName,m_szTargSpecies,m_pszTrackTitle);
		m_NumPackedBaseAlleles+=1;	
		}
	else
		{	// initial allocation and header not required
			// but need to allocate additional memory for this new chrom?
		if (((size_t)m_pChromSNPs->ChromLen + 0x0ffff) > m_AllocPackedBaseAllelesMem) 
			{
			memreq = ((size_t)m_pChromSNPs->ChromLen + 0x0ffff);  // reallocating a little extra to minimize chances of a realloc required if a subsequent larger chrom processed
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "GenPBAChrom: Increasing memory allocation to %zd bytes", memreq);
#ifdef _WIN32
			pPackedBaseAlleles = (uint8_t*)realloc(m_pPackedBaseAlleles, memreq);
			if (pPackedBaseAlleles == nullptr)
				{
#else
				pPackedBaseAlleles = (uint8_t*)mremap(m_pPackedBaseAlleles, m_AllocPackedBaseAllelesMem, memreq, MREMAP_MAYMOVE);
				if (pPackedBaseAlleles == MAP_FAILED)
				{
#endif
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "OutputSNPs: Memory reallocation PackedBaseAlleles to %zd bytes from %zd failed - (%d) %s", memreq, m_AllocPackedBaseAllelesMem, errno, strerror(errno));
				return(eBSFerrMem);
				}
			m_pPackedBaseAlleles = pPackedBaseAlleles;
			m_AllocPackedBaseAllelesMem = memreq;
			}
		m_NumPackedBaseAlleles = 0;	// no header
		}
	
	pPackedBaseAlleles = &m_pPackedBaseAlleles[m_NumPackedBaseAlleles];

	*pPackedBaseAlleles = (uint8_t)strlen(szChromName);
	strcpy((char *)&pPackedBaseAlleles[1],szChromName);
	pPackedBaseAlleles += 2 + *pPackedBaseAlleles;
	*(uint32_t *)pPackedBaseAlleles = m_pChromSNPs->ChromLen;
	pPackedBaseAlleles += 4;
	m_NumPackedBaseAlleles = (uint32_t)(pPackedBaseAlleles - m_pPackedBaseAlleles);
	pSNP = &m_pChromSNPs->Cnts[0];
	InitialiseWIGSpan();
	for (Loci = 0; Loci < m_pChromSNPs->ChromLen; Loci++, pSNP++,pPackedBaseAlleles++,m_NumPackedBaseAlleles++)
		{
		// iterate over each base and score according to proportion of total at this loci
		PackedBaseAlleles = 0;			// assume no coverage
		uint32_t Coverage = pSNP->NumNonRefBases + pSNP->NumRefBases - pSNP->NonRefBaseCnts[4];	// coverage excludes indeterminate bases - 'N'
		if((Rslt=AccumWIGCnts(m_pChromSNPs->ChromID,Loci+1,Coverage))!=eBSFSuccess) // WIG loci start from 1
			return(Rslt);
		double AlleleProp;
		if(Coverage > 0)
			{
			// PBA: Allele A in bits 7.6, C in bits 5.4, G in bits 3.2, T in bits 1.0
			for(int BaseIdx = 0; BaseIdx < 4; BaseIdx++)
				{
				PackedBaseAlleles <<= 2;
				if(BaseIdx == pSNP->RefBase)
					AlleleProp = pSNP->NumRefBases / (double)Coverage;
				else
					AlleleProp = pSNP->NonRefBaseCnts[BaseIdx] / (double)Coverage;
				if(Coverage >= 5)
					{
					if(AlleleProp >= cScorePBA3MinProp)
						PackedBaseAlleles |= 0x03;
					else
						if(AlleleProp >= cScorePBA2MinProp)
							PackedBaseAlleles |= 0x02;
						else
							if(AlleleProp >= cScorePBA1MinProp)
								PackedBaseAlleles |= 0x01;
					}
				else
					{
					if(AlleleProp >= cScorePBA2MinLCProp)
						PackedBaseAlleles |= 0x02;
					else
						if(AlleleProp >= cScorePBA1MinLCProp)
							PackedBaseAlleles |= 0x01;
					}
				}
			}
		*pPackedBaseAlleles = PackedBaseAlleles;
		}
	CompleteWIGSpan(true);
	if(m_hPackedBaseAllelesFile != -1)
		{
		if(m_NumPackedBaseAlleles > 0)
			{
			if(!CUtility::RetryWrites(m_hPackedBaseAllelesFile, m_pPackedBaseAlleles, m_NumPackedBaseAlleles))
				{
				gDiagnostics.DiagOut(eDLInfo,gszProcName,"RetryWrites() error");
				return(eBSFerrWrite);
				}
			m_NumPackedBaseAlleles = 0;
			}
		}
	return(eBSFSuccess);
	}	// completed PBA processing

	// NOTE: set a floor on the global (whole chromosome) sequencing error rate
	GlobalSeqErrRate = max(cMinSeqErrRate, (double)m_pChromSNPs->TotMismatch / (double)(1 + m_pChromSNPs->TotMatch + m_pChromSNPs->TotMismatch));

	pSNPWinR = &m_pChromSNPs->Cnts[0];
	LocalBkgndRateWinFlank = cSNPBkgndRateWindow / 2;
	LocalBkgndRateWindow = (LocalBkgndRateWinFlank * 2) + 1;
	LocalTotMismatches = 0;
	LocalTotMatches = 0;
	CoverageSegLen = 0;
	StartCoverageSegLoci = 0;

	for (Loci = 0; Loci < min(LocalBkgndRateWindow, m_pChromSNPs->ChromLen); Loci++, pSNPWinR++)
		{
		LocalTotMismatches += pSNPWinR->NumNonRefBases;
		LocalTotMatches += pSNPWinR->NumRefBases;
		}

	pLociPValues = m_pLociPValues;
	m_NumLociPValues = 0;
	pSNP = &m_pChromSNPs->Cnts[0];

	pSNPWinL = pSNP;
	LineLen = 0;

	PrevScaledLociCoverage = 0;
	CoverageSegLen = 0;
	StartCoverageSegLoci = 0;
	m_MaxDiSNPSep = min(cDfltMaxDiSNPSep, (int32_t)m_pChromSNPs->MeanReadLen);
	InitialiseWIGSpan();
	for (Loci = 0; Loci < m_pChromSNPs->ChromLen; Loci++, pSNP++)
		{
		// determine background expected error rate from window surrounding the current loci
		if (Loci > LocalBkgndRateWinFlank && (Loci + LocalBkgndRateWinFlank) < m_pChromSNPs->ChromLen)
			{
			// need to ensure that LocalTotMismatches and LocalTotMismatches will never underflow
			if (LocalTotMismatches >= pSNPWinL->NumNonRefBases)
				LocalTotMismatches -= pSNPWinL->NumNonRefBases;
			else
				LocalTotMismatches = 0;
			if (LocalTotMatches >= pSNPWinL->NumRefBases)
				LocalTotMatches -= pSNPWinL->NumRefBases;
			else
				LocalTotMatches = 0;
			LocalTotMismatches += pSNPWinR->NumNonRefBases;
			LocalTotMatches += pSNPWinR->NumRefBases;
			pSNPWinL += 1;
			pSNPWinR += 1;
			}

		TotBases = pSNP->NumNonRefBases + pSNP->NumRefBases;
		if (TotBases > 0)
			{
			m_LociBasesCovered += 1;
			m_LociBasesCoverage += TotBases;
			}

		AccumWIGCnts(m_pChromSNPs->ChromID,Loci,TotBases);

		if (TotBases < m_MinSNPreads) // if not meeting threshold for SNP calling coverage at current loci then iterate onto next loci after recording loci base classification
			continue;

		if (m_hSNPCentsfile != -1)
			{
			// get 4bases up/dn stream from loci with SNP and use these to inc centroid counts of from/to counts
			if (Loci >= cSNPCentfFlankLen && Loci < (m_pChromSNPs->ChromLen - cSNPCentfFlankLen))
				{
				m_pSfxArray->GetSeq(m_pChromSNPs->ChromID, Loci - (uint32_t)cSNPCentfFlankLen, SNPFlanks, cSNPCentroidLen);
				pSNPFlank = &SNPFlanks[cSNPCentroidLen - 1];
				SNPCentroidIdx = 0;
				for (SNPFlankIdx = 0; SNPFlankIdx < cSNPCentroidLen; SNPFlankIdx++, pSNPFlank--)
					{
					Base = *pSNPFlank & 0x07;
					if (Base > eBaseT)
						break;
					SNPCentroidIdx |= (Base << (SNPFlankIdx * 2));
					}
				if (SNPFlankIdx == cSNPCentroidLen)
					m_pSNPCentroids[SNPCentroidIdx].NumInsts += 1;
				}
			}


		if (pSNP->NumNonRefBases < cMinSNPreads)
			continue;
		Proportion = (double)pSNP->NumNonRefBases / TotBases;
		if (Proportion < m_SNPNonRefPcnt)	// needs to be at least this proportion of non-ref bases to be worth exploring as being SNP
			continue;

		// needing to allocate more memory? NOTE: allowing small safety margin of 10 tsLociPValues
		if (m_AllocLociPValuesMem < (sizeof(tsLociPValues) * (m_NumLociPValues + 10)))
			{
			size_t memreq = m_AllocLociPValuesMem + (cAllocLociPValues * sizeof(tsLociPValues));
#ifdef _WIN32
			pLociPValues = (tsLociPValues*)realloc(m_pLociPValues, memreq);
			if (pLociPValues == nullptr)
				{
#else
				pLociPValues = (tsLociPValues*)mremap(m_pLociPValues, m_AllocLociPValuesMem, memreq, MREMAP_MAYMOVE);
				if (pLociPValues == MAP_FAILED)
				{
#endif
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "OutputSNPs: Memory reallocation to %zd bytes failed - %s", memreq, strerror(errno));
				return(eBSFerrMem);
				}
			m_pLociPValues = pLociPValues;
			m_AllocLociPValuesMem = memreq;
			pLociPValues = &m_pLociPValues[m_NumLociPValues];
			}



		if (pSNP->NumNonRefBases <= LocalTotMismatches)
			LocTMM = LocalTotMismatches - pSNP->NumNonRefBases;
		else
			LocTMM = 0;

		if (pSNP->NumRefBases < LocalTotMatches)
			LocTM = LocalTotMatches - pSNP->NumRefBases;
		else
			LocTM = 0;


		if ((LocTMM + LocTM) == 0)
			LocalSeqErrRate = GlobalSeqErrRate;
		else
			{
			LocalSeqErrRate = (double)LocTMM / (double)(LocTMM + LocTM);
			if (LocalSeqErrRate < GlobalSeqErrRate)
				LocalSeqErrRate = GlobalSeqErrRate;
			}
		if (LocalSeqErrRate > cMaxBkgdNoiseThres)	// don't bother attempting to call if the background is too noisy
			continue;
		
		// accepting as being a putative SNP
		if(m_bXCSVFrameShifts)
			{
			// determine frame shifted codon counts for aligned sequences which overlay the SNP site
			memset(pLociPValues->FrameShiftedCodons,0,sizeof(pLociPValues->FrameShiftedCodons));

			NumCovReads = FrameShiftCodons(m_pChromSNPs->ChromID,	// SNP loci is on this chromosome
											Loci,					// returned codons are frame shifted relative to this loci
											&pLociPValues->FrameShiftedCodons[0][0]);	// where to return cnts of frame shifted codons (set of codon counts [3][64])

			// determine the canonical target reference codon in each frame shift
			int RelLoci;
			int FrameShift;
			int FrameCodonIdx;
			etSeqBase LociBase;
			for(FrameShift = 0; FrameShift < 3; FrameShift++)
				{
				pLociPValues->FrameShiftedRefCodons[FrameShift] = -1;
				FrameCodonIdx = 0;
				for(RelLoci = -2; RelLoci <= 0; RelLoci++)
					{
					if((LociBase = (m_pSfxArray->GetBase(m_pChromSNPs->ChromID, Loci + FrameShift + RelLoci) & 0x0ff)) > eBaseT)
						break;
					FrameCodonIdx <<= 2;
					FrameCodonIdx |= LociBase;
					}
				if(LociBase <= eBaseT)
					pLociPValues->FrameShiftedRefCodons[FrameShift] = FrameCodonIdx;
				}
			}
		
		// if outputting as marker sequence then get SNP up/dnstream sequence and report
		int MarkerStartLoci;
		int MarkerSeqIdx;
		int AllelicIdx;
		tsSNPcnts* pMarkerBase;
		etSeqBase MarkerSequence[2000];
		etSeqBase* pMarkerSeq;
		int TotMarkerLociBases;
		double MarkerLociBaseProportion;
		int MarkerLen;
		int NumPolymorphicSites;
		if (m_hMarkerFile != -1)			// output marker sequences? 
			{
			// ensure putative marker sequence would be completely contained within the chromosome
			if (Loci < (uint32_t)m_Marker5Len)
				continue;
			if ((Loci + m_Marker3Len) >= m_pChromSNPs->ChromLen)
				continue;
			TotMarkerLociBases = pSNP->NumNonRefBases + pSNP->NumRefBases;
			MarkerLociBaseProportion = (double)pSNP->NumNonRefBases / TotMarkerLociBases;
			if (MarkerLociBaseProportion < 0.5)
				continue;

			NumPolymorphicSites = 0;
			MarkerLen = 1 + m_Marker5Len + m_Marker3Len;
			MarkerStartLoci = Loci - m_Marker5Len;
			pMarkerSeq = MarkerSequence;
			pMarkerBase = &m_pChromSNPs->Cnts[MarkerStartLoci];
			// check there are alignments covering the complete putative marker sequence
			// and that at any loci covered by the marker has a significant allelic base
			for (MarkerSeqIdx = 0; MarkerSeqIdx < MarkerLen; MarkerSeqIdx++, pMarkerBase++, pMarkerSeq++)
				{
				if ((TotMarkerLociBases = pMarkerBase->NumNonRefBases + pMarkerBase->NumRefBases) < m_MinSNPreads)	// must be at least enough reads covering to have confidence in base call
					break;
				MarkerLociBaseProportion = (double)pMarkerBase->NumNonRefBases / TotMarkerLociBases;
				if (MarkerLociBaseProportion <= m_MarkerPolyThres)													// if no more than polymorphic threshold then can simply accept RefBase
					{
					if (MarkerLociBaseProportion > 0.1)
						NumPolymorphicSites += 1;
					*pMarkerSeq = CSeqTrans::MapBase2Ascii(pMarkerBase->RefBase);
					continue;
					}
				// need to find a major allelic base - base must account for very high proportion of counts
				for (AllelicIdx = 0; AllelicIdx < 5; AllelicIdx++)
					if (pMarkerBase->NonRefBaseCnts[AllelicIdx] > 0 && (MarkerLociBaseProportion = ((double)pMarkerBase->NonRefBaseCnts[AllelicIdx] / TotMarkerLociBases)) >= (1.0 - m_MarkerPolyThres))
						{
						if (MarkerLociBaseProportion < 0.9)
							NumPolymorphicSites += 1;
						*pMarkerSeq = CSeqTrans::MapBase2Ascii(AllelicIdx);
						break;
						}
				if (AllelicIdx == 5)
					break;
				}
			if (MarkerSeqIdx != MarkerLen)			// only reporting SNPs which are consistent with reported markers
				continue;
			char SNPbase;
			char RefBase;
			RefBase = CSeqTrans::MapBase2Ascii(pSNP->RefBase);
			SNPbase = MarkerSequence[m_Marker5Len];
			if (RefBase == SNPbase)					// double check that the reference base is not being called as being the SNP base
				continue;
			MarkerSequence[MarkerLen] = '\0';

			// accepted marker sequence
			m_MarkerID += 1;
			pLociPValues->MarkerID = m_MarkerID;
			pLociPValues->NumPolymorphicSites = NumPolymorphicSites;
			// >MarkerNNN  Chrom StartLoci|MarkerLen|SNPLoci|Marker5Len,SNPbase|RefBase|NumPolymorphicSites
			LineLen += sprintf(&m_pszLineBuff[LineLen], ">Marker%d %s %d|%d|%d|%d|%c|%c|%d\n%s\n",
				m_MarkerID, szChromName, MarkerStartLoci, MarkerLen, Loci, m_Marker5Len, SNPbase, RefBase, NumPolymorphicSites, MarkerSequence);

			if ((LineLen + cMaxSeqLen) > cAllocLineBuffSize)
			{
				CUtility::RetryWrites(m_hMarkerFile, m_pszLineBuff, LineLen);
				LineLen = 0;
			}
		}

		if (m_hMarkerFile == -1)
			{
			pLociPValues->MarkerID = 0;
			pLociPValues->NumPolymorphicSites = 0;
			}
		PValue = 1.0 - Stats.Binomial(TotBases, pSNP->NumNonRefBases, LocalSeqErrRate);
		pLociPValues->PValue = PValue;
		pLociPValues->Loci = Loci;
		pLociPValues->Rank = 0;
		pLociPValues->LocalBkGndSubRate = LocalSeqErrRate;
		pLociPValues->LocalReads = LocTMM + LocTM;
		pLociPValues->LocalSubs = LocTMM;
		pLociPValues->NumReads = TotBases;
		pLociPValues->SNPcnts = *pSNP;
		pLociPValues->NumSubs = pSNP->NumNonRefBases;
		pLociPValues += 1;
		m_NumLociPValues += 1;
		}

	if (m_hMarkerFile != -1 && LineLen)
		{
		CUtility::RetryWrites(m_hMarkerFile, m_pszLineBuff, LineLen);
		LineLen = 0;
		}


	if (m_NumLociPValues == 0)
		{
		if(m_hWIGSpansFile != -1)
			{
			if(CoverageSegLen > 0 && PrevScaledLociCoverage > 0)
				{
				m_CovSegBuffIdx+=sprintf(&m_pszCovSegBuff[m_CovSegBuffIdx],"\n%s\t%u\t%u\tScaled:%u\t%u",szChromName,StartCoverageSegLoci,StartCoverageSegLoci+CoverageSegLen,PrevScaledLociCoverage,PrevScaledLociCoverage);
				CoverageSegLen = 0;
				}
			if(m_CovSegBuffIdx > 0)
				{
				if(!CUtility::RetryWrites(m_hWIGSpansFile, m_pszCovSegBuff, m_CovSegBuffIdx))
					{
					gDiagnostics.DiagOut(eDLInfo,gszProcName,"RetryWrites() error");
					return(eBSFerrWrite);
					}
				m_CovSegBuffIdx = 0;
				}
			}

		return(eBSFSuccess);
		}

	if (m_NumLociPValues > 1)
		m_mtqsort.qsort(m_pLociPValues, m_NumLociPValues, sizeof(tsLociPValues), SortLociPValues);
	pLociPValues = m_pLociPValues;
	NumSNPs = 0;
	for (Idx = 0; Idx < (int)m_NumLociPValues; Idx++, pLociPValues++)
		{
		AdjPValue = ((Idx + 1) / (double)m_NumLociPValues) * m_QValue;
		if (pLociPValues->PValue >= AdjPValue)
			break;
		NumSNPs += 1;
		pLociPValues->Rank = Idx + 1;
		}
	m_NumLociPValues = NumSNPs;
	if (m_NumLociPValues > 1)
		m_mtqsort.qsort(m_pLociPValues, m_NumLociPValues, sizeof(tsLociPValues), SortPValuesLoci);

	DiSNPBuffIdx = 0;
	TriSNPBuffIdx = 0;
	CurDiSNPLoci = 0;
	PrevDiSNPLoci = -1;
	CurTriSNPLoci = 0;
	PrevTriSNPLoci = -1;
	FirstTriSNPLoci = -1;
	TotNumDiSNPs = 0;
	TotNumTriSNPs = 0;
	LineLen = 0;
	pLociPValues = m_pLociPValues;
	for (Idx = 0; Idx < (int)m_NumLociPValues; Idx++, pLociPValues++)
	{
		m_TotNumSNPs += 1;
		RelRank = max(1, (int32_t)(999 - ((999 * pLociPValues->Rank) / m_NumLociPValues)));
		if (m_FMode == eFMbed)
		{
			LineLen += sprintf(&m_pszLineBuff[LineLen], "%s\t%d\t%d\tSNP_%d\t%d\t+\n",
				szChromName, pLociPValues->Loci, pLociPValues->Loci + 1, m_TotNumSNPs, RelRank);
		}
		else   // else could be either CSV or VCF
		{
			if (m_bSNPsVCF)
			{
				char szALTs[100];
				char szAltFreq[100];
				int AltOfs;
				int SNPPhred;
				int AltFreqOfs;
				int AltIdx;
				uint32_t CntsThres;		// only reporting cnts which are at least 10% of the highest non-ref base counts. 
										// otherwise too many noise cnt bases are reported 

				CntsThres = 0;
				for (AltIdx = 0; AltIdx < eBaseN; AltIdx++)
				{
					if (AltIdx == pLociPValues->SNPcnts.RefBase)
						continue;
					if (pLociPValues->SNPcnts.NonRefBaseCnts[AltIdx] > CntsThres)
						CntsThres = pLociPValues->SNPcnts.NonRefBaseCnts[AltIdx];
				}
				CntsThres = max((CntsThres + 5) / 10, (uint32_t)1);
				AltOfs = 0;
				AltFreqOfs = 0;
				for (AltIdx = 0; AltIdx < eBaseN; AltIdx++)
				{
					if (AltIdx == pLociPValues->SNPcnts.RefBase)
						continue;
					if (pLociPValues->SNPcnts.NonRefBaseCnts[AltIdx] >= CntsThres)
					{
						if (AltOfs > 0)
						{
							szALTs[AltOfs++] = ',';
							szAltFreq[AltFreqOfs++] = ',';
						}
						szALTs[AltOfs++] = CSeqTrans::MapBase2Ascii(AltIdx);
						szALTs[AltOfs] = '\0';
						AltFreqOfs += sprintf(&szAltFreq[AltFreqOfs], "%1.4f", (double)pLociPValues->SNPcnts.NonRefBaseCnts[AltIdx] / pLociPValues->NumReads);
					}
				}
				if (pLociPValues->PValue < 0.0000000001)
					SNPPhred = 100;
				else
					SNPPhred = (int)(0.5 + (10.0 * log10(1.0 / pLociPValues->PValue)));
				LineLen += sprintf(&m_pszLineBuff[LineLen], "%s\t%u\tSNP%d\t%c\t%s\t%d\tPASS\tAF=%s;DP=%d\n",
					szChromName, pLociPValues->Loci + 1, m_TotNumSNPs, CSeqTrans::MapBase2Ascii(pLociPValues->SNPcnts.RefBase),
					szALTs, SNPPhred, szAltFreq, pLociPValues->NumReads);
				}
			else
				{
				// for consistency now including refbase counts as if nonref - totals across all nonref bases will sum to be same as numreads covering the SNP loci
				pLociPValues->SNPcnts.NonRefBaseCnts[pLociPValues->SNPcnts.RefBase] = pLociPValues->NumReads - pLociPValues->NumSubs;
				LineLen += sprintf(&m_pszLineBuff[LineLen], "%d,\"SNP\",\"%s\",\"%s\",%d,%d,1,\"+\",%d,%f,%d,%d,\"%c\",%d,%d,%d,%d,%d,%f,%d,%d,%d,%d",
					m_TotNumSNPs, m_szTargSpecies, szChromName, pLociPValues->Loci, pLociPValues->Loci, RelRank, pLociPValues->PValue,
					pLociPValues->NumReads, pLociPValues->NumSubs,
					CSeqTrans::MapBase2Ascii(pLociPValues->SNPcnts.RefBase),
					pLociPValues->SNPcnts.NonRefBaseCnts[0], pLociPValues->SNPcnts.NonRefBaseCnts[1], pLociPValues->SNPcnts.NonRefBaseCnts[2], pLociPValues->SNPcnts.NonRefBaseCnts[3], pLociPValues->SNPcnts.NonRefBaseCnts[4],
					pLociPValues->LocalBkGndSubRate, pLociPValues->LocalReads, pLociPValues->LocalSubs, pLociPValues->MarkerID, pLociPValues->NumPolymorphicSites);
				if(m_bXCSVFrameShifts)
					{
					char szCodon[4];
					for(int FrameIdx = 0; FrameIdx < 3; FrameIdx++)
						{
						int Codon = pLociPValues->FrameShiftedRefCodons[FrameIdx];
						if(Codon < 0)
							strcpy(szCodon, "NNN");
						else
							for(int Base = 0; Base < 3; Base++)
								{
								switch(Codon & 0x030) {
									case 0x00: szCodon[Base] = 'A'; break;
									case 0x10: szCodon[Base] = 'C'; break;
									case 0x20: szCodon[Base] = 'G'; break;
									case 0x30: szCodon[Base] = 'T'; break;
									}
								Codon <<= 2;
								}
						szCodon[3] = '\0';
						LineLen += sprintf(&m_pszLineBuff[LineLen], ",\"%s\"", szCodon);
						for(int CodonIdx = 0; CodonIdx < 64; CodonIdx++)
							LineLen += sprintf(&m_pszLineBuff[LineLen],",%d", pLociPValues->FrameShiftedCodons[FrameIdx][CodonIdx]);
						}
					}
				LineLen += sprintf(&m_pszLineBuff[LineLen],"\n");
				}
			}
		if ((LineLen + cMaxSeqLen + 1) > cAllocLineBuffSize)
			{
			CUtility::RetryWrites(m_hSNPfile, m_pszLineBuff, LineLen);
			LineLen = 0;
			}

		if (gProcessingID > 0)
			{
			sMonoSNP.MonoSnpPID = m_TotNumSNPs;						// SNP instance, processing instance unique
			strcpy(sMonoSNP.szElType, "SNP");						// SNP type
			strcpy(sMonoSNP.szSpecies, m_szTargSpecies);			// SNP located for alignments against this target/species assembly	
			strcpy(sMonoSNP.szChrom, szChromName);					// SNP is on this chrom
			sMonoSNP.StartLoci = pLociPValues->Loci;				// offset (0..N) at which SNP located
			sMonoSNP.EndLoci = pLociPValues->Loci;					// offset (0..N) at which SNP located - allowing for future polymorphic varation covering multiple bases
			sMonoSNP.Len = 1;										// polymorphic variation is of this length
			sMonoSNP.szStrand[0] = '+'; sMonoSNP.szStrand[1] = '\0'; // SNP relative to this strand
			sMonoSNP.Rank = RelRank;								// ranking confidence in thisSNP - min 0, max 1000
			sMonoSNP.PValue = pLociPValues->PValue;					// probability of this SNP being a false positive
			sMonoSNP.Bases = pLociPValues->NumReads;				// total number of bases aligning over the SNP loci
			sMonoSNP.Mismatches = pLociPValues->NumSubs;			// aligned bases were aligning with this many mismatches
			sMonoSNP.szRefBase[0] = CSeqTrans::MapBase2Ascii(pLociPValues->SNPcnts.RefBase); sMonoSNP.szRefBase[1] = '\0';			// target sequence base at the SNP locai
			sMonoSNP.MMBaseA = pLociPValues->SNPcnts.NonRefBaseCnts[0];			// this many mismatched bases were A
			sMonoSNP.MMBaseC = pLociPValues->SNPcnts.NonRefBaseCnts[1];			// this many mismatched bases were C
			sMonoSNP.MMBaseG = pLociPValues->SNPcnts.NonRefBaseCnts[2];			// this many mismatched bases were G
			sMonoSNP.MMBaseT = pLociPValues->SNPcnts.NonRefBaseCnts[3];			// this many mismatched bases were T
			sMonoSNP.MMBaseN = pLociPValues->SNPcnts.NonRefBaseCnts[4];			// this many mismatched bases were N
			sMonoSNP.BackgroundSubRate = pLociPValues->LocalBkGndSubRate;		// background substitution rate within a window centered at SNP loci
			sMonoSNP.TotWinBases = pLociPValues->LocalReads;					// total number of bases within centeredwindow
			sMonoSNP.TotWinMismatches = pLociPValues->LocalSubs;				// total number of mismatched bases within centered window
			sMonoSNP.MarkerID = pLociPValues->MarkerID;							// marker identifier
			sMonoSNP.NumPolymorphicSites = pLociPValues->NumPolymorphicSites;	// number of polymorphic sites within marker
			gSQLiteSummaries.AddMonoSNP(gExperimentID, gProcessingID, &sMonoSNP);
			}

#ifdef _DISNPS_
		if (!m_bIsSOLiD && m_hDiSNPfile != -1)
			{
			// try to find all reads which are overlapping this SNP plus the prev within m_MaxDiSNPSep (max 300bp)
			tsReadHit* pCurOverlappingRead;
			uint8_t PrevDiSNPBase;
			uint8_t CurDiSNPBase;
			uint8_t FirstTriSNPBase;
			uint8_t PrevTriSNPBase;
			uint8_t CurTriSNPBase;

			int NumHaplotypes;
			int HaplotypeCntThres;
			int NumReadsOverlapping;
			int NumReadsAntisense;
			int DiSNPIdx;
			int DiSNPCnts[64];

			Loci = pLociPValues->Loci;
			CurDiSNPLoci = Loci;
			pSNP = &m_pChromSNPs->Cnts[Loci];
			if (PrevDiSNPLoci != -1 && CurDiSNPLoci > 0 && ((CurDiSNPLoci - PrevDiSNPLoci) <= m_MaxDiSNPSep))
				{
				NumReadsOverlapping = 0;
				NumReadsAntisense = 0;
				memset(DiSNPCnts, 0, sizeof(DiSNPCnts));
				memset(PrevSNPs, 0, sizeof(PrevSNPs));
				PrevSNPs[0].RefBase = (&m_pChromSNPs->Cnts[PrevDiSNPLoci])->RefBase;
				PrevSNPs[1].RefBase = pSNP->RefBase;

				while ((pCurOverlappingRead = IterateReadsOverlapping(false, m_pChromSNPs, PrevDiSNPLoci, CurDiSNPLoci)) != nullptr)
					{
					// get bases at both SNP loci for current overlapping read
					PrevDiSNPBase = AdjAlignSNPBase(pCurOverlappingRead, m_pChromSNPs->ChromID, PrevDiSNPLoci);
					if (PrevDiSNPBase > eBaseT)
						continue;
					CurDiSNPBase = AdjAlignSNPBase(pCurOverlappingRead, m_pChromSNPs->ChromID, CurDiSNPLoci);
					if (CurDiSNPBase > eBaseT)
						continue;
					PrevSNPs[0].NonRefBaseCnts[PrevDiSNPBase] += 1;
					PrevSNPs[1].NonRefBaseCnts[CurDiSNPBase] += 1;
					if (PrevDiSNPBase == PrevSNPs[0].RefBase)
						PrevSNPs[0].NumRefBases += 1;
					else
						PrevSNPs[0].NumNonRefBases += 1;
					if (CurDiSNPBase == PrevSNPs[1].RefBase)
						PrevSNPs[1].NumRefBases += 1;
					else
						PrevSNPs[1].NumNonRefBases += 1;
					NumReadsOverlapping += 1;
					if (pCurOverlappingRead->HitLoci.Hit.Seg[0].Strand == '-')
						NumReadsAntisense += 1;
					DiSNPIdx = ((PrevDiSNPBase & 0x03) << 2) | (CurDiSNPBase & 0x03);
					DiSNPCnts[DiSNPIdx] += 1;
					}
				NumHaplotypes = 0;	// haplotype calling very simplistic very simplistic - calling haplotype if was called a SNP loci and allele coverage >= than 10% of total counts (min 5)
				if (NumReadsOverlapping >= m_MinSNPreads)
					{
					HaplotypeCntThres = max(5, (NumReadsOverlapping+5) / 10); // if fewer than 5 reads overlapping with same allele then simply too noisy to use as a haplotype!
					for (DiSNPIdx = 0; DiSNPIdx < 16; DiSNPIdx++)
						if (DiSNPCnts[DiSNPIdx] >= HaplotypeCntThres)
							NumHaplotypes += 1;
						else
							DiSNPCnts[DiSNPIdx] = 0;
					}

				if(NumHaplotypes >= 2)		// only reporting if at least 2 DiSNP localised haplotypes
					{
					TotNumDiSNPs += 1;
					DiSNPBuffIdx += sprintf(&szDiSNPs[DiSNPBuffIdx], "%d,\"DiSNPs\",\"%s\",\"%s\",", TotNumDiSNPs, m_szTargSpecies, szChromName);

					if (gProcessingID > 0)
						{
						sDiSNP.DiSnpPID = TotNumDiSNPs;		// SNP instance, processing instance unique
						strcpy(sDiSNP.szElType, "DiSNPs");		// SNP type
						strcpy(sDiSNP.szSpecies, m_szTargSpecies);		// SNP located for alignments againts this target/species assembly	
						strcpy(sDiSNP.szChrom, szChromName);		// SNP is on this chrom
						}

					char SNPrefBases[2];
					int RefBasesIdx;
					for (RefBasesIdx = 0; RefBasesIdx <= 1; RefBasesIdx++)
						{
						switch (PrevSNPs[RefBasesIdx].RefBase) {
							case eBaseA:
								SNPrefBases[RefBasesIdx] = 'a';
								break;
							case eBaseC:
								SNPrefBases[RefBasesIdx] = 'c';
								break;
							case eBaseG:
								SNPrefBases[RefBasesIdx] = 'g';
								break;
							case eBaseT:
								SNPrefBases[RefBasesIdx] = 't';
								break;
							case eBaseN:
								SNPrefBases[RefBasesIdx] = 'n';
								break;
							};
						switch (RefBasesIdx) {
							case 0:
								DiSNPBuffIdx += sprintf(&szDiSNPs[DiSNPBuffIdx], "%d,", PrevDiSNPLoci);
								if (gProcessingID > 0)
									{
									sDiSNP.SNP1Loci = PrevDiSNPLoci;
									sDiSNP.szSNP1RefBase[0] = SNPrefBases[RefBasesIdx];
									sDiSNP.szSNP1RefBase[1] = '\0';
									sDiSNP.SNP1BaseAcnt = PrevSNPs[RefBasesIdx].NonRefBaseCnts[0];
									sDiSNP.SNP1BaseCcnt = PrevSNPs[RefBasesIdx].NonRefBaseCnts[1];
									sDiSNP.SNP1BaseGcnt = PrevSNPs[RefBasesIdx].NonRefBaseCnts[2];
									sDiSNP.SNP1BaseTcnt = PrevSNPs[RefBasesIdx].NonRefBaseCnts[3];
									sDiSNP.SNP1BaseNcnt = 0;
									}
								break;
							case 1:
								DiSNPBuffIdx += sprintf(&szDiSNPs[DiSNPBuffIdx], "%d,", CurDiSNPLoci);
								if (gProcessingID > 0)
									{
									sDiSNP.SNP2Loci = CurDiSNPLoci;
									sDiSNP.szSNP2RefBase[0] = SNPrefBases[RefBasesIdx];
									sDiSNP.szSNP2RefBase[1] = '\0';
									sDiSNP.SNP2BaseAcnt = PrevSNPs[RefBasesIdx].NonRefBaseCnts[0];
									sDiSNP.SNP2BaseCcnt = PrevSNPs[RefBasesIdx].NonRefBaseCnts[1];
									sDiSNP.SNP2BaseGcnt = PrevSNPs[RefBasesIdx].NonRefBaseCnts[2];
									sDiSNP.SNP2BaseTcnt = PrevSNPs[RefBasesIdx].NonRefBaseCnts[3];
									sDiSNP.SNP2BaseNcnt = 0;
									}
								break;
							}

						DiSNPBuffIdx += sprintf(&szDiSNPs[DiSNPBuffIdx], "\"%c\",%d,%d,%d,%d,0,", SNPrefBases[RefBasesIdx], PrevSNPs[RefBasesIdx].NonRefBaseCnts[0], PrevSNPs[RefBasesIdx].NonRefBaseCnts[1], PrevSNPs[RefBasesIdx].NonRefBaseCnts[2], PrevSNPs[RefBasesIdx].NonRefBaseCnts[3]);
						}

					DiSNPBuffIdx += sprintf(&szDiSNPs[DiSNPBuffIdx], "%d,%d,%d", NumReadsOverlapping, NumReadsAntisense, NumHaplotypes);
					if (gProcessingID > 0)
						{
						sDiSNP.Depth = NumReadsOverlapping;
						sDiSNP.Antisense = NumReadsAntisense;
						sDiSNP.Haplotypes = NumHaplotypes;
						}


					for (DiSNPIdx = 0; DiSNPIdx < 16; DiSNPIdx++)
						{
						DiSNPBuffIdx += sprintf(&szDiSNPs[DiSNPBuffIdx], ",%d", DiSNPCnts[DiSNPIdx]);
						if (gProcessingID > 0)
							sDiSNP.HaplotypeCnts[DiSNPIdx] = DiSNPCnts[DiSNPIdx];
						}
					DiSNPBuffIdx += sprintf(&szDiSNPs[DiSNPBuffIdx], "\n");
					if (gProcessingID > 0)
						gSQLiteSummaries.AddDiSNP(gExperimentID, gProcessingID, &sDiSNP);

					if ((DiSNPBuffIdx + 200) > sizeof(szDiSNPs))
						{
						if(!CUtility::RetryWrites(m_hDiSNPfile, szDiSNPs, DiSNPBuffIdx))
							{
							gDiagnostics.DiagOut(eDLInfo,gszProcName,"RetryWrites() error");
							return(eBSFerrWrite);
							}
						DiSNPBuffIdx = 0;
						}
					}
				}
			// DiSNP processing completed

			// starting on TriSNPs
			CurTriSNPLoci = Loci;
			if (FirstTriSNPLoci != -1 && PrevTriSNPLoci > 0 && CurTriSNPLoci > 0 && ((CurTriSNPLoci - FirstTriSNPLoci) <= m_MaxDiSNPSep))
				{
				NumReadsOverlapping = 0;
				NumReadsAntisense = 0;
				memset(DiSNPCnts, 0, sizeof(DiSNPCnts));

				memset(PrevSNPs, 0, sizeof(PrevSNPs));
				PrevSNPs[0].RefBase = (&m_pChromSNPs->Cnts[FirstTriSNPLoci])->RefBase;
				PrevSNPs[1].RefBase = (&m_pChromSNPs->Cnts[PrevTriSNPLoci])->RefBase;
				PrevSNPs[2].RefBase = pSNP->RefBase;

				while ((pCurOverlappingRead = IterateReadsOverlapping(true, m_pChromSNPs, FirstTriSNPLoci, CurTriSNPLoci)) != nullptr)
					{
					// get bases at all three SNP loci
					FirstTriSNPBase = AdjAlignSNPBase(pCurOverlappingRead, m_pChromSNPs->ChromID, FirstTriSNPLoci);
					if (FirstTriSNPBase > eBaseT)
						continue;
					PrevTriSNPBase = AdjAlignSNPBase(pCurOverlappingRead, m_pChromSNPs->ChromID, PrevTriSNPLoci);
					if (PrevTriSNPBase > eBaseT)
						continue;
					CurTriSNPBase = AdjAlignSNPBase(pCurOverlappingRead, m_pChromSNPs->ChromID, CurTriSNPLoci);
					if (CurTriSNPBase > eBaseT)
						continue;

					PrevSNPs[0].NonRefBaseCnts[FirstTriSNPBase] += 1;
					PrevSNPs[1].NonRefBaseCnts[PrevTriSNPBase] += 1;
					PrevSNPs[2].NonRefBaseCnts[CurTriSNPBase] += 1;
					if (FirstTriSNPBase == PrevSNPs[0].RefBase)
						PrevSNPs[0].NumRefBases += 1;
					else
						PrevSNPs[0].NumNonRefBases += 1;
					if (PrevTriSNPBase == PrevSNPs[1].RefBase)
						PrevSNPs[1].NumRefBases += 1;
					else
						PrevSNPs[1].NumNonRefBases += 1;
					if (CurTriSNPBase == PrevSNPs[2].RefBase)
						PrevSNPs[2].NumRefBases += 1;
					else
						PrevSNPs[2].NumNonRefBases += 1;

					NumReadsOverlapping += 1;
					if (pCurOverlappingRead->HitLoci.Hit.Seg[0].Strand == '-')
						NumReadsAntisense += 1;
					DiSNPIdx = ((FirstTriSNPBase & 0x03) << 4) | ((PrevTriSNPBase & 0x03) << 2) | (CurTriSNPBase & 0x03);
					DiSNPCnts[DiSNPIdx] += 1;
					}

				NumHaplotypes = 0;	// haplotype calling very simplistic very simplistic - calling haplotype if was called a SNP loci and allele coverage >= than 10% of total counts (min 5)
				if (NumReadsOverlapping >= m_MinSNPreads)
					{
					HaplotypeCntThres = max(5, (NumReadsOverlapping+5) / 10); // if fewer than 5 reads overlapping with same allele then simply too noisy to use as a haplotype!
					for (DiSNPIdx = 0; DiSNPIdx < 64; DiSNPIdx++)
						if (DiSNPCnts[DiSNPIdx] >= HaplotypeCntThres)
							NumHaplotypes += 1;
						else
							DiSNPCnts[DiSNPIdx] = 0;
					}

				if(NumHaplotypes >= 3)		// only reporting if at least 3 TriSNP localised haplotypes
					{
					TotNumTriSNPs += 1;

					TriSNPBuffIdx += sprintf(&szTriSNPs[TriSNPBuffIdx], "%d,\"TriSNPs\",\"%s\",\"%s\",", TotNumTriSNPs, m_szTargSpecies, szChromName);
					if (gProcessingID > 0)
						{
						sTriSNP.TriSnpPID = TotNumTriSNPs;				// SNP instance, processing instance unique
						strcpy(sTriSNP.szElType, "TriSNPs");			// SNP type
						strcpy(sTriSNP.szSpecies, m_szTargSpecies);		// SNP located for alignments against this target/species assembly	
						strcpy(sTriSNP.szChrom, szChromName);			// SNP is on this chrom
						}
					char SNPrefBases[3];
					int RefBasesIdx;

					for (RefBasesIdx = 0; RefBasesIdx <= 2; RefBasesIdx++)
						{
						switch (PrevSNPs[RefBasesIdx].RefBase) {
							case eBaseA:
								SNPrefBases[RefBasesIdx] = 'a';
								break;
							case eBaseC:
								SNPrefBases[RefBasesIdx] = 'c';
								break;
							case eBaseG:
								SNPrefBases[RefBasesIdx] = 'g';
								break;
							case eBaseT:
								SNPrefBases[RefBasesIdx] = 't';
								break;
							case eBaseN:
								SNPrefBases[RefBasesIdx] = 'n';
								break;
							};
					switch (RefBasesIdx) {
						case 0:
							TriSNPBuffIdx += sprintf(&szTriSNPs[TriSNPBuffIdx], "%d,", FirstTriSNPLoci);
							if (gProcessingID > 0)
								{
								sTriSNP.SNP1Loci = FirstTriSNPLoci;
								sTriSNP.szSNP1RefBase[0] = SNPrefBases[RefBasesIdx];
								sTriSNP.szSNP1RefBase[1] = '\0';
								sTriSNP.SNP1BaseAcnt = PrevSNPs[RefBasesIdx].NonRefBaseCnts[0];
								sTriSNP.SNP1BaseCcnt = PrevSNPs[RefBasesIdx].NonRefBaseCnts[1];
								sTriSNP.SNP1BaseGcnt = PrevSNPs[RefBasesIdx].NonRefBaseCnts[2];
								sTriSNP.SNP1BaseTcnt = PrevSNPs[RefBasesIdx].NonRefBaseCnts[3];
								sTriSNP.SNP1BaseNcnt = 0;
								}
							break;
						case 1:
							TriSNPBuffIdx += sprintf(&szTriSNPs[TriSNPBuffIdx], "%d,", PrevTriSNPLoci);
							if (gProcessingID > 0)
								{
								sTriSNP.SNP2Loci = PrevTriSNPLoci;
								sTriSNP.szSNP2RefBase[0] = SNPrefBases[RefBasesIdx];
								sTriSNP.szSNP2RefBase[1] = '\0';
								sTriSNP.SNP2BaseAcnt = PrevSNPs[RefBasesIdx].NonRefBaseCnts[0];
								sTriSNP.SNP2BaseCcnt = PrevSNPs[RefBasesIdx].NonRefBaseCnts[1];
								sTriSNP.SNP2BaseGcnt = PrevSNPs[RefBasesIdx].NonRefBaseCnts[2];
								sTriSNP.SNP2BaseTcnt = PrevSNPs[RefBasesIdx].NonRefBaseCnts[3];
								sTriSNP.SNP2BaseNcnt = 0;
								}
							break;
						case 2:
							TriSNPBuffIdx += sprintf(&szTriSNPs[TriSNPBuffIdx], "%d,", CurTriSNPLoci);
							if (gProcessingID > 0)
								{
								sTriSNP.SNP3Loci = CurTriSNPLoci;
								sTriSNP.szSNP3RefBase[0] = SNPrefBases[RefBasesIdx];
								sTriSNP.szSNP3RefBase[1] = '\0';
								sTriSNP.SNP3BaseAcnt = PrevSNPs[RefBasesIdx].NonRefBaseCnts[0];
								sTriSNP.SNP3BaseCcnt = PrevSNPs[RefBasesIdx].NonRefBaseCnts[1];
								sTriSNP.SNP3BaseGcnt = PrevSNPs[RefBasesIdx].NonRefBaseCnts[2];
								sTriSNP.SNP3BaseTcnt = PrevSNPs[RefBasesIdx].NonRefBaseCnts[3];
								sTriSNP.SNP3BaseNcnt = 0;
								}
							break;
						}

					TriSNPBuffIdx += sprintf(&szTriSNPs[TriSNPBuffIdx], "\"%c\",%d,%d,%d,%d,0,", SNPrefBases[RefBasesIdx], PrevSNPs[RefBasesIdx].NonRefBaseCnts[0], PrevSNPs[RefBasesIdx].NonRefBaseCnts[1], PrevSNPs[RefBasesIdx].NonRefBaseCnts[2], PrevSNPs[RefBasesIdx].NonRefBaseCnts[3]);
					}

				TriSNPBuffIdx += sprintf(&szTriSNPs[TriSNPBuffIdx], "%d,%d,%d", NumReadsOverlapping, NumReadsAntisense, NumHaplotypes);
				if (gProcessingID > 0)
					{
					sTriSNP.Depth = NumReadsOverlapping;
					sTriSNP.Antisense = NumReadsAntisense;
					sTriSNP.Haplotypes = NumHaplotypes;
					}
				for (DiSNPIdx = 0; DiSNPIdx < 64; DiSNPIdx++)
					{
					TriSNPBuffIdx += sprintf(&szTriSNPs[TriSNPBuffIdx], ",%d", DiSNPCnts[DiSNPIdx]);
					if (gProcessingID > 0)
						sTriSNP.HaplotypeCnts[DiSNPIdx] = DiSNPCnts[DiSNPIdx];
					}
				TriSNPBuffIdx += sprintf(&szTriSNPs[TriSNPBuffIdx], "\n");
				if (gProcessingID > 0)
					gSQLiteSummaries.AddTriSNP(gExperimentID, gProcessingID, &sTriSNP);
				if ((TriSNPBuffIdx + 500) > sizeof(szTriSNPs))
					{
					CUtility::RetryWrites(m_hTriSNPfile, szTriSNPs, TriSNPBuffIdx);
					TriSNPBuffIdx = 0;
					}
				}
			}
		PrevDiSNPLoci = CurDiSNPLoci;
		FirstTriSNPLoci = PrevTriSNPLoci;
		PrevTriSNPLoci = CurTriSNPLoci;
		}
#endif

	if (m_hSNPCentsfile != -1)
		{
		// get 4bases up/dn stream from loci with SNP and use these to inc centroid counts of from/to counts
		if (pLociPValues->Loci >= cSNPCentfFlankLen && pLociPValues->Loci < (m_pChromSNPs->ChromLen - cSNPCentfFlankLen))
			{
			m_pSfxArray->GetSeq(m_pChromSNPs->ChromID, pLociPValues->Loci - (uint32_t)cSNPCentfFlankLen, SNPFlanks, cSNPCentroidLen);
			pSNPFlank = &SNPFlanks[cSNPCentroidLen - 1];
			SNPCentroidIdx = 0;
			for (SNPFlankIdx = 0; SNPFlankIdx < cSNPCentroidLen; SNPFlankIdx++, pSNPFlank--)
				{
				Base = *pSNPFlank & 0x07;
				if (Base > eBaseT)
					break;
				SNPCentroidIdx |= (Base << (SNPFlankIdx * 2));
				}
			if (SNPFlankIdx != cSNPCentroidLen)
				continue;

			pSNP = &m_pChromSNPs->Cnts[pLociPValues->Loci];
			pCentroid = &m_pSNPCentroids[SNPCentroidIdx];
			pCentroid->CentroidID = SNPCentroidIdx;
			pCentroid->RefBaseCnt += pSNP->NumRefBases;
			pCentroid->NonRefBaseCnts[0] += pSNP->NonRefBaseCnts[0];
			pCentroid->NonRefBaseCnts[1] += pSNP->NonRefBaseCnts[1];
			pCentroid->NonRefBaseCnts[2] += pSNP->NonRefBaseCnts[2];
			pCentroid->NonRefBaseCnts[3] += pSNP->NonRefBaseCnts[3];
			pCentroid->NonRefBaseCnts[4] += pSNP->NonRefBaseCnts[4];
			pCentroid->NumSNPs += 1;
			}
		}
	}
CompleteWIGSpan(true);

if (m_hDiSNPfile != -1 && DiSNPBuffIdx > 0)
	{
	if(!CUtility::RetryWrites(m_hDiSNPfile, szDiSNPs, DiSNPBuffIdx))
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"RetryWrites() error");
		return(eBSFerrWrite);
		}
	DiSNPBuffIdx = 0;
	}
if (m_hTriSNPfile != -1 && TriSNPBuffIdx > 0)
	{
	if(!CUtility::RetryWrites(m_hTriSNPfile, szTriSNPs, TriSNPBuffIdx))
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"RetryWrites() error");
		return(eBSFerrWrite);
		}
	DiSNPBuffIdx = 0;
	}

if (LineLen)
	if(!CUtility::RetryWrites(m_hSNPfile, m_pszLineBuff, LineLen))
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"RetryWrites() error");
		return(eBSFerrWrite);
		}
return(eBSFSuccess);
}



int
CKAligner::ProcessSNPs(void)
{
int Rslt;
int LineLen;

uint32_t SeqIdx;
etSeqBase ReadSeq[cMaxFastQSeqLen+10000];	// to hold sequence (sans quality scores) for current read
etSeqBase AssembSeq[cMaxFastQSeqLen+ 10000];	// to hold targeted genome assembly sequence

etSeqBase TargBases[3];
etSeqBase ReadBase;
tsSNPcnts *pSNP;
uint8_t *pSeqVal;
etSeqBase *pReadSeq;
etSeqBase *pAssembSeq;
tsSegLoci *pSeg;
tsReadHit *pReadHit;
tBSFEntryID PrevTargEntry;
uint32_t ChromLen;
uint32_t HitLoci;
uint32_t MatchLen;
uint32_t PrevMMChromID;
uint32_t PrevMMLoci;

if(!m_bPackedBaseAlleles)
	{
	if(m_FMode == eFMbed)
		{
		LineLen = sprintf(m_pszLineBuff,"track type=bed name=\"%s_SNPs\" description=\"%s SNPs\"\n",m_pszTrackTitle,m_pszTrackTitle);
		CUtility::RetryWrites(m_hSNPfile,m_pszLineBuff,LineLen);
		LineLen = 0;
		}
	else			// else must be either CSV or VCF
		{
		if(m_bSNPsVCF)
			{
			LineLen = sprintf(m_pszLineBuff,"##fileformat=VCFv4.1\n##source=kalign%s\n##reference=%s\n##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n",
														kit4bversion,m_pszSfxFile);
			LineLen += sprintf(&m_pszLineBuff[LineLen],"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
			}
		else
			{
			LineLen = sprintf(m_pszLineBuff,"\"SNP_ID\",\"ElType\",\"Species\",\"Chrom\",\"StartLoci\",\"EndLoci\",\"Len\",\"Strand\",\"Rank\",\"PValue\",\"Bases\",\"Mismatches\",\"RefBase\",\"MMBaseA\",\"MMBaseC\",\"MMBaseG\",\"MMBaseT\",\"MMBaseN\",\"BackgroundSubRate\",\"TotWinBases\",\"TotWinMismatches\",\"MarkerID\",\"NumPolymorphicSites\"");
			if(m_bXCSVFrameShifts)
				{
				char szCodon[4];
				for(int FrameShift = 0; FrameShift < 3; FrameShift++)
					{
					LineLen += sprintf(&m_pszLineBuff[LineLen], ",RefFrame:%d", FrameShift);
					for(int CodonCnts = 0; CodonCnts < 64; CodonCnts++)
						{
						int Codon = CodonCnts;
						for(int Base = 0; Base < 3 ; Base++)
							{
							switch(Codon & 0x030) {
								case 0x00: szCodon[Base] = 'A'; break;
								case 0x10: szCodon[Base] = 'C'; break;
								case 0x20: szCodon[Base] = 'G'; break;
								case 0x30: szCodon[Base] = 'T'; break;
								}
							Codon <<= 2;
							}
						szCodon[3] = '\0';
						LineLen += sprintf(&m_pszLineBuff[LineLen], ",\"RFS%d: %s\"", FrameShift,szCodon);
						}
					}
				}
			LineLen += sprintf(&m_pszLineBuff[LineLen],"\n");
			}
		CUtility::RetryWrites(m_hSNPfile,m_pszLineBuff,LineLen);
		LineLen = 0;
		}

	if(m_hWIGSpansFile != -1)
		{
		LineLen = sprintf(m_pszLineBuff,"track type=wiggle_0 name=\"Coverage\" description=\"Alignment Segment Coverage\" useScore=1\n");
		if(!CUtility::RetryWrites(m_hWIGSpansFile,m_pszLineBuff,LineLen))
			{
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"RetryWrites() error");
			return(eBSFerrWrite);
			}
		LineLen = 0;
		}

	if(m_hDiSNPfile != -1)
		{
		int Idx;
		char szDiSNPs[3];
		LineLen = sprintf(m_pszLineBuff,"\"DiSNPs_ID\",\"ElType\",\"Species\",\"Chrom\",\"SNP1Loci\",\"SNP1RefBase\",\"SNP1BaseAcnt\",\"SNP1BaseCcnt\",\"SNP1BaseGcnt\",\"SNP1BaseTcnt\",\"SNP1BaseNcnt\",\"SNP2Loci\",\"SNP2RefBase\",\"SNP2BaseAcnt\",\"SNP2BaseCcnt\",\"SNP2BaseGcnt\",\"SNP2BaseTcnt\",\"SNP2BaseNcnt\",\"Depth\",\"Antisense\",\"Haplotypes\"");
		for(Idx = 0; Idx < 16; Idx++)
			{
			switch(Idx & 0x03) {
				case 0: szDiSNPs[1] = 'a'; break;
				case 1: szDiSNPs[1] = 'c'; break;
				case 2: szDiSNPs[1] = 'g'; break;
				case 3: szDiSNPs[1] = 't'; break;
				}
			switch((Idx >> 2) & 0x03) {
				case 0: szDiSNPs[0] = 'a'; break;
				case 1: szDiSNPs[0] = 'c'; break;
				case 2: szDiSNPs[0] = 'g'; break;
				case 3: szDiSNPs[0] = 't'; break;
				}
			szDiSNPs[2] = '\0';
			LineLen += sprintf(&m_pszLineBuff[LineLen],",\"%s\"",szDiSNPs);
			}
		LineLen += sprintf(&m_pszLineBuff[LineLen],"\n");
		if(!CUtility::RetryWrites(m_hDiSNPfile,m_pszLineBuff,LineLen))
			{
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"RetryWrites() error");
			return(eBSFerrWrite);
			}
		LineLen = 0;
		}

	if(m_hTriSNPfile != -1)
		{
		int Idx;
		char szTriSNPs[4];
		LineLen = sprintf(m_pszLineBuff,"\"TriSNPs_ID\",\"ElType\",\"Species\",\"Chrom\",\"SNP1Loci\",\"SNP1RefBase\",\"SNP1BaseAcnt\",\"SNP1BaseCcnt\",\"SNP1BaseGcnt\",\"SNP1BaseTcnt\",\"SNP1BaseNcnt\",\"SNP2Loci\",\"SNP2RefBase\",\"SNP2BaseAcnt\",\"SNP2BaseCcnt\",\"SNP2BaseGcnt\",\"SNP2BaseTcnt\",\"SNP2BaseNcnt\",\"SNP3Loci\",\"SNP3RefBase\",\"SNP3BaseAcnt\",\"SNP3BaseCcnt\",\"SNP3BaseGcnt\",\"SNP3BaseTcnt\",\"SNP3BaseNcnt\",\"Depth\",\"Antisense\",\"Haplotypes\"");
		for(Idx = 0; Idx < 64; Idx++)
			{
			switch(Idx & 0x03) {
				case 0: szTriSNPs[2] = 'a'; break;
				case 1: szTriSNPs[2] = 'c'; break;
				case 2: szTriSNPs[2] = 'g'; break;
				case 3: szTriSNPs[2] = 't'; break;
				}
			switch((Idx >> 2) & 0x03) {
				case 0: szTriSNPs[1] = 'a'; break;
				case 1: szTriSNPs[1] = 'c'; break;
				case 2: szTriSNPs[1] = 'g'; break;
				case 3: szTriSNPs[1] = 't'; break;
				}
			switch((Idx >> 4) & 0x03) {
				case 0: szTriSNPs[0] = 'a'; break;
				case 1: szTriSNPs[0] = 'c'; break;
				case 2: szTriSNPs[0] = 'g'; break;
				case 3: szTriSNPs[0] = 't'; break;
				}
			szTriSNPs[3] = '\0';
			LineLen += sprintf(&m_pszLineBuff[LineLen],",\"%s\"",szTriSNPs);
			}
		LineLen += sprintf(&m_pszLineBuff[LineLen],"\n");
		CUtility::RetryWrites(m_hTriSNPfile,m_pszLineBuff,LineLen);
		LineLen = 0;
		}

	// need to check that there are accepted aligned reads to process for SNPs!!!
	if(m_hSNPCentsfile != -1)
		{
		if(m_pSNPCentroids == nullptr)
			{
			int CentroidIdx;
			tsSNPCentroid *pCentroid;
			if((m_pSNPCentroids = new tsSNPCentroid[cSNPCentroidEls + 16])==nullptr)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessSNPs: Memory allocation of %d SNP centroid elements failed",cSNPCentroidEls + 16);
				Reset(false);
				return(eBSFerrMem);
				}
			memset(m_pSNPCentroids,0,sizeof(tsSNPCentroid) * (cSNPCentroidEls+16));
			pCentroid = m_pSNPCentroids;
			for(CentroidIdx = 1; CentroidIdx <= cSNPCentroidEls; CentroidIdx++,pCentroid++)
				pCentroid->CentroidID = CentroidIdx;
			}
		}
	}

pReadHit = nullptr;
LineLen = 0;
PrevTargEntry = 0;

m_TotNumSNPs = 0;
while((pReadHit = IterSortedReads(pReadHit))!=nullptr)
	{
	if(pReadHit->NAR == eNARAccepted)
		{
		if(pReadHit->HitLoci.Hit.FlgInDel || pReadHit->HitLoci.Hit.FlgSplice)
			continue;

		pSeg = &pReadHit->HitLoci.Hit.Seg[0];
		if(pSeg->ChromID != (uint32_t)PrevTargEntry)
			{
			if(m_pChromSNPs != nullptr)
				{
				// this is where the SNPs for the previously processed chrom need to saved off as new chrom is about to be processed
				m_pChromSNPs->MeanReadLen = (uint32_t)(((m_pChromSNPs->TotReadLen + m_pChromSNPs->NumReads - 1) / m_pChromSNPs->NumReads));
				if((Rslt=OutputSNPs())!=eBSFSuccess)
					{
					Reset(false);
					return(Rslt);
					}
				}

			PrevTargEntry = pSeg->ChromID;
			ChromLen = m_pSfxArray->GetSeqLen(PrevTargEntry);
			if(m_pChromSNPs == nullptr || (m_pChromSNPs != nullptr && (ChromLen + 16) > m_pChromSNPs->AllocChromLen))
				{
				if(m_pChromSNPs != nullptr)
					{
					delete []m_pChromSNPs;
					m_pChromSNPs = nullptr;
					}
				size_t AllocSize = sizeof(tsChromSNPs) + ((ChromLen + 16) * sizeof(tsSNPcnts));
				if((m_pChromSNPs = (tsChromSNPs *)new uint8_t[AllocSize])==nullptr)
					{
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessSNPs: Memory allocation of %zd bytes - %s",(int64_t)AllocSize,strerror(errno));
					Reset(false);
					return(eBSFerrMem);
					}
				m_pChromSNPs->AllocChromLen = ChromLen + 16;
				}
			memset(&m_pChromSNPs->Cnts,0,((ChromLen + 16) * sizeof(tsSNPcnts)));
			m_pChromSNPs->ChromLen = (uint32_t)ChromLen;
			m_pChromSNPs->ChromID = pSeg->ChromID;
			m_pChromSNPs->TotMatch = 0;
			m_pChromSNPs->TotMismatch = 0;
			m_pChromSNPs->MeanReadLen = 0;
			m_pChromSNPs->NumReads = 0;
			m_pChromSNPs->TotReadLen = 0;
			m_pChromSNPs->AdjacentSNPs[0].StartLoci = 0;
			m_pChromSNPs->AdjacentSNPs[0].EndLoci = 0;
			m_pChromSNPs->AdjacentSNPs[0].pFirstIterReadHit = nullptr;
			m_pChromSNPs->AdjacentSNPs[0].pPrevIterReadHit = 0;
			m_pChromSNPs->AdjacentSNPs[1].StartLoci = 0;
			m_pChromSNPs->AdjacentSNPs[1].EndLoci = 0;
			m_pChromSNPs->AdjacentSNPs[1].pFirstIterReadHit = nullptr;
			m_pChromSNPs->AdjacentSNPs[1].pPrevIterReadHit = 0;
			m_pChromSNPs->pFirstReadHit = nullptr;
			m_pChromSNPs->pLastReadHit = nullptr;
			PrevTargEntry = m_pChromSNPs->ChromID;
			PrevMMChromID = 0;
			PrevMMLoci = -1;
			}

		// get target genome sequence
		if(m_bIsSOLiD)
			{
			MatchLen = AdjHitLen(pSeg);
			HitLoci = AdjStartLoci(pSeg);
			HitLoci += 1;
			MatchLen -= 1;
			m_pSfxArray->GetColorspaceSeq(pSeg->ChromID,
									HitLoci,
									AssembSeq,MatchLen);	// get colorspace sequence


			}
		else
			{
			// get target assembly sequence for entry starting at offset and of length len
			MatchLen = AdjHitLen(pSeg);
			HitLoci = AdjStartLoci(pSeg);
			uint32_t RetSeqLen;
			if((RetSeqLen = m_pSfxArray->GetSeq(pSeg->ChromID,HitLoci,AssembSeq,MatchLen)) != MatchLen)
				continue;
			pAssembSeq = AssembSeq;
			for(SeqIdx = 0; SeqIdx < MatchLen; SeqIdx++,pAssembSeq++)
				*pAssembSeq = *pAssembSeq & 0x07;
			}

			// get accepted aligned read sequence
		pSeqVal = &pReadHit->Read[pReadHit->DescrLen+1];
		pSeqVal += pSeg->ReadOfs + pSeg->TrimLeft;
		pReadSeq = ReadSeq;

		if(m_bIsSOLiD)
			{
			// convert read sequence into colorspace
			uint8_t PrvBase = *pSeqVal & 0x07;
			for(SeqIdx = 1; SeqIdx <= MatchLen; SeqIdx++,pReadSeq++,pSeqVal++)
				{
				*pReadSeq = SOLiDmap[PrvBase][pSeqVal[1] & 0x07];
				PrvBase = pSeqVal[1] & 0x07;
				}
			// reverse, not complement, sequence if hit was onto '-' strand
			if(pSeg->Strand == '-')
				CSeqTrans::ReverseSeq(MatchLen,ReadSeq);
			}
		else
			{
			for(SeqIdx = 0; SeqIdx < MatchLen; SeqIdx++,pReadSeq++,pSeqVal++)
				*pReadSeq = *pSeqVal & 0x07;
			if(pSeg->Strand == '-')
				CSeqTrans::ReverseComplement(MatchLen,ReadSeq);
			}

		// double check not about to update snp counts past the expected chrom length
		if((HitLoci + MatchLen) > ChromLen)
			{
			if((MatchLen = (int)ChromLen - HitLoci) < 10)
				continue;
			}

		if(m_pChromSNPs->pFirstReadHit == nullptr)
			m_pChromSNPs->pFirstReadHit = pReadHit;
		m_pChromSNPs->pLastReadHit = pReadHit;
		m_pChromSNPs->TotReadLen += MatchLen;
		m_pChromSNPs->NumReads += 1;

		// now iterate read bases and if mismatch then update appropriate counts
		pSNP = &m_pChromSNPs->Cnts[HitLoci];
		pAssembSeq = &AssembSeq[0];
		pReadSeq = &ReadSeq[0];
		uint32_t Loci = HitLoci;
		bool bPairMM = false;
		int SeqMM = 0;
		for(SeqIdx = 0; SeqIdx < MatchLen; SeqIdx++, Loci++,pReadSeq++, pAssembSeq++,pSNP++)
			{
			if(*pAssembSeq >= eBaseN || (m_bIsSOLiD && *pReadSeq >= eBaseN) || *pReadSeq > eBaseN)
				{
				SeqMM += 1;
				continue;
				}

			if(m_bIsSOLiD)		// in colorspace, unpaired mismatches assumed to be sequencer errors and simply sloughed when identifying SNPs
				{
				if(Loci == 0)	// too problematic with SNPs in colorspace at the start of the target sequence, simply slough
					continue;

				if(!bPairMM && *pAssembSeq != *pReadSeq)
					{
					if(SeqIdx < (1+MatchLen))
						{
						if(pReadSeq[1] == pAssembSeq[1])
							{
							pSNP->NumRefBases += 1;
							m_pChromSNPs->TotMatch += 1;
							SeqMM += 1;
							continue;
							}
						}

					// get the previous target sequence base and use this + read colorspace space to derive the mismatch in basespace
					if(pSeg->ChromID != PrevMMChromID || Loci != PrevMMLoci)
						{
						PrevMMChromID = pSeg->ChromID;
						PrevMMLoci = Loci;
						m_pSfxArray->GetSeq(pSeg->ChromID,Loci-1,&TargBases[0],2);
						if(TargBases[0] > eBaseN)
							TargBases[0] = eBaseN;
						if(TargBases[1] > eBaseN)
							TargBases[1] = eBaseN;
						pSNP->RefBase = TargBases[1];
						}

					ReadBase = *pReadSeq;
					if(ReadBase > eBaseT)
						ReadBase = eBaseN;

					if(SeqMM == 0)
						ReadBase = SOLiDmap[TargBases[0]][ReadBase];
					else
						ReadBase = eBaseN;

					// sometimes it seems that a colorspace read may have had a sequencing error earlier in the read
					// or some mismatch such that the current loci mismatches in colorspace but matches in basespace
					// these strange bases are treated as though they are undefined and accrue counts as being eBaseN's
					if(ReadBase == pSNP->RefBase)
						ReadBase = eBaseN;
					pSNP->NonRefBaseCnts[ReadBase] += 1;
					pSNP->NumNonRefBases += 1;
					m_pChromSNPs->TotMismatch += 1;
					bPairMM = true;
					SeqMM += 1;
					}
				else
					{
					pSNP->NumRefBases += 1;
					m_pChromSNPs->TotMatch += 1;
					bPairMM = false;
					SeqMM = 0;
					}
				}
			else				// in basespace any mismatch is counted as a NonRefCnt
				{
				ReadBase = *pReadSeq & 0x07;
				TargBases[0] = *pAssembSeq & 0x07;

				pSNP->RefBase = TargBases[0];
				if(TargBases[0] == ReadBase)
					{
					pSNP->NumRefBases += 1;
					m_pChromSNPs->TotMatch += 1;
					}
				else
					{
					if(ReadBase > eBaseT)
						ReadBase = eBaseN;
					pSNP->NonRefBaseCnts[ReadBase] += 1;
					pSNP->NumNonRefBases += 1;
					m_pChromSNPs->TotMismatch += 1;
					}
				}
			}
		}
	}

if(m_pChromSNPs != nullptr)
	{
	m_pChromSNPs->MeanReadLen = (uint32_t)(((m_pChromSNPs->TotReadLen + m_pChromSNPs->NumReads - 1) / m_pChromSNPs->NumReads));
	if((Rslt=OutputSNPs())!=eBSFSuccess)
		{
		Reset(false);
		return(Rslt);
		}
	}

if(!m_bPackedBaseAlleles)
	{
	if(m_hSNPfile != -1)
		{
	#ifdef _WIN32
		_commit(m_hSNPfile);
	#else
		fsync(m_hSNPfile);
	#endif
		close(m_hSNPfile);
		m_hSNPfile = -1;
		}

	if(m_hDiSNPfile != -1)
		{
	#ifdef _WIN32
		_commit(m_hDiSNPfile);
	#else
		fsync(m_hDiSNPfile);
	#endif
		close(m_hDiSNPfile);
		m_hDiSNPfile = -1;
		}

	if(m_hTriSNPfile != -1)
		{
	#ifdef _WIN32
		_commit(m_hTriSNPfile);
	#else
		fsync(m_hTriSNPfile);
	#endif
		close(m_hTriSNPfile);
		m_hTriSNPfile = -1;
		}

	if(m_hWIGSpansFile != -1)
		{
	#ifdef _WIN32
		_commit(m_hWIGSpansFile);
	#else
		fsync(m_hWIGSpansFile);
	#endif
		close(m_hWIGSpansFile);
		m_hWIGSpansFile = -1;
		}
	

if(m_hSNPCentsfile != -1)
	{
	// report on the SNP centroid distributions
	int SNPCentroidIdx;
	int CentroidSeq;
	uint8_t Bases[cSNPCentroidLen];
	int BaseIdx;
	tsSNPCentroid *pCentroid;
	char szCentroids[4096];
	int BuffIdx;

	BuffIdx = sprintf(szCentroids,"\"CentroidID\",\"Seq\",\"NumInsts\",\"NumSNPs\",\"RefBase\",\"RefBaseCnt\",\"BaseA\",\"BaseC\",\"BaseG\",\"BaseT\",\"BaseN\"\n");
	pCentroid = m_pSNPCentroids;
	for(SNPCentroidIdx = 0; SNPCentroidIdx < cSNPCentroidEls; SNPCentroidIdx++, pCentroid++)
		{
		CentroidSeq = SNPCentroidIdx;
		for(BaseIdx = cSNPCentroidLen-1; BaseIdx >= 0; BaseIdx--)
			{
			Bases[BaseIdx] = CentroidSeq & 0x03;
			CentroidSeq >>= 2;
			}

		BuffIdx += sprintf(&szCentroids[BuffIdx],"%d,\"%s\",%d,%d,\"%c\",%d,%d,%d,%d,%d,%d\n",
							SNPCentroidIdx+1,CSeqTrans::MapSeq2Ascii(Bases,cSNPCentroidLen),pCentroid->NumInsts,pCentroid->NumSNPs,CSeqTrans::MapBase2Ascii(Bases[cSNPCentfFlankLen]),
							pCentroid->RefBaseCnt,pCentroid->NonRefBaseCnts[0],pCentroid->NonRefBaseCnts[1],pCentroid->NonRefBaseCnts[2],pCentroid->NonRefBaseCnts[3],pCentroid->NonRefBaseCnts[4]);


		if(BuffIdx + 200 > sizeof(szCentroids))
			{
			CUtility::RetryWrites(m_hSNPCentsfile,szCentroids,BuffIdx);
			BuffIdx = 0;
			}
		}
	if(BuffIdx)
		CUtility::RetryWrites(m_hSNPCentsfile,szCentroids,BuffIdx);
	}

	if(m_hSNPCentsfile != -1)
		{
	#ifdef _WIN32
		_commit(m_hSNPCentsfile);
	#else
		fsync(m_hSNPCentsfile);
	#endif
		close(m_hSNPCentsfile);
		m_hSNPCentsfile = -1;
		}

	if(m_hMarkerFile != -1)
		{
	#ifdef _WIN32
		_commit(m_hMarkerFile);
	#else
		fsync(m_hMarkerFile);
	#endif
		close(m_hMarkerFile);
		m_hMarkerFile = -1;
		}
	}
else
	{
	if(m_hPackedBaseAllelesFile != -1)
		{
#ifdef _WIN32
		_commit(m_hPackedBaseAllelesFile);
#else
		fsync(m_hPackedBaseAllelesFile);
#endif
		close(m_hPackedBaseAllelesFile);
		m_hPackedBaseAllelesFile = -1;
		}
	}


if(m_pChromSNPs != nullptr)
	{
	delete []m_pChromSNPs;
	m_pChromSNPs = nullptr;
	}
return(eBSFSuccess);
}

//--- the following function ProcessSiteProbabilites() is targeted for use in RNA-seq processing
// which should help in identifying differentially expressed transcripts
int
CKAligner::ProcessSiteProbabilites(int RelSiteStartOfs)	// offset the site octamer by this relative start offset (read start base == 0)
{
int LineLen;

int SeqIdx;
int SiteIdx;
etSeqBase AssembSeq[cMaxFastQSeqLen+1];	// to hold genome assembly sequence for current read loci

tsOctSitePrefs *pSitePrefs;
tsOctSitePrefs *pSitePref;
etSeqBase *pAssembSeq;
tsReadHit *pReadHit;
tsSegLoci *pSeg;
tBSFEntryID PrevTargEntry;
uint32_t HitLoci;
uint32_t PrevLoci;
uint32_t MatchLen;
uint32_t CurChromLen;

int TotOccs;
int TotSites;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing for alignment site probabilities...");

TotOccs = 0;
TotSites = 0;
memset(m_OctSitePrefs,0,sizeof(m_OctSitePrefs));
pSitePrefs = m_OctSitePrefs[0];
pSitePref = m_OctSitePrefs[1];
for(SiteIdx = 0; SiteIdx < cNumOctamers; SiteIdx++,pSitePref++,pSitePrefs++)
	{
	pSitePref->Octamer = SiteIdx;
	pSitePrefs->Octamer = SiteIdx;
	}
pReadHit = nullptr;
LineLen = 0;
PrevTargEntry = 0;
PrevLoci = -1;
CurChromLen = -1;
// iterate all accepted as aligned reads
// reads are assumed to have been sorted by chrom and start loci
while((pReadHit = IterSortedReads(pReadHit))!=nullptr)
	{
	if(pReadHit->NAR == eNARAccepted)
		{
		if(pReadHit->HitLoci.Hit.FlgInDel || pReadHit->HitLoci.Hit.FlgSplice)
			continue;

		pSeg = &pReadHit->HitLoci.Hit.Seg[0];
		if(pSeg->ChromID != (uint32_t)PrevTargEntry)
			{
			PrevTargEntry = pSeg->ChromID;
			PrevLoci = -1;
			CurChromLen = m_pSfxArray->GetSeqLen(pSeg->ChromID);
			}

		// get target assembly sequence for entry starting at MatchLoci, offset by RelSiteStartOfs, and of octamer length
		MatchLen = (uint32_t)pSeg->MatchLen;
		HitLoci = (uint32_t)pSeg->MatchLoci;

		if(pSeg->Strand == '+')
			HitLoci += RelSiteStartOfs;
		else
			{
			HitLoci += pSeg->MatchLen - 1;
			HitLoci -= RelSiteStartOfs;
			HitLoci -= 7;
			}

		if(HitLoci < 0)						// force octamer to start/end within the chrom sequence
			HitLoci = 0;
		else
			if((HitLoci + 8) >= CurChromLen)
				HitLoci = CurChromLen - 9;

		m_pSfxArray->GetSeq(pSeg->ChromID,HitLoci,AssembSeq,8);
		pAssembSeq = AssembSeq;
		for(SeqIdx = 0; SeqIdx < 8; SeqIdx++,pAssembSeq++)
			*pAssembSeq = *pAssembSeq & 0x07;

		if(pSeg->Strand == '-')
			{
			pSitePrefs = m_OctSitePrefs[1];
			CSeqTrans::ReverseComplement(8,AssembSeq);
			}
		else
			pSitePrefs = m_OctSitePrefs[0];
		pAssembSeq = &AssembSeq[0];

		// generate site index
		SiteIdx = 0;
		for(SeqIdx = 0; SeqIdx < 8; SeqIdx++)
			{
			if(*pAssembSeq > eBaseT)
				break;
			SiteIdx <<= 2;
			SiteIdx |= *pAssembSeq++;
			}
		if(SeqIdx != 8)
			continue;
		pReadHit->SiteIdx = SiteIdx;
		pSitePref = &pSitePrefs[SiteIdx];

		pSitePref->NumOccs += 1;
		TotOccs += 1;
		if(HitLoci != PrevLoci)
			{
			pSitePref->NumSites += 1;
			PrevLoci = HitLoci;
			TotSites += 1;
			}
		}
	}

// now to generate the relative abundance scores
// firstly determine the mean number of read alignment occurrences for each octamer
// find the top 64 (~ 0.1%) with highest mean occurrences
// all octamers are then scaled to this mean of the top 0.1%
pSitePrefs = m_OctSitePrefs[0];
pSitePref = m_OctSitePrefs[1];
for(SiteIdx = 0; SiteIdx < cNumOctamers; SiteIdx++,pSitePrefs++,pSitePref++)
	{
	if(pSitePrefs->NumSites >= 1)
		pSitePrefs->RelScale = (double)pSitePrefs->NumOccs/pSitePrefs->NumSites;
	else
		pSitePrefs->RelScale = 0.0;
	if(pSitePref->NumSites >= 1)
		pSitePref->RelScale = (double)pSitePref->NumOccs/pSitePref->NumSites;
	else
		pSitePref->RelScale = 0;
	}

// now sort ascending  by RelScale so top 0.1% can be determined and their mean determined
m_mtqsort.qsort(m_OctSitePrefs[0],cNumOctamers,sizeof(tsOctSitePrefs),SortSiteRelScale);
m_mtqsort.qsort(m_OctSitePrefs[1],cNumOctamers,sizeof(tsOctSitePrefs),SortSiteRelScale);

double TopWatsonMean = 0.0;
double TopCrickMean = 0.0;
pSitePrefs = &m_OctSitePrefs[0][0x0ffc0];
pSitePref = &m_OctSitePrefs[1][0x0ffc0];
for(SiteIdx = 0x0ffc0; SiteIdx < cNumOctamers; SiteIdx++,pSitePrefs++,pSitePref++)
	{
	TopWatsonMean += pSitePrefs->RelScale;
	pSitePrefs->RelScale = 1.0;
	TopCrickMean +=	pSitePref->RelScale;
	pSitePref->RelScale = 1.0;
	}
TopWatsonMean /= 64;
TopCrickMean /= 64;

// top means known now set normalisation scale factors
pSitePrefs = m_OctSitePrefs[0];
pSitePref = m_OctSitePrefs[1];
for(SiteIdx = 0; SiteIdx < 0x0ffc0; SiteIdx++,pSitePrefs++,pSitePref++)
	{
	if(pSitePrefs->RelScale > 0.0)
		pSitePrefs->RelScale = max(0.0001,pSitePrefs->RelScale / TopWatsonMean);
	if(pSitePref->RelScale)
		pSitePref->RelScale = max(0.0001,pSitePref->RelScale / TopCrickMean);
	}

// restore m_OctSitePrefs to original octamer ascending order
m_mtqsort.qsort(m_OctSitePrefs[0],cNumOctamers,sizeof(tsOctSitePrefs),SortSiteRelOctamer);
m_mtqsort.qsort(m_OctSitePrefs[1],cNumOctamers,sizeof(tsOctSitePrefs),SortSiteRelOctamer);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Completed alignment site probabilities");
return(eBSFSuccess);
}



char *
CKAligner::Octamer2Txt(int Octamer) // Report on site octamer site preferencing distribution
{
static char szOctamer[0x09];	// to contain '\0' terminated octamer as nucleotide bases 'a'..'t'
char *pChr;
int Idx;
pChr = &szOctamer[8];
*pChr-- = '\0';
for(Idx = 0; Idx < 8; Idx++,pChr--)
	{
	switch(Octamer & 0x03) {
		case 0:
			*pChr = 'a';
			break;
		case 1:
			*pChr = 'c';
			break;
		case 2:
			*pChr = 'g';
			break;
		case 3:
			*pChr = 't';
			break;
		}
	Octamer >>= 2;
	}
return(szOctamer);
return(nullptr);
}

int
CKAligner::WriteSitePrefs(void)
{
char szBuff[0x03fff];
int BuffIdx;
int SiteIdx;
tsOctSitePrefs *pSite;
BuffIdx = sprintf(szBuff,"\"Id\",\"Strand\",\"Octamer\",\"TotalHits\",\"UniqueLoci\",\"RelScale\"\n");
pSite = m_OctSitePrefs[0];
for(SiteIdx = 0; SiteIdx < 0x0ffff; SiteIdx++,pSite++)
	{
	BuffIdx += sprintf(&szBuff[BuffIdx],"%d,\"+\",\"%s\",%d,%d,%1.3f\n",SiteIdx+1,Octamer2Txt(pSite->Octamer),pSite->NumOccs,pSite->NumSites,pSite->RelScale);
	if(BuffIdx + 200 > sizeof(szBuff))
		{
		CUtility::RetryWrites(m_hSitePrefsFile,szBuff,BuffIdx);
		BuffIdx = 0;
		}
	}
pSite = m_OctSitePrefs[1];
for(SiteIdx = 0; SiteIdx < 0x0ffff; SiteIdx++,pSite++)
	{
	BuffIdx += sprintf(&szBuff[BuffIdx],"%d,\"-\",\"%s\",%d,%d,%1.3f\n",SiteIdx+1,Octamer2Txt(pSite->Octamer),pSite->NumOccs,pSite->NumSites,pSite->RelScale);
	if(BuffIdx + 200 > sizeof(szBuff))
		{
		CUtility::RetryWrites(m_hSitePrefsFile,szBuff,BuffIdx);
		BuffIdx = 0;
		}
	}

if(BuffIdx > 0)
	{
	CUtility::RetryWrites(m_hSitePrefsFile,szBuff,BuffIdx);
	BuffIdx = 0;
	}
return(eBSFSuccess);
}

int
CKAligner::LoadReads(char *pszRdsFile)	// file containing preprocessed reads (genreads output)
{
int Rslt;
size_t memreq;
int RdLen;
tsRawReadV5 *pReadV5;					// current preprocessed read being processed if V5
tsRawReadV6 *pReadV6;					// current preprocessed read being processed if V6

uint8_t *pReadBuff;						// alloc'd to buffer the preprocessed reads from disk
tsReadHit *pReadHit;					// current read hit
int BuffLen;
int BuffOfs;
char *pChr;
char Chr;
// check if file name contains any wildcard chars
pChr = pszRdsFile;
while(Chr = *pChr++)
	if(Chr == '*' || Chr == '?' || Chr == '[' || Chr == ']')
		return(eBSFerrNotBioseq);

#ifdef _WIN32
m_hInFile = open(pszRdsFile, O_READSEQ ); // file access is normally sequential..
#else
m_hInFile = open64(pszRdsFile, O_READSEQ ); // file access is normally sequential..
#endif

if(m_hInFile == -1)					// check if file open succeeded
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open %s - %s",pszRdsFile,strerror(errno));
	return(eBSFerrOpnFile);
	}

// expecting a preprocessed .rds file as input, header processing will confirm!
if((Rslt=Disk2Hdr(pszRdsFile))!=eBSFSuccess)
	{
	close(m_hInFile);
	m_hInFile = -1;
	return((teBSFrsltCodes)Rslt);
	}

if((pReadBuff = new uint8_t [cRdsBuffAlloc])==nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Memory allocation of %d bytes - %s",cRdsBuffAlloc,strerror(errno));
	close(m_hInFile);
	m_hInFile = -1;
	return(eBSFerrMem);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Reads file '%s' generator version: %d",pszRdsFile,m_FileHdr.Version);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Contains total of %d reads from %d source reads with duplicate sequences %s, PMode was %d, Quality was %d, %d bases 5' trimmed, %d bases 3' trimmed",
			m_FileHdr.NumRds,m_FileHdr.OrigNumReads,m_FileHdr.FlagsK ? "retained":"removed", m_FileHdr.PMode,m_FileHdr.QMode,m_FileHdr.Trim5,m_FileHdr.Trim3);
if(m_FileHdr.FlagsPR)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Reads were processed as %s",m_FileHdr.FlagsPR ? "paired end reads" : "single end reads");

gDiagnostics.DiagOutMsgOnly(eDLInfo,"Reads were processed from %d files",m_FileHdr.NumFiles);
char *pszSrcFile = (char *)m_FileHdr.FileNames;
for(BuffOfs=0;BuffOfs<m_FileHdr.NumFiles;BuffOfs++)
	{
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Source file: '%s'",pszSrcFile);
	pszSrcFile += strlen(pszSrcFile) + 1;
	}

// ensure there is at least one read...
if(m_FileHdr.NumRds == 0 || m_FileHdr.TotReadsLen == 0)
	{
	gDiagnostics.DiagOut(eDLWarn,gszProcName,"Nothing to do, '%s' contains no reads...",pszRdsFile);
	delete []pReadBuff;
	close(m_hInFile);
	m_hInFile = -1;
	return(eBSFerrOpnFile);
	}

// initial allocation of memory to hold all pre-processed reads plus a little safety margin (10000) bytes)
memreq = ((size_t)m_FileHdr.NumRds * sizeof(tsReadHit)) + (size_t)m_FileHdr.TotReadsLen + 10000;
AcquireExclusiveLock();
#ifdef _WIN32
m_pReadHits = (tsReadHit *) malloc(memreq);	// initial and perhaps the only allocation

if(m_pReadHits == nullptr)
	{
	ReleaseExclusiveLock();
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Memory allocation of %zd bytes - %s",(int64_t)memreq,strerror(errno));
	delete []pReadBuff;
	close(m_hInFile);
	m_hInFile = -1;
	return(eBSFerrMem);
	}
#else
// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
m_pReadHits = (tsReadHit *)mmap(nullptr,memreq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
if(m_pReadHits == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Memory allocation of %zd bytes through mmap()  failed - %s",(int64_t)memreq,strerror(errno));
	m_pReadHits = nullptr;
	ReleaseExclusiveLock();
	delete pReadBuff;
	close(m_hInFile);
	m_hInFile = -1;
	return(eBSFerrMem);
	}
#endif

m_AllocdReadHitsMem = memreq;
m_UsedReadHitsMem = 0;
m_FinalReadID = 0;
m_NumReadsLoaded = 0;
m_NumDescrReads = 0;
ReleaseExclusiveLock();

// iterate each read sequence starting from the first
lseek(m_hInFile,(long)m_FileHdr.RdsOfs,SEEK_SET);
BuffLen = 0;
BuffOfs = 0;

int SizeOfRawRead = m_FileHdr.Version == 5 ? sizeof(tsRawReadV5) : sizeof(tsRawReadV6);
int CurReadLen;
int CurDescrLen;

while((RdLen = read(m_hInFile,&pReadBuff[BuffLen],cRdsBuffAlloc - BuffLen)) > 0)
	{
	BuffLen += RdLen;
	BuffOfs = 0;
	while((BuffLen - BuffOfs) >=  SizeOfRawRead)
		{
		if(m_FileHdr.Version == 5)
			{
			pReadV5 = (tsRawReadV5 *)&pReadBuff[BuffOfs];
			CurDescrLen = (int)pReadV5->DescrLen;
			CurReadLen = (int)pReadV5->ReadLen;
			}
		else
			{
			pReadV6 = (tsRawReadV6 *)&pReadBuff[BuffOfs];
			CurDescrLen = (int)pReadV6->DescrLen;
			CurReadLen = (int)pReadV6->ReadLen;
			}

		if((int)(CurDescrLen + CurReadLen + SizeOfRawRead) > (BuffLen - BuffOfs))
			break;
		BuffOfs += CurDescrLen + CurReadLen + SizeOfRawRead;
			// pReadV5/V6 now pts at a read
			// shouldn't but is there a need to allocate more memory?
		if(m_UsedReadHitsMem + (sizeof(tsReadHit) + CurReadLen + CurDescrLen) >= (m_AllocdReadHitsMem - 5000))
			{
			AcquireExclusiveLock();
			memreq = m_AllocdReadHitsMem + ((sizeof(tsReadHit) + (size_t)cDfltReadLen) * cReadsHitReAlloc);
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadReads: Needing memory re-allocation to %zd bytes from %zd",(int64_t)m_AllocdReadHitsMem,(int64_t)memreq);

#ifdef _WIN32
			pReadHit = (tsReadHit *) realloc(m_pReadHits,memreq);
#else
			pReadHit = (tsReadHit *)mremap(m_pReadHits,m_AllocdReadHitsMem,memreq,MREMAP_MAYMOVE);
			if(pReadHit == MAP_FAILED)
				pReadHit = nullptr;
#endif
			if(pReadHit == nullptr)
				{
				ReleaseExclusiveLock();
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Memory re-allocation to %zd bytes - %s",(int64_t)(m_AllocdReadHitsMem + ((sizeof(tsReadHit) + cDfltReadLen)*cReadsHitReAlloc)),strerror(errno));
				delete []pReadBuff;
				close(m_hInFile);
				m_hInFile = -1;
				return(eBSFerrMem);
				}
			m_pReadHits = pReadHit;
			m_AllocdReadHitsMem = memreq;
			ReleaseExclusiveLock();
			}

		pReadHit = (tsReadHit *)((uint8_t *)m_pReadHits + m_UsedReadHitsMem);
		m_UsedReadHitsMem += sizeof(tsReadHit) + CurReadLen + CurDescrLen;
		memset(pReadHit,0,sizeof(tsReadHit));
		pReadHit->HitLoci.Hit.Seg[0].Strand = '?';
		pReadHit->ReadLen = CurReadLen;
		pReadHit->DescrLen = CurDescrLen;

		if(m_FileHdr.Version == 5)
			{
			pReadHit->ReadID = pReadV5->ReadID;
			pReadHit->PairReadID = pReadV5->PairReadID;
			pReadHit->NumReads = pReadV5->NumReads;

			memmove(pReadHit->Read,pReadV5->Read,CurDescrLen + 1 + CurReadLen);
			m_FinalReadID = pReadV5->ReadID;
			}
		else
			{
			pReadHit->ReadID = pReadV6->ReadID;
			pReadHit->PairReadID = pReadV6->PairReadID;
			pReadHit->NumReads = pReadV6->NumReads;
			memmove(pReadHit->Read,pReadV6->Read,CurDescrLen + 1 + CurReadLen);
			m_FinalReadID = pReadV6->ReadID;
			}

		m_NumDescrReads += 1;

		// processing threads are only updated with actual number of loaded reads every 100000 reads so as
		// to minimise disruption to the actual aligner threads which will also be serialised through m_hMtxIterReads
		if(m_NumDescrReads > 0 && !(m_NumDescrReads % 100000))
			{
			AcquireSerialise();
			m_FinalReadID = m_NumDescrReads;
			m_NumReadsLoaded = m_NumDescrReads;
			ReleaseSerialise();
			}
		}
	BuffLen -= BuffOfs;
	if(BuffLen)
		memmove(pReadBuff,&pReadBuff[BuffOfs],BuffLen);
	}
delete []pReadBuff;
close(m_hInFile);
m_hInFile = -1;
if(m_NumDescrReads != m_NumReadsLoaded)
	{
	AcquireSerialise();
	m_FinalReadID = m_NumDescrReads;
	m_NumReadsLoaded = m_NumDescrReads;
	ReleaseSerialise();
	}
return(m_NumDescrReads);
}

int
CKAligner::CreateMutexes(void)
{
if(m_bMutexesCreated)
	return(eBSFSuccess);

#ifdef _WIN32
InitializeSRWLock(&m_hRwLock);
#else
if(pthread_rwlock_init (&m_hRwLock,nullptr)!=0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to create rwlock");
	return(eBSFerrInternal);
	}
#endif

#ifdef _WIN32
if((m_hMtxIterReads = CreateMutex(nullptr,false,nullptr))==nullptr)
	{
#else
if(pthread_mutex_init (&m_hMtxIterReads,nullptr)!=0)
	{
	pthread_rwlock_destroy(&m_hRwLock);
#endif
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to create mutex");
	return(eBSFerrInternal);
	}

#ifdef _WIN32
if((m_hMtxMHReads = CreateMutex(nullptr,false,nullptr))==nullptr)
	{
#else
if(pthread_mutex_init (&m_hMtxMHReads,nullptr)!=0)
	{
#endif
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to create mutex");
#ifdef _WIN32
	CloseHandle(m_hMtxIterReads);
#else
	pthread_rwlock_destroy(&m_hRwLock);
	pthread_mutex_destroy(&m_hMtxIterReads);
#endif
	return(eBSFerrInternal);
	}
if(m_MLMode != eMLdefault)
	{
#ifdef _WIN32
	if((m_hMtxMultiMatches = CreateMutex(nullptr,false,nullptr))==nullptr)
		{
		CloseHandle(m_hMtxIterReads);
		CloseHandle(m_hMtxMHReads);
#else
	if(pthread_mutex_init (&m_hMtxMultiMatches,nullptr)!=0)
		{
		pthread_mutex_destroy(&m_hMtxIterReads);
		pthread_mutex_destroy(&m_hMtxMHReads);
		pthread_rwlock_destroy(&m_hRwLock);
#endif
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to create mutex");
		return(eBSFerrInternal);
		}
	}
m_bMutexesCreated = true;
return(eBSFSuccess);
}

void
CKAligner::DeleteMutexes(void)
{
if(!m_bMutexesCreated)
	return;
#ifdef _WIN32
CloseHandle(m_hMtxIterReads);
CloseHandle(m_hMtxMHReads);
if(m_MLMode != eMLdefault)
	CloseHandle(m_hMtxMultiMatches);
#else
pthread_mutex_destroy(&m_hMtxIterReads);
pthread_mutex_destroy(&m_hMtxMHReads);
pthread_rwlock_destroy(&m_hRwLock);
if(m_MLMode != eMLdefault)
	pthread_mutex_destroy(&m_hMtxMultiMatches);
#endif
m_bMutexesCreated = false;
}

#ifdef _WIN32
unsigned __stdcall KCoredApproxThread(void * pThreadPars)
#else
void *KCoredApproxThread(void * pThreadPars)
#endif
{
int Rslt;
tsThreadMatchPars *pPars = (tsThreadMatchPars *)pThreadPars; // makes it easier not having to deal with casts!
CKAligner *pKAligner = (CKAligner *)pPars->pThis;
Rslt = pKAligner->ProcCoredApprox(pPars);
pPars->Rslt = Rslt;
#ifdef _WIN32
_endthreadex(Rslt < 0 ? 1 : 0);
return(Rslt < 0 ? 1 : 0);
#else
pthread_exit(&pPars->Rslt);
#endif
}

// LocateCoredApprox
// Locates all cored approximates
int
CKAligner::LocateCoredApprox(int MinEditDist,	// any matches must have at least this edit distance to the next best match
							int MaxSubs)		// maximum number of substitutions allowed per 100bp of accepted aligned read length
{
int Rslt;
int CurBlockID;							// current suffix block being processed
tBSFEntryID CurChromID;				    // current suffix array entry being processed
uint32_t TotNumReadsProc;					// total number of reads processed
uint32_t PlusHits;
uint32_t MinusHits;
uint32_t ChimericHits;

uint32_t CurReadsProcessed;
uint32_t PrevReadsProcessed;
uint32_t CurReadsLoaded;
uint32_t PrevReadsLoaded;
int MaxNumSlides;

int ThreadIdx;
tsThreadMatchPars *pWorkerThreads;

if((pWorkerThreads = new tsThreadMatchPars [m_NumThreads+1])==nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Fatal: unable to allocate memory for %d tsThreadMatchPars", m_NumThreads);
	Reset(false);
	return(eBSFerrMem);
	}

m_PerThreadAllocdIdentNodes = cMaxNumIdentNodes;
m_TotAllocdIdentNodes = m_PerThreadAllocdIdentNodes * m_NumThreads;
if((m_pAllocsIdentNodes = new tsIdentNode [m_TotAllocdIdentNodes])==nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to allocate memory for %d tsIdentNodes",m_TotAllocdIdentNodes);
	delete []pWorkerThreads;
	Reset(false);
	return(eBSFerrMem);
	}

if(m_PEproc != ePEdefault)
	m_MaxMLPEmatches = cMaxMLPEmatches;
else
	m_MaxMLPEmatches = 0;

if((m_pAllocsMultiHitLoci = new tsHitLoci [m_NumThreads * (max(m_MaxMLPEmatches,m_MaxMLmatches) + cPriorityExacts)])==nullptr)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to allocate memory for %d tsHitLoci",m_NumThreads * (m_MaxMLmatches + cPriorityExacts));
	delete []pWorkerThreads;
	Reset(false);
	return(eBSFerrMem);
	}

if(m_MLMode == eMLall)
	{
	if((m_pAllocsMultiHitBuff = new uint8_t [m_NumThreads * cReadHitBuffLen])==nullptr)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to allocate memory for %d multihit record buffering",m_NumThreads * cReadHitBuffLen);
		delete []pWorkerThreads;
		Reset(false);
		return(eBSFerrMem);
		}
	}
else
	m_pAllocsMultiHitBuff = nullptr;

// load single SfxBlock, expected to contain all chromosomes, and process all reads against that block
PlusHits = 0;
MinusHits = 0;
ChimericHits = 0;
TotNumReadsProc = 0;
PrevReadsProcessed = 0;
PrevReadsLoaded = 0;

CurChromID = 0;
CurBlockID = 1;
if((Rslt=m_pSfxArray->SetTargBlock(CurBlockID))<0)
	{
	while(m_pSfxArray->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pSfxArray->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to load genome assembly suffix array");
	delete []pWorkerThreads;
	return(Rslt);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Genome assembly suffix array loaded");

	// determine minimum core length from targeted sequence length
	// core length is a balance between sensitivity and throughput
	// reducing core size has a relatively minor effect on sensitivity but significantly reduces throughput
	// large genomes require larger cores, more sensitive alignments require smaller cores
m_BlockTotSeqLen = m_pSfxArray->GetTotSeqsLen();

int AutoCoreLen = 1;
while (m_BlockTotSeqLen >>= 2)
	AutoCoreLen++;
AutoCoreLen -= 1;
m_MinCoreLen = max(cKAMinCoreLen, AutoCoreLen);

// MaxNumSlides is per 100bp of read length
// more slides enables higher sensitivity but negatively impacts on alignment throughput
switch(m_PMode) {
	case ePMUltraSens:				// ultra sensitive - much slower
		m_MinCoreLen -= 2;
		MaxNumSlides = 9;			
		break;
	case ePMMoreSens:				// more sensitive - slower
		m_MinCoreLen -= 1;
		MaxNumSlides = 8;
		break;
	case ePMdefault:				// default processing mode
		MaxNumSlides = 8;
		break;
	default:						// less sensitive but quicker
		m_MinCoreLen += 2;
		MaxNumSlides = 6;
		break;
	}
m_MaxNumSlides = MaxNumSlides;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Now aligning with minimum core size of %dbp...\n",m_MinCoreLen);
m_pSfxArray->InitialiseCoreKMers(m_MinCoreLen);

m_ThreadCoredApproxRslt = 0;
ResetThreadedIterReads();
memset(pWorkerThreads,0,sizeof(tsThreadMatchPars) * m_NumThreads);
tsThreadMatchPars *pWorkerThread;
pWorkerThread = pWorkerThreads;
for(ThreadIdx = 0; ThreadIdx < m_NumThreads; ThreadIdx++, pWorkerThread++)
	{
	pWorkerThread->ThreadIdx = ThreadIdx + 1;
	pWorkerThread->pThis = this;
	pWorkerThread->NumIdentNodes = m_PerThreadAllocdIdentNodes;
	pWorkerThread->pIdentNodes = &m_pAllocsIdentNodes[m_PerThreadAllocdIdentNodes * ThreadIdx];
	pWorkerThread->pMultiHits = &m_pAllocsMultiHitLoci[(max(m_MaxMLPEmatches, m_MaxMLmatches) + cPriorityExacts) * ThreadIdx];
	pWorkerThread->CurBlockID = CurBlockID;
	pWorkerThread->MinEditDist = MinEditDist;
	pWorkerThread->MaxSubs = MaxSubs;
	pWorkerThread->AlignStrand = m_AlignStrand;
	pWorkerThread->microInDelLen = m_microInDelLen;
	pWorkerThread->SpliceJunctLen = m_SpliceJunctLen;
	pWorkerThread->MinCoreLen = m_MinCoreLen;
	pWorkerThread->MaxNumSlides = MaxNumSlides;
	pWorkerThread->MinChimericLen = m_MinChimericLen;
	if(m_MLMode == eMLall)
		pWorkerThread->pszOutBuff = &m_pAllocsMultiHitBuff[cReadHitBuffLen * ThreadIdx];
	else
		pWorkerThread->pszOutBuff = nullptr;
	pWorkerThread->OutBuffIdx = 0;
#ifdef _WIN32
	pWorkerThread->threadHandle = (HANDLE)_beginthreadex(nullptr,0x0fffff,KCoredApproxThread, pWorkerThread,0,&pWorkerThread->threadID);
#else
	pWorkerThread->threadRslt =	pthread_create (&pWorkerThread->threadID , nullptr , KCoredApproxThread , pWorkerThread);
#endif
	}

// allow threads a few seconds to startup
#ifdef _WIN32
	Sleep(5000);
#else
	sleep(5);
#endif

uint32_t ReportProgressSecs;
ReportProgressSecs = 60;
if(m_SampleNthRawRead > 1)
	ReportProgressSecs = 30;

// let user know that kit4b is working hard...
ApproxNumReadsProcessed(&PrevReadsProcessed,&PrevReadsLoaded);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: %u reads aligned from %u loaded",PrevReadsProcessed,PrevReadsLoaded);

// wait for all threads to have completed
pWorkerThread = pWorkerThreads;
for(ThreadIdx = 0; ThreadIdx < m_NumThreads; ThreadIdx++, pWorkerThread++)
	{
#ifdef _WIN32
	while(WAIT_TIMEOUT == WaitForSingleObject(pWorkerThread->threadHandle, (DWORD)ReportProgressSecs * 1000))
		{
		ApproxNumReadsProcessed(&CurReadsProcessed,&CurReadsLoaded);
		if(CurReadsProcessed > PrevReadsProcessed || CurReadsLoaded > PrevReadsLoaded)
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: %u reads aligned from %u loaded",CurReadsProcessed,CurReadsLoaded);
		PrevReadsProcessed = CurReadsProcessed;
		PrevReadsLoaded = CurReadsLoaded;
		}
	CloseHandle(pWorkerThread->threadHandle);
#else
	struct timespec ts;
	int JoinRlt;
	clock_gettime(CLOCK_REALTIME, &ts);
	ts.tv_sec += ReportProgressSecs;
	while((JoinRlt = pthread_timedjoin_np(pWorkerThread->threadID, nullptr, &ts)) != 0)
		{
		ApproxNumReadsProcessed(&CurReadsProcessed,&CurReadsLoaded);
		if(CurReadsProcessed > PrevReadsProcessed || CurReadsLoaded > PrevReadsLoaded)
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: %u reads aligned from %u loaded",CurReadsProcessed,CurReadsLoaded);
		PrevReadsProcessed = CurReadsProcessed;
		PrevReadsLoaded = CurReadsLoaded;
		ts.tv_sec += ReportProgressSecs;
		}

#endif
	PlusHits += pWorkerThread->PlusHits;
	MinusHits += pWorkerThread->MinusHits;
	ChimericHits += pWorkerThread->ChimericHits;
	TotNumReadsProc += pWorkerThread->NumReadsProc;

	if(pWorkerThread->OutBuffIdx != 0)
		{
#ifdef _WIN32
		WaitForSingleObject(m_hMtxMultiMatches,INFINITE);
#else
		pthread_mutex_lock(&m_hMtxMultiMatches);
#endif
		if((cAllocLineBuffSize - m_szLineBuffIdx) < (int)(pWorkerThread->OutBuffIdx + ((cMaxFastQSeqLen * 2) + 1024)))
			{
			if(!m_bgzOutFile)
				CUtility::RetryWrites(m_hOutFile,m_pszLineBuff,m_szLineBuffIdx);
			else
				CUtility::SafeWrite_gz(m_gzOutFile,m_pszLineBuff,m_szLineBuffIdx);
			m_szLineBuffIdx = 0;
			}
		memmove(&m_pszLineBuff[m_szLineBuffIdx], pWorkerThread->pszOutBuff, pWorkerThread->OutBuffIdx);
		m_szLineBuffIdx += pWorkerThread->OutBuffIdx;
#ifdef _WIN32
		ReleaseMutex(m_hMtxMultiMatches);
#else
		pthread_mutex_unlock(&m_hMtxMultiMatches);
#endif
		pWorkerThread->OutBuffIdx = 0;
		}
	}


// pickup the read loader thread, if the reads processing threads all finished then the loader thread should also have finished
#ifdef _WIN32
if(m_hThreadLoadReads != nullptr)
	{
	while(WAIT_TIMEOUT == WaitForSingleObject(m_hThreadLoadReads, 5000))
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: waiting");
		}
	CloseHandle(m_hThreadLoadReads);
	m_hThreadLoadReads = nullptr;
	}
#else
if(m_ThreadLoadReadsID != 0)
	{
	struct timespec ts;
	int JoinRlt;
	clock_gettime(CLOCK_REALTIME, &ts);
	ts.tv_sec += 5;
	while((JoinRlt = pthread_timedjoin_np(m_ThreadLoadReadsID, nullptr, &ts)) != 0)
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: waiting");
		ts.tv_sec += 60;
		}
	m_ThreadLoadReadsID = 0;
	}
#endif

// Checking here that the reads were all loaded w/o any major dramas!
if(m_ThreadLoadReadsRslt < 0 || m_ThreadCoredApproxRslt < 0)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: Early terminated");
	Reset(false);
	return(m_ThreadLoadReadsRslt < 0 ? m_ThreadLoadReadsRslt : m_ThreadCoredApproxRslt);
	}
ApproxNumReadsProcessed(&CurReadsProcessed,&CurReadsLoaded);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Alignment of %u from %u loaded completed",CurReadsProcessed,CurReadsLoaded);

m_PerThreadAllocdIdentNodes = 0;
m_TotAllocdIdentNodes = 0;
if(m_pAllocsIdentNodes != nullptr)
	{
	delete []m_pAllocsIdentNodes;
	m_pAllocsIdentNodes = nullptr;
	}
if(m_pAllocsMultiHitLoci != nullptr)
	{
	delete []m_pAllocsMultiHitLoci;
	m_pAllocsMultiHitLoci = nullptr;
	}

if(m_pAllocsMultiHitBuff != nullptr)
	{
	delete []m_pAllocsMultiHitBuff;
	m_pAllocsMultiHitBuff = nullptr;
	}

if((m_FMode == eFMsam || m_FMode == eFMsamAll) && m_MLMode == eMLall)
	m_pszLineBuff[m_szLineBuffIdx] = '\0';

if(m_szLineBuffIdx > 0)
	{
	if(!m_bgzOutFile)
		CUtility::RetryWrites(m_hOutFile,m_pszLineBuff,m_szLineBuffIdx);
	else
		CUtility::SafeWrite_gz(m_gzOutFile,m_pszLineBuff,m_szLineBuffIdx);
	m_szLineBuffIdx = 0;
	}

delete []pWorkerThreads;
return(eBSFSuccess);
}

int
CKAligner::AlignRead(tsReadHit* pReadHit,
				tsThreadMatchPars* pPars)
{
char szPriorityChromName[100];
int PriorityChromID;

int SeqIdx;
int NumNs;
uint8_t* pSeqVal;
etSeqBase* pSeq;
int Rslt;
int HitRslt;

int ProbeLen;
int RandIdx;
int HitIdx;
tsHitLoci* pHit;
int MaxTotMM;
int CoreLen;
int CoreDelta;
int MaxNumSlides;
int MaxML = max(m_MaxMLPEmatches,m_MaxMLmatches);

// functionalise starts
pReadHit->NAR = eNARNoHit;			// assume unable to align read
pReadHit->FlgPEAligned = 0;
pReadHit->LowHitInstances = 0;
pReadHit->LowMMCnt = 0;
pReadHit->NumHits = 0;
pPars->NumReadsProc += 1;

// get sequence for read and remove any packed quality values
pSeqVal = &pReadHit->Read[pReadHit->DescrLen + 1];
pSeq = pPars->Sequence;
NumNs = 0;
int MaxNsSeq = 0;		// maximum allowed for this sequence
if (m_MaxNs)
	MaxNsSeq = max(((pReadHit->ReadLen * m_MaxNs) / 100), m_MaxNs);

for (SeqIdx = 0; SeqIdx < pReadHit->ReadLen; SeqIdx++, pSeq++, pSeqVal++)
	{
	if ((*pSeq = (*pSeqVal & 0x07)) > eBaseN)
		break;
	*pSeq &= ~cRptMskFlg;
	if (*pSeq == eBaseN)
		{
		if (++NumNs > MaxNsSeq)
			break;
		}
	}
if (SeqIdx != pReadHit->ReadLen) // if too many 'N's...
	{
	pPars->pPrevReadHit = nullptr;
	pReadHit->NAR = eNARNs;
	pPars->PrevMatchLen = 0;
	pPars->NumSloughedNs += 1;
	return(eHRSeqErrs);
	}

if (m_bIsSOLiD)	// if SOLiD colorspace then need to convert read back into colorspace before attempting to locate
	{
	pSeq = pPars->Sequence;
	uint8_t PrvBase = *pSeq;
	for (SeqIdx = 1; SeqIdx < pReadHit->ReadLen; SeqIdx++, pSeq++)
		{
		*pSeq = SOLiDmap[PrvBase][pSeq[1]];
		PrvBase = pSeq[1];
		}
	ProbeLen = pReadHit->ReadLen;
	pPars->MatchLen = ProbeLen - 1;
	}
else
	{
	ProbeLen = pReadHit->ReadLen;
	pPars->MatchLen = ProbeLen;
	}

	// note: MaxSubs is specified by user as being per 100bp of read length, e.g. if user specified '-s5' and a read is 200bp then 
	// 10 mismatches will be allowed for that specific read
MaxTotMM = pPars->MaxSubs == 0 ? 0 : max(1, (int)(0.5 + (pPars->MatchLen * pPars->MaxSubs) / 100.0));

if (MaxTotMM > cMaxTotAllowedSubs)		// irrespective of length allow at most this many subs
	MaxTotMM = cMaxTotAllowedSubs;

	// The window core length is set to be read length / (subs+1) for minimum Hamming difference of 1, and
	// to be read length / (subs+2) for minimum Hamming difference of 2
	// The window core length is clamped to be at least m_MinCoreLen
CoreLen = max(m_MinCoreLen, pPars->MatchLen / (pPars->MinEditDist == 1 ? MaxTotMM + 1 : MaxTotMM + 2));
MaxNumSlides = max(1, ((pPars->MaxNumSlides * ProbeLen) + 99) / 100);
CoreDelta = max(ProbeLen / MaxNumSlides - 1, CoreLen);

int LowHitInstances;
int LowMMCnt;
int NxtLowMMCnt;

LowHitInstances = pReadHit->LowHitInstances;
LowMMCnt = pReadHit->LowMMCnt;
NxtLowMMCnt = pReadHit->NxtLowMMCnt;

	// for priority alignments to known reference sequences (example would be cDNA transcripts assembled with a de Novo assembly) then
	// a) align for exact matches allowing say 10 multiloci hits
	// b) iterate these multiloci hits and discard any not aligning to a ref sequence
	// c) process those aligning to reference as if these were uniquely aligning
int RefExacts = cPriorityExacts;
bool bProcNorm;

if (m_pPriorityRegionBED != nullptr)
	RefExacts = cPriorityExacts;
else
	RefExacts = 0;

	// Heuristic is that if read sequence is identical to previously processed sequence then
	// simply reuse previously AlignReads() hit results - saves a lot of processing time
if (pPars->pPrevReadHit == nullptr || pPars->bForceNewAlignment || pPars->MatchLen != pPars->PrevMatchLen || memcmp(pPars->Sequence, pPars->PrevSequence, pPars->MatchLen))
	{
	memset(pPars->pMultiHits, 0, sizeof(tsHitLoci));
	bProcNorm = true;
	if (RefExacts > 0)
		{
		HitRslt = m_pSfxArray->AlignReads(0,						// flags indicating if lower levels need to do any form of extended processing with this specific read...
				pReadHit->ReadID,				// identifies this read
				pPars->MinChimericLen,			// minimum chimeric length as a percentage (0 to disable, otherwise 50..99) of probe sequence
				MaxTotMM,						// max number of mismatches allowed
				CoreLen,						// core window length
				CoreDelta,						// core window offset increment (1..n)
				MaxNumSlides,					// limit on number of times core window can be moved or slide to right over read per 100bp of read length
				pPars->MinCoreLen,				// minimum core length allowed
				pPars->MinEditDist,				// minimum (1..n) mismatch difference between the best and next best core alignment
				pPars->AlignStrand,				// watson, crick or both?
				0,								// microInDel length maximum
				0,								// maximum splice junction length when aligning RNAseq reads
				&LowHitInstances,				// Out number of match instances for lowest number of mismatches thus far for this read
				&LowMMCnt,						// Out lowest number of mismatches thus far for this read
				&NxtLowMMCnt,					// Out next to lowest number of mismatches thus far for this read
				pPars->Sequence,				// probe
				pPars->MatchLen,				// probe length
				MaxML + RefExacts,				// (IN) process for at most this number of hits
				pPars->pMultiHits,				// where to return the loci for each hit by current read
				pPars->NumIdentNodes,			// memory has been allocated by caller for holding up to this many tsIdentNodes
				pPars->pIdentNodes);			// memory allocated by caller for holding tsIdentNodes

		if (HitRslt == eHRhits && m_pPriorityRegionBED != nullptr)
			{
			int NumInPriorityRegions;
			int NumInNonPriorityRegions;
			tsHitLoci* pPriorityHits;

			NumInPriorityRegions = 0;
			NumInNonPriorityRegions = 0;
			pHit = pPars->pMultiHits;
			pPriorityHits = pHit;
			for (HitIdx = 0; HitIdx < min(LowHitInstances, MaxML + RefExacts); HitIdx++, pHit++)
				{
				// check if hit loci within region designated as being a priority exact matching region
				if (m_pSfxArray->GetIdentName(pHit->Seg->ChromID, sizeof(szPriorityChromName), szPriorityChromName) != eBSFSuccess)
					{
					NumInNonPriorityRegions += 1;
					continue;
					}
				if ((PriorityChromID = m_pPriorityRegionBED->LocateChromIDbyName(szPriorityChromName)) < 1)
					{
					NumInNonPriorityRegions += 1;
					continue;
					}
				if (!m_pPriorityRegionBED->InAnyFeature(PriorityChromID, (int)pHit->Seg->MatchLoci, (int)(pHit->Seg->MatchLoci + pHit->Seg->MatchLen - 1)))
					{
					NumInNonPriorityRegions += 1;
					continue;
					}
				if (NumInNonPriorityRegions > 0)
					*pPriorityHits++ = *pHit;
				NumInPriorityRegions += 1;
				}
			if (NumInPriorityRegions > 0 && (m_bClampMaxMLmatches || NumInPriorityRegions <= MaxML))
				{
				if (NumInPriorityRegions > MaxML)
					NumInPriorityRegions = MaxML;
				LowHitInstances = NumInPriorityRegions;
				bProcNorm = false;
				}
			}
		if (bProcNorm)
			{
			LowHitInstances = pReadHit->LowHitInstances;
			LowMMCnt = pReadHit->LowMMCnt;
			NxtLowMMCnt = pReadHit->NxtLowMMCnt;
			}
		}
	else
		bProcNorm = true;

	if (bProcNorm)
		{
		if (m_bLocateBestMatches)
			{
			HitRslt =						// < 0 if errors, 0 if no matches, 1..MaxHits, or MaxHits+1 if additional matches have been sloughed
					m_pSfxArray->LocateBestMatches(pReadHit->ReadID,			// identifies this read
						MaxTotMM,			        // return matches having at most this number of mismatches
						CoreLen,					// core window length
						CoreDelta,					// core window offset increment (1..n)
						MaxNumSlides,				// max number of times to slide core on each strand
						pPars->AlignStrand,			// watson, crick or both?
						pPars->Sequence, pPars->MatchLen,
						MaxML,						// process for at most this many hits by current read
						&LowHitInstances,			// returned number of match instances in pHits
						pPars->pMultiHits,			// where to return the loci for each hit by current read
						pPars->MaxIter,					// max allowed iterations per subsegmented sequence when matching that subsegment
						pPars->NumIdentNodes,		// memory has been allocated by caller for holding up to this many tsIdentNodes
						pPars->pIdentNodes);		// memory allocated by caller for holding tsIdentNodes
				if (HitRslt == 0)
					HitRslt = eHRnone;
				else
					if (HitRslt >= 1)
						HitRslt = eHRhits;
			}
		else
			HitRslt = m_pSfxArray->AlignReads(0,	// flags indicating if lower levels need to do any form of extended processing with this specific read...
					pReadHit->ReadID,				// identifies this read
					pPars->MinChimericLen,			// minimum chimeric length as a percentage (0 to disable, otherwise 50..99) of probe sequence
					MaxTotMM, CoreLen, CoreDelta,
					MaxNumSlides,					// limit on number of times core window can be moved or slide to right over read per 100bp of read length
					pPars->MinCoreLen,				// minimum core length allowed
					pPars->MinEditDist,
					pPars->AlignStrand,				// watson, crick or both?
					pPars->microInDelLen,			// microInDel length maximum
					pPars->SpliceJunctLen,			// maximum splice junction length when aligning RNAseq reads
					&LowHitInstances,
					&LowMMCnt,
					&NxtLowMMCnt,
					pPars->Sequence, pPars->MatchLen,
					MaxML,					// process for at most this many hits by current read
					pPars->pMultiHits,				// where to return the loci for each hit by current read
					pPars->NumIdentNodes,
					pPars->pIdentNodes);
		}

		// user may be interested in multihits upto m_MaxMLmatches limit even if there were actually many more so in this case remap the hit result
	if (LowHitInstances > MaxML)
		LowHitInstances = MaxML + 1;
	if (m_bClampMaxMLmatches && HitRslt == eHRHitInsts)
		{
		LowHitInstances = MaxML;
		HitRslt = eHRhits;
		pPars->bForceNewAlignment = true;
		}
	pPars->PrevHitRslt = HitRslt;
	pPars->PrevLowHitInstances = LowHitInstances;
	pPars->PrevLowMMCnt = LowMMCnt;
	pPars->PrevNxtLowMMCnt = NxtLowMMCnt;
	pPars->PrevMatchLen = pPars->MatchLen;
	pPars->pPrevReadHit = pReadHit;
	memmove(pPars->PrevSequence, pPars->Sequence, pPars->MatchLen);
	}
else		// else, identical sequence, reuse previous hit results from AlignReads
	{
	HitRslt = pPars->PrevHitRslt;
	LowHitInstances = pPars->PrevLowHitInstances;
	LowMMCnt = pPars->PrevLowMMCnt;
	NxtLowMMCnt = pPars->PrevNxtLowMMCnt;
	}

	// if SOLiD colorspace then need to normalise the hits back to as if aligned in basespace
if (m_bIsSOLiD)
	{
	switch (HitRslt) {
		case eHRHitInsts:
			if (m_MLMode != eMLall)
				break;
		case eHRhits:
			pHit = pPars->pMultiHits;
			for (HitIdx = 0; HitIdx < min(LowHitInstances, MaxML); HitIdx++, pHit++)
				{
				pHit->Seg[0].MatchLoci -= 1;
				pHit->Seg[0].MatchLen += 1;
				if (pHit->FlgInDel == 1 || pHit->FlgSplice == 1)
					pHit->Seg[1].ReadOfs += 1;
				}
			break;
		default:
			break;
		}
	}

	// ensure that TrimLeft/Right/Mismatches are consistent with the fact that trimming is a post alignment phase
	// if alignment was the result of a chimeric alignment then accept the flank trimming
if (HitRslt == eHRHitInsts || HitRslt == eHRhits)
	{
	pHit = pPars->pMultiHits;
	for (HitIdx = 0; HitIdx < min(LowHitInstances, MaxML); HitIdx++, pHit++)
		{
		if (pHit->FlgChimeric != 1)
			{
			pHit->Seg[0].TrimLeft = 0;
			pHit->Seg[0].TrimRight = 0;
			}
		pHit->Seg[0].TrimMismatches = pHit->Seg[0].Mismatches;
		if (pHit->FlgInDel == 1 || pHit->FlgSplice == 1)
			{
			pHit->FlgChimeric = 0;
			pHit->Seg[1].TrimLeft = 0;
			pHit->Seg[1].TrimRight = 0;
			pHit->Seg[1].TrimMismatches = pHit->Seg[1].Mismatches;
			}
		}
	}

Rslt = 0;
switch (HitRslt) {
	case eHRnone:			// no change or no hits
		pReadHit->NAR = eNARNoHit;
		pReadHit->LowMMCnt = 0;
		pReadHit->NumHits = 0;
		pReadHit->LowHitInstances = 0;
		pReadHit->HitLoci.Hit.FlgChimeric = 0;
		if (m_FMode == eFMsamAll)
			{
			pReadHit->LowMMCnt = (int8_t)0;
			if ((Rslt = WriteHitLoci(pPars, pReadHit, 0, pPars->pMultiHits)) < 0)
				break;
			}
		pPars->NumNonAligned += 1;
		pPars->bForceNewAlignment = false;
		break;

	case eHRhits:			// MMDelta criteria met and within the max allowed number of hits
		if (LowHitInstances > 1)
			pPars->bForceNewAlignment = true;
		else
			pPars->bForceNewAlignment = false;

		pReadHit->NAR = eNARAccepted;
		if (m_MLMode == eMLall)
			{
			int NumWriteHitLoci;
			// report each hit here...
			pReadHit->LowMMCnt = (int8_t)LowMMCnt;
			if ((Rslt = NumWriteHitLoci = WriteHitLoci(pPars, pReadHit, LowHitInstances, pPars->pMultiHits)) < 0)
				break;
			if (NumWriteHitLoci)
				{
				pPars->NumAcceptedAsAligned += 1;
				pPars->NumLociAligned += NumWriteHitLoci;
				if (NumWriteHitLoci == 1)
					pPars->TotAcceptedAsUniqueAligned += 1;
				else
					pPars->TotAcceptedAsMultiAligned += 1;
				}
			HitRslt = eHRHitInsts;
			break;
			}

		pPars->NumAcceptedAsAligned += 1;
		pPars->NumLociAligned += LowHitInstances;

		if (LowHitInstances == 1)
			pPars->TotAcceptedAsUniqueAligned += 1;
		else
			pPars->TotAcceptedAsMultiAligned += 1;

		if(LowHitInstances > 0)
			m_MultiHitDist[LowHitInstances - 1] += 1;

		if(m_PEproc == ePEdefault)   // not PE processing
			{
			if (m_MLMode == eMLrand)
				RandIdx = MaxML == 1 ? 0 : (rand() % LowHitInstances);
			else
				RandIdx = 0;
			if (m_MLMode <= eMLrand || LowHitInstances == 1)
				{
				if ((m_MLMode == eMLdist && LowHitInstances == 1) || m_MLMode != eMLdist)
					{
					if (LowHitInstances == 1)	// was this a unique hit?
						pReadHit->HitLoci.FlagHL = (int)eHLunique;
					else
						pReadHit->HitLoci.FlagHL = (int)eHLrandom;
					LowHitInstances = 1;			// currently just accepting one
					pReadHit->NumHits = 1;
					pReadHit->HitLoci.Hit = pPars->pMultiHits[RandIdx];
					pReadHit->HitLoci.FlagSegs = (pReadHit->HitLoci.Hit.FlgInDel == 1 || pReadHit->HitLoci.Hit.FlgSplice == 1) ? 1 : 0;
					}
				else
					{
					pReadHit->NAR = eNARMultiAlign;
					pReadHit->NumHits = 0;
					}
				}
			else
				if (LowHitInstances == 1)
					{
					pReadHit->NAR = eNARAccepted;
					pReadHit->NumHits = 1;
					}
				else
					{
					pReadHit->NAR = eNARMultiAlign;
					pReadHit->NumHits = 0;
					}
			}
		else   // else is PE processing
			{
			if(LowHitInstances == 1)
				{
				pReadHit->NAR = eNARAccepted;
				pReadHit->NumHits = 1;
				pReadHit->HitLoci.Hit = pPars->pMultiHits[0];
				pReadHit->HitLoci.FlagSegs = (pReadHit->HitLoci.Hit.FlgInDel == 1 || pReadHit->HitLoci.Hit.FlgSplice == 1) ? 1 : 0;
				}
			else
				{
				pReadHit->NAR = eNARMultiAlign;
				pReadHit->NumHits = LowHitInstances;
				}
			}

		pReadHit->LowHitInstances = (int16_t)LowHitInstances;
		pReadHit->LowMMCnt = (int8_t)LowMMCnt;
		pReadHit->NxtLowMMCnt = (int8_t)NxtLowMMCnt;
		// handling multiple loci aligned reads, make a copy of each multihit loci
		if ((m_PEproc != ePEdefault && LowHitInstances > 1) || m_MLMode > eMLrand)
			{
			int HitIdx;
			tsHitLoci* pHit;
			tsReadHit* pMHit;
			pHit = pPars->pMultiHits;
			pMHit = &pPars->HitReads[pPars->HitReadsOfs];
			for (HitIdx = 0; HitIdx < LowHitInstances; HitIdx++, pMHit++, pHit++)
				{
				memcpy(pMHit, pReadHit, sizeof(tsReadHit));
				pMHit->FlgPEAligned = false;
				pMHit->HitLoci.Hit = *pHit;
				pMHit->HitLoci.FlagSegs = (pMHit->HitLoci.Hit.FlgInDel == 1 || pMHit->HitLoci.Hit.FlgSplice == 1) ? 1 : 0;
				pMHit->HitLoci.FlagMH = LowHitInstances > 1 ? 1 : 0;
				pMHit->HitLoci.FlagMHA = 0;
				}
			
			if(m_PEproc == ePEdefault)
				if ((Rslt = AddMHitReads(LowHitInstances, &pPars->HitReads[pPars->HitReadsOfs])) < 0)		// pts to array of hit loci
					break;
			pPars->HitReadsOfs = HitIdx;
			}
		// finish handling multiple aligned reads
		break;

	case eHRMMDelta:			// same or a new LowMMCnt unique hit but MMDelta criteria not met
		pPars->NumNotAcceptedDelta += 1;
		memset(&pReadHit->HitLoci.Hit.Seg[0], 0, sizeof(tsSegLoci));
		pReadHit->NAR = eNARMMDelta;
		pReadHit->NumHits = 0;
		pReadHit->HitLoci.Hit.Seg[0].Strand = '?';
		pReadHit->HitLoci.Hit.BisBase = eBaseN;
		pReadHit->HitLoci.Hit.Seg[0].MatchLen = (uint16_t)ProbeLen;
		pReadHit->LowHitInstances = (int16_t)LowHitInstances;
		pReadHit->LowMMCnt = (int8_t)LowMMCnt;
		pReadHit->NxtLowMMCnt = (int8_t)NxtLowMMCnt;
		pPars->bForceNewAlignment = false;
		break;

	case eHRHitInsts:			// same or new LowMMCnt but now simply too many multiple hit instances, treat as none-aligned
		pReadHit->NAR = eNARMultiAlign;
		if (m_MLMode == eMLall && m_FMode == eFMsamAll)
			{
			pReadHit->LowMMCnt = (int8_t)0;
			if ((Rslt = WriteHitLoci(pPars, pReadHit, 0, pPars->pMultiHits)) < 0)
				break;
			}
		pPars->NumNonAligned += 1;
		pPars->bForceNewAlignment = false;
		break;

	case eHRRMMDelta:			// reduced NxtLowMMCnt only
		pReadHit->NxtLowMMCnt = (int8_t)NxtLowMMCnt;
		pPars->bForceNewAlignment = false;
		break;
	}
	
if (Rslt < 0)
	{
	AcquireSerialise();
	m_ThreadCoredApproxRslt = Rslt;
	ReleaseSerialise();
	return(eHRFatalError);
	}


if (HitRslt == eHRHitInsts)
	{
	pReadHit->NumHits = 0;
	pReadHit->NAR = eNARMultiAlign;
	memset(&pReadHit->HitLoci.Hit.Seg[0], 0, sizeof(tsSegLoci));
	pReadHit->HitLoci.Hit.Seg[0].Strand = '?';
	pReadHit->HitLoci.Hit.BisBase = eBaseN;
	pReadHit->HitLoci.Hit.Seg[0].MatchLen = (uint16_t)ProbeLen;
	pReadHit->LowHitInstances = (int16_t)LowHitInstances;
	pReadHit->LowMMCnt = (int8_t)LowMMCnt;
	pReadHit->NxtLowMMCnt = (int8_t)NxtLowMMCnt;
	}

	// if SOLiD colorspace then need to restore any hits to what they were immediately following return from AlignReads
	// as when SOLiD processing MatchLoci/MatchLen and ReadOfs would have all been normalised as if basespace processing and if
	// a duplicate of this read sequence follows then hits will be simply reused w/o the overhead of a call to AlignReads
if (m_bIsSOLiD)
	{
	switch (HitRslt) {
		case eHRHitInsts:
			if (m_MLMode != eMLall)
				break;
		case eHRhits:
			pHit = pPars->pMultiHits;
			for (HitIdx = 0; HitIdx < min(LowHitInstances, MaxML); HitIdx++, pHit++)
				{
				pHit->Seg[0].MatchLoci += 1;
				pHit->Seg[0].MatchLen -= 1;
				if (pHit->FlgInDel == 1 || pHit->FlgSplice == 1)
					pHit->Seg[1].ReadOfs -= 1;
				}
			break;
		default:
			break;
		}
	}
return(HitRslt);
}



int
CKAligner::ProcCoredApprox(tsThreadMatchPars *pPars)
{
int Rslt;
int PE1HitRslt;
int PE2HitRslt;
int ReadsHitIdx;
int ChkPE1Ofs;
int ChkPE2Ofs;
bool bMultiAligned;
bool bMultiAlignedAccepted;

int MaxML = max(m_MaxMLmatches,m_MaxMLPEmatches);

tsReadHit* pPE1ReadHit;					// current PE1 or SE read being processed
tsReadHit* pPE2ReadHit;					// current PE2 read being processed
tsReadsHitBlock *pReadsHitBlock;		// block of reads for this thread to process

int PE1ReadMHOfs;
int PE2ReadMHOfs;

tsReadHit *pPE1ReadMH;					// current multihit PE1 read being processed
tsReadHit* pPE2ReadMH;					// current multihit PE2 read being processed
tsReadHitLoci PE1ProvHitLoci;				// provisional PE1 hit loci if multiloci hits
tsReadHitLoci PE2ProvHitLoci;				// provisional PE2 hit loci if multiloci hits


if((pReadsHitBlock = (tsReadsHitBlock *)calloc(1,sizeof(tsReadsHitBlock)))==nullptr)
	return(-1);
pReadsHitBlock->MaxReads = cMaxReadsPerBlock;

// iterate each read sequence starting from the first
// assumes that reads will have been sorted by ReadID
pPars->PlusHits = 0;
pPars->MinusHits = 0;
pPars->ChimericHits = 0;
pPars->NumReadsProc = 0;
pPE1ReadHit = nullptr;
pPE2ReadHit = nullptr;

pReadsHitBlock->NumReads = 0;
Rslt = 0;
pPars->MaxIter = m_pSfxArray->GetMaxIter();
pPars->bForceNewAlignment = true;

// m_hRwLock will be released and regained within ThreadedIterReads so need to always have acquired a read lock before calling ThreadedIterReads
AcquireSharedLock();
while(ThreadedIterReads(pReadsHitBlock))
	{
	Rslt = 0;
	for(ReadsHitIdx = 0; ReadsHitIdx < pReadsHitBlock->NumReads; ReadsHitIdx++)
		{
		pPE1ReadHit = pReadsHitBlock->pReadHits[ReadsHitIdx];
		pPars->ReadsHitIdx = ReadsHitIdx;
		pPars->bIsPE2= false;
		PE1ReadMHOfs = pPars->HitReadsOfs = 0;
		PE1HitRslt = AlignRead(pPE1ReadHit, pPars);
		if(PE1HitRslt == eHRFatalError)
			{
			ReleaseSharedLock();
			free(pReadsHitBlock);
			return(-1);
			}
		if(m_PEproc != ePEdefault)
			{
			pPE2ReadHit = pReadsHitBlock->pReadHits[++ReadsHitIdx];
			pPars->ReadsHitIdx = ReadsHitIdx;
			pPars->bIsPE2 = true;
			PE2ReadMHOfs = pPars->HitReadsOfs;
			PE2HitRslt = AlignRead(pPE2ReadHit, pPars);
			if (PE2HitRslt == eHRFatalError)
				{
				ReleaseSharedLock();
				free(pReadsHitBlock);
				return(-1);
				}
			if (PE1HitRslt != eHRhits || PE2HitRslt != eHRhits)
				continue;

			if(pPE1ReadHit->LowHitInstances == 1 && pPE2ReadHit->LowHitInstances == 1)
				continue;

				
			if(pPE1ReadHit->LowHitInstances >= m_MaxMLPEmatches || pPE2ReadHit->LowHitInstances >= m_MaxMLPEmatches)
				continue;

			bMultiAligned = false;
			bMultiAlignedAccepted = false;
			for(ChkPE1Ofs = 0; !(bMultiAligned && !bMultiAlignedAccepted) &&  ChkPE1Ofs < pPE1ReadHit->LowHitInstances; ChkPE1Ofs++)
				{
				if (pPE1ReadHit->LowHitInstances == 1)
					pPE1ReadMH = pPE1ReadHit;
				else
					pPE1ReadMH = &pPars->HitReads[PE1ReadMHOfs + ChkPE1Ofs];
				pPE1ReadMH->NumHits = 1;
				for (ChkPE2Ofs = 0; ChkPE2Ofs < pPE2ReadHit->LowHitInstances; ChkPE2Ofs++)
					{
					if (pPE2ReadHit->LowHitInstances == 1)
						pPE2ReadMH = pPE2ReadHit;
					else
						pPE2ReadMH = &pPars->HitReads[PE2ReadMHOfs+ ChkPE2Ofs];
					pPE2ReadMH->NumHits = 1;

					int SeqFragLen = AcceptProvPE(m_PairMinLen, m_PairMaxLen, m_bPairStrand, pPE1ReadMH, pPE2ReadMH);
					if(SeqFragLen > 0)
						{
						if(!bMultiAligned)
							{
							PE1ProvHitLoci = pPE1ReadMH->HitLoci;
							PE2ProvHitLoci = pPE2ReadMH->HitLoci;
							bMultiAligned = true;
							bMultiAlignedAccepted = true;
							}
						else
							{
							bMultiAlignedAccepted = false;
							break;
							}
						}
					}
				}

			if(bMultiAlignedAccepted)
				{
				pPE1ReadHit->HitLoci = PE1ProvHitLoci;
				pPE1ReadHit->NAR = eNARAccepted;
				pPE1ReadHit->NumHits = 1;
				pPE2ReadHit->HitLoci = PE2ProvHitLoci;
				pPE2ReadHit->NAR = eNARAccepted;
				pPE2ReadHit->NumHits = 1;
				}
			}
		}
	}
ReleaseSharedLock();
AcquireSerialise();
m_NumSloughedNs += pPars->NumSloughedNs;
m_TotNonAligned += pPars->NumNonAligned;
m_TotAcceptedAsUniqueAligned += pPars->TotAcceptedAsUniqueAligned;
m_TotAcceptedAsMultiAligned += pPars->TotAcceptedAsMultiAligned;
m_TotAcceptedAsAligned += pPars->NumAcceptedAsAligned;
m_TotLociAligned += pPars->NumLociAligned;
m_TotNotAcceptedDelta += pPars->NumNotAcceptedDelta;
m_TotAcceptedHitInsts += pPars->NumAcceptedHitInsts;

if(m_MLMode != eMLall)
	{
	int Idx;
	for(Idx = 0; Idx < MaxML; Idx++)
		m_MultiHitDist[Idx] += pPars->MultiHitDist[Idx];
	}
free(pReadsHitBlock);
ReleaseSerialise();
return(1);
}

// LocateRead
tsReadHit *
CKAligner::LocateRead(uint32_t ReadID) // Locate read with requested ReadID
{
int Rslt;
tsReadHit *pProbe;
int Lo,Mid,Hi;	// search limits

if(m_ppReadHitsIdx == nullptr || m_AllocdReadHitsIdx < m_NumReadsLoaded)
	return(nullptr);

Lo = 0;
Hi = m_NumReadsLoaded-1;
while(Hi >= Lo) {
	Mid = (Hi + Lo)/2;
	pProbe = m_ppReadHitsIdx[Mid];
	Rslt = pProbe->ReadID - ReadID;
	if(Rslt > 0)
		{
		Hi = Mid - 1;
		continue;
		}
	if(Rslt < 0)
		{
		Lo = Mid + 1;
		continue;
		}
	return(pProbe);
	}
return(nullptr);
}


int
CKAligner::AddMHitReads(uint32_t NumHits,	// number of multimatches loci in pHits
		tsReadHit *pHits)		// pts to array of hit loci
{
size_t memreq;
tsReadHit *pDstHits;
// ensure actually processing multihits
if(m_MLMode <= eMLrand)
	return(0);					// silently slough these hits

AcquireSerialiseMH();

if((m_AllocdMultiHits - m_NumMultiHits) < (NumHits+1000))	// need to realloc? -- added 1000 to provide a little safety margin
	{
	memreq = (m_AllocdMultiHits + cAllocMultihits) * sizeof(tsReadHit);
#ifdef _WIN32
	pDstHits = (tsReadHit *) realloc(m_pMultiHits,memreq);
#else
	pDstHits = (tsReadHit *)mremap(m_pMultiHits,m_AllocdMultiHitsMem,memreq,MREMAP_MAYMOVE);
	if(pDstHits == MAP_FAILED)
		pDstHits = nullptr;
#endif
	if(pDstHits == nullptr)
		{
		ReleaseSerialiseMH();
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddMHitReads: Memory re-allocation to %zd bytes - %s",memreq,strerror(errno));
		return(eBSFerrMem);
		}
	m_pMultiHits = pDstHits;
	m_AllocdMultiHitsMem = memreq;
	m_AllocdMultiHits += cAllocMultihits;
	}
pDstHits = &m_pMultiHits[m_NumMultiHits];
memmove(pDstHits,pHits,sizeof(tsReadHit) * NumHits);
m_NumMultiHits += NumHits;
if(NumHits == 1)
	m_NumUniqueMultiHits += 1;
else
	m_NumProvMultiAligned += 1;
ReleaseSerialiseMH();
return((int)NumHits);
}


// ResetThreadedIterReads
void
CKAligner::ResetThreadedIterReads(void) // must be called by master thread prior to worker threads calling ThreadedIterReads()
{
m_NumReadsProc = 0;
m_NxtReadProcOfs = 0;
m_ProcessingStartSecs = gStopWatch.ReadUSecs();
}


uint32_t		// Returns the number of reads thus far loaded and processed for alignment
CKAligner::ApproxNumReadsProcessed(uint32_t *pNumProcessed,uint32_t *pNumLoaded)
{
uint32_t NumReadsProc;
AcquireSerialise();
NumReadsProc = m_NumReadsProc;
if(pNumProcessed != nullptr)
	*pNumProcessed = NumReadsProc;
if(pNumLoaded != nullptr)
	*pNumLoaded = m_NumReadsLoaded;
ReleaseSerialise();
return(NumReadsProc);
}

// ThreadedIterReads
// Iterates over all reads
// The body of this function is serialised
bool	// returns false if no more reads avail for processing by calling thread
CKAligner::ThreadedIterReads(tsReadsHitBlock *pRetBlock)
{
uint32_t NumReadsLeft;
uint32_t MaxReads2Proc;
uint32_t AdjReadsPerBlock;
tsReadHit *pCurReadHit;
pRetBlock->NumReads = 0;

AdjReadsPerBlock = cMaxReadsPerBlock;
if(m_SampleNthRawRead > 1)
	AdjReadsPerBlock = min((uint32_t)100,AdjReadsPerBlock/m_SampleNthRawRead);

ReleaseSharedLock();
while(1) {
	AcquireSerialise();
	AcquireSharedLock();
	if(m_bAllReadsLoaded || ((m_NumReadsLoaded - m_NumReadsProc) >= (uint32_t)min(AdjReadsPerBlock,(uint32_t)pRetBlock->MaxReads)) || m_ThreadCoredApproxRslt < 0)
		break;

	ReleaseSharedLock();
	ReleaseSerialise();
#ifdef _WIN32
	Sleep(2000);			// must have caught up to the reads loader, allow it some breathing space to parse and load some more reads...
#else
	sleep(2);
#endif
	}

if(m_pReadHits == nullptr ||
	m_ThreadCoredApproxRslt < 0 ||
	(m_bAllReadsLoaded && (m_LoadReadsRslt != eBSFSuccess)) ||
	m_bAllReadsLoaded && (m_NumReadsLoaded == 0 || m_NumReadsProc == m_NumReadsLoaded)) // if all reads have been loaded and all processed then time to move onto next processing phase
	{
	pRetBlock->NumReads = 0;
	pRetBlock->pReadHits[0] = nullptr;
	ReleaseSerialise();
	return(false);
	}

// adjust pRetBlock->MaxReads according to the number of reads remaining and threads still processing these reads
// idea is to maximise the number of threads still processing when most reads have been processed so that
// the last thread processing doesn't end up with a large block of reads needing lengthly processing
NumReadsLeft = m_NumReadsLoaded - m_NumReadsProc;
if(NumReadsLeft < AdjReadsPerBlock/4)	// if < cMaxReadsPerBlock/4 yet to be processed then give it all to the one thread
	MaxReads2Proc = NumReadsLeft;
else
	{
	MaxReads2Proc = min((uint32_t)pRetBlock->MaxReads,10 + (NumReadsLeft / (uint32_t)m_NumThreads));
	// assume PE processing so ensure MaxReads2Proc is a multiple of 2
	MaxReads2Proc &= ~0x01;
	}
MaxReads2Proc = min(MaxReads2Proc,NumReadsLeft);
if(!m_NumReadsProc)
	m_NxtReadProcOfs = 0;
pCurReadHit = (tsReadHit *)((uint8_t *)m_pReadHits + m_NxtReadProcOfs);

while(MaxReads2Proc)
	{
	pRetBlock->pReadHits[pRetBlock->NumReads++] = pCurReadHit;
	pCurReadHit = (tsReadHit *)((uint8_t *)pCurReadHit + sizeof(tsReadHit) + pCurReadHit->ReadLen + pCurReadHit->DescrLen);
	MaxReads2Proc -= 1;
	}

m_NumReadsProc += pRetBlock->NumReads;
m_NxtReadProcOfs = (size_t)((uint8_t *)pCurReadHit - (uint8_t *)m_pReadHits);

ReleaseSerialise();
return(true);
}


// IterReads
// use to iterate over reads returning ptr to next read following the current read
// Reads need not be sorted
// To start from first read then pass in nullptr as pCurReadHit
// Returns nullptr if all read hits have been iterated
tsReadHit *
CKAligner::IterReads(tsReadHit *pCurReadHit) // to start from first read then pass in nullptr as pCurReadHit
{
tsReadHit *pNxtReadHit = nullptr;
if(pCurReadHit == nullptr)
	pNxtReadHit = m_pReadHits;
else
	if(pCurReadHit->ReadID != m_FinalReadID)
		pNxtReadHit = (tsReadHit *)((uint8_t *)pCurReadHit + sizeof(tsReadHit) + pCurReadHit->ReadLen + pCurReadHit->DescrLen);
return(pNxtReadHit);
}

// IterSortedReads
// use to iterate over sorted reads returning ptr to next read following the current read
// To start from first read then pass in nullptr as pCurReadHit
// Returns nullptr if all read hits have been iterated
tsReadHit *
CKAligner::IterSortedReads(tsReadHit *pCurReadHit)
{
tsReadHit *pNxtReadHit = nullptr;
if(pCurReadHit == nullptr)
	pNxtReadHit = m_ppReadHitsIdx[0];
else
	if(pCurReadHit->ReadHitIdx < m_NumReadsLoaded)
		pNxtReadHit = m_ppReadHitsIdx[pCurReadHit->ReadHitIdx];
return(pNxtReadHit);
}

tsReadHit *		// returned read which overlaps StartLoci and EndLoci, nullptr if no read located
CKAligner::IterateReadsOverlapping(bool bTriSNPs, // false if iterating DiSNPs, true if iterating TriSNPs
						tsChromSNPs *pChromSNPs, // processing SNPs on this chromosome
						int StartLoci,				// returned reads are required to overlap both this starting and
						int EndLoci)				// this ending loci
{
bool bNewStartEnd;
int CurStart;
int CurEnd;
tsAdjacentSNPs *pAdjacentSNPs;
tsReadHit *pNxtReadHit = nullptr;

pAdjacentSNPs = bTriSNPs ? &pChromSNPs->AdjacentSNPs[1] : &pChromSNPs->AdjacentSNPs[0];
bNewStartEnd = false;
if(pAdjacentSNPs->pFirstIterReadHit == nullptr ||				// only nullptr if iterations just starting for SNPs on this chrom
	StartLoci <  pAdjacentSNPs->StartLoci || EndLoci <  pAdjacentSNPs->EndLoci)	// or if StartLoci/EndLoci now 5' to previous so needing to search from start
	{
	pNxtReadHit = pChromSNPs->pFirstReadHit;
	pAdjacentSNPs->pFirstIterReadHit = pChromSNPs->pFirstReadHit;
	pAdjacentSNPs->pPrevIterReadHit = nullptr;
	pAdjacentSNPs->StartLoci = StartLoci;
	pAdjacentSNPs->EndLoci = EndLoci;
	bNewStartEnd = true;
	}
else
	{
	if(StartLoci == pAdjacentSNPs->StartLoci && EndLoci == pAdjacentSNPs->EndLoci) // same start/end loci as previous so must be iterating same start/end loci
		{
		if(pAdjacentSNPs->pPrevIterReadHit == pChromSNPs->pLastReadHit)		// check if last returned was last on chrom
			return(nullptr);

		pNxtReadHit = m_ppReadHitsIdx[pAdjacentSNPs->pPrevIterReadHit->ReadHitIdx]; // recommence search from next
		}
	else   // else starting new Start/EndLoci which is 3' to previous
		{
		pNxtReadHit = pAdjacentSNPs->pFirstIterReadHit;
		pAdjacentSNPs->pPrevIterReadHit = pNxtReadHit;
		pAdjacentSNPs->StartLoci = StartLoci;
		pAdjacentSNPs->EndLoci = EndLoci;
		bNewStartEnd = true;
		}	
	}

	// have the current read, iterate forward returning reads until reads align past the StartLoci
do {
	if(pNxtReadHit == nullptr)
		return(nullptr);
	if(pNxtReadHit->NAR == eNARAccepted &&
		!(pNxtReadHit->HitLoci.Hit.FlgInDel || pNxtReadHit->HitLoci.Hit.FlgSplice))
		{
		if(pNxtReadHit->HitLoci.Hit.Seg[0].ChromID !=  pChromSNPs->ChromID) // double check read is aligning on to expected chrom
			return(nullptr);

		CurStart = AdjStartLoci(&pNxtReadHit->HitLoci.Hit.Seg[0]); 
		CurEnd = AdjEndLoci(&pNxtReadHit->HitLoci.Hit.Seg[0]);

		if(CurStart <= StartLoci && CurEnd >= EndLoci)   // this read overlaps both start and end?
			{
			pAdjacentSNPs->pPrevIterReadHit = pNxtReadHit;
			if(bNewStartEnd)		
				pAdjacentSNPs->pFirstIterReadHit = pNxtReadHit;
			return(pNxtReadHit);
			}
		else
			if(CurStart > StartLoci)
				return(nullptr);
		}

	}
while((pNxtReadHit != pChromSNPs->pLastReadHit) && (pNxtReadHit = m_ppReadHitsIdx[pNxtReadHit->ReadHitIdx])!=nullptr);

return(nullptr);
}


// returns counts for all 3 frame shifted codon combination counts for read alignments at SNP loci
// frame shifts start at -2,-1, 0 relative to SNP loci
// codon combinations are relative to sense strand
int
CKAligner::FrameShiftCodons(uint32_t ChromID,					// SNP loci is on this chromosome
				 int Loci,				// returned codons are frame shifted relative to this loci (-2..0,-1..1,0..2)
				 int *pCodons)	// where to return cnts of frame shifted codons (set of codon counts [3][64])
{
int NumReads;
int InFrameCodons[3][64];
int FrameShift;
int FrameCodonIdx;
int RelLoci;
etSeqBase LociBase;
tsReadHit* pCurReadHit = nullptr;
if(pCodons == nullptr)
	return(0);

memset(pCodons,0,sizeof(InFrameCodons));
memset(InFrameCodons, 0, sizeof(InFrameCodons));

if((pCurReadHit = Locate1stReadOverlapLoci(ChromID,Loci))==nullptr)
	return(0);
NumReads = 0;
do
	{
	if(pCurReadHit->NAR != eNARAccepted || pCurReadHit->HitLoci.Hit.FlgInDel || pCurReadHit->HitLoci.Hit.FlgSplice)
		continue;
	NumReads += 1;
	// for each frame shift then get codon combination
	for(FrameShift = 0; FrameShift < 3; FrameShift++)
		{
		FrameCodonIdx = 0;
		for(RelLoci = -2; RelLoci <= 0; RelLoci++)
			{
			if((LociBase = AdjAlignSNPBase(pCurReadHit, ChromID, Loci + FrameShift + RelLoci)) > eBaseT)
				break;
			FrameCodonIdx <<= 2;
			FrameCodonIdx |= LociBase;
			}
		if(LociBase <= eBaseT)
			InFrameCodons[FrameShift][FrameCodonIdx] += 1;
		}
	}
while((pCurReadHit = IterReadsOverlapLoci(ChromID, Loci, pCurReadHit))!=nullptr);
memcpy(pCodons,InFrameCodons,sizeof(InFrameCodons));
return(NumReads);
}

tsReadHit *				// returned lowest sorted read which overlaps ChromID.Loci, nullptr if non overlaps
CKAligner::Locate1stReadOverlapLoci(uint32_t ChromID,				// loci is on this chromosome
				 int Loci)							// returned read is lowest sorted read which overlaps this loci
{
tsReadHit *pCurRead;

int CurStart;
int CurEnd;
int64_t MidIdx;
int64_t HiIdx = m_NumReadsLoaded - 1;
int64_t LoIdx = 0;
int Start1st = 0;
tsReadHit* p1stOverlapRead = nullptr;
while(HiIdx >= LoIdx) {
	MidIdx = (HiIdx + LoIdx)/2;
	pCurRead = m_ppReadHitsIdx[MidIdx];
	if(pCurRead->HitLoci.Hit.Seg[0].ChromID == ChromID)
		{
		// read aligned to requested chrom
		CurStart = AdjStartLoci(&pCurRead->HitLoci.Hit.Seg[0]);
		CurEnd = AdjEndLoci(&pCurRead->HitLoci.Hit.Seg[0]);
		if(Loci >= CurStart && Loci <= CurEnd)	// have an overlap?
			{
			Start1st = CurStart;
			p1stOverlapRead = pCurRead;	// record the overlap, might be the only!
			while(MidIdx--)     // doing a backward linear search ...
				{
				pCurRead = m_ppReadHitsIdx[MidIdx];
				CurStart = AdjStartLoci(&pCurRead->HitLoci.Hit.Seg[0]);
				CurEnd = AdjEndLoci(&pCurRead->HitLoci.Hit.Seg[0]);
				if(Loci >= CurStart && Loci <= CurEnd)	// have a new overlap?
					{
					Start1st = CurStart;
					p1stOverlapRead = pCurRead;	// record the overlap, might be the first!
					}
				if((Start1st - CurStart) > (int)m_MaxReadsLen)
					break;	
				}
			return(p1stOverlapRead);
			}

		if(CurStart > Loci)		// CurStart after loci then reads sorted higher can't overlap - places an upper limit on search bounds
			HiIdx = MidIdx - 1;
		else                    // CurStart is <= loci but CurEnd is also < target loci
			{
			if((Loci - CurStart) > (int)m_MaxReadsLen)
				LoIdx = MidIdx + 1;
			else
				LoIdx += 1;
			}
		continue;
		}

	// need to locate chromosome before can start looking for read overlaps onto loci
	if(pCurRead->HitLoci.Hit.Seg[0].ChromID < ChromID)
		LoIdx = MidIdx+1;
	else
		HiIdx = MidIdx - 1;
	};

return(p1stOverlapRead);
}

tsReadHit*													// located read, or nullptr if unable to locate any more reads overlapping loci
CKAligner::IterReadsOverlapLoci(uint32_t ChromID,				// loci is on this chromosome
						   int Loci,						// returned reads are required to overlap this loci
						   tsReadHit* pPrevRead)			// any previously returned read which overlapped loci, nullptr if no read previously returned
{
bool bNewChrom;
int CurStart;
int CurEnd;
tsReadHit* pNxtReadHit = nullptr;

if(pPrevRead->NAR != eNARAccepted)
	pPrevRead = nullptr;

bNewChrom = (pPrevRead == nullptr || pPrevRead->HitLoci.Hit.Seg[0].ChromID != ChromID) ? true : false;
if(pPrevRead != nullptr)
	{
	CurStart = AdjStartLoci(&pPrevRead->HitLoci.Hit.Seg[0]);
	if(pPrevRead->HitLoci.Hit.Seg[0].ChromID < ChromID || Loci < CurStart)
		pPrevRead = nullptr;
	}
while((pNxtReadHit = IterSortedReads(pPrevRead))!=nullptr)
	{
	if(pNxtReadHit->NAR == eNARAccepted &&
	   !(pNxtReadHit->HitLoci.Hit.FlgInDel || pNxtReadHit->HitLoci.Hit.FlgSplice))
		{
		if(pNxtReadHit->HitLoci.Hit.Seg[0].ChromID != ChromID)	// double check read is aligning on to expected chrom
			if(bNewChrom)		// loci is on a different chromosome to that previously iterated so keep iterating 
				{
				if((pPrevRead = pNxtReadHit) == nullptr)
					return(nullptr);
				continue;
				}
			else
				return(nullptr);
		bNewChrom = false;
		CurStart = AdjStartLoci(&pNxtReadHit->HitLoci.Hit.Seg[0]);
		CurEnd = AdjEndLoci(&pNxtReadHit->HitLoci.Hit.Seg[0]);

		if(CurStart <= Loci && CurEnd >= Loci)					// this read overlaps loci?
			return(pNxtReadHit);
		if(CurStart > Loci)
			return(nullptr);
		}
	if((pPrevRead = pNxtReadHit) == nullptr)
		break;
	}
return(nullptr);
}


// NumDnUniques
// Returns number of unique loci to which reads map which are downstream of the specified current read
int											// returned number of unique reads
CKAligner::NumDnUniques(tsReadHit *pCurReadHit,		// current read
				int WinLen,					// only interested in unique reads starting within this window
				bool bStrandDep)			// if true then unique loci reads must be on current read stand
{
int NumUniques;
uint32_t CurChromID;
uint8_t CurStrand;
int CurStart;
int NxtStart;
int PrvNxtStart;
uint32_t NxtReadIdx;
tsReadHit *pNxtReadHit = nullptr;

if (pCurReadHit == nullptr)									// if nullptr then start from first
	{
	NxtReadIdx = 0;
	pCurReadHit = m_ppReadHitsIdx[NxtReadIdx];
	}
else
	if((NxtReadIdx = pCurReadHit->ReadHitIdx) == m_NumReadsLoaded)
		return(0);

// have the next read, iterate forward counting the uniques until loci outside of window
CurChromID = pCurReadHit->HitLoci.Hit.Seg[0].ChromID;
CurStart = AdjStartLoci(&pCurReadHit->HitLoci.Hit.Seg[0]);
PrvNxtStart = CurStart;
CurStrand = pCurReadHit->HitLoci.Hit.Seg[0].Strand;
NumUniques = 0;
do {
	pNxtReadHit = m_ppReadHitsIdx[NxtReadIdx];
	NxtReadIdx = pNxtReadHit->ReadHitIdx;
	if(CurChromID != pNxtReadHit->HitLoci.Hit.Seg[0].ChromID)	// must be on same chromosome
		return(NumUniques);
	if(pNxtReadHit->NAR != eNARAccepted)						// skip over any read not uniquely aligned
		continue;
	NxtStart = AdjStartLoci(&pNxtReadHit->HitLoci.Hit.Seg[0]);

	if((CurStart + WinLen) < NxtStart)					// outside of window?
		return(NumUniques);

	if(bStrandDep && CurStrand != pNxtReadHit->HitLoci.Hit.Seg[0].Strand) // has to be on same strand?
		continue;

	if(NxtStart != PrvNxtStart)
		{
		NumUniques += 1;
		PrvNxtStart = NxtStart;
		}
	}
while(NxtReadIdx != m_NumReadsLoaded);
return(NumUniques);
}


// NumUpUniques
// Returns number of unique loci to which reads map which are upstream of the specified current read
int											// returned number of unique reads
CKAligner::NumUpUniques(tsReadHit *pCurReadHit,		// current read
				int WinLen,					// only interested in unique reads starting within this window
				bool bStrandDep)			// if true then unique loci reads must be on current read stand
{
int NumUniques;
uint32_t CurChromID;
uint8_t CurStrand;
int CurStart;
int NxtStart;
int PrvNxtStart;
uint32_t NxtReadIdx;
tsReadHit *pNxtReadHit = nullptr;

if(pCurReadHit == nullptr || (NxtReadIdx = pCurReadHit->ReadHitIdx)==1) // if nullptr or first then
	return(0);								// can't be any upstream from the first!

CurChromID = pCurReadHit->HitLoci.Hit.Seg[0].ChromID;
CurStart = AdjStartLoci(&pCurReadHit->HitLoci.Hit.Seg[0]);
PrvNxtStart = CurStart;
CurStrand = pCurReadHit->HitLoci.Hit.Seg[0].Strand;
NumUniques = 0;
do {
	pNxtReadHit = m_ppReadHitsIdx[NxtReadIdx-1];
	if(CurChromID != pNxtReadHit->HitLoci.Hit.Seg[0].ChromID)		// must be on same chromosome
		return(NumUniques);
	if(pNxtReadHit->NAR != eNARAccepted)						// skip over any read not uniquely aligned
		continue;
	NxtStart = AdjStartLoci(&pNxtReadHit->HitLoci.Hit.Seg[0]);

	if(CurStart > WinLen && (CurStart - WinLen) > NxtStart)					// outside of window?
		return(NumUniques);

	if(bStrandDep && CurStrand != pNxtReadHit->HitLoci.Hit.Seg[0].Strand) // has to be on same strand?
		continue;

	if(NxtStart != PrvNxtStart)
		{
		NumUniques += 1;
		PrvNxtStart = NxtStart;
		}
	}
while(--NxtReadIdx > 0);
return(NumUniques);
}

int
CKAligner::SortReadHits(etReadsSortMode SortMode,		// sort mode required
				bool bSeqSorted,			// used to optimise eRSMSeq processing, if it is known that reads are already sorted in sequence order (loaded from pre-processed .rds file)
				bool bForce)				// if true then force sort
{
tsReadHit *pReadHit;
uint32_t Idx;

if(!bForce && SortMode == m_CurReadsSortMode && m_ppReadHitsIdx != nullptr && m_AllocdReadHitsIdx >= m_NumReadsLoaded)
	return(eBSFSuccess);					// if already in requested mode

if(m_ppReadHitsIdx == nullptr || m_AllocdReadHitsIdx < m_NumReadsLoaded)
	{
	if(m_ppReadHitsIdx != nullptr)
		{
		delete []m_ppReadHitsIdx;
		m_ppReadHitsIdx = nullptr;
		m_AllocdReadHitsIdx = 0;
		}
	if((m_ppReadHitsIdx = (tsReadHit **) new tsReadHit * [m_NumReadsLoaded])==nullptr)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"SortReadHits: Memory reads index allocation for %d ptrs - %s",m_NumReadsLoaded,strerror(errno));
		return(eBSFerrMem);
		}
	m_AllocdReadHitsIdx = m_NumReadsLoaded;
	}

if(SortMode != eRSMSeq || bSeqSorted)
	{
	pReadHit = nullptr;
	tsReadHit **pIdx = m_ppReadHitsIdx;
	while((pReadHit = IterReads(pReadHit))!=nullptr)
		*pIdx++ = pReadHit;
	}
else
	{
	pReadHit = m_pReadHits;
	size_t BuffOfs = 0;
	for(Idx = 0; Idx < m_NumReadsLoaded; Idx++)
		{
		m_ppReadHitsIdx[Idx] = pReadHit;
		BuffOfs += sizeof(tsReadHit) + pReadHit->ReadLen + pReadHit->DescrLen;
		pReadHit = (tsReadHit *)((char *)m_pReadHits+BuffOfs);
		}
	}

switch(SortMode) {
	case eRSMReadID:
		m_mtqsort.qsort(m_ppReadHitsIdx,m_NumReadsLoaded,sizeof(tsReadHit *),SortReadIDs);
		break;
	case eRSMPairReadID:
		m_mtqsort.qsort(m_ppReadHitsIdx,m_NumReadsLoaded,sizeof(tsReadHit *),SortPairReadIDs);
		break;
	case eRSMHitMatch:
		m_mtqsort.qsort(m_ppReadHitsIdx,m_NumReadsLoaded,sizeof(tsReadHit *),SortHitMatch);
		break;

	case eRSMPEHitMatch:
		m_mtqsort.qsort(m_ppReadHitsIdx,m_NumReadsLoaded,sizeof(tsReadHit *),SortPEHitMatch);
		break;		

	case eRSMSeq:
		if(!bSeqSorted)
			m_mtqsort.qsort(m_ppReadHitsIdx,m_NumReadsLoaded,sizeof(tsReadHit *),SortReadSeqs);
		break;
	default:
		break;
	}

// m_ppReadHitsIdx now in requested order, assign sequentially incrementing ReadHitIdx to the reads
for(Idx = 1; Idx <= m_NumReadsLoaded; Idx++)
	m_ppReadHitsIdx[Idx-1]->ReadHitIdx = Idx;

m_CurReadsSortMode = SortMode;
return(eBSFSuccess);
}


// SortReadIDs
// Sort reads by ascending read identifiers
int
CKAligner::SortReadIDs(const void *arg1, const void *arg2)
{
tsReadHit *pEl1 = *(tsReadHit **)arg1;
tsReadHit *pEl2 = *(tsReadHit **)arg2;

if(pEl1->ReadID < pEl2->ReadID )
		return(-1);
if(pEl1->ReadID > pEl2->ReadID )
	return(1);
return(0);
}

// SortPairReadIDs
// Sort paired reads reads by ascending PairReadID identifiers
int
CKAligner::SortPairReadIDs(const void *arg1, const void *arg2)
{
tsReadHit *pEl1 = *(tsReadHit **)arg1;
tsReadHit *pEl2 = *(tsReadHit **)arg2;
uint32_t El1ID = pEl1->PairReadID & 0x7fffffff;
uint32_t El2ID = pEl2->PairReadID & 0x7fffffff;

if(El1ID < El2ID)
	return(-1);
if(El1ID > El2ID)
	return(1);
if(pEl1->PairReadID & 0x80000000)		// if same PairReadID then the 3' read is after the 5' read
	return(1);
return(0);
}

// SortPEHitmatch
// Sort by NAR eNARAccepted with FlgPEAligned, then ascending chrom,  then PairReadID, then 5' forward read
int
CKAligner::SortPEHitMatch(const void *arg1, const void *arg2)
{
tsReadHit *pEl1 = *(tsReadHit **)arg1;
tsReadHit *pEl2 = *(tsReadHit **)arg2;
uint32_t El1ID;
uint32_t El2ID;

// order by eNARAccepted.FlgPEAligned
if((pEl1->NAR == eNARAccepted && pEl1->FlgPEAligned) && !(pEl2->NAR == eNARAccepted && pEl2->FlgPEAligned))
	return(-1);
if(!(pEl1->NAR == eNARAccepted && pEl1->FlgPEAligned) && (pEl2->NAR == eNARAccepted && pEl2->FlgPEAligned))
	return(1);
if(!(pEl1->NAR == eNARAccepted && pEl1->FlgPEAligned &&  pEl2->NAR == eNARAccepted && pEl2->FlgPEAligned))
	return(0);

// sort by Chrom ascending
if(pEl1->HitLoci.Hit.Seg[0].ChromID < pEl2->HitLoci.Hit.Seg[0].ChromID)
	return(-1);
if(pEl1->HitLoci.Hit.Seg[0].ChromID > pEl2->HitLoci.Hit.Seg[0].ChromID)
	return(1);

// sort by PairReadID ascending
El1ID = pEl1->PairReadID & 0x7fffffff;
El2ID = pEl2->PairReadID & 0x7fffffff;
if(El1ID < El2ID)
	return(-1);
if(El1ID > El2ID)
	return(1);
// same PairReadID, order by 5' followed by 3' read
if(pEl1->PairReadID & 0x80000000)		
	return(1);
return(0);
}


// SortHitmatch
// Sort by ascending read NAR,NumHits(1,0,2,3..), chrom, loci, len, strand, LowMMCnt
int
CKAligner::SortHitMatch(const void *arg1, const void *arg2)
{
tsReadHit *pEl1 = *(tsReadHit **)arg1;
tsReadHit *pEl2 = *(tsReadHit **)arg2;

if(pEl1->NAR < pEl2->NAR)
	return(-1);
if(pEl1->NAR > pEl2->NAR)
	return(1);
if(pEl1->NumHits == 1 && pEl2->NumHits != 1)	// hit counts are to be sorted 1,0,2,3,...
	return(-1);
if(pEl1->NumHits != 1 && pEl2->NumHits == 1)
	return(1);
if(pEl1->NumHits != 1 && pEl2->NumHits != 1)
	{
	if(pEl1->NumHits < pEl2->NumHits)
		return(-1);
	if(pEl1->NumHits > pEl2->NumHits)
		return(1);
	return(0);
	}

if(pEl1->HitLoci.Hit.Seg[0].ChromID < pEl2->HitLoci.Hit.Seg[0].ChromID)
		return(-1);
if(pEl1->HitLoci.Hit.Seg[0].ChromID > pEl2->HitLoci.Hit.Seg[0].ChromID)
	return(1);
if(AdjStartLoci(&pEl1->HitLoci.Hit.Seg[0]) < AdjStartLoci(&pEl2->HitLoci.Hit.Seg[0]))
		return(-1);
if(AdjStartLoci(&pEl1->HitLoci.Hit.Seg[0]) > AdjStartLoci(&pEl2->HitLoci.Hit.Seg[0]))
	return(1);

if(AdjHitLen(&pEl1->HitLoci.Hit.Seg[0]) < AdjHitLen(&pEl2->HitLoci.Hit.Seg[0]))
		return(-1);
if(AdjHitLen(&pEl1->HitLoci.Hit.Seg[0]) > AdjHitLen(&pEl2->HitLoci.Hit.Seg[0]))
	return(1);

if(pEl1->HitLoci.Hit.Seg[0].Strand < pEl2->HitLoci.Hit.Seg[0].Strand)
		return(-1);
if(pEl1->HitLoci.Hit.Seg[0].Strand > pEl2->HitLoci.Hit.Seg[0].Strand)
	return(1);
if(pEl1->LowMMCnt < pEl2->LowMMCnt)
		return(-1);
if(pEl1->LowMMCnt > pEl2->LowMMCnt)
	return(1);
return(0);
}

// SortMultiHits
// Sort by ascending ChromID, HitLoci, HitLen, Hitmismatches, Strand, then ReadID
int
CKAligner::SortMultiHits(const void *arg1, const void *arg2)
{
tsReadHit *pEl1 = (tsReadHit *)arg1;
tsReadHit *pEl2 = (tsReadHit *)arg2;

if(pEl1->HitLoci.Hit.Seg[0].ChromID < pEl2->HitLoci.Hit.Seg[0].ChromID )
		return(-1);
if(pEl1->HitLoci.Hit.Seg[0].ChromID > pEl2->HitLoci.Hit.Seg[0].ChromID )
	return(1);
if(AdjStartLoci(&pEl1->HitLoci.Hit.Seg[0]) < AdjStartLoci(&pEl2->HitLoci.Hit.Seg[0]))
		return(-1);
if(AdjStartLoci(&pEl1->HitLoci.Hit.Seg[0]) > AdjStartLoci(&pEl2->HitLoci.Hit.Seg[0]))
	return(1);
if(AdjHitLen(&pEl1->HitLoci.Hit.Seg[0]) < AdjHitLen(&pEl2->HitLoci.Hit.Seg[0]))
	return(-1);
if(AdjHitLen(&pEl1->HitLoci.Hit.Seg[0]) > AdjHitLen(&pEl2->HitLoci.Hit.Seg[0]))
	return(1);

if(pEl1->HitLoci.Hit.Seg[0].Mismatches < pEl2->HitLoci.Hit.Seg[0].Mismatches )
	return(-1);
if(pEl1->HitLoci.Hit.Seg[0].Mismatches > pEl2->HitLoci.Hit.Seg[0].Mismatches )
	return(1);

if(pEl1->HitLoci.Hit.Seg[0].Strand < pEl2->HitLoci.Hit.Seg[0].Strand )
		return(-1);
if(pEl1->HitLoci.Hit.Seg[0].Strand > pEl2->HitLoci.Hit.Seg[0].Strand )
	return(1);

if(pEl1->ReadID < pEl2->ReadID )
	return(-1);
if(pEl1->ReadID > pEl2->ReadID )
	return(1);
return(0);
}

// SortMultiHitReadIDs
// sort by ascending ReadID then by descending scores
// to break tied scores, further sorted by
int
CKAligner::SortMultiHitReadIDs(const void *arg1, const void *arg2)
{
tsReadHit *pEl1 = (tsReadHit *)arg1;
tsReadHit *pEl2 = (tsReadHit *)arg2;

if(pEl1->ReadID < pEl2->ReadID)
		return(-1);
if(pEl1->ReadID > pEl2->ReadID)
	return(1);

if(pEl1->HitLoci.Hit.Score > pEl2->HitLoci.Hit.Score)
		return(-1);
if(pEl1->HitLoci.Hit.Score < pEl2->HitLoci.Hit.Score)
	return(1);

if(pEl1->HitLoci.Hit.Seg[0].ChromID < pEl2->HitLoci.Hit.Seg[0].ChromID )
		return(-1);
if(pEl1->HitLoci.Hit.Seg[0].ChromID > pEl2->HitLoci.Hit.Seg[0].ChromID )
	return(1);

if(AdjHitLen(&pEl1->HitLoci.Hit.Seg[0]) < AdjHitLen(&pEl2->HitLoci.Hit.Seg[0]))
	return(-1);
if(AdjHitLen(&pEl1->HitLoci.Hit.Seg[0]) > AdjHitLen(&pEl2->HitLoci.Hit.Seg[0]))
	return(1);

if(pEl1->HitLoci.Hit.Seg[0].Mismatches < pEl2->HitLoci.Hit.Seg[0].Mismatches )
	return(-1);
if(pEl1->HitLoci.Hit.Seg[0].Mismatches > pEl2->HitLoci.Hit.Seg[0].Mismatches )
	return(1);

if(AdjStartLoci(&pEl1->HitLoci.Hit.Seg[0]) < AdjStartLoci(&pEl2->HitLoci.Hit.Seg[0]))
		return(-1);
if(AdjStartLoci(&pEl1->HitLoci.Hit.Seg[0]) > AdjStartLoci(&pEl2->HitLoci.Hit.Seg[0]))
	return(1);

if(pEl1->HitLoci.Hit.Seg[0].Strand < pEl2->HitLoci.Hit.Seg[0].Strand )
		return(-1);
if(pEl1->HitLoci.Hit.Seg[0].Strand > pEl2->HitLoci.Hit.Seg[0].Strand )
	return(1);

return(0);
}

// SortSegJuncts
// sort by ascending chrom, start, end
int
CKAligner::SortSegJuncts(const void *arg1, const void *arg2)
{
tsSegJuncts *pEl1 = (tsSegJuncts *)arg1;
tsSegJuncts *pEl2 = (tsSegJuncts *)arg2;

if(pEl1->ChromID < pEl2->ChromID)
		return(-1);
if(pEl1->ChromID > pEl2->ChromID)
	return(1);

if(pEl1->Starts < pEl2->Starts)
		return(-1);
if(pEl1->Starts > pEl2->Starts)
	return(1);
if(pEl1->Ends < pEl2->Ends)
		return(-1);
if(pEl1->Ends > pEl2->Ends)
	return(1);
return(0);
}

int
CKAligner::SortReadSeqs(const void *arg1, const void *arg2)
{
int Idx;
uint8_t *pSeq1;
uint8_t *pSeq2;

tsReadHit *pEl1 = *(tsReadHit **)arg1;
tsReadHit *pEl2 = *(tsReadHit **)arg2;

pSeq1 = &pEl1->Read[pEl1->DescrLen+1];
pSeq2 = &pEl2->Read[pEl2->DescrLen+1];

for(Idx = 0; Idx < pEl1->ReadLen; Idx++, pSeq1++, pSeq2++)
	{
	if((*pSeq1 & 0x07) < (*pSeq2 & 0x07))
		return(-1);
	if((*pSeq1 & 0x07) > (*pSeq2 & 0x07))
		return(1);
	if(Idx >= pEl2->ReadLen)
		return(1);
	}
if(pEl1->ReadLen < pEl2->ReadLen)
	return(-1);
return(0);
}

int
CKAligner::SortLociPValues(const void *arg1, const void *arg2)
{
tsLociPValues *pEl1 = (tsLociPValues *)arg1;
tsLociPValues *pEl2 = (tsLociPValues *)arg2;
if(pEl1->PValue < pEl2->PValue)
		return(-1);
if(pEl1->PValue > pEl2->PValue)
	return(1);
return(0);
}

int
CKAligner::SortPValuesLoci(const void *arg1, const void *arg2)
{
tsLociPValues *pEl1 = (tsLociPValues *)arg1;
tsLociPValues *pEl2 = (tsLociPValues *)arg2;
if(pEl1->Loci < pEl2->Loci)
		return(-1);
if(pEl1->Loci > pEl2->Loci)
	return(1);
return(0);
}


int
CKAligner::SortSiteRelScale(const void *arg1, const void *arg2)
{
tsOctSitePrefs *pEl1 = (tsOctSitePrefs *)arg1;
tsOctSitePrefs *pEl2 = (tsOctSitePrefs *)arg2;
if(pEl1->RelScale < pEl2->RelScale)
		return(-1);
if(pEl1->RelScale > pEl2->RelScale)
	return(1);
return(0);
}

int
CKAligner::SortSiteRelOctamer(const void *arg1, const void *arg2)
{
tsOctSitePrefs *pEl1 = (tsOctSitePrefs *)arg1;
tsOctSitePrefs *pEl2 = (tsOctSitePrefs *)arg2;
if(pEl1->Octamer < pEl2->Octamer)
		return(-1);
if(pEl1->Octamer > pEl2->Octamer)
	return(1);
return(0);
}

int
CKAligner::SortConstraintLoci(const void *arg1, const void *arg2)
{
tsConstraintLoci *pEl1 = (tsConstraintLoci *)arg1;
tsConstraintLoci *pEl2 = (tsConstraintLoci *)arg2;
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

void
CKAligner::AcquireSerialise(void)  // mutually exclusive lock - only a single thread can own this lock
{
#ifdef _WIN32
WaitForSingleObject(m_hMtxIterReads,INFINITE);
#else
pthread_mutex_lock(&m_hMtxIterReads);
#endif
}

void
CKAligner::ReleaseSerialise(void)  // release lock allowing other threads compete for ownership
{
#ifdef _WIN32
ReleaseMutex(m_hMtxIterReads);
#else
pthread_mutex_unlock(&m_hMtxIterReads);
#endif
}

void
CKAligner::AcquireSerialiseMH(void)
{
#ifdef _WIN32
WaitForSingleObject(m_hMtxMHReads,INFINITE);
#else
pthread_mutex_lock(&m_hMtxMHReads);
#endif
}

void
CKAligner::ReleaseSerialiseMH(void)
{
#ifdef _WIN32
ReleaseMutex(m_hMtxMHReads);
#else
pthread_mutex_unlock(&m_hMtxMHReads);
#endif
}

void
CKAligner::AcquireSharedLock(void)  // acquire SRW lock ownership in shared mode with other threads
{
#ifdef _WIN32
AcquireSRWLockShared(&m_hRwLock); // Acquires an SRW lock in shared mode for calling thread
#else
pthread_rwlock_rdlock(&m_hRwLock);
#endif
}

void
CKAligner::AcquireExclusiveLock(void)   // acquire exclusive R/W lock for calling thread
{
int Rslt = 0;
int SpinCnt = 500;
int BackoffMS = 5;

while(1)
	{
	AcquireSerialise();     // serialise access to locks
#ifdef _WIN32
	if(TryAcquireSRWLockExclusive(&m_hRwLock))   // attempt to acquire an exclusive lock
		break;
#else
	Rslt = pthread_rwlock_trywrlock(&m_hRwLock);
	if(Rslt == 0 || Rslt == EDEADLK)
		break;
#endif
	ReleaseSerialise();   // allow other threads a chance to release shared or their exclusive lock
	if(SpinCnt -= 1)      
		continue;
	CUtility::SleepMillisecs(BackoffMS);
	SpinCnt = 100;
	if(BackoffMS < 500)
		BackoffMS += 2;
	}
}

void
CKAligner::ReleaseExclusiveLock(void)   // release exclusive R/W lock previously acquired by calling thread
{
#ifdef _WIN32
ReleaseSRWLockExclusive(&m_hRwLock); // Releases SRW lock from exclusive mode
#else
pthread_rwlock_unlock(&m_hRwLock);
#endif
ReleaseSerialise();   // allow other threads to gain shared or exclusive locks
}

void
CKAligner::ReleaseSharedLock(void)  // release SRW shared lock ownership
{
#ifdef _WIN32
ReleaseSRWLockShared(&m_hRwLock); // Releases SRW lock from shared mode for calling thread
#else
pthread_rwlock_unlock(&m_hRwLock);
#endif
}

int
CKAligner::ProcLoadReadFiles(tsLoadReadsThreadPars *pPars)
{
teBSFrsltCodes Rslt;
int *pRslt = pPars->pRslt;
int Idx;
char *pszInfile;

int NumInputFilesProcessed;

AcquireExclusiveLock();
m_bAllReadsLoaded = false;
m_LoadReadsRslt = eBSFSuccess;
ReleaseExclusiveLock();

// first try to load as pre-processed reads '.rds' (as generated by 'simreads' or 'genreads')
// if unable to load as a '.rds' then try to load as raw fasta/fastq

if(((Rslt = (teBSFrsltCodes)LoadReads(m_ppszPE1InputFiles[0])) < eBSFSuccess) && Rslt != eBSFerrNotBioseq )
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to load reads from %s",m_ppszPE1InputFiles[0]);
	AcquireExclusiveLock();
	*pRslt = Rslt;
	m_bAllReadsLoaded = true;
	ReleaseExclusiveLock();
	return(Rslt);
	}

if(Rslt != eBSFerrNotBioseq)
	{
	if(Rslt == 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"There were no reads in %s",m_ppszPE1InputFiles[0]);
		Rslt = eBSFerrNoEntries;
		}
	else
		Rslt = eBSFSuccess;
	AcquireExclusiveLock();
	*pRslt = Rslt;
	m_bAllReadsLoaded = true;
	ReleaseExclusiveLock();
	return(Rslt);
	}

if(m_hInFile != -1)
	{
	close(m_hInFile);
	m_hInFile = -1;
	}

AcquireExclusiveLock();
m_bAllReadsLoaded = false;
memset(&m_FileHdr,0,sizeof(tsBSFRdsHdr));
m_FileHdr.Magic[0] = 'b'; m_FileHdr.Magic[1] = 'i'; m_FileHdr.Magic[2] = 'o'; m_FileHdr.Magic[3] = 'r';
m_FileHdr.Version = cBSFRdsVersion;
m_FileHdr.FlagsK = 1;
m_FileHdr.FlagsCS = m_bIsSOLiD;
m_FileHdr.FlagsPR = m_PEproc == ePEdefault ? 0 : 1;
m_FileHdr.PMode = (uint8_t)0;
m_FileHdr.QMode = m_QMethod;
m_FileHdr.Trim5 = m_Trim5;
m_FileHdr.Trim3 = m_Trim3;
m_pReadHits = nullptr;
m_AllocdReadHitsMem = 0;
m_NumReadsLoaded = 0;
m_UsedReadHitsMem = 0;
m_FinalReadID = 0;
m_NumDescrReads = 0;
ReleaseExclusiveLock();

if(m_FileHdr.FlagsPR && (m_NumPE1InputFiles < 1 || (m_NumPE1InputFiles != m_NumPE2InputFiles)))
	{
	Rslt = eBSFerrParams;
	AcquireExclusiveLock();
	*pRslt = Rslt;
	m_bAllReadsLoaded = true;
	ReleaseExclusiveLock();
	return(Rslt);
	}

NumInputFilesProcessed = 0;
if(!m_FileHdr.FlagsPR)
	{
	CSimpleGlob glob(SG_GLOB_FULLSORT);
	for(Idx = 0; Idx < m_NumPE1InputFiles; Idx++)
		{
		glob.Init();
		if(glob.Add(m_ppszPE1InputFiles[Idx]) < SG_SUCCESS)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to glob '%s",m_ppszPE1InputFiles[Idx]);
			Rslt = eBSFerrOpnFile;
			AcquireExclusiveLock();
			*pRslt = Rslt;
			m_bAllReadsLoaded = true;
			m_LoadReadsRslt = Rslt;
			ReleaseExclusiveLock();
			return(Rslt);	// treat as though unable to open file
			}

		if(glob.FileCount() <= 0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to locate any source raw reads file matching '%s",m_ppszPE1InputFiles[Idx]);
			continue;
			}

		Rslt = eBSFSuccess;
		for (int FileID = 0; Rslt >= eBSFSuccess &&  FileID < glob.FileCount(); ++FileID)
			{
			pszInfile = glob.File(FileID);

			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading and parsing reads from raw sequence file '%s'\n",pszInfile);
			Rslt = LoadRawReads(false,NumInputFilesProcessed+1,pszInfile,nullptr);
			if(Rslt != eBSFSuccess)
				{
				if(m_TermBackgoundThreads == 0)
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Load failed for input raw sequence file '%s'\n",pszInfile);
				AcquireExclusiveLock();
				*pRslt = Rslt;
				m_bAllReadsLoaded = true;
				m_LoadReadsRslt = Rslt;
				ReleaseExclusiveLock();
				return(Rslt);
				}
			NumInputFilesProcessed += 1;
			}
		}
	}
else
	{
	Rslt = eBSFSuccess;
	NumInputFilesProcessed = 0;
	for(int FileID = 0; Rslt >= eBSFSuccess &&  FileID < m_NumPE1InputFiles; ++FileID)
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading and parsing reads from 5' PE1 raw sequence file '%s'",m_ppszPE1InputFiles[FileID]);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading and parsing reads from 3' PE2 raw sequence file '%s'",m_ppszPE2InputFiles[FileID]);

		Rslt = LoadRawReads(true,NumInputFilesProcessed+1,m_ppszPE1InputFiles[FileID],m_ppszPE2InputFiles[FileID]);
		if(Rslt != eBSFSuccess)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Load failed for input raw sequence PE files '%s' and '%s'\n",m_ppszPE1InputFiles[FileID],m_ppszPE2InputFiles[FileID]);
			AcquireExclusiveLock();
			*pRslt = Rslt;
			m_bAllReadsLoaded = true;
			m_LoadReadsRslt = Rslt;
			ReleaseExclusiveLock();
			return(Rslt);
			}
		NumInputFilesProcessed += 2;
		}
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"%d raw sequence file%s parsed and reads loaded for aligning", NumInputFilesProcessed, NumInputFilesProcessed == 1 ? " was" : "s were");
if(NumInputFilesProcessed == 0)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Nothing to do, no raw sequence files to be filtered");
	AcquireExclusiveLock();
	*pRslt = Rslt;
	m_bAllReadsLoaded = true;
	m_LoadReadsRslt = eBSFerrNoEntries;
	ReleaseExclusiveLock();
	return(Rslt);
	}

// can now update header with final numbers
AcquireExclusiveLock();
m_FileHdr.NumRds = m_NumDescrReads;
m_FileHdr.OrigNumReads = m_NumDescrReads;
m_FileHdr.TotReadsLen = m_DataBuffOfs;
m_FinalReadID = m_NumDescrReads;
m_NumReadsLoaded = m_NumDescrReads;
m_LoadReadsRslt = m_NumReadsLoaded > 0 ? eBSFSuccess : eBSFerrNoEntries;
m_bAllReadsLoaded = true;
*pRslt = Rslt;
ReleaseExclusiveLock();
return(Rslt);
}



int
CKAligner::AddEntry(bool bPEProcessing,	// true if entries are from a PE dataset, false if SE dataset
			bool bIsPairRead,	// false if SE or PE1, true if this is the paired read PE2
		 uint32_t PairReadID,		// identifies partner of this read if paired read processing (0 if no partner read)
		 uint8_t FileID,			// identifies file from which this read was parsed
		 int DescrLen,			// length of following descriptor
		 char *pszReadDescr,	// copy of descriptor, used to pair reads with matching descriptors
		 int ReadLen,			// length of following read
		 uint8_t *pszReadBuff)	// packed read + phred score
{
uint8_t *pTmpAlloc;
tsReadHit *pReadHit;
size_t memreq;

if(m_pReadHits == nullptr)
	{
	memreq = cDataBuffAlloc;
	AcquireExclusiveLock();
#ifdef _WIN32
	m_pReadHits = (tsReadHit *) malloc((size_t)memreq);
	if(m_pReadHits == nullptr)
		{
		ReleaseExclusiveLock();
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddEntry: Memory allocation of %zd bytes failed",(int64_t)memreq);
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	m_pReadHits = (tsReadHit *)mmap(nullptr,(size_t)memreq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pReadHits == MAP_FAILED)
		{
		m_pReadHits = nullptr;
		ReleaseExclusiveLock();
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddEntry: Memory allocation of %zd bytes through mmap()  failed",(int64_t)memreq,strerror(errno));
		return(eBSFerrMem);
		}
#endif
	m_AllocdReadHitsMem = memreq;
	m_DataBuffOfs = 0;
	ReleaseExclusiveLock();
	}

// need to allocate more memory? NOTE: allowing margin of 16K
if((m_AllocdReadHitsMem - m_DataBuffOfs) < (sizeof(tsReadHit) +  ReadLen + DescrLen + 0x03fff))
	{
	memreq = m_AllocdReadHitsMem + cDataBuffAlloc;
	AcquireExclusiveLock();
#ifdef _WIN32
	pTmpAlloc = (uint8_t *) realloc(m_pReadHits,memreq);
	if(pTmpAlloc == nullptr)
		{
#else
	pTmpAlloc = (uint8_t *)mremap(m_pReadHits,m_AllocdReadHitsMem,memreq,MREMAP_MAYMOVE);
	if(pTmpAlloc == MAP_FAILED)
		{
		pTmpAlloc = nullptr;
#endif

		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddEntry: Memory reallocation to %zd bytes failed - %s",memreq,strerror(errno));
		ReleaseExclusiveLock();
		return(eBSFerrMem);
		}
	m_pReadHits = (tsReadHit *)pTmpAlloc;
	m_AllocdReadHitsMem = memreq;
	ReleaseExclusiveLock();
	}

pReadHit = (tsReadHit *)((uint8_t *)m_pReadHits + m_DataBuffOfs);
m_DataBuffOfs += sizeof(tsReadHit) + ReadLen + DescrLen;

memset(pReadHit,0,sizeof(tsReadHit));
pReadHit->PrevSizeOf = m_PrevSizeOf;
pReadHit->ReadID = ++m_NumDescrReads;
pReadHit->PairReadID = PairReadID;
if(bIsPairRead)
	pReadHit->PairReadID |= 0x80000000;
pReadHit->NumReads = 1;
pReadHit->HitLoci.Hit.Seg[0].Strand = '?';
pReadHit->ReadLen = ReadLen;
pReadHit->DescrLen = (uint8_t)DescrLen;
if(DescrLen > 0)
	memcpy((char *)&pReadHit->Read[0],pszReadDescr,DescrLen+1);
else
	pReadHit->Read[0] = '\0';
memmove(&pReadHit->Read[DescrLen+1],pszReadBuff,ReadLen);
m_PrevSizeOf = (uint32_t)sizeof(tsReadHit) + ReadLen + DescrLen;

// processing threads are only updated with actual number of loaded reads every 250K reads so as
// to minimise disruption to the actual aligner threads which will also be serialised through m_hMtxIterReads
uint32_t RptDiff = 250000;
if(m_SampleNthRawRead > 1)
	RptDiff = 1 + (RptDiff/m_SampleNthRawRead);
  
if((!bPEProcessing || (bPEProcessing && bIsPairRead)) && m_NumDescrReads > 0 && (m_NumDescrReads - m_NumReadsLoaded) >= RptDiff)
	{
	AcquireSerialise();
	m_FinalReadID = m_NumDescrReads;
	m_NumReadsLoaded = m_NumDescrReads;
	ReleaseSerialise();
	}
return(eBSFSuccess);
}


double												// returned prob of read being error free
CKAligner::GenProbErrFreeRead(int QSSchema,			// guestimated scoring schema - 0: no scoring, 1: Solexa, 2: Illumina 1.3+, 3: Illumina 1.5+, 4: Illumina 1.8+ or could be Sanger 
				   int ReadLen,						// read length
				   uint8_t *pQScores)					// Phred quality scores
{
int SeqOfs;
int Score;
double ProbErr;
double ProbNoReadErr;
if(QSSchema == 0 || ReadLen < 1 || pQScores == nullptr || pQScores[0] == '\0')		// if can't score then have to assume the best - read is sequencing base call error free
	return(1.0);

ProbNoReadErr = 1.0;
for (SeqOfs = 0; SeqOfs < ReadLen; SeqOfs++, pQScores++)
	{
	switch(QSSchema) { // quality scoring schema 
		case 4: // Illumina 1.8+ or Sanger. Sanger is '!' (0) to 'J' (41) and Illumina is '#' (2) to 'J' (41)
			Score = (int)(*pQScores - (uint8_t)'!');
			break;
		default:		// all other scoring schemas have 40 at 'h'
			if(*pQScores <= '@')
				Score = 0;
			else
				Score = (int)(*pQScores - (uint8_t)'@');
			break;
			}
	if(Score < 0)			// force scores to always be in the range 0..41 --- shouldn't be outside 0..41 but some qscores have been observed to be outside the expected range!
		Score = 0;
	else
		if(Score > 41)
			Score = 41;

	ProbErr = 1.0 / pow(10.0,(double)Score/10.0);
	ProbNoReadErr *= (1.0 - ProbErr); 
	}
if(ProbNoReadErr < 0.0)
	ProbNoReadErr = 0.0;
else
	if(ProbNoReadErr > 1.0)
		ProbNoReadErr = 1.0;
return(ProbNoReadErr);
}

teBSFrsltCodes
CKAligner::LoadRawReads(bool bIsPairReads,	// true if paired end processing - PE1 reads in pszPE1File and PE2 reads in pszPE2File
		  int FileID,						// uniquely identifies source file for PE1, FileID + 1 uniquely identifies PE2 file
		  char *pszPE1File,					// process PE1 reads from this file
		  char *pszPE2File)					// optionally process PE2 reads from this file
{
static int FileNamesOfs = 0;
teBSFrsltCodes Rslt;
int Idx;
bool bIsFastq;

uint8_t *pReadBuff;
uint8_t *pQualBuff;
uint8_t Qphred;

int PE1NumDescrReads;
int PE1DescrLen;
uint8_t szPE1DescrBuff[cMaxKADescrLen+1];
int PE1ReadLen;
int PE1QualLen;
uint8_t* pszPE1ReadBuff = nullptr;
uint8_t *pszPE1QualBuff = nullptr;
int PE1NumReadsAccepted;
int PE1NumInvalValues;
int PE1NumUnsupportedBases;
int PE1NumUnderlength;
int PE1NumOverlength;

int PE2NumDescrReads;
int PE2DescrLen;
uint8_t szPE2DescrBuff[cMaxKADescrLen+1];
int PE2ReadLen;

int PE2QualLen;
uint8_t* pszPE2ReadBuff = nullptr;
uint8_t *pszPE2QualBuff = nullptr;
int PE2NumReadsAccepted;
int PE2NumInvalValues;
int PE2NumUnsupportedBases;
int PE2NumUnderlength;
int PE2NumOverlength;

int ContamLen5PE1;
int ContamLen3PE1;
int ContamLen5PE2;
int ContamLen3PE2;

int NumContamLen5PE1;
int NumContamLen3PE1;
int NumContamLen5PE2;
int NumContamLen3PE2;

uint32_t NxtToSample;

CFasta PE1Fasta;
CFasta PE2Fasta;

int32_t EstScoreSchema;				// guestimated scoring schema - 0: no scoring, 1: Solexa, 2: Illumina 1.3+, 3: Illumina 1.5+, 4: Illumina 1.8+ or could be Sanger 
int32_t EstSeqLen;
int32_t EstDescrLen;
uint32_t EstNumSeqs;
int64_t ReqAllocSize;

// try and guestimate memory requirements so these can be allocated upfront - will be realloc'd if guestimate is wrong
if((EstNumSeqs = (teBSFrsltCodes)PE1Fasta.FastaEstSizes(pszPE1File,nullptr,nullptr,&EstDescrLen,nullptr,&EstSeqLen,&EstScoreSchema)) == 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to estimate memory requirements for file '%s'",pszPE1File);
	return(eBSFerrOpnFile);
	}

if(m_SampleNthRawRead > 1)
	EstNumSeqs = 1 + (EstNumSeqs / m_SampleNthRawRead);	// ensure at least 1 read will be sampled!

// ReqAllocSize assumes -
// a) EstNumSeqs is not accurate so adds another 100000 to reduce chance of subsequent realloc required
// b) descriptors will be trimmed to 1st whitespace
ReqAllocSize = (int64_t)(EstNumSeqs + 100000) * (int64_t)(sizeof(tsReadHit) + EstSeqLen + max(20,(EstDescrLen/2)));	
if(bIsPairReads)
	ReqAllocSize *= 2;

uint32_t PairReadID;

PairReadID = (m_NumDescrReads/2) + 1;				// if bIsPairReads then start paired reads identifiers from this value and increment after each read processed
if((Rslt=(teBSFrsltCodes)PE1Fasta.Open(pszPE1File,true))!=eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Load: Unable to open '%s' [%s] %s",pszPE1File,PE1Fasta.ErrText((teBSFrsltCodes)Rslt),PE1Fasta.GetErrMsg());
	return(Rslt);
	}
if(!m_FileHdr.NumFiles)
	FileNamesOfs = 0;

m_FileHdr.NumFiles += 1;
if((FileNamesOfs + strlen(pszPE1File) + 1) < sizeof(m_FileHdr.FileNames))
	{
	strcpy((char *)&m_FileHdr.FileNames[FileNamesOfs],pszPE1File);
	FileNamesOfs += (int)strlen(pszPE1File) + 1;
	}

if(m_bIsSOLiD == true && m_bIsSOLiD != PE1Fasta.IsSOLiD())
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Load: Aligning for %s but reads file '%s' is in %s",m_bIsSOLiD ? "Colorspace" : "Basespace",pszPE1File,m_bIsSOLiD ? "Basespace" : "Colorspace");
	PE1Fasta.Close();
	return(eBSFerrCvrtType);
	}

bIsFastq = PE1Fasta.IsFastq();

if(bIsPairReads)	
	{
	if((Rslt=(teBSFrsltCodes)PE2Fasta.Open(pszPE2File,true))!=eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Load: Unable to open '%s' [%s] %s",pszPE2File,PE2Fasta.ErrText((teBSFrsltCodes)Rslt),PE2Fasta.GetErrMsg());
		PE1Fasta.Close();
		return(Rslt);
		}

	m_FileHdr.NumFiles += 1;
	if((FileNamesOfs + strlen(pszPE2File) + 1) < sizeof(m_FileHdr.FileNames))
		{
		strcpy((char *)&m_FileHdr.FileNames[FileNamesOfs],pszPE2File);
		FileNamesOfs += (int)strlen(pszPE2File) + 1;
		}

	if (m_bIsSOLiD == true && m_bIsSOLiD != PE1Fasta.IsSOLiD())
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Load: Aligning for %s but reads file '%s' is in %s",m_bIsSOLiD ? "Colorspace" : "Basespace",pszPE2File,m_bIsSOLiD ? "Basespace" : "Colorspace");
		PE1Fasta.Close();
		PE2Fasta.Close();
		return(eBSFerrCvrtType);
		}

	if(bIsFastq != PE2Fasta.IsFastq())
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Load: Paired end files not of same type");
		PE1Fasta.Close();
		PE2Fasta.Close();
		return(eBSFerrCvrtType);
		}
	}

AcquireExclusiveLock();
if(m_pReadHits == nullptr)
	{
	m_AllocdReadHitsMem = (size_t)ReqAllocSize;
#ifdef _WIN32
	m_pReadHits = (tsReadHit *) malloc((size_t)m_AllocdReadHitsMem);
	if(m_pReadHits == nullptr)
		{
		ReleaseExclusiveLock();
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Concatenated sequences memory allocation of %zd bytes - %s",(int64_t)m_AllocdReadHitsMem,strerror(errno));
		PE1Fasta.Close();
		m_AllocdReadHitsMem = 0;
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	m_pReadHits = (tsReadHit *)mmap(nullptr,m_AllocdReadHitsMem, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pReadHits == MAP_FAILED)
		{
		ReleaseExclusiveLock();
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Concatenated sequences memory of %zd bytes through mmap()  failed - %s",(int64_t)m_AllocdReadHitsMem,strerror(errno));
		PE1Fasta.Close();
		m_pReadHits = nullptr;
		return(eBSFerrMem);
		}
#endif
	m_DataBuffOfs = 0;

	}
else
	{
	uint8_t *pDstSeq;
	size_t memreq;
	if((m_DataBuffOfs + ReqAllocSize + 0x0fffff) >= m_AllocdReadHitsMem)		// 1M as a small safety margin!
		{
		memreq = (size_t)(m_AllocdReadHitsMem + ReqAllocSize);
#ifdef _WIN32
		pDstSeq = (uint8_t *) realloc(m_pReadHits,memreq);
#else
		pDstSeq = (uint8_t *)mremap(m_pReadHits,m_AllocdReadHitsMem,memreq,MREMAP_MAYMOVE);
		if(pDstSeq == MAP_FAILED)
			pDstSeq = nullptr;
#endif
		if(pDstSeq == nullptr)
			{
			ReleaseExclusiveLock();
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddSeq: Memory re-allocation to %zd bytes - %s",memreq,strerror(errno));
			PE1Fasta.Close();
			return(eBSFerrMem);
			}
		m_AllocdReadHitsMem = memreq;
		m_pReadHits = (tsReadHit *)pDstSeq;
		}
	}
ReleaseExclusiveLock();

PE1NumUnsupportedBases = 0;
PE1NumDescrReads = 0;
PE1NumReadsAccepted = 0;
PE1NumInvalValues = 0;
PE1NumUnderlength = 0;
PE1NumOverlength = 0;
PE2NumUnsupportedBases = 0;
PE2NumDescrReads = 0;
PE2NumReadsAccepted = 0;
PE2NumInvalValues = 0;
PE2NumUnderlength = 0;
PE2NumOverlength = 0;

NumContamLen5PE1 = 0;
NumContamLen3PE1 = 0;
NumContamLen5PE2 = 0;
NumContamLen3PE2 = 0;

NxtToSample = m_SampleNthRawRead;
pszPE1ReadBuff = new uint8_t[cMaxReadLen+1000];
pszPE1QualBuff = new uint8_t[cMaxReadLen+1000];
if (bIsPairReads)
	{
	pszPE2ReadBuff = new uint8_t[cMaxReadLen + 1000];
	pszPE2QualBuff = new uint8_t[cMaxReadLen + 1000];
	}

while((Rslt = (teBSFrsltCodes)(PE1ReadLen = PE1Fasta.ReadSequence(pszPE1ReadBuff, cMaxReadLen-1,true,false))) > eBSFSuccess)
	{
	if(m_TermBackgoundThreads != 0)	// need to immediately self-terminate?
		break;

	PE1NumDescrReads += 1;
	if(PE1ReadLen == eBSFFastaDescr)		// just read a descriptor line which would be as expected for multifasta or fastq
		{
		PE1DescrLen = PE1Fasta.ReadDescriptor((char *)szPE1DescrBuff,sizeof(szPE1DescrBuff)-1);
		szPE1DescrBuff[cMaxKADescrLen-1] = '\0';
		PE1ReadLen = PE1Fasta.ReadSequence(pszPE1ReadBuff, cMaxReadLen-1); // an observed issue is that some fastq triming toolkits can trim out the complete sequence leaving the descriptor imeediately followed by another descriptor with no sequence present 
		if(PE1ReadLen < 0 || PE1ReadLen == eBSFFastaDescr || PE1ReadLen >= cMaxReadLen)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Problem parsing PE1 sequence after %d reads parsed (pre-filtering) from file: '%s'",PE1NumDescrReads,pszPE1File);
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Last PE1 descriptor parsed: %s",szPE1DescrBuff);
			if(Rslt < eBSFSuccess)
				{
				while(PE1Fasta.NumErrMsgs())
					gDiagnostics.DiagOut(eDLFatal, gszProcName, PE1Fasta.GetErrMsg());
				}
			PE1Fasta.Close();
			if(bIsPairReads)
				PE2Fasta.Close();
			if (pszPE1ReadBuff != nullptr)
				delete []pszPE1ReadBuff;
			if (pszPE2ReadBuff != nullptr)
				delete[]pszPE2ReadBuff;
			if (pszPE1QualBuff != nullptr)
				delete []pszPE1QualBuff;
			if (pszPE2QualBuff != nullptr)
				delete[]pszPE2QualBuff;
			return(eBSFerrParse);
			}

		// if paired end processing then also load PE2 read
		if(bIsPairReads)
			{
			if((Rslt = (teBSFrsltCodes)(PE2ReadLen = PE2Fasta.ReadSequence(pszPE2ReadBuff, cMaxReadLen-1,true,false))) <= eBSFSuccess)
				{
				gDiagnostics.DiagOut(eDLWarn,gszProcName,"Successfully parsed %d PE1 descriptors (pre-filtering) from file: '%s'",PE1NumDescrReads,pszPE1File);
				gDiagnostics.DiagOut(eDLWarn,gszProcName,"Only able to parse corresponding %d PE2 descriptors (pre-filtering) from file: '%s'",PE2NumDescrReads,pszPE2File);
				if(PE2ReadLen == 0)
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected matching numbers of PE1 and PE2 (pre-filtering) descriptors");
				else
					{
					if(Rslt < eBSFSuccess)
						{
						while(PE2Fasta.NumErrMsgs())
							gDiagnostics.DiagOut(eDLFatal, gszProcName, PE2Fasta.GetErrMsg());
						}
					}
				if(PE1NumDescrReads)
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"Last PE1 descriptor parsed: %s",szPE1DescrBuff);
				if(PE2NumDescrReads)
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"Last PE2 descriptor parsed: %s",szPE2DescrBuff);
				PE1Fasta.Close();
				PE2Fasta.Close();
				if (pszPE1ReadBuff != nullptr)
					delete[]pszPE1ReadBuff;
				if (pszPE2ReadBuff != nullptr)
					delete[]pszPE2ReadBuff;
				if (pszPE1QualBuff != nullptr)
					delete[]pszPE1QualBuff;
				if (pszPE2QualBuff != nullptr)
					delete[]pszPE2QualBuff;
				return(eBSFerrParse);
				}

			if(m_TermBackgoundThreads != 0)	// need to immediately self-terminate?
				break;

			PE2NumDescrReads += 1;
			if(PE2ReadLen == eBSFFastaDescr)		// just read a descriptor line which would be as expected for multifasta or fastq
				{
				PE2DescrLen = PE2Fasta.ReadDescriptor((char *)szPE2DescrBuff,sizeof(szPE2DescrBuff)-1);
				szPE2DescrBuff[cMaxKADescrLen-1] = '\0';
				PE2ReadLen = PE2Fasta.ReadSequence(pszPE2ReadBuff, cMaxReadLen-1);
				if(PE2ReadLen < 0  || PE1ReadLen == eBSFFastaDescr  || PE1ReadLen >= cMaxReadLen)
					{
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"Problem PE2 parsing sequence after %d reads parsed (pre-filtering) from file: '%s'",PE2NumDescrReads,pszPE2File);
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"Last descriptor parsed: %s",szPE2DescrBuff);
					PE1Fasta.Close();
					PE2Fasta.Close();
					if (pszPE1ReadBuff != nullptr)
						delete[]pszPE1ReadBuff;
					if (pszPE2ReadBuff != nullptr)
						delete[]pszPE2ReadBuff;
					if (pszPE1QualBuff != nullptr)
						delete[]pszPE1QualBuff;
					if (pszPE2QualBuff != nullptr)
						delete[]pszPE2QualBuff;
					return(eBSFerrParse);
					}
				}
			else
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"raw sequence file '%s' processing error: %s ",pszPE2File,
													Rslt > 0 ? "over length read" : "not a multifasta short reads or fastq file");
				PE1Fasta.Close();
				if(bIsPairReads)
					PE2Fasta.Close();
				if (pszPE1ReadBuff != nullptr)
					delete[]pszPE1ReadBuff;
				if (pszPE2ReadBuff != nullptr)
					delete[]pszPE2ReadBuff;
				if (pszPE1QualBuff != nullptr)
					delete[]pszPE1QualBuff;
				if (pszPE2QualBuff != nullptr)
					delete[]pszPE2QualBuff;
				return(eBSFerrParse);
				}
			}

		if(m_SampleNthRawRead > 1)
			{
			NxtToSample += 1;
			if(m_SampleNthRawRead > NxtToSample)
				continue;
			NxtToSample = 0;
			}


		if(m_pContaminants != nullptr)
			{
			// currently treating any contaminant match errors as if simply there was no overlap - should really report these!!!
			if((ContamLen5PE1 = m_pContaminants->MatchContaminants(eAOF5PE1Targ,1,m_Trim5+1,PE1ReadLen,pszPE1ReadBuff)) <= m_Trim5)
				ContamLen5PE1 = 0;
			else
				{
				if(ContamLen5PE1 > m_Trim5)
					ContamLen5PE1 -= m_Trim5;
				}
			if((ContamLen3PE1 = m_pContaminants->MatchContaminants(eAOF3PE1Targ,1,m_Trim3+1,PE1ReadLen,pszPE1ReadBuff)) <= m_Trim3)
				ContamLen3PE1 = 0;
			else
				{
				if(ContamLen3PE1 > m_Trim3)
					ContamLen3PE1 -= m_Trim3;
				}

			if(bIsPairReads)
				{
				if((ContamLen5PE2 = m_pContaminants->MatchContaminants(eAOF5PE2Targ,1,m_Trim5+1,PE2ReadLen,pszPE2ReadBuff)) <= m_Trim5)
					ContamLen5PE2 = 0;
				else
					{
					if(ContamLen5PE2 > m_Trim5)
						ContamLen5PE2 -= m_Trim5;
					}
				if((ContamLen3PE2 = m_pContaminants->MatchContaminants(eAOF3PE2Targ,1,m_Trim3+1,PE2ReadLen,pszPE2ReadBuff)) <= m_Trim3)
					ContamLen3PE2 = 0;
				else
					{
					if(ContamLen3PE2 > m_Trim3)
						ContamLen3PE2 -= m_Trim3;
					}
				}
			else
				{
				ContamLen5PE2 = 0;
				ContamLen3PE2 = 0;
				}
			}
		else
			{
			ContamLen5PE1 = 0;
			ContamLen3PE1 = 0;
			ContamLen5PE2 = 0;
			ContamLen3PE2 = 0;
			}

		// ensure would still have a sequence of at least m_MinAcceptReadLen after any end trims were applied
		if((m_Trim5 + m_Trim3 + ContamLen5PE1 + ContamLen3PE1 + (int)m_MinAcceptReadLen) > PE1ReadLen)
			{
			PE1NumUnderlength += 1;
			if(PE1NumUnderlength <= 10)
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Load: under length (%d) sequence in '%s' after end trims has been sloughed..",PE1ReadLen,pszPE1File);
			continue;
			}

		// ensure would have a sequence of no more than m_MaxAcceptReadLen after any end trims were applied
		if((m_Trim5 + m_Trim3 + ContamLen5PE1 + ContamLen3PE1 + (int)m_MaxAcceptReadLen) < PE1ReadLen)
			{
			PE1NumOverlength += 1;
			if(PE1NumOverlength <= 10)
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Load: over length (%d) sequence in '%s' after end trims has been sloughed..",PE1ReadLen,pszPE1File);
			continue;
			}

		if(bIsPairReads)
			{
			if((m_Trim5 + m_Trim3 + ContamLen5PE2 + ContamLen3PE2  + (int)m_MinAcceptReadLen) > PE2ReadLen)
				{
				PE2NumUnderlength += 1;
				if(PE2NumUnderlength <= 10)
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"Load: under length (%d) sequence in '%s' after end trims has been sloughed..",PE2ReadLen,pszPE2File);
				continue;
				}
			if((m_Trim5 + m_Trim3 + ContamLen5PE2 + ContamLen3PE2  + (int)m_MaxAcceptReadLen) < PE2ReadLen)
				{
				PE2NumOverlength += 1;
				if(PE2NumOverlength <= 10)
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"Load: over length (%d) sequence in '%s' after end trims has been sloughed..",PE2ReadLen,pszPE2File);
				continue;
				}
			}

		if(bIsFastq && m_QMethod != eFQIgnore)
			{
			PE1QualLen = PE1Fasta.ReadQValues((char *)pszPE1QualBuff, cMaxReadLen-1);
			if(PE1QualLen != PE1ReadLen)		// must be same...
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Load: quality length (%d) not same as read length (%d) for '%s' entry in file '%s'",PE1QualLen,PE1ReadLen,szPE1DescrBuff,pszPE1File);
				PE1Fasta.Close();
				if(bIsPairReads)
					PE2Fasta.Close();
				if (pszPE1ReadBuff != nullptr)
					delete[]pszPE1ReadBuff;
				if (pszPE2ReadBuff != nullptr)
					delete[]pszPE2ReadBuff;
				if (pszPE1QualBuff != nullptr)
					delete[]pszPE1QualBuff;
				if (pszPE2QualBuff != nullptr)
					delete[]pszPE2QualBuff;
				return(eBSFerrParse);
				}
			// normalise the quality score to be in range 0..15 (needs to fit into 4 bits!)
			pQualBuff = pszPE1QualBuff;
			for(Idx = 0; Idx < PE1ReadLen; Idx++,pQualBuff++)
				{
				switch(m_QMethod) {
					case eFQIgnore:		// simply treat as the minimum phred
						Qphred = 0;
						break;

					case eFQSanger:		// Qphred = -10 log10(P), where Qphred is in range 0..93; Illumina 1.8+ is essentially the same as Sanger
						if(*pQualBuff < 33 || *pQualBuff >= 126)
							{
							if(!PE1NumInvalValues++)
								gDiagnostics.DiagOut(eDLWarn,gszProcName,"Load: quality value of '%c' outside range expected in read '%s' from file '%s' for Sanger",*(char *)pQualBuff,szPE1DescrBuff,pszPE1File);
							if(*pQualBuff < 33)
								*pQualBuff = 33;
							else
								*pQualBuff = 125;
							}
						Qphred = *pQualBuff - 33;	// Sanger encodes into ascii starting from decimal 33 '!'
						if(Qphred > 40)				// clamp at phred equiv to 0.0001
							Qphred = 40;
						break;

					case eFQIllumia:	//Qphred = -10 log10(P), where Qphred is in range 0..63
						if(*pQualBuff < 64 || *pQualBuff >= 126)
							{
							if(!PE1NumInvalValues++)
								gDiagnostics.DiagOut(eDLWarn,gszProcName,"Load: quality value of '%c' outside range expected in read '%s' from file '%s' for Illumina 1.3+",*(char *)pQualBuff,szPE1DescrBuff,pszPE1File);
							if(*pQualBuff < 64)
								*pQualBuff = 64;
							else
								*pQualBuff = 125;
							}

						Qphred = *pQualBuff - 64;	//Illumia encodes into ascii starting from decimal 64
						if(Qphred > 40)				// clamp at phred equiv to 0.0001
							Qphred = 40;
						break;

					case eFQSolexa:		// SolexaQ = -10 log10(P/(1-P)), where SolexaQ is in range -5 to 62, note the negative value
										// negative values will result if P > 0.5
										// $Q = 10 * log(1 + 10 ** (ord(SolexaQphred) - 64) / 10.0)) / log(10);
										// once Qphred is over about 15 then essentially same as Sanger and Illumina 1.3+ so
										// is it worth doing the full conversion????
						if(*pQualBuff < 59 || *pQualBuff >= 126)
							{
							if(!PE1NumInvalValues++)
								gDiagnostics.DiagOut(eDLWarn,gszProcName,"Load: quality value of '%c' outside range expected in read '%s' from file '%s' for Solexa/Illumina pre 1.3",*(char *)pQualBuff,szPE1DescrBuff,pszPE1File);
							if(*pQualBuff < 64)
								*pQualBuff = 64;
							else
								*pQualBuff = 125;
							}
						Qphred = *pQualBuff - 59;	// Solexa/Illumina encodes into ascii starting from decimal 59
						Qphred = (uint8_t)(10 * log(1 + pow(10.0,((double)Qphred/10.0) / log(10.0))));	//
						if(Qphred > 40)				// clamp at phred equiv to 0.0001
							Qphred = 40;
						break;
					}
				*pQualBuff = (uint8_t)((((uint32_t)Qphred+2)*15)/40);
				}

			// pack the read and quality, read into the low order bits 0..3, quality into bits 4..7
			pQualBuff = pszPE1QualBuff;
			pReadBuff = pszPE1ReadBuff;
			for(Idx = 0; Idx < PE1ReadLen; Idx++,pQualBuff++,pReadBuff++)
				pszPE1ReadBuff[Idx] |= *pQualBuff << 4;

			if(bIsPairReads)
				{
				PE2QualLen = PE2Fasta.ReadQValues((char *)pszPE2QualBuff, cMaxReadLen-1);
				if(PE2QualLen != PE2ReadLen)		// must be same...
					{
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"Load: quality length (%d) not same as read length (%d) for '%s' entry in file '%s'",PE2QualLen,PE2ReadLen,szPE2DescrBuff,pszPE2File);
					PE1Fasta.Close();
					PE2Fasta.Close();
					if (pszPE1ReadBuff != nullptr)
						delete[]pszPE1ReadBuff;
					if (pszPE2ReadBuff != nullptr)
						delete[]pszPE2ReadBuff;
					if (pszPE1QualBuff != nullptr)
						delete[]pszPE1QualBuff;
					if (pszPE2QualBuff != nullptr)
						delete[]pszPE2QualBuff;
					return(eBSFerrParse);
					}
				// normalise the quality score to be in range 0..15 (needs to fit into 4 bits!)
				pQualBuff = pszPE2QualBuff;
				for(Idx = 0; Idx < PE2ReadLen; Idx++,pQualBuff++)
					{
					switch(m_QMethod) {
						case eFQIgnore:		// simply treat as the minimum phred
							Qphred = 0;
							break;

						case eFQSanger:		// Qphred = -10 log10(P), where Qphred is in range 0..93
							if(*pQualBuff < 33 || *pQualBuff >= 126)
								{
								if(!PE2NumInvalValues++)
									gDiagnostics.DiagOut(eDLWarn,gszProcName,"Load: quality value of '%c' outside range expected in read '%s' from file '%s' for Sanger",*(char *)pQualBuff,szPE2DescrBuff,pszPE2File);
								if(*pQualBuff < 33)
									*pQualBuff = 33;
								else
									*pQualBuff = 125;
								}
							Qphred = *pQualBuff - 33;	// Sanger encodes into ascii starting from decimal 33 '!'
							if(Qphred > 40)				// clamp at phred equiv to 0.0001
								Qphred = 40;
							break;

						case eFQIllumia:	//Qphred = -10 log10(P), where Qphred is in range 0..63
							if(*pQualBuff < 64 || *pQualBuff >= 126)
								{
								if(!PE2NumInvalValues++)
									gDiagnostics.DiagOut(eDLWarn,gszProcName,"Load: quality value of '%c' outside range expected in read '%s' from file '%s' for Illumina 1.3+",*(char *)pQualBuff,szPE2DescrBuff,pszPE2File);
								if(*pQualBuff < 64)
									*pQualBuff = 64;
								else
									*pQualBuff = 125;
								}

							Qphred = *pQualBuff - 64;	//Illumia encodes into ascii starting from decimal 64
							if(Qphred > 40)				// clamp at phred equiv to 0.0001
								Qphred = 40;
							break;

						case eFQSolexa:		// SolexaQ = -10 log10(P/(1-P)), where SolexaQ is in range -5 to 62, note the negative value
											// negative values will result if P > 0.5
											// $Q = 10 * log(1 + 10 ** (ord(SolexaQphred) - 64) / 10.0)) / log(10);
											// once Qphred is over about 15 then essentially same as Sanger and Illumina 1.3+ so
											// is it worth doing the full conversion????
							if(*pQualBuff < 59 || *pQualBuff >= 126)
								{
								if(!PE2NumInvalValues++)
									gDiagnostics.DiagOut(eDLWarn,gszProcName,"Load: quality value of '%c' outside range expected in read '%s' from file '%s' for Solexa/Illumina pre 1.3",*(char *)pQualBuff,szPE2DescrBuff,pszPE2File);
								if(*pQualBuff < 64)
									*pQualBuff = 64;
								else
									*pQualBuff = 125;
								}
							Qphred = *pQualBuff - 59;	// Solexa/Illumina encodes into ascii starting from decimal 59
							Qphred = (uint8_t)(10 * log(1 + pow(10.0,((double)Qphred/10.0) / log(10.0))));	//
							if(Qphred > 40)				// clamp at phred equiv to 0.0001
								Qphred = 40;
							break;
						}
					*pQualBuff = (uint8_t)((((uint32_t)Qphred+2)*15)/40);
					}

					// pack the read and quality, read into the low order bits 0..3, quality into bits 4..7
				pQualBuff = pszPE2QualBuff;
				pReadBuff = pszPE2ReadBuff;
				for(Idx = 0; Idx < PE2ReadLen; Idx++,pQualBuff++,pReadBuff++)
					pszPE2ReadBuff[Idx] |= *pQualBuff << 4;
				}
			}

		// apply any end trims, note that because quality scores packed into same byte as the base then end trims also trim the quality scores...
		if(m_Trim5 || ContamLen5PE1)
			{
			PE1ReadLen -= (m_Trim5 + ContamLen5PE1);
			memmove(pszPE1ReadBuff,&pszPE1ReadBuff[m_Trim5 + ContamLen5PE1],PE1ReadLen);
			}
		if(m_Trim3 || ContamLen3PE1)
			PE1ReadLen -= (m_Trim3 + ContamLen3PE1);
		if(ContamLen5PE1 > 0)
			NumContamLen5PE1 += 1;
		if(ContamLen3PE1 > 0)
			NumContamLen3PE1 += 1;


		// truncate descriptors at 1st whitespace
		for(Idx = 0; Idx < cMaxDescrIDLen-1; Idx++)
			{
			if(szPE1DescrBuff[Idx] == '\0' || isspace(szPE1DescrBuff[Idx]))
				break;
			}
		szPE1DescrBuff[Idx] = '\0';
		PE1DescrLen = Idx;

		if(bIsPairReads)
			{
			// apply any end trims
			if(m_Trim5 || ContamLen5PE2)
				{
				PE2ReadLen -= (m_Trim5 + ContamLen5PE2);
				memmove(pszPE2ReadBuff,&pszPE2ReadBuff[m_Trim5 + ContamLen5PE2],PE2ReadLen);
				}
			if(m_Trim3 || ContamLen3PE2)
				PE2ReadLen -= (m_Trim3 + ContamLen3PE2);
			if(ContamLen5PE2 > 0)
				NumContamLen5PE2 += 1;
			if(ContamLen3PE2 > 0)
				NumContamLen3PE2 += 1;

			// truncate descriptors at 1st whitespace unless fasta was generated by simulating reads in which case
			// the whole descriptor is retained as where the read was simulated from is of interest
			for(Idx = 0; Idx < cMaxDescrIDLen-1; Idx++)
				{
				if(szPE2DescrBuff[Idx] == '\0' || isspace(szPE2DescrBuff[Idx]))
					break;
				}
			szPE2DescrBuff[Idx] = '\0';
			PE2DescrLen = Idx;
			}

		if(m_TermBackgoundThreads != 0)	// need to immediately self-terminate?
			break;

		if(!bIsPairReads)
			{
			if((Rslt=(teBSFrsltCodes)AddEntry(bIsPairReads,false,0,FileID,PE1DescrLen,(char *)szPE1DescrBuff,PE1ReadLen,pszPE1ReadBuff))!=eBSFSuccess)
				break;
			PE1NumReadsAccepted += 1;
			}
		else
			{
			if((Rslt=(teBSFrsltCodes)AddEntry(bIsPairReads,false,PairReadID,FileID,PE1DescrLen,(char *)szPE1DescrBuff,PE1ReadLen,pszPE1ReadBuff))!=eBSFSuccess)
				break;
			PE1NumReadsAccepted += 1;
			if((Rslt=(teBSFrsltCodes)AddEntry(bIsPairReads,true,PairReadID,FileID+1,PE2DescrLen,(char *)szPE2DescrBuff,PE2ReadLen,pszPE2ReadBuff))!=eBSFSuccess)
				break;
			PE2NumReadsAccepted += 1;
			PairReadID += 1;
			}
		continue;
		}
	else
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"raw sequence file '%s' processing error: %s ",pszPE1File,
											Rslt > 0 ? "over length read" : "not a multifasta short reads or fastq file");
		PE1Fasta.Close();
		if(bIsPairReads)
			PE2Fasta.Close();
		if (pszPE1ReadBuff != nullptr)
			delete[]pszPE1ReadBuff;
		if (pszPE2ReadBuff != nullptr)
			delete[]pszPE2ReadBuff;
		if (pszPE1QualBuff != nullptr)
			delete[]pszPE1QualBuff;
		if (pszPE2QualBuff != nullptr)
			delete[]pszPE2QualBuff;
		return(eBSFerrParse);
		}
	}

if(Rslt != eBSFSuccess)
	{
	if(m_TermBackgoundThreads == 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Errors processing file: %s ",pszPE1File);
		while(PE1Fasta.NumErrMsgs())
			gDiagnostics.DiagOut(eDLFatal,gszProcName,PE1Fasta.GetErrMsg());
		}
	PE1Fasta.Close();
	if(bIsPairReads)
		PE2Fasta.Close();
	if (pszPE1ReadBuff != nullptr)
		delete[]pszPE1ReadBuff;
	if (pszPE2ReadBuff != nullptr)
		delete[]pszPE2ReadBuff;
	if (pszPE1QualBuff != nullptr)
		delete[]pszPE1QualBuff;
	if (pszPE2QualBuff != nullptr)
		delete[]pszPE2QualBuff;
	return(Rslt);
	}
PE1Fasta.Close();
if(bIsPairReads)
	PE2Fasta.Close();
if (pszPE1ReadBuff != nullptr)
delete[]pszPE1ReadBuff;
if (pszPE2ReadBuff != nullptr)
delete[]pszPE2ReadBuff;
if (pszPE1QualBuff != nullptr)
delete[]pszPE1QualBuff;
if (pszPE2QualBuff != nullptr)
delete[]pszPE2QualBuff;
if(m_TermBackgoundThreads != 0)	// need to immediately self-terminate?
	return(eBSErrSession);

AcquireSerialise();
m_FinalReadID = m_NumDescrReads;
m_NumReadsLoaded = m_NumDescrReads;
ReleaseSerialise();

gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadReads: Total of %1.9d reads parsed and loaded from %s",PE1NumReadsAccepted,pszPE1File);
if(PE1NumInvalValues > 0)
	gDiagnostics.DiagOut(eDLWarn,gszProcName,"Load: total of %d unexpected quality values read from file '%s'",PE1NumInvalValues,pszPE1File);
if(PE1NumUnsupportedBases > 0)
	gDiagnostics.DiagOut(eDLWarn,gszProcName,"Load: total of %d unsupported bases read from file '%s'",PE1NumUnsupportedBases,pszPE1File);
if(PE1NumUnderlength > 0)
	gDiagnostics.DiagOut(eDLWarn,gszProcName,"Load: total of %d under length sequences sloughed from file '%s'",PE1NumUnderlength,pszPE1File);
if(PE1NumOverlength > 0)
	gDiagnostics.DiagOut(eDLWarn,gszProcName,"Load: total of %d over length sequences sloughed from file '%s'",PE1NumOverlength,pszPE1File);
if(m_pContaminants != nullptr)
	{
	gDiagnostics.DiagOut(eDLWarn,gszProcName,"Load: total of %d sequences PE1 sequences were 5' contaminate trimmed",NumContamLen5PE1);
	gDiagnostics.DiagOut(eDLWarn,gszProcName,"Load: total of %d sequences PE1 sequences were 3' contaminate trimmed",NumContamLen3PE1);
	}

if(bIsPairReads)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadReads: Total of %1.9d reads parsed and loaded from %s",PE2NumReadsAccepted,pszPE2File);
	if(PE2NumInvalValues > 0)
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"Load: total of %d unexpected quality values read from file '%s'",PE2NumInvalValues,pszPE2File);
	if(PE2NumUnsupportedBases > 0)
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"Load: total of %d unsupported bases read from file '%s'",PE2NumUnsupportedBases,pszPE2File);
	if(PE2NumUnderlength > 0)
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"Load: total of %d under length sequences sloughed from file '%s'",PE2NumUnderlength,pszPE2File);
	if(PE2NumOverlength > 0)
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"Load: total of %d over length sequences sloughed from file '%s'",PE2NumOverlength,pszPE2File);
	if(m_pContaminants != nullptr)
		{
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"Load: total of %d sequences PE1 sequences were 5' contaminant trimmed",NumContamLen5PE2);
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"Load: total of %d sequences PE1 sequences were 3' contaminant trimmed",NumContamLen3PE2);
		}
	}
return(eBSFSuccess);
}





