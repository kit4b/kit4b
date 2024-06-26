// usimdiffexpr.cpp : Defines the entry point for the console application
// simulates differentially expressed transcripts
// transcript lengths are from BED transcriptomes

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

const int cDfltNumReads = 20;			// default number of million reads 
const int cMinNumReads =  1;			// minimum number of million reads
const int cMaxNumReads =  100;			// maximum number of million reads

const int cDfltVaryReads = 10;			// default reads per replicate variation percentage
const int cMaxVaryReads = 30;			// default reads per replicate variation percentage

const int cDfltWhiteNoise = 100;		// default white noise to apply to read acceptance thresholds each time read sampled for count
const int cMaxWhiteNoise = 900;			// max white noise to apply to read acceptance thresholds each time read sampled for count

const int cMaxNumReps  = 100;			// maximum number of replicates supported

const double cDfltShotNoise = 0.0;			// mean shot noise per transcript as proportion of counts
const double cMaxShotNoise = 0.5;			// max mean shot noise per transcript as proportion of counts

// processing modes
typedef enum TAG_ePMode {
	ePMStandard,					// default - standard uniform read sampling probabilities
	ePMrandom,						// random read sampling probabilities
	ePMrandprof,					// non-linear profiled random read sampling probabilities
	ePMplaceholder					// used to set the enumeration range
	} etPMode;

// output format
typedef enum TAG_eFMode {
	eFMcsv,					// default - CSV 
	eFMtab,					// tab delimited
	eFMplaceholder			// used to set the enumeration range
	} etFMode;


int
Process(etPMode PMode,		// processing mode
		etFMode FMode,		// output format
		int WhiteNoise,		// white noise, randomly select up to this value to add/remove to read threshold every time the read is tested for accepting as a count
		int NumReads,		// number of reads required
		int VaryReads,		// vary number of reads by up to this percentage of NumReads
		int VaryTrans,		// after half replicates then vary this percentage of transcripts have differential expression
		int NumReps,		// number of replicates
		int NumThreads,		// number of worker threads to use
		char *pszInFile,	// generated counts are for these features or genes in this BED file
		char *pszOutFile);	// output simulated transcript counts to this file

int RunSamplingThreads(int NumThreads, int Num2Sample, int RepIdxStart, int NumReps);

static int SortSimReadEls(const void *arg1, const void *arg2);

CStopWatch gStopWatch;
CDiagnostics gDiagnostics;				// for writing diagnostics messages to log file
char gszProcName[_MAX_FNAME];			// process name

#pragma pack(1)

#ifdef _WIN32
SRWLOCK m_hCntsRwLock;
unsigned __stdcall SamplingThread(void * pThreadPars);
#else
pthread_rwlock_t m_hCntsRwLock;
void *SamplingThread(void * pThreadPars);
#endif

typedef struct TAG_tsReadEl {
	int RandCnt;		// dual purpose, used when randomly sorting read order, and to hold number of times this read was accepted as a count
	int BaseThres;		// base acceptance threshold for this read
	int SampThres;		// BaseThres +/- noise component representing sequencing bias noise          
	int TransID;		// read count is will be attributed to this transcript
} tsReadEl;

typedef struct TAG_tsTransCnt {
	int TransID;		// identifies this transcript
	int Score;			// score for this transcript
	int Thres;			// reads on this transcript will be sampled if above this threshold (0..1000000)
	int TranscribedLen; // length of this transcript
	int NumReads;		// number of reads associated to this transcript
	int Counts[cMaxNumReps];	// total counts for each replicate of this transcript
} tsTransCnt;


typedef struct TAG_sSamplingThreadPars {
#ifdef _WIN32
	HANDLE threadHandle;	// handle as returned by _beginthreadex()
	unsigned int threadID;	// identifier as set by _beginthreadex()
#else
	int threadRslt;			// result as returned by pthread_create ()
	pthread_t threadID;		// identifier as set by pthread_create ()
#endif
	int ThreadIdx;			// index of this thread 0..n
	int Rslt;				// processing result, 1 for success
	etPMode PMode;			// processing mode
	int RepIdxStart;		// which replicate to start from (0..N)
	int NumReps;			// number of replicates to process
	int NumReqSamples;		// number of samples to be generated by this worker thread
	int WhiteNoise;		    // white noise, randomly select up to this value to add/remove to read threshold every time the read is tested for accepting as a count
	int NumReadEls;
	tsReadEl *pReadEls;
	tsTransCnt *pTransCnts;
	uint64_t NumTestedSamples;// number of samples tested
	int NumAcceptedSamples;	// number of actually accepted samples

} tsSamplingThreadPars;
#pragma pack()

#ifdef _WIN32
int _tmain(int argc, char* argv[])
{
// determine my process name
_splitpath(argv[0],NULL,NULL,gszProcName,NULL);
#else
int 
main(int argc, const char** argv)
{
// determine my process name
CUtility::splitpath((char *)argv[0],NULL,gszProcName);
#endif
int iScreenLogLevel;		// level of screen diagnostics
int iFileLogLevel;			// level of file diagnostics
char szLogFile[_MAX_PATH];	// write diagnostics to this file
int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure

int PMode;					// processing mode
int FMode;					// output format - csv or tab delimited
int NumReads;				// number of reads required
int VaryReads;				// vary total reads per replicate by up to this percentage
int VaryTrans;				// after half replicates then vary this proportion of transcripts
int NumReps;				// number of replicates
int NumberOfProcessors;		// number of installed CPUs
int NumThreads;				// number of threads (0 defaults to number of CPUs)

int WhiteNoise;				// induce this added shot noise as a proportion of transcript counts

char szInFile[_MAX_PATH]; // input BED feature or gene file containing transcript
char szOutFile[_MAX_PATH];	// output simulated transcript counts to this file


// command line args
struct arg_lit  *help    = arg_lit0("h","help",                "print this help and exit");
struct arg_lit  *version = arg_lit0("v","version,ver",			"print version information and exit");
struct arg_int *FileLogLevel=arg_int0("f", "FileLogLevel",		"<int>","Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
struct arg_file *LogFile = arg_file0("F","log","<file>",		"diagnostics log file");

struct arg_int *pmode = arg_int0("m","mode","<int>",		    "processing mode: 0 - uniform, 1 - linear random, 2 - non-linear profiled read sampling probabilities (0 - default)");
struct arg_int *fmode = arg_int0("M","format","<int>",		    "output format: 0 - csv, 1 - tab delimited (default: 0)");
struct arg_int *numreads = arg_int0("n","ncounts","<int>",	    "total number of counts or simulated aligned reads in millions, max of 100M (default = 50)");
struct arg_int *varyreads = arg_int0("R","rcounts","<int>",	    "randomly vary total counts between replicates by up to this percentage (0..30) of nominal total (default = 10)");

struct arg_int *varytrans = arg_int0("e","trans","<int>",	    "percentage (0..99) of transcripts with differential expression (default 0)");


struct arg_int *numreps = arg_int0("r","nreplicates","<int>",   "total number of replicates, max of 100 (default = 2)");
struct arg_int *whitenoise = arg_int0("w","whitenoise","<int>",	"white noise to apply to read acceptance thresholds between replicates representing change in sampling biases (default = 100)");
struct arg_file *infile = arg_file1("i","infile","<file>",		"use features or genes in this BED file for defining transcriptome isoform lengths");
struct arg_file *outfile = arg_file1("o","outfile","<file>",	"output transcript counts to this file");
struct arg_int *threads = arg_int0("T","threads","<int>",		"number of processing threads 0..128 (defaults to 0 which sets threads to number of CPU cores)");


struct arg_end *end = arg_end(20);

void *argtable[] = {help,version,FileLogLevel,LogFile,
					pmode,fmode,numreads,varyreads,varytrans,numreps, whitenoise,
					infile,outfile,threads,
					end};

char **pAllArgs;
int argerrors;
argerrors = CUtility::arg_parsefromfile(argc,(char **)argv,&pAllArgs);
if(argerrors >= 0)
	argerrors = arg_parse(argerrors,pAllArgs,argtable);

/* special case: '--help' takes precedence over error reporting */
if (help->count > 0)
        {
		printf("\n%s the Kanga simulated differential transcript generator, Version %s\nOptions ---\n", gszProcName,kit4bversion);
        arg_print_syntax(stdout,argtable,"\n");
        arg_print_glossary(stdout,argtable,"  %-25s %s\n");
		printf("\nNote: Parameters can be entered into a parameter file, one parameter per line.");
		printf("\n      To invoke this parameter file then precede its name with '@'");
		printf("\n      e.g. %s @myparams.txt\n",gszProcName);
		printf("\nPlease report any issues regarding usage of %s to https://github.com/kit4b/kit4b/issues\n\n",gszProcName);
		exit(1);
        }

    /* special case: '--version' takes precedence error reporting */
if (version->count > 0)
        {
		printf("\n%s Version %s\n",gszProcName,kit4bversion);
		exit(1);
        }

if (!argerrors)
	{
	if(FileLogLevel->count && !LogFile->count)
		{
		printf("\nError: FileLogLevel '-f%d' specified but no logfile '-F<logfile>'",FileLogLevel->ival[0]);
		exit(1);
		}

	iScreenLogLevel = iFileLogLevel = FileLogLevel->count ? FileLogLevel->ival[0] : eDLInfo;
	if(iFileLogLevel < eDLNone || iFileLogLevel > eDLDebug)
		{
		printf("\nError: FileLogLevel '-l%d' specified outside of range %d..%d",iFileLogLevel,eDLNone,eDLDebug);
		exit(1);
		}
	if(LogFile->count)
		{
		strncpy(szLogFile,LogFile->filename[0],_MAX_PATH);
		szLogFile[_MAX_PATH-1] = '\0';
		}
	else
		{
		iFileLogLevel = eDLNone;
		szLogFile[0] = '\0';
		}

		// now that log parameters have been parsed then initialise diagnostics log system
	if(!gDiagnostics.Open(szLogFile,(etDiagLevel)iScreenLogLevel,(etDiagLevel)iFileLogLevel,true))
		{
		printf("\nError: Unable to start diagnostics subsystem\n");
		if(szLogFile[0] != '\0')
			printf(" Most likely cause is that logfile '%s' can't be opened/created\n",szLogFile);
		exit(1);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Version: %s",kit4bversion);


	PMode = (etPMode)(pmode->count ? pmode->ival[0] : ePMStandard);
	if(PMode < ePMStandard || PMode >= ePMplaceholder)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Processing mode '-m%d' specified outside of range %d..%d",PMode,(int)ePMStandard,(int)ePMplaceholder-1);
		exit(1);
		}


	FMode = (etFMode)(fmode->count ? fmode->ival[0] : eFMcsv);
	if(FMode < eFMcsv || FMode >= eFMplaceholder)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Output format '-m%d' specified outside of range %d..%d",PMode,(int)eFMcsv,(int)eFMplaceholder-1);
		exit(1);
		}

#ifdef _WIN32
	SYSTEM_INFO SystemInfo;
	GetSystemInfo(&SystemInfo);
	NumberOfProcessors = SystemInfo.dwNumberOfProcessors;
#else
	NumberOfProcessors = sysconf(_SC_NPROCESSORS_CONF);
#endif

	int MaxAllowedThreads = min(cMaxWorkerThreads,NumberOfProcessors);	// limit to be at most cMaxWorkerThreads
	if((NumThreads = threads->count ? threads->ival[0] : MaxAllowedThreads)==0)
		NumThreads = MaxAllowedThreads;
	if(NumThreads < 0 || NumThreads > MaxAllowedThreads)
		{
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"Warning: Number of threads '-T%d' specified was outside of range %d..%d",NumThreads,1,MaxAllowedThreads);
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"Warning: Defaulting number of threads to %d",MaxAllowedThreads);
		NumThreads = MaxAllowedThreads;
		}

	NumReads = numreads->count ? numreads->ival[0] : cDfltNumReads;
	if(NumReads < cMinNumReads || NumReads > cMaxNumReads)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Number of reads '-n%d' (millions) specified outside of range %d..%d",NumReads,cMinNumReads,cMaxNumReads);
		exit(1);
		}

	VaryReads = varyreads->count ? varyreads->ival[0] : cDfltVaryReads;
	if(VaryReads < 0 || VaryReads > cMaxVaryReads)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Variation in total reads between replicates '-R%d' specified outside of range 0..%d",VaryReads,cMaxVaryReads);
		exit(1);
		}

	VaryTrans = varytrans->count ? varytrans->ival[0] : 0;
	if(VaryTrans < 0 || VaryTrans > 99)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Variation in transcript differential expression '-e%d' specified outside of range 0..99",VaryTrans);
		exit(1);
		}

	NumReps = numreps->count ? numreps->ival[0] : 2;
	if(NumReps < 2 || NumReps > cMaxNumReps)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Number of replicates '-r%d' specified outside of range 2..%d",NumReps,cMaxNumReps);
		exit(1);
		}

	WhiteNoise = whitenoise->count ? whitenoise->ival[0] : cDfltWhiteNoise;
	if(WhiteNoise < 0 || WhiteNoise >= cMaxWhiteNoise)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: White noise '-w%d' specified outside of range 0 to %d",WhiteNoise,cMaxWhiteNoise);
		exit(1);
		}


	strcpy(szInFile,infile->filename[0]);
	strcpy(szOutFile,outfile->filename[0]);

// show user current resource limits
#ifndef _WIN32
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Resources: %s",CUtility::ReportResourceLimits());
#endif

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing parameters:");

	const char *pszDescr;
	switch(PMode) {
		case ePMStandard:				
			pszDescr = "uniform";
			break;
		case ePMrandom:
			pszDescr = "random";
			break;
		case ePMrandprof:
			pszDescr = "non-linear profiled random read sampling probabilities";
			break;
		}
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"processing mode is : %s read sampling probabilities",pszDescr);

	switch(FMode) {
		case eFMcsv:				
			pszDescr = "CSV";
			break;
		case eFMtab:				
			pszDescr = "tab delimited";
			break;
		}


	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output format is : '%s'",pszDescr);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Total number of read counts (millions) required: %d", NumReads);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Vary total number of reads per replicate by up to this percentage: %d", VaryReads);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Percentage of transcripts with differential expression: %d", VaryTrans);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Number of replicates: %d", NumReps);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"white noise to apply to read acceptance thresholds between replicates representing change in sampling biases: %d",WhiteNoise);

	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Input transcriptome BED feature or gene file: '%s'",szInFile);
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Output file: '%s'",szOutFile);

    gDiagnostics.DiagOutMsgOnly(eDLInfo,"number of threads : %d",NumThreads);
#ifdef _WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	gStopWatch.Start();
	Rslt = Process((etPMode)PMode,(etFMode)FMode,WhiteNoise,NumReads,VaryReads,VaryTrans,NumReps,NumThreads,szInFile,szOutFile);
	gStopWatch.Stop();
	Rslt = Rslt >=0 ? 0 : 1;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Exit code: %d Total processing time: %s",Rslt,gStopWatch.Read());
	exit(Rslt);
	}
else
	{
	printf("\n%s Kanga simulated differential transcript generator, Version %s\n",gszProcName,kit4bversion);
	arg_print_errors(stdout,end,gszProcName);
	arg_print_syntax(stdout,argtable,"\nUse '-h' to view option and parameter usage\n");
	exit(1);
	}
return 0;
}


etPMode m_PMode;				// processing mode
etFMode m_FMode;				// requested output format
int m_WhiteNoise;					// white noise, randomly select up to this value to add/remove to read threshold every time the read is tested for accepting as a count

int m_hOutFile;					// output results file handle
int m_NumReadEls;				// number of reads for which counts are to be generated
tsReadEl *m_pReadEls;			// allocated to hold all reads

int m_NumFeatures;				// number of features in BED file and number of transcripts
tsTransCnt *m_pTransCnts;		// used to hold transcript counts
CBEDfile *m_pBEDFile;			// holds gene/features representing transcripts

uint64_t m_NumTestedSamples;		// number of sample reads tested
uint32_t m_NumAcceptedSamples;    // number of sample reads accepted as counts 

int m_TotTranscribedLen;		// total transcribed length of all gene/features


void
Reset(void)
{
if(m_hOutFile != -1)
	{
	close(m_hOutFile);
	m_hOutFile = -1;
	}
if(m_pReadEls != NULL)
	{
	delete m_pReadEls;
	m_pReadEls = NULL;
	}
if(m_pTransCnts != NULL)
	{
	delete m_pTransCnts;
	m_pTransCnts = NULL;
	}
if(m_pBEDFile != NULL)
	{
	delete m_pBEDFile;
	m_pBEDFile = NULL;
	}
m_NumReadEls = 0;
m_NumFeatures = 0;
m_TotTranscribedLen = 0;
m_NumTestedSamples = 0;
m_NumAcceptedSamples = 0;
}

void
Init(void)
{
m_pReadEls = NULL;
m_pTransCnts = NULL;
m_pBEDFile = NULL;
m_hOutFile = -1;
Reset();
}



int
Process(etPMode PMode,		// processing mode
		etFMode FMode,		// output format
		int WhiteNoise,		// white noise, randomly select up to this value to add/remove to read threshold every time the read is tested for accepting as a count
		int NumReads,		// number of reads required
		int VaryReads,		// vary number of reads by up to this percentage of NumReads
		int VaryTrans,		// after half replicates then vary this percentage of transcripts have differential expression
		int NumReps,		// number of replicates
		int NumThreads,		// number of worker threads to use
		char *pszInFile,	// generated counts are for these features or genes in this BED file
		char *pszOutFile)	// output simulated transcript counts to this file
{
int Rslt;
int Num2Vary;
int Num2Sample;
char szBuff[32000];
int BuffIdx;
int *pCnt;
char szFeatName[128];
tsReadEl *pReadEl;
tsTransCnt *pTransCnt;
int ElIdx;
int NumReadsCurTrans;
int CurFeatLen;
int TotRemainFeatLen;
int RepIdx;
int CurFeatureID;

int NumTransVaried;
int NumTrans2Vary;
bool bUpReg;
int	NumUpRegulated;
int	NumDnRegulated;

Init();

m_PMode = PMode;
m_FMode = FMode;
m_WhiteNoise = WhiteNoise;

// seed random number generator
uint32_t RandSeed;
#ifdef _WIN32
RandSeed = (uint32_t)(uint64_t)GetTickCount();   // number of millisecs since system was started
#else
struct tms Times;
RandSeed = (uint32_t)times(&Times);			// number of clock ticks since an arbitrary point in the past - it's only a seed so not too concerned as to when this point was
#endif
if(!RandSeed)								// just in case
	RandSeed = (uint32_t)time(NULL);
TRandomCombined<CRandomMother,CRandomMersenne> RGseeds(RandSeed);
// generator can take a number of iterations to settle down before generating reasonable psuedo-random values so run it for a million cycles
for(RandSeed = 0; RandSeed < 100000; RandSeed++)
	RGseeds.IRandom(0,1000000);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Opening BED file '%s'",pszInFile);
// open feature file
if((m_pBEDFile = new CBEDfile) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBEDfile object");
	return(eBSFerrObj);
	}
if((Rslt = m_pBEDFile->Open(pszInFile))!=eBSFSuccess)
	{
	while(m_pBEDFile->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pBEDFile->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open feature/gene BED file '%s'",pszInFile);
	Reset();
	return(Rslt);
	}
if((m_NumFeatures = m_pBEDFile->GetNumFeatures()) == 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Feature/gene BED file '%s' is empty",pszInFile);
	Reset();
	return(Rslt);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"BED file '%s' contains %d features or genes",pszInFile,m_NumFeatures);


gDiagnostics.DiagOut(eDLInfo,gszProcName,"Creating/truncating output counts file '%s'",pszOutFile);
#ifdef _WIN32
m_hOutFile = open(pszOutFile,O_CREATETRUNC );
#else
if((m_hOutFile = open(pszOutFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
    if(ftruncate(m_hOutFile,0)!=0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",pszOutFile,strerror(errno));
			Reset();
			return(eBSFerrCreateFile);
			}
#endif

if(m_hOutFile < 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output file '%s'",pszOutFile);
	Reset();
	return(eBSFerrCreateFile);
	}

if(FMode == eFMcsv)
	{
	BuffIdx = sprintf(szBuff,"\"TransID\",\"Transcript\",\"TransFoldChange\",\"CtrlTransThres\",\"ExprTransThres\",\"TransLen\",\"TransReads\",\"CtrlMean\",\"ExprMean\",\"CtrlExprFoldChange\",\"MaxRepsFoldChange\"");
	for(RepIdx = 0; RepIdx < NumReps; RepIdx++)
		if(VaryTrans > 0 && RepIdx >= NumReps/2)
			BuffIdx += sprintf(&szBuff[BuffIdx],",\"Expr:%d\"",1 + RepIdx - (NumReps/2));
		else
			BuffIdx += sprintf(&szBuff[BuffIdx],",\"Ctrl:%d\"",RepIdx+1);
	BuffIdx += sprintf(&szBuff[BuffIdx],"\n");
	}
else
	{
	BuffIdx = 0;
	for(RepIdx = 0; RepIdx < NumReps; RepIdx++)
		if(VaryTrans > 0 && RepIdx >= NumReps/2)
			BuffIdx += sprintf(&szBuff[BuffIdx],"\tExpr:%d",1 + RepIdx - (NumReps/2));
		else
			BuffIdx += sprintf(&szBuff[BuffIdx],"\tCtrl:%d",RepIdx+1);
	BuffIdx += sprintf(&szBuff[BuffIdx],"\n");
	}
if(write(m_hOutFile,szBuff,BuffIdx) != BuffIdx)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Errors whilst writing %d bytes - '%s' - %s",BuffIdx, pszOutFile, strerror(errno));
	Reset();
	return(eBSFerrFileAccess);
	}
BuffIdx = 0;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Initialising ready to start sampling...");


// allocate for required number of transcripts
if((m_pTransCnts = new tsTransCnt [m_NumFeatures]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Memory allocation for %d transcripts failed",m_NumFeatures);
	Reset();
	return(eBSFerrMem);
	}

// allocate for required number of reads plus max allowance for the variation
NumReads *= 1000000;		// reads were specified as millions
if(VaryReads > 0)
	Num2Vary = (int)(((uint64_t)NumReads * VaryReads)/100);
else
	Num2Vary = 0;
if((m_pReadEls = new tsReadEl [NumReads + Num2Vary]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadFasta: Memory allocation for %d read counts failed",NumReads+Num2Vary);
	Reset();
	return(eBSFerrMem);
	}
m_NumReadEls = NumReads + Num2Vary;

// determine transcript lengths and total length
// needed so that reads can be uniformly distributed along the transcripts
m_TotTranscribedLen = 0;
CurFeatureID = 0;
pTransCnt = m_pTransCnts;
while((CurFeatureID = m_pBEDFile->GetNextFeatureID(CurFeatureID)) > 0)
	{
	memset(pTransCnt,0,sizeof(*pTransCnt));
	pTransCnt->TransID = CurFeatureID;
	pTransCnt->Score = min(1000,m_pBEDFile->GetFeatScore(CurFeatureID));
	if(pTransCnt->Score == 0)						// if no score in BED file then randomly associate a score in the range 1..1000
		pTransCnt->Score = RGseeds.IRandom(1,1000);
	pTransCnt->Thres = pTransCnt->Score;		
	pTransCnt->TranscribedLen = max(1,m_pBEDFile->GetTranscribedLen(CurFeatureID));
	m_TotTranscribedLen += pTransCnt->TranscribedLen;
	pTransCnt += 1;
	}

// can now start allocating reads to features
// allocation is in proportion to either feature length or randomly
// if proportional to feature length then feature X gets featlen/totfeatlen * totreads
NumReadsCurTrans = 0;
TotRemainFeatLen = m_TotTranscribedLen;
pReadEl = m_pReadEls;
pTransCnt = m_pTransCnts;
for(ElIdx = 0; ElIdx < m_NumReadEls; ElIdx++, pReadEl++)
	{
	if(NumReadsCurTrans == 0)		// 0 if starting to allocate for a new transcript
		{
		CurFeatLen = pTransCnt->TranscribedLen;
		CurFeatureID = pTransCnt->TransID;
    	NumReadsCurTrans = max(1,(int)(((uint64_t)(m_NumReadEls - ElIdx) * CurFeatLen)/TotRemainFeatLen));
		pTransCnt->NumReads = NumReadsCurTrans;
		TotRemainFeatLen -= CurFeatLen;
		pTransCnt += 1;
		}
	switch(PMode) {
		case ePMStandard:			// default - standard uniform read sampling probabilities
			pReadEl->BaseThres = 500;
			break;
		case ePMrandom:				// random read sampling probabilities
			pReadEl->BaseThres = RGseeds.IRandom(1,1000);
			break;
		case ePMrandprof:					// non-linear profiled random read sampling probabilities
			pReadEl->BaseThres = max(1,(int)((1000.0 / (double)RGseeds.IRandom(1,1000))));
			break;
		}
	pReadEl->SampThres = pReadEl->BaseThres;
	pReadEl->TransID = CurFeatureID;
	NumReadsCurTrans -= 1;
	}

// need to serialise access to counts updating between threads, do this with an exclusive read/write lock
#ifdef _WIN32
InitializeSRWLock(&m_hCntsRwLock);
#else
if(pthread_rwlock_init (&m_hCntsRwLock,NULL)!=0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to create rwlock");
	return(eBSFerrInternal);
	}
#endif

// at long last, for each replicate, start sampling the reads with replacement
for(RepIdx = 0; RepIdx < NumReps;RepIdx++)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Initialising for replicate %d, randomising read element order...",RepIdx+1);

	// now sort so that read positions are randomly distributed
	// object is to reduce the impact of any random generator biases when later randomly selecting a read 
	// and comparing a randomly generated value against that read's threshold
	pReadEl = m_pReadEls;
	for(ElIdx = 0; ElIdx < m_NumReadEls; ElIdx++, pReadEl++)
		pReadEl->RandCnt = RGseeds.IRandom(1,m_NumReadEls * 10);	// upper limit much larger than m_NumReadEls so as to reduce chances of duplicate RandCnt values  
	qsort(m_pReadEls,m_NumReadEls,sizeof(tsReadEl),SortSimReadEls);
		pReadEl = m_pReadEls;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Randomisation of read element order completed");
	
	// reset counts for each read and also apply any read bias noise
	for(ElIdx = 0; ElIdx < m_NumReadEls; ElIdx++, pReadEl++)
		{
		pReadEl->RandCnt = 0;										// RandCnt also used to hold number of cnts resulting from sampling of this read
		if(WhiteNoise > 0)
			pReadEl->SampThres = RGseeds.IRandom(max(1,pReadEl->BaseThres-WhiteNoise),min(1000,pReadEl->BaseThres+WhiteNoise));
		}

	// if user doing control and experiment then change the transcript abundance through up/down on the individual transcript thresholds
	if(VaryTrans > 0 && RepIdx == NumReps/2)
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Starting experimental replicates at replicate %d, Up/Down regulating %d%% transcripts...",RepIdx+1,VaryTrans);
		NumTrans2Vary = (99 + (m_NumFeatures * VaryTrans)) / 100;
		NumTransVaried = 0;
		NumUpRegulated = 0;
		NumDnRegulated = 0;
		while(NumTransVaried < NumTrans2Vary) { 
			ElIdx = RGseeds.IRandom(0,m_NumFeatures-1);
			pTransCnt = &m_pTransCnts[ElIdx];
			if(pTransCnt->Thres != pTransCnt->Score)	// if already varied then don't re-vary
				continue;
			
			// simulate at least a 10% change, or 10 which ever is larger, in threshold
			int MinDelta = max(10,pTransCnt->Score / 10);

			if(pTransCnt->Score < (MinDelta + 10))
				bUpReg = true;
			else
				if((pTransCnt->Score + MinDelta) > 990)
					bUpReg = false;
				else
					bUpReg = RGseeds.IRandom(1,100000) > 50000 ? false : true; // using 1..100000 instead of 0..1 just in case there is a significant bias with such a restricted range

			if(bUpReg)
				{
				pTransCnt->Thres = RGseeds.IRandom(pTransCnt->Score + MinDelta,999);
				NumUpRegulated += 1;
				}
			else
				{
				pTransCnt->Thres = RGseeds.IRandom(1, pTransCnt->Score - MinDelta);
				NumDnRegulated += 1;
				}
			NumTransVaried += 1;
			}
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Up regulated %d, down regulated %d transcripts...",NumUpRegulated,NumDnRegulated);
		}

	pReadEl = m_pReadEls;
	if(VaryReads == 0)
		Num2Sample = m_NumReadEls;
	else
		Num2Sample = RGseeds.IRandom(NumReads - Num2Vary,NumReads + Num2Vary);

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Sampling for total of %d counts on replicate %d...",Num2Sample,RepIdx+1);
	m_NumTestedSamples = 0;
	m_NumAcceptedSamples = 0;
	RunSamplingThreads(NumThreads,Num2Sample,RepIdx,1);
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Sampling on replicate %d completed",RepIdx+1);
	}

#ifndef _WIN32
pthread_rwlock_destroy(&m_hCntsRwLock);
#endif

double CtrlMean;
double ExprMean;
double CntsFoldChange;
double TransFoldChange;
int RepMin;
int RepMax;
double RepsCntsFoldChange;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Replicates all sampled now writing replicate counts to '%s'",pszOutFile);
pTransCnt = m_pTransCnts;
BuffIdx = 0;
for(CurFeatureID = 0; CurFeatureID < m_NumFeatures; CurFeatureID++, pTransCnt++)
	{
	m_pBEDFile->GetFeature(pTransCnt->TransID,szFeatName);
	switch(FMode) {
		case eFMcsv:
			// what fold change between the control and experiment replicates means
			// note that a fudge is used if the denominator is ever 0 by accepting the numerator as the fold change
			pCnt = pTransCnt->Counts;
			CtrlMean = 0;
			ExprMean = 0;
			RepMin = -1;
			RepMax = -1;
			for(RepIdx = 0; RepIdx < NumReps; RepIdx++,pCnt++)
				{
				if(RepIdx < NumReps/2)
					CtrlMean += (double)*pCnt;
				else
					ExprMean += (double)*pCnt;
				if(*pCnt > RepMax)
					RepMax = *pCnt;
				if(RepMin == -1 || *pCnt < RepMin)
					RepMin = *pCnt;
				}

			CtrlMean /= (double)(NumReps/2);
			ExprMean /= (double)(NumReps - (NumReps/2));
				// a fudge to prevent divide by  zero issues
			if(CtrlMean == 0.0)
				CntsFoldChange = ExprMean;
			else
				if(ExprMean == 0.0)
					CntsFoldChange = -1 * CtrlMean;
			if(CtrlMean > 0.0 && ExprMean > 0.0)
				{
				CntsFoldChange = ExprMean/CtrlMean;
				if(CntsFoldChange < 1.0)
					CntsFoldChange = -1.0/CntsFoldChange;
				}

			// now for the transcript fold change
			if(pTransCnt->Score == 0)
				TransFoldChange = (double)pTransCnt->Thres;
			else
				if(pTransCnt->Thres == 0)
					TransFoldChange = -1.0 * pTransCnt->Score;
			if(pTransCnt->Score > 0 && pTransCnt->Thres > 0)
				{
				TransFoldChange = (double)pTransCnt->Thres/pTransCnt->Score;
				if(TransFoldChange < 1.0)
					TransFoldChange = -1.0/TransFoldChange;
				}

			// now for the worst case replicate counts fold change
			if(RepMin == 0)
				RepsCntsFoldChange = (double)RepMax;
			else
				RepsCntsFoldChange = (double)RepMax/RepMin;
			
			BuffIdx += sprintf(&szBuff[BuffIdx],"%d,\"%s\",%1.3f,%d,%d,%d,%d,",pTransCnt->TransID,szFeatName,TransFoldChange,pTransCnt->Score,pTransCnt->Thres,pTransCnt->TranscribedLen,pTransCnt->NumReads);
			BuffIdx += sprintf(&szBuff[BuffIdx],"%1.3f,%1.3f,%1.3f,%1.3f",CtrlMean,ExprMean,CntsFoldChange,RepsCntsFoldChange);

			pCnt = pTransCnt->Counts;
			for(RepIdx = 0; RepIdx < NumReps; RepIdx++,pCnt++)
				BuffIdx += sprintf(&szBuff[BuffIdx],",%d",*pCnt);
			break;
		case eFMtab:
			BuffIdx += sprintf(&szBuff[BuffIdx],"%s",szFeatName);
			pCnt = pTransCnt->Counts;
			for(RepIdx = 0; RepIdx < NumReps; RepIdx++,pCnt++)
				BuffIdx += sprintf(&szBuff[BuffIdx],"\t%d",*pCnt);
			break;
		}

	BuffIdx += sprintf(&szBuff[BuffIdx],"\n");
	if((BuffIdx + 1024) > sizeof(szBuff))
		{
		if(write(m_hOutFile,szBuff,BuffIdx) != BuffIdx)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Errors whilst writing %d bytes - '%s' - %s",BuffIdx, pszOutFile, strerror(errno));
			Reset();
			return(eBSFerrFileAccess);
			}
		BuffIdx = 0;
		}
	}
if(BuffIdx > 0)
	{
	if(write(m_hOutFile,szBuff,BuffIdx) != BuffIdx)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Errors whilst writing %d bytes - '%s' - %s",BuffIdx, pszOutFile, strerror(errno));
		Reset();
		return(eBSFerrFileAccess);
		}
	}
close(m_hOutFile);
m_hOutFile = -1;
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Simulation completed");

Reset();
return(0);
}

void
ShowProgress(void)
{
uint64_t NumTestedSamples;
uint32_t NumAcceptedSamples;
#ifdef _WIN32
		AcquireSRWLockExclusive(&m_hCntsRwLock);
#else
		pthread_rwlock_wrlock(&m_hCntsRwLock);
#endif
		NumTestedSamples = m_NumTestedSamples;
		NumAcceptedSamples = m_NumAcceptedSamples;
#ifdef _WIN32
		ReleaseSRWLockExclusive(&m_hCntsRwLock);
#else
		pthread_rwlock_unlock(&m_hCntsRwLock);
#endif
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress, from %zd sampled reads have accepted %d as counts...",NumTestedSamples,NumAcceptedSamples);
}

int
RunSamplingThreads(int NumThreads, int Num2Sample, int RepIdxStart, int NumReps)
{
int ThreadIdx;
tsSamplingThreadPars SamplingThreads[cMaxWorkerThreads];

memset(SamplingThreads,0,sizeof(SamplingThreads));

for(ThreadIdx = 0; ThreadIdx < NumThreads; ThreadIdx++)
	{
	SamplingThreads[ThreadIdx].ThreadIdx = ThreadIdx + 1;
	SamplingThreads[ThreadIdx].NumReqSamples = Num2Sample / (NumThreads - ThreadIdx);
	Num2Sample -= SamplingThreads[ThreadIdx].NumReqSamples;
	SamplingThreads[ThreadIdx].RepIdxStart = RepIdxStart;
	SamplingThreads[ThreadIdx].NumReps = NumReps;
	SamplingThreads[ThreadIdx].NumReadEls = m_NumReadEls;
	SamplingThreads[ThreadIdx].pReadEls = m_pReadEls;
	SamplingThreads[ThreadIdx].pTransCnts = m_pTransCnts;
	SamplingThreads[ThreadIdx].PMode = m_PMode;
	SamplingThreads[ThreadIdx].WhiteNoise = m_WhiteNoise;

#ifdef _WIN32
	SamplingThreads[ThreadIdx].threadHandle = (HANDLE)_beginthreadex(NULL,0x0fffff,SamplingThread,&SamplingThreads[ThreadIdx],0,&SamplingThreads[ThreadIdx].threadID);
#else
	SamplingThreads[ThreadIdx].threadRslt =	pthread_create (&SamplingThreads[ThreadIdx].threadID , NULL , SamplingThread , &SamplingThreads[ThreadIdx] );
#endif
	}

for(ThreadIdx = 0; ThreadIdx < NumThreads; ThreadIdx++)
	{
#ifdef _WIN32
	while(WAIT_TIMEOUT == WaitForSingleObject( SamplingThreads[ThreadIdx].threadHandle, 60000))
		{
		ShowProgress();
		}
	CloseHandle( SamplingThreads[ThreadIdx].threadHandle);
#else
	struct timespec ts;
	int JoinRlt;
	clock_gettime(CLOCK_REALTIME, &ts);
	ts.tv_sec += 60;
	while((JoinRlt = pthread_timedjoin_np(SamplingThreads[ThreadIdx].threadID, NULL, &ts)) != 0)
		{
		ShowProgress();
		ts.tv_sec += 60;
		}
#endif
	}
ShowProgress();
return(eBSFSuccess);
}


#ifdef _WIN32
unsigned __stdcall SamplingThread(void * pThreadPars)
#else
void *SamplingThread(void * pThreadPars)
#endif
{
tsSamplingThreadPars *pPars = (tsSamplingThreadPars *)pThreadPars; // makes it easier not having to deal with casts!
int Num2Sample = pPars->NumReqSamples;
etPMode PMode =  pPars->PMode;
int NumReadEls = pPars->NumReadEls;
tsReadEl *pReadEls = pPars->pReadEls;
tsTransCnt *pTransCnts =  pPars->pTransCnts;
int RepStartIdx = pPars->RepIdxStart;
int NumReps = pPars->NumReps;
int	NumAcceptedSamples;
int	DeltaAcceptedSamples;
int CurThres;
int ElIdx;
int RepIdx;
uint64_t NumTestedSamples;
uint64_t DeltaTestedSamples;
tsTransCnt *pTransCnt;
tsReadEl *pReadEl;

// seed random number generator
uint32_t RandSeed;
#ifdef _WIN32
RandSeed = (uint32_t)(uint64_t)GetTickCount();   // number of millisecs since system was started
#else
struct tms Times;
RandSeed = (uint32_t)times(&Times);			// number of clock ticks since an arbitrary point in the past - it's only a seed so not too concerned as to when this point was
#endif
RandSeed += (uint32_t)pPars->threadID + (uint32_t)(uint64_t)&RandSeed;
if(!RandSeed)								// just in case
	RandSeed = (uint32_t)time(NULL);
TRandomCombined<CRandomMother,CRandomMersenne> RGseeds(RandSeed);
for(ElIdx = 0; ElIdx < 10000; ElIdx++)
	RGseeds.IRandom(0,1000000);
NumAcceptedSamples = 0;
NumTestedSamples = 0;
DeltaAcceptedSamples = 0;
DeltaTestedSamples = 0;
for(RepIdx = RepStartIdx; RepIdx < (RepStartIdx + NumReps); RepIdx++)
	{
	while(NumAcceptedSamples < Num2Sample)
		{
		NumTestedSamples += 1;
		DeltaTestedSamples += 1;
		// let foreground know how we are doing
		if(NumTestedSamples > 0 && !(NumTestedSamples % 100000))
			{
#ifdef _WIN32
			AcquireSRWLockExclusive(&m_hCntsRwLock);
#else
			pthread_rwlock_wrlock(&m_hCntsRwLock);
#endif
			m_NumTestedSamples += DeltaTestedSamples;
			m_NumAcceptedSamples += DeltaAcceptedSamples;
#ifdef _WIN32
			ReleaseSRWLockExclusive(&m_hCntsRwLock);
#else
			pthread_rwlock_unlock(&m_hCntsRwLock);
#endif
			DeltaTestedSamples = 0;
			DeltaAcceptedSamples = 0;
			}

		ElIdx = RGseeds.IRandom(0,NumReadEls-1);
		pReadEl = &pReadEls[ElIdx];
		pTransCnt = &pTransCnts[pReadEl->TransID-1];	// first try at the transcript level 
		CurThres = RGseeds.IRandom(0,1000);
		if(pTransCnt->Thres < CurThres)
			continue;
		CurThres = RGseeds.IRandom(0,1000);		
		if(pReadEl->SampThres < CurThres)
			continue;

#ifdef _WIN32
		AcquireSRWLockExclusive(&m_hCntsRwLock);
#else
		pthread_rwlock_wrlock(&m_hCntsRwLock);
#endif
		pTransCnt->Counts[RepIdx] += 1;
		pReadEl->RandCnt += 1;
#ifdef _WIN32
		ReleaseSRWLockExclusive(&m_hCntsRwLock);
#else
		pthread_rwlock_unlock(&m_hCntsRwLock);
#endif
		NumAcceptedSamples += 1;
		DeltaAcceptedSamples += 1;
		}
	}

#ifdef _WIN32
AcquireSRWLockExclusive(&m_hCntsRwLock);
#else
pthread_rwlock_wrlock(&m_hCntsRwLock);
#endif
m_NumTestedSamples += DeltaTestedSamples;
m_NumAcceptedSamples += DeltaAcceptedSamples;
#ifdef _WIN32
ReleaseSRWLockExclusive(&m_hCntsRwLock);
#else
pthread_rwlock_unlock(&m_hCntsRwLock);
#endif
pPars->NumAcceptedSamples = NumAcceptedSamples;
pPars->NumTestedSamples = NumTestedSamples;
pPars->Rslt = 1;
#ifdef _WIN32
_endthreadex(0);
return(eBSFSuccess);
#else
pthread_exit(&pPars->Rslt);
#endif
}

// SortSimReadEls
// Sort simulated reads by Rand
int
SortSimReadEls(const void *arg1, const void *arg2)
{
tsReadEl *pEl1 = (tsReadEl *)arg1;
tsReadEl *pEl2 = (tsReadEl *)arg2;
if(pEl1->RandCnt < pEl2->RandCnt)
	return(-1);
if(pEl1->RandCnt > pEl2->RandCnt)
	return(1);
return(0);
}