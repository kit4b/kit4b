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
#include "SNPs2pgSNPs.h"

typedef struct TAG_sCodonXlate
	{
	uint8_t Idx;		// AminoAcid tRNA complement, aaa (0) ... ttt (63)
	uint8_t SynGroup;	// synonymous AminoAcid grouping - 0..20 including stop codon
	double HumanFreqPerK; // frequency (per 1000) at which codon observed in human
	char AminoAcid;		// single char AminoAcid representation
	bool bWobble;		// true if a wobble codon, with no cognate tRNA in human
	} tsCodonXlate;

typedef struct TAG_sAminoAcid
	{
	uint8_t SynGroup;	// synonymous amino acid grouping - 0..20 including stop codon
	char AminoAcid;		// single char symbol
	char szPeptide[4];	// three char symbol
	} tsAminoAcid;

static tsAminoAcid gAminoAcids[] = {
	{0,'A',"Ala"},
	{1,'C',"Cys"},
	{2,'D',"Asp"},
	{3,'E',"Glu"},
	{4,'F',"Phe"},
	{5,'G',"Gly"},
	{6,'H',"His"},
	{7,'I',"Ile"},
	{8,'K',"Lys"},
	{9,'L',"Leu"},
	{10,'M',"Met"},
	{11,'N',"Asn"},
	{12,'P',"Pro"},
	{13,'Q',"Gln"},
	{14,'R',"Arg"},
	{15,'S',"Ser"},
	{16,'T',"Thr"},
	{17,'V',"Val"},
	{18,'W',"Trp"},
	{19,'Y',"Tyr"},
	{20,'*',"Stp"}};


static tsCodonXlate gXlateAminoAcids[] = {		// all possible codon AminoAcids (*) slow with no cognate tRNA in human!
	{0,8,24.0,'K',false},	// aaa
	{1,11,19.5,'N',true},	// aac *
	{2,8,32.9,'K',false},	// aag
	{3,11,16.7,'N',false},	// aat

	{4,16,14.8,'T',false},	// aca
	{5,16,19.2,'T',false},	// acc
	{6,16,6.2,'T',false},	// acg	
	{7,16,12.8,'T',false},	// act

	{8,14,11.5,'R',false},	// aga
	{9,15,19.4,'S',false},	// agc 
	{10,14,11.4,'R',false},	// agg
	{11,15,11.9,'S',true},	// agt *

	{12,7,7.1,'I',false},	// ata
	{13,7,21.4,'I',false},	// atc
	{14,10,22.3,'M',false},	// atg	
	{15,7,15.7,'I',false},	// att

	{16,13,11.8,'Q',false},	// caa
	{17,6,14.9,'H',false},	// cac 
	{18,13,34.6,'Q',false},	// cag
	{19,6,10.4,'H',true},	// cat *

	{20,12,16.7,'P',false},	// cca
	{21,12,20.0,'P',true},	// ccc *
	{22,12,7.0,'P',false},	// ccg	
	{23,12,17.3,'P',false},	// cct

	{24,14,6.3,'R',false},	// cga
	{25,14,10.9,'R',true},	// cgc *
	{26,14,11.9,'R',false},	// cgg
	{27,14,4.7,'R',false},	// cgt

	{28,9,6.9,'L',false},	// cta
	{29,9,19.4,'L',true},	// ctc *
	{30,9,40.3,'L',false},	// ctg	
	{31,9,12.8,'L',false},	// ctt

	{32,3,29.0,'E',false},	// gaa
	{33,2,26.0,'D',false},	// gac 
	{34,3,40.8,'E',false},	// gag
	{35,2,22.3,'D',true},	// gat *

	{36,0,16.0,'A',false},	// gca
	{37,0,28.5,'A',false},	// gcc
	{38,0,7.6,'A',false},	// gcg	
	{39,0,18.6,'A',false},	// gct

	{40,5,16.3,'G',false},	// gga
	{41,5,22.8,'G',false},	// ggc 
	{42,5,16.4,'G',false},	// ggg
	{43,5,10.8,'G',false},	// ggt

	{44,17,7.0,'V',false},	// gta
	{45,17,14.6,'V',true},	// gtc *
	{46,17,28.9,'V',false},	// gtg	
	{47,17,10.9,'V',false},	// gtt

	{48,20,0.7,'*',false},	// taa
	{49,19,15.6,'Y',false},	// tac 
	{50,20,0.5,'*',false},	// tag
	{51,19,12.0,'Y',false},	// tat

	{52,15,11.7,'S',false},	// tca
	{53,15,17.4,'S',true},	// tcc *
	{54,15,4.5,'S',false},	// tcg	
	{55,15,14.6,'S',false},	// tct

	{56,20,1.3,'*',false},	// tga
	{57,1,12.2,'C',false},	// tgc 
	{58,18,12.8,'W',false},	// tgg
	{59,1,9.9,'C',false},	// tgt *

	{60,9,7.2,'L',false},	// tta
	{61,4,20.4,'F',false},	// ttc
	{62,9,12.6,'L',false},	// ttg	
	{63,4,16.9,'F',true}	// ttt *
	};

int Process(eModepgSSNP Mode,				// processing mode
			bool bAllowInferenced,			// true if both original kalign SNP calls and snpmarker inferenced SNP calls to be used, default is for kalign SNPs only
			double LocalSeqErrRate,			// local sequencing error rate in a 100bp window centered around the putative SNP site
			int MinCoverage,				// must be at least this coverage at SNP site
			double MinAlleleProp,			// putative allele must be at least this proportion of total site read coverage
			double PValueThres,				// only accept SNP alleles which have a PValue <= this threshold
			char* pszTrackName,				// track name
			char* pszAssemblyName,			// UCSC assembly name - for SARS-CoV-2 it is "wuhCor1"
			char* pszExperimentDescr,		// describes experiment
			etSSetOp SetOp,					// set operation on SetA and/or SetB				
			int NumInSetA,					// number of species/cultivars/isolates in SetA
			char** ppszSetA,				// names of those cspecies/cultivars/isolates which are in SetA
			int NumInSetB,					// number of species/cultivars/isolates in SetB
			char **ppszSetB,				// names of those cspecies/cultivars/isolates in SetB
			char *pszGFFFile,				// general feature format file - identifies start/end loci of features in assembly, i.e. transcripts etc.
			char* pszSNPFile,				// load SNP calls from this CSV file, can be either SNPs generated by kalign or snpmarkers
			char* pszOutFile);				// output SNPs to this UCSC Personal Genome SNP format file, or VCF4.1 if extension is '.vcf'



#ifdef _WIN32
int SNPs2pgSNPs(int argc, char* argv[])
{
	// determine my process name
_splitpath(argv[0], NULL, NULL, gszProcName, NULL);
#else
int
SNPs2pgSNPs(int argc, char** argv)
{
	// determine my process name
	CUtility::splitpath((char*)argv[0], NULL, gszProcName);
#endif
	int Len;
	int iFileLogLevel;			// level of file diagnostics
	int iScreenLogLevel;		// level of file diagnostics
	char szLogFile[_MAX_PATH];	// write diagnostics to this file
	int Rslt = 0;   			// function result code >= 0 represents success, < 0 on failure
	
	int NumberOfProcessors;		// number of installed CPUs

	eModepgSSNP PMode;			// processing mode
	bool bAllowInferenced;		// true if both original kalign SNP calls and snpmarker inferenced SNP calls to be used, default is for kalign SNPs only 
	double LocalSeqErrRate;		// default local sequencing error rate in a 100bp window centered around the putative SNP site
	int MinCoverage;			// must be at least this coverage at SNP site
	double MinAlleleProp;			// putative allele must be at least this proportion of total site read coverage
	double PValueThres;				// only accept SNP alleles which have a PValue <= this threshold
	char szInSNPsFile[_MAX_PATH];	// processing this kalign SNP calls or snpmarker CSV file
	char szGFFFile[_MAX_PATH];			// general feature format file - identifies start/end loci of features in assembly, i.e. transcripts etc.
	char szOutpgSNPsFile[_MAX_PATH];		// output in UCSC Personal Genome SNP format to this file

	int ChkAIdx;					// used when checking for duplicate names in sets
	int ChkBIdx;					// used when checking for duplicate names in sets
	etSSetOp SetOp;					// set operation on SetA and/or SetB				
	int NumInSetA;					// number of species/cultivars/isolates in SetA
	char* pszSetA[cMaxSSetMembers];	// names of those cspecies/cultivars/isolates which are in SetA
	int NumInSetB;					// number of species/cultivars/isolates in SetB
	char* pszSetB[cMaxSSetMembers];	// names of those cspecies/cultivars/isolates in SetB

	char szTrackName[cMaxDatasetSpeciesChrom + 1];	// UCSC track name
	char szAssemblyName[cMaxDatasetSpeciesChrom + 1];	// UCSC assembly name
	char szExperimentDescr[cMaxDatasetSpeciesChrom + 1];	// describes experiment


	struct arg_lit* help = arg_lit0("h", "help", "print this help and exit");
	struct arg_lit* version = arg_lit0("v", "version,ver", "print version information and exit");
	struct arg_int* FileLogLevel = arg_int0("f", "FileLogLevel", "<int>", "Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug");
	struct arg_file* LogFile = arg_file0("F", "log", "<file>", "diagnostics log file");

	struct arg_int* pmode = arg_int0("m", "mode", "<int>", "processing mode: 0 kalign SNP calls, 1: snpmarkers calls");

	struct arg_lit* allowinferenced = arg_lit0("A", "allowinferenced", "allow snpmarker inferenced SNP sites, default is use kalign SNP sites only");
	struct arg_file* insnps = arg_file1("i", "insnps", "<file>", "input called SNPs or snpmarkers from this CSV file");
	struct arg_file* gff = arg_file0("g", "gtf", "<file>", "input transcriptomic features from this virus GTF file, WARNING: targeting sense strand only transcripts, with a max of 100 CDS features");
	struct arg_int* setop = arg_int0("V", "setop", "<int>", "SetA/B operator 0: none, 1: (A union B), 2: (A intersect B), 3: complement of (A union B), 4: complement of (A intersect B) 5: (A subtract B)");
	struct arg_str* seta = arg_strn("s", "seta", "<seta>", 0, cMaxSSetMembers, "name of snpmarkers species/cultivar/strain in SetA");
	struct arg_str* setb = arg_strn("S", "setb", "<setb>", 0, cMaxSSetMembers, "name of snpmarkers species/cultivar/strain in SetB");

	struct arg_file* outpgsnps = arg_file1("o", "outpgsnps", "<file>", "output in UCSC pgSNP format to this file (if extn is '.vcf' then format will be VCF 4.1)");
	struct arg_str* experimentdescr = arg_str1("e", "experiment", "<str>", "UCSC pgSNP experiment description");
	struct arg_str* assemblyname = arg_str1("a", "assembly", "<str>", "UCSC assembly name");
	struct arg_str* trackname = arg_str1("t", "track", "<str>", "UCSC track name");

	struct arg_dbl* locerrrate = arg_dbl0("l", "locerrrate", "<dbl>", "default local sequencing error rate in a 100bp window centered around the putative SNP site (default 0.02, range 0.001 .. 0.25)");
	struct arg_dbl* pvaluethres = arg_dbl0("p", "pvalue", "<dbl>", "SNP PValue threshold (default 0.10, range 0.005..0.25)");
	struct arg_dbl* minalleleprop = arg_dbl0("P", "minalleleprop", "<dbl>", "SNP minimum allele proportion of loci site coverage (default 0.025, range 0.01..0.95)");
	struct arg_int* mincoverage = arg_int0("c", "mincoverage", "<int>", "SNP minimum loci site coverage (default 20, range 1..10000)");

	struct arg_end* end = arg_end(200);

	void* argtable[] = { help,version,FileLogLevel,LogFile,
						pmode,allowinferenced,locerrrate,pvaluethres,minalleleprop,mincoverage,insnps,gff,setop,seta,setb,outpgsnps,experimentdescr,assemblyname,trackname,end };

	char** pAllArgs;
	int argerrors;
	argerrors = CUtility::arg_parsefromfile(argc, (char**)argv, &pAllArgs);
	if (argerrors >= 0)
		argerrors = arg_parse(argerrors, pAllArgs, argtable);

	/* special case: '--help' takes precedence over error reporting */
	if (help->count > 0)
	{
		printf("\nCAUTION: This processing module is specifically targeting analysis of virus (e.g. SARS-CoV-2) sized assembly SNP calls");
		printf("\n%s %s %s, Version %s\nOptions ---\n", gszProcName, gpszSubProcess->pszName, gpszSubProcess->pszFullDescr, kit4bversion);
		arg_print_syntax(stdout, argtable, "\n");
		arg_print_glossary(stdout, argtable, "  %-25s %s\n");
		printf("\nNote: Parameters can be entered into a parameter file, one parameter per line.");
		printf("\n      To invoke this parameter file then precede its name with '@'");
		printf("\n      e.g. %s %s @myparams.txt\n", gszProcName, gpszSubProcess->pszName);
		printf("\nPlease report any issues regarding usage of %s at https://github.com/kit4b/issues\n\n", gszProcName);
		return(1);
	}

	/* special case: '--version' takes precedence error reporting */
	if (version->count > 0)
	{
		printf("\n%s %s Version %s\n", gszProcName, gpszSubProcess->pszName, kit4bversion);
		return(1);
	}

	if (!argerrors)
	{
		if (FileLogLevel->count && !LogFile->count)
		{
			printf("\nError: FileLogLevel '-f%d' specified but no logfile '-F<logfile>\n'", FileLogLevel->ival[0]);
			exit(1);
		}

		iScreenLogLevel = iFileLogLevel = FileLogLevel->count ? FileLogLevel->ival[0] : eDLInfo;
		if (iFileLogLevel < eDLNone || iFileLogLevel > eDLDebug)
		{
			printf("\nError: FileLogLevel '-l%d' specified outside of range %d..%d\n", iFileLogLevel, eDLNone, eDLDebug);
			exit(1);
		}

		if (LogFile->count)
		{
			strncpy(szLogFile, LogFile->filename[0], _MAX_PATH);
			szLogFile[_MAX_PATH - 1] = '\0';
		}
		else
		{
			iFileLogLevel = eDLNone;
			szLogFile[0] = '\0';
		}

		// now that log parameters have been parsed then initialise diagnostics log system
		if (!gDiagnostics.Open(szLogFile, (etDiagLevel)iScreenLogLevel, (etDiagLevel)iFileLogLevel, true))
		{
			printf("\nError: Unable to start diagnostics subsystem\n");
			if (szLogFile[0] != '\0')
				printf(" Most likely cause is that logfile '%s' can't be opened/created\n", szLogFile);
			exit(1);
		}

		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Subprocess %s Version %s starting", gpszSubProcess->pszName, kit4bversion);
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "CAUTION: This processing module is specifically targeting analysis of virus (e.g. SARS-CoV-2) sized assemblies SNP calls");
		gExperimentID = 0;
		gProcessID = 0;
		gProcessingID = 0;
		szExperimentDescr[0] = '\0';

		PMode = pmode->count ? (eModepgSSNP)pmode->ival[0] : eMpgSSNPKalign;
		if (PMode < eMpgSSNPKalign || PMode > eMpgSSNPmarkers)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Processing mode '-m%d' specified outside of range %d..%d\n", PMode, eMpgSSNPKalign, (int)eMpgSSNPmarkers);
			exit(1);
			}

		bAllowInferenced = allowinferenced->count ? true : false;

		strncpy(szExperimentDescr, experimentdescr->sval[0], sizeof(szExperimentDescr) - 1);
		szExperimentDescr[sizeof(szExperimentDescr) - 1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(szExperimentDescr);
		CUtility::ReduceWhitespace(szExperimentDescr);
		if ((Len =(int)strlen(szExperimentDescr)) < 3 || Len > 80)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Expected UCSC experiment description length to be in range of 3..80\n");;
			exit(1);
			}

		strncpy(szAssemblyName, assemblyname->sval[0], sizeof(szAssemblyName) - 1);
		szAssemblyName[sizeof(szAssemblyName) - 1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(szAssemblyName);
		CUtility::ReduceWhitespace(szAssemblyName);
		if ((Len = (int)strlen(szAssemblyName)) < 3 || Len > 50)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Expected UCSC assembly name length to be in range of 3..50\n");
			exit(1);
			}

		strncpy(szTrackName, trackname->sval[0], sizeof(szTrackName) - 1);
		szTrackName[sizeof(szTrackName) - 1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(szTrackName);
		CUtility::ReduceWhitespace(szTrackName);
		if ((Len=(int)strlen(szTrackName)) < 3 || Len > 50)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error: Expected UCSC track name length to be in range of 3..50\n");
			exit(1);
			}

		LocalSeqErrRate = locerrrate->count ? locerrrate->dval[0] : 0.02;
		if(LocalSeqErrRate <= 0.001)
			LocalSeqErrRate = 0.001;
		else
			if(LocalSeqErrRate > 0.25)
				LocalSeqErrRate = 0.25;

		PValueThres = pvaluethres->count ? pvaluethres->dval[0] : 0.10;
		if(PValueThres < 0.005)
			PValueThres = 0.005;
		else
			if(PValueThres > 0.25)
				PValueThres = 0.25;

		MinAlleleProp = minalleleprop->count ? minalleleprop->dval[0] : 0.025;
		if(MinAlleleProp <= 0.01)
			MinAlleleProp = 0.01;
		else
			if(MinAlleleProp > 0.99)
				MinAlleleProp = 0.95;

		MinCoverage = mincoverage->count ? mincoverage->ival[0] : 20;
		if(MinCoverage <= 1)
			MinCoverage = 1;
		else
			if(MinCoverage > 10000)
				MinCoverage = 10000;

		// show user current resource limits
#ifndef _WIN32
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Resources: %s", CUtility::ReportResourceLimits());
#endif

#ifdef _WIN32
		SYSTEM_INFO SystemInfo;
		GetSystemInfo(&SystemInfo);
		NumberOfProcessors = SystemInfo.dwNumberOfProcessors;
#else
		NumberOfProcessors = sysconf(_SC_NPROCESSORS_CONF);
#endif

		strcpy(szInSNPsFile, insnps->filename[0]);
		CUtility::TrimQuotedWhitespcExtd(szInSNPsFile);
		if (szInSNPsFile[0] == '\0')
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "No input SNPs file specified");
			exit(1);
			}

		if(gff->count)
			{
			strcpy(szGFFFile, gff->filename[0]);
			CUtility::TrimQuotedWhitespcExtd(szGFFFile);
			if(szGFFFile[0] == '\0')
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "GTF file specified but after whitespace trimming no file specified");
				exit(1);
				}
			}
		else
			szGFFFile[0] = '\0';

		strncpy(szOutpgSNPsFile, outpgsnps->filename[0],sizeof(szOutpgSNPsFile));
		szOutpgSNPsFile[sizeof(szOutpgSNPsFile) - 1] = '\0';
		CUtility::TrimQuotedWhitespcExtd(szOutpgSNPsFile);
		if (szOutpgSNPsFile[0] == '\0')
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "No UCSC Personal Genome or VCF SNP format output file specified");
			exit(1);
			}
		

		SetOp = setop->count ? (etSSetOp)setop->ival[0] : eSSONone;
		if (SetOp < eSSONone || SetOp >= eSSOPlaceHolder)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Set operation must be in range 0..%d", eSSOPlaceHolder-1);
			exit(1);
			}

		memset(pszSetA,0,sizeof(pszSetA));
		NumInSetA = 0;
		if (SetOp != eSSONone && PMode == eMpgSSNPmarkers && seta->count)
			{ 
			for (int Idx = 0; NumInSetA < cMaxSSetMembers && Idx < seta->count; Idx++)
				{
				pszSetA[Idx] = NULL;
				if (pszSetA[NumInSetA] == NULL)
					pszSetA[NumInSetA] = new char[_MAX_PATH];
				strncpy(pszSetA[NumInSetA], seta->sval[Idx], _MAX_PATH);
				pszSetA[NumInSetA][_MAX_PATH - 1] = '\0';
				CUtility::TrimQuotedWhitespcExtd(pszSetA[NumInSetA]);
				if (pszSetA[NumInSetA][0] != '\0')
					NumInSetA++;
				}
			}

		if(SetOp != eSSONone && NumInSetA == 0)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Set operation requested but SetA is empty");
			exit(1);
			}

		if(NumInSetA)
			{
			// check for duplicate names
			for (ChkAIdx = 0; ChkAIdx < NumInSetA - 1; ChkAIdx++)
				{
				for (ChkBIdx = ChkAIdx + 1; ChkBIdx < NumInSetA; ChkBIdx++)
					{
					if (!stricmp(pszSetA[ChkAIdx], pszSetA[ChkBIdx]))
						{
						// have a duplicate name in Set
						gDiagnostics.DiagOut(eDLFatal, gszProcName, "SetA set contains duplicate name '%s'", pszSetA[ChkAIdx]);
						exit(1);
						}
					}
				}
			}

		memset(pszSetB, 0, sizeof(pszSetB));
		NumInSetB = 0;
		if (SetOp != eSSONone && PMode == eMpgSSNPmarkers && setb->count)
			{
			for (int Idx = 0; NumInSetB < cMaxSSetMembers && Idx < setb->count; Idx++)
				{
				pszSetB[Idx] = NULL;
				if (pszSetB[NumInSetB] == NULL)
					pszSetB[NumInSetB] = new char[_MAX_PATH];
				strncpy(pszSetB[NumInSetB], setb->sval[Idx], _MAX_PATH);
				pszSetB[NumInSetB][_MAX_PATH - 1] = '\0';
				CUtility::TrimQuotedWhitespcExtd(pszSetB[NumInSetB]);
				if (pszSetB[NumInSetB][0] != '\0')
					NumInSetB++;
				}
			}

		if (NumInSetB)
			{
			// check for duplicate names
			for (ChkAIdx = 0; ChkAIdx < NumInSetB - 1; ChkAIdx++)
				{
				for (ChkBIdx = ChkAIdx + 1; ChkBIdx < NumInSetB; ChkBIdx++)
					{
						if (!stricmp(pszSetB[ChkAIdx], pszSetB[ChkBIdx]))
						{
							// have a duplicate name in set
							gDiagnostics.DiagOut(eDLFatal, gszProcName, "SetB contains duplicate name '%s'", pszSetB[ChkAIdx]);
							exit(1);
						}
					}
				}

		// check for names common to both SetA and SetB
			for (ChkAIdx = 0; ChkAIdx < NumInSetA; ChkAIdx++)
				{
				for (ChkBIdx = 0; ChkBIdx < NumInSetB; ChkBIdx++)
					{
					if (!stricmp(pszSetA[ChkAIdx], pszSetB[ChkBIdx]))
						{
						gDiagnostics.DiagOut(eDLFatal, gszProcName, "SetA and SetB contain common name '%s'", pszSetA[ChkAIdx]);
						exit(1);
						}
					}
				}
			}

		switch (SetOp) {
			case eSSONone:						// no set operation
				break;
			case eSSOUnion:						// union of SetA and SetB
				if (!NumInSetB)
					{
					gDiagnostics.DiagOut(eDLWarn, gszProcName, "Set union operation requested but SetB is empty");
					gDiagnostics.DiagOut(eDLWarn, gszProcName, "Treating SetA as the union");
					}
				break;
			case eSSOIntersect:					// intersect of SetA and SetB
				if(!NumInSetB)
					{
					gDiagnostics.DiagOut(eDLFatal, gszProcName, "Set intersect operation requested but SetB is empty");
					exit(1);
					}
				break;
			case eSSOCplUnion:					// complement of union SetA or SetB
				if (!NumInSetB)
					{
					gDiagnostics.DiagOut(eDLWarn, gszProcName, "Set complement union operation requested but SetB is empty");
					gDiagnostics.DiagOut(eDLWarn, gszProcName, "Treating SetA as the union");
					}
				break;

			case eSSOCplIntersect:				// complement of intersect SetA and SetB
				if (!NumInSetB)
					{
					gDiagnostics.DiagOut(eDLFatal, gszProcName, "Set complement intersect operation requested but SetB is empty");
					exit(1);
					}
				break;

			case eSSOSubtract:					// A subtract B
				if (!NumInSetB)
					{
					gDiagnostics.DiagOut(eDLWarn, gszProcName, "Set subtract union operation requested but SetB is empty");
					gDiagnostics.DiagOut(eDLWarn, gszProcName, "Using all of SetA");
					}
				break;
			}

		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Processing parameters:");
		const char* pszDescr;
		switch (PMode) {
			case eMpgSSNPKalign:
				pszDescr = "input SNPs are those generated by kalign";
				break;

			case eMpgSSNPmarkers:
				pszDescr = "input SNPs are those generated by snpmarkers";
				break;
			}

		gDiagnostics.DiagOutMsgOnly(eDLInfo, "SNP to UCSC Personal Genome SNP or VCF4.1 format conversion : '%s'", pszDescr);

		gDiagnostics.DiagOutMsgOnly(eDLInfo, "Using SNP sites from : '%s'", bAllowInferenced ? "both 'kalign' called and 'snpmarkers' inferenced" : "only 'kaligned' called");

		if(szGFFFile[0] != '\0')
			{
			gDiagnostics.DiagOutMsgOnly(eDLInfo,"Load transcriptome features from: '%s'", szGFFFile);
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "WARNING: targeting sense strand only transcripts, with a max of 100 CDS features");
			}
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "Local sequencing error rate : '%f'", LocalSeqErrRate);
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "SNP minimum loci site coverage : '%d'", MinCoverage);
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "SNP minimum allele proportion of loci site coverage : '%f'", MinAlleleProp);

		gDiagnostics.DiagOutMsgOnly(eDLInfo, "SNP allele PValue threshold : '%f'", PValueThres);
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "UCSC assembly name : '%s'", szAssemblyName);
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "UCSC track name : '%s'", szTrackName);
		gDiagnostics.DiagOutMsgOnly(eDLInfo, "UCSC pgSNP experiment description : '%s'", szExperimentDescr);

		gDiagnostics.DiagOutMsgOnly(eDLInfo, "Input SNPs file : '%s'", szInSNPsFile);

		// check file extension, if '.vcf' then generate output formated for vcf instead of the default pgSNP format
		int Len;
		if ((Len = (int)strlen(szOutpgSNPsFile)) >= 5 && !stricmp(&szOutpgSNPsFile[Len - 4], ".vcf"))
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "Output to '%s' is in VCF 4.1 format", szOutpgSNPsFile);
		else
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "Output to '%s' is in UCSC pgSNP format", szOutpgSNPsFile);

		switch (SetOp) {
			case eSSONone:
				pszDescr = "no Set operation";
				break;
			case eSSOUnion:
				pszDescr = "union of SetA and SetB";
				break;
			case eSSOIntersect:
				pszDescr = "intersect of SetA and SetB";
				break;
			case eSSOCplUnion:
				pszDescr = "complement of union SetA and SetB";
				break;
			case eSSOCplIntersect:
				pszDescr = "complement of intersect SetA and SetB";
				break;
			case eSSOSubtract:
				pszDescr = "SetA subtract SetB";
				break;
			}

		gDiagnostics.DiagOutMsgOnly(eDLInfo, "Set operation : '%s'", pszDescr);

		if(NumInSetA)
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "Number of species/cultivars/isolates in SetA : '%d'", NumInSetA);

		if (NumInSetB)
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "Number of species/cultivars/isolates in SetB : '%d'", NumInSetB);


		if (szExperimentDescr[0] != '\0')
			gDiagnostics.DiagOutMsgOnly(eDLInfo, "Experiment description: %s", szExperimentDescr);

#ifdef _WIN32
		SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
		gStopWatch.Start();
		Rslt = 0;
		Rslt = Process(PMode,					// processing mode
					bAllowInferenced,			// true if both original kalign SNP calls and snpmarker inferenced SNP calls to be used, default is for kalign SNPs only
					LocalSeqErrRate,			// local sequencing error rate in a 100bp window centered around the putative SNP site
					MinCoverage,				// must be at least this coverage at SNP site
					MinAlleleProp,				// putative allele must be at least this proportion of total site read coverage
					PValueThres,				// only accept SNP alleles which have a PValue <= this threshold
					szTrackName,				// track name
					szAssemblyName,				// UCSC assembly name - for SARS-CoV-2 it is "wuhCor1"
					szExperimentDescr,			// describes experiment
					SetOp,						// set operation on SetA and/or SetB				
					NumInSetA,					// number of species/cultivars/isolates in SetA
					pszSetA,					// names of those cspecies/cultivars/isolates which are in SetA
					NumInSetB,					// number of species/cultivars/isolates in SetB
					pszSetB,					// names of those cspecies/cultivars/isolates in SetB
					szGFFFile,					// general feature format file - identifies start/end loci of features in assembly, i.e. transcripts etc.
					szInSNPsFile,				// load SNP calls from these CSV files (can be either SNPs generated by kalign or snpmarkers
					szOutpgSNPsFile);			// output SNPs to this UCSC Personal Genome SNP format file
		Rslt = Rslt >= 0 ? 0 : 1;
		gStopWatch.Stop();

		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Exit code: %d Total processing time: %s", Rslt, gStopWatch.Read());
		exit(Rslt);
	}
	else
	{
		printf("\n%s %s %s, Version %s\n", gszProcName, gpszSubProcess->pszName, gpszSubProcess->pszFullDescr, kit4bversion);
		arg_print_errors(stdout, end, gszProcName);
		arg_print_syntax(stdout, argtable, "\nUse '-h' to view option and parameter usage\n");
		exit(1);
	}
return 0;
}

int Process(eModepgSSNP Mode,				// processing mode
	bool bAllowInferenced,			// true if both original kalign SNP calls and snpmarker inferenced SNP calls to be used, default is for kalign SNPs only
	double LocalSeqErrRate,			// local sequencing error rate in a 100bp window centered around the putative SNP site
	int MinCoverage,				// must be at least this coverage at SNP site
	double MinAlleleProp,			// putative allele must be at least this proportion of total site read coverage
	double PValueThres,				// only accept SNP alleles which have a PValue <= this threshold
	char* pszTrackName,				// track name
	char* pszAssemblyName,			// UCSC assembly name - for SARS-CoV-2 it is "wuhCor1"
	char* pszExperimentDescr,		// describes experiment
	etSSetOp SetOp,					// set operation on SetA and SetB	
	int NumInSetA,					// number of species/cultivars/isolates in SetA
	char** ppszSetA,				// names of those cspecies/cultivars/isolates which are in SetA
	int NumInSetB,					// number of species/cultivars/isolates in SetB
	char** ppszSetB,				// names of those cspecies/cultivars/isolates in SetB
	char* pszGFFFile,				// general feature format file - identifies start/end loci of features in assembly, i.e. transcripts etc.
	char* pszSNPFile,				// load SNP calls from this CSV file, can be either SNPs generated by kalign or snpmarkers
	char* pszOutFile)				// output SNPs to this UCSC Personal Genome SNP format file
{
int Rslt;
CSNPs2pgSNPs *pSNPs2pgSNPs;
if ((pSNPs2pgSNPs = new CSNPs2pgSNPs) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to instantiate instance of CSNPs2pgSNPs");
	return(eBSFerrInternal);
	}

Rslt = pSNPs2pgSNPs->Process(Mode, bAllowInferenced, LocalSeqErrRate,MinCoverage, MinAlleleProp,PValueThres, pszTrackName, pszAssemblyName,pszExperimentDescr, SetOp, NumInSetA, ppszSetA,NumInSetB, ppszSetB, pszGFFFile, pszSNPFile, pszOutFile);
if(pSNPs2pgSNPs != NULL)
	delete pSNPs2pgSNPs;
return(Rslt);
}


CSNPs2pgSNPs::CSNPs2pgSNPs(void)
{
m_pSNPSites = NULL;
m_pCSV = NULL;
m_pszOutBuff = NULL;
m_pszSiteDistBuff = NULL;
m_pGapCntDist = NULL;
m_pSharedSNPsDist = NULL;
m_pIndividualSNPsDist = NULL;
m_pGFFFile = NULL;
m_pFeatures = NULL;
m_pIsolateFeatSNPs = NULL;
m_pMatrix = NULL;
m_pSimilaritiesMatrix = NULL;
m_hOutFile = -1;
m_hOutSiteDist = -1; 
Reset();
}


CSNPs2pgSNPs::~CSNPs2pgSNPs(void)
{
if(m_pGFFFile != NULL)
	delete m_pGFFFile;
if(m_pCSV != NULL)
	delete m_pCSV;
if(m_hOutFile !=-1)
	close(m_hOutFile);
if (m_pszOutBuff != NULL)
	delete []m_pszOutBuff;
if(m_pszSiteDistBuff != NULL)
	delete []m_pszSiteDistBuff;
if(m_pMatrix != NULL)
	delete []m_pMatrix;
if(m_pSimilaritiesMatrix == NULL)
	delete []m_pSimilaritiesMatrix;
if(m_pSNPSites != NULL)
	delete m_pSNPSites;
if(m_pIsolateFeatSNPs != NULL)
	{
#ifdef _WIN32
	free(m_pIsolateFeatSNPs);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pIsolateFeatSNPs != MAP_FAILED)
		munmap(m_pIsolateFeatSNPs, m_AllocdIsolateFeatSNPsMem);
#endif
	}
if(m_pGapCntDist != NULL)
	delete []m_pGapCntDist;
if(m_pSharedSNPsDist != NULL)
	delete []m_pSharedSNPsDist;
if(m_pIndividualSNPsDist != NULL)
	delete []m_pIndividualSNPsDist;
if(m_pFeatures != NULL)
	delete []m_pFeatures;
}

void
CSNPs2pgSNPs::Reset(void)
{
if(m_pGFFFile != NULL)
	{
	delete m_pGFFFile;
	m_pGFFFile = NULL;
	}

if(m_pFeatures != NULL)
	{
	delete []m_pFeatures;
	m_pFeatures = NULL;
	}



if(m_pCSV != NULL)
	{
	delete m_pCSV;
	m_pCSV = NULL;
	}
if (m_hOutFile != -1)
	{
	if (m_OutBuffOffs && m_pszOutBuff != NULL)
		CUtility::RetryWrites(m_hOutFile, m_pszOutBuff, m_OutBuffOffs);
#ifdef _WIN32
	_commit(m_hOutFile);
#else
	fsync(m_hOutFile);
#endif
	close(m_hOutFile);
	m_hOutFile = -1;
	}

if(m_hOutSiteDist != -1)
	{
	if(m_SiteDistOffs && m_pszSiteDistBuff != NULL)
		CUtility::RetryWrites(m_hOutSiteDist, m_pszSiteDistBuff, m_SiteDistOffs);
#ifdef _WIN32
	_commit(m_hOutSiteDist);
#else
	fsync(m_hOutSiteDist);
#endif
	close(m_hOutSiteDist);
	m_hOutSiteDist = -1;
	}



if(m_pszOutBuff != NULL)
	{
	delete []m_pszOutBuff;
	m_pszOutBuff = NULL;
	}

if(m_pszSiteDistBuff != NULL)
	{
	delete[]m_pszSiteDistBuff;
	m_pszSiteDistBuff = NULL;
	}

if(m_pSNPSites != NULL)
	{
	delete m_pSNPSites;
	m_pSNPSites = NULL;
	}

if(m_pGapCntDist != NULL)
	{
	delete[]m_pGapCntDist;
	m_pGapCntDist = NULL;
	}
if(m_pSharedSNPsDist != NULL)
	{
	delete[]m_pSharedSNPsDist;
	m_pSharedSNPsDist = NULL;
	}
if(m_pIndividualSNPsDist != NULL)
	{
	delete[]m_pIndividualSNPsDist;
	m_pIndividualSNPsDist = NULL;
	}

if(m_pIsolateFeatSNPs != NULL)
	{
#ifdef _WIN32
	free(m_pIsolateFeatSNPs);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pIsolateFeatSNPs != MAP_FAILED)
		munmap(m_pIsolateFeatSNPs, m_AllocdIsolateFeatSNPsMem);
#endif
	m_pIsolateFeatSNPs = NULL;
	}

if(m_pMatrix != NULL)
	{
	delete[]m_pMatrix;
	m_pMatrix = NULL;
	}

if(m_pSimilaritiesMatrix != NULL)
	{
	delete []m_pSimilaritiesMatrix;
	m_pSimilaritiesMatrix = NULL;
	}

m_SortOrder = 0;
m_MatrixCols = 0;
m_MatrixRows = 0;
m_SimilaritiesMatrixCols = 0;
m_SimilaritiesMatrixRows = 0;

m_NumUniqueSiteLoci = 0;
m_NumCultivars = 0;

m_NumChromNames = 0;
m_szChromNames[0] = '\0';
m_szChromIdx[0] = 0;
m_NxtszChromIdx = 0;

m_NumReadsetNames = 0;
m_szReadsetNames[0] = '\0';
m_szReadsetIdx[0] = 0;
m_NxtszReadsetIdx = 0;

m_szGFFFile[0] = '\0';
m_AllocdOutBuff = 0;
m_OutBuffOffs = 0;
m_szInSNPsFile[0] = '\0';
m_szOutSiteDistFile[0] = '\0';

m_AllocdSNPSites = 0;
m_NumSNPs = 0;
m_SNPID = 0;

m_MaxSNPgap = 0;
m_MaxSharedSNPs = 0;
m_SiteDistOffs = 0;
m_AllocdSiteDistBuff = 0;

m_IsolateFeatSNPs = 0;
m_AllocdIsolateFeatSNPs = 0;
m_AllocdIsolateFeatSNPsMem = 0;

m_NumFeatures = 0;
m_AllocFeatures = 0;
m_bRptHdr = true;
}



int 
CSNPs2pgSNPs::ProcessSnpmarkersSNPs(void)
{
	int Rslt;
	int NumFields;
	uint32_t CultivarIdx;
	tsSCultivar* pCultivar;
	uint32_t PrevSNPloci;
	uint32_t EstNumSNPs;
	uint32_t NumSNPsParsed;
	uint32_t RowNumber;
	bool bMarkerFormat;
	int ExpNumFields;
	uint32_t NumSpeciesWithCnts;
	tsSNPSSite SNPSite;
	tsSNPSSite* pSNPSite;

	uint32_t NumDiracs;
	uint32_t NumNonDiracs;
	uint32_t DiracCnts[5];

	bool bRprtHdr = true;
	EstNumSNPs = m_pCSV->EstNumRows();
	NumSNPsParsed = 0;
	RowNumber = 0;
	ExpNumFields = 0;
	m_NumCultivars = 0;
	bMarkerFormat = false;
	pSNPSite = &SNPSite;
	PrevSNPloci = 0xffffffff;
	NumDiracs = 0;
	NumNonDiracs = 0;
	while((Rslt = m_pCSV->NextLine()) > 0)				// onto next line containing fields
		{
		memset(DiracCnts,0,sizeof(DiracCnts));
		RowNumber += 1;
		NumFields = m_pCSV->GetCurFields();
		if(NumFields == cAlignNumSSNPfields || NumFields == cAlignNumSSNPXfields || NumFields < cSSNPMarkerNumFields)
		{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Number of fields (%d) at row %d does not match expected fields for a snpmarkers file '%s'", NumFields, RowNumber, m_szInSNPsFile, RowNumber);
			Reset();
			return(eBSFerrFieldCnt);
		}
		if(ExpNumFields && ExpNumFields != NumFields)
		{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Inconsistency in number of fields, previously %d but %d fields parsed from '%s' near line %d", ExpNumFields, NumFields, m_szInSNPsFile, RowNumber);
			Reset();
			return(eBSFerrFieldCnt);
		}
		if(!ExpNumFields)
		{
			m_NumCultivars = (NumFields - 4) / 9;
			if(((m_NumCultivars * 9) + 4) != NumFields)
			{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Invalid number (%d) of fields parsed from '%s' near line %d, was this file generated by 'ngskit4b kalign/snpmarkers'?", NumFields, m_szInSNPsFile, RowNumber);
				Reset();
				return(eBSFerrFieldCnt);
			}
			if(m_NumCultivars > cMaxCultivars)
			{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Too many cultivars (%d) in '%s', max allowed is %d", m_NumCultivars, m_szInSNPsFile, cMaxCultivars);
				Reset();
				return(eBSFerrFieldCnt);
			}
			bMarkerFormat = true;
			ExpNumFields = NumFields;
		}

		if(RowNumber == 1)
		{
			if(m_pCSV->IsLikelyHeaderLine())
			{
				if(!NumSNPsParsed && bMarkerFormat)
				{
				// parse out target assembly name against which alignments were made and the cultivar names
					char* pSrc;
					char* pDst;
					char Chr;
					int Len;

					m_pCSV->GetText(1, &pSrc);
					pDst = m_TargAssemblyName;
					Len = 0;
					while(Len < sizeof(m_TargAssemblyName) - 1 && (Chr = *pSrc++) && Chr != ':')
					{
						*pDst++ = Chr;
						Len++;
						*pDst = '\0';
					}
					memset(m_Cultivars, 0, sizeof(m_Cultivars));
					for(CultivarIdx = 0; CultivarIdx < m_NumCultivars; CultivarIdx++)
					{
						m_pCSV->GetText(5 + (CultivarIdx * 9), &pSrc);
						pDst = m_Cultivars[CultivarIdx].szName;
						Len = 0;
						while(Len < sizeof(m_Cultivars[CultivarIdx].szName) - 1 && (Chr = *pSrc++) && Chr != ':')
						{
							*pDst++ = Chr;
							Len++;
							*pDst = '\0';
						}
						if(!m_NumInSetA)
							m_Cultivars[CultivarIdx].flgInSetA = false;
						else
						{
							int WIdx;
							for(WIdx = 0; WIdx < m_NumInSetA; WIdx++)
							{
								if(!stricmp(m_ppszSetA[WIdx], m_Cultivars[CultivarIdx].szName))
									break;
							}
							if(WIdx == m_NumInSetA)
								m_Cultivars[CultivarIdx].flgInSetA = false;
							else
								m_Cultivars[CultivarIdx].flgInSetA = true;
						}
						if(!m_NumInSetB)
							m_Cultivars[CultivarIdx].flgInSetB = false;
						else
						{
							int WIdx;
							for(WIdx = 0; WIdx < m_NumInSetB; WIdx++)
							{
								if(!stricmp(m_ppszSetB[WIdx], m_Cultivars[CultivarIdx].szName))
									break;
							}
							if(WIdx == m_NumInSetB)
								m_Cultivars[CultivarIdx].flgInSetB = false;
							else
								m_Cultivars[CultivarIdx].flgInSetB = true;
						}
					}
				}
			}
			else
			{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Expected CSV file '%s' first line to be a header line with  fully quoted field names", m_szInSNPsFile);
				Reset();
				return(eBSFerrFieldCnt);
			}
			continue;
		}
		memset(pSNPSite, 0, sizeof(*pSNPSite));
		strcpy(pSNPSite->szSeqReadSet, m_ReadsetIdentifier);
		if((Rslt = pSNPSite->ReadsetID = AddReadset(pSNPSite->szSeqReadSet)) < 1)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Failed generating a unique readset name identifier '%s'", pSNPSite->szSeqReadSet);
			Reset();
			return(Rslt);
			}

		char* pszTxt;
		m_pCSV->GetText(1, &pszTxt);			// get ":TargSeq"
		strncpy(pSNPSite->szChrom, pszTxt, sizeof(pSNPSite->szChrom));
		pSNPSite->szChrom[sizeof(pSNPSite->szChrom) - 1] = '\0';
		if((Rslt = pSNPSite->ChromID = AddChrom(pSNPSite->szChrom)) < 1)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Failed generating a unique chromosome name identifier '%s'", pSNPSite->szChrom);
			Reset();
			return(Rslt);
			}
		m_pCSV->GetInt(2, (int*)&pSNPSite->SNPLoci);			// get "Loci"
		m_pCSV->GetText(3, &pszTxt);			// get "TargBase"
		switch(*pszTxt) {
			case 'a': case 'A':
				pSNPSite->RefBase = eBaseA;
				break;
			case 'c': case 'C':
				pSNPSite->RefBase = eBaseC;
				break;
			case 'g': case 'G':
				pSNPSite->RefBase = eBaseG;
				break;
			case 't': case 'T': case 'u': case 'U':	// U in case RNA alignments..
				pSNPSite->RefBase = eBaseT;
				break;
			case 'n': case 'N':				// unlikely to have a SNP against an indeterminate base but you never know...
				pSNPSite->RefBase = eBaseN;
				break;
			default:
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Expected SNP TargBase ('%s') to be one of 'ACGTN' in CSV file '%s' near line %d", pszTxt, m_szInSNPsFile, RowNumber);
				Reset();
				return(eBSFerrFieldCnt);
		}
		m_pCSV->GetInt(4, (int *)&NumSpeciesWithCnts);			// get "NumSpeciesWithCnts"
		if(NumSpeciesWithCnts < 1 || NumSpeciesWithCnts > m_NumCultivars)
		{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Expected NumSpeciesWithCnts (%d) to be between 1 and %d in CSV file '%s' near line %d", NumSpeciesWithCnts, m_NumCultivars, m_szInSNPsFile, RowNumber);
			Reset();
			return(eBSFerrNumRange);
		}

		int NumNotInferenced = 0;
		int FieldIdx = 5;
		pCultivar = m_Cultivars;
		for(CultivarIdx = 0; CultivarIdx < m_NumCultivars; CultivarIdx++, pCultivar++)
		{
			m_pCSV->GetText(FieldIdx++, &pszTxt);			// get SNP call source, either directly called 'S' or inferenced 'I'
			if(*pszTxt == 'I' || *pszTxt == 'i')
				pCultivar->flgSNPInferenced = true;
			else
			{
				pCultivar->flgSNPInferenced = false;
				NumNotInferenced += 1;
			}
			m_pCSV->GetText(FieldIdx++, &pszTxt);		// get cultivar ":Base"
			switch(*pszTxt) {
				case 'a': case 'A':
					pCultivar->CalledBase = eBaseA;
					break;
				case 'c': case 'C':
					pCultivar->CalledBase = eBaseC;
					break;
				case 'g': case 'G':
					pCultivar->CalledBase = eBaseG;
					break;
				case 't': case 'T': case 'u': case 'U':	// U in case RNA alignments..
					pCultivar->CalledBase = eBaseT;
					break;
				case 'n': case 'N':				// unlikely to have a SNP against an indeterminate base but you never know...
					pCultivar->CalledBase = eBaseN;
					break;
				default:
					gDiagnostics.DiagOut(eDLFatal, gszProcName, "Expected SNP TargBase ('%s') to be one of 'ACGTN' in CSV file '%s' near line %d", pszTxt, m_szInSNPsFile, RowNumber);
					Reset();
					return(eBSFerrFieldCnt);
			}
			m_pCSV->GetDouble(FieldIdx++, &pCultivar->Score);			// get ":LSER"
			m_pCSV->GetUint(FieldIdx++, &pCultivar->TotBaseCnts);		// get ":BaseCntTot"
			m_pCSV->GetUint(FieldIdx++, &pCultivar->BaseCnts[0]);		// get ":BaseCntA"
			m_pCSV->GetUint(FieldIdx++, &pCultivar->BaseCnts[1]);		// get ":BaseCntC"
			m_pCSV->GetUint(FieldIdx++, &pCultivar->BaseCnts[2]);		// get ":BaseCntG"
			m_pCSV->GetUint(FieldIdx++, &pCultivar->BaseCnts[3]);		// get ":BaseCntT"
			m_pCSV->GetUint(FieldIdx++, &pCultivar->BaseCnts[4]);		// get ":BaseCntN"
		}

		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Evaluating SNP alleles at %s:%d with %d non-inferenced sites", pSNPSite->szChrom, pSNPSite->SNPLoci, NumNotInferenced);

		// have all cultivars and their counts at the SNP loci
		// iterate the cultivars and for every cultivar determine which bases are allelic and count occurrences over all cultivars
		bool bHaveSNP;
		double PValue;
		int NumInSetASnps;
		int NumInSetBSnps;


		bHaveSNP = false;
		NumInSetASnps = 0;
		NumInSetBSnps = 0;
		memset(pSNPSite->BaseCnts, 0, sizeof(pSNPSite->BaseCnts));
		pSNPSite->ClassifySite = 0x00;
		pCultivar = m_Cultivars;
		for(CultivarIdx = 0; CultivarIdx < m_NumCultivars; CultivarIdx++, pCultivar++)
		{
			if(pCultivar->flgSNPInferenced)
				pCultivar->NumSitesInferenced += 1;
			else
				pCultivar->NumSitesCalled += 1;

			switch(m_SetOp) {
				case eSSONone:					// no Set operation, all (universe) isolates were SNP called
					break;						// just accept!

				case eSSOUnion:					// union of SetA and SetB
					if(!(pCultivar->flgInSetA || pCultivar->flgInSetB)) // no interest in the universe of all SNPs, just SetA and SetB
						continue;
					break;

				case eSSOIntersect:				// intersect of SetA and SetB
					if(!(pCultivar->flgInSetA && pCultivar->flgInSetB))
						continue;
					break;

				case eSSOCplUnion:				// complement of union SetA and SetB
					if(pCultivar->flgInSetA || pCultivar->flgInSetB)
						continue;
					break;						// have a SNP which is not in either Set

				case eSSOCplIntersect:			// complement of intersect SetA and SetB, all (universe) isolates were processed for SNPs
					if(pCultivar->flgInSetA && pCultivar->flgInSetB) // looking for universe SNPs which are not in both Sets
						continue;
					break;						// have a universe SNP which is not in both Sets

				case eSSOSubtract:				// SetA subtract SetB
					if(!pCultivar->flgInSetA)
						continue;
					break;

				default:
					break;
			}

			if(pCultivar->flgSNPInferenced && !m_bAllowInferenced)
				continue;

			pSNPSite->TotBaseCnts += 1;

			pCultivar->TotMismatches = 0;
			pCultivar->TotBaseCnts = 0;
			for(int Idx = 0; Idx < 4; Idx++)
				{
				pCultivar->TotBaseCnts += pCultivar->BaseCnts[Idx];
				if(Idx == pSNPSite->RefBase)
					continue;
				pCultivar->TotMismatches += pCultivar->BaseCnts[Idx];
				}


			if(pCultivar->TotBaseCnts < m_MinCoverage)		// have to have at least this many bases covering the loci
				continue;

			if(m_MinAlleleProp > 0.0 && (pCultivar->TotBaseCnts == 0 || (pCultivar->TotMismatches / (double)pCultivar->TotBaseCnts) < m_MinAlleleProp))
				continue;

			PValue = 1.0 - m_Stats.Binomial(pCultivar->TotBaseCnts, pCultivar->TotMismatches, pCultivar->Score);
			if(PValue > m_PValueThres)
				continue;

			bool bIsDirac = false;
			bHaveSNP = true;
			for(int Idx = 0; Idx < 4; Idx++)
				{
				if(Idx == pSNPSite->RefBase)
					continue;
				if((pCultivar->BaseCnts[Idx] / (double)pCultivar->TotBaseCnts) >= m_DiracThres)
					{
					bIsDirac = true;
					DiracCnts[Idx] += 1;
					}

				if(m_MinAlleleProp > 0.0 && (pCultivar->BaseCnts[Idx] == 0 || (((double)pCultivar->BaseCnts[Idx] * 4) / (double)pCultivar->TotBaseCnts) < m_MinAlleleProp))
					continue;

				pSNPSite->BaseCnts[Idx].Cnts += 1;
				PValue = min(0.5, 1.0 - m_Stats.Binomial(pCultivar->TotBaseCnts, pCultivar->BaseCnts[Idx], pCultivar->Score));
				pSNPSite->BaseCnts[Idx].QScore += PValue;
				}
			if(bIsDirac)
				NumDiracs++;
			else
				NumNonDiracs++;

			if(pCultivar->flgInSetA)
				NumInSetASnps++;
			if(pCultivar->flgInSetB)
				NumInSetBSnps++;
			pSNPSite->TotMismatches += 1;
			pCultivar->NumSitesAccepted += 1;
			}

		if(bHaveSNP)
			{
			for(int Idx = 0; Idx < 4; Idx++)
				if(pSNPSite->BaseCnts[Idx].Cnts > 1)
					pSNPSite->BaseCnts[Idx].QScore /= pSNPSite->BaseCnts[Idx].Cnts;

			// Note: some of these set operations have been prefiltered such that
			// only isolates relevant to the operation have been processed
			switch(m_SetOp) {
				case eSSONone:					// no Set operation, all (universe) isolates were SNP called
					break;						// just accept!

				case eSSOUnion:					// union of SetA and SetB, only isolates marked as in SetA or SetB were SNP called
					if(!(NumInSetASnps || NumInSetBSnps))	// require at least one SNP in either Set
						bHaveSNP = false;
					break;						// have at least one so accept!

				case eSSOIntersect:					// intersect of SetA and SetB, only isolates marked as in SetA or SetB were SNP called
					if(!(NumInSetASnps && NumInSetBSnps)) // require at least one SNP in both Sets
						bHaveSNP = false;;
					break;						// SNPs in both so accept!

				case eSSOCplUnion:				// complement of union SetA and SetB, all (universe) isolates were SNP called
					if(NumInSetASnps || NumInSetBSnps) // looking for universe SNPs which are not in either Set
						bHaveSNP = false;;
					break;						// have a SNP which is not in either Set


				case eSSOCplIntersect:			// complement of intersect SetA and SetB, all (universe) isolates were processed for SNPs
					if(NumInSetASnps && NumInSetBSnps) // looking for universe SNPs which are not in both Sets
						bHaveSNP = false;;
					break;						// have a universe SNP which is not in both Sets

				case eSSOSubtract:				// SNPs in A but which are not in B, only isolates marked as in SetA or SetB were processed for SNPs
					if(!NumInSetASnps || NumInSetBSnps)
						bHaveSNP = false;;
					break;						// have SNPs in SetA but none in SetB
				}
			}

		if(bHaveSNP)	// if still accepted as being a called SNP
			{
			// classify site loci
			for(int Idx = 0; Idx < 4; Idx++)
				{
				if(Idx == pSNPSite->RefBase)
					continue;
				if((pCultivar->BaseCnts[Idx] / (double)pCultivar->TotBaseCnts) >= m_DiracThres)
					pSNPSite->ClassifySite = 0x03;
				else
					if(pSNPSite->ClassifySite == 0 && (pCultivar->BaseCnts[Idx] / (double)pCultivar->TotBaseCnts) >= (m_DiracThres/2))
						pSNPSite->ClassifySite = 0x02;
				}
			if(pSNPSite->ClassifySite == 0x00)	// if not previously classified then can only be a minor allele
				pSNPSite->ClassifySite = 0x01;	
			if(PrevSNPloci != 0xffffffff)
				{
				uint32_t CurGap;
				CurGap = pSNPSite->SNPLoci - PrevSNPloci;
				if(CurGap < 1)					// must be at least one otherwise upstream generation of loci failed
					{
					gDiagnostics.DiagOut(eDLFatal, gszProcName, "Overlap in loci discovered at loci: %u...",pSNPSite->SNPLoci);
					return(eBSFerrLocField);
					}
				if(CurGap > cMaxGapBetweenSNPS)
					CurGap = cMaxGapBetweenSNPS;
				m_pGapCntDist[CurGap - 1] += 1;
				if(m_MaxSNPgap < CurGap)
					m_MaxSNPgap = CurGap;
				}
			PrevSNPloci = pSNPSite->SNPLoci;

			uint32_t CurShared = pSNPSite->TotMismatches;
			if(CurShared > cMaxSharedSNPS)
				CurShared = cMaxSharedSNPS;
			m_pSharedSNPsDist[CurShared - 1] += 1;
			if(m_MaxSharedSNPs < CurShared)
				m_MaxSharedSNPs = CurShared;

			memcpy(pSNPSite->DiracCnts,DiracCnts,sizeof(pSNPSite->DiracCnts));
			m_NumSNPs++;
			pSNPSite->bSNPPlaceholder = false;
			}
		else
			{
			pSNPSite->bSNPPlaceholder = true;
			pSNPSite->ClassifySite = 0;
			}

		pSNPSite->SNPId = ++m_SNPID;		// all SNPs have an identifier even if only a placeholder
		if(m_bRptHdr && m_NumSNPs)
			{
			Report(true, pSNPSite);
			m_bRptHdr = false;
			}
		Report(false, pSNPSite);
	}

if((Rslt = WriteOutFile(true)) < eBSFSuccess)
	return(Rslt);


if(m_SiteDistOffs)
	{
	CUtility::RetryWrites(m_hOutSiteDist, m_pszSiteDistBuff, m_SiteDistOffs);
	m_SiteDistOffs = 0;
	}

#ifdef _WIN32
_commit(m_hOutSiteDist);
#else
fsync(m_hOutSiteDist);
#endif
close(m_hOutSiteDist);
m_hOutSiteDist = -1;



	m_MaxIndividualSNPs = 0;
	pCultivar = m_Cultivars;
	for(CultivarIdx = 0; CultivarIdx < m_NumCultivars; CultivarIdx++, pCultivar++)
		{
		if(pCultivar->NumSitesAccepted >= cMaxIndividualSNPS)
			pCultivar->NumSitesAccepted = cMaxIndividualSNPS;
		m_pIndividualSNPsDist[pCultivar->NumSitesAccepted] += 1;
		if(pCultivar->NumSitesAccepted > m_MaxIndividualSNPs)
			m_MaxIndividualSNPs = pCultivar->NumSitesAccepted;
		}
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Discovered %u unique SNP loci with longest gap between loci of %u", m_NumSNPs, m_MaxSNPgap);
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Maximum number of accepted SNP loci in any individual cultivar was %u", m_MaxIndividualSNPs);
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Maximum number of accepted SNP loci shared between cultivars was %u", m_MaxSharedSNPs);
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Proportion of individual cultivars SNPs characterised as dirac is: %f", NumDiracs/ (double)((uint64_t)NumDiracs + NumNonDiracs));
return(0);
}



int
CSNPs2pgSNPs::ProcessKalignSNPs(void)
{
int Rslt;
int NumFields;
UINT32 EstNumSNPs;
UINT32 NumSNPsParsed;
UINT32 RowNumber;
bool bMarkerFormat;
int ExpNumFields;
int BaseCntsIdx;
tsSBaseCnts* pBaseCnts;
int CntRef;
double PValue;
bool bHaveSNP;
tsSNPSSite SNPSite;
tsSNPSSite *pSNPSite;

EstNumSNPs = m_pCSV->EstNumRows();
NumSNPsParsed = 0;
RowNumber = 0;
ExpNumFields = 0;
bMarkerFormat = false;
while ((Rslt = m_pCSV->NextLine()) > 0)				// onto next line containing fields
	{
	memset(&SNPSite,0,sizeof(SNPSite));
	pSNPSite = &SNPSite;
	RowNumber += 1;
	NumFields = m_pCSV->GetCurFields();
	if (ExpNumFields && ExpNumFields != NumFields)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Inconsistency in number of fields, previously %d but %d fields parsed from '%s' near line %d", ExpNumFields, NumFields, m_szInSNPsFile, RowNumber);
		Reset();
		return(eBSFerrFieldCnt);
		}
	if (!ExpNumFields)
		{
		if (!(NumFields == cAlignNumSSNPfields || NumFields == cAlignNumSSNPXfields))								// must be exactly this many if 'ngskit4b kalign' format
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "%d fields parsed from '%s' near line %d, expected %d, was file generated by 'kit4b kalign'?", NumFields, m_szInSNPsFile, RowNumber, cAlignNumSSNPfields);
			Reset();
			return(eBSFerrFieldCnt);
			}
		bMarkerFormat = false;
		ExpNumFields = NumFields;
		}

	if (RowNumber == 1)
		{
		if (m_pCSV->IsLikelyHeaderLine())
			continue;
		else
			gDiagnostics.DiagOut(eDLWarn, gszProcName, "Expected CSV file '%s' first line to be a header line with  fully quoted field names", m_szInSNPsFile); // found instances where header line was missing when aligned on amazon ec2
		}

	if (!NumSNPsParsed)
		{
		// parse out target assembly name against which alignments were made and the cultivar names
		char* pSrc;
		char* pDst;
		char Chr;
		int Len;

		m_pCSV->GetText(3, &pSrc);
		pDst = m_TargAssemblyName;
		Len = 0;
		while (Len < sizeof(m_TargAssemblyName) - 1 && (Chr = *pSrc++))
			{
			*pDst++ = Chr;
			Len++;
			*pDst = '\0';
			}
		memset(&m_Cultivars[0], 0, sizeof(m_Cultivars[0]));
		strcpy(m_Cultivars[0].szName, m_TargAssemblyName);
		}

	char* pszTxt;

	memset(pSNPSite, 0, sizeof(SNPSite));
	strcpy(SNPSite.szSeqReadSet, m_ReadsetIdentifier);
	pBaseCnts = &pSNPSite->BaseCnts[0];
	for (BaseCntsIdx = 0; BaseCntsIdx < 5; BaseCntsIdx++, pBaseCnts++)
		pBaseCnts->Base = (etSeqBase)BaseCntsIdx;


	int SNPlen;
	m_pCSV->GetInt(7, &SNPlen);			// check that we are processing SNPs!
	if (SNPlen != 1)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Expected SNP CSV file '%s' to only contain 1 base SNPs, 'Len' = %d near line %d", m_szInSNPsFile, SNPlen, RowNumber);
		Reset();
		return(eBSFerrFieldCnt);
		}
	m_pCSV->GetUint(11, &pSNPSite->TotBaseCnts);		// get "Bases"
	if(pSNPSite->TotBaseCnts < m_MinCoverage)
		continue;
	pSNPSite->ReadsetID = m_ReadsetID;
	m_pCSV->GetUint(1, &pSNPSite->ReadsetSiteId);
	m_pCSV->GetText(4, &pszTxt);			// get "Chrom"
	strncpy(pSNPSite->szChrom, pszTxt, sizeof(pSNPSite->szChrom));
	pSNPSite->szChrom[sizeof(pSNPSite->szChrom) - 1] = '\0';
	if((Rslt = pSNPSite->ChromID = AddChrom(pSNPSite->szChrom)) < 1)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Failed generating a unique chromosome name identifier '%s'", pSNPSite->szChrom);
		Reset();
		return(Rslt);
		}

	m_pCSV->GetUint(5, &pSNPSite->SNPLoci);			// get "StartLoci"

	m_pCSV->GetUint(12, &pSNPSite->TotMismatches);		// get "Mismatches"
	CntRef = pSNPSite->TotBaseCnts - pSNPSite->TotMismatches;
	m_pCSV->GetUint(14, &pSNPSite->BaseCnts[0].Cnts);	// get "MMBaseA"
	m_pCSV->GetUint(15, &pSNPSite->BaseCnts[1].Cnts);	// get "MMBaseC"
	m_pCSV->GetUint(16, &pSNPSite->BaseCnts[2].Cnts);	// get "MMBaseG"
	m_pCSV->GetUint(17, &pSNPSite->BaseCnts[3].Cnts);	// get "MMBaseT"
	m_pCSV->GetUint(18, &pSNPSite->BaseCnts[4].Cnts);	// get "MMBaseN"
	m_pCSV->GetText(13, &pszTxt);		// get "RefBase"
	switch (*pszTxt) {
		case 'a': case 'A':
			pSNPSite->RefBase = eBaseA;
			pSNPSite->BaseCnts[0].Cnts = CntRef;
			break;
		case 'c': case 'C':
			pSNPSite->RefBase = eBaseC;
			pSNPSite->BaseCnts[1].Cnts = CntRef;
			break;
		case 'g': case 'G':
			pSNPSite->RefBase = eBaseG;
			pSNPSite->BaseCnts[2].Cnts = CntRef;
			break;
		case 't': case 'T': case 'u': case 'U':	// U in case RNA alignments..
			pSNPSite->RefBase = eBaseT;
			pSNPSite->BaseCnts[3].Cnts = CntRef;
			break;
		case 'n': case 'N':				// unlikely to have a SNP against an indeterminate base but you never know...
			pSNPSite->RefBase = eBaseN;
			pSNPSite->BaseCnts[4].Cnts = CntRef;
			break;
		default:
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Expected SNP RefBase ('%s') to be one of 'ACGTN' in CSV file '%s' near line %d", pszTxt, m_szInSNPsFile, RowNumber);
			Reset();
			return(eBSFerrFieldCnt);
		}
	m_pCSV->GetDouble(19, &pSNPSite->LocalSeqErrRate);	// get "BackgroundSubRate"

	if(m_MinAlleleProp > 0.0 && (pSNPSite->TotBaseCnts == 0 || (pSNPSite->TotMismatches / (double)pSNPSite->TotBaseCnts) < m_MinAlleleProp))
		continue;

	PValue = 1.0 - m_Stats.Binomial(pSNPSite->TotBaseCnts, pSNPSite->TotMismatches, pSNPSite->LocalSeqErrRate);
	if(PValue > m_PValueThres)
		continue;

	bHaveSNP = false;
	pSNPSite->ClassifySite = 0;
	for (int Idx = 0; Idx < 4; Idx++)
		{
		pSNPSite->OrigBaseCnts[Idx] = pSNPSite->BaseCnts[Idx].Cnts;
		if(m_MinAlleleProp > 0.0 && (pSNPSite->BaseCnts[Idx].Cnts == 0 || (((double)pSNPSite->BaseCnts[Idx].Cnts * 4) / (double)pSNPSite->TotBaseCnts) < m_MinAlleleProp))
			{
			pSNPSite->BaseCnts[Idx].Cnts = 0;
			pSNPSite->BaseCnts[Idx].QScore = 0.0;
			continue;
			}

		PValue = 1.0 - m_Stats.Binomial(pSNPSite->TotBaseCnts, pSNPSite->BaseCnts[Idx].Cnts, pSNPSite->LocalSeqErrRate);
		if(PValue > m_PValueThres)
			{
			pSNPSite->BaseCnts[Idx].Cnts = 0;
			pSNPSite->BaseCnts[Idx].QScore = 0.0;
			}
		else
			{
			if(Idx != pSNPSite->RefBase)   // only non-ref alleles can be characterised
				{
				if((pSNPSite->BaseCnts[Idx].Cnts / (double)pSNPSite->TotBaseCnts) >= m_DiracThres)
					{
					pSNPSite->DiracCnts[Idx] = 1;
					pSNPSite->ClassifySite = 0x03;
					}
				else
					if(pSNPSite->ClassifySite < 0x02)
						{
						if((pSNPSite->BaseCnts[Idx].Cnts / (double)pSNPSite->TotBaseCnts) >= (m_DiracThres/2))
							pSNPSite->ClassifySite = 0x02;
						else
							pSNPSite->ClassifySite = 0x01;
						}
				}
			pSNPSite->BaseCnts[Idx].QScore = PValue; 
			if(Idx != pSNPSite->RefBase)
				bHaveSNP = true;
			}
		}
	if(!bHaveSNP)
		{
		pSNPSite->ClassifySite = 0x00;
		continue;
		}
	NumSNPsParsed++;

	tsSummaryFeatCnts* pFeature = m_pFeatures;
	eSeqBase Base;
	int FeatureIdx;
	int FrameShift;
	int FieldIdx;
	int PolypeptideCnts[21];
	char szRefCodon[4];
	char *pszRefCodon;
	int RefCodonIdx;
	int AlignCodonIdx;
	int AlignCodons[64];

	pSNPSite->bInFeature = false;
	pSNPSite->FeatureIdx = 0;
	if(m_NumFeatures > 0)
		{
				// locate the feature in which this SNP is located
		for(FeatureIdx = 0; FeatureIdx < m_NumFeatures; FeatureIdx++, pFeature++)
			{
			if(stricmp(pSNPSite->szChrom, pFeature->szChrom) && stricmp(pSNPSite->szChrom, "NC_045512.2"))
				continue;
			if(pSNPSite->SNPLoci >= pFeature->FeatStart && pSNPSite->SNPLoci <= pFeature->FeatEnd)
				{
				memset(&pFeature->SiteFeatCnts,0,sizeof(pFeature->SiteFeatCnts));
				break;
				}
			}
		if(FeatureIdx < m_NumFeatures)
			{
			pSNPSite->bInFeature = true;
			pSNPSite->FeatureIdx = FeatureIdx;
			pFeature->TotSNPs += 1;
			if((pSNPSite->ClassifySite & 0x03) == 3)
				{
				pFeature->SiteFeatCnts.IsDirac = true;
				for(FieldIdx = 0; FieldIdx <= 4; FieldIdx++)
					pFeature->DiracCnts[FieldIdx] += pSNPSite->DiracCnts[FieldIdx];
				}
			if(NumFields == cAlignNumSSNPXfields)	// if SNP file is an extended SNP file then can utilise the additional codon frame shift fields present in the extension
				{
				// which frame shift to use?
				FrameShift = (pSNPSite->SNPLoci - pFeature->FeatStart) % 3;	// need to determine frame relative to feature start loci
				if(FrameShift == 1)
					FrameShift = 0;
				else
					if(FrameShift == 0)
						FrameShift = 1;
				FieldIdx = cAlignNumSSNPfields + 1 + (65 * FrameShift);
				m_pCSV->GetText(FieldIdx++, &pszRefCodon);
				strncpy(szRefCodon, pszRefCodon,sizeof(szRefCodon));
				RefCodonIdx = 0;
				for(int Idx = 0; Idx <= 2; Idx++)
					{
					switch(szRefCodon[Idx])
						{
							case 'a': case 'A':
								Base = eBaseA;
								break;

							case 'c': case 'C':
								Base = eBaseC;
								break;

							case 'g': case 'G':
								Base = eBaseG;
								break;

							case 't': case 'T':
								Base = eBaseT;
								break;

							case 'n': case 'N':
								break;
						}
					RefCodonIdx = (RefCodonIdx <<= 2) | Base;
					}
				memset(PolypeptideCnts,0,sizeof(PolypeptideCnts));
				int SumCodonCnts = 0;
				for(AlignCodonIdx = 0; AlignCodonIdx < 64; AlignCodonIdx++)
					{
					m_pCSV->GetInt(FieldIdx++, &AlignCodons[AlignCodonIdx]);
					SumCodonCnts += AlignCodons[AlignCodonIdx];
					}


				bool bIsSynonymous = false;
				bool bIsNonSynonymous = false;
				bool bIsWobble = false;
				pFeature->SiteFeatCnts.RefSynGroup = gXlateAminoAcids[RefCodonIdx].SynGroup;
				for(AlignCodonIdx = 0; AlignCodonIdx < 64; AlignCodonIdx++)
					{
					// what count threshold should be used before treating as noise and thus counts to be ignored???
					// here I'm requiring at least 2 codon instances covering site with a PValue <= 2*m_PValueThres which is slightly better than guessing
					if(AlignCodons[AlignCodonIdx] < 2)
						continue;
					PValue = 1.0 - m_Stats.Binomial(SumCodonCnts, AlignCodons[AlignCodonIdx], pSNPSite->LocalSeqErrRate);
					if(PValue > (m_PValueThres * 2))
						continue;

					// accepting as codon change
					pFeature->FromAminoAcidChanges[gXlateAminoAcids[RefCodonIdx].SynGroup] += 1;
					pFeature->ToAminoAcidChanges[gXlateAminoAcids[RefCodonIdx].SynGroup][gXlateAminoAcids[AlignCodonIdx].SynGroup] += 1;
					pFeature->SiteFeatCnts.ToAminoAcidChanges[gXlateAminoAcids[AlignCodonIdx].SynGroup] += AlignCodons[AlignCodonIdx];

					PolypeptideCnts[gXlateAminoAcids[AlignCodonIdx].SynGroup] += AlignCodons[AlignCodonIdx];
					if(gXlateAminoAcids[AlignCodonIdx].SynGroup == gXlateAminoAcids[RefCodonIdx].SynGroup)
						{
						bIsSynonymous = true;
						pFeature->NumSiteSynonCodons++;
						pFeature->SiteFeatCnts.TotSynonymous += AlignCodons[AlignCodonIdx];
						// attempt to see if can pickup on codon bias adaptation towards humans more abundant tRNAs
						if(gXlateAminoAcids[AlignCodonIdx].HumanFreqPerK > gXlateAminoAcids[RefCodonIdx].HumanFreqPerK)
							{
							pFeature->NumPosBiasedCodons++;
							pFeature->SiteFeatCnts.NumPosBiasedCodons += AlignCodons[AlignCodonIdx];
							}
						else
							if(gXlateAminoAcids[AlignCodonIdx].HumanFreqPerK < gXlateAminoAcids[RefCodonIdx].HumanFreqPerK)
								{
								pFeature->NumNegBiasedCodons++;
								pFeature->SiteFeatCnts.NumNegBiasedCodons += AlignCodons[AlignCodonIdx];
								}
						}
					else
						{
						if(gXlateAminoAcids[AlignCodonIdx].bWobble)
							{
							bIsWobble = true;
							pFeature->SiteFeatCnts.TotWobble += AlignCodons[AlignCodonIdx];
							}
						else
							{
							bIsNonSynonymous = true;
							pFeature->SiteFeatCnts.TotNonSynonymous += AlignCodons[AlignCodonIdx];
							}
						}
					}

				if(bIsSynonymous && !(bIsNonSynonymous || bIsWobble))
					pFeature->TotExclSynonymous += 1;
				if(bIsNonSynonymous && !(bIsSynonymous || bIsWobble))
					pFeature->TotExclNonSynonymous += 1;
				if(bIsWobble && !(bIsSynonymous || bIsNonSynonymous))
					pFeature->TotExclWobble += 1;

				if(bIsSynonymous)
					pFeature->TotSynonymous += 1;
				if(bIsNonSynonymous)
					pFeature->TotNonSynonymous += 1;
				if(bIsWobble)
					pFeature->TotWobble += 1;
				}
			}
		}

	m_SNPID+=1;
	pSNPSite->SNPId = m_NumSNPs = m_SNPID;

	if(m_bRptHdr)
		Report(true, pSNPSite);
	m_bRptHdr = false;
	Report(false, pSNPSite);

	tsIsolateFeatSNPs* pIsolateFeatSNPs;
	if(m_pIsolateFeatSNPs == NULL)
		{
		size_t memreq = (size_t)(sizeof(tsIsolateFeatSNPs) * cAllocNumIsolateFeatSNPs);

#ifdef _WIN32
		m_pIsolateFeatSNPs = (tsIsolateFeatSNPs*)malloc(memreq);	// initial and perhaps the only allocation
		if(m_pIsolateFeatSNPs == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "IsolateFeatSNPs: Initial memory allocation of %lld bytes - %s", (INT64)memreq, strerror(errno));
			Reset();
			return(eBSFerrMem);
			}
#else
		// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
		m_pIsolateFeatSNPs = (tsIsolateFeatSNPs *)mmap(NULL, memreq, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
		if(m_pIsolateFeatSNPs == MAP_FAILED)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "IsolateFeatSNPs: Memory allocation of %lld bytes through mmap()  failed - %s", (INT64)memreq, strerror(errno));
			m_pIsolateFeatSNPs = NULL;
			Reset();
			return(eBSFerrMem);
			}
#endif
		m_AllocdIsolateFeatSNPsMem = memreq;
		m_IsolateFeatSNPs = 0;
		m_AllocdIsolateFeatSNPs = cAllocNumIsolateFeatSNPs;
		}
	else
		{
		if((m_IsolateFeatSNPs + 100) > m_AllocdIsolateFeatSNPs)
			{
			size_t memreq = m_AllocdIsolateFeatSNPsMem + (cAllocNumIsolateFeatSNPs * sizeof(tsIsolateFeatSNPs));
				
#ifdef _WIN32
			pIsolateFeatSNPs = (tsIsolateFeatSNPs*)realloc(m_pIsolateFeatSNPs, memreq);
#else
			pIsolateFeatSNPs = (tsIsolateFeatSNPs*)mremap(m_pIsolateFeatSNPs, m_AllocdIsolateFeatSNPsMem, memreq, MREMAP_MAYMOVE);
			if(pIsolateFeatSNPs == MAP_FAILED)
				pIsolateFeatSNPs = NULL;
#endif
			if(pIsolateFeatSNPs == NULL)
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "IsolateFeatSNPs: Memory re-allocation to %lld bytes - %s", memreq, strerror(errno));
				return(eBSFerrMem);
				}
			m_pIsolateFeatSNPs = pIsolateFeatSNPs;
			m_AllocdIsolateFeatSNPs += cAllocNumIsolateFeatSNPs;
			m_AllocdIsolateFeatSNPsMem = memreq;
			}
		}
	pIsolateFeatSNPs = &m_pIsolateFeatSNPs[m_IsolateFeatSNPs++];
	memset(pIsolateFeatSNPs,0,sizeof(tsIsolateFeatSNPs));
	if(m_NumFeatures > 0 && pFeature != NULL)
		memcpy(&pIsolateFeatSNPs->SiteFeatCnts,pFeature,sizeof(tsSummaryFeatCnts));
	memcpy(&pIsolateFeatSNPs->SNPSSite, pSNPSite,sizeof(tsSNPSSite));
	}
return(NumSNPsParsed);
}

int		// returned chrom identifier, < 1 if unable to accept this chromosome name
CSNPs2pgSNPs::AddChrom(char* pszChrom) // associate unique identifier with this chromosome name
{
int ChromNameIdx;
int ChromNameLen;

// iterate over all known chroms in case this chrom to add is a duplicate
for(ChromNameIdx = 0; ChromNameIdx < m_NumChromNames; ChromNameIdx++)
	if(!stricmp(pszChrom, &m_szChromNames[m_szChromIdx[ChromNameIdx]]))
		return(ChromNameIdx + 1);

// chrom is not a duplicate
ChromNameLen = (int)strlen(pszChrom);
if((m_NxtszChromIdx + ChromNameLen + 1) > (int)sizeof(m_szChromNames))
	return(eBSFerrMaxEntries);
if(m_NumChromNames == cMaxChromNames)
	return(eBSFerrMaxEntries);

m_szChromIdx[m_NumChromNames] = m_NxtszChromIdx;
strcpy(&m_szChromNames[m_NxtszChromIdx], pszChrom);
m_NxtszChromIdx += ChromNameLen + 1;
return(++m_NumChromNames);
}

char* 
CSNPs2pgSNPs::LocateChrom(int ChromID)
{
if(ChromID < 1 || ChromID > m_NumChromNames)
	return(NULL);
return(&m_szChromNames[m_szChromIdx[ChromID-1]]);
}

int		// returned readset identifier, < 1 if unable to accept this readset name
CSNPs2pgSNPs::AddReadset(char* pszReadset) // associate unique identifier with this readset name
{
	int ReadsetIdx;
	int ReadsetNameLen;

	// iterate over all known readset in case this readset to add is a duplicate
	for(ReadsetIdx = 0; ReadsetIdx < m_NumReadsetNames; ReadsetIdx++)
		if(!stricmp(pszReadset, &m_szReadsetNames[m_szReadsetIdx[ReadsetIdx]]))
			return(ReadsetIdx + 1);

	// readset is not a duplicate
	ReadsetNameLen = (int)strlen(pszReadset);
	if((m_NxtszReadsetIdx + ReadsetNameLen + 1) > sizeof(m_szReadsetNames))
		return(eBSFerrMaxEntries);
	if(m_NumReadsetNames == cMaxReadsetNames)
		return(eBSFerrMaxEntries);

	m_szReadsetIdx[m_NumReadsetNames] = m_NxtszReadsetIdx;
	strcpy(&m_szReadsetNames[m_NxtszReadsetIdx], pszReadset);
	m_NxtszReadsetIdx += ReadsetNameLen + 1;
	return(++m_NumReadsetNames);
}

char* // returned ptr to readset name 
CSNPs2pgSNPs::LocateReadset(int ReadsetID) // readset name identifier
{
	if(ReadsetID < 1 || ReadsetID > m_NumReadsetNames)
		return(NULL);
	return(&m_szReadsetNames[m_szReadsetIdx[ReadsetID - 1]]);
}



// following function is experimental, enabling reporting in three different file formats
// pgSNP, VCF4.1, or CSV
int
CSNPs2pgSNPs::Report(bool bHeader,		// if true then header line(s) to be generated, otherwise SNPs or alleles
					 tsSNPSSite* pSNPSite)
{
int Rslt;
char Base;
double SumQScores;
int NumAlleles;
int Depth;
int PercentScore;
char szChrom[cMaxDatasetSpeciesChrom + 1];

if(bHeader)
	{
	int StartLoci = pSNPSite->SNPLoci - 500;
	if(StartLoci < 0)
		StartLoci = 0;
	int EndLoci = pSNPSite->SNPLoci + 500;
	switch (m_ReportFormat) {
		case eRMFpgSNP:						// report in pgSNP format
			// a special case: UCSC references the reference SARS-CoV-2 as being NC_O45512v2,not NC_O45512.2
			// so change NC_O45512.2 to be NC_O45512v2 when reporting as pgSNP format
			if(!stricmp(pSNPSite->szChrom, "NC_045512.2"))
				strcpy(szChrom, "NC_045512v2");
			else
				strcpy(szChrom, pSNPSite->szChrom);
			m_OutBuffOffs = sprintf(m_pszOutBuff, "track type=pgSnp visibility=3 db=%s name=\"%s\" description=\"%s\"\n", m_SpecAssemblyName, m_szTrackName, m_szDescription);
			m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], "browser position %s:%d-%d", szChrom, StartLoci, EndLoci);
			break;

		case eRMFvcf:						// report in VCF 4.1 format
			m_OutBuffOffs = sprintf(m_pszOutBuff, "##fileformat=VCFv4.1\n##source=ngskit4b%s\n##reference=%s\n##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n",
										kit4bversion, m_szInSNPsFile);
			m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
			break;

		case eRMFcsv:						// report in CSV format
			m_OutBuffOffs = sprintf(m_pszOutBuff, "\"Readset\",\"Site ID\",\"Site Chrom\",\"Site Loci\",\"SNP ID\",\"Ref Base\",\"CntA\",\"CntC\",\"CntG\",\"CntT\",\"CntN\"");
			m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], ",\"Alleles\",\"Depth\"");

			m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], ",\"Feature Chrom\",\"Feature Name\",\"Feature Start Loci\",\"Feature End Loci\"");
			m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], ",\"Synonymous\",\"NonSynonymous\",\"Wobble\"");
			m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], ",\"NumPosFreqBiasedCodons\",\"NumNegFreqBiasedCodons\"");
			m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], ",\"IsDirac\",\"Ref Peptide\"");
			for(int ToPeptideIdx = 0; ToPeptideIdx <= 20; ToPeptideIdx++)
				m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], ",\"To:%s\"", gAminoAcids[ToPeptideIdx].szPeptide);
			m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], "\n");
			break;						
		}
	if((Rslt = WriteOutFile(true)) < eBSFSuccess)
		return(Rslt);
	}
else
	{
	if(m_hOutSiteDist != -1 && m_pszSiteDistBuff != NULL)
		{
		if((m_SiteDistOffs + 500) >= m_AllocdSiteDistBuff)
			{
			CUtility::RetryWrites(m_hOutSiteDist, m_pszSiteDistBuff, m_SiteDistOffs);
			m_SiteDistOffs = 0;
			}

		if(pSNPSite->bSNPPlaceholder)
			{
			m_SiteDistOffs += sprintf(&m_pszSiteDistBuff[m_SiteDistOffs], "%u,\"%s\",%u,\"N\",", pSNPSite->SNPId, pSNPSite->szChrom, pSNPSite->SNPLoci);
			m_SiteDistOffs += sprintf(&m_pszSiteDistBuff[m_SiteDistOffs], "0,0,");
			m_SiteDistOffs += sprintf(&m_pszSiteDistBuff[m_SiteDistOffs], "0,0,0,0,0,");
			m_SiteDistOffs += sprintf(&m_pszSiteDistBuff[m_SiteDistOffs], "0,0,0,0,0\n");
			return(eBSFSuccess);
			}

		switch(pSNPSite->RefBase)
			{
				case 0: Base = 'A'; break;
				case 1: Base = 'C'; break;
				case 2: Base = 'G'; break;
				case 3: Base = 'T'; break;
				case 4: Base = 'N'; break;
			}
		m_SiteDistOffs += sprintf(&m_pszSiteDistBuff[m_SiteDistOffs],"%u,\"%s\",%u,\"%c\",", pSNPSite->SNPId, pSNPSite->szChrom, pSNPSite->SNPLoci,Base);
		m_SiteDistOffs += sprintf(&m_pszSiteDistBuff[m_SiteDistOffs], "%u,%u,", pSNPSite->TotBaseCnts, pSNPSite->TotMismatches);
		m_SiteDistOffs += sprintf(&m_pszSiteDistBuff[m_SiteDistOffs], "%u,%u,%u,%u,%u,", pSNPSite->BaseCnts[0].Cnts, pSNPSite->BaseCnts[1].Cnts, pSNPSite->BaseCnts[2].Cnts, pSNPSite->BaseCnts[3].Cnts, pSNPSite->BaseCnts[4].Cnts);
		m_SiteDistOffs += sprintf(&m_pszSiteDistBuff[m_SiteDistOffs],"%u,%u,%u,%u,%u\n",pSNPSite->DiracCnts[0], pSNPSite->DiracCnts[1], pSNPSite->DiracCnts[2], pSNPSite->DiracCnts[3], pSNPSite->DiracCnts[4]);
		}

	if((Rslt = WriteOutFile()) < eBSFSuccess)
		return(Rslt);
	switch (m_ReportFormat) {
		case eRMFpgSNP:						// report in pgSNP format
			// a special case: UCSC references the reference SARS-CoV-2 as being NC_O45512v2,not NC_O45512.2
			// so change NC_O45512.2 to be NC_O45512v2 when reporting as pgSNP format

			if(!stricmp(pSNPSite->szChrom, "NC_045512.2"))
				strcpy(szChrom, "NC_045512v2");
			else
				strcpy(szChrom, pSNPSite->szChrom);
			m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], "\n%s\t%d\t%d\t", szChrom, pSNPSite->SNPLoci, pSNPSite->SNPLoci + 1);
			NumAlleles = 0;
			for (int Idx = 0; Idx < 4; Idx++)
				{
				if (pSNPSite->BaseCnts[Idx].Cnts > 0)
					{
					if (m_pszOutBuff[m_OutBuffOffs - 1] != '\t')
						m_pszOutBuff[m_OutBuffOffs++] = '/';
					switch (Idx) {
						case 0: Base = 'A'; break;
						case 1: Base = 'C'; break;
						case 2: Base = 'G'; break;
						case 3: Base = 'T'; break;
						case 4: Base = 'N'; break;
						}
					NumAlleles++;
					m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], "%c", Base);
					}
				}
			m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], "\t%d\t", NumAlleles);
			for (int Idx = 0; Idx < 4; Idx++)
				{
				if (pSNPSite->BaseCnts[Idx].Cnts > 0)
					{
					if (m_pszOutBuff[m_OutBuffOffs - 1] != '\t')
						m_pszOutBuff[m_OutBuffOffs++] = ',';
					m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], "%d", pSNPSite->BaseCnts[Idx].Cnts);
					}
				}
			m_pszOutBuff[m_OutBuffOffs++] = '\t';

			// following is gross misuse of scores! But I want to let the viewer know the percentages of all isolates/cultivars/whatever at the SNP loci sharing at least one allelic variant
			// so first the proportion ( * 10)  of all sites is reported
			PercentScore = max(1, (int)(0.5 + ((double)(pSNPSite->TotMismatches * 100.0) / pSNPSite->TotBaseCnts)));
			m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], "\t%d", PercentScore);
			// next the population size
			if(NumAlleles >= 2)
				m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], ",%d", pSNPSite->TotBaseCnts);
			// next the number of samples which have at least one allelic variation
			if(NumAlleles >= 3)
				m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], ",%d", pSNPSite->TotMismatches);
			if(NumAlleles == 4)
				m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], ",0");
			break;

		case eRMFvcf:						// report in VCF 4.1 format
			NumAlleles = 0;
			SumQScores = 0;
			Depth = 0;
			if(!stricmp(pSNPSite->szChrom, "NC_045512.2"))
				strcpy(szChrom, "NC_045512v2");
			else
				strcpy(szChrom, pSNPSite->szChrom);
			m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], "%s\t%u\tSNP%d\t%c\t", szChrom, pSNPSite->SNPLoci + 1, pSNPSite->SNPId, CSeqTrans::MapBase2Ascii(pSNPSite->RefBase));
			
			for (int Idx = 0; Idx < 4; Idx++)
				{
				Depth += pSNPSite->BaseCnts[Idx].Cnts;
				if(Idx == pSNPSite->RefBase)
					continue;
				if (pSNPSite->BaseCnts[Idx].Cnts > 0)
					{
					if (m_pszOutBuff[m_OutBuffOffs - 1] != '\t')
						m_pszOutBuff[m_OutBuffOffs++] = ',';
					switch (Idx) {
						case 0: Base = 'A'; break;
						case 1: Base = 'C'; break;
						case 2: Base = 'G'; break;
						case 3: Base = 'T'; break;
						case 4: Base = 'N'; break;
						}
					NumAlleles++;
					m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], "%c", Base);
					SumQScores += min(99.0, ((m_PValueThres - pSNPSite->BaseCnts[Idx].QScore) * (100 / m_PValueThres)));
					}
				}
			SumQScores /= NumAlleles;
			m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], "\t%d\tPASS\tAF=", (int)SumQScores);
			for (int Idx = 0; Idx < 4; Idx++)
				{
				if (Idx == pSNPSite->RefBase)
					continue;
				if (pSNPSite->BaseCnts[Idx].Cnts > 0)
					{
					if (m_pszOutBuff[m_OutBuffOffs - 1] != '=')
						m_pszOutBuff[m_OutBuffOffs++] = ',';
					m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], "%.4f", pSNPSite->BaseCnts[Idx].Cnts/(double)Depth);
					}
				}
			m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs],";DP=%d\n", Depth);
			break;

		case eRMFcsv:						// report in CSV format
			NumAlleles = 0;
			SumQScores = 0;
			Depth = 0;
			m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], "\"%s\",%u,\"%s\",%u,%d,%c", pSNPSite->szSeqReadSet,pSNPSite->ReadsetSiteId, pSNPSite->szChrom, pSNPSite->SNPLoci, pSNPSite->SNPId, CSeqTrans::MapBase2Ascii(pSNPSite->RefBase));
			
			for(int Idx = 0; Idx <= 4; Idx++)
				m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], ",%u", pSNPSite->OrigBaseCnts[Idx]);
			m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], ",");

			for(int Idx = 0; Idx < 4; Idx++)
				{
				Depth += pSNPSite->BaseCnts[Idx].Cnts;
				if(Idx == pSNPSite->RefBase)
					continue;
				if(pSNPSite->BaseCnts[Idx].Cnts > 0)
					{
					if(m_pszOutBuff[m_OutBuffOffs - 1] != ',')
						m_pszOutBuff[m_OutBuffOffs++] = '|';
					switch(Idx)
						{
							case 0: Base = 'A'; break;
							case 1: Base = 'C'; break;
							case 2: Base = 'G'; break;
							case 3: Base = 'T'; break;
							case 4: Base = 'N'; break;
						}
					NumAlleles++;
					m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], "%c:%d:%f", Base, pSNPSite->BaseCnts[Idx].Cnts, pSNPSite->BaseCnts[Idx].QScore);
					}
				}
			m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], ",%d", Depth);
			if(m_pFeatures != NULL && m_NumFeatures)
				{
				tsSummaryFeatCnts *pFeature;
				if(pSNPSite->bInFeature)
					{
					pFeature = &m_pFeatures[pSNPSite->FeatureIdx];
					m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], ",\"%s\",\"%s\",%u,%u", pFeature->szChrom, pFeature->szFeatName, pFeature->FeatStart, pFeature->FeatEnd);
					m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], ",%u,%u,%u", pFeature->SiteFeatCnts.TotSynonymous, pFeature->SiteFeatCnts.TotNonSynonymous, pFeature->SiteFeatCnts.TotWobble);
					m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], ",%u,%u", pFeature->SiteFeatCnts.NumPosBiasedCodons, pFeature->SiteFeatCnts.NumNegBiasedCodons);
					m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], ",%u", pFeature->SiteFeatCnts.IsDirac ? 1 : 0);
					m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs],",\"%s\"", gAminoAcids[pFeature->SiteFeatCnts.RefSynGroup].szPeptide);
					for(int ToPeptideIdx = 0; ToPeptideIdx <= 20; ToPeptideIdx++)
						m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], ",%u", pFeature->SiteFeatCnts.ToAminoAcidChanges[ToPeptideIdx]);
					}
				else
					{
					m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], ",\"-\",\"-\",0,0");
					m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], ",0,0,0,0,0,0");
					m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], "0,\"-\"");
					for(int ToPeptideIdx = 0; ToPeptideIdx <= 20; ToPeptideIdx++)
						m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], ",0");
					}
				}

			m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], "\n");
			if((Rslt = WriteOutFile()) < eBSFSuccess)
				return(Rslt);
			break;
		}
	}
return(eBSFSuccess);
}

int
CSNPs2pgSNPs::ReportFeatSNPcnts(bool bHeader,	// header required
						int NumFeatures,		// number of features
						tsSummaryFeatCnts* pFeatures)	// array of feature SNP counts and codons
{
int Rslt;
int FeatIdx;
int FromPeptideIdx;
int ToPeptideIdx;

if(NumFeatures == 0 || pFeatures == NULL)	// may be no features to be reported on!
	return(eBSFSuccess);

m_OutBuffOffs = 0;
if(bHeader)
	{
	m_OutBuffOffs = sprintf(m_pszOutBuff, "\"Chrom\",\"Feature\",\"Start Loci\",\"End Loci\"");
	m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], ",\"SNPS\",\"Synonymous\",\"NonSynonymous\",\"Wobble\",\"Excl.Synonymous\",\"Excl.NonSynonymous\",\"Excl.Wobble\"");
	m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], ",\"SiteSynonymous\",\"NumPosFreqBiasedCodons\",\"NumNegFreqBiasedCodons\"");
	m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], ",\"DiracA\",\"DiracC\",\"DiracG\",\"DiracT\"");
	for(FromPeptideIdx = 0; FromPeptideIdx <= 20; FromPeptideIdx++)
		{
		m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs],",\"From:%s\"", gAminoAcids[FromPeptideIdx].szPeptide);
		for(ToPeptideIdx = 0; ToPeptideIdx <= 20; ToPeptideIdx++)
			m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], ",\"To:%s\"", gAminoAcids[ToPeptideIdx].szPeptide);
		if((Rslt = WriteOutFile())<eBSFSuccess)
			return(Rslt);
		}
	m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], "\n");
	}

tsSummaryFeatCnts* pFeature = pFeatures;

for(FeatIdx = 0; FeatIdx < NumFeatures; FeatIdx++, pFeature++)
	{
	if((Rslt = WriteOutFile()) < eBSFSuccess)
		return(Rslt);
	m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs],"\"%s\",\"%s\",%u,%u",pFeature->szChrom, pFeature->szFeatName, pFeature->FeatStart, pFeature->FeatEnd);
	m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs],",%u,%u,%u,%u", pFeature->TotSNPs,pFeature->TotSynonymous,pFeature->TotNonSynonymous,pFeature->TotWobble);
	m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], ",%u,%u,%u", pFeature->TotExclSynonymous, pFeature->TotExclNonSynonymous, pFeature->TotExclWobble);
	m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], ",%u,%u,%u", pFeature->NumSiteSynonCodons,pFeature->NumPosBiasedCodons, pFeature->NumNegBiasedCodons);
	m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], ",%u,%u,%u,%u", pFeature->DiracCnts[0], pFeature->DiracCnts[1], pFeature->DiracCnts[2], pFeature->DiracCnts[3]);
	for(FromPeptideIdx = 0; FromPeptideIdx <= 20; FromPeptideIdx++)
		{
		if((Rslt = WriteOutFile()) < eBSFSuccess)
			return(Rslt);
		m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], ",%u", pFeature->FromAminoAcidChanges[FromPeptideIdx]);
		for(ToPeptideIdx = 0; ToPeptideIdx <= 20; ToPeptideIdx++)
			m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs],",%u",pFeature->ToAminoAcidChanges[FromPeptideIdx][ToPeptideIdx]);
		}
	m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], "\n");
	if((Rslt = WriteOutFile(true)) < eBSFSuccess)
		return(Rslt);
	}

return(eBSFSuccess);
}


int
CSNPs2pgSNPs::WriteOutFile(bool bForce)	// true if write is forced
{
if(m_hOutFile == -1)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "WriteOutFile: no opened file");
	Reset();
	return(eBSFerrClosed);
	}
if(bForce || ((m_OutBuffOffs + 200000) >= m_AllocdOutBuff))
	{
	if(m_OutBuffOffs)
		{
		if(CUtility::RetryWrites(m_hOutFile, m_pszOutBuff, m_OutBuffOffs) == false)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "WriteOutFile: write to '%s' failed", m_szOutFile);
			close(m_hOutFile);
			m_hOutFile = -1;
			Reset();
			return(eBSFerrFileAccess);
			}
		m_OutBuffOffs = 0;
		}
	}
return(eBSFSuccess);
}

int
CSNPs2pgSNPs::CreateOutFile(char *pszExtn)
{
// close any previously opened output file
CloseOutFile();

// add extension to szOutFilePrefix to form szOutFile
strcpy(m_szOutFile,m_szOutFilePrefix);
strcat(m_szOutFile, pszExtn);
#ifdef _WIN32
if((m_hOutFile = open(m_szOutFile, (_O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC), (_S_IREAD | _S_IWRITE))) == -1)
#else
if((m_hOutFile = open(m_szOutFile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE)) == -1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to create/truncate output file - '%s' - %s", m_szOutFile, strerror(errno));
	Reset();
	return(eBSFerrOpnFile);
	}
return(eBSFSuccess);
}

void
CSNPs2pgSNPs::CloseOutFile(void)			// closes output file if currently opened
{
if(m_hOutFile != -1)
	{
	if(m_OutBuffOffs && m_pszOutBuff != NULL)
		CUtility::RetryWrites(m_hOutFile, m_pszOutBuff, m_OutBuffOffs);

#ifdef _WIN32
	_commit(m_hOutFile);
#else
	fsync(m_hOutFile);
#endif
	close(m_hOutFile);
	m_hOutFile = -1;
	}
m_OutBuffOffs = 0;
}

int 
CSNPs2pgSNPs::Process(eModepgSSNP Mode,				// processing mode
					bool bAllowInferenced,		// true if both original kalign SNP calls and snpmarker inferenced SNP calls to be used, default is for kalign SNPs only
					double LocalSeqErrRate,			// local sequencing error rate in a 100bp window centered around the putative SNP site
					int MinCoverage,				// must be at least this coverage at SNP site
					double MinAlleleProp,			// putative allele must be at least this proportion of total site read coverage
					double PValueThres,				// only accept SNP alleles which have a PValue <= this threshold
					char *pszTrackName,				// track name
					char *pszAssemblyName,			// UCSC assembly name - for SARS-CoV-2 it is "wuhCor1"
					char* pszExperimentDescr,		// describes experiment
					etSSetOp SetOp,					// set operation on SetA and/or SetB				
					int NumInSetA,					// number of species/cultivars/isolates in SetA
					char** ppszSetA,				// names of those cspecies/cultivars/isolates which are in SetA
					int NumInSetB,					// number of species/cultivars/isolates in SetB
					char** ppszSetB,				// names of those cspecies/cultivars/isolates in SetB
					char* pszGFFFile,				// general feature format file - identifies start/end loci of features in assembly, i.e. transcripts etc.
					char* pszSNPFile,				// load SNP calls from this CSV file, can be either SNPs generated by kalign or snpmarkers
					char* pszOutFile)				// output SNPs to this UCSC Personal Genome SNP format file
{
int Rslt;
int Len;
Reset();
strcpy(m_szInSNPsFile, pszSNPFile);
strcpy(m_SpecAssemblyName, pszAssemblyName);
strcpy(m_szDescription,pszExperimentDescr);
strcpy(m_szTrackName, pszTrackName);
if(pszGFFFile != NULL && pszGFFFile[0] != '\0')
	strcpy(m_szGFFFile, pszGFFFile);
else
	m_szGFFFile[0] = '\0';
m_bAllowInferenced = bAllowInferenced;
m_LocalSeqErrRate = LocalSeqErrRate;
m_PValueThres = PValueThres;
m_MinCoverage = MinCoverage;
m_MinAlleleProp = MinAlleleProp;
m_SetOp = SetOp;
m_NumInSetA = NumInSetA;
m_ppszSetA = ppszSetA;
m_NumInSetB = NumInSetB;
m_ppszSetB = ppszSetB;
m_DiracThres = cDfltDiracThres;

if(m_szGFFFile[0] != '\0')
	{
	if((m_pGFFFile = new CGFFFile) == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to instantiate CGFFFile");
		Reset();
		return(eBSFerrObj);
		}
	if((Rslt = m_pGFFFile->Open(pszGFFFile)) != eBSFSuccess)
		{
		while(m_pGFFFile->NumErrMsgs())
			gDiagnostics.DiagOut(eDLFatal, gszProcName, m_pGFFFile->GetErrMsg());
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to open input GFF3 file '%s'", pszGFFFile);
		Reset();
		return(Rslt);
		}

	if((m_pFeatures = new tsSummaryFeatCnts[cAllocGeneFeatures])==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to allocate memory for gene feature counts");
		Reset();
		return(eBSFerrMem);
		}
	m_AllocFeatures = cAllocGeneFeatures;
	m_NumFeatures = 0;
	tsSummaryFeatCnts *pGeneCnts = m_pFeatures;
	char *pFeature;
	int StartAt;
	int EndAt;
	char * pszTmp;
	while(m_NumFeatures < m_AllocFeatures && (Rslt = m_pGFFFile->NextRecordOfType(eGGFFany)) > 0)
		{
		pFeature = m_pGFFFile->GetFeature();
		if(stricmp(pFeature, (char*)"CDS"))
			continue;
		if(m_pGFFFile->GetStrand() != '+')
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Sense strand only CDS features supported, pszGFFFile '%s'", pszGFFFile);
			Reset();
			return(eBSFerrFeature);
			}

		memset(pGeneCnts,0,sizeof(tsSummaryFeatCnts));
		pszTmp = m_pGFFFile->GetNoteValue((char *)"gene_id");
		strncpy(pGeneCnts->szFeatName, pszTmp, sizeof(pGeneCnts->szFeatName));
		CUtility::TrimQuotes(pGeneCnts->szFeatName);
		pszTmp = m_pGFFFile->GetSeqName();
		strncpy(pGeneCnts->szChrom, pszTmp, sizeof(pGeneCnts->szChrom));
		StartAt = m_pGFFFile->GetStart();
		EndAt = m_pGFFFile->GetEnd() + 3;	// CDS size does not include stop codon
		pGeneCnts->FeatStart = StartAt;
		pGeneCnts->FeatEnd = EndAt;
		pGeneCnts++;
		m_NumFeatures++;
		}
	m_pGFFFile->Close();
	delete m_pGFFFile;
	m_pGFFFile = NULL;
	if(m_NumFeatures > m_AllocFeatures)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Too many CDS gene features in pszGFFFile '%s', hard limit of %d", pszGFFFile,cAllocGeneFeatures);
		Reset();
		return(eBSFerrFeature);
		}
	}


if ((m_pszOutBuff = new char[cAllocLineBuffSize]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to allocate memory for output line buffering -- %s", strerror(errno));
	return(eBSFerrMem);
	}
m_AllocdOutBuff = cAllocLineBuffSize;

if((m_pGapCntDist = new uint32_t[cMaxGapBetweenSNPS+1]) == NULL)
{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to allocate memory for inter SNP gap distributions -- %s", strerror(errno));
	Reset();
	return(eBSFerrMem);
}
memset(m_pGapCntDist, 0, sizeof(uint32_t) * (cMaxGapBetweenSNPS + 1));


if((m_pSharedSNPsDist = new uint32_t[cMaxSharedSNPS+1]) == NULL)
{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to allocate memory for shared SNP distributions -- %s", strerror(errno));
	Reset();
	return(eBSFerrMem);
}
memset(m_pSharedSNPsDist, 0, sizeof(uint32_t) * (cMaxSharedSNPS + 1));


if((m_pIndividualSNPsDist = new uint32_t[cMaxIndividualSNPS+1]) == NULL)
{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to allocate memory for individual SNP distributions -- %s", strerror(errno));
	Reset();
	return(eBSFerrMem);
}
memset(m_pIndividualSNPsDist, 0, sizeof(uint32_t) * (cMaxIndividualSNPS + 1));


// check file extension, if '.vcf' then generate output formated for vcf, if '.csv' then report format for csv, instead of the default pgSNP format
strcpy(m_szOutFilePrefix, pszOutFile);
m_ReportFormat = eRMFsnp::eRMFpgSNP;	// default output format
Len = (int)strlen(pszOutFile);
if(Len >= 8)
	{
	if(!stricmp(&pszOutFile[Len - 6], ".pgSNP"))
		m_szOutFilePrefix[Len - 6] = '\0';
	}
if(Len >= 5)
	{
	if(!stricmp(&pszOutFile[Len-4],".vcf"))
		{
		m_ReportFormat = eRMFsnp::eRMFvcf;
		m_szOutFilePrefix[Len-4] = '\0';
		}
	else
		if(!stricmp(&pszOutFile[Len - 4], ".csv"))
			{
			m_ReportFormat = eRMFsnp::eRMFcsv;
			m_szOutFilePrefix[Len - 4] = '\0';
			}
	}

// create/open the primary output file ready for writing
switch(m_ReportFormat) {
	case eRMFcsv:
		Rslt = CreateOutFile((char *)".csv");
		break;
	case eRMFvcf:
		Rslt = CreateOutFile((char *)".vcf");
		break;
	default:
		Rslt = CreateOutFile((char *)".pgSNP");
		break;
	}
if(Rslt < eBSFSuccess)
	{
	Reset();
	return(Rslt);
	}

if ((m_pCSV = new CCSVFile) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to instantiate CCSVfile");
	Reset();
	return(eBSFerrObj);
	}
if(Mode == eMpgSSNPmarkers)
	m_pCSV->SetMaxFields((cMaxCultivars *9) + 4);	// number of fields is unknown so need to allow for a worst case
else
	m_pCSV->SetMaxFields(cAlignNumSSNPXfields);


Rslt = eBSFSuccess;
m_ReadsetID = 0;
m_bRptHdr = true;
char *pszCurSNPFile;
CSimpleGlob glob(SG_GLOB_FULLSORT);
switch(Mode) {
	case eMpgSSNPKalign:
		glob.Init();
		if(glob.Add(pszSNPFile) < SG_SUCCESS)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to glob '%s", pszSNPFile);
			Reset();
			return(eBSFerrOpnFile);	// treat as though unable to open file
			}

		if(glob.FileCount() <= 0)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to locate any SNP file matching '%s", pszSNPFile);
			Reset();
			return(eBSFerrOpnFile);	// treat as though unable to open file
			}

		Rslt = eBSFSuccess;
		for(int FileID = 0; Rslt >= eBSFSuccess && FileID < glob.FileCount(); ++FileID)
			{
			pszCurSNPFile = glob.File(FileID);

			gDiagnostics.DiagOut(eDLInfo, gszProcName, "Processing SNP file '%s'\n", pszCurSNPFile);
			strcpy(m_szInSNPsFile, pszCurSNPFile);
			if((Rslt = m_pCSV->Open(pszCurSNPFile)) != eBSFSuccess)
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to open file: '%s'", pszCurSNPFile);
				Reset();
				return(Rslt);
				}

			// readset name is simply the input SNP file sans path and file name extension if '.snp.csv' or just '.csv' - very simplistic but SNP CSVs do not contain sequencing readset names!
			CUtility::splitpath(pszCurSNPFile,NULL, m_ReadsetIdentifier);
			if((Len = (int)strlen(m_ReadsetIdentifier)) >= 10 && !stricmp(&m_ReadsetIdentifier[Len - 9], ".snps.csv"))
				m_ReadsetIdentifier[Len - 9] = '\0';
			else
				if((Len = (int)strlen(m_ReadsetIdentifier)) >= 5 && !stricmp(&m_ReadsetIdentifier[Len - 4], ".csv"))
					m_ReadsetIdentifier[Len - 3] = '\0';
			if((Rslt = m_ReadsetID = AddReadset(m_ReadsetIdentifier)) < eBSFSuccess)
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to add readset name: '%s'", m_ReadsetIdentifier);
				Reset();
				return(Rslt);
				}
			Rslt = ProcessKalignSNPs();
			if(m_pCSV != NULL)
				m_pCSV->Close();
			if(Rslt < eBSFSuccess)
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Failed processing SNP file '%s'\n", pszCurSNPFile);
				Reset();
				return(Rslt);
				}
			if(Rslt == 0)			// 0 only if file was parsed successfully but contained no SNPs
				{
				gDiagnostics.DiagOut(eDLWarn, gszProcName, "Processed SNP file '%s' contained no accepted SNPs", pszCurSNPFile);
				m_NumReadsetNames -= 1;
				}
			}
		CloseOutFile();					// primary output file
		Rslt = 0;
		if(m_NumReadsetNames == 0)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "No SNP file contained any accepted SNPs");
			Reset();
			return(eBSFerrLocField);
			}
		m_NumCultivars = m_NumReadsetNames;
		// load and report loci isolate matrix
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Initialising loci isolate matrix ...");
		if((Rslt = LoadLociIsolateMatrix()) < 0)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Initialising loci isolate matrix failed");
			Reset();
			return(Rslt);
			}
		
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Reporting loci isolate matrix ...");
		if((Rslt = ReportLociIsolateMatrix((char *)".ReadsetLociMatrix.csv"))<0)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Reporting loci isolate matrix failed");
			Reset();
			return(Rslt);
			}

		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Reporting readset similarities matrix ...");
		if((Rslt = ReportSimilaritiesMatrix((char *)".ReadsetSimilaritiesMatrix.csv")) < 0)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Reporting readset similarities matrix failed");
			Reset();
			return(Rslt);
			}

		// now identifying the shared SNP site loci distributions ..
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Identifying shared SNP site loci distributions ...");
		if((Rslt = IdentifySharedSiteLociDist((char *)".SharedReadsetLociDist.csv")) < 0)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Identifying shared SNP site loci distributions failed");
			Reset();
			return(Rslt);
			}


		if(m_NumFeatures)
			{
			if((Rslt = CreateOutFile((char *)".Feats.csv")) < eBSFSuccess)
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Creating output file failed");
				Reset();
				return(Rslt);
				}
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "Reporting shared SNP site loci distributions ...");
			Rslt = ReportFeatSNPcnts(true,m_NumFeatures,m_pFeatures);
			CloseOutFile();
			if(Rslt < 0)
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Reporting shared SNP site loci distributions failed");
				Reset();
				return(Rslt);
				}
			}
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Reporting distributions completed");
		break;

	case eMpgSSNPmarkers:
		if((Rslt = m_pCSV->Open(pszSNPFile)) != eBSFSuccess)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to open file: %s", pszSNPFile);
			Reset();
			return(Rslt);
			}
			// readset name is simply the input SNP file sans path and file name extension if it is '.csv' - very simplistic but snpmarkers CSVs do not contain sequencing readset names!
		CUtility::splitpath(pszSNPFile, NULL, m_ReadsetIdentifier);
		if((Len = (int)strlen(m_ReadsetIdentifier)) >= 5 && !stricmp(&m_ReadsetIdentifier[Len - 4], ".csv"))
			m_ReadsetIdentifier[Len-5] = '\0';
		
		if((m_pszSiteDistBuff = new char[cAllocLineBuffSize]) == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to allocate memory for site distribution line buffering -- %s", strerror(errno));
			Reset();
			return(eBSFerrMem);
			}
		m_AllocdSiteDistBuff = cAllocLineBuffSize;

		strcpy(m_szOutSiteDistFile, pszOutFile);
		strcat(m_szOutSiteDistFile, ".sitedist.csv");
#ifdef _WIN32
		if((m_hOutSiteDist = open(m_szOutSiteDistFile, (_O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC), (_S_IREAD | _S_IWRITE))) == -1)
#else
		if((m_hOutSiteDist = open(m_szOutSiteDistFile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE)) == -1)
#endif
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to create/truncate output file - '%s' - %s", m_szOutSiteDistFile, strerror(errno));
			Reset();
			return(eBSFerrOpnFile);
			}


		// header is generated here rather than in reporting function as actual SNP reporting may occur after place holder SNP loci have been reported ....
		m_SiteDistOffs = sprintf(m_pszSiteDistBuff, "\"ID\",\"Chrom\",\"Loci\",\"RefBase\",\"CovBases\",\"AlleleBases\",");
		m_SiteDistOffs += sprintf(&m_pszSiteDistBuff[m_SiteDistOffs], "\"AlleleA\",\"AlleleC\",\"AlleleG\",\"AlleleT\",\"AlleleN\",");
		m_SiteDistOffs += sprintf(&m_pszSiteDistBuff[m_SiteDistOffs], "\"DiracA\",\"DiracC\",\"DiracG\",\"DiracT\",\"DiracN\"\n");

		Rslt = ProcessSnpmarkersSNPs();

		if(Rslt >= eBSFSuccess)
			{
			if((Rslt=CreateOutFile((char *)".GapDist.csv"))< eBSFSuccess)
				{
				Reset();
				return(Rslt);
				}

			m_OutBuffOffs = sprintf(m_pszOutBuff, "\"Len\",\"Instances\"\n");
			for(uint32_t Idx = 0; Idx <= m_MaxSNPgap; Idx++)
				m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], "%d,%d\n", Idx + 1, m_pGapCntDist[Idx]);

			CloseOutFile();


			if((Rslt = CreateOutFile((char *)".SharedReadsetLociDist.csv")) < eBSFSuccess)
				{
				Reset();
				return(Rslt);
				}

			m_OutBuffOffs = sprintf(m_pszOutBuff, "\"Shared\",\"Instances\"\n");
			for(uint32_t Idx = 0; Idx <= m_MaxSharedSNPs; Idx++)
				m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], "%d,%d\n", Idx + 1, m_pSharedSNPsDist[Idx]);
			CloseOutFile();

			if((Rslt = CreateOutFile((char *)".IndividualSnpDist.csv")) < eBSFSuccess)
				{
				Reset();
				return(Rslt);
				}

			m_OutBuffOffs = sprintf(m_pszOutBuff, "\"SNPs\",\"Instances\"\n");
			for(uint32_t Idx = 0; Idx <= m_MaxIndividualSNPs; Idx++)
				m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], "%d,%d\n", Idx, m_pIndividualSNPsDist[Idx]);
			CloseOutFile();
		}
	break;
	}

Reset();
return(Rslt);
}

uint32_t					// returned number of unique sites
CSNPs2pgSNPs::GenNumUniqueSiteLoci(void)
{
uint32_t NumUniqueSiteLoci;
uint32_t CurSNPLoci;
uint32_t CurSNPChromID;
uint32_t FeatSNPIdx;
tsIsolateFeatSNPs* pFeatSNP;

if(m_IsolateFeatSNPs < 2)
	return(m_IsolateFeatSNPs);

// sort chrom.loci.readset ascending
if(m_SortOrder != 1)	// 0: unsorted, 1: SortSNPChromLociReadset, 2: SortSNPReadsetChromLoci
	{
	qsort(m_pIsolateFeatSNPs, m_IsolateFeatSNPs, sizeof(tsIsolateFeatSNPs), SortSNPChromLociReadset);
	m_SortOrder = 1;
	}

// determine number of unique site loci in target genome
NumUniqueSiteLoci = 0;
CurSNPLoci = (uint32_t)-1;
CurSNPChromID = (uint32_t)-1;
pFeatSNP = m_pIsolateFeatSNPs;
for(FeatSNPIdx = 0; FeatSNPIdx < m_IsolateFeatSNPs; FeatSNPIdx++, pFeatSNP++)
	{
	if(pFeatSNP->SNPSSite.ChromID == CurSNPChromID && pFeatSNP->SNPSSite.SNPLoci == CurSNPLoci)
		continue;
	CurSNPLoci = pFeatSNP->SNPSSite.SNPLoci;
	CurSNPChromID = pFeatSNP->SNPSSite.ChromID;
	NumUniqueSiteLoci++;
	}
return(NumUniqueSiteLoci);
}

int
CSNPs2pgSNPs::LoadLociIsolateMatrix(void)				// generate a matrix (rows - isolates, columns - loci) with cells absence/presence of sites at that loci in that isolate
{
tsIsolateFeatSNPs* pFeatSNP;
uint32_t * pCurMatrixCell;
uint32_t *pCurLoci;
uint32_t FeatSNPIdx;
uint32_t ReadsetIdx;
uint32_t CurReadsetID;
uint32_t CurSNPLoci;
uint32_t LociIdx;

if(m_pMatrix != NULL)
	{
	delete []m_pMatrix;
	m_pMatrix = NULL;
	}
if(m_NumUniqueSiteLoci == 0)
	m_NumUniqueSiteLoci = GenNumUniqueSiteLoci();

m_pMatrix = new uint32_t[(size_t)(m_NumCultivars + 1) * (m_NumUniqueSiteLoci+1)];
memset(m_pMatrix,0, (size_t)(m_NumCultivars + 1) * (m_NumUniqueSiteLoci + 1) * sizeof(uint32_t));
m_MatrixRows = m_NumCultivars + 1;
m_MatrixCols = m_NumUniqueSiteLoci + 1;
// sort chrom.loci.readset ascending so can fill in unique loci as first row
if(m_SortOrder != 1)	// 0: unsorted, 1: SortSNPChromLociReadset, 2: SortSNPReadsetChromLoci
	{
	qsort(m_pIsolateFeatSNPs, m_IsolateFeatSNPs, sizeof(tsIsolateFeatSNPs), SortSNPChromLociReadset);
	m_SortOrder = 1;
	}

pFeatSNP = m_pIsolateFeatSNPs;
CurSNPLoci = (uint32_t)-1;
LociIdx = 1;
for(FeatSNPIdx = 0; FeatSNPIdx < m_IsolateFeatSNPs; FeatSNPIdx++, pFeatSNP++)
	{
	if(CurSNPLoci != pFeatSNP->SNPSSite.SNPLoci)
		{
		m_pMatrix[LociIdx++] = pFeatSNP->SNPSSite.SNPLoci;
		CurSNPLoci = pFeatSNP->SNPSSite.SNPLoci;
		}
	}

// sort readset.chrom.loci ascending so can process by readset 
if(m_SortOrder != 2)	// 0: unsorted, 1: SortSNPChromLociReadset, 2: SortSNPReadsetChromLoci
	{
	qsort(m_pIsolateFeatSNPs, m_IsolateFeatSNPs, sizeof(tsIsolateFeatSNPs), SortSNPReadsetChromLoci);
	m_SortOrder = 2;
	}
pFeatSNP = m_pIsolateFeatSNPs;
ReadsetIdx = 0;
CurReadsetID = (uint32_t)-1;
for(FeatSNPIdx = 0; FeatSNPIdx < m_IsolateFeatSNPs; FeatSNPIdx++, pFeatSNP++)
	{
	if(CurReadsetID != pFeatSNP->SNPSSite.ReadsetID)
		{
		ReadsetIdx += m_NumUniqueSiteLoci + 1;
		CurReadsetID = pFeatSNP->SNPSSite.ReadsetID;
		m_pMatrix[ReadsetIdx] = CurReadsetID;
		pCurLoci = &m_pMatrix[1];
		pCurMatrixCell = &m_pMatrix[ReadsetIdx + 1];
		}

	while(pFeatSNP->SNPSSite.SNPLoci != *pCurLoci)
		{
		pCurMatrixCell++;
		pCurLoci++;
		}
	*pCurMatrixCell = (uint32_t)pFeatSNP->SNPSSite.ClassifySite;
	}
return(eBSFSuccess);
}


int
CSNPs2pgSNPs::ReportLociIsolateMatrix(char *pszFileExtn)	// report loci vs isolate matrix using this file extension, will generate loci isolate matrix if not already generated
{
int Rslt;
uint32_t LociIdx;
uint32_t ReadsetIdx;

if(m_pMatrix == NULL || m_NumUniqueSiteLoci == 0)
	if((Rslt=LoadLociIsolateMatrix())<eBSFSuccess)
		{
		Reset();
		return(Rslt);
		}
if((Rslt = CreateOutFile(pszFileExtn)) < eBSFSuccess)
	{
	Reset();
	return(Rslt);
	}

// first the heading, heading fields are simply each of the unique loci at which a SNP was present
uint32_t* pCurMatrixCell;
pCurMatrixCell = &m_pMatrix[1];
for(LociIdx = 1; LociIdx <= m_NumUniqueSiteLoci; LociIdx++, pCurMatrixCell++)
	m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs],",\"Loci:%u\"", *pCurMatrixCell);
m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], "\n");
// now each row
for(ReadsetIdx = 1; ReadsetIdx <= m_NumCultivars; ReadsetIdx++)
	{
	m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], "\"%s\"", LocateReadset(*pCurMatrixCell++)); // this the readset name
	for(LociIdx = 0; LociIdx < m_NumUniqueSiteLoci; LociIdx++, pCurMatrixCell++)		 // these are the matrix cells which contain SNP site classifications
		m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], ",%u", *pCurMatrixCell & 0x03);
	m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], "\n");
	if((Rslt = WriteOutFile()) < eBSFSuccess)
		return(Rslt);
	}
CloseOutFile();
return(eBSFSuccess);
}

int	
CSNPs2pgSNPs::IdentifySharedSiteLociDist(char *pszFileExtn)		// determine the SNP site loci sharing distributions and associate a PValue
{
int Rslt;
uint32_t FromReadsetID;
uint32_t ToReadsetID;
uint32_t FromFeatSNPIdx;

uint32_t CurSNPLoci;
uint32_t CurSNPChromID;
uint32_t NumSharedLoci;
double PValue;
double ExpRate;
tsIsolateFeatSNPs *pFromFeatSNP;
tsIsolateFeatSNPs* pToFeatSNP;
pFromFeatSNP = m_pIsolateFeatSNPs;
tsIsolateFeatSNPs *pLastFeatSNPloci = &m_pIsolateFeatSNPs[m_IsolateFeatSNPs - 1];
FromReadsetID = 0;
ToReadsetID = 0;

if((Rslt = CreateOutFile(pszFileExtn)) < eBSFSuccess)
	{
	Reset();
	return(Rslt);
	}
m_OutBuffOffs = sprintf(m_pszOutBuff,"\"SNP Loci\",\"PValue\",\"Instances\",\"Dirac Instances\",\"Readsets\"\n");
// determine PValues for observing numbers of shared loci between isolates
// this requires an estimate of the expected mutation rate
// to derive this I am firstly taking the sum of all observed mutations over all isolates and the sum of all sites over all isolates
// at which at least one of the isolates had a mutation at that site. In later processing, I iterate each of the sites showing at least one mutation
// and then subtract that site contribution to the overall sums count. The amended counts (sum of all mutations divided by sum of all sites) are used
// as the expected (background frequency) and the binomial distribution is used to calculate the PValue for each mutated loci along the gRNA.
//
// determine number of unique site loci in target genome
if(!m_NumUniqueSiteLoci)
	m_NumUniqueSiteLoci = GenNumUniqueSiteLoci();

if(m_SortOrder != 1)	// 0: unsorted, 1: SortSNPChromLociReadset, 2: SortSNPReadsetChromLoci
	{
	qsort(m_pIsolateFeatSNPs, m_IsolateFeatSNPs, sizeof(tsIsolateFeatSNPs), SortSNPChromLociReadset);
	m_SortOrder = 1;
	}

int32_t SumSites = m_NumCultivars * m_NumUniqueSiteLoci;
int32_t SumMutations = m_IsolateFeatSNPs;
int32_t AcceptedPValues;
int32_t NumSharedDiracLoci;
AcceptedPValues = 0;
NumSharedLoci = 0;
NumSharedDiracLoci = 0;
pFromFeatSNP = m_pIsolateFeatSNPs;
CurSNPLoci = pFromFeatSNP->SNPSSite.SNPLoci;
CurSNPChromID = pFromFeatSNP->SNPSSite.ChromID;
for(FromFeatSNPIdx = 0; FromFeatSNPIdx <= m_IsolateFeatSNPs; FromFeatSNPIdx++, pFromFeatSNP++)
	{
	if(FromFeatSNPIdx < m_IsolateFeatSNPs && pFromFeatSNP->SNPSSite.ChromID == CurSNPChromID && pFromFeatSNP->SNPSSite.SNPLoci == CurSNPLoci)
		{
		NumSharedLoci += 1;
		if(pFromFeatSNP->SNPSSite.DiracCnts[0] || pFromFeatSNP->SNPSSite.DiracCnts[1] || pFromFeatSNP->SNPSSite.DiracCnts[2] || pFromFeatSNP->SNPSSite.DiracCnts[3])
			NumSharedDiracLoci += 1;
		if(NumSharedLoci == 1)
			pToFeatSNP = pFromFeatSNP;
		continue;
		}
	ExpRate = ((double)SumMutations - NumSharedLoci) / ((double)SumSites - m_NumCultivars);
	PValue = 1.0 - m_Stats.Binomial(m_NumCultivars,NumSharedLoci, ExpRate);

	m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs],"%u,%.6f,%u,%u,%u\n",CurSNPLoci,PValue,NumSharedLoci,NumSharedDiracLoci,m_NumCultivars);
	if((Rslt = WriteOutFile()) < eBSFSuccess)
		return(Rslt);
	do
		{
		pToFeatSNP->Shared = NumSharedLoci;
		}
	while(++pToFeatSNP != pFromFeatSNP);
	CurSNPLoci = pFromFeatSNP->SNPSSite.SNPLoci;
	CurSNPChromID = pFromFeatSNP->SNPSSite.ChromID;
	if(pFromFeatSNP->SNPSSite.DiracCnts[0] || pFromFeatSNP->SNPSSite.DiracCnts[1] || pFromFeatSNP->SNPSSite.DiracCnts[2] || pFromFeatSNP->SNPSSite.DiracCnts[3])
		NumSharedDiracLoci = 1;
	else
		NumSharedDiracLoci = 0;
	NumSharedLoci = 1;
	pToFeatSNP = pFromFeatSNP;
	}
CloseOutFile();
return(eBSFSuccess);
}


int
CSNPs2pgSNPs::LoadSimilaritiesMatrix(void)				// generate a matrix (readsets x readsets) with cells containing scores
{
int Rslt;
double RowPValue;
double ColPValue;
double Score;

uint32_t* pRowMatrixCell;
uint32_t* pColMatrixCell;

uint32_t RowReadsetIdx;

uint32_t ColReadsetIdx;


uint32_t LociIdx;

if(m_pMatrix == NULL || m_NumUniqueSiteLoci == 0)
	if((Rslt = LoadLociIsolateMatrix()) < eBSFSuccess){
		Reset();
		return(Rslt);
		}

if(m_pSimilaritiesMatrix != NULL)
	{
	delete[]m_pSimilaritiesMatrix;
	m_pSimilaritiesMatrix = NULL;
	}

m_SimilaritiesMatrixCols = m_NumCultivars;
m_SimilaritiesMatrixRows = m_NumCultivars;
m_pSimilaritiesMatrix = new double[(size_t)m_SimilaritiesMatrixCols * m_SimilaritiesMatrixRows];
memset(m_pSimilaritiesMatrix, 0, (size_t)m_SimilaritiesMatrixCols * m_SimilaritiesMatrixRows * sizeof(double));

// for each readset sum number of SNPs in that readset - assuming random distribution then can determine the probability of observing SNP at any loci
int SumReadsetSNPs;
int RowSumReadsetSNPs;
int ColSumReadsetSNPs;
for(RowReadsetIdx = 1; RowReadsetIdx <= m_SimilaritiesMatrixRows; RowReadsetIdx++)
	{
	pRowMatrixCell = &m_pMatrix[RowReadsetIdx * m_MatrixCols];
	if(RowReadsetIdx != *pRowMatrixCell++)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Inconsistency between Similarities and Isolates matrix ReadsetIDs");
		Reset();
		return(-1);
		}
	if(!(RowReadsetIdx % 500))
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Processing isolate %u",RowReadsetIdx);
	SumReadsetSNPs = 0;
	for(LociIdx = 1; LociIdx < m_MatrixCols; LociIdx++)
		SumReadsetSNPs += *pRowMatrixCell++;
	for(ColReadsetIdx = 1; ColReadsetIdx <= m_SimilaritiesMatrixCols; ColReadsetIdx++) 
		{
		RowSumReadsetSNPs = SumReadsetSNPs;
		RowPValue = (double)RowSumReadsetSNPs/(m_MatrixCols-1);
		pColMatrixCell = &m_pMatrix[ColReadsetIdx * m_MatrixCols];
		
		if(ColReadsetIdx != *pColMatrixCell++)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Inconsistency between Similarities and Isolates matrix ReadsetIDs");
			Reset();
			return(-1);
			}
		ColSumReadsetSNPs = 0;
		for(LociIdx = 1; LociIdx < m_MatrixCols; LociIdx++)
			ColSumReadsetSNPs += *pColMatrixCell++;
		ColPValue = (double)ColSumReadsetSNPs / (m_MatrixCols-1);

		// iterate base x base and score
		Score = 1.0;
		pRowMatrixCell = &m_pMatrix[(RowReadsetIdx * m_MatrixCols) + 1];
		pColMatrixCell = &m_pMatrix[(ColReadsetIdx * m_MatrixCols) + 1];
		for(LociIdx = 1; LociIdx < m_MatrixCols; LociIdx++, pRowMatrixCell++, pColMatrixCell++)
			{
			if(*pRowMatrixCell == 0 && *pColMatrixCell == 0) // both match the reference
				continue;
			if(*pRowMatrixCell == 0 && *pColMatrixCell != 0) // column has a SNP
				{
				Score *= 1.0 - ColPValue;
				ColPValue = (double)(--ColSumReadsetSNPs) / (m_MatrixCols - 1);
				continue;
				}

			if(*pRowMatrixCell != 0 && *pColMatrixCell == 0) // row has a SNP
				{
				Score *= 1.0 - RowPValue;
				RowPValue = (double)(--RowSumReadsetSNPs) / (m_MatrixCols);
				continue;
				}

			if(*pRowMatrixCell != 0 && *pColMatrixCell != 0) // both row and column have SNP
				{
				ColPValue = (double)(--ColSumReadsetSNPs) / (m_MatrixCols);
				RowPValue = (double)(--RowSumReadsetSNPs) / (m_MatrixCols);
				continue;
				}
			}
		m_pSimilaritiesMatrix[((RowReadsetIdx - 1) * m_SimilaritiesMatrixCols) + (ColReadsetIdx - 1)] = Score;
		}
	}
gDiagnostics.DiagOut(eDLInfo, gszProcName, "Processed %u isolates",m_SimilaritiesMatrixRows);
return(eBSFSuccess);
}

int
CSNPs2pgSNPs::ReportSimilaritiesMatrix(char* pszFileExtn)	// report readset vs readset  similarities matrix using this file extension, will generate readset vs readset matrix if not already generated
{
int Rslt;
uint32_t ColReadsetIdx;
uint32_t RowReadsetIdx;
double* pCurMatrixCell;

if(m_pSimilaritiesMatrix == NULL)
	if((Rslt = LoadSimilaritiesMatrix()) < eBSFSuccess)
		{
		Reset();
		return(Rslt);
		}
if((Rslt = CreateOutFile(pszFileExtn)) < eBSFSuccess)
	{
	Reset();
	return(Rslt);
	}

// first the heading, heading fields are simply each of the unique loci at which a SNP was present
for(ColReadsetIdx = 1; ColReadsetIdx <= m_NumCultivars; ColReadsetIdx++)
	m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], ",\"%s\"", LocateReadset(ColReadsetIdx));
m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], "\n");
// now each row
pCurMatrixCell = m_pSimilaritiesMatrix;
for(RowReadsetIdx = 1; RowReadsetIdx <= m_NumCultivars; RowReadsetIdx++)
	{
	m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], "\"%s\"", LocateReadset(RowReadsetIdx)); // this the readset name
	for(ColReadsetIdx = 0; ColReadsetIdx < m_NumCultivars; ColReadsetIdx++, pCurMatrixCell++)		 // these are the matrix cells containing scores
		m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], ",%.9f", *pCurMatrixCell);
	m_OutBuffOffs += sprintf(&m_pszOutBuff[m_OutBuffOffs], "\n");
	if((Rslt = WriteOutFile()) < eBSFSuccess)
		return(Rslt);
	}
CloseOutFile();
return(eBSFSuccess);

}

// SortSNPChromLoci
// Sort tsIsolateFeatSNPs SNP loci by ascending chrom, loci, readset
int
CSNPs2pgSNPs::SortSNPChromLociReadset(const void* arg1, const void* arg2)
	{
	tsIsolateFeatSNPs* pEl1 = (tsIsolateFeatSNPs*)arg1;
	tsIsolateFeatSNPs* pEl2 = (tsIsolateFeatSNPs*)arg2;

	if(pEl1->SNPSSite.ChromID < pEl2->SNPSSite.ChromID)
		return(-1);
	if(pEl1->SNPSSite.ChromID > pEl2->SNPSSite.ChromID)
		return(1);

	if(pEl1->SNPSSite.SNPLoci < pEl2->SNPSSite.SNPLoci)
		return(-1);
	if(pEl1->SNPSSite.SNPLoci > pEl2->SNPSSite.SNPLoci)
		return(1);

	if(pEl1->SNPSSite.ReadsetID < pEl2->SNPSSite.ReadsetID)
		return(-1);
	if(pEl1->SNPSSite.ReadsetID > pEl2->SNPSSite.ReadsetID)
		return(1);

	return(0);
	}

// SortSNPLoci
// Sort tsIsolateFeatSNPs SNP loci by ascending readset, chrom, loci
int
CSNPs2pgSNPs::SortSNPReadsetChromLoci(const void* arg1, const void* arg2)
	{
	tsIsolateFeatSNPs* pEl1 = (tsIsolateFeatSNPs*)arg1;
	tsIsolateFeatSNPs* pEl2 = (tsIsolateFeatSNPs*)arg2;

	if(pEl1->SNPSSite.ReadsetID < pEl2->SNPSSite.ReadsetID)
		return(-1);
	if(pEl1->SNPSSite.ReadsetID > pEl2->SNPSSite.ReadsetID)
		return(1);

	if(pEl1->SNPSSite.ChromID < pEl2->SNPSSite.ChromID)
		return(-1);
	if(pEl1->SNPSSite.ChromID > pEl2->SNPSSite.ChromID)
		return(1);

	if(pEl1->SNPSSite.SNPLoci < pEl2->SNPSSite.SNPLoci)
		return(-1);
	if(pEl1->SNPSSite.SNPLoci > pEl2->SNPSSite.SNPLoci)
		return(1);
	return(0);
	}