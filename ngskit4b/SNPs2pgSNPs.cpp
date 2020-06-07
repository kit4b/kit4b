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
	{3,'E',"Gle"},
	{4,'F',"Phe"},
	{5,'G',"Gly"},
	{6,'H',"His"},
	{7,'I',"Lle"},
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
	struct arg_dbl* minalleleprop = arg_dbl0("P", "minalleleprop", "<dbl>", "SNP minimum allele proportion of loci site coverage (default 0.025, range 0.01..0.25)");
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
			if(MinAlleleProp > 0.25)
				MinAlleleProp = 0.25;

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

		strcpy(szOutpgSNPsFile, outpgsnps->filename[0]);
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
			gDiagnostics.DiagOutMsgOnly(eDLInfo, gszProcName, "Output to '%s' is in VCF 4.1 format", szOutpgSNPsFile);
		else
			gDiagnostics.DiagOutMsgOnly(eDLInfo, gszProcName, "Output to '%s' is in UCSC pgSNP format", szOutpgSNPsFile);

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
Rslt = eBSFSuccess;

Rslt = pSNPs2pgSNPs->Process(Mode, bAllowInferenced, LocalSeqErrRate,MinCoverage, MinAlleleProp,PValueThres, pszTrackName, pszAssemblyName,pszExperimentDescr, SetOp, NumInSetA, ppszSetA,NumInSetB, ppszSetB, pszGFFFile, pszSNPFile, pszOutFile);
if(pSNPs2pgSNPs != NULL)
	delete pSNPs2pgSNPs;
return(Rslt);
}


CSNPs2pgSNPs::CSNPs2pgSNPs(void)
{
m_pSNPSites = NULL;
m_pCSV = NULL;
m_pszLineBuff = NULL;
m_pszSiteDistBuff = NULL;
m_pGapCntDist = NULL;
m_pSharedSNPsDist = NULL;
m_pIndividualSNPsDist = NULL;
m_pGFFFile = NULL;
m_pFeatures = NULL;
m_hOutpgSNPs = -1;
m_hOutSiteDist = -1; 
m_hOutFeatures = -1;
Reset();
}


CSNPs2pgSNPs::~CSNPs2pgSNPs(void)
{
if(m_pGFFFile != NULL)
	delete m_pGFFFile;
if(m_pCSV != NULL)
	delete m_pCSV;
if(m_hOutpgSNPs!=-1)
	close(m_hOutpgSNPs);
if(m_hOutFeatures!=-1)
	close(m_hOutFeatures);
if (m_pszLineBuff != NULL)
	delete []m_pszLineBuff;
if(m_pszSiteDistBuff != NULL)
	delete []m_pszSiteDistBuff;
if(m_pSNPSites != NULL)
	delete m_pSNPSites;
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
if (m_hOutpgSNPs != -1)
	{
	if (m_LineBuffOffs && m_pszLineBuff != NULL)
		CUtility::SafeWrite(m_hOutpgSNPs, m_pszLineBuff, m_LineBuffOffs);
#ifdef _WIN32
	_commit(m_hOutpgSNPs);
#else
	fsync(m_hOutpgSNPs);
#endif
	close(m_hOutpgSNPs);
	m_hOutpgSNPs = -1;
	}

if(m_hOutFeatures != -1)
	{
#ifdef _WIN32
	_commit(m_hOutFeatures);
#else
	fsync(m_hOutFeatures);
#endif
	close(m_hOutFeatures);
	m_hOutFeatures = -1;
	}

if(m_hOutSiteDist != -1)
	{
	if(m_SiteDistOffs && m_pszSiteDistBuff != NULL)
		CUtility::SafeWrite(m_hOutSiteDist, m_pszSiteDistBuff, m_SiteDistOffs);
#ifdef _WIN32
	_commit(m_hOutSiteDist);
#else
	fsync(m_hOutSiteDist);
#endif
	close(m_hOutSiteDist);
	m_hOutSiteDist = -1;
	}



if(m_pszLineBuff != NULL)
	{
	delete []m_pszLineBuff;
	m_pszLineBuff = NULL;
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

m_szGFFFile[0] = '\0';
m_AllocdLineBuff = 0;
m_LineBuffOffs = 0;
m_szInSNPsFile[0] = '\0';
m_szOutpgSNPsFile[0] = '\0';
m_szOutSiteDistFile[0] = '\0';

m_AllocdSNPSites = 0;
m_NumSNPs = 0;
m_SNPID = 0;

m_MaxSNPgap = 0;
m_MaxSharedSNPs = 0;
m_SiteDistOffs = 0;
m_AllocdSiteDistBuff = 0;

m_NumFeatures = 0;
m_AllocFeatures = 0;
m_bRptHdr = true;
}



int 
CSNPs2pgSNPs::ProcessSnpmarkersSNPs(void)
{
	int Rslt;
	int NumFields;
	int CultivarIdx;
	tsSCultivar* pCultivar;
	uint32_t PrevSNPloci;
	uint32_t EstNumSNPs;
	uint32_t NumSNPsParsed;
	uint32_t RowNumber;
	bool bMarkerFormat;
	int ExpNumFields;
	CStats Stats;
	int NumSpeciesWithCnts;
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
		char* pszTxt;
		m_pCSV->GetText(1, &pszTxt);			// get ":TargSeq"
		strncpy(pSNPSite->szChrom, pszTxt, sizeof(pSNPSite->szChrom));
		pSNPSite->szChrom[sizeof(pSNPSite->szChrom) - 1] = '\0';
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
		m_pCSV->GetInt(4, &NumSpeciesWithCnts);			// get "NumSpeciesWithCnts"
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

			PValue = 1.0 - Stats.Binomial(pCultivar->TotBaseCnts, pCultivar->TotMismatches, pCultivar->Score);
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
				PValue = min(0.5, 1.0 - Stats.Binomial(pCultivar->TotBaseCnts, pCultivar->BaseCnts[Idx], pCultivar->Score));
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

		if(bHaveSNP)
			{
			if(PrevSNPloci != 0xffffffff)
				{
				uint32_t CurGap;
				CurGap = pSNPSite->SNPLoci - PrevSNPloci;
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
			pSNPSite->bSNPPlaceholder = true;

		pSNPSite->SNPId = ++m_SNPID;		// all SNPs have an identifier even if only a placeholder
		if(m_bRptHdr && m_NumSNPs)
			{
			Report(true, pSNPSite);
			m_bRptHdr = false;
			}
		Report(false, pSNPSite);
	}

if (m_LineBuffOffs)
	{
	CUtility::SafeWrite(m_hOutpgSNPs, m_pszLineBuff, m_LineBuffOffs);
	m_LineBuffOffs = 0;
	}

#ifdef _WIN32
	_commit(m_hOutpgSNPs);
#else
	fsync(m_hOutpgSNPs);
#endif
	close(m_hOutpgSNPs);
	m_hOutpgSNPs = -1;


if(m_SiteDistOffs)
	{
	CUtility::SafeWrite(m_hOutSiteDist, m_pszSiteDistBuff, m_SiteDistOffs);
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

CStats Stats;

EstNumSNPs = m_pCSV->EstNumRows();
NumSNPsParsed = 0;
RowNumber = 0;
ExpNumFields = 0;
m_NumCultivars = 0;
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
		m_NumCultivars = 1;
		bMarkerFormat = false;
		ExpNumFields = NumFields;
		}

	if (RowNumber == 1)
		{
		if (m_pCSV->IsLikelyHeaderLine())
			continue;
		else
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Expected CSV file '%s' first line to be a header line with  fully quoted field names", m_szInSNPsFile);
			Reset();
			return(eBSFerrFieldCnt);
			}
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
		strcpy(m_Cultivars[0].szName, (char*)"Unknown");
		}

	char* pszTxt;

	memset(pSNPSite, 0, sizeof(SNPSite));
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

	m_pCSV->GetText(4, &pszTxt);			// get "Chrom"
	strncpy(pSNPSite->szChrom, pszTxt, sizeof(pSNPSite->szChrom));
	pSNPSite->szChrom[sizeof(pSNPSite->szChrom) - 1] = '\0';
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

	PValue = 1.0 - Stats.Binomial(pSNPSite->TotBaseCnts, pSNPSite->TotMismatches, pSNPSite->LocalSeqErrRate);
	if(PValue > m_PValueThres)
		continue;

	bHaveSNP = false;
	for (int Idx = 0; Idx < 4; Idx++)
		{
		if(m_MinAlleleProp > 0.0 && (pSNPSite->BaseCnts[Idx].Cnts == 0 || (((double)pSNPSite->BaseCnts[Idx].Cnts * 4) / (double)pSNPSite->TotBaseCnts) < m_MinAlleleProp))
			{
			pSNPSite->BaseCnts[Idx].Cnts = 0;
			pSNPSite->BaseCnts[Idx].QScore = 0.0;
			continue;
			}

		PValue = 1.0 - Stats.Binomial(pSNPSite->TotBaseCnts, pSNPSite->BaseCnts[Idx].Cnts, pSNPSite->LocalSeqErrRate);
		if(PValue > m_PValueThres)
			{
			pSNPSite->BaseCnts[Idx].Cnts = 0;
			pSNPSite->BaseCnts[Idx].QScore = 0.0;
			}
		else
			{
			if((pSNPSite->BaseCnts[Idx].Cnts / (double)pSNPSite->TotBaseCnts) >= m_DiracThres)
				pSNPSite->DiracCnts[Idx] = 1;
			pSNPSite->BaseCnts[Idx].QScore = PValue; 
			if(Idx != pSNPSite->RefBase)
				bHaveSNP = true;
			}
		}
	if(!bHaveSNP)
		continue;

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
	if(m_NumFeatures > 0)
		{
				// locate the feature in which this SNP is located
		for(FeatureIdx = 0; FeatureIdx < m_NumFeatures; FeatureIdx++, pFeature++)
			{
			if(stricmp(pSNPSite->szChrom, pFeature->szChrom) && stricmp(pSNPSite->szChrom, "NC_045512.2"))
				continue;
			if(pSNPSite->SNPLoci >= pFeature->FeatStart && pSNPSite->SNPLoci <= pFeature->FeatEnd)
				break;
			}
		if(FeatureIdx < m_NumFeatures)
			{
			pFeature->TotSNPs += 1;
			for(FieldIdx = 0; FieldIdx <= 4; FieldIdx++)
				pFeature->DiracCnts[FieldIdx] += pSNPSite->DiracCnts[FieldIdx];
		
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

				for(AlignCodonIdx = 0; AlignCodonIdx < 64; AlignCodonIdx++)
					{
					// what count threshold should be used before treating as noise and thus counts to be ignored???
					// here I'm requiring at least 2 codon instances covering site with a PValue <= 2*m_PValueThres which is slightly better than guessing
					if(AlignCodons[AlignCodonIdx] < 2)
						continue;
					PValue = 1.0 - Stats.Binomial(SumCodonCnts, AlignCodons[AlignCodonIdx], pSNPSite->LocalSeqErrRate);
					if(PValue > (m_PValueThres * 2))
						continue;

					// accepting as codon change
					pFeature->FromAminoAcidChanges[gXlateAminoAcids[RefCodonIdx].SynGroup] += 1;
					pFeature->ToAminoAcidChanges[gXlateAminoAcids[RefCodonIdx].SynGroup][gXlateAminoAcids[AlignCodonIdx].SynGroup] += 1;

					PolypeptideCnts[gXlateAminoAcids[AlignCodonIdx].SynGroup] += AlignCodons[AlignCodonIdx];
					if(gXlateAminoAcids[AlignCodonIdx].SynGroup == gXlateAminoAcids[RefCodonIdx].SynGroup)
						{
						bIsSynonymous = true;
						pFeature->NumSiteSynonCodons++;
						// attempt to see if can pickup on codon bias adaptation towards humans more abundant tRNAs
						if(gXlateAminoAcids[AlignCodonIdx].HumanFreqPerK > gXlateAminoAcids[RefCodonIdx].HumanFreqPerK)
							pFeature->NumPosBiasedCodons++;
						else
							if(gXlateAminoAcids[AlignCodonIdx].HumanFreqPerK < gXlateAminoAcids[RefCodonIdx].HumanFreqPerK)
								pFeature->NumNegBiasedCodons++;
						}
					else
						{
						if(gXlateAminoAcids[AlignCodonIdx].bWobble)
							bIsWobble = true;
						else
							bIsNonSynonymous = true;
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
	}
if (m_LineBuffOffs)
	{
	CUtility::SafeWrite(m_hOutpgSNPs, m_pszLineBuff, m_LineBuffOffs);
	m_LineBuffOffs = 0;
	}
#ifdef _WIN32
_commit(m_hOutpgSNPs);
#else
fsync(m_hOutpgSNPs);
#endif
return(0);
}

// following function is experimental, enabling reporting in three different file formats
// pgSNP, VCF4.1, or CSV
int
CSNPs2pgSNPs::Report(bool bHeader,		// if true then header line(s) to be generated, otherwise SNPs or alleles
					 tsSNPSSite* pSNPSite)
{
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
			m_LineBuffOffs = sprintf(m_pszLineBuff, "track type=pgSnp visibility=3 db=%s name=\"%s\" description=\"%s\"\n", m_SpecAssemblyName, m_szTrackName, m_szDescription);
			m_LineBuffOffs += sprintf(&m_pszLineBuff[m_LineBuffOffs], "browser position %s:%d-%d", szChrom, StartLoci, EndLoci);
			break;

		case eRMFvcf:						// report in VCF 4.1 format
			m_LineBuffOffs = sprintf(m_pszLineBuff, "##fileformat=VCFv4.1\n##source=ngskit4b%s\n##reference=%s\n##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n",
										kit4bversion, m_szInSNPsFile);
			m_LineBuffOffs += sprintf(&m_pszLineBuff[m_LineBuffOffs], "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
			break;

		}
	CUtility::SafeWrite(m_hOutpgSNPs, m_pszLineBuff, m_LineBuffOffs);
	m_LineBuffOffs = 0;
	}
else
	{
	if(m_hOutSiteDist != -1 && m_pszSiteDistBuff != NULL)
		{
		if((m_SiteDistOffs + 500) >= m_AllocdSiteDistBuff)
			{
			CUtility::SafeWrite(m_hOutSiteDist, m_pszSiteDistBuff, m_SiteDistOffs);
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

	if ((m_LineBuffOffs + 500) >= m_AllocdLineBuff)
		{
		CUtility::SafeWrite(m_hOutpgSNPs, m_pszLineBuff, m_LineBuffOffs);
		m_LineBuffOffs = 0;
		}
	switch (m_ReportFormat) {
		case eRMFpgSNP:						// report in pgSNP format
			// a special case: UCSC references the reference SARS-CoV-2 as being NC_O45512v2,not NC_O45512.2
			// so change NC_O45512.2 to be NC_O45512v2 when reporting as pgSNP format

			if(!stricmp(pSNPSite->szChrom, "NC_045512.2"))
				strcpy(szChrom, "NC_045512v2");
			else
				strcpy(szChrom, pSNPSite->szChrom);
			m_LineBuffOffs += sprintf(&m_pszLineBuff[m_LineBuffOffs], "\n%s\t%d\t%d\t", szChrom, pSNPSite->SNPLoci, pSNPSite->SNPLoci + 1);
			NumAlleles = 0;
			for (int Idx = 0; Idx < 4; Idx++)
				{
				if (pSNPSite->BaseCnts[Idx].Cnts > 0)
					{
					if (m_pszLineBuff[m_LineBuffOffs - 1] != '\t')
						m_pszLineBuff[m_LineBuffOffs++] = '/';
					switch (Idx) {
						case 0: Base = 'A'; break;
						case 1: Base = 'C'; break;
						case 2: Base = 'G'; break;
						case 3: Base = 'T'; break;
						case 4: Base = 'N'; break;
						}
					NumAlleles++;
					m_LineBuffOffs += sprintf(&m_pszLineBuff[m_LineBuffOffs], "%c", Base);
					}
				}
			m_LineBuffOffs += sprintf(&m_pszLineBuff[m_LineBuffOffs], "\t%d\t", NumAlleles);
			for (int Idx = 0; Idx < 4; Idx++)
				{
				if (pSNPSite->BaseCnts[Idx].Cnts > 0)
					{
					if (m_pszLineBuff[m_LineBuffOffs - 1] != '\t')
						m_pszLineBuff[m_LineBuffOffs++] = ',';
					m_LineBuffOffs += sprintf(&m_pszLineBuff[m_LineBuffOffs], "%d", pSNPSite->BaseCnts[Idx].Cnts);
					}
				}
			m_pszLineBuff[m_LineBuffOffs++] = '\t';

			// following is gross misuse of scores! But I want to let the viewer know the percentages of all isolates/cultivars/whatever at the SNP loci sharing at least one allelic variant
			// so first the proportion ( * 10)  of all sites is reported
			PercentScore = max(1, (int)(0.5 + ((double)(pSNPSite->TotMismatches * 100.0) / pSNPSite->TotBaseCnts)));
			m_LineBuffOffs += sprintf(&m_pszLineBuff[m_LineBuffOffs], "\t%d", PercentScore);
			// next the population size
			if(NumAlleles >= 2)
				m_LineBuffOffs += sprintf(&m_pszLineBuff[m_LineBuffOffs], ",%d", pSNPSite->TotBaseCnts);
			// next the number of samples which have at least one allelic variation
			if(NumAlleles >= 3)
				m_LineBuffOffs += sprintf(&m_pszLineBuff[m_LineBuffOffs], ",%d", pSNPSite->TotMismatches);
			if(NumAlleles == 4)
				m_LineBuffOffs += sprintf(&m_pszLineBuff[m_LineBuffOffs], ",0");
			break;

		case eRMFvcf:						// report in VCF 4.1 format
			NumAlleles = 0;
			SumQScores = 0;
			Depth = 0;
			if(!stricmp(pSNPSite->szChrom, "NC_045512.2"))
				strcpy(szChrom, "NC_045512v2");
			else
				strcpy(szChrom, pSNPSite->szChrom);
			m_LineBuffOffs += sprintf(&m_pszLineBuff[m_LineBuffOffs], "%s\t%u\tSNP%d\t%c\t", szChrom, pSNPSite->SNPLoci + 1, pSNPSite->SNPId, CSeqTrans::MapBase2Ascii(pSNPSite->RefBase));
			for (int Idx = 0; Idx < 4; Idx++)
				{
				Depth += pSNPSite->BaseCnts[Idx].Cnts;
				if(Idx == pSNPSite->RefBase)
					continue;
				if (pSNPSite->BaseCnts[Idx].Cnts > 0)
					{
					if (m_pszLineBuff[m_LineBuffOffs - 1] != '\t')
						m_pszLineBuff[m_LineBuffOffs++] = ',';
					switch (Idx) {
						case 0: Base = 'A'; break;
						case 1: Base = 'C'; break;
						case 2: Base = 'G'; break;
						case 3: Base = 'T'; break;
						case 4: Base = 'N'; break;
						}
					NumAlleles++;
					m_LineBuffOffs += sprintf(&m_pszLineBuff[m_LineBuffOffs], "%c", Base);
					SumQScores += min(99.0, ((m_PValueThres - pSNPSite->BaseCnts[Idx].QScore) * (100 / m_PValueThres)));
					}
				}
			SumQScores /= NumAlleles;
			m_LineBuffOffs += sprintf(&m_pszLineBuff[m_LineBuffOffs], "\t%d\tPASS\tAF=", (int)SumQScores);
			for (int Idx = 0; Idx < 4; Idx++)
				{
				if (Idx == pSNPSite->RefBase)
					continue;
				if (pSNPSite->BaseCnts[Idx].Cnts > 0)
					{
					if (m_pszLineBuff[m_LineBuffOffs - 1] != '=')
						m_pszLineBuff[m_LineBuffOffs++] = ',';
					m_LineBuffOffs += sprintf(&m_pszLineBuff[m_LineBuffOffs], "%.4f", pSNPSite->BaseCnts[Idx].Cnts/(double)Depth);
					}
				}
			m_LineBuffOffs += sprintf(&m_pszLineBuff[m_LineBuffOffs],";DP=%d\n", Depth);
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
int FeatIdx;
int BufIdx;
int FromPeptideIdx;
int ToPeptideIdx;
char szLineBuff[4000];

if(m_hOutFeatures == -1 || NumFeatures == 0 || pFeatures == NULL)
	return(0);
BufIdx = 0;
if(bHeader)
	{
	BufIdx = sprintf(szLineBuff, "\"Chrom\",\"Feature\",\"Start Loci\",\"End Loci\"");
	BufIdx += sprintf(&szLineBuff[BufIdx], ",\"SNPS\",\"Synonymous\",\"NonSynonymous\",\"Wobble\",\"Excl.Synonymous\",\"Excl.NonSynonymous\",\"Excl.Wobble\"");
	BufIdx += sprintf(&szLineBuff[BufIdx], ",\"SiteSynonymous\",\"NumPosFreqBiasedCodons\",\"NumNegFreqBiasedCodons\"");
	BufIdx += sprintf(&szLineBuff[BufIdx], ",\"DiracA\",\"DiracC\",\"DiracG\",\"DiracT\"");
	for(FromPeptideIdx = 0; FromPeptideIdx <= 20; FromPeptideIdx++)
		{
		BufIdx += sprintf(&szLineBuff[BufIdx],",\"From:%s\"", gAminoAcids[FromPeptideIdx].szPeptide);
		for(ToPeptideIdx = 0; ToPeptideIdx <= 20; ToPeptideIdx++)
			BufIdx += sprintf(&szLineBuff[BufIdx], ",\"To:%s\"", gAminoAcids[ToPeptideIdx].szPeptide);
		if((BufIdx + 2000) > sizeof(szLineBuff))
			{
			CUtility::SafeWrite(m_hOutFeatures, szLineBuff, BufIdx);
			BufIdx = 0;
			}
		}
	BufIdx += sprintf(&szLineBuff[BufIdx], "\n");
	}

tsSummaryFeatCnts* pFeature = pFeatures;

for(FeatIdx = 0; FeatIdx < NumFeatures; FeatIdx++, pFeature++)
	{
	if((BufIdx + 2000) > sizeof(szLineBuff))
		{
		CUtility::SafeWrite(m_hOutFeatures, szLineBuff, BufIdx);
		BufIdx = 0;
		}
	BufIdx += sprintf(&szLineBuff[BufIdx],"\"%s\",\"%s\",%u,%u",pFeature->szChrom, pFeature->szFeatName, pFeature->FeatStart, pFeature->FeatEnd);
	BufIdx += sprintf(&szLineBuff[BufIdx],",%u,%u,%u,%u", pFeature->TotSNPs,pFeature->TotSynonymous,pFeature->TotNonSynonymous,pFeature->TotWobble);
	BufIdx += sprintf(&szLineBuff[BufIdx], ",%u,%u,%u", pFeature->TotExclSynonymous, pFeature->TotExclNonSynonymous, pFeature->TotExclWobble);
	BufIdx += sprintf(&szLineBuff[BufIdx], ",%u,%u,%u", pFeature->NumSiteSynonCodons,pFeature->NumPosBiasedCodons, pFeature->NumNegBiasedCodons);
	BufIdx += sprintf(&szLineBuff[BufIdx], ",%u,%u,%u,%u", pFeature->DiracCnts[0], pFeature->DiracCnts[1], pFeature->DiracCnts[2], pFeature->DiracCnts[3]);
	for(FromPeptideIdx = 0; FromPeptideIdx <= 20; FromPeptideIdx++)
		{
		if((BufIdx + 2000) > sizeof(szLineBuff))
			{
			CUtility::SafeWrite(m_hOutFeatures, szLineBuff, BufIdx);
			BufIdx = 0;
			}
		BufIdx += sprintf(&szLineBuff[BufIdx], ",%u", pFeature->FromAminoAcidChanges[FromPeptideIdx]);
		for(ToPeptideIdx = 0; ToPeptideIdx <= 20; ToPeptideIdx++)
			BufIdx += sprintf(&szLineBuff[BufIdx],",%u",pFeature->ToAminoAcidChanges[FromPeptideIdx][ToPeptideIdx]);
		}
	BufIdx += sprintf(&szLineBuff[BufIdx], "\n");
	CUtility::SafeWrite(m_hOutFeatures, szLineBuff, BufIdx);
	BufIdx = 0;
	}

return(BufIdx);
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
strcpy(m_szOutpgSNPsFile, pszOutFile);
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


if ((m_pszLineBuff = new char[cAllocLineBuffSize]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to allocate memory for output line buffering -- %s", strerror(errno));
	return(eBSFerrMem);
	}
m_AllocdLineBuff = cAllocLineBuffSize;

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
m_ReportFormat = eRMFsnp::eRMFpgSNP;
if((Len=(int)strlen(pszOutFile)) >= 5 && !stricmp(&pszOutFile[Len-4],".vcf"))
	m_ReportFormat = eRMFsnp::eRMFvcf;



#ifdef _WIN32
if ((m_hOutpgSNPs = open(pszOutFile, (_O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC), (_S_IREAD | _S_IWRITE))) == -1)
#else
if ((m_hOutpgSNPs = open(pszOutFile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE)) == -1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to create/truncate output file - '%s' - %s", pszOutFile, strerror(errno));
	Reset();
	return(eBSFerrOpnFile);
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

if(Mode == eMpgSSNPKalign)
	{
	strcpy(m_szOutFeatures, pszOutFile);
	strcat(m_szOutFeatures,".feats.csv");
#ifdef _WIN32
	if((m_hOutFeatures = open(m_szOutFeatures, (_O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC), (_S_IREAD | _S_IWRITE))) == -1)
#else
	if((m_hOutFeatures = open(m_szOutFeatures, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE)) == -1)
#endif
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to create/truncate output file - '%s' - %s", m_szOutFeatures, strerror(errno));
		Reset();
		return(eBSFerrOpnFile);
		}
	}

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
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to open file: %s", pszCurSNPFile);
				Reset();
				return(Rslt);
				}
			Rslt = ProcessKalignSNPs();
			if(m_pCSV != NULL)
				m_pCSV->Close();
			if(Rslt != eBSFSuccess)
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Failed processing SNP file '%s'\n", pszCurSNPFile);
				Reset();
				return(Rslt);
				}
			}
		if(m_hOutpgSNPs != -1)
			{
			if(m_LineBuffOffs && m_pszLineBuff != NULL)
				CUtility::SafeWrite(m_hOutpgSNPs, m_pszLineBuff, m_LineBuffOffs);
		#ifdef _WIN32
			_commit(m_hOutpgSNPs);
		#else
			fsync(m_hOutpgSNPs);
		#endif
			close(m_hOutpgSNPs);
			m_hOutpgSNPs = -1;
			}
		if(m_NumFeatures && m_hOutFeatures != -1)
			ReportFeatSNPcnts(true,m_NumFeatures,m_pFeatures);
		if(m_hOutFeatures != -1)
			{
	#ifdef _WIN32
			_commit(m_hOutFeatures);
	#else
			fsync(m_hOutFeatures);
	#endif
			close(m_hOutFeatures);
			m_hOutFeatures = -1;
			}
		break;

	case eMpgSSNPmarkers:
		if((Rslt = m_pCSV->Open(pszSNPFile)) != eBSFSuccess)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to open file: %s", pszSNPFile);
			Reset();
			return(Rslt);
			}
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
			strcpy(m_szOutpgSNPsFile, pszOutFile);
			strcat(m_szOutpgSNPsFile, ".gapdist.csv");
#ifdef _WIN32
			if((m_hOutpgSNPs = open(m_szOutpgSNPsFile, (_O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC), (_S_IREAD | _S_IWRITE))) == -1)
#else
			if((m_hOutpgSNPs = open(m_szOutpgSNPsFile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE)) == -1)
#endif
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to create/truncate output file - '%s' - %s", m_szOutpgSNPsFile, strerror(errno));
				Reset();
				return(eBSFerrOpnFile);
				}
			m_LineBuffOffs = sprintf(m_pszLineBuff, "\"Len\",\"Instances\"\n");
			for(uint32_t Idx = 0; Idx <= m_MaxSNPgap; Idx++)
				m_LineBuffOffs += sprintf(&m_pszLineBuff[m_LineBuffOffs], "%d,%d\n", Idx + 1, m_pGapCntDist[Idx]);

			CUtility::SafeWrite(m_hOutpgSNPs, m_pszLineBuff, m_LineBuffOffs);
#ifdef _WIN32
			_commit(m_hOutpgSNPs);
#else
			fsync(m_hOutpgSNPs);
#endif
			close(m_hOutpgSNPs);
			m_hOutpgSNPs = -1;

			strcpy(m_szOutpgSNPsFile, pszOutFile);
			strcat(m_szOutpgSNPsFile, ".sharedsnpdist.csv");
#ifdef _WIN32
			if((m_hOutpgSNPs = open(m_szOutpgSNPsFile, (_O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC), (_S_IREAD | _S_IWRITE))) == -1)
#else
			if((m_hOutpgSNPs = open(m_szOutpgSNPsFile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE)) == -1)
#endif
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to create/truncate output file - '%s' - %s", m_szOutpgSNPsFile, strerror(errno));
				Reset();
				return(eBSFerrOpnFile);
				}
			m_LineBuffOffs = sprintf(m_pszLineBuff, "\"Shared\",\"Instances\"\n");
			for(uint32_t Idx = 0; Idx <= m_MaxSharedSNPs; Idx++)
				m_LineBuffOffs += sprintf(&m_pszLineBuff[m_LineBuffOffs], "%d,%d\n", Idx + 1, m_pSharedSNPsDist[Idx]);
			CUtility::SafeWrite(m_hOutpgSNPs, m_pszLineBuff, m_LineBuffOffs);
#ifdef _WIN32
			_commit(m_hOutpgSNPs);
#else
			fsync(m_hOutpgSNPs);
#endif
			close(m_hOutpgSNPs);
			m_hOutpgSNPs = -1;

			strcpy(m_szOutpgSNPsFile, pszOutFile);
			strcat(m_szOutpgSNPsFile, ".individualsnpdist.csv");
#ifdef _WIN32
			if((m_hOutpgSNPs = open(m_szOutpgSNPsFile, (_O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC), (_S_IREAD | _S_IWRITE))) == -1)
#else
			if((m_hOutpgSNPs = open(m_szOutpgSNPsFile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE)) == -1)
#endif
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: Unable to create/truncate output file - '%s' - %s", m_szOutpgSNPsFile, strerror(errno));
				Reset();
				return(eBSFerrOpnFile);
				}
			m_LineBuffOffs = sprintf(m_pszLineBuff, "\"SNPs\",\"Instances\"\n");
			for(uint32_t Idx = 0; Idx <= m_MaxIndividualSNPs; Idx++)
				m_LineBuffOffs += sprintf(&m_pszLineBuff[m_LineBuffOffs], "%d,%d\n", Idx, m_pIndividualSNPsDist[Idx]);
			CUtility::SafeWrite(m_hOutpgSNPs, m_pszLineBuff, m_LineBuffOffs);
#ifdef _WIN32
			_commit(m_hOutpgSNPs);
#else
			fsync(m_hOutpgSNPs);
#endif
			close(m_hOutpgSNPs);
			m_hOutpgSNPs = -1;
		}
	break;
	}


Reset();
return(Rslt);
}